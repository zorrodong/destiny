library(Matrix)
library(ggplot2)

if (!'env_playground' %in% search()) {
	attach(new.env(), name = 'env_playground')
	assign('import_from', function(pkg, names) {
		for (name in names) {
			val <- get(name, asNamespace(pkg), inherits = FALSE)
			assign(name, val, as.environment('env_playground'))
		}
	}, as.environment('env_playground'))
}


# gene_relevance ----


import_from('destiny', c('eigenvectors', 'eigenvalues', 'dataset', 'find_knn'))
import_from('destiny', c('dataset_extract_doublematrix', 'get_smoothing'))
import_from('matrixStats', 'rowMedians')
import_from('smoother', 'smth.gaussian')

gene_relevance <- function(coords, exprs, ..., k = 20L, dims = 1:2, distance = NULL, smooth = TRUE, verbose = FALSE) {
	if (is(coords, 'DiffusionMap')) {
		dm <- coords
		coords <- eigenvectors(dm)
		exprs <- dataset_extract_doublematrix(dataset(dm))
		weights <- eigenvalues(dm)[dims]
	}
	
	distance <- match.arg(distance, c('euclidean', 'cosine', 'rankcor'))
	coords_used <- coords[, dims, drop = FALSE]
	n_dims <- ncol(coords_used) # has to be defined early for `weights` default argument.
	smooth <- get_smoothing(smooth)
	
	if (is.null(colnames(exprs))) stop('The expression matrix columns need to be named but are NULL')
	if (n_dims != length(weights)) stop(n_dims, ' dimensions, but ', length(weights), ' weights were provided')
	
	nn_index <- find_knn(exprs, k, distance = distance)$index
	
	k <- ncol(nn_index)
	n_cells <- nrow(coords_used)
	n_genes <- ncol(exprs)
	partials <- array(
		NA,
		dim = c(n_cells, n_genes, n_dims),
		dimnames = list(NULL, colnames(exprs), if (is.character(dims)) dims else colnames(coords_used)))
	
	# a very small value to subtract from the differential
	small <- min(exprs[exprs != 0]) / length(exprs[exprs == 0])
	if (verbose) cat('Calculating expression differential\n')
	gene_differential <- function(expr_gene) {
		# Compute change in expression
		# Do not compute if reference is zero, could be drop-out
		expr_masked <- expr_gene
		expr_masked[expr_masked == 0] <- small
		differential_expr <- apply(nn_index, 2, function(nn) expr_gene[nn] - expr_masked)
		differential_expr[differential_expr == 0] <- NA  # Cannot evaluate partial
		#stopifnot(identical(dim(differential_expr), c(n_cells, k)))
		differential_expr
	}
	differential_exprs <- apply(exprs, 2L, gene_differential)
	#stopifnot(identical(dim(differential_exprs), c(n_cells * k, n_genes)))
	# apply only handles returning vectors, so we have to reshape the return value
	dim(differential_exprs) <- c(n_cells, k, n_genes)
	dimnames(differential_exprs)[[3L]] <- colnames(exprs)
	#stopifnot(identical(gene_differential(exprs[, 1L]), differential_exprs[, , 1L]))
	
	for (d in seq_len(n_dims)) {
		# Compute partial derivatives in direction of current dimension
		
		if (verbose) cat('Calculating partial derivatives of dimension ', d, '/', n_dims, '\n')
		# We could optionaly add normalization by max(coords_used[, d]) - min(coords_used[, d])
		differential_coord <- apply(nn_index, 2L, function(nn) coords_used[nn, d] - coords_used[, d])
		
		partials_unweighted <- apply(differential_exprs, 3L, function(grad_gene_exprs) {
			# Compute median of difference quotients to NN
			difference_quotients <- differential_coord / grad_gene_exprs
			# Only compute differential if at least two observations are present!
			stable_cells <- rowSums(!is.na(difference_quotients)) >= 2L
			ifelse(stable_cells, rowMedians(difference_quotients, na.rm = TRUE), NA)
		})
		colnames(partials_unweighted) <- colnames(exprs)
		
		if (!any(is.na(smooth))) {
			order_coor <- order(coords_used[, d])
			order_orig <- order(order_coor)
			partials_unweighted <- apply(partials_unweighted, 2L, function(gene_exprs) {
				ordered <- gene_exprs[order_coor]
				smoothed <- smth.gaussian(ordered, smooth[[1L]], smooth[[2L]], tails = TRUE)
				smoothed[order_orig]
			})
			colnames(partials_unweighted) <- colnames(exprs)
		}
		
		partials[, , d] <- weights[[d]] * partials_unweighted
	}
	
	# Compute norm over partial derivates: Frobenius
	partials_norm <- apply(partials, c(1, 2), function(z) sqrt(sum(z^2, na.rm = TRUE)))
	
	# Find outlier cells: Not in NN of more than 1 other cell
	# Remove these as they tend to receive very large norms
	#outliers <- sapply(seq_len(n_cells), function(cell) sum(nn_index == cell) > 1)
	#partials_norm[, outliers] <- NA
	
	# Prepare output
	colnames(partials_norm) <- colnames(partials)
	new(
		'GeneRelevance',
		coords = coords,
		exprs = exprs,
		partials = partials,
		partials_norm = partials_norm,
		nn_index = nn_index,
		dims = dims,
		smooth_window = smooth[[1L]],
		smooth_alpha  = smooth[[2L]],
		distance = distance)
}

# plot_gene_relevance ----


import_from('destiny', 'get_coords')
import_from('utils', 'head')

plot_gene_relevance <- function(relevance_map, ..., iter_smooth = 2L, genes = 5L, dims = 1:2, pal = palette()) {
	partials_norm <- relevance_map@partials_norm
	coords <- get_coords(relevance_map, dims)
	
	if (is.character(genes)) {
		found <- sapply(genes, function(id) length(grep(id, colnames(partials_norm))) > 0)
		gene_ids <- genes[found]
	} else if (length(genes) == 1L) {
		n_genes <- min(genes, ncol(relevance_map@exprs), na.rm = TRUE)
		# gene with max norm for each cell
		genes_max <- colnames(partials_norm)[apply(partials_norm, 1L, function(cell) which.max(cell))]
		counts <- as.data.frame(table(genes_max), stringsAsFactors = FALSE)
		n_genes <- min(n_genes, nrow(counts))
		gene_ids <- counts[order(counts$Freq, decreasing = TRUE)[1:n_genes], 'genes_max']
	} else {
		gene_ids <- colnames(partials_norm)[genes]
	}
	if (is.function(pal)) pal <- pal(length(gene_ids))
	
	num_top <- min(5L, length(gene_ids))
	top_n <- apply(partials_norm, 1L, function(cell) {
		idxs <- head(order(cell, decreasing = TRUE), num_top)
		names <- colnames(partials_norm)[idxs]
		txt <- sprintf('%s. %s (%.3f)', seq_len(num_top), names, cell[idxs])
		paste(txt, collapse = '\n')
	})
	
	# Plot a single map with cells coloured by gene which has 
	# the highest differential norm of all genes considered.
	
	max_gene_idx <- apply(partials_norm[, gene_ids, drop = FALSE], 1L, which.max)
	max_gene_idx[sapply(max_gene_idx, length) == 0] <- NA
	max_gene <- gene_ids[unlist(max_gene_idx)]
	# Label smoothing through graph structure
	for (i in seq_len(iter_smooth)) {
		max_gene <- apply(relevance_map@nn_index, 1, function(cell) {
			max_genes_nn <- unique(max_gene[cell])
			max_genes_nn_hist <- sapply(max_genes_nn, function(gene) sum(gene == max_gene[cell], na.rm = TRUE))
			names(max_genes_nn_hist)[which.max(max_genes_nn_hist)]
		})
	}
	# Add more than two DC and return data frame so that user
	# can easily rebuild relevance map on other DC combination than 1 and 2.
	rel_map_data <- cbind(as.data.frame(coords), Gene = factor(max_gene, levels = gene_ids), TopN = top_n)
	
	d1 <- colnames(coords)[[1]]
	d2 <- colnames(coords)[[2]]
	rel_map <- ggplot(rel_map_data, aes_string(x = d1, y = d2, colour = 'Gene', text = 'TopN')) +
		geom_point(alpha = .8) + 
		scale_color_manual(values = pal) +
		ggtitle(sprintf('Gene relevance map'))
	
	rel_map$ids <- gene_ids
	
	rel_map
}


# Test code ----


import_from('destiny', 'DiffusionMap')
import_from('plotly', 'ggplotly')

data(guo_norm)
dm_guo <- DiffusionMap(guo_norm)
gr_guo <- gene_relevance(dm_guo, smooth = FALSE)
gm_guo <- plot_gene_relevance(gr_guo, iter_smooth = 0, genes = NA)
print(ggplotly(gm_guo + scale_colour_hue()))


# playground ----

gr_guo@exprs

