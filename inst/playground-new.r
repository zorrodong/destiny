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


# GR plotting ----


import_from('destiny', c('eigenvectors', 'eigenvalues', 'dataset', 'find_knn'))
import_from('destiny', c('dataset_extract_doublematrix', 'get_smoothing', 'get_louvain_clusters'))

plot_gr <- function(dm, gr, n_gene_clusts = 500L, dims = 1:2, k = 30, distance = c('euclidean', 'cosine', 'rankcor')) {
	#d <- dataset_extract_doublematrix(dataset(dm))
	distance <- match.arg(distance)
	
	cell_clust_label <- get_louvain_clusters(dm@transitions)
	cell_clusts <- seq_len(max(cell_clust_label))
	
	best_genes <- lapply(cell_clusts, function(cc) {
		chunk <- gr@exprs[cell_clust_label == cc, ]
		dists <- find_knn(t(chunk), k, distance = distance)$dist_mat
		adj <- sparseMatrix(dists@i, p = dists@p, x = 1/dists@x, dims = dim(dists), symmetric = TRUE, index1 = FALSE)
		
		# no resolution settable, so using one of the hierarchical: cluster_walktrap, cluster_fast_greedy, cluster_edge_betweenness
		graph <- igraph::graph_from_adjacency_matrix(adj, 'undirected', weighted = TRUE)
		gene_clust_label <- as.integer(unclass(igraph::cut_at(igraph::cluster_walktrap(graph), n_gene_clusts)))
		gene_clusts <- seq_len(max(gene_clust_label))
		
		gene_clust_scores <- lapply(gene_clusts, function(cg) {
			gene_names <- Biobase::featureNames(gr)[gene_clust_label == cg]
			scores <- gr@partials_norm[cell_clust_label == cc, gene_clust_label == cg, drop = FALSE]
			gene_scores <- apply(scores, 2L, mean)
			sort(setNames(gene_scores, gene_names), decreasing = TRUE)
		})
		gene_clust_score_means <- sapply(gene_clusts, function(cg) mean(gene_clust_scores[[cg]]))
		names(gene_clust_scores[[which.max(gene_clust_score_means)]])
	})
	
	plot_data <- data.frame(gr@coords, Cluster = sapply(best_genes[cell_clust_label], paste, collapse = ', '))
	gg <- ggplot(plot_data, aes(DC1, DC2, colour = Cluster)) + geom_point()
	gg$ids <- best_genes
	gg
}


# Colors ----

import_from('RColorBrewer', 'brewer.pal')

prettier_colors <- c(
	setNames(brewer.pal(8, 'Dark2')[c(1,4,6,8)], c('turquoise', 'magenta', 'yellow', 'black')),
	setNames(brewer.pal(9, 'Set1')[-c(4,6)], c('red', 'blue', 'green', 'darkorange', 'brown', 'pink', 'grey'))
)
scale_colour_cluster <- function(...) scale_color_manual(values = prettier_colors, ...)
scale_fill_cluster <- function(...) scale_fill_manual(values = prettier_colors, ...)

ggplot(tibble::enframe(prettier_colors), aes(name, 1, fill = value)) +
	geom_col() + coord_fixed() + scale_fill_identity() + labs(x = NULL, y = NULL) +
	theme_void() + theme(axis.text.x = element_text())


# Test code ----


import_from('destiny', c('DiffusionMap', 'gene_relevance'))
import_from('plotly', c('ggplotly', 'export'))
import_from('glue', 'glue')

scial <- readRDS('~/Analysis/destiny-paper/scialdone.rds')
scial_blood <- subset(scial, rep(TRUE, nrow(scial)), scial$cluster %in% c('darkorange', 'brown'))

do_it <- function(dat, name, distnc) {
	dm_scial <- DiffusionMap(dat, distance = distnc)
	gr_scial <- gene_relevance(dm_scial, smooth = FALSE)
	featureNames(gr_scial) <- SummarizedExperiment::rowData(dat)$Symbol
	gm_scial <- plot_gr(dm_scial, gr_scial)
	gm_scial_ly <- ggplotly(gm_scial + scale_colour_hue())
	print(gm_scial_ly)
	
	dir.create('mail', showWarnings = FALSE)
	base <- glue('mail/{name} {distnc}')
	export(ggplotly(plot(dm_scial, 1:2, col_by = 'cluster') + scale_fill_cluster()), glue('{base}.png'))
	export(gm_scial_ly, glue('{base} gene relevance v2.png'))
	writeChar(jsonlite::toJSON(setNames(gm_scial$ids, 1:4), pretty = TRUE), glue('{base} gene relevance v2.json'))
}

for (distnc in c('euclidean', 'cosine', 'rankcor')) {
	do_it(scial,       'normal', distnc)
	do_it(scial_blood, 'blood',  distnc)
}
