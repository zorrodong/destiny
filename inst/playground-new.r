library(Matrix)
library(ggplot2)
library(patchwork)

import_env <- 'env_playground'
if (!import_env %in% search()) {
	attach(new.env(), name = import_env)
	assign('import_from', function(pkg, names) {
		for (name in names) {
			val <- get(name, asNamespace(pkg), inherits = FALSE)
			assign(name, val, as.environment(import_env))
		}
	}, as.environment(import_env))
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
	gg <- ggplot(plot_data, aes(DC1, DC2, colour = Cluster, text = )) + geom_point()
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
	theme_void() + theme(axis.text.x = element_text(size = 8, angle = 0, debug = FALSE))


# Test code ----


import_from('destiny', c('DiffusionMap', 'gene_relevance', 'plot_gene_relevance'))
import_from('plotly', c('ggplotly', 'export'))
import_from('glue', 'glue')
import_from('SummarizedExperiment', 'rowData')
import_from('magrittr', '%>%')

scial <- readRDS('~/Analysis/destiny-paper/scialdone.rds')
# scial_blood <- subset(scial, rep(TRUE, nrow(scial)), scial$cluster %in% c('yellow', 'brown'))

# do_it <- function(dat, name, distnc) {
# 	dm <- DiffusionMap(dat, distance = distnc)
# 	gr <- gene_relevance(dm, smooth = FALSE)
# 	featureNames(gr) <- SummarizedExperiment::rowData(dat)$Symbol
# 	gm <- plot_gr(dm, gr)
# 	dm_ly <- ggplotly(plot(dm, 1:2, col_by = 'cluster') + scale_fill_cluster())
# 	gm_ly <- ggplotly(gm + scale_colour_hue())
# 	print(dm_ly)
# 	print(gm_ly)
# 	
# 	dir.create('mail', showWarnings = FALSE)
# 	base <- glue('mail/{name} {distnc}')
# 	export(dm_ly, glue('{base}.png'))
# 	export(gm_ly, glue('{base} gene relevance v2.png'))
# 	writeChar(jsonlite::toJSON(setNames(gm$ids, seq_len(min(length(gm$ids), 4))), pretty = TRUE), glue('{base} gene relevance v2.json'))
# 	
# 	gm
# }
# 
# if (!('maps' %in% ls())) maps <- list()
# for (distnc in c('rankcor', 'cosine', 'euclidean')) {
# if (!distnc %in% names(maps)) maps[[distnc]] <- list(
# 	normal = do_it(scial,     'normal', distnc),
# 	blood = do_it(scial_blood, 'blood', distnc)
# )
# }

if (!'dm_scial' %in% ls()) dm_scial <- DiffusionMap(scial, distance = 'rankcor')
if (!'gr_scial' %in% ls()) gr_scial <- gene_relevance(dm_scial)
if (!'gr_scial_unsmooth' %in% ls()) gr_scial_unsmooth <- gene_relevance(dm_scial, smooth = FALSE)
featureNames(gr_scial) <- featureNames(gr_scial_unsmooth) <- rowData(scial)$Symbol
gg_scial <- plot(dm_scial, 1:2, col_by = 'cluster') + scale_fill_cluster()
gm_scial <- plot_gene_relevance(gr_scial) + ggtitle('GR smoothed partials+F1')
gm_scial_unsmooth <- plot_gene_relevance(gr_scial_unsmooth, genes = 8) + ggtitle('GR smoothed F1')
gm_scial_louvain <- plot_gr(dm_scial, gr_scial_unsmooth) + ggtitle('GR louvain clustering')

print(ggplotly(gg_scial))
print(ggplotly(gm_scial))
print(ggplotly(gm_scial_unsmooth))
print(ggplotly(gm_scial_louvain))
ggsave(
	'mail/normal rankcor gene relevance v1.pdf',
	(plot(gr_scial, gm_scial$ids) | (gm_scial / gg_scial)) + plot_layout(widths = c(3, 1)),
	width = 12)
ggsave(
	'mail/normal rankcor gene relevance v1 unsmoothed.pdf',
	(plot(gr_scial, gm_scial_unsmooth$ids) | (gm_scial_unsmooth / gg_scial)) + plot_layout(widths = c(3, 1)),
	width = 12)
ggsave(
	'mail/normal rankcor gene relevance v2 louvain.pdf',
	(plot(gr_scial, unlist(gm_scial_louvain$ids)) | (gm_scial_louvain / gg_scial)) + plot_layout(widths = c(3, 1)),
	width = 12)

# Just turquoise ----

gr_turquoise <- subset(scial, rep(TRUE, nrow(scial)), scial$cluster %in% 'turquoise') %>% DiffusionMap(distance = 'rankcor') %>% gene_relevance()
featureNames(gr_turquoise) <- rowData(scial)$Symbol
print(plot(gr_turquoise, c('T', 'Mesp1', 'Evx1', 'Id3', 'Wfdc2')) / plot_gene_relevance(gr_turquoise))

# differential map ----

import_from('destiny', c('geom_voronoi', 'scale_fill_cube_helix'))
import_from('magrittr', '%>%')
import_from('purrr', c('map', 'pmap', 'imap', 'reduce'))
import_from('dplyr', c('arrange', 'select'))

plot_differential_maps <- function(genes, name, gr) {
	g_missing <- setdiff(genes, featureNames(gr))
	if (length(g_missing) > 0) message('Missing genes: ', paste(g_missing, collapse = ', '))
	genes <- intersect(genes, featureNames(gr))
	if (length(genes) == 0) {
		warning('No genes left')
		return()
	}
	dtm <- destiny:::differential_map(gr, genes, 1:2)
	dtm$scatters %>% arrange(Expression) %>% ggplot(aes(DC1, DC2)) +
		geom_voronoi(aes(fill = PartialsNorm)) +
		geom_point(aes(colour = Expression), shape = 20) +
		scale_colour_viridis_c() + scale_fill_cube_helix(discrete = FALSE, reverse = TRUE) +
		facet_wrap(~ Gene) + ggtitle(name) +
		scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
		theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank())
}

plot_all_differential_maps <- function(gr, gene_lists, path) {
	gm <- plot_gene_relevance(gr)
	p <- c(setNames(gene_lists, seq_along(gene_lists)), list(determined = head(gm$ids))) %>%
		imap(plot_differential_maps, gr = gr) %>%
		reduce(`+`) + plot_layout(ncol = 2L)
	ggsave(path, p, scale = 2)
	p
}
	

list(
	genes = list(gm_scial$ids, c('T', 'Mesp1', 'Evx1', 'Id3', 'Wfdc2')),
  name = c('Found genes', 'Turquoise cluster genes')
) %>%
	purrr::pmap(plot_differential_maps) %>%
	Reduce(`/`, .) %>%
	print()

gene_lists <-
	list(pseudospace = 4, pseudotime = 5) %>%
	purrr::imap(function(fig, name)
		lapply(1:3, function(n)
			readxl::read_xlsx(paste0('mail/nature18633-s', fig, '.xlsx'), n, skip = 1) %>%
			.$gene.name %>%
			head()
    )
	)

subsets_gr <- function(ss) {
	scial_subset <- scial %>% subset(select = cluster %in% ss)
	gr <- scial_subset %>% DiffusionMap(distance = 'rankcor') %>% gene_relevance()
	featureNames(gr) <- rowData(scial_subset)$Symbol
	gr
}


# pseudospace = 4, blue, mesodermal progenitor
# pseudotime = 7/8, yellow/brown, blod development

plot_all_differential_maps(gr_scial, gene_lists$pseudospace, 'mail/meso-pseudospace.pdf') %>% print()
plot_all_differential_maps(gr_scial, gene_lists$pseudotime,  'mail/blood-pseudotime.pdf') %>% print()

gr_scial_meso  <- subsets_gr(c('blue'))
gr_scial_blood <- subsets_gr(c('yellow', 'brown'))

plot_all_differential_maps(gr_scial_meso,  gene_lists$pseudospace, 'mail/meso-pseudospace-only.pdf') %>% print()
plot_all_differential_maps(gr_scial_blood, gene_lists$pseudotime,  'mail/blood-pseudotime-only.pdf')  %>% print()


# umap ----

import_from('destiny', c('find_knn', 'gene_relevance', 'plot_gene_relevance'))
import_from('dplyr', c('bind_cols'))

dists <- scial %>% assay('logcounts') %>% t() %>% find_knn(15) %>% .$dist_mat
if (!'tumap' %in% ls()) tumap <- uwot::tumap(dists) %>% `colnames<-`(paste0('umap', 1:2))
gg_umap <- list(tumap, colData(scial)) %>%
	lapply(as.data.frame) %>%
	bind_cols() %>%
	ggplot(aes(umap1, umap2, colour = cluster)) +
	  geom_point() +
		scale_colour_cluster() + ggtitle('UMAP')

if (!'gr_tumap' %in% ls()) gr_tumap <- gene_relevance(tumap, scial %>% assay('logcounts') %>% t() %>% as.matrix())
featureNames(gr_tumap) <- rowData(scial)$Symbol
gm_umap <- plot_gene_relevance(gr_tumap)

print(gg_umap / gm_umap)
ggsave('mail/umap.pdf', gg_umap / gm_umap, height = 10)
