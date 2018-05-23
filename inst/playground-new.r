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
	ggplot(plot_data, aes(DC1, DC2, colour = Cluster)) + geom_point()
}


# Test code ----


import_from('destiny', c('DiffusionMap', 'gene_relevance'))
import_from('plotly', 'ggplotly')

scial <- readRDS('~/Analysis/destiny-paper/scialdone.rds')
if (!'dm_scial' %in% ls()) dm_scial <- DiffusionMap(scial)
if (!'gr_scial' %in% ls()) {
	gr_scial <- gene_relevance(dm_scial, smooth = FALSE)
	featureNames(gr_scial) <- Biobase::fData(scial)$Symbol
}
gm_scial <- plot_gr(dm_scial, gr_scial)
print(ggplotly(gm_scial + scale_colour_hue()))
