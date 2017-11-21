library(mygene)
library(destiny)
library(dplyr)
library(biomaRt)

mouse <- useMart('ensembl', 'mmusculus_gene_ensembl')

gene_symbols <- function(ids)
	queryMany(ids, scopes = 'ensembl.gene', fields = 'symbol', returnall = TRUE)$response$symbol

gene_descr <- function(syms) {
	getBM(attributes = c('mgi_symbol', 'description'), filters = 'mgi_symbol', values = syms, mart = mouse) %>%
		with(setNames(description, mgi_symbol))
}

data(guo_norm)
dm <- DiffusionMap(guo_norm)

#phil
gr_phil <- gene_relevance(dm)

gr_plot_phil <- plot(gr_phil)
gene_descr(gr_plot_phil$ids) %>% as.list

#david
gr_david <- david_compute_gene_relevance(dm = dm)

gr_plot_david <- david_plot_gene_relevance_map(dm = gr_david)
gene_descr(gr_plot_david$ids) %>% as.list

# combined
gridExtra::grid.arrange(gr_plot_phil, gr_plot_david, ncol = 2)
