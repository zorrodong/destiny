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
gr <- gene_relevance(dm)

(gr_plot <- plot(gr))
gene_descr(gr_plot$ids) %>% as.list
