% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene-relevance-methods.r
\docType{methods}
\name{Gene Relevance methods}
\alias{Gene Relevance methods}
\alias{featureNames,GeneRelevance-method}
\alias{featureNames<-,GeneRelevance,characterOrFactor-method}
\alias{dataset,GeneRelevance-method}
\alias{dataset<-,GeneRelevance-method}
\title{Gene Relevance methods}
\usage{
\S4method{featureNames}{GeneRelevance}(object)


  \S4method{featureNames}{GeneRelevance,characterOrFactor}(object) <- value

\S4method{dataset}{GeneRelevance}(object)

\S4method{dataset}{GeneRelevance}(object) <- value
}
\arguments{
\item{object}{\code{\link{GeneRelevance}} object}

\item{value}{A text vector (\code{\link{character}} or \code{\link{factor}})}
}
\description{
\code{featureNames <- ...} Can be used to set the gene names used for plotting
(e.g. if the data contains hardly readably gene or transcript IDs)
}
\examples{
data(guo_norm)
gr <- gene_relevance(DiffusionMap(guo_norm))
featureNames(gr)[[37]] <- 'Id2 (suppresses differentiation)'
# now plot it with the changed gene name(s)

}
\seealso{
\code{\link{gene_relevance}}, \link{Gene Relevance plotting}
}
