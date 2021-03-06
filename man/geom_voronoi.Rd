% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom-voronoi.r
\docType{data}
\name{geom_voronoi}
\alias{geom_voronoi}
\alias{StatVoronoi}
\title{Voronoi tesselation}
\format{An object of class \code{StatVoronoi} (inherits from \code{Stat}, \code{ggproto}, \code{gg}) of length 3.}
\usage{
geom_voronoi(mapping = NULL, data = NULL, stat = "voronoi",
  position = "identity", na.rm = FALSE, eps = 1e-09,
  normalize = FALSE, expand = 0, show.legend = NA,
  inherit.aes = TRUE, ...)

StatVoronoi
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[=aes]{aes()}} or
\code{\link[=aes_]{aes_()}}. If specified and \code{inherit.aes = TRUE} (the
default), it is combined with the default mapping at the top level of the
plot. You must supply \code{mapping} if there is no plot mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[=ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[=fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data.}

\item{stat}{The statistical transformation to use on the data for this
layer, as a string.}

\item{position}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{eps}{A value of epsilon used in testing whether a quantity is zero,
mainly in the context of whether points are collinear.
If anomalous errors arise, it is possible that these may averted
by adjusting the value of eps upward or downward.}

\item{normalize}{Should coordinates be normalized prior to calculations.
If \code{x} and \code{y} are in wildly different ranges
it can lead to tesselation and triangulation that seems off
when plotted without \code{\link[ggplot2]{coord_fixed}()}.
Normalization of coordinates solves this.
The coordinates are transformed back after calculations.}

\item{expand}{How much to expand the convex hull around the points.
Default: 0.1; \eqn{expand \times max(span(data$x), span(data$y))}}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[=borders]{borders()}}.}

\item{...}{Other arguments passed on to \code{\link[=layer]{layer()}}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour = "red"} or \code{size = 3}. They may also be parameters
to the paired geom/stat.}
}
\description{
This stat and geom allows you to display voronoi tesselation as polygons.
The computations are based on the \code{\link[deldir]{deldir}()} function.
}
\section{Aesthetics}{

Understands the following aesthetics.
Required aesthetics are in bold:

\itemize{
\item \strong{x}
\item \strong{y}
\item alpha
\item color
\item fill
\item linetype
\item size
}
}

\examples{
library(ggplot2)
ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_voronoi(aes(fill = Species)) +
  geom_point()

# Difference of normalize = TRUE
ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_voronoi(aes(fill = Species), normalize = TRUE) +
  geom_point()
}
\keyword{datasets}
