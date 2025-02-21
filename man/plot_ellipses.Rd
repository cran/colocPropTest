% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plot_ellipses}
\alias{plot_ellipses}
\title{draw two ellipses}
\usage{
plot_ellipses(
  b1,
  vb1,
  b2,
  vb2,
  legend = c("inferred", "observed"),
  include_origin = FALSE,
  ...
)
}
\arguments{
\item{b1}{ellipse 1 centre (2d)}

\item{vb1}{ellipse 1 vcov matrix}

\item{b2}{ellipse 2 centre (2d)}

\item{vb2}{ellipse 2 vcov matrix}

\item{legend}{character vector length 2 naming ellipse 1 and 2}

\item{include_origin}{if TRUE, ensure plot includes (0,0)}

\item{...}{arguments passed to plot()}
}
\value{
draw ellipses on current graphics device
}
\description{
draw two ellipses
}
\examples{
plot_ellipses(b1=c(5,5), vb1=diag(2),
              b2=c(2,2), vb2=matrix( c(1,0.5,0.5,1), 2, 2 ),
              legend=c("circle", "ellipse"),
              include.origin=TRUE) 
}
\author{
Chris Wallace
}
