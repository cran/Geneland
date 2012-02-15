\name{PlotTessellation}
\alias{PlotTessellation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maps of posterior probability  of membership}
\description{Plots maps of posterior probabilities  of population
  membership for each population}
\usage{
PlotTessellation(coordinates, path.mcmc,printit,path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
   columns and one line per individual.}
   \item{path.mcmc}{Character : Path to output files directory }
  \item{printit}{Logical : if TRUE, figures are also printed}
  \item{path}{Character : Path to directory where figures
      should be printed}
}

\author{Gilles Guillot  }


\seealso{ \code{\link{PostProcessChain}}}
