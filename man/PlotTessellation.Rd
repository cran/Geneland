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

\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }


\seealso{ \code{\link{PostProcessChain}}}


\keyword{}
