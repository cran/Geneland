\name{Plotntile}
\alias{Plotntile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of number of tiles along the MCMC run}
\description{Plot of number of tiles in the Poisson-Voronoi tessellation
  along the MCMC run}

\usage{
Plotntile(path.mcmc,burnin,printit=FALSE,file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path.mcmc}{ Path to output files directory }
   \item{printit}{Logical : if TRUE, figures are also printed}
  \item{file}{Character : Path to file where figures
    should be printed}
  \item{burnin}{An integer: number of saved iterations to discard for
    the representation of the histogram of the chain}
}

\author{Gilles Guillot  }

