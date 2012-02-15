\name{Plotnpop}
\alias{Plotnpop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of number of populations along the MCMC run}
\description{Gives a plot of the number of populations along the MCMC run
}
\usage{
Plotnpop(path.mcmc,burnin,printit=FALSE,file,format="pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path.mcmc}{ Path to output files directory }
  \item{printit}{Logical : if TRUE, figures are also printed}
  \item{file}{Character : Path to file where figures
    should be printed}
  \item{format}{format of the output file, should be either \code{"ps"}
    or \code{"pdf"}}
  \item{burnin}{An integer: number of saved iterations to discard for
    the representation of the histogram of the chain}
}

\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, Mortier, Estoup, Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, 2005.}
\author{Gilles Guillot  }
 