\name{PostProcessChain}
\alias{PostProcessChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computation for  maps of posterior probability of population membership}
\description{
 Computes posterior probabilities of population membership for each pixel
 of the spatial domain. 
}
\usage{
PostProcessChain(coordinates,
path.mcmc,nxdom, nydom,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.}

  \item{path.mcmc}{ Path to output files directory }

  \item{nxdom}{Number of pixel for discretization of the spatial domain
    in the horizontal direction}

  \item{nydom}{Number of pixel for discretization of the spatial domain
    in the vertical direction}

  \item{burnin}{Number of iterations of the chain to throw away.
    WARNING : this argument should be given the number of stored
    iterations (and not the number of computed iterations which differ
    if \code{thinning} !=1). If you have
    \code{nit}=100000 and \code{thinning}=100, then only 1000 iterations
    are stored. Then \code{burnin}=10 will throw away 10 stored
    iterations, namely 100*10 computed iterations.}
}

\value{

  \item{Posterior probability of population membership for each pixel:}{
    They  are
    written in an ascii file called \file{proba.pop.membership.txt}.
    Two first columns are coordinates of pixels then 
    one column per population (thus \code{npopmax} values are computed
  for each pixel). 
    Images in each column of \file{proba.pop.membership.txt} are stored
    starting from the bottom left pixel.
    First line of \file{proba.pop.membership.txt}  = bottom left pixel ,
    second line of \file{proba.pop.membership.txt} = upward neighboor of
    the previous pixel, etc...}
  
  

  \item{Posterior probability of population membership for each
    individual:}{ They are written in a file named
    \file{proba.pop.membership.indiv.txt}. One line per individuals and
    \code{2+npopmax} columns per individual. Two first columns are
    spatial coordinates.
  }
  
  

  \item{Label of modal population for pixels and individuals:}{
    They are written in files named \file{modal.pop.txt} and
    \file{modal.pop.indiv.txt}
    respectively. 
    See the example section  of function \code{\link{MCMC}} to see
    how they can be added in a plot.}
}
\references{G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
  model for landscape genetics. Genetics, 170, 1261-1280, 2005.

  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
  genetics. Molecular Ecology Notes, 5, 712-715, 2005.}

\author{Gilles Guillot  }




\seealso{\code{\link{PlotTessellation}}}

