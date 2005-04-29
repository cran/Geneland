\name{PostProcessChain}
\alias{PostProcessChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computation for  maps of posterior probability of population membership}
\description{
 Computes posterior probabilities of population membership for each pixel
 of the spatial domain. 
}
\usage{
PostProcessChain(coordinates,genotypes,allele.numbers,
path.mcmc,nxdom, nydom,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.}
  \item{genotypes}{Genotypes of individuals. A matrix with one line per
    individual and 2 columns per locus}
  \item{allele.numbers}{A vector of integer containing the number of possible
    allele for each locus}
  \item{path.mcmc}{ Path to output files directory }
  \item{nxdom}{Number of pixel for discretization of the spatial domain
    in the horizontal direction}
  \item{nydom}{Number of pixel for discretization of the spatial domain
    in the vertical direction}
  \item{burnin}{Number of iterations of the chain to throw away.
    WARNING : this argument should be given the number of stored
    iterations (and not the number of computed iterations which differ
    if \code{burnin} !=1). If you have
    \code{nit}=100000 and \code{thinning}=100, then only 1000 iterations
    are stored. Then \code{burnin}=10 will throw away 10 stored
    iterations, namely 100*10 computed iterations.}
}

\value{
  \item{Posterior probability of population membership for each pixel:}{
    They  are
    written in an ascii file called \file{proba.pop.membership.txt}
    (one column per population,
    \code{npopmax} values are computed for each pixel.
    Images in each column of \file{proba.pop.membership.txt} are stored
    column-wise starting from the bottom left pixel.
    First line of \file{proba.pop.membership.txt}  = bottom left pixel ,
    second line of \file{proba.pop.membership.txt} = upward neighboor of
    the previous pixel, etc...)
    
    Another file called \file{proba.pop.membership.perm.txt} tries
    to get rid of label switching
    issues by labelling the population according to a fixed constraint.
    (This has proved to be usefull with a small number of loci, (e.g. nloc=3), for well
    differentiated populations.)}
  
  
  \item{Posterior probability of population membership for each
    individual:}{ They are written in a file named
    \file{proba.pop.membership.indiv.txt}.
  }
  
  
  \item{Label of modal population for individuals:}{
    They are written in a file named \file{modal.pop.indiv.txt}. 
    See the example section  of function \code{\link{mcmcFmodel}} to see
    how they can be added in a plot.}
}
\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, Mortier, Estoup, Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, 2005.}
\author{Gilles Guillot  }




\seealso{\code{\link{PlotTessellation}}}


\keyword{}
