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
    iterations, namely 1000*10 computed iterations.}
}
\details{ \file{path.data} directory should contain at least three files
  named   exactly :
  \file{coordinates.txt}, \file{genotypes.txt} and
  \file{allele.numbers.txt} and containing respectively the spatial
  coordinates, the genotypes and the number of alleles per locus.

    See format of files simulated by \code{\link{simFmodel}} for an
  example. 
  
}
\value{
  Posterior probability of population membership for each pixel are
  written in an ascii file called \file{proba.pop.membership.txt}
  (one column per population,
  \code{npopmax} values are computed for each pixel, 
  pixels are stored  row-wise starting from bottom left).

  Another file called \file{proba.pop.membership.perm.txt} tries
  to get rid of label switching
  issues by labelling the population according to a fixed constraint.
  (This has proved to be usefull with a small number of loci, (e.g. nloc=3), for well
  differentiated populations.)
}
\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }




\seealso{\code{\link{PlotTessellation}}}


\keyword{}
