\name{PostProcessChain}
\alias{PostProcessChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computation for  maps of posterior probability of population membership}
\description{
 Computes posterior probabilities of population membership for each pixel
 of the spatial domain. 
}
\usage{
PostProcessChain(repdat, repmcmc,nxdom, nydom,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repdat}{ Path to input files directory }
  \item{repmcmc}{ Path to output files directory }
  \item{nxdom}{Number of pixel for discretization of the spatial domain
    in the horizontal direction}
  \item{nydom}{Number of pixel for discretization of the spatial domain
    in the vertical direction}
  \item{burnin}{Number of iterations of the chain to throwaway}
}
\details{ \code{repdat} directory should contain at least three files
  named   exactly :
  \code{coordinates.txt}, \code{genotypes.txt} and
  \code{allele.numbers.txt} and containing respectively the spatial
  coordinates, the genotypes and the number of alleles per locus.

    See format of files simulated by \code{\link{simFmodel}} for an
  example. 
  
}
\value{
  Posterior probability of population membership for each pixel are
  written in an ascii file called \code{proba.pop.membership.txt}
  (one column per population,
  \code{nclassmax} values are computed for each pixel, 
  pixels are stored  row-wise starting from bottom left).

  Another file called \code{proba.pop.membership.perm.txt} tries
  to get rid of label switching
  issues by labelling the population according to a fixed constraint.
  (This has has proved to be usefull with a small number of loci, (e.g. nloc=3), for well
  differentiated populations.)
}
\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }




\seealso{\code{\link{PlotTessellation}}}


\keyword{}
