\name{simIBD}
\alias{simIBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of georeferenced genotypes under an IBD + barrier
  model}

\description{Simulates coordinates and genotypes for a \code{npop}
  populations.
  Each population is supposed to be under an Isolation by Distance model
  and different populations are supposed to be separated by impermeable
  barriers. The barriers are given by a Poisson-Voronoi tessellation.}
\usage{
simIBD(nindiv,
coord.indiv = NULL,
coord.lim,
rate = NULL,
number.nuclei = NULL,
coord.nuclei = NULL,
color.nuclei = NULL,
allele.numbers,
model,
param,
npop,
seed.coord = NULL,
seed.tess = NULL,
seed.freq = NULL,
give.tess.grid=FALSE,
give.freq.grid = FALSE,
npix = NULL,
comp.Fst = FALSE,
comp.Dsigma2=FALSE,
comp.diff=FALSE,
width,
plot.pairs.borders=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nindiv}{Number of indivuals}
  \item{coord.indiv}{Coordinates of the individuals}
  \item{coord.lim}{Limits of the geographical domain. The domain is
    supposed to be rectangular and the limits are given as  (abs min, abs
    max, ord min, ord max)}
  \item{rate}{Rate of the Poisson process governing the hidden tessellation}
  \item{number.nuclei}{Number of nuclei in the tessellation (if given, 
    then \code{rate} is ignored)}
  \item{coord.nuclei}{Coordinates of the nuclei (the number of
    coordinates of the nuclei
    given here as a matrix has to comply with \code{number.nuclei} )}
  \item{color.nuclei}{Population membeship of the nuclei: a vector of
    integer of length \code{number of nuclei} whose values are between 1
  and \code{npop}}
  \item{allele.numbers}{A vecor giving the number of alleles observed at
    each locus}
  \item{model}{Model of spatial covariance function used for the
    underlying Gaussian fields (see documentation of the package
    \code{RandomFields}
  for details)}
  \item{param}{Parametersof spatial covariance function used for the
    underlying Gaussian fields (see documentation of the package
    \code{RandomFields}
  for details)}
  \item{npop}{Number of Populations}
  \item{seed.coord}{Random seed to initialise the simulation of the
    coordinates
  (mostly for debugging)}
  \item{seed.tess}{Random seed to initialise the simulation of the
    tessellation  (mostly for debugging)}
  \item{seed.freq}{ Random seed to initialise the simulation of the
    frequencies  (mostly for debugging)}
  \item{give.freq.grid}{Logical to tell whether frequencies  on a grid are
    also returned}
  \item{give.tess.grid}{Logical to tell whether population memberships of pixels  on a grid are also returned}
  \item{npix}{A vector of two integers telling how many horizontal and
    vertical pixel should contain the grid for the graphical
    representations}
  \item{comp.Fst}{Logical to tell whether Fst, Fis and Fit should be
    computed}
  \item{comp.Dsigma2}{Logical to tell whether IBD index Dsgma2 should be
    computed}
  \item{comp.diff}{Logical to tell whether the local differentiation across the barriers
    should be computed}
  \item{width}{Real number specifying the width around the barrier in
    the computation of its local differentiation}
  \item{plot.pairs.borders}{Logical to tell whether the pairs of
      individuals coming into the computation of the differentiation of the barriers
      should be plotted}
}


\value{A list whose components can be seen using \code{summary}}

\references{}
\author{Gilles Guillot}
\note{ }
\seealso{}
\examples{}

\keyword{}

