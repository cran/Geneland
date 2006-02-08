\name{simFmodel}
\alias{simFmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation of multi-locus genetic data from the spatial F-model}
\description{
  Simulates multi-locus genotypes and spatial coordinates for 
  individuals belonging to some spatially organised populations.
}
\usage{
simFmodel(nindiv, coordinates = NULL, coord.lim, number.nuclei,
coord.nuclei = NULL, color.nuclei = NULL,nall, npop, freq.model="Dirichlet",drift=NULL, seed=NULL, plots = FALSE, ploth = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nindiv}{Integer:  Number of individuals}
  \item{coordinates}{ Matrix (2 rows, \code{nindiv} columns) of spatial
    coordinates of individuals}
  \item{coord.lim}{ Vector of limits of spatial domain to be considered
    (x min, x max, y min, y max) }
  \item{number.nuclei}{ Integer: number of nuclei in the Voronoi tessellation}
  \item{coord.nuclei}{Coordinates of nuclei of Voronoi tessellation}
  \item{color.nuclei}{Population labels of the nuclei (vector of integers of size \code{number.nuclei})}
  \item{nall}{ Vector of integers giving number of alleles at each locus}
  \item{npop}{Number of populations }
  \item{freq.model}{model for frequencies:"Falush" or "Dirichlet" }
  \item{drift}{Vector (of size \code{npop}) of drift factors between 0
  and 1 (only for the Falush model)}
  \item{seed}{Seed to initialize the random number generator}
  \item{plots}{Logical: if TRUE, spatial coordinates are ploted}
  \item{ploth}{Logical: if TRUE, barplots for allele frequencies are ploted}
}
\details{ \code{number.nuclei} uniform i.i.d points are randomly spread on
  the rectangular domain. These points generates the so called Voronoi
  tessellation of the domain in \code{number.nuclei} polygonal sub-domains. Each
  polygon is given a color uniformly on {1, \code{npop}}. The union
  of polygons of the color k gives the domain of population k.
Then \code{nindiv} uniform i.i.d  points are randomly spread on the
domain and stand for the locations of individuals.
Allele frequencies in the ancestral population are sampled from
independent Dirichlet D(1,...,1).
Allele frequencies in the present time population are drawn from
Dirichlet distrubution whose parameters depend on drift factors
\code{drift} and allele frequencies in the ancestral population.
Individual genotypes in each population are drawn from the allele
frequencies of the corresponding population assuming Hardy-Weinberg
equilibrium and linkage equilibrium.}

\value{A list of variables involved in the simulation. The elements of
  this list are: 
coordinates, genotypes, allele.numbers, number.nuclei, coord.nuclei, color.nuclei, 
frequencies,  ancestral.frequencies, drifts, index.nearest.nucleus}

\references{Falush D., M. Stephens, J.K. Pritchard, 
  Inference of population structure using multilocus genotype data:
  Linked loci and correlated allele frequencies, Genetics, pp 1567-1587,
  vol 164, 2003

 G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
  model for landscape genetics. Genetics, 170, 1261-1280, 2005.

  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
  genetics. Molecular Ecology Notes, 5, 712-715, 2005.}

\author{Gilles Guillot  }




\seealso{ Function \code{\link{mcmcFmodel}} for inference}

\keyword{}