\name{simFmodel}
\alias{simFmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation of multi-locus genetic data from the spatial F-model}
\description{
  Simulates multi-locus genotypes and spatial coordinates for 
  individuals belonging to some spatially organised populations.
}
\usage{
simFmodel(nindiv, s = NULL, slim, npp, u = NULL, c = NULL,
nloc, nall, nclass, drift, seed=NULL, plots = TRUE, ploth = TRUE, write = FALSE, repout=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nindiv}{Integer:  Number of individuals}
  \item{s}{ Matrix (2 rows, \code{nindiv} columns) of spatial
    coordinates of individuals}
  \item{slim}{ Vector of limits of spatial domain to be considered (abs
    min, abs max, ord min, ord max) }
  \item{npp}{ Integer: number of nuclei in the Voronoi tessellation}
  \item{u}{Coordinates of nuclei of Voronoi tessellation}
  \item{c}{Population labels of the nuclei (vector of integers of size \code{npp})}
  \item{nloc}{ number of loci}
  \item{nall}{ Vector of integers (of size \code{nloc}) giving number of alleles at each locus}
  \item{nclass}{Number of populations }
  \item{drift}{Vector (of size \code{nclass}) of drift factors between 0
  and 1}
  \item{seed}{Seed to initialize the random number generator}
  \item{plots}{Logical: if TRUE, spatial coordinates are ploted}
  \item{ploth}{Logical: if TRUE, barplots for allele frequencies are ploted}
  \item{write}{Logical: if TRUE, data are written to directory \code{repout}}
  \item{repout}{Path to output directory}
}
\details{ \code{npp} uniform i.i.d points are randomly spread on
  the rectangular domain. These points generates the so called Voronoi
  tessellation of the domain in \code{npp} polygonal sub-domains. Each
  polygon is given a color uniformly on {1, \code{nclass}}. The union
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
equilibrium and linkage equilibrium. }

\value{Simulated values are written in files in \code{repout} directory}

\references{Falush D., M. Stephens, J.K. Pritchard, 
  Inference of population structure using multilocus genotype data:
  Linked loci and correlated allele frequencies, Genetics, pp 1567-1587,
  vol 164, 2003

  A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }




\seealso{ Function \code{\link{mcmcFmodel}} for inference}

\keyword{}