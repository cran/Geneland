\name{nullify}
\alias{nullify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulates null alleles}
\description{Simulates null alleles}
\usage{
nullify(genotypes, nall.null = 1, nloc.null)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{ a matrix of genotypes as produced by simFmodel and simIBD}
  \item{nall.null}{number of null alleles on each locus}
  \item{nloc.null}{number of loci carrying null alleles}
}
\details{For diploid data only}

\value{A list with component: 
  \item{genotypes }{The new genotypes after alteration.}
}

\author{Gilles Guillot}

