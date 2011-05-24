\name{EstimateFreqNA}
\alias{EstimateFreqNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate frequencies of null alleles.}
\description{Estimate frequencies of null alleles at each locus.}


\usage{
EstimateFreqNA(path.mcmc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path.mcmc}{path to directory containing stored MCMC computations}
}

\value{A vector of length \code{nloc} (the number of loci) whose entries
  are estimated frequencies of null alleles.
}

\author{Gilles Guillot}

\references{  Gilles Guillot, Filipe Santos  and Arnaud Estoup,
  Analysing georeferenced population genetics data with
  Geneland: a new algorithm to deal with null alleles and a friendly graphical user interface
Bioinformatics 2008 24(11):1406-1407.}
