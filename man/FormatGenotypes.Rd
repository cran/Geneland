\name{FormatGenotypes}
\alias{FormatGenotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Formatting file of genotype data.}
\description{Takes genotype data as a matrix with
  one line per individual and two columns per locus,
  with alleles coded by integers (number of replications for
  micro-satellites data).
  Build a new matrix with alleles acodes as consecutive integer.
  If allele j has 7 alleles they will be coded as 1,2,...7.
}
\usage{
FormatGenotypes(genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{A matrix with
  one line per individual and two columns per locus,
  with alleles coded by integers}
}

\value{
  A matrix with
  one line per individual and two columns per locus,
  with alleles coded by integers}


\author{Gilles Guillot}


\keyword{}