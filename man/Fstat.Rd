\name{Fstat}
\alias{Fstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{F statistics}
\description{Computes F statistics according to Weir and Cockerham's 
  estimators.}

\usage{
Fstat(genotypes,allele.numbers, path.mcmc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{genotypes}{Genotypes of individuals. A matrix with one line per
    individual and 2 columns per locus}
  \item{allele.numbers}{A vector of integer containing the number of possible
    allele for each locus}
  \item{path.mcmc}{ Path to output files directory }
}
}
\value{A list with components 
  \item{Total.Fit }{A real number estimating the total Fit}
  \item{Pairwise.Fis }{A matrix of real numbers estimating the pairwise Fis}
  \item{Pairwise.Fst }{A matrix of real numbers estimating the pairwise Fst}
  \item{Pairwise.Fit }{A matrix of real numbers estimating the pairwise Fit}
}

\author{Arnaud Estoup for original code in Turbo Pascal.

  Translation in Fortran and interface with R by Gilles Guillot}

\references{ Weir, B.S. and C.C. Cockerham, 
  Estimating F-statistics for the analysis of population structure,
  Evolution, 1984,
  vol. 38, 1358-1370.
}

\keyword{}