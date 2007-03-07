\name{Fstat.output}
\alias{Fstat.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes F statistics on the output of an inference by MCMC simulation}
\description{Computes F statistics according to Weir and Cockerham's 
  estimators.    Missing values are allowed but as of today, the NA code is only treated as an extra
  allele which might bias the result.}

\usage{
Fstat.output(genotypes,path.mcmc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{genotypes}{Genotypes of individuals. A matrix with one line per
    individual and 2 columns per locus}
  \item{path.mcmc}{ Path to output files directory }
}
}
\value{A list with components 
  \item{Pairwise.Fis }{A matrix of real numbers estimating the pairwise Fis}
  \item{Pairwise.Fst }{A matrix of real numbers estimating the pairwise Fst}
}

\author{Arnaud Estoup for original code in Turbo Pascal.

  Translation in Fortran and interface with R by Gilles Guillot}

\references{ Weir, B.S. and C.C. Cockerham, 
  Estimating F-statistics for the analysis of population structure,
  Evolution, 1984,
  vol. 38, 1358-1370.
}

\keyword{}
