\name{Fstat}
\alias{Fstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes F statistics on a population genetics dataset given as R
  object}
\description{Computes F statistics according to Weir and Cockerham's 
  estimators.
  Missing values are allowed but as of today, the NA code is only treated as an extra
  allele which alters the result.}

\usage{
Fstat(genotypes,npop,pop.mbrship)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{genotypes}{Genotypes of individuals. A matrix with one line per
    individual and 2 columns per locus}
  \item{npop}{total number of popluation present in the dataset}
  \item{pop.mbrship}{Vector of integers giving the population membership
     for each individual  }
}
\value{A list with components 
  \item{Pairwise.Fis }{A vector of estimations of within-population Fis}
  \item{Pairwise.Fst }{A matrix of estimations the pairwise population
    Fst. (Only the upper triangular part is returned).}
}

\author{Arnaud Estoup for original code in Turbo Pascal.

  Translation in Fortran and interface with R by Gilles Guillot}

\references{ Weir, B.S. and C.C. Cockerham, 
  Estimating F-statistics for the analysis of population structure,
  Evolution, 1984,
  vol. 38, 1358-1370.
}

\keyword{}
