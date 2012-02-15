\name{Fstat}
\alias{Fstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes F statistics on a population genetics dataset given as R
  object}
\description{Computes F statistics according to Weir and Cockerham's 
  estimators.
  Missing values are allowed and accounting for in the computation of
  Fst.
  The presence of missing values involves a downward bias in the
  computation of Fis.
  This function should not be used on haploid data.
}


\usage{
Fstat(genotypes,npop,pop.mbrship,ploidy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{genotypes}{Diploid codominant genotype data. A matrix with one line per
    individual and 2 columns per locus}
  \item{npop}{total number of popluation present in the dataset}
  \item{pop.mbrship}{Vector of integers giving the population membership
    for each individual}
   \item{ploidy}{Integer: 1 or 2 (default is 2) under development. Do
     not use for haploid data.}
  }

\value{A list with components 
  \item{Fis }{A vector of estimations of within-population Fis}
  \item{Fst }{A matrix of estimations the pairwise population
    Fst. (Only the upper triangular part is returned).}
}

\author{Arnaud Estoup for original code in Turbo Pascal.
  Translation in Fortran and interface with R by Gilles Guillot}

\references{ Weir, B.S. and C.C. Cockerham, 
  Estimating F-statistics for the analysis of population structure,
  Evolution, 1984,
  vol. 38, 1358-1370.
}
