\name{PlotFreq}
\alias{PlotFreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the trace of allele frequencies in present time population
  and optionnaly print it.}

\description{Plot of allele frequencies in present time populations number
  \code{ipop},
  for allele \code{iall}  of locus \code{iloc}}

\usage{
PlotFreq(repdat, repmcmc,ipop,iloc,iall,printit=FALSE,path=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{repdat}{ Path to input files directory }
   \item{repmcmc}{ Path to output files directory }
   \item{ipop} {Integer number : index of population}
   \item{iloc}{Integer number : index of locus}
   \item{iall}{Integer number : index of allele}
   \item{printit}{Logical : if TRUE, figures are also printed}
   \item{path}{Character : Path to directory where figures
     should be printed}
  }




\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005}
\author{Gilles Guillot  }

\keyword{}




