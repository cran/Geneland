\name{PlotFreqA}
\alias{PlotFreqA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots allele frequencies in ancestral population along the MCMC run
  and optionally prints it}

\description{Plot frequency of  allele \code{iall}  of locus \code{iloc}
  in ancestral population}

\usage{
PlotFreqA(repdat, repmcmc,iloc,iall,printit=FALSE,path=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{repdat}{ Path to input files directory }
   \item{repmcmc}{ Path to output files directory }
   \item{iloc}{Integer number : index of locus}
   \item{iall}{Integer number : index of allele}
   \item{printit}{Logical : if TRUE, figures are also printed}
   \item{path}{Character : Path to directory where figures
    should be printed}
 }


\references{A spatial statistical model for landscape genetics, Guillot,
  Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }



\keyword{}