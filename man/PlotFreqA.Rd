\name{PlotFreqA}
\alias{PlotFreqA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots allele frequencies in ancestral population along the MCMC run
  and optionally prints it}

\description{Plot frequency of  allele \code{iall}  of locus \code{iloc}
  in ancestral population}

\usage{
PlotFreqA(path.mcmc,iloc,iall,printit=FALSE,path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{path.mcmc}{ Path to output files directory }
   \item{iloc}{Integer number : index of locus}
   \item{iall}{Integer number : index of allele}
   \item{printit}{Logical : if TRUE, figures are also printed}
   \item{path}{Character : Path to directory where figures
    should be printed}
 }


\author{Gilles Guillot  }


