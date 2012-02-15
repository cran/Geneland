\name{PlotFreq}
\alias{PlotFreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the trace of allele frequencies in present time population
  and optionnaly print it.}

\description{Plot of allele frequencies in present time populations number
  \code{ipop},
  for allele \code{iall}  of locus \code{iloc}}

\usage{
PlotFreq(path.mcmc,ipop,iloc,iall,printit=FALSE,path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{path.mcmc}{ Path to output files directory }
   \item{ipop}{Integer number : index of population}
   \item{iloc}{Integer number : index of locus}
   \item{iall}{Integer number : index of allele. If \code{MCMC}
     was launched with option \code{filter.null.alleles=TRUE}, an extra
     fictive allele standing for putative null alleles is created. It
     estimated frequency can be also plotted. If there is say, 5 alleles
   at a locus, the estimated frequency of null alleles can be seen
   invoking \code{PlotFreq} with \code{iall=6}.}
   \item{printit}{Logical : if TRUE, figures are also printed}
   \item{path}{Character : Path to directory where figures
     should be printed}
  }


\author{Gilles Guillot  }



