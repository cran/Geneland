\name{PosteriorMode}
\alias{PosteriorMode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Map of mode of posterior distribution of population membership}
\description{Plots a map giving  the modal population as a
  color code for each pixel.}

\usage{
function(path.data,path.mcmc,write=FALSE,plotit=TRUE,
                          printit=FALSE,file=NULL,main.title=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path.data}{ Path to input files directory }
  \item{path.mcmc}{ Path to output files directory }
  \item{write}{Logical: if TRUE the modal class of each individual is
    written in a plain ascci file named posterior.mode.txt}
  \item{plotit}{Logical: if TRUE the map is plotted}
  \item{printit}{Logical : if TRUE, figures are also printed}
  \item{file}{Character : Path to file where figures
    should be printed}
  \item{main.title}{Character : Title to appear on top of the graphic}
}

\references{A spatial statistical model for landscape genetics, Guillot,
  Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.

}
\author{Gilles Guillot  }

\keyword{}
