\name{PosteriorMode}
\alias{PosteriorMode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Map of mode of posterior distribution of population membership}
\description{Plots a map giving  the modal population as a
  color code for each pixel.}

\usage{
PosteriorMode(coordinates,path.mcmc,plotit=TRUE,format="pdf",new.dev=TRUE,
                          printit=FALSE,file,main.title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.}
  \item{path.mcmc}{ Path to output files directory }
  \item{plotit}{Logical: if TRUE the map is plotted}
  \item{printit}{Logical : if TRUE, figures are also printed}
  \item{format}{\code{"ps"} or \code{"pdf"}}
  \item{file}{Character : Path to file where figures
    should be printed}
  \item{main.title}{Character : Title to appear on top of the graphic}
   \item{new.dev}{Logical. Open a new graphical device if TRUE  and \code{printit}
    is FALSE}
}

\references{G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
  model for landscape genetics. Genetics, 170, 1261-1280, 2005.

  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
  genetics. Molecular Ecology Notes, 5, 712-715, 2005.}


\author{Gilles Guillot  }
