\name{show.simdata}
\alias{show.simdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical display of data simulated by simdata}
\description{
  Graphical display of data simulated by simdata}

\usage{
show.simdata(dataset,
plot.coord = FALSE,
file.plot.coord,
plot.tess = FALSE,
file.plot.tess,
plot.freq.grid = FALSE,
file.plot.freq.grid,
loc.grid = 1,
plot.freq.indiv = FALSE,
file.plot.freq.indiv,
loc.indiv=1,
zlim.freq=c(0,1),
plot.gen = FALSE,
file.plot.gen)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{An R object produced by function \code{simdata}}
  \item{plot.coord}{Logical}
  \item{file.plot.coord}{Character string. Path to the file where the
    graphic should be stored}
  \item{plot.tess}{ Logical}
  \item{file.plot.tess}{Character string. Path to the file where the
    graphic should be stored}
  \item{plot.freq.grid}{Logical}
  \item{file.plot.freq.grid}{Character string. Path to the file where the
    graphic should be stored}
  \item{loc.grid}{Vector of integers giving indices of loci for which the
    map of frequencies is required}
  \item{plot.freq.indiv}{Logical}
  \item{file.plot.freq.indiv}{Character string. Path to the file where the
    graphic should be stored}
  \item{loc.indiv}{Vector of integers giving indices of loci for which the
    plot of frequencies at sites of individuals is required}
  \item{plot.gen}{Logical}
  \item{file.plot.gen}{Character string. Path to the file where the
    graphic should be stored}
  \item{zlim.freq}{A vector of two integers giving the limit of the
    values for the image of fields of frequencies}
}

\value{
No value returned
}
\references{Guillot G. and Estoup A. }
\author{Gilles Guillot}

