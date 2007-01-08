\name{show.simIBD}
\alias{show.simIBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical display of data simulated by simIBD}
\description{
  Graphical display of data simulated by simIBD}
}
\usage{
show.simIBD(dataset,
plot.coord = FALSE,
plot.tess = FALSE,
plot.freq.grid = FALSE,
plot.freq.indiv = FALSE,
loc.grid = 1:nloc,
all.grid = dataset$allele.numbers,
zlim.freq=c(0,1),
plot.gen = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{An R object produced by function \code{sim.IBD}}
  \item{plot.coord}{Logical}
  \item{plot.tess}{ Logical}
  \item{plot.freq.grid}{Logical}
  \item{plot.freq.indiv}{Logical}
  \item{loc.grid}{Vector of integers giving indices of loci for which the
    field of frequencies is required}
  \item{all.grid}{Vector of integesr giving indices of alleles for which the
    field of frequencies is required}
  \item{plot.gen}{Logical}
  \item{zlim.freq}{A vector of two integers giving the limit of the
    values for the image of fields of frequencies}
}
\details{

}
\value{
No value returned
}
\references{Guillot G. and Estoup A. }
\author{Gilles Guillot}
\note{}
\examples{}
\keyword{}

