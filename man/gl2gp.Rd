\name{gl2gp}
\alias{gl2gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Geneland to Genepop conversion}
\description{Takes the matrices  of  coordinates and genotypes in the Geneland format
  and writes it as an ascci file in the Genepop format}
\usage{
gl2gp(coordinates, genotypes, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{ Matrix of coordinates}
  \item{genotypes}{MAtrix of genotypes}
  \item{file}{ Character string giving the path to the file
    where the data in Genepop format should be written}
}

\value{An ascci file is written. This file may need to be
process in order to convert eol character from unix style to dos style. }

\author{Gilles Guillot and Tanguy de l'Argentaye}


