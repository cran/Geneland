\name{MCMC}
\alias{MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov Chain Monte-Carlo inference under various models}
\description{ Markov Chain Monte-Carlo inference under various models}
\usage{
MCMC(
      coordinates=NULL,
      genotypes,## 
                 ploidy=2,
                 dominance="Codominant",
                 allele.numbers,
                                        # path to output directory
                 path.mcmc,
                                        # hyper-prior parameters
                 rate.max,delta.coord=0,shape1=2,shape2=20,
                 npopmin=1,npopinit,npopmax,
                                        # dimensions
                 nb.nuclei.max,
                                        # options in mcmc computations
                 nit,thinning=1,freq.model="Correlated",
                 varnpop=TRUE,
                 spatial=TRUE,
                 jcf=TRUE,
                 filter.null.alleles=TRUE,
                 prop.update.cell=0.1,
                 write.rate.Poisson.process=FALSE,
                 write.number.nuclei=TRUE,
                 write.number.pop=TRUE,
                 write.coord.nuclei=TRUE,
                 write.color.nuclei=TRUE,
                 write.freq=TRUE,
                 write.ancestral.freq=TRUE,
                 write.drifts=TRUE,
                 write.logposterior=TRUE,
                 write.loglikelihood=TRUE,
                 write.true.coord=TRUE,
                 write.size.pop=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.}
  \item{genotypes}{Genotypes of individuals. A matrix with dimension as
    described below.


    For diploid organisms at codominant markers: a matrix with one line per
    individual and two  columns per locus.

    
    For diploid organisms at dominant markers: a matrix with one line per
    individual and one column per locus. Presence/absence of a band should be
    coded as 0/1 (0 for absence / 1 for presence).


    For haploid organisms: a matrix with one line per
    individual and one  column per locus.}
    

  \item{ploidy}{An integer equal to 1 for haploid data and 2 for diploid
    data. }
  \item{dominance}{A character string that should be either "Codominant"
    (default) for codominant markers such as microsatellites and SNPs
    or "Dominant" for AFLPs.}
  \item{allele.numbers}{Vector of integers giving the number of observed
  alleles at each locus. Only for debugging. Should be left unspecified.}
  \item{path.mcmc}{ Path to output files directory. It seems that the
    path should be given in the Unix style even under Windows (use \/
    instead of \\).
    This path *has to*  end with a slash (\/)
    (e.g. path.mcmc="/home/me/Geneland-stuffs/")
  }
  \item{rate.max}{Maximum rate of Poisson process  (real number >0).
  Setting \code{rate.max} equal to the number of individuals in the
  dataset has proved to be efficient in many cases.}
  \item{delta.coord}{ Parameter prescribing the amount of unctertainty attached
    to spatial coordinates. If \code{delta.coord}=0 spatial coordinates are
    consiered as true coordinates, if \code{delta.coord}>0 it is assumed that observed
  coordinates are true coordinates blurred by an additive noise uniform
  on a square of side \code{delta.coord} centered on 0.}
\item{shape1}{First parameter in the Beta(shape1,shape2) prior
  distribution of the drift parameters in the Correlated model.}
\item{shape2}{Second parameter in the Beta(shape1,shape2) prior
  distribution of the drift parameters in the Correlated model.}
  \item{npopmin}{Minimum number of populations (integer >=1) }
  \item{npopinit}{Initial number of populations
    ( integer sucht that
    \code{npopmin} =< \code{npopinit} =< \code{npopmax})}
  \item{npopmax}{Maximum number of populations (integer >=
    \code{npopinit}).
    There is no obvious rule to select \code{npopmax},
    it should be set to a value larger than any value that
    you can reasonably expect for your data.}
  \item{nb.nuclei.max}{Integer: Maximum number of nuclei in the
    Poisson-Voronoi tessellation. A good guess consists in setting this
    value equal to \code{3*rate.max}. Lower values
    can also be used in order to speed up computations. The relevance of
    the value set can be
    checked by inspection of the MCMC run. The number of tiles should not
    go too close to \code{nb.nuclei.max}. If it does, you should re-run your
    chain  with a larger value for \code{nb.nuclei.max}. In case of use
    of the option \code{SPATIAL=FALSE}, \code{nb.nuclei.max} should be
    set equal to the number of individuals.}
  \item{nit}{Number of MCMC iterations}
  \item{thinning}{Number of MCMC iterations between two writing steps (if \code{thinning}=1, all
    states are saved whereas if e.g. \code{thinning}=10 only each 10 iteration is saved)}
  \item{freq.model}{Character: "Correlated" or "Uncorrelated" (model for
    frequencies). 
  See also details in detail section of \code{\link{Geneland}} help page.}
  \item{varnpop}{Logical: if TRUE the number of class is treated as
    unknown and will vary along the MCMC inference, if FALSE it will be
    fixed to the initial value \code{npopinit}. 
    \code{varnpop = TRUE} *should not* be used in conjunction with
    \code{freq.model = "Correlated"} as in this case it seems that large numbers
    of populations are not penalized enough and there is a serious risk
    of inferring spurious sub-populations.}
    
  \item{spatial}{Logical: if TRUE the colored Poisson-Voronoi
    tessellation is used as a prior for the spatial organisation of
    populations. If FALSE, all clustering receive equal prior
    probability. In this case spatial information (i.e coordinates)
    are not used  and the locations of  the nuclei are initialized and
    kept fixed at the locations of individuals.}
  \item{jcf}{Logical: if true update of c and f are performed jointly}
  \item{filter.null.alleles}{Logical: if TRUE, tries to filter out null
    alleles. An extra fictive null allele is created at each locus coding
  for all putative null allele. Its frequency  is estimated and can be
  viewed with function \code{PlotFreq}. This option is available only
  with \code{freq.model="Uncorrelated"}.}
\item{prop.update.cell}{Integer between 0 and 1. Proportion of cell updated. For
  debugging only.}
\item{write.rate.Poisson.process}{Logical: if TRUE (default) write rate
  of Poisson process simulated by MCMC}
\item{write.number.nuclei}{Logical: if TRUE (default) write number of nuclei simulated by MCMC }
\item{write.number.pop}{Logical: if TRUE (default) write number of
  populations simulated by MCMC }
\item{write.coord.nuclei}{Logical: if TRUE (default) write coordinates
  of nuclei simulated by MCMC  }
\item{write.color.nuclei}{Logical: if TRUE (default) write color of
  nuclei simulated by MCMC  }
\item{write.freq}{Logical: if TRUE (default is FALSE) write allele
  frequencies simulated by MCMC  }
\item{write.ancestral.freq}{Logical: if TRUE (default is FALSE) write
  ancestral allele frequencies simulated by MCMC }
\item{write.drifts}{Logical: if TRUE (default is FALSE) write drifts simulated by MCMC }
\item{write.logposterior}{Logical: if TRUE (default is FALSE) write
  logposterior simulated by MCMC }
\item{write.loglikelihood}{Logical: if TRUE (default is FALSE) write
  loglikelihood simulated by MCMC }
\item{write.true.coord}{Logical: if TRUE (default is FALSE) write true
  spatial coordinates simulated by MCMC }
\item{write.size.pop}{Logical: if TRUE (default is FALSE) write size of
  populations simulated by MCMC }
}


\value{
Successive states of all blocks of parameters are written in files
contained in \code{path.mcmc} and named after the type of parameters they contain.
}
\section{Storage format}{
All parameters processed by function \code{\link{MCMC}} are
written in the directory  specified by \file{path.mcmc} as follows:

\item File \file{population.numbers.txt} contains values of the number of
populations (\code{nit} lines, one line per iteration of the MCMC algorithm)

\item File \file{nuclei.numbers.txt} contains the number of points in the Poisson
point process generating the Voronoi tessellation

\item File \file{color.nuclei.txt} contains vectors of integers of
length \code{nb.nuclei.max} coding the class membership of each Voronoi tile.
Vectors of class membership for successive states of the chain are
concatenated in one column. Some entries of the vector containing
clas membership for a current state may have missing values as the
actual number of polygon may be smaller that the maximum number allowed
\code{nb.nuclei.max}. 
This file has \code{nb.nuclei.max*chain/thinning} lines

\item File \file{coord.nuclei.txt} contains coordinates of points in the Poisson
point process generating the Voronoi tessellation. It has
\code{nb.nuclei.max*chain/thinning} lines
and two columns (hor. and vert. coordinates).

\item File \file{drifts.txt} contains the drift factors for each
population, (one column per population). 

\item File \file{ancestral.frequencies.txt} contains allele frequencies in ancestral
population. Each line contains all frequencies of the current state.
The file has \code{nit} lines.
In each line, values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first).

\item File \file{frequencies.txt} contains allele frequencies of present time
populations. Column xx contains frequencies of population numer xx.
In each column values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first), and
values of successive iterations are pasted.
The file has \code{nallmax*nloc*nit/thinning} lines where \code{nallmax} 
is the maximum number of alleles over all loci.

\item File \file{Poisson.process.rate.txt} contains rates of Poisson
process

\item File \file{hidden.coord.txt} contains the coordinates of each
individual as updated along the chain if those given as input are not
considered as exact coordinates (which is specified by 
\code{delta.coord} to a non zero value).

\item File \file{log.likelihood.txt} contains log-likelihood of data
for the current state of parameters of the Markov chain.

\item File \file{log.posterior.density.txt} contains log of posterior probability
(up to marginal density of data) of the
current state of parameters in the Markov chain.


}

\references{ G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
  model for landscape genetics. Genetics, 170, 1261-1280, 2005.

  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
  genetics. Molecular Ecology Notes, 5, 712-715, 2005.

   Gilles Guillot, Filipe Santos  and Arnaud Estoup,
  Analysing georeferenced population genetics data with
  Geneland: a new algorithm to deal with null alleles and a friendly graphical user interface
Bioinformatics 2008 24(11):1406-1407.}


\author{Gilles Guillot  }



\seealso{ \code{\link{simFmodel}}}

\examples{


}

\keyword{}
