\name{mcmcFmodel}
\alias{mcmcFmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inference in a  spatial statistical model  }
\description{ Markov Chain Monte-Carlo inference in the spatial F-model}
\usage{
mcmcFmodel(coordinates,genotypes,allele.numbers,
path.mcmc, rate.max, delta.coord, npopmin, npopinit,
npopmax, nb.nuclei.max, nit, thinning, freq.model, varnpop, spatial)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.}
  \item{genotypes}{Genotypes of individuals. A matrix with one line per
    individual and 2 columns per locus}
  \item{allele.numbers}{A vector of integer containing the number of possible
    allele for each locus}
  \item{path.mcmc}{ Path to output files directory. It seems that the
    path should be given in the Unix style even under Windows (use \/
    instead of \\).
    This path *have to*  end with a slash (\/)
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
    go too close from \code{nb.nuclei.max}. If it does, you should re-run your
    chain  with a larger value for \code{nb.nuclei.max}.}
  \item{nit}{Number of MCMC iterations}
  \item{thinning}{Number of MCMC iterations between two writing steps (if \code{thinning}=1, all
    states are saved whereas if e.g. \code{thinning}=10 only each 10 iteration is saved)}
  \item{freq.model}{Character: "Falush" or "Dirichlet" (model for
    frequencies). 
  See also details in detail section of \code{\link{Geneland}} help page.}
  \item{varnpop}{Logical: if TRUE the number of class is treated as
    unknown and will vary along the MCMC inference, if FALSE it will be
    fixed to the initial value \code{npopinit}. 
    \code{varnpop = TRUE} *should not* be used in conjunction with
    \code{freq.model = "Falush"} as in this case it seems that large numbers
    of populations are not penalized enough and there is a serious risk
    of inferring spurious sub-populations.}
    
  \item{spatial}{Logical: if TRUE the colored Poisson-Voronoi
    tessellation is used as a prior for the spatial organisation of
    populations. If FALSE, all clustering receive equal prior
    probability. In this case spatial information (i.e coordinates)
    are not used  and the locations of  the nuclei are initialized and
    kept fixed at the locations of individuals.}
   
}


\value{
Successive states of all blocks of parameters are written in files
contained in \code{path.mcmc} and named after the type of parameters they contain.
}
\section{Storage format}{
All parameters processed by function \code{\link{mcmcFmodel}} are
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
is the maximum numer of alleles over all loci.

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

\references{A spatial statistical model for landscape genetics, Guillot,
Estoup, Mortier, Cosson, Genetics, 2005
  
 Guillot, G., Geneland : A program for landscape genetics.
 Molecular Ecology  Notes, submited.}
\author{Gilles Guillot  }



\seealso{ \code{\link{simFmodel}}}

\examples{
  # Below is a complete sequence 
  # of commands using Geneland functions

  # we assume that Geneland is installed
  # and loaded by library("Geneland")



  #  define a place for MCMC outputs
path.mcmc = paste(tempdir(),"/",sep="")

 # Simulation of a dataset made of 2 populations
 # 2 loci and 2 alleles per locus
sim = simFmodel(nindiv=100,
          coord.lim=c(0,1,0,1),
          number.nuclei=2,
          nloc=5,
          nall=c(5,5,5,5,5),
          npop=2,
          drift=c(.3,.3),
          plots=FALSE,
          ploth=FALSE,
          seed=123)


  # First run of MCMC algorithm
  # in order to get the posterior mode of the number of populations
mcmcFmodel(sim$coordinates,sim$genotypes,sim$allele.numbers,
           path.mcmc=path.mcmc,
           rate.max=100,
           delta.coord=0,
           npopmin=1,
           npopinit=5,
           npopmax=10,
           nb.nuclei.max=200,
           nit=10000,
           thinning=10,
           freq.model="Dirichlet",
           varnpop=TRUE,
           spatial=TRUE)

  # Trace of number of populations
  # Should display a mode at 2
Plotnpop(path.mcmc)


 # Then re-run the chain with fixed number of populations
mcmcFmodel(sim$coordinates,sim$genotypes,sim$allele.numbers,
           path.mcmc=path.mcmc,
           rate.max=100,
           delta.coord=0,
           npopmin=1,
           npopinit=2,
           npopmax=2,
           nb.nuclei.max=200,
           nit=5000,
           thinning=10,
           freq.model="Dirichlet",
           varnpop=FALSE,
           spatial=TRUE)

   # Post-processing the chain 
PostProcessChain(sim$coordinates,sim$genotypes,sim$allele.numbers,
                  path.mcmc=path.mcmc,
                  nxdom=50,
                  nydom=50,
                  burnin=0)

   # Plots allele frequencies of allele #1 at locus  #1 in sub-population  #1
 PlotFreq(sim$allele.numbers,
path.mcmc=path.mcmc,ipop=1,iloc=1,iall=1)

   # Map of posterior probabilites
   # of population membership
 PlotTessellation(sim$coordinates,path.mcmc=path.mcmc)


   # Map of posterior mode 
   # of population membership
 PosteriorMode(sim$coordinates,path.mcmc=path.mcmc,
               write=FALSE,plotit=TRUE)




}

\keyword{}
