\name{mcmcFmodel}
\alias{mcmcFmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inference in a  spatial statistical model  }
\description{ Markov Chain Monte-Carlo inference in the spatial F-model}
\usage{
mcmcFmodel(repdat, repmcmc, lambdamax, dt, nclassmin, nclassinit,
nclassmax, nppmax, nchain, stepw, model, varnclass, spatial)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repdat}{ Path to input files directory }
  \item{repmcmc}{ Path to output files directory }
  \item{lambdamax}{Maximum rate of Poisson process  (real number >0).
  Setting \code{lambdamax} equal to the number of individuals in the
  dataset has proved to be efficient in many cases.}
  \item{dt}{ Parameter prescribing the amount of unctertainty attached
    to spatial coordinates. If \code{dt}=0 spatial coordinates are
    consiered as true coordinates, if \code{dt}>0 it is assumed that observed
  coordinates are true coordinates blurred by an additive noise uniform
  on a square of side \code{dt} centered on 0.}
  \item{nclassmin}{Minimum number of populations (integer >=1) }
  \item{nclassinit}{Initial number of populations
    ( integer sucht that
    \code{nclassmin} =< \code{nclassinit} =< \code{nclassmax})}
  \item{nclassmax}{Maximum number of populations (integer >=
    \code{nclassinit}).
    There is no obvious rule to select \code{nclassmax},
    it should be set to a value larger than any value that
    you can reasonably expect for your data.}
  \item{nppmax}{Integer: Maximum number of nuclei in the
    Poisson-Voronoi tessellation. A good guess consists in setting this
    value equal to \code{10*lambdamax}. The relevance of this rule can be
    checked by inspection of the MCMC run. The number of tiles should not
    go too close from \code{nppmax}. If it does, you should re-run your
    chain  with a larger value for \code{nppmax}}.
  \item{nchain}{Number of MCMC iterations}
  \item{stepw}{Number of MCMC iterations between two writing steps (if \code{stepw}=1, all
    states are saved whereas if e.g. \code{stepw}=10 only each 10 iteration is saved)}
  \item{model}{Character: "Falush" or "Dirichlet" (model for
    frequencies). 
  See also details in detail section of \code{\link{Geneland}} help page.}
  \item{varnclass}{Logical: if TRUE the number of class is treated as
    unknown and will vary along the MCMC inference, if FALSE it will be
    fixed to the initial value \code{nclassinit}. 
    \code{varnclass = TRUE} *should not* be used in conjunction with
    \code{model = "Falush"} as in this case it seems that large numbers
    of populations are not penalized enough and there is a serious risk
    of inferring spurious sub-populations.}
    
  \item{spatial}{Logical: if TRUE the colored Poisson-Voronoi
    tessellation is used as a prior for the spatial organisation of
    populations. If FALSE, all clustering receive equal prior
    probability. In this case spatial information (i.e coordinates)
    are not used  and the locations of  the nuclei are initialized and
    kept fixed at the locations of individuals.}
   
}
\details{ \code{repdat} directory should contain at least three files
  named exactly :   \code{ coord.txt},  \code{genotypes.txt} and
  \code{allele.numbers.txt}   and containing respectively the spatial
  coordinates, the genotypes and the number of alleles per locus.
  
  These files should be plain ascci files with no header.
  \code{coordinates.txt} should have one line per individual
  and two columns (x, y coordinates) missing value are not allowed.
  
  \code{genotypes.txt}  should have one line per individual
  and two columns per locus,
  alleles should be coded as integer between 1 and the maximum number of
  alleles for this locus, missing genotypic values should be coded as -999.
  If your initial file as the required number of lines and columns but
  with alleles coded with non-consecutive integers, then the function
  \code{\link{FormatGenotypes}} should be used first. 
  
  \code{allele.numbers.txt} should have one line per locus an one column
  giving number of alleles for consecutive loci.

  
}


\value{
Successive states of all blocks of parameters are written in files
contained in \code{repout} and named after the type of parameters they contain.
}
\section{Storage format}{
All parameters processed by function \code{\link{mcmcFmodel}} are
written in the directory  specified by \code{repout} as follows:

\item File \code{population.numbers.txt} contains values of the number of
populations (\code{nchain} lines, one line per iteration of the MCMC algorithm)

\item File \code{nuclei.numbers.txt} contains the number of points in the Poisson
point process generating the Voronoi tessellation

\item File \code{color.nuclei.txt} contains vectors of integers of
length \code{nppmax} coding the class membership of each Voronoi tile.
Vectors of class membership for successive states of the chain are
concatenated in one column. Some entries of the vector containing
clas membership for a current state may have missing values as the
actual number of polygon may be smaller that the maximum number allowed
\code{nppmax}. 
This file has \code{nppmax*chain/stepw} lines

\item File \code{coord.nuclei.txt} contains coordinates of points in the Poisson
point process generating the Voronoi tessellation. It has
\code{nppmax*chain/stepw} lines
and two columns (hor. and vert. coordinates).

\item File \code{drifts.txt} contains the drift factors for each
population, (one column per population). 

\item File \code{ancestral.frequencies.txt} contains allele frequencies in ancestral
population. Each line contains all frequencies of the current state.
The file has \code{nchain} lines.
In each line, values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first).

\item File \code{frequencies.txt} contains allele frequencies of present time
populations. Column xx contains frequencies of population numer xx.
In each column values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first), and
values of successive iterations are pasted.
The file has \code{nallmax*nloc*nchain/stepw} lines where \code{nallmax} 
is the maximum numer of alleles over all loci.

\item File \code{Poisson.process.rate.txt} contains rates of Poisson process

\item File \code{log.likelihood.txt} contains log-likelihood of data
for the current state of parameters of the Markov chain.

\item File \code{log.posterior.density.txt} contains log of posterior probability
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



  #  first look for a place to write
  #  simulated data and MCMC outputs
if(.Platform$OS == "unix"){
repdat = "/tmp/"
repmcmc= "/tmp/"
}

if(.Platform$OS == "windows"){
repdat = "/temp/"
repmcmc= "/temp/"
}

 # Simulation of a dataset made of 2 populations
 # 2 loci and 2 alleles per locus
simFmodel(nindiv=100,
          slim=c(0,1,0,1),
          npp=2,
          nloc=5,
          nall=c(5,5,5,5,5),
          nclass=2,
          drift=c(.3,.3),
          plots=FALSE,
          ploth=FALSE,
          seed=12345,
          write=TRUE,
          repout=repdat)


  # First run of MCMC algorithm
  # in order to get the posterior mode of the number of populations
mcmcFmodel(repdat=repdat,
           repmcmc=repmcmc,
           lambdamax=100,
           dt=0,
           nclassmin=1,
           nclassinit=5,
           nclassmax=10,
           nppmax=200,
           nchain=5000,
           stepw=10,
           model="Dirichlet",
           varnclass=TRUE,
           spatial=TRUE)

 # Trace of number of populations
 # Should display a mode at 2
Plotnclass(repdat,repmcmc)


 # Then re-run the chain with fixed number of populations
mcmcFmodel(repdat=repdat,
           repmcmc=repmcmc,
           lambdamax=100,
           dt=0,
           nclassmin=1,
           nclassinit=2,
           nclassmax=2,
           nppmax=200,
           nchain=5000,
           stepw=10,
           model="Dirichlet",
           varnclass=FALSE,
           spatial=TRUE)

  # Post-processing the chain 
PostProcessChain(repdat=repdat,
                 repmcmc=repmcmc,
                 nxdom=50,
                 nydom=50,
                 burnin=0)

  # Plots allele frequencies of allele #1 at locus  #1 in sub-population  #1
PlotFreq(repdat=repdat,repmcmc=repmcmc,ipop=1,iloc=1,iall=1)

  # Map of posterior probabilites
  # of population membership
PlotTessellation(repdat=repdat,repmcmc=repmcmc)


  # Map of posterior mode 
  # of population membership
PosteriorMode(repdat=repdat,repmcmc=repmcmc,
              write=FALSE,plotit=TRUE)




}

\keyword{}
