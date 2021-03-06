\name{MCMC}
\alias{MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov Chain Monte-Carlo inference of clusters from genotpype data}
\description{ Markov Chain Monte-Carlo inference of clusters from genotpype data}
\usage{
          MCMC(
 ## input data
                 coordinates=NULL, # spatial coordinates
                 geno.dip.codom=NULL, # diploid codominant markers
                                      # one line per indiv.
                                      # two column per marker
                 geno.dip.dom=NULL, # diploid dominant markers
                                    # one line per indiv.
                                    # one column per marker
                 geno.hap=NULL, # haploid
                                # one line per indiv.
                                # one column per marker
                 qtc, # quantitative continuous variables
                 qtd, # quantitative discrete variables
                 ql, # qualitative variables
                 ## path to output directory
                 path.mcmc,
                 ## hyper-prior parameters
                 rate.max,delta.coord=0,shape1=2,shape2=20,
                 npopmin=1,npopinit,npopmax,
                 ## dimensions
                 nb.nuclei.max,
                 ## mcmc computations options
                 nit,
                 thinning=1,
                 freq.model="Uncorrelated",
                 varnpop=TRUE,
                 spatial=TRUE,
                 jcf=TRUE,
                 filter.null.alleles=TRUE,
                 prop.update.cell=0.1,
                 ## writing mcmc output files options 
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
                 write.size.pop=FALSE,
                 write.mean.quanti=TRUE,
                 write.sd.quanti=TRUE,
                 write.betaqtc=FALSE,
                 miss.loc=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{coordinates}{Spatial coordinates of individuals. A matrix with 2
    columns and one line per individual.
    }

  \item{geno.dip.codom}{Genotypes  for diploid data with codominant markers.
    A matrix with one line per individual and two  columns per locus.
  Note that the object has to be of type matrix not table. This can be
    forced by function \code{as.matrix}.}
    
 \item{geno.dip.dom}{Genotypes for diploid data with dominant
   markers. A matrix with one line per individual and one  column per
   locus. Presence/absence of a band should be
    coded as 0/1 (0 for absence / 1 for presence). Dominant and
   codominant markers can be analyzed jointly by passing variables to arguments 
   geno.dip.codom and geno.dip.dom.
 Haploid data and diploid dominant data can not be analyzed jointly in
 the current version.
 Note that the object has to be of type matrix not table. This can be
    forced by function \code{as.matrix}.}

 \item{geno.hap}{Genotypes of haploid data.
   A matrix with one line per individual and one  column per
   locus.     Dominant diploid data and haploid data 
   can be analyzed jointly (e.g. to analyse microsatelite data or SNP
   data together with mtDNA.
 Haploid data and diploid dominant data can not be analyzed jointly in
   the current version.
 Note that the object has to be of type matrix not table. This can be
    forced by function \code{as.matrix}.}
    

  \item{qtc}{A matrix of continuous quantitative  phenotypic variables. One line per individual and one  column per
   phenotypic variable. 
 Note that the object has to be of type matrix not table. This can be
    forced by function \code{as.matrix}. }

  \item{qtd}{A matrix of discrete quantitative phenotypic variables. NOT IMPLEMENTED YET}

   \item{ql}{A matrix of categorical phenotypic variables. NOT IMPLEMENTED YET}

 
\item{path.mcmc}{ Path to output files directory. It seems that the
    path should be given in the Unix style even under Windows (use \/
    instead of \\).
    This path *has to*  end with a slash (\/)
    (e.g. path.mcmc="/home/me/Geneland-stuffs/")
  }

  \item{rate.max}{Maximum rate of Poisson process  (real number >0).
  Setting \code{rate.max} equal to the number of individuals in the
  dataset has proved to be efficient in many cases.}

\item{delta.coord}{ Parameter prescribing the amount of uncertainty attached
    to spatial coordinates. If \code{delta.coord}=0 spatial coordinates are
    considered as true coordinates, if \code{delta.coord}>0 it is assumed that observed
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

\item{write.mean.quanti}{Logical: if TRUE (default is FALSE) write
    means of quantitatives variables in the various groups simulated by MCMC}

\item{write.sd.quanti}{Logical: if TRUE (default is FALSE) write
    standard deviations of quantitatives variables in the various groups simulated by MCMC}

  \item{write.betaqtc}{Logical: if TRUE (default is FALSE) write
    hyper-parameter beta of distribution of quantitatives variables simulated by MCMC}

\item{miss.loc}{A matrix with \code{nindiv} lines and \code{nloc}
  columns of 0 or 1. For each individual, at each locus it says if the
  locus is genuinely missing (no attempt to measure it). This info is
    used under the option \code{filterNA=TRUE} do decide how a double
    missing value should be treated (genuine missing data
    or double null allele).}
}



\value{
Successive states of all blocks of parameters are written in files
contained in \code{path.mcmc} and named after the type of parameters they contain.
}

\section{Storage format}{
All parameters processed by function \code{\link{MCMC}} are
written in the directory  specified by \file{path.mcmc} as follows:
\itemize{
\item File \file{population.numbers.txt} contains values of the number of
populations (\code{nit} lines, one line per iteration of the MCMC
algorithm).


\item File \file{population.numbers.txt} contains values of the number of
populations (\code{nit} lines, one line per iteration of the MCMC algorithm).

\item  File \file{nuclei.numbers.txt} contains the number of points in the Poisson
point process generating the Voronoi tessellation.

\item  File \file{color.nuclei.txt} contains vectors of integers of
length \code{nb.nuclei.max} coding the class membership of each Voronoi tile.
Vectors of class membership for successive states of the chain are
concatenated in one column. Some entries of the vector containing
clas membership for a current state may have missing values as the
actual number of polygon may be smaller that the maximum number allowed
\code{nb.nuclei.max}. This file has \code{nb.nuclei.max*chain/thinning} lines.

\item  File \file{coord.nuclei.txt} contains coordinates of points in the Poisson
point process generating the Voronoi tessellation. It has
\code{nb.nuclei.max*chain/thinning} lines
and two columns (hor. and vert. coordinates).

\item  File \file{drifts.txt} contains the drift factors for each
population, (one column per population).

\item  File \file{ancestral.frequencies.txt} contains allele frequencies in ancestral
population. Each line contains all frequencies of the current state.
The file has \code{nit} lines.
In each line, values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first).

\item  File \file{frequencies.txt}contains allele frequencies of present time
populations. Column xx contains frequencies of population numer xx.
In each column values of allele frequencies are stored by increasing
allele index and and locus index (allele index varying first), and
values of successive iterations are pasted.
The file has \code{nallmax*nloc*nit/thinning} lines where \code{nallmax} 
is the maximum number of alleles over all loci.

\item  File \file{Poisson.process.rate.txt} contains rates of Poisson
process.

\item  File \file{hidden.coord.txt} contains the coordinates of each
individual as updated along the chain if those given as input are not
considered as exact coordinates (which is specified by 
\code{delta.coord} to a non zero value).

\item  File \file{log.likelihood.txt} contains log-likelihood of data
for the current state of parameters of the Markov chain.

\item  File \file{log.posterior.density.txt} contains log of posterior probability
(up to marginal density of data) of the
current state of parameters in the Markov chain.

}}

\references{
  \itemize{
    \item G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
    model for landscape genetics. Genetics, 170, 1261-1280, 2005.
    
    \item  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
    genetics. Molecular Ecology Notes, 5, 712-715, 2005.

    
    \item   Gilles Guillot, Filipe Santos  and Arnaud Estoup,
    Analysing georeferenced population genetics data with
    Geneland: a new algorithm to deal with null alleles and a friendly graphical user interface
    Bioinformatics 2008 24(11):1406-1407.
    
    \item G. Guillot.  Inference of structure in subdivided populations at
    low levels of genetic differentiation. The correlated allele
    frequencies 
    model revisited. Bioinformatics, 24:2222-2228, 2008
  
    \item G. Guillot and F. Santos  A computer program to simulate
    multilocus 
    genotype data with spatially auto-correlated allele frequencies.
    Molecular Ecology Resources, 2009

    \item G. Guillot, R. Leblois, A. Coulon, A. Frantz  Statistical
    methods in spatial genetics, Molecular Ecology, 2009.

    \item B. Guedj and G. Guillot. Estimating the location and shape of hybrid
  zones. Molecular Ecology Resources, 11(6) 1119-1123, 2011

   \item G. Guillot, S. Renaud, R. Ledevin, J. Michaux and J. Claude. A
  Unifying Model for the Analysis of Phenotypic, Genetic and Geographic
  Data. Systematic Biology, to appear, 2012.

}}


\author{Gilles Guillot  }



\seealso{ \code{\link{simFmodel}}}
