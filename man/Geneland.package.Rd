\name{Geneland-package}

\alias{Geneland-package}

\alias{Geneland}

\docType{package}

\title{Simulation and inference of a 
  spatial statistical model for landscape genetics.}

\description{
  The main function \code{\link{mcmcFmodel}} 
  takes geo-referenced individual multilocus genetic data and tries to detect
  population structure, i.e sub-populations making use of both spatial and genetic information.
}
}
\details{
  \tabular{ll}{
    Package: \tab Geneland\cr
    Type: \tab Package\cr
    Version: \tab 1.0.5\cr
    Date: \tab 2006-01-10\cr
    License: \tab GPL\cr
  }
  
  The overall population is viewed as a co-existence of sub-populations at 
  Hardy-Weiberg equilibrium. The sub-populations are supposed to be spatially 
  organised through the so-called 
  colored Poisson-Voronoi tessellation. 
  Allele frequencies are assumed to be drawn from the F-model 
  as described in (Falush 2003) although 
  the particular case of independent Dirichlet allele frequencies, 
  as described by Pritchard (2000),  is also handled. 
  Individuals within populations are assumed to be randomly 
  located and  Hardy-Weinberg equilibrium 
  and linkage equilibrium are assumed.
  
  The main purpose of the program is to perform Bayesian inference of all the parameters involved 
  through Markov Chain Monte-Carlo simulation.
  This is achievied by the function   \code{\link{mcmcFmodel}}.
  Function \code{\link{PostProcessChain}} read some output files of
  \code{\link{mcmcFmodel}}
  and computes some statistics suitable to print maps of inferred
  populations.
  
  See \code{Storage format} section in \code{\link{mcmcFmodel} help page.}

  The following functions are provided by the package:

  
  \code{\link{simFmodel}}: simulation from the prior of the spatial
  F-model
  
  \code{\link{mcmcFmodel}}: Full Bayesian Markov Chain Monte Carlo
  inference of parameters in the spatial F-model
  
  \code{\link{PostProcessChain}}: Post-procesing of MCMC output
  for maps of posterior probability of populations subdomains
    
  \code{\link{PlotTessellation}}: Graphical display of inferred
  sub-domains

  The following  functions are very basic and are only intended to be an
  aid for those not familiar with R. Most probably you may want to use
  directly the output files of \code{\link{mcmcFmodel}} and
  \code{\link{PostProcessChain}} to print your own figures.
  
  \code{\link{PlotDrift}}: Graphical display of drift factors along MCMC
  run
  
   \code{\link{PlotFreqA}}: Graphical display of allele frequencies in the
   ancestral population along MCMC run
   
   \code{\link{PlotFreq}}: Graphical display of allele frequencies in the
   present time population along MCMC run
   
   \code{\link{Plotnpop}}: Graphical display of number of populations
   along MCMC run
   
   \code{\link{Plotntile}}: Graphical display of number of tiles along
   MCMC run
   
   \code{\link{PosteriorMode}}: Computation and/or graphical display of
   mode in the posterior distribution of class membership at each pixel
   
   \code{\link{Fstat}}: Computations of pairwise F statistics between
   inferred subpopulations 

   \code{\link{FormatGenotypes}}: Internal function. 
   Transform a file of genotypes into
   a format suitable for function \code{\link{mcmcFmodel}}

   \code{\link{setplot}}: Internal function

}
\author{Gilles Guillot

  \url{www.inapg.inra.fr/ens_rech/mathinfo/personnel/guillot/welcome.html}

  See also
  
  \url{www.inapg.inra.fr/ens_rech/mathinfo/personnel/guillot/Geneland.html}
}

\references{

  On the implementation of mixture models in population genetics:
  
- J.K. Pritchard, M. Stephens and P. Donnelly,
Inference of population structure using multilocus genotype data,
Genetics, pp 945-959 vol. 155, 2000

- Falush D., M. Stephens and  J.K. Pritchard, 
Inference of population structure using multilocus genotype data:
Linked loci and correlated allele frequencies, Genetics, pp 1567-1587,
vol 164, 2003

On the implementation of variable dimension MCMC algorihtm in
population genetics:

- Corander, J.C., Waldmann, P. and Sillanpaa, M.J.,
Bayesian analysis of genetic differentiation between populations,
Genetics, 2003, 163, 367-374

- Corander, J.C.,  P. Waldmann, P. Martinen and  M.J. Sillanpaa, ,
BAPS2: Enhanced possibilities for the analysis of genetic population structure,
Bioinformatics, vol. 20,number 15, 2004

On the use of Voronoi tessellations in population genetics :

- Dupanloup, I., Schneider, S. and Excoffier, L.,
A simulated annealing approach to define genetic structure of populations,
Molecular Ecology, 2002, 11, 2571-2581

On the model (and sub-models) implemented in \code{Geneland}

- G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
model for landscape genetics. Genetics, 170, 1261-1280, 2005.


- G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
genetics. Molecular Ecology Notes, 5, 712-715, 2005.

}
}
\keyword{ }

\seealso{}

\examples{}
