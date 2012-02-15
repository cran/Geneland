\name{Geneland-package}

\alias{Geneland-package}

\alias{Geneland}

\docType{package}

\title{Simulation and inference for subdivided populations }

\description{
  Detect population structure (i.e sub-populations),
  making use of genetic (and optionally geographic) information.
}


\author{Arnaud Estoup, Gilles Guillot, Filipe Santos

  \url{http://www2.imm.dtu.dk/~gigu/Geneland/}

  See also
  
  \url{http://www2.imm.dtu.dk/~gigu/Geneland/}
}

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
    }
  
}


\details{
 

  
  The main purpose of the program is to perform Bayesian inference of all the parameters involved 
  through Markov Chain Monte-Carlo simulation.
  This is achievied by the function   \code{\link{MCMC}}.
  Function \code{\link{PostProcessChain}} read some output files of
  \code{\link{MCMC}}
  and computes some statistics suitable to print maps of inferred
  populations.
  
  See \code{Storage format} section in \code{\link{MCMC} help page.}

  The following functions are provided by the package:

  
  \code{\link{simFmodel}}: simulation from the prior of the spatial
  F-model

  \code{\link{simdata}}: Simulation of georeferenced genotypes under an IBD + barrier
  model

   \code{\link{show.simdata}}:  Graphical display of data simulated by simdata
  
  \code{\link{MCMC}}: Full Bayesian Markov Chain Monte Carlo
  inference of parameters in the spatial F-model
  
  \code{\link{PostProcessChain}}: Post-procesing of MCMC output
  for maps of posterior probability of populations subdomains
    
  \code{\link{PlotTessellation}}: Graphical display of inferred
  sub-domains

  The following  functions are very basic and are only intended to be an
  aid for those not familiar with R. Most probably you may want to use
  directly the output files of \code{\link{MCMC}} and
  \code{\link{PostProcessChain}} to print your own figures.
  
  \code{\link{PlotDrift}}: Graphical display of drift factors along MCMC
  run
  
   \code{\link{PlotFreqA}}: Graphical display of allele frequencies in the
   ancestral population along MCMC run
   
   \code{\link{PlotFreq}}: Graphical display of allele frequencies in the
   present time population along MCMC run
   
   \code{\link{Plotnpop}}: Graphical display of number of populations
   along MCMC run
   
%    \code{\link{Plotntile}}: Graphical display of number of tiles along
%    MCMC run
   
%    \code{\link{PosteriorMode}}: Computation and/or graphical display of
%    mode in the posterior distribution of class membership at each pixel
   
%    \code{\link{Fstat}}: Computations of pairwise F statistics between
%    known subpopulations of a dataset

%     \code{\link{Fstat.output}}: Computations of pairwise F statistics between
%    inferred subpopulations 

%     \code{\link{nullify}}: Alter genotpyes of a dataset by simulating
%     null alleles.
    
%    \code{\link{FormatGenotypes}}: Internal function. 
%    Transform a file of genotypes into
%    a format suitable for function \code{\link{MCMC}}

 \tabular{ll}{
    Package: \tab Geneland\cr
    Type: \tab Package\cr
    License: \tab GPL\cr
  }
  
   
}