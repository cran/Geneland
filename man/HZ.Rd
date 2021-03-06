\name{HZ}
\alias{HZ}

\title{Estimate parameters of the hybrid zone model by MCMC simulation}
\description{Estimate parameters of a hybrid zone model by MCMC
  simulation. The function does not currently accept more than one genotype matrix.
}
\usage{
            HZ(coordinates,
               geno.dip.codom=NULL,
               geno.dip.dom=NULL,
               geno.hap=NULL,
               dist.IC=NULL,
               allele.freq=NULL,
               ncluster=NULL,
               cluster.indiv=NULL,
               path.mcmc.noadm=NULL,
               a.init=NULL,
               b.init=NULL,
               c.init=1,
               a.max=10,
               b.max=NULL,
               c.max=1,
               estimate.a=TRUE,
               estimate.b=TRUE,
               estimate.c=FALSE,
               common.param=TRUE,
               nit,
               thinning,
               path.mcmc.adm=NULL)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with
  2}

  \item{geno.dip.codom}{Genotypes  for diploid data with codominant markers.
    A matrix with one line per individual and two  columns per locus.
    }
  
  \item{geno.dip.dom}{Genotypes for diploid data with dominant
    markers. A matrix with one line per individual and one  column per
    locus. Presence/absence of a band should be
    coded as 0/1 (0 for absence / 1 for presence). Dominant and
    codominant markers can be analyzed jointly by passing variables to arguments 
    geno.dip.codom and geno.dip.dom.
    Haploid data and diploid dominant data can not be analyzed jointly in
    the current version.}
  
  \item{geno.hap}{Genotypes of haploid data.
    A matrix with one line per individual and one  column per
    locus.     Dominant diploid data and haploid data 
    can be analyzed jointly (e.g. to analyse microsatelite data or SNP
    data together with mtDNA.
    Haploid data and diploid dominant data can not be analyzed jointly in
    the current version.}

 \item{dist.IC}{A matrix with \code{nindiv} lines and \code{ncluster}
   column. The parameter \code{ncluster} being the number of clusters,
   most often 2 when there is a single hybrid zone. 
   If \code{dist.IC} is missing, the user has to provide instead the
   path to a directory storing results from a no-admixture MCMC run. }
 
  \item{allele.freq}{An array with \code{ncluster} x \code{nloc} x \code{nalmax}
  If missing, it will be estimated from the output from the MCMC run to
  estimate clusters. }

  \item{ncluster}{Number of clusters. If missing, the user has to provide instead the
   path to a directory storing results from a no-admixture MCMC run}

  \item{cluster.indiv}{Cluster membership of individuals. A numeric
   vector with integer values (maximum values being the total number of
   clusters) }
 
  \item{path.mcmc.noadm}{ Path to output files directory of the previous
   Geneland    no-admixture run. It seems that the
    path should be given in the Unix style even under Windows (use \/
    instead of \\).
    This path *has to*  end with a slash (\/)
    (e.g. path.mcmc="/home/me/Geneland-noadmixture/")  }

  \item{a.init}{A numerical value to use as fixed or initial value for the \code{a} parameter}
  
  \item{b.init}{A numerical value to use as fixed or initial value for the \code{b} parameter}
  
  \item{c.init}{A numerical value to use as fixed or initial value for the \code{c} parameter}

  \item{a.max}{Maximum value allowed along MCMC simulation for parameter
    \code{a} (default is 1)}

  \item{b.max}{Maximum value allowed along MCMC simulation for parameter
    \code{b} (default is a small fraction of the study area diameter)}
  
  \item{c.max}{Maximum value allowed along MCMC simulation for parameter
    \code{c} (default is 1)}
  
  \item{estimate.a}{Logical. If TRUE, parameter \code{a} is estimated, if
    FALSE it is left at the initial value.}
  
  \item{estimate.b}{Logical. If TRUE, parameter \code{b} is estimated, if
    FALSE it is left at the initial value.}

  \item{estimate.c}{Logical. If TRUE, parameter \code{c} is estimated, if
    FALSE it is left at the initial value.}

  \item{common.param}{If TRUE, parameters \code{a}, \code{b} and \code{c}
    are common to all clusters. If FALSE, the program attemps to estimate
    cluster-specific values.}

  \item{nit}{Number of MCMC iteration.}

  \item{thinning}{Number of MCMC iterations between two writing steps (if \code{thinning}=1, all
    states are saved whereas if e.g. \code{thinning}=10 only each 10
    iteration is saved)}

  \item{path.mcmc.adm}{ Path to output files directory for the admixture model. It seems that the
    path should be given in the Unix style even under Windows (use \/ instead of \\).
    This path *has to*  end with a slash (\/)    (e.g. path.mcmc="/home/me/Geneland-admixture/")
  }
}

\value{No object is returned. All outputs are stored in ascii file
  located in the \code{path.mcmc} directory}

\references{B. Guedj and G. Guillot, A Bayesian model for infering
  hybrid zones. 
%% ~put references to the literature/web site here ~
}
\author{G. Guillot}

%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hybrid zone}
\keyword{Admixture}
\keyword{MCMC}
\keyword{Spatial genetics}