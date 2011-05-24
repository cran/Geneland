\name{show.estimate.hz}

\alias{show.estimate.hz}

\title{Show estimate of parameters of the hybrid zone model by MCMC simulation}
\description{Show estimate of parameters of a hybrid zone model.}
\usage{
show.estimate.hz(coordinates,
   path.mcmc.adm,
   burnin,
   angle=0,
   plot.distruct=TRUE,
   plot.mcmc=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{Spatial coordinates of individuals. A matrix with
  2}

  \item{path.mcmc.adm}{ Path to output files directory for the admixture model. It seems that the
    path should be given in the Unix style even under Windows (use \/ instead of \\).
    This path *has to*  end with a slash (\/)    (e.g. path.mcmc="/home/me/Geneland-admixture/")
  }
  \item{burnin}{burnin}
  \item{angle}{angle of the axis on which sampling sites are projected
    (in radian)}
  \item{plot.distruct}{Logical}
  \item{plot.mcmc}{Logical}
}

\value{No object is returned. All outputs are stored in ascii files 
  located in the \code{path.mcmc} directory}

\references{B. Guedj and G. Guillot, 
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