"PostProcessChain" <-
function(repdat,# path to data directory
                             repmcmc, # path to MCMC output directory
                             nxdom,nydom,# resolution
                             burnin # number of iterations of the chain to throw away
                             )
  {
    files <- paste(repdat,"coordinates.txt",sep="")
    filez <- paste(repdat,"genotypes.txt",sep="")
    filenall <- paste(repdat,"allele.numbers.txt",sep="")

    filenclass <- paste(repmcmc,"populations.numbers.txt",sep="")
    filenpp <- paste(repmcmc,"nuclei.numbers.txt",sep="")
    fileu <- paste(repmcmc,"coord.nuclei.txt",sep="")
    filec <- paste(repmcmc,"color.nuclei.txt",sep="")
    filef <- paste(repmcmc,"frequencies.txt",sep="")
    
    filefperm <- paste(repmcmc,"fperm.txt",sep="")
    filedom <- paste(repmcmc,"proba.pop.membership.txt",sep="")
    filedomperm <- paste(repmcmc,"proba.pop.membership.perm.txt",sep="")
    
    nindiv <- nrow(read.table(filez))
    nloc <- ncol(read.table(filez))/2

    param <- as.matrix(read.table(paste(repmcmc,"parameters.txt",sep="")))
    dt <- as.numeric(param[4,3])
    nclassmax <- as.numeric(param[7,3])
    nppmax <- as.numeric(param[10,3])
    nchain <- as.numeric(param[11,3])
    stepw <- as.numeric(param[12,3])
    
    dom <- matrix(nr=nxdom*nydom,nc=nclassmax,data=0)
    domperm <- matrix(nr=nxdom*nydom,nc=nclassmax,data=0)
    coorddom <- matrix(nr=2,nc=nxdom*nydom,data=-999)
    indvois <- numeric(nxdom*nydom)
    distvois <- numeric(nxdom*nydom)
    f11 <- numeric(nclassmax)
    orderf11 <- numeric(nclassmax)
    nallmax <- max(scan(filenall))
    s <- matrix(nr=2,nc=nindiv,data=-999)
    u <- matrix(nr=2,nc=nppmax,data=-999)
    c <- rep(times=nppmax,-999)
    nall <- scan(filenall)
    
    out.res<- .Fortran(name="postprocesschain",
                       PACKAGE="Geneland",
                       as.integer(nindiv),
                       as.integer(nxdom),
                       as.integer(nydom),
                       as.integer(burnin),
                       as.integer(nclassmax),
                       as.integer(nppmax),
                       as.integer(nloc),
                       as.integer(nindiv),
                       as.integer(nloc),
                       as.character(filenall),
                       as.integer(nallmax),
                       as.character(files),
                       as.single(dt),
                       as.integer(nchain/stepw),
                       as.character(filenpp),
                       as.character(filenclass),
                       as.character(fileu),
                       as.character(filec),
                       as.character(filef),
                       as.character(filefperm),
                       as.character(filedom),
                       as.character(filedomperm),
                       as.single(s),
                       as.single(u),
                       as.integer(c),
                       as.single(dom),
                       as.single(domperm),
                       as.single(coorddom),
                       as.integer(indvois),
                       as.single(distvois),
                       as.integer(nall),
                       as.single(f11),
                       as.single(orderf11))
    param <- c(paste("nxdom :",nxdom),
               paste("nydom :",nydom))
    write.table(param,file=paste(repmcmc,"postprocess.parameters.txt",sep=""),
                quote=FALSE,row.name=FALSE,col.name=FALSE)
  }

