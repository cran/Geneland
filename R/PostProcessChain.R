"PostProcessChain" <-
function(path.data,# path to data directory
                             path.mcmc, # path to MCMC output directory
                             nxdom,nydom,# resolution
                             burnin # number of iterations of the chain to throw away
                             )
  {
    files <- paste(path.data,"coordinates.txt",sep="")
    filez <- paste(path.data,"genotypes.txt",sep="")
    filenall <- paste(path.data,"allele.numbers.txt",sep="")

    filenpop <- paste(path.mcmc,"populations.numbers.txt",sep="")
    filenpp <- paste(path.mcmc,"nuclei.numbers.txt",sep="")
    fileu <- paste(path.mcmc,"coord.nuclei.txt",sep="")
    filec <- paste(path.mcmc,"color.nuclei.txt",sep="")
    filef <- paste(path.mcmc,"frequencies.txt",sep="")
    
    filefperm <- paste(path.mcmc,"fperm.txt",sep="")
    filedom <- paste(path.mcmc,"proba.pop.membership.txt",sep="")
    filedomperm <- paste(path.mcmc,"proba.pop.membership.perm.txt",sep="")
    
    nindiv <- nrow(read.table(filez))
    nloc <- ncol(read.table(filez))/2

    param <- as.matrix(read.table(paste(path.mcmc,"parameters.txt",sep="")))
    delta.coord <- as.numeric(param[4,3])
    npopmax <- as.numeric(param[7,3])
    nb.nuclei.max <- as.numeric(param[10,3])
    nit <- as.numeric(param[11,3])
    thinning <- as.numeric(param[12,3])
    
    dom <- matrix(nr=nxdom*nydom,nc=npopmax,data=0)
    domperm <- matrix(nr=nxdom*nydom,nc=npopmax,data=0)
    coorddom <- matrix(nr=2,nc=nxdom*nydom,data=-999)
    indvois <- numeric(nxdom*nydom)
    distvois <- numeric(nxdom*nydom)
    f11 <- numeric(npopmax)
    orderf11 <- numeric(npopmax)
    nallmax <- max(scan(filenall))
    s <- matrix(nr=2,nc=nindiv,data=-999)
    u <- matrix(nr=2,nc=nb.nuclei.max,data=-999)
    c <- rep(times=nb.nuclei.max,-999)
    nall <- scan(filenall)
    
    out.res<- .Fortran(name="postprocesschain",
                       PACKAGE="Geneland",
                       as.integer(nindiv),
                       as.integer(nxdom),
                       as.integer(nydom),
                       as.integer(burnin),
                       as.integer(npopmax),
                       as.integer(nb.nuclei.max),
                       as.integer(nloc),
                       as.integer(nindiv),
                       as.integer(nloc),
                       as.character(filenall),
                       as.integer(nallmax),
                       as.character(files),
                       as.single(delta.coord),
                       as.integer(nit/thinning),
                       as.character(filenpp),
                       as.character(filenpop),
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
    write.table(param,file=paste(path.mcmc,"postprocess.parameters.txt",sep=""),
                quote=FALSE,row.name=FALSE,col.name=FALSE)
  }

