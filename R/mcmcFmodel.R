"mcmcFmodel" <-
function(# path to input directory
                 path.data,
                 # path to output directory
                 path.mcmc,
                 # hyper-prior parameters
                 rate.max,delta.coord,npopmin,npopinit,npopmax,
                 # dimensions
                 nb.nuclei.max,
                 # options in mcmc computations
                 nit,thinning,freq.model,varnpop,spatial)
  {
    
                                        # input files
    files <- paste(path.data,"coordinates.txt",sep="")
    filez <- paste(path.data,"genotypes.txt",sep="")
    filenall <- paste(path.data,"allele.numbers.txt",sep="")
    
                                        # output files
    filenpp <- paste(path.mcmc,"nuclei.numbers.txt",sep="")
    fileu <- paste(path.mcmc,"coord.nuclei.txt",sep="")
    filec <- paste(path.mcmc,"color.nuclei.txt",sep="")
    filef <- paste(path.mcmc,"frequencies.txt",sep="")
    filefa <- paste(path.mcmc,"ancestral.frequencies.txt",sep="")
    filedrift <- paste(path.mcmc,"drifts.txt",sep="")
    filenpop <- paste(path.mcmc,"populations.numbers.txt",sep="")
    filelambda <- paste(path.mcmc,"Poisson.process.rate.txt",sep="")
    filelpp <- paste(path.mcmc,"log.posterior.density.txt",sep="")
    filell <- paste(path.mcmc,"log.likelihood.txt",sep="")
    
    nindiv <- nrow(read.table(filez))
    nloc <- ncol(read.table(filez))/2
    
    npp <- ifelse(spatial==1,
                  1+floor(runif(1)*nb.nuclei.max/2), #avoid large init for npp w.r.t lambda
                  nindiv)
    nallmax <- max(scan(filenall,quiet=TRUE))
    s <- matrix(nr=2,nc=nindiv,data=-999)
    z <- matrix(nr=nindiv,nc=nloc*2,data=-999)
    u <- matrix(nr=2,nc=nb.nuclei.max,data=-999)
    utemp <- matrix(nr=2,nc=nb.nuclei.max,data=-999)
    c <- rep(times=nb.nuclei.max,-999)
    ctemp <- rep(times=nb.nuclei.max,-999)
    t <-matrix(nr=2,nc=nindiv,data=-999)
    ttemp <-matrix(nr=2,nc=nindiv,data=-999)
    f <- array(dim=c(npopmax,nloc,nallmax),data=-999)
    ftemp <- array(dim=c(npopmax,nloc,nallmax),data=-999)
    nall <- scan(filenall,quiet=TRUE)
    fa <- array(dim=c(nloc,nallmax),data=-999)
    drift <- rep(-999,npopmax)
    drifttemp <- rep(-999,npopmax)
    indcell <- rep(times=nindiv,-999)
    indcelltemp <- rep(times=nindiv,-999)
    distcell <- rep(times=nindiv,-999)
    distcelltemp <- rep(times=nindiv,-999)
    n <-  array(dim=c(npopmax,nloc,nallmax),data=-999)
    ntemp <-  array(dim=c(npopmax,nloc,nallmax),data=-999)
    a <- rep(times=nallmax,-999)
    ptemp <- rep(times=nallmax,-999)
    effcl <- rep(times=npopmax,-999)
    iclv <- rep(times=npopmax,-999)
    cellclass <- rep(times=nb.nuclei.max,-999)
    listcell <- rep(times=nb.nuclei.max,-999)
    fmodel <- ifelse(freq.model=="Falush",1,0) # Falush or Dirichlet model
    kfix <- 1-as.integer(varnpop)
    spatial <- as.integer(spatial)
    

    out.res<- .Fortran(name="mcmc",
                       PACKAGE="Geneland",
                       as.character(files),
                       as.character(filez),
                       as.character(filenall),
                       as.character(filelambda),
                       as.character(filenpp),
                       as.character(fileu),
                       as.character(filec),
                       as.character(filef),
                       as.character(filefa),
                       as.character(filedrift),
                       as.character(filenpop),
                       as.character(filelpp),
                       as.character(filell),
                       as.single(rate.max),
                       as.single(delta.coord),
                       as.integer(nit),
                       as.integer(thinning),
                       as.integer(nindiv),
                       as.integer(nloc),
                       as.integer(nloc*2),
                       as.integer(nallmax),
                       as.integer(npp),
                       as.integer(nb.nuclei.max),
                       as.integer(npopinit),
                       as.integer(npopmin),
                       as.integer(npopmax),
                       as.single(s),
                       as.integer(z),
                       as.single(t),
                       as.single(ttemp),
                       as.single(u),
                       as.single(utemp),
                       as.integer(c),
                       as.integer(ctemp),
                       as.single(f),
                       as.single(ftemp),
                       as.single(fa),
                       as.single(drift),
                       as.single(drifttemp),
                       as.single(nall),
                       as.integer(indcell),
                       as.integer(indcelltemp),
                       as.single(distcell),
                       as.single(distcelltemp),
                       as.integer(n),
                       as.integer(ntemp),
                       as.single(a),
                       as.single(ptemp),
                       as.integer(effcl),
                       as.integer(iclv),
                       as.integer(cellclass),
                       as.integer(listcell),
                       as.integer(fmodel),
                       as.integer(kfix),
                       as.integer(spatial))

                                        # write parameters of the run in an ascii file
    param <- c(paste("path.data :",path.data),
               paste("path.mcmc :",path.mcmc),
               paste("rate.max :",rate.max),
               paste("delta.coord :",delta.coord),
               paste("npopmin :",npopmin),
               paste("npopinit :",npopinit),
               paste("npopmax :",npopmax),
               paste("nindiv :",nindiv),
               paste("nloc :",nloc),
               paste("nb.nuclei.max :",nb.nuclei.max),
               paste("nit :",nit),
               paste("thinning :",thinning),
               paste("freq.model :",freq.model),
               paste("varnpop :",varnpop),
               paste("spatial :",spatial))
    write.table(param,file=paste(path.mcmc,"parameters.txt",sep=""),
                quote=FALSE,row.name=FALSE,col.name=FALSE)
#    NULL               
                   
}

