"mcmcFmodel" <-
function(# path to input directory
                 repdat,
                 # path to output directory
                 repmcmc,
                 # hyper-prior parameters
                 lambdamax,dt,nclassmin,nclassinit,nclassmax,
                 # dimensions
                 nppmax,
                 # options in mcmc computations
                 nchain,stepw,model,varnclass,spatial)
  {
    
                                        # input files
    files <- paste(repdat,"coordinates.txt",sep="")
    filez <- paste(repdat,"genotypes.txt",sep="")
    filenall <- paste(repdat,"allele.numbers.txt",sep="")
    
                                        # output files
    filenpp <- paste(repmcmc,"nuclei.numbers.txt",sep="")
    fileu <- paste(repmcmc,"coord.nuclei.txt",sep="")
    filec <- paste(repmcmc,"color.nuclei.txt",sep="")
    filef <- paste(repmcmc,"frequencies.txt",sep="")
    filefa <- paste(repmcmc,"ancestral.frequencies.txt",sep="")
    filedrift <- paste(repmcmc,"drifts.txt",sep="")
    filenclass <- paste(repmcmc,"populations.numbers.txt",sep="")
    filelambda <- paste(repmcmc,"Poisson.process.rate.txt",sep="")
    filelpp <- paste(repmcmc,"log.posterior.density.txt",sep="")
    filell <- paste(repmcmc,"log.likelihood.txt",sep="")
    
    nindiv <- nrow(read.table(filez))
    nloc <- ncol(read.table(filez))/2
    
    npp <- ifelse(spatial==1,
                  1+floor(runif(1)*nppmax/2), #avoid large init for npp w.r.t lambda
                  nindiv)
    nallmax <- max(scan(filenall))
    s <- matrix(nr=2,nc=nindiv,data=-999)
    z <- matrix(nr=nindiv,nc=nloc*2,data=-999)
    u <- matrix(nr=2,nc=nppmax,data=-999)
    utemp <- matrix(nr=2,nc=nppmax,data=-999)
    c <- rep(times=nppmax,-999)
    ctemp <- rep(times=nppmax,-999)
    t <-matrix(nr=2,nc=nindiv,data=-999)
    ttemp <-matrix(nr=2,nc=nindiv,data=-999)
    f <- array(dim=c(nclassmax,nloc,nallmax),data=-999)
    ftemp <- array(dim=c(nclassmax,nloc,nallmax),data=-999)
    nall <- scan(filenall)
    fa <- array(dim=c(nloc,nallmax),data=-999)
    drift <- rep(-999,nclassmax)
    drifttemp <- rep(-999,nclassmax)
    indcell <- rep(times=nindiv,-999)
    indcelltemp <- rep(times=nindiv,-999)
    distcell <- rep(times=nindiv,-999)
    distcelltemp <- rep(times=nindiv,-999)
    n <-  array(dim=c(nclassmax,nloc,nallmax),data=-999)
    ntemp <-  array(dim=c(nclassmax,nloc,nallmax),data=-999)
    a <- rep(times=nallmax,-999)
    ptemp <- rep(times=nallmax,-999)
    effcl <- rep(times=nclassmax,-999)
    iclv <- rep(times=nclassmax,-999)
    cellclass <- rep(times=nppmax,-999)
    listcell <- rep(times=nppmax,-999)
    fmodel <- ifelse(model=="Falush",1,0) # Falush or Dirichlet model
    kfix <- 1-as.integer(varnclass)
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
                       as.character(filenclass),
                       as.character(filelpp),
                       as.character(filell),
                       as.single(lambdamax),
                       as.single(dt),
                       as.integer(nchain),
                       as.integer(stepw),
                       as.integer(nindiv),
                       as.integer(nloc),
                       as.integer(nloc*2),
                       as.integer(nallmax),
                       as.integer(npp),
                       as.integer(nppmax),
                       as.integer(nclassinit),
                       as.integer(nclassmin),
                       as.integer(nclassmax),
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
    param <- c(paste("repdat :",repdat),
               paste("repmcmc :",repmcmc),
               paste("lambdamax :",lambdamax),
               paste("dt :",dt),
               paste("nclassmin :",nclassmin),
               paste("nclassinit :",nclassinit),
               paste("nclassmax :",nclassmax),
               paste("nindiv :",nindiv),
               paste("nloc :",nloc),
               paste("nppmax :",nppmax),
               paste("nchain :",nchain),
               paste("stepw :",stepw),
               paste("model :",model),
               paste("varnclass :",varnclass),
               paste("spatial :",spatial))
    write.table(param,file=paste(repmcmc,"parameters.txt",sep=""),
                quote=FALSE,row.name=FALSE,col.name=FALSE)
    NULL               
                   
}

