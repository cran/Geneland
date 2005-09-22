"mcmcFmodel" <-
function(
                                        # input data
                 coordinates,genotypes,ploidy=2,
                 # path to output directory
                 path.mcmc,
                 # hyper-prior parameters
                 rate.max,delta.coord,npopmin,npopinit,npopmax,
                 # dimensions
                 nb.nuclei.max,
                 # options in mcmc computations
                 nit,thinning=1,freq.model="Dirichlet",varnpop,spatial=TRUE)
  {
    ## check parameters
    if((nit %% thinning) != 0)
      {
        stop('nit/thinning is not an integer')
      }
    if(file.exists(path.mcmc) != TRUE)
      {
        stop('The output directory does not exist')
      }
    

    ##
    ## Reformatting data
    if(ploidy == 1)
      {
        ## diploidize the data
        data.tmp <- matrix(nrow=nrow(genotypes),
                            ncol=ncol(genotypes)*2)
        data.tmp[,seq(1,ncol(genotypes)*2-1,2)] <- genotypes
        data.tmp[,seq(2,ncol(genotypes)*2,2)] <- genotypes
        genotypes <- data.tmp
      }
    data.tmp <- FormatGenotypes(as.matrix(genotypes))
    coordinates <- as.matrix(coordinates)
    genotypes <- data.tmp$genotypes
    allele.numbers <- data.tmp$allele.numbers

    

    ##
    ## Initializing variables
    
                                        # output files
    filenpp <- paste(path.mcmc,"nuclei.numbers.txt",sep="")
    fileu <- paste(path.mcmc,"coord.nuclei.txt",sep="")
    filec <- paste(path.mcmc,"color.nuclei.txt",sep="")
    filef <- paste(path.mcmc,"frequencies.txt",sep="")
    filefa <- paste(path.mcmc,"ancestral.frequencies.txt",sep="")
    filedrift <- paste(path.mcmc,"drifts.txt",sep="")
    filenpop <- paste(path.mcmc,"populations.numbers.txt",sep="")
    filelambda <- paste(path.mcmc,"Poisson.process.rate.txt",sep="")
    filet <-  paste(path.mcmc,"hidden.coord.txt",sep="")
    filelpp <- paste(path.mcmc,"log.posterior.density.txt",sep="")
    filell <- paste(path.mcmc,"log.likelihood.txt",sep="")


    nindiv <- nrow(genotypes)
    nloc <- length(allele.numbers)
    
    npp <- ifelse(spatial==1,
                  1+floor(runif(1)*nb.nuclei.max/2), #avoid large init for npp w.r.t lambda
                  nindiv)
    nallmax <- max(allele.numbers)
    u <- matrix(nr=2,nc=nb.nuclei.max,data=-999)
    utemp <- matrix(nr=2,nc=nb.nuclei.max,data=-999)
    c <- rep(times=nb.nuclei.max,-999)
    ctemp <- rep(times=nb.nuclei.max,-999)
    t <-matrix(nr=2,nc=nindiv,data=-999)
    ttemp <-matrix(nr=2,nc=nindiv,data=-999)
    f <- array(dim=c(npopmax,nloc,nallmax),data=-999)
    ftemp <- array(dim=c(npopmax,nloc,nallmax),data=-999)
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

    ##
    ## Starting MCMC computations

    out.res<- .Fortran(name="mcmc",
                       PACKAGE="Geneland",
                       as.double(t(coordinates)),
                       as.integer(genotypes),
                       as.integer(allele.numbers),
                       as.integer(ploidy),
                       as.character(filelambda),
                       as.character(filenpp),
                       as.character(fileu),
                       as.character(filec),
                       as.character(filef),
                       as.character(filefa),
                       as.character(filedrift),
                       as.character(filenpop),
                       as.character(filet),
                       as.character(filelpp),
                       as.character(filell),
                       as.double(rate.max),
                       as.double(delta.coord),
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
                       as.double(t),
                       as.double(ttemp),
                       as.double(u),
                       as.double(utemp),
                       as.integer(c),
                       as.integer(ctemp),
                       as.double(f),
                       as.double(ftemp),
                       as.double(fa),
                       as.double(drift),
                       as.double(drifttemp),
                       as.integer(indcell),
                       as.integer(indcelltemp),
                       as.double(distcell),
                       as.double(distcelltemp),
                       as.integer(n),
                       as.integer(ntemp),
                       as.double(a),
                       as.double(ptemp),
                       as.integer(effcl),
                       as.integer(iclv),
                       as.integer(cellclass),
                       as.integer(listcell),
                       as.integer(fmodel),
                       as.integer(kfix),
                       as.integer(spatial))

                                        # write parameters of the present run in an ascii file
    param <- c(#paste("path.mcmc :",path.mcmc),
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

