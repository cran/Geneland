`PlotFreqA` <-
function(genotypes,path.mcmc,iloc,iall,printit=FALSE,path)
  {
   
    
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])
    ploidy <- as.numeric(param[param[,1]=="ploidy",3])
    if(ploidy == 1)
      {
        ## diploidize the data
        data.tmp <- matrix(nrow=nrow(genotypes),
                           ncol=ncol(genotypes)*2)
        data.tmp[,seq(1,ncol(genotypes)*2-1,2)] <- genotypes
        data.tmp[,seq(2,ncol(genotypes)*2,2)] <- genotypes
        genotypes <- data.tmp
       }
    allele.numbers <- FormatGenotypes(genotypes)$allele.numbers
    nloc <- length(allele.numbers)
    filefa <-  paste(path.mcmc,"ancestral.frequencies.txt",sep="")
    fa <- as.matrix(read.table(filefa))
    get(getOption("device"))()
    plot(fa[,(iloc-1)*max(allele.numbers)+iall],
         xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
         ylab=paste("Frequency of allele",
           iall,"at locus",iloc),type="l")
    title(main=ifelse(iall==1,
            paste("Allele freq. in ancestral pop. at locus",iloc),
            ""))
    if(printit==T)
      {
        postscript(file=paste(path,"freq.ancestral.pop.loc",iloc,".ps",sep=""))
        plot(fa[,(iloc-1)*max(allele.numbers)+iall],
             xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
             ylab=paste("Frequency of allele",
               iall,"at locus",iloc),type="l")
        title(main=ifelse(iall==1,
                paste("Allele freq. in ancestral pop. at locus",iloc),
                ""))
        dev.off()
      }

  }

