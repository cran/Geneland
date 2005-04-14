"PlotFreqA" <-
function(allele.numbers,path.mcmc,iloc,iall,printit=FALSE,path=NULL)
  {
    coordinates <- as.matrix(coordinates)
    genotypes <- as.matrix(genotypes)
    allele.numbers <- as.matrix(allele.numbers)
    
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])

    nloc <- length(allele.numbers)
    filefa <-  paste(path.mcmc,"ancestral.frequencies.txt",sep="")
    fa <- as.matrix(read.table(filefa))
    X11()
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

