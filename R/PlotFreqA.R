"PlotFreqA" <-
function(repdat,repmcmc,iloc,iall,printit=FALSE,path=NULL)
  {
    fileparam <- paste(repmcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    stepw <- as.numeric(param[12,3])

     
    nall <- scan(paste(repdat,"allele.numbers.txt",sep=""))
    nloc <- length(nall)
    filefa <-  paste(repmcmc,"ancestral.frequencies.txt",sep="")
    fa <- as.matrix(read.table(filefa))
    X11()
    plot(fa[,(iloc-1)*max(nall)+iall],
         xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""),
         ylab=paste("Frequency of allele",
           iall,"at locus",iloc),type="l")
    title(main=ifelse(iall==1,
            paste("Allele freq. in ancestral pop. at locus",iloc),
            ""))
    if(printit==T)
      {
        postscript(file=paste(path,"freq.ancestral.pop.loc",iloc,".ps",sep=""))
        plot(fa[,(iloc-1)*max(nall)+iall],
             xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""),
             ylab=paste("Frequency of allele",
               iall,"at locus",iloc),type="l")
        title(main=ifelse(iall==1,
                paste("Allele freq. in ancestral pop. at locus",iloc),
                ""))
        dev.off()
      }

  }

