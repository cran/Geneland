"PlotDrift" <-
function(repmcmc,printit=FALSE,file=NULL)
  {
    fileparam <- paste(repmcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    stepw <- as.numeric(param[12,3])
    
    filedrift <-  paste(repmcmc,"drifts.txt",sep="")
    drift <- as.matrix(read.table(filedrift))
    nclassmax <- ncol(drift)
    if(printit == TRUE){
      postscript(file=file)
      par(mfrow=c(nclassmax,1))
      for(iclass in 1:nclassmax)
        {
          plot(drift[,iclass],
               xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""),
               ylab=paste("Drift of population",iclass))
        }
      dev.off()
    }
    X11()
    par(mfrow=c(nclassmax,1))
    for(iclass in 1:nclassmax)
      {
        plot(drift[,iclass],
             xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""),
             ylab=paste("Drift of population",iclass))
      }
  }

