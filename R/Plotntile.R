"Plotntile" <-
function(repmcmc,printit=FALSE,file=NULL)
  {
    fileparam <- paste(repmcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    stepw <- as.numeric(param[12,3])
    
    filenpp <- paste(repmcmc,"nuclei.numbers.txt",sep="")
    npp <- scan(filenpp)
    if(printit==T){
      postscript(file)
      par(mfrow=c(1,2))
      plot(npp,type="l",
           ylab="Number of Tiles",
           xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""))
      hist(npp,plot=T,prob=T,xlab="Number of tiles along the chain")
      dev.off()
    }
    else{
      par(mfrow=c(1,2))
      plot(npp,type="l",ylab="Number of Tiles",
           xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""))
      hist(npp,plot=T,prob=T,xlab="Number of tiles along the chain")
    }
  }

