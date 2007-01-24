`Plotntile` <-
function(path.mcmc,printit=FALSE,file)
  {
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])
    
    filenpp <- paste(path.mcmc,"nuclei.numbers.txt",sep="")
    npp <- scan(filenpp)
    if(printit==T){
      postscript(file)
      par(mfrow=c(1,2))
      plot(npp,type="l",
           ylab="Number of Tiles",
           xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""))
      hist(npp,plot=T,prob=T,xlab="Number of tiles along the chain")
      dev.off()
    }
    else{
      par(mfrow=c(1,2))
      plot(npp,type="l",ylab="Number of Tiles",
           xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""))
      hist(npp,plot=T,prob=T,xlab="Number of tiles along the chain")
    }
  }

