"Plotnclass" <-
function(repmcmc,printit=FALSE,file=NULL)
  {
    fileparam <- paste(repmcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    stepw <- as.numeric(param[12,3])
    
    filenclass <- paste(repmcmc,"populations.numbers.txt",sep="")
    nclass <- scan(filenclass)
    #layout(mat=matrix(nr=1,nc=2,data=1:2),width=c(4,2))
    if(printit==TRUE){
      postscript(file)
      par(mfrow=c(1,2))
      plot(nclass,type="l",ylab="Number of classes",
           xlab=paste("Index of MCMC iteration"," (x ",stepw,")",sep=""),
           ylim=c(1,max(nclass)+0.5))
      hist(nclass,plot=TRUE,prob=TRUE,breaks=seq(.5,max(nclass)+0.5,1),
           xlab="Number of populations along the chain")
    dev.off()
    }
    else{
      par(mfrow=c(1,2))
      plot(nclass,type="l",ylab="Number of classes",xlab="Index of MCMC iterations",
           ylim=c(1,max(nclass)+0.5))
      hist(nclass,plot=TRUE,prob=TRUE,breaks=seq(.5,max(nclass)+0.5,1),
           xlab="Number of populations along the chain")
    }
  }

