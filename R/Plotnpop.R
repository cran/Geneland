"Plotnpop" <-
function(path.mcmc,printit=FALSE,file=NULL)
  {
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])
    
    filenpop <- paste(path.mcmc,"populations.numbers.txt",sep="")
    npop <- scan(filenpop)
    #layout(mat=matrix(nr=1,nc=2,data=1:2),width=c(4,2))
    if(printit==TRUE){
      postscript(file)
      par(mfrow=c(1,2))
      plot(npop,type="l",ylab="Number of classes",
           xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
           ylim=c(1,max(npop)+0.5))
      hist(npop,plot=TRUE,prob=TRUE,breaks=seq(.5,max(npop)+0.5,1),
           xlab="Number of populations along the chain")
    dev.off()
    }
    else{
      par(mfrow=c(1,2))
      plot(npop,type="l",ylab="Number of classes",xlab="Index of MCMC iterations",
           ylim=c(1,max(npop)+0.5))
      hist(npop,plot=TRUE,prob=TRUE,breaks=seq(.5,max(npop)+0.5,1),
           xlab="Number of populations along the chain")
    }
  }

