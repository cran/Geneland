"PlotFreq" <-
function(path.data,path.mcmc,ipop,iloc,iall,printit=FALSE,path=NULL)
  {
    cat(paste("Reading file of frequencies ",
                path.mcmc,
                "frequencies.txt",sep=""))
                                        # get informations about the MCMC run 
    fileparam <- paste(path.mcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    nit <- as.numeric(param[param[,1]=="nit",3])
    thinning <-  as.numeric(param[param[,1]=="thinning",3])

    filef <-  paste(path.mcmc,"frequencies.txt",sep="")
    f <- as.matrix(read.table(filef))
    npopmax <- ncol(f)
    nall <- scan(paste(path.data,"allele.numbers.txt",sep=""))
    nloc <- length(nall)
                                        # extract frequencies
                                        # from messy matrix f
    sub1 <- rep(FALSE,times=iall-1)
    sub2 <- TRUE
    sub3 <- rep(FALSE,times=max(nall)-iall)
    sub <- c(sub1,sub2,sub3)
    sub1 <- rep(FALSE,(iloc-1)*max(nall))
    sub2 <- sub
    sub3 <- rep(FALSE,(nloc-iloc)*max(nall))
    sub <- rep(c(sub1,sub2,sub3),times=nit/thinning)
    plot(f[sub,ipop],
         xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
         ylab=paste("Frequency of allele",
           iall,"at locus",iloc),type="l")
    title(main=ifelse(iall==1,
            paste("Allele frequencies in population",ipop,
                  "for locus",iloc),""))
              
    if(printit==TRUE)
      {
        postscript(file=paste(path,"freq.pop",ipop,".loc",iloc,".ps",sep=""))
        par(mfrow=c(ceiling(sqrt(nall[iloc])),ceiling(sqrt(nall[iloc]))))
                                        # extract frequencies
                                        # in this messy tabular
        sub1 <- rep(FALSE,times=iall-1)
        sub2 <- TRUE
        sub3 <- rep(FALSE,times=max(nall)-iall)
        sub <- c(sub1,sub2,sub3)
        sub1 <- rep(FALSE,(iloc-1)*max(nall))
        sub2 <- sub
        sub3 <- rep(FALSE,(nloc-iloc)*max(nall))
        sub <- rep(c(sub1,sub2,sub3),times=nit/thinning)
        plot(f[sub,ipop],
             xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
             ylab=paste("Frequency of allele",
               iall,"at locus",iloc),type="l")
        title(main=ifelse(iall==1,
                paste("Allele frequencies in population",ipop,
                      "for locus",iloc),""))
        dev.off()
      }
  }

