"PosteriorMode" <-
function(repdat,repmcmc,write=FALSE,plotit=TRUE,
                          printit=FALSE,file=NULL)
  {
                                        # get informations about the MCMC run 
    fileparam <- paste(repmcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    dt <-  as.numeric(param[param[,1]=="dt",3])
    nclassmax <-  as.numeric(param[param[,1]=="nclassmax",3])

    param.postprocess <- as.matrix(read.table(paste(repmcmc,"postprocess.parameters.txt",sep="")))
    nxdom <- as.numeric(param.postprocess[1,3])
    nydom <- as.numeric(param.postprocess[2,3])

    files <- paste(repdat,"coordinates.txt",sep="/")
    s <- as.matrix(read.table(files))
    filedom <- paste(repmcmc,"proba.pop.membership.txt",sep="/")
    dom.post <- as.matrix(read.table(filedom))
    
    
    s[,1] <- s[,1] - min(s[,1])
    s[,2] <- s[,2] - min(s[,2])
                                        # rounding sites coordinates
    xlim <- c(min(s[,1])-dt/2,max(s[,1])+dt/2)
    ylim <- c(min(s[,2])-dt/2,max(s[,2])+dt/2)
    Dx <- (xlim[2]-xlim[1])/(nxdom-1)
    Dy <- (ylim[2]-ylim[1])/(nydom-1)
    s.discr <- s
    s.discr[,1] <- floor(s[,1]/Dx)*Dx
    s.discr[,2] <- floor(s[,2]/Dy)*Dy
    
    is <- (s.discr[,1]-xlim[1])/Dx + 1
    js <- (s.discr[,2]-ylim[1])/Dy + 1

    ks <- (is-1)*nydom + js
    map <- numeric(length(ks))
    for(k in 1:length(ks))
      {map[k] <- order(dom.post[ks[k],],decreasing=TRUE)[1]}
    
    filepm <- paste(repmcmc,"posterior.mode.txt",sep="/")
    if(write) write.table(map,file=filepm,quote=FALSE,row.name=FALSE,col.name=FALSE)

    map.dom <- t(apply(dom.post,1,order))[,nclassmax]
    
    if(plotit) {
      X11()
      setplot(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
              seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom))
      frame <- max(max(s[,1])-min(s[,1]),max(s[,2])-min(s[,2]))/40
      image(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
            seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
            matrix(map.dom,nr=nxdom,nc=nydom,byrow=TRUE),
            xlab="x coordinates ",ylab="y coordinates",
            main="",cex=1.5,cex.lab=1.5,col=terrain.colors(nclassmax),
            xlim=c(min(s[,1]-dt/2-frame),max(s[,1]+dt/2+frame)),
            ylim=c(min(s[,2]-dt/2-frame),max(s[,2]+dt/2+frame)),
            sub="Map of posterior mode of population memmbership")
    }
    if(printit){
      postscript(file)
      setplot(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
              seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom))
      frame <- max(max(s[,1])-min(s[,1]),max(s[,2])-min(s[,2]))/40
      image(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
            seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
            matrix(map.dom,nr=nxdom,nc=nydom,byrow=TRUE),
            xlab="x coordinates ",ylab="y coordinates",
            main="",cex=1.5,cex.lab=1.5,col=terrain.colors(nclassmax),
            xlim=c(min(s[,1]-dt/2-frame),max(s[,1]+dt/2+frame)),
            ylim=c(min(s[,2]-dt/2-frame),max(s[,2]+dt/2+frame)),
            sub="Map of posterior mode of population memmbership")
      dev.off()
    }
  }

