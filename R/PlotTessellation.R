"PlotTessellation" <-
function(repdat,repmcmc,printit=FALSE,path=NULL)
  {

    param.postprocess <- as.matrix(read.table(paste(repmcmc,
                                                    "postprocess.parameters.txt",
                                                    sep="")))
    nxdom <- as.numeric(param.postprocess[1,3])
    nydom <- as.numeric(param.postprocess[2,3])

    param <- as.matrix(read.table(paste(repmcmc,"parameters.txt",sep="")))
    dt <- as.numeric(param[4,3])
    
    files <-  paste(repdat,"coordinates.txt",sep="")
    s <- as.matrix(read.table(files))
    filedom <- paste(repmcmc,"proba.pop.membership.txt",sep="")
    dom.post <- as.matrix(read.table(filedom))
    nclassmax <- ncol(dom.post)
    for(iclass in 1:nclassmax)
      {
        X11()
        setplot(s[,1],s[,2])
                                        #par(omi=c(0,0,0,0))
        image(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
              seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
              matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
              xlab="x coordinates",ylab="y coordinates",
              main=paste("Map of posterior probability to belong to class ",
                iclass))
        contour(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
                seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                add=TRUE)
        points(s[,1],s[,2],pch=16)
        if(printit==TRUE) {
          postscript(file=paste(path,"/mapproba.pop",iclass,".ps",sep=""))
          setplot(s[,1],s[,2])
                                        #par(omi=c(0,0,0,0))
          image(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
                seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                xlab="x coordinates",ylab="y coordinates",
                main=paste("Map of posterior probability to belong to class ",
                  iclass))
          contour(seq(min(s[,1]-dt/2),max(s[,1]+dt/2),length=nxdom),
                  seq(min(s[,2]-dt/2),max(s[,2]+dt/2),length=nydom),
                  matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                  add=TRUE)
          points(s[,1],s[,2],pch=16)
          dev.off()
        }
      }
  }

