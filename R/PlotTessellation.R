"PlotTessellation" <-
function(coordinates,path.mcmc,printit=FALSE,path=NULL)
  {

    param.postprocess <- as.matrix(read.table(paste(path.mcmc,
                                                    "postprocess.parameters.txt",
                                                    sep="")))
    nxdom <- as.numeric(param.postprocess[1,3])
    nydom <- as.numeric(param.postprocess[2,3])

    param <- as.matrix(read.table(paste(path.mcmc,"parameters.txt",sep="")))
    delta.coord <- as.numeric(param[3,3])
    
    s <- coordinates
    filedom <- paste(path.mcmc,"proba.pop.membership.txt",sep="")
    dom.post <- as.matrix(read.table(filedom))
    npopmax <- ncol(dom.post)
    for(iclass in 1:npopmax)
      {
        X11()
        setplot(s[,1],s[,2])
                                        #par(omi=c(0,0,0,0))
        image(seq(min(s[,1]-delta.coord/2),max(s[,1]+delta.coord/2),length=nxdom),
              seq(min(s[,2]-delta.coord/2),max(s[,2]+delta.coord/2),length=nydom),
              matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
              xlab="x coordinates",ylab="y coordinates",
              main=paste("Map of posterior probability to belong to class ",
                iclass))
        contour(seq(min(s[,1]-delta.coord/2),max(s[,1]+delta.coord/2),length=nxdom),
                seq(min(s[,2]-delta.coord/2),max(s[,2]+delta.coord/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                add=TRUE)
        points(s[,1],s[,2],pch=16)
        if(printit==TRUE) {
          postscript(file=paste(path,"/mapproba.pop",iclass,".ps",sep=""))
          setplot(s[,1],s[,2])
                                        #par(omi=c(0,0,0,0))
          image(seq(min(s[,1]-delta.coord/2),max(s[,1]+delta.coord/2),length=nxdom),
                seq(min(s[,2]-delta.coord/2),max(s[,2]+delta.coord/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                xlab="x coordinates",ylab="y coordinates",
                main=paste("Map of posterior probability to belong to class ",
                  iclass))
          contour(seq(min(s[,1]-delta.coord/2),max(s[,1]+delta.coord/2),length=nxdom),
                  seq(min(s[,2]-delta.coord/2),max(s[,2]+delta.coord/2),length=nydom),
                  matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                  add=TRUE)
          points(s[,1],s[,2],pch=16)
          dev.off()
        }
      }
  }

