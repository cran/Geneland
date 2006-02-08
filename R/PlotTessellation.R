"PlotTessellation" <-
function(coordinates,path.mcmc,printit=FALSE,path=NULL)
  {
    coordinates <- as.matrix(coordinates)
      
    param.postprocess <- as.matrix(read.table(paste(path.mcmc,
                                                    "postprocess.parameters.txt",
                                                    sep="")))
    nxdom <- as.numeric(param.postprocess[1,3])
    nydom <- as.numeric(param.postprocess[2,3])

    param <- as.matrix(read.table(paste(path.mcmc,"parameters.txt",sep="")))
    delta.coord <- as.numeric(param[param[,1]=="delta.coord",3])
    
    filedom <- paste(path.mcmc,"proba.pop.membership.txt",sep="")
    dom.post <- as.matrix(read.table(filedom))
    npopmax <- ncol(dom.post)
    for(iclass in 1:npopmax)
      {
        get(getOption("device"))()
        setplot(coordinates[,1],coordinates[,2])
                                        #par(omi=c(0,0,0,0))
        image(seq(min(coordinates[,1]-delta.coord/2),max(coordinates[,1]+delta.coord/2),length=nxdom),
              seq(min(coordinates[,2]-delta.coord/2),max(coordinates[,2]+delta.coord/2),length=nydom),
              matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
              xlab="x coordinates",ylab="y coordinates")
        title(main=paste("Map of posterior probability to belong to class ",
                iclass))
        contour(seq(min(coordinates[,1]-delta.coord/2),max(coordinates[,1]+delta.coord/2),length=nxdom),
                seq(min(coordinates[,2]-delta.coord/2),max(coordinates[,2]+delta.coord/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                add=TRUE)
        points(coordinates[,1],coordinates[,2],pch=16)
        if(printit==TRUE) {
          postscript(file=paste(path,"/mapproba.pop",iclass,".ps",sep=""))
          setplot(coordinates[,1],coordinates[,2])
                                        #par(omi=c(0,0,0,0))
          image(seq(min(coordinates[,1]-delta.coord/2),max(coordinates[,1]+delta.coord/2),length=nxdom),
                seq(min(coordinates[,2]-delta.coord/2),max(coordinates[,2]+delta.coord/2),length=nydom),
                matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                xlab="x coordinates",ylab="y coordinates")
          title(main=paste("Map of posterior probability to belong to class ",
                  iclass))
          contour(seq(min(coordinates[,1]-delta.coord/2),max(coordinates[,1]+delta.coord/2),length=nxdom),
                  seq(min(coordinates[,2]-delta.coord/2),max(coordinates[,2]+delta.coord/2),length=nydom),
                  matrix(dom.post[,iclass],nr=nxdom,nc=nydom,byrow=TRUE),
                  add=TRUE)
          points(coordinates[,1],coordinates[,2],pch=16)
          dev.off()
        }
      }
  }

