"show.simIBD" <-
function(dataset,
                     plot.coord=FALSE,
                       plot.tess=FALSE,
                       plot.freq.grid=FALSE,
                       plot.freq.indiv=FALSE,
                       loc.grid=1:nloc,
                       all.grid=dataset$allele.numbers,
                        zlim.freq=c(0,1),
                       plot.gen=FALSE)
  {
    setplot(seq(dataset$coord.lim[1],dataset$coord.lim[2],(dataset$coord.lim[2]-dataset$coord.lim[1])/10),
            seq(dataset$coord.lim[3],dataset$coord.lim[4],(dataset$coord.lim[4]-dataset$coord.lim[3])/10))
                
    nindiv <- nrow(dataset$genotypes)
    nloc <- length(dataset$allele.numbers)
   if(plot.coord==TRUE)
      { 
        #X11()
        plot(dataset$coord.indiv[1,],dataset$coord.indiv[2,],
             xlab="x dataset$coordinates",ylab="y dataset$coordinates")
        points(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],col=2)
        text(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],1:dataset$number.nuclei,pos=1,
             col=2,pch=2,cex=1.2)
        text(dataset$coord.indiv[1,],dataset$coord.indiv[2,],
             dataset$color.nuclei[dataset$nearest.nucleus.indiv],pos=2)
        title(main="Location of individuals")
      }
   
   if(plot.tess==TRUE)
     {
       #X11()
       image(seq(from=dataset$coord.lim[1],to=dataset$coord.lim[2],length=dataset$npix[1]),
             seq(from=dataset$coord.lim[3],to=dataset$coord.lim[4],length=dataset$npix[2]),
             matrix(nr=dataset$npix[1],nc=dataset$npix[2],
                    dataset$color.nuclei[dataset$nearest.nucleus.grid],
                    byrow=F),
             xlab="",ylab="",col=terrain.colors(dataset$npop))
       points(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],col=1,pch=".",cex=1.5)
       #text(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],dataset$color.nuclei,pos=2,cex=2)
       title(sub="Partition of the domain according to population membership")
      }
   
    if(plot.freq.grid==TRUE)
      {
        for(iloc in loc.grid)
          {
            for(iall in 1:all.grid[iloc])
              {
                #X11()
                FF <- rep(-999,prod(dataset$npix))
                for(ipop in 1:dataset$npop)
                  {
                    ff <- dataset$freq.grid[ipop,,iloc,iall]
                    cc <- dataset$color.nuclei[dataset$nearest.nucleus.grid]
                    FF[cc==ipop] <- ff[cc==ipop]
                  }
                image(seq(from=dataset$coord.lim[1],to=dataset$coord.lim[2],length=dataset$npix[1]),
                      seq(from=dataset$coord.lim[3],to=dataset$coord.lim[4],length=dataset$npix[2]),
                      matrix(nr=dataset$npix[1],nc=dataset$npix[2],FF,byrow=F),
                      col=heat.colors(500),
                      xlab="",ylab="",
                      #breaks=seq(min(dataset$freq.grid),max(dataset$freq.grid),length=501),
                      zlim=zlim.freq)
               ##  contour(seq(from=dataset$coord.lim[1],to=dataset$coord.lim[2],length=dataset$npix[1]),
##                       seq(from=dataset$coord.lim[3],to=dataset$coord.lim[4],length=dataset$npix[2]),
##                       matrix(nr=dataset$npix[2],nc=dataset$npix[1],FF,byrow=F),add=T)
                
                #points(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],pch=".")
                #text(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],dataset$color.nuclei,cex=2)
                title(paste("Frequencies of allele #",iall,"at locus #",iloc))
              }
          }
      }
    if(plot.freq.indiv==TRUE)
      {
        for(iloc in 1: loc.grid)
          {
            for(iall in 1:all.grid[iloc])
              {
                FF <- rep(NA,nindiv)
                for(ipop in 1:dataset$npop)
                  {
                    ff <- dataset$freq.indiv[ipop,,iloc,iall]
                    cc <- dataset$color.nuclei[nearest.nucleus.indiv]
                    FF[cc==ipop] <- ff[cc==ipop]
                  }
                                        #X11()
                look <- as.image(x=t(dataset$coord.indiv),Z=FF)
                image.plot(look,main=paste("Field of frequencies for locus #",
                                  iloc,"allele #",iall))
                points(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],col=2,cex=2,lwd=3)
                                        #text(dataset$coord.nuclei[1,],dataset$coord.nuclei[2,],dataset$color.nuclei,col=2,pos=1)
                
                                        #text(dataset$coord.indiv[1,],dataset$coord.indiv[2,],dataset$genotypes[,2*(iloc-1)+iall],
                                        #     cex=1.5,lwd=2)
              }
          }
      }

   
   if(plot.gen==TRUE)
     {
       nindiv <- ncol(dataset$coord.indiv)
       for(iloc in 1:loc.grid)
         {
           #X11()
           if(dataset$allele.numbers[iloc]==2){
             aa <- ((dataset$genotypes[,2*(iloc-1)+1]==dataset$genotypes[,2*(iloc-1)+2])&
                    dataset$genotypes[,2*(iloc-1)+1]==1)
             AA <- ((dataset$genotypes[,2*(iloc-1)+1]==dataset$genotypes[,2*(iloc-1)+2])&
                    dataset$genotypes[,2*(iloc-1)+1]==2)
             aA <- !aa & !AA
             print(dataset$coord.indiv[1,])
             print(dataset$coord.indiv[2,])
             plot(dataset$coord.indiv[1,],dataset$coord.indiv[2,],type="n",cex=1,lwd=1,
                  xlab="",ylab="")
             title(paste("Genotypes at locus #",iloc))
             points(dataset$coord.indiv[1,aa],dataset$coord.indiv[2,aa],pch=15,col="red")
             points(dataset$coord.indiv[1,AA],dataset$coord.indiv[2,AA],pch=16,col="cyan")
             points(dataset$coord.indiv[1,aA],dataset$coord.indiv[2,aA],pch=17,col="green")
                    
           }else
           {
             text(dataset$coord.indiv[1,],dataset$coord.indiv[2,],dataset$genotypes[,2*(iloc-1)+1],
                  cex=1,lwd=2,pos=2,adj=0.1)
             text(dataset$coord.indiv[1,],dataset$coord.indiv[2,],dataset$genotypes[,2*(iloc-1)+2],
                  cex=1,lwd=2,pos=4,col=2,adj=0.1)
           }
         }
     }
  }

