"simFmodel" <-
function(nindiv,
                      coordinates=NULL,
                      coord.lim,
                      number.nuclei,
                      coord.nuclei=NULL,
                      color.nuclei=NULL,
                      nloc,
                      nall,
                      npop,
                      drift,
                      seed=NULL,
                      plots=TRUE,
                      ploth=TRUE)
#  ,
#                      write=FALSE,
#                      repout=NULL)
  {
    if(!is.null(seed)) { set.seed(seed)}
    
                                        # spatial coord of indviduals
    if(is.null(coordinates))
      {
        coordinates <- rbind(runif(min=coord.lim[1],max=coord.lim[2],nindiv),
                   runif(min=coord.lim[3],max=coord.lim[4],nindiv))
      }
    
                                        # centers and colors of tiles
    if(is.null(coord.nuclei))
      {
        coord.nuclei <-  rbind(runif(min=0,max=1,number.nuclei),
                    runif(min=0,max=1,number.nuclei))
        color.nuclei <- numeric(number.nuclei)
        for(ipp in 1:number.nuclei)color.nuclei[ipp] <- rdiscr(rep(1/npop,npop))
                                        # avoid to have only one pop
                                        # assuming we simulate at least two pop
        if(number.nuclei==2) color.nuclei <- 1:2
      } 
    ppvois <- numeric(nindiv)
    for(iindiv in 1:nindiv)
      {
        k <- 1
        dd <- 1e+9
        for(ipp in 1:number.nuclei)
          {
            ddnew <- (coord.nuclei[1,ipp]-coordinates[1,iindiv])^2+(coord.nuclei[2,ipp]-coordinates[2,iindiv])^2
            if(ddnew < dd)
              {
                k <- ipp
                dd <- ddnew
              }
          }
        ppvois[iindiv] <- k
      }

                                        # alleles frequencies in ancestral population
    fa <- matrix(nr=nloc,nc=max(nall),data=-999)
    for(iloc in 1:nloc)
          {
            fa[iloc,1:nall[iloc]] <- rexp(n=nall[iloc])
            fa[iloc,1:nall[iloc]] <- fa[iloc,1:nall[iloc]] /
              sum(fa[iloc,1:nall[iloc]])
          }

    
                                        # alleles frequencies in present time population
    freq <- array(dim=c(npop,nloc,max(nall)),data=-999)
    for(iclass in 1:npop)
      {
        for(iloc in 1:nloc)
          {
            freq[iclass,iloc,1:nall[iloc]] <- rgamma(n=nall[iloc],
                                                   scale=1,
                                                   shape=fa[iloc,(1:nall[iloc])]*
                                                   (1-drift[iclass])/drift[iclass])
            freq[iclass,iloc,1:nall[iloc]] <- freq[iclass,iloc,1:nall[iloc]] /
              sum(freq[iclass,iloc,1:nall[iloc]])
          }
      }
    
                                        # impose condition on  freqs :  f111 > f211
    if(npop > 1)
      {
        if(freq[1,1,1] < freq[2,1,1])
          {
            tmp <- freq[2,,]
            freq[2,,] <- freq[1,,]
            freq[1,,] <- tmp
          }
      }
          
    
                                        # genotypes     
    z <- matrix(nr=nindiv,nc=nloc*2)
    for(iclass in 1:npop)
      {
        for(iindiv in (1:nindiv)[color.nuclei[ppvois]==iclass] )
          {
            for(iloc in 1:nloc)
              {
                z[iindiv,2*iloc-1] <- rdiscr(freq[iclass,iloc,1:nall[iloc]])
                z[iindiv,2*iloc]   <- rdiscr(freq[iclass,iloc,1:nall[iloc]])
              }
          }
      }
    
    if(plots==TRUE)
      { 
        X11()
        plot(coordinates[1,],coordinates[2,],
             xlab="x coordinates",ylab="y coordinates")
        points(coord.nuclei[1,],coord.nuclei[2,],col=2)
        text(coord.nuclei[1,],coord.nuclei[2,],1:number.nuclei,pos=1,col=2,pch=2,cex=1.2)
        text(coordinates[1,],coordinates[2,],ppvois,pos=1)
      }
    if(ploth==TRUE)
      {
        X11()
        par(mfrow=c(floor(sqrt(nloc)+1),
                  floor(sqrt(nloc))))
            #par(mfrow=c(1,nloc))
            for(iloc in 1:nloc)
              {
                plot(1:nall[iloc],fa[iloc,1:nall[iloc]],
                       type="h",col=2,xlab="",axes=FALSE,
                     sub=paste("Locus",iloc),ylim=c(0,1),
                     main="Frequencies in ancestral population")
              }
        #X11()
        #par(mfrow=c(npop,nloc))
        for(iclass in 1:npop)
          {
            X11()
            par(mfrow=c(floor(sqrt(nloc)+1),
                  floor(sqrt(nloc))))
            #par(mfrow=c(1,nloc))
            for(iloc in 1:nloc)
              {
                hist(c(z[color.nuclei[ppvois]==iclass,2*(iloc-1)+1],
                       z[color.nuclei[ppvois]==iclass,2*(iloc-1)+2]),
                     breaks=seq(.5,nall[iloc]+.5,1),
                     prob=TRUE,main=paste("Histogram of freq. in pop.",iclass,
                              ", locus",iloc),
                     xlab="",ylim=c(0,1),axes=F)
                points(1:nall[iloc],freq[iclass,iloc,1:nall[iloc]],
                       type="h",col=2)
              }
          }
      }
    list(coordinates=t(coordinates),
         genotypes=z,
         allele.numbers=nall,
         coord.nuclei=t(coord.nuclei),
         color.nuclei=c,
         frequencies=freq,
         ancestral.frequencies=fa,
         drifts=drift,
         index.nearest.nucleus=ppvois)
    
  }

