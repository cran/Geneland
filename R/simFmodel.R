"simFmodel" <-
function(nindiv,
                      s=NULL,
                      slim,
                      npp,
                      u=NULL,
                      c=NULL,
                      nloc,
                      nall,
                      npop,
                      drift,
                      seed=NULL,
                      plots=TRUE,
                      ploth=TRUE,
                      write=FALSE,
                      repout=NULL)
  {
    if(!is.null(seed)) { set.seed(seed)}
                                        # path to output files
    files <- paste(repout,"coordinates.txt",sep="")
    filez <- paste(repout,"genotypes.txt",sep="")
    filenall <- paste(repout,"allele.numbers.txt",sep="")
    fileu <- paste(repout,"coord.nuclei.txt",sep="")
    filec <- paste(repout,"color.nuclei.txt",sep="")
    filef <- paste(repout,"frequencies.",sep="")
    filefa <- paste(repout,"ancestral.frequencies.txt",sep="")
    filedrift <- paste(repout,"drifts.txt",sep="")
    filenall <- paste(repout,"allele.numbers.txt",sep="")
    fileppvois <- paste(repout,"index.nearest.nucleus.txt",sep="")
    
                                        # spatial coord of indviduals
    if(is.null(s))
      {
        s <- rbind(runif(min=slim[1],max=slim[2],nindiv),
                   runif(min=slim[3],max=slim[4],nindiv))
      }
    
                                        # centers and colors of tiles
    if(is.null(u))
      {
        u <-  rbind(runif(min=0,max=1,npp),
                    runif(min=0,max=1,npp))
        c <- numeric(npp)
        for(ipp in 1:npp)c[ipp] <- rdiscr(rep(1/npop,npop))
                                        # avoid to have only one pop
                                        # assuming we simulate at least two pop
        if(npp==2) c <- 1:2
      } 
    ppvois <- numeric(nindiv)
    for(iindiv in 1:nindiv)
      {
        k <- 1
        dd <- 1e+9
        for(ipp in 1:npp)
          {
            ddnew <- (u[1,ipp]-s[1,iindiv])^2+(u[2,ipp]-s[2,iindiv])^2
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
        for(iindiv in (1:nindiv)[c[ppvois]==iclass] )
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
        plot(s[1,],s[2,],
             xlab="x coordinates",ylab="y coordinates")
        points(u[1,],u[2,],col=2)
        text(u[1,],u[2,],1:npp,pos=1,col=2,pch=2,cex=1.2)
        text(s[1,],s[2,],ppvois,pos=1)
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
                hist(c(z[c[ppvois]==iclass,2*(iloc-1)+1],
                       z[c[ppvois]==iclass,2*(iloc-1)+2]),
                     breaks=seq(.5,nall[iloc]+.5,1),
                     prob=TRUE,main=paste("Histogram of freq. in pop.",iclass,
                              ", locus",iloc),
                     xlab="",ylim=c(0,1),axes=F)
                points(1:nall[iloc],freq[iclass,iloc,1:nall[iloc]],
                       type="h",col=2)
              }
          }
      }
    if(write==TRUE)
      {
        write.table(floor(t(s)*1e+3)*1e-3,
                    files,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(z,
                    filez,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(nall,
                    filenall,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(floor(t(u)*1e+3)*1e-3,
                    fileu,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(c,
                    filec,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(nall,
                    filenall,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(fa,
                    filefa,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        write.table(floor(drift*1e+3)*1e-3,
                    filedrift,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
        for(iclass in 1:npop)
          {
            write.table(freq[iclass,,],
                        quote=FALSE,
                        col.names=FALSE,
                        row.names=FALSE,
                        paste(filef,"class",iclass,".txt",sep=""))
          }
         write.table(ppvois,
                    fileppvois,
                    quote=FALSE,
                    col.names=FALSE,
                    row.names=FALSE)
      }
  }

