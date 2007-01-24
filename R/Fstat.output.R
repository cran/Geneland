`Fstat.output` <-
function(genotypes,allele.numbers,path.mcmc)
  {
    fileparam <- paste(path.mcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    filepm <- paste(path.mcmc,"posterior.mode.txt",sep="")

    nindiv <- nrow(genotypes)
    npop <- as.numeric(param[param[,1]=="npopmax",3])
    map <- scan(filepm)

     if(sum(is.na(genotypes)) != 0)  warning("Genotypes contain missing values which might bias computations")


                                         # Pairwise computations
     Fis=matrix(nr=npop,nc=npop,-999)
     Fst=matrix(nr=npop,nc=npop,-999)
     Fit=matrix(nr=npop,nc=npop,-999)

    for(iclass1 in 1:(npop-1))
      {
        for(iclass2 in (iclass1+1):npop)
          {
            sub1 <- map==iclass1
            sub2 <- map==iclass2
            if((sum(sub1)!=0)  & (sum(sub2)!=0)){
              ztmp <- genotypes[sub1 | sub2,]
              nindivtmp <- nrow(ztmp)
              maptmp <- map[sub1 | sub2]
              maptmp[maptmp==iclass1] <- 1
              maptmp[maptmp==iclass2] <- 2
              tabindiv <- matrix(nr=nindivtmp,nc=2,data=-999)
              kk <- numeric(2)
              effcl <- table(maptmp)
              nb.nuclei.max <- nindivtmp
              nloc <- length(allele.numbers)
              nloc2 <- 2*nloc
              Fistmp <- Fsttmp <- Fittmp <- -999
              res<- .Fortran(name="fstae",
                             PACKAGE="Geneland",
                             as.integer(nindivtmp),
                             as.integer(nb.nuclei.max),
                             as.integer(nloc),
                             as.integer(nloc2),
                             as.integer(allele.numbers),
                             as.integer(2),
                             as.integer(effcl),
                             as.integer(ztmp),
                             as.integer(maptmp),
                             as.integer(tabindiv),
                             as.integer(kk),
                             as.single(Fistmp),
                             as.single(Fsttmp),
                             as.single(Fittmp))
              Fis[iclass1,iclass2] <- res[[12]][1]
              Fst[iclass1,iclass2] <- res[[13]][1]
              Fit[iclass1,iclass2] <- res[[14]][1]
            }
          }
      }
    #list(Total.Fis=res.glob[[12]][1],
         #Total.Fst=res.glob[[13]][1],
         #Total.Fit=res.glob[[14]][1],
    list(Pairwise.Fis=Fis,
         Pairwise.Fst=Fst,
         Pairwise.Fit=Fit)
  }

