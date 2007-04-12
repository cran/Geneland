`Fstat.output` <-
function(genotypes,path.mcmc)
  {
    fileparam <- paste(path.mcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    filepm <- paste(path.mcmc,"modal.pop.indiv.txt",sep="")
    
    allele.numbers <- FormatGenotypes(genotypes)$allele.numbers

    nindiv <- nrow(genotypes)
    npop <- as.numeric(param[param[,1]=="npopmax",3])
    pop.mbrship <- scan(filepm)

     if(sum(is.na(genotypes)) != 0)  warning("Genotypes contain missing values which might bias computations")


                                         # Pairwise computations
     Fis=rep(-999,npop)
     Fst=matrix(nr=npop,nc=npop,-999)

    if(npop > 1)
      {
        for(iclass1 in 1:(npop-1))
          {
            for(iclass2 in (iclass1+1):npop)
              {
                sub1 <- pop.mbrship==iclass1
                sub2 <- pop.mbrship==iclass2
                if((sum(sub1)!=0)  & (sum(sub2)!=0))
                  {
                    print("coucou")
                    ztmp <- genotypes[sub1 | sub2,]
                    nindivtmp <- nrow(ztmp)
                    pop.mbrshiptmp <- pop.mbrship[sub1 | sub2]
                    pop.mbrshiptmp[pop.mbrshiptmp==iclass1] <- 1
                    pop.mbrshiptmp[pop.mbrshiptmp==iclass2] <- 2
                    tabindiv <- matrix(nr=nindivtmp,nc=2,data=-999)
                    kk <- numeric(2)
                    effcl <- table(pop.mbrshiptmp)
                    nloc <- length(allele.numbers)
                    nloc2 <- 2*nloc
                    Fistmp <- Fsttmp <- Fittmp <- -999
                    ## print("Calling Fortran function fstae")
                    res<- .Fortran(name="ggfst",
                                   PACKAGE="Geneland",
                                   as.integer(nindivtmp),
                                   as.integer(nloc),
                                   as.integer(nloc2),
                                   as.integer(allele.numbers),
                                   as.integer(2),
                                   as.integer(effcl),
                                   as.integer(ztmp),
                                   as.integer(pop.mbrshiptmp),
                                   as.integer(tabindiv),
                                   as.integer(kk),
                                   as.single(Fistmp),
                                   as.single(Fsttmp),
                                   as.single(Fittmp))
                    Fst[iclass1,iclass2] <- res[[12]][1]
                  }
              }
          }
      }
    for(iclass1 in 1:npop)
      {
        sub1 <- pop.mbrship==iclass1
        if((sum(sub1)!=0))
          {
                                        print("coucou")
            ztmp <-  rbind(genotypes[sub1,],genotypes[sub1,])
            nindivtmp <- nrow(ztmp)
            pop.mbrshiptmp <- c(rep(1,sum(sub1)),rep(2,sum(sub1)))
            tabindiv <- matrix(nr=nindivtmp,nc=2,data=-999)
            kk <- numeric(2)
            effcl <- table(pop.mbrshiptmp)
            nloc <- length(allele.numbers)
            nloc2 <- 2*nloc
            Fistmp <- Fsttmp <- Fittmp <- -999
            ## print("Calling Fortran function fstae")
            res<- .Fortran(name="ggfst",
                           PACKAGE="Geneland",
                           as.integer(nindivtmp),
                           as.integer(nloc),
                           as.integer(nloc2),
                           as.integer(allele.numbers),
                           as.integer(2),
                           as.integer(effcl),
                           as.integer(ztmp),
                           as.integer(pop.mbrshiptmp),
                           as.integer(tabindiv),
                           as.integer(kk),
                           as.single(Fistmp),
                           as.single(Fsttmp),
                           as.single(Fittmp))
            Fis[iclass1] <- res[[11]][1]
          }
      }
    list(Fis=Fis,
         Fst=Fst)
  }

