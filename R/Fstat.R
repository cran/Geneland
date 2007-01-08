"Fstat" <-
function(genotypes,allele.numbers,npop,pop.mbrship)
  {
    #dyn.load("~/mathlib/Fst/myfst.so")

    #print("debut fonction R Fstat ")
    #print(genotypes)
    #print(allele.numbers)
    #print(npop)
    #print(pop.mbrship)
     
    nindiv <- nrow(genotypes)
    
                                        # Pairwise computations
    Fis=matrix(nr=npop,nc=npop,-999)
    Fst=matrix(nr=npop,nc=npop,-999)
    Fit=matrix(nr=npop,nc=npop,-999)

    if(npop > 1){
    for(iclass1 in 1:(npop-1))
      {
        for(iclass2 in (iclass1+1):npop)
          {
            sub1 <- pop.mbrship==iclass1
            sub2 <- pop.mbrship==iclass2
            if((sum(sub1)!=0)  & (sum(sub2)!=0))
              {
                                        #print("coucou")
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
                res<- .Fortran(name="fstae",
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
                Fis[iclass1,iclass2] <- res[[11]][1]
                Fst[iclass1,iclass2] <- res[[12]][1]
                Fit[iclass1,iclass2] <- res[[13]][1]
              }
          }
      }
  }
  
    if(npop == 1)
      {
         #sub1 <- pop.mbrship==1
         ztmp <- rbind(genotypes,genotypes)
         nindivtmp <- nrow(ztmp)
         pop.mbrshiptmp <- c(rep(1,nindivtmp/2),rep(2,nindivtmp/2))
         #pop.mbrshiptmp[pop.mbrshiptmp==iclass1] <- 1
         #pop.mbrshiptmp[pop.mbrshiptmp==iclass2] <- 2
         tabindiv <- matrix(nr=nindivtmp,nc=2,data=-999)
         kk <- numeric(2)
         effcl <- table(pop.mbrshiptmp)
         nloc <- length(allele.numbers)
         nloc2 <- 2*nloc
         Fistmp <- Fsttmp <- Fittmp <- -999
##         print("Calling Fortran function fstae")
         res<- .Fortran(name="fstae",
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
         Fis[1,1] <- res[[11]][1]
         Fst[1,1] <- res[[12]][1]
         Fit[1,1] <- res[[13]][1]
       }
    list(Fis=Fis,
         Fst=Fst,
         Fit=Fit)
  }

