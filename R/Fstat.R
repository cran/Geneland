"Fstat" <-
function(repdat,repmcmc)
  {
    fileparam <- paste(repmcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    filez <- paste(repdat,"genotypes.txt",sep="")
    filenall <- paste(repdat,"allele.numbers.txt",sep="")
    filepm <- paste(repmcmc,"posterior.mode.txt",sep="")

    z <- as.matrix(read.table(filez))
    nindiv <- nrow(z)
    nall <- scan(filenall)
    nclass <- as.numeric(param[param[,1]=="nclassmax",3])
    map <- scan(filepm)


#                                         # Global computation
#     tabindiv <- matrix(nr=nindiv,nc=nclass,data=-999)
#     kk <- numeric(nclass)
#     effcl <- table(map)
#     nppmax <- nindiv
#     nloc <- length(nall)
#     nloc2 <- 2*nloc
#     Total.Fis <- Total.Fst <- Total.Fit <- -999
#     res.glob <- .Fortran(name="fstat",
#                          PACKAGE="Geneland",
#                          as.integer(nindiv),
#                          as.integer(nppmax),
#                          as.integer(nloc),
#                          as.integer(nloc2),
#                          as.integer(nall),
#                          as.integer(nclass),
#                          as.integer(effcl),
#                          as.integer(z),
#                          as.integer(map),
#                          as.integer(tabindiv),
#                          as.integer(kk),
#                          as.single(Total.Fis),
#                          as.single(Total.Fst),
#                          as.single(Total.Fit))


                                         # Pairwise computations
     Fis=matrix(nr=nclass,nc=nclass,-999)
     Fst=matrix(nr=nclass,nc=nclass,-999)
     Fit=matrix(nr=nclass,nc=nclass,-999)

    for(iclass1 in 1:(nclass-1))
      {
        for(iclass2 in (iclass1+1):nclass)
          {
            sub1 <- map==iclass1
            sub2 <- map==iclass2
            if((sum(sub1)!=0)  & (sum(sub2)!=0)){
              ztmp <- z[sub1 | sub2,]
              nindivtmp <- nrow(ztmp)
              maptmp <- map[sub1 | sub2]
              maptmp[maptmp==iclass1] <- 1
              maptmp[maptmp==iclass2] <- 2
              tabindiv <- matrix(nr=nindivtmp,nc=2,data=-999)
              kk <- numeric(2)
              effcl <- table(maptmp)
              nppmax <- nindivtmp
              nloc <- length(nall)
              nloc2 <- 2*nloc
              Fistmp <- Fsttmp <- Fittmp <- -999
              res<- .Fortran(name="fstat",
                             PACKAGE="Geneland",
                             as.integer(nindivtmp),
                             as.integer(nppmax),
                             as.integer(nloc),
                             as.integer(nloc2),
                             as.integer(nall),
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

