"FormatGenotypes" <-
function(genotypes)
  {
    formatted <- genotypes
                                        # vector whose entries
                                        # will be the numbers
                                        # of allele per locus
    nall <- numeric(ncol(genotypes)/2)
    
    for(ic in 1:(ncol(genotypes)/2))
      {
        print(ic)
        ens <- sort(unique(c(genotypes[,2*ic-1],genotypes[,2*ic])))
        max.ens <- length(ens)
        for(il in 1:(dim(genotypes)[1]))
          {
            formatted[il,2*ic-1] <- ifelse(is.na(genotypes[il,2*ic-1]),
                                           NA,
                                           (1:max.ens)[genotypes[il,2*ic-1]==ens])
            
            formatted[il,2*ic]   <- ifelse(is.na(genotypes[il,2*ic]),
                                           NA,
                                           (1:max.ens)[genotypes[il,2*ic]==ens])
          }
        nall[ic] <- max.ens
      }
    genotypes <- as.matrix(genotypes)
    genotypes[is.na(genotypes)] <- -999
    list(genotypes=formatted,allele.numbers=nall)
  }

