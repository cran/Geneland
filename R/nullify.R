"nullify" <-
function(genotypes,nall.null=1,nloc.null)
  {
    ## select ,nall.null allele at random among observed alleles at loci in 1:nloc.null
    ## and simulate a null allele
    ## return the empirical frequency of this allele
    z <- genotypes
    freq <- rep(0,nloc.null)
    for(iloc in 1:nloc.null)
      {
        for(iall in 1:nall.null)
            {
              zz <- z[,c(2*iloc-1,2*iloc)]
              ina <- sample(1:length(unique(c(zz[,1],zz[,2]))),2)
              freq[iloc] <- freq[iloc] + sum(zz==ina[iall])
              for(iindiv in 1:nrow(genotypes))
                {
                  if(sum(is.na(z[iindiv,c(2*iloc-1,2*iloc)]))==0)
                    {
                      if(sum(z[iindiv,c(2*iloc-1,2*iloc)]==ina[iall])==2)
                        {
                          zz[iindiv,] <- c(NA,NA)
                        }else{
                          if(z[iindiv,2*iloc-1]==ina[iall]) zz[iindiv,1] <- z[iindiv,2*iloc]
                          if(z[iindiv,2*iloc]==ina[iall]) zz[iindiv,2] <- z[iindiv,2*iloc-1]
                        }
                    }else{
                      if(sum(is.na(z[iindiv,c(2*iloc-1,2*iloc)]))==2)
                        {
                          zz[iindiv,] <- z[iindiv,]
                        }
                      if(is.na(z[iindiv,2*iloc]))
                         {
                           if(z[iindiv,2*iloc-1]==ina[iall]) zz[iindiv,1] <- NA
                         }
                      if(is.na(z[iindiv,2*iloc-1]))
                        {                         
                          if(z[iindiv,2*iloc]==ina[iall]) zz[iindiv,2] <- NA
                        }
                      }
                }
              z[,c(2*iloc-1,iloc)] <- zz
            }
        freq[iloc] <- freq[iloc]/nall.null
      }
    list(genotypes=z,freq=freq)
  }

