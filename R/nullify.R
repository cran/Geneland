`nullify` <-
function (genotypes, nall.null = 1, nloc.null) 
{
    z <- genotypes
    for (iloc in 1:nloc.null) {
        ina <- sample(1:length(unique(as.vector(z[, c(2 * iloc - 
            1, 2 * iloc)]))), nall.null)
        for (iall in 1:nall.null) {
            zz <- z[, c(2 * iloc - 1, 2 * iloc)]
            for (iindiv in 1:nrow(genotypes)) {
                if (sum(is.na(z[iindiv, c(2 * iloc - 1, 2 * iloc)])) == 
                  0) {
                  if (sum(z[iindiv, c(2 * iloc - 1, 2 * iloc)] == 
                    ina[iall]) == 2) {
                    zz[iindiv, ] <- c(NA, NA)
                  }
                  else {
                    if (z[iindiv, 2 * iloc - 1] == ina[iall]) 
                      zz[iindiv, 1] <- z[iindiv, 2 * iloc]
                    if (z[iindiv, 2 * iloc] == ina[iall]) 
                      zz[iindiv, 2] <- z[iindiv, 2 * iloc - 1]
                  }
                }
                else {
                  if (sum(is.na(z[iindiv, c(2 * iloc - 1, 2 * 
                    iloc)])) == 2) {
                    zz[iindiv, ] <- z[iindiv, c(2 * iloc - 1, 
                      2 * iloc)]
                  }
                  else {
                    if (is.na(z[iindiv, 2 * iloc])) {
                      if (z[iindiv, 2 * iloc - 1] == ina[iall]) 
                        zz[iindiv, 1] <- NA
                    }
                    if (is.na(z[iindiv, 2 * iloc - 1])) {
                      if (z[iindiv, 2 * iloc] == ina[iall]) 
                        zz[iindiv, 2] <- NA
                    }
                  }
                }
            }
            z[, c(2 * iloc - 1, 2 * iloc)] <- zz
        }
    }
    list(genotypes = z)
}
