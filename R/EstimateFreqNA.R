`EstimateFreqNA` <-
function (genotypes, path.mcmc) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "")
    param <- as.matrix(read.table(fileparam))
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    ploidy <- as.numeric(param[param[, 1] == "ploidy", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    if (!filter.null.alleles) 
        stop("The null alleles option has not been selected for MCMC inference")
    if (ploidy == 1) 
        stop("Null alleles can not be filtered with haploid data")
    allele.numbers <- FormatGenotypes(genotypes)$allele.numbers
    allele.numbers <- allele.numbers + 1
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    f <- as.matrix(read.table(filef))
    npopmax <- ncol(f)
    nloc <- length(allele.numbers)
    res <- matrix(nrow = nloc, ncol = 2, data = 0)
    for (ipop in 1:npopmax) {
        for (iloc in 1:nloc) {
            iall <- allele.numbers[iloc]
            sub1 <- rep(FALSE, times = iall - 1)
            sub2 <- TRUE
            sub3 <- rep(FALSE, times = max(allele.numbers) - 
                iall)
            sub <- c(sub1, sub2, sub3)
            sub1 <- rep(FALSE, (iloc - 1) * max(allele.numbers))
            sub2 <- sub
            sub3 <- rep(FALSE, (nloc - iloc) * max(allele.numbers))
            sub <- rep(c(sub1, sub2, sub3), times = nit/thinning)
            ff <- f[sub, ipop]
            ff <- ff[ff > 0]
            res[iloc, 1] <- res[iloc, 1] + sum(ff)
            res[iloc, 2] <- res[iloc, 2] + length(ff)
        }
    }
    for (iloc in 1:nloc) {
        res[iloc, 1] <- res[iloc, 1]/res[iloc, 2]
    }
    rownames(res) <- paste("Locus", 1:nloc)
    res[, 1]
}
