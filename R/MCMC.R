MCMC <-
function (coordinates = NULL, genotypes, ploidy = 2, dominance = "Codominant", 
    allele.numbers, path.mcmc, rate.max, delta.coord = 0, shape1 = 2, 
    shape2 = 20, npopmin = 1, npopinit, npopmax, nb.nuclei.max, 
    nit, thinning = 1, freq.model = "Uncorrelated", varnpop = TRUE, 
    spatial = TRUE, jcf = TRUE, filter.null.alleles = TRUE, prop.update.cell = 0.1, 
    write.rate.Poisson.process = FALSE, write.number.nuclei = TRUE, 
    write.number.pop = TRUE, write.coord.nuclei = TRUE, write.color.nuclei = TRUE, 
    write.freq = TRUE, write.ancestral.freq = TRUE, write.drifts = TRUE, 
    write.logposterior = TRUE, write.loglikelihood = TRUE, write.true.coord = TRUE, 
    write.size.pop = FALSE, miss.loc = NULL) 
{
    if (substring(path.mcmc, first = nchar(path.mcmc), last = nchar(path.mcmc)) != 
        "/") 
        stop(" path.mcmc has to end with /")
    if ((nit%%thinning) != 0) 
        stop("nit/thinning is not an integer")
    if (missing(npopinit)) 
        npopinit <- npopmax
    if (npopinit > npopmax) 
        stop("npopinit > npopmax")
    if (npopinit < npopmin) 
        stop("npopinit < npopmin")
    if ((freq.model != "Correlated") & (freq.model != "Uncorrelated")) {
        stop(paste("Error:", freq.model, "is not a frequency model. Check spelling (case sensitive) "))
    }
    if (is.null(coordinates)) {
        if (spatial) {
            stop("Please give spatial coordinates of individuals or set argument spatial to FALSE")
        }
        else {
            nindiv <- nrow(genotypes)
            n.int <- ceiling(sqrt(nindiv))
            x <- rep(seq(from = 0, to = 1, length = n.int), n.int)
            y <- rep(seq(from = 0, to = 1, length = n.int), n.int)
            y <- as.vector(t(matrix(nr = n.int, nc = n.int, y, 
                byrow = FALSE)))
            coordinates <- cbind(x, y)[1:nindiv, ]
        }
    }
    else {
        if (ncol(coordinates) != 2) 
            stop("matrix of coordinates does not have 2 columns")
    }
    if (ploidy == 1 & filter.null.alleles) 
        stop("Null alleles can not be filtered with haploid data")
    if (ploidy == 1 & dominance == "Dominant") 
        stop("Option ploidy=1 (haploid) and dominance='Dominant' not compatible")
    if (dominance != "Codominant") {
        if (dominance != "Dominant") {
            print(paste(dominance, " is not a valid argument."))
            print("See on-line help for details and check spelling (case sensitive)")
            stop()
        }
    }
    nindiv <- nrow(genotypes)
    if (is.null(miss.loc)) 
        miss.loc <- matrix(nr = nindiv, nc = ifelse(ploidy == 
            2, ncol(genotypes)/2, ncol(genotypes)), data = 0)
    if (missing(rate.max)) 
        rate.max <- nrow(genotypes)
    if (missing(nb.nuclei.max)) {
        nb.nuclei.max <- ifelse(spatial, 2 * nrow(genotypes), 
            nrow(genotypes))
    }
    if (!spatial & (nb.nuclei.max < nrow(genotypes))) {
        stop("With the option spatial=FALSE, nb.nuclei.max should be at least equal to the number of individuals")
    }
    if (spatial & (nb.nuclei.max < 2 * rate.max)) {
        stop("nb.nuclei.max is too small as compared to rate.max")
    }
    output.files <- c(write.rate.Poisson.process, write.number.nuclei, 
        write.number.pop, write.coord.nuclei, write.color.nuclei, 
        write.freq, write.ancestral.freq, write.drifts, write.logposterior, 
        write.loglikelihood, write.true.coord, write.size.pop)
    if (!is.matrix(coordinates)) 
        coordinates <- as.matrix(coordinates)
    if (!is.matrix(genotypes)) 
        genotypes <- as.matrix(genotypes)
    if (ploidy == 1 | dominance == "Dominant") {
        data.tmp <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) * 
            2)
        data.tmp[, seq(1, ncol(genotypes) * 2 - 1, 2)] <- genotypes
        data.tmp[, seq(2, ncol(genotypes) * 2, 2)] <- genotypes
        genotypes <- data.tmp
    }
    data.tmp <- FormatGenotypes(as.matrix(genotypes))
    genotypes <- data.tmp$genotypes
    if (missing(allele.numbers)) {
        allele.numbers <- data.tmp$allele.numbers
    }
    print("In R function MCMC, the number of observed alleles are:")
    print(allele.numbers)
    if (sum(allele.numbers == 1) > 0) {
        stop("Some of the markers do not display any polymorphism")
    }
    print(paste("There are ", sum(is.na(genotypes)), "missing values"))
    nchar.path <- nchar(path.mcmc)
    nindiv <- nrow(genotypes)
    nloc <- length(allele.numbers)
    if (filter.null.alleles & dominance == "Dominant") {
        filter.null.alleles <- FALSE
        print("Dominant markers -> Option regarding null allele set to FALSE")
    }
    if (filter.null.alleles & dominance == "Codominant") {
        allele.numbers <- allele.numbers + 1
    }
    if (dominance == "Dominant") 
        dom <- 1
    if (dominance == "Codominant") 
        dom <- 0
    nalmax <- max(allele.numbers)
    npp <- -999
    u <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    utemp <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    ctemp <- rep(times = nb.nuclei.max, -999)
    t <- matrix(nr = 2, nc = nindiv, data = -999)
    ttemp <- matrix(nr = 2, nc = nindiv, data = -999)
    f <- array(dim = c(npopmax, nloc, nalmax), data = -999)
    ftemp <- array(dim = c(npopmax, nloc, nalmax), data = -999)
    fa <- array(dim = c(nloc, nalmax), data = -999)
    drift <- rep(-999, npopmax)
    drifttemp <- rep(-999, npopmax)
    indcell <- rep(times = nindiv, -999)
    indcelltemp <- rep(times = nindiv, -999)
    distcell <- rep(times = nindiv, -999)
    distcelltemp <- rep(times = nindiv, -999)
    xlim <- ylim <- rep(-999, times = 2)
    n <- array(dim = c(npopmax, nloc, nalmax), data = -999)
    ntemp <- array(dim = c(npopmax, nloc, nalmax), data = -999)
    a <- rep(times = nalmax, -999)
    ptemp <- rep(times = nalmax, -999)
    effcl <- rep(times = npopmax, -999)
    iclv <- rep(times = npopmax, -999)
    cellclass <- rep(times = nb.nuclei.max, -999)
    listcell <- rep(times = nb.nuclei.max, -999)
    fmodel <- ifelse(freq.model == "Correlated", 1, 0)
    kfix <- 1 - as.integer(varnpop)
    full.cond.y <- matrix(nr = nalmax, nc = 2, 0)
    geno.999 <- genotypes
    geno.999[is.na(genotypes)] <- -999
    true.genotypes <- geno.999
    out.res <- .Fortran(name = "mcmcgld", PACKAGE = "Geneland", 
        as.double(t(coordinates)), as.integer(geno.999), as.integer(allele.numbers), 
        as.integer(ploidy), as.integer(dom), as.character(path.mcmc), 
        as.integer(nchar.path), as.double(rate.max), as.double(delta.coord), 
        as.double(shape1), as.double(shape2), as.integer(nit), 
        as.integer(thinning), as.integer(filter.null.alleles), 
        as.integer(nindiv), as.integer(nloc), as.integer(nloc * 
            2), as.integer(nalmax), as.integer(npp), as.integer(nb.nuclei.max), 
        as.integer(npopinit), as.integer(npopmin), as.integer(npopmax), 
        as.double(t), as.double(ttemp), as.double(u), as.double(utemp), 
        as.integer(c), as.integer(ctemp), as.double(f), as.double(ftemp), 
        as.double(fa), as.double(drift), as.double(drifttemp), 
        as.integer(indcell), as.integer(indcelltemp), as.double(distcell), 
        as.double(distcelltemp), as.double(xlim), as.double(ylim), 
        as.integer(n), as.integer(ntemp), as.double(a), as.double(ptemp), 
        as.integer(cellclass), as.integer(listcell), as.integer(fmodel), 
        as.integer(kfix), as.integer(spatial), as.integer(jcf), 
        as.integer(true.genotypes), as.double(full.cond.y), as.integer(output.files), 
        as.double(prop.update.cell), as.integer(miss.loc))
    param <- c(paste("ploidy :", ploidy), paste("dominance :", 
        dominance), paste("rate.max :", rate.max), paste("delta.coord :", 
        delta.coord), paste("npopmin :", npopmin), paste("npopinit :", 
        npopinit), paste("npopmax :", npopmax), paste("nindiv :", 
        nindiv), paste("nloc :", nloc), paste("nb.nuclei.max :", 
        nb.nuclei.max), paste("nit :", nit), paste("thinning :", 
        thinning), paste("freq.model :", freq.model), paste("varnpop :", 
        varnpop), paste("filter.null.alleles :", filter.null.alleles), 
        paste("spatial :", spatial))
    write.table(param, file = paste(path.mcmc, "parameters.txt", 
        sep = ""), quote = FALSE, row.name = FALSE, col.name = FALSE)
}
