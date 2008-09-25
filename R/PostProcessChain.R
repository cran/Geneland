`PostProcessChain` <-
function (coordinates = NULL, genotypes, path.mcmc, nxdom, nydom, 
    burnin) 
{
    print("Reading MCMC parameter file")
    param <- as.matrix(read.table(paste(path.mcmc, "parameters.txt", 
        sep = "")))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    nb.nuclei.max <- as.numeric(param[param[, 1] == "nb.nuclei.max", 
        3])
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    ploidy <- as.numeric(param[param[, 1] == "ploidy", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    if (ploidy == 1) {
        data.tmp <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) * 
            2)
        data.tmp[, seq(1, ncol(genotypes) * 2 - 1, 2)] <- genotypes
        data.tmp[, seq(2, ncol(genotypes) * 2, 2)] <- genotypes
        genotypes <- data.tmp
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
    coordinates <- as.matrix(coordinates)
    data.tmp <- FormatGenotypes(as.matrix(genotypes))
    genotypes <- data.tmp$genotypes
    allele.numbers <- data.tmp$allele.numbers
    if (filter.null.alleles) {
        allele.numbers <- allele.numbers + 1
    }
    nindiv <- nrow(genotypes)
    nloc <- length(allele.numbers)
    filenpop <- paste(path.mcmc, "populations.numbers.txt", sep = "")
    filenpp <- paste(path.mcmc, "nuclei.numbers.txt", sep = "")
    fileu <- paste(path.mcmc, "coord.nuclei.txt", sep = "")
    filec <- paste(path.mcmc, "color.nuclei.txt", sep = "")
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    fileperm <- paste(path.mcmc, "perm.txt", sep = "")
    filedom <- paste(path.mcmc, "proba.pop.membership.txt", sep = "")
    filedomperm <- paste(path.mcmc, "proba.pop.membership.perm.txt", 
        sep = "")
    print("Estimating number of populations")
    npop <- scan(paste(path.mcmc, "populations.numbers.txt", 
        sep = ""))
    npop.est <- order(hist(npop[-(1:burnin)], breaks = seq(0.5, 
        npopmax + 0.5, 1), plot = FALSE)$counts, decreasing = TRUE)[1]
    lpd <- scan(paste(path.mcmc, "log.posterior.density.txt", 
        sep = ""))
    subK <- (npop == npop.est) & ((1:(nit/thinning)) > burnin)
    pivot <- order(lpd[subK], decreasing = TRUE)[1]
    dom <- matrix(nr = nxdom * nydom, nc = npopmax, data = 0)
    domperm <- matrix(nr = nxdom * nydom, nc = npopmax, data = 0)
    coorddom <- matrix(nr = 2, nc = nxdom * nydom, data = -999)
    indvois <- numeric(nxdom * nydom)
    distvois <- numeric(nxdom * nydom)
    orderf <- orderftmp <- 1:npopmax
    nalmax <- max(allele.numbers)
    u <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    xlim <- ylim <- rep(-999, 2)
    f <- fpiv <- array(dim = c(npopmax, nloc, nalmax), -999)
    ninrub = 0
    print("Calling Fortran function postprocesschain2")
    out.res <- .Fortran(name = "postprocesschain2", PACKAGE = "Geneland", 
        as.integer(nxdom), as.integer(nydom), as.integer(burnin), 
        as.integer(ninrub), as.integer(npopmax), as.integer(nb.nuclei.max), 
        as.integer(nindiv), as.integer(nloc), as.integer(allele.numbers), 
        as.integer(nalmax), as.double(xlim), as.double(ylim), 
        as.double(delta.coord), as.integer(nit), as.integer(thinning), 
        as.character(filenpop), as.character(filenpp), as.character(fileu), 
        as.character(filec), as.character(filef), as.character(fileperm), 
        as.character(filedom), as.double(t(coordinates)), as.double(u), 
        as.integer(c), as.double(f), as.integer(pivot), as.double(fpiv), 
        as.double(dom), as.double(coorddom), as.integer(indvois), 
        as.double(distvois), as.integer(orderf), as.integer(orderftmp), 
        as.integer(npop.est))
    print("End of Fortran function postprocesschain2")
    coord.grid <- read.table(filedom)[, 1:2]
    pmbr <- as.matrix(read.table(filedom)[, -(1:2)])
    pmp <- rep(NA, dim(pmbr)[1])
    for (ipix in 1:(dim(pmbr)[1])) {
        pmp[ipix] <- order(pmbr[ipix, ], decreasing = TRUE)[1]
    }
    write.table(cbind(coord.grid, pmp), file = paste(path.mcmc, 
        "modal.pop.txt", sep = ""), quote = FALSE, row.name = FALSE, 
        col.name = FALSE)
    indvois <- numeric(nindiv)
    distvois <- numeric(nindiv)
    u <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    pmp <- matrix(nr = nindiv, nc = npopmax, data = 0)
    out.res <- .Fortran(name = "pppmindiv2", PACKAGE = "Geneland", 
        as.integer(nindiv), as.double(t(coordinates)), as.integer(npopmax), 
        as.integer(nb.nuclei.max), as.integer(indvois), as.double(distvois), 
        as.double(u), as.integer(c), as.double(pmp), as.character(filenpop), 
        as.character(filenpp), as.character(fileu), as.character(filec), 
        as.character(fileperm), as.integer(nit), as.integer(thinning), 
        as.integer(burnin), as.integer(orderf), as.integer(npop.est), 
        as.integer(pivot))
    pmp <- matrix(nr = nindiv, nc = npopmax, data = out.res[[9]])
    mod.pop.indiv <- t(apply(pmp, 1, order))[, npopmax]
    write.table(cbind(coordinates, pmp), file = paste(path.mcmc, 
        "proba.pop.membership.indiv.txt", sep = ""), quote = FALSE, 
        row.name = FALSE, col.name = FALSE)
    write.table(cbind(coordinates, mod.pop.indiv), file = paste(path.mcmc, 
        "modal.pop.indiv.txt", sep = ""), quote = FALSE, row.name = FALSE, 
        col.name = FALSE)
    param <- c(paste("nxdom :", nxdom), paste("nydom :", nydom), 
        paste("burnin :", burnin))
    write.table(param, file = paste(path.mcmc, "postprocess.parameters.txt", 
        sep = ""), quote = FALSE, row.name = FALSE, col.name = FALSE)
}
