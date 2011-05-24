PostProcessMultChain <-
function (coordinates = NULL, genotypes, ploidy, path.all, nrun, 
    nxdom, nydom, burnin) 
{
    print("Reading parameters")
    path.mcmc <- paste(path.all, "1/", sep = "/")
    param <- as.matrix(read.table(paste(path.mcmc, "parameters.txt", 
        sep = "")))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    nb.nuclei.max <- as.numeric(param[param[, 1] == "nb.nuclei.max", 
        3])
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    spatial <- as.logical(param[param[, 1] == "spatial", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    if (ploidy == 1) {
        data.tmp <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) * 
            2)
        data.tmp[, seq(1, ncol(genotypes) * 2 - 1, 2)] <- genotypes
        data.tmp[, seq(2, ncol(genotypes) * 2, 2)] <- genotypes
        genotypes <- data.tmp
    }
    print("defining dummy coordinates if coord are missing")
    nindiv <- nrow(genotypes)
    if (is.null(coordinates)) {
        if (spatial) {
            stop("Please give spatial coordinates of individuals or set argument spatial to FALSE")
        }
        else {
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
        if (nrow(coordinates) != nindiv) {
            print(paste("number of individuals in data matrix =", 
                nindiv))
            print(paste("number of individuals in coordinate matrix =", 
                nrow(coordinates)))
            stop("Number of rows in coordinate matrix and data matrices do not match ")
        }
    }
    coordinates <- as.matrix(coordinates)
    data.tmp <- FormatGenotypes(as.matrix(genotypes), ploidy = ploidy)
    genotypes <- data.tmp$genotypes
    allele.numbers <- data.tmp$allele.numbers
    if (filter.null.alleles) {
        allele.numbers <- allele.numbers + 1
    }
    nloc <- length(allele.numbers)
    npop <- rep(NA, nit/thinning)
    for (irun in 1:nrun) {
        path.mcmc <- paste(path.all, irun, "/", sep = "/")
        npop[(1:(nit/thinning)) + (irun - 1) * nit/thinning] <- scan(paste(path.mcmc, 
            "populations.numbers.txt", sep = ""), quiet = TRUE)
    }
    npop.est <- order(hist(npop[-(1:burnin)], breaks = seq(0.5, 
        npopmax + 0.5, 1), plot = FALSE)$counts, decreasing = TRUE)[1]
    lpd <- rep(NA, nit/thinning)
    for (irun in 1:nrun) {
        path.mcmc <- paste(path.all, irun, "/", sep = "/")
        lpd[(1:(nit/thinning)) + (irun - 1) * nit/thinning] <- scan(paste(path.mcmc, 
            "log.posterior.density.txt", sep = ""), quiet = TRUE)
    }
    pivot <- order(lpd, decreasing = TRUE)[1]
    irun.piv <- 1 + floor((pivot - 1)/(nit/thinning))
    char.irun.piv <- paste(irun.piv)
    tmp <- character(255)
    substr(tmp, 1, nchar(path.all)) <- path.all
    nchar.path.all <- nchar(path.all)
    nchar.irun.piv <- nchar(char.irun.piv)
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
    out.res <- .Fortran(name = "postprocessmultchain", PACKAGE = "Geneland", 
        as.integer(nxdom), as.integer(nydom), as.integer(burnin), 
        as.integer(ninrub), as.integer(npopmax), as.integer(nb.nuclei.max), 
        as.integer(nindiv), as.integer(nloc), as.integer(allele.numbers), 
        as.integer(nalmax), as.double(xlim), as.double(ylim), 
        as.double(delta.coord), as.integer(nit), as.integer(thinning), 
        as.integer(nrun), as.character(path.all), as.integer(nchar.path.all), 
        as.double(t(coordinates)), as.double(u), as.integer(c), 
        as.double(f), as.integer(pivot), as.character(char.irun.piv), 
        as.integer(nchar.irun.piv), as.double(fpiv), as.double(dom), 
        as.double(coorddom), as.integer(indvois), as.double(distvois), 
        as.integer(orderf), as.integer(orderftmp), as.integer(npop.est))
    filedom <- paste(path.all, "proba.pop.membership.txt", sep = "")
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
    nchar.path.all <- nchar(path.all)
    out.res <- .Fortran(name = "pppmindivmultchain", PACKAGE = "Geneland", 
        as.integer(nindiv), as.double(t(coordinates)), as.integer(npopmax), 
        as.integer(nb.nuclei.max), as.integer(indvois), as.double(distvois), 
        as.double(u), as.integer(c), as.double(pmp), as.character(path.all), 
        as.integer(nchar.path.all), as.integer(nrun), as.integer(nit), 
        as.integer(thinning), as.integer(burnin), as.integer(orderf), 
        as.integer(npop.est), as.integer(pivot))
    pmp <- matrix(nr = nindiv, nc = npopmax, data = out.res[[9]])
    mod.pop.indiv <- numeric(nindiv)
    for (iindiv in 1:nindiv) {
        mod.pop.indiv[iindiv] <- order(pmp[iindiv, ], decreasing = TRUE)[1]
    }
    write.table(cbind(coordinates, pmp), file = paste(path.all, 
        "proba.pop.membership.indiv.txt", sep = ""), quote = FALSE, 
        row.name = FALSE, col.name = FALSE)
    write.table(cbind(coordinates, mod.pop.indiv), file = paste(path.all, 
        "modal.pop.indiv.txt", sep = ""), quote = FALSE, row.name = FALSE, 
        col.name = FALSE)
    param <- rbind(param, paste("nrun :", nrun))
    write.table(param, file = paste(path.all, "parameters.txt", 
        sep = ""), quote = FALSE, row.name = FALSE, col.name = FALSE)
    param <- c(paste("nxdom :", nxdom), paste("nydom :", nydom))
    write.table(param, file = paste(path.all, "postprocess.parameters.txt", 
        sep = ""), quote = FALSE, row.name = FALSE, col.name = FALSE)
}
