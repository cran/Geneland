`PostProcessChain` <-
function (coordinates, genotypes, path.mcmc, nxdom, nydom, burnin) 
{
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
    if (ploidy == 1) {
        data.tmp <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) * 
            2)
        data.tmp[, seq(1, ncol(genotypes) * 2 - 1, 2)] <- genotypes
        data.tmp[, seq(2, ncol(genotypes) * 2, 2)] <- genotypes
        genotypes <- data.tmp
    }
    coordinates <- as.matrix(coordinates)
    data.tmp <- FormatGenotypes(as.matrix(genotypes))
    genotypes <- data.tmp$genotypes
    allele.numbers <- data.tmp$allele.numbers
    filenpop <- paste(path.mcmc, "populations.numbers.txt", sep = "")
    filenpp <- paste(path.mcmc, "nuclei.numbers.txt", sep = "")
    fileu <- paste(path.mcmc, "coord.nuclei.txt", sep = "")
    filec <- paste(path.mcmc, "color.nuclei.txt", sep = "")
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    filefperm <- paste(path.mcmc, "fperm.txt", sep = "")
    filedom <- paste(path.mcmc, "proba.pop.membership.txt", sep = "")
    filedomperm <- paste(path.mcmc, "proba.pop.membership.perm.txt", 
        sep = "")
    nindiv <- nrow(genotypes)
    nloc <- length(allele.numbers)
    dom <- matrix(nr = nxdom * nydom, nc = npopmax, data = 0)
    domperm <- matrix(nr = nxdom * nydom, nc = npopmax, data = 0)
    coorddom <- matrix(nr = 2, nc = nxdom * nydom, data = -999)
    indvois <- numeric(nxdom * nydom)
    distvois <- numeric(nxdom * nydom)
    f11 <- numeric(npopmax)
    orderf11 <- numeric(npopmax)
    nallmax <- max(allele.numbers)
    u <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    out.res <- .Fortran(name = "postprocesschain", PACKAGE = "Geneland", 
        as.integer(nxdom), as.integer(nydom), as.integer(burnin), 
        as.integer(npopmax), as.integer(nb.nuclei.max), as.integer(nloc), 
        as.integer(nindiv), as.integer(nloc), as.integer(nallmax), 
        as.single(delta.coord), as.integer(nit/thinning), as.character(filenpp), 
        as.character(filenpop), as.character(fileu), as.character(filec), 
        as.character(filef), as.character(filefperm), as.character(filedom), 
        as.character(filedomperm), as.single(t(coordinates)), 
        as.single(u), as.integer(c), as.single(dom), as.single(domperm), 
        as.single(coorddom), as.integer(indvois), as.single(distvois), 
        as.single(f11), as.single(orderf11))
    param <- c(paste("nxdom :", nxdom), paste("nydom :", nydom))
    write.table(param, file = paste(path.mcmc, "postprocess.parameters.txt", 
        sep = ""), quote = FALSE, row.name = FALSE, col.name = FALSE)
    indvois <- numeric(nindiv)
    distvois <- numeric(nindiv)
    u <- matrix(nr = 2, nc = nb.nuclei.max, data = -999)
    c <- rep(times = nb.nuclei.max, -999)
    pmp <- matrix(nr = nindiv, nc = npopmax, data = 0)
    out.res <- .Fortran(name = "pppmindiv", PACKAGE = "Geneland", 
        as.integer(nindiv), as.single(t(coordinates)), as.integer(npopmax), 
        as.integer(nb.nuclei.max), as.integer(indvois), as.single(distvois), 
        as.single(u), as.integer(c), as.single(pmp), as.character(filenpp), 
        as.character(fileu), as.character(filec), as.integer(nit/thinning), 
        as.integer(burnin))
    pmp <- matrix(nr = nindiv, nc = npopmax, data = out.res[[9]])
    mod.pop.indiv <- numeric(nindiv)
    for (iindiv in 1:nindiv) {
        mod.pop.indiv[iindiv] <- order(pmp[iindiv, ], decreasing = TRUE)[1]
    }
    write.table(cbind(coordinates, pmp), file = paste(path.mcmc, 
        "proba.pop.membership.indiv.txt", sep = ""), quote = FALSE, 
        row.name = FALSE, col.name = FALSE)
    write.table(cbind(coordinates, mod.pop.indiv), file = paste(path.mcmc, 
        "modal.pop.indiv.txt", sep = ""), quote = FALSE, row.name = FALSE, 
        col.name = FALSE)
}
