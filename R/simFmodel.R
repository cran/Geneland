`simFmodel` <-
function (nindiv, coordinates, coord.lim, number.nuclei, coord.nuclei, 
    color.nuclei, nall, npop, freq.model = "Dirichlet", drift, 
    seed, plots = FALSE, ploth = FALSE) 
{
    if (!missing(seed)) {
        set.seed(seed)
    }
    if (missing(coordinates)) {
        coordinates <- rbind(runif(min = coord.lim[1], max = coord.lim[2], 
            nindiv), runif(min = coord.lim[3], max = coord.lim[4], 
            nindiv))
    }
    if (missing(coord.nuclei)) {
        coord.nuclei <- rbind(runif(min = 0, max = 1, number.nuclei), 
            runif(min = 0, max = 1, number.nuclei))
        color.nuclei <- sample(x = 1:npop, size = number.nuclei, 
            replace = TRUE)
    }
    ppvois <- numeric(nindiv)
    for (iindiv in 1:nindiv) {
        k <- 1
        dd <- 1e+09
        for (ipp in 1:number.nuclei) {
            ddnew <- (coord.nuclei[1, ipp] - coordinates[1, iindiv])^2 + 
                (coord.nuclei[2, ipp] - coordinates[2, iindiv])^2
            if (ddnew < dd) {
                k <- ipp
                dd <- ddnew
            }
        }
        ppvois[iindiv] <- k
    }
    nloc <- length(nall)
    if (freq.model == "Falush") {
        if (missing(drift)) {
            drift <- rbeta(shape1 = 2, shape2 = 20, npop)
        }
        fa <- matrix(nr = nloc, nc = max(nall), data = -999)
        for (iloc in 1:nloc) {
            fa[iloc, 1:nall[iloc]] <- rexp(n = nall[iloc])
            fa[iloc, 1:nall[iloc]] <- fa[iloc, 1:nall[iloc]]/sum(fa[iloc, 
                1:nall[iloc]])
        }
    }
    freq <- array(dim = c(npop, nloc, max(nall)), data = -999)
    for (iclass in 1:npop) {
        for (iloc in 1:nloc) {
            if (freq.model == "Dirichlet") {
                freq[iclass, iloc, 1:nall[iloc]] <- rgamma(n = nall[iloc], 
                  scale = 1, shape = 1)
                freq[iclass, iloc, 1:nall[iloc]] <- freq[iclass, 
                  iloc, 1:nall[iloc]]/sum(freq[iclass, iloc, 
                  1:nall[iloc]])
            }
            if (freq.model == "Falush") {
                freq[iclass, iloc, 1:nall[iloc]] <- rgamma(n = nall[iloc], 
                  scale = 1, shape = fa[iloc, (1:nall[iloc])] * 
                    (1 - drift[iclass])/drift[iclass])
                freq[iclass, iloc, 1:nall[iloc]] <- freq[iclass, 
                  iloc, 1:nall[iloc]]/sum(freq[iclass, iloc, 
                  1:nall[iloc]])
            }
        }
    }
    if (npop > 1) {
        if (freq[1, 1, 1] < freq[2, 1, 1]) {
            tmp <- freq[2, , ]
            freq[2, , ] <- freq[1, , ]
            freq[1, , ] <- tmp
        }
    }
    z <- matrix(nr = nindiv, nc = nloc * 2)
    for (iclass in 1:npop) {
        subclass <- (1:nindiv)[color.nuclei[ppvois] == iclass]
        for (iloc in 1:nloc) {
            z[subclass, c(2 * iloc - 1, 2 * iloc)] <- sample(x = 1:nall[iloc], 
                size = 2 * length(subclass), prob = freq[iclass, 
                  iloc, 1:nall[iloc]], replace = TRUE)
        }
    }
    if (plots == TRUE) {
        get(getOption("device"))()
        plot(coordinates[1, ], coordinates[2, ], xlab = "x coordinates", 
            ylab = "y coordinates")
        points(coord.nuclei[1, ], coord.nuclei[2, ], col = color.nuclei, 
            cex = 2, pch = 16)
        text(coord.nuclei[1, ], coord.nuclei[2, ], 1:number.nuclei, 
            pos = 1, col = 2, pch = 2, cex = 1.2)
        text(coordinates[1, ], coordinates[2, ], ppvois, pos = 1)
    }
    if (ploth == TRUE) {
        if (freq.model == "Falush") {
            get(getOption("device"))()
            par(mfrow = c(floor(sqrt(nloc) + 1), floor(sqrt(nloc))))
            for (iloc in 1:nloc) {
                plot(1:nall[iloc], fa[iloc, 1:nall[iloc]], type = "h", 
                  col = 2, xlab = "", axes = FALSE, sub = paste("Locus", 
                    iloc), ylim = c(0, 1), main = "Frequencies in ancestral population")
            }
        }
        get(getOption("device"))()
        for (iclass in 1:npop) {
            get(getOption("device"))()
            par(mfrow = c(floor(sqrt(nloc) + 1), floor(sqrt(nloc))))
            for (iloc in 1:nloc) {
                hist(c(z[color.nuclei[ppvois] == iclass, 2 * 
                  (iloc - 1) + 1], z[color.nuclei[ppvois] == 
                  iclass, 2 * (iloc - 1) + 2]), breaks = seq(0.5, 
                  nall[iloc] + 0.5, 1), prob = TRUE, main = paste("Histogram of freq. in pop.", 
                  iclass, ", locus", iloc), xlab = "", ylim = c(0, 
                  1), axes = FALSE)
                points(1:nall[iloc], freq[iclass, iloc, 1:nall[iloc]], 
                  type = "h", col = 2)
            }
        }
    }
    if (freq.model == "Dirichlet") {
        res <- list(coordinates = t(coordinates), genotypes = z, 
            allele.numbers = nall, number.nuclei = number.nuclei, 
            coord.nuclei = t(coord.nuclei), color.nuclei = color.nuclei, 
            frequencies = freq, index.nearest.nucleus = ppvois)
        return(res)
    }
    if (freq.model == "Falush") {
        res <- list(coordinates = t(coordinates), genotypes = z, 
            allele.numbers = nall, number.nuclei = number.nuclei, 
            coord.nuclei = t(coord.nuclei), color.nuclei = color.nuclei, 
            frequencies = freq, ancestral.frequencies = fa, drifts = drift, 
            index.nearest.nucleus = ppvois)
        return(res)
    }
}
