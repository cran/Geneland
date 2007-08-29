`PosteriorMode` <-
function (coordinates, path.mcmc, plotit = TRUE, printit = FALSE, 
    file, main.title = "") 
{
    coordinates <- as.matrix(coordinates)
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "")
    param <- as.matrix(read.table(fileparam))
    delta.coord <- as.numeric(param[param[, 1] == "delta.coord", 
        3])
    npopmax <- as.numeric(param[param[, 1] == "npopmax", 3])
    param.postprocess <- as.matrix(read.table(paste(path.mcmc, 
        "postprocess.parameters.txt", sep = "")))
    nxdom <- as.numeric(param.postprocess[1, 3])
    nydom <- as.numeric(param.postprocess[2, 3])
    s <- coordinates
    filedom <- paste(path.mcmc, "proba.pop.membership.txt", sep = "/")
    dom.post <- as.matrix(read.table(filedom))[, -(1:2)]
    coord.grid <- as.matrix(read.table(filedom))[, (1:2)]
    s[, 1] <- s[, 1] - min(s[, 1])
    s[, 2] <- s[, 2] - min(s[, 2])
    xlim <- c(min(s[, 1]) - delta.coord/2, max(s[, 1]) + delta.coord/2)
    ylim <- c(min(s[, 2]) - delta.coord/2, max(s[, 2]) + delta.coord/2)
    Dx <- (xlim[2] - xlim[1])/(nxdom - 1)
    Dy <- (ylim[2] - ylim[1])/(nydom - 1)
    s.discr <- s
    s.discr[, 1] <- floor(s[, 1]/Dx) * Dx
    s.discr[, 2] <- floor(s[, 2]/Dy) * Dy
    is <- s.discr[, 1]/Dx + 1
    js <- s.discr[, 2]/Dy + 1
    ks <- (is - 1) * nydom + js
    map <- numeric(length(ks))
    for (k in 1:length(ks)) {
        map[k] <- order(dom.post[ks[k], ], decreasing = TRUE)[1]
    }
    map.dom <- t(apply(dom.post, 1, order))[, npopmax]
    if (plotit) {
        get(getOption("device"))()
        frame <- max(max(coordinates[, 1]) - min(coordinates[, 
            1]), max(coordinates[, 2]) - min(coordinates[, 2]))/40
        image(seq(min(coordinates[, 1] - delta.coord/2), max(coordinates[, 
            1] + delta.coord/2), length = nxdom), seq(min(coordinates[, 
            2] - delta.coord/2), max(coordinates[, 2] + delta.coord/2), 
            length = nydom), matrix(map.dom, nr = nxdom, nc = nydom, 
            byrow = TRUE), xlab = "", ylab = "", main = "", cex = 1.5, 
            cex.lab = 1.5, col = terrain.colors(npopmax), xlim = c(min(coordinates[, 
                1] - delta.coord/2 - frame), max(coordinates[, 
                1] + delta.coord/2 + frame)), ylim = c(min(coordinates[, 
                2] - delta.coord/2 - frame), max(coordinates[, 
                2] + delta.coord/2 + frame)), asp = 1)
        title(sub = "Posterior mode of population membership")
        title(main = main.title, pch = 16)
    }
    if (printit) {
        postscript(file)
        frame <- max(max(coordinates[, 1]) - min(coordinates[, 
            1]), max(coordinates[, 2]) - min(coordinates[, 2]))/40
        image(seq(min(coordinates[, 1] - delta.coord/2), max(coordinates[, 
            1] + delta.coord/2), length = nxdom), seq(min(coordinates[, 
            2] - delta.coord/2), max(coordinates[, 2] + delta.coord/2), 
            length = nydom), matrix(map.dom, nr = nxdom, nc = nydom, 
            byrow = TRUE), xlab = "", ylab = "", main = "", cex = 1.5, 
            cex.lab = 1.5, col = terrain.colors(npopmax), xlim = c(min(coordinates[, 
                1] - delta.coord/2 - frame), max(coordinates[, 
                1] + delta.coord/2 + frame)), ylim = c(min(coordinates[, 
                2] - delta.coord/2 - frame), max(coordinates[, 
                2] + delta.coord/2 + frame)), asp = 1)
        points(coordinates, pch = 16)
        title(main = main.title, sub = "Posterior mode of population membership")
        dev.off()
    }
}
