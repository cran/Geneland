`Plotnpop` <-
function (path.mcmc, burnin, printit = FALSE, file) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filenpop <- paste(path.mcmc, "populations.numbers.txt", sep = "")
    npop <- scan(filenpop)
    if (burnin > 0) {
        sub <- -(1:burnin)
    }
    else {
        sub <- 1:length(npop)
    }
    if (printit == TRUE) {
        postscript(file)
        par(mfrow = c(1, 2))
        plot((1:length(npop)) * thinning, npop, type = "l", ylab = "Number of classes", 
            xlab = paste("Index of MCMC iteration", "\n Whole chain", 
                sep = ""), ylim = c(1, max(npop) + 0.5))
        hist(npop[sub], plot = TRUE, prob = TRUE, breaks = seq(0.5, 
            max(npop) + 0.5, 1), xlab = paste("Nb. of pop. along the chain \n(after a burnin of ", 
            burnin, "x", thinning, "i t.)", sep = ""))
        dev.off()
    }
    else {
        par(mfrow = c(1, 2))
        plot((1:length(npop)) * thinning, npop, type = "l", ylab = "Number of classes", 
            xlab = paste("Index of MCMC iteration", "\n Whole chain", 
                sep = ""), ylim = c(1, max(npop) + 0.5))
        hist(npop[sub], plot = TRUE, prob = TRUE, breaks = seq(0.5, 
            max(npop) + 0.5, 1), xlab = paste("Nb. of pop. along the chain \n(after a burnin of ", 
            burnin, "x", thinning, " it.)", sep = ""))
    }
}
