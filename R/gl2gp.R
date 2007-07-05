`gl2gp` <-
function (coordinates, genotypes, file) 
{
    nindiv <- nrow(genotypes)
    nloc <- ncol(genotypes)/2
    write.table("Header line", file, quote = F, col.name = F, 
        row.name = F)
    for (iloc in 1:nloc) {
        write.table(paste("loc", iloc), file, quote = F, col.name = F, 
            row.name = F, append = T)
    }
    for (iindiv in 1:nindiv) {
        write.table("Pop", file, quote = F, col.name = F, row.name = F, 
            append = T)
        string <- c(round(coordinates[, iindiv], digits = 3), 
            ",")
        for (iloc in 1:nloc) {
            all1 <- genotypes[iindiv, 2 * iloc - 1]
            all2 <- genotypes[iindiv, 2 * iloc]
            char1 <- "00"
            if (all1 < 10) {
                substring(char1, 2, 2) <- as.character(all1)
            }
            if ((all1 > 9) & (all1 < 100)) {
                substring(char1, 1, 2) <- as.character(all1)
            }
            char2 <- "00"
            if (all2 < 10) {
                substring(char2, 2, 2) <- as.character(all2)
            }
            if ((all2 > 9) & (all2 < 100)) {
                substring(char2, 1, 2) <- as.character(all2)
            }
            char <- paste(char1, char2, sep = "")
            string <- c(string, char)
        }
        if (iindiv < nindiv) {
            write.table(x = t(string), file = file, quote = F, 
                col.name = F, row.name = F, append = T)
        }
        if (iindiv == nindiv) {
            write.table(x = t(string), file = file, eol = "", 
                quote = F, col.name = F, row.name = F, append = T)
        }
    }
}
