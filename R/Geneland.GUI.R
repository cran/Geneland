`Geneland.GUI` <-
function (lib.loc = NULL) 
{
    require(tcltk)
    tt <- tktoplevel()
    tkwm.title(tt, "Geneland - Graphical Interface")
    tkwm.geometry(tt, "+100+100")
    image1 <- tclVar()
    tcl("image", "create", "photo", image1, file = system.file("images/geneland6.gif", 
        package = "Geneland", lib.loc = lib.loc))
    imgAsLabel <- tklabel(tt, image = image1)
    tkpack(imgAsLabel)
    tkfocus(tt)
    tcl("after", "3000", "destroy", tt)
    require(Geneland)
    idb.dataset <- 0
    globalcoordinates <- NULL
    globalgenotypes <- 0
    globallabels <- NA
    burnin <- tclVar(0)
    falush <- 0
    LogState <- tclVar(0)
    labelcoordtext <- tclVar("")
    labelgenotext <- tclVar("")
    imageconfigure <- tclVar()
    imagerun <- tclVar()
    imagerun2 <- tclVar()
    imagefmodel <- tclVar()
    imagefstat <- tclVar()
    imageplot <- tclVar()
    imagepostprocess <- tclVar()
    imagedraw <- tclVar()
    imageok <- tclVar()
    imagepleasewait <- tclVar()
    imageibd <- tclVar()
    imageconvert <- tclVar()
    imageline <- tclVar()
    tcl("image", "create", "photo", imageconfigure, file = system.file("images/icon-configure.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagerun, file = system.file("images/icon-run.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagerun2, file = system.file("images/icon-run2.gif", 
        package = "Geneland"))
    tcl("image", "create", "photo", imagefmodel, file = system.file("images/icon-fmodel.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagefstat, file = system.file("images/icon-fstat.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imageplot, file = system.file("images/icon-output.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagepostprocess, file = system.file("images/icon-postprocess.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagedraw, file = system.file("images/icon-draw.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imageok, file = system.file("images/icon-ok.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imageibd, file = system.file("images/icon-ibd.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imagepleasewait, file = system.file("images/pleasewait.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imageconvert, file = system.file("images/icon-convert.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tcl("image", "create", "photo", imageline, file = system.file("images/line.gif", 
        package = "Geneland", lib.loc = lib.loc))
    tkwait.window(tt)
    coordinatesfile <- tclVar("")
    genotypefile <- tclVar("")
    outputdir <- tclVar("")
    labelsfile <- tclVar("")
    advanced <- tclVar(0)
    sep1 <- tclVar("White space")
    sep2 <- tclVar("White space")
    md <- tclVar("NA")
    processors <- tclVar(1)
    cluster <- 0
    usecluster <- FALSE
    tt <- tktoplevel()
    tkwm.title(tt, "Geneland - Graphical Interface")
    tkwm.geometry(tt, "750x600+100+100")
    tkwm.resizable(tt, 0, 0)
    ttconf <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttrun <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttinit <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttpost <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttsimf <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttplot <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttfstat <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttibd <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttplot2 <- tkframe(tt, borderwidth = 2, relief = "sunken")
    ttpan <- tkframe(tt)
    numberofdigits <- function(number) {
        i <- 0
        while (number >= 1) {
            i <- i + 1
            number <- number/10
        }
        return(i)
    }
    matrix2str <- function(mat) {
        if (is.null(mat)) 
            str = "NULL"
        else {
            str = "rbind("
            for (i in 1:NROW(mat)) {
                str = paste(str, "c(", mat[i, 1], sep = "")
                for (j in 2:NCOL(mat)) {
                  str = paste(str, ",", as.character(mat[i, j]), 
                    sep = "")
                }
                str = paste(str, "),", sep = "")
            }
            str = strtrim(str, nchar(str) - 1)
            str = paste(str, ")", sep = "")
        }
        return(str)
    }
    fadvanced <- function() {
        if (tclvalue(advanced) == 1) {
            tkconfigure(labelsimulation, state = "normal")
            tkconfigure(buttonsimfmodel, state = "normal")
            tkconfigure(buttonibd, state = "normal")
            tkconfigure(buttonplot2, state = "normal")
            tcl(toolsMenu, "entryconfigure", "1", state = "normal")
        }
        else {
            tkconfigure(labelsimulation, state = "disabled")
            tkconfigure(buttonsimfmodel, state = "disabled")
            tkconfigure(buttonibd, state = "disabled")
            tkconfigure(buttonplot2, state = "disabled")
            tcl(toolsMenu, "entryconfigure", "1", state = "disabled")
        }
        run()
    }
    helpWindow <- function() {
        help.start("Geneland")
        tkfocus(tt)
    }
    creditsWindow <- function() {
        ttcredits <- tktoplevel(parent = .TkRoot)
        tkwm.title(ttcredits, "Credits")
        label1.widget <- tklabel(ttcredits, text = "Authors:")
        label2.widget <- tklabel(ttcredits, text = " ")
        label3.widget <- tklabel(ttcredits, text = "Gilles Guillot:")
        label4.widget <- tklabel(ttcredits, text = "Statistical Fortran and R code")
        label5.widget <- tklabel(ttcredits, text = "")
        label6.widget <- tklabel(ttcredits, text = "Filipe Santos:")
        label7.widget <- tklabel(ttcredits, text = " Graphical interface (code and design)")
        label8.widget <- tklabel(ttcredits, text = "")
        label9.widget <- tklabel(ttcredits, text = "Arnaud Estoup:")
        label10.widget <- tklabel(ttcredits, text = " Graphical interface (design and test)")
        tkgrid(label1.widget, row = 1, column = 1, sticky = "w")
        tkgrid(label2.widget, row = 2, column = 1, sticky = "w")
        tkgrid(label3.widget, row = 3, column = 1, sticky = "w")
        tkgrid(label4.widget, row = 3, column = 2, sticky = "w")
        tkgrid(label5.widget, row = 4, column = 1, sticky = "w")
        tkgrid(label6.widget, row = 5, column = 1, sticky = "w")
        tkgrid(label7.widget, row = 5, column = 2, sticky = "w")
        tkgrid(label8.widget, row = 6, column = 1, sticky = "w")
        tkgrid(label9.widget, row = 7, column = 1, sticky = "w")
        tkgrid(label10.widget, row = 7, column = 2, sticky = "w")
    }
    configure <- function() {
        label1.widget <- tklabel(ttconf, text = tclvalue(coordinatesfile), 
            width = 40, justify = "left")
        tkconfigure(label1.widget, textvariable = coordinatesfile)
        getcoordinatesfile <- function() {
            tclvalue(coordinatesfile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *}", 
                title = "Choose coordinate file"))
            if (tclvalue(coordinatesfile) != "") {
                if (tclvalue(sep1) == "White space") 
                  tclvalue(sep1) <- ""
                globalcoordinates <<- as.matrix(read.table(tclvalue(coordinatesfile), 
                  sep = tclvalue(sep1)))
                Log(paste("as.matrix(read.table(", tclvalue(coordinatesfile), 
                  "),sep=\"", tclvalue(sep1), "\")", sep = ""), 
                  "[SUCCESS] ")
                if (tclvalue(sep1) == "") 
                  tclvalue(sep1) <- "White space"
                tclvalue(labelcoordtext) = "Coordinate: File loaded"
            }
            else {
                globalcoordinates <<- NULL
                tclvalue(labelcoordtext) = "Coordinate: Data unloaded"
            }
            tkfocus(tt)
        }
        button1.widget <- tkbutton(ttconf, text = "Coordinate file", 
            command = getcoordinatesfile, width = 15, justify = "left")
        tkgrid(button1.widget, row = 2, column = 1, sticky = "we")
        tkgrid(label1.widget, row = 2, column = 2, sticky = "we")
        label2.widget <- tklabel(ttconf, text = tclvalue(genotypefile), 
            width = 40, justify = "left")
        tkconfigure(label2.widget, textvariable = genotypefile)
        getgenotypefile <- function() {
            tclvalue(genotypefile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *}", 
                title = "Choose genotype file"))
            if (tclvalue(genotypefile) != "") {
                if (tclvalue(sep2) == "White space") 
                  tclvalue(sep2) <- ""
                globalgenotypes <<- as.matrix(read.table(tclvalue(genotypefile), 
                  sep = tclvalue(sep2), na.strings = md))
                Log(paste("as.matrix(read.table(", tclvalue(genotypefile), 
                  "),sep=\"", tclvalue(sep2), "\",na.strings=", 
                  md, "\")", sep = ""), "[SUCCESS] ")
                if (tclvalue(sep2) == "") 
                  tclvalue(sep2) <- "White space"
                tclvalue(labelgenotext) = "Genotype:  File loaded"
            }
            else {
                globalgenotypes <<- 0
                tclvalue(labelgenotext) = "Genotype:  Data unloaded"
            }
            tkfocus(tt)
        }
        button2.widget <- tkbutton(ttconf, text = "Genotype file", 
            command = getgenotypefile, width = 15, justify = "left")
        tkgrid(button2.widget, row = 3, column = 1, sticky = "we")
        tkgrid(label2.widget, row = 3, column = 2, sticky = "we")
        label3.widget <- tklabel(ttconf, text = tclvalue(outputdir), 
            width = 40, justify = "left")
        tkconfigure(label3.widget, textvariable = outputdir)
        getoutputdir <- function() {
            tclvalue(outputdir) <- tclvalue(tkchooseDirectory(parent = tt, 
                title = "Please choose an output directory"))
            if (tclvalue(outputdir) != "") {
                tcl("regsub", "-all", "\\\\", tclvalue(outputdir), 
                  "/", outputdir)
                tcl("append", outputdir, "/")
                auxblink <<- 2
                tkconfigure(extralabel.widget, text = paste("Output Directory : ", 
                  tclvalue(outputdir), sep = ""))
            }
            tkfocus(tt)
        }
        button3.widget <- tkbutton(ttconf, text = "Output directory", 
            command = getoutputdir, width = 15, justify = "left")
        tkgrid(button3.widget, row = 4, column = 1, sticky = "we")
        tkgrid(label3.widget, row = 4, column = 2, sticky = "we")
        label5.widget <- tklabel(ttconf, font = "*-Courier--i-normal--12-*", 
            foreground = "blue", text = "-------------Optional-------------", 
            width = 40, justify = "left")
        tkgrid(label5.widget, row = 5, column = 1, columnspan = 2, 
            sticky = "we")
        label4.widget <- tklabel(ttconf, text = tclvalue(labelsfile), 
            width = 40, justify = "left")
        tkconfigure(label4.widget, textvariable = labelsfile)
        getlabelsfile <- function() {
            tclvalue(labelsfile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *}", 
                title = "Choose individual label file"))
            globallabels <<- as.matrix(read.table(tclvalue(labelsfile)))
            tkfocus(tt)
        }
        button4.widget <- tkbutton(ttconf, text = "Individual label file", 
            command = getlabelsfile, width = 15, justify = "left")
        tkgrid(button4.widget, row = 6, column = 1, sticky = "we")
        tkgrid(label4.widget, row = 6, column = 2, sticky = "we")
    }
    run <- function() {
        tclvalue(burnin) <<- 0
        testnumberpop <- tclVar(0)
        RunmcmcFmodel <- function() {
            if (tclvalue(ploidy) == "Haploid") 
                dploidy <- 1
            else dploidy <- 2
            if (tclvalue(advanced) != 1) 
                tclvalue(npopinit) <- tclvalue(npopmax)
            varnpop <- TRUE
            if (as.integer(tclvalue(npopmin)) == as.integer(tclvalue(npopmax))) 
                varnpop <- FALSE
            tttry <- tktoplevel(tt)
            tkgrab(tttry)
            tkwm.geometry(tttry, "+200+200")
            tkwm.title(tttry, "wait")
            warn <- tklabel(tttry, image = imagepleasewait)
            tkpack(warn)
            tkfocus(tttry)
            tcl("update")
            if (tclvalue(testnumberpop) == 0) {
                print("Starting...")
                Sys.sleep(0.5)
                tcl("update")
                err <- try(MCMC(coordinates = globalcoordinates, 
                  genotypes = globalgenotypes, ploidy = dploidy, 
                  path.mcmc = tclvalue(outputdir), rate.max = as.numeric(tclvalue(rate)), 
                  delta.coord = as.numeric(tclvalue(delta)), 
                  npopmin = as.numeric(tclvalue(npopmin)), npopinit = as.numeric(tclvalue(npopinit)), 
                  npopmax = as.numeric(tclvalue(npopmax)), nb.nuclei.max = as.numeric(tclvalue(nuclei)), 
                  nit = as.numeric(tclvalue(nit)), thinning = as.numeric(tclvalue(thinning)), 
                  freq.model = tclvalue(freq), shape1 = as.numeric(tclvalue(gshape1)), 
                  shape2 = as.numeric(tclvalue(gshape2)), varnpop = varnpop, 
                  spatial = as.logical(tclvalue(spatial)), jcf = as.logical(tclvalue(jcf)), 
                  filter.null.alleles = as.logical(tclvalue(null))), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  Log(paste("MCMC(coordinates=", matrix2str(globalcoordinates), 
                    ",genotypes=", matrix2str(globalgenotypes), 
                    ",ploidy=", dploidy, ",path.mcmc=\"", tclvalue(outputdir), 
                    "\",rate.max=", as.numeric(tclvalue(rate)), 
                    ",delta.coord=", as.numeric(tclvalue(delta)), 
                    ",npopmin=", as.numeric(tclvalue(npopmin)), 
                    ",npopinit=", as.numeric(tclvalue(npopinit)), 
                    ",npopmax=", as.numeric(tclvalue(npopmax)), 
                    ",nb.nuclei.max=", as.numeric(tclvalue(nuclei)), 
                    ",nit=", as.numeric(tclvalue(nit)), ",thinning=", 
                    as.numeric(tclvalue(thinning)), ",freq.model=\"", 
                    tclvalue(freq), "\",shape1=", shape1 = as.numeric(tclvalue(gshape1)), 
                    ",shape2=", shape2 = as.numeric(tclvalue(gshape2)), 
                    ",varnpop=", varnpop, ",spatial=", as.logical(tclvalue(spatial)), 
                    ",jcf=", as.logical(tclvalue(jcf)), ",filter.null.alleles =", 
                    as.logical(tclvalue(null)), ")", sep = ""), 
                    "[FAILED] ")
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  Log(paste("MCMC(coordinates=", matrix2str(globalcoordinates), 
                    ",genotypes=", matrix2str(globalgenotypes), 
                    ",ploidy=", dploidy, ",path.mcmc=\"", tclvalue(outputdir), 
                    "\",rate.max=", as.numeric(tclvalue(rate)), 
                    ",delta.coord=", as.numeric(tclvalue(delta)), 
                    ",npopmin=", as.numeric(tclvalue(npopmin)), 
                    ",npopinit=", as.numeric(tclvalue(npopinit)), 
                    ",npopmax=", as.numeric(tclvalue(npopmax)), 
                    ",nb.nuclei.max=", as.numeric(tclvalue(nuclei)), 
                    ",nit=", as.numeric(tclvalue(nit)), ",thinning=", 
                    as.numeric(tclvalue(thinning)), ",freq.model=\"", 
                    tclvalue(freq), "\",shape1=", shape1 = as.numeric(tclvalue(gshape1)), 
                    ",shape2=", shape2 = as.numeric(tclvalue(gshape2)), 
                    ",varnpop=", varnpop, ",spatial=", as.logical(tclvalue(spatial)), 
                    ",jcf=", as.logical(tclvalue(jcf)), ",filter.null.alleles =", 
                    as.logical(tclvalue(null)), ")", sep = ""), 
                    "[SUCCESS] ")
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                  if (tclvalue(freq) == "Correlated") 
                    falush <<- 1
                  else falush <<- 0
                }
            }
            else {
                probs <<- c()
                pops <<- c()
                runs <<- c()
                txtleft <- c()
                txtmidle <- c()
                txtright <- c()
                reburnvalue <- tclVar(0)
                rbValue <- tclVar(0)
                auxoutputdir <- outputdir
                totaltime <- 0
                runtime <- 0
                defineOutdir <- function(n) {
                  outputdir <<- tclVar(paste(tclvalue(auxoutputdir), 
                    as.character(n), "/", sep = ""))
                  tkconfigure(extralabel.widget, text = paste("Output Directory : ", 
                    tclvalue(outputdir), sep = ""))
                }
                makebutton <- function(n) {
                  rbload <- tclVar()
                  rbload <- tkradiobutton(tttextpop, variable = rbValue, 
                    value = n, borderwidth = 2, command = function() defineOutdir(n))
                  tkwindow.create(load, "end", window = rbload)
                }
                Sort <- function() {
                  sorted <- order(probs, decreasing = TRUE)
                  new <- rbind(probs[sorted], pops[sorted], runs[sorted])
                  tkdelete(left, "1.0", "end")
                  tkdelete(midle, "1.0", "end")
                  tkdelete(right, "1.0", "end")
                  tkdelete(load, "1.0", "end")
                  txtleft <<- c()
                  txtmidle <<- c()
                  txtright <<- c()
                  for (i in 1:length(runs)) {
                    ltext = tklabel(tttextpop, text = as.character(new[3, 
                      i]), border = 2)
                    tkwindow.create(left, "end", window = ltext)
                    tkinsert(left, "end", "\n")
                    lmiddle = tklabel(tttextpop, text = new[2, 
                      i], border = 2)
                    tkwindow.create(midle, "end", window = lmiddle)
                    tkinsert(midle, "end", "\n")
                    lright = tklabel(tttextpop, text = new[1, 
                      i], border = 2)
                    tkwindow.create(right, "end", window = lright)
                    tkinsert(right, "end", "\n")
                    makebutton(new[3, i])
                    tkinsert(load, "end", "\n")
                    txtleft <<- c(txtleft, as.character(new[3, 
                      i]))
                    txtmidle <<- c(txtmidle, new[2, i])
                    txtright <<- c(txtright, new[1, i])
                  }
                }
                Reburn <- function() {
                  if (as.numeric(tclvalue(reburnvalue)) >= as.numeric(tclvalue(nit))/as.numeric(tclvalue(thinning))) {
                    tclvalue(reburnvalue) <- as.numeric(tclvalue(nit))/as.numeric(tclvalue(thinning)) - 
                      1
                    tkconfigure(burnentry, textvariable = reburnvalue)
                  }
                  tkdelete(left, "1.0", "end")
                  tkdelete(midle, "1.0", "end")
                  tkdelete(right, "1.0", "end")
                  tkdelete(load, "1.0", "end")
                  txtleft <<- c()
                  txtmidle <<- c()
                  txtright <<- c()
                  probs <<- c()
                  pops <<- c()
                  runs <<- c()
                  for (i in 1:as.numeric(tclvalue(ntestpop))) {
                    tempoutputdir <- paste(tclvalue(auxoutputdir), 
                      as.character(i), "/", sep = "")
                    ltext = tklabel(tttextpop, text = as.character(i), 
                      border = 3)
                    tkwindow.create(left, "end", window = ltext)
                    tkinsert(left, "end", "\n")
                    runs <<- c(runs, i)
                    file <- try(scan(paste(tempoutputdir, "populations.numbers.txt", 
                      sep = "")), silent = TRUE)
                    if (tclvalue(reburnvalue) != 0) 
                      dist <- hist(file[-(1:as.numeric(tclvalue(reburnvalue)))], 
                        plot = FALSE, breaks = seq(0.5, max(file) + 
                          0.5, 1))
                    else dist <- hist(file, plot = FALSE, breaks = seq(0.5, 
                      max(file) + 0.5, 1))
                    firsttime <- 0
                    straux <- ""
                    for (j in 1:length(dist$counts)) {
                      if (dist$counts[j] == max(dist$counts)) {
                        if (firsttime == 0) {
                          straux <- as.character(dist$mids[j])
                          firsttime <- 1
                        }
                        else {
                          straux <- paste(straux, " and ", sep = "")
                          straux <- paste(straux, as.character(dist$mids[j]), 
                            sep = "")
                        }
                      }
                    }
                    straux <- paste(straux, " ( ", sep = "")
                    straux <- paste(straux, as.character(as.double(max(dist$density) * 
                      100)), sep = "")
                    straux <- paste(straux, " %) ", sep = "")
                    lmiddle = tklabel(tttextpop, text = straux, 
                      border = 3)
                    tkwindow.create(midle, "end", window = lmiddle)
                    pops <<- c(pops, straux)
                    tkinsert(midle, "end", "\n")
                    file <- try(scan(paste(tempoutputdir, "log.posterior.density.txt", 
                      sep = "")), silent = TRUE)
                    if (tclvalue(reburnvalue) != 0) 
                      mpd <- mean(file[-(1:as.numeric(tclvalue(reburnvalue)))])
                    else mpd <- mean(file)
                    lright = tklabel(tttextpop, text = mpd, border = 3)
                    tkwindow.create(right, "end", window = lright)
                    tkinsert(right, "end", "\n")
                    probs <<- c(probs, mpd)
                    makebutton(i)
                    tkinsert(load, "end", "\n")
                    txtleft <<- c(txtleft, as.character(i))
                    txtmidle <<- c(txtmidle, straux)
                    txtright <<- c(txtright, mpd)
                  }
                }
                Output <- function() {
                  outfile <- tclvalue(tkgetSaveFile(filetypes = "{{.txt} *.txt}", 
                    title = "Save to file"))
                  if (outfile != "") {
                    zz1 <- file(outfile, "w")
                    for (i in 1:length(runs)) {
                      cat(as.character(txtleft[i]), "\t", file = zz1)
                      cat(as.character(txtmidle[i]), "\t", file = zz1)
                      cat(as.character(txtright[i]), "\n", file = zz1)
                    }
                    close(zz1)
                    tkmessageBox(message = paste("Wrote: ", outfile, 
                      sep = ""), icon = "info", type = "ok", 
                      parent = tt)
                  }
                }
                tttextpop <- tktoplevel(parent = .TkRoot)
                tkgrab(tttextpop)
                tkwm.title(tttextpop, "Multiple runs for inferring the number of populations")
                tkfocus(tttextpop)
                tkwait.visibility(tttextpop)
                left <- tktext(tttextpop)
                midle <- tktext(tttextpop)
                right <- tktext(tttextpop)
                load <- tktext(tttextpop)
                posx <- tclVar("")
                posy <- tclVar("")
                yscr <- tkscrollbar(tttextpop, repeatinterval = 5, 
                  command = function(...) {
                    tkyview(midle, ...)
                    tkyview(left, ...)
                    tkyview(right, ...)
                  })
                ltop <- tktext(tttextpop, font = tkfont.create(family = "courrier"), 
                  height = 1, width = 15, wrap = "none")
                mtop <- tktext(tttextpop, font = tkfont.create(family = "courrier"), 
                  height = 1, width = 30, wrap = "none")
                rtop <- tktext(tttextpop, font = tkfont.create(family = "courrier"), 
                  height = 1, width = 30, wrap = "none")
                loadtop <- tktext(tttextpop, font = tkfont.create(family = "courrier"), 
                  height = 1, width = 15, wrap = "none")
                tkconfigure(left, font = tkfont.create(family = "courier"), 
                  wrap = "none", width = 15, yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(midle, as.double(...))
                    tkyview.moveto(right, as.double(...))
                    tkyview.moveto(load, as.double(...))
                  })
                tkconfigure(midle, font = tkfont.create(family = "courier"), 
                  wrap = "none", width = 30, yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(left, as.double(...))
                    tkyview.moveto(right, as.double(...))
                    tkyview.moveto(load, as.double(...))
                  })
                tkconfigure(right, font = tkfont.create(family = "courier"), 
                  wrap = "none", width = 30, yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(left, as.double(...))
                    tkyview.moveto(midle, as.double(...))
                    tkyview.moveto(load, as.double(...))
                  })
                tkconfigure(load, font = tkfont.create(family = "courier"), 
                  wrap = "none", width = 15, yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(left, as.double(...))
                    tkyview.moveto(midle, as.double(...))
                    tkyview.moveto(right, as.double(...))
                  })
                sortbutton <- tkbutton(tttextpop, text = "Sort by\nposterior probability", 
                  command = Sort)
                burnbutton <- tkbutton(tttextpop, text = "Recalculate\n with burnin", 
                  command = Reburn)
                burnentry <- tkentry(tttextpop, width = "10", 
                  textvariable = reburnvalue)
                timelabel.widget <- tklabel(tttextpop, text = "...", 
                  foreground = "blue")
                outputbutton <- tkbutton(tttextpop, text = "Save to file", 
                  command = Output)
                tkinsert(ltop, "end", "Run")
                tkinsert(mtop, "end", "Number of populations")
                tkinsert(rtop, "end", "Average log posterior probability")
                tkinsert(loadtop, "end", "Select a run")
                tkgrid(ltop, row = 1, column = 1)
                tkgrid(mtop, row = 1, column = 2, columnspan = 2)
                tkgrid(rtop, row = 1, column = 4)
                tkgrid(loadtop, row = 1, column = 5)
                tkgrid(left, row = 2, column = 1)
                tkgrid(midle, row = 2, column = 2, columnspan = 2)
                tkgrid(right, row = 2, column = 4)
                tkgrid(load, row = 2, column = 5)
                tkgrid(yscr, row = 2, column = 6, sticky = "ns")
                tkgrid(sortbutton, row = 3, column = 4)
                tkgrid(burnentry, row = 3, column = 3, sticky = "w")
                tkgrid(burnbutton, row = 3, column = 2, sticky = "e")
                tkgrid(outputbutton, row = 3, column = 1)
                tkgrid(timelabel.widget, row = 4, column = 1, 
                  columnspan = 4)
                print("Starting...")
                initialtime <- as.numeric(Sys.time(), "secs")
                mr <- function(i) {
                  tempoutputdir <- paste(tclvalue(outputdir), 
                    as.character(i), "/", sep = "")
                  dir.create(tempoutputdir, showWarnings = FALSE)
                  Log(paste("dir.create(", tempoutputdir, ", showWarnings = FALSE)", 
                    sep = ""), "[SUCCESS]")
                  tcl("update")
                  Sys.sleep(0.5)
                  err <- try(MCMC(coordinates = globalcoordinates, 
                    genotypes = globalgenotypes, ploidy = dploidy, 
                    path.mcmc = tempoutputdir, rate.max = as.numeric(tclvalue(rate)), 
                    delta.coord = as.numeric(tclvalue(delta)), 
                    npopmin = as.numeric(tclvalue(npopmin)), 
                    npopinit = as.numeric(tclvalue(npopinit)), 
                    npopmax = as.numeric(tclvalue(npopmax)), 
                    nb.nuclei.max = as.numeric(tclvalue(nuclei)), 
                    nit = as.numeric(tclvalue(nit)), thinning = as.numeric(tclvalue(thinning)), 
                    freq.model = tclvalue(freq), shape1 = as.numeric(tclvalue(gshape1)), 
                    shape2 = as.numeric(tclvalue(gshape2)), varnpop = varnpop, 
                    spatial = as.logical(tclvalue(spatial)), 
                    jcf = as.logical(tclvalue(jcf)), filter.null.alleles = as.logical(tclvalue(null))), 
                    silent = TRUE)
                  print("Done")
                  if (class(err) == "try-error") {
                    Log(paste("MCMC(coordinates=", matrix2str(globalcoordinates), 
                      ",genotypes=", matrix2str(globalgenotypes), 
                      ",ploidy=", dploidy, ",path.mcmc=\"", tempoutputdir, 
                      "\",rate.max=", as.numeric(tclvalue(rate)), 
                      ",delta.coord=", as.numeric(tclvalue(delta)), 
                      ",npopmin=", as.numeric(tclvalue(npopmin)), 
                      ",npopinit=", as.numeric(tclvalue(npopinit)), 
                      ",npopmax=", as.numeric(tclvalue(npopmax)), 
                      ",nb.nuclei.max=", as.numeric(tclvalue(nuclei)), 
                      ",nit=", as.numeric(tclvalue(nit)), ",thinning=", 
                      as.numeric(tclvalue(thinning)), ",freq.model=\"", 
                      tclvalue(freq), "\",shape1=", shape1 = as.numeric(tclvalue(gshape1)), 
                      ",shape2=", shape2 = as.numeric(tclvalue(gshape2)), 
                      ",varnpop=", varnpop, ",spatial=", as.logical(tclvalue(spatial)), 
                      ",jcf=", as.logical(tclvalue(jcf)), ",filter.null.alleles =", 
                      as.logical(tclvalue(null)), ")", sep = ""), 
                      "[FAILED] ")
                    tkinsert(left, "end", as.character(i))
                    tkinsert(left, "end", "\n")
                    tkinsert(midle, "end", "failed\n")
                    tkinsert(right, "end", "failed\n")
                    runs <- c(runs, i)
                    pops <<- c(pops, NA)
                    probs <<- c(probs, NA)
                  }
                  else {
                    Log(paste("MCMC(coordinates=", matrix2str(globalcoordinates), 
                      ",genotypes=", matrix2str(globalgenotypes), 
                      ",ploidy=", dploidy, ",path.mcmc=\"", tempoutputdir, 
                      "\",rate.max=", as.numeric(tclvalue(rate)), 
                      ",delta.coord=", as.numeric(tclvalue(delta)), 
                      ",npopmin=", as.numeric(tclvalue(npopmin)), 
                      ",npopinit=", as.numeric(tclvalue(npopinit)), 
                      ",npopmax=", as.numeric(tclvalue(npopmax)), 
                      ",nb.nuclei.max=", as.numeric(tclvalue(nuclei)), 
                      ",nit=", as.numeric(tclvalue(nit)), ",thinning=", 
                      as.numeric(tclvalue(thinning)), ",freq.model=\"", 
                      tclvalue(freq), "\",shape1=", shape1 = as.numeric(tclvalue(gshape1)), 
                      ",shape2=", shape2 = as.numeric(tclvalue(gshape2)), 
                      ",varnpop=", varnpop, ",spatial=", as.logical(tclvalue(spatial)), 
                      ",jcf=", as.logical(tclvalue(jcf)), ",filter.null.alleles =", 
                      as.logical(tclvalue(null)), ")", sep = ""), 
                      "[SUCCESS] ")
                    ltext = tklabel(tttextpop, text = as.character(i), 
                      border = 3)
                    tkwindow.create(left, "end", window = ltext)
                    tkinsert(left, "end", "\n")
                    runs <<- c(runs, i)
                    file <- try(scan(paste(tempoutputdir, "populations.numbers.txt", 
                      sep = "")), silent = TRUE)
                    dist <- hist(file, plot = FALSE, breaks = seq(0.5, 
                      max(file) + 0.5, 1))
                    firsttime <- 0
                    straux <- ""
                    for (j in 1:length(dist$counts)) {
                      if (dist$counts[j] == max(dist$counts)) {
                        if (firsttime == 0) {
                          straux <- as.character(dist$mids[j])
                          firsttime <- 1
                        }
                        else {
                          straux <- paste(straux, " and ", sep = "")
                          straux <- paste(straux, as.character(dist$mids[j]), 
                            sep = "")
                        }
                      }
                    }
                    straux <- paste(straux, " ( ", sep = "")
                    straux <- paste(straux, as.character(as.double(max(dist$density) * 
                      100)), sep = "")
                    straux <- paste(straux, " %) ", sep = "")
                    lmiddle = tklabel(tttextpop, text = straux, 
                      border = 3)
                    tkwindow.create(midle, "end", window = lmiddle)
                    pops <<- c(pops, straux)
                    tkinsert(midle, "end", "\n")
                    file <- try(scan(paste(tempoutputdir, "log.posterior.density.txt", 
                      sep = "")), silent = TRUE)
                    mpd <- mean(file)
                    lright = tklabel(tttextpop, text = mpd, border = 3)
                    tkwindow.create(right, "end", window = lright)
                    tkinsert(right, "end", "\n")
                    probs <<- c(probs, mpd)
                    if (i == 1) {
                      runtime <<- as.numeric(Sys.time(), "secs") - 
                        initialtime
                      totaltime <<- runtime * as.integer(tclvalue(ntestpop))
                    }
                    makebutton(i)
                    tkinsert(load, "end", "\n")
                    changetotime <- function() {
                      seconds <- (totaltime - runtime)%%60
                      aux <- (totaltime - runtime)%/%60
                      minutes <- aux%%60
                      aux <- aux%/%60
                      hours <- aux%%24
                      aux <- aux%/%24
                      str <- ""
                      if (aux != 0) 
                        str <- paste(str, as.integer(aux), " day(s), ", 
                          sep = "")
                      if (hours != 0) 
                        str <- paste(str, as.integer(hours), 
                          " hour(s), ", sep = "")
                      if (minutes != 0) 
                        str <- paste(str, as.integer(minutes), 
                          " minute(s), ", sep = "")
                      str <- paste(str, as.integer(seconds), 
                        " second(s)", sep = "")
                      return(str)
                    }
                    tkconfigure(timelabel.widget, text = paste("about ", 
                      changetotime(), " remaining", sep = ""))
                    totaltime <<- totaltime - runtime
                    txtleft <- c(txtleft, as.character(i))
                    txtmidle <- c(txtmidle, straux)
                    txtright <- c(txtright, mpd)
                  }
                  tkyview.moveto(left, 1)
                }
                mrcluster <- function(i, outdir, crate.max, cdelta.coord, 
                  cnpopmin, cnpopinit, cnpopmax, cnb.nuclei.max, 
                  cnit, cthinning, cfreq.model, cshape1, cshape2, 
                  cspatial, cjcf, cfilter.null.alleles) {
                  tempoutputdir <- paste(outdir, as.character(i), 
                    "/", sep = "")
                  dir.create(tempoutputdir, showWarnings = FALSE)
                  require("Geneland")
                  Sys.sleep(0.5)
                  err <- try(MCMC(coordinates = globalcoordinates, 
                    genotypes = globalgenotypes, ploidy = dploidy, 
                    path.mcmc = tempoutputdir, rate.max = crate.max, 
                    delta.coord = cdelta.coord, npopmin = cnpopmin, 
                    npopinit = cnpopinit, npopmax = cnpopmax, 
                    nb.nuclei.max = cnb.nuclei.max, nit = cnit, 
                    thinning = cthinning, freq.model = cfreq.model, 
                    shape1 = cshape1, shape2 = cshape2, varnpop = varnpop, 
                    spatial = cspatial, jcf = cjcf, filter.null.alleles = cfilter.null.alleles), 
                    silent = TRUE)
                  zz1 <- file(paste(tempoutputdir, "ClusterLog.txt", 
                    sep = ""), "a")
                  cat(err, "\n", file = zz1)
                  close(zz1)
                  print("Done")
                }
                mrafter <- function(i) {
                  tempoutputdir <- paste(tclvalue(outputdir), 
                    as.character(i), "/", sep = "")
                  ltext = tklabel(tttextpop, text = as.character(i), 
                    border = 3)
                  tkwindow.create(left, "end", window = ltext)
                  tkinsert(left, "end", "\n")
                  runs <<- c(runs, i)
                  file <- try(scan(paste(tempoutputdir, "populations.numbers.txt", 
                    sep = "")), silent = TRUE)
                  dist <- hist(file, plot = FALSE, breaks = seq(0.5, 
                    max(file) + 0.5, 1))
                  firsttime <- 0
                  straux <- ""
                  for (j in 1:length(dist$counts)) {
                    if (dist$counts[j] == max(dist$counts)) {
                      if (firsttime == 0) {
                        straux <- as.character(dist$mids[j])
                        firsttime <- 1
                      }
                      else {
                        straux <- paste(straux, " and ", sep = "")
                        straux <- paste(straux, as.character(dist$mids[j]), 
                          sep = "")
                      }
                    }
                  }
                  straux <- paste(straux, " ( ", sep = "")
                  straux <- paste(straux, as.character(as.double(max(dist$density) * 
                    100)), sep = "")
                  straux <- paste(straux, " %) ", sep = "")
                  lmiddle = tklabel(tttextpop, text = straux, 
                    border = 3)
                  tkwindow.create(midle, "end", window = lmiddle)
                  pops <<- c(pops, straux)
                  tkinsert(midle, "end", "\n")
                  file <- try(scan(paste(tempoutputdir, "log.posterior.density.txt", 
                    sep = "")), silent = TRUE)
                  mpd <- mean(file)
                  lright = tklabel(tttextpop, text = mpd, border = 3)
                  tkwindow.create(right, "end", window = lright)
                  tkinsert(right, "end", "\n")
                  probs <<- c(probs, mpd)
                  if (i == 1) {
                    runtime <<- as.numeric(Sys.time(), "secs") - 
                      initialtime
                    totaltime <<- runtime * as.integer(tclvalue(ntestpop))
                  }
                  makebutton(i)
                  tkinsert(load, "end", "\n")
                  changetotime <- function() {
                    seconds <- (totaltime - runtime)%%60
                    aux <- (totaltime - runtime)%/%60
                    minutes <- aux%%60
                    aux <- aux%/%60
                    hours <- aux%%24
                    aux <- aux%/%24
                    str <- ""
                    if (aux != 0) 
                      str <- paste(str, as.integer(aux), " day(s), ", 
                        sep = "")
                    if (hours != 0) 
                      str <- paste(str, as.integer(hours), " hour(s), ", 
                        sep = "")
                    if (minutes != 0) 
                      str <- paste(str, as.integer(minutes), 
                        " minute(s), ", sep = "")
                    str <- paste(str, as.integer(seconds), " second(s)", 
                      sep = "")
                    return(str)
                  }
                  tkconfigure(timelabel.widget, text = paste("about ", 
                    changetotime(), " remaining", sep = ""))
                  totaltime <<- totaltime - runtime
                  txtleft <- c(txtleft, as.character(i))
                  txtmidle <- c(txtmidle, straux)
                  txtright <- c(txtright, mpd)
                  tkyview.moveto(left, 1)
                }
                tcl("update")
                if (usecluster) {
                  tkconfigure(timelabel.widget, text = "Parallel processing...")
                  tcl("update")
                  clusterApply(cluster, 1:tclvalue(ntestpop), 
                    mrcluster, tclvalue(outputdir), as.numeric(tclvalue(rate)), 
                    as.numeric(tclvalue(delta)), as.numeric(tclvalue(npopmin)), 
                    as.numeric(tclvalue(npopinit)), as.numeric(tclvalue(npopmax)), 
                    as.numeric(tclvalue(nuclei)), as.numeric(tclvalue(nit)), 
                    cthinning <- as.numeric(tclvalue(thinning)), 
                    tclvalue(freq), as.numeric(tclvalue(gshape1)), 
                    as.numeric(tclvalue(gshape2)), as.logical(tclvalue(spatial)), 
                    as.logical(tclvalue(jcf)), as.logical(tclvalue(null)))
                  for (i in 1:as.numeric(tclvalue(ntestpop))) mrafter(i)
                }
                else {
                  for (i in 1:as.numeric(tclvalue(ntestpop))) mr(i)
                }
                tkconfigure(timelabel.widget, text = "Done")
                if (tclvalue(freq) == "Correlated") 
                  falush <<- 1
                else falush <<- 0
                tkdestroy(tttry)
                tkfocus(tttextpop)
                tkgrab("release", tttextpop)
            }
        }
        ploidylabel.widget <- tklabel(ttrun, text = "Ploidy:")
        ploidy <- tclVar("Diploid")
        wdiplody <- .Tk.subwin(ttrun)
        ploidyoptionmenu.widget <- tcl("tk_optionMenu", wdiplody, 
            ploidy, "Diploid", "Haploid")
        tkgrid(ploidylabel.widget, row = 2, column = 1, sticky = "w")
        tkgrid(wdiplody, row = 2, column = 3, columnspan = 3, 
            sticky = "w")
        rate <- tclVar(100)
        rate.widget <- tkentry(ttrun, width = "23", textvariable = rate)
        ratelabel.widget <- tklabel(ttrun, text = "Maximum rate of poisson process:")
        delta <- tclVar(0)
        delta.widget <- tkentry(ttrun, width = "23", textvariable = delta)
        deltalabel.widget <- tklabel(ttrun, text = "Uncertainty on coordinate:")
        tkgrid(deltalabel.widget, row = 4, column = 1, sticky = "w")
        tkgrid(delta.widget, row = 4, column = 3, columnspan = 3, 
            sticky = "w")
        poplabel.widget <- tklabel(ttrun, text = "Number of populations:")
        npopmin <- tclVar(1)
        npopinit <- tclVar(1)
        npopmax <- tclVar(1)
        wmin <- .Tk.subwin(ttrun)
        winit <- .Tk.subwin(ttrun)
        wmax <- .Tk.subwin(ttrun)
        actualizaspin <- function() {
            tkconfigure(min, "-to", as.numeric(tclvalue(npopmax)))
            tkconfigure(max, "-from", as.numeric(tclvalue(npopmin)))
            tkconfigure(init, "-from", as.numeric(tclvalue(npopmin)), 
                "-to", as.numeric(tclvalue(npopmax)))
            if (as.numeric(tclvalue(npopinit)) > as.numeric(tclvalue(npopmax))) {
                tclvalue(npopinit) <- tclvalue(npopmax)
                tkconfigure(init, "-textvariable", npopinit)
            }
            else if (as.numeric(tclvalue(npopinit)) < as.numeric(tclvalue(npopmin))) {
                tclvalue(npopinit) <- tclvalue(npopmin)
                tkconfigure(init, "-textvariable", npopinit)
            }
        }
        min <- tcl("spinbox", wmin, "-textvariable", npopmin, 
            "-width", 5, "-increment", 1, "-from", 1, "-to", 
            as.numeric(tclvalue(npopmax)), "-command", actualizaspin)
        init <- tcl("spinbox", winit, "-textvariable", npopinit, 
            "-width", 5, "-increment", 1, "-from", as.numeric(tclvalue(npopmin)), 
            "-to", as.numeric(tclvalue(npopmax)), "-command", 
            actualizaspin)
        max <- tcl("spinbox", wmax, "-textvariable", npopmax, 
            "-width", 5, "-increment", 1, "-from", as.numeric(tclvalue(npopmin)), 
            "-to", 1e+05, "-command", actualizaspin)
        npopminlabel.widget <- tklabel(ttrun, text = "pop min")
        npopinitlabel.widget <- tklabel(ttrun, text = "pop init")
        npopmaxlabel.widget <- tklabel(ttrun, text = "pop max")
        tkgrid(npopminlabel.widget, row = 5, column = 3, sticky = "w")
        tkgrid(npopinitlabel.widget, row = 5, column = 4, sticky = "w")
        tkgrid(npopmaxlabel.widget, row = 5, column = 5, sticky = "w")
        tkgrid(poplabel.widget, row = 6, column = 1, sticky = "w")
        tkgrid(wmin, row = 6, column = 3, sticky = "w")
        tkgrid(winit, row = 6, column = 4, sticky = "w")
        tkgrid(wmax, row = 6, column = 5, sticky = "w")
        nuclei <- tclVar(300)
        nuclei.widget <- tkentry(ttrun, width = "23", textvariable = nuclei)
        nucleilabel.widget <- tklabel(ttrun, text = "Maximum number of nuclei:")
        nit <- tclVar()
        nit.widget <- tkentry(ttrun, width = "23", textvariable = nit)
        nitlabel.widget <- tklabel(ttrun, text = "Number of iterations:")
        tkgrid(nitlabel.widget, row = 8, column = 1, sticky = "w")
        tkgrid(nit.widget, row = 8, column = 3, columnspan = 3, 
            sticky = "w")
        thinning <- tclVar()
        thinning.widget <- tkentry(ttrun, width = "23", textvariable = thinning)
        thinninglabel.widget <- tklabel(ttrun, text = "Thinning:")
        tkgrid(thinninglabel.widget, row = 9, column = 1, sticky = "w")
        tkgrid(thinning.widget, row = 9, column = 3, columnspan = 3, 
            sticky = "w")
        freqlabel.widget <- tklabel(ttrun, text = "Allele frequency model:")
        freq <- tclVar("Uncorrelated")
        wfreq <- .Tk.subwin(ttrun)
        freqoptionmenu.widget <- tcl("tk_optionMenu", wfreq, 
            freq, "Correlated", "Uncorrelated")
        gshape1 <- tclVar(2)
        gshape2 <- tclVar(18)
        shape1.widget <- tkentry(ttrun, width = "6", textvariable = gshape1)
        shape2.widget <- tkentry(ttrun, width = "6", textvariable = gshape2)
        shapelabel.widget <- tklabel(ttrun, text = "Drift factor prior-gamma[a,b] :")
        spatiallabel.widget <- tklabel(ttrun, text = "Spatial model:")
        spatial <- tclVar("TRUE")
        wspatial <- .Tk.subwin(ttrun)
        spatialoptionmenu.widget <- tcl("tk_optionMenu", wspatial, 
            spatial, "FALSE", "TRUE")
        tkgrid(spatiallabel.widget, row = 12, column = 1, sticky = "w")
        tkgrid(wspatial, row = 12, column = 3, columnspan = 3, 
            sticky = "w")
        nulllabel.widget <- tklabel(ttrun, text = "Null allele model:")
        null <- tclVar("TRUE")
        wnull <- .Tk.subwin(ttrun)
        nulloptionmenu.widget <- tcl("tk_optionMenu", wnull, 
            null, "FALSE", "TRUE")
        tkgrid(nulllabel.widget, row = 13, column = 1, sticky = "w")
        tkgrid(wnull, row = 13, column = 3, columnspan = 3, sticky = "w")
        jcf <- tclVar("FALSE")
        labelspace0 <- tklabel(ttrun, text = " ")
        tkgrid(labelspace0, row = 14, column = 1)
        labelspace <- tklabel(ttrun, image = imageline)
        tkgrid(labelspace, row = 15, column = 1, columnspan = 6, 
            sticky = "news")
        labelspace1 <- tklabel(ttrun, text = " ")
        tkgrid(labelspace1, row = 16, column = 1)
        activatetestpop <- function() {
            if (tclvalue(testnumberpop) == 1) 
                tkconfigure(testpop, state = "normal")
            else tkconfigure(testpop, state = "disable")
        }
        ttestpoplabel.widget <- tklabel(ttrun, text = "Multiple independent runs:")
        ttestpoplabelyes.widget <- tklabel(ttrun, text = "No")
        ttestpoplabelno.widget <- tklabel(ttrun, text = "Yes")
        testpopyes.widget <- tkradiobutton(ttrun, command = activatetestpop, 
            variable = testnumberpop, value = 0, selectcolor = "white")
        testpopno.widget <- tkradiobutton(ttrun, command = activatetestpop, 
            variable = testnumberpop, value = 1, selectcolor = "white")
        ntestpop <- tclVar(1)
        wtestpop <- .Tk.subwin(ttrun)
        testpop <- tcl("spinbox", wtestpop, "-textvariable", 
            ntestpop, "-width", 5, "-increment", 1, "-from", 
            1, "-to", 999, state = "disable")
        tkconfigure(testpop, state = "disable")
        tkgrid(ttestpoplabel.widget, row = 17, column = 1, columnspan = 2, 
            rowspan = 2, sticky = "w")
        tkgrid(ttestpoplabelyes.widget, row = 17, column = 3, 
            sticky = "w")
        tkgrid(ttestpoplabelno.widget, row = 17, column = 4, 
            sticky = "w")
        tkgrid(testpopyes.widget, row = 18, column = 3, sticky = "w")
        tkgrid(testpopno.widget, row = 18, column = 4, sticky = "w")
        tkgrid(wtestpop, row = 18, column = 5, sticky = "w")
        labelspace2 <- tklabel(ttrun, text = " ")
        tkgrid(labelspace2, row = 19, column = 1)
        labelspace3 <- tklabel(ttrun, text = "             ")
        labelspace4 <- tklabel(ttrun, text = "             ")
        nextbutton <- tkbutton(ttrun, image = imagerun2, text = "RUN >>", 
            command = RunmcmcFmodel)
        tkgrid(nextbutton, row = 20, column = 3, columnspan = 3, 
            sticky = "e")
        tkfocus(ttrun)
        tkgrid(ratelabel.widget, row = 3, column = 1, sticky = "w")
        tkgrid(rate.widget, row = 3, column = 3, columnspan = 3, 
            sticky = "w")
        tkgrid(nucleilabel.widget, row = 7, column = 1, sticky = "w")
        tkgrid(nuclei.widget, row = 7, column = 3, columnspan = 3, 
            sticky = "w")
        tkgrid(freqlabel.widget, row = 10, column = 1, sticky = "w")
        tkgrid(wfreq, row = 10, column = 3, columnspan = 3, sticky = "w")
        tkgrid(shapelabel.widget, row = 11, column = 1, sticky = "w")
        tkgrid(shape1.widget, row = 11, column = 3, sticky = "w")
        tkgrid(shape2.widget, row = 11, column = 4, sticky = "w")
        if (tclvalue(advanced) == 1) {
            tkconfigure(ratelabel.widget, state = "normal")
            tkconfigure(rate.widget, state = "normal")
            tkconfigure(nucleilabel.widget, state = "normal")
            tkconfigure(nuclei.widget, state = "normal")
            tkconfigure(shapelabel.widget, state = "normal")
            tkconfigure(shape1.widget, state = "normal")
            tkconfigure(shape2.widget, state = "normal")
            tkgrid.remove(labelspace3)
            tkgrid(winit, row = 6, column = 4, sticky = "w")
            tkconfigure(npopinitlabel.widget, state = "normal")
        }
        else {
            tkconfigure(ratelabel.widget, state = "disable")
            tkconfigure(rate.widget, state = "disable")
            tkconfigure(nucleilabel.widget, state = "disable")
            tkconfigure(nuclei.widget, state = "disable")
            tkconfigure(shapelabel.widget, state = "disable")
            tkconfigure(shape1.widget, state = "disable")
            tkconfigure(shape2.widget, state = "disable")
            tkconfigure(winit, relief = "raised")
            tkconfigure(npopinitlabel.widget, state = "disable")
            tkgrid.remove(winit)
            tkgrid(labelspace3, row = 6, column = 4, sticky = "w")
        }
    }
    initialimage <- function() {
        imgAsLabel <- tklabel(ttinit, image = image1, bg = "white")
        tkgrid(imgAsLabel, sticky = "news")
        notice <- tklabel(ttinit, text = "Geneland is loaded\n\n* Please *\n\nRegister on http://folk.uio.no/gillesg/Geneland/register.php")
        tkbind(notice, "<Button-1>", function() {
            browseURL("http://folk.uio.no/gillesg/Geneland/register.php")
        })
        tkgrid(notice, sticky = "news")
    }
    postproc <- function() {
        RunpostprocessChain <- function() {
            tttry <- tktoplevel(parent = .TkRoot)
            tkgrab(tttry)
            tkwm.geometry(tttry, "+200+200")
            tkwm.title(tttry, "wait")
            warn <- tklabel(tttry, image = imagepleasewait)
            tkpack(warn)
            tkfocus(tttry)
            tcl("update")
            Sys.sleep(0.5)
            print("Starting...")
            err <- try(PostProcessChain(coordinates = globalcoordinates, 
                genotypes = globalgenotypes, path.mcmc = tclvalue(outputdir), 
                nxdom = as.numeric(tclvalue(nxdom)), nydom = as.numeric(tclvalue(nydom)), 
                burnin = as.numeric(tclvalue(burnin))), silent = TRUE)
            tkdestroy(tttry)
            print("Done.")
            if (class(err) == "try-error") {
                Log(paste("PostProcessChain(coordinates=", matrix2str(globalcoordinates), 
                  ",genotypes=", matrix2str(globalgenotypes), 
                  ",path.mcmc=\"", tclvalue(outputdir), "\",nxdom=", 
                  as.numeric(tclvalue(nxdom)), ",nydom=", as.numeric(tclvalue(nydom)), 
                  ",burnin=", as.numeric(tclvalue(burnin)), ")", 
                  sep = ""), "[FAILED] ")
                tkmessageBox(message = err, icon = "error", type = "ok", 
                  parent = tt)
            }
            else {
                tkmessageBox(message = "Terminated with success", 
                  type = "ok", parent = tt)
                Log(paste("PostProcessChain(coordinates=", matrix2str(globalcoordinates), 
                  ",genotypes=", matrix2str(globalgenotypes), 
                  ",path.mcmc=\"", tclvalue(outputdir), "\",nxdom=", 
                  as.numeric(tclvalue(nxdom)), ",nydom=", as.numeric(tclvalue(nydom)), 
                  ",burnin=", as.numeric(tclvalue(burnin)), ")", 
                  sep = ""), "[SUCCESS] ")
            }
        }
        Nullalleles <- function() {
            Output <- function() {
                outfile <- tclvalue(tkgetSaveFile(filetypes = "{{.txt} *.txt}", 
                  title = "Save to file"))
                if (outfile != "") {
                  zz1 <- file(outfile, "w")
                  for (i in 1:length(err)) {
                    cat(as.character(names(err[i])), "\t", file = zz1)
                    cat(as.character(err[i]), "\n", file = zz1)
                  }
                  close(zz1)
                  tkmessageBox(message = paste("Wrote: ", outfile, 
                    sep = ""), icon = "info", type = "ok", parent = tt)
                }
            }
            print("Starting...")
            err <- try(EstimateFreqNA(genotypes = globalgenotypes, 
                path.mcmc = tclvalue(outputdir)), silent = TRUE)
            if (class(err) == "try-error") {
                Log(paste("EstimateFreqNA(genotypes=", matrix2str(globalgenotypes), 
                  ",path.mcmc=\"", tclvalue(outputdir), "\")", 
                  sep = ""), "[FAILED] ")
                tkmessageBox(message = err, icon = "error", type = "ok", 
                  parent = tt)
            }
            else {
                Log(paste("EstimateFreqNA(genotypes=", matrix2str(globalgenotypes), 
                  ",path.mcmc=\"", tclvalue(outputdir), "\")", 
                  sep = ""), "[SUCCESS] ")
                tkgrid.remove(nabutton)
                tttextna <- tktoplevel(parent = .TkRoot)
                tkwm.title(tttextna, "Estimated frequency of null alelles")
                posx <- tclVar("")
                posy <- tclVar("")
                left <- tktext(tttextna, width = 30)
                txt <- tktext(tttextna, width = 30)
                yscr <- tkscrollbar(tttextna, repeatinterval = 5, 
                  command = function(...) {
                    tkyview(txt, ...)
                    tkyview(left, ...)
                  })
                tkconfigure(txt, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(left, as.double(...))
                  })
                tkconfigure(left, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(txt, as.double(...))
                  })
                auxtxt <- ""
                auxleft <- ""
                for (i in 1:length(err)) {
                  auxleft <- paste(auxleft, names(err[i]), "\n", 
                    sep = "")
                  auxtxt <- paste(auxtxt, as.character(err[i]), 
                    "\n", sep = "")
                }
                tkinsert(left, "end", auxleft)
                tkinsert(txt, "end", auxtxt)
                tkgrid(txt, row = 2, column = 2)
                tkgrid(left, row = 2, column = 1)
                tkgrid(yscr, row = 2, column = 3, sticky = "ns")
                outputbutton <- tkbutton(tttextna, text = "Save to file", 
                  command = Output)
                tkgrid(outputbutton, row = 3, column = 1, columnspan = 2)
                print(err)
                tcl("update")
            }
            print("Done.")
        }
        ndomlabel.widget <- tklabel(ttpost, text = "Number of pixels in the spatial domain:")
        tkgrid(ndomlabel.widget, row = 1, column = 1, columnspan = 2, 
            sticky = "w")
        nxdom <- tclVar(50)
        nydom <- tclVar(50)
        nxdomlabel.widget <- tklabel(ttpost, text = "X")
        nydomlabel.widget <- tklabel(ttpost, text = "Y")
        tkgrid(nxdomlabel.widget, row = 2, column = 1, sticky = "w")
        tkgrid(nydomlabel.widget, row = 2, column = 2, sticky = "w")
        nxdom.widget <- tkentry(ttpost, width = "20", textvariable = nxdom)
        nydom.widget <- tkentry(ttpost, width = "20", textvariable = nydom)
        tkgrid(nxdom.widget, row = 3, column = 1, sticky = "w")
        tkgrid(nydom.widget, row = 3, column = 2, sticky = "w")
        burnin.widget <- tkentry(ttpost, width = "20", textvariable = burnin)
        burninlabel.widget <- tklabel(ttpost, text = "Burnin:")
        tkgrid(burninlabel.widget, row = 4, column = 1, sticky = "w")
        tkgrid(burnin.widget, row = 4, column = 2, sticky = "w")
        labelspace <- tklabel(ttpost, text = " ")
        tkgrid(labelspace, row = 5, column = 1)
        nabutton <- tkbutton(ttpost, text = "Click here to estimate the frequency of null alleles", 
            command = Nullalleles)
        tkgrid(nabutton, row = 6, column = 1, columnspan = 2, 
            sticky = "w")
        labelspace2 <- tklabel(ttpost, text = " ")
        tkgrid(labelspace2, row = 7, column = 1)
        nextbutton <- tkbutton(ttpost, image = imagerun2, text = "RUN >>", 
            command = RunpostprocessChain)
        tkgrid(nextbutton, row = 8, column = 2, sticky = "e")
        tkfocus(ttpost)
    }
    GraficalIBD <- function() {
        if (length(idb.dataset) == 1) {
            tkmessageBox(message = "First simulate some data", 
                icon = "error", type = "ok", parent = tt)
        }
        else {
            DrawShowIBD <- function() {
                vect1 <- c()
                vec1 <- unlist(strsplit(tclvalue(loc.grid), ","))
                for (i in 1:length(vec1)) vect1[i] <- as.numeric(vec1[i])
                if (tclvalue(plot.coord.ps) == "") 
                  file.plot.coord <- NA
                else file.plot.coord <- tclvalue(plot.coord.ps)
                if (tclvalue(plot.tess.ps) == "") 
                  file.plot.tess <- NA
                else file.plot.tess <- tclvalue(plot.tess.ps)
                if (tclvalue(plot.freq.grid.ps) == "") 
                  file.plot.freq.grid <- NA
                else file.plot.freq.grid <- tclvalue(plot.freq.grid.ps)
                if (tclvalue(plot.freq.indiv.ps) == "") 
                  file.plot.freq.indiv <- NA
                else file.plot.freq.indiv <- tclvalue(plot.freq.indiv.ps)
                if (tclvalue(plot.gen.ps) == "") 
                  file.plot.gen <- NA
                else file.plot.gen <- tclvalue(plot.gen.ps)
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                Sys.sleep(0.5)
                print("Starting...")
                err <- try(show.simdata(dataset = idb.dataset, 
                  plot.coord = as.logical(tclvalue(plot.coord)), 
                  file.plot.coord = file.plot.coord, plot.tess = as.logical(tclvalue(plot.tess)), 
                  file.plot.tess = file.plot.tess, plot.freq.grid = as.logical(tclvalue(plot.freq.grid)), 
                  file.plot.freq.grid = file.plot.freq.grid, 
                  plot.freq.indiv = as.logical(tclvalue(plot.freq.indiv)), 
                  file.plot.freq.indiv = file.plot.freq.indiv, 
                  loc.grid = vect1, loc.indiv = as.numeric(tclvalue(loc.indiv)), 
                  zlim.freq = c(as.numeric(tclvalue(zlimin)), 
                    as.numeric(tclvalue(zlimax))), plot.gen = as.logical(tclvalue(plot.gen)), 
                  file.plot.gen = file.plot.gen), silent = TRUE)
                tkdestroy(tttry)
                print("Done")
                if (class(err) == "try-error") {
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                }
            }
            ttshowibd <- tktoplevel()
            tkwm.title(ttshowibd, "Graphical Display of data simulated by IBD")
            plot.coord.ps <- tclVar("")
            plot.tess.ps <- tclVar("")
            plot.freq.grid.ps <- tclVar("")
            plot.freq.indiv.ps <- tclVar("")
            plot.gen.ps <- tclVar("")
            getprintfile <- function() {
                printfile <- tclVar()
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttshowibd)
                return(printfile)
            }
            plot.coord.psbutton.widget <- tkbutton(ttshowibd, 
                text = "Save File", command = function() {
                  plot.coord.ps <<- getprintfile()
                  tkconfigure(plot.coord.pslabel.widget, text = tclvalue(plot.coord.ps))
                }, width = 15)
            plot.coord.pslabel.widget <- tklabel(ttshowibd, text = tclvalue(plot.coord.ps), 
                width = 50)
            plot.tess.psbutton.widget <- tkbutton(ttshowibd, 
                text = "Save File", command = function() {
                  plot.tess.ps <<- getprintfile()
                  tkconfigure(plot.tess.pslabel.widget, text = tclvalue(plot.tess.ps))
                }, width = 15)
            plot.tess.pslabel.widget <- tklabel(ttshowibd, text = tclvalue(plot.tess.ps), 
                width = 50)
            plot.freq.grid.psbutton.widget <- tkbutton(ttshowibd, 
                text = "Save File", command = function() {
                  plot.freq.grid.ps <<- getprintfile()
                  tkconfigure(plot.freq.grid.pslabel.widget, 
                    text = tclvalue(plot.freq.grid.ps))
                }, width = 15)
            plot.freq.grid.pslabel.widget <- tklabel(ttshowibd, 
                text = tclvalue(plot.freq.grid.ps), width = 50)
            plot.freq.indiv.psbutton.widget <- tkbutton(ttshowibd, 
                text = "Save File", command = function() {
                  plot.freq.indiv.ps <<- getprintfile()
                  tkconfigure(plot.freq.indiv.pslabel.widget, 
                    text = tclvalue(plot.freq.indiv.ps))
                }, width = 15)
            plot.freq.indiv.pslabel.widget <- tklabel(ttshowibd, 
                text = tclvalue(plot.freq.indiv.ps), width = 50)
            plot.gen.psbutton.widget <- tkbutton(ttshowibd, text = "Save File", 
                command = function() {
                  plot.gen.ps <<- getprintfile()
                  tkconfigure(plot.gen.pslabel.widget, text = tclvalue(plot.gen.ps))
                }, width = 15)
            plot.gen.pslabel.widget <- tklabel(ttshowibd, text = tclvalue(plot.gen.ps), 
                width = 50)
            plot.coordlabel.widget <- tklabel(ttshowibd, text = "Plot coordinate of individuals:")
            plot.coord <- tclVar("FALSE")
            wplot.coord <- .Tk.subwin(ttshowibd)
            plot.coordoptionmenu.widget <- tcl("tk_optionMenu", 
                wplot.coord, plot.coord, "FALSE", "TRUE")
            plot.tesslabel.widget <- tklabel(ttshowibd, text = "Plot tessellation:")
            plot.tess <- tclVar("FALSE")
            wplot.tess <- .Tk.subwin(ttshowibd)
            plot.tessoptionmenu.widget <- tcl("tk_optionMenu", 
                wplot.tess, plot.tess, "FALSE", "TRUE")
            plot.freq.gridlabel.widget <- tklabel(ttshowibd, 
                text = "Plot allele frequencies for all pixels:")
            plot.freq.grid <- tclVar("FALSE")
            wplot.freq.grid <- .Tk.subwin(ttshowibd)
            plot.freq.gridoptionmenu.widget <- tcl("tk_optionMenu", 
                wplot.freq.grid, plot.freq.grid, "FALSE", "TRUE")
            plot.freq.indivlabel.widget <- tklabel(ttshowibd, 
                text = "Plot allele frequencies at individual sites:")
            plot.freq.indiv <- tclVar("FALSE")
            wplot.freq.indiv <- .Tk.subwin(ttshowibd)
            plot.freq.indivoptionmenu.widget <- tcl("tk_optionMenu", 
                wplot.freq.indiv, plot.freq.indiv, "FALSE", "TRUE")
            loc.grid <- tclVar(1)
            loc.grid.widget <- tkentry(ttshowibd, width = "20", 
                textvariable = loc.grid)
            loc.gridlabel.widget <- tklabel(ttshowibd, text = "Indices of loci (E.g: 1,2,4,...):")
            loc.indiv <- tclVar(1)
            loc.indiv.widget <- tkentry(ttshowibd, width = "20", 
                textvariable = loc.indiv)
            loc.indivlabel.widget <- tklabel(ttshowibd, text = "Indices of loci (E.g: 1,3,6,...):")
            zlimin <- tclVar(0)
            zlimax <- tclVar(1)
            plot.genlabel.widget <- tklabel(ttshowibd, text = "Plot genotype:")
            plot.gen <- tclVar("FALSE")
            wplot.gen <- .Tk.subwin(ttshowibd)
            plot.genoptionmenu.widget <- tcl("tk_optionMenu", 
                wplot.gen, plot.gen, "FALSE", "TRUE")
            tkgrid(plot.coordlabel.widget, row = 2, column = 1, 
                sticky = "w")
            tkgrid(wplot.coord, row = 2, column = 2, sticky = "w")
            tkgrid(plot.coord.psbutton.widget, row = 2, column = 3, 
                sticky = "e")
            tkgrid(plot.coord.pslabel.widget, row = 2, column = 4, 
                sticky = "w")
            tkgrid(plot.tesslabel.widget, row = 3, column = 1, 
                sticky = "w")
            tkgrid(wplot.tess, row = 3, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(plot.tess.psbutton.widget, row = 3, column = 3, 
                sticky = "e")
            tkgrid(plot.tess.pslabel.widget, row = 3, column = 4, 
                sticky = "w")
            tkgrid(plot.freq.gridlabel.widget, row = 4, column = 1, 
                sticky = "w")
            tkgrid(wplot.freq.grid, row = 4, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(plot.freq.grid.psbutton.widget, row = 4, column = 3, 
                sticky = "e")
            tkgrid(plot.freq.grid.pslabel.widget, row = 4, column = 4, 
                sticky = "w")
            tkgrid(plot.freq.indivlabel.widget, row = 6, column = 1, 
                sticky = "w")
            tkgrid(wplot.freq.indiv, row = 6, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(plot.freq.indiv.psbutton.widget, row = 6, 
                column = 3, sticky = "e")
            tkgrid(plot.freq.indiv.pslabel.widget, row = 6, column = 4, 
                sticky = "w")
            tkgrid(loc.gridlabel.widget, row = 5, column = 1, 
                sticky = "w")
            tkgrid(loc.grid.widget, row = 5, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(loc.indivlabel.widget, row = 7, column = 1, 
                sticky = "w")
            tkgrid(loc.indiv.widget, row = 7, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(plot.genlabel.widget, row = 8, column = 1, 
                sticky = "w")
            tkgrid(wplot.gen, row = 8, column = 2, columnspan = 2, 
                sticky = "w")
            tkgrid(plot.gen.psbutton.widget, row = 8, column = 3, 
                sticky = "e")
            tkgrid(plot.gen.pslabel.widget, row = 8, column = 4, 
                sticky = "w")
            labelspace <- tklabel(ttshowibd, text = " ")
            tkgrid(labelspace, row = 9, column = 1)
            nextbutton <- tkbutton(ttshowibd, image = imagedraw, 
                text = "Draw >>", command = DrawShowIBD)
            tkgrid(nextbutton, row = 10, column = 2, columnspan = 2, 
                sticky = "e")
            tkfocus(ttshowibd)
        }
    }
    plot <- function() {
        Drift <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            DrawDrift <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(PlotDrift(path.mcmc = tclvalue(outputdir), 
                    printit = TRUE, file = tclvalue(printfile)), 
                    silent = TRUE)
                }
                else {
                  err <- try(PlotDrift(path.mcmc = tclvalue(outputdir), 
                    printit = FALSE), silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotDrift(path.mcmc=\"", tclvalue(outputdir), 
                      "\",printit=TRUE,file=\"", tclvalue(printfile), 
                      "\")", sep = ""), "[FAILED] ")
                  }
                  else {
                    Log(paste("PlotDrift(path.mcmc=\"", tclvalue(outputdir), 
                      "\",printit=FALSE", ")", sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotDrift(path.mcmc=\"", tclvalue(outputdir), 
                      "\",printit=TRUE,file=\"", tclvalue(printfile), 
                      "\")", sep = ""), "[SUCCESS] ")
                  }
                  else {
                    Log(paste("PlotDrift(path.mcmc=\"", tclvalue(outputdir), 
                      "\",printit=FALSE", ")", sep = ""), "[SUCCESS] ")
                  }
                }
            }
            ttdrift <- tktoplevel()
            tkwm.title(ttdrift, "Drift factors")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttdrift)
            }
            alellelabel.widget <- tklabel(ttdrift, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttdrift, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttdrift, text = "Save File", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttdrift, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            labelspace <- tklabel(ttdrift, text = " ")
            tkgrid(labelspace, row = 2, column = 1)
            nextbutton <- tkbutton(ttdrift, image = imagedraw, 
                text = "Draw >>", command = DrawDrift)
            tkgrid(nextbutton, row = 3, column = 2, sticky = "e")
            tkfocus(ttdrift)
        }
        Freq <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            ipop <- tclVar(1)
            iloc <- tclVar(1)
            iall <- tclVar(1)
            Drawfreq <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(PlotFreq(genotypes = globalgenotypes, 
                    path.mcmc = tclvalue(outputdir), ipop = as.numeric(tclvalue(ipop)), 
                    iloc = as.numeric(tclvalue(iloc)), iall = as.numeric(tclvalue(iall)), 
                    printit = TRUE, path = tclvalue(printfile)), 
                    silent = TRUE)
                }
                else {
                  err <- try(PlotFreq(genotypes = globalgenotypes, 
                    path.mcmc = tclvalue(outputdir), ipop = as.numeric(tclvalue(ipop)), 
                    iloc = as.numeric(tclvalue(iloc)), iall = as.numeric(tclvalue(iall)), 
                    printit = FALSE), silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotFreq(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",ipop=", 
                      as.numeric(tclvalue(ipop)), ",iloc=", as.numeric(tclvalue(iloc)), 
                      ",iall=", as.numeric(tclvalue(iall)), ",printit=TRUE,path=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  else {
                    Log(paste("PlotFreq(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",ipop=", 
                      as.numeric(tclvalue(ipop)), ",iloc=", as.numeric(tclvalue(iloc)), 
                      ",iall=", as.numeric(tclvalue(iall)), ",printit=FALSE", 
                      ")", sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotFreq(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",ipop=", 
                      as.numeric(tclvalue(ipop)), ",iloc=", as.numeric(tclvalue(iloc)), 
                      ",iall=", as.numeric(tclvalue(iall)), ",printit=TRUE,path=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                  else {
                    Log(paste("PlotFreq(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",ipop=", 
                      as.numeric(tclvalue(ipop)), ",iloc=", as.numeric(tclvalue(iloc)), 
                      ",iall=", as.numeric(tclvalue(iall)), ",printit=FALSE)", 
                      sep = ""), "[SUCCESS] ")
                  }
                }
            }
            ttfreq <- tktoplevel()
            tkwm.title(ttfreq, "Frequencies in population")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkchooseDirectory(title = "Please choose an output directory"))
                tcl("regsub", "-all", "\\\\", tclvalue(printfile), 
                  "/", printfile)
                tcl("append", printfile, "/")
                tkfocus(ttfreq)
            }
            alellelabel.widget <- tklabel(ttfreq, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttfreq, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttfreq, text = "Save Directory", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttfreq, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            ipop.widget <- tkentry(ttfreq, width = "20", textvariable = ipop)
            ipoplabel.widget <- tklabel(ttfreq, text = "Index of populations:")
            tkgrid(ipoplabel.widget, row = 3, column = 1, sticky = "w")
            tkgrid(ipop.widget, row = 3, column = 2, sticky = "w")
            iloc.widget <- tkentry(ttfreq, width = "20", textvariable = iloc)
            iloclabel.widget <- tklabel(ttfreq, text = "Index of locus:")
            tkgrid(iloclabel.widget, row = 4, column = 1, sticky = "w")
            tkgrid(iloc.widget, row = 4, column = 2, sticky = "w")
            iall.widget <- tkentry(ttfreq, width = "20", textvariable = iall)
            ialllabel.widget <- tklabel(ttfreq, text = "Index of allele:")
            tkgrid(ialllabel.widget, row = 5, column = 1, sticky = "w")
            tkgrid(iall.widget, row = 5, column = 2, sticky = "w")
            labelspace <- tklabel(ttfreq, text = " ")
            tkgrid(labelspace, row = 6, column = 1)
            nextbutton <- tkbutton(ttfreq, image = imagedraw, 
                text = "Draw >>", command = Drawfreq)
            tkgrid(nextbutton, row = 7, column = 2, sticky = "e")
            tkfocus(ttfreq)
        }
        FreqA <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            iloc <- tclVar(1)
            iall <- tclVar(1)
            Drawfreqa <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(PlotFreqA(genotypes = globalgenotypes, 
                    path.mcmc = tclvalue(outputdir), iloc = as.numeric(tclvalue(iloc)), 
                    iall = as.numeric(tclvalue(iall)), printit = TRUE, 
                    path = tclvalue(printfile)), silent = TRUE)
                }
                else {
                  err <- try(PlotFreqA(genotypes = globalgenotypes, 
                    path.mcmc = tclvalue(outputdir), iloc = as.numeric(tclvalue(iloc)), 
                    iall = as.numeric(tclvalue(iall)), printit = FALSE), 
                    silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotFreqA(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",iloc=", 
                      as.numeric(tclvalue(iloc)), ",iall=", as.numeric(tclvalue(iall)), 
                      ",printit=TRUE,path=\"", tclvalue(printfile), 
                      "\")", sep = ""), "[FAILED] ")
                  }
                  else {
                    Log(paste("PlotFreqA(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",iloc=", 
                      as.numeric(tclvalue(iloc)), ",iall=", as.numeric(tclvalue(iall)), 
                      ",printit=FALSE)", sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotFreqA(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",iloc=", 
                      as.numeric(tclvalue(iloc)), ",iall=", as.numeric(tclvalue(iall)), 
                      ",printit=TRUE,path=\"", tclvalue(printfile), 
                      "\")", sep = ""), "[SUCCESS] ")
                  }
                  else {
                    Log(paste("PlotFreqA(genotypes=", matrix2str(globalgenotypes), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",iloc=", 
                      as.numeric(tclvalue(iloc)), ",iall=", as.numeric(tclvalue(iall)), 
                      ",printit=FALSE)", sep = ""), "[SUCCESS] ")
                  }
                }
            }
            ttfreqa <- tktoplevel()
            tkwm.title(ttfreqa, "Frequencies in ancestral population")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkchooseDirectory(title = "Please choose an output directory"))
                tcl("regsub", "-all", "\\\\", tclvalue(printfile), 
                  "/", printfile)
                tcl("append", printfile, "/")
                tkfocus(ttfreqa)
            }
            alellelabel.widget <- tklabel(ttfreqa, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttfreqa, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttfreqa, text = "Save Directory", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttfreqa, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            iloc.widget <- tkentry(ttfreqa, width = "20", textvariable = iloc)
            iloclabel.widget <- tklabel(ttfreqa, text = "Index of locus:")
            tkgrid(iloclabel.widget, row = 3, column = 1, sticky = "w")
            tkgrid(iloc.widget, row = 3, column = 2, sticky = "w")
            iall.widget <- tkentry(ttfreqa, width = "20", textvariable = iall)
            ialllabel.widget <- tklabel(ttfreqa, text = "Index of allele:")
            tkgrid(ialllabel.widget, row = 4, column = 1, sticky = "w")
            tkgrid(iall.widget, row = 4, column = 2, sticky = "w")
            labelspace <- tklabel(ttfreqa, text = " ")
            tkgrid(labelspace, row = 5, column = 1)
            nextbutton <- tkbutton(ttfreqa, image = imagedraw, 
                text = "Draw >>", command = Drawfreqa)
            tkgrid(nextbutton, row = 6, column = 2, sticky = "e")
            tkfocus(ttfreqa)
        }
        Tessellation <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            Drawtessellation <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(PlotTessellation(coordinates = globalcoordinates, 
                    path.mcmc = tclvalue(outputdir), printit = TRUE, 
                    path = tclvalue(printfile)), silent = TRUE)
                }
                else {
                  err <- try(PlotTessellation(coordinates = globalcoordinates, 
                    path.mcmc = tclvalue(outputdir), printit = FALSE), 
                    silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotTessellation(coordinates=", 
                      matrix2str(globalcoordinates), ",path.mcmc=\"", 
                      tclvalue(outputdir), "\",printit=TRUE,path=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  else {
                    Log(paste("PlotTessellation(coordinates=", 
                      matrix2str(globalcoordinates), ",path.mcmc=\"", 
                      tclvalue(outputdir), "\",printit=FALSE)", 
                      sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PlotTessellation(coordinates=", 
                      matrix2str(globalcoordinates), ",path.mcmc=\"", 
                      tclvalue(outputdir), "\",printit=TRUE,path=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                  else {
                    Log(paste("PlotTessellation(coordinates=", 
                      matrix2str(globalcoordinates), ",path.mcmc=\"", 
                      tclvalue(outputdir), "\",printit=FALSE)", 
                      sep = ""), "[SUCCESS] ")
                  }
                }
            }
            tttessellation <- tktoplevel()
            tkwm.title(tttessellation, "Tessellation")
            printit <- tclVar(0)
            printfile <- tclVar("")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkchooseDirectory(title = "Please choose an output directory"))
                tcl("regsub", "-all", "\\\\", tclvalue(printfile), 
                  "/", printfile)
                tcl("append", printfile, "/")
                tkfocus(tttessellation)
            }
            alellelabel.widget <- tklabel(tttessellation, text = "Save to file?")
            cb2.widget <- tkcheckbutton(tttessellation, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(tttessellation, text = "Save Directory", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(tttessellation, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            labelspace <- tklabel(tttessellation, text = " ")
            tkgrid(labelspace, row = 2, column = 1)
            nextbutton <- tkbutton(tttessellation, image = imagedraw, 
                text = "Draw >>", command = Drawtessellation)
            tkgrid(nextbutton, row = 3, column = 2, sticky = "e")
            tkfocus(tttessellation)
        }
        Npop <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            Drawnpop <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(Plotnpop(path.mcmc = tclvalue(outputdir), 
                    burnin = as.numeric(tclvalue(burnin)), printit = TRUE, 
                    file = tclvalue(printfile)), silent = TRUE)
                }
                else {
                  err <- try(Plotnpop(path.mcmc = tclvalue(outputdir), 
                    burnin = as.numeric(tclvalue(burnin)), printit = FALSE), 
                    silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("Plotnpop(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  else {
                    Log(paste("Plotnpop(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=FALSE)", 
                      sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("Plotnpop(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                  else {
                    Log(paste("Plotnpop(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=FALSE)", 
                      sep = ""), "[SUCCESS] ")
                  }
                }
            }
            ttnpop <- tktoplevel()
            tkwm.title(ttnpop, "Number of populations")
            printit <- tclVar(0)
            printfile <- tclVar("")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttnpop)
            }
            alellelabel.widget <- tklabel(ttnpop, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttnpop, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttnpop, text = "Save File", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttnpop, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            burns <- tkentry(ttnpop, width = "20", textvariable = burnin)
            burns.widget <- tklabel(ttnpop, text = "Burnin:")
            tkgrid(burns, row = 2, column = 2, sticky = "w")
            tkgrid(burns.widget, row = 2, column = 1, sticky = "w")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            labelspace <- tklabel(ttnpop, text = " ")
            tkgrid(labelspace, row = 3, column = 1)
            nextbutton <- tkbutton(ttnpop, image = imagedraw, 
                text = "Draw >>", command = Drawnpop)
            tkgrid(nextbutton, row = 4, column = 2, sticky = "e")
            tkfocus(ttnpop)
        }
        Dpost <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            Drawdpost <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                file <- try(scan(paste(tclvalue(outputdir), "log.posterior.density.txt", 
                  sep = "")), silent = TRUE)
                if (class(file) == "try-error") {
                  tkmessageBox(message = "File hasn't been created or bad output path", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(burnin) != 0) 
                    lpp <- file[-(1:as.numeric(tclvalue(burnin)))]
                  else lpp <- file
                  mean.lpp <- mean(lpp)
                  if (tclvalue(printit) == 1) {
                    postscript(tclvalue(printfile))
                    plot.default(x = as.vector(lpp), type = "l", 
                      ylab = "Log posterior density")
                    title(main = paste("Posterior density of model (values in log)\nMean=", 
                      mean.lpp))
                    dev.off()
                  }
                  else {
                    plot.default(x = as.vector(lpp), type = "l", 
                      ylab = "Log posterior density")
                    title(main = paste("Posterior density of model (values in log),burnin=", 
                      tclvalue(burnin), "\nMean=", mean.lpp))
                  }
                  tkdestroy(tttry)
                  print("Done.")
                  if (class(file) == "try-error") {
                    if (tclvalue(printit) == 1) {
                      Log(paste("postscript(\"", tclvalue(printfile), 
                        "\")", sep = ""))
                      Log(paste("plot.default(x=", as.vector(file), 
                        ",type=\"l\",ylab=\"Log posterior density\")", 
                        sep = ""), "[FAILED] ")
                      Log(paste("title(main=paste(\"Posterior density of model (values in log)\\nMean=\"", 
                        mean.lpp, "\"))", sep = ""), "[FAILED] ")
                      Log("dev.off()")
                    }
                    else {
                      Log(paste("plot.default(x=", as.vector(file), 
                        ",type=\"l\",ylab=\"Log posterior density\")", 
                        sep = ""), "[FAILED] ")
                      Log(paste("title(main=paste(\"Posterior density of model (values in log)\\nMean=\"", 
                        mean.lpp, "\"))", sep = ""), "[FAILED] ")
                    }
                    tkmessageBox(message = file, icon = "error", 
                      type = "ok", parent = tt)
                  }
                  else {
                    if (tclvalue(printit) == 1) {
                      Log(paste("plot.default(x=", as.vector(file), 
                        ",type=\"l\",ylab=\"Log posterior density\")", 
                        sep = ""), "[SUCCESS] ")
                      Log(paste("title(main=paste(\"mean log posterior density =\",", 
                        mean.lpp, "))", sep = ""), "[SUCCESS] ")
                    }
                    else {
                      Log(paste("plot.default(x=", as.vector(file), 
                        ",type=\"l\",ylab=\"Log posterior density\")", 
                        sep = ""), "[SUCCESS] ")
                      Log(paste("title(main=paste(\"mean log posterior density =\",", 
                        mean.lpp, "))", sep = ""), "[SUCCESS] ")
                    }
                  }
                }
            }
            ttdpost <- tktoplevel()
            tkwm.title(ttdpost, "Model global density values over the MCMC")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttdpost)
            }
            alellelabel.widget <- tklabel(ttdpost, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttdpost, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttdpost, text = "Save File", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttdpost, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            burns <- tkentry(ttdpost, width = "20", textvariable = burnin)
            burns.widget <- tklabel(ttdpost, text = "Burnin:")
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            tkgrid(burns.widget, row = 2, column = 1, sticky = "w")
            tkgrid(burns, row = 2, column = 2, sticky = "w")
            labelspace <- tklabel(ttdpost, text = " ")
            tkgrid(labelspace, row = 3, column = 1)
            nextbutton <- tkbutton(ttdpost, image = imagedraw, 
                text = "Draw >>", command = Drawdpost)
            tkgrid(nextbutton, row = 4, column = 2, sticky = "e")
            tkfocus(ttdpost)
        }
        Ntile <- function() {
            printit <- tclVar(0)
            printfile <- tclVar("")
            Drawntile <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(Plotntile(path.mcmc = tclvalue(outputdir), 
                    burnin = as.numeric(tclvalue(burnin)), printit = TRUE, 
                    file = tclvalue(printfile)), silent = TRUE)
                }
                else {
                  err <- try(Plotntile(path.mcmc = tclvalue(outputdir), 
                    burnin = as.numeric(tclvalue(burnin)), printit = FALSE), 
                    silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("Plotntile(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  else {
                    Log(paste("Plotntile(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE)", 
                      sep = ""), "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("Plotntile(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                  else {
                    Log(paste("Plotntile(path.mcmc=\"", tclvalue(outputdir), 
                      "\",burnin=", tclvalue(burnin), ",printit=TRUE)", 
                      sep = ""), "[SUCCESS] ")
                  }
                }
            }
            ttntile <- tktoplevel()
            tkwm.title(ttntile, "Number of tiles")
            printit <- tclVar(0)
            printfile <- tclVar("")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttntile)
            }
            alellelabel.widget <- tklabel(ttntile, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttntile, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttntile, text = "Save File", 
                command = getprintfile, width = 15)
            filelabel.widget <- tklabel(ttntile, textvariable = printfile, 
                width = 50)
            tkdeselect(cb2.widget)
            tkconfigure(printbutton.widget, state = "disable")
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            labelspace <- tklabel(ttntile, text = " ")
            tkgrid(labelspace, row = 2, column = 1)
            nextbutton <- tkbutton(ttntile, image = imagedraw, 
                text = "Draw >>", command = Drawntile)
            tkgrid(nextbutton, row = 3, column = 2, sticky = "e")
            tkfocus(ttntile)
        }
        PosteriorM <- function() {
            ttposm <- tktoplevel()
            tkwm.title(ttposm, "Map of populations")
            write <- tclVar("FALSE")
            plotit <- tclVar("TRUE")
            maintitle <- tclVar("")
            Drawposm <- function() {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (tclvalue(printit) == 1) {
                  err <- try(PosteriorMode(coordinates = globalcoordinates, 
                    path.mcmc = tclvalue(outputdir), plotit = as.logical(tclvalue(plotit)), 
                    printit = TRUE, file = tclvalue(printfile), 
                    main.title = tclvalue(maintitle)), silent = TRUE)
                }
                else {
                  err <- try(PosteriorMode(coordinates = globalcoordinates, 
                    path.mcmc = tclvalue(outputdir), plotit = as.logical(tclvalue(plotit)), 
                    printit = FALSE, file = "", main.title = as.character(tclvalue(maintitle))), 
                    silent = TRUE)
                }
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PosteriorMode(coordinates=", matrix2str(globalcoordinates), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",plotit=", 
                      as.logical(tclvalue(plotit)), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\",main.title=\"", 
                      tclvalue(maintitle), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  else {
                    Log(paste("PosteriorMode(coordinates=", matrix2str(globalcoordinates), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",plotit=", 
                      as.logical(tclvalue(plotit)), ",printit=FALSE,file=\"", 
                      tclvalue(printfile), "\",main.title=\"", 
                      tclvalue(maintitle), "\")", sep = ""), 
                      "[FAILED] ")
                  }
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  if (tclvalue(printit) == 1) {
                    Log(paste("PosteriorMode(coordinates=", matrix2str(globalcoordinates), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",plotit=", 
                      as.logical(tclvalue(plotit)), ",printit=TRUE,file=\"", 
                      tclvalue(printfile), "\",main.title=\"", 
                      tclvalue(maintitle), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                  else {
                    Log(paste("PosteriorMode(coordinates=", matrix2str(globalcoordinates), 
                      ",path.mcmc=\"", tclvalue(outputdir), "\",plotit=", 
                      as.logical(tclvalue(plotit)), ",printit=FALSE,file=\"", 
                      tclvalue(printfile), "\",main.title=\"", 
                      tclvalue(maintitle), "\")", sep = ""), 
                      "[SUCCESS] ")
                  }
                }
            }
            printit <- tclVar(0)
            printfile <- tclVar("")
            setprint <- function() {
                if (tclvalue(printit) == 1) {
                  tkconfigure(printbutton.widget, state = "normal")
                }
                else {
                  tkconfigure(printbutton.widget, state = "disable")
                  tclvalue(printfile) <- ""
                }
            }
            getprintfile <- function() {
                tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes = "{{.ps} *.ps}"))
                tkfocus(ttposm)
            }
            alellelabel.widget <- tklabel(ttposm, text = "Save to file?")
            cb2.widget <- tkcheckbutton(ttposm, command = setprint, 
                variable = printit, onvalue = 1, offvalue = 0)
            printbutton.widget <- tkbutton(ttposm, text = "Save File", 
                command = getprintfile, width = 15, state = "disable")
            filelabel.widget <- tklabel(ttposm, textvariable = printfile, 
                width = 50)
            maintitle.widget <- tkentry(ttposm, width = "20", 
                textvariable = maintitle)
            maintitlelabel.widget <- tklabel(ttposm, text = "Graph title:")
            labelspace <- tklabel(ttposm, text = " ")
            tkgrid(labelspace, row = 6, column = 1)
            tkgrid(alellelabel.widget, row = 1, column = 1, sticky = "w")
            tkgrid(cb2.widget, row = 1, column = 2, sticky = "w")
            tkgrid(printbutton.widget, row = 1, column = 3, sticky = "w")
            tkgrid(filelabel.widget, row = 1, column = 4, sticky = "w")
            tkgrid(maintitlelabel.widget, row = 2, column = 1, 
                sticky = "w")
            tkgrid(maintitle.widget, row = 2, column = 2, sticky = "w")
            nextbutton <- tkbutton(ttposm, image = imagedraw, 
                text = "Draw >>", command = Drawposm)
            tkgrid(nextbutton, row = 7, column = 2, sticky = "e")
            tkfocus(ttposm)
        }
        if (falush == 0) {
            buttondrift <- tkbutton(ttplot, width = 30, text = "Drift factors", 
                command = Drift, state = "disabled")
            buttonfreqA <- tkbutton(ttplot, width = 30, text = "Frequencies in ancestral population", 
                command = FreqA, state = "disabled")
        }
        else {
            buttondrift <- tkbutton(ttplot, width = 30, text = "Drift factors", 
                command = Drift)
            buttonfreqA <- tkbutton(ttplot, width = 30, text = "Frequencies in ancestral population", 
                command = FreqA)
        }
        buttonfreq <- tkbutton(ttplot, width = 30, text = "Frequencies in population", 
            command = Freq)
        buttontessellation <- tkbutton(ttplot, width = 30, text = "Map of proba. of pop. membership", 
            command = Tessellation)
        buttonnpop <- tkbutton(ttplot, width = 30, text = "Number of populations", 
            command = Npop)
        buttonntile <- tkbutton(ttplot, width = 30, text = "Number of tiles", 
            command = Ntile)
        buttonposm <- tkbutton(ttplot, width = 30, text = "Map of population membership", 
            command = PosteriorM)
        buttondpost <- tkbutton(ttplot, width = 30, text = "Posterior density of model", 
            command = Dpost)
        labelfigures <- tklabel(ttplot, text = "-Graphics-", 
            font = "*-Times-bold-i-normal--20-*", foreground = "blue")
        labeltables <- tklabel(ttplot, text = "-Tables-", font = "*-Times-bold-i-normal--20-*", 
            foreground = "blue")
        labelspace1 <- tklabel(ttplot, text = " ")
        labelspace2 <- tklabel(ttplot, text = " ")
        labelinfo <- tklabel(ttplot, text = paste(tclvalue(outputdir), 
            sep = ""))
        proba.pop.membership <- tklabel(ttplot, text = "Posterior probability of population membership for each pixel")
        tkbind(proba.pop.membership, "<Any-Enter>", function() {
            tkconfigure(proba.pop.membership, foreground = "blue")
        })
        tkbind(proba.pop.membership, "<Any-Leave>", function() {
            tkconfigure(proba.pop.membership, foreground = "black")
        })
        tkbind(proba.pop.membership, "<Button-1>", function() {
            Showtextproba("proba.pop.membership.txt")
        })
        proba.pop.membership.ind <- tklabel(ttplot, text = "Posterior probability of population membership for each individual")
        tkbind(proba.pop.membership.ind, "<Any-Enter>", function() {
            tkconfigure(proba.pop.membership.ind, foreground = "blue")
        })
        tkbind(proba.pop.membership.ind, "<Any-Leave>", function() {
            tkconfigure(proba.pop.membership.ind, foreground = "black")
        })
        tkbind(proba.pop.membership.ind, "<Button-1>", function() {
            Showtextindiv("proba.pop.membership.indiv.txt")
        })
        modal.pop.ind <- tklabel(ttplot, text = "Label of modal population for each individual")
        tkbind(modal.pop.ind, "<Any-Enter>", function() {
            tkconfigure(modal.pop.ind, foreground = "blue")
        })
        tkbind(modal.pop.ind, "<Any-Leave>", function() {
            tkconfigure(modal.pop.ind, foreground = "black")
        })
        tkbind(modal.pop.ind, "<Button-1>", function() {
            Showtextmodal("modal.pop.indiv.txt")
        })
        model.global.density <- tklabel(ttplot, text = "Posterior density of model")
        tkbind(model.global.density, "<Any-Enter>", function() {
            tkconfigure(model.global.density, foreground = "blue")
        })
        tkbind(model.global.density, "<Any-Leave>", function() {
            tkconfigure(model.global.density, foreground = "black")
        })
        tkbind(model.global.density, "<Button-1>", function() {
            Showtextposterior("log.posterior.density.txt")
        })
        labelfstat <- tklabel(ttplot, text = "-F statistics-", 
            font = "*-Times-bold-i-normal--20-*", foreground = "blue")
        buttonfstat <- tkbutton(ttplot, width = 20, text = "Fst and Fis", 
            font = "*-Helvetica-bold-i-*", command = Gfstat)
        tkgrid(labelfigures, row = 1, column = 1, sticky = "w")
        tkgrid(buttondrift, row = 5, column = 2, sticky = "w")
        tkgrid(buttonfreq, row = 3, column = 2, sticky = "w")
        tkgrid(buttonfreqA, row = 4, column = 2, sticky = "w")
        tkgrid(buttontessellation, row = 4, column = 1, sticky = "w")
        tkgrid(buttonnpop, row = 2, column = 1, sticky = "w")
        tkgrid(buttonntile, row = 2, column = 2, sticky = "w")
        tkgrid(buttonposm, row = 3, column = 1, sticky = "w")
        tkgrid(buttondpost, row = 5, column = 1, sticky = "w")
        tkgrid(labelspace1, row = 6, column = 1, sticky = "w")
        tkgrid(labeltables, row = 7, column = 1, sticky = "w")
        tkgrid(proba.pop.membership, row = 8, column = 1, columnspan = 2, 
            sticky = "w")
        tkgrid(proba.pop.membership.ind, row = 9, column = 1, 
            columnspan = 2, sticky = "w")
        tkgrid(modal.pop.ind, row = 10, column = 1, columnspan = 2, 
            sticky = "w")
        tkgrid(model.global.density, row = 11, column = 1, columnspan = 2, 
            sticky = "w")
        tkgrid(labelspace2, row = 12, column = 1, sticky = "w")
        tkgrid(labelfstat, row = 13, column = 1, sticky = "w")
        tkgrid(buttonfstat, row = 14, column = 1, sticky = "w")
    }
    plot2 <- function() {
        pfstat <- function() {
            if (length(idb.dataset) == 1) {
                tkmessageBox(message = "First simulate some data", 
                  icon = "error", type = "ok", parent = tt)
            }
            else {
                tttext <- tktoplevel(parent = .TkRoot)
                tkwm.title(tttext, "F Statistics")
                yscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  command = function(...) tkyview(txt, ...))
                xscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  orient = "horizontal", command = function(...) tkxview(txt, 
                    ...))
                txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) tkset(yscr, 
                    ...), xscrollcommand = function(...) tkset(xscr, 
                    ...))
                tkinsert(txt, "end", "Pairwise Fst\n\n\n")
                for (i in 1:nrow(idb.dataset$Fst)) {
                  for (j in 1:ncol(idb.dataset$Fst)) {
                    tkinsert(txt, "end", idb.dataset$Fst[i, j])
                    tkinsert(txt, "end", "\t")
                  }
                  tkinsert(txt, "end", "\n")
                }
                tkinsert(txt, "end", "\n\n\nFis\n\n\n")
                for (j in 1:length(idb.dataset$Fis)) {
                  tkinsert(txt, "end", idb.dataset$Fis[j])
                  tkinsert(txt, "end", "\t")
                }
                tkinsert(txt, "end", "\n")
                tkgrid(txt, row = 1, column = 1)
                tkgrid(yscr, row = 1, column = 2, sticky = "ns")
                tkgrid(xscr, row = 2, column = 1, sticky = "we")
            }
        }
        pdsigma <- function() {
            if (length(idb.dataset) == 1) {
                tkmessageBox(message = "First simulate some data", 
                  icon = "error", type = "ok", parent = tt)
            }
            else {
                tttext <- tktoplevel(parent = .TkRoot)
                tkwm.title(tttext, "Dsigma2")
                yscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  command = function(...) tkyview(txt, ...))
                xscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  orient = "horizontal", command = function(...) tkxview(txt, 
                    ...))
                txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) tkset(yscr, 
                    ...), xscrollcommand = function(...) tkset(xscr, 
                    ...))
                tkinsert(txt, "end", "Dsigma2\n\n\n")
                for (j in 1:length(idb.dataset$Dsigma2)) {
                  tkinsert(txt, "end", idb.dataset$Dsigma2[j])
                  tkinsert(txt, "end", "\t")
                }
                tkinsert(txt, "end", "\n")
                tkgrid(txt, row = 1, column = 1)
                tkgrid(yscr, row = 1, column = 2, sticky = "ns")
                tkgrid(xscr, row = 2, column = 1, sticky = "we")
            }
        }
        pdiff <- function() {
            if (length(idb.dataset) == 1) {
                tkmessageBox(message = "First simulate some data", 
                  icon = "error", type = "ok", parent = tt)
            }
            else {
                print(idb.dataset$diff.W)
                print(idb.dataset$diff.B)
                tttext <- tktoplevel(parent = .TkRoot)
                tkwm.title(tttext, "Variability of allele frequency")
                yscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  command = function(...) tkyview(txt, ...))
                xscr <- tkscrollbar(tttext, repeatinterval = 5, 
                  orient = "horizontal", command = function(...) tkxview(txt, 
                    ...))
                txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) tkset(yscr, 
                    ...), xscrollcommand = function(...) tkset(xscr, 
                    ...))
                tkinsert(txt, "end", "Differentiation within populations\n\n\n")
                for (j in 1:length(idb.dataset$diff.W)) {
                  tkinsert(txt, "end", idb.dataset$diff.W[j])
                  tkinsert(txt, "end", "\t")
                }
                tkinsert(txt, "end", "\n\n\nDifferentiation around barrier between populations\n\n\n")
                for (i in 1:nrow(idb.dataset$diff.B)) {
                  for (j in 1:ncol(idb.dataset$diff.B)) {
                    tkinsert(txt, "end", idb.dataset$diff.B[i, 
                      j])
                    tkinsert(txt, "end", "\t")
                  }
                  tkinsert(txt, "end", "\n")
                }
                tkinsert(txt, "end", "\n")
                tkgrid(txt, row = 1, column = 1)
                tkgrid(yscr, row = 1, column = 2, sticky = "ns")
                tkgrid(xscr, row = 2, column = 1, sticky = "we")
            }
        }
        buttonshowibd <- tkbutton(ttplot2, width = 30, text = "Show simulated data", 
            command = GraficalIBD)
        labelfigures <- tklabel(ttplot2, text = "-Graphics-", 
            font = "*-Times-bold-i-normal--20-*", foreground = "blue")
        labeltables <- tklabel(ttplot2, text = "-Tables-", font = "*-Times-bold-i-normal--20-*", 
            foreground = "blue")
        labelspace1 <- tklabel(ttplot2, text = " ")
        labelfstat <- tklabel(ttplot2, text = "-F statistic-", 
            font = "*-Times-bold-i-normal--20-*", foreground = "blue")
        buttonfstat <- tkbutton(ttplot2, width = 30, text = "Fst and Fis", 
            command = pfstat)
        labelspace2 <- tklabel(ttplot2, text = " ")
        labeldsigma <- tklabel(ttplot2, text = "-Dsigma2-", font = "*-Times-bold-i-normal--20-*", 
            foreground = "blue")
        buttondsigma <- tkbutton(ttplot2, width = 30, text = "Dsigma2", 
            command = pdsigma)
        labelspace3 <- tklabel(ttplot2, text = " ")
        labeldiff <- tklabel(ttplot2, text = "-Differentiation-", 
            font = "*-Times-bold-i-normal--20-*", foreground = "blue")
        buttondiff <- tkbutton(ttplot2, width = 30, text = "Within and between populations", 
            command = pdiff)
        labelspace4 <- tklabel(ttplot2, text = "                                      ")
        tkgrid(labelfigures, row = 1, column = 1, sticky = "w")
        tkgrid(labelspace4, row = 1, column = 2, sticky = "w")
        tkgrid(buttonshowibd, row = 2, column = 1, sticky = "w")
        tkgrid(labelspace1, row = 3, column = 1, sticky = "w")
        tkgrid(labelfstat, row = 4, column = 1, sticky = "w")
        tkgrid(buttonfstat, row = 5, column = 1, sticky = "w")
        tkgrid(labelspace2, row = 6, column = 1, sticky = "w")
        tkgrid(labeldiff, row = 10, column = 1, sticky = "w")
        tkgrid(buttondiff, row = 11, column = 1, sticky = "w")
    }
    Gfstat <- function() {
        Runfstat <- function() {
            ShowtextFis <- function() {
                file <- try(read.table(paste(tclvalue(outputdir), 
                  "Fis.txt", sep = "")), silent = TRUE)
                if (class(file) == "try-error") {
                  tkmessageBox(message = "File hasn't been created or bad output path", 
                    type = "ok", parent = tt)
                }
                else {
                  tttext <- tktoplevel(parent = .TkRoot)
                  tkwm.title(tttext, "FIS")
                  posx <- tclVar("")
                  posy <- tclVar("")
                  yscr <- tkscrollbar(tttext, repeatinterval = 5, 
                    command = function(...) {
                      tkyview(txt, ...)
                      tkyview(left, ...)
                    })
                  xscr <- tkscrollbar(tttext, repeatinterval = 5, 
                    orient = "horizontal", command = function(...) {
                      tkxview(txt, ...)
                    })
                  txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                    wrap = "none", yscrollcommand = function(...) {
                      tkset(yscr, ...)
                      tkyview.moveto(left, as.double(...))
                    }, xscrollcommand = function(...) {
                      tkset(xscr, ...)
                    })
                  left <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                    width = 5, wrap = "none", yscrollcommand = function(...) {
                      tkset(yscr, ...)
                      tkyview.moveto(txt, as.double(...))
                    })
                  col <- NCOL(file)
                  row <- NROW(file)
                  for (i in 1:row) {
                    for (j in 1:col) {
                      if (as.double(file[i, j]) == -999) 
                        tkinsert(txt, "end", sprintf("%.3f", 
                          as.double(file[i, j])))
                      else tkinsert(txt, "end", sprintf("%.5f", 
                        as.double(file[i, j])))
                      tkinsert(txt, "end", "\t\t")
                    }
                    tkinsert(txt, "end", "\n")
                  }
                  for (i in 1:row) {
                    tkinsert(left, "end", paste("pop", as.character(i), 
                      sep = ""))
                    tkinsert(left, "end", "\n")
                  }
                  tkdestroy(tttry)
                  tkgrid(txt, row = 2, column = 2)
                  tkgrid(left, row = 2, column = 1)
                  tkgrid(yscr, row = 2, column = 3, sticky = "ns")
                  tkgrid(xscr, row = 3, column = 2, sticky = "we")
                }
            }
            ShowtextFst <- function() {
                file <- try(read.table(paste(tclvalue(outputdir), 
                  "Fst.txt", sep = "")), silent = TRUE)
                if (class(file) == "try-error") {
                  tkmessageBox(message = "File hasn't been created or bad output path", 
                    type = "ok", parent = tt)
                }
                else {
                  tttext <- tktoplevel(parent = .TkRoot)
                  tkwm.title(tttext, "FST")
                  posx <- tclVar("")
                  posy <- tclVar("")
                  yscr <- tkscrollbar(tttext, repeatinterval = 5, 
                    command = function(...) {
                      tkyview(txt, ...)
                      tkyview(left, ...)
                    })
                  xscr <- tkscrollbar(tttext, repeatinterval = 5, 
                    orient = "horizontal", command = function(...) {
                      tkxview(txt, ...)
                      tkxview(top, ...)
                    })
                  txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                    wrap = "none", yscrollcommand = function(...) {
                      tkset(yscr, ...)
                      tkyview.moveto(left, as.double(...))
                    }, xscrollcommand = function(...) {
                      tkset(xscr, ...)
                      tkxview.moveto(top, as.double(...))
                    })
                  top <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                    height = 1, wrap = "none", xscrollcommand = function(...) {
                      tkset(xscr, ...)
                      tkxview.moveto(txt, as.double(...))
                    })
                  left <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                    width = 5, wrap = "none", yscrollcommand = function(...) {
                      tkset(yscr, ...)
                      tkyview.moveto(txt, as.double(...))
                    })
                  col <- NCOL(file)
                  row <- NROW(file)
                  for (i in 1:row) {
                    for (j in 1:col) {
                      if (as.double(file[i, j]) == -999) 
                        tkinsert(txt, "end", sprintf("%.2f", 
                          as.double(file[i, j])))
                      else tkinsert(txt, "end", sprintf("%.5f", 
                        as.double(file[i, j])))
                      tkinsert(txt, "end", "\t\t")
                    }
                    tkinsert(txt, "end", "\n")
                  }
                  for (i in 1:row) {
                    tkinsert(left, "end", paste("pop", as.character(i), 
                      sep = ""))
                    tkinsert(left, "end", "\n")
                  }
                  for (i in 1:col) {
                    tkinsert(top, "end", paste("pop", as.character(i), 
                      sep = ""))
                    tkinsert(top, "end", "\t\t")
                  }
                  tkdestroy(tttry)
                  tkgrid(top, row = 1, column = 2)
                  tkgrid(txt, row = 2, column = 2)
                  tkgrid(left, row = 2, column = 1)
                  tkgrid(yscr, row = 2, column = 3, sticky = "ns")
                  tkgrid(xscr, row = 3, column = 2, sticky = "we")
                }
            }
            tttry <- tktoplevel(parent = .TkRoot)
            tkgrab(tttry)
            tkwm.geometry(tttry, "+200+200")
            tkwm.title(tttry, "wait")
            warn <- tklabel(tttry, image = imagepleasewait)
            tkpack(warn)
            tkfocus(tttry)
            tcl("update")
            print("Starting...")
            Sys.sleep(0.5)
            param <- as.matrix(read.table(paste(tclvalue(outputdir), 
                "parameters.txt", sep = "")))
            nploidy <- as.numeric(param[param[, 1] == "ploidy", 
                3])
            if (as.numeric(nploidy) == 1) {
                tkmessageBox(message = "FST for haploid data:\nNot implemented", 
                  icon = "error", type = "ok", parent = tt)
                tkdestroy(tttry)
                print("Done.")
            }
            else {
                err <- try(Fstat.output(coordinates = NULL, genotypes = globalgenotypes, 
                  ploidy = nploidy, burnin = NULL, path.mcmc = tclvalue(outputdir)), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  Log(paste("Fstat.output(coordinates=NULL,genotypes=", 
                    matrix2str(globalgenotypes), ",ploidy=", 
                    nploidy, ",burnin=NULL,path.mcmc=\"", tclvalue(outputdir), 
                    "\")", sep = ""), "[FAILED] ")
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  Log(paste("Fstat.output(coordinates=", matrix2str(globalcoordinates), 
                    ",genotypes=", matrix2str(globalgenotypes), 
                    ",ploidy=", nploidy, ",burnin=NULL,path.mcmc=\"", 
                    tclvalue(outputdir), "\")", sep = ""), "[SUCCESS] ")
                  if (tclvalue(sep1) == "White space") 
                    tclvalue(sep1) <- " "
                  if (tclvalue(sep2) == "White space") 
                    tclvalue(sep2) <- " "
                  write.table(err$Fis, file = paste(tclvalue(outputdir), 
                    "Fis.txt", sep = ""), row.names = FALSE, 
                    col.names = FALSE)
                  write.table(err$Fst, file = paste(tclvalue(outputdir), 
                    "Fst.txt", sep = ""), row.names = FALSE, 
                    col.names = FALSE)
                  ShowtextFis()
                  ShowtextFst()
                  if (tclvalue(sep1) == " ") 
                    tclvalue(sep1) <- "White space"
                  if (tclvalue(sep2) == " ") 
                    tclvalue(sep2) <- "White space"
                }
            }
        }
        Runfstat()
    }
    SimnonIBD <- function() {
        runnonibd <- function() {
            if (tclvalue(sep1) == "White space") 
                tclvalue(sep1) <- " "
            if (tclvalue(sep2) == "White space") 
                tclvalue(sep2) <- " "
            vect1 <- c()
            vec1 <- unlist(strsplit(tclvalue(nall), ","))
            for (i in 1:length(vec1)) vect1[i] <- as.numeric(vec1[i])
            if (tclvalue(nloc) != length(vect1)) {
                tkmessageBox(message = "Number of locus must be equal to the length of the number of alleles per locus", 
                  icon = "error", type = "ok", parent = tt)
            }
            else {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (is.null(globalcoordinates) || tclvalue(prevcoord) == 
                  0) 
                  idb.dataset <<- try(simdata(nindiv = as.numeric(tclvalue(nindiv)), 
                    coord.lim = c(as.numeric(tclvalue(absmin)), 
                      as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                      as.numeric(tclvalue(ordmax))), number.nuclei = as.numeric(tclvalue(nuclei)), 
                    allele.numbers = vect1, IBD = FALSE, npop = as.numeric(tclvalue(npop)), 
                    give.freq.grid = as.logical(tclvalue(freq.grid)), 
                    give.tess.grid = as.logical(tclvalue(tess.grid)), 
                    npix = c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                    comp.Fst = as.logical(tclvalue(comp.Fst))), 
                    silent = TRUE)
                else idb.dataset <<- try(simdata(nindiv = as.numeric(tclvalue(nindiv)), 
                  coord.indiv = globalcoordinates, coord.lim = c(as.numeric(tclvalue(absmin)), 
                    as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                    as.numeric(tclvalue(ordmax))), number.nuclei = as.numeric(tclvalue(nuclei)), 
                  allele.numbers = vect1, IBD = FALSE, npop = as.numeric(tclvalue(npop)), 
                  give.freq.grid = as.logical(tclvalue(freq.grid)), 
                  give.tess.grid = as.logical(tclvalue(tess.grid)), 
                  npix = c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                  comp.Fst = as.logical(tclvalue(comp.Fst))), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(idb.dataset) == "try-error") {
                  tkmessageBox(message = idb.dataset, icon = "error", 
                    type = "ok", parent = tt)
                  idb.dataset <<- 0
                }
                else {
                  if (is.null(globalcoordinates) || tclvalue(prevcoord) == 
                    0) 
                    Log(paste("simdata(nindiv=", as.numeric(tclvalue(nindiv)), 
                      ",coord.lim=", c(as.numeric(tclvalue(absmin)), 
                        as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                        as.numeric(tclvalue(ordmax))), ",number.nuclei=", 
                      as.numeric(tclvalue(nuclei)), ",allele.numbers=", 
                      vect1, ",IBD=FALSE,npop=", as.numeric(tclvalue(npop)), 
                      ",give.freq.grid=", as.logical(tclvalue(freq.grid)), 
                      ",give.tess.grid=", as.logical(tclvalue(tess.grid)), 
                      ",npix=", c(as.numeric(tclvalue(npixh)), 
                        as.numeric(tclvalue(npixv))), ",comp.Fst=", 
                      as.logical(tclvalue(comp.Fst)), ")", sep = ""), 
                      "[SUCCESS] ")
                  else Log(paste("simdata(nindiv=", as.numeric(tclvalue(nindiv)), 
                    "coord.indiv=", matrix2str(globalcoordinates), 
                    ",coord.lim=", c(as.numeric(tclvalue(absmin)), 
                      as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                      as.numeric(tclvalue(ordmax))), ",number.nuclei=", 
                    as.numeric(tclvalue(nuclei)), ",allele.numbers=", 
                    vect1, ",IBD=FALSE,npop=", as.numeric(tclvalue(npop)), 
                    ",give.freq.grid=", as.logical(tclvalue(freq.grid)), 
                    ",give.tess.grid=", as.logical(tclvalue(tess.grid)), 
                    ",npix=", c(as.numeric(tclvalue(npixh)), 
                      as.numeric(tclvalue(npixv))), ",comp.Fst=", 
                    as.logical(tclvalue(comp.Fst)), ")", sep = ""), 
                    "[SUCCESS] ")
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                  globalcoordinates <<- idb.dataset$coord.indiv
                  tclvalue(labelcoordtext) <- "Coordinate:    Simulated panmictic data loaded"
                  globalgenotypes <<- idb.dataset$genotypes
                  tclvalue(labelgenotext) <- "Genotype:      Simulated panmictic data loaded"
                  if (tclvalue(save) == 1) {
                    auxcoord <- tclVar()
                    tclvalue(auxcoord) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                      initialdir = tclvalue(outputdir), title = "Save coordinate file to:"))
                    auxgen <- tclVar()
                    tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                      initialdir = tclvalue(outputdir), title = "Save genotype file to:"))
                    write.table(idb.dataset$coord.indiv, file = tclvalue(auxcoord), 
                      sep = tclvalue(sep1), row.names = FALSE, 
                      col.names = FALSE)
                    write.table(idb.dataset$genotypes, file = tclvalue(auxgen), 
                      sep = tclvalue(sep2), row.names = FALSE, 
                      col.names = FALSE)
                  }
                }
            }
            if (tclvalue(sep1) == " ") 
                tclvalue(sep1) <- "White space"
            if (tclvalue(sep2) == " ") 
                tclvalue(sep2) <- "White space"
        }
        nindiv = tclVar("0")
        nindiv.widget <- tkentry(ttsimf, width = "20", textvariable = nindiv)
        nindivlabel.widget <- tklabel(ttsimf, text = "Number of individuals:")
        coordxlabel.widget <- tklabel(ttsimf, text = "Limits of geographical domain:")
        absmin <- tclVar(0)
        absmax <- tclVar(1)
        absmin.widget <- tkentry(ttsimf, width = "8", textvariable = absmin)
        absmax.widget <- tkentry(ttsimf, width = "8", textvariable = absmax)
        abslabel.widget <- tklabel(ttsimf, text = "   abs (min|max) :")
        ordmin <- tclVar(0)
        ordmax <- tclVar(1)
        ordmin.widget <- tkentry(ttsimf, width = "8", textvariable = ordmin)
        ordmax.widget <- tkentry(ttsimf, width = "8", textvariable = ordmax)
        ordlabel.widget <- tklabel(ttsimf, text = "   ord (min|max) :")
        nuclei = tclVar("0")
        nuclei.widget <- tkentry(ttsimf, width = "20", textvariable = nuclei)
        nucleilabel.widget <- tklabel(ttsimf, text = "Number of nuclei in tessellation:")
        nloc <- tclVar(0)
        nloc.widget <- tkentry(ttsimf, width = "20", textvariable = nloc)
        nloclabel.widget <- tklabel(ttsimf, text = "Number of loci:")
        nall <- tclVar()
        nall.widget <- tkentry(ttsimf, width = "20", textvariable = nall)
        nalllabel.widget <- tklabel(ttsimf, text = "Number of alleles per locus (E.g: 10,3,8,..):")
        npop <- tclVar("")
        npop.widget <- tkentry(ttsimf, width = "20", textvariable = npop)
        npoplabel.widget <- tklabel(ttsimf, text = "Number of populations:")
        freq.gridlabel.widget <- tklabel(ttsimf, text = "Return frequencies on grid:")
        freq.grid <- tclVar("FALSE")
        wfreq.grid <- .Tk.subwin(ttsimf)
        freq.gridoptionmenu.widget <- tcl("tk_optionMenu", wfreq.grid, 
            freq.grid, "FALSE", "TRUE")
        tess.gridlabel.widget <- tklabel(ttsimf, text = "Return population membership on grid:")
        tess.grid <- tclVar("FALSE")
        wtess.grid <- .Tk.subwin(ttsimf)
        tess.gridoptionmenu.widget <- tcl("tk_optionMenu", wtess.grid, 
            tess.grid, "FALSE", "TRUE")
        npixh <- tclVar(50)
        npixv <- tclVar(50)
        npixh.widget <- tkentry(ttsimf, width = "8", textvariable = npixh)
        npixv.widget <- tkentry(ttsimf, width = "8", textvariable = npixv)
        npixlabel.widget <- tklabel(ttsimf, text = "Number of pixels for representation (hor|ver):")
        comp.Fstlabel.widget <- tklabel(ttsimf, text = "Compute F statistics:")
        comp.Fst <- tclVar("FALSE")
        wcomp.Fst <- .Tk.subwin(ttsimf)
        comp.Fstoptionmenu.widget <- tcl("tk_optionMenu", wcomp.Fst, 
            comp.Fst, "FALSE", "TRUE")
        prevcoord <- tclVar(0)
        prevcoordlabel.widget <- tklabel(ttsimf, text = "Use loaded coordinate file:")
        prevcoord.widget <- tkcheckbutton(ttsimf, variable = save, 
            onvalue = 1, offvalue = 0)
        save <- tclVar(0)
        savelabel.widget <- tklabel(ttsimf, text = "Save coordinate and genotype files:")
        save.widget <- tkcheckbutton(ttsimf, variable = save, 
            onvalue = 1, offvalue = 0)
        tkgrid(nindivlabel.widget, row = 1, column = 1, sticky = "w")
        tkgrid(nindiv.widget, row = 1, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(coordxlabel.widget, row = 2, column = 1, sticky = "w")
        tkgrid(abslabel.widget, row = 3, column = 1, sticky = "w")
        tkgrid(absmin.widget, row = 3, column = 2, sticky = "w")
        tkgrid(absmax.widget, row = 3, column = 3, sticky = "w")
        tkgrid(ordlabel.widget, row = 4, column = 1, sticky = "w")
        tkgrid(ordmin.widget, row = 4, column = 2, sticky = "w")
        tkgrid(ordmax.widget, row = 4, column = 3, sticky = "w")
        tkgrid(nucleilabel.widget, row = 5, column = 1, sticky = "w")
        tkgrid(nuclei.widget, row = 5, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(nloclabel.widget, row = 6, column = 1, sticky = "w")
        tkgrid(nloc.widget, row = 6, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(nalllabel.widget, row = 7, column = 1, sticky = "w")
        tkgrid(nall.widget, row = 7, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(npoplabel.widget, row = 10, column = 1, sticky = "w")
        tkgrid(npop.widget, row = 10, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(freq.gridlabel.widget, row = 11, column = 1, sticky = "w")
        tkgrid(wfreq.grid, row = 11, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(tess.gridlabel.widget, row = 12, column = 1, sticky = "w")
        tkgrid(wtess.grid, row = 12, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(npixlabel.widget, row = 13, column = 1, sticky = "w")
        tkgrid(npixh.widget, row = 13, column = 2, sticky = "w")
        tkgrid(npixv.widget, row = 13, column = 3, sticky = "w")
        tkgrid(comp.Fstlabel.widget, row = 14, column = 1, sticky = "w")
        tkgrid(wcomp.Fst, row = 14, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(prevcoordlabel.widget, row = 15, column = 1, columnspan = 2, 
            sticky = "w")
        tkgrid(prevcoord.widget, row = 15, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(savelabel.widget, row = 19, column = 1, sticky = "w")
        tkgrid(save.widget, row = 19, column = 2, columnspan = 2, 
            sticky = "w")
        labelspace <- tklabel(ttsimf, text = " ")
        tkgrid(labelspace, row = 20, column = 1)
        nextbutton <- tkbutton(ttsimf, image = imagerun2, text = "RUN >>", 
            command = runnonibd)
        tkgrid(nextbutton, row = 21, column = 2, columnspan = 2, 
            sticky = "e")
    }
    Convert <- function() {
        ttcon <- tktoplevel()
        filename <- tclVar("")
        tkwm.title(ttcon, "Convert loaded data into genepop format")
        gltgp <- function() {
            if (tclvalue(filename) == "" | is.null(globalcoordinates) | 
                length(globalgenotypes) == 1) {
                tkmessageBox(message = "You must define filename, coordinate file and genotype file", 
                  icon = "error", type = "ok")
            }
            else {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                err <- try(gl2gp(coordinates = globalcoordinates, 
                  genotypes = globalgenotypes, file = tclvalue(filename)), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  Log(paste("gl2gp(coordinates=", matrix2str(globalcoordinates), 
                    ",genotypes=", matrix2str(globalgenotypes), 
                    ",file=", tclvalue(filename), ")", sep = ""), 
                    "[FAILED] ")
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                }
                else {
                  Log(paste("gl2gp(coordinates=", matrix2str(globalcoordinates), 
                    ",genotypes=", matrix2str(globalgenotypes), 
                    ",file=", tclvalue(filename), ")", sep = ""), 
                    "[SUCCESS] ")
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                }
            }
        }
        setgenepopfile <- function() {
            tclvalue(filename) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                title = "Choose a File"))
        }
        genepopbutton.widget <- tkbutton(ttcon, text = "Choose file name", 
            command = setgenepopfile, width = 15)
        filelabel.widget <- tklabel(ttcon, textvariable = filename, 
            width = 50)
        tkgrid(genepopbutton.widget, row = 1, column = 1, sticky = "w")
        tkgrid(filelabel.widget, row = 1, column = 2, sticky = "w")
        labelspace <- tklabel(ttcon, text = " ")
        tkgrid(labelspace, row = 2, column = 1)
        nextbutton <- tkbutton(ttcon, image = imageconvert, text = "RUN >>", 
            command = gltgp)
        tkgrid(nextbutton, row = 3, column = 2, sticky = "e")
        tkfocus(ttcon)
    }
    SimIBD <- function() {
        runibd <- function() {
            if (tclvalue(sep1) == "White space") 
                tclvalue(sep1) <- " "
            if (tclvalue(sep2) == "White space") 
                tclvalue(sep2) <- " "
            vect1 <- c()
            vec1 <- unlist(strsplit(tclvalue(nall), ","))
            for (i in 1:length(vec1)) vect1[i] <- as.numeric(vec1[i])
            if (tclvalue(nloc) != length(vect1)) {
                tkmessageBox(message = "Number of locus must be equal to the length of the number of alleles per locus", 
                  icon = "error", type = "ok", parent = tt)
            }
            else {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                if (is.null(globalcoordinates) || tclvalue(prevcoord) == 
                  0) 
                  idb.dataset <<- try(simdata(nindiv = as.numeric(tclvalue(nindiv)), 
                    coord.lim = c(as.numeric(tclvalue(absmin)), 
                      as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                      as.numeric(tclvalue(ordmax))), number.nuclei = as.numeric(tclvalue(nuclei)), 
                    allele.numbers = vect1, IBD = TRUE, beta = as.numeric(tclvalue(beta)), 
                    npop = as.numeric(tclvalue(npop)), give.freq.grid = as.logical(tclvalue(freq.grid)), 
                    give.tess.grid = as.logical(tclvalue(tess.grid)), 
                    npix = c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                    comp.Fst = as.logical(tclvalue(comp.Fst)), 
                    comp.Dsigma2 = as.logical(tclvalue(comp.Dsigma2)), 
                    comp.diff = as.logical(tclvalue(comp.height)), 
                    width = as.numeric(tclvalue(hwidth)), plot.pairs.borders = as.logical(tclvalue(plot.pairs.borders))), 
                    silent = TRUE)
                else idb.dataset <<- try(simdata(nindiv = as.numeric(tclvalue(nindiv)), 
                  coord.indiv = globalcoordinates, coord.lim = c(as.numeric(tclvalue(absmin)), 
                    as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                    as.numeric(tclvalue(ordmax))), number.nuclei = as.numeric(tclvalue(nuclei)), 
                  allele.numbers = vect1, IBD = TRUE, beta = as.numeric(tclvalue(beta)), 
                  npop = as.numeric(tclvalue(npop)), give.freq.grid = as.logical(tclvalue(freq.grid)), 
                  give.tess.grid = as.logical(tclvalue(tess.grid)), 
                  npix = c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                  comp.Fst = as.logical(tclvalue(comp.Fst)), 
                  comp.Dsigma2 = as.logical(tclvalue(comp.Dsigma2)), 
                  comp.diff = as.logical(tclvalue(comp.height)), 
                  width = as.numeric(tclvalue(hwidth)), plot.pairs.borders = as.logical(tclvalue(plot.pairs.borders))), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(idb.dataset) == "try-error") {
                  tkmessageBox(message = idb.dataset, icon = "error", 
                    type = "ok", parent = tt)
                  idb.dataset <<- 0
                }
                else {
                  if (is.null(globalcoordinates) || tclvalue(prevcoord) == 
                    0) 
                    Log(paste("simdata(nindiv=", as.numeric(tclvalue(nindiv)), 
                      ",coord.lim=", c(as.numeric(tclvalue(absmin)), 
                        as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                        as.numeric(tclvalue(ordmax))), ",number.nuclei=", 
                      as.numeric(tclvalue(nuclei)), ",allele.numbers=", 
                      vect1, ",IBD=TRUE,beta=", as.numeric(tclvalue(beta)), 
                      ",npop=", as.numeric(tclvalue(npop)), ",give.freq.grid=", 
                      as.logical(tclvalue(freq.grid)), ",give.tess.grid=", 
                      as.logical(tclvalue(tess.grid)), ",npix=", 
                      c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                      ",comp.Fst=", as.logical(tclvalue(comp.Fst)), 
                      ",comp.Dsigma2=", as.logical(tclvalue(comp.Dsigma2)), 
                      ",comp.diff=", as.logical(tclvalue(comp.height)), 
                      ",width=", as.numeric(tclvalue(hwidth)), 
                      ",plot.pairs.borders=", as.logical(tclvalue(plot.pairs.borders)), 
                      ")", sep = ""), "[SUCCESS] ")
                  else Log(paste("simdata(nindiv=", as.numeric(tclvalue(nindiv)), 
                    ",coord.indiv=", matrix2str(globalcoordinates), 
                    ",coord.lim=", c(as.numeric(tclvalue(absmin)), 
                      as.numeric(tclvalue(absmax)), as.numeric(tclvalue(ordmin)), 
                      as.numeric(tclvalue(ordmax))), ",number.nuclei=", 
                    as.numeric(tclvalue(nuclei)), ",allele.numbers=", 
                    vect1, ",IBD=TRUE,beta=", as.numeric(tclvalue(beta)), 
                    ",npop=", as.numeric(tclvalue(npop)), ",give.freq.grid=", 
                    as.logical(tclvalue(freq.grid)), ",give.tess.grid=", 
                    as.logical(tclvalue(tess.grid)), ",npix=", 
                    c(as.numeric(tclvalue(npixh)), as.numeric(tclvalue(npixv))), 
                    ",comp.Fst=", as.logical(tclvalue(comp.Fst)), 
                    ",comp.Dsigma2=", as.logical(tclvalue(comp.Dsigma2)), 
                    ",comp.diff=", as.logical(tclvalue(comp.height)), 
                    ",width=", as.numeric(tclvalue(hwidth)), 
                    ",plot.pairs.borders=", as.logical(tclvalue(plot.pairs.borders)), 
                    ")", sep = ""), "[SUCCESS] ")
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                  globalcoordinates <<- idb.dataset$coord.indiv
                  tclvalue(labelcoordtext) <- "Coordinate:    Simulated IBD data loaded"
                  globalgenotypes <<- idb.dataset$genotypes
                  tclvalue(labelgenotext) <- "Genotype:       Simulated IBD data loaded"
                  if (tclvalue(save) == 1) {
                    auxcoord <- tclVar()
                    tclvalue(auxcoord) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                      initialdir = tclvalue(outputdir), title = "Save coordinate file to:"))
                    auxgen <- tclVar()
                    tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                      initialdir = tclvalue(outputdir), title = "Save genotype file to:"))
                    write.table(idb.dataset$coord.indiv, file = tclvalue(auxcoord), 
                      sep = tclvalue(sep1), row.names = FALSE, 
                      col.names = FALSE)
                    write.table(idb.dataset$genotypes, file = tclvalue(auxgen), 
                      sep = tclvalue(sep2), row.names = FALSE, 
                      col.names = FALSE)
                  }
                }
            }
            if (tclvalue(sep1) == " ") 
                tclvalue(sep1) <- "White space"
            if (tclvalue(sep2) == " ") 
                tclvalue(sep2) <- "White space"
        }
        nindiv = tclVar("0")
        nindiv.widget <- tkentry(ttibd, width = "20", textvariable = nindiv)
        nindivlabel.widget <- tklabel(ttibd, text = "Number of individuals:")
        coordxlabel.widget <- tklabel(ttibd, text = "Limits of geographical domain:")
        absmin <- tclVar(0)
        absmax <- tclVar(1)
        absmin.widget <- tkentry(ttibd, width = "8", textvariable = absmin)
        absmax.widget <- tkentry(ttibd, width = "8", textvariable = absmax)
        abslabel.widget <- tklabel(ttibd, text = "   abs (min|max) :")
        ordmin <- tclVar(0)
        ordmax <- tclVar(1)
        ordmin.widget <- tkentry(ttibd, width = "8", textvariable = ordmin)
        ordmax.widget <- tkentry(ttibd, width = "8", textvariable = ordmax)
        ordlabel.widget <- tklabel(ttibd, text = "   ord (min|max) :")
        nuclei = tclVar("0")
        nuclei.widget <- tkentry(ttibd, width = "20", textvariable = nuclei)
        nucleilabel.widget <- tklabel(ttibd, text = "Number of nuclei in tessellation:")
        nloc <- tclVar(0)
        nloc.widget <- tkentry(ttibd, width = "20", textvariable = nloc)
        nloclabel.widget <- tklabel(ttibd, text = "Number of loci:")
        nall <- tclVar()
        nall.widget <- tkentry(ttibd, width = "20", textvariable = nall)
        nalllabel.widget <- tklabel(ttibd, text = "Number of alleles per locus (E.g: 10,3,8,..):")
        beta <- tclVar("")
        beta.widget <- tkentry(ttibd, width = "20", textvariable = beta)
        betalabel.widget <- tklabel(ttibd, text = "Spatial correlation parameter for frequencies:")
        npop <- tclVar("")
        npop.widget <- tkentry(ttibd, width = "20", textvariable = npop)
        npoplabel.widget <- tklabel(ttibd, text = "Number of populations:")
        freq.gridlabel.widget <- tklabel(ttibd, text = "Return frequencies on grid:")
        freq.grid <- tclVar("FALSE")
        wfreq.grid <- .Tk.subwin(ttibd)
        freq.gridoptionmenu.widget <- tcl("tk_optionMenu", wfreq.grid, 
            freq.grid, "FALSE", "TRUE")
        tess.gridlabel.widget <- tklabel(ttibd, text = "Return population membership on grid:")
        tess.grid <- tclVar("FALSE")
        wtess.grid <- .Tk.subwin(ttibd)
        tess.gridoptionmenu.widget <- tcl("tk_optionMenu", wtess.grid, 
            tess.grid, "FALSE", "TRUE")
        npixh <- tclVar(50)
        npixv <- tclVar(50)
        npixh.widget <- tkentry(ttibd, width = "8", textvariable = npixh)
        npixv.widget <- tkentry(ttibd, width = "8", textvariable = npixv)
        npixlabel.widget <- tklabel(ttibd, text = "Number of pixels for representation (hor|ver):")
        comp.Fstlabel.widget <- tklabel(ttibd, text = "Compute F statistics:")
        comp.Fst <- tclVar("FALSE")
        wcomp.Fst <- .Tk.subwin(ttibd)
        comp.Fstoptionmenu.widget <- tcl("tk_optionMenu", wcomp.Fst, 
            comp.Fst, "FALSE", "TRUE")
        comp.Dsigma2label.widget <- tklabel(ttibd, text = "Compute IBD index Dsigma2:")
        comp.Dsigma2 <- tclVar("FALSE")
        wcomp.Dsigma2 <- .Tk.subwin(ttibd)
        comp.Dsigma2optionmenu.widget <- tcl("tk_optionMenu", 
            wcomp.Dsigma2, comp.Dsigma2, "FALSE", "TRUE")
        comp.heightlabel.widget <- tklabel(ttibd, text = "Index of allele freq. variation:")
        comp.height <- tclVar("FALSE")
        wcomp.height <- .Tk.subwin(ttibd)
        comp.heightoptionmenu.widget <- tcl("tk_optionMenu", 
            wcomp.height, comp.height, "FALSE", "TRUE")
        hwidth <- tclVar(0.1)
        hwidth.widget <- tkentry(ttibd, width = "20", textvariable = hwidth)
        hwidthlabel.widget <- tklabel(ttibd, text = "Width around the barriers:")
        plot.pairs.borders <- tclVar("FALSE")
        prevcoord <- tclVar(0)
        prevcoordlabel.widget <- tklabel(ttibd, text = "Use loaded coordinate file:")
        prevcoord.widget <- tkcheckbutton(ttibd, variable = save, 
            onvalue = 1, offvalue = 0)
        save <- tclVar(0)
        savelabel.widget <- tklabel(ttibd, text = "Save coordinate and genotype files:")
        save.widget <- tkcheckbutton(ttibd, variable = save, 
            onvalue = 1, offvalue = 0)
        tkgrid(nindivlabel.widget, row = 1, column = 1, sticky = "w")
        tkgrid(nindiv.widget, row = 1, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(coordxlabel.widget, row = 2, column = 1, sticky = "w")
        tkgrid(abslabel.widget, row = 3, column = 1, sticky = "w")
        tkgrid(absmin.widget, row = 3, column = 2, sticky = "w")
        tkgrid(absmax.widget, row = 3, column = 3, sticky = "w")
        tkgrid(ordlabel.widget, row = 4, column = 1, sticky = "w")
        tkgrid(ordmin.widget, row = 4, column = 2, sticky = "w")
        tkgrid(ordmax.widget, row = 4, column = 3, sticky = "w")
        tkgrid(nucleilabel.widget, row = 5, column = 1, sticky = "w")
        tkgrid(nuclei.widget, row = 5, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(nloclabel.widget, row = 6, column = 1, sticky = "w")
        tkgrid(nloc.widget, row = 6, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(nalllabel.widget, row = 7, column = 1, sticky = "w")
        tkgrid(nall.widget, row = 7, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(betalabel.widget, row = 9, column = 1, columnspan = 2, 
            sticky = "w")
        tkgrid(beta.widget, row = 9, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(npoplabel.widget, row = 10, column = 1, sticky = "w")
        tkgrid(npop.widget, row = 10, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(freq.gridlabel.widget, row = 11, column = 1, sticky = "w")
        tkgrid(wfreq.grid, row = 11, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(tess.gridlabel.widget, row = 12, column = 1, sticky = "w")
        tkgrid(wtess.grid, row = 12, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(npixlabel.widget, row = 13, column = 1, sticky = "w")
        tkgrid(npixh.widget, row = 13, column = 2, sticky = "w")
        tkgrid(npixv.widget, row = 13, column = 3, sticky = "w")
        tkgrid(comp.Fstlabel.widget, row = 14, column = 1, sticky = "w")
        tkgrid(wcomp.Fst, row = 14, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(comp.Dsigma2label.widget, row = 15, column = 1, 
            sticky = "w")
        tkgrid(wcomp.Dsigma2, row = 15, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(comp.heightlabel.widget, row = 16, column = 1, 
            sticky = "w")
        tkgrid(wcomp.height, row = 16, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(hwidthlabel.widget, row = 17, column = 1, sticky = "w")
        tkgrid(hwidth.widget, row = 17, column = 2, columnspan = 2, 
            sticky = "w")
        tkgrid(prevcoordlabel.widget, row = 19, column = 1, sticky = "w")
        tkgrid(prevcoord.widget, row = 19, column = 2, sticky = "w")
        tkgrid(savelabel.widget, row = 20, column = 1, sticky = "w")
        tkgrid(save.widget, row = 20, column = 2, columnspan = 2, 
            sticky = "w")
        labelspace <- tklabel(ttibd, text = " ")
        tkgrid(labelspace, row = 21, column = 1)
        nextbutton <- tkbutton(ttibd, image = imagerun2, text = "RUN >>", 
            command = runibd)
        tkgrid(nextbutton, row = 22, column = 2, columnspan = 2, 
            sticky = "e")
    }
    Nullify <- function() {
        ttnul <- tktoplevel()
        tkwm.title(ttnul, "Simulate genotype with null alleles from loaded dataset")
        gltgp <- function() {
            if (length(globalgenotypes) == 0) {
                tkmessageBox(message = "You must define genotype file", 
                  icon = "error", type = "ok")
            }
            else {
                tttry <- tktoplevel(parent = .TkRoot)
                tkgrab(tttry)
                tkwm.geometry(tttry, "+200+200")
                tkwm.title(tttry, "wait")
                warn <- tklabel(tttry, image = imagepleasewait)
                tkpack(warn)
                tkfocus(tttry)
                tcl("update")
                print("Starting...")
                Sys.sleep(0.5)
                err <- try(nullify(genotypes = globalgenotypes, 
                  nall.null = as.integer(tclvalue(nall)), nloc.null = as.integer(tclvalue(nloc))), 
                  silent = TRUE)
                tkdestroy(tttry)
                print("Done.")
                if (class(err) == "try-error") {
                  tkmessageBox(message = err, icon = "error", 
                    type = "ok", parent = tt)
                  Log(paste("nullify(genotypes=", matrix2str(globalgenotypes), 
                    ",nall.null=", as.integer(tclvalue(nall)), 
                    ",nloc.null=", as.integer(tclvalue(nloc)), 
                    ")", sep = ""), "[FAILED] ")
                }
                else {
                  Log(paste("nullify(genotypes=", matrix2str(globalgenotypes), 
                    ",nall.null=", as.integer(tclvalue(nall)), 
                    ",nloc.null=", as.integer(tclvalue(nloc)), 
                    ")", sep = ""), "[SUCCESS] ")
                  tkmessageBox(message = "Terminated with success", 
                    type = "ok", parent = tt)
                  globalgenotypes <<- err$genotypes
                  if (tclvalue(save) == 1) {
                    if (tclvalue(sep2) == "White space") 
                      tclvalue(sep2) <- " "
                    auxgen <- tclVar()
                    tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes = "{{All files} *}", 
                      initialdir = tclvalue(outputdir), title = "Save genotype file to:"))
                    write.table(err$genotypes, file = tclvalue(auxgen), 
                      sep = tclvalue(sep2), row.names = FALSE, 
                      col.names = FALSE)
                    if (tclvalue(sep2) == " ") 
                      tclvalue(sep2) <- "White space"
                  }
                  tclvalue(labelgenotext) <- paste(tclvalue(labelgenotext), 
                    "with null alleles")
                }
            }
        }
        nall <- tclVar(1)
        nallentry.widget <- tkentry(ttnul, textvariable = nall, 
            width = 15)
        nalllabel.widget <- tklabel(ttnul, text = "Number of null alleles on each locus:")
        tkgrid(nalllabel.widget, row = 1, column = 1, sticky = "w")
        tkgrid(nallentry.widget, row = 1, column = 2, sticky = "w")
        nloc <- tclVar(1)
        nlocentry.widget <- tkentry(ttnul, textvariable = nloc, 
            width = 15)
        nloclabel.widget <- tklabel(ttnul, text = "Number of loci with null alleles:")
        tkgrid(nloclabel.widget, row = 2, column = 1, sticky = "w")
        tkgrid(nlocentry.widget, row = 2, column = 2, sticky = "w")
        save <- tclVar(0)
        savelabel.widget <- tklabel(ttnul, text = "Save genotype file:")
        save.widget <- tkcheckbutton(ttnul, variable = save, 
            onvalue = 1, offvalue = 0)
        tkgrid(savelabel.widget, row = 3, column = 1, sticky = "w")
        tkgrid(save.widget, row = 3, column = 2, sticky = "w")
        labelspace <- tklabel(ttnul, text = " ")
        tkgrid(labelspace, row = 4, column = 1)
        nextbutton <- tkbutton(ttnul, image = imageok, text = "RUN >>", 
            command = gltgp)
        tkgrid(nextbutton, row = 5, column = 2, sticky = "e")
        tkfocus(ttnul)
    }
    Reset <- function() {
        auxblink <<- 0
        idb.dataset <<- 0
        globalcoordinates <<- NULL
        globalgenotypes <<- 0
        globallabels <<- NA
        tclvalue(coordinatesfile) <<- ""
        tclvalue(genotypefile) <<- ""
        tclvalue(outputdir) <<- ""
        tclvalue(advanced) <<- 0
        tclvalue(burnin) <<- 0
        tclvalue(labelcoordtext) = ""
        tclvalue(labelgenotext) = ""
        configure()
        run()
        postproc()
        plot()
        SimnonIBD()
        SimIBD()
        plot2()
    }
    Showtextproba <- function(filename) {
        tttry <- tktoplevel(parent = .TkRoot)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry, "wait")
        warn <- tklabel(tttry, image = imagepleasewait)
        tkpack(warn)
        tkwait.visibility(tttry)
        tkfocus(tttry)
        tcl("update")
        Sys.sleep(0.5)
        file <- try(read.table(paste(tclvalue(outputdir), filename, 
            sep = "")), silent = TRUE)
        if (class(file) == "try-error") {
            tkmessageBox(message = "File hasn't been created or bad output path", 
                type = "ok", parent = tt)
        }
        else {
            file2 <- try(readLines(paste(tclvalue(outputdir), 
                "postprocess.parameters.txt", sep = "")), silent = TRUE)
            tttext <- tktoplevel(parent = .TkRoot)
            tkwm.title(tttext, filename)
            posx <- tclVar("")
            posy <- tclVar("")
            left <- tktext(tttext)
            txt <- tktext(tttext)
            top <- tktext(tttext)
            yscr <- tkscrollbar(tttext, repeatinterval = 5, command = function(...) {
                tkyview(txt, ...)
                tkyview(left, ...)
            })
            xscr <- tkscrollbar(tttext, repeatinterval = 5, orient = "horizontal", 
                command = function(...) {
                  tkxview(txt, ...)
                  tkxview(top, ...)
                })
            tkconfigure(txt, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(left, as.double(...))
                }, xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(top, as.double(...))
                })
            row <- NROW(file)
            col <- NCOL(file)
            tkconfigure(top, font = tkfont.create(family = "courrier"), 
                height = 1, wrap = "none", xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(txt, as.double(...))
                })
            auxtxt <- ""
            if (class(file2) == "try-error") {
                tkconfigure(left, font = tkfont.create(family = "courrier"), 
                  wrap = "none", width = numberofdigits(row + 
                    6), yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(txt, as.double(...))
                  })
                for (i in 1:row) {
                  auxtxt <- paste(auxtxt, "pixel ", sep = "")
                  auxtxt <- paste(auxtxt, i, sep = "")
                  auxtxt <- paste(auxtxt, "\n", sep = "")
                }
            }
            else {
                tkconfigure(left, font = tkfont.create(family = "courrier"), 
                  wrap = "none", yscrollcommand = function(...) {
                    tkset(yscr, ...)
                    tkyview.moveto(txt, as.double(...))
                  })
                aux1 <- as.integer(substr(file2[1], 9, 15))
                aux2 <- as.integer(substr(file2[2], 9, 15))
                for (i in 1:aux1) {
                  for (j in 1:aux2) {
                    auxtxt <- paste(auxtxt, "pixel ", sep = "")
                    auxtxt <- paste(auxtxt, i, sep = "")
                    auxtxt <- paste(auxtxt, "x", sep = "")
                    auxtxt <- paste(auxtxt, j, sep = "")
                    auxtxt <- paste(auxtxt, "\n", sep = "")
                  }
                }
                tkconfigure(left, width = (numberofdigits(aux1) + 
                  numberofdigits(aux2) + 7))
            }
            tkinsert(left, "end", auxtxt)
            auxtxt <- ""
            tkinsert(top, "end", "x coordinate\ty coordinate\t")
            for (j in 3:col) {
                auxtxt <- paste(auxtxt, paste("pop", as.character(j - 
                  2), sep = ""), sep = "")
                auxtxt <- paste(auxtxt, "\t\t", sep = "")
            }
            tkinsert(top, "end", auxtxt)
            auxtxt <- ""
            for (i in 1:row) {
                for (j in 1:col) {
                  if (as.double(file[i, j]) < 0) 
                    auxtxt <- paste(auxtxt, sprintf("%.4f", as.double(file[i, 
                      j])), sep = "")
                  else auxtxt <- paste(auxtxt, sprintf("%.5f", 
                    as.double(file[i, j])), sep = "")
                  auxtxt <- paste(auxtxt, "\t\t", sep = "")
                }
                auxtxt <- paste(auxtxt, "\n", sep = "")
            }
            tkinsert(txt, "end", auxtxt)
            tkdestroy(tttry)
            tkgrid(top, row = 1, column = 2)
            tkgrid(txt, row = 2, column = 2)
            tkgrid(left, row = 2, column = 1)
            tkgrid(yscr, row = 2, column = 3, sticky = "ns")
            tkgrid(xscr, row = 3, column = 2, sticky = "we")
        }
    }
    Showtextindiv <- function(filename) {
        tttry <- tktoplevel(parent = .TkRoot)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry, "wait")
        warn <- tklabel(tttry, image = imagepleasewait)
        tkpack(warn)
        tkwait.visibility(tttry)
        tkfocus(tttry)
        tcl("update")
        Sys.sleep(0.5)
        file <- try(read.table(paste(tclvalue(outputdir), filename, 
            sep = "")), silent = TRUE)
        if (class(file) == "try-error") {
            tkmessageBox(message = "File hasn't been created or bad output path", 
                type = "ok", parent = tt)
        }
        else {
            tttext <- tktoplevel(parent = .TkRoot)
            tkwm.title(tttext, filename)
            posx <- tclVar("")
            posy <- tclVar("")
            txt <- tktext(tttext)
            left <- tktext(tttext)
            top <- tktext(tttext)
            yscr <- tkscrollbar(tttext, repeatinterval = 5, command = function(...) {
                tkyview(txt, ...)
                tkyview(left, ...)
            })
            xscr <- tkscrollbar(tttext, repeatinterval = 5, orient = "horizontal", 
                command = function(...) {
                  tkxview(txt, ...)
                  tkxview(top, ...)
                })
            tkconfigure(txt, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(left, as.double(...))
                }, xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(top, as.double(...))
                })
            row <- NROW(file)
            col <- NCOL(file)
            tkconfigure(top, font = tkfont.create(family = "courrier"), 
                height = 1, wrap = "none", xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(txt, as.double(...))
                })
            tkconfigure(left, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(txt, as.double(...))
                })
            auxtxt <- ""
            auxtxt <- paste(auxtxt, "x coordinate\ty coordinate", 
                sep = "")
            auxtxt <- paste(auxtxt, "\t", sep = "")
            for (j in 1:(col - 2)) {
                auxtxt <- paste(auxtxt, paste("pop", as.character(j), 
                  sep = ""), sep = "")
                auxtxt <- paste(auxtxt, "\t\t", sep = "")
            }
            tkinsert(top, "end", auxtxt)
            auxtxt <- ""
            for (i in 1:row) {
                for (j in 1:col) {
                  auxtxt <- paste(auxtxt, sprintf("%.5f", as.double(file[i, 
                    j])), sep = "")
                  auxtxt <- paste(auxtxt, "\t\t", sep = "")
                }
                auxtxt <- paste(auxtxt, "\n", sep = "")
            }
            tkinsert(txt, "end", auxtxt)
            auxtxt <- ""
            if (length(globallabels) == 1) {
                for (i in 1:row) {
                  auxtxt <- paste(auxtxt, "ind ", sep = "")
                  auxtxt <- paste(auxtxt, i, sep = "")
                  auxtxt <- paste(auxtxt, "\n", sep = "")
                }
                tkconfigure(left, width = numberofdigits(row) + 
                  4)
            }
            else {
                size <- 0
                for (i in 1:row) {
                  auxtxt <- paste(auxtxt, globallabels[i], sep = "")
                  auxtxt <- paste(auxtxt, "\n", sep = "")
                  if (size < nchar(globallabels[i], "width")) 
                    size <- nchar(globallabels[i], "width")
                }
                tkconfigure(left, width = size)
            }
            tkinsert(left, "end", auxtxt)
            tkdestroy(tttry)
            tkgrid(top, row = 1, column = 2)
            tkgrid(txt, row = 2, column = 2)
            tkgrid(left, row = 2, column = 1)
            tkgrid(yscr, row = 2, column = 3, sticky = "ns")
            tkgrid(xscr, row = 3, column = 2, sticky = "we")
        }
    }
    Showtextmodal <- function(filename) {
        tttry <- tktoplevel(parent = .TkRoot)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry, "wait")
        warn <- tklabel(tttry, image = imagepleasewait)
        tkpack(warn)
        tkwait.visibility(tttry)
        tkfocus(tttry)
        tcl("update")
        Sys.sleep(0.5)
        file <- try(read.table(paste(tclvalue(outputdir), filename, 
            sep = "")), silent = TRUE)
        if (class(file) == "try-error") {
            tkmessageBox(message = "File hasn't been created or bad output path", 
                type = "ok", parent = tt)
        }
        else {
            tttext <- tktoplevel(parent = .TkRoot)
            tkwm.title(tttext, filename)
            posx <- tclVar("")
            posy <- tclVar("")
            txt <- tktext(tttext)
            left <- tktext(tttext)
            top <- tktext(tttext)
            yscr <- tkscrollbar(tttext, repeatinterval = 5, command = function(...) {
                tkyview(txt, ...)
                tkyview(left, ...)
            })
            xscr <- tkscrollbar(tttext, repeatinterval = 5, orient = "horizontal", 
                command = function(...) {
                  tkxview(txt, ...)
                  tkxview(top, ...)
                })
            tkconfigure(txt, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(left, as.double(...))
                }, xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(top, as.double(...))
                })
            row <- NROW(file)
            col <- NCOL(file)
            tkconfigure(top, font = tkfont.create(family = "courrier"), 
                height = 1, wrap = "none", xscrollcommand = function(...) {
                  tkset(xscr, ...)
                  tkxview.moveto(txt, as.double(...))
                })
            tkconfigure(left, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(txt, as.double(...))
                })
            auxtxt <- ""
            auxtxt <- paste(auxtxt, "x coordinate\t", sep = "")
            auxtxt <- paste(auxtxt, "y coordinate\t", sep = "")
            auxtxt <- paste(auxtxt, "most likely population", 
                sep = "")
            tkinsert(top, "end", auxtxt)
            auxtxt <- ""
            for (i in 1:row) {
                for (j in 1:col) {
                  if (j == 3) 
                    auxtxt <- paste(auxtxt, file[i, j], sep = "")
                  else auxtxt <- paste(auxtxt, sprintf("%.8f", 
                    as.double(file[i, j])), sep = "")
                  auxtxt <- paste(auxtxt, "\t", sep = "")
                }
                auxtxt <- paste(auxtxt, "\n", sep = "")
            }
            tkinsert(txt, "end", auxtxt)
            auxtxt <- ""
            if (length(globallabels) == 1) {
                for (i in 1:row) {
                  auxtxt <- paste(auxtxt, "ind ", sep = "")
                  auxtxt <- paste(auxtxt, i, sep = "")
                  auxtxt <- paste(auxtxt, "\n", sep = "")
                }
                tkconfigure(left, width = numberofdigits(row) + 
                  4)
            }
            else {
                size <- 0
                for (i in 1:row) {
                  auxtxt <- paste(auxtxt, globallabels[i], sep = "")
                  auxtxt <- paste(auxtxt, "\n", sep = "")
                  if (size < nchar(globallabels[i], "width")) 
                    size <- nchar(globallabels[i], "width")
                }
                tkconfigure(left, width = size)
            }
            tkinsert(left, "end", auxtxt)
            tkdestroy(tttry)
            tkgrid(top, row = 1, column = 2)
            tkgrid(txt, row = 2, column = 2)
            tkgrid(left, row = 2, column = 1)
            tkgrid(yscr, row = 2, column = 3, sticky = "ns")
            tkgrid(xscr, row = 3, column = 2, sticky = "we")
        }
    }
    Showtextposterior <- function(filename) {
        tttry <- tktoplevel(parent = .TkRoot)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry, "wait")
        warn <- tklabel(tttry, image = imagepleasewait)
        tkpack(warn)
        tkwait.visibility(tttry)
        tkfocus(tttry)
        tcl("update")
        Sys.sleep(0.5)
        file <- try(scan(paste(tclvalue(outputdir), filename, 
            sep = "")), silent = TRUE)
        if (class(file) == "try-error") {
            tkmessageBox(message = "File hasn't been created or bad output path", 
                type = "ok", parent = tt)
        }
        else {
            tttext <- tktoplevel(parent = .TkRoot)
            tkwm.title(tttext, "Posterior density of model (values in log)")
            posx <- tclVar("")
            posy <- tclVar("")
            yscr <- tkscrollbar(tttext, repeatinterval = 5, command = function(...) {
                tkyview(txt, ...)
                tkyview(left, ...)
            })
            xscr <- tkscrollbar(tttext, repeatinterval = 5, orient = "horizontal", 
                command = function(...) {
                  tkxview(txt, ...)
                })
            txt <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(left, as.double(...))
                }, xscrollcommand = function(...) {
                  tkset(xscr, ...)
                })
            if (tclvalue(burnin) != 0) 
                row <- NROW(file[-(1:as.numeric(tclvalue(burnin)))])
            else row <- NROW(file)
            left <- tktext(tttext, font = tkfont.create(family = "courrier"), 
                wrap = "none", yscrollcommand = function(...) {
                  tkset(yscr, ...)
                  tkyview.moveto(txt, as.double(...))
                })
            if (tclvalue(burnin) != 0) 
                mean.lpp <- mean(file[-(1:as.numeric(tclvalue(burnin)))])
            else mean.lpp <- mean(file)
            auxtxt <- "\nMean= "
            auxtxt <- paste(auxtxt, mean.lpp, sep = "")
            auxtxt <- paste(auxtxt, "\n\nPoint values along the chain\n\n", 
                sep = "")
            for (i in 1:row) {
                auxtxt <- paste(auxtxt, file[i], sep = "")
                auxtxt <- paste(auxtxt, "\n", sep = "")
            }
            tkinsert(txt, "end", auxtxt)
            auxtxt <- paste("\nBurnin=", tclvalue(burnin), "\n----\n\n", 
                sep = "")
            for (i in 1:row) {
                auxtxt <- paste(auxtxt, "Sample ", sep = "")
                auxtxt <- paste(auxtxt, i, sep = "")
                auxtxt <- paste(auxtxt, "\n", sep = "")
            }
            tkconfigure(left, width = numberofdigits(row) + 7)
            tkinsert(left, "end", auxtxt)
            tkdestroy(tttry)
            tkgrid(txt, row = 2, column = 2)
            tkgrid(left, row = 2, column = 1)
            tkgrid(yscr, row = 2, column = 3, sticky = "ns")
            tkgrid(xscr, row = 3, column = 2, sticky = "we")
        }
    }
    Log <- function(line, state) {
        if (tclvalue(LogState) == 1) {
            if (tclvalue(outputdir) == "") 
                tkmessageBox(message = "Output Directory is not set. Log will not be written.", 
                  icon = "error", type = "ok", parent = tt)
            else {
                zz1 <- file(paste(tclvalue(outputdir), "ExecLog.txt", 
                  sep = ""), "a")
                cat(line, "\n", file = zz1)
                close(zz1)
                zz2 <- file(paste(tclvalue(outputdir), "StatesLog.txt", 
                  sep = ""), "a")
                cat(state, line, "\n", file = zz2)
                close(zz2)
            }
        }
    }
    Parallel <- function() {
        Setparallel <- function() {
            if (tclvalue(processors) < 1) {
                tkmessageBox(message = "Geneland doesn't support abacus", 
                  icon = "error", type = "ok", parent = ttpara)
                usecluster <<- FALSE
                tkconfigure(nextbutton, state = "normal")
                tkconfigure(cancelbutton, state = "disabled")
                tkconfigure(info.label, text = "Parallel mode = OFF")
            }
            else if (tclvalue(processors) == 1) {
                tkmessageBox(message = "For parallel processing it is advised to use at least two nodes", 
                  icon = "error", type = "ok", parent = ttpara)
                usecluster <<- FALSE
                tkconfigure(nextbutton, state = "normal")
                tkconfigure(cancelbutton, state = "disabled")
                tkconfigure(info.label, text = "Parallel mode = OFF")
            }
            else {
                cluster <<- try(makeCluster(as.numeric(tclvalue(processors)), 
                  type = tclvalue(pmethod)), silent = TRUE)
                if (class(cluster) == "try-error") {
                  tkmessageBox(message = paste("Error, please read snow package documentation.\n The received error message was:", 
                    cluster, sep = "\n"), icon = "error", type = "ok", 
                    parent = ttpara)
                  usecluster <<- FALSE
                  tkconfigure(nextbutton, state = "normal")
                  tkconfigure(cancelbutton, state = "disabled")
                  tkconfigure(info.label, text = "Parallel mode = OFF")
                }
                else {
                  usecluster <<- TRUE
                  tkconfigure(nextbutton, state = "disabled")
                  tkconfigure(cancelbutton, state = "normal")
                  tkconfigure(info.label, text = "Parallel mode = ON")
                }
            }
        }
        if (require(snow) == FALSE) 
            tkmessageBox(message = "Snow not found. Install it before using this feature", 
                icon = "error", type = "ok", parent = ttpara)
        else {
            ttpara <- tktoplevel(parent = .TkRoot)
            tkwm.geometry(ttpara, "+200+100")
            tkwm.title(ttpara, "Parallel options")
            processors.label <- tklabel(ttpara, text = "Number of nodes:")
            processors.entry <- tkentry(ttpara, width = "10", 
                textvariable = processors)
            labelspace1 <- tklabel(ttpara, text = " ")
            method.label <- tklabel(ttpara, text = "Parallelization method:")
            pmethod <- tclVar("MPI")
            method.label1 <- tklabel(ttpara, text = "MPI")
            method.label2 <- tklabel(ttpara, text = "PVM")
            method.radio1 <- tkradiobutton(ttpara, variable = pmethod, 
                value = "MPI", selectcolor = "white")
            method.radio2 <- tkradiobutton(ttpara, variable = pmethod, 
                value = "PVM", selectcolor = "white")
            labelspace2 <- tklabel(ttpara, text = " ")
            info.label <- tklabel(ttpara, text = "", foreground = "blue")
            labelspace3 <- tklabel(ttpara, text = " ")
            nextbutton <- tkbutton(ttpara, text = "Start", command = Setparallel)
            cancelbutton <- tkbutton(ttpara, text = "Stop", command = function() {
                stopCluster(cluster)
                usecluster <<- FALSE
                tkconfigure(nextbutton, state = "normal")
                tkconfigure(cancelbutton, state = "disable")
                tkconfigure(info.label, text = "Parallel mode = OFF")
            })
            if (usecluster) {
                tkconfigure(nextbutton, state = "disabled")
                tkconfigure(cancelbutton, state = "normal")
                tkconfigure(info.label, text = "Parallel mode = ON")
            }
            else {
                tkconfigure(cancelbutton, state = "disabled")
                tkconfigure(nextbutton, state = "normal")
                tkconfigure(info.label, text = "Parallel mode = OFF")
            }
            tkgrid(processors.label, row = 1, column = 1, sticky = "e")
            tkgrid(processors.entry, row = 1, column = 2, columnspan = 3, 
                sticky = "e")
            tkgrid(labelspace1, row = 2, column = 1, columnspan = 4, 
                sticky = "e")
            tkgrid(method.label, row = 3, column = 1, rowspan = 2, 
                sticky = "ns")
            tkgrid(method.label1, row = 3, column = 3, sticky = "we")
            tkgrid(method.label2, row = 3, column = 4, sticky = "we")
            tkgrid(method.radio1, row = 4, column = 3, sticky = "we")
            tkgrid(method.radio2, row = 4, column = 4, sticky = "we")
            tkgrid(labelspace2, row = 5, column = 1, columnspan = 4, 
                sticky = "e")
            tkgrid(info.label, row = 5, column = 1, columnspan = 4, 
                sticky = "we")
            tkgrid(labelspace3, row = 7, column = 1, columnspan = 4, 
                sticky = "e")
            tkgrid(nextbutton, row = 8, column = 1, columnspan = 2, 
                sticky = "we")
            tkgrid(cancelbutton, row = 8, column = 3, columnspan = 2, 
                sticky = "we")
            tkfocus(ttpara)
        }
    }
    initialimage()
    topMenu <- tkmenu(tt)
    tkconfigure(tt, menu = topMenu)
    fileMenu <- tkmenu(topMenu, tearoff = FALSE)
    coordinatesMenu <- tkmenu(topMenu, tearoff = FALSE)
    genotypesMenu <- tkmenu(topMenu, tearoff = FALSE)
    missingMenu <- tkmenu(topMenu, tearoff = FALSE)
    toolsMenu <- tkmenu(topMenu, tearoff = FALSE)
    helpMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(fileMenu, "checkbutton", label = "Advanced options", 
        variable = advanced, selectcolor = "blue", onvalue = 1, 
        offvalue = 0, command = function() fadvanced())
    tkadd(coordinatesMenu, "radiobutton", label = "White space", 
        value = "", variable = sep1, selectcolor = "blue")
    tkadd(coordinatesMenu, "radiobutton", label = ",", value = ",", 
        variable = sep1, selectcolor = "blue")
    tkadd(coordinatesMenu, "radiobutton", label = ";", value = ";", 
        variable = sep1, selectcolor = "blue")
    tkadd(genotypesMenu, "radiobutton", label = "White space", 
        value = "", variable = sep2, selectcolor = "blue")
    tkadd(genotypesMenu, "radiobutton", label = ",", value = ",", 
        variable = sep2, selectcolor = "blue")
    tkadd(genotypesMenu, "radiobutton", label = ";", value = ";", 
        variable = sep2, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "0", value = "0", 
        variable = md, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "00", value = "00", 
        variable = md, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "000", value = "000", 
        variable = md, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "A", value = "A", 
        variable = md, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "NA", value = "NA", 
        variable = md, selectcolor = "blue")
    tkadd(missingMenu, "radiobutton", label = "NAD", value = "NAD", 
        variable = md, selectcolor = "blue")
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "cascade", label = "Missing data symbol", 
        menu = missingMenu)
    tkadd(fileMenu, "cascade", label = "Coordinate file values separator", 
        menu = coordinatesMenu)
    tkadd(fileMenu, "cascade", label = "Genotype file values separator", 
        menu = genotypesMenu)
    tkadd(toolsMenu, "command", label = "Convert to Genepop files", 
        command = function() Convert())
    tkadd(toolsMenu, "command", label = "Simulate null alleles", 
        state = "disabled", command = function() Nullify())
    tkadd(toolsMenu, "command", label = "Use parallel processing", 
        command = function() Parallel())
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "checkbutton", label = "Create log file", 
        variable = LogState, offvalue = 0, onvalue = 1, selectcolor = "blue", 
        command = function() {
            if (tclvalue(LogState) == 1) {
                if (tclvalue(outputdir) == "") {
                  tkmessageBox(message = "Output Directory is not set. Log will not be written.", 
                    icon = "error", type = "ok", parent = tt)
                }
            }
        })
    tkadd(fileMenu, "command", label = "Reset values", command = function() Reset())
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Quit", command = function() tkdestroy(tt))
    tkadd(topMenu, "cascade", label = "Menu", menu = fileMenu)
    tkadd(topMenu, "cascade", label = "Tools", menu = toolsMenu)
    tkadd(helpMenu, "command", label = "Help", command = function() helpWindow())
    tkadd(helpMenu, "command", label = "Credits", command = function() creditsWindow())
    tkadd(topMenu, "cascade", label = "Help", menu = helpMenu)
    labelinference <- tklabel(ttpan, text = "-Inference-", font = "*-Times-bold-i-normal--20-*", 
        foreground = "blue")
    labelsimulation <- tklabel(ttpan, text = "-Simulation-", 
        state = "disabled", font = "*-Times-bold-i-normal--20-*", 
        foreground = "blue")
    labelspace <- tklabel(ttpan, text = " ")
    buttonconf <- tkbutton(ttpan, image = imageconfigure, text = "Configure", 
        command = function() {
            configure()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttpost)
            tkgrid.remove(ttplot)
            tkgrid(ttconf, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonrun <- tkbutton(ttpan, image = imagerun, text = "Run", 
        command = function() {
            run()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttpost)
            tkgrid.remove(ttplot)
            tkgrid(ttrun, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonpostprocess <- tkbutton(ttpan, image = imagepostprocess, 
        text = "Postprocess", command = function() {
            postproc()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttplot)
            tkgrid(ttpost, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonsimfmodel <- tkbutton(ttpan, image = imagefmodel, text = "F-model", 
        state = "disabled", command = function() {
            SimnonIBD()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttplot)
            tkgrid(ttsimf, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonplot <- tkbutton(ttpan, image = imageplot, text = "Plot", 
        command = function() {
            plot()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttpost)
            tkgrid(ttplot, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonibd <- tkbutton(ttpan, image = imageibd, text = "IBD", 
        state = "disabled", command = function() {
            SimIBD()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttplot)
            tkgrid.remove(ttplot2)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttpost)
            tkgrid(ttibd, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    buttonplot2 <- tkbutton(ttpan, image = imageplot, text = "Plot2", 
        state = "disabled", command = function() {
            plot2()
            tkgrid.remove(ttinit)
            tkgrid.remove(ttconf)
            tkgrid.remove(ttsimf)
            tkgrid.remove(ttibd)
            tkgrid.remove(ttplot)
            tkgrid.remove(ttfstat)
            tkgrid.remove(ttrun)
            tkgrid.remove(ttpost)
            tkgrid(ttplot2, row = 1, column = 2, sticky = "we", 
                pady = 10)
        })
    tkgrid(labelinference, row = 1, column = 1, sticky = "w", 
        padx = 10, pady = 10)
    tkgrid(buttonconf, row = 2, column = 1, sticky = "we", padx = 10)
    tkgrid(buttonrun, row = 3, column = 1, sticky = "we", padx = 10)
    tkgrid(buttonpostprocess, row = 4, column = 1, sticky = "we", 
        padx = 10)
    tkgrid(buttonplot, row = 5, column = 1, sticky = "we", padx = 10)
    tkgrid(labelspace, row = 6, column = 1, sticky = "w", padx = 10)
    tkgrid(labelsimulation, row = 7, column = 1, sticky = "w", 
        padx = 10)
    tkgrid(buttonsimfmodel, row = 8, column = 1, sticky = "we", 
        padx = 10)
    tkgrid(buttonplot2, row = 10, column = 1, sticky = "we", 
        padx = 10)
    coordownlabel.widget <- tklabel(tt, textvariable = labelcoordtext, 
        foreground = "blue")
    genodownlabel.widget <- tklabel(tt, textvariable = labelgenotext, 
        foreground = "blue")
    auxblink <- 1
    extralabel.widget <- tklabel(tt, text = "Please configure output directory", 
        foreground = "blue")
    tkbind(extralabel.widget, "<Button-1>", function() {
        configure()
        tkgrid.remove(ttinit)
        tkgrid.remove(ttsimf)
        tkgrid.remove(ttibd)
        tkgrid.remove(ttfstat)
        tkgrid.remove(ttplot2)
        tkgrid.remove(ttrun)
        tkgrid.remove(ttpost)
        tkgrid.remove(ttplot)
        tkgrid(ttconf, row = 1, column = 2, sticky = "we")
    })
    blink <- function() {
        if (auxblink == 1) {
            error <- try(tkconfigure(extralabel.widget, text = ""), 
                silent = TRUE)
            if (class(error) == "try-error") 
                tkdestroy(tt)
            auxblink <<- 0
        }
        else if (auxblink == 0) {
            error <- try(tkconfigure(extralabel.widget, text = "Please configure output directory"), 
                silent = TRUE)
            if (class(error) == "try-error") 
                tkdestroy(tt)
            auxblink <<- 1
        }
        try(tcl("after", 1000, blink), silent = TRUE)
    }
    blink()
    tkgrid(ttpan, row = 1, column = 1, sticky = "we")
    tkgrid(ttinit, row = 1, column = 2, sticky = "we", pady = 10)
    tkgrid(coordownlabel.widget, row = 2, column = 1, columnspan = 2, 
        sticky = "w", padx = 3)
    tkgrid(genodownlabel.widget, row = 3, column = 1, columnspan = 2, 
        sticky = "w", padx = 3)
    tkgrid(extralabel.widget, row = 4, column = 1, columnspan = 2, 
        sticky = "w", padx = 3)
    tkgrid.columnconfigure(tt, 1, minsize = 200)
    tkgrid.columnconfigure(tt, 2, minsize = 500)
    tkgrid.rowconfigure(tt, 1, minsize = 520)
}
