`Geneland.GUI` <-
function(lib.loc=NULL) {

  
  
  #first the initial image
  require(tcltk)
  tt <- tktoplevel()
  tkwm.title(tt,"Geneland - Graphical Interface")
  tkwm.geometry(tt, "+100+100")
  
  image1 <- tclVar()
  tcl("image","create","photo",image1,file=system.file("images/geneland6.gif", package="Geneland",lib.loc=lib.loc))
  imgAsLabel <- tklabel(tt,image=image1)
  tkpack(imgAsLabel)
  tkfocus(tt)
  tcl("after","3000","destroy",tt)
  
  require(Geneland)
  #idb global object
  idb.dataset <- 0

  globalcoordinates <- 0
  globalgenotypes <- 0
  globalallele.numbers <- 0

  falush <- 0
  
  labelcoordtext <- tclVar("")
  labelgenotext <- tclVar("")

  
  #images for pretty interface
  imageconfigure <- tclVar()
  imagerun  <- tclVar()
  imagerun2  <- tclVar()
  imagefmodel  <- tclVar()
  imagefstat  <- tclVar()
  imageplot  <- tclVar()
  imagepostprocess  <- tclVar()
  imagedraw  <- tclVar()
  imageok  <- tclVar()
  imagepleasewait  <- tclVar()
  imageibd <- tclVar()
  imageconvert <- tclVar()
  #liblo
  
    #  imagecancel  <- tclVar()
  tcl("image","create","photo",imageconfigure,file=system.file("images/icon-configure.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagerun,file=system.file("images/icon-run.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagerun2,file=system.file("images/icon-run2.gif", package="Geneland"))
  tcl("image","create","photo",imagefmodel,file=system.file("images/icon-fmodel.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagefstat,file=system.file("images/icon-fstat.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imageplot,file=system.file("images/icon-output.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagepostprocess,file=system.file("images/icon-postprocess.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagedraw,file=system.file("images/icon-draw.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imageok,file=system.file("images/icon-ok.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imageibd,file=system.file("images/icon-ibd.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imagepleasewait,file=system.file("images/pleasewait.gif", package="Geneland",lib.loc=lib.loc))
  tcl("image","create","photo",imageconvert,file=system.file("images/icon-convert.gif", package="Geneland",lib.loc=lib.loc))
  
  tkwait.window(tt)
    # and now tcharaaaa

    #global variables here, please
  
  coordinatesfile <- tclVar("")
  genotypefile <- tclVar("")  
  outputdir <- tclVar("")
  advanced <- tclVar(0)
  sep1 <- tclVar("White space")
  sep2 <- tclVar("White space")
  
  #the "code" that really matters




  tt <- tktoplevel()
  tkwm.title(tt,"Geneland - Graphical Interface")
  tkwm.geometry(tt, "750x600+100+100") 
  tkwm.resizable(tt,0,0)

  
  ttconf <- tkframe(tt,borderwidth=2,relief="sunken")
  ttrun <- tkframe(tt,borderwidth=2,relief="sunken")
  ttinit <- tkframe(tt,borderwidth=2,relief="sunken")
  ttpost <- tkframe(tt,borderwidth=2,relief="sunken")
  ttsimf <- tkframe(tt,borderwidth=2,relief="sunken")
  ttplot <- tkframe(tt,borderwidth=2,relief="sunken")
  ttfstat <- tkframe(tt,borderwidth=2,relief="sunken")
  ttibd <- tkframe(tt,borderwidth=2,relief="sunken")
  ttplot2 <- tkframe(tt,borderwidth=2,relief="sunken")
  ttpan <- tkframe(tt)

  #tkoptionmenu <- function(...) tcl("tk_optionMenu", ...)
  #tkspinbox <- function(...) tcl("spinbox", ...)


  
  fadvanced <- function(){
    if(tclvalue(advanced)==1){
      tkconfigure(labelsimulation,state="normal")
      tkconfigure(buttonsimfmodel,state="normal")
      tkconfigure(buttonibd,state="normal")
      tkconfigure(buttonplot2,state="normal")
      #tkmenu.entryconfigure(toolsMenu,"1",state="normal")
      tcl(toolsMenu,"entryconfigure","1",state="normal")

    }
    else{
      tkconfigure(labelsimulation,state="disabled")
      tkconfigure(buttonsimfmodel,state="disabled")
      tkconfigure(buttonibd,state="disabled")
      tkconfigure(buttonplot2,state="disabled")
      tcl(toolsMenu,"entryconfigure","1",state="disabled")
    }
      
    run() 
  }
  
  
    
     # -----------------------------------------------
     # HELLO | HELP STARTS HERE
     #------------------------------------



  helpWindow  <- function() {
    help.start(Geneland)
    tkfocus(tt)
  }

     # -----------------------------------------------
     # HELLO | CREDITS STARTS HERE
     #------------------------------------


    creditsWindow  <- function() {

      ttcredits <- tktoplevel(parent=.TkRoot)
      tkwm.title(ttcredits,"Credits")

      label1.widget <- tklabel(ttcredits,text="Authors:")
      label2.widget <- tklabel(ttcredits,text=" ")
      label3.widget <- tklabel(ttcredits,text="Gilles Guillot:")
      label4.widget <- tklabel(ttcredits,text=" Fortran and R Code")
      label5.widget <- tklabel(ttcredits,text="")
      label6.widget <- tklabel(ttcredits,text="Filipe Santos:")
      label7.widget <- tklabel(ttcredits,text=" Graphical Interface (code and design)")
      label8.widget <- tklabel(ttcredits,text="")
      label9.widget <- tklabel(ttcredits,text="Arnaud Estoup:")
      label10.widget <- tklabel(ttcredits,text=" Graphical Interface (design and test)")

      tkgrid(label1.widget,row=1,column=1,sticky="w")
      tkgrid(label2.widget,row=2,column=1,sticky="w")
      tkgrid(label3.widget,row=3,column=1,sticky="w")
      tkgrid(label4.widget,row=3,column=2,sticky="w")
      tkgrid(label5.widget,row=4,column=1,sticky="w")
      tkgrid(label6.widget,row=5,column=1,sticky="w")
      tkgrid(label7.widget,row=5,column=2,sticky="w")
      tkgrid(label8.widget,row=6,column=1,sticky="w")
      tkgrid(label9.widget,row=7,column=1,sticky="w")
      tkgrid(label10.widget,row=7,column=2,sticky="w")
      
      
      #tkmessageBox(title="Credits",message="Authors:\n\t Gilles Guillot\n\n\t Filipe Santos",icon="info",type="ok")
      #tkfocus(tt)
    }





     # -----------------------------------------------
     # HELLO | CONFIGURE STARTS HERE
     #------------------------------------


  configure<-function(){
    
    #load files
    label1.widget <-tklabel(ttconf,text=tclvalue(coordinatesfile),width=40)
    tkconfigure(label1.widget,textvariable=coordinatesfile)
    getcoordinatesfile <- function()  {
      tclvalue(coordinatesfile) <- tclvalue(tkgetOpenFile(filetypes="{{All files} *}",title="Choose Coordinates File"))
      if(tclvalue(coordinatesfile)!=""){
        if(tclvalue(sep1)=="White space")
          tclvalue(sep1) <- ""
        globalcoordinates <<- read.table(tclvalue(coordinatesfile),sep= tclvalue(sep1))
        if(tclvalue(sep1)=="")
          tclvalue(sep1) <- "White space"
        tclvalue(labelcoordtext)="Coordinates: File loaded"
        #tkconfigure(coordownlabel.widget,text="Coordinates: File loaded")
      } else {
        globalcoordinates <<- 0
        tclvalue(labelcoordtext)="Coordinates: Data unloaded"
        #tkconfigure(coordownlabel.widget,text="Coordinates: File unloaded")
      }
      tkfocus(tt)
      }
    
    button1.widget <- tkbutton(ttconf,text="Coordinates File",command=getcoordinatesfile,width=15)
    
    
    tkgrid(button1.widget,row=2,column=1,sticky="we")
    tkgrid(label1.widget,row=2,column=2,sticky="we")
    
    label2.widget <- tklabel(ttconf,text=tclvalue(genotypefile),width=40)
    tkconfigure(label2.widget,textvariable=genotypefile)
    getgenotypefile <- function()  {
      tclvalue(genotypefile) <- tclvalue(tkgetOpenFile(filetypes="{{All files} *}",title="Choose Genotypes File"))
      
      if(tclvalue(genotypefile)!=""){
        if(tclvalue(sep2)=="White space")
          tclvalue(sep2) <- ""
        globalgenotypes <<- read.table(tclvalue(genotypefile),sep= tclvalue(sep2))
        auxallele <- FormatGenotypes(globalgenotypes)
        globalgenotypes <<- auxallele$genotypes
        globalallele.numbers <<- auxallele$allele.numbers
        if(tclvalue(sep2)=="")
          tclvalue(sep2) <- "White space"
        tclvalue(labelgenotext)="Genotype:    File loaded"
        #tkconfigure(genodownlabel.widget,text="Genotype:    File loaded")
        
      } else {
        globalgenotypes <<- 0
        globalallele.numbers <<- 0
        tclvalue(labelgenotext)="Genotype:    Data unloaded"
        #tkconfigure(genodownlabel.widget,text="Genotype:    File unloaded")
      }
      tkfocus(tt)
    }
    button2.widget <- tkbutton(ttconf,text="Genotype File",command=getgenotypefile,width=15)
    
    tkgrid(button2.widget,row=3,column=1,sticky="we")
    tkgrid(label2.widget,row=3,column=2,sticky="we")
    
    
    label3.widget <- tklabel(ttconf,text=tclvalue(outputdir),width=40)
    tkconfigure(label3.widget,textvariable=outputdir)
    getoutputdir <- function()  {
      tclvalue(outputdir) <- tclvalue(tkchooseDirectory(parent=tt,title="Please choose an output directory"))
      if(tclvalue(outputdir)!=""){
        tcl("regsub","-all","\\\\",tclvalue(outputdir),"/",outputdir)
        tcl("append",outputdir,"/")
        auxblink<<-2
        tkconfigure(extralabel.widget,text="")
      }
      tkfocus(tt)
    }
    button3.widget <- tkbutton(ttconf,text="Output Directory",command=getoutputdir,width=15)
   
    tkgrid(button3.widget,row=4,column=1,sticky="we")
    tkgrid(label3.widget,row=4,column=2,sticky="we")
  }

     # -----------------------------------------------
     # HELLO | RUN STARTS HERE
     #------------------------------------

  
  run<-function(){
    
    
    
    RunmcmcFmodel <- function(){
      
      if(tclvalue(ploidy)=="Haploid")
        dploidy <- 1
      else
        dploidy <- 2
      
#      vect<- c()
      
 #     vec<- unlist(strsplit(tclvalue(allelenumber),","))
 #     for (i in 1:length(vec)) vect[i]<-as.numeric(vec[i])
      
      
      tttry <- tktoplevel()
      tkgrab(tttry)
      tkwm.geometry(tttry, "+200+200")
      tkwm.title(tttry,"wait")
      warn<-tklabel(tttry,image=imagepleasewait)
      tkpack(warn)

      tkfocus(tttry)

         
                
      print("Starting...")
      Sys.sleep(0.1)

      #if(tclvalue(cb1state)==0){
      #  err <- try(mcmcFmodel(coordinates=globalcoordinates,genotypes=globalgenotypes,ploidy=dploidy,path.mcmc=tclvalue(outputdir),rate.max=as.numeric(tclvalue(rate)),delta.coord=as.numeric(tclvalue(delta)),npopmin=as.numeric(tclvalue(npopmin)),npopinit=as.numeric(tclvalue(npopinit)),npopmax=as.numeric(tclvalue(npopmax)),nb.nuclei.max=as.numeric(tclvalue(nuclei)),nit=as.numeric(tclvalue(nit)),thinning=as.numeric(tclvalue(thinning)),freq.model=tclvalue(freq),varnpop=as.logical(tclvalue(varnpop)),spatial=as.logical(tclvalue(spatial)),jcf=as.logical(tclvalue(jcf))),silent=TRUE)
        
        #mcmcFmodel(coordinates=globalcoordinates,genotypes=globalgenotypes,ploidy=dploidy,path.mcmc=tclvalue(outputdir),rate.max=as.numeric(tclvalue(rate)),delta.coord=as.numeric(tclvalue(delta)),npopmin=as.numeric(tclvalue(npopmin)),npopinit=as.numeric(tclvalue(npopinit)),npopmax=as.numeric(tclvalue(npopmax)),nb.nuclei.max=as.numeric(tclvalue(nuclei)),nit=as.numeric(tclvalue(nit)),thinning=as.numeric(tclvalue(thinning)),freq.model=tclvalue(freq),varnpop=as.logical(tclvalue(varnpop)),spatial=as.logical(tclvalue(spatial)),jcf=as.logical(tclvalue(jcf)))
      #}
      #else{
      err <- try(mcmcFmodel(coordinates=globalcoordinates,genotypes=globalgenotypes,ploidy=dploidy,path.mcmc=tclvalue(outputdir),rate.max=as.numeric(tclvalue(rate)),delta.coord=as.numeric(tclvalue(delta)),npopmin=as.numeric(tclvalue(npopmin)),npopinit=as.numeric(tclvalue(npopinit)),npopmax=as.numeric(tclvalue(npopmax)),nb.nuclei.max=as.numeric(tclvalue(nuclei)),nit=as.numeric(tclvalue(nit)),thinning=as.numeric(tclvalue(thinning)),freq.model=tclvalue(freq),varnpop=as.logical(tclvalue(varnpop)),spatial=as.logical(tclvalue(spatial)),jcf=as.logical(tclvalue(jcf))),silent=TRUE)
      #}
      tkdestroy(tttry)
      print("Done.")
    # validate <- 0
      if (class(err) == "try-error") 
        tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
      else{
        tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        if(tclvalue(freq)=="Falush")
          falush <<- 1
        else
          falush <<- 0
      }
    }
    #and now the pick of the allele.numbers

    #print(globalallele.numbers)
    
    #strallele <- paste(globalallele.numbers, sep = ",", collapse =",")#globalallele.numbers[0]
    #for (i in 1:length(globalallele.numbers)) strallele <-"," +strallele   
    #print(strallele)
    
    #cb1state <- tclVar()
    #allelenumber <- tclVar(strallele)
    #annoyingnumber <- 0
    #setallelenumber <- function()  {
    #  annoyingnumber <<- annoyingnumber + 1
    #  if(annoyingnumber==20){
    #    tkmessageBox(title="Hey you",message="For God sake. Stop that.",icon="info",type="ok")
    #  }
    #  else if(annoyingnumber==30){
    #      tkmessageBox(title="",message="I will bless you.",icon="info",type="ok")
    #    }
    #  else if(annoyingnumber>30){
    #    tkmessageBox(title="",message="There is no evolution.",icon="info",type="ok")
    #  }
    #  if(tclvalue(cb1state)==1){
    #    tkconfigure(entry1.widget,state="normal")
    #  }else{
    #    tkconfigure(entry1.widget,state="disable")
    #  }
    #  
    #}
    #entry1.widget <-tkentry(ttrun,width="20",textvariable=allelenumber,state="disable")
    #cb1.widget <- tkcheckbutton(ttrun,command=setallelenumber,variable=cb1state,onvalue=1, offvalue=0,selectcolor="blue")
    #tkdeselect(cb1.widget)
    #alellelabel.widget <- tklabel(ttrun,text="Number of alleles per locus (E.g: n,m,..):")



    #ploidy
    ploidylabel.widget <- tklabel(ttrun,text="Ploidy:")
    ploidy <- tclVar("Diploid")
    wdiplody <- .Tk.subwin(ttrun)
    ploidyoptionmenu.widget <- tcl("tk_optionMenu",wdiplody,ploidy,"Diploid","Haploid")
    #ploidycombo.widget <-
    #tk("tk_optionMenu","menubuttonoptions","haploid","diploid")
    tkgrid(ploidylabel.widget,row=2,column=1,sticky="w")
    tkgrid(wdiplody,row=2,column=3,sticky="w")

    #rate max
    rate <- tclVar(100)
    rate.widget <-tkentry(ttrun,width="20",textvariable=rate)
    ratelabel.widget <- tklabel(ttrun,text="Maximum rate of poisson process:")




    #delta.coord
    delta <- tclVar(0)
    delta.widget <-tkentry(ttrun,width="20",textvariable=delta)
    deltalabel.widget <- tklabel(ttrun,text="Uncertainty on coordinates:")
    
    tkgrid(deltalabel.widget,row=4,column=1,sticky="w")
    tkgrid(delta.widget,row=4,column=3,sticky="w")

    #populations
    npopmin <- tclVar("1")
    npopinit <- tclVar("1") 
    npopmax <- tclVar("1")
    
    wmin  <- .Tk.subwin(ttrun)
    winit <- .Tk.subwin(ttrun)
    wmax <- .Tk.subwin(ttrun)
    
    actualizaspin<-function(){
      tkconfigure(min,"-to",tclvalue(npopmax))
      tkconfigure(max,"-from",tclvalue(npopmin))
      tkconfigure(init,"-from",tclvalue(npopmin),"-to",tclvalue(npopmax))
      if(tclvalue(npopinit)>tclvalue(npopmax)){
        tclvalue(npopinit)<-tclvalue(npopmax)
        tkconfigure(init,"-textvariable",npopinit)
      }
      else if(tclvalue(npopinit)<tclvalue(npopmin)){
        tclvalue(npopinit)<-tclvalue(npopmin)
        tkconfigure(init,"-textvariable",npopinit)
        }
    }
    
    min<-tcl("spinbox",wmin,"-textvariable",npopmin,"-width",5,"-increment",1,"-from",1,"-to",tclvalue(npopmax),"-command",actualizaspin)
    init<-tcl("spinbox", winit,"-textvariable",npopinit,"-width",5,"-increment",1,"-from",tclvalue(npopmin),"-to",tclvalue(npopmax),"-command",actualizaspin)
    max<-tcl("spinbox",wmax,"-textvariable",npopmax,"-width",5,"-increment",1,"-from",tclvalue(npopmin),"-to",100000,"-command",actualizaspin)


    
    npopminlabel.widget <- tklabel(ttrun,text="pop min")
    npopinitlabel.widget <- tklabel(ttrun,text="pop init")
    npopmaxlabel.widget <- tklabel(ttrun,text="pop max")
    
    tkgrid(npopminlabel.widget,row=5,column=1,sticky="e")
    tkgrid(npopinitlabel.widget,row=5,column=2,sticky="w")
    tkgrid(npopmaxlabel.widget,row=5,column=3,sticky="w")
    
    tkgrid(wmin,row=6,column=1,sticky="e")
    tkgrid(winit,row=6,column=2,sticky="w")
    tkgrid(wmax,row=6,column=3,sticky="w")

    # nb.nuclei.max

    nuclei <- tclVar(300)
    nuclei.widget <-tkentry(ttrun,width="20",textvariable=nuclei)
    nucleilabel.widget <- tklabel(ttrun,text="Maximum number of nuclei:")
    
    
    
    # nit

    nit <- tclVar()
    nit.widget <-tkentry(ttrun,width="20",textvariable=nit)
    nitlabel.widget <- tklabel(ttrun,text="Number of iterations:")
    
    tkgrid(nitlabel.widget,row=8,column=1,sticky="w")
    tkgrid(nit.widget,row=8,column=3,sticky="w")
    
    # thinning

    thinning <- tclVar()
    thinning.widget <-tkentry(ttrun,width="20",textvariable=thinning)
    thinninglabel.widget <- tklabel(ttrun,text="Thinning:")
    
    tkgrid(thinninglabel.widget,row=9,column=1,sticky="w")
    tkgrid(thinning.widget,row=9,column=3,sticky="w")

    #freq.model

    freqlabel.widget <- tklabel(ttrun,text="Frequency Model:")
    freq <- tclVar("Dirichlet")
    wfreq <- .Tk.subwin(ttrun)
    freqoptionmenu.widget <-  tcl("tk_optionMenu",wfreq,freq,"Falush","Dirichlet")
    #ploidycombo.widget <-
    #tk("tk_optionMenu","menubuttonoptions","haploid","diploid")


    #varnpop

    varnpoplabel.widget <- tklabel(ttrun,text="Unknown number of populations:")
    varnpop <- tclVar("TRUE")
    wvarnpop <- .Tk.subwin(ttrun)
    varnpopoptionmenu.widget <-  tcl("tk_optionMenu", wvarnpop,varnpop,"FALSE","TRUE")
    #ploidycombo.widget <-
    #tk("tk_optionMenu","menubuttonoptions","haploid","diploid")
    tkgrid(varnpoplabel.widget,row=11,column=1,sticky="w")
    tkgrid(wvarnpop,row=11,column=3,sticky="w")
    
    #spatial

    spatiallabel.widget <- tklabel(ttrun,text="Spatial model:")
    spatial <- tclVar("TRUE")
    wspatial <- .Tk.subwin(ttrun)
    spatialoptionmenu.widget <-  tcl("tk_optionMenu",wspatial,spatial,"FALSE","TRUE")
    #ploidycombo.widget <-
    #tk("tk_optionMenu","menubuttonoptions","haploid","diploid")
    tkgrid(spatiallabel.widget,row=12,column=1,sticky="w")
    tkgrid(wspatial,row=12,column=3,sticky="w")
    
    #jcf

#    jcflabel.widget <- tklabel(ttrun,text="Joint update of colors and frequencies:")
    jcf <- tclVar("FALSE")
#    wjcf <- .Tk.subwin(ttrun)
#    jcfoptionmenu.widget <-  tcl("tk_optionMenu",wjcf,jcf,"FALSE","TRUE")
    #ploidycombo.widget <-
    #tk("tk_optionMenu","menubuttonoptions","haploid","diploid")

    labelspace <-tklabel(ttrun,text=" ")
    tkgrid(labelspace,row=13,column=1)

    #next

    nextbutton <- tkbutton(ttrun,image=imagerun2,text="RUN >>",command=RunmcmcFmodel)
    tkgrid(nextbutton,row=14,column=3,sticky="e")

    tkfocus(ttrun)

#    tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
#    tkgrid(cb1.widget,row=1,column=2,sticky="w")
#    tkgrid(entry1.widget,row=1,column=3,sticky="w")
    tkgrid(ratelabel.widget,row=3,column=1,sticky="w")
    tkgrid(rate.widget,row=3,column=3,sticky="w")
    tkgrid(nucleilabel.widget,row=7,column=1,sticky="w")
    tkgrid(nuclei.widget,row=7,column=3,sticky="w")
 #  tkgrid(jcflabel.widget,row=13,column=1,sticky="w")
 #   tkgrid(wjcf,row=13,column=3,sticky="w")
    tkgrid(freqlabel.widget,row=10,column=1,sticky="w")
    tkgrid(wfreq,row=10,column=3,sticky="w")
    
      
    if(tclvalue(advanced)==1){
      #tkconfigure(alellelabel.widget,state="normal")
      #tkconfigure(cb1.widget,state="normal")
        #tkconfigure(entry1.widget,state="normal")
      tkconfigure(ratelabel.widget,state="normal")
      tkconfigure(rate.widget,state="normal")
      tkconfigure(nucleilabel.widget,state="normal")
      tkconfigure(nuclei.widget,state="normal")
  #    tkconfigure(jcflabel.widget,state="normal")
 #     tkconfigure(wjcf,state="normal")
      tkconfigure(freqlabel.widget,state="normal")
      tkconfigure(wfreq,state="normal")
      
    }
    else{
      
      #tkconfigure(alellelabel.widget,state="disable")
      #tkconfigure(cb1.widget,state="disable")
      #tkconfigure(entry1.widget,state="disable")
      tkconfigure(ratelabel.widget,state="disable")
      tkconfigure(rate.widget,state="disable")
      tkconfigure(nucleilabel.widget,state="disable")
      tkconfigure(nuclei.widget,state="disable")
 #     tkconfigure(jcflabel.widget,state="disable")
 #     tkconfigure(wjcf,state="disable")
      tkconfigure(freqlabel.widget,state="disable")
      tkconfigure(wfreq,state="disable")
      
      
    }
  }
  
  initialimage <-function(){
    imgAsLabel <- tklabel(ttinit,image=image1,bg="white")
    tkgrid(imgAsLabel,sticky="news")
    }

     # -----------------------------------------------
     # HELLO | POSTPROCESSING STARTS HERE
     #------------------------------------

  postproc <-function(){
    
    RunpostprocessChain <- function(){
      
      tttry <- tktoplevel(parent=.TkRoot)
      tkgrab(tttry)
      tkwm.geometry(tttry, "+200+200")
      tkwm.title(tttry,"wait")
      warn<-tklabel(tttry,image=imagepleasewait)
      tkpack(warn)
      
      
      tkfocus(tttry)
      
      Sys.sleep(0.1)

      
      print("Starting...")
     
      
      err <- try(PostProcessChain(coordinates=globalcoordinates,genotypes=globalgenotypes,path.mcmc=tclvalue(outputdir),nxdom=as.numeric(tclvalue(nxdom)),nydom=as.numeric(tclvalue(nydom)),burnin=as.numeric(tclvalue(burnin))),silent=TRUE)
      tkdestroy(tttry)
      
      print("Done.")
      if (class(err) == "try-error") 
        tkmessageBox(message=err,icon="error",type="ok",parent=tt)               
      else
        tkmessageBox(message="Terminated with success",type="ok",parent=tt)
    }

    ndomlabel.widget <- tklabel(ttpost,text="Number of pixels in the spatial domain:")
    tkgrid(ndomlabel.widget,row=1,column=1,columnspan=2,sticky="w")

      #nxdom  and nydom
    nxdom <- tclVar(50)
    nydom <- tclVar(50)
    nxdomlabel.widget <- tklabel(ttpost,text="X")
    nydomlabel.widget <- tklabel(ttpost,text="Y")
    tkgrid(nxdomlabel.widget,row=2,column=1,sticky="w")
    tkgrid(nydomlabel.widget,row=2,column=2,sticky="w")
    
    nxdom.widget <-tkentry(ttpost,width="20",textvariable=nxdom)
    nydom.widget <-tkentry(ttpost,width="20",textvariable=nydom)
    
    tkgrid(nxdom.widget,row=3,column=1,sticky="w")
    tkgrid(nydom.widget,row=3,column=2,sticky="w")

      #burnin
    
    burnin <- tclVar(0)
    burnin.widget <-tkentry(ttpost,width="20",textvariable=burnin)
    burninlabel.widget <- tklabel(ttpost,text="Burnin:")
    
    tkgrid(burninlabel.widget,row=4,column=1,sticky="w")
    tkgrid(burnin.widget,row=4,column=2,sticky="w")

    labelspace <-tklabel(ttpost,text=" ")
    tkgrid(labelspace,row=5,column=1)
    
    nextbutton <- tkbutton(ttpost,image=imagerun2,text="RUN >>",command=RunpostprocessChain)
    tkgrid(nextbutton,row=6,column=2,sticky="e")
    
    tkfocus(ttpost)

  }

     # -----------------------------------------------
     # HELLO | ShowIBD STARTS HERE
     #------------------------------------

  
  GraficalIBD <- function(){

    if (length(idb.dataset)==1){
      
      tkmessageBox(message="First simulate some data",icon="error",type="ok",parent=tt)
      
    }
    else{

      DrawShowIBD <- function(){

        vect1<- c()
        #vect2 <- c()
        
        vec1<-unlist(strsplit(tclvalue(loc.grid),","))
        for (i in 1:length(vec1)) vect1[i]<-as.numeric(vec1[i])
        #vec2<-unlist(strsplit(tclvalue(all.grid),","))
        #for (i in 1:length(vec2)) vect2[[i]]<-as.numeric(vec2[i])

        if(tclvalue(plot.coord.ps) == "")
          file.plot.coord<-NA
        else
          file.plot.coord<-tclvalue(plot.coord.ps)
        
        if(tclvalue(plot.tess.ps) == "")
          file.plot.tess<-NA
        else
          file.plot.tess<-tclvalue(plot.tess.ps)
        
                    
        if(tclvalue(plot.freq.grid.ps) == "")
          file.plot.freq.grid<-NA
        else
         file.plot.freq.grid<-tclvalue(plot.freq.grid.ps)
        
          
        if(tclvalue(plot.freq.indiv.ps) == "")
         file.plot.freq.indiv<-NA
        else
          file.plot.freq.indiv<-tclvalue(plot.freq.indiv.ps)
        
        if(tclvalue(plot.gen.ps) == "")
          file.plot.gen <- NA
        else
          file.plot.gen<-tclvalue(plot.gen.ps)
        
        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        
      
        tkfocus(tttry)
        
        Sys.sleep(0.1)

      
        print("Starting...")
        err <- try(show.simdata(dataset=idb.dataset,plot.coord=as.logical(tclvalue(plot.coord)),file.plot.coord=file.plot.coord,plot.tess=as.logical(tclvalue(plot.tess)),file.plot.tess=file.plot.tess,plot.freq.grid=as.logical(tclvalue(plot.freq.grid)),file.plot.freq.grid=file.plot.freq.grid,plot.freq.indiv=as.logical(tclvalue(plot.freq.indiv)),file.plot.freq.indiv=file.plot.freq.indiv,loc.grid=vect1,loc.indiv=as.numeric(tclvalue(loc.indiv)),zlim.freq=c(as.numeric(tclvalue(zlimin)),as.numeric(tclvalue(zlimax))),plot.gen=as.logical(tclvalue(plot.gen)),file.plot.gen=file.plot.gen),silent=TRUE)
        tkdestroy(tttry)
        print("Done")

        if (class(err) == "try-error") 
          tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
        else
          tkmessageBox(message="Terminated with success",type="ok",parent=tt)
      }
      
      ttshowibd <- tktoplevel()
      tkwm.title(ttshowibd,"Graphical Display of data simulated by IBD")

#      printit <- tclVar(0)
      plot.coord.ps  <- tclVar("")
      plot.tess.ps<- tclVar("")
      plot.freq.grid.ps<- tclVar("")
      plot.freq.indiv.ps<- tclVar("")
      plot.gen.ps<- tclVar("")
      

      #setprint <- function()  {
      #  
      #  if(tclvalue(printit)==1){
      #    tkconfigure(printbutton.widget,state="normal")
      #  }else{
      #    tkconfigure(printbutton.widget,state="disable")
      #    tclvalue(printfile) <- ""
      #  }
      #}
      
      getprintfile <- function()  {
        printfile <- tclVar()
        tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes="{{.ps} *.ps}"))
        tkfocus(ttshowibd)
        return(printfile)
      }
      
      #alellelabel.widget <- tklabel(ttshowibd,text="Save to directory?")
      #cb2.widget <- tkcheckbutton(ttshowibd,command=setprint,variable=printit,onvalue=1, offvalue=0)
      plot.coord.psbutton.widget <- tkbutton(ttshowibd,text="Save File",command=function(){plot.coord.ps<<-getprintfile();tkconfigure( plot.coord.pslabel.widget,text=tclvalue(plot.coord.ps))},width=15)
      plot.coord.pslabel.widget <- tklabel(ttshowibd,text=tclvalue(plot.coord.ps),width=50)
      
      plot.tess.psbutton.widget <- tkbutton(ttshowibd,text="Save File",command=function(){plot.tess.ps<<-getprintfile();tkconfigure( plot.tess.pslabel.widget,text=tclvalue(plot.tess.ps))},width=15)
      plot.tess.pslabel.widget <- tklabel(ttshowibd,text=tclvalue(plot.tess.ps),width=50)

      plot.freq.grid.psbutton.widget <- tkbutton(ttshowibd,text="Save File",command=function(){plot.freq.grid.ps<<-getprintfile();tkconfigure( plot.freq.grid.pslabel.widget,text=tclvalue(plot.freq.grid.ps))},width=15)
      plot.freq.grid.pslabel.widget <- tklabel(ttshowibd,text=tclvalue(plot.freq.grid.ps),width=50)

      plot.freq.indiv.psbutton.widget <- tkbutton(ttshowibd,text="Save File",command=function(){plot.freq.indiv.ps<<-getprintfile();tkconfigure(plot.freq.indiv.pslabel.widget ,text=tclvalue(plot.freq.indiv.ps))},width=15)
      plot.freq.indiv.pslabel.widget <- tklabel(ttshowibd,text=tclvalue(plot.freq.indiv.ps),width=50)
      
      plot.gen.psbutton.widget <- tkbutton(ttshowibd,text="Save File",command=function(){plot.gen.ps<<-getprintfile();tkconfigure(plot.gen.pslabel.widget,text=tclvalue(plot.gen.ps))},width=15)
      plot.gen.pslabel.widget <- tklabel(ttshowibd,text=tclvalue(plot.gen.ps),width=50)
      
      
       #plot.coord
      plot.coordlabel.widget <- tklabel(ttshowibd,text="Plot coordinates of individuals:")
      plot.coord <- tclVar("FALSE")
      wplot.coord <- .Tk.subwin(ttshowibd)
      plot.coordoptionmenu.widget <-  tcl("tk_optionMenu",wplot.coord,plot.coord,"FALSE","TRUE")
      
       #plot.tess
      plot.tesslabel.widget <- tklabel(ttshowibd,text="Plot tessellation:")
      plot.tess <- tclVar("FALSE")
      wplot.tess <- .Tk.subwin(ttshowibd)
      plot.tessoptionmenu.widget <-  tcl("tk_optionMenu",wplot.tess,plot.tess,"FALSE","TRUE")
    
       #plot.freq.grid
      plot.freq.gridlabel.widget <- tklabel(ttshowibd,text="Plot allele frequencies for all pixels:")
      plot.freq.grid <- tclVar("FALSE")
      wplot.freq.grid <- .Tk.subwin(ttshowibd)
      plot.freq.gridoptionmenu.widget <-  tcl("tk_optionMenu",wplot.freq.grid,plot.freq.grid,"FALSE","TRUE")
    
       #plot.freq.indiv
      plot.freq.indivlabel.widget <- tklabel(ttshowibd,text="Plot allele frequencies at individual sites:")
      plot.freq.indiv <- tclVar("FALSE")
      wplot.freq.indiv <- .Tk.subwin(ttshowibd)
      plot.freq.indivoptionmenu.widget <-  tcl("tk_optionMenu", wplot.freq.indiv,plot.freq.indiv,"FALSE","TRUE")
    

       #loc.grid
      loc.grid <- tclVar(1)
      loc.grid.widget <-tkentry(ttshowibd,width="20",textvariable=loc.grid)
      loc.gridlabel.widget <- tklabel(ttshowibd,text="Indices of loci (E.g: 1,2,4,...):")

       #all.grid
#      all.grid <- tclVar(0)
#      all.grid.widget <-tkentry(ttshowibd,width="20",textvariable=all.grid)
#      all.gridlabel.widget <- tklabel(ttshowibd,text="Indices of alleles (E.g: n,m,..):")

      #loc.indiv
      loc.indiv <- tclVar(1)
      loc.indiv.widget <-tkentry(ttshowibd,width="20",textvariable=loc.indiv)
      loc.indivlabel.widget <- tklabel(ttshowibd,text="Indices of loci (E.g: 1,3,6,...):")
      
      #zlim.freq=c(0,1)

      zlimin <- tclVar(0)
      zlimax <- tclVar(1)
      #zlimin.widget <-tkentry(ttshowibd,width="10",textvariable=zlimin,state="disabled")
      #zlimax.widget <-tkentry(ttshowibd,width="10",textvariable=zlimax,state="disabled")
      #zlimlabel.widget <- tklabel(ttshowibd,text="Z Limit:",state="disabled")
      
      
      #plot.gen
      plot.genlabel.widget <- tklabel(ttshowibd,text="Plot genotypes:")
      plot.gen <- tclVar("FALSE")
      wplot.gen <- .Tk.subwin(ttshowibd)
      plot.genoptionmenu.widget <-  tcl("tk_optionMenu", wplot.gen,plot.gen,"FALSE","TRUE")

#      tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
#      tkgrid(cb2.widget,row=1,column=2,sticky="w")
      
      tkgrid(plot.coordlabel.widget,row=2,column=1,sticky="w")
      tkgrid(wplot.coord,row=2,column=2,sticky="w")
      tkgrid(plot.coord.psbutton.widget,row=2,column=3,sticky="e")
      tkgrid(plot.coord.pslabel.widget,row=2,column=4,sticky="w")
      tkgrid(plot.tesslabel.widget,row=3,column=1,sticky="w")
      tkgrid(wplot.tess,row=3,column=2,columnspan=2,sticky="w")
      tkgrid(plot.tess.psbutton.widget,row=3,column=3,sticky="e")
      tkgrid(plot.tess.pslabel.widget,row=3,column=4,sticky="w")
      tkgrid(plot.freq.gridlabel.widget,row=4,column=1,sticky="w")
      tkgrid(wplot.freq.grid,row=4,column=2,columnspan=2,sticky="w")
      tkgrid(plot.freq.grid.psbutton.widget,row=4,column=3,sticky="e")
      tkgrid(plot.freq.grid.pslabel.widget,row=4,column=4,sticky="w")
      tkgrid(plot.freq.indivlabel.widget,row=6,column=1,sticky="w")
      tkgrid(wplot.freq.indiv,row=6,column=2,columnspan=2,sticky="w")
      tkgrid(plot.freq.indiv.psbutton.widget,row=6,column=3,sticky="e")
      tkgrid(plot.freq.indiv.pslabel.widget,row=6,column=4,sticky="w")
      tkgrid(loc.gridlabel.widget,row=5,column=1,sticky="w")
      tkgrid(loc.grid.widget,row=5,column=2,columnspan=2,sticky="w")
#      tkgrid(all.gridlabel.widget,row=6,column=1,sticky="w")
#      tkgrid(all.grid.widget,row=6,column=2,sticky="w")
      tkgrid(loc.indivlabel.widget,row=7,column=1,sticky="w")
      tkgrid(loc.indiv.widget,row=7,column=2,columnspan=2,sticky="w")
      #tkgrid(zlimlabel.widget,row=7,column=1,sticky="w")
      #tkgrid(zlimin.widget,row=7,column=2,sticky="w")
      #tkgrid(zlimax.widget,row=7,column=3,sticky="w")
      tkgrid(plot.genlabel.widget,row=8,column=1,sticky="w")
      tkgrid(wplot.gen,row=8,column=2,columnspan=2,sticky="w")
      tkgrid(plot.gen.psbutton.widget,row=8,column=3,sticky="e")
      tkgrid(plot.gen.pslabel.widget,row=8,column=4,sticky="w")
      labelspace <-tklabel(ttshowibd,text=" ")
      tkgrid(labelspace,row=9,column=1) 
      
      nextbutton <- tkbutton(ttshowibd,image=imagedraw,text="Draw >>",command=DrawShowIBD)
      tkgrid(nextbutton,row=10,column=2,columnspan=2,sticky="e")
      
      tkfocus(ttshowibd)
    }
  }
  
   
     # -----------------------------------------------
     # HELLO | PLOT STARTS HERE
     #------------------------------------

  
  plot <- function(){
      
      
      # -----------------------------------------------
      # DRIFT
      #------------------------------------

    Drift <- function(){
      
      printit <- tclVar(0)
      printfile <- tclVar("")
      
      DrawDrift <- function(){
        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        tkfocus(tttry)
        
        print("Starting...")
        
        Sys.sleep(0.1)
        
        if(tclvalue(printit)==1){
          err <- try(PlotDrift(path.mcmc=tclvalue(outputdir),printit=TRUE,file=tclvalue(printfile)),silent=TRUE)
        }else{
          err <- try(PlotDrift(path.mcmc=tclvalue(outputdir),printit=FALSE),silent=TRUE)
        }
        tkdestroy(tttry)
        print("Done.")
        if (class(err) == "try-error") 
          tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
        else
          tkmessageBox(message="Terminated with success",type="ok",parent=tt)
      }

      
      #Drift buttons
      ttdrift <- tktoplevel()
      
      tkwm.title(ttdrift,"Drift factors")
      
      setprint <- function()  {
        
        if(tclvalue(printit)==1){
          tkconfigure(printbutton.widget,state="normal")
        }else{
          tkconfigure(printbutton.widget,state="disable")
          tclvalue(printfile) <- ""
        }
      }
      
      getprintfile <- function()  {
        tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes="{{.ps} *.ps}"))
        tkfocus(ttdrift)
      }
      
      alellelabel.widget <- tklabel(ttdrift,text="Save to file?")
      cb2.widget <- tkcheckbutton(ttdrift,command=setprint,variable=printit,onvalue=1, offvalue=0)
      printbutton.widget <- tkbutton(ttdrift,text="Save File",command=getprintfile,width=15)
      filelabel.widget <- tklabel(ttdrift,textvariable=printfile,width=50)
      
      
      tkdeselect(cb2.widget)
      tkconfigure(printbutton.widget,state="disable")
      
      tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
      tkgrid(cb2.widget,row=1,column=2,sticky="w")
      tkgrid(printbutton.widget,row=1,column=3,sticky="w")
      tkgrid(filelabel.widget,row=1,column=4,sticky="w")
      
      labelspace <-tklabel(ttdrift,text=" ")
      tkgrid(labelspace,row=2,column=1)
      
      nextbutton <- tkbutton(ttdrift,image=imagedraw,text="Draw >>",command=DrawDrift)
      tkgrid(nextbutton,row=3,column=2,sticky="e")
      
      tkfocus(ttdrift)
      
      
    }

      # -----------------------------------------------
      # FREQ
      #------------------------------------

    Freq <- function(){
      
      printit <- tclVar(0)
      printfile <- tclVar("")
      
      ipop <- tclVar(1)
      iloc <- tclVar(1)
      iall <- tclVar(1)

      strallele <- paste(globalallele.numbers, sep = ",", collapse =",")#globalallele.numbers[0]
      
      
      allelenumbers<-tclVar(strallele)
      
      Drawfreq <- function(){
        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        tkfocus(tttry)
        
        vect1<- c()
        vec1<-unlist(strsplit(tclvalue(allelenumbers),","))
        for (i in 1:length(vec1)) vect1[i]<-as.numeric(vec1[i])
        
          
        print("Starting...")
        
        
        Sys.sleep(0.1)
        
        if(tclvalue(printit)==1){
          
            err <- try(PlotFreq(allele.numbers=vect1,path.mcmc=tclvalue(outputdir),ipop=as.numeric(tclvalue(ipop)),iloc=as.numeric(tclvalue(iloc)),iall=as.numeric(tclvalue(iall)),printit=TRUE,path=tclvalue(printfile)),silent=TRUE)
          }else{
            err <- try(PlotFreq(allele.numbers=vect1,path.mcmc=tclvalue(outputdir),ipop=as.numeric(tclvalue(ipop)),iloc=as.numeric(tclvalue(iloc)),iall=as.numeric(tclvalue(iall)),printit=FALSE),silent=TRUE)
          }
        tkdestroy(tttry)
        print("Done.")
        if (class(err) == "try-error") 
          tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
        else
          tkmessageBox(message="Terminated with success",type="ok",parent=tt)
      }
      
      
      
      ttfreq <- tktoplevel()
      tkwm.title(ttfreq,"Frequencies in population")
      


        setprint <- function()  {

          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }

        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkchooseDirectory(title="Please choose an output directory"))
          tcl("regsub","-all","\\\\",tclvalue(printfile),"/",printfile)
          tcl("append",printfile,"/")
          tkfocus(ttfreq)
        }

        alellelabel.widget <- tklabel(ttfreq,text="Save to file?")
        cb2.widget <- tkcheckbutton(ttfreq,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(ttfreq,text="Save Directory",command=getprintfile,width=15)
        filelabel.widget <- tklabel(ttfreq,textvariable=printfile,width=50)


        tkdeselect(cb2.widget)
        tkconfigure(printbutton.widget,state="disable")

        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")


        #allelenumbers
        #allelenumbers.widget <-tkentry(ttfreq,width="20",textvariable=allelenumbers)
        #allelenumberslabel.widget <- tklabel(ttfreq,text="Number of alleles per locus (E.g: n,m,..):")

        #tkgrid(allelenumberslabel.widget,row=2,column=1,sticky="w")
        #tkgrid(allelenumbers.widget,row=2,column=2,sticky="w")

        #ipop
        ipop.widget <-tkentry(ttfreq,width="20",textvariable=ipop)
        ipoplabel.widget <- tklabel(ttfreq,text="Index of populations:")

        tkgrid(ipoplabel.widget,row=3,column=1,sticky="w")
        tkgrid(ipop.widget,row=3,column=2,sticky="w")
        #iloc
        iloc.widget <-tkentry(ttfreq,width="20",textvariable=iloc)
        iloclabel.widget <- tklabel(ttfreq,text="Index of locus:")

        tkgrid(iloclabel.widget,row=4,column=1,sticky="w")
        tkgrid(iloc.widget,row=4,column=2,sticky="w")

        #iall
        iall.widget <-tkentry(ttfreq,width="20",textvariable=iall)
        ialllabel.widget <- tklabel(ttfreq,text="Index of allele:")

        tkgrid(ialllabel.widget,row=5,column=1,sticky="w")
        tkgrid(iall.widget,row=5,column=2,sticky="w")

        labelspace <-tklabel(ttfreq,text=" ")
        tkgrid(labelspace,row=6,column=1)
        

        nextbutton <- tkbutton(ttfreq,image=imagedraw,text="Draw >>",command=Drawfreq)
        tkgrid(nextbutton,row=7,column=2,sticky="e")

        tkfocus(ttfreq)


      }




      # -----------------------------------------------
      # FREQA
      #------------------------------------

      FreqA <- function(){

        printit <- tclVar(0)
        printfile <- tclVar("")

        iloc <- tclVar(1)
        iall <- tclVar(1)
        #strallele <- paste(globalallele.numbers, sep = ",", collapse =",")#globalallele.numbers[0]

        #allelenumbers<-tclVar(strallele)

        Drawfreqa <- function(){

          #vect1<- c()
          #vec1<-unlist(strsplit(tclvalue(allelenumbers),","))
          #for (i in 1:length(vec1)) vect1[i]<-as.numeric(vec1[i])
          
          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)
          print("Starting...")
       

          
          Sys.sleep(0.1)
          
          if(tclvalue(printit)==1){
            err <- try(PlotFreqA(genotypes=globalgenotypes,path.mcmc=tclvalue(outputdir),iloc=as.numeric(tclvalue(iloc)),iall=as.numeric(tclvalue(iall)),printit=TRUE,path=tclvalue(printfile)),silent=TRUE)
          }else{
            err <- try(PlotFreqA(allele.numbers=vect1,path.mcmc=tclvalue(outputdir),iloc=as.numeric(tclvalue(iloc)),iall=as.numeric(tclvalue(iall)),printit=FALSE),silent=TRUE)
          }
          tkdestroy(tttry)
          print("Done.")
          if (class(err) == "try-error") 
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
          else
            tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        }



        ttfreqa <- tktoplevel()
        tkwm.title(ttfreqa,"Frequencies in ancestral population")



        setprint <- function()  {

          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }

        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkchooseDirectory(title="Please choose an output directory"))
          tcl("regsub","-all","\\\\",tclvalue(printfile),"/",printfile)
          tcl("append",printfile,"/")

          tkfocus(ttfreqa)
        }

        alellelabel.widget <- tklabel(ttfreqa,text="Save to file?")
        cb2.widget <- tkcheckbutton(ttfreqa,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(ttfreqa,text="Save Directory",command=getprintfile,width=15)
        filelabel.widget <- tklabel(ttfreqa,textvariable=printfile,width=50)


        tkdeselect(cb2.widget)
        tkconfigure(printbutton.widget,state="disable")

        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")


        #allelenumbers
        #allelenumbers.widget <-tkentry(ttfreqa,width="20",textvariable=allelenumbers)
        #allelenumberslabel.widget <- tklabel(ttfreqa,text="Number of alleles per locus (E.g: n,m,..):")

        #tkgrid(allelenumberslabel.widget,row=2,column=1,sticky="w")
        #tkgrid(allelenumbers.widget,row=2,column=2,sticky="w")

        #iloc
        iloc.widget <-tkentry(ttfreqa,width="20",textvariable=iloc)
        iloclabel.widget <- tklabel(ttfreqa,text="Index of locus:")

        tkgrid(iloclabel.widget,row=3,column=1,sticky="w")
        tkgrid(iloc.widget,row=3,column=2,sticky="w")

        #iall
        iall.widget <-tkentry(ttfreqa,width="20",textvariable=iall)
        ialllabel.widget <- tklabel(ttfreqa,text="Index of allele:")

        tkgrid(ialllabel.widget,row=4,column=1,sticky="w")
        tkgrid(iall.widget,row=4,column=2,sticky="w")

        labelspace <-tklabel(ttfreqa,text=" ")
        tkgrid(labelspace,row=5,column=1)

        nextbutton <- tkbutton(ttfreqa,image=imagedraw,text="Draw >>",command=Drawfreqa)
        tkgrid(nextbutton,row=6,column=2,sticky="e")

        tkfocus(ttfreqa)


      }

      # -----------------------------------------------
      # TESSELLATION
      #------------------------------------



      Tessellation <- function(){

        printit <- tclVar(0)
        printfile <- tclVar("")


        Drawtessellation <- function(){

    
          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)
          
          print("Starting...")
          
          Sys.sleep(0.1)
          
          if(tclvalue(printit)==1){
            err <- try(PlotTessellation(coordinates=globalcoordinates,path.mcmc=tclvalue(outputdir),printit=TRUE,path=tclvalue(printfile)),silent=TRUE)
          }else{
            err <- try(PlotTessellation(coordinates=globalcoordinates,path.mcmc=tclvalue(outputdir),printit=FALSE),silent=TRUE)
          }
          tkdestroy(tttry)
          print("Done.")
          if (class(err) == "try-error")
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)
          else
            tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        }


        tttessellation <- tktoplevel()
        tkwm.title(tttessellation,"Tessellation")

        printit <- tclVar(0)
        printfile <- tclVar("")

        setprint <- function()  {

          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }

        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkchooseDirectory(title="Please choose an output directory"))
          tcl("regsub","-all","\\\\",tclvalue(printfile),"/",printfile)
          tcl("append",printfile,"/")
          tkfocus(tttessellation)
        }

        alellelabel.widget <- tklabel(tttessellation,text="Save to file?")
        cb2.widget <- tkcheckbutton(tttessellation,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(tttessellation,text="Save Directory",command=getprintfile,width=15)
        filelabel.widget <- tklabel(tttessellation,textvariable=printfile,width=50)


        tkdeselect(cb2.widget)
        tkconfigure(printbutton.widget,state="disable")

        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")

        labelspace <-tklabel(tttessellation,text=" ")
        tkgrid(labelspace,row=2,column=1)

        nextbutton <- tkbutton(tttessellation,image=imagedraw,text="Draw >>",command=Drawtessellation)
        tkgrid(nextbutton,row=3,column=2,sticky="e")

        tkfocus(tttessellation)


      }





      # -----------------------------------------------
      # NPOP
      #------------------------------------
      Npop <- function(){


        printit <- tclVar(0)
        printfile <- tclVar("")

        Drawnpop <- function(){
          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)
          print("Starting...")
          
          Sys.sleep(0.1)
                    
          if(tclvalue(printit)==1){
            err <- try(Plotnpop(path.mcmc=tclvalue(outputdir),printit=TRUE,file=tclvalue(printfile)),silent=TRUE)
          }else{
            err <- try(Plotnpop(path.mcmc=tclvalue(outputdir),printit=FALSE),silent=TRUE)
          }
          tkdestroy(tttry)
          print("Done.")
          if (class(err) == "try-error") 
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
          else
            tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        }
        
       

        ttnpop <- tktoplevel()
        tkwm.title(ttnpop,"Number of populations")


        printit <- tclVar(0)
        printfile <- tclVar("")

        setprint <- function()  {

          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }

        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes="{{.ps} *.ps}"))
          tkfocus(ttnpop)
        }

        alellelabel.widget <- tklabel(ttnpop,text="Save to file?")
        cb2.widget <- tkcheckbutton(ttnpop,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(ttnpop,text="Save File",command=getprintfile,width=15)
        filelabel.widget <- tklabel(ttnpop,textvariable=printfile,width=50)


        tkdeselect(cb2.widget)
        tkconfigure(printbutton.widget,state="disable")

        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")

        labelspace <-tklabel(ttnpop,text=" ")
        tkgrid(labelspace,row=2,column=1)

        nextbutton <- tkbutton(ttnpop,image=imagedraw,text="Draw >>",command=Drawnpop)
        tkgrid(nextbutton,row=3,column=2,sticky="e")

        tkfocus(ttnpop)


      }

      # -----------------------------------------------
      # NTILE
      #------------------------------------


      
      Ntile <- function(){


        printit <- tclVar(0)
        printfile <- tclVar("")
        
        Drawntile <- function(){
          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)

          print("Starting...")
          Sys.sleep(0.1)
                    
          if(tclvalue(printit)==1){
            err <- try(Plotntile(path.mcmc=tclvalue(outputdir),printit=TRUE,file=tclvalue(printfile)),silent=TRUE)
          }else{
            err <- try(Plotntile(path.mcmc=tclvalue(outputdir),printit=FALSE),silent=TRUE)
          }
          tkdestroy(tttry)
          print("Done.")
          if (class(err) == "try-error") 
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
          else
            tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        }
        

        
        ttntile <- tktoplevel()
        tkwm.title(ttntile,"Number of tiles")
        
        
        printit <- tclVar(0)
        printfile <- tclVar("")
        
        setprint <- function()  {

          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }

        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes="{{.ps} *.ps}"))
          tkfocus(ttntile)
        }

        alellelabel.widget <- tklabel(ttntile,text="Save to file?")
        cb2.widget <- tkcheckbutton(ttntile,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(ttntile,text="Save File",command=getprintfile,width=15)
        filelabel.widget <- tklabel(ttntile,textvariable=printfile,width=50)


        tkdeselect(cb2.widget)
        tkconfigure(printbutton.widget,state="disable")

        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")

        labelspace <-tklabel(ttntile,text=" ")
        tkgrid(labelspace,row=2,column=1)
        
        nextbutton <- tkbutton(ttntile,image=imagedraw,text="Draw >>",command=Drawntile)
        tkgrid(nextbutton,row=3,column=2,sticky="e")

        tkfocus(ttntile)


      }
      
      # -----------------------------------------------
      # POSTERIOR MODE
      #------------------------------------
      
      
      PosteriorM <- function(){
        

        ttposm <- tktoplevel()
        tkwm.title(ttposm,"Map of populations")

        
        
        write <- tclVar("FALSE")
        plotit <- tclVar("TRUE")
        maintitle <- tclVar("")


        Drawposm <- function(){
          

          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)
          
          print("Starting...")

          Sys.sleep(0.1)
          
          
          if(tclvalue(printit)==1){
            err <- try(PosteriorMode(coordinates=globalcoordinates,path.mcmc=tclvalue(outputdir),write=as.logical(tclvalue(write)),plotit=as.logical(tclvalue(plotit)),printit=TRUE,file=tclvalue(printfile),main.title=tclvalue(maintitle)),silent=TRUE)
          }else{
            err <- try(PosteriorMode(coordinates=globalcoordinates,path.mcmc=tclvalue(outputdir),write=as.logical(tclvalue(write)),plotit=as.logical(tclvalue(plotit)),printit=FALSE,file="",main.title=as.character(tclvalue(maintitle))),silent=TRUE)
          }
          tkdestroy(tttry)
          print("Done.")
          if (class(err) == "try-error") 
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
          else
            tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        }
        
        
        printit <- tclVar(0)
        printfile <- tclVar("")
        
        
        setprint <- function()  {
          
          if(tclvalue(printit)==1){
            tkconfigure(printbutton.widget,state="normal")
          }else{
            tkconfigure(printbutton.widget,state="disable")
            tclvalue(printfile) <- ""
          }
        }
        
        getprintfile <- function()  {
          tclvalue(printfile) <- tclvalue(tkgetSaveFile(filetypes="{{.ps} *.ps}"))
          tkfocus(ttposm)
        }
        
        alellelabel.widget <- tklabel(ttposm,text="Save to file?")
        cb2.widget <- tkcheckbutton(ttposm,command=setprint,variable=printit,onvalue=1, offvalue=0)
        printbutton.widget <- tkbutton(ttposm,text="Save File",command=getprintfile,width=15,state="disable")
        filelabel.widget <- tklabel(ttposm,textvariable=printfile,width=50)
        
        writelabel.widget <- tklabel(ttposm,text="Write in ASCII(readable) format:")
        wwrite <- .Tk.subwin(ttposm)
        writeoptionmenu.widget <-  tcl("tk_optionMenu",wwrite,write,"FALSE","TRUE")

        plotitlabel.widget <- tklabel(ttposm,text="Plot the map:")
        wplotit <- .Tk.subwin(ttposm)
        plotitoptionmenu.widget <-  tcl("tk_optionMenu",wplotit,plotit,"FALSE","TRUE")


        maintitle.widget <-tkentry(ttposm,width="20",textvariable=maintitle)
        maintitlelabel.widget <- tklabel(ttposm,text="Graph Title:")

        labelspace <-tklabel(ttposm,text=" ")
        tkgrid(labelspace,row=6,column=1)
        
        tkgrid(alellelabel.widget,row=1,column=1,sticky="w")
        tkgrid(cb2.widget,row=1,column=2,sticky="w")
        tkgrid(printbutton.widget,row=1,column=3,sticky="w")
        tkgrid(filelabel.widget,row=1,column=4,sticky="w")
        tkgrid(maintitlelabel.widget,row=2,column=1,sticky="w")
        tkgrid(maintitle.widget,row=2,column=2,sticky="w")
        tkgrid(writelabel.widget,row=4,column=1,sticky="w")
        tkgrid(wwrite,row=4,column=2,sticky="w")
        tkgrid(plotitlabel.widget,row=5,column=1,sticky="w")
        tkgrid(wplotit,row=5,column=2,sticky="w")

        nextbutton <- tkbutton(ttposm,image=imagedraw,text="Draw >>",command=Drawposm)
        tkgrid(nextbutton,row=7,column=2,sticky="e")

        tkfocus(ttposm)

      }

    
    
    
    
    if(falush==0){
      buttondrift <- tkbutton(ttplot,width=30,text="Drift factors",command=Drift,state="disabled")
      buttonfreqA <- tkbutton(ttplot,width=30,text="Frequencies in ancestral population",command=FreqA,state="disabled")
    }
    else {
  buttondrift <- tkbutton(ttplot,width=30,text="Drift factors",command=Drift)
      buttonfreqA <- tkbutton(ttplot,width=30,text="Frequencies in ancestral population",command=FreqA)
    }
    buttonfreq <- tkbutton(ttplot,width=30,text="Frequencies in population",command=Freq)
   
    buttontessellation <- tkbutton(ttplot,width=30,text="Map of proba. of pop. membership",command=Tessellation)
    buttonnpop <- tkbutton(ttplot,width=30,text="Number of populations",command=Npop)
    buttonntile <- tkbutton(ttplot,width=30,text="Number of tiles",command=Ntile)
    buttonposm <- tkbutton(ttplot,width=30,text="Map of population membership",command=PosteriorM)

   
    
    labelfigures <-tklabel(ttplot,text="-Graphics-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    labeltables <-tklabel(ttplot,text="-Tables-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    labelspace1 <-tklabel(ttplot,text=" ")
    labelspace2 <-tklabel(ttplot,text=" ")
    labelinfo <-tklabel(ttplot,text=paste(tclvalue(outputdir),sep=""))
    
    
    proba.pop.membership <- tklabel(ttplot,text="Posterior probability of population membership for each pixel")
    tkbind(proba.pop.membership,"<Any-Enter>",function() { tkconfigure(proba.pop.membership,foreground="blue")} )
    tkbind(proba.pop.membership,"<Any-Leave>",function() { tkconfigure(proba.pop.membership,foreground="black")} )
    tkbind(proba.pop.membership, "<Button-1>",function() { Showtext("proba.pop.membership.txt")} )

    proba.pop.membership.ind <- tklabel(ttplot,text="Posterior probability of population membership for each individual")
    tkbind(proba.pop.membership.ind,"<Any-Enter>",function() { tkconfigure(proba.pop.membership.ind,foreground="blue")} )
    tkbind(proba.pop.membership.ind,"<Any-Leave>",function() { tkconfigure(proba.pop.membership.ind,foreground="black")} )
    tkbind(proba.pop.membership.ind, "<Button-1>",function() { Showtext("proba.pop.membership.indiv.txt")} )
    
    modal.pop.ind <- tklabel(ttplot,text="Label of modal population for each individual")
    tkbind(modal.pop.ind,"<Any-Enter>",function() { tkconfigure(modal.pop.ind,foreground="blue")} )
    tkbind(modal.pop.ind,"<Any-Leave>",function() { tkconfigure(modal.pop.ind,foreground="black")} )
    tkbind(modal.pop.ind, "<Button-1>",function() { Showtext("modal.pop.indiv.txt")} )
   
    
    labelfstat <-tklabel(ttplot,text="-F statistics-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    buttonfstat <- tkbutton(ttplot,width=20,text="Fst and Fis",font="*-Helvetica-bold-i-*",command=Gfstat)
    
    tkgrid(labelfigures,row=1,column=1,sticky="w")
    tkgrid(buttondrift,row=5,column=2,sticky="w")
    tkgrid(buttonfreq,row=3,column=2,sticky="w")
    tkgrid(buttonfreqA,row=4,column=2,sticky="w")
    tkgrid(buttontessellation,row=4,column=1,sticky="w")
    tkgrid(buttonnpop,row=2,column=1,sticky="w")
    tkgrid(buttonntile,row=2,column=2,sticky="w")
    tkgrid(buttonposm,row=3,column=1,sticky="w")

    tkgrid(labelspace1,row=6,column=1,sticky="w")
    tkgrid(labeltables,row=7,column=1,sticky="w")
    tkgrid(proba.pop.membership,row=8,column=1,columnspan=2,sticky="w")
    tkgrid(proba.pop.membership.ind,row=9,column=1,columnspan=2,sticky="w")
    tkgrid(modal.pop.ind,row=10,column=1,columnspan=2,sticky="w")

    tkgrid(labelspace2,row=11,column=1,sticky="w")
    tkgrid(labelfstat,row=12,column=1,sticky="w")
    tkgrid(buttonfstat,row=13,column=1,sticky="w")
  }


       # -----------------------------------------------
     # HELLO | PLOT2 STARTS HERE
     #------------------------------------

  
  plot2 <- function(){
    
  
     
      # -----------------------------------------------
      # pfstat
      #------------------------------------
      
    pfstat <- function(){
      if (length(idb.dataset)==1){
      
        tkmessageBox(message="First simulate some data",icon="error",type="ok",parent=tt)
      
      }
      else{
        tttext <- tktoplevel(parent=.TkRoot)
        tkwm.title(tttext,"F Statistics")
        yscr <- tkscrollbar(tttext, repeatinterval=5,command=function(...)tkyview(txt,...))
        xscr <- tkscrollbar(tttext, repeatinterval=5,orient="horizontal",command=function(...)tkxview(txt,...))
        txt <- tktext(tttext,font="courier", wrap="none",yscrollcommand=function(...)tkset(yscr,...),xscrollcommand=function(...)tkset(xscr,...))
        tkinsert(txt,"end","Pairwise Fst\n\n\n")
        
        for (i in 1:nrow(idb.dataset$Fst)){
          for (j in 1:ncol(idb.dataset$Fst)){
            tkinsert(txt,"end",idb.dataset$Fst[i,j])
            tkinsert(txt,"end","\t")
          }
          tkinsert(txt,"end","\n")
        }
        tkinsert(txt,"end","\n\n\nFis\n\n\n")
        for (j in 1:length(idb.dataset$Fis)){
          tkinsert(txt,"end",idb.dataset$Fis[j])
          tkinsert(txt,"end","\t")
        }
        tkinsert(txt,"end","\n")
        
        tkgrid( txt,row=1,column=1)
        tkgrid( yscr,row=1,column=2, sticky="ns")
        tkgrid( xscr,row=2,column=1, sticky="we")
      }
    }
    
      # -----------------------------------------------
      # pdsigma
      #------------------------------------
  
    pdsigma <- function(){
      if (length(idb.dataset)==1){
      
        tkmessageBox(message="First simulate some data",icon="error",type="ok",parent=tt)
      
      }
      else{
        tttext <- tktoplevel(parent=.TkRoot)
        tkwm.title(tttext,"Dsigma2")
        yscr <- tkscrollbar(tttext, repeatinterval=5,command=function(...)tkyview(txt,...))
        xscr <- tkscrollbar(tttext, repeatinterval=5,orient="horizontal",command=function(...)tkxview(txt,...))
        txt <- tktext(tttext,font="courier", wrap="none",yscrollcommand=function(...)tkset(yscr,...),xscrollcommand=function(...)tkset(xscr,...))
        tkinsert(txt,"end","Dsigma2\n\n\n")
        
        for (j in 1:length(idb.dataset$Dsigma2)){
          tkinsert(txt,"end",idb.dataset$Dsigma2[j])
          tkinsert(txt,"end","\t")
        }
        tkinsert(txt,"end","\n")
        
        tkgrid( txt,row=1,column=1)
        tkgrid( yscr,row=1,column=2, sticky="ns")
        tkgrid( xscr,row=2,column=1, sticky="we")
      
   
      } 
    }

      # -----------------------------------------------
      # dsigma
      #------------------------------------
      
    pdiff <- function(){
      if (length(idb.dataset)==1){
      
        tkmessageBox(message="First simulate some data",icon="error",type="ok",parent=tt)
        
      }
      else{

        print(idb.dataset$diff.W)
        print(idb.dataset$diff.B)
        
        tttext <- tktoplevel(parent=.TkRoot)
        tkwm.title(tttext,"Variability of allele frequency")
        yscr <- tkscrollbar(tttext, repeatinterval=5,command=function(...)tkyview(txt,...))
        xscr <- tkscrollbar(tttext, repeatinterval=5,orient="horizontal",command=function(...)tkxview(txt,...))
        txt <- tktext(tttext,font="courier", wrap="none",yscrollcommand=function(...)tkset(yscr,...),xscrollcommand=function(...)tkset(xscr,...))
        tkinsert(txt,"end","Differentiation within populations\n\n\n")
        
        for (j in 1:length(idb.dataset$diff.W)){
          tkinsert(txt,"end",idb.dataset$diff.W[j])
          tkinsert(txt,"end","\t")
        }
        
        tkinsert(txt,"end","\n\n\nDifferentiation around barrier between populations\n\n\n")
        for (i in 1:nrow(idb.dataset$diff.B)){
          for (j in 1:ncol(idb.dataset$diff.B)){
            tkinsert(txt,"end",idb.dataset$diff.B[i,j])
            
            tkinsert(txt,"end","\t")
          }
          tkinsert(txt,"end","\n")
        }


        tkinsert(txt,"end","\n")


        
        tkgrid( txt,row=1,column=1)
        tkgrid( yscr,row=1,column=2, sticky="ns")
        tkgrid( xscr,row=2,column=1, sticky="we")
        
      }
    }
    
    buttonshowibd <- tkbutton(ttplot2,width=30,text="Show Simulated Data",command=GraficalIBD)
    
    labelfigures <-tklabel(ttplot2,text="-Graphics-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    labeltables <-tklabel(ttplot2,text="-Tables-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    labelspace1 <-tklabel(ttplot2,text=" ")

    labelfstat <-tklabel(ttplot2,text="-F statistic-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    #buttonfstat <- tkbutton(ttplot2,width=20,text="Fst and Fis",font="*-Helvetica-bold-i-*",command=pfstat)
    buttonfstat <- tkbutton(ttplot2,width=30,text="Fst and Fis",command=pfstat)

    labelspace2 <-tklabel(ttplot2,text=" ")

    labeldsigma <-tklabel(ttplot2,text="-Dsigma-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    #buttondsigma <- tkbutton(ttplot2,width=20,text="Dsigma2",font="*-Helvetica-bold-i-*",command=pdsigma)
    buttondsigma <- tkbutton(ttplot2,width=30,text="Dsigma2",command=pdsigma)
    labelspace3 <-tklabel(ttplot2,text=" ")

    labeldiff <-tklabel(ttplot2,text="-Differentiations-",font="*-Times-bold-i-normal--20-*",foreground="blue")
    #buttondiff <- tkbutton(ttplot2,width=20,text="Within and Between Individuals",font="*-Helvetica-bold-i-*",command=pdiff)
    buttondiff <- tkbutton(ttplot2,width=30,text="Within and Between Populations",command=pdiff)
    labelspace4 <-tklabel(ttplot2,text="                                      ")
    
    tkgrid(labelfigures,row=1,column=1,sticky="w")
    tkgrid(labelspace4,row=1,column=2,sticky="w")
    tkgrid(buttonshowibd,row=2,column=1,sticky="w")
    tkgrid(labelspace1,row=3,column=1,sticky="w")
    tkgrid(labelfstat,row=4,column=1,sticky="w")
    tkgrid(buttonfstat,row=5,column=1,sticky="w")
    tkgrid(labelspace2,row=6,column=1,sticky="w")
    tkgrid(labeldsigma,row=7,column=1,sticky="w")
    tkgrid(buttondsigma,row=8,column=1,sticky="w")
    tkgrid(labelspace3,row=9,column=1,sticky="w")
    tkgrid(labeldiff,row=10,column=1,sticky="w")
    tkgrid(buttondiff,row=11,column=1,sticky="w")
    
  }


     # -----------------------------------------------
     # HELLO | NFSTAT STARTS HERE
     #------------------------------------




    Nfstat <- function(){
      ttfstat <- tktoplevel()
      tkwm.title(ttnfstat,"F statistics")



      Runfstat <- function(){
        
        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        tkfocus(tttry)
        
        print("Starting...")
        Sys.sleep(0.1)
        
        
        err2 <- try(Fstat(genotypes=globalgenotypes,npop=as.integer(tclvalue(npop)),pop.mbrship=as.integer(tclvalue(pop.mbr))))
        tkdestroy(tttry)
        print("Done.")
        if (class(err2) == "try-error") 
          tkmessageBox(message=err2,icon="error",type="ok",parent=tt)              
        else{
          if(tclvalue(sep1)=="White space")
            tclvalue(sep1) <- " "
          if(tclvalue(sep2)=="White space")
            tclvalue(sep2) <- " " 
          tkmessageBox(message="Terminated with success",type="ok",parent=tt) 
          
          print(err2)
          
          if(tclvalue(sep1)==" ")
            tclvalue(sep1) <- "White space"
          if(tclvalue(sep2)==" ")
            tclvalue(sep2) <- "White space"
          
        }
      }

      #genotypes
      printbutton.widget <- tkbutton(ttnfstat,text="Save File",command=getprintfile,width=15)
      filelabel.widget <- tklabel(ttnfstat,textvariable=printfile,width=50)
      

      
      #npop
      npop<-tclVar(0)
      npoplabel.widget <- tklabel(ttnfstat,text="Number of populations",state="disabled")
      npop.widget <-tkentry(ttnfstat,width="20",textvariable=npop,state="disabled")
      
      tkgrid(npoplabel.widget,row=3,column=1,sticky="w")
      tkgrid(npop.widget,row=3,column=2,sticky="w")
      
            
      nextbutton <- tkbutton(ttnfstat,image=imagerun2,text="RUN >>",command=Runfstat)
      tkgrid(nextbutton,row=5,column=2,sticky="e")
    }
  
  
  
     # -----------------------------------------------
     # HELLO | FSTAT STARTS HERE
     #------------------------------------




    Gfstat <- function(){

#      ttfstat <- tktoplevel()

#      tkwm.title(ttfstat,"F statistics")



      Runfstat <- function(){


         #vect<- c()
         
         #vec<-unlist(strsplit(tclvalue(allelenumbers),","))
         #for (i in 1:length(vec)) vect[i]<-as.numeric(vec[i])

        # allele<- tclVar("")

        # tcl("append",allele,"c","(",tclvalue(allelenumbers),")")
        #if(tclvalue(cb)==1){
          tttry <- tktoplevel(parent=.TkRoot)
          tkgrab(tttry)
          tkwm.geometry(tttry, "+200+200")
          tkwm.title(tttry,"wait")
          warn<-tklabel(tttry,image=imagepleasewait)
          tkpack(warn)
          tkfocus(tttry)
          
          print("Starting...")
          Sys.sleep(0.1)
          
          
          err <- try(Fstat.output(genotypes=globalgenotypes,path.mcmc=tclvalue(outputdir)),silent=TRUE)
          tkdestroy(tttry)
          print("Done.")
#          print (err)
          if (class(err) == "try-error") 
            tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
          else{
            if(tclvalue(sep1)=="White space")
              tclvalue(sep1) <- " "
            if(tclvalue(sep2)=="White space")
              tclvalue(sep2) <- " " 
            tkmessageBox(message="Terminated with success",type="ok",parent=tt) 
            
            
            write.table(err$Fis,file=paste(tclvalue(outputdir),"Fis.txt",sep=""),row.names=FALSE,col.names=FALSE)
            write.table(err$Fst,file=paste(tclvalue(outputdir),"Fst.txt",sep=""),row.names=FALSE,col.names=FALSE)
            
            Showtext("Fis.txt")
            Showtext("Fst.txt")
            
            
            
            if(tclvalue(sep1)==" ")
              tclvalue(sep1) <- "White space"
            if(tclvalue(sep2)==" ")
              tclvalue(sep2) <- "White space"
            
          }
        }
        #}
        
        #if(tclvalue(cb2)==1){
        #  tttry <- tktoplevel(parent=.TkRoot)
        #  tkgrab(tttry)
        #  tkwm.geometry(tttry, "+200+200")
        #  tkwm.title(tttry,"wait")
        #  warn<-tklabel(tttry,image=imagepleasewait)
        #  tkpack(warn)
        #  tkfocus(tttry)
        #  
        #  print("Starting...")
        #  Sys.sleep(0.1)
        #  
        #  
        #  err2 <- try(Fstat(genotypes=globalgenotypes,npop=as.integer(tclvalue(npop)),pop.mbrship=as.integer(tclvalue(pop.mbr))))
        #  tkdestroy(tttry)
        #  print("Done.")
        #  if (class(err2) == "try-error") 
        #    tkmessageBox(message=err2,icon="error",type="ok",parent=tt)              
        #  else{
        #    if(tclvalue(sep1)=="White space")
        #      tclvalue(sep1) <- " "
        #    if(tclvalue(sep2)=="White space")
        #      tclvalue(sep2) <- " " 
        #    tkmessageBox(message="Terminated with success",type="ok",parent=tt) 
        #    
        #    print(err2)
        #    
        #    if(tclvalue(sep1)==" ")
        #      tclvalue(sep1) <- "White space"
        #    if(tclvalue(sep2)==" ")
        #      tclvalue(sep2) <- "White space"
        #    
        #  }
        #}
      #}
      
      #viewothers <- function(){
      #  
      #  if(tclvalue(cb2)==0){
      #    tkconfigure(npoplabel.widget,state="disabled")
      #    tkconfigure(npop.widget,state="disabled")
      #    tkconfigure(pop.mbrlabel.widget,state="disabled")
      #    tkconfigure(pop.mbr.widget,state="disabled")
      #  }
      #  else{
      #    tkconfigure(npoplabel.widget,state="normal")
      #    tkconfigure(npop.widget,state="normal")
      #    tkconfigure(pop.mbrlabel.widget,state="normal")
      #    tkconfigure(pop.mbr.widget,state="normal")
      #  }
      #}
      
      
      #cb<-tclVar(0)
      #cb2<-tclVar(0)
      
      #cblabel.widget <- tklabel(ttfstat,text="old")
      #cb.widget <- tkcheckbutton(ttfstat,variable=cb,onvalue=1,offvalue=0)

      #cb2label.widget <- tklabel(ttfstat,text="new")
      #cb2.widget <- tkcheckbutton(ttfstat,variable=cb2,onvalue=1,offvalue=0,command=viewothers)

      
      #npop
      #  npop<-tclVar(0)
      #  npoplabel.widget <- tklabel(ttfstat,text="Number of populations",state="disabled")
      #  npop.widget <-tkentry(ttfstat,width="20",textvariable=npop,state="disabled")

      #pop.mbrship
      #  pop.mbr <-tclVar(0)
      #pop.mbrlabel.widget <- tklabel(ttfstat,text="Population membership",state="disabled")
      #pop.mbr.widget <-tkentry(ttfstat,width="20",textvariable=pop.mbr,state="disabled")

      #tkgrid(cblabel.widget,row=1,column=1,sticky="w")
      #tkgrid(cb.widget,row=1,column=2,sticky="w")
      #tkgrid(cb2label.widget,row=2,column=1,sticky="w")
      #tkgrid(cb2.widget,row=2,column=2,sticky="w")
      #tkgrid(npoplabel.widget,row=3,column=1,sticky="w")
      #tkgrid(npop.widget,row=3,column=2,sticky="w")
      #tkgrid(pop.mbrlabel.widget,row=4,column=1,sticky="w")
      #tkgrid(pop.mbr.widget,row=4,column=2,sticky="w")
       #strallele <- paste(globalallele.numbers, sep = ",", collapse =",")#globalallele.numbers[0]
       #allelenumbers
       #allelenumbers <- tclVar(strallele)
       #allelenumbers.widget <-tkentry(ttfstat,width="20",textvariable=allelenumbers)
       #allelenumberslabel.widget <- tklabel(ttfstat,text="Number of alleles per locus (E.g: n,m,..):")

       #tkgrid(allelenumberslabel.widget,row=1,column=1,sticky="w")
       #tkgrid(allelenumbers.widget,row=1,column=2,sticky="w")
      
      #nextbutton <- tkbutton(ttfstat,image=imagerun2,text="RUN >>",command=Runfstat)
      #tkgrid(nextbutton,row=5,column=2,sticky="e")
      
       
       
       
       Runfstat()
       

     }

     # -----------------------------------------------
     # HELLO | NON-IBD STARTS HERE
     #------------------------------------

     
   SimnonIBD <- function(){

     #this is going to be complicated



     runnonibd <- function(){

       if(tclvalue(sep1)=="White space")
         tclvalue(sep1) <- " "
       if(tclvalue(sep2)=="White space")
        tclvalue(sep2) <- " " 
       

       vect1<- c()
       #vect2 <- c()
       
       vec1<-unlist(strsplit(tclvalue(nall),","))
       for (i in 1:length(vec1)) vect1[i]<-as.numeric(vec1[i])
       
       #vec2<-unlist(strsplit(tclvalue(param),","))
       #for (i in 1:length(vec2)) vect2[[i]]<-as.numeric(vec2[i])
       
       if(tclvalue(nloc)!=length(vect1)){
         tkmessageBox(message="Number of locus must be equal to the length of the number of alleles per locus",icon="error",type="ok",parent=tt)
       }
       else{
         #beta <- as.numeric(tclvalue(absmax))
         #if(beta < as.numeric(tclvalue(ordmax)))
         #  beta <- as.numeric(tclvalue(ordmax))
         #beta <- beta*100
         tttry <- tktoplevel(parent=.TkRoot)
         tkgrab(tttry)
         tkwm.geometry(tttry, "+200+200")
         tkwm.title(tttry,"wait")
         warn<-tklabel(tttry,image=imagepleasewait)
         tkpack(warn)
         tkfocus(tttry)
         
         print("Starting...")
         Sys.sleep(0.1)
         if(length(globalcoordinates)==1 || tclvalue(prevcoord)==0)
           idb.dataset <<- try(simdata( nindiv=as.numeric(tclvalue(nindiv)),coord.lim=c(as.numeric(tclvalue(absmin)),as.numeric(tclvalue(absmax)),as.numeric(tclvalue(ordmin)),as.numeric(tclvalue(ordmax))),number.nuclei=as.numeric(tclvalue(nuclei)),allele.numbers=vect1,IBD=FALSE,npop=as.numeric(tclvalue(npop)),give.freq.grid=as.logical(tclvalue(freq.grid)),give.tess.grid=as.logical(tclvalue(tess.grid)),npix=c(as.numeric(tclvalue(npixh)),as.numeric(tclvalue(npixv))),comp.Fst=as.logical(tclvalue(comp.Fst))),silent=TRUE)
         else
           idb.dataset <<- try(simdata( nindiv=as.numeric(tclvalue(nindiv)),coord.indiv=globalcoordinates,coord.lim=c(as.numeric(tclvalue(absmin)),as.numeric(tclvalue(absmax)),as.numeric(tclvalue(ordmin)),as.numeric(tclvalue(ordmax))),number.nuclei=as.numeric(tclvalue(nuclei)),allele.numbers=vect1,IBD=FALSE,npop=as.numeric(tclvalue(npop)),give.freq.grid=as.logical(tclvalue(freq.grid)),give.tess.grid=as.logical(tclvalue(tess.grid)),npix=c(as.numeric(tclvalue(npixh)),as.numeric(tclvalue(npixv))),comp.Fst=as.logical(tclvalue(comp.Fst))),silent=TRUE)
         tkdestroy(tttry)
         print("Done.")
         if (class(idb.dataset) == "try-error"){
           tkmessageBox(message=idb.dataset,icon="error",type="ok",parent=tt)              
           idb.dataset <<- 0
         }
         else{
           tkmessageBox(message="Terminated with success",type="ok",parent=tt)
           #print(idb.dataset)
           #print(idb.dataset$Fst)
           #print(idb.dataset$Fis)
           #print(idb.dataset$Fit)
           globalcoordinates <<- t(idb.dataset$coord.indiv)
           #tkconfigure(coordownlabel.widget,text="Coordinates:    Simulated panmitic data loaded")
           tclvalue(labelcoordtext) <- "Coordinates:    Simulated panmitic data loaded"
           globalgenotypes <<- idb.dataset$genotypes
           #tkconfigure(genodownlabel.widget,text="Genotype:       Simulated panmitic data loaded")
           tclvalue(labelgenotext) <- "Genotype:      Simulated panmitic data loaded"
           globalallele.numbers <<- idb.dataset$allele.numbers

           #print(idb.dataset$coord.indiv)
           
           if(tclvalue(save)==1){
             
             auxcoord <-tclVar()
             tclvalue(auxcoord) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save coordinates file to:"))
             auxgen <-tclVar()
             tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save genotypes file to:"))
             
             write.table(t(idb.dataset$coord.indiv),file=tclvalue(auxcoord),sep=tclvalue(sep1),row.names=FALSE,col.names=FALSE)
             write.table(idb.dataset$genotypes,file = tclvalue(auxgen),sep=tclvalue(sep2),row.names=FALSE,col.names=FALSE)
           }
         }
       }
       if(tclvalue(sep1)==" ")
         tclvalue(sep1) <- "White space"
       if(tclvalue(sep2)==" ")
         tclvalue(sep2) <- "White space"
       
     }
     
     #nindiv
     nindiv=tclVar("0")
     nindiv.widget <-tkentry(ttsimf,width="20",textvariable=nindiv)
     nindivlabel.widget <- tklabel(ttsimf,text="Number of individuals:")
     
     #coord.lim

     coordxlabel.widget <- tklabel(ttsimf,text="Limits of geographical domain:")
     
     absmin <- tclVar(0)
     absmax <- tclVar(1)
     absmin.widget <-tkentry(ttsimf,width="8",textvariable=absmin)
     absmax.widget <-tkentry(ttsimf,width="8",textvariable=absmax)
     abslabel.widget <- tklabel(ttsimf,text="   abs (min|max) :")
     
     ordmin <- tclVar(0)
     ordmax <- tclVar(1)
     ordmin.widget <-tkentry(ttsimf,width="8",textvariable=ordmin)
     ordmax.widget <-tkentry(ttsimf,width="8",textvariable=ordmax)
     ordlabel.widget <- tklabel(ttsimf,text="   ord (min|max) :")

     #rate
     #rate=tclVar("0")
     #rate.widget <-tkentry(ttsimf,width="20",textvariable=rate)
     #ratelabel.widget <- tklabel(ttsimf,text="Rate of the Poisson process:")

     # number.nuclei 

     nuclei=tclVar("0")
     nuclei.widget <-tkentry(ttsimf,width="20",textvariable=nuclei)
     nucleilabel.widget <- tklabel(ttsimf,text="Number of nuclei in tessellation:")
     
     #coord.nuclei

     #cuclei <- tclVar(0)
     #cuclei.widget <-tkentry(ttsimf,width="20",textvariable=cuclei)
     #cucleilabel.widget <- tklabel(ttsimf,text="Coordinates of nuclei of Voronoi tessellation:")

     #color.nuclei

     #coclei <- tclVar(0)
     #coclei.widget <-tkentry(ttsimf,width="20",textvariable=coclei)
     #cocleilabel.widget <- tklabel(ttsimf,text="Population labels of the nuclei (E.g: n,m,..):")

     #allele.numbers

     nloc <- tclVar(0)
     nloc.widget <-tkentry(ttsimf,width="20",textvariable=nloc)
     nloclabel.widget <- tklabel(ttsimf,text="Number of loci:")
     
     
     #allele.numbers

     nall <- tclVar()
     nall.widget <-tkentry(ttsimf,width="20",textvariable=nall)
     nalllabel.widget <- tklabel(ttsimf,text="Number of alleles per locus (E.g: 10,3,8,..):")
     
     #model

#     modellabel.widget <- tklabel(ttsimf,text="Spatial covariance model:")
#     model <- tclVar("besset")
#     wmodel <- .Tk.subwin(ttsimf)
#     modeloptionmenu.widget <-  tcl("tk_optionMenu",wmodel,model,"besset","cauchy","cauchytbm","circular","constant","cone","cubic","cutoff","dampedcosine","exponential","fractalB","FD","fractgauss","gauss","gencauchy","gengneiting","gneiting","hyperbolic","iacocesare","lgd1","mastein","nsst","nsst2","nugget","penta","power","qexponential","spherical","stable","Stein","steinst1","wave","whittlematern")

     #param

#     param <- tclVar("")
#     param.widget <-tkentry(ttsimf,width="20",textvariable=param)
#     paramlabel.widget <- tklabel(ttsimf,text="Parameters (E.g:mean, variance, nugget, scale, ...):")

     #Beta
#     beta <- tclVar("")
#     beta.widget <-tkentry(ttsimf,width="20",textvariable=beta)
#     betalabel.widget <- tklabel(ttsimf,text="Spacial correlation parameter:")
     
     #npop
     npop <- tclVar("")
     npop.widget <-tkentry(ttsimf,width="20",textvariable=npop)
     npoplabel.widget <- tklabel(ttsimf,text="Number of populations:")

     #freq.grid
     freq.gridlabel.widget <- tklabel(ttsimf,text="Return frequencies on grid:")
     freq.grid <- tclVar("FALSE")
     wfreq.grid <- .Tk.subwin(ttsimf)
     freq.gridoptionmenu.widget <-  tcl("tk_optionMenu",wfreq.grid,freq.grid,"FALSE","TRUE")

     #tess.grid
     tess.gridlabel.widget <- tklabel(ttsimf,text="Return population membership on grid:")
     tess.grid <- tclVar("FALSE")
     wtess.grid <- .Tk.subwin(ttsimf)
     tess.gridoptionmenu.widget <-  tcl("tk_optionMenu", wtess.grid,tess.grid,"FALSE","TRUE")

     #npix

     npixh <- tclVar(50)
     npixv <- tclVar(50)
     npixh.widget <-tkentry(ttsimf,width="8",textvariable=npixh)
     npixv.widget <-tkentry(ttsimf,width="8",textvariable=npixv)
     npixlabel.widget <- tklabel(ttsimf,text="Number of pixels for representation (hor|ver):")
     
     #comp.Fst
     comp.Fstlabel.widget <- tklabel(ttsimf,text="Compute F statistics:")
     comp.Fst <- tclVar("FALSE")
     wcomp.Fst <- .Tk.subwin(ttsimf)
     comp.Fstoptionmenu.widget <-  tcl("tk_optionMenu", wcomp.Fst,comp.Fst,"FALSE","TRUE")
     
     #comp.Dsigma2
 #    comp.Dsigma2label.widget <- tklabel(ttsimf,text="Compute SIMF index Sigma2:")
 #    comp.Dsigma2 <- tclVar("FALSE")
 #    wcomp.Dsigma2 <- .Tk.subwin(ttsimf)
 #    comp.Dsigma2optionmenu.widget <-  tcl("tk_optionMenu",wcomp.Dsigma2,comp.Dsigma2,"FALSE","TRUE")

     #comp.height
#     comp.heightlabel.widget <- tklabel(ttsimf,text="Compute the \"height\" of the barriers:")
#     comp.height <- tclVar("FALSE")
#     wcomp.height <- .Tk.subwin(ttsimf)
#     comp.heightoptionmenu.widget <-  tcl("tk_optionMenu", wcomp.height,comp.height,"FALSE","TRUE")

     
     #width

#     hwidth <- tclVar(1)
#     hwidth.widget <-tkentry(ttsimf,width="20",textvariable=hwidth)
#     hwidthlabel.widget <- tklabel(ttsimf,text="Width around the barriers:")

     #plot.pairs.borders
#     plot.pairs.borderslabel.widget <- tklabel(ttsimf,text="Plot pairs of individuals:")
#     plot.pairs.borders <- tclVar("FALSE")
#     wplot.pairs.borders <- .Tk.subwin(ttsimf)
#     plot.pairs.bordersoptionmenu.widget <-  tcl("tk_optionMenu", wplot.pairs.borders,plot.pairs.borders,"FALSE","TRUE")

     prevcoord <- tclVar(0)
     
     prevcoordlabel.widget <- tklabel(ttsimf,text="Use loaded coordinates file:")
     prevcoord.widget <- tkcheckbutton(ttsimf,variable=save,onvalue=1, offvalue=0)


     
     #savefiles
     save <- tclVar(0)
     savelabel.widget <- tklabel(ttsimf,text="Save coordinates and genotypes files:")
     save.widget <- tkcheckbutton(ttsimf,variable=save,onvalue=1, offvalue=0)
     
     tkgrid(nindivlabel.widget,row=1,column=1,sticky="w")
     tkgrid(nindiv.widget,row=1,column=2,columnspan=2,sticky="w")
     tkgrid(coordxlabel.widget,row=2,column=1,sticky="w")
     tkgrid(abslabel.widget,row=3,column=1,sticky="w")
     tkgrid(absmin.widget,row=3,column=2,sticky="w")
     tkgrid(absmax.widget,row=3,column=3,sticky="w")
     tkgrid(ordlabel.widget,row=4,column=1,sticky="w")
     tkgrid(ordmin.widget,row=4,column=2,sticky="w")
     tkgrid(ordmax.widget,row=4,column=3,sticky="w")
     #tkgrid(ratelabel.widget,row=5,column=1,sticky="w")
     #tkgrid(rate.widget,row=5,column=2,columnspan=2,sticky="w")
     tkgrid(nucleilabel.widget,row=5,column=1,sticky="w")
     tkgrid(nuclei.widget,row=5,column=2,columnspan=2,sticky="w")
     tkgrid(nloclabel.widget,row=6,column=1,sticky="w")
     tkgrid(nloc.widget,row=6,column=2,columnspan=2,sticky="w")
     tkgrid(nalllabel.widget,row=7,column=1,sticky="w")
     tkgrid(nall.widget,row=7,column=2,columnspan=2,sticky="w")
##     tkgrid(modellabel.widget,row=8,column=1,sticky="w")
##     tkgrid(wmodel,row=8,column=2,columnspan=2,sticky="w")
##     tkgrid(paramlabel.widget,row=9,column=1,sticky="w")
##     tkgrid(param.widget,row=9,column=2,columnspan=2,sticky="w")
#     tkgrid(betalabel.widget,row=9,column=1,columnspan=2,sticky="w")
#     tkgrid(beta.widget,row=9,column=2,columnspan=2,sticky="w")
     tkgrid(npoplabel.widget,row=10,column=1,sticky="w")
     tkgrid(npop.widget,row=10,column=2,columnspan=2,sticky="w")
     tkgrid(freq.gridlabel.widget,row=11,column=1,sticky="w")
     tkgrid(wfreq.grid,row=11,column=2,columnspan=2,sticky="w")
     tkgrid(tess.gridlabel.widget,row=12,column=1,sticky="w")
     tkgrid(wtess.grid,row=12,column=2,columnspan=2,sticky="w")
     tkgrid(npixlabel.widget,row=13,column=1,sticky="w")
     tkgrid(npixh.widget,row=13,column=2,sticky="w")
     tkgrid(npixv.widget,row=13,column=3,sticky="w")
     tkgrid(comp.Fstlabel.widget,row=14,column=1,sticky="w")
     tkgrid(wcomp.Fst,row=14,column=2,columnspan=2,sticky="w")
#     tkgrid(comp.Dsigma2label.widget,row=15,column=1,sticky="w")
#     tkgrid(wcomp.Dsigma2,row=15,column=2,columnspan=2,sticky="w")
#     tkgrid(comp.heightlabel.widget,row=16,column=1,sticky="w")
#     tkgrid(wcomp.height,row=16,column=2,columnspan=2,sticky="w")
#     tkgrid(hwidthlabel.widget,row=17,column=1,sticky="w")
#     tkgrid(hwidth.widget,row=17,column=2,columnspan=2,sticky="w")
#     tkgrid(plot.pairs.borderslabel.widget,row=18,column=1,sticky="w")
#     tkgrid(wplot.pairs.borders,row=18,column=2,columnspan=2,sticky="w")
     tkgrid(prevcoordlabel.widget,row=15,column=1,columnspan=2,sticky="w")
     tkgrid(prevcoord.widget,row=15,column=2,columnspan=2,sticky="w")

     tkgrid(savelabel.widget,row=19,column=1,sticky="w")
     tkgrid(save.widget,row=19,column=2,columnspan=2,sticky="w")
     
     labelspace <-tklabel(ttsimf,text=" ")
     tkgrid(labelspace,row=20,column=1)
     
     nextbutton <- tkbutton(ttsimf,image=imagerun2,text="RUN >>",command=runnonibd)
     tkgrid(nextbutton,row=21,column=2,columnspan=2,sticky="e")
   }


  simFModel <- function(){

    
    Runsimf <- function(){
      
      if(tclvalue(sep1)=="White space")
        tclvalue(sep1) <- ""
      if(tclvalue(sep2)=="White space")
        tclvalue(sep2) <- ""
      
       # vec1<-as.vector(strsplit(tclvalue(coclei),","))
       #  for (i in 1:length(vec1)) vec1[[i]]<-as.numeric(vec1[[i]])

      vect2<- c()
      vect3<- c()
      
      vec2<-unlist(strsplit(tclvalue(nall),","))
      for (i in 1:length(vec2)) vect2[i]<-as.numeric(vec2[i])
      
      vec3<-unlist(strsplit(tclvalue(drift),","))
      for (i in 1:length(vec3)) vect3[[i]]<-as.numeric(vec3[i])

      tttry <- tktoplevel(parent=.TkRoot)
      tkgrab(tttry)
      tkwm.geometry(tttry, "+200+200")
      tkwm.title(tttry,"wait")
      warn<-tklabel(tttry,image=imagepleasewait)
      tkpack(warn)
      tkfocus(tttry)
      
      print("Starting...")
      Sys.sleep(0.1)

      if(length(globalcoordinates)==1 || tclvalue(prevcoord)==0){
        #try(sim <- simFmodel(nindiv=as.numeric(tclvalue(nindiv)),coordinates=NULL,coord.lim=c(as.numeric(tclvalue(coordminx)),as.numeric(tclvalue(coordmaxx)),as.numeric(tclvalue(coordminy)),as.numeric(tclvalue(coordmaxy))),number.nuclei=as.numeric(tclvalue(nuclei)),coord.nuclei=NULL,color.nuclei=NULL,nall=vect2,npop=as.numeric(tclvalue(npop)),freq.model=tclvalue(freq),drift=vect3,seed=as.numeric(tclvalue(seed)),plots =as.logical(tclvalue(plots)),ploth =as.logical(tclvalue(ploth))),silent = TRUE)
        # And for history....
        # sim <- simFmodel(nindiv=as.numeric(tclvalue(nindiv)),coordinates=NULL,coord.lim=c(as.numeric(tclvalue(coordminx)),as.numeric(tclvalue(coordmaxx)),as.numeric(tclvalue(coordminy)),as.numeric(tclvalue(coordmaxy))),number.nuclei=as.numeric(tclvalue(nuclei)),coord.nuclei=NULL,color.nuclei=NULL,nall=vect2,npop=as.numeric(tclvalue(npop)),freq.model=tclvalue(freq),drift=vect3,seed=as.numeric(tclvalue(seed)),plots =as.logical(tclvalue(plots)),ploth =as.logical(tclvalue(ploth)))
        sim <- try(simFmodel(nindiv=as.numeric(tclvalue(nindiv)),coord.lim=c(as.numeric(tclvalue(coordminx)),as.numeric(tclvalue(coordmaxx)),as.numeric(tclvalue(coordminy)),as.numeric(tclvalue(coordmaxy))),number.nuclei=as.numeric(tclvalue(nuclei)),nall=vect2,npop=as.numeric(tclvalue(npop)),freq.model=tclvalue(freq),drift=vect3,seed=as.numeric(tclvalue(seed)),plots =as.logical(tclvalue(plots)),ploth =as.logical(tclvalue(ploth))),silent=TRUE)
      }
      else  {
        #try(sim <- simFmodel(nindiv=as.numeric(tclvalue(nindiv)),coordinates=globalcoordinates,coord.lim=c(as.numeric(tclvalue(coordminx)),as.numeric(tclvalue(coordmaxx)),as.numeric(tclvalue(coordminy)),as.numeric(tclvalue(coordmaxy))),number.nuclei=as.numeric(tclvalue(nuclei)),coord.nuclei=NULL,color.nuclei=NULL,nall=vect2,npop=as.numeric(tclvalue(npop)),freq.model=tclvalue(freq),drift=vect3,seed=as.numeric(tclvalue(seed)),plots =as.logical(tclvalue(plots)),ploth =as.logical(tclvalue(ploth))),silent = TRUE)
        sim <- try(simFmodel(nindiv=as.numeric(tclvalue(nindiv)),coordinates=globalcoordinates,coord.lim=c(as.numeric(tclvalue(coordminx)),as.numeric(tclvalue(coordmaxx)),as.numeric(tclvalue(coordminy)),as.numeric(tclvalue(coordmaxy))),number.nuclei=as.numeric(tclvalue(nuclei)),nall=vect2,npop=as.numeric(tclvalue(npop)),freq.model=tclvalue(freq),drift=vect3,seed=as.numeric(tclvalue(seed)),plots =as.logical(tclvalue(plots)),ploth =as.logical(tclvalue(ploth))),silent=TRUE)
      }
        
      tkdestroy(tttry)
      print("Done.")
      if (class(sim) == "try-error") 
        tkmessageBox(message=sim,icon="error",type="ok",parent=tt)              
      else{
        
        tkmessageBox(message="Terminated with success",type="ok",parent=tt)
        globalcoordinates <<- sim$coordinates
        tclvalue(labelcoordtext) <- "Coordinates:    Simulated non-IBD data loaded"
        #tkconfigure(coordownlabel.widget,text="Coordinates:    Simulated non-IBD data loaded")
        globalgenotypes <<- sim$genotypes
        #tkconfigure(genodownlabel.widget,text="Genotype:       Simulated non-IBD data loaded")
        tclvalue(labelgenotext) <- "Genotype:    Simulated non-IBD data loaded"

        #allele.numbers <<- sim$allele
        globalallele.numbers <<- sim$allele
      
        if(tclvalue(save)==1){
          
          auxcoord <-tclVar()
          tclvalue(auxcoord) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save coordinates file to:"))
          auxgen <-tclVar()
          tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save genotypes file to:"))
          
          
          write.table(sim$coordinates,file=tclvalue(auxcoord),sep=tclvalue(sep1),row.names=FALSE,col.names=FALSE)
          write.table(sim$genotypes,file = tclvalue(auxgen),sep=tclvalue(sep2),row.names=FALSE,col.names=FALSE)
        }
      }
      if(tclvalue(sep1)=="")
        tclvalue(sep1) <- "White space"
      if(tclvalue(sep2)=="")
        tclvalue(sep2) <- "White space"
      
    }
    
    #previou loaded coordinates
    prevcoord <- tclVar(0)
    
    prevcoordlabel.widget <- tklabel(ttsimf,text="Use loaded coordinates file:")
    prevcoord.widget <- tkcheckbutton(ttsimf,variable=save,onvalue=1, offvalue=0)

    
    
      #nindiv

    nindiv <- tclVar(0)
    nindiv.widget <-tkentry(ttsimf,width="20",textvariable=nindiv)
    nindivlabel.widget <- tklabel(ttsimf,text="Number of individuals:")
    
      #coord.lim

    coordxlabel.widget <- tklabel(ttsimf,text="Limits of spatial domain:")
    
    coordminx <- tclVar(0)
    coordmaxx <- tclVar(1)
    coordminx.widget <-tkentry(ttsimf,width="8",textvariable=coordminx)
    coordmaxx.widget <-tkentry(ttsimf,width="8",textvariable=coordmaxx)
    coordxxlabel.widget <- tklabel(ttsimf,text="X (min|max) :")
    
    coordminy <- tclVar(0)
    coordmaxy <- tclVar(1)
    coordminy.widget <-tkentry(ttsimf,width="8",textvariable=coordminy)
    coordmaxy.widget <-tkentry(ttsimf,width="8",textvariable=coordmaxy)
    coordxylabel.widget <- tklabel(ttsimf,text="Y (min|max) :")

      #number.nuclei

    nuclei <- tclVar(0)
    nuclei.widget <-tkentry(ttsimf,width="20",textvariable=nuclei)
    nucleilabel.widget <- tklabel(ttsimf,text="Number of nuclei in Voronoi tessellation:")

      #coord.nuclei

      #cuclei <- tclVar(0)
      #cuclei.widget <-tkentry(ttsimf,width="20",textvariable=cuclei)
      #cucleilabel.widget <- tklabel(ttsimf,text="Coordinates of nuclei of Voronoi tessellation:")

      #color.nuclei

      #coclei <- tclVar(0)
      #coclei.widget <-tkentry(ttsimf,width="20",textvariable=coclei)
      #cocleilabel.widget <- tklabel(ttsimf,text="Population labels of the nuclei (E.g: n,m,..):")

      #nall

    nall <- tclVar(0)
    nall.widget <-tkentry(ttsimf,width="20",textvariable=nall)
    nalllabel.widget <- tklabel(ttsimf,text="Number of alleles at each locus (Eg: A1,A2,..):")
    

      #populations
    npop <- tclVar("1")
    wpop  <- .Tk.subwin(ttsimf)
    pop<-tcl("spinbox", wpop,"-textvariable",npop,"-width",5,"-increment",1,"-from",1,"-to",100000)
    npoplabel.widget <- tklabel(ttsimf,text="Number of populations:")
    
      #freq.model

    freqlabel.widget <- tklabel(ttsimf,text="Frequency Model:")
    freq <- tclVar("Dirichlet")
    wfreq <- .Tk.subwin(ttsimf)
    freqoptionmenu.widget <-  tcl("tk_optionMenu",wfreq,freq,"Falush","Dirichlet")

      #drift

    drift <- tclVar(0)
    drift.widget <-tkentry(ttsimf,width="20",textvariable=drift)
    driftlabel.widget <- tklabel(ttsimf,text="Drift factors for Falush Model (E.g: n,m,..):")
    
      #seed

    seed <- tclVar(0)
    seed.widget <-tkentry(ttsimf,width="20",textvariable=seed)
    seedlabel.widget <- tklabel(ttsimf,text="Initial seed:")
    
      #plots

    plotslabel.widget <- tklabel(ttsimf,text="Plot spatial coord:")
    plots <- tclVar("TRUE")
    wplots <- .Tk.subwin(ttsimf)
    plotsoptionmenu.widget <-  tcl("tk_optionMenu",wplots,plots,"FALSE","TRUE")
    
      #ploth

    plothlabel.widget <- tklabel(ttsimf,text="Barplots for allele frequencies:")
    ploth <- tclVar("TRUE")
    wploth <- .Tk.subwin(ttsimf)
    plothoptionmenu.widget <-  tcl("tk_optionMenu",wploth,ploth,"FALSE","TRUE")

    #save files to disk
    save <- tclVar(0)
    savelabel.widget <- tklabel(ttsimf,text="Save coordinates and genotypes files?")
    save.widget <- tkcheckbutton(ttsimf,variable=save,onvalue=1, offvalue=0)
    #print(tclvalue(save)
  
    tkgrid(nindivlabel.widget,row=2,column=1,sticky="w")
    tkgrid(nindiv.widget,row=2,column=2,columnspan=2,sticky="w")
    tkgrid(coordxlabel.widget,row=3,column=1,sticky="w")
    tkgrid(coordxxlabel.widget,row=4,column=1,sticky="w")
    tkgrid(coordminx.widget,row=4,column=2,sticky="w")
    tkgrid(coordmaxx.widget,row=4,column=3,sticky="w")
    tkgrid(coordxylabel.widget,row=5,column=1,sticky="w")
    tkgrid(coordminy.widget,row=5,column=2,sticky="w")
    tkgrid(coordmaxy.widget,row=5,column=3,sticky="w")
    tkgrid(nucleilabel.widget,row=6,column=1,sticky="w")
    tkgrid(nuclei.widget,row=6,column=2,columnspan=2,sticky="w")
    #tkgrid(cucleilabel.widget,row=6,column=1,sticky="w")
    #tkgrid(cuclei.widget,row=6,column=2,columnspan=2,sticky="w")
    #tkgrid(cocleilabel.widget,row=7,column=1,sticky="w")
    #tkgrid(coclei.widget,row=7,column=2,columnspan=2,sticky="w")
    tkgrid(nalllabel.widget,row=8,column=1,sticky="w")
    tkgrid(nall.widget,row=8,column=2,columnspan=2,sticky="w")
    tkgrid(npoplabel.widget,row=9,column=1,sticky="w")
    tkgrid(wpop,row=9,column=2,columnspan=2,sticky="w")
    tkgrid(npoplabel.widget,row=10,column=1,sticky="w")
    tkgrid(wpop,row=10,column=2,columnspan=2,sticky="w")
    tkgrid(freqlabel.widget,row=11,column=1,sticky="w")
    tkgrid(wfreq,row=11,column=2,columnspan=2,sticky="w")
    tkgrid(driftlabel.widget,row=12,column=1,sticky="w")
    tkgrid(drift.widget,row=12,column=2,columnspan=2,sticky="w")
    tkgrid(seedlabel.widget,row=13,column=1,sticky="w")
    tkgrid(seed.widget,row=13,column=2,columnspan=2,sticky="w")
    tkgrid(plotslabel.widget,row=14,column=1,sticky="w")
    tkgrid(wplots,row=14,column=2,columnspan=2,sticky="w")
    tkgrid(plothlabel.widget,row=15,column=1,sticky="w")
    tkgrid(wploth,row=15,column=2,columnspan=2,sticky="w")
    tkgrid(prevcoordlabel.widget,row=16,column=1,sticky="w")
    tkgrid(prevcoord.widget,row=16,column=2,columnspan=2,sticky="w")
    tkgrid(savelabel.widget,row=17,column=1,sticky="w")
    tkgrid(save.widget,row=17,column=2,columnspan=2,sticky="w")

    nextbutton <- tkbutton(ttsimf,image=imagerun2,text="RUN >>",command=Runsimf)
    tkgrid(nextbutton,row=18,column=2,columnspan=2,sticky="e")
    
  }

     # -----------------------------------------------
     # HELLO | GL2GP STARTS HERE
     #------------------------------------

  Convert <- function(){
    
    ttcon<-tktoplevel()
    filename<- tclVar("")

    tkwm.title(ttcon,"Convert loaded data into genepop format")
    
    gltgp <- function(){
      
      if(tclvalue(filename)=="" | length(globalcoordinates)==1 | length(globalgenotypes)==1){
        tkmessageBox(message="You must define filename, coordinates file and genotypes file",icon="error",type="ok")
        
      }
      else{

        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        tkfocus(tttry)
        
        print("Starting...")
        Sys.sleep(0.1)
        
        
        err <- try(gl2gp(coordinates=t(globalcoordinates),genotypes=globalgenotypes,file=tclvalue(filename)),silent=TRUE)
        tkdestroy(tttry)
        print("Done.")
        if (class(err) == "try-error") 
          tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
        else
          tkmessageBox(message="Terminated with success",type="ok",parent=tt)
      }
      
    }
    
    setgenepopfile <- function(){
      tclvalue(filename) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",title="Choose a File"))   
    }
    

    genepopbutton.widget <- tkbutton(ttcon,text="Choose file name",command=setgenepopfile,width=15)
    filelabel.widget <- tklabel(ttcon,textvariable=filename,width=50)
    
    
    tkgrid(genepopbutton.widget,row=1,column=1,sticky="w")
    tkgrid(filelabel.widget,row=1,column=2,sticky="w")

    labelspace <-tklabel(ttcon,text=" ")
    tkgrid(labelspace,row=2,column=1)
    
    nextbutton <- tkbutton(ttcon,image=imageconvert,text="RUN >>",command=gltgp)
    tkgrid(nextbutton,row=3,column=2,sticky="e")      
    
    tkfocus(ttcon)
    
  }


  
  
     # -----------------------------------------------
     # HELLO | SIMIBD STARTS HERE
     #------------------------------------

   SimIBD <- function(){

     #this is going to be complicated



     runibd <- function(){

       if(tclvalue(sep1)=="White space")
         tclvalue(sep1) <- " "
       if(tclvalue(sep2)=="White space")
        tclvalue(sep2) <- " " 
       
       vect1<- c()
       #vect2 <- c()
       
       vec1<-unlist(strsplit(tclvalue(nall),","))
       for (i in 1:length(vec1)) vect1[i]<-as.numeric(vec1[i])
       
       #vec2<-unlist(strsplit(tclvalue(param),","))
       #for (i in 1:length(vec2)) vect2[[i]]<-as.numeric(vec2[i])
       
       if(tclvalue(nloc)!=length(vect1)){
         tkmessageBox(message="Number of locus must be equal to the length of the number of alleles per locus",icon="error",type="ok",parent=tt)
       }
       else{
         tttry <- tktoplevel(parent=.TkRoot)
         tkgrab(tttry)
         tkwm.geometry(tttry, "+200+200")
         tkwm.title(tttry,"wait")
         warn<-tklabel(tttry,image=imagepleasewait)
         tkpack(warn)
         tkfocus(tttry)
         
         print("Starting...")
         Sys.sleep(0.1)
         if(length(globalcoordinates)==1 || tclvalue(prevcoord)==0)
           idb.dataset <<- try(simdata( nindiv=as.numeric(tclvalue(nindiv)),coord.lim=c(as.numeric(tclvalue(absmin)),as.numeric(tclvalue(absmax)),as.numeric(tclvalue(ordmin)),as.numeric(tclvalue(ordmax))),number.nuclei=as.numeric(tclvalue(nuclei)),allele.numbers=vect1,IBD=TRUE,beta=as.numeric(tclvalue(beta)),npop=as.numeric(tclvalue(npop)),give.freq.grid=as.logical(tclvalue(freq.grid)),give.tess.grid=as.logical(tclvalue(tess.grid)),npix=c(as.numeric(tclvalue(npixh)),as.numeric(tclvalue(npixv))),comp.Fst=as.logical(tclvalue(comp.Fst)),comp.Dsigma2=as.logical(tclvalue(comp.Dsigma2)),comp.diff=as.logical(tclvalue(comp.height)),width=as.numeric(tclvalue(hwidth)),plot.pairs.borders=as.logical(tclvalue(plot.pairs.borders))),silent=TRUE)
         else
           idb.dataset <<- try(simdata( nindiv=as.numeric(tclvalue(nindiv)),coord.indiv=globalcoordinates,coord.lim=c(as.numeric(tclvalue(absmin)),as.numeric(tclvalue(absmax)),as.numeric(tclvalue(ordmin)),as.numeric(tclvalue(ordmax))),number.nuclei=as.numeric(tclvalue(nuclei)),allele.numbers=vect1,IBD=TRUE,beta=as.numeric(tclvalue(beta)),npop=as.numeric(tclvalue(npop)),give.freq.grid=as.logical(tclvalue(freq.grid)),give.tess.grid=as.logical(tclvalue(tess.grid)),npix=c(as.numeric(tclvalue(npixh)),as.numeric(tclvalue(npixv))),comp.Fst=as.logical(tclvalue(comp.Fst)),comp.Dsigma2=as.logical(tclvalue(comp.Dsigma2)),comp.diff=as.logical(tclvalue(comp.height)),width=as.numeric(tclvalue(hwidth)),plot.pairs.borders=as.logical(tclvalue(plot.pairs.borders))),silent=TRUE)
         tkdestroy(tttry)
         print("Done.")
         if (class(idb.dataset) == "try-error"){
           tkmessageBox(message=idb.dataset,icon="error",type="ok",parent=tt)              
           idb.dataset <<- 0
         }
         else{
         #  print("FST,FIS,FIT")
         #  print(idb.dataset$Fst)
         #  print(idb.dataset$Fis)
         #  print(idb.dataset$Fit)
           
         #  print("DSIGMA")
         #  print(idb.dataset$Dsigma2)

         #  print("diff.B,diff.W")
         #  print(idb.dataset$diff.B)
         #  print(idb.dataset$diff.W)
           
           tkmessageBox(message="Terminated with success",type="ok",parent=tt)
           globalcoordinates <<- t(idb.dataset$coord.indiv)
           #tkconfigure(coordownlabel.widget,text="Coordinates:    Simulated IBD data loaded")
           tclvalue(labelcoordtext) <- "Coordinates:    Simulated IBD data loaded"
           globalgenotypes <<- idb.dataset$genotypes
           #tkconfigure(genodownlabel.widget,text="Genotype:       Simulated IBD data loaded")
           tclvalue(labelgenotext)  <- "Genotype:       Simulated IBD data loaded"
           globalallele.numbers <<- idb.dataset$allele.numbers
           
           if(tclvalue(save)==1){
             
             auxcoord <-tclVar()
             tclvalue(auxcoord) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save coordinates file to:"))
             auxgen <-tclVar()
             tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save genotypes file to:"))
             
             
             write.table(t(idb.dataset$coord.indiv),file=tclvalue(auxcoord),sep=tclvalue(sep1),row.names=FALSE,col.names=FALSE)
             write.table(idb.dataset$genotypes,file = tclvalue(auxgen),sep=tclvalue(sep2),row.names=FALSE,col.names=FALSE)
           }
         }
         
       }
       if(tclvalue(sep1)==" ")
         tclvalue(sep1) <- "White space"
       if(tclvalue(sep2)==" ")
         tclvalue(sep2) <- "White space"
       
     }
     
     #nindiv
     nindiv=tclVar("0")
     nindiv.widget <-tkentry(ttibd,width="20",textvariable=nindiv)
     nindivlabel.widget <- tklabel(ttibd,text="Number of individuals:")
     
     #coord.lim

     coordxlabel.widget <- tklabel(ttibd,text="Limits of geographical domain:")
     
     absmin <- tclVar(0)
     absmax <- tclVar(1)
     absmin.widget <-tkentry(ttibd,width="8",textvariable=absmin)
     absmax.widget <-tkentry(ttibd,width="8",textvariable=absmax)
     abslabel.widget <- tklabel(ttibd,text="   abs (min|max) :")
     
     ordmin <- tclVar(0)
     ordmax <- tclVar(1)
     ordmin.widget <-tkentry(ttibd,width="8",textvariable=ordmin)
     ordmax.widget <-tkentry(ttibd,width="8",textvariable=ordmax)
     ordlabel.widget <- tklabel(ttibd,text="   ord (min|max) :")

     #rate
     #rate=tclVar("0")
     #rate.widget <-tkentry(ttibd,width="20",textvariable=rate)
     #ratelabel.widget <- tklabel(ttibd,text="Rate of the Poisson process:")

     # number.nuclei 

     nuclei=tclVar("0")
     nuclei.widget <-tkentry(ttibd,width="20",textvariable=nuclei)
     nucleilabel.widget <- tklabel(ttibd,text="Number of nuclei in tessellation:")
     
     #coord.nuclei

     #cuclei <- tclVar(0)
     #cuclei.widget <-tkentry(ttibd,width="20",textvariable=cuclei)
     #cucleilabel.widget <- tklabel(ttibd,text="Coordinates of nuclei of Voronoi tessellation:")

     #color.nuclei

     #coclei <- tclVar(0)
     #coclei.widget <-tkentry(ttibd,width="20",textvariable=coclei)
     #cocleilabel.widget <- tklabel(ttibd,text="Population labels of the nuclei (E.g: n,m,..):")

     #allele.numbers

     nloc <- tclVar(0)
     nloc.widget <-tkentry(ttibd,width="20",textvariable=nloc)
     nloclabel.widget <- tklabel(ttibd,text="Number of loci:")
     
     
     #allele.numbers

     nall <- tclVar()
     nall.widget <-tkentry(ttibd,width="20",textvariable=nall)
     nalllabel.widget <- tklabel(ttibd,text="Number of alleles per locus (E.g: 10,3,8,..):")
     
     #model

#     modellabel.widget <- tklabel(ttibd,text="Spatial covariance model:")
#     model <- tclVar("besset")
#     wmodel <- .Tk.subwin(ttibd)
#     modeloptionmenu.widget <-  tcl("tk_optionMenu",wmodel,model,"besset","cauchy","cauchytbm","circular","constant","cone","cubic","cutoff","dampedcosine","exponential","fractalB","FD","fractgauss","gauss","gencauchy","gengneiting","gneiting","hyperbolic","iacocesare","lgd1","mastein","nsst","nsst2","nugget","penta","power","qexponential","spherical","stable","Stein","steinst1","wave","whittlematern")

     #param

#     param <- tclVar("")
#     param.widget <-tkentry(ttibd,width="20",textvariable=param)
#     paramlabel.widget <- tklabel(ttibd,text="Parameters (E.g:mean, variance, nugget, scale, ...):")

     #Beta
     beta <- tclVar("")
     beta.widget <-tkentry(ttibd,width="20",textvariable=beta)
     betalabel.widget <- tklabel(ttibd,text="Spatial correlation parameter for frequencies:")
     
     #npop
     npop <- tclVar("")
     npop.widget <-tkentry(ttibd,width="20",textvariable=npop)
     npoplabel.widget <- tklabel(ttibd,text="Number of populations:")

     #freq.grid
     freq.gridlabel.widget <- tklabel(ttibd,text="Return frequencies on grid:")
     freq.grid <- tclVar("FALSE")
     wfreq.grid <- .Tk.subwin(ttibd)
     freq.gridoptionmenu.widget <-  tcl("tk_optionMenu",wfreq.grid,freq.grid,"FALSE","TRUE")

     #tess.grid
     tess.gridlabel.widget <- tklabel(ttibd,text="Return population membership on grid:")
     tess.grid <- tclVar("FALSE")
     wtess.grid <- .Tk.subwin(ttibd)
     tess.gridoptionmenu.widget <-  tcl("tk_optionMenu", wtess.grid,tess.grid,"FALSE","TRUE")

     #npix

     npixh <- tclVar(50)
     npixv <- tclVar(50)
     npixh.widget <-tkentry(ttibd,width="8",textvariable=npixh)
     npixv.widget <-tkentry(ttibd,width="8",textvariable=npixv)
     npixlabel.widget <- tklabel(ttibd,text="Number of pixels for representation (hor|ver):")
     
     #comp.Fst
     comp.Fstlabel.widget <- tklabel(ttibd,text="Compute F statistics:")
     comp.Fst <- tclVar("FALSE")
     wcomp.Fst <- .Tk.subwin(ttibd)
     comp.Fstoptionmenu.widget <-  tcl("tk_optionMenu", wcomp.Fst,comp.Fst,"FALSE","TRUE")
     
     #comp.Dsigma2
     comp.Dsigma2label.widget <- tklabel(ttibd,text="Compute IBD index Dsigma2:")
     comp.Dsigma2 <- tclVar("FALSE")
     wcomp.Dsigma2 <- .Tk.subwin(ttibd)
     comp.Dsigma2optionmenu.widget <-  tcl("tk_optionMenu",wcomp.Dsigma2,comp.Dsigma2,"FALSE","TRUE")

     #comp.height
     comp.heightlabel.widget <- tklabel(ttibd,text="Index of variability of allele freq.:")
     comp.height <- tclVar("FALSE")
     wcomp.height <- .Tk.subwin(ttibd)
     comp.heightoptionmenu.widget <-  tcl("tk_optionMenu", wcomp.height,comp.height,"FALSE","TRUE")

     
     #width

     hwidth <- tclVar(0.1)
     hwidth.widget <-tkentry(ttibd,width="20",textvariable=hwidth)
     hwidthlabel.widget <- tklabel(ttibd,text="Width around the barriers:")

     #plot.pairs.borders
#     plot.pairs.borderslabel.widget <- tklabel(ttibd,text="Plot pairs of individuals:")
     plot.pairs.borders <- tclVar("FALSE")
#     wplot.pairs.borders <- .Tk.subwin(ttibd)
#     plot.pairs.bordersoptionmenu.widget <-  tcl("tk_optionMenu", wplot.pairs.borders,plot.pairs.borders,"FALSE","TRUE")

     prevcoord <- tclVar(0)
     
     prevcoordlabel.widget <- tklabel(ttibd,text="Use loaded coordinates file:")
     prevcoord.widget <- tkcheckbutton(ttibd,variable=save,onvalue=1, offvalue=0)


     
     #savefiles
     save <- tclVar(0)
     savelabel.widget <- tklabel(ttibd,text="Save coordinates and genotypes files:")
     save.widget <- tkcheckbutton(ttibd,variable=save,onvalue=1, offvalue=0)
     
     tkgrid(nindivlabel.widget,row=1,column=1,sticky="w")
     tkgrid(nindiv.widget,row=1,column=2,columnspan=2,sticky="w")
     tkgrid(coordxlabel.widget,row=2,column=1,sticky="w")
     tkgrid(abslabel.widget,row=3,column=1,sticky="w")
     tkgrid(absmin.widget,row=3,column=2,sticky="w")
     tkgrid(absmax.widget,row=3,column=3,sticky="w")
     tkgrid(ordlabel.widget,row=4,column=1,sticky="w")
     tkgrid(ordmin.widget,row=4,column=2,sticky="w")
     tkgrid(ordmax.widget,row=4,column=3,sticky="w")
     #tkgrid(ratelabel.widget,row=5,column=1,sticky="w")
     #tkgrid(rate.widget,row=5,column=2,columnspan=2,sticky="w")
     tkgrid(nucleilabel.widget,row=5,column=1,sticky="w")
     tkgrid(nuclei.widget,row=5,column=2,columnspan=2,sticky="w")
     tkgrid(nloclabel.widget,row=6,column=1,sticky="w")
     tkgrid(nloc.widget,row=6,column=2,columnspan=2,sticky="w")
     tkgrid(nalllabel.widget,row=7,column=1,sticky="w")
     tkgrid(nall.widget,row=7,column=2,columnspan=2,sticky="w")
##     tkgrid(modellabel.widget,row=8,column=1,sticky="w")
##     tkgrid(wmodel,row=8,column=2,columnspan=2,sticky="w")
##     tkgrid(paramlabel.widget,row=9,column=1,sticky="w")
##     tkgrid(param.widget,row=9,column=2,columnspan=2,sticky="w")
     tkgrid(betalabel.widget,row=9,column=1,columnspan=2,sticky="w")
     tkgrid(beta.widget,row=9,column=2,columnspan=2,sticky="w")
     tkgrid(npoplabel.widget,row=10,column=1,sticky="w")
     tkgrid(npop.widget,row=10,column=2,columnspan=2,sticky="w")
     tkgrid(freq.gridlabel.widget,row=11,column=1,sticky="w")
     tkgrid(wfreq.grid,row=11,column=2,columnspan=2,sticky="w")
     tkgrid(tess.gridlabel.widget,row=12,column=1,sticky="w")
     tkgrid(wtess.grid,row=12,column=2,columnspan=2,sticky="w")
     tkgrid(npixlabel.widget,row=13,column=1,sticky="w")
     tkgrid(npixh.widget,row=13,column=2,sticky="w")
     tkgrid(npixv.widget,row=13,column=3,sticky="w")
     tkgrid(comp.Fstlabel.widget,row=14,column=1,sticky="w")
     tkgrid(wcomp.Fst,row=14,column=2,columnspan=2,sticky="w")
     tkgrid(comp.Dsigma2label.widget,row=15,column=1,sticky="w")
     tkgrid(wcomp.Dsigma2,row=15,column=2,columnspan=2,sticky="w")
     tkgrid(comp.heightlabel.widget,row=16,column=1,sticky="w")
     tkgrid(wcomp.height,row=16,column=2,columnspan=2,sticky="w")
     tkgrid(hwidthlabel.widget,row=17,column=1,sticky="w")
     tkgrid(hwidth.widget,row=17,column=2,columnspan=2,sticky="w")
     #tkgrid(plot.pairs.borderslabel.widget,row=18,column=1,sticky="w")
     #tkgrid(wplot.pairs.borders,row=18,column=2,columnspan=2,sticky="w")
     tkgrid(prevcoordlabel.widget,row=19,column=1,sticky="w")
     tkgrid(prevcoord.widget,row=19,column=2,sticky="w")
     tkgrid(savelabel.widget,row=20,column=1,sticky="w")
     tkgrid(save.widget,row=20,column=2,columnspan=2,sticky="w")

     labelspace <-tklabel(ttibd,text=" ")
     tkgrid(labelspace,row=21,column=1)
     
     nextbutton <- tkbutton(ttibd,image=imagerun2,text="RUN >>",command=runibd)
     tkgrid(nextbutton,row=22,column=2,columnspan=2,sticky="e")
   }


     # -----------------------------------------------
     # HELLO | NULLIFY STARTS HERE
     #------------------------------------

  Nullify <- function(){

    ttnul<-tktoplevel()
    
    tkwm.title(ttnul,"Simulate genotypes with null alleles from loaded dataset")
    
    gltgp <- function(){
      
      
      if(length(globalgenotypes)==0){
        tkmessageBox(message="You must define genotypes file",icon="error",type="ok")
        
      }
      else{
        
        tttry <- tktoplevel(parent=.TkRoot)
        tkgrab(tttry)
        tkwm.geometry(tttry, "+200+200")
        tkwm.title(tttry,"wait")
        warn<-tklabel(tttry,image=imagepleasewait)
        tkpack(warn)
        tkfocus(tttry)
        
        print("Starting...")
        Sys.sleep(0.1)
        
        err <- try(nullify(genotypes=globalcoordinates,nall.null=as.integer(tclvalue(nall)),nloc.null=as.integer(tclvalue(nloc))),silent=TRUE)
        tkdestroy(tttry)
        print("Done.")
        if (class(err) == "try-error") 
          tkmessageBox(message=err,icon="error",type="ok",parent=tt)              
        else{
          tkmessageBox(message="Terminated with success",type="ok",parent=tt)
          globalgenotypes <<- err$genotypes
          if(tclvalue(save)==1){
            if(tclvalue(sep2)=="White space")
              tclvalue(sep2) <- " " 
            auxgen <-tclVar()
            tclvalue(auxgen) <- tclvalue(tkgetSaveFile(filetypes="{{All files} *}",initialdir=tclvalue(outputdir),title="Save genotypes file to:"))
            write.table(err$genotypes,file = tclvalue(auxgen),sep=tclvalue(sep2),row.names=FALSE,col.names=FALSE)
            if(tclvalue(sep2)==" ")
              tclvalue(sep2) <- "White space"
          }
          tclvalue(labelgenotext) <- paste(tclvalue(labelgenotext),"with null alleles")
        }
       
      }
    }
    #nall.null
    nall<- tclVar(1)
    
    nallentry.widget <- tkentry(ttnul,textvariable=nall,width=15)
    nalllabel.widget <- tklabel(ttnul,text="Number of null alleles on each locus:")
    
    tkgrid(nalllabel.widget,row=1,column=1,sticky="w")
    tkgrid(nallentry.widget,row=1,column=2,sticky="w")
    
    #nloc.null
    nloc<- tclVar(1)
    
    nlocentry.widget <- tkentry(ttnul,textvariable=nloc,width=15)
    nloclabel.widget <- tklabel(ttnul,text="Number of loci with null alleles:")
    
    tkgrid(nloclabel.widget,row=2,column=1,sticky="w")
    tkgrid(nlocentry.widget,row=2,column=2,sticky="w")

    save <- tclVar(0)
    savelabel.widget <- tklabel(ttnul,text="Save genotypes files:")
    save.widget <- tkcheckbutton(ttnul,variable=save,onvalue=1, offvalue=0)

    tkgrid(savelabel.widget,row=3,column=1,sticky="w")
    tkgrid(save.widget,row=3,column=2,sticky="w")
    labelspace <-tklabel(ttnul,text=" ")
    tkgrid(labelspace,row=4,column=1)
    
    nextbutton <- tkbutton(ttnul,image=imageok,text="RUN >>",command=gltgp)
    tkgrid(nextbutton,row=5,column=2,sticky="e")      
    
    tkfocus(ttnul)

  }
     # -----------------------------------------------
     # HELLO | RESET STARTS HERE
     #------------------------------------

  Reset <- function(){
    auxblink<<-0
    idb.dataset <<- 0
    globalcoordinates <<- 0
    globalgenotypes <<- 0
    allele.numbers <<- 0
    
    tclvalue(coordinatesfile) <<- ""
    tclvalue(genotypefile) <<- "" 
    tclvalue(outputdir) <<- ""
    tclvalue(advanced) <<- 0
    
    tkconfigure(coordownlabel.widget,text="")
    tkconfigure(genodownlabel.widget,text="")


    configure()
    run()
    postproc()
    plot()
    Gfstat()
    simFModel()
    SimIBD()
    
    
  }

     # -----------------------------------------------
     # HELLO | ShowTEXT STARTS HERE
     #------------------------------------

  Showtext <- function(filename){
    
    #tttry <- tktoplevel()
    #tkgrab(tttry)
    #tkwm.geometry(tttry, "+200+200")
    #tkwm.title(tttry,"wait")
    #warn<-tklabel(tttry,image=imagepleasewait)
    #tkpack(warn)
    #tkfocus(tttry)
    #Sys.sleep(0.1)
    
    
    file <- try(readLines(paste(tclvalue(outputdir),filename,sep="")),silent=TRUE)
    
    if (class(file) == "try-error"){ 
     # tkdestroy(tttry)
      tkmessageBox(message="File hasn't been created or bad output path",type="ok",parent=tt)              
    }
    else{
      #print(file)
      
      tttext <- tktoplevel(parent=.TkRoot)
      tkwm.title(tttext,filename)
      yscr <- tkscrollbar(tttext, repeatinterval=5,command=function(...)tkyview(txt,...))
      xscr <- tkscrollbar(tttext, repeatinterval=5,orient="horizontal",command=function(...)tkxview(txt,...))
      txt <- tktext(tttext,font="courier", wrap="none",yscrollcommand=function(...)tkset(yscr,...),xscrollcommand=function(...)tkset(xscr,...))
      #print(file)
      #print(length(file))
      for (i in 1:length(file)){
        tkinsert(txt,"end",file[i])
        tkinsert(txt,"end","\n")
      }
      #tkdestroy(tttry)
      tkgrid( txt,row=1,column=1)
      tkgrid( yscr,row=1,column=2, sticky="ns")
      tkgrid( xscr,row=2,column=1, sticky="we")
    }
  }

  
  initialimage()    
  
  topMenu <- tkmenu(tt)
  
  tkconfigure(tt,menu=topMenu)
  fileMenu <- tkmenu(topMenu,tearoff=FALSE)
  coordinatesMenu <- tkmenu(topMenu,tearoff=FALSE)
  genotypesMenu <- tkmenu(topMenu,tearoff=FALSE)
  toolsMenu <- tkmenu(topMenu,tearoff=FALSE)
  helpMenu <- tkmenu(topMenu,tearoff=FALSE)
  
  tkadd(fileMenu,"checkbutton",label="Advanced Options",variable=advanced,selectcolor="blue",onvalue=1,offvalue=0,command=function() fadvanced())
  tkadd(coordinatesMenu,"radiobutton",label="White space",value="",variable=sep1,selectcolor="blue")
  tkadd(coordinatesMenu,"radiobutton",label=",",value=",",variable=sep1,selectcolor="blue")
  tkadd(coordinatesMenu,"radiobutton",label=";",value=";",variable=sep1,selectcolor="blue")
  tkadd(genotypesMenu,"radiobutton",label="White space",value="",variable=sep2,selectcolor="blue")
  tkadd(genotypesMenu,"radiobutton",label=",",value=",",variable=sep2,selectcolor="blue")
  tkadd(genotypesMenu,"radiobutton",label=";",value=";",variable=sep2,selectcolor="blue")
  tkadd(fileMenu,"separator")
  tkadd(fileMenu,"cascade",label="Coordinates File Values Separator",menu=coordinatesMenu)
  tkadd(fileMenu,"cascade",label="Genotypes File Values Separator",menu=genotypesMenu)
  tkadd(toolsMenu,"command",label="Convert to Genepop Files",command=function() Convert())
  tkadd(toolsMenu,"command",label="Simulate null alleles",state="disabled",command=function() Nullify())
  tkadd(fileMenu,"separator")
  tkadd(fileMenu,"command",label="Reset Values",command=function() Reset())
  tkadd(fileMenu,"separator")
  tkadd(fileMenu,"command",label="Quit",command=function() tkdestroy(tt))
  tkadd(topMenu,"cascade",label="Menu",menu=fileMenu)
  tkadd(topMenu,"cascade",label="Tools",menu=toolsMenu)
  tkadd(helpMenu,"command",label="Help",command=function() helpWindow())
  tkadd(helpMenu,"command",label="Credits",command=function() creditsWindow())
  tkadd(topMenu,"cascade",label="Help",menu=helpMenu)
  
  labelinference <-tklabel(ttpan,text="-Inference-",font="*-Times-bold-i-normal--20-*",foreground="blue")
  labelsimulation <-tklabel(ttpan,text="-Simulation-",state="disabled",font="*-Times-bold-i-normal--20-*",foreground="blue")
  labelspace <-tklabel(ttpan,text=" ")
  buttonconf <- tkbutton(ttpan,image=imageconfigure,text="Configure",command=function() {configure();tkgrid.remove(ttinit);tkgrid.remove(ttsimf);tkgrid.remove(ttibd);tkgrid.remove(ttfstat);tkgrid.remove(ttplot2);tkgrid.remove(ttrun);tkgrid.remove(ttpost);tkgrid.remove(ttplot);tkgrid(ttconf,row=1,column=2,sticky="we")})
  buttonfstat <- tkbutton(ttpan,image=imagefstat,text="Fstat",command=function(){ Gfstat();tkgrid.remove(ttinit);tkgrid.remove(ttsimf);tkgrid.remove(ttconf);tkgrid.remove(ttrun);tkgrid.remove(ttibd);tkgrid.remove(ttplot2);tkgrid.remove(ttpost);tkgrid.remove(ttplot);tkgrid(ttfstat,row=1,column=2,sticky="we")})
  buttonrun <- tkbutton(ttpan,image=imagerun,text="Run",command=function(){ run();tkgrid.remove(ttinit);tkgrid.remove(ttsimf);tkgrid.remove(ttconf);tkgrid.remove(ttfstat);tkgrid.remove(ttibd);tkgrid.remove(ttplot2);tkgrid.remove(ttpost);tkgrid.remove(ttplot);tkgrid(ttrun,row=1,column=2,sticky="we")})
  buttonpostprocess <- tkbutton(ttpan,image=imagepostprocess,text="Postprocess",command=function(){ postproc();tkgrid.remove(ttinit);tkgrid.remove(ttsimf);tkgrid.remove(ttibd);tkgrid.remove(ttfstat);tkgrid.remove(ttplot2);tkgrid.remove(ttconf);tkgrid.remove(ttrun);tkgrid.remove(ttplot);tkgrid(ttpost,row=1,column=2,sticky="we")})
  buttonsimfmodel <- tkbutton(ttpan,image=imagefmodel,text="F-model",state="disabled",command=function(){ SimnonIBD();tkgrid.remove(ttinit);tkgrid.remove(ttconf);tkgrid.remove(ttibd);tkgrid.remove(ttfstat);tkgrid.remove(ttplot2);tkgrid.remove(ttrun);tkgrid.remove(ttplot);tkgrid(ttsimf,row=1,column=2,sticky="we")})
  buttonplot <- tkbutton(ttpan,image=imageplot,text="Plot",command=function(){ plot();tkgrid.remove(ttinit);tkgrid.remove(ttconf);tkgrid.remove(ttsimf);tkgrid.remove(ttibd);tkgrid.remove(ttfstat);tkgrid.remove(ttplot2);tkgrid.remove(ttrun);tkgrid.remove(ttpost);tkgrid(ttplot,row=1,column=2,sticky="we")})
  buttonibd <- tkbutton(ttpan,image=imageibd,text="IBD",state="disabled",command=function(){ SimIBD();tkgrid.remove(ttinit);tkgrid.remove(ttconf);tkgrid.remove(ttsimf);tkgrid.remove(ttplot);tkgrid.remove(ttplot2);tkgrid.remove(ttfstat);tkgrid.remove(ttrun);tkgrid.remove(ttpost);tkgrid(ttibd,row=1,column=2,sticky="we")})
  buttonplot2 <- tkbutton(ttpan,image=imageplot,text="Plot2",state="disabled",command=function(){ plot2();tkgrid.remove(ttinit);tkgrid.remove(ttconf);tkgrid.remove(ttsimf);tkgrid.remove(ttibd);tkgrid.remove(ttplot);tkgrid.remove(ttfstat);tkgrid.remove(ttrun);tkgrid.remove(ttpost);tkgrid(ttplot2,row=1,column=2,sticky="we")})
  
  tkgrid(labelinference,row=1,column=1,sticky="w")
  tkgrid(buttonconf,row=2,column=1,sticky="we")
  tkgrid(buttonrun,row=3,column=1,sticky="we")
  tkgrid(buttonpostprocess,row=4,column=1,sticky="we")
  tkgrid(buttonplot,row=5,column=1,sticky="we")
  #tkgrid(buttonfstat,row=6,column=1,sticky="we")
  tkgrid(labelspace,row=6,column=1,sticky="w")
  tkgrid(labelsimulation,row=7,column=1,sticky="w")
  tkgrid(buttonsimfmodel,row=8,column=1,sticky="we")
  tkgrid(buttonibd,row=9,column=1,sticky="we")
  tkgrid(buttonplot2,row=10,column=1,sticky="we")


  
  coordownlabel.widget <- tklabel(tt,textvariable=labelcoordtext,foreground="blue")
  genodownlabel.widget <- tklabel(tt,textvariable=labelgenotext,foreground="blue")
  
  auxblink <- 1
  extralabel.widget <- tklabel(tt,text="Please configure output dir",foreground="blue")
  blink <- function(){
    if (auxblink==1){
      try(tkconfigure(extralabel.widget,text=""),silent = TRUE)
      auxblink<<-0
    }
    else if (auxblink==0){
      try(tkconfigure(extralabel.widget,text="Please configure output dir"),silent = TRUE)
      auxblink<<-1
    }
    try(tcl("after",1000,blink),silent = TRUE)
  }
  blink()
  
  tkgrid(ttpan,row=1,column=1,sticky="we")  
  tkgrid(ttinit,row=1,column=2,sticky="e")  
  tkgrid(coordownlabel.widget,row=2,column=1,columnspan=2,sticky="w",padx=3)
  tkgrid(genodownlabel.widget,row=3,column=1,columnspan=2,sticky="w",padx=3)
  tkgrid(extralabel.widget,row=4,column=1,columnspan=2,sticky="w",padx=3)
  tkgrid.columnconfigure(tt,1,minsize=200)
  tkgrid.columnconfigure(tt,2,minsize=500)
  tkgrid.rowconfigure(tt,1,minsize=450)
  
}

