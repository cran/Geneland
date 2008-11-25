      subroutine mcmcgld(s,z,nal,ploidy,dom,path,nchpath,lambdamax,
     &     dt,shape1,shape2,nit,thinning,filtna,
     &     nindiv,nloc,nloc2,nalmax,npp,nppmax,
     &     npop,npopmin,npopmax,
     &     t,ttmp,u,utmp,c,ctmp,f,ftmp,fa,drift,drifttmp,
     &     indcell,indcelltmp,
     &     distcell,distcelltmp,xlim,ylim,n,ntmp,a,ptmp,
     &     cellpop,listcell,fmodel,kfix,spatial,jcf,
     &     y,fcy,ofiles,pudcel,missloc) 
      implicit none 

*     data
      integer nindiv,nloc,nloc2,
     &     nal,nalmax,z,ploidy,dom,jcf,nchpath,missloc
      double precision s

*     hyper parameters
      double precision lambdamax,dt,shape1,shape2

*     parameters
      integer npp,nppmax,npop,npopmin,npopmax,c,ctmp
      double precision lambda,u,utmp,f,t,fa,drift,ftmp,drifttmp,alpha

*     computing options
      integer nit,thinning,fmodel,kfix,spatial,filtna,ofiles,nudcel
      double precision du,pudcel

*     variables de travail
      integer iloc,iindiv,iit,ipp,ipop,ial,indcell,indcelltmp,
     &     n,cellpop,listcell,cellpophost,ntmp,nn,y,iud,nud
      double precision ptmp,xlim,ylim,ggrunif,rpostlamb,
     &     distcell,distcelltmp,a,ttmp,lpp,ll,fcy,pct
      character*255 path,filef,filenpp,filelambda,filenpop,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet,filesize
 
*     dimensionnement 
      dimension s(2,nindiv),t(2,nindiv),z(nindiv,nloc2),
     &     u(2,nppmax),utmp(2,nppmax),c(nppmax),ctmp(nppmax),
     &     f(npopmax,nloc,nalmax),xlim(2),ylim(2),
     &     nal(nloc),indcell(nindiv),indcelltmp(nindiv),
     &     distcell(nindiv),ttmp(2,nindiv),
     &     distcelltmp(nindiv),n(npopmax,nloc,nalmax),
     &     ntmp(npopmax,nloc,nalmax),
     &     a(nalmax),ptmp(nalmax),
     &     fa(nloc,nalmax),drift(npopmax),
     &     ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     cellpop(nppmax),listcell(nppmax),cellpophost(nppmax),
     &     y(nindiv,nloc2),fcy(nalmax,2),ofiles(12),missloc(nindiv,nloc)


      call intpr('***************************************',-1,0,0)
      call intpr('***    Starting MCMC simulation     ***',-1,0,0)
      call intpr('***************************************',-1,0,0)

*     init RNG
      call rndstart()
   
 
*     look for smallest rectangle enclosing the spatial domain
      call limit(nindiv,s,xlim,ylim,dt)


*     Ouverture des fichiers pour l'ecriture des sorties
      filelambda = path(1:nchpath) // "Poisson.process.rate.txt"
      filenpp = path(1:nchpath) // "nuclei.numbers.txt"
      filenpop = path(1:nchpath) // "populations.numbers.txt"
      fileu = path(1:nchpath) // "coord.nuclei.txt"
      filec = path(1:nchpath) // "color.nuclei.txt"
      filef = path(1:nchpath) // "frequencies.txt"
      filefa = path(1:nchpath) // "ancestral.frequencies.txt"
      filedrift = path(1:nchpath) // "drifts.txt"
      filelpp = path(1:nchpath) // "log.posterior.density.txt"
      filell = path(1:nchpath) // "log.likelihood.txt"
      filet = path(1:nchpath) // "hidden.coord.txt"
      filesize = path(1:nchpath) // "size.pop.txt"

c 100  format(f7.3,' %')
 1000 format (2(1x,e15.8,1x))
 2000 format (300(1x,e15.8,1x))

      if(ofiles(1) .eq.1) then
         open(9,file=filelambda)
      endif
      if(ofiles(2) .eq.1) then
         open(10,file=filenpp)
      endif
      if(ofiles(3) .eq.1) then
         open(11,file=filenpop)
      endif   
      if(ofiles(4) .eq.1) then
         open(12,file=fileu)
      endif
      if(ofiles(5) .eq.1) then
         open(13,file=filec)
      endif
      if(ofiles(6) .eq.1) then
         open(14,file=filef)
      endif
      if(ofiles(7) .eq.1) then
         open(15,file=filefa) 
      endif
      if(ofiles(8) .eq.1) then
         open(16,file=filedrift)
      endif
      if(ofiles(9) .eq.1) then
         open(17,file=filelpp)
      endif
      if(ofiles(10) .eq.1) then
         open(18,file=filell) 
      endif
      if(ofiles(1) .eq.1) then
         open(19,file=filet) 
      endif
      if(ofiles(12) .eq.1) then
         open(20,file=filesize)
      endif

       call intpr('Output files opened:',-1,0,0)

 
************************
*     Initialization
************************
*     parameter for the Dirichlet model for allele freq. 
*     (not used if Falush .eq.1)
      alpha = 1

c$$$      npp = 2
c$$$      u(1,1) = .25
c$$$      u(2,1) = .5
c$$$      u(1,2) = .75
c$$$      u(2,2) = .5           
c$$$      c(1) = 1
c$$$      c(2) = 2

      
c$$$      npp = 4
c$$$      u(1,1) = 4
c$$$      u(2,1) = 3
c$$$      u(1,2) = 8
c$$$      u(2,2) = 3
c$$$      u(1,3) = 4
c$$$      u(2,3) = 5
c$$$      u(1,4) = 8
c$$$      u(2,4) = 5
c$$$      c(1) = 1 
c$$$      c(2) = 2 
c$$$      c(3) = 2
c$$$      c(4) = 2

c$$$      npp = 20
c$$$      do ipp = 1,(npp/2)
c$$$         u(1,ipp) = 4
c$$$         u(2,ipp) = 2.5+2.5*dble(ipp)/(npp/2)
c$$$      enddo
c$$$      do ipp = (npp/2)+1,npp
c$$$         u(1,ipp) = 8
c$$$         u(2,ipp) = 2.5+2.5*(dble(ipp)-(npp/2))/(npp/2)
c$$$      enddo
c$$$      do ipp = 1,ipp
c$$$         c(ipp) = 1
c$$$      enddo

      du = dsqrt((xlim(2)-xlim(1))*(ylim(2)-ylim(1))/dble(nindiv))
      lambda = lambdamax*ggrunif(0.d0,1.d0) 


      if(spatial .eq. 1) then
         npp = 1 + idint(dint(lambda))
      else
         npp = nindiv
      endif
      call intpr('npp initialised:',-1,0,0)


      if(spatial .eq. 1) then 
         call rprioru(npp,nppmax,xlim,ylim,u)
      else 
         do iindiv=1,nindiv
            u(1,iindiv) = s(1,iindiv)
            u(2,iindiv) = s(2,iindiv) 
         enddo
      endif 
      call intpr('u initialised:',-1,0,0)

      call rpriorc(npp,nppmax,npop,c)
      call intpr('c initialised:',-1,0,0)

      do iindiv=1,nindiv
         t(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         t(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
      enddo
      call intpr('t initialised:',-1,0,0)

      call calccell(nindiv,t,npp,nppmax,u,indcell,distcell)
      call intpr('Voronoi cells initialised:',-1,0,0)

      call rpriordrift(npop,npopmax,drift,fmodel,shape1,shape2)
      call intpr('drift initialised:',-1,0,0)

      call rpriorfa(nloc,nloc,nal,nalmax,fa,fmodel,ptmp)
      call intpr('fa initialised:',-1,0,0)

      call rpostf2(npop,npopmax,nloc,nloc,nal,nalmax,f,fa,drift,
     &     nindiv,nloc2,y,nppmax,c,indcell,n,a,ptmp,ploidy)
      call intpr('f initialised:',-1,0,0)
      


      call intpr('All parameters initialised:',-1,0,0)
     

c$$$      write(*,*) 'nindiv =', nindiv,"\n"
c$$$      write(*,*) 'nloc =',nloc   ,"\n"
c$$$      write(*,*) 'nal =',nal ,"\n"
c$$$      write(*,*) 'nalmax =', nalmax,"\n"
c$$$      write(*,*) 'z =', z,"\n"
c$$$      write(*,*) 's =', s,"\n"
c$$$      write(*,*) 'lambda=',lambda,"\n"
c$$$      write(*,*) 'lambdamax=',lambdamax  ,"\n"
c$$$      write(*,*) 'npp =', npp ,"\n"
c$$$      write(*,*) 'nppmax =', nppmax   ,"\n"  
c$$$      write(*,*) 'npopmin =',npopmin  ,"\n"
c$$$      write(*,*) 'npopmax =', npopmax ,"\n" 
c$$$      write(*,*) 'dt =',dt,"\n" 
c$$$      write(*,*) 'shape1 =',shape1,"\n" 
c$$$      write(*,*) 'shape2 =',shape2,"\n" 
c$$$      write(*,*) 'filtna =',filtna,"\n" 
c$$$      write(*,*) 'fmodel=',fmodel,"\n" 
c$$$      write(*,*) 'kfix =',kfix,"\n" 
c$$$      write(*,*) 'spatial =',spatial,"\n" 
c$$$      write(*,*) 'jcf =',jcf,"\n" 
c$$$      write(*,*) 'pudcel =',pudcel,"\n" 
c$$$      write(*,*) 'c =', c,"\n"
c$$$      write(*,*) 'ctmp =',ctmp ,"\n"
c$$$      write(*,*) 'u =',u  ,"\n"
c$$$      write(*,*) 'utmp =',utmp ,"\n"
c$$$      write(*,*) 'f =', f  ,"\n"
c$$$      write(*,*) 'ftmp =', ftmp,"\n"
c$$$      write(*,*) 'nit =', nit ,"\n" 
c$$$      write(*,*) 'thinning =',thinning  ,"\n"
c$$$      write(*,*) 'indcell =', indcell ,"\n"
c$$$      write(*,*) 'distcell =', distcell,"\n"
c$$$      write(*,*) 'indcelltmp =', indcelltmp ,"\n"
c$$$      write(*,*) 'distcelltmp =',distcelltmp ,"\n"
      write(*,*) 'xlim=', xlim 
      write(*,*) 'ylim=', ylim 
************************  
* mise a jour iterative  
************************

      call intpr('Percentage of computations:',-1,0,0)


c      write(*,*) 'Starting updates'
      do iit=1,nit
         if(mod(iit,thinning) .eq. 0) then

            pct = dble(iit)/dble(nit)*100.
            call dblepr('                     ',-1,pct,1)

            if(ofiles(1) .eq.1) then
               write(9,*) lambda
            endif
            if(ofiles(2) .eq.1) then
               write(10,*) npp
            endif
            if(ofiles(3) .eq.1) then
               write(11,*) npop
            endif
            if(ofiles(4) .eq.1) then
               write(12,1000) 
     &              (sngl(u(1,ipp)),sngl(u(2,ipp)), ipp=1,nppmax)
            endif
            if(ofiles(5) .eq.1) then
               do ipp=1,nppmax
                  write(13,*) c(ipp)
               enddo
            endif
            if(ofiles(6) .eq.1) then
               do iloc=1,nloc
                  do ial=1,nalmax
                     write(14,2000) (sngl(f(ipop,iloc,ial)),
     &                    ipop=1,npopmax)
                  enddo
               enddo  
            endif
            if(ofiles(7) .eq.1) then
               write(15,2000) 
     &              ((sngl(fa(iloc,ial)),ial=1,nalmax),iloc=1,nloc)
            endif
            if(ofiles(8) .eq.1) then
               write(16,2000) (sngl(drift(ipop)),ipop=1,npopmax)
            endif
            if(ofiles(9) .eq.1) then
               write(17,*) lpp(lambdamax,lambda,y,npop,npp,drift,f,fa,c,
     &              nppmax,nindiv,nloc2,npopmax,nloc,nalmax,
     &              indcell,nal,fmodel,xlim,ylim,shape1,shape2)    
            endif
            if(ofiles(10) .eq.1) then
               write(18,*) ll(y,nindiv,nloc,nloc2,npopmax,
     &              nalmax,nppmax,c,f,indcell)
            endif
            if(ofiles(11) .eq.1) then
               if(dt .gt. 1.d-300) then 
                  write(19,1000) (sngl(t(1,iindiv)),sngl(t(2,iindiv)),
     &                 iindiv=1,nindiv)
               endif
            endif
            if(ofiles(12) .eq.1) then
*     counting nb of individuals in each pop
               do ipop = 1,npopmax
                  nn = 0
                  do iindiv = 1,nindiv
                     if(c(indcell(iindiv)) .eq. ipop) then
                        nn = nn + 1
                     endif
                  enddo
                  write(20,*) nn
               enddo
            endif
         endif

*     nb of cells updated
         nudcel = max0(1,min0(idint(dint(dble(npp)*pudcel)),npp))

*******************
*     update lambda
c          write(*,*) 'update lambda'
          if(spatial .eq. 1) then
             lambda = rpostlamb(lambdamax,npp)
          endif

          if(fmodel .eq. 1) then 
*******************
*     update drift
c              write(*,*) 'update drift'
             call  upddrift(npop,npopmax,nloc,nalmax,nal,
     &            f,fa,drift,shape1,shape2)
*******************
*     update fa
c             write(*,*) 'update fa'
             call updfa(npop,npopmax,nloc,nalmax,nal,
     &            f,fa,drift) 
          endif
******************* 
*     update c and f 
c          nud =  1 + idint(ggrunif(0.d0,1.d0)*dble(npp))
c          write(*,*) 'nud=',nud
c          do iud = 1,nud
c     write(*,*) 'update f'  
             if(jcf .eq. 1) then 
c     write(*,*) 'update c and f' 
                if(fmodel .eq. 0) then
*     joint update of c anf f
                   call  udcf(npop,npopmax,f,nloc,nloc,nloc2,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,
     &                  c,ctmp,a,ptmp,ftmp,y,n,ntmp,ploidy,alpha,nudcel)
                else
                   call udcf2(npop,npopmax,f,fa,drift,nloc,nloc,nloc2,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,
     &                  c,ctmp,a,ptmp,ftmp,y,n,ntmp,ploidy,nudcel)
                endif
             else          
                call rpostf2(npop,npopmax,nloc,nloc,nal,nalmax,
     &               f,fa,drift,
     &               nindiv,nloc2,y,nppmax,c,indcell,
     &               n,a,ptmp,ploidy)
                call updc(npp,nppmax,c,ctmp,y,nindiv,nloc,
     &               nloc,nloc2,nalmax,npop,npopmax,f,indcell,ploidy,
     &               nudcel)
                
             endif
c          enddo

          
          if(spatial .eq. 1) then
*******************
*     update u et mise a jour de indcell et distcell
             call updurw(npp,nppmax,c,u,y,nindiv,nloc,nloc,
     &            nloc2,nalmax,npopmax,f,indcell,distcell,
     &            indcelltmp,distcelltmp,t,xlim,ylim,du,ploidy,nudcel)

*******************             
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.d-300) then 
c                write(*,*) 'update t'
                call updt(npp,nppmax,nindiv,
     &               nloc,nloc,nloc2,nalmax,npopmax,
     &               t,ttmp,dt,s,c,indcell,distcell,
     &               indcelltmp,distcelltmp,u,y,f,ploidy)
             endif
*******************
*     birth/death des points du pp
c             write(*,*) 'update npp'
             call bdpp(nindiv,u,c,utmp,ctmp,
     &            npop,npopmax,nloc,nloc,
     &            nloc2,nalmax,npp,nppmax,y,f,t,xlim,ylim,indcell,
     &            distcell,indcelltmp,distcelltmp,lambda,ploidy)
          
          endif

*******************
*     birth/death de pop
          if(kfix .eq. 0) then 
             if(dble(iit)/dble(nit) .ge. 0.1) then 
c     write(*,*) 'update npop'
                if(fmodel .eq. 0) then 
                   call smd(npop,npopmin,npopmax,f,fa,drift,
     &                  nloc,nloc2,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &                  a,ptmp,ftmp,drifttmp,y,cellpop,listcell,
     &                  cellpophost,n,ntmp,ploidy)
                else
                   call sm2(npop,npopmin,npopmax,f,fa,drift,
     &                  nloc,nloc2,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &                  a,ptmp,ftmp,drifttmp,y,cellpop,listcell,
     &                  cellpophost,n,ntmp,ploidy,shape1,shape2)
                endif
             endif
          endif

*******************
*     update true unobserved genotypes
*     in the case where null alleles are suspected
         if(filtna .eq. 1) then
c            call udy2(nindiv,nloc,nloc2,nal,nalmax,y,z,
c     &           npopmax,f,fcy,npop)
            call udyNA(nindiv,nloc,nloc2,nal,nalmax,nppmax,y,z,
     &     c,indcell,npopmax,f,fcy,npop,missloc) 
         endif    
*     in the case of dominant data
         if(dom .eq. 1) then
            call udyDOM(nindiv,nloc,nloc2,nal,nalmax,nppmax,y,z,
     &     c,indcell,npopmax,f,fcy,npop)
         endif                  

      enddo
*     end of main loop

*     closing all files
      if(ofiles(1) .eq.1) then
         close(9)
      endif
      if(ofiles(2) .eq.1) then
         close(10)
      endif
      if(ofiles(3) .eq.1) then
         close(11)
      endif   
      if(ofiles(4) .eq.1) then
         close(12)
      endif
      if(ofiles(5) .eq.1) then
         close(13)
      endif
      if(ofiles(6) .eq.1) then
         close(14)
      endif
      if(ofiles(7) .eq.1) then
         close(15)
      endif
      if(ofiles(8) .eq.1) then
         close(16)
      endif
      if(ofiles(9) .eq.1) then
         close(17)
      endif
      if(ofiles(10) .eq.1) then
         close(18)
      endif
      if(ofiles(1) .eq.1) then
         close(19)
      endif
      if(ofiles(12) .eq.1) then
         close(20)
      endif

      call intpr('************************************',-1,0,0)
      call intpr('***    End of MCMC simulation    ***',-1,0,0)
      call intpr('************************************',-1,0,0)

c$$$      write(*,*) 'nindiv =', nindiv,"\n"
c$$$      write(*,*) 'nloc =',nloc   ,"\n"
c$$$      write(*,*) 'nloc2 =', nloc2,"\n"
c$$$      write(*,*) 'nal =',nal ,"\n"
c$$$      write(*,*) 'nalmax =', nalmax,"\n"
c$$$      write(*,*) 'z =', z,"\n"
c$$$      write(*,*) 's =', s,"\n"
c$$$      write(*,*) 'lambda=',lambda,"\n"
c$$$      write(*,*) 'lambdamax=',lambdamax  ,"\n"
c$$$      write(*,*) 'npp =', npp ,"\n"
c$$$      write(*,*) 'nppmax =', nppmax   ,"\n"  
c$$$      write(*,*) 'npopmin =',npopmin  ,"\n"
c$$$      write(*,*) 'npopmax =', npopmax ,"\n" 
c$$$      write(*,*) 'c =', c,"\n"
c$$$      write(*,*) 'ctmp =',ctmp ,"\n"
c$$$      write(*,*) 'u =',u  ,"\n"
c$$$      write(*,*) 'utmp =',utmp ,"\n"
c$$$      write(*,*) 'f =', f  ,"\n"
c$$$      write(*,*) 'ftmp =', ftmp,"\n"
c$$$      write(*,*) 'nit =', nit ,"\n" 
c$$$      write(*,*) 'thinning =',thinning  ,"\n"
c$$$      write(*,*) 'indcell =', indcell ,"\n"
c$$$      write(*,*) 'distcell =', distcell,"\n"
c$$$      write(*,*) 'indcelltmp =', indcelltmp ,"\n"
c$$$      write(*,*) 'distcelltmp =',distcelltmp ,"\n"
c$$$      write(*,*) 'xlim=', xlim 
c$$$      write(*,*) 'ylim=', ylim 


      call rndend()
      end subroutine mcmcgld


c$$$
c$$$*************************************************************************      
c$$$*     update matrix of true genotypes
c$$$*     if given genotypes are true corrupted by the presence
c$$$*     of null alleles
c$$$C CONTIENT UNE ERREUR 
c$$$C TOUS LES INDIVS SONT MIS A JOUR POUR CHAQUE POP
c$$$      subroutine udy2(nindiv,nloc,nloc2,nal,nalmax,y,z,
c$$$     &     npopmax,f,fcy,npop)
c$$$      implicit none
c$$$      integer nindiv,nloc,nloc2,nal,nalmax,y,z,npopmax,npop,nppmax
c$$$      double precision f,fcy
c$$$      dimension nal(nloc),y(nindiv,nloc2),z(nindiv,nloc2),
c$$$     &     f(npopmax,nloc,nalmax),fcy(nalmax,2)
c$$$
c$$$      integer iindiv,iloc,ial1,ipop,ial2,ial,yy,alpha
c$$$      double precision u,ggrunif,sp
c$$$
c$$$c      write(*,*) 'debut udy2'
c$$$
c$$$      do iloc = 1,nloc
c$$$         do ipop = 1,npop
c$$$*     computes posterior proba of true genotypes 
c$$$*     given allele freq and observed genotypes (= true genotypes 
c$$$*     blurred by null alleles)
c$$$            call postpyNA(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
c$$$*     sample y
c$$$            do iindiv = 1,nindiv
c$$$*     only for indiv in pop ipop 
c$$$*     only for indiv with ambigous genotype (homozygous)
c$$$               if(z(iindiv,2*iloc-1) .eq. z(iindiv,2*iloc)) then
c$$$*     case doubly missing data
c$$$                  if((z(iindiv,2*iloc-1) .eq. -999)) then 
c$$$                     y(iindiv,2*iloc-1) = nal(iloc)
c$$$                     y(iindiv,2*iloc)   = nal(iloc)
c$$$                  else
c$$$*     case homozygous (non missing data)
c$$$                     u=ggrunif(0.d0,1.d0)
c$$$                     alpha = z(iindiv,2*iloc-1)
c$$$                     if(u .le. fcy(alpha,1)) then 
c$$$                        y(iindiv,2*iloc-1) = alpha
c$$$                        y(iindiv,2*iloc)   = alpha
c$$$                     else
c$$$                        y(iindiv,2*iloc-1) = alpha
c$$$                        y(iindiv,2*iloc)   = nal(iloc)
c$$$                     endif
c$$$                  endif
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$c      write(*,*) 'fin udy2'
c$$$      end subroutine udy2       



*************************************************************************      
*     update matrix of true genotypes
*     if given genotypes are true corrupted by the presence
*     of null alleles
      subroutine udyNA(nindiv,nloc,nloc2,nal,nalmax,nppmax,y,z,
     &     c,indcell,npopmax,f,fcy,npop,missloc)
      implicit none
      integer nindiv,nloc,nloc2,nal,nalmax,y,z,npopmax,npop,nppmax,
     &     c,indcell,missloc
      double precision f,fcy
      dimension nal(nloc),y(nindiv,nloc2),z(nindiv,nloc2),
     &     f(npopmax,nloc,nalmax),fcy(nalmax,2),c(nppmax),
     &     indcell(nindiv),missloc(nindiv,nloc)
      integer iindiv,iloc,ial1,ipop,ial2,ial,yy,alpha
      double precision u,ggrunif,sp
c      write(*,*) 'debut udyNA'
      do iloc = 1,nloc
         do ipop = 1,npop
*     computes posterior proba of true genotypes 
*     given allele freq and observed genotypes (i.e. true genotypes 
*     blurred by null alleles)
            call postpyNA(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
*     sample y
            do iindiv = 1,nindiv
*     only for indiv in pop ipop 
               if(c(indcell(iindiv)) .eq. ipop) then
*     only for indiv with ambigous genotype (homozygous)
                  if(z(iindiv,2*iloc-1) .eq. z(iindiv,2*iloc)) then
*     case doubly missing data NOT at a missing locus
                     if((z(iindiv,2*iloc-1) .eq. -999) .and.
     &                    missloc(iindiv,iloc) .eq. 0) then 
                        y(iindiv,2*iloc-1) = nal(iloc)
                        y(iindiv,2*iloc)   = nal(iloc)
                     else   
*     case homozygous (non missing data)
                        u=ggrunif(0.d0,1.d0)
                        alpha = z(iindiv,2*iloc-1)
                        if(u .le. fcy(alpha,1)) then 
                           y(iindiv,2*iloc-1) = alpha
                           y(iindiv,2*iloc)   = alpha
                        else
                           y(iindiv,2*iloc-1) = alpha
                           y(iindiv,2*iloc)   = nal(iloc)
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
c      write(*,*) 'fin udyNA'
      end subroutine udyNA


*     computes posterior proba of true genotypes 
*     given allele freq and observed genotypes (obs. = true genotypes 
*     corrupted by null alleles)
*     if presence of null allele is assumed, an extra allele 
*     is assumed and information relative to this allele 
*     is stored in the last non empty entry of f,fcy,...
      subroutine postpyNA(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
      implicit none
      integer ipop,iloc,npopmax,nloc,nalmax,nal
      double precision f,fcy
      dimension f(npopmax,nloc,nalmax),nal(nloc),fcy(nalmax,2)
      integer ial
c$$$      write(*,*) 'debut postpyNA'
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax 
c$$$      write(*,*) 'fcy=',fcy
c$$$      write(*,*) 'f=',f
*     visit all "genuine" alleles
      do ial = 1,nal(iloc)-1
*     proba to have a genuine homoziguous ial,ial
            fcy(ial,1) = f(ipop,iloc,ial)/
     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
*     proba to have a false homoziguous
c            fcy(ial,2) = 2*f(ipop,iloc,nal(iloc))/
c     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
      enddo
c      write(*,*) 'fin postpyNA'
      end subroutine postpyNA



*************************************************************************      
*     update the matrix of true genotypes
*     in the case of dominant markers such as AFLP
      subroutine udyDOM(nindiv,nloc,nloc2,nal,nalmax,nppmax,y,z,
     &     c,indcell,npopmax,f,fcy,npop)
      implicit none
      integer nindiv,nloc,nloc2,nal,nalmax,y,z,npopmax,npop,nppmax,
     &     c,indcell
      double precision f,fcy
      dimension nal(nloc),y(nindiv,nloc2),z(nindiv,nloc2),
     &     f(npopmax,nloc,nalmax),fcy(nalmax,2),c(nppmax),
     &     indcell(nindiv)

      integer iindiv,iloc,ial1,ipop,ial2,ial,yy,alpha
      double precision u,ggrunif,p

c      write(*,*) 'debut udyDOM'

      do iloc = 1,nloc
         do ipop = 1,npop
*     computes posterior proba of genuine homozigous
*     given allele freq and data 
c            call postpyDOM(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
            p = f(ipop,iloc,2)/(f(ipop,iloc,2) + 2*f(ipop,iloc,1))

*     sample y
            do iindiv = 1,nindiv
*     only for indiv in pop ipop 
               if(c(indcell(iindiv)) .eq. ipop) then
*     only for indiv with ambigous obs. (presence of a band)
                  if((z(iindiv,2*iloc-1) .eq. z(iindiv,2*iloc)) .and. 
     &                 (z(iindiv,2*iloc-1) .eq. 2)) then
                     u=ggrunif(0.d0,1.d0)
                     if(u .le. p) then 
                        y(iindiv,2*iloc-1) = 2
                        y(iindiv,2*iloc)   = 2
                     else
                        y(iindiv,2*iloc-1) = 2
                        y(iindiv,2*iloc)   = 1
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
c      write(*,*) 'fin udyDOM'
      end subroutine udyDOM


c$$$*     computes posterior proba of true genotypes 
c$$$*     given allele freq and data 
c$$$*     in the case of dominant markers such as AFLP
c$$$      subroutine postpyDOM(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
c$$$      implicit none
c$$$      integer ipop,iloc,npopmax,nloc,nalmax,nal
c$$$      double precision f,fcy
c$$$      dimension f(npopmax,nloc,nalmax),
c$$$     &     nal(nloc),fcy(nalmax,2)
c$$$      integer ial
c$$$c$$$      write(*,*) 'debut postpyDOM'
c$$$c$$$      write(*,*) 'nal=',nal
c$$$c$$$      write(*,*) 'nalmax=',nalmax 
c$$$c$$$      write(*,*) 'fcy=',fcy
c$$$c$$$      write(*,*) 'f=',f
c$$$*     allele 2 (ial=2) codes the presence of a band 
c$$$*     proba to have a genuine homoziguous 
c$$$      fcy(ial,1) = f(ipop,iloc,2)/
c$$$     &     (f(ipop,iloc,2) + 2*f(ipop,iloc,1))
c$$$*     proba to have a false homoziguous
c$$$c      fcy(ial,2) = 2*f(ipop,iloc,1)/
c$$$c     &     (f(ipop,iloc,2) + 2*f(ipop,iloc,1))
c$$$c      write(*,*) 'fin postpyDOM'
c$$$      end subroutine postpyDOM  




*     Limites du rectangle contenant les coordonnees
      subroutine limit(nindiv,s,xlim,ylim,dt)
      implicit none
      integer nindiv
      double precision s(2,nindiv),xlim(2),ylim(2),dt
      integer iindiv
      xlim(1) = 1.d+300
      xlim(2) = -1.d+300
      ylim(1) = 1.d+300
      ylim(2) = -1.d+300
      do iindiv=1,nindiv
         xlim(1) = dmin1(s(1,iindiv),xlim(1))
         xlim(2) = dmax1(s(1,iindiv),xlim(2))
         ylim(1) = dmin1(s(2,iindiv),ylim(1))
         ylim(2) = dmax1(s(2,iindiv),ylim(2))
      enddo
      xlim(1) = xlim(1) - dt*.5
      xlim(2) = xlim(2) + dt*.5
      ylim(1) = ylim(1) - dt*.5
      ylim(2) = ylim(2) + dt*.5
      end


*     points uniformes dans [0,1]x[0,1]
      subroutine rprioru(npp,nppmax,xlim,ylim,u)
      implicit none
      integer npp,nppmax
      double precision u(2,nppmax),ggrunif,xlim(2),ylim(2)
      integer i
c      call intpr('Begin init u ',-1,0,0)
c      write(*,*) 'npp=',npp
c      write(*,*) 'nppmax=',nppmax
c      write(*,*) 'xlim=',xlim
c      write(*,*) 'ylim=',ylim
c      write(*,*) 'u=',u
      do i=1,npp
         u(1,i) = xlim(1)+(xlim(2)-xlim(1))*ggrunif(0.d0,1.d0)
         u(2,i) = ylim(1)+(ylim(2)-ylim(1))*ggrunif(0.d0,1.d0)
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1, nppmax
            u(1,i) = -999.
            u(2,i) = -999.
         enddo
      endif
c      call intpr('End init u ',-1,0,0)
      end

*     affectation dans les pops selon une loi uniforme
      subroutine rpriorc(npp,nppmax,npop,c)
      implicit none
      integer npp,nppmax,npop,c(nppmax)
      double precision ggrunif
      integer i
      do i=1,npp
         c(i) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1,nppmax
            c(i) = -999
         enddo
      endif
      end

********************************************************************
*     init de la dérive 
*     selon un prior beta(shape1,shape2)
      subroutine rpriordrift(npop,npopmax,drift,fmodel,shape1,shape2)
      implicit none
      integer npop,npopmax,fmodel
      double precision drift(npopmax),shape1,shape2
      integer ipop
      double precision ggrbet
      if(fmodel .eq. 0) then
         do ipop=1,npop
            drift(ipop) = 0.5d0
         enddo
      else
         do ipop=1,npop
            drift(ipop) = ggrbet(shape1,shape2)
         enddo
      endif
      if(npopmax .gt. npop) then
         do ipop=npop+1,npopmax
            drift(ipop) = -999
         enddo
      endif
      end subroutine rpriordrift

*     simulation d'une Dirichlet(1,...,1)
*     (p(1),...,p(k)) uniforme dans l'ensemble {p(1)+...+p(k)=1}
      subroutine dirichlet1(nal,nalmax,p)
      implicit none
      integer nal,nalmax
      double precision p(nalmax)
      integer i
      double precision s,ggrexp
      s = 0.
      do i=1,nal
         p(i) = ggrexp(1.d0)
         s = s + p(i)
      enddo
      do i=1,nal
         p(i) =  p(i)/s
      enddo
      if(nalmax .gt. nal) then
         do i=nal+1,nalmax
            p(i) =  -1
         enddo
      endif
      end

*     simulation d'une Dirichlet(a1,...,an)
      subroutine dirichlet(n,nmax,a,p)
      implicit none
      integer n,nmax
      double precision a(nmax),p(nmax)
      integer i
      double precision s,ggrgam
c      write(*,*) 'debut dirichlet'
      s = 0.
      do i=1,n
         p(i) = 0.
         do while(p(i) .lt. 1d-300) 
c            p(i) = ggrgam(1.d0,a(i))
            p(i) = ggrgam(a(i),1.d0)
         enddo
         s = s + p(i)
      enddo
      do i=1,n
         p(i) =  p(i)/s
      enddo
      if(nmax .gt. n) then
         do i=n+1,nmax
            p(i) =  -1
         enddo
      endif
c      write(*,*) 'fin dirichlet'
      end

      
      subroutine rank(n,nmax,x,p)
      implicit none 
      integer n,nmax,p(nmax)
      double precision x(nmax)
      integer i,j
      p(1) = 1
      do i=2,n
         p(i) = 1
         do j=1,i-1
            if(x(i) .lt. x(j)) then 
               p(i) = p(i) + 1
            else 
               p(j) = p(j) + 1
            endif
         enddo
      enddo      
      end

*     from numerical recipe p 233
      subroutine indexx(n,nmax,arrin,indx)
      dimension arrin(nmax),indx(nmax),arrtmp(nmax)

      do j=1,n
         arrtmp(j) = - arrin(j)
      enddo
      do j=1,n
         indx(j) = j
      enddo
      if(nmax .gt. n) then 
         do j=n+1, nmax
            indx(j) = j
         enddo
      endif
      if(n .eq. 1) return
      l=n/2+1
      ir = n
 10   continue
      if(l .gt. 1) then
         l = l-1
         indxt = indx(l)
         q = arrtmp(indxt)
      else
         indxt = indx(ir)
         q = arrtmp(indxt)
         indx(ir) = indx(1)
         ir = ir -1
         if(ir .eq. 1) then 
            indx(1) = indxt
            return
         endif
      endif
      i = l
      j = l+l
 20   if(j .le. ir) then
         if(j .lt. ir) then
            if(arrtmp(indx(j)) .lt. arrtmp(indx(j+1))) j = j+1
         endif
         if(q .lt. arrtmp(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j+j
         else
            j = ir+1
         endif
         go to 20
      endif
      indx(i) = indxt
      go to 10
      end


*     tirage des frequences dans toutes les pops
*     a tous les locus
      subroutine rpriorf(npop,npopmax,nloc,nlocmax,nal,nalmax,f,
     &     ptmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),nalmax
      double precision f(npopmax,nlocmax,nalmax)
      integer k,l,i
      double precision ptmp(nalmax)
c      call intpr('in rpriorf',-1,0,0)
c      call intpr('npop=',-1,npop,1)
c      call intpr('nloc=',-1,nloc,1)
c      call intpr('nalmax=',-1,nalmax,1)

      do  k=1,npop
c           call intpr('k=',-1,k,1)
         do l=1,nloc
c             call intpr('l=',-1,l,1)
            call dirichlet1(nal(l),nalmax,ptmp)
            do i=1,nalmax
c                call intpr('i=',-1,i,1)
               f(k,l,i)  = ptmp(i)
            enddo
         enddo
      enddo
c      call intpr('in rpriorf non dummy entries done',-1,0,0)

      if(npopmax .gt. npop) then
         do k=npop+1, npopmax
            do l=1,nloc
               do i=1,nalmax
                  f(k,l,i)  = -999
               enddo
            enddo
         enddo
      endif
      end subroutine rpriorf

*     tirage des frequences dans la pop ancestrale 
*     a tous les locus
      subroutine rpriorfa(nloc,nlocmax,nal,nalmax,fa,fmodel,ptmp)
      implicit none
      integer nloc,nlocmax,nal(nlocmax),nalmax,fmodel
      double precision fa(nlocmax,nalmax)
      integer l,i
      double precision ptmp(nalmax)
c      call intpr('in rpriorfa',-1,0,0)
      if(fmodel .eq. 0) then
         do l=1,nloc
            do i=1,nalmax
               fa(l,i)  = 1
            enddo
         enddo
      else
        do l=1,nloc
           call dirichlet1(nal(l),nalmax,ptmp)
           do i=1,nalmax
               fa(l,i)  = ptmp(i)
            enddo
         enddo 
      endif
      end subroutine rpriorfa


*     Mise a jour gibbsienne de f
*     prior p(f) Dirichlet(1,...,1) 
*     p(f|...) Dirichlet(1+ n1,..., 1+np)
*     ni = nbre d'alleles observes
      subroutine rpostf(npop,npopmax,nloc,nlocmax,nal,nalmax,f,
     &     nindiv,nlocmax2,z,nppmax,c,indcell,n,a,f11tmp)
      implicit none 
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),nalmax,
     &     nindiv,nlocmax2,nppmax,c(nppmax),
     &     indcell(nindiv)
      double precision f(npopmax,nlocmax,nalmax)
      integer ipop,iloc,iindiv,ial,n(npopmax,nlocmax,nalmax),
     &     z(nindiv,nlocmax2)
      double precision a(nalmax),f11tmp(nalmax)

*     comptage des effectifs
      do ipop = 1,npop
         do iloc = 1,nloc
            do ial =1,nal(iloc)
               n(ipop,iloc,ial)=0
            enddo
         enddo
      enddo
      do iindiv = 1,nindiv
         do iloc = 1,nloc
            if(z(iindiv,2*iloc-1) .ne. -999) then
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) + 1 
            endif
            if(z(iindiv,2*iloc) .ne. -999) then 
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) + 1 
            endif
         enddo
      enddo
*     tirage Dirichlet
      do ipop = 1,npop
         do iloc = 1,nloc
            do ial =1,nal(iloc)
               a(ial) = 1+dble(n(ipop,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,f11tmp)
            do ial =1,nal(iloc)
               f(ipop,iloc,ial) = f11tmp(ial)
            enddo
         enddo
      enddo
      end



*     Mise a jour gibbsienne de f
*     prior p(f) Dirichlet
*     et une paramétrisation style F-model
*     p(f|...)  est aussi Dirichlet 
*     ni = nbre d'alleles observes
*     (Cf Falush P. 26)
      subroutine rpostf2 (npop,npopmax,nloc,nlocmax,nal,nalmax,
     &     f,fa,drift,
     &     nindiv,nlocmax2,z,nppmax,c,indcell,n,a,ptmp,
     &     ploidy)
      implicit none 
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),nalmax,
     &     nindiv,nlocmax2,nppmax,c(nppmax),
     &     indcell(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax),fa(nlocmax,nalmax),
     &     drift(npopmax)
      integer ipop,iloc,iindiv,ial,n(npopmax,nlocmax,nalmax),
     &     z(nindiv,nlocmax2)
      double precision a(nalmax),ptmp(nalmax)

c$$$      write(*,*) ''
c$$$      write(*,*) ''
c$$$      write(*,*) 'dans rpostf2'
c$$$      write(*,*) 'npop=',npop
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nlocmax=',nlocmax
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax
c$$$      write(*,*) 'f=',f
c$$$      write(*,*) 'fa=',fa
c$$$      write(*,*) 'drift=',drift
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nlocmax2=',nlocmax2
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'npp=',npp
c$$$      write(*,*) 'nppmax=',nppmax
c$$$      write(*,*) 'c=',c
c$$$      write(*,*) 'indcell=',indcell
c$$$      write(*,*) 'n=',n
c$$$      write(*,*) 'a=',a
c$$$      write(*,*) 'ptmp=',ptmp
      
*     comptage des effectifs
      do ipop = 1,npop
         do iloc = 1,nloc
            do ial =1,nal(iloc)
               n(ipop,iloc,ial)=0
            enddo
         enddo
      enddo
      do iindiv = 1,nindiv
         do iloc = 1,nloc
            if(z(iindiv,2*iloc-1) .ne. -999) then
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) + 1 
            endif
            if(z(iindiv,2*iloc) .ne. -999) then 
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) + 1 
            endif
         enddo
      enddo
*     tirage Dirichlet
      do ipop = 1,npop
c         write(*,*) 'ipop=', ipop
         do iloc = 1,nloc
c            write(*,*) 'iloc=',iloc
            do ial =1,nal(iloc)
c               write(*,*) 'ial=',ial
               if(ploidy .eq. 1) then
                  a(ial) = fa(iloc,ial)*(1/drift(ipop)-1) + 
     &              dble(n(ipop,iloc,ial))/2
               endif
               if(ploidy .eq. 2) then
                  a(ial) = fa(iloc,ial)*(1/drift(ipop)-1) + 
     &              dble(n(ipop,iloc,ial))
               endif
c               write(*,*) 'a(',ial,')=',a(ial) 
c               write(*,*) 'fa(',iloc,',',ial,')=', fa(iloc,ial)
c               write(*,*) 'drift(',ipop,')=',drift(ipop)
c               write(*,*) 'n(',ipop,',',iloc,',',ial,')=', 
c     &              n(ipop,iloc,ial)
            enddo
c            write(*,*) 'fa=',fa
c            write(*,*) 'drift=',drift
c            write(*,*) 'a=',a
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial =1,nal(iloc)
               f(ipop,iloc,ial) = ptmp(ial)
            enddo
         enddo
      enddo
      end subroutine rpostf2





*
*     Mise à jour M-H des freq allelique de la pop ancestrale 
*     (F-Model de Falush et al.)
*     fa admet un prior Dirichlet(1,...,1)
      subroutine updfa(npop,npopmax,nlocmax,nalmax,nal,f,fa,drift)
      implicit none 
      integer npop,npopmax,nlocmax,nalmax,nal(nlocmax)
      double precision f(npopmax,nlocmax,nalmax),fa(nlocmax,nalmax),
     &     drift(npopmax)
      integer iloc,ial1,ial2,ipop
      double precision delta,ggrunif,ggrnorm,sigdelta,fa1,fa2,ratio,
     &     lratio,gglgamfn,q,u
      parameter(sigdelta = 0.05) 
      
*     boucle sur les loci
      do iloc=1,nlocmax

*     tirage des deux formes alleliques dont les freq seront 
*     mises à jour 
         ial1 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
         ial2 = ial1 

         do while(ial2 .eq. ial1)
c            write(*,*) 'dans le while'
            ial2 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
         enddo

*     tirage de l'increment
         delta = ggrnorm(0.d0,1.d0)*sigdelta

*     perturbation des deux freq
         fa1 = fa(iloc,ial1) + delta
         fa2 = fa(iloc,ial2) - delta
         if(((fa1 .gt. 1d-300) .and. (1-fa1 .gt. 1d-300)) .and.
     &      ((fa2 .gt. 1d-300) .and. (1-fa2 .gt. 1d-300))) then 
*     calcul du log du ratio 
            lratio = 0.
            do ipop = 1,npop
               q = (1.d0-drift(ipop))/drift(ipop)
               lratio = lratio 
     &              + gglgamfn(fa(iloc,ial1)*q)-gglgamfn(fa1*q)
     &              + gglgamfn(fa(iloc,ial2)*q)-gglgamfn(fa2*q)
     &              + delta*q
     &              *log(f(ipop,iloc,ial1)/f(ipop,iloc,ial2))
            enddo
            lratio = dmin1(0.d0,lratio) 
            ratio = dexp(lratio)

c$$$            write(*,*) 'delta=', delta
c$$$            write(*,*) 'fa1=',fa1
c$$$            write(*,*) 'fa2=',fa2
c$$$            write(*,*) 'q=',q
c$$$            write(*,*) 'fa(iloc,ial2)*q=',fa(iloc,ial2)*q
c$$$            write(*,*) 'gamma(fa(iloc,ial1)*q)',gamma(fa(iloc,ial1)*q)
c$$$            write(*,*) 'ratio=',ratio
c$$$            write(*,*) 'delta*q=',delta*q
c$$$            write(*,*) ''

            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               fa(iloc,ial1) = fa1 
               fa(iloc,ial2) = fa2
            endif
         endif 
      enddo
      end subroutine updfa
      

*
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. beta sur chaque composante
      subroutine upddrift(npop,npopmax,nlocmax,nalmax,nal,
     &     f,fa,drift,shape1,shape2)
      implicit none 
      integer npop,npopmax,nlocmax,nalmax,nal(nlocmax)
      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax)
      integer ipop,iloc,ial
      double precision dtmp,q,qtmp,sigdelta,ratio,lratio,shape1,shape2,
     &     sall,ggrnorm,gglgamfn,u,ggrunif
c      parameter(sigdelta = 0.01)
      sigdelta = 0.5*shape1/(shape1+shape2)

*     boucle sur les pops
      do ipop=1,npop
*     proposition nouvelle valeur
         dtmp = drift(ipop) + ggrnorm(0.d0,1.d0)*sigdelta
         q = (1-drift(ipop))/drift(ipop)
         qtmp = (1-dtmp)/dtmp
         if((dtmp .gt. 1d-300 ) .and. (1-dtmp .gt. 1d-300)) then 

*     calcul du log du ratio
c     prior uniforme
            lratio = 0 
c     prior beta(shape1,shape2)
            lratio = (shape1-1)*dlog(dtmp/drift(ipop)) + 
     &           (shape2-1)*dlog((1-dtmp)/(1-drift(ipop)))
            do iloc=1,nlocmax
               sall = 0.
               do ial = 1,nal(iloc)
                  sall = sall + gglgamfn(fa(iloc,ial)*q)-
     &                 gglgamfn(fa(iloc,ial)*qtmp) +
     &                 fa(iloc,ial)*(qtmp-q)*dlog(f(ipop,iloc,ial))
               enddo
               lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
            enddo
 
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               drift(ipop) = dtmp 
            endif
         endif
      enddo
      end subroutine upddrift

c$$$*
c$$$*     Mise à jour M-H du vecteur de dérives génétiques 
c$$$*     prior indep. uniforme sur chaque composante
c$$$      subroutine upddrift2(npop,npopmax,nlocmax,nalmax,nal,
c$$$     &     f,fa,drift)  
c$$$      implicit none 
c$$$      integer npop,npopmax,nlocmax,nalmax,nal(nlocmax)
c$$$      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
c$$$     &     fa(nlocmax,nalmax)
c$$$      integer ipop,iloc,ial
c$$$      double precision d,q,qtmp,sigdelta,ratio,lratio,alpha,sall,ggrnorm,
c$$$     &     gglgamfn,u,ggrunif
c$$$      parameter(sigdelta = 0.01,alpha=5000) 
c$$$
c$$$*     boucle sur les pops
c$$$      do ipop=1,npop
c$$$*     proposition nouvelle valeur
c$$$         d = drift(ipop) + ggrnorm()*sigdelta
c$$$         q = (1-drift(ipop))/drift(ipop)
c$$$         qtmp = (1-d)/d
c$$$         if((d .gt. 1d-300) .and. (1-d .gt. 1d-300)) then 
c$$$
c$$$*     calcul du log du ratio
c$$$            lratio = 0 
c$$$c     decommenter la ligne suivante pour avoir un prior exponentiel tronqué
c$$$c     sinon le prior est uniforme
c$$$c            lratio = -alpha*(d-drift(ipop))
c$$$            do iloc=1,nlocmax
c$$$               sall = 0.
c$$$               do ial = 1,nal(iloc)
c$$$                  sall = sall + gglgamfn(fa(iloc,ial)*q)-
c$$$     &                 gglgamfn(fa(iloc,ial)*qtmp) +
c$$$     &                 fa(iloc,ial)*(qtmp-q)*dlog(f(ipop,iloc,ial))
c$$$               enddo
c$$$               lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
c$$$            enddo
c$$$
c$$$
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            ratio = exp(lratio)
c$$$            u = ggrunif(0.d0,1.d0)
c$$$            if(u .le. ratio) then 
c$$$               drift(ipop) = d 
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$      end subroutine upddrift2



*     recherche la cellule de chaque individu
*     stockage des indices dans indcell
*     stockage des carres des distances dans distcell
      subroutine calccell(nindiv,s,npp,nppmax,u,
     &     indcell,distcell)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv)
      double precision distcell(nindiv)
      double precision s(2,nindiv),u(2,nppmax)
      integer iindiv,ipp
      double precision d
c      write(*,*) 's=',s
c      write(*,*) 'u=',u
      do iindiv=1,nindiv
         indcell(iindiv) = -999
         distcell(iindiv) = 1.d+300
         do ipp=1,npp
            d = (s(1,iindiv)-u(1,ipp))**2 + (s(2,iindiv)-u(2,ipp))**2
c           write(*,*) 'd=',d
            if( d .lt. distcell(iindiv) ) then 
               indcell(iindiv) = ipp
               distcell(iindiv) = d
c               write(*,*) 'ipp=',ipp
            endif
         enddo
      enddo
c      write(*,*) 'indcell =',indcell
      end




*     mise a jour de indcell et distcell
*     apres le deplacement d'un point de u (celui d'indice j)
      subroutine vormove(nindiv,s,npp,nppmax,u,
     &     indcell,distcell,indcelltmp,distcelltmp,j)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),j
      double precision s(2,nindiv),u(2,nppmax),distcell(nindiv),
     &     distcelltmp(nindiv),d
      integer iindiv,ipp

C       write(6,*) 'debut de  vormove'
C       write(6,*) 'j =', j
C       write(6,*) 'indcell',indcell
C       write(6,*)'distcell',distcell
C       write(6,*) 'indcelltmp',indcelltmp
C       write(6,*)'distcelltmp',distcelltmp

      do iindiv=1,nindiv
         if(indcell(iindiv) .eq. j) then 
*     pour les indiv qui etaient dans la cellule j on cherche
*     la nouvelle cellule
            d = 3.e+37
            indcelltmp(iindiv) = -999
            distcelltmp(iindiv) = 3.e+37
            do ipp=1,npp
               d= (s(1,iindiv)-u(1,ipp))**2+(s(2,iindiv)-u(2,ipp))**2
               if( d .lt. distcelltmp(iindiv) ) then 
                  indcelltmp(iindiv) = ipp
                  distcelltmp(iindiv) = d
               endif
            enddo
*     pour les autres indiv on regarde si le nouveau uj s'est intercale
         else
            d = (s(1,iindiv)-u(1,j))**2+(s(2,iindiv)-u(2,j))**2
            if(d .lt. distcell(iindiv)) then
               indcelltmp(iindiv) = j
               distcelltmp(iindiv) = d
            else
               indcelltmp(iindiv) = indcell(iindiv)
               distcelltmp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo

c$$$         call calccell(nindiv,s,npp,nppmax,u,
c$$$     &        indcelltmp2,distcelltmp2)
c$$$         do iindiv=1,nindiv
c$$$            if((indcelltmp2(iindiv) .ne. indcelltmp(iindiv)) .or. 
c$$$     &         (distcelltmp2(iindiv) .ne. distcelltmp(iindiv)))  then
c$$$               write(6,*) 'fin de  vormove'
c$$$               write(6,*) 'j =', j
c$$$               write(6,*) 'iindiv=',iindiv
c$$$               write(6,*) 'indcell',indcell
c$$$               write(6,*)'distcell',distcell
c$$$               write(6,*) 'indcelltmp',indcelltmp
c$$$               write(6,*)'distcelltmp',distcelltmp
c$$$               write(6,*)'indceltmp2',indcelltmp2
c$$$               write(6,*)'distceltmp2',distcelltmp2
c$$$               stop
c$$$            endif
c$$$         enddo
C          write(6,*) 'fin de  vormove'
C          write(6,*) 'j =', j
C          write(6,*) 'indcell',indcell
C          write(6,*)'distcell',distcell
C          write(6,*) 'indcelltmp',indcelltmp
C          write(6,*)'distcelltmp',distcelltmp
      end



*     mise a jour de indcell et distcell
*     apres naissance d'un point de u 
      subroutine voradd(s,utmp,
     &     indcell,distcell,indcelltmp,distcelltmp,
     &     nindiv,npp,nppmax)
      implicit none 
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),iindiv
      double precision s(2,nindiv),distcell(nindiv),
     &     distcelltmp(nindiv),d,utmp(2,nppmax)
      
      do iindiv =1,nindiv
*     est-ce que le nouveau point s'est intercale ?
         d = (s(1,iindiv)-utmp(1,npp+1))**2+
     &        (s(2,iindiv)-utmp(2,npp+1))**2
         if(d .lt. distcell(iindiv)) then 
            distcelltmp(iindiv) = d
            indcelltmp(iindiv) = npp+1
         else
            distcelltmp(iindiv) = distcell(iindiv)
            indcelltmp(iindiv) = indcell(iindiv) 
         endif
      enddo
      end 

 


     
*     mise a jour de indcell et distcell
*     apres mort d'un point de u 
      subroutine vorrem(s,utmp,ipprem,
     &     indcell,distcell,indcelltmp,distcelltmp,
     &     nindiv,npp,nppmax)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),ipprem,iindiv
      double precision s(2,nindiv),utmp(2,nppmax),
     &     distcell(nindiv),distcelltmp(nindiv),d
      integer ipp
      
      do iindiv =1,nindiv
*     est-ce que le site courant dependait de la cellule disparue ?
         if(indcell(iindiv) .eq. ipprem) then
*     si oui on recherche sa nouvelle cellule parmi celles qui restent
*     (les nouvelles)
            distcelltmp(iindiv) = 3.e+37
            do ipp=1,npp-1
               d = (s(1,iindiv)-utmp(1,ipp))**2+
     &              (s(2,iindiv)-utmp(2,ipp))**2
               if( d .lt. distcelltmp(iindiv) ) then 
                  indcelltmp(iindiv) = ipp
                  distcelltmp(iindiv) = d
               endif
            enddo
         else
            if(indcell(iindiv) .lt. ipprem) then
               indcelltmp(iindiv) = indcell(iindiv)
               distcelltmp(iindiv) = distcell(iindiv)
            else
               indcelltmp(iindiv) = indcell(iindiv) - 1
               distcelltmp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo
      end


*     Tirage de lambda selon p(lambda|m) 
*     avec 
*     p(m|lambda) Poisson translatee
*     p(lambda) uniforme dans [0,lambdamax]
*     p(lambda|m) gamma tronquee
      double precision function rpostlamb(lambdamax,m)
      implicit none 
      double precision lambdamax,ggrgam
      integer m
*      write(*,*) 'beg rpostlamb'
      rpostlamb = lambdamax + 1
      do while(rpostlamb .gt. lambdamax)
c         rpostlamb = ggrgam(1.d0,dble(m))
         rpostlamb = ggrgam(dble(m),1.d0)
      enddo
*      write(*,*) 'end rpostlamb'
      end


*     
*     Mise a jour de c sans modif de npp
      subroutine updc(npp,nppmax,c,ctmp,z,nindiv,nloc,
     &     nlocmax,nlocmax2,nalmax,npop,npopmax,f,indcell,ploidy,nudcel)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,npop,nalmax,npopmax,z(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax)
      integer ipp,ctmp(nppmax),iud
      double precision ggrunif,r,alpha,ratio,ggrbinom,bern
      
      do ipp=1,npp
         ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif
c     write(*,*) ''
      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         r = ratio(z,f,c,ctmp,indcell,indcell,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'c=',c
c         write(*,*) 'ctmp=',ctmp
c         write(*,*) 'r=',r
         alpha = dmin1(1.d0,r)
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
      end subroutine updc





***********************************************************************
*     joint update of c and f under the Dirichlet model
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf(npop,npopmax,f,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,z,n,ntmp,ploidy,alpha,nudcel)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),z(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax),
     &      ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax)
      integer ipop,ipp,ipop1,ipop2,iloc,ial
      double precision ggrunif,lrpf,lratio,ratio,llr6
      integer iipp
      integer n1,n2,ntmp1,ntmp2,iud
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern
      double precision alpha,lrp
c      write(*,*) 'begin udcf'
c      write(*,*) 'npop =',npop
c      write(*,*) 'npp=',npp

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp) 
         ipop2 = ctmp(ipp)
*     counting alleles for states c and ctmp
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,z,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
c         write(*,*) 'in udcf c=',c(1),c(2)
c         write(*,*) '     ctmp=',ctmp(1),ctmp(2)

         bern = 0 
         if(ipop1 .ne. ipop2) then 
            lratio = 0
            do iloc=1,nloc
               n1 = 0
               n2 = 0
               ntmp1 = 0
               ntmp2 = 0
               do ial=1,nal(iloc)
                  lratio = lratio 
     &                 - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                 - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                  n1 = n1 + n(ipop1,iloc,ial)
                  n2 = n2 + n(ipop2,iloc,ial)
                  ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                  ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
               enddo
               lratio = lratio + 
     &              gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &              + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
            enddo

c            write(*,*) 'in udcf lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            bern = ggrbinom(1.d0,ratio)
         endif
         
         if((bern .eq. 1) .or. (ipop1 .eq. ipop2)) then
*     sample new frequencies
c            write(*,*) 'accept move in udcf'
            call samplef(npop,npopmax,nloc,nlocmax,
     &           nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
            c(ipp) = ctmp(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                  f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
               enddo
            enddo
         else 
            ctmp(ipp) = c(ipp)
          endif
      enddo
c      write(*,*) 'end udcf'
      end subroutine udcf
***********************************************************************


c$$$
c$$$
c$$$
c$$$***********************************************************************
c$$$*     joint update of c and f under the Dirichlet model
c$$$*     update of c by splitting a pop
c$$$*     new f is proposed according to full conditionnal pi(f*|c*,z)
c$$$      subroutine udcfsplit(npop,npopmax,f,nloc,nloc2,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
c$$$     &     a,ptmp,ftmp,z,n,ntmp,ploidy,alpha,cellpop,listcell)
c$$$      implicit none
c$$$      integer npop,npopmax,nloc,nal(nloc),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
c$$$     &     nloc2,c(nppmax),ctmp(nppmax),z(nindiv,nloc2),
c$$$     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
c$$$     &     ploidy,cellpop(nppmax),listcell(nppmax)
c$$$      double precision f(npopmax,nloc,nalmax),ftmp(npopmax,nloc,nalmax),
c$$$     &     a(nalmax),ptmp(nalmax)
c$$$      integer ipop,ipp,ipop1,ipop2,iloc,ial
c$$$      double precision ggrunif,lrpf,lratio,ratio,llr6
c$$$      double precision bern,ggrbinom
c$$$      integer iipp
c$$$      integer nu1,nu2,nu,n1,n2,ntmp1,ntmp2
c$$$      double precision junk,termf9bis,gglgamfn
c$$$      double precision alpha,lrp
c$$$c      write(*,*) 'debut udcfsplit'
c$$$* if npop > 1
c$$$
c$$$*     choix de la pop qui split and fake init for ipop2
c$$$      ipop1 = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$      ipop2 = ipop1
c$$$c      write(*,*) 'ipop1=',ipop1
c$$$*     recherche des cellules affectees a cette pop
c$$$      call who(c,ipop1,npp,nppmax,cellpop,nu1)
c$$$c      write(*,*) 'nu1=',nu1
c$$$      bern = 1
c$$$      if(nu1.gt. 0) then
c$$$*     tirage du nombre de cellules reallouees
c$$$         nu = idint(dint(dble(nu1+1)*ggrunif(0.d0,1.d0)))
c$$$c         write(*,*) 'nu=',nu
c$$$         if(nu .gt. 0) then
c$$$*     tirage des cellules reallouees
c$$$            call sample2(cellpop,nppmax,nu,nu1,listcell)
c$$$*     choix de la pop hote
c$$$            do while(ipop2 .eq. ipop1)
c$$$               ipop2 = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            enddo
c$$$c            write(*,*) 'ipop2=',ipop2
c$$$*     proposition de reallocation dans la pop ipop2
c$$$            call split(ipop2,c,ctmp,nppmax,nu,listcell)
c$$$         else
c$$$            do ipp = 1,nppmax
c$$$               ctmp(ipp) = c(ipp)
c$$$            enddo
c$$$         endif
c$$$      else
c$$$         do ipp = 1,nppmax
c$$$            ctmp(ipp) = c(ipp)
c$$$         enddo
c$$$c     write(*,*) 'ipop2=',ipop2
c$$$      endif
c$$$*     counting alleles for states c and ctmp
c$$$      call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &     nppmax,nal,nalmax,z,n,indcell,c,ploidy)
c$$$      call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &     nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
c$$$c     write(*,*) 'alleles counted'
c$$$      lratio = 0
c$$$      do iloc=1,nloc
c$$$         n1 = 0
c$$$         n2 = 0
c$$$         ntmp1 = 0
c$$$         ntmp2 = 0
c$$$         do ial=1,nal(iloc)
c$$$            lratio = lratio 
c$$$     &           - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
c$$$     &           - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
c$$$     &           + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
c$$$     &           + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
c$$$            n1 = n1 + n(ipop1,iloc,ial)
c$$$            n2 = n2 + n(ipop2,iloc,ial)
c$$$            ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
c$$$            ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
c$$$         enddo
c$$$         lratio = lratio + 
c$$$     &        gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
c$$$     &        + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
c$$$     &        - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
c$$$     &        - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
c$$$      enddo
c$$$c      write(*,*) 'lratio=',lratio
c$$$c      write(*,*) 'lRqc=', gglgamfn(dble(nu2+1))
c$$$c     &     - gglgamfn(dble(nu1+nu2+2)) - gglgamfn(dble(nu1-nu+1))
c$$$      lratio = lratio + gglgamfn(dble(nu1+2)) + gglgamfn(dble(nu2+1))
c$$$     &     - gglgamfn(dble(nu1+nu2+2)) - gglgamfn(dble(nu1-nu+1))
c$$$      lratio = dmin1(0.d0,lratio)
c$$$      ratio = dexp(lratio)
c$$$      bern = ggrbinom(1.d0,ratio)
c$$$c      if(bern .eq. 1) write(*,*) 'accept in udcf split'
c$$$      if(bern .eq. 1) then 
c$$$         call samplef(npop,npopmax,nloc,nloc,
c$$$     &        nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
c$$$         do iloc=1,nloc
c$$$            do ial=1,nal(iloc)
c$$$               f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
c$$$               f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
c$$$            enddo
c$$$         enddo  
c$$$         do ipp = 1,npp
c$$$            c(ipp) = ctmp(ipp)
c$$$         enddo
c$$$      endif
c$$$c      write(*,*) 'end udcfsplit'
c$$$      end subroutine udcfsplit
c$$$***********************************************************************
c$$$






***********************************************************************
*     joint update of c and f under the CFM
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf2(npop,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,z,n,ntmp,ploidy,nudcel)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),z(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iud
      double precision alpha,ggrunif,lratio,lTf
      double precision bern,ggrbinom

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp)
         ipop2 = ctmp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,z,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
         
         bern = 0
         if(ipop1 .ne. ipop2) then 
            lratio = lTf(ipop1,ntmp,fa,drift,npopmax,nloc,nal,nalmax) + 
     &           lTf(ipop2,ntmp,fa,drift,npopmax,nloc,nal,nalmax) - 
     &           lTf(ipop1,n,fa,drift,npopmax,nloc,nal,nalmax) - 
     &           lTf(ipop2,n,fa,drift,npopmax,nloc,nal,nalmax)
            
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
         endif
c         write(*,*) 'bern=',bern

          if((bern .eq. 1) .or. (ipop1 .eq. ipop2)) then
            call samplef2(npop,npopmax,nloc,nlocmax,
     &           nal,nalmax,ipop1,ipop2,f,ftmp,
     &           fa,drift,a,ptmp,ntmp)
            c(ipp) = ctmp(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                  f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
               enddo
            enddo
         else 
            ctmp(ipp) = c(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
                  ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
               enddo
            enddo
         endif        
      enddo
      end subroutine udcf2
***********************************************************************



*     Modification de u
*     composante par composante 
*     avec proposal uniforme sur un carre de cote du 
*     centre sur le point courant (random walk)
      subroutine updurw(npp,nppmax,c,u,z,nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,f,indcell,distcell,
     &     indcelltmp,distcelltmp,
     &     s,xlim,ylim,du,ploidy,nudcel)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,z(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel
      double precision u(2,nppmax),f(npopmax,nlocmax,nalmax),
     &     distcell(nindiv),s(2,nindiv),xlim(2),ylim(2),du
      integer ipp,iindiv,indcelltmp(nindiv),iud
      double precision utmp(2,nppmax),ggrunif,r,alpha,
     &     distcelltmp(nindiv),surf,surftmp,dx,dy,ratio
      double precision bern,ggrbinom
*     initialisation du tableau tmporaire
      do ipp=1,npp
         utmp(1,ipp) = u(1,ipp)
         utmp(2,ipp) = u(2,ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            utmp(1,ipp) = -999.
            utmp(2,ipp) = -999.
         enddo
      endif

c      write(*,*) 'npp=', npp
c      write(*,*) 'u=', u

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     proposition d un deplacement d un point de u
         utmp(1,ipp) = max(u(1,ipp)-du/2.,xlim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(1,ipp)+du/2.,xlim(2))-max(u(1,ipp)-du/2.,xlim(1)))
         utmp(2,ipp) = max(u(2,ipp)-du/2.,ylim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(2,ipp)+du/2.,ylim(2))-max(u(2,ipp)-du/2.,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
         surf = (dx+du/2.)*(dy+du/2.)
         dx = min(du/2.,utmp(1,ipp)-xlim(1),xlim(2)-utmp(1,ipp))
         dy = min(du/2.,utmp(2,ipp)-ylim(1),ylim(2)-utmp(2,ipp))
         surftmp = (dx+du/2.)*(dy+du/2.)

*     modif de indcell et distcell
         call vormove(nindiv,s,npp,nppmax,utmp,
     &        indcell,distcell,indcelltmp,distcelltmp,ipp)

c         write(*,*) 'apres vormove'



         r = ratio(z,f,c,c,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'r=',r
         r = r*surf/surftmp

         alpha = dmin1(1.d0,r)
c         write(*,*) 'alpha=',alpha
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            u(1,ipp) = utmp(1,ipp)
            u(2,ipp) = utmp(2,ipp)
            do iindiv=1,nindiv
               indcell(iindiv) = indcelltmp(iindiv)
               distcell(iindiv) = distcelltmp(iindiv)
            enddo
         else 
            utmp(1,ipp) = u(1,ipp)
            utmp(2,ipp) = u(2,ipp)
         endif
      enddo
      end subroutine updurw




*
*     mise a jour de t    
* 
      subroutine updt(npp,nppmax,nindiv,
     &     nloc,nlocmax,nlocmax2,nalmax,npopmax,
     &     t,ttmp,dt,s,c,indcell,distcell,indcelltmp,distcelltmp,
     &     u,z,f,ploidy)
      implicit none 
      integer npp,nppmax,nindiv,nloc,nlocmax,nlocmax2,nalmax,
     &     npopmax,c(nppmax),indcell(nindiv),z(nindiv,nlocmax2),
     &     ploidy
      double precision t(2,nindiv),s(2,nindiv),distcell(nindiv),
     &     u(2,nppmax),f(npopmax,nlocmax,nalmax),dt
      integer iindiv,ipp,accept,indcelltmp(nindiv)
      double precision ggrunif,d,ttmp(2,nindiv),r,alpha,
     &     distcelltmp(nindiv),ratio
      double precision bern,ggrbinom

*     initialisation
      do iindiv = 1,nindiv
         ttmp(1,iindiv) = t(1,iindiv)
         ttmp(2,iindiv) = t(2,iindiv)
         indcelltmp(iindiv) = indcell(iindiv)
         distcelltmp(iindiv) = distcell(iindiv)
      enddo

c      do iindiv = 1,nindiv
      iindiv= 1 + idint(dint(dble(nindiv)*ggrunif(0.d0,1.d0)))
*     proposition d'une modif de t
         ttmp(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         ttmp(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)

*     modif de indcell et distcell
         distcelltmp(iindiv) = 3.e+37
         do ipp = 1,npp
            d = (ttmp(1,iindiv)-u(1,ipp))**2+
     &           (ttmp(2,iindiv)-u(2,ipp))**2
            if(d .lt. distcelltmp(iindiv)) then 
               indcelltmp(iindiv)  = ipp
               distcelltmp(iindiv) = d
            endif
         enddo

*     proba d'acceptation
         if(indcelltmp(iindiv) .ne. indcell(iindiv)) then 
            r = ratio(z,f,c,c,indcell,indcelltmp,
     &           npopmax,nlocmax,nalmax,nindiv,nloc,
     &           nlocmax2,nppmax,ploidy)
         else 
            r = 1.
         endif
         alpha = dmin1(1.d0,r)
         accept = ggrbinom(1.d0,alpha)
*     mise a jour en cas d'acceptation
         if(accept .eq. 1) then 
            indcell(iindiv) = indcelltmp(iindiv)
            distcell(iindiv) = distcelltmp(iindiv)
            t(1,iindiv) = ttmp(1,iindiv) 
            t(2,iindiv) = ttmp(2,iindiv)
         endif
c      enddo
      end subroutine updt






*
*     naissance ou mort d'une cellule
*     avec prior Poisson(lambda) tronquée :   1 < m < nppmax
      subroutine bdpp(nindiv,u,c,utmp,ctmp,npop,npopmax,
     &     nloc,nlocmax,nlocmax2,nalmax,npp,nppmax,z,f,s,xlim,ylim,
     &     indcell,distcell,indcelltmp,distcelltmp,lambda,ploidy)
      implicit none 
      integer nindiv,nloc,nlocmax,nlocmax2,
     &     npop,npopmax,
     &     nalmax,npp,nppmax,z(nindiv,nlocmax2),c(nppmax),
     &     indcell(nindiv),ploidy
      double precision u(2,nindiv),f(npopmax,nlocmax,nalmax),xlim(2),
     &     ylim(2),s(2,nindiv),distcell(nindiv),lambda

      integer ctmp(nppmax),indcelltmp(nindiv),ipp,npptmp,
     &     accept,iindiv,ipprem
      double precision utmp(2,nppmax),distcelltmp(nindiv),ggrunif,
     &     ratio,r,alpha,ggrbinom,b
      
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utmp (1,ipp) = u(1,ipp)
               utmp (2,ipp) = u(2,ipp)
               ctmp(ipp) = c(ipp)
            enddo
            npptmp = npp + 1
            ctmp(npptmp) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            utmp(1,npptmp) = xlim(1)+(xlim(2)-xlim(1))*
     &           ggrunif(0.d0,1.d0)
            utmp(2,npptmp) = ylim(1)+(ylim(2)-ylim(1))*
     &           ggrunif(0.d0,1.d0)
            if(nppmax .gt. npptmp) then
               do ipp=npptmp+1,nppmax
                  ctmp(ipp) = -999
                  utmp(1,ipp) = -999.
                  utmp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,utmp,indcell,distcell,indcelltmp,
     &           distcelltmp,nindiv,npp,nppmax)
            r = ratio(z,f,c,ctmp,indcell,indcelltmp,npopmax,nlocmax,
     &           nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy)
            r = r*lambda/dble(npp+1)
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npptmp
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      else
*     mort
         if(npp .ne. 1) then 
            ipprem = 1+ aint(dble(npp)*ggrunif(0.d0,1.d0))
            if(ipprem .ne. 1) then 
               do ipp = 1,ipprem-1
                  utmp (1,ipp) = u(1,ipp)
                  utmp (2,ipp) = u(2,ipp)
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            if(ipprem .ne. npp) then 
               do ipp = ipprem,npp-1
                  utmp (1,ipp) = u(1,ipp+1)
                  utmp (2,ipp) = u(2,ipp+1)
                  ctmp(ipp) = c(ipp+1)
               enddo
            endif
            do ipp=npp,nppmax
               utmp (1,ipp) = -999.
               utmp (2,ipp) = -999.
               ctmp(ipp) = -999
            enddo

            call vorrem(s,utmp,ipprem,indcell,distcell,
     &           indcelltmp,distcelltmp,nindiv,npp,nppmax)

            r = ratio(z,f,c,ctmp,indcell,indcelltmp,npopmax,nlocmax,
     &           nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy)
            r = r*dble(npp)/lambda
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npp-1
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      endif
      end subroutine bdpp






*     calcul du ratio p(z|theta*)/p(z|theta)
*     ca ne depend pas de lambda
      double precision function ratio(z,f,c,ctmp,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,z(nindiv,nlocmax2),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp

c$$$      write(*,*) 'debut de ratio'
c$$$      write(*,*) 'indcell=',indcell
c$$$      write(*,*) 'indcelltmp=',indcelltmp
c$$$      write(*,*) 'c=',c
c$$$      write(*,*) 'ctmp=',ctmp


      ratio = 1.
      do iindiv=1,nindiv
c         write(*,*) 'iindiv=', iindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
C         write(*,*) 'indcell=',indcell
C          write(*,*) 'indcelltmp=',indcelltmp
C          write(*,*) 'c=',c
C          write(*,*) 'ctmp=',ctmp
C          write(*,*) 'ipop=',ipop
C          write(*,*) 'ipoptmp=',ipoptmp

         do iloc=1,nloc
c             write(*,*) 'iloc=',iloc
c            write(6,*) 'z=',z(iindiv,2*iloc-1)
c            write(6,*) 'z=',z(iindiv,2*iloc)
            ial1 = z(iindiv,2*iloc-1)
            ial2 = z(iindiv,2*iloc)
c            ratio = ratio*
c     &           (f(ipoptmp,iloc,ial1)/f(ipop,iloc,ial1))*
c     &           (f(ipoptmp,iloc,ial2)/f(ipop,iloc,ial2))
            if(ial1 .ne. -999) then 
c               write(*,*) f(ipoptmp,iloc,ial1)
c               write(*,*) f(ipop,iloc,ial1)
               ratio = ratio*
     &              (f(ipoptmp,iloc,ial1)/f(ipop,iloc,ial1))
               
            endif
            if(ial2 .ne. -999) then 
c               write(*,*) f(ipoptmp,iloc,ial2)
c               write(*,*) f(ipop,iloc,ial2)
               ratio = ratio*
     &              (f(ipoptmp,iloc,ial2)/f(ipop,iloc,ial2))
            endif
c            write(*,*) 'ratio =',ratio
         enddo
      enddo
      if(ploidy .eq. 1) then 
         ratio = dsqrt(ratio)
      endif
c      write(*,*) 'fin de ratio'
      end function ratio


************************************************************************
*     Indice des cellules dans une pop
      subroutine who(c,ipop,npp,nppmax,cellpop,
     &     ncellpop)
      implicit none
      integer npp,nppmax,c(nppmax),ipop,cellpop(nppmax),ncellpop
      integer ipp,ii
c      write(*,*) 'who'
      ii = 1
      ncellpop = 0
      do ipp=1,npp
         if(c(ipp) .eq. ipop) then
           cellpop(ii) = ipp
           ncellpop = ncellpop + 1
           ii = ii + 1
        endif
      enddo
      if(nppmax .gt. ncellpop) then
         do ipp= ncellpop+1, nppmax
            cellpop(ipp) = -999
         enddo
      endif
      end subroutine who


************************************************************************
*     Tirage de nu cellules parmi ncellpop cellules
*      
      subroutine sample(cellpop,nppmax,nu,ncellpop,listcell)
      implicit none
      integer nppmax,cellpop(nppmax),nu,ncellpop,listcell(nppmax)
      integer isamp,ii,jj
      double precision ggrunif
c      write(*,*) 'sample'

*     init
      ii = 1 + idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
      listcell(1) = cellpop(ii)
      if(nu .gt. 1) then
         do isamp = 2,nu
c             write(*,*) 'cellpop=',cellpop
*     translation
            if(ii .eq. 1) then 
               do jj=2,ncellpop
                  cellpop(jj-1) = cellpop(jj)
               enddo
            else
               if(ii .ne. ncellpop) then
                  do jj=ii+1,ncellpop
                     cellpop(jj-1) = cellpop(jj)
                  enddo
               endif
            endif
            cellpop(ncellpop-isamp+1) = -999
*     tirage parmi les ncellpop-isamp cellules restantes
            ii = 1 + 
     &           idint(dint(dble(ncellpop-isamp)*ggrunif(0.d0,1.d0)))
            listcell(isamp) = cellpop(ii)
         enddo
      endif
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
c      write(*,*) 'listcell=',listcell
      end subroutine sample

***********************************************************************
*
*     Tirage de nu cellules parmi ncellpop cellules
*      version corrigée de sample apres un bug 
*     trouvé en septembre 2005 à Göteborg
      subroutine sample2(cellpop,nppmax,nu,ncellpop,listcell)
      implicit none
      integer nppmax,cellpop(nppmax),nu,ncellpop,listcell(nppmax)
      integer isamp,ii,jj
      double precision ggrunif
c      call rndstart()
c      write(*,*) 'sample2'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
*     init
      ii = 1 + idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
      listcell(1) = cellpop(ii)
c      write(*,*) 'listcell(1)=',listcell(1)
      if(nu .gt. 1) then
         do isamp = 2,nu
c             write(*,*) 'cellpop=',cellpop
*     translation
            if(ii .eq. 1) then 
               do jj=2,ncellpop
                  cellpop(jj-1) = cellpop(jj)
               enddo
            else
               if(ii .ne. ncellpop) then
                  do jj=ii+1,ncellpop
                     cellpop(jj-1) = cellpop(jj)
                  enddo
               endif
            endif
            cellpop(ncellpop-(isamp-1)+1) = -999
c             write(*,*) 'cellpop=',cellpop
*     tirage parmi les ncellpop-isamp cellules restantes
            ii = 1 + 
     &           idint(dint(dble(ncellpop-isamp)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'ii=',ii
            listcell(isamp) = cellpop(ii)
c            write(*,*) 'listcell(isamp)=',listcell(isamp)
         enddo
      endif
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
c      write(*,*) 'listcell=',listcell
c      call rndend()
      end subroutine sample2


*******************************************************************
*     split d'une pop en deux
*     reallocation de nu cellules dont les indices
*     sont dans listcell
*     dans la pop ipop
      subroutine split(ipop,c,ctmp,nppmax,nu,listcell)
      implicit none
      integer ipop,nppmax,c(nppmax),ctmp(nppmax),nu,
     &     listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de split'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipop=',ipop
      do ipp=1,nppmax
         ctmp(ipp) = c(ipp)
      enddo
      if(nu .gt. 0) then
         do ii=1,nu
            ctmp(listcell(ii)) = ipop
         enddo
      endif
c      write(*,*)'c=',c
c      write(*,*)'ctmp=',ctmp
c      write(*,*) 'fin de split'
      end subroutine split

*******************************************************************
*     merge de deux  pops en une : 
*     reallocation des nu cellules de la pop ipoprem 
*     dont les indices sont dans listcell
*     dans la pop ipophost
      subroutine merging(ipoprem,ipophost,
     &     c,ctmp,nppmax,nu,listcell)
      implicit none
      integer ipoprem,ipophost,nppmax,
     &     c(nppmax),ctmp(nppmax),nu,listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de merge'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipoprem=',ipoprem
      do ipp=1,nppmax
         ctmp(ipp) = c(ipp)
      enddo
      if(ipoprem .gt. ipophost) then 
         if(nu .gt. 0) then
            do ii=1,nu
               ctmp(listcell(ii)) = ipophost
            enddo
         endif
      else
         if(nu .gt. 0) then
            do ii=1,nu
               ctmp(listcell(ii)) = ipophost - 1
            enddo
         endif
      endif
      do ipp=1,nppmax
         if(c(ipp) .gt. ipoprem) ctmp(ipp) = c(ipp)-1
      enddo
c      write(*,*)'c=',c
c      write(*,*)'ctmp=',ctmp
c      write(*,*) 'fin de merge'
      end subroutine merging

 

******************************************************************
*     Mise a jour de c et f en cas d acceptation d'un split/merge
      subroutine accept5(nppmax,npopmax,nlocmax,nalmax,
     &     nal,c,ctmp,f,ftmp,drift,drifttmp)
      implicit none
      integer nppmax,npopmax,nlocmax,nalmax,
     &     nal(nlocmax),c(nppmax),ctmp(nppmax)
      double precision f(npopmax,nlocmax,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     drift(npopmax),drifttmp(npopmax)
      integer ipop,iloc,ial,ipp
c      write(*,*) 'debut de accept5'
c      write(*,*) 'f=',f
c      write(*,*) 'ftmp=',ftmp
      do ipp=1,nppmax
         c(ipp) = ctmp(ipp)
      enddo
      do ipop = 1,npopmax
         do iloc= 1,nlocmax
            do ial=1,nal(iloc)
               f(ipop,iloc,ial) = ftmp(ipop,iloc,ial)
            enddo
         enddo
         drift(ipop) = drifttmp(ipop)
      enddo

c      write(*,*) 'f=',f
c      write(*,*) 'fin de accept5'
      end subroutine accept5

*
*     coefficients du binome C_n^p
*
      double precision function bico(n,p)
      implicit none
      integer n,p
      double precision gglgamfn
      bico = dexp(gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
     &     gglgamfn(dble(n-p+1)))
c      write(*,*) 'in bico '
c$$$      write(*,*) 'n=', n
c$$$      write(*,*) 'p=', p
c$$$      write(*,*) 'gglgamfn(dble(n+1))=',gglgamfn(dble(n+1))
c$$$      write(*,*) 'gglgamfn(dble(p+1))=',gglgamfn(dble(p+1))
c$$$      write(*,*) 'gglgamfn(dble(n-p+1)))=',gglgamfn(dble(n-p+1))
c$$$      write(*,*) 'dexp()=',dexp(gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
c$$$     &     gglgamfn(dble(n-p+1)))
c$$$      write(*,*) 'bico =', nint(dexp(gglgamfn(dble(n+1))-
c$$$     &     gglgamfn(dble(p+1))-
c$$$     &     gglgamfn(dble(n-p+1))))
c      write(*,*) 'bico =',bico
      
      end function bico



*****************************************************************
*     ln du coefficient du binome C_n^p
*
      double precision function lbico(n,p)
      implicit none
      integer n,p
      double precision gglgamfn
      lbico = gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
     &     gglgamfn(dble(n-p+1))
      end function lbico




*
*     log du ratio des vraisemblances dans bdpop6
*
      double precision function llr6(z,f,ftmp,c,ctmp,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,z(nindiv,nlocmax2),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp

      llr6 = 0

*     log du rapport des vraisemblances
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
         do iloc=1,nloc
            ial1 = z(iindiv,2*iloc-1)
            ial2 = z(iindiv,2*iloc)
            if(ial1 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftmp(ipoptmp,iloc,ial1)) - 
     &              dlog(f(ipop,iloc,ial1))
            endif
            if(ial2 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftmp(ipoptmp,iloc,ial2)) - 
     &              dlog(f(ipop,iloc,ial2))
            endif
            
c$$$            if(ipoptmp .eq. 3) then 
c$$$               write(*,*) 'ipoptmp=',ipoptmp
c$$$               write(*,*) 'iinidiv=',iindiv
c$$$               write(*,*) 'c=', c(indcell(iindiv))
c$$$               write(*,*) 'ctmp=', ctmp(indcelltmp(iindiv))
c$$$               write(*,*) 'z=',z(iindiv,1),z(iindiv,2)
c$$$               write(*,*) 'llr6 =',llr6
c$$$               write(*,*) 'ftmp(ipoptmp,iloc,ial1)=',
c$$$     &              ftmp(ipoptmp,iloc,ial1) 
c$$$               write(*,*) 'f(ipop,iloc,ial1)=',
c$$$     &              f(ipop,iloc,ial1)
c$$$               write(*,*) 'ftmp(ipoptmp,iloc,ial2)=',
c$$$     &              ftmp(ipoptmp,iloc,ial2)
c$$$               write(*,*) 'f(ipop,iloc,ial2)=',
c$$$     &              f(ipop,iloc,ial2)
c$$$            endif

         enddo
      enddo
      if(ploidy .eq. 1) llr6 = 0.5d0*llr6 
      end function llr6


*
*     comptage des alleles dans chaque pop pour c
*
      subroutine countn(nindiv,nlocmax,nlocmax2,npopmax,
     &     nppmax,nal,nalmax,z,n,indcell,c,ploidy)
      implicit none
      integer  nindiv,nlocmax,nlocmax2,npopmax,nppmax,nalmax,ploidy,
     &     z(nindiv,nlocmax2),nal(nlocmax),
     &     n(npopmax,nlocmax,nalmax),c(nppmax),indcell(nindiv)
      integer ipop,iloc,ial,iindiv
*     init du tableau
      do ipop = 1,npopmax
         do iloc = 1,nlocmax
            do ial =1,nal(iloc)
               n(ipop,iloc,ial)=0
            enddo
         enddo
      enddo
*     comptage
      do iindiv = 1,nindiv
         do iloc = 1,nlocmax
            if(z(iindiv,2*iloc-1) .ne. -999) then
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1))+ 1 
            endif
            if(ploidy .eq. 2) then
               if(z(iindiv,2*iloc) .ne. -999) then 
                  n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) = 
     &                 n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) + 1 
               endif
            endif
         enddo
      enddo

      end subroutine countn



*
*     log du ratio (prob cond. complete)/prior
*     pour les frequences
*     dans un split de la pop ipop
      double precision function lrf(ipop,npopmax,nlocmax,nal,nalmax,
     &     f,fa,drift,n)
      implicit none
      integer ipop,npopmax,nlocmax,nal(nlocmax),nalmax,
     &     n(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),
     &     drift(npopmax)
      integer iloc,ial,nn
      double precision ss,gglgamfn,q

      lrf = 0.
      q = (1-drift(ipop))/drift(ipop)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ipop,iloc,ial))) +
     &           (1 - fa(iloc,ial) * q - dble(n(ipop,iloc,ial)))*
     &           dlog(f(ipop,iloc,ial))
            nn = nn + n(ipop,iloc,ial)
         enddo
c         write(*,*) 'nn=',nn
         lrf = lrf + gglgamfn(dble(nal(iloc))) -
     &        gglgamfn(q + nn) + ss

      enddo
      end function lrf



***********************************************************************
*
*     Naissance et mort de pop avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdpop9(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),z(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termf9
      double precision b,bern,ggrbinom      
       do ipp=1,nppmax
          cellpop(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ggrbinom(1.d0,0.5d0)
       
       if(b .eq. 1) then
          if(npop .lt. npopmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la pop qui split
             isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
             
*     recherche des cellules affectees a cette pop
             call who(c,isplit,npp,nppmax,cellpop,ncellpop)
             
             if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
                nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
                if(nu .gt. 0) then
                   
*     tirage des cellules reallouees
                   call sample2(cellpop,nppmax,nu,ncellpop,
     &                  listcell)
                   
*     proposition de reallocation dans la pop npop+1
                   call split(npop+1,c,ctmp,nppmax,nu,
     &                  listcell)
                else 
                   do ipp = 1,nppmax
                      ctmp(ipp) = c(ipp)
                   enddo
                endif
             else
                nu = 0
                do ipp = 1,nppmax
                   ctmp(ipp) = c(ipp)
                enddo
             endif

*     comptage des alleles sur chaque locus pour c puis ctmp
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nal,nalmax,z,n,indcell,c,ploidy)
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
             
*     proposition nouvelle freq et derive 
c     call addfreq5(isplit,npop,npopmax,nloc,nlocmax,
c     &     nal,nalmax,f,ftmp,fa,drift,drifttmp,a,ptmp)
             call addfreq7(npop,npopmax,nloc,nlocmax,
     &            nal,nalmax,isplit,
     &            f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
             
*     calcul du log du ratio
*     terme des vraisemblances
             lratio =  llr6(z,f,ftmp,c,ctmp,indcell,
     &            indcell,npopmax,nlocmax,nalmax,
     &            nindiv,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
c$$$             lratio = lratio + termfsplit(isplit,npop,npopmax,
c$$$     &            nlocmax,nal,nalmax,
c$$$     &            f,ftmp,n,ntmp,fa,drift,drifttmp) 
             lratio = lratio 
     &           + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,isplit) 
     &            -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &            drifttmp,isplit) 
     &            -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &            drifttmp,npop+1) 

             
*     terme des proposal sur c
             lratio = lratio + dlog(2*dble(ncellpop+1)) + 
     &            lbico(ncellpop,nu) - dlog(dble(npop+1)) 
             
*     terme des priors sur c
             lratio = lratio + 
     &            dble(npp)*(dlog(dble(npop)) - 
     &            dlog(dble(npop+1)))
             
             lratio = dmin1(0.d0,lratio)
             alpha = dexp(lratio)
             bern = ggrbinom(1.d0,alpha)

c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctmp=',ctmp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',ial,')=',
c$$$     &                           n(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ntmp(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',ial,')=',
c$$$     &                       f(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ftmp(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftmp,c,ctmp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nalmax,
c$$$     &                  nindiv,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,npop,npopmax,
c$$$     &                  nlocmax,nal,nalmax,
c$$$     &                  f,ftmp,n,ntmp,fa,drift,drifttmp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(npop)) - 
c$$$     &                  dlog(dble(npop+1)))
c                    write(*,*) 'alpha=',alpha 

             
             if(bern .eq. 1) then
                call accept5(nppmax,npopmax,nlocmax,
     &               nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                npop = npop + 1
             endif
          endif

*     merge
      else
         if(npop .gt. npopmin) then 
c             write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + 
     &              idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo

*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
 
*     next line corrected by Gilles on 05/01/08
*            if(ncellpop .gt. 0) then
*     proposition de reallocation dans la pop ipophost
               call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &              ncellpop,cellpop)

*     comptage des alleles sur chaque locus pour c puis ctmp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,z,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(ipoprem,ipophost,npop,npopmax,
c     &              nloc,nlocmax,nal,nalmax,f,ftmp,drift,drifttmp,
c     &              a,fa)
               call remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftmp,c,ctmp,indcell,
     &              indcell,npopmax,nlocmax,nalmax,
     &              nindiv,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
c$$$               lratio = lratio + 
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nal,nalmax,
c$$$     &              f,ftmp,n,ntmp,fa,drift,drifttmp)
            lratio = lratio  
     &     + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,ipophost) 
     &      + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,ipoprem) 
     & -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &              drifttmp,ipophost)


*     terme des proposal sur c
               lratio = lratio + dlog(dble(npop)) - 
     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop-1)))
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
      
         
c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctmp=',ctmp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',ial,')=',
c$$$     &                       n(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ntmp(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',ial,')=',
c$$$     &                       f(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ftmp(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftmp,c,ctmp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nalmax,
c$$$     &              nindiv,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nal,nalmax,
c$$$     &              f,ftmp,n,ntmp,fa,drift,drifttmp)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(npop)) - 
c$$$     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(npop)) - 
c$$$     &              dlog(dble(npop-1)))
c               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,npopmax,nlocmax,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop - 1
               endif
*            endif
         endif
      endif
      end subroutine bdpop9




***********************************************************************
*     split/merge populations in the spatial D-model
*     changes from bdpop7bis:
*     - process populations whatever the number of tiles or individuals
*       they have
      subroutine bdpop9bis(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),z(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termf9bis,
     &     gglgamfn
      integer b,iloc
      double precision bern,ggrbinom

c      write(*,*) ''

      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'isplit=',isplit

*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)

            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees 
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then                 

*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)
 
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,
     &                 listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)

            
*     proposition nouvelle freq et derive 
            call addfreq7bis(npop,npopmax,nloc,
     &           nlocmax,nal,nalmax,isplit,
     &           f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)

*     calcul du log du ratio
*     terme des vraisemblances
c            write(*,*) 'calcul du log du ratio'
            lratio =  llr6(z,f,ftmp,c,ctmp,indcell,
     &           indcell,npopmax,nlocmax,nalmax,
     &           nindiv,nloc,nlocmax2,nppmax,ploidy)

*     term proposal freq.
            lratio = lratio 
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,isplit) 
     &         - termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,isplit) 
     &         - termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,npop+1)  

* term prior freq
            do iloc = 1,nloc
               lratio = lratio + gglgamfn(dble(nal(iloc)))
            enddo

*     terme des proposal sur c
            lratio = lratio + dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))

            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)

            if(bern .eq. 1) then
               call accept5(nppmax,npopmax,nlocmax,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif

      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 
     &              + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)

*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)

*     propostion du nouveau tableau de freq et de derives
            call remfreq7bis(ipoprem,ipophost,
     &           npop,npopmax,nloc,nlocmax,nal,
     &           nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,
     &           ntmp)
            
*     calcul du log du ratio  
*     terme des vraisemblances
            lratio =  llr6(z,f,ftmp,c,ctmp,indcell,
     &           indcell,npopmax,nlocmax,nalmax,
     &           nindiv,nloc,nlocmax2,nppmax,ploidy)

            lratio = lratio  
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,ipophost) 
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,ipoprem) 
     &       -termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,ipophost)

* term prior freq
            do iloc = 1,nloc
               lratio = lratio - gglgamfn(dble(nal(iloc)))
            enddo
     
*     terme des proposal sur c
            lratio = lratio + dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))

            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha) 

            if(bern .eq. 1) then
               call accept5(nppmax,npopmax,nlocmax,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif 
      
      end subroutine bdpop9bis

 


*
*     ajoute une pop 
*     dans  le tableau des dérives selon le prior
*     et dans le tableau des frequences 
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq7(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,fa,
     &     drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax),ggrunif
      double precision bern,ggrbinom

c      write(*,*) 'debut addfreq7'
c      write(*,*) 'fa=',fa
c      write(*,*) 'drift=',drift
c      write(*,*) 'ntmp=',ntmp

*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         drifttmp(ipop) = drift(ipop)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttmp(isplit) = ggrunif(0.d0,1.d0)
      drifttmp(npop+1) = ggrunif(0.d0,1.d0)

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftmp(',isplit,1,1,')=',ftmp(isplit,1,1)
c$$$      write(*,*) 'ftmp(',isplit,1,2,')=',ftmp(isplit,1,2)
c$$$      write(*,*) 'f(',npop+1,1,1,')=',f(npop+1,1,1)
c$$$      write(*,*) 'f(',npop+1,1,2,')=',f(npop+1,1,2)
c$$$      write(*,*) 'ftmp(',npop+1,1,1,')=',ftmp(npop+1,1,1)
c$$$      write(*,*) 'ftmp(',npop+1,1,2,')=',ftmp(npop+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftmp(',isplit,2,1,')=',ftmp(isplit,2,1)
c$$$      write(*,*) 'ftmp(',isplit,2,2,')=',ftmp(isplit,2,2)
c$$$      write(*,*) 'f(',npop+1,2,1,')=',f(npop+1,2,1)
c$$$      write(*,*) 'f(',npop+1,2,2,')=',f(npop+1,2,2)
c$$$      write(*,*) 'ftmp(',npop+1,2,1,')=',ftmp(npop+1,2,1)
c$$$      write(*,*) 'ftmp(',npop+1,2,2,')=',ftmp(npop+1,2,2)
c$$$

      end subroutine addfreq7


*
*     ajoute une pop 
*     dans le tableau des frequences 
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq8(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,fa,
     &     drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax),ggrunif
*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo
*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo
      end subroutine addfreq8

*
*     ajoute une pop 
*     dans  le tableau des frequences  selon le prior
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     et une valeur 0.5d0 dans le tableau des dérives 
*     (sans modifier les tableaux en entrée)
*     pour court-circuiter le F-model
      subroutine addfreq7bis(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,
     &     fa,drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax)


*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         drifttmp(ipop) = drift(ipop)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttmp(isplit) = 0.5d0
      drifttmp(npop+1) = 0.5d0

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftmp(',isplit,1,1,')=',ftmp(isplit,1,1)
c$$$      write(*,*) 'ftmp(',isplit,1,2,')=',ftmp(isplit,1,2)
c$$$      write(*,*) 'f(',npop+1,1,1,')=',f(npop+1,1,1)
c$$$      write(*,*) 'f(',npop+1,1,2,')=',f(npop+1,1,2)
c$$$      write(*,*) 'ftmp(',npop+1,1,1,')=',ftmp(npop+1,1,1)
c$$$      write(*,*) 'ftmp(',npop+1,1,2,')=',ftmp(npop+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftmp(',isplit,2,1,')=',ftmp(isplit,2,1)
c$$$      write(*,*) 'ftmp(',isplit,2,2,')=',ftmp(isplit,2,2)
c$$$      write(*,*) 'f(',npop+1,2,1,')=',f(npop+1,2,1)
c$$$      write(*,*) 'f(',npop+1,2,2,')=',f(npop+1,2,2)
c$$$      write(*,*) 'ftmp(',npop+1,2,1,')=',ftmp(npop+1,2,1)
c$$$      write(*,*) 'ftmp(',npop+1,2,2,')=',ftmp(npop+1,2,2)
c$$$

      end subroutine addfreq7bis



*     enleve une pop des tableau des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux pops
*     tirage d'une derive selon prior 
*     sans modifier des tableaux en entrée
      subroutine remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax),ggrunif

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
         drifttmp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pop
      drifttmp(ipophost) = ggrunif(0.d0,1.d0)
      
*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftmp(',ipophost,1,1,')=',ftmp(ipophost,1,1)
c$$$      write(*,*) 'ftmp(',ipophost,1,2,')=',ftmp(ipophost,1,2)
c$$$      write(*,*) 'ftmp(',ipophost,2,1,')=',ftmp(ipophost,2,1)
c$$$      write(*,*) 'ftmp(',ipophost,2,2,')=',ftmp(ipophost,2,2)

      end subroutine remfreq7



******************************************************************
*     enleve une pop des tableau des frequences
*     tirage d'une freq selon posterior apres un merge de deux pops
*     sans modifier des tableaux en entrée
      subroutine remfreq8(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax),ggrunif

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
      enddo

*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo
      end subroutine remfreq8



***********************************************
*     enleve une pop des tableau des derives
*     sans modifier des tableaux en entrée
      subroutine remdrift(ipoprem,ipophost,npop,npopmax,drift,drifttmp,
     &     shape1,shape2)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax
      double precision drift(npopmax),drifttmp(npopmax),shape1,shape2
      integer ipop
      double precision ggrbet
*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         drifttmp(ipop) = -999
      enddo
*     derive pour la nouvelle pop
      drifttmp(ipophost) =  ggrbet(shape1,shape2)
      end subroutine remdrift


***********************************************
*     enleve une pop des tableau des derives
*     sans modifier des tableaux en entrée
*     d* = (d1+d2)/2
      subroutine remdrift2(ipoprem,ipophost,npop,npopmax,drift,drifttmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax
      double precision drift(npopmax),drifttmp(npopmax)
      integer ipop
      double precision ggrbet
*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         drifttmp(ipop) = -999
      enddo
*     derive pour la nouvelle pop
      drifttmp(ipophost) = (drift(ipophost) + drift(ipoprem))/2
      end subroutine remdrift2



*****************************************************************
*     enleve une pop des tableaux des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux pops
*     la nouvelle derive est mise à 0.5d0
*     (sans modifier les tableaux en entrée)
*     c'est pour court-circuiter le F-model 
      subroutine remfreq7bis(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax)

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
         drifttmp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pop
      drifttmp(ipophost) = 0.5d0
      
*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftmp(',ipophost,1,1,')=',ftmp(ipophost,1,1)
c$$$      write(*,*) 'ftmp(',ipophost,1,2,')=',ftmp(ipophost,1,2)
c$$$      write(*,*) 'ftmp(',ipophost,2,1,')=',ftmp(ipophost,2,1)
c$$$      write(*,*) 'ftmp(',ipophost,2,2,')=',ftmp(ipophost,2,2)

      end subroutine remfreq7bis



************************************************************************
*     terme des freq dans le log ratio pour un split 
      double precision function termfsplit(isplit,npop,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     isplit,npop
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfsplit  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(isplit,iloc,ial))) +
     &           dble(n(isplit,iloc,ial))*
     &           dlog(f(isplit,iloc,ial))
            nn = nn + n(isplit,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttmp(isplit))/drifttmp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(isplit,iloc,ial))) +
     &           dble(ntmp(isplit,iloc,ial))*
     &           dlog(ftmp(isplit,iloc,ial))
            nn = nn + ntmp(isplit,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttmp(npop+1))/drifttmp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(npop+1,iloc,ial))) +
     &           dble(ntmp(npop+1,iloc,ial))*
     &           dlog(ftmp(npop+1,iloc,ial))
            nn = nn + ntmp(npop+1,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit - tt
      end function termfsplit



*******************************************************
*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5d0
      double precision function termfsplitbis(isplit,npop,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     isplit,npop
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfsplitbis  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(isplit,iloc,ial))) +
     &           dble(n(isplit,iloc,ial))*
     &           dlog(f(isplit,iloc,ial))
            nn = nn + n(isplit,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttmp(isplit))/drifttmp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(isplit,iloc,ial))) +
     &           dble(ntmp(isplit,iloc,ial))*
     &           dlog(ftmp(isplit,iloc,ial))
            nn = nn + ntmp(isplit,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttmp(npop+1))/drifttmp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(npop+1,iloc,ial))) +
     &           dble(ntmp(npop+1,iloc,ial))*
     &           dlog(ftmp(npop+1,iloc,ial))
            nn = nn + ntmp(npop+1,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

      end function termfsplitbis


*********************************************************************
*     terme des freq dans le log ratio pour un split 
      double precision function termfmerge(ihost,irem,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     ihost,irem
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfmerge  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ihost,iloc,ial))) +
     &           dble(n(ihost,iloc,ial))*
     &           dlog(f(ihost,iloc,ial))
            nn = nn + n(ihost,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(irem,iloc,ial))) +
     &           dble(n(irem,iloc,ial))*
     &           dlog(f(irem,iloc,ial))
            nn = nn + n(irem,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttmp(ihost))/drifttmp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(ihost,iloc,ial))) +
     &           dble(ntmp(ihost,iloc,ial))*
     &           dlog(ftmp(ihost,iloc,ial))
            nn = nn + ntmp(ihost,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge - tt
      end function termfmerge


*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5d0
      double precision function termfmergebis(ihost,irem,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     ihost,irem
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfmergebis  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ihost,iloc,ial))) +
     &           dble(n(ihost,iloc,ial))*
     &           dlog(f(ihost,iloc,ial))
            nn = nn + n(ihost,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(irem,iloc,ial))) +
     &           dble(n(irem,iloc,ial))*
     &           dlog(f(irem,iloc,ial))
            nn = nn + n(irem,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttmp(ihost))/drifttmp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(ihost,iloc,ial))) +
     &           dble(ntmp(ihost,iloc,ial))*
     &           dlog(ftmp(ihost,iloc,ial))
            nn = nn + ntmp(ihost,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis - tt
      end function termfmergebis



*******************************************************
*     term from proposal of frequencies in a split in bdpop9bis
      double precision function termf9(npopmax,nloc,nal,nalmax,n,
     &     f,fa,drift,ipop)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ipop,ipop2
      double precision f,fa,drift
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     f(npopmax,nloc,nalmax),fa(nloc,nalmax),drift(npopmax)
      integer iloc,ial,nn
      double precision gglgamfn,q
*     next line corrected by gilles on 5/01/08
*      q = drift(ipop) / (1- drift(ipop))
      q = (1-drift(ipop)) / drift(ipop)
      termf9 = 0.d0
      do iloc=1,nloc
         nn = 0
         do ial=1,nal(iloc)
            termf9 = termf9 + 
     &           dble(n(ipop,iloc,ial))*dlog(f(ipop,iloc,ial)) - 
     &           gglgamfn(q*fa(iloc,ial)+dble(n(ipop,iloc,ial))) + 
     &           gglgamfn(q*fa(iloc,ial))
            nn = nn + n(ipop,iloc,ial)
         enddo
*     next line corrected by gilles on 5/01/08
*         termf9 = termf9 + 
*     &        gglgamfn(q+dble(nn)) -  gglgamfn(dble(nn))
         termf9 = termf9 + 
     &        gglgamfn(q+dble(nn)) - gglgamfn(q)
      enddo
      end function termf9

 
*******************************************************
*     term from proposal of frequencies in a split in bdpop9bis
      double precision function termf9bis(npopmax,nloc,nal,nalmax,n,
     &     f,ipop)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ipop,ipop2
      double precision f
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     f(npopmax,nloc,nalmax)
      integer iloc,ial,nn
      double precision gglgamfn
      termf9bis = 0.d0
      do iloc=1,nloc
         nn = 0
         do ial=1,nal(iloc)
            termf9bis = termf9bis + 
     &           dble(n(ipop,iloc,ial))*dlog(f(ipop,iloc,ial)) - 
     &           gglgamfn(1+dble(n(ipop,iloc,ial))) 
            nn = nn + n(ipop,iloc,ial)
         enddo
         termf9bis = termf9bis + 
     &        gglgamfn(dble(nal(iloc)+nn))
      enddo
      end function termf9bis

 

********************************************************
*     log vraisemblance
      double precision function ll(z,nindiv,nlocmax,nlocmax2,npopmax,
     &     nalmax,nppmax,c,f,indcell)
      implicit none
      integer nindiv,nlocmax,nlocmax2,npopmax,
     &     z(nindiv,nlocmax2),nppmax,c(nppmax),nalmax,
     &     indcell(nindiv)
      double precision f(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop
      ll = 0
      do iindiv = 1,nindiv
         ipop = c(indcell(iindiv))
         do iloc = 1,nlocmax
            ial1 = z(iindiv,2*iloc-1)
            ial2 = z(iindiv,2*iloc)
            if(ial1 .ne. -999) then 
               ll = ll + dlog(f(ipop,iloc,ial1))
            endif
            if(ial2 .ne. -999) then 
               ll = ll + dlog(f(ipop,iloc,ial2))
            endif
            if((ial1 .ne. ial2))  then 
               ll = ll + dlog(2.d0)
            endif
         enddo
      enddo
      end function ll

 
*******************************************************************
*     log de la proba a posteriori du vecteur de parametres
      double precision function lpp(lambdamax,lambda,z,npop,npp,drift,f,
     &     fa,c,nppmax,
     &     nindiv,nlocmax2,npopmax,nlocmax,nalmax,indcell,nal,
     &     fmodel,xlim,ylim,shape1,shape2)
      implicit none
      integer nindiv,nlocmax2,npop,npopmax,nlocmax,nalmax,
     &     npp,nppmax,z(nindiv,nlocmax2),indcell(nindiv),
     &     c(nppmax),nal(nlocmax),fmodel
      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),lambdamax,lambda,xlim(2),ylim(2)
      integer ipp,ipop,iloc,ial
      double precision gglgamfn,shape1,shape2,ll,lg
*     contrib npop
      lpp = lpp - dlog(dble(npopmax))
*     contrib lambda
      lpp = - dlog(lambdamax)
*     contrib npp
      lpp = lpp - lambda + dble(npp)*dlog(lambda) - 
     &     gglgamfn(dble(npp+1))
*     contrib c
      lpp = lpp - dble(npp)*dlog(dble(npop))
*     contrib u
      lpp = lpp - dble(npp)*dlog((xlim(2)-xlim(1))*(ylim(2)-ylim(1)))
*     contrib freq
      if(fmodel .eq. 0) then
         lg = 0
         do iloc = 1,nlocmax
            lg = lg + gglgamfn(dble(nal(iloc)))
         enddo
         lpp= lpp + dble(npop) * lg
      endif
      if(fmodel .eq. 1) then
*     contrib ancestral freq
         lg = 0
         do iloc = 1,nlocmax
            lg = lg + gglgamfn(dble(nal(iloc)))
         enddo
         lpp= lpp + lg
*     contrib drifts
         do ipop = 1,npop
            lpp = lpp + 
     &           gglgamfn(shape1+shape2) - gglgamfn(shape1) 
     &           - gglgamfn(shape2) + (shape1-1)*dlog(drift(ipop)) + 
     &           (shape2-1)*dlog(1-drift(ipop))
         enddo
*     contrib freq
         do ipop = 1,npop
            do iloc = 1,nlocmax
               lpp = lpp + gglgamfn((1-drift(ipop))/drift(ipop))
               do ial=1,nal(iloc)
                  lpp = lpp -gglgamfn(fa(iloc,ial)*
     &                 (1-drift(ipop))/drift(ipop)) + 
     &                 (fa(iloc,ial)* 
     &                 (1-drift(ipop))/drift(ipop)-1)*
     &              dlog(f(ipop,iloc,ial))
               enddo
            enddo
         enddo
      endif
      
*     contrib likelihood
      lpp = lpp + ll(z,nindiv,nlocmax,nlocmax2,npopmax,
     &     nalmax,nppmax,c,f,indcell)

      end function lpp





***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,z)
*     drifts and ancestral are not changed
      subroutine samplef(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,ntmp,
     &     ipop1,ipop2
      double precision ftmp,a
      dimension nal(nlocmax),
     &     ntmp(npopmax,nloc,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax)
      integer iloc,ial,ipop
      double precision ptmp
      double precision alpha

c      write(*,*) 'debut sample f'

*     new freq for pop ipop1
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = alpha + dble(ntmp(ipop1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

*     for pop ipop2
      if(ipop2 .ne. ipop1) then
         do iloc=1,nloc
            do ial = 1,nal(iloc)
               a(ial) = 1. + dble(ntmp(ipop2,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop2,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
c      write(*,*) 'end  sample f'
      end subroutine samplef 


***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,z)
*     drifts and ancestral are not changed
*     if CFM for allele frequencies
      subroutine samplef2(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,ipop1,ipop2,f,ftmp,
     &     fa,drift,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,ntmp,
     &     ipop1,ipop2
      double precision f,drift,ftmp,fa,a
      dimension nal(nlocmax),
     &     ntmp(npopmax,nloc,nalmax),
     &     f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),a(nalmax),ptmp(nalmax)
      integer iloc,ial,ipop
      double precision ptmp

*     init ftmp 
c$$$      do ipop = 1,npop
c$$$         do iloc=1,nloc
c$$$            do ial=1,nal(iloc)
c$$$               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
c$$$            enddo 
c$$$         enddo
c$$$      enddo

      do iloc=1,nloc
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
         enddo 
      enddo
      do iloc=1,nloc
         do ial=1,nal(iloc)
            ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
         enddo 
      enddo


*     new freq
*     for pop ipop1
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drift(ipop1))/
     &           drift(ipop1)+dble(ntmp(ipop1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

*     for pop ipop2
      if(ipop2 .ne. ipop1) then
         do iloc=1,nloc
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drift(ipop2))/
     &              drift(ipop2)+dble(ntmp(ipop2,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop2,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      end subroutine samplef2 



***********************************************************************
*     log of ratio of proposals in a joint update of c and f
*     (in subroutine udcf)
      double precision function lrpf(npopmax,nloc,nal,nalmax,n,
     &     ntmp,f,ftmp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ntmp,ipop1,ipop2
      double precision f,ftmp
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     ntmp(npopmax,nloc,nalmax),f(npopmax,nloc,nalmax),
     &     ftmp(npopmax,nloc,nalmax)
      integer iloc,ial,n1,n2,ntmp1,ntmp2
      double precision gglgamfn
      double precision alpha
      alpha = 1

      lrpf = 0.d0
      do iloc=1,nloc
c     write(*,*) 'iloc=',iloc
         n1 = 0
         n2 = 0
         ntmp1 = 0
         ntmp2 = 0
         do ial=1,nal(iloc)
c     write(*,*) 'ial=',ial
            lrpf = lrpf + 
     &     (n(ipop1,iloc,ial)+alpha-1)*dlog(f(ipop1,iloc,ial)) 
     &   + (n(ipop2,iloc,ial)+alpha-1)*dlog(f(ipop2,iloc,ial))
     &   - (ntmp(ipop1,iloc,ial)+alpha-1)*dlog(ftmp(ipop1,iloc,ial)) 
     &   - (ntmp(ipop2,iloc,ial)+alpha-1)*dlog(ftmp(ipop2,iloc,ial)) 
     &           - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &           - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &           + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &           + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
            n1 = n1 + n(ipop1,iloc,ial)
            n2 = n2 + n(ipop2,iloc,ial)
            ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
            ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
c     write(*,*) 'lrpf=',lrpf
         enddo
         lrpf = lrpf + 
     &          gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &        + gglgamfn(alpha*dble(nal(iloc))+dble(n2)) 
     &        - gglgamfn(alpha*dble(nal(iloc))+dble(ntmp1)) 
     &        - gglgamfn(alpha*dble(nal(iloc))+dble(ntmp2))
c     write(*,*) 'lrpf=',lrpf
      enddo
c$$$      if(ipop1 .ne. ipop2) then
c$$$         do iloc=1,nloc
c$$$c         write(*,*) 'iloc=',iloc
c$$$            n1 = 0
c$$$            n2 = 0
c$$$            ntmp1 = 0
c$$$            ntmp2 = 0
c$$$            do ial=1,nal(iloc)
c$$$c     write(*,*) 'ial=',ial
c$$$               lrpf = lrpf + 
c$$$     &              n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
c$$$     &            + n(ipop2,iloc,ial)*dlog(f(ipop2,iloc,ial))
c$$$     &             - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
c$$$     &             - ntmp(ipop2,iloc,ial)*dlog(ftmp(ipop2,iloc,ial)) 
c$$$     &              - gglgamfn(1+dble(n(ipop1,iloc,ial)))
c$$$     &              - gglgamfn(1+dble(n(ipop2,iloc,ial)))
c$$$     &              + gglgamfn(1+dble(ntmp(ipop1,iloc,ial)))
c$$$     &              + gglgamfn(1+dble(ntmp(ipop2,iloc,ial)))
c$$$               n1 = n1 + n(ipop1,iloc,ial)
c$$$               n2 = n2 + n(ipop2,iloc,ial)
c$$$               ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
c$$$               ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$            enddo
c$$$            lrpf = lrpf + gglgamfn(dble(nal(iloc)+n1)) 
c$$$     &           + gglgamfn(dble(nal(iloc)+n2)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp1)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp2))
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$         enddo
c$$$      else
c$$$         do iloc=1,nloc
c$$$c         write(*,*) 'iloc=',iloc
c$$$            n1 = 0
c$$$            ntmp1 = 0
c$$$            do ial=1,nal(iloc)
c$$$c     write(*,*) 'ial=',ial
c$$$               lrpf = lrpf + 
c$$$     &              n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
c$$$     &             - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
c$$$     &             - gglgamfn(1+dble(n(ipop1,iloc,ial)))
c$$$     &             + gglgamfn(1+dble(ntmp(ipop1,iloc,ial)))
c$$$               n1 = n1 + n(ipop1,iloc,ial)
c$$$               ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$            enddo
c$$$            lrpf = lrpf + gglgamfn(dble(nal(iloc)+n1)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp1)) 
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$         enddo
c$$$      endif                     
      end 




***********************************************************************
*
*     log of ratio of proposals x priors in a joint update of c and f
*     (in subroutine udcf2)
*     
      double precision function lrppf(npopmax,nloc,nal,nalmax,n,
     &     ntmp,fa,drift,f,ftmp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ntmp,ipop1,ipop2
      double precision fa,drift,f,ftmp
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     ntmp(npopmax,nloc,nalmax),f(npopmax,nloc,nalmax),
     &     fa(nloc,nalmax),ftmp(npopmax,nloc,nalmax),
     &     drift(npopmax)
      integer iloc,ial,n1,n2,ntmp1,ntmp2
      double precision gglgamfn,q1,q2
      lrppf = 0.d0
      q1 = drift(ipop1)/(1-drift(ipop1))
      q2 = drift(ipop2)/(1-drift(ipop2))
      
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntmp1 = 0
         ntmp2 = 0
         do ial=1,nal(iloc)
            lrppf = lrppf + 
     &         n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
     &       + n(ipop2,iloc,ial)*dlog(f(ipop2,iloc,ial) )
     &       - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
     &       - ntmp(ipop2,iloc,ial)*dlog(ftmp(ipop2,iloc,ial)) 
     &       + gglgamfn(fa(iloc,ial)*q1+
     &           dble(ntmp(ipop1,iloc,ial)))
     &       + gglgamfn(fa(iloc,ial)*q2+
     &           dble(ntmp(ipop2,iloc,ial)))
     &       - gglgamfn(fa(iloc,ial)*q1+
     &        dble(n(ipop1,iloc,ial)))
     &       - gglgamfn(fa(iloc,ial)*q2+
     &        dble(n(ipop2,iloc,ial)))
            n1 = n1 + n(ipop1,iloc,ial)
            n2 = n2 + n(ipop2,iloc,ial)
            ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
            ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
         enddo
         lrppf = lrppf + gglgamfn(q1+dble(n1)) 
     &           + gglgamfn(q2+dble(n2)) 
     &           - gglgamfn(q1+dble(ntmp1)) 
     &           - gglgamfn(q2+dble(ntmp2))
      enddo
      end function lrppf



************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     MH ratio does not depend on proposed frequencies
*     D-model
      subroutine smd(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),z(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTfd
      integer b,iloc
      double precision llr6,termf9bis,gglgamfn
      double precision bern,ggrbinom
      
c      write(*,*) 'debut smd' 
c      write(*,*) 'npop =',npop

      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 

c            write(*,*) 'naissance'
            
*     split
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
C     ligne suivante corrigee le 17/01/08
C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c               write(*,*) 'nu=',nu
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c                   write(*,*) 'listcell=',listcell
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
            
*     calcul du log du ratio
*     terme des proposal sur c
            lratio = dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
            
*     terme des priors sur c
            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))
            
c            write(*,*) 'in smd   c= ',c(1),c(2)
c            write(*,*) '         ctmp= ',ctmp(1),ctmp(2)
c            write(*,*) '       Rc=',lratio 
            
*     terme des frequences
            lratio = lratio +
     &           lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  +
     &           lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) -
     &           lTfd(isplit,n,npopmax,nloc,nal,nalmax) 

c            write(*,*) 'in smd Tf=',
c     &           lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  
c     &           + lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) 
c     &           - lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'in smd lratio=',lratio 

c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
             
            if(bern .eq. 1) then
*     proposition nouvelle freq et derive 
               call addfreq7bis(npop,npopmax,nloc,nloc,
     &              nal,nalmax,isplit,
     &              f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif
         
*     merge
      else
         if(npop .gt. npopmin) then 
c            write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
            
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
            
*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c            write(*,*) 'terme en c lratio =',lratio
*     term des freq
            lratio = lratio + 
     &           lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax) -  
     &           lTfd(ipophost,n,npopmax,nloc,nal,nalmax) - 
     &           lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'terme en f=',
c     &           lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax)  
c     &           - lTfd(ipophost,n,npopmax,nloc,nal,nalmax) 
c     &           - lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'lratio =',lratio
             lratio = dmin1(0.d0,lratio)
             alpha  = dexp(lratio)
             bern = ggrbinom(1.d0,alpha)      
                          
             if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
                call remfreq7bis(ipoprem,ipophost,
     &               npop,npopmax,nloc,nloc,nal,
     &               nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
                call accept5(nppmax,npopmax,nloc,
     &               nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                npop = npop - 1
             endif
          endif
       endif

c       write(*,*) 'fin smd' 
       end subroutine smd
      

************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     CFM
*     d* drawn form prior
      subroutine sm(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy,shape1,shape2)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),z(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax),shape1,shape2
      integer ipoprem,ipp,isplit,ipop,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTf
      integer b,iloc
      double precision gglgamfn,ggrbet,rr
      double precision bern,ggrbinom

c         write(*,*) 'debut sm' 
      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
C     ligne suivante corrigee le 17/01/08
C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
*     nouvelles derives      
            do ipop = 1,npop
               drifttmp(ipop) = drift(ipop)
            enddo
            drifttmp(isplit) = ggrbet(shape1,shape2)
            drifttmp(npop+1) = ggrbet(shape1,shape2)
*     terme des proposal sur c
            lratio = dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))
*     terme des frequences
            lratio = lratio + 
     &          lTf(isplit,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  +
     &          lTf(npop+1,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) -
     &          lTf(isplit,n,fa,drift,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
            if(bern .eq. 1) then
*     proposition nouvelle freq 
               call addfreq8(npop,npopmax,nloc,nloc,
     &              nal,nalmax,isplit,
     &              f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif
         
*     merge
      else
         if(npop .gt. npopmin) then 
c            write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
*     enleve une pop dans le tableau tmporaire des derives            
            call remdrift(ipoprem,ipophost,npop,npopmax,drift,
     &           drifttmp,shape1,shape2)
*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c     write(*,*) 'term en c lratio =',lratio
*     term des freq
            lratio = lratio + 
     &         lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) -  
     &         lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) - 
     &         lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'term en f=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drifttmp(ipophost)=',drifttmp(ipophost)
c$$$            write(*,*) 'term en f*_k1=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drift(ipophost)=',drift(ipophost)
c$$$             write(*,*) 'term en f_k1=',
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$             write(*,*) 'drift(ipoprem)=',drift(ipoprem)
c$$$            write(*,*) 'term en f_k2=',
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'lratio =',lratio
c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)      
                          
            if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
               call remfreq8(ipoprem,ipophost,
     &              npop,npopmax,nloc,nloc,nal,
     &              nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif
      end subroutine sm
      

************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     CFM
*     (d1*,d2*) = d1-u,d2+u
      subroutine sm2(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy,shape1,shape2)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),z(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax),shape1,shape2
      integer ipoprem,ipp,isplit,ipop,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTf
      integer b,iloc
      double precision gglgamfn,rr,deltad
      double precision bern,ggrbinom,ggrnorm

      deltad = shape1/(shape1+shape2)
C      deltad = .2

c         write(*,*) 'debut sm' 
      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
C     ligne suivante corrigee le 17/01/08
C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
*     nouvelles derives      
            do ipop = 1,npop
               drifttmp(ipop) = drift(ipop)
            enddo
c     next line modified 29/02/08
c            rr = ggrunif(0.d0,1.d0)*deltad
            rr = ggrnorm(0.d0,1.d0)*deltad
            drifttmp(isplit) = drift(isplit) - rr
            drifttmp(npop+1) = drift(isplit) + rr
            if(((drifttmp(isplit)   .gt. 1d-300) .and. 
     &          (1-drifttmp(isplit) .gt. 1d-300)).and.
     &         ((drifttmp(npop+1)   .gt. 1d-300) .and. 
     &          (1-drifttmp(npop+1) .gt. 1d-300))) then

*     terme des proposal sur c
               lratio = dlog(2*dble(ncellpop+1)) + 
     &              lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
               lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop+1)))
*     terme des frequences
               lratio = lratio + 
     &         lTf(isplit,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  +
     &         lTf(npop+1,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) -
     &         lTf(isplit,n,fa,drift,npopmax,nloc,nal,nalmax) 
*     term proposal drift
               lratio = lratio + dlog(2*deltad) 
     &              + .5*rr**2 + dlog(dsqrt(2*3.141593d0))
*     term prior drift
               lratio = lratio 
     &              + (shape1-1)*dlog(drifttmp(isplit))
     &              + (shape2-1)*dlog(1-drifttmp(isplit))
     &              + (shape1-1)*dlog(drifttmp(npop+1))
     &              + (shape2-1)*dlog(1-drifttmp(npop+1))
     &              - (shape1-1)*dlog(drift(isplit))
     &              - (shape2-1)*dlog(1-drift(isplit))
     &              + gglgamfn(shape1+shape2)
     &              - gglgamfn(shape1)-gglgamfn(shape2)
c            write(*,*) 'lratio=',lratio
            
               lratio = dmin1(0.d0,lratio)
               alpha = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)
             
               if(bern .eq. 1) then
*     proposition nouvelle freq 
                  call addfreq8(npop,npopmax,nloc,nloc,
     &                 nal,nalmax,isplit,
     &                 f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
                  call accept5(nppmax,npopmax,nloc,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop + 1
               endif
            endif
         endif
         
*     merge
      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
c            if(abs(drift(ipoprem)-drift(ipophost)) .lt. 2*deltad) then
*     on range dans la pop d'indice le plus petit
               if(ipophost .gt. ipoprem) then
                  ii = ipophost
                  ipophost = ipoprem
                  ipoprem = ii
               endif
               
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'drift=',drift

*     recherche des cellules qui vont etre reallouees
               call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
               
*     recherche des cellules de la pop hote
               call who(c,ipophost,npp,nppmax,cellpophost,
     &              ncellpophost)
               
*     proposition de reallocation dans la pop ipophost
               call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &              ncellpop,cellpop)
               
*     comptage des alleles sur chaque locus pour c puis ctmp
               call countn(nindiv,nloc,nloc2,npopmax,
     &              nppmax,nal,nalmax,z,n,indcell,c,ploidy)
               call countn(nindiv,nloc,nloc2,npopmax,
     &              nppmax,nal,nalmax,z,ntmp,indcell,ctmp,ploidy)
               
*     enleve une pop dans le tableau tmporaire des derives            
               call remdrift2(ipoprem,ipophost,npop,npopmax,drift,
     &           drifttmp)
               rr = (drift(ipophost) - drift(ipoprem))/(2*deltad)
c$$$               write(*,*) 'drifttmp=',drifttmp

*     calcul du log du ratio  
*     terme des proposal sur c
               lratio = dlog(dble(npop)) - 
     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 
            
*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop-1)))
c$$$               write(*,*) 'term en c lratio =',lratio
*     term des freq
               lratio = lratio + 
     &        lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) -  
     &        lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) -
     &        lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
*     terme proposal d
               lratio = lratio - dlog(2*deltad)
     &              - .5*rr**2 - dlog(dsqrt(2*3.141593d0))
c$$$               write(*,*) 'term en d lratio =',- dlog(2*deltad)
c$$$     &              - .5*rr**2 - dlog(dsqrt(2*3.141593d0))

*     term prior drift
               lratio = lratio 
     &              + (shape1-1)*dlog(drifttmp(ipophost))
     &              + (shape2-1)*dlog(1-drifttmp(ipophost))
     &              - (shape1-1)*dlog(drift(ipophost))
     &              - (shape2-1)*dlog(1-drift(ipophost))
     &              - (shape1-1)*dlog(drift(ipoprem))
     &              - (shape2-1)*dlog(1-drift(ipoprem))
     &              - gglgamfn(shape1+shape2)
     &              + gglgamfn(shape1)+gglgamfn(shape2)               

c$$$            write(*,*) 'term en f=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drifttmp(ipophost)=',drifttmp(ipophost)
c$$$            write(*,*) 'term en f*_k1=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drift(ipophost)=',drift(ipophost)
c$$$             write(*,*) 'term en f_k1=',
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$             write(*,*) 'drift(ipoprem)=',drift(ipoprem)
c$$$            write(*,*) 'term en f_k2=',
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'lratio =',lratio
 
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
                          
               if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
                  call remfreq8(ipoprem,ipophost,
     &                 npop,npopmax,nloc,nloc,nal,
     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
                  call accept5(nppmax,npopmax,nloc,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop - 1
               endif
c            endif
         endif
       endif
       end subroutine sm2
      


***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under CFM
      double precision function lTf(ipop,n,fa,drift,npopmax,nloc,nal,
     &     nalmax)
      implicit none
      integer ipop,n,npopmin,npopmax,nloc,nal,nalmax
      double precision fa,drift
      dimension n(npopmax,nloc,nalmax),fa(nloc,nalmax),nal(nloc), 
     &     drift(npopmax)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTf = 0
      do iloc = 1,nloc
         sumn = 0
         do ial = 1,nal(iloc)
            lTf = lTf + gglgamfn(dble(n(ipop,iloc,ial)) + 
     &           fa(iloc,ial)*(1-drift(ipop))/drift(ipop)) - 
     &           gglgamfn(fa(iloc,ial)*(1-drift(ipop))/drift(ipop))
            sumn = sumn + n(ipop,iloc,ial)
         enddo
         lTf = lTf + gglgamfn((1-drift(ipop))/drift(ipop)) - 
     &        gglgamfn(sumn+(1-drift(ipop))/drift(ipop))
      enddo
      end function lTf



***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under the Dirichlet model
      double precision function lTfd(ipop,n,npopmax,nloc,nal,nalmax)
      implicit none
      integer ipop,n,npopmin,npopmax,nloc,nal,nalmax
      dimension n(npopmax,nloc,nalmax),nal(nloc)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTfd = 0
      do iloc = 1,nloc
         sumn = 0
         do ial = 1,nal(iloc)
            lTfd = lTfd + gglgamfn(dble(n(ipop,iloc,ial)) + 1)
            sumn = sumn + n(ipop,iloc,ial)
         enddo
         lTfd = lTfd + gglgamfn(dble(nal(iloc))) - 
     &        gglgamfn(sumn+dble(nal(iloc)))
      enddo
      end function lTfd


***********************************************************************
      subroutine  postprocesschain2(nxdommax,nydommax,burnin,ninrub,
     &     npopmax,nppmax,nindiv,nloc,nal,nalmax,xlim,ylim,dt,nit,
     &     thinning,filenpop,filenpp,fileu,filec,filef,fileperm,filedom,
     &     s,u,c,f,pivot,fpiv,dom,coorddom,indcel,distcel,
     &     order,ordertmp,npopest)
      implicit none
      character*255 fileu,filec,filenpp,filenpop,filedom,filef,fileperm      
      integer nit,thinning,npp,npop,iit,nindiv,nxdommax,
     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indcel,
     &     ipop,nloc,nal,nalmax,ijunk,order,ordertmp,ipopperm,burnin,
     &     ninrub,npopest,nnit,iloc,ial,pivot
      double precision s,u,xlim,ylim,coorddom,dom,domperm,distcel,dt,
     &     f,fpiv
      integer iitsub
*     dimensionnement 
      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
     &     domperm(nxdommax*nydommax,npopmax),
     &     coorddom(2,nxdommax*nydommax),indcel(nxdommax*nydommax),
     &     distcel(nxdommax*nydommax),order(npopmax),ordertmp(npopmax),
     &     f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),nal(nloc)
      open(9,file=filenpop)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      open(13,file=filef)
      open(14,file=fileperm)
      open(15,file=filedom)

c      write(6,*) 'debut postproc order=',order

c      write(6,*) 'npopest=', npopest

*     coordonnées de la grille 
      call limit(nindiv,s,xlim,ylim,dt)
      idom = 1
      do ixdom =1,nxdommax
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nydommax
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           float(ixdom-1)*(xlim(2) - xlim(1))/float(nxdommax-1)
            coorddom(2,idom) = ylim(1) +
     &           float(iydom-1)*(ylim(2) - ylim(1))/float(nydommax-1)
            do ipop=1,npopmax
               dom(idom,ipop) = 0.
               domperm(idom,ipop) = 0.
            enddo
            idom = idom + 1
         enddo
      enddo

****************
*     read frequencies for pivot state   
c      write(6,*) 'look for pivot state'
      do iit=1,pivot
         do iloc=1,nloc
c            write(6,*) 'iloc=',iloc
            do ial=1,nalmax
c               write(6,*) 'ial=',ial
               read(13,*) (fpiv(ipop,iloc,ial),ipop=1,npopmax)
            enddo
         enddo
      enddo
      rewind 13
c$$$      iit = 1
c$$$      iitsub = 0
c$$$      do while(iitsub .lt. pivot)
c$$$c         write(6,*) 'iit=',iit
c$$$         read(9,*) npop
c$$$         do iloc=1,nloc
c$$$c     write(6,*) 'iloc=',iloc
c$$$            do ial=1,nalmax
c$$$c     write(6,*) 'ial=',ial
c$$$               read(13,*) (fpiv(ipop,iloc,ial),ipop=1,npopmax)
c$$$            enddo
c$$$         enddo
c$$$         if((npop .eq. 
c$$$            iitsub = iitsub + 1
c$$$c            write(6,*) 'iitsub=',iitsub
c$$$         endif  
c$$$         iit = iit + 1
c$$$      enddo
c$$$c      write(*,*) 'piv = ',iitsub
c$$$c      write(*,*) 'en fortran fpiv=',fpiv
c$$$      rewind 9
c$$$      rewind 13                 

      
**************
*     relabel wrt to pivot or take pivot as estimator
c      write(6,*) 'relabel'
      nnit = 0
      do iit=1,int(float(nit)/float(thinning))
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         do iloc=1,nloc
c             write(6,*) 'iloc=',iloc
            do ial=1,nalmax
c                write(6,*) 'ial=',ial
               read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
            enddo
         enddo  
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
c            write(6,*) 'avant relab order=',order
            nnit = nnit + 1 
            if(npopest .lt. 10) then 
               call Relabel(npopmax,nloc,nalmax,nal,npopest,f,fpiv,
     &              order,ordertmp)
c            write(6,*) 'apre relab order=',order
               write(14,*) (order(ipop),ipop = 1,npopmax)
            endif
c            write(6,*) 'iit=',iit
c            write(6,*) 'pivot=', pivot
            if((npopest .lt. 10) .or. (nnit .eq. pivot)) then 
               call calccell(nxdommax*nydommax,coorddom,
     &              npp,nppmax,u,indcel,distcel)
               do idom=1,nxdommax*nydommax
                  ipop = order(c(indcel(idom)))
c                  write(*,*) 'ipop=',ipop
                  dom(idom,ipop) = dom(idom,ipop) + 1.
               enddo
            endif
         endif
      enddo
      if(npopest .lt. 10) then 
         do idom=1,nxdommax*nydommax
            do ipop=1,npopmax
               dom(idom,ipop) = dom(idom,ipop)/float(nnit)
            enddo
         enddo
      endif

 2000 format (1000(e15.5,1x))
      do idom=1,nxdommax*nydommax
         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
     &        (dom(idom,ipop), ipop=1,npopmax)
      enddo
c      write(*,*) coorddom
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)                 
      end subroutine postprocesschain2


c$$$***********************************************************************
c$$$      subroutine  postprocesschain3(nxdommax,nydommax,burnin,ninrub,
c$$$     &     npopmax,nppmax,nindiv,nloc,nal,nalmax,xlim,ylim,dt,nit,
c$$$     &     thinning,filenpop,filenpp,fileu,filec,filef,fileperm,filedom,
c$$$     &     filemeanf,
c$$$     &     s,u,c,f,pivot,fpiv,dom,coorddom,indcel,distcel,
c$$$     &     order,ordertmp,npopest,meanf)
c$$$      implicit none
c$$$      character*255 fileu,filec,filenpp,filenpop,filedom,filef,fileperm,
c$$$     &     filemeanf
c$$$      integer nit,thinning,npp,npop,iit,nindiv,nxdommax,
c$$$     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indcel,
c$$$     &     ipop,nloc,nal,nalmax,ijunk,order,ordertmp,ipopperm,burnin,
c$$$     &     ninrub,npopest,nnit,iloc,ial,pivot
c$$$      double precision s,u,xlim,ylim,coorddom,dom,domperm,distcel,dt,
c$$$     &     f,fpiv,meanf
c$$$
c$$$      integer iitsub
c$$$*     dimensionnement 
c$$$      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
c$$$     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
c$$$     &     domperm(nxdommax*nydommax,npopmax),
c$$$     &     coorddom(2,nxdommax*nydommax),indcel(nxdommax*nydommax),
c$$$     &     distcel(nxdommax*nydommax),order(npopmax),ordertmp(npopmax),
c$$$     &     f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),nal(nloc),
c$$$     &     meanf(npopmax,nloc,nalmax)
c$$$      open(9,file=filenpop)
c$$$      open(10,file=filenpp)
c$$$      open(11,file=fileu)
c$$$      open(12,file=filec)
c$$$      open(13,file=filef)
c$$$      open(14,file=fileperm)
c$$$      open(15,file=filedom)
c$$$      open(16,file=filemeanf)
c$$$
c$$$c      write(6,*) 'debut postproc order=',order
c$$$
c$$$c      write(6,*) 'npopest=', npopest
c$$$      do ipop = 1,npopmax
c$$$         do iloc = 1,nloc
c$$$            do ial = 1,nalmax
c$$$               meanf(ipop,iloc,ial) = 0
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$*     coordonnées de la grille 
c$$$      call limit(nindiv,s,xlim,ylim,dt)
c$$$      idom = 1
c$$$      do ixdom =1,nxdommax
c$$$c         write(6,*) 'ixdom=',ixdom
c$$$         do iydom=1,nydommax
c$$$c            write(6,*) 'iydom=',iydom
c$$$            coorddom(1,idom) = xlim(1) + 
c$$$     &           float(ixdom-1)*(xlim(2) - xlim(1))/float(nxdommax-1)
c$$$            coorddom(2,idom) = ylim(1) +
c$$$     &           float(iydom-1)*(ylim(2) - ylim(1))/float(nydommax-1)
c$$$            do ipop=1,npopmax
c$$$               dom(idom,ipop) = 0.
c$$$               domperm(idom,ipop) = 0.
c$$$            enddo
c$$$            idom = idom + 1
c$$$         enddo
c$$$      enddo
c$$$
c$$$****************
c$$$*     read frequencies for pivot state   
c$$$c      write(6,*) 'look for pivot state'
c$$$      do iit=1,pivot
c$$$         do iloc=1,nloc
c$$$c            write(6,*) 'iloc=',iloc
c$$$            do ial=1,nalmax
c$$$c               write(6,*) 'ial=',ial
c$$$               read(13,*) (fpiv(ipop,iloc,ial),ipop=1,npopmax)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      rewind 13
c$$$
c$$$
c$$$**************
c$$$*     relabel wrt to pivot or take pivot as estimator
c$$$c      write(6,*) 'relabel'
c$$$      nnit = 0
c$$$      do iit=1,int(float(nit)/float(thinning))
c$$$         read(9,*) npop
c$$$         read(10,*) npp
c$$$         do ipp=1,nppmax
c$$$            read(11,*) u(1,ipp),u(2,ipp)
c$$$            read(12,*) c(ipp)
c$$$         enddo
c$$$         do iloc=1,nloc
c$$$            do ial=1,nalmax
c$$$               read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
c$$$            enddo
c$$$         enddo  
c$$$         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
c$$$c            write(6,*) 'avant relab order=',order
c$$$            nnit = nnit + 1 
c$$$            if(npopest .lt. 10) then 
c$$$               call Relabel(npopmax,nloc,nalmax,nal,npopest,f,fpiv,
c$$$     &              order,ordertmp)
c$$$c            write(6,*) 'apre relab order=',order
c$$$               write(14,*) (order(ipop),ipop = 1,npopmax)
c$$$            endif
c$$$c            write(6,*) 'iit=',iit
c$$$c            write(6,*) 'pivot=', pivot
c$$$            if((npopest .lt. 10) .or. (nnit .eq. pivot)) then 
c$$$               call calccell(nxdommax*nydommax,coorddom,
c$$$     &              npp,nppmax,u,indcel,distcel)
c$$$               do idom=1,nxdommax*nydommax
c$$$                  ipop = order(c(indcel(idom)))
c$$$c                  write(*,*) 'ipop=',ipop
c$$$                  dom(idom,ipop) = dom(idom,ipop) + 1.
c$$$               enddo
c$$$*     increment estimated allele frequencies
c$$$               do ipop = 1,npopmax
c$$$                  do iloc = 1,nloc
c$$$                     do ial = 1,nalmax
c$$$                        meanf(ipop,iloc,ial) =  meanf(ipop,iloc,ial) +  
c$$$     &                       f(order(ipop),iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$      if(npopest .lt. 10) then 
c$$$         do idom=1,nxdommax*nydommax
c$$$            do ipop=1,npopmax
c$$$               dom(idom,ipop) = dom(idom,ipop)/float(nnit)
c$$$            enddo
c$$$         enddo
c$$$         do ipop = 1,npopmax
c$$$            do iloc = 1,nloc
c$$$               do ial = 1,nalmax
c$$$                  meanf(ipop,iloc,ial) = meanf(ipop,iloc,ial) / 
c$$$     &                 float(nnit)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$
c$$$ 2000 format (1000(e15.5,1x))
c$$$      do idom=1,nxdommax*nydommax
c$$$         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
c$$$     &        (dom(idom,ipop), ipop=1,npopmax)
c$$$      enddo
c$$$      
c$$$ 3000 format (300(1x,e15.8,1x))
c$$$      do iloc=1,nloc
c$$$         do ial=1,nalmax
c$$$            write(16,3000) (sngl(meanf(ipop,iloc,ial)),ipop=1,npopmax)
c$$$         enddo
c$$$      enddo  
c$$$      
c$$$c     write(*,*) coorddom
c$$$      close(9)
c$$$      close(10)
c$$$      close(11)
c$$$      close(12)
c$$$      close(13)
c$$$      close(14)
c$$$      close(15)
c$$$      close(16)
c$$$      end subroutine postprocesschain3


***********************************************************************
      subroutine  postprocessmultchain(nxdommax,nydommax,burnin,ninrub,
     &     npopmax,nppmax,nindiv,nloc,nal,nalmax,xlim,ylim,dt,nit,
     &     thinning,nrun,pathall,nchpathall,
     &     s,u,c,f,pivot,chirunpiv,nchirunpiv,fpiv,dom,coorddom,indcel,
     &     distcel,order,ordertmp,npopest)
      implicit none
      character*255 fileu,filec,filenpp,filenpop,filedom,filef,fileperm,
     &     pathall,path,chirunpiv,chirun,fileftmp     
      integer nit,thinning,npp,npop,iit,nindiv,nxdommax,
     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indcel,
     &     ipop,nloc,nal,nalmax,ijunk,order,ordertmp,ipopperm,burnin,
     &     ninrub,npopest,nnit,iloc,ial,pivot,nrun,irun,nchpathall,
     &     irunpiv,nchirunpiv,resirun,nchpath
      double precision s,u,xlim,ylim,coorddom,dom,domperm,distcel,dt,
     &     f,fpiv
      integer iitsub
*     dimensionnement 
      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
     &     domperm(nxdommax*nydommax,npopmax),
     &     coorddom(2,nxdommax*nydommax),indcel(nxdommax*nydommax),
     &     distcel(nxdommax*nydommax),order(npopmax),ordertmp(npopmax),
     &     f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),nal(nloc)

c      write(6,*) 'debut postproc order=',order
      fileperm = pathall(1:nchpathall) // "/perm.txt"
      filedom = pathall(1:nchpathall) // "/proba.pop.membership.txt"
      open(14,file=fileperm)
      open(15,file=filedom)
******************************
*     coordonnées de la grille 
      call limit(nindiv,s,xlim,ylim,dt)
      idom = 1
      do ixdom =1,nxdommax
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nydommax
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           float(ixdom-1)*(xlim(2) - xlim(1))/float(nxdommax-1)
            coorddom(2,idom) = ylim(1) +
     &           float(iydom-1)*(ylim(2) - ylim(1))/float(nydommax-1)
            do ipop=1,npopmax
               dom(idom,ipop) = 0.
               domperm(idom,ipop) = 0.
            enddo
            idom = idom + 1
         enddo
      enddo

****************************************
*     read frequencies for pivot state  
*     index of the run containing pivot state
c      write(*,*) 'en fortran pathall=',pathall
c      write(*,*) 'en fortran pivot=',pivot
c      write(*,*) 'en fortran nchpathall=',nchpathall
c      write(*,*) 'en fortran chirunpiv=',chirunpiv
c      write(*,*) 'en fortran nchirunpiv=',nchirunpiv
c      nchpathall = len_trim(pathall)
c      nchirunpiv =  len_trim(chirunpiv)
      irunpiv = 1 + 
     &     int(aint(float(pivot-1)/(float(nit)/float(thinning))))
      fileftmp = pathall(1:nchpathall) // chirunpiv(1:nchirunpiv)
c      write(*,*) 'en fortran fileftmp=',fileftmp
      filef = fileftmp(1:(nchpathall+nchirunpiv)) // "/frequencies.txt"
c      write(*,*) 'en fortran filef=',filef
      open(13,file=filef)
*     index of saved iteration containing pivot  in this file
      iit = pivot - (irunpiv-1)*float(nit)/float(thinning)
      do iitsub = 1,iit
         do iloc=1,nloc
c     write(6,*) 'iloc=',iloc
            do ial=1,nalmax
c     write(6,*) 'ial=',ial
               read(13,*) (fpiv(ipop,iloc,ial),ipop=1,npopmax)
            enddo
         enddo
      enddo
      close(13)


      nnit = 0
      do iit=1,int(float(nrun*nit)/float(thinning))
         irun = 1 + aint(float(iit-1)/(float(nit)/float(thinning)))
         if(mod(iit-1,int(float(nit)/float(thinning))) .eq. 0) then 
c            write(*,*) 'iit=',iit
c            write(*,*) 'irun=',irun
*     define path to file and  open files (Ref.: book Maryse Ain, p.340) 
c           write(*,*) 'opening files'
c            write(*,*) 'pathall=',pathall
            resirun = irun
            if(irun .gt. 999) then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/1000)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 1000*int(aint(float(irun)/1000))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/100)) + ichar('0'))
               nchpath = nchpath + 1
               resirun = resirun - 100*int(aint(float(resirun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 99) .and. (irun .le. 999))then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/100)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 100*int(aint(float(irun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               resirun = resirun - 10*int(aint(float(resirun)/10))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 9) .and. (irun .le. 99)) then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/10)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 10*int(aint(float(irun)/10))
               path = path(1:nchpath) // 
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if(irun .le. 9) then 
               path = pathall(1:nchpathall) // char(irun + ichar('0'))
               nchpath = nchpathall + 1
            endif

            filenpp = path(1:nchpath) // "/nuclei.numbers.txt"
            filenpop = path(1:nchpath)// "/populations.numbers.txt"
            fileu = path(1:nchpath)   // "/coord.nuclei.txt"
            filec = path(1:nchpath)   // "/color.nuclei.txt"
            filef = path(1:nchpath)   // "/frequencies.txt"
            open(9,file=filenpop)
            open(10,file=filenpp)
            open(11,file=fileu)
            open(12,file=filec)
            open(13,file=filef)
         endif


*     read and relabel  wrt to pivot or take pivot as estimator
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         do iloc=1,nloc
            do ial=1,nalmax
               read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
            enddo
         enddo  
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then 
               call Relabel(npopmax,nloc,nalmax,nal,npopest,f,fpiv,
     &              order,ordertmp)
               write(14,*) (order(ipop),ipop = 1,npopmax)
            endif
            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
               call calccell(nxdommax*nydommax,coorddom,
     &              npp,nppmax,u,indcel,distcel)
               do idom=1,nxdommax*nydommax
                  ipop = order(c(indcel(idom)))
                  dom(idom,ipop) = dom(idom,ipop) + 1.
               enddo
            endif
         endif
      
         if(mod(iit,int(float(nit)/float(thinning))) .eq. 0) then 
            close(9)
            close(10)
            close(11)
            close(12)
            close(13)
         endif
      enddo
      if(npopest .lt. 10) then 
         do idom=1,nxdommax*nydommax
            do ipop=1,npopmax
               dom(idom,ipop) = dom(idom,ipop)/float(nnit)
            enddo
         enddo
      endif

 2000 format (1000(e15.5,1x))
      do idom=1,nxdommax*nydommax
         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
     &        (dom(idom,ipop), ipop=1,npopmax)
      enddo  
      close(14)
      close(15)     
      
      end subroutine postprocessmultchain




*********************************************************************
*     posterior probability of population membership for individuals
*
      subroutine  pppmindiv2(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,filenpop,filenpp,fileu,filec,
     &     fileperm,nit,thinning,burnin,order,npopest,pivot)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npop,npp,c,
     &     nit,burnin,npopest,order(npopmax),thinning,pivot
      double precision pmp,distcell,u,s

      integer iit,ipp,iindiv,ipop,nnit,iitsub
      character*255 filenpop,fileu,filec,filenpp,fileperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
c$$$      write(6,*) '      **********************************************'
c$$$      write(6,*) '      *  Computing posterior probabilities          '
c$$$      write(6,*) '      *  of population membership for individuals   '
c$$$      write(6,*) '      **********************************************'
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'
      open(9,file=filenpop)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      open(13,file=fileperm)
      nnit = 0 
      do iit=1,int(float(nit)/float(thinning))
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then
               read(13,*) (order(ipop),ipop=1,npopmax)
            endif
            if((npopest .lt. 10) .or. (nnit .eq. pivot)) then 
               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
               do iindiv=1,nindiv
                  ipop = order(c(indcell(iindiv)))
                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
               enddo
            endif
         endif
      enddo
      if(npopest .lt. 10) then
         do iindiv=1,nindiv 
            do ipop=1,npopmax
               pmp(iindiv,ipop) = pmp(iindiv,ipop)/float(nnit)
            enddo
         enddo
      endif
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)

c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      end subroutine  pppmindiv2



*********************************************************************
*     posterior probability of population membership for individuals
*
      subroutine  pppmindivmultchain(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,pathall,nchpathall,nrun,
     &     nit,thinning,burnin,order,npopest,pivot)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npop,npp,c,nrun,
     &     nit,burnin,npopest,order(npopmax),thinning,pivot,
     &     nchpathall
      double precision pmp,distcell,u,s
      character*255 pathall

      integer iit,ipp,iindiv,ipop,nnit,iitsub,irun,nchpath,resirun
      character*255 path,filenpop,fileu,filec,filenpp,fileperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
c$$$      write(6,*) '      **********************************************'
c$$$      write(6,*) '      *  Computing posterior probabilities          '
c$$$      write(6,*) '      *  of population membership for individuals   '
c$$$      write(6,*) '      **********************************************'
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      fileperm = pathall(1:nchpathall) // "/perm.txt"
      open(13,file=fileperm)

      nnit = 0
      do iit=1,int(float(nrun*nit)/float(thinning))
         irun = 1 + aint(float(iit-1)/(float(nit)/float(thinning)))
         if(mod(iit-1,int(float(nit)/float(thinning))) .eq. 0) then 
c            write(*,*) 'iit=',iit
c            write(*,*) 'irun=',irun
*     open files
c           write(*,*) 'opening files'
c            write(*,*) 'pathall=',pathall
*     cf book M Ain, p. 340 
            resirun = irun
            if(irun .gt. 999)then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/1000)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 1000*int(aint(float(irun)/1000))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/100)) + ichar('0'))
               nchpath = nchpath + 1
               resirun = resirun - 100*int(aint(float(resirun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 99) .and. (irun .le. 999))then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/100)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 100*int(aint(float(irun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               resirun = resirun - 10*int(aint(float(resirun)/10))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 9) .and. (irun .le. 99)) then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/10)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 10*int(aint(float(irun)/10))
               path = path(1:nchpath) // 
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if(irun .le. 9) then 
               path = pathall(1:nchpathall) // char(irun + ichar('0'))
               nchpath = nchpathall + 1
            endif
            filenpp = path(1:nchpath) // "/nuclei.numbers.txt"
            filenpop = path(1:nchpath) // "/populations.numbers.txt"
            fileu = path(1:nchpath) // "/coord.nuclei.txt"
            filec = path(1:nchpath) // "/color.nuclei.txt"           
            open(9,file=filenpop)
            open(10,file=filenpp)
            open(11,file=fileu)
            open(12,file=filec)
         endif
*     do the real job now
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then
               read(13,*) (order(ipop),ipop=1,npopmax)
            endif
            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
               do iindiv=1,nindiv
                  ipop = order(c(indcell(iindiv)))
                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
               enddo
            endif
         endif
*     close files
         if(mod(iit,int(float(nit)/float(thinning))) .eq. 0) then 
            close(9)
            close(10)
            close(11)
            close(12)
         endif
      enddo
      close(13)
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$      open(9,file=filenpop)
c$$$      open(10,file=filenpp)
c$$$      open(11,file=fileu)
c$$$      open(12,file=filec)
c$$$      open(13,file=fileperm)
c$$$      nnit = 0 
c$$$      do iit=1,int(float(nit)/float(thinning))
c$$$         read(9,*) npop
c$$$         read(10,*) npp
c$$$         do ipp=1,nppmax
c$$$            read(11,*) u(1,ipp),u(2,ipp)
c$$$            read(12,*) c(ipp)
c$$$         enddo
c$$$         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
c$$$            nnit = nnit + 1 
c$$$            if(npopest .lt. 10) then
c$$$               read(13,*) (order(ipop),ipop=1,npopmax)
c$$$            endif
c$$$            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
c$$$               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
c$$$               do iindiv=1,nindiv
c$$$                  ipop = order(c(indcell(iindiv)))
c$$$                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
c$$$               enddo
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$      if(npopest .lt. 10) then
c$$$         do iindiv=1,nindiv 
c$$$            do ipop=1,npopmax
c$$$               pmp(iindiv,ipop) = pmp(iindiv,ipop)/float(nnit)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$      close(9)
c$$$      close(10)
c$$$      close(11)
c$$$      close(12)
c$$$      close(13)
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCc

c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'
      end subroutine  pppmindivmultchain


************************************************
      Subroutine Relabel(npopmax,nloc,nalmax,nal,npop,f,fpiv,order,
     &     ordertmp)
*     find partition that minimizes scalar product between f and fpiv
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order,ordertmp
      double precision f,fpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),ordertmp(npopmax),nal(nloc)
      double precision sp,sptmp,spf
      integer ipop
      Integer*4 I,I1,J,G,H

      sp = 0 
      do ipop=1,npop
         ordertmp(ipop) = ipop
      enddo
      If (npop.Gt.1) Go To 10
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
90    Return
10    Continue
      I=npop-2
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      G=ordertmp(npop-1)
      H=ordertmp(npop)
      If (G .eq. H) Go To 20
      ordertmp(npop)=G
      ordertmp(npop-1)=H
C      Call Sum(npop,order)
c      Write(*,*) 'order=',order
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      ordertmp(npop-1)=G
      ordertmp(npop)=H
20    Continue
      If (I.Eq.0) Go To 90
      H=ordertmp(I)
      I1=I+1
      Do 30 J=I1,npop
      If (ordertmp(J) .Le. H) Go To 30
      ordertmp(I)=ordertmp(J)
      ordertmp(J)=H
      Go To 10
30    Continue
31    Continue
      Do 40 J=I1,npop
      ordertmp(J-1)=ordertmp(J)
40    Continue
      ordertmp(npop)=H
      I=(I-1)
      Go To 20
      End Subroutine Relabel

************************************************************************
*     scalar product of two arrays of frequencies
      double precision function spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,
     &     order)
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order
      double precision f,fpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),nal(nloc)
      integer ipop,iloc,ial
      spf = 0
      do ipop=1,npop
         do iloc = 1,nloc
            do ial=1,nal(iloc)
               spf = spf + f(order(ipop),iloc,ial)*fpiv(ipop,iloc,ial)
            enddo
         enddo
      enddo
      end function spf

      
c$$$ 
c$$$
c$$$************************************************
c$$$      SUBROUTINE PERMUT (N,E)
c$$$C<title> CALCULATES ALL PERMUTATIONS OF AN ARRAY (E1,.....EN) </title>
c$$$C=====IN LEXICOGRAPHIC ORDER WITHOUT REPETITION.
c$$$      INTEGER N
c$$$      INTEGER E
c$$$      DIMENSION E(N)
c$$$C
c$$$C  ARGUMENTS:
c$$$C   N:NUMBER OF ELEMENTS TO PERMUTE
c$$$C   E:COMPONENTS OF VECTOR E ARE THE NUMBERS TO BE PERMUTED,
c$$$C     THEY MUST BE ORDERED SO,THAT E(I-1) <= E(I),
c$$$C     THE ORIGINAL ORDER WILL BE RESTORED.
c$$$C  SUM IS A ROUTINE TO BE CALLED BY CALL SUM(N,E) AFTER EACH
c$$$C      PERMUTATION TO ACT ON IT.
c$$$C
c$$$      INTEGER*4 G,H
c$$$      IF (N.GT.1) GO TO 10
c$$$c      CALL SUM(N,E)
c$$$      write(*,*) 'E=',E
c$$$90    RETURN
c$$$10    CONTINUE
c$$$      I=N-2
c$$$c      CALL SUM(N,E)
c$$$      write(*,*) 'E=',E
c$$$      G=E(N-1)
c$$$      H=E(N)
c$$$      IF (G.EQ.H) GO TO 20
c$$$      E(N)=G
c$$$      E(N-1)=H
c$$$c      CALL SUM(N,E)
c$$$      write(*,*) 'E=',E
c$$$      E(N-1)=G
c$$$      E(N)=H
c$$$20    CONTINUE
c$$$      IF (I.EQ.0) GO TO 90
c$$$      H=E(I)
c$$$      I1=I+1
c$$$      DO 30 J=I1,N
c$$$      IF (E(J) .LE. H) GO TO 30
c$$$      E(I)=E(J)
c$$$      E(J)=H
c$$$      GO TO 10
c$$$30    CONTINUE
c$$$31    CONTINUE
c$$$      DO 40 J=I1,N
c$$$      E(J-1)=E(J)
c$$$40    CONTINUE
c$$$      E(N)=H
c$$$      I=(I-1)
c$$$      GO TO 20
c$$$      END
c$$$      
c$$$ 
c$$$C      
c$$$
c$$$
c$$$C     
