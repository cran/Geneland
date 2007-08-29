***********************************************************************  
*     Simulation selon le posterior des parametres 
*     du modele de melange de genotypes par une Voronoi cachee
*
*     Routines de simulations  
*     de ~/mathlib/randlib-1.3.f/src/all.f
*
***********************************************************************
     
   
      subroutine mcmc(s,z,nall,ploidy,filelambda,
     &     filenpp,fileu,filec,filef,filefa,filedrift,
     &     filenpop,filet,filelpp,filell,filesize,lambdamax,
     &     dt,nit,thinning,
     &     nindiv,nloc,nloc2,nallmax,npp,nppmax,
     &     npop,npopmin,npopmax,
     &     t,ttmp,u,utmp,c,ctmp,f,ftmp,fa,drift,drifttmp,
     &     indcell,indcelltmp,
     &     distcell,distcelltmp,xlim,ylim,n,ntmp,a,ptmp,
     &     cellpop,listcell,fmodel,kfix,spatial,jcf,seed1,seed2) 
      implicit none 

*     les donnees
      integer nindiv,nloc,nloc2,
     &     nall,nallmax,z,ploidy,jcf,seed1,seed2
      real s

*     parametres des priors 
      real lambdamax,dt

*     les parametres a inferer (et valeurs maximales)
      integer npp,nppmax,npop,npopmin,npopmax,
     &     c,ctmp
      real lambda,u,utmp,f,t,fa,drift,ftmp,drifttmp

*     variables de travail
      integer iloc,iindiv,ichain,nit,
     &     ipp,ipop,iall,thinning,indcell,indcelltmp,
     &     n,cellpop,listcell,
     &     cellpophost,ntmp,fmodel,kfix,spatial
      real ptmp,xlim,ylim,ranf,rpostlamb,
     &     distcell,du,distcelltmp,a,ttmp,lpp,ll
      character*255 filef,filenpp,
     &     filelambda,filenpop,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet,filesize

      integer nn
 
*     dimensionnement 
      dimension s(2,nindiv),t(2,nindiv),z(nindiv,nloc2),
     &     u(2,nppmax),utmp(2,nppmax),c(nppmax),ctmp(nppmax),
     &     f(npopmax,nloc,nallmax),xlim(2),ylim(2),
     &     nall(nloc),indcell(nindiv),indcelltmp(nindiv),
     &     distcell(nindiv),ttmp(2,nindiv),
     &     distcelltmp(nindiv),n(npopmax,nloc,nallmax),
     &     ntmp(npopmax,nloc,nallmax),
     &     a(nallmax),ptmp(nallmax),
     &     fa(nloc,nallmax),drift(npopmax),
     &     ftmp(npopmax,nloc,nallmax),drifttmp(npopmax),
     &     cellpop(nppmax),listcell(nppmax),cellpophost(nppmax)



      write(6,*) '          *****************************'
      write(6,*) '          ***    MCMC inference     ***'
      write(6,*) '          *****************************'

 


      call setall(seed1,seed2) 
   
      nloc = nloc
  
 
*     look for smallest rectangle enclosing the spatial domain
      call limit(nindiv,s,xlim,ylim,dt)
c      write(*,*) 'fin de limit'
c      write(*,*) 's=',s
c      write(*,*) 'z=',z
c      write(*,*) 'nall=',nall


*     Ouverture des fichiers pour l'ecriture des sorties
      open(9,file=filelambda)
      open(10,file=filenpp)
      open(11,file=filenpop)
      open(12,file=fileu)
      open(13,file=filec)
      open(14,file=filef)
      open(15,file=filefa) 
      open(16,file=filedrift)
      open(17,file=filelpp)
      open(18,file=filell) 
      open(19,file=filet)
      open(20,file=filesize)

c       write(*,*) 'fin de l ouverture'



 
************************
*     Initialization
************************
      du = sqrt((xlim(2)-xlim(1))*(ylim(2)-ylim(1))/nindiv)
      lambda = lambdamax*ranf()
      if(spatial .eq. 1) then
         npp = 1 + int(aint(lambda))
      else
         npp = nindiv
      endif


      do iindiv=1,nindiv
         t(1,iindiv) = s(1,iindiv) + dt*(ranf()-.5)
         t(2,iindiv) = s(2,iindiv) + dt*(ranf()-.5)
      enddo
         
      if(spatial .eq. 1) then 
         call rprioru(npp,nppmax,xlim,ylim,u)
      else 
         do iindiv=1,nindiv
            u(1,iindiv) = s(1,iindiv)
            u(2,iindiv) = s(2,iindiv) 
         enddo
      endif  

      call rpriorc(npp,nppmax,npop,c)
      call calccell(nindiv,t,npp,nppmax,u,indcell,distcell)
      call rpriordrift(npop,npopmax,drift,fmodel)
      call rpriorfa(nloc,nloc,nall,nallmax,fa,fmodel,ptmp)
      call rpriorf(npop,npopmax,nloc,nloc,nall,nallmax,f,ptmp)
 

c$$$      write(*,*) 'nindiv =', nindiv,"\n"
c$$$      write(*,*) 'nloc =',nloc   ,"\n"
c$$$      write(*,*) 'nall =',nall ,"\n"
c$$$      write(*,*) 'nallmax =', nallmax,"\n"
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



************************  
* mise a jour iterative  
************************
 
c      write(*,*) 'Starting updates'
      do ichain=1,nit
*     ecriture dans les fichiers (tous les thinning)
         if(mod(ichain,thinning) .eq. 0) then
 100        format(f7.3,' %')
            write(6,100)float(ichain)/float(nit)*100.
            write(9,*) lambda
            write(10,*) npp
            write(11,*) npop
 1000       format (2(1x,e15.8,1x))
            
            write(12,1000) (u(1,ipp),u(2,ipp), ipp=1,nppmax)
            do ipp=1,nppmax
               write(13,*) c(ipp)
            enddo
            
 2000       format (300(1x,e15.8,1x))
            do iloc=1,nloc
               do iall=1,nallmax
                  write(14,2000) (f(ipop,iloc,iall),
     &                 ipop=1,npopmax)
               enddo
            enddo  
 3000       format (300(1x,e15.8,1x))    
            write(15,3000) ((fa(iloc,iall),iall=1,nallmax),iloc=1,nloc)
            write(16,2000) (drift(ipop),ipop=1,npopmax)
            write(17,*) lpp(lambdamax,lambda,z,npop,npp,drift,f,fa,c,
     &           nppmax,nindiv,nloc2,npopmax,nloc,nallmax,
     &           indcell,nall,fmodel,xlim,ylim)           
            write(18,*) ll(z,nindiv,nloc,nloc2,npopmax,
     &           nallmax,nppmax,c,f,indcell)
            if(dt .gt. 1.e-30) then 
               write(19,1000) (t(1,iindiv),t(2,iindiv),iindiv=1,nindiv)
            endif

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




*     update lambda
*          write(*,*) 'update lambda'
          if(spatial .eq. 1) then
             lambda = rpostlamb(lambdamax,npp)
          endif

          if(fmodel .eq. 1) then 
*     update drift
*              write(*,*) 'update drift'
             call  upddrift(npop,npopmax,nloc,nallmax,nall,
     &            f,fa,drift)
*     update fa
*             write(*,*) 'update fa'
             call updfa(npop,npopmax,nloc,nallmax,nall,
     &            f,fa,drift) 
          endif
 
*     update f 
*          write(*,*) 'update f'  
          if(jcf .eq. 1) then 
c     write(*,*) 'update c and f' 
             if(fmodel .eq. 0) then 
*     joint update of c anf f
                call  udcf(npop,npopmax,f,fa,drift,
     &               nloc,nloc,nloc2,
     &               nall,nallmax,indcell,nindiv,npp,nppmax,
     &               c,ctmp,a,ptmp,ftmp,z,n,ntmp,ploidy)
             else
                call udcf2(npop,npopmax,f,fa,drift,
     &               nloc,nloc,nloc2,
     &               nall,nallmax,indcell,nindiv,npp,nppmax,
     &               c,ctmp,a,ptmp,ftmp,z,n,ntmp,ploidy)
             endif
          else          
             call rpostf2(npop,npopmax,nloc,nloc,nall,nallmax,
     &            f,fa,drift,
     &            nindiv,nloc2,z,nppmax,c,indcell,
     &            n,a,ptmp,ploidy)
             call updc(npp,nppmax,c,ctmp,z,nindiv,nloc,
     &         nloc,nloc2,nallmax,npop,npopmax,f,indcell,ploidy)

          endif
 

          if(spatial .eq. 1) then 
*     update u et mise a jour de indcell et distcell
c$$$             write(*,*) 'before update u'
c$$$             write(*,*) 'xlim=', xlim 
             call updurw(npp,nppmax,c,u,z,nindiv,nloc,nloc,
     &            nloc2,nallmax,npopmax,f,indcell,distcell,
     &            indcelltmp,distcelltmp,t,xlim,ylim,du,ploidy)
c$$$             write(*,*) 'xlim=', xlim
c$$$             write(*,*) 'after update u'
             
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.e-30) then 
*                write(*,*) 'update t'
                call updt(npp,nppmax,nindiv,
     &               nloc,nloc,nloc2,nallmax,npopmax,
     &               t,ttmp,dt,s,c,indcell,distcell,
     &               indcelltmp,distcelltmp,u,z,f,ploidy)
             endif

*     birth/death des points du pp
*             write(*,*) 'update npp'
             call bdpp(nindiv,u,c,utmp,ctmp,
     &            npop,npopmax,nloc,nloc,
     &            nloc2,nallmax,npp,nppmax,z,f,t,xlim,ylim,indcell,
     &            distcell,indcelltmp,distcelltmp,lambda,ploidy)
          endif


*     birth/death de pop
         if(kfix .eq. 0) then 
*            write(*,*) 'update npop'
* split/merge avec prop de f selon cond. complete dans les deux sens
            if(fmodel .eq. 0) then 
*     avec drift=0.5 et fa=1 pour court-circuiter le F-model
               call bdpop9bis(npop,npopmin,npopmax,f,fa,drift,
     &              nloc,nloc,nloc2,
     &              nall,nallmax,indcell,nindiv,npp,nppmax,c,
     &              ctmp,a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &              cellpophost,n,ntmp,ploidy)
            else
               call bdpop9(npop,npopmin,npopmax,f,fa,drift,
     &              nloc,nloc,nloc2,
     &              nall,nallmax,indcell,nindiv,npp,nppmax,c,
     &              ctmp,a,ptmp,ftmp,drifttmp,z,cellpop,listcell,
     &              cellpophost,n,ntmp,ploidy)
            endif
         endif

    

      enddo
       
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)

      close(20)
       
      write(6,*) '          ************************************'
      write(6,*) '          ***    End of MCMC inference     ***'
      write(6,*) '          ************************************'

c$$$      write(*,*) 'nindiv =', nindiv,"\n"
c$$$      write(*,*) 'nloc =',nloc   ,"\n"
c$$$      write(*,*) 'nloc2 =', nloc2,"\n"
c$$$      write(*,*) 'nall =',nall ,"\n"
c$$$      write(*,*) 'nallmax =', nallmax,"\n"
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



      end
      
c      include 'sub.f'
      

 
