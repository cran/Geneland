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
     &     dt,nchain,stepw,
     &     nindiv,nlocmax,nlocmax2,nallmax,npp,nppmax,
     &     npop,npopmin,npopmax,
     &     t,ttemp,u,utemp,c,ctemp,f,ftemp,fa,drift,drifttemp,
     &     indcell,indcelltemp,
     &     distcell,distcelltemp,n,ntemp,a,ptemp,
     &     cellpop,listcell,fmodel,kfix,spatial,jcf,seed1,seed2) 
      implicit none 

*     les donnees
      integer nindiv,nloc,nlocmax,nlocmax2,
     &     nall,nallmax,z,ploidy,jcf,seed1,seed2
      real s

*     parametres des priors 
      real lambdamax,dt

*     les parametres a inferer (et valeurs maximales)
      integer npp,nppmax,npop,npopmin,npopmax,
     &     c,ctemp
      real lambda,u,utemp,f,t,fa,drift,ftemp,drifttemp

*     variables de travail
      integer iloc,iindiv,ichain,nchain,
     &     ipp,ipop,iall,stepw,indcell,indcelltemp,
     &     n,cellpop,listcell,
     &     cellpophost,ntemp,fmodel,kfix,spatial
      real ptemp,xlim(2),ylim(2),ranf,rpostlamb,
     &     distcell,du,distcelltemp,a,ttemp,lpp,ll
      character*256 filef,filenpp,
     &     filelambda,filenpop,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet,filesize

      integer nn
 
*     dimensionnement 
      dimension s(2,nindiv),t(2,nindiv),z(nindiv,nlocmax2),
     &     u(2,nppmax),utemp(2,nppmax),c(nppmax),ctemp(nppmax),
     &     f(npopmax,nlocmax,nallmax),
     &     nall(nlocmax),indcell(nindiv),indcelltemp(nindiv),
     &     distcell(nindiv),ttemp(2,nindiv),
     &     distcelltemp(nindiv),n(npopmax,nlocmax,nallmax),
     &     ntemp(npopmax,nlocmax,nallmax),
     &     a(nallmax),ptemp(nallmax),
     &     fa(nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     cellpop(nppmax),listcell(nppmax),cellpophost(nppmax)



      write(6,*) '          *****************************'
      write(6,*) '          ***    MCMC inference     ***'
      write(6,*) '          *****************************'

 


      call setall(seed1,seed2) 
   
      nloc = nlocmax
  
 
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     for debugging only
c      npp = 2
CCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     for debugging only
c$$$      c(1) = 1 
c$$$      c(2) = 1
c$$$      u(1,1) = 0.25
c$$$      u(2,1) = 0.5  
c$$$      u(1,2) = 0.75
c$$$      u(2,2) = 0.5
c      write(*,*) 'f=',f
c      write(*,*) 'ftemp=',ftemp
c      write(*,*) 'c=', c
c      write(*,*) 'u=',(u(1,ipp),ipp=1,nppmax)
c      write(*,*) 'u=',(u(2,ipp),ipp=1,nppmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCC


      call calccell(nindiv,t,npp,nppmax,u,indcell,distcell)
c      write(*,*) 'indcell =',indcell
      
      call rpriordrift(npop,npopmax,drift,fmodel)

      call rpriorfa(nloc,nlocmax,nall,nallmax,fa,fmodel,ptemp)
       
      call rpriorf(npop,npopmax,nloc,nlocmax,nall,nallmax,f,ptemp)
 
C
  
************************  
* mise a jour iterative  
************************
 
c      write(*,*) 'Starting updates'
      do ichain=1,nchain
*     ecriture dans les fichiers (tous les stepw)
          if(mod(ichain,stepw) .eq. 0) then
 100         format(f7.3,' %')
             write(6,100)float(ichain)/float(nchain)*100.
             write(9,*) lambda
             write(10,*) npp
             write(11,*) npop
 1000        format (2(1x,e15.8,1x))
             write(12,1000) (u(1,ipp),u(2,ipp), ipp=1,nppmax)
             do ipp=1,nppmax
                write(13,*) c(ipp)
             enddo
c$$$ 2000        format (2(e10.5,1x))
c$$$             write(14,2000) (((f(ipop,iloc,iall),ipop=1,npopmax),
c$$$     &            iall=1,nallmax),iloc=1,nloc)
 2000        format (300(1x,e15.8,1x))
             do iloc=1,nlocmax
                do iall=1,nallmax
                   write(14,2000) (f(ipop,iloc,iall),
     &                  ipop=1,npopmax)
                enddo
             enddo  
 3000        format (300(1x,e15.8,1x))    
             write(15,3000) ((fa(iloc,iall),iall=1,nallmax),iloc=1,nloc)
             write(16,2000) (drift(ipop),ipop=1,npopmax)
             write(17,*) lpp(lambda,z,npop,npp,drift,f,fa,c,nppmax,
     &            nindiv,nlocmax2,npopmax,nlocmax,nallmax,indcell,
     &            nall,fmodel)           
             write(18,*) ll(z,nindiv,nlocmax,nlocmax2,npopmax,
     &     nallmax,nppmax,c,f,indcell)
             if(dt .gt. 1.e-30) then 
                write(19,1000) (t(1,iindiv),t(2,iindiv),iindiv=1,nindiv)
             endif

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
             call  upddrift(npop,npopmax,nlocmax,nallmax,nall,
     &            f,fa,drift)
*     update fa
*             write(*,*) 'update fa'
             call updfa(npop,npopmax,nlocmax,nallmax,nall,
     &            f,fa,drift) 
          endif
 
*     update f 
*          write(*,*) 'update f'  
          if(jcf .eq. 1) then 
c     write(*,*) 'update c and f' 
             if(fmodel .eq. 0) then 
*     joint update of c anf f
                call  udcf(npop,npopmax,f,fa,drift,
     &               nloc,nlocmax,nlocmax2,
     &               nall,nallmax,indcell,nindiv,npp,nppmax,
     &               c,ctemp,a,ptemp,ftemp,z,n,ntemp,ploidy)
             else
                call udcf2(npop,npopmax,f,fa,drift,
     &               nloc,nlocmax,nlocmax2,
     &               nall,nallmax,indcell,nindiv,npp,nppmax,
     &               c,ctemp,a,ptemp,ftemp,z,n,ntemp,ploidy)
             endif
          else          
             call rpostf2(npop,npopmax,nloc,nlocmax,nall,nallmax,
     &            f,fa,drift,
     &            nindiv,nlocmax2,z,nppmax,c,indcell,
     &            n,a,ptemp,ploidy)
             call updc(npp,nppmax,c,ctemp,z,nindiv,nloc,
     &         nlocmax,nlocmax2,nallmax,npop,npopmax,f,indcell,ploidy)

          endif
 

          if(spatial .eq. 1) then 
*     update u et mise a jour de indcell et distcell
*             write(*,*) 'update u'
             call updurw(npp,nppmax,c,u,z,nindiv,nloc,nlocmax,
     &            nlocmax2,nallmax,npopmax,f,indcell,distcell,
     &            indcelltemp,distcelltemp,t,xlim,ylim,du,ploidy)
             
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.e-30) then 
*                write(*,*) 'update t'
                call updt(npp,nppmax,nindiv,
     &               nloc,nlocmax,nlocmax2,nallmax,npopmax,
     &               t,ttemp,dt,s,c,indcell,distcell,
     &               indcelltemp,distcelltemp,u,z,f,ploidy)
             endif

*     birth/death des points du pp
*             write(*,*) 'update npp'
             call bdpp(nindiv,u,c,utemp,ctemp,
     &            npop,npopmax,nloc,nlocmax,
     &            nlocmax2,nallmax,npp,nppmax,z,f,t,xlim,ylim,indcell,
     &            distcell,indcelltemp,distcelltemp,lambda,ploidy)
          endif


*     birth/death de pop
         if(kfix .eq. 0) then 
*            write(*,*) 'update npop'
* split/merge avec prop de f selon cond. complete dans les deux sens
            if(fmodel .eq. 0) then 
*     avec drift=0.5 et fa=1 pour court-circuiter le F-model
               call bdpop9bis(npop,npopmin,npopmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &              cellpophost,n,ntemp,ploidy)
            else
               call bdpop9(npop,npopmin,npopmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &              cellpophost,n,ntemp,ploidy)
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
      end
      
c      include 'sub.f'
      

 
