***********************************************************************  
*     Simulation selon le posterior des parametres 
*     du modele de melange de genotypes par une Voronoi cachee
*
*     Routines de simulations 
*     de ~/mathlib/randlib-1.3.f/src/all.f
*
*
***********************************************************************


      subroutine mcmc(s,z,nall,filelambda,
     &     filenpp,fileu,filec,filef,filefa,filedrift,
     &     filenclass,filet,filelpp,filell,lambdamax,
     &     dt,nchain,stepw,
     &     nindivmax,nlocmax,nlocmax2,nallmax,npp,nppmax,
     &     nclass,nclassmin,nclassmax,
     &     t,ttemp,u,utemp,c,ctemp,f,ftemp,fa,drift,drifttemp,
     &     indcell,indcelltemp,
     &     distcell,distcelltemp,n,ntemp,a,ptemp,effcl,iclv,
     &     cellclass,listcell,fmodel,kfix,spatial) 
      implicit none 

*     les donnees
      integer nindiv,nindivmax,nloc,nlocmax,nlocmax2,
     &     nall,nallmax,z
      real s

*     parametres des priors 
      real lambdamax,dt

*     les parametres a inferer (et valeurs maximales)
      integer npp,nppmax,nclass,nclassmin,nclassmax,
     &     c,ctemp
      real lambda,u,utemp,f,t,fa,drift,ftemp,drifttemp

*     variables de travail
      integer iloc,iindiv,ichain,nchain,
     &     ipp,iclass,iall,stepw,indcell,indcelltemp,
     &     n,effcl,iclv,cellclass,listcell,
     &     cellclasshost,ntemp,fmodel,kfix,spatial
      real ptemp,xlim(2),ylim(2),ranf,rpostlamb,
     &     distcell,du,distcelltemp,a,ttemp,lpp,ll
      character*200 files,filez,filef,filenall,filenpp,
     &     filelambda,filenclass,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet
 
*     dimensionnement 
      dimension s(2,nindivmax),t(2,nindivmax),z(nindivmax,nlocmax2),
     &     u(2,nppmax),utemp(2,nppmax),c(nppmax),ctemp(nppmax),
     &     f(nclassmax,nlocmax,nallmax),
     &     nall(nlocmax),indcell(nindivmax),indcelltemp(nindivmax),
     &     distcell(nindivmax),ttemp(2,nindivmax),
     &     distcelltemp(nindivmax),n(nclassmax,nlocmax,nallmax),
     &     ntemp(nclassmax,nlocmax,nallmax),
     &     a(nallmax),ptemp(nallmax),effcl(nclassmax),iclv(nclassmax),
     &     fa(nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     cellclass(nppmax),listcell(nppmax),cellclasshost(nppmax)



      write(6,*) '          *****************************'
      write(6,*) '          ***    MCMC inference     ***'
      write(6,*) '          *****************************'


      nindiv = nindivmax
      nloc = nlocmax

**************************
*     read data
**************************
C       open(10,file=files)
C       do iindiv=1,nindiv
C c         write(*,*) 'coucou'
C          read(10,*) s(1,iindiv), s(2,iindiv)
C       enddo
C       close(10)
      
C       call limit(nindiv,nindivmax,s,xlim,ylim,dt)
C c      xlim(1) = 0
C c      xlim(2) = 1
C c      ylim(1) = 0
C c      ylim(2) = 1
        
C       open(10,file=filenall)
C       do iloc=1,nloc
C          read(10,*) nall(iloc)
C       enddo
C       close(10)


C       open(10,file=filez)
C       do iindiv=1,nindiv
C          read(10,*) (z(iindiv,iloc),iloc=1,2*nloc)
C       enddo
C       close(10)


*     look for smallest rectangle enclosing the spatial domain
      call limit(nindiv,nindivmax,s,xlim,ylim,dt)
c      write(*,*) 'fin de limit'
c      write(*,*) 's=',s
c      write(*,*) 'z=',z
c      write(*,*) 'nall=',nall


*     Ouverture des fichiers pour l'ecriture des sorties
      open(9,file=filelambda)
      open(10,file=filenpp)
      open(11,file=filenclass)
      open(12,file=fileu)
      open(13,file=filec)
      open(14,file=filef)
      open(15,file=filefa)
      open(16,file=filedrift)
      open(17,file=filelpp)
      open(18,file=filell)
      open(19,file=filet)

c       write(*,*) 'fin de l ouverture'

************************
*     Initialization
************************
      du = sqrt((xlim(2)-xlim(1))*(ylim(2)-ylim(1))/nindiv)
      lambda = lambdamax*ranf()

c$$$      npp = 1+ignpoi(lambda)
c$$$      lambda = lambdamax
c$$$      nclass = rpriornclass(mu,nclassmin,nclassmax)

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

      call rpriorc(npp,nppmax,nclass,c)

      call calccell(nindiv,nindivmax,t,npp,nppmax,u,indcell,distcell)
      
      call rpriordrift(nclass,nclassmax,drift,fmodel)

      call rpriorfa(nloc,nlocmax,nall,nallmax,fa,fmodel)
       
      call rpriorf(nclass,nclassmax,nloc,nlocmax,nall,nallmax,f)

      


************************
* mise a jour iterative 
************************

c      write(*,*) 'Starting updates'
      do ichain=1,nchain
*     ecriture dans les fichiers (tous les stepw)
          if(mod(ichain,stepw) .eq. 0) then
 100         format(f7.3,' %')
*             write(*,*) ''
             write(6,100)float(ichain)/float(nchain)*100.
             write(9,*) lambda
             write(10,*) npp
             write(11,*) nclass
 1000        format (2(1x,e15.8,1x))
             write(12,1000) (u(1,ipp),u(2,ipp), ipp=1,nppmax)
             do ipp=1,nppmax
                write(13,*) c(ipp)
             enddo
c$$$ 2000        format (2(e10.5,1x))
c$$$             write(14,2000) (((f(iclass,iloc,iall),iclass=1,nclassmax),
c$$$     &            iall=1,nallmax),iloc=1,nloc)
 2000        format (300(1x,e15.8,1x))
             do iloc=1,nlocmax
                do iall=1,nallmax
                   write(14,2000) (f(iclass,iloc,iall),
     &                  iclass=1,nclassmax)
                enddo
             enddo
 3000        format (300(1x,e15.8,1x))
             write(15,3000) ((fa(iloc,iall),iall=1,nallmax),iloc=1,nloc)
             write(16,2000) (drift(iclass),iclass=1,nclassmax)
             write(17,*) lpp(lambda,z,nclass,npp,drift,f,fa,c,nppmax,
     &            nindivmax,nlocmax2,nclassmax,nlocmax,nallmax,indcell,
     &            nall,fmodel)
             write(18,*) ll(z,nindivmax,nlocmax,nlocmax2,nall,nclassmax,
     &     nallmax,nppmax,c,f,indcell)
             if(dt .gt. 1.e-30) then 
                write(19,1000) (t(1,iindiv),t(2,iindiv),iindiv=1,nindiv)
             endif
          endif


*     update lambda
*          write(*,*) 'update lambda'
          lambda = rpostlamb(lambdamax,npp)

          if(fmodel .eq. 1) then 
*     update drift
*              write(*,*) 'update drift'
             call  upddrift(nclass,nclassmax,nlocmax,nallmax,nall,
     &            f,fa,drift)
*     update fa
*             write(*,*) 'update fa'
             call updfa(nclass,nclassmax,nlocmax,nallmax,nall,
     &            f,fa,drift)
          endif

*     update f
*          write(*,*) 'update f'
          call rpostf2(nclass,nclassmax,nloc,nlocmax,nall,nallmax,
     &         f,fa,drift,
     &         nindiv,nindivmax,nlocmax2,z,npp,nppmax,c,indcell,
     &         n,a,ptemp)
          
*     update c 
*          write(*,*) 'update c'
          call updc(npp,nppmax,c,ctemp,z,nindiv,nindivmax,nloc,nlocmax,
     &         nlocmax2,nallmax,nclass,nclassmax,f,indcell)


          if(spatial .eq. 1) then 
*     update u et mise a jour de indcell et distcell
*             write(*,*) 'update u'
             call updurw(npp,nppmax,c,u,z,nindiv,nindivmax,nloc,nlocmax,
     &            nlocmax2,nallmax,nclass,nclassmax,f,indcell,distcell,
     &            indcelltemp,distcelltemp,t,xlim,ylim,du)
             
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.e-30) then 
*                write(*,*) 'update t'
                call updt(npp,nppmax,nindiv,
     &               nindivmax,nloc,nlocmax,nlocmax2,nallmax,nclassmax,
     &               t,ttemp,dt,s,c,indcell,distcell,
     &               indcelltemp,distcelltemp,u,z,f)
             endif

*     birth/death des points du pp
*             write(*,*) 'update npp'
             call bdpp(nindiv,nindivmax,u,c,utemp,ctemp,
     &            nclass,nclassmax,nloc,nlocmax,
     &            nlocmax2,nallmax,npp,nppmax,z,f,t,xlim,ylim,indcell,
     &            distcell,indcelltemp,distcelltemp,lambda)
          endif


*     birth/death de classes vides
         if(kfix .eq. 0) then 
*            write(*,*) 'update nclass'
* split/merge avec prop de f selon cond. complete dans les deux sens
            if(fmodel .eq. 0) then 
*     avec drift=0.5 et fa=1 pour court-circuiter le F-model
               call bdclass7bis(nclass,nclassmin,nclassmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &              cellclasshost,n,ntemp)
            else
               call bdclass7(nclass,nclassmin,nclassmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &              cellclasshost,n,ntemp)
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

       
       write(6,*) '          ************************************'
       write(6,*) '          ***    End of MCMC inference     ***'
       write(6,*) '          ************************************'
      end
      
c      include 'sub.f'
      

