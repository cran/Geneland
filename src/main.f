***********************************************************************  
*     Simulation selon le posterior des parametres 
*     du modele de melange de genotypes par une Voronoi cachee
*
*     Routines de simulations 
*     de ~/mathlib/randlib-1.3.f/src/all.f
*
*
***********************************************************************
 

      subroutine mcmc(s,z,nall,ploidy,filelambda,
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
     &     nall,nallmax,z,ploidy
      double precision s

*     parametres des priors 
      double precision lambdamax,dt

*     les parametres a inferer (et valeurs maximales)
      integer npp,nppmax,nclass,nclassmin,nclassmax,
     &     c,ctemp
      double precision lambda,u,utemp,f,t,fa,drift,ftemp,drifttemp

*     variables de travail
      integer iloc,iindiv,ichain,nchain,
     &     ipp,iclass,iall,stepw,indcell,indcelltemp,
     &     n,effcl,iclv,cellclass,listcell,
     &     cellclasshost,ntemp,fmodel,kfix,spatial
      double precision ptemp,xlim(2),ylim(2),ggrunif,rpostlamb,
     &     distcell,du,distcelltemp,a,ttemp,lpp,ll
      character*256 files,filez,filef,filenall,filenpp,
     &     filelambda,filenclass,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet
      character*1 prchar
 
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

      parameter(prchar='*')



      call intpr(' ',-1,':',0) 
      call intpr(' ',-1,':',0) 
      call intpr('************************************',-1,'*',0)
      call intpr('***  Starting MCMC inference     ***',-1,'*',0)
      call intpr('************************************',-1,'*',0)

************************************************************************
*     init RNG
      call rndstart()

      nindiv = nindivmax
      nloc = nlocmax


*     look for smallest rectangle enclosing the spatial domain
      call limit(nindiv,nindivmax,s,xlim,ylim,dt)


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


************************
*     Initialization
************************
      du = dsqrt((xlim(2)-xlim(1))*(ylim(2)-ylim(1))/dble(nindiv))
      lambda = lambdamax*ggrunif(0.d0,1.d0)

      do iindiv=1,nindiv
         t(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         t(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
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

      call rpriorfa(nloc,nlocmax,nall,nallmax,fa,fmodel,ptemp)
       
      call rpriorf(nclass,nclassmax,nloc,nlocmax,nall,nallmax,f,ptemp)


      
************************
* mise a jour iterative 
************************
      call intpr(' ',-1,':',0) 
      call intpr(' ',-1,':',0) 
      call intpr('Percentage of iterations:',-1,':',0) 
      call intpr('-------------------------',-1,':',0) 

      do ichain=1,nchain
*     ecriture dans les fichiers (tous les stepw)
          if(mod(ichain,stepw) .eq. 0) then
 100         format(f7.3,' %')
*             write(*,*) ''
c             write(6,100)float(ichain)/float(nchain)*100.

             call realpr('', -1,float(ichain)/float(nchain)*100,1)

             write(9,*) lambda
             write(10,*) npp
             write(11,*) nclass
 1000        format (2(1x,d15.8,1x))
             write(12,1000) (u(1,ipp),u(2,ipp), ipp=1,nppmax)
             do ipp=1,nppmax
                write(13,*) c(ipp)
             enddo
c$$$ 2000        format (2(e10.5,1x))
c$$$             write(14,2000) (((f(iclass,iloc,iall),iclass=1,nclassmax),
c$$$     &            iall=1,nallmax),iloc=1,nloc)
 2000        format (300(1x,d15.8,1x))
             do iloc=1,nlocmax
                do iall=1,nallmax
                   write(14,2000) (f(iclass,iloc,iall),
     &                  iclass=1,nclassmax)
                enddo
             enddo
 3000        format (300(1x,d15.8,1x))
             write(15,3000) ((fa(iloc,iall),iall=1,nallmax),iloc=1,nloc)
             write(16,2000) (drift(iclass),iclass=1,nclassmax)
             write(17,*) lpp(lambda,z,nclass,npp,drift,f,fa,c,nppmax,
     &            nindivmax,nlocmax2,nclassmax,nlocmax,nallmax,indcell,
     &            nall,fmodel)
             write(18,*) ll(z,nindivmax,nlocmax,nlocmax2,nall,nclassmax,
     &     nallmax,nppmax,c,f,indcell)
             if(dt .gt. 1.d-30) then 
                write(19,1000) (t(1,iindiv),t(2,iindiv),iindiv=1,nindiv)
             endif
          endif


*     update lambda
c          write(*,*) 'update lambda'

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
c          write(*,*) 'update f'
          call rpostf2(nclass,nclassmax,nloc,nlocmax,nall,nallmax,
     &         f,fa,drift,
     &         nindiv,nindivmax,nlocmax2,z,npp,nppmax,c,indcell,
     &         n,a,ptemp,ploidy)
          
*     update c 
c          write(*,*) 'update c'
c          call updc(npp,nppmax,c,ctemp,z,nindiv,nindivmax,nloc,nlocmax,
c     &         nlocmax2,nallmax,nclass,nclassmax,f,indcell)

*     joint update of c anf f
c          write(*,*) 'update c and f'
          if(fmodel .eq. 0) then 
             call  udcf(nclass,nclassmin,nclassmax,f,fa,drift,
     &            nloc,nlocmax,nlocmax2,
     &            nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,
     &            c,ctemp,a,ptemp,ftemp,drifttemp,z,n,ntemp,ploidy)
          else
             call udcf2(nclass,nclassmin,nclassmax,f,fa,drift,
     &            nloc,nlocmax,nlocmax2,
     &            nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,
     &            c,ctemp,a,ptemp,ftemp,drifttemp,z,n,ntemp,ploidy)
          endif


          if(spatial .eq. 1) then 
*     update u et mise a jour de indcell et distcell
c             write(*,*) 'update u'
             call updurw(npp,nppmax,c,u,z,nindiv,nindivmax,nloc,nlocmax,
     &            nlocmax2,nallmax,nclass,nclassmax,f,indcell,distcell,
     &            indcelltemp,distcelltemp,t,xlim,ylim,du,ploidy)
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.d-200) then 
*                write(*,*) 'update t'
                call updt(npp,nppmax,nindiv,
     &               nindivmax,nloc,nlocmax,nlocmax2,nallmax,nclassmax,
     &               t,ttemp,dt,s,c,indcell,distcell,
     &               indcelltemp,distcelltemp,u,z,f,ploidy)
             endif

*     birth/death des points du pp
c             write(*,*) 'update npp'
             call bdpp(nindiv,nindivmax,u,c,utemp,ctemp,
     &            nclass,nclassmax,nloc,nlocmax,
     &            nlocmax2,nallmax,npp,nppmax,z,f,t,xlim,ylim,indcell,
     &            distcell,indcelltemp,distcelltemp,lambda,ploidy)

          endif


*     birth/death de pop
         if(kfix .eq. 0) then 
c            write(*,*) 'begin update nclass'
* split/merge avec prop de f selon cond. complete dans les deux sens
            if(fmodel .eq. 0) then 
*     avec drift=0.5 et fa=1 pour court-circuiter le F-model
               call bdclass8bis(nclass,nclassmin,nclassmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &              cellclasshost,n,ntemp,ploidy)
            else
               call bdclass8(nclass,nclassmin,nclassmax,f,fa,drift,
     &              nloc,nlocmax,nlocmax2,
     &              nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,
     &              ctemp,a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &              cellclasshost,n,ntemp,ploidy)
            endif
c            write(*,*) 'end update nclass'
         endif

      enddo
      call rndend()  
      
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

       
      call intpr('***************************************',-1,'*',0)
      call intpr('***    End of MCMC computations     ***',-1,'*',0)
      call intpr('***************************************',-1,'*',0)

c      write(6,*) '          ************************************'
c      write(6,*) '          ***    End of MCMC inference     ***'
c      write(6,*) '          ************************************'
       
      end

c      include './sub.f'

