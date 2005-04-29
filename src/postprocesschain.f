      subroutine  postprocesschain(nindivmax,nxdommax,nydommax,burn,
     &     nclassmax,nppmax,
     &     nlocmax,nindiv,nloc,nall,nallmax,dt,nchain,filenpp,
     &     filenclass,fileu,filec,filef,filefperm,filedom,filedomperm,
     &     s,u,c,dom,domperm,coorddom,indvois,distvois,f,orderf)
      implicit none

      character*256 files,fileu,filec,filenpp,
     &     filenclass,filedom,filedomperm,filenall,filef,filefperm
      
      integer nchain,npp,nclass,ichain,nindiv,nindivmax,nxdommax,
     &     nydommax,nclassmax,ipp,nppmax,c,ixdom,iydom,idom,indvois,
     &     iclass,nloc, nlocmax,nall,nallmax,ijunk,orderf,iindiv,iloc,
     &     iclassperm,burn
      real s,u,xlim(2),ylim(2),coorddom,dom,domperm,distvois,f,dt

*     dimensionnement 
*     PENSER AU FORMAT D ECRITURE DES FREQUENCES 
*     QUI DEPEND DE NCLASSMAX
      dimension s(2,nindivmax),u(2,nppmax),c(nppmax),
     &     dom(nxdommax*nydommax,nclassmax),
     &     domperm(nxdommax*nydommax,nclassmax),
     &     coorddom(2,nxdommax*nydommax),indvois(nxdommax*nydommax),
     &     distvois(nxdommax*nydommax),nall(nlocmax),f(nclassmax),
     &     orderf(nclassmax)
      

**************************
*     lecture des données
**************************
      write(6,*) '      *****************************************'
      write(6,*) '      *  Computing posterior probabilities     '
      write(6,*) '      *  of population membership for pixels   '
      write(6,*) '      *****************************************'
      write(6,*) ''
      write(6,*) ''
      write(6,*) ''
      write(6,*) ''

      call limit(nindiv,nindivmax,s,xlim,ylim,dt)


*     coordonnées de la grille 
      idom = 1
      do ixdom =1,nxdommax
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nydommax
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           float(ixdom-1)*(xlim(2) - xlim(1))/float(nxdommax-1)
            coorddom(2,idom) = ylim(1) +
     &           float(iydom-1)*(ylim(2) - ylim(1))/float(nydommax-1)
            do iclass=1,nclassmax
               dom(idom,iclass) = 0.
               domperm(idom,iclass) = 0.
            enddo
            idom = idom + 1
         enddo
      enddo

      open(9,file=filenclass)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      open(13,file=filef)
      open(14,file=filefperm)
      open(15,file=filedom)
      open(16,file=filedomperm)


      do ichain=1,nchain
 100     format(f7.3,' %')
c         write(6,100)float(ichain)/float(nchain)*100.
         
         read(9,*) nclass
         
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         
         read(13,*) f
         call indexx(nclass,nclassmax,f,orderf)
         
         write(14,2000) (f(orderf(iclass)),iclass=1,nclassmax)
         
         do ijunk=2,nallmax*nloc
            read(13,*) f
            write(14,2000) (f(orderf(iclass)),iclass=1,nclassmax)
         enddo
         
         
         
         call calccell(nxdommax*nydommax,nxdommax*nydommax,coorddom,
     &        npp,nppmax,u,indvois,distvois)
         
         if(ichain .gt. burn)  then
            do idom=1,nxdommax*nydommax
               iclass = c(indvois(idom))
               iclassperm = orderf(c(indvois(idom)))
               dom(idom,iclass) =  dom(idom,iclass) + 1.
               domperm(idom,iclassperm) = domperm(idom,iclassperm) + 1.
            enddo
         endif
      enddo


      do idom=1,nxdommax*nydommax
         do iclass=1,nclassmax
            dom(idom,iclass) = dom(idom,iclass)/float(nchain-burn)
            domperm(idom,iclass) = domperm(idom,iclass)/
     &           float(nchain-burn)
         enddo
      enddo

 2000 format (1000(f8.3,1x))

      do idom=1,nxdommax*nydommax
         write(15,2000) (dom(idom,iclass), iclass=1,nclassmax)
         write(16,2000) (domperm(idom,iclass), iclass=1,nclassmax)
      enddo
      
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
 
      end subroutine postprocesschain

*
*     posterior probability of population membership for individuals
*
      subroutine  pppmindiv(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,filenpp,fileu,filec,nit,
     &     burn)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npp,c,
     &     nit,burn
      real pmp,distcell,u,s

      integer iit,ipp,iindiv,nppcur,ccur,ipop
      real xlim(2),ylim(2),ucur,dt
      character*256 files,fileu,filec,filenpp,
     &     filenclass,filedom,filedomperm,filenall,filef,filefperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
      write(6,*) '      **********************************************'
      write(6,*) '      *  Computing posterior probabilities          '
      write(6,*) '      *  of population membership for individuals   '
      write(6,*) '      **********************************************'

      
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)

*     sequentially processes states of the chain
      do iit=1,nit
10000    format(f7.3,' %')
c         write(6,10000)float(iit)/float(nit)*100.
         
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         
         if(iit .gt. burn)  then
            call calccell(nindiv,nindiv,s,npp,nppmax,u,indcell,distcell)
            do iindiv=1,nindiv
               ipop = c(indcell(iindiv))
               pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
            enddo
         endif 
      enddo
      
c      write(*,*) 'pmp=',pmp
      do iindiv=1,nindiv 
         do ipop=1,npopmax
            pmp(iindiv,ipop) = pmp(iindiv,ipop)/float(nit-burn)
         enddo
      enddo
c      write(*,*) 'pmp=',pmp
      close(10)
      close(11)
      close(12)
      end subroutine  pppmindiv


      
 
 



