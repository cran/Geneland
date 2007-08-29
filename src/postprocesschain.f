      subroutine  postprocesschain(nxdommax,nydommax,burn,
     &     npopmax,nppmax,
     &     nlocmax,nindiv,nloc,nallmax,xlim,ylim,dt,nchain,filenpp,
     &     filenpop,fileu,filec,filef,filefperm,filedom,filedomperm,
     &     s,u,c,dom,domperm,coorddom,indvois,distvois,f,orderf)
      implicit none

      character*255 fileu,filec,filenpp,
     &     filenpop,filedom,filedomperm,filef,filefperm
      
      integer nchain,npp,npop,ichain,nindiv,nxdommax,
     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indvois,
     &     ipop,nloc, nlocmax,nallmax,ijunk,orderf,ipopperm,burn
      real s,u,xlim,ylim,coorddom,dom,domperm,distvois,f,dt

*     dimensionnement 
      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
     &     domperm(nxdommax*nydommax,npopmax),
     &     coorddom(2,nxdommax*nydommax),indvois(nxdommax*nydommax),
     &     distvois(nxdommax*nydommax),f(npopmax),
     &     orderf(npopmax)
      

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

      call limit(nindiv,s,xlim,ylim,dt)


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
            do ipop=1,npopmax
               dom(idom,ipop) = 0.
               domperm(idom,ipop) = 0.
            enddo
            idom = idom + 1
         enddo
      enddo

      open(9,file=filenpop)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      open(13,file=filef)
      open(14,file=filefperm)
      open(15,file=filedom)
      open(16,file=filedomperm)


      do ichain=1,nchain
c 100     format(f7.3,' %')
c         write(6,100)float(ichain)/float(nchain)*100.
         
         read(9,*) npop
         
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         
         read(13,*) f
         call indexx(npop,npopmax,f,orderf)
         
         write(14,2000) (f(orderf(ipop)),ipop=1,npopmax)
         
         do ijunk=2,nallmax*nloc
            read(13,*) f
            write(14,2000) (f(orderf(ipop)),ipop=1,npopmax)
         enddo
         
         
         
         call calccell(nxdommax*nydommax,coorddom,
     &        npp,nppmax,u,indvois,distvois)
         
         if(ichain .gt. burn)  then
            do idom=1,nxdommax*nydommax
               ipop = c(indvois(idom))
               ipopperm = orderf(c(indvois(idom)))
               dom(idom,ipop) =  dom(idom,ipop) + 1.
               domperm(idom,ipopperm) = domperm(idom,ipopperm) + 1.
            enddo
         endif
      enddo


      do idom=1,nxdommax*nydommax
         do ipop=1,npopmax
            dom(idom,ipop) = dom(idom,ipop)/float(nchain-burn)
            domperm(idom,ipop) = domperm(idom,ipop)/
     &           float(nchain-burn)
         enddo
      enddo

 2000 format (1000(f15.3,1x))

      do idom=1,nxdommax*nydommax
         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
     &        (dom(idom,ipop), ipop=1,npopmax)
         write(16,2000)  coorddom(1,idom),  coorddom(2,idom), 
     &        (domperm(idom,ipop), ipop=1,npopmax)
      enddo
      
c      write(*,*) coorddom

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

      integer iit,ipp,iindiv,ipop
      character*255 fileu,filec,filenpp,filenpop
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
      write(6,*) '      **********************************************'
      write(6,*) '      *  Computing posterior probabilities          '
      write(6,*) '      *  of population membership for individuals   '
      write(6,*) '      **********************************************'
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
c$$$      write(*,*) 'burn=',burn,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)

*     sequentially processes states of the chain
     
      do iit=1,nit
c10000    format(f7.3,' %')
c         write(6,10000)float(iit)/float(nit)*100.
         
c$$$         write(*,*) 'nit=',nit
c$$$         write(*,*) 'iit=',iit
c$$$         write(*,*) 'npp=',npp, '\n' 

         read(10,*) npp
         do ipp=1,nppmax
c             write(*,*) 'ipp=',ipp
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         
         if(iit .gt. burn)  then
            call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
            do iindiv=1,nindiv
               ipop = c(indcell(iindiv))
               pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
            enddo
         endif 

c$$$         write(*,*) 'nit=',nit
c$$$         write(*,*) 'iit=',iit
c$$$         write(*,*) 'npp=',npp, '\n'

      enddo
     
c      write(*,*) 'sortie de la boucle'
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
c$$$      write(*,*) 'burn=',burn,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      end subroutine  pppmindiv


      
 
 



