*     limites du rectangle contenant les coordonnees
      subroutine limit(nindiv,nindivmax,s,xlim,ylim,dt)
      implicit none
      integer nindiv,nindivmax
      real s(2,nindivmax),xlim(2),ylim(2),dt
      integer iindiv
      xlim(1) = 1.e+30
      xlim(2) = -1.e+30
      ylim(1) = 1.e+30
      ylim(2) = -1.e+30
      do iindiv=1,nindiv
         xlim(1) = amin1(s(1,iindiv),xlim(1))
         xlim(2) = amax1(s(1,iindiv),xlim(2))
         ylim(1) = amin1(s(2,iindiv),ylim(1))
         ylim(2) = amax1(s(2,iindiv),ylim(2))
      enddo
      xlim(1) = xlim(1) - dt*.5
      xlim(2) = xlim(2) + dt*.5
      ylim(1) = ylim(1) - dt*.5
      ylim(2) = ylim(2) + dt*.5
      end

c$$$*     Poisson translatee de 1 tronquee
c$$$      integer function rpriornpop(mu,npopmin,npopmax)
c$$$      implicit none
c$$$      integer npopmin,npopmax
c$$$      real mu
c$$$      integer n,ignpoi
c$$$      n = npopmax +1
c$$$      do while((n .lt. npopmin) .or. (n .gt. npopmax))
c$$$         n = 1+ ignpoi(mu)
c$$$      enddo
c$$$      rpriornpop = n
c$$$      end
      

*     points uniformes dans [0,1]x[0,1]
      subroutine rprioru(npp,nppmax,xlim,ylim,u)
      implicit none
      integer npp,nppmax
      real u(2,nppmax),ranf,xlim(2),ylim(2)
      integer i
      do i=1,npp
         u(1,i) = xlim(1)+(xlim(2)-xlim(1))*ranf()
         u(2,i) = ylim(1)+(ylim(2)-ylim(1))*ranf()
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1, nppmax
            u(1,i) = -999.
            u(2,i) = -999.
         enddo
      endif
      end

*     affectation dans les popes selon une loi uniforme
      subroutine rpriorc(npp,nppmax,npop,c)
      implicit none
      integer npp,nppmax,npop,c(nppmax)
      real ranf
      integer i
      do i=1,npp
         c(i) = 1+ int(aint(float(npop)*ranf()))
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1,nppmax
            c(i) = -999
         enddo
      endif
      end

*     init de la dérive 
*     selon un prior uniforme sur [0,1]
      subroutine rpriordrift(npop,npopmax,drift,fmodel)
      implicit none
      integer npop,npopmax,fmodel
      real drift(npopmax)
      integer ipop
      real ranf
      if(fmodel .eq. 0) then
         do ipop=1,npop
            drift(ipop) = 0.5
         enddo
      else
         do ipop=1,npop
            drift(ipop) = ranf()
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
      subroutine dirichlet1(nall,nallmax,p)
      implicit none
      integer nall,nallmax
      real p(nallmax)
      integer i
      real s,genexp
      s = 0.
      do i=1,nall
         p(i) = genexp(1.)
         s = s + p(i)
      enddo
      do i=1,nall
         p(i) =  p(i)/s
      enddo
      if(nallmax .gt. nall) then
         do i=nall+1,nallmax
            p(i) =  -1
         enddo
      endif
      end

*     simulation d'une Dirichlet(a1,...,an)
      subroutine dirichlet(n,nmax,a,p)
      implicit none
      integer n,nmax
      real a(nmax),p(nmax)
      integer i
      real s,gengam
      s = 0.
      do i=1,n
         p(i) = 0.
         do while(p(i) .lt. 1e-37) 
            p(i) = gengam(1.,a(i))
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
      end

      
      subroutine rank(n,nmax,x,p)
      implicit none 
      integer n,nmax,p(nmax)
      real x(nmax)
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


*     tirage des frequences dans toutes les popes
*     a tous les locus
      subroutine rpriorf(npop,npopmax,nloc,nlocmax,nall,nallmax,f,
     &     ptemp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),nallmax
      real f(npopmax,nlocmax,nallmax)
      integer k,l,i
      real ptemp(nallmax)
      do  k=1,npop
         do l=1,nloc
            call dirichlet1(nall(l),nallmax,ptemp)
            do i=1,nallmax
               f(k,l,i)  = ptemp(i)
            enddo
         enddo
      enddo
      if(npopmax .gt. npop) then
         do k=npop+1, npopmax
            do l=1,nloc
               do i=1,nallmax
                  f(k,l,i)  = -999
               enddo
            enddo
         enddo
      endif
      end subroutine rpriorf

*     tirage des frequences dans la pop ancestrale 
*     a tous les locus
      subroutine rpriorfa(nloc,nlocmax,nall,nallmax,fa,fmodel,ptemp)
      implicit none
      integer nloc,nlocmax,nall(nlocmax),nallmax,fmodel
      real fa(nlocmax,nallmax)
      integer l,i
      real ptemp(nallmax)
      if(fmodel .eq. 0) then
         do l=1,nloc
            do i=1,nallmax
               fa(l,i)  = 1
            enddo
         enddo
      else
        do l=1,nloc
           call dirichlet1(nall(l),nallmax,ptemp)
           do i=1,nallmax
               fa(l,i)  = ptemp(i)
            enddo
         enddo 
      endif
      end subroutine rpriorfa


*     Mise a jour gibbsienne de f
*     prior p(f) Dirichlet(1,...,1) 
*     p(f|...) Dirichlet(1+ n1,..., 1+np)
*     ni = nbre d'alleles observes
*     (Cf Green 95, p738) 
*     (et on fait comme si on avait 2*n individus haploides) 
      subroutine rpostf(npop,npopmax,nloc,nlocmax,nall,nallmax,f,
     &     nindiv,nindivmax,nlocmax2,z,nppmax,c,indcell,n,a,f11temp)
      implicit none 
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),nallmax,
     &     nindiv,nindivmax,nlocmax2,nppmax,c(nppmax),
     &     indcell(nindivmax)
      real f(npopmax,nlocmax,nallmax)
      integer ipop,iloc,iindiv,iall,n(npopmax,nlocmax,nallmax),
     &     z(nindivmax,nlocmax2)
      real a(nallmax),f11temp(nallmax)

*     comptage des effectifs
      do ipop = 1,npop
         do iloc = 1,nloc
            do iall =1,nall(iloc)
               n(ipop,iloc,iall)=0
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
            do iall =1,nall(iloc)
               a(iall) = 1+float(n(ipop,iloc,iall))
            enddo
            call dirichlet(nall(iloc),nallmax,a,f11temp)
            do iall =1,nall(iloc)
               f(ipop,iloc,iall) = f11temp(iall)
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
      subroutine rpostf2 (npop,npopmax,nloc,nlocmax,nall,nallmax,
     &     f,fa,drift,
     &     nindiv,nindivmax,nlocmax2,z,nppmax,c,indcell,n,a,ptemp,
     &     ploidy)
      implicit none 
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),nallmax,
     &     nindiv,nindivmax,nlocmax2,nppmax,c(nppmax),
     &     indcell(nindivmax),ploidy
      real f(npopmax,nlocmax,nallmax),fa(nlocmax,nallmax),
     &     drift(npopmax)
      integer ipop,iloc,iindiv,iall,n(npopmax,nlocmax,nallmax),
     &     z(nindivmax,nlocmax2)
      real a(nallmax),ptemp(nallmax)

c$$$      write(*,*) ''
c$$$      write(*,*) ''
c$$$      write(*,*) 'dans rpostf2'
c$$$      write(*,*) 'npop=',npop
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nlocmax=',nlocmax
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'f=',f
c$$$      write(*,*) 'fa=',fa
c$$$      write(*,*) 'drift=',drift
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nindivmax=',nindivmax
c$$$      write(*,*) 'nlocmax2=',nlocmax2
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'npp=',npp
c$$$      write(*,*) 'nppmax=',nppmax
c$$$      write(*,*) 'c=',c
c$$$      write(*,*) 'indcell=',indcell
c$$$      write(*,*) 'n=',n
c$$$      write(*,*) 'a=',a
c$$$      write(*,*) 'ptemp=',ptemp
      
*     comptage des effectifs
      do ipop = 1,npop
         do iloc = 1,nloc
            do iall =1,nall(iloc)
               n(ipop,iloc,iall)=0
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
            do iall =1,nall(iloc)
c               write(*,*) 'iall=',iall
               if(ploidy .eq. 1) then
                  a(iall) = fa(iloc,iall)*(1/drift(ipop)-1) + 
     &              float(n(ipop,iloc,iall))/2
               endif
               if(ploidy .eq. 2) then
                  a(iall) = fa(iloc,iall)*(1/drift(ipop)-1) + 
     &              float(n(ipop,iloc,iall))
               endif
c               write(*,*) 'a(',iall,')=',a(iall) 
c               write(*,*) 'fa(',iloc,',',iall,')=', fa(iloc,iall)
c               write(*,*) 'drift(',ipop,')=',drift(ipop)
c               write(*,*) 'n(',ipop,',',iloc,',',iall,')=', 
c     &              n(ipop,iloc,iall)
            enddo
c            write(*,*) 'fa=',fa
c            write(*,*) 'drift=',drift
c            write(*,*) 'a=',a
            call dirichlet(nall(iloc),nallmax,a,ptemp)
            do iall =1,nall(iloc)
               f(ipop,iloc,iall) = ptemp(iall)
            enddo
         enddo
      enddo
      end subroutine rpostf2





*
*     Mise à jour M-H des freq allelique de la pop ancestrale 
*     (F-Model de Falush et al.)
*     fa admet un prior Dirichlet(1,...,1)
      subroutine updfa(npop,npopmax,nlocmax,nallmax,nall,f,fa,drift)
      implicit none 
      integer npop,npopmax,nlocmax,nallmax,nall(nlocmax)
      real f(npopmax,nlocmax,nallmax),fa(nlocmax,nallmax),
     &     drift(npopmax)
      integer iloc,iall1,iall2,ipop
      real delta,ranf,mysnorm,sigdelta,fa1,fa2,ratio,lratio,algama,q,u
      parameter(sigdelta = 0.05) 
      
*     boucle sur les loci
      do iloc=1,nlocmax

*     tirage des deux formes alleliques dont les freq seront 
*     mises à jour 
         iall1 = 1+ int(aint(float(nall(iloc))*ranf()))
         iall2 = iall1 

         do while(iall2 .eq. iall1)
c            write(*,*) 'dans le while'
            iall2 = 1+ int(aint(float(nall(iloc))*ranf()))
         enddo

*     tirage de l'increment
         delta = mysnorm()*sigdelta

*     perturbation des deux freq
         fa1 = fa(iloc,iall1) + delta
         fa2 = fa(iloc,iall2) - delta
         if(((fa1 .gt. 1e-37) .and. (1-fa1 .gt. 1e-37)) .and.
     &      ((fa2 .gt. 1e-37) .and. (1-fa2 .gt. 1e-37))) then 
*     calcul du log du ratio 
            lratio = 0.
            do ipop = 1,npop
               q = (1.-drift(ipop))/drift(ipop)
               lratio = lratio 
     &              + algama(fa(iloc,iall1)*q)-algama(fa1*q)
     &              + algama(fa(iloc,iall2)*q)-algama(fa2*q)
     &              + delta*q
     &              *log(f(ipop,iloc,iall1)/f(ipop,iloc,iall2))
            enddo
            lratio = amin1(0.,lratio) 
            ratio = exp(lratio)

c$$$            write(*,*) 'delta=', delta
c$$$            write(*,*) 'fa1=',fa1
c$$$            write(*,*) 'fa2=',fa2
c$$$            write(*,*) 'q=',q
c$$$            write(*,*) 'fa(iloc,iall2)*q=',fa(iloc,iall2)*q
c$$$            write(*,*) 'gamma(fa(iloc,iall1)*q)',gamma(fa(iloc,iall1)*q)
c$$$            write(*,*) 'ratio=',ratio
c$$$            write(*,*) 'delta*q=',delta*q
c$$$            write(*,*) ''

            u = ranf()
            if(u .le. ratio) then 
               fa(iloc,iall1) = fa1 
               fa(iloc,iall2) = fa2
            endif
         endif 
      enddo
      end subroutine updfa
      

*
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. uniforme sur chaque composante
      subroutine upddrift(npop,npopmax,nlocmax,nallmax,nall,
     &     f,fa,drift)
      implicit none 
      integer npop,npopmax,nlocmax,nallmax,nall(nlocmax)
      real drift(npopmax),f(npopmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax)
      integer ipop,iloc,iall
      real d,q,qtemp,sigdelta,ratio,lratio,shape1,shape2,sall,mysnorm,
     &     algama,u,ranf
      parameter(sigdelta = 0.01,shape1=2.,shape2=20.) 

*     boucle sur les popes
      do ipop=1,npop
*     proposition nouvelle valeur
         d = drift(ipop) + mysnorm()*sigdelta
         q = (1-drift(ipop))/drift(ipop)
         qtemp = (1-d)/d
         if((d .gt. 1e-37) .and. (1-d .gt. 1e-37)) then 

*     calcul du log du ratio
c     prior uniforme
            lratio = 0 
c     prior beta(shape1,shape2)
            lratio = (shape1-1)*alog(d/drift(ipop)) + 
     &           (shape2-1)*alog((1-d)/(1-drift(ipop)))
            do iloc=1,nlocmax
               sall = 0.
               do iall = 1,nall(iloc)
                  sall = sall + algama(fa(iloc,iall)*q)-
     &                 algama(fa(iloc,iall)*qtemp) +
     &                 fa(iloc,iall)*(qtemp-q)*alog(f(ipop,iloc,iall))
               enddo
               lratio = lratio + sall + (algama(qtemp)-algama(q))
            enddo
 

            lratio = amin1(0.,lratio)
            ratio = exp(lratio)
            u = ranf()
            if(u .le. ratio) then 
               drift(ipop) = d 
            endif
         endif
      enddo
      end subroutine upddrift
*
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. uniforme sur chaque composante
      subroutine upddrift2(npop,npopmax,nlocmax,nallmax,nall,
     &     f,fa,drift)  
      implicit none 
      integer npop,npopmax,nlocmax,nallmax,nall(nlocmax)
      real drift(npopmax),f(npopmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax)
      integer ipop,iloc,iall
      real d,q,qtemp,sigdelta,ratio,lratio,alpha,sall,mysnorm,
     &     algama,u,ranf
      parameter(sigdelta = 0.01,alpha=5000) 

*     boucle sur les popes
      do ipop=1,npop
*     proposition nouvelle valeur
         d = drift(ipop) + mysnorm()*sigdelta
         q = (1-drift(ipop))/drift(ipop)
         qtemp = (1-d)/d
         if((d .gt. 1e-37) .and. (1-d .gt. 1e-37)) then 

*     calcul du log du ratio
            lratio = 0 
c     decommenter la ligne suivante pour avoir un prior exponentiel tronqué
c     sinon le prior est uniforme
c            lratio = -alpha*(d-drift(ipop))
            do iloc=1,nlocmax
               sall = 0.
               do iall = 1,nall(iloc)
                  sall = sall + algama(fa(iloc,iall)*q)-
     &                 algama(fa(iloc,iall)*qtemp) +
     &                 fa(iloc,iall)*(qtemp-q)*alog(f(ipop,iloc,iall))
               enddo
               lratio = lratio + sall + (algama(qtemp)-algama(q))
            enddo


            lratio = amin1(0.,lratio)
            ratio = exp(lratio)
            u = ranf()
            if(u .le. ratio) then 
               drift(ipop) = d 
            endif
         endif
      enddo
      end subroutine upddrift2




*     recherche la cellule de chaque individu
*     stockage des indices dans indcell
*     stockage des carres des distances dans distcell
      subroutine calccell(nindiv,nindivmax,s,npp,nppmax,u,
     &     indcell,distcell)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax)
      real distcell(nindivmax)
      real s(2,nindivmax),u(2,nppmax)
      integer iindiv,ipp
      real d
      do iindiv=1,nindiv
         indcell(iindiv) = -999
         distcell(iindiv) = 1.e+36
         do ipp=1,npp
            d = (s(1,iindiv)-u(1,ipp))**2 +
     &           (s(2,iindiv)-u(2,ipp))**2
            if( d .lt. distcell(iindiv) ) then 
               indcell(iindiv) = ipp
               distcell(iindiv) = d
            endif
         enddo
      enddo
      if(nindivmax .gt. nindiv) then 
         do iindiv=nindiv+1,nindivmax
            indcell(iindiv) = -999
            distcell(iindiv) = -999.
         enddo
      endif
      end




*     mise a jour de indcell et distcell
*     apres le deplacement d'un point de u (celui d'indice j)
      subroutine vormove(nindiv,nindivmax,s,npp,nppmax,u,
     &     indcell,distcell,indcelltemp,distcelltemp,j)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),j
      real s(2,nindivmax),u(2,nppmax),distcell(nindivmax),
     &     distcelltemp(nindivmax),d
      integer iindiv,ipp

C       write(6,*) 'debut de  vormove'
C       write(6,*) 'j =', j
C       write(6,*) 'indcell',indcell
C       write(6,*)'distcell',distcell
C       write(6,*) 'indcelltemp',indcelltemp
C       write(6,*)'distcelltemp',distcelltemp

      do iindiv=1,nindiv
         if(indcell(iindiv) .eq. j) then 
*     pour les indiv qui etaient dans la cellule j on cherche
*     la nouvelle cellule
            d = 3.e+37
            indcelltemp(iindiv) = -999
            distcelltemp(iindiv) = 3.e+37
            do ipp=1,npp
               d= (s(1,iindiv)-u(1,ipp))**2+(s(2,iindiv)-u(2,ipp))**2
               if( d .lt. distcelltemp(iindiv) ) then 
                  indcelltemp(iindiv) = ipp
                  distcelltemp(iindiv) = d
               endif
            enddo
*     pour les autres indiv on regarde si le nouveau uj s'est intercale
         else
            d = (s(1,iindiv)-u(1,j))**2+(s(2,iindiv)-u(2,j))**2
            if(d .lt. distcell(iindiv)) then
               indcelltemp(iindiv) = j
               distcelltemp(iindiv) = d
            else
               indcelltemp(iindiv) = indcell(iindiv)
               distcelltemp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo
         if(nindivmax .gt. nindiv) then
            do iindiv=nindiv+1,nindivmax
               indcelltemp(iindiv) = -999
               distcelltemp(iindiv) = -999.
            enddo
         endif
         

c$$$         call calccell(nindiv,nindivmax,s,npp,nppmax,u,
c$$$     &        indcelltemp2,distcelltemp2)
c$$$         do iindiv=1,nindiv
c$$$            if((indcelltemp2(iindiv) .ne. indcelltemp(iindiv)) .or. 
c$$$     &         (distcelltemp2(iindiv) .ne. distcelltemp(iindiv)))  then
c$$$               write(6,*) 'fin de  vormove'
c$$$               write(6,*) 'j =', j
c$$$               write(6,*) 'iindiv=',iindiv
c$$$               write(6,*) 'indcell',indcell
c$$$               write(6,*)'distcell',distcell
c$$$               write(6,*) 'indcelltemp',indcelltemp
c$$$               write(6,*)'distcelltemp',distcelltemp
c$$$               write(6,*)'indceltmp2',indcelltemp2
c$$$               write(6,*)'distceltmp2',distcelltemp2
c$$$               stop
c$$$            endif
c$$$         enddo
C          write(6,*) 'fin de  vormove'
C          write(6,*) 'j =', j
C          write(6,*) 'indcell',indcell
C          write(6,*)'distcell',distcell
C          write(6,*) 'indcelltemp',indcelltemp
C          write(6,*)'distcelltemp',distcelltemp
      end



*     mise a jour de indcell et distcell
*     apres naissance d'un point de u 
      subroutine voradd(s,utemp,
     &     indcell,distcell,indcelltemp,distcelltemp,
     &     nindiv,nindivmax,npp,nppmax)
      implicit none 
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),iindiv
      real s(2,nindivmax),distcell(nindivmax),
     &     distcelltemp(nindivmax),d,utemp(2,nppmax)
      
      do iindiv =1,nindiv
*     est-ce que le nouveau point s'est intercale ?
         d = (s(1,iindiv)-utemp(1,npp+1))**2+
     &        (s(2,iindiv)-utemp(2,npp+1))**2
         if(d .lt. distcell(iindiv)) then 
            distcelltemp(iindiv) = d
            indcelltemp(iindiv) = npp+1
         else
            distcelltemp(iindiv) = distcell(iindiv)
            indcelltemp(iindiv) = indcell(iindiv) 
         endif
      enddo
      if(nindivmax .gt. nindiv) then
         do iindiv=nindiv+1,nindivmax
            indcelltemp(iindiv) = -999
            distcelltemp(iindiv) = -999.
         enddo
      endif
      end 

 


     
*     mise a jour de indcell et distcell
*     apres mort d'un point de u 
      subroutine vorrem(s,utemp,ipprem,
     &     indcell,distcell,indcelltemp,distcelltemp,
     &     nindiv,nindivmax,npp,nppmax)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),ipprem,iindiv
      real s(2,nindivmax),utemp(2,nppmax),
     &     distcell(nindivmax),distcelltemp(nindivmax),d
      integer ipp
      
      do iindiv =1,nindiv
*     est-ce que le site courant dependait de la cellule disparue ?
         if(indcell(iindiv) .eq. ipprem) then
*     si oui on recherche sa nouvelle cellule parmi celles qui restent
*     (les nouvelles)
            distcelltemp(iindiv) = 3.e+37
            do ipp=1,npp-1
               d = (s(1,iindiv)-utemp(1,ipp))**2+
     &              (s(2,iindiv)-utemp(2,ipp))**2
               if( d .lt. distcelltemp(iindiv) ) then 
                  indcelltemp(iindiv) = ipp
                  distcelltemp(iindiv) = d
               endif
            enddo
         else
            if(indcell(iindiv) .lt. ipprem) then
               indcelltemp(iindiv) = indcell(iindiv)
               distcelltemp(iindiv) = distcell(iindiv)
            else
               indcelltemp(iindiv) = indcell(iindiv) - 1
               distcelltemp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo
      if(nindivmax .gt. nindiv) then
         do iindiv=nindiv+1,nindivmax
            indcelltemp(iindiv) = -999
            distcelltemp(iindiv) = -999.
         enddo
      endif
      end




*
*     Tirage de lambda selon p(lambda|m) 
*     avec 
*     p(m|lambda) Poisson translatee
*     p(lambda) uniforme dans [0,lambdamax]
*     p(lambda|m) gamma tronquee
      real function rpostlamb(lambdamax,m)
      implicit none 
      real lambdamax,gengam
      integer m
*      write(*,*) 'beg rpostlamb'
      rpostlamb = lambdamax + 1
      do while(rpostlamb .gt. lambdamax)
         rpostlamb = gengam(1.,float(m))
      enddo
*      write(*,*) 'end rpostlamb'
      end







*     
*     Mise a jour de c sans modif de npp
*
      subroutine updc(npp,nppmax,c,ctemp,z,nindiv,nindivmax,nloc,
     &     nlocmax,nlocmax2,nallmax,npop,npopmax,f,indcell,ploidy)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,npop,nallmax,npopmax,z(nindivmax,nlocmax2),
     &     indcell(nindivmax),ploidy
      real f(npopmax,nlocmax,nallmax)
      integer ipp,ctemp(nppmax),ignbin,bern
      real ranf,r,alpha,ratio

      do ipp=1,npp
         ctemp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctemp(ipp) = -999
         enddo
      endif
      do ipp=1,npp
         ctemp(ipp) = 1 + int(aint(float(npop)*ranf()))
         r = ratio(z,f,c,ctemp,indcell,indcell,
     &     npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)

         alpha = amin1(1.,r)
         bern = ignbin(1,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
         else 
            ctemp(ipp) = c(ipp)
         endif
      enddo
      end subroutine updc





***********************************************************************
*     joint update of c and f under the Dirichlet model
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf(npop,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,z,n,ntemp,ploidy)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer ipop,ipp,ipop1,ipop2,iloc,iall,ignbin
      real alpha,ranf,lrpf,lratio,llr6
      integer bern

*     init. temp. vector of population membership
      do ipp=1,npp
          ctemp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctemp(ipp) = -999
         enddo
      endif

*     init temp freq
      do ipop = 1,npop
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
            enddo 
         enddo
      enddo

      do ipp=1,npp
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctemp(ipp) = 1 + int(aint(float(npop)*ranf()))
         ipop1 = c(ipp)
         ipop2 = ctemp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nall,nallmax,z,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)


c         write(*,*) 'n=',n
c         write(*,*) 'ntemp=',ntemp

*     sample new frequencies
         call samplef(npop,npopmax,nloc,nlocmax,
     &     nall,nallmax,ipop1,ipop2,f,ftemp,
     &     fa,drift,a,ptemp,ntemp)

c         write(*,*) 'f=',f
c         write(*,*) 'ftemp=',ftemp



*     compute M-H ratio
*     likelihood
         lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,npopmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)

c         write(*,*) 'lratio =',lratio

*     contrib proposals
         lratio = lratio + lrpf(npopmax,nloc,nall,nallmax,n,ntemp,
     &        f,ftemp,ipop1,ipop2)

c          write(*,*) 'lratio =',lratio

         lratio = amin1(0.,lratio)
         alpha = exp(lratio)

c         write(*,*) 'alpha=',alpha

         bern = ignbin(1,alpha)

c         write(*,*) 'bern=',bern

         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  f(ipop1,iloc,iall) = ftemp(ipop1,iloc,iall)
                  f(ipop2,iloc,iall) = ftemp(ipop2,iloc,iall)
               enddo
            enddo
         else 
            ctemp(ipp) = c(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop1,iloc,iall) = f(ipop1,iloc,iall)
                  ftemp(ipop2,iloc,iall) = f(ipop2,iloc,iall)
               enddo
            enddo
         endif
      enddo
      end subroutine udcf
***********************************************************************







***********************************************************************
*     joint update of c and f under the Falush model
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf2(npop,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,z,n,ntemp,ploidy)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer ipop,ipp,ipop1,ipop2,iloc,iall
      real alpha,ranf,lrppf,lratio,llr6
      integer ignbin, bern


*     init. temp. vector of population membership
      do ipp=1,npp
          ctemp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctemp(ipp) = -999
         enddo
      endif

*     init temp freq
      do ipop = 1,npop
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
            enddo 
         enddo
      enddo

      do ipp=1,npp
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctemp(ipp) = 1 + int(aint(float(npop)*ranf()))
         ipop1 = c(ipp)
         ipop2 = ctemp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nall,nallmax,z,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)
         
*     sample new frequencies
         call samplef(npop,npopmax,nloc,nlocmax,
     &        nall,nallmax,ipop1,ipop2,f,ftemp,
     &        fa,drift,a,ptemp,ntemp)
         
*     compute M-H ratio
*     likelihood
         lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &        indcell,npopmax,nlocmax,nallmax,
     &        nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
         
*     contrib prior and proposals
         lratio = lratio + lrppf(npopmax,nloc,nall,nallmax,n,
     &     ntemp,fa,drift,f,ftemp,ipop1,ipop2)

         lratio = amin1(0.,lratio)
         alpha = exp(lratio)
         bern = ignbin(1,alpha)
         
c         write(*,*) 'bern=',bern

         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  f(ipop1,iloc,iall) = ftemp(ipop1,iloc,iall)
                  f(ipop2,iloc,iall) = ftemp(ipop2,iloc,iall)
               enddo
            enddo
         else 
            ctemp(ipp) = c(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop1,iloc,iall) = f(ipop1,iloc,iall)
                  ftemp(ipop2,iloc,iall) = f(ipop2,iloc,iall)
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
      subroutine updurw(npp,nppmax,c,u,z,nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,nallmax,npopmax,f,indcell,distcell,
     &     indcelltemp,distcelltemp,
     &     s,xlim,ylim,du,ploidy)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,nallmax,npopmax,z(nindivmax,nlocmax2),
     &     indcell(nindivmax),ploidy
      real u(2,nppmax),f(npopmax,nlocmax,nallmax),
     &     distcell(nindivmax),s(2,nindivmax),xlim(2),ylim(2),du
      integer ipp,iindiv,ignbin,bern,indcelltemp(nindivmax)
      real utemp(2,nppmax),ranf,r,alpha,distcelltemp(nindivmax),
     &     surf,surftemp,dx,dy,ratio

*     initialisation du tableau temporaire
      do ipp=1,npp
         utemp(1,ipp) = u(1,ipp)
         utemp(2,ipp) = u(2,ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            utemp(1,ipp) = -999.
            utemp(2,ipp) = -999.
         enddo
      endif

c      write(*,*) 'npp=', npp
c      write(*,*) 'u=', u

      do ipp=1,npp
*     proposition d un deplacement d un point de u
         utemp(1,ipp) = max(u(1,ipp)-du/2.,xlim(1)) + ranf()*
     &        (min(u(1,ipp)+du/2.,xlim(2))-max(u(1,ipp)-du/2.,xlim(1)))
         utemp(2,ipp) = max(u(2,ipp)-du/2.,ylim(1)) + ranf()*
     &        (min(u(2,ipp)+du/2.,ylim(2))-max(u(2,ipp)-du/2.,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
         surf = (dx+du/2.)*(dy+du/2.)
         dx = min(du/2.,utemp(1,ipp)-xlim(1),xlim(2)-utemp(1,ipp))
         dy = min(du/2.,utemp(2,ipp)-ylim(1),ylim(2)-utemp(2,ipp))
         surftemp = (dx+du/2.)*(dy+du/2.)

*     modif de indcell et distcell
         call vormove(nindiv,nindivmax,s,npp,nppmax,utemp,
     &        indcell,distcell,indcelltemp,distcelltemp,ipp)

c         write(*,*) 'apres vormove'



         r = ratio(z,f,c,c,indcell,indcelltemp,
     &     npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'r=',r
         r = r*surf/surftemp

         alpha = amin1(1.,r)
c         write(*,*) 'alpha=',alpha
         bern = ignbin(1,alpha)
         if(bern .eq. 1) then
            u(1,ipp) = utemp(1,ipp)
            u(2,ipp) = utemp(2,ipp)
            do iindiv=1,nindiv
               indcell(iindiv) = indcelltemp(iindiv)
               distcell(iindiv) = distcelltemp(iindiv)
            enddo
         else 
            utemp(1,ipp) = u(1,ipp)
            utemp(2,ipp) = u(2,ipp)
         endif
      enddo
      end subroutine updurw




*
*     mise a jour de t    
* 
      subroutine updt(npp,nppmax,nindiv,
     &     nindivmax,nloc,nlocmax,nlocmax2,nallmax,npopmax,
     &     t,ttemp,dt,s,c,indcell,distcell,indcelltemp,distcelltemp,
     &     u,z,f,ploidy)
      implicit none 
      integer npp,nppmax,nindiv,nindivmax,nloc,nlocmax,nlocmax2,nallmax,
     &     npopmax,c(nppmax),indcell(nindivmax),z(nindivmax,nlocmax2),
     &     ploidy
      real t(2,nindivmax),s(2,nindivmax),distcell(nindivmax),
     &     u(2,nppmax),f(npopmax,nlocmax,nallmax),dt
      integer iindiv,ipp,accept,ignbin,indcelltemp(nindivmax)
      real ranf,d,ttemp(2,nindivmax),r,alpha,distcelltemp(nindivmax),
     &     ratio

*     initialisation
      do iindiv = 1,nindiv
         ttemp(1,iindiv) = t(1,iindiv)
         ttemp(2,iindiv) = t(2,iindiv)
         indcelltemp(iindiv) = indcell(iindiv)
         distcelltemp(iindiv) = distcell(iindiv)
      enddo
      if(nindivmax .gt. nindiv) then 
         do iindiv = nindiv+1,nindivmax
            ttemp(1,iindiv) = -999.
            ttemp(2,iindiv) = -999.
         enddo
      endif

      do iindiv = 1,nindiv
*     proposition d'une modif de t
         ttemp(1,iindiv) = s(1,iindiv) + dt*(ranf()-.5)
         ttemp(2,iindiv) = s(2,iindiv) + dt*(ranf()-.5)

*     modif de indcell et distcell
         distcelltemp(iindiv) = 3.e+37
         do ipp = 1,npp
            d = (ttemp(1,iindiv)-u(1,ipp))**2+
     &           (ttemp(2,iindiv)-u(2,ipp))**2
            if(d .lt. distcelltemp(iindiv)) then 
               indcelltemp(iindiv)  = ipp
               distcelltemp(iindiv) = d
            endif
         enddo

*     proba d'acceptation
         if(indcelltemp(iindiv) .ne. indcell(iindiv)) then 
            r = ratio(z,f,c,c,indcell,indcelltemp,
     &           npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,
     &           nlocmax2,nppmax,ploidy)
         else 
            r = 1.
         endif
         alpha = amin1(1.,r)
         accept = ignbin(1,alpha)
*     mise a jour en cas d'acceptation
         if(accept .eq. 1) then 
            indcell(iindiv) = indcelltemp(iindiv)
            distcell(iindiv) = distcelltemp(iindiv)
            t(1,iindiv) = ttemp(1,iindiv) 
            t(2,iindiv) = ttemp(2,iindiv)
         endif
      enddo
      end subroutine updt






*
*     naissance ou mort d'une cellule
*     avec prior Poisson(lambda) tronquée :   1 < m < nppmax
      subroutine bdpp(nindiv,nindivmax,u,c,utemp,ctemp,npop,npopmax,
     &     nloc,nlocmax,nlocmax2,nallmax,npp,nppmax,z,f,s,xlim,ylim,
     &     indcell,distcell,indcelltemp,distcelltemp,lambda,ploidy)
      implicit none 
      integer nindiv,nindivmax,nloc,nlocmax,nlocmax2,
     &     npop,npopmax,
     &     nallmax,npp,nppmax,z(nindivmax,nlocmax2),c(nppmax),
     &     indcell(nindivmax),ploidy
      real u(2,nindivmax),f(npopmax,nlocmax,nallmax),xlim(2),
     &     ylim(2),s(2,nindivmax),distcell(nindivmax),lambda

      integer ignbin,b,ctemp(nppmax),indcelltemp(nindivmax),ipp,npptemp,
     &     accept,iindiv,ipprem
      real utemp(2,nppmax),distcelltemp(nindivmax),ranf,
     &     ratio,r,alpha
      
*     naissance ou mort ?
      b = ignbin(1,0.5)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utemp (1,ipp) = u(1,ipp)
               utemp (2,ipp) = u(2,ipp)
               ctemp(ipp) = c(ipp)
            enddo
            npptemp = npp + 1
            ctemp(npptemp) = 1+ int(aint(float(npop)*ranf()))
            utemp(1,npptemp) = xlim(1)+(xlim(2)-xlim(1))*ranf()
            utemp(2,npptemp) = ylim(1)+(ylim(2)-ylim(1))*ranf()
            if(nppmax .gt. npptemp) then
               do ipp=npptemp+1,nppmax
                  ctemp(ipp) = -999
                  utemp(1,ipp) = -999.
                  utemp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,utemp,indcell,distcell,indcelltemp,
     &           distcelltemp,nindiv,nindivmax,npp,nppmax)
            r = ratio(z,f,c,ctemp,indcell,indcelltemp,npopmax,nlocmax,
     &           nallmax,nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
            r = r*lambda/float(npp+1)
            alpha = amin1(1.,r)
            accept = ignbin(1,alpha)
            if(accept .eq. 1) then 
               npp = npptemp
               do iindiv=1,nindivmax
                  indcell(iindiv) = indcelltemp(iindiv)
                  distcell(iindiv) = distcelltemp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utemp(1,ipp)
                  u (2,ipp) = utemp(2,ipp)
                  c(ipp) = ctemp(ipp)
               enddo
            endif
         endif
      else
*     mort
         if(npp .ne. 1) then 
            ipprem = 1+ aint(float(npp)*ranf())
            if(ipprem .ne. 1) then 
               do ipp = 1,ipprem-1
                  utemp (1,ipp) = u(1,ipp)
                  utemp (2,ipp) = u(2,ipp)
                  ctemp(ipp) = c(ipp)
               enddo
            endif
            if(ipprem .ne. npp) then 
               do ipp = ipprem,npp-1
                  utemp (1,ipp) = u(1,ipp+1)
                  utemp (2,ipp) = u(2,ipp+1)
                  ctemp(ipp) = c(ipp+1)
               enddo
            endif
            do ipp=npp,nppmax
               utemp (1,ipp) = -999.
               utemp (2,ipp) = -999.
               ctemp(ipp) = -999
            enddo

            call vorrem(s,utemp,ipprem,indcell,distcell,
     &           indcelltemp,distcelltemp,nindiv,nindivmax,npp,nppmax)

            r = ratio(z,f,c,ctemp,indcell,indcelltemp,npopmax,nlocmax,
     &           nallmax,nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
            r = r*float(npp)/lambda
            alpha = amin1(1.,r)
            accept = ignbin(1,alpha)
            if(accept .eq. 1) then 
               npp = npp-1
               do iindiv=1,nindivmax
                  indcell(iindiv) = indcelltemp(iindiv)
                  distcell(iindiv) = distcelltemp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utemp(1,ipp)
                  u (2,ipp) = utemp(2,ipp)
                  c(ipp) = ctemp(ipp)
               enddo
            endif
         endif
      endif
      end subroutine bdpp






*     calcul du ratio p(z|theta*)/p(z|theta)
*     ca ne depend pas de lambda
      real function ratio(z,f,c,ctemp,indcell,indcelltemp,
     &     npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
     &     indcell(nindivmax),indcelltemp(nindivmax),ploidy
      real f(npopmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,ipop,ipoptemp

c      write(*,*) 'debut de ratio'
c      write(*,*) 'indcell=',indcell
c      write(*,*) 'indcelltemp=',indcelltemp
c      write(*,*) 'c=',c
c      write(*,*) 'ctemp=',ctemp


      ratio = 1.
      do iindiv=1,nindiv
c         write(*,*) 'iindiv=', iindiv
         ipop = c(indcell(iindiv))
         ipoptemp = ctemp(indcelltemp(iindiv))
C         write(*,*) 'indcell=',indcell
C          write(*,*) 'indcelltemp=',indcelltemp
C          write(*,*) 'c=',c
C          write(*,*) 'ctemp=',ctemp
C          write(*,*) 'ipop=',ipop
C          write(*,*) 'ipoptemp=',ipoptemp

         do iloc=1,nloc
c             write(*,*) 'iloc=',iloc
c            write(6,*) 'z=',z(iindiv,2*iloc-1)
c            write(6,*) 'z=',z(iindiv,2*iloc)
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
c            ratio = ratio*
c     &           (f(ipoptemp,iloc,iall1)/f(ipop,iloc,iall1))*
c     &           (f(ipoptemp,iloc,iall2)/f(ipop,iloc,iall2))
            if(iall1 .ne. -999) then 
c               write(*,*) f(ipoptemp,iloc,iall1)
c               write(*,*) f(ipop,iloc,iall1)
               ratio = ratio*
     &              (f(ipoptemp,iloc,iall1)/f(ipop,iloc,iall1))
               
            endif
            if(iall2 .ne. -999) then 
c               write(*,*) f(ipoptemp,iloc,iall2)
c               write(*,*) f(ipop,iloc,iall2)
               ratio = ratio*
     &              (f(ipoptemp,iloc,iall2)/f(ipop,iloc,iall2))
            endif
         enddo
      enddo
      if(ploidy .eq. 1) then 
         ratio = sqrt(ratio)
      endif
c      write(*,*) 'fin de ratio'
      end function ratio



c$$$
c$$$*     calcul du ratio p(z|theta*)/p(z|theta)
c$$$*     quand f est  modifié
c$$$      real function ratiobd(z,f,ftemp,c,ctemp,indcell,indcelltemp,
c$$$     &     npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
c$$$     &     nppmax)
c$$$      implicit none
c$$$      integer npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
c$$$     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
c$$$     &     indcell(nindivmax),indcelltemp(nindivmax)
c$$$      real f(npopmax,nlocmax,nallmax),ftemp(npopmax,nlocmax,nallmax)
c$$$      integer iindiv,iloc,iall1,iall2,ipop,ipoptemp
c$$$
c$$$      ratiobd = 1.
c$$$      do iindiv=1,nindiv
c$$$         ipop = c(indcell(iindiv))
c$$$         ipoptemp = ctemp(indcelltemp(iindiv))
c$$$         do iloc=1,nloc
c$$$            iall1 = z(iindiv,2*iloc-1)
c$$$            iall2 = z(iindiv,2*iloc)
c$$$            if(iall1 .ne. -999) then 
c$$$               ratiobd = ratiobd*
c$$$     &              (ftemp(ipoptemp,iloc,iall1)/f(ipop,iloc,iall1))
c$$$            endif
c$$$            if(iall2 .ne. -999) then 
c$$$               ratiobd = ratiobd*
c$$$     &              (ftemp(ipoptemp,iloc,iall2)/f(ipop,iloc,iall2))
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$c      write(*,*) 'f=',f
c$$$c      write(*,*) 'ftemp=',ftemp
c$$$c      write(*,*) 'c=',c
c$$$c      write(*,*) 'ctemp=',ctemp
c$$$c      write(*,*) 'ratiobd=',ratiobd
c$$$      end function ratiobd      C





*
*     Indice des cellules dans une pope
*
      subroutine who(c,ipop,npp,nppmax,cellpop,
     &     ncellpop)
      implicit none
      integer npp,nppmax,c(nppmax),ipop,cellpop(nppmax),
     &     ncellpop
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
      real ranf
c      write(*,*) 'sample'

*     init
      ii = 1 + int(aint(float(ncellpop)*ranf()))
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
            ii = 1 + int(aint(float(ncellpop-isamp)*ranf()))
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
      real ranf
c      call rndstart()
c      write(*,*) 'sample2'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
*     init
      ii = 1 + int(aint(float(ncellpop)*ranf()))
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
     &           int(aint(float(ncellpop-isamp)*ranf()))
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
*     split d'une pope en deux
*     reallocation de nu cellules dont les indices
*     sont dans listcell
*     dans la pope ipop
*
      subroutine split(ipop,c,ctemp,nppmax,nu,listcell)
      implicit none
      integer ipop,nppmax,c(nppmax),ctemp(nppmax),nu,
     &     listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de split'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipop=',ipop
      do ipp=1,nppmax
         ctemp(ipp) = c(ipp)
      enddo
      if(nu .gt. 0) then
         do ii=1,nu
            ctemp(listcell(ii)) = ipop
         enddo
      endif
c      write(*,*)'c=',c
c      write(*,*)'ctemp=',ctemp
c      write(*,*) 'fin de split'
      end subroutine split


*     merge de deux  popes en une : 
*     reallocation des nu cellules de la pope ipoprem 
*     dont les indices sont dans listcell
*     dans la pope ipophost
*
      subroutine merging(ipoprem,ipophost,
     &     c,ctemp,nppmax,nu,listcell)
      implicit none
      integer ipoprem,ipophost,nppmax,
     &     c(nppmax),ctemp(nppmax),nu,listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de merge'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipoprem=',ipoprem
      do ipp=1,nppmax
         ctemp(ipp) = c(ipp)
      enddo
      if(ipoprem .gt. ipophost) then 
         if(nu .gt. 0) then
            do ii=1,nu
               ctemp(listcell(ii)) = ipophost
            enddo
         endif
      else
         if(nu .gt. 0) then
            do ii=1,nu
               ctemp(listcell(ii)) = ipophost - 1
            enddo
         endif
      endif
      do ipp=1,nppmax
         if(c(ipp) .gt. ipoprem) ctemp(ipp) = c(ipp)-1
      enddo
c      write(*,*)'c=',c
c      write(*,*)'ctemp=',ctemp
c      write(*,*) 'fin de merge'
      end subroutine merging

 




*
*     Mise a jour de c et f en cas d acceptation d'un split/merge
*
      subroutine accept5(nppmax,npopmax,nlocmax,nallmax,
     &     nall,c,ctemp,f,ftemp,drift,drifttemp)
      implicit none
      integer nppmax,npopmax,nlocmax,nallmax,
     &     nall(nlocmax),c(nppmax),ctemp(nppmax)
      real f(npopmax,nlocmax,nallmax),
     &     ftemp(npopmax,nlocmax,nallmax),
     &     drift(npopmax),drifttemp(npopmax)
      integer ipop,iloc,iall,ipp
c      write(*,*) 'debut de accept5'
c      write(*,*) 'f=',f
c      write(*,*) 'ftemp=',ftemp
      do ipp=1,nppmax
         c(ipp) = ctemp(ipp)
      enddo
      do ipop = 1,npopmax
         do iloc= 1,nlocmax
            do iall=1,nall(iloc)
               f(ipop,iloc,iall) = ftemp(ipop,iloc,iall)
            enddo
         enddo
         drift(ipop) = drifttemp(ipop)
      enddo

c      write(*,*) 'f=',f
c      write(*,*) 'fin de accept5'
      end subroutine accept5

*
*     coefficients du binome C_n^p
*
      real function bico(n,p)
      implicit none
      integer n,p
      real algama
      bico = exp(algama(float(n+1))-algama(float(p+1))-
     &     algama(float(n-p+1)))
c      write(*,*) 'in bico '
c$$$      write(*,*) 'n=', n
c$$$      write(*,*) 'p=', p
c$$$      write(*,*) 'algama(float(n+1))=',algama(float(n+1))
c$$$      write(*,*) 'algama(float(p+1))=',algama(float(p+1))
c$$$      write(*,*) 'algama(float(n-p+1)))=',algama(float(n-p+1))
c$$$      write(*,*) 'exp()=',exp(algama(float(n+1))-algama(float(p+1))-
c$$$     &     algama(float(n-p+1)))
c$$$      write(*,*) 'bico =', nint(exp(algama(float(n+1))-
c$$$     &     algama(float(p+1))-
c$$$     &     algama(float(n-p+1))))
c      write(*,*) 'bico =',bico
      
      end function bico



*****************************************************************
*     ln du coefficient du binome C_n^p
*
      real function lbico(n,p)
      implicit none
      integer n,p
      real algama
      lbico = algama(float(n+1))-algama(float(p+1))-
     &     algama(float(n-p+1))
      end function lbico




*
*     log du ratio des vraisemblances dans bdpop6
*
      real function llr6(z,f,ftemp,c,ctemp,indcell,indcelltemp,
     &     npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
     &     indcell(nindivmax),indcelltemp(nindivmax),ploidy
      real f(npopmax,nlocmax,nallmax),ftemp(npopmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,ipop,ipoptemp

      llr6 = 0

*     log du rapport des vraisemblances
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptemp = ctemp(indcelltemp(iindiv))
         do iloc=1,nloc
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
               llr6 = llr6 + 
     &              alog(ftemp(ipoptemp,iloc,iall1)) - 
     &              alog(f(ipop,iloc,iall1))
            endif
            if(iall2 .ne. -999) then 
               llr6 = llr6 + 
     &              alog(ftemp(ipoptemp,iloc,iall2)) - 
     &              alog(f(ipop,iloc,iall2))
            endif
            
c$$$            if(ipoptemp .eq. 3) then 
c$$$               write(*,*) 'ipoptemp=',ipoptemp
c$$$               write(*,*) 'iinidiv=',iindiv
c$$$               write(*,*) 'c=', c(indcell(iindiv))
c$$$               write(*,*) 'ctemp=', ctemp(indcelltemp(iindiv))
c$$$               write(*,*) 'z=',z(iindiv,1),z(iindiv,2)
c$$$               write(*,*) 'llr6 =',llr6
c$$$               write(*,*) 'ftemp(ipoptemp,iloc,iall1)=',
c$$$     &              ftemp(ipoptemp,iloc,iall1) 
c$$$               write(*,*) 'f(ipop,iloc,iall1)=',
c$$$     &              f(ipop,iloc,iall1)
c$$$               write(*,*) 'ftemp(ipoptemp,iloc,iall2)=',
c$$$     &              ftemp(ipoptemp,iloc,iall2)
c$$$               write(*,*) 'f(ipop,iloc,iall2)=',
c$$$     &              f(ipop,iloc,iall2)
c$$$            endif

         enddo
      enddo
      if(ploidy .eq. 1) llr6 = 0.5*llr6 
      end function llr6


*
*     comptage des alleles dans chaque pop pour c
*
      subroutine countn(nindiv,nlocmax,nlocmax2,npopmax,
     &     nppmax,nall,nallmax,z,n,indcell,c,ploidy)
      implicit none
      integer  nindiv,nlocmax,nlocmax2,npopmax,nppmax,nallmax,ploidy,
     &     z(nindiv,nlocmax2),nall(nlocmax),
     &     n(npopmax,nlocmax,nallmax),c(nppmax),indcell(nindiv)
      integer ipop,iloc,iall,iindiv
*     init du tableau
      do ipop = 1,npopmax
         do iloc = 1,nlocmax
            do iall =1,nall(iloc)
               n(ipop,iloc,iall)=0
            enddo
         enddo
      enddo
c      write(*,*) 'dans countn'
c      write(*,*) 'c=',c
c      write(*,*) 'indcell=',indcell
c      write(*,*) 'n=',n
*     comptage
      do iindiv = 1,nindiv
         do iloc = 1,nlocmax
            if(z(iindiv,2*iloc-1) .ne. -999) then
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1))+ 1 
c               write(*,*) 'n=',
c     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) 
            endif
            if(ploidy .eq. 2) then
               if(z(iindiv,2*iloc) .ne. -999) then 
                  n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) = 
     &                 n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) + 1 
c     write(*,*) 'n=',
c     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) 
               endif
            endif
         enddo
      enddo
c      write(*,*) 'n=',n
      end subroutine countn



*
*     log du ratio (prob cond. complete)/prior
*     pour les frequences
*     dans un split de la pope ipop
      real function lrf(ipop,npopmax,nlocmax,nall,nallmax,
     &     f,fa,drift,n)
      implicit none
      integer ipop,npopmax,nlocmax,nall(nlocmax),nallmax,
     &     n(npopmax,nlocmax,nallmax)
      real f(npopmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),
     &     drift(npopmax)
      integer iloc,iall,nn
      real ss,algama,q

      lrf = 0.
      q = (1-drift(ipop))/drift(ipop)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q +
     &           float(n(ipop,iloc,iall))) +
     &           (1 - fa(iloc,iall) * q - float(n(ipop,iloc,iall)))*
     &           alog(f(ipop,iloc,iall))
            nn = nn + n(ipop,iloc,iall)
         enddo
c         write(*,*) 'nn=',nn
         lrf = lrf + algama(float(nall(iloc))) -
     &        algama(q + nn) + ss

      enddo
      end function lrf



      
*
*     Naissance et mort de popes avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdpop7(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &     cellpophost,n,ntemp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer b,ignbin,
     &     ipoprem,ipp,bern,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      real alpha,ranf,lbico,lratio,llr6,termfsplit,termfmerge
      
       do ipp=1,nppmax
          cellpop(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ignbin(1,0.5)
       
       if(b .eq. 1) then
          if(npop .lt. npopmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la pope qui split
             isplit = 1 + int(aint(float(npop)*ranf()))

*     recherche des cellules affectees a cette pope
             call who(c,isplit,npp,nppmax,cellpop,ncellpop)
             if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
                nu = int(aint(float(ncellpop)*ranf()))
                if(nu .gt. 0) then

*     tirage des cellules reallouees
                   call sample(cellpop,nppmax,nu,ncellpop,listcell)

*     proposition de reallocation dans la pope npop+1
                   call split(npop+1,c,ctemp,nppmax,nu,
     &                  listcell)

*     comptage des alleles sur chaque locus pour c puis ctemp
                   call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &                  nppmax,nall,nallmax,z,n,indcell,c,ploidy)
                   call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &                  nppmax,nall,nallmax,z,ntemp,indcell,ctemp,
     &                  ploidy)

*     proposition nouvelle freq et derive 
c                    call addfreq5(isplit,npop,npopmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
                   call addfreq7(npop,npopmax,nloc,nlocmax,
     &     nall,nallmax,isplit,f,ftemp,
     &     fa,drift,drifttemp,a,ptemp,ntemp)

*     calcul du log du ratio
*     terme des vraisemblances
                   lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,npopmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
                   lratio = lratio + termfsplit(isplit,npop,npopmax,
     &                  nlocmax,nall,nallmax,
     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)

*     terme des proposal sur c
                   lratio = lratio + alog(2*float(ncellpop+1)) + 
     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 

*     terme des priors sur c
                   lratio = lratio + 
     &                  float(npp)*(alog(float(npop)) - 
     &                  alog(float(npop+1)))

                   lratio = amin1(0.,lratio)
                   alpha = exp(lratio)
                   bern = ignbin(1,alpha)

c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                           n(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,npop,npopmax,
c$$$     &                  nlocmax,nall,nallmax,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  alog(2*float(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   float(npp)*(alog(float(npop)) - 
c$$$     &                  alog(float(npop+1)))
c$$$                   write(*,*) 'alpha=',alpha 


                   if(bern .eq. 1) then
                   call accept5(nppmax,npopmax,nlocmax,
     &                  nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                   npop = npop + 1
                endif
               endif
            endif 
         endif

*     merge
      else
         if(npop .gt. npopmin) then 
c             write(*,*) 'mort'
*     tirage de la pope qui meurt
            ipoprem = 1 + int(aint(float(npop)*ranf()))
            
*     tirage de la pope hote
            ipophost = 1 + int(aint(float(npop)*ranf()))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + int(aint(float(npop)*ranf()))
            enddo

*     on range dans la pope d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pope hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
               
            if(ncellpop .gt. 0) then
*     proposition de reallocation dans la pope ipophost
               call merging(ipoprem,ipophost,c,ctemp,nppmax,
     &              ncellpop,cellpop)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(ipoprem,ipophost,npop,npopmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,npopmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmerge(ipophost,ipoprem,
     &              npopmax,nlocmax,
     &              nall,nallmax,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp)

*     terme des proposal sur c
               lratio = lratio + alog(float(npop)) - 
     &              alog(2*float(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
               lratio = lratio + 
     &              float(npp)*(alog(float(npop)) - 
     &              alog(float(npop-1)))
               lratio = amin1(0.,lratio)
               alpha  = exp(lratio)
               bern = ignbin(1,alpha)      
      
         
c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                       n(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npop,npopmax,nlocmax,
c$$$     &              nall,nallmax,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$               write(*,*) ' term prop c=',
c$$$     &              alog(float(npop)) - 
c$$$     &              alog(2*float(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              float(npp)*(alog(float(npop)) - 
c$$$     &              alog(float(npop-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,npopmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  npop = npop - 1
               endif
            endif
         endif
      endif
      end subroutine bdpop7

*
*     Naissance et mort de popes avec réallocations 
*     (split/merge)
*     drift* reste à 0.5
*     pour court-circuiter le F-model
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
* 
      subroutine bdpop7bis(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &     cellpophost,n,ntemp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer b,ignbin,
     &     ipoprem,ipp,bern,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      real alpha,ranf,lbico,lratio,llr6,
     &     termfsplitbis,termfmergebis
      
       do ipp=1,nppmax
          cellpop(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ignbin(1,0.5)
       
       if(b .eq. 1) then
          if(npop .lt. npopmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la pope qui split
             isplit = 1 + int(aint(float(npop)*ranf()))

*     recherche des cellules affectees a cette pope
             call who(c,isplit,npp,nppmax,cellpop,ncellpop)
             if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
                nu = int(aint(float(ncellpop)*ranf()))
                if(nu .gt. 0) then

*     tirage des cellules reallouees
                   call sample(cellpop,nppmax,nu,ncellpop,listcell)

*     proposition de reallocation dans la pope npop+1
                   call split(npop+1,c,ctemp,nppmax,nu,
     &                  listcell)

*     comptage des alleles sur chaque locus pour c puis ctemp
                   call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &                  nppmax,nall,nallmax,z,n,indcell,c,ploidy)
                   call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &                  nppmax,nall,nallmax,z,ntemp,indcell,ctemp,
     &                  ploidy)

*     proposition nouvelle freq et derive 
c                    call addfreq5(isplit,npop,npopmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
                   call addfreq7bis(npop,npopmax,nloc,
     &     nlocmax,nall,nallmax,isplit,
     &     f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)

*     calcul du log du ratio
*     terme des vraisemblances
                   lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,npopmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
                   lratio = lratio + termfsplitbis(isplit,npop,
     &                  npopmax,nlocmax,nall,nallmax,
     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)

*     terme des proposal sur c
                   lratio = lratio + alog(2*float(ncellpop+1)) + 
     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 

*     terme des priors sur c
                   lratio = lratio + 
     &                  float(npp)*(alog(float(npop)) - 
     &                  alog(float(npop+1)))

                   lratio = amin1(0.,lratio)
                   alpha = exp(lratio)
                   bern = ignbin(1,alpha)

c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                           n(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,npop,npopmax,
c$$$     &                  nlocmax,nall,nallmax,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  alog(2*float(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   float(npp)*(alog(float(npop)) - 
c$$$     &                  alog(float(npop+1)))
c$$$                   write(*,*) 'alpha=',alpha 


                   if(bern .eq. 1) then
                   call accept5(nppmax,npopmax,nlocmax,
     &                  nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                   npop = npop + 1
                endif
               endif
            endif 
         endif

*     merge
      else
         if(npop .gt. npopmin) then 
c             write(*,*) 'mort'
*     tirage de la pope qui meurt
            ipoprem = 1 + int(aint(float(npop)*ranf()))
            
*     tirage de la pope hote
            ipophost = 1 + int(aint(float(npop)*ranf()))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + int(aint(float(npop)*ranf()))
            enddo

*     on range dans la pope d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pope hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
               
            if(ncellpop .gt. 0) then
*     proposition de reallocation dans la pope ipophost
               call merging(ipoprem,ipophost,c,ctemp,nppmax,
     &              ncellpop,cellpop)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(ipoprem,ipophost,npop,npopmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7bis(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,npopmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmergebis(ipophost,ipoprem,
     &              npopmax,nlocmax,
     &              nall,nallmax,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp)

*     terme des proposal sur c
               lratio = lratio + alog(float(npop)) - 
     &              alog(2*float(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
               lratio = lratio + 
     &              float(npp)*(alog(float(npop)) - 
     &              alog(float(npop-1)))
               lratio = amin1(0.,lratio)
               alpha  = exp(lratio)
               bern = ignbin(1,alpha)       
      
         
c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                       n(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,nlocmax2,
c$$$     &              nall,nallmax,nindiv,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$               write(*,*) ' term prop c=',
c$$$     &              alog(float(npop)) - 
c$$$     &              alog(2*float(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              float(npp)*(alog(float(npop)) - 
c$$$     &              alog(float(npop-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,npopmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  npop = npop - 1
               endif
            endif
         endif
      endif
      end subroutine bdpop7bis







***********************************************************************
*
*     Naissance et mort de popes avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdpop8(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &     cellpophost,n,ntemp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      real alpha,ranf,lbico,lratio,llr6,termfsplit,
     &     termfmerge
      integer b,ignbin,bern
      
       do ipp=1,nppmax
          cellpop(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ignbin(1,0.5e0)
       
       if(b .eq. 1) then
          if(npop .lt. npopmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la pope qui split
             isplit = 1 + int(aint(float(npop)*ranf()))
             
*     recherche des cellules affectees a cette pope
             call who(c,isplit,npp,nppmax,cellpop,ncellpop)
             
             if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
                nu = int(aint(float(ncellpop)*ranf()))
                if(nu .gt. 0) then
                   
*     tirage des cellules reallouees
                   call sample2(cellpop,nppmax,nu,ncellpop,
     &                  listcell)
                   
*     proposition de reallocation dans la pope npop+1
                   call split(npop+1,c,ctemp,nppmax,nu,
     &                  listcell)
                else 
                   do ipp = 1,nppmax
                      ctemp(ipp) = c(ipp)
                   enddo
                endif
             else
                nu = 0
                do ipp = 1,nppmax
                   ctemp(ipp) = c(ipp)
                enddo
             endif

*     comptage des alleles sur chaque locus pour c puis ctemp
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nall,nallmax,z,n,indcell,c,ploidy)
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)
             
*     proposition nouvelle freq et derive 
c     call addfreq5(isplit,npop,npopmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
             call addfreq7(npop,npopmax,nloc,nlocmax,
     &            nall,nallmax,isplit,
     &            f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)
             
*     calcul du log du ratio
*     terme des vraisemblances
             lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &            indcell,npopmax,nlocmax,nallmax,
     &            nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
             lratio = lratio + termfsplit(isplit,npop,npopmax,
     &            nlocmax,nall,nallmax,
     &            f,ftemp,n,ntemp,fa,drift,drifttemp)
             
*     terme des proposal sur c
             lratio = lratio + alog(2*float(ncellpop+1)) + 
     &            lbico(ncellpop,nu) - alog(float(npop+1)) 
             
*     terme des priors sur c
             lratio = lratio + 
     &            float(npp)*(alog(float(npop)) - 
     &            alog(float(npop+1)))
             
             lratio = amin1(0.e0,lratio)
             alpha = exp(lratio)
             bern = ignbin(1,alpha)

c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                           n(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,npop,npopmax,
c$$$     &                  nlocmax,nall,nallmax,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  alog(2*float(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   float(npp)*(alog(float(npop)) - 
c$$$     &                  alog(float(npop+1)))
c$$$                   write(*,*) 'alpha=',alpha 

             
             if(bern .eq. 1) then
                call accept5(nppmax,npopmax,nlocmax,
     &               nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                npop = npop + 1
             endif
          endif

*     merge
      else
         if(npop .gt. npopmin) then 
c             write(*,*) 'mort'
*     tirage de la pope qui meurt
            ipoprem = 1 + int(aint(float(npop)*ranf()))
            
*     tirage de la pope hote
            ipophost = 1 + int(aint(float(npop)*ranf()))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + 
     &              int(aint(float(npop)*ranf()))
            enddo

*     on range dans la pope d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pope hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
               
            if(ncellpop .gt. 0) then
*     proposition de reallocation dans la pope ipophost
               call merging(ipoprem,ipophost,c,ctemp,nppmax,
     &              ncellpop,cellpop)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(ipoprem,ipophost,npop,npopmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,npopmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmerge(ipophost,ipoprem,
     &              npopmax,nlocmax,
     &              nall,nallmax,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp)

*     terme des proposal sur c
               lratio = lratio + alog(float(npop)) - 
     &              alog(2*float(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
               lratio = lratio + 
     &              float(npp)*(alog(float(npop)) - 
     &              alog(float(npop-1)))
               lratio = amin1(0.e0,lratio)
               alpha  = exp(lratio)
               bern = ignbin(1,alpha)      
      
         
c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                       n(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nall,nallmax,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$               write(*,*) ' term prop c=',
c$$$     &              alog(float(npop)) - 
c$$$     &              alog(2*float(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              float(npp)*(alog(float(npop)) - 
c$$$     &              alog(float(npop-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,npopmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  npop = npop - 1
               endif
            endif
         endif
      endif
      end subroutine bdpop8




***********************************************************************
*     split/merge populations in the spatial D-model
*     changes from bdpop7bis:
*     - process populations whatever the number of tiles or individuals
*       they have
      subroutine bdpop8bis(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellpop,listcell,
     &     cellpophost,n,ntemp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(npopmax,nloc,nallmax),ntemp(npopmax,nloc,nallmax),
     &     ploidy
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &      ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      real alpha,ranf,lbico,lratio,llr6,
     &     termfsplitbis,termfmergebis
      integer b,ignbin,bern


      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ignbin(1,0.5)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pope qui split
            isplit = 1 + int(aint(float(npop)*ranf()))
c            write(*,*) 'isplit=',isplit

*     recherche des cellules affectees a cette pope
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)

            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees 
               nu = int(aint(float(ncellpop+1)*ranf()))
               if(nu .gt. 0) then                 

*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)

 
*     proposition de reallocation dans la pope npop+1
                  call split(npop+1,c,ctemp,nppmax,nu,
     &                 listcell)
c                  write(*,*) 'apres split'
c                  write(*,*) 'z(67,9)=',z(67,9)
c                  write(*,*) 'z(67,10)=',z(67,10)
               else 
                  do ipp = 1,nppmax
                     ctemp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctemp(ipp) = c(ipp)
               enddo
            endif
            
*     comptage des alleles sur chaque locus pour c puis ctemp
c            write(*,*) 'comptage des alleles'
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nall,nallmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)
c            write(*,*) 'apres count'
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)
            
*     proposition nouvelle freq et derive 
c            write(*,*) 'ajoutage des freq'
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)
            call addfreq7bis(npop,npopmax,nloc,
     &           nlocmax,nall,nallmax,isplit,
     &           f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)

*     calcul du log du ratio
*     terme des vraisemblances
c            write(*,*) 'calcul du log du ratio'
            lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &           indcell,npopmax,nlocmax,nallmax,
     &           nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
c            write(*,*) 'lratio =',lratio

*     terme des freq.
            lratio = lratio + termfsplitbis(isplit,npop,
     &           npopmax,nlocmax,nall,nallmax,
     &           f,ftemp,n,ntemp,fa,drift,drifttemp)
c            write(*,*) 'lratio =',lratio
*     terme des proposal sur c
            lratio = lratio + alog(2*float(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - alog(float(npop+1)) 
c            write(*,*) 'lratio =',lratio
*     terme des priors sur c
            lratio = lratio + 
     &           float(npp)*(alog(float(npop)) - 
     &           alog(float(npop+1)))
c*     Poisson prior on npop
c     &           -alog(float(npop+1)
c            write(*,*) 'lratio =',lratio

            lratio = amin1(0.e0,lratio)
            alpha = exp(lratio)
            bern = ignbin(1,alpha)
c$$$
c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                           n(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplitbis(isplit,npop,npopmax,
c$$$     &                  nlocmax,nall,nallmax,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  alog(2*float(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - alog(float(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   float(npp)*(alog(float(npop)) - 
c$$$     &                  alog(float(npop+1)))
c$$$                   write(*,*) 'alpha=',alpha 
            
            if(bern .eq. 1) then
c               write(*,*) 'accept split'
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)

               call accept5(nppmax,npopmax,nlocmax,
     &              nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
               npop = npop + 1
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
            endif
         endif

      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
*     tirage de la pope qui meurt
            ipoprem = 1 + int(aint(float(npop)*ranf()))
            
*     tirage de la pope hote
            ipophost = 1 + int(aint(float(npop)*ranf()))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 
     &              + int(aint(float(npop)*ranf()))
            enddo
            
*     on range dans la pope d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pope hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)

*     proposition de reallocation dans la pope ipophost
            call merging(ipoprem,ipophost,c,ctemp,nppmax,
     &           ncellpop,cellpop)
            
*     comptage des alleles sur chaque locus pour c puis ctemp
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nall,nallmax,z,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nall,nallmax,z,ntemp,indcell,ctemp,ploidy)

*     propostion du nouveau tableau de freq et de derives
            call remfreq7bis(ipoprem,ipophost,
     &           npop,npopmax,nloc,nlocmax,nall,
     &           nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,
     &           ntemp)
            
*     calcul du log du ratio  
*     terme des vraisemblances
            lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &           indcell,npopmax,nlocmax,nallmax,
     &           nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
c            write(*,*) 'lratio =',lratio
*     terme des freq.
            lratio = lratio + 
     &           termfmergebis(ipophost,ipoprem,
     &           npopmax,nlocmax,
     &           nall,nallmax,
     &           f,ftemp,n,ntemp,fa,drift,drifttemp)
c            write(*,*) 'lratio =',lratio
*     terme des proposal sur c
            lratio = lratio + alog(float(npop)) - 
     &           alog(2*float(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
c            write(*,*) 'lratio =',lratio
*     terme des priors sur c
            lratio = lratio + 
     &           float(npp)*(alog(float(npop)) - 
     &           alog(float(npop-1)))
c*     Poisson prior on npop
c     &           + alog(float(npop))
c            write(*,*) 'lratio =',lratio

            lratio = amin1(0.e0,lratio)
            alpha  = exp(lratio)
            bern = ignbin(1,alpha) 

c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',iall,')=',
c$$$     &                       n(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ntemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',iall,')=',
c$$$     &                       f(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',ipop,',',iloc,',',iall,')=',
c$$$     &                           ftemp(ipop,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmergebis(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nall,nallmax,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp)
c$$$               write(*,*) ' term prop c=',
c$$$     &              alog(float(npop)) - 
c$$$     &              alog(2*float(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              float(npp)*(alog(float(npop)) - 
c$$$     &              alog(float(npop-1)))
c$$$               write(*,*) 'alpha=',alpha 
  
            if(bern .eq. 1) then
c               write(*,*) 'accept merge'
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
               call accept5(nppmax,npopmax,nlocmax,
     &              nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
               npop = npop - 1
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
            endif
         endif
      endif 
      
      end subroutine bdpop8bis

 




*
*     ajoute une pope 
*     dans  le tableau des dérives selon le prior
*     et dans le tableau des frequences 
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq7(npop,npopmax,nloc,nlocmax,
     &     nall,nallmax,isplit,f,ftemp,fa,
     &     drift,drifttemp,a,ptemp,ntemp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,ntemp(npopmax,nloc,nallmax),
     &     isplit
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax),a(nallmax)
      integer iloc,iall,ipop
      real ptemp(nallmax),ranf


*     remplissage de f et drift pour le popes pre-existantes
      do ipop = 1,npop
         drifttemp(ipop) = drift(ipop)
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttemp(isplit) = ranf()
      drifttemp(npop+1) = ranf()

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(isplit))/
     &           drifttemp(isplit)+float(ntemp(isplit,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(isplit,iloc,iall)  = ptemp(iall)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(npop+1))/
     &           drifttemp(npop+1)+float(ntemp(npop+1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(npop+1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftemp(',isplit,1,1,')=',ftemp(isplit,1,1)
c$$$      write(*,*) 'ftemp(',isplit,1,2,')=',ftemp(isplit,1,2)
c$$$      write(*,*) 'f(',npop+1,1,1,')=',f(npop+1,1,1)
c$$$      write(*,*) 'f(',npop+1,1,2,')=',f(npop+1,1,2)
c$$$      write(*,*) 'ftemp(',npop+1,1,1,')=',ftemp(npop+1,1,1)
c$$$      write(*,*) 'ftemp(',npop+1,1,2,')=',ftemp(npop+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftemp(',isplit,2,1,')=',ftemp(isplit,2,1)
c$$$      write(*,*) 'ftemp(',isplit,2,2,')=',ftemp(isplit,2,2)
c$$$      write(*,*) 'f(',npop+1,2,1,')=',f(npop+1,2,1)
c$$$      write(*,*) 'f(',npop+1,2,2,')=',f(npop+1,2,2)
c$$$      write(*,*) 'ftemp(',npop+1,2,1,')=',ftemp(npop+1,2,1)
c$$$      write(*,*) 'ftemp(',npop+1,2,2,')=',ftemp(npop+1,2,2)
c$$$

      end subroutine addfreq7

*
*     ajoute une pope 
*     dans  le tableau des frequences  selon le prior
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     et une valeur 0.5 dans le tableau des dérives 
*     (sans modifier les tableaux en entrée)
*     pour court-circuiter le F-model
      subroutine addfreq7bis(npop,npopmax,nloc,nlocmax,
     &     nall,nallmax,isplit,f,ftemp,
     &     fa,drift,drifttemp,a,ptemp,ntemp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,ntemp(npopmax,nloc,nallmax),
     &     isplit
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax),a(nallmax)
      integer iloc,iall,ipop
      real ptemp(nallmax)


*     remplissage de f et drift pour le popes pre-existantes
      do ipop = 1,npop
         drifttemp(ipop) = drift(ipop)
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttemp(isplit) = 0.5
      drifttemp(npop+1) = 0.5

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(isplit))/
     &           drifttemp(isplit)+float(ntemp(isplit,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(isplit,iloc,iall)  = ptemp(iall)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(npop+1))/
     &           drifttemp(npop+1)+float(ntemp(npop+1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(npop+1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftemp(',isplit,1,1,')=',ftemp(isplit,1,1)
c$$$      write(*,*) 'ftemp(',isplit,1,2,')=',ftemp(isplit,1,2)
c$$$      write(*,*) 'f(',npop+1,1,1,')=',f(npop+1,1,1)
c$$$      write(*,*) 'f(',npop+1,1,2,')=',f(npop+1,1,2)
c$$$      write(*,*) 'ftemp(',npop+1,1,1,')=',ftemp(npop+1,1,1)
c$$$      write(*,*) 'ftemp(',npop+1,1,2,')=',ftemp(npop+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftemp(',isplit,2,1,')=',ftemp(isplit,2,1)
c$$$      write(*,*) 'ftemp(',isplit,2,2,')=',ftemp(isplit,2,2)
c$$$      write(*,*) 'f(',npop+1,2,1,')=',f(npop+1,2,1)
c$$$      write(*,*) 'f(',npop+1,2,2,')=',f(npop+1,2,2)
c$$$      write(*,*) 'ftemp(',npop+1,2,1,')=',ftemp(npop+1,2,1)
c$$$      write(*,*) 'ftemp(',npop+1,2,2,')=',ftemp(npop+1,2,2)
c$$$

      end subroutine addfreq7bis



*     enleve une pope des tableau des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux popes
*     tirage d'une derive selon prior 
*     sans modifier des tableaux en entrée
      subroutine remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,ntemp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,nindiv,
     &     ntemp(npopmax,nlocmax,nallmax)
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax),a(nallmax) 
      integer ipop,iloc,iall
      real ptemp(nallmax),ranf

*     supprime la pope ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop-1,iloc,iall) = f(ipop,iloc,iall)
               enddo 
            enddo
            drifttemp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
               enddo 
            enddo
            drifttemp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do iall=1,nall(iloc)
                     ftemp(ipop-1,iloc,iall) = f(ipop,iloc,iall)
                  enddo 
               enddo
               drifttemp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = -999
            enddo 
         enddo
         drifttemp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pope
      drifttemp(ipophost) = ranf()
      
*     frequences pour la nouvelle pope
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(ipophost))/
     &           drifttemp(ipophost)+
     &           float(ntemp(ipophost,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(ipophost,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftemp(',ipophost,1,1,')=',ftemp(ipophost,1,1)
c$$$      write(*,*) 'ftemp(',ipophost,1,2,')=',ftemp(ipophost,1,2)
c$$$      write(*,*) 'ftemp(',ipophost,2,1,')=',ftemp(ipophost,2,1)
c$$$      write(*,*) 'ftemp(',ipophost,2,2,')=',ftemp(ipophost,2,2)

      end subroutine remfreq7

*     enleve une pope des tableaux des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux popes
*     la nouvelle derive est mise à 0.5
*     (sans modifier les tableaux en entrée)
*     c'est pour court-circuiter le F-model 
      subroutine remfreq7bis(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,f,ftemp,drift,drifttemp,fa,a,ptemp,ntemp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,nindiv,
     &     ntemp(npopmax,nlocmax,nallmax)
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax),a(nallmax) 
      integer ipop,iloc,iall
      real ptemp(nallmax)

*     supprime la pope ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop-1,iloc,iall) = f(ipop,iloc,iall)
               enddo 
            enddo
            drifttemp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
               enddo 
            enddo
            drifttemp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do iall=1,nall(iloc)
                     ftemp(ipop-1,iloc,iall) = f(ipop,iloc,iall)
                  enddo 
               enddo
               drifttemp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = -999
            enddo 
         enddo
         drifttemp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pope
      drifttemp(ipophost) = 0.5
      
*     frequences pour la nouvelle pope
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(ipophost))/
     &           drifttemp(ipophost)+
     &           float(ntemp(ipophost,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(ipophost,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftemp(',ipophost,1,1,')=',ftemp(ipophost,1,1)
c$$$      write(*,*) 'ftemp(',ipophost,1,2,')=',ftemp(ipophost,1,2)
c$$$      write(*,*) 'ftemp(',ipophost,2,1,')=',ftemp(ipophost,2,1)
c$$$      write(*,*) 'ftemp(',ipophost,2,2,')=',ftemp(ipophost,2,2)

      end subroutine remfreq7bis



*
*     terme des freq dans le log ratio pour un split 
      real function termfsplit(isplit,npop,npopmax,
     &     nlocmax,nall,nallmax,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp)
      implicit none
      integer npopmax,nlocmax,nallmax,
     &     n(npopmax,nlocmax,nallmax),
     &     ntemp(npopmax,nlocmax,nallmax),nall(nlocmax),
     &     isplit,npop
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      real q,algama,ss,tt

      
      termfsplit  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(isplit,iloc,iall))) +
     &           float(n(isplit,iloc,iall))*
     &           alog(f(isplit,iloc,iall))
            nn = nn + n(isplit,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfsplit = termfsplit + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttemp(isplit))/drifttemp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(isplit,iloc,iall))) +
     &           float(ntemp(isplit,iloc,iall))*
     &           alog(ftemp(isplit,iloc,iall))
            nn = nn + ntemp(isplit,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfsplit = termfsplit - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttemp(npop+1))/drifttemp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(npop+1,iloc,iall))) +
     &           float(ntemp(npop+1,iloc,iall))*
     &           alog(ftemp(npop+1,iloc,iall))
            nn = nn + ntemp(npop+1,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfsplit = termfsplit - tt



      end function termfsplit



*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5
      real function termfsplitbis(isplit,npop,npopmax,
     &     nlocmax,nall,nallmax,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp)
      implicit none
      integer npopmax,nlocmax,nallmax,
     &     n(npopmax,nlocmax,nallmax),
     &     ntemp(npopmax,nlocmax,nallmax),nall(nlocmax),
     &     isplit,npop
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      real q,algama,ss,tt

      
      termfsplitbis  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(isplit,iloc,iall))) +
     &           float(n(isplit,iloc,iall))*
     &           alog(f(isplit,iloc,iall))
            nn = nn + n(isplit,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttemp(isplit))/drifttemp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(isplit,iloc,iall))) +
     &           float(ntemp(isplit,iloc,iall))*
     &           alog(ftemp(isplit,iloc,iall))
            nn = nn + ntemp(isplit,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttemp(npop+1))/drifttemp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(npop+1,iloc,iall))) +
     &           float(ntemp(npop+1,iloc,iall))*
     &           alog(ftemp(npop+1,iloc,iall))
            nn = nn + ntemp(npop+1,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

      end function termfsplitbis




*
*     terme des freq dans le log ratio pour un split 
      real function termfmerge(ihost,irem,npopmax,
     &     nlocmax,nall,nallmax,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp)
      implicit none
      integer npopmax,nlocmax,nallmax,
     &     n(npopmax,nlocmax,nallmax),
     &     ntemp(npopmax,nlocmax,nallmax),nall(nlocmax),
     &     ihost,irem
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      real q,algama,ss,tt

      
      termfmerge  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(ihost,iloc,iall))) +
     &           float(n(ihost,iloc,iall))*
     &           alog(f(ihost,iloc,iall))
            nn = nn + n(ihost,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(irem,iloc,iall))) +
     &           float(n(irem,iloc,iall))*
     &           alog(f(irem,iloc,iall))
            nn = nn + n(irem,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttemp(ihost))/drifttemp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(ihost,iloc,iall))) +
     &           float(ntemp(ihost,iloc,iall))*
     &           alog(ftemp(ihost,iloc,iall))
            nn = nn + ntemp(ihost,iloc,iall)
         enddo
         tt = algama(q + float(nn)) - algama(q) + ss
      enddo
      termfmerge = termfmerge - tt
      end function termfmerge


*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5
      real function termfmergebis(ihost,irem,npopmax,
     &     nlocmax,nall,nallmax,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp)
      implicit none
      integer npopmax,nlocmax,nallmax,
     &     n(npopmax,nlocmax,nallmax),
     &     ntemp(npopmax,nlocmax,nallmax),nall(nlocmax),
     &     ihost,irem
      real f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),drifttemp(npopmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      real q,algama,ss,tt

      
      termfmergebis  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(ihost,iloc,iall))) +
     &           float(n(ihost,iloc,iall))*
     &           alog(f(ihost,iloc,iall))
            nn = nn + n(ihost,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(n(irem,iloc,iall))) +
     &           float(n(irem,iloc,iall))*
     &           alog(f(irem,iloc,iall))
            nn = nn + n(irem,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttemp(ihost))/drifttemp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + algama(fa(iloc,iall) * q) -
     &           algama(fa(iloc,iall) * q +
     &           float(ntemp(ihost,iloc,iall))) +
     &           float(ntemp(ihost,iloc,iall))*
     &           alog(ftemp(ihost,iloc,iall))
            nn = nn + ntemp(ihost,iloc,iall)
         enddo
         tt = algama(float(nall(iloc)) + float(nn)) - 
     &        algama(float(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis - tt
      end function termfmergebis

*
*     log vraisemblance
*
      real function ll(z,nindivmax,nlocmax,nlocmax2,npopmax,
     &     nallmax,nppmax,c,f,indcell)
      implicit none
      integer nindivmax,nlocmax,nlocmax2,npopmax,
     &     z(nindivmax,nlocmax2),nppmax,c(nppmax),nallmax,
     &     indcell(nindivmax)
      real f(npopmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,ipop

      ll = 0
      do iindiv = 1,nindivmax
         ipop = c(indcell(iindiv))
         do iloc = 1,nlocmax
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
               ll = ll + alog(f(ipop,iloc,iall1))
            endif
            if(iall2 .ne. -999) then 
               ll = ll + alog(f(ipop,iloc,iall2))
               if((iall1 .ne. iall2) .and. (iall1 .ne. -999))  then 
                  ll = ll + alog(2.)
               endif
            endif
         enddo
      enddo
      end function ll

 
*
*     log de la proba a posteriori du vecteur de parametres
*
      real function lpp(lambda,z,npop,npp,drift,f,fa,c,nppmax,
     &     nindivmax,nlocmax2,npopmax,nlocmax,nallmax,indcell,nall,
     &     fmodel)
      implicit none
      integer nindivmax,nlocmax2,npop,npopmax,nlocmax,nallmax,
     &     npp,nppmax,z(nindivmax,nlocmax2),indcell(nindivmax),
     &     c(nppmax),nall(nlocmax),fmodel
      real drift(npopmax),f(npopmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),lambda
      integer ipp,ipop,iloc,iall
      real algama,shape1,shape2,ll

      parameter(shape1=2.,shape2=20.)

      lpp = -lambda + npp*alog(lambda) 
      do ipp=1,npp
         lpp = lpp - alog(float(ipp))
      enddo
      lpp = lpp -npp*alog(float(npop))
      
      if(fmodel .eq. 1) then
         do ipop = 1,npop
            lpp = lpp + 
     &           algama(shape1+shape2) - algama(shape1) - algama(shape2)
     &           + (shape1-1)*alog(drift(ipop)) + 
     &           (shape2-1)*alog(1-drift(ipop))
            do iloc = 1,nlocmax
               lpp = lpp + algama((1-drift(ipop))/drift(ipop))
               do iall=1,nall(iloc)
                  lpp = lpp -algama(fa(iloc,iall)*
     &                 (1-drift(ipop))/drift(ipop)) + 
     &                 (fa(iloc,iall)* 
     &                 (1-drift(ipop))/drift(ipop)-1)*
     &              f(ipop,iloc,iall)
               enddo
            enddo
         enddo
      endif
      
      lpp = lpp + ll(z,nindivmax,nlocmax,nlocmax2,npopmax,
     &     nallmax,nppmax,c,f,indcell)

      end function lpp





***********************************************************************
*     
*     sample freq from full contitionnal pi(f|u,c,z)
*     drifts and ancestral are not changed
      subroutine samplef(npop,npopmax,nloc,nlocmax,
     &     nall,nallmax,ipop1,ipop2,f,ftemp,
     &     fa,drift,a,ptemp,ntemp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nall,
     &     nallmax,ntemp,
     &     ipop1,ipop2
      real f,drift,ftemp,fa,a
      dimension nall(nlocmax),
     &     ntemp(npopmax,nloc,nallmax),
     &     f(npopmax,nlocmax,nallmax),drift(npopmax),
     &     ftemp(npopmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),a(nallmax),ptemp(nallmax)
      integer iloc,iall,ipop
      real ptemp


*     init ftemp 
      do ipop = 1,npop
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(ipop,iloc,iall) = f(ipop,iloc,iall)
            enddo 
         enddo
      enddo


*     new freq
*     for pop ipop1
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drift(ipop1))/
     &           drift(ipop1)+float(ntemp(ipop1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(ipop1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

*     for pop ipop2
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drift(ipop2))/
     &           drift(ipop2)+float(ntemp(ipop2,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(ipop2,iloc,iall)  = ptemp(iall)
         enddo
      enddo
      end subroutine  samplef 



***********************************************************************
*
*     log of ratio of proposals in a joint update of c and f
*     (in subroutine udcf)
      real function lrpf(npopmax,nloc,nall,nallmax,n,
     &     ntemp,f,ftemp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nall,nallmax,n,ntemp,ipop1,ipop2
      real f,ftemp
      dimension nall(nloc),n(npopmax,nloc,nallmax),
     &     ntemp(npopmax,nloc,nallmax),f(npopmax,nloc,nallmax),
     &     ftemp(npopmax,nloc,nallmax)
      integer iloc,iall,n1,n2,ntemp1,ntemp2
      real algama
      lrpf = 0.e0
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntemp1 = 0
         ntemp2 = 0
         do iall=1,nall(iloc)
            lrpf = lrpf + 
     &         n(ipop1,iloc,iall)*alog(f(ipop1,iloc,iall)) 
     &       + n(ipop2,iloc,iall)*alog(f(ipop2,iloc,iall) )
     &       - ntemp(ipop1,iloc,iall)*alog(ftemp(ipop1,iloc,iall)) 
     &       - ntemp(ipop2,iloc,iall)*alog(ftemp(ipop2,iloc,iall)) 
     &       - algama(1+float(n(ipop1,iloc,iall)))
     &       - algama(1+float(n(ipop2,iloc,iall)))
     &       + algama(1+float(ntemp(ipop1,iloc,iall)))
     &       + algama(1+float(ntemp(ipop2,iloc,iall)))
            n1 = n1 + n(ipop1,iloc,iall)
            n2 = n2 + n(ipop2,iloc,iall)
            ntemp1 = ntemp1 + ntemp(ipop1,iloc,iall)
            ntemp2 = ntemp2 + ntemp(ipop2,iloc,iall)
         enddo
         lrpf = lrpf + algama(float(nall(iloc)+n1)) 
     &           + algama(float(nall(iloc)+n2)) 
     &           - algama(float(nall(iloc)+ntemp1)) 
     &           - algama(float(nall(iloc)+ntemp2))
      enddo
      end 




***********************************************************************
*
*     log of ratio of proposals x priors in a joint update of c and f
*     (in subroutine udcf2)
*     
      real function lrppf(npopmax,nloc,nall,nallmax,n,
     &     ntemp,fa,drift,f,ftemp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nall,nallmax,n,ntemp,ipop1,ipop2
      real fa,drift,f,ftemp
      dimension nall(nloc),n(npopmax,nloc,nallmax),
     &     ntemp(npopmax,nloc,nallmax),f(npopmax,nloc,nallmax),
     &     fa(nloc,nallmax),ftemp(npopmax,nloc,nallmax),
     &     drift(npopmax)
      integer iloc,iall,n1,n2,ntemp1,ntemp2
      real algama,q1,q2
      lrppf = 0.e0
      q1 = drift(ipop1)/(1-drift(ipop1))
      q2 = drift(ipop2)/(1-drift(ipop2))
      
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntemp1 = 0
         ntemp2 = 0
         do iall=1,nall(iloc)
            lrppf = lrppf + 
     &         n(ipop1,iloc,iall)*alog(f(ipop1,iloc,iall)) 
     &       + n(ipop2,iloc,iall)*alog(f(ipop2,iloc,iall) )
     &       - ntemp(ipop1,iloc,iall)*alog(ftemp(ipop1,iloc,iall)) 
     &       - ntemp(ipop2,iloc,iall)*alog(ftemp(ipop2,iloc,iall)) 
     &       + algama(fa(iloc,iall)*q1+
     &           float(ntemp(ipop1,iloc,iall)))
     &       + algama(fa(iloc,iall)*q2+
     &           float(ntemp(ipop2,iloc,iall)))
     &       - algama(fa(iloc,iall)*q1+
     &        float(n(ipop1,iloc,iall)))
     &       - algama(fa(iloc,iall)*q2+
     &        float(n(ipop2,iloc,iall)))
            n1 = n1 + n(ipop1,iloc,iall)
            n2 = n2 + n(ipop2,iloc,iall)
            ntemp1 = ntemp1 + ntemp(ipop1,iloc,iall)
            ntemp2 = ntemp2 + ntemp(ipop2,iloc,iall)
         enddo
         lrppf = lrppf + algama(q1+float(n1)) 
     &           + algama(q2+float(n2)) 
     &           - algama(q1+float(ntemp1)) 
     &           - algama(q2+float(ntemp2))
      enddo
      end 



c      include './randlib-1.3.f'
c      include './algama.f'
   
