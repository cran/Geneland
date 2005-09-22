***********************************************************************
*     limites du rectangle contenant les coordonnees
      subroutine limit(nindiv,nindivmax,s,xlim,ylim,dt)
      implicit none
      integer nindiv,nindivmax
      double precision s(2,nindivmax),xlim(2),ylim(2),dt
      integer iindiv
      xlim(1) = 1.d+30
      xlim(2) = -1.d+30
      ylim(1) = 1.d+30
      ylim(2) = -1.d+30
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
***********************************************************************





***********************************************************************
*     points uniformes dans [0,1]x[0,1]
      subroutine rprioru(npp,nppmax,xlim,ylim,u)
      implicit none
      integer npp,nppmax
      double precision u(2,nppmax),ggrunif,xlim(2),ylim(2)
      integer i
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
      end





***********************************************************************
*     affectation dans les classes selon une loi uniforme
      subroutine rpriorc(npp,nppmax,nclass,c)
      implicit none
      integer npp,nppmax,nclass,c(nppmax)
      double precision ggrunif
      integer i
      do i=1,npp
         c(i) = 1+ int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1,nppmax
            c(i) = -999
         enddo
      endif
      end


***********************************************************************
*     init de la dérive 
*     selon un prior uniforme sur [0,1]
      subroutine rpriordrift(nclass,nclassmax,drift,fmodel)
      implicit none
      integer nclass,nclassmax,fmodel
      double precision drift(nclassmax)
      integer iclass
      double precision ggrunif
      if(fmodel .eq. 0) then
         do iclass=1,nclass
            drift(iclass) = 0.5d0
         enddo
      else
         do iclass=1,nclass
            drift(iclass) = ggrunif(0.d0,1.d0)
         enddo
      endif
      if(nclassmax .gt. nclass) then
         do iclass=nclass+1,nclassmax
            drift(iclass) = -999
         enddo
      endif
      end subroutine rpriordrift


***********************************************************************
*     simulation d'une Dirichlet(1,...,1)
*     (p(1),...,p(k)) uniforme dans l'ensemble {p(1)+...+p(k)=1}
      subroutine dirichlet1(nall,nallmax,p)
      implicit none
      integer nall,nallmax
      double precision p(nallmax)
      integer i
      double precision s,ggrexp
      s = 0.d0
      do i=1,nall
         p(i) = ggrexp(1.d0)
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


***********************************************************************
*     simulation d'une Dirichlet(a1,...,an)
      subroutine dirichlet(n,nmax,a,p)
      implicit none
      integer n,nmax
      double precision a(nmax),p(nmax)
      integer i
      double precision s,ggrgam
      s = 0.d0
      do i=1,n
         p(i) = 0.d0
         do while(p(i) .lt. 1.d-38) 
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
      end

 



***********************************************************************     
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





***********************************************************************
*     from numerical recipes p 233
      subroutine indexx(n,nmax,arrin,indx)
      implicit none
      integer n,nmax,indx,indxt
      double precision arrin,arrtmp
      dimension arrin(nmax),indx(nmax),arrtmp(nmax)
      integer i,j,l,ir,q

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



***********************************************************************
*     tirage des frequences dans toutes les classes
*     a tous les locus
      subroutine rpriorf(nclass,nclassmax,nloc,nlocmax,nall,nallmax,f,
     &     ptemp)
      implicit none
      integer nclass,nclassmax,nloc,nlocmax,nall(nlocmax),nallmax
      double precision f(nclassmax,nlocmax,nallmax)
      integer k,l,i
      double precision ptemp(nallmax)
      do  k=1,nclass
         do l=1,nloc
            call dirichlet1(nall(l),nallmax,ptemp)
            do i=1,nallmax
               f(k,l,i)  = ptemp(i)
            enddo
         enddo
      enddo
      if(nclassmax .gt. nclass) then
         do k=nclass+1, nclassmax
            do l=1,nloc
               do i=1,nallmax
                  f(k,l,i)  = -999
               enddo
            enddo
         enddo
      endif
      end subroutine rpriorf


***********************************************************************
*     tirage des frequences dans la pop ancestrale 
*     a tous les locus
      subroutine rpriorfa(nloc,nlocmax,nall,nallmax,fa,fmodel,ptemp)
      implicit none
      integer nloc,nlocmax,nall(nlocmax),nallmax,fmodel
      double precision fa(nlocmax,nallmax)
      integer l,i
      double precision ptemp(nallmax)
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



***********************************************************************
*     Mise a jour gibbsienne de f
*     prior p(f) Dirichlet(1,...,1) 
*     p(f|...) Dirichlet(1+ n1,..., 1+np)
*     ni = nbre d'alleles observes
*     (Cf Green 95, p738) 
*     (et on fait comme si on avait 2*n individus haploides) 
      subroutine rpostf(nclass,nclassmax,nloc,nlocmax,nall,nallmax,f,
     &     nindiv,nindivmax,nlocmax2,z,npp,nppmax,c,indcell,n,a,f11temp)
      implicit none 
      integer nclass,nclassmax,nloc,nlocmax,nall(nlocmax),nallmax,
     &     nindiv,nindivmax,nlocmax2,npp,nppmax,c(nppmax),
     &     indcell(nindivmax)
      double precision f(nclassmax,nlocmax,nallmax)
      integer iclass,iloc,iindiv,iall,n(nclassmax,nlocmax,nallmax),
     &     z(nindivmax,nlocmax2)
      double precision a(nallmax),f11temp(nallmax)

*     comptage des effectifs
      do iclass = 1,nclass
         do iloc = 1,nloc
            do iall =1,nall(iloc)
               n(iclass,iloc,iall)=0
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
      do iclass = 1,nclass
         do iloc = 1,nloc
            do iall =1,nall(iloc)
               a(iall) = 1+dble(n(iclass,iloc,iall))
            enddo
            call dirichlet(nall(iloc),nallmax,a,f11temp)
            do iall =1,nall(iloc)
               f(iclass,iloc,iall) = f11temp(iall)
            enddo
         enddo
      enddo
      end




***********************************************************************
*     Mise a jour gibbsienne de f
*     prior p(f) F-model
*     n = nbre d'alleles observes
      subroutine rpostf2 (nclass,nclassmax,nloc,nlocmax,nall,nallmax,
     &     f,fa,drift,
     &     nindiv,nindivmax,nlocmax2,z,npp,nppmax,c,indcell,n,a,ptemp,
     &     ploidy)
      implicit none 
      integer nclass,nclassmax,nloc,nlocmax,nall(nlocmax),nallmax,
     &     nindiv,nindivmax,nlocmax2,npp,nppmax,c(nppmax),
     &     indcell(nindivmax),ploidy
      double precision f(nclassmax,nlocmax,nallmax),fa(nlocmax,nallmax),
     &     drift(nclassmax)
      integer iclass,iloc,iindiv,iall,n(nclassmax,nlocmax,nallmax),
     &     z(nindivmax,nlocmax2)
      double precision a(nallmax),ptemp(nallmax)

c$$$      write(*,*) ''
c$$$      write(*,*) ''
c$$$      write(*,*) 'dans rpostf2'
c$$$      write(*,*) 'nclass=',nclass
c$$$      write(*,*) 'nclassmax=',nclassmax
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
      do iclass = 1,nclass
         do iloc = 1,nloc
            do iall =1,nall(iloc)
               n(iclass,iloc,iall)=0
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
      do iclass = 1,nclass
c         write(*,*) 'iclass=', iclass
         do iloc = 1,nloc
c            write(*,*) 'iloc=',iloc
            do iall =1,nall(iloc)
c               write(*,*) 'iall=',iall
               if(ploidy .eq. 1) then
                  a(iall) = fa(iloc,iall)*(1/drift(iclass)-1) + 
     &                 dble(n(iclass,iloc,iall))/2
               endif
               if(ploidy .eq. 2) then
                  a(iall) = fa(iloc,iall)*(1/drift(iclass)-1) + 
     &                 dble(n(iclass,iloc,iall))
               endif
c               write(*,*) 'a(',iall,')=',a(iall) 
c               write(*,*) 'fa(',iloc,',',iall,')=', fa(iloc,iall)
c               write(*,*) 'drift(',iclass,')=',drift(iclass)
c               write(*,*) 'n(',iclass,',',iloc,',',iall,')=', 
c     &              n(iclass,iloc,iall)
            enddo
c            write(*,*) 'fa=',fa
c            write(*,*) 'drift=',drift
c            write(*,*) 'a=',a
            call dirichlet(nall(iloc),nallmax,a,ptemp)
            do iall =1,nall(iloc)
               f(iclass,iloc,iall) = ptemp(iall)
            enddo
         enddo
      enddo
      end subroutine rpostf2






***********************************************************************
*
*     Mise à jour M-H des freq allelique de la pop ancestrale 
*     (F-Model de Falush et al.)
*     fa admet un prior Dirichlet(1,...,1)
      subroutine updfa(nclass,nclassmax,nlocmax,nallmax,nall,f,fa,drift)
      implicit none 
      integer nclass,nclassmax,nlocmax,nallmax,nall(nlocmax)
      double precision f(nclassmax,nlocmax,nallmax),fa(nlocmax,nallmax),
     &     drift(nclassmax)
      integer iloc,iall1,iall2,iclass
      double precision delta,ggrunif,ggrnorm,sigdelta,fa1,fa2,ratio,
     &     lratio,gglgammafn,q,u
      parameter(sigdelta = 0.05) 
      
*     boucle sur les loci
      do iloc=1,nlocmax

*     tirage des deux formes alleliques dont les freq seront 
*     mises à jour 
         iall1 = 1+ int(dint(dble(nall(iloc))*ggrunif(0.d0,1.d0)))
         iall2 = iall1 

         do while(iall2 .eq. iall1)
c            write(*,*) 'dans le while'
            iall2 = 1+ int(dint(dble(nall(iloc))*ggrunif(0.d0,1.d0)))
         enddo

*     tirage de l'increment
         delta = ggrnorm(0.d0,1.d0)*sigdelta

*     perturbation des deux freq
         fa1 = fa(iloc,iall1) + delta
         fa2 = fa(iloc,iall2) - delta
         if(((fa1 .gt. 1.d-38) .and. (1-fa1 .gt. 1.d-38)) .and.
     &      ((fa2 .gt. 1.d-38) .and. (1-fa2 .gt. 1.d-38))) then 
*     calcul du log du ratio 
            lratio = 0.d0
            do iclass = 1,nclass
               q = (1.d0-drift(iclass))/drift(iclass)
               lratio = lratio 
     &              + gglgammafn(fa(iloc,iall1)*q)-gglgammafn(fa1*q)
     &              + gglgammafn(fa(iloc,iall2)*q)-gglgammafn(fa2*q)
     &              + delta*q
     &              *log(f(iclass,iloc,iall1)/f(iclass,iloc,iall2))
            enddo
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)

c$$$            write(*,*) 'delta=', delta
c$$$            write(*,*) 'fa1=',fa1
c$$$            write(*,*) 'fa2=',fa2
c$$$            write(*,*) 'q=',q
c$$$            write(*,*) 'fa(iloc,iall2)*q=',fa(iloc,iall2)*q
c$$$            write(*,*) 'gamma(fa(iloc,iall1)*q)',gamma(fa(iloc,iall1)*q)
c$$$            write(*,*) 'ratio=',ratio
c$$$            write(*,*) 'delta*q=',delta*q
c$$$            write(*,*) ''

            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               fa(iloc,iall1) = fa1 
               fa(iloc,iall2) = fa2
            endif
         endif 
      enddo
      end subroutine updfa
      


***********************************************************************
*
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. uniforme sur chaque composante
      subroutine upddrift(nclass,nclassmax,nlocmax,nallmax,nall,
     &     f,fa,drift)
      implicit none 
      integer nclass,nclassmax,nlocmax,nallmax,nall(nlocmax)
      double precision drift(nclassmax),f(nclassmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax)
      integer iclass,iloc,iall
      double precision d,q,qtemp,sigdelta,ratio,lratio,shape1,shape2,
     &     sall,ggrnorm,
     &     gglgammafn,u,ggrunif
      parameter(sigdelta = 0.01,shape1=2.d0,shape2=20.d0) 

*     boucle sur les classes
      do iclass=1,nclass
*     proposition nouvelle valeur
         d = drift(iclass) + ggrnorm(0.d0,1.d0)*sigdelta
         q = (1-drift(iclass))/drift(iclass)
         qtemp = (1-d)/d
         if((d .gt. 1.d-38) .and. (1-d .gt. 1.d-38)) then 

*     calcul du log du ratio
c     prior uniforme
            lratio = 0 
c     prior beta(shape1,shape2)
            lratio = (shape1-1)*dlog(d/drift(iclass)) + 
     &           (shape2-1)*dlog((1-d)/(1-drift(iclass)))
            do iloc=1,nlocmax
               sall = 0.d0
               do iall = 1,nall(iloc)
                  sall = sall + gglgammafn(fa(iloc,iall)*q)-
     &                 gglgammafn(fa(iloc,iall)*qtemp) +
     &                 fa(iloc,iall)*(qtemp-q)*dlog(f(iclass,iloc,iall))
               enddo
               lratio = lratio +sall+ (gglgammafn(qtemp)-gglgammafn(q))
            enddo
 

            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               drift(iclass) = d 
            endif
         endif
      enddo
      end subroutine upddrift



***********************************************************************
*
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. uniforme sur chaque composante
      subroutine upddrift2(nclass,nclassmax,nlocmax,nallmax,nall,
     &     f,fa,drift)  
      implicit none 
      integer nclass,nclassmax,nlocmax,nallmax,nall(nlocmax)
      double precision drift(nclassmax),f(nclassmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax)
      integer iclass,iloc,iall
      double precision d,q,qtemp,sigdelta,ratio,lratio,alpha,sall,
     &     ggrnorm,gglgammafn,u,ggrunif
      parameter(sigdelta = 0.01,alpha=5000) 

*     boucle sur les classes
      do iclass=1,nclass
*     proposition nouvelle valeur
         d = drift(iclass) + ggrnorm(0.d0,1.d0)*sigdelta
         q = (1-drift(iclass))/drift(iclass)
         qtemp = (1-d)/d
         if((d .gt. 1.d-38) .and. (1-d .gt. 1.d-38)) then 

*     calcul du log du ratio
            lratio = 0 
c     decommenter la ligne suivante pour avoir un prior exponentiel tronqué
c     sinon le prior est uniforme
c            lratio = -alpha*(d-drift(iclass))
            do iloc=1,nlocmax
               sall = 0.d0
               do iall = 1,nall(iloc)
                  sall = sall + gglgammafn(fa(iloc,iall)*q)-
     &                 gglgammafn(fa(iloc,iall)*qtemp) +
     &                 fa(iloc,iall)*(qtemp-q)*dlog(f(iclass,iloc,iall))
               enddo
               lratio = lratio +sall+ (gglgammafn(qtemp)-gglgammafn(q))
            enddo


            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               drift(iclass) = d 
            endif
         endif
      enddo
      end subroutine upddrift2





***********************************************************************
*     recherche la cellule de chaque individu
*     stockage des indices dans indcell
*     stockage des carres des distances dans distcell
      subroutine calccell(nindiv,nindivmax,s,npp,nppmax,u,
     &     indcell,distcell)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax)
      double precision distcell(nindivmax)
      double precision s(2,nindivmax),u(2,nppmax)
      integer iindiv,ipp
      double precision d
      do iindiv=1,nindiv
         indcell(iindiv) = -999
         distcell(iindiv) = 1.d+36
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





***********************************************************************
*     mise a jour de indcell et distcell
*     apres le deplacement d'un point de u (celui d'indice j)
      subroutine vormove(nindiv,nindivmax,s,npp,nppmax,u,
     &     indcell,distcell,indcelltemp,distcelltemp,j)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),j
      double precision s(2,nindivmax),u(2,nppmax),distcell(nindivmax),
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




***********************************************************************
*     mise a jour de indcell et distcell
*     apres naissance d'un point de u 
      subroutine voradd(s,u,c,utemp,ctemp,
     &     indcell,distcell,indcelltemp,distcelltemp,
     &     nindiv,nindivmax,npp,nppmax)
      implicit none 
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),iindiv,c(nppmax),ctemp(nppmax)
      double precision s(2,nindivmax),u(2,nppmax),distcell(nindivmax),
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

 


     

***********************************************************************
*     mise a jour de indcell et distcell
*     apres mort d'un point de u 
      subroutine vorrem(s,u,c,utemp,ctemp,ipprem,
     &     indcell,distcell,indcelltemp,distcelltemp,
     &     nindiv,nindivmax,npp,nppmax)
      implicit none
      integer nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     indcelltemp(nindivmax),ipprem,iindiv,c(nppmax),ctemp(nppmax)
      double precision s(2,nindivmax),u(2,nppmax),utemp(2,nppmax),
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





***********************************************************************
*
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
         rpostlamb = ggrgam(dble(m),1.d0)
      enddo
*      write(*,*) 'end rpostlamb'
      end






*****************************************************************
*     
*     Mise a jour de c sans modif de npp
*
      subroutine updc(npp,nppmax,c,ctemp,z,nindiv,nindivmax,nloc,
     &     nlocmax,nlocmax2,nallmax,nclass,nclassmax,f,indcell,ploidy)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,nclass,nallmax,nclassmax,z(nindivmax,nlocmax2),
     &     indcell(nindivmax),ploidy
      double precision f(nclassmax,nlocmax,nallmax)
      integer ipp,ctemp(nppmax)
      double precision ggrunif,r,alpha,ratio,ggrbinom,bern
      integer iclass,iloc,iall

c      write(*,*) 'debut de updc'

c      write(*,*) 'c=',c

      do ipp=1,npp
         ctemp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctemp(ipp) = -999
         enddo
      endif
      do ipp=1,npp

         ctemp(ipp) = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
c         write(*,*) 'ctemp=',ctemp
c         write(*,*) 'avant ratio'
         r = ratio(z,f,c,ctemp,indcell,indcell,
     &     nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'apres ratio'
c         write(*,*) 'r=',r
         alpha = dmin1(1.d0,r)
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
         else 
            ctemp(ipp) = c(ipp)
         endif
      enddo

c      write(*,*) 'fin de updc'
      end subroutine updc





c$$$*****************************************************************
c$$$*     Modification de u
c$$$*     composante par composante 
c$$$*     avec proposal uniforme sur le domaine (independence sampler)
c$$$      subroutine upduis(npp,nppmax,c,u,z,nindiv,nindivmax,nloc,nlocmax,
c$$$     &     nlocmax2,nallmax,nclass,nclassmax,f,indcell,distcell,
c$$$     &     s,xlim,ylim)
c$$$      implicit none 
c$$$      integer npp, nppmax,c(nppmax),nindiv,nindivmax,nloc,nlocmax,
c$$$     &     nlocmax2,nclass,nallmax,nclassmax,z(nindivmax,nlocmax2),
c$$$     &     indcell(nindivmax)
c$$$      double precision u(2,nppmax),f(nclassmax,nlocmax,nallmax),
c$$$     &     distcell(nindivmax),s(2,nindivmax),xlim(2),ylim(2)
c$$$      integer ipp,iclass,iclasstemp,iall1,iall2,iindiv,
c$$$     &     iloc,indcelltemp(nindivmax)
c$$$      double precision utemp(2,nppmax),ggrunif,r,alpha,ggrbinom,bern,
c$$$     &     distcelltemp(nindivmax)
c$$$
c$$$      do ipp=1,npp
c$$$         utemp(1,ipp) = u(1,ipp)
c$$$         utemp(1,ipp) = u(1,ipp)
c$$$      enddo
c$$$      if(nppmax .gt. npp) then
c$$$         do ipp=npp+1,nppmax
c$$$            utemp(1,ipp) = -999.
c$$$            utemp(1,ipp) = -999.
c$$$         enddo
c$$$      endif
c$$$
c$$$      do ipp=1,npp
c$$$         utemp(1,ipp) = xlim(1)+(xlim(2)-xlim(1))*ggrunif(0.d0,1.d0)
c$$$         utemp(2,ipp) = ylim(1)+(ylim(2)-ylim(1))*ggrunif(0.d0,1.d0)
c$$$         call vormove(nindiv,nindivmax,s,npp,nppmax,utemp,
c$$$     &        indcell,distcell,indcelltemp,distcelltemp,ipp)
c$$$         r =1.d0
c$$$         do iindiv=1,nindiv
c$$$            iclass = c(indcell(iindiv))
c$$$            iclasstemp = c(indcelltemp(iindiv))
c$$$            do iloc=1,nloc
c$$$               iall1 = z(iindiv,2*iloc-1)
c$$$               iall2 = z(iindiv,2*iloc)
c$$$               r = r*
c$$$     &              (f(iclasstemp,iloc,iall1)/f(iclass,iloc,iall1))*
c$$$     &              (f(iclasstemp,iloc,iall2)/f(iclass,iloc,iall2))
c$$$            enddo
c$$$         enddo
c$$$         alpha = dmin1(1.d0,r)
c$$$         bern = ggrbinom(1.d0,alpha)
c$$$         if(bern .eq. 1) then
c$$$            u(1,ipp) = utemp(1,ipp)
c$$$            u(2,ipp) = utemp(2,ipp)
c$$$            do iindiv=1,nindiv
c$$$               indcell(iindiv) = indcelltemp(iindiv)
c$$$               distcell(iindiv) = distcelltemp(iindiv)
c$$$            enddo
c$$$         else 
c$$$            utemp(1,ipp) = u(1,ipp)
c$$$            utemp(2,ipp) = u(2,ipp)
c$$$         endif
c$$$      enddo
c$$$      end subroutine upduis






***********************************************************************
*     Modification de u
*     composante par composante 
*     avec proposal uniforme sur un carre de cote du 
*     centre sur le point courant (random walk)
      subroutine updurw(npp,nppmax,c,u,z,nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,nallmax,nclass,nclassmax,f,indcell,distcell,
     &     indcelltemp,distcelltemp,
     &     s,xlim,ylim,du,ploidy)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nindivmax,nloc,nlocmax,
     &     nlocmax2,nclass,nallmax,nclassmax,z(nindivmax,nlocmax2),
     &     indcell(nindivmax),ploidy
      double precision u(2,nppmax),f(nclassmax,nlocmax,nallmax),
     &     distcell(nindivmax),s(2,nindivmax),xlim(2),ylim(2),du
      integer ipp,iindiv,indcelltemp(nindivmax)
      double precision utemp(2,nppmax),ggrunif,r,alpha,ggrbinom,bern,
     &     distcelltemp(nindivmax),surf,surftemp,dx,dy,ratio

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
c         write(*,*) 'u(1,ipp)=', u(1,ipp)
c         write(*,*) 'u(2,ipp)=', u(2,ipp)
*     proposition d un deplacement d un point de u
         utemp(1,ipp) = max(u(1,ipp)-du/2.d0,xlim(1)) + 
     &        ggrunif(0.d0,1.d0)*
     &        (min(u(1,ipp)+du/2.d0,xlim(2))-
     &         max(u(1,ipp)-du/2.d0,xlim(1)))
         utemp(2,ipp) = max(u(2,ipp)-du/2.d0,ylim(1)) + 
     &        ggrunif(0.d0,1.d0)*
     &        (min(u(2,ipp)+du/2.d0,ylim(2))-
     &        max(u(2,ipp)-du/2.d0,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.d0,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.d0,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
c         write(*,*) 'dx=',dx
c         write(*,*) 'dy=',dy
         surf = (dx+du/2.d0)*(dy+du/2.d0)
         dx = min(du/2.d0,utemp(1,ipp)-xlim(1),xlim(2)-utemp(1,ipp))
         dy = min(du/2.d0,utemp(2,ipp)-ylim(1),ylim(2)-utemp(2,ipp))
         surftemp = (dx+du/2.d0)*(dy+du/2.d0)

*     modif de indcell et distcell
         call vormove(nindiv,nindivmax,s,npp,nppmax,utemp,
     &        indcell,distcell,indcelltemp,distcelltemp,ipp)

c         write(*,*) 'apres vormove'



         r = ratio(z,f,c,c,indcell,indcelltemp,
     &     nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'r=',r
         r = r*surf/surftemp
c         write(*,*) 'surf=',surf
c         write(*,*) 'surftemp=',surftemp
         alpha = dmin1(1.d0,r)
c         write(*,*) 'alpha=',alpha
         bern = ggrbinom(1.d0,alpha)
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





***********************************************************************
*
*     mise a jour de t    
* 
      subroutine updt(npp,nppmax,nindiv,
     &     nindivmax,nloc,nlocmax,nlocmax2,nallmax,nclassmax,
     &     t,ttemp,dt,s,c,indcell,distcell,indcelltemp,distcelltemp,
     &     u,z,f,ploidy)
      implicit none 
      integer npp,nppmax,nindiv,nindivmax,nloc,nlocmax,nlocmax2,nallmax,
     &     nclassmax,c(nppmax),indcell(nindivmax),z(nindivmax,nlocmax2),
     &     ploidy
      double precision t(2,nindivmax),s(2,nindivmax),
     &     distcell(nindivmax),u(2,nppmax),f(nclassmax,nlocmax,nallmax),
     &     dt
      integer iindiv,ipp,indcelltemp(nindivmax)
      double precision ggrunif,d,ttemp(2,nindivmax),r,alpha,
     &     distcelltemp(nindivmax),ratio,ggrbinom,accept

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
         ttemp(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         ttemp(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)

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
     &           nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,
     &           nlocmax2,nppmax,ploidy)
         else 
            r = 1.d0
         endif
         alpha = dmin1(1.d0,r)
         accept = ggrbinom(1.d0,alpha)
*     mise a jour en cas d'acceptation
         if(accept .eq. 1) then 
            indcell(iindiv) = indcelltemp(iindiv)
            distcell(iindiv) = distcelltemp(iindiv)
            t(1,iindiv) = ttemp(1,iindiv) 
            t(2,iindiv) = ttemp(2,iindiv)
         endif
      enddo
      end subroutine updt






***********************************************************************
*
*     naissance ou mort d'une cellule
*     avec prior Poisson(lambda) tronquée :   1 < m < nppmax
      subroutine bdpp(nindiv,nindivmax,u,c,utemp,ctemp,nclass,nclassmax,
     &     nloc,nlocmax,nlocmax2,nallmax,npp,nppmax,z,f,s,xlim,ylim,
     &     indcell,distcell,indcelltemp,distcelltemp,lambda,ploidy)
      implicit none 
      integer nindiv,nindivmax,nloc,nlocmax,nlocmax2,
     &     nclass,nclassmax,
     &     nallmax,npp,nppmax,z(nindivmax,nlocmax2),c(nppmax),
     &     indcell(nindivmax),ploidy
      double precision u(2,nindivmax),f(nclassmax,nlocmax,nallmax),
     &     xlim(2),ylim(2),s(2,nindivmax),distcell(nindivmax),lambda

      integer b,ctemp(nppmax),indcelltemp(nindivmax),ipp,npptemp,
     &     iindiv,ipprem
      double precision utemp(2,nppmax),distcelltemp(nindivmax),ggrunif,
     &     ratio,r,accept,alpha,ggrbinom
      
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utemp (1,ipp) = u(1,ipp)
               utemp (2,ipp) = u(2,ipp)
               ctemp(ipp) = c(ipp)
            enddo
            npptemp = npp + 1
            ctemp(npptemp) = 1+ 
     &           int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            utemp(1,npptemp) = xlim(1)+(xlim(2)-xlim(1))*
     &           ggrunif(0.d0,1.d0)
            utemp(2,npptemp) = ylim(1)+(ylim(2)-ylim(1))*
     &           ggrunif(0.d0,1.d0)
            if(nppmax .gt. npptemp) then
               do ipp=npptemp+1,nppmax
                  ctemp(ipp) = -999
                  utemp(1,ipp) = -999.
                  utemp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,u,c,utemp,ctemp,indcell,distcell,indcelltemp,
     &           distcelltemp,nindiv,nindivmax,npp,nppmax)
            r = ratio(z,f,c,ctemp,indcell,indcelltemp,nclassmax,nlocmax,
     &           nallmax,nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
            r = r*lambda/dble(npp+1)
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
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
            ipprem = 1+ dint(dble(npp)*ggrunif(0.d0,1.d0))
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

            call vorrem(s,u,c,utemp,ctemp,ipprem,indcell,distcell,
     &           indcelltemp,distcelltemp,nindiv,nindivmax,npp,nppmax)

            r = ratio(z,f,c,ctemp,indcell,indcelltemp,nclassmax,nlocmax,
     &           nallmax,nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
            r = r*dble(npp)/lambda
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
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







***********************************************************************
*     calcul du ratio p(z|theta*)/p(z|theta)
*     ca ne depend pas de lambda
      double precision function ratio(z,f,c,ctemp,indcell,indcelltemp,
     &     nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
     &     indcell(nindivmax),indcelltemp(nindivmax),ploidy
      double precision f(nclassmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,iclass,iclasstemp,iall

c      write(*,*) 'debut de ratio'
c      write(*,*) 'indcell=',indcell
c      write(*,*) 'indcelltemp=',indcelltemp
c      write(*,*) 'c=',c
c      write(*,*) 'ctemp=',ctemp


      ratio = 1.d0
      do iindiv=1,nindiv
c         write(*,*) 'iindiv=', iindiv
         iclass = c(indcell(iindiv))
         iclasstemp = ctemp(indcelltemp(iindiv))
C         write(*,*) 'indcell=',indcell
C          write(*,*) 'indcelltemp=',indcelltemp
C          write(*,*) 'c=',c
C          write(*,*) 'ctemp=',ctemp
C          write(*,*) 'iclass=',iclass
C          write(*,*) 'iclasstemp=',iclasstemp

         do iloc=1,nloc
c             write(*,*) 'iloc=',iloc
c            write(6,*) 'z=',z(iindiv,2*iloc-1)
c            write(6,*) 'z=',z(iindiv,2*iloc)
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
c               write(*,*) f(iclasstemp,iloc,iall1)
c               write(*,*) f(iclass,iloc,iall1)
               ratio = ratio*
     &              (f(iclasstemp,iloc,iall1)/f(iclass,iloc,iall1))
               
            endif
            if(iall2 .ne. -999) then 
c               write(*,*) f(iclasstemp,iloc,iall2)
c               write(*,*) f(iclass,iloc,iall2)
               ratio = ratio*
     &              (f(iclasstemp,iloc,iall2)/f(iclass,iloc,iall2))
            endif
            if(ratio .lt. 0) then
               write(*,*) 'iloc=',iloc
               write(*,*) 'iindiv=',iindiv
               write(*,*) 'iall1=',iall1
               write(*,*) 'iall2=',iall2
               write(*,*) 
     &              'f(',iclass,',',iloc,',',iall1,')=',
     &              f(iclass,iloc,iall1)
               write(*,*) 
     &              'f(',iclass,',',iloc,',',iall2,')=',
     &              f(iclass,iloc,iall2)
               write(*,*) 
     &              'f(',iclasstemp,',',iloc,',',iall1,')=',
     &              f(iclasstemp,iloc,iall1)
               write(*,*) 
     &              'f(',iclasstemp,',',iloc,',',iall2,')=',
     &              f(iclasstemp,iloc,iall2)
               stop
            endif
         enddo
      enddo
c      write(*,*) 'fin de ratio'
c$$$      if(ratio .lt. 0) then 
c$$$         write(*,*) 'indcell=',indcell
c$$$         write(*,*) 'indcelltemp=',indcelltemp
c$$$         write(*,*) 'c=',c
c$$$         write(*,*) 'ctemp=',ctemp
c$$$         do iclass=1,nclassmax
c$$$            do iloc=1,nlocmax
c$$$               do iall=1,nallmax
c$$$                  write(*,*) 
c$$$     &                 'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                 f(iclass,iloc,iall)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif
      if(ploidy .eq. 1) then 
         ratio = dsqrt(ratio)
      endif
      end function ratio








********************************************************************
*     calcul du ratio p(z|theta*)/p(z|theta)
*     quand f est  modifié
      double precision function ratiobd(z,f,ftemp,c,ctemp,indcell,
     &     indcelltemp,nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,
     &     nlocmax2,nppmax)
      implicit none
      integer nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
     &     indcell(nindivmax),indcelltemp(nindivmax)
      double precision f(nclassmax,nlocmax,nallmax),
     &     ftemp(nclassmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,iclass,iclasstemp

      ratiobd = 1.d0
      do iindiv=1,nindiv
         iclass = c(indcell(iindiv))
         iclasstemp = ctemp(indcelltemp(iindiv))
         do iloc=1,nloc
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
               ratiobd = ratiobd*
     &              (ftemp(iclasstemp,iloc,iall1)/f(iclass,iloc,iall1))
            endif
            if(iall2 .ne. -999) then 
               ratiobd = ratiobd*
     &              (ftemp(iclasstemp,iloc,iall2)/f(iclass,iloc,iall2))
            endif
         enddo
      enddo
c      write(*,*) 'f=',f
c      write(*,*) 'ftemp=',ftemp
c      write(*,*) 'c=',c
c      write(*,*) 'ctemp=',ctemp
c      write(*,*) 'ratiobd=',ratiobd
      end function ratiobd






***********************************************************************
*
*     Indice des cellules dans une classe
*
      subroutine who(c,iclass,npp,nppmax,cellclass,
     &     ncellclass)
      implicit none
      integer npp,nppmax,c(nppmax),iclass,cellclass(nppmax),
     &     ncellclass
      integer ipp,ii
c      write(*,*) 'who'
      ii = 1
      ncellclass = 0
      do ipp=1,npp
         if(c(ipp) .eq. iclass) then
           cellclass(ii) = ipp
           ncellclass = ncellclass + 1
           ii = ii + 1
        endif
      enddo
      if(nppmax .gt. ncellclass) then
         do ipp= ncellclass+1, nppmax
            cellclass(ipp) = -999
         enddo
      endif
      end subroutine who



***********************************************************************
*
*     Tirage de nu cellules parmi ncellclass cellules
*      
      subroutine sample(cellclass,nppmax,nu,ncellclass,listcell)
      implicit none
      integer nppmax,cellclass(nppmax),nu,ncellclass,listcell(nppmax)
      integer isamp,ii,jj
      double precision ggrunif
c      write(*,*) 'sample'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellclass=',ncellclass
c      write(*,*) 'cellclass=',cellclass
*     init
      ii = 1 + int(dint(dble(ncellclass)*ggrunif(0.d0,1.d0)))
      listcell(1) = cellclass(ii)
c      write(*,*) 'listcell(1)=',listcell(1)
      if(nu .gt. 1) then
         do isamp = 2,nu
c             write(*,*) 'cellclass=',cellclass
*     translation
            if(ii .eq. 1) then 
               do jj=2,ncellclass
                  cellclass(jj-1) = cellclass(jj)
               enddo
            else
               if(ii .ne. ncellclass) then
                  do jj=ii+1,ncellclass
                     cellclass(jj-1) = cellclass(jj)
                  enddo
               endif
            endif
            cellclass(ncellclass-(isamp-1)+1) = -999
c             write(*,*) 'cellclass=',cellclass
*     tirage parmi les ncellclass-isamp cellules restantes
            ii = 1 + 
     &           int(dint(dble(ncellclass-isamp)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'ii=',ii
            listcell(isamp) = cellclass(ii)
c            write(*,*) 'listcell(isamp)=',listcell(isamp)
         enddo
      endif
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellclass=',ncellclass
c      write(*,*) 'cellclass=',cellclass
c      write(*,*) 'listcell=',listcell
      end subroutine sample



***********************************************************************
*     split d'une classe en deux
*     double precisionlocation de nu cellules dont les indices
*     sont dans listcell
*     dans la classe iclass
*
      subroutine split(iclass,nclass,c,ctemp,nppmax,nu,listcell)
      implicit none
      integer iclass,nclass,nppmax,c(nppmax),ctemp(nppmax),nu,
     &     listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de split'
c      write(*,*) 'nu=',nu
c      write(*,*) 'iclass=',iclass
c      write(*,*) 'listcell=',listcell
      do ipp=1,nppmax
         ctemp(ipp) = c(ipp)
      enddo
      if(nu .gt. 0) then
         do ii=1,nu
            ctemp(listcell(ii)) = iclass
         enddo
      endif
c      write(*,*)'c=',c
c      write(*,*)'ctemp=',ctemp
c      write(*,*) 'fin de split'
      end subroutine split



***********************************************************************
*     merge de deux  classes en une : 
*     double precisionlocation des nu cellules de la classe iclassrem 
*     dont les indices sont dans listcell
*     dans la classe iclasshost
*
      subroutine merging(iclassrem,iclasshost,
     &     nclass,c,ctemp,nppmax,nu,listcell)
      implicit none
      integer iclassrem,iclasshost,nclass,nppmax,
     &     c(nppmax),ctemp(nppmax),nu,listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de merge'
c      write(*,*) 'nu=',nu
c      write(*,*) 'iclassrem=',iclassrem
      do ipp=1,nppmax
         ctemp(ipp) = c(ipp)
      enddo
      if(iclassrem .gt. iclasshost) then 
         if(nu .gt. 0) then
            do ii=1,nu
               ctemp(listcell(ii)) = iclasshost
            enddo
         endif
      else
         if(nu .gt. 0) then
            do ii=1,nu
               ctemp(listcell(ii)) = iclasshost - 1
            enddo
         endif
      endif
      do ipp=1,nppmax
         if(c(ipp) .gt. iclassrem) ctemp(ipp) = c(ipp)-1
      enddo
c      write(*,*)'c=',c
c      write(*,*)'ctemp=',ctemp
c      write(*,*) 'fin de merge'
      end subroutine merging






***********************************************************************
*
*     Mise a jour de c et f en cas d acceptation d'un split/merge
*
      subroutine accept5(nppmax,nclassmax,nlocmax,nallmax,
     &     nall,c,ctemp,f,ftemp,drift,drifttemp)
      implicit none
      integer nppmax,nclassmax,nlocmax,nallmax,
     &     nall(nlocmax),c(nppmax),ctemp(nppmax)
      double precision f(nclassmax,nlocmax,nallmax),
     &     ftemp(nclassmax,nlocmax,nallmax),
     &     drift(nclassmax),drifttemp(nclassmax)
      integer iclass,iloc,iall,ipp
c      write(*,*) 'debut de accept5'
c      write(*,*) 'f=',f
c      write(*,*) 'ftemp=',ftemp
      do ipp=1,nppmax
         c(ipp) = ctemp(ipp)
      enddo
      do iclass = 1,nclassmax
         do iloc= 1,nlocmax
            do iall=1,nall(iloc)
               f(iclass,iloc,iall) = ftemp(iclass,iloc,iall)
            enddo
         enddo
         drift(iclass) = drifttemp(iclass)
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
      double precision gglgammafn
      bico = dexp(gglgammafn(dble(n+1))-gglgammafn(dble(p+1))-
     &     gglgammafn(dble(n-p+1)))
c      write(*,*) 'in bico '
c$$$      write(*,*) 'n=', n
c$$$      write(*,*) 'p=', p
c$$$      write(*,*) 'gglgammafn(dble(n+1))=',gglgammafn(dble(n+1))
c$$$      write(*,*) 'gglgammafn(dble(p+1))=',gglgammafn(dble(p+1))
c$$$      write(*,*) 'gglgammafn(dble(n-p+1)))=',gglgammafn(dble(n-p+1))
c$$$      write(*,*) 'dexp()=',dexp(gglgammafn(dble(n+1))-gglgammafn(dble(p+1))-
c$$$     &     gglgammafn(dble(n-p+1)))
c$$$      write(*,*) 'bico =', nint(dexp(gglgammafn(dble(n+1))-
c$$$     &     gglgammafn(dble(p+1))-
c$$$     &     gglgammafn(dble(n-p+1))))
c      write(*,*) 'bico =',bico
      
      end function bico


***********************************************************************
*
*     ln du coefficient du binome C_n^p
*
      double precision function lbico(n,p)
      implicit none
      integer n,p
      double precision gglgammafn
      lbico = gglgammafn(dble(n+1))-gglgammafn(dble(p+1))-
     &     gglgammafn(dble(n-p+1))
      end function lbico





***********************************************************************
*
*     log du ratio des vraisemblances dans bdclass6
*
      double precision function llr6(z,f,ftemp,c,ctemp,indcell,
     &     indcelltemp,
     &     nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer nclassmax,nlocmax,nallmax,nindiv,nindivmax,nloc,nlocmax2,
     &     nppmax,z(nindivmax,nlocmax2),c(nppmax),ctemp(nppmax),
     &     indcell(nindivmax),indcelltemp(nindivmax),ploidy
      double precision f(nclassmax,nlocmax,nallmax),
     &     ftemp(nclassmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,iclass,iclasstemp

      llr6 = 0

*     log du rapport des vraisemblances
      do iindiv=1,nindiv
         iclass = c(indcell(iindiv))
         iclasstemp = ctemp(indcelltemp(iindiv))
         do iloc=1,nloc
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftemp(iclasstemp,iloc,iall1)) - 
     &              dlog(f(iclass,iloc,iall1))
            endif
            if(iall2 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftemp(iclasstemp,iloc,iall2)) - 
     &              dlog(f(iclass,iloc,iall2))
            endif
            
c$$$            if(iclasstemp .eq. 3) then 
c$$$               write(*,*) 'iclasstemp=',iclasstemp
c$$$               write(*,*) 'iinidiv=',iindiv
c$$$               write(*,*) 'c=', c(indcell(iindiv))
c$$$               write(*,*) 'ctemp=', ctemp(indcelltemp(iindiv))
c$$$               write(*,*) 'z=',z(iindiv,1),z(iindiv,2)
c$$$               write(*,*) 'llr6 =',llr6
c$$$               write(*,*) 'ftemp(iclasstemp,iloc,iall1)=',
c$$$     &              ftemp(iclasstemp,iloc,iall1) 
c$$$               write(*,*) 'f(iclass,iloc,iall1)=',
c$$$     &              f(iclass,iloc,iall1)
c$$$               write(*,*) 'ftemp(iclasstemp,iloc,iall2)=',
c$$$     &              ftemp(iclasstemp,iloc,iall2)
c$$$               write(*,*) 'f(iclass,iloc,iall2)=',
c$$$     &              f(iclass,iloc,iall2)
c$$$            endif

         enddo
      enddo
      if(ploidy .eq. 1) llr6 = llr6 / 2.d0
      end function llr6



***********************************************************************
*
*     comptage des alleles dans chaque classe pour c
*
      subroutine countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &     nppmax,nall,nallmax,z,n,indcell,c)
      implicit none
      integer  nindiv,nlocmax,nlocmax2,nclassmax,nppmax,nallmax,
     &     nclass,z(nindiv,nlocmax2),nall(nlocmax),
     &     n(nclassmax,nlocmax,nallmax),c(nppmax),indcell(nindiv)
      integer iclass,iloc,iall,iindiv
*     init du tableau
      do iclass = 1,nclassmax
         do iloc = 1,nlocmax
            do iall =1,nall(iloc)
               n(iclass,iloc,iall)=0
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
            if(z(iindiv,2*iloc) .ne. -999) then 
               n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) = 
     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc)) + 1 
c                write(*,*) 'n=',
c     &              n(c(indcell(iindiv)),iloc,z(iindiv,2*iloc-1)) 
            endif
         enddo 
      enddo
c      write(*,*) 'n=',n
c      do iclass=1,nclass
c         do iloc=1,nlocmax
c            write(*,*) 'n(',iclass,iloc,')=',
c     &           (n(iclass,iloc,iall),iall=1,nall(iloc))
c         enddo
c      enddo
      
      end subroutine countn




***********************************************************************
*
*     log du ratio (prob cond. complete)/prior
*     pour les frequences
*     dans un split de la classe iclass
      double precision function lrf(iclass,nclassmax,nlocmax,nall,
     &     nallmax,f,fa,drift,n)
      implicit none
      integer iclass,nclassmax,nlocmax,nall(nlocmax),nallmax,
     &     n(nclassmax,nlocmax,nallmax)
      double precision f(nclassmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),
     &     drift(nclassmax)
      integer iloc,iall,nn
      double precision ss,gglgammafn,q

      lrf = 0.d0
      q = (1-drift(iclass))/drift(iclass)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q +
     &           dble(n(iclass,iloc,iall))) +
     &           (1 - fa(iloc,iall) * q - dble(n(iclass,iloc,iall)))*
     &           dlog(f(iclass,iloc,iall))
            nn = nn + n(iclass,iloc,iall)
         enddo
c         write(*,*) 'nn=',nn
         lrf = lrf + gglgammafn(dble(nall(iloc))) -
     &        gglgammafn(q + nn) + ss

      enddo
      end function lrf



      

***********************************************************************
*
*     Naissance et mort de classes avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdclass7(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &     cellclasshost,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclassrem,ipp,isplit,
     &     cellclass(nppmax),ncellclass,nu,listcell(nppmax),
     &     iclasshost,ncellclasshost,cellclasshost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termfsplit,
     &     termfmerge,b,ggrbinom,bern
      
       do ipp=1,nppmax
          cellclass(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ggrbinom(1.d0,0.5d0)
       
       if(b .eq. 1) then
          if(nclass .lt. nclassmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la classe qui split
             isplit = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))

*     recherche des cellules affectees a cette classe
             call who(c,isplit,npp,nppmax,cellclass,ncellclass)
             if(ncellclass .gt. 0) then
*     tirage du nombre de cellules double precisionlouees
                nu = int(dint(dble(ncellclass)*ggrunif(0.d0,1.d0)))
                if(nu .gt. 0) then

*     tirage des cellules double precisionlouees
                   call sample(cellclass,nppmax,nu,ncellclass,listcell)

*     proposition de double precisionlocation dans la classe nclass+1
                   call split(nclass+1,nclass,c,ctemp,nppmax,nu,
     &                  listcell)

*     comptage des alleles sur chaque locus pour c puis ctemp
                   call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &                  nppmax,nall,nallmax,z,n,indcell,c)
                   call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &                  nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     proposition nouvelle freq et derive 
c                    call addfreq5(isplit,nclass,nclassmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
                   call addfreq7(nclass,nclassmax,nindiv,nloc,nlocmax,
     &     nlocmax2,nall,nallmax,nppmax,indcell,z,isplit,ctemp,f,ftemp,
     &     fa,drift,drifttemp,a,ptemp,ntemp)

*     calcul du log du ratio
*     terme des vraisemblances
                   lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,nclassmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
                   lratio = lratio + termfsplit(isplit,nclass,nclassmax,
     &                  nlocmax,nlocmax2,nall,nallmax,nindiv,
     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)

*     terme des proposal sur c
                   lratio = lratio + dlog(2*dble(ncellclass+1)) + 
     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 

*     terme des priors sur c
                   lratio = lratio + 
     &                  dble(npp)*(dlog(dble(nclass)) - 
     &                  dlog(dble(nclass+1)))

                   lratio = dmin1(0.d0,lratio)
                   alpha = dexp(lratio)
                   bern = ggrbinom(1.d0,alpha)

c$$$                   write(*,*) 'nclass=',nclass
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellclass=',ncellclass  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                           n(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,nclassmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,nclass,nclassmax,
c$$$     &                  nlocmax,nlocmax2,nall,nallmax,nindiv,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellclass+1)) + 
c$$$     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(nclass)) - 
c$$$     &                  dlog(dble(nclass+1)))
c$$$                   write(*,*) 'alpha=',alpha 


                   if(bern .eq. 1) then
                   call accept5(nppmax,nclassmax,nlocmax,
     &                  nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                   nclass = nclass + 1
                endif
               endif
            endif 
         endif

*     merge
      else
         if(nclass .gt. nclassmin) then 
c             write(*,*) 'mort'
*     tirage de la classe qui meurt
            iclassrem = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            
*     tirage de la classe hote
            iclasshost = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            do while(iclasshost .eq. iclassrem)
               iclasshost = 1 + 
     &              int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            enddo

*     on range dans la classe d'indice le plus petit
            if(iclasshost .gt. iclassrem) then
               ii = iclasshost
               iclasshost = iclassrem
               iclassrem = ii
            endif
            
*     recherche des cellules qui vont etre double precisionlouees
            call who(c,iclassrem,npp,nppmax,cellclass,ncellclass)
            
*     recherche des cellules de la classe hote
            call who(c,iclasshost,npp,nppmax,cellclasshost,
     &           ncellclasshost)
               
            if(ncellclass .gt. 0) then
*     proposition de double precisionlocation dans la classe iclasshost
               call merging(iclassrem,iclasshost,nclass,c,ctemp,nppmax,
     &              ncellclass,cellclass)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,n,indcell,c)
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(iclassrem,iclasshost,nclass,nclassmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7(iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &     nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,nclassmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmerge(iclasshost,iclassrem,
     &              nclass,nclassmax,nlocmax,nlocmax2,
     &              nall,nallmax,nindiv,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)

*     terme des proposal sur c
               lratio = lratio + dlog(dble(nclass)) - 
     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
     &              lbico(ncellclass+ncellclasshost,ncellclass) 

*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(nclass)) - 
     &              dlog(dble(nclass-1)))
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
      
         
c$$$               write(*,*) 'nclass=',nclass
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'iclassrem=',iclassrem
c$$$               write(*,*) 'cellclass=',cellclass
c$$$               write(*,*) 'iclasshost=',iclasshost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                       n(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,nclassmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(iclasshost,iclassrem,
c$$$     &              nclass,nclassmax,nlocmax,nlocmax2,
c$$$     &              nall,nallmax,nindiv,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(nclass)) - 
c$$$     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
c$$$     &              lbico(ncellclass+ncellclasshost,ncellclass) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(nclass)) - 
c$$$     &              dlog(dble(nclass-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,nclassmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  nclass = nclass - 1
               endif
            endif
         endif
      endif
      end subroutine bdclass7



***********************************************************************
*
*     Naissance et mort de classes avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdclass8(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &     cellclasshost,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclassrem,ipp,isplit,
     &     cellclass(nppmax),ncellclass,nu,listcell(nppmax),
     &     iclasshost,ncellclasshost,cellclasshost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termfsplit,
     &     termfmerge,b,ggrbinom,bern
      
       do ipp=1,nppmax
          cellclass(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ggrbinom(1.d0,0.5d0)
       
       if(b .eq. 1) then
          if(nclass .lt. nclassmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la classe qui split
             isplit = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
             
*     recherche des cellules affectees a cette classe
             call who(c,isplit,npp,nppmax,cellclass,ncellclass)
             
             if(ncellclass .gt. 0) then
*     tirage du nombre de cellules double precisionlouees
                nu = int(dint(dble(ncellclass)*ggrunif(0.d0,1.d0)))
                if(nu .gt. 0) then
                   
*     tirage des cellules double precisionlouees
                   call sample(cellclass,nppmax,nu,ncellclass,listcell)
                   
*     proposition de double precisionlocation dans la classe nclass+1
                   call split(nclass+1,nclass,c,ctemp,nppmax,nu,
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
             call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &            nppmax,nall,nallmax,z,n,indcell,c)
             call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &            nppmax,nall,nallmax,z,ntemp,indcell,ctemp)
             
*     proposition nouvelle freq et derive 
c     call addfreq5(isplit,nclass,nclassmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
             call addfreq7(nclass,nclassmax,nindiv,nloc,nlocmax,
     &            nlocmax2,nall,nallmax,nppmax,indcell,z,isplit,ctemp,
     &            f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)
             
*     calcul du log du ratio
*     terme des vraisemblances
             lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &            indcell,nclassmax,nlocmax,nallmax,
     &            nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
             lratio = lratio + termfsplit(isplit,nclass,nclassmax,
     &            nlocmax,nlocmax2,nall,nallmax,nindiv,
     &            f,ftemp,n,ntemp,fa,drift,drifttemp,z)
             
*     terme des proposal sur c
             lratio = lratio + dlog(2*dble(ncellclass+1)) + 
     &            lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
             
*     terme des priors sur c
             lratio = lratio + 
     &            dble(npp)*(dlog(dble(nclass)) - 
     &            dlog(dble(nclass+1)))
             
             lratio = dmin1(0.d0,lratio)
             alpha = dexp(lratio)
             bern = ggrbinom(1.d0,alpha)

c$$$                   write(*,*) 'nclass=',nclass
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellclass=',ncellclass  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                           n(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,nclassmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,nclass,nclassmax,
c$$$     &                  nlocmax,nlocmax2,nall,nallmax,nindiv,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellclass+1)) + 
c$$$     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(nclass)) - 
c$$$     &                  dlog(dble(nclass+1)))
c$$$                   write(*,*) 'alpha=',alpha 

             
             if(bern .eq. 1) then
                call accept5(nppmax,nclassmax,nlocmax,
     &               nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                nclass = nclass + 1
             endif
          endif

*     merge
      else
         if(nclass .gt. nclassmin) then 
c             write(*,*) 'mort'
*     tirage de la classe qui meurt
            iclassrem = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            
*     tirage de la classe hote
            iclasshost = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            do while(iclasshost .eq. iclassrem)
               iclasshost = 1 + 
     &              int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            enddo

*     on range dans la classe d'indice le plus petit
            if(iclasshost .gt. iclassrem) then
               ii = iclasshost
               iclasshost = iclassrem
               iclassrem = ii
            endif
            
*     recherche des cellules qui vont etre double precisionlouees
            call who(c,iclassrem,npp,nppmax,cellclass,ncellclass)
            
*     recherche des cellules de la classe hote
            call who(c,iclasshost,npp,nppmax,cellclasshost,
     &           ncellclasshost)
               
            if(ncellclass .gt. 0) then
*     proposition de double precisionlocation dans la classe iclasshost
               call merging(iclassrem,iclasshost,nclass,c,ctemp,nppmax,
     &              ncellclass,cellclass)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,n,indcell,c)
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(iclassrem,iclasshost,nclass,nclassmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7(iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &     nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,nclassmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmerge(iclasshost,iclassrem,
     &              nclass,nclassmax,nlocmax,nlocmax2,
     &              nall,nallmax,nindiv,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)

*     terme des proposal sur c
               lratio = lratio + dlog(dble(nclass)) - 
     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
     &              lbico(ncellclass+ncellclasshost,ncellclass) 

*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(nclass)) - 
     &              dlog(dble(nclass-1)))
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
      
         
c$$$               write(*,*) 'nclass=',nclass
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'iclassrem=',iclassrem
c$$$               write(*,*) 'cellclass=',cellclass
c$$$               write(*,*) 'iclasshost=',iclasshost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                       n(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,nclassmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(iclasshost,iclassrem,
c$$$     &              nclass,nclassmax,nlocmax,nlocmax2,
c$$$     &              nall,nallmax,nindiv,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(nclass)) - 
c$$$     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
c$$$     &              lbico(ncellclass+ncellclasshost,ncellclass) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(nclass)) - 
c$$$     &              dlog(dble(nclass-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,nclassmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  nclass = nclass - 1
               endif
            endif
         endif
      endif
      end subroutine bdclass8







******************************************************************
*     Naissance et mort de classes avec réallocations 
*     (split/merge)
*     drift* reste à 0.5d0
*     pour court-circuiter le F-model
*     proposition de f* selon conditionnelle complète 
*     dans les deux sens
*
      subroutine bdclass7bis(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &     cellclasshost,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclassrem,ipp,isplit,
     &     cellclass(nppmax),ncellclass,nu,listcell(nppmax),
     &     iclasshost,ncellclasshost,cellclasshost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,
     &     termfsplitbis,termfmergebis,b,ggrbinom,bern
      
       do ipp=1,nppmax
          cellclass(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ggrbinom(1.d0,0.5d0)
       
       if(b .eq. 1) then
          if(nclass .lt. nclassmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la classe qui split
             isplit = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))

*     recherche des cellules affectees a cette classe
             call who(c,isplit,npp,nppmax,cellclass,ncellclass)
             if(ncellclass .gt. 0) then
*     tirage du nombre de cellules double precisionlouees
                nu = int(dint(dble(ncellclass)*ggrunif(0.d0,1.d0)))
                if(nu .gt. 0) then

*     tirage des cellules double precisionlouees
                   call sample(cellclass,nppmax,nu,ncellclass,listcell)

*     proposition de double precisionlocation dans la classe nclass+1
                   call split(nclass+1,nclass,c,ctemp,nppmax,nu,
     &                  listcell)

*     comptage des alleles sur chaque locus pour c puis ctemp
                   call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &                  nppmax,nall,nallmax,z,n,indcell,c)
                   call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &                  nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     proposition nouvelle freq et derive 
c                    call addfreq5(isplit,nclass,nclassmax,nloc,nlocmax,
c     &     nall,nallmax,f,ftemp,fa,drift,drifttemp,a,ptemp)
                   call addfreq7bis(nclass,nclassmax,nindiv,nloc,
     &     nlocmax,nlocmax2,nall,nallmax,nppmax,indcell,z,isplit,ctemp,
     &     f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)

*     calcul du log du ratio
*     terme des vraisemblances
                   lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,nclassmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
                   lratio = lratio + termfsplitbis(isplit,nclass,
     &                  nclassmax,nlocmax,nlocmax2,nall,nallmax,nindiv,
     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)

*     terme des proposal sur c
                   lratio = lratio + dlog(2*dble(ncellclass+1)) + 
     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 

*     terme des priors sur c
                   lratio = lratio + 
     &                  dble(npp)*(dlog(dble(nclass)) - 
     &                  dlog(dble(nclass+1)))

                   lratio = dmin1(0.d0,lratio)
                   alpha = dexp(lratio)
                   bern = ggrbinom(1.d0,alpha)

c$$$                   write(*,*) 'nclass=',nclass
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellclass=',ncellclass  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                           n(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,nclassmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,nclass,nclassmax,
c$$$     &                  nlocmax,nlocmax2,nall,nallmax,nindiv,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellclass+1)) + 
c$$$     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(nclass)) - 
c$$$     &                  dlog(dble(nclass+1)))
c$$$                   write(*,*) 'alpha=',alpha 


                   if(bern .eq. 1) then
                   call accept5(nppmax,nclassmax,nlocmax,
     &                  nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                   nclass = nclass + 1
                endif
               endif
            endif 
         endif

*     merge
      else
         if(nclass .gt. nclassmin) then 
c             write(*,*) 'mort'
*     tirage de la classe qui meurt
            iclassrem = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            
*     tirage de la classe hote
            iclasshost = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            do while(iclasshost .eq. iclassrem)
               iclasshost = 1 + 
     &              int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            enddo

*     on range dans la classe d'indice le plus petit
            if(iclasshost .gt. iclassrem) then
               ii = iclasshost
               iclasshost = iclassrem
               iclassrem = ii
            endif
            
*     recherche des cellules qui vont etre double precisionlouees
            call who(c,iclassrem,npp,nppmax,cellclass,ncellclass)
            
*     recherche des cellules de la classe hote
            call who(c,iclasshost,npp,nppmax,cellclasshost,
     &           ncellclasshost)
               
            if(ncellclass .gt. 0) then
*     proposition de double precisionlocation dans la classe iclasshost
               call merging(iclassrem,iclasshost,nclass,c,ctemp,nppmax,
     &              ncellclass,cellclass)

*     comptage des alleles sur chaque locus pour c puis ctemp
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,n,indcell,c)
               call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &              nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(iclassrem,iclasshost,nclass,nclassmax,
c     &              nloc,nlocmax,nall,nallmax,f,ftemp,drift,drifttemp,
c     &              a,fa)
               call remfreq7bis(iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &     nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,ntemp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &              indcell,nclassmax,nlocmax,nallmax,
     &              nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
               lratio = lratio + 
     &              termfmergebis(iclasshost,iclassrem,
     &              nclass,nclassmax,nlocmax,nlocmax2,
     &              nall,nallmax,nindiv,
     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)

*     terme des proposal sur c
               lratio = lratio + dlog(dble(nclass)) - 
     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
     &              lbico(ncellclass+ncellclasshost,ncellclass) 

*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(nclass)) - 
     &              dlog(dble(nclass-1)))
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
      
         
c$$$               write(*,*) 'nclass=',nclass
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'iclassrem=',iclassrem
c$$$               write(*,*) 'cellclass=',cellclass
c$$$               write(*,*) 'iclasshost=',iclasshost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                       n(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,nclassmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(iclasshost,iclassrem,
c$$$     &              nclass,nclassmax,nlocmax,nlocmax2,
c$$$     &              nall,nallmax,nindiv,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(nclass)) - 
c$$$     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
c$$$     &              lbico(ncellclass+ncellclasshost,ncellclass) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(nclass)) - 
c$$$     &              dlog(dble(nclass-1)))
c$$$               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,nclassmax,nlocmax,
     &                 nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
                  nclass = nclass - 1
               endif
            endif
         endif
      endif
      end subroutine bdclass7bis






***********************************************************************
*     split/merge populations in the spatial D-model
*     changes from bdclass7bis:
*     - process populations whatever the number of tiles or individuals
*       they have
      subroutine bdclass8bis(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,cellclass,listcell,
     &     cellclasshost,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclassrem,ipp,isplit,
     &     cellclass(nppmax),ncellclass,nu,listcell(nppmax),
     &     iclasshost,ncellclasshost,cellclasshost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,
     &     termfsplitbis,termfmergebis,b,ggrbinom,bern
      integer iclass,iloc,iall

      do ipp=1,nppmax
         cellclass(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(nclass .lt. nclassmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la classe qui split
            isplit = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'isplit=',isplit

*     recherche des cellules affectees a cette classe
            call who(c,isplit,npp,nppmax,cellclass,ncellclass)

            if(ncellclass .gt. 0) then
*     tirage du nombre de cellules double precisionlouees 
               nu = int(dint(dble(ncellclass+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then                 

*     tirage des cellules double precisionlouees
                  call sample(cellclass,nppmax,nu,ncellclass,listcell)


*     proposition de double precisionlocation dans la classe nclass+1
                  call split(nclass+1,nclass,c,ctemp,nppmax,nu,
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
            call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &           nppmax,nall,nallmax,z,n,indcell,c)
            call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &           nppmax,nall,nallmax,z,ntemp,indcell,ctemp)
c            write(*,*) 'apres count'
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)
            
*     proposition nouvelle freq et derive 
c            write(*,*) 'ajoutage des freq'
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)
            call addfreq7bis(nclass,nclassmax,nindiv,nloc,
     &           nlocmax,nlocmax2,nall,nallmax,nppmax,indcell,z,isplit,
     &           ctemp,f,ftemp,fa,drift,drifttemp,a,ptemp,ntemp)
c            write(*,*) 'z(67,9)=',z(67,9)
c            write(*,*) 'z(67,10)=',z(67,10)

*     calcul du log du ratio
*     terme des vraisemblances
c            write(*,*) 'calcul du log du ratio'
            lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &           indcell,nclassmax,nlocmax,nallmax,
     &           nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
c            write(*,*) 'lratio =',lratio

*     terme des freq.
            lratio = lratio + termfsplitbis(isplit,nclass,
     &           nclassmax,nlocmax,nlocmax2,nall,nallmax,nindiv,
     &           f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c            write(*,*) 'lratio =',lratio
*     terme des proposal sur c
            lratio = lratio + dlog(2*dble(ncellclass+1)) + 
     &           lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
c            write(*,*) 'lratio =',lratio
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(nclass)) - 
     &           dlog(dble(nclass+1)))
c            write(*,*) 'lratio =',lratio

            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
c$$$
c$$$                   write(*,*) 'nclass=',nclass
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellclass=',ncellclass  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctemp=',ctemp 
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                           n(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do iclass=1,nclassmax
c$$$                      do iloc=1,nlocmax
c$$$                         do iall=1,nall(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &                  indcell,nclassmax,nlocmax,nallmax,
c$$$     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplitbis(isplit,nclass,nclassmax,
c$$$     &                  nlocmax,nlocmax2,nall,nallmax,nindiv,
c$$$     &                  f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellclass+1)) + 
c$$$     &                  lbico(ncellclass,nu) - dlog(dble(nclass+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(nclass)) - 
c$$$     &                  dlog(dble(nclass+1)))
c$$$                   write(*,*) 'alpha=',alpha 
            
            if(bern .eq. 1) then
c               write(*,*) 'accept split'
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)

               call accept5(nppmax,nclassmax,nlocmax,
     &              nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
               nclass = nclass + 1
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
            endif
         endif

      else
         if(nclass .gt. nclassmin) then 
c      write(*,*) 'mort'
*     tirage de la classe qui meurt
            iclassrem = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            
*     tirage de la classe hote
            iclasshost = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            do while(iclasshost .eq. iclassrem)
               iclasshost = 1 
     &              + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la classe d'indice le plus petit
            if(iclasshost .gt. iclassrem) then
               ii = iclasshost
               iclasshost = iclassrem
               iclassrem = ii
            endif
            
*     recherche des cellules qui vont etre double precisionlouees
            call who(c,iclassrem,npp,nppmax,cellclass,ncellclass)
            
*     recherche des cellules de la classe hote
            call who(c,iclasshost,npp,nppmax,cellclasshost,
     &           ncellclasshost)

*     proposition de double precisionlocation dans la classe iclasshost
            call merging(iclassrem,iclasshost,nclass,c,ctemp,nppmax,
     &           ncellclass,cellclass)
            
*     comptage des alleles sur chaque locus pour c puis ctemp
            call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &           nppmax,nall,nallmax,z,n,indcell,c)
            call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &           nppmax,nall,nallmax,z,ntemp,indcell,ctemp)

*     propostion du nouveau tableau de freq et de derives
            call remfreq7bis(iclassrem,iclasshost,
     &           nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &           nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,
     &           ntemp)
            
*     calcul du log du ratio  
*     terme des vraisemblances
            lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &           indcell,nclassmax,nlocmax,nallmax,
     &           nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
c            write(*,*) 'lratio =',lratio
*     terme des freq.
            lratio = lratio + 
     &           termfmergebis(iclasshost,iclassrem,
     &           nclass,nclassmax,nlocmax,nlocmax2,
     &           nall,nallmax,nindiv,
     &           f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c            write(*,*) 'lratio =',lratio
*     terme des proposal sur c
            lratio = lratio + dlog(dble(nclass)) - 
     &           dlog(2*dble(ncellclass+ncellclasshost+1)) -
     &           lbico(ncellclass+ncellclasshost,ncellclass) 
c            write(*,*) 'lratio =',lratio
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(nclass)) - 
     &           dlog(dble(nclass-1)))
c            write(*,*) 'lratio =',lratio

            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha) 

c$$$               write(*,*) 'nclass=',nclass
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'iclassrem=',iclassrem
c$$$               write(*,*) 'cellclass=',cellclass
c$$$               write(*,*) 'iclasshost=',iclasshost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctemp=',ctemp 
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',iclass,',',iloc,',',iall,')=',
c$$$     &                       n(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ntemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',iclass,',',iloc,',',iall,')=',
c$$$     &                       f(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do iclass=1,nclassmax
c$$$                  do iloc=1,nlocmax
c$$$                     do iall=1,nall(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftemp(',iclass,',',iloc,',',iall,')=',
c$$$     &                           ftemp(iclass,iloc,iall)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(z,f,ftemp,c,ctemp,indcell,
c$$$     &              indcell,nclassmax,nlocmax,nallmax,
c$$$     &              nindiv,nindivmax,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmergebis(iclasshost,iclassrem,
c$$$     &              nclass,nclassmax,nlocmax,nlocmax2,
c$$$     &              nall,nallmax,nindiv,
c$$$     &              f,ftemp,n,ntemp,fa,drift,drifttemp,z)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(nclass)) - 
c$$$     &              dlog(2*dble(ncellclass+ncellclasshost+1)) -
c$$$     &              lbico(ncellclass+ncellclasshost,ncellclass) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(nclass)) - 
c$$$     &              dlog(dble(nclass-1)))
c$$$               write(*,*) 'alpha=',alpha 
  
            if(bern .eq. 1) then
c               write(*,*) 'accept merge'
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
               call accept5(nppmax,nclassmax,nlocmax,
     &              nallmax,nall,c,ctemp,f,ftemp,drift,drifttemp)
               nclass = nclass - 1
c               write(*,*) 'z(67,9)=',z(67,9)
c               write(*,*) 'z(67,10)=',z(67,10)
            endif
         endif
      endif 
      
      end subroutine bdclass8bis




**********************************************************************
*     ajoute une classe 
*     dans  le tableau des dérives selon le prior
*     et dans le tableau des frequences 
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq7(nclass,nclassmax,nindiv,nloc,nlocmax,nlocmax2,
     &     nall,nallmax,nppmax,indcell,z,isplit,ctemp,f,ftemp,fa,
     &     drift,drifttemp,a,ptemp,ntemp)
      implicit none
      integer nclass,nclassmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,ctemp(nppmax),ntemp(nclassmax,nloc,nallmax),
     &     nindiv,indcell(nindiv),z(nindiv,nlocmax2),isplit
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax),a(nallmax)
      integer iloc,iall,iclass
      double precision ptemp(nallmax),ggrunif


*     remplissage de f et drift pour les classes pre-existantes
      do iclass = 1,nclass
         drifttemp(iclass) = drift(iclass)
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttemp(isplit) = ggrunif(0.d0,1.d0)
      drifttemp(nclass+1) = ggrunif(0.d0,1.d0)

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(isplit))/
     &           drifttemp(isplit)+dble(ntemp(isplit,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(isplit,iloc,iall)  = ptemp(iall)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(nclass+1))/
     &           drifttemp(nclass+1)+dble(ntemp(nclass+1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(nclass+1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftemp(',isplit,1,1,')=',ftemp(isplit,1,1)
c$$$      write(*,*) 'ftemp(',isplit,1,2,')=',ftemp(isplit,1,2)
c$$$      write(*,*) 'f(',nclass+1,1,1,')=',f(nclass+1,1,1)
c$$$      write(*,*) 'f(',nclass+1,1,2,')=',f(nclass+1,1,2)
c$$$      write(*,*) 'ftemp(',nclass+1,1,1,')=',ftemp(nclass+1,1,1)
c$$$      write(*,*) 'ftemp(',nclass+1,1,2,')=',ftemp(nclass+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftemp(',isplit,2,1,')=',ftemp(isplit,2,1)
c$$$      write(*,*) 'ftemp(',isplit,2,2,')=',ftemp(isplit,2,2)
c$$$      write(*,*) 'f(',nclass+1,2,1,')=',f(nclass+1,2,1)
c$$$      write(*,*) 'f(',nclass+1,2,2,')=',f(nclass+1,2,2)
c$$$      write(*,*) 'ftemp(',nclass+1,2,1,')=',ftemp(nclass+1,2,1)
c$$$      write(*,*) 'ftemp(',nclass+1,2,2,')=',ftemp(nclass+1,2,2)
c$$$
      end subroutine addfreq7






***********************************************************************
*
*     ajoute une classe 
*     dans  le tableau des frequences  selon le prior
*     selon la conditionnelle complète pour les deux nouveaux groupes
*     et une valeur 0.5d0 dans le tableau des dérives 
*     (sans modifier les tableaux en entrée)
*     pour court-circuiter le F-model
      subroutine addfreq7bis(nclass,nclassmax,nindiv,nloc,nlocmax,
     &     nlocmax2,nall,nallmax,nppmax,indcell,z,isplit,ctemp,f,ftemp,
     &     fa,drift,drifttemp,a,ptemp,ntemp)
      implicit none
      integer nclass,nclassmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,ctemp(nppmax),ntemp(nclassmax,nloc,nallmax),
     &     nindiv,indcell(nindiv),z(nindiv,nlocmax2),isplit
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax),a(nallmax)
      integer iloc,iall,iclass
      double precision ptemp(nallmax)


*     remplissage de f et drift pour le classes pre-existantes
      do iclass = 1,nclass
         drifttemp(iclass) = drift(iclass)
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
            enddo 
         enddo
      enddo
 
*     nouvelles derives
      drifttemp(isplit) = 0.5d0
      drifttemp(nclass+1) = 0.5d0

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(isplit))/
     &           drifttemp(isplit)+dble(ntemp(isplit,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(isplit,iloc,iall)  = ptemp(iall)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(nclass+1))/
     &           drifttemp(nclass+1)+dble(ntemp(nclass+1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(nclass+1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftemp(',isplit,1,1,')=',ftemp(isplit,1,1)
c$$$      write(*,*) 'ftemp(',isplit,1,2,')=',ftemp(isplit,1,2)
c$$$      write(*,*) 'f(',nclass+1,1,1,')=',f(nclass+1,1,1)
c$$$      write(*,*) 'f(',nclass+1,1,2,')=',f(nclass+1,1,2)
c$$$      write(*,*) 'ftemp(',nclass+1,1,1,')=',ftemp(nclass+1,1,1)
c$$$      write(*,*) 'ftemp(',nclass+1,1,2,')=',ftemp(nclass+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftemp(',isplit,2,1,')=',ftemp(isplit,2,1)
c$$$      write(*,*) 'ftemp(',isplit,2,2,')=',ftemp(isplit,2,2)
c$$$      write(*,*) 'f(',nclass+1,2,1,')=',f(nclass+1,2,1)
c$$$      write(*,*) 'f(',nclass+1,2,2,')=',f(nclass+1,2,2)
c$$$      write(*,*) 'ftemp(',nclass+1,2,1,')=',ftemp(nclass+1,2,1)
c$$$      write(*,*) 'ftemp(',nclass+1,2,2,')=',ftemp(nclass+1,2,2)
c$$$

      end subroutine addfreq7bis




***********************************************************************
*     enleve une classe des tableau des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux classes
*     tirage d'une derive selon prior 
*     sans modifier des tableaux en entrée
      subroutine remfreq7(iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &     nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,ntemp)
      implicit none
      integer iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,c(nppmax),nindiv,z(nindiv,nlocmax2),
     &     ntemp(nclassmax,nlocmax,nallmax)
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax),a(nallmax) 
      integer iclass,iloc,iall
      double precision ptemp(nallmax),ggrunif

*     supprime la classe iclassrem
      if(iclassrem .eq. 1) then
         do iclass =2,nclass
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass-1,iloc,iall) = f(iclass,iloc,iall)
               enddo 
            enddo
            drifttemp(iclass-1) = drift(iclass)
         enddo
      else
         do iclass =1,iclassrem-1
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
               enddo 
            enddo
            drifttemp(iclass) = drift(iclass)
         enddo
         if(iclassrem .ne. nclass) then
            do iclass =iclassrem+1,nclass
               do iloc=1,nloc
                  do iall=1,nall(iloc)
                     ftemp(iclass-1,iloc,iall) = f(iclass,iloc,iall)
                  enddo 
               enddo
               drifttemp(iclass-1) = drift(iclass)
            enddo
         endif
      endif
      do iclass=nclass,nclassmax
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = -999
            enddo 
         enddo
         drifttemp(iclass) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle classe
      drifttemp(iclasshost) = ggrunif(0.d0,1.d0)
      
*     frequences pour la nouvelle classe
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(iclasshost))/
     &           drifttemp(iclasshost)+
     &           dble(ntemp(iclasshost,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(iclasshost,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',iclasshost,1,1,')=',f(iclasshost,1,1)
c$$$      write(*,*) 'f(',iclasshost,1,2,')=',f(iclasshost,1,2)
c$$$      write(*,*) 'f(',iclasshost,2,1,')=',f(iclasshost,2,1)
c$$$      write(*,*) 'f(',iclasshost,2,2,')=',f(iclasshost,2,2)
c$$$
c$$$      write(*,*) 'ftemp(',iclasshost,1,1,')=',ftemp(iclasshost,1,1)
c$$$      write(*,*) 'ftemp(',iclasshost,1,2,')=',ftemp(iclasshost,1,2)
c$$$      write(*,*) 'ftemp(',iclasshost,2,1,')=',ftemp(iclasshost,2,1)
c$$$      write(*,*) 'ftemp(',iclasshost,2,2,')=',ftemp(iclasshost,2,2)

      end subroutine remfreq7





***********************************************************************
*     enleve une classe des tableaux des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux classes
*     la nouvelle derive est mise à 0.5d0
*     (sans modifier les tableaux en entrée)
*     c'est pour court-circuiter le F-model 
      subroutine remfreq7bis(iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,nppmax,
     &     nallmax,f,ftemp,drift,drifttemp,c,nindiv,z,fa,a,ptemp,ntemp)
      implicit none
      integer iclassrem,iclasshost,
     &     nclass,nclassmax,nloc,nlocmax,nlocmax2,nall(nlocmax),
     &     nallmax,nppmax,c(nppmax),nindiv,z(nindiv,nlocmax2),
     &     ntemp(nclassmax,nlocmax,nallmax)
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax),a(nallmax) 
      integer iclass,iloc,iall
      double precision ptemp(nallmax)

*     supprime la classe iclassrem
      if(iclassrem .eq. 1) then
         do iclass =2,nclass
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass-1,iloc,iall) = f(iclass,iloc,iall)
               enddo 
            enddo
            drifttemp(iclass-1) = drift(iclass)
         enddo
      else
         do iclass =1,iclassrem-1
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
               enddo 
            enddo
            drifttemp(iclass) = drift(iclass)
         enddo
         if(iclassrem .ne. nclass) then
            do iclass =iclassrem+1,nclass
               do iloc=1,nloc
                  do iall=1,nall(iloc)
                     ftemp(iclass-1,iloc,iall) = f(iclass,iloc,iall)
                  enddo 
               enddo
               drifttemp(iclass-1) = drift(iclass)
            enddo
         endif
      endif
      do iclass=nclass,nclassmax
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = -999
            enddo 
         enddo
         drifttemp(iclass) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle classe
      drifttemp(iclasshost) = 0.5d0
      
*     frequences pour la nouvelle classe
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drifttemp(iclasshost))/
     &           drifttemp(iclasshost)+
     &           dble(ntemp(iclasshost,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(iclasshost,iloc,iall)  = ptemp(iall)
         enddo
      enddo

c$$$      write(*,*) 'f(',iclasshost,1,1,')=',f(iclasshost,1,1)
c$$$      write(*,*) 'f(',iclasshost,1,2,')=',f(iclasshost,1,2)
c$$$      write(*,*) 'f(',iclasshost,2,1,')=',f(iclasshost,2,1)
c$$$      write(*,*) 'f(',iclasshost,2,2,')=',f(iclasshost,2,2)
c$$$
c$$$      write(*,*) 'ftemp(',iclasshost,1,1,')=',ftemp(iclasshost,1,1)
c$$$      write(*,*) 'ftemp(',iclasshost,1,2,')=',ftemp(iclasshost,1,2)
c$$$      write(*,*) 'ftemp(',iclasshost,2,1,')=',ftemp(iclasshost,2,1)
c$$$      write(*,*) 'ftemp(',iclasshost,2,2,')=',ftemp(iclasshost,2,2)

      end subroutine remfreq7bis




***********************************************************************
*
*     terme des freq dans le log ratio pour un split 
      double precision function termfsplit(isplit,nclass,nclassmax,
     &     nlocmax,nlocmax2,nall,nallmax,nindiv,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp,z)
      implicit none
      integer nclassmax,nlocmax,nlocmax2,nallmax,nindiv,
     &     n(nclassmax,nlocmax,nallmax),
     &     ntemp(nclassmax,nlocmax,nallmax),nall(nlocmax),
     &     z(nindiv,nlocmax2),isplit,nclass
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      double precision q,gglgammafn,ss,tt

      
      termfsplit  = 0.d0
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(isplit,iloc,iall))) +
     &           dble(n(isplit,iloc,iall))*
     &           dlog(f(isplit,iloc,iall))
            nn = nn + n(isplit,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfsplit = termfsplit + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttemp(isplit))/drifttemp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(isplit,iloc,iall))) +
     &           dble(ntemp(isplit,iloc,iall))*
     &           dlog(ftemp(isplit,iloc,iall))
            nn = nn + ntemp(isplit,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfsplit = termfsplit - tt

*    ln( pi[ f_{nclass+1}^* |fa, drift] / pi[f_{nclass+1}^*| ...] ) 
      tt = 0
      q = (1-drifttemp(nclass+1))/drifttemp(nclass+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(nclass+1,iloc,iall))) +
     &           dble(ntemp(nclass+1,iloc,iall))*
     &           dlog(ftemp(nclass+1,iloc,iall))
            nn = nn + ntemp(nclass+1,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfsplit = termfsplit - tt



      end function termfsplit




***********************************************************************
*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(iclass) =0.5d0
      double precision function termfsplitbis(isplit,nclass,nclassmax,
     &     nlocmax,nlocmax2,nall,nallmax,nindiv,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp,z)
      implicit none
      integer nclassmax,nlocmax,nlocmax2,nallmax,nindiv,
     &     n(nclassmax,nlocmax,nallmax),
     &     ntemp(nclassmax,nlocmax,nallmax),nall(nlocmax),
     &     z(nindiv,nlocmax2),isplit,nclass
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      double precision q,gglgammafn,ss,tt

      
      termfsplitbis  = 0.d0
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(isplit,iloc,iall))) +
     &           dble(n(isplit,iloc,iall))*
     &           dlog(f(isplit,iloc,iall))
            nn = nn + n(isplit,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttemp(isplit))/drifttemp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(isplit,iloc,iall))) +
     &           dble(ntemp(isplit,iloc,iall))*
     &           dlog(ftemp(isplit,iloc,iall))
            nn = nn + ntemp(isplit,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

*    ln( pi[ f_{nclass+1}^* |fa, drift] / pi[f_{nclass+1}^*| ...] ) 
      tt = 0
      q = (1-drifttemp(nclass+1))/drifttemp(nclass+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(nclass+1,iloc,iall))) +
     &           dble(ntemp(nclass+1,iloc,iall))*
     &           dlog(ftemp(nclass+1,iloc,iall))
            nn = nn + ntemp(nclass+1,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

      end function termfsplitbis





***********************************************************************
*
*     terme des freq dans le log ratio pour un split 
      double precision function termfmerge(ihost,irem,nclass,nclassmax,
     &     nlocmax,nlocmax2,nall,nallmax,nindiv,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp,z)
      implicit none
      integer nclassmax,nlocmax,nlocmax2,nallmax,nindiv,
     &     n(nclassmax,nlocmax,nallmax),
     &     ntemp(nclassmax,nlocmax,nallmax),nall(nlocmax),
     &     z(nindiv,nlocmax2),ihost,irem,nclass
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      double precision q,gglgammafn,ss,tt

      
      termfmerge  = 0.d0
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(ihost,iloc,iall))) +
     &           dble(n(ihost,iloc,iall))*
     &           dlog(f(ihost,iloc,iall))
            nn = nn + n(ihost,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(irem,iloc,iall))) +
     &           dble(n(irem,iloc,iall))*
     &           dlog(f(irem,iloc,iall))
            nn = nn + n(irem,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttemp(ihost))/drifttemp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(ihost,iloc,iall))) +
     &           dble(ntemp(ihost,iloc,iall))*
     &           dlog(ftemp(ihost,iloc,iall))
            nn = nn + ntemp(ihost,iloc,iall)
         enddo
         tt = gglgammafn(q + dble(nn)) - gglgammafn(q) + ss
      enddo
      termfmerge = termfmerge - tt
      end function termfmerge





***********************************************************************
*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(iclass) =0.5d0
      double precision function termfmergebis(ihost,irem,nclass,
     &     nclassmax,nlocmax,nlocmax2,nall,nallmax,nindiv,
     &     f,ftemp,n,ntemp,fa,drift,drifttemp,z)
      implicit none
      integer nclassmax,nlocmax,nlocmax2,nallmax,nindiv,
     &     n(nclassmax,nlocmax,nallmax),
     &     ntemp(nclassmax,nlocmax,nallmax),nall(nlocmax),
     &     z(nindiv,nlocmax2),ihost,irem,nclass
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     fa(nlocmax,nallmax)
      integer iloc, iall,nn
      double precision q,gglgammafn,ss,tt

      
      termfmergebis  = 0.d0
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(ihost,iloc,iall))) +
     &           dble(n(ihost,iloc,iall))*
     &           dlog(f(ihost,iloc,iall))
            nn = nn + n(ihost,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(n(irem,iloc,iall))) +
     &           dble(n(irem,iloc,iall))*
     &           dlog(f(irem,iloc,iall))
            nn = nn + n(irem,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttemp(ihost))/drifttemp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do iall = 1,nall(iloc)
            ss = ss + gglgammafn(fa(iloc,iall) * q) -
     &           gglgammafn(fa(iloc,iall) * q +
     &           dble(ntemp(ihost,iloc,iall))) +
     &           dble(ntemp(ihost,iloc,iall))*
     &           dlog(ftemp(ihost,iloc,iall))
            nn = nn + ntemp(ihost,iloc,iall)
         enddo
         tt = gglgammafn(dble(nall(iloc)) + dble(nn)) - 
     &        gglgammafn(dble(nall(iloc))) + ss
      enddo
      termfmergebis = termfmergebis - tt
      end function termfmergebis





***********************************************************************
*
*     log vraisemblance
*
      double precision function ll(z,nindivmax,nlocmax,nlocmax2,nall,
     &     nclassmax,nallmax,nppmax,c,f,indcell)
      implicit none
      integer nindivmax,nlocmax,nlocmax2,nall(nlocmax),nclassmax,
     &     z(nindivmax,nlocmax2),nppmax,c(nppmax),nallmax,
     &     indcell(nindivmax)
      double precision f(nclassmax,nlocmax,nallmax)
      integer iindiv,iloc,iall1,iall2,iclass

      ll = 0
      do iindiv = 1,nindivmax
         iclass = c(indcell(iindiv))
         do iloc = 1,nlocmax
            iall1 = z(iindiv,2*iloc-1)
            iall2 = z(iindiv,2*iloc)
            if(iall1 .ne. -999) then 
               ll = ll + dlog(f(iclass,iloc,iall1))
            endif
            if(iall2 .ne. -999) then 
               ll = ll + dlog(f(iclass,iloc,iall2))
               if((iall1 .ne. iall2) .and. (iall1 .ne. -999))  then 
                  ll = ll + dlog(2.d0)
               endif
            endif
         enddo
      enddo
      end function ll



*********************************************************************** 
*
*     log de la proba a posteriori du vecteur de parametres
*
      double precision function lpp(lambda,z,nclass,npp,drift,f,fa,c,
     &     nppmax,nindivmax,nlocmax2,nclassmax,nlocmax,nallmax,indcell,
     &     nall,fmodel)
      implicit none
      integer nindivmax,nlocmax2,nclass,nclassmax,nlocmax,nallmax,
     &     npp,nppmax,z(nindivmax,nlocmax2),indcell(nindivmax),
     &     c(nppmax),nall(nlocmax),fmodel
      double precision drift(nclassmax),f(nclassmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),lambda
      integer ipp,iclass,iloc,iall
      double precision gglgammafn,shape1,shape2,ll

      parameter(shape1=2.d0,shape2=20.d0)

      lpp = -lambda + npp*dlog(lambda) 
      do ipp=1,npp
         lpp = lpp - dlog(dble(ipp))
      enddo
      lpp = lpp -npp*dlog(dble(nclass))
      
      if(fmodel .eq. 1) then
         do iclass = 1,nclass
            lpp = lpp + 
     &           gglgammafn(shape1+shape2) - gglgammafn(shape1) - 
     &           gglgammafn(shape2)
     &           + (shape1-1)*dlog(drift(iclass)) + 
     &           (shape2-1)*dlog(1-drift(iclass))
            do iloc = 1,nlocmax
               lpp = lpp + gglgammafn((1-drift(iclass))/drift(iclass))
               do iall=1,nall(iloc)
                  lpp = lpp -gglgammafn(fa(iloc,iall)*
     &                 (1-drift(iclass))/drift(iclass)) + 
     &                 (fa(iloc,iall)* 
     &                 (1-drift(iclass))/drift(iclass)-1)*
     &              f(iclass,iloc,iall)
               enddo
            enddo
         enddo
      endif
      
      lpp = lpp + ll(z,nindivmax,nlocmax,nlocmax2,nall,nclassmax,
     &     nallmax,nppmax,c,f,indcell)

      end function lpp






***********************************************************************
*     joint update of c and f under the Dirichlet model
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclass,ipp,isplit,iclass1,iclass2,iloc,iall
      double precision alpha,ggrunif,lrpf,lratio,llr6,b,ggrbinom,bern

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
      do iclass = 1,nclass
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
            enddo 
         enddo
      enddo

      do ipp=1,npp
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctemp(ipp) = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
         iclass1 = c(ipp)
         iclass2 = ctemp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &        nppmax,nall,nallmax,z,n,indcell,c)
         call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &        nppmax,nall,nallmax,z,ntemp,indcell,ctemp)


c         write(*,*) 'n=',n
c         write(*,*) 'ntemp=',ntemp

*     sample new frequencies
         call samplef(nclass,nclassmax,nindiv,nloc,nlocmax,nlocmax2,
     &     nall,nallmax,nppmax,z,iclass1,iclass2,ctemp,f,ftemp,
     &     fa,drift,a,ptemp,ntemp)

c         write(*,*) 'f=',f
c         write(*,*) 'ftemp=',ftemp



*     compute M-H ratio
*     likelihood
         lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &                  indcell,nclassmax,nlocmax,nallmax,
     &                  nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)

c         write(*,*) 'lratio =',lratio

*     contrib proposals
         lratio = lratio + lrpf(nclassmax,nloc,nall,nallmax,n,ntemp,
     &        f,ftemp,iclass1,iclass2)

c          write(*,*) 'lratio =',lratio

         lratio = dmin1(0.d0,lratio)
         alpha = dexp(lratio)

c         write(*,*) 'alpha=',alpha

         bern = ggrbinom(1.d0,alpha)

c         write(*,*) 'bern=',bern

         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  f(iclass1,iloc,iall) = ftemp(iclass1,iloc,iall)
                  f(iclass2,iloc,iall) = ftemp(iclass2,iloc,iall)
               enddo
            enddo
         else 
            ctemp(ipp) = c(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass1,iloc,iall) = f(iclass1,iloc,iall)
                  ftemp(iclass2,iloc,iall) = f(iclass2,iloc,iall)
               enddo
            enddo
         endif
      enddo
      end subroutine udcf


***********************************************************************
*     joint update of c and f under the Falush model
*     single component update of c
*     new f is proposed according to full conditionnal pi(f*|c*,z)
      subroutine udcf2(nclass,nclassmin,nclassmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nall,nallmax,indcell,nindiv,nindivmax,npp,nppmax,c,ctemp,
     &     a,ptemp,ftemp,drifttemp,z,n,ntemp,ploidy)
      implicit none
      integer nclass,nclassmin,nclassmax,nloc,nlocmax,nall(nlocmax),
     &     nallmax,nindiv,nindivmax,npp,nppmax,indcell(nindivmax),
     &     nlocmax2,c(nppmax),ctemp(nppmax),z(nindivmax,nlocmax2),
     &     n(nclassmax,nloc,nallmax),ntemp(nclassmax,nloc,nallmax),
     &     ploidy
      double precision f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &      ftemp(nclassmax,nlocmax,nallmax),drifttemp(nclassmax),
     &     a(nallmax),ptemp(nallmax),fa(nlocmax,nallmax)
      integer iclass,ipp,isplit,iclass1,iclass2,iloc,iall
      double precision alpha,ggrunif,lrppf,lratio,llr6,b,ggrbinom,bern


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
      do iclass = 1,nclass
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
            enddo 
         enddo
      enddo

      do ipp=1,npp
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctemp(ipp) = 1 + int(dint(dble(nclass)*ggrunif(0.d0,1.d0)))
         iclass1 = c(ipp)
         iclass2 = ctemp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &        nppmax,nall,nallmax,z,n,indcell,c)
         call countn(nindiv,nlocmax,nlocmax2,nclass,nclassmax,
     &        nppmax,nall,nallmax,z,ntemp,indcell,ctemp)
         
*     sample new frequencies
         call samplef(nclass,nclassmax,nindiv,nloc,nlocmax,nlocmax2,
     &        nall,nallmax,nppmax,z,iclass1,iclass2,ctemp,f,ftemp,
     &        fa,drift,a,ptemp,ntemp)
         
*     compute M-H ratio
*     likelihood
         lratio =  llr6(z,f,ftemp,c,ctemp,indcell,
     &        indcell,nclassmax,nlocmax,nallmax,
     &        nindiv,nindivmax,nloc,nlocmax2,nppmax,ploidy)
         
*     contrib prior and proposals
         lratio = lratio + lrppf(nclassmax,nloc,nall,nallmax,n,
     &     ntemp,fa,drift,f,ftemp,iclass1,iclass2)

         lratio = dmin1(0.d0,lratio)
         alpha = dexp(lratio)
         bern = ggrbinom(1.d0,alpha)
         
c         write(*,*) 'bern=',bern

         if(bern .eq. 1) then
            c(ipp) = ctemp(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  f(iclass1,iloc,iall) = ftemp(iclass1,iloc,iall)
                  f(iclass2,iloc,iall) = ftemp(iclass2,iloc,iall)
               enddo
            enddo
         else 
            ctemp(ipp) = c(ipp)
            do iloc=1,nloc
               do iall=1,nall(iloc)
                  ftemp(iclass1,iloc,iall) = f(iclass1,iloc,iall)
                  ftemp(iclass2,iloc,iall) = f(iclass2,iloc,iall)
               enddo
            enddo
         endif        
      enddo
      end subroutine udcf2





***********************************************************************
*
*     log of ratio of proposals in a joint update of c and f
*     (in subroutine udcf)
*     
      double precision function lrpf(nclassmax,nloc,nall,nallmax,n,
     &     ntemp,f,ftemp,iclass1,iclass2)
      implicit none 
      integer nclassmax,nloc,nall,nallmax,n,ntemp,iclass1,iclass2
      double precision f,ftemp
      dimension nall(nloc),n(nclassmax,nloc,nallmax),
     &     ntemp(nclassmax,nloc,nallmax),f(nclassmax,nloc,nallmax),
     &     ftemp(nclassmax,nloc,nallmax)
      integer iloc,iall,n1,n2,ntemp1,ntemp2
      double precision gglgammafn
      lrpf = 0.d0
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntemp1 = 0
         ntemp2 = 0
         do iall=1,nall(iloc)
            lrpf = lrpf + 
     &         n(iclass1,iloc,iall)*dlog(f(iclass1,iloc,iall)) 
     &       + n(iclass2,iloc,iall)*dlog(f(iclass2,iloc,iall) )
     &       - ntemp(iclass1,iloc,iall)*dlog(ftemp(iclass1,iloc,iall)) 
     &       - ntemp(iclass2,iloc,iall)*dlog(ftemp(iclass2,iloc,iall)) 
     &       - gglgammafn(1+dble(n(iclass1,iloc,iall)))
     &       - gglgammafn(1+dble(n(iclass2,iloc,iall)))
     &       + gglgammafn(1+dble(ntemp(iclass1,iloc,iall)))
     &       + gglgammafn(1+dble(ntemp(iclass2,iloc,iall)))
            n1 = n1 + n(iclass1,iloc,iall)
            n2 = n2 + n(iclass2,iloc,iall)
            ntemp1 = ntemp1 + ntemp(iclass1,iloc,iall)
            ntemp2 = ntemp2 + ntemp(iclass2,iloc,iall)
         enddo
         lrpf = lrpf + gglgammafn(dble(nall(iloc)+n1)) 
     &           + gglgammafn(dble(nall(iloc)+n2)) 
     &           - gglgammafn(dble(nall(iloc)+ntemp1)) 
     &           - gglgammafn(dble(nall(iloc)+ntemp2))
      enddo
      end 





***********************************************************************
*
*     log of ratio of proposals x priors in a joint update of c and f
*     (in subroutine udcf2)
*     
      double precision function lrppf(nclassmax,nloc,nall,nallmax,n,
     &     ntemp,fa,drift,f,ftemp,iclass1,iclass2)
      implicit none 
      integer nclassmax,nloc,nall,nallmax,n,ntemp,iclass1,iclass2
      double precision fa,drift,f,ftemp
      dimension nall(nloc),n(nclassmax,nloc,nallmax),
     &     ntemp(nclassmax,nloc,nallmax),f(nclassmax,nloc,nallmax),
     &     fa(nloc,nallmax),ftemp(nclassmax,nloc,nallmax),
     &     drift(nclassmax)
      integer iloc,iall,n1,n2,ntemp1,ntemp2
      double precision gglgammafn,q1,q2
      lrppf = 0.d0
      q1 = drift(iclass1)/(1-drift(iclass1))
      q2 = drift(iclass2)/(1-drift(iclass2))
      
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntemp1 = 0
         ntemp2 = 0
         do iall=1,nall(iloc)
            lrppf = lrppf + 
     &         n(iclass1,iloc,iall)*dlog(f(iclass1,iloc,iall)) 
     &       + n(iclass2,iloc,iall)*dlog(f(iclass2,iloc,iall) )
     &       - ntemp(iclass1,iloc,iall)*dlog(ftemp(iclass1,iloc,iall)) 
     &       - ntemp(iclass2,iloc,iall)*dlog(ftemp(iclass2,iloc,iall)) 
     &       + gglgammafn(fa(iloc,iall)*q1+
     &           dble(ntemp(iclass1,iloc,iall)))
     &       + gglgammafn(fa(iloc,iall)*q2+
     &           dble(ntemp(iclass2,iloc,iall)))
     &       - gglgammafn(fa(iloc,iall)*dble(q1)+
     &        dble(n(iclass1,iloc,iall)))
     &       - gglgammafn(fa(iloc,iall)*dble(q2)+
     &        dble(n(iclass2,iloc,iall)))
            n1 = n1 + n(iclass1,iloc,iall)
            n2 = n2 + n(iclass2,iloc,iall)
            ntemp1 = ntemp1 + ntemp(iclass1,iloc,iall)
            ntemp2 = ntemp2 + ntemp(iclass2,iloc,iall)
         enddo
         lrppf = lrppf + gglgammafn(q1+dble(n1)) 
     &           + gglgammafn(q2+dble(n2)) 
     &           - gglgammafn(q1+dble(ntemp1)) 
     &           - gglgammafn(q2+dble(ntemp2))
      enddo
      end 





***********************************************************************
*     
*     sample freq from full contitionnal pi(f|u,c,z)
*     drifts and ancestral are not changed
*
      subroutine samplef(nclass,nclassmax,nindiv,nloc,nlocmax,nlocmax2,
     &     nall,nallmax,nppmax,z,iclass1,iclass2,ctemp,f,ftemp,
     &     fa,drift,a,ptemp,ntemp)
      implicit none
      integer nclass,nclassmax,nloc,nlocmax,nlocmax2,nall,
     &     nallmax,nppmax,ctemp,ntemp,
     &     nindiv,z,iclass1,iclass2
      double precision f,drift,ftemp,fa,a
      dimension nall(nlocmax),ctemp(nppmax),
     &     ntemp(nclassmax,nloc,nallmax),
     &     z(nindiv,nlocmax2),
     &     f(nclassmax,nlocmax,nallmax),drift(nclassmax),
     &     ftemp(nclassmax,nlocmax,nallmax),
     &     fa(nlocmax,nallmax),a(nallmax),ptemp(nallmax)
      integer iloc,iall,iclass
      double precision ptemp,ggrunif


*     init ftemp 
      do iclass = 1,nclass
         do iloc=1,nloc
            do iall=1,nall(iloc)
               ftemp(iclass,iloc,iall) = f(iclass,iloc,iall)
            enddo 
         enddo
      enddo


*     new freq
*     for pop iclass1
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drift(iclass1))/
     &           drift(iclass1)+dble(ntemp(iclass1,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(iclass1,iloc,iall)  = ptemp(iall)
         enddo
      enddo

*     for pop iclass2
      do iloc=1,nloc
         do iall = 1,nall(iloc)
            a(iall) = fa(iloc,iall)*(1-drift(iclass2))/
     &           drift(iclass2)+dble(ntemp(iclass2,iloc,iall))
         enddo
         call dirichlet(nall(iloc),nallmax,a,ptemp)
         do iall=1,nall(iloc)
            ftemp(iclass2,iloc,iall)  = ptemp(iall)
         enddo
      enddo
      end subroutine  samplef




***********************************************************************
*
*     compute posterior probabilities of pop. membership
*     for pixels and individuals 
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
      double precision s,u,xlim(2),ylim(2),coorddom,dom,domperm,
     &     distvois,f,dt

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
      call intpr('*****************************************',-1,'*',0)
      call intpr('*  Computing posterior probabilities     ',-1,'*',0)
      call intpr('*  of population membership for pixels   ',-1,'*',0)
      call intpr('*****************************************',-1,'*',0)
      call intpr(' ',-1,':',0) 
      call intpr(' ',-1,':',0) 

      call limit(nindiv,nindivmax,s,xlim,ylim,dt)


*     coordonnées de la grille 
      idom = 1
      do ixdom =1,nxdommax
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nydommax
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           dble(ixdom-1)*(xlim(2) - xlim(1))/dble(nxdommax-1)
            coorddom(2,idom) = ylim(1) +
     &           dble(iydom-1)*(ylim(2) - ylim(1))/dble(nydommax-1)
            do iclass=1,nclassmax
               dom(idom,iclass) = 0.d0
               domperm(idom,iclass) = 0.d0
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
c         write(6,100)dble(ichain)/dble(nchain)*100.
         
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
               dom(idom,iclass) =  dom(idom,iclass) + 1.d0
               domperm(idom,iclassperm) = domperm(idom,iclassperm)+ 1.d0
            enddo
         endif
      enddo


      do idom=1,nxdommax*nydommax
         do iclass=1,nclassmax
            dom(idom,iclass) = dom(idom,iclass)/dble(nchain-burn)
            domperm(idom,iclass) = domperm(idom,iclass)/
     &           dble(nchain-burn)
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



********************************************************************
*     posterior probability of population membership for individuals
*
      subroutine  pppmindiv(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,filenpp,fileu,filec,nit,
     &     burn)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npp,c,
     &     nit,burn
      double precision pmp,distcell,u,s

      integer iit,ipp,iindiv,nppcur,ccur,ipop
      double precision xlim(2),ylim(2),ucur,dt
      character*256 files,fileu,filec,filenpp,
     &     filenclass,filedom,filedomperm,filenall,filef,filefperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)

      call intpr('*******************************************',-1,'*',0)
      call intpr('*  Computing posterior probabilities       ',-1,'*',0)
      call intpr('*  of population membership for individuals',-1,'*',0)
      call intpr('*******************************************',-1,'*',0)


      
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)

*     sequentially processes states of the chain
      do iit=1,nit
10000    format(f7.3,' %')
c         write(6,10000)dble(iit)/dble(nit)*100.
         
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         
         if(iit .gt. burn)  then
            call calccell(nindiv,nindiv,s,npp,nppmax,u,indcell,distcell)
            do iindiv=1,nindiv
               ipop = c(indcell(iindiv))
               pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.d0
            enddo
         endif 
      enddo
      
c      write(*,*) 'pmp=',pmp
      do iindiv=1,nindiv 
         do ipop=1,npopmax
            pmp(iindiv,ipop) = pmp(iindiv,ipop)/dble(nit-burn)
         enddo
      enddo
c      write(*,*) 'pmp=',pmp
      close(10)
      close(11)
      close(12)
      end subroutine  pppmindiv


      
 
*     Calcul du Fst d'apres programme en Turbo Pascal d'Arnaud Estoup
*     le vecteur c contient la variable de classe
*     les effectifs des classes doivent etre donnees en entree

      subroutine fstat(nindiv,nppmax,nloc,nloc2,nall,nclass,effcl,z,c,
     &     tabindiv,kk,Fistot,Fsttot,Fittot)
      implicit none
      integer nindiv,nppmax,nloc,nloc2,nall(nloc),nclass,
     &     z(nindiv,nloc2),c(nppmax),effcl(nclass)
      real Fsttot,Fittot,Fistot
      integer iloc,iclass,iall,iindiv
      integer g1,g2,k,tabindiv(nindiv,nclass),kk(nclass)
      real s1,s2,s3,s1l,s2l,s3l,ni,sni,sni2,sniA,sniAA,s2A,nA,AA,
     &     nc,MSG,MSI,MSP,s2G,s2I,s2P,Fst,Fit,Fis

*     Recherche des indices des indiv de chaque classe 
      do iclass=1,nclass
         kk(iclass) = 1
      enddo
      do iindiv=1,nindiv
         tabindiv(kk(c(iindiv)),c(iindiv)) = iindiv
         kk(c(iindiv)) = kk(c(iindiv)) + 1
      enddo
c      write(*,*) (tabindiv(iindiv,1),iindiv=1,nindiv)
c      write(*,*) (tabindiv(iindiv,2),iindiv=1,nindiv)
       
      s1 = 0.
      s2 = 0.
      s3 = 0.
      do iloc=1,nloc
c         write(*,*) 'iloc=', iloc
         s1l = 0.
         s2l = 0.
         s3l = 0.
         do iall=1,nall(iloc)
c            write(*,*) 'iall=', iall
            sni = 0.
            sni2 = 0.
            sniA = 0.
            sniAA = 0.
            s2A = 0.
c            k = 1
            do iclass=1,nclass
c               write(*,*) 'iclass=',iclass
               ni = effcl(iclass)
c                write(*,*) 'ni=',ni
               nA = 0.
               AA = 0.
c               do iindiv=k,(k+ni-1)
               do iindiv=1,ni 
                  k = tabindiv(iindiv,iclass)
c                  write(*,*) 'iindiv=',iindiv
                  g1 = z(k,2*(iloc-1)+1)
                  g2 = z(k,2*(iloc-1)+2)
                  if((g1 .eq. iall) .and. (g2 .eq. iall)) then 
                     AA = AA + 1.
                  endif
                  if(g1 .eq. iall) then 
                     nA = nA + 1. 
                  endif
                  if(g2 .eq. iall) then 
                     nA = nA + 1. 
                  endif
c                  write(*,*) 'variables=',AA,nA
               enddo
c               k = k+ni
               sniA = sniA + nA
               sniAA = sniAA + AA
               sni = sni + ni 
               sni2 = sni2 + ni*ni
               s2A = s2A + nA*nA/(2*ni)
c               write(*,*) 'variables=',sniA,sniAA,sni,sni2,s2A
            enddo
            nc = (sni-sni2/sni)/(nclass-1)
            MSG = (0.5*sniA-sniAA)/sni
            MSI = (0.5*sniA+sniAA-s2A)/(sni-nclass)
            MSP = (s2A-0.5*(sniA**2)/sni)/(nclass-1)
            s2G = MSG
            s2I = 0.5*(MSI-MSG)
            s2P = (MSP-MSI)/(2*nc)
            s1l = s1l + s2P
            s2l = s2l + s2P + s2I
            s3l = s3l + s2P + s2I + s2G
c            write(*,*) 'variables=',nc,MSG,MSI,s2G,s2I,s2P,s1l,s2l,s3l
         enddo
         Fst = s1l/s3l
         Fit = s2l/s3l
         Fis = (Fit-Fst)/(1-Fst)
         s1 = s1 + s1l
         s2 = s2 + s2l
         s3 = s3 + s3l
c         write(*,*) 'variables=',Fst,Fit,Fis,s1,s2,s3
      enddo
      Fsttot = s1/s3
      Fittot = S2/s3
      Fistot = (Fittot-Fsttot)/(1-Fsttot)
      end
      




































 



