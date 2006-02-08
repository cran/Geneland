*     Calcul du Fst d'apres programme en Turbo Pascal d'Arnaud Estoup
*     le vecteur c contient la variable de classe
*     les effectifs des classes doivent etre donnees en entree

      subroutine fstae(nindiv,nppmax,nloc,nloc2,nall,npop,effcl,z,c,
     &     tabindiv,kk,Fistot,Fsttot,Fittot)
      implicit none
      integer nindiv,nppmax,nloc,nloc2,nall(nloc),npop,
     &     z(nindiv,nloc2),c(nppmax),effcl(npop)
      real Fsttot,Fittot,Fistot
      integer iloc,ipop,iall,iindiv
      integer g1,g2,k,tabindiv(nindiv,npop),kk(npop)
      real s1,s2,s3,s1l,s2l,s3l,ni,sni,sni2,sniA,sniAA,s2A,nA,AA,
     &     nc,MSG,MSI,MSP,s2G,s2I,s2P,Fst,Fit,Fis

*     Recherche des indices des indiv de chaque pope 
      do ipop=1,npop
         kk(ipop) = 1
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
            do ipop=1,npop
c               write(*,*) 'ipop=',ipop
               ni = effcl(ipop)
c                write(*,*) 'ni=',ni
               nA = 0.
               AA = 0.
c               do iindiv=k,(k+ni-1)
               do iindiv=1,ni 
                  k = tabindiv(iindiv,ipop)
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
            nc = (sni-sni2/sni)/(npop-1)
            MSG = (0.5*sniA-sniAA)/sni
            MSI = (0.5*sniA+sniAA-s2A)/(sni-npop)
            MSP = (s2A-0.5*(sniA**2)/sni)/(npop-1)
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
      




































