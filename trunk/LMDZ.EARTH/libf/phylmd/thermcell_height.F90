      SUBROUTINE thermcell_height(ngrid,nlay,lalim,lmin,linter,lmix,  &
     &           zw2,zlev,lmax,zmax,zmax0,zmix,wmax,lev_out)                            

!-----------------------------------------------------------------------------
!thermcell_height: calcul des caracteristiques du thermique: zmax,wmax,zmix
!-----------------------------------------------------------------------------
      IMPLICIT NONE
#include "iniprint.h"
#include "thermcell.h"

      INTEGER ig,l
      INTEGER ngrid,nlay
      INTEGER lalim(ngrid),lmin(ngrid)
      INTEGER lmix(ngrid)
      REAL linter(ngrid)
      integer lev_out                           ! niveau pour les print

      REAL zw2(ngrid,nlay+1)
      REAL zlev(ngrid,nlay+1)

      REAL wmax(ngrid)
      INTEGER lmax(ngrid)
      REAL zmax(ngrid)
      REAL zmax0(ngrid)
      REAL zmix(ngrid)
      REAL num(ngrid)
      REAL denom(ngrid)

      REAL zlevinter(ngrid)

!calcul de la hauteur max du thermique
      do ig=1,ngrid
         lmax(ig)=lalim(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lalim(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo

! On traite le cas particulier qu'il faudrait éviter ou le thermique
! atteind le haut du modele ...
      do ig=1,ngrid
      if ( zw2(ig,nlay) > 1.e-10 ) then
          print*,'WARNING !!!!! W2 thermiques non nul derniere couche '
          lmax(ig)=nlay
      endif
      enddo

! pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
             lmax(ig)=1
             lmin(ig)=1
             lalim(ig)=1
         endif
      enddo 
!    
! Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
                  print*,'pb2 zw2<0'
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

!   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo

      if (iflag_thermals_ed.ge.1) then

         num(:)=0.
         denom(:)=0.
         do ig=1,ngrid
          do l=1,nlay
             num(ig)=num(ig)+zw2(ig,l)*zlev(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
             denom(ig)=denom(ig)+zw2(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
          enddo
       enddo
       do ig=1,ngrid
       if (denom(ig).gt.1.e-10) then
          zmax(ig)=2.*num(ig)/denom(ig)
          zmax0(ig)=zmax(ig)
       endif 
       enddo

       else

      do  ig=1,ngrid
! calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*  &
     &    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)  &
     &    -zlev(ig,lmax(ig)))
!pour le cas ou on prend tjs lmin=1
!       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,1))
       zmax0(ig)=zmax(ig)
      enddo


      endif
!endif iflag_thermals_ed
!
! def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1) then
! test 
              if (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
     &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))  &
     &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
     &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).gt.1e-10)  &
     &        then
!             
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
     &        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)  &
     &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
     &        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))  &
     &        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
     &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))  &
     &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
     &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
              else
              zmix(ig)=zlev(ig,lmix(ig))
              print*,'pb zmix'
              endif
          else 
              zmix(ig)=0.
          endif
!test
         if ((zmax(ig)-zmix(ig)).le.0.) then
            zmix(ig)=0.9*zmax(ig)
!            print*,'pb zmix>zmax'
         endif
      enddo
!
! calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,nlay
            if (zmix(ig).ge.zlev(ig,l).and.  &
     &          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
!
      return 
      end
