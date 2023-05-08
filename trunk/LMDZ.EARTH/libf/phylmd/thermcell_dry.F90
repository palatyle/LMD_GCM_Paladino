!
! $Id: thermcell_dry.F90 1403 2010-07-01 09:02:53Z fairhead $
!
       SUBROUTINE thermcell_dry(ngrid,nlay,zlev,pphi,ztv,alim_star,  &
     &                            lalim,lmin,zmax,wmax,lev_out)

!--------------------------------------------------------------------------
!thermcell_dry: calcul de zmax et wmax du thermique sec
! Calcul de la vitesse maximum et de la hauteur maximum pour un panache
! ascendant avec une fonction d'alimentation alim_star et sans changement 
! de phase.
! Le calcul pourrait etre sans doute simplifier.
! La temperature potentielle virtuelle dans la panache ascendant est
! la temperature potentielle virtuelle pondérée par alim_star.
!--------------------------------------------------------------------------

       IMPLICIT NONE
#include "YOMCST.h"       
#include "iniprint.h"
       INTEGER l,ig

       INTEGER ngrid,nlay
       REAL zlev(ngrid,nlay+1)
       REAL pphi(ngrid,nlay)
       REAl ztv(ngrid,nlay)
       REAL alim_star(ngrid,nlay)
       INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print

       REAL zmax(ngrid)
       REAL wmax(ngrid)

!variables locales
       REAL zw2(ngrid,nlay+1)
       REAL f_star(ngrid,nlay+1)
       REAL ztva(ngrid,nlay+1)
       REAL wmaxa(ngrid)
       REAL wa_moy(ngrid,nlay+1)
       REAL linter(ngrid),zlevinter(ngrid)
       INTEGER lmix(ngrid),lmax(ngrid),lmin(ngrid)
      CHARACTER (LEN=20) :: modname='thermcell_dry'
      CHARACTER (LEN=80) :: abort_message

!initialisations
       do ig=1,ngrid
          do l=1,nlay+1
             zw2(ig,l)=0.
             wa_moy(ig,l)=0.
          enddo
       enddo
       do ig=1,ngrid
          do l=1,nlay
             ztva(ig,l)=ztv(ig,l)
          enddo
       enddo
       do ig=1,ngrid
          wmax(ig)=0.
          wmaxa(ig)=0.
       enddo
!calcul de la vitesse a partir de la CAPE en melangeant thetav


! Calcul des F^*, integrale verticale de E^*
       f_star(:,1)=0.
       do l=1,nlay
          f_star(:,l+1)=f_star(:,l)+alim_star(:,l)
       enddo

! niveau (reel) auquel zw2 s'annule FH :n'etait pas initialise
       linter(:)=0.

! couche la plus haute concernee par le thermique. 
       lmax(:)=1

! Le niveau linter est une variable continue qui se trouve dans la couche
! lmax

       do l=1,nlay-2
         do ig=1,ngrid
            if (l.eq.lmin(ig).and.lalim(ig).gt.1) then

!------------------------------------------------------------------------
!  Calcul de la vitesse en haut de la premiere couche instable.
!  Premiere couche du panache thermique
!------------------------------------------------------------------------

               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)  &
     &                     *(zlev(ig,l+1)-zlev(ig,l))  &
     &                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))

!------------------------------------------------------------------------
! Tant que la vitesse en bas de la couche et la somme du flux de masse
! et de l'entrainement (c'est a dire le flux de masse en haut) sont
! positifs, on calcul
! 1. le flux de masse en haut  f_star(ig,l+1)
! 2. la temperature potentielle virtuelle dans la couche ztva(ig,l)
! 3. la vitesse au carr� en haut zw2(ig,l+1)
!------------------------------------------------------------------------

            else if (zw2(ig,l).ge.1e-10) then

               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+alim_star(ig,l)  &
     &                    *ztv(ig,l))/f_star(ig,l+1)
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+  &
     &                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
     &                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
! determination de zmax continu par interpolation lineaire
!------------------------------------------------------------------------

            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
!               stop'On tombe sur le cas particulier de thermcell_dry'
!               print*,'On tombe sur le cas particulier de thermcell_dry'
                zw2(ig,l+1)=0.
                linter(ig)=l+1
                lmax(ig)=l
            endif

            if (zw2(ig,l+1).lt.0.) then
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmax(ig)=l
            endif

               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))

            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
       if (prt_level.ge.1) print*,'fin calcul zw2'
!
! Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
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
      do  ig=1,ngrid
! calcul de zlevinter
          zlevinter(ig)=zlev(ig,lmax(ig)) + &
     &    (linter(ig)-lmax(ig))*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
           zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

      RETURN
      END
