      subroutine thermcell_dv2(ngrid,nlay,ptimestep,fm,entr,masse  &
     &    ,fraca,larga  &
     &    ,u,v,du,dv,ua,va,lev_out)
      implicit none

#include "iniprint.h"
!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
! Vectorisation, FH : 2010/03/08
!
!=======================================================================


      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)
      integer lev_out                           ! niveau pour les print

      real qa(ngrid,nlay),detr(ngrid,nlay),zf,zf2
      real wvd(ngrid,nlay+1),wud(ngrid,nlay+1)
      real gamma0(ngrid,nlay+1),gamma(ngrid,nlay+1)
      real ue(ngrid,nlay),ve(ngrid,nlay)
      LOGICAL ltherm(ngrid,nlay)
      real dua(ngrid,nlay),dva(ngrid,nlay)
      integer iter

      integer ig,k,nlarga0

!-------------------------------------------------------------------------

!   calcul du detrainement
!---------------------------

!      print*,'THERMCELL DV2 OPTIMISE 3'

      nlarga0=0.

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

!   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
         ue(ig,1)=u(ig,1)
         ve(ig,1)=v(ig,1)
      enddo

      IF(prt_level>9)WRITE(lunout,*)                                    &
     &      'WARNING on initialise gamma(1:ngrid,1)=0.'
      gamma(1:ngrid,1)=0.
      do k=2,nlay
         do ig=1,ngrid
            ltherm(ig,k)=(fm(ig,k+1)+detr(ig,k))*ptimestep > 1.e-5*masse(ig,k)
            if(ltherm(ig,k).and.larga(ig)>0.) then
               gamma0(ig,k)=masse(ig,k)  &
     &         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )  &
     &         *0.5/larga(ig)  &
     &         *1.
            else
               gamma0(ig,k)=0.
            endif
            if (ltherm(ig,k).and.larga(ig)<=0.) nlarga0=nlarga0+1
         enddo
      enddo

      gamma(:,:)=0.

      do k=2,nlay

         do ig=1,ngrid
            if (ltherm(ig,k)) then
               dua(ig,k)=ua(ig,k-1)-u(ig,k-1)
               dva(ig,k)=va(ig,k-1)-v(ig,k-1)
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               ue(ig,k)=u(ig,k)
               ve(ig,k)=v(ig,k)
            endif
         enddo


! Debut des iterations
!----------------------
do iter=1,5
         do ig=1,ngrid
! Pour memoire : calcul prenant en compte la fraction reelle
!              zf=0.5*(fraca(ig,k)+fraca(ig,k+1))
!              zf2=1./(1.-zf)
! Calcul avec fraction infiniement petite
               zf=0.
               zf2=1.

!  la première fois on multiplie le coefficient de freinage
!  par le module du vent dans la couche en dessous.
!  Mais pourquoi donc ???
               if (ltherm(ig,k)) then
!   On choisit une relaxation lineaire.
!                 gamma(ig,k)=gamma0(ig,k)
!   On choisit une relaxation quadratique.
                  gamma(ig,k)=gamma0(ig,k)*sqrt(dua(ig,k)**2+dva(ig,k)**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1)  &
     &               +(zf2*entr(ig,k)+gamma(ig,k))*u(ig,k))  &
     &               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2  &
     &                 +gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1)  &
     &               +(zf2*entr(ig,k)+gamma(ig,k))*v(ig,k))  &
     &               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2  &
     &                 +gamma(ig,k))
!                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua(ig,k),dva(ig,k)
                  dua(ig,k)=ua(ig,k)-u(ig,k)
                  dva(ig,k)=va(ig,k)-v(ig,k)
                  ue(ig,k)=(u(ig,k)-zf*ua(ig,k))*zf2
                  ve(ig,k)=(v(ig,k)-zf*va(ig,k))*zf2
               endif
         enddo
! Fin des iterations
!--------------------
enddo

      enddo ! k=2,nlay


! Calcul du flux vertical de moment dans l'environnement.
!---------------------------------------------------------
      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*ue(ig,k)
            wvd(ig,k)=fm(ig,k)*ve(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

! calcul des tendances.
!-----------------------
      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k)  &
     &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
     &               -wud(ig,k)+wud(ig,k+1))  &
     &               /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k)  &
     &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
     &               -wvd(ig,k)+wvd(ig,k+1))  &
     &               /masse(ig,k)
         enddo
      enddo


! Sorties eventuelles.
!----------------------

   if(prt_level.GE.10) then
      do k=1,nlay
         do ig=1,ngrid
           print*,'th_dv2 ig k gamma entr detr ua ue va ve wud wvd masse',ig,k,gamma(ig,k), &
     &   entr(ig,k),detr(ig,k),ua(ig,k),ue(ig,k),va(ig,k),ve(ig,k),wud(ig,k),wvd(ig,k),wud(ig,k+1),wvd(ig,k+1), &
     &   masse(ig,k)
         enddo
      enddo
   endif
!
     if (nlarga0>0) then
          print*,'WARNING !!!!!! DANS THERMCELL_DV2 '
          print*,nlarga0,' points pour lesquels laraga=0. dans un thermique'
          print*,'Il faudrait decortiquer ces points'
     endif

      return
      end
