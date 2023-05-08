      subroutine thermcell_dtke(ngrid,nlay,nsrf,ptimestep,fm0,entr0,  &
     &           rg,pplev,tke)
      implicit none

#include "iniprint.h"
!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
!=======================================================================

      integer ngrid,nlay,nsrf

      real ptimestep
      real masse0(ngrid,nlay),fm0(ngrid,nlay+1),pplev(ngrid,nlay+1)
      real entr0(ngrid,nlay),rg
      real tke(ngrid,nlay,nsrf)
      real detr0(ngrid,nlay)


      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay)
      real q(ngrid,nlay)
      integer lev_out                           ! niveau pour les print

      real qa(ngrid,nlay),detr(ngrid,nlay),wqd(ngrid,nlay+1)

      real zzm

      integer ig,k
      integer isrf


      lev_out=0


      if (prt_level.ge.1) print*,'Q2 THERMCEL_DQ 0'

!   calcul du detrainement
      do k=1,nlay
         detr0(:,k)=fm0(:,k)-fm0(:,k+1)+entr0(:,k)
         masse0(:,k)=(pplev(:,k)-pplev(:,k+1))/RG
      enddo


! Decalage vertical des entrainements et detrainements.
      masse(:,1)=0.5*masse0(:,1)
      entr(:,1)=0.5*entr0(:,1)
      detr(:,1)=0.5*detr0(:,1)
      fm(:,1)=0.
      do k=1,nlay-1
         masse(:,k+1)=0.5*(masse0(:,k)+masse0(:,k+1))
         entr(:,k+1)=0.5*(entr0(:,k)+entr0(:,k+1))
         detr(:,k+1)=0.5*(detr0(:,k)+detr0(:,k+1))
         fm(:,k+1)=fm(:,k)+entr(:,k)-detr(:,k)
      enddo
      fm(:,nlay+1)=0.

!   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo



do isrf=1,nsrf

   q(:,:)=tke(:,:,isrf)

    if (1==1) then
      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.  &
     &         1.e-5*masse(ig,k)) then
         qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
     &         /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
            if (qa(ig,k).lt.0.) then
!               print*,'qa<0!!!'
            endif
            if (q(ig,k).lt.0.) then
!               print*,'q<0!!!'
            endif
         enddo
      enddo

! Calcul du flux subsident

      do k=2,nlay
         do ig=1,ngrid
            wqd(ig,k)=fm(ig,k)*q(ig,k)
            if (wqd(ig,k).lt.0.) then
!               print*,'wqd<0!!!'
            endif
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo
     

! Calcul des tendances
      do k=1,nlay
         do ig=1,ngrid
            q(ig,k)=q(ig,k)+(detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)  &
     &               -wqd(ig,k)+wqd(ig,k+1))  &
     &               *ptimestep/masse(ig,k)
         enddo
      enddo

 endif

   tke(:,:,isrf)=q(:,:)

enddo

      return
      end
