!
! $Header$
!
      SUBROUTINE thermcell_init(ngrid,nlay,ztv,zlay,zlev,  &
     &                  lalim,lmin,alim_star,alim_star_tot,lev_out)

!----------------------------------------------------------------------
!thermcell_init: calcul du profil d alimentation du thermique
!----------------------------------------------------------------------
      IMPLICIT NONE
#include "iniprint.h"
#include "thermcell.h"

      INTEGER l,ig
!arguments d entree
      INTEGER ngrid,nlay
      REAL ztv(ngrid,nlay)
      REAL zlay(ngrid,nlay)
      REAL zlev(ngrid,nlay+1)
!arguments de sortie
      INTEGER lalim(ngrid)
      INTEGER lmin(ngrid)
      REAL alim_star(ngrid,nlay)
      REAL alim_star_tot(ngrid)
      integer lev_out                           ! niveau pour les print
      
      REAL zzalim(ngrid)
!CR: ponderation entrainement des couches instables
!def des alim_star tels que alim=f*alim_star      


      write(lunout,*)'THERM INIT V20C '

      alim_star_tot(:)=0.
      alim_star(:,:)=0.
      lmin(:)=1
      lalim(:)=1

      do l=1,nlay-1
         do ig=1,ngrid
            if (ztv(ig,l)> ztv(ig,l+1) .and. ztv(ig,1)>=ztv(ig,l) ) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1)) 
               lalim(:)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
            endif
         enddo
      enddo
      do l=1,nlay
         do ig=1,ngrid 
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo
      alim_star_tot(:)=1.

      return
      end  
