!
! $Header$
!
      SUBROUTINE thermcell_closure(ngrid,nlay,r_aspect,ptimestep,rho,  &
     &   zlev,lalim,alim_star,f_star,zmax,wmax,f,lev_out)

!-------------------------------------------------------------------------
!thermcell_closure: fermeture, determination de f
!
! Modification 7 septembre 2009
! 1. On enleve alim_star_tot des arguments pour le recalculer et etre ainis
! coherent avec l'integrale au numerateur.
! 2. On ne garde qu'une version des couples wmax,zmax et wmax_sec,zmax_sec
! l'idee etant que le choix se fasse a l'appel de thermcell_closure
! 3. Vectorisation en mettant les boucles en l l'exterieur avec des if
!-------------------------------------------------------------------------
      IMPLICIT NONE

#include "iniprint.h"
#include "thermcell.h"
INTEGER ngrid,nlay
INTEGER ig,k       
REAL r_aspect,ptimestep
integer lev_out                           ! niveau pour les print

INTEGER lalim(ngrid)
REAL alim_star(ngrid,nlay)
REAL f_star(ngrid,nlay+1)
REAL rho(ngrid,nlay)
REAL zlev(ngrid,nlay)
REAL zmax(ngrid)
REAL wmax(ngrid)
REAL zdenom(ngrid)
REAL alim_star2(ngrid)
REAL f(ngrid)

REAL alim_star_tot(ngrid)
INTEGER llmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*,'THERMCELL CLOSURE 26E'

alim_star2(:)=0.
alim_star_tot(:)=0.
f(:)=0.

! Indice vertical max (max de lalim) atteint par les thermiques sur le domaine
llmax=1
do ig=1,ngrid
   if (lalim(ig)>llmax) llmax=lalim(ig)
enddo


! Calcul des integrales sur la verticale de alim_star et de
!   alim_star^2/(rho dz)
do k=1,llmax-1
   do ig=1,ngrid
      if (k<lalim(ig)) then
         alim_star2(ig)=alim_star2(ig)+alim_star(ig,k)**2  &
&                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
         alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,k)
      endif
   enddo
enddo


do ig=1,ngrid
   if (alim_star2(ig)>1.e-10) then
      f(ig)=wmax(ig)*alim_star_tot(ig)/  &
&     (max(500.,zmax(ig))*r_aspect*alim_star2(ig))
   endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TESTS POUR UNE NOUVELLE FERMETURE DANS LAQUELLE ALIM_STAR NE SERAIT
! PAS NORMALISE
!           f(ig)=f(ig)*f_star(ig,2)/(f_star(ig,lalim(ig)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end
