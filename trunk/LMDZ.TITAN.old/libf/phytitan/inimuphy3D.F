
************************************************************************
*      inimphycst()
*      initialise les constante utilisée dans la microphysique.
*
*      Oui je sais certaine constantes, voir toutes les definitions
*      faites ici devraient etre dans YOMCST.h mais la j'ai pas le 
*      courage pour le moment !
************************************************************************
       subroutine inimphycst()
       implicit none
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"
#include "YOMCST.h"

         pi = RPI

         nav = RNAVO
         rgp = R
         kbz = RKBOL

         rtit = RA
         g0   = RG

         mch4  = 16.e-3      ! kg/mol !!!!
         mc2h6 = 30.e-3
         mc2h2 = 26.e-3
         mar   = 36.4e-3
         mn2   = 28.e-3
         mair  = 28.e-3      ! l'air c'est du N2 ^^
  
         rhol      = 1.e+3     ! rho aerosols
         rhoi_ch4  = 425.
         rhoi_c2h6 = 544.6
         rhoi_c2h2 = 615.0 

         mtetach4  = 0.92
         mtetac2h6 = 0.92
         mtetac2h2 = 0.92

         return
       end
************************************************************************
*      Routine calculant/initialisant la grille verticale
*      ainsi que tous les variables communes aux routines microphysique.
*
*      CONVENTION :
*      suffixe b    : bord des couches.
*      SANS suffixe : centre des couches.
************************************************************************
       subroutine inimuphy(ihor,tplev,tplay,tzlev,tzlay,tpt)
       implicit none
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"
*
*      INPUT :
*
       integer ihor
       real tplev(nz+1),tplay(nz),tzlev(nz+1),tzlay(nz),tpt(nz)
*
*      LOCAL
*
       integer i,j
       real plev(nz+1),play(nz),zlev(nz+1),zlay(nz),pt(nz)


       if (nz.ne.llm) then
         print*,"<nz.ne.llm> dans inimuphys !"
         STOP
       endif

*-----------------------
*      INVERSION DES TABLEAUX
*-----------------------

      do j=1,nz+1 
        plev(j)=tplev(nz-j+2)
        zlev(j)=tzlev(nz-j+2)
      enddo
      do j=1,nz
        play(j)=tplay(nz-j+1)
        zlay(j)=tzlay(nz-j+1)
        pt(j) = tpt(nz-j+1)
      enddo
*
*      Calcul de la grille verticale :
*      z,dz,zb,dzb
*      
       do i=1,nz-1         
         dz(i)=zlay(i)-zlay(i+1)
         dzb(i)=zlev(i)-zlev(i+1)
       enddo
       dz(nz)=dz(nz-1)                 ! ARBITRAIRE  
       dzb(nz)=zlev(nz)-zlev(nz+1)
       z=zlay
       zb=zlev
       p=play
       pb=plev

*
*      conditions speciale pour le dernier niveau de pression
*      on crée un niveau tres tres tres fin mais qui existe quand meme !)
       pb(1) = plev(2)*1.e-7


c Interpolation de tb a partir de t
c*******************************************************
       t=pt
       do i=1,nz-1
         tb(i+1)=(t(i)+t(i+1))/2.  
       enddo
       tb(1)=t(1)
       tb(nz+1)=(t(nz)-t(nz-1))*.5+t(nz)

         return
       end
