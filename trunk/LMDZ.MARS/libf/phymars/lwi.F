      subroutine lwi (ig0,kdlon,kflev
     .                ,psi,zdblay,pdp
     .                ,newpcolc )


      use dimradmars_mod, only: ndlo2, ndlon, nflev, nir
      use yomlw_h, only: gcp, nlaylte, xi
      USE comcstfi_h, ONLY: g, cpp
      USE time_phylmdz_mod, ONLY: dtphys
      implicit none

#include "comg1d.h"
#include "callkeys.h"
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                             -   lwi    -    
C
C     PURPOSE:       Shema semi - implicite 
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c************************************************************************
c
c        0.    Declarations
c              ------------
c
c-------------------------------------------------------------------------
c        0.1   Arguments
c              ---------
c
 
      integer ig0,kdlon,kflev

      real    psi(ndlo2,kflev)
     .     ,  zdblay(ndlo2,nir,kflev)
     .     ,  pdp(ndlo2,kflev)


      real    newpcolc(ndlo2,kflev)

c-------------------------------------------------------------------------
c        0.2   local arrays
c              ------------
c
      real    di(ndlon,nflev) 
     .      , hi(ndlon,nflev) 
     .      , bi(ndlon,nflev) 

      real    ci(ndlon,nflev) 
     .      , ai(ndlon,nflev) 
      real    deltat

      real   semit, denom

      integer i, jl

c************************************************************************
c
c        1.    Initialisations
c              ---------------
c
c-----------------------------------------------------------------------
 
        deltat = dtphys * iradia
c       print*,'SEMI = ',semi, '(expl:0  semi-implicite:0.5  impl:1)'
        semit = semi * deltat
c       semi = 0.

c       print*,'dtphys,iradia,deltat,semit:',dtphys,iradia,deltat,semit
c       print*,'g,cpp',g,cpp


c************************************************************************
c
c        2.    
c              ---------------
c
c-------------------------------------------------------------------------
c        2.1   Calcul des di
c              -------------
c


      do i = 1 , nlaylte-1
        do jl = 1 , kdlon
c     -------------------
      di(jl,i) =  1 + semit * (g / pdp(jl,i) / cpp) * (
     .    ( xi(ig0+jl,1,i,nlaylte+1)
     .    + xi(ig0+jl,1,i,i+1)
     .    + xi(ig0+jl,1,i,i-1) )
     .    *    zdblay(jl,1,i)
     .  + ( xi(ig0+jl,2,i,nlaylte+1)
     .    + xi(ig0+jl,2,i,i+1)
     .    + xi(ig0+jl,2,i,i-1) )
     .    *    zdblay(jl,2,i)
     .     )
c     -------------------
        enddo
      enddo

c couche nlaylte
c ------------
c      , on enleve i,i+1 sinon on a 2 fois le cooling2space

      do jl = 1 , kdlon
c     -------------------
      di(jl,nlaylte) =  1 + semit * (g / pdp(jl,nlaylte) / cpp) * (
     .    ( xi(ig0+jl,1,nlaylte,nlaylte+1)
     .    + xi(ig0+jl,1,nlaylte,nlaylte-1) )
     .    *    zdblay(jl,1,nlaylte)
     .  + ( xi(ig0+jl,2,nlaylte,nlaylte+1)
     .    + xi(ig0+jl,2,nlaylte,nlaylte-1) )
     .    *    zdblay(jl,2,nlaylte)
     .     )
c     -------------------
      enddo

c-------------------------------------------------------------------------
c        2.2   Calcul des hi
c              -------------
c

      do i = 1 , nlaylte-1
        do jl = 1 , kdlon
c     -------------------
      hi(jl,i) =    - semit * (g / pdp(jl,i) / cpp) *
     .            (    xi(ig0+jl,1,i,i+1) * zdblay(jl,1,i+1)   
     .               + xi(ig0+jl,2,i,i+1) * zdblay(jl,2,i+1)   )
c     -------------------
        enddo
      enddo

c-------------------------------------------------------------------------
c        2.3   Calcul des bi
c              -------------
c


      do i = 2 , nlaylte
        do jl = 1 , kdlon
c     -------------------
      bi(jl,i) =   - semit * (g / pdp(jl,i) / cpp) * 
     .           (     xi(ig0+jl,1,i,i-1) * zdblay(jl,1,i-1)   
     .               + xi(ig0+jl,2,i,i-1) * zdblay(jl,2,i-1)   )
c     -------------------
        enddo
      enddo


c couche 1
c --------
c  tant qu'on a pas un calcul propre de zdblay(0) qui tienne compte de 
c    la discontinuite de temperature au sol , on met  b1 = 0


      do jl = 1 , kdlon
        bi(jl,1) = 0 
      enddo

c-------------------------------------------------------------------------
c        2.4   
c              -------------
c

c couche nlaylte
c ------------

      do jl = 1 , kdlon
c     -------------------
      ci(jl,nlaylte) = (gcp * psi(jl,nlaylte) / pdp(jl,nlaylte))
     .                   / di(jl,nlaylte)

      ai(jl,nlaylte) = - bi(jl,nlaylte) / di(jl,nlaylte)
c     -------------------
      enddo



      do i = nlaylte-1 , 1 , -1
        do jl = 1 , kdlon
c     -------------------
      denom = di(jl,i) + hi(jl,i) * ai(jl,i+1)

      ci(jl,i) = (  gcp * psi(jl,i) / pdp(jl,i)
     .             - hi(jl,i) * ci(jl,i+1)  )  / denom
 
      ai(jl,i) = -bi(jl,i) / denom
c     -------------------
        enddo
      enddo


c-------------------------------------------------------------------------
c        2.5   
c              -------------
c

c couche 1
c -------
      do jl = 1 , kdlon
        newpcolc(jl,1) = ci(jl,1)
      enddo


      do i = 2 , nlaylte
        do jl = 1 , kdlon
           newpcolc(jl,i) = ci(jl,i) + ai(jl,i) * newpcolc(jl,i-1)
        enddo
      enddo



c-------------------------------------------------------------------------
      RETURN
      END
