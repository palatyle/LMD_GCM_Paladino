      subroutine lwxb (ig0,kdlon,kflev
     .                ,emis
     .                ,aer_t,co2_u,co2_up)

c----------------------------------------------------------------------
c     LWXB   computes transmission function and exchange coefficients
c                          for boundaries 
c                          (co2 / aerosols)
c                       (bands 1 and 2 of co2) 
c----------------------------------------------------------------------
c 
c              |---|---|---|---|---|---|---|---|
c   kflev+1    |***|***|***|***|***|***|***| 0 |  (space)
c              |---|---|---|---|---|---|---|---|
c    kflev     |***|   |   |   |   |   | 0 |***|
c              |---|---|---|---|---|---|---|---|
c     ...      |***|   |   |   |   | 0 |   |***|
c              |---|---|---|---|---|---|---|---|
c      4       |***|   |   |   | 0 |   |   |***|
c              |---|---|---|---|---|---|---|---|
c      3       |***|   |   | 0 |   |   |   |***|
c              |---|---|---|---|---|---|---|---|
c      2       |***|   | 0 |   |   |   |   |***|
c              |---|---|---|---|---|---|---|---|
c      1       |***| 0 |   |   |   |   |   |***|
c              |---|---|---|---|---|---|---|---|
c      0       | 0 |***|***|***|***|***|***|***|  (ground)
c              |---|---|---|---|---|---|---|---|
c                0   1   2   3   4  ...  k |k+1
c             (ground)                    (space)
c
c  (*)  xi computed in this subroutine
c----------------------------------------------------------------------

      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "callkeys.h"
 
#include "yomlw.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer kdlon          ! part of ngrid
      integer kflev          ! part of nalyer

      real emis (ndlo2)                  ! surface emissivity
      real aer_t (ndlo2,nuco2,kflev+1)   ! transmission (aer)
      real co2_u (ndlo2,nuco2,kflev+1)   ! absorber amounts (co2)
      real co2_up (ndlo2,nuco2,kflev+1)  ! idem scaled by the pressure (co2)

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer ja,jl,jk,ig0

      real zt_co2 (ndlon,nuco2)
      real zt_aer (ndlon,nuco2)
      real zu (ndlon,nuco2)
      real zup (ndlon,nuco2)
c                                          2 for ground(1) and space(2)
      real trans (ndlon,nuco2,2,0:nflev+1)
      real ksi (ndlon,nuco2,2,0:nflev+1)
c                                                       only for space
      real trans_emis (ndlon,nuco2,0:nflev+1)
      real ksi_emis (ndlon,nuco2,0:nflev+1)

c*************************************************************************
c         1.0   Transmissions 
c               -------------
c----------------------------------------------------------------------
c         1.1   Direct Transmission 
c               -------------------

c                                                                 space
c                                                                 -----
      do jk = 1 , nlaylte+1

        do ja = 1 , nuco2
          do jl = 1 , kdlon
            zu(jl,ja)  = co2_u(jl,ja,jk)
            zup(jl,ja) = co2_up(jl,ja,jk)
            zt_aer(jl,ja) = aer_t(jl,ja,jk)
          enddo
        enddo

        call lwtt(kdlon,zu,zup,nuco2,zt_co2)

        do ja = 1 , nuco2
          do jl = 1 , kdlon
            trans(jl,ja,2,jk)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo

      enddo
c                                                                 ground
c                                                                 -----
      do jk = 1 , nlaylte+1

        do ja = 1 , nuco2
          do jl = 1 , kdlon
            zu(jl,ja) =  co2_u(jl,ja,1)  - co2_u(jl,ja,jk)
            zup(jl,ja) = co2_up(jl,ja,1) - co2_up(jl,ja,jk)
            zt_aer(jl,ja) = aer_t(jl,ja,1) /aer_t(jl,ja,jk)
          enddo
        enddo

        call lwtt(kdlon,zu,zup,nuco2,zt_co2)

        do ja = 1 , nuco2
          do jl = 1 , kdlon
            trans(jl,ja,1,jk)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo

      enddo

c----------------------------------------------------------------------
c         1.2   Transmission with reflexion
c               ---------------------------

c                                                                 space
c                                                                 -----
      do jk = 1 , nlaylte+1

        do ja = 1 , nuco2
          do jl = 1 , kdlon

            zu(jl,ja) =   2 * co2_u(jl,ja,1)  - co2_u(jl,ja,jk)
            zup(jl,ja) =  2 * co2_up(jl,ja,1) - co2_up(jl,ja,jk)
            zt_aer(jl,ja) = aer_t(jl,ja,1) 
     .                    * aer_t(jl,ja,1)
     .                    / aer_t(jl,ja,jk)

          enddo
        enddo

        call lwtt(kdlon,zu,zup,nuco2,zt_co2)

        do ja = 1 , nuco2
          do jl = 1 , kdlon
            trans_emis(jl,ja,jk)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo

      enddo

c*************************************************************************
c         2.0   Exchange Coefficiants
c               ---------------------

      do jk = 1 , nlaylte
        do ja = 1 , nuco2
          do jl = 1 , kdlon

c-------------------------------------------------------------------------
c        2.1    colling to space  (from layer 1,nlaylte toward "layer" nlaylte+1)
c               ----------------


      ksi(jl,ja,2,jk) = trans(jl,ja,2,jk+1)
     .                - trans(jl,ja,2,jk)

      ksi_emis(jl,ja,jk) = trans_emis(jl,ja,jk)
     .                   - trans_emis(jl,ja,jk+1) 

      xi(ig0+jl,ja,jk,nlaylte+1)= ksi(jl,ja,2,jk) 
     .                        + ksi_emis(jl,ja,jk)* (1 - emis(jl))

c                                                         ksi Reciprocity
c                                                         ---------------
      xi(ig0+jl,ja,nlaylte+1,jk)      = xi(ig0+jl,ja,jk,nlaylte+1)

c-------------------------------------------------------------------------
c        2.2    echange with ground  (from "layer" 0 toward layers 1,nlaylte)
c               -------------------


      ksi(jl,ja,1,jk) = trans(jl,ja,1,jk) 
     .                - trans(jl,ja,1,jk+1)

      xi(ig0+jl,ja,0,jk) = ksi(jl,ja,1,jk) * emis(jl)

c                                                         ksi Reciprocity
c                                                         ---------------
      xi(ig0+jl,ja,jk,0) = xi(ig0+jl,ja,0,jk)

c-------------------------------------------------------------------------
          enddo
        enddo
      enddo

c-------------------------------------------------------------------------
c       2.3     echange ground-space  (from "layer" 0 toward "layer" nlaylte+1)
c               ----------------------

c Is not used because we use sigma T4 for the ground budget in physiq.F

      do ja = 1 , nuco2
        do jl = 1 , kdlon

      ksi(jl,ja,1,nlaylte+1) = trans(jl,ja,1,nlaylte+1) 
      xi(ig0+jl,ja,0,nlaylte+1) = ksi(jl,ja,1,nlaylte+1) * emis(jl)

c                                                         ksi Reciprocity
c                                                         ---------------
      xi(ig0+jl,ja,nlaylte+1,0) = xi(ig0+jl,ja,0,nlaylte+1)

        enddo
      enddo

c-------------------------------------------------------------------------
      return
      end
