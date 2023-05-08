      subroutine lwxd (ig0,kdlon,kflev,emis
     .                ,aer_t,co2_u,co2_up)

c----------------------------------------------------------------------
c     LWXD   computes transmission function and exchange coefficiants
c                        for distant layers
c                          (co2 / aerosols)
c                       (bands 1 and 2 of co2) 
c----------------------------------------------------------------------
c  
c              |---|---|---|---|---|---|---|---|
c   kflev+1    |   |   |   |   |   |   |   | 0 |  (space)
c              |---|---|---|---|---|---|---|---|
c    kflev     |   |***|***|***|***|   | 0 |   |
c              |---|---|---|---|---|---|---|---|
c     ...      |   |***|***|***|   | 0 |   |   |
c              |---|---|---|---|---|---|---|---|
c      4       |   |***|***|   | 0 |   |***|   | 
c              |---|---|---|---|---|---|---|---|
c      3       |   |***|   | 0 |   |***|***|   |
c              |---|---|---|---|---|---|---|---|
c      2       |   |   | 0 |   |   |***|***|   |
c              |---|---|---|---|---|---|---|---|
c      1       |   | 0 |   |   |***|***|***|   |
c              |---|---|---|---|---|---|---|---|
c      0       | 0 |   |   |***|***|***|***|   |  (ground)
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
 
#include "yomlw.h"
#include "callkeys.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer ig0
      integer kdlon      ! part of ngrid
      integer kflev      ! part of nalyer
 
      real emis (ndlo2)                  ! surface emissivity
      real aer_t (ndlo2,nuco2,kflev+1)   ! transmission (aer)
      real co2_u (ndlo2,nuco2,kflev+1)   ! absorber amounts (co2)
      real co2_up (ndlo2,nuco2,kflev+1)  ! idem scaled by the pressure (co2)

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer ja,jl,jk,jkk,ndim
      parameter(ndim = ndlon*nuco2*(nflev+2)*(nflev+2))


      real zu (ndlon,nuco2)
      real zup (ndlon,nuco2)
      real zt_co2 (ndlon,nuco2)
      real zt_aer (ndlon,nuco2)

      real ksi (ndlon,nuco2,0:nflev+1,0:nflev+1)
      real ksi_emis (ndlon,nuco2,0:nflev+1,0:nflev+1)
      real trans (ndlon,nuco2,0:nflev+1,0:nflev+1)
      real trans_emis (ndlon,nuco2,0:nflev+1,0:nflev+1)

c----------------------------------------------------------------------
      call zerophys(ndim,ksi_emis)
c----------------------------------------------------------------------
c         1.0   Transmission functions
c               ----------------------

c----------------------------------------------------------------------
c        1.1   Direct transmission 
c              -------------------

      do jk = 1 , nlaylte+1
        do jkk = jk , nlaylte+1

          do ja = 1 , nuco2
            do jl = 1 , kdlon
c                                                                   co2
c                                                                   ---
              zu(jl,ja) =  co2_u(jl,ja,jk)  - co2_u(jl,ja,jkk)
              zup(jl,ja) = co2_up(jl,ja,jk) - co2_up(jl,ja,jkk)
c                                                                   aer
c                                                                   ---
              zt_aer(jl,ja)= aer_t(jl,ja,jk)
     .                      /aer_t(jl,ja,jkk)

            enddo
          enddo

          call lwtt(kdlon,zu,zup,nuco2,zt_co2)
c                                                            co2 and aer
c                                                            -----------
          do ja = 1 , nuco2
            do jl = 1 , kdlon
              trans(jl,ja,jk,jkk) = zt_co2(jl,ja) * zt_aer(jl,ja)
            enddo
          enddo
c                                                       trans reciprocity
c                                                       -----------------
          do ja = 1 , nuco2
            do jl = 1 , kdlon
              trans(jl,ja,jkk,jk) = trans(jl,ja,jk,jkk) 
c           if (trans(jl,ja,jk,jkk) .LT. 0 ) then
c           print*,'trans bande',ja,jk,jkk,trans(jl,ja,jk,jkk)
c           endif
c           if (trans(jl,ja,jk,jkk) .GT. 1) then
c           print*,'trans bande',ja,jk,jkk,trans(jl,ja,jk,jkk)
c           trans(jl,ja,jk,jkk)=1
c           print*,'trans bande',ja,jk,jkk,trans(jl,ja,jk,jkk)
c           endif

            enddo
          enddo
     
        enddo
      enddo

c----------------------------------------------------------------------
c         1.2   Transmission with reflexion
c               ---------------------------

      do jk = 1 , nlaylte+1
        do jkk = jk , nlaylte+1

      if (callemis) then
          do ja = 1 , nuco2
            do jl = 1 , kdlon
c                                                                   co2
c                                                                   ---
              zu(jl,ja)  = 2 * co2_u(jl,ja,1)  - co2_u(jl,ja,jk)
     .                                         - co2_u(jl,ja,jkk)
              zup(jl,ja) = 2 * co2_up(jl,ja,1) - co2_up(jl,ja,jk)
     .                                         - co2_up(jl,ja,jkk)
c                                                                   aer
c                                                                   ---
                zt_aer(jl,ja)  = aer_t(jl,ja,1)
     .                         * aer_t(jl,ja,1)
     .                         / aer_t(jl,ja,jk)
     .                         / aer_t(jl,ja,jkk)
            enddo
          enddo

          call lwtt(kdlon,zu,zup,nuco2,zt_co2)
c                                                            co2 and aer
c                                                            -----------
          do ja = 1 , nuco2
            do jl = 1 , kdlon
              trans_emis(jl,ja,jk,jkk) = zt_co2(jl,ja)
     .                                 * zt_aer(jl,ja)
            enddo
          enddo

      else

          do ja = 1 , nuco2
            do jl = 1 , kdlon
              trans_emis(jl,ja,jk,jkk) = 1.
            enddo
          enddo

      endif
c                                                       trans reciprocity
c                                                       -----------------
          do ja = 1 , nuco2
            do jl = 1 , kdlon
              trans_emis(jl,ja,jkk,jk) = trans_emis(jl,ja,jk,jkk)
c         if (trans_emis(jl,ja,jk,jkk) .LT. 0 
c    .                .OR.  trans_emis(jl,ja,jk,jkk) .GT. 1) then
c     print*,'trans_emis bande',ja,jk,jkk,trans_emis(jl,ja,jk,jkk)
c         endif
            enddo
          enddo

        enddo
      enddo

c----------------------------------------------------------------------
c         2.0   Exchange Coefficiants
c               ---------------------

      do jk = 1 , nlaylte-2
        do jkk = jk+2 , nlaylte
          do ja = 1 , nuco2
            do jl = 1 , kdlon

      ksi(jl,ja,jk,jkk) = 
     .            trans(jl,ja,jk+1,jkk)   - trans(jl,ja,jk,jkk)
     .          - trans(jl,ja,jk+1,jkk+1) + trans(jl,ja,jk,jkk+1)

      ksi_emis(jl,ja,jk,jkk) = 
     .   trans_emis(jl,ja,jk,jkk)   - trans_emis(jl,ja,jk+1,jkk)
     . - trans_emis(jl,ja,jk,jkk+1) + trans_emis(jl,ja,jk+1,jkk+1)

c       if (ksi(jl,ja,jk,jkk) .LT. 0 ) then
c           print*,'ksi bande',ja,jk,jkk,ksi(jl,ja,jk,jkk)
c           ksi(jl,ja,jk,jkk)=0
c           print*,'ksi bande',ja,jk,jkk,ksi(jl,ja,jk,jkk)
c       endif
c       if (ksi(jl,ja,jk,jkk) .GT. 1) then
c           print*,'ksi bande',ja,jk,jkk,ksi(jl,ja,jk,jkk)
c           ksi(jl,ja,jk,jkk)=1
c           print*,'ksi bande',ja,jk,jkk,ksi(jl,ja,jk,jkk)
c       endif

c       if (ksi_emis(jl,ja,jk,jkk) .LT. 0 
c    .                .OR.  ksi_emis(jl,ja,jk,jkk) .GT. 1) then
c     print*,'ksi_emis bande',ja,jk,jkk,ksi_emis(jl,ja,jk,jkk)
c       endif

      xi(ig0+jl,ja,jk,jkk) = ksi(jl,ja,jk,jkk)
     .      + ksi_emis(jl,ja,jk,jkk) * (1 - emis(jl))

c                                                        ksi reciprocity
c                                                        ---------------
      ksi(jl,ja,jkk,jk)      = ksi(jl,ja,jk,jkk)
      ksi_emis(jl,ja,jkk,jk) = ksi_emis(jl,ja,jk,jkk)
      xi(ig0+jl,ja,jkk,jk)   = xi(ig0+jl,ja,jk,jkk)

            enddo
          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         2.1   Save xi_emis for neighbours (lwxn.F)
c               -----------------------------------

      do jk = 1 , nlaylte-1
        do ja = 1 , nuco2
          do jl = 1 , kdlon

c     ksi_emis(jl,ja,jk,jk+1) = 
c    .   trans_emis(jl,ja,jk,jk+1)   - trans_emis(jl,ja,jk+1,jk+1)
c    . - trans_emis(jl,ja,jk,jk+2) + trans_emis(jl,ja,jk+1,jk+2)

      xi_emis(ig0+jl,ja,jk) =
     .                 ksi_emis(jl,ja,jk,jk+1) * (1-emis(jl))

          enddo
        enddo
      enddo

c----------------------------------------------------------------------
      return
      end
