!
! $Id: addfi.F 1446 2010-10-22 09:27:25Z emillour $
!
      SUBROUTINE addfi(pdt, leapf, forward,
     S          pucov, pvcov, pteta, pq   , pps ,
     S          pdufi, pdvfi, pdhfi,pdqfi, pdpfi  )

      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : planet_type
      USE comconst_mod, ONLY: kappa
      IMPLICIT NONE
c
c=======================================================================
c
c    Addition of the physical tendencies
c
c    Interface :
c    -----------
c
c      Input :
c      -------
c      pdt                    time step of integration
c      leapf                  logical
c      forward                logical
c      pucov(ip1jmp1,llm)     first component of the covariant velocity
c      pvcov(ip1ip1jm,llm)    second component of the covariant velocity
c      pteta(ip1jmp1,llm)     potential temperature
c      pts(ip1jmp1,llm)       surface temperature
c      pdufi(ip1jmp1,llm)     |
c      pdvfi(ip1jm,llm)       |   respective
c      pdhfi(ip1jmp1)         |      tendencies  (unit/s)
c      pdtsfi(ip1jmp1)        |
c
c      Output :
c      --------
c      pucov
c      pvcov
c      ph
c      pts
c
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
c    Arguments :
c    -----------
c
      REAL,INTENT(IN) :: pdt ! time step for the integration (s)
c
      REAL,INTENT(INOUT) :: pvcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(INOUT) :: pucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(INOUT) :: pteta(ip1jmp1,llm) ! potential temperature
      REAL,INTENT(INOUT) :: pq(ip1jmp1,llm,nqtot) ! tracers
      REAL,INTENT(INOUT) :: pps(ip1jmp1) ! surface pressure (Pa)
c respective tendencies (.../s) to add
      REAL,INTENT(IN) :: pdvfi(ip1jm,llm)
      REAL,INTENT(IN) :: pdufi(ip1jmp1,llm)
      REAL,INTENT(IN) :: pdqfi(ip1jmp1,llm,nqtot)
      REAL,INTENT(IN) :: pdhfi(ip1jmp1,llm)
      REAL,INTENT(IN) :: pdpfi(ip1jmp1)
c
      LOGICAL,INTENT(IN) :: leapf,forward ! not used
c
c
c    Local variables :
c    -----------------
c
      REAL xpn(iim),xps(iim),tpn,tps
      INTEGER j,k,iq,ij
      REAL,PARAMETER :: qtestw = 1.0e-15
      REAL,PARAMETER :: qtestt = 1.0e-40

      REAL SSUM
c
c-----------------------------------------------------------------------

      DO k = 1,llm
         DO j = 1,ip1jmp1
            pteta(j,k)= pteta(j,k) + pdhfi(j,k) * pdt
         ENDDO
      ENDDO

      DO  k    = 1, llm
       DO  ij   = 1, iim
         xpn(ij) = aire(   ij   ) * pteta(  ij    ,k)
         xps(ij) = aire(ij+ip1jm) * pteta(ij+ip1jm,k)
       ENDDO
       tpn      = SSUM(iim,xpn,1)/ apoln
       tps      = SSUM(iim,xps,1)/ apols

       DO ij   = 1, iip1
         pteta(   ij   ,k)  = tpn
         pteta(ij+ip1jm,k)  = tps
       ENDDO
      ENDDO
!***********************
! Correction on teta due to surface pressure changes
      DO k = 1,llm
        DO j = 1,ip1jmp1
           pteta(j,k)= pteta(j,k)*(1+pdpfi(j)*pdt/pps(j))**kappa
        ENDDO
      ENDDO
!***********************

      DO k = 1,llm
         DO j = iip2,ip1jm
            pucov(j,k)= pucov(j,k) + pdufi(j,k) * pdt
         ENDDO
      ENDDO

      DO k = 1,llm
         DO j = 1,ip1jm
            pvcov(j,k)= pvcov(j,k) + pdvfi(j,k) * pdt
         ENDDO
      ENDDO

c
      DO j = 1,ip1jmp1
         pps(j) = pps(j) + pdpfi(j) * pdt
      ENDDO
 
      if (planet_type=="earth") then
      ! earth case, special treatment for first 2 tracers (water)
       DO iq = 1, 2
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestw )
            ENDDO
         ENDDO
       ENDDO

       DO iq = 3, nqtot
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
       ENDDO
      else
      ! general case, treat all tracers equally)
       DO iq = 1, nqtot
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
       ENDDO
      endif ! of if (planet_type=="earth")

      DO  ij   = 1, iim
        xpn(ij) = aire(   ij   ) * pps(  ij     )
        xps(ij) = aire(ij+ip1jm) * pps(ij+ip1jm )
      ENDDO
      tpn      = SSUM(iim,xpn,1)/apoln
      tps      = SSUM(iim,xps,1)/apols

      DO ij   = 1, iip1
        pps (   ij     )  = tpn
        pps ( ij+ip1jm )  = tps
      ENDDO


      DO iq = 1, nqtot
        DO  k    = 1, llm
          DO  ij   = 1, iim
            xpn(ij) = aire(   ij   ) * pq(  ij    ,k,iq)
            xps(ij) = aire(ij+ip1jm) * pq(ij+ip1jm,k,iq)
          ENDDO
          tpn      = SSUM(iim,xpn,1)/apoln
          tps      = SSUM(iim,xps,1)/apols

          DO ij   = 1, iip1
            pq (   ij   ,k,iq)  = tpn
            pq (ij+ip1jm,k,iq)  = tps
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
