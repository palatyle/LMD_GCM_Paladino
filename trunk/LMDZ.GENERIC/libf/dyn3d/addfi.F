      SUBROUTINE addfi(nq, pdt, leapf, forward,
     S          pucov, pvcov, pteta, pq   , pps , pmasse ,
     S          pdufi, pdvfi, pdhfi,pdqfi, pdpfi  )

      USE comvert_mod, ONLY: ap,bp     

      IMPLICIT NONE
c
c=======================================================================
c
c    Addition of the physical tendencies
c    Modif special Mars :
c       - no special treatment of iq = 1, 2, FF 2003
c       - recompute mass after incrementation of ps (EM,FF 2008)
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
c      pdhfi(ip1jmp1)         |      tendencies
c      pdtsfi(ip1jmp1)        |
c
c      Output :
c      --------
c      pucov
c      pvcov
c      ph
c      pts
c      pmasse(ip1jmp1,llm)    ! mass
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
      INTEGER nq

      REAL pdt
c
      REAL pvcov(ip1jm,llm),pucov(ip1jmp1,llm)
      REAL pteta(ip1jmp1,llm),pq(ip1jmp1,llm,nq),pps(ip1jmp1)
c
      REAL pdvfi(ip1jm,llm),pdufi(ip1jmp1,llm)
      REAL pdqfi(ip1jmp1,llm,nq),pdhfi(ip1jmp1,llm),pdpfi(ip1jmp1)

      real pmasse(ip1jmp1,llm)
c
      LOGICAL leapf,forward
c
c
c    Local variables :
c    -----------------
c
      REAL xpn(iim),xps(iim),tpn,tps
      INTEGER j,k,iq,ij
      REAL qtestw, qtestt
      PARAMETER ( qtestw = 1.0e-15 )
      PARAMETER ( qtestt = 1.0e-30 )

      real p(ip1jmp1,llmp1) ! pressure

      REAL SSUM
      EXTERNAL SSUM
c
c-----------------------------------------------------------------------
! increment potential temperature
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
c

! increment covariant zonal wind
      DO k = 1,llm
         DO j = iip2,ip1jm
            pucov(j,k)= pucov(j,k) + pdufi(j,k) * pdt
         ENDDO
      ENDDO

! increment covariant meridional wind
      DO k = 1,llm
         DO j = 1,ip1jm
            pvcov(j,k)= pvcov(j,k) + pdvfi(j,k) * pdt
         ENDDO
      ENDDO

c
! increment surface pressure
      DO j = 1,ip1jmp1
         pps(j) = pps(j) + pdpfi(j) * pdt
      ENDDO

c     DO iq = 1, 2 ! special Mars: no special treatment for water
c        DO k = 1,llm
c           DO j = 1,ip1jmp1
c              pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
c              pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestw )
c           ENDDO
c        ENDDO
c     ENDDO

! increment tracers
      DO iq = 1, nq
         DO k = 1,llm
            DO j = 1,ip1jmp1
!               pq(j,k,iq)=  pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt ) ! forbid negative tracer values
            ENDDO
         ENDDO
      ENDDO
!      print*,'MAJOR TRACER MOD in addfi for test'
                                       

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

! recompute mass: (to be synchronous with update of ps)
      CALL pression ( ip1jmp1, ap, bp, pps, p )
      CALL massdair (     p  , pmasse         )
 


       DO iq = 1, nq
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
