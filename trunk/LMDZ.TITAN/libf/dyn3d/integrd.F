      SUBROUTINE integrd
     $  (  nq,vcovm1,ucovm1,tetam1,psm1,massem1,
     $     dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps,masse,phis,finvmaold )

      USE comvert_mod, ONLY: ap,bp
      USE logic_mod, ONLY: leapf
      USE serre_mod, ONLY: alphax
      USE temps_mod, ONLY: dt

      IMPLICIT NONE


c=======================================================================
cccccccccccccccccccccccccccccccccccccccccccc
c
!Mars       VERSION MARTIENNE de integrd.F
c
c   ..   modification de l'integration de  q   . 26/04/94 ..
c   ....   Si shema Van-leer pour advection de q , on n'integre pas  q
c      car q  a ete deja integre   dans "tracvl.F" appele par vanleer   ...
cccccccccccccccccccccccccccccccccccccccccccc

c
c   Auteur:  P. Le Van
c   -------
c
c   objet:
c   ------
c
c   Incrementation des tendances dynamiques
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c   Arguments:
c   ----------

      INTEGER nq

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL q(ip1jmp1,llm,nq)
      REAL ps(ip1jmp1),masse(ip1jmp1,llm),phis(ip1jmp1)

      REAL vcovm1(ip1jm,llm),ucovm1(ip1jmp1,llm)
      REAL tetam1(ip1jmp1,llm),psm1(ip1jmp1),massem1(ip1jmp1,llm)

      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dteta(ip1jmp1,llm),dp(ip1jmp1)
      REAL dq(ip1jmp1,llm,nq), finvmaold(ip1jmp1,llm)

c   Local:
c   ------

      REAL vscr( ip1jm ),uscr( ip1jmp1 ),hscr( ip1jmp1 ),pscr(ip1jmp1)
      REAL massescr( ip1jmp1,llm ), finvmasse(ip1jmp1,llm)
      REAL p(ip1jmp1,llmp1)
      REAL tpn,tps,tppn(iim),tpps(iim)
      REAL qpn,qps,qppn(iim),qpps(iim)
      REAL deltap( ip1jmp1,llm )

      INTEGER  l,ij,iq

      EXTERNAL  filtreg,massdair,pression
      EXTERNAL  SCOPY 
      REAL SSUM
      EXTERNAL SSUM

c-----------------------------------------------------------------------

      DO  l = 1,llm
        DO  ij = 1,iip1
         ucov(    ij    , l) = 0.
         ucov( ij +ip1jm, l) = 0.
         uscr(     ij      ) = 0.
         uscr( ij +ip1jm   ) = 0.
        ENDDO
      ENDDO


c    ............    integration  de       ps         ..............

      CALL SCOPY(ip1jmp1*llm, masse, 1, massescr, 1)

      DO 2 ij = 1,ip1jmp1
       pscr (ij)    = ps(ij)
       ps (ij)      = psm1(ij) + dt * dp(ij)

   2  CONTINUE
c
      DO ij = 1,ip1jmp1
        IF( ps(ij).LT.0. ) THEN
         PRINT*,' Au point ij = ',ij, ' , pression sol neg. ', ps(ij)
         write(*,*)'psm1(ij)=',psm1(ij),' dp(ij)=',dp(ij),
     &             'dp(ij)*dt=',dp(ij)*dt
         STOP' dans integrd'
        ENDIF
      ENDDO
c
      IF( alphax.NE.0. )   THEN
         DO  ij    = 1, iim
          tppn(ij) = aire(   ij   ) * ps(  ij    )
          tpps(ij) = aire(ij+ip1jm) * ps(ij+ip1jm)
         ENDDO
          tpn      = SSUM(iim,tppn,1)/apoln
          tps      = SSUM(iim,tpps,1)/apols
         DO ij   = 1, iip1
          ps(   ij   )  = tpn
          ps(ij+ip1jm)  = tps
         ENDDO
      ENDIF
c
c  ... Calcul  de la nouvelle masse d'air au dernier temps integre t+1 ...
c
      CALL pression ( ip1jmp1, ap, bp, ps, p )
      CALL massdair (     p  , masse         )

      CALL SCOPY( ijp1llm  , masse, 1, finvmasse,  1      )
      CALL filtreg( finvmasse, jjp1, llm, -2, 2, .TRUE., 1  )
c

c    ............   integration  de  ucov, vcov,  h     ..............

      DO 10 l = 1,llm

      DO 4 ij = iip2,ip1jm
      uscr( ij )   =  ucov( ij,l )
      ucov( ij,l ) = ucovm1( ij,l ) + dt * du( ij,l )
   4  CONTINUE

      DO 5 ij = 1,ip1jm
      vscr( ij )   =  vcov( ij,l )
      vcov( ij,l ) = vcovm1( ij,l ) + dt * dv( ij,l )
   5  CONTINUE

      DO 6 ij = 1,ip1jmp1
      hscr( ij )    =  teta(ij,l)
      teta ( ij,l ) = tetam1(ij,l) *  massem1(ij,l) / masse(ij,l) 
     $                + dt * dteta(ij,l) / masse(ij,l)
   6  CONTINUE

c   ....  Calcul de la valeur moyenne, unique  aux poles pour  teta    ......
c
c
      DO  ij   = 1, iim
        tppn(ij) = aire(   ij   ) * teta(  ij    ,l)
        tpps(ij) = aire(ij+ip1jm) * teta(ij+ip1jm,l)
      ENDDO
        tpn      = SSUM(iim,tppn,1)/apoln
        tps      = SSUM(iim,tpps,1)/apols

      DO ij   = 1, iip1
        teta(   ij   ,l)  = tpn
        teta(ij+ip1jm,l)  = tps
      ENDDO
c

      IF(leapf)  THEN
         CALL SCOPY ( ip1jmp1, uscr(1), 1, ucovm1(1, l), 1 )
         CALL SCOPY (   ip1jm, vscr(1), 1, vcovm1(1, l), 1 )
         CALL SCOPY ( ip1jmp1, hscr(1), 1, tetam1(1, l), 1 )
      END IF

  10  CONTINUE


c
c   .......  integration de   q   ......
c
c
c     .....   FIN  de l'integration  de   q    .......

c    .................................................................


      IF( leapf )  THEN
         CALL SCOPY (    ip1jmp1 ,  pscr   , 1,   psm1  , 1 )
         CALL SCOPY ( ip1jmp1*llm, massescr, 1,  massem1, 1 )
      END IF

      RETURN
      END
