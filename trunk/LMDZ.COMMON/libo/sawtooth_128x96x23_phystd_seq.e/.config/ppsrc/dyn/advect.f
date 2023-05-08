










!
! $Header$
!
      SUBROUTINE advect(ucov,vcov,teta,w,massebx,masseby,du,dv,dteta)

      USE comconst_mod, ONLY: daysec
      USE logic_mod, ONLY: conser
      USE ener_mod, ONLY: gtot

      IMPLICIT NONE
c=======================================================================
c
c   Auteurs:  P. Le Van , Fr. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   *************************************************************
c   .... calcul des termes d'advection vertic.pour u,v,teta,q ...
c   *************************************************************
c        ces termes sont ajoutes a du,dv,dteta et dq .
c  Modif F.Forget 03/94 : on retire q de advect
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ip1jm,llm)
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm)
      REAL,INTENT(IN) :: teta(ip1jmp1,llm)
      REAL,INTENT(IN) :: massebx(ip1jmp1,llm)
      REAL,INTENT(IN) :: masseby(ip1jm,llm)
      REAL,INTENT(IN) :: w(ip1jmp1,llm)
      REAL,INTENT(INOUT) :: dv(ip1jm,llm)
      REAL,INTENT(INOUT) :: du(ip1jmp1,llm)
      REAL,INTENT(INOUT) :: dteta(ip1jmp1,llm)

c   Local:
c   ------

      REAL uav(ip1jmp1,llm),vav(ip1jm,llm),wsur2(ip1jmp1)
      REAL unsaire2(ip1jmp1), ge(ip1jmp1)
      REAL deuxjour, ww, gt, uu, vv

      INTEGER  ij,l

      REAL      SSUM

c-----------------------------------------------------------------------
c   2. Calculs preliminaires:
c   -------------------------

      IF (conser)  THEN
         deuxjour = 2. * daysec

         DO ij   = 1, ip1jmp1
           unsaire2(ij) = unsaire(ij) * unsaire(ij)
         ENDDO
      END IF


c------------------  -yy ----------------------------------------------
c   .  Calcul de     u

      DO  l=1,llm
         DO    ij     = iip2, ip1jmp1
            uav(ij,l) = 0.25 * ( ucov(ij,l) + ucov(ij-iip1,l) )
         ENDDO
         DO    ij     = iip2, ip1jm
            uav(ij,l) = uav(ij,l) + uav(ij+iip1,l)
         ENDDO
         DO      ij         = 1, iip1
            uav(ij      ,l) = 0.
            uav(ip1jm+ij,l) = 0.
         ENDDO
      ENDDO

c------------------  -xx ----------------------------------------------
c   .  Calcul de     v

      DO  l=1,llm
         DO    ij   = 2, ip1jm
          vav(ij,l) = 0.25 * ( vcov(ij,l) + vcov(ij-1,l) )
         ENDDO
         DO    ij   = 1,ip1jm,iip1
          vav(ij,l) = vav(ij+iim,l)
         ENDDO
         DO    ij   = 1, ip1jm-1
          vav(ij,l) = vav(ij,l) + vav(ij+1,l)
         ENDDO
         DO    ij       = 1, ip1jm, iip1
          vav(ij+iim,l) = vav(ij,l)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------

c
      DO l = 1, llmm1


c       ......   calcul de  - w/2.    au niveau  l+1   .......

        DO ij   = 1, ip1jmp1
          wsur2( ij ) = - 0.5 * w( ij,l+1 )
        ENDDO


c     .....................     calcul pour  du     ..................

        DO ij = iip2 ,ip1jm-1
          ww        = wsur2 (  ij  )     + wsur2( ij+1 ) 
          uu        = 0.5 * ( ucov(ij,l) + ucov(ij,l+1) )
          du(ij,l)  = du(ij,l)   -ww*(uu-uav(ij,l))/massebx(ij,l)
          du(ij,l+1)= du(ij,l+1) +ww*(uu-uav(ij,l+1))/massebx(ij,l+1)
        ENDDO

c     .....  correction pour  du(iip1,j,l)  ........
c     .....     du(iip1,j,l)= du(1,j,l)   .....

CDIR$ IVDEP
        DO ij   = iip1 +iip1, ip1jm, iip1
          du( ij, l  ) = du( ij -iim, l  )
          du( ij,l+1 ) = du( ij -iim,l+1 )
        ENDDO

c     .................    calcul pour   dv      .....................

        DO ij = 1, ip1jm
          ww        = wsur2( ij+iip1 )   + wsur2( ij )
          vv        = 0.5 * ( vcov(ij,l) + vcov(ij,l+1) )
          dv(ij,l)  = dv(ij, l ) - ww*(vv-vav(ij,l))/masseby(ij,l)
          dv(ij,l+1)= dv(ij,l+1) + ww*(vv-vav(ij,l+1))/masseby(ij,l+1)
        ENDDO

c

c     ............................................................
c     ...............    calcul pour   dh      ...................
c     ............................................................

c                       ---z
c       calcul de  - d( teta  * w )      qu'on ajoute a   dh
c                   ...............

        DO ij = 1, ip1jmp1
         ww            = wsur2(ij) * (teta(ij,l) + teta(ij,l+1) )
         dteta(ij, l ) = dteta(ij, l )  -  ww
         dteta(ij,l+1) = dteta(ij,l+1)  +  ww
        ENDDO

        IF( conser)  THEN
          DO ij = 1,ip1jmp1
            ge(ij)   = wsur2(ij) * wsur2(ij) * unsaire2(ij)
          ENDDO
          gt       = SSUM( ip1jmp1,ge,1 )
          gtot(l)  = deuxjour * SQRT( gt/ip1jmp1 )
        END IF

      ENDDO ! of DO l = 1, llmm1
 
      END
