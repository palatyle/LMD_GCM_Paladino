!
! $Id: $
!
      SUBROUTINE caldyn
     $ (itau,ucov,vcov,teta,ps,masse,pk,pkf,tsurpk,phis ,
     $  phi,conser,du,dv,dteta,dp,w,pbaru,pbarv,time )

      USE comvert_mod, ONLY: ap,bp

      IMPLICIT NONE

!=======================================================================
!
!  Auteur :  P. Le Van
!
!   Objet:
!   ------
!
!   Calcul des tendances dynamiques.
!
! Modif 04/93 F.Forget
!=======================================================================

!-----------------------------------------------------------------------
!   0. Declarations:
!   ----------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

!   Arguments:
!   ----------

      LOGICAL,INTENT(IN) :: conser ! triggers printing some diagnostics
      INTEGER,INTENT(IN) :: itau ! time step index
      REAL,INTENT(IN) :: vcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ip1jmp1,llm) ! potential temperature
      REAL,INTENT(IN) :: ps(ip1jmp1) ! surface pressure
      REAL,INTENT(IN) :: phis(ip1jmp1) ! geopotential at the surface
      REAL,INTENT(IN) :: pk(ip1jmp1,llm) ! Exner at mid-layer
      REAL,INTENT(IN) :: pkf(ip1jmp1,llm) ! filtered Exner
      REAL,INTENT(IN) :: tsurpk(ip1jmp1,llm) ! cpp * temperature / pk
      REAL,INTENT(IN) :: phi(ip1jmp1,llm) ! geopotential
      REAL,INTENT(OUT) :: masse(ip1jmp1,llm) ! air mass
      REAL,INTENT(OUT) :: dv(ip1jm,llm) ! tendency on vcov
      REAL,INTENT(OUT) :: du(ip1jmp1,llm) ! tendency on ucov
      REAL,INTENT(OUT) :: dteta(ip1jmp1,llm) ! tenddency on teta
      REAL,INTENT(OUT) :: dp(ip1jmp1) ! tendency on ps
      REAL,INTENT(OUT) :: w(ip1jmp1,llm) ! vertical velocity
      REAL,INTENT(OUT) :: pbaru(ip1jmp1,llm) ! mass flux in the zonal direction
      REAL,INTENT(OUT) :: pbarv(ip1jm,llm) ! mass flux in the meridional direction
      REAL,INTENT(IN) :: time ! current time

!   Local:
!   ------

      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL ang(ip1jmp1,llm),p(ip1jmp1,llmp1)
      REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm),psexbarxy(ip1jm)
      REAL vorpot(ip1jm,llm)
      REAL ecin(ip1jmp1,llm),convm(ip1jmp1,llm)
      REAL bern(ip1jmp1,llm)
      REAL massebxy(ip1jm,llm)
      REAL temp(ip1jmp1,llm)

      INTEGER   ij,l

!-----------------------------------------------------------------------
!   Compute dynamical tendencies:
!--------------------------------

      ! compute contravariant winds ucont() and vcont
      CALL covcont  ( llm    , ucov    , vcov , ucont, vcont        )
      ! compute pressure p()
      CALL pression ( ip1jmp1, ap      , bp   ,  ps  , p            )
      ! compute psexbarxy() XY-area weighted-averaged surface pressure (what for?)
      CALL psextbar (   ps   , psexbarxy                            )
      ! compute mass in each atmospheric mesh: masse()
      CALL massdair (    p   , masse                                )
      ! compute X and Y-averages of mass, massebx() and masseby()
      CALL massbar  (   masse, massebx , masseby                    )
      ! compute XY-average of mass, massebxy()
      call massbarxy(   masse, massebxy                             )
      ! compute mass fluxes pbaru() and pbarv()
      CALL flumass  ( massebx, masseby , vcont, ucont ,pbaru, pbarv )
      ! compute dteta() , horizontal converging flux of theta
      CALL dteta1   (   teta , pbaru   , pbarv, dteta               )
      ! compute convm(), horizontal converging flux of mass
      CALL convmas  (   pbaru, pbarv   , convm                      )

      ! compute pressure variation due to mass convergence
      DO ij =1, ip1jmp1
         dp( ij ) = convm( ij,1 ) / airesurg( ij )
      ENDDO

      ! compute vertical velocity w()
      CALL vitvert ( convm  , w                                  )
      ! compute potential vorticity vorpot()
      CALL tourpot ( vcov   , ucov  , massebxy  , vorpot         )
      ! compute rotation induced du() and dv()
      CALL dudv1   ( vorpot , pbaru , pbarv     , du     , dv    )
      ! compute kinetic energy ecin()
      CALL enercin ( vcov   , ucov  , vcont     , ucont  , ecin  )
      ! compute Bernouilli function bern()
      CALL bernoui ( ip1jmp1, llm   , phi       , ecin   , bern  )
      ! compute and add du() and dv() contributions from Bernouilli and pressure
      CALL dudv2   ( tsurpk , pkf   , bern      , du     , dv    )


      DO l=1,llm
         DO ij=1,ip1jmp1
            ang(ij,l) = ucov(ij,l) + constang(ij)
      ENDDO
      ENDDO

      ! compute vertical advection contributions to du(), dv() and dteta()
      CALL advect( ang, vcov, teta, w, massebx, masseby, du, dv,dteta ) 

!  WARNING probleme de peridocite de dv sur les PC/linux. Pb d'arrondi 
!          probablement. Observe sur le code compile avec pgf90 3.0-1 

      DO l = 1, llm
         DO ij = 1, ip1jm, iip1
           IF( dv(ij,l).NE.dv(ij+iim,l) )  THEN
!         PRINT *,'!!!ATTENTION!!! probleme de periodicite sur vcov',  
!    ,   ' dans caldyn'
!         PRINT *,' l,  ij = ', l, ij, ij+iim,dv(ij+iim,l),dv(ij,l)
          dv(ij+iim,l) = dv(ij,l)
          endif
         enddo
      enddo

!-----------------------------------------------------------------------
!   Output some control variables:
!---------------------------------

      IF( conser )  THEN
        CALL sortvarc
     $ (itau,ucov,tsurpk,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov)

      ENDIF

      END
