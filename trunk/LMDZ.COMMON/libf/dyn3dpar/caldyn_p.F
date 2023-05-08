!
! $Id: $
!
#undef DEBUG_IO
c#define DEBUG_IO

      SUBROUTINE caldyn_p
     $ (itau,ucov,vcov,teta,ps,masse,pk,pkf,tsurpk,phis ,
     $  phi,conser,du,dv,dteta,dp,w,pbaru,pbarv,time )
      USE parallel_lmdz
      USE Write_Field_p
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

      REAL,SAVE :: vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL,SAVE :: ang(ip1jmp1,llm)
      REAL,SAVE :: p(ip1jmp1,llmp1)
      REAL,SAVE :: massebx(ip1jmp1,llm),masseby(ip1jm,llm)
      REAL,SAVE :: psexbarxy(ip1jm)
      REAL,SAVE :: vorpot(ip1jm,llm)
      REAL,SAVE :: ecin(ip1jmp1,llm)
      REAL,SAVE :: bern(ip1jmp1,llm)
      REAL,SAVE :: massebxy(ip1jm,llm)
      REAL,SAVE :: convm(ip1jmp1,llm)
!      REAL,SAVE :: temp(ip1jmp1,llm)
      INTEGER   ij,l,ijb,ije,ierr

!-----------------------------------------------------------------------
!   Compute dynamical tendencies:
!--------------------------------

      ! compute contravariant winds ucont() and vcont
      CALL covcont_p  ( llm    , ucov    , vcov , ucont, vcont        )
      ! compute pressure p()
      CALL pression_p ( ip1jmp1, ap      , bp   ,  ps  , p            )
cym      CALL psextbar (   ps   , psexbarxy                          )
c$OMP BARRIER
      ! compute mass in each atmospheric mesh: masse()
      CALL massdair_p (    p   , masse                                )
      ! compute X and Y-averages of mass, massebx() and masseby()
      CALL massbar_p  (   masse, massebx , masseby                    )
      ! compute XY-average of mass, massebxy()
      call massbarxy_p(   masse, massebxy                             )
      ! compute mass fluxes pbaru() and pbarv()
      CALL flumass_p  ( massebx, masseby , vcont, ucont ,pbaru, pbarv )
      ! compute dteta() , horizontal converging flux of theta
      CALL dteta1_p   (   teta , pbaru   , pbarv, dteta               )
      ! compute convm(), horizontal converging flux of mass
      CALL convmas1_p  (   pbaru, pbarv   , convm                      )
c$OMP BARRIER      
      CALL convmas2_p  (   convm                      )
c$OMP BARRIER
#ifdef DEBUG_IO
c$OMP BARRIER
c$OMP MASTER
      call WriteField_p('ucont',reshape(ucont,(/iip1,jmp1,llm/)))
      call WriteField_p('vcont',reshape(vcont,(/iip1,jjm,llm/)))
      call WriteField_p('p',reshape(p,(/iip1,jmp1,llmp1/)))
      call WriteField_p('masse',reshape(masse,(/iip1,jmp1,llm/)))
      call WriteField_p('massebx',reshape(massebx,(/iip1,jmp1,llm/)))
      call WriteField_p('masseby',reshape(masseby,(/iip1,jjm,llm/)))
      call WriteField_p('massebxy',reshape(massebxy,(/iip1,jjm,llm/)))
      call WriteField_p('pbaru',reshape(pbaru,(/iip1,jmp1,llm/)))
      call WriteField_p('pbarv',reshape(pbarv,(/iip1,jjm,llm/)))
      call WriteField_p('dteta',reshape(dteta,(/iip1,jmp1,llm/)))
      call WriteField_p('convm',reshape(convm,(/iip1,jmp1,llm/)))
c$OMP END MASTER
c$OMP BARRIER
#endif      

c$OMP BARRIER
c$OMP MASTER
      ijb=ij_begin
      ije=ij_end
      ! compute pressure variation due to mass convergence
      DO ij =ijb, ije
         dp( ij ) = convm( ij,1 ) / airesurg( ij )
      ENDDO
c$OMP END MASTER
c$OMP BARRIER
c$OMP FLUSH

      ! compute vertical velocity w()
      CALL vitvert_p ( convm  , w                                  )
      ! compute potential vorticity vorpot()
      CALL tourpot_p ( vcov   , ucov  , massebxy  , vorpot         )
      ! compute rotation induced du() and dv()
      CALL dudv1_p   ( vorpot , pbaru , pbarv     , du     , dv    )

#ifdef DEBUG_IO      
c$OMP BARRIER
c$OMP MASTER
      call WriteField_p('w',reshape(w,(/iip1,jmp1,llm/)))
      call WriteField_p('vorpot',reshape(vorpot,(/iip1,jjm,llm/)))
      call WriteField_p('du',reshape(du,(/iip1,jmp1,llm/)))
      call WriteField_p('dv',reshape(dv,(/iip1,jjm,llm/)))
c$OMP END MASTER
c$OMP BARRIER
#endif      
      
      ! compute kinetic energy ecin()
      CALL enercin_p ( vcov   , ucov  , vcont     , ucont  , ecin  )
      ! compute Bernouilli function bern()
      CALL bernoui_p ( ip1jmp1, llm   , phi       , ecin   , bern  )
      ! compute and add du() and dv() contributions from Bernouilli and pressure
      CALL dudv2_p   ( tsurpk , pkf   , bern      , du     , dv    )

#ifdef DEBUG_IO
c$OMP BARRIER
c$OMP MASTER
      call WriteField_p('ecin',reshape(ecin,(/iip1,jmp1,llm/)))
      call WriteField_p('bern',reshape(bern,(/iip1,jmp1,llm/)))
      call WriteField_p('du',reshape(du,(/iip1,jmp1,llm/)))
      call WriteField_p('dv',reshape(dv,(/iip1,jjm,llm/)))
      call WriteField_p('pkf',reshape(pkf,(/iip1,jmp1,llm/)))
c$OMP END MASTER
c$OMP BARRIER
#endif
      
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      
      if (pole_nord) ijb=ij_begin
      if (pole_sud) ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
         DO ij=ijb,ije
            ang(ij,l) = ucov(ij,l) + constang(ij)
        ENDDO
      ENDDO
c$OMP END DO

      ! compute vertical advection contributions to du(), dv() and dteta()
      CALL advect_new_p(ang,vcov,teta,w,massebx,masseby,du,dv,dteta) 

C  WARNING probleme de peridocite de dv sur les PC/linux. Pb d'arrondi 
C          probablement. Observe sur le code compile avec pgf90 3.0-1 
      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l = 1, llm
         DO ij = ijb, ije, iip1
           IF( dv(ij,l).NE.dv(ij+iim,l) )  THEN
c         PRINT *,'!!!ATTENTION!!! probleme de periodicite sur vcov',  
c    ,   ' dans caldyn'
c         PRINT *,' l,  ij = ', l, ij, ij+iim,dv(ij+iim,l),dv(ij,l)
          dv(ij+iim,l) = dv(ij,l)
          endif
         enddo
      enddo
c$OMP END DO NOWAIT      
c-----------------------------------------------------------------------
c   Sorties eventuelles des variables de controle:
c   ----------------------------------------------

      IF( conser )  THEN
c ym ---> exige communication collective ( aussi dans advect)
        CALL sortvarc
     $ (itau,ucov,tsurpk,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov)

      ENDIF

      END
