










!
! $Id: sortvarc.F 2083 2014-07-09 14:43:31Z emillour $
!
      SUBROUTINE sortvarc
     $(itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time ,
     $ vcov )

      USE control_mod, ONLY: resetvarc 
      USE comconst_mod, ONLY: daysec,dtvr,rad,g,omeg
      USE logic_mod, ONLY: read_start
      USE ener_mod, ONLY: etot,ptot,ztot,stot,ang,
     .			etot0,ptot0,ztot0,stot0,ang0,
     .			rmsdpdt,rmsv
      IMPLICIT NONE


c=======================================================================
c
c   Auteur:    P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   sortie des variables de controle
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      INCLUDE "dimensions.h"
      INCLUDE "paramet.h"
      INCLUDE "comgeom.h"
      INCLUDE "iniprint.h"

c   Arguments:
c   ----------

      INTEGER,INTENT(IN) :: itau
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm)
      REAL,INTENT(IN) :: teta(ip1jmp1,llm)
      REAL,INTENT(IN) :: masse(ip1jmp1,llm)
      REAL,INTENT(IN) :: vcov(ip1jm,llm)
      REAL,INTENT(IN) :: ps(ip1jmp1)
      REAL,INTENT(IN) :: phis(ip1jmp1)
      REAL,INTENT(IN) :: vorpot(ip1jm,llm)
      REAL,INTENT(IN) :: phi(ip1jmp1,llm)
      REAL,INTENT(IN) :: bern(ip1jmp1,llm)
      REAL,INTENT(IN) :: dp(ip1jmp1)
      REAL,INTENT(IN) :: time
      REAL,INTENT(IN) :: pk(ip1jmp1,llm)

c   Local:
c   ------

      REAL vor(ip1jm),bernf(ip1jmp1,llm),ztotl(llm)
      REAL etotl(llm),stotl(llm),rmsvl(llm),angl(llm),ge(ip1jmp1)
      REAL cosphi(ip1jm),omegcosp(ip1jm)
      REAL dtvrs1j,rjour,heure,radsg,radomeg
      REAL massebxy(ip1jm,llm)
      INTEGER  l, ij, imjmp1

      REAL       SSUM
      LOGICAL,SAVE :: firstcal=.true.
      CHARACTER(LEN=*),PARAMETER :: modname="sortvarc"

c-----------------------------------------------------------------------
! Ehouarn: when no initialization fields from file, resetvarc should be
!          set to false
       if (firstcal) then
         if (.not.read_start) then
           resetvarc=.true.
         endif
       endif

       dtvrs1j   = dtvr/daysec
       rjour     = REAL( INT( itau * dtvrs1j ))
       heure     = ( itau*dtvrs1j-rjour ) * 24.
       imjmp1    = iim * jjp1
       IF(ABS(heure - 24.).LE.0.0001 ) heure = 0.
c
       CALL massbarxy ( masse, massebxy )

c   .....  Calcul  de  rmsdpdt  .....

       ge(:)=dp(:)*dp(:)

       rmsdpdt = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)
c
       rmsdpdt = daysec* 1.e-2 * SQRT(rmsdpdt/imjmp1) 

       CALL SCOPY( ijp1llm,bern,1,bernf,1 )
       CALL filtreg(bernf,jjp1,llm,-2,2,.TRUE.,1)

c   .....  Calcul du moment  angulaire   .....

       radsg    = rad /g
       radomeg  = rad * omeg
c
       DO ij=iip2,ip1jm
          cosphi( ij ) = COS(rlatu((ij-1)/iip1+1))
          omegcosp(ij) = radomeg   * cosphi(ij)
       ENDDO

c  ...  Calcul  de l'energie,de l'enstrophie,de l'entropie et de rmsv  .

       DO l=1,llm
          DO ij = 1,ip1jm
             vor(ij)=vorpot(ij,l)*vorpot(ij,l)*massebxy(ij,l)
          ENDDO
          ztotl(l)=(SSUM(ip1jm,vor,1)-SSUM(jjm,vor,iip1))

          DO ij = 1,ip1jmp1
             ge(ij)= masse(ij,l)*(phis(ij)+teta(ij,l)*pk(ij,l)  +
     s        bernf(ij,l)-phi(ij,l))
          ENDDO
          etotl(l) = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO   ij   = 1, ip1jmp1
             ge(ij) = masse(ij,l)*teta(ij,l)
          ENDDO
          stotl(l)= SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO ij=1,ip1jmp1
             ge(ij)=masse(ij,l)*AMAX1(bernf(ij,l)-phi(ij,l),0.)
          ENDDO
          rmsvl(l)=2.*(SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1))

          DO ij =iip2,ip1jm
             ge(ij)=(ucov(ij,l)/cu(ij)+omegcosp(ij))*masse(ij,l) *
     *               cosphi(ij)
          ENDDO
          angl(l) = rad *
     s    (SSUM(ip1jm-iip1,ge(iip2),1)-SSUM(jjm-1,ge(iip2),iip1))
      ENDDO

          DO ij=1,ip1jmp1
            ge(ij)= ps(ij)*aire(ij)
          ENDDO
      ptot  = SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1)
      etot  = SSUM(     llm, etotl, 1 )
      ztot  = SSUM(     llm, ztotl, 1 )
      stot  = SSUM(     llm, stotl, 1 )
      rmsv  = SSUM(     llm, rmsvl, 1 )
      ang   = SSUM(     llm,  angl, 1 )

      IF (firstcal.and.resetvarc) then
         WRITE(lunout,3500) itau, rjour, heure, time
         WRITE(lunout,*) trim(modname),
     &     ' WARNING!!! Recomputing initial values of : '
         WRITE(lunout,*) 'ptot,rmsdpdt,etot,ztot,stot,rmsv,ang'
         WRITE(lunout,*) ptot,rmsdpdt,etot,ztot,stot,rmsv,ang
         etot0 = etot
         ptot0 = ptot
         ztot0 = ztot
         stot0 = stot
         ang0  = ang
      END IF

      ! compute relative changes in etot,... (except if 'reference' values
      ! are zero, which can happen when using iniacademic)
      if (etot0.ne.0) then
        etot= etot/etot0
      else
        etot=1.
      endif
      rmsv= SQRT(rmsv/ptot)
      if (ptot0.ne.0) then
        ptot= ptot/ptot0
      else
        ptot=1.
      endif
      if (ztot0.ne.0) then
        ztot= ztot/ztot0
      else
        ztot=1.
      endif
      if (stot0.ne.0) then
        stot= stot/stot0
      else
        stot=1.
      endif
      if (ang0.ne.0) then
        ang = ang /ang0
      else
        ang=1.
      endif

      firstcal = .false.

      WRITE(lunout,3500) itau, rjour, heure, time
      WRITE(lunout,4000) ptot,rmsdpdt,etot,ztot,stot,rmsv,ang

3500   FORMAT(10("*"),4x,'pas',i7,5x,'jour',f9.0,'heure',f5.1,4x 
     *   ,'date',f14.4,4x,10("*"))
4000   FORMAT(10x,'masse',4x,'rmsdpdt',7x,'energie',2x,'enstrophie'
     * ,2x,'entropie',3x,'rmsv',4x,'mt.ang',/,'GLOB  '
     .  ,f10.6,e13.6,5f10.3/
     * )
      END

