!
! $Id: sortvarc0.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE sortvarc0
     $(itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time ,
     $ vcov)

      USE comconst_mod, ONLY: daysec,dtvr,rad,g,omeg
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0,rmsv,rmsdpdt

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

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c   Arguments:
c   ----------

      INTEGER itau
      REAL ucov(ip1jmp1,llm),teta(ip1jmp1,llm),masse(ip1jmp1,llm)
      REAL vcov(ip1jm,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)
      REAL vorpot(ip1jm,llm)
      REAL phi(ip1jmp1,llm),bern(ip1jmp1,llm)
      REAL dp(ip1jmp1)
      REAL time
      REAL pk(ip1jmp1,llm)

c   Local:
c   ------

      REAL vor(ip1jm),bernf(ip1jmp1,llm),ztotl(llm)
      REAL etotl(llm),stotl(llm),rmsvl(llm),angl(llm),ge(ip1jmp1)
      REAL cosphi(ip1jm),omegcosp(ip1jm)
      REAL dtvrs1j,rjour,heure,radsg,radomeg
      REAL rday, massebxy(ip1jm,llm)
      INTEGER  l, ij, imjmp1

      REAL       SSUM
      integer  ismin,ismax

c-----------------------------------------------------------------------

       dtvrs1j   = dtvr/daysec
       rjour     = REAL( INT( itau * dtvrs1j ))
       heure     = ( itau*dtvrs1j-rjour ) * 24.
       imjmp1    = iim * jjp1
       IF(ABS(heure - 24.).LE.0.0001 ) heure = 0.
c
       CALL massbarxy ( masse, massebxy )

c   .....  Calcul  de  rmsdpdt  .....

       ge=dp*dp

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
      ptot0  = SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1)
      etot0  = SSUM(     llm, etotl, 1 )
      ztot0  = SSUM(     llm, ztotl, 1 )
      stot0  = SSUM(     llm, stotl, 1 )
      rmsv   = SSUM(     llm, rmsvl, 1 )
      ang0   = SSUM(     llm,  angl, 1 )

      rday = REAL(INT (time ))
c
      PRINT 3500, itau, rday, heure, time
      PRINT *, ptot0,etot0,ztot0,stot0,ang0

3500   FORMAT(10("*"),4x,'pas',i7,5x,'jour',f5.0,'heure',f5.1,4x 
     *   ,'date',f10.5,4x,10("*"))
      RETURN
      END

