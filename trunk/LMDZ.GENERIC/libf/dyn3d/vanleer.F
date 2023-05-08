      SUBROUTINE vanleer(numvanle,iapp_tracvl,nq,q,pbaru,pbarv ,
     *                     p ,masse, dq ,  iadv1, teta, pk      )
c
      USE comconst_mod, ONLY: dtvr

      IMPLICIT NONE
c
c     Auteurs:   F.Hourdin , P.Le Van, F.Forget, F.Codron  
c
c=======================================================================
c
c       Shema de  Van Leer
c       Version simplifiée pour Mars
c       (L'original bugue a cause des histoires d'eau terrestre !!)
c       FF 2003
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"

c   Arguments:
c   ----------
      INTEGER nq, numvanle, iapp_tracvl, iadv1 
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm),masse(ip1jmp1,llm)
      REAL p( ip1jmp1,llmp1),q( ip1jmp1,llm,nq),dq( ip1jmp1,llm,nq )
      REAL teta( ip1jmp1,llm),pk( ip1jmp1,llm)
c  ..................................................................
c
c   Local:
c   ------

!      EXTERNAL  tracvl,minmaxq, qminimum
      INTEGER ij,l, iq, iapptrac
      REAL finmasse(ip1jmp1,llm), dtvrtrac
cc

        CALL SCOPY( nq * ijp1llm, q, 1, dq, 1 )

c   advection

      CALL tracvl( numvanle,iapp_tracvl,nq,pbaru,pbarv,p , masse,q  ,
     *                      iapptrac, iadv1, teta ,pk              )

cc     ....  Calcul  de deltap  qu'on stocke dans finmasse   ...
c
          DO l = 1, llm
           DO ij = 1, ip1jmp1
             finmasse(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
          ENDDO

          CALL SCOPY   ( ip1jmp1*llm, masse, 1, finmasse,       1 )
          CALL filtreg ( finmasse ,  jjp1,  llm, -2, 2, .TRUE., 1 )
c
c   *****  Calcul de dq pour l'eau , pour le passer a la physique ******

          dtvrtrac = iapp_tracvl * dtvr
c
           DO iq = 1 , nq ! modif special mars : 2 devient nq
            DO l = 1 , llm
             DO ij = 1,ip1jmp1
             dq(ij,l,iq) = ( q(ij,l,iq) - dq(ij,l,iq) ) * finmasse(ij,l)
     *                               /  dtvrtrac
             ENDDO
            ENDDO
           ENDDO


      RETURN
      END
