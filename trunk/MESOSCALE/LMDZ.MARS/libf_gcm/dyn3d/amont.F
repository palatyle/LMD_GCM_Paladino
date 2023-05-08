      SUBROUTINE amont (nq,iq,q,w, pbaru, pbarv,dq)
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van / F.Forget
c   -------
c
c   ********************************************************************
c   Transport d'un traceur q par shema amont 3D
c   ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"
#include "logic.h"
#include "comvert.h"

c   Arguments:
c   ----------
      INTEGER nq,iq  
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm,nq), dq( ip1jmp1,llm,nq )
      REAL w(ip1jmp1,llm)

c   Local:
c   ------

      INTEGER   l,ij

      REAL qbyv( ip1jm,llm ), qbxu( ip1jmp1,llm )
      REAL ww
      REAL dqh( ip1jmp1,llm)

      EXTERNAL     convflu

c     ----------------------------------------------------------------
c     --        calcul transport horiz. pour  dq                    --
c     -                                                              -
c     -         shema ' amont '                                      -
c     ----------------------------------------------------------------
c
      IF( iq.GT.2 )   THEN
       DO   l = 1, llm
        DO ij = 1, ip1jmp1
          dq(ij,l,iq) = 0.
        ENDDO
       ENDDO
            RETURN
      ENDIF
 
c
      IF( forward.OR.leapf )   THEN
 
        DO 10 l = 1, llm
 
        DO 6  ij     = iip2, ip1jm - 1
           IF( pbaru(ij,l).GT.0.  )   THEN
            qbxu( ij,l ) =  pbaru( ij,l )  * q(  ij  ,l,iq)
           ELSE
            qbxu( ij,l ) =  pbaru( ij,l )  * q(ij +1 ,l,iq)
           ENDIF
   6    CONTINUE
 
c     ..... correction  pour  qbxu(iip1,j,l)   .....
c     ...   qbxu(iip1,j,l)= qbxu(1,j,l)  ...
 
CDIR$ IVDEP
        DO  7   ij   = iip1 +iip1, ip1jm, iip1
        qbxu( ij,l ) = qbxu( ij - iim, l )
   7    CONTINUE
 
        DO  8  ij     = 1, ip1jm
           IF( pbarv(ij,l).GT.0.  )   THEN
            qbyv( ij,l ) =  pbarv( ij,l )  * q(ij+iip1, l,iq)
           ELSE
            qbyv( ij,l ) =  pbarv( ij,l )  * q(   ij  , l,iq)
           ENDIF
   8    CONTINUE
c
  10    CONTINUE
 
c     stockage dans  dqh de la convergence horiz.du flux d'humidite  .
c     ------------------------------------------------------------- 
c
                  CALL convflu(qbxu, qbyv, llm, dqh )


c ---------------------------------------------------------------
c   .... calcul des termes d'advection vertic.pour q. Shema amont
c ---------------------------------------------------------------

c     calcul de  - d( q   * w )    qu'on ajoute a  dqh pour calculer dq
c

       DO 20 l = 1,llmm1
c
         DO 11 ij = 1,ip1jmp1
           IF( w(ij,l+1).LT.0.  )   THEN
            ww  =   - w( ij,l+1 ) * q( ij,  l  ,iq)
           ELSE
            ww  =   - w( ij,l+1 ) * q( ij, l +1,iq)
           ENDIF

          dq (ij, l ,iq ) = dqh(ij, l )   -  ww
          dqh(ij,l+1    ) = dqh(ij,l+1)   +  ww
  11     CONTINUE
c
  20   CONTINUE 

c
c       special dq (ij,llm) (qui n'a pas ete rempli ! ) 
c
        DO  ij = 1,ip1jmp1
          dq( ij,llm,iq ) = dqh( ij,llm )
        END DO
           
      END IF

      RETURN
      END
