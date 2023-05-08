      SUBROUTINE amont_qsat (nq,iq,q,teta,pk,w, pbaru, pbarv,dq)
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van / F.Forget
c   -------
c
c   ********************************************************************
c   Transport d'un traceur q par shema amont 3D  avec contrainte thermi-
c           -que
c   ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"
#include "logic.h"
#include "comvert.h"
#include "comconst.h"
#include "comgeom.h"
#include "serre.h"

c   Arguments:
c   ----------
      INTEGER nq,iq  
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm,nq), dq( ip1jmp1,llm,nq )
      REAL teta(ip1jmp1,llm),pk(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)

c   Local:
c   ------

      INTEGER   l,ij

      REAL qbyv( ip1jm,llm ), qbxu( ip1jmp1,llm )
      REAL ww
      REAL dqh( ip1jmp1,llm)

      REAL ptarg,pdelarg,foeew,zdelta
      REAL rtt,retv,r3les,r3ies,r4les,r4ies,unskappa,play
      REAL tempe(ip1jmp1),qsat(ip1jmp1,llm),r2es

      EXTERNAL     convflu

       FOEEW ( PTARG,PDELARG ) = EXP (
     *          (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT)
     * / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )
C


c     ----------------------------------------------------------------
c     --        calcul transport horiz. pour  dq                    --
c     -                                                              -
c     -         shema ' amont '                                      -
c     ----------------------------------------------------------------
 
        r2es  = 380.11733 
        r3les = 17.269
        r3ies = 21.875
        r4les = 35.86
        r4ies = 7.66
        retv = 0.6077667
        rtt  = 273.16
 
        unskappa =  1./ kappa

        DO l = 1, llm
         DO ij = 1, ip1jmp1
          tempe(ij) = teta(ij,l) * pk(ij,l) /cpp
         ENDDO
         DO ij = 1, ip1jmp1
          zdelta = MAX( 0., SIGN(1., rtt - tempe(ij)) )
          play   =  preff * (pk(ij,l)/cpp) ** unskappa
          qsat(ij,l) = MIN(0.5, r2es* FOEEW(tempe(ij),zdelta) / play )
          qsat(ij,l) = qsat(ij,l) / ( 1. - retv * qsat(ij,l) )
         ENDDO
        ENDDO

        DO 10 l = 1, llm

        DO 6  ij     = iip2, ip1jm - 1
         IF( pbaru(ij,l).GT.0  )  THEN
          qbxu( ij,l ) = pbaru( ij,l ) * MIN( q(ij,l,iq),qsat(ij +1,l))
         ELSE
          qbxu( ij,l ) = pbaru( ij,l ) * MIN( q(ij+1,l,iq),qsat(ij,l) )
         ENDIF
       
   6    CONTINUE
 
c     ..... correction  pour  qbxu(iip1,j,l)   .....
c     ...   qbxu(iip1,j,l)= qbxu(1,j,l)  ...
 
CDIR$ IVDEP
        DO  7   ij   = iip1 +iip1, ip1jm, iip1
        qbxu( ij,l ) = qbxu( ij - iim, l )
   7    CONTINUE
 
        DO  8  ij     = 1, ip1jm
         IF( pbarv(ij,l).GT.0  )  THEN
          qbyv(ij,l) = pbarv( ij,l ) * MIN( q(ij+iip1,l,iq),qsat(ij,l) )
         ELSE
          qbyv(ij,l) = pbarv( ij,l ) * MIN( q(ij,l,iq),qsat(ij+iip1,l) )
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

c calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dqh pour calculer dq
c

       DO 20 l = 1,llmm1
c
         DO 11 ij = 1,ip1jmp1
           IF( w(ij,l+1).LT.0. )   THEN
             ww =  -w(ij,l+1) * MIN( q(ij, l ,iq),qsat(ij ,l+1) )
           ELSE
             ww =  -w(ij,l+1) * MIN( q(ij,l+1,iq),qsat(ij , l ) )
           ENDIF

cc     .... Modif . P. Le Van  (26/09/95) ....

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
           
      RETURN
      END
