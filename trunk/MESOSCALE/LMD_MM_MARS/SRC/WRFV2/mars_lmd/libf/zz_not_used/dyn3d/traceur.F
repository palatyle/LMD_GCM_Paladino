      SUBROUTINE traceur (iq,iadv,q,teta,pk,w, pbaru, pbarv,dq )
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  F.Forget / A.Armengaud
c   -------
c
c   ....  Modif  P.Le Van  ,  29/04/97  et 01/06/98    ....
c
c   ********************************************************************
c   Transport des traceurs (dont l'humidite q(1)) par differents shemas
c   ********************************************************************
c     iq,q,pbaru et pbarv sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c     iadv = 1 :  Schema de transport type "humidite specifique  LMD"
c                 (shema amont (horiz),moy arithmetique (vertic))
c
c     iadv = 2  :  Shema amont
c
c
c     iadv = 0  : Indice montrant que q n'est pas vraiment un 
c                traceur, mais une variable pronostique attachee
c                (par exemple : une pente du shemas de pente)
c=======================================================================


#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"

c   Arguments:
c   ----------
      INTEGER iq,iadv(iq) 
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q( ip1jmp1,llm,iq), dq( ip1jmp1,llm,iq)
      REAL w(ip1jmp1,llm),teta(ip1jmp1,llm),pk(ip1jmp1,llm)

c   Local:
c   ------

      LOGICAL     qsat
      PARAMETER ( qsat = .TRUE. )

c
cc ******     Si  iadv(iq) =  2   et   qsat = .TRUE.  ,   
cc        ,  on choisit  le shema  amont avec le q sature   ********
c

      EXTERNAL adv_h2o, amont, amont_qsat


       IF( iadv(iq).GT.2 )  THEN
         PRINT *,' Erreur dans le choix de iadv  pour la routine 
     *  traceur .  Corrigez . '
         STOP
       ENDIF

c   ----------------------------------------------------------------
c   Schema de transport type "humidite specifique  LMD" 
c   (shema amont (horiz),moy arithmetique (vertic))
c   ----------------------------------------------------------------
         IF (iadv(iq).EQ.1) THEN
              CALL adv_h2o (iq,iq,q,w, pbaru, pbarv,dq)       
c   ----------------------------------------------------------------
c   Shema amont (dans les 3 dimensions)
c   ----------------------------------------------------------------
         ELSE IF (iadv(iq).EQ.2) THEN
           IF( qsat )   THEN
              CALL amont_qsat ( iq,iq,q,teta,pk,w, pbaru, pbarv,dq )
           ELSE
              CALL amont      ( iq,iq,q,        w, pbaru, pbarv,dq )
           ENDIF
c   ----------------------------------------------------------------

         ENDIF


      RETURN
      END
