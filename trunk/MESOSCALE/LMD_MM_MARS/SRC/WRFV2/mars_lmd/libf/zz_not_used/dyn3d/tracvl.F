      SUBROUTINE tracvl(numvanle,iapp_tracvl,nq,pbaru,pbarv ,
     *                    p, masse , q, iapptrac, iadv1, teta, pk  )
c
c     Auteur :  F. Hourdin
c
c
ccc   ..   Modif. P. Le Van  ( 20/12/97 )  ...
c                 F. Codron     (10/99)

c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom.h"

c     .... Arguments  ....
c
      INTEGER numvanle, nq, iapp_tracvl, iapptrac, iadv1

      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL q(ip1jmp1,llm,nq),masse(ip1jmp1,llm)
      REAL p( ip1jmp1,llmp1 ),teta(ip1jmp1,llm)
      REAL pk(ip1jmp1,llm)

c     ....  var. locales  .....
c
      REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),zdp(ip1jmp1)

      REAL pbarug(ip1jmp1,llm),pbarvg(ip1jm,llm),wg(ip1jmp1,llm)


      INTEGER iadvtr, numvan
      INTEGER ij,l,iq
      REAL zdpmin, zdpmax
      EXTERNAL  minmax
      SAVE iadvtr, massem, pbaruc, pbarvc, numvan
      DATA iadvtr/0/

      numvan = numvanle

      IF(iadvtr.EQ.0) THEN
         CALL initial0(ijp1llm,pbaruc)
         CALL initial0(ijmllm,pbarvc)
      ENDIF

c   accumulation des flux de masse horizontaux
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN

         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
ccc         CALL filtreg ( massem ,jjp1, llm,-2, 2, .TRUE., 1 )
c
      ENDIF

      iadvtr   = iadvtr+1
      iapptrac = iadvtr


c   Test pour savoir si on advecte a ce pas de temps
      IF ( iadvtr.EQ.iapp_tracvl ) THEN

cc   ..  Modif P.Le Van  ( 20/12/97 )  ....
cc

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.

        CALL groupe( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )


c  test sur l'eventuelle creation de valeurs negatives de la masse
         DO l=1,llm-1
            DO ij = iip2+1,ip1jm
              zdp(ij) = pbarug(ij-1,l) - pbarug(ij,l)
     s                  - pbarvg(ij-iip1,l) + pbarvg(ij,l)
     s                  + wg(ij,l+1) - wg(ij,l)
            ENDDO
            CALL SCOPY( jjm -1 ,zdp(iip1+iip1),iip1,zdp(iip2),iip1 )
            DO ij = iip2,ip1jm
               zdp(ij)= zdp(ij)*dtvr/ massem(ij,l) 
            ENDDO 


            CALL minmax ( ip1jm-iip1, zdp(iip2), zdpmin,zdpmax )

            IF(MAX(ABS(zdpmin),ABS(zdpmax)).GT.0.5) THEN
            PRINT*,'WARNING DP/P l=',l,'  MIN:',zdpmin,
     s        '   MAX:', zdpmax
            ENDIF

         ENDDO

c   Advection proprement dite.
c
c   test sur iadv1 pour le schema de vapeur d'eau
c
         IF (numvanle.EQ.1.AND.iadv1.EQ.4) THEN
           CALL vlspltqs( q(1,1,1), 2., massem, wg ,
     *                 pbarug,pbarvg,dtvr,p,pk,teta )
           numvan = 2
         ENDIF

         DO iq = numvan, nq
          CALL vlsplt( q(1,1,iq), 2. ,massem,wg,pbarug,pbarvg,dtvr )
         ENDDO
         iadvtr=0

c   on reinitialise a zero les flux de masse cumules.

      ENDIF ! if iadvtr.EQ.iapp_tracvl

      RETURN
      END
