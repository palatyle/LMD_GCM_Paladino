      SUBROUTINE fluxstoke(pbaru,pbarv,masse,teta,phi,phis)
c
c     Auteur :  F. Hourdin
c
c
ccc   ..   Modif. P. Le Van  ( 20/12/97 )  ...
c
      USE comvert_mod, ONLY: presnivs
      USE comconst_mod, ONLY: dtvr,pi

      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "tracstoke.h"


      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL masse(ip1jmp1,llm),teta(ip1jmp1,llm),phi(ip1jmp1,llm)
      REAL phis(ip1jmp1)

      REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),tetac(ip1jmp1,llm),phic(ip1jmp1,llm)

      REAL pbarug(ip1jmp1,llm),pbarvg(iip1,jjm,llm),wg(ip1jmp1,llm)

      REAL pbarvst(iip1,jjp1,llm)


      INTEGER iadvtr 
      INTEGER ij,l,irec,i,j
 
      SAVE iadvtr, massem,pbaruc,pbarvc,irec
      SAVE phic,tetac
      logical first
      save first
      data first/.true./
      DATA iadvtr/0/

      if(first) then
#ifdef CRAY
         CALL ASSIGN("assign -N ieee -F null f:fluxmass")
#endif
         open(47,file='fluxmass',form='unformatted',
     s        access='direct',recl=4*(6*ijp1llm))
         irec=1
         first=.false.

         open(77,file='fluxmass.ctl',status='unknown',form='formatted')

      endif


      IF(iadvtr.EQ.0) THEN
         CALL initial0(ijp1llm,phic)
         CALL initial0(ijp1llm,tetac)
         CALL initial0(ijp1llm,pbaruc)
         CALL initial0(ijmllm,pbarvc)
      ENDIF

c   accumulation des flux de masse horizontaux
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
            tetac(ij,l) = tetac(ij,l) + teta(ij,l)
            phic(ij,l) = phic(ij,l) + phi(ij,l)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN
         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
      ENDIF

      iadvtr   = iadvtr+1


c   Test pour savoir si on advecte a ce pas de temps
      IF ( iadvtr.EQ.istdyn ) THEN

c    normalisation
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l)/float(istdyn)
            tetac(ij,l) = tetac(ij,l)/float(istdyn)
            phic(ij,l) = phic(ij,l)/float(istdyn)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l)/float(istdyn)
         ENDDO
      ENDDO

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.

        CALL groupe( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

        do l=1,llm
           do j=1,jjm
              do i=1,iip1
                 pbarvst(i,j,l)=pbarvg(i,j,l)
              enddo
           enddo
           do i=1,iip1
              pbarvst(i,jjp1,l)=0.
           enddo
        enddo

         iadvtr=0

         irec=irec+1
         write(47,rec=1) float(irec),dtvr,float(istdyn),
     s    float(iim),float(jjm),float(llm),rlonu,rlonv,rlatu,rlatv
     s    ,aire,phis
         write(47,rec=irec) massem,pbarug,pbarvst,wg,tetac,phic

c   on reinitialise a zero les flux de masse cumules.

      write(77,'(a4,2x,a40)')
     &       'DSET ','^fluxmass'

      write(77,'(a12)') 'UNDEF 1.0E30'
      write(77,'(a5,1x,a40)') 'TITLE ','Titre a voir'
      call formcoord(77,iip1,rlonv,180./pi,.false.,'XDEF')
      call formcoord(77,jjp1,rlatu,180./pi,.true.,'YDEF')
      call formcoord(77,llm,presnivs,1.,.false.,'ZDEF')
      write(77,'(a4,i10,a30)')
     &       'TDEF ',irec,' LINEAR 02JAN1987 1DY '
      write(77,'(a4,2x,i5)') 'VARS',6
      write(77,1000) 'masse',llm,99,'masse    '
      write(77,1000) 'pbaru',llm,99,'pbaru    '
      write(77,1000) 'pbarv',llm,99,'pbarv    '
      write(77,1000) 'w    ',llm,99,'w        '
      write(77,1000) 'teta ',llm,99,'teta     '
      write(77,1000) 'phi  ',llm,99,'phi      '
      write(77,'(a7)') 'ENDVARS'

1000  format(a5,3x,i4,i3,1x,a39)



      ENDIF ! if iadvtr.EQ.istdyn

      RETURN
      END
