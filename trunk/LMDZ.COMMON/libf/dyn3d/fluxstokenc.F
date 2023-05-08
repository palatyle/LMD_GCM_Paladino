!
! $Id: fluxstokenc.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )
#ifdef CPP_IOIPSL
! This routine is designed to work with ioipsl

       USE IOIPSL
c
c     Auteur :  F. Hourdin
c
c
ccc   ..   Modif. P. Le Van  ( 20/12/97 )  ...
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "tracstoke.h"
#include "iniprint.h"

      REAL time_step,t_wrt, t_ops
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL masse(ip1jmp1,llm),teta(ip1jmp1,llm),phi(ip1jmp1,llm)
      REAL phis(ip1jmp1)

      REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),tetac(ip1jmp1,llm),phic(ip1jmp1,llm)

      REAL pbarug(ip1jmp1,llm),pbarvg(iip1,jjm,llm),wg(ip1jmp1,llm)

      REAL pbarvst(iip1,jjp1,llm),zistdyn
	real dtcum

      INTEGER iadvtr,ndex(1) 
      integer nscal
      real tst(1),ist(1),istp(1)
      INTEGER ij,l,irec,i,j,itau
      INTEGER, SAVE :: fluxid, fluxvid,fluxdid
 
      SAVE iadvtr, massem,pbaruc,pbarvc,irec
      SAVE phic,tetac
      logical first
      save first
      data first/.true./
      DATA iadvtr/0/


c AC initialisations
      pbarug(:,:)   = 0.
      pbarvg(:,:,:) = 0.
      wg(:,:)       = 0.
      

      if(first) then

	CALL initfluxsto( 'fluxstoke',
     .  time_step,istdyn* time_step,istdyn* time_step,
     .  fluxid,fluxvid,fluxdid) 
	
	ndex(1) = 0
        call histwrite(fluxid, 'phis', 1, phis, iip1*jjp1, ndex)
        call histwrite(fluxid, 'aire', 1, aire, iip1*jjp1, ndex)
	
	ndex(1) = 0
        nscal = 1
        tst(1) = time_step
        call histwrite(fluxdid, 'dtvr', 1, tst, nscal, ndex)
        ist(1)=istdyn
        call histwrite(fluxdid, 'istdyn', 1, ist, nscal, ndex)
        istp(1)= istphy
        call histwrite(fluxdid, 'istphy', 1, istp, nscal, ndex)
	
	first = .false.

      endif


      IF(iadvtr.EQ.0) THEN
         phic(:,:)=0
         tetac(:,:)=0
         pbaruc(:,:)=0
         pbarvc(:,:)=0
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
            pbaruc(ij,l) = pbaruc(ij,l)/REAL(istdyn)
            tetac(ij,l) = tetac(ij,l)/REAL(istdyn)
            phic(ij,l) = phic(ij,l)/REAL(istdyn)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l)/REAL(istdyn)
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
	write(lunout,*)'ITAU auquel on stoke les fluxmasses',itau
	
	call histwrite(fluxid, 'masse', itau, massem,
     .               iip1*jjp1*llm, ndex)
	
	call histwrite(fluxid, 'pbaru', itau, pbarug,
     .               iip1*jjp1*llm, ndex)
	
	call histwrite(fluxvid, 'pbarv', itau, pbarvg,
     .               iip1*jjm*llm, ndex)
	
        call histwrite(fluxid, 'w' ,itau, wg, 
     .             iip1*jjp1*llm, ndex) 
	
	call histwrite(fluxid, 'teta' ,itau, tetac, 
     .             iip1*jjp1*llm, ndex) 
	
	call histwrite(fluxid, 'phi' ,itau, phic, 
     .             iip1*jjp1*llm, ndex) 
	
C

      ENDIF ! if iadvtr.EQ.istdyn

#else
      write(lunout,*)
     & 'fluxstokenc: Needs IOIPSL to function'
#endif
! of #ifdef CPP_IOIPSL
      RETURN
      END
