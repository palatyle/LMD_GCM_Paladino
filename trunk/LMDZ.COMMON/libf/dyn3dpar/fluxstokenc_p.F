!
! $Id: fluxstokenc_p.F 1454 2010-11-18 12:01:24Z fairhead $
!
      SUBROUTINE fluxstokenc_p(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )
#ifdef CPP_IOIPSL
! This routine is designed to work with ioipsl

       USE IOIPSL
       USE parallel_lmdz
       USE misc_mod
       USE mod_hallo
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

      REAL,SAVE :: pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),tetac(ip1jmp1,llm),phic(ip1jmp1,llm)

      REAL pbarug(ip1jmp1,llm),pbarvg(iip1,jjm,llm),wg(ip1jmp1,llm)

      REAL pbarvst(iip1,jjp1,llm),zistdyn
	real dtcum

      INTEGER iadvtr,ndex(1) 
      integer nscal
      real tst(1),ist(1),istp(1)
      INTEGER ij,l,irec,i,j,itau
      INTEGER,SAVE :: fluxid, fluxvid,fluxdid
 
      SAVE iadvtr, massem,irec
      SAVE phic,tetac
      logical first
      save first
      data first/.true./
      DATA iadvtr/0/
      integer :: ijb,ije,jjb,jje,jjn
      type(Request) :: Req

c AC initialisations
      pbarug(:,:)   = 0.
cym      pbarvg(:,:,:) = 0.
cym      wg(:,:)       = 0.

c$OMP MASTER

      if(first) then

	CALL initfluxsto_p( 'fluxstoke',
     .  time_step,istdyn* time_step,istdyn* time_step,
     .  fluxid,fluxvid,fluxdid) 
	
        ijb=ij_begin
        ije=ij_end
        jjn=jj_nb

	ndex(1) = 0
        call histwrite(fluxid, 'phis', 1, phis(ijb:ije),
     .	               iip1*jjn, ndex)
        call histwrite(fluxid, 'aire', 1, aire(ijb:ije),
     .                 iip1*jjn, ndex)
	
	ndex(1) = 0
        nscal = 1
        
	if (mpi_rank==0) then
          tst(1) = time_step
          call histwrite(fluxdid, 'dtvr', 1, tst, nscal, ndex)
          ist(1)=istdyn
          call histwrite(fluxdid, 'istdyn', 1, ist, nscal, ndex)
          istp(1)= istphy
          call histwrite(fluxdid, 'istphy', 1, istp, nscal, ndex)
	endif
	first = .false.

      endif


      IF(iadvtr.EQ.0) THEN
cym         CALL initial0(ijp1llm,phic)
cym        CALL initial0(ijp1llm,tetac)
cym         CALL initial0(ijp1llm,pbaruc)
cym         CALL initial0(ijmllm,pbarvc)
        ijb=ij_begin
        ije=ij_end
        phic(ijb:ije,1:llm)=0
	tetac(ijb:ije,1:llm)=0
	pbaruc(ijb:ije,1:llm)=0
	
	IF (pole_sud) ije=ij_end-iip1
	pbarvc(ijb:ije,1:llm)=0
      ENDIF

c   accumulation des flux de masse horizontaux
      ijb=ij_begin
      ije=ij_end
      
      DO l=1,llm
         DO ij = ijb,ije
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
            tetac(ij,l) = tetac(ij,l) + teta(ij,l)
            phic(ij,l) = phic(ij,l) + phi(ij,l)
         ENDDO
       ENDDO
      
      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1
	
      DO l=1,llm
         DO ij = ijb,ije
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN
cym         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
        ijb=ij_begin
        ije=ij_end 
	massem(ijb:ije,1:llm)=masse(ijb:ije,1:llm)
      ENDIF

      iadvtr   = iadvtr+1

c$OMP END MASTER
c$OMP BARRIER
c   Test pour savoir si on advecte a ce pas de temps
      IF ( iadvtr.EQ.istdyn ) THEN
c$OMP MASTER
c    normalisation
      ijb=ij_begin
      ije=ij_end 

      DO l=1,llm
         DO ij = ijb,ije
            pbaruc(ij,l) = pbaruc(ij,l)/REAL(istdyn)
            tetac(ij,l) = tetac(ij,l)/REAL(istdyn)
            phic(ij,l) = phic(ij,l)/REAL(istdyn)
         ENDDO
      ENDDO

      ijb=ij_begin
      ije=ij_end 
      if (pole_sud) ije=ij_end-iip1      
      
      DO l=1,llm
          DO ij = ijb,ije
            pbarvc(ij,l) = pbarvc(ij,l)/REAL(istdyn)
         ENDDO
      ENDDO

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.
c$OMP END MASTER
c$OMP BARRIER 
        call Register_Hallo(pbaruc,ip1jmp1,llm,1,1,1,1,Req)
	call Register_Hallo(pbarvc,ip1jm,llm,1,1,1,1,Req)
        call SendRequest(Req)
c$OMP BARRIER
        call WaitRequest(Req)
c$OMP BARRIER
c$OMP MASTER
        CALL groupe_p( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )
        
        jjb=jj_begin
	jje=jj_end
	if (pole_sud) jje=jj_end-1
	
        do l=1,llm
           do j=jjb,jje
              do i=1,iip1
                 pbarvst(i,j,l)=pbarvg(i,j,l)
              enddo
           enddo
	 enddo
	 
	 if (pole_sud) then
           do i=1,iip1
              pbarvst(i,jjp1,l)=0.
           enddo
        endif
      
         iadvtr=0
	write(lunout,*)'ITAU auquel on stoke les fluxmasses',itau
	
        ijb=ij_begin
	ije=ij_end
	jjn=jj_nb
	
	call histwrite(fluxid, 'masse', itau, massem(ijb:ije,:),
     .               iip1*jjn*llm, ndex)
	
	call histwrite(fluxid, 'pbaru', itau, pbarug(ijb:ije,:),
     .               iip1*jjn*llm, ndex)
	
        jjb=jj_begin
	jje=jj_end
	jjn=jj_nb
	if (pole_sud) then
	  jje=jj_end-1
	  jjn=jj_nb-1
	endif
	
	call histwrite(fluxvid, 'pbarv', itau, pbarvg(:,jjb:jje,:),
     .               iip1*jjn*llm, ndex)
	
        ijb=ij_begin
	ije=ij_end
	jjn=jj_nb
	
        call histwrite(fluxid, 'w' ,itau, wg(ijb:ije,:), 
     .             iip1*jjn*llm, ndex) 
	
	call histwrite(fluxid, 'teta' ,itau, tetac(ijb:ije,:), 
     .             iip1*jjn*llm, ndex) 
	
	call histwrite(fluxid, 'phi' ,itau, phic(ijb:ije,:), 
     .             iip1*jjn*llm, ndex) 
	
C
c$OMP END MASTER
      ENDIF ! if iadvtr.EQ.istdyn

#else
      write(lunout,*)
     & 'fluxstokenc: Needs IOIPSL to function'
#endif
! of #ifdef CPP_IOIPSL
      RETURN
      END
