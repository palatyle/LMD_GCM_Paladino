      SUBROUTINE SFLUXV(IPRINT,IG,dist_sol,falbe,icld)

      use dimphy
      USE TGMDAT_MOD, ONLY: UBARI,UBARV,UBAR0
      USE TGMDAT_MOD, ONLY: F0PI
      IMPLICIT NONE
#include "dimensions.h"
#include "comorbit.h"

c  ASTUCE POUR EVITER klon... EN ATTENDANT MIEUX
      INTEGER   ngrid
      PARAMETER (ngrid=(jjm-1)*iim+2)  ! = klon
c
      INTEGER IG,IPRINT,icld
      real dist_sol,falbe(ngrid)

      INTEGER NLAYER,NLEVEL,NSPECV,NSPC1V
      PARAMETER (NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECV=24,NSPC1V=25)
      INTEGER NT,NTERM(NSPECV),J,K

      REAL FUW(NLEVEL),FDW(NLEVEL)
      REAL DT0(NLAYER),T0(NLEVEL),WB0(NLAYER),CO0(NLAYER)
      REAL BTOP, BSURF
      REAL ATERM(4,NSPECV),BTERM(4,NSPECV)
      REAL PEXPON(NSPECV), SOLARF(NSPECV)
      REAL  DTAUV(ngrid,NLAYER,NSPECV,4)
     &     ,TAUV (ngrid,NLEVEL,NSPECV,4)
     &     ,WBARV(ngrid,NLAYER,NSPECV,4) 
     &     ,COSBV(ngrid,NLAYER,NSPECV,4)
     &     ,DTAUVP(ngrid,NLAYER,NSPECV,4)
     &     ,TAUVP(ngrid,NLEVEL,NSPECV,4)
     &     ,WBARVP(ngrid,NLAYER,NSPECV,4)
     &     ,COSBVP(ngrid,NLAYER,NSPECV,4)
      REAL BWNV(NSPC1V),WNOV(NSPECV)
     &     ,DWNV(NSPECV),WLNV(NSPECV)
      REAL  FNETV(ngrid,NLEVEL),     
     &      FUPV(ngrid,NLEVEL,NSPECV),
     &      FDV(ngrid,NLEVEL,NSPECV),
     &      FMNETV(ngrid,NLEVEL),
     &      FMUPV(NLEVEL),FMDV(NLEVEL)

      COMMON /VISGAS/SOLARF,NTERM,PEXPON,
     &         ATERM,BTERM

      COMMON /OPTICV/  DTAUV
     &                ,TAUV 
     &                ,WBARV 
     &                ,COSBV
     &                ,DTAUVP
     &                ,TAUVP
     &                ,WBARVP
     &                ,COSBVP

      COMMON /SPECTV/ BWNV,WNOV
     &               ,DWNV,WLNV

      COMMON /FLUXvV/ FNETV,     
     &               FUPV,
     &               FDV,
     &               FMNETV


* ON NE FAIT PAS LE CALCUL POUR TOUS LES IG EN MEME TEMPS
* IG EST EN ARGUMENT...et SFLUXV EST APPELLEE NGRIDMX FOIS!

C ZERO THE NET FLUXES
      DO 212 J=1,NLEVEL
      FNETV(ig,J)=-0.
      FMNETV(ig,J)=-0.
  212 CONTINUE
C
C WE NOW ENTER A MAJOR LOOP OVER SPECRAL INTERVALS IN THE VISIBLE
C AND OVER THE HORIZONTAL GRIDS
C TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL
C
C ***************************************************************


      DO 500 K=1,NSPECV          ! #2
C ZERO THE SPECTRAL FLUXES IN ANTCIPATION OF SUMMING OVER NTERMS

       DO 214 J=1,NLEVEL         ! #3
       FUPV(ig,J,K)=0.
       FDV(ig,J,K)=0.
 214   CONTINUE
C
C SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
      F0PI=SOLARF(K)*(p_elips/dist_sol)**2.
      BTOP=0.0
C
C LOOP OVER THE NTERMS BEGINING HERE
      DO 912 NT=1,NTERM(K)
      IF (ICLD.eq.1) THEN
        BSURF=0.+ falbe(ig)*UBAR0*F0PI*EXP(-TAUV(ig,NLEVEL,K,NT)/UBAR0)
      ELSE
        BSURF=0.+ falbe(ig)*UBAR0*F0PI*EXP(-TAUVP(ig,NLEVEL,K,NT)/UBAR0)
      ENDIF
C
C* WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C  CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
C
C FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO
C RETURN FLUXES FOR A GIVEN NT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C  USE DT0,T0,WB0,CO0 INSTEAD OF DTAUV(ig,1,K,NT)..etc...

       IF (ICLD.EQ.1) THEN
         DO  J=1,NLAYER        
           DT0(J)=DTAUV(ig,J,K,NT)
           T0(J) =TAUV(ig,J,K,NT)
           WB0(J)=WBARV(ig,J,K,NT)
           CO0(J)=COSBV(ig,J,K,NT)
         ENDDO
         T0(NLEVEL)=TAUV(ig,NLEVEL,K,NT)
       ELSE
         DO  J=1,NLAYER        
           DT0(J)=DTAUVP(ig,J,K,NT)
           T0(J) =TAUVP(ig,J,K,NT)
           WB0(J)=WBARVP(ig,J,K,NT)
           CO0(J)=COSBVP(ig,J,K,NT)
         ENDDO
         T0(NLEVEL)=TAUVP(ig,NLEVEL,K,NT)

       ENDIF

c       PRINT*,'entree gfluxv #: ',ig,K
c        write(*,*) (DT0(J),J=1,NLAYER)
c        print*,'---'
c        write(*,*) (T0(J),J=1,NLEVEL)
c        print*,'---'
c        write(*,*) (WB0(J),J=1,NLAYER)
c        print*,'---'
c        write(*,*) (CO0(J),J=1,NLAYER)
c        print*,'UBAR0 ',UBAR0
c      print*,NLEVEL,WNOV(K),F0PI,falbe(ig),BTOP,BSURF
       FUW = 0.0
       FDW = 0.0
       FMUPV=0.0
       FMDV= 0.0
       
      CALL GFLUXV(NLEVEL,WNOV(K),DT0,T0,
     & WB0,CO0,F0PI,falbe(ig),BTOP,BSURF,FUW,FDW,FMUPV,
     &    FMDV,IPRINT)
c       PRINT*,'sortie gfluxv #: ',ig,K
c        print*,'UBAR0 ',UBAR0

C NOW CALCULTE THE CUMULATIVE VISIBLE NET FLUX

             DO 300 J=1,NLEVEL         !<------------
      FMNETV(ig,J)=FMNETV(ig,J)+( FMUPV(J)-FMDV(J) )*ATERM(NT,K)
      FNETV(ig,J)=FNETV(ig,J)+( FUW(J)-FDW(J) )*ATERM(NT,K)

C AND THE SPECTRAL FLUXES SUMMED OVER THE NTERMS
      FUPV(ig,J,K)=FUPV(ig,J,K)+FUW(J)*ATERM(NT,K)
      FDV(ig,J,K)=FDV(ig,J,K)+FDW(J)*ATERM(NT,K)
  300         CONTINUE                !<--------------
C
C
  912 CONTINUE
  500 CONTINUE

C *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****
      RETURN
      END
