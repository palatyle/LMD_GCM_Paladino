!
! $Header$
!
      SUBROUTINE histo_o500_pctau(nbreg,pct_ocean,w,histo,histoW,nhisto)
      USE dimphy
      IMPLICIT none

      INTEGER :: ij, k, l, nw
      INTEGER :: nreg, nbreg
cym#include "dimensions.h"
cym#include "dimphy.h"
      INTEGER, PARAMETER :: kmax=8, lmax=8
      INTEGER, PARAMETER :: kmaxm1=kmax-1, lmaxm1=lmax-1
      INTEGER, PARAMETER :: iwmax=40

      INTEGER, dimension(klon) :: iw
      REAL, dimension(klon) :: w
      REAL, PARAMETER :: wmin=-200., pas_w=10.
      REAL, dimension(kmaxm1,lmaxm1,iwmax,nbreg) :: histoW, nhisto
      REAL, dimension(klon,kmaxm1,lmaxm1) :: histo

!     LOGICAL, dimension(klon,nbreg) :: pct_ocean
      INTEGER, dimension(klon,nbreg) :: pct_ocean

! initialisation
      histoW(:,:,:,:)=0.      
      nhisto(:,:,:,:)=0.   
!   
!calcul de l'histogramme de chaque regime dynamique
      DO nreg=1, nbreg
       DO ij=1, klon
        iw(ij) = int((w(ij)-wmin)/pas_w) +1
c       IF(pct_ocean(ij,nreg)) THEN
c       IF(pct_ocean(ij,nreg).EQ.1) THEN
         IF(iw(ij).GE.1.AND.iw(ij).LE.iwmax) THEN 
          DO l=1, lmaxm1
           DO k=1, kmaxm1
            IF(histo(ij,k,l).GT.0.) THEN
             histoW(k,l,iw(ij),nreg) = histoW(k,l,iw(ij),nreg) 
     &       + histo(ij,k,l)*pct_ocean(ij,nreg)
             nhisto(k,l,iw(ij),nreg)= nhisto(k,l,iw(ij),nreg) + 
     &       pct_ocean(ij,nreg)
            ENDIF
           ENDDO !k
          ENDDO !l
c        ELSE IF (iw(ij).LE.0.OR.iw(ij).GT.iwmax) THEN !iw
c         PRINT*,'ij,iw=',ij,iw(ij)
         ENDIF !iw
c       ENDIF !pct_ocean
       ENDDO !ij
!normalisation
       DO nw=1, iwmax
        DO l=1, lmaxm1
         DO k=1, kmaxm1
          IF(nhisto(k,l,nw,nreg).NE.0.) THEN
           histoW(k,l,nw,nreg) = 100.*histoW(k,l,nw,nreg)
     &     /nhisto(k,l,nw,nreg)
c          PRINT*,'k,l,nw,nreg,histoW',k,l,nw,nreg,
c    &     histoW(k,l,nw,nreg)
          ENDIF
         ENDDO !k
        ENDDO !l
       ENDDO !nw
      ENDDO !nreg

      RETURN
      END
