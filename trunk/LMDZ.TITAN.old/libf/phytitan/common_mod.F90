      MODULE common_mod
! Variables sauvegardees - anciens common
!======================================================================
!
!
!======================================================================
! Declaration des variables
      USE dimphy
      implicit none
#include "dimensions.h"

! ancien rnuabar
      REAL,    ALLOCATABLE, SAVE :: rmcbar(:,:),xfbar(:,:,:)
      INTEGER, ALLOCATABLE, SAVE :: ncount(:,:)
!$OMP THREADPRIVATE(rmcbar,xfbar)
!$OMP THREADPRIVATE(ncount)

! Common de aerprod.h... aerprod et htoh2
! Etaient utilises par calchim et la microphysique pour echanger 
! des variables pour la production des aerosols par la chimie
! et la chimie heterogene de H a la surface des aerosols
! Ca a disparu de la microphysique... 
! A remettre si on veut reutiliser cette fonctionnalite

! ancien aerprod
      INTEGER, SAVE :: utilaer(16)
      REAL,    ALLOCATABLE, SAVE :: maer(:,:,:),prodaer(:,:,:)
      REAL,    ALLOCATABLE, SAVE :: csn(:,:,:),csh(:,:,:)
!$OMP THREADPRIVATE(utilaer)
!$OMP THREADPRIVATE(maer,prodaer,csn,csh)
! ancien htoh2
      REAL,    ALLOCATABLE, SAVE :: psurfhaze(:,:)
!$OMP THREADPRIVATE(psurfhaze)

! ancien titan_for.h
      INTEGER, PARAMETER :: NLEV=llm+70,NC=44,ND=54,NR=377
!$OMP THREADPRIVATE(NLEV,NC,ND,NR)
!!!  doivent etre en accord avec titan.h
! pour l'UV (650 niveaux de 2 km)
      INTEGER, PARAMETER :: NLRT=650

! ancien diagmuphy.h
!     toutes les variables
!     diagnostiques sorties de la microphysique.

! ---- flux de glace (1:CH4 / 2:C2H6 / 3:C2H2)
      REAL, ALLOCATABLE, SAVE :: flxesp_i(:,:,:)
!$OMP THREADPRIVATE(flxesp_i)
! ---- taux sedimentation gouttes, aerosols sec
      REAL, ALLOCATABLE, SAVE :: tau_drop(:,:),tau_aer(:,:,:)
!$OMP THREADPRIVATE(tau_drop,tau_aer)
! ---- Production de glace (negatif si disparition)
      REAL, ALLOCATABLE, SAVE :: solesp(:,:,:)
!$OMP THREADPRIVATE(solesp)
! ---- Evaporation CH4
      REAL, ALLOCATABLE, SAVE :: evapch4(:)
!$OMP THREADPRIVATE(evapch4)
! ---- occurences des nuages
      REAL, ALLOCATABLE, SAVE :: occcld_m(:,:,:)
!$OMP THREADPRIVATE(occcld_m)
! ---- occcld sert a obtenir les opacités/extinction des nuages (proxy)
      REAL, ALLOCATABLE, SAVE :: occcld(:,:)
!$OMP THREADPRIVATE(occcld)
! ---- saturation CH4,C2H6,C2H2
      REAL, ALLOCATABLE, SAVE :: satch4(:,:),satc2h6(:,:),satc2h2(:,:)
!$OMP THREADPRIVATE(satch4,satc2h6,satc2h2)
! ---- precipitations (CH4, C2H6, C2H2, noyaux, aerosols)
      REAL, ALLOCATABLE, SAVE :: precip(:,:)
!$OMP THREADPRIVATE(precip)
! ---- rayon moyen des gouttes
      REAL, ALLOCATABLE, SAVE :: rmcloud(:,:)
!$OMP THREADPRIVATE(rmcloud)

! Anciently /TAUD/
      REAL, ALLOCATABLE, SAVE :: TauHID(:,:,:) ! cumulative Haze   IR  opacity
      REAL, ALLOCATABLE, SAVE :: TauCID(:,:,:) ! cumulative Clouds IR  opacity
      REAL, ALLOCATABLE, SAVE :: TauGID(:,:,:) ! cumulative Gas    IR  opacity
      REAL, ALLOCATABLE, SAVE :: TauHVD(:,:,:) ! cumulative Haze   Vis opacity
      REAL, ALLOCATABLE, SAVE :: TauCVD(:,:,:) ! cumulative Clouds Vis opacity
      REAL, ALLOCATABLE, SAVE :: TauGVD(:,:,:) ! cumulative Gas    Vis opacity
!$OMP THREADPRIVATE(TauHID,TauCID,TauGID)
!$OMP THREADPRIVATE(TauHVD,TauCVD,TauGVD)
! besoin en plus en l'absence de racommon_h
      INTEGER,PARAMETER :: NSPECI=46,NSPECV=24

CONTAINS

!======================================================================
SUBROUTINE common_init
use dimphy
IMPLICIT NONE
#include "microtab.h"

      ALLOCATE(rmcbar(klon,klev),xfbar(klon,klev,4))
      ALLOCATE(ncount(klon,klev))

      ALLOCATE(maer(klon,klev,4),prodaer(klon,klev,4))
      ALLOCATE(csn(klon,klev,4),csh(klon,klev,4))
      ALLOCATE(psurfhaze(klon,klev))

      ALLOCATE(flxesp_i(klon,klev,3),tau_drop(klon,klev))
      ALLOCATE(tau_aer(klon,klev,nrad),solesp(klon,klev,3))
      ALLOCATE(evapch4(klon),occcld_m(klon,klev,12))
      ALLOCATE(occcld(klon,klev),satch4(klon,klev))
      ALLOCATE(satc2h6(klon,klev),satc2h2(klon,klev))
      ALLOCATE(precip(klon,5),rmcloud(klon,klev))

      ALLOCATE(TauHID(klon,klev,NSPECI))
      ALLOCATE(TauCID(klon,klev,NSPECI))
      ALLOCATE(TauGID(klon,klev,NSPECI))
      ALLOCATE(TauHVD(klon,klev,NSPECV))
      ALLOCATE(TauCVD(klon,klev,NSPECV))
      ALLOCATE(TauGVD(klon,klev,NSPECV))

END SUBROUTINE common_init

!======================================================================
      END MODULE common_mod
