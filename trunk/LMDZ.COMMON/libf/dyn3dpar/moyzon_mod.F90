      MODULE moyzon_mod
! Moyennes zonales pour transmission a la physique 
!======================================================================
! Specifique a Titan
!
!======================================================================
! Declaration des variables

      REAL,ALLOCATABLE,SAVE :: zplevbar_mpi(:,:),zplaybar_mpi(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfibar_mpi(:,:),zqfibar_mpi(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zphibar_mpi(:,:),zphisbar_mpi(:)
      REAL,ALLOCATABLE,SAVE :: zpkbar_mpi(:,:),ztetabar_mpi(:,:)

      REAL,ALLOCATABLE,SAVE :: zplevbar(:,:),zplaybar(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfibar(:,:),zqfibar(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zphibar(:,:),zphisbar(:)
      REAL,ALLOCATABLE,SAVE :: zzlevbar(:,:),zzlaybar(:,:)
!$OMP THREADPRIVATE(zplevbar,zplaybar,ztfibar,zqfibar)
!$OMP THREADPRIVATE(zphibar,zphisbar,zzlevbar,zzlaybar)

! pmoy: global averaged pressure...
! tmoy: global averaged temperature...
! put here to be transfered to Titan routines...
! to be changed...
      REAL,ALLOCATABLE,SAVE :: plevmoy(:),playmoy(:)
      REAL,ALLOCATABLE,SAVE :: zlevmoy(:),zlaymoy(:),phimoy(:)
      REAL,ALLOCATABLE,SAVE :: tmoy(:),tetamoy(:),pkmoy(:)
      INTEGER,ALLOCATABLE,SAVE :: klat(:)

CONTAINS

!======================================================================
SUBROUTINE moyzon_init(klon,llm,nqtot)
#ifdef CPP_PHYS
!! This routine needs physics
IMPLICIT NONE
       INTEGER :: klon,llm,nqtot

      ALLOCATE(zplevbar_mpi(klon,llm+1),zplaybar_mpi(klon,llm))
      ALLOCATE(zphibar_mpi(klon,llm),zphisbar_mpi(klon))
      ALLOCATE(ztfibar_mpi(klon,llm),zqfibar_mpi(klon,llm,nqtot))
      ALLOCATE(zpkbar_mpi(klon,llm),ztetabar_mpi(klon,llm))
#endif
END SUBROUTINE moyzon_init

!======================================================================
SUBROUTINE moyzon_init_omp(nlon,llm,nqtot)
#ifdef CPP_PHYS
IMPLICIT NONE

      INTEGER :: nlon,llm,nqtot

      ALLOCATE(zplevbar(nlon,llm+1),zplaybar(nlon,llm))
      ALLOCATE(zphibar(nlon,llm),zphisbar(nlon))
      ALLOCATE(ztfibar(nlon,llm),zqfibar(nlon,llm,nqtot))
      ALLOCATE(zzlevbar(nlon,llm+1),zzlaybar(nlon,llm))
#endif
END SUBROUTINE moyzon_init_omp

!======================================================================
SUBROUTINE moyzon(nlev,var,varbar)

IMPLICIT NONE
#include "dimensions.h"
#include "paramet.h"

      INTEGER :: nlev
      REAL,dimension(iip1,nlev) :: var
      REAL,dimension(nlev)      :: varbar

      INTEGER :: i

      varbar(:) = 0.
      do i=1,iim
        varbar(:)=varbar(:)+var(i,:)/iim
      enddo

      return
END SUBROUTINE moyzon

!======================================================================
      END MODULE moyzon_mod
