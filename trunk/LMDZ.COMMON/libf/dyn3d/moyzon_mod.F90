      MODULE moyzon_mod
! Moyennes zonales pour transmission a la physique 
!======================================================================
! Specifique a Titan
!
!======================================================================
! Declaration des variables

      REAL,ALLOCATABLE,SAVE :: zplevbar(:,:),zplaybar(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfibar(:,:),zqfibar(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zphibar(:,:),zphisbar(:)
      REAL,ALLOCATABLE,SAVE :: zpkbar(:,:),ztetabar(:,:)
      REAL,ALLOCATABLE,SAVE :: zzlevbar(:,:),zzlaybar(:,:)

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

      ALLOCATE(zplevbar(klon,llm+1),zplaybar(klon,llm))
      ALLOCATE(zphibar(klon,llm),zphisbar(klon))
      ALLOCATE(ztfibar(klon,llm),zqfibar(klon,llm,nqtot))
      ALLOCATE(zpkbar(klon,llm),ztetabar(klon,llm))
      ALLOCATE(zzlevbar(klon,llm+1),zzlaybar(klon,llm))
#else
    include "iniprint.h"
    write(lunout,*) "moyzon_init Error: should only be called with physics"
    stop
#endif
END SUBROUTINE moyzon_init

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
