!
! $Id: phys_state_var_mod.F90 1670 2012-10-17 08:42:04Z idelkadi $
!
      MODULE phys_state_var_mod
! Variables sauvegardees pour le startphy.nc
!======================================================================
!
!
!======================================================================
! Declaration des variables
      USE dimphy
!      INTEGER, SAVE :: radpas
!!$OMP THREADPRIVATE(radpas)
!      REAL, SAVE :: dtime
!!$OMP THREADPRIVATE(dtime)

      REAL, ALLOCATABLE, SAVE :: ftsol(:)
!$OMP THREADPRIVATE(ftsol)
      REAL, ALLOCATABLE, SAVE :: ftsoil(:,:)
!$OMP THREADPRIVATE(ftsoil)
      REAL, ALLOCATABLE, SAVE :: falbe(:)
!$OMP THREADPRIVATE(falbe)

!clesphy0 param physiq
!
! Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
!
      REAL, ALLOCATABLE, SAVE :: zmea(:), zstd(:), zsig(:), zgam(:)
!$OMP THREADPRIVATE(zmea, zstd, zsig, zgam)
      REAL, ALLOCATABLE, SAVE :: zthe(:), zpic(:), zval(:)
!$OMP THREADPRIVATE(zthe, zpic, zval)
!     REAL tabcntr0(100)
      REAL, ALLOCATABLE, SAVE :: rugoro(:)
!$OMP THREADPRIVATE(rugoro)
      REAL, ALLOCATABLE, SAVE :: t_ancien(:,:), q_ancien(:,:)
!$OMP THREADPRIVATE(t_ancien, q_ancien)
      REAL, ALLOCATABLE, SAVE :: u_ancien(:,:), v_ancien(:,:)
!$OMP THREADPRIVATE(u_ancien, v_ancien)
      LOGICAL, SAVE :: ancien_ok
!$OMP THREADPRIVATE(ancien_ok)
! pressure level
      REAL,ALLOCATABLE,SAVE :: zuthe(:),zvthe(:)
!$OMP THREADPRIVATE(zuthe,zvthe)
!
! heat : chauffage solaire
! heat0: chauffage solaire ciel clair
! cool : refroidissement infrarouge
! cool0 : refroidissement infrarouge ciel clair
! sollwdown : downward LW flux at surface
! sollwdownclr : downward CS LW flux at surface
! toplwdown : downward CS LW flux at TOA
! toplwdownclr : downward CS LW flux at TOA
! swnet,swdn,lwdn: + downward
! lwnet,swup,lwup: + upward
      REAL,ALLOCATABLE,SAVE :: swnet(:,:),swup(:,:),swdn(:,:)   
!$OMP THREADPRIVATE(swnet,swup,swdn)
      REAL,ALLOCATABLE,SAVE :: lwnet(:,:),lwup(:,:),lwdn(:,:)
!$OMP THREADPRIVATE(lwnet,lwup,lwdn)
      REAL,ALLOCATABLE,SAVE :: heat(:,:)   
!$OMP THREADPRIVATE(heat)
      REAL,ALLOCATABLE,SAVE :: heat0(:,:)
!$OMP THREADPRIVATE(heat0)
      REAL,ALLOCATABLE,SAVE :: cool(:,:)
!$OMP THREADPRIVATE(cool)
      REAL,ALLOCATABLE,SAVE :: cool0(:,:)
!$OMP THREADPRIVATE(cool0)
      REAL,ALLOCATABLE,SAVE :: dtrad(:,:)   
!$OMP THREADPRIVATE(dtrad)
      REAL,ALLOCATABLE,SAVE :: topsw(:), toplw(:)
!$OMP THREADPRIVATE(topsw,toplw)
      REAL, ALLOCATABLE, SAVE :: solsw(:), sollw(:)
!$OMP THREADPRIVATE(solsw, sollw)
      REAL, ALLOCATABLE, SAVE :: radsol(:)
!$OMP THREADPRIVATE(radsol)
      REAL,ALLOCATABLE,SAVE :: sollwdown(:)
!$OMP THREADPRIVATE(sollwdown)
      REAL,ALLOCATABLE,SAVE :: sollwdownclr(:)
!$OMP THREADPRIVATE(sollwdownclr)
      REAL,ALLOCATABLE,SAVE :: toplwdown(:)
!$OMP THREADPRIVATE(toplwdown)
      REAL,ALLOCATABLE,SAVE :: toplwdownclr(:)
!$OMP THREADPRIVATE(toplwdownclr)
      REAL,ALLOCATABLE,SAVE :: topsw0(:),toplw0(:),solsw0(:),sollw0(:)
!$OMP THREADPRIVATE(topsw0,toplw0,solsw0,sollw0)
      REAL,save,allocatable :: dlw(:)  ! derivee infra rouge
      REAL,save,allocatable :: fder(:) ! Derive de flux (sensible et latente) 
!$OMP THREADPRIVATE(dlw,fder)

!
! Parametres pour le cycle du methane:
!
      REAL,save,allocatable :: resch4(:) ! surface reservoir CH4
!$OMP THREADPRIVATE(resch4)

CONTAINS

!======================================================================
SUBROUTINE phys_state_var_init

IMPLICIT NONE
#include "dimsoil.h"

      ALLOCATE(ftsol(klon))            ! temperature de surface
      ALLOCATE(ftsoil(klon,nsoilmx))   ! temperature dans le sol
      ALLOCATE(falbe(klon))            ! albedo

!  Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
!
!zmea(:)   ! orographie moyenne
!zstd(:)   ! deviation standard de l'OESM
!zsig(:)   ! pente de l'OESM
!zgam(:)   ! anisotropie de l'OESM
!zthe(:)   ! orientation de l'OESM
!zpic(:)   ! Maximum de l'OESM
!zval(:)   ! Minimum de l'OESM
!rugoro(:) ! longueur de rugosite de l'OESM
      ALLOCATE(zmea(klon), zstd(klon), zsig(klon), zgam(klon))
      ALLOCATE(zthe(klon), zpic(klon), zval(klon))
      ALLOCATE(rugoro(klon))

      ALLOCATE(t_ancien(klon,klev), q_ancien(klon,klev))
      ALLOCATE(u_ancien(klon,klev), v_ancien(klon,klev))

      ALLOCATE(zuthe(klon),zvthe(klon))
!
      ALLOCATE(swnet(klon,klev+1), lwnet(klon,klev+1)) 
      ALLOCATE(swup(klon,klev+1), lwup(klon,klev+1))
      ALLOCATE(swdn(klon,klev+1), lwdn(klon,klev+1))
      ALLOCATE(heat(klon,klev), heat0(klon,klev)) 
      ALLOCATE(cool(klon,klev), cool0(klon,klev))
      ALLOCATE(dtrad(klon,klev))
      ALLOCATE(topsw(klon), toplw(klon))
      ALLOCATE(solsw(klon), sollw(klon))
      ALLOCATE(radsol(klon))  ! bilan radiatif au sol calcule par code radiatif
      ALLOCATE(sollwdown(klon), sollwdownclr(klon))
      ALLOCATE(toplwdown(klon), toplwdownclr(klon))
      ALLOCATE(topsw0(klon),toplw0(klon),solsw0(klon),sollw0(klon))
      ALLOCATE(dlw(klon), fder(klon))
      
      ALLOCATE(resch4(klon))

END SUBROUTINE phys_state_var_init

!======================================================================
SUBROUTINE phys_state_var_end

IMPLICIT NONE

      deallocate(ftsol, ftsoil, falbe)
      deallocate(zmea, zstd, zsig, zgam)
      deallocate(zthe, zpic, zval)
      deallocate(rugoro, t_ancien, q_ancien)
      deallocate(        u_ancien, v_ancien)
      deallocate(zuthe, zvthe)
      deallocate(swnet, lwnet) 
      deallocate(swup, lwup)
      deallocate(swdn, lwdn)
      deallocate(heat, heat0) 
      deallocate(cool, cool0)
      deallocate(dtrad)
      deallocate(solsw, sollw, radsol)
      deallocate(topsw, toplw)
      deallocate(sollwdown, sollwdownclr)
      deallocate(toplwdown, toplwdownclr)
      deallocate(topsw0,toplw0,solsw0,sollw0)
      deallocate(dlw, fder)

      deallocate(resch4)

END SUBROUTINE phys_state_var_end

      END MODULE phys_state_var_mod
