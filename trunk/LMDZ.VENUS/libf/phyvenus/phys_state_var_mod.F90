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
      USE turb_mod
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
! composition for upper atmosphere
      REAL,ALLOCATABLE,SAVE :: co2vmr_gcm(:,:)   
!$OMP THREADPRIVATE(co2vmr_gcm)
      REAL,ALLOCATABLE,SAVE :: ovmr_gcm(:,:)   
!$OMP THREADPRIVATE(ovmr_gcm)
      REAL,ALLOCATABLE,SAVE :: n2vmr_gcm(:,:)   
!$OMP THREADPRIVATE(n2vmr_gcm)
      REAL,ALLOCATABLE,SAVE :: nvmr_gcm(:,:)   
!$OMP THREADPRIVATE(nvmr_gcm)
      REAL,ALLOCATABLE,SAVE :: covmr_gcm(:,:)   
!$OMP THREADPRIVATE(covmr_gcm)

! photochemistry and microphysics
      REAL,ALLOCATABLE,SAVE :: d_tr_chem(:,:,:), d_tr_sed(:,:,:)   
!$OMP THREADPRIVATE(d_tr_chem,d_tr_sed)
      INTEGER,ALLOCATABLE,SAVE :: iter(:,:)   
!$OMP THREADPRIVATE(iter)

! Tendencies due to radiative scheme   [K/s]
! heat : chauffage solaire
! heat0: chauffage solaire ciel clair
! cool : refroidissement infrarouge
! cool0 : refroidissement infrarouge ciel clair
      REAL,ALLOCATABLE,SAVE :: heat(:,:)   
!$OMP THREADPRIVATE(heat)
      REAL,ALLOCATABLE,SAVE :: heat0(:,:)
!$OMP THREADPRIVATE(heat0)
      REAL,ALLOCATABLE,SAVE :: cool(:,:)
!$OMP THREADPRIVATE(cool)
      REAL,ALLOCATABLE,SAVE :: cool0(:,:)
!$OMP THREADPRIVATE(cool0)
     REAL,ALLOCATABLE,SAVE :: dtsw(:,:)   
!$OMP THREADPRIVATE(dtsw)
     REAL,ALLOCATABLE,SAVE :: dtlw(:,:)   
!$OMP THREADPRIVATE(dtlw)
     REAL,ALLOCATABLE,SAVE :: d_t_rad(:,:),d_t_euv(:,:)   
!$OMP THREADPRIVATE(d_t_rad,d_t_euv)
     REAL,ALLOCATABLE,SAVE :: d_t_nirco2(:,:),d_t_nlte(:,:)
!$OMP THREADPRIVATE(d_t_nirco2,d_t_nlte)

! Case for Newtonian cooling (physideal=.true.)
     REAL,ALLOCATABLE,SAVE :: zt_eq(:,:)   
!$OMP THREADPRIVATE(zt_eq)

! Fluxes (W/m2)
! sollwdown : downward LW flux at surface
! sollwdownclr : downward CS LW flux at surface
! toplwdown : downward CS LW flux at TOA
! toplwdownclr : downward CS LW flux at TOA
      REAL,ALLOCATABLE,SAVE :: swnet(:,:)   
!$OMP THREADPRIVATE(swnet)
      REAL,ALLOCATABLE,SAVE :: lwnet(:,:)   
!$OMP THREADPRIVATE(lwnet)
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
CONTAINS

!======================================================================
SUBROUTINE phys_state_var_init(nqmax)

IMPLICIT NONE
#include "dimsoil.h"

      integer :: nqmax

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

      ALLOCATE(ovmr_gcm(klon,klev),covmr_gcm(klon,klev))
      ALLOCATE(co2vmr_gcm(klon,klev),n2vmr_gcm(klon,klev))
      ALLOCATE(nvmr_gcm(klon,klev))

      ALLOCATE(d_tr_chem(klon,klev,nqmax),d_tr_sed(klon,klev,nqmax))
      ALLOCATE(iter(klon,klev))

      ALLOCATE(heat(klon,klev), heat0(klon,klev)) 
      ALLOCATE(cool(klon,klev), cool0(klon,klev))
      ALLOCATE(dtlw(klon,klev),dtsw(klon,klev))
      ALLOCATE(d_t_rad(klon,klev),d_t_euv(klon,klev))
      ALLOCATE(d_t_nirco2(klon,klev),d_t_nlte(klon,klev))
      ALLOCATE(zt_eq(klon,klev))

      ALLOCATE(swnet(klon,klev+1), lwnet(klon,klev+1)) 
      ALLOCATE(topsw(klon), toplw(klon))
      ALLOCATE(solsw(klon), sollw(klon))
      ALLOCATE(radsol(klon))  ! bilan radiatif au sol calcule par code radiatif
      ALLOCATE(sollwdown(klon), sollwdownclr(klon))
      ALLOCATE(toplwdown(klon), toplwdownclr(klon))
      ALLOCATE(topsw0(klon),toplw0(klon),solsw0(klon),sollw0(klon))
      ALLOCATE(dlw(klon), fder(klon))
      allocate(sens(klon))      
      allocate(q2(klon,klev+1))
      allocate(l0(klon))
      allocate(wstar(klon))
      allocate(yustar(klon))
      allocate(tstar(klon))
      allocate(hfmax_th(klon))
      allocate(zmax_th(klon))

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

      deallocate(ovmr_gcm,covmr_gcm)
      deallocate(co2vmr_gcm,n2vmr_gcm)
      deallocate(nvmr_gcm)

      deallocate(d_tr_chem,d_tr_sed)
      deallocate(iter)

      deallocate(heat, heat0) 
      deallocate(cool, cool0)
      deallocate(dtlw,dtsw)
      deallocate(d_t_rad,d_t_euv)
      deallocate(d_t_nirco2,d_t_nlte)
      deallocate(zt_eq)

      deallocate(swnet, lwnet) 
      deallocate(solsw, sollw, radsol)
      deallocate(topsw, toplw)
      deallocate(sollwdown, sollwdownclr)
      deallocate(toplwdown, toplwdownclr)
      deallocate(topsw0,toplw0,solsw0,sollw0)
      deallocate(dlw, fder)
      deallocate(sens)
      deallocate(q2)
      deallocate(l0)
      deallocate(wstar)
      deallocate(yustar)
      deallocate(tstar)
      deallocate(hfmax_th)
      deallocate(zmax_th)
END SUBROUTINE phys_state_var_end

      END MODULE phys_state_var_mod
