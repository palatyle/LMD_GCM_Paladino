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
      USE dimphy, only : klon,klev
      USE comsoil_h, only : nsoilmx
      use comsaison_h, only: mu0, fract
      use radcommon_h, only: glat
      use slab_ice_h, only : noceanmx
      USE radinc_h, only : L_NSPECTI, L_NSPECTV,naerkind
      use surfdat_h, only: phisfi, albedodat,  &
                        zmea, zstd, zsig, zgam, zthe
      use turb_mod, only: q2,sensibFlux,wstar,ustar,tstar,hfmax_th,zmax_th

      real,allocatable,dimension(:,:),save :: ztprevious ! Previous loop Atmospheric Temperature (K)
! Useful for Dynamical Heating calculation.
      real,allocatable,dimension(:,:),save :: zuprevious
!$OMP THREADPRIVATE(ztprevious,zuprevious)

      real, dimension(:),allocatable,save ::  tsurf                ! Surface temperature (K).
      real, dimension(:,:),allocatable,save ::  tsoil              ! Sub-surface temperatures (K).
      real, dimension(:,:),allocatable,save :: albedo              ! Surface Spectral albedo. By MT2015.
      real, dimension(:),allocatable,save :: albedo_equivalent     ! Spectral Mean albedo.
      real, dimension(:),allocatable,save :: albedo_snow_SPECTV    ! Snow Spectral albedo.
      real, dimension(:),allocatable,save :: albedo_co2_ice_SPECTV ! CO2 ice Spectral albedo.
!$OMP THREADPRIVATE(tsurf,tsoil,albedo,albedo_equivalent,albedo_snow_SPECTV,albedo_co2_ice_SPECTV)

      real,dimension(:),allocatable,save :: albedo_bareground ! Bare Ground Albedo. By MT 2015.
      real,dimension(:),allocatable,save :: rnat              ! Defines the type of the grid (ocean,continent,...). By BC.
!$OMP THREADPRIVATE(albedo_bareground,rnat)

      real,dimension(:),allocatable,save :: emis        ! Thermal IR surface emissivity.
      real,dimension(:,:),allocatable,save :: dtrad     ! Net atmospheric radiative heating rate (K.s-1).
      real,dimension(:),allocatable,save :: fluxrad_sky ! Radiative flux from sky absorbed by surface (W.m-2).
      real,dimension(:),allocatable,save :: fluxrad     ! Net radiative surface flux (W.m-2).
      real,dimension(:),allocatable,save :: capcal      ! Surface heat capacity (J m-2 K-1).
      real,dimension(:),allocatable,save :: fluxgrd     ! Surface conduction flux (W.m-2).
      real,dimension(:,:),allocatable,save :: qsurf     ! Tracer on surface (e.g. kg.m-2).
!$OMP THREADPRIVATE(emis,dtrad,fluxrad_sky,fluxrad,capcal,fluxgrd,qsurf)

      ! FOR DIAGNOSTIC :

      real,dimension(:),allocatable,save :: fluxsurf_lw     ! Incident Long Wave (IR) surface flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxsurf_sw     ! Incident Short Wave (stellar) surface flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxsurfabs_sw  ! Absorbed Short Wave (stellar) flux by the surface (W.m-2).
!$OMP THREADPRIVATE(fluxsurf_lw,fluxsurf_sw,fluxsurfabs_sw)

      real,dimension(:),allocatable,save :: fluxtop_lw      ! Outgoing LW (IR) flux to space (W.m-2).
      real,dimension(:),allocatable,save :: fluxabs_sw      ! Absorbed SW (stellar) flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxtop_dn      ! Incoming SW (stellar) radiation at the top of the atmosphere (W.m-2).
      real,dimension(:),allocatable,save :: fluxdyn         ! Horizontal heat transport by dynamics (W.m-2).
!$OMP THREADPRIVATE(fluxtop_lw,fluxabs_sw,fluxtop_dn,fluxdyn)

      real,dimension(:,:),allocatable,save :: OLR_nu        ! Outgoing LW radiation in each band (Normalized to the band width (W/m2/cm-1)).
      real,dimension(:,:),allocatable,save :: OSR_nu        ! Outgoing SW radiation in each band (Normalized to the band width (W/m2/cm-1)).
      real,dimension(:,:),allocatable,save :: zdtlw         ! LW heating tendencies (K/s).
      real,dimension(:,:),allocatable,save :: zdtsw         ! SW heating tendencies (K/s).
      !real,dimension(:),allocatable,save :: sensibFlux      ! Turbulent flux given by the atmosphere to the surface (W.m-2).
!$OMP THREADPRIVATE(OLR_nu,OSR_nu,zdtlw,zdtsw)

      real,dimension(:,:,:),allocatable,save :: int_dtauv   ! VI optical thickness of layers within narrowbands for diags ().
      real,dimension(:,:,:),allocatable,save :: int_dtaui   ! IR optical thickness of layers within narrowbands for diags ().
!$OMP THREADPRIVATE(int_dtaui,int_dtauv) 

      real,allocatable,dimension(:),save :: tau_col ! Total Aerosol Optical Depth.
      real,allocatable,save :: hice(:) ! Oceanic Ice height. by BC
!$OMP THREADPRIVATE(tau_col,hice)

      real,allocatable,dimension(:,:),save :: cloudfrac  ! Fraction of clouds (%).
      real,allocatable,dimension(:),save :: totcloudfrac ! Column fraction of clouds (%).
!$OMP THREADPRIVATE(cloudfrac,totcloudfrac)

      real,allocatable,dimension(:,:),save :: qsurf_hist
      real,allocatable,dimension(:,:,:),save :: nueffrad ! Aerosol effective radius variance. By RW
!$OMP THREADPRIVATE(qsurf_hist,nueffrad)

      real,allocatable,dimension(:),save :: ice_initial
      real,allocatable,dimension(:),save :: ice_min
!$OMP THREADPRIVATE(ice_initial,ice_min)

      real, dimension(:),allocatable,save ::  pctsrf_sic
      real, dimension(:,:),allocatable,save :: tslab
      real, dimension(:),allocatable,save ::  tsea_ice
      real, dimension(:),allocatable,save :: sea_ice
      real, dimension(:),allocatable,save :: zmasq
      integer, dimension(:),allocatable,save ::knindex
      real,allocatable,dimension(:,:,:),save :: reffrad
!$OMP THREADPRIVATE(pctsrf_sic,tslab,tsea_ice,sea_ice,zmasq,knindex,reffrad)
      
CONTAINS

!======================================================================
SUBROUTINE phys_state_var_init(nqtot)

IMPLICIT NONE

        integer,intent(in) :: nqtot

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
        print*,'klon',klon,'klev',klev
        ALLOCATE(phisfi(klon))
        ALLOCATE(tsurf(klon))
        ALLOCATE(tsoil(klon,nsoilmx))
        ALLOCATE(albedo(klon,L_NSPECTV))
        ALLOCATE(albedo_equivalent(klon))
        ALLOCATE(albedo_snow_SPECTV(L_NSPECTV))
        ALLOCATE(albedo_co2_ice_SPECTV(L_NSPECTV))
        ALLOCATE(albedo_bareground(klon))
        ALLOCATE(rnat(klon))
        ALLOCATE(emis(klon))
        ALLOCATE(dtrad(klon,klev))
        ALLOCATE(fluxrad_sky(klon))
        ALLOCATE(fluxrad(klon))
        ALLOCATE(capcal(klon))
        ALLOCATE(fluxgrd(klon))
        ALLOCATE(qsurf(klon,nqtot))
        ALLOCATE(q2(klon,klev+1))
        ALLOCATE(ztprevious(klon,klev))
        ALLOCATE(zuprevious(klon,klev))
        ALLOCATE(cloudfrac(klon,klev))
        ALLOCATE(totcloudfrac(klon))
        ALLOCATE(hice(klon))
        ALLOCATE(qsurf_hist(klon,nqtot))
        ALLOCATE(reffrad(klon,klev,naerkind))
        ALLOCATE(nueffrad(klon,klev,naerkind))
        ALLOCATE(ice_initial(klon))
        ALLOCATE(ice_min(klon))
        ALLOCATE(fluxsurf_lw(klon))
        ALLOCATE(fluxsurf_sw(klon))
        ALLOCATE(fluxsurfabs_sw(klon))
        ALLOCATE(fluxtop_lw(klon))
        ALLOCATE(fluxabs_sw(klon))
        ALLOCATE(fluxtop_dn(klon))
        ALLOCATE(fluxdyn(klon))
        ALLOCATE(OLR_nu(klon,L_NSPECTI))
        ALLOCATE(OSR_nu(klon,L_NSPECTV))
        ALLOCATE(int_dtaui(klon,klev,L_NSPECTI))
        ALLOCATE(int_dtauv(klon,klev,L_NSPECTV))
        ALLOCATE(sensibFlux(klon))
        ALLOCATE(zdtlw(klon,klev))
        ALLOCATE(zdtsw(klon,klev))
        ALLOCATE(tau_col(klon))
        ALLOCATE(pctsrf_sic(klon))
        ALLOCATE(tslab(klon,noceanmx))
        ALLOCATE(tsea_ice(klon))
        ALLOCATE(sea_ice(klon))
        ALLOCATE(zmasq(klon))
        ALLOCATE(knindex(klon))
        ! This is defined in comsaison_h
        ALLOCATE(mu0(klon))
        ALLOCATE(fract(klon))
         ! This is defined in radcommon_h
        ALLOCATE(glat(klon))
        ALLOCATE(albedodat(klon))
        ALLOCATE(zmea(klon))
        ALLOCATE(zstd(klon))
        ALLOCATE(zsig(klon))
        ALLOCATE(zgam(klon))
        ALLOCATE(zthe(klon))
        !allocate(l0(klon))
        allocate(wstar(klon))
        allocate(ustar(klon))
        allocate(tstar(klon))
        allocate(hfmax_th(klon))
        allocate(zmax_th(klon))

END SUBROUTINE phys_state_var_init

!======================================================================
SUBROUTINE phys_state_var_end

IMPLICIT NONE

        DEALLOCATE(tsurf)
        DEALLOCATE(tsoil)
        DEALLOCATE(albedo)
        DEALLOCATE(albedo_equivalent)
        DEALLOCATE(albedo_snow_SPECTV)
        DEALLOCATE(albedo_co2_ice_SPECTV)
        DEALLOCATE(albedo_bareground)
        DEALLOCATE(rnat)
        DEALLOCATE(emis)
        DEALLOCATE(dtrad)
        DEALLOCATE(fluxrad_sky)
        DEALLOCATE(fluxrad)
        DEALLOCATE(capcal)

        DEALLOCATE(fluxgrd)
        DEALLOCATE(qsurf)
        DEALLOCATE(q2)
        DEALLOCATE(ztprevious)
        DEALLOCATE(zuprevious)
        DEALLOCATE(cloudfrac)
        DEALLOCATE(totcloudfrac)
        DEALLOCATE(hice)
        DEALLOCATE(qsurf_hist)
        DEALLOCATE(reffrad)
        DEALLOCATE(nueffrad)
        DEALLOCATE(ice_initial)
        DEALLOCATE(ice_min)
        DEALLOCATE(fluxsurf_lw)
        DEALLOCATE(fluxsurf_sw)
        DEALLOCATE(fluxsurfabs_sw)
        DEALLOCATE(fluxtop_lw)
        DEALLOCATE(fluxabs_sw)
        DEALLOCATE(fluxtop_dn)
        DEALLOCATE(fluxdyn)
        DEALLOCATE(OLR_nu)
        DEALLOCATE(OSR_nu)
        DEALLOCATE(int_dtaui)
        DEALLOCATE(int_dtauv)
        DEALLOCATE(sensibFlux)
        DEALLOCATE(zdtlw)
        DEALLOCATE(zdtsw)
        DEALLOCATE(tau_col)
        DEALLOCATE(pctsrf_sic)
        DEALLOCATE(tslab)
        DEALLOCATE(tsea_ice)
        DEALLOCATE(sea_ice)
        DEALLOCATE(zmasq)
        DEALLOCATE(knindex)
        DEALLOCATE(mu0)
        DEALLOCATE(fract)
        DEALLOCATE(glat)
        DEALLOCATE(phisfi)
        DEALLOCATE(albedodat)
        DEALLOCATE(zmea)
        DEALLOCATE(zstd)
        DEALLOCATE(zsig)
        DEALLOCATE(zgam)
        DEALLOCATE(zthe)
        !deallocate(l0)
        deallocate(wstar)
        deallocate(ustar)
        deallocate(tstar)
        deallocate(hfmax_th)
        deallocate(zmax_th)


END SUBROUTINE phys_state_var_end

      END MODULE phys_state_var_mod
