MODULE update_inputs_physiq_mod

IMPLICIT NONE

CHARACTER(len=20),save,allocatable,dimension(:) :: traceurs ! tracer names

CONTAINS

!SUBROUTINE update_inputs_physiq_time
!SUBROUTINE update_inputs_physiq_tracers
!SUBROUTINE update_inputs_physiq_constants
!SUBROUTINE update_inputs_physiq_geom
!SUBROUTINE update_inputs_physiq_surf
!SUBROUTINE update_inputs_physiq_soil
!SUBROUTINE update_inputs_physiq_turb
!SUBROUTINE update_inputs_physiq_rad
!SUBROUTINE update_inputs_physiq_slope

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!c  Returns solar longitude, Ls (in deg.), from day number (in sol),
!c  where sol=0=Ls=0 at the northern hemisphere spring equinox


!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_time(&
            JULYR,JULDAY,GMT,&
            elaps,&
            lct_input,lon_input,ls_input,&
            MY)

  USE variables_mod, only: JD_cur,JH_cur_split,phour_ini
  !! JD_cur <> pday ! Julian day
  !! JH_cur_split <> ptime ! Julian hour (fraction of day)

  implicit none

  INTEGER, INTENT(IN) :: JULDAY, JULYR
  REAL, INTENT(IN) :: GMT,elaps,lon_input,ls_input,lct_input
  REAL,INTENT(OUT) :: MY

    !
    ! specified
    !
  IF (JULYR .ne. 9999) THEN
    JH_cur_split = (GMT + elaps/421200.) !! universal time (0<JH_cur_split<1): JH_cur_split=0.5 at 12:00 UT 
    JH_cur_split = MODULO(JH_cur_split,24.)   !! the two arguments of MODULO must be of the same type 
    JH_cur_split = JH_cur_split / 24.
    JD_cur = (JULDAY - 1 + INT((421200.*GMT+elaps)/10108800.))
    JD_cur = MODULO(int(JD_cur),2)
    MY = (JULYR-2000) + (10108800.*(JULDAY - 1)+421200.*GMT+elaps)/20217600.
    MY = INT(MY)
  ELSE
    JH_cur_split = lct_input - lon_input / 15. + elaps/(421200.)
    JH_cur_split = MODULO(JH_cur_split,24.)
    JH_cur_split = JH_cur_split / 24.
    print*,'JH_cur_split',JH_cur_split
    JD_cur = INT((3700*(lct_input - lon_input / 15.) +elaps)/88800)
    JD_cur = MODULO(int(JD_cur),669)
    MY = 2024
 ENDIF

END SUBROUTINE update_inputs_physiq_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_tracers(nq,MARS_MODE)

  use infotrac_phy, only: nqtot,tname
  use logic_mod, only: iflag_trac
  implicit none
  INTEGER, INTENT(IN) :: nq,MARS_MODE

  IF (.not.ALLOCATED(tname)) ALLOCATE(tname(nq))
  nqtot=nq

  IF (MARS_MODE .EQ. 34) THEN
    tname(1)="co2"
    tname(2)="co"
    tname(3)="h2"
    tname(4)="h2o"
    tname(5)="o1d"
    tname(6)="o"
    tname(7)="o2"
    tname(8)="o2dg"
    tname(9)="o3"
    tname(10)="h"
    tname(11)="oh"
    tname(12)="ho2"
    tname(13)="h2o2"
    tname(14)="cl"
    tname(15)="clo"
    tname(16)="cl2"
    tname(17)="hcl"
    tname(18)="hocl"
    tname(19)="clco"
    tname(20)="clco3"
    tname(21)="cocl2"
    tname(22)="s"
    tname(23)="so"
    tname(24)="so2"
    tname(25)="so3"
    tname(26)="s2o2"
    tname(27)="ocs"
    tname(28)="hso3"
    tname(29)="h2so4"
    tname(30)="s2"
    tname(31)="clso2"
    tname(32)="oscl"
    tname(33)="h2oliq"
    tname(34)="h2so4liq"
    iflag_trac=1
  ELSE
    tname(:)='zolbxs'
  ENDIF

END SUBROUTINE update_inputs_physiq_tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_constants

   !USE module_model_constants
   !use comcstfi_h, only: omeg,mugaz
   !use planete_h,  only: year_day,periheli,aphelie, &
   !                       peri_day,obliquit,emin_turb, &
   !                       lmixmin
   !use surfdat_h,  only: emissiv,albedice,iceradius, &
   !                      emisice,dtemisice, &
   !                       z0_default
   !use comsoil_h, only: volcapa

   implicit none

   !! comcstfi_h
   !omeg = womeg                
   !mugaz = wmugaz  
   !print*,"check: omeg,mugaz"
   !print*,omeg,mugaz
   !! planet_h

END SUBROUTINE update_inputs_physiq_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_geom( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,ngrid,nlayer,&
            DDX,DDY,MSFT,&
            lat_input, lon_input,&
            XLAT,XLONG)

   ! in WRF (share)
   USE module_model_constants, only: DEGRAD
   ! in LMD (phymars)
   !use comgeomfi_h, only: ini_fillgeom
   ! in LMD (phy_common)
   USE mod_grid_phy_lmdz, ONLY: init_grid_phy_lmdz
   USE geometry_mod, ONLY: latitude,latitude_deg,&
                           longitude,longitude_deg,&
                           cell_area,dx,dy

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,ngrid,nlayer
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     MSFT,XLAT,XLONG
   REAL, INTENT(IN) :: ddx,ddy
   REAL, INTENT(IN) :: lat_input, lon_input
   INTEGER :: i,j,subs
   REAL,DIMENSION(ngrid) :: plon, plat, parea

   DO j = jps,jpe
   DO i = ips,ipe

    !-----------------------------------!
    ! 1D subscript for physics "cursor" !
    !-----------------------------------!
    subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

    !----------------------------------------!
    ! Surface of each part of the grid (m^2) !  
    !----------------------------------------!
    !parea(subs) = dx*dy                           !! 1. idealized cases - computational grid
    parea(subs) = (ddx/msft(i,j))*(ddy/msft(i,j))    !! 2. WRF map scale factors - assume that msfx=msfy (msf=covariance)
    !parea(subs)=dx*dy/msfu(i,j)                   !! 3. special for Mercator GCM-like simulations

    !---------------------------------------------!
    ! Mass-point latitude and longitude (radians) !
    !---------------------------------------------!
    IF (JULYR .ne. 9999) THEN
     plat(subs) = XLAT(i,j)*DEGRAD
     plon(subs) = XLONG(i,j)*DEGRAD
    ELSE
     !!! IDEALIZED CASE
     IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION lat: ',lat_input
     IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION lon: ',lon_input
     plat(subs) = lat_input*DEGRAD
     plon(subs) = lon_input*DEGRAD
    ENDIF 

   ENDDO
   ENDDO

   !! FILL GEOMETRICAL ARRAYS !!
   !call ini_fillgeom(ngrid,plat,plon,parea)

   !!! ----------------------------------------------------------
   !!! --- initializing geometry in phy_common
   !!! --- (this is quite planet-independent)
   !!! ----------------------------------------------------------
   ! initialize mod_grid_phy_lmdz
   CALL init_grid_phy_lmdz(1,1,ipe-ips+1,jpe-jps+1,nlayer)
   ! fill in geometry_mod variables
   ! ... copy over local grid longitudes and latitudes
   ! ... partly what is done in init_geometry
   IF(.not.ALLOCATED(longitude)) ALLOCATE(longitude(ngrid))
   IF(.not.ALLOCATED(longitude_deg)) ALLOCATE(longitude_deg(ngrid))
   IF(.not.ALLOCATED(latitude)) ALLOCATE(latitude(ngrid))
   IF(.not.ALLOCATED(latitude_deg)) ALLOCATE(latitude_deg(ngrid))
   IF(.not.ALLOCATED(cell_area)) ALLOCATE(cell_area(ngrid))
   IF (.not.ALLOCATED(dx)) ALLOCATE(dx(ngrid))
   IF (.not.ALLOCATED(dy)) ALLOCATE(dy(ngrid))
   longitude(:) = plon(:)
   latitude(:) = plat(:)
   longitude_deg(:) = plon(:)/DEGRAD
   latitude_deg(:) = plat(:)/DEGRAD
   cell_area(:) = parea(:)
   !!! ----------------------------------------------------------

END SUBROUTINE update_inputs_physiq_geom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_surf( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,MARS_MODE,&
            M_ALBEDO,CST_AL,&
            M_TSURF,M_EMISS,M_CO2ICE,&
            M_GW,M_Z0,CST_Z0,&
            M_H2OICE,&
            phisfi_val)

   use phys_state_var_mod, only : falbe,ftsol 

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,MARS_MODE
   INTEGER :: i,j,subs,nlast
   REAL, INTENT(IN  ) :: CST_AL, phisfi_val, CST_Z0
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     M_ALBEDO,M_TSURF,M_EMISS,M_CO2ICE,M_H2OICE,M_Z0
   REAL, DIMENSION( ims:ime, 5, jms:jme ), INTENT(IN   )  :: M_GW  

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !---------------------!
     ! Ground geopotential !
     !---------------------!
     !phisfi(subs) = phisfi_val

     !---------------!
     ! Ground albedo !
     !---------------!
     IF (JULYR .ne. 9999) THEN 
      IF (CST_AL == 0) THEN
       falbe(subs)=M_ALBEDO(i,j)
      ELSE 
       falbe(subs)=CST_AL
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT ALBEDO ', falbe
      ENDIF
     ELSE
      IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION albedo: ', CST_AL
      falbe(subs)=CST_AL
     ENDIF

     !-----------------------------------------!
     ! Gravity wave parametrization            !
     ! NB: usually 0 in mesoscale applications !
     !-----------------------------------------!
      !zmea(subs)=0.
      !zstd(subs)=0.
      !zsig(subs)=0.
      !zgam(subs)=0.
      !zthe(subs)=0.


     !-----------------------------------------------!
     ! Ground temperature, emissivity, CO2 ice cover !
     !-----------------------------------------------!
     ftsol(subs) = M_TSURF(i,j)
     !emis(subs) = M_EMISS(i,j)

     !-------------------!
     ! Tracer at surface !
     !-------------------!
     !qsurf(subs,:)=0. ! default case


   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   !print*,"check: phisfi",phisfi(1),phisfi(nlast)
   print*,"check: albedodat",falbe(1),falbe(nlast)
   print*,"check: tsurf",ftsol(1),ftsol(nlast)
   !print*,"check: emis",emis(1),emis(nlast)


END SUBROUTINE update_inputs_physiq_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_soil( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,nsoil,&
            M_TI,CST_TI,&
            M_ISOIL,M_DSOIL,&
            M_TSOIL,M_TSURF)

   !use comsoil_h, only: tsoil
   use phys_state_var_mod, only : ftsoil

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,nsoil
   INTEGER :: i,j,k,subs,nlast
   REAL, INTENT(IN  ) :: CST_TI 
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     M_TI, M_TSURF
   REAL, DIMENSION( ims:ime, nsoil, jms:jme ), INTENT(IN)  :: &
     M_TSOIL, M_ISOIL, M_DSOIL
   REAL :: inertiedat_val
   REAL :: lay1,alpha

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !------------------------!
     ! Deep soil temperatures ! 
     !------------------------!
     IF (M_TSOIL(i,1,j) .gt. 0. .and. JULYR .ne. 9999) THEN
       ftsoil(subs,:)=M_TSOIL(i,:,j)
     ELSE
       IF ((i == ips) .AND. (j == jps)) PRINT *,'no tsoil.set it tsurf.'
       do k=1,nsoil
        ftsoil(subs,k) = M_TSURF(i,j)
       enddo
     ENDIF

   ENDDO
   ENDDO
   !print*,'M_TSURF',M_TSURF(i,j)
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   !print*,"check: inertiedat",inertiedat(1,:),inertiedat(nlast,:)
   !print*,"check: mlayer",mlayer(:)
   !print*,"check: layer",layer(:)
   print*,"check: tsoil",ftsoil(1,:),ftsoil(nlast,:)

END SUBROUTINE update_inputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_turb( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,&
            RESTART,isles,&
            M_Q2,M_WSTAR)

   use turb_mod, only: q2,wstar,turb_resolved

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER :: i,j,subs,nlast      
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: M_WSTAR
   REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(IN) :: M_Q2
   LOGICAL, INTENT(IN ) :: RESTART,isles

   !! to know if this is turbulence-resolving run or not
   turb_resolved = isles

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !PBL wind variance
     IF (.not. restart) THEN
      q2(subs,:) = 1.E-6      
      wstar(subs)=0.
     ELSE
      q2(subs,:)=M_Q2(i,:,j)
      wstar(subs)=M_WSTAR(i,j)
     ENDIF

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: q2",q2(1,1),q2(nlast,kme+1)
   print*,"check: wstar",wstar(1),wstar(nlast)

END SUBROUTINE update_inputs_physiq_turb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            RESTART,&
            M_FLUXRAD)

   !use dimradmars_mod, only: fluxrad
   use phys_state_var_mod, only : radsol

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER :: i,j,subs,nlast      
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: M_FLUXRAD
   LOGICAL, INTENT(IN ) :: RESTART

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     ! fluxrad_save
     IF (.not. restart) THEN
      radsol(subs)=0.
     ELSE
      radsol(subs)=M_FLUXRAD(i,j)
     ENDIF
     !! et fluxrad_sky ???!???

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: fluxrad",radsol(1),radsol(nlast)

END SUBROUTINE update_inputs_physiq_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_slope( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,&
            SLPX,SLPY)

   !USE module_model_constants, only: DEGRAD
   !USE slope_mod, ONLY: theta_sl, psi_sl

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: SLPX,SLPY
   INTEGER :: i,j,subs,nlast

   !! Nothing is done: no slope package (yet) in Venus' physics

END SUBROUTINE update_inputs_physiq_slope

END MODULE update_inputs_physiq_mod

