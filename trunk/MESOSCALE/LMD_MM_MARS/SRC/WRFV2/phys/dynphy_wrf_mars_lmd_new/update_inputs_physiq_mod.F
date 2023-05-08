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
      real function ls2sol(ls)

!c  Returns solar longitude, Ls (in deg.), from day number (in sol),
!c  where sol=0=Ls=0 at the northern hemisphere spring equinox

      implicit none

!c  Arguments:
      real ls

!c  Local:
      double precision xref,zx0,zteta,zz
!c	xref: mean anomaly, zteta: true anomaly, zx0: eccentric anomaly
      double precision year_day 
      double precision peri_day,timeperi,e_elips
      double precision pi,degrad 
      parameter (year_day=668.6d0) ! number of sols in a amartian year
!c      data peri_day /485.0/
      parameter (peri_day=485.35d0) ! date (in sols) of perihelion
!c  timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
      parameter (timeperi=1.90258341759902d0)
      parameter (e_elips=0.0934d0)  ! eccentricity of orbit
      parameter (pi=3.14159265358979d0)
      parameter (degrad=57.2957795130823d0)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            ls2sol = 0.0
         else
            ls2sol = year_day
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      ls2sol = zz*year_day + peri_day
      if (ls2sol.lt.0.0) ls2sol = ls2sol + year_day
      if (ls2sol.ge.year_day) ls2sol = ls2sol - year_day

      return
      end function ls2sol
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

  IF (JULYR .ne. 9999) THEN
    !
    ! specified
    !
    JH_cur_split = (GMT + elaps/3700.) !! universal time (0<JH_cur_split<1): JH_cur_split=0.5 at 12:00 UT 
    JH_cur_split = MODULO(JH_cur_split,24.)   !! the two arguments of MODULO must be of the same type 
    JH_cur_split = JH_cur_split / 24.
    JD_cur = (JULDAY - 1 + INT((3700*GMT+elaps)/88800))
    JD_cur = MODULO(int(JD_cur),669)
    MY = (JULYR-2000) + (88800.*(JULDAY - 1)+3700.*GMT+elaps)/59496000.
    MY = INT(MY)
  ELSE
    !
    ! idealized
    !
    JH_cur_split = lct_input - lon_input / 15. + elaps/3700.
    JH_cur_split = MODULO(JH_cur_split,24.)
    JH_cur_split = JH_cur_split / 24.
    JD_cur = floor(ls2sol(ls_input)) + INT((3700*(lct_input - lon_input / 15.) + elaps)/88800)
    JD_cur = MODULO(int(JD_cur),669)
    MY = 2024
    !day_ini = floor(ls2sol(ls_input)) !! JD_cur at firstcall is day_ini
  ENDIF
  print *,'** Mars ** TIME IS', JD_cur, JH_cur_split*24.

END SUBROUTINE update_inputs_physiq_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_tracers(nq,MARS_MODE)

  INTEGER, INTENT(IN) :: nq,MARS_MODE

  ALLOCATE(traceurs(nq))

  !!! name of tracers to mimic entries in tracer.def 
  !!! ----> IT IS IMPORTANT TO KEEP SAME ORDERING AS IN THE REGISTRY !!!!
  !!!                               SAME NAMING   AS IN THE LMD PHYSICS !!!!
  SELECT CASE (MARS_MODE)
    CASE(0,10)  
      traceurs(nq) = 'co2'
    CASE(1)    ! scalar:qh2o,qh2o_ice
      traceurs(1)  = 'h2o_vap'
      traceurs(2)  = 'h2o_ice'
    CASE(2)    ! scalar:qdust
      traceurs(1)  = 'dust01'
    CASE(3)    ! scalar:qdust,qdustn
      traceurs(1)  = 'dust_mass'
      traceurs(2)  = 'dust_number'
    CASE(11)   ! scalar:qh2o,qh2o_ice,qdust,qdustn
      traceurs(1)  = 'h2o_vap'
      traceurs(2)  = 'h2o_ice'
      traceurs(3)  = 'dust_mass'
      traceurs(4)  = 'dust_number'
    CASE(12)
      traceurs(1)  = 'h2o_vap'
      traceurs(2)  = 'h2o_ice'
      traceurs(3)  = 'dust_mass'
      traceurs(4)  = 'dust_number'
      traceurs(5)  = 'ccn_mass'
      traceurs(6)  = 'ccn_number'
    CASE(20) 
      traceurs(1) = 'qtrac1'
    CASE(32) ! scalar:qh2o,qh2o_ice,qdust,qdustn,qccn,qccnn,qco2,qco2_ice,qccn_co2,qccnn_co2
      traceurs(1)   = 'h2o_vap'
      traceurs(2)   = 'h2o_ice'
      traceurs(3)   = 'dust_mass'
      traceurs(4)   = 'dust_number'
      traceurs(5)   = 'ccn_mass'
      traceurs(6)   = 'ccn_number'
      traceurs(7)   = 'co2'
      traceurs(8)   = 'co2_ice'
      traceurs(9)   = 'ccnco2_mass'
      traceurs(10)  = 'ccnco2_number'
    CASE(42)
      traceurs(1)  = 'co2'
      traceurs(2)  = 'co'
      traceurs(3)  = 'o'
      traceurs(4)  = 'o1d'
      traceurs(5)  = 'o2'
      traceurs(6)  = 'o3'
      traceurs(7)  = 'h'
      traceurs(8)  = 'h2'
      traceurs(9)  = 'oh'
      traceurs(10)  = 'ho2'
      traceurs(11)  = 'h2o2'
      traceurs(12)  = 'ch4'
      traceurs(13)  = 'n2'
      traceurs(14)  = 'ar'
      traceurs(15)  = 'h2o_ice'
      traceurs(16)  = 'h2o_vap'
      traceurs(17)  = 'dust_mass'
      traceurs(18)  = 'dust_number'
  END SELECT

  !!---------------------!!
  !! OUTPUT FOR CHECKING !!
  !!---------------------!!
  print*,"check: traceurs",traceurs

END SUBROUTINE update_inputs_physiq_tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_constants

   USE module_model_constants
   use comcstfi_h, only: omeg,mugaz
   use planete_h,  only: year_day,periheli,aphelie, &
                          peri_day,obliquit,emin_turb, &
                          lmixmin
   use surfdat_h,  only: emissiv,albedice,iceradius, &
                          emisice,dtemisice, &
                          z0_default
   use comsoil_h, only: volcapa

   !! comcstfi_h
   omeg = womeg                
   mugaz = wmugaz  
   print*,"check: omeg,mugaz"
   print*,omeg,mugaz
   !! planet_h
   year_day = wyear_day
   periheli = wperiheli
   aphelie = waphelie
   peri_day = wperi_day
   obliquit = wobliquit
   emin_turb = wemin_turb
   lmixmin = wlmixmin
   print*,"check: year_day,periheli,aphelie,peri_day,obliquit"
   print*,year_day,periheli,aphelie,peri_day,obliquit
   print*,"check: emin_turb,lmixmin"
   print*,emin_turb,lmixmin
   !! surfdat_h
   emissiv = wemissiv
   emisice(1) = wemissiceN
   emisice(2) = wemissiceS
   albedice(1) = walbediceN
   albedice(2) = walbediceS
   iceradius(1) = wiceradiusN
   iceradius(2) = wiceradiusS
   dtemisice(1) = wdtemisiceN
   dtemisice(2) = wdtemisiceS
   z0_default = wz0
   print*,"check: z0def,emissiv,emisice,albedice,iceradius,dtemisice"
   print*,z0_default,emissiv,emisice,albedice,iceradius,dtemisice
   !! comsoil_h
   volcapa = wvolcapa
   print*,"check: volcapa",volcapa

END SUBROUTINE update_inputs_physiq_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_geom( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,ngrid,nlayer,&
            DX,DY,MSFT,&
            lat_input, lon_input,&
            XLAT,XLONG)

   ! in WRF (share)
   USE module_model_constants, only: DEGRAD
   ! in LMD (phymars)
   use comgeomfi_h, only: ini_fillgeom
   ! in LMD (phy_common)
   USE mod_grid_phy_lmdz, ONLY: init_grid_phy_lmdz
   USE geometry_mod, ONLY: latitude,latitude_deg,&
                           longitude,longitude_deg,&
                           cell_area

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,ngrid,nlayer
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     MSFT,XLAT,XLONG
   REAL, INTENT(IN) :: dx,dy
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
    parea(subs) = (dx/msft(i,j))*(dy/msft(i,j))    !! 2. WRF map scale factors - assume that msfx=msfy (msf=covariance)
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
   call ini_fillgeom(ngrid,plat,plon,parea)

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

   use surfdat_h, only: phisfi, albedodat,  &
                        zmea, zstd, zsig, zgam, zthe, z0, &
                        tsurf, co2ice, emis, qsurf

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
     phisfi(subs) = phisfi_val

     !---------------!
     ! Ground albedo !
     !---------------!
     IF (JULYR .ne. 9999) THEN 
      IF (CST_AL == 0) THEN
       albedodat(subs)=M_ALBEDO(i,j)
      ELSE 
       albedodat(subs)=CST_AL
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT ALBEDO ', albedodat
      ENDIF
     ELSE
      IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION albedo: ', CST_AL
      albedodat(subs)=CST_AL
     ENDIF

     !-----------------------------------------!
     ! Gravity wave parametrization            !
     ! NB: usually 0 in mesoscale applications !
     !-----------------------------------------!
     IF (JULYR .ne. 9999) THEN
      zmea(subs)=M_GW(i,1,j)
      zstd(subs)=M_GW(i,2,j)
      zsig(subs)=M_GW(i,3,j)
      zgam(subs)=M_GW(i,4,j)
      zthe(subs)=M_GW(i,5,j)
     ELSE
      IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION GWdrag OFF'
      zmea(subs)=0.
      zstd(subs)=0.
      zsig(subs)=0.
      zgam(subs)=0.
      zthe(subs)=0.
     ENDIF

     !----------------------------!
     ! Variable surface roughness !
     !----------------------------!
     IF (JULYR .ne. 9999) THEN
      IF (CST_Z0 == 0) THEN
       z0(subs) = M_Z0(i,j)
      ELSE
       z0(subs) = CST_Z0
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT SURF ROUGHNESS (m) ',CST_Z0
      ENDIF
     ELSE
      IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION z0 (m) ', CST_Z0
      z0(subs)=CST_Z0
     ENDIF
     !!!! ADDITIONAL SECURITY. THIS MIGHT HAPPEN WITH OLD INIT FILES.
     IF (z0(subs) == 0.) THEN
       IF ( (i == ips) .AND. (j == jps) ) PRINT *, 'WELL, z0 is 0, this is no good. Setting to old defaults value 0.01 m'
       z0(subs) = 0.01
     ENDIF
     !!!! ADDITIONAL SECURITY. INTERP+SMOOTH IN GEOGRID MIGHT YIELD NEGATIVE Z0 !!!
     IF (z0(subs) < 0.) THEN
       PRINT *, 'WELL, z0 is NEGATIVE, this is impossible. better stop here.'
       PRINT *, 'advice --> correct interpolation / smoothing of z0 in WPS'
       PRINT *, '           -- or check the constant value set in namelist.input'
       STOP
     ENDIF

     !-----------------------------------------------!
     ! Ground temperature, emissivity, CO2 ice cover !
     !-----------------------------------------------!
     tsurf(subs) = M_TSURF(i,j)
     emis(subs) = M_EMISS(i,j)
     co2ice(subs) = M_CO2ICE(i,j)

     !-------------------!
     ! Tracer at surface !
     !-------------------!
     qsurf(subs,:)=0. ! default case
     SELECT CASE (MARS_MODE)
       CASE(1)
       qsurf(subs,2)=M_H2OICE(i,j)    !! logique avec noms(2) = 'h2o_ice' defini ci-dessus
                                      !! ----- retrocompatible ancienne physique
                                      !! ----- [H2O ice is last tracer in qsurf in LMD physics]
       CASE(2)   
       qsurf(subs,1)=0.                !! not coupled with lifting for the moment [non remobilise]
       !CASE(3)
       !qsurf(subs,1)=q_prof(1,1)                !!! temporaire, a definir       
       !qsurf(subs,2)=q_prof(1,2)     
       CASE(11)
       qsurf(subs,2)=M_H2OICE(i,j)    !! logique avec noms(2) = 'h2o_ice' defini ci-dessus
       qsurf(subs,3)=0.                !! not coupled with lifting for the moment [non remobilise]
       CASE(12)
       qsurf(subs,2)=M_H2OICE(i,j)    !! logique avec noms(2) = 'h2o_ice' defini ci-dessus
       qsurf(subs,3)=0.                !! not coupled with lifting for the moment [non remobilise]
     END SELECT

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: phisfi",phisfi(1),phisfi(nlast)
   print*,"check: albedodat",albedodat(1),albedodat(nlast)
   print*,"check: zmea",zmea(1),zmea(nlast)
   print*,"check: zstd",zstd(1),zstd(nlast)
   print*,"check: zsig",zsig(1),zsig(nlast)
   print*,"check: zgam",zgam(1),zgam(nlast)
   print*,"check: zthe",zthe(1),zthe(nlast)
   print*,"check: z0",z0(1),z0(nlast)
   print*,"check: tsurf",tsurf(1),tsurf(nlast)
   print*,"check: emis",emis(1),emis(nlast)
   print*,"check: co2ice",co2ice(1),co2ice(nlast)
   print*,"check: qsurf",qsurf(1,:),qsurf(nlast,:)

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

   use comsoil_h, only: inertiedat,mlayer,layer,tsoil

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,nsoil
   INTEGER :: i,j,subs,nlast,k
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

     !-----------------!
     ! Thermal Inertia !
     !-----------------!
     IF (JULYR .ne. 9999) THEN
      IF (CST_TI == 0) THEN
       inertiedat_val=M_TI(i,j)
      ELSE
       inertiedat_val=CST_TI
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT THERMAL INERTIA ', inertiedat_val
      ENDIF
     ELSE
      IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION inertia: ',CST_TI
      inertiedat_val=CST_TI
     ENDIF
     !inertiedat(subs) = inertiedat_val
     !--pb de dimensions???!!???
     IF (JULYR .ne. 9999) THEN
       inertiedat(subs,:)=M_ISOIL(i,:,j) !! verifier que cest bien hires TI en surface
       mlayer(0:nsoil-1)=M_DSOIL(i,:,j)
     ELSE
        IF ( nsoil .lt. 18 ) THEN
            PRINT *,'** Mars ** WRONG NUMBER OF SOIL LAYERS. SHOULD BE 18 AND IT IS ',nsoil
            STOP
        ENDIF
        IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION isoil and dsoil standard'
        do k=1,nsoil
         inertiedat(subs,k) = inertiedat_val
         !mlayer(k-1) = sqrt(887.75/3.14)*((2.**(k-0.5))-1.) * inertiedat_val / wvolcapa    ! old setting
         mlayer(k-1) = 2.E-4 * (2.**(k-0.5-1.))                                            ! new gcm settings
        enddo
     ENDIF
     IF ( (i == ips) .AND. (j == jps) ) &
         PRINT *,'** Mars ** TI and depth profiles are',inertiedat(subs,:),mlayer(0:nsoil-1)

     !!!!!!!!!!!!!!!!! DONE in soil_setting.F 
     ! 1.5 Build layer(); following the same law as mlayer()
     ! Assuming layer distribution follows mid-layer law:
     ! layer(k)=lay1*alpha**(k-1)
     lay1=sqrt(mlayer(0)*mlayer(1))
     alpha=mlayer(1)/mlayer(0)
     do k=1,nsoil
       layer(k)=lay1*(alpha**(k-1))
     enddo

     !------------------------!
     ! Deep soil temperatures ! 
     !------------------------!
     IF (M_TSOIL(i,1,j) .gt. 0. .and. JULYR .ne. 9999) THEN
       tsoil(subs,:)=M_TSOIL(i,:,j)
     ELSE
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** no tsoil. set it to tsurf.'
       do k=1,nsoil
        tsoil(subs,k) = M_TSURF(i,j)
       enddo
     ENDIF

   ENDDO
   ENDDO

   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: inertiedat",inertiedat(1,:),inertiedat(nlast,:)
   print*,"check: mlayer",mlayer(:)
   print*,"check: layer",layer(:)
   print*,"check: tsoil",tsoil(1,:),tsoil(nlast,:)

END SUBROUTINE update_inputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_turb( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,&
            RESTART,isles,&
            M_Q2,M_WSTAR)

   use turb_mod, only: q2,wstar,turb_resolved

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
   print*,"check: q2",q2(1,:),q2(nlast,:)
   print*,"check: wstar",wstar(1),wstar(nlast)

END SUBROUTINE update_inputs_physiq_turb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            RESTART,&
            M_FLUXRAD)

   use dimradmars_mod, only: fluxrad

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
      fluxrad(subs)=0.
     ELSE
      fluxrad(subs)=M_FLUXRAD(i,j)
     ENDIF
     !! et fluxrad_sky ???!???

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: fluxrad",fluxrad(1),fluxrad(nlast)

END SUBROUTINE update_inputs_physiq_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_slope( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,&
            SLPX,SLPY)

   USE module_model_constants, only: DEGRAD
   USE slope_mod, ONLY: theta_sl, psi_sl

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: SLPX,SLPY
   INTEGER :: i,j,subs,nlast

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !-------------------!
     ! Slope inclination !
     !-------------------!
     IF (JULYR .ne. 9999) THEN
       theta_sl(subs)=atan(sqrt( (1000.*SLPX(i,j))**2 + (1000.*SLPY(i,j))**2 ))
       theta_sl(subs)=theta_sl(subs)/DEGRAD
     ELSE
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'IDEALIZED SIMULATION: theta_sl = 0'
       theta_sl(subs)=0.
     ENDIF

     !-------------------------------------------!
     ! Slope orientation; 0 is north, 90 is east !
     !-------------------------------------------!
     IF (JULYR .ne. 9999) THEN
       psi_sl(subs)=-90.*DEGRAD-atan(SLPY(i,j)/SLPX(i,j))
       if (SLPX(i,j) .ge. 0.) then
         psi_sl(subs)=psi_sl(subs)-180.*DEGRAD
       endif
       psi_sl(subs)=360.*DEGRAD+psi_sl(subs)
       psi_sl(subs)=psi_sl(subs)/DEGRAD
       psi_sl(subs) = MODULO(psi_sl(subs)+180.,360.)
     ELSE
       IF ( (i == ips) .AND. (j == jps) ) PRINT *,'IDEALIZED SIMULATION: psi_sl = 0'
       psi_sl(subs) = 0.
     ENDIF

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: theta_sl",theta_sl(1),theta_sl(nlast)
   print*,"check: psi_sl",psi_sl(1),psi_sl(nlast)

END SUBROUTINE update_inputs_physiq_slope

END MODULE update_inputs_physiq_mod

