MODULE update_inputs_physiq_mod

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_time(&
            JULYR,JULDAY,GMT,&
            elaps,&
            lct_input,lon_input,ls_input,&
            MY)
  USE variables_mod, only: JD_cur,JH_cur_split,phour_ini
  use callkeys_mod, only : tlocked
  INTEGER, INTENT(IN) :: JULDAY, JULYR
  REAL, INTENT(IN) :: GMT,elaps,lon_input,ls_input,lct_input
  REAL,INTENT(OUT) :: MY
  REAL :: sec,nsec
  
  IF (JULYR .ne. 9999) THEN
    if (tlocked .eqv. .false.) THEN
      JH_cur_split = (GMT + elaps/3600.) !! universal time (0<JH_cur_split<1): JH_cur_split=0.5 at 12:00 UT 
      JH_cur_split = MODULO(JH_cur_split,24.)   !! the two arguments of MODULO must be of the same type 
      JH_cur_split = JH_cur_split / 24.
      JD_cur = (JULDAY - 1 + INT((3600.0*GMT+elaps)/86400))
      JD_cur = MODULO(int(JD_cur),2)
      MY = (JULYR-2000) + (86400*(JULDAY - 1)+3600.0*GMT+elaps)/31968000
      MY = INT(MY)
    ELSE
      JH_cur_split = (GMT)! + elaps/420000.) !! universal time (0<JH_cur_split<1): JH_cur_split=0.5 at 12:00 UT 
      JH_cur_split = MODULO(JH_cur_split,24.)   !! the two arguments of MODULO must be of the same type 
      JH_cur_split = JH_cur_split / 24.
      JD_cur = (JULDAY - 1 + INT((3600*GMT)/86400))!+elaps)/1.008e7))
      JD_cur = MODULO(int(JD_cur),2)
      MY = (JULYR-2000) + (86400*(JULDAY - 1)+3600*GMT+elaps)/31968000
      MY = INT(MY)
    ENDIF
  ELSE
    if (tlocked .eqv. .false.) THEN
      JH_cur_split = lct_input - lon_input / 15. + elaps/1500.0
      JD_cur =  INT((sec*(lct_input - lon_input / 15.) + elaps)/36000)
      !ptime = lct_input - lon_input / 15. + elaps/3600.
      !pday =  INT((3600*(lct_input - lon_input / 15.) + elaps)/86400)
    ELSE
      JH_cur_split = lct_input - lon_input / 15. !+ elaps/1500.0
      !pday =  INT((sec*(lct_input - lon_input / 15.)+ elaps)/36000)
      JD_cur =  INT((sec*(lct_input - lon_input / 15.))/3600)
      !print*,'ptime',ptime
      !print*,'pday',pday 
    !pday =  INT((3600*(lct_input - lon_input / 15.) + elaps)/86400)
      JH_cur_split = MODULO(ptime,24.)
      JH_cur_split = JH_cur_split / 24.
      JD_cur = MODULO(int(pday),365)
      MY = 2024
    ENDIF
  ENDIF

END SUBROUTINE update_inputs_physiq_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_tracers(nq,MARS_MODE)

  use tracer_h, only: noms,nqtot
  !use phys_state_var_mod , only : qsurf
  !USE dimphy, only : klon,klev 
  INTEGER, INTENT(IN) :: nq,MARS_MODE
  INTEGER :: i,k
!  CHARACTER(len=20), INTENT(INOUT) :: tname(nq) ! tracer names
  logical :: end_of_file

  !! TRACERS
  IF (.not.ALLOCATED(noms)) ALLOCATE(noms(nq)) !! est fait dans initracer normalement
      
                 !! tableau dans tracer_mod.F90
  print*,'nq',nq
  print*,'MARS_MODE',MARS_MODE
  nqtot=nq
  if (MARS_MODE .eq. 1) THEN
    noms(1)="h2o_vap"
    noms(2)="h2o_ice"
  else if (MARS_MODE .eq. 22)  then 
    noms(1)="co2_vap"
    noms(2)="co2_ice"
  else
    noms(:)="zolbxs"
  endif
  !noms(:)=tname(:)
  print*,'noms',noms
  !!---------------------!!
  !! OUTPUT FOR CHECKING !!
  !!---------------------!!

END SUBROUTINE update_inputs_physiq_tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_constants

   !USE module_model_constants
   !use comcstfi_h, only: omeg,mugaz
   !use planete_h,  only: year_day,periheli,aphelie, &
   !                       peri_day,obliquit,emin_turb, &
   !                       lmixmin
   use planete_mod, only: year_day, periastr, apoastr, peri_day,&
                       obliquit, z0, lmixmin, emin_turb
   use surfdat_h,  only: emissiv,iceradius, &
                         emisice,dtemisice
   !                       z0_default
   !use comsoil_h, only: volcapa
   use comcstfi_mod, only: rad, omeg, g, mugaz, rcp, cpp, r
   !! comcstfi_h
   use phys_state_var_mod, only :cloudfrac,totcloudfrac,hice,rnat,pctsrf_sic,tsea_ice
   !use iniorbit
   !use time_phylmdz_mod, only: dtphys, daysec,day_ini
  
  
     !open(17,file='controle.txt',form='formatted',status='old')
     !rewind(17)
     !read(17,*)
     !read(17,*) 
     !read(17,*) day_ini !(tab0+3)
     !read(17,*) 
     !read(17,*) !tab0+5)
     !read(17,*) omeg !(tab0+6)
     !read(17,*) !(tab0+7)
     !read(17,*)  !(tab0+8)
     !read(17,*)  !(tab0+9)
     !read(17,*) daysec
     !read(17,*) dtphys !tab0+11)
     !read(17,*)
     !read(17,*)
     !read(17,*) year_day !(tab0+14)
     !read(17,*) periastr !tab0+15)
     !read(17,*) apoastr !tab0+16)
     !read(17,*) peri_day !tab0+17)
     !read(17,*) obliquit !tab0+18)
     !read(17,*) z0
     !read(17,*)
     !read(17,*)
     !read(17,*)
     !read(17,*)
     !read(17,*) emisice(1)
     !read(17,*) emisice(2)
     !read(17,*) emissiv
     !read(17,*)
     !read(17,*)
     !read(17,*)
     !read(17,*)
     !read(17,*) iceradius(1)
     !read(17,*) iceradius(2)
     !read(17,*) dtemisice(1)
     !read(17,*) dtemisice(2)
     !close(17) 
     !cpp=(8.314511/(mugaz/1000.0))/rcp
     !print*,'cpp',cpp
     !print*,'g',g

     !emissiv(:)=EMIS
     !cloudfrac(:,:)=0.5
     !totcloudfrac(:)=0.5
     !hice(:)=0.
     !rnat(:)=0.
     !pctsrf_sic(:)=0.
     !tsea_ice(:)=0.
     !qsurf(:,:) = 0. 
     !print*,'iceradius',iceradius,'dtemisice',dtemisice
     !print*,'apoastr,periastr,year_day,peri_day,obliq',apoastr,periastr,year_day,peri_day,obliquit
     !print*,'emissiv',emissiv
 
     ! SUBROUTINE iniorbit(apoastr,periastr,year_day,peri_day,obliq)

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
   USE module_model_constants, only: DEGRAD,p0
   ! in LMD (phymars)
   !use comgeomfi_h, only: ini_fillgeom
   ! in LMD (phy_common)
   USE mod_grid_phy_lmdz, ONLY: init_grid_phy_lmdz
   USE geometry_mod, ONLY: latitude,latitude_deg,&
                           longitude,longitude_deg,&
                           cell_area
   use comdiurn_h, only: sinlat, coslat, sinlon, coslon
   !use planetwide_mod, only: planetwide_sumval
   use comgeomfi_h, only: totarea, totarea_planet
   use planete_mod, only: ini_planete_mod

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,ngrid,nlayer
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     MSFT,XLAT,XLONG
   REAL, INTENT(IN) :: dx,dy
   REAL, INTENT(IN) :: lat_input, lon_input
   INTEGER :: i,j,subs,ig,k
   REAL,DIMENSION(ngrid) :: plon, plat, parea
   REAL, DIMENSION(nlayer+1) :: znw
   REAL*8, DIMENSION(nlayer+1) :: apdyn,bpdyn
   REAL :: ptop
   REAL*8 :: PP0
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
   longitude(:) = plon(:)
   latitude(:) = plat(:)
   longitude_deg(:) = plon(:)/DEGRAD
   latitude_deg(:) = plat(:)/DEGRAD
   cell_area(:) = parea(:)
   !call planetwide_sumval(parea,totarea_planet)
   !print*,'parea',parea(1)
   !totarea=SSUM(ngrid,parea,1)
   totarea=ngrid*parea(1)
   !totarea_planet=SSUM(ngrid,parea,1)
   totarea_planet=ngrid*parea(1)

   IF (.not.ALLOCATED(sinlat)) ALLOCATE(sinlat(ngrid))
   IF (.not.ALLOCATED(coslat)) ALLOCATE(coslat(ngrid))
   IF (.not.ALLOCATED(sinlon)) ALLOCATE(sinlon(ngrid))
   IF (.not.ALLOCATED(coslon)) ALLOCATE(coslon(ngrid))
   DO ig=1,ngrid
     sinlat(ig)=sin(plat(ig))
     coslat(ig)=cos(plat(ig))
     sinlon(ig)=sin(plon(ig))
     coslon(ig)=cos(plon(ig))
   ENDDO
   
   open(unit=12,file='levels',form='formatted',status='old')
   rewind(12)
   DO k=1, nlayer
   read(12,*) znw(k)
   !write(6,*) 'read level ', k,grid%znw(k)
   ENDDO
   close(12)
   
   ptop=0.5
   apdyn=ptop*(1-znw)
   bpdyn=znw
   PP0=p0
   CALL ini_planete_mod(nlayer,PP0,apdyn,bpdyn)
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
                       zmea, zstd, zsig, zgam, zthe
   use slab_ice_h, only: noceanmx
   use planete_mod, only: z0
   use phys_state_var_mod, only : tsurf, emis, qsurf,tslab

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,MARS_MODE
   INTEGER :: i,j,subs,nlast,iq
   REAL, INTENT(IN  ) :: CST_AL, phisfi_val, CST_Z0
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     M_ALBEDO,M_TSURF,M_EMISS,M_CO2ICE,M_H2OICE,M_Z0
   REAL, DIMENSION( ims:ime, 5, jms:jme ), INTENT(IN   )  :: M_GW  

   !print*,'ALLOCATED(phisfi)',ALLOCATED(phisfi)
   !print*,'size phisfi',size(phisfi)
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
!     print*,'size phisfi',size(phisfi)
     !print*,'phisfi',phisfi(subs)
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
     zmea(subs)=0.
     zstd(subs)=0.
     zsig(subs)=0.
     zgam(subs)=0.
     zthe(subs)=0.

     !----------------------------!
     ! Variable surface roughness !
     !----------------------------!
     z0 = CST_Z0
     !IF (JULYR .ne. 9999) THEN
     ! IF (CST_Z0 == 0) THEN
     !  z0(subs) = M_Z0(i,j)
     ! ELSE
     !  z0(subs) = CST_Z0
     !  IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT SURF ROUGHNESS (m) ',CST_Z0
     ! ENDIF
     !ELSE
     ! IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION z0 (m) ', CST_Z0
     ! z0(subs)=CST_Z0
     !ENDIF
     !!!!! ADDITIONAL SECURITY. THIS MIGHT HAPPEN WITH OLD INIT FILES.
     !IF (z0(subs) == 0.) THEN
     !  IF ( (i == ips) .AND. (j == jps) ) PRINT *, 'WELL, z0 is 0, this is no good. Setting to old defaults value 0.01 m'
     !  z0(subs) = 0.01
     !ENDIF
     !!!!! ADDITIONAL SECURITY. INTERP+SMOOTH IN GEOGRID MIGHT YIELD NEGATIVE Z0 !!!
     !IF (z0(subs) < 0.) THEN
     !  PRINT *, 'WELL, z0 is NEGATIVE, this is impossible. better stop here.'
     !  PRINT *, 'advice --> correct interpolation / smoothing of z0 in WPS'
     !  PRINT *, '           -- or check the constant value set in namelist.input'
     !  STOP
     !ENDIF

     !-----------------------------------------------!
     ! Ground temperature, emissivity, CO2 ice cover !
     !-----------------------------------------------!
     tsurf(subs) = M_TSURF(i,j)
     emis(subs) = M_EMISS(i,j)     
     !do i=1,noceanmx
     tslab(subs,:)=tsurf(subs)
     !enddo
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
   print*,"check: z0",z0
   print*,"check: tsurf",tsurf(1),tsurf(nlast)
   print*,"check: emis",emis(1),emis(nlast)
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

   use comsoil_h, only: inertiedat,mlayer,layer,volcapa
   use phys_state_var_mod, only: tsoil

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,nsoil
   INTEGER :: i,j,subs,nlast
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
         PRINT *,'** Mars ** TI and depth profiles are',inertiedat(subs,:)!,mlayer(0:nsoil-1)

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
        !print*,'M_TSURF(i,j)',M_TSURF(1,:)
        !print*,'size M_TSURF',size(M_TSURF)
        !print*,'size tsoil',size(tsoil)
        tsoil(subs,k) = M_TSURF(i,j)
        !print*,'tsoil(subs,k)',tsoil(subs,k)
       enddo
     ENDIF

   ENDDO
   ENDDO

   volcapa=1.e6   

   print*,'zolbxs'
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
   !use phys_state_var_mod, only : q2,

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER :: i,j,subs,nlast      
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: M_WSTAR
   REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(IN) :: M_Q2
   LOGICAL, INTENT(IN ) :: RESTART,isles

   !! to know if this is turbulence-resolving run or not

   turb_resolved = isles
   print*,'isles',isles

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
      q2(subs,:)=M_Q2(i,:,j)!
      !q2(subs,:) = 1.e-3
      wstar(subs)=M_WSTAR(i,j)
     ENDIF

   ENDDO
   ENDDO
  
   !!---------------------!!
   !! OUTPUT FOR CHECKING !!
   !!---------------------!!
   nlast = (ipe-ips+1)*(jpe-jps+1)
   print*,"check: q2",q2(1,1)!,q2(nlast,:)
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
   use phys_state_var_mod, only : fluxrad

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

   !USE module_model_constants, only: DEGRAD
   !USE slope_mod, ONLY: theta_sl, psi_sl

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: SLPX,SLPY
   INTEGER :: i,j,subs,nlast

END SUBROUTINE update_inputs_physiq_slope

END MODULE update_inputs_physiq_mod
