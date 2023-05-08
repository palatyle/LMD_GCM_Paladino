!WRF:MODEL_LAYER:CONSTANTS
!
!PLANET MARS

 MODULE module_model_constants
   USE module_configure
   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.
   character(len=15) :: planet
   REAL    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical parameter used within the model.

   !! added in WRFV3
   real    , parameter :: mwdry        = 43.49 ! molecular weight of dry air (g/mole)
!
!-------------------------------

   REAL    , PARAMETER :: asselin      = .025
!   REAL    , PARAMETER :: asselin      = .0
   REAL    , PARAMETER :: cb           = 25.

   REAL    , PARAMETER :: XLV0         = 3.15E6
   REAL    , PARAMETER :: XLV1         = 2370.
   REAL    , PARAMETER :: XLS0         = 2.905E6
   REAL    , PARAMETER :: XLS1         = 259.532

   REAL    , PARAMETER :: XLS          = 2.85E6
   REAL    , PARAMETER :: XLV          = 2.5E6
   REAL    , PARAMETER :: XLF          = 3.50E5

   REAL    , PARAMETER :: rhowater     = 1000.
   REAL    , PARAMETER :: rhosnow      = 100.
   REAL    , PARAMETER :: rhoair0      = 0.02   !Earth Surface density: 1.217 kg/m3
                                                !Mars Surface density: 0.020 kg/m3
   REAL    , PARAMETER :: DEGRAD       = 3.1415926/180.
   REAL    , PARAMETER :: DPD          = 360./669. ! longitude solaire d'un jour?

   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  KARMAN=0.4  ! von karman constant 
   REAL    , PARAMETER ::  STBOLT=5.67051E-8 ! stefan-boltzmann constant

   !! NB: c_s and c_k in WRFV2 moved to module_diffusion_em.F
   !! NB: no need in WRFV3 because set in namelist
   REAL    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
!   REAL    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL    , PARAMETER ::  w_alpha = 2.0 ! strength m/s/s
!   REAL    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number
   REAL    , PARAMETER ::  w_beta  = 0.8 ! activation cfl number

!! in WRFV3
!   REAL    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
!   REAL    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number
!!
       REAL , PARAMETER ::  pq0=379.90516
       REAL , PARAMETER ::  epsq2=0.2
!      REAL , PARAMETER ::  epsq2=0.02
       REAL , PARAMETER ::  a2=17.2693882
       REAL , PARAMETER ::  a3=273.16
       REAL , PARAMETER ::  a4=35.86
       REAL , PARAMETER ::  epsq=1.e-12
!#if ( NMM_CORE == 1 )
       REAL , PARAMETER ::  climit=1.e-20
       REAL , PARAMETER ::  cm1=2937.4
       REAL , PARAMETER ::  cm2=4.9283
       REAL , PARAMETER ::  cm3=23.5518
!       REAL , PARAMETER ::  defc=8.0
!       REAL , PARAMETER ::  defm=32.0
       REAL , PARAMETER ::  defc=0.0
       REAL , PARAMETER ::  defm=99999.0
       REAL , PARAMETER ::  epsfc=1./1.05
       REAL , PARAMETER ::  epswet=0.0
       REAL , PARAMETER ::  fcdif=1./3.
!       REAL , PARAMETER ::  fcm=0.003
       REAL , PARAMETER ::  fcm=0.00003
!       REAL , PARAMETER ::  fcm=0.0
       REAL , PARAMETER ::  p400=40000.0
       REAL , PARAMETER ::  phitp=15000.0
       REAL , PARAMETER ::  pi2=2.*3.1415926
       REAL , PARAMETER ::  plbtm=105000.0
       REAL , PARAMETER ::  plomd=64200.0
       REAL , PARAMETER ::  pmdhi=35000.0
       REAL , PARAMETER ::  q2ini=0.50
       REAL , PARAMETER ::  rhcrit_land=0.75
       REAL , PARAMETER ::  rhcrit_sea=0.80
       REAL , PARAMETER ::  rlag=14.8125
       REAL , PARAMETER ::  rlx=0.90
       REAL , PARAMETER ::  scq2=50.0
       REAL , PARAMETER ::  slopht=0.001
       REAL , PARAMETER ::  tlc=2.*0.703972477
       REAL , PARAMETER ::  wa=0.15
       REAL , PARAMETER ::  wght=0.35
       REAL , PARAMETER ::  wpc=0.075
       REAL , PARAMETER ::  z0land=0.10
!       REAL , PARAMETER ::  z0max=0.01
       REAL , PARAMETER ::  z0max=0.008
       REAL , PARAMETER ::  z0sea=0.001
       REAL, PARAMETER :: r_v=461.6 ! gas constant for water vapor
!#endif

!!!! needed in WRFV3
   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL , PARAMETER :: OBLIQUITY = 23.5
   REAL , PARAMETER :: ECCENTRICITY = 0.014
   REAL , PARAMETER :: SEMIMAJORAXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL , PARAMETER :: EQUINOX_FRACTION= 0.0
!!!! needed in WRFV3

   !  These are the physical constants planet dependent used within the model.

     real :: g
     REAL :: r_d        
     REAL :: cp         
     !REAL :: cpp
     REAL :: cv               
     REAL :: cpv         
     REAL :: cvv          
     REAL :: cvpm        
     REAL :: cliq        
     REAL :: cice         
     REAL :: psat         
     REAL :: rcv          
     REAL :: rcp       
     REAL :: rovg        
     REAL :: c2     
     REAL :: EOMEG ! angular rotation rate (rad.s-1)  

     REAL :: p1000mb     
     REAL :: t0        
     REAL :: p0          
     REAL :: cpovcv       
     REAL :: cvovcp      
     REAL :: rvovrd     
     REAL ::  gma
     REAL ::  EP_1
     REAL ::  EP_2
     REAL ::  rfcp
     REAL ::  p608
     REAL :: rad
     REAL :: reradius  
     Real :: molmas !Molecular mass

    REAL  :: wdaysec        ! duree du sol (s)  ~88775 s
    REAL  :: wmugaz         ! Masse molaire de l'atm (g.mol-1) ~43.49
    REAL  :: womeg          ! omega (rad.s-1)
    REAL  :: wyear_day      ! Duree de l'annee (sols) ~668.6
    REAL  :: wperiheli      ! Dist.min. soleil-mars (Mkm) ~206.66      
    REAL  :: waphelie       ! Dist.max. soleil-mars (Mkm) ~249.22
    REAL  :: wperi_day      ! Date du perihelie (sols depuis printemps)
    REAL  :: wobliquit      ! Obliquite de la planete (deg) ~25.2
    REAL  :: wz0            ! surface roughness (m) ~0.01
    REAL  :: wlmixmin       ! longueur de melange ~100
    REAL  :: wemin_turb     ! energie minimale ~1.e-8 
    REAL  :: wemissiv       ! Emissivite du sol martien ~.95
    REAL  :: wemissiceN     !  Emissivite calotte nord
    REAL  :: wemissiceS     ! Emissivite calotte sud
    REAL  :: walbediceN     ! Albedo calotte nord !0.5
    REAL  :: walbediceS     ! Albedo calotte sud  !0.5
    REAL  :: wiceradiusN    ! mean scat radius of CO2 snow (north)
    REAL  :: wiceradiusS    ! mean scat radius of CO2 snow (south)
    REAL  :: wdtemisiceN    ! time scale for snow metamorphism (north) !2  
    REAL  :: wdtemisiceS    ! time scale for snow metamorphism (south) !2 
    REAL  :: wvolcapa       ! volumetric capacity of soil (new soil model)
    REAL  :: z_scale        ! scale height 

    REAL :: TT00
    REAL :: nu
 CONTAINS
  SUBROUTINE init_module_model_constants
    IMPLICIT NONE
    Call initial_config
    planet=model_config_rec%planet
    !write(*,*) "unknown planet type", planet
    !stop
    
   !  These are the physical constants used within the model
    
    IF ( planet == "mars" ) then 
! JM NOTE -- can we name this grav instead?
     g = 3.72  ! acceleration due to gravity (m {s}^-2)

!#if ( NMM_CORE == 1 )
!    REAL :: r_d          = 192.  ! gas constant m2 s-2 K-1
!    REAL :: cp           = 844.6   ! r= 8.314511E+0 *1000.E+0/mugaz
!#else
!    REAL :: r_d          = 192.
!    REAL :: cp           = 844.6
!#endif

     r_d          = 191.
     cp           = 744.5
     cliq         = 4190.
     cice         = 2106.
     psat         = 610.78
     t0           = 220. ! earth : 300 ... à remplacer par 220 ?
     p0           = 610.
     p1000mb      = 610.
     reradius     = 1./3397200.	! Rayon de mars (m)  ~3397200 m
     EOMEG        =7.0721E-5 ! angular rotation rate (rad.s-1)
     wdaysec       = 88800.     ! duree du sol (s)  ~88775 s
     wmugaz        = 43.49      ! Masse molaire de l'atm (g.mol-1) ~43.49
     womeg         = 7.0721E-5  ! omega (rad.s-1)
     wyear_day     = 669.       ! Duree de l'annee (sols) ~668.6
     wperiheli     = 206.66     ! Dist.min. soleil-mars (Mkm) ~206.66      
     waphelie      = 249.22     ! Dist.max. soleil-mars (Mkm) ~249.22
     wperi_day     = 485.       ! Date du perihelie (sols depuis printemps)
     wobliquit     = 25.2       ! Obliquite de la planete (deg) ~25.2
     wz0           = 1.e-2      ! surface roughness (m) ~0.01
     wlmixmin      = 30.        ! longueur de melange ~100
     wemin_turb    = 1.e-6      ! energie minimale ~1.e-8 
     wemissiv      = 0.95       ! Emissivite du sol martien ~.95
     wemissiceN    = 0.95       ! Emissivite calotte nord
     wemissiceS    = 0.95       ! Emissivite calotte sud
     walbediceN    = 0.65        ! Albedo calotte nord !0.5
     walbediceS    = 0.65        ! Albedo calotte sud  !0.5
     wiceradiusN   = 100.e-6    ! mean scat radius of CO2 snow (north)
     wiceradiusS   = 100.e-6    ! mean scat radius of CO2 snow (south)
     wdtemisiceN   = 0.4  ! time scale for snow metamorphism (north) !2  
     wdtemisiceS   = 0.4  ! time scale for snow metamorphism (south) !2 
#ifdef NEWPHYS
     wvolcapa      = 1.e6       ! volumetric capacity of soil (new soil model)
#endif
     TT00=0.
     nu=0.
!#endif
    ELSE if ( planet == "venus") then 
     g            = 8.87
     r_d          = 191.84383904727036
     !cpp          = 1000.
     cp           = 1000.
     cliq         = 4190.
     psat         = 610.78
     t0           = 735. ! earth : 300
     p0           = 92.e5
     p1000mb      = 92.e5
     reradius     = 1./6051800.	! Rayon de mars (m)  ~6051800 m
     rad          = 6051800.
     wmugaz       = 43.49 
     womeg        = 0.2992549E-06
     wdaysec      = 0.1008707E08
     TT00=460.
     nu=0.35
    ELSE if (( planet == "prescribed").or.(planet == "generic")) then
     open(17,file='planet_constant',form='formatted',status='old')
     rewind(17)
     read(17,*) g
     read(17,*) cp
     read(17,*) molmas
     read(17,*) t0
     read(17,*) p0
     read(17,*) rad
     close(17)
     r_d=8.314/(molmas*1e-3)
     p1000mb=p0
     reradius=1./rad
     TT00=0.
     nu=0.
    Else 
     write(*,*) "unknown planet type:", planet
     stop
    endif
     cv           = cp-r_d
     cpv          = 3.9*r_v
     cvv          = cpv-r_v
     cvpm         = -cv/cp
     rcv          = r_d/cv
     rcp          = r_d/cp
     rovg         = r_d/g
     c2           = cp * rcv
     cpovcv       = cp/(cp-r_d)
     cvovcp       = 1./cpovcv
     rvovrd       = r_v/r_d
     gma          = -r_d*(1.-rcp)*0.5 
     EP_1         = r_v/r_d-1.
     EP_2         = r_d/r_v
     rfcp         =0.25/cp
     p608         =rvovrd-1.
     z_scale      =(r_d*t0/g)*1e-3
    write(*,*) "planet : ", planet
    write(*,*) "g = ", g
    write(*,*) "cp = ", cp
    write(*,*) "r_d = ", r_d
    write(*,*) "t0 = ", t0
    write(*,*) "p0 = ", p0
    write(*,*) "T0 = ", TT00
    write(*,*) "nu = ", nu
  END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
