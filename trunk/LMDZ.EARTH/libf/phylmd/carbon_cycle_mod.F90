MODULE carbon_cycle_mod
! Controle module for the carbon CO2 tracers :
!   - Identification
!   - Get concentrations comming from coupled model or read from file to tracers
!   - Calculate new RCO2 for radiation scheme
!   - Calculate new carbon flux for sending to coupled models (PISCES and ORCHIDEE)
!
! Author : Josefine GHATTAS, Patricia CADULE

  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: carbon_cycle_init, carbon_cycle

! Variables read from parmeter file physiq.def
  LOGICAL, PUBLIC :: carbon_cycle_tr        ! 3D transport of CO2 in the atmosphere, parameter read in conf_phys
!$OMP THREADPRIVATE(carbon_cycle_tr)
  LOGICAL, PUBLIC :: carbon_cycle_cpl       ! Coupling of CO2 fluxes between LMDZ/ORCHIDEE and LMDZ/OCEAN(PISCES) 
!$OMP THREADPRIVATE(carbon_cycle_cpl)

  LOGICAL :: carbon_cycle_emis_comp=.FALSE. ! Calculation of emission compatible
!$OMP THREADPRIVATE(carbon_cycle_emis_comp)

  LOGICAL :: RCO2_inter  ! RCO2 interactive : if true calculate new value RCO2 for the radiation scheme
!$OMP THREADPRIVATE(RCO2_inter)

! Scalare values when no transport, from physiq.def
  REAL :: fos_fuel_s  ! carbon_cycle_fos_fuel dans physiq.def
!$OMP THREADPRIVATE(fos_fuel_s)
  REAL :: emis_land_s ! not yet implemented
!$OMP THREADPRIVATE(emis_land_s)

  REAL :: airetot     ! Total area of the earth surface
!$OMP THREADPRIVATE(airetot)

  INTEGER :: ntr_co2  ! Number of tracers concerning the carbon cycle
!$OMP THREADPRIVATE(ntr_co2)

! fco2_ocn_day : flux CO2 from ocean for 1 day (cumulated) [gC/m2/d]. Allocation and initalization done in cpl_mod
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_ocn_day 
!$OMP THREADPRIVATE(fco2_ocn_day)

  REAL, DIMENSION(:), ALLOCATABLE :: fco2_land_day   ! flux CO2 from land for 1 day (cumulated)  [gC/m2/d]
!$OMP THREADPRIVATE(fco2_land_day)
  REAL, DIMENSION(:), ALLOCATABLE :: fco2_lu_day     ! Emission from land use change for 1 day (cumulated) [gC/m2/d]
!$OMP THREADPRIVATE(fco2_lu_day)

  REAL, DIMENSION(:,:), ALLOCATABLE :: dtr_add       ! Tracer concentration to be injected 
!$OMP THREADPRIVATE(dtr_add)

! Following 2 fields will be allocated and initialized in surf_land_orchidee
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_land_inst  ! flux CO2 from land at one time step
!$OMP THREADPRIVATE(fco2_land_inst)
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: fco2_lu_inst    ! Emission from land use change at one time step
!$OMP THREADPRIVATE(fco2_lu_inst)

! Calculated co2 field to be send to the ocean via the coupler and to ORCHIDEE 
  REAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: co2_send ! Field allocated in phyetat0
!$OMP THREADPRIVATE(co2_send)


  TYPE, PUBLIC ::   co2_trac_type
     CHARACTER(len = 8) :: name       ! Tracer name in tracer.def
     INTEGER            :: id         ! Index in total tracer list, tr_seri
     CHARACTER(len=30)  :: file       ! File name
     LOGICAL            :: cpl        ! True if this tracers is coupled from ORCHIDEE or PISCES. 
                                      ! False if read from file.
     INTEGER            :: updatefreq ! Frequence to inject in second
     INTEGER            :: readstep   ! Actual time step to read in file
     LOGICAL            :: updatenow  ! True if this tracer should be updated this time step
  END TYPE co2_trac_type
  INTEGER,PARAMETER :: maxco2trac=5  ! Maximum number of different CO2 fluxes
  TYPE(co2_trac_type), DIMENSION(maxco2trac) :: co2trac

CONTAINS
  
  SUBROUTINE carbon_cycle_init(tr_seri, pdtphys, aerosol, radio)
! This subroutine is called from traclmdz_init, only at first timestep.
! - Read controle parameters from .def input file
! - Search for carbon tracers and set default values
! - Allocate variables
! - Test for compatibility

    USE dimphy
    USE comgeomphy
    USE mod_phys_lmdz_transfert_para
    USE infotrac
    USE IOIPSL
    USE surface_data, ONLY : ok_veget, type_ocean
    USE phys_cal_mod, ONLY : mth_len

    IMPLICIT NONE
    INCLUDE "clesphys.h"
    INCLUDE "iniprint.h"
 
! Input argument
    REAL,DIMENSION(klon,klev,nbtr),INTENT(IN) :: tr_seri ! Concentration Traceur [U/KgA]  
    REAL,INTENT(IN)                           :: pdtphys ! length of time step in physiq (sec)

! InOutput arguments
    LOGICAL,DIMENSION(nbtr), INTENT(INOUT) :: aerosol
    LOGICAL,DIMENSION(nbtr), INTENT(INOUT) :: radio

! Local variables
    INTEGER               :: ierr, it, iiq, itc
    INTEGER               :: teststop



! 1) Read controle parameters from .def input file
! ------------------------------------------------
    ! Read fosil fuel value if no transport
    IF (.NOT. carbon_cycle_tr) THEN
       fos_fuel_s = 0.
       CALL getin ('carbon_cycle_fos_fuel',fos_fuel_s)
       WRITE(lunout,*) 'carbon_cycle_fos_fuel = ', fos_fuel_s 
    END IF


    ! Read parmeter for calculation compatible emission
    IF (.NOT. carbon_cycle_tr) THEN
       carbon_cycle_emis_comp=.FALSE.
       CALL getin('carbon_cycle_emis_comp',carbon_cycle_emis_comp)
       WRITE(lunout,*) 'carbon_cycle_emis_comp = ',carbon_cycle_emis_comp
       IF (carbon_cycle_emis_comp) THEN
          CALL abort_gcm('carbon_cycle_init', 'carbon_cycle_emis_comp option not yet implemented!!',1)
       END IF
    END IF

    ! Read parameter for interactive calculation of the CO2 value for the radiation scheme
    RCO2_inter=.FALSE.
    CALL getin('RCO2_inter',RCO2_inter)
    WRITE(lunout,*) 'RCO2_inter = ', RCO2_inter
    IF (RCO2_inter) THEN
       WRITE(lunout,*) 'RCO2 will be recalculated once a day'
       WRITE(lunout,*) 'RCO2 initial = ', RCO2
    END IF


! 2) Search for carbon tracers and set default values
! ---------------------------------------------------
    itc=0
    DO it=1,nbtr
       iiq=niadv(it+2)
       
       SELECT CASE(tname(iiq))
       CASE("fCO2_ocn")
          itc = itc + 1
          co2trac(itc)%name='fCO2_ocn'
          co2trac(itc)%id=it
          co2trac(itc)%file='fl_co2_ocean.nc'
          IF (carbon_cycle_cpl .AND. type_ocean=='couple') THEN 
             co2trac(itc)%cpl=.TRUE.
             co2trac(itc)%updatefreq = 86400 ! Once a day as the coupling with OASIS/PISCES
          ELSE
             co2trac(itc)%cpl=.FALSE.
             co2trac(itc)%updatefreq = 86400*mth_len ! Once a month
          END IF
       CASE("fCO2_land")
          itc = itc + 1
          co2trac(itc)%name='fCO2_land'
          co2trac(itc)%id=it
          co2trac(itc)%file='fl_co2_land.nc'
          IF (carbon_cycle_cpl .AND. ok_veget) THEN 
             co2trac(itc)%cpl=.TRUE.
             co2trac(itc)%updatefreq = INT(pdtphys) ! Each timestep as the coupling with ORCHIDEE
          ELSE
             co2trac(itc)%cpl=.FALSE.
!             co2trac(itc)%updatefreq = 10800   ! 10800sec = 3H
             co2trac(itc)%updatefreq = 86400*mth_len ! Once a month
          END IF
       CASE("fCO2_land_use")
          itc = itc + 1
          co2trac(itc)%name='fCO2_land_use'
          co2trac(itc)%id=it
          co2trac(itc)%file='fl_co2_land_use.nc'
          IF (carbon_cycle_cpl .AND. ok_veget) THEN 
             co2trac(it)%cpl=.TRUE.
             co2trac(itc)%updatefreq = INT(pdtphys) ! Each timestep as the coupling with ORCHIDEE
          ELSE
             co2trac(itc)%cpl=.FALSE.
             co2trac(itc)%updatefreq = 10800   ! 10800sec = 3H
          END IF
       CASE("fCO2_fos_fuel")
          itc = itc + 1
          co2trac(itc)%name='fCO2_fos_fuel'
          co2trac(itc)%id=it
          co2trac(itc)%file='fossil_fuel.nc'
          co2trac(itc)%cpl=.FALSE.       ! This tracer always read from file
!         co2trac(itc)%updatefreq = 86400  ! 86400sec = 24H Cadule case
          co2trac(itc)%updatefreq = 86400*mth_len ! Once a month
       CASE("fCO2_bbg")
          itc = itc + 1
          co2trac(itc)%name='fCO2_bbg'
          co2trac(itc)%id=it
          co2trac(itc)%file='fl_co2_bbg.nc'
          co2trac(itc)%cpl=.FALSE.       ! This tracer always read from file
          co2trac(itc)%updatefreq = 86400*mth_len ! Once a month
       CASE("fCO2")
          ! fCO2 : One tracer transporting the total CO2 flux
          itc = itc + 1
          co2trac(itc)%name='fCO2'
          co2trac(itc)%id=it
          co2trac(itc)%file='fl_co2.nc'
          IF (carbon_cycle_cpl) THEN 
             co2trac(itc)%cpl=.TRUE.
          ELSE
             co2trac(itc)%cpl=.FALSE.
          END IF
          co2trac(itc)%updatefreq = 86400
          ! DOES THIS WORK ???? Problematic due to implementation of the coupled fluxes...
          CALL abort_gcm('carbon_cycle_init','transport of total CO2 has to be implemented and tested',1)
       END SELECT
    END DO

    ! Total number of carbon CO2 tracers
    ntr_co2 = itc 
    
    ! Definition of control varaiables for the tracers
    DO it=1,ntr_co2
       aerosol(co2trac(it)%id) = .FALSE.
       radio(co2trac(it)%id)   = .FALSE.
    END DO
    
    ! Vector indicating which timestep to read for each tracer
    ! Always start read in the beginning of the file
    co2trac(:)%readstep = 0
   

! 3) Allocate variables
! ---------------------
    ! Allocate vector for storing fluxes to inject
    ALLOCATE(dtr_add(klon,maxco2trac), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('carbon_cycle_init', 'pb in allocation 11',1)       
    
    ! Allocate variables for cumulating fluxes from ORCHIDEE
    IF (RCO2_inter) THEN
       IF (.NOT. carbon_cycle_tr .AND. carbon_cycle_cpl) THEN
          ALLOCATE(fco2_land_day(klon), stat=ierr)
          IF (ierr /= 0) CALL abort_gcm('carbon_cycle_init', 'pb in allocation 2',1)
          fco2_land_day(1:klon) = 0.
          
          ALLOCATE(fco2_lu_day(klon), stat=ierr)
          IF (ierr /= 0) CALL abort_gcm('carbon_cycle_init', 'pb in allocation 3',1)
          fco2_lu_day(1:klon)   = 0.
       END IF
    END IF


! 4) Test for compatibility
! -------------------------
!    IF (carbon_cycle_cpl .AND. type_ocean/='couple') THEN
!       WRITE(lunout,*) 'Coupling with ocean model is needed for carbon_cycle_cpl'
!       CALL abort_gcm('carbon_cycle_init', 'coupled ocean is needed for carbon_cycle_cpl',1)
!    END IF
!
!    IF (carbon_cycle_cpl .AND..NOT. ok_veget) THEN
!       WRITE(lunout,*) 'Coupling with surface land model ORCHDIEE is needed for carbon_cycle_cpl'
!       CALL abort_gcm('carbon_cycle_init', 'ok_veget is needed for carbon_cycle_cpl',1)
!    END IF

    ! Compiler test : following should never happen
    teststop=0
    DO it=1,teststop
       CALL abort_gcm('carbon_cycle_init', 'Entering loop from 1 to 0',1)
    END DO

    IF (ntr_co2==0) THEN
       ! No carbon tracers found in tracer.def. It is not possible to do carbon cycle 
       WRITE(lunout,*) 'No carbon tracers found in tracer.def. Not ok with carbon_cycle_tr and/or carbon_cycle_cp'
       CALL abort_gcm('carbon_cycle_init', 'No carbon tracers found in tracer.def',1)
    END IF
    
! 5) Calculate total area of the earth surface
! --------------------------------------------
    CALL reduce_sum(SUM(airephy),airetot)
    CALL bcast(airetot)

  END SUBROUTINE carbon_cycle_init


  SUBROUTINE carbon_cycle(nstep, pdtphys, pctsrf, tr_seri, source)
! Subroutine for injection of co2 in the tracers
!
! - Find out if it is time to update
! - Get tracer from coupled model or from file
! - Calculate new RCO2 value for the radiation scheme
! - Calculate CO2 flux to send to ocean and land models (PISCES and ORCHIDEE)

    USE infotrac
    USE dimphy
    USE mod_phys_lmdz_transfert_para
    USE phys_cal_mod, ONLY : mth_cur, mth_len
    USE phys_cal_mod, ONLY : day_cur
    USE comgeomphy

    IMPLICIT NONE

    INCLUDE "clesphys.h"
    INCLUDE "indicesol.h"
    INCLUDE "iniprint.h"
    INCLUDE "YOMCST.h"

! In/Output arguments
    INTEGER,INTENT(IN) :: nstep      ! time step in physiq
    REAL,INTENT(IN)    :: pdtphys    ! length of time step in physiq (sec)
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf            ! Surface fraction
    REAL, DIMENSION(klon,klev,nbtr), INTENT(INOUT)  :: tr_seri ! All tracers
    REAL, DIMENSION(klon,nbtr), INTENT(INOUT)       :: source  ! Source for all tracers

! Local variables
    INTEGER :: it
    LOGICAL :: newmonth ! indicates if a new month just started
    LOGICAL :: newday   ! indicates if a new day just started
    LOGICAL :: endday   ! indicated if last time step in a day

    REAL, PARAMETER :: fact=1.E-15/2.12  ! transformation factor from gC/m2/day => ppm/m2/day
    REAL, DIMENSION(klon) :: fco2_tmp
    REAL :: sumtmp
    REAL :: delta_co2_ppm
    

! 1) Calculate logicals indicating if it is a new month, new day or the last time step in a day (end day)
! -------------------------------------------------------------------------------------------------------

    newday = .FALSE.; endday = .FALSE.; newmonth = .FALSE.

    IF (MOD(nstep,INT(86400./pdtphys))==1) newday=.TRUE.
    IF (MOD(nstep,INT(86400./pdtphys))==0) endday=.TRUE.
    IF (newday .AND. day_cur==1) newmonth=.TRUE.

! 2)  For each carbon tracer find out if it is time to inject (update)
! --------------------------------------------------------------------
    DO it = 1, ntr_co2
       IF ( MOD(nstep,INT(co2trac(it)%updatefreq/pdtphys)) == 1 ) THEN
          co2trac(it)%updatenow = .TRUE.
       ELSE
          co2trac(it)%updatenow = .FALSE.
       END IF
    END DO

! 3) Get tracer update
! --------------------------------------
    DO it = 1, ntr_co2
       IF ( co2trac(it)%updatenow ) THEN
          IF ( co2trac(it)%cpl ) THEN
             ! Get tracer from coupled model
             SELECT CASE(co2trac(it)%name)
             CASE('fCO2_land')     ! from ORCHIDEE
                dtr_add(:,it) = fco2_land_inst(:)*pctsrf(:,is_ter)*fact ! [ppm/m2/day]
             CASE('fCO2_land_use') ! from ORCHIDEE
                dtr_add(:,it) = fco2_lu_inst(:)  *pctsrf(:,is_ter)*fact ! [ppm/m2/day]
             CASE('fCO2_ocn')      ! from PISCES
                dtr_add(:,it) = fco2_ocn_day(:)  *pctsrf(:,is_oce)*fact ! [ppm/m2/day]
             CASE DEFAULT
                WRITE(lunout,*) 'Error with tracer ',co2trac(it)%name
                CALL abort_gcm('carbon_cycle', 'No coupling implemented for this tracer',1)
             END SELECT
          ELSE
             ! Read tracer from file
             co2trac(it)%readstep = co2trac(it)%readstep + 1 ! increment time step in file
! Patricia   CALL read_map2D(co2trac(it)%file,'fco2',co2trac(it)%readstep,.FALSE.,dtr_add(:,it))
             CALL read_map2D(co2trac(it)%file,'fco2',co2trac(it)%readstep,.TRUE.,dtr_add(:,it))

             ! Converte from kgC/m2/h to kgC/m2/s
             dtr_add(:,it) = dtr_add(:,it)/3600
             ! Add individual treatment of values read from file
             SELECT CASE(co2trac(it)%name)
             CASE('fCO2_land')
                dtr_add(:,it) = dtr_add(:,it) *pctsrf(:,is_ter)
             CASE('fCO2_land_use')
                dtr_add(:,it) = dtr_add(:,it) *pctsrf(:,is_ter)
             CASE('fCO2_ocn')
                dtr_add(:,it) = dtr_add(:,it) *pctsrf(:,is_oce)
! Patricia :
!             CASE('fCO2_fos_fuel')
!                dtr_add(:,it) = dtr_add(:,it)/mth_len
!                co2trac(it)%readstep = 0 ! Always read same value for fossil fuel(Cadule case)
             END SELECT
          END IF
       END IF
    END DO

! 4) Update co2 tracers : 
!    Loop over all carbon tracers and add source
! ------------------------------------------------------------------
    IF (carbon_cycle_tr) THEN
       DO it = 1, ntr_co2
          IF (.FALSE.) THEN
             tr_seri(1:klon,1,co2trac(it)%id) = tr_seri(1:klon,1,co2trac(it)%id) + dtr_add(1:klon,it)
             source(1:klon,co2trac(it)%id) = 0.
          ELSE
             source(1:klon,co2trac(it)%id) = dtr_add(1:klon,it)
          END IF
       END DO
    END IF


! 5) Calculations for new CO2 value for the radiation scheme(instead of reading value from .def)
! ----------------------------------------------------------------------------------------------
    IF (RCO2_inter) THEN
       ! Cumulate fluxes from ORCHIDEE at each timestep
       IF (.NOT. carbon_cycle_tr .AND. carbon_cycle_cpl) THEN
          IF (newday) THEN ! Reset cumulative variables once a day 
             fco2_land_day(1:klon) = 0.
             fco2_lu_day(1:klon)   = 0.
          END IF
          fco2_land_day(1:klon) = fco2_land_day(1:klon) + fco2_land_inst(1:klon) ![gC/m2/day]
          fco2_lu_day(1:klon)   = fco2_lu_day(1:klon)   + fco2_lu_inst(1:klon)   ![gC/m2/day]
       END IF

       ! At the end of a new day, calculate a mean scalare value of CO2
       ! JG : Ici on utilise uniquement le traceur du premier couche du modele. Est-ce que c'est correcte ? 
       IF (endday) THEN

          IF (carbon_cycle_tr) THEN
             ! Sum all co2 tracers to get the total delta CO2 flux
             fco2_tmp(:) = 0.
             DO it = 1, ntr_co2
                fco2_tmp(1:klon) = fco2_tmp(1:klon) + tr_seri(1:klon,1,co2trac(it)%id)
             END DO
             
          ELSE IF (carbon_cycle_cpl) THEN ! no carbon_cycle_tr
             ! Sum co2 fluxes comming from coupled models and parameter for fossil fuel
             fco2_tmp(1:klon) = fos_fuel_s + ((fco2_lu_day(1:klon) + fco2_land_day(1:klon))*pctsrf(1:klon,is_ter) &
                  + fco2_ocn_day(:)*pctsrf(:,is_oce)) * fact
          END IF

          ! Calculate a global mean value of delta CO2 flux
          fco2_tmp(1:klon) = fco2_tmp(1:klon) * airephy(1:klon)
          CALL reduce_sum(SUM(fco2_tmp),sumtmp)
          CALL bcast(sumtmp)
          delta_co2_ppm = sumtmp/airetot
          
          ! Add initial value for co2_ppm and delta value
          co2_ppm = co2_ppm0 + delta_co2_ppm
          
          ! Transformation of atmospheric CO2 concentration for the radiation code
          RCO2 = co2_ppm * 1.0e-06  * 44.011/28.97 
          
          WRITE(lunout,*) 'RCO2 is now updated! RCO2 = ', RCO2
       END IF ! endday

    END IF ! RCO2_inter


! 6) Calculate CO2 flux to send to ocean and land models : PISCES and ORCHIDEE         
! ----------------------------------------------------------------------------
    IF (carbon_cycle_cpl) THEN

       IF (carbon_cycle_tr) THEN
          ! Sum all co2 tracers to get the total delta CO2 flux at first model layer
          fco2_tmp(:) = 0.
          DO it = 1, ntr_co2
             fco2_tmp(1:klon) = fco2_tmp(1:klon) + tr_seri(1:klon,1,co2trac(it)%id)
          END DO
          co2_send(1:klon) = fco2_tmp(1:klon) + co2_ppm0
       ELSE
          ! Send a scalare value in 2D variable to ocean and land model (PISCES and ORCHIDEE)
          co2_send(1:klon) = co2_ppm
       END IF

    END IF

  END SUBROUTINE carbon_cycle
  
END MODULE carbon_cycle_mod
