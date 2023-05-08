program kcm1d

  use infotrac, only: nqtot
  use radinc_h,      only: NAERKIND
  use radcommon_h,   only: L_NSPECTI, L_NSPECTV, sigma, glat
  use watercommon_h, only: mH2O
  use ioipsl_getincom, only: getin
  use comsaison_h, only: mu0, fract, dist_star
  use planete_mod
  use callkeys_mod, only: check_cpp_match, pceil, tstrat, tracer
  use inifis_mod, only: inifis
  use comcstfi_mod
  use mod_grid_phy_lmdz, only : regular_lonlat
  use regular_lonlat_mod, only: init_regular_lonlat
  use physics_distribution_mod, only: init_physics_distribution
  use regular_lonlat_mod, only: init_regular_lonlat
  use mod_interface_dyn_phys, only: init_interface_dyn_phys
  use geometry_mod, only: init_geometry
  use dimphy, only : init_dimphy

  implicit none

  !==================================================================
  !     
  !     Purpose
  !     -------
  !     Run the universal model radiative transfer once in a 1D column.
  !     Useful for climate evolution studies etc.
  !     
  !     It can be compiled with a command like (e.g. for 25 layers):
  !                                  "makegcm -p std -d 25 kcm1d"
  !     It requires the files "callphys.def", "gases.def"
  !     "traceur.def", and "run.def" with a line "INCLUDEDEF=callphys.def"
  !
  !     Authors
  !     -------
  !     - R. Wordsworth
  !     - updated by M. Turbet (June 2017)
  !
  !==================================================================

#include "dimensions.h"
!#include "dimphys.h"

  ! --------------------------------------------------------------
  !  Declarations
  ! --------------------------------------------------------------

  integer nlayer,nlevel,nq
  integer ilay,ilev,iq,iw,iter
  real play(llm)     ! Pressure at the middle of the layers [Pa]
  real zlay(llm)     ! Altitude at middle of the layers [km]
  real plev(llm+1)   ! Intermediate pressure levels [Pa]
  real temp(llm)     ! temperature at the middle of the layers [K]
  real,allocatable :: q(:,:)   ! tracer mixing ratio [kg/kg]
  real,allocatable :: vmr(:,:) ! tracer mixing ratio [mol/mol]
  real,allocatable :: qsurf(:)        ! tracer surface budget [kg/kg] ????
  real psurf,psurf_n,tsurf      

  real emis, albedo, albedo_equivalent
  real albedo_wv(1,L_NSPECTV)

  real muvar(llm+1)

  real dtsw(llm) ! heating rate (K/s) due to SW
  real dtlw(llm) ! heating rate (K/s) due to LW
  real fluxsurf_lw   ! incident LW flux to surf (W/m2)
  real fluxtop_lw    ! outgoing LW flux to space (W/m2)
  real fluxsurf_sw   ! incident SW flux to surf (W/m2)
  real fluxsurfabs_sw ! absorbed SW flux by the surf (W/m2)
  real fluxabs_sw    ! SW flux absorbed by planet (W/m2)
  real fluxtop_dn    ! incident top of atmosphere SW flux (W/m2)

  ! not used
  real cloudfrac(llm)
  real totcloudfrac
  real tau_col

  real dTstrat
  real aerosol(llm,naerkind) ! aerosol tau (kg/kg)
  real OLR_nu(1,L_NSPECTI)
  real OSR_nu(1,L_NSPECTV)
  real Eatmtot

  integer ierr
  logical firstcall,lastcall,global1d

  character*20,allocatable :: nametrac(:)   ! name of the tracer (no need for adv trac common)
  
  real :: latitude(1), longitude(1), cell_area(1), phisfi(1)
  
  ! --------------
  ! Initialisation
  ! --------------

  pi=2.E+0*asin(1.E+0)

  cloudfrac(:)  = 0.0
  totcloudfrac  = 0.0


  nlayer=llm
  nlevel=nlayer+1

  !! this is defined in comsaison_h
  ALLOCATE(mu0(1))
  ALLOCATE(fract(1))
  ALLOCATE(glat(1))

  !  Load parameters from "run.def" 
  ! -------------------------------

  ! check if 'kcm1d.def' file is around (otherwise reading parameters
  ! from callphys.def via getin() routine won't work.)
  open(90,file='kcm1d.def',status='old',form='formatted',&
       iostat=ierr)
  if (ierr.ne.0) then
     write(*,*) 'Cannot find required file "kcm1d.def"'
     write(*,*) '  (which should contain some input parameters'
     write(*,*) '   along with the following line:'
     write(*,*) '   INCLUDEDEF=callphys.def'
     write(*,*) '   )'
     write(*,*) ' ... might as well stop here ...'
     stop
  else
     close(90)
  endif
  
! now, run.def is needed anyway. so we create a dummy temporary one
! for ioipsl to work. if a run.def is already here, stop the
! program and ask the user to do a bit of cleaning
  open(90,file='run.def',status='old',form='formatted',& 
       iostat=ierr)
  if (ierr.eq.0) then
     close(90)
     write(*,*) 'There is already a run.def file.'
     write(*,*) 'This is not compatible with 1D runs.'
     write(*,*) 'Please remove the file and restart the run.'
     write(*,*) 'Runtime parameters are supposed to be in kcm1d.def'
     stop
  else
     call system('touch run.def')
     call system("echo 'INCLUDEDEF=callphys.def' >> run.def")
     call system("echo 'INCLUDEDEF=kcm1d.def' >> run.def")
  endif
    
! check if we are going to run with or without tracers
  write(*,*) "Run with or without tracer transport ?"
  tracer=.false. ! default value
  call getin("tracer",tracer)
  write(*,*) " tracer = ",tracer



  ! Tracer initialisation
  ! ---------------------
  if (tracer) then
     ! load tracer names from file 'traceur.def'
     open(90,file='traceur.def',status='old',form='formatted',&
          iostat=ierr)
     if (ierr.eq.0) then
        write(*,*) "kcm1d: Reading file traceur.def"
        ! read number of tracers:
        read(90,*,iostat=ierr) nq
        if (ierr.ne.0) then
           write(*,*) "kcm1d: error reading number of tracers"
           write(*,*) "   (first line of traceur.def) "
           stop
        endif
        nqtot=nq
        ! allocate arrays which depend on number of tracers
        allocate(nametrac(nq))
        allocate(q(nlayer,nq))
        allocate(vmr(nlayer,nq))
        allocate(qsurf(nq))

        do iq=1,nq
           ! minimal version, just read in the tracer names, 1 per line
           read(90,*,iostat=ierr) nametrac(iq)
           write(*,*) 'tracer here is', nametrac(iq)
           if (ierr.ne.0) then
              write(*,*) 'kcm1d: error reading tracer names...'
              stop
           endif
        enddo !of do iq=1,nq
     endif

     call initracer(1,nq,nametrac)

  endif
  



  global1d = .false. ! default value
  call getin("global1d",global1d)
  write(*,*) " global1d = ",global1d
  if(.not.global1d)then
     print*,'Error, kcm1d must have global1d=.true. in kcm1d.def!'
     stop
  end if

  psurf_n=100000. ! default value for psurf
  print*,'Dry surface pressure (Pa)?'
  call getin("psurf",psurf_n)
  write(*,*) " psurf = ",psurf_n

! OK. now that run.def has been read once -- any variable is in memory.
! so we can dump the dummy run.def
  call system("rm -rf run.def")

  tsurf=300.0 ! default value for tsurf
  print*,'Surface temperature (K)?'
  call getin("tref",tsurf)
  write(*,*) " tsurf = ",tsurf

  g=10.0 ! default value for g
  print*,'Gravity ?'
  call getin("g",g)
  write(*,*) " g = ",g
  glat(1)=g

  periastr = 1.0
  apoastr  = 1.0
  print*,'Periastron (AU)?'
  call getin("periastr",periastr)
  write(*,*) "periastron = ",periastr
  dist_star = periastr 
  fract     = 0.5
  print*,'Apoastron (AU)?'
  call getin("apoastr",apoastr)
  write(*,*) "apoastron = ",apoastr

  albedo=0.2 ! default value for albedo
  print*,'Albedo of bare ground?'
  call getin("albedo",albedo)
  write(*,*) " albedo = ",albedo
  do iw=1,L_NSPECTV
     albedo_wv(1,iw)=albedo
  enddo

  emis=1.0 ! default value for emissivity
  PRINT *,'Emissivity of bare ground ?'
  call getin("emis",emis)
  write(*,*) " emis = ",emis

  pceil=100.0 ! Pascals
  PRINT *,'Ceiling pressure (Pa) ?'
  call getin("pceil",pceil)
  write(*,*) " pceil = ", pceil

  check_cpp_match = .false.
  call getin("check_cpp_match",check_cpp_match)
  if (check_cpp_match) then
     print*,"In 1D modeling, check_cpp_match is supposed to be F"
     print*,"Please correct callphys.def"
     stop
  endif


!!! GEOGRAPHICAL INITIALIZATIONS
     !!! AREA. useless in 1D
  cell_area(1)=1.E+0 !JL+EM to have access to the area in the diagfi.nc files.
     !!! GEOPOTENTIAL. useless in 1D because control by surface pressure
  phisfi(1)=0.E+0
     !!! LATITUDE. can be set. 
  latitude=0 ! default value for latitude
  PRINT *,'latitude (in degrees) ?'
  call getin("latitude",latitude)
  write(*,*) " latitude = ",latitude
  latitude=latitude*pi/180.E+0
     !!! LONGITUDE. useless in 1D.
  longitude=0.E+0
  longitude=longitude*pi/180.E+0

  rad=6400000.
  !rad = -99999.
  !PRINT *,'PLANETARY RADIUS in m ?'
  !call getin("rad",rad)
  ! Planetary  radius is needed to compute shadow of the rings
  !IF (rad.eq.-99999. .and. rings_shadow .eqv. .true.) THEN
  !   PRINT *,"STOP. I NEED rad IN RCM1D.DEF."
  !   STOP
  !ELSE
  !   PRINT *,"--> rad = ",rad
  !ENDIF



  call init_physics_distribution(regular_lonlat,4,1,1,1,nlayer,1)
  call init_interface_dyn_phys
  CALL init_regular_lonlat(1,1,longitude,latitude,(/0.,0./),(/0.,0./))
  call init_geometry(1,longitude,latitude,(/0.,0.,0.,0./),(/0.,0.,0.,0./),cell_area)
  call init_dimphy(1,nlayer) ! Initialize dimphy module

!!! CALL INIFIS
!!! - read callphys.def
!!! - calculate sine and cosine of longitude and latitude
!!! - calculate mugaz and cp 
!!! NB: some operations are being done dummily in inifis in 1D
!!! - physical constants: nevermind, things are done allright below
!!! - physical frequency: nevermind, in inifis this is a simple print

  cpp=1.d-7 !JL because we divide by cpp in inifis, there may be a more elegant solution
  CALL inifis(1,llm,nq,0,86400.0,1,1.0,latitude,longitude,cell_area,rad,g,r,cpp)

  !write(*,*) 'BASE 1'
  !write(*,*) 1,llm,nq,0,86400.0,1,1.0,latitude,longitude,cell_area,rad,g,r,cpp

  do iq=1,nq
     do ilay=1,nlayer
        q(ilay,iq) = 0.
     enddo
  enddo
  
  do iq=1,nq
     qsurf(iq) = 0.
  enddo

  firstcall = .true.
  lastcall  = .false.

  iter    = 1
  Tstrat  = 150.0
  dTstrat = 1000.0

  ! ---------
  ! Run model
  ! ---------
  !do
  psurf = psurf_n

     !    Create vertical profiles
  call kcmprof_fn(nlayer,psurf,qsurf(1),tsurf,    &
       tstrat,play,plev,zlay,temp,q(:,1),muvar(1))
         

  !write(*,*) 'BASE 2'
  !write(*,*) nlayer,psurf,qsurf(1),tsurf,    &
  !           tstrat
  !write(*,*) play
  !write(*,*) plev
  !write(*,*) zlay
  !write(*,*) temp
  !write(*,*) q(:,1),muvar(1)

         
               
     call iniaerosol()
     

     

     !    Run radiative transfer
     call callcorrk(1,nlayer,q,nq,qsurf,                  &
          albedo_wv,albedo_equivalent,                    &
          emis,mu0,plev,play,temp,                        &
          tsurf,fract,dist_star,aerosol,muvar,            &
          dtlw,dtsw,fluxsurf_lw,fluxsurf_sw,              & 
          fluxsurfabs_sw,fluxtop_lw,                      &
          fluxabs_sw,fluxtop_dn,OLR_nu,OSR_nu,tau_col,    &
          cloudfrac,totcloudfrac,.false.,firstcall,lastcall)

     !write(*,*) 'BASE 3'
     !write(*,*) 1,nlayer,'A',q,'B',nq,'C',qsurf,                      &
     !     albedo_wv,'D',albedo_equivalent,'E',                    &
     !     emis,'F',mu0,'G',plev,'H',play,'I',temp,'J',                        &
     !     tsurf,'K',fract,'L',dist_star,'M',aerosol,'N',muvar,'O',            &
     !     dtlw,'P',dtsw,'Q',fluxsurf_lw,'R',fluxsurf_sw,'S',              & 
     !     fluxsurfabs_sw,'T',fluxtop_lw,'U',                      &
     !     fluxabs_sw,'V',fluxtop_dn,'W',OLR_nu,'X',OSR_nu,'Y',tau_col,'Z',    &
     !     cloudfrac,'A1',totcloudfrac,'A2',.false.,'A3',firstcall,'A4',lastcall







     !    Iterate stratospheric temperature
     print*,'Tstrat = ',Tstrat
     dTstrat = Tstrat
     !Tstrat  = Tsurf*(0.2786*(psurf/100000.)**(-1.123))**0.25 
     ! skin temperature (gray approx.) using analytic pure H2 expression
     !Tstrat  = (fluxabs_sw/(2*sigma))**0.25 ! skin temperature (gray approx.)
     Tstrat  = (fluxtop_lw/(2*sigma))**0.25 ! skin temperature (gray approx.)
     dTstrat = dTstrat-Tstrat

     !if(abs(dTstrat).lt.1.0)then
     !   print*,'dTstrat = ',dTstrat
     !   exit
     !endif

     !iter=iter+1
     !if(iter.eq.100)then
     !   print*,'Stratosphere failed to converge after'
     !   print*,'100 iteration, aborting run.'
     !   call abort
     !endif

  !end do


  ! Run radiative transfer one last time to get OLR,OSR
  firstcall=.false.
  lastcall=.true.
  call callcorrk(1,nlayer,q,nq,qsurf,                          &
       albedo_wv,albedo_equivalent,emis,mu0,plev,play,temp,    &
       tsurf,fract,dist_star,aerosol,muvar,                    &
       dtlw,dtsw,fluxsurf_lw,fluxsurf_sw,fluxsurfabs_sw,       &
       fluxtop_lw, fluxabs_sw,fluxtop_dn,OLR_nu,OSR_nu,        &
       tau_col,                                                &
       cloudfrac,totcloudfrac,.false.,firstcall,lastcall)


     !write(*,*) 'BASE 4'
     !write(*,*) 1,nlayer,'A',q,'B',nq,'C',qsurf,                      &
     !     albedo_wv,'D',albedo_equivalent,'E',                    &
     !     emis,'F',mu0,'G',plev,'H',play,'I',temp,'J',                        &
     !     tsurf,'K',fract,'L',dist_star,'M',aerosol,'N',muvar,'O',            &
     !     dtlw,'P',dtsw,'Q',fluxsurf_lw,'R',fluxsurf_sw,'S',              & 
     !     fluxsurfabs_sw,'T',fluxtop_lw,'U',                      &
     !     fluxabs_sw,'V',fluxtop_dn,'W',OLR_nu,'X',OSR_nu,'Y',tau_col,'Z',    &
     !     cloudfrac,'A1',totcloudfrac,'A2',.false.,'A3',firstcall,'A4',lastcall







  ! Calculate total atmospheric energy
  Eatmtot=0.0
  !  call calcenergy_kcm(nlayer,tsurf,temp,play,plev,qsurf,&
  !     q(:,1),muvar,Eatmtot)

  ! ------------------------
  ! Save data to ascii files
  ! ------------------------

  print*,'Saving profiles...'
  open(115,file='profpres.out',form='formatted')
  open(116,file='proftemp.out',form='formatted')
  open(117,file='profztab.out',form='formatted')
  open(118,file='profqvar.out',form='formatted')
  open(119,file='profvvar.out',form='formatted')

  write(115,*) psurf
  write(116,*) tsurf
  write(117,*) 0.0
  write(118,*) qsurf(1)
  write(119,*) qsurf(1)*(muvar(1)/mH2O)
  do ilay=1,nlayer
     vmr(ilay,1) = q(ilay,1)*(muvar(ilay+1)/mH2O)
     write(115,*) play(ilay)
     write(116,*) temp(ilay)
     write(117,*) zlay(ilay)
     write(118,*) q(ilay,1)
     write(119,*) vmr(ilay,1)
  enddo
  close(115)
  close(116)
  close(117)
  close(118)
  close(119)

  print*, tsurf,psurf,fluxtop_dn,fluxabs_sw,fluxtop_lw 

  print*,'Saving scalars...'
  open(116,file='surf_vals.out')
  write(116,*) tsurf,psurf,fluxtop_dn,         &
       fluxabs_sw,fluxtop_lw 
  close(116)
  open(111,file='ene_vals.out')
  write(111,*) tsurf,psurf,Eatmtot,Tstrat
  close(111)

  print*,'Saving spectra...'
  open(117,file='OLRnu.out')
  do iw=1,L_NSPECTI
     write(117,*) OLR_nu(1,iw)
  enddo
  close(117)

  open(127,file='OSRnu.out')
  do iw=1,L_NSPECTV
     write(127,*) OSR_nu(1,iw)
  enddo
  close(127)  

end program kcm1d
