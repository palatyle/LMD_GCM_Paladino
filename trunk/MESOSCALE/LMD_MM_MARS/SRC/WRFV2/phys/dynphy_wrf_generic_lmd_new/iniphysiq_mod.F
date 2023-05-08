MODULE iniphysiq_mod

CONTAINS

subroutine iniphysiq(ngrid,nlayer,nq,piphysiq,&
                     punjours, pdayref, &
                     prad,pg,pr,pcpp,iflag_phys)

!use control_mod, only: nday
!use surf_heat_transp_mod, only: ini_surf_heat_transp
!use infotrac, only : nqtot ! number of advected tracers
!USE comvert_mod, ONLY: ap,bp,preff
use inifis_mod, only: inifis
use ioipsl_getin_p_mod, only: getin_p

!use inigeomphy_mod, only: inigeomphy
!use geometry_mod, only: cell_area, & ! physics grid area (m2)
!                        longitude, & ! longitudes (rad)
!                        latitude ! latitudes (rad)
! necessary to get klon_omp
!USE mod_phys_lmdz_para, ONLY: klon_omp ! number of columns (on local omp grid)
USE mod_phys_lmdz_para, ONLY: Init_phys_lmdz_para
USE dimphy, ONLY: init_dimphy
USE phys_state_var_mod
!use planete_mod, only: year_day, periastr, apoastr, peri_day,&
!                       obliquit, z0, lmixmin, emin_turb
!                       init_planete_mod
use planete_mod
use time_phylmdz_mod, only: dtphys, daysec,day_ini
use planete_mod, only: year_day, periastr, apoastr, peri_day,&
                   obliquit, z0, lmixmin, emin_turb
use surfdat_h,  only: emissiv,iceradius, &
                    emisice,dtemisice
use comcstfi_mod, only: omeg,mugaz
use tracer_h, only: nqtot
!use comm_wrf, only : allocate_comm_wrf

implicit none

INCLUDE 'mpif.h'

REAL,intent(in) :: prad
REAL,intent(in) :: pg
REAL,intent(in) :: pr
REAL,intent(in) :: pcpp
REAL,intent(in) :: punjours
!DOUBLE PRECISION,intent(in) :: ptimestep

!real,intent(in) :: prad ! radius of the planet (m)
!real,intent(in) :: pg ! gravitational acceleration (m/s2)
!real,intent(in) :: pr ! ! reduced gas constant R/mu
!real,intent(in) :: pcpp ! specific heat Cp
!real,intent(in) :: punjours ! length (in s) of a standard day [daysec]
integer,intent(in) :: pdayref ! reference day of for the simulation [day_ini]
integer,intent(in) :: iflag_phys ! type of physics to be called

integer :: nday=0 ! this is dummy for mesoscale (in dyn3d/control_mod)

integer,intent(in) :: ngrid ! number of physics columns for this MPI process
integer,intent(in) :: nlayer ! number of atmospheric layers
integer,intent(in) :: nq ! number of tracers
!real,intent(in) :: phour_ini   ! start time (fraction of day) of the run
!0=<phour_ini<1
real,intent(in) :: piphysiq   ! call physics every piphysiq dynamical timesteps
!real :: latitude(ngrid),longitude(ngrid),cell_area(ngrid)
!  real,intent(in) :: prefff ! reference surface pressure (Pa)
!  real,intent(in) :: apf(nlayer+1) ! hybrid coordinate at interfaces
!  real,intent(in) :: bpf(nlayer+1) 
logical :: ok_slab_ocean
real*8 :: lat(ngrid),long(ngrid),cellarea(ngrid)
REAL*8 :: pprad,ppg,ppr,ppcpp,ppunjours
REAL*8 :: dummy 
  ! the common part for all planetary physics
  !------------------------------------------
  ! --> initialize physics distribution, global fields and geometry
  ! (i.e. things in phy_common or dynphy_lonlat)

  ! the distinct part for all planetary physics (ie. things in phystd)
  !------------------------------------------

CALL Init_phys_lmdz_para(1,1,1,MPI_COMM_WORLD)

!call phys_state_var_init
print*,'ngrid',ngrid,'nlayer',nlayer
call init_dimphy(ngrid,nlayer)
call phys_state_var_init(nqtot)
! copy over preff , ap() and bp()
!call ini_planete_mod(nlayer,prefff,apf,bpf)

! for slab ocean, copy over some arrays
ok_slab_ocean=.false. ! default value
!call getin_p("ok_slab_ocean",ok_slab_ocean)
!if (ok_slab_ocean) then
!  call ini_surf_heat_transp(ip1jm,ip1jmp1,unsairez,fext,unsaire, &
!                            cu,cuvsurcv,cv,cvusurcu,aire,apoln,apols, &
!                            aireu,airev)
!endif

dummy=1.
lat(:)=0.
long(:)=0.
cellarea(:)=1.
print*,'pg',pg
!ppunjours=punjours
ppunjours=1.
pprad=prad
ppg=pg
ppr=pr
ppcpp=pcpp
call inifis(ngrid,nlayer,nq,pdayref,ppunjours,nday,dummy, &
            lat,long,cellarea,pprad,ppg,ppr,ppcpp)

     open(17,file='controle.txt',form='formatted',status='old')
     rewind(17)
     read(17,*)
     read(17,*)
     read(17,*) day_ini !(tab0+3)
     read(17,*)
     read(17,*) !tab0+5)
     read(17,*) omeg !(tab0+6)
     read(17,*) !(tab0+7)
     read(17,*) mugaz
     read(17,*)  !(tab0+9)
     read(17,*) daysec
     read(17,*) dtphys !tab0+11)
     read(17,*)
     read(17,*)
     read(17,*) year_day !(tab0+14)
     read(17,*) periastr !tab0+15)
     read(17,*) apoastr !tab0+16)
     read(17,*) peri_day !tab0+17)
     read(17,*) obliquit !tab0+18)
     read(17,*) z0
     read(17,*)
     read(17,*)
     read(17,*)
     read(17,*)
     read(17,*) emisice(1)
     read(17,*) emisice(2)
     read(17,*) emissiv
     read(17,*)
     read(17,*)
     read(17,*)
     read(17,*)
     read(17,*) iceradius(1)
     read(17,*) iceradius(2)
     read(17,*) dtemisice(1)
     read(17,*) dtemisice(2)
     close(17)
     !print*,'g',g

     !emissiv(:)=EMIS
     cloudfrac(:,:)=0.5
     totcloudfrac(:)=0.5
     hice(:)=0.
     rnat(:)=0.
     pctsrf_sic(:)=0.
     tsea_ice(:)=0.
     !qsurf(:,:) = 0. 
     print*,'check'
     print*,'iceradius',iceradius,'dtemisice',dtemisice
     print*,'apoastr,periastr,year_day,peri_day,obliq',apoastr,periastr,year_day,peri_day,obliquit
     print*,'emissiv',emissiv
     print*,'mugaz',mugaz
     
end subroutine iniphysiq


END MODULE iniphysiq_mod
