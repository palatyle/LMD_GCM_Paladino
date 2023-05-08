program fft

! SL 01/2010:
! This program reads 4D (lon-lat-alt-time) fields recast in log P coordinates
!
! it computes fft of temperature, zonal and merid winds from high-frequency outputs:
!
! fftaT -- 4D -- FFT in amplitude of temperature field (K)
! fftau -- 4D -- FFT in amplitude of zonal wind (m s-1)
! fftav -- 4D -- FFT in amplitude of meridional wind (m s-1)
! ulf   -- 4D -- low  freq part of zonal wind perturbation uprim (m s-1)
! ubf   -- 4D -- band freq part of zonal wind perturbation uprim (m s-1)
! uhf   -- 4D -- high freq part of zonal wind perturbation uprim (m s-1)
! vlf   -- 4D -- low  freq part of meridional wind perturbation vprim (m s-1)
! vbf   -- 4D -- band freq part of meridional wind perturbation vprim (m s-1)
! vhf   -- 4D -- high freq part of meridional wind perturbation vprim (m s-1)
! wlf   -- 4D -- low  freq part of vertical wind perturbation wprim (Pa s-1)
! wbf   -- 4D -- band freq part of vertical wind perturbation wprim (Pa s-1)
! whf   -- 4D -- high freq part of vertical wind perturbation wprim (Pa s-1)
! Tlf   -- 4D -- low  freq part of temperature perturbation Tprim (K)
! Tbf   -- 4D -- band freq part of temperature perturbation Tprim (K)
! Thf   -- 4D -- high freq part of temperature perturbation Tprim (K)
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - pressure vertical coordinate
! - atmospheric temperature
! - zonal, meridional and vertical winds
!
! We use the FFTW library:   http://www.fftw.org
! These routines are in C, but also include Fortran interfaces.
!
! Convention: qbar  <=> zonal average    / qstar = q - qbar
!             qmean <=> temporal average / qprim = q - qmean
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  FILTRES
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  low  frequencies: qlf= lowfreq(qprim)
!  band frequencies: qbf=bandfreq(qprim)
!  high frequencies: qhf=highfreq(qprim)
!
!  Les frequences seuils sont ajustables dans filter.h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (name_P.nc)
character (len=128) :: outfile1,outfile2,outfile3,outfile4 ! output file names

character (len=64) :: text ! to store some text
integer infid ! NetCDF input file ID
integer outfid1,outfid2,outfid3,outfid4 ! NetCDF output files ID
integer lon_dimid1,lat_dimid1,alt_dimid1,time_dimid1 ! NetCDF dimension IDs
integer lon_dimid2,lat_dimid2,alt_dimid2,time_dimid2 ! NetCDF dimension IDs
integer lon_dimid3,lat_dimid3,alt_dimid3,time_dimid3 ! NetCDF dimension IDs
integer lon_dimid4,lat_dimid4,alt_dimid4,time_dimid4 ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid
integer              :: datashape1d ! shape of 1D datasets
integer,dimension(4) :: datashape4d ! shape of 4D datasets

real :: miss_val ! special "missing value" to specify missing data
real,parameter :: miss_val_def=-9.99e+33 ! default value for "missing value"
real :: pi
real,dimension(:),allocatable :: lon ! longitude
integer lonlength ! # of grid points along longitude
real,dimension(:),allocatable :: lat ! latitude
integer latlength ! # of grid points along latitude
real,dimension(:),allocatable :: plev ! Pressure levels (Pa)
integer altlength ! # of grid point along altitude (of input datasets)
real,dimension(:),allocatable :: time ! time
real,dimension(:),allocatable :: freq ! frequencies of the FFT (only timelength/2+1 values)
integer timelength ! # of points along time
real,dimension(:,:,:,:),allocatable :: temp ! atmospheric temperature
real,dimension(:,:,:,:),allocatable :: vitu ! zonal wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitw ! vertical wind (in Pa/s)

!!! output variables
real,dimension(:,:,:,:),allocatable :: fftaT ! FFT in amplitude of temperature (K)
real,dimension(:,:,:,:),allocatable :: fftau ! FFT in amplitude of zonal wind (m s-1)
real,dimension(:,:,:,:),allocatable :: fftav ! FFT in amplitude of meridional wind (m s-1)
real,dimension(:,:,:,:),allocatable :: fftaw ! FFT in amplitude of vertical wind (Pa s-1)
real,dimension(:,:,:,:),allocatable :: ulf ! low  freq part of zonal wind perturbation uprim (m s-1)
real,dimension(:,:,:,:),allocatable :: ubf ! band freq part of zonal wind perturbation uprim (m s-1)
real,dimension(:,:,:,:),allocatable :: uhf ! high freq part of zonal wind perturbation uprim (m s-1)
real,dimension(:,:,:,:),allocatable :: vlf ! low  freq part of meridional wind perturbation vprim (m s-1)
real,dimension(:,:,:,:),allocatable :: vbf ! band freq part of meridional wind perturbation vprim (m s-1)
real,dimension(:,:,:,:),allocatable :: vhf ! high freq part of meridional wind perturbation vprim (m s-1)
real,dimension(:,:,:,:),allocatable :: wlf ! low  freq part of vertical wind perturbation vprim (Pa s-1)
real,dimension(:,:,:,:),allocatable :: wbf ! band freq part of vertical wind perturbation vprim (Pa s-1)
real,dimension(:,:,:,:),allocatable :: whf ! high freq part of vertical wind perturbation vprim (Pa s-1)
real,dimension(:,:,:,:),allocatable :: Tlf ! low  freq part of temperature perturbation Tprim (K)
real,dimension(:,:,:,:),allocatable :: Tbf ! band freq part of temperature perturbation Tprim (K)
real,dimension(:,:,:,:),allocatable :: Thf ! high freq part of temperature perturbation Tprim (K)

! local variables
real,dimension(:,:,:,:),allocatable :: uprim
real,dimension(:,:,:,:),allocatable :: vprim
real,dimension(:,:,:,:),allocatable :: wprim
real,dimension(:,:,:,:),allocatable :: Tprim 

! lon,lat,alt
real,dimension(:,:,:),allocatable :: umean 
real,dimension(:,:,:),allocatable :: vmean 
real,dimension(:,:,:),allocatable :: wmean 
real,dimension(:,:,:),allocatable :: Tmean 

! for FFTW routines
real,dimension(:),allocatable :: wndow
double precision,dimension(:),allocatable :: var,fltvar
double complex,dimension(:),allocatable :: fftvar,fltfft
double complex,dimension(:),allocatable :: filtrelf,filtrebf,filtrehf
integer   :: M_fft
integer*8 :: planf,planb


integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
logical flagfft
logical :: lmdflag

! Tuning parameters
real :: fcoup1,fcoup2,width
real :: fcoup1tmp,fcoup2tmp,widthtmp
logical,dimension(4) :: ok_out
character (len=1) :: ok_outtmp

include "planet.h"

#include <fftw3.f> 

!===============================================================================
! 1. Input parameters
!===============================================================================

pi = 2.*asin(1.)
miss_val = miss_val_def

write(*,*) ""
write(*,*) "You are working on the atmosphere of ",planet

! initialisation
!----------------

! Par defaut

! Define the filters
! Low  cutting frequency, in Hz    : fcoup1
fcoup1=2.5e-6
! High cutting frequency, in Hz    : fcoup2
fcoup2=6.5e-6
! Half-width of the filters, in Hz : width
width=4.e-7
! Outputs (U,     V,      W,     T)
ok_out=(/.true.,.true.,.false.,.true./)

print*,"Low  cutting frequency, in Hz ?"
print*,"between 1e-5 and 1e-7, 0 for default => 2.5e-6"
read(*,*) fcoup1tmp
if ((fcoup1tmp.lt.1e-5).and.(fcoup1tmp.gt.1e-7)) fcoup1=fcoup1tmp
print*,"=",fcoup1
print*,"High cutting frequency, in Hz ?"
print*,"between 1e-5 and 1e-7, 0 for default => 6.5e-6"
read(*,*) fcoup2tmp
if ((fcoup2tmp.lt.1e-5).and.(fcoup2tmp.gt.1e-7)) fcoup2=fcoup2tmp
print*,"=",fcoup2
print*,"Half-width of the filters, in Hz ?"
print*,"between 1e-6 and 1e-8, 0 for default => 4e-7)"
read(*,*) widthtmp
if ((widthtmp.lt.1e-6).and.(widthtmp.gt.1e-8)) width=widthtmp
print*,"=",width
!width = 3./time(timelength)

! Outputs 
print*,"Output of zonal wind ? (y or n, default is y)"
read(*,'(a1)') ok_outtmp
if (ok_outtmp.eq."n") ok_out(1)=.false.
print*,"=",ok_out(1)
print*,"Output of meridional wind ? (y or n, default is y)"
read(*,'(a1)') ok_outtmp
if (ok_outtmp.eq."n") ok_out(2)=.false.
print*,"=",ok_out(2)
print*,"Output of vertical wind ? (y or n, default is n)"
read(*,'(a1)') ok_outtmp
if (ok_outtmp.eq."y") ok_out(3)=.true.
print*,"=",ok_out(3)
print*,"Output of temperature ? (y or n, default is y)"
read(*,'(a1)') ok_outtmp
if (ok_outtmp.eq."n") ok_out(4)=.false.
print*,"=",ok_out(4)

!===============================================================================
! 1.1 Input file
!===============================================================================

write(*,*) ""
write(*,*) " Program valid for files with pressure axis (*_P.nc)"
write(*,*) "Enter input file name:"

read(*,'(a128)') infile
write(*,*) ""

! open input file

ierr = NF_OPEN(infile,NF_NOWRITE,infid)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening file ',trim(infile)
   stop ""
endif

!===============================================================================
! 1.2 Get grids in lon,lat,alt(pressure),time
!===============================================================================

call get_iddim(infid,lat_varid,latlength,lon_varid,lonlength,&
                     alt_varid,altlength,time_varid,timelength,lmdflag )

allocate(lon(lonlength))
ierr=NF_GET_VAR_REAL(infid,lon_varid,lon)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading longitude"

allocate(lat(latlength))
ierr=NF_GET_VAR_REAL(infid,lat_varid,lat)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading lat"

allocate(plev(altlength))
ierr=NF_GET_VAR_REAL(infid,alt_varid,plev)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading altitude (ie pressure levels)"

allocate(time(timelength))
ierr=NF_GET_VAR_REAL(infid,time_varid,time)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading time"

!===============================================================================
! 1.3 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
outfile1=infile(1:len_trim(infile)-3)//"_UFFT.nc"
outfile2=infile(1:len_trim(infile)-3)//"_VFFT.nc"
outfile3=infile(1:len_trim(infile)-3)//"_WFFT.nc"
outfile4=infile(1:len_trim(infile)-3)//"_TFFT.nc"
write(*,*) "Output file names are: "
if (ok_out(1)) write(*,*) trim(outfile1)
if (ok_out(2)) write(*,*) trim(outfile2)
if (ok_out(3)) write(*,*) trim(outfile3)
if (ok_out(4)) write(*,*) trim(outfile4)


!===============================================================================
! 2.1 Store needed fields 
!===============================================================================

!===============================================================================
! 2.1.1 Atmospheric temperature
!===============================================================================
if (ok_out(4)) then
allocate(temp(lonlength,latlength,altlength,timelength))

text="temp"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for t instead... "
  text="t"
  call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) then
    print*,"Error: Failed to get temperature ID"
    ok_out(4)=.false.
  endif
endif
if (ierr2.ne.NF_NOERR) then
  print*,"Error: Failed reading temperature"
  ok_out(4)=.false.
endif
endif !ok_out(4)

!===============================================================================
! 2.1.2 Winds
!===============================================================================
! zonal wind vitu (in m/s)
if (ok_out(1)) then
allocate(vitu(lonlength,latlength,altlength,timelength))

text="vitu"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitu,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  print*,"Error: Failed to get vitu ID"
  ok_out(1)=.false.
endif
if (ierr2.ne.NF_NOERR) then
  print*,"Error: Failed reading zonal wind"
  ok_out(1)=.false.
endif
endif !ok_out(1)

! meridional wind vitv (in m/s)
if (ok_out(2)) then
allocate(vitv(lonlength,latlength,altlength,timelength))

text="vitv"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitv,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  print*,"Error: Failed to get vitv ID"
  ok_out(2)=.false.
endif
if (ierr2.ne.NF_NOERR) then
  print*,"Error: Failed reading meridional wind"
  ok_out(2)=.false.
endif
endif !ok_out(2)

! vertical wind vitw (in Pa/s)
if (ok_out(3)) then
allocate(vitw(lonlength,latlength,altlength,timelength))

text="vitw"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitw,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  print*,"Error: Failed to get vitw ID"
  ok_out(3)=.false.
endif
if (ierr2.ne.NF_NOERR) then
  print*,"Error: Failed reading vertical wind"
  ok_out(3)=.false.
endif
endif !ok_out(3)

!===============================================================================
! 2.2 Computations 
!===============================================================================

print*,"debut calcul"

!===============================================================================
! 2.2.1 FFT and filtering
!===============================================================================

! allocations
!-------------
if (ok_out(1)) then
allocate(fftau(lonlength,latlength,altlength,timelength))
allocate(uprim(lonlength,latlength,altlength,timelength))
allocate(ulf(lonlength,latlength,altlength,timelength))
allocate(ubf(lonlength,latlength,altlength,timelength))
allocate(uhf(lonlength,latlength,altlength,timelength))
endif !ok_out(1)
if (ok_out(2)) then
allocate(fftav(lonlength,latlength,altlength,timelength))
allocate(vprim(lonlength,latlength,altlength,timelength))
allocate(vlf(lonlength,latlength,altlength,timelength))
allocate(vbf(lonlength,latlength,altlength,timelength))
allocate(vhf(lonlength,latlength,altlength,timelength))
endif !ok_out(2)
if (ok_out(3)) then
allocate(fftaw(lonlength,latlength,altlength,timelength))
allocate(wprim(lonlength,latlength,altlength,timelength))
allocate(wlf(lonlength,latlength,altlength,timelength))
allocate(wbf(lonlength,latlength,altlength,timelength))
allocate(whf(lonlength,latlength,altlength,timelength))
endif !ok_out(3)
if (ok_out(4)) then
allocate(fftaT(lonlength,latlength,altlength,timelength))
allocate(Tprim(lonlength,latlength,altlength,timelength))
allocate(Tlf(lonlength,latlength,altlength,timelength))
allocate(Tbf(lonlength,latlength,altlength,timelength))
allocate(Thf(lonlength,latlength,altlength,timelength))
endif !ok_out(4)

! lon,lat,alt
if (ok_out(1)) allocate(umean(lonlength,latlength,altlength))
if (ok_out(2)) allocate(vmean(lonlength,latlength,altlength))
if (ok_out(3)) allocate(wmean(lonlength,latlength,altlength))
if (ok_out(4)) allocate(Tmean(lonlength,latlength,altlength))

! time / frequencies
allocate(freq(timelength))
allocate(wndow(timelength))
allocate(var(timelength))
allocate(fltvar(timelength))
M_fft = timelength/2
allocate(fftvar(M_fft+1))
allocate(fltfft(M_fft+1))
allocate(filtrelf(M_fft+1))
allocate(filtrebf(M_fft+1))
allocate(filtrehf(M_fft+1))

! intermediates
!-----------------

if (ok_out(1)) call moytim(lonlength,latlength,altlength,timelength,miss_val,vitu,umean)
if (ok_out(2)) call moytim(lonlength,latlength,altlength,timelength,miss_val,vitv,vmean)
if (ok_out(3)) call moytim(lonlength,latlength,altlength,timelength,miss_val,vitw,wmean)
if (ok_out(4)) call moytim(lonlength,latlength,altlength,timelength,miss_val,temp,Tmean)

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
   do itim=1,timelength
if (ok_out(1)) then
    if ((vitu(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (umean(ilon,ilat,ilev)    .ne.miss_val)) then
  uprim(ilon,ilat,ilev,itim) = vitu(ilon,ilat,ilev,itim)-umean(ilon,ilat,ilev)
    else
  uprim(ilon,ilat,ilev,itim) = miss_val
    endif
endif !ok_out(1)
if (ok_out(2)) then
    if ((vitv(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (vmean(ilon,ilat,ilev)    .ne.miss_val)) then
  vprim(ilon,ilat,ilev,itim) = vitv(ilon,ilat,ilev,itim)-vmean(ilon,ilat,ilev)
    else
  vprim(ilon,ilat,ilev,itim) = miss_val
    endif
endif !ok_out(2)
if (ok_out(3)) then
    if ((vitw(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (wmean(ilon,ilat,ilev)    .ne.miss_val)) then
  wprim(ilon,ilat,ilev,itim) = vitw(ilon,ilat,ilev,itim)-wmean(ilon,ilat,ilev)
    else
  wprim(ilon,ilat,ilev,itim) = miss_val
    endif
endif !ok_out(3)
if (ok_out(4)) then
    if ((temp(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (Tmean(ilon,ilat,ilev)    .ne.miss_val)) then
  Tprim(ilon,ilat,ilev,itim) = temp(ilon,ilat,ilev,itim)-Tmean(ilon,ilat,ilev)
    else
  Tprim(ilon,ilat,ilev,itim) = miss_val
    endif
endif !ok_out(4)
   enddo
  enddo
 enddo
enddo ! lonlength

! fft intermediates
!-------------

! Define the frequencies
do itim=1,M_fft+1
  freq(itim) = (itim-1)/(timelength*(time(2)-time(1)))
enddo
do itim=M_fft+2,timelength
  freq(itim) = 0.
enddo

! Define the window (triangle)
do itim=1,timelength
! N window:
!  wndow(itim)= 1.
! triangulaire de moyenne = 1
  wndow(itim)= 2.*(1. - abs(real(itim-0.5-M_fft)/real(M_fft)))
enddo

do itim=1,M_fft+1
  if (freq(itim).lt.(fcoup1-width)) then
     filtrelf(itim) = 1.
  elseif (freq(itim).gt.(fcoup1+width)) then
     filtrelf(itim) = 0.
  else
     filtrelf(itim) = (1.+sin(pi*(fcoup1-freq(itim))/(2.*width)))/2.
  endif
  if (freq(itim).lt.(fcoup2-width)) then
     filtrehf(itim) = 0.
  elseif (freq(itim).gt.(fcoup2+width)) then
     filtrehf(itim) = 1.
  else
     filtrehf(itim) = (1.-sin(pi*(fcoup2-freq(itim))/(2.*width)))/2.
  endif
  filtrebf(itim) = (1.-filtrelf(itim))*(1.-filtrehf(itim))
enddo


! fft and filtering
!-------------

!---FFTW routines
call dfftw_plan_dft_r2c_1d(planf,timelength,var,fftvar,FFTW_MEASURE)
call dfftw_plan_dft_c2r_1d(planb,timelength,fltfft,fltvar,FFTW_MEASURE)
!---

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength

! For zonal wind field
if (ok_out(1)) then

   flagfft=.true.
   do itim=1,timelength
     if (uprim(ilon,ilat,ilev,itim).eq.miss_val) flagfft=.false.
   enddo

   if (flagfft) then

! 1/ windowing to improve spectral analysis
      var(:)=uprim(ilon,ilat,ilev,:)*wndow(:)
! 2/ FFT computation
!---FFTW routines
      call dfftw_execute_dft_r2c(planf,var,fftvar)
!---
! 3/ Amplitude of the FFT, for spectral analysis
      fftau(ilon,ilat,ilev,1)=abs(fftvar(1))/M_fft
      do itim=2,M_fft
       fftau(ilon,ilat,ilev,itim) = abs(fftvar(itim))/M_fft
      enddo
      fftau(ilon,ilat,ilev,M_fft+1)=abs(fftvar(M_fft+1))/M_fft
      do itim=M_fft+2,timelength
       fftau(ilon,ilat,ilev,itim) = 0.
      enddo

! 4/ filtering the FFT in three regions
! filtering + normalisation (low freq)
      fltfft(:) = fftvar(:)*filtrelf(:)/timelength
! 5/ backward FFT for each region
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
! 6/ reverse the windowing
      ulf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (band freq)
      fltfft(:) = fftvar(:)*filtrebf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      ubf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (high freq)
      fltfft(:) = fftvar(:)*filtrehf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      uhf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

   else
     fftau(ilon,ilat,ilev,itim) = miss_val
       ulf(ilon,ilat,ilev,itim) = miss_val
       ubf(ilon,ilat,ilev,itim) = miss_val
       uhf(ilon,ilat,ilev,itim) = miss_val
   endif ! flagfft

endif !ok_out(1)

! For meridional wind wind field
if (ok_out(2)) then

   flagfft=.true.
   do itim=1,timelength
     if (vprim(ilon,ilat,ilev,itim).eq.miss_val) flagfft=.false.
   enddo

   if (flagfft) then

! 1/ windowing to improve spectral analysis
      var(:)=vprim(ilon,ilat,ilev,:)*wndow(:)
! 2/ FFT computation
!---FFTW routines
      call dfftw_execute_dft_r2c(planf,var,fftvar)
!---
! 3/ Amplitude of the FFT, for spectral analysis
      fftav(ilon,ilat,ilev,1)=abs(fftvar(1))/M_fft
      do itim=2,M_fft
       fftav(ilon,ilat,ilev,itim) = abs(fftvar(itim))/M_fft
      enddo
      fftav(ilon,ilat,ilev,M_fft+1)=abs(fftvar(M_fft+1))/M_fft
      do itim=M_fft+2,timelength
       fftav(ilon,ilat,ilev,itim) = 0.
      enddo

! 4/ filtering the FFT in three regions
! filtering + normalisation (low freq)
      fltfft(:) = fftvar(:)*filtrelf(:)/timelength
! 5/ backward FFT for each region
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
! 6/ reverse the windowing
      vlf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (band freq)
      fltfft(:) = fftvar(:)*filtrebf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      vbf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (high freq)
      fltfft(:) = fftvar(:)*filtrehf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      vhf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

   else
     fftav(ilon,ilat,ilev,itim) = miss_val
       vlf(ilon,ilat,ilev,itim) = miss_val
       vbf(ilon,ilat,ilev,itim) = miss_val
       vhf(ilon,ilat,ilev,itim) = miss_val
   endif ! flagfft

endif !ok_out(2)

! For vertical wind wind field
if (ok_out(3)) then

   flagfft=.true.
   do itim=1,timelength
     if (wprim(ilon,ilat,ilev,itim).eq.miss_val) flagfft=.false.
   enddo

   if (flagfft) then

! 1/ windowing to improve spectral analysis
      var(:)=wprim(ilon,ilat,ilev,:)*wndow(:)
! 2/ FFT computation
!---FFTW routines
      call dfftw_execute_dft_r2c(planf,var,fftvar)
!---
! 3/ Amplitude of the FFT, for spectral analysis
      fftaw(ilon,ilat,ilev,1)=abs(fftvar(1))/M_fft
      do itim=2,M_fft
       fftaw(ilon,ilat,ilev,itim) = abs(fftvar(itim))/M_fft
      enddo
      fftaw(ilon,ilat,ilev,M_fft+1)=abs(fftvar(M_fft+1))/M_fft
      do itim=M_fft+2,timelength
       fftaw(ilon,ilat,ilev,itim) = 0.
      enddo

! 4/ filtering the FFT in three regions
! filtering + normalisation (low freq)
      fltfft(:) = fftvar(:)*filtrelf(:)/timelength
! 5/ backward FFT for each region
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
! 6/ reverse the windowing
      wlf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (band freq)
      fltfft(:) = fftvar(:)*filtrebf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      wbf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (high freq)
      fltfft(:) = fftvar(:)*filtrehf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      whf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

   else
     fftaw(ilon,ilat,ilev,itim) = miss_val
       wlf(ilon,ilat,ilev,itim) = miss_val
       wbf(ilon,ilat,ilev,itim) = miss_val
       whf(ilon,ilat,ilev,itim) = miss_val
   endif ! flagfft

endif !ok_out(3)

! For temperature field
if (ok_out(4)) then

   flagfft=.true.
   do itim=1,timelength
     if (Tprim(ilon,ilat,ilev,itim).eq.miss_val) flagfft=.false.
   enddo

   if (flagfft) then

! 1/ windowing to improve spectral analysis
      var(:)=Tprim(ilon,ilat,ilev,:)*wndow(:)
! 2/ FFT computation
!---FFTW routines
      call dfftw_execute_dft_r2c(planf,var,fftvar)
!---
! 3/ Amplitude of the FFT, for spectral analysis
      fftaT(ilon,ilat,ilev,1)=abs(fftvar(1))/M_fft
      do itim=2,M_fft
       fftaT(ilon,ilat,ilev,itim) = abs(fftvar(itim))/M_fft
      enddo
      fftaT(ilon,ilat,ilev,M_fft+1)=abs(fftvar(M_fft+1))/M_fft
      do itim=M_fft+2,timelength
       fftaT(ilon,ilat,ilev,itim) = 0.
      enddo

! 4/ filtering the FFT in three regions
! filtering + normalisation (low freq)
      fltfft(:) = fftvar(:)*filtrelf(:)/timelength
! 5/ backward FFT for each region
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
! 6/ reverse the windowing
      Tlf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (band freq)
      fltfft(:) = fftvar(:)*filtrebf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      Tbf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

! filtering + normalisation (high freq)
      fltfft(:) = fftvar(:)*filtrehf(:)/timelength
!---FFTW routines
      call dfftw_execute_dft_c2r(planb,fltfft,fltvar)
!---
      Thf(ilon,ilat,ilev,:) = fltvar(:)/wndow(:)

   else
     fftaT(ilon,ilat,ilev,itim) = miss_val
       Tlf(ilon,ilat,ilev,itim) = miss_val
       Tbf(ilon,ilat,ilev,itim) = miss_val
       Thf(ilon,ilat,ilev,itim) = miss_val
   endif ! flagfft

endif !ok_out(4)

  enddo
 enddo
enddo ! lonlength

!---FFTW routines
call dfftw_destroy_plan(planf)
call dfftw_destroy_plan(planb)
!---

print*,"End of computations"

!===============================================================================
! 3. Create output files
!===============================================================================

! Create output files
if (ok_out(1)) then
ierr=NF_CREATE(outfile1,NF_CLOBBER,outfid1)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile1
  stop
endif
endif !ok_out(1)

if (ok_out(2)) then
ierr=NF_CREATE(outfile2,NF_CLOBBER,outfid2)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile2
  stop
endif
endif !ok_out(2)

if (ok_out(3)) then
ierr=NF_CREATE(outfile3,NF_CLOBBER,outfid3)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile3
  stop
endif
endif !ok_out(3)

if (ok_out(4)) then
ierr=NF_CREATE(outfile4,NF_CLOBBER,outfid4)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile4
  stop
endif
endif !ok_out(4)

!===============================================================================
! 3.1. Define and write dimensions
!===============================================================================

if (ok_out(1)) &
call write_dim(outfid1,lonlength,latlength,altlength,timelength, &
    lon,lat,plev,time,lon_dimid1,lat_dimid1,alt_dimid1,time_dimid1)
if (ok_out(2)) &
call write_dim(outfid2,lonlength,latlength,altlength,timelength, &
    lon,lat,plev,time,lon_dimid2,lat_dimid2,alt_dimid2,time_dimid2)
if (ok_out(3)) &
call write_dim(outfid3,lonlength,latlength,altlength,timelength, &
    lon,lat,plev,time,lon_dimid3,lat_dimid3,alt_dimid3,time_dimid3)
if (ok_out(4)) &
call write_dim(outfid4,lonlength,latlength,altlength,timelength, &
    lon,lat,plev,time,lon_dimid4,lat_dimid4,alt_dimid4,time_dimid4)

!===============================================================================
! 3.2. Define and write variables
!===============================================================================

if (ok_out(1)) then

datashape4d(1)=lon_dimid1
datashape4d(2)=lat_dimid1
datashape4d(3)=alt_dimid1
datashape4d(4)=time_dimid1
datashape1d   =time_dimid1

call write_var1d(outfid1,datashape1d,timelength,&
                "freq      ", "FFT frequencies     ","s-1       ",miss_val,&
                 freq )

call write_var4d(outfid1,datashape4d,lonlength,latlength,altlength,timelength,&
                 "fftau     ", "FFT ampl of vitu    ","m s-1     ",miss_val,&
                  fftau )

call write_var4d(outfid1,datashape4d,lonlength,latlength,altlength,timelength,&
                 "ulf       ", "low freq part vitu  ","m s-1     ",miss_val,&
                  ulf )

call write_var4d(outfid1,datashape4d,lonlength,latlength,altlength,timelength,&
                 "ubf       ", "band freq part vitu ","m s-1     ",miss_val,&
                  ubf )

call write_var4d(outfid1,datashape4d,lonlength,latlength,altlength,timelength,&
                 "uhf       ", "high freq part vitu ","m s-1     ",miss_val,&
                  uhf )
endif !ok_out(1)

if (ok_out(2)) then

datashape4d(1)=lon_dimid2
datashape4d(2)=lat_dimid2
datashape4d(3)=alt_dimid2
datashape4d(4)=time_dimid2
datashape1d   =time_dimid2

call write_var1d(outfid2,datashape1d,timelength,&
                "freq      ", "FFT frequencies     ","s-1       ",miss_val,&
                 freq )

call write_var4d(outfid2,datashape4d,lonlength,latlength,altlength,timelength,&
                 "fftav     ", "FFT ampl of vitv    ","m s-1     ",miss_val,&
                  fftav )

call write_var4d(outfid2,datashape4d,lonlength,latlength,altlength,timelength,&
                 "vlf       ", "low freq part vitv  ","m s-1     ",miss_val,&
                  vlf )

call write_var4d(outfid2,datashape4d,lonlength,latlength,altlength,timelength,&
                 "vbf       ", "band freq part vitv ","m s-1     ",miss_val,&
                  vbf )

call write_var4d(outfid2,datashape4d,lonlength,latlength,altlength,timelength,&
                 "vhf       ", "high freq part vitv ","m s-1     ",miss_val,&
                  vhf )
endif !ok_out(2)

if (ok_out(3)) then

datashape4d(1)=lon_dimid3
datashape4d(2)=lat_dimid3
datashape4d(3)=alt_dimid3
datashape4d(4)=time_dimid3
datashape1d   =time_dimid3

call write_var1d(outfid3,datashape1d,timelength,&
                "freq      ", "FFT frequencies     ","s-1       ",miss_val,&
                 freq )

call write_var4d(outfid3,datashape4d,lonlength,latlength,altlength,timelength,&
                 "fftaw     ", "FFT ampl of vitw    ","Pa s-1    ",miss_val,&
                  fftaw )

call write_var4d(outfid3,datashape4d,lonlength,latlength,altlength,timelength,&
                 "wlf       ", "low freq part vitw  ","Pa s-1    ",miss_val,&
                  wlf )

call write_var4d(outfid3,datashape4d,lonlength,latlength,altlength,timelength,&
                 "wbf       ", "band freq part vitw ","Pa s-1    ",miss_val,&
                  wbf )

 call write_var4d(outfid3,datashape4d,lonlength,latlength,altlength,timelength,&
                 "whf       ", "high freq part vitw ","Pa s-1    ",miss_val,&
                  whf )
endif !ok_out(3)

if (ok_out(4)) then

datashape4d(1)=lon_dimid4
datashape4d(2)=lat_dimid4
datashape4d(3)=alt_dimid4
datashape4d(4)=time_dimid4
datashape1d   =time_dimid4

call write_var1d(outfid4,datashape1d,timelength,&
                "freq      ", "FFT frequencies     ","s-1       ",miss_val,&
                 freq )

call write_var4d(outfid4,datashape4d,lonlength,latlength,altlength,timelength,&
                 "fftaT     ", "FFT ampl of temp    ","K         ",miss_val,&
                  fftaT )

call write_var4d(outfid4,datashape4d,lonlength,latlength,altlength,timelength,&
                 "tlf       ", "low freq part temp  ","K         ",miss_val,&
                  Tlf )

call write_var4d(outfid4,datashape4d,lonlength,latlength,altlength,timelength,&
                 "tbf       ", "band freq part temp ","K         ",miss_val,&
                  Tbf )

call write_var4d(outfid4,datashape4d,lonlength,latlength,altlength,timelength,&
                 "thf       ", "high freq part temp ","K         ",miss_val,&
                  Thf )
endif !ok_out(4)

!!!! Close output files
if (ok_out(1)) then
ierr=NF_CLOSE(outfid1)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile1
endif !ok_out(1)

if (ok_out(2)) then
ierr=NF_CLOSE(outfid2)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile2
endif !ok_out(2)

if (ok_out(3)) then
ierr=NF_CLOSE(outfid3)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile3
endif !ok_out(3)

if (ok_out(4)) then
ierr=NF_CLOSE(outfid4)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile4
endif !ok_out(4)


end program
