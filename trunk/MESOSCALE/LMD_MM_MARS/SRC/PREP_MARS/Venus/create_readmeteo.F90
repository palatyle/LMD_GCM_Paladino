program create_readmeteo

implicit none
include "netcdf.inc"

!------------------------------------------------------------------------!
! create_readmeteo generates the file 'readmeteo.def'                    !
! ... must be run prior to readmeteo                                     !
! ... the user will be asked a few questions                             !
!                                                                        !
! A. Spiga - 01/08/2007                                                  !
!------------------------------------------------------------------------!

INTEGER, PARAMETER :: MONTHS_PER_YEAR = 100
INTEGER, PARAMETER :: mday = 24 
INTEGER :: start_day,init,i,month,day
INTEGER :: ierr,nid,nvarid
INTEGER :: n,start,my,interval_subs,subs
INTEGER :: start_hour,interval,hour,inc_hour
INTEGER :: no_please
INTEGER :: timedim,timelen
INTEGER :: inc
REAL, DIMENSION(100) :: param
REAL, DIMENSION(:), ALLOCATABLE :: time


!
! Init
!
interval = 0

!mday(:)=24
!
! Open input NETCDF file
!
write(*,*) "Scanning netcdf file ..."
ierr=NF_OPEN ("input_diagfi.nc",NF_NOWRITE,nid)
!ierr=NF90_OPEN("input_diagfi.nc",NF_NOWRITE,nid)
IF (ierr.NE.NF_NOERR) THEN
      write(*,*)'**** Please create a symbolic link called input_diagfi.nc'
      CALL ABORT
ENDIF

ierr=NF_INQ_DIMID(nid,"time_counter",timedim)
IF (ierr .NE. NF_NOERR) THEN
    ierr=NF_INQ_DIMID(nid,"time",timedim)
ENDIF
ierr=NF_INQ_DIMLEN(nid,timedim,timelen)


!
! Get starting time
!
ALLOCATE(time(timelen))
ierr = NF_INQ_VARID (nid, "time_counter",nvarid)
IF (ierr .NE. NF_NOERR) THEN
        ierr = NF_INQ_VARID (nid, "time",nvarid)
           IF (ierr .NE. NF_NOERR) THEN
                PRINT *, "Error: Readmeteo <Time> not found"
                stop
           ENDIF
ENDIF
#ifdef NC_DOUBLE
ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif

ierr = NF_INQ_VARID (nid,"controle",nvarid)
IF (ierr .NE. NF_NOERR) THEN
        PRINT *, "Error: Readmeteo <ps> not found"
        !stop
ENDIF
#ifdef NC_DOUBLE
ierr = NF_GET_VAR_DOUBLE(nid, nvarid, param)
#else
ierr = NF_GET_VAR_REAL(nid, nvarid, param)
#endif

write(*,*) "Done."


!
! fill the include if no info can be found in the diagfi
!
include "fix_no_info.inc"


! beware, param(4) is the day reference of start and startfi
! ...have to add time(1) to get the real starting date in diagfi
!start_day=floor(param(4)+time(1))		
!start_hour=nint((param(4)-floor(param(4))+time(1))*13)	! starting hour	
!start_hour=MOD(start_hour,24)   
!IF (interval .eq. 0) interval=nint(time(1)*24)	 ! interval between each time subscript	

start_day=0
start_hour=0
interval=864 ! in hours, 864 HH = 1 MM = 1/100 Vd
PRINT *,'*****************'
PRINT *,'GCM data file starts at sol ',start_day,'and hour',start_hour
PRINT *,'GCM data interval is ',interval,'hours'

month=0
day=0
!interval=864
!
! User defined parameters
!
write(*,*) "This will erase readmeteo.def file"
write(*,*) "  NB: Default mode is to generate "
write(*,*) "      a WPS-compatible file for each timestep "
write(*,*) "      included in the diagfi "
write(*,*) "Total number of generated files: ",timelen 
write(*,*) "Continue ? 0 if no, 1 if yes, 99 for settings"
read(*,*) no_please
if (no_please == 0) stop
if (no_please == 1) then
        my=1000
        n=timelen
        start=0
        interval_subs=1
else
        write(*,*) "-- GCM data file information --"
        write(*,*) "Starting Martian year ? ex: 24,25,26..."
        read(*,*) my
        my=1000+my
        write(*,*) "-- WRF data file information --"
        write(*,*) "How many files do you want to create ? at least 2, max is",timelen
        read(*,*) n
        IF (n == timelen) THEN
                start=1
                interval_subs=1
        ELSE
                write(*,*) "Time subscript you want to begin with (first is 1) ? max is",timelen-n+1
                read(*,*) start
                write(*,*) "Time subscript interval you want to get ? no more than",(timelen-start+1)/n
                read(*,*) interval_subs
        END IF
endif



!
! Info
!
PRINT *,'-------'
PRINT *,'run readmeteo with the command line:'
PRINT *,'readmeteo.exe < readmeteo.def'


!
! Generate readmeteo.def
!
OPEN(6,file='readmeteo.def',status='replace',form='formatted')
! files
IF (n < 10) THEN
        write(6,fmt='(I1)') n
ELSE IF (n < 100) THEN
        write(6,fmt='(I2)') n
ELSE IF (n < 1000) THEN
        write(6,fmt='(I3)') n
ELSE
        write(6,fmt='(I4)') n        
END IF
! subscript in GCM data file
DO i=1,n
        subs=start+(i-1)*interval_subs
        IF (subs < 10) THEN
                write(6,fmt='(I1)') subs 
        ELSE IF (subs < 100) THEN 
                write(6,fmt='(I2)') subs
        ELSE IF (subs < 1000) THEN
                write(6,fmt='(I3)') subs
        ELSE 
                write(6,fmt='(I4)') subs   
        END IF        
END DO
! WRF time reference
hour=interval
inc_hour=interval*interval_subs
IF (hour >= 36) day=day+INT(hour/36)
hour=MOD(hour,36)
IF (day >= mday) THEN
      inc=INT(day/mday)
      day=day-inc*mday
      month=month+inc
END IF
IF (month >= MONTHS_PER_YEAR) THEN
      my=my+1
      month=0
END IF
DO i=1,n
        write(6,fmt='(I4)') my
        IF (month < 10) THEN
                write(6,fmt='(I1,I1)') 0,month 
        ELSE 
                write(6,fmt='(I2)') month
        END IF  
        IF (day < 10) THEN
                write(6,fmt='(I1,I1)') 0,day
        ELSE
                write(6,fmt='(I2)') day
        END IF
        IF (hour < 10) THEN
                write(6,fmt='(I1,I1)') 0,hour
        ELSE
                write(6,fmt='(I2)') hour
        END IF
        write(6,fmt='(A1)') 'y'
IF (hour+inc_hour >= 36) day=day+INT((hour+inc_hour)/36)
hour=MOD(hour+inc_hour,36)
IF (day >= mday) THEN
        inc=INT(day/mday)
        day=day-inc*mday
        month=month+inc
END IF        
IF (month >= MONTHS_PER_YEAR) THEN
        my=my+1
        month=0
END IF        
END DO
close(6)

END
