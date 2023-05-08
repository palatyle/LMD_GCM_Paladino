PROGRAM convdate

IMPLICIT NONE

!
! pgf90 time.F makecal.F90
!

INTEGER :: gcm_day, init_oui

INTEGER :: wrf_month,wrf_day

INTEGER :: init,i
INTEGER, PARAMETER :: MONTHS_PER_YEAR = 12
INTEGER, PARAMETER :: mday(MONTHS_PER_YEAR)   &
           = (/61,66,66,65,60,54,50,46,47,47,51,56/)
REAL :: ls

CHARACTER(2) :: cday,cmonth

!!PRINT *, '     init --','   GCM sol --', '   GCM ls --', '     WRF month --', '   WRF day --'
!!PRINT *, '     init --','   GCM sol --', '   GCM ls --', '     WRF'
PRINT *, '-- sol -- ', ' -- ls -- ', ' -- month -- ', ' -- day --'

DO gcm_day=0,668

init_oui=0

init=gcm_day+1    !!+1 sinon on decale tout
DO i=1,MONTHS_PER_YEAR
        wrf_month=i
        init=init-mday(i)
        IF (init <= 0) EXIT
END DO

wrf_day=init+mday(wrf_month)

IF (MODULO(gcm_day+1,3) == 0) init_oui=99

CALL sol2ls(real(gcm_day),ls)

Write( cday, '(i0)' )  wrf_day
Write( cmonth, '(i0)' )  wrf_month
IF (wrf_day < 10) cday='0'//cday
IF (wrf_month < 10) cmonth='0'//cmonth

IF (ls < 1e-3) ls=0.

!!PRINT *, init_oui, gcm_day, ls, wrf_day, wrf_month, '2024'//'-'//cmonth//'-'//cday//'_'//'00'//':00:00'
!!PRINT *, init_oui, gcm_day, ls, '2024'//'-'//cmonth//'-'//cday//'_'//'00'//':00:00'
PRINT *, gcm_day, ls, cmonth//'  '//cday

ENDDO

END PROGRAM convdate
