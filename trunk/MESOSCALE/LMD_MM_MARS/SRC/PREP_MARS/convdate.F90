PROGRAM convdate

IMPLICIT NONE

!
! compile: pgf90 convdate.F90 -o convdate 
! use: echo 352 | convdate
!


INTEGER :: gcm_day

INTEGER :: wrf_month,wrf_day

INTEGER :: init,i
INTEGER, PARAMETER :: MONTHS_PER_YEAR = 12
INTEGER, PARAMETER :: mday(MONTHS_PER_YEAR)   &
           = (/61,66,66,65,60,54,50,46,47,47,51,56/)

print *, 'which day ?'
read *, gcm_day

PRINT *,'*****************'
print *, 'GCM sol is ... (starting 0)', gcm_day
PRINT *,'*****************'

!
! Find WRF month and day
!

IF (gcm_day >= 669) THEN    !! gcm_day commence au jour 0
        PRINT *,'out of bounds ! martian year is 669 sols !'
        gcm_day=MOD(gcm_day,669)
        !STOP
ENDIF


init=gcm_day+1    !!+1 sinon on decale tout
DO i=1,MONTHS_PER_YEAR
        wrf_month=i
        init=init-mday(i)
        IF (init <= 0) EXIT
END DO

PRINT *,'corresponding WRF month is ',wrf_month
wrf_day=init+mday(wrf_month)
PRINT *,'corresponding WRF day is ',wrf_day
PRINT *,'*****************'


END PROGRAM convdate
