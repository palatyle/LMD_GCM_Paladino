subroutine cirs_haze(press,wno,taeros,taeroscat,cbar)
IMPLICIT NONE

real,intent(in)   :: press,wno
real,intent(inout):: taeros,taeroscat,cbar

!---------------------------
! Attention: taeros est en m-1
! mais la valeur semble devoir etre multipliee par lambda (en m) pour
! obtenir une valeur comparable a celle fournies par Sandrine...
! a tirer au clair...
!---------------------------

real         :: taerosold
logical,save :: firstcall=.true.

if (firstcall) then
   print*,"CIRS HAZE"
   firstcall=.false.
endif

if (wno.eq.600.) then
 print*,press,wno,taeros,taeroscat,cbar
endif

taerosold = taeros

! modif de taeros

! taeroscat est modifie proportionnellement a taeros

if (taerosold.ne.0.) then
   taeroscat = taeroscat/taerosold*taeros
endif

! Je maintiens le cbar du calcul microphysique

end subroutine cirs_haze
