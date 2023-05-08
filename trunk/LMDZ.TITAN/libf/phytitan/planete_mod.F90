MODULE planete_mod
  IMPLICIT NONE
  
  REAL,SAVE :: apoastr ! maximum star-planet distance (AU)
  REAL,SAVE :: periastr ! minimum star-planet distance (AU)
  REAL,SAVE :: year_day ! length of year (sols)
  REAL,SAVE :: peri_day ! date of periastron (sols since N. spring)
  REAL,SAVE :: obliquit ! Obliquity of the planet (deg)
!$OMP THREADPRIVATE(apoastr,periastr,year_day,peri_day,obliquit)
  REAL,SAVE :: nres ! tidal resonance ratio
  REAL,SAVE :: z0 ! surface roughness (m)
  REAL,SAVE :: lmixmin ! mixing length
  REAL,SAVE :: emin_turb ! minimal energy
!$OMP THREADPRIVATE(nres,z0,lmixmin,emin_turb)
  REAL,SAVE :: coefvis
  REAL,SAVE :: coefir
  REAL,SAVE :: timeperi
  REAL,SAVE :: e_elips
  REAL,SAVE :: p_elips
!$OMP THREADPRIVATE(coefvis,coefir,timeperi,e_elips,p_elips)
  
  REAL,SAVE :: preff ! reference surface pressure (Pa)	!read by master
  REAL,SAVE,ALLOCATABLE :: ap(:) ! hybrid coordinate at layer interface	!read by master
  REAL,SAVE,ALLOCATABLE :: bp(:) ! hybrid coordinate at layer interface 	!read by master
!$OMP THREADPRIVATE(preff,ap,bp)

  CONTAINS
  
  subroutine ini_planete_mod(nlayer,preff_dyn,ap_dyn,bp_dyn)
  
  implicit none
  integer,intent(in) :: nlayer ! number of atmospheric layers
  real,intent(in) :: preff_dyn ! reference surface pressure (Pa)
  real,intent(in) :: ap_dyn(nlayer+1) ! hybrid coordinate at interfaces
  real,intent(in) :: bp_dyn(nlayer+1) ! hybrid coordinate at interfaces
  
  allocate(ap(nlayer+1))
  allocate(bp(nlayer+1))
  
  preff=preff_dyn
  ap(:)=ap_dyn(:)
  bp(:)=bp_dyn(:)
  
  end subroutine ini_planete_mod
  
END MODULE planete_mod
