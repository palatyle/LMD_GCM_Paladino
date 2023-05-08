MODULE variables_mod

IMPLICIT NONE

REAL*8 :: JD_cur !pday !JD_cur ! Julian day
REAL*8 :: JH_cur_split !ptime !JH_cur_split ! Julian hour (fraction of day)
REAL*8 :: zdt_split !ptimestep !zdt_split ! time step over which the physics are evaluated
REAL*8 :: phour_ini ! start time (fraction of day) of the run 0=<phour_ini<1
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zplev_omp !pplev !zplev_omp(klon,llm+1) ! interlayer pressure (Pa)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zplay_omp !pplay !zplay_omp(klon,llm) ! mid-layer pressure (Pa)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zpk_omp!(klon,llm)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zphi_omp !pphi !zphi_omp(klon,llm) ! geopotential at midlayer
REAL*8,DIMENSION(:),ALLOCATABLE :: zphis_omp!(klon) ! surface geopotential
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zufi_omp !pu !zufi_omp(klon,llm) ! zonal wind (m/s)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zvfi_omp !pv !zvfi_omp(klon,llm) ! meridional wind (m/s)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zrfi_omp!(klon,llm) ! relative wind vorticity, in s-1
REAL*8,DIMENSION(:,:),ALLOCATABLE :: ztfi_omp !pt !ztfi_omp(klon,llm) ! temperature (K)
REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: zqfi_omp !pq !zqfi_omp(klon,llm,nqtot) ! tracers (*/kg of air)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: flxwfi_omp !flxw !flxwfi_omp(klon,llm) ! Vertical mass flux on lower mesh interfaces (kg/s) 
! tendencies (in */s) from the physics:
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zdufi_omp !pdu !zdufi_omp(klon,llm) ! tendency on zonal winds
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zdvfi_omp !pdv !zdvfi_omp(klon,llm) ! tendency on meridional winds
REAL*8,DIMENSION(:,:),ALLOCATABLE :: zdtfi_omp !pdt !zdtfi_omp(klon,llm) ! tendency on temperature
REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: zdqfi_omp !pdq !zdqfi_omp(klon,llm,nqtot) ! tendency on tracers
REAL*8,DIMENSION(:),ALLOCATABLE :: zdpsrf_omp !pdpsrf !zdpsrf_omp(klon) ! tendency on surface pressure

CONTAINS

SUBROUTINE allocate_interface(klon,llm,nqtot)

  IMPLICIT NONE
 
  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
 
  IF (.NOT.ALLOCATED(zplev_omp)) ALLOCATE(zplev_omp(klon,llm+1))
  IF (.NOT.ALLOCATED(zplay_omp)) ALLOCATE(zplay_omp(klon,llm))
  IF (.NOT.ALLOCATED(zpk_omp)) ALLOCATE(zpk_omp(klon,llm))
  IF (.NOT.ALLOCATED(zphi_omp)) ALLOCATE(zphi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zphis_omp)) ALLOCATE(zphis_omp(klon))
  IF (.NOT.ALLOCATED(zufi_omp)) ALLOCATE(zufi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zvfi_omp)) ALLOCATE(zvfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zrfi_omp)) ALLOCATE(zrfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(ztfi_omp)) ALLOCATE(ztfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zqfi_omp)) ALLOCATE(zqfi_omp(klon,llm,nqtot))
  IF (.NOT.ALLOCATED(flxwfi_omp)) ALLOCATE(flxwfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zdufi_omp)) ALLOCATE(zdufi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zdvfi_omp)) ALLOCATE(zdvfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zdtfi_omp)) ALLOCATE(zdtfi_omp(klon,llm))
  IF (.NOT.ALLOCATED(zdqfi_omp)) ALLOCATE(zdqfi_omp(klon,llm,nqtot))
  IF (.NOT.ALLOCATED(zdpsrf_omp)) ALLOCATE(zdpsrf_omp(klon))

END SUBROUTINE allocate_interface

SUBROUTINE deallocate_interface

  DEALLOCATE(zplev_omp)
  DEALLOCATE(zplay_omp)
  DEALLOCATE(zpk_omp)
  DEALLOCATE(zphi_omp)
  DEALLOCATE(zphis_omp)
  DEALLOCATE(zufi_omp)
  DEALLOCATE(zvfi_omp)
  DEALLOCATE(zrfi_omp)
  DEALLOCATE(ztfi_omp)
  DEALLOCATE(zqfi_omp)
  DEALLOCATE(flxwfi_omp)
  DEALLOCATE(zdufi_omp)
  DEALLOCATE(zdvfi_omp)
  DEALLOCATE(zdtfi_omp)
  DEALLOCATE(zdqfi_omp)
  DEALLOCATE(zdpsrf_omp)

END SUBROUTINE deallocate_interface

END MODULE variables_mod
