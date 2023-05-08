!
! $Header$
!
MODULE limit_read_mod
!
! This module reads the fichier "limit.nc" containing fields for surface forcing.
!
! Module subroutines :
!  limit_read_frac    : call limit_read_tot and return the fractions
!  limit_read_rug_alb : return rugosity and albedo, if coupled ocean call limit_read_tot first
!  limit_read_sst     : return sea ice temperature   
!  limit_read_tot     : read limit.nc and store the fields in local modules variables
!
  IMPLICIT NONE

  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE, PRIVATE :: pctsrf
!$OMP THREADPRIVATE(pctsrf)
  REAL, ALLOCATABLE, DIMENSION(:),   SAVE, PRIVATE :: rugos
!$OMP THREADPRIVATE(rugos)
  REAL, ALLOCATABLE, DIMENSION(:),   SAVE, PRIVATE :: albedo
!$OMP THREADPRIVATE(albedo)  
  REAL, ALLOCATABLE, DIMENSION(:),   SAVE, PRIVATE :: sst
!$OMP THREADPRIVATE(sst)  
  LOGICAL,SAVE :: read_continents=.FALSE.
!$OMP THREADPRIVATE(read_continents) 

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Public subroutines :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE limit_read_frac(itime, dtime, jour, pctsrf_new, is_modified)
!
! This subroutine is called from "change_srf_frac" for case of 
! ocean=force or from ocean_slab_frac for ocean=slab.
! The fraction for all sub-surfaces at actual time step is returned.

    USE dimphy
    INCLUDE "indicesol.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN) :: itime   ! time step
    INTEGER, INTENT(IN) :: jour    ! current day
    REAL   , INTENT(IN) :: dtime   ! length of time step
  
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon,nbsrf), INTENT(OUT) :: pctsrf_new  ! sub surface fractions
    LOGICAL, INTENT(OUT)                     :: is_modified ! true if pctsrf is modified at this time step

! End declaration
!****************************************************************************************

! 1) Read file limit.nc
    CALL limit_read_tot(itime, dtime, jour, is_modified)

! 2) Return the fraction read in limit_read_tot
    pctsrf_new(:,:) = pctsrf(:,:)
    
  END SUBROUTINE limit_read_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE limit_read_rug_alb(itime, dtime, jour, &
       knon, knindex, &
       rugos_out, alb_out)
!
! This subroutine is called from surf_land_bucket. 
! The flag "ok_veget" must can not be true. If coupled run, "ocean=couple"
! then this routine will call limit_read_tot.
!
    USE dimphy
    USE surface_data

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN) :: itime                     ! numero du pas de temps courant
    INTEGER, INTENT(IN) :: jour                      ! jour a lire dans l'annee
    REAL   , INTENT(IN) :: dtime                     ! pas de temps de la physique (en s)
    INTEGER, INTENT(IN) :: knon                      ! nomber of points on compressed grid
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex  ! grid point number for compressed grid
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT) :: rugos_out
    REAL, DIMENSION(klon), INTENT(OUT) :: alb_out
    
! Local variables
!****************************************************************************************
    INTEGER :: i
    LOGICAL :: is_modified
!****************************************************************************************

    IF (type_ocean == 'couple') THEN
       ! limit.nc has not yet been read. Do it now!
       CALL limit_read_tot(itime, dtime, jour, is_modified)
    END IF

    DO i=1,knon
       rugos_out(i) = rugos(knindex(i))
       alb_out(i)  = albedo(knindex(i))
    END DO

  END SUBROUTINE limit_read_rug_alb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE limit_read_sst(knon, knindex, sst_out)
!
! This subroutine returns the sea surface temperature already read from limit.nc.
!
    USE dimphy, ONLY : klon

    INTEGER, INTENT(IN)                  :: knon     ! nomber of points on compressed grid
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex  ! grid point number for compressed grid
    REAL, DIMENSION(klon), INTENT(OUT)   :: sst_out

    INTEGER :: i

    DO i = 1, knon
       sst_out(i) = sst(knindex(i))
    END DO

  END SUBROUTINE limit_read_sst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Private subroutine :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE limit_read_tot(itime, dtime, jour, is_modified)
!
! Read everything needed from limit.nc
!
! 0) Initialize
! 1) Open the file limit.nc, if it is time
! 2) Read fraction, if not type_ocean=couple
! 3) Read sea surface temperature, if not type_ocean=couple
! 4) Read albedo and rugosity for land surface, only in case of no vegetation model
! 5) Close file and distribuate variables to all processus

    USE dimphy
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
    USE surface_data, ONLY : type_ocean, ok_veget
    USE netcdf

    IMPLICIT NONE
    
    INCLUDE "indicesol.h"

! In- and ouput arguments
!****************************************************************************************
    INTEGER, INTENT(IN) :: itime   ! numero du pas de temps courant
    INTEGER, INTENT(IN) :: jour    ! jour a lire dans l'annee
    REAL   , INTENT(IN) :: dtime   ! pas de temps de la physique (en s)

    LOGICAL, INTENT(OUT) :: is_modified  ! true if pctsrf is modified at this time step

! Locals variables with attribute SAVE
!****************************************************************************************
! frequence de lecture des conditions limites (en pas de physique) 
    INTEGER,SAVE                              :: lmt_pas
!$OMP THREADPRIVATE(lmt_pas) 
    LOGICAL, SAVE                             :: first_call=.TRUE.
!$OMP THREADPRIVATE(first_call)    
! Locals variables
!****************************************************************************************
    INTEGER                                   :: nid, nvarid
    INTEGER                                   :: ii, ierr
    INTEGER, DIMENSION(2)                     :: start, epais
    REAL, DIMENSION(klon_glo,nbsrf)           :: pct_glo  ! fraction at global grid
    REAL, DIMENSION(klon_glo)                 :: sst_glo  ! sea-surface temperature at global grid
    REAL, DIMENSION(klon_glo)                 :: rug_glo  ! rugosity at global grid
    REAL, DIMENSION(klon_glo)                 :: alb_glo  ! albedo at global grid
    CHARACTER(len=20)                         :: modname='limit_read_mod'     

! End declaration
!****************************************************************************************

!****************************************************************************************
! 0) Initialization
!
!****************************************************************************************
    IF (first_call) THEN
       ! calculate number of time steps for one day
       lmt_pas = NINT(86400./dtime * 1.0)
       
       ! Allocate module save variables
       IF ( type_ocean /= 'couple' ) THEN
          ALLOCATE(pctsrf(klon,nbsrf), sst(klon), stat=ierr)
          IF (ierr /= 0) CALL abort_gcm(modname, 'PB in allocating pctsrf and sst',1)
       END IF

       IF ( .NOT. ok_veget ) THEN
          ALLOCATE(rugos(klon), albedo(klon), stat=ierr)
          IF (ierr /= 0) CALL abort_gcm(modname, 'PB in allocating rugos and albedo',1)
       END IF

       first_call=.FALSE.
    ENDIF
  
!****************************************************************************************
! 1) Open the file limit.nc if it is the right moment to read, once a day.
!    The file is read only by the master thread of the master mpi process(is_mpi_root)
!
!****************************************************************************************

    is_modified = .FALSE.
    IF (MOD(itime-1, lmt_pas) == 0) THEN   ! time to read
       is_modified = .TRUE.
!$OMP MASTER  ! Only master thread
       IF (is_mpi_root) THEN ! Only master processus

          ierr = NF90_OPEN ('limit.nc', NF90_NOWRITE, nid)
          IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,&
               'Pb d''ouverture du fichier de conditions aux limites',1)
          
          ! La tranche de donnees a lire:
          start(1) = 1
          start(2) = jour
          epais(1) = klon_glo
          epais(2) = 1


!****************************************************************************************
! 2) Read fraction if not type_ocean=couple
!
!****************************************************************************************

          IF ( type_ocean /= 'couple') THEN
!
! Ocean fraction
             ierr = NF90_INQ_VARID(nid, 'FOCE', nvarid)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname, 'Le champ <FOCE> est absent',1)
             
             ierr = NF90_GET_VAR(nid,nvarid,pct_glo(:,is_oce),start,epais)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <FOCE>' ,1)
!
! Sea-ice fraction
             ierr = NF90_INQ_VARID(nid, 'FSIC', nvarid)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <FSIC> est absent',1)

             ierr = NF90_GET_VAR(nid,nvarid,pct_glo(:,is_sic),start,epais)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <FSIC>' ,1)


! Read land and continentals fraction only if asked for
             IF (read_continents .OR. itime == 1) THEN
!
! Land fraction
                ierr = NF90_INQ_VARID(nid, 'FTER', nvarid)
                IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <FTER> est absent',1)
                
                ierr = NF90_GET_VAR(nid,nvarid,pct_glo(:,is_ter),start,epais)
                IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <FTER>',1)
!
! Continentale ice fraction
                ierr = NF90_INQ_VARID(nid, 'FLIC', nvarid)
                IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <FLIC> est absent',1)

                ierr = NF90_GET_VAR(nid,nvarid,pct_glo(:,is_lic),start,epais)
                IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <FLIC>',1)
             END IF

          END IF ! type_ocean /= couple

!****************************************************************************************
! 3) Read sea-surface temperature, if not coupled ocean
!
!****************************************************************************************
          IF ( type_ocean /= 'couple') THEN

             ierr = NF90_INQ_VARID(nid, 'SST', nvarid)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <SST> est absent',1)

             ierr = NF90_GET_VAR(nid,nvarid,sst_glo,start,epais)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <SST>',1)
          
          END IF

!****************************************************************************************
! 4) Read albedo and rugosity for land surface, only in case of no vegetation model
!
!****************************************************************************************

          IF (.NOT. ok_veget) THEN
!
! Read albedo
             ierr = NF90_INQ_VARID(nid, 'ALB', nvarid)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <ALB> est absent',1)

             ierr = NF90_GET_VAR(nid,nvarid,alb_glo,start,epais)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <ALB>',1)
!
! Read rugosity
             ierr = NF90_INQ_VARID(nid, 'RUG', nvarid)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Le champ <RUG> est absent',1)

             ierr = NF90_GET_VAR(nid,nvarid,rug_glo,start,epais)
             IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Lecture echouee pour <RUG>',1)

          END IF

!****************************************************************************************
! 5) Close file and distribuate variables to all processus
!
!****************************************************************************************
          ierr = NF90_CLOSE(nid)
          IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Pb when closing file', 1)
       ENDIF ! is_mpi_root

!$OMP END MASTER
!$OMP BARRIER

       IF ( type_ocean /= 'couple') THEN
          CALL Scatter(sst_glo,sst)
          CALL Scatter(pct_glo(:,is_oce),pctsrf(:,is_oce))
          CALL Scatter(pct_glo(:,is_sic),pctsrf(:,is_sic))
          IF (read_continents .OR. itime == 1) THEN
             CALL Scatter(pct_glo(:,is_ter),pctsrf(:,is_ter))
             CALL Scatter(pct_glo(:,is_lic),pctsrf(:,is_lic))
          END IF
       END IF

       IF (.NOT. ok_veget) THEN
          CALL Scatter(alb_glo, albedo)
          CALL Scatter(rug_glo, rugos)
       END IF

    ENDIF ! time to read

  END SUBROUTINE limit_read_tot


END MODULE limit_read_mod
