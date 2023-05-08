! $Header$

SUBROUTINE limit_slab(itime, dtime, jour, lmt_bils, lmt_foce, diff_sst)

  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE netcdf 

  IMPLICIT NONE

  INCLUDE "indicesol.h"
  INCLUDE "temps.h"
  INCLUDE "clesphys.h"
  INCLUDE "dimensions.h"

! In- and ouput arguments
!****************************************************************************************
  INTEGER, INTENT(IN) :: itime   ! numero du pas de temps courant
  INTEGER, INTENT(IN) :: jour    ! jour a lire dans l'annee
  REAL   , INTENT(IN) :: dtime   ! pas de temps de la physique (en s)
  REAL, DIMENSION(klon), INTENT(OUT) :: lmt_bils, lmt_foce, diff_sst

! Locals variables with attribute SAVE
!****************************************************************************************
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bils_save, foce_save
!$OMP THREADPRIVATE(bils_save, foce_save)

! Locals variables
!****************************************************************************************
  INTEGER                  :: lmt_pas   
  INTEGER                  :: nvarid, nid, ierr, i
  INTEGER, DIMENSION(2)    :: start, epais 
  REAL, DIMENSION(klon_glo):: bils_glo, foce_glo, sst_l_glo, sst_lp1_glo, diff_sst_glo
  CHARACTER (len = 20)     :: modname = 'limit_slab'

! End declaration
!****************************************************************************************

  ! calculate number of time steps for one day
  lmt_pas = NINT(86400./dtime)
  
  IF (MOD(itime-1, lmt_pas) == 0) THEN   ! time to read
     !$OMP MASTER  ! Only master thread
     IF (is_mpi_root) THEN ! Only master processus
        print*,'in limit_slab time to read, itime=',itime
        
        ierr = NF90_OPEN ('limit_slab.nc', NF90_NOWRITE, nid)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,&
             'Pb in opening file limit_slab.nc',1)
        
        ! La tranche de donnees a lire:
        start(1) = 1
        start(2) = jour
        epais(1) = klon_glo
        epais(2) = 1

!****************************************************************************************
! 2) Read bils and ocean fraction
!
!****************************************************************************************
!
! Read bils_glo
        ierr = NF90_INQ_VARID(nid, 'BILS_OCE', nvarid)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'The variable <BILS_OCE> is abstent',1)

        ierr = NF90_GET_VAR(nid,nvarid,bils_glo,start,epais)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Reading of <BILS_OCE> failed',1)
!
! Read foce_glo
        ierr = NF90_INQ_VARID(nid, 'FOCE', nvarid)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'The variable <FOCE> is abstent',1)

        ierr = NF90_GET_VAR(nid,nvarid,foce_glo,start,epais)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Reading of <FOCE> failed',1)
!
! Read sst_glo for this day
        ierr = NF90_INQ_VARID(nid, 'SST', nvarid)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'The variable <SST> is abstent',1)

        ierr = NF90_GET_VAR(nid,nvarid,sst_l_glo,start,epais)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Reading of <SST> failed',1)

! Read sst_glo for one day ahead
        start(2) = jour + 1
        IF (start(2) > 360) start(2)=1
        ierr = NF90_GET_VAR(nid,nvarid,sst_lp1_glo,start,epais)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Reading of <SST> day+1 failed',1)

! Calculate difference in temperature between this day and one ahead
        DO i=1, klon_glo-1
           diff_sst_glo(i) = sst_lp1_glo(i) - sst_l_glo(i)
        END DO
        diff_sst_glo(klon_glo) = sst_lp1_glo(klon_glo) - sst_l_glo(1)

!****************************************************************************************
! 5) Close file and distribuate variables to all processus
!
!****************************************************************************************
        ierr = NF90_CLOSE(nid)
        IF (ierr /= NF90_NOERR) CALL abort_gcm(modname,'Pb when closing file', 1)
     ENDIF ! is_mpi_root

!$OMP END MASTER
       
     IF (.NOT. ALLOCATED(bils_save)) THEN
        ALLOCATE(bils_save(klon), foce_save(klon), stat=ierr)
        IF (ierr /= 0) CALL abort_gcm('limit_slab', 'pb in allocation',1)
     END IF

     CALL Scatter(bils_glo, bils_save)
     CALL Scatter(foce_glo, foce_save)
     CALL Scatter(diff_sst_glo, diff_sst)
     
  ELSE ! not time to read
     diff_sst(:) = 0.
  ENDIF ! time to read

  lmt_bils(:) = bils_save(:)
  lmt_foce(:) = foce_save(:)
  
END SUBROUTINE limit_slab
