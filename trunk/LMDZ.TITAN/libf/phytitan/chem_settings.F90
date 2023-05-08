SUBROUTINE chem_settings(nid,ngrid,nlayer,indextime)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Author : Jan Vatant d'Ollone (2018)
! ------
!
! Purpose : + Intialise upper atmosphere chemistry pressure grid
! -------   and composition fields.
!           + This subroutine is called in phyetat0 and reads 
!           from a NetCDF "startfi.nc" file.
!           + The presence of pressure grid is compulsory in the
!           file but not composition fields. The presence of the
!           1st field is tested and then we assume there's either
!           no one or all of the 44 chemistry scheme species.
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE comchem_h
USE iostart, only: get_field, get_var, inquire_dimension_length
USE callkeys_mod, only : callchim

USE vertical_layers_mod, only: presnivs, pseudoalt

IMPLICIT NONE

!======================================================================
!  INPUTS :
!  --------
  INTEGER,INTENT(IN)  :: nid		! Input Netcdf file ID 
  INTEGER,INTENT(IN)  :: ngrid		! # of horizontal grid points
  INTEGER,INTENT(IN)  :: nlayer 	! # of vertical layers
  INTEGER,INTENT(IN)  :: indextime	! position on time axis
!======================================================================
! local variables:
  REAL :: phi0, phi

  INTEGER :: ierr	! status (returned by NetCDF functions)
  INTEGER :: nvarid	! ID of NetCDF variable
  INTEGER :: dimid	! ID of NetCDF dimension
  
  INTEGER :: ilay, iq
  
  LOGICAL ::  found
  
  CHARACTER(LEN=13) :: modname
!======================================================================

  ! 0. Start by reading how many layers of upper chemistry there are
  
  nlaykim_up=inquire_dimension_length("upper_chemistry_layers")
  
  ! 1. Allocates arrays in comchem_h
  
  CALL ini_comchem_h(ngrid)

  ! 2. Load upper chemistry pressure grid
  
  CALL get_var("preskim",preskim,found)
  IF (.NOT.found) THEN
    CALL abort_physic(modname,"Failed loading <preskim>",1)
  ENDIF
  WRITE(*,*) "chem_settings: Upper chemistry pressure grid <preskim> range:", &
               maxval(preskim), minval(preskim) 
               
  ! 3. Compute others chemistry grid
  
  ! a. Total pressure grid (0->1300km)
  DO ilay=1,nlayer ! GCM levels
    preskim_tot(ilay) = presnivs(ilay)
  ENDDO
  DO ilay=1,nlaykim_up ! Upper chemistry
    preskim_tot(ilay+nlayer) = preskim(ilay)
  ENDDO
  
  ! b. Pseudo-altitudes ( TBD - hydrostatic equilibrium or read somewhere ?)          
 
  ! 4. Inquire ( and load ) upper chemistry composition fields
 
  CALL get_field("H_up",ykim_up(1,:,:),found,indextime)
  IF (.NOT.found) THEN
  
    ! We assume we can't do uncomplete chemistry
    WRITE(*,*) "chem_settings: No upper chemistry fields."
    
    IF ( callchim ) THEN ! if chemistry we must have the upper fields !
      CALL abort_physic(modname,"Failed loading uppper chemistry fields, whereas callchim set to true !",1)
    ENDIF
    
    DEALLOCATE(ykim_up) ! it will be useless
    DEALLOCATE(ykim_tot) ! it will be useless
    
  ELSE
  
    WRITE(*,*) "chem_settings: H in upper atmosphere <H_up> range:", &
               minval(ykim_up(1,:,:)), maxval(ykim_up(1,:,:))
               
    ! Load others fields if first one found only as we assume we can't do uncomplete chemistry
    ! NB : We assume a given order of the 44 chemistry species !!
    ! ( H=1, H2=2 ..., C4N2=44) -> cf comchem_h
    
    DO iq=2,nkim
      CALL get_field(trim(cnames(iq))//"_up",ykim_up(iq,:,:),found,indextime)
      IF (.NOT.found) THEN
        CALL abort_physic(modname,"Failed loading <"//trim(cnames(iq))//"_up>",1)
      ENDIF
      WRITE(*,*) "chem_settings: "//trim(cnames(iq))//" in upper atmosphere <"//trim(cnames(iq))//"_up> range:", &
                 minval(ykim_up(iq,:,:)), maxval(ykim_up(iq,:,:))
    ENDDO
      
  ENDIF ! of if H_up found
 
END SUBROUTINE chem_settings
