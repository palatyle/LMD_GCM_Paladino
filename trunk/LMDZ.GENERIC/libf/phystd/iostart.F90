MODULE iostart

    IMPLICIT NONE
PRIVATE
    INTEGER,SAVE :: nid_start ! NetCDF file identifier for startfi.nc file
    INTEGER,SAVE :: nid_restart ! NetCDF file identifier for restartfi.nc file
!$OMP THREADPRIVATE(nid_start,nid_restart)
    
    ! restartfi.nc file dimension identifiers: (see open_restartphy())
    INTEGER,SAVE :: idim1 ! "index" dimension
    INTEGER,SAVE :: idim2 ! "physical_points" dimension
    INTEGER,SAVE :: idim3 ! "subsurface_layers" dimension
    INTEGER,SAVE :: idim4 ! "nlayer_plus_1" dimension
    INTEGER,SAVE :: idim5 ! "number_of_advected_fields" dimension
    INTEGER,SAVE :: idim6 ! "nlayer" dimension
    INTEGER,SAVE :: idim7 ! "Time" dimension
    INTEGER,SAVE :: idim8 ! "ocean_layers" dimension
    INTEGER,SAVE :: timeindex ! current time index (for time-dependent fields)
!$OMP THREADPRIVATE(idim1,idim2,idim3,idim4,idim5,idim6,idim7,timeindex)
    INTEGER,PARAMETER :: length=100 ! size of tab_cntrl array
    
    INTERFACE get_field
      MODULE PROCEDURE Get_field_r1,Get_field_r2,Get_field_r3
    END INTERFACE get_field
    
    INTERFACE get_var
      MODULE PROCEDURE get_var_r0,Get_var_r1,Get_var_r2,Get_var_r3
    END INTERFACE get_var

    INTERFACE put_field
      MODULE PROCEDURE put_field_r1,put_field_r2,put_field_r3
    END INTERFACE put_field

    INTERFACE put_var
      MODULE PROCEDURE put_var_r0,put_var_r1,put_var_r2,put_var_r3
    END INTERFACE put_var

    PUBLIC nid_start, length
    PUBLIC get_field,get_var,put_field,put_var
    PUBLIC inquire_dimension, inquire_dimension_length
    PUBLIC inquire_field, inquire_field_ndims
    PUBLIC open_startphy,close_startphy,open_restartphy,close_restartphy
    
CONTAINS

  SUBROUTINE open_startphy(filename)
  USE netcdf, only: NF90_OPEN, NF90_NOERR, NF90_NOWRITE, nf90_strerror
  USE mod_phys_lmdz_para, only: is_master, bcast
  IMPLICIT NONE
    CHARACTER(LEN=*) :: filename
    INTEGER          :: ierr

    IF (is_master) THEN
      ierr = NF90_OPEN (filename, NF90_NOWRITE, nid_start)
      IF (ierr.NE.NF90_NOERR) THEN
        write(*,*)'open_startphy: problem opening file '//trim(filename)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
    ENDIF
    
    CALL bcast(nid_start) ! tell all procs about nid_start
  
  END SUBROUTINE open_startphy

  SUBROUTINE close_startphy
  USE netcdf, only: NF90_CLOSE
  USE mod_phys_lmdz_para, only: is_master
  IMPLICIT NONE
    INTEGER          :: ierr

    IF (is_master) THEN
        ierr = NF90_CLOSE (nid_start)
    ENDIF

  END SUBROUTINE close_startphy


  FUNCTION inquire_field(Field_name)
  ! check if a given field is present in the input file
  USE netcdf, only: NF90_INQ_VARID, NF90_NOERR
  USE mod_phys_lmdz_para, only: is_master, bcast
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Field_name
    LOGICAL :: inquire_field
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_master) THEN
      ierr=NF90_INQ_VARID(nid_start,Field_name,varid)
      IF (ierr==NF90_NOERR) THEN
        Inquire_field=.TRUE.
      ELSE
        Inquire_field=.FALSE.
      ENDIF
    ENDIF

    CALL bcast(inquire_field)

  END FUNCTION inquire_field


  FUNCTION inquire_field_ndims(Field_name)
  ! give the number of dimensions of "Field_name" stored in the input file 
  USE netcdf, only: nf90_inq_varid, nf90_inquire_variable, &
                    NF90_NOERR, nf90_strerror
  USE mod_phys_lmdz_para, only: is_master, bcast
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Field_name
    INTEGER :: inquire_field_ndims
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_master) THEN
      ierr=nf90_inq_varid(nid_start,Field_name,varid)
      ierr=nf90_inquire_variable(nid_start,varid,&
                                  ndims=inquire_field_ndims)
      IF (ierr.NE.NF90_NOERR) THEN
        write(*,*)'inquire_field_ndims: problem obtaining ndims of '&
                  //trim(field_name)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
    ENDIF

    CALL bcast(inquire_field_ndims)

  END FUNCTION inquire_field_ndims


  FUNCTION inquire_dimension(Field_name)
  ! check if a given dimension is present in the input file
  USE netcdf, only: nf90_inq_dimid, NF90_NOERR
  USE mod_phys_lmdz_para, only: is_master, bcast
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Field_name
    LOGICAL :: inquire_dimension
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_master) THEN
      ierr=NF90_INQ_DIMID(nid_start,Field_name,varid)
      IF (ierr==NF90_NOERR) THEN
        Inquire_dimension=.TRUE.
      ELSE
        Inquire_dimension=.FALSE.
      ENDIF
    ENDIF

    CALL bcast(inquire_dimension)

  END FUNCTION inquire_dimension

  FUNCTION inquire_dimension_length(Field_name)
  ! give the length of the "Field_name" dimension stored in the input file
  USE netcdf, only: nf90_inquire_dimension, nf90_inq_dimid, &
                    NF90_NOERR, nf90_strerror
  USE mod_phys_lmdz_para, only: is_master, bcast
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Field_name
    INTEGER :: inquire_dimension_length
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_master) THEN
      ierr=nf90_inq_dimid(nid_start,Field_name,varid)
      ierr=nf90_inquire_dimension(nid_start,varid,&
                                  len=inquire_dimension_length)
      IF (ierr.NE.NF90_NOERR) THEN
        write(*,*)'inquire_field_length: problem obtaining length of '&
                  //trim(field_name)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
    ENDIF

    CALL bcast(inquire_dimension_length)

  END FUNCTION inquire_dimension_length



  SUBROUTINE Get_Field_r1(field_name,field,found,timeindex)
  ! For a surface field
  use mod_grid_phy_lmdz, only: klon_glo ! number of atmospheric columns (full grid)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 
    INTEGER,INTENT(IN),OPTIONAL    :: timeindex ! time index of sought data

    integer :: edges(4), corners(4)

    edges(1)=klon_glo
    edges(2:4)=1
    corners(1:4)=1
    if (PRESENT(timeindex)) then
      corners(2)=timeindex
    endif

    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,1,corners,edges,found)
    ELSE
      CALL Get_field_rgen(field_name,field,1,corners,edges)
    ENDIF
      
  END SUBROUTINE Get_Field_r1
  
  SUBROUTINE Get_Field_r2(field_name,field,found,timeindex)
  ! For a "3D" horizontal-vertical field
  use mod_grid_phy_lmdz, only: klon_glo ! number of atmospheric columns (full grid)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:,:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 
    INTEGER,INTENT(IN),OPTIONAL    :: timeindex ! time index of sought data

    integer :: edges(4), corners(4)

    edges(1)=klon_glo
    edges(2)=size(field,2)
    edges(3:4)=1
    corners(1:4)=1
    if (PRESENT(timeindex)) then
      corners(3)=timeindex
    endif
    
    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,size(field,2),&
                          corners,edges,found)
    ELSE
      CALL Get_field_rgen(field_name,field,size(field,2),&
                          corners,edges)
    ENDIF

      
  END SUBROUTINE Get_Field_r2
  
  SUBROUTINE Get_Field_r3(field_name,field,found,timeindex)
  ! for a "4D" field surf/alt/??
  use mod_grid_phy_lmdz, only: klon_glo ! number of atmospheric columns (full grid)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:,:,:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 
    INTEGER,INTENT(IN),OPTIONAL    :: timeindex ! time index of sought data

    integer :: edges(4), corners(4)

    edges(1)=klon_glo
    edges(2)=size(field,2)
    edges(3)=size(field,3)
    edges(4)=1
    corners(1:4)=1
    if (PRESENT(timeindex)) then
      corners(4)=timeindex
    endif
    
    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,size(field,2)*size(field,3),&
                          corners,edges,found)
    ELSE
      CALL Get_field_rgen(field_name,field,size(field,2)*size(field,3),&
                          corners,edges)
    ENDIF
      
  END SUBROUTINE Get_Field_r3
  
  SUBROUTINE Get_field_rgen(field_name,field,field_size, &
                            corners,edges,found)
  USE netcdf
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: Field_name
    INTEGER          :: field_size
    REAL             :: field(klon,field_size)
    INTEGER,INTENT(IN) :: corners(4)
    INTEGER,INTENT(IN) :: edges(4)
    LOGICAL,OPTIONAL :: found
    
    REAL    :: field_glo(klon_glo,field_size)
    LOGICAL :: tmp_found
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_master) THEN
  
      ierr=NF90_INQ_VARID(nid_start,Field_name,varid)
      
      IF (ierr==NF90_NOERR) THEN
        CALL body(field_glo)
        tmp_found=.TRUE.
      ELSE
        tmp_found=.FALSE.
      ENDIF
    
    ENDIF
    
    CALL bcast(tmp_found)

    IF (tmp_found) THEN
      CALL scatter(field_glo,field)
    ENDIF
    
    IF (PRESENT(found)) THEN
      found=tmp_found
    ELSE
      IF (.NOT. tmp_found) THEN
        PRINT*, 'get_field_rgen: Field <'//field_name//'> not found'
        CALL abort
      ENDIF
    ENDIF
 
    
    CONTAINS
     
     SUBROUTINE body(field_glo)
       REAL :: field_glo(klon_glo*field_size)
         ierr=NF90_GET_VAR(nid_start,varid,field_glo,corners,edges)
         IF (ierr/=NF90_NOERR) THEN
           ! La variable exist dans le fichier mais la lecture a echouee. 
           PRINT*, 'get_field_rgen: Failed reading <'//field_name//'>'

!           IF (field_name=='CLWCON' .OR. field_name=='RNEBCON' .OR. field_name=='RATQS') THEN
!              ! Essaye de lire le variable sur surface uniqument, comme fait avant
!              field_glo(:)=0.
!              ierr=NF90_GET_VAR(nid_start,varid,field_glo(1:klon_glo))
!              IF (ierr/=NF90_NOERR) THEN
!                 PRINT*, 'phyetat0: Lecture echouee aussi en 2D pour <'//field_name//'>'
!                 CALL abort
!              ELSE
!                 PRINT*, 'phyetat0: La variable <'//field_name//'> lu sur surface seulement'!, selon ancien format, le reste mis a zero'
!              END IF
!           ELSE
              CALL abort
!           ENDIF
         ENDIF

     END SUBROUTINE body

  END SUBROUTINE Get_field_rgen


  SUBROUTINE get_var_r0(var_name,var,found)
  ! Get a scalar from input file
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var
    LOGICAL,OPTIONAL,INTENT(OUT) :: found

    REAL                         :: varout(1)
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,varout,size(varout),found)
    ELSE
      CALL Get_var_rgen(var_name,varout,size(varout))
    ENDIF
    var=varout(1)
 
  END SUBROUTINE get_var_r0

  SUBROUTINE get_var_r1(var_name,var,found)
  ! Get a vector from input file
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var(:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r1

  SUBROUTINE get_var_r2(var_name,var,found)
  ! Get a 2D field from input file
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(OUT)             :: var(:,:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r2

  SUBROUTINE get_var_r3(var_name,var,found)
  ! Get a 3D field frominput file
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var(:,:,:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r3

  SUBROUTINE Get_var_rgen(var_name,var,var_size,found)
  USE netcdf
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: var_name
    INTEGER          :: var_size
    REAL             :: var(var_size)
    LOGICAL,OPTIONAL :: found
    
    LOGICAL :: tmp_found
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_mpi_root .AND. is_omp_root) THEN
  
      ierr=NF90_INQ_VARID(nid_start,var_name,varid)
      
      IF (ierr==NF90_NOERR) THEN
        ierr=NF90_GET_VAR(nid_start,varid,var)
        IF (ierr/=NF90_NOERR) THEN
          PRINT*, 'phyetat0: Failed loading <'//trim(var_name)//'>'
          CALL abort
        ENDIF
        tmp_found=.TRUE.
      ELSE
        tmp_found=.FALSE.
      ENDIF
    
    ENDIF
    
    CALL bcast(tmp_found)

    IF (tmp_found) THEN
      CALL bcast(var)
    ENDIF
    
    IF (PRESENT(found)) THEN
      found=tmp_found
    ELSE
      IF (.NOT. tmp_found) THEN
        PRINT*, 'phyetat0: Variable <'//trim(var_name)//'> not found'
        CALL abort
      ENDIF
    ENDIF

  END SUBROUTINE Get_var_rgen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE open_restartphy(filename)
  USE netcdf, only: NF90_CREATE, NF90_CLOBBER, NF90_64BIT_OFFSET, &
                    NF90_NOERR, nf90_strerror, &
                    NF90_PUT_ATT, NF90_GLOBAL, NF90_DEF_DIM, &
                    NF90_UNLIMITED, NF90_ENDDEF, &
                    NF90_WRITE, NF90_OPEN
  USE mod_phys_lmdz_para, only: is_master
  USE mod_grid_phy_lmdz, only: klon_glo
  USE dimphy, only: klev, klevp1
  USE tracer_h, only: nqtot
  USE comsoil_h, only: nsoilmx
  USE slab_ice_h, only: noceanmx

  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER                     :: ierr
    LOGICAL,SAVE :: already_created=.false.
!$OMP THREADPRIVATE(already_created)
    
    IF (is_master) THEN
      if (.not.already_created) then
        ! At the very first call, create the file
        ierr=NF90_CREATE(filename,IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), &
                          nid_restart)
        IF (ierr/=NF90_NOERR) THEN
          write(*,*)'open_restartphy: problem creating file '//trim(filename)
          write(*,*)trim(nf90_strerror(ierr))
          CALL ABORT
        ENDIF
        already_created=.true.
      else
        ! Just open the file
        ierr=NF90_OPEN(filename,NF90_WRITE,nid_restart)
        IF (ierr/=NF90_NOERR) THEN
          write(*,*)'open_restartphy: problem opening file '//trim(filename)
          write(*,*)trim(nf90_strerror(ierr))
          CALL ABORT
        ENDIF
        return
      endif ! of if (.not.already_created)

      ierr=NF90_PUT_ATT(nid_restart,NF90_GLOBAL,"title",&
                        "Physics start file")
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem writing title '
        write(*,*)trim(nf90_strerror(ierr))
      ENDIF

      ierr=NF90_DEF_DIM(nid_restart,"index",length,idim1)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining index dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
      
      ierr=NF90_DEF_DIM(nid_restart,"physical_points",klon_glo,idim2)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining physical_points dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
      
      ierr=NF90_DEF_DIM(nid_restart,"subsurface_layers",nsoilmx,idim3)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining subsurface_layers dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
      
      ierr=NF90_DEF_DIM(nid_restart,"nlayer_plus_1",klevp1,idim4)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining nlayer_plus_1 dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
      
      if (nqtot>0) then
        ! only define a tracer dimension if there are tracers
        ierr=NF90_DEF_DIM(nid_restart,"number_of_advected_fields",nqtot,idim5)
        IF (ierr/=NF90_NOERR) THEN
          write(*,*)'open_restartphy: problem defining number_of_advected_fields dimension '
          write(*,*)trim(nf90_strerror(ierr))
          CALL ABORT
        ENDIF
      endif

      ierr=NF90_DEF_DIM(nid_restart,"nlayer",klev,idim6)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining nlayer dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
      
      ierr=NF90_DEF_DIM(nid_restart,"Time",NF90_UNLIMITED,idim7)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining Time dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF

      ierr=NF90_DEF_DIM(nid_restart,"ocean_layers",noceanmx,idim8)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem defining oceanic layer dimension '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF


      ierr=NF90_ENDDEF(nid_restart)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'open_restartphy: problem ending definition mode '
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
    ENDIF

  END SUBROUTINE open_restartphy

  SUBROUTINE close_restartphy
  USE netcdf, only: NF90_CLOSE
  USE mod_phys_lmdz_para, only: is_master
  IMPLICIT NONE
    INTEGER          :: ierr

    IF (is_master) THEN
      ierr = NF90_CLOSE (nid_restart)
    ENDIF
 
  END SUBROUTINE close_restartphy

  SUBROUTINE put_field_r1(field_name,title,field,time)
  ! For a surface field
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:)
  REAL,OPTIONAL,INTENT(IN)       :: time
  
  IF (present(time)) THEN
    ! if timeindex is present, it is a time-dependent variable
    CALL put_field_rgen(field_name,title,field,1,time)
  ELSE
    CALL put_field_rgen(field_name,title,field,1)
  ENDIF
  
  END SUBROUTINE put_field_r1

  SUBROUTINE put_field_r2(field_name,title,field,time)
  ! For a "3D" horizontal-vertical field
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:,:)
  REAL,OPTIONAL,INTENT(IN)       :: time
  
  IF (present(time)) THEN
    ! if timeindex is present, it is a time-dependent variable
    CALL put_field_rgen(field_name,title,field,size(field,2),time)
  ELSE
    CALL put_field_rgen(field_name,title,field,size(field,2))
  ENDIF
  
  END SUBROUTINE put_field_r2

  SUBROUTINE put_field_r3(field_name,title,field,time)
  ! For a "4D" field surf/alt/??
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:,:,:)
  REAL,OPTIONAL,INTENT(IN)       :: time
  
  IF (present(time)) THEN
    ! if timeindex is present, it is a time-dependent variable
    CALL put_field_rgen(field_name,title,field,size(field,2)*size(field,3),&
                        time)
  ELSE  
    CALL put_field_rgen(field_name,title,field,size(field,2)*size(field,3))
  ENDIF
  
  END SUBROUTINE put_field_r3
  
  SUBROUTINE put_field_rgen(field_name,title,field,field_size,time)
  USE netcdf
  USE dimphy
  USE comsoil_h, only: nsoilmx
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE slab_ice_h, only: noceanmx
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  INTEGER,INTENT(IN)             :: field_size
  REAL,INTENT(IN)                :: field(klon,field_size)
  REAL,OPTIONAL,INTENT(IN)       :: time
  
  REAL                           :: field_glo(klon_glo,field_size)
  INTEGER                        :: ierr
  INTEGER                        :: nvarid
  INTEGER                        :: idim
   
    CALL gather(field,field_glo)
    
    IF (is_master) THEN

      IF (field_size==1) THEN
        ! input is a 1D "surface field" array
        if (.not.present(time)) then ! for a time-independent field
          ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                            (/idim2/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                            (/idim2/),nvarid)
#endif
          if (ierr.ne.NF90_NOERR) then
            write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
            write(*,*)trim(nf90_strerror(ierr))
          endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr=NF90_ENDDEF(nid_restart)
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo)
        else
          ! check if the variable has already been defined:
          ierr=NF90_INQ_VARID(nid_restart,field_name,nvarid)
          if (ierr/=NF90_NOERR) then ! variable not found, define it
            ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                              (/idim2,idim7/),nvarid)
#else
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                              (/idim2,idim7/),nvarid)
#endif
            if (ierr.ne.NF90_NOERR) then
              write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
              write(*,*)trim(nf90_strerror(ierr))
            endif
            IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
            ierr=NF90_ENDDEF(nid_restart)
          endif
          ! Write the variable
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo,&
                            start=(/1,timeindex/))
        endif ! of if (.not.present(timeindex))

      ELSE IF (field_size==klev) THEN
        ! input is a 2D "atmospheric field" array
        if (.not.present(time)) then ! for a time-independent field
          ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                            (/idim2,idim6/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                            (/idim2,idim6/),nvarid)
#endif
          if (ierr.ne.NF90_NOERR) then
            write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
            write(*,*)trim(nf90_strerror(ierr))
          endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr=NF90_ENDDEF(nid_restart)
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo)
        else
          ! check if the variable has already been defined:
          ierr=NF90_INQ_VARID(nid_restart,field_name,nvarid)
          if (ierr/=NF90_NOERR) then ! variable not found, define it
            ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                              (/idim2,idim6,idim7/),nvarid)
#else
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                              (/idim2,idim6,idim7/),nvarid)
#endif
            if (ierr.ne.NF90_NOERR) then
              write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
              write(*,*)trim(nf90_strerror(ierr))
            endif
            IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
            ierr=NF90_ENDDEF(nid_restart)
          endif
          ! Write the variable
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo,&
                            start=(/1,1,timeindex/))
        endif ! of if (.not.present(time))

      ELSE IF (field_size==klevp1) THEN
        ! input is a 2D "interlayer atmospheric field" array
        if (.not.present(time)) then ! for a time-independent field
          ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                            (/idim2,idim4/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                            (/idim2,idim4/),nvarid)
#endif
          if (ierr.ne.NF90_NOERR) then
            write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
            write(*,*)trim(nf90_strerror(ierr))
          endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr = NF90_ENDDEF(nid_restart)
          ierr = NF90_PUT_VAR(nid_restart,nvarid,field_glo)
        else
          ! check if the variable has already been defined:
          ierr=NF90_INQ_VARID(nid_restart,field_name,nvarid)
          if (ierr/=NF90_NOERR) then ! variable not found, define it
            ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                              (/idim2,idim4,idim7/),nvarid)
#else
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                              (/idim2,idim4,idim7/),nvarid)
#endif
            if (ierr.ne.NF90_NOERR) then
              write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
              write(*,*)trim(nf90_strerror(ierr))
            endif
            IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
            ierr=NF90_ENDDEF(nid_restart)
          endif
          ! Write the variable
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo,&
                            start=(/1,1,timeindex/))
        endif ! of if (.not.present(timeindex))

      ELSE IF (field_size==nsoilmx) THEN
        ! input is a 2D "subsurface field" array
        if (.not.present(time)) then ! for a time-independent field
          ierr = NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                            (/idim2,idim3/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                            (/idim2,idim3/),nvarid)
#endif
          if (ierr.ne.NF90_NOERR) then
            write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
            write(*,*)trim(nf90_strerror(ierr))
          endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr = NF90_ENDDEF(nid_restart)
          ierr = NF90_PUT_VAR(nid_restart,nvarid,field_glo)
        else
          ! check if the variable has already been defined:
          ierr=NF90_INQ_VARID(nid_restart,field_name,nvarid)
          if (ierr/=NF90_NOERR) then ! variable not found, define it
            ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                              (/idim2,idim3,idim7/),nvarid)
#else
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                              (/idim2,idim3,idim7/),nvarid)
#endif
           if (ierr.ne.NF90_NOERR) then
              write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
              write(*,*)trim(nf90_strerror(ierr))
            endif
            IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
            ierr=NF90_ENDDEF(nid_restart)
          endif
          ! Write the variable
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo,&
                            start=(/1,1,timeindex/))

        endif ! of if (.not.present(time))

      ELSE IF (field_size==noceanmx) THEN
        ! input is a 2D "oceanic field" array
        if (.not.present(time)) then ! for a time-independent field
          ierr = NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                            (/idim2,idim8/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                            (/idim2,idim8/),nvarid)
#endif
          if (ierr.ne.NF90_NOERR) then
            write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
            write(*,*)trim(nf90_strerror(ierr))
          endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr = NF90_ENDDEF(nid_restart)
          ierr = NF90_PUT_VAR(nid_restart,nvarid,field_glo)
        else
          ! check if the variable has already been defined:
          ierr=NF90_INQ_VARID(nid_restart,field_name,nvarid)
          if (ierr/=NF90_NOERR) then ! variable not found, define it
            ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_DOUBLE,&
                              (/idim2,idim8,idim7/),nvarid)
#else
            ierr=NF90_DEF_VAR(nid_restart,field_name,NF90_FLOAT,&
                              (/idim2,idim8,idim7/),nvarid)
#endif
           if (ierr.ne.NF90_NOERR) then
              write(*,*)"put_field_rgen error: failed to define "//trim(field_name)
              write(*,*)trim(nf90_strerror(ierr))
            endif
            IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
            ierr=NF90_ENDDEF(nid_restart)
          endif
          ! Write the variable
          ierr=NF90_PUT_VAR(nid_restart,nvarid,field_glo,&
                            start=(/1,1,timeindex/))

        endif ! of if (.not.present(time))


      ELSE
        PRINT *, "Error phyredem(put_field_rgen) : wrong dimension for ",trim(field_name)
        write(*,*) "  field_size =",field_size
        CALL ABORT
      ENDIF

      ! Check the writting of field to file went OK
      if (ierr.ne.NF90_NOERR) then
        write(*,*) " Error phyredem(put_field_rgen) : failed writing ",trim(field_name)
        write(*,*)trim(nf90_strerror(ierr))
        call abort
      endif

    ENDIF ! of IF (is_master)
    
  END SUBROUTINE put_field_rgen  
  
  SUBROUTINE put_var_r0(var_name,title,var)
  ! Put a scalar in file
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var
     REAL                        :: varin(1)
     
     varin(1)=var
     
     CALL put_var_rgen(var_name,title,varin,size(varin))

  END SUBROUTINE put_var_r0


  SUBROUTINE put_var_r1(var_name,title,var)
  ! Put a vector in file
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r1
 
  SUBROUTINE put_var_r2(var_name,title,var)
  ! Put a 2D field in file
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:,:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r2     
  
  SUBROUTINE put_var_r3(var_name,title,var)
  ! Put a 3D field in file
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:,:,:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r3

  SUBROUTINE put_var_rgen(var_name,title,var,var_size)
  USE netcdf, only: NF90_REDEF, NF90_DEF_VAR, NF90_ENDDEF, NF90_PUT_VAR, &
                    NF90_FLOAT, NF90_DOUBLE, &
                    NF90_PUT_ATT, NF90_NOERR, nf90_strerror, &
                    nf90_inq_dimid, nf90_inquire_dimension, NF90_INQ_VARID
  USE comsoil_h, only: nsoilmx
  USE mod_phys_lmdz_para, only: is_master
  USE slab_ice_h, only: noceanmx
  IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     INTEGER,INTENT(IN)          :: var_size
     REAL,INTENT(IN)             :: var(var_size)
     
     INTEGER :: ierr
     INTEGER :: nvarid
     INTEGER :: idim1d
     logical,save :: firsttime=.true.
!$OMP THREADPRIVATE(firsttime)
         
    IF (is_master) THEN

      IF (var_name=="Time") THEN
        ! Very specific case of "Time" variable
        if (firsttime) then
          ! Create the "Time variable"
          ierr=NF90_REDEF(nid_restart)
#ifdef NC_DOUBLE
          ierr=NF90_DEF_VAR(nid_restart,var_name,NF90_DOUBLE,&
                            (/idim7/),nvarid)
#else
          ierr=NF90_DEF_VAR(nid_restart,var_name,NF90_FLOAT,&
                            (/idim7/),nvarid)
#endif
          IF (LEN_TRIM(title) > 0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
          ierr=NF90_ENDDEF(nid_restart)
          
          firsttime=.false.
        endif
        ! Append "time" value
        ! get current length of "Time" dimension
        ierr=nf90_inq_dimid(nid_restart,var_name,idim1d)
        ierr=nf90_inquire_dimension(nid_restart,idim1d,len=timeindex)
        timeindex=timeindex+1
        ierr=NF90_INQ_VARID(nid_restart,var_name,nvarid)
        ierr=NF90_PUT_VAR(nid_restart,nvarid,var,&
                           start=(/timeindex/))
        IF (ierr/=NF90_NOERR) THEN
          write(*,*)'put_var_rgen: problem writing Time'
          write(*,*)trim(nf90_strerror(ierr))
          CALL ABORT
        ENDIF
        return ! nothing left to do
      ELSEIF (var_size==length) THEN
        ! We know it is a "controle" kind of 1D array
        idim1d=idim1
      ELSEIF (var_size==nsoilmx) THEN
        ! We know it is an  "mlayer" kind of 1D array
        idim1d=idim3
      ELSEIF (var_size==noceanmx) THEN
        ! We know it is an  "mlayer" kind of 1D array
        idim1d=idim8
      ELSE 
        PRINT *, "put_var_rgen error : wrong dimension"
        write(*,*) "  var_size =",var_size
        CALL abort

      ENDIF ! of IF (var_size==length) THEN

      ! Swich to NetCDF define mode
      ierr=NF90_REDEF (nid_restart)
      ! Define the variable
#ifdef NC_DOUBLE
      ierr=NF90_DEF_VAR(nid_restart,var_name,NF90_DOUBLE,(/idim1d/),nvarid)
#else
      ierr=NF90_DEF_VAR(nid_restart,var_name,NF90_FLOAT,(/idim1d/),nvarid)
#endif
      ! Add a "title" attribute
      IF (LEN_TRIM(title)>0) ierr=NF90_PUT_ATT(nid_restart,nvarid,"title",title)
      ! Swich out of define mode
      ierr=NF90_ENDDEF(nid_restart)
      ! Write variable to file
      ierr=NF90_PUT_VAR(nid_restart,nvarid,var)
      IF (ierr/=NF90_NOERR) THEN
        write(*,*)'put_var_rgen: problem writing '//trim(var_name)
        write(*,*)trim(nf90_strerror(ierr))
        CALL ABORT
      ENDIF
    ENDIF ! of IF (is_master)
    
  END SUBROUTINE put_var_rgen     

END MODULE iostart
