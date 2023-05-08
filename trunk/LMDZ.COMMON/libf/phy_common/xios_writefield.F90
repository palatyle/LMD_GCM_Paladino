MODULE xios_writefield_mod
#ifdef CPP_XIOS
  PRIVATE
  
  INTEGER,PARAMETER :: size_name=255
  INTEGER, PARAMETER :: MaxWriteField = 100
  CHARACTER(len=size_name), SAVE ::  FieldName(MaxWriteField) 
!$OMP THREADPRIVATE(FieldName)
  INTEGER,                  SAVE  :: FieldIt (MaxWriteField)
!$OMP THREADPRIVATE(FieldIt)
  INTEGER,SAVE :: NbField = 0
!$OMP THREADPRIVATE(NbField) 

  INTERFACE xios_writefield
    MODULE PROCEDURE xios_writefield2d,xios_writefield3d
  END INTERFACE

  LOGICAL,SAVE :: output_native_grid
!$OMP THREADPRIVATE( output_native_grid)  

  INTEGER,SAVE :: ni_glo , nj_glo
!$OMP THREADPRIVATE( ni_glo, nj_glo)  

   
  PUBLIC xios_writefield

CONTAINS
  
  FUNCTION NameId(name_in)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)      :: name_in
    INTEGER                          :: nameId
    CHARACTER(LEN=size_name)         :: name
    INTEGER :: n  
    
    name=name_in
     
    DO n=1,NbField
      IF (name==fieldName(n)) THEN
        nameId=n
        RETURN
      ENDIF
    ENDDO
      
    nameId=0
    RETURN
  END FUNCTION NameId
     
 
 
  SUBROUTINE xios_writefield2d(field,name_in)
  USE dimphy
  USE mod_phys_lmdz_para
  USE xios
  USE print_control_mod, ONLY:  lunout
  USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo, grid_type, unstructured, regular_lonlat
  USE wxios, ONLY: wxios_domain_param_unstructured, wxios_domain_param,  wxios_set_context
  USE ioipsl
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: field
  CHARACTER(LEN=*)                 :: name_in
  CHARACTER(LEN=size_name)         :: name 
  TYPE(xios_context) :: xios_ctx
  TYPE(xios_domaingroup) :: domain_definition
  TYPE(xios_domaingroup) :: my_domaingroup
  TYPE(xios_domain)      :: my_domain
  TYPE(xios_domain)      :: my_regular_domain
  TYPE(xios_filegroup)   :: file_definition
  TYPE(xios_file)        :: my_file
  TYPE(xios_field)       :: my_field
  TYPE(xios_fieldgroup)  :: field_definition
  TYPE(xios_generate_rectilinear_domain) :: generate_domain 
  TYPE(xios_interpolate_domain) :: interpolate_domain
  INTEGER                :: id
  REAL,DIMENSION(klon_mpi) :: field_mpi  
  REAL :: Field2d(nbp_lon,jj_nb)
  
    
    IF (size(field,1) /= klon_omp) THEN
      WRITE(lunout,*) 'xios_writefield :: '//FieldName//' is not on the model grid'
      RETURN
    ENDIF
  
    name = TRIM(ADJUSTL(name_in))
    id=nameId(name)
    
    IF (id/=0) THEN
      IF (is_omp_master) THEN
        CALL xios_get_handle('context_lmdz_'//TRIM(name), xios_ctx)    !Récupération
        CALL xios_set_current_context(xios_ctx)            !Activation
      ENDIF
    ELSE
      output_native_grid=.FALSE.
      ni_glo=0
      nj_glo=0
      IF (is_master) CALL getin("xios_writefield_nlon",ni_glo)
      IF (is_master) CALL getin("xios_writefield_nlat",nj_glo)
      CALL bcast(ni_glo)
      CALL bcast(nj_glo)
      IF (ni_glo==0 .OR. nj_glo==0) output_native_grid=.TRUE.

      IF (is_omp_master) THEN
        CALL xios_context_initialize('context_lmdz_'//TRIM(name), COMM_LMDZ_PHY)
        CALL xios_get_handle('context_lmdz_'//TRIM(name), xios_ctx)    !Récupération
        CALL xios_set_current_context(xios_ctx)            !Activation
      
        CALL xios_define_calendar("D360")
        CALL xios_set_start_date(xios_date(2000,1,1,0,0,0))
        CALL xios_set_time_origin(xios_date(2000,1,1,0,0,0))
        CALL xios_set_timestep(xios_second)
      
        CALL xios_get_handle("domain_definition",domain_definition)
        CALL xios_add_child(domain_definition,my_domaingroup,"domaingroup")
        CALL xios_add_child(my_domaingroup,my_domain,"domain")
        IF (grid_type==unstructured .AND. .NOT. output_native_grid ) THEN
          CALL xios_add_child(domain_definition, my_regular_domain, "regular_domain")
          CALL xios_set_attr(my_regular_domain,ni_glo=ni_glo, nj_glo=nj_glo, type="rectilinear")
          CALL xios_add_child(my_regular_domain,generate_domain)
          CALL xios_set_attr(generate_domain,lon_start=-180., lat_start=90., lat_end=-90.)
          CALL xios_add_child(my_regular_domain,interpolate_domain)
        ENDIF
          
      ENDIF
      
      IF (grid_type==regular_lonlat) THEN
        CALL wxios_domain_param("domain")
      ELSE IF (grid_type==unstructured) THEN
        CALL wxios_domain_param_unstructured("domaingroup")
      ENDIF       
      
      IF (is_omp_master) THEN
        CALL xios_get_handle("file_definition",file_definition)
        CALL xios_add_child(file_definition, my_file)
        CALL xios_set_attr(my_file,name=TRIM(name),output_freq=xios_timestep,sync_freq=xios_timestep,type="one_file")
        CALL xios_get_handle("field_definition",field_definition) 
        CALL xios_add_child(field_definition, my_field, TRIM(name))
        CALL xios_set_attr(my_field,domain_ref="domain",operation="instant")
        CALL xios_add_child(my_file, my_field)
        CALL xios_set_attr(my_field,field_ref=TRIM(name))
        IF (grid_type==unstructured .AND. .NOT. output_native_grid) CALL xios_set_attr(my_field,domain_ref="regular_domain")

        CALL xios_close_context_definition()
      
        NbField=NbField+1
        FieldName(NbField)=TRIM(name)
        FieldIt(NbField)=0
        id=NbField
      ENDIF
    ENDIF
 
    CALL Gather_omp(field,field_mpi)

    IF (is_omp_master) THEN
      FieldIt(id)=FieldIt(id)+1
      CALL xios_update_calendar(FieldIt(id)) 

      IF (grid_type==regular_lonlat) THEN
        CALL grid1Dto2D_mpi(field_mpi,Field2d)   
        CALL xios_send_field(TRIM(name), Field2d)
      ELSE IF (grid_type==unstructured) THEN
        CALL xios_send_field(TRIM(name), field_mpi)
      ENDIF
    ENDIF
    
    CALL wxios_set_context()
  
  END SUBROUTINE  xios_writefield2d

  SUBROUTINE xios_writefield3d(field,name_in)
  USE dimphy
  USE mod_phys_lmdz_para
  USE xios
  USE print_control_mod, ONLY:  lunout
  USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo, grid_type, unstructured, regular_lonlat
  USE wxios, ONLY: wxios_domain_param_unstructured, wxios_domain_param,  wxios_set_context
  USE ioipsl
  IMPLICIT NONE
  REAL, DIMENSION(:,:), INTENT(IN) :: field
  CHARACTER(LEN=*)                 :: name_in
  CHARACTER(LEN=size_name)         :: name 
  TYPE(xios_context) :: xios_ctx
  TYPE(xios_domaingroup) :: domain_definition
  TYPE(xios_domaingroup) :: my_domaingroup
  TYPE(xios_domain)      :: my_domain
  TYPE(xios_domain)      :: my_regular_domain
  TYPE(xios_filegroup)   :: file_definition
  TYPE(xios_file)        :: my_file
  TYPE(xios_field)       :: my_field
  TYPE(xios_fieldgroup)  :: field_definition
  TYPE(xios_generate_rectilinear_domain) :: generate_domain 
  TYPE(xios_interpolate_domain) :: interpolate_domain
  TYPE(xios_axisgroup)   :: axis_definition
  TYPE(xios_axis)        :: my_axis
  INTEGER                :: id
  REAL,DIMENSION(klon_mpi,size(field,2)) :: field_mpi  
  REAL :: Field2d(nbp_lon,jj_nb,size(field,2))
  REAL :: axis_value(size(field,2))
  INTEGER :: i
    
    IF (size(field,1) /= klon_omp) THEN
      WRITE(lunout,*) 'xios_writefield :: '//FieldName//' is not on the model grid'
      RETURN
    ENDIF
  
    name = TRIM(ADJUSTL(name_in))
    id=nameId(name)
    
    IF (id/=0) THEN
   
      IF (is_omp_master) CALL xios_get_handle('context_lmdz_'//TRIM(name), xios_ctx)    !Récupération
      IF (is_omp_master) CALL xios_set_current_context(xios_ctx)            !Activation
   
    ELSE
      output_native_grid=.FALSE.
      ni_glo=0
      nj_glo=0
      IF (is_master) CALL getin("xios_writefield_nlon",ni_glo)
      IF (is_master) CALL getin("xios_writefield_nlat",nj_glo)
      CALL bcast(ni_glo)
      CALL bcast(nj_glo)
      IF (ni_glo==0 .OR. nj_glo==0) output_native_grid=.TRUE.

      IF (is_omp_master) THEN
        CALL xios_context_initialize('context_lmdz_'//TRIM(name), COMM_LMDZ_PHY)
        CALL xios_get_handle('context_lmdz_'//TRIM(name), xios_ctx)    !Récupération
        CALL xios_set_current_context(xios_ctx)            !Activation
      
        CALL xios_define_calendar("D360")
        CALL xios_set_start_date(xios_date(2000,1,1,0,0,0))
        CALL xios_set_time_origin(xios_date(2000,1,1,0,0,0))
        CALL xios_set_timestep(xios_second)
      
        CALL xios_get_handle("domain_definition",domain_definition)
        CALL xios_add_child(domain_definition,my_domaingroup,"domaingroup")
        CALL xios_add_child(my_domaingroup,my_domain,"domain")
        IF (grid_type==unstructured .AND. .NOT. output_native_grid) THEN
          CALL xios_add_child(domain_definition, my_regular_domain, "regular_domain")
          CALL xios_set_attr(my_regular_domain,ni_glo=ni_glo, nj_glo=nj_glo, type="rectilinear")
          CALL xios_add_child(my_regular_domain,generate_domain)
          CALL xios_set_attr(generate_domain,lon_start=-180., lat_start=90., lat_end=-90.)
          CALL xios_add_child(my_regular_domain,interpolate_domain)
        ENDIF
        CALL xios_get_handle("axis_definition",axis_definition)
        CALL xios_add_child(axis_definition,my_axis,"axis")
        axis_value=(/(i,i=1,size(field,2))/)
        CALL xios_set_attr(my_axis,name="z",n_glo=size(field,2),value=axis_value, unit="level" )        
      ENDIF
      
      IF (grid_type==regular_lonlat) THEN
        CALL wxios_domain_param("domain")
      ELSE IF (grid_type==unstructured) THEN
        CALL wxios_domain_param_unstructured("domaingroup")
      ENDIF       
      
      IF (is_omp_master) THEN
        CALL xios_get_handle("file_definition",file_definition)
        CALL xios_add_child(file_definition, my_file)
        CALL xios_set_attr(my_file,name=TRIM(name),output_freq=xios_timestep,sync_freq=xios_timestep,type="one_file")
        CALL xios_get_handle("field_definition",field_definition) 
        CALL xios_add_child(field_definition, my_field, TRIM(name))
        CALL xios_set_attr(my_field,domain_ref="domain",axis_ref="axis", operation="instant")
        CALL xios_add_child(my_file, my_field)
        CALL xios_set_attr(my_field,field_ref=TRIM(name))
        IF (grid_type==unstructured .AND. .NOT. output_native_grid) CALL xios_set_attr(my_field,domain_ref="regular_domain")

        CALL xios_close_context_definition()
      
        NbField=NbField+1
        FieldName(NbField)=TRIM(name)
        FieldIt(NbField)=0
        id=NbField
      ENDIF
    ENDIF
 
    CALL Gather_omp(field,field_mpi)

    IF (is_omp_master) THEN
      FieldIt(id)=FieldIt(id)+1
      CALL xios_update_calendar(FieldIt(id)) 

      IF (grid_type==regular_lonlat) THEN
        CALL grid1Dto2D_mpi(field_mpi,Field2d)   
        CALL xios_send_field(TRIM(name), Field2d)
      ELSE IF (grid_type==unstructured) THEN
        CALL xios_send_field(TRIM(name), field_mpi)
      ENDIF
    ENDIF
    
    CALL wxios_set_context()
  
  END SUBROUTINE  xios_writefield3d

#endif
END MODULE
