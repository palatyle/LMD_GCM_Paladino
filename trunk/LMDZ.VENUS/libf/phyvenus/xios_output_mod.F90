MODULE xios_output_mod

 IMPLICIT NONE
 
 INTEGER,PRIVATE,SAVE :: time_it=0 ! store number of iterations with calls to XIOS since start
! does not need to be threadprivate; managed by omp master

 CHARACTER(LEN=*), PARAMETER :: context_id= "LMDZ" ! same as in context_lmdz_physics.xml
 
#ifdef CPP_XIOS

 INTERFACE send_xios_field
    MODULE PROCEDURE histwrite0d_xios,histwrite2d_xios,histwrite3d_xios
 END INTERFACE
 

CONTAINS

  SUBROUTINE initialize_xios_output(day,timeofday,dtphys,daysec,&
                                    presnivs,pseudoalt)
!  USE mod_phys_lmdz_para, only: gather, bcast, &
!                                jj_nb, jj_begin, jj_end, ii_begin, ii_end, &
!                                mpi_size, mpi_rank, klon_mpi, &
!                                is_sequential, is_south_pole_dyn
  USE mod_phys_lmdz_para, ONLY: jj_nb, jj_begin, jj_end, ii_begin, ii_end, &
                                mpi_size, mpi_rank, klon_mpi, &
                                is_sequential, is_south_pole_dyn
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, grid_type, unstructured
  USE print_control_mod, ONLY: lunout, prt_level
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
  USE nrtype, ONLY: pi
#ifdef CPP_XIOS
  USE xios
#endif
  USE wxios, ONLY: wxios_domain_param, wxios_domain_param_unstructured, wxios_closedef
  IMPLICIT NONE
  
  REAL,INTENT(IN) :: day ! Number of elapsed sols since reference Ls=0.
  REAL,INTENT(IN) :: timeofday ! "Universal time", given as fraction of sol (e.g.: 0.5 for noon).
  REAL,INTENT(IN) :: dtphys ! physics time step (s)
  REAL,INTENT(IN) :: daysec ! lengthof a standard day (s)
  REAL,INTENT(IN) :: presnivs(:) ! vertical grid approximate pressure (Pa)
  REAL,INTENT(IN) :: pseudoalt(:) ! vertical grid approximate altitude (km)
  
  
  INTEGER :: data_ibegin, data_iend
  TYPE(xios_duration) :: timestep
  TYPE(xios_date) :: time_origin
  TYPE(xios_date) :: start_date
  
!$OMP BARRIER
!$OMP MASTER

    ! 1. Declare available vertical axes to be used in output files:
    IF (prt_level>=10) WRITE(lunout,*) "initialize_xios_output: call xios_set_axis_attr for presnivs"
    CALL xios_set_axis_attr("presnivs", n_glo=size(presnivs), value=presnivs,&
                            unit="Pa",positive="down")
    IF (prt_level>=10) WRITE(lunout,*) "initialize_xios_output: call xios_set_axis_attr for altitude"
    CALL xios_set_axis_attr("altitude", n_glo=size(pseudoalt), value=pseudoalt,&
                            unit="km",positive="up")
    
    ! 2. Declare horizontal domain
    ! Set values for the mask:
!    IF (mpi_rank == 0) THEN
!        data_ibegin = 0
!    ELSE 
!        data_ibegin = ii_begin - 1
!    END IF

!    IF (mpi_rank == mpi_size-1) THEN
!        data_iend = nbp_lon
!    ELSE
!        data_iend = ii_end + 1
!    END IF

!    if (prt_level>=10) then
!      write(lunout,*) "initialize_xios_output: mpirank=",mpi_rank," iibegin=",ii_begin , " ii_end=",ii_end," jjbegin=",jj_begin," jj_nb=",jj_nb," jj_end=",jj_end
!      write(lunout,*) "initialize_xios_output: mpirank=",mpi_rank," nbp_lon=",nbp_lon," nbp_lat=",nbp_lat
!      write(lunout,*) "initialize_xios_output: mpirank=",mpi_rank," data_ibegin=",data_ibegin," data_iend=",data_iend
!      write(lunout,*) "initialize_xios_output: mpirank=",mpi_rank," data_ibegin=",data_ibegin," data_iend=",data_iend
!      write(lunout,*) "initialize_xios_output: mpirank=",mpi_rank," is_south_pole=",is_south_pole_dyn
!    endif

!$OMP END MASTER
!$OMP BARRIER
    ! Initialize the XIOS domain corresponding to this process:
    if (prt_level>=10) write(lunout,*) "initialize_xios_output: call wxios_domain_param"
!    CALL wxios_domain_param("dom_glo", is_sequential, nbp_lon, jj_nb, nbp_lon, nbp_lat, &
!                            1, nbp_lon, ii_begin, ii_end, jj_begin, jj_end,             &
!                            klon_mpi+2*(nbp_lon-1), data_ibegin, data_iend,             &
! VENUS IS SEEN UPSIDE DOWN, SO CENTRAL SYMMETRY TO PUT NORTH UP AGAIN
!                           -1.*lat_reg*(180./pi), -1.*lon_reg*(180./pi),                &
!                            is_south_pole_dyn,mpi_rank)

    IF (grid_type==unstructured) THEN
      CALL wxios_domain_param_unstructured("dom_glo",.true.)
    ELSE
      CALL wxios_domain_param("dom_glo",.true.)
    ENDIF

!$OMP MASTER
    ! 3. Declare calendar and time step
    if (prt_level>=10) then
     write(lunout,*) "initialize_xios_output: build calendar"
    endif
    timestep%second=nint(dtphys)
    if (nint(dtphys).ne.dtphys) then
      write(*,*) "initialize_xios_output: warning physics timestep is not an integer!"
    endif
    if (nint(daysec).ne.daysec) then
      write(*,*) "initialize_xios_output: warning day length is not an integer!"
    endif
    ! Important: do no operations involving dates and calendars
    ! before defining the calendar!
    CALL xios_define_calendar(type="user_defined", &
                              timestep=timestep, &
                              day_length=nint(daysec), &
                              month_lengths=[2]) ! one month, 2 days long

    ! time origin of the simulation (default: 1st year/1st month/1st day, Ls=0)
    time_origin=xios_date(1,1,1,0,0,0)
    CALL xios_set_time_origin(time_origin=time_origin)
    if (prt_level>=10) then
     write(lunout,*) "initialize_xios_output: time_origin=",time_origin
    endif

    ! Now define the start time of this simulation
    ! NB: we substract dtphys because we want to set the origin of the time axis
    start_date=time_origin+xios_duration(0,0,day,0,0,timeofday*daysec-dtphys)
    call xios_set_start_date(start_date=start_date)
    if (prt_level>=10) then
     write(lunout,*) "initialize_xios_output: start_date=",start_date
    endif
    
    ! 4. Finalize the context:
    if (prt_level>=10) write(*,*) "initialize_xios_output: call wxios_closedef"
    CALL wxios_closedef()

!$OMP END MASTER
!$OMP BARRIER
  
  END SUBROUTINE initialize_xios_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE finalize_xios_output
  USE xios
  IMPLICIT NONE
!$OMP BARRIER    
!$OMP MASTER
    CALL xios_context_finalize
!$OMP END MASTER    
!$OMP BARRIER    
  
  END SUBROUTINE finalize_xios_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE update_xios_timestep
  USE xios
  IMPLICIT NONE
    CALL set_xios_context
!$OMP MASTER
    time_it=time_it+1
    CALL xios_update_calendar(time_it)
!$OMP END MASTER    
  END SUBROUTINE update_xios_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE set_xios_context
  USE XIOS
  IMPLICIT NONE
    TYPE(xios_context) :: ctx_hdl

!$OMP MASTER
    CALL xios_get_handle(context_id,ctx_hdl)
    CALL xios_set_current_context(ctx_hdl)
!$OMP END MASTER    
  END SUBROUTINE set_xios_context

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE histwrite0d_xios(field_name,field)
  USE xios, ONLY: xios_send_field
  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE
  
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    REAL, INTENT(IN) :: field
    
    IF (prt_level >= 10) WRITE(lunout,*)'Begin histrwrite0d_xios ',trim(field_name)
    
!$OMP MASTER
    CALL xios_send_field(field_name,field)
!$OMP END MASTER
    
    IF (prt_level >= 10) WRITE(lunout,*)'End histrwrite0d_xios ',trim(field_name)
    
  END SUBROUTINE histwrite0d_xios

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE histwrite2d_xios(field_name,field)
  USE dimphy, only: klon
  USE mod_phys_lmdz_para, only: gather_omp, grid1Dto2D_mpi, &
                                jj_nb, klon_mpi
  USE xios, only: xios_send_field
  USE print_control_mod, ONLY: prt_level, lunout
  USE mod_grid_phy_lmdz, ONLY: nbp_lon
  IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: field_name
    REAL, DIMENSION(:), INTENT(IN) :: field
      
    REAL,DIMENSION(klon_mpi) :: buffer_omp
    REAL :: Field2d(nbp_lon,jj_nb)

    IF (prt_level >= 10) WRITE(lunout,*)'Begin histrwrite2d_xios ',trim(field_name)

    IF (SIZE(field)/=klon) CALL abort_physic('iophy::histwrite2d_xios','Field first DIMENSION not equal to klon',1)
    
    CALL Gather_omp(field,buffer_omp)    
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,Field2d)
    
    CALL xios_send_field(field_name, Field2d)
!$OMP END MASTER   

    IF (prt_level >= 10) WRITE(lunout,*)'End histrwrite2d_xios ',trim(field_name)
  END SUBROUTINE histwrite2d_xios

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE histwrite3d_xios(field_name, field)
  USE dimphy, only: klon, klev
  USE mod_phys_lmdz_para, only: gather_omp, grid1Dto2D_mpi, &
                                jj_nb, klon_mpi
  USE xios, only: xios_send_field
  USE print_control_mod, ONLY: prt_level,lunout
  USE mod_grid_phy_lmdz, ONLY: nbp_lon

  IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: field_name
    REAL, DIMENSION(:,:), INTENT(IN) :: field ! --> field(klon,:)

    REAL,DIMENSION(klon_mpi,SIZE(field,2)) :: buffer_omp
    REAL :: Field3d(nbp_lon,jj_nb,SIZE(field,2))
    INTEGER :: ip, n, nlev

  IF (prt_level >= 10) write(lunout,*)'Begin histrwrite3d_xios ',trim(field_name)

    !Et on.... écrit
    IF (SIZE(field,1)/=klon) CALL abort_physic('iophy::histwrite3d','Field first DIMENSION not equal to klon',1)
    nlev=SIZE(field,2)


    CALL Gather_omp(field,buffer_omp)
!$OMP MASTER
    CALL grid1Dto2D_mpi(buffer_omp,field3d)

    CALL xios_send_field(field_name, Field3d(:,:,1:nlev))
!$OMP END MASTER   

    IF (prt_level >= 10) write(lunout,*)'End histrwrite3d_xios ',trim(field_name)
  END SUBROUTINE histwrite3d_xios

#endif

END MODULE xios_output_mod
