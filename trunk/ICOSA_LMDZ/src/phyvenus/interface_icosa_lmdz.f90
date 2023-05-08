MODULE interface_icosa_lmdz_mod

  USE field_mod, ONLY: t_field
  USE transfert_mod, ONLY: t_message 
  
 
  TYPE(t_message),SAVE :: req_u
  TYPE(t_message),SAVE :: req_dps0, req_dulon0, req_dulat0, req_dTemp0, req_dq0

  TYPE(t_field),POINTER,SAVE :: f_p(:) 
  TYPE(t_field),POINTER,SAVE :: f_pks(:)  
  TYPE(t_field),POINTER,SAVE :: f_pk(:)  
  TYPE(t_field),POINTER,SAVE :: f_p_layer(:)    
  TYPE(t_field),POINTER,SAVE :: f_theta(:)   
  TYPE(t_field),POINTER,SAVE :: f_phi(:)   
  TYPE(t_field),POINTER,SAVE :: f_Temp(:)    
  TYPE(t_field),POINTER,SAVE :: f_ulon(:)    
  TYPE(t_field),POINTER,SAVE :: f_ulat(:)   
  TYPE(t_field),POINTER,SAVE :: f_dulon(:)
  TYPE(t_field),POINTER,SAVE :: f_dulat(:)
  TYPE(t_field),POINTER,SAVE :: f_dTemp(:)
  TYPE(t_field),POINTER,SAVE :: f_dq(:)
  TYPE(t_field),POINTER,SAVE :: f_dps(:)
  TYPE(t_field),POINTER,SAVE :: f_duc(:)
  TYPE(t_field),POINTER,SAVE :: f_bounds_lon(:)
  TYPE(t_field),POINTER,SAVE :: f_bounds_lat(:)

  INTEGER,SAVE :: start_clock
  INTEGER,SAVE :: stop_clock
  INTEGER,SAVE :: count_clock=0
  
  REAL,SAVE :: day_length ! length of a day (s)
!  REAL,SAVE :: year_length ! length of a year (s)
  INTEGER,SAVE :: start_day ! reference sol value at beginning of the run wrt Ls=0
  
  INTEGER,SAVE :: nbp_phys
  INTEGER,SAVE :: nbp_phys_glo
  
  CHARACTER(len=30),SAVE,ALLOCATABLE :: tname(:) ! tracer names
  CHARACTER(len=33),SAVE,ALLOCATABLE :: ttext(:) ! tracer long name for diagnostics
  REAL,SAVE :: pday ! number of ellapsed sols since Ls=0
  REAL,SAVE :: ptime ! "universal time" as fraction of sol (e.g. 0.5 for noon)


CONTAINS

  SUBROUTINE initialize_physics
  USE distrib_icosa_lmdz_mod, ONLY : init_distrib_icosa_lmdz, transfer_icosa_to_lmdz
! from dynamico
  USE domain_mod
  USE dimensions
  USE mpi_mod
  USE mpipara
  USE disvert_mod
  USE xios_mod
  USE time_mod , init_time_icosa=> init_time 
  USE transfert_mod
  
! from LMDZ
  USE mod_grid_phy_lmdz, ONLY : unstructured
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp
  USE transfert_mod
  USE physics_distribution_mod, ONLY : init_physics_distribution
   
  
  IMPLICIT NONE
  INTEGER  :: ind,i,j,ij,pos
  REAL(rstd),POINTER :: bounds_lon(:,:)
  REAL(rstd),POINTER :: bounds_lat(:,:)
  
  REAL(rstd),ALLOCATABLE :: latfi(:)
  REAL(rstd),ALLOCATABLE :: lonfi(:)
  REAL(rstd),ALLOCATABLE :: airefi(:)
  REAL(rstd),ALLOCATABLE :: bounds_latfi(:,:)
  REAL(rstd),ALLOCATABLE :: bounds_lonfi(:,:)
  REAL(rstd) :: pseudoalt(llm)

  INTEGER :: nbp_phys, nbp_phys_glo
  
!$OMP PARALLEL
    CALL allocate_field(f_bounds_lon,field_t,type_real,6)
    CALL allocate_field(f_bounds_lat,field_t,type_real,6)
    CALL allocate_field(f_p,field_t,type_real,llm+1,name="p_in")
    CALL allocate_field(f_pks,field_t,type_real)
    CALL allocate_field(f_pk,field_t,type_real,llm)
    CALL allocate_field(f_p_layer,field_t,type_real,llm,name="p_layer_in")
    CALL allocate_field(f_theta,field_t,type_real,llm)
    CALL allocate_field(f_phi,field_t,type_real,llm,name="phi_in")
    CALL allocate_field(f_Temp,field_t,type_real,llm,name="Temp_in")
    CALL allocate_field(f_ulon,field_t,type_real,llm,name="ulon_in")
    CALL allocate_field(f_ulat,field_t,type_real,llm,name="ulat_in")
    CALL allocate_field(f_dulon,field_t,type_real,llm,name="dulon_out")
    CALL allocate_field(f_dulat,field_t,type_real,llm,name="dulat_out")
    CALL allocate_field(f_dTemp,field_t,type_real,llm,name="dTemp_out")
    CALL allocate_field(f_dq,field_t,type_real,llm,nqtot,name="dq_out")
    CALL allocate_field(f_dps,field_t,type_real,name="dps_out")
    CALL allocate_field(f_duc,field_t,type_real,3,llm)    

    CALL init_message(f_dps,req_i0,req_dps0)
    CALL init_message(f_dulon,req_i0,req_dulon0)
    CALL init_message(f_dulat,req_i0,req_dulat0)
    CALL init_message(f_dTemp,req_i0,req_dTemp0)
    CALL init_message(f_dq,req_i0,req_dq0)
!$OMP END PARALLEL    

    nbp_phys=0
    DO ind=1,ndomain
      CALL swap_dimensions(ind)
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          IF (domain(ind)%own(i,j)) nbp_phys=nbp_phys+1
        ENDDO
      ENDDO
    ENDDO
    

!initialize LMDZ5 physic mpi decomposition
    CALL MPI_ALLREDUCE(nbp_phys,nbp_phys_glo,1,MPI_INTEGER,MPI_SUM,comm_icosa,ierr)
    CALL init_physics_distribution(unstructured, 6, nbp_phys, 1, nbp_phys_glo, llm, comm_icosa)
    
    DO ind=1,ndomain
        CALL swap_dimensions(ind)
        CALL swap_geometry(ind)
        bounds_lon=f_bounds_lon(ind)
        bounds_lat=f_bounds_lat(ind)
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            CALL xyz2lonlat(xyz_v(ij+z_rup,:), bounds_lon(ij,1), bounds_lat(ij,1))
            CALL xyz2lonlat(xyz_v(ij+z_up,:), bounds_lon(ij,2), bounds_lat(ij,2))
            CALL xyz2lonlat(xyz_v(ij+z_lup,:), bounds_lon(ij,3), bounds_lat(ij,3))
            CALL xyz2lonlat(xyz_v(ij+z_ldown,:), bounds_lon(ij,4), bounds_lat(ij,4))
            CALL xyz2lonlat(xyz_v(ij+z_down,:), bounds_lon(ij,5), bounds_lat(ij,5))
            CALL xyz2lonlat(xyz_v(ij+z_rdown,:), bounds_lon(ij,6), bounds_lat(ij,6))
         ENDDO
       ENDDO            
    ENDDO
          
!$OMP PARALLEL
    CALL initialize_physics_omp
!$OMP END PARALLEL            

    CALL xios_set_context    


     

  END SUBROUTINE initialize_physics


  SUBROUTINE initialize_physics_omp
  USE distrib_icosa_lmdz_mod, ONLY : init_distrib_icosa_lmdz, transfer_icosa_to_lmdz
! from dynamico
  USE domain_mod
  USE dimensions
  USE mpi_mod
  USE mpipara
  USE disvert_mod
  USE xios_mod
  USE time_mod , ONLY: init_time_icosa=> init_time, dt, itaumax, itau_physics
  USE omp_para

! from LMDZ
  USE mod_grid_phy_lmdz, ONLY : unstructured
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp
  USE time_phylmdz_mod, ONLY: init_time_lmdz => init_time
  USE transfert_mod
  USE physics_distribution_mod, ONLY : init_physics_distribution
  USE dimphy, ONLY: init_dimphy
  USE geometry_mod, ONLY : init_geometry
  USE vertical_layers_mod, ONLY : init_vertical_layers
!  USE planete_mod, ONLY: ini_planete_mod
  USE cpdet_phy_mod, ONLY: init_cpdet_phy
  USE infotrac_phy, ONLY: init_infotrac_phy
 
  USE netcdf  
  
  IMPLICIT NONE



  INTEGER  :: ind,i,j,k,ij,pos
  REAL(rstd),POINTER :: bounds_lon(:,:)
  REAL(rstd),POINTER :: bounds_lat(:,:)
  
  REAL(rstd),ALLOCATABLE :: latfi(:)
  REAL(rstd),ALLOCATABLE :: lonfi(:)
  REAL(rstd),ALLOCATABLE :: airefi(:)
  REAL(rstd),ALLOCATABLE :: bounds_latfi(:,:)
  REAL(rstd),ALLOCATABLE :: bounds_lonfi(:,:)
  REAL(rstd),ALLOCATABLE :: ind_cell_glo(:)
  REAL(rstd),ALLOCATABLE :: dx(:)
  REAL(rstd),ALLOCATABLE :: dy(:)

  REAL(rstd) :: pseudoalt(llm)
  REAL(rstd) :: aps(llm)
  REAL(rstd) :: bps(llm)
  real(rstd) :: scaleheight

  ! For Cp(T)
  REAL(rstd) :: nu_venus
  REAL(rstd) :: t0_venus

  ! Calendar related stuff
  REAL(rstd) :: ptimestep
  REAL(rstd) :: run_length  
  INTEGER :: annee_ref  
  INTEGER :: day_ref    
  INTEGER :: day_ini
  INTEGER :: day_end
  INTEGER :: raz_date

  ! Tracer related information
  INTEGER :: iflag_trac
  CHARACTER(len=4)              :: type_trac
!  CHARACTER(len=30),ALLOCATABLE :: tname(:)    ! tracer short name for restart and diagnostics
!  CHARACTER(len=33),ALLOCATABLE :: ttext(:)     ! tracer long name for diagnostics
  TYPE(t_field),POINTER,SAVE    :: f_ind_cell_glo(:)
  
  INTEGER :: iflag_phys    
  INTEGER :: nq

  !! to get values from startphy.nc controle array
  !! ---------------------------------------------
  logical :: startphy_file
  ! NetCDF stuff
  integer :: status ! NetCDF return code
  integer :: ncid ! NetCDF file ID
  integer :: varid ! NetCDF variable ID
  real :: tab_cntrl(100)

    CALL init_distrib_icosa_lmdz
    
    ALLOCATE(latfi(klon_omp))
    ALLOCATE(lonfi(klon_omp))
    ALLOCATE(airefi(klon_omp))
    ALLOCATE(bounds_latfi(klon_omp,6))
    ALLOCATE(bounds_lonfi(klon_omp,6))
    ALLOCATE(ind_cell_glo(klon_omp))
    ALLOCATE(dx(klon_omp))
    ALLOCATE(dy(klon_omp))

    CALL transfer_icosa_to_lmdz(geom%lat_i,latfi)
    CALL transfer_icosa_to_lmdz(geom%lon_i,lonfi)
    CALL transfer_icosa_to_lmdz(f_bounds_lat,bounds_latfi)
    CALL transfer_icosa_to_lmdz(f_bounds_lon,bounds_lonfi)
    CALL transfer_icosa_to_lmdz(geom%Ai,airefi)

    CALL allocate_field(f_ind_cell_glo,field_t,type_real)
    
    DO ind=1,ndomain
      IF (.NOT. assigned_domain(ind)  .OR. .NOT. is_omp_level_master ) CYCLE
      CALL swap_dimensions(ind)
      CALL swap_geometry(ind)
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          ij=(j-1)*iim+i
          f_ind_cell_glo(ind)%rval2d(ij)=domain(ind)%assign_cell_glo(i,j)
        ENDDO
      ENDDO
    ENDDO

      
    CALL transfer_icosa_to_lmdz(f_ind_cell_glo,ind_cell_glo)
    CALL deallocate_field(f_ind_cell_glo)
      
              
    ! Initialize dimphy module
    CALL init_dimphy(klon_omp,llm)

    ! Dummy initializations for dx(),dy() In principle these are not used
    ! in the physics; but this should be checked further...
    dx(:)=1 ; dy(:)=1
    CALL init_geometry(klon_omp,lonfi,latfi,bounds_lonfi,bounds_latfi,&
                       airefi,INT(ind_cell_glo),dx,dy)

    scaleheight=scale_height/1000. ! Atmospheric scale height (km)
    aps(1:llm)=0.5*(ap(1:llm)+ap(2:llm+1))
    bps(1:llm)=0.5*(bp(1:llm)+bp(2:llm+1))
    pseudoalt(:)=-scaleheight*log(presnivs(:)/preff)
    CALL init_vertical_layers(llm,preff,scaleheight,ap,bp,aps,bps,presnivs,pseudoalt)

    ! Initialize planet_mod (quite redundant wrt vertical_levels...)
!    CALL ini_planete_mod(llm,preff,ap,bp)

    ! Initialize tracer names, numbers, etc. for physics

! init tracers model for standard lmdz case
!$OMP MASTER
    ALLOCATE(tname(nqtot))
    ALLOCATE(ttext(nqtot))
!$OMP END MASTER
!$OMP BARRIER

! read tname() from traceur.def file
    IF (is_mpi_root) THEN
!$OMP MASTER
    OPEN(unit=42,file="traceur.def",form="formatted",status="old",iostat=ierr)
    IF (ierr==0) THEN
      READ(42,*) nq ! should be the same as nqtot
      IF (nq /= nqtot) THEN
        WRITE(*,*) "Error: number of tracers in tracer.def should match nqtot!"
        WRITE(*,*) "       will just use nqtot=",nqtot," tracers"
      ENDIF
      DO i=1,nqtot
        READ(42,*) j,k,tname(i)
        ttext(i)=trim(tname(i))//"VL1"
      ENDDO
      CLOSE(42)
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
    ENDIF ! of (is_mpi_root)

    DO i=1,nqtot
      CALL bcast(tname(i))
      CALL bcast(ttext(i))
    ENDDO

   ! Get/set some constants for the physics


    startphy_file=.true.
    CALL getin('startphy_file',startphy_file)

    ! value in physics daysec=10087066.76s, important for solar radiation, 
    ! for the rest ptime/timeofday is what matters. To well estimate them,
    ! day_length must be multiple of dt*itau_physics. Error order:1e-4.
    day_length=10087000   
    CALL getin('day_length',day_length)
    
    run_length=day_length ! default
    CALL getin('run_length',run_length)
    
    raz_date=0 ! default: 0: no change in date
    CALL getin('raz_date',raz_date)

    IF (startphy_file) THEN
      ! Read in some information from the startphy.nc file

      IF (is_mpi_root) THEN
!$OMP MASTER      
      status=nf90_open('startphy.nc',NF90_NOWRITE,ncid)
      if (status.ne.nf90_noerr) then
        write(*,*)"Failed to open startphy.nc"
        write(*,*)trim(nf90_strerror(status))
        stop
      endif

      status=nf90_inq_varid(ncid,"controle",varid)
      if (status.ne.nf90_noerr) then
        write(*,*)"Failed to find controle variable"
        write(*,*)trim(nf90_strerror(status))
        stop
      endif

      status=nf90_get_var(ncid,varid,tab_cntrl)
      ! extract needed variables from tab_cntrl 
      day_ini=tab_cntrl(13)
      annee_ref=tab_cntrl(14)
      ! check if tab_cntrl(1), the stored physics time step
      ! is the same as the current physics time step (within roundoff precision)
      if (abs(((itau_physics*dt)-tab_cntrl(1))/(itau_physics*dt))<=1.e-8) then
        ! Everything OK
        ptime=modulo((tab_cntrl(15)*tab_cntrl(1))/day_length,1.0)
      else ! unless raz_date == 1 , we have a problem
        if (raz_date==1) then
          ! we reset date to midnight at lon=0 in the physics
          ptime=0.0
        else
          write(*,*)"Error: physics time step in startphy.nc is different"
          write(*,*)"       from what is specified via run_icosa.def"
          write(*,*)"       From run_icosa.def:    ",itau_physics*dt
          write(*,*)"       From startphy.nc file: ",tab_cntrl(1)
          write(*,*)"       You must reset date to midnight at lon=0"
          write(*,*)"       by specifying raz_date=1 in your def file"
          call abort
        endif
      endif

      status=nf90_close(ncid)
!$OMP END MASTER      
!$OMP BARRIER
      ENDIF ! of !IF (is_mpi_root)
      
      CALL bcast(day_ini)
      CALL bcast(annee_ref)
      CALL bcast(ptime)

    ELSE
      ! required information that is not in tab_cntrl
      ! has to be default or read from def files
      day_ini=0
      annee_ref=1
      CALL getin('annee_ref',annee_ref)
      ptime=0.
    ENDIF

    day_end=day_ini+nint(run_length/day_length)

    ! Other required values which have to be read from def files
    day_ref=1
    CALL getin('day_ref',day_ref)
    iflag_trac=0
    CALL getin('iflag_trac',iflag_trac)


    ! Initialize some physical constants
    CALL suphec

    ! Initialize cpdet_phy module
    nu_venus=0.35
    t0_venus=460.
    CALL init_cpdet_phy(cpp,nu_venus,t0_venus)
  
    ! Initialize some "temporal and calendar" related variables
    ptimestep=itau_physics*dt
    CALL init_time_lmdz(annee_ref,day_ref,day_ini,day_end,ptimestep)
  
    ! Initialize tracers in physics
    CALL init_infotrac_phy(iflag_trac,nqtot,tname,ttext)
    
    
    ! Initializations of some module variables
    
!$OMP MASTER
    ! initialize pday
    pday = day_ini
!$OMP END MASTER
!$OMP BARRIER
  
  END SUBROUTINE  initialize_physics_omp 
  
  


  SUBROUTINE physics
  USE icosa
  USE time_mod
  USE disvert_mod
  USE transfert_mod
  USE mpipara
  USE xios_mod
  USE wxios
  USE trace
  USE distrib_icosa_lmdz_mod, ONLY : transfer_icosa_to_lmdz, transfer_lmdz_to_icosa
  USE physics_external_mod, ONLY : it, f_phis, f_ps, f_theta_rhodz, f_u, f_wflux, f_q 
  USE write_field_mod
  USE checksum_mod
! from LMDZ
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp
  USE geometry_mod, ONLY : cell_area
  USE physiq_mod, ONLY: physiq
  IMPLICIT NONE
  
    REAL(rstd),POINTER :: phis(:)
    REAL(rstd),POINTER :: ps(:)
    REAL(rstd),POINTER :: theta_rhodz(:,:,:)
    REAL(rstd),POINTER :: u(:,:)
    REAL(rstd),POINTER :: wflux(:,:)
    REAL(rstd),POINTER :: q(:,:,:)
    REAL(rstd),POINTER :: p(:,:)
    REAL(rstd),POINTER :: pks(:)
    REAL(rstd),POINTER :: pk(:,:)
    REAL(rstd),POINTER :: p_layer(:,:)
    REAL(rstd),POINTER :: theta(:,:)
    REAL(rstd),POINTER :: phi(:,:)
    REAL(rstd),POINTER :: Temp(:,:)
    REAL(rstd),POINTER :: ulon(:,:)
    REAL(rstd),POINTER :: ulat(:,:)
    REAL(rstd),POINTER :: dulon(:,:)
    REAL(rstd),POINTER :: dulat(:,:)
    REAL(rstd),POINTER :: dTemp(:,:)
    REAL(rstd),POINTER :: dq(:,:,:)
    REAL(rstd),POINTER :: dps(:)
    REAL(rstd),POINTER :: duc(:,:,:)


    INTEGER :: ind,l
    
    REAL(rstd),ALLOCATABLE,SAVE :: ps_phy(:)
!$OMP THREADPRIVATE(ps_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: p_phy(:,:)
!$OMP THREADPRIVATE(p_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: p_layer_phy(:,:)
!$OMP THREADPRIVATE(p_layer_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: pk_phy(:,:)
!$OMP THREADPRIVATE(p_layer_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: Temp_phy(:,:)
!$OMP THREADPRIVATE(Temp_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: phis_phy(:)
!$OMP THREADPRIVATE(phis_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: phi_phy(:,:)
!$OMP THREADPRIVATE(phi_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: ulon_phy(:,:)
!$OMP THREADPRIVATE(ulon_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: ulat_phy(:,:)
!$OMP THREADPRIVATE(ulat_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: q_phy(:,:,:)
!$OMP THREADPRIVATE(q_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: wflux_phy(:,:)
!$OMP THREADPRIVATE(wflux_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: dulon_phy(:,:)
!$OMP THREADPRIVATE(dulon_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: dulat_phy(:,:)
!$OMP THREADPRIVATE(dulat_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: dTemp_phy(:,:)
!$OMP THREADPRIVATE(dTemp_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: dq_phy(:,:,:)
!$OMP THREADPRIVATE(dq_phy)
    REAL(rstd),ALLOCATABLE,SAVE :: dps_phy(:)
!$OMP THREADPRIVATE(dps_phy)
    REAL(rstd)   :: dtphy 
    LOGICAL      :: debut
    LOGICAL      :: lafin
    LOGICAL,SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)
    REAL(rstd),ALLOCATABLE,SAVE :: plevmoy(:)
!$OMP THREADPRIVATE(plevmoy)
    REAL(rstd),ALLOCATABLE,SAVE :: tmoy(:)
!$OMP THREADPRIVATE(tmoy)
    
    IF(first) THEN 
      debut=.TRUE.
    ELSE 
      debut=.FALSE.
    ENDIF


    IF(it-itau0>=itaumax) THEN 
      lafin=.TRUE.
    ELSE 
      lafin=.FALSE.
    ENDIF

    IF (first) THEN
      first=.FALSE.
      CALL init_message(f_u,req_e1_vect,req_u)
      ALLOCATE(ps_phy(klon_omp))
      ALLOCATE(p_phy(klon_omp,llm+1))
      ALLOCATE(p_layer_phy(klon_omp,llm))
      ALLOCATE(pk_phy(klon_omp,llm))
      ALLOCATE(Temp_phy(klon_omp,llm))
      ALLOCATE(phis_phy(klon_omp))
      ALLOCATE(phi_phy(klon_omp,llm))
      ALLOCATE(ulon_phy(klon_omp,llm))
      ALLOCATE(ulat_phy(klon_omp,llm))
      ALLOCATE(q_phy(klon_omp,llm,nqtot))
      ALLOCATE(wflux_phy(klon_omp,llm))
      ALLOCATE(dulon_phy(klon_omp,llm))
      ALLOCATE(dulat_phy(klon_omp,llm))
      ALLOCATE(dTemp_phy(klon_omp,llm))
      ALLOCATE(dq_phy(klon_omp,llm,nqtot))
      ALLOCATE(dps_phy(klon_omp))
      ALLOCATE(plevmoy(llm+1))
      ALLOCATE(tmoy(llm))
!$OMP BARRIER
    ENDIF


!$OMP MASTER        
!    CALL update_calendar(it)
!$OMP END MASTER
!$OMP BARRIER
    dtphy=itau_physics*dt
    
    
    
    CALL transfert_message(f_u,req_u)
    
    DO ind=1,ndomain
      CALL swap_dimensions(ind)
      IF (assigned_domain(ind)) THEN
        CALL swap_geometry(ind)
      
        phis=f_phis(ind)
        ps=f_ps(ind)
        theta_rhodz=f_theta_rhodz(ind)
        u=f_u(ind)
        q=f_q(ind)
        wflux=f_wflux(ind)
        p=f_p(ind)
        pks=f_pks(ind)
        pk=f_pk(ind)
        p_layer=f_p_layer(ind)
        theta=f_theta(ind)
        phi=f_phi(ind)
        Temp=f_Temp(ind)
        ulon=f_ulon(ind)
        ulat=f_ulat(ind)
            
        CALL grid_icosa_to_physics

      ENDIF
    ENDDO
   
!$OMP BARRIER
!$OMP MASTER
    CALL SYSTEM_CLOCK(start_clock)
!$OMP END MASTER
    CALL trace_start("physic")
!    CALL trace_off()


!    CALL writeField("p_in",f_p)
!    CALL writeField("p_layer_in",f_p_layer)
!    CALL writeField("phi_in",f_phi)
!    CALL writeField("phis_in",f_phis)
!    CALL writeField("ulon_in",f_ulon)
!    CALL writeField("ulat_in",f_ulat)
!    CALL writeField("Temp_in",f_Temp)
!    CALL writeField("q_in",f_q)
!    CALL writeField("wflux_in",f_wflux)

!    CALL checksum(f_p)
!    CALL checksum(f_p_layer)
!    CALL checksum(f_phi)
!    CALL checksum(f_phis)
!    CALL checksum(f_ulon)
!    CALL checksum(f_ulat)
!    CALL checksum(f_Temp)
!    CALL checksum(f_q)
!    CALL checksum(f_wflux)

    CALL transfer_icosa_to_lmdz(f_p      , p_phy)
    CALL transfer_icosa_to_lmdz(f_p_layer, p_layer_phy)
    CALL transfer_icosa_to_lmdz(f_pk     , pk_phy)
    CALL transfer_icosa_to_lmdz(f_phi    , phi_phy)
    CALL transfer_icosa_to_lmdz(f_phis   , phis_phy )
    CALL transfer_icosa_to_lmdz(f_ulon   , ulon_phy )
    CALL transfer_icosa_to_lmdz(f_ulat   , ulat_phy)
    CALL transfer_icosa_to_lmdz(f_Temp   , Temp_phy)
    CALL transfer_icosa_to_lmdz(f_q      , q_phy)
    CALL transfer_icosa_to_lmdz(f_wflux  , wflux_phy)

    DO l=1,llm
      ! Warning: In the physics, vertical flux convention is positive if downwards!
      wflux_phy(:,l)= - wflux_phy(:,l)*cell_area(:)
      ! Compute relative geopotential
      phi_phy(:,l)=phi_phy(:,l)-phis_phy(:)
    ENDDO
    
    CALL wxios_set_context()

    ! Update pday and ptime to send to physics
!$OMP MASTER
    ptime=ptime+dtphy/day_length
    IF (ptime >= 1.) THEN
      ptime=ptime-1
      pday=pday+1
    ENDIF
!$OMP END MASTER
!$OMP BARRIER    

! NB: arguments plevmoy(:) planet-averaged mean pressure (Pa) at interlayers
!     and tmoy(:) planet-averaged mean temperature (K) at mid-layers
!     are an issue. For now we just set them to unphysical values
    plevmoy(:)=-999
    tmoy(:)=-999
    ! Ehouarn test: add some noise to ulon_phy and vlon_phy if they are zero
    IF(minval(ulon_phy).eq.maxval(ulon_phy)) ulon_phy(:,:)=1.e-15
    IF(minval(ulat_phy).eq.maxval(ulat_phy)) ulat_phy(:,:)=-1.e-15
    CALL physiq(klon_omp, llm, nqtot, &
                debut, lafin, &
                pday, ptime, dtphy, &
                p_phy, p_layer_phy, pk_phy, &
                phi_phy, phis_phy, presnivs, &
                ulon_phy, ulat_phy, Temp_phy, q_phy, &
                wflux_phy, plevmoy, tmoy, &
                dulon_phy, dulat_phy, dTemp_phy, dq_phy, dps_phy)
    
    CALL transfer_lmdz_to_icosa(dulon_phy, f_dulon )
    CALL transfer_lmdz_to_icosa(dulat_phy, f_dulat )
    CALL transfer_lmdz_to_icosa(dTemp_phy, f_dTemp )
    CALL transfer_lmdz_to_icosa(dq_phy   , f_dq )
    CALL transfer_lmdz_to_icosa(dps_phy  , f_dps )
 
!    CALL writeField("dulon_out",f_dulon)
!    CALL writeField("dulat_out",f_dulat)
!    CALL writeField("dTemp_out",f_dTemp)
!    CALL writeField("dq_out",f_dq)
!    CALL writeField("dps_out",f_dps)

!    CALL checksum(f_dulon)
!    CALL checksum(f_dulat)
!    CALL checksum(f_dTemp)
!    CALL checksum(f_dq)
!    CALL checksum(f_dps)
    
    CALL send_message(f_dps,req_dps0)
    CALL send_message(f_dulon,req_dulon0)
    CALL send_message(f_dulat,req_dulat0)
    CALL send_message(f_dTemp,req_dTemp0)
    CALL send_message(f_dq,req_dq0)

    CALL wait_message(req_dps0)
    CALL wait_message(req_dulon0)
    CALL wait_message(req_dulat0)
    CALL wait_message(req_dTemp0)
    CALL wait_message(req_dq0)


!    CALL trace_on()
    CALL trace_end("physic")
!$OMP MASTER
    CALL SYSTEM_CLOCK(stop_clock)
    count_clock=count_clock+stop_clock-start_clock
!$OMP END MASTER

!$OMP BARRIER                       

    DO ind=1,ndomain
      CALL swap_dimensions(ind)
      IF (assigned_domain(ind)) THEN
        CALL swap_geometry(ind)

        theta_rhodz=f_theta_rhodz(ind)
        u=f_u(ind)
        q=f_q(ind)
        ps=f_ps(ind)
        dulon=f_dulon(ind)
        dulat=f_dulat(ind)
        Temp=f_temp(ind)
        dTemp=f_dTemp(ind)
        dq=f_dq(ind)
        dps=f_dps(ind)
        duc=f_duc(ind)
        p=f_p(ind)
        pks=f_pks(ind)
        pk=f_pk(ind)
      
        CALL grid_physics_to_icosa
      ENDIF
    ENDDO

!$OMP BARRIER
    CALL xios_set_context    
   
 
  CONTAINS

    SUBROUTINE grid_icosa_to_physics
    USE pression_mod
    USE exner_mod
    USE theta2theta_rhodz_mod
    USE geopotential_mod
    USE wind_mod
    USE omp_para
    IMPLICIT NONE
    
    REAL(rstd) :: uc(3)
    INTEGER :: i,j,ij,l
    

! compute pression

      DO    l    = ll_begin,ll_endp1
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            p(ij,l) = ap(l) + bp(l) * ps(ij)
          ENDDO
        ENDDO
      ENDDO

!$OMP BARRIER

! compute exner
       
       IF (is_omp_first_level) THEN
         DO j=jj_begin,jj_end
            DO i=ii_begin,ii_end
               ij=(j-1)*iim+i
               pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
            ENDDO
         ENDDO
       ENDIF

       ! 3D : pk
       DO l = ll_begin,ll_end
          DO j=jj_begin,jj_end
             DO i=ii_begin,ii_end
                ij=(j-1)*iim+i
                pk(ij,l) = cpp * ((.5/preff)*(p(ij,l)+p(ij,l+1))) ** kappa
             ENDDO
          ENDDO
       ENDDO

!$OMP BARRIER

!   compute theta, temperature and pression at layer
    DO    l    = ll_begin, ll_end
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          ij=(j-1)*iim+i
          theta(ij,l) = theta_rhodz(ij,l,1) / ((p(ij,l)-p(ij,l+1))/g)
          Temp(ij,l) = theta(ij,l) * pk(ij,l) / cpp
          p_layer(ij,l)=preff*(pk(ij,l)/cpp)**(1./kappa) 
        ENDDO
      ENDDO
    ENDDO


!!! Compute geopotential
       
  ! for first layer
  IF (is_omp_first_level) THEN
    DO j=jj_begin,jj_end
      DO i=ii_begin,ii_end
        ij=(j-1)*iim+i
        phi( ij,1 ) = phis( ij ) + theta(ij,1) * ( pks(ij) - pk(ij,1) )
      ENDDO
    ENDDO
  ENDIF
!!-> implicit flush on phi(:,1)
          
  ! for other layers
  DO l = ll_beginp1, ll_end
    DO j=jj_begin,jj_end
      DO i=ii_begin,ii_end
        ij=(j-1)*iim+i
        phi(ij,l) =  0.5 * ( theta(ij,l)  + theta(ij,l-1) )  & 
                         * (  pk(ij,l-1) -  pk(ij,l)    )
      ENDDO
    ENDDO
  ENDDO       

!$OMP BARRIER


  IF (is_omp_first_level) THEN
    DO l = 2, llm
      DO j=jj_begin,jj_end
! ---> Bug compilo intel ici en openmp
! ---> Couper la boucle
       IF (j==jj_end+1) PRINT*,"this message must not be printed"
        DO i=ii_begin,ii_end
          ij=(j-1)*iim+i
          phi(ij,l) = phi(ij,l)+ phi(ij,l-1)
        ENDDO
      ENDDO
    ENDDO
! --> IMPLICIT FLUSH on phi --> non 
  ENDIF 

! compute wind centered lon lat compound
    DO l=ll_begin,ll_end
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          ij=(j-1)*iim+i
          uc(:)=1/Ai(ij)*                                                                                                &
                        ( ne(ij,right)*u(ij+u_right,l)*le(ij+u_right)*((xyz_v(ij+z_rdown,:)+xyz_v(ij+z_rup,:))/2-centroid(ij,:))  &
                         + ne(ij,rup)*u(ij+u_rup,l)*le(ij+u_rup)*((xyz_v(ij+z_rup,:)+xyz_v(ij+z_up,:))/2-centroid(ij,:))          &
                         + ne(ij,lup)*u(ij+u_lup,l)*le(ij+u_lup)*((xyz_v(ij+z_up,:)+xyz_v(ij+z_lup,:))/2-centroid(ij,:))          &
                         + ne(ij,left)*u(ij+u_left,l)*le(ij+u_left)*((xyz_v(ij+z_lup,:)+xyz_v(ij+z_ldown,:))/2-centroid(ij,:))    &
                         + ne(ij,ldown)*u(ij+u_ldown,l)*le(ij+u_ldown)*((xyz_v(ij+z_ldown,:)+xyz_v(ij+z_down,:))/2-centroid(ij,:))&
                         + ne(ij,rdown)*u(ij+u_rdown,l)*le(ij+u_rdown)*((xyz_v(ij+z_down,:)+xyz_v(ij+z_rdown,:))/2-centroid(ij,:)))
          ulon(ij,l)=sum(uc(:)*elon_i(ij,:))
          ulat(ij,l)=sum(uc(:)*elat_i(ij,:)) 
        ENDDO
      ENDDO
    ENDDO

!$OMP BARRIER
    END SUBROUTINE grid_icosa_to_physics


    SUBROUTINE grid_physics_to_icosa
    USE theta2theta_rhodz_mod
    USE omp_para
    IMPLICIT NONE
      INTEGER :: i,j,ij,l,iq
          
      DO l=ll_begin,ll_end
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            duc(ij,:,l)=dulon(ij,l)*elon_i(ij,:)+dulat(ij,l)*elat_i(ij,:)
          ENDDO
        ENDDO
      ENDDO

      DO l=ll_begin,ll_end
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            u(ij+u_right,l) = u(ij+u_right,l) + dtphy * sum( 0.5*(duc(ij,:,l) + duc(ij+t_right,:,l))*ep_e(ij+u_right,:) )
            u(ij+u_lup,l) = u(ij+u_lup,l) + dtphy * sum( 0.5*(duc(ij,:,l) + duc(ij+t_lup,:,l))*ep_e(ij+u_lup,:) )
            u(ij+u_ldown,l) = u(ij+u_ldown,l) + dtphy*sum( 0.5*(duc(ij,:,l) + duc(ij+t_ldown,:,l))*ep_e(ij+u_ldown,:) )
          ENDDO
        ENDDO
      ENDDO          

      DO l=ll_begin,ll_end
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            Temp(ij,l)=Temp(ij,l)+ dtphy * dTemp(ij,l)
          ENDDO
        ENDDO
      ENDDO          
      
      DO iq=1,nqtot
        DO l=ll_begin,ll_end
          DO j=jj_begin,jj_end
            DO i=ii_begin,ii_end
              ij=(j-1)*iim+i
              q(ij,l,iq)=q(ij,l,iq)+ dtphy * dq(ij,l,iq)
            ENDDO
          ENDDO
        ENDDO 
      ENDDO

!$OMP BARRIER
      
     IF (is_omp_first_level) THEN
       DO j=jj_begin,jj_end
         DO i=ii_begin,ii_end
           ij=(j-1)*iim+i
           ps(ij)=ps(ij)+ dtphy * dps(ij)
          ENDDO
       ENDDO
     ENDIF

!     CALL compute_temperature2theta_rhodz(ps,Temp,theta_rhodz,0)

! compute pression
!$OMP BARRIER
      DO    l    = ll_begin,ll_endp1
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            p(ij,l) = ap(l) + bp(l) * ps(ij)
          ENDDO
        ENDDO
      ENDDO

!$OMP BARRIER

! compute exner
       
       IF (is_omp_first_level) THEN
         DO j=jj_begin,jj_end
            DO i=ii_begin,ii_end
               ij=(j-1)*iim+i
               pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
            ENDDO
         ENDDO
       ENDIF

       ! 3D : pk
       DO l = ll_begin,ll_end
          DO j=jj_begin,jj_end
             DO i=ii_begin,ii_end
                ij=(j-1)*iim+i
                pk(ij,l) = cpp * ((.5/preff)*(p(ij,l)+p(ij,l+1))) ** kappa
             ENDDO
          ENDDO
       ENDDO

!$OMP BARRIER

!   compute theta, temperature and pression at layer
    DO    l    = ll_begin, ll_end
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          ij=(j-1)*iim+i
          theta_rhodz(ij,l,1) = temp(ij,l) * ((p(ij,l)-p(ij,l+1))/g) / (pk(ij,l) / cpp )
        ENDDO
      ENDDO
    ENDDO
    
    END SUBROUTINE grid_physics_to_icosa



  END SUBROUTINE physics





END MODULE interface_icosa_lmdz_mod
