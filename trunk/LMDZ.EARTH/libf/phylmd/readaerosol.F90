! $Id: readaerosol.F90 1484 2011-02-09 15:44:57Z jghattas $
!
MODULE readaerosol_mod

  REAL, SAVE :: not_valid=-333.

CONTAINS

SUBROUTINE readaerosol(name_aero, type, filename, iyr_in, klev_src, pt_ap, pt_b, pt_out, psurf, load)

!****************************************************************************************
! This routine will read the aersosol from file. 
!
! Read a year data with get_aero_fromfile depending on aer_type : 
! - actuel   : read year 1980
! - preind   : read natural data
! - scenario : read one or two years and do eventually linare time interpolation
!
! Return pointer, pt_out, to the year read or result from interpolation
!****************************************************************************************
  USE dimphy

  IMPLICIT NONE

 INCLUDE "iniprint.h"

  ! Input arguments
  CHARACTER(len=7), INTENT(IN) :: name_aero
  CHARACTER(len=*), INTENT(IN) :: type  ! actuel, annuel, scenario or preind
  CHARACTER(len=8), INTENT(IN) :: filename
  INTEGER, INTENT(IN)          :: iyr_in

  ! Output
  INTEGER, INTENT(OUT)            :: klev_src
  REAL, POINTER, DIMENSION(:)     :: pt_ap        ! Pointer for describing the vertical levels      
  REAL, POINTER, DIMENSION(:)     :: pt_b         ! Pointer for describing the vertical levels      
  REAL, POINTER, DIMENSION(:,:,:) :: pt_out       ! The massvar distributions, DIMENSION(klon, klev_src, 12)
  REAL, DIMENSION(klon,12), INTENT(OUT) :: psurf  ! Surface pression for 12 months
  REAL, DIMENSION(klon,12), INTENT(OUT) :: load   ! Aerosol mass load in each column for 12 months

  ! Local variables
  CHARACTER(len=4)                :: cyear
  REAL, POINTER, DIMENSION(:,:,:) :: pt_2
  REAL, DIMENSION(klon,12)        :: psurf2, load2
  REAL                            :: p0           ! Reference pressure
  INTEGER                         :: iyr1, iyr2, klev_src2
  INTEGER                         :: it, k, i
  LOGICAL, PARAMETER              :: lonlyone=.FALSE.

!****************************************************************************************
! Read data depending on aer_type
!
!****************************************************************************************

  IF (type == 'actuel') THEN
! Read and return data for year 1980
!****************************************************************************************
     cyear='1980'
     ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
     

  ELSE IF (type == 'preind') THEN
! Read and return data from file with suffix .nat
!****************************************************************************************     
     cyear='.nat'
     ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
     
  ELSE IF (type == 'annuel') THEN
! Read and return data from scenario annual files
!****************************************************************************************     
     WRITE(cyear,'(I4)') iyr_in
     WRITE(lunout,*) 'get_aero 3 iyr_in=', iyr_in,'   ',cyear
     ! get_aero_fromfile returns pt_out allocated and initialized with data for nbr_tsteps month 
     ! pt_out has dimensions (klon, klev_src, 12)
     CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
     
  ELSE IF (type == 'scenario') THEN
! Read data depending on actual year and interpolate if necessary
!****************************************************************************************
     IF (iyr_in .LT. 1850) THEN
        cyear='.nat'
        WRITE(lunout,*) 'get_aero 1 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
        
     ELSE IF (iyr_in .GE. 2100) THEN
        cyear='2100'
        WRITE(lunout,*) 'get_aero 2 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
        
     ELSE
        ! Read data from 2 decades and interpolate to actual year
        ! a) from actual 10-yr-period
        IF (iyr_in.LT.1900) THEN
           iyr1 = 1850
           iyr2 = 1900
        ELSE IF (iyr_in.GE.1900.AND.iyr_in.LT.1920) THEN
           iyr1 = 1900
           iyr2 = 1920
        ELSE 
           iyr1 = INT(iyr_in/10)*10
           iyr2 = INT(1+iyr_in/10)*10
        ENDIF
        
        WRITE(cyear,'(I4)') iyr1
        WRITE(lunout,*) 'get_aero 3 iyr_in=', iyr_in,'   ',cyear
        ! get_aero_fromfile returns pt_out allocated and initialized with data for 12 month 
        ! pt_out has dimensions (klon, klev_src, 12)
        CALL get_aero_fromfile(name_aero, cyear, filename, klev_src, pt_ap, pt_b, p0, pt_out, psurf, load)
        
        ! If to read two decades:
        IF (.NOT.lonlyone) THEN 
           
           ! b) from the next following one
           WRITE(cyear,'(I4)') iyr2
           WRITE(lunout,*) 'get_aero 4 iyr_in=', iyr_in,'   ',cyear
           
           NULLIFY(pt_2)
           ! get_aero_fromfile returns pt_2 allocated and initialized with data for 12 month 
           ! pt_2 has dimensions (klon, klev_src, 12)
           CALL get_aero_fromfile(name_aero, cyear, filename, klev_src2, pt_ap, pt_b, p0, pt_2, psurf2, load2)
           ! Test for same number of vertical levels
           IF (klev_src /= klev_src2) THEN
              WRITE(lunout,*) 'Two aerosols files with different number of vertical levels is not allowded'
              CALL abort_gcm('readaersosol','Error in number of vertical levels',1)
           END IF
           
           ! Linare interpolate to the actual year:
           DO it=1,12
              DO k=1,klev_src
                 DO i = 1, klon
                    pt_out(i,k,it) = &
                         pt_out(i,k,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                         (pt_out(i,k,it) - pt_2(i,k,it))
                 END DO
              END DO

              DO i = 1, klon
                 psurf(i,it) = &
                      psurf(i,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                      (psurf(i,it) - psurf2(i,it))

                 load(i,it) = &
                      load(i,it) - REAL(iyr_in-iyr1)/REAL(iyr2-iyr1) * &
                      (load(i,it) - load2(i,it))
              END DO
           END DO

           ! Deallocate pt_2 no more needed
           DEALLOCATE(pt_2)
           
        END IF ! lonlyone
     END IF ! iyr_in .LT. 1850

  ELSE
     WRITE(lunout,*)'This option is not implemented : aer_type = ', type,' name_aero=',name_aero
     CALL abort_gcm('readaerosol','Error : aer_type parameter not accepted',1)
  END IF ! type


END SUBROUTINE readaerosol


  SUBROUTINE get_aero_fromfile(varname, cyr, filename, klev_src, pt_ap, pt_b, p0, pt_year, psurf_out, load_out)
!****************************************************************************************
! Read 12 month aerosol from file and distribute to local process on physical grid. 
! Vertical levels, klev_src, may differ from model levels if new file format.
!
! For mpi_root and master thread :
! 1) Open file 
! 2) Find vertical dimension klev_src
! 3) Read field month by month
! 4) Close file  
! 5) Transform the global field from 2D(iim, jjp+1) to 1D(klon_glo)
!     - Also the levels and the latitudes have to be inversed
!
! For all processes and threads :
! 6) Scatter global field(klon_glo) to local process domain(klon)
! 7) Test for negative values
!****************************************************************************************

    USE netcdf
    USE dimphy
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
    USE iophy, ONLY : io_lon, io_lat

    IMPLICIT NONE
      
    INCLUDE "dimensions.h"      
    INCLUDE "iniprint.h"

! Input argumets
    CHARACTER(len=7), INTENT(IN)          :: varname
    CHARACTER(len=4), INTENT(IN)          :: cyr
    CHARACTER(len=8), INTENT(IN)          :: filename

! Output arguments
    INTEGER, INTENT(OUT)                  :: klev_src     ! Number of vertical levels in file
    REAL, POINTER, DIMENSION(:)           :: pt_ap        ! Pointer for describing the vertical levels      
    REAL, POINTER, DIMENSION(:)           :: pt_b         ! Pointer for describing the vertical levels      
    REAL                                  :: p0           ! Reference pressure value
    REAL, POINTER, DIMENSION(:,:,:)       :: pt_year      ! Pointer-variabale from file, 12 month, grid : klon,klev_src
    REAL, DIMENSION(klon,12), INTENT(OUT) :: psurf_out    ! Surface pression for 12 months
    REAL, DIMENSION(klon,12), INTENT(OUT) :: load_out     ! Aerosol mass load in each column
    INTEGER                               :: nbr_tsteps   ! number of month in file read

! Local variables
    CHARACTER(len=30)     :: fname
    CHARACTER(len=30)     :: cvar
    INTEGER               :: ncid, dimid, varid
    INTEGER               :: imth, i, j, k, ierr
    REAL                  :: npole, spole
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: varmth
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: varyear       ! Global variable read from file, 12 month
    REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: varyear_glo1D !(klon_glo, klev_src, 12)
    REAL, ALLOCATABLE, DIMENSION(:)       :: varktmp

    REAL, DIMENSION(iim,jjm+1,12)         :: psurf_glo2D   ! Surface pression for 12 months on dynamics global grid
    REAL, DIMENSION(klon_glo,12)          :: psurf_glo1D   ! -"- on physical global grid
    REAL, DIMENSION(iim,jjm+1,12)         :: load_glo2D    ! Load for 12 months on dynamics global grid
    REAL, DIMENSION(klon_glo,12)          :: load_glo1D    ! -"- on physical global grid
    REAL, DIMENSION(iim,jjm+1)            :: vartmp
    REAL, DIMENSION(iim)                  :: lon_src              ! longitudes in file
    REAL, DIMENSION(jjm+1)                :: lat_src, lat_src_inv ! latitudes in file
    LOGICAL                               :: new_file             ! true if new file format detected
    LOGICAL                               :: invert_lat           ! true if the field has to be inverted for latitudes


    ! Deallocate pointers
    IF (ASSOCIATED(pt_ap)) DEALLOCATE(pt_ap)
    IF (ASSOCIATED(pt_b))  DEALLOCATE(pt_b)

    IF (is_mpi_root .AND. is_omp_root) THEN

! 1) Open file 
!****************************************************************************************
! Add suffix to filename
       fname = trim(filename)//cyr//'.nc'
  
       WRITE(lunout,*) 'reading variable ',TRIM(varname),' in file ', TRIM(fname)
       CALL check_err( nf90_open(TRIM(fname), NF90_NOWRITE, ncid) )

! Test for equal longitudes and latitudes in file and model
!****************************************************************************************
       ! Read and test longitudes
       CALL check_err( nf90_inq_varid(ncid, 'lon', varid) )
       CALL check_err( nf90_get_var(ncid, varid, lon_src(:)) )
       
       IF (maxval(ABS(lon_src - io_lon)) > 0.001) THEN
          WRITE(lunout,*) 'Problem in longitudes read from file : ',TRIM(fname)
          WRITE(lunout,*) 'longitudes in file ', TRIM(fname),' : ', lon_src
          WRITE(lunout,*) 'longitudes in model :', io_lon
          
          CALL abort_gcm('get_aero_fromfile', 'longitudes are not the same in file and model',1)
       END IF

       ! Read and test latitudes
       CALL check_err( nf90_inq_varid(ncid, 'lat', varid) )
       CALL check_err( nf90_get_var(ncid, varid, lat_src(:)) )

       ! Invert source latitudes
       DO j = 1, jjm+1
          lat_src_inv(j) = lat_src(jjm+1 +1 -j)
       END DO

       IF (maxval(ABS(lat_src - io_lat)) < 0.001) THEN
          ! Latitudes are the same
          invert_lat=.FALSE.
       ELSE IF (maxval(ABS(lat_src_inv - io_lat)) < 0.001) THEN
          ! Inverted source latitudes correspond to model latitudes
          WRITE(lunout,*) 'latitudes will be inverted for file : ',TRIM(fname)
          invert_lat=.TRUE.
       ELSE
          WRITE(lunout,*) 'Problem in latitudes read from file : ',TRIM(fname)
          WRITE(lunout,*) 'latitudes in file ', TRIM(fname),' : ', lat_src      
          WRITE(lunout,*) 'latitudes in model :', io_lat
          CALL abort_gcm('get_aero_fromfile', 'latitudes do not correspond between file and model',1)
       END IF


! 2) Check if old or new file is avalabale.
!    New type of file should contain the dimension 'lev'
!    Old type of file should contain the dimension 'PRESNIVS'
!****************************************************************************************
       ierr = nf90_inq_dimid(ncid, 'lev', dimid) 
       IF (ierr /= NF90_NOERR) THEN
          ! Coordinate axe lev not found. Check for presnivs.
          ierr = nf90_inq_dimid(ncid, 'PRESNIVS', dimid)
          IF (ierr /= NF90_NOERR) THEN
             ! Dimension PRESNIVS not found either
             CALL abort_gcm('get_aero_fromfile', 'dimension lev or presnivs not in file',1)
          ELSE 
             ! Old file found
             new_file=.FALSE.
             WRITE(lunout,*) 'Vertical interpolation for ',TRIM(varname),' will not be done'
          END IF
       ELSE
          ! New file found
          new_file=.TRUE.
          WRITE(lunout,*) 'Vertical interpolation for ',TRIM(varname),' will be done'
       END IF
       
! 2) Find vertical dimension klev_src
!****************************************************************************************
       CALL check_err( nf90_inquire_dimension(ncid, dimid, len = klev_src) )
       
     ! Allocate variables depending on the number of vertical levels
       ALLOCATE(varmth(iim, jjm+1, klev_src), varyear(iim, jjm+1, klev_src, 12), stat=ierr)
       IF (ierr /= 0) CALL abort_gcm('get_aero_fromfile', 'pb in allocation 1',1)

       ALLOCATE(pt_ap(klev_src), pt_b(klev_src), varktmp(klev_src), stat=ierr)
       IF (ierr /= 0) CALL abort_gcm('get_aero_fromfile', 'pb in allocation 2',1)

! 3) Read all variables from file
!    There is 2 options for the file structure :
!    new_file=TRUE  : read varyear, ps, pt_ap and pt_b
!    new_file=FALSE : read varyear month by month
!****************************************************************************************

       IF (new_file) THEN
! ++) Check number of month in file opened
!**************************************************************************************************
       ierr = nf90_inq_dimid(ncid, 'TIME',dimid)
       CALL check_err( nf90_inquire_dimension(ncid, dimid, len = nbr_tsteps) )
!       IF (nbr_tsteps /= 12 .AND. nbr_tsteps /= 14) THEN
       IF (nbr_tsteps /= 12 ) THEN
         CALL abort_gcm('get_aero_fromfile', 'not the right number of months in aerosol file read (should be 12 for the moment)',1)
       ENDIF

! ++) Read the aerosol concentration month by month and concatenate to total variable varyear
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, TRIM(varname), varid) )
          
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, varyear(:,:,:,:)) )
          
! ++) Read surface pression, 12 month in one variable
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "ps", varid) )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, psurf_glo2D) )
          
! ++) Read mass load, 12 month in one variable
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "load_"//TRIM(varname), varid) )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, load_glo2D) )
          
! ++) Read ap
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "ap", varid) )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, pt_ap) )

! ++) Read b
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "b", varid) )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, pt_b) )

! ++) Read p0 : reference pressure
!****************************************************************************************
          ! Get variable id
          CALL check_err( nf90_inq_varid(ncid, "p0", varid) )
          ! Get the variable
          CALL check_err( nf90_get_var(ncid, varid, p0) )
          

       ELSE  ! old file

! ++) Read the aerosol concentration month by month and concatenate to total variable varyear
!****************************************************************************************
          DO imth=1, 12
             IF (imth.EQ.1) THEN
                cvar=TRIM(varname)//'JAN'
             ELSE IF (imth.EQ.2) THEN
                cvar=TRIM(varname)//'FEB'
             ELSE IF (imth.EQ.3) THEN
                cvar=TRIM(varname)//'MAR'
             ELSE IF (imth.EQ.4) THEN
                cvar=TRIM(varname)//'APR'
             ELSE IF (imth.EQ.5) THEN
                cvar=TRIM(varname)//'MAY'
             ELSE IF (imth.EQ.6) THEN
                cvar=TRIM(varname)//'JUN'
             ELSE IF (imth.EQ.7) THEN
                cvar=TRIM(varname)//'JUL'
             ELSE IF (imth.EQ.8) THEN
                cvar=TRIM(varname)//'AUG'
             ELSE IF (imth.EQ.9) THEN
                cvar=TRIM(varname)//'SEP'
             ELSE IF (imth.EQ.10) THEN
                cvar=TRIM(varname)//'OCT'
             ELSE IF (imth.EQ.11) THEN
                cvar=TRIM(varname)//'NOV'
             ELSE IF (imth.EQ.12) THEN
                cvar=TRIM(varname)//'DEC'
             END IF
             
             ! Get variable id
             CALL check_err( nf90_inq_varid(ncid, TRIM(cvar), varid) )
             
             ! Get the variable
             CALL check_err( nf90_get_var(ncid, varid, varmth) )
             
             ! Store in variable for the whole year
             varyear(:,:,:,imth)=varmth(:,:,:)
             
          END DO
          
          ! Putting dummy 
          psurf_glo2D(:,:,:) = not_valid
          load_glo2D(:,:,:)  = not_valid
          pt_ap(:) = not_valid
          pt_b(:)  = not_valid

       END IF

! 4) Close file  
!****************************************************************************************
       CALL check_err( nf90_close(ncid) )
     

! 5) Transform the global field from 2D(iim, jjp+1) to 1D(klon_glo)
!****************************************************************************************
! Test if vertical levels have to be inversed

       IF ((pt_b(1) < pt_b(klev_src)) .OR. .NOT. new_file) THEN
!          WRITE(lunout,*) 'Vertical axis in file ',TRIM(fname), ' needs to be inverted'
!          WRITE(lunout,*) 'before pt_ap = ', pt_ap
!          WRITE(lunout,*) 'before pt_b = ', pt_b
          
          ! Inverse vertical levels for varyear 
          DO imth=1, 12
             varmth(:,:,:) = varyear(:,:,:,imth) ! use varmth temporarly
             DO k=1, klev_src
                DO j=1, jjm+1
                   DO i=1,iim 
                      varyear(i,j,k,imth) = varmth(i,j,klev_src+1-k)
                   END DO
                END DO
             END DO
          END DO
           
          ! Inverte vertical axes for pt_ap and pt_b
          varktmp(:) = pt_ap(:)
          DO k=1, klev_src
             pt_ap(k) = varktmp(klev_src+1-k)
          END DO

          varktmp(:) = pt_b(:)
          DO k=1, klev_src
             pt_b(k) = varktmp(klev_src+1-k)
          END DO
          WRITE(lunout,*) 'after pt_ap = ', pt_ap
          WRITE(lunout,*) 'after pt_b = ', pt_b

       ELSE 
          WRITE(lunout,*) 'Vertical axis in file ',TRIM(fname), ' is ok, no vertical inversion is done'
          WRITE(lunout,*) 'pt_ap = ', pt_ap
          WRITE(lunout,*) 'pt_b = ', pt_b
       END IF

!     - Invert latitudes if necessary
       DO imth=1, 12
          IF (invert_lat) THEN

             ! Invert latitudes for the variable
             varmth(:,:,:) = varyear(:,:,:,imth) ! use varmth temporarly
             DO k=1,klev_src
                DO j=1,jjm+1
                   DO i=1,iim
                      varyear(i,j,k,imth) = varmth(i,jjm+1+1-j,k)
                   END DO
                END DO
             END DO
             
             ! Invert latitudes for surface pressure
             vartmp(:,:) = psurf_glo2D(:,:,imth)
             DO j=1, jjm+1
                DO i=1,iim
                   psurf_glo2D(i,j,imth)= vartmp(i,jjm+1+1-j)
                END DO
             END DO
             
             ! Invert latitudes for the load
             vartmp(:,:) = load_glo2D(:,:,imth)
             DO j=1, jjm+1
                DO i=1,iim
                   load_glo2D(i,j,imth)= vartmp(i,jjm+1+1-j)
                END DO
             END DO
          END IF ! invert_lat
             
          ! Do zonal mead at poles and distribut at whole first and last latitude
          DO k=1, klev_src
             npole=0.  ! North pole, j=1
             spole=0.  ! South pole, j=jjm+1         
             DO i=1,iim
                npole = npole + varyear(i,1,k,imth)
                spole = spole + varyear(i,jjm+1,k,imth)
             END DO
             npole = npole/REAL(iim)
             spole = spole/REAL(iim)
             varyear(:,1,    k,imth) = npole
             varyear(:,jjm+1,k,imth) = spole
          END DO
       END DO ! imth
       
       ALLOCATE(varyear_glo1D(klon_glo, klev_src, 12), stat=ierr)
       IF (ierr /= 0) CALL abort_gcm('get_aero_fromfile', 'pb in allocation 3',1)
       
       ! Transform from 2D to 1D field
       CALL grid2Dto1D_glo(varyear,varyear_glo1D)
       CALL grid2Dto1D_glo(psurf_glo2D,psurf_glo1D)
       CALL grid2Dto1D_glo(load_glo2D,load_glo1D)

    ELSE
      ALLOCATE(varyear_glo1D(0,0,0))        
    END IF ! is_mpi_root .AND. is_omp_root

!$OMP BARRIER
  
! 6) Distribute to all processes
!    Scatter global field(klon_glo) to local process domain(klon)
!    and distribute klev_src to all processes
!****************************************************************************************

    ! Distribute klev_src
    CALL bcast(klev_src)

    ! Allocate and distribute pt_ap and pt_b
    IF (.NOT. ASSOCIATED(pt_ap)) THEN  ! if pt_ap is allocated also pt_b is allocated
       ALLOCATE(pt_ap(klev_src), pt_b(klev_src), stat=ierr)
       IF (ierr /= 0) CALL abort_gcm('get_aero_fromfile', 'pb in allocation 4',1)
    END IF
    CALL bcast(pt_ap)
    CALL bcast(pt_b)

    ! Allocate space for output pointer variable at local process
    IF (ASSOCIATED(pt_year)) DEALLOCATE(pt_year)
    ALLOCATE(pt_year(klon, klev_src, 12), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('get_aero_fromfile', 'pb in allocation 5',1)

    ! Scatter global field to local domain at local process
    CALL scatter(varyear_glo1D, pt_year)
    CALL scatter(psurf_glo1D, psurf_out)
    CALL scatter(load_glo1D,  load_out)

! 7) Test for negative values
!****************************************************************************************
    IF (MINVAL(pt_year) < 0.) THEN
       WRITE(lunout,*) 'Warning! Negative values read from file :', fname
    END IF

  END SUBROUTINE get_aero_fromfile


  SUBROUTINE check_err(status)
    USE netcdf
    IMPLICIT NONE

    INCLUDE "iniprint.h"
    INTEGER, INTENT (IN) :: status

    IF (status /= NF90_NOERR) THEN
       WRITE(lunout,*) 'Error in get_aero_fromfile ',status
       CALL abort_gcm('get_aero_fromfile',trim(nf90_strerror(status)),1)
    END IF

  END SUBROUTINE check_err


END MODULE readaerosol_mod
