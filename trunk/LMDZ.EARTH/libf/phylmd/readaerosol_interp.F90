! $Id$
!
SUBROUTINE readaerosol_interp(id_aero, itap, pdtphys, r_day, first, pplay, paprs, t_seri, mass_out, pi_mass_out, load_src)
!
! This routine will return the mass concentration at actual day(mass_out) and 
! the pre-industrial values(pi_mass_out) for aerosol corresponding to "id_aero".
! The mass concentrations for all aerosols are saved in this routine but each
! call to this routine only treats the aerosol "id_aero".
!
! 1) Read in data for the whole year, only at first time step
! 2) Interpolate to the actual day, only at new day
! 3) Interpolate to the model vertical grid (target grid), only at new day
! 4) Test for negative mass values

  USE ioipsl
  USE dimphy, ONLY : klev,klon
  USE mod_phys_lmdz_para, ONLY : mpi_rank  
  USE readaerosol_mod
  USE aero_mod, ONLY : naero_spc, name_aero
  USE write_field_phy
  USE phys_cal_mod

  IMPLICIT NONE

  INCLUDE "YOMCST.h"
  INCLUDE "chem.h"      
  INCLUDE "temps.h"      
  INCLUDE "clesphys.h"
  INCLUDE "iniprint.h"
  INCLUDE "dimensions.h"
  INCLUDE "comvert.h"
!
! Input:
!****************************************************************************************
  INTEGER, INTENT(IN)                    :: id_aero! Identity number for the aerosol to treat
  INTEGER, INTENT(IN)                    :: itap   ! Physic step count
  REAL, INTENT(IN)                       :: pdtphys! Physic day step
  REAL, INTENT(IN)                       :: r_day  ! Day of integration
  LOGICAL, INTENT(IN)                    :: first  ! First model timestep 
  REAL, DIMENSION(klon,klev), INTENT(IN) :: pplay  ! pression at model mid-layers
  REAL, DIMENSION(klon,klev+1),INTENT(IN):: paprs  ! pression between model layers
  REAL, DIMENSION(klon,klev), INTENT(IN) :: t_seri ! air temperature
!      
! Output:      
!****************************************************************************************
  REAL, INTENT(OUT) :: mass_out(klon,klev)    ! Mass of aerosol (monthly mean data,from file) [ug AIBCM/m3]
  REAL, INTENT(OUT) :: pi_mass_out(klon,klev) ! Mass of preindustrial aerosol (monthly mean data,from file) [ug AIBCM/m3]
  REAL, INTENT(OUT) :: load_src(klon) ! Load of aerosol (monthly mean data,from file) [kg/m3]
!      
! Local Variables:
!****************************************************************************************
  INTEGER                         :: i, k, ierr
  INTEGER                         :: iday, iyr, lmt_pas
!  INTEGER                         :: im, day1, day2, im2
  INTEGER                         :: im, im2
  REAL                            :: day1, day2
  INTEGER                         :: pi_klev_src ! Only for testing purpose
  INTEGER, SAVE                   :: klev_src    ! Number of vertical levles in source field
!$OMP THREADPRIVATE(klev_src)

  REAL                              :: zrho      ! Air density [kg/m3]
  REAL                              :: volm      ! Volyme de melange [kg/kg]
  REAL, DIMENSION(klon)             :: psurf_day, pi_psurf_day
  REAL, DIMENSION(klon)             :: pi_load_src  ! Mass load at source grid
  REAL, DIMENSION(klon)             :: load_tgt, load_tgt_test
  REAL, DIMENSION(klon,klev)        :: delp ! pressure difference in each model layer

  REAL, ALLOCATABLE, DIMENSION(:,:)            :: pplay_src ! pression mid-layer at source levels
  REAL, ALLOCATABLE, DIMENSION(:,:)            :: tmp1, tmp2  ! Temporary variables
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:), SAVE  :: var_year    ! VAR in right dimension for the total year
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:), SAVE  :: pi_var_year ! pre-industrial VAR, -"-
!$OMP THREADPRIVATE(var_year,pi_var_year)
  REAL, ALLOCATABLE, DIMENSION(:,:,:),SAVE     :: var_day     ! VAR interpolated to the actual day and model grid
  REAL, ALLOCATABLE, DIMENSION(:,:,:),SAVE     :: pi_var_day  ! pre-industrial VAR, -"-
!$OMP THREADPRIVATE(var_day,pi_var_day)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE    :: psurf_year, pi_psurf_year ! surface pressure for the total year
!$OMP THREADPRIVATE(psurf_year, pi_psurf_year)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE    :: load_year, pi_load_year   ! load in the column for the total year
!$OMP THREADPRIVATE(load_year, pi_load_year)

  REAL, DIMENSION(:,:,:), POINTER   :: pt_tmp      ! Pointer allocated in readaerosol
  REAL, POINTER, DIMENSION(:), SAVE :: pt_ap, pt_b ! Pointer for describing the vertical levels 
!$OMP THREADPRIVATE(pt_ap, pt_b)
  INTEGER, SAVE                     :: nbr_tsteps ! number of time steps in file read
  REAL, DIMENSION(14), SAVE         :: month_len, month_start, month_mid
!$OMP THREADPRIVATE(nbr_tsteps, month_len, month_start, month_mid)
  REAL                              :: jDay

  LOGICAL            :: lnewday      ! Indicates if first time step at a new day
  LOGICAL            :: OLDNEWDAY
  LOGICAL,SAVE       :: vert_interp  ! Indicates if vertical interpolation will be done
  LOGICAL,SAVE       :: debug=.FALSE.! Debugging in this subroutine
!$OMP THREADPRIVATE(vert_interp, debug)
  CHARACTER(len=8)      :: type
  CHARACTER(len=8)      :: filename


!****************************************************************************************
! Initialization
!
!****************************************************************************************

! Calculation to find if it is a new day

  IF(mpi_rank == 0 .AND. debug )then
     PRINT*,'CONTROL PANEL REGARDING TIME STEPING'
  ENDIF

  ! Use phys_cal_mod
  iday= day_cur
  iyr = year_cur
  im  = mth_cur

!  iday = INT(r_day)
!  iyr  = iday/360
!  iday = iday-iyr*360         ! day of the actual year
!  iyr  = iyr + annee_ref      ! year of the run   
!  im   = iday/30 +1           ! the actual month
  CALL ymds2ju(iyr, im, iday, 0., jDay)
!   CALL ymds2ju(iyr, im, iday-(im-1)*30, 0., jDay)


  IF(MOD(itap-1,NINT(86400./pdtphys)) == 0)THEN
     lnewday=.TRUE.
  ELSE
     lnewday=.FALSE.
  ENDIF

  IF(mpi_rank == 0 .AND. debug)then
     ! 0.02 is about 0.5/24, namly less than half an hour
     OLDNEWDAY = (r_day-REAL(iday) < 0.02)
     ! Once per day, update aerosol fields
     lmt_pas = NINT(86400./pdtphys)
     PRINT*,'r_day-REAL(iday) =',r_day-REAL(iday) 
     PRINT*,'itap =',itap
     PRINT*,'pdtphys =',pdtphys
     PRINT*,'lmt_pas =',lmt_pas
     PRINT*,'iday =',iday
     PRINT*,'r_day =',r_day
     PRINT*,'day_cur =',day_cur
     PRINT*,'mth_cur =',mth_cur
     PRINT*,'year_cur =',year_cur
     PRINT*,'NINT(86400./pdtphys) =',NINT(86400./pdtphys)
     PRINT*,'MOD(0,1) =',MOD(0,1)
     PRINT*,'lnewday =',lnewday
     PRINT*,'OLDNEWDAY =',OLDNEWDAY
  ENDIF

  IF (.NOT. ALLOCATED(var_day)) THEN
     ALLOCATE( var_day(klon, klev, naero_spc), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 1',1)
     ALLOCATE( pi_var_day(klon, klev, naero_spc), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 2',1)

     ALLOCATE( psurf_year(klon, 12, naero_spc), pi_psurf_year(klon, 12, naero_spc), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 3',1)

     ALLOCATE( load_year(klon, 12, naero_spc), pi_load_year(klon, 12, naero_spc), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 4',1)

     lnewday=.TRUE.

     NULLIFY(pt_ap)
     NULLIFY(pt_b)
  END IF

!****************************************************************************************
! 1) Read in data : corresponding to the actual year and preindustrial data. 
!    Only for the first day of the year.
!
!****************************************************************************************
  IF ( (first .OR. iday==0) .AND. lnewday ) THEN 
     NULLIFY(pt_tmp)

     ! Reading values corresponding to the closest year taking into count the choice of aer_type. 
     ! For aer_type=scenario interpolation between 2 data sets is done in readaerosol.
     ! If aer_type=mix1 or mix2, the run type and file name depends on the aerosol.
     IF (aer_type=='preind' .OR. aer_type=='actuel' .OR. aer_type=='annuel' .OR. aer_type=='scenario') THEN
        ! Standard case
        filename='aerosols'
        type=aer_type
     ELSE IF (aer_type == 'mix1') THEN
        ! Special case using a mix of decenal sulfate file and annual aerosols(all aerosols except sulfate)
        IF (name_aero(id_aero) == 'SO4') THEN
           filename='so4.run '
           type='scenario'
        ELSE
           filename='aerosols'
           type='annuel'
        END IF
     ELSE  IF (aer_type == 'mix2') THEN
        ! Special case using a mix of decenal sulfate file and natrual aerosols
        IF (name_aero(id_aero) == 'SO4') THEN
           filename='so4.run '
           type='scenario'
        ELSE
           filename='aerosols'
           type='preind'
        END IF
     ELSE
        CALL abort_gcm('readaerosol_interp', 'this aer_type not supported',1)
     END IF

     CALL readaerosol(name_aero(id_aero), type, filename, iyr, klev_src, pt_ap, pt_b, pt_tmp, &
          psurf_year(:,:,id_aero), load_year(:,:,id_aero))
     IF (.NOT. ALLOCATED(var_year)) THEN
        ALLOCATE(var_year(klon, klev_src, 12, naero_spc), stat=ierr)
        IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 5',1)
     END IF
     var_year(:,:,:,id_aero) = pt_tmp(:,:,:)

     ! Reading values corresponding to the preindustrial concentrations.
     type='preind'
     CALL readaerosol(name_aero(id_aero), type, filename, iyr, pi_klev_src, pt_ap, pt_b, pt_tmp, &
          pi_psurf_year(:,:,id_aero), pi_load_year(:,:,id_aero))

     ! klev_src must be the same in both files. 
     ! Also supposing pt_ap and pt_b to be the same in the 2 files without testing. 
     IF (pi_klev_src /= klev_src) THEN
        WRITE(lunout,*) 'Error! All forcing files for the same aerosol must have the same vertical dimension'
        WRITE(lunout,*) 'Aerosol : ', name_aero(id_aero)
        CALL abort_gcm('readaerosol_interp','Differnt vertical axes in aerosol forcing files',1)
     END IF

     IF (.NOT. ALLOCATED(pi_var_year)) THEN
        ALLOCATE(pi_var_year(klon, klev_src, 12, naero_spc), stat=ierr)
        IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 6',1)
     END IF
     pi_var_year(:,:,:,id_aero) = pt_tmp(:,:,:)
    
     IF (debug) THEN
        CALL writefield_phy('var_year_jan',var_year(:,:,1,id_aero),klev_src)
        CALL writefield_phy('var_year_dec',var_year(:,:,12,id_aero),klev_src)
        CALL writefield_phy('psurf_src',psurf_year(:,:,id_aero),1)
        CALL writefield_phy('pi_psurf_src',pi_psurf_year(:,:,id_aero),1)
        CALL writefield_phy('load_year_src',load_year(:,:,id_aero),1)
        CALL writefield_phy('pi_load_year_src',pi_load_year(:,:,id_aero),1)
     END IF

     ! Pointer no more useful, deallocate. 
     DEALLOCATE(pt_tmp)

     ! Test if vertical interpolation will be needed.
     IF (psurf_year(1,1,id_aero)==not_valid .OR. pi_psurf_year(1,1,id_aero)==not_valid ) THEN
        ! Pressure=not_valid indicates old file format, see module readaerosol
        vert_interp = .FALSE.

        ! If old file format, both psurf_year and pi_psurf_year must be not_valid
        IF (  psurf_year(1,1,id_aero) /= pi_psurf_year(1,1,id_aero) ) THEN
           WRITE(lunout,*) 'Warning! All forcing files for the same aerosol must have the same structure'
           CALL abort_gcm('readaerosol_interp', 'The aerosol files have not the same format',1)
        END IF
        
        IF (klev /= klev_src) THEN
           WRITE(lunout,*) 'Old format of aerosol file do not allowed vertical interpolation'
           CALL abort_gcm('readaerosol_interp', 'Old aerosol file not possible',1)
        END IF

     ELSE 
        vert_interp = .TRUE.
     END IF

!    Calendar initialisation
!
     DO i = 2, 13
       month_len(i) = REAL(ioget_mon_len(year_cur, i-1))
       CALL ymds2ju(year_cur, i-1, 1, 0.0, month_start(i))
     ENDDO
     month_len(1) = REAL(ioget_mon_len(year_cur-1, 12))
     CALL ymds2ju(year_cur-1, 12, 1, 0.0, month_start(1))
     month_len(14) = REAL(ioget_mon_len(year_cur+1, 1))
     CALL ymds2ju(year_cur+1, 1, 1, 0.0, month_start(14))
     month_mid(:) = month_start (:) + month_len(:)/2.
     
     if (debug) then
       write(lunout,*)' month_len = ',month_len
       write(lunout,*)' month_mid = ',month_mid
     endif

  END IF  ! IF ( (first .OR. iday==0) .AND. lnewday ) THEN 
  
!****************************************************************************************
! - 2) Interpolate to the actual day.
! - 3) Interpolate to the model vertical grid.
!
!****************************************************************************************

  IF (lnewday) THEN ! only if new day
!****************************************************************************************
! 2) Interpolate to the actual day
! 
!****************************************************************************************
    ! Find which months and days to use for time interpolation
     nbr_tsteps = 12
     IF (nbr_tsteps == 12) then
       IF (jDay < month_mid(im+1)) THEN
          im2=im-1
          day2 = month_mid(im2+1)
          day1 = month_mid(im+1)
          IF (im2 <= 0) THEN
             ! the month is january, thus the month before december
             im2=12
          END IF
       ELSE
          ! the second half of the month
          im2=im+1
          day2 = month_mid(im+1)
          day1 = month_mid(im2+1)
          IF (im2 > 12) THEN
             ! the month is december, the following thus january
             im2=1
          ENDIF
       END IF
     ELSE IF (nbr_tsteps == 14) then
       im = im + 1
       IF (jDay < month_mid(im)) THEN
          ! in the first half of the month use month before and actual month
          im2=im-1
          day2 = month_mid(im2)
          day1 = month_mid(im)
       ELSE
          ! the second half of the month
          im2=im+1
          day2 = month_mid(im)
          day1 = month_mid(im2)
       END IF
     ELSE
       CALL abort_gcm('readaerosol_interp', 'number of months undefined',1)
     ENDIF
     if (debug) then
       write(lunout,*)' jDay, day1, day2, im, im2 = ', jDay, day1, day2, im, im2
     endif

 
     ! Time interpolation, still on vertical source grid
     ALLOCATE(tmp1(klon,klev_src), tmp2(klon,klev_src),stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 7',1)

     ALLOCATE(pplay_src(klon,klev_src), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('readaerosol_interp', 'pb in allocation 8',1)
     

     DO k=1,klev_src
        DO i=1,klon 
           tmp1(i,k) = &
                var_year(i,k,im2,id_aero) - (jDay-day2)/(day1-day2) * &
                (var_year(i,k,im2,id_aero) - var_year(i,k,im,id_aero))
           
           tmp2(i,k) = &
                pi_var_year(i,k,im2,id_aero) - (jDay-day2)/(day1-day2) * &
                (pi_var_year(i,k,im2,id_aero) - pi_var_year(i,k,im,id_aero))
        END DO
     END DO

     ! Time interpolation for pressure at surface, still on vertical source grid
     DO i=1,klon 
        psurf_day(i) = &
             psurf_year(i,im2,id_aero) - (jDay-day2)/(day1-day2) * &
             (psurf_year(i,im2,id_aero) - psurf_year(i,im,id_aero))
        
        pi_psurf_day(i) = &
             pi_psurf_year(i,im2,id_aero) - (jDay-day2)/(day1-day2) * &
             (pi_psurf_year(i,im2,id_aero) - pi_psurf_year(i,im,id_aero))
     END DO

     ! Time interpolation for the load, still on vertical source grid
     DO i=1,klon 
        load_src(i) = &
             load_year(i,im2,id_aero) - (jDay-day2)/(day1-day2) * &
             (load_year(i,im2,id_aero) - load_year(i,im,id_aero))
        
        pi_load_src(i) = &
             pi_load_year(i,im2,id_aero) - (jDay-day2)/(day1-day2) * &
             (pi_load_year(i,im2,id_aero) - pi_load_year(i,im,id_aero))
     END DO

!****************************************************************************************
! 3) Interpolate to the model vertical grid (target grid)
!
!****************************************************************************************

     IF (vert_interp) THEN

        ! - Interpolate variable tmp1 (on source grid) to var_day (on target grid)
        !********************************************************************************
        ! a) calculate pression at vertical levels for the source grid using the
        !    hybrid-sigma coordinates ap and b and the surface pressure, variables from file.
        DO k = 1, klev_src
           DO i = 1, klon
              pplay_src(i,k)= pt_ap(k) + pt_b(k)*psurf_day(i)
           END DO
        END DO
        
        IF (debug) THEN
           CALL writefield_phy('psurf_day_src',psurf_day(:),1)
           CALL writefield_phy('pplay_src',pplay_src(:,:),klev_src)
           CALL writefield_phy('pplay',pplay(:,:),klev)
           CALL writefield_phy('day_src',tmp1,klev_src)
           CALL writefield_phy('pi_day_src',tmp2,klev_src)
        END IF

        ! b) vertical interpolation on pressure leveles
        CALL pres2lev(tmp1(:,:), var_day(:,:,id_aero), klev_src, klev, pplay_src, pplay, &
             1, klon, .FALSE.)
        
        IF (debug) CALL writefield_phy('day_tgt',var_day(:,:,id_aero),klev)
        
        ! c) adjust to conserve total aerosol mass load in the vertical pillar
        !    Calculate the load in the actual pillar and compare with the load
        !    read from aerosol file.
        
        ! Find the pressure difference in each model layer
        DO k = 1, klev
           DO i = 1, klon
              delp(i,k) = paprs(i,k) - paprs (i,k+1)
           END DO
        END DO

        ! Find the mass load in the actual pillar, on target grid
        load_tgt(:) = 0.
        DO k= 1, klev
           DO i = 1, klon
              zrho = pplay(i,k)/t_seri(i,k)/RD       ! [kg/m3]
              volm = var_day(i,k,id_aero)*1.E-9/zrho ! [kg/kg]
              load_tgt(i) = load_tgt(i) + 1/RG * volm *delp(i,k)
           END DO
        END DO
        
        ! Adjust, uniform
        DO k = 1, klev
           DO i = 1, klon
              var_day(i,k,id_aero) = var_day(i,k,id_aero)*load_src(i)/load_tgt(i) 
           END DO
        END DO
        
        IF (debug) THEN
           load_tgt_test(:) = 0.
           DO k= 1, klev
              DO i = 1, klon
                 zrho = pplay(i,k)/t_seri(i,k)/RD       ! [kg/m3]
                 volm = var_day(i,k,id_aero)*1.E-9/zrho ! [kg/kg]
                 load_tgt_test(i) = load_tgt_test(i) + 1/RG * volm*delp(i,k)
              END DO
           END DO
           
           CALL writefield_phy('day_tgt2',var_day(:,:,id_aero),klev)
           CALL writefield_phy('load_tgt',load_tgt(:),1)
           CALL writefield_phy('load_tgt_test',load_tgt_test(:),1)
           CALL writefield_phy('load_src',load_src(:),1)
        END IF

        ! - Interpolate variable tmp2 (source grid) to pi_var_day (target grid)
        !********************************************************************************
        ! a) calculate pression at vertical levels at source grid    
        DO k = 1, klev_src
           DO i = 1, klon
              pplay_src(i,k)= pt_ap(k) + pt_b(k)*pi_psurf_day(i)
           END DO
        END DO

        IF (debug) THEN
           CALL writefield_phy('pi_psurf_day_src',pi_psurf_day(:),1)
           CALL writefield_phy('pi_pplay_src',pplay_src(:,:),klev_src)
        END IF

        ! b) vertical interpolation on pressure leveles
        CALL pres2lev(tmp2(:,:), pi_var_day(:,:,id_aero), klev_src, klev, pplay_src, pplay, &
             1, klon, .FALSE.)

        IF (debug) CALL writefield_phy('pi_day_tgt',pi_var_day(:,:,id_aero),klev)

        ! c) adjust to conserve total aerosol mass load in the vertical pillar
        !    Calculate the load in the actual pillar and compare with the load
        !    read from aerosol file.

        ! Find the load in the actual pillar, on target grid
        load_tgt(:) = 0.
        DO k = 1, klev
           DO i = 1, klon
              zrho = pplay(i,k)/t_seri(i,k)/RD          ! [kg/m3]
              volm = pi_var_day(i,k,id_aero)*1.E-9/zrho ! [kg/kg]
              load_tgt(i) = load_tgt(i) + 1/RG * volm * delp(i,k)
           END DO
        END DO

        DO k = 1, klev
           DO i = 1, klon
              pi_var_day(i,k,id_aero) = pi_var_day(i,k,id_aero)*pi_load_src(i)/load_tgt(i)
           END DO
        END DO

        IF (debug) THEN
           load_tgt_test(:) = 0.
           DO k = 1, klev
              DO i = 1, klon
                 zrho = pplay(i,k)/t_seri(i,k)/RD          ! [kg/m3]
                 volm = pi_var_day(i,k,id_aero)*1.E-9/zrho ! [kg/kg]
                 load_tgt_test(i) = load_tgt_test(i) + 1/RG * volm * delp(i,k)
              END DO
           END DO
           CALL writefield_phy('pi_day_tgt2',pi_var_day(:,:,id_aero),klev)
           CALL writefield_phy('pi_load_tgt',load_tgt(:),1)
           CALL writefield_phy('pi_load_tgt_test',load_tgt_test(:),1)
           CALL writefield_phy('pi_load_src',pi_load_src(:),1)
        END IF


     ELSE   ! No vertical interpolation done

        var_day(:,:,id_aero)    = tmp1(:,:)
        pi_var_day(:,:,id_aero) = tmp2(:,:)

     END IF ! vert_interp


     ! Deallocation
     DEALLOCATE(tmp1, tmp2, pplay_src, stat=ierr)

!****************************************************************************************
! 4) Test for negative mass values
!
!****************************************************************************************
     IF (MINVAL(var_day(:,:,id_aero)) < 0.) THEN
        DO k=1,klev
           DO i=1,klon 
              ! Test for var_day
              IF (var_day(i,k,id_aero) < 0.) THEN
                 IF (jDay-day2 < 0.) WRITE(lunout,*) 'jDay-day2=',jDay-day2
                 IF (var_year(i,k,im2,id_aero) - var_year(i,k,im,id_aero) < 0.) THEN
                    WRITE(lunout,*) trim(name_aero(id_aero)),'(i,k,im2)-', &
                         trim(name_aero(id_aero)),'(i,k,im)=',           &
                         var_year(i,k,im2,id_aero) - var_year(i,k,im,id_aero)
                 END IF
                 WRITE(lunout,*) 'stop for aerosol : ',name_aero(id_aero)
                 WRITE(lunout,*) 'day1, day2, jDay = ', day1, day2, jDay 
                 CALL abort_gcm('readaerosol_interp','Error in interpolation 1',1)
              END IF
           END DO 
        END DO
     END IF

     IF (MINVAL(pi_var_day(:,:,id_aero)) < 0. ) THEN
        DO k=1, klev
           DO i=1,klon
              ! Test for pi_var_day
              IF (pi_var_day(i,k,id_aero) < 0.) THEN
                 IF (jDay-day2 < 0.) WRITE(lunout,*) 'jDay-day2=',jDay-day2
                 IF (pi_var_year(i,k,im2,id_aero) - pi_var_year(i,k,im,id_aero) < 0.) THEN
                    WRITE(lunout,*) trim(name_aero(id_aero)),'(i,k,im2)-', &
                         trim(name_aero(id_aero)),'(i,k,im)=',           &
                         pi_var_year(i,k,im2,id_aero) - pi_var_year(i,k,im,id_aero)
                 END IF
                 
                 WRITE(lunout,*) 'stop for aerosol : ',name_aero(id_aero)
                 CALL abort_gcm('readaerosol_interp','Error in interpolation 2',1)
              END IF
           END DO
        END DO
     END IF

  END IF ! lnewday

!****************************************************************************************
! Copy output from saved variables
!
!****************************************************************************************

  mass_out(:,:)    = var_day(:,:,id_aero) 
  pi_mass_out(:,:) = pi_var_day(:,:,id_aero)
  
END SUBROUTINE readaerosol_interp
