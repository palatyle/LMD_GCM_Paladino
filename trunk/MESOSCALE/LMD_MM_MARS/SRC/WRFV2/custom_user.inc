!---------------------------------------------
!
! custom_user.inc
!
! TOUTE MODIFICATION IMPLIQUE DE RECOMPILER
! EN N'OUBLIANT PAS DE PRECEDER LA
! COMPILATION DE :
! "touch dyn_em/module_initialize_real.F"
!
!---------------------------------------------
! cas où config_flags%adjust_heights est true
!---------------------------------------------
print *,  '**********************************'
print *,  'PLEASE NOTE YOU ARE IN CUSTOM MODE'
print *,  '**********************************'


!!
!! stuff for LES
!!
!
!	 print *, 'setting uniform values and profiles'
!        print *, 'u'
!        print *, grid%em_u_gc(its+1,:,jts+1)
!        print *, 'v'
!        print *, grid%em_v_gc(its+1,:,jts+1)
!        print *, 't'
!        print *, grid%em_t_gc(its+1,:,jts+1)
!        print *, 'p'
!        print *, grid%em_rh_gc(its+1,:,jts+1)
!        print *, 'geop'
!        print *, grid%em_ght_gc(its+1,:,jts+1)
!        print *, 'albedo'
!        print *, grid%u10(its+1,jts+1)
!        print *, 'thermal inertia'
!        print *, grid%v10(its+1,jts+1)
!        print *, 'topography'
!        print *, grid%ht(its+1,jts+1)
!        print *, 'toposoil'
!        print *, grid%toposoil(its+1,jts+1)
!        print *, 'surface temperature'
!        print *, grid%tsk(its+1,jts+1)
!        print *, 'surface pressure'
!        print *, grid%psfc(its+1,jts+1)
!        print *, grid%em_psfc_gc(its+1,jts+1)
!
!        DO j = jts, MIN(jte,jde-1)
!         DO i = its, MIN(ite,ide-1)
!                grid%em_u_gc(i,:,j)=grid%em_u_gc(its+1,:,jts+1)
!                grid%em_v_gc(i,:,j)=grid%em_v_gc(its+1,:,jts+1)
!                grid%em_t_gc(i,:,j)=grid%em_t_gc(its+1,:,jts+1)
!                grid%em_rh_gc(i,:,j)=grid%em_rh_gc(its+1,:,jts+1)
!                grid%em_ght_gc(i,:,j) = grid%em_ght_gc(its+1,:,jts+1)
!                grid%u10(i,j) = grid%u10(its+1,jts+1)
!                grid%v10(i,j) = grid%v10(its+1,jts+1)
!                grid%ht(i,j) = grid%ht(its+1,jts+1)
!                grid%toposoil(i,j) = grid%toposoil(its+1,jts+1)
!                grid%tsk(i,j) = grid%tsk(its+1,jts+1)
!                grid%psfc(i,j) = grid%psfc(its+1,jts+1)
!                grid%em_psfc_gc(i,j) = grid%em_psfc_gc(its+1,jts+1)
!          ENDDO
!        ENDDO
!
!!!
!!!
!print *, 'u divided by 10'
!grid%em_u_gc=grid%em_u_gc/10.
!grid%em_v_gc=grid%em_v_gc/10.
!!!
!!!


!!****************************
!! CONSTANT THERMAL INERTIA
!!****************************
print *, 'fixed thermal inertia'

DO j = jts, MIN(jte,jde-1)
DO i = its, MIN(ite,ide-1)

       grid%v10(i,j) = 85.    ! thermal inertia

ENDDO
ENDDO




!!!******************************
!!! IDEALIZED CRATER SIMULATION
!!!******************************
!print *, 'fixed albedo, thermal inertia, zonal and meridional wind'
!
!DO j = jts, MIN(jte,jde-1)
!DO i = its, MIN(ite,ide-1)
!
!	grid%u10(i,j) = 0.19	! albedo
!	grid%v10(i,j) = 315.	! thermal inertia
!
!	!! INPUT: LIKE A "SOUNDING"
!	grid%em_u_gc(i,:,j)=grid%em_u_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))	! zonal wind
!	grid%em_v_gc(i,:,j)=0.
!!!grid%em_v_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))	! meridional wind
!	!! INPUT: LIKE A "SOUNDING"
!
!ENDDO
!ENDDO
!
!!!!!INVERSE WIND
!!!grid%em_u_gc = -grid%em_u_gc
!
!
!!! FIX for the STAGGERED SPECIFICITY
!grid%em_u_gc(MIN(ite,ide-1)+1,:,:)=grid%em_u_gc(MIN(ite,ide-1),:,:)
!grid%em_v_gc(:,:,MIN(jte,jde-1)+1)=grid%em_v_gc(:,:,MIN(jte,jde-1))
!
!
!!! CHECK
!print *, 'wind profile'
!print *, 'took at ...', int((MIN(ite,ide-1)-its)/2), int((MIN(jte,jde-1)-jts)/2)
!print *, '--zonal'
!print *, grid%em_u_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))
!print *, '--meridional'
!print *, grid%em_v_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))
!

!!
!! Later, sigma levels are changed to avoid extrapolation
!! ... so no change in the pressure arrays is needed
!!



!!grid%u10 = grid%u10*0. + 0.19	! albedo
!!grid%v10 = grid%v10*0. + 315.	! thermal inertia
!!grid%em_u_gc = grid%em_u_gc*0. - 15.	! zonal wind
!!grid%em_v_gc = grid%em_v_gc*0. 		! meridional wind
!!
!!grid%em_u_gc(:,1,:)=0.
!!grid%em_v_gc(:,1,:)=0.


!! TEST TRACERS
!grid%scalar(its+10,1:10,jts+10,:)=10.


!!!!!!!!!!!!!!!!!!!
!!! READ PROFILE !!
!!!!!!!!!!!!!!!!!!!
!
!open(unit=10,file='input_sounding',form='formatted',status='old')
!rewind(10)
!read(10,*) grid%em_u_gc(1,:,1)
!
!DO j = jts, MIN(jte,jde-1)
!DO i = its, MIN(ite,ide-1)
!
!!!        grid%u10(i,j) = 0.19    ! albedo
!!!        grid%v10(i,j) = 315.    ! thermal inertia
!
!        grid%u10(i,j) = 0.29    ! albedo
!        grid%v10(i,j) = 100.    ! thermal inertia
!
!        grid%em_u_gc(i,:,j)=grid%em_u_gc(1,:,1)     ! zonal wind
!        grid%em_v_gc(i,:,j)=0.			    ! meridional wind
!
!ENDDO
!ENDDO
!
!
!!! FIX for the STAGGERED SPECIFICITY
!grid%em_u_gc(MIN(ite,ide-1)+1,:,:)=grid%em_u_gc(MIN(ite,ide-1),:,:)
!grid%em_v_gc(:,:,MIN(jte,jde-1)+1)=grid%em_v_gc(:,:,MIN(jte,jde-1))
!
!
!print *, 'wind profile'
!print *, '--zonal'
!print *, grid%em_u_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))
!print *, '--meridional'
!print *, grid%em_v_gc(int((MIN(ite,ide-1)-its)/2),:,int((MIN(jte,jde-1)-jts)/2))




