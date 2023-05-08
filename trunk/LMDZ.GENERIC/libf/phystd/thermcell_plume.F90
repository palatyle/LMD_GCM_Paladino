!
!
!
SUBROUTINE thermcell_plume(ngrid,nlay,nq,ptimestep,                           &
                           ztv,zhl,zqt,zql,zlev,pplev,zpopsk,                 &
                           detr_star,entr_star,f_star,                        &
                           ztva,zhla,zqta,zqla,zqsa,                          &
                           zw2,lbot,lmin)
      
      
!===============================================================================
!  Purpose: calcule les valeurs de qt, thetal et w dans l ascendance
!  
!  Nota Bene
!     ql   means "non-gaseous water mass mixing ratio" (liquid and solid)
!     qv   means "vapor mass mixing ratio"
!     qt   means "total water mass mixing ratio"
!     TP   means "potential temperature"
!     TRPV means "virtual potential temperature with latent heat release"  
!     TPV  means "virtual potential temperature"
!     TR   means "temperature with latent heat release"
!===============================================================================
      
      USE print_control_mod, ONLY: prt_level
      USE watercommon_h, ONLY: RLvCp, RETV, Psat_water
      USE tracer_h, ONLY: igcm_h2o_vap
      USE thermcell_mod
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: nq
      
      INTEGER, INTENT(in) :: lbot(ngrid)              ! First considered layer
      
      REAL, INTENT(in) :: ptimestep
      REAL, INTENT(in) :: zlev(ngrid,nlay+1)          ! Levels altitude
      REAL, INTENT(in) :: pplev(ngrid,nlay+1)         ! Levels pressure
      REAL, INTENT(in) :: zpopsk(ngrid,nlay)          ! Exner function
      
      REAL, INTENT(in) :: ztv(ngrid,nlay)             ! TRPV environment
      REAL, INTENT(in) :: zhl(ngrid,nlay)             ! TP   environment
      REAL, INTENT(in) :: zqt(ngrid,nlay)             ! qt   environment
      REAL, INTENT(in) :: zql(ngrid,nlay)             ! ql   environment
      
!     Outputs:
!     --------
      
      INTEGER, INTENT(out) :: lmin(ngrid)             ! Plume bottom level (first unstable level)
      
      REAL, INTENT(out) :: detr_star(ngrid,nlay)      ! Normalized detrainment
      REAL, INTENT(out) :: entr_star(ngrid,nlay)      ! Normalized entrainment
      REAL, INTENT(out) :: f_star(ngrid,nlay+1)       ! Normalized mass flux
      
      REAL, INTENT(out) :: ztva(ngrid,nlay)           ! TRPV plume (after mixing)
      REAL, INTENT(out) :: zhla(ngrid,nlay)           ! TP   plume (after mixing)
      REAL, INTENT(out) :: zqla(ngrid,nlay)           ! ql   plume (after mixing)
      REAL, INTENT(out) :: zqta(ngrid,nlay)           ! qt   plume (after mixing)
      REAL, INTENT(out) :: zqsa(ngrid,nlay)           ! qsat plume (after mixing)
      REAL, INTENT(out) :: zw2(ngrid,nlay+1)          ! w    plume (after mixing)
      
!     Local:
!     ------
      
      INTEGER ig, l, k
      INTEGER l_start
      
      REAL ztva_est(ngrid,nlay)                       ! TRPV plume (before mixing)
      REAL zqla_est(ngrid,nlay)                       ! ql   plume (before mixing)
      REAL zta_est(ngrid,nlay)                        ! TR   plume (before mixing)
      REAL zqsa_est(ngrid)                            ! qsat plume (before mixing)
      REAL zw2_est(ngrid,nlay+1)                      ! w    plume (before mixing)
      
      REAL zta(ngrid,nlay)                            ! TR   plume (after mixing)
      
      REAL zbuoy(ngrid,nlay)                          ! Plume buoyancy
      REAL ztemp(ngrid)                               ! Temperature to compute saturation vapor pressure
      REAL zdz                                        ! Layers heights
      REAL ztv2(ngrid,nlay)                           ! ztv + d_temp * Dirac(l=linf)
      
      REAL zdw2                                       ! 
      REAL zw2fact                                    ! 
      REAL zw2m                                       ! Average vertical velocity between two successive levels
      REAL gamma                                      ! Plume acceleration term (to compute vertical velocity)
      REAL test                                       ! 
      
      REAL psat                                       ! Dummy argument for Psat_water()
      
      LOGICAL active(ngrid)                           ! If the plume is active (speed and incoming mass flux > 0)
      LOGICAL activetmp(ngrid)                        ! If the plume is active (active=true and outgoing mass flux > 0)
      
!===============================================================================
! Initialization
!===============================================================================
      
      ztva(:,:) = ztv(:,:)                                                       ! ztva is set to TPV environment
      zhla(:,:) = zhl(:,:)                                                       ! zhla is set to TP  environment
      zqta(:,:) = zqt(:,:)                                                       ! zqta is set to qt  environment
      zqla(:,:) = zql(:,:)                                                       ! zqla is set to ql  environment
      zqsa(:,:) = 0.
      zw2(:,:)  = 0.
      
      ztva_est(:,:) = ztv(:,:)                                                   ! ztva_est is set to TPV environment
      zqla_est(:,:) = zql(:,:)                                                   ! zqla_est is set to ql  environment
      zqsa_est(:)   = 0.
      zw2_est(:,:)  = 0.
      
      zbuoy(:,:) = 0.
      
      f_star(:,:)    = 0.
      detr_star(:,:) = 0.
      entr_star(:,:) = 0.
      
      lmin(:) = lbot(:)
      
      ztv2(:,:)    = ztv(:,:)
      ztv2(:,linf) = ztv(:,linf) + d_temp
      
      active(:) = .false.
      
      l_start = nlay
      
!===============================================================================
! First layer computation
!===============================================================================
      
      DO ig=1,ngrid
         l = lbot(ig)
         l_start = MIN(l_start, lbot(ig)+1)
         DO WHILE (.not.active(ig).and.(pplev(ig,l+1) > pres_limit).and.(l < nlay))
            zbuoy(ig,l) = RG * (ztv2(ig,l) - ztv2(ig,l+1)) / ztv2(ig,l+1)
            IF (zbuoy(ig,l) > 0.) THEN
               lmin(ig) = l
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: entrainement and mass flux initial values are set to 1. The physical value
!     will be computed thanks to the closure relation in thermcell_closure.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               entr_star(ig,l) = 1.
               f_star(ig,l+1) = 1.
               zdz = zlev(ig,l+1) - zlev(ig,l)
               zw2fact = 2. * fact_epsilon * zdz / (1. + betalpha)
               zdw2 = 2. * afact * zbuoy(ig,l) * zdz / (1. + betalpha)
               zw2_est(ig,l+1) = exp(-zw2fact) * zdw2
               zw2(ig,l+1) = zw2_est(ig,l+1)
               active(ig) = .true.
            ENDIF
            l = l + 1
         ENDDO
      ENDDO
      
!===============================================================================
! Thermal plumes computations
!===============================================================================
      
      DO l=l_start,nlay-1
         
!-------------------------------------------------------------------------------
! Is thermal plume (still) active ?
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            active(ig) = (zw2(ig,l) > 1.e-9).and.(f_star(ig,l) > 1.e-9)
         ENDDO
         
!-------------------------------------------------------------------------------
! Latent heat release (before mixing)
!-------------------------------------------------------------------------------
         
         ztemp(:) = zpopsk(:,l) * zhla(:,l-1)
         
         DO ig=1,ngrid
            IF (active(ig)) THEN
               CALL Psat_water(ztemp(ig), pplev(ig,l), psat, zqsa_est(ig))
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Vertical speed (before mixing)
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (active(ig)) THEN
               zqla_est(ig,l) = MAX(0.,zqta(ig,l-1) - zqsa_est(ig))              ! zqla_est is set to ql   plume
               zta_est(ig,l)  = zhla(ig,l-1) * zpopsk(ig,l)                   &  ! zta_est  is set to TR   plume
               &              + RLvCp * zqla_est(ig,l)
               ztva_est(ig,l) = zta_est(ig,l) / zpopsk(ig,l)                  &  ! ztva_est is set to TRPV plume
               &              * (1. + RETV * (zqta(ig,l-1)-zqla_est(ig,l)) - zqla_est(ig,l))
               
               zbuoy(ig,l) = RG * (ztva_est(ig,l) - ztv(ig,l)) / ztv(ig,l)
               zdz = zlev(ig,l+1) - zlev(ig,l)
               
               zw2fact = 2. * fact_epsilon * zdz / (1. + betalpha)
               zdw2 = afact * zbuoy(ig,l) / fact_epsilon
               zw2_est(ig,l+1) = MAX(0., exp(-zw2fact) * (zw2_est(ig,l) - zdw2) + zdw2)
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Mass flux, entrainment and detrainment
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (active(ig)) THEN
               
               zdz = zlev(ig,l+1) - zlev(ig,l)
               zw2m = (zw2_est(ig,l+1) + zw2(ig,l)) / 2.
               gamma = afact * zbuoy(ig,l) - fact_epsilon * zw2m
               
               IF (zw2m > 0.) THEN
                  test = gamma / zw2m - nu
               ELSE
                  test = 0.
                  print *, 'WARNING: vertical speed is negative while plume is active!'
                  print *, 'ig,l', ig, l
                  print *, 'zw2m', zw2m
               ENDIF
               
               IF (test > 0.) THEN
                  detr_star(ig,l) = zdz * f_star(ig,l) * nu
                  entr_star(ig,l) = zdz * f_star(ig,l) * (betalpha * gamma / zw2m + nu) / (betalpha + 1)
               ELSE
                  detr_star(ig,l) = zdz * f_star(ig,l) * ((betalpha + 1) * nu - betalpha * gamma / zw2m)
                  entr_star(ig,l) = zdz * f_star(ig,l) * nu
               ENDIF
               
               f_star(ig,l+1) = f_star(ig,l) + entr_star(ig,l) - detr_star(ig,l)
               
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Mixing between thermal plume and environment
!-------------------------------------------------------------------------------
         
         activetmp(:) = active(:).and.(f_star(:,l+1) > 1.e-9)
         
         DO ig=1,ngrid
            IF (activetmp(ig)) THEN
               zhla(ig,l) = (f_star(ig,l) * zhla(ig,l-1)                      &  ! zhla is set to TP in plume (mixed)
               &          + entr_star(ig,l) * zhl(ig,l))                      &
               &          / (f_star(ig,l+1) + detr_star(ig,l))
               zqta(ig,l) = (f_star(ig,l) * zqta(ig,l-1) +                    &  ! zqta is set to qt in plume (mixed)
               &          + entr_star(ig,l) * zqt(ig,l))                      &
               &          / (f_star(ig,l+1) + detr_star(ig,l))
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Latent heat release (after mixing)
!-------------------------------------------------------------------------------
         
         ztemp(:) = zpopsk(:,l) * zhla(:,l)
         
         DO ig=1,ngrid
            IF (activetmp(ig)) THEN
               CALL Psat_water(ztemp(ig), pplev(ig,l), psat, zqsa(ig,l))
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Vertical speed (after mixing)
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (activetmp(ig)) THEN
               zqla(ig,l) = MAX(0.,zqta(ig,l) - zqsa(ig,l))                      ! zqla is set to ql   plume (mixed)
               zta(ig,l)  = zhla(ig,l) * zpopsk(ig,l)                         &  ! ztva is set to TR   plume (mixed)
               &          + RLvCp * zqla(ig,l)
               ztva(ig,l) = zta(ig,l) / zpopsk(ig,l)                          &  ! ztva is set to TRPV plume (mixed)
               &          * (1. + RETV*(zqta(ig,l)-zqla(ig,l)) - zqla(ig,l))
               
               zbuoy(ig,l) = RG * (ztva(ig,l) - ztv(ig,l)) / ztv(ig,l)
               zdz = zlev(ig,l+1) - zlev(ig,l)
               
               zw2fact = 2. * fact_epsilon * zdz / (1. + betalpha)
               zdw2 = afact * zbuoy(ig,l) / fact_epsilon
               zw2(ig,l+1) = MAX(0., exp(-zw2fact) * (zw2(ig,l) - zdw2) + zdw2)
            ENDIF
         ENDDO
         
      ENDDO
      
      
RETURN
END
