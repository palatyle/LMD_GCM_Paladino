










!
!
!
SUBROUTINE thermcell_env(ngrid,nlay,nq,pq,pt,pu,pv,pplay,pplev,               &
                         zqt,zql,zt,ztv,zhl,zu,zv,zpopsk,zqs)
      
      
!===============================================================================
!  Purpose: calcul des caracteristiques de l'environnement necessaires au calcul
!           des proprietes dans le thermique.
!  
!  Modif 2019/04 (AB alexandre.boissinot@lmd.jussieu.fr)
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
      USE thermcell_mod, ONLY: RKAPPA
      USE watercommon_h, ONLY: RLvCp, RETV, Psat_water
      USE tracer_h, ONLY: igcm_h2o_vap, igcm_h2o_ice
      USE callkeys_mod, ONLY: water
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: nq
      
      REAL, INTENT(in) :: pq(ngrid,nlay,nq)           ! Large scale water
      REAL, INTENT(in) :: pt(ngrid,nlay)              ! Large scale temperature
      REAL, INTENT(in) :: pu(ngrid,nlay)              ! Large scale zonal wind
      REAL, INTENT(in) :: pv(ngrid,nlay)              ! Large scale meridional wind
      REAL, INTENT(in) :: pplay(ngrid,nlay)           ! Layers pressure
      REAL, INTENT(in) :: pplev(ngrid,nlay+1)         ! Levels pressure
      REAL, INTENT(in) :: zpopsk(ngrid,nlay)          ! Exner function
      
!     Outputs:
!     --------
      
      REAL, INTENT(out) :: zt(ngrid,nlay)             ! T    environment
      REAL, INTENT(out) :: ztv(ngrid,nlay)            ! TRPV environment
      REAL, INTENT(out) :: zhl(ngrid,nlay)            ! TP   environment
      REAL, INTENT(out) :: zu(ngrid,nlay)             ! u    environment
      REAL, INTENT(out) :: zv(ngrid,nlay)             ! v    environment
      REAL, INTENT(out) :: zqt(ngrid,nlay)            ! qt   environment
      REAL, INTENT(out) :: zql(ngrid,nlay)            ! ql   environment
      REAL, INTENT(out) :: zqs(ngrid,nlay)            ! qsat environment
      
!     Local:
!     ------
      
      INTEGER ig, k
      
      REAL psat                                       ! Dummy argument for Psat_water()
      
!===============================================================================
! Initialization
!===============================================================================
      
      zu(:,:) = pu(:,:)
      zv(:,:) = pv(:,:)
      
      zhl(:,:) = pt(:,:) / zpopsk(:,:)
      
      zqt(:,:) = 0.
      zql(:,:) = 0.
      
!===============================================================================
! Condensation and latent heat release
!===============================================================================
      
      IF (water) THEN
         
         zqt(:,:) = pq(:,:,igcm_h2o_vap)
         
         DO k=1,nlay
            DO ig=1,ngrid
               CALL Psat_water(pt(ig,k), pplev(ig,k), psat, zqs(ig,k))
            ENDDO
         ENDDO
         
         DO k=1,nlay
            DO ig=1,ngrid
               zql(ig,k) = max(0.,pq(ig,k,igcm_h2o_vap) - zqs(ig,k))
               zt(ig,k) = pt(ig,k) + RLvCp * zql(ig,k)
               ztv(ig,k) = zt(ig,k) / zpopsk(ig,k)                            &
               &         * (1. + RETV * (zqt(ig,k)-zql(ig,k)) - zql(ig,k))
            ENDDO
         ENDDO
         
      ELSE
         
         zt(:,:) = pt(:,:)
         ztv(:,:) = pt(:,:) / zpopsk(:,:)
         
      ENDIF
      
      
RETURN
END
