










!
!
!
SUBROUTINE thermcell_main(ngrid,nlay,nq,ptimestep,firstcall,                  &
                          pplay,pplev,pphi,zpopsk,                            &
                          pu,pv,pt,pq,                                        &
                          pduadj,pdvadj,pdtadj,pdqadj,                        &
                          fm_tot,entr_tot,detr_tot,zw2_tot,fraca)
      
      
!===============================================================================
!   Auteurs: Frederic Hourdin, Catherine Rio, Anne Mathieu
!   Version du 09.02.07
!   Calcul du transport vertical dans la couche limite en presence
!   de "thermiques" explicitement representes avec processus nuageux
!
!   Reecriture a partir d'un listing papier a Habas, le 14/02/00
!
!   le thermique est suppose homogene et dissipe par melange avec
!   son environnement. la longueur l_mix controle l'efficacite du
!   melange
!
!   Le calcul du transport des differentes especes se fait en prenant
!   en compte:
!     1. un flux de masse montant
!     2. un flux de masse descendant
!     3. un entrainement
!     4. un detrainement
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!    Introduction of an implicit computation of vertical advection in
!    the environment of thermal plumes in thermcell_dq
!    impl = 0 : explicit ; impl = 1 : implicit ; impl =-1 : old version
!    controled by iflag_thermals =
!       15, 16 run with impl=-1 : numerical convergence with NPv3
!       17, 18 run with impl=1  : more stable
!    15 and 17 correspond to the activation of the stratocumulus "bidouille"
! 
! Major changes 2018-19 (AB alexandre.boissinot@lmd.jussieu.fr)
!    New detr and entre formulae (no longer alimentation)
!    lmin can be greater than 1
!    Mix every tracer
!    Can stack verticaly multiple plumes (it makes thermcell_dv2 unusable for the moment)
! 
!===============================================================================
      
      USE thermcell_mod
      USE print_control_mod, ONLY: prt_level
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: nq
      
      REAL, INTENT(in) :: ptimestep
      REAL, INTENT(in) :: pplay(ngrid,nlay)           ! Layer pressure
      REAL, INTENT(in) :: pplev(ngrid,nlay+1)         ! Level pressure
      REAL, INTENT(in) :: pphi(ngrid,nlay)            ! Geopotential
      REAL, INTENT(in) :: zpopsk(ngrid,nlay)          ! Exner function
      
      REAL, INTENT(in) :: pu(ngrid,nlay)              ! Zonal wind
      REAL, INTENT(in) :: pv(ngrid,nlay)              ! Meridional wind
      REAL, INTENT(in) :: pt(ngrid,nlay)              ! Temperature
      REAL, INTENT(in) :: pq(ngrid,nlay,nq)           ! Tracers mass mixing ratio
      
      LOGICAL, INTENT(in) :: firstcall
      
!     Outputs:
!     --------
      
      REAL, INTENT(out) :: pduadj(ngrid,nlay)         ! u convective variations
      REAL, INTENT(out) :: pdvadj(ngrid,nlay)         ! v convective variations
      REAL, INTENT(out) :: pdtadj(ngrid,nlay)         ! t convective variations
      REAL, INTENT(out) :: pdqadj(ngrid,nlay,nq)      ! q convective variations
      
      REAL, INTENT(inout) :: fm_tot(ngrid,nlay+1)     ! Total mass flux
      REAL, INTENT(inout) :: entr_tot(ngrid,nlay)     ! Total entrainment
      REAL, INTENT(inout) :: detr_tot(ngrid,nlay)     ! Total detrainment
      
      REAL, INTENT(out) :: fraca(ngrid,nlay+1)        ! Updraft fraction
      REAL, INTENT(out) :: zw2_tot(ngrid,nlay+1)      ! Total plume vertical speed
      
!     Local:
!     ------
      
      INTEGER ig, k, l, iq
      INTEGER lmax(ngrid)                             ! Highest layer reached by the plume
      INTEGER lmix(ngrid)                             ! Layer in which plume vertical speed is maximal
      INTEGER lmin(ngrid)                             ! First unstable layer
      
      REAL zmix(ngrid)                                ! Altitude of maximal vertical speed
      REAL zmax(ngrid)                                ! Maximal altitudes where plumes are active
      REAL zmin(ngrid)                                ! Minimal altitudes where plumes are active
      
      REAL zlay(ngrid,nlay)                           ! Layers altitudes
      REAL zlev(ngrid,nlay+1)                         ! Levels altitudes
      REAL rho(ngrid,nlay)                            ! Layers densities
      REAL rhobarz(ngrid,nlay)                        ! Levels densities
      REAL masse(ngrid,nlay)                          ! Layers masses
      
      REAL zu(ngrid,nlay)                             ! u    environment
      REAL zv(ngrid,nlay)                             ! v    environment
      REAL zt(ngrid,nlay)                             ! TR   environment
      REAL zqt(ngrid,nlay)                            ! qt   environment
      REAL zql(ngrid,nlay)                            ! ql   environment
      REAL zhl(ngrid,nlay)                            ! TP   environment
      REAL ztv(ngrid,nlay)                            ! TRPV environment
      REAL zqs(ngrid,nlay)                            ! qsat environment
      
      REAL zua(ngrid,nlay)                            ! u    plume
      REAL zva(ngrid,nlay)                            ! v    plume
      REAL zqla(ngrid,nlay)                           ! qv   plume
      REAL zqta(ngrid,nlay)                           ! qt   plume
      REAL zhla(ngrid,nlay)                           ! TP   plume
      REAL ztva(ngrid,nlay)                           ! TRPV plume
      REAL zqsa(ngrid,nlay)                           ! qsat plume
      
      REAL zqa(ngrid,nlay,nq)                         ! q    plume (ql=0, qv=qt)
      
      REAL f_star(ngrid,nlay+1)                       ! Normalized mass flux
      REAL entr_star(ngrid,nlay)                      ! Normalized entrainment
      REAL detr_star(ngrid,nlay)                      ! Normalized detrainment
      
      REAL f(ngrid)                                   ! Mass flux norm
      REAL fm(ngrid,nlay+1)                           ! Mass flux
      REAL entr(ngrid,nlay)                           ! Entrainment
      REAL detr(ngrid,nlay)                           ! Detrainment
      
      REAL zw2(ngrid,nlay+1)                          ! Plume vertical speed
      REAL wmax(ngrid)                                ! Maximal vertical speed
      REAL zdthladj(ngrid,nlay)                       ! Potential temperature variations
      REAL dummy(ngrid,nlay)                          ! Dummy argument for thermcell_dq()
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER lbot(ngrid)
      LOGICAL re_tpm
      INTEGER while_loop_counter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!===============================================================================
! Initialization
!===============================================================================
      
      fm_tot(:,:) = 0.
      entr_tot(:,:) = 0.
      detr_tot(:,:) = 0.
      zw2_tot(:,:) = 0.
      
      pduadj(:,:) = 0.0
      pdvadj(:,:) = 0.0
      pdtadj(:,:) = 0.0
      pdqadj(:,:,:) = 0.0
      
      zdthladj(:,:) = 0.0
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      re_tpm = .true.
      lbot(:) = linf
      while_loop_counter = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!===============================================================================
! Environment settings
!===============================================================================
      
!-------------------------------------------------------------------------------
! Calcul de T,q,ql a partir de Tl et qt dans l environnement
!-------------------------------------------------------------------------------
      
      CALL thermcell_env(ngrid,nlay,nq,pq,pt,pu,pv,pplay,pplev,               &
      &                  zqt,zql,zt,ztv,zhl,zu,zv,zpopsk,zqs)
      
!-------------------------------------------------------------------------------
! Levels and layers altitudes
!-------------------------------------------------------------------------------
      
      DO l=2,nlay
         zlev(:,l) = 0.5 * (pphi(:,l) + pphi(:,l-1)) / RG
      ENDDO
      
      zlev(:,1) = 0.
      zlev(:,nlay+1) = (2. * pphi(:,nlay) - pphi(:,nlay-1)) / RG
      
      DO l=1,nlay
         zlay(:,l) = pphi(:,l) / RG
      ENDDO
      
!-------------------------------------------------------------------------------
! Levels and layers densities
!-------------------------------------------------------------------------------
      
      rho(:,:) = pplay(:,:) / (zpopsk(:,:) * RD * ztv(:,:))
      
      rhobarz(:,1) = rho(:,1)
      IF (prt_level.ge.10) THEN
         print *, 'WARNING: density in the first layer is equal to density at the first level!'
      ENDIF
      
      DO l=2,nlay
         rhobarz(:,l) = 0.5 * (rho(:,l) + rho(:,l-1))
      ENDDO
      
!-------------------------------------------------------------------------------
! Layers masses
!-------------------------------------------------------------------------------
      
      DO l=1,nlay
         masse(:,l) = (pplev(:,l) - pplev(:,l+1)) / RG
      ENDDO
      
!===============================================================================
! Explicative schemes
!===============================================================================

!-------------------------------------------------------------------------------
! Thermal plume variables
!-------------------------------------------------------------------------------
      
!           top of the model     
!     ===========================
!                                
!     ---------------------------
!                                       _
!     ----- F_lmax+1=0 ------zmax        !     lmax                                |
!     ------F_lmax>0-------------         |
!                                         |
!     ---------------------------         |
!                                         |
!     ---------------------------         |
!                                         |
!     ------------------wmax,zmix         |
!     lmix                                |
!     ---------------------------         |
!                                         |
!     ---------------------------         |
!                                         | E, D
!     ---------------------------         |
!                                         |
!     --------------------------- rhobarz, f_star, fm, zw2, fraca
!         zt, zu, zv, zo, rho             |
!     ---------------------------         |
!                                         |
!     ---------------------------         |
!                                         |
!     ---------------------------         |
!                                         |
!     ------F_lmin+1>0-----------         |
!     lmin                                |
!     ----- F_lmin=0 ------------       _/
!                                
!     ---------------------------
!                                
!     ===========================
!         bottom of the model    
      
!-------------------------------------------------------------------------------
! Zoom on layers k and k-1
!-------------------------------------------------------------------------------
      
!     |     /|\                  |                          |
!     |----  |  F_k+1 -----------|--------------------------|   level k+1
!     |      |  w_k+1         |                             |
!     |                     --|--> D_k                      |
!     |                       |                             |   layer k
!     |                    <--|--  E_k                      |
!     |     /|\               |                             |
!     |----  |  F_k ----------|-----------------------------|   level k
!     |      |  w_k        |                                |
!     |                  --|--> D_k-1                       |
!     |                    |                                |   layer k-1
!     |                 <--|--  E_k-1                       |
!     |     /|\            |                                |
!     |----  |  F_k-1 -----|--------------------------------|   level k-1
!            |  w_k-1                                        
!     0                  fraca                              1
!      \__________________/ \______________________________/
!          plume (fraca)          environment (1-fraca)   
      
!===============================================================================
! Thermal plumes computation
!===============================================================================
      
      DO WHILE (re_tpm.and.(while_loop_counter<nlay))
         while_loop_counter = while_loop_counter + 1
         
!-------------------------------------------------------------------------------
! Thermal plumes speeds, normalized fluxes, tracers and temperatures
!-------------------------------------------------------------------------------
         
         CALL thermcell_plume(ngrid,nlay,nq,ptimestep,                        &
         &                    ztv,zhl,zqt,zql,zlev,pplev,zpopsk,              &
         &                    detr_star,entr_star,f_star,                     &
         &                    ztva,zhla,zqta,zqla,zqsa,                       &
         &                    zw2,lbot,lmin)
         
!-------------------------------------------------------------------------------
! Thermal plumes characteristics: zmax, zmix, wmax
!-------------------------------------------------------------------------------
         
! AB: Careful, zw2 became its square root in thermcell_height!
         CALL thermcell_height(ngrid,nlay,lmin,lmix,lmax,                     &
         &                     zlev,zmin,zmix,zmax,zw2,wmax,f_star)
         
!-------------------------------------------------------------------------------
! Closure
!-------------------------------------------------------------------------------
         
         CALL thermcell_closure(ngrid,nlay,ptimestep,rho,zlev,                &
         &                      lmax,entr_star,zmin,zmax,wmax,f)
         
! FH: Test valable seulement en 1D mais pas genant
         IF (.not. (f(1).ge.0.) ) THEN
            print *, 'ERROR: mass flux norm is not positive!'
            print *, 'f =', f(1)
            CALL abort
         ENDIF
         
!-------------------------------------------------------------------------------
! Mass fluxes
!-------------------------------------------------------------------------------
         
         CALL thermcell_flux(ngrid,nlay,ptimestep,masse,                      &
         &                   lmin,lmax,entr_star,detr_star,                   &
         &                   f,rhobarz,zlev,zw2,fm,entr,detr)
         
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
         
         re_tpm = .false.
         DO ig=1,ngrid
            IF(lmax(ig) > lmin(ig)) THEN
               lbot(ig) = lmax(ig)
               re_tpm = .true.
            ELSE
               lbot(ig) = nlay
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Thermal plumes stacking
!-------------------------------------------------------------------------------
         
         zw2_tot(:,:)  = zw2_tot(:,:)  + zw2(:,:)
         entr_tot(:,:) = entr_tot(:,:) + entr(:,:)
         detr_tot(:,:) = detr_tot(:,:) + detr(:,:)
         fm_tot(:,:)   = fm_tot(:,:)   + fm(:,:)
         
      ENDDO
      
!===============================================================================
! Updraft fraction
!===============================================================================
      
      DO ig=1,ngrid
         fraca(ig,1) = 0.
         fraca(ig,nlay+1) = 0.
      ENDDO
      
      DO l=2,nlay
         DO ig=1,ngrid
            IF (zw2_tot(ig,l) > 0.) THEN
               fraca(ig,l) = fm_tot(ig,l) / (rhobarz(ig,l) * zw2_tot(ig,l))
            ELSE
               fraca(ig,l) = 0.
            ENDIF
         ENDDO
      ENDDO
      
!===============================================================================
! Vertical transport
!===============================================================================
      
!-------------------------------------------------------------------------------
! Calcul du transport vertical de la temperature potentielle
!-------------------------------------------------------------------------------
      
      CALL thermcell_dq(ngrid,nlay,ptimestep,fm_tot,entr_tot,detr_tot,        &
      &                 masse,zhl,zdthladj,dummy)
      
      DO l=1,nlay
         DO ig=1,ngrid
            pdtadj(ig,l) = zdthladj(ig,l) * zpopsk(ig,l)  
         ENDDO
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul du transport vertical des traceurs
!-------------------------------------------------------------------------------
      
      DO iq=1,nq
         CALL thermcell_dq(ngrid,nlay,ptimestep,fm_tot,entr_tot,detr_tot,     &
         &                 masse,pq(:,:,iq),pdqadj(:,:,iq),zqa(:,:,iq))
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul du transport vertical du moment horizontal
!-------------------------------------------------------------------------------
      
! AB: Careful, thermcell_dv2 wasn't checked! It is not sure that it works
!     correctly with the plumes stacking (zmin, wmax doesn't make sense).
      IF (dvimpl) THEN
         CALL thermcell_dv2(ngrid,nlay,ptimestep,fm_tot,entr_tot,detr_tot,    &
         &                  masse,fraca,zmax,zmin,pu,pv,pduadj,pdvadj,zua,zva)
      ELSE
         CALL thermcell_dq(ngrid,nlay,ptimestep,fm_tot,entr_tot,detr_tot,     &
         &                 masse,zu,pduadj,zua)
         CALL thermcell_dq(ngrid,nlay,ptimestep,fm_tot,entr_tot,detr_tot,     &
         &                 masse,zv,pdvadj,zva)
      ENDIF
      
      
RETURN
END
