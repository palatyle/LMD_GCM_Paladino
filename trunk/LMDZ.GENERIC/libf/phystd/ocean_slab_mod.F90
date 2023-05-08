!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ocean_slab_mod.F90,v 1.3 2008-02-04 16:24:28 fairhead Exp $
!
MODULE ocean_slab_mod
!
! This module is used for both surface ocean and sea-ice when using the slab ocean,
! "ocean=slab".
!


      use slab_ice_h
      use watercommon_h, only: T_h2O_ice_liq
      use surf_heat_transp_mod
      implicit none





  
  !LOGICAL, PRIVATE, SAVE  :: ok_slab_sic,ok_slab_heaT_h2O_ice_liqBS
  !!$OMP THREADPRIVATE(ok_slab_sic,ok_slab_heaT_h2O_ice_liqBS)
  !INTEGER, PRIVATE, SAVE                           :: slab_ekman, slab_cadj
  !!$OMP THREADPRIVATE(slab_ekman,slab_cadj)
  INTEGER, PRIVATE, SAVE                           :: lmt_pas, julien, idayvrai
  !$OMP THREADPRIVATE(lmt_pas,julien,idayvrai)
  INTEGER, PRIVATE, SAVE                           :: cpl_pas
  !$OMP THREADPRIVATE(cpl_pas)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: tmp_seaice
  !$OMP THREADPRIVATE(tmp_seaice)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: tmp_tslab_loc
  !$OMP THREADPRIVATE(tmp_tslab_loc)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: slab_bils
  !$OMP THREADPRIVATE(slab_bils)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE , SAVE  :: lmt_bils
  !$OMP THREADPRIVATE(lmt_bils)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE :: tmp_pctsrf_slab
  !$OMP THREADPRIVATE(tmp_pctsrf_slab)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: tmp_tslab
  !$OMP THREADPRIVATE(tmp_tslab)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: tmp_radsol
  !$OMP THREADPRIVATE(tmp_radsol)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: dt_hdiff
  !$OMP THREADPRIVATE(dt_hdiff)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: dt_ekman
  !$OMP THREADPRIVATE(dt_ekman)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: tmp_flux_o, tmp_flux_g
  !$OMP THREADPRIVATE(tmp_flux_o,tmp_flux_g)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE   :: slabh
  !$OMP THREADPRIVATE(slabh)
  LOGICAL, PRIVATE, SAVE                           :: check = .FALSE.
  !$OMP THREADPRIVATE(check)

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_init(ngrid,dtime, tslab_rst, seaice_rst, pctsrf_rst)

      use slab_ice_h




! Input variables
!****************************************************************************************
        
    integer,intent(in) :: ngrid ! number of atmospherci columns
    REAL, INTENT(IN)                         :: dtime
! Variables read from restart file
    REAL, DIMENSION(ngrid,noceanmx), INTENT(IN)  :: tslab_rst         
    REAL, DIMENSION(ngrid), INTENT(IN)        :: seaice_rst
    REAL, DIMENSION(ngrid), INTENT(IN) :: pctsrf_rst


! Local variables
!****************************************************************************************
    INTEGER                :: error
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'ocean_slab_intit'


    print*,'************************'
    print*,'SLAB OCEAN est actif, prenez precautions !'
    print*,'************************'    

! Lecture des parametres:
!    CALL getpar('ok_slab_sic',.true.,ok_slab_sic,'glace de mer dans slab')
!    CALL getpar('slab_ekman',0,slab_ekman,'type transport Ekman pour slab (0=rien)')
!    CALL getpar('slab_cadj',1,slab_cadj,'type ajustement conv slab 2 couches')

! Allocate variables initialize from restart fields
      ALLOCATE(tmp_tslab(ngrid,noceanmx), stat = error)
      IF (error /= 0) THEN
         abort_message='Pb allocation tmp_tslab'
         CALL abort_physic(modname,abort_message,1)
      ENDIF
      tmp_tslab(:,:) = tslab_rst(:,:)
      ALLOCATE(tmp_tslab_loc(ngrid,noceanmx), stat = error)
      IF (error /= 0) THEN
         abort_message='Pb allocation tmp_tslab_loc'
         CALL abort_physic(modname,abort_message,1)
      ENDIF
      tmp_tslab_loc(:,:) = tslab_rst(:,:)

    ALLOCATE(tmp_seaice(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_seaice'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    tmp_seaice(:) = seaice_rst(:)

    ALLOCATE(tmp_pctsrf_slab(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_pctsrf_slab'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    tmp_pctsrf_slab(:) = pctsrf_rst(:)
    
! Allocate some other variables internal in module mod_oceanslab
    ALLOCATE(tmp_radsol(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_radsol'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

    ALLOCATE(tmp_flux_o(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_flux_o'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    
    ALLOCATE(tmp_flux_g(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation tmp_flux_g'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

! a mettre un slab_bils aussi en force !!!
    ALLOCATE(slab_bils(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation slab_bils'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    slab_bils(:) = 0.0   

    ALLOCATE(dt_hdiff(ngrid,noceanmx), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation dt_hdiff'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    dt_hdiff = 0.0   

    ALLOCATE(dt_ekman(ngrid,noceanmx), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation dt_hdiff'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    dt_ekman = 0.0   


    ALLOCATE(lmt_bils(ngrid), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation lmt_bils'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    lmt_bils(:) = 0.0

    ALLOCATE(slabh(noceanmx), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation slabh'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    slabh(1)=50.
    IF (noceanmx.GE.2) slabh(2)=150.
    IF (noceanmx.GE.3) slabh(3)=2800.


!       CALL init_masquv




  END SUBROUTINE ocean_slab_init
!
!****************************************************************************************
!
 
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_ice(dtime, &
       ngrid, knindex, tsurf, radsol, &
       precip_snow, snow, evap, &
       fluxsens,tsurf_new, pctsrf_sic, &
       taux_slab,tauy_slab,icount)

      use slab_ice_h
      use callkeys_mod, only: ok_slab_sic


! Input arguments  
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: ngrid
    INTEGER, DIMENSION(ngrid), INTENT(IN) :: knindex
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(ngrid), INTENT(IN)    :: tsurf
    REAL, DIMENSION(ngrid), INTENT(IN)    :: taux_slab
    REAL, DIMENSION(ngrid), INTENT(IN)    :: tauy_slab
    REAL, DIMENSION(ngrid), INTENT(IN)    :: evap, fluxsens
    REAL, DIMENSION(ngrid), INTENT(IN)    :: precip_snow
    REAL, DIMENSION(ngrid), INTENT(IN)    :: radsol
    INTEGER, INTENT(IN)                     :: icount


!In/Output arguments
!****************************************************************************************l
    REAL, DIMENSION(ngrid), INTENT(INOUT)          :: snow

!    REAL, DIMENSION(ngridmx), INTENT(INOUT)          :: agesno


! Output arguments
!****************************************************************************************
!    REAL, DIMENSION(ngridmx), INTENT(OUT)            :: alb1_new  ! new albedo in visible SW interval
!    REAL, DIMENSION(ngridmx), INTENT(OUT)            :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(ngrid), INTENT(OUT)            :: tsurf_new     
    REAL, DIMENSION(ngrid), INTENT(OUT)            :: pctsrf_sic

! Local variables
!****************************************************************************************
    INTEGER                                 :: i
    REAL, DIMENSION(ngrid)                :: cal, beta, dif_grnd, capsol
    REAL, DIMENSION(ngrid)                :: alb_neig, alb_ice!, tsurf_temp
    REAL, DIMENSION(ngrid)                :: soilcap, soilflux
    REAL, DIMENSION(ngrid)                :: zfra
    REAL                                    :: snow_evap, fonte 
    REAL, DIMENSION(ngrid,noceanmx)        :: tslab
    REAL, DIMENSION(ngrid)                :: seaice,tsea_ice ! glace de mer (kg/m2)
    REAL, DIMENSION(ngrid)                :: pctsrf_new


!*****************************************************************************


! Initialization of output variables
!    alb1_new(:) = 0.0

!******************************************************************************
! F. Codron: 3 cas 
! -Glace interactive, quantité seaice : T glace suit modèle simple
! -Pas de glace: T oce peut descendre en dessous de 0°C sans geler
!*****************************************************************************

       pctsrf_new(:)=tmp_pctsrf_slab(:)
       tmp_radsol(:)=radsol(:)
       tmp_flux_o(:)=fluxsens(:)

    DO i = 1, ngrid
       tsurf_new(i) = tsurf(knindex(i))
       seaice(i)=tmp_seaice(knindex(i))



       IF (pctsrf_new(knindex(i)) < EPSFRA) THEN
          snow(i) = 0.0
          tsurf_new(i) = T_h2O_ice_liq - 1.8
          !IF (soil_model) tsoil(i,:) = T_h2O_ice_liq -1.8
       ENDIF
    ENDDO
    tmp_flux_g(:) = 0.0
    tsea_ice(:)=T_h2O_ice_liq-1.8
! Calculs flux glace/mer et glace/air
    IF (ok_slab_sic) THEN   
!*********************************
! Calcul de beta, heat capacity
!********************************* 
!    CALL calbeta(dtime, is_sic, knon, snow, qsol, beta, capsol, dif_grnd)
    
!    IF ((soil_model)) THEN 

!    ELSE 
!       dif_grnd = 1.0 / tau_gl
!       cal = RCPD * calice
!       WHERE (snow > 0.0) cal = RCPD * calsno 
!    ENDIF
!    tsurf_temp = tsurf_new
!    beta = 1.0


! **********************************************
! Evolution avec glace interactive:
! *************************************
!      tsurf_new=tsurf_temp
      DO i = 1, ngrid
        IF (pctsrf_new(knindex(i)) .GT. epsfra) THEN
! *************************************
! + Calcul Flux entre glace et océan 
! *************************************
          tmp_flux_g(knindex(i))=(tsurf_new(i)-(T_h2O_ice_liq-1.8)) &
            * ice_cond*ice_den/seaice(i)

! ****************************************
! + Calcul nouvelle température de la glace 
! ****************************************
          tsurf_new(i)=tsurf_new(i)+2*(radsol(i)+fluxsens(i) &
             -tmp_flux_g(knindex(i))) &
             /(ice_cap*seaice(i))*dtime

! ***************************************
! + Precip and evaporation of snow and ice
! ***************************************
! Add precip
          IF (precip_snow(i).GT.0.) THEN
            snow(i)=snow(i)+precip_snow(i)*dtime
          ENDIF
! Evaporation 
          snow_evap=0.
          IF (evap(i).GT.0.) THEN
            snow_evap = MIN (snow(i) / dtime, evap(i))
            snow(i) = snow(i) - snow_evap*dtime
            snow(i) = MAX(0.0, snow(i))
            seaice(i)=MAX(0.0,seaice(i)-(evap(i)-snow_evap)*dtime)
          ENDIF
        
! *****************************************
! + Fonte neige & glace par le haut
! *****************************************
! Snow melt
          IF ((snow(i) > epsfra) .AND. (tsurf_new(i) >= T_h2O_ice_liq)) THEN
            fonte=MIN(MAX((tsurf_new(i)-T_h2O_ice_liq)*seaice(i) &
                  /2.*ice_cap/ice_lat,0.0),snow(i))
            snow(i) = MAX(0.,(snow(i)-fonte))
            tsurf_new(i)=tsurf_new(i)-fonte*2*ice_lat/(seaice(i)*ice_cap)
          ENDIF
          snow(i)=MIN(snow(i),3000.)
! Ice melt
          IF ((seaice(i) > epsfra) .AND. (tsurf_new(i) >= T_h2O_ice_liq)) THEN       
             fonte=seaice(i)*ice_cap*(tsurf_new(i)-T_h2O_ice_liq) &
                  /(2*ice_lat+ice_cap*1.8)
             CALL icemelt(fonte,pctsrf_new(knindex(i)),seaice(i))
             tmp_flux_g(knindex(i))=tmp_flux_g(knindex(i)) &
                                 +fonte*ice_lat/dtime
             tsurf_new(i)=T_h2O_ice_liq
          ENDIF  
          tmp_seaice(knindex(i))=seaice(i)
        ENDIF!(pctsrf_new(knindex(i)) .GT. epsfra)
          tsea_ice(knindex(i))=tsurf_new(i)

      ENDDO
    ENDIF
! ******************************************
! CALL interfoce:
! cumul/calcul nouvelle T_h2O_ice_liqce
! fonte/formation de glace en dessous
! ******************************************   
    tslab = tmp_tslab
    
    CALL interfoce_slab(ngrid, dtime, &
         tmp_radsol, tmp_flux_o, tmp_flux_g, tmp_pctsrf_slab, &
         tslab, tsea_ice, pctsrf_new,taux_slab,tauy_slab,icount)
    
    tmp_pctsrf_slab(:)=pctsrf_new(:)
!    tmp_pctsrf_slab(:,is_oce)=1.0-tmp_pctsrf_slab(:) &
!            -tmp_pctsrf_slab(:,is_lic)-tmp_pctsrf_slab(:,is_ter)

    DO i=1, ngrid
       tmp_tslab(knindex(i),:)=tslab(knindex(i),:)
       seaice(i)=tmp_seaice(knindex(i))
       tsurf_new(i)=tsea_ice(knindex(i))

    ENDDO

! ****************************
! calcul new albedo
! ****************************
! Albedo neige
!    CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
!    WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
! Fraction neige (hauteur critique 45kg/m2~15cm)
!    zfra(1:knon) = MAX(0.0,MIN(1.0,snow(1:knon)/45.0))
! Albedo glace
!    IF (ok_slab_sicOBS) THEN
!       alb_ice=0.6
!    ELSE
!       alb_ice(1:knon)=alb_ice_max-(alb_ice_max-alb_ice_min) &
!         *exp(-seaice(1:knon)/h_alb_ice)
!    ENDIF
!Albedo final
!    alb1_new(1 : knon) = alb_neig(1 : knon) *zfra(1:knon) + & 
!         alb_ice * (1.0-zfra(1:knon))
!    alb2_new(:) = alb1_new(:)

        
! ********************************
! Return the fraction of sea-ice
! ********************************
    pctsrf_sic(:) =  tmp_pctsrf_slab(:)


  END SUBROUTINE ocean_slab_ice


! 
!***************************************************************************
!
  SUBROUTINE interfoce_slab(ngrid, dtime, &
       radsol, fluxo, fluxg, pctsrf, &
       tslab, tsea_ice, pctsrf_slab, &
       taux_slab, tauy_slab,icount)
!
! Cette routine calcule la temperature d'un slab ocean, la glace de mer 
! et les pourcentages de la maille couverte par l'ocean libre et/ou 
! la glace de mer pour un "slab" ocean de 50m
!
! Conception: Laurent Li
! Re-ecriture + adaptation LMDZ4: I. Musat
! Transport, nouveau modèle glace, 2 couches: F.Codron
!
! input:
!   klon         nombre T_h2O_ice_liqtal de points de grille
!   itap         numero du pas de temps
!   dtime        pas de temps de la physique (en s)
!   ijour        jour dans l'annee en cours
!   radsol       rayonnement net au sol (LW + SW)
!   fluxo        flux turbulent (sensible + latent) sur les mailles oceaniques 
!   fluxg        flux de conduction entre la surface de la glace de mer et l'ocean
!   pctsrf       tableau des pourcentages de surface de chaque maille
! output: 
!   tslab        temperature de l'ocean libre
!   tsea_ice         temperature de la glace (surface)
!   pctsrf_slab  "pourcentages" (valeurs entre 0. et 1.) surfaces issus du slab

    use slab_ice_h 
    use callkeys_mod, only: ok_slab_sic,ok_slab_heat_transp

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                       :: ngrid,icount
!    INTEGER, INTENT(IN)                       :: itap
    REAL, INTENT(IN)                          :: dtime       ! not used
!    INTEGER, INTENT(IN)                       :: ijour
    REAL, DIMENSION(ngrid), INTENT(IN)         :: radsol
    REAL, DIMENSION(ngrid), INTENT(IN)         :: fluxo
    REAL, DIMENSION(ngrid), INTENT(IN)         :: fluxg
    REAL, DIMENSION(ngrid), INTENT(IN)  :: pctsrf
    REAL, DIMENSION(ngrid), INTENT(IN)  :: taux_slab
    REAL, DIMENSION(ngrid), INTENT(IN)  :: tauy_slab

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(ngrid,noceanmx), INTENT(OUT) :: tslab
    REAL, DIMENSION(ngrid), INTENT(INOUT)       :: pctsrf_slab
    REAL, DIMENSION(ngrid), INTENT(INOUT)       :: tsea_ice

! Local variables
!****************************************************************************************
    REAL                    :: fonte,t_cadj
    INTEGER                 :: i,k,kt,kb
    REAL                    :: zz, za, zb
    REAL                    :: cyang ! capacite calorifique slab J/(m2 K)
    REAL, PARAMETER         :: tau_conv=432000. ! temps ajust conv (5 jours)
    REAL, DIMENSION(ngrid,noceanmx) :: dtdiff_loc,dtekman_loc



!***************************************************
! Capacite calorifique de la couche de surface
    cyang=capcalocean!slabh(1)*4.228e+06
    cpl_pas=1!4*iradia

!
! lecture du bilan au sol lmt_bils issu d'une simulation forcee en debut de journee
!!    julien = MOD(ijour,360)
!!    IF (ok_slab_heaT_h2O_ice_liqBS .and. (MOD(itap,lmt_pas) .EQ. 1)) THEN  
!!       ! 1er pas de temps de la journee
!!       idayvrai = ijour
!!       CALL condsurf(julien,idayvrai, lmt_bils)
!!    ENDIF

! ----------------------------------------------
! A chaque pas de temps: cumul flux de chaleur
! ----------------------------------------------
! bilan du flux de chaleur dans l'ocean :
!
    DO i = 1, ngrid
       za = radsol(i) + fluxo(i)
       zb = fluxg(i)
!       IF(((1-pctsrf(i)).GT.epsfra).OR. &
!            (pctsrf(i).GT.epsfra)) THEN
          slab_bils(i)=slab_bils(i)+(za*(1-pctsrf(i)) &
               +zb*pctsrf(i))/ cpl_pas
!       ENDIF

    ENDDO !klon

! ---------------------------------------------
! T_h2O_ice_liqus les cpl_pas: update de tslab et seaice
! ---------------------------------------------

    IF (mod(icount-1,cpl_pas).eq.0) THEN !fin de journee
!
! ---------------------------------------------
! Transport de chaleur par circulation
! Decoupage par pas de temps pour stabilite numérique 
! (diffusion, schema centre pour advection)
! ---------------------------------------------
      dt_hdiff(:,:)=0.
      dt_ekman(:,:)=0.

      IF (ok_slab_heat_transp) THEN
!       DO i=1,cpl_pas
!  Transport diffusif
!         IF (ok_soil_hdiff) THEN
             CALL divgrad_phy(ngrid,noceanmx,tmp_tslab_loc,dtdiff_loc)
             dtdiff_loc=dtdiff_loc*soil_hdiff*50./SUM(slabh)!*100.
 !           dtdiff_loc(:,1)=dtdiff_loc(:,1)*soil_hdiff*50./SUM(slabh)*0.8
 !           dtdiff_loc(:,2)=dtdiff_loc(:,2)*soil_hdiff*50./SUM(slabh)*0.2
             dt_hdiff=dt_hdiff+dtdiff_loc
             tmp_tslab_loc=tmp_tslab_loc+dtdiff_loc*dtime
!         END IF

! Calcul  de transport par Ekman

          CALL slab_ekman2(ngrid,taux_slab,tauy_slab,tslab,dtekman_loc)


!        END SELECT
!        IF (slab_ekman.GT.0) THEN
          do k=1,noceanmx
            dtekman_loc(:,k)=dtekman_loc(:,k)/(slabh(k)*1000.)!*0.
          enddo
          dt_ekman(:,:)=dt_ekman(:,:)+dtekman_loc(:,:)
          tmp_tslab_loc=tmp_tslab_loc+dtekman_loc*dtime
!        ENDIF
!      ENDDO ! time splitting 1,cpl_pas     
!      IF (slab_ekman.GT.0) THEN
!	  taux_slab(:)=0.
!	  tauy_slab(:)=0.
!      ENDIF

!       print*, 'slab_bils=',slab_bils(1)


      ENDIF!(ok_slab_heat_transp)

      DO i = 1, ngrid
!      IF (((1-pctsrf(i)).GT.epsfra).OR. &
!             (pctsrf(i).GT.epsfra)) THEN
! ---------------------------------------------
! Ajout des flux de chaleur dans l'océan
! ---------------------------------------

!print*, 'T_h2O_ice_liqcean1=',tmp_tslab_loc(i,1)
        tmp_tslab_loc(i,1) = tmp_tslab_loc(i,1) + &
             slab_bils(i)/cyang*dtime*cpl_pas

!print*, 'capcalocean=',capcalocean
!print*, 'cyang=',cyang
!print*, 'dT_h2O_ice_liqcean=',slab_bils(i)/cyang*dtime
!print*, 'T_h2O_ice_liqcean2=',tmp_tslab_loc(i,1)


! on remet l'accumulation a 0
        slab_bils(i) = 0.
! ---------------------------------------------
! Glace interactive 
! ---------------------------------------------
        IF (ok_slab_sic) THEN
! Fondre la glace si température > 0. 
! -----------------------------------
          IF ((tmp_tslab_loc(i,1).GT.T_h2O_ice_liq-1.8) .AND. (tmp_seaice(i).GT.0.0)) THEN
            fonte=(tmp_tslab_loc(i,1)-T_h2O_ice_liq+1.8)*cyang &
                /(ice_lat+ice_cap/2*(T_h2O_ice_liq-1.8-tsea_ice(i)))
            CALL icemelt(fonte,pctsrf_slab(i),tmp_seaice(i))
                 tmp_tslab_loc(i,1)=T_h2O_ice_liq-1.8+fonte*ice_lat/cyang
          ENDIF 
! fabriquer de la glace si congelation atteinte:
! ----------------------------------------------
          IF (tmp_tslab_loc(i,1).LT.(T_h2O_ice_liq-1.8)) THEN

            IF (tmp_seaice(i).LT.h_ice_min) THEN
! Creation glace nouvelle
!             IF (pctsrf_slab(i).LT.(epsfra)) THEN
                 fonte=(T_h2O_ice_liq-1.8-tmp_tslab_loc(i,1))*cyang/ice_lat
              IF (fonte.GT.h_ice_min*ice_frac_min) THEN
                tmp_seaice(i)=MIN(h_ice_thin,fonte/ice_frac_min)
                tmp_seaice(i)=MAX(tmp_seaice(i),fonte/ice_frac_max)
!             IF (fonte.GT.0) THEN
!               tmp_seaice(i)=fonte
                tsea_ice(i)=T_h2O_ice_liq-1.8
                pctsrf_slab(i)=(1-pctsrf_slab(i))*fonte/tmp_seaice(i)
!                pctsrf_slab(i)=1.
                tmp_tslab_loc(i,1)=T_h2O_ice_liq-1.8
              ENDIF
            ELSE  
! glace déjà présente 
! Augmenter epaisseur
              fonte=(T_h2O_ice_liq-1.8-tmp_tslab_loc(i,1))*cyang &
                   /(ice_lat+ice_cap/2*(T_h2O_ice_liq-1.8-tsea_ice(i)))           
              zz=tmp_seaice(i)
              tmp_seaice(i)=MAX(zz,MIN(h_ice_thick,fonte+zz))
! Augmenter couverture (oce libre et h>h_thick) 
              za=1-pctsrf_slab(i)
              zb=pctsrf_slab(i)
              fonte=fonte*za+MAX(0.0,fonte+zz-tmp_seaice(i))*zb
              pctsrf_slab(i)=MIN(zb+fonte/tmp_seaice(i), &
                                    (za+zb)*ice_frac_max)
              fonte=MAX(0.0,fonte-(pctsrf_slab(i)-zb)*tmp_seaice(i))
! Augmenter epaisseur si couverture complete
              tmp_seaice(i)=tmp_seaice(i)+fonte/pctsrf_slab(i)
              tmp_tslab_loc(i,1) = T_h2O_ice_liq-1.8
            ENDIF ! presence glace
          ENDIF !congelation
! vérifier limites de hauteur de glace
          IF(tmp_seaice(i).GT.h_ice_min) THEN
            tmp_seaice(i) = MIN(tmp_seaice(i),h_ice_max)
          ELSE
            tmp_seaice(i) = 0. 
            pctsrf_slab(i)=0.
          ENDIF! (tmp_seaice(i).GT.h_ice_min)
       
        ENDIF !(ok_slab_sic) Glace interactive
! ----------------------------------
! Ajustement convectif si > 1 layers
! ----------------------------------

        IF ((noceanmx.GT.1)) THEN
          DO kt=1,noceanmx-1
            kb=kt           
            DO k=kt+1,noceanmx !look for instability
              IF (tmp_tslab_loc(i,k).GT.tmp_tslab_loc(i,kt)) THEN
                kb=k
              ENDIF
            ENDDO
            IF (kb.GT.kt) THEN !ajust conv
!               IF (slab_cadj.EQ.1) THEN
!             :   t_cadj=SUM(tmp_tslab_loc(i,kt:kb)*slabh(kt:kb))/SUM(slabh(kt:kb))
!                DO k=kt,kb
!                  tmp_tslab_loc(i,k)=t_cadj
!                ENDDO
!              ELSEIF (slab_cadj.EQ.2) THEN
                t_cadj=(tmp_tslab_loc(i,kb)-tmp_tslab_loc(i,kt))*dtime/tau_conv*cpl_pas
                tmp_tslab_loc(i,kt)=tmp_tslab_loc(i,kt)+t_cadj
                tmp_tslab_loc(i,kb)=tmp_tslab_loc(i,kb)-t_cadj*slabh(kt)/slabh(kb)
              !ENDIF 
            ENDIF
          ENDDO 
        ENDIF !ajust 2 layers
!      ENDIF !pctsrf
      ENDDO !klon

! On met a jour la temperature et la glace
    tslab  = tmp_tslab_loc


    ENDIF !(mod(icount-1,cpl_pas).eq.0)

  END SUBROUTINE interfoce_slab
!
!******************************************************************************  
  SUBROUTINE icemelt(fonte,pctsrf,seaice)

      use slab_ice_h



     REAL, INTENT(INOUT)  :: pctsrf
     REAL , INTENT(INOUT)   ::fonte,seaice !kg/m2
     REAL :: hh !auxilliary
     REAL :: ff !auxilliary


! ice gt h_ice_thick: decrease thickness up T_h2O_ice_liq h_ice_thick
     IF (seaice.GT.h_ice_thick) THEN
        hh=seaice
!        ff=fonte
        seaice=max(h_ice_thick,hh-fonte)    
        fonte=max(0.0,fonte+h_ice_thick-hh)

!        seaice=max(0.,hh-fonte)    
!        fonte=max(0.0,fonte-(seaice-hh))
!     IF (seaice.LT.epsfra) THEN
!        pctsrf=0.
!        seaice=0.
!        fonte=ff-hh
!     ENDIF

     ENDIF
! ice gt h_ice_thin: partially decrease thickness
     IF ((seaice.GE.h_ice_thin).AND.(fonte.GT.0.0)) THEN
       hh=seaice
       seaice=MAX(hh-0.6*fonte,h_ice_thin)
       fonte=MAX(0.0,fonte-hh+seaice)
     ENDIF
! use rest T_h2O_ice_liq decrease area
       hh=pctsrf
       pctsrf=MIN(hh,MAX(0.0,hh-fonte/seaice))
       fonte=MAX(0.0,fonte-(hh-pctsrf)*seaice)

    END SUBROUTINE icemelt

!****************************************************************************************
!
  SUBROUTINE ocean_slab_final!(tslab_rst, seaice_rst)

! This subroutine will send T_h2O_ice_liq phyredem the variables concerning the slab 
! ocean that should be written T_h2O_ice_liq restart file.

!****************************************************************************************

!    REAL, DIMENSION(ngridmx,noceanmx), INTENT(OUT) :: tslab_rst
!    REAL, DIMENSION(ngridmx), INTENT(OUT) :: seaice_rst

!****************************************************************************************
! Set the output variables
!    tslab_rst(:,:)  = tmp_tslab(:,:)
!    tslab_rst(:)  = tmp_tslab_loc(:)
!    seaice_rst(:) = tmp_seaice(:)

! Deallocation of all variables in module
    DEALLOCATE(tmp_tslab, tmp_tslab_loc, tmp_pctsrf_slab)
    DEALLOCATE(tmp_seaice, tmp_radsol, tmp_flux_o, tmp_flux_g)
    DEALLOCATE(slab_bils, lmt_bils)
    DEALLOCATE(dt_hdiff)

  END SUBROUTINE ocean_slab_final
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_get_vars(ngrid,tslab_loc, seaice_loc, flux_o_loc, flux_g_loc, &
       dt_hdiff_loc,dt_ekman_loc)
  
! "Get some variables from module ocean_slab_mod"
! This subroutine prints variables T_h2O_ice_liq a external routine

    INTEGER, INTENT(IN)                     :: ngrid
    REAL, DIMENSION(ngrid,noceanmx), INTENT(OUT) :: tslab_loc
    REAL, DIMENSION(ngrid), INTENT(OUT) :: seaice_loc
    REAL, DIMENSION(ngrid), INTENT(OUT) :: flux_o_loc
    REAL, DIMENSION(ngrid), INTENT(OUT) :: flux_g_loc
    REAL, DIMENSION(ngrid,noceanmx), INTENT(OUT) :: dt_hdiff_loc
    REAL, DIMENSION(ngrid,noceanmx), INTENT(OUT) :: dt_ekman_loc
    INTEGER :: i


! Set the output variables
    tslab_loc=0.
    dt_hdiff_loc=0.
    dt_ekman_loc=0.
    tslab_loc  = tmp_tslab
    seaice_loc(:) = tmp_seaice(:)
    flux_o_loc(:) = tmp_flux_o(:)
    flux_g_loc(:) = tmp_flux_g(:)
    DO i=1,noceanmx
      dt_hdiff_loc(:,i) = dt_hdiff(:,i)*slabh(i)*1000.*4228. !Convert en W/m2
      dt_ekman_loc(:,i) = dt_ekman(:,i)*slabh(i)*1000.*4228. 
    ENDDO
  
  

  END SUBROUTINE ocean_slab_get_vars
!
!****************************************************************************************
!
END MODULE ocean_slab_mod
