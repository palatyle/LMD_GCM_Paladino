! $Id: phys_output_mod.F90 1473 2011-01-12 15:07:43Z fairhead $
!
! Abderrahmane 12 2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Ecreture des Sorties du modele dans les fichiers Netcdf :
! histmth.nc : moyennes mensuelles
! histday.nc : moyennes journalieres
! histhf.nc  : moyennes toutes les 3 heures
! histins.nc : valeurs instantanees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE phys_output_mod 

  IMPLICIT NONE

  private histdef2d, histdef3d, conf_physoutputs


   integer, parameter                           :: nfiles = 5
   logical, dimension(nfiles), save             :: clef_files
   integer, dimension(nfiles), save             :: lev_files
   integer, dimension(nfiles), save             :: nid_files
!!$OMP THREADPRIVATE(clef_files, lev_files,nid_files)
 
   integer, dimension(nfiles), private, save :: nhorim, nvertm
   integer, dimension(nfiles), private, save :: nvertap, nvertbp, nvertAlt
!   integer, dimension(nfiles), private, save :: nvertp0
   real, dimension(nfiles), private, save                :: zoutm
   real,                    private, save                :: zdtime
   CHARACTER(len=20), dimension(nfiles), private, save   :: type_ecri
!$OMP THREADPRIVATE(nhorim, nvertm, zoutm,zdtime,type_ecri)

!   integer, save                     :: nid_hf3d 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition pour chaque variable du niveau d ecriture dans chaque fichier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/ histmth, histday, histhf, histins /),'!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, private:: levmin(nfiles) = 1
  integer, private:: levmax(nfiles)

  TYPE ctrl_out
   integer,dimension(5) :: flag
   character(len=20)     :: name
  END TYPE ctrl_out

!!! Comosentes de la coordonnee sigma-hybride
!!! Ap et Bp
  type(ctrl_out),save :: o_Ahyb         = ctrl_out((/ 1, 1, 1, 1, 1 /), 'Ap')
  type(ctrl_out),save :: o_Bhyb         = ctrl_out((/ 1, 1, 1, 1, 1 /), 'Bp')
  type(ctrl_out),save :: o_Alt          = ctrl_out((/ 1, 1, 1, 1, 1 /), 'Alt')

!!! 1D
  type(ctrl_out),save :: o_phis         = ctrl_out((/ 1, 1, 10, 5, 1 /), 'phis') 
  type(ctrl_out),save :: o_aire         = ctrl_out((/ 1, 1, 10,  10, 1 /),'aire')
  type(ctrl_out),save :: o_contfracATM  = ctrl_out((/ 10, 1,  1, 10, 10 /),'contfracATM')
  type(ctrl_out),save :: o_contfracOR   = ctrl_out((/ 10, 1,  1, 10, 10 /),'contfracOR')
  type(ctrl_out),save :: o_aireTER      = ctrl_out((/ 10, 10, 1, 10, 10 /),'aireTER')
  
!!! 2D

  type(ctrl_out),save :: o_flat         = ctrl_out((/ 5, 1, 10, 5, 1 /),'flat')
  type(ctrl_out),save :: o_slp          = ctrl_out((/ 1, 1, 1, 10, 1 /),'slp')
  type(ctrl_out),save :: o_tsol         = ctrl_out((/ 1, 1, 1, 5, 1 /),'tsol')
  type(ctrl_out),save :: o_t2m          = ctrl_out((/ 1, 1, 1, 5, 1 /),'t2m')
  type(ctrl_out),save :: o_t2m_min      = ctrl_out((/ 1, 1, 10, 10, 10 /),'t2m_min')
  type(ctrl_out),save :: o_t2m_max      = ctrl_out((/ 1, 1, 10, 10, 10 /),'t2m_max')
  type(ctrl_out),save,dimension(4) :: o_t2m_srf      = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'t2m_ter'), &
                                                 ctrl_out((/ 10, 4, 10, 10, 10 /),'t2m_lic'), &
                                                 ctrl_out((/ 10, 4, 10, 10, 10 /),'t2m_oce'), &
                                                 ctrl_out((/ 10, 4, 10, 10, 10 /),'t2m_sic') /)

  type(ctrl_out),save :: o_wind10m      = ctrl_out((/ 1, 1, 1, 10, 10 /),'wind10m')
  type(ctrl_out),save :: o_wind10max    = ctrl_out((/ 10, 1, 10, 10, 10 /),'wind10max')
  type(ctrl_out),save :: o_sicf         = ctrl_out((/ 1, 1, 10, 10, 10 /),'sicf')
  type(ctrl_out),save :: o_q2m          = ctrl_out((/ 1, 1, 1, 5, 1 /),'q2m')
  type(ctrl_out),save :: o_u10m         = ctrl_out((/ 1, 1, 1, 5, 1 /),'u10m')
  type(ctrl_out),save :: o_v10m         = ctrl_out((/ 1, 1, 1, 5, 1 /),'v10m')
  type(ctrl_out),save :: o_psol         = ctrl_out((/ 1, 1, 1, 5, 1 /),'psol')
  type(ctrl_out),save :: o_qsurf        = ctrl_out((/ 1, 10, 10, 10, 10 /),'qsurf')

  type(ctrl_out),save,dimension(4) :: o_u10m_srf     = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'u10m_ter'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'u10m_lic'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'u10m_oce'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'u10m_sic') /)

  type(ctrl_out),save,dimension(4) :: o_v10m_srf     = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'v10m_ter'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'v10m_lic'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'v10m_oce'), &
                                              ctrl_out((/ 10, 4, 10, 10, 10 /),'v10m_sic') /)

  type(ctrl_out),save :: o_qsol         = ctrl_out((/ 1, 10, 10, 10, 10 /),'qsol')

  type(ctrl_out),save :: o_ndayrain     = ctrl_out((/ 1, 10, 10, 10, 10 /),'ndayrain')
  type(ctrl_out),save :: o_precip       = ctrl_out((/ 1, 1, 1, 5, 10 /),'precip')
  type(ctrl_out),save :: o_plul         = ctrl_out((/ 1, 1, 1, 10, 10 /),'plul')

  type(ctrl_out),save :: o_pluc         = ctrl_out((/ 1, 1, 1, 5, 10 /),'pluc')
  type(ctrl_out),save :: o_snow         = ctrl_out((/ 1, 1, 10, 5, 10 /),'snow') 
  type(ctrl_out),save :: o_evap         = ctrl_out((/ 1, 1, 10, 10, 10 /),'evap')
  type(ctrl_out),save,dimension(4) :: o_evap_srf     = (/ ctrl_out((/ 1, 1, 10, 10, 10 /),'evap_ter'), &
                                           ctrl_out((/ 1, 1, 10, 10, 10 /),'evap_lic'), &
                                           ctrl_out((/ 1, 1, 10, 10, 10 /),'evap_oce'), &
                                           ctrl_out((/ 1, 1, 10, 10, 10 /),'evap_sic') /)
  type(ctrl_out),save :: o_msnow       = ctrl_out((/ 1, 10, 10, 10, 10 /),'msnow')
  type(ctrl_out),save :: o_fsnow       = ctrl_out((/ 1, 10, 10, 10, 10 /),'fsnow')

  type(ctrl_out),save :: o_tops         = ctrl_out((/ 1, 1, 10, 10, 10 /),'tops')
  type(ctrl_out),save :: o_tops0        = ctrl_out((/ 1, 5, 10, 10, 10 /),'tops0')
  type(ctrl_out),save :: o_topl         = ctrl_out((/ 1, 1, 10, 5, 10 /),'topl')
  type(ctrl_out),save :: o_topl0        = ctrl_out((/ 1, 5, 10, 10, 10 /),'topl0')
  type(ctrl_out),save :: o_SWupTOA      = ctrl_out((/ 1, 4, 10, 10, 10 /),'SWupTOA')
  type(ctrl_out),save :: o_SWupTOAclr   = ctrl_out((/ 1, 4, 10, 10, 10 /),'SWupTOAclr')
  type(ctrl_out),save :: o_SWdnTOA      = ctrl_out((/ 1, 4, 10, 10, 10 /),'SWdnTOA')
  type(ctrl_out),save :: o_SWdnTOAclr   = ctrl_out((/ 1, 4, 10, 10, 10 /),'SWdnTOAclr')
  type(ctrl_out),save :: o_nettop       = ctrl_out((/ 1, 4, 10, 10, 10 /),'nettop')

  type(ctrl_out),save :: o_SWup200      = ctrl_out((/ 1, 10, 10, 10, 10 /),'SWup200')
  type(ctrl_out),save :: o_SWup200clr   = ctrl_out((/ 10, 1, 10, 10, 10 /),'SWup200clr')
  type(ctrl_out),save :: o_SWdn200      = ctrl_out((/ 1, 10, 10, 10, 10 /),'SWdn200')
  type(ctrl_out),save :: o_SWdn200clr   = ctrl_out((/ 10, 1, 10, 10, 10 /),'SWdn200clr')

! arajouter
!  type(ctrl_out),save :: o_LWupTOA     = ctrl_out((/ 1, 4, 10, 10, 10 /),'LWupTOA')
!  type(ctrl_out),save :: o_LWupTOAclr  = ctrl_out((/ 1, 4, 10, 10, 10 /),'LWupTOAclr')
!  type(ctrl_out),save :: o_LWdnTOA     = ctrl_out((/ 1, 4, 10, 10, 10 /),'LWdnTOA')
!  type(ctrl_out),save :: o_LWdnTOAclr  = ctrl_out((/ 1, 4, 10, 10, 10 /),'LWdnTOAclr')

  type(ctrl_out),save :: o_LWup200      = ctrl_out((/ 1, 10, 10, 10, 10 /),'LWup200')
  type(ctrl_out),save :: o_LWup200clr   = ctrl_out((/ 1, 10, 10, 10, 10 /),'LWup200clr')
  type(ctrl_out),save :: o_LWdn200      = ctrl_out((/ 1, 10, 10, 10, 10 /),'LWdn200')
  type(ctrl_out),save :: o_LWdn200clr   = ctrl_out((/ 1, 10, 10, 10, 10 /),'LWdn200clr')
  type(ctrl_out),save :: o_sols         = ctrl_out((/ 1, 1, 10, 10, 10 /),'sols')
  type(ctrl_out),save :: o_sols0        = ctrl_out((/ 1, 5, 10, 10, 10 /),'sols0')
  type(ctrl_out),save :: o_soll         = ctrl_out((/ 1, 1, 10, 10, 10 /),'soll')
  type(ctrl_out),save :: o_soll0        = ctrl_out((/ 1, 5, 10, 10, 10 /),'soll0')
  type(ctrl_out),save :: o_radsol       = ctrl_out((/ 1, 1, 10, 10, 10 /),'radsol')
  type(ctrl_out),save :: o_SWupSFC      = ctrl_out((/ 1, 4, 10, 5, 10 /),'SWupSFC')
  type(ctrl_out),save :: o_SWupSFCclr   = ctrl_out((/ 1, 4, 10, 10, 10 /),'SWupSFCclr')
  type(ctrl_out),save :: o_SWdnSFC      = ctrl_out((/ 1, 1, 10, 5, 10 /),'SWdnSFC') 
  type(ctrl_out),save :: o_SWdnSFCclr   = ctrl_out((/ 1, 4, 10, 5, 10 /),'SWdnSFCclr')
  type(ctrl_out),save :: o_LWupSFC      = ctrl_out((/ 1, 4, 10, 10, 10 /),'LWupSFC')
  type(ctrl_out),save :: o_LWupSFCclr   = ctrl_out((/ 1, 4, 10, 5, 10 /),'LWupSFCclr')
  type(ctrl_out),save :: o_LWdnSFC      = ctrl_out((/ 1, 4, 10, 5, 10 /),'LWdnSFC')
  type(ctrl_out),save :: o_LWdnSFCclr   = ctrl_out((/ 1, 4, 10, 5, 10 /),'LWdnSFCclr')
  type(ctrl_out),save :: o_bils         = ctrl_out((/ 1, 2, 10, 5, 10 /),'bils')
  type(ctrl_out),save :: o_sens         = ctrl_out((/ 1, 1, 10, 5, 10 /),'sens')
  type(ctrl_out),save :: o_fder         = ctrl_out((/ 1, 2, 10, 10, 10 /),'fder')
  type(ctrl_out),save :: o_ffonte       = ctrl_out((/ 1, 10, 10, 10, 10 /),'ffonte')
  type(ctrl_out),save :: o_fqcalving    = ctrl_out((/ 1, 10, 10, 10, 10 /),'fqcalving')
  type(ctrl_out),save :: o_fqfonte      = ctrl_out((/ 1, 10, 10, 10, 10 /),'fqfonte')

  type(ctrl_out),save :: o_taux         = ctrl_out((/ 1, 10, 10, 10, 10 /),'taux')
  type(ctrl_out),save :: o_tauy         = ctrl_out((/ 1, 10, 10, 10, 10 /),'tauy')
  type(ctrl_out),save,dimension(4) :: o_taux_srf     = (/ ctrl_out((/ 1, 4, 10, 10, 10 /),'taux_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'taux_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'taux_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'taux_sic') /)

  type(ctrl_out),save,dimension(4) :: o_tauy_srf     = (/ ctrl_out((/ 1, 4, 10, 10, 10 /),'tauy_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tauy_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tauy_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tauy_sic') /)


  type(ctrl_out),save,dimension(4) :: o_pourc_srf    = (/ ctrl_out((/ 1, 4, 10, 10, 10 /),'pourc_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'pourc_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'pourc_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'pourc_sic') /)     

  type(ctrl_out),save,dimension(4) :: o_fract_srf    = (/ ctrl_out((/ 1, 4, 10, 10, 10 /),'fract_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'fract_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'fract_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'fract_sic') /)

  type(ctrl_out),save,dimension(4) :: o_tsol_srf     = (/ ctrl_out((/ 1, 4, 10, 10, 10 /),'tsol_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tsol_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tsol_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 10, 10 /),'tsol_sic') /)

  type(ctrl_out),save,dimension(4) :: o_sens_srf     = (/ ctrl_out((/ 1, 4, 10, 1, 10 /),'sens_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'sens_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'sens_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'sens_sic') /)

  type(ctrl_out),save,dimension(4) :: o_lat_srf      = (/ ctrl_out((/ 1, 4, 10, 1, 10 /),'lat_ter'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'lat_lic'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'lat_oce'), &
                                                 ctrl_out((/ 1, 4, 10, 1, 10 /),'lat_sic') /)

  type(ctrl_out),save,dimension(4) :: o_flw_srf      = (/ ctrl_out((/ 1, 10, 10, 10, 10 /),'flw_ter'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'flw_lic'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'flw_oce'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'flw_sic') /)
                                                 
  type(ctrl_out),save,dimension(4) :: o_fsw_srf      = (/ ctrl_out((/ 1, 10, 10, 10, 10 /),'fsw_ter'), &
                                                  ctrl_out((/ 1, 10, 10, 10, 10 /),'fsw_lic'), &
                                                  ctrl_out((/ 1, 10, 10, 10, 10 /),'fsw_oce'), &
                                                  ctrl_out((/ 1, 10, 10, 10, 10 /),'fsw_sic') /)

  type(ctrl_out),save,dimension(4) :: o_wbils_srf    = (/ ctrl_out((/ 1, 10, 10, 10, 10 /),'wbils_ter'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'wbils_lic'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'wbils_oce'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'wbils_sic') /)

  type(ctrl_out),save,dimension(4) :: o_wbilo_srf    = (/ ctrl_out((/ 1, 10, 10, 10, 10 /),'wbilo_ter'), &
                                                     ctrl_out((/ 1, 10, 10, 10, 10 /),'wbilo_lic'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'wbilo_oce'), &
                                                 ctrl_out((/ 1, 10, 10, 10, 10 /),'wbilo_sic') /)


  type(ctrl_out),save :: o_cdrm         = ctrl_out((/ 1, 10, 10, 10, 10 /),'cdrm')
  type(ctrl_out),save :: o_cdrh         = ctrl_out((/ 1, 10, 10, 1, 10 /),'cdrh')
  type(ctrl_out),save :: o_cldl         = ctrl_out((/ 1, 1, 10, 10, 10 /),'cldl')
  type(ctrl_out),save :: o_cldm         = ctrl_out((/ 1, 1, 10, 10, 10 /),'cldm')
  type(ctrl_out),save :: o_cldh         = ctrl_out((/ 1, 1, 10, 10, 10 /),'cldh')
  type(ctrl_out),save :: o_cldt         = ctrl_out((/ 1, 1, 2, 5, 10 /),'cldt')
  type(ctrl_out),save :: o_cldq         = ctrl_out((/ 1, 1, 10, 10, 10 /),'cldq')
  type(ctrl_out),save :: o_lwp          = ctrl_out((/ 1, 5, 10, 10, 10 /),'lwp')
  type(ctrl_out),save :: o_iwp          = ctrl_out((/ 1, 5, 10, 10, 10 /),'iwp')
  type(ctrl_out),save :: o_ue           = ctrl_out((/ 1, 10, 10, 10, 10 /),'ue')
  type(ctrl_out),save :: o_ve           = ctrl_out((/ 1, 10, 10, 10, 10 /),'ve')
  type(ctrl_out),save :: o_uq           = ctrl_out((/ 1, 10, 10, 10, 10 /),'uq')
  type(ctrl_out),save :: o_vq           = ctrl_out((/ 1, 10, 10, 10, 10 /),'vq')
 
  type(ctrl_out),save :: o_cape         = ctrl_out((/ 1, 10, 10, 10, 10 /),'cape')
  type(ctrl_out),save :: o_pbase        = ctrl_out((/ 1, 5, 10, 10, 10 /),'pbase')
  type(ctrl_out),save :: o_ptop         = ctrl_out((/ 1, 5, 10, 10, 10 /),'ptop')
  type(ctrl_out),save :: o_fbase        = ctrl_out((/ 1, 10, 10, 10, 10 /),'fbase')
  type(ctrl_out),save :: o_prw          = ctrl_out((/ 1, 1, 10, 10, 10 /),'prw')

  type(ctrl_out),save :: o_s_pblh       = ctrl_out((/ 1, 10, 10, 10, 1 /),'s_pblh')
  type(ctrl_out),save :: o_s_pblt       = ctrl_out((/ 1, 10, 10, 10, 1 /),'s_pblt')
  type(ctrl_out),save :: o_s_lcl        = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_lcl')
  type(ctrl_out),save :: o_s_capCL      = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_capCL')
  type(ctrl_out),save :: o_s_oliqCL     = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_oliqCL')
  type(ctrl_out),save :: o_s_cteiCL     = ctrl_out((/ 1, 10, 10, 10, 1 /),'s_cteiCL')
  type(ctrl_out),save :: o_s_therm      = ctrl_out((/ 1, 10, 10, 10, 1 /),'s_therm')
  type(ctrl_out),save :: o_s_trmb1      = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_trmb1')
  type(ctrl_out),save :: o_s_trmb2      = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_trmb2')
  type(ctrl_out),save :: o_s_trmb3      = ctrl_out((/ 1, 10, 10, 10, 10 /),'s_trmb3')

  type(ctrl_out),save :: o_slab_bils    = ctrl_out((/ 1, 1, 10, 10, 10 /),'slab_bils_oce')

  type(ctrl_out),save :: o_ale_bl       = ctrl_out((/ 1, 1, 1, 10, 10 /),'ale_bl')
  type(ctrl_out),save :: o_alp_bl       = ctrl_out((/ 1, 1, 1, 10, 10 /),'alp_bl')
  type(ctrl_out),save :: o_ale_wk       = ctrl_out((/ 1, 1, 1, 10, 10 /),'ale_wk')
  type(ctrl_out),save :: o_alp_wk       = ctrl_out((/ 1, 1, 1, 10, 10 /),'alp_wk')

  type(ctrl_out),save :: o_ale          = ctrl_out((/ 1, 1, 1, 10, 10 /),'ale')
  type(ctrl_out),save :: o_alp          = ctrl_out((/ 1, 1, 1, 10, 10 /),'alp')
  type(ctrl_out),save :: o_cin          = ctrl_out((/ 1, 1, 1, 10, 10 /),'cin')
  type(ctrl_out),save :: o_wape         = ctrl_out((/ 1, 1, 1, 10, 10 /),'wape')


! Champs interpolles sur des niveaux de pression ??? a faire correctement
                                              
  type(ctrl_out),save,dimension(6) :: o_uSTDlevs     = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'u850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'u700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'u500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'u200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'u50'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'u10') /) 
                                                     

  type(ctrl_out),save,dimension(6) :: o_vSTDlevs     = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'v850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'v700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'v500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'v200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'v50'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'v10') /)

  type(ctrl_out),save,dimension(6) :: o_wSTDlevs     = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'w850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'w700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'w500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'w200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'w50'), & 
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'w10') /)

  type(ctrl_out),save,dimension(6) :: o_tSTDlevs     = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'t850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'t700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'t500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'t200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'t50'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'t10') /)

  type(ctrl_out),save,dimension(6) :: o_qSTDlevs     = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'q850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'q700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'q500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'q200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'q50'), & 
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'q10') /)

  type(ctrl_out),save,dimension(6) :: o_zSTDlevs   = (/ ctrl_out((/ 1, 1, 3, 10, 10 /),'z850'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'z700'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'z500'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'z200'), &
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'z50'), & 
                                                     ctrl_out((/ 1, 1, 3, 10, 10 /),'z10') /)


  type(ctrl_out),save :: o_t_oce_sic    = ctrl_out((/ 1, 10, 10, 10, 10 /),'t_oce_sic')

  type(ctrl_out),save :: o_weakinv      = ctrl_out((/ 10, 1, 10, 10, 10 /),'weakinv')
  type(ctrl_out),save :: o_dthmin       = ctrl_out((/ 10, 1, 10, 10, 10 /),'dthmin')
  type(ctrl_out),save,dimension(4) :: o_u10_srf      = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'u10_ter'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'u10_lic'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'u10_oce'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'u10_sic') /)

  type(ctrl_out),save,dimension(4) :: o_v10_srf      = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'v10_ter'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'v10_lic'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'v10_oce'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'v10_sic') /)

  type(ctrl_out),save :: o_cldtau       = ctrl_out((/ 10, 5, 10, 10, 10 /),'cldtau')                     
  type(ctrl_out),save :: o_cldemi       = ctrl_out((/ 10, 5, 10, 10, 10 /),'cldemi')
  type(ctrl_out),save :: o_rh2m         = ctrl_out((/ 5, 5, 10, 10, 10 /),'rh2m')
  type(ctrl_out),save :: o_rh2m_min     = ctrl_out((/ 10, 5, 10, 10, 10 /),'rh2m_min')
  type(ctrl_out),save :: o_rh2m_max     = ctrl_out((/ 10, 5, 10, 10, 10 /),'rh2m_max')
  type(ctrl_out),save :: o_qsat2m       = ctrl_out((/ 10, 5, 10, 10, 10 /),'qsat2m')
  type(ctrl_out),save :: o_tpot         = ctrl_out((/ 10, 5, 10, 10, 10 /),'tpot')
  type(ctrl_out),save :: o_tpote        = ctrl_out((/ 10, 5, 10, 10, 10 /),'tpote')
  type(ctrl_out),save :: o_tke          = ctrl_out((/ 4, 10, 10, 10, 10 /),'tke ')
  type(ctrl_out),save :: o_tke_max      = ctrl_out((/ 4, 10, 10, 10, 10 /),'tke_max')

  type(ctrl_out),save,dimension(4) :: o_tke_srf      = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_ter'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_lic'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_oce'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_sic') /)

  type(ctrl_out),save,dimension(4) :: o_tke_max_srf  = (/ ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_max_ter'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_max_lic'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_max_oce'), &
                                                     ctrl_out((/ 10, 4, 10, 10, 10 /),'tke_max_sic') /)

  type(ctrl_out),save :: o_kz           = ctrl_out((/ 4, 10, 10, 10, 10 /),'kz')
  type(ctrl_out),save :: o_kz_max       = ctrl_out((/ 4, 10, 10, 10, 10 /),'kz_max')
  type(ctrl_out),save :: o_SWnetOR      = ctrl_out((/ 10, 10, 2, 10, 10 /),'SWnetOR')
  type(ctrl_out),save :: o_SWdownOR     = ctrl_out((/ 10, 10, 2, 10, 10 /),'SWdownOR')
  type(ctrl_out),save :: o_LWdownOR     = ctrl_out((/ 10, 10, 2, 10, 10 /),'LWdownOR')

  type(ctrl_out),save :: o_snowl        = ctrl_out((/ 10, 1, 10, 10, 10 /),'snowl')
  type(ctrl_out),save :: o_cape_max     = ctrl_out((/ 10, 1, 10, 10, 10 /),'cape_max')
  type(ctrl_out),save :: o_solldown     = ctrl_out((/ 10, 1, 10, 10, 10 /),'solldown')

  type(ctrl_out),save :: o_dtsvdfo      = ctrl_out((/ 10, 10, 10, 10, 10 /),'dtsvdfo')
  type(ctrl_out),save :: o_dtsvdft      = ctrl_out((/ 10, 10, 10, 10, 10 /),'dtsvdft')
  type(ctrl_out),save :: o_dtsvdfg      = ctrl_out((/ 10, 10, 10, 10, 10 /),'dtsvdfg')
  type(ctrl_out),save :: o_dtsvdfi      = ctrl_out((/ 10, 10, 10, 10, 10 /),'dtsvdfi')
  type(ctrl_out),save :: o_rugs         = ctrl_out((/ 10, 10, 10, 10, 10 /),'rugs')

  type(ctrl_out),save :: o_topswad      = ctrl_out((/ 2, 10, 10, 10, 10 /),'topswad')
  type(ctrl_out),save :: o_topswai      = ctrl_out((/ 2, 10, 10, 10, 10 /),'topswai')
  type(ctrl_out),save :: o_solswad      = ctrl_out((/ 2, 10, 10, 10, 10 /),'solswad')
  type(ctrl_out),save :: o_solswai      = ctrl_out((/ 2, 10, 10, 10, 10 /),'solswai')

  type(ctrl_out),save,dimension(10) :: o_tausumaero  = (/ ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_ASBCM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_ASPOMM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_ASSO4M'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_CSSO4M'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_SSSSM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_ASSSM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_CSSSM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_CIDUSTM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_AIBCM'), &
                                                     ctrl_out((/ 4, 4, 10, 10, 10 /),'OD550_AIPOMM') /)

  type(ctrl_out),save :: o_od550aer         = ctrl_out((/ 4, 4, 10, 10, 10 /),'od550aer')
  type(ctrl_out),save :: o_od865aer         = ctrl_out((/ 4, 4, 10, 10, 10 /),'od865aer')
  type(ctrl_out),save :: o_absvisaer        = ctrl_out((/ 4, 4, 10, 10, 10 /),'absvisaer')
  type(ctrl_out),save :: o_od550lt1aer      = ctrl_out((/ 4, 4, 10, 10, 10 /),'od550lt1aer')

  type(ctrl_out),save :: o_sconcso4         = ctrl_out((/ 4, 4, 10, 10, 10 /),'sconcso4')
  type(ctrl_out),save :: o_sconcoa          = ctrl_out((/ 4, 4, 10, 10, 10 /),'sconcoa')
  type(ctrl_out),save :: o_sconcbc          = ctrl_out((/ 4, 4, 10, 10, 10 /),'sconcbc')
  type(ctrl_out),save :: o_sconcss          = ctrl_out((/ 4, 4, 10, 10, 10 /),'sconcss')
  type(ctrl_out),save :: o_sconcdust        = ctrl_out((/ 4, 4, 10, 10, 10 /),'sconcdust')
  type(ctrl_out),save :: o_concso4          = ctrl_out((/ 4, 4, 10, 10, 10 /),'concso4')
  type(ctrl_out),save :: o_concoa           = ctrl_out((/ 4, 4, 10, 10, 10 /),'concoa')
  type(ctrl_out),save :: o_concbc           = ctrl_out((/ 4, 4, 10, 10, 10 /),'concbc')
  type(ctrl_out),save :: o_concss           = ctrl_out((/ 4, 4, 10, 10, 10 /),'concss')
  type(ctrl_out),save :: o_concdust         = ctrl_out((/ 4, 4, 10, 10, 10 /),'concdust')
  type(ctrl_out),save :: o_loadso4          = ctrl_out((/ 4, 4, 10, 10, 10 /),'loadso4')
  type(ctrl_out),save :: o_loadoa           = ctrl_out((/ 4, 4, 10, 10, 10 /),'loadoa')
  type(ctrl_out),save :: o_loadbc           = ctrl_out((/ 4, 4, 10, 10, 10 /),'loadbc')
  type(ctrl_out),save :: o_loadss           = ctrl_out((/ 4, 4, 10, 10, 10 /),'loadss')
  type(ctrl_out),save :: o_loaddust         = ctrl_out((/ 4, 4, 10, 10, 10 /),'loaddust')


  type(ctrl_out),save :: o_swtoaas_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoaas_nat')
  type(ctrl_out),save :: o_swsrfas_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfas_nat')
  type(ctrl_out),save :: o_swtoacs_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoacs_nat')
  type(ctrl_out),save :: o_swsrfcs_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfcs_nat')

  type(ctrl_out),save :: o_swtoaas_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoaas_ant')
  type(ctrl_out),save :: o_swsrfas_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfas_ant')
  type(ctrl_out),save :: o_swtoacs_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoacs_ant')
  type(ctrl_out),save :: o_swsrfcs_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfcs_ant')
  type(ctrl_out),save :: o_swtoacf_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoacf_nat')
  type(ctrl_out),save :: o_swsrfcf_nat      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfcf_nat')
  type(ctrl_out),save :: o_swtoacf_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoacf_ant')
  type(ctrl_out),save :: o_swsrfcf_ant      = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfcf_ant')
  type(ctrl_out),save :: o_swtoacf_zero     = ctrl_out((/ 4, 4, 10, 10, 10 /),'swtoacf_zero')
  type(ctrl_out),save :: o_swsrfcf_zero     = ctrl_out((/ 4, 4, 10, 10, 10 /),'swsrfcf_zero')

  type(ctrl_out),save :: o_cldncl          = ctrl_out((/ 4, 4, 10, 10, 10 /),'cldncl')
  type(ctrl_out),save :: o_reffclwtop      = ctrl_out((/ 4, 4, 10, 10, 10 /),'reffclwtop')
  type(ctrl_out),save :: o_cldnvi          = ctrl_out((/ 4, 4, 10, 10, 10 /),'cldnvi')
  type(ctrl_out),save :: o_lcc             = ctrl_out((/ 4, 4, 10, 10, 10 /),'lcc')


!!!!!!!!!!!!!!!!!!!!!! 3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  type(ctrl_out),save :: o_ec550aer     = ctrl_out((/ 4,  4, 10, 10, 1 /),'ec550aer')
  type(ctrl_out),save :: o_lwcon        = ctrl_out((/ 2, 5, 10, 10, 1 /),'lwcon')
  type(ctrl_out),save :: o_iwcon        = ctrl_out((/ 2, 5, 10, 10, 10 /),'iwcon')
  type(ctrl_out),save :: o_temp         = ctrl_out((/ 2, 3, 4, 10, 1 /),'temp')
  type(ctrl_out),save :: o_theta        = ctrl_out((/ 2, 3, 4, 10, 1 /),'theta')
  type(ctrl_out),save :: o_ovap         = ctrl_out((/ 2, 3, 4, 10, 1 /),'ovap')
  type(ctrl_out),save :: o_ovapinit         = ctrl_out((/ 2, 3, 10, 10, 1 /),'ovapinit')
  type(ctrl_out),save :: o_wvapp        = ctrl_out((/ 2, 10, 10, 10, 10 /),'wvapp')
  type(ctrl_out),save :: o_geop         = ctrl_out((/ 2, 3, 10, 10, 1 /),'geop')
  type(ctrl_out),save :: o_vitu         = ctrl_out((/ 2, 3, 4, 5, 1 /),'vitu')
  type(ctrl_out),save :: o_vitv         = ctrl_out((/ 2, 3, 4, 5, 1 /),'vitv')
  type(ctrl_out),save :: o_vitw         = ctrl_out((/ 2, 3, 10, 5, 1 /),'vitw')
  type(ctrl_out),save :: o_pres         = ctrl_out((/ 2, 3, 10, 10, 1 /),'pres')
  type(ctrl_out),save :: o_paprs        = ctrl_out((/ 2, 3, 10, 10, 1 /),'paprs')
  type(ctrl_out),save :: o_mass        = ctrl_out((/ 2, 3, 10, 10, 1 /),'mass')

  type(ctrl_out),save :: o_rneb         = ctrl_out((/ 2, 5, 10, 10, 1 /),'rneb')
  type(ctrl_out),save :: o_rnebcon      = ctrl_out((/ 2, 5, 10, 10, 1 /),'rnebcon')
  type(ctrl_out),save :: o_rhum         = ctrl_out((/ 2, 5, 10, 10, 10 /),'rhum')
  type(ctrl_out),save :: o_ozone        = ctrl_out((/ 2, 10, 10, 10, 10 /),'ozone')
  type(ctrl_out),save :: o_ozone_light  = ctrl_out((/ 2, 10, 10, 10, 10 /),'ozone_daylight')
  type(ctrl_out),save :: o_upwd         = ctrl_out((/ 2, 10, 10, 10, 10 /),'upwd')
  type(ctrl_out),save :: o_dtphy        = ctrl_out((/ 2, 10, 10, 10, 10 /),'dtphy')
  type(ctrl_out),save :: o_dqphy        = ctrl_out((/ 2, 10, 10, 10, 10 /),'dqphy')
  type(ctrl_out),save :: o_pr_con_l     = ctrl_out((/ 2, 10, 10, 10, 10 /),'pr_con_l')
  type(ctrl_out),save :: o_pr_con_i     = ctrl_out((/ 2, 10, 10, 10, 10 /),'pr_con_i')
  type(ctrl_out),save :: o_pr_lsc_l     = ctrl_out((/ 2, 10, 10, 10, 10 /),'pr_lsc_l')
  type(ctrl_out),save :: o_pr_lsc_i     = ctrl_out((/ 2, 10, 10, 10, 10 /),'pr_lsc_i')
  type(ctrl_out),save :: o_re           = ctrl_out((/ 5, 10, 10, 10, 10 /),'re')
  type(ctrl_out),save :: o_fl           = ctrl_out((/ 5, 10, 10, 10, 10 /),'fl')
  type(ctrl_out),save :: o_scdnc        =ctrl_out((/ 4,  4, 10, 10, 1 /),'scdnc')
  type(ctrl_out),save :: o_reffclws     =ctrl_out((/ 4,  4, 10, 10, 1 /),'reffclws')
  type(ctrl_out),save :: o_reffclwc     =ctrl_out((/ 4,  4, 10, 10, 1 /),'reffclwc')
  type(ctrl_out),save :: o_lcc3d        =ctrl_out((/ 4,  4, 10, 10, 1 /),'lcc3d')
  type(ctrl_out),save :: o_lcc3dcon     =ctrl_out((/ 4,  4, 10, 10, 1 /),'lcc3dcon')
  type(ctrl_out),save :: o_lcc3dstra    =ctrl_out((/ 4,  4, 10, 10, 1 /),'lcc3dstra')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(ctrl_out),save,dimension(4) :: o_albe_srf     = (/ ctrl_out((/ 3, 4, 10, 1, 10 /),'albe_ter'), &
                                                     ctrl_out((/ 3, 4, 10, 1, 10 /),'albe_lic'), &
                                                     ctrl_out((/ 3, 4, 10, 1, 10 /),'albe_oce'), &
                                                     ctrl_out((/ 3, 4, 10, 1, 10 /),'albe_sic') /) 

  type(ctrl_out),save,dimension(4) :: o_ages_srf     = (/ ctrl_out((/ 10, 10, 10, 10, 10 /),'ages_ter'), &
                                                     ctrl_out((/ 3, 10, 10, 10, 10 /),'ages_lic'), &
                                                     ctrl_out((/ 10, 10, 10, 10, 10 /),'ages_oce'), &
                                                     ctrl_out((/ 3, 10, 10, 10, 10 /),'ages_sic') /)

  type(ctrl_out),save,dimension(4) :: o_rugs_srf     = (/ ctrl_out((/ 3, 4, 10, 10, 10 /),'rugs_ter'), &
                                                     ctrl_out((/ 3, 4, 10, 10, 10 /),'rugs_lic'), &
                                                     ctrl_out((/ 3, 4, 10, 10, 10 /),'rugs_oce'), &
                                                     ctrl_out((/ 3, 4, 10, 10, 10 /),'rugs_sic') /)

  type(ctrl_out),save :: o_alb1         = ctrl_out((/ 3, 10, 10, 10, 10 /),'albs')
  type(ctrl_out),save :: o_alb2       = ctrl_out((/ 3, 10, 10, 10, 10 /),'albslw')

  type(ctrl_out),save :: o_clwcon       = ctrl_out((/ 4, 10, 10, 10, 10 /),'clwcon')
  type(ctrl_out),save :: o_Ma           = ctrl_out((/ 4, 10, 10, 10, 10 /),'Ma')
  type(ctrl_out),save :: o_dnwd         = ctrl_out((/ 4, 10, 10, 10, 10 /),'dnwd')
  type(ctrl_out),save :: o_dnwd0        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dnwd0')
  type(ctrl_out),save :: o_mc           = ctrl_out((/ 4, 10, 10, 10, 10 /),'mc')
  type(ctrl_out),save :: o_ftime_con    = ctrl_out((/ 4, 10, 10, 10, 10 /),'ftime_con')
  type(ctrl_out),save :: o_dtdyn        = ctrl_out((/ 4, 10, 10, 10, 1 /),'dtdyn')
  type(ctrl_out),save :: o_dqdyn        = ctrl_out((/ 4, 10, 10, 10, 1 /),'dqdyn')
  type(ctrl_out),save :: o_dudyn        = ctrl_out((/ 4, 10, 10, 10, 1 /),'dudyn')  !AXC
  type(ctrl_out),save :: o_dvdyn        = ctrl_out((/ 4, 10, 10, 10, 1 /),'dvdyn')  !AXC
  type(ctrl_out),save :: o_dtcon        = ctrl_out((/ 4, 5, 10, 10, 10 /),'dtcon')
  type(ctrl_out),save :: o_ducon        = ctrl_out((/ 4, 10, 10, 10, 10 /),'ducon')
  type(ctrl_out),save :: o_dqcon        = ctrl_out((/ 4, 5, 10, 10, 10 /),'dqcon')
  type(ctrl_out),save :: o_dtwak        = ctrl_out((/ 4, 5, 10, 10, 10 /),'dtwak')
  type(ctrl_out),save :: o_dqwak        = ctrl_out((/ 4, 5, 10, 10, 10 /),'dqwak')
  type(ctrl_out),save :: o_wake_h       = ctrl_out((/ 4, 5, 10, 10, 10 /),'wake_h')
  type(ctrl_out),save :: o_wake_s       = ctrl_out((/ 4, 5, 10, 10, 10 /),'wake_s')
  type(ctrl_out),save :: o_wake_deltat  = ctrl_out((/ 4, 5, 10, 10, 10 /),'wake_deltat')
  type(ctrl_out),save :: o_wake_deltaq  = ctrl_out((/ 4, 5, 10, 10, 10 /),'wake_deltaq')
  type(ctrl_out),save :: o_wake_omg     = ctrl_out((/ 4, 5, 10, 10, 10 /),'wake_omg')
  type(ctrl_out),save :: o_Vprecip      = ctrl_out((/ 10, 10, 10, 10, 10 /),'Vprecip')
  type(ctrl_out),save :: o_ftd          = ctrl_out((/ 4, 5, 10, 10, 10 /),'ftd')
  type(ctrl_out),save :: o_fqd          = ctrl_out((/ 4, 5, 10, 10, 10 /),'fqd')
  type(ctrl_out),save :: o_dtlsc        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtlsc')
  type(ctrl_out),save :: o_dtlschr      = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtlschr')
  type(ctrl_out),save :: o_dqlsc        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dqlsc')
  type(ctrl_out),save :: o_dtvdf        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtvdf')
  type(ctrl_out),save :: o_dqvdf        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dqvdf')
  type(ctrl_out),save :: o_dteva        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dteva')
  type(ctrl_out),save :: o_dqeva        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dqeva')
  type(ctrl_out),save :: o_ptconv       = ctrl_out((/ 4, 10, 10, 10, 10 /),'ptconv')
  type(ctrl_out),save :: o_ratqs        = ctrl_out((/ 4, 10, 10, 10, 10 /),'ratqs')
  type(ctrl_out),save :: o_dtthe        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtthe')
  type(ctrl_out),save :: o_f_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'f_th')
  type(ctrl_out),save :: o_e_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'e_th')
  type(ctrl_out),save :: o_w_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'w_th')
  type(ctrl_out),save :: o_lambda_th    = ctrl_out((/ 10, 10, 10, 10, 10 /),'lambda_th')
  type(ctrl_out),save :: o_ftime_th     = ctrl_out((/ 10, 10, 10, 10, 10 /),'ftime_th')
  type(ctrl_out),save :: o_q_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'q_th')
  type(ctrl_out),save :: o_a_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'a_th')
  type(ctrl_out),save :: o_d_th         = ctrl_out((/ 4, 10, 10, 10, 10 /),'d_th')
  type(ctrl_out),save :: o_f0_th        = ctrl_out((/ 4, 10, 10, 10, 10 /),'f0_th')
  type(ctrl_out),save :: o_zmax_th      = ctrl_out((/ 4, 10, 10, 10, 10 /),'zmax_th')
  type(ctrl_out),save :: o_dqthe        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dqthe')
  type(ctrl_out),save :: o_dtajs        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtajs')
  type(ctrl_out),save :: o_dqajs        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dqajs')
  type(ctrl_out),save :: o_dtswr        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtswr')
  type(ctrl_out),save :: o_dtsw0        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtsw0')
  type(ctrl_out),save :: o_dtlwr        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtlwr')
  type(ctrl_out),save :: o_dtlw0        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtlw0')
  type(ctrl_out),save :: o_dtec         = ctrl_out((/ 4, 10, 10, 10, 10 /),'dtec')
  type(ctrl_out),save :: o_duvdf        = ctrl_out((/ 4, 10, 10, 10, 10 /),'duvdf')
  type(ctrl_out),save :: o_dvvdf        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dvvdf')
  type(ctrl_out),save :: o_duoro        = ctrl_out((/ 4, 10, 10, 10, 10 /),'duoro')
  type(ctrl_out),save :: o_dvoro        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dvoro')
  type(ctrl_out),save :: o_dulif        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dulif')
  type(ctrl_out),save :: o_dvlif        = ctrl_out((/ 4, 10, 10, 10, 10 /),'dvlif')
  type(ctrl_out),save,allocatable :: o_trac(:)
    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Ouverture des fichier et definition des variable de sortie !!!!!!!!
!! histbeg, histvert et histdef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  SUBROUTINE phys_output_open(jjmp1,nlevSTD,clevSTD,nbteta, &
       ctetaSTD,dtime, ok_veget, &
       type_ocean, iflag_pbl,ok_mensuel,ok_journe, &
       ok_hf,ok_instan,ok_LES,ok_ade,ok_aie, read_climoz, &
       new_aod, aerosol_couple)   


  USE iophy 
  USE dimphy
  USE infotrac
  USE ioipsl
  USE mod_phys_lmdz_para
  USE aero_mod, only : naero_spc,name_aero

  IMPLICIT NONE
  include "dimensions.h"
  include "temps.h"
  include "indicesol.h"
  include "clesphys.h"
  include "thermcell.h"
  include "comvert.h"

  integer                               :: jjmp1
  integer                               :: nbteta, nlevSTD, radpas
  logical                               :: ok_mensuel, ok_journe, ok_hf, ok_instan
  logical                               :: ok_LES,ok_ade,ok_aie
  logical                               :: new_aod, aerosol_couple
  integer, intent(in)::  read_climoz ! read ozone climatology
  !     Allowed values are 0, 1 and 2
  !     0: do not read an ozone climatology
  !     1: read a single ozone climatology that will be used day and night
  !     2: read two ozone climatologies, the average day and night
  !     climatology and the daylight climatology

  real                                  :: dtime
  integer                               :: idayref
  real                                  :: zjulian
  real, dimension(klev)                 :: Ahyb, Bhyb, Alt
  character(len=4), dimension(nlevSTD)  :: clevSTD
  integer                               :: nsrf, k, iq, iiq, iff, i, j, ilev
  integer                               :: naero
  logical                               :: ok_veget
  integer                               :: iflag_pbl
  CHARACTER(len=4)                      :: bb2
  CHARACTER(len=2)                      :: bb3
  character(len=6)                      :: type_ocean
  CHARACTER(len=3)                      :: ctetaSTD(nbteta)
  real, dimension(nfiles)               :: ecrit_files
  CHARACTER(len=20), dimension(nfiles)  :: phys_out_filenames
  INTEGER, dimension(iim*jjmp1)         ::  ndex2d
  INTEGER, dimension(iim*jjmp1*klev)    :: ndex3d
  integer                               :: imin_ins, imax_ins
  integer                               :: jmin_ins, jmax_ins
  integer, dimension(nfiles)            :: phys_out_levmin, phys_out_levmax
  integer, dimension(nfiles)            :: phys_out_filelevels
  CHARACTER(len=20), dimension(nfiles)  :: type_ecri_files, phys_out_filetypes
  character(len=20), dimension(nfiles)  :: chtimestep   = (/ 'DefFreq', 'DefFreq','DefFreq', 'DefFreq', 'DefFreq' /)
  logical, dimension(nfiles)            :: phys_out_filekeys

!!!!!!!!!! stockage dans une region limitee pour chaque fichier !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 entre [phys_out_lonmin,phys_out_lonmax] et [phys_out_latmin,phys_out_latmax]

  logical, dimension(nfiles), save  :: phys_out_regfkey       = (/ .false., .false., .false., .false., .false. /)
  real, dimension(nfiles), save     :: phys_out_lonmin        = (/ -180., -180., -180., -180., -180. /)
  real, dimension(nfiles), save     :: phys_out_lonmax        = (/ 180., 180., 180., 180., 180. /)
  real, dimension(nfiles), save     :: phys_out_latmin        = (/ -90., -90., -90., -90., -90. /)
  real, dimension(nfiles), save     :: phys_out_latmax        = (/ 90., 90., 90., 90., 90. /)

!IM definition dynamique flag o_trac pour sortie traceurs
  INTEGER :: nq
  CHARACTER(len=8) :: solsym(nqtot)

   print*,'Debut phys_output_mod.F90'
! Initialisations (Valeurs par defaut

   if (.not. allocated(o_trac)) ALLOCATE(o_trac(nqtot))

   levmax = (/ klev, klev, klev, klev, klev /)

   phys_out_filenames(1) = 'histmth'
   phys_out_filenames(2) = 'histday'
   phys_out_filenames(3) = 'histhf'
   phys_out_filenames(4) = 'histins'
   phys_out_filenames(5) = 'histLES'

   type_ecri(1) = 'ave(X)'
   type_ecri(2) = 'ave(X)'
   type_ecri(3) = 'ave(X)'
   type_ecri(4) = 'inst(X)'
   type_ecri(5) = 'ave(X)'

   clef_files(1) = ok_mensuel
   clef_files(2) = ok_journe
   clef_files(3) = ok_hf
   clef_files(4) = ok_instan
   clef_files(5) = ok_LES

   lev_files(1) = lev_histmth
   lev_files(2) = lev_histday
   lev_files(3) = lev_histhf
   lev_files(4) = lev_histins
   lev_files(5) = lev_histLES


   ecrit_files(1) = ecrit_mth
   ecrit_files(2) = ecrit_day
   ecrit_files(3) = ecrit_hf
   ecrit_files(4) = ecrit_ins
   ecrit_files(5) = ecrit_LES
 
!! Lectures des parametres de sorties dans physiq.def

   call getin('phys_out_regfkey',phys_out_regfkey)
   call getin('phys_out_lonmin',phys_out_lonmin)
   call getin('phys_out_lonmax',phys_out_lonmax)
   call getin('phys_out_latmin',phys_out_latmin)
   call getin('phys_out_latmax',phys_out_latmax)
     phys_out_levmin(:)=levmin(:)
   call getin('phys_out_levmin',levmin)
     phys_out_levmax(:)=levmax(:)
   call getin('phys_out_levmax',levmax)
   call getin('phys_out_filenames',phys_out_filenames)
     phys_out_filekeys(:)=clef_files(:)
   call getin('phys_out_filekeys',clef_files)
     phys_out_filelevels(:)=lev_files(:)
   call getin('phys_out_filelevels',lev_files)
   call getin('phys_out_filetimesteps',chtimestep)
     phys_out_filetypes(:)=type_ecri(:)
   call getin('phys_out_filetypes',type_ecri)

   type_ecri_files(:)=type_ecri(:)

   print*,'phys_out_lonmin=',phys_out_lonmin
   print*,'phys_out_lonmax=',phys_out_lonmax
   print*,'phys_out_latmin=',phys_out_latmin
   print*,'phys_out_latmax=',phys_out_latmax
   print*,'phys_out_filenames=',phys_out_filenames
   print*,'phys_out_filetypes=',type_ecri
   print*,'phys_out_filekeys=',clef_files
   print*,'phys_out_filelevels=',lev_files

!!!!!!!!!!!!!!!!!!!!!!! Boucle sur les fichiers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel de histbeg et histvert pour creer le fichier et les niveaux verticaux !!
! Appel des histbeg pour definir les variables (nom, moy ou inst, freq de sortie ..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 zdtime = dtime         ! Frequence ou l on moyenne

! Calcul des Ahyb, Bhyb et Alt
         do k=1,klev
          Ahyb(k)=(ap(k)+ap(k+1))/2.
          Bhyb(k)=(bp(k)+bp(k+1))/2.
          Alt(k)=log(preff/presnivs(k))*8.
         enddo
!          if(prt_level.ge.1) then
           print*,'Ap Hybrid = ',Ahyb(1:klev)
           print*,'Bp Hybrid = ',Bhyb(1:klev)
           print*,'Alt approx des couches pour une haut d echelle de 8km = ',Alt(1:klev)
!          endif
 DO iff=1,nfiles

    IF (clef_files(iff)) THEN

      if ( chtimestep(iff).eq.'DefFreq' ) then
! Par defaut ecrit_files = (ecrit_mensuel ecrit_jour ecrit_hf ...)*86400.
        ecrit_files(iff)=ecrit_files(iff)*86400.
      else
        call convers_timesteps(chtimestep(iff),ecrit_files(iff)) 
      endif
       print*,'ecrit_files(',iff,')= ',ecrit_files(iff)

      zoutm(iff) = ecrit_files(iff) ! Frequence ou l on ecrit en seconde

      idayref = day_ref
      CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)

!!!!!!!!!!!!!!!!! Traitement dans le cas ou l'on veut stocker sur un domaine limite !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (phys_out_regfkey(iff)) then

        imin_ins=1
        imax_ins=iim
        jmin_ins=1
        jmax_ins=jjmp1

! correction abderr        
        do i=1,iim
           print*,'io_lon(i)=',io_lon(i)
           if (io_lon(i).le.phys_out_lonmin(iff)) imin_ins=i
           if (io_lon(i).le.phys_out_lonmax(iff)) imax_ins=i+1
        enddo

        do j=1,jjmp1
            print*,'io_lat(j)=',io_lat(j)
            if (io_lat(j).ge.phys_out_latmin(iff)) jmax_ins=j+1
            if (io_lat(j).ge.phys_out_latmax(iff)) jmin_ins=j
        enddo

        print*,'On stoke le fichier histoire numero ',iff,' sur ', &
         imin_ins,imax_ins,jmin_ins,jmax_ins
         print*,'longitudes : ', &
         io_lon(imin_ins),io_lon(imax_ins), &
         'latitudes : ', &
         io_lat(jmax_ins),io_lat(jmin_ins)

 CALL histbeg(phys_out_filenames(iff),iim,io_lon,jjmp1,io_lat, &
              imin_ins,imax_ins-imin_ins+1, &
              jmin_ins,jmax_ins-jmin_ins+1, &
              itau_phy,zjulian,dtime,nhorim(iff),nid_files(iff))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else
 CALL histbeg_phy(phys_out_filenames(iff),itau_phy,zjulian,dtime,nhorim(iff),nid_files(iff))
       endif
 
      CALL histvert(nid_files(iff), "presnivs", "Vertical levels", "Pa", &
           levmax(iff) - levmin(iff) + 1, &
           presnivs(levmin(iff):levmax(iff)), nvertm(iff),"down")

!!!!!!!!!!!!! Traitement des champs 3D pour histhf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! A Revoir plus tard !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IF (iff.eq.3.and.lev_files(iff).ge.4) THEN
!          CALL histbeg_phy("histhf3d",itau_phy, &
!     &                     zjulian, dtime, &
!     &                     nhorim, nid_hf3d)

!         CALL histvert(nid_hf3d, "presnivs", &
!     &                 "Vertical levels", "mb", &
!     &                 klev, presnivs/100., nvertm)
!          ENDIF
!
!!!! Composentes de la coordonnee sigma-hybride 
   CALL histvert(nid_files(iff), "Ahyb","Ahyb comp of Hyb Cord ", "Pa", &
                 levmax(iff) - levmin(iff) + 1,Ahyb,nvertap(iff))

   CALL histvert(nid_files(iff), "Bhyb","Bhyb comp of Hyb Cord", " ", &
                 levmax(iff) - levmin(iff) + 1,Bhyb,nvertbp(iff))

   CALL histvert(nid_files(iff), "Alt","Height approx for scale heigh of 8km at levels", "Km", &
                 levmax(iff) - levmin(iff) + 1,Alt,nvertAlt(iff))

!   CALL histvert(nid_files(iff), "preff","Reference pressure", "Pa", &
!                 1,preff,nvertp0(iff))
!!! Champs 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL histdef2d(iff,o_phis%flag,o_phis%name,"Surface geop.height", "m2/s2") 
   type_ecri(1) = 'once'
   type_ecri(2) = 'once'
   type_ecri(3) = 'once'
   type_ecri(4) = 'once'
   type_ecri(5) = 'once'
 CALL histdef2d(iff,o_aire%flag,o_aire%name,"Grid area", "-")
 CALL histdef2d(iff,o_contfracATM%flag,o_contfracATM%name,"% sfce ter+lic", "-")
   type_ecri(:) = type_ecri_files(:)

!!! Champs 2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CALL histdef2d(iff,o_contfracOR%flag,o_contfracOR%name,"% sfce terre OR", "-" )
 CALL histdef2d(iff,o_aireTER%flag,o_aireTER%name,"Grid area CONT", "-" )
 CALL histdef2d(iff,o_flat%flag,o_flat%name, "Latent heat flux", "W/m2")
 CALL histdef2d(iff,o_slp%flag,o_slp%name, "Sea Level Pressure", "Pa" )
 CALL histdef2d(iff,o_tsol%flag,o_tsol%name, "Surface Temperature", "K")
 CALL histdef2d(iff,o_t2m%flag,o_t2m%name, "Temperature 2m", "K" )
   type_ecri(1) = 't_min(X)'
   type_ecri(2) = 't_min(X)'
   type_ecri(3) = 't_min(X)'
   type_ecri(4) = 't_min(X)'
   type_ecri(5) = 't_min(X)' 
 CALL histdef2d(iff,o_t2m_min%flag,o_t2m_min%name, "Temp 2m min", "K" )
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)' 
 CALL histdef2d(iff,o_t2m_max%flag,o_t2m_max%name, "Temp 2m max", "K" )
   type_ecri(:) = type_ecri_files(:)
 CALL histdef2d(iff,o_wind10m%flag,o_wind10m%name, "10-m wind speed", "m/s")
 CALL histdef2d(iff,o_wind10max%flag,o_wind10max%name, "10m wind speed max", "m/s")
 CALL histdef2d(iff,o_sicf%flag,o_sicf%name, "Sea-ice fraction", "-" )
 CALL histdef2d(iff,o_q2m%flag,o_q2m%name, "Specific humidity 2m", "kg/kg")
 CALL histdef2d(iff,o_u10m%flag,o_u10m%name, "Vent zonal 10m", "m/s" )
 CALL histdef2d(iff,o_v10m%flag,o_v10m%name, "Vent meridien 10m", "m/s")
 CALL histdef2d(iff,o_psol%flag,o_psol%name, "Surface Pressure", "Pa" ) 
 CALL histdef2d(iff,o_qsurf%flag,o_qsurf%name, "Surface Air humidity", "kg/kg")

  if (.not. ok_veget) then
 CALL histdef2d(iff,o_qsol%flag,o_qsol%name, "Soil watter content", "mm" )
  endif

 CALL histdef2d(iff,o_ndayrain%flag,o_ndayrain%name, "Number of dayrain(liq+sol)", "-")
 CALL histdef2d(iff,o_precip%flag,o_precip%name, "Precip Totale liq+sol", "kg/(s*m2)" )
 CALL histdef2d(iff,o_plul%flag,o_plul%name, "Large-scale Precip.", "kg/(s*m2)") 
 CALL histdef2d(iff,o_pluc%flag,o_pluc%name, "Convective Precip.", "kg/(s*m2)")
 CALL histdef2d(iff,o_snow%flag,o_snow%name, "Snow fall", "kg/(s*m2)" )
 CALL histdef2d(iff,o_msnow%flag,o_msnow%name, "Surface snow amount", "kg/m2" )
 CALL histdef2d(iff,o_fsnow%flag,o_fsnow%name, "Surface snow area fraction", "-" )
 CALL histdef2d(iff,o_evap%flag,o_evap%name, "Evaporat", "kg/(s*m2)" ) 
 CALL histdef2d(iff,o_tops%flag,o_tops%name, "Solar rad. at TOA", "W/m2") 
 CALL histdef2d(iff,o_tops0%flag,o_tops0%name, "CS Solar rad. at TOA", "W/m2")
 CALL histdef2d(iff,o_topl%flag,o_topl%name, "IR rad. at TOA", "W/m2" )
 CALL histdef2d(iff,o_topl0%flag,o_topl0%name, "IR rad. at TOA", "W/m2")
 CALL histdef2d(iff,o_SWupTOA%flag,o_SWupTOA%name, "SWup at TOA", "W/m2") 
 CALL histdef2d(iff,o_SWupTOAclr%flag,o_SWupTOAclr%name, "SWup clear sky at TOA", "W/m2")
 CALL histdef2d(iff,o_SWdnTOA%flag,o_SWdnTOA%name, "SWdn at TOA", "W/m2" )
 CALL histdef2d(iff,o_SWdnTOAclr%flag,o_SWdnTOAclr%name, "SWdn clear sky at TOA", "W/m2") 
 CALL histdef2d(iff,o_nettop%flag,o_nettop%name, "Net dn radiatif flux at TOA", "W/m2")
 CALL histdef2d(iff,o_SWup200%flag,o_SWup200%name, "SWup at 200mb", "W/m2" ) 
 CALL histdef2d(iff,o_SWup200clr%flag,o_SWup200clr%name, "SWup clear sky at 200mb", "W/m2")
 CALL histdef2d(iff,o_SWdn200%flag,o_SWdn200%name, "SWdn at 200mb", "W/m2" )
 CALL histdef2d(iff,o_SWdn200clr%flag,o_SWdn200clr%name, "SWdn clear sky at 200mb", "W/m2")
 CALL histdef2d(iff,o_LWup200%flag,o_LWup200%name, "LWup at 200mb", "W/m2") 
 CALL histdef2d(iff,o_LWup200clr%flag,o_LWup200clr%name, "LWup clear sky at 200mb", "W/m2")
 CALL histdef2d(iff,o_LWdn200%flag,o_LWdn200%name, "LWdn at 200mb", "W/m2") 
 CALL histdef2d(iff,o_LWdn200clr%flag,o_LWdn200clr%name, "LWdn clear sky at 200mb", "W/m2")
 CALL histdef2d(iff,o_sols%flag,o_sols%name, "Solar rad. at surf.", "W/m2")
 CALL histdef2d(iff,o_sols0%flag,o_sols0%name, "Solar rad. at surf.", "W/m2")
 CALL histdef2d(iff,o_soll%flag,o_soll%name, "IR rad. at surface", "W/m2")  
 CALL histdef2d(iff,o_radsol%flag,o_radsol%name, "Rayonnement au sol", "W/m2")
 CALL histdef2d(iff,o_soll0%flag,o_soll0%name, "IR rad. at surface", "W/m2") 
 CALL histdef2d(iff,o_SWupSFC%flag,o_SWupSFC%name, "SWup at surface", "W/m2")
 CALL histdef2d(iff,o_SWupSFCclr%flag,o_SWupSFCclr%name, "SWup clear sky at surface", "W/m2")
 CALL histdef2d(iff,o_SWdnSFC%flag,o_SWdnSFC%name, "SWdn at surface", "W/m2")
 CALL histdef2d(iff,o_SWdnSFCclr%flag,o_SWdnSFCclr%name, "SWdn clear sky at surface", "W/m2")
 CALL histdef2d(iff,o_LWupSFC%flag,o_LWupSFC%name, "Upwd. IR rad. at surface", "W/m2")
 CALL histdef2d(iff,o_LWdnSFC%flag,o_LWdnSFC%name, "Down. IR rad. at surface", "W/m2")
 CALL histdef2d(iff,o_LWupSFCclr%flag,o_LWupSFCclr%name, "CS Upwd. IR rad. at surface", "W/m2")
 CALL histdef2d(iff,o_LWdnSFCclr%flag,o_LWdnSFCclr%name, "Down. CS IR rad. at surface", "W/m2")
 CALL histdef2d(iff,o_bils%flag,o_bils%name, "Surf. total heat flux", "W/m2")
 CALL histdef2d(iff,o_sens%flag,o_sens%name, "Sensible heat flux", "W/m2")
 CALL histdef2d(iff,o_fder%flag,o_fder%name, "Heat flux derivation", "W/m2")
 CALL histdef2d(iff,o_ffonte%flag,o_ffonte%name, "Thermal flux for snow melting", "W/m2")
 CALL histdef2d(iff,o_fqcalving%flag,o_fqcalving%name, "Ice Calving", "kg/m2/s") 
 CALL histdef2d(iff,o_fqfonte%flag,o_fqfonte%name, "Land ice melt", "kg/m2/s") 

 CALL histdef2d(iff,o_taux%flag,o_taux%name, "Zonal wind stress","Pa")
 CALL histdef2d(iff,o_tauy%flag,o_tauy%name, "Meridional wind stress","Pa")

     DO nsrf = 1, nbsrf
 CALL histdef2d(iff,o_pourc_srf(nsrf)%flag,o_pourc_srf(nsrf)%name,"% "//clnsurf(nsrf),"%")
 CALL histdef2d(iff,o_fract_srf(nsrf)%flag,o_fract_srf(nsrf)%name,"Fraction "//clnsurf(nsrf),"1")
 CALL histdef2d(iff,o_taux_srf(nsrf)%flag,o_taux_srf(nsrf)%name,"Zonal wind stress"//clnsurf(nsrf),"Pa")
 CALL histdef2d(iff,o_tauy_srf(nsrf)%flag,o_tauy_srf(nsrf)%name,"Meridional wind stress "//clnsurf(nsrf),"Pa")
 CALL histdef2d(iff,o_tsol_srf(nsrf)%flag,o_tsol_srf(nsrf)%name,"Temperature "//clnsurf(nsrf),"K")
 CALL histdef2d(iff,o_u10m_srf(nsrf)%flag,o_u10m_srf(nsrf)%name,"Vent Zonal 10m "//clnsurf(nsrf),"m/s")
 CALL histdef2d(iff,o_evap_srf(nsrf)%flag,o_evap_srf(nsrf)%name,"evaporation at surface "//clnsurf(nsrf),"kg/(s*m2)")
 CALL histdef2d(iff,o_v10m_srf(nsrf)%flag,o_v10m_srf(nsrf)%name,"Vent meredien 10m "//clnsurf(nsrf),"m/s")
 CALL histdef2d(iff,o_t2m_srf(nsrf)%flag,o_t2m_srf(nsrf)%name,"Temp 2m "//clnsurf(nsrf),"K")
 CALL histdef2d(iff,o_sens_srf(nsrf)%flag,o_sens_srf(nsrf)%name,"Sensible heat flux "//clnsurf(nsrf),"W/m2")
 CALL histdef2d(iff,o_lat_srf(nsrf)%flag,o_lat_srf(nsrf)%name,"Latent heat flux "//clnsurf(nsrf),"W/m2")
 CALL histdef2d(iff,o_flw_srf(nsrf)%flag,o_flw_srf(nsrf)%name,"LW "//clnsurf(nsrf),"W/m2")
 CALL histdef2d(iff,o_fsw_srf(nsrf)%flag,o_fsw_srf(nsrf)%name,"SW "//clnsurf(nsrf),"W/m2")
 CALL histdef2d(iff,o_wbils_srf(nsrf)%flag,o_wbils_srf(nsrf)%name,"Bilan sol "//clnsurf(nsrf),"W/m2" )
 CALL histdef2d(iff,o_wbilo_srf(nsrf)%flag,o_wbilo_srf(nsrf)%name,"Bilan eau "//clnsurf(nsrf),"kg/(m2*s)")
  if (iflag_pbl>1 .and. lev_files(iff).gt.10 ) then
 CALL histdef2d(iff,o_tke_srf(nsrf)%flag,o_tke_srf(nsrf)%name,"Max Turb. Kinetic Energy "//clnsurf(nsrf),"-")
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)'
 CALL histdef2d(iff,o_tke_max_srf(nsrf)%flag,o_tke_max_srf(nsrf)%name,"Max Turb. Kinetic Energy "//clnsurf(nsrf),"-")
   type_ecri(:) = type_ecri_files(:)
  endif
 CALL histdef2d(iff,o_albe_srf(nsrf)%flag,o_albe_srf(nsrf)%name,"Albedo VIS surf. "//clnsurf(nsrf),"-")
 CALL histdef2d(iff,o_rugs_srf(nsrf)%flag,o_rugs_srf(nsrf)%name,"Surface roughness "//clnsurf(nsrf),"m")
 CALL histdef2d(iff,o_ages_srf(nsrf)%flag,o_ages_srf(nsrf)%name,"Snow age", "day")
END DO

IF (new_aod .AND. (.NOT. aerosol_couple)) THEN
 IF (ok_ade.OR.ok_aie) THEN

  CALL histdef2d(iff,o_od550aer%flag,o_od550aer%name, "Total aerosol optical depth at 550nm", "-")
  CALL histdef2d(iff,o_od865aer%flag,o_od865aer%name, "Total aerosol optical depth at 870nm", "-")
  CALL histdef2d(iff,o_absvisaer%flag,o_absvisaer%name, "Absorption aerosol visible optical depth", "-")
  CALL histdef2d(iff,o_od550lt1aer%flag,o_od550lt1aer%name, "Fine mode optical depth", "-")
  CALL histdef2d(iff,o_sconcso4%flag,o_sconcso4%name,"Surface Concentration of Sulfate ","kg/m3")
  CALL histdef2d(iff,o_sconcoa%flag,o_sconcoa%name,"Surface Concentration of Organic Aerosol ","kg/m3")
  CALL histdef2d(iff,o_sconcbc%flag,o_sconcbc%name,"Surface Concentration of Black Carbon ","kg/m3")
  CALL histdef2d(iff,o_sconcss%flag,o_sconcss%name,"Surface Concentration of Sea Salt ","kg/m3")
  CALL histdef2d(iff,o_sconcdust%flag,o_sconcdust%name,"Surface Concentration of Dust ","kg/m3")
  CALL histdef3d(iff,o_concso4%flag,o_concso4%name,"Concentration of Sulfate ","kg/m3")
  CALL histdef3d(iff,o_concoa%flag,o_concoa%name,"Concentration of Organic Aerosol ","kg/m3")
  CALL histdef3d(iff,o_concbc%flag,o_concbc%name,"Concentration of Black Carbon ","kg/m3")
  CALL histdef3d(iff,o_concss%flag,o_concss%name,"Concentration of Sea Salt ","kg/m3")
  CALL histdef3d(iff,o_concdust%flag,o_concdust%name,"Concentration of Dust ","kg/m3")
  CALL histdef2d(iff,o_loadso4%flag,o_loadso4%name,"Column Load of Sulfate ","kg/m2")
  CALL histdef2d(iff,o_loadoa%flag,o_loadoa%name,"Column Load of Organic Aerosol ","kg/m2")
  CALL histdef2d(iff,o_loadbc%flag,o_loadbc%name,"Column Load of Black Carbon ","kg/m2")
  CALL histdef2d(iff,o_loadss%flag,o_loadss%name,"Column Load of Sea Salt ","kg/m2")
  CALL histdef2d(iff,o_loaddust%flag,o_loaddust%name,"Column Load of Dust ","kg/m2")

  DO naero = 1, naero_spc
  CALL histdef2d(iff,o_tausumaero(naero)%flag,o_tausumaero(naero)%name,"Aerosol Optical depth at 550 nm "//name_aero(naero),"1")
  END DO
 ENDIF
ENDIF

 IF (ok_ade) THEN
  CALL histdef2d(iff,o_topswad%flag,o_topswad%name, "ADE at TOA", "W/m2")
  CALL histdef2d(iff,o_solswad%flag,o_solswad%name, "ADE at SRF", "W/m2")

 CALL histdef2d(iff,o_swtoaas_nat%flag,o_swtoaas_nat%name, "Natural aerosol radiative forcing all-sky at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfas_nat%flag,o_swsrfas_nat%name, "Natural aerosol radiative forcing all-sky at SRF", "W/m2")
 CALL histdef2d(iff,o_swtoacs_nat%flag,o_swtoacs_nat%name, "Natural aerosol radiative forcing clear-sky at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfcs_nat%flag,o_swsrfcs_nat%name, "Natural aerosol radiative forcing clear-sky at SRF", "W/m2")

 CALL histdef2d(iff,o_swtoaas_ant%flag,o_swtoaas_ant%name, "Anthropogenic aerosol radiative forcing all-sky at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfas_ant%flag,o_swsrfas_ant%name, "Anthropogenic aerosol radiative forcing all-sky at SRF", "W/m2")
 CALL histdef2d(iff,o_swtoacs_ant%flag,o_swtoacs_ant%name, "Anthropogenic aerosol radiative forcing clear-sky at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfcs_ant%flag,o_swsrfcs_ant%name, "Anthropogenic aerosol radiative forcing clear-sky at SRF", "W/m2")

 IF (.NOT. aerosol_couple) THEN 
 CALL histdef2d(iff,o_swtoacf_nat%flag,o_swtoacf_nat%name, "Natural aerosol impact on cloud radiative forcing at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfcf_nat%flag,o_swsrfcf_nat%name, "Natural aerosol impact on cloud radiative forcing  at SRF", "W/m2")
 CALL histdef2d(iff,o_swtoacf_ant%flag,o_swtoacf_ant%name, "Anthropogenic aerosol impact on cloud radiative forcing at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfcf_ant%flag,o_swsrfcf_ant%name, "Anthropogenic aerosol impact on cloud radiative forcing at SRF", "W/m2")
 CALL histdef2d(iff,o_swtoacf_zero%flag,o_swtoacf_zero%name, "Cloud radiative forcing (allsky-clearsky fluxes) at TOA", "W/m2")
 CALL histdef2d(iff,o_swsrfcf_zero%flag,o_swsrfcf_zero%name, "Cloud radiative forcing (allsky-clearsky fluxes) at SRF", "W/m2")
 ENDIF

 ENDIF

 IF (ok_aie) THEN
  CALL histdef2d(iff,o_topswai%flag,o_topswai%name, "AIE at TOA", "W/m2")
  CALL histdef2d(iff,o_solswai%flag,o_solswai%name, "AIE at SFR", "W/m2")
!Cloud droplet number concentration
  CALL histdef3d(iff,o_scdnc%flag,o_scdnc%name, "Cloud droplet number concentration","m-3")
  CALL histdef2d(iff,o_cldncl%flag,o_cldncl%name, "CDNC at top of liquid water cloud", "m-3")
  CALL histdef3d(iff,o_reffclws%flag,o_reffclws%name, "Stratiform Cloud Droplet Effective Radius","m")
  CALL histdef3d(iff,o_reffclwc%flag,o_reffclwc%name, "Convective Cloud Droplet Effective Radius","m")
  CALL histdef2d(iff,o_cldnvi%flag,o_cldnvi%name, "Column Integrated Cloud Droplet Number", "m-2")
  CALL histdef3d(iff,o_lcc3d%flag,o_lcc3d%name, "Cloud liquid fraction","1")
  CALL histdef3d(iff,o_lcc3dcon%flag,o_lcc3dcon%name, "Convective cloud liquid fraction","1")
  CALL histdef3d(iff,o_lcc3dstra%flag,o_lcc3dstra%name, "Stratiform cloud liquid fraction","1")
  CALL histdef2d(iff,o_lcc%flag,o_lcc%name, "Cloud liquid fraction at top of cloud","1")
  CALL histdef2d(iff,o_reffclwtop%flag,o_reffclwtop%name, "Droplet effective radius at top of liquid water cloud", "m")
 ENDIF


 CALL histdef2d(iff,o_alb1%flag,o_alb1%name, "Surface VIS albedo", "-")
 CALL histdef2d(iff,o_alb2%flag,o_alb2%name, "Surface Near IR albedo", "-")
 CALL histdef2d(iff,o_cdrm%flag,o_cdrm%name, "Momentum drag coef.", "-")
 CALL histdef2d(iff,o_cdrh%flag,o_cdrh%name, "Heat drag coef.", "-" )
 CALL histdef2d(iff,o_cldl%flag,o_cldl%name, "Low-level cloudiness", "-")
 CALL histdef2d(iff,o_cldm%flag,o_cldm%name, "Mid-level cloudiness", "-")
 CALL histdef2d(iff,o_cldh%flag,o_cldh%name, "High-level cloudiness", "-")
 CALL histdef2d(iff,o_cldt%flag,o_cldt%name, "Total cloudiness", "-")
 CALL histdef2d(iff,o_cldq%flag,o_cldq%name, "Cloud liquid water path", "kg/m2")
 CALL histdef2d(iff,o_lwp%flag,o_lwp%name, "Cloud water path", "kg/m2")
 CALL histdef2d(iff,o_iwp%flag,o_iwp%name, "Cloud ice water path", "kg/m2" )
 CALL histdef2d(iff,o_ue%flag,o_ue%name, "Zonal energy transport", "-")
 CALL histdef2d(iff,o_ve%flag,o_ve%name, "Merid energy transport", "-")
 CALL histdef2d(iff,o_uq%flag,o_uq%name, "Zonal humidity transport", "-")
 CALL histdef2d(iff,o_vq%flag,o_vq%name, "Merid humidity transport", "-")

     IF(iflag_con.GE.3) THEN ! sb
 CALL histdef2d(iff,o_cape%flag,o_cape%name, "Conv avlbl pot ener", "J/kg")
 CALL histdef2d(iff,o_pbase%flag,o_pbase%name, "Cld base pressure", "Pa")
 CALL histdef2d(iff,o_ptop%flag,o_ptop%name, "Cld top pressure", "Pa")
 CALL histdef2d(iff,o_fbase%flag,o_fbase%name, "Cld base mass flux", "kg/m2/s")
 CALL histdef2d(iff,o_prw%flag,o_prw%name, "Precipitable water", "kg/m2")
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)'
 CALL histdef2d(iff,o_cape_max%flag,o_cape_max%name, "CAPE max.", "J/kg")
   type_ecri(:) = type_ecri_files(:)
 CALL histdef3d(iff,o_upwd%flag,o_upwd%name, "saturated updraft", "kg/m2/s")
 CALL histdef3d(iff,o_Ma%flag,o_Ma%name, "undilute adiab updraft", "kg/m2/s")
 CALL histdef3d(iff,o_dnwd%flag,o_dnwd%name, "saturated downdraft", "kg/m2/s")
 CALL histdef3d(iff,o_dnwd0%flag,o_dnwd0%name, "unsat. downdraft", "kg/m2/s")
 CALL histdef3d(iff,o_mc%flag,o_mc%name, "Convective mass flux", "kg/m2/s")
   type_ecri(1) = 'inst(X)'
   type_ecri(2) = 'inst(X)'
   type_ecri(3) = 'inst(X)'
   type_ecri(4) = 'inst(X)'
   type_ecri(5) = 'inst(X)'
 CALL histdef2d(iff,o_ftime_con%flag,o_ftime_con%name, "Fraction of time convection Occurs", " ")
   type_ecri(:) = type_ecri_files(:)
     ENDIF !iflag_con .GE. 3

 CALL histdef2d(iff,o_s_pblh%flag,o_s_pblh%name, "Boundary Layer Height", "m")
 CALL histdef2d(iff,o_s_pblt%flag,o_s_pblt%name, "t at Boundary Layer Height", "K")
 CALL histdef2d(iff,o_s_lcl%flag,o_s_lcl%name, "Condensation level", "m")
 CALL histdef2d(iff,o_s_therm%flag,o_s_therm%name, "Exces du thermique", "K")
!IM : Les champs suivants (s_oliqCL, s_cteiCL, s_trmb1, s_trmb2, s_trmb3) ne sont pas definis dans HBTM.F
!CALL histdef2d(iff,o_s_capCL%flag,o_s_capCL%name, "Conv avlbl pot enerfor ABL", "J/m2" )
!CALL histdef2d(iff,o_s_oliqCL%flag,o_s_oliqCL%name, "Liq Water in BL", "kg/m2")
!CALL histdef2d(iff,o_s_cteiCL%flag,o_s_cteiCL%name, "Instability criteria(ABL)", "K")
!CALL histdef2d(iff,o_s_trmb1%flag,o_s_trmb1%name, "deep_cape(HBTM2)", "J/m2")
!CALL histdef2d(iff,o_s_trmb2%flag,o_s_trmb2%name, "inhibition (HBTM2)", "J/m2")
!CALL histdef2d(iff,o_s_trmb3%flag,o_s_trmb3%name, "Point Omega (HBTM2)", "m")

! Champs interpolles sur des niveaux de pression

   type_ecri(1) = 'inst(X)'
   type_ecri(2) = 'inst(X)'
   type_ecri(3) = 'inst(X)'
   type_ecri(4) = 'inst(X)'
   type_ecri(5) = 'inst(X)'

! Attention a reverifier

        ilev=0        
        DO k=1, nlevSTD
!     IF(k.GE.2.AND.k.LE.12) bb2=clevSTD(k)
     bb2=clevSTD(k)
     IF(bb2.EQ."850".OR.bb2.EQ."700".OR.bb2.EQ."500".OR.bb2.EQ."200".OR.bb2.EQ."50".OR.bb2.EQ."10")THEN
      ilev=ilev+1
      print*,'ilev k bb2 flag name ',ilev,k, bb2,o_uSTDlevs(ilev)%flag,o_uSTDlevs(ilev)%name
 CALL histdef2d(iff,o_uSTDlevs(ilev)%flag,o_uSTDlevs(ilev)%name,"Zonal wind "//bb2//"hPa", "m/s")
 CALL histdef2d(iff,o_vSTDlevs(ilev)%flag,o_vSTDlevs(ilev)%name,"Meridional wind "//bb2//"hPa", "m/s")
 CALL histdef2d(iff,o_wSTDlevs(ilev)%flag,o_wSTDlevs(ilev)%name,"Vertical wind "//bb2//"hPa", "Pa/s")
 CALL histdef2d(iff,o_zSTDlevs(ilev)%flag,o_zSTDlevs(ilev)%name,"Geopotential height "//bb2//"hPa", "m")
 CALL histdef2d(iff,o_qSTDlevs(ilev)%flag,o_qSTDlevs(ilev)%name,"Specific humidity "//bb2//"hPa", "kg/kg" )
 CALL histdef2d(iff,o_tSTDlevs(ilev)%flag,o_tSTDlevs(ilev)%name,"Temperature "//bb2//"hPa", "K")
     ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR."500".OR.bb2.EQ."200".OR.bb2.EQ."50".OR.bb2.EQ."10")
       ENDDO
   type_ecri(:) = type_ecri_files(:)

 CALL histdef2d(iff,o_t_oce_sic%flag,o_t_oce_sic%name, "Temp mixte oce-sic", "K")

 IF (type_ocean=='slab') & 
     CALL histdef2d(iff,o_slab_bils%flag, o_slab_bils%name,"Bilan au sol sur ocean slab", "W/m2")

! Couplage conv-CL
 IF (iflag_con.GE.3) THEN
    IF (iflag_coupl.EQ.1) THEN
 CALL histdef2d(iff,o_ale_bl%flag,o_ale_bl%name, "ALE BL", "m2/s2") 
 CALL histdef2d(iff,o_alp_bl%flag,o_alp_bl%name, "ALP BL", "m2/s2") 
    ENDIF
 ENDIF !(iflag_con.GE.3)

 CALL histdef2d(iff,o_weakinv%flag,o_weakinv%name, "Weak inversion", "-")
 CALL histdef2d(iff,o_dthmin%flag,o_dthmin%name, "dTheta mini", "K/m")
 CALL histdef2d(iff,o_rh2m%flag,o_rh2m%name, "Relative humidity at 2m", "%" )
   type_ecri(1) = 't_min(X)'
   type_ecri(2) = 't_min(X)'
   type_ecri(3) = 't_min(X)'
   type_ecri(4) = 't_min(X)'
   type_ecri(5) = 't_min(X)'
 CALL histdef2d(iff,o_rh2m_min%flag,o_rh2m_min%name, "Min Relative humidity at 2m", "%" )
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)'
 CALL histdef2d(iff,o_rh2m_max%flag,o_rh2m_max%name, "Max Relative humidity at 2m", "%" )
   type_ecri(:) = type_ecri_files(:)
 CALL histdef2d(iff,o_qsat2m%flag,o_qsat2m%name, "Saturant humidity at 2m", "%")
 CALL histdef2d(iff,o_tpot%flag,o_tpot%name, "Surface air potential temperature", "K")
 CALL histdef2d(iff,o_tpote%flag,o_tpote%name, "Surface air equivalent potential temperature", "K")
 CALL histdef2d(iff,o_SWnetOR%flag,o_SWnetOR%name, "Sfce net SW radiation OR", "W/m2")
 CALL histdef2d(iff,o_SWdownOR%flag,o_SWdownOR%name, "Sfce incident SW radiation OR", "W/m2")
 CALL histdef2d(iff,o_LWdownOR%flag,o_LWdownOR%name, "Sfce incident LW radiation OR", "W/m2")
 CALL histdef2d(iff,o_snowl%flag,o_snowl%name, "Solid Large-scale Precip.", "kg/(m2*s)")

 CALL histdef2d(iff,o_solldown%flag,o_solldown%name, "Down. IR rad. at surface", "W/m2")
 CALL histdef2d(iff,o_dtsvdfo%flag,o_dtsvdfo%name, "Boundary-layer dTs(o)", "K/s")
 CALL histdef2d(iff,o_dtsvdft%flag,o_dtsvdft%name, "Boundary-layer dTs(t)", "K/s")
 CALL histdef2d(iff,o_dtsvdfg%flag,o_dtsvdfg%name, "Boundary-layer dTs(g)", "K/s")
 CALL histdef2d(iff,o_dtsvdfi%flag,o_dtsvdfi%name, "Boundary-layer dTs(g)", "K/s")
 CALL histdef2d(iff,o_rugs%flag,o_rugs%name, "rugosity", "-" )

! Champs 3D:
 CALL histdef3d(iff,o_ec550aer%flag,o_ec550aer%name, "Extinction at 550nm", "m^-1")
 CALL histdef3d(iff,o_lwcon%flag,o_lwcon%name, "Cloud liquid water content", "kg/kg")
 CALL histdef3d(iff,o_iwcon%flag,o_iwcon%name, "Cloud ice water content", "kg/kg")
 CALL histdef3d(iff,o_temp%flag,o_temp%name, "Air temperature", "K" )
 CALL histdef3d(iff,o_theta%flag,o_theta%name, "Potential air temperature", "K" )
 CALL histdef3d(iff,o_ovap%flag,o_ovap%name, "Specific humidity", "kg/kg" )
 CALL histdef3d(iff,o_ovapinit%flag,o_ovapinit%name, "Specific humidity (begin of timestep)", "kg/kg" )
 CALL histdef3d(iff,o_geop%flag,o_geop%name, "Geopotential height", "m2/s2")
 CALL histdef3d(iff,o_vitu%flag,o_vitu%name, "Zonal wind", "m/s" )
 CALL histdef3d(iff,o_vitv%flag,o_vitv%name, "Meridional wind", "m/s" )
 CALL histdef3d(iff,o_vitw%flag,o_vitw%name, "Vertical wind", "Pa/s" )
 CALL histdef3d(iff,o_pres%flag,o_pres%name, "Air pressure", "Pa" )
 CALL histdef3d(iff,o_paprs%flag,o_paprs%name, "Air pressure Inter-Couches", "Pa" )
 CALL histdef3d(iff,o_mass%flag,o_mass%name, "Masse Couches", "kg/m2" )
 CALL histdef3d(iff,o_rneb%flag,o_rneb%name, "Cloud fraction", "-")
 CALL histdef3d(iff,o_rnebcon%flag,o_rnebcon%name, "Convective Cloud Fraction", "-")
 CALL histdef3d(iff,o_rhum%flag,o_rhum%name, "Relative humidity", "-")
 CALL histdef3d(iff,o_ozone%flag,o_ozone%name, "Ozone mole fraction", "-")
 if (read_climoz == 2) &
      CALL histdef3d(iff,o_ozone_light%flag,o_ozone_light%name, &
      "Daylight ozone mole fraction", "-")
 CALL histdef3d(iff,o_dtphy%flag,o_dtphy%name, "Physics dT", "K/s")
 CALL histdef3d(iff,o_dqphy%flag,o_dqphy%name, "Physics dQ", "(kg/kg)/s")
 CALL histdef3d(iff,o_cldtau%flag,o_cldtau%name, "Cloud optical thickness", "1")
 CALL histdef3d(iff,o_cldemi%flag,o_cldemi%name, "Cloud optical emissivity", "1")
!IM: bug ?? dimensionnement variables (klon,klev+1) pmflxr, pmflxs, prfl, psfl
 CALL histdef3d(iff,o_pr_con_l%flag,o_pr_con_l%name, "Convective precipitation lic", " ")
 CALL histdef3d(iff,o_pr_con_i%flag,o_pr_con_i%name, "Convective precipitation ice", " ")
 CALL histdef3d(iff,o_pr_lsc_l%flag,o_pr_lsc_l%name, "Large scale precipitation lic", " ")
 CALL histdef3d(iff,o_pr_lsc_i%flag,o_pr_lsc_i%name, "Large scale precipitation ice", " ")
!Cloud droplet effective radius
 CALL histdef3d(iff,o_re%flag,o_re%name, "Cloud droplet effective radius","um")
 CALL histdef3d(iff,o_fl%flag,o_fl%name, "Denominator of Cloud droplet effective radius"," ")
!FH Sorties pour la couche limite
     if (iflag_pbl>1) then
 CALL histdef3d(iff,o_tke%flag,o_tke%name, "TKE", "m2/s2")
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)'
 CALL histdef3d(iff,o_tke_max%flag,o_tke_max%name, "TKE max", "m2/s2")
   type_ecri(:) = type_ecri_files(:)
     endif

 CALL histdef3d(iff,o_kz%flag,o_kz%name, "Kz melange", "m2/s")
   type_ecri(1) = 't_max(X)'
   type_ecri(2) = 't_max(X)'
   type_ecri(3) = 't_max(X)'
   type_ecri(4) = 't_max(X)'
   type_ecri(5) = 't_max(X)'
 CALL histdef3d(iff,o_kz_max%flag,o_kz_max%name, "Kz melange max", "m2/s" )
   type_ecri(:) = type_ecri_files(:)
 CALL histdef3d(iff,o_clwcon%flag,o_clwcon%name, "Convective Cloud Liquid water content", "kg/kg") 
 CALL histdef3d(iff,o_dtdyn%flag,o_dtdyn%name, "Dynamics dT", "K/s")
 CALL histdef3d(iff,o_dqdyn%flag,o_dqdyn%name, "Dynamics dQ", "(kg/kg)/s")
 CALL histdef3d(iff,o_dudyn%flag,o_dudyn%name, "Dynamics dU", "m/s2")
 CALL histdef3d(iff,o_dvdyn%flag,o_dvdyn%name, "Dynamics dV", "m/s2")
 CALL histdef3d(iff,o_dtcon%flag,o_dtcon%name, "Convection dT", "K/s")
 CALL histdef3d(iff,o_ducon%flag,o_ducon%name, "Convection du", "m/s2")
 CALL histdef3d(iff,o_dqcon%flag,o_dqcon%name, "Convection dQ", "(kg/kg)/s")

! Wakes
 IF(iflag_con.EQ.3) THEN
 IF (iflag_wake == 1) THEN
   CALL histdef2d(iff,o_ale_wk%flag,o_ale_wk%name, "ALE WK", "m2/s2")
   CALL histdef2d(iff,o_alp_wk%flag,o_alp_wk%name, "ALP WK", "m2/s2")
   CALL histdef2d(iff,o_ale%flag,o_ale%name, "ALE", "m2/s2")
   CALL histdef2d(iff,o_alp%flag,o_alp%name, "ALP", "W/m2")
   CALL histdef2d(iff,o_cin%flag,o_cin%name, "Convective INhibition", "m2/s2")
   CALL histdef2d(iff,o_wape%flag,o_WAPE%name, "WAPE", "m2/s2")
   CALL histdef2d(iff,o_wake_h%flag,o_wake_h%name, "wake_h", "-")
   CALL histdef2d(iff,o_wake_s%flag,o_wake_s%name, "wake_s", "-")
   CALL histdef3d(iff,o_dtwak%flag,o_dtwak%name, "Wake dT", "K/s")
   CALL histdef3d(iff,o_dqwak%flag,o_dqwak%name, "Wake dQ", "(kg/kg)/s")
   CALL histdef3d(iff,o_wake_deltat%flag,o_wake_deltat%name, "wake_deltat", " ")
   CALL histdef3d(iff,o_wake_deltaq%flag,o_wake_deltaq%name, "wake_deltaq", " ")
   CALL histdef3d(iff,o_wake_omg%flag,o_wake_omg%name, "wake_omg", " ")
 ENDIF
   CALL histdef3d(iff,o_Vprecip%flag,o_Vprecip%name, "precipitation vertical profile", "-")
   CALL histdef3d(iff,o_ftd%flag,o_ftd%name, "tend temp due aux descentes precip", "-")
   CALL histdef3d(iff,o_fqd%flag,o_fqd%name,"tend vap eau due aux descentes precip", "-")
 ENDIF !(iflag_con.EQ.3)

 CALL histdef3d(iff,o_dtlsc%flag,o_dtlsc%name, "Condensation dT", "K/s")
 CALL histdef3d(iff,o_dtlschr%flag,o_dtlschr%name,"Large-scale condensational heating rate","K/s")
 CALL histdef3d(iff,o_dqlsc%flag,o_dqlsc%name, "Condensation dQ", "(kg/kg)/s")
 CALL histdef3d(iff,o_dtvdf%flag,o_dtvdf%name, "Boundary-layer dT", "K/s")
 CALL histdef3d(iff,o_dqvdf%flag,o_dqvdf%name, "Boundary-layer dQ", "(kg/kg)/s") 
 CALL histdef3d(iff,o_dteva%flag,o_dteva%name, "Reevaporation dT", "K/s")
 CALL histdef3d(iff,o_dqeva%flag,o_dqeva%name, "Reevaporation dQ", "(kg/kg)/s")
 CALL histdef3d(iff,o_ptconv%flag,o_ptconv%name, "POINTS CONVECTIFS", " ")
 CALL histdef3d(iff,o_ratqs%flag,o_ratqs%name, "RATQS", " ")
 CALL histdef3d(iff,o_dtthe%flag,o_dtthe%name, "Dry adjust. dT", "K/s")

if(iflag_thermals.gt.1) THEN
 CALL histdef3d(iff,o_f_th%flag,o_f_th%name, "Thermal plume mass flux", "K/s")
 CALL histdef3d(iff,o_e_th%flag,o_e_th%name,"Thermal plume entrainment","K/s")
 CALL histdef3d(iff,o_w_th%flag,o_w_th%name,"Thermal plume vertical velocity","m/s")
 CALL histdef3d(iff,o_lambda_th%flag,o_lambda_th%name,"Thermal plume vertical velocity","m/s")
 CALL histdef2d(iff,o_ftime_th%flag,o_ftime_th%name,"Fraction of time Shallow convection occurs"," ")
 CALL histdef3d(iff,o_q_th%flag,o_q_th%name, "Thermal plume total humidity", "kg/kg")
 CALL histdef3d(iff,o_a_th%flag,o_a_th%name, "Thermal plume fraction", "")
 CALL histdef3d(iff,o_d_th%flag,o_d_th%name, "Thermal plume detrainment", "K/s")
!IM endif !iflag_thermals.gt.1
 CALL histdef2d(iff,o_f0_th%flag,o_f0_th%name, "Thermal closure mass flux", "K/s")
 CALL histdef2d(iff,o_zmax_th%flag,o_zmax_th%name, "Thermal plume height", "K/s")
 CALL histdef3d(iff,o_dqthe%flag,o_dqthe%name, "Dry adjust. dQ", "(kg/kg)/s")
endif !iflag_thermals.gt.1
 CALL histdef3d(iff,o_dtajs%flag,o_dtajs%name, "Dry adjust. dT", "K/s")
 CALL histdef3d(iff,o_dqajs%flag,o_dqajs%name, "Dry adjust. dQ", "(kg/kg)/s")
 CALL histdef3d(iff,o_dtswr%flag,o_dtswr%name, "SW radiation dT", "K/s")
 CALL histdef3d(iff,o_dtsw0%flag,o_dtsw0%name, "CS SW radiation dT", "K/s")
 CALL histdef3d(iff,o_dtlwr%flag,o_dtlwr%name, "LW radiation dT", "K/s")
 CALL histdef3d(iff,o_dtlw0%flag,o_dtlw0%name, "CS LW radiation dT", "K/s")
 CALL histdef3d(iff,o_dtec%flag,o_dtec%name, "Cinetic dissip dT", "K/s")
 CALL histdef3d(iff,o_duvdf%flag,o_duvdf%name, "Boundary-layer dU", "m/s2")
 CALL histdef3d(iff,o_dvvdf%flag,o_dvvdf%name, "Boundary-layer dV", "m/s2")

     IF (ok_orodr) THEN
 CALL histdef3d(iff,o_duoro%flag,o_duoro%name, "Orography dU", "m/s2")
 CALL histdef3d(iff,o_dvoro%flag,o_dvoro%name, "Orography dV", "m/s2")
     ENDIF

     IF (ok_orolf) THEN
 CALL histdef3d(iff,o_dulif%flag,o_dulif%name, "Orography dU", "m/s2")
 CALL histdef3d(iff,o_dvlif%flag,o_dvlif%name, "Orography dV", "m/s2")
     ENDIF


    IF (nqtot>=3) THEN
     DO iq=3,nqtot  
       iiq=niadv(iq)
       o_trac(iq-2) = ctrl_out((/ 4, 5, 1, 1, 1 /),tname(iiq))
       CALL histdef3d (iff, o_trac(iq-2)%flag,o_trac(iq-2)%name,'Tracer '//ttext(iiq), "-" )
     ENDDO
    ENDIF

        CALL histend(nid_files(iff))

         ndex2d = 0
         ndex3d = 0

         ENDIF ! clef_files

         ENDDO !  iff
     print*,'Fin phys_output_mod.F90'
      end subroutine phys_output_open

      SUBROUTINE histdef2d (iff,flag_var,nomvar,titrevar,unitvar)
      
       use ioipsl
       USE dimphy
       USE mod_phys_lmdz_para

       IMPLICIT NONE
       
       include "dimensions.h"
       include "temps.h"
       include "indicesol.h"
       include "clesphys.h"

       integer                          :: iff
       integer, dimension(nfiles)       :: flag_var
       character(len=20)                 :: nomvar
       character(len=*)                 :: titrevar
       character(len=*)                 :: unitvar

       real zstophym

       if (type_ecri(iff)=='inst(X)'.OR.type_ecri(iff)=='once') then
         zstophym=zoutm(iff)
       else
         zstophym=zdtime
       endif

! Appel a la lecture des noms et niveau d'ecriture des variables dans output.def
       call conf_physoutputs(nomvar,flag_var)
       
       if ( flag_var(iff)<=lev_files(iff) ) then
 call histdef (nid_files(iff),nomvar,titrevar,unitvar, &
               iim,jj_nb,nhorim(iff), 1,1,1, -99, 32, &
               type_ecri(iff), zstophym,zoutm(iff))                
       endif                      
      end subroutine histdef2d

      SUBROUTINE histdef3d (iff,flag_var,nomvar,titrevar,unitvar)

       use ioipsl
       USE dimphy
       USE mod_phys_lmdz_para

       IMPLICIT NONE

       include "dimensions.h"
       include "temps.h"
       include "indicesol.h"
       include "clesphys.h"

       integer                          :: iff
       integer, dimension(nfiles)       :: flag_var
       character(len=20)                 :: nomvar
       character(len=*)                 :: titrevar
       character(len=*)                 :: unitvar

       real zstophym

! Appel a la lecture des noms et niveau d'ecriture des variables dans output.def
       call conf_physoutputs(nomvar,flag_var)

       if (type_ecri(iff)=='inst(X)'.OR.type_ecri(iff)=='once') then
         zstophym=zoutm(iff)
       else
         zstophym=zdtime
       endif

       if ( flag_var(iff)<=lev_files(iff) ) then
          call histdef (nid_files(iff), nomvar, titrevar, unitvar, &
               iim, jj_nb, nhorim(iff), klev, levmin(iff), &
               levmax(iff)-levmin(iff)+1, nvertm(iff), 32, type_ecri(iff), &
               zstophym, zoutm(iff))
       endif
      end subroutine histdef3d

      SUBROUTINE conf_physoutputs(nam_var,flag_var)
!!! Lecture des noms et niveau de sortie des variables dans output.def
!   en utilisant les routines getin de IOIPSL  
       use ioipsl

       IMPLICIT NONE

       include 'iniprint.h'

       character(len=20)                :: nam_var
       integer, dimension(nfiles)      :: flag_var

        IF(prt_level>10) WRITE(lunout,*)'Avant getin: nam_var flag_var ',nam_var,flag_var(:)
        call getin('flag_'//nam_var,flag_var)
        call getin('name_'//nam_var,nam_var)
        IF(prt_level>10) WRITE(lunout,*)'Apres getin: nam_var flag_var ',nam_var,flag_var(:)

      END SUBROUTINE conf_physoutputs

      SUBROUTINE convers_timesteps(str,timestep)

        use ioipsl

        IMPLICIT NONE

        character(len=20)   :: str
        character(len=10)   :: type
        integer             :: ipos,il
        real                :: ttt,xxx,timestep,dayseconde
        parameter (dayseconde=86400.)
        include "temps.h"
        include "comconst.h"

        ipos=scan(str,'0123456789.',.true.)
!  
        il=len_trim(str)
        print*,ipos,il
        read(str(1:ipos),*) ttt
        print*,ttt
        type=str(ipos+1:il)


        if ( il == ipos ) then
        type='day'
        endif

        if ( type == 'day'.or.type == 'days'.or.type == 'jours'.or.type == 'jour' ) timestep = ttt * dayseconde
        if ( type == 'mounths'.or.type == 'mth'.or.type == 'mois' ) then
           print*,'annee_ref,day_ref mon_len',annee_ref,day_ref,ioget_mon_len(annee_ref,day_ref)
           timestep = ttt * dayseconde * ioget_mon_len(annee_ref,day_ref)
        endif
        if ( type == 'hours'.or.type == 'hr'.or.type == 'heurs') timestep = ttt * dayseconde / 24.
        if ( type == 'mn'.or.type == 'minutes'  ) timestep = ttt * 60.
        if ( type == 's'.or.type == 'sec'.or.type == 'secondes'   ) timestep = ttt
        if ( type == 'TS' ) timestep = dtphys

        print*,'type =      ',type
        print*,'nb j/h/m =  ',ttt
        print*,'timestep(s)=',timestep

        END SUBROUTINE convers_timesteps

END MODULE phys_output_mod

