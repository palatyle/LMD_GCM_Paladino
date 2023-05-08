MODULE callkeys_mod
IMPLICIT NONE  

      logical,save :: callrad,corrk,calldifv,UseTurbDiff
!$OMP THREADPRIVATE(callrad,corrk,calldifv,UseTurbDiff)
      !! Paladino
      logical,save :: calladj,calltherm,co2cond,callsoil,callvolcano
      character(50),save :: volc_name,input_key
      character(15),save :: atmos_type
      character(300),save :: output_dir
      real,save :: rho_volc 
!$OMP THREADPRIVATE(calladj,calltherm,co2cond,callsoil,callvolcano)
      logical,save :: season,diurnal,tlocked,rings_shadow,lwrite
!$OMP THREADPRIVATE(season,diurnal,tlocked,rings_shadow,lwrite)
      logical,save :: callstats,calleofdump
!$OMP THREADPRIVATE(callstats,calleofdump)
      logical,save :: callgasvis,continuum,H2Ocont_simple,graybody
!$OMP THREADPRIVATE(callgasvis,continuum,H2Ocont_simple,graybody)
      logical,save :: strictboundcorrk                                     
!$OMP THREADPRIVATE(strictboundcorrk)

      logical,save :: enertest
      logical,save :: nonideal
      logical,save :: meanOLR
      logical,save :: specOLR
      logical,save :: kastprof
      logical,save :: diagdtau
!$OMP THREADPRIVATE(enertest,nonideal,meanOLR,kastprof,diagdtau)
      logical,save :: newtonian
      logical,save :: check_cpp_match
      logical,save :: force_cpp
      logical,save :: testradtimes
      logical,save :: rayleigh
!$OMP THREADPRIVATE(newtonian,check_cpp_match,force_cpp,testradtimes,rayleigh)
      logical,save :: stelbbody
      logical,save :: ozone
      logical,save :: nearco2cond
      logical,save :: tracer
      logical,save :: mass_redistrib
!$OMP THREADPRIVATE(stelbbody,ozone,nearco2cond,tracer,mass_redistrib)
      logical,save :: varactive
      logical,save :: varfixed
      logical,save :: radfixed
      logical,save :: sedimentation
!$OMP THREADPRIVATE(varactive,varfixed,radfixed,sedimentation)
      logical,save :: water,watercond,waterrain
!$OMP THREADPRIVATE(water,watercond,waterrain)
      logical,save :: aeroco2,aeroh2o,aeroh2so4,aeroback2lay
!$OMP THREADPRIVATE(aeroco2,aeroh2o,aeroh2so4,aeroback2lay)
      logical,save :: aeronh3, aeroaurora
!$OMP THREADPRIVATE(aeronh3,aeroaurora)
      logical,save :: aerofixco2,aerofixh2o
!$OMP THREADPRIVATE(aerofixco2,aerofixh2o)
      logical,save :: hydrology
      logical,save :: sourceevol
      logical,save :: CLFvarying
      logical,save :: nosurf
      logical,save :: oblate
!$OMP THREADPRIVATE(hydrology,sourceevol,CLFvarying,nosurf,oblate)
      logical,save :: ok_slab_ocean
      logical,save :: ok_slab_sic
      logical,save :: ok_slab_heat_transp
      logical,save :: albedo_spectral_mode
!$OMP THREADPRIVATE(ok_slab_ocean,ok_slab_sic,ok_slab_heat_transp,albedo_spectral_mode)
      logical,save :: photochem
      logical,save :: haze
!$OMP THREADPRIVATE(photochem)

      integer,save :: versH2H2cia
      integer,save :: iddist
      integer,save :: iaervar
      integer,save :: iradia
      integer,save :: startype
!$OMP THREADPRIVATE(versH2H2cia,iddist,iaervar,iradia,startype)

      real,save :: tplanckmin
      real,save :: tplanckmax
      real,save :: dtplanck
!$OMP THREADPRIVATE(tplanckmin,tplanckmax,dtplanck)
      real,save :: topdustref
      real,save :: Nmix_co2
      real,save :: dusttau
      real,save :: Fat1AU
      real,save :: stelTbb
!$OMP THREADPRIVATE(topdustref,Nmix_co2,dusttau,Fat1AU,stelTbb)
      real,save :: Tstrat
      real,save :: tplanet
      real,save :: obs_tau_col_tropo
      real,save :: obs_tau_col_strato
!$OMP THREADPRIVATE(Tstrat,tplanet,obs_tau_col_tropo,obs_tau_col_strato)
      character(64),save :: optprop_back2lay_vis
      character(64),save :: optprop_back2lay_ir
      real,save :: pres_bottom_tropo
      real,save :: pres_top_tropo
      real,save :: pres_bottom_strato
      real,save :: pres_top_strato
!$OMP THREADPRIVATE(pres_bottom_tropo,pres_top_tropo,pres_bottom_strato,pres_top_strato)
      real,save :: size_tropo
      real,save :: size_strato
      real,save :: satval
      real,save :: CLFfixval
      real,save :: n2mixratio
!$OMP THREADPRIVATE(size_tropo,size_strato,satval,CLFfixval,n2mixratio)
      real,save :: size_nh3_cloud
      real,save :: pres_nh3_cloud
      real,save :: tau_nh3_cloud
!$OMP THREADPRIVATE(size_nh3_cloud, pres_nh3_cloud, tau_nh3_cloud)
      real,save :: co2supsat
      real,save :: pceil
      real,save :: albedosnow
      real,save :: albedoco2ice
      real,save :: maxicethick
!$OMP THREADPRIVATE(co2supsat,pceil,albedosnow,albedoco2ice,maxicethick)
      real,save :: Tsaldiff
      real,save :: tau_relax
      real,save :: cloudlvl
      real,save :: icetstep
      real,save :: intheat
!$OMP THREADPRIVATE(Tsaldiff,tau_relax,cloudlvl,icetstep,intheat)
      real,save :: flatten
      real,save :: Rmean
      real,save :: J2
      real,save :: MassPlanet
!$OMP THREADPRIVATE(flatten,Rmean,J2,MassPlanet)
      real,save :: surfalbedo
      real,save :: surfemis
!$OMP THREADPRIVATE(surfalbedo,surfemis)

      logical,save :: iscallphys=.false.!existence of callphys.def
!$OMP THREADPRIVATE(iscallphys)

      ! do we read a startphy.nc file (default=.true.)
      logical,save :: startphy_file=.true. 
!$OMP THREADPRIVATE(startphy_file)

END MODULE callkeys_mod
