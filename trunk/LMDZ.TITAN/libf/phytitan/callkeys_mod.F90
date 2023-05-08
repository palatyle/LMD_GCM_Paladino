MODULE callkeys_mod
IMPLICIT NONE  

      logical,save :: callrad,corrk,calldifv,UseTurbDiff
!$OMP THREADPRIVATE(callrad,corrk,calldifv,UseTurbDiff)
      logical,save :: calladj,callsoil
!$OMP THREADPRIVATE(calladj,callsoil)
      logical,save :: season,diurnal,tlocked,rings_shadow,lwrite
!$OMP THREADPRIVATE(season,diurnal,tlocked,rings_shadow,lwrite)
      logical,save :: callstats,calleofdump
!$OMP THREADPRIVATE(callstats,calleofdump)
      logical,save :: callgasvis,continuum,graybody
!$OMP THREADPRIVATE(callgasvis,continuum,graybody)
      logical,save :: strictboundcorrk
!$OMP THREADPRIVATE(strictboundcorrk)
      logical,save :: corrk_recombin
!$OMP_THREADPRIVATE(corrk_recombin)
      logical,save :: seashaze,uncoupl_optic_haze
!$OMP THREADPRIVATE(seashaze,uncoupl_optic_haze)

      logical,save :: callchim, callmufi, callclouds
!$OMP THREADPRIVATE(callchim,callmufi,callclouds)
      logical,save :: global1d
!$OMP THREADPRIVATE(global1d)
      logical,save :: enertest
      logical,save :: nonideal
      logical,save :: meanOLR
      logical,save :: specOLR
      logical,save :: diagdtau
!$OMP THREADPRIVATE(enertest,nonideal,meanOLR,specOLR,diagdtau)
      logical,save :: newtonian
      logical,save :: check_cpp_match
      logical,save :: force_cpp
      logical,save :: testradtimes
      logical,save :: rayleigh
!$OMP THREADPRIVATE(newtonian,check_cpp_match,force_cpp,testradtimes,rayleigh)
      logical,save :: stelbbody
      logical,save :: ozone
      logical,save :: tracer
      logical,save :: mass_redistrib
!$OMP THREADPRIVATE(stelbbody,ozone,tracer,mass_redistrib)
      logical,save :: nosurf
      logical,save :: oblate
!$OMP THREADPRIVATE(nosurf,oblate)
      logical,save :: eff_gz
!$OMP THREADPRIVATE(eff_gz)
      
      integer,save :: ichim
!$OMP THREADPRIVATE(ichim)
      integer,save :: versH2H2cia
      integer,save :: iddist
      integer,save :: iradia
      integer,save :: startype
!$OMP THREADPRIVATE(versH2H2cia,iddist,iradia,startype)
      
      real,save :: p_prod, tx_prod, rc_prod
      real,save :: air_rad
!$OMP THREADPRIVATE(p_prod,tx_prod,rc_prod,air_rad)
      
      real,save :: szangle
!$OMP THREADPRIVATE(szangle)
      real,save :: Fat1AU
      real,save :: stelTbb
!$OMP THREADPRIVATE(Fat1AU,stelTbb)
      real,save :: pceil
!$OMP THREADPRIVATE(pceil)
      real,save :: tau_relax
      real,save :: intheat
!$OMP THREADPRIVATE(tau_relax,intheat)
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
