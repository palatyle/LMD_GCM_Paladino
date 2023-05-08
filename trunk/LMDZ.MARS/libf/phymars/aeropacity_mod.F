      MODULE aeropacity_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE aeropacity(ngrid,nlayer,nq,zday,pplay,pplev,ls,
     &    pq,tauscaling,tauref,tau,taucloudtes,aerosol,dsodust,reffrad,
     &    QREFvis3d,QREFir3d,omegaREFir3d,
     &    totstormfract,clearatm,dsords,dsotop,
     &    alpha_hmons,nohmons,
     &    clearsky,totcloudfrac)
                                                         
! to use  'getin'
      USE ioipsl_getincom, only: getin
      use tracer_mod, only: noms, igcm_h2o_ice, igcm_dust_mass,
     &                      igcm_dust_submicron, rho_dust, rho_ice,
     &                      nqdust, igcm_stormdust_mass,
     &                      igcm_topdust_mass
      use geometry_mod, only: latitude ! grid point latitudes (rad)
      use comgeomfi_h, only: sinlat ! sines of grid point latitudes
#ifdef DUSTSTORM
      use geometry_mod, only: longitude
      use tracer_mod, only: r3n_q, ref_r0, igcm_dust_number
#endif
      use planete_h
      USE comcstfi_h
      use dimradmars_mod, only: naerkind, name_iaer, 
     &            iaerdust,tauvis,
     &            iaer_dust_conrath,iaer_dust_doubleq,
     &            iaer_dust_submicron,iaer_h2o_ice,
     &            iaer_stormdust_doubleq,
     &            iaer_topdust_doubleq
      USE calcstormfract_mod
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing aerosol optical depth in each gridbox.
c
c   author: F.Forget 
c   ------
c   update F. Montmessin (water ice scheme) 
c      and S. Lebonnois (12/06/2003) compatibility dust/ice/chemistry
c   update J.-B. Madeleine 2008-2009:
c       - added 3D scattering by aerosols;
c       - dustopacity transferred from physiq.F to callradite.F,
c           and renamed into aeropacity.F;
c   update E. Millour, march 2012:
c         - reference pressure is now set to 610Pa (not 700Pa)
c   
c   input:
c   ----- 
c   ngrid             Number of gridpoint of horizontal grid
c   nlayer            Number of layer
c   nq                Number of tracer
c   zday                  Date (time since Ls=0, in martian days)
c   ls                Solar longitude (Ls) , radian
c   pplay,pplev       pressure (Pa) in the middle and boundary of each layer
c   pq                Dust mixing ratio (used if tracer =T and active=T).
c   reffrad(ngrid,nlayer,naerkind)  Aerosol effective radius
c   QREFvis3d(ngrid,nlayer,naerkind) \ 3d extinction coefficients
c   QREFir3d(ngrid,nlayer,naerkind)  / at reference wavelengths;
c   omegaREFir3d(ngrid,nlayer,naerkind)  / at reference wavelengths;
c
c   output:
c   -------
c   tauref       Prescribed mean column optical depth at 610 Pa 
c   tau          Column total visible dust optical depth at each point
c   aerosol      aerosol(ig,l,1) is the dust optical
c                depth in layer l, grid point ig
c   taualldust   CW17 total opacity for all dust scatterer stormdust included
c
c=======================================================================
      include "callkeys.h"

c-----------------------------------------------------------------------
c
c    Declarations :
c    --------------
c
c    Input/Output
c    ------------
      INTEGER, INTENT(IN) ::  ngrid,nlayer,nq
      REAL, INTENT(IN) :: ls,zday
      REAL, INTENT(IN) ::  pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL, INTENT(IN) ::  pq(ngrid,nlayer,nq)
      REAL, INTENT(OUT) :: tauref(ngrid)
      REAL, INTENT(OUT) :: tau(ngrid,naerkind)
      REAL, INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind)
      REAL, INTENT(INOUT) ::  dsodust(ngrid,nlayer)
      REAL, INTENT(INOUT) ::  dsords(ngrid,nlayer) !dso of stormdust
      REAL, INTENT(INOUT) ::  dsotop(ngrid,nlayer) !dso of topdust  
      REAL, INTENT(INOUT) :: reffrad(ngrid,nlayer,naerkind)
      REAL, INTENT(IN) :: QREFvis3d(ngrid,nlayer,naerkind)
      REAL, INTENT(IN) :: QREFir3d(ngrid,nlayer,naerkind)
      REAL, INTENT(IN) ::  omegaREFir3d(ngrid,nlayer,naerkind)
      LOGICAL, INTENT(IN) :: clearatm
      REAL, INTENT(IN) :: totstormfract(ngrid)
      LOGICAL, INTENT(IN) :: nohmons
      REAL, INTENT(IN) :: alpha_hmons(ngrid)
      REAL, INTENT(OUT) ::  tauscaling(ngrid) ! Scaling factor for qdust and Ndust
      REAL,INTENT(IN) :: totcloudfrac(ngrid) ! total cloud fraction
      LOGICAL,INTENT(IN) :: clearsky ! true for part without clouds,false for part with clouds (total or sub-grid clouds)
c
c    Local variables :
c    -----------------
      REAL CLFtot ! total cloud fraction
      real expfactor 
      INTEGER l,ig,iq,i,j
      INTEGER iaer           ! Aerosol index
      real topdust(ngrid)
      real zlsconst, zp
      real taueq,tauS,tauN
c     Mean Qext(vis)/Qext(ir) profile
      real msolsir(nlayer,naerkind)
c     Mean Qext(ir)/Qabs(ir) profile
      real mqextsqabs(nlayer,naerkind)
c     Variables used when multiple particle sizes are used
c       for dust or water ice particles in the radiative transfer
c       (see callradite.F for more information).
      REAL taudusttmp(ngrid)! Temporary dust opacity used before scaling
      REAL taubackdusttmp(ngrid)! Temporary background dust opacity used before scaling
      REAL taualldust(ngrid)! dust opacity all dust
      REAL taudust(ngrid)! dust opacity dust doubleq
      REAL taustormdust(ngrid)! dust opacity stormdust doubleq
      REAL taustormdusttmp(ngrid)! dust opacity stormdust doubleq before tauscaling
      REAL taudustvis(ngrid) ! Dust opacity after scaling
      REAL taudusttes(ngrid) ! Dust opacity at IR ref. wav. as
                               !   "seen" by the GCM.
      REAL taucloudvis(ngrid)! Cloud opacity at visible
                               !   reference wavelength
      REAL taucloudtes(ngrid)! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)
      REAL topdust0(ngrid)

#ifdef DUSTSTORM
!! Local dust storms
      logical localstorm        ! =true to create a local dust storm
      real taulocref,ztoploc,radloc,lonloc,latloc  ! local dust storm parameters
      real reffstorm, yeah
      REAL ray(ngrid) ! distance from dust storm center
      REAL tauuser(ngrid) ! opacity perturbation due to dust storm
      REAL more_dust(ngrid,nlayer,2) ! Mass mixing ratio perturbation due to the dust storm
      REAL int_factor(ngrid) ! useful factor to compute mmr perturbation
      real l_top ! layer of the storm's top
      REAL zalt(ngrid, nlayer) ! useful factor to compute l_top 
#endif

c   local saved variables
c   ---------------------

c     Level under which the dust mixing ratio is held constant
c       when computing the dust opacity in each layer
c       (this applies when doubleq and active are true)
      INTEGER, PARAMETER :: cstdustlevel0 = 7
      INTEGER, SAVE      :: cstdustlevel

      LOGICAL,SAVE :: firstcall=.true.

! indexes of water ice and dust tracers:
      INTEGER,SAVE :: i_ice=0  ! water ice
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
      CHARACTER(LEN=20) :: txt ! to temporarly store text
      CHARACTER(LEN=1) :: txt2 ! to temporarly store text
! indexes of dust scatterers:
      INTEGER,SAVE :: naerdust ! number of dust scatterers

! initializations
      tau(1:ngrid,1:naerkind)=0

! identify tracers

      !! AS: firstcall OK absolute
      IF (firstcall) THEN
        ! identify scatterers that are dust
        naerdust=0
        DO iaer=1,naerkind
          txt=name_iaer(iaer)
        ! CW17: choice tauscaling for stormdust or not
          IF ((txt(1:4).eq."dust").OR.(txt(1:5).eq."storm")
     &         .OR.(txt(1:3).eq."top")) THEN !MV19: topdust tracer
            naerdust=naerdust+1
            iaerdust(naerdust)=iaer
          ENDIF
        ENDDO
        ! identify tracers which are dust
        i=0
        DO iq=1,nq
          txt=noms(iq)
          IF (txt(1:4).eq."dust") THEN
          i=i+1
          nqdust(i)=iq
          ENDIF
        ENDDO

        IF (water.AND.activice) THEN
          i_ice=igcm_h2o_ice
          write(*,*) "aeropacity: i_ice=",i_ice
        ENDIF

c       typical profile of solsir and (1-w)^(-1):
c       --- purely for diagnostics and printing
        msolsir(1:nlayer,1:naerkind)=0
        mqextsqabs(1:nlayer,1:naerkind)=0
        WRITE(*,*) "Typical profiles of Qext(vis)/Qext(IR)"
        WRITE(*,*) "  and Qext(IR)/Qabs(IR):"
        DO iaer = 1, naerkind ! Loop on aerosol kind
          WRITE(*,*) "Aerosol # ",iaer
          DO l=1,nlayer
            DO ig=1,ngrid
              msolsir(l,iaer)=msolsir(l,iaer)+
     &              QREFvis3d(ig,l,iaer)/
     &              QREFir3d(ig,l,iaer)
              mqextsqabs(l,iaer)=mqextsqabs(l,iaer)+
     &              (1.E0-omegaREFir3d(ig,l,iaer))**(-1)
            ENDDO
            msolsir(l,iaer)=msolsir(l,iaer)/REAL(ngrid)
            mqextsqabs(l,iaer)=mqextsqabs(l,iaer)/REAL(ngrid)
          ENDDO
          WRITE(*,*) "solsir: ",msolsir(:,iaer)
          WRITE(*,*) "Qext/Qabs(IR): ",mqextsqabs(:,iaer)
        ENDDO

!       load value of tauvis from callphys.def (if given there,
!       otherwise default value read from starfi.nc file will be used)
        call getin("tauvis",tauvis)

        IF (freedust.or.rdstorm) THEN ! if rdstorm no need to held opacity constant at the first levels
          cstdustlevel = 1
        ELSE
          cstdustlevel = cstdustlevel0 !Opacity in the first levels is held constant to 
                                       !avoid unrealistic values due to constant lifting
        ENDIF


#ifndef DUSTSTORM
        firstcall=.false.
#endif

      END IF ! end of if firstcall

c     Vertical column optical depth at "odpref" Pa 
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(freedust) THEN
         tauref(:) = 0. ! tauref is computed after, instead of being forced

      ELSE IF(iaervar.eq.1) THEN 
         do ig=1, ngrid
          tauref(ig)=max(tauvis,1.e-9) ! tauvis=cste (set in callphys.def
                                       ! or read in starfi
        end do
      ELSE IF (iaervar.eq.2) THEN   ! << "Viking" Scenario>>

        tauref(1) = 0.7+.3*cos(ls+80.*pi/180.) ! like seen by VL1
        do ig=2,ngrid
          tauref(ig) = tauref(1)
        end do

      ELSE IF (iaervar.eq.3) THEN  ! << "MGS" scenario >>

        taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls-4.363)))**14
        tauS= 0.1 +(0.5-0.1)  *(cos(0.5*(ls-4.363)))**14
        tauN = 0.1
c	   if (peri_day.eq.150) then
c	     tauS=0.1
c	     tauN=0.1 +(0.5-0.1)  *(cos(0.5*(ls+pi-4.363)))**14
c	     taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls+pi-4.363)))**14
c           endif
        do ig=1,ngrid
          if (latitude(ig).ge.0) then
          ! Northern hemisphere
            tauref(ig)= tauN +
     &      (taueq-tauN)*0.5*(1+tanh((45-latitude(ig)*180./pi)*6/60))
          else
          ! Southern hemisphere
            tauref(ig)= tauS +
     &      (taueq-tauS)*0.5*(1+tanh((45+latitude(ig)*180./pi)*6/60))
          endif
        enddo ! of do ig=1,ngrid
      ELSE IF (iaervar.eq.5) THEN   ! << Escalier Scenario>>
c         tauref(1) = 0.2
c         if ((ls.ge.210.*pi/180.).and.(ls.le.330.*pi/180.))
c    &                              tauref(1) = 2.5
        tauref(1) = 2.5
        if ((ls.ge.30.*pi/180.).and.(ls.le.150.*pi/180.))
     &                              tauref(1) = .2

        do ig=2,ngrid
          tauref(ig) = tauref(1)
        end do
      ELSE IF ((iaervar.ge.6).and.(iaervar.le.8)) THEN
      ! clim, cold or warm synthetic scenarios
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE IF ((iaervar.ge.24).and.(iaervar.le.34))
     &     THEN  ! << MY... dust scenarios >>
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE IF ((iaervar.eq.4).or.
     &         ((iaervar.ge.124).and.(iaervar.le.126))) THEN
       ! "old" TES assimation dust scenario (values at 700Pa in files!)
        call read_dust_scenario(ngrid,nlayer,zday,pplev,tauref)
      ELSE
        stop 'problem with iaervar in aeropacity.F'
      ENDIF

c -----------------------------------------------------------------
c Computing the opacity in each layer
c -----------------------------------------------------------------

      DO iaer = 1, naerkind ! Loop on aerosol kind
c     --------------------------------------------
        aerkind: SELECT CASE (name_iaer(iaer))
c==================================================================
        CASE("dust_conrath") aerkind      ! Typical dust profile
c==================================================================

c       Altitude of the top of the dust layer
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        zlsconst=SIN(ls-2.76)
        if (iddist.eq.1) then
          do ig=1,ngrid
             topdust(ig)=topdustref         ! constant dust layer top
          end do

        else if (iddist.eq.2) then          ! "Viking" scenario
          do ig=1,ngrid
            ! altitude of the top of the aerosol layer (km) at Ls=2.76rad:
            ! in the Viking year scenario
            topdust0(ig)=60. -22.*sinlat(ig)**2
            topdust(ig)=topdust0(ig)+18.*zlsconst
          end do

        else if(iddist.eq.3) then         !"MGS" scenario
          do ig=1,ngrid
            topdust(ig)=60.+18.*zlsconst
     &                -(32+18*zlsconst)*sin(latitude(ig))**4
     &                 - 8*zlsconst*(sin(latitude(ig)))**5
          end do
        endif

c       Optical depth in each layer :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(iddist.ge.1) then

          expfactor=0.
          DO l=1,nlayer
            DO ig=1,ngrid
c             Typical mixing ratio profile 
              if(pplay(ig,l).gt.odpref
     $                        /(988.**(topdust(ig)/70.))) then
                zp=(odpref/pplay(ig,l))**(70./topdust(ig))
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)
              else    
                expfactor=1.e-3
              endif
c             Vertical scaling function
              aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1)) *
     &          expfactor *
     &          QREFvis3d(ig,l,iaer) / QREFvis3d(ig,1,iaer)
            ENDDO
          ENDDO

        else if(iddist.eq.0) then   
c         old dust vertical distribution function (pollack90)
          DO l=1,nlayer
             DO ig=1,ngrid
                zp=odpref/pplay(ig,l)
                aerosol(ig,l,1)= tauref(ig)/odpref *
     s           (pplev(ig,l)-pplev(ig,l+1))
     s           *max( exp(.03*(1.-max(zp,1.))) , 1.E-3 )
             ENDDO
          ENDDO
        end if

c==================================================================
        CASE("dust_doubleq") aerkind! Two-moment scheme for background dust
c        (transport of mass and number mixing ratio)
c==================================================================

          DO l=1,nlayer
            IF (l.LE.cstdustlevel) THEN
c           Opacity in the first levels is held constant to 
c             avoid unrealistic values due to constant lifting:
              DO ig=1,ngrid
	      ! OPTICAL DEPTH for the computation of tauref, 
	      ! which is to be compared with tauref_scenario
	      ! => visible wavelength
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &          ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &          pq(ig,cstdustlevel,igcm_dust_mass) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ! DENSITY SCALED OPACITY :
	      ! GCM output to be compared with observations
	      ! => infrared wavelength
                dsodust(ig,l) =
     &          (  0.75 * QREFir3d(ig,cstdustlevel,iaer) /
     &          ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &          pq(ig,cstdustlevel,igcm_dust_mass)
              ENDDO
            ELSE
              DO ig=1,ngrid
              ! OPTICAL DEPTH for the computation of tauref, 
	      ! which is to be compared with tauref_scenario
	      ! => visible wavelength 
	        aerosol(ig,l,iaer) =
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_mass) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ! DENSITY SCALED OPACITY :
	      ! GCM output to be compared with observations
	      ! => infrared wavelength
                dsodust(ig,l) =
     &          (  0.75 * QREFir3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_mass)
              ENDDO
            ENDIF
          ENDDO

c==================================================================
        CASE("dust_submicron") aerkind   ! Small dust population
c==================================================================

          DO l=1,nlayer
            IF (l.LE.cstdustlevel) THEN
c           Opacity in the first levels is held constant to 
c             avoid unrealistic values due to constant lifting:
              DO ig=1,ngrid
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &          ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &          pq(ig,cstdustlevel,igcm_dust_submicron) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ELSE
              DO ig=1,ngrid
                aerosol(ig,l,iaer) = 
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_submicron) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
              ENDDO
            ENDIF
          ENDDO

c==================================================================
        CASE("h2o_ice") aerkind             ! Water ice crystals
c==================================================================

c       1. Initialization
        aerosol(1:ngrid,1:nlayer,iaer) = 0.
        taucloudvis(1:ngrid) = 0.
        taucloudtes(1:ngrid) = 0.
c       2. Opacity calculation
        ! NO CLOUDS
        IF (clearsky) THEN
            aerosol(1:ngrid,1:nlayer,iaer) =1.e-9
        ! CLOUDSs
        ELSE ! else (clearsky)
          DO ig=1, ngrid
            DO l=1,nlayer
              aerosol(ig,l,iaer) = max(1E-20,
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_ice * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,i_ice) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
     &                              )
              taucloudvis(ig) = taucloudvis(ig) + aerosol(ig,l,iaer)
              taucloudtes(ig) = taucloudtes(ig) + aerosol(ig,l,iaer)*
     &          QREFir3d(ig,l,iaer) / QREFvis3d(ig,l,iaer) *
     &          ( 1.E0 - omegaREFir3d(ig,l,iaer) )
            ENDDO
          ENDDO
          ! SUB-GRID SCALE CLOUDS
          IF (CLFvarying) THEN
             DO ig=1, ngrid
                DO l=1,nlayer-1
                   CLFtot  = max(totcloudfrac(ig),0.01) 
                   aerosol(ig,l,iaer)=
     &                    aerosol(ig,l,iaer)/CLFtot
                   aerosol(ig,l,iaer) = 
     &                    max(aerosol(ig,l,iaer),1.e-9)
                ENDDO
             ENDDO
!          ELSE ! else (CLFvarying)
!             DO ig=1, ngrid
!                DO l=1,nlayer-1 ! to stop the rad tran bug
!                   CLFtot  = CLFfixval
!                   aerosol(ig,l,iaer)=
!     &                    aerosol(ig,l,iaer)/CLFtot
!                   aerosol(ig,l,iaer) = 
!     &                    max(aerosol(ig,l,iaer),1.e-9)
!                ENDDO
!             ENDDO
          ENDIF ! end (CLFvarying)             
        ENDIF ! end (clearsky)

c==================================================================
        CASE("stormdust_doubleq") aerkind ! CW17 : Two-moment scheme for 
c       stormdust  (transport of mass and number mixing ratio) 
c==================================================================
c       aerosol is calculated twice : once within the dust storm (clearatm=false)
c       and once in the part of the mesh without dust storm (clearatm=true)
        aerosol(1:ngrid,1:nlayer,iaer) = 0.
        IF (clearatm) THEN  ! considering part of the mesh without storm
          aerosol(1:ngrid,1:nlayer,iaer)=1.e-25
        ELSE  ! part of the mesh with concentred dust storm
          DO l=1,nlayer
             IF (l.LE.cstdustlevel) THEN
c          Opacity in the first levels is held constant to 
c           avoid unrealistic values due to constant lifting:
               DO ig=1,ngrid
	       ! OPTICAL DEPTH for the computation of tauref, 
	       ! which is to be compared with tauref_scenario
	       ! => visible wavelength
                 aerosol(ig,l,iaer) = 
     &           (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &           ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &           pq(ig,cstdustlevel,igcm_stormdust_mass) *
     &           ( pplev(ig,l) - pplev(ig,l+1) ) / g
               ! DENSITY SCALED OPACITY :
	       ! GCM output to be compared with observations
	       ! => infrared wavelength
                 dsords(ig,l) =
     &           (  0.75 * QREFir3d(ig,cstdustlevel,iaer) /
     &           ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &           pq(ig,cstdustlevel,igcm_stormdust_mass)
               ENDDO
             ELSE
               DO ig=1,ngrid
	       ! OPTICAL DEPTH for the computation of tauref, 
	       ! which is to be compared with tauref_scenario
	       ! => visible wavelength
                 aerosol(ig,l,iaer) =
     &           (  0.75 * QREFvis3d(ig,l,iaer) /
     &           ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &           pq(ig,l,igcm_stormdust_mass) *
     &           ( pplev(ig,l) - pplev(ig,l+1) ) / g
               ! DENSITY SCALED OPACITY :
	       ! GCM output to be compared with observations
	       ! => infrared wavelength
                 dsords(ig,l) =
     &           (  0.75 * QREFir3d(ig,l,iaer) /
     &           ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &           pq(ig,l,igcm_stormdust_mass)
               ENDDO
             ENDIF
          ENDDO
        ENDIF
c==================================================================
        CASE("topdust_doubleq") aerkind ! MV18 : Two-moment scheme for 
c       topdust  (transport of mass and number mixing ratio) 
c==================================================================
c       aerosol is calculated twice : once "above" the sub-grid mountain (nohmons=false)
c       and once in the part of the mesh without the sub-grid mountain (nohmons=true)
        aerosol(1:ngrid,1:nlayer,iaer) = 0.
        IF (nohmons) THEN  ! considering part of the mesh without storm
          aerosol(1:ngrid,1:nlayer,iaer)=1.e-25
        ELSE  ! part of the mesh with concentred dust storm
          DO l=1,nlayer
             IF (l.LE.cstdustlevel) THEN
c          Opacity in the first levels is held constant to 
c           avoid unrealistic values due to constant lifting:
               DO ig=1,ngrid
                  aerosol(ig,l,iaer) = 
     &           (  0.75 * QREFvis3d(ig,cstdustlevel,iaer) /
     &           ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &           pq(ig,cstdustlevel,igcm_topdust_mass) *
     &           ( pplev(ig,l) - pplev(ig,l+1) ) / g
               ! DENSITY SCALED OPACITY :
	       ! GCM output to be compared with observations
	       ! => infrared wavelength
                 dsotop(ig,l) =
     &           (  0.75 * QREFir3d(ig,cstdustlevel,iaer) /
     &           ( rho_dust * reffrad(ig,cstdustlevel,iaer) )  ) *
     &           pq(ig,cstdustlevel,igcm_topdust_mass)
               ENDDO
             ELSE
               DO ig=1,ngrid
               ! OPTICAL DEPTH for the computation of tauref, 
	       ! which is to be compared with tauref_scenario
	       ! => visible wavelength
	         aerosol(ig,l,iaer) =
     &           (  0.75 * QREFvis3d(ig,l,iaer) /
     &           ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &           pq(ig,l,igcm_topdust_mass) *
     &           ( pplev(ig,l) - pplev(ig,l+1) ) / g
               ! DENSITY SCALED OPACITY :
	       ! GCM output to be compared with observations
	       ! => infrared wavelength
                 dsotop(ig,l) =
     &           (  0.75 * QREFir3d(ig,l,iaer) /
     &           ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &           pq(ig,l,igcm_topdust_mass)
               ENDDO
             ENDIF
          ENDDO
        ENDIF
c==================================================================
        END SELECT aerkind
c     -----------------------------------
      ENDDO ! iaer (loop on aerosol kind)

c -----------------------------------------------------------------
c Rescaling each layer to reproduce the choosen (or assimilated)
c   dust extinction opacity at visible reference wavelength, which
c   is originally scaled to an equivalent odpref Pa pressure surface.
c -----------------------------------------------------------------


#ifdef DUSTSTORM
c -----------------------------------------------------------------
! Calculate reference opacity without perturbation
c -----------------------------------------------------------------
      IF (firstcall) THEN
        DO iaer=1,naerdust
          DO l=1,nlayer
            DO ig=1,ngrid
              tauref(ig) = tauref(ig) +
     &                    aerosol(ig,l,iaerdust(iaer))
            ENDDO
          ENDDO
        ENDDO
        tauref(:) = tauref(:) * odpref / pplev(:,1)

c--------------------------------------------------
c Get parameters of the opacity perturbation
c--------------------------------------------------
        iaer=1  ! just change dust

        write(*,*) "Add a local storm ?"
        localstorm=.true. ! default value
        call getin("localstorm",localstorm)
        write(*,*) " localstorm = ",localstorm

        IF (localstorm) THEN
          WRITE(*,*) "********************"
          WRITE(*,*) "ADDING A LOCAL STORM"
          WRITE(*,*) "********************"

          write(*,*) "ref opacity of local dust storm"
              taulocref = 4.25 ! default value
              call getin("taulocref",taulocref)
              write(*,*) " taulocref = ",taulocref

          write(*,*) "target altitude of local storm (km)"
              ztoploc = 10.0 ! default value
              call getin("ztoploc",ztoploc)
              write(*,*) " ztoploc = ",ztoploc

          write(*,*) "radius of dust storm (degree)"
              radloc = 0.5 ! default value
              call getin("radloc",radloc)
              write(*,*) " radloc = ",radloc

          write(*,*) "center longitude of storm (deg)"
              lonloc = 25.0 ! default value
              call getin("lonloc",lonloc)
              write(*,*) " lonloc = ",lonloc

          write(*,*) "center latitude of storm (deg)"
              latloc = -2.5 ! default value
              call getin("latloc",latloc)
              write(*,*) " latloc = ",latloc
        
          write(*,*) "reff storm (mic) 0. for background"
              reffstorm = 0.0 ! default value
              call getin("reffstorm",reffstorm)
              write(*,*) " reffstorm = ",reffstorm

!! LOOP: modify opacity
      DO ig=1,ngrid

      !! distance to the center:
      ray(ig)=SQRT((latitude(ig)*180./pi-latloc)**2 +
     &          (longitude(ig)*180./pi -lonloc)**2)

      !! transition factor for storm
      !! factor is hardcoded -- increase it to steepen
      yeah = (TANH(2.+(radloc-ray(ig))*10.)+1.)/2.

      !! new opacity field
      !! -- add an opacity set to taulocref
      !! -- the additional reference opacity will
      !!      thus be taulocref*odpref/pplev
      tauuser(ig)=max( tauref(ig) * pplev(ig,1) /odpref , 
     &          taulocref * yeah )

      !! compute l_top
          DO l=1,nlayer
            zalt(ig,l) = LOG( pplev(ig,1)/pplev(ig,l) )
     &                      / g / 44.01
     &                    * 8.31 * 210.
                IF (     (ztoploc .lt. zalt(ig,l)  )
     &          .and. (ztoploc .gt. zalt(ig,l-1)) ) l_top=l-1
          ENDDO

     !! change reffrad if ever needed
      IF (reffstorm .gt. 0.) THEN
          DO l=1,nlayer
             IF (l .lt. l_top+1) THEN
                reffrad(ig,l,iaer) = max( reffrad(ig,l,iaer), reffstorm
     &          * 1.e-6 * yeah )
             ENDIF
          ENDDO
      ENDIF

      ENDDO
!! END LOOP

      !! compute perturbation in each layer (equation 8 in Spiga et al. JGR 2013)
      DO ig=1,ngrid
          int_factor(ig)=0.
          DO l=1,nlayer
             IF (l .lt. l_top+1) THEN
                      int_factor(ig) =
     &                int_factor(ig) +
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
             ENDIF
          ENDDO
          DO l=1, nlayer
          !! Mass mixing ratio perturbation due to local dust storm in each layer
          more_dust(ig,l,1)=
     &                     (tauuser(ig)-(tauref(ig)
     &                      * pplev(ig,1) /odpref)) /
     &                      int_factor(ig)
          more_dust(ig,l,2)=
     &                     (tauuser(ig)-(tauref(ig) *
     &                     pplev(ig,1) /odpref))
     &                      / int_factor(ig) *
     &                     ((ref_r0/reffrad(ig,l,iaer))**3)
     &                      * r3n_q 
          ENDDO
      ENDDO

      !! quantity of dust for each layer with the addition of the perturbation
      DO l=1, l_top
          pq(:,l,igcm_dust_mass)= pq(:,l,igcm_dust_mass)
     .            + more_dust(:,l,1)
          pq(:,l,igcm_dust_number)= pq(:,l,igcm_dust_number)
     .            + more_dust(:,l,2)
      ENDDO
      ENDIF !! IF (localstorm)
      tauref(:)=0.
      ENDIF !! IF (firstcall)
#endif

      IF (freedust) THEN
          tauscaling(:) = 1.
c        opacity obtained with stormdust
        IF (rdstorm) THEN
           taustormdusttmp(1:ngrid)=0.
           DO l=1,nlayer
             DO ig=1,ngrid
                taustormdusttmp(ig) = taustormdusttmp(ig)+ 
     &            aerosol(ig,l,iaerdust(2))
             ENDDO
           ENDDO
           !opacity obtained with background dust only
           taubackdusttmp(1:ngrid)=0.  
           DO l=1,nlayer
             DO ig=1,ngrid
                taubackdusttmp(ig) = taubackdusttmp(ig)+
     &           aerosol(ig,l,iaerdust(1))
             ENDDO
           ENDDO
        ENDIF !rdsstorm
      ELSE
c       Temporary scaling factor
        taudusttmp(1:ngrid)=0.
        DO iaer=1,naerdust
          DO l=1,nlayer
            DO ig=1,ngrid
c             Scaling factor
              taudusttmp(ig) = taudusttmp(ig) + 
     &                         aerosol(ig,l,iaerdust(iaer))
            ENDDO
          ENDDO
        ENDDO

c       Saved scaling factor
        DO ig=1,ngrid
            tauscaling(ig) = tauref(ig) * 
     &                       pplev(ig,1) / odpref / taudusttmp(ig)
        ENDDO

      ENDIF ! IF (freedust)

c     Opacity computation
      DO iaer=1,naerdust
        DO l=1,nlayer
          DO ig=1,ngrid
            aerosol(ig,l,iaerdust(iaer)) = max(1E-20,
     &                aerosol(ig,l,iaerdust(iaer))* tauscaling(ig))
          ENDDO
        ENDDO
      ENDDO

      IF (freedust) THEN
        ! tauref has been initialized to 0 before.
        DO iaer=1,naerdust
          DO l=1,nlayer
            DO ig=1,ngrid
#ifdef DUSTSTORM
      !! recalculate opacity because storm perturbation has been added
      IF (firstcall) THEN
              aerosol(ig,l,iaer) =
     &          (  0.75 * QREFvis3d(ig,l,iaer) /
     &          ( rho_dust * reffrad(ig,l,iaer) )  ) *
     &          pq(ig,l,igcm_dust_mass) *
     &          ( pplev(ig,l) - pplev(ig,l+1) ) / g
      ENDIF
#endif
!              tauref(ig) = tauref(ig) + 
!     &                    aerosol(ig,l,iaerdust(iaer))
c      MV19: tauref must ALWAYS contain the opacity of all dust tracers
       IF (name_iaer(iaerdust(iaer)).eq."dust_doubleq") THEN
              tauref(ig) = tauref(ig) +
     &  (  0.75 * QREFvis3d(ig,l,iaerdust(iaer)) /
     &  ( rho_dust * reffrad(ig,l,iaerdust(iaer)) )  ) *
     &  pq(ig,l,igcm_dust_mass) *
     &  ( pplev(ig,l) - pplev(ig,l+1) ) / g
       ELSE IF (name_iaer(iaerdust(iaer)).eq."stormdust_doubleq") THEN
              tauref(ig) = tauref(ig) +
     &  (  0.75 * QREFvis3d(ig,l,iaerdust(iaer)) /
     &  ( rho_dust * reffrad(ig,l,iaerdust(iaer)) )  ) *
     &  pq(ig,l,igcm_stormdust_mass) *
     &  ( pplev(ig,l) - pplev(ig,l+1) ) / g 
       ELSE IF (name_iaer(iaerdust(iaer)).eq."topdust_doubleq") THEN
              tauref(ig) = tauref(ig) +
     &  (  0.75 * QREFvis3d(ig,l,iaerdust(iaer)) /
     &  ( rho_dust * reffrad(ig,l,iaerdust(iaer)) )  ) *
     &  pq(ig,l,igcm_topdust_mass) *
     &  ( pplev(ig,l) - pplev(ig,l+1) ) / g
       ENDIF

            ENDDO
          ENDDO
        ENDDO
        tauref(:) = tauref(:) * odpref / pplev(:,1)
      ENDIF

c -----------------------------------------------------------------
c Column integrated visible optical depth in each point
c -----------------------------------------------------------------
      DO iaer=1,naerkind
        do l=1,nlayer
           do ig=1,ngrid
             tau(ig,iaer) = tau(ig,iaer) + aerosol(ig,l,iaer)
           end do
        end do
      ENDDO

c     for diagnostics: opacity for all dust scatterers stormdust included 
      taualldust(1:ngrid)=0.
      DO iaer=1,naerdust
        DO l=1,nlayer
          DO ig=1,ngrid
            taualldust(ig) = taualldust(ig) + 
     &                         aerosol(ig,l,iaerdust(iaer))
          ENDDO
        ENDDO
      ENDDO
     
      IF (rdstorm) THEN
 
c     for diagnostics: opacity for dust in background only  
       taudust(1:ngrid)=0.
        DO l=1,nlayer
         DO ig=1,ngrid
           taudust(ig) = taudust(ig) + 
     &                       aerosol(ig,l,iaer_dust_doubleq)
         ENDDO
        ENDDO

c     for diagnostics: opacity for dust in storm only  
       taustormdust(1:ngrid)=0.
        DO l=1,nlayer
         DO ig=1,ngrid
           taustormdust(ig) = taustormdust(ig) + 
     &                       aerosol(ig,l,iaer_stormdust_doubleq)
         ENDDO
        ENDDO
 
      ENDIF
      

#ifdef DUSTSTORM
      IF (firstcall) THEN
        firstcall=.false.
      ENDIF
#endif


c -----------------------------------------------------------------
c aerosol/X for stormdust to prepare calculation of radiative transfer
c -----------------------------------------------------------------
      IF (rdstorm) THEN
         DO l=1,nlayer
           DO ig=1,ngrid
               ! stormdust: opacity relative to the storm fraction (stormdust/x)
               aerosol(ig,l,iaer_stormdust_doubleq) = 
     &           aerosol(ig,l,iaer_stormdust_doubleq)/totstormfract(ig)
           ENDDO
         ENDDO
      ENDIF 

c -----------------------------------------------------------------
c aerosol/X for topdust to prepare calculation of radiative transfer
c -----------------------------------------------------------------
      IF (slpwind) THEN
        DO ig=1,ngrid
          IF (alpha_hmons(ig) .gt. 0.) THEN
            DO l=1,nlayer
             ! topdust: opacity relative to the storm fraction (topdust/x)
              aerosol(ig,l,iaer_topdust_doubleq) =
     &        aerosol(ig,l,iaer_topdust_doubleq)/alpha_hmons(ig)
            ENDDO
          ENDIF
        ENDDO 
      ENDIF

      END SUBROUTINE aeropacity
      
      END MODULE aeropacity_mod
