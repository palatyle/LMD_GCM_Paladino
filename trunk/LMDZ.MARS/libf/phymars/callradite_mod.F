      MODULE callradite_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE callradite(icount,ngrid,nlayer,nq,zday,ls,pq,albedo,
     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
     $     dtlw,dtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,
     $     fluxtop_sw,tauref,tau,aerosol,dsodust,tauscaling,
     $     taucloudtes,rdust,rice,nuice,co2ice,rstormdust,rtopdust,
     $     totstormfract,clearatm,dsords,dsotop,alpha_hmons,nohmons,
     $     clearsky,totcloudfrac)

      use aeropacity_mod, only: aeropacity
      use updatereffrad_mod, only: updatereffrad
      use dimradmars_mod, only: ndomainsz, nflev, nsun, nir
      use dimradmars_mod, only: naerkind, name_iaer,
     &            iaer_dust_conrath,iaer_dust_doubleq,
     &            iaer_dust_submicron,iaer_h2o_ice,
     &            iaer_stormdust_doubleq,iaer_topdust_doubleq
      use yomlw_h, only: gcp, nlaylte
      use comcstfi_h, only: g,cpp
      use time_phylmdz_mod, only: daysec
      use lwmain_mod, only: lwmain
      use swmain_mod, only: swmain
      IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Subroutine designed to call the main canonic
c   radiative transfer subroutine "lwmain" et "swmain"
c   to compute radiative heating and cooling rate and
c   radiative fluxes to the surface.
c
c   These calculations are only valid on the part of the atmosphere
c   where Local Thermal Equilibrium (LTE) is verified. In practice
c   the calculations are only performed for the first "nlaylte"
c   parameters (nlaylte is calculated by subroutine "nlthermeq"
c   and stored in module "yomlw_h").
c
c   The purpose of this subroutine is to:
c      1) Make some initial calculation at first call
c      2) Split the calculation in several sub-grid
c        ("sub-domain") to save memory and
c        be able run on a workstation at high resolution
c        The sub-grid size is defined in dimradmars_mod
c      3) Compute the 3D scattering parameters depending on the
c        size distribution of the different tracers (added by JBM)
c      4) call "lwmain" and "swmain"
c
c
c   authors:   
c   ------
c   Francois Forget / Christophe Hourdin / J.-B. Madeleine (2009)
c
c
c   3D scattering scheme user's guide (J.-B. Madeleine)
c   ---------------------------------
c
c   This routine has been modified to take into account 3D, time
c   dependent scattering properties of the aerosols.
c---- The look-up tables that contain the scattering parameters
c   of a given tracer, for different sizes, are read by SUAER.F90.
c   The names of the corresponding ASCII files have to be set in
c   this subroutine (file_id variable), and files must be in the
c   directory specified in datafile_mod. Please make sure that the
c   ASCII files are correctly written, and that the range
c   of particle sizes is consistent with what you would expect.
c---- SUAER.F90 is in charge of reading the ASCII files and averaging
c   the scattering parameters in each GCM channel, using the three last
c   equations of Forget et al. 1998 (GRL 25, No.7, p.1105-1108).
c---- These look-up tables, loaded during the firstcall, are then
c   constantly used by the subroutine "aeroptproperties.F" to compute,
c   online, the 3D scattering parameters, based on the size distribution
c   (reffrad and nueffrad) of the different tracers, in each grid box.
c   These 3D size distributions are loaded by the "updatereffrad.F"
c   subroutine. A log-normal distribution is then assumed in
c   "aeroptproperties.F", along with a Gauss-Legendre integration.
c---- The optical depth at the visible reference wavelength (set in
c   SUAER.F90, after the file_id variable) is then computed by
c   the subroutine "aeropacity.F", by using the size and spatial
c   distribution of the corresponding tracer. This connection has to
c   be implemented in "aeropacity.F" when adding a new tracer. To do so,
c   one can use equation 2 of Forget et al. 1998 (Icarus 131, p.302-316).
c---- The resulting variables "aerosol", "QVISsQREF3d", "omegaVIS3d" and
c   "gVIS3d" (same in the infrared) are finally used by lwmain.F and 
c   swmain.F to solve the radiative transfer equation.
c
c   changes:
c   -------
c
c   > SRL 7/2000
c   
c   This version has been modified to only calculate radiative tendencies
c   over layers 1..NFLEV (set in dimradmars_mod).  Returns zero for higher
c   layers, if any.
c   In other routines, nlayer -> nflev.
c   Routines affected: lwflux, lwi, lwmain, lwxb, lwxd, lwxn.
c
c   > J.-B. Madeleine 10W12
c   This version uses the variable's splitting, which can be usefull
c     when performing very high resolution simulation like LES.
c
c   ----------
c   Here, solar band#1 is spectral interval between "long1vis" and "long2vis"
c   set in dimradmars_mod 
c   Here, solar band#2 is spectral interval between "long2vis" and "long3vis"
c   set in dimradmars_mod
c
c   input:
c   ----- 
c   icount                counter of call to subroutine physic by gcm
c   ngrid                 number of gridpoint of horizontal grid
c   nlayer                Number of layer
c   nq                    Number of tracer
c   ls                    Solar longitude (Ls) , radian
c   zday                  Date (time since Ls=0, in martian days)
c   pq(ngrid,nlayer,nq)   Advected fields
c
c   albedo (ngrid,2)      hemispheric surface albedo
c                         albedo (i,1) : mean albedo for solar band#1 
c                                        (see below)
c                         albedo (i,2) : mean albedo for solar band#2
c                                        (see below)
c   emis                  Thermal IR surface emissivity (no unit)
c   mu0(ngrid)           cos of solar zenith angle
c                           (=1 when sun at zenith)
c   pplay(ngrid,nlayer)    pressure (Pa) in the middle of each layer
c   pplev(ngrid,nlayer+1)  pressure (Pa) at boundaries of each layer
c   pt(ngrid,nlayer)       atmospheric temperature in each layer (K)
c   tsurf(ngrid)           surface temperature (K)
c   fract(ngrid)         day fraction of the time interval 
c                          =1 during the full day ; =0 during the night
c   declin                 latitude of subsolar point
c   dist_sol               sun-Mars distance (AU)
c   igout                  coordinate of analysed point for debugging
c   reffrad(ngrid,nlayer,naerkind)  Aerosol effective radius
c   nueffrad(ngrid,nlayer,naerkind) Aerosol effective variance

c
c  output:
c  -------
c dtlw (ngrid,nlayer)       longwave (IR) heating rate (K/s)
c dtsw(ngrid,nlayer)        shortwave (Solar) heating rate (K/s)
c fluxsurf_lw(ngrid)        surface downward flux tota LW (thermal IR) (W.m-2)
c fluxsurf_sw(ngrid,1)      surface downward flux SW for solar band#1 (W.m-2)
c fluxsurf_sw(ngrid,2)      surface downward flux SW for solar band#2 (W.m-2)
c
c fluxtop_lw(ngrid)         outgoing upward flux tota LW (thermal IR) (W.m-2)
c fluxtop_sw(ngrid,1)       outgoing upward flux SW for solar band#1 (W.m-2)
c fluxtop_sw(ngrid,2)       outgoing upward flux SW for solar band#2 (W.m-2)

c   tauref       Prescribed mean column optical depth at 610 Pa 
c   tau          Column total visible dust optical depth at each point
c   aerosol(ngrid,nlayer,naerkind)    aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars_h , in each layer, for one of
c                         the "naerkind" kind of aerosol optical
c                         properties.
c=======================================================================
c
c    Declarations :
c    -------------
c
      include "callkeys.h"

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER,INTENT(IN) :: icount        
      INTEGER,INTENT(IN) :: ngrid,nlayer,nq 
      INTEGER,INTENT(IN) :: igout

      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq)
      REAL,INTENT(INOUT) :: tauscaling(ngrid) ! Conversion factor for 
                               ! qdust and Ndust
      REAL,INTENT(IN) :: albedo(ngrid,2),emis(ngrid)
      REAL,INTENT(IN) :: ls,zday

      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL,INTENT(IN) :: pt(ngrid,nlayer)
      REAL,INTENT(IN) :: tsurf(ngrid)
      REAL,INTENT(IN) :: dist_sol,mu0(ngrid),fract(ngrid)
      REAL,INTENT(OUT) :: dtlw(ngrid,nlayer),dtsw(ngrid,nlayer)
      REAL,INTENT(OUT) :: fluxsurf_lw(ngrid), fluxtop_lw(ngrid)
      REAL,INTENT(OUT) :: fluxsurf_sw(ngrid,2), fluxtop_sw(ngrid,2)

      REAL,INTENT(OUT) :: tauref(ngrid), tau(ngrid,naerkind)
      REAL,INTENT(OUT) :: taucloudtes(ngrid)! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)
      REAL,INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind)
      REAL,INTENT(INOUT) :: dsodust(ngrid,nlayer)
      REAL,INTENT(OUT) :: rdust(ngrid,nlayer)  ! Dust geometric mean radius (m)
      REAL,INTENT(OUT) :: rice(ngrid,nlayer)   ! Ice geometric mean radius (m)
      REAL,INTENT(OUT) :: nuice(ngrid,nlayer)  ! Estimated effective variance
      REAL,INTENT(IN) :: co2ice(ngrid)           ! co2 ice surface layer (kg.m-2)

c     rocket dust storm
      LOGICAL,INTENT(IN) :: clearatm ! true for background dust
      REAL,INTENT(IN) :: totstormfract(ngrid) ! dust storm mesh fraction
      REAL,INTENT(OUT) :: rstormdust(ngrid,nlayer)  ! Storm dust geometric mean radius (m)
      REAL,INTENT(INOUT) :: dsords(ngrid,nlayer) ! density scaled opacity for rocket dust storm dust
      
c     entrainment by slope wind
      LOGICAL, INTENT(IN) :: nohmons ! true for background dust 
      REAL, INTENT(IN) :: alpha_hmons(ngrid) ! sub-grid scale topography mesh fraction
      REAL,INTENT(OUT) :: rtopdust(ngrid,nlayer)  ! Topdust geometric mean radius (m)
      REAL,INTENT(INOUT) :: dsotop(ngrid,nlayer) ! density scaled opacity for topmons dust
      
c     sub-grid scale water ice clouds
      LOGICAL,INTENT(IN) :: clearsky
      REAL,INTENT(IN) :: totcloudfrac(ngrid)

c
c    Local variables :
c    -----------------

      INTEGER j,l,ig,n,ich
      INTEGER aer_count,iaer
      INTEGER jd,ig0,nd

      real  cste_mars ! solar constant on Mars (Wm-2)
      REAL ptlev(ngrid,nlayer+1)

      INTEGER :: ndomain

c     Thermal IR net radiative budget (W m-2)

      real znetrad(ndomainsz,nflev)

      real zfluxd_sw(ndomainsz,nflev+1,2)
      real zfluxu_sw(ndomainsz,nflev+1,2)

      REAL zplev(ndomainsz,nflev+1)
      REAL zztlev(ndomainsz,nflev+1)
      REAL zplay(ndomainsz,nflev)
      REAL zt(ndomainsz,nflev)
      REAL zaerosol(ndomainsz,nflev,naerkind)
      REAL zalbedo(ndomainsz,2)
      REAL zdp(ndomainsz,nflev)
      REAL zdt0(ndomainsz)

      REAL zzdtlw(ndomainsz,nflev)
      REAL zzdtsw(ndomainsz,nflev)
      REAL zzflux(ndomainsz,6)
      real zrmuz

      REAL :: zQVISsQREF3d(ndomainsz,nflev,nsun,naerkind)
      REAL :: zomegaVIS3d(ndomainsz,nflev,nsun,naerkind)
      REAL :: zgVIS3d(ndomainsz,nflev,nsun,naerkind)

      REAL :: zQIRsQREF3d(ndomainsz,nflev,nir,naerkind)
      REAL :: zomegaIR3d(ndomainsz,nflev,nir,naerkind)
      REAL :: zgIR3d(ndomainsz,nflev,nir,naerkind)

c     Aerosol size distribution
      REAL :: reffrad(ngrid,nlayer,naerkind)
      REAL :: nueffrad(ngrid,nlayer,naerkind)
c     Aerosol optical properties
      REAL :: QVISsQREF3d(ngrid,nlayer,nsun,naerkind)
      REAL :: omegaVIS3d(ngrid,nlayer,nsun,naerkind)
      REAL :: gVIS3d(ngrid,nlayer,nsun,naerkind)

      REAL :: QIRsQREF3d(ngrid,nlayer,nir,naerkind)
      REAL :: omegaIR3d(ngrid,nlayer,nir,naerkind)
      REAL :: gIR3d(ngrid,nlayer,nir,naerkind)

      REAL :: QREFvis3d(ngrid,nlayer,naerkind)
      ! QREFvis3d : Extinction efficiency at the VISible reference wavelength
      REAL :: QREFir3d(ngrid,nlayer,naerkind)
      ! QREFir3d : Extinction efficiency at the InfraRed reference wavelength

      REAL :: omegaREFvis3d(ngrid,nlayer,naerkind)
      REAL :: omegaREFir3d(ngrid,nlayer,naerkind)

c   local saved variables
c   ---------------------

      real zco2   ! volume fraction of CO2 in Mars atmosphere
      DATA zco2/0.95/
      SAVE zco2

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

c----------------------------------------------------------------------

c     Initialisation
c     --------------

! compute ndomain
! AS: moved out of firstcall to allow nesting+evoluting domain
! ------------------------------------------------------------
      ndomain= (ngrid-1) / ndomainsz + 1

      IF (firstcall) THEN

         write(*,*) 'Splitting radiative calculations: ',
     $              ' ngrid,ndomainsz,ndomain',
     $                ngrid,ndomainsz,ndomain

c        Assign a number to the different scatterers
c        -------------------------------------------

         iaer_dust_conrath=0
         iaer_dust_doubleq=0
         iaer_dust_submicron=0
         iaer_h2o_ice=0
         iaer_stormdust_doubleq=0
         iaer_topdust_doubleq=0

         aer_count=0
         if (.NOT.active) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."dust_conrath") then
               iaer_dust_conrath = iaer
               aer_count = aer_count + 1
             endif
           enddo
         endif
         if (doubleq.AND.active) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."dust_doubleq") then
               iaer_dust_doubleq = iaer
               aer_count = aer_count + 1
             endif
           enddo
         endif
         if (submicron.AND.active) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."dust_submicron") then
               iaer_dust_submicron = iaer
               aer_count = aer_count + 1
             endif
           enddo
         endif
         if (water.AND.activice) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."h2o_ice") then
               iaer_h2o_ice = iaer
               aer_count = aer_count + 1
             endif
           enddo
         endif
         if (rdstorm.AND.active) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."stormdust_doubleq") then
               iaer_stormdust_doubleq = iaer
               aer_count = aer_count + 1
             endif
           enddo
         end if
         if (slpwind.AND.active) then
           do iaer=1,naerkind
             if (name_iaer(iaer).eq."topdust_doubleq") then
               iaer_topdust_doubleq = iaer
               aer_count = aer_count + 1
             endif
           enddo
         end if

c        Check that we identified all tracers:
         if (aer_count.ne.naerkind) then
           write(*,*) "callradite: found only ",aer_count," scatterers"
           write(*,*) "               expected ",naerkind
           write(*,*) "please make sure that the number of"
           write(*,*) "scatterers in scatterers.h, the names"
           write(*,*) "in callradite.F, and the flags in"
           write(*,*) "callphys.def are all consistent!"
           do iaer=1,naerkind
             write(*,*)'      ',iaer,' ',trim(name_iaer(iaer))
           enddo
           stop
         else
           write(*,*) "callradite: found all scatterers, namely:"
           do iaer=1,naerkind
             write(*,*)'      ',iaer,' ',trim(name_iaer(iaer))
           enddo
         endif
c        -------------------------------------------

         gcp = g/cpp

c        Loading the optical properties in external look-up tables:
         CALL SUAER
!         CALL SULW ! this step is now done in ini_yomlw_h

         if (ngrid .EQ. 1) then
           if (ndomainsz .NE. 1) then
             print*
             print*,'ATTENTION !!!'
             print*,'pour tourner en 1D, '
             print*,'fixer ndomainsz=1 dans phymars/dimradmars_h'
             print*
             call exit(1)
           endif
         endif

         firstcall=.false.
      END IF

c     Computing aerosol optical properties and opacity
c     ------------------------------------------------

c     Updating aerosol size distributions:
      CALL updatereffrad(ngrid,nlayer,
     &                rdust,rstormdust,rtopdust,rice,nuice,
     &                reffrad,nueffrad,
     &                pq,tauscaling,tau,pplay)

c     Computing 3D scattering parameters:
      CALL aeroptproperties(ngrid,nlayer,reffrad,nueffrad,
     &                      QVISsQREF3d,omegaVIS3d,gVIS3d,
     &                      QIRsQREF3d,omegaIR3d,gIR3d,
     &                      QREFvis3d,QREFir3d,
     &                      omegaREFvis3d,omegaREFir3d)

c     Computing aerosol optical depth in each layer:
      CALL aeropacity(ngrid,nlayer,nq,zday,pplay,pplev,ls,
     &    pq,tauscaling,tauref,tau,taucloudtes,aerosol,dsodust,reffrad,
     &    QREFvis3d,QREFir3d,omegaREFir3d,
     &    totstormfract,clearatm,dsords,dsotop,
     &    alpha_hmons,nohmons,
     &    clearsky,totcloudfrac)

c     Starting loop on sub-domain
c     ----------------------------

      DO jd=1,ndomain
        ig0=(jd-1)*ndomainsz
        if (jd.eq.ndomain) then
         nd=ngrid-ig0
        else
         nd=ndomainsz
        endif

c       Spliting input variable in sub-domain input variables
c       ---------------------------------------------------

        do l=1,nlaylte
         do ig = 1,nd
           do iaer = 1, naerkind
             do ich = 1, nsun
               zQVISsQREF3d(ig,l,ich,iaer) = 
     &                           QVISsQREF3d(ig0+ig,l,ich,iaer)
               zomegaVIS3d(ig,l,ich,iaer) = 
     &                           omegaVIS3d(ig0+ig,l,ich,iaer)
               zgVIS3d(ig,l,ich,iaer) = 
     &                           gVIS3d(ig0+ig,l,ich,iaer)
             enddo
             do ich = 1, nir
               zQIRsQREF3d(ig,l,ich,iaer) = 
     &                           QIRsQREF3d(ig0+ig,l,ich,iaer)
               zomegaIR3d(ig,l,ich,iaer) = 
     &                           omegaIR3d(ig0+ig,l,ich,iaer)
               zgIR3d(ig,l,ich,iaer) = 
     &                           gIR3d(ig0+ig,l,ich,iaer)
             enddo
           enddo
         enddo
        enddo

        do l=1,nlaylte+1
         do ig = 1,nd
          zplev(ig,l) = pplev(ig0+ig,l)
         enddo
        enddo

        do l=1,nlaylte
         do ig = 1,nd
          zplay(ig,l) = pplay(ig0+ig,l)
          zt(ig,l) = pt(ig0+ig,l)
c         Thickness of each layer (Pa) :
          zdp(ig,l)= pplev(ig0+ig,l) - pplev(ig0+ig,l+1)
         enddo
        enddo

        do n=1,naerkind
          do l=1,nlaylte
            do ig=1,nd
              zaerosol(ig,l,n) = aerosol(ig0+ig,l,n)
            enddo
          enddo
        enddo

        do j=1,2
          do ig = 1,nd
           zalbedo(ig,j) = albedo(ig0+ig,j)
          enddo
        enddo

c       Intermediate  levels: (computing tlev)
c       ---------------------------------------
c       Extrapolation for the air temperature above the surface
        DO ig=1,nd
              zztlev(ig,1)=zt(ig,1)+
     s        (zplev(ig,1)-zplay(ig,1))*
     s        (zt(ig,1)-zt(ig,2))/(zplay(ig,1)-zplay(ig,2))

              zdt0(ig) = tsurf(ig0+ig) - zztlev(ig,1)
        ENDDO

        DO l=2,nlaylte
         DO ig=1, nd
               zztlev(ig,l)=0.5*(zt(ig,l-1)+zt(ig,l)) 
         ENDDO
        ENDDO

        DO ig=1, nd
           zztlev(ig,nlaylte+1)=zt(ig,nlaylte)
        ENDDO


c       Longwave ("lw") radiative transfer (= thermal infrared)
c       -------------------------------------------------------
        call lwmain (ig0,icount,nd,nflev
     .        ,zdp,zdt0,emis(ig0+1),zplev,zztlev,zt
     .        ,zaerosol,zzdtlw
     .        ,fluxsurf_lw(ig0+1),fluxtop_lw(ig0+1)
     .        ,znetrad
     &        ,zQIRsQREF3d,zomegaIR3d,zgIR3d
     &        ,co2ice(ig0+1))

c       Shortwave ("sw") radiative transfer (= solar radiation)
c       -------------------------------------------------------
c          Mars solar constant (W m-2)
c          1370 W.m-2 is the solar constant at 1 AU.
           cste_mars=1370./(dist_sol*dist_sol)

           call swmain ( nd, nflev,
     S     cste_mars, zalbedo,
     S     mu0(ig0+1), zdp, zplev, zaerosol, fract(ig0+1),
     S     zzdtsw, zfluxd_sw, zfluxu_sw,
     &     zQVISsQREF3d,zomegaVIS3d,zgVIS3d)

c       ------------------------------------------------------------
c       Un-spliting output variable from sub-domain input variables
c       ------------------------------------------------------------

        do l=1,nlaylte
         do ig = 1,nd
          dtlw(ig0+ig,l) = zzdtlw(ig,l)
          dtsw(ig0+ig,l) = zzdtsw(ig,l)
         enddo
        enddo

        do l=1,nlaylte+1
         do ig = 1,nd
          ptlev(ig0+ig,l) = zztlev(ig,l)
         enddo
        enddo

        do ig = 1,nd
          fluxsurf_sw(ig0+ig,1) = zfluxd_sw(ig,1,1)
          fluxsurf_sw(ig0+ig,2) = zfluxd_sw(ig,1,2)
          fluxtop_sw(ig0+ig,1) = zfluxu_sw(ig,nlaylte+1,1)
          fluxtop_sw(ig0+ig,2) = zfluxu_sw(ig,nlaylte+1,2)
        enddo

      ENDDO         !   (boucle jd=1, ndomain)

c     Zero tendencies for any remaining layers between nlaylte and nlayer
      if (nlayer.gt.nlaylte) then
         do l = nlaylte+1, nlayer
            do ig = 1, ngrid
               dtlw(ig, l) = 0.
               dtsw(ig, l) = 0.
            enddo
         enddo
      endif

c     Output for debugging if lwrite=T
c     --------------------------------
c     Write all nlayer layers, even though only nlaylte layers may have
c     non-zero tendencies.

         IF(lwrite) THEN
            PRINT*,'Diagnotique for the radiation'
            PRINT*,'albedo, emissiv, mu0,fract,fluxsurf_lw,fluxsurf_sw'
            PRINT*,albedo(igout,1),emis(igout),mu0(igout),
     s           fract(igout), fluxsurf_lw(igout),
     $     fluxsurf_sw(igout,1)+fluxsurf_sw(igout,2)
            PRINT*,'Tlay Tlev Play Plev dT/dt SW dT/dt LW (K/s)'
            PRINT*,'daysec',daysec
            DO l=1,nlayer
               PRINT*,pt(igout,l),ptlev(igout,l),
     s         pplay(igout,l),pplev(igout,l),
     s         dtsw(igout,l),dtlw(igout,l)
            ENDDO
         ENDIF


      END SUBROUTINE callradite

      END MODULE callradite_mod 
