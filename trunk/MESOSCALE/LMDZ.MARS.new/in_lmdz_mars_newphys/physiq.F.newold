      SUBROUTINE physiq(ngrid,nlayer,nq,
     $            firstcall,lastcall,
     $            pday,ptime,ptimestep,
     $            pplev,pplay,pphi,
     $            pu,pv,pt,pq,
     $            pw,
     $            pdu,pdv,pdt,pdq,pdpsrf,tracerdyn)


      IMPLICIT NONE
c=======================================================================
c
c   subject:
c   --------
c
c   Organisation of the physical parametrisations of the LMD 
c   martian atmospheric general circulation model.
c
c   The GCM can be run without or with tracer transport
c   depending on the value of Logical "tracer" in file  "callphys.def"
c   Tracers may be water vapor, ice OR chemical species OR dust particles
c
c   SEE comments in initracer.F about numbering of tracer species...
c
c   It includes:
c
c      1. Initialization:
c      1.1 First call initializations
c      1.2 Initialization for every call to physiq
c      1.2.5 Compute mean mass and cp, R and thermal conduction coeff.
c      2. Compute radiative transfer tendencies
c         (longwave and shortwave) for CO2 and dust.
c      3. Gravity wave and subgrid scale topography drag :
c      4. Vertical diffusion (turbulent mixing):
c      5. Convective adjustment
c      6. Condensation and sublimation of carbon dioxide.
c      7.  TRACERS :
c       7a. water and water ice
c       7b. call for photochemistry when tracers are chemical species
c       7c. other scheme for tracer (dust) transport (lifting, sedimentation)
c       7d. updates (CO2 pressure variations, surface budget)
c      8. Contribution to tendencies due to thermosphere
c      9. Surface and sub-surface temperature calculations
c     10. Write outputs :
c           - "startfi", "histfi" (if it's time)
c           - Saving statistics (if "callstats = .true.")
c           - Dumping eof (if "calleofdump = .true.")
c           - Output any needed variables in "diagfi" 
c      10. Diagnostic: mass conservation of tracers
c 
c   author: 
c   ------- 
c           Frederic Hourdin	15/10/93
c           Francois Forget		1994
c           Christophe Hourdin	02/1997 
c           Subroutine completly rewritten by F.Forget (01/2000)
c           Introduction of the photochemical module: S. Lebonnois (11/2002)
c           Introduction of the thermosphere module: M. Angelats i Coll (2002)
c           Water ice clouds: Franck Montmessin (update 06/2003)
c           
c
c
c   arguments:
c   ----------
c
c   input:
c   ------
c    ecri                  period (in dynamical timestep) to write output
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    nq                    Number of advected fields
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c    ptime                 Universal time (0<ptime<1): ptime=0.5 at 12:00 UT
c    ptimestep             timestep (s)
c    pplay(ngrid,nlayer)   Pressure at the middle of the layers (Pa)
c    pplev(ngrid,nlayer+1) intermediate pressure levels (pa)
c    pphi(ngrid,nlayer)    Geopotential at the middle of the layers (m2s-2)
c    pu(ngrid,nlayer)      u component of the wind (ms-1)
c    pv(ngrid,nlayer)      v component of the wind (ms-1)
c    pt(ngrid,nlayer)      Temperature (K)
c    pq(ngrid,nlayer,nq)   Advected fields
c    pudyn(ngrid,nlayer)    \ 
c    pvdyn(ngrid,nlayer)     \ Dynamical temporal derivative for the
c    ptdyn(ngrid,nlayer)     / corresponding variables
c    pqdyn(ngrid,nlayer,nq) /
c    pw(ngrid,?)           vertical velocity
c
c   output:
c   -------
c
c    pdu(ngrid,nlayermx)        \
c    pdv(ngrid,nlayermx)         \  Temporal derivative of the corresponding
c    pdt(ngrid,nlayermx)         /  variables due to physical processes.
c    pdq(ngrid,nlayermx)        /
c    pdpsrf(ngrid)             /
c    tracerdyn                 call tracer in dynamical part of GCM ?

c
c=======================================================================
c
c    0.  Declarations :
c    ------------------

#include "dimensions.h"
#include "dimphys.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "control.h"
#include "dimradmars.h"
#include "comg1d.h"
#include "tracer.h"
#include "nlteparams.h"

#include "chimiedata.h"
#include "watercap.h"
#include "fisice.h"
#include "param.h"
#include "param_v3.h"
#include "conc.h"

#include "netcdf.inc"



c Arguments :
c -----------

c   inputs:
c   -------
      INTEGER ngrid,nlayer,nq
      REAL ptimestep
      REAL pplev(ngridmx,nlayer+1),pplay(ngridmx,nlayer)
      REAL pphi(ngridmx,nlayer)
      REAL pu(ngridmx,nlayer),pv(ngridmx,nlayer)
      REAL pt(ngridmx,nlayer),pq(ngridmx,nlayer,nq)
      REAL pw(ngridmx,nlayer) !Mars pvervel transmit par dyn3d
      REAL zh(ngridmx,nlayermx)      ! potential temperature (K)
      LOGICAL firstcall,lastcall

      REAL pday
      REAL ptime 
      logical tracerdyn

c   outputs:
c   --------
c     physical tendencies
      REAL pdu(ngridmx,nlayer),pdv(ngridmx,nlayer)
      REAL pdt(ngridmx,nlayer),pdq(ngridmx,nlayer,nq)
      REAL pdpsrf(ngridmx) ! surface pressure tendency


c Local saved variables:
c ----------------------
c     aerosol (dust or ice) extinction optical depth  at reference wavelength 
c     "longrefvis" set in dimradmars.h , for one of the "naerkind"  kind of
c      aerosol optical properties  :
      REAL aerosol(ngridmx,nlayermx,naerkind) 

      INTEGER day_ini  ! Initial date of the run (sol since Ls=0) 
      INTEGER icount     ! counter of calls to physiq during the run.
      REAL tsurf(ngridmx)            ! Surface temperature (K)
      REAL tsoil(ngridmx,nsoilmx)    ! sub-surface temperatures (K)
      REAL co2ice(ngridmx)           ! co2 ice surface layer (kg.m-2)  
      REAL albedo(ngridmx,2)         ! Surface albedo in each solar band
      REAL emis(ngridmx)             ! Thermal IR surface emissivity
      REAL dtrad(ngridmx,nlayermx)   ! Net atm. radiative heating rate (K.s-1)
      REAL fluxrad_sky(ngridmx)      ! rad. flux from sky absorbed by surface (W.m-2)
      REAL fluxrad(ngridmx)          ! Net radiative surface flux (W.m-2)
      REAL capcal(ngridmx)           ! surface heat capacity (J m-2 K-1)
      REAL fluxgrd(ngridmx)          ! surface conduction flux (W.m-2)
      REAL qsurf(ngridmx,nqmx)       ! tracer on surface (e.g. kg.m-2)
      REAL q2(ngridmx,nlayermx+1)    ! Turbulent Kinetic Energy 
      INTEGER ig_vl1                 ! Grid Point near VL1   (for diagnostic) 

      SAVE day_ini, icount
      SAVE aerosol, tsurf,tsoil
      SAVE co2ice,albedo,emis, q2
      SAVE capcal,fluxgrd,dtrad,fluxrad,fluxrad_sky,qsurf
      SAVE ig_vl1

      REAL stephan   
      DATA stephan/5.67e-08/  ! Stephan Boltzman constant
      SAVE stephan

c Local variables :
c -----------------


      CHARACTER*80 fichier 
      INTEGER l,ig,ierr,igout,iq, tapphys
      INTEGER iqmin                     ! Used if iceparty engaged

      REAL fluxsurf_lw(ngridmx)      !incident LW (IR) surface flux (W.m-2)
      REAL fluxsurf_sw(ngridmx,2)    !incident SW (solar) surface flux (W.m-2)
      REAL fluxtop_lw(ngridmx)       !Outgoing LW (IR) flux to space (W.m-2)
      REAL fluxtop_sw(ngridmx,2)     !Outgoing SW (solar) flux to space (W.m-2)
c     for clear area (uncovered by clouds) :
      REAL clsurf_lw(ngridmx)      !incident LW (IR) surface flux (W.m-2) 
      REAL clsurf_sw(ngridmx,2)    !incident SW (solar) surface flux (W.m-2)
      REAL cltop_lw(ngridmx)       !Outgoing LW (IR) flux to space (W.m-2)
      REAL cltop_sw(ngridmx,2)     !Outgoing SW (solar) flux to space (W.m-2)
      REAL tauref(ngridmx)           ! Reference column optical depth at 700 Pa
                                     ! (used if active=F) 
      REAL tau(ngridmx,naerkind)     ! Column dust optical depth at each point
      REAL zls                       !  solar longitude (rad)
      REAL zday                      ! date (time since Ls=0, in martian days)
      REAL zzlay(ngridmx,nlayermx)   ! altitude at the middle of the layers
      REAL zzlev(ngridmx,nlayermx+1) ! altitude at layer boundaries
      REAL latvl1,lonvl1             ! Viking Lander 1 point (for diagnostic)


c     Tendancies due to various processes:
      REAL dqsurf(ngridmx,nqmx)
      REAL zdtlw(ngridmx,nlayermx)     ! (K/s)
      REAL zdtsw(ngridmx,nlayermx)     ! (K/s)
      REAL cldtlw(ngridmx,nlayermx)     ! (K/s) LW heating rate for clear area
      REAL cldtsw(ngridmx,nlayermx)     ! (K/s) SW heating rate for clear area
      REAL zdtnirco2(ngridmx,nlayermx) ! (K/s)
      REAL zdtnlte(ngridmx,nlayermx)   ! (K/s)
      REAL zdtsurf(ngridmx)            ! (K/s)
      REAL zdtcloud(ngridmx,nlayermx)
      REAL zdvdif(ngridmx,nlayermx),zdudif(ngridmx,nlayermx)  ! (m.s-2)
      REAL zdhdif(ngridmx,nlayermx), zdtsdif(ngridmx)         ! (K/s)
      REAL zdvadj(ngridmx,nlayermx),zduadj(ngridmx,nlayermx)  ! (m.s-2)
      REAL zdhadj(ngridmx,nlayermx)                           ! (K/s)
      REAL zdtgw(ngridmx,nlayermx)                            ! (K/s)
      REAL zdugw(ngridmx,nlayermx),zdvgw(ngridmx,nlayermx)    ! (m.s-2)
      REAL zdtc(ngridmx,nlayermx),zdtsurfc(ngridmx)
      REAL zdvc(ngridmx,nlayermx),zduc(ngridmx,nlayermx)

      REAL zdqdif(ngridmx,nlayermx,nqmx), zdqsdif(ngridmx,nqmx)
      REAL zdqsed(ngridmx,nlayermx,nqmx), zdqssed(ngridmx,nqmx)
      REAL zdqdev(ngridmx,nlayermx,nqmx), zdqsdev(ngridmx,nqmx)
      REAL zdqadj(ngridmx,nlayermx,nqmx)
      REAL zdqc(ngridmx,nlayermx,nqmx)
      REAL zdqcloud(ngridmx,nlayermx,nqmx)
      REAL zdqscloud(ngridmx,nqmx)
      REAL zdqchim(ngridmx,nlayermx,nqmx)
      REAL zdqschim(ngridmx,nqmx)

      REAL zdteuv(ngridmx,nlayermx)    ! (K/s)
      REAL zdtconduc(ngridmx,nlayermx) ! (K/s)
      REAL zdumolvis(ngridmx,nlayermx)
      REAL zdvmolvis(ngridmx,nlayermx)
      real zdqmoldiff(ngridmx,nlayermx,nqmx)

c     Local variable for local intermediate calcul:
      REAL zflubid(ngridmx)
      REAL zplanck(ngridmx),zpopsk(ngridmx,nlayermx)
      REAL zdum1(ngridmx,nlayermx)
      REAL zdum2(ngridmx,nlayermx)
      REAL ztim1,ztim2,ztim3, z1,z2
      REAL ztime_fin
      REAL zdh(ngridmx,nlayermx)
      REAL pclc_min ! minimum of the cloud fraction over the whole domain
      INTEGER length
      PARAMETER (length=100)

c local variables only used for diagnostic (output in file "diagfi" or "stats")
c -----------------------------------------------------------------------------
      REAL ps(ngridmx), zt(ngridmx,nlayermx)
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx)
      REAL zq(ngridmx,nlayermx,nqmx)
      REAL fluxtop_sw_tot(ngridmx), fluxsurf_sw_tot(ngridmx)
      character*2 str2
      character*5 str5
      real zdtdif(ngridmx,nlayermx), zdtadj(ngridmx,nlayermx)
      real reff(ngridmx,nlayermx) ! effective dust radius (used if doubleq=T)
      real qtot1,qtot2 ! total aerosol mass
      integer igmin, lmin
      logical tdiag


      REAL zplev(ngrid,nlayermx+1),zplay(ngrid,nlayermx)
      REAL zstress(ngrid), cd
      real hco2(nqmx),tmean, zlocal(nlayermx)
      real rho(ngridmx,nlayermx) ! density
      real vmr(ngridmx,nlayermx) ! volume mixing ratio

      REAL time_phys
     
c=======================================================================

c 1. Initialisation:
c -----------------

c  1.1   Initialisation only at first call
c  ---------------------------------------
      IF (firstcall) THEN

c        variables set to 0
c        ~~~~~~~~~~~~~~~~~~
         call zerophys(ngrid*nlayer*naerkind,aerosol)
         call zerophys(ngrid*nlayer,dtrad)
         call zerophys(ngrid,fluxrad)

c        read startfi 
c        ~~~~~~~~~~~~

! Read netcdf initial physical parameters.
         CALL phyetat0 ("startfi.nc",0,0,
     &         nsoilmx,nq,
     &         day_ini,time_phys,
     &         tsurf,tsoil,emis,q2,qsurf,co2ice)

         if (pday.ne.day_ini) then
           write(*,*) "PHYSIQ: ERROR: bad synchronization between ",
     &                "physics and dynamics"
           write(*,*) "dynamics day: ",pday
           write(*,*) "physics day:  ",day_ini
           stop
         endif

         write (*,*) 'In physiq day_ini =', day_ini

c        Initialize albedo and orbital calculation
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL surfini(ngrid,co2ice,qsurf,albedo)
         CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)

c        initialize soil 
c        ~~~~~~~~~~~~~~~
         IF (callsoil) THEN
            CALL soil(ngrid,nsoilmx,firstcall,inertiedat,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
         ELSE
            PRINT*,'WARNING! Thermal conduction in the soil turned off'
            DO ig=1,ngrid
               capcal(ig)=1.e5
               fluxgrd(ig)=0.
            ENDDO
         ENDIF
         icount=1


c        initialize tracers
c        ~~~~~~~~~~~~~~~~~~
         tracerdyn=tracer
         IF (tracer) THEN
            CALL initracer(qsurf,co2ice)
         ENDIF  ! end tracer

c        Determining gridpoint near Viking Lander 1 (used for diagnostic only)
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if(ngrid.ne.1) then
           latvl1= 22.27 
           lonvl1= -47.94 
           ig_vl1= 1+ int( (1.5-(latvl1-90.)*jjm/180.)  -2 )*iim +
     &              int(1.5+(lonvl1+180)*iim/360.)
           write(*,*) 'Viking Lander 1 GCM point: lat,lon',
     &              lati(ig_vl1)*180/pi, long(ig_vl1)*180/pi
         end if 

c        Initialize thermospheric parameters
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (callthermos) call param_read

c        Initialize R and Cp as constant

         if (.not.callthermos .and. .not.photochem) then
                 do l=1,nlayermx
                  do ig=1,ngridmx
                   rnew(ig,l)=r
                   cpnew(ig,l)=cpp
                   mmean(ig,l)=mugaz
                   enddo
                  enddo  
         endif         
                   
      ENDIF        !  (end of "if firstcall")

c ---------------------------------------------------
c 1.2   Initializations done at every physical timestep:
c ---------------------------------------------------
c
      IF (ngrid.NE.ngridmx) THEN
         PRINT*,'STOP in PHYSIQ'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'ngridmx   = ',ngridmx
         STOP
      ENDIF

c     Initialize various variables
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call zerophys(ngrid*nlayer, pdv)
      call zerophys(ngrid*nlayer, pdu)
      call zerophys(ngrid*nlayer, pdt)
      call zerophys(ngrid*nlayer*nq, pdq)
      call zerophys(ngrid, pdpsrf)
      call zerophys(ngrid, zflubid)
      call zerophys(ngrid, zdtsurf)
      call zerophys(ngrid*nq, dqsurf)
      igout=ngrid/2+1 


      zday=pday+ptime ! compute time, in sols (and fraction thereof)

c     Compute Solar Longitude (Ls) :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (season) then
         call solarlong(zday,zls)
      else
         call solarlong(float(day_ini),zls)
      end if

c     Compute geopotential at interlayers
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ponderation des altitudes au niveau des couches en dp/p

      DO l=1,nlayer
         DO ig=1,ngrid
            zzlay(ig,l)=pphi(ig,l)/g
         ENDDO
      ENDDO
      DO ig=1,ngrid
         zzlev(ig,1)=0.
         zzlev(ig,nlayer+1)=1.e7    ! dummy top of last layer above 10000 km...
      ENDDO
      DO l=2,nlayer
         DO ig=1,ngrid
            z1=(pplay(ig,l-1)+pplev(ig,l))/(pplay(ig,l-1)-pplev(ig,l))
            z2=(pplev(ig,l)+pplay(ig,l))/(pplev(ig,l)-pplay(ig,l))
            zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
         ENDDO
      ENDDO


!     Potential temperature calculation not the same in physiq and dynamic

c     Compute potential temperature
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO l=1,nlayer
         DO ig=1,ngrid 
            zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c    1.2.5 Compute mean mass, cp, and R
c    --------------------------------

      if(photochem.or.callthermos) then
         call concentrations(pplay,pt,pdt,pq,pdq,ptimestep)
      endif
c-----------------------------------------------------------------------
c    2. Compute radiative tendencies :
c------------------------------------


      IF (callrad) THEN
         IF( MOD(icount-1,iradia).EQ.0) THEN

c          Local Solar zenith angle
c          ~~~~~~~~~~~~~~~~~~~~~~~~
           CALL orbite(zls,dist_sol,declin)

           IF(diurnal) THEN
               ztim1=SIN(declin)
               ztim2=COS(declin)*COS(2.*pi*(zday-.5))
               ztim3=-COS(declin)*SIN(2.*pi*(zday-.5))

               CALL solang(ngrid,sinlon,coslon,sinlat,coslat,
     s         ztim1,ztim2,ztim3, mu0,fract)

           ELSE
               CALL mucorr(ngrid,declin, lati, mu0, fract,10000.,rad)
           ENDIF

c          NLTE cooling from CO2 emission
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           IF(callnlte) CALL nltecool(ngrid,nlayer,pplay,pt,zdtnlte)

c          Find number of layers for LTE radiation calculations
           IF(MOD(iphysiq*(icount-1),day_step).EQ.0)
     &          CALL nlthermeq(ngrid,nlayer,pplev,pplay)


c          Atmospheric dust opacity and aerosol distribution:
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           CALL dustopacity(ngrid,nlayer,nq,zday,pplay,pplev,zls,pq,
     $      tauref,tau,aerosol)
         
c          Call main radiative transfer scheme
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          Transfer through dust and CO2, except NIR CO2 absorption

           CALL callradite(icount,ngrid,nlayer,aerosol,albedo,
     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
     $     zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,fluxtop_sw)

c          CO2 near infrared absorption
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           call zerophys(ngrid*nlayer,zdtnirco2)
           if (callnirco2) then
              call nirco2abs (ngrid,nlayer,pplay,dist_sol,
     .                       mu0,fract,declin, zdtnirco2)
           endif

c          Radiative flux from the sky absorbed by the surface (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,ngrid
               fluxrad_sky(ig)=emis(ig)*fluxsurf_lw(ig)
     $         +fluxsurf_sw(ig,1)*(1.-albedo(ig,1))
     $         +fluxsurf_sw(ig,2)*(1.-albedo(ig,2))
           ENDDO


c          Net atmospheric radiative heating rate (K.s-1)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF(callnlte) THEN
              CALL blendrad(ngrid, nlayer, pplay,
     &             zdtsw, zdtlw, zdtnirco2, zdtnlte, dtrad)
           ELSE
              DO l=1,nlayer
                 DO ig=1,ngrid
                    dtrad(ig,l)=zdtsw(ig,l)+zdtlw(ig,l)
     &                          +zdtnirco2(ig,l)
                  ENDDO
              ENDDO
           ENDIF



        ENDIF ! of if(mod(icount-1,iradia).eq.0)



c    Transformation of the radiative tendencies:
c    -------------------------------------------

c          Net radiative surface flux (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~
c
           DO ig=1,ngrid
               zplanck(ig)=tsurf(ig)*tsurf(ig)
               zplanck(ig)=emis(ig)*
     $         stephan*zplanck(ig)*zplanck(ig)
               fluxrad(ig)=fluxrad_sky(ig)-zplanck(ig)
           ENDDO


         DO l=1,nlayer
            DO ig=1,ngrid
               pdt(ig,l)=pdt(ig,l)+dtrad(ig,l)
            ENDDO
         ENDDO

      ENDIF ! of IF (callrad)

c-----------------------------------------------------------------------
c    3. Gravity wave and subgrid scale topography drag :
c    -------------------------------------------------


      IF(calllott)THEN

        CALL calldrag_noro(ngrid,nlayer,ptimestep,
     &                 pplay,pplev,pt,pu,pv,zdtgw,zdugw,zdvgw)
 
        DO l=1,nlayer
          DO ig=1,ngrid
            pdv(ig,l)=pdv(ig,l)+zdvgw(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdugw(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtgw(ig,l)
          ENDDO
        ENDDO
      ENDIF

c-----------------------------------------------------------------------
c    4. Vertical diffusion (turbulent mixing):
c    -----------------------------------------
c
      IF (calldifv) THEN


         DO ig=1,ngrid
            zflubid(ig)=fluxrad(ig)+fluxgrd(ig)
         ENDDO

         CALL zerophys(ngrid*nlayer,zdum1)
         CALL zerophys(ngrid*nlayer,zdum2)
         do l=1,nlayer
            do ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            enddo
         enddo

c        Calling vdif (Martian version WITH CO2 condensation)
         CALL vdifc(ngrid,nlayer,nq,co2ice,zpopsk,
     $        ptimestep,capcal,lwrite,
     $        pplay,pplev,zzlay,zzlev,z0,
     $        pu,pv,zh,pq,tsurf,emis,qsurf,
     $        zdum1,zdum2,zdh,pdq,zflubid,
     $        zdudif,zdvdif,zdhdif,zdtsdif,q2,
     &        zdqdif,zdqsdif)

         DO l=1,nlayer
            DO ig=1,ngrid
               pdv(ig,l)=pdv(ig,l)+zdvdif(ig,l)
               pdu(ig,l)=pdu(ig,l)+zdudif(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhdif(ig,l)*zpopsk(ig,l)

               zdtdif(ig,l)=zdhdif(ig,l)*zpopsk(ig,l) ! for diagnostic only

            ENDDO
         ENDDO

         DO ig=1,ngrid
            zdtsurf(ig)=zdtsurf(ig)+zdtsdif(ig)
         ENDDO

         if (tracer) then 
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdif(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
           DO iq=1, nq
              DO ig=1,ngrid
                 dqsurf(ig,iq)=dqsurf(ig,iq) + zdqsdif(ig,iq)
              ENDDO
           ENDDO
         end if ! of if (tracer)

      ELSE    
         DO ig=1,ngrid
            zdtsurf(ig)=zdtsurf(ig)+
     s      (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
         ENDDO
      ENDIF ! of IF (calldifv)


c-----------------------------------------------------------------------
c   5. Dry convective adjustment:
c   -----------------------------

      IF(calladj) THEN

         DO l=1,nlayer
            DO ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            ENDDO
         ENDDO
         CALL zerophys(ngrid*nlayer,zduadj)
         CALL zerophys(ngrid*nlayer,zdvadj)
         CALL zerophys(ngrid*nlayer,zdhadj)

         CALL convadj(ngrid,nlayer,nq,ptimestep,
     $                pplay,pplev,zpopsk,
     $                pu,pv,zh,pq,
     $                pdu,pdv,zdh,pdq,
     $                zduadj,zdvadj,zdhadj,
     $                zdqadj)

         DO l=1,nlayer
            DO ig=1,ngrid
               pdu(ig,l)=pdu(ig,l)+zduadj(ig,l)
               pdv(ig,l)=pdv(ig,l)+zdvadj(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhadj(ig,l)*zpopsk(ig,l)

               zdtadj(ig,l)=zdhadj(ig,l)*zpopsk(ig,l) ! for diagnostic only
            ENDDO
         ENDDO

         if(tracer) then 
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqadj(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         end if
      ENDIF ! of IF(calladj)

c-----------------------------------------------------------------------
c   6. Carbon dioxide condensation-sublimation:
c   -------------------------------------------

      IF (callcond) THEN
         CALL newcondens(ngrid,nlayer,nq,ptimestep,
     $              capcal,pplay,pplev,tsurf,pt,
     $              pphi,pdt,pdu,pdv,zdtsurf,pu,pv,pq,pdq,
     $              co2ice,albedo,emis,
     $              zdtc,zdtsurfc,pdpsrf,zduc,zdvc,zdqc,
     $	            fluxsurf_sw) 

         DO l=1,nlayer
           DO ig=1,ngrid
             pdt(ig,l)=pdt(ig,l)+zdtc(ig,l)
             pdv(ig,l)=pdv(ig,l)+zdvc(ig,l)
             pdu(ig,l)=pdu(ig,l)+zduc(ig,l)
           ENDDO
         ENDDO
         DO ig=1,ngrid
            zdtsurf(ig) = zdtsurf(ig) + zdtsurfc(ig)
         ENDDO

         IF (tracer) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqc(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         ENDIF ! of IF (tracer)

      ENDIF  ! of IF (callcond)

c-----------------------------------------------------------------------
c   7. Specific parameterizations for tracers 
c:   -----------------------------------------

      if (tracer) then 

c   7a. Water and ice
c     ---------------

c        ---------------------------------------
c        Water ice condensation in the atmosphere
c        ----------------------------------------
         IF (water) THEN
           call watercloud(ngrid,nlayer, ptimestep,
     &                pplev,pplay,pdpsrf,zzlev,zzlay, pt,pdt,
     &                pq,pdq,zdqcloud,qsurf,zdqscloud,zdtcloud,
     &                nq,naerkind,tau,icount,zls)

           if (activice) then
c Temperature variation due to latent heat release
           DO l=1,nlayer
             DO ig=1,ngrid
               pdt(ig,l)=pdt(ig,l)+zdtcloud(ig,l)
             ENDDO
           ENDDO
           endif

           IF (iceparty) then
             iqmin=nq-1
           ELSE
             iqmin=nq
           ENDIF

           DO iq=iqmin,nq
             DO l=1,nlayer
                DO ig=1,ngrid
                   pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqcloud(ig,l,iq)
                ENDDO
             ENDDO
             DO ig=1,ngrid
                dqsurf(ig,iq)= dqsurf(ig,iq) + zdqscloud(ig,iq)
             ENDDO
           ENDDO

         END IF  ! of IF (water)

c   7b. Chemical species
c     ------------------

c        --------------
c        photochemistry :
c        --------------
         IF (photochem .or. thermochem) then
          call calchim(ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,
     .      zzlay,zday,pq,pdq,zdqchim,zdqschim,zdqcloud,zdqscloud)
           
c Photochemistry includes condensation of H2O2

           do iq=nqchem_min,nq
            if (noms(iq).eq."h2o2") then
             DO l=1,nlayer
               DO ig=1,ngrid
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqchim(ig,l,iq)
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqcloud(ig,l,iq)
               ENDDO
             ENDDO
            else
             DO l=1,nlayer
               DO ig=1,ngrid
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqchim(ig,l,iq)
               ENDDO
             ENDDO
            endif
           ENDDO
           do iq=nqchem_min,nq
            if (noms(iq).eq."h2o2") then
               DO ig=1,ngrid
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqschim(ig,iq)
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqscloud(ig,iq)
               ENDDO
            else
               DO ig=1,ngrid
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqschim(ig,iq)
               ENDDO
            endif
           ENDDO

         END IF  ! of IF (photochem.or.thermochem)

c   7c. Aerosol particles
c     -------------------

c        ----------
c        Dust devil :
c        ----------
         IF(callddevil) then 
           call dustdevil(ngrid,nlayer,nq, pplev,pu,pv,pt, tsurf,q2,
     &                zdqdev,zdqsdev)
 
           if (dustbin.ge.1) then
              do iq=1,nq
                 DO l=1,nlayer
                    DO ig=1,ngrid
                       pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdev(ig,l,iq)
                    ENDDO
                 ENDDO
              enddo
              do iq=1,nq
                 DO ig=1,ngrid
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqsdev(ig,iq)
                 ENDDO
              enddo
           endif  ! of if (dustbin.ge.1)

         END IF ! of IF (callddevil)

c        ------------- 
c        Sedimentation :   acts also on water ice
c        ------------- 
         IF (sedimentation) THEN 
           call zerophys(ngrid*nlayer*nq, zdqsed)
           call zerophys(ngrid*nq, zdqssed)

           if(doubleq) then
            call callsedim2q(ngrid,nlayer, ptimestep,
     &                pplev,zzlev, pt,
     &                pq, pdq, zdqsed, zdqssed,nq)
           else
            call callsedim(ngrid,nlayer, ptimestep,
     &                pplev,zzlev, pt,
     &                pq, pdq, zdqsed, zdqssed,nq)
           end if


           DO iq=1, nq
             DO l=1,nlayer
               DO ig=1,ngrid
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqsed(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO
           DO iq=1, nq
             DO ig=1,ngrid
                dqsurf(ig,iq)= dqsurf(ig,iq) + zdqssed(ig,iq)
             ENDDO
           ENDDO
         END IF   ! of IF (sedimentation)

c   7d. Updates
c     ---------

        DO iq=1, nq
          DO ig=1,ngrid

c       ---------------------------------
c       Updating tracer budget on surface
c       ---------------------------------
            qsurf(ig,iq)=qsurf(ig,iq)+ptimestep*dqsurf(ig,iq)

          ENDDO  ! (ig)
        ENDDO    ! (iq)

      endif !  of if (tracer) 


c-----------------------------------------------------------------------
c   8. THERMOSPHERE CALCULATION
c-----------------------------------------------------------------------

      if (callthermos) then
        call thermosphere(pplev,pplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)
c           do iq=nqchem_min,nq
c           write(*,*) 'thermo iq,pq',iq,pq(690,1,iq)
c           enddo

        DO l=1,nlayer
          DO ig=1,ngrid
            dtrad(ig,l)=dtrad(ig,l)+zdteuv(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtconduc(ig,l)
     &                         +zdteuv(ig,l)
            pdv(ig,l)=pdv(ig,l)+zdvmolvis(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdumolvis(ig,l)
            DO iq=1, nq
              pdq(ig,l,iq)=pdq(ig,l,iq)+zdqmoldiff(ig,l,iq)
            ENDDO
          ENDDO
        ENDDO


      endif
c-----------------------------------------------------------------------
c   9. Surface  and sub-surface soil temperature
c-----------------------------------------------------------------------
c
c
c   9.1 Increment Surface temperature:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO ig=1,ngrid
         tsurf(ig)=tsurf(ig)+ptimestep*zdtsurf(ig) 
      ENDDO

c  Prescribe a cold trap at south pole (except at high obliquity !!)
c  Temperature at the surface is set there to be the temperature
c  corresponding to equilibrium temperature between phases of CO2

      IF (tracer.AND.water.AND.ngridmx.NE.1) THEN
         if (caps.and.(obliquit.lt.27.)) then
           ! NB: Updated surface pressure, at grid point 'ngrid', is
           !     ps(ngrid)=pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep
           tsurf(ngrid)=1./(1./136.27-r/5.9e+5*alog(0.0095*
     &                     (pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep)))
         endif
c       -------------------------------------------------------------
c       Change of surface albedo (set to 0.4) in case of ground frost
c       everywhere except on the north permanent cap and in regions
c       covered by dry ice. 
c              ALWAYS PLACE these lines after newcondens !!!
c       -------------------------------------------------------------
         do ig=1,ngrid

c       -------------- Version temporaire fit TES 2008 ------------
         if (co2ice(ig).eq.0.and.qsurf(ig,nqmx).gt.0.005) then
              albedo(ig,1)=0.45                                     
              albedo(ig,2)=0.45                                     
          endif

c         if (co2ice(ig).eq.0.and.qsurf(ig,nqmx).gt.0.005) then
c           if (.not.watercaptag(ig)) then
c             albedo(ig,1)=0.4
c             albedo(ig,2)=0.4
c           endif
c         endif
c       -------------- version Francois ---------------
c          if (co2ice(ig).eq.0.and.
c    &        ((qsurf(ig,nqmx).gt.0.005).or.(watercaptag(ig)))) then
c              albedo(ig,1)=max(0.4,albedodat(ig))
c              albedo(ig,2)=albedo(ig,1)
c          endif
c       ---------------------------------------------
         enddo  ! of do ig=1,ngrid
      ENDIF  ! of IF (tracer.AND.water.AND.ngridmx.NE.1)

c
c   9.2 Compute soil temperatures and subsurface heat flux:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (callsoil) THEN
         CALL soil(ngrid,nsoilmx,.false.,inertiedat,
     &          ptimestep,tsurf,tsoil,capcal,fluxgrd)
      ENDIF

c-----------------------------------------------------------------------
c  10. Write output files
c  ----------------------

c    -------------------------------
c    Dynamical fields incrementation
c    -------------------------------
c (FOR OUTPUT ONLY : the actual model integration is performed in the dynamics)
      ! temperature, zonal and meridional wind
      DO l=1,nlayer
        DO ig=1,ngrid
          zt(ig,l)=pt(ig,l)  + pdt(ig,l)*ptimestep
          zu(ig,l)=pu(ig,l)  + pdu(ig,l)*ptimestep
          zv(ig,l)=pv(ig,l)  + pdv(ig,l)*ptimestep
        ENDDO
      ENDDO

      ! tracers
      DO iq=1, nq
        DO l=1,nlayer
          DO ig=1,ngrid
            zq(ig,l,iq)=pq(ig,l,iq) +pdq(ig,l,iq)*ptimestep
          ENDDO
        ENDDO
      ENDDO

      ! surface pressure
      DO ig=1,ngrid
          ps(ig)=pplev(ig,1) + pdpsrf(ig)*ptimestep
      ENDDO

      ! pressure
      DO l=1,nlayer
        DO ig=1,ngrid
             zplev(ig,l)=pplev(ig,l)/pplev(ig,1)*ps(ig)
             zplay(ig,l)=pplay(ig,l)/pplev(ig,1)*ps(ig)
        ENDDO
      ENDDO

      ! Density 
      DO l=1,nlayer
         DO ig=1,ngrid
            rho(ig,l) = zplay(ig,l)/(rnew(ig,l)*zt(ig,l))
         ENDDO
      ENDDO

c    Compute surface stress : (NB: z0 is a common in planete.h)
c     DO ig=1,ngrid
c        cd = (0.4/log(zzlay(ig,1)/z0))**2
c        zstress(ig) = rho(ig,1)*cd*(zu(ig,1)**2 + zv(ig,1)**2)
c     ENDDO

c     Sum of fluxes in solar spectral bands (for output only)
      DO ig=1,ngrid
	     fluxtop_sw_tot(ig)=fluxtop_sw(ig,1) + fluxtop_sw(ig,2)
	     fluxsurf_sw_tot(ig)=fluxsurf_sw(ig,1) + fluxsurf_sw(ig,2)
      ENDDO
c ******* TEST ******************************************************
      ztim1 = 999
      DO l=1,nlayer
        DO ig=1,ngrid
           if (pt(ig,l).lt.ztim1) then
               ztim1 = pt(ig,l)
               igmin = ig
               lmin = l 
           end if
        ENDDO
      ENDDO
      if(min(pt(igmin,lmin),zt(igmin,lmin)).lt.70.) then	
        write(*,*) 'PHYSIQ: stability WARNING :'
        write(*,*) 'pt, zt Tmin = ', pt(igmin,lmin), zt(igmin,lmin),
     &              'ig l =', igmin, lmin
      end if
c *******************************************************************

c     ---------------------
c     Outputs to the screen 
c     ---------------------

      IF (lwrite) THEN
         PRINT*,'Global diagnostics for the physics'
         PRINT*,'Variables and their increments x and dx/dt * dt'
         WRITE(*,'(a6,a10,2a15)') 'Ts','dTs','ps','dps'
         WRITE(*,'(2f10.5,2f15.5)')
     s   tsurf(igout),zdtsurf(igout)*ptimestep,
     s   pplev(igout,1),pdpsrf(igout)*ptimestep
         WRITE(*,'(a4,a6,5a10)') 'l','u','du','v','dv','T','dT'
         WRITE(*,'(i4,6f10.5)') (l,
     s   pu(igout,l),pdu(igout,l)*ptimestep,
     s   pv(igout,l),pdv(igout,l)*ptimestep,
     s   pt(igout,l),pdt(igout,l)*ptimestep,
     s   l=1,nlayer)
      ENDIF ! of IF (lwrite)


      IF (ngrid.NE.1) THEN
         print*,'Ls =',zls*180./pi,
     &   ' tauref(700 Pa,lat=0) =',tauref(ngrid/2),
     &   ' tau(Viking1) =',tau(ig_vl1,1)

c        -------------------------------------------------------------------
c        Writing NetCDF file  "RESTARTFI" at the end of the run
c        -------------------------------------------------------------------
c        Note: 'restartfi' is stored just before dynamics are stored
c              in 'restart'. Between now and the writting of 'restart',
c              there will have been the itau=itau+1 instruction and
c              a reset of 'time' (lastacll = .true. when itau+1= itaufin)
c              thus we store for time=time+dtvr

         IF(lastcall) THEN
            ztime_fin = ptime + ptimestep/(float(iphysiq)*daysec) 
            write(*,*)'PHYSIQ: for physdem ztime_fin =',ztime_fin
            call physdem1("restartfi.nc",long,lati,nsoilmx,nq,
     .              ptimestep,pday,
     .              ztime_fin,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .              area,albedodat,inertiedat,zmea,zstd,zsig,
     .              zgam,zthe)
         ENDIF

c        -----------------------------------------------------------------
c        Saving statistics :
c        -----------------------------------------------------------------
c        ("stats" stores and accumulates 8 key variables in file "stats.nc"
c        which can later be used to make the statistic files of the run:
c        "stats")          only possible in 3D runs !

         
         IF (callstats) THEN


            call wstats(ngrid,"ps","Surface pressure","K",2,ps)
            call wstats(ngrid,"tsurf","Surface temperature","K",2,tsurf)
            call wstats(ngrid,"co2ice","CO2 ice cover",
     .                  "kg.m-2",2,co2ice)
c            call wstats(ngrid,"emis","Surface emissivity","w.m-1",2,
c    .                  emis)
             call wstats(ngrid,"fluxsurf_lw",
     .                 "Thermal IR radiative flux to surface","W.m-2",2,
     .                 fluxsurf_lw)
             call wstats(ngrid,"fluxsurf_sw",
     .                  "Solar radiative flux to surface","W.m-2",2,
     .                   fluxsurf_sw_tot)
             call wstats(ngrid,"fluxtop_lw",
     .                  "Thermal IR radiative flux to space","W.m-2",2,
     .                  fluxtop_lw)
             call wstats(ngrid,"fluxtop_sw",
     .                  "Solar radiative flux to space","W.m-2",2,
     .                  fluxtop_sw_tot)
            call wstats(ngrid,"dod","Dust optical depth"," ",2,tau)

            call wstats(ngrid,"temp","Atmospheric temperature","K",3,zt)
            call wstats(ngrid,"u","Zonal (East-West) wind","m.s-1",3,zu)
            call wstats(ngrid,"v","Meridional (North-South) wind",
     .                  "m.s-1",3,zv)
            call wstats(ngrid,"w","Vertical (down-up) wind",
     .                  "m.s-1",3,pw)
            call wstats(ngrid,"rho","Atmospheric density","none",3,rho)
             call wstats(ngrid,"q2",
     .               "Boundary layer eddy kinetic energy","m2.s-2",3,q2)

           if (tracer) then
            if (water) then
               vmr=zq(1:ngridmx,1:nlayermx,nqmx)*mugaz/mmol(nqmx)
               call wstats(ngrid,"vmr_h2ovapor",
     .                    "H2O vapor volume mixing ratio","mol/mol",
     .                    3,vmr)
               if (iceparty) then
                  vmr=zq(1:ngridmx,1:nlayermx,nqmx-1)*mugaz/mmol(nqmx-1)
                  call wstats(ngrid,"vmr_h2oice",
     .                       "H2O ice volume mixing ratio","mol/mol",
     .                       3,vmr)
               endif
            endif

             if (thermochem.or.photochem) then
                do iq=1,nq
                   if ((noms(iq).eq."o").or.(noms(iq).eq."co2").or.
     .               (noms(iq).eq."co").or.(noms(iq).eq."n2").or.
     .               (noms(iq).eq."h2").or.
     .               (noms(iq).eq."o3")) then
                        do l=1,nlayer
                          do ig=1,ngrid
                            vmr(ig,l)=zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
                          end do
                        end do
                        call wstats(ngrid,"vmr_"//trim(noms(iq)),
     .                     "Volume mixing ratio","mol/mol",3,vmr)
                   endif
                enddo
             endif
           endif !tracer 

            IF(lastcall) THEN
              write (*,*) "Writing stats..."
              call mkstats(ierr)
            ENDIF
          ENDIF !if callstats

c        (Store EOF for Mars Climate database software)
         IF (calleofdump) THEN
            CALL eofdump(ngrid, nlayer, zu, zv, zt, rho, ps)
         ENDIF


c        ----------------------------------------------------------------------
c        output in netcdf file "DIAGFI", containing any variable for diagnostic
c        (output with  period "ecritphy", set in "run.def")
c        ----------------------------------------------------------------------
c        WRITEDIAGFI can ALSO be called from any other subroutines
c        for any variables !!
       call WRITEDIAGFI(ngrid,"emis","Surface emissivity","w.m-1",2,
     .                  emis)
        call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",2,
     .                  tsurf)
        call WRITEDIAGFI(ngrid,"ps","surface pressure","K",2,ps)
        call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness","kg.m-2",2,
     .                  co2ice)
c       call WRITEDIAGFI(ngrid,"fluxsurf_lw","fluxsurf_lw","W.m-2",2,
c    .       fluxsurf_lw)
c       call WRITEDIAGFI(ngrid,"fluxsurf_sw","fluxsurf_sw","W.m-2",2,
c    .       fluxsurf_sw_tot)
c       call WRITEDIAGFI(ngrid,"fluxtop_lw","fluxtop_lw","W.m-2",2,
c    .       fluxtop_lw)
c       call WRITEDIAGFI(ngrid,"fluxtop_sw","fluxtop_sw","W.m-2",2,
c    .       fluxtop_sw_tot)
        call WRITEDIAGFI(ngrid,"temp","temperature","K",3,zt)
c        call WRITEDIAGFI(ngrid,"tau","tau"," ",2,tau)
        call WRITEDIAGFI(ngrid,"u","Zonal wind","m.s-1",3,zu)
        call WRITEDIAGFI(ngrid,"v","Meridional wind","m.s-1",3,zv)
        call WRITEDIAGFI(ngrid,"w","Vertical wind","m.s-1",3,pw)
c       call WRITEDIAGFI(ngrid,"rho","density","none",3,rho)
c        call WRITEDIAGFI(ngrid,"q2","q2","kg.m-3",3,q2)
c       call WRITEDIAGFI(ngridm,'Teta','T potentielle','K',3,zh)
c       call WRITEDIAGFI(ngridm,'Pression','Pression','Pa',3,pplay)
        call WRITEDIAGFI(ngrid,"tsoil","soil temperature","K",3,tsoil)


c       OUTPUT of tracer mass mixing ratio and surface value :    
       if (tracer) then
c         (for photochemistry, outputs are done in calchim)
          do iq=1,nqmx
            write(str2(1:2),'(i2.2)') iq
            call WRITEDIAGFI(ngridmx,'q'//str2,noms(iq),
     &                         'kg/kg',3,zq(1,1,iq))
            call WRITEDIAGFI(ngridmx,'qsurf'//str2,noms(iq),
     &                     'kg.m-2',2,qsurf(1,iq))
          end do
        end if
c ***************************************************

c  Outputs for H2O
       if (tracer) then
        if (activice) then
c         call WRITEDIAGFI(ngridmx,'tauice','tau','SI',2,tau(1,2))
c         call WRITEDIAGFI(ngridmx,'sw_htrt','sw heat. rate',
c     &                   'w.m-2',3,zdtsw)
c         call WRITEDIAGFI(ngridmx,'lw_htrt','lw heat. rate',
c     &                   'w.m-2',3,zdtlw)
        endif  !(activice)

        if (water.and..not.photochem) then
           iq=nq
c          write(str2(1:2),'(i2.2)') iq
c          call WRITEDIAGFI(ngridmx,'dqs'//str2,'dqscloud',
c     &                    'kg.m-2',2,zdqscloud(1,iq))
c          call WRITEDIAGFI(ngridmx,'dqch'//str2,'var chim',
c     &                    'kg/kg',3,zdqchim(1,1,iq))
c          call WRITEDIAGFI(ngridmx,'dqd'//str2,'var dif',
c     &                    'kg/kg',3,zdqdif(1,1,iq))
c          call WRITEDIAGFI(ngridmx,'dqa'//str2,'var adj',
c     &                    'kg/kg',3,zdqadj(1,1,iq))
c          call WRITEDIAGFI(ngridmx,'dqc'//str2,'var c',
c     &                    'kg/kg',3,zdqc(1,1,iq))
        endif  !(water.and..not.photochem)

c       if (iceparty) then
c          iq=nq-1
c          write(str2(1:2),'(i2.2)') iq
c          call WRITEDIAGFI(ngridmx,'q'//str2,'mix. ratio',
c    &                     'kg/kg',3,zq(1,1,iq))
c       endif  !(iceparty)
       endif

c  Outputs for dust tracers

         if (tracer.and.(dustbin.ne.0)) then
            call WRITEDIAGFI(ngridmx,'tau','taudust','SI',2,tau(1,1))
            if (doubleq) then
               call WRITEDIAGFI(ngridmx,'qsurf','qsurf',
     &                         'kg.m-2',2,qsurf(1,1))
               call WRITEDIAGFI(ngridmx,'Nsurf','N particles',
     &                         'N.m-2',2,qsurf(1,2))
               call WRITEDIAGFI(ngridmx,'dqsdev','ddevil lift',
     &                         'kg.m-2.s-1',2,zdqsdev(1,1))
               call WRITEDIAGFI(ngridmx,'dqssed','sedimentation',
     &                         'kg.m-2.s-1',2,zdqssed(1,1))
               do  l=1,nlayer
                  do ig=1, ngrid
                     reff(ig,l)= ref_r0 *
     &               (r3n_q*pq(ig,l,1)/max(pq(ig,l,2),0.01))**(1./3.)
                     reff(ig,l)=min(max(reff(ig,l),1.e-10),500.e-6)
                  end do
               end do
               call WRITEDIAGFI(ngridmx,'reff','reff','m',3,reff)
            else
               do iq=1,dustbin
                  write(str2(1:2),'(i2.2)') iq
                  call WRITEDIAGFI(ngridmx,'q'//str2,'mix. ratio',
     &                            'kg/kg',3,zq(1,1,iq))
                  call WRITEDIAGFI(ngridmx,'qsurf'//str2,'qsurf',
     &                            'kg.m-2',2,qsurf(1,iq))
               end do
            endif ! (doubleq)
         end if  ! (tracer.and.(dustbin.ne.0))


      ELSE     ! if(ngrid.eq.1)

         print*,'Ls =',zls*180./pi,
     &  '  tauref(700 Pa) =',tauref,' local tau =',tau(1,1)
c      ----------------------------------------------------------------------
c      Output in grads file "g1d" (ONLY when using testphys1d)
c      (output at every X physical timestep)
c      ----------------------------------------------------------------------
c
c        CALL writeg1d(ngrid,1,fluxsurf_lw,'Fs_ir','W.m-2')
c        CALL writeg1d(ngrid,1,tsurf,'tsurf','K')
         CALL writeg1d(ngrid,1,ps,'ps','Pa')
         
         CALL writeg1d(ngrid,nlayer,zt,'T','K')
c        CALL writeg1d(ngrid,nlayer,pu,'u','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pv,'v','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pw,'w','m.s-1')

         if(tracer) then
c           CALL writeg1d(ngrid,1,tau,'tau','SI')
            do iq=1,nq
              CALL writeg1d(ngrid,nlayer,zq(1,1,iq),noms(iq),'kg/kg') 
            end do
         end if

         zlocal(1)=-log(pplay(1,1)/pplev(1,1))* Rnew(1,1)*zt(1,1)/g

         do l=2,nlayer-1
            tmean=zt(1,l)
            if(zt(1,l).ne.zt(1,l-1))
     &        tmean=(zt(1,l)-zt(1,l-1))/log(zt(1,l)/zt(1,l-1))
              zlocal(l)= zlocal(l-1)
     &        -log(pplay(1,l)/pplay(1,l-1))*rnew(1,l)*tmean/g
         enddo
         zlocal(nlayer)= zlocal(nlayer-1)-
     &                   log(pplay(1,nlayer)/pplay(1,nlayer-1))*
     &                   rnew(1,nlayer)*tmean/g

c         if(tracer) then
c         do l=2,nlayer
c            do iq=1,nq
c               hco2(iq)=zq(1,l,iq)/zq(1,l-1,iq)
c               hco2(iq)=-(zlocal(l)-zlocal(l-1))/log(hco2(iq))/1000.
c            enddo
c            write(225,*) l,pt(1,l),(hco2(iq),iq=1,6)
c            write(226,*) l,(hco2(iq),iq=7,13)
c         enddo
c         endif

      END IF       ! if(ngrid.ne.1)


      icount=icount+1
      RETURN
      END
