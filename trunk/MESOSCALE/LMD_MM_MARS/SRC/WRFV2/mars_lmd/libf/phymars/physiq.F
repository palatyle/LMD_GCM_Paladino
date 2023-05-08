      SUBROUTINE physiq(ngrid,nlayer,nq,
     $            firstcall,lastcall,
     $            pday,ptime,ptimestep,
     $            pplev,pplay,pphi,
     $            pu,pv,pt,pq,
     $            pw,
     $            pdu,pdv,pdt,pdq,pdpsrf,tracerdyn,
     $            wtsurf,wtsoil,wemis,wq2,wqsurf,wco2ice,
     $            wday_ini,
     $            output_tab2d, output_tab3d,
     $            flag_LES)


      IMPLICIT NONE
c=======================================================================
c
c	CAREFUL: THIS IS A VERSION TO BE USED WITH WRF !!!
c
c	... CHECK THE ****WRF lines
c
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
c      1. Initialisation:
c      1.5 Calculation of mean mass and cp, R and thermal conduction coeff.
c      2. Calcul of the radiative tendencies : radiative transfer
c         (longwave and shortwave) for CO2 and dust.
c      3. Gravity wave and subgrid scale topography drag :
c      4. Vertical diffusion (turbulent mixing):
c      5. convective adjustment
c      6.  TRACERS :
c       6a. water and water ice
c       6b. call for photochemistry when tracers are chemical species
c       6c. other scheme for tracer (dust) transport (lifting, sedimentation)
c       6d. updates (CO2 pressure variations, surface budget)
c      7.  condensation and sublimation of carbon dioxide.
c      8. Surface and sub-surface temperature calculations
c      9. Writing output files :
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
c	    WRF version: Aymeric Spiga (01-03/2007)
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
c
c    ****WRF 
c	day_ini,tsurf,tsoil,emis,q2,qsurf,co2ice are inputs 
c		and locally saved variables
c    			(no need to call phyetat0)
c
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

#include "slope.h"
#include "wrf_output_2d.h"
#include "wrf_output_3d.h"


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
c ****WRF
      INTEGER wday_ini  
      REAL wtsurf(ngridmx)  ! input only ay firstcall - output            
      REAL wtsoil(ngridmx,nsoilmx)    
      REAL wco2ice(ngridmx)             
      REAL wemis(ngridmx)             
      REAL wqsurf(ngridmx,nqmx)       
      REAL wq2(ngridmx,nlayermx+1) 
      REAL output_tab2d(ngridmx,n2d)
      REAL output_tab3d(ngridmx,nlayer,n3d)
      REAL sl_ls, sl_lct, sl_lat, sl_tau, sl_alb, sl_the, sl_psi
      REAL sl_fl0, sl_flu
      REAL sl_ra, sl_di0
      REAL sky
      REAL sensheat(ngridmx)    !! pour LES avec isfflx!=0
      REAL ustar(ngridmx)    !! pour LES avec isfflx!=0
      LOGICAL flag_LES     !! pour LES avec isfflx!=0
      REAL qsurflast(ngridmx) !! pour diagnostics
c ****WRF
      REAL pday
      REAL ptime 
      logical tracerdyn
   

c   outputs:
c   --------
c     physical tendencies
      REAL pdu(ngridmx,nlayer),pdv(ngridmx,nlayer)
      REAL pdt(ngridmx,nlayer),pdq(ngridmx,nlayer,nq)
      REAL pdpsrf(ngridmx)

!!!!!!!TEST TEST
!      REAL spdu(ngridmx,nlayermx)
!      REAL spdv(ngridmx,nlayermx)
!      REAL spdt(ngridmx,nlayermx)
!      REAL spdq(ngridmx,nlayermx,nqmx)
!      REAL spdpsrf(ngridmx)
!      SAVE spdu
!      SAVE spdv
!      SAVE spdt
!      SAVE spdq
!      SAVE spdpsrf
!      LOGICAL nocalculphy
!!!!!!!TEST TEST      


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
      REAL fluxrad_save(ngridmx)     ! Net radiative surface flux (W.m-2)
      REAL capcal(ngridmx)           ! surface heat capacity (J m-2 K-1)
      REAL fluxgrd(ngridmx)          ! surface conduction flux (W.m-2)
      REAL qsurf(ngridmx,nqmx)       ! tracer on surface (e.g. kg.m-2)
      REAL q2(ngridmx,nlayermx+1)    ! Turbulent Kinetic Energy 
      INTEGER ig_vl1                 ! Grid Point near VL1   (for diagnostic) 

      SAVE day_ini, icount
      SAVE aerosol, tsurf,tsoil
      SAVE co2ice,albedo,emis, q2
      SAVE capcal,fluxgrd,dtrad,fluxrad_save, qsurf
      SAVE ig_vl1

      REAL stephan   
      DATA stephan/5.67e-08/  ! Stephan Boltzman constant
      SAVE stephan

c Local variables :
c -----------------

      REAL fluxrad(ngridmx)     ! Net radiative surface flux (W.m-2)

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
      real hco2(nqmx),tmean, zlocal(nlayermx)
      real rho(ngridmx,nlayermx) ! density
      real vmr(ngridmx,nlayermx) ! volume mixing ratio

      REAL time_phys

!!! WRF for retrocompatibility with newphys
      REAL tauTES(ngridmx)      ! column optical depth at 825 cm-1
      REAL qsurfice(ngridmx)    ! pour diagnostics
!!! WRF for retrocompatibility with newphys
     
c=======================================================================




c 1. Initialisation:
c -----------------
 

c     Initialisation only at first call
c     ---------------------------------
      IF(firstcall) THEN

c        variables set to 0
c        ~~~~~~~~~~~~~~~~~~
         call zerophys(ngrid*nlayer*naerkind,aerosol)
         call zerophys(ngrid*nlayer,dtrad)
         call zerophys(ngrid,fluxrad)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ****WRF
c
c	No need to use startfi.nc
c		> part of the job of phyetat0 is done in inifis
c		> remaining initializations are passed here from the WRF variables
c		> beware, some operations were done by phyetat0 (ex: tracers)
c			> if any problems, look in phyetat0
c
      tsurf(:)=wtsurf(:)
      PRINT*,'check: tsurf ',tsurf(1),tsurf(ngridmx)
      tsoil(:,:)=wtsoil(:,:)
      PRINT*,'check: tsoil ',tsoil(1,1),tsoil(ngridmx,nsoilmx)
      emis(:)=wemis(:)
      PRINT*,'check: emis ',emis(1),emis(ngridmx)
      q2(:,:)=wq2(:,:)
      PRINT*,'check: q2 ',q2(1,1),q2(ngridmx,nlayermx+1)
      qsurf(:,:)=wqsurf(:,:)
      PRINT*,'check: qsurf ',qsurf(1,1),qsurf(ngridmx,nqmx)
      co2ice(:)=wco2ice(:)
      PRINT*,'check: co2 ',co2ice(1),co2ice(ngridmx)
      day_ini=wday_ini

c	artificially filling dyn3d/control.h is also required
c	> iphysiq is put in WRF to be set easily (cf ptimestep)
c	> day_step is simply deduced:
c
      day_step=daysec/ptimestep
      PRINT*,'Call to LMD physics:',day_step,' per Martian day'
c
      iphysiq=ptimestep
c
      ecritphy=8.e18  !! a dummy low frequency
      PRINT*,'Write LMD physics each:',ecritphy,' seconds'
              !!PRINT*,ecri_phys  
              !!PRINT*,float(ecri_phys) ...
              !!renvoient tous deux des nombres absurdes 
              !!pourtant callkeys.h est inclus ...
              !!
              !!donc ecritphys est passe en argument ... 
      PRINT*,'Write LMD physics each:',ecritphy,' seconds'
c
c ****WRF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         if (pday.ne.day_ini) then
           write(*,*) "***ERROR Pb de synchronisation entre phys et dyn"
           write(*,*) "jour dynamique: ",pday
           write(*,*) "jour physique: ",day_ini
           stop
         endif

         write (*,*) 'In physic day_ini =', day_ini



c        Initialize albedo and orbital calculation
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL surfini(ngrid,co2ice,qsurf,albedo)
         CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)

c        initialisation soil 
c        ~~~~~~~~~~~~~~~~~~~
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


c        initialisation pour les traceurs
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tracerdyn=tracer
         IF (tracer) THEN
            CALL initracer(qsurf,co2ice)
         ENDIF  ! end tracer


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ****WRF
c
c nosense in mesoscale modeling
c
cc        Determining gridpoint near Viking Lander 1 (used for diagnostic only)
cc        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         if(ngrid.ne.1) then
c           latvl1= 22.27 
c           lonvl1= -47.94 
c           ig_vl1= 1+ int( (1.5-(latvl1-90.)*jjm/180.)  -2 )*iim +
c     &              int(1.5+(lonvl1+180)*iim/360.)
c           write(*,*) 'Viking Lander 1 GCM point: lat,lon',
c     &              lati(ig_vl1)*180/pi, long(ig_vl1)*180/pi
c         end if 
c
c ****WRF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c        Initializing thermospheric parameters
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (callthermos) call param_read

c        Initializing R and Cp as constant

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

!!!!!!!!!!!!!!!!TEST TEST
!      IF (nocalculphy) THEN
!
!            write(*,*) 'tendencies are not recalculated !'  
!            pdu(:,:)=spdu(:,:) 
!            pdv(:,:)=spdv(:,:)
!            pdt(:,:)=spdt(:,:)
!            pdq(:,:,:)=spdq(:,:,:)
!            pdpsrf(:)=spdpsrf(:)
!            write(*,*) pdu(100,10), pdv(100,10), pdt(100,10)  
!      ELSE      
!!!!!!!!!!!!!!!!TEST TEST


c    ------------------------------------------
c    Initialisations at every physical timestep:
c    ------------------------------------------
c
      IF (ngrid.NE.ngridmx) THEN
         PRINT*,'STOP in PHYSIQ'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'ngridmx   = ',ngridmx
         STOP
      ENDIF

      zday=pday+ptime



c     Computing Solar Longitude (Ls) :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (season) then 
         PRINT*,'day',zday
         CALL solarlong(zday,zls)
      else
         PRINT*,'day_ini',day_ini
         call solarlong(float(day_ini),zls) 
      end if


c     Initializing various variable
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call zerophys(ngrid*nlayer, pdv)
      call zerophys(ngrid*nlayer, pdu)
      call zerophys(ngrid*nlayer, pdt)
      call zerophys(ngrid*nlayer*nq, pdq)
      call zerophys(ngrid, pdpsrf)
      call zerophys(ngrid, zflubid)
      call zerophys(ngrid, zdtsurf)
      call zerophys(ngrid*nq, dqsurf)
      igout=ngrid/2+1 

c     computing geopotentiel at interlayer levels
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

c     Computing potential temperature
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO l=1,nlayer
         DO ig=1,ngrid 
            zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c    1.5 Calculation of mean mass, cp, and R
c    ---------------------------------------

      if(photochem.or.callthermos) then
         call concentrations(pplay,pt,pdt,pq,pdq,ptimestep)
      endif
c-----------------------------------------------------------------------
c    2. Calcul of the radiative tendencies :
c    ---------------------------------------


      IF(callrad) THEN
         IF( MOD(icount-1,iradia).EQ.0) THEN

           write (*,*) 'call radiative transfer'

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
           IF(MOD(iphysiq*(icount-1),day_step).EQ.0) THEN          
                CALL nlthermeq(ngrid,nlayer,pplev,pplay)
           ENDIF


c          Atmospheric dust opacity and aerosol distribution:
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc*****WRF
         CALL meso_dustopacity(ngrid,nlayer,nq,
     $          zday,pplay,pplev,zls,pq,
     $      tauref,tau,aerosol)

cc*****WRF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        
c          Calling main radiative transfer scheme
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          Transfer through dust and CO2, except NIR CO2 absorption

c ----------
c partie rajoutee par Franck, commentee pour l'instant
c ----------
c
c           if (ngridmx.eq.1) pclc(1)=1. !TEST for 1D simulation
c
c           pclc_min=1.
c           if (activice.and.naerkind.gt.1) then
c             do ig=1,ngrid
c               pclc_min=min(pclc_min,pclc(ig))
c             enddo
c           endif
c
c           IF(activice.and.naerkind.gt.1.and.pclc_min.lt.1) THEN
c          Computing radiative tendencies accounting for a cloudy area (0< pclc(ngrid) <1)
c          pclc is set in initracer (prescribed for the moment).
c          two steps : 1/rad. tend. without clouds (aerosol(*,*,2)=0.)
c          ~~~~~~~~~   2/rad. tend. with clouds
c                      3/final tendencies=average of 1/ and 2/ weighted by the
cloud area (pclc)
c
c
c 1/
c           call zerophys(nlayer*ngrid,aerosol(1,1,2))  !remettre
c           CALL callradite(icount,ngrid,nlayer,aerosol,albedo,
c     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
c     $     cldtlw,cldtsw,clsurf_lw,clsurf_sw,cltop_lw,cltop_sw)
c
c 2/
c           CALL h2oiceopacity(ngrid,nlayer,nq,pplay,pplev,pt,pq,
c     $     tau,aerosol,zls,co2ice)
c
c           CALL callradite(icount,ngrid,nlayer,aerosol,albedo,
c     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
c     $     zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,fluxtop_sw)
c
c 3/
c           do l=1,nlayer
c             do ig=1,ngrid
c             zdtlw(ig,l)=(1.-pclc(ig))*cldtlw(ig,l)+pclc(ig)*zdtlw(ig,l)
c             zdtsw(ig,l)=(1.-pclc(ig))*cldtsw(ig,l)+pclc(ig)*zdtsw(ig,l)
c             enddo
c           enddo
c           do ig=1,ngrid
c             fluxsurf_lw(ig)=(1.-pclc(ig))*clsurf_lw(ig)+
c     $        pclc(ig)*fluxsurf_lw(ig)
c             fluxtop_lw(ig)=(1.-pclc(ig))*cltop_lw(ig)+
c     $        pclc(ig)*fluxtop_lw(ig)
c
c             fluxsurf_sw(ig,1)=(1.-pclc(ig))*clsurf_sw(ig,1)+
c     $        pclc(ig)*fluxsurf_sw(ig,1)
c             fluxsurf_sw(ig,2)=(1.-pclc(ig))*clsurf_sw(ig,2)+
c     $        pclc(ig)*fluxsurf_sw(ig,2)
c
c             fluxtop_sw(ig,1)=(1.-pclc(ig))*cltop_sw(ig,1)+
c     $        pclc(ig)*fluxtop_sw(ig,1)
c             fluxtop_sw(ig,2)=(1.-pclc(ig))*cltop_sw(ig,2)+
c     $        pclc(ig)*fluxtop_sw(ig,2)
c           enddo
c
c           ELSE
c
c          Atmospheric water ice opacity and particle distribution:
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c           if (activice.and.naerkind.gt.1)
c     &      CALL h2oiceopacity(ngrid,nlayer,nq,pplay,pplev,pt,pq,
c     &      tau,aerosol,zls,co2ice)
c
c
c ----------
c  fin partie rajoutee par Franck (plus ENDIF ci-dessous)
c ----------


!          DO ig=1,ngrid
!            DO l=1,nlayer
!                   IF ( (pplev(ig,l+1) - pplev(ig,l) ) .gt. 0. )
!     .                PRINT*,'pb1 ',
!     .                        pplev(ig,l),
!     .                        pplev(ig,l+1),
!     .                        pplev(ig,l+1)-pplev(ig,l),
!     .                        l,
!     .                        ig
!            ENDDO
!         ENDDO

           CALL callradite(icount,ngrid,nlayer,aerosol,albedo,
     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
     $     zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,fluxtop_sw)

!           DO ig=1,ngrid
!            zdtlw(ig)=zdtlw(1)
!            zdtsw(ig)=zdtsw(1)
!            fluxsurf_lw(ig)=fluxsurf_lw(1)
!            fluxsurf_sw(ig,1)=fluxsurf_sw(1,1)
!            fluxsurf_sw(ig,2)=fluxsurf_sw(1,2)
!            fluxtop_lw(ig)=fluxtop_lw(1)
!            fluxtop_sw(ig,1)=fluxtop_sw(1,1) 
!            fluxtop_sw(ig,2)=fluxtop_sw(1,2)
!           ENDDO

c ----------
c           ENDIF ! end of condition on the cloudy fraction
c ----------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc
ccccc PARAM SLOPE : Insolation (direct + scattered)
ccccc
      DO ig=1,ngrid  
        sl_the = theta_sl(ig)
        IF (sl_the .ne. 0.) THEN
         ztim1=fluxsurf_sw(ig,1)+fluxsurf_sw(ig,2)
          DO l=1,2
           sl_lct = ptime*24. + 180.*long(ig)/pi/15.
           sl_ra  = pi*(1.0-sl_lct/12.)
           sl_lat = 180.*lati(ig)/pi
           sl_tau = tau(ig,1)
           sl_alb = albedo(ig,l)
           sl_psi = psi_sl(ig)
           sl_fl0 = fluxsurf_sw(ig,l)
           sl_di0 = 0.
           if (mu0(ig) .gt. 0.) then
            sl_di0 = mu0(ig)*(exp(-sl_tau/mu0(ig)))
            sl_di0 = sl_di0*1370./dist_sol/dist_sol
            sl_di0 = sl_di0/ztim1
            sl_di0 = fluxsurf_sw(ig,l)*sl_di0
           endif
           ! sait-on jamais (a cause des arrondis)
           if (sl_fl0 .lt. sl_di0) sl_di0=sl_fl0
     !!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL param_slope( mu0(ig), declin, sl_ra, sl_lat, 
     &            sl_tau, sl_alb, 
     &            sl_the, sl_psi, sl_di0, sl_fl0, sl_flu)
     !!!!!!!!!!!!!!!!!!!!!!!!!!
           fluxsurf_sw(ig,l) = sl_flu
                !!      sl_ls = 180.*zls/pi
                !!      sl_lct = ptime*24. + 180.*long(ig)/pi/15.
                !!      sl_lat = 180.*lati(ig)/pi
                !!      sl_tau = tau(ig,1)
                !!      sl_alb = albedo(ig,l)
                !!      sl_the = theta_sl(ig)
                !!      sl_psi = psi_sl(ig)
                !!      sl_fl0 = fluxsurf_sw(ig,l)
                !!      CALL param_slope_full(sl_ls, sl_lct, sl_lat, 
                !!     &                   sl_tau, sl_alb, 
                !!     &                   sl_the, sl_psi, sl_fl0, sl_flu)
          ENDDO
          !!! compute correction on IR flux as well
          sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
          fluxsurf_lw(ig)= fluxsurf_lw(ig)*sky
        ENDIF    
      ENDDO
ccccc
ccccc PARAM SLOPE
ccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


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
               fluxrad_save(ig)=emis(ig)*fluxsurf_lw(ig)
     $         +fluxsurf_sw(ig,1)*(1.-albedo(ig,1))
     $         +fluxsurf_sw(ig,2)*(1.-albedo(ig,2))

            !print*,'RAD ', fluxrad_save(ig)
            !print*,'LW ', emis(ig)*fluxsurf_lw(ig)
            !print*,'SW ', fluxsurf_sw(ig,2)*(1.-albedo(ig,2))

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

           !PRINT*,'zdtsw',zdtsw
           !PRINT*,'zdtlw',zdtlw
           !PRINT*,'zdtnirco2',zdtnirco2
           !PRINT*,'dtrad',dtrad


        ENDIF ! mod(icount-1,iradia).eq.0

c    Transformation of the radiative tendencies:
c    -----------------------------------------------

c          Net radiative surface flux (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,ngrid
               zplanck(ig)=tsurf(ig)*tsurf(ig)
               zplanck(ig)=emis(ig)*
     $         stephan*zplanck(ig)*zplanck(ig)
               fluxrad(ig)=fluxrad_save(ig)-zplanck(ig)
cccc
cccc param slope
cccc
               sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
               fluxrad(ig)=fluxrad(ig)+(1.-sky)*zplanck(ig)
cccc
cccc
cccc
           ENDDO


         DO l=1,nlayer
            DO ig=1,ngrid
               pdt(ig,l)=pdt(ig,l)+dtrad(ig,l)
            ENDDO
         ENDDO

      ENDIF



!c-----------------------------------------------------------------------
!c    3. Gravity wave and subgrid scale topography drag :
!c    -------------------------------------------------
!
!
!      IF(calllott)THEN
!
!        CALL calldrag_noro(ngrid,nlayer,ptimestep,
!     &                 pplay,pplev,pt,pu,pv,zdtgw,zdugw,zdvgw)
! 
!        DO l=1,nlayer
!          DO ig=1,ngrid
!            pdv(ig,l)=pdv(ig,l)+zdvgw(ig,l)
!            pdu(ig,l)=pdu(ig,l)+zdugw(ig,l)
!            pdt(ig,l)=pdt(ig,l)+zdtgw(ig,l)
!          ENDDO
!        ENDDO
!      ENDIF

c-----------------------------------------------------------------------
c    4. Vertical diffusion (turbulent mixing):
c    -----------------------------------------
c

      IF(calldifv) THEN

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

         DO ig=1,ngrid
          !! sensible heat flux in W/m2
          sensheat(ig) = zflubid(ig)-capcal(ig)*zdtsdif(ig)
          !! u star in similarity theory in m/s
          ustar(ig) = 0.4
     .               * sqrt( pu(ig,1)*pu(ig,1) + pv(ig,1)*pv(ig,1) )
     .               / log( 1.E+0 + zzlay(ig,1)/z0 )
         ENDDO   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! LES LES 
       IF (flag_LES) THEN        

         write (*,*) '************************************************' 
         write (*,*) '** LES mode: the difv part is only used to'
         write (*,*) '**  provide HFX and UST to the dynamics'
         write (*,*) '** NB: - dudif, dvdif, dhdif, dqdif are set to 0'
         write (*,*) '**     - tsurf is updated'     
         write (*,*) '************************************************'

!
          DO l=1,nlayer
            zdvdif(ig,l) = 0.
            zdudif(ig,l) = 0.
            zdhdif(ig,l) = 0.
            DO iq=1, nq
              zdqdif(ig,l,iq) = 0.
              zdqsdif(ig,iq) = 0. !! sortir de la boucle
            ENDDO 
          ENDDO
!
!         ENDDO
         !write (*,*) 'RAD ',fluxrad(igout)
         !write (*,*) 'GRD ',fluxgrd(igout)
         !write (*,*) 'dTs/dt ',capcal(igout)*zdtsurf(igout)
      ENDIF
!!! LES LES        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!!!gros caca : sans chaleur sensible
!         DO ig=1,ngrid
!            zdtsurf(ig)=zdtsurf(ig)+
!     s      (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
!         ENDDO


         if(tracer) then 
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
         end if

      ELSE    
         DO ig=1,ngrid
            zdtsurf(ig)=zdtsurf(ig)+
     s      (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (flag_LES) THEN
            write(*,*) 'LES mode !' 
            write(*,*) 'Please set calldifv to T in callphys.def'
            STOP
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF


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
      ENDIF

c-----------------------------------------------------------------------
c   6. Carbon dioxide condensation-sublimation:
c   -------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! get the actual co2 seasonal cap from Titus observations
      !!! (inherited from GCM only at first step, but not a big deal)
      IF (tituscap) THEN
         CALL geticecover( ngrid, 180.*zls/pi,
     .                  180.*long/pi, 180.*lati/pi, co2ice )
         co2ice(:) = co2ice(:)*10000.
         emis(:) = 0.95 ! so that points outside the cap are indeed at 0.95
                        ! avoid unwanted patchiness from GCM initial state
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(callcond) THEN
         CALL newcondens(ngrid,nlayer,nq,ptimestep,
     $              capcal,pplay,pplev,tsurf,pt,
     $              pphi,pdt,pdu,pdv,zdtsurf,pu,pv,pq,pdq,
     $              co2ice,albedo,emis,
     $              zdtc,zdtsurfc,pdpsrf,zduc,zdvc,zdqc,
     $	            fluxsurf_sw) 

         !! Titus cap security. 
         !! -- neglect pressure changes caused by sublimation/condensation
         IF (tituscap) THEN 
           DO ig=1,ngrid
             pdpsrf(ig) = 0.
           ENDDO
         ENDIF

         DO l=1,nlayer
           DO ig=1,ngrid
             pdt(ig,l)=pdt(ig,l)+zdtc(ig,l)
             pdv(ig,l)=pdv(ig,l)+zdvc(ig,l)
             pdu(ig,l)=pdu(ig,l)+zduc(ig,l)
           ENDDO
         ENDDO
         DO ig=1,ngrid
            zdtsurf(ig) = zdtsurf(ig) + zdtsurfc(ig)
            ps(ig)=pplev(ig,1) + pdpsrf(ig)*ptimestep
         ENDDO

         IF(tracer) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqc(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         END IF !(tracer)

      ENDIF  !(callcond) 

c        print*,'condens ok'
c-----------------------------------------------------------------------
c   7. Specific parameterizations for tracers 
c:   -----------------------------------------

      if(tracer) then 

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

         END IF  ! (water)

c   7b. Chemical species
c     ------------------

!c        --------------
!c        photochemistry :
!c        --------------
!         IF(photochem .or. thermochem) then
!          call calchim(ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,
!     .      zzlay,zday,pq,pdq,zdqchim,zdqschim,zdqcloud,zdqscloud)
!           
!c Photochemistry includes condensation of H2O2
!
!           do iq=nqchem_min,nq
!            if (noms(iq).eq."h2o2") then
!             DO l=1,nlayer
!               DO ig=1,ngrid
!                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqchim(ig,l,iq)
!                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqcloud(ig,l,iq)
!               ENDDO
!             ENDDO
!            else
!             DO l=1,nlayer
!               DO ig=1,ngrid
!                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqchim(ig,l,iq)
!               ENDDO
!             ENDDO
!            endif
!           ENDDO
!           do iq=nqchem_min,nq
!            if (noms(iq).eq."h2o2") then
!               DO ig=1,ngrid
!                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqschim(ig,iq)
!                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqscloud(ig,iq)
!               ENDDO
!            else
!               DO ig=1,ngrid
!                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqschim(ig,iq)
!               ENDDO
!            endif
!           ENDDO
!
!         END IF  ! (photochem.or.thermochem)
!c        print*,'photochem ok'

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
           endif  ! (test sur dustbin)

         END IF

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
         END IF   ! (sedimentation)

c        print*,'sedim ok'

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

      END IF    ! (tracer) 

c        print*,'tracers ok'


!c-----------------------------------------------------------------------
!c   8.5 THERMOSPHERE CALCULATION
!c-----------------------------------------------------------------------
!
!      if (callthermos) then
!        call thermosphere(pplev,pplay,dist_sol,
!     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
!     &     pt,pq,pu,pv,pdt,pdq,
!     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)
!c           do iq=nqchem_min,nq
!c           write(*,*) 'thermo iq,pq',iq,pq(690,1,iq)
!c           enddo
!
!        DO l=1,nlayer
!          DO ig=1,ngrid
!            dtrad(ig,l)=dtrad(ig,l)+zdteuv(ig,l)
!            pdt(ig,l)=pdt(ig,l)+zdtconduc(ig,l)
!     &                         +zdteuv(ig,l)
!            pdv(ig,l)=pdv(ig,l)+zdvmolvis(ig,l)
!            pdu(ig,l)=pdu(ig,l)+zdumolvis(ig,l)
!            DO iq=1, nq
!              pdq(ig,l,iq)=pdq(ig,l,iq)+zdqmoldiff(ig,l,iq)
!            ENDDO
!          ENDDO
!        ENDDO
!
!
!      endif

c-----------------------------------------------------------------------
c   8. Surface  and sub-surface soil temperature
c-----------------------------------------------------------------------
c
c
c   Surface temperature incrementation :
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO ig=1,ngrid
         tsurf(ig)=tsurf(ig)+ptimestep*zdtsurf(ig) 
      ENDDO


c  Prescription piege froid au pole sud (Except at high obliquity !!)
c  temperature en surface = temperature equilibre de phases du CO2
       
      IF (tracer.AND.water.AND.ngridmx.NE.1) THEN
         !if (caps.and.(obliquit.lt.27.)) then
         !  tsurf(ngrid)=1/(1/136.27-r/5.9e+5*alog(0.0095*ps(ngrid)))
         !endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! note WRF MESOSCALE AYMERIC -- mot cle "caps"
!!!!! watercaptag n'est plus utilise que dans vdifc
!!!!! ... pour que la sublimation ne soit pas stoppee 
!!!!! ... dans la calotte permanente nord si qsurf=0
!!!!! on desire garder cet effet regle par caps=T
!!!!! on a donc commente "if (caps.and.(obliquit.lt.27.))" ci-dessus
!!!!! --- remplacer ces lignes par qqch de plus approprie
!!!!!      si on s attaque a la calotte polaire sud
!!!!! pas d'autre occurrence majeure du mot-cle "caps"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


c        print*,'tracer, water and 3D ok'

c
c   soil temperatures and subsurface heat flux:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (callsoil) THEN
         CALL soil(ngrid,nsoilmx,.false.,inertiedat,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
      ENDIF

c        print*,'soil ok'





c-----------------------------------------------------------------------
c  9. Writing output files
c  ------------------------

c    -------------------------------
c    Dynamical fields incrementation
c    -------------------------------
c (FOR OUTPUT ONLY : the actual model integration is performed in the dynamics)

      DO l=1,nlayer
        DO ig=1,ngrid
          zt(ig,l)=pt(ig,l)  + pdt(ig,l)*ptimestep
          zu(ig,l)=pu(ig,l)  + pdu(ig,l)*ptimestep
          zv(ig,l)=pv(ig,l)  + pdv(ig,l)*ptimestep
        ENDDO
      ENDDO
      DO iq=1, nq
        DO l=1,nlayer
          DO ig=1,ngrid
            zq(ig,l,iq)=pq(ig,l,iq) +pdq(ig,l,iq)*ptimestep
          ENDDO
        ENDDO
      ENDDO
      DO ig=1,ngrid
          ps(ig)=pplev(ig,1) + pdpsrf(ig)*ptimestep !already in 7
      ENDDO
      DO l=1,nlayer
        DO ig=1,ngrid
             zplev(ig,l)=pplev(ig,l)/pplev(ig,1)*ps(ig)
             zplay(ig,l)=pplay(ig,l)/pplev(ig,1)*ps(ig)
        ENDDO
      ENDDO
      ! Density calculation
      DO l=1,nlayer
         DO ig=1,ngrid
            rho(ig,l) = zplay(ig,l)/(rnew(ig,l)*zt(ig,l))
         ENDDO
      ENDDO

c     Sum of fluxes in solar spectral bands (for output only)
      DO ig=1,ngrid
            fluxtop_sw_tot(ig)=fluxtop_sw(ig,1) + fluxtop_sw(ig,2)
            fluxsurf_sw_tot(ig)=fluxsurf_sw(ig,1) + fluxsurf_sw(ig,2)
      ENDDO
c ******* TEMPORAIRE ************************************************
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
        write(*,*) 'stability WARNING :' 
        write(*,*) 'pt, zt Tmin = ', pt(igmin,lmin), zt(igmin,lmin),
     &              'ig l =', igmin, lmin
      end if
c *******************************************************************

c     --------------------
c     Output on the screen 
c     --------------------


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
      ENDIF

cc****WRF
cc      IF (ngrid.NE.1) THEN
c	PRINT*,'check - values after update at ngrid/2+1'
c         WRITE(*,'(a4,a6,2a10)') 'l','u','v','T'
c         WRITE(*,'(i4,3f10.5)') (l,
c     s   zu(igout,l),
c     s   zv(igout,l),
c     s   zt(igout,l),
c     s   l=1,nlayer)
c
cc****WRF

         print*,'Ls =',zls*180./pi,
     &   ' tauref(700 Pa,lat=0) =',tauref(ngrid/2)
c     &   ' tau(Viking1) =',tau(ig_vl1,1)


c        -------------------------------------------------------------------
c        Writing NetCDF file  "RESTARTFI" at the end of the run
c        -------------------------------------------------------------------
c        Remarque : On  stocke restarfi 
c        juste avant que la dynamique ne le soit dans restart.
c        entre maintenant et l'ecriture de restart,
c        on aura itau = itau +1 et remise a jour de time.
c        (lastcall = .true. lorsque itau+1 = itaufin)
c        Donc on stocke  avec time = time + dtvr

         IF(lastcall) THEN
            ztime_fin = ptime + ptimestep/(float(iphysiq)*daysec) 
            write(*,*)'pour physdem ztime_fin =',ztime_fin
            call physdem1("restartfi.nc",long,lati,nsoilmx,nq,
     .              ptimestep,pday,
     .              ztime_fin,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .              area,albedodat,inertiedat,zmea,zstd,zsig,
     .              zgam,zthe)
         ENDIF


ccc**************** WRF OUTPUT **************************
ccc**************** WRF OUTPUT **************************
ccc**************** WRF OUTPUT **************************
      do ig=1,ngrid
         wtsurf(ig) = tsurf(ig)    !! surface temperature
         wco2ice(ig) = co2ice(ig)  !! co2 ice 

         !!! TEMP TEMP TEMP TEMP TEMP TEMP TEMP
         !!! specific to WRF WRF WRF
         !!! just to output water ice on surface
         !!! [it might not be water ice on surface but OK]
         !!! uncomment the Registry entry
         qsurflast(ig) = qsurf(ig,nqmx)

      enddo
c
c THIS INCLUDE IS AUTOMATICALLY GENERATED FROM REGISTRY
c
#include "fill_save.inc"
c
ccc**************** WRF OUTPUT **************************
ccc**************** WRF OUTPUT **************************
ccc**************** WRF OUTPUT **************************


cc-----------------------------------
cc you can still use meso_WRITEDIAGFI (e.g. for debugging purpose), 
cc though this is not the default strategy now
cc-----------------------------------
cc please use cudt in namelist.input to set frequency of outputs
cc----------------------------------- 
cc BEWARE: if at least one call to meso_WRITEDIAGFI is performed,
cc cudt cannot be 0 - otherwise you'll get a "Floating exception"
cc-----------------------------------         
!      call meso_WRITEDIAGFI(ngrid,"tauref",
!     .  "tauref","W.m-2",2,
!     .       tauref)
!      call meso_WRITEDIAGFI(ngrid,"dtrad",
!     .  "dtrad","W.m-2",2,
!     .       dtrad)
c      call meso_WRITEDIAGFI(ngrid,"tsurf",
c     .  "tsurf","K",2,
c     .       tsurf)
c
!      call meso_WRITEDIAGFI(ngrid,"zt",
!     .  "zt","W.m-2",3,
!     .       zt)
!      call meso_WRITEDIAGFI(ngrid,"zdtlw",
!     .  "zdtlw","W.m-2",2,
!     .       zdtlw)
!      call meso_WRITEDIAGFI(ngrid,"zdtsw",
!     .  "zdtsw","W.m-2",2,
!     .       zdtsw)


      icount=icount+1

!!!!!!!!!!!!!TEST TEST
!        ENDIF
!!!!!!!!!!!!!TEST TEST      

      write(*,*) 'now, back to the dynamical core...'
      RETURN
      END
