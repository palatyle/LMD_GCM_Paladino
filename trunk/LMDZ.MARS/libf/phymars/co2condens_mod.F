      MODULE co2condens_mod

      IMPLICIT NONE

      logical, save :: scavco2cond = .false. ! flag for using scavenging_by_co2
      
      CONTAINS

      SUBROUTINE co2condens(ngrid,nlayer,nq,ptimestep,
     $                  pcapcal,pplay,pplev,ptsrf,pt,
     $                  pphi,pdt,pdu,pdv,pdtsrf,pu,pv,pq,pdq,
     $                  piceco2,psolaralb,pemisurf,rdust,
     $                  pdtc,pdtsrfc,pdpsrf,pduc,pdvc,pdqc,
     $                  fluxsurf_sw,zls,
     $                  zdqssed_co2,pcondicea_co2microp,
     $                  zdtcloudco2,pdqsc)
                                                   
       use tracer_mod, only: noms, igcm_h2o_ice,
     &                      igcm_dust_mass, igcm_dust_number,
     &                      igcm_ccn_mass, igcm_ccn_number
       use surfdat_h, only: emissiv, phisfi
       use geometry_mod, only: latitude, ! grid point latitudes (rad)
     &                         longitude_deg, latitude_deg
       use planete_h, only: obliquit
       use comcstfi_h, only: cpp, g, r, pi 
       
#ifndef MESOSCALE
       USE vertical_layers_mod, ONLY: ap, bp
#endif
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Condensation/sublimation of CO2 ice on the ground and in the
c   atmosphere 
c   (Scheme described in Forget et al., Icarus, 1998)
c
c   author:   Francois Forget     1994-1996 ; updated 1996 -- 2018
c   ------
c             adapted to external CO2 ice clouds scheme by Deborah Bardet (2018) '
c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c
      include "callkeys.h"

c-----------------------------------------------------------------------
c    Arguments :
c    ---------
      INTEGER,INTENT(IN) :: ngrid  ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of vertical layers
      INTEGER,INTENT(IN) :: nq  ! number of tracers

      REAL,INTENT(IN) :: ptimestep ! physics timestep (s)
      REAL,INTENT(IN) :: pcapcal(ngrid)
      REAL,INTENT(IN) :: pplay(ngrid,nlayer) !mid-layer pressure (Pa)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: ptsrf(ngrid) ! surface temperature (K)
      REAL,INTENT(IN) :: pt(ngrid,nlayer) ! atmospheric temperature (K)
      REAL,INTENT(IN) :: pphi(ngrid,nlayer) ! geopotential (m2.s-2)
      REAL,INTENT(IN) :: pdt(ngrid,nlayer) ! tendency on temperature from
                                           ! previous physical processes (K/s)
      REAL,INTENT(IN) :: pdu(ngrid,nlayer) ! tendency on zonal wind (m/s2)
                                           ! from previous physical processes
      REAL,INTENT(IN) :: pdv(ngrid,nlayer) ! tendency on meridional wind (m/s2)
                                           ! from previous physical processes
      REAL,INTENT(IN) :: pdtsrf(ngrid) ! tendency on surface temperature from
                                       ! previous physical processes (K/s)
      REAL,INTENT(IN) :: pu(ngrid,nlayer) ! zonal wind (m/s)
      REAL,INTENT(IN) :: pv(ngrid,nlayer) ! meridional wind (m/s)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq) ! tracers (../kg_air)
      REAL,INTENT(IN) :: pdq(ngrid,nlayer,nq) ! tendency on tracers from
                                              ! previous physical processes

      REAL,INTENT(IN) :: zdqssed_co2(ngrid)  ! CO2 flux at the surface  (kg.m-2.s-1)
      REAL,INTENT(IN) :: pcondicea_co2microp(ngrid,nlayer)! tendency due to CO2 condensation (kg/kg.s-1)
      REAL,INTENT(IN) :: zdtcloudco2(ngrid,nlayer) ! tendency on temperature due to latent heat                 

      REAL,INTENT(INOUT) :: piceco2(ngrid) ! CO2 ice on the surface (kg.m-2)
      REAL,INTENT(INOUT) :: psolaralb(ngrid,2) ! albedo of the surface
      REAL,INTENT(INOUT) :: pemisurf(ngrid) ! emissivity of the surface
      REAL,INTENT(IN) :: rdust(ngrid,nlayer) ! dust effective radius
      
      ! tendencies due to CO2 condensation/sublimation:
      REAL,INTENT(OUT) :: pdtc(ngrid,nlayer) ! tendency on temperature (K/s)
      REAL,INTENT(OUT) :: pdtsrfc(ngrid) ! tendency on surface temperature (K/s)
      REAL,INTENT(OUT) :: pdpsrf(ngrid) ! tendency on surface pressure (Pa/s)
      REAL,INTENT(OUT) :: pduc(ngrid,nlayer) ! tendency on zonal wind (m.s-2)
      REAL,INTENT(OUT) :: pdvc(ngrid,nlayer) ! tendency on meridional wind (m.s-2)
      REAL,INTENT(OUT) :: pdqc(ngrid,nlayer,nq) ! tendency on tracers
      REAL,INTENT(OUT) :: pdqsc(ngrid,nq) ! tendency on surface tracers
      
      ! added to calculate flux dependent albedo:
      REAL,intent(in) :: fluxsurf_sw(ngrid,2)
      real,intent(in) :: zls ! solar longitude (rad)

c
c    Local variables :
c    -----------------

      INTEGER i,j
      INTEGER l,ig,iq,icap
      REAL zt(ngrid,nlayer)
      REAL zcpi
      REAL ztcond (ngrid,nlayer+1) ! CO2 condensation temperature (atm)
      REAL ztcondsol(ngrid) ! CO2 condensation temperature (surface)
      REAL zdiceco2(ngrid)
      REAL zcondicea(ngrid,nlayer) ! condensation rate in layer  l  (kg/m2/s)
      REAL zcondices(ngrid) ! condensation rate on the ground (kg/m2/s)
      REAL zfallice(ngrid,nlayer+1) ! amount of ice falling from layer l (kg/m2/s)
      REAL zfallheat
      REAL zmflux(nlayer+1)
      REAL zu(nlayer),zv(nlayer)
      REAL zqc(nlayer,nq),zq1(nlayer)
      REAL ztsrf(ngrid)
      REAL ztc(nlayer), ztm(nlayer+1) 
      REAL zum(nlayer+1) , zvm(nlayer+1)
      REAL zqm(nlayer+1,nq),zqm1(nlayer+1)
      REAL masse(nlayer),w(nlayer+1)
      REAL Sm(nlayer),Smq(nlayer,nq),mixmas,qmix
      REAL availco2
      LOGICAL condsub(ngrid)

      real :: emisref(ngrid)
      
      REAL zdq_scav(ngrid,nlayer,nq) ! tendancy due to scavenging by co2
      REAL zq(ngrid,nlayer,nq)

c variable speciale diagnostique
      real tconda1(ngrid,nlayer)
      real tconda2(ngrid,nlayer)
c     REAL zdiceco2a(ngrid) ! for diagnostic only
      real zdtsig (ngrid,nlayer)
      real zdt (ngrid,nlayer)
      real vmr_co2(ngrid,nlayer) ! co2 volume mixing ratio
! improved_ztcond flag: If set to .true. (AND running with a 'co2' tracer)
! then condensation temperature is computed using partial pressure of CO2
      logical,parameter :: improved_ztcond=.true.

c   local saved variables
      integer,save :: ico2 ! index of CO2 tracer
      real,save :: qco2,mmean
      real,parameter :: latcond=5.9e5 ! (J/kg) Latent heat of solid CO2 ice
      real,parameter :: tcond1mb=136.27 ! condensation temperature (K) at 1 mbar
      real,parameter :: cpice=1000. ! (J.kg-1.K-1) specific heat of CO2 ice
      REAL,SAVE :: acond,bcond,ccond
      real,save :: m_co2, m_noco2, A , B

      LOGICAL,SAVE :: firstcall = .true. !,firstcall2=.true.

c D.BARDET: for debug 
      real ztc3D(ngrid,nlayer)
      REAL ztm3D(ngrid,nlayer)
      REAL zmflux3D(ngrid,nlayer)
c----------------------------------------------------------------------

c   Initialisation
c   --------------
c
      ! AS: firstcall OK absolute
      IF (firstcall) THEN
         
         bcond=1./tcond1mb
         ccond=cpp/(g*latcond)
         acond=r/latcond

         firstcall=.false.
         write(*,*) 'CO2condens: improved_ztcond=',improved_ztcond
         PRINT*,'In co2condens:Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
         PRINT*,'acond,bcond,ccond',acond,bcond,ccond

         ico2=0

         if (tracer) then
c          Prepare Special treatment if one of the tracer is CO2 gas
           do iq=1,nq
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)   
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)   
c               Compute A and B coefficient use to compute
c               mean molecular mass Mair defined by
c               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
c               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             endif
           enddo
         end if
      ENDIF ! of IF (firstcall)
      zcpi=1./cpp

c
c======================================================================
c    Calcul of CO2 condensation sublimation 
c    ============================================================
c  
c    Used variable :
c       piceco2(ngrid)   :  amount of co2 ice on the ground (kg/m2)
c       zcondicea(ngrid,l):  condensation rate in layer  l  (kg/m2/s)
c       zcondices(ngrid):  condensation rate on the ground (kg/m2/s)
c       zfallice(ngrid,l):amount of ice falling from layer l (kg/m2/s)
c                           
c       pdtc(ngrid,nlayer)  : dT/dt due to cond/sub
c
c
c     Tendencies set to 0
c     -------------------------------------
      zcondicea(1:ngrid,1:nlayer) = 0.
      zfallice(1:ngrid,1:nlayer+1) = 0.
      pduc(1:ngrid,1:nlayer)  = 0
      pdvc(1:ngrid,1:nlayer)  = 0
      pdtc(1:ngrid,1:nlayer) = 0.
      pdqsc(1:ngrid,1:nq) = 0
         
      pdqc(1:ngrid,1:nlayer,1:nq)  = 0

      zcondices(1:ngrid) = 0.
      pdtsrfc(1:ngrid) = 0.
      pdpsrf(1:ngrid) = 0.
      condsub(1:ngrid) = .false.
      zdiceco2(1:ngrid) = 0.

      zfallheat=0
      
      zdq_scav(:,:,:)=0.

c     Update tendencies from previous processes
c     -------------------------------------
      DO l=1,nlayer
         DO ig=1,ngrid
            zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep
	    do iq=1,nq
             zq(ig,l,iq)=pq(ig,l,iq)+pdq(ig,l,iq)*ptimestep
            enddo
         ENDDO
      ENDDO
      
c     *************************
c     ATMOSPHERIC CONDENSATION
c     *************************

c     Compute CO2 Volume mixing ratio
c     -------------------------------
      if (improved_ztcond.and.(ico2.ne.0)) then
         DO l=1,nlayer
            DO ig=1,ngrid
              qco2=pq(ig,l,ico2)+pdq(ig,l,ico2)*ptimestep
c             Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
              mmean=1/(A*qco2 +B)
              vmr_co2(ig,l) = qco2*mmean/m_co2 
            ENDDO
         ENDDO
      else
         DO l=1,nlayer
            DO ig=1,ngrid
              vmr_co2(ig,l)=0.95
            ENDDO
         ENDDO
      endif

      IF (.NOT. co2clouds) then
c     forecast of atmospheric temperature zt and frost temperature ztcond
c     --------------------------------------------------------------------

      DO l=1,nlayer
         DO ig=1,ngrid
!            ztcond(ig,l)=1./(bcond-acond*log(.0095*pplay(ig,l)))
            if (pplay(ig,l).ge.1e-4) then
              ztcond(ig,l)=
     &         1./(bcond-acond*log(.01*vmr_co2(ig,l)*pplay(ig,l)))
            else
              ztcond(ig,l)=0.0 !mars Monica
            endif
         ENDDO
      ENDDO

      ztcond(:,nlayer+1)=ztcond(:,nlayer)
 
c      Condensation/sublimation in the atmosphere
c      ------------------------------------------
c      (calcul of zcondicea , zfallice and pdtc)
c
      DO l=nlayer , 1, -1
         DO ig=1,ngrid
           pdtc(ig,l)=0.
           IF((zt(ig,l).LT.ztcond(ig,l)).or.(zfallice(ig,l+1).gt.0))THEN
               condsub(ig)=.true.
               IF (zfallice(ig,l+1).gt.0) then  
                 zfallheat=zfallice(ig,l+1)*
     &           (pphi(ig,l+1)-pphi(ig,l) +
     &          cpice*(ztcond(ig,l+1)-ztcond(ig,l)))/latcond
               ELSE
                   zfallheat=0.
               ENDIF
               pdtc(ig,l)=(ztcond(ig,l) - zt(ig,l))/ptimestep
               zcondicea(ig,l)=(pplev(ig,l)-pplev(ig,l+1))
     &                        *ccond*pdtc(ig,l)- zfallheat
c              Case when the ice from above sublimes entirely
c              """""""""""""""""""""""""""""""""""""""""""""""
               IF (zfallice(ig,l+1).lt.- zcondicea(ig,l)) then
                  pdtc(ig,l)=(-zfallice(ig,l+1)+zfallheat)/
     &              (ccond*(pplev(ig,l)-pplev(ig,l+1)))
                  zcondicea(ig,l)= -zfallice(ig,l+1)
               END IF

               zfallice(ig,l) = zcondicea(ig,l)+zfallice(ig,l+1)
            END IF
         ENDDO
      ENDDO
      
      if (scavco2cond) then
        call scavenging_by_co2(ngrid,nlayer,nq,ptimestep,pplev,zq,
     &                       rdust,zcondicea,zfallice,zdq_scav,pdqsc)
      endif
      
      ELSE ! if co2 clouds
        DO ig=1,ngrid 
         zfallice(ig,1) = zdqssed_co2(ig)
        ENDDO
        DO l=nlayer , 1, -1
            DO ig=1,ngrid 
         zcondicea(ig,l) = pcondicea_co2microp(ig,l)* 
     &                        (pplev(ig,l) - pplev(ig,l+1))/g
            ENDDO
        ENDDO
 
      ENDIF ! end of if co2clouds

      call WRITEdiagfi(ngrid,"pdtc_atm",
     &         "temperature tendency due to CO2 condensation",
     &         " ",3,pdtc) 
     
       call WRITEdiagfi(ngrid,"zcondicea",
     &         "",
     &         " ",3,zcondicea)         
      
       call WRITEdiagfi(ngrid,"zfallice",
     &         "",
     &         " ",2,zfallice(ngrid,1))
     

c     *************************
c       SURFACE  CONDENSATION
c     *************************

c     forecast of ground temperature ztsrf and frost temperature ztcondsol
c     --------------------------------------------------------------------
      DO ig=1,ngrid
         ztcondsol(ig)=
     &          1./(bcond-acond*log(.01*vmr_co2(ig,1)*pplev(ig,1)))
         ztsrf(ig) = ptsrf(ig) + pdtsrf(ig)*ptimestep
      ENDDO

c
c      Condensation/sublimation on the ground 
c      --------------------------------------
c      (compute zcondices and pdtsrfc)
c
      DO ig=1,ngrid
         IF(latitude(ig).lt.0) THEN
           ! Southern hemisphere
            icap=2
         ELSE
           ! Northern hemisphere
            icap=1
         ENDIF
       
c
c        Loop on where we have  condensation/ sublimation 
         IF ((ztsrf(ig) .LT. ztcondsol(ig)) .OR.   ! ground cond 
     $       (zfallice(ig,1).NE.0.) .OR.           ! falling snow
     $      ((ztsrf(ig) .GT. ztcondsol(ig)) .AND.  ! ground sublim.
     $      ((piceco2(ig)+zfallice(ig,1)*ptimestep) .NE. 0.))) THEN
            condsub(ig) = .true. 

            IF (zfallice(ig,1).gt.0) then  
                 zfallheat=zfallice(ig,1)*
     &           (pphi(ig,1)- phisfi(ig) +
     &           cpice*(ztcond(ig,1)-ztcondsol(ig)))/latcond
            ELSE
                 zfallheat=0.
            ENDIF

c           condensation or partial sublimation of CO2 ice
c           """""""""""""""""""""""""""""""""""""""""""""""
            zcondices(ig)=pcapcal(ig)*(ztcondsol(ig)-ztsrf(ig)) 
     &      /(latcond*ptimestep)         - zfallheat
            pdtsrfc(ig) = (ztcondsol(ig) - ztsrf(ig))/ptimestep
            zdiceco2(ig) = zcondices(ig) + zfallice(ig,1)

#ifdef MESOSCALE
      print*, "not enough CO2 tracer in 1st layer to condense"
      print*, ">>> to be implemented in the mesoscale case"
      print*, "because this uses ap levels..."
#else
c           If there is not enough CO2 tracer in 1st layer to condense
c           """"""""""""""""""""""""""""""""""""""""""""""""""""""
            IF(ico2.ne.0) then
c              Available CO2 tracer in layer 1 at end of timestep (kg/m2)
                availco2= pq(ig,1,ico2)*((ap(1)-ap(2))+
     &          (bp(1)-bp(2))*(pplev(ig,1)/g-zdiceco2(ig)*ptimestep))

               IF ((zcondices(ig) + zcondicea(ig,1))*ptimestep
     &           .gt.availco2) then
                   zcondices(ig) = availco2/ptimestep -zcondicea(ig,1)
                   zdiceco2(ig) = zcondices(ig) + zfallice(ig,1)
                   pdtsrfc(ig)=(latcond/pcapcal(ig))*
     &                          (zcondices(ig)+zfallheat)
               ENDIF 
            ENDIF 
#endif

c           If the entire CO2 ice layer sublimes
c           """"""""""""""""""""""""""""""""""""""""""""""""""""
c           (including what has just condensed in the atmosphere)

            IF((piceco2(ig)/ptimestep+zfallice(ig,1)).LE.
     &          -zcondices(ig))THEN
               zcondices(ig) = -piceco2(ig)/ptimestep - zfallice(ig,1) 
               pdtsrfc(ig)=(latcond/pcapcal(ig))*
     &         (zcondices(ig)+zfallheat)
               zdiceco2(ig) = zcondices(ig) + zfallice(ig,1)
            END IF

c           Changing CO2 ice amount and pressure :
c           """"""""""""""""""""""""""""""""""""

            piceco2(ig) = piceco2(ig) + zdiceco2(ig)*ptimestep
            pdpsrf(ig) = -zdiceco2(ig)*g

            IF(ABS(pdpsrf(ig)*ptimestep).GT.pplev(ig,1)) THEN
               PRINT*,'STOP in condens'
               PRINT*,'condensing more than total mass'
               PRINT*,'Grid point ',ig
               PRINT*,'Longitude(degrees): ',longitude_deg(ig)
               PRINT*,'Latitude(degrees): ',latitude_deg(ig)
               PRINT*,'Ps = ',pplev(ig,1)
               PRINT*,'d Ps = ',pdpsrf(ig)
               STOP
            ENDIF
         END IF ! if there is condensation/sublimation
      ENDDO ! of DO ig=1,ngrid

c  ********************************************************************
c  Surface albedo and emissivity of the surface below the snow (emisref)
c  ********************************************************************

!     Check that amont of CO2 ice is not problematic
      DO ig=1,ngrid
           if(.not.piceco2(ig).ge.0.) THEN
              if(piceco2(ig).le.-5.e-8) print*,
     $         'WARNING co2condens piceco2(',ig,')=', piceco2(ig)
               piceco2(ig)=0.
           endif
      ENDDO
      
!     Set albedo and emissivity of the surface
!     ----------------------------------------
      CALL albedocaps(zls,ngrid,piceco2,psolaralb,emisref)

! set pemisurf() to emissiv when there is bare surface (needed for co2snow)
      DO ig=1,ngrid
        if (piceco2(ig).eq.0) then
          pemisurf(ig)=emissiv
        endif
      ENDDO

!         firstcall2=.false.
c ***************************************************************
c  Correction to account for redistribution between sigma or hybrid 
c  layers when changing surface pressure (and warming/cooling
c  of the CO2 currently changing phase).
c *************************************************************

      DO ig=1,ngrid
        if (condsub(ig)) then
           do l=1,nlayer
             ztc(l)  =zt(ig,l)   +pdtc(ig,l)  *ptimestep
             zu(l)   =pu(ig,l)   +pdu( ig,l)  *ptimestep
             zv(l)   =pv(ig,l)   +pdv( ig,l)  *ptimestep
            do iq=1,nq
             zqc(l,iq)=zq(ig,l,iq)+zdq_scav(ig,l,iq)*ptimestep ! zdq_scav=0 if watercloud=false
            enddo
           enddo

c  Mass fluxes through the sigma levels (kg.m-2.s-1)  (>0 when up)
c  """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

            zmflux(1) = -zcondices(ig)
            DO l=1,nlayer
             zmflux(l+1) = zmflux(l) -zcondicea(ig,l)
#ifndef MESOSCALE
     &        + (bp(l)-bp(l+1))*(zfallice(ig,1)-zmflux(1))
c zmflux set to 0 if very low to avoid: top layer is disappearing in v1ld  
          if (abs(zmflux(l+1)).lt.1E-13.OR.bp(l+1).eq.0.) zmflux(l+1)=0.
#else
          if (abs(zmflux(l+1)).lt.1E-13) zmflux(l+1)=0.
#endif
            END DO

#ifdef MESOSCALE
         print*,"absurd mass set because bp not available"
         print*,"TO BE FIXED"
#else
c Mass of each layer at the end of timestep
c ----------------------------------------- 
            DO l=1,nlayer
              masse(l)=( pplev(ig,l) - pplev(ig,l+1) +  
     &                 (bp(l)-bp(l+1))*pdpsrf(ig)*ptimestep)/g
            END DO
#endif

c  Corresponding fluxes for T,U,V,Q
c  """"""""""""""""""""""""""""""""

c           averaging operator for TRANSPORT  
c           """"""""""""""""""""""""""""""""
c           Value transfert at the surface interface when condensation
c           sublimation:
            ztm(1) = ztsrf(ig) + pdtsrfc(ig)*ptimestep
            zum(1) = 0  
            zvm(1) = 0  
            do iq=1,nq
              zqm(1,iq)=0. ! most tracer do not condense !
            enddo
c           Special case if one of the tracer is CO2 gas
            if (ico2.ne.0) zqm(1,ico2)=1. ! flux is 100% CO2

c           Van Leer scheme:
            DO l=1,nlayer+1
                w(l)=-zmflux(l)*ptimestep
            END DO
            call vl1d(nlayer,ztc,2.,masse,w,ztm) 
            call vl1d(nlayer,zu ,2.,masse,w,zum) 
            call vl1d(nlayer,zv ,2.,masse,w,zvm) 
            do iq=1,nq
             do l=1,nlayer
              zq1(l)=zqc(l,iq)
             enddo
             zqm1(1)=zqm(1,iq)
             call vl1d(nlayer,zq1,2.,masse,w,zqm1)
             do l=2,nlayer
              zqc(l,iq)=zq1(l)
              zqm(l,iq)=zqm1(l)
             enddo
            enddo

c           Surface condensation affects low winds
            if (zmflux(1).lt.0) then 
                zum(1)= zu(1) *  (w(1)/masse(1))
                zvm(1)= zv(1) *  (w(1)/masse(1))
                if (w(1).gt.masse(1)) then ! ensure numerical stability
                  zum(1)= (zu(1)-zum(2))*masse(1)/w(1) + zum(2)
                  zvm(1)= (zv(1)-zvm(2))*masse(1)/w(1) + zvm(2)
                end if
            end if
                    
            ztm(nlayer+1)= ztc(nlayer) ! should not be used, but... 
            zum(nlayer+1)= zu(nlayer)  ! should not be used, but...
            zvm(nlayer+1)= zv(nlayer)  ! should not be used, but...
            do iq=1,nq
             zqm(nlayer+1,iq)= zqc(nlayer,iq)
            enddo

#ifdef MESOSCALE
!!!! AS: This part must be commented in the mesoscale model
!!!! AS: ... to avoid instabilities.
!!!! AS: you have to compile with -DMESOSCALE to do so
#else  
c           Tendencies on T, U, V, Q 
c           """"""""""""""""""""""""
            DO l=1,nlayer
               IF(.not. co2clouds) THEN
c             Tendencies on T
                zdtsig(ig,l) = (1/masse(l)) *
     &        ( zmflux(l)*(ztm(l) - ztc(l)) 
     &        - zmflux(l+1)*(ztm(l+1) - ztc(l))
     &        + zcondicea(ig,l)*(ztcond(ig,l)-ztc(l))  )
               ELSE
                zdtsig(ig,l) = (1/masse(l)) *
     &        ( zmflux(l)*(ztm(l) - ztc(l)) 
     &        - zmflux(l+1)*(ztm(l+1) - ztc(l)))               
               ENDIF
c D.BARDET: for diagnotics 
                zmflux3D(ig,l)=zmflux(l)
                ztm3D(ig,l)=ztm(l)
                ztc3D(ig,l)=ztc(l)
                
                pdtc(ig,l) =  pdtc(ig,l) + zdtsig(ig,l)

c             Tendencies on U
                pduc(ig,l)   = (1/masse(l)) *
     &        ( zmflux(l)*(zum(l) - zu(l))
     &        - zmflux(l+1)*(zum(l+1) - zu(l)) )


c             Tendencies on V
                pdvc(ig,l)   = (1/masse(l)) *
     &        ( zmflux(l)*(zvm(l) - zv(l))
     &        - zmflux(l+1)*(zvm(l+1) - zv(l)) )

            END DO

#endif

c           Tendencies on Q
            do iq=1,nq
!              if (noms(iq).eq.'co2') then 
              if (iq.eq.ico2) then 
c               SPECIAL Case when the tracer IS CO2 :
                DO l=1,nlayer
                  pdqc(ig,l,iq)= (1/masse(l)) *
     &           ( zmflux(l)*(zqm(l,iq) - zqc(l,iq))
     &           - zmflux(l+1)*(zqm(l+1,iq) - zqc(l,iq)) 
     &           + zcondicea(ig,l)*(zqc(l,iq)-1.) )
                END DO
              else
                DO l=1,nlayer
                  pdqc(ig,l,iq)= (1/masse(l)) *
     &           ( zmflux(l)*(zqm(l,iq) - zqc(l,iq))
     &           - zmflux(l+1)*(zqm(l+1,iq) - zqc(l,iq)) 
     &           + zcondicea(ig,l)*zqc(l,iq) )
     
                  pdqc(ig,l,iq)=pdqc(ig,l,iq)+zdq_scav(ig,l,iq) ! zdq_scav=0 if watercloud=false
                END DO
              end if
            enddo

       end if ! if (condsub)
      END DO  ! loop on ig 

c ***************************************************************
c   CO2 snow / clouds scheme
c ***************************************************************

      call co2snow(ngrid,nlayer,ptimestep,emisref,condsub,pplev,
     &        zcondicea,zcondices,zfallice,pemisurf)

c ***************************************************************
c Ecriture des diagnostiques
c ***************************************************************

c     DO l=1,nlayer
c        DO ig=1,ngrid
c         Taux de cond en kg.m-2.pa-1.s-1
c          tconda1(ig,l)=zcondicea(ig,l)/(pplev(ig,l)-pplev(ig,l+1))
c          Taux de cond en kg.m-3.s-1
c          tconda2(ig,l)=tconda1(ig,l)*pplay(ig,l)*g/(r*pt(ig,l))
c        END DO
c     END DO
c     call WRITEDIAGFI(ngrid,'tconda1',
c    &'Taux de condensation CO2 atmospherique /Pa',
c    & 'kg.m-2.Pa-1.s-1',3,tconda1)
c     call WRITEDIAGFI(ngrid,'tconda2',
c    &'Taux de condensation CO2 atmospherique /m',
c    & 'kg.m-3.s-1',3,tconda2)

! output falling co2 ice in 1st layer:
!      call WRITEDIAGFI(ngrid,'fallice',
!     &'Precipitation of co2 ice',
!     & 'kg.m-2.s-1',2,zfallice(1,1))

#ifndef MESOSCALE 
! Extra special case for surface temperature tendency pdtsrfc:
! we want to fix the south pole temperature to CO2 condensation temperature
         if (caps.and.(obliquit.lt.27.)) then
           ! check if last grid point is the south pole
           if (abs(latitude(ngrid)-(-pi/2.)).lt.1.e-5) then
           ! NB: Updated surface pressure, at grid point 'ngrid', is
           !     ps(ngrid)=pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep
!             write(*,*) "co2condens: South pole: latitude(ngrid)=",
!     &                                           latitude(ngrid)
             ztcondsol(ngrid)=
     &          1./(bcond-acond*log(.01*vmr_co2(ngrid,1)*
     &                    (pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep)))
             pdtsrfc(ngrid)=(ztcondsol(ngrid)-ztsrf(ngrid))/ptimestep
           endif
         endif
#endif

      END SUBROUTINE co2condens



c *****************************************************************
      SUBROUTINE vl1d(nlayer,q,pente_max,masse,w,qm)
c
c    
c     Operateur de moyenne inter-couche pour calcul de transport type
c     Van-Leer " pseudo amont " dans la verticale
c    q,w sont des arguments d'entree  pour le s-pg ....
c    masse : masse de la couche Dp/g
c    w : masse d'atm ``transferee'' a chaque pas de temps (kg.m-2)
c    pente_max = 2 conseillee
c
c
c   --------------------------------------------------------------------
      IMPLICIT NONE

c
c
c
c   Arguments:
c   ----------
      integer nlayer
      real masse(nlayer),pente_max
      REAL q(nlayer),qm(nlayer+1)
      REAL w(nlayer+1)
c
c      Local 
c   ---------
c
      INTEGER l
c
      real dzq(nlayer),dzqw(nlayer),adzqw(nlayer),dzqmax
      real sigw, Mtot, MQtot
      integer m 
c     integer ismax,ismin 


c    On oriente tout dans le sens de la pression 
c     W > 0 WHEN DOWN !!!!!!!!!!!!!

      do l=2,nlayer
            dzqw(l)=q(l-1)-q(l)
            adzqw(l)=abs(dzqw(l))
      enddo

      do l=2,nlayer-1
            if(dzqw(l)*dzqw(l+1).gt.0.) then
                dzq(l)=0.5*(dzqw(l)+dzqw(l+1))
            else
                dzq(l)=0.
            endif
            dzqmax=pente_max*min(adzqw(l),adzqw(l+1))
            dzq(l)=sign(min(abs(dzq(l)),dzqmax),dzq(l))
      enddo

         dzq(1)=0.
         dzq(nlayer)=0.

       do l = 1,nlayer-1

c         Regular scheme (transfered mass < layer mass)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(l+1).gt.0. .and. w(l+1).le.masse(l+1)) then
             sigw=w(l+1)/masse(l+1)
             qm(l+1)=(q(l+1)+0.5*(1.-sigw)*dzq(l+1))
          else if(w(l+1).le.0. .and. -w(l+1).le.masse(l)) then
             sigw=w(l+1)/masse(l)
             qm(l+1)=(q(l)-0.5*(1.+sigw)*dzq(l))

c         Extended scheme (transfered mass > layer mass)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if(w(l+1).gt.0.) then
             m=l+1
             Mtot = masse(m)
             MQtot = masse(m)*q(m)
             do while ((m.lt.nlayer).and.(w(l+1).gt.(Mtot+masse(m+1))))
                m=m+1
                Mtot = Mtot + masse(m)
                MQtot = MQtot + masse(m)*q(m)
             end do
             if (m.lt.nlayer) then
                sigw=(w(l+1)-Mtot)/masse(m+1)
                qm(l+1)= (1/w(l+1))*(MQtot + (w(l+1)-Mtot)*
     &          (q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
             else
                w(l+1) = Mtot
                qm(l+1) = Mqtot / Mtot
                write(*,*) 'top layer is disapearing !'
                stop
             end if
          else      ! if(w(l+1).lt.0) 
             m = l-1 
             Mtot = masse(m+1)
             MQtot = masse(m+1)*q(m+1)
             if (m.gt.0) then ! because some compilers will have problems
                              ! evaluating masse(0)
              do while ((m.gt.0).and.(-w(l+1).gt.(Mtot+masse(m))))
                m=m-1
                Mtot = Mtot + masse(m+1)
                MQtot = MQtot + masse(m+1)*q(m+1)
                if (m.eq.0) exit
              end do
             endif
             if (m.gt.0) then
                sigw=(w(l+1)+Mtot)/masse(m)
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)* 
     &          (q(m)-0.5*(1.+sigw)*dzq(m))  )
             else
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)*qm(1))
             end if   
          end if
       enddo

c     boundary conditions (not used in co2condens !!)
c         qm(nlayer+1)=0.
c         if(w(1).gt.0.) then
c            qm(1)=q(1)
c         else 
c           qm(1)=0.
c         end if

      END SUBROUTINE vl1d
          
c *****************************************************************
      SUBROUTINE scavenging_by_co2(ngrid,nlayer,nq,ptimestep,pplev,pq,
     &                          rdust,pcondicea,pfallice,pdq_scav,pdqsc)
                      
c
c    
c     Calcul de la quantite de poussiere lessivee par les nuages de CO2
c     
c   --------------------------------------------------------------------
      use tracer_mod, only: nqmx, igcm_h2o_vap, igcm_h2o_ice,
     &                      igcm_dust_mass, igcm_dust_number,
     &                      igcm_ccn_mass, igcm_ccn_number,
     &                      rho_dust, nuice_sed, nuice_ref,r3n_q
      use comcstfi_h, only: g
      
      IMPLICIT NONE
      include "callkeys.h" ! for the flags water, microphys and freedust
c
c
c   Arguments:
      INTEGER,INTENT(IN) :: ngrid  ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of vertical layers
      INTEGER,INTENT(IN) :: nq  ! number of tracers
      REAL,INTENT(IN) :: ptimestep ! physics timestep (s)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq)
      REAL,INTENT(IN) :: rdust(ngrid,nlayer) ! dust effective radius
      REAL,INTENT(IN) :: pcondicea(ngrid,nlayer) ! condensation rate in layer  l  (kg/m2/s)
      REAL,INTENT(IN) :: pfallice(ngrid,nlayer+1) ! amount of ice falling from layer l (kg/m2/s)
      
      REAL,INTENT(OUT) :: pdq_scav(ngrid,nlayer,nq) ! tendancy due to scavenging by co2
      REAL,INTENT(OUT) :: pdqsc(ngrid,nq) ! tendency on surface tracers
      
c   Locals:
      INTEGER l,ig
      REAL scav_ratio_dust, scav_ratio_wice ! ratio of the dust/water ice mass mixing ratios in condensing CO2 ice and in air
      REAL scav_dust_mass(nlayer+1) ! dust flux (mass) scavenged towards the lower layer (kg/m2/s) (POSITIVE WHEN DOWNWARD)
      REAL scav_dust_number(nlayer+1) ! dust flux (number) scavenged towards the lower layer (kg/m2/s) (POSITIVE WHEN DOWNWARD)
      REAL scav_ccn_mass(nlayer+1) ! ccn flux (mass) scavenged towards the lower layer
      REAL scav_ccn_number(nlayer+1) ! ccn flux (number) scavenged towards the lower layer
      REAL scav_h2o_ice(nlayer+1) ! water ice flux (mass) scavenged towards the lower layer
      REAL massl ! mass of the layer l at point ig (kg/m2)
      
c   Initialization:
      scav_ratio_dust = 100 !1 !10 !100 !1000
      scav_ratio_wice = scav_ratio_dust
      pdq_scav(:,:,:)=0.
      
      DO ig=1,ngrid
        scav_dust_mass(nlayer+1)=0.
        scav_dust_number(nlayer+1)=0.
        scav_ccn_mass(nlayer+1)=0.
        scav_ccn_number(nlayer+1)=0.
        scav_h2o_ice(nlayer+1)=0.
	
        DO l=nlayer , 1, -1
	 massl=(pplev(ig,l)-pplev(ig,l+1))/g
         IF(pcondicea(ig,l).GT.0.)THEN ! if CO2 condenses and traps dust/water ice
           ! Calculation of the tendencies
	   if (freedust) then
	     pdq_scav(ig,l,igcm_dust_mass)=-pq(ig,l,igcm_dust_mass)
     &                                     /ptimestep*(1-exp(
     &              -scav_ratio_dust*pcondicea(ig,l)*ptimestep/massl))
             
	     pdq_scav(ig,l,igcm_dust_number)=pdq_scav(ig,l,igcm_dust_mass)
     &                                    *r3n_q/rdust(ig,l)
           endif
	   if (freedust.AND.microphys) then
             pdq_scav(ig,l,igcm_ccn_mass)=-pq(ig,l,igcm_ccn_mass)
     &                                     /ptimestep*(1-exp(
     &              -scav_ratio_wice*pcondicea(ig,l)*ptimestep/massl))
	     pdq_scav(ig,l,igcm_ccn_number)=pdq_scav(ig,l,igcm_ccn_mass)
     &                                    *r3n_q/rdust(ig,l)
           endif
           if (water) then
             pdq_scav(ig,l,igcm_h2o_ice)=-pq(ig,l,igcm_h2o_ice)
     &                                     /ptimestep*(1-exp(
     &              -scav_ratio_wice*pcondicea(ig,l)*ptimestep/massl))
           endif
     
         ELSE IF(pcondicea(ig,l).LT.0.)THEN ! if CO2 sublimates and releases dust/water ice
	   ! Calculation of the tendencies
	   if (freedust) then
	     pdq_scav(ig,l,igcm_dust_mass)=-pcondicea(ig,l)/massl*
     &                              scav_dust_mass(l+1)/pfallice(ig,l+1)
           
	     pdq_scav(ig,l,igcm_dust_number)=-pcondicea(ig,l)/massl*
     &                            scav_dust_number(l+1)/pfallice(ig,l+1)
           endif
           if (freedust.AND.microphys) then
	     pdq_scav(ig,l,igcm_ccn_mass)=-pcondicea(ig,l)/massl*
     &                              scav_ccn_mass(l+1)/pfallice(ig,l+1)
           
	     pdq_scav(ig,l,igcm_ccn_number)=-pcondicea(ig,l)/massl*
     &                            scav_ccn_number(l+1)/pfallice(ig,l+1)
           endif
	   if (water) then
	     pdq_scav(ig,l,igcm_h2o_ice)=-pcondicea(ig,l)/massl*
     &                              scav_h2o_ice(l+1)/pfallice(ig,l+1)
           endif
  
         END IF
	 ! Calculation of the scavenged dust/wice flux towards the lower layers
	 if (freedust) then
	   scav_dust_mass(l)=-pdq_scav(ig,l,igcm_dust_mass)*massl
     &                     +scav_dust_mass(l+1)
         
	   scav_dust_number(l)=-pdq_scav(ig,l,igcm_dust_number)*massl
     &                     +scav_dust_number(l+1)
         endif
	 if (freedust.AND.microphys) then
           scav_ccn_mass(l)=-pdq_scav(ig,l,igcm_ccn_mass)*massl
     &                     +scav_ccn_mass(l+1)
         
           scav_ccn_number(l)=-pdq_scav(ig,l,igcm_ccn_number)*massl
     &                     +scav_dust_number(l+1)
         endif
	 if (water) then
           scav_h2o_ice(l)=-pdq_scav(ig,l,igcm_h2o_ice)*massl
     &                     +scav_h2o_ice(l+1)
         endif
	 
       ENDDO
       ! Calculation of the surface tendencies
       pdqsc(ig,igcm_dust_mass)=0.
       pdqsc(ig,igcm_dust_number)=0.
       
       if (freedust) then
         pdqsc(ig,igcm_dust_mass)=pdqsc(ig,igcm_dust_mass)
     &                           +scav_dust_mass(1)
         pdqsc(ig,igcm_dust_number)=pdqsc(ig,igcm_dust_number)
     &                             +scav_dust_number(1)
       endif
       if (freedust.AND.microphys) then
         pdqsc(ig,igcm_dust_mass)=pdqsc(ig,igcm_dust_mass)
     &                           +scav_ccn_mass(1)
         pdqsc(ig,igcm_dust_number)=pdqsc(ig,igcm_dust_number)
     &                             +scav_ccn_number(1)
       endif
       if (water) then
         pdqsc(ig,igcm_h2o_ice)=scav_h2o_ice(1)
       endif
      ENDDO
      
      END SUBROUTINE scavenging_by_co2
      
      END MODULE co2condens_mod
