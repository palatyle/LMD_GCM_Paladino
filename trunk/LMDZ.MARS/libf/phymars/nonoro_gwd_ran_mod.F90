MODULE nonoro_gwd_ran_mod

IMPLICIT NONE

REAL,allocatable,save :: du_nonoro_gwd(:,:) ! Zonal wind tendency due to GWD 
REAL,allocatable,save :: dv_nonoro_gwd(:,:) ! Meridional wind tendency due to GWD
REAL,ALLOCATABLE,SAVE :: east_gwstress(:,:) ! Profile of eastward stress
REAL,ALLOCATABLE,SAVE :: west_gwstress(:,:) ! Profile of westward stress

CONTAINS

      SUBROUTINE NONORO_GWD_RAN(ngrid,nlayer,DTIME, pp,  &
                  zmax_therm, pt, pu, pv, pdt, pdu, pdv, &
                  zustr,zvstr,d_t, d_u, d_v)

    !--------------------------------------------------------------------------------
    ! Parametrization of the momentum flux deposition due to a discrete
    ! number of gravity waves. 
    ! F. Lott
    ! Version 14, Gaussian distribution of the source
    ! LMDz model online version      
    ! ADAPTED FOR VENUS /  F. LOTT + S. LEBONNOIS
    ! Version adapted on 03/04/2013:
    !      - input flux compensated in the deepest layers
    !                           
    ! ADAPTED FOR MARS     G.GILLI     02/2016
    !        Revision with F.Forget    06/2016  Variable EP-flux according to
    !                                           PBL variation (max velocity thermals)
    ! UPDATED              D.BARDET    01/2020  - reproductibility of the
    !                                           launching altitude calculation 
    !                                           - wave characteristic
    !                                           calculation using MOD
    !                                           - adding east_gwstress and
    !                                           west_gwstress variables    
    !---------------------------------------------------------------------------------

      use comcstfi_h, only: g, pi, cpp, r
      USE ioipsl_getin_p_mod, ONLY : getin_p
      use assert_m, only : assert
      use vertical_layers_mod, only : presnivs
      implicit none

      include "yoegwd.h" 
      include "callkeys.h"

      CHARACTER (LEN=20) :: modname='flott_gwd_rando'
      CHARACTER (LEN=80) :: abort_message


    ! 0. DECLARATIONS:

    ! 0.1 INPUTS
    INTEGER, intent(in):: ngrid, nlayer 
    REAL, intent(in):: DTIME ! Time step of the Physics
    REAL, intent(in):: zmax_therm(ngrid) ! altitude of max velocity thermals (m) 
    REAL, intent(in):: pp(ngrid,nlayer)   ! Pressure at full levels
    REAL, intent(in):: pt(ngrid,nlayer)   ! Temp at full levels 
    REAL, intent(in):: pu(ngrid,nlayer),pv(ngrid,nlayer) ! Hor winds at full levels
    REAL,INTENT(in) :: pdt(ngrid,nlayer) ! tendency on temperature
    REAL,INTENT(in) :: pdu(ngrid,nlayer) ! tendency on zonal wind
    REAL,INTENT(in) :: pdv(ngrid,nlayer) ! tendency on meridional wind

    ! 0.2 OUTPUTS
    REAL, intent(out):: zustr(ngrid), zvstr(ngrid) ! Surface Stresses
    REAL, intent(out):: d_t(ngrid, nlayer)        ! Tendency on Temp.
    REAL, intent(out):: d_u(ngrid, nlayer), d_v(ngrid, nlayer) ! tendency on winds

    ! O.3 INTERNAL ARRAYS
    REAL :: TT(ngrid, nlayer)   ! Temp at full levels 
    REAL :: UU(ngrid, nlayer) , VV(ngrid, nlayer) ! Hor winds at full levels
    REAL :: BVLOW(ngrid)
    REAL :: DZ

    INTEGER II, JJ, LL

    ! 0.3.0 TIME SCALE OF THE LIFE CYCLE OF THE WAVES PARAMETERIZED

    REAL, parameter:: DELTAT = 24. * 3600.
    
    ! 0.3.1 GRAVITY-WAVES SPECIFICATIONS

    INTEGER, PARAMETER:: NK = 2 ! number of horizontal wavenumbers
    INTEGER, PARAMETER:: NP = 2 ! directions (eastward and westward) phase speed
    INTEGER, PARAMETER:: NO = 2 ! absolute values of phase speed
    INTEGER, PARAMETER:: NW = NK * NP * NO ! Total numbers of gravity waves
    INTEGER JK, JP, JO, JW
    INTEGER, PARAMETER:: NA = 5 ! number of realizations to get the phase speed
    REAL, parameter:: kmax = 7.e-4 ! Max horizontal wavenumber
    REAL, parameter:: kmin = 2.e-5 ! Min horizontal wavenumber
    REAL, parameter:: cmax = 30.   ! Max horizontal absolute phase velocity
    REAL, parameter:: cmin = 1.    ! Min horizontal absolute phase velocity    
    REAL CPHA                   ! absolute PHASE VELOCITY frequency
    REAL ZK(NW, ngrid)           ! Horizontal wavenumber amplitude
    REAL ZP(NW, ngrid)           ! Horizontal wavenumber angle        
    REAL ZO(NW, ngrid)           ! Absolute frequency

    REAL intr_freq_m(nw, ngrid)          ! Waves Intr. freq. at the 1/2 lev below the full level (previous name: ZOM)
    REAL intr_freq_p(nw, ngrid)          ! Waves Intr. freq. at the 1/2 lev above the full level (previous name: ZOP)
    REAL wwm(nw, ngrid)                  ! Wave EP-fluxes at the 1/2 level below the full level
    REAL wwp(nw, ngrid)                  ! Wave EP-fluxes at the 1/2 level above the full level
    REAL u_epflux_p(nw, ngrid)           ! Partial zonal flux (=for each wave) at the 1/2 level above the full level (previous name: RUWP)
    REAL v_epflux_p(nw, ngrid)           ! Partial meridional flux (=for each wave) at the 1/2 level above the full level (previous name: RVWP)
    REAL u_epflux_tot(ngrid, nlayer + 1) ! Total zonal flux (=for all waves (nw)) at the 1/2 level above the full level (3D) (previous name: RUW)  
    REAL v_epflux_tot(ngrid, nlayer + 1) ! Total meridional flux (=for all waves (nw)) at the 1/2 level above the full level (3D) (previous name: RVW) 
    REAL epflux_0(nw, ngrid)             ! Fluxes at launching level (previous name: RUW0)
    REAL, save :: epflux_max             ! Max EP flux value at launching altitude (previous name: RUWMAX)
    INTEGER LAUNCH      ! Launching altitude
    REAL, parameter:: xlaunch = 0.4      ! Control the launching altitude
    REAL, parameter:: zmaxth_top = 8000. ! Top of convective layer (approx.)


    REAL PREC(ngrid)
    REAL PRMAX ! Maximum value of PREC, and for which our linear formula


    ! 0.3.2 PARAMETERS OF WAVES DISSIPATIONS
    REAL, parameter:: sat   = 1.     ! saturation parameter
    REAL, parameter:: rdiss = 1.     ! coefficient of dissipation
    REAL, parameter:: zoisec = 1.e-6 ! security for intrinsic freguency

    ! 0.3.3 Background flow at 1/2 levels and vertical coordinate
    REAL H0bis(ngrid, nlayer)          ! Characteristic Height of the atmosphere
    REAL, save:: H0                 ! Characteristic Height of the atmosphere
    REAL, parameter:: pr = 250      ! Reference pressure [Pa]
    REAL, parameter:: tr = 220.     ! Reference temperature [K]
    REAL ZH(ngrid, nlayer + 1)         ! Log-pressure altitude (constant H0)
    REAL ZHbis(ngrid, nlayer + 1)      ! Log-pressure altitude (varying H)
    REAL UH(ngrid, nlayer + 1)         ! zonal wind at 1/2 levels
    REAL VH(ngrid, nlayer + 1)         ! meridional wind at 1/2 levels
    REAL PH(ngrid, nlayer + 1)         ! Pressure at 1/2 levels
    REAL, parameter:: psec = 1.e-6  ! Security to avoid division by 0 pressure
    REAL BV(ngrid, nlayer + 1)         ! Brunt Vaisala freq. (BVF) at 1/2 levels
    REAL, parameter:: bvsec = 1.e-5 ! Security to avoid negative BV  
    REAL HREF(nlayer + 1)             ! Reference altitude for launching alt.


! COSMETICS TO DIAGNOSE EACH WAVES CONTRIBUTION.
    logical,save :: output=.false.
 ! CAUTION ! IF output is .true. THEN change NO to 10 at least !
    character*14 outform
    character*2  str2
    integer      ieq

    REAL RAN_NUM_1,RAN_NUM_2,RAN_NUM_3


    LOGICAL,SAVE :: firstcall = .true.


   !-----------------------------------------------------------------
    !  1. INITIALISATIONS

     IF (firstcall) THEN
        write(*,*) "nonoro_gwd_ran: FLott non-oro GW scheme is active!"
        epflux_max = 7.E-7 ! Mars' value !!
        call getin_p("nonoro_gwd_epflux_max", epflux_max)
        write(*,*) "nonoro_gwd_ran: epflux_max=", epflux_max
        ! Characteristic vertical scale height
        H0 = r * tr / g
        ! Control
        if (deltat .LT. dtime) THEN
             call abort_physic("nonoro_gwd_ran","gwd random: deltat lower than dtime!",1)
        endif
        if (nlayer .LT. nw) THEN
             call abort_physic("nonoro_gwd_ran","gwd random: nlayer lower than nw!",1)
        endif
        firstcall = .false.
     ENDIF

    gwd_convective_source=.false.

    ! Compute current values of temperature and winds
    tt(:,:)=pt(:,:)+dtime*pdt(:,:)
    uu(:,:)=pu(:,:)+dtime*pdu(:,:)
    vv(:,:)=pv(:,:)+dtime*pdv(:,:)



    !  2. EVALUATION OF THE BACKGROUND FLOW AT SEMI-LEVELS
    !-------------------------------------------------------------

    !Online output
    if (output) OPEN(11,file="impact-gwno.dat")

    ! Pressure and Inv of pressure, Temperature / at 1/2 level
    DO LL = 2, nlayer
       PH(:, LL) = EXP((LOG(PP(:, LL)) + LOG(PP(:, LL - 1))) / 2.)
    end DO

    PH(:, nlayer + 1) = 0. 
    PH(:, 1) = 2. * PP(:, 1) - PH(:, 2)

    ! Launching altitude

    !Pour revenir a la version non reproductible en changeant le nombre de
    !process
    ! Reprend la formule qui calcule PH en fonction de PP=play
    DO LL = 2, nlayer
       HREF(LL) = EXP((LOG(presnivs(LL))+ LOG(presnivs(LL - 1))) / 2.)
    end DO
    HREF(nlayer + 1) = 0.
    HREF(1) = 2. * presnivs(1) - HREF(2)

    LAUNCH=0
    DO LL = 1, nlayer
       IF (HREF(LL) / HREF(1) > XLAUNCH) LAUNCH = LL
    ENDDO
 
    if (output) print*, " WE ARE IN FLOTT GW SCHEME "
    
    ! Log pressure vert. coordinate
    DO LL = 1, nlayer + 1 
       ZH(:, LL) = H0 * LOG(PR / (PH(:, LL) + PSEC))
    end DO

    if (output) then
    ! altitude above surface
       ZHbis(:,1) = 0.    
       DO LL = 2, nlayer + 1 
          H0bis(:, LL-1) = r * TT(:, LL-1) / g 
          ZHbis(:, LL) = ZHbis(:, LL-1) &
           + H0bis(:, LL-1)*(PH(:, LL-1)-PH(:,LL))/PP(:, LL-1)
       end DO
    endif

    ! Winds and BV frequency
    DO LL = 2, nlayer
       UH(:, LL) = 0.5 * (UU(:, LL) + UU(:, LL - 1)) ! Zonal wind
       VH(:, LL) = 0.5 * (VV(:, LL) + VV(:, LL - 1)) ! Meridional wind
       ! GG test	
       !print*, 'TT, UH, VH, ZH at launch', TT(ngrid/2,LAUNCH), UH(ngrid/2,LAUNCH),VH(ngrid/2, LAUNCH), ZH(ngrid/2,LAUNCH)	
       ! BVSEC: BV Frequency
       BV(:, LL) = 0.5 * (TT(:, LL) + TT(:, LL - 1)) &
            * r**2 / cpp / H0**2 + (TT(:, LL) &
            - TT(:, LL - 1)) / (ZH(:, LL) - ZH(:, LL - 1)) * r / H0
       BV(:,LL) =SQRT(MAX(BVSEC,BV(:,LL)))
    end DO
       !GG test
       !print*, 'BV freq in flott_gwnoro:',LAUNCH,  BV(ngrid/2, LAUNCH)  

    BV(:, 1) = BV(:, 2)
    UH(:, 1) = 0.
    VH(:, 1) = 0.
    BV(:, nlayer + 1) = BV(:, nlayer)
    UH(:, nlayer + 1) = UU(:, nlayer)
    VH(:, nlayer + 1) = VV(:, nlayer)


    ! 3. WAVES CHARACTERISTICS CHOSEN RANDOMLY
    !-------------------------------------------

    ! The mod function of here a weird arguments
    ! are used to produce the waves characteristics
    ! in a stochastic way

    DO JW = 1, NW
             !  Angle
             DO II = 1, ngrid
                ! Angle (0 or PI so far)
                RAN_NUM_1=MOD(TT(II, JW) * 10., 1.)
                RAN_NUM_2= MOD(TT(II, JW) * 100., 1.)
                ZP(JW, II) = (SIGN(1., 0.5 - RAN_NUM_1) + 1.) &
                     * PI / 2.
                ! Horizontal wavenumber amplitude
                ZK(JW, II) = KMIN + (KMAX - KMIN) *RAN_NUM_2
                ! Horizontal phase speed
                CPHA = 0.
                DO JJ = 1, NA
                    RAN_NUM_3=MOD(TT(II, JW+3*JJ)**2, 1.)
                    CPHA = CPHA + &
                    CMAX*2.*(RAN_NUM_3 -0.5)*SQRT(3.)/SQRT(NA*1.)
                END DO
                IF (CPHA.LT.0.)  THEN
                   CPHA = -1.*CPHA
                   ZP(JW,II) = ZP(JW,II) + PI
                ENDIF
       !Online output: linear
                if (output) CPHA = CMIN + (CMAX - CMIN) * (JO-1)/(NO-1)
                ! Intrinsic frequency
                ZO(JW, II) = CPHA * ZK(JW, II)
                ! Intrinsic frequency  is imposed
                    ZO(JW, II) = ZO(JW, II)      &
                  + ZK(JW, II) * COS(ZP(JW, II)) * UH(II, LAUNCH) &
                  + ZK(JW, II) * SIN(ZP(JW, II)) * VH(II, LAUNCH)
                ! Momentum flux at launch lev 
                ! epflux_0(JW, II) = epflux_max / REAL(NW) &
                epflux_0(JW, II) = epflux_max &
                     * MOD(100. * (UU(II, JW)**2 + VV(II, JW)**2), 1.)
       !Online output: fixed to max
                if (output) epflux_0(JW, II) = epflux_max
             ENDDO
   end DO

    ! 4. COMPUTE THE FLUXES
    !--------------------------

    !  4.1  Vertical velocity at launching altitude to ensure 
    !       the correct value to the imposed fluxes.
    !
    DO JW = 1, NW
       ! Evaluate intrinsic frequency at launching altitude:
       intr_freq_p(JW, :) = ZO(JW, :) &
            - ZK(JW, :) * COS(ZP(JW, :)) * UH(:, LAUNCH) &
            - ZK(JW, :) * SIN(ZP(JW, :)) * VH(:, LAUNCH) 
    end DO

    IF (gwd_convective_source) THEN
         DO JW = 1, NW
       ! VERSION WITH CONVECTIVE SOURCE (designed for Earth)

       ! Vertical velocity at launch level, value to ensure the
       ! imposed mmt flux factor related to the convective forcing:
       ! precipitations.

       ! tanh limitation to values above prmax:
!       WWP(JW, :) = epflux_0(JW, :) &
!            * (r / cpp / H0 * RLVTT * PRMAX * TANH(PREC(:) / PRMAX))**2
!       Here, we neglected the kinetic energy providing of the thermodynamic
!       phase change

!

       ! Factor related to the characteristics of the waves:
            WWP(JW, :) = WWP(JW, :) * ZK(JW, :)**3 / KMIN / BVLOW(:)  &
                 / MAX(ABS(intr_freq_p(JW, :)), ZOISEC)**3

      ! Moderation by the depth of the source (dz here):
            WWP(JW, :) = WWP(JW, :) &
                 * EXP(- BVLOW(:)**2 / MAX(ABS(intr_freq_p(JW, :)), ZOISEC)**2 &
                 * ZK(JW, :)**2 * DZ**2)

      ! Put the stress in the right direction:
            u_epflux_p(JW, :) = intr_freq_p(JW, :) / MAX(ABS(intr_freq_p(JW, :)), ZOISEC)**2 &
                 * BV(:, LAUNCH) * COS(ZP(JW, :)) * WWP(JW, :)**2
            v_epflux_p(JW, :) = intr_freq_p(JW, :) / MAX(ABS(intr_freq_p(JW, :)), ZOISEC)**2 &
                 * BV(:, LAUNCH) * SIN(ZP(JW, :)) * WWP(JW, :)**2
         end DO
    ELSE ! VERSION WITHOUT CONVECTIVE SOURCE
       ! Vertical velocity at launch level, value to ensure the imposed
       ! mom flux:
         DO JW = 1, NW
       ! WW is directly a flux, here, not vertical velocity anymore
            WWP(JW, :) = epflux_0(JW,:)
            u_epflux_p(JW, :) = COS(ZP(JW, :)) * SIGN(1., intr_freq_p(JW, :)) * epflux_0(JW, :)
            v_epflux_p(JW, :) = SIN(ZP(JW, :)) * SIGN(1., intr_freq_p(JW, :)) * epflux_0(JW, :)

         end DO
    ENDIF
    !  4.2 Initial flux at launching altitude

    u_epflux_tot(:, LAUNCH) = 0
    v_epflux_tot(:, LAUNCH) = 0
    DO JW = 1, NW
       u_epflux_tot(:, LAUNCH) = u_epflux_tot(:, LAUNCH) + u_epflux_p(JW, :)
       v_epflux_tot(:, LAUNCH) = v_epflux_tot(:, LAUNCH) + v_epflux_p(JW, :)
    end DO

    !  4.3 Loop over altitudes, with passage from one level to the
    !      next done by i) conserving the EP flux, ii) dissipating
    !      a little, iii) testing critical levels, and vi) testing
    !      the breaking.

    !Online output
    if (output) then
        ieq=ngrid/2+1
        write(str2,'(i2)') NW+2
        outform="("//str2//"(E12.4,1X))"
        WRITE(11,outform) ZH(IEQ, 1) / 1000., ZHbis(IEQ, 1) / 1000., &
               (ZO(JW, IEQ)/ZK(JW, IEQ)*COS(ZP(JW, IEQ)), JW = 1, NW)
    endif

    DO LL = LAUNCH, nlayer - 1


       !  W(KB)ARNING: ALL THE PHYSICS IS HERE (PASSAGE FROM ONE LEVEL
       ! TO THE NEXT)
       DO JW = 1, NW
          intr_freq_m(JW, :) = intr_freq_p(JW, :)
          WWM(JW, :) = WWP(JW, :)
          ! Intrinsic Frequency
          intr_freq_p(JW, :) = ZO(JW, :) - ZK(JW, :) * COS(ZP(JW,:)) * UH(:, LL + 1) &
               - ZK(JW, :) * SIN(ZP(JW,:)) * VH(:, LL + 1) 

          WWP(JW, :) = MIN( & 
       ! No breaking (Eq.6)
               WWM(JW, :) & 
      ! Dissipation (Eq. 8):
               * EXP(- RDISS * PR / (PH(:, LL + 1) + PH(:, LL)) &
               * ((BV(:, LL + 1) + BV(:, LL)) / 2.)**3 &
               / MAX(ABS(intr_freq_p(JW, :) + intr_freq_m(JW, :)) / 2., ZOISEC)**4 &
               * ZK(JW, :)**3 * (ZH(:, LL + 1) - ZH(:, LL))), &
       ! Critical levels (forced to zero if intrinsic frequency changes sign)
               MAX(0., SIGN(1., intr_freq_p(JW, :) * intr_freq_m(JW, :))) &
       ! Saturation (Eq. 12)
               * ABS(intr_freq_p(JW, :))**3 /BV(:, LL+1) & 
               * EXP(-ZH(:, LL + 1)/H0) * SAT**2*KMIN**2/ZK(JW, :)**4)  
       end DO

       ! END OF W(KB)ARNING
       ! Evaluate EP-flux from Eq. 7 and 
       ! Give the right orientation to the stress
       DO JW = 1, NW
          u_epflux_p(JW, :) = SIGN(1.,intr_freq_p(JW, :)) * COS(ZP(JW, :)) *  WWP(JW, :)
          v_epflux_p(JW, :) = SIGN(1.,intr_freq_p(JW, :)) * SIN(ZP(JW, :)) *  WWP(JW, :)
       end DO
 
       u_epflux_tot(:, LL + 1) = 0.
       v_epflux_tot(:, LL + 1) = 0.

       DO JW = 1, NW
          u_epflux_tot(:, LL + 1) = u_epflux_tot(:, LL + 1) + u_epflux_p(JW, :) 
          v_epflux_tot(:, LL + 1) = v_epflux_tot(:, LL + 1) + v_epflux_p(JW, :) 
          EAST_GWSTRESS(:, LL)=EAST_GWSTRESS(:, LL)+MAX(0.,u_epflux_p(JW,:))/FLOAT(NW)
          WEST_GWSTRESS(:, LL)=WEST_GWSTRESS(:, LL)+MIN(0.,u_epflux_p(JW,:))/FLOAT(NW)
       end DO
       !Online output
       if (output) then
         do JW=1,NW
            if(u_epflux_p(JW, IEQ).gt.0.) then
              u_epflux_p(JW, IEQ) = max(u_epflux_p(JW, IEQ), 1.e-99)
            else
              u_epflux_p(JW, IEQ) = min(u_epflux_p(JW, IEQ), -1.e-99)
            endif
         enddo
                   WRITE(11,outform) ZH(IEQ, LL+1) / 1000., &
                                  ZHbis(IEQ, LL+1) / 1000., &
                                  (u_epflux_p(JW, IEQ), JW = 1, NW)
       endif

    end DO

    ! 5 CALCUL DES TENDANCES:
    !------------------------

    ! 5.1 Rectification des flux au sommet et dans les basses couches:

! Attention, ici c'est le total sur toutes les ondes...

    u_epflux_tot(:, nlayer + 1) = 0.
    v_epflux_tot(:, nlayer + 1) = 0.

    ! Here, big change compared to FLott version:
    ! We compensate (u_epflux_tot(:, LAUNCH), ie total emitted upward flux
    !  over the layers max(1,LAUNCH-3) to LAUNCH-1
    DO LL = 1, max(1,LAUNCH-3)
      u_epflux_tot(:, LL) = 0.
      v_epflux_tot(:, LL) = 0.
    end DO
    DO LL = max(2,LAUNCH-2), LAUNCH-1
       u_epflux_tot(:, LL) = u_epflux_tot(:, LL - 1) + u_epflux_tot(:, LAUNCH) * &
            (PH(:,LL)-PH(:,LL-1)) / (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
       v_epflux_tot(:, LL) = v_epflux_tot(:, LL - 1) + v_epflux_tot(:, LAUNCH) * &
            (PH(:,LL)-PH(:,LL-1)) / (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
       EAST_GWSTRESS(:,LL) = EAST_GWSTRESS(:, LL - 1) + &
            EAST_GWSTRESS(:, LAUNCH) * (PH(:,LL)-PH(:,LL-1))/ &
            (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
       WEST_GWSTRESS(:,LL) = WEST_GWSTRESS(:, LL - 1) + &
            WEST_GWSTRESS(:, LAUNCH) * (PH(:,LL)-PH(:,LL-1))/ &
            (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
    end DO
    ! This way, the total flux from GW is zero, but there is a net transport
    ! (upward) that should be compensated by circulation 
    ! and induce additional friction at the surface 

    !Online output
    if (output) then
       DO LL = 1, nlayer - 1
           WRITE(11,*) ZHbis(IEQ, LL)/1000.,u_epflux_tot(IEQ,LL)
       end DO
       CLOSE(11)
       stop
    endif


    ! 5.2 AR-1 RECURSIVE FORMULA (13) IN VERSION 4
    !---------------------------------------------
    DO LL = 1, nlayer
       d_u(:, LL) = G * (u_epflux_tot(:, LL + 1) - u_epflux_tot(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
       d_v(:, LL) = G * (v_epflux_tot(:, LL + 1) - v_epflux_tot(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
    ENDDO
    d_t(:,:) = 0.

    ! 5.3 Update tendency of wind with the previous (and saved) values
    !-----------------------------------------------------------------
    d_u(:,:) = DTIME/DELTAT/REAL(NW) * d_u(:,:)                       &
             + (1.-DTIME/DELTAT) * du_nonoro_gwd(:,:)
    d_v(:,:) = DTIME/DELTAT/REAL(NW) * d_v(:,:)                       &
             + (1.-DTIME/DELTAT) * dv_nonoro_gwd(:,:)
    du_nonoro_gwd(:,:) = d_u(:,:)
    dv_nonoro_gwd(:,:) = d_v(:,:)

    ! Cosmetic: evaluation of the cumulated stress

    ZUSTR(:) = 0.
    ZVSTR(:) = 0.
    DO LL = 1, nlayer
       ZUSTR(:) = ZUSTR(:) + D_U(:, LL) / g * (PH(:, LL + 1) - PH(:, LL))
       ZVSTR(:) = ZVSTR(:) + D_V(:, LL) / g * (PH(:, LL + 1) - PH(:, LL))
    ENDDO

  END SUBROUTINE NONORO_GWD_RAN

! ========================================================
! Subroutines used to allocate/deallocate module variables       
! ========================================================

  SUBROUTINE ini_nonoro_gwd_ran(ngrid,nlayer)

  IMPLICIT NONE

      INTEGER, INTENT (in) :: ngrid  ! number of atmospheric columns
      INTEGER, INTENT (in) :: nlayer ! number of atmospheric layers

         allocate(du_nonoro_gwd(ngrid,nlayer))
         allocate(dv_nonoro_gwd(ngrid,nlayer))
         allocate(east_gwstress(ngrid,nlayer))
         east_gwstress(:,:)=0
         allocate(west_gwstress(ngrid,nlayer))
         west_gwstress(:,:)=0

  END SUBROUTINE ini_nonoro_gwd_ran
! ----------------------------------
  SUBROUTINE end_nonoro_gwd_ran

  IMPLICIT NONE

         if (allocated(du_nonoro_gwd)) deallocate(du_nonoro_gwd)
         if (allocated(dv_nonoro_gwd)) deallocate(dv_nonoro_gwd)
         if (allocated(east_gwstress)) deallocate(east_gwstress)
         if (allocated(west_gwstress)) deallocate(west_gwstress)

  END SUBROUTINE end_nonoro_gwd_ran

END MODULE nonoro_gwd_ran_mod
