      SUBROUTINE FLOTT_GWD_RAN(NLON,NLEV,DTIME, pp, pn2,  &
                  tt,uu,vv, plevmoy, &
                  zustr,zvstr,d_t, d_u, d_v)

    !----------------------------------------------------------------------
    ! Parametrization of the momentum flux deposition due to a discrete
    ! number of gravity waves. 
    ! F. Lott
    ! Version 14, Gaussian distribution of the source
    ! LMDz model online version      
    ! ADAPTED FOR VENUS /  F. LOTT + S. LEBONNOIS
    ! Version adapted on 03/04/2013:
    !      - input flux compensated in the deepest layers
    !---------------------------------------------------------------------

      use dimphy
      implicit none

#include "YOEGWD.h"
#include "YOMCST.h"

    ! 0. DECLARATIONS:

    ! 0.1 INPUTS
    INTEGER, intent(in):: NLON, NLEV 
    REAL, intent(in):: DTIME ! Time step of the Physics
    REAL, intent(in):: pp(NLON, NLEV)   ! Pressure at full levels
! VENUS ATTENTION: CP VARIABLE PN2 CALCULE EN AMONT DES PARAMETRISATIONS
    REAL, intent(in):: pn2(NLON,NLEV)   ! N2 (BV^2) at 1/2 levels
    REAL, intent(in):: TT(NLON, NLEV)   ! Temp at full levels 

    REAL, intent(in):: UU(NLON, NLEV) , VV(NLON, NLEV)
    ! Hor winds at full levels
    REAL, intent(in) :: plevmoy(NLEV+1) ! press (Pa) at interlayers, at klon/2+1

    ! 0.2 OUTPUTS
    REAL, intent(out):: zustr(NLON), zvstr(NLON) ! Surface Stresses
    REAL, intent(inout):: d_t(NLON, NLEV)        ! Tendency on Temp.

    REAL, intent(inout):: d_u(NLON, NLEV), d_v(NLON, NLEV)
    ! Tendencies on winds

    ! O.3 INTERNAL ARRAYS

    INTEGER II, LL

    ! 0.3.0 TIME SCALE OF THE LIFE CYCLE OF THE WAVES PARAMETERIZED

    REAL DELTAT

    ! 0.3.1 GRAVITY-WAVES SPECIFICATIONS

!VENUS 
    INTEGER, PARAMETER:: NK = 2, NP = 2, NO = 2, NW = NK * NP * NO
    !Online output: change NO
!    INTEGER, PARAMETER:: NK = 1, NP = 2, NO = 11, NW = NK * NP * NO
    INTEGER JK, JP, JO, JW
    REAL KMIN, KMAX ! Min and Max horizontal wavenumbers
    REAL CMIN, CMAX ! Min and Max absolute ph. vel.
    REAL CPHA ! absolute PHASE VELOCITY frequency
    REAL ZK(NW, KLON) ! Horizontal wavenumber amplitude
    REAL ZP(NW)       ! Horizontal wavenumber angle        
    REAL ZO(NW, KLON) ! Absolute frequency      !

    ! Waves Intr. freq. at the 1/2 lev surrounding the full level
    REAL ZOM(NW, KLON), ZOP(NW, KLON)

    ! Wave EP-fluxes at the 2 semi levels surrounding the full level
    REAL WWM(NW, KLON), WWP(NW, KLON)
    ! Fluxes X and Y for each waves at 1/2 Levels
    REAL RUWP(NW, KLON), RVWP(NW, KLON)
    REAL RUW(KLON, KLEV + 1) ! Flux x at semi levels
    REAL RVW(KLON, KLEV + 1) ! Flux y at semi levels

    REAL RUW0(NW, KLON) ! Fluxes at launching level
    REAL RUWMAX         ! Max value
    INTEGER LAUNCH      ! Launching altitude
    REAL XLAUNCH        ! Control the launching altitude

    ! 0.3.2 PARAMETERS OF WAVES DISSIPATIONS

    REAL SAT  ! saturation parameter
    REAL RDISS, ZOISEC ! COEFF DE DISSIPATION, SECURITY FOR INTRINSIC FREQ

    ! 0.3.3 BACKGROUND FLOW AT 1/2 LEVELS AND VERTICAL COORDINATE

    REAL H0bis(KLON, KLEV) ! Characteristic Height of the atmosphere
    REAL H0 ! Characteristic Height of the atmosphere
    REAL PR, TR ! Reference Pressure and Temperature

    REAL ZH(KLON, KLEV + 1) ! Log-pressure altitude (constant H0)
    REAL ZHbis(KLON, KLEV + 1) ! Log-pressure altitude (varying H)

    REAL UH(KLON, KLEV + 1), VH(KLON, KLEV + 1) ! Winds at 1/2 levels
    REAL PH(KLON, KLEV + 1) ! Pressure at 1/2 levels
    REAL PSEC ! Security to avoid division by 0 pressure
    REAL BV(KLON, KLEV + 1) ! Brunt Vaisala freq. (BVF) at 1/2 levels
    REAL BVSEC ! Security to avoid negative BVF

! COSMETICS TO DIAGNOSE EACH WAVES CONTRIBUTION.
    logical output
    data output/.false./
!    data output/.true./
 ! CAUTION ! IF output is .true. THEN change NO to 10 at least !
    character*14 outform
    character*2  str2
    integer      ieq

! ON CONSERVE LA MEMOIRE un certain temps AVEC UN SAVE
    real,save,allocatable :: d_u_sav(:,:),d_v_sav(:,:)
    LOGICAL firstcall
    SAVE firstcall
    DATA firstcall/.true./

      REAL ALEAS
      EXTERNAL ALEAS

    !-----------------------------------------------------------------
    !  1. INITIALISATIONS

      IF (firstcall) THEN
        allocate(d_u_sav(NLON,NLEV),d_v_sav(NLON,NLEV))
        d_u_sav = 0.
        d_v_sav = 0.
        firstcall=.false.
      ENDIF
    
    !    1.1 Basic parameter

    !  PARAMETERS CORRESPONDING TO V3:
    RUWMAX = 0.005      ! Max EP-Flux at Launch altitude
    SAT    = 0.85       ! Saturation parameter: Sc in (12)
    RDISS  = 0.1        ! Diffusion parameter 

    DELTAT=24.*3600.    ! Time scale of the waves (first introduced in 9b)

!!!! TEST GG Values corresponding to min/max horizontal wavel 50-500 km
!(similar to observations)
    KMIN = 1.E-5        ! Min horizontal wavenumber
!    KMIN = 6.3E-6       ! Min horizontal wavenumber
    KMAX = 1.E-4        ! Max horizontal wavenumber
    !Online output: one value only
    if (output) then
      KMIN = 3E-5
      KMAX = 3E-5
    endif
    CMIN = 1.           ! Min phase velocity
    CMAX = 61.          ! Max phase speed velocity
!    XLAUNCH=0.6         ! Parameter that control launching altitude
    XLAUNCH=5e-3        ! Value for top of cloud convective region

!    PR = 9.2e6          ! Reference pressure    ! VENUS!!
    PR = 5e5            ! Reference pressure    ! VENUS: cloud layer
    TR = 300.           ! Reference Temperature ! VENUS: cloud layer
    H0 = RD * TR / RG   ! Characteristic vertical scale height

    BVSEC  = 1.E-5      ! Security to avoid negative BVF  
    PSEC   = 1.E-8      ! Security to avoid division by 0 pressure
    ZOISEC = 1.E-8      ! Security FOR 0 INTRINSIC FREQ

    IF(DELTAT.LT.DTIME)THEN
       PRINT *,'GWD RANDO: DELTAT LT DTIME!'
       STOP
    ENDIF


    IF (NLEV < NW) THEN
       PRINT *, 'YOU WILL HAVE PROBLEM WITH RANDOM NUMBERS'
       PRINT *, 'FLOTT GWD STOP'
       STOP 1
    ENDIF

    !  2. EVALUATION OF THE BACKGROUND FLOW AT SEMI-LEVELS
    !-------------------------------------------------------------

    !Online output
    if (output) OPEN(11,file="impact-gwno.dat")

    ! Pressure and Inv of pressure, Temperature / at 1/2 level
    DO LL = 2, KLEV
       PH(:, LL) = EXP((LOG(PP(:, LL)) + LOG(PP(:, LL - 1))) / 2.)
    end DO

    PH(:, KLEV + 1) = 0. 
    PH(:, 1) = 2. * PP(:, 1) - PH(:, 2)

    ! Launching altitude

    DO LL = 1, NLEV
       IF (plevmoy(LL) / plevmoy(1) > XLAUNCH) LAUNCH = LL
    ENDDO
! test
!    print*,"launch=",LAUNCH
!    print*,"launch p,N2=",plevmoy(LAUNCH),pn2(nlon/2+1,LAUNCH)

    ! Log pressure vert. coordinate
    DO LL = 1, KLEV + 1 
       ZH(:, LL) = H0 * LOG(PR / (PH(:, LL) + PSEC))
    end DO

    if (output) then
    ! altitude above surface
       ZHbis(:,1) = 0.    
       DO LL = 2, KLEV + 1 
          H0bis(:, LL-1) = RD * TT(:, LL-1) / RG 
          ZHbis(:, LL) = ZHbis(:, LL-1) &
           + H0bis(:, LL-1)*(PH(:, LL-1)-PH(:,LL))/PP(:, LL-1)
       end DO
    endif

    ! Winds and BV frequency
    DO LL = 2, KLEV
       UH(:, LL) = 0.5 * (UU(:, LL) + UU(:, LL - 1)) ! Zonal wind
       VH(:, LL) = 0.5 * (VV(:, LL) + VV(:, LL - 1)) ! Meridional wind
       ! BVSEC: BV Frequency
! VENUS ATTENTION: CP VARIABLE PSTAB CALCULE EN AMONT DES PARAMETRISATIONS
       BV(:, LL) = MAX(BVSEC,SQRT(pn2(:,LL)))
    end DO
    BV(:, 1) = BV(:, 2)
    UH(:, 1) = 0.
    VH(:, 1) = 0.
    BV(:, KLEV + 1) = BV(:, KLEV)
    UH(:, KLEV + 1) = UU(:, KLEV)
    VH(:, KLEV + 1) = VV(:, KLEV)


    ! 3. WAVES CHARACTERISTICS CHOSEN RANDOMLY
    !-------------------------------------------

    ! The mod function of here a weird arguments
    ! are used to produce the waves characteristics
    ! in a stochastic way

!! A REVOIR: 
!! - utilisation de MOD ou bien de aleas ?
!! - distribution gaussienne des CPHA ? (avec signe ZP qui est ajuste apres)

    JW = 0
    DO JP = 1, NP
       DO JK = 1, NK
          DO JO = 1, NO
             JW = JW + 1
             !  Angle
             ZP(JW) = 2. * RPI * REAL(JP - 1) / REAL(NP) 
             DO II = 1, KLON
                ! Horizontal wavenumber amplitude
!                ZK(JW, II) = KMIN + (KMAX - KMIN) * MOD(TT(II, JW) * 100., 1.)
                ZK(JW, II) = KMIN + (KMAX - KMIN) * ALEAS(0.)
                ! Horizontal phase speed
!                CPHA = CMIN + (CMAX - CMIN) * MOD(TT(II, JW)**2, 1.)
                CPHA = CMIN + (CMAX - CMIN) * ALEAS(0.)
       !Online output: linear
                if (output) CPHA = CMIN + (CMAX - CMIN) * (JO-1)/(NO-1)
                ! Intrinsic frequency
                ZO(JW, II) = CPHA * ZK(JW, II)
                ! Intrinsic frequency  is imposed
                    ZO(JW, II) = ZO(JW, II)      &
                  + ZK(JW, II) * COS(ZP(JW)) * UH(II, LAUNCH) &
                  + ZK(JW, II) * SIN(ZP(JW)) * VH(II, LAUNCH)
                ! Momentum flux at launch lev 
                ! RUW0(JW, II) = RUWMAX / REAL(NW) &
                RUW0(JW, II) = RUWMAX &
!                     * MOD(100. * (UU(II, JW)**2 + VV(II, JW)**2), 1.)
                     * ALEAS(0.)
       !Online output: fixed to max
                if (output) RUW0(JW, II) = RUWMAX
             ENDDO
          end DO
       end DO
    end DO

    ! 4. COMPUTE THE FLUXES
    !--------------------------

    !  4.1  Vertical velocity at launching altitude to ensure 
    !       the correct value to the imposed fluxes.
    !
    DO JW = 1, NW

       ! Evaluate intrinsic frequency at launching altitude:
       ZOP(JW, :) = ZO(JW, :) &
            - ZK(JW, :) * COS(ZP(JW)) * UH(:, LAUNCH) &
            - ZK(JW, :) * SIN(ZP(JW)) * VH(:, LAUNCH) 
       ! WW is directly a flux, here, not vertical velocity anymore
       WWP(JW, :) = RUW0(JW,:)
       RUWP(JW, :) = COS(ZP(JW)) * SIGN(1., ZOP(JW, :)) * RUW0(JW, :)
       RVWP(JW, :) = SIN(ZP(JW)) * SIGN(1., ZOP(JW, :)) * RUW0(JW, :)

    end DO

    !  4.2 Initial flux at launching altitude

    RUW(:, LAUNCH) = 0
    RVW(:, LAUNCH) = 0
    DO JW = 1, NW
       RUW(:, LAUNCH) = RUW(:, LAUNCH) + RUWP(JW, :)
       RVW(:, LAUNCH) = RVW(:, LAUNCH) + RVWP(JW, :)
    end DO

    !  4.3 Loop over altitudes, with passage from one level to the
    !      next done by i) conserving the EP flux, ii) dissipating
    !      a little, iii) testing critical levels, and vi) testing
    !      the breaking.

    !Online output
    if (output) then
        ieq=nlon/2+1
        write(str2,'(i2)') NW+2
        outform="("//str2//"(E12.4,1X))"
        WRITE(11,outform) ZH(IEQ, 1) / 1000., ZHbis(IEQ, 1) / 1000., &
               (ZO(JW, IEQ)/ZK(JW, IEQ)*COS(ZP(JW)), JW = 1, NW)
    endif

    DO LL = LAUNCH, KLEV - 1


       !  W(KB)ARNING: ALL THE PHYSICS IS HERE (PASSAGE FROM ONE LEVEL
       ! TO THE NEXT)
       DO JW = 1, NW
          ZOM(JW, :) = ZOP(JW, :)
          WWM(JW, :) = WWP(JW, :)
          ! Intrinsic Frequency
          ZOP(JW, :) = ZO(JW, :) - ZK(JW, :) * COS(ZP(JW)) * UH(:, LL + 1) &
               - ZK(JW, :) * SIN(ZP(JW)) * VH(:, LL + 1) 

          WWP(JW, :) = MIN( & 
       ! No breaking (Eq.6)
               WWM(JW, :) & 
       ! Dissipation (Eq. 8):
               * EXP(- RDISS * PR / (PH(:, LL + 1) + PH(:, LL)) &
               * ((BV(:, LL + 1) + BV(:, LL)) / 2.)**3 &
               / MAX(ABS(ZOP(JW, :) + ZOM(JW, :)) / 2., ZOISEC)**4 &
               * ZK(JW, :)**3 * (ZH(:, LL + 1) - ZH(:, LL))), &
       ! Critical levels (forced to zero if intrinsic frequency changes sign)
               MAX(0., SIGN(1., ZOP(JW, :) * ZOM(JW, :))) &
       ! Saturation (Eq. 12)
               * ABS(ZOP(JW, :))**3 /BV(:, LL+1) & 
               * EXP(-ZH(:, LL + 1)/H0) * SAT**2*KMIN**2/ZK(JW, :)**4)  
       end DO

       ! END OF W(KB)ARNING
       ! Evaluate EP-flux from Eq. 7 and 
       ! Give the right orientation to the stress

       DO JW = 1, NW
          RUWP(JW, :) = SIGN(1.,ZOP(JW, :))*COS(ZP(JW))*WWP(JW, :)
          RVWP(JW, :) = SIGN(1.,ZOP(JW, :))*SIN(ZP(JW))*WWP(JW, :)
       end DO
       !
       RUW(:, LL + 1) = 0.
       RVW(:, LL + 1) = 0.

       DO JW = 1, NW
          RUW(:, LL + 1) = RUW(:, LL + 1) + RUWP(JW, :) 
          RVW(:, LL + 1) = RVW(:, LL + 1) + RVWP(JW, :) 
       end DO
       !Online output
       if (output) then
         do JW=1,NW
            if(RUWP(JW, IEQ).gt.0.) then
              RUWP(JW, IEQ) = max(RUWP(JW, IEQ), 1.e-99)
            else
              RUWP(JW, IEQ) = min(RUWP(JW, IEQ), -1.e-99)
            endif
         enddo
                   WRITE(11,outform) ZH(IEQ, LL+1) / 1000., &
                                  ZHbis(IEQ, LL+1) / 1000., &
                                  (RUWP(JW, IEQ), JW = 1, NW)
       endif

    end DO

    ! 5 CALCUL DES TENDANCES:
    !------------------------

    ! 5.1 Rectification des flux au sommet et dans les basses couches:
! MODIF SL
 
! Attention, ici c'est le total sur toutes les ondes...

    RUW(:, KLEV + 1) = 0.
    RVW(:, KLEV + 1) = 0.

    ! Here, big change compared to FLott version:
    ! We compensate (RUW(:, LAUNCH), ie total emitted upward flux
    !  over the layers max(1,LAUNCH-3) to LAUNCH-1
    DO LL = 1, max(1,LAUNCH-3)
      RUW(:, LL) = 0.
      RVW(:, LL) = 0.
    end DO
    DO LL = max(2,LAUNCH-2), LAUNCH-1
       RUW(:, LL) = RUW(:, LL - 1) + RUW(:, LAUNCH) * &
            (PH(:,LL)-PH(:,LL-1)) / (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
       RVW(:, LL) = RVW(:, LL - 1) + RVW(:, LAUNCH) * &
            (PH(:,LL)-PH(:,LL-1)) / (PH(:,LAUNCH)-PH(:,max(1,LAUNCH-3)))
    end DO
    ! This way, the total flux from GW is zero, but there is a net transport
    ! (upward) that should be compensated by circulation 
    ! and induce additional friction at the surface 

    !Online output
    if (output) then 
       DO LL = 1, KLEV - 1
           WRITE(11,*) ZHbis(IEQ, LL)/1000.,RUW(IEQ,LL)
       end DO
       CLOSE(11)
       stop
    endif

    ! AR-1 RECURSIVE FORMULA (13) IN VERSION 4
    DO LL = 1, KLEV
       d_u(:, LL) = RG * (RUW(:, LL + 1) - RUW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
       d_v(:, LL) = RG * (RVW(:, LL + 1) - RVW(:, LL)) &
            / (PH(:, LL + 1) - PH(:, LL)) * DTIME
    ENDDO
        d_t = 0.
    ! ON CONSERVE LA MEMOIRE un certain temps AVEC UN SAVE
        d_u = DTIME/DELTAT/REAL(NW) * d_u + (1.-DTIME/DELTAT) * d_u_sav
        d_v = DTIME/DELTAT/REAL(NW) * d_v + (1.-DTIME/DELTAT) * d_v_sav
	d_u_sav = d_u
	d_v_sav = d_v

    ! Cosmetic: evaluation of the cumulated stress

    ZUSTR(:) = 0.
    ZVSTR(:) = 0.
    DO LL = 1, KLEV
       ZUSTR(:) = ZUSTR(:) + D_U(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))
       ZVSTR(:) = ZVSTR(:) + D_V(:, LL) / RG * (PH(:, LL + 1) - PH(:, LL))
    ENDDO

  END SUBROUTINE FLOTT_GWD_RAN

!===================================================================
!===================================================================
!===================================================================
!===================================================================

      FUNCTION ALEAS (R)
!***BEGIN PROLOGUE  ALEAS
!***PURPOSE  Generate a uniformly distributed random number.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  L6A21
!***TYPE      SINGLE PRECISION (ALEAS-S)
!***KEYWORDS  FNLIB, ALEAS NUMBER, SPECIAL FUNCTIONS, UNIFORM
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!      This pseudo-random number generator is portable among a wide
! variety of computers.  RAND(R) undoubtedly is not as good as many
! readily available installation dependent versions, and so this
! routine is not recommended for widespread usage.  Its redeeming
! feature is that the exact same random numbers (to within final round-
! off error) can be generated from machine to machine.  Thus, programs
! that make use of random numbers can be easily transported to and
! checked in a new environment.
!
!      The random numbers are generated by the linear congruential
! method described, e.g., by Knuth in Seminumerical Methods (p.9),
! Addison-Wesley, 1969.  Given the I-th number of a pseudo-random
! sequence, the I+1 -st number is generated from
!             X(I+1) = (A*X(I) + C) MOD M,
! where here M = 2**22 = 4194304, C = 1731 and several suitable values
! of the multiplier A are discussed below.  Both the multiplier A and
! random number X are represented in double precision as two 11-bit
! words.  The constants are chosen so that the period is the maximum
! possible, 4194304.
!
!      In order that the same numbers be generated from machine to
! machine, it is necessary that 23-bit integers be reducible modulo
! 2**11 exactly, that 23-bit integers be added exactly, and that 11-bit
! integers be multiplied exactly.  Furthermore, if the restart option
! is used (where R is between 0 and 1), then the product R*2**22 =
! R*4194304 must be correct to the nearest integer.
!
!      The first four random numbers should be .0004127026,
! .6750836372, .1614754200, and .9086198807.  The tenth random number
! is .5527787209, and the hundredth is .3600893021 .  The thousandth
! number should be .2176990509 .
!
!      In order to generate several effectively independent sequences
! with the same generator, it is necessary to know the random number
! for several widely spaced calls.  The I-th random number times 2**22,
! where I=K*P/8 and P is the period of the sequence (P = 2**22), is
! still of the form L*P/8.  In particular we find the I-th random
! number multiplied by 2**22 is given by
! I   =  0  1*P/8  2*P/8  3*P/8  4*P/8  5*P/8  6*P/8  7*P/8  8*P/8
! RAND=  0  5*P/8  2*P/8  7*P/8  4*P/8  1*P/8  6*P/8  3*P/8  0
! Thus the 4*P/8 = 2097152 random number is 2097152/2**22.
!
!      Several multipliers have been subjected to the spectral test
! (see Knuth, p. 82).  Four suitable multipliers roughly in order of
! goodness according to the spectral test are
!    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5
!    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5
!    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5
!    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1
!
!      In the table below LOG10(NU(I)) gives roughly the number of
! random decimal digits in the random numbers considered I at a time.
! C is the primary measure of goodness.  In both cases bigger is better.
!
!                   LOG10 NU(I)              C(I)
!       A       I=2  I=3  I=4  I=5    I=2  I=3  I=4  I=5
!
!    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
!    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
!    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
!    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
!   Best
!    Possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
!
!             Input Argument --
! R      If R=0., the next random number of the sequence is generated.
!        If R .LT. 0., the last generated number will be returned for
!          possible use in a restart procedure.
!        If R .GT. 0., the sequence of random numbers will start with
!          the seed R mod 1.  This seed is also returned as the value of
!          RAND provided the arithmetic is done exactly.
!
!             Output Value --
! RAND   a pseudo-random number between 0. and 1.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  RAND
      SAVE IA1, IA0, IA1MA0, IC, IX1, IX0
      DATA IA1, IA0, IA1MA0 /1536, 1029, 507/
      DATA IC /1731/
      DATA IX1, IX0 /0, 0/
!***FIRST EXECUTABLE STATEMENT  RAND
!
!           A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
!                   + IA0*IX0) + IA0*IX0
!
      IF (R.EQ.0.) THEN
       IY0 = IA0*IX0
       IY1 = IA1*IX1 + IA1MA0*(IX0-IX1) + IY0
       IY0 = IY0 + IC
       IX0 = MOD (IY0, 2048)
       IY1 = IY1 + (IY0-IX0)/2048
       IX1 = MOD (IY1, 2048)
      ENDIF

      IF (R.GT.0.) THEN
       IX1 = MOD(R,1.)*4194304. + 0.5
       IX0 = MOD (IX1, 2048)
       IX1 = (IX1-IX0)/2048
      ENDIF

      ALEAS = IX1*2048 + IX0
      ALEAS = ALEAS / 4194304.
      RETURN

      END


