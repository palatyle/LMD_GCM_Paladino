!
! $Id: hines_gwd.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE HINES_GWD(NLON,NLEV,DTIME,paphm1x, papm1x,
     I      rlat,tx,ux,vx,
     O      zustrhi,zvstrhi,
     O      d_t_hin, d_u_hin, d_v_hin)

C ########################################################################
C Parametrization of the momentum flux deposition due to a broad band 
C spectrum of gravity waves, following Hines (1997a,b), as coded by 
C McLANDRESS (1995). Modified by McFARLANE and MANZINI (1995-1997) 
C                 MAECHAM model stand alone version
C ########################################################################
C 
C
         USE dimphy
         implicit none

cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOEGWD.h"
#include "YOMCST.h"

      INTEGER NAZMTH
      PARAMETER(NAZMTH=8)

C     INPUT ARGUMENTS.
C     ----- ----------
C
C  - 2D
C  PAPHM1   : HALF LEVEL PRESSURE (T-DT)
C  PAPM1    : FULL LEVEL PRESSURE (T-DT)
C  PTM1     : TEMPERATURE (T-DT)
C  PUM1     : ZONAL WIND (T-DT)
C  PVM1     : MERIDIONAL WIND (T-DT)
C

C     REFERENCE.
C     ----------
C         SEE MODEL DOCUMENTATION
C
C     AUTHOR.
C     -------
C
C      N. MCFARLANE   DKRZ-HAMBURG   MAY 1995
C      STAND ALONE E. MANZINI MPI-HAMBURG FEBRUARY 1997
C
C      BASED ON A COMBINATION OF THE OROGRAPHIC SCHEME BY N.MCFARLANE 1987
C      AND THE HINES SCHEME AS CODED BY C. MCLANDRESS 1995.                       
C
C
C
cym      INTEGER KLEVM1
C
      REAL PAPHM1(klon,klev+1), PAPM1(klon,KLEV)  
      REAL PTM1(klon,KLEV), PUM1(klon,KLEV), PVM1(klon,KLEV)
      REAL PRFLUX(klon)
C1
C1
C1
      REAL RLAT(klon),COSLAT(KLON)
C 
      REAL TH(klon,KLEV),
     2     UTENDGW(klon,KLEV), VTENDGW(klon,KLEV), 
     3     PRESSG(klon),
     4     UHS(klon,KLEV),     VHS(klon,KLEV), ZPR(klon)

C     * VERTICAL POSITIONING ARRAYS.

      REAL SGJ(klon,KLEV),     SHJ(klon,KLEV),    
     1     SHXKJ(klon,KLEV),   DSGJ(klon,KLEV)

C     * LOGICAL SWITCHES TO CONTROL ROOF DRAG, ENVELOP GW DRAG AND
C     * HINES' DOPPLER SPREADING EXTROWAVE GW DRAG.
C     * LOZPR IS TRUE FOR ZPR ENHANCEMENT


C     * WORK ARRAYS.

      REAL M_ALPHA(klon,KLEV,NAZMTH),     V_ALPHA(klon,KLEV,NAZMTH),
     1     SIGMA_ALPHA(klon,KLEV,NAZMTH), 
     1     SIGSQH_ALPHA(klon,KLEV,NAZMTH),
     2     DRAG_U(klon,KLEV),   DRAG_V(klon,KLEV),  FLUX_U(klon,KLEV),
     3     FLUX_V(klon,KLEV),   HEAT(klon,KLEV),    DIFFCO(klon,KLEV),
     4     BVFREQ(klon,KLEV),   DENSITY(klon,KLEV), SIGMA_T(klon,KLEV),
     5     VISC_MOL(klon,KLEV), ALT(klon,KLEV),      
     6     SIGSQMCW(klon,KLEV,NAZMTH), 
     6     SIGMATM(klon,KLEV), 
     7     AK_ALPHA(klon,NAZMTH),       K_ALPHA(klon,NAZMTH),
     8     MMIN_ALPHA(klon,NAZMTH),     I_ALPHA(klon,NAZMTH),
     9     RMSWIND(klon), BVFBOT(klon), DENSBOT(klon)
      REAL  SMOOTHR1(klon,KLEV), SMOOTHR2(klon,KLEV)
      REAL  SIGALPMC(klon,KLEV,NAZMTH)      
      REAL  F2MOD(klon,KLEV)
      
C     * THES ARE THE INPUT PARAMETERS FOR HINES ROUTINE AND
C     * ARE SPECIFIED IN ROUTINE HINES_SETUP. SINCE THIS IS CALLED
C     * ONLY AT FIRST CALL TO THIS ROUTINE THESE VARIABLES MUST BE SAVED
C     * FOR USE AT SUBSEQUENT CALLS. THIS CAN BE AVOIDED BY CALLING
C     * HINES_SETUP IN MAIN PROGRAM AND PASSING THE PARAMETERS AS
C     * SUBROUTINE ARGUEMENTS.
C

      REAL    RMSCON
      INTEGER NMESSG, IPRINT, ILRMS
      INTEGER IFL
C
      INTEGER  NAZ,ICUTOFF,NSMAX,IHEATCAL
      REAL  SLOPE,F1,F2,F3,F5,F6,KSTAR(KLON),ALT_CUTOFF,SMCO
C
C    PROVIDED AS INPUT
C
      integer nlon,nlev

      real dtime
      real paphm1x(nlon,nlev+1), papm1x(nlon,nlev)
      real ux(nlon,nlev), vx(nlon,nlev), tx(nlon,nlev)
c
c     VARIABLES FOR OUTPUT
c

      real d_t_hin(nlon,nlev),d_u_hin(nlon,nlev),d_v_hin(nlon,nlev)
      real zustrhi(nlon),zvstrhi(nlon) 

C
C     * LOGICAL SWITCHES TO CONTROL PRECIP ENHANCEMENT AND
C     * HINES' DOPPLER SPREADING EXTROWAVE GW DRAG.
C     * LOZPR IS TRUE FOR ZPR ENHANCEMENT
C  
      LOGICAL LOZPR, LORMS(klon)
C
C  LOCAL PARAMETERS TO MAKE THINGS WORK (TEMPORARY VARIABLE)

      REAL RHOH2O,ZPCONS,RGOCP,ZLAT,DTTDSF,RATIO,HSCAL
      INTEGER I,J,L,JL,JK,LE,LREF,LREFP,LEVBOT
C
C  DATA PARAMETERS NEEDED, EXPLAINED LATER

      REAL V0,VMIN,DMPSCAL,TAUFAC,HMIN,APIBT,CPART,FCRIT
      REAL PCRIT,PCONS
      INTEGER IPLEV,IERROR
   
C
C      
C     PRINT *,' IT IS STARTED HINES GOING ON...'
C
C
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
C
      d_t_hin(:,:)=0.
      
      RHOH2O=1000.    
      ZPCONS = (1000.*86400.)/RHOH2O
cym      KLEVM1=KLEV-1
C

        do jl=kidia,kfdia
        PAPHM1(JL,1) = paphm1x(JL,klev+1)
          do jk=1,klev      
          le=klev+1-jk
          PAPHM1(JL,JK+1) =  paphm1x(JL,le) 
          PAPM1(JL,JK) = papm1x(JL,le)
          PTM1(JL,JK) = tx(JL,le)
          PUM1(JL,JK) = ux(JL,le)
          PVM1(JL,JK) = vx(JL,le)
          enddo
        enddo
C
  100 CONTINUE
C
C    Define constants and arrays needed for the ccc/mam gwd scheme
C    *Constants:

      RGOCP=RD/RCPD
      LREFP=KLEV-1
      LREF=KLEV-2
C1
C1    *Arrays
C1
      DO 2101 JK=1,KLEV
      DO 2102 JL=KIDIA,KFDIA
      SHJ(JL,JK)=PAPM1(JL,JK)/PAPHM1(JL,klev+1)
      SGJ(JL,JK)=PAPM1(JL,JK)/PAPHM1(JL,klev+1)
      DSGJ(JL,JK)=(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))/PAPHM1(JL,klev+1)
      SHXKJ(JL,JK)=(PAPM1(JL,JK)/PAPHM1(JL,klev+1))**RGOCP 
      TH(JL,JK)= PTM1(JL,JK)
 2102 CONTINUE
 2101 CONTINUE    
      
CC
      DO 211 JL=KIDIA,KFDIA
      PRESSG(JL)=PAPHM1(JL,klev+1)
  211 CONTINUE
C
C
      DO 301 JL=KIDIA,KFDIA
      PRFLUX(JL) = 0.0
      ZPR(JL)=ZPCONS*PRFLUX(JL)
      ZLAT=(RLAT(JL)/180.)*RPI
      COSLAT(Jl)=COS(ZLAT)    
  301 CONTINUE
C
C
  400 CONTINUE
C  
C
C
C
*/#########################################################################
*/
*/
C
C     * AUG. 14/95 - C. MCLANDRESS.
C     * SEP.    95   N. MCFARLANE.
C
C     * THIS ROUTINE CALCULATES THE HORIZONTAL WIND TENDENCIES
C     * DUE TO MCFARLANE'S OROGRAPHIC GW DRAG SCHEME, HINES'
C     * DOPPLER SPREAD SCHEME FOR "EXTROWAVES" AND ADDS ON
C     * ROOF DRAG. IT IS BASED ON THE ROUTINE GWDFLX8.
C
C     * LREFP IS THE INDEX OF THE MODEL LEVEL BELOW THE REFERENCE LEVEL
C     * I/O ARRAYS PASSED FROM MAIN.
C     * (PRESSG = SURFACE PRESSURE)
C
C
C
C
C     * CONSTANTS VALUES DEFINED IN DATA STATEMENT ARE :
C     * VMIN     = MIMINUM WIND IN THE DIRECTION OF REFERENCE LEVEL
C     *            WIND BEFORE WE CONSIDER BREAKING TO HAVE OCCURED.
C     * DMPSCAL  = DAMPING TIME FOR GW DRAG IN SECONDS.
C     * TAUFAC   = 1/(LENGTH SCALE).
C     * HMIN     = MIMINUM ENVELOPE HEIGHT REQUIRED TO PRODUCE GW DRAG.
C     * V0       = VALUE OF WIND THAT APPROXIMATES ZERO.


      DATA    VMIN  /    5.0 /, V0       / 1.E-10 /,
     1        TAUFAC/  5.E-6 /, HMIN     /   40000. /,
     3        DMPSCAL  / 6.5E+6 /, APIBT / 1.5708 /,
     4        CPART /    0.7 /, FCRIT    / 1. /

C     * HINES EXTROWAVE GWD CONSTANTS DEFINED IN DATA STATEMENT ARE:
C     * RMSCON = ROOT MEAN SQUARE GRAVITY WAVE WIND AT LOWEST LEVEL (M/S).
C     * NMESSG  = UNIT NUMBER FOR PRINTED MESSAGES.
C     * IPRINT  = 1 TO DO PRINT OUT SOME HINES ARRAYS.
C     * IFL     = FIRST CALL FLAG TO HINES_SETUP ("SAVE" IT)
C     * PCRIT = CRITICAL VALUE OF ZPR (MM/D)
C     * IPLEV = LEVEL OF APPLICATION OF PRCIT
C     * PCONS = FACTOR OF ZPR ENHANCEMENT
C

      DATA PCRIT / 5. /, PCONS / 4.75 /

      IPLEV = LREFP-1
C
      DATA    RMSCON  / 1.00 /
     1        IPRINT   /  2  /, NMESSG  /   6   /
      DATA    IFL / 0 /
C
      LOZPR = .FALSE.
C
C-----------------------------------------------------------------------
C
C
C
C     * SET ERROR FLAG

      IERROR = 0

C     * SPECIFY VARIOUS PARAMETERS FOR HINES ROUTINE AT VERY FIRST CALL.
C     * (NOTE THAT ARRAY K_ALPHA IS SPECIFIED SO MAKE SURE THAT
C     * IT IS NOT OVERWRITTEN LATER ON).
C
        CALL HINES_SETUP (NAZ,SLOPE,F1,F2,F3,F5,F6,KSTAR,
     1                    ICUTOFF,ALT_CUTOFF,SMCO,NSMAX,IHEATCAL,
     2                   K_ALPHA,IERROR,NMESSG,klon,NAZMTH,COSLAT)
        IF (IERROR.NE.0)  GO TO 999
C
C     * START GWD CALCULATIONS.

      LREF  = LREFP-1

C
      DO 105 J=1,NAZMTH
      DO 105 L=1,KLEV
      DO 105 I=kidia,klon
        SIGSQMCW(I,L,J) = 0.
  105 CONTINUE
c


C     * INITIALIZE NECESSARY ARRAYS.
C
      DO 120 L=1,KLEV
      DO 120 I=KIDIA,KFDIA
        UTENDGW(I,L) = 0.
        VTENDGW(I,L) = 0.

        UHS(I,L) = 0.
        VHS(I,L) = 0.

 120  CONTINUE
C
C     * IF USING HINES SCHEME THEN CALCULATE B V FREQUENCY AT ALL POINTS
C     * AND SMOOTH BVFREQ.

        DO 130 L=2,KLEV
        DO 130 I=KIDIA,KFDIA
          DTTDSF=(TH(I,L)/SHXKJ(I,L)-TH(I,L-1)/
     1            SHXKJ(I,L-1))/(SHJ(I,L)-SHJ(I,L-1))
          DTTDSF=MIN(DTTDSF, -5./SGJ(I,L))
          BVFREQ(I,L)=SQRT(-DTTDSF*SGJ(I,L)*(SGJ(I,L)**RGOCP)/RD)
     1                     *RG/PTM1(I,L)
  130   CONTINUE
        DO 135 L=1,KLEV
        DO 135 I=KIDIA,KFDIA
          IF(L.EQ.1)                        THEN
            BVFREQ(I,L) = BVFREQ(I,L+1)
          ENDIF
          IF(L.GT.1)                        THEN
            RATIO=5.*LOG(SGJ(I,L)/SGJ(I,L-1))
            BVFREQ(I,L) = (BVFREQ(I,L-1) + RATIO*BVFREQ(I,L))
     1                       /(1.+RATIO)
          ENDIF
  135   CONTINUE
C
C
  300 CONTINUE

C     * CALCULATE GW DRAG DUE TO HINES' EXTROWAVES
C     * SET MOLECULAR VISCOSITY TO A VERY SMALL VALUE.
C     * IF THE MODEL TOP IS GREATER THAN 100 KM THEN THE ACTUAL
C     * VISCOSITY COEFFICIENT COULD BE SPECIFIED HERE.

      DO 310 L=1,KLEV
      DO 310 I=KIDIA,KFDIA
         VISC_MOL(I,L) = 1.5E-5
         DRAG_U(I,L) = 0.
         DRAG_V(I,L) = 0.
         FLUX_U(I,L) = 0.
         FLUX_V(I,L) = 0.
         HEAT(I,L)   = 0.
         DIFFCO(I,L) = 0.
 310  CONTINUE

C     * ALTITUDE AND DENSITY AT BOTTOM.

      DO 330 I=KIDIA,KFDIA
         HSCAL = RD * PTM1(I,KLEV) / RG
         DENSITY(I,KLEV) = SGJ(I,KLEV) * PRESSG(I) / (RG*HSCAL)
         ALT(I,KLEV) = 0.
  330 CONTINUE

C     * ALTITUDE AND DENSITY AT REMAINING LEVELS.

      DO 340 L=KLEV-1,1,-1
      DO 340 I=KIDIA,KFDIA
         HSCAL = RD * PTM1(I,L) / RG
         ALT(I,L) = ALT(I,L+1) + HSCAL * DSGJ(I,L) / SGJ(I,L)
         DENSITY(I,L) = SGJ(I,L) * PRESSG(I) / (RG*HSCAL)
  340 CONTINUE

C
C     * INITIALIZE SWITCHES FOR HINES GWD CALCULATION
C
      ILRMS = 0
C
      DO 345 I=KIDIA,KFDIA 
      LORMS(I) = .FALSE.
  345 CONTINUE 
C
C
C     * DEFILE BOTTOM LAUNCH LEVEL
C
      LEVBOT = IPLEV
C
C     * BACKGROUND WIND MINUS VALUE AT BOTTOM LAUNCH LEVEL.
C
      DO 351 L=1,LEVBOT
      DO 351 I=KIDIA,KFDIA 
      UHS(I,L) = PUM1(I,L) - PUM1(I,LEVBOT)
      VHS(I,L) = PVM1(I,L) - PVM1(I,LEVBOT)
  351 CONTINUE
C
C     * SPECIFY ROOT MEAN SQUARE WIND AT BOTTOM LAUNCH LEVEL.
C
       DO 355 I=KIDIA,KFDIA 
       RMSWIND(I) = RMSCON
  355  CONTINUE

      IF (LOZPR) THEN
        DO 350 I=KIDIA,KFDIA 
        IF (ZPR(I) .GT. PCRIT) THEN
          RMSWIND(I) = RMSCON
     >                +( (ZPR(I)-PCRIT)/ZPR(I) )*PCONS
        ENDIF
  350   CONTINUE
      ENDIF
C
      DO 356 I=KIDIA,KFDIA 
      IF (RMSWIND(I) .GT. 0.0) THEN
      ILRMS = ILRMS+1
      LORMS(I) = .TRUE.
      ENDIF
  356 CONTINUE
C
C     * CALCULATE GWD (NOTE THAT DIFFUSION COEFFICIENT AND
C     * HEATING RATE ONLY CALCULATED IF IHEATCAL = 1).
C
      IF ( ILRMS.GT.0 )       THEN                    
C
      CALL HINES_EXTRO0 (DRAG_U,DRAG_V,HEAT,DIFFCO,FLUX_U,FLUX_V,
     1                   UHS,VHS,BVFREQ,DENSITY,VISC_MOL,ALT,
     2                   RMSWIND,K_ALPHA,M_ALPHA,V_ALPHA,
     3                   SIGMA_ALPHA,SIGSQH_ALPHA,AK_ALPHA,
     4                   MMIN_ALPHA,I_ALPHA,SIGMA_T,DENSBOT,BVFBOT,
     5                   1,IHEATCAL,ICUTOFF,IPRINT,NSMAX,
     6                   SMCO,ALT_CUTOFF,KSTAR,SLOPE,
     7                   F1,F2,F3,F5,F6,NAZ,SIGSQMCW,SIGMATM,
     8                   KIDIA,klon,1,LEVBOT,KLON,KLEV,NAZMTH,
     9                   LORMS,SMOOTHR1,SMOOTHR2,
     9                   SIGALPMC,F2MOD)

C     * ADD ON HINES' GWD TENDENCIES TO OROGRAPHIC TENDENCIES AND
C     * APPLY HINES' GW DRAG ON (UROW,VROW) WORK ARRAYS.

      DO 360 L=1,KLEV
      DO 360 I=KIDIA,KFDIA
         UTENDGW(I,L) = UTENDGW(I,L) + DRAG_U(I,L)
         VTENDGW(I,L) = VTENDGW(I,L) + DRAG_V(I,L)
  360 CONTINUE
C

C     * END OF HINES CALCULATIONS.
C
      ENDIF
C
  500 CONTINUE


C-----------------------------------------------------------------------
C
        do jl=kidia,kfdia
          zustrhi(jl)=FLUX_U(jl,1)
          zvstrhi(jl)=FLUX_v(jl,1)
          do jk=1,klev
          le=klev-jk+1
          d_u_hin(jl,JK) =  UTENDGW(jl,le) * dtime
          d_v_hin(jl,JK) =  VTENDGW(jl,le) * dtime
          enddo
        enddo

c     PRINT *,'UTENDGW:',UTENDGW

C     PRINT *,' HINES HAS BEEN COMPLETED (LONG ISNT IT...)'

      RETURN
 999  CONTINUE

C     * IF ERROR DETECTED THEN ABORT.

      WRITE (NMESSG,6000)
      WRITE (NMESSG,6010) IERROR
 6000 FORMAT (/' EXECUTION ABORTED IN GWDOREXV')
 6010 FORMAT ('     ERROR FLAG =',I4)

C
      RETURN
      END
*/
*/


      SUBROUTINE HINES_EXTRO0 (DRAG_U,DRAG_V,HEAT,DIFFCO,FLUX_U,FLUX_V,
     1                         VEL_U,VEL_V,BVFREQ,DENSITY,VISC_MOL,ALT,
     2                         RMSWIND,K_ALPHA,M_ALPHA,V_ALPHA,
     3                         SIGMA_ALPHA,SIGSQH_ALPHA,AK_ALPHA,
     4                         MMIN_ALPHA,I_ALPHA,SIGMA_T,DENSB,BVFB,
     5                         IORDER,IHEATCAL,ICUTOFF,IPRINT,NSMAX,
     6                         SMCO,ALT_CUTOFF,KSTAR,SLOPE,
     7                         F1,F2,F3,F5,F6,NAZ,SIGSQMCW,SIGMATM,
     8                         IL1,IL2,LEV1,LEV2,NLONS,NLEVS,NAZMTH,
     9                         LORMS,SMOOTHR1,SMOOTHR2,
     9                         SIGALPMC,F2MOD)

       implicit none
C
C  Main routine for Hines' "extrowave" gravity wave parameterization based
C  on Hines' Doppler spread theory. This routine calculates zonal
C  and meridional components of gravity wave drag, heating rates
C  and diffusion coefficient on a longitude by altitude grid.
C  No "mythical" lower boundary region calculation is made so it
C  is assumed that lowest level winds are weak (i.e, approximately zero).
C
C  Aug. 13/95 - C. McLandress
C  SEPT. /95  - N.McFarlane
C
C  Modifications:
C
C  Output arguements:
C
C     * DRAG_U = zonal component of gravity wave drag (m/s^2).
C     * DRAG_V = meridional component of gravity wave drag (m/s^2).
C     * HEAT   = gravity wave heating (K/sec).
C     * DIFFCO = diffusion coefficient (m^2/sec)
C     * FLUX_U = zonal component of vertical momentum flux (Pascals)
C     * FLUX_V = meridional component of vertical momentum flux (Pascals)
C
C  Input arguements:
C
C     * VEL_U      = background zonal wind component (m/s).
C     * VEL_V      = background meridional wind component (m/s).
C     * BVFREQ     = background Brunt Vassala frequency (radians/sec).
C     * DENSITY    = background density (kg/m^3) 
C     * VISC_MOL   = molecular viscosity (m^2/s)
C     * ALT        = altitude of momentum, density, buoyancy levels (m)
C     *              (NOTE: levels ordered so that ALT(I,1) > ALT(I,2), etc.)
C     * RMSWIND   = root mean square gravity wave wind at lowest level (m/s).
C     * K_ALPHA    = horizontal wavenumber of each azimuth (1/m).
C     * IORDER	   = 1 means vertical levels are indexed from top down 
C     *              (i.e., highest level indexed 1 and lowest level NLEVS);
C     *           .NE. 1 highest level is index NLEVS.
C     * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
C     * IPRINT     = 1 to print out various arrays.
C     * ICUTOFF    = 1 to exponentially damp GWD, heating and diffusion 
C     *              arrays above ALT_CUTOFF; otherwise arrays not modified.
C     * ALT_CUTOFF = altitude in meters above which exponential decay applied.
C     * SMCO       = smoothing factor used to smooth cutoff vertical 
C     *              wavenumbers and total rms winds in vertical direction
C     *              before calculating drag or heating
C     *              (SMCO >= 1 ==> 1:SMCO:1 stencil used).
C     * NSMAX      = number of times smoother applied ( >= 1),
C     *            = 0 means no smoothing performed.
C     * KSTAR      = typical gravity wave horizontal wavenumber (1/m).
C     * SLOPE      = slope of incident vertical wavenumber spectrum
C     *              (SLOPE must equal 1., 1.5 or 2.).
C     * F1 to F6   = Hines's fudge factors (F4 not needed since used for
C     *              vertical flux of vertical momentum).
C     * NAZ        = actual number of horizontal azimuths used.
C     * IL1        = first longitudinal index to use (IL1 >= 1).
C     * IL2        = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1       = index of first level for drag calculation.
C     * LEV2       = index of last level for drag calculation 
C     *              (i.e., LEV1 < LEV2 <= NLEVS).
C     * NLONS      = number of longitudes.
C     * NLEVS      = number of vertical levels.
C     * NAZMTH     = azimuthal array dimension (NAZMTH >= NAZ).
C 
C  Work arrays.
C
C     * M_ALPHA      = cutoff vertical wavenumber (1/m).
C     * V_ALPHA      = wind component at each azimuth (m/s) and if IHEATCAL=1
C     *                holds vertical derivative of cutoff wavenumber.
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * SIGMA_T      = total rms horizontal wind (m/s).
C     * AK_ALPHA     = spectral amplitude factor at each azimuth 
C     *                (i.e.,{AjKj}) in m^4/s^2.
C     * I_ALPHA      = Hines' integral.
C     * MMIN_ALPHA   = minimum value of cutoff wavenumber.
C     * DENSB        = background density at bottom level.
C     * BVFB         = buoyancy frequency at bottom level and
C     *                work array for ICUTOFF = 1.
C
C     * LORMS       = .TRUE. for drag computation 
C
      INTEGER  NAZ, NLONS, NLEVS, NAZMTH, IL1, IL2, LEV1, LEV2
      INTEGER  ICUTOFF, NSMAX, IORDER, IHEATCAL, IPRINT
      REAL  KSTAR(NLONS), F1, F2, F3, F5, F6, SLOPE
      REAL  ALT_CUTOFF, SMCO
      REAL  DRAG_U(NLONS,NLEVS),   DRAG_V(NLONS,NLEVS) 
      REAL  HEAT(NLONS,NLEVS),     DIFFCO(NLONS,NLEVS)
      REAL  FLUX_U(NLONS,NLEVS),   FLUX_V(NLONS,NLEVS)
      REAL  VEL_U(NLONS,NLEVS),    VEL_V(NLONS,NLEVS)
      REAL  BVFREQ(NLONS,NLEVS),   DENSITY(NLONS,NLEVS)
      REAL  VISC_MOL(NLONS,NLEVS), ALT(NLONS,NLEVS)
      REAL  RMSWIND(NLONS),       BVFB(NLONS),   DENSB(NLONS)
      REAL  SIGMA_T(NLONS,NLEVS), SIGSQMCW(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH), SIGMATM(NLONS,NLEVS)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH), V_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  AK_ALPHA(NLONS,NAZMTH),      K_ALPHA(NLONS,NAZMTH)
      REAL  MMIN_ALPHA(NLONS,NAZMTH) ,   I_ALPHA(NLONS,NAZMTH)
      REAL  SMOOTHR1(NLONS,NLEVS), SMOOTHR2(NLONS,NLEVS)
      REAL  SIGALPMC(NLONS,NLEVS,NAZMTH)
      REAL  F2MOD(NLONS,NLEVS)
C
      LOGICAL LORMS(NLONS)
C
C  Internal variables.
C
      INTEGER  LEVBOT, LEVTOP, I, N, L, LEV1P, LEV2M
      INTEGER  ILPRT1, ILPRT2
C----------------------------------------------------------------------- 
C
C     PRINT *,' IN HINES_EXTRO0'
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
C
C  Index of lowest altitude level (bottom of drag calculation).
C
      LEVBOT = LEV2
      LEVTOP = LEV1
      IF (IORDER.NE.1)  THEN
      write(6,1)
   1  format(2x,' error: IORDER NOT ONE! ')
      END IF
C
C  Buoyancy and density at bottom level.
C
      DO 10 I = IL1,IL2
        BVFB(I)  = BVFREQ(I,LEVBOT)
        DENSB(I) = DENSITY(I,LEVBOT)
 10   CONTINUE
C
C  initialize some variables
C
      DO 20 N = 1,NAZ
      DO 20 L=LEV1,LEV2
      DO 20 I=IL1,IL2
      M_ALPHA(I,L,N) = 0.0
  20  CONTINUE
      DO 21 L=LEV1,LEV2
      DO 21 I=IL1,IL2
      SIGMA_T(I,L) = 0.0
  21  CONTINUE
      DO 22 N = 1,NAZ
      DO 22 I=IL1,IL2
      I_ALPHA(I,N) = 0.0
  22  CONTINUE 
C
C  Compute azimuthal wind components from zonal and meridional winds.
C
      CALL HINES_WIND ( V_ALPHA, 
     ^                  VEL_U, VEL_V, NAZ,
     ^                  IL1, IL2, LEV1, LEV2, NLONS, NLEVS, NAZMTH )
C
C  Calculate cutoff vertical wavenumber and velocity variances.
C
      CALL HINES_WAVNUM ( M_ALPHA, SIGMA_ALPHA, SIGSQH_ALPHA, SIGMA_T,
     ^                    AK_ALPHA, V_ALPHA, VISC_MOL, DENSITY, DENSB,
     ^                    BVFREQ, BVFB, RMSWIND, I_ALPHA, MMIN_ALPHA,
     ^                    KSTAR, SLOPE, F1, F2, F3, NAZ, LEVBOT,
     ^                    LEVTOP,IL1,IL2,NLONS,NLEVS,NAZMTH, SIGSQMCW,
     ^                    SIGMATM,LORMS,SIGALPMC,F2MOD)
C
C  Smooth cutoff wavenumbers and total rms velocity in the vertical 
C  direction NSMAX times, using FLUX_U as temporary work array.
C   
      IF (NSMAX.GT.0)  THEN
        DO 80 N = 1,NAZ
          DO 81 L=LEV1,LEV2
          DO 81 I=IL1,IL2
          SMOOTHR1(I,L) = M_ALPHA(I,L,N)
 81       CONTINUE 
             CALL VERT_SMOOTH (SMOOTHR1, 
     ^                       SMOOTHR2, SMCO, NSMAX,
     ^                       IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
        DO 83 L=LEV1,LEV2
        DO 83 I=IL1,IL2
        M_ALPHA(I,L,N) = SMOOTHR1(I,L)
 83     CONTINUE
 80     CONTINUE
        CALL VERT_SMOOTH ( SIGMA_T, 
     ^                     SMOOTHR2, SMCO, NSMAX,
     ^                     IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
      END IF
C
C  Calculate zonal and meridional components of the
C  momentum flux and drag.
C
      CALL HINES_FLUX ( FLUX_U, FLUX_V, DRAG_U, DRAG_V, 
     ^                  ALT, DENSITY, DENSB, M_ALPHA, 
     ^                  AK_ALPHA, K_ALPHA, SLOPE, NAZ,
     ^                  IL1, IL2, LEV1, LEV2, NLONS, NLEVS, NAZMTH,
     ^                  LORMS)
C
C  Cutoff drag above ALT_CUTOFF, using BVFB as temporary work array.
C
      IF (ICUTOFF.EQ.1)  THEN		
        CALL HINES_EXP ( DRAG_U, 
     ^                   BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                   IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
        CALL HINES_EXP ( DRAG_V, 
     ^                   BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                   IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
      END IF   
C
C  Print out various arrays for diagnostic purposes.
C
      IF (IPRINT.EQ.1)  THEN
        ILPRT1 = 15
        ILPRT2 = 16
        CALL HINES_PRINT ( FLUX_U, FLUX_V, DRAG_U, DRAG_V, ALT,
     ^                     SIGMA_T, SIGMA_ALPHA, V_ALPHA, M_ALPHA,
     ^                     1, 1, 6, ILPRT1, ILPRT2, LEV1, LEV2,
     ^                     NAZ, NLONS, NLEVS, NAZMTH)
      END IF
C
C  If not calculating heating rate and diffusion coefficient then finished.
C
      IF (IHEATCAL.NE.1)  RETURN
C
C  Calculate vertical derivative of cutoff wavenumber (store
C  in array V_ALPHA) using centered differences at interior gridpoints
C  and one-sided differences at first and last levels.
C 
      DO 130 N = 1,NAZ
        DO 100 L = LEV1P,LEV2M
          DO 90 I = IL1,IL2
            V_ALPHA(I,L,N) = ( M_ALPHA(I,L+1,N) - M_ALPHA(I,L-1,N) ) 
     ^                       / ( ALT(I,L+1) - ALT(I,L-1) )
  90      CONTINUE
 100    CONTINUE
        DO 110 I = IL1,IL2
          V_ALPHA(I,LEV1,N) = ( M_ALPHA(I,LEV1P,N) - M_ALPHA(I,LEV1,N) ) 
     ^                       / ( ALT(I,LEV1P) - ALT(I,LEV1) )
 110    CONTINUE
        DO 120 I = IL1,IL2
          V_ALPHA(I,LEV2,N) = ( M_ALPHA(I,LEV2,N) - M_ALPHA(I,LEV2M,N) ) 
     ^                       / ( ALT(I,LEV2) - ALT(I,LEV2M) )
 120    CONTINUE
 130  CONTINUE
C
C  Heating rate and diffusion coefficient.
C
      CALL HINES_HEAT ( HEAT, DIFFCO, 
     ^                  M_ALPHA, V_ALPHA, AK_ALPHA, K_ALPHA, 
     ^                  BVFREQ, DENSITY, DENSB, SIGMA_T, VISC_MOL, 
     ^                  KSTAR, SLOPE, F2, F3, F5, F6, NAZ, 
     ^                  IL1, IL2, LEV1, LEV2, NLONS, NLEVS, NAZMTH)
C
C  Finished.
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_WAVNUM (M_ALPHA,SIGMA_ALPHA,SIGSQH_ALPHA,SIGMA_T,
     1                         AK_ALPHA,V_ALPHA,VISC_MOL,DENSITY,DENSB,
     2                         BVFREQ,BVFB,RMS_WIND,I_ALPHA,MMIN_ALPHA,
     3                         KSTAR,SLOPE,F1,F2,F3,NAZ,LEVBOT,LEVTOP,
     4                         IL1,IL2,NLONS,NLEVS,NAZMTH,SIGSQMCW,
     5                         SIGMATM,LORMS,SIGALPMC,F2MOD)
C
C  This routine calculates the cutoff vertical wavenumber and velocity
C  variances on a longitude by altitude grid for the Hines' Doppler 
C  spread gravity wave drag parameterization scheme.
C  NOTE: (1) only values of four or eight can be used for # azimuths (NAZ).
C        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (SLOPE). 
C
C  Aug. 10/95 - C. McLandress
C
C  Output arguements:
C
C     * M_ALPHA      = cutoff wavenumber at each azimuth (1/m).
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * SIGMA_T      = total rms horizontal wind (m/s).
C     * AK_ALPHA     = spectral amplitude factor at each azimuth 
C     *                (i.e.,{AjKj}) in m^4/s^2.
C
C  Input arguements:
C
C     * V_ALPHA  = wind component at each azimuth (m/s). 
C     * VISC_MOL = molecular viscosity (m^2/s)
C     * DENSITY  = background density (kg/m^3).
C     * DENSB    = background density at model bottom (kg/m^3).
C     * BVFREQ   = background Brunt Vassala frequency (radians/sec).
C     * BVFB     = background Brunt Vassala frequency at model bottom.
C     * RMS_WIND = root mean square gravity wave wind at lowest level (m/s).
C     * KSTAR    = typical gravity wave horizontal wavenumber (1/m).
C     * SLOPE    = slope of incident vertical wavenumber spectrum
C     *            (SLOPE = 1., 1.5 or 2.).
C     * F1,F2,F3 = Hines's fudge factors.
C     * NAZ      = actual number of horizontal azimuths used (4 or 8).
C     * LEVBOT   = index of lowest vertical level.
C     * LEVTOP   = index of highest vertical level 
C     *            (NOTE: if LEVTOP < LEVBOT then level index 
C     *             increases from top down).
C     * IL1      = first longitudinal index to use (IL1 >= 1).
C     * IL2      = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS    = number of longitudes.
C     * NLEVS    = number of vertical levels.
C     * NAZMTH   = azimuthal array dimension (NAZMTH >= NAZ).
C
C     * LORMS       = .TRUE. for drag computation 
C
C  Input work arrays:
C
C     * I_ALPHA    = Hines' integral at a single level.
C     * MMIN_ALPHA = minimum value of cutoff wavenumber.
C
      INTEGER  NAZ, LEVBOT, LEVTOP, IL1, IL2, NLONS, NLEVS, NAZMTH
      REAL  SLOPE, KSTAR(NLONS), F1, F2, F3
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGALPMC(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQMCW(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_T(NLONS,NLEVS)
      REAL  SIGMATM(NLONS,NLEVS)
      REAL  AK_ALPHA(NLONS,NAZMTH)
      REAL  V_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  VISC_MOL(NLONS,NLEVS)
      REAL  F2MOD(NLONS,NLEVS)
      REAL  DENSITY(NLONS,NLEVS),  DENSB(NLONS)
      REAL  BVFREQ(NLONS,NLEVS),   BVFB(NLONS),  RMS_WIND(NLONS)
      REAL  I_ALPHA(NLONS,NAZMTH), MMIN_ALPHA(NLONS,NAZMTH)
C
      LOGICAL LORMS(NLONS)
C
C Internal variables.
C
      INTEGER  I, L, N, LSTART, LEND, LINCR, LBELOW
      REAL  M_SUB_M_TURB, M_SUB_M_MOL, M_TRIAL
      REAL  VISC, VISC_MIN, AZFAC, SP1

cc      REAL  N_OVER_M(1000), SIGFAC(1000)

      REAL  N_OVER_M(NLONS), SIGFAC(NLONS)
      DATA  VISC_MIN / 1.E-10 / 
C-----------------------------------------------------------------------     
C

C     PRINT *,'IN HINES_WAVNUM'
      SP1 = SLOPE + 1.
C
C  Indices of levels to process.
C
      IF (LEVBOT.GT.LEVTOP)  THEN
        LSTART = LEVBOT - 1     
        LEND   = LEVTOP         
        LINCR  = -1
      ELSE
      write(6,1)
   1  format(2x,' error: IORDER NOT ONE! ')
      END IF
C
C  Use horizontal isotropy to calculate azimuthal variances at bottom level.
C
      AZFAC = 1. / REAL(NAZ)
      DO 20 N = 1,NAZ
        DO 10 I = IL1,IL2
          SIGSQH_ALPHA(I,LEVBOT,N) = AZFAC * RMS_WIND(I)**2
 10     CONTINUE
 20   CONTINUE
C
C  Velocity variances at bottom level.
C
      CALL HINES_SIGMA ( SIGMA_T, SIGMA_ALPHA, 
     ^                   SIGSQH_ALPHA, NAZ, LEVBOT, 
     ^                   IL1, IL2, NLONS, NLEVS, NAZMTH)
c
      CALL HINES_SIGMA ( SIGMATM, SIGALPMC, 
     ^                   SIGSQMCW, NAZ, LEVBOT, 
     ^                   IL1, IL2, NLONS, NLEVS, NAZMTH) 
C
C  Calculate cutoff wavenumber and spectral amplitude factor 
C  at bottom level where it is assumed that background winds vanish
C  and also initialize minimum value of cutoff wavnumber.
C
      DO 40 N = 1,NAZ
        DO 30 I = IL1,IL2
        IF (LORMS(I)) THEN
          M_ALPHA(I,LEVBOT,N) =  BVFB(I) / 
     ^                           ( F1 * SIGMA_ALPHA(I,LEVBOT,N) 
     ^                           + F2 * SIGMA_T(I,LEVBOT) )
          AK_ALPHA(I,N)   = SIGSQH_ALPHA(I,LEVBOT,N) 
     ^                      / ( M_ALPHA(I,LEVBOT,N)**SP1 / SP1 )
          MMIN_ALPHA(I,N) = M_ALPHA(I,LEVBOT,N)
        ENDIF
 30     CONTINUE
 40   CONTINUE
C
C  Calculate quantities from the bottom upwards, 
C  starting one level above bottom.
C
      DO 150 L = LSTART,LEND,LINCR
C
C  Level beneath present level.
C
        LBELOW = L - LINCR 
C
C  Calculate N/m_M where m_M is maximum permissible value of the vertical
C  wavenumber (i.e., m > m_M are obliterated) and N is buoyancy frequency.
C  m_M is taken as the smaller of the instability-induced 
C  wavenumber (M_SUB_M_TURB) and that imposed by molecular viscosity
C  (M_SUB_M_MOL). Since variance at this level is not yet known
C  use value at level below.
C
        DO 50 I = IL1,IL2
        IF (LORMS(I)) THEN
c
        F2MFAC=SIGMATM(I,LBELOW)**2
        F2MOD(I,LBELOW) =1.+ 2.*F2MFAC
     ^                      / ( F2MFAC+SIGMA_T(I,LBELOW)**2 )
c
         VISC = AMAX1 ( VISC_MOL(I,L), VISC_MIN )
         M_SUB_M_TURB = BVFREQ(I,L) 
     ^                 / ( F2 *F2MOD(I,LBELOW)*SIGMA_T(I,LBELOW))
         M_SUB_M_MOL = (BVFREQ(I,L)*KSTAR(I)/VISC)**0.33333333/F3
          IF (M_SUB_M_TURB .LT. M_SUB_M_MOL)  THEN
            N_OVER_M(I) = F2 *F2MOD(I,LBELOW)*SIGMA_T(I,LBELOW)
          ELSE
            N_OVER_M(I) = BVFREQ(I,L) / M_SUB_M_MOL 
          END IF
        ENDIF
  50    CONTINUE
C
C  Calculate cutoff wavenumber at this level.
C
        DO 70 N = 1,NAZ
          DO 60 I = IL1,IL2
          IF (LORMS(I)) THEN
C
C  Calculate trial value (since variance at this level is not yet known
C  use value at level below). If trial value is negative or if it exceeds 
C  minimum value (not permitted) then set it to minimum value. 
C                                                                      
            M_TRIAL = BVFB(I) / ( F1 * ( SIGMA_ALPHA(I,LBELOW,N)+  
     ^       SIGALPMC(I,LBELOW,N)) + N_OVER_M(I) + V_ALPHA(I,L,N) )
            IF (M_TRIAL.LE.0. .OR. M_TRIAL.GT.MMIN_ALPHA(I,N))  THEN
              M_TRIAL = MMIN_ALPHA(I,N)
            END IF
            M_ALPHA(I,L,N) = M_TRIAL
C
C  Reset minimum value of cutoff wavenumber if necessary.
C
            IF (M_ALPHA(I,L,N) .LT. MMIN_ALPHA(I,N))  THEN
              MMIN_ALPHA(I,N) = M_ALPHA(I,L,N)
            END IF
C
          ENDIF
 60       CONTINUE
 70     CONTINUE
C
C  Calculate the Hines integral at this level.
C
        CALL HINES_INTGRL ( I_ALPHA, 
     ^                      V_ALPHA, M_ALPHA, BVFB, SLOPE, NAZ, 
     ^                      L, IL1, IL2, NLONS, NLEVS, NAZMTH,
     ^                      LORMS )

C
C  Calculate the velocity variances at this level.
C
        DO 80 I = IL1,IL2
          SIGFAC(I) = DENSB(I) / DENSITY(I,L) 
     ^                * BVFREQ(I,L) / BVFB(I) 
 80     CONTINUE
        DO 100 N = 1,NAZ
          DO 90 I = IL1,IL2
            SIGSQH_ALPHA(I,L,N) = SIGFAC(I) * AK_ALPHA(I,N) 
     ^                            * I_ALPHA(I,N)
  90      CONTINUE
 100    CONTINUE
        CALL HINES_SIGMA ( SIGMA_T, SIGMA_ALPHA, 
     ^                     SIGSQH_ALPHA, NAZ, L, 
     ^                     IL1, IL2, NLONS, NLEVS, NAZMTH )
c
        CALL HINES_SIGMA ( SIGMATM, SIGALPMC, 
     ^                     SIGSQMCW, NAZ, L, 
     ^                     IL1, IL2, NLONS, NLEVS, NAZMTH )
C
C  End of level loop.
C
 150  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_WIND (V_ALPHA,VEL_U,VEL_V,
     1                       NAZ,IL1,IL2,LEV1,LEV2,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the azimuthal horizontal background wind components 
C  on a longitude by altitude grid for the case of 4 or 8 azimuths for
C  the Hines' Doppler spread GWD parameterization scheme.
C
C  Aug. 7/95 - C. McLandress
C
C  Output arguement:
C
C     * V_ALPHA   = background wind component at each azimuth (m/s). 
C     *             (note: first azimuth is in eastward direction
C     *              and rotate in counterclockwise direction.)
C
C  Input arguements:
C
C     * VEL_U     = background zonal wind component (m/s).
C     * VEL_V     = background meridional wind component (m/s).
C     * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
C     * IL1       = first longitudinal index to use (IL1 >= 1).
C     * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1      = first altitude level to use (LEV1 >=1). 
C     * LEV2      = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS     = number of longitudes.
C     * NLEVS     = number of vertical levels.
C     * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Constants in DATA statements.
C
C     * COS45 = cosine of 45 degrees. 		
C     * UMIN  = minimum allowable value for zonal or meridional 
C     *         wind component (m/s).
C
C  Subroutine arguements.
C
      INTEGER  NAZ, IL1, IL2, LEV1, LEV2
      INTEGER  NLONS, NLEVS, NAZMTH
      REAL  V_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  VEL_U(NLONS,NLEVS), VEL_V(NLONS,NLEVS)
C
C  Internal variables.
C
      INTEGER  I, L
      REAL U, V, COS45, UMIN
C
      DATA  COS45 / 0.7071068 /
      DATA  UMIN / 0.001 /
C-----------------------------------------------------------------------     
C
C  Case with 4 azimuths.
C

C      PRINT *,'IN HINES_WIND'
      IF (NAZ.EQ.4)  THEN
        DO 20 L = LEV1,LEV2
          DO 10 I = IL1,IL2
            U = VEL_U(I,L)
            V = VEL_V(I,L)
            IF (ABS(U) .LT. UMIN)  U = UMIN 
            IF (ABS(V) .LT. UMIN)  V = UMIN 
            V_ALPHA(I,L,1) = U 
            V_ALPHA(I,L,2) = V
            V_ALPHA(I,L,3) = - U
            V_ALPHA(I,L,4) = - V
 10       CONTINUE
 20     CONTINUE
      END IF
C
C  Case with 8 azimuths.
C
      IF (NAZ.EQ.8)  THEN
        DO 40 L = LEV1,LEV2
          DO 30 I = IL1,IL2
            U = VEL_U(I,L)
            V = VEL_V(I,L)
            IF (ABS(U) .LT. UMIN)  U = UMIN  
            IF (ABS(V) .LT. UMIN)  V = UMIN  
            V_ALPHA(I,L,1) = U 
            V_ALPHA(I,L,2) = COS45 * ( V + U )
            V_ALPHA(I,L,3) = V
            V_ALPHA(I,L,4) = COS45 * ( V - U )
            V_ALPHA(I,L,5) = - U
            V_ALPHA(I,L,6) = - V_ALPHA(I,L,2)
            V_ALPHA(I,L,7) = - V
            V_ALPHA(I,L,8) = - V_ALPHA(I,L,4)
 30       CONTINUE
 40     CONTINUE
      END IF
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_FLUX (FLUX_U,FLUX_V,DRAG_U,DRAG_V,ALT,DENSITY,
     1                       DENSB,M_ALPHA,AK_ALPHA,K_ALPHA,SLOPE,
     2                       NAZ,IL1,IL2,LEV1,LEV2,NLONS,NLEVS,NAZMTH,
     3                       LORMS)
C
C  Calculate zonal and meridional components of the vertical flux 
C  of horizontal momentum and corresponding wave drag (force per unit mass)
C  on a longitude by altitude grid for the Hines' Doppler spread 
C  GWD parameterization scheme.
C  NOTE: only 4 or 8 azimuths can be used.
C
C  Aug. 6/95 - C. McLandress
C
C  Output arguements:
C
C     * FLUX_U = zonal component of vertical momentum flux (Pascals)
C     * FLUX_V = meridional component of vertical momentum flux (Pascals)
C     * DRAG_U = zonal component of drag (m/s^2).
C     * DRAG_V = meridional component of drag (m/s^2).
C
C  Input arguements:
C
C     * ALT       = altitudes (m).
C     * DENSITY   = background density (kg/m^3).
C     * DENSB     = background density at bottom level (kg/m^3).
C     * M_ALPHA   = cutoff vertical wavenumber (1/m).
C     * AK_ALPHA  = spectral amplitude factor (i.e., {AjKj} in m^4/s^2).
C     * K_ALPHA   = horizontal wavenumber (1/m).
C     * SLOPE     = slope of incident vertical wavenumber spectrum.
C     * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
C     * IL1       = first longitudinal index to use (IL1 >= 1).
C     * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1      = first altitude level to use (LEV1 >=1). 
C     * LEV2      = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS     = number of longitudes.
C     * NLEVS     = number of vertical levels.
C     * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).
C
C     * LORMS       = .TRUE. for drag computation 
C
C  Constant in DATA statement.
C
C     * COS45 = cosine of 45 degrees. 		
C
C  Subroutine arguements.
C
      INTEGER  NAZ, IL1, IL2, LEV1, LEV2
      INTEGER  NLONS, NLEVS, NAZMTH
      REAL  SLOPE
      REAL  FLUX_U(NLONS,NLEVS), FLUX_V(NLONS,NLEVS)
      REAL  DRAG_U(NLONS,NLEVS), DRAG_V(NLONS,NLEVS)
      REAL  ALT(NLONS,NLEVS),    DENSITY(NLONS,NLEVS), DENSB(NLONS)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  AK_ALPHA(NLONS,NAZMTH), K_ALPHA(NLONS,NAZMTH)
C
      LOGICAL LORMS(NLONS)
C
C  Internal variables.
C
      INTEGER  I, L, LEV1P, LEV2M
      REAL  COS45, PROD2, PROD4, PROD6, PROD8, DENDZ, DENDZ2
      DATA  COS45 / 0.7071068 /   
C-----------------------------------------------------------------------
C
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
      LEV2P = LEV2 + 1
C
C  Sum over azimuths for case where SLOPE = 1.
C
      IF (SLOPE.EQ.1.)  THEN
C
C  Case with 4 azimuths.
C
        IF (NAZ.EQ.4)  THEN
          DO 20 L = LEV1,LEV2
            DO 10 I = IL1,IL2
              FLUX_U(I,L) = AK_ALPHA(I,1)*K_ALPHA(I,1)*M_ALPHA(I,L,1)
     ^                    - AK_ALPHA(I,3)*K_ALPHA(I,3)*M_ALPHA(I,L,3)
              FLUX_V(I,L) = AK_ALPHA(I,2)*K_ALPHA(I,2)*M_ALPHA(I,L,2)
     ^                    - AK_ALPHA(I,4)*K_ALPHA(I,4)*M_ALPHA(I,L,4)
 10         CONTINUE
 20       CONTINUE
        END IF
C
C  Case with 8 azimuths.
C
        IF (NAZ.EQ.8)  THEN
          DO 40 L = LEV1,LEV2
            DO 30 I = IL1,IL2
              PROD2 = AK_ALPHA(I,2)*K_ALPHA(I,2)*M_ALPHA(I,L,2)
              PROD4 = AK_ALPHA(I,4)*K_ALPHA(I,4)*M_ALPHA(I,L,4)
              PROD6 = AK_ALPHA(I,6)*K_ALPHA(I,6)*M_ALPHA(I,L,6)
              PROD8 = AK_ALPHA(I,8)*K_ALPHA(I,8)*M_ALPHA(I,L,8)
              FLUX_U(I,L) = 
     ^                AK_ALPHA(I,1)*K_ALPHA(I,1)*M_ALPHA(I,L,1)
     ^              - AK_ALPHA(I,5)*K_ALPHA(I,5)*M_ALPHA(I,L,5)
     ^              + COS45 * ( PROD2 - PROD4 - PROD6 + PROD8 )
              FLUX_V(I,L) = 
     ^                AK_ALPHA(I,3)*K_ALPHA(I,3)*M_ALPHA(I,L,3)
     ^              - AK_ALPHA(I,7)*K_ALPHA(I,7)*M_ALPHA(I,L,7)
     ^              + COS45 * ( PROD2 + PROD4 - PROD6 - PROD8 )
 30         CONTINUE
 40       CONTINUE
        END IF
C
      END IF
C
C  Sum over azimuths for case where SLOPE not equal to 1.
C
      IF (SLOPE.NE.1.)  THEN
C
C  Case with 4 azimuths.
C
        IF (NAZ.EQ.4)  THEN
          DO 60 L = LEV1,LEV2
            DO 50 I = IL1,IL2
              FLUX_U(I,L) = 
     ^               AK_ALPHA(I,1)*K_ALPHA(I,1)*M_ALPHA(I,L,1)**SLOPE
     ^             - AK_ALPHA(I,3)*K_ALPHA(I,3)*M_ALPHA(I,L,3)**SLOPE
              FLUX_V(I,L) = 
     ^               AK_ALPHA(I,2)*K_ALPHA(I,2)*M_ALPHA(I,L,2)**SLOPE
     ^             - AK_ALPHA(I,4)*K_ALPHA(I,4)*M_ALPHA(I,L,4)**SLOPE
 50         CONTINUE
 60       CONTINUE
        END IF
C
C  Case with 8 azimuths.
C
        IF (NAZ.EQ.8)  THEN
          DO 80 L = LEV1,LEV2
            DO 70 I = IL1,IL2
              PROD2 = AK_ALPHA(I,2)*K_ALPHA(I,2)*M_ALPHA(I,L,2)**SLOPE
              PROD4 = AK_ALPHA(I,4)*K_ALPHA(I,4)*M_ALPHA(I,L,4)**SLOPE
              PROD6 = AK_ALPHA(I,6)*K_ALPHA(I,6)*M_ALPHA(I,L,6)**SLOPE
              PROD8 = AK_ALPHA(I,8)*K_ALPHA(I,8)*M_ALPHA(I,L,8)**SLOPE
              FLUX_U(I,L) = 
     ^                AK_ALPHA(I,1)*K_ALPHA(I,1)*M_ALPHA(I,L,1)**SLOPE
     ^              - AK_ALPHA(I,5)*K_ALPHA(I,5)*M_ALPHA(I,L,5)**SLOPE
     ^              + COS45 * ( PROD2 - PROD4 - PROD6 + PROD8 )
              FLUX_V(I,L) = 
     ^                AK_ALPHA(I,3)*K_ALPHA(I,3)*M_ALPHA(I,L,3)**SLOPE
     ^              - AK_ALPHA(I,7)*K_ALPHA(I,7)*M_ALPHA(I,L,7)**SLOPE
     ^              + COS45 * ( PROD2 + PROD4 - PROD6 - PROD8 )
 70         CONTINUE
 80       CONTINUE
        END IF
C
      END IF
C
C  Calculate flux from sum.
C
      DO 100 L = LEV1,LEV2
        DO 90 I = IL1,IL2
          FLUX_U(I,L) = FLUX_U(I,L) * DENSB(I) / SLOPE
          FLUX_V(I,L) = FLUX_V(I,L) * DENSB(I) / SLOPE
  90    CONTINUE
 100  CONTINUE
C
C  Calculate drag at intermediate levels using centered differences 
C      
      DO 120 L = LEV1P,LEV2M
        DO 110 I = IL1,IL2
        IF (LORMS(I)) THEN
ccc       DENDZ2 = DENSITY(I,L) * ( ALT(I,L+1) - ALT(I,L-1) )
          DENDZ2 = DENSITY(I,L) * ( ALT(I,L-1) - ALT(I,L) ) 
ccc       DRAG_U(I,L) = - ( FLUX_U(I,L+1) - FLUX_U(I,L-1) ) / DENDZ2
          DRAG_U(I,L) = - ( FLUX_U(I,L-1) - FLUX_U(I,L) ) / DENDZ2
ccc       DRAG_V(I,L) = - ( FLUX_V(I,L+1) - FLUX_V(I,L-1) ) / DENDZ2
          DRAG_V(I,L) = - ( FLUX_V(I,L-1) - FLUX_V(I,L) ) / DENDZ2
          
        ENDIF
 110    CONTINUE
 120  CONTINUE
C
C  Drag at first and last levels using one-side differences.
C 
      DO 130 I = IL1,IL2
      IF (LORMS(I)) THEN
        DENDZ = DENSITY(I,LEV1) * ( ALT(I,LEV1) - ALT(I,LEV1P) ) 
        DRAG_U(I,LEV1) =  FLUX_U(I,LEV1)  / DENDZ
        DRAG_V(I,LEV1) =  FLUX_V(I,LEV1)  / DENDZ
      ENDIF
 130  CONTINUE
      DO 140 I = IL1,IL2
      IF (LORMS(I)) THEN
        DENDZ = DENSITY(I,LEV2) * ( ALT(I,LEV2M) - ALT(I,LEV2) )
        DRAG_U(I,LEV2) = - ( FLUX_U(I,LEV2M) - FLUX_U(I,LEV2) ) / DENDZ
        DRAG_V(I,LEV2) = - ( FLUX_V(I,LEV2M) - FLUX_V(I,LEV2) ) / DENDZ
      ENDIF
 140  CONTINUE
      IF (NLEVS .GT. LEV2) THEN
      DO 150 I = IL1,IL2
      IF (LORMS(I)) THEN
        DENDZ = DENSITY(I,LEV2P) * ( ALT(I,LEV2) - ALT(I,LEV2P) )
        DRAG_U(I,LEV2P) = -  FLUX_U(I,LEV2)  / DENDZ
        DRAG_V(I,LEV2P) = - FLUX_V(I,LEV2)  / DENDZ
      ENDIF
150   CONTINUE
      ENDIF
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_HEAT (HEAT,DIFFCO,M_ALPHA,DMDZ_ALPHA,
     1                       AK_ALPHA,K_ALPHA,BVFREQ,DENSITY,DENSB,
     2                       SIGMA_T,VISC_MOL,KSTAR,SLOPE,F2,F3,F5,F6,
     3                       NAZ,IL1,IL2,LEV1,LEV2,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the gravity wave induced heating and 
C  diffusion coefficient on a longitude by altitude grid for  
C  the Hines' Doppler spread gravity wave drag parameterization scheme.
C
C  Aug. 6/95 - C. McLandress
C
C  Output arguements:
C
C     * HEAT   = gravity wave heating (K/sec).
C     * DIFFCO = diffusion coefficient (m^2/sec)
C
C  Input arguements:
C
C     * M_ALPHA     = cutoff vertical wavenumber (1/m).
C     * DMDZ_ALPHA  = vertical derivative of cutoff wavenumber.
C     * AK_ALPHA    = spectral amplitude factor of each azimuth 
C                     (i.e., {AjKj} in m^4/s^2).
C     * K_ALPHA     = horizontal wavenumber of each azimuth (1/m).
C     * BVFREQ      = background Brunt Vassala frequency (rad/sec).
C     * DENSITY     = background density (kg/m^3).
C     * DENSB       = background density at bottom level (kg/m^3).
C     * SIGMA_T     = total rms horizontal wind (m/s).
C     * VISC_MOL    = molecular viscosity (m^2/s).
C     * KSTAR       = typical gravity wave horizontal wavenumber (1/m).
C     * SLOPE       = slope of incident vertical wavenumber spectrum.
C     * F2,F3,F5,F6 = Hines's fudge factors.
C     * NAZ         = actual number of horizontal azimuths used.
C     * IL1         = first longitudinal index to use (IL1 >= 1).
C     * IL2         = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1        = first altitude level to use (LEV1 >=1). 
C     * LEV2        = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS       = number of longitudes.
C     * NLEVS       = number of vertical levels.
C     * NAZMTH      = azimuthal array dimension (NAZMTH >= NAZ).
C
      INTEGER  NAZ, IL1, IL2, LEV1, LEV2, NLONS, NLEVS, NAZMTH
      REAL  KSTAR(NLONS), SLOPE, F2, F3, F5, F6
      REAL  HEAT(NLONS,NLEVS), DIFFCO(NLONS,NLEVS)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH), DMDZ_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  AK_ALPHA(NLONS,NAZMTH), K_ALPHA(NLONS,NAZMTH)
      REAL  BVFREQ(NLONS,NLEVS), DENSITY(NLONS,NLEVS),  DENSB(NLONS) 
      REAL  SIGMA_T(NLONS,NLEVS), VISC_MOL(NLONS,NLEVS)
C
C Internal variables.
C
      INTEGER  I, L, N
      REAL  M_SUB_M_TURB, M_SUB_M_MOL, M_SUB_M, HEATNG
      REAL  VISC, VISC_MIN, CPGAS, SM1
C
C specific heat at constant pressure
C
      DATA  CPGAS / 1004. / 
C             
C minimum permissible viscosity
C
      DATA  VISC_MIN / 1.E-10 /       
C-----------------------------------------------------------------------     
C
C  Initialize heating array.
C
      DO 20 L = 1,NLEVS
        DO 10 I = 1,NLONS
          HEAT(I,L) = 0.
  10    CONTINUE
  20  CONTINUE
C
C  Perform sum over azimuths for case where SLOPE = 1.
C
      IF (SLOPE.EQ.1.)  THEN
        DO 50 N = 1,NAZ
          DO 40 L = LEV1,LEV2
            DO 30 I = IL1,IL2
              HEAT(I,L) = HEAT(I,L) + AK_ALPHA(I,N) * K_ALPHA(I,N) 
     ^                    * DMDZ_ALPHA(I,L,N) 
 30         CONTINUE
 40       CONTINUE
 50     CONTINUE
      END IF
C
C  Perform sum over azimuths for case where SLOPE not 1.
C
      IF (SLOPE.NE.1.)  THEN
        SM1 = SLOPE - 1.
        DO 80 N = 1,NAZ
          DO 70 L = LEV1,LEV2
            DO 60 I = IL1,IL2
              HEAT(I,L) = HEAT(I,L) + AK_ALPHA(I,N) * K_ALPHA(I,N) 
     ^                    * M_ALPHA(I,L,N)**SM1 * DMDZ_ALPHA(I,L,N) 
 60         CONTINUE
 70       CONTINUE
 80     CONTINUE
      END IF
C
C  Heating and diffusion.
C
      DO 100 L = LEV1,LEV2
        DO 90 I = IL1,IL2
C
C  Maximum permissible value of cutoff wavenumber is the smaller 
C  of the instability-induced wavenumber (M_SUB_M_TURB) and 
C  that imposed by molecular viscosity (M_SUB_M_MOL).
C
          VISC    = AMAX1 ( VISC_MOL(I,L), VISC_MIN )
          M_SUB_M_TURB = BVFREQ(I,L) / ( F2 * SIGMA_T(I,L) )
          M_SUB_M_MOL  = (BVFREQ(I,L)*KSTAR(I)/VISC)**0.33333333/F3
          M_SUB_M      = AMIN1 ( M_SUB_M_TURB, M_SUB_M_MOL )
C
          HEATNG = - HEAT(I,L) * F5 * BVFREQ(I,L) / M_SUB_M 
     ^               * DENSB(I) / DENSITY(I,L) 
          DIFFCO(I,L) = F6 * HEATNG**0.33333333 / M_SUB_M**1.33333333
          HEAT(I,L)   = HEATNG / CPGAS
C
 90     CONTINUE
 100  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_SIGMA (SIGMA_T,SIGMA_ALPHA,SIGSQH_ALPHA,
     1                        NAZ,LEV,IL1,IL2,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the total rms and azimuthal rms horizontal 
C  velocities at a given level on a longitude by altitude grid for 
C  the Hines' Doppler spread GWD parameterization scheme.
C  NOTE: only four or eight azimuths can be used.
C
C  Aug. 7/95 - C. McLandress
C
C  Output arguements:
C
C     * SIGMA_T      = total rms horizontal wind (m/s).
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C
C  Input arguements:
C
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * NAZ       = actual number of horizontal azimuths used (must be 4 or 8).
C     * LEV       = altitude level to process.
C     * IL1       = first longitudinal index to use (IL1 >= 1).
C     * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS     = number of longitudes.
C     * NLEVS     = number of vertical levels.
C     * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Subroutine arguements.
C
      INTEGER  LEV, NAZ, IL1, IL2
      INTEGER  NLONS, NLEVS, NAZMTH
      REAL  SIGMA_T(NLONS,NLEVS)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
C
C  Internal variables.
C
      INTEGER  I, N
      REAL  SUM_EVEN, SUM_ODD 
C-----------------------------------------------------------------------     
C
C  Calculate azimuthal rms velocity for the 4 azimuth case.
C
      IF (NAZ.EQ.4)  THEN
        DO 10 I = IL1,IL2
          SIGMA_ALPHA(I,LEV,1) = SQRT ( SIGSQH_ALPHA(I,LEV,1)
     ^                                + SIGSQH_ALPHA(I,LEV,3) )
          SIGMA_ALPHA(I,LEV,2) = SQRT ( SIGSQH_ALPHA(I,LEV,2)
     ^                                + SIGSQH_ALPHA(I,LEV,4) )
          SIGMA_ALPHA(I,LEV,3) = SIGMA_ALPHA(I,LEV,1)
          SIGMA_ALPHA(I,LEV,4) = SIGMA_ALPHA(I,LEV,2)
 10     CONTINUE
      END IF
C
C  Calculate azimuthal rms velocity for the 8 azimuth case.
C
      IF (NAZ.EQ.8)  THEN
        DO 20 I = IL1,IL2
          SUM_ODD  = ( SIGSQH_ALPHA(I,LEV,1) 
     ^                 + SIGSQH_ALPHA(I,LEV,3) 
     ^                 + SIGSQH_ALPHA(I,LEV,5) 
     ^                 + SIGSQH_ALPHA(I,LEV,7) ) / 2.
          SUM_EVEN = ( SIGSQH_ALPHA(I,LEV,2) 
     ^                 + SIGSQH_ALPHA(I,LEV,4)
     ^                 + SIGSQH_ALPHA(I,LEV,6) 
     ^                 + SIGSQH_ALPHA(I,LEV,8) ) / 2.
          SIGMA_ALPHA(I,LEV,1) = SQRT ( SIGSQH_ALPHA(I,LEV,1) 
     ^                           + SIGSQH_ALPHA(I,LEV,5) + SUM_EVEN )
          SIGMA_ALPHA(I,LEV,2) = SQRT ( SIGSQH_ALPHA(I,LEV,2) 
     ^                           + SIGSQH_ALPHA(I,LEV,6) + SUM_ODD )
          SIGMA_ALPHA(I,LEV,3) = SQRT ( SIGSQH_ALPHA(I,LEV,3) 
     ^                           + SIGSQH_ALPHA(I,LEV,7) + SUM_EVEN )
          SIGMA_ALPHA(I,LEV,4) = SQRT ( SIGSQH_ALPHA(I,LEV,4) 
     ^                           + SIGSQH_ALPHA(I,LEV,8) + SUM_ODD )
          SIGMA_ALPHA(I,LEV,5) = SIGMA_ALPHA(I,LEV,1)
          SIGMA_ALPHA(I,LEV,6) = SIGMA_ALPHA(I,LEV,2)
          SIGMA_ALPHA(I,LEV,7) = SIGMA_ALPHA(I,LEV,3)
          SIGMA_ALPHA(I,LEV,8) = SIGMA_ALPHA(I,LEV,4)
 20     CONTINUE
      END IF
C
C  Calculate total rms velocity.
C
      DO 50 I = IL1,IL2
        SIGMA_T(I,LEV) = 0.
 50   CONTINUE
      DO 70 N = 1,NAZ
        DO 60 I = IL1,IL2
          SIGMA_T(I,LEV) = SIGMA_T(I,LEV) + SIGSQH_ALPHA(I,LEV,N)
 60     CONTINUE
 70   CONTINUE
      DO 80 I = IL1,IL2
        SIGMA_T(I,LEV) = SQRT ( SIGMA_T(I,LEV) )
 80   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------     
      END

      SUBROUTINE HINES_INTGRL (I_ALPHA,V_ALPHA,M_ALPHA,BVFB,SLOPE,
     1                         NAZ,LEV,IL1,IL2,NLONS,NLEVS,NAZMTH,
     2                         LORMS)
C
C  This routine calculates the vertical wavenumber integral
C  for a single vertical level at each azimuth on a longitude grid
C  for the Hines' Doppler spread GWD parameterization scheme.
C  NOTE: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
C        (2) the integral is written in terms of the product QM
C            which by construction is always less than 1. Series
C            solutions are used for small |QM| and analytical solutions
C            for remaining values.
C
C  Aug. 8/95 - C. McLandress
C
C  Output arguement:
C
C     * I_ALPHA = Hines' integral.
C
C  Input arguements:
C
C     * V_ALPHA = azimuthal wind component (m/s). 
C     * M_ALPHA = azimuthal cutoff vertical wavenumber (1/m).
C     * BVFB    = background Brunt Vassala frequency at model bottom.
C     * SLOPE   = slope of initial vertical wavenumber spectrum 
C     *           (must use SLOPE = 1., 1.5 or 2.)
C     * NAZ     = actual number of horizontal azimuths used.
C     * LEV     = altitude level to process.
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical levels.
C     * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).
C
C     * LORMS       = .TRUE. for drag computation 
C
C  Constants in DATA statements:
C
C     * QMIN = minimum value of Q_ALPHA (avoids indeterminant form of integral)
C     * QM_MIN = minimum value of Q_ALPHA * M_ALPHA (used to avoid numerical
C     *          problems).
C
      INTEGER  LEV, NAZ, IL1, IL2, NLONS, NLEVS, NAZMTH
      REAL  I_ALPHA(NLONS,NAZMTH)
      REAL  V_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  BVFB(NLONS), SLOPE
C
      LOGICAL LORMS(NLONS)
C
C  Internal variables.
C
      INTEGER  I, N
      REAL  Q_ALPHA, QM, SQRTQM, Q_MIN, QM_MIN
C
      DATA  Q_MIN / 1.0 /, QM_MIN / 0.01 /
C-----------------------------------------------------------------------     
C
C  For integer value SLOPE = 1.
C
      IF (SLOPE .EQ. 1.)  THEN
C
        DO 20 N = 1,NAZ
          DO 10 I = IL1,IL2
          IF (LORMS(I)) THEN
C
            Q_ALPHA = V_ALPHA(I,LEV,N) / BVFB(I)
            QM = Q_ALPHA * M_ALPHA(I,LEV,N)
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. 0.)  THEN
                I_ALPHA(I,N) = M_ALPHA(I,LEV,N)**2 / 2.
              ELSE
                I_ALPHA(I,N) = ( QM**2/2. + QM**3/3. + QM**4/4.
     ^                           + QM**5/5. ) / Q_ALPHA**2
              END IF
            ELSE
              I_ALPHA(I,N) = - ( ALOG(1.-QM) + QM ) / Q_ALPHA**2
            END IF
C
          ENDIF
 10       CONTINUE
 20     CONTINUE
C
      END IF
C
C  For integer value SLOPE = 2.
C
      IF (SLOPE .EQ. 2.)  THEN
C
        DO 40 N = 1,NAZ
          DO 30 I = IL1,IL2
          IF (LORMS(I)) THEN
C
            Q_ALPHA = V_ALPHA(I,LEV,N) / BVFB(I)
            QM = Q_ALPHA * M_ALPHA(I,LEV,N)
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. 0.)  THEN
                I_ALPHA(I,N) = M_ALPHA(I,LEV,N)**3 / 3.
              ELSE
                I_ALPHA(I,N) = ( QM**3/3. + QM**4/4. + QM**5/5. 
     ^                           + QM**6/6. ) / Q_ALPHA**3
              END IF
            ELSE
              I_ALPHA(I,N) = - ( ALOG(1.-QM) + QM + QM**2/2.) 
     ^                         / Q_ALPHA**3
            END IF
C
          ENDIF
 30       CONTINUE
 40     CONTINUE
C
      END IF
C
C  For real value SLOPE = 1.5
C
      IF (SLOPE .EQ. 1.5)  THEN
C
        DO 60 N = 1,NAZ
          DO 50 I = IL1,IL2
          IF (LORMS(I)) THEN
C
            Q_ALPHA = V_ALPHA(I,LEV,N) / BVFB(I)
            QM = Q_ALPHA * M_ALPHA(I,LEV,N)       
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. 0.)  THEN
                I_ALPHA(I,N) = M_ALPHA(I,LEV,N)**2.5 / 2.5
              ELSE
                I_ALPHA(I,N) = ( QM/2.5 + QM**2/3.5 
     ^                           + QM**3/4.5 + QM**4/5.5 ) 
     ^                           * M_ALPHA(I,LEV,N)**1.5 / Q_ALPHA
              END IF
            ELSE
              QM     = ABS(QM)
              SQRTQM = SQRT(QM)
              IF (Q_ALPHA .GE. 0.)  THEN
                I_ALPHA(I,N) = ( ALOG( (1.+SQRTQM)/(1.-SQRTQM) )
     ^                          -2.*SQRTQM*(1.+QM/3.) ) / Q_ALPHA**2.5
              ELSE
                I_ALPHA(I,N) = 2. * ( ATAN(SQRTQM) + SQRTQM*(QM/3.-1.) )
     ^                          / ABS(Q_ALPHA)**2.5
              END IF
            END IF
C
          ENDIF
 50       CONTINUE
 60     CONTINUE
C
      END IF
C
C  If integral is negative (which in principal should not happen) then
C  print a message and some info since execution will abort when calculating
C  the variances.
C
c      DO 80 N = 1,NAZ
c        DO 70 I = IL1,IL2
c          IF (I_ALPHA(I,N).LT.0.)  THEN
c            WRITE (6,*) 
c            WRITE (6,*) '******************************'
c            WRITE (6,*) 'Hines integral I_ALPHA < 0 '
c            WRITE (6,*) '  longitude I=',I
c            WRITE (6,*) '  azimuth   N=',N
c            WRITE (6,*) '  level   LEV=',LEV
c            WRITE (6,*) '  I_ALPHA =',I_ALPHA(I,N)
c            WRITE (6,*) '  V_ALPHA =',V_ALPHA(I,LEV,N)
c            WRITE (6,*) '  M_ALPHA =',M_ALPHA(I,LEV,N)
c            WRITE (6,*) '  Q_ALPHA =',V_ALPHA(I,LEV,N) / BVFB(I)
c            WRITE (6,*) '  QM      =',V_ALPHA(I,LEV,N) / BVFB(I) 
c     ^                                * M_ALPHA(I,LEV,N)
c            WRITE (6,*) '******************************'
c          END IF
c 70     CONTINUE
c 80   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_SETUP (NAZ,SLOPE,F1,F2,F3,F5,F6,KSTAR,
     1                        ICUTOFF,ALT_CUTOFF,SMCO,NSMAX,IHEATCAL,
     2                       K_ALPHA,IERROR,NMESSG,NLONS,NAZMTH,COSLAT)
C
C  This routine specifies various parameters needed for the
C  the Hines' Doppler spread gravity wave drag parameterization scheme.
C
C  Aug. 8/95 - C. McLandress
C
C  Output arguements:
C
C     * NAZ        = actual number of horizontal azimuths used
C     *              (code set up presently for only NAZ = 4 or 8).
C     * SLOPE      = slope of incident vertical wavenumber spectrum
C     *              (code set up presently for SLOPE 1., 1.5 or 2.).
C     * F1         = "fudge factor" used in calculation of trial value of
C     *              azimuthal cutoff wavenumber M_ALPHA (1.2 <= F1 <= 1.9).
C     * F2         = "fudge factor" used in calculation of maximum
C     *              permissible instabiliy-induced cutoff wavenumber 
C     *              M_SUB_M_TURB (0.1 <= F2 <= 1.4).
C     * F3         = "fudge factor" used in calculation of maximum 
C     *              permissible molecular viscosity-induced cutoff wavenumber 
C     *              M_SUB_M_MOL (0.1 <= F2 <= 1.4).
C     * F5         = "fudge factor" used in calculation of heating rate
C     *              (1 <= F5 <= 3).
C     * F6         = "fudge factor" used in calculation of turbulent 
C     *              diffusivity coefficient.
C     * KSTAR      = typical gravity wave horizontal wavenumber (1/m)
C     *              used in calculation of M_SUB_M_TURB.
C     * ICUTOFF    = 1 to exponentially damp off GWD, heating and diffusion 
C     *              arrays above ALT_CUTOFF; otherwise arrays not modified.
C     * ALT_CUTOFF = altitude in meters above which exponential decay applied.
C     * SMCO       = smoother used to smooth cutoff vertical wavenumbers
C     *              and total rms winds before calculating drag or heating.
C     *              (==> a 1:SMCO:1 stencil used; SMCO >= 1.).
C     * NSMAX      = number of times smoother applied ( >= 1),
C     *            = 0 means no smoothing performed.
C     * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
C     *            = 0 means only drag and flux calculated.
C     * K_ALPHA    = horizontal wavenumber of each azimuth (1/m) which
C     *              is set here to KSTAR.
C     * IERROR     = error flag.
C     *            = 0 no errors.
C     *            = 10 ==> NAZ > NAZMTH
C     *            = 20 ==> invalid number of azimuths (NAZ must be 4 or 8).
C     *            = 30 ==> invalid slope (SLOPE must be 1., 1.5 or 2.).
C     *            = 40 ==> invalid smoother (SMCO must be >= 1.)
C
C  Input arguements:
C
C     * NMESSG  = output unit number where messages to be printed.
C     * NLONS   = number of longitudes.
C     * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).
C
      INTEGER  NAZ, NLONS, NAZMTH, IHEATCAL, ICUTOFF
      INTEGER  NMESSG, NSMAX, IERROR
      REAL  KSTAR(NLONS), SLOPE, F1, F2, F3, F5, F6, ALT_CUTOFF, SMCO
      REAL  K_ALPHA(NLONS,NAZMTH),COSLAT(NLONS)
      REAL  KSMIN, KSMAX
C
C Internal variables.
C
      INTEGER  I, N
C-----------------------------------------------------------------------     
C
C  Specify constants.
C
      NAZ   = 8
      SLOPE = 1.
      F1    = 1.5 
      F2    = 0.3 
      F3    = 1.0 
      F5    = 3.0 
      F6    = 1.0       
      KSMIN = 1.E-5       
      KSMAX = 1.E-4       
      DO I=1,NLONS
         KSTAR(I) = KSMIN/( COSLAT(I)+(KSMIN/KSMAX) )      
      ENDDO
      ICUTOFF    = 1   
      ALT_CUTOFF = 105.E3
      SMCO       = 2.0 
c      SMCO       = 1.0 
      NSMAX      = 5
c      NSMAX      = 2
      IHEATCAL   = 0 
C
C  Print information to output file.
C
c      WRITE (NMESSG,6000)
c 6000 FORMAT (/' Subroutine HINES_SETUP:')
c      WRITE (NMESSG,*)  '  SLOPE = ', SLOPE
c      WRITE (NMESSG,*)  '  NAZ = ', NAZ
c      WRITE (NMESSG,*)  '  F1,F2,F3  = ', F1, F2, F3
c      WRITE (NMESSG,*)  '  F5,F6     = ', F5, F6
c      WRITE (NMESSG,*)  '  KSTAR     = ', KSTAR
c     >           ,'  COSLAT     = ', COSLAT
c      IF (ICUTOFF .EQ. 1)  THEN
c        WRITE (NMESSG,*) '  Drag exponentially damped above ',
c     &                       ALT_CUTOFF/1.E3
c     END IF
c      IF (NSMAX.LT.1 )  THEN
c        WRITE (NMESSG,*) '  No smoothing of cutoff wavenumbers, etc'
c      ELSE
c        WRITE (NMESSG,*) '  Cutoff wavenumbers and sig_t smoothed:'
c        WRITE (NMESSG,*) '    SMCO  =', SMCO
c        WRITE (NMESSG,*) '    NSMAX =', NSMAX
c     END IF
C
C  Check that things are setup correctly and log error if not
C
      IERROR = 0
      IF (NAZ .GT. NAZMTH)                                  IERROR = 10
      IF (NAZ.NE.4 .AND. NAZ.NE.8)                          IERROR = 20
      IF (SLOPE.NE.1. .AND. SLOPE.NE.1.5 .AND. SLOPE.NE.2.) IERROR = 30
      IF (SMCO .LT. 1.)                                     IERROR = 40
C
C  Use single value for azimuthal-dependent horizontal wavenumber.
C
      DO 20 N = 1,NAZ
        DO 10 I = 1,NLONS
          K_ALPHA(I,N) = KSTAR(I)
 10     CONTINUE
 20   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_PRINT (FLUX_U,FLUX_V,DRAG_U,DRAG_V,ALT,SIGMA_T,
     1                        SIGMA_ALPHA,V_ALPHA,M_ALPHA,
     2                        IU_PRINT,IV_PRINT,NMESSG,
     3                        ILPRT1,ILPRT2,LEVPRT1,LEVPRT2,
     4                        NAZ,NLONS,NLEVS,NAZMTH)
C
C  Print out altitude profiles of various quantities from
C  Hines' Doppler spread gravity wave drag parameterization scheme.
C  (NOTE: only for NAZ = 4 or 8). 
C
C  Aug. 8/95 - C. McLandress
C
C  Input arguements:
C
C     * IU_PRINT = 1 to print out values in east-west direction.
C     * IV_PRINT = 1 to print out values in north-south direction.
C     * NMESSG   = unit number for printed output.
C     * ILPRT1   = first longitudinal index to print.
C     * ILPRT2   = last longitudinal index to print.
C     * LEVPRT1  = first altitude level to print.
C     * LEVPRT2  = last altitude level to print.
C
      INTEGER  NAZ, ILPRT1, ILPRT2, LEVPRT1, LEVPRT2
      INTEGER  NLONS, NLEVS, NAZMTH
      INTEGER  IU_PRINT, IV_PRINT, NMESSG
      REAL  FLUX_U(NLONS,NLEVS), FLUX_V(NLONS,NLEVS)
      REAL  DRAG_U(NLONS,NLEVS), DRAG_V(NLONS,NLEVS)
      REAL  ALT(NLONS,NLEVS), SIGMA_T(NLONS,NLEVS)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  V_ALPHA(NLONS,NLEVS,NAZMTH), M_ALPHA(NLONS,NLEVS,NAZMTH)
C
C  Internal variables.
C
      INTEGER  N_EAST, N_WEST, N_NORTH, N_SOUTH
      INTEGER  I, L
C-----------------------------------------------------------------------
C
C  Azimuthal indices of cardinal directions.
C
      N_EAST = 1
      IF (NAZ.EQ.4)  THEN
        N_WEST  = 3       
        N_NORTH = 2
        N_SOUTH = 4       
      ELSE IF (NAZ.EQ.8)  THEN
        N_WEST  = 5       
        N_NORTH = 3
        N_SOUTH = 7       
      END IF
C
C  Print out values for range of longitudes.
C
      DO 100 I = ILPRT1,ILPRT2
C
C  Print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
C
        IF (IU_PRINT.EQ.1)  THEN
          WRITE (NMESSG,*) 
          WRITE (NMESSG,6001) I
          WRITE (NMESSG,6005) 
 6001     FORMAT ( 'Hines GW (east-west) at longitude I =',I3)
 6005     FORMAT (15x,' U ',2x,'sig_E',2x,'sig_T',3x,'m_E',
     &            4x,'m_W',4x,'fluxU',5x,'gwdU')
          DO 10 L = LEVPRT1,LEVPRT2
            WRITE (NMESSG,6701) ALT(I,L)/1.E3, V_ALPHA(I,L,N_EAST), 
     &                          SIGMA_ALPHA(I,L,N_EAST), SIGMA_T(I,L),
     &                          M_ALPHA(I,L,N_EAST)*1.E3, 
     &                          M_ALPHA(I,L,N_WEST)*1.E3,
     &                          FLUX_U(I,L)*1.E5, DRAG_U(I,L)*24.*3600.
  10      CONTINUE
 6701     FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
        END IF
C
C  Print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
C
        IF (IV_PRINT.EQ.1)  THEN
          WRITE(NMESSG,*) 
          WRITE(NMESSG,6002) I
 6002     FORMAT ( 'Hines GW (north-south) at longitude I =',I3)
          WRITE(NMESSG,6006) 
 6006     FORMAT (15x,' V ',2x,'sig_N',2x,'sig_T',3x,'m_N',
     &            4x,'m_S',4x,'fluxV',5x,'gwdV')
          DO 20 L = LEVPRT1,LEVPRT2
            WRITE (NMESSG,6701) ALT(I,L)/1.E3, V_ALPHA(I,L,N_NORTH), 
     &                          SIGMA_ALPHA(I,L,N_NORTH), SIGMA_T(I,L),
     &                          M_ALPHA(I,L,N_NORTH)*1.E3, 
     &                          M_ALPHA(I,L,N_SOUTH)*1.E3,
     &                          FLUX_V(I,L)*1.E5, DRAG_V(I,L)*24.*3600.
 20       CONTINUE
        END IF
C
 100  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE HINES_EXP (DATA,DATA_ZMAX,ALT,ALT_EXP,IORDER,
     1                      IL1,IL2,LEV1,LEV2,NLONS,NLEVS)
C
C  This routine exponentially damps a longitude by altitude array 
C  of data above a specified altitude.
C
C  Aug. 13/95 - C. McLandress
C
C  Output arguements:
C
C     * DATA = modified data array.
C
C  Input arguements:
C
C     * DATA    = original data array.
C     * ALT     = altitudes.
C     * ALT_EXP = altitude above which exponential decay applied.
C     * IORDER	= 1 means vertical levels are indexed from top down 
C     *           (i.e., highest level indexed 1 and lowest level NLEVS);
C     *           .NE. 1 highest level is index NLEVS.
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1    = first altitude level to use (LEV1 >=1). 
C     * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical
C
C  Input work arrays:
C
C     * DATA_ZMAX = data values just above altitude ALT_EXP.
C
      INTEGER  IORDER, IL1, IL2, LEV1, LEV2, NLONS, NLEVS
      REAL  ALT_EXP
      REAL  DATA(NLONS,NLEVS), DATA_ZMAX(NLONS), ALT(NLONS,NLEVS)
C
C Internal variables.
C
      INTEGER  LEVBOT, LEVTOP, LINCR, I, L
      REAL  HSCALE
      DATA  HSCALE / 5.E3 /
C-----------------------------------------------------------------------     
C
C  Index of lowest altitude level (bottom of drag calculation).
C
      LEVBOT = LEV2
      LEVTOP = LEV1
      LINCR  = 1
      IF (IORDER.NE.1)  THEN
        LEVBOT = LEV1
        LEVTOP = LEV2
        LINCR  = -1
      END IF
C
C  Data values at first level above ALT_EXP.
C
      DO 20 I = IL1,IL2
        DO 10 L = LEVTOP,LEVBOT,LINCR
          IF (ALT(I,L) .GE. ALT_EXP)  THEN
            DATA_ZMAX(I) = DATA(I,L) 
          END IF	   
 10     CONTINUE
 20   CONTINUE
C
C  Exponentially damp field above ALT_EXP to model top at L=1.
C
      DO 40 L = 1,LEV2 
        DO 30 I = IL1,IL2
          IF (ALT(I,L) .GE. ALT_EXP)  THEN
            DATA(I,L) = DATA_ZMAX(I) * EXP( (ALT_EXP-ALT(I,L))/HSCALE )
          END IF
 30     CONTINUE
 40   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END

      SUBROUTINE VERT_SMOOTH (DATA,WORK,COEFF,NSMOOTH,
     1                        IL1,IL2,LEV1,LEV2,NLONS,NLEVS)
C
C  Smooth a longitude by altitude array in the vertical over a
C  specified number of levels using a three point smoother. 
C
C  NOTE: input array DATA is modified on output!
C
C  Aug. 3/95 - C. McLandress
C
C  Output arguement:
C
C     * DATA    = smoothed array (on output).
C
C  Input arguements:
C
C     * DATA    = unsmoothed array of data (on input).
C     * WORK    = work array of same dimension as DATA.
C     * COEFF   = smoothing coefficient for a 1:COEFF:1 stencil.
C     *           (e.g., COEFF = 2 will result in a smoother which
C     *           weights the level L gridpoint by two and the two 
C     *           adjecent levels (L+1 and L-1) by one).
C     * NSMOOTH = number of times to smooth in vertical.
C     *           (e.g., NSMOOTH=1 means smoothed only once, 
C     *           NSMOOTH=2 means smoothing repeated twice, etc.)
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1    = first altitude level to use (LEV1 >=1). 
C     * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical levels.
C
C  Subroutine arguements.
C
      INTEGER  NSMOOTH, IL1, IL2, LEV1, LEV2, NLONS, NLEVS
      REAL  COEFF
      REAL  DATA(NLONS,NLEVS), WORK(NLONS,NLEVS)
C
C  Internal variables.
C
      INTEGER  I, L, NS, LEV1P, LEV2M
      REAL  SUM_WTS
C-----------------------------------------------------------------------     
C
C  Calculate sum of weights.
C
      SUM_WTS = COEFF + 2.
C
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
C
C  Smooth NSMOOTH times
C
      DO 50 NS = 1,NSMOOTH
C
C  Copy data into work array.
C
        DO 20 L = LEV1,LEV2
          DO 10 I = IL1,IL2
            WORK(I,L) = DATA(I,L)
 10       CONTINUE
 20     CONTINUE
C
C  Smooth array WORK in vertical direction and put into DATA.
C
        DO 40 L = LEV1P,LEV2M
          DO 30 I = IL1,IL2
            DATA(I,L) = ( WORK(I,L+1) + COEFF*WORK(I,L) + WORK(I,L-1) ) 
     &                    / SUM_WTS 
 30       CONTINUE
 40     CONTINUE
C
 50   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END




      
      
