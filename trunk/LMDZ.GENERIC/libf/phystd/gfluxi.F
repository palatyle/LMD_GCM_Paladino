      SUBROUTINE GFLUXI(NLL,TLEV,NW,DW,DTAU,TAUCUM,W0,COSBAR,UBARI,
     *                  RSF,BTOP,BSURF,FTOPUP,FMIDP,FMIDM)
      
      use radinc_h
      use radcommon_h, only: planckir
      use comcstfi_mod, only: pi
      
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
!  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
!  THE LEVELS.  THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
!  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS  
!  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
!  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
!  FROM TOP TO BOTTOM.  SEE C.P. MCKAY, TGM NOTES.
!  THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY 
!  VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
!
! NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
! TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
! WAVEN          = WAVELENGTH FOR THE COMPUTATION
! DW             = WAVENUMBER INTERVAL
! DTAU(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
! W0(NLEVEL)     = SINGLE SCATTERING ALBEDO
! COSBAR(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
! UBARI          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
! RSF            = SURFACE REFLECTANCE
! BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
! BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
! FP(NLEVEL)     = UPWARD FLUX AT LEVELS
! FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
! FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
! FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
!-----------------------------------------------------------------------
      
      INTEGER NLL, NLAYER, L, NW, NT, NT2
      REAL*8  TERM, CPMID, CMMID
      REAL*8  PLANCK
      REAL*8  EM,EP
      REAL*8  COSBAR(L_NLAYRAD), W0(L_NLAYRAD), DTAU(L_NLAYRAD)
      REAL*8  TAUCUM(L_LEVELS), DTAUK
      REAL*8  TLEV(L_LEVELS)
      REAL*8  WAVEN, DW, UBARI, RSF
      REAL*8  BTOP, BSURF, FMIDP(L_NLAYRAD), FMIDM(L_NLAYRAD)
      REAL*8  B0(L_NLAYRAD)
      REAL*8  B1(L_NLAYRAD)
      REAL*8  ALPHA(L_NLAYRAD)
      REAL*8  LAMDA(L_NLAYRAD),XK1(L_NLAYRAD),XK2(L_NLAYRAD)
      REAL*8  GAMA(L_NLAYRAD),CP(L_NLAYRAD),CM(L_NLAYRAD)
      REAL*8  CPM1(L_NLAYRAD),CMM1(L_NLAYRAD),E1(L_NLAYRAD)
      REAL*8  E2(L_NLAYRAD)
      REAL*8  E3(L_NLAYRAD)
      REAL*8  E4(L_NLAYRAD)
      REAL*8  FTOPUP, FLUXUP, FLUXDN
      REAL*8 :: TAUMAX = L_TAUMAX

! AB : variables for interpolation
      REAL*8 C1
      REAL*8 C2
      REAL*8 P1
      REAL*8 P2
      
!=======================================================================
!     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
      
      NLAYER = L_NLAYRAD

      DO L=1,L_NLAYRAD-1

!-----------------------------------------------------------------------
! There is a problem when W0 = 1
!         open(888,file='W0')
!           if ((W0(L).eq.0.).or.(W0(L).eq.1.)) then
!             write(888,*) W0(L), L, 'gfluxi'
!           endif
! Prevent this with an if statement:
!-----------------------------------------------------------------------
         if (W0(L).eq.1.D0) then
            W0(L) = 0.99999D0
         endif
         
         ALPHA(L) = SQRT( (1.0D0-W0(L))/(1.0D0-W0(L)*COSBAR(L)) )
         LAMDA(L) = ALPHA(L)*(1.0D0-W0(L)*COSBAR(L))/UBARI
         
         NT    = int(TLEV(2*L)*NTfac)   - NTstart+1
         NT2   = int(TLEV(2*L+2)*NTfac) - NTstart+1
         
! AB : PLANCKIR(NW,NT) is replaced by P1, the linear interpolation result for a temperature NT
! AB : idem for PLANCKIR(NW,NT2) and P2
         C1 = TLEV(2*L) * NTfac - int(TLEV(2*L) * NTfac)
         C2 = TLEV(2*L+2)*NTfac - int(TLEV(2*L+2)*NTfac)
         P1 = (1.0D0 - C1) * PLANCKIR(NW,NT) + C1 * PLANCKIR(NW,NT+1)
         P2 = (1.0D0 - C2) * PLANCKIR(NW,NT2) + C2 * PLANCKIR(NW,NT2+1)
         B1(L) = (P2 - P1) / DTAU(L)
         B0(L) = P1
      END DO
      
!     Take care of special lower layer
      
      L        = L_NLAYRAD

      if (W0(L).eq.1.) then
          W0(L) = 0.99999D0
      end if
      
      ALPHA(L) = SQRT( (1.0D0-W0(L))/(1.0D0-W0(L)*COSBAR(L)) )
      LAMDA(L) = ALPHA(L)*(1.0D0-W0(L)*COSBAR(L))/UBARI
      
      ! Tsurf is used for 1st layer source function
      ! -- same results for most thin atmospheres
      ! -- and stabilizes integrations
      NT    = int(TLEV(2*L+1)*NTfac) - NTstart+1
      !! For deep, opaque, thick first layers (e.g. Saturn)
      !! what is below works much better, not unstable, ...
      !! ... and actually fully accurate because 1st layer temp (JL) 
      !NT    = int(TLEV(2*L)*NTfac) - NTstart+1
      !! (or this one yields same results
      !NT    = int( (TLEV(2*L)+TLEV(2*L+1))*0.5*NTfac ) - NTstart+1
      
      NT2   = int(TLEV(2*L)*NTfac)   - NTstart+1
      
! AB : PLANCKIR(NW,NT) is replaced by P1, the linear interpolation result for a temperature NT
! AB : idem for PLANCKIR(NW,NT2) and P2
      C1 = TLEV(2*L+1)*NTfac - int(TLEV(2*L+1)*NTfac)
      C2 = TLEV(2*L) * NTfac - int(TLEV(2*L) * NTfac)
      P1 = (1.0D0 - C1) * PLANCKIR(NW,NT) + C1 * PLANCKIR(NW,NT+1)
      P2 = (1.0D0 - C2) * PLANCKIR(NW,NT2) + C2 * PLANCKIR(NW,NT2+1)
      B1(L) = (P1 - P2) / DTAU(L)
      B0(L) = P2
      
      DO L=1,L_NLAYRAD
         GAMA(L) = (1.0D0-ALPHA(L))/(1.0D0+ALPHA(L))
         TERM    = UBARI/(1.0D0-W0(L)*COSBAR(L))
         
! CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
! AT THE TOP OF THE LAYER, THAT IS ZERO OPTICAL DEPTH
         
         CPM1(L) = B0(L)+B1(L)*TERM
         CMM1(L) = B0(L)-B1(L)*TERM
         
! CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
! BOTTOM OF THE LAYER.  THAT IS AT DTAU OPTICAL DEPTH.
! JL18 put CP and CM after the calculation of CPM1 and CMM1 to avoid unecessary calculations. 
         
         CP(L) = CPM1(L) +B1(L)*DTAU(L) 
         CM(L) = CMM1(L) +B1(L)*DTAU(L) 
      END DO
      
! NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
! FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
! WARNING IF DTAU(J) IS GREATER THAN ABOUT 35 (VAX)
! WE CLIP IT TO AVOID OVERFLOW.
      
      DO L=1,L_NLAYRAD
        EP    = EXP( MIN((LAMDA(L)*DTAU(L)),TAUMAX)) ! CLIPPED EXPONENTIAL
        EM    = 1.0D0/EP
        E1(L) = EP+GAMA(L)*EM
        E2(L) = EP-GAMA(L)*EM
        E3(L) = GAMA(L)*EP+EM
        E4(L) = GAMA(L)*EP-EM
      END DO
      
!      B81=BTOP  ! RENAME BEFORE CALLING DSOLVER - used to be to set
!      B82=BSURF ! them to real*8 - but now everything is real*8
!      R81=RSF   ! so this may not be necessary

! DOUBLE PRECISION TRIDIAGONAL SOLVER
      
      CALL DSOLVER(NLAYER,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
     *             BSURF,RSF,XK1,XK2)
      
! NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.
      
      DO L=1,L_NLAYRAD-1
         DTAUK = TAUCUM(2*L+1)-TAUCUM(2*L)
         EP    = EXP(MIN(LAMDA(L)*DTAUK,TAUMAX)) ! CLIPPED EXPONENTIAL 
         EM    = 1.0D0/EP
         TERM  = UBARI/(1.D0-W0(L)*COSBAR(L))
         
! CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
! BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH
         
         CPMID    = B0(L)+B1(L)*DTAUK +B1(L)*TERM
         CMMID    = B0(L)+B1(L)*DTAUK -B1(L)*TERM
         FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
         FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
         
! FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT
         
         FMIDP(L) = FMIDP(L)*PI
         FMIDM(L) = FMIDM(L)*PI
      END DO
      
! And now, for the special bottom layer

      L    = L_NLAYRAD

      EP   = EXP(MIN((LAMDA(L)*DTAU(L)),TAUMAX)) ! CLIPPED EXPONENTIAL 
      EM   = 1.0D0/EP
      TERM = UBARI/(1.D0-W0(L)*COSBAR(L))

! CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
! BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

      CPMID    = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
      CMMID    = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM
      FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
      FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
 
! FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT
      
      FMIDP(L) = FMIDP(L)*PI
      FMIDM(L) = FMIDM(L)*PI
      
! FLUX AT THE PTOP LEVEL
      
      EP   = 1.0D0
      EM   = 1.0D0
      TERM = UBARI/(1.0D0-W0(1)*COSBAR(1))
      
! CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
! BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH
      
      CPMID  = B0(1)+B1(1)*TERM
      CMMID  = B0(1)-B1(1)*TERM
      
      FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
      FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID
      
! FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT
      
      FTOPUP = (FLUXUP-FLUXDN)*PI
      
      
      RETURN
      END
