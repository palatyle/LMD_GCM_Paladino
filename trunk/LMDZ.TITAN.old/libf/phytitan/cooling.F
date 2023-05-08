      SUBROUTINE COOLING(NG,NL,PRESS,TEMP,Z,Q0,zlwup,zlwdn,pfluxi,icld)

c=======================================================================
c
c   Author :  C. P. Mc Kay            01/02/91
c   ------
c
c   Object :
c   --------
c
C THIS SUBROUTINE RETURNS THE COOLING RATE IN TITAN'S ATMOSPHERE
C INPUTS ARE PRESS(BARS), TEMP(K), Z(KM)
C OUTPUT IS: Q(K/SEC)C
C
C COOLING RATE COMPUTED NEGLECTING SCATTERING.
C THE TRICK OF THIS ROUTINE IS THAT IT READS IN THE OPACITIES
C FOR EACH LAYER AT EACH WAVENUMBER IN THE SPECTRAL DOMAIN
C THESE OPACITIES ARE HELD CONSTANT WITH TEMPERATURE AND TIME.
c
c   Interface:
c   ----------
c
c   Arguments:
c   ----------
c
c      input:
c      ------
c
c      nl               number of levels
c      press(nl)        pressure levels (layers)
c      temp(nl)         temperature (layers)
c      z(nl)            altitude   (m, levels)
c
c      output:
c      -------
c
c      q0(nl-1)         radiative cooling in K/sec
c      zlwup(nl)         up fluxes,   (+) upward
c      zlwdn(nl)         down fluxes, (+) downward
c      pfluxi          IR descendant a la surface (+ vers le bas)
c
c   Commons:
c   --------
c
c     COMMON/IRTAUS/dtaui(nlayer,nspeci)
c     infrared opacities of the differents layers for differents
c     spectral ranges. This common is initialized by radtitan.
c   
c     COMMON /PLANT/ CSUBP,F0PI
c     This common is initialized by tgmdat.
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   ------------

      use dimphy
      use tgmdat_mod, only: CSUBP,F0PI
      IMPLICIT NONE
#include "dimensions.h"
#include "YOMCST.h"
#include "clesphys.h"
      INTEGER NLAYER,NSPECI,NSPC1I
      PARAMETER(NLAYER=llm)
      PARAMETER (NSPECI=46,NSPC1I=47)

c  ASTUCE POUR EVITER klon... EN ATTENDANT MIEUX
      INTEGER   ngrid
      PARAMETER (ngrid=(jjm-1)*iim+2)  ! = klon
c
c   Arguments:
c   ----------

      INTEGER NG,NL,icld
      REAL PRESS(NG,NL),TEMP(NG,NL)
      REAL Z(NG,NL),Q0(NG,NL-1)
      REAL zlwup(NG,NL),zlwdn(NG,NL),UBARI2
      real pfluxi(NG)


c   Common:
c   -------

C DTAU IS PASSED EN-MASS, SO ITS DEMENSIONS ARE CRITICAL
      REAL dtaui(ngrid,NLAYER,NSPECI)
      REAL dtauip(ngrid,NLAYER,NSPECI)
      COMMON /IRTAUS/ dtaui,dtauip

c   Local:
c   ------

      REAL WNOI(NSPECI),DWNI(NSPECI)   ! SPECTAL INTERVALS
      REAL B0(ngrid,llm+1)
      REAL EM(ngrid,llm+1)
      REAL DW,WAVEN,TJ,BSURF,QOUT,QIN,eff_g,COLDEN

      INTEGER ig,K,J,I,L

c     EXTERNAL PLNCK
      REAL PLNCK,zz1,zz2,zz3,zz4,WAVNUM,Xtest

      REAL FNETIS(ngrid,llm+1),FNETI(ngrid,llm+1)
      REAL FDIS(ngrid,llm+1,nspeci),FUPIS(ngrid,llm+1,nspeci)
      REAL FDI(ngrid,llm+1), FUPI(ngrid,llm+1)

c   Data:
c   -----

      REAL RHOP,UBARI
      DATA RHOP/1.E4/      ! CONVERSION FROM PRESSURE TO MASS
      DATA UBARI/0.5/      ! MEAN COSINE FOR 2-STREAM
      DATA WNOI/
     &    11.500,    20.000,    31.250,    50.000,    75.000,
     &   100.000,   125.000,   150.000,   175.000,   200.000,
     &   225.000,   250.000,   275.000,   300.000,   325.000,
     &   350.000,   375.000,   400.000,   425.000,   450.000,
     &   475.000,   500.000,   525.000,   550.000,   575.000,
     &   600.000,   628.750,   662.838,   681.757,   683.919,
     &   686.541,   689.623,   692.704,   695.786,   715.141,
     &   733.836,   735.597,   737.358,   739.119,   742.720,
     &   748.160,   753.600,   834.560,   917.333,   926.400,
     &   935.466/
      DATA DWNI/
     &     7.000,    10.000,    12.500,    25.000,    25.000,
     &    25.000,    25.000,    25.000,    25.000,    25.000,
     &    25.000,    25.000,    25.000,    25.000,    25.000,
     &    25.000,    25.000,    25.000,    25.000,    25.000,
     &    25.000,    25.000,    25.000,    25.000,    25.000,
     &    25.000,    32.500,    35.676,     2.162,     2.162,
     &     3.082,     3.082,     3.082,     3.082,    35.629,
     &     1.761,     1.761,     1.761,     1.761,     5.440,
     &     5.440,     5.440,   156.480,     9.067,     9.067,
     &     9.067/


      save RHOP,UBARI,WNOI,DWNI

      REAL effg    ! effg est une fonction(z en m)

c-----------------------------------------------------------------------

c   Initialisations:
c   ----------------

      UBARI2=1./1.66
      UBARI2=UBARI

C ZERO THE FLUXES
         Q0    = 0.0
         zlwup = 0.0 
         zlwdn = 0.0 

c-----------------------------------------------------------------------
C WE NOW ENTER A MAJOR LOOP OVER SPECRAL INTERVALS IN THE INFRARED
C TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL
c-----------------------------------------------------------------------

       DO 2000 K=1,NSPECI    ! *** START OF SPECTRAL LOOP

c-----------------------------------------------------------------------
C SET UP ALTITIDUE PARAMETERS

          WAVEN=WNOI(K)
          DW=DWNI(K)
          zz1=DW/(2.*2)
          EM = 0.
          B0 = 0.

          DO J=1,NL-1
             DO ig=1,NG
                TJ=TEMP(ig,J)

    
C  Modif: in-lining de la fonction planck pour vectorisation
C               B0(ig,J)=PLNCK(WAVEN,TJ,DW)
C     FUNCTION PLNCK(WAV,T,DW)
C* PLNCK FUNCTION RETURNS B IN CGS UNITS, ERGS CM-2 WAVENUMBER-1
C* WAVNUM IS WAVENUMBER IN CM-1
C* T IS IN KELVIN
                PLNCK=0.
                DO I=-2,2,1
                   WAVNUM=WAVEN + I*zz1
                   zz2=EXP(-1.4388 * WAVNUM/TEMP(ig,J))
                   zz3=WAVNUM*WAVNUM*WAVNUM
                   PLNCK=PLNCK+1.191E-5* zz3*zz2/(1.-zz2)
                ENDDO
                B0(ig,J)=.2*PLNCK
             ENDDO

             IF (ICLD.EQ.1) THEN
               DO ig=1,NG
                 zz4=EXP(-DTAUI(ig,J,K)/UBARI2)
                 EM(ig,J)=zz4
               ENDDO
             ELSE
               DO ig=1,NG
                 zz4=EXP(-DTAUIP(ig,J,K)/UBARI2)
                 EM(ig,J)=zz4
               ENDDO
             ENDIF
          ENDDO

c-----------------------------------------------------------------------
C CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
C OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY

           FDI  =0.
           FDIS =0.
           FUPI =0.
           FUPIS=0.

        DO 2220 J=1,NL-1
           DO 2230 ig=1,NG
              FDI(ig,J+1) = FDI(ig,J)*EM(ig,J) + 2.*RPI*UBARI*
     &        B0(ig,J)*(1.-EM(ig,J))
              FDIS(ig,J+1,K) = FDIS(ig,J,K)*EM(ig,J) + 2.*RPI*UBARI*
     &        B0(ig,J)*(1.-EM(ig,J))
2230       CONTINUE
2220    CONTINUE
c     write(*,*)
c     write(*,*) 'cooling : EM  =' ,
c    & ((EM(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : B0  =' ,
c    & ((B0(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : FDI =' ,
c    & ((FDI(i,l),l=1,nl),i=1,ngrid)

c-----------------------------------------------------------------------
C UPWARD FLUXES: SURFACE EMISSIONS

        DO 2310 ig=1,NG
          PLNCK=0.
          DO I=-2,2,1
             WAVNUM=WAVEN + I*zz1
             zz2=EXP(-1.4388 * WAVNUM/TEMP(ig,NL))
             zz3=WAVNUM*WAVNUM*WAVNUM
             PLNCK=PLNCK+1.191E-5* zz3*zz2/(1.-zz2)
          ENDDO
c          BSURF=PLNCK( WAVEN, TEMP(ig,NL), DW)
           BSURF=.2*PLNCK*emis
        FUPI(ig,NL)   =BSURF*2.*RPI*UBARI+(1-emis)*FDI(ig,NL)
        FUPIS(ig,NL,K)=BSURF*2.*RPI*UBARI+(1-emis)*FDIS(ig,NL,K)
2310    CONTINUE
c     write(*,*)
c     write(*,*) 'cooling : FUPI/NL =' ,
c    & ((FUPI(i,l),l=nl,nl),i=1,NG)
c     write(*,*)
c     write(*,*) 'cooling : FDI/NL =' ,
c    & ((FDI(i,l),l=nl,nl),i=1,NG)

        DO 2320 J=NL-1,1,-1
           DO 2330 ig=1,NG
              FUPI(ig,J) = FUPI(ig,J+1)*EM(ig,J) + 2.*RPI*UBARI*
     &        B0(ig,J)*(1.-EM(ig,J))
              FUPIS(ig,J,K) = FUPIS(ig,J+1,K)*EM(ig,J)+2.*RPI*UBARI*
     &        B0(ig,J)*(1.-EM(ig,J))
2330       CONTINUE
2320    CONTINUE
c     write(*,*)
c     write(*,*) 'cooling : EM  =' ,
c    & ((EM(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : B0  =' ,
c    & ((B0(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : FUPI =' ,
c    & ((FUPI(i,l),l=1,nl),i=1,ngrid)

c compute the downward IR flux at the surface:
c 
          DO 3520 ig=1,NG
             pfluxi(ig)=pfluxi(ig)+ DWNI(K)*FDI(ig,NL)
3520      CONTINUE

c compute the up (+ upward) and down (+ downward) IR fluxes:
c 
          DO J=1,NL
          DO ig=1,NG
             zlwup(ig,J)= zlwup(ig,J)+ DWNI(K)*FUPI(ig,J)
             zlwdn(ig,J)= zlwdn(ig,J)+ DWNI(K)*FDI(ig,J)
          ENDDO
          ENDDO
	  
          DO 3210 J=1,NL-1
             DO 3220 ig=1,NG
                QOUT=FUPI(ig,J) + FDI(ig,J+1)   ! OUT OF LAYER
                QIN =FDI(ig,J)  + FUPI(ig,J+1)  ! INTO LAYER
                Q0(ig,J)=Q0(ig,J)+(QOUT-QIN)*DWNI(K)
3220         CONTINUE
3210      CONTINUE

c     write(*,*)
c     write(*,*) 'cooling/loop : FUPI =' ,
c    & ((FUPI(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : FDI  =' ,
c    & ((FDI(i,l),l=1,nl),i=1,ngrid)
c     write(*,*)
c     write(*,*) 'cooling : Q0 =' ,
c    & ((Q0(i,l),l=1,nl-1),i=1,ngrid)


c-----------------------------------------------------------------------

2000  CONTINUE ! *** END SPECTRAL INTERVAL COMPUTATIONS


c-----------------------------------------------------------------------

c   convertion erg/cm2 -> J/m2
      DO 3550 ig=1,NG
         pfluxi(ig)  = 1.e-3*pfluxi(ig)
	 zlwup(ig,:) = 1.e-3*zlwup(ig,:)
	 zlwdn(ig,:) = 1.e-3*zlwdn(ig,:)
3550  CONTINUE

c     PRINT*,'flux IR'
c     WRITE(*,'(8e10.2)') pfluxi

C COMPUTE THE BASELINE COOLING RATE

       DO 3000 J=1,NL-1
C TURN THE Q'S INTO TIMESCALES.....
          DO 3300 ig=1,NG
          COLDEN = RHOP*(PRESS(ig,J+1)-PRESS(ig,J))/effg(Z(ig,J))
c            Q0(J) = (COLDEN * CSUBP )/Q0(J)
             Q0(ig,J) = Q0(ig,J) / (COLDEN*CSUBP) 
3300      CONTINUE
3000  CONTINUE

c     write(*,*)
c     write(*,*) 'cooling/end : Q0 =' 
c     write(*,*) ((Q0(k,l)*1e7,l=1,nl-1),k=1,ngrid)
c-----------------------------------------------------------------------

      RETURN
      END
