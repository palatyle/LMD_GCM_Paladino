










      SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,RSFI,WNOI,DWNI,
     *                  COSBI,WBARI,NFLUXTOPI,NFLUXTOPI_nu,
     *                  FMNETI,fluxupi,fluxdni,fluxupi_nu,
     *                  FZEROI,TAUGSURF)
      
      use radinc_h
      use radcommon_h, only: planckir, tlimit,sigma, gweight
      use comcstfi_mod, only: pi
      
      implicit none
      
      integer NLEVRAD, L, NW, NG, NTS, NTT
      
      real*8 TLEV(L_LEVELS), PLEV(L_LEVELS)
      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 FMNETI(L_NLAYRAD)
      real*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI)
      real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 NFLUXTOPI
      real*8 NFLUXTOPI_nu(L_NSPECTI)
      real*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI)
      real*8 FTOPUP
      
      real*8 UBARI, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
      real*8 PLANCK, PLTOP
      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
      real*8 FZEROI(L_NSPECTI)
      real*8 taugsurf(L_NSPECTI,L_NGAUSS-1), fzero
      
      real*8 fup_tmp(L_NSPECTI),fdn_tmp(L_NSPECTI)
      real*8 PLANCKSUM,PLANCKREF
      
! AB : variables for interpolation
      REAL*8 C1
      REAL*8 C2
      REAL*8 P1
      
!======================================================================C
      
      NLEVRAD = L_NLEVRAD
      
! ZERO THE NET FLUXES
      NFLUXTOPI = 0.0D0
      
      DO NW=1,L_NSPECTI
        NFLUXTOPI_nu(NW) = 0.0D0
        DO L=1,L_NLAYRAD
           FLUXUPI_nu(L,NW) = 0.0D0
           fup_tmp(nw)=0.0D0
           fdn_tmp(nw)=0.0D0
        END DO
      END DO
      
      DO L=1,L_NLAYRAD
        FMNETI(L)  = 0.0D0
        FLUXUPI(L) = 0.0D0
        FLUXDNI(L) = 0.0D0
      END DO
      
! WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
! TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL
      
      TTOP  = TLEV(2)  ! JL12 why not (1) ???
      TSURF = TLEV(L_LEVELS)

      NTS   = int(TSURF*NTfac)-NTstart+1
      NTT   = int(TTOP *NTfac)-NTstart+1

!JL12 corrects the surface planck function so that its integral is equal to sigma Tsurf^4
!JL12 this ensure that no flux is lost due to:
!JL12          -truncation of the planck function at high/low wavenumber
!JL12          -numerical error during first spectral integration
!JL12          -discrepancy between Tsurf and NTS/NTfac
      PLANCKSUM = 0.d0
      PLANCKREF = TSURF * TSURF
      PLANCKREF = sigma * PLANCKREF * PLANCKREF
      
      DO NW=1,L_NSPECTI
! AB : PLANCKIR(NW,NTS) is replaced by P1, the linear interpolation result for a temperature TSURF
         C1 = TSURF * NTfac - int(TSURF * NTfac)
         P1 = (1.0D0 - C1) * PLANCKIR(NW,NTS) + C1 * PLANCKIR(NW,NTS+1)
         PLANCKSUM = PLANCKSUM + P1 * DWNI(NW)
      ENDDO
      
      PLANCKSUM = PLANCKREF / (PLANCKSUM * Pi)
!JL12
      
      DO 501 NW=1,L_NSPECTI
! SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS
! AB : PLANCKIR(NW,NTS) is replaced by P1, the linear interpolation result for a temperature TSURF
! AB : idem for PLANCKIR(NW,NTT) and PLTOP
         C1 = TSURF * NTfac - int(TSURF * NTfac)
         C2 = TTOP  * NTfac - int(TTOP  * NTfac)
         P1 = (1.0D0 - C1) * PLANCKIR(NW,NTS) + C1 * PLANCKIR(NW,NTS+1)
         BSURF = (1. - RSFI) * P1 * PLANCKSUM
         PLTOP = (1.0D0 - C2) * PLANCKIR(NW,NTT) + C2*PLANCKIR(NW,NTT+1)
         
! If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
! special Gauss point at the end.
         FZERO = FZEROI(NW)
         
         IF(FZERO.ge.0.99) goto 40
         
         DO NG=1,L_NGAUSS-1
            
            if(TAUGSURF(NW,NG).lt. TLIMIT) then
               fzero = fzero + (1.0D0-FZEROI(NW))*GWEIGHT(NG)
               goto 30
            end if
            
! SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
! CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
! OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
            
!            TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
            TAUTOP = TAUCUMI(2,NW,NG)
            BTOP   = (1.0D0-EXP(-TAUTOP/UBARI))*PLTOP
            
! WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
! CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
! WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
            
            CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)
         
! NOW CALCULATE THE CUMULATIVE IR NET FLUX
            NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)
     *                * (1.0D0-FZEROI(NW))
            
! and same thing by spectral band... (RDW)
            NFLUXTOPI_nu(NW) = NFLUXTOPI_nu(NW) + FTOPUP * DWNI(NW)
     *                       * GWEIGHT(NG) * (1.0D0-FZEROI(NW))
            
            DO L=1,L_NLEVRAD-1
!           CORRECT FOR THE WAVENUMBER INTERVALS
               FMNETI(L)  = FMNETI(L) + (FMUPI(L)-FMDI(L)) * DWNI(NW)
     *                    * GWEIGHT(NG)*(1.0D0-FZEROI(NW))
               FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)
     *                    * (1.0D0-FZEROI(NW))
               FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)
     *                    * (1.0D0-FZEROI(NW))
!         and same thing by spectral band... (RW)
               FLUXUPI_nu(L,NW) = FLUXUPI_nu(L,NW) + FMUPI(L)*DWNI(NW)
     *                          * GWEIGHT(NG) * (1.0D0 - FZEROI(NW))
            END DO
            
   30       CONTINUE
         
         END DO       !End NGAUSS LOOP
         
   40    CONTINUE
         
! SPECIAL 17th Gauss point
         NG     = L_NGAUSS
         
!         TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
         TAUTOP = TAUCUMI(2,NW,NG)
         BTOP   = (1.0D0-EXP(-TAUTOP/UBARI))*PLTOP
         
! WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
! CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
! WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
         
         CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)
         
! NOW CALCULATE THE CUMULATIVE IR NET FLUX
         NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO
         
!         and same thing by spectral band... (RW)
         NFLUXTOPI_nu(NW) = NFLUXTOPI_nu(NW)
     *      +FTOPUP*DWNI(NW)*FZERO
         
         DO L=1,L_NLEVRAD-1
! CORRECT FOR THE WAVENUMBER INTERVALS
            FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
            FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
            FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO
! and same thing by spectral band... (RW)
            FLUXUPI_nu(L,NW) = FLUXUPI_nu(L,NW)
     *                       + FMUPI(L) * DWNI(NW) * FZERO
         END DO
         
  501 CONTINUE      !End Spectral Interval LOOP
! *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****
      
      RETURN
      END
