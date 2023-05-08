










      SUBROUTINE SFLUXV(DTAUV,TAUV,TAUCUMV,RSFV,DWNV,WBARV,COSBV,
     *                  UBAR0,STEL,NFLUXTOPV,FLUXTOPVDN,
     *                  NFLUXOUTV_nu,NFLUXGNDV_nu,
     *                  FMNETV,FLUXUPV,FLUXDNV,FZEROV,taugsurf)

      use radinc_h
      use radcommon_h, only: tlimit, gweight

      implicit none

      real*8 FMNETV(L_NLAYRAD)
      real*8 TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8 TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8 DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS), DWNV(L_NSPECTV)
      real*8 FMUPV(L_NLAYRAD), FMDV(L_NLAYRAD)
      real*8 COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 STEL(L_NSPECTV)
      real*8 FLUXUPV(L_NLAYRAD), FLUXDNV(L_NLAYRAD)
      real*8 NFLUXTOPV, FLUXUP, FLUXDN,FLUXTOPVDN
      real*8 NFLUXOUTV_nu(L_NSPECTV)
      real*8 NFLUXGNDV_nu(L_NSPECTV)

      integer L, NG, NW, NG1,k
      real*8 ubar0, f0pi, btop, bsurf, taumax, eterm
      real*8 rsfv(L_NSPECTV) ! Spectral dependency added by MT2015.
      real*8 FZEROV(L_NSPECTV)

      real*8 DIFFV, DIFFVT
      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1), fzero

C======================================================================C

      TAUMAX = L_TAUMAX

C     ZERO THE NET FLUXES

      NFLUXTOPV = 0.0
      FLUXTOPVDN = 0.0

      DO NW=1,L_NSPECTV
         NFLUXOUTV_nu(NW)=0.0
         NFLUXGNDV_nu(NW)=0.0
      END DO

      DO L=1,L_NLAYRAD
        FMNETV(L)  = 0.0
        FLUXUPV(L) = 0.0
        FLUXDNV(L) = 0.0
      END DO

      DIFFVT = 0.0

C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE VISIBLE
C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      DO 500 NW=1,L_NSPECTV
      
        F0PI = STEL(NW)

        FZERO = FZEROV(NW)
        IF(FZERO.ge.0.99) goto 40
        DO NG=1,L_NGAUSS-1

          if(TAUGSURF(NW,NG) .lt. TLIMIT) then

            fzero = fzero + (1.0-FZEROV(NW))*GWEIGHT(NG)

            goto 30
          end if

C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE

          BTOP  = 0.0
          !BSURF = 0./0. ! why was this here?
          BSURF = 0.
C         LOOP OVER THE NTERMS BEGINNING HERE
 

!      FACTOR    = 1.0D0 - WDEL(1)*CDEL(1)**2
!      TAU(1)    = TDEL(1)*FACTOR


          ETERM = MIN(TAUV(L_NLEVRAD,NW,NG),TAUMAX)
          BSURF = RSFV(NW)*UBAR0*STEL(NW)*EXP(-ETERM/UBAR0)

C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
C 
C         FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
C         RETURN FLUXES FOR A GIVEN NT


          CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),
     *                WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV(NW),    
     *                BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN)

C         NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

          NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*GWEIGHT(NG)*
     *                          (1.0-FZEROV(NW))
          FLUXTOPVDN = FLUXTOPVDN+FLUXDN*GWEIGHT(NG)*
     *                          (1.0-FZEROV(NW))
          DO L=1,L_NLAYRAD
            FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*
     *                           GWEIGHT(NG)*(1.0-FZEROV(NW))
            FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*GWEIGHT(NG)*
     *                   (1.0-FZEROV(NW))
            FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*GWEIGHT(NG)*
     *                   (1.0-FZEROV(NW))
          END DO

c     band-resolved flux leaving TOA (RDW)
          NFLUXOUTV_nu(NW) = NFLUXOUTV_nu(NW)
     *      +FLUXUP*GWEIGHT(NG)*(1.0-FZEROV(NW))

c     band-resolved flux at ground (RDW)
          NFLUXGNDV_nu(NW) = NFLUXGNDV_nu(NW)
     *      +FMDV(L_NLAYRAD)*GWEIGHT(NG)*(1.0-FZEROV(NW))


C         THE DIFFUSE COMPONENT OF THE DOWNWARD STELLAR FLUX

          DIFFVT = DIFFVT + DIFFV*GWEIGHT(NG)*(1.0-FZEROV(NW))

   30     CONTINUE 

        END DO   ! the Gauss loop 

   40   continue 
C       Special 17th Gauss point

        NG = L_NGAUSS

C       SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
 
        BTOP = 0.0

C       LOOP OVER THE NTERMS BEGINNING HERE
 
        ETERM = MIN(TAUV(L_NLEVRAD,NW,NG),TAUMAX)
        BSURF = RSFV(NW)*UBAR0*STEL(NW)*EXP(-ETERM/UBAR0)


C       WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C       CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C       WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
C 
C       FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
C       RETURN FLUXES FOR A GIVEN NT

        CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),
     *              WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV(NW),
     *              BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN)


C       NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

        NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*FZERO
        FLUXTOPVDN = FLUXTOPVDN+FLUXDN*FZERO
        DO L=1,L_NLAYRAD
          FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*FZERO
          FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*FZERO
          FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*FZERO
        END DO

c     band-resolved flux leaving TOA (RDW)
          NFLUXOUTV_nu(NW) = NFLUXOUTV_nu(NW)
     *      +FLUXUP*FZERO

c     band-resolved flux at ground (RDW)
          NFLUXGNDV_nu(NW) = NFLUXGNDV_nu(NW)+FMDV(L_NLAYRAD)*FZERO


C       THE DIFFUSE COMPONENT OF THE DOWNWARD STELLAR FLUX

        DIFFVT = DIFFVT + DIFFV*FZERO


  500 CONTINUE


C     *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****


      RETURN
      END
