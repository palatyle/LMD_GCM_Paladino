      SUBROUTINE optcv_1pt3(zqaer_1pt,rcdb,xfrb,ioptv,IPRINT)

      use dimphy
      USE TGMDAT_MOD, ONLY: RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,
     &                      RCLOUD,FARGON
      IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "clesphys.h"

      integer nlayer, nlevel, nspeci, nspc1i, nspecv, nspc1v,nterm
      real z,press, den, temp, ch4, xn2, h2, ar, xmu, gas1,
     &     colden, c2h2, c2h6, hcn, radius, xnumb, reali,
     &     ximgi, realv, ximgv, rcldi, xicldi, rcldv, xicldv, rcldi2,
     &     xicldi2, rcldv2, xicldv2
      real bwni, wnoi, dwni, wlni,prod

      integer k, j,inq,nt

      real tbar, pbar, bmu, coef1, effg, taeros, taeroscat, cbar,
     &     qext, qsct, qabs, qbar, xmono, xrule, deltaz, 
     &     cnbar, qextc, qsctc, qabsc, qbarc, taugas, pnn,
     &     pcc, pcn, phn, kgas, u, ig, tau2, tlimit,
     &     solarf, pexpon, aterm, bterm, bwnv, wnov, dwnv,
     &     wlnv, v, rayon, vrat, dr, dv, taerosm1, deltazm1, 
     &     taeroscatm1, tauray

      PARAMETER(NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECI=46,NSPC1I=47,NSPECV=24,NSPC1V=25)

c   Arguments:
c   ---------
      integer IPRINT,ioptv
C ioptv: premier appel, on ne calcule qu'une fois les QM et QF
* nrad dans microtab.h
      real   zqaer_1pt(NLAYER,2*nrad)
#include "optcv_1pt.h"
c   ---------

      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)

      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),AR(NLEVEL)
     & ,XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)

      COMMON /VISGAS/SOLARF(NSPECV),NTERM(NSPECV),PEXPON(NSPECV),
     &         ATERM(4,NSPECV),BTERM(4,NSPECV)

      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER)
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)

      COMMON /CLOUD/ 
     &               RCLDI(NSPECI), XICLDI(NSPECI)
     &             , RCLDV(NSPECV), XICLDV(NSPECV)
     &             , RCLDI2(NSPECI), XICLDI2(NSPECI)
     &             , RCLDV2(NSPECV), XICLDV2(NSPECV)

      COMMON /SPECTV/ BWNV(NSPC1V),WNOV(NSPECV) 
     &               ,DWNV(NSPECV),WLNV(NSPECV)

* nrad dans microtab.h
      COMMON /part/ v(nrad),rayon(nrad),vrat,dr(nrad),dv(nrad)

      REAL QF1(nrad,NSPECV),QF2(nrad,NSPECV)
      REAL QF3(nrad,NSPECV),QF4(nrad,NSPECV)
      REAL QM1(nrad,NSPECV),QM2(nrad,NSPECV)
      REAL QM3(nrad,NSPECV),QM4(nrad,NSPECV)

c---- NUAGES
      real TNUEXT,TNUSCAT
      real   rcdb(NLAYER), xfrb(NLAYER,4)

      save qf1,qf2,qf3,qf4,qm1,qm2,qm3,qm4


C*
C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISIBLE
C IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE VIS
C LAYER: WBAR, DTAU, COSBAR
C LEVEL: TAU
C
C ZERO THE COLUMN OPTICAL DEPTHS OF EACH TYPE
C ??FLAG? THE OPTICAL DEPTH OF THE TOP OF THE MODEL
C MAY NOT BE ZERO.

c******* DEBUT DES BOUCLES ************************
      DO 100 K=1,NSPECV         !b! BOUCLE SUR LAMBDA

      TAURV_1pt(K)=0.
      TAUHV_1pt(K)=0.            ! INTEGRATED TAU.......INITIALIZATION.
      TAUCV_1pt(K)=0.            ! Rayleigh, Haze, Cloud, Gas
      TAUGV_1pt(K)=0.            !   sca,    abs,  abs  , abs

      DO 100 J=1,NLAYER         !a! BOUCLE SUR L"ALTITUDE

C #1:                   HAZE
c---------------------------

c     CALL THE MIE CODE TO GIVE THE AEROSOL PROPERTIES
c     USE XFRAC FOR FRACTAL AEROSOLS PROPERTIES AT LAMBDA < 2. um 




c                    /\
c                   /  \
c                  /    \
c                 / _O   \
c                / |/     \
c               /  / \     \
c              /   |\ \/\   \
c             /    || /  \   \
c             ----------------
c            |     WARNING    |
c            |    SLOW DOWN   |
c             ---------------- 




c*********** EN TRAVAUX ***************************
 
         TAEROS=0.
         TAEROSCAT=0.
         CBAR=0.

c       print*,"rayon=",rayon
c       print*,"RF=",RF

      DO inq=1,nrad         !BOUCLE SUR LES TAILLE D"AEROSOLS


            IF (rayon(inq).lt.RF(inq)) THEN    ! aerosols spheriques

              
            if(ioptv.eq.0.and.J.eq.1) then
c                  CALL XMIE(rayon(inq)*1.e6,REALV(K),XIMGV(K),
c    &             QEXT,QSCT,QABS,QBAR,WNOV(K))

                CALL CMIE(1.E-2/WNOV(K),REALV(K),XIMGV(K),rayon(inq),
     &          QEXT,QSCT,QABS,QBAR)

c       print*,'inq=',inq,' QM1=',QM1(inq,K),' QEXT=',QEXT

              QM1(inq,K)=QEXT
              QM2(inq,K)=QSCT
              QM3(inq,K)=QABS
              QM4(inq,K)=QBAR
            endif

       TAEROS=QM1(inq,K)*zqaer_1pt(NLAYER+1-J,inq)*1.e-4+TAEROS
       TAEROSCAT=QM2(inq,K)*zqaer_1pt(NLAYER+1-J,inq)*1.e-4+TAEROSCAT
       CBAR=CBAR+QM4(inq,K)*QM2(inq,K)*zqaer_1pt(NLAYER+1-J,inq)*1.e-4

            ELSE                        ! aerosols fractals

               XMONO=(rayon(inq)/RF(inq))**3.
               XRULE=1. 

            if(XMONO.gt.16384./1.5) then 
             XRULE=(XMONO/16384.) 
             XMONO=16384.
            endif 

            if(ioptv.eq.0.and.J.eq.1) then

c        CALL OPTFRAC(XMONO,10000./WNOV(K)
c     &                        ,QEXT,QSCT,QABS,QBAR)

        CALL CFFFV11(1.e-2/WNOV(K),REALV(K),XIMGV(K),RF(inq),2.
     &   ,XMONO,QSCT,QEXT,QABS,QBAR)


              QF1(inq,K)=QEXT*XRULE
              QF2(inq,K)=QSCT*XRULE
              QF3(inq,K)=QABS*XRULE
              QF4(inq,K)=QBAR

c       print*,'inq=',inq,' QF1=',QF1(inq,K),' QEXT=',QEXT,' XRULE=',XRULE
                   
            endif

        TAEROS=QF1(inq,K)*zqaer_1pt(NLAYER+1-J,inq)+TAEROS
        TAEROSCAT=QF2(inq,K)*zqaer_1pt(NLAYER+1-J,inq)+TAEROSCAT
        CBAR=CBAR+QF4(inq,K)*QF2(inq,K)*zqaer_1pt(NLAYER+1-J,inq)

           ENDIF 

       ENDDO    ! nrad


       if (TAEROSCAT.ne.0.) CBAR=CBAR/TAEROSCAT

        DELTAZ=Z(J)-Z(J+1)

c --------------------------------------------------------------------
c profil brume Pascal: fit T (sauf tropopause) et albedo
c -------------------
        if( cutoff.eq.1) then
         IF(PRESS(J).gt.9.e-3) THEN
          TAEROS=TAEROSM1*DELTAZ/DELTAZM1*0.85
          TAEROSCAT=TAEROSCATM1*DELTAZ/DELTAZM1*0.85
c         TAEROS=0.
c         TAEROSCAT=0.
         ENDIF

         IF(PRESS(J).gt.1.e-1) THEN
          TAEROS=TAEROSM1*DELTAZ/DELTAZM1*1.15
          TAEROSCAT=TAEROSCATM1*DELTAZ/DELTAZM1*1.15
c         TAEROS=0.
c         TAEROSCAT=0.
         ENDIF
        endif !cutoff=1

c profil brume pour fit T (y compris tropopause), mais ne fit plus albedo...
c -----------------------
        if( cutoff.eq.2) then
         IF(PRESS(J).gt.1.e-1) THEN
          TAEROS=0.
          TAEROSCAT=0.
         ENDIF
        endif !cutoff=2
c --------------------------------------------------------------------

         TAEROSM1=TAEROS
         TAEROSCATM1=TAEROSCAT
         DELTAZM1=DELTAZ


       IF (TAEROSCAT.le.0.) CBAR=0. 

c     print*, 'HERE, CIRS AEROSOLS'
c     call cirs_haze(PRESS(J),WNOV(K),TAEROS,TAEROSCAT,CBAR)

1699  FORMAT(a3,2I3,3(ES15.7,1X))

c*********** EN TRAVAUX ***************************

C #2:                   RAYLEIGH
c-------------------------------

C RAYLEIGH SCATTERING STRAIGHT FROM HANSEN AND TRAVIS...SEE NOTES
C RATIOED BY THE LAYER COLUMN NUMBER TO THE TOTAL
C COLUMN NUMBER ON EARTH. CM-2
C THIS IS THE SCATTERING BY THE ATMOSPHERE

      TAURAY=(COLDEN(J)*28.9/(XMU(J)*1013.25))*
     &(.008569/WLNV(K)**4)*(1.+.0113/WLNV(K)**2+.00013/WLNV(K)**4)


C #3:                   CLOUD
c----------------------------
C NEXT COMPUTE TAU CLOUD
c 
c  Menu special :
c  On utilise ici une look-up table afin de calculer
c  les proprietes optique des nuages.
c  Le principe est le suivant :
c  La look-up table contient les proprietes optique d'une goutte
c  de methane pur de 3 um.
c  On approxime les proprietes optiques pour une goutte de rayon r a
c  de la table.
c
       TNUEXT=0.
       TNUSCAT=0.
       CNBAR=0.
       IF (clouds.eq.1) THEN

         CALL getoptcld(1.E-2/WNOV(K),rcdb(nlayer+1-J),
     &                  QEXTC,QSCTC,QABSC,QBARC) 
         TNUEXT=0.
         TNUSCAT=0.
         CNBAR=0.
         IF (rcdb(nlayer+1-J).gt.1.1e-10) THEN
           TNUEXT =QEXTC/xnuf*SUM(zqaer_1pt(NLAYER+1-J,nrad+1:2*nrad))
           TNUSCAT=QSCTC/xnuf*SUM(zqaer_1pt(NLAYER+1-J,nrad+1:2*nrad))
           CNBAR  =QBARC 
         ENDIF
           IF(TNUSCAT.GE.0.8*TNUEXT) TNUSCAT=0.8*TNUEXT
       ENDIF  ! Cond. CLD

       TAUCV_1pt(K)=TAUCV_1pt(K)+TNUEXT
       TAUCVD_1pt(J,K)=TAUCV_1pt(K)

       TAURV_1pt(K)=TAURV_1pt(K)+TAURAY
       TAUGVD_1pt(J,K)=TAURV_1pt(K)

       TAUHV_1pt(K)=TAUHV_1pt(K)+TAEROS ! INTEGRATED Quant.
       TAUHVD_1pt(J,K)=TAUHV_1pt(K)



C #4:                  TAUGAS
C----------------------------

C LOOP OVER THE NTERMS
C THIS IS THE ABSORPTION BY THE ATMOSPHERE (METHANE)


       DO 909 NT=1,NTERM(K)
         TAUGAS=COLDEN(J)*GAS1(J)*BTERM(NT,K)*
     &   (   (PRESS(J+1) + PRESS(J))*.5  )**PEXPON(K)


*  COSBV ET COSBVP
*-----------------

         IF(TAEROSCAT+TNUSCAT+TAURAY .ne. 0.) THEN
           COSBV_1pt(J,K,NT)=(CBAR*TAEROSCAT + CNBAR*TNUSCAT)
     &     /(TAEROSCAT+TNUSCAT+TAURAY) !CBAR_RAY=0.
         ELSE
           COSBV_1pt(J,K,NT)=0.
         ENDIF

         IF(TAEROSCAT+TAURAY .ne. 0.) THEN
           COSBVP_1pt(J,K,NT)=(CBAR*TAEROSCAT)
     &     /(TAEROSCAT+TAURAY) !CBAR_RAY=0.
         ELSE
           COSBVP_1pt(J,K,NT)=0.
         ENDIF

*  DTAUV ET DTAUVP
*-----------------

         DTAUV_1pt(J,K,NT) =TAUGAS+TAEROS+TAURAY+TNUEXT !TAU_ABS_METH
         DTAUVP_1pt(J,K,NT)=TAUGAS+TAEROS+TAURAY       !TAU_ABS_METH

         TAUGV_1pt(K)=TAUGV_1pt(K)+TAUGAS*ATERM(NT,K) !INTEG.

*  WBARV ET WBARVP
*-----------------

         IF(TAUGAS+TAEROS+TAURAY+TNUEXT .ne.  0.) THEN
           WBARV_1pt(J,K,NT)=(TAEROSCAT+TAURAY*0.9999999 + TNUSCAT)
     &     /(TAUGAS+TAEROS+TAURAY+TNUEXT)
         ELSE
           WBARV_1pt(J,K,NT)=0.
         ENDIF

         IF(TAUGAS+TAEROS+TAURAY .ne.  0.) THEN
           WBARVP_1pt(J,K,NT)=(TAEROSCAT+TAURAY*0.9999999 )
     &     /(TAUGAS+TAEROS+TAURAY)
         ELSE
           WBARVP_1pt(J,K,NT)=0.
         ENDIF
    
 909   CONTINUE
 
       TAUGVD_1pt(J,K)=TAUGVD_1pt(J,K)+TAUGV_1pt(K)

 100  CONTINUE

       ioptv=1

c HERE END OF THE LOOPS *******
c******************************
         
C TOTAL EXTINCTION OPTICAL DEPTHS
          DO 119 K=1,NSPECV
C LOOP OVER NTERMS
           DO 119 NT=1,NTERM(K)
           TAUV_1pt(1,K,NT)=0.0
           TAUVP_1pt(1,K,NT)=0.0
             DO 119 J=1,NLAYER
             TAUV_1pt(J+1,K,NT)=TAUV_1pt(J,K,NT)+DTAUV_1pt(J,K,NT)
             TAUVP_1pt(J+1,K,NT)=TAUVP_1pt(J,K,NT)+DTAUVP_1pt(J,K,NT)
 119     CONTINUE


c       print*,'SETUP'
c      do i=1,NSPECV
c      print*,WLNV(i)
c       do j=1,NLAYER+1
c       print*,Z(j),TAUV(1,j,i,1),WBARV(1,j,i,1),COSBV(1,j,i,1)
c       enddo
c      enddo
c
c     IF (IPRINT .GT. 1) THEN
c           NT=1
c     IF (2 .GT. 1) THEN
c          WRITE (6,120)
c 120      FORMAT(///'  OPTICAL CONSTANTS IN THE VISIBLE (@EQUATOR) ')
c          DO 200 K=1,NSPECV
c          WRITE (6,190)
c          WRITE (6,210)K,WLNV(K),WNOV(K),BWNV(K)
c    &    ,BWNV(K)+DWNV(K),DWNV(K)
c          WRITE (6,230)REALV(K),XIMGV(K)
c          DO 195 J=1,NLAYER,NLAYER
c          WRITE (6,220)XNUMB(J), WBARV_1pt(J,K,NT),COSBV_1pt(J,K,NT)
c    &      ,DTAUV_1pt(J,K,NT),TAUV_1pt(J,K,NT)
c 195      CONTINUE
c          WRITE (6,240) TAUV_1pt(NLEVEL,K,NT)
c 200      CONTINUE
c     END IF

c  210 FORMAT(1X,I3,F10.3,F10.2,F10.2,'-',F8.2,F10.3)
c  190 FORMAT(1X//'  SNUM  MICRONS   WAVENU   INTERVAL    DELTA-WN')
c  230 FORMAT(1X,'NREAL(LAYER)= ',1PE10.3,' NIMG(LAYER)= ',E10.3/
c     &' #AEROSOLS   WBAR  COSBAR       DTAU     TAU'
c     & ,9X,'RAY     GAS    AEROSOL')
c  220 FORMAT(8(1X,F9.3))
c  240 FORMAT(41X,F9.3)

      RETURN
      END
