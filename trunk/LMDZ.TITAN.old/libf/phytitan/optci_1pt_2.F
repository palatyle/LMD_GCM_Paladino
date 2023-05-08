      SUBROUTINE optci_1pt2(zqaer_1pt,rcdb,xfrb,iopti,IPRINT)

      use dimphy
      USE TGMDAT_MOD, ONLY: RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,
     &                      RCLOUD,FARGON
      USE TGMDAT_MOD, ONLY: RGAS
#include "dimensions.h"
#include "microtab.h"
#include "numchimrad.h"
#include "clesphys.h"

      PARAMETER(NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECI=46,NSPC1I=47,NSPECV=24,NSPC1V=25)

c   Arguments:
c   ---------
      integer IPRINT,iopti
C iopti: premier appel, on ne calcule qu'une fois les QM et QF
* nrad dans microtab.h
      real   zqaer_1pt(NLAYER,2*nrad)
#include "optci_1pt.h"
c   ---------

      COMMON /ATM/ Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)

      COMMON /GASS/ CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),AR(NLEVEL)
     & ,XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)

      COMMON /STRATO/ C2H2(NLAYER),C2H6(NLAYER)
      COMMON /STRAT2/ HCN(NLAYER)

      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER)
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)

      COMMON /CLOUD/
     &               RCLDI(NSPECI), XICLDI(NSPECI)
     &             , RCLDV(NSPECV), XICLDV(NSPECV)
     &             , RCLDI2(NSPECI), XICLDI2(NSPECI)
     &             , RCLDV2(NSPECV), XICLDV2(NSPECV)

      COMMON /SPECTI/ BWNI(NSPC1I), WNOI(NSPECI),
     &                DWNI(NSPECI), WLNI(NSPECI)

      COMMON /part/v,rayon,vrat,dr,dv

      DIMENSION PROD(NLEVEL)
* nrad dans microtab.h
      real v(nrad),rayon(nrad),vrat,dr(nrad),dv(nrad)
      real xv1(klev,nspeci),xv2(klev,nspeci)
      real xv3(klev,nspeci)
      REAL QF1(nrad,NSPECI),QF2(nrad,NSPECI)
      REAL QF3(nrad,NSPECI),QF4(nrad,NSPECI)
      REAL QM1(nrad,NSPECI),QM2(nrad,NSPECI)
      REAL QM3(nrad,NSPECI),QM4(nrad,NSPECI)
      REAL QC1(nrad,NSPECI),QC2(nrad,NSPECI)
      REAL QC3(nrad,NSPECI),QC4(nrad,NSPECI)
      real emu
      REAL TAEROSM1(NSPECI),TAEROSCATM1(NSPECI),DELTAZM1(NSPECI)

c ---- nuages      
      REAL TNUAGE,TNUAGESCAT
      REAL rcdb(nlayer),xfrb(nlayer,4)

      save qf1,qf2,qf3,qf4,qm1,qm2,qm3,qm4,qc1,qc2,qc3,qc4


C THE PRESSURE INDUCED TRANSITIONS ARE FROM REGIS
C THE LAST SEVENTEEN INTERVALS ARE THE BANDS FROM GNF.
C
C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
C IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
C LAYER: WBAR, DTAU, COSBAR
C LEVEL: TAU
C

        DO 80 K=1,NSPECI
           TAUHI_1pt(K)=0.
           TAUCI_1pt(K)=0.
           TAUGI_1pt(K)=0.
 80     CONTINUE

c************************************************************************
          DO 100 J=1,NLAYER    ! BOUCLE SUR L'ALTITUDE
c************************************************************************

C SET UP THE COEFFICIENT TO REDUCE MASS PATH TO STP ...SEE NOTES
C T0 =273.15   PO=1.01325 BAR

        TBAR=0.5*(TEMP(J)+TEMP(J+1))
        PBAR=SQRT(PRESS(J)*PRESS(J+1))
        BMU=0.5*(XMU(J+1)+XMU(J))
c attention ici, Z en km doit etre passe en m
        COEF1=RGAS*273.15**2*.5E5* (PRESS(J+1)**2 - PRESS(J)**2)
     & /(1.01325**2 *EFFG(Z(J)*1000.)*TBAR*BMU)

      IF (IPRINT .GT. 9) WRITE(6,21) J,EFFG(Z(J)*1000.),TBAR,BMU,COEF1
 21   FORMAT(' J, EFFG, TBAR, BMU, COEF1,: ',I3,1P6E10.3)

c------------------------------------------------------------------------
         DO 101 K=1,NSPECI     ! BOUCLE SUR LES L.D'O
c------------------------------------------------------------------------


C #1:             HAZE
C---------------------

C FIRST COMPUTE TAU AEROSOL


c
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


      DO inq=1,nrad         !BOUCLE SUR LES TAILLE D"AEROSOLS


      IF (WNOI(K).lt.wco) THEN    ! lamda > 56 um

           if (iopti.eq.0) then

c          CALL XMIE(rayon(inq)*1.e6,REALI(K),XIMGI(K),
c    &     QEXT,QSCT,QABS,QBAR,WNOI(K))

           CALL CMIE(1.E-2/WNOI(K),REALI(K),XIMGI(K),rayon(inq),
     &     QEXT,QSCT,QABS,QBAR)


           QM1(inq,K)=QEXT
           QM2(inq,K)=QSCT
           QM3(inq,K)=QABS
           QM4(inq,K)=QBAR

          endif         ! end iopti


      TAEROS=QM1(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4+TAEROS
      TAEROSCAT=QM2(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4+TAEROSCAT
      CBAR=CBAR+QM4(inq,K)*QM2(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4


         ELSE                           ! 0.2 < lambda < 56 um


            if(rayon(inq).lt.RF(inq)) THEN

              if (iopti.eq.0) then

                   CALL XMIE(rayon(inq)*1.e6,REALI(K),XIMGI(K),
     &             QEXT,QSCT,QABS,QBAR,WNOI(K))

              QM1(inq,K)=QEXT
              QM2(inq,K)=QSCT
              QM3(inq,K)=QABS
              QM4(inq,K)=QBAR
              endif         ! end iopti


        TAEROS=QM1(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4+TAEROS
        TAEROSCAT=QM2(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4+TAEROSCAT
        CBAR=CBAR+QM4(inq,K)*QM2(inq,K)*zqaer_1pt(nlayer+1-J,inq)*1.e-4

           else

               XMONO=(rayon(inq)/RF(inq))**3.
               XRULE=1.

            if(XMONO.gt.16384./1.5) then
             XRULE=(XMONO/16384.)
             XMONO=16384.
            endif


             if (iopti.eq.0) then

       CALL CFFFV11(1.e-2/WNOI(K),REALI(K),XIMGI(K),RF(inq),2.
     &                ,XMONO,QSCT,QEXT,QABS,QBAR)

              QF1(inq,K)=QEXT*XRULE
              QF2(inq,K)=QSCT*XRULE
              QF3(inq,K)=QABS*XRULE
              QF4(inq,K)=QBAR
             endif         ! end iopti

        TAEROS=QF1(inq,K)*zqaer_1pt(nlayer+1-J,inq)+TAEROS
        TAEROSCAT=QF2(inq,K)*zqaer_1pt(nlayer+1-J,inq)+TAEROSCAT
        CBAR=CBAR+QF4(inq,K)*QF2(inq,K)*zqaer_1pt(nlayer+1-J,inq)

           endif

               IF(TAEROS.LT.1.e-10) TAEROS=1.e-10

         ENDIF
       ENDDO             ! FIN DE LA BOUCLE SUR nrad

       IF(TAEROSCAT.le.0.) then
        CBAR=0.
       ELSE
        CBAR=CBAR/TAEROSCAT
       ENDIF

        DELTAZ=Z(J)-Z(J+1)

c --------------------------------------------------------------------
c profil brume Pascal: fit T (sauf tropopause) et albedo
c -------------------
        if( cutoff.eq.1) then
         IF(PRESS(J).gt.9.e-3) THEN
          TAEROS=TAEROSM1(K)*DELTAZ/DELTAZM1(K)*0.85
          TAEROSCAT=TAEROSCATM1(K)*DELTAZ/DELTAZM1(K)*0.85
c         TAEROS=0.
c         TAEROSCAT=0.
         ENDIF

         IF(PRESS(J).gt.1.e-1) THEN
          TAEROS=TAEROSM1(K)*DELTAZ/DELTAZM1(K)*1.15
          TAEROSCAT=TAEROSCATM1(K)*DELTAZ/DELTAZM1(K)*1.15
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

         TAEROSM1(K)=TAEROS
         TAEROSCATM1(K)=TAEROSCAT
         DELTAZM1(K)=DELTAZ


      IF(TAEROSCAT.le.0.) CBAR=0.

c     print*,'HERE, MCKAY AEROSOLS IR'
c     TAEROS=xv1(j,k)
c     TAEROSCAT=xv2(j,k)
c     CBAR=xv3(j,k)

c*********** EN TRAVAUX ***************************

C #2:         CLOUD
C------------------
C NEXT COMPUTE TAU CLOUD
c
c  Menu special :
c  Afin d'eviter la surcharge de calcul on ne calcule les 
c  propriétes optiques des nuages qu'une seule fois
c  avec un rayon de particule effectif de 3um et une composition
c  de goutte : 90% CH4 / 10% NOYAUX
c  Puis on ajute les section efficace par la surface reelle de 
c  la goutte.
c
c  ---> A TESTER !!!!
c
        TNUAGE=0.
        TNUAGESCAT=0.
        CNBAR=0.
        IF (clouds.eq.1) THEN
          IF (iopti.eq.0) THEN !--> au premier appel
            QEXTC=0.
            QSCTC=0.
            QABSC=0.
            CBARC=0.
            DO inq=1,nrad         !BOUCLE SUR LES NQMX TAILLE D"AEROSOLS
              QC1(inq,K)=0.
              QC2(inq,K)=0.
              QC3(inq,K)=0.
              QC4(inq,K)=0.
            ENDDO
** OPTICAL CONSTANT : MIXING RULES
** Fraction volumique fixe :
** 10% noyaux.
** 90% methane.
            XNR = 0.5 * REALI(K)
     &          + 0.5 * RCLDI(K)
            XNI = 0.5 * XIMGI(K)
     &          + 0.5 * XICLDI(K)
**
**   Efficacite : particule de 3um de rayon
            CALL CMIE(1.E-2/WNOI(K),XNR,XNI,3.e-6,
     &                QEXTC,QSCTC,QABSC,CBARC)
**
**   ATTENTION CE SONT DES EFFICACITES : il faut les x par la surface REELLE de la goutte.
            DO inq=1,nrad
              QC1(inq,K)=QEXTC/xnuf
              QC2(inq,K)=QSCTC/xnuf
              QC3(inq,K)=QABSC/xnuf
              QC4(inq,K)=CBARC
            ENDDO
          ENDIF   ! iopti = 0

c ----- On ne calcule les constante optiques que si Rgoutte > 1e-10
          IF (rcdb(nlayer+1-J).gt.1.1e-10) THEN
            DO inq=1,nrad        
              TNUAGE=QC1(inq,K)*(rcdb(nlayer+1-J)/3.e-6)**2.*1.e-4*
     &               zqaer_1pt(nlayer+1-J,inq+nrad) +
     &               TNUAGE
              TNUAGESCAT=QC2(inq,K)*(rcdb(nlayer+1-J)/3.e-6)**2.*
     &                   1.e-4*zqaer_1pt(nlayer+1-J,inq+nrad) +
     &                   TNUAGESCAT
              CNBAR=QC4(inq,K)*QC2(inq,K)*
     &              (rcdb(nlayer+1-J)/3.e-6)**2.*
     &              1.e-4*zqaer_1pt(nlayer+1-J,inq+nrad) +
     &              CNBAR
            ENDDO
          ENDIF

          IF(TNUAGESCAT.EQ.0.) THEN
            CNBAR=0.
          ELSE
            CNBAR=CNBAR/TNUAGESCAT
          ENDIF

        ENDIF    ! Cond CLD
c
C #3:          GAZ
C------------------

C NOW COMPUTE TAUGAS DUE TO THE PIA TERM ONLY FOR LAMDA LT 940
       TAUGAS=0.0
       IF (WNOI(K) .LT. 940. ) THEN
                 CALL PIA(K,TBAR,PNN,PCC,PCN,PHN)
C HERE IS WHERE WE COULD SCALE THE PIA COEFFICEINTS TO FIT DATA
C BASED ON REGIS' NOTES. ---TGM HAS THIS ADJUST IN IT AS DEFAULT
                 PCN=PCN*MIN(1.75 , AMAX1(1.0,WNOI(K)/200.))
C***REPLACE ABOVE WITH: PCN=PCN*1.25*MIN(1.75 , AMAX1(1.0,WNOI(K)/200.))
C 1.25 FACTOR (NOT FROM DATA) SUGGESTED BY TOON et al. (1988)
                 TAUGAS=COEF1*
     &           (XN2(J)*XN2(J)*PNN + CH4(J)*CH4(J)*PCC
     &           + XN2(J)*CH4(J)*PCN + XN2(J)*H2(J)*PHN)
            IF (J .EQ. NLAYER .AND. IPRINT .GT. 9)
     &          WRITE (6,22) WNOI(K),TAUGAS,XN2(J),CH4(J),H2(J),
     &          TBAR, PNN,PCC,PCN, PHN,
     &          XN2(J)*XN2(J)*PNN , CH4(J)*CH4(J)*PCC ,
     &          XN2(J)*CH4(J)*PCN , XN2(J)*H2(J)*PHN
 22             FORMAT(1X,1P8E10.2)
       ENDIF

       IF (K .GT. 28) THEN
                KGAS=K-28
C     ??FLAG? HERE MUST BE WATCHED CAREFULLY
                     U=COLDEN(J)*6.02204E23/BMU
                     if((ylellouch).or.(.not.hcnrad)) then
                       CALL GAS2_NOHCN(J, KGAS,TBAR,PBAR,U,TAU2)
                     else
                       CALL GAS2(J, KGAS,TBAR,PBAR,U,TAU2)
                     endif
                     TAUGAS=TAUGAS+TAU2
       ENDIF
C

      DTAUI_1pt(J,K)=TAUGAS+TAEROS+TNUAGE
      DTAUIP_1pt(J,K)=TAUGAS+TAEROS

      TAUHI_1pt(K)=TAUHI_1pt(K) + TAEROS
      TAUHID_1pt(J,K)=TAUHI_1pt(K)

      TAUGI_1pt(K)=TAUGI_1pt(K) + TAUGAS
      TAUGID_1pt(J,K)=TAUGI_1pt(K)

      TAUCI_1pt(K)=TAUCI_1pt(K) + TNUAGE
      TAUCID_1pt(J,K)=TAUCI_1pt(K)
 
C ??FLAG? SERIOUS PROBLEM WITH THE CODE HERE!

      TLIMIT=1.E-16


      IF (TAEROSCAT + TNUAGESCAT .GT. 0.) THEN
         COSBI_1pt(J,K)=(CBAR*TAEROSCAT + CNBAR*TNUAGESCAT )
     &                     /(TAEROSCAT + TNUAGESCAT)
      ELSE
         COSBI_1pt(J,K)=0.0
      ENDIF

      IF (TAEROSCAT  .GT. 0.) THEN
         COSBIP_1pt(J,K)=(CBAR*TAEROSCAT)
     &                     /(TAEROSCAT)
      ELSE
         COSBIP_1pt(J,K)=0.0
      ENDIF

*---------

      IF (DTAUI_1pt(J,K) .GT.  TLIMIT) THEN
          WBARI_1pt(J,K)=(TAEROSCAT+TNUAGESCAT) /DTAUI_1pt(J,K)
      ELSE
         WBARI_1pt(J,K)=0.0
         DTAUI_1pt(J,K)=TLIMIT
      ENDIF

      IF (DTAUIP_1pt(J,K) .GT.  TLIMIT) THEN
          WBARIP_1pt(J,K)=(TAEROSCAT) /DTAUIP_1pt(J,K)
      ELSE
         WBARIP_1pt(J,K)=0.0
         DTAUIP_1pt(J,K)=TLIMIT
      ENDIF


c     IF (IPRINT .GT. 9)
c    & WRITE(6,73)J,K,TAUGAS,TAEROS,QEXT,QSCT
  73           FORMAT(2I3,1P8E10.3)
 

c------------------------------------------------------------------------
 101  CONTINUE   ! FIN BOUCLE L D'O
c------------------------------------------------------------------------
 
      iopti=1

c************************************************************************
 100  CONTINUE   ! FIN BOUCLE ALTITUDE
c************************************************************************
 
        DO 119 K=1,NSPECI
           TAUI_1pt(1,K)=0.0
           TAUIP_1pt(1,K)=0.0
        DO 119 J=1,NLAYER
           TAUI_1pt(J+1,K)=TAUI_1pt(J,K)+DTAUI_1pt(J,K)
           TAUIP_1pt(J+1,K)=TAUIP_1pt(J,K)+DTAUIP_1pt(J,K)
 119    CONTINUE

c      IF (IPRINT .GT. 2) THEN
c          WRITE (6,120)
c  120      FORMAT(///'  OPTICAL CONSTANTS IN THE INFRARED')

c        DO 200 K=1,NSPECI           ! #2
c          WRITE (6,190)
c          WRITE (6,210)K,WLNI(K),WNOI(K),BWNI(K)
c    &    ,BWNI(K)+DWNI(K),DWNI(K)
c          WRITE (6,230)REALI(K),XIMGI(K)

c        DO 195 J=1,NLAYER         !   #3
c          WRITE (6,220)XNUMB(J), WBARI_1pt(J,K),COSBI_1pt(J,K)
c    &                          , DTAUI_1pt(J,K),TAUI_1pt(J,K)
c 195    CONTINUE

c 200    CONTINUE

c          END IF


c  210 FORMAT(1X,I3,F10.3,F10.2,F10.2,'-',F8.2,F10.3)
c  190 FORMAT(1X//'  SNUM  MICRONS   WAVENU   INTERVAL    DELTA-WN')
c  230 FORMAT(1X,'NREAL(LAYER)= ',1PE10.3,' NIMG(LAYER)= ',E10.3/
c     &' #AEROSOLS   WBAR  COSBAR DTAU TAU')
c  220 FORMAT(5(1X,G9.3))
c  240 FORMAT(41X,G9.3)

      RETURN
      END
