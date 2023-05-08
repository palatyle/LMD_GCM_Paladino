      SUBROUTINE optcv_1pt(zqaer_1pt,rcdb,xfrb,ioptv,IPRINT)


      use dimphy
#include "dimensions.h"
#include "microtab.h"
#include "clesphys.h"

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
      REAL QC1(nrad,NSPECV),QC2(nrad,NSPECV)
      REAL QC3(nrad,NSPECV),QC4(nrad,NSPECV)

c---- NUAGES
      real TNUABS,TNUSCAT
      real   rcdb(NLAYER)
      real     xfrb(NLAYER,4)

      save qf1,qf2,qf3,qf4,qm1,qm2,qm3,qm4
 
      integer ilat,jalt
      common/toto/ilat,jalt

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
        jalt=j
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


       CBAR=CBAR/TAEROSCAT

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

c      if (IPRINT.eq.1) then
c      if (k.eq.NSPECV/2) then   
c       write(*,1699) '@VI',K,J,TAEROS,TAEROSCAT,CBAR
c       write(*,1699) '@  ',K,J,QF1(1,K),QF2(1,K),zqaer_1pt(NLAYER+1-J,1)
c       write(*,1699) '@  ',K,J,QF1(3,K),QF2(3,K),zqaer_1pt(NLAYER+1-J,3)
c       write(*,1699) '@  ',K,J,QF1(5,K),QF2(5,K),zqaer_1pt(NLAYER+1-J,5)
c       write(*,1699) '@  ',K,J,QF1(7,K),QF2(7,K),zqaer_1pt(NLAYER+1-J,7)
c       write(*,1699) '@  ',K,J,QF1(9,K),QF2(9,K),zqaer_1pt(NLAYER+1-J,9)
c       print*
c      endif
c      endif

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

c       PRINT*,WLNV(K)
c      COLX=0.
c      COLP=0.
c      COLT=0.
c     DO IU=1,NLAYER
c      COLP=COLDEN(IU)*1.e+1*1.35+COLP
c     TAURAY=(COLDEN(IU)*28.9/(XMU(IU)*1013.25))*
c    & (.008569/WLNV(K)**4)*(1.+.0113/WLNV(K)**2
c    & +.00013/WLNV(K)**4)
c      COLT=COLT+TAURAY
c      COLX=COLDEN(IU)*1.e+1/(1.E5*28./22.4E3)*1.e-1*0.0933e-1+COLX
c                           |    
c                           |    
c           g/cm2->kg/m2    |  m2/kg   
c      Print*,IU, tauray, 
c    &   COLDEN(IU)*1.e+1/(1.E5*28./22.4E3)*1.e-1*0.543e-1
c     ENDDO
c       PRINT*,COLP,' PRESSURE AT GROUND;'
c       PRINT*,COLX,' TAU_GAS AT GROUND;'
c       print*,colt,colx,' COLT, COLX'
c      STOP
							
c       DZ=Z(J)-Z(J+1)
c     PRINT*, Z(J),WLNV(K),
c    &(28.9/(XMU(J)*1013.25))*(.008569/WLNV(K)**4)*
c    &(1.+.0113/WLNV(K)**2+.00013/WLNV(K)**4)
c    & ,COLDEN(J)/DZ/100000.,
c    &(28.9/(XMU(J)*1013.25))*(.008569/WLNV(K)**4)*
c    &(1.+.0113/WLNV(K)**2+.00013/WLNV(K)**4)
c    & *COLDEN(J)/DZ/100000.
    
 

C #3:                   CLOUD
c----------------------------

C NEXT COMPUTE TAU CLOUD

       IF (clouds.eq.0) THEN
         CNBAR=0.
         TNUSCAT=0.
         TNUABS=0.
         TBNUABS=0.
       ELSE
         CNBAR=0.
         TNUSCAT=0.
         TNUABS=0.
         TBNUABS=0.
         QEXTC=0.
         QSCTC=0.
         QABSC=0.
         CBARC=0.

         do inq=1,nrad
           QC1(INQ,k)=0.
           QC2(INQ,k)=0.
           QC3(INQ,k)=0.
           QC4(INQ,k)=0.
         enddo

         IF (rcdb(nlayer+1-J).gt.1.1e-10) THEN

** OPTICAL CONSTANT : MIXING RULES

           XNR=xfrb(nlayer+1-J,1)*REALV(K)                   !
     &        +xfrb(nlayer+1-J,2)*RCLDV(K)                   !
     &        +xfrb(nlayer+1-J,3)*RCLDV2(K)                  !
     &        +xfrb(nlayer+1-J,4)*RCLDV2(K)                  !

           XNI=xfrb(nlayer+1-J,1)*XIMGV(K)
     &        +xfrb(nlayer+1-J,2)*XICLDV(K)
     &        +xfrb(nlayer+1-J,3)*XICLDV2(K)
     &        +xfrb(nlayer+1-J,4)*XICLDV2(K)

** OPTICAL CONSTANT : LIQUID DROP = THOLIN
           IF(xfrb(nlayer+1-J,1).ge.0.01) THEN
             XNI=XIMGV(K)
             XNR=REALV(K)
           ENDIF

           IF (XNI.gt.1.e-10  .and. XNR.gt.1.00) THEN
             CALL CMIE(1.E-2/WNOV(K),XNR,XNI,
     &       rcdb(nlayer+1-J),
     &       QEXTC,QSCTC,QABSC,CBARC)
           ELSE
             PRINT*,' WARNING XNR/XNI in optcv: ',XNR,XNI
             QEXTC=0.
             QSCTC=0.
             QABSC=0.
             CBARC=0.
             STOP
           ENDIF
         ELSE
           QEXTC=0.
           QSCTC=0.
           QABSC=0.
           CBARC=0.
         ENDIF

         DO inq=1,nrad

           QC1(INQ,k)=QEXTC/xnuf
           QC2(INQ,k)=QSCTC/xnuf
           QC3(INQ,k)=QABSC/xnuf
           QC4(INQ,k)=CBARC

           TNUABS=QC1(inq,K)*zqaer_1pt(NLAYER+1-J,inq+nrad)*1.e-4
     &           +TNUABS

           TNUSCAT=QC2(inq,K)*zqaer_1pt(NLAYER+1-J,inq+nrad)*1.e-4
     &            +TNUSCAT

           CNBAR=QC4(inq,K)*QC2(inq,K)*
     &           zqaer_1pt(NLAYER+1-J,inq+nrad)*1.e-4 + CNBAR

         ENDDO

         IF(TNUSCAT.EQ.0) CNBAR=0.
         IF(TNUSCAT.NE.0.) CNBAR=CNBAR/TNUSCAT


       ENDIF  ! Cond. CLD

       TAURV_1pt(K)=TAURV_1pt(K)+TAURAY
       TAUGVD_1pt(J,K)=TAURV_1pt(K)

       TAUHV_1pt(K)=TAUHV_1pt(K)+TAEROS ! INTEGRATED Quant.
       TAUHVD_1pt(J,K)=TAUHV_1pt(K)

       TAUCV_1pt(K)=TAUCV_1pt(K)+TNUABS
       TAUCVD_1pt(J,K)=TAUCV_1pt(K)


C #4:                  TAUGAS
C----------------------------

C LOOP OVER THE NTERMS
C THIS IS THE ABSORPTION BY THE ATMOSPHERE (METHANE)


      DO 909 NT=1,NTERM(K)
      TAUGAS=COLDEN(J)*GAS1(J)*BTERM(NT,K)*
     &  (   (PRESS(J+1) + PRESS(J))*.5  )**PEXPON(K)


*  COSBV ET COSBVP
*-----------------

      IF(TAEROSCAT+TNUSCAT+TAURAY .ne. 0.) THEN
         COSBV_1pt(J,K,NT)=(CBAR*TAEROSCAT + CNBAR*TNUSCAT)
     &  /(TAEROSCAT+TNUSCAT+TAURAY) !CBAR_RAY=0.
      ELSE
         COSBV_1pt(J,K,NT)=0.
      ENDIF

      IF(TAEROSCAT+TAURAY .ne. 0.) THEN
         COSBVP_1pt(J,K,NT)=(CBAR*TAEROSCAT)
     &  /(TAEROSCAT+TAURAY) !CBAR_RAY=0.
      ELSE
         COSBVP_1pt(J,K,NT)=0.
      ENDIF

*  DTAUV ET DTAUVP
*-----------------

      DTAUV_1pt(J,K,NT) =TAUGAS+TAEROS+TAURAY+TNUABS !TAU_ABS_METH
      DTAUVP_1pt(J,K,NT)=TAUGAS+TAEROS+TAURAY       !TAU_ABS_METH

      TAUGV_1pt(K)=TAUGV_1pt(K)+TAUGAS*ATERM(NT,K) !INTEG.

*  WBARV ET WBARVP
*-----------------

      IF(TAUGAS+TAEROS+TAURAY+TNUABS .ne.  0.) THEN
         WBARV_1pt(J,K,NT)=(TAEROSCAT+TAURAY*0.9999999 + TNUSCAT)
     & /(TAUGAS+TAEROS+TAURAY+TNUABS)
      ELSE
         WBARV_1pt(J,K,NT)=0.
      ENDIF

      IF(TAUGAS+TAEROS+TAURAY .ne.  0.) THEN
         WBARVP_1pt(J,K,NT)=(TAEROSCAT+TAURAY*0.9999999 )
     & /(TAUGAS+TAEROS+TAURAY)
      ELSE
         WBARVP_1pt(J,K,NT)=0.
      ENDIF

 909  CONTINUE
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
c          WRITE(6,*) 'latitude:',ig
c          DO 200 K=1,NSPECV
c          WRITE (6,190)
c          WRITE (6,210)K,WLNV(K),WNOV(K),BWNV(K)
c    &    ,BWNV(K)+DWNV(K),DWNV(K)
c          WRITE (6,230)REALV(K),XIMGV(K)
c          DO 195 J=1,NLAYER,NLAYER
C RECALCULATE FOR PRINT OUT ONLY, ONLY FIRST NTERM AT ig=12 (EQUATOR)
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

       if (IPRINT.eq.1) stop

      RETURN
      END
