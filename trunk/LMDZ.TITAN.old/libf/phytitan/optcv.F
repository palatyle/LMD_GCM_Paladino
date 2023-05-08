      SUBROUTINE OPTCV(qaer,nmicro,IPRINT)

      use dimphy
      use infotrac_phy, only: nqtot
      use common_mod, only:rmcbar,xfbar,ncount,TauHVD,TauCVD,TauGVD
      USE TGMDAT_MOD, ONLY: RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,
     &                      RCLOUD,FARGON
#include "dimensions.h"
#include "microtab.h"
#include "clesphys.h"

c   Argument:
c   ---------
      REAL    qaer(klon,klev,nqtot)
      integer nmicro
c   ---------

c  ASTUCE POUR EVITER klon... EN ATTENDANT MIEUX
      INTEGER   ngrid
      PARAMETER (ngrid=(jjm-1)*iim+2)  ! = klon
c
      PARAMETER(NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECI=46,NSPC1I=47,NSPECV=24,NSPC1V=25)

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

      COMMON /TAUS/   TAUHI(ngrid,NSPECI), TAUCI(ngrid,NSPECI)
     &               ,TAUGI(ngrid,NSPECI), TAURV(ngrid,NSPECV)
     &               ,TAUHV(ngrid,NSPECV) ,TAUCV(ngrid,NSPECV)
     &               ,TAUGV(ngrid,NSPECV)

      COMMON /OPTICV/ DTAUV(ngrid,NLAYER,NSPECV,4)
     &               ,TAUV(ngrid,NLEVEL,NSPECV,4)
     &               ,WBARV(ngrid,NLAYER,NSPECV,4)
     &               ,COSBV(ngrid,NLAYER,NSPECV,4)
     &               ,DTAUVP(ngrid,NLAYER,NSPECV,4)
     &               ,TAUVP(ngrid,NLEVEL,NSPECV,4)
     &               ,WBARVP(ngrid,NLAYER,NSPECV,4)
     &               ,COSBVP(ngrid,NLAYER,NSPECV,4)

      COMMON /SPECTV/ BWNV(NSPC1V),WNOV(NSPECV) 
     &               ,DWNV(NSPECV),WLNV(NSPECV)

      COMMON /part/ v(nrad),rayon(nrad),vrat,dr(nrad),dv(nrad)

      REAL xv1(klev,NSPECV)
      REAL xv2(klev,NSPECV)
      REAL xv3(klev,NSPECV)

      REAL QF1(nrad,NSPECV),QF2(nrad,NSPECV)
      REAL QF3(nrad,NSPECV),QF4(nrad,NSPECV)
      REAL QM1(nrad,NSPECV),QM2(nrad,NSPECV)
      REAL QM3(nrad,NSPECV),QM4(nrad,NSPECV)

      save qf1,qf2,qf3,qf4,qm1,qm2,qm3,qm4
 
      integer ioptv,iwarning     ! ioptv: premier appel, une seule boucle sur les l.d'o.
      integer ig_,seulmtunpt
      save ioptv,iwarning,seulmtunpt
      data ioptv,iwarning,seulmtunpt/0,0,0/

      real   zqaer_1pt(NLAYER,2*nrad)
#include "optcv_1pt.h"

      character*100 dummy
      real   dummy2,dummy3

C*
C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISIBLE
C IT CALCULATES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE VIS
C LAYER: WBAR, DTAU, COSBAR
C LEVEL: TAU
C
       sum=0.
       PRINT*,'OPTCV'
       print*,'ATTENTION, TAU UNIFORME DANS OPTCV'

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c INITIALISATIONS UNE SEULE FOIS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (ioptv.eq.0) then

c verif pour taille zqaer_1pt, sachant que si microfi=0 et nqtot=1, 
c il faut quand meme qu'on lise la look-up table de dim nrad=10
c et si microfi=1, on doit avoir nmicro=nrad (dans microtab.h)
c 
c Nouvelle verif pour nuages : 
c La condition ci-dessus n'est plus realisable !
c nmicro comprend maintenant aussi des glaces
c Donc on teste juste que nmicro soit > 2*nrad (ou nrad si on ne fait pas de nuages)
       if (microfi.ge.1) then
         if ((clouds.eq.1).and.(nmicro.lt.2*nrad)) then
           print*,"OPTCV :"
           print*,"clouds = 1 MAIS nmicro < 2*nrad"
           print*,"Probleme pour zqaer_1pt dans optcv."
           stop
         endif
         if ((clouds.eq.0).and.(nmicro.lt.nrad)) then
           print*,"OPTCV :"
           print*,"nmicro < nrad"
           print*,"Probleme pour zqaer_1pt dans optcv."
           stop
         endif
       endif
      
      DO 130 K=1,NSPECV
C LETS USE THE OPTICAL CONSTANTS FOR THOLIN
c     CALL THOLIN(WLNV(K),TNR,TNI)
      CALL THOLIN_CVD(WLNV(K),TNR,TNI)
      REALV(K)=TNR
      XIMGV(K)=TNI*FHVIS
C BUT WE NOW USE THE GEOMETRIC ALBEDO FITTED RESULTS
C      XIMGV(K)=FITEDT(WLNV(K))
C      XIMGV(K)=FITEDN(WLNV(K))
C THE CLOUD IS CLEAR IN THE VISIBLE
      CALL LIQCH4(WLNV(K),TNR,TNI)
      RCLDV(K)=TNR
      XICLDV(K)=TNI
      CALL LIQC2H6(WLNV(K),TNR,TNI)
      RCLDV2(K)=TNR
      XICLDV2(K)=TNI
 130  CONTINUE
C
c      open (unit=1,file='xsetupv')
c       do j=1,nspecv
c        read(1,*) a
c        do i=1,klev
c            read(1,*) xv1(i,j),xv2(i,j),xv3(i,j)
c        enddo
c       enddo
c       close(1)

c DEBUG
c       print*,"wnov=",WNOV

      endif    ! fin initialisations premier appel

c******* DEBUT DES BOUCLE GRILLE ************************
c     PRINT*, 'AEROSOLS EN VISIBLE'

      DO 101 ig=1,klon       !c! BOUCLE SUR GRILLE HORIZONTALE

        if (microfi.ge.1) then
           do iq=1,2*nrad
             if (clouds.eq.0.and.iq.gt.nrad) then
                zqaer_1pt(:,iq)=0.
             else
               do j=1,NLAYER
                  zqaer_1pt(j,iq)=qaer(ig,j,iq)
               enddo
             endif
           enddo
        else
         if (ig.eq.1)  then
c initialisation zqaer_1pt a partir d une look-up table (uniforme en ig)
c boucle sur nrad=10
           open(10,file="qaer_eq_1d.dat")
           do iq=1,15
             read(10,'(A100)') dummy
           enddo
           do j=NLAYER,1,-1
             read(10,*) dummy2,dummy3,(zqaer_1pt(j,iq),iq=1,nrad)
           enddo
           close(10)
         endif
        endif
	
c        if ((ig.eq.klon/2).or.(microfi.eq.0))  then
c       print*,"Q01=",zqaer_1pt(:,1)
c       print*,"Q05=",zqaer_1pt(:,5)
c       print*,"Q10=",zqaer_1pt(:,10)
c       stop
c        endif
	
        iout=0
c       if ((microfi.eq.0).or.(ig.eq.klon/2)) iout=1
        if (seulmtunpt.eq.0) then
          call optcv_1pt3(zqaer_1pt,rmcbar(ig,:),xfbar(ig,:,:),
     &                   ioptv,IPRINT)
           ioptv = 1
	endif

c Pas de microphysique, ni de composition variable: un seul passage
c dans optcv_1pt.
        if ((microfi.eq.0).and.(ylellouch)) then
	   seulmtunpt = 1
	endif
	
        COSBV(ig,:,:,:)= MAX(MIN(COSBV_1pt(:,:,:),0.999999),1e-6)
        WBARV(ig,:,:,:)= MAX(MIN(WBARV_1pt(:,:,:),0.999999),1e-6)
        DTAUV(ig,:,:,:)= DTAUV_1pt(:,:,:) 
        TAUV(ig,:,:,:) = TAUV_1pt(:,:,:) 

        COSBVP(ig,:,:,:)= MAX(MIN(COSBVP_1pt(:,:,:),0.999999),1e-6)
        WBARVP(ig,:,:,:)= MAX(MIN(WBARVP_1pt(:,:,:),0.999999),1e-6)
        DTAUVP(ig,:,:,:)= DTAUVP_1pt(:,:,:) 
        TAUVP(ig,:,:,:) = TAUVP_1pt(:,:,:) 

        TAUHV(ig,:)    = TAUHV_1pt(:) 
        TAUCV(ig,:)    = TAUCV_1pt(:) 
        TAURV(ig,:)    = TAURV_1pt(:) 
        TAUGV(ig,:)    = TAUGV_1pt(:) 

        TauHVD(ig,:,:) = TAUHVD_1pt(:,:) 
        TauCVD(ig,:,:) = TAUCVD_1pt(:,:) 
        TauGVD(ig,:,:) = TAUGVD_1pt(:,:) 

c DEBUG
c     if(ig.eq.(ngrid/2+16)) then
c         print*,ig,'/',KLON,':'
c         print*,'TauHVD_1',TAUHVD(ig,1,:)
c         print*,'TauGVD_1',TAUGVD(ig,1,:)
c         print*,'TauHVD_50',TAUHVD(ig,50,:)
c         print*,'TauGVD_50',TAUGVD(ig,50,:)
c     stop
c     endif

 101  CONTINUE

c FIN BOUCLE GRILLE     *******
c******************************
         
       PRINT*, 'FIN OPTCV'
      RETURN
      END
