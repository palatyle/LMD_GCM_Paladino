      SUBROUTINE OPTCI(ykim,qaer,nmicro,IPRINT)
      use dimphy
      use infotrac_phy, only: nqtot
      use common_mod, only:rmcbar,xfbar,ncount,TauHID,TauCID,TauGID
      USE TGMDAT_MOD, ONLY: RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,
     &                      RCLOUD,FARGON
#include "dimensions.h"
#include "microtab.h"
#include "numchimrad.h"
#include "clesphys.h"

c   Arguments:
c   ---------
      REAL    ykim(klon,klev,nqtot)
      real    qaer(klon,klev,nqtot)
      integer nmicro,IPRINT
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

      COMMON /STRATO/ C2H2(NLAYER),C2H6(NLAYER)
      COMMON /STRAT2/ HCN(NLAYER)

      COMMON /AERSOL/ RADIUS(NLAYER), XNUMB(NLAYER)
     & , REALI(NSPECI), XIMGI(NSPECI), REALV(NSPECV), XIMGV(NSPECV)

      COMMON /CLOUD/
     &               RCLDI(NSPECI), XICLDI(NSPECI)
     &             , RCLDV(NSPECV), XICLDV(NSPECV)
     &             , RCLDI2(NSPECI), XICLDI2(NSPECI)
     &             , RCLDV2(NSPECV), XICLDV2(NSPECV)

      COMMON /TAUS/   TAUHI(ngrid,NSPECI),TAUCI(ngrid,NSPECI),
     &                TAUGI(ngrid,NSPECI),TAURV(ngrid,NSPECV),
     &                TAUHV(ngrid,NSPECV),TAUCV(ngrid,NSPECV),
     &                TAUGV(ngrid,NSPECV)

      COMMON /OPTICI/ DTAUI(ngrid,NLAYER,NSPECI)
     &               ,TAUI (ngrid,NLEVEL,NSPECI)
     &               ,WBARI(ngrid,NLAYER,NSPECI)
     &               ,COSBI(ngrid,NLAYER,NSPECI)
     &               ,DTAUIP(ngrid,NLAYER,NSPECI)
     &               ,TAUIP(ngrid,NLEVEL,NSPECI)
     &               ,WBARIP(ngrid,NLAYER,NSPECI)
     &               ,COSBIP(ngrid,NLAYER,NSPECI)

      COMMON /SPECTI/ BWNI(NSPC1I), WNOI(NSPECI),
     &                DWNI(NSPECI), WLNI(NSPECI)

      REAL DTAUP(ngrid,NLAYER,NSPECI),DTAUPP(ngrid,NLAYER,NSPECI)
      COMMON /IRTAUS/ DTAUP,DTAUPP

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
      real emu
      REAL TAEROSM1(NSPECI),TAEROSCATM1(NSPECI),DELTAZM1(NSPECI)
      
      save qf1,qf2,qf3,qf4,qm1,qm2,qm3,qm4

      integer iopti,iwarning     ! iopti: premier appel, une seule boucle sur les l.d'o.
      integer ig,seulmtunpt,iout
      save iopti,iwarning,seulmtunpt
      data iopti,iwarning,seulmtunpt/0,0,0/

      real   zqaer_1pt(NLAYER,2*nrad)
#include "optci_1pt.h"

      character*100 dummy
      real   dummy2,dummy3

C THE PRESSURE INDUCED TRANSITIONS ARE FROM REGIS
C THE LAST SEVENTEEN INTERVALS ARE THE BANDS FROM GNF.
C
C THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
C IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
C LAYER: WBAR, DTAU, COSBAR
C LEVEL: TAU
C
       print*,'START OPTCI'

c Diagnostic eventuellement:
c      if (nmicro.gt.0) then
c      sum=0.
c      do nng=2,klon
c        do i=1,klev
c         do j=1,nmicro
c          print*,'j,rj',j,rayon(j)
c          print*,'paer',qaer(nng,i,j)
c           sum=sum+qaer(nng,i,j)*rayon(j)**3.*1.3333*3.1415*1000.
c         enddo
c        enddo
c        enddo
c      print*,sum/(klon-1),'SOMME COLONNE/OPTCI'
c      endif


c      do inq=1,nrad
c          print*,inq,rayon(inq),vrat,qaer(12,25,inq)
c      enddo
             
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c INITIALISATIONS UNE SEULE FOIS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (iopti.eq.0) then

c verif pour taille zqaer_1pt, sachant que si microfi=0 et nqtot=1, 
c il faut quand meme qu on lise la look-up table de dim nrad=10
c et si microfi=1, on doit avoir nmicro=nrad (dans microtab.h)
c
c Nouvelle verif pour nuages : 
c La condition ci-dessus n'est plus realisable !
c nmicro comprend maintenant aussi des glaces
c Donc on teste juste que nmicro soit > 2*nrad (ou nrad si on ne fait pas de nuages)
       if (microfi.ge.1) then
         if ((clouds.eq.1).and.(nmicro.lt.2*nrad)) then
           print*,"OPTCI :"
           print*,"clouds = 1 MAIS nmicro < 2*nrad"
           print*,"Probleme pour zqaer_1pt dans optci."
           stop
         endif
         if ((clouds.eq.0).and.(nmicro.lt.nrad)) then
           print*,"OPTCI :"
           print*,"nmicro < nrad"
           print*,"Probleme pour zqaer_1pt dans optci."
           stop
         endif
       endif

      DO 420 K=1,NSPECI
C LETS USE THE THOLIN OPTICAL CONSTANTS FOR THE HAZE.
c         CALL THOLIN(WLNI(K),TNR,TNI)
          CALL THOLIN_CVD(WLNI(K),TNR,TNI)
          REALI(K)=TNR
          XIMGI(K)=TNI*FHIR
C SET UP THE OPTICAL CONSTANTS FOR THE CLOUD
          CALL LIQCH4(WLNI(K),TNR,TNI)
          RCLDI(K)=TNR
          XICLDI(K)=TNI
          CALL LIQC2H6(WLNI(K),TNR,TNI)
          RCLDI2(K)=TNR
          XICLDI2(K)=TNI
 420  CONTINUE

c DEBUG
c       print*,"wnoi=",WNOI

C
C ZERO ALL OPTICAL DEPTHS.
C ??FLAG? FOR SOME APPLCIATIONS THE TOP OPACITY MAY NOT VANISH

c      open  (unit=1,file='xsetupi')
c      do i=1,klev
c       read(1,*) a
c        do j=1,nspeci
c            read(1,*) xv1(i,j),xv2(i,j),xv3(i,j)
c        enddo
c       enddo
c       close(1)

      endif    ! fin initialisations premier appel

c************************************************************************
c************************************************************************
      DO 79 ig=1,klon      ! BOUCLE SUR GRILLE HORIZONTALE    
c      print*,'ig NEW optci',ig
c************************************************************************
c************************************************************************

        if (.not.ylellouch) then
        
            XN2(1) = ykim(ig,1,iradn2)
            CH4(1) = ykim(ig,1,iradch4)
             H2(1) = ykim(ig,1,iradh2)
            do j=2,nlayer
               XN2(j) = (ykim(ig,j,iradn2)+ykim(ig,j-1,iradn2))/2.
               CH4(j) = (ykim(ig,j,iradch4)+ykim(ig,j-1,iradch4))/2.
                H2(j) = (ykim(ig,j,iradh2)+ykim(ig,j-1,iradh2))/2.
            enddo
            XN2(nlevel) = ykim(ig,nlayer,iradn2)
            CH4(nlevel) = ykim(ig,nlayer,iradch4)
             H2(nlevel) = ykim(ig,nlayer,iradh2)     

            do j=1,nlayer
               emu = ( xmu(j) + xmu(j+1) )/2.
               C2H2(j) = ykim(ig,j,iradc2h2) * 26./emu
               C2H6(j) = ykim(ig,j,iradc2h6) * 30./emu
                HCN(j) = ykim(ig,j,iradhcn ) * 27./emu
            enddo
                  
        endif

c     if ((.not.ylellouch).and.(ig.eq.klon/2)) then
c        print*,' LAYER      C2H2         C2H6       HCN masmix ratios'
c        do j=1,nlayer
c            print*,j,C2H2(j),C2H6(j),HCN(j)
c        enddo
c     endif   

        if (microfi.ge.1) then
          do iq=1,2*nrad
c           si on ne fait pas de nuages on ne copie que les nrad premieres valeurs.
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
           zqaer_1pt = 0.
c initialisation zqaer_1pt a partir d une look-up table (uniforme en ig)
c boucle sur nrad=10 (dans microtab.h)
           open(10,file="qaer_eq_1d.dat")
           do iq=1,15
             read(10,'(A100)') dummy
           enddo
           do j=NLAYER,1,-1
             read(10,*) dummy2,dummy3,(zqaer_1pt(j,iq),iq=1,nrad)
           enddo
           close(10)
c ici, les tableaux definissant la structure des aerosols sont
c remplis: rf,df(nq),rayon(nq,)v(nq)......
           call rdf()
         endif
        endif

c        if ((ig.eq.klon/2).or.(microfi.eq.0))  then
c       print*,"Q01=",zqaer_1pt(:,1)
c       print*,"Q05=",zqaer_1pt(:,5)
c       print*,"Q10=",zqaer_1pt(:,10)
c       stop
c        endif
	
        iout=0
c       if ((microfi.eq.0).or.(ig.eq.(klon/2+16))) iout=1
        if (seulmtunpt.eq.0) then
          call optci_1pt3(zqaer_1pt,rmcbar(ig,:),xfbar(ig,:,:),
     &                   iopti,iout)
           iopti = 1
	endif

c Pas de microphysique, ni de composition variable: un seul passage
c dans optci_1pt.
        if ((microfi.eq.0).and.(ylellouch)) then
	   seulmtunpt = 1
	endif
	
        COSBI(ig,:,:)  = MAX(MIN(COSBI_1pt(:,:),0.999999),1e-6) 
        WBARI(ig,:,:)  = MAX(MIN(WBARI_1pt(:,:),0.999999),1e-6) 
        DTAUI(ig,:,:)  = DTAUI_1pt(:,:) 
        TAUI(ig,:,:)   = TAUI_1pt(:,:) 

        COSBIP(ig,:,:)  = MAX(MIN(COSBIP_1pt(:,:),0.999999),1e-6) 
        WBARIP(ig,:,:)  = MAX(MIN(WBARIP_1pt(:,:),0.999999),1e-6) 
        DTAUIP(ig,:,:)  = DTAUIP_1pt(:,:) 
        TAUIP(ig,:,:)   = TAUIP_1pt(:,:) 

        TAUHI(ig,:)    = TAUHI_1pt(:) 
        TAUCI(ig,:)    = TAUCI_1pt(:) 
        TAUGI(ig,:)    = TAUGI_1pt(:) 

        TauHID(ig,:,:) = TAUHID_1pt(:,:) 
        TauCID(ig,:,:) = TAUCID_1pt(:,:) 
        TauGID(ig,:,:) = TAUGID_1pt(:,:) 

c DEBUG
c     if(ig.eq.(ngrid/2+16)) then
c         print*,ig,'/',KLON,':'
c         print*,'TauHID_1',TAUHID(ig,1,:)
c         print*,'TauGID_1',TAUGID(ig,1,:)
c         print*,'TauHID_50',TAUHID(ig,50,:)
c         print*,'TauGID_50',TAUGID(ig,50,:)
c         print*,"DTAUI_1=",DTAUI(ig,1,:)
c         print*,"DTAUI_50=",DTAUI(ig,50,:)
c         print*,'cosBI_1',COSBI(ig,1,:)
c         print*,'cosBI_50',COSBI(ig,50,:)
c         print*,'WBARI_1',WBARI(ig,1,:)
c         print*,'WBARI_50',WBARI(ig,50,:)
c         stop
c     endif

c************************************************************************
c************************************************************************
  79  CONTINUE   ! FIN BOUCLE GRILLE HORIZONTALE
c************************************************************************
c************************************************************************
C THIS ROUTINE HAS ALREADY SET THE DTAUI(J,K) VALUES BUT MUST BE PASSED
        DO 225 IG=1,klon
         DO 220 J=1,NLAYER
          DO 230 K=1,NSPECI
              DTAUP(IG,J,K)=DTAUI(IG,J,K)
              DTAUPP(IG,J,K)=DTAUIP(IG,J,K)
230       CONTINUE
220      CONTINUE
225     CONTINUE

      print*, 'FIN OPTCI'

      RETURN
      END

