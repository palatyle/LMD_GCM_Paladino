       SUBROUTINE RADTITAN(p,nq,nmicro,ycomp,qaer)

c=======================================================================
c
c   Authors: C.P. Mc Kay  01/02/91
c   -------
c
c   Object:  Computation of the solar and infra-red
c   -------  Opacities    (dans des common...)
c
c ON TITAN       ADAPTED FROM BEST.FOR FEB 91
c                           C.P. McKAY
c
c   Arguments:
c   ----------
c
c      Input:
c      ------
c
c        p(klon,nl)    pressure (level)
c        nq       nombre de traceurs
c        nmicro   nombre de traceurs microphysiques
c        ycomp(klon,nlayer,nq) 
c
c      Output:
c      -------
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      USE infotrac_phy, ONLY: tname
      use dimphy
      USE geometry_mod, ONLY: latitude ! in radians
      USE optcld, only : iniqcld
      use moyzon_mod, only:plevmoy
      USE TGMDAT_MOD, ONLY: RHCH4,FH2,FHAZE,FHVIS,FHIR,TAUFAC,
     &                      RCLOUD,FARGON
      USE TGMDAT_MOD, ONLY: RHOP
      IMPLICIT NONE
#include "dimensions.h"
#include "clesphys.h"
#include "microtab.h"
#include "numchimrad.h"
#include "YOMCST.h"

c Pour le CRAY, les block data doivent etre declares external
c pour etre pris en compte
      EXTERNAL TGMDAT

      INTEGER NLEVEL,NLAYER,NSPECI,NSPC1I,NSPECV,NSPC1V,NSPV
      PARAMETER(NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECI=46,NSPC1I=47,NSPECV=24,NSPC1V=25)
      PARAMETER (NSPV=21)  ! LDO POUR CALCUL ALBEDO

c
c  ASTUCE POUR EVITER klon... EN ATTENDANT MIEUX
      INTEGER   ngrid
      PARAMETER (ngrid=(jjm-1)*iim+2)  ! = klon
c

c   Arguments:
c   ----------

      INTEGER nq,nmicro  

      REAL p(klon,nlevel)
      REAL ycomp(klon,nlayer,nq)
      REAL qaer(klon,klev,nq)

c   Local:
c   ------

      INTEGER I,J,IG,K,IPRINT
      INTEGER IPREM
      LOGICAL notfirstcall
      SAVE IPREM,notfirstcall
      data notfirstcall/.false./

      REAL emu,somcoslat,coslat(ngrid)
 
      REAL PCH4, effg,FH2L,RHCH4L,SSUM    ! effg est une fonction(z)

c   COMMONS for interface with local subroutines:
c   ---------------------------------------------

      REAL DZED(NLAYER)
      REAL Z(NLEVEL),PRESS(NLEVEL),DEN(NLEVEL),TEMP(NLEVEL)
      REAL  CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),AR(NLEVEL)
      REAL  XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)
      REAL  C2H2(NLAYER),C2H6(NLAYER),HCN(NLAYER)

      COMMON /VERTICAL/ DZED

      COMMON /ATM/ Z,PRESS
     &            ,DEN,TEMP


      COMMON /GASS/ CH4,XN2
     &              ,H2,AR
     &              ,XMU,GAS1
     &              ,COLDEN

      COMMON /STRATO/ C2H2,C2H6
      COMMON /STRAT2/ HCN

c-----------------------------------------------------------------------
c   1. Initialisations:
c   -------------------


C IPRINT CONTOLS OUTPUT AMOUNT:0=IRREDUCIBLE OUTPUT,LESS THAN 1 PAGE
C PER RUN, 0=MINIMAL OUTPUT, 1=BACKGROUND ATM AND SPEC; 10=FULL DEBUG
      IPRINT=1

C&&
      FHAZE=0.3
C&&
       if(iprem.eq.0) then
         TAUFAC=0
c xvis et xir lus dans physiq.def  (ancien fichier initpar)
       FHVIS= xvis
       FHIR = xir
c      on initialise le paquet optcld
       if (clouds.eq.1) call iniqcld()
       iprem=1
       endif

c-----------------------------------------------------------------------
c   2. Calcul of the atmospheric profile:
c   -------------------------------------

       print*,'dans radtitan ',klon
       print*,notfirstcall
       IF(notfirstcall) GOTO 300           !F au premier appel!
       print*,notfirstcall

c   pression moyenne globale
c   passage au pressions en bar avec indice 1 au sommet.
c   (similaire zp dans radlwsw)
      DO 210 J=2,NLEVEL
         PRESS(J)=plevmoy(NLEVEL+1-j)*1.e-5
210   CONTINUE
      PRESS(1) = PRESS(2)*0.001

c  a cause du tableau predefini dans lell.F (et lell_light.F)
c     IF(press(nlevel-1).GE.1.44) then
      IF(press(nlevel-1).GE.1.48) then 
           STOP'pression au sol trop grande'
          PRINT*,'pression au sol trop grande'
      endif

c      PRESS(nlevel)=1.48
c      XCORR=1.48/PRESS(nlevel)
c     DO 211 J=1,NLEVEL
c        PRESS(J)=XCORR*PRESS(J)
c11   CONTINUE

c *********************************************************
c + 20/1/00: S.Lebonnois: model with chemistry
c   ++ 22/07/02: ajout HCN ++
c *********************************************************
      if (ylellouch) then
c------------------------------------------------------
c  initialisation de l'atmosphere et de la composition
c------------------------------------------------------
          CALL LELL(NLEVEL,Z,RHCH4L,FH2L,FARGON,TEMP,PRESS,DEN,XMU,
     &              CH4,H2,XN2,AR,IPRINT)
 
          print*,'LELLOUCH'
          do i=1,55
             print*,z(i),PRESS(i)
          enddo
C
C
C    NOW CALCULATE THE LAYER AVERAGE GAS MIXING RATIOS.
          CALL GASSES(IPRINT)

      else
c------------------------------------------------------
c  initialisation seulement de l'atmosphere
c------------------------------------------------------
          CALL LELL_LIGHT(NLEVEL,Z,FARGON,TEMP,PRESS,DEN,XMU,
     &              CH4,H2,XN2,AR,IPRINT)
 
          print*,'LELLOUCH LIGHT'
          do i=1,55
             print*,z(i),PRESS(i)
          enddo

c ++ remplace gasses.F ++

          do i=1,nq
             if (tname(i).eq."CH4") then
                iradch4=i
             elseif (tname(i).eq."C2H2") then
                iradc2h2=i
             elseif (tname(i).eq."C2H6") then
                iradc2h6=i
             elseif (tname(i).eq."HCN") then
                iradhcn=i
             elseif (tname(i).eq."N2") then
                iradn2=i
             elseif (tname(i).eq."H2") then
                iradh2=i
             endif 
          enddo
          
c          print*,iradch4,iradc2h2,iradc2h6,iradhcn,iradn2,iradh2
          
          print*,' ALT   CH4 mass mixing ratio '
          
          somcoslat=0.
          do j=1,klon
	    coslat(j) = cos(latitude(j))
            somcoslat=somcoslat+coslat(j)
          enddo
          do i=1,nlayer
c attention ici, Z en km doit etre passe en m
             colden(i)=rhop*(press(i+1)-press(i))/effg(z(i)*1000.)
             gas1(i)=0.
             emu=(xmu(i+1)+xmu(i))/2.
             do j=1,klon
               gas1(i) = gas1(i) + 
     $            coslat(j)/somcoslat*ycomp(j,i,iradch4)*(16./emu)
             enddo
             print*,z(i),gas1(i)
          enddo
          
          RHCH4=0.
          do j=1,klon
             RHCH4 = RHCH4 + coslat(j)/somcoslat*ycomp(j,nlayer,iradch4)
          enddo
          RHCH4 = RHCH4*press(nlevel)/PCH4(temp(nlevel))
          print*,'RHCH4 = ',RHCH4

      endif

c *********************************************************

C
C CALL A ROUTINE THAT SETS UP THE IR SPECTRAL INTERVALS
      CALL SETSPI(IPRINT)
      CALL SETSPV(IPRINT)
C SET UP PIA COEFFICIENTS
      CALL SETPIA(IPRINT,1)

      IF (TAUFAC .GT. 0.)  CALL  CLD(IPRINT)

C
C CALL A SUBROUTINE THAT SETS UP THE OPTICAL PROPERTIES IN THE
C  INFRARED. AND THEN IN THE VISIBLE.

C  NOW, THIS COMPUTATION IS DONE FOR EACH VALUE OF klon
C  AND AT EACH CALL OF THE PHYSICS

      print*,'aerosol/gas/cloud properties'

      CALL OPTCI(ycomp,qaer,nmicro,IPRINT)        ! #1
      print*,'On sort de optci'

C  NOW, THIS COMPUTATION IS DONE FOR EACH VALUE OF klon
C  INFRARED. AND THEN IN THE VISIBLE.

        CALL OPTCV(qaer,nmicro,IPRINT)        ! #2

        do j=1,NLAYER
        DZED(j)=Z(J)-Z(J+1)
        enddo
       
c      print*,wlnv
c      print*,""
c      print*,wlni
c      stop

300   CONTINUE  ! fin notfirstcall
       
       
c -----------------------------
c on ne recalcule pas optci si microfi=0 et compo lellouch
c -----------------------------
      IF ((MICROFI.ge.1).or.(.not.ylellouch)) THEN  
      IF(notfirstcall)    THEN  !F au 1er appel   T aux autres appels!!
       print*,'aerosol/gas/cloud properties'
       CALL OPTCI(ycomp,qaer,nmicro,IPRINT)        ! #1
      ENDIF
      ENDIF
      
c ni optcv si microfi=0

      IF (MICROFI.ge.1) THEN  
      IF(notfirstcall)    THEN  !F au 1er appel   T aux autres appels!!
       print*,'aerosol/gas/cloud properties'
       CALL OPTCV(qaer,nmicro,IPRINT)        ! #2
      ENDIF
      ENDIF
     
c -----------------------------
         if (klon.eq.1) then 
           ig=1
         else
           ig=klon/2
         endif
c       print*,"DTAUI(equateur,:,1)=",DTAUI(ig,:,1)
c       print*,"DTAUI(equateur,:,10)=",DTAUI(ig,:,10)
c       print*,"DTAUI(equateur,:,NSPECI)=",DTAUI(ig,:,NSPECI)
c       print*,"DTAUV(equateur,:,1,2)=",DTAUV(ig,:,1,2)
c       print*,"DTAUV(equateur,:,10,2)=",DTAUV(ig,:,10,2)
c       print*,"DTAUV(equateur,:,NSPECV,2)=",DTAUV(ig,:,NSPECV,2)
c       stop

      notfirstcall=.true.

      RETURN
 191  FORMAT(F8.2,1P10E10.2)
 192  FORMAT(a8,1P10E10.2)
      END
