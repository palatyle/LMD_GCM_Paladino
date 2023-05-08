!
! $Id: cv3_cine.F 1403 2010-07-01 09:02:53Z fairhead $
!
        SUBROUTINE cv3_cine(nloc,ncum,nd,icb,inb
     :                      ,pbase,plcl,p,ph,tv,tvp
     :                      ,cina,cinb,plfc)

***************************************************************
*                                                             *
* CV3_CINE                                                    *
*                                                             *
*                                                             *
* written by   :   Frederique Cheruy                          *
* vectorization:   Jean-Yves Grandpeix, 19/06/2003, 11.54.43  *
* modified by :                                               *
***************************************************************
*
      implicit none
c
#include "YOMCST.h"
#include "cvthermo.h"
#include "cv3param.h"
c input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real pbase(nloc),plcl(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real tv(nloc,nd),tvp(nloc,nd)
c
c output
      real cina(nloc),cinb(nloc),plfc(nloc)
c
c local variables
      integer il,i,j,k
      integer itop(nloc),ineg(nloc),ilow(nloc)
      integer ifst(nloc),isublcl(nloc)
      logical lswitch(nloc),lswitch1(nloc),lswitch2(nloc)
      logical exist_lfc(nloc)
      real dpmax
      real deltap,dcin
      real buoylcl(nloc),tvplcl(nloc),tvlcl(nloc)
      real p0(nloc)
      real buoyz(nloc), buoy(nloc,nd)
c
c-------------------------------------------------------------
c     Initialization
c-------------------------------------------------------------
      do il = 1,ncum
       cina(il) = 0.
       cinb(il) = 0.
      enddo
c
c--------------------------------------------------------------
c      Recompute buoyancies
c--------------------------------------------------------------
      DO k = 1,nd
        DO il = 1,ncum
!      print*,'tvp tv=',tvp(il,k),tv(il,k)
          buoy(il,k) = tvp(il,k) - tv(il,k)
        ENDDO
      ENDDO
c---------------------------------------------------------------
c
c   calcul de la flottabilite a LCL (Buoylcl)
c     ifst = first P-level above lcl
c     isublcl = highest P-level below lcl.
c---------------------------------------------------------------
c
      do il = 1,ncum
       TVPlcl(il) = TVP(il,1)*(Plcl(il)/P(il,1))**(2./7.) !For dry air, R/Cp=2/7
      enddo
c
      do il = 1,ncum
       IF (Plcl(il) .GT. P(il,icb(il))) THEN
        ifst(il) = icb(il)
        isublcl(il) = icb(il)-1
       ELSE
        ifst(il) = icb(il)+1
        isublcl(il) = icb(il)
       ENDIF
      enddo
c
      do il = 1,ncum
       TVlcl(il)=TV(il,ifst(il)-1)+(TV(il,ifst(il))-TV(il,ifst(il)-1))
     $   *(Plcl(il)-P(il,ifst(il)-1))/(P(il,ifst(il))-P(il,ifst(il)-1))
      enddo
c
      do il = 1,ncum
        BUOYlcl(il) = TVPlcl(il)-TVlcl(il)
      enddo
c
c---------------------------------------------------------------
c premiere couche contenant un  niveau de flotabilite positive
c et premiere couche contenant un  niveau de flotabilite negative
c  au dessus du niveau de condensation
c---------------------------------------------------------------
      do il = 1,ncum
        itop(il) =nl-1
        ineg(il) = nl-1
        exist_lfc(il) = .FALSE.
      enddo
      do 100 k=nl-1,1,-1
       do 110 il=1,ncum
        if (k .ge. ifst(il)) then
         if (buoy(il,k) .gt. 0.) then
          itop(il)=k
          exist_lfc(il) = .TRUE.
         else
          ineg(il)=k
         endif
        endif
110    continue
100   continue
c
c---------------------------------------------------------------
c When there is no positive buoyancy level, set Plfc, Cina and Cinb
c to arbitrary extreme values.
c---------------------------------------------------------------
      DO il = 1,ncum
       IF (.NOT.exist_lfc(il)) THEN
         Plfc(il) = 1.111
         Cinb(il) = -1111.
         Cina(il) = -1112.
       ENDIF
      ENDDO
c
c
c---------------------------------------------------------------
c -- Two cases : BUOYlcl >= 0 and BUOYlcl < 0.
c---------------------------------------------------------------
C
C--------------------
C -- 1.0 BUOYlcl >=0.
C--------------------
c
      DPMAX = 50.
      DO il = 1,ncum
        lswitch1(il)=BUOYlcl(il) .GE. 0. .AND. exist_lfc(il)
        lswitch(il) = lswitch1(il)
      ENDDO
c
c 1.1 No inhibition case
c ----------------------
C   If buoyancy is positive at LCL and stays positive over a large enough
C pressure interval (=DPMAX), inhibition is set to zero,
C
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        IF (P(il,ineg(il)) .LT. P(il,icb(il))-DPmax) THEN
          PLFC(il) = Plcl(il)
          Cina(il) = 0.
          Cinb(il) = 0.
        ENDIF
      ENDIF
      ENDDO
c
c 1.2 Upper inhibition only case
c ------------------------------
      DO il = 1,ncum
        lswitch2(il)= P(il,ineg(il)) .GE. P(il,icb(il))-DPmax
        lswitch(il) = lswitch1(il) .AND. lswitch2(il)
      ENDDO
c
      DO il = 1,ncum
      IF (lswitch(il)) THEN
          Cinb(il) = 0.
c
c 1.2.1  Calcul de la pression du niveau de flot. nulle juste au-dessus de LCL
c ---------------------------------------------------------------------------
         IF (ineg(il) .GT. isublcl(il)+1) THEN
C In order to get P0, one may interpolate linearly buoyancies
C  between P(ineg) and P(ineg-1).
        P0(il)=(buoy(il,ineg(il))*P(il,ineg(il)-1)
     $         -buoy(il,ineg(il)-1)*P(il,ineg(il)))
     :           / (buoy(il,ineg(il))-buoy(il,ineg(il)-1))
         ELSE
C In order to get P0, one has to interpolate between P(ineg) and Plcl.
        P0(il) = (BUOY(il,ineg(il))*Plcl(il)-BUOYlcl(il)*P(il,ineg(il)))
     $          /(BUOY(il,ineg(il))     -BUOYlcl(il))
         ENDIF
      ENDIF
      ENDDO
c
c 1.2.2 Recompute itop (=1st layer with positive buoyancy above ineg)
c -------------------------------------------------------------------
      do il = 1,ncum
      IF (lswitch(il)) THEN
        itop(il) =nl-1
      ENDIF
      enddo
c
      do  k=nl,1,-1
       do  il=1,ncum
       IF (lswitch(il)) THEN
        if (k .ge. ineg(il) .and. buoy(il,k) .gt. 0) then
         itop(il)=k
        endif
       ENDIF
       enddo
      enddo
c
c 1.2.3 Computation of PLFC
c -------------------------
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        PLFC(il)=(buoy(il,itop(il))*P(il,itop(il)-1)
     $           -buoy(il,itop(il)-1)*P(il,itop(il)))
     $           / (buoy(il,itop(il))-buoy(il,itop(il)-1))
      ENDIF
      ENDDO
c
c 1.2.4 Computation of CINA
c -------------------------
c
C   Upper part of CINA : integral from P(itop-1) to Plfc
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P(il,itop(il)-1)-Plfc(il)
        dcin = RD*BUOY(il,itop(il)-1)*deltap
     $        / (P(il,itop(il)-1)+Plfc(il))
        CINA(il) = min(0.,dcin)
      ENDIF
      ENDDO
c
C   Middle part of CINA : integral from P(ineg) to P(itop-1)
      DO k = 1,nl
        DO il = 1,ncum
        IF (lswitch(il)) THEN
          IF (k .GE. ineg(il) .AND. k .LE. itop(il)-2) THEN
           deltap = P(il,k)-P(il,k+1)
           dcin = 0.5*RD*(BUOY(il,k)+BUOY(il,k+1))*deltap/PH(il,k+1)
           CINA(il) = CINA(il) + min(0.,dcin)
          ENDIF
        ENDIF
        ENDDO
      ENDDO
c
C   Lower part of CINA : integral from P0 to P(ineg)
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P0(il)-P(il,ineg(il))
        dcin = RD*BUOY(il,ineg(il))*deltap/(P(il,ineg(il))+P0(il))
        CINA(il) = CINA(il) + min(0.,dcin)
      ENDIF
      ENDDO
c
C
C ------------------
C -- 2.0 BUOYlcl <0.
C ------------------
C
      DO il = 1,ncum
        lswitch1(il)=BUOYlcl(il) .LT. 0. .AND. exist_lfc(il)
        lswitch(il) = lswitch1(il)
      ENDDO
c
c 2.0.1 Premiere  couche ou la flotabilite est negative au dessus du sol
c ----------------------------------------------------
c    au cas ou elle existe  sinon ilow=1 (nk apres)
c      on suppose que la parcelle part de la premiere couche
c
      DO il = 1,ncum
      IF (lswitch(il)) THEN
       ilow(il)=1
      ENDIF
      ENDDO
c
      do 200 k=nl,1,-1
        DO il = 1,ncum
        IF (lswitch(il) .AND. k .LE.icb(il)-1) THEN
         if(buoy(il,k).lt. 0.) then
           ilow(il) = k
          endif
        ENDIF
        ENDDO
 200  continue

c 2.0.2  Calcul de la pression du niveau de flot. nulle sous le nuage
c ----------------------------------------------------
      DO il = 1,ncum
      IF (lswitch(il)) THEN
       if(ilow(il).gt. 1) then
         P0(il)=(buoy(il,ilow(il))*P(il,ilow(il)-1)
     $          -buoy(il,ilow(il)-1)*P(il,ilow(il)))
     :            / (buoy(il,ilow(il))-buoy(il,ilow(il)-1))
         BUOYz(il) = 0.
       else
         P0(il) = P(il,1)
         BUOYz(il) = BUOY(il,1)
       endif
      ENDIF
      ENDDO
c
C 2.1. Computation of CINB
C -----------------------
c
      DO il = 1,ncum
        lswitch2(il)= (isublcl(il) .EQ. 1 .AND. ilow(il) .EQ. 1)
     $                  .OR.(isublcl(il) .EQ. ilow(il)-1)
        lswitch(il) = lswitch1(il) .AND. lswitch2(il)
      ENDDO
cc      IF (    (isublcl .EQ. 1 .AND. ilow .EQ. 1)
cc     $    .OR.(isublcl .EQ. ilow-1)) THEN
c
c 2.1.1 First case : Plcl just above P0
c -------------------------------------
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P0(il)-Plcl(il)
        dcin = RD*(BUOYz(il)+BUOYlcl(il))*deltap/(P0(il)+Plcl(il))
        CINB(il) = min(0.,dcin)
      ENDIF
      ENDDO
c
      DO il = 1,ncum
        lswitch(il) = lswitch1(il) .AND. .NOT. lswitch2(il)
      ENDDO
cc      ELSE
c
c 2.1.2 Second case : there is at least one P-level between P0 and Plcl
c ---------------------------------------------------------------------
c
C   Lower part of CINB : integral from P0 to P(ilow)
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P0(il)-P(il,ilow(il))
        dcin = RD*(BUOYz(il)+BUOY(il,ilow(il)))*deltap
     $         /(P0(il)+P(il,ilow(il)))
        CINB(il) = min(0.,dcin)
      ENDIF
      ENDDO
c
c
C  Middle part of CINB : integral from P(ilow) to P(isublcl)
cc      DO k = ilow,isublcl-1
      DO k = 1,nl
        DO il = 1,ncum
        IF (lswitch(il)
     $   .AND. k .GE. ilow(il) .AND. k .LE. isublcl(il)-1) THEN
          deltap = P(il,k)-P(il,k+1)
          dcin = 0.5*RD*(BUOY(il,k)+BUOY(il,k+1))*deltap/PH(il,k+1)
          CINB(il) = CINB(il) + min(0.,dcin)
        ENDIF
        ENDDO
      ENDDO
c
C  Upper part of CINB : integral from P(isublcl) to Plcl
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P(il,isublcl(il)) - Plcl(il)
        dcin = RD*(BUOY(il,isublcl(il))+BUOYlcl(il))*deltap
     $         /(P(il,isublcl(il))+Plcl(il))
        CINB(il) = CINB(il)+min(0.,dcin)
      ENDIF
      ENDDO
C
c
cc      ENDIF
c
C 2.2 Computation of CINA
c ---------------------
c
      DO il = 1,ncum
        lswitch2(il)= Plcl(il) .GT. P(il,itop(il)-1)
        lswitch(il) = lswitch1(il) .AND. lswitch2(il)
      ENDDO
c
c 2.2.1 FIrst case : Plcl > P(itop-1)
C ---------------------------------
C In order to get Plfc, one may interpolate linearly buoyancies
C  between P(itop) and P(itop-1).
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        PLFC(il)=(buoy(il,itop(il))*P(il,itop(il)-1)
     $           -buoy(il,itop(il)-1)*P(il,itop(il)))
     $           / (buoy(il,itop(il))-buoy(il,itop(il)-1))
      ENDIF
      ENDDO
c
C   Upper part of CINA : integral from P(itop-1) to Plfc
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = P(il,itop(il)-1)-Plfc(il)
        dcin = RD*BUOY(il,itop(il)-1)*deltap
     $         /(P(il,itop(il)-1)+Plfc(il))
        CINA(il) = min(0.,dcin)
      ENDIF
      ENDDO
c
C   Middle part of CINA : integral from P(icb+1) to P(itop-1)
      DO k = 1,nl
        DO il = 1,ncum
        IF (lswitch(il)
     $     .AND. k .GE. icb(il)+1 .AND. k .LE. itop(il)-2) THEN
          deltap = P(il,k)-P(il,k+1)
          dcin = 0.5*RD*(BUOY(il,k)+BUOY(il,k+1))*deltap/PH(il,k+1)
          CINA(il) = CINA(il) + min(0.,dcin)
        ENDIF
        ENDDO
      ENDDO
c
C   Lower part of CINA : integral from Plcl to P(icb+1)
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        IF (Plcl(il) .GT. P(il,icb(il))) THEN
          IF (icb(il) .LT. itop(il)-1) THEN
            deltap = P(il,icb(il))-P(il,icb(il)+1)
            dcin = 0.5*RD*(BUOY(il,icb(il))+BUOY(il,icb(il)+1))
     $                   *deltap/PH(il,icb(il)+1)
            CINA(il) = CINA(il)+min(0.,dcin)
          ENDIF
c
          deltap = Plcl(il)-P(il,icb(il))
          dcin = RD*(BUOYlcl(il)+BUOY(il,icb(il)))
     $              *deltap/(Plcl(il)+P(il,icb(il)))
          CINA(il) = CINA(il)+min(0.,dcin)
        ELSE
          deltap = Plcl(il)-P(il,icb(il)+1)
          dcin = RD*(BUOYlcl(il)+BUOY(il,icb(il)+1))
     $             *deltap/(Plcl(il)+P(il,icb(il)+1))
          CINA(il) = CINA(il)+min(0.,dcin)
        ENDIF
      ENDIF
      ENDDO
c
      DO il = 1,ncum
        lswitch(il) = lswitch1(il) .AND. .NOT. lswitch2(il)
      ENDDO
cc      ELSE
c
c 2.2.2 Second case : Plcl lies between P(itop-1) and P(itop);
C ----------------------------------------------------------
C In order to get Plfc, one has to interpolate between P(itop) and Plcl.
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        PLFC(il) =
     $    (BUOY(il,itop(il))*Plcl(il)-BUOYlcl(il)*P(il,itop(il)))
     $          /(BUOY(il,itop(il))     -BUOYlcl(il))
      ENDIF
      ENDDO
c
      DO il = 1,ncum
      IF (lswitch(il)) THEN
        deltap = Plcl(il)-Plfc(il)
        dcin = RD*BUOYlcl(il)*deltap/(Plcl(il)+Plfc(il))
        CINA(il) = min(0.,dcin)
      ENDIF
      ENDDO
cc      ENDIF
c


      RETURN
      END
