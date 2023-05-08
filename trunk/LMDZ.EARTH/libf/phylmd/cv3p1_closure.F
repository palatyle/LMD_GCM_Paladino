!
! $Id: cv3p1_closure.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE cv3p1_closure(nloc,ncum,nd,icb,inb
     :                      ,pbase,plcl,p,ph,tv,tvp,buoy
     :                      ,Supmax,ok_inhib,Ale,Alp
     o                      ,sig,w0,ptop2,cape,cin,m,iflag,coef
     :                      ,Plim1,Plim2,asupmax,supmax0
     :                      ,asupmaxmin,cbmf)

*
***************************************************************
*                                                             *
* CV3P1_CLOSURE                                               *
*                  Ale & Alp Closure of Convect3              *
*                                                             *
* written by   :   Kerry Emanuel                              *
* vectorization:   S. Bony                                    *
* modified by :    Jean-Yves Grandpeix, 18/06/2003, 19.32.10  *
*                  Julie Frohwirth,     14/10/2005  17.44.22  *
***************************************************************
*
      implicit none

#include "cvthermo.h"
#include "cv3param.h"
#include "YOMCST2.h"
#include "YOMCST.h"
#include "conema3.h"
#include "iniprint.h"

c input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real pbase(nloc),plcl(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real tv(nloc,nd),tvp(nloc,nd), buoy(nloc,nd)
      real Supmax(nloc,nd)
      logical ok_inhib ! enable convection inhibition by dryness
      real Ale(nloc),Alp(nloc)

c input/output:
      real sig(nloc,nd), w0(nloc,nd), ptop2(nloc)

c output:
      real cape(nloc),cin(nloc)
      real m(nloc,nd)
      real Plim1(nloc),Plim2(nloc)
      real asupmax(nloc,nd),supmax0(nloc)
      real asupmaxmin(nloc)
      integer iflag(nloc)
c
c local variables:
      integer il, i, j, k, icbmax, i0(nloc)
      real deltap, fac, w, amu
      real rhodp
      real Pbmxup
      real dtmin(nloc,nd), sigold(nloc,nd)
      real coefmix(nloc,nd)
      real pzero(nloc),ptop2old(nloc)
      real cina(nloc),cinb(nloc)
      integer ibeg(nloc)
      integer nsupmax(nloc)
      real supcrit,temp(nloc,nd)
      real P1(nloc),Pmin(nloc),plfc(nloc)
      real asupmax0(nloc)
      logical ok(nloc)
      real siglim(nloc,nd),wlim(nloc,nd),mlim(nloc,nd)
      real wb2(nloc)
      real cbmflim(nloc),cbmf1(nloc),cbmfmax(nloc),cbmf(nloc)
      real cbmflast(nloc)
      real coef(nloc)
      real xp(nloc),xq(nloc),xr(nloc),discr(nloc),b3(nloc),b4(nloc)
      real theta(nloc),bb(nloc)
      real term1,term2,term3
      real alp2(nloc) ! Alp with offset
      real wb,sigmax
      data wb /2./, sigmax /0.1/

      CHARACTER (LEN=20) :: modname='cv3p1_closure'
      CHARACTER (LEN=80) :: abort_message
c
c      print *,' -> cv3p1_closure, Ale ',ale(1)
c

c -------------------------------------------------------
c -- Initialization
c -------------------------------------------------------

c
c
      do il = 1,ncum
       alp2(il) = max(alp(il),1.e-5)
cIM 
       alp2(il) = max(alp(il),1.e-12)
      enddo
c
      PBMXUP=50.    ! PBMXUP+PBCRIT = cloud depth above which mixed updraughts
c                     exist (if any)

       if(prt_level.GE.20)
     . print*,'cv3p1_param nloc ncum nd icb inb nl',nloc,ncum,nd,
     . icb(nloc),inb(nloc),nl
      do k=1,nl
       do il=1,ncum
        m(il,k)=0.0
       enddo
      enddo

c -------------------------------------------------------
c -- Reset sig(i) and w0(i) for i>inb and i<icb
c -------------------------------------------------------

c update sig and w0 above LNB:

      do 100 k=1,nl-1
       do 110 il=1,ncum
        if ((inb(il).lt.(nl-1)).and.(k.ge.(inb(il)+1)))then
         sig(il,k)=beta*sig(il,k)
     :            +2.*alpha*buoy(il,inb(il))*ABS(buoy(il,inb(il)))
         sig(il,k)=AMAX1(sig(il,k),0.0)
         w0(il,k)=beta*w0(il,k)
        endif
 110   continue
 100  continue

c      if(prt.level.GE.20) print*,'cv3p1_param apres 100'
c compute icbmax:

      icbmax=2
      do 200 il=1,ncum
        icbmax=MAX(icbmax,icb(il))
 200  continue
!     if(prt.level.GE.20) print*,'cv3p1_param apres 200'

c update sig and w0 below cloud base:

      do 300 k=1,icbmax
       do 310 il=1,ncum
        if (k.le.icb(il))then
         sig(il,k)=beta*sig(il,k)-2.*alpha*buoy(il,icb(il))
     $                                    *buoy(il,icb(il))
         sig(il,k)=amax1(sig(il,k),0.0)
         w0(il,k)=beta*w0(il,k)
        endif
310    continue
300    continue
       if(prt_level.GE.20) print*,'cv3p1_param apres 300'
c -------------------------------------------------------------
c -- Reset fractional areas of updrafts and w0 at initial time
c -- and after 10 time steps of no convection
c -------------------------------------------------------------

      do 400 k=1,nl-1
       do 410 il=1,ncum
        if (sig(il,nd).lt.1.5.or.sig(il,nd).gt.12.0)then
         sig(il,k)=0.0
         w0(il,k)=0.0
        endif
 410   continue
 400  continue
      if(prt_level.GE.20) print*,'cv3p1_param apres 400'
c
c -------------------------------------------------------------
Cjyg1
C --  Calculate adiabatic ascent top pressure (ptop)
c -------------------------------------------------------------
C
c
cc 1. Start at first level where precipitations form
      do il = 1,ncum
        Pzero(il) = Plcl(il)-PBcrit
      enddo
c
cc 2. Add offset
      do il = 1,ncum
        Pzero(il) = Pzero(il)-PBmxup
      enddo
      do il=1,ncum
         ptop2old(il)=ptop2(il)
      enddo
c
      do il = 1,ncum
cCR:c est quoi ce 300??
        P1(il) = Pzero(il)-300.
      enddo

c    compute asupmax=abs(supmax) up to lnm+1

      DO il=1,ncum
        ok(il)=.true.
        nsupmax(il)=inb(il)
      ENDDO

      DO i = 1,nl
        DO il = 1,ncum
        IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
        IF (P(il,i) .LE. Pzero(il) .and.
     $       supmax(il,i) .lt. 0 .and. ok(il)) THEN
           nsupmax(il)=i
           ok(il)=.false.
        ENDIF    ! end IF (P(i) ...
        ENDIF    ! end IF (icb+1 le i le inb)
        ENDDO
      ENDDO

      if(prt_level.GE.20) print*,'cv3p1_param apres 2.'
      DO i = 1,nl
        DO il = 1,ncum
        asupmax(il,i)=abs(supmax(il,i))
        ENDDO
      ENDDO

c
        DO il = 1,ncum
        asupmaxmin(il)=10.
        Pmin(il)=100.
!IM ??
        asupmax0(il)=0.
        ENDDO

cc 3.  Compute in which level is Pzero

cIM bug      i0 = 18 
       DO il = 1,ncum
        i0(il) = nl
       ENDDO

       DO i = 1,nl
        DO il = 1,ncum
         IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
           IF (P(il,i) .LE. Pzero(il) .AND. P(il,i) .GE. P1(il)) THEN
            IF (Pzero(il) .GT. P(il,i) .AND.
     $           Pzero(il) .LT. P(il,i-1)) THEN
             i0(il) = i
            ENDIF
           ENDIF
          ENDIF
        ENDDO
       ENDDO
       if(prt_level.GE.20) print*,'cv3p1_param apres 3.'

cc 4.  Compute asupmax at Pzero

       DO i = 1,nl
        DO il = 1,ncum
         IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
           IF (P(il,i) .LE. Pzero(il) .AND. P(il,i) .GE. P1(il)) THEN
             asupmax0(il) = 
     $             ((Pzero(il)-P(il,i0(il)-1))*asupmax(il,i0(il))
     $             -(Pzero(il)-P(il,i0(il)))*asupmax(il,i0(il)-1))
     $             /(P(il,i0(il))-P(il,i0(il)-1))
           ENDIF
         ENDIF
        ENDDO
       ENDDO


      DO i = 1,nl
        DO il = 1,ncum
         IF (P(il,i) .EQ. Pzero(il)) THEN
           asupmax(i,il) = asupmax0(il)
         ENDIF
        ENDDO
      ENDDO
      if(prt_level.GE.20) print*,'cv3p1_param apres 4.'

cc 5. Compute asupmaxmin, minimum of asupmax

      DO i = 1,nl
        DO il = 1,ncum
        IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
        IF (P(il,i) .LE. Pzero(il) .AND. P(il,i) .GE. P1(il)) THEN
          IF (asupmax(il,i) .LT. asupmaxmin(il)) THEN
            asupmaxmin(il)=asupmax(il,i)
            Pmin(il)=P(il,i)
          ENDIF
        ENDIF
        ENDIF
        ENDDO
      ENDDO

      DO il = 1,ncum
!IM
        if(prt_level.GE.20) THEN
         print*,'cv3p1_closure il asupmax0 asupmaxmin',il,asupmax0(il),
     $ asupmaxmin(il) ,Pzero(il),Pmin(il)
        endif
          IF (asupmax0(il) .LT. asupmaxmin(il)) THEN
             asupmaxmin(il) = asupmax0(il)
             Pmin(il) = Pzero(il)
          ENDIF
      ENDDO
      if(prt_level.GE.20) print*,'cv3p1_param apres 5.' 

c
c   Compute Supmax at Pzero
c
      DO i = 1,nl
        DO il = 1,ncum
        IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
        IF (P(il,i) .LE. Pzero(il)) THEN
         Supmax0(il) = ((P(il,i  )-Pzero(il))*aSupmax(il,i-1)
     $             -(P(il,i-1)-Pzero(il))*aSupmax(il,i  ))
     $             /(P(il,i)-P(il,i-1))
         GO TO 425
        ENDIF    ! end IF (P(i) ...
        ENDIF    ! end IF (icb+1 le i le inb)
        ENDDO
      ENDDO

425   continue
      if(prt_level.GE.20) print*,'cv3p1_param apres 425.'

cc 6. Calculate ptop2
c
      DO il = 1,ncum
        IF (asupmaxmin(il) .LT. Supcrit1) THEN
          Ptop2(il) = Pmin(il)
        ENDIF

        IF (asupmaxmin(il) .GT. Supcrit1
     $ .AND. asupmaxmin(il) .LT. Supcrit2) THEN
          Ptop2(il) = Ptop2old(il)
        ENDIF

        IF (asupmaxmin(il) .GT. Supcrit2) THEN
            Ptop2(il) =  Ph(il,inb(il))
        ENDIF
      ENDDO
c
      if(prt_level.GE.20) print*,'cv3p1_param apres 6.'

cc 7. Compute multiplying factor for adiabatic updraught mass flux
c
c
      IF (ok_inhib) THEN
c
      DO i = 1,nl
        DO il = 1,ncum
         IF (i .le. nl) THEN
         coefmix(il,i) = (min(ptop2(il),ph(il,i))-ph(il,i))
     $                  /(ph(il,i+1)-ph(il,i))
         coefmix(il,i) = min(coefmix(il,i),1.)
         ENDIF
        ENDDO
      ENDDO
c
c
      ELSE   ! when inhibition is not taken into account, coefmix=1
c

c
      DO i = 1,nl
        DO il = 1,ncum
         IF (i .le. nl) THEN
         coefmix(il,i) = 1.
         ENDIF
        ENDDO
      ENDDO
c
      ENDIF  ! ok_inhib
      if(prt_level.GE.20) print*,'cv3p1_param apres 7.'
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c

Cjyg2
C
c==========================================================================
C
c
c -------------------------------------------------------------
c -- Calculate convective inhibition (CIN)
c -------------------------------------------------------------

c      do i=1,nloc
c      print*,'avant cine p',pbase(i),plcl(i)
c      enddo
c     do j=1,nd
c     do i=1,nloc
c      print*,'avant cine t',tv(i),tvp(i)
c     enddo
c     enddo
      CALL cv3_cine (nloc,ncum,nd,icb,inb
     :                      ,pbase,plcl,p,ph,tv,tvp
     :                      ,cina,cinb,plfc)
c
      DO il = 1,ncum
        cin(il) = cina(il)+cinb(il)
      ENDDO
      if(prt_level.GE.20) print*,'cv3p1_param apres cv3_cine'
c -------------------------------------------------------------
c --Update buoyancies to account for Ale
c -------------------------------------------------------------
c
      CALL cv3_buoy (nloc,ncum,nd,icb,inb
     :                      ,pbase,plcl,p,ph,Ale,Cin
     :                      ,tv,tvp
     :                      ,buoy )
      if(prt_level.GE.20) print*,'cv3p1_param apres cv3_buoy'

c -------------------------------------------------------------
c -- Calculate convective available potential energy (cape),
c -- vertical velocity (w), fractional area covered by
c -- undilute updraft (sig), and updraft mass flux (m)
c -------------------------------------------------------------

      do 500 il=1,ncum
       cape(il)=0.0
 500  continue

c compute dtmin (minimum buoyancy between ICB and given level k):

      do k=1,nl
       do il=1,ncum
         dtmin(il,k)=100.0
       enddo
      enddo

      do 550 k=1,nl
       do 560 j=minorig,nl
        do 570 il=1,ncum
          if ( (k.ge.(icb(il)+1)).and.(k.le.inb(il)).and.
     :         (j.ge.icb(il)).and.(j.le.(k-1)) )then
           dtmin(il,k)=AMIN1(dtmin(il,k),buoy(il,j))
          endif
 570     continue
 560   continue
 550  continue

c the interval on which cape is computed starts at pbase :

      do 600 k=1,nl
       do 610 il=1,ncum

        if ((k.ge.(icb(il)+1)).and.(k.le.inb(il))) then

         deltap = MIN(pbase(il),ph(il,k-1))-MIN(pbase(il),ph(il,k))
         cape(il)=cape(il)+rrd*buoy(il,k-1)*deltap/p(il,k-1)
         cape(il)=AMAX1(0.0,cape(il))
         sigold(il,k)=sig(il,k)


cjyg       Coefficient coefmix limits convection to levels where a sufficient
c          fraction of mixed draughts are ascending.
         siglim(il,k)=coefmix(il,k)*alpha1*dtmin(il,k)*ABS(dtmin(il,k))
         siglim(il,k)=amax1(siglim(il,k),0.0)
         siglim(il,k)=amin1(siglim(il,k),0.01)
cc         fac=AMIN1(((dtcrit-dtmin(il,k))/dtcrit),1.0)
         fac = 1.
         wlim(il,k)=fac*SQRT(cape(il))
         amu=siglim(il,k)*wlim(il,k)
         rhodp = 0.007*p(il,k)*(ph(il,k)-ph(il,k+1))/tv(il,k)
         mlim(il,k)=amu*rhodp
c         print*, 'siglim ', k,siglim(1,k)
        endif

 610   continue
 600  continue
      if(prt_level.GE.20) print*,'cv3p1_param apres 600'

      do 700 il=1,ncum
!IM beg
        if(prt_level.GE.20) THEN
         print*,'cv3p1_closure il icb mlim ph ph+1 ph+2',il,
     $icb(il),mlim(il,icb(il)+1),ph(il,icb(il)),
     $ph(il,icb(il)+1),ph(il,icb(il)+2)
        endif

        if (icb(il)+1.le.inb(il)) then
!IM end
       mlim(il,icb(il))=0.5*mlim(il,icb(il)+1)
     :             *(ph(il,icb(il))-ph(il,icb(il)+1))
     :             /(ph(il,icb(il)+1)-ph(il,icb(il)+2))
!IM beg
        endif !(icb(il.le.inb(il))) then
!IM end
 700  continue
      if(prt_level.GE.20) print*,'cv3p1_param apres 700'

cjyg1
c------------------------------------------------------------------------
cc     Correct mass fluxes so that power used to overcome CIN does not
cc     exceed Power Available for Lifting (PAL).
c------------------------------------------------------------------------
c
      do il = 1,ncum
       cbmflim(il) = 0.
       cbmf(il) = 0.
      enddo
c
cc 1. Compute cloud base mass flux of elementary system (Cbmf0=Cbmflim)
c
      do k= 1,nl
       do il = 1,ncum
!old       IF (k .ge. icb(il) .and. k .le. inb(il)) THEN
!IM        IF (k .ge. icb(il)+1 .and. k .le. inb(il)) THEN
       IF (k .ge. icb(il) .and. k .le. inb(il)         !cor jyg
     $     .and. icb(il)+1 .le. inb(il)) THEN          !cor jyg
         cbmflim(il) = cbmflim(il)+MLIM(il,k)
        ENDIF
       enddo
      enddo
      if(prt_level.GE.20) print*,'cv3p1_param apres cbmflim'

cc 1.5 Compute cloud base mass flux given by Alp closure (Cbmf1), maximum
cc     allowed mass flux (Cbmfmax) and final target mass flux (Cbmf)
cc     Cbmf is set to zero if Cbmflim (the mass flux of elementary cloud) is
c--    exceedingly small.
c
      DO il = 1,ncum
        wb2(il) = sqrt(2.*max(Ale(il)+cin(il),0.))
      ENDDO
c
      DO il = 1,ncum
       cbmf1(il) = alp2(il)/(0.5*wb*wb-Cin(il))
       if(cbmf1(il).EQ.0.AND.alp2(il).NE.0.) THEN
        write(lunout,*)
     &  'cv3p1_closure cbmf1=0 and alp NE 0 il alp2 alp cin ',il,
     . alp2(il),alp(il),cin(il)
        abort_message = ''
        CALL abort_gcm (modname,abort_message,1)
       endif
       cbmfmax(il) = sigmax*wb2(il)*100.*p(il,icb(il))
     :              /(rrd*tv(il,icb(il)))
      ENDDO
c
      DO il = 1,ncum
       IF (cbmflim(il) .gt. 1.e-6) THEN
cATTENTION TEST CR
c         if (cbmfmax(il).lt.1.e-12) then
        cbmf(il) = min(cbmf1(il),cbmfmax(il))
c         else
c         cbmf(il) = cbmf1(il)
c         endif
c        print*,'cbmf',cbmf1(il),cbmfmax(il)
       ENDIF
      ENDDO
      if(prt_level.GE.20) print*,'cv3p1_param apres cbmflim_testCR'
c
cc 2. Compute coefficient and apply correction
c
      do il = 1,ncum
       coef(il) = (cbmf(il)+1.e-10)/(cbmflim(il)+1.e-10)
      enddo
      if(prt_level.GE.20) print*,'cv3p1_param apres coef_plantePLUS'
c
      DO k = 1,nl
        do il = 1,ncum
         IF ( k .ge. icb(il)+1 .AND. k .le. inb(il)) THEN
         amu=beta*sig(il,k)*w0(il,k)+
     :   (1.-beta)*coef(il)*siglim(il,k)*wlim(il,k)
         w0(il,k) = wlim(il,k)
         w0(il,k) =max(w0(il,k),1.e-10)
         sig(il,k)=amu/w0(il,k)
         sig(il,k)=min(sig(il,k),1.)
cc         amu = 0.5*(SIG(il,k)+sigold(il,k))*W0(il,k)
         M(il,k)=AMU*0.007*P(il,k)*(PH(il,k)-PH(il,k+1))/TV(il,k)
         ENDIF
        enddo
      ENDDO
cjyg2
      DO il = 1,ncum
       w0(il,icb(il))=0.5*w0(il,icb(il)+1)
       m(il,icb(il))=0.5*m(il,icb(il)+1)
     $       *(ph(il,icb(il))-ph(il,icb(il)+1))
     $       /(ph(il,icb(il)+1)-ph(il,icb(il)+2))
       sig(il,icb(il))=sig(il,icb(il)+1)
       sig(il,icb(il)-1)=sig(il,icb(il))
      ENDDO
      if(prt_level.GE.20) print*,'cv3p1_param apres w0_sig_M'
c
cc 3. Compute final cloud base mass flux and set iflag to 3 if
cc    cloud base mass flux is exceedingly small and is decreasing (i.e. if
cc    the final mass flux (cbmflast) is greater than the target mass flux
cc    (cbmf)).
c
      do il = 1,ncum
       cbmflast(il) = 0.
      enddo
c
      do k= 1,nl
       do il = 1,ncum
        IF (k .ge. icb(il) .and. k .le. inb(il)) THEN
 !IMpropo??      IF ((k.ge.(icb(il)+1)).and.(k.le.inb(il))) THEN
         cbmflast(il) = cbmflast(il)+M(il,k)
        ENDIF
       enddo
      enddo
c
      do il = 1,ncum
       IF (cbmflast(il) .lt. 1.e-6 .and.
     $     cbmflast(il) .ge. cbmf(il)) THEN
         iflag(il) = 3
       ENDIF
      enddo
c
      do k= 1,nl
       do il = 1,ncum
        IF (iflag(il) .ge. 3) THEN
         M(il,k) = 0.
         sig(il,k) = 0.
         w0(il,k) = 0.
        ENDIF
       enddo
      enddo
      if(prt_level.GE.20) print*,'cv3p1_param apres iflag'
c
cc 4. Introduce a correcting factor for coef, in order to obtain an effective
cc    sigdz larger in the present case (using cv3p1_closure) than in the old
cc    closure (using cv3_closure).
      if (1.eq.0) then
       do il = 1,ncum 
cc      coef(il) = 2.*coef(il)
        coef(il) = 5.*coef(il)
       enddo
c version CVS du ..2008
      else
       if (iflag_cvl_sigd.eq.0) then
ctest pour verifier qu on fait la meme chose qu avant: sid constant
        coef(1:ncum)=1.
       else
        coef(1:ncum) = min(2.*coef(1:ncum),5.)
        coef(1:ncum) = max(2.*coef(1:ncum),0.2)
       endif
      endif
c
      if(prt_level.GE.20) print*,'cv3p1_param FIN'
       return
       end


