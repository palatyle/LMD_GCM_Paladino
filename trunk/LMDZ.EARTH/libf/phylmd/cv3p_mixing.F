      SUBROUTINE cv3p_mixing(nloc,ncum,nd,na,ntra,icb,nk,inb
     :                    ,ph,t,rr,rs,u,v,tra,h,lv,qnk
     :                    ,unk,vnk,hp,tv,tvp,ep,clw,sig
     :                    ,ment,qent,hent,uent,vent,nent
     :                   ,sij,elij,supmax,ments,qents,traent)
***************************************************************
*                                                             *
* CV3P_MIXING : compute mixed draught properties and,         *
*               within a scaling factor, mixed draught        *
*               mass fluxes.                                  *
* written by  : VTJ Philips,JY Grandpeix, 21/05/2003, 09.14.15*
* modified by :                                               *
***************************************************************
*
      implicit none
c
#include "cvthermo.h"
#include "cv3param.h"
#include "YOMCST2.h"

c inputs:
      integer ncum, nd, na, ntra, nloc
      integer icb(nloc), inb(nloc), nk(nloc)
      real sig(nloc,nd)
      real qnk(nloc),unk(nloc),vnk(nloc)
      real ph(nloc,nd+1)
      real t(nloc,nd), rr(nloc,nd), rs(nloc,nd)
      real u(nloc,nd), v(nloc,nd)
      real tra(nloc,nd,ntra) ! input of convect3
      real lv(nloc,na)
      real h(nloc,na)  !liquid water static energy of environment
      real hp(nloc,na) !liquid water static energy of air shed from adiab. asc.
      real tv(nloc,na), tvp(nloc,na), ep(nloc,na), clw(nloc,na)

c outputs:
      real ment(nloc,na,na), qent(nloc,na,na)
      real uent(nloc,na,na), vent(nloc,na,na)
      real sij(nloc,na,na), elij(nloc,na,na)
      real supmax(nloc,na)     ! Highest mixing fraction of mixed updraughts
                               ! with the sign of (h-hp)
      real traent(nloc,nd,nd,ntra)
      real ments(nloc,nd,nd), qents(nloc,nd,nd)
      real sigij(nloc,nd,nd)
      real hent(nloc,nd,nd)
      integer nent(nloc,nd)

c local variables:
      integer i, j, k, il, im, jm
      integer num1, num2
      real rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
      real alt, delp, delm
      real Qmixmax(nloc), Rmixmax(nloc), SQmRmax(nloc)
      real Qmixmin(nloc), Rmixmin(nloc), SQmRmin(nloc)
      real signhpmh(nloc)
      real Sx, Scrit2
      integer Jx
      real smid(nloc), sjmin(nloc), sjmax(nloc)
      real Sbef(nloc), Sup(nloc), Smin(nloc)
      real asij(nloc), smax(nloc), scrit(nloc)
      real csum(nloc,nd)
      real awat
      logical lwork(nloc)
c
      REAL amxupcrit, df, ff
      INTEGER nstep
C
c --   Mixing probability distribution functions
c
      real Qcoef1,Qcoef2,QFF,QFFF,Qmix,Rmix,Qmix1,Rmix1,Qmix2,Rmix2,F
      Qcoef1(F) = tanh(F/gammas)
      Qcoef2(F) = ( tanh(F/gammas) + gammas *
     $            log(cosh((1.- F)/gammas)/cosh(F/gammas)))
      QFF(F) = Max(Min(F,1.),0.)
      QFFF(F) = Min(QFF(F),scut)
      Qmix1(F) = ( tanh((QFF(F) - Fmax)/gammas)+Qcoef1max )/
     $           Qcoef2max
      Rmix1(F) = ( gammas*log(cosh((QFF(F)-Fmax)/gammas))
     1             +QFF(F)*Qcoef1max ) / Qcoef2max
      Qmix2(F) = -Log(1.-QFFF(F))/scut
      Rmix2(F) = (QFFF(F)+(1.-QFF(F))*Log(1.-QFFF(F)))/scut
      Qmix(F) = qqa1*Qmix1(F) + qqa2*Qmix2(F)
      Rmix(F) = qqa1*Rmix1(F) + qqa2*Rmix2(F)
C
      INTEGER ifrst
      DATA ifrst/0/
C

c=====================================================================
c --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
c=====================================================================
c
c -- Initialize mixing PDF coefficients
      IF (ifrst .EQ. 0) THEN
        ifrst = 1
        Qcoef1max = Qcoef1(Fmax)
        Qcoef2max = Qcoef2(Fmax)
c
      ENDIF
c

c ori        do 360 i=1,ncum*nlp
        do 361 j=1,nl
        do 360 i=1,ncum
          nent(i,j)=0
c in convect3, m is computed in cv3_closure
c ori          m(i,1)=0.0
 360    continue
 361    continue

c ori      do 400 k=1,nlp
c ori       do 390 j=1,nlp
      do 400 j=1,nl
       do 390 k=1,nl
          do 385 i=1,ncum
            qent(i,k,j)=rr(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
            elij(i,k,j)=0.0
            hent(i,k,j)=0.0
            ment(i,k,j)=0.0
            sij(i,k,j)=0.0
 385      continue
 390    continue
 400  continue

      do k=1,ntra
       do j=1,nd  ! instead nlp
        do i=1,nd ! instead nlp
         do il=1,ncum
            traent(il,i,j,k)=tra(il,j,k)
         enddo
        enddo
       enddo
      enddo

c=====================================================================
c --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
c --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
c --- FRACTION (sij)
c=====================================================================

      do 750 i=minorig+1, nl

       do 710 j=minorig,nl
        do 700 il=1,ncum
         if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
     :      (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then

          rti=qnk(il)-ep(il,i)*clw(il,i)
          bf2=1.+lv(il,j)*lv(il,j)*rs(il,j)/(rrv*t(il,j)*t(il,j)*cpd)
          anum=h(il,j)-hp(il,i)+(cpv-cpd)*t(il,j)*(rti-rr(il,j))
          denom=h(il,i)-hp(il,i)+(cpd-cpv)*(rr(il,i)-rti)*t(il,j)
          dei=denom
          if(abs(dei).lt.0.01)dei=0.01
          sij(il,i,j)=anum/dei
          sij(il,i,i)=1.0
          altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
          altem=altem/bf2
          cwat=clw(il,j)*(1.-ep(il,j))
          stemp=sij(il,i,j)
          if((stemp.lt.0.0.or.stemp.gt.1.0.or.altem.gt.cwat)
     :                 .and.j.gt.i)then
           anum=anum-lv(il,j)*(rti-rs(il,j)-cwat*bf2)
           denom=denom+lv(il,j)*(rr(il,i)-rti)
           if(abs(denom).lt.0.01)denom=0.01
           sij(il,i,j)=anum/denom
           altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
           altem=altem-(bf2-1.)*cwat
          end if
              if(sij(il,i,j).gt.0.0)then
ccc                 ment(il,i,j)=m(il,i)
                 ment(il,i,j)=1.
                 elij(il,i,j)=altem
                 elij(il,i,j)=amax1(0.0,elij(il,i,j))
                 nent(il,i)=nent(il,i)+1
              endif

         sij(il,i,j)=amax1(0.0,sij(il,i,j))
         sij(il,i,j)=amin1(1.0,sij(il,i,j))
         endif ! new
 700   continue
 710  continue

c
c   ***   if no air can entrain at level i assume that updraft detrains  ***
c   ***   at that level and calculate detrained air flux and properties  ***
c

c@      do 170 i=icb(il),inb(il)

      do 740 il=1,ncum
      if ((i.ge.icb(il)).and.(i.le.inb(il))
     :                  .and.(nent(il,i).eq.0)) then
c@      if(nent(il,i).eq.0)then
ccc      ment(il,i,i)=m(il,i)
      ment(il,i,i)=1.
      qent(il,i,i)=qnk(il)-ep(il,i)*clw(il,i)
      uent(il,i,i)=unk(il)
      vent(il,i,i)=vnk(il)
      elij(il,i,i)=clw(il,i)*(1.-ep(il,i))
      sij(il,i,i)=0.0
      end if
 740  continue
 750  continue

      do j=1,ntra
       do i=minorig+1,nl
        do il=1,ncum
         if (i.ge.icb(il) .and. i.le.inb(il)
     :                    .and. nent(il,i).eq.0) then
          traent(il,i,i,j)=tra(il,nk(il),j)
         endif
        enddo
       enddo
      enddo

      do 100 j=minorig,nl
      do 101 i=minorig,nl
      do 102 il=1,ncum
      if ((j.ge.(icb(il)-1)).and.(j.le.inb(il))
     :    .and.(i.ge.icb(il)).and.(i.le.inb(il)))then
       sigij(il,i,j)=sij(il,i,j)
      endif
 102  continue
 101  continue
 100  continue
c@      enddo

c@170   continue

c=====================================================================
c   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
c   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
c=====================================================================

      call zilch(csum,nloc*nd)

      do il=1,ncum
       lwork(il) = .FALSE.
      enddo

      DO 789 i=minorig+1,nl

      num1=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) num1=num1+1
      enddo
      if (num1.le.0) goto 789


      do 781 il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) then
        lwork(il)=(nent(il,i).ne.0)
        Signhpmh(il) = sign(1.,hp(il,i)-h(il,i))
        qp=qnk(il)-ep(il,i)*clw(il,i)
        anum=h(il,i)-hp(il,i)-lv(il,i)*(qp-rs(il,i))
     :           +(cpv-cpd)*t(il,i)*(qp-rr(il,i))
        denom=h(il,i)-hp(il,i)+lv(il,i)*(rr(il,i)-qp)
     :           +(cpd-cpv)*t(il,i)*(rr(il,i)-qp)
        if(abs(denom).lt.0.01)denom=0.01
        scrit(il)=min(anum/denom,1.)
        alt=qp-rs(il,i)+scrit(il)*(rr(il,i)-qp)
c
cjyg1    Find maximum of SIJ for J>I, if any, and new critical value Scrit2
c        such that : Sij > Scrit2  => mixed draught will detrain at J<I
c                    Sij < Scrit2  => mixed draught will detrain at J>I
c
       Sx = 0.
       Jx = 0.
       Sbef(il) = max(0.,signhpmh(il))
       DO j = i+1,inb(il)
         IF (Sbef(il) .LT. Sij(il,i,j)) THEN
           Sx = max(Sij(il,i,j),Sx)
           Jx = J
         ENDIF
         Sbef(il) = Sij(il,i,j)
       ENDDO
c
       Scrit2 = min(Scrit(il),Sx)*max(0.,-signhpmh(il))
     :         +Scrit(il)*max(0.,signhpmh(il))
c
       Scrit(il) = Scrit2
c
cjyg    Correction pour la nouvelle logique; la correction pour ALT
c       est un peu au hazard
       if(scrit(il).le.0.0)scrit(il)=0.0
       if(alt.le.0.0) scrit(il)=1.0
C
        smax(il)=0.0
        asij(il)=0.0
       Sup(il)=0.     ! upper S-value reached by descending draughts
       endif
781   continue

      do 175 j=minorig,nl

      num2=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) num2=num2+1
      enddo
      if (num2.le.0) goto 175

c -----------------------------------------------
         IF (j .GT. i) THEN
c -----------------------------------------------
      do 782 il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then
       if(sij(il,i,j).gt.0.0)then
           Smid(il)=min(Sij(il,i,j),Scrit(il))
           Sjmax(il)=Smid(il)
           Sjmin(il)=Smid(il)
           IF (Smid(il) .LT. Smin(il) .AND.
     1                         Sij(il,i,j+1) .LT. Smid(il)) THEN
             Smin(il)=Smid(il)
             Sjmax(il)=min( (Sij(il,i,j+1)+Sij(il,i,j))/2. ,
     1                  Sij(il,i,j) ,
     1                  Scrit(il) )
             Sjmin(il)=max( (Sbef(il)+Sij(il,i,j))/2. ,
     1                  Sij(il,i,j) )
             Sjmin(il)=min(Sjmin(il),Scrit(il))
             Sbef(il) = Sij(il,i,j)
           ENDIF
      endif
      endif
782   continue
c -----------------------------------------------
         ELSE IF (j .EQ. i) THEN
c -----------------------------------------------
      do 783 il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then
       if(sij(il,i,j).gt.0.0)then
           Smid(il) = 1.
           Sjmin(il) = max((Sij(il,i,j-1)+Smid(il))/2.,Scrit(il))
     1                                         *max(0.,-signhpmh(il))
     1            +min((Sij(il,i,j+1)+Smid(il))/2.,Scrit(il))
     1                                         *max(0., signhpmh(il))
           Sjmin(il) = max(Sjmin(il),Sup(il))
           Sjmax(il) = 1.
c
c-           preparation des variables Scrit, Smin et Sbef pour la partie j>i
           Scrit(il) = min(Sjmin(il),Sjmax(il),Scrit(il))

           Smin(il) = 1.
           Sbef(il) = max(0.,signhpmh(il))
           Supmax(il,i) = sign(Scrit(il),-signhpmh(il))
      endif
      endif
783   continue
c -----------------------------------------------
         ELSE IF ( j .LT. i) THEN
c -----------------------------------------------
      do 784 il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then
       if(sij(il,i,j).gt.0.0)then
           Smid(il)=max(Sij(il,i,j),Scrit(il))
           Sjmax(il) = Smid(il)
           Sjmin(il) = Smid(il)
           IF (Smid(il) .GT. Smax(il) .AND.
     1                          Sij(il,i,j+1) .GT. Smid(il)) THEN
             Smax(il) = Smid(il)
             Sjmax(il) = max( (Sij(il,i,j+1)+Sij(il,i,j))/2. ,
     1                                               Sij(il,i,j) )
             Sjmax(il) = max(Sjmax(il),Scrit(il))
             Sjmin(il) = min( (Sbef(il)+Sij(il,i,j))/2. ,
     1                                               Sij(il,i,j) )
             Sjmin(il) = max(Sjmin(il),Scrit(il))
             Sbef(il) = Sij(il,i,j)
           ENDIF
          IF (abs(Sjmin(il)-Sjmax(il)) .GT. 1.e-10) Sup(il)=
     1                            max(Sjmin(il),Sjmax(il),Sup(il))
      endif
      endif
784   continue
c -----------------------------------------------
         END IF
c -----------------------------------------------
c
c
      do il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then
       if(sij(il,i,j).gt.0.0)then
         rti=qnk(il)-ep(il,i)*clw(il,i)
         Qmixmax(il)=Qmix(Sjmax(il))
         Qmixmin(il)=Qmix(Sjmin(il))
         Rmixmax(il)=Rmix(Sjmax(il))
         Rmixmin(il)=Rmix(Sjmin(il))
         SQmRmax(il)= Sjmax(il)*Qmix(Sjmax(il))-Rmix(Sjmax(il))
         SQmRmin(il)= Sjmin(il)*Qmix(Sjmin(il))-Rmix(Sjmin(il))
c
         Ment(il,i,j) = abs(Qmixmax(il)-Qmixmin(il))*Ment(il,i,j)
c
c    Sigij(i,j) is the 'true' mixing fraction of mixture Ment(i,j)
         IF (abs(Qmixmax(il)-Qmixmin(il)) .GT. 1.e-10) THEN
           Sigij(il,i,j) =
     :           (SQmRmax(il)-SQmRmin(il))/(Qmixmax(il)-Qmixmin(il))
         ELSE
           Sigij(il,i,j) = 0.
         ENDIF
c
c --    Compute Qent, uent, vent according to the true mixing fraction
        Qent(il,i,j) = (1.-Sigij(il,i,j))*rti
     :               + Sigij(il,i,j)*rr(il,i)
        uent(il,i,j) = (1.-Sigij(il,i,j))*unk(il)
     :               + Sigij(il,i,j)*u(il,i)
        vent(il,i,j) = (1.-Sigij(il,i,j))*vnk(il)
     :               + Sigij(il,i,j)*v(il,i)
c
c--     Compute liquid water static energy of mixed draughts
c       IF (j .GT. i) THEN
c        awat=elij(il,i,j)-(1.-ep(il,j))*clw(il,j)
c        awat=amax1(awat,0.0)
c       ELSE
c        awat = 0.
c       ENDIF
c       Hent(il,i,j) = (1.-Sigij(il,i,j))*HP(il,i)
c    :         + Sigij(il,i,j)*H(il,i)
c    :         + (LV(il,j)+(cpd-cpv)*t(il,j))*awat
cIM 301008 beg
        Hent(il,i,j) = (1.-Sigij(il,i,j))*HP(il,i)
     :         + Sigij(il,i,j)*H(il,i)

        Elij(il,i,j) = Qent(il,i,j)-rs(il,j)
        Elij(il,i,j) = Elij(il,i,j)
     :    + ((h(il,j)-Hent(il,i,j))*rs(il,j)*LV(il,j)
     :    / ((cpd*(1.-Qent(il,i,j))+Qent(il,i,j)*cpv)
     :    * rrv*t(il,j)*t(il,j)))
        Elij(il,i,j) = Elij(il,i,j)
     :    / (1.+LV(il,j)*LV(il,j)*rs(il,j)
     :    / ((cpd*(1.-Qent(il,i,j))+Qent(il,i,j)*cpv)
     :    * rrv*t(il,j)*t(il,j)))

        Elij(il,i,j) = max(elij(il,i,j),0.)

        Elij(il,i,j) = min(elij(il,i,j),Qent(il,i,j))

        IF (j .GT. i) THEN
         awat=elij(il,i,j)-(1.-ep(il,j))*clw(il,j)
         awat=amax1(awat,0.0)
        ELSE
         awat = 0.
        ENDIF

c        print *,h(il,j)-hent(il,i,j),LV(il,j)*rs(il,j)/(cpd*rrv*t(il,j)*
c    :         t(il,j))

        Hent(il,i,j) =  Hent(il,i,j)+(LV(il,j)+(cpd-cpv)*t(il,j))
     :         * awat
cIM 301008 end
c
c      print *,'mix : i,j,hent(il,i,j),sigij(il,i,j) ',
c     :               i,j,hent(il,i,j),sigij(il,i,j)
c
c --      ASij is the integral of P(F) over the relevant F interval
         ASij(il) = ASij(il)
     1               + abs(Qmixmax(il)*(1.-Sjmax(il))+Rmixmax(il)
     1                    -Qmixmin(il)*(1.-Sjmin(il))-Rmixmin(il))
c
      endif
      endif
      enddo
       do k=1,ntra
         do il=1,ncum
          if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
     :       (j.ge.(icb(il)-1)).and.(j.le.inb(il))
     :      .and. lwork(il) ) then
          if(sij(il,i,j).gt.0.0)then
            traent(il,i,j,k)=sigij(il,i,j)*tra(il,i,k)
     :            +(1.-sigij(il,i,j))*tra(il,nk(il),k)
          endif
          endif
         enddo
       enddo
c
c --    If I=J (detrainement and entrainement at the same level), then only the
c --    adiabatic ascent part of the mixture is considered
        IF (I .EQ. J) THEN
      do il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and.
     :      j.ge.(icb(il)-1) .and. j.le.inb(il)
     :      .and. lwork(il) ) then
       if(sij(il,i,j).gt.0.0)then
          rti=qnk(il)-ep(il,i)*clw(il,i)
ccc          Ment(il,i,i) = m(il,i)*abs(Qmixmax(il)*(1.-Sjmax(il))
          Ment(il,i,i) = abs(Qmixmax(il)*(1.-Sjmax(il))
     1                    +Rmixmax(il)
     1                    -Qmixmin(il)*(1.-Sjmin(il))-Rmixmin(il))
          Qent(il,i,i) = rti
          uent(il,i,i) = unk(il)
          vent(il,i,i) = vnk(il)
          Hent(il,i,i) = hp(il,i)
          Elij(il,i,i) = clw(il,i)*(1.-ep(il,i))
          Sigij(il,i,i) = 0.
      endif
      endif
      enddo
       do k=1,ntra
         do il=1,ncum
          if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
     :       (j.ge.(icb(il)-1)).and.(j.le.inb(il))
     :      .and. lwork(il) ) then
          if(sij(il,i,j).gt.0.0)then
            traent(il,i,i,k)=tra(il,nk(il),k)
          endif
          endif
         enddo
       enddo
c
        ENDIF
c
175   continue

      do il=1,ncum
       if (i.ge.icb(il).and.i.le.inb(il).and.lwork(il)) then
        asij(il)=amax1(1.0e-16,asij(il))
        asij(il)=1.0/asij(il)
        csum(il,i)=0.0
       endif
      enddo

      do 180 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         ment(il,i,j)=ment(il,i,j)*asij(il)
        endif
       enddo
180   continue

      do 197 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :   .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         csum(il,i)=csum(il,i)+ment(il,i,j)
        endif
       enddo
197   continue

      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :     .and. csum(il,i).lt.1. ) then
ccc     :     .and. csum(il,i).lt.m(il,i) ) then
        nent(il,i)=0
ccc        ment(il,i,i)=m(il,i)
        ment(il,i,i)=1.
        qent(il,i,i)=qnk(il)-ep(il,i)*clw(il,i)
        uent(il,i,i)=unk(il)
        vent(il,i,i)=vnk(il)
        elij(il,i,i)=clw(il,i)*(1.-ep(il,i))
        sij(il,i,i)=0.0
       endif
      enddo ! il

      do j=1,ntra
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
     :     .and. csum(il,i).lt.1. ) then
ccc     :     .and. csum(il,i).lt.m(il,i) ) then
         traent(il,i,i,j)=tra(il,nk(il),j)
        endif
       enddo
      enddo
c
789   continue
c
      return
      end

