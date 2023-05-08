      SUBROUTINE PIAN2(NU,NT,XROT,XTRA,T)
C
C   THIS ROUTINE WAS FIRST WRITTEN BY REGIS COURTIN IN 1980 :
C   THE PARAMETERS WERE DETERMINED BY THE MEASUREMENTS OF
C   BUONTEMPO ET AL (1975), JOURNAL OF CHEMICAL PHYSICS 63, 2570.
C   LAST MODIFIED BY REGIS COURTIN (APRIL 23,1986) :
C   THE QUADS1 PARAMETER WAS ADJUSTED TO FIT THE MEASUREMENTS OF
C   DAGG ET AL (1985), CANADIAN JOURNAL OF PHYSICS 63, 625.
C
      IMPLICIT REAL (A-H,O-Z)
      REAL NU,I8,NSQR,K
      real c
      DIMENSION T(100)
      DIMENSION XROT(100),XTRA(100)
      DIMENSION XJ(32),G(32),E(32),WW(32),CG(32),RHO(100,32),
     $I8(100),QUADS1(100),TAU1(100),TAU2(100),A2(100)
      save ww,a2,pi,c,hck,tau1,tau2,cg,rho,c1,k,asqr,quads1,i8
      DATA JMAX1/30/,JMAX2/32/,PI/3.141592654/,H/1.05459E-27/,
     $C/2.997925E10/,K/1.38054E-16/,HCK/1.43892/,XNLOS/2.687E19/,
     $ASQR/3.3124/,EPSK/71.4/,PW1/.639/,PW2/.423/,
     $DSDT/-.0204/
      DO 8 IT=1,NT
      X=4.*EPSK/T(IT)
      I8(IT)=FI8(X)
    8 CONTINUE
      DO 10 J=1,JMAX2
      XJ(J)=FLOAT(J-1)
      G(J)=6.
      IF (MOD(J,2).EQ.0) G(J)=3.
      E(J)=AZENER(0.,XJ(J))
   10 CONTINUE
      DO 20 J=1,JMAX1
      CG(J)=1.5*(XJ(J)+1.)*(XJ(J)+2.)/(2.*XJ(J)+3.)
      WW(J)=2.*PI*C*(E(J+2)-E(J))
   20 CONTINUE
      NSQR=(XNLOS*1.0E-24)**2
      C1=4.*PI*PI*NSQR/H/C
      DO 90 IT=1,NT
      QUADS1(IT)=27.72+DSDT*T(IT)
      TAU1(IT)=9.797E-12/(T(IT)**PW1)
      TAU2(IT)=1.518E-12/(T(IT)**PW2)
      SUM=0.
      DO 50 J=1,JMAX2
      ARG=HCK*E(J)/T(IT)
      RHO(IT,J)=EXP(-ARG)
      SUM=SUM+(2.*XJ(J)+1.)*G(J)*RHO(IT,J)
   50 CONTINUE
      DO 55 J=1,JMAX2
      RHO(IT,J)=G(J)*RHO(IT,J)/SUM
   55 CONTINUE
      A2(IT)=0.
      DO 60 J=1,JMAX1
      A2(IT)=A2(IT)+XJ(J)*(XJ(J)+1.)*(2.*XJ(J)+1.)*RHO(IT,J)/
     $  ((2.*XJ(J)-1.)*(2.*XJ(J)+3.))
   60 CONTINUE
   90 CONTINUE
c     print*,'A2 dans pia...',a2
c     print*,'var dans pia...',
c    s    pi,c,hck,tau1,tau2,cg,rho,c1,k,asqr,quads1,i8

      RETURN
c!!!!  ENTRY OPAN2(NU,NT,XROT,XTRA)
      ENTRY OPAN2(NU,NT,XROT,XTRA,T)
c     print*,'A2 dans opa...',a2
c     print*,'var dans opa...',
c    s    pi,c,hck,tau1,tau2,cg,rho,c1,k,asqr,quads1,i8
      W=2.*PI*C*NU
      DO 200 IT=1,NT
      CZ=HCK*NU/T(IT)
      CW=W*(1.-EXP(-CZ))
      XTRA(IT)=A2(IT)*GAMMB(W,T(IT),TAU1(IT),TAU2(IT))
      XROT(IT)=0.
      DO 125 J=1,JMAX1
      XROT(IT)=XROT(IT)+CG(J)*
     $(RHO(IT,J  )*GAMMB(W-WW(J),T(IT),TAU1(IT),TAU2(IT))
     $+RHO(IT,J+2)*GAMMB(W+WW(J),T(IT),TAU1(IT),TAU2(IT)))
  125 CONTINUE
      COEF=CW*C1*K*ASQR*QUADS1(IT)*I8(IT)
      XROT(IT)=XROT(IT)*COEF
      XTRA(IT)=XTRA(IT)*COEF
  200 CONTINUE
      RETURN
      END
