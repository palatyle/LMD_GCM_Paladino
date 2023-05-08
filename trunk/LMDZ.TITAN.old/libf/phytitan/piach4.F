      SUBROUTINE PIACH4(NU,NT,XROT,XTRA,T)
      IMPLICIT REAL (A-H,O-Z)
      REAL NU,I10,NSQR,K
      DIMENSION T(100)
      DIMENSION FOX(60),XJ(23),G(23),E(23),WW(3,23),RHO(100,23),
     $I10(100),XROT(100),XTRA(100)

      save ww,a2,pi,c,hck,tau1,tau2,cg,rho,c1,k,asqr,quads1,i10,xj

      DATA FOX /1.00,1.00,4.12,1.00,0.52,0.81,1.00,0.70,0.75,0.80,
     $0.69,1.24,0.44,1.29,1.05,0.56,0.79,0.79,0.89,1.16,1.27,0.86,
     $0.95,1.00,0.72,1.16,0.92,1.00,0.92,1.18,0.90,1.18,1.09,0.85,
     $0.96,0.93,0.82,1.18,1.17,0.95,0.96,1.00,0.85,1.06,0.90,15*1.00/
      DATA JMAX1/20/,JMAX2/23/,PI/3.141592654/,H/1.05459E-27/,
     $C/2.997925E10/,K /1.38054E-16/,HCK/1.43892/,XNLOS/2.687E19/,          276.
     $ASQR/6.7081/,EPSK/148.6/,B/5.25/,ANORM/4.797E17/
      DO 8 IT=1,NT
      X=4.*EPSK/T(IT)
      I10(IT)=FI10(X)
    8 CONTINUE
      DO 10 J=1,JMAX2
      XJ(J)=FLOAT(J-1)
      E(J)=B*XJ(J)*(XJ(J)+1.)
   10 CONTINUE
      G(1)=5.
      G(2)=3.
      G(3)=5.
      G(4)=11.
      G(5)=13.
      G(6)=11.
!!! correction au nez
      DO 15 J=7,JMAX2
   15 G(J)=G(J-6)+16.
      DO 20 J=1,JMAX1
      DO 20 I=1,3
      WW(I,J)=2.*PI*C*(E(J+I)-E(J))
   20 CONTINUE
      NSQR=(XNLOS*1.0E-24)**2
c     C1=64./63.*PI*PI*NSQR/H/C/ANORM
      C1=64./63.*PI*PI*NSQR/H/C
      QUADS1=11.10
      TAU1=9.00E-14
      TAU2=13.60E-14
      DO 90 IT=1,NT
      SUM=0.
      DO 50 J=1,JMAX2
      ARG=HCK*E(J)/T(IT)
      RHO(IT,J)=EXP(-ARG)
      SUM=SUM+(2.*XJ(J)+1.)*G(J)*RHO(IT,J)
   50 CONTINUE
      DO 55 J=1,JMAX2
      RHO(IT,J)=RHO(IT,J)/SUM
   55 CONTINUE
   90 CONTINUE
      RETURN
      ENTRY OPACH4(NU,NT,XROT,XTRA,T)
      W=2.*PI*C*NU
      DO 200 IT=1,NT
      XROT(IT)=0.
      XTRA(IT)=0.
      CZ=HCK*NU/T(IT)
      CW=W*(1.-EXP(-CZ))
      DO 125 J=1,JMAX1
      XTRA(IT)=XTRA(IT)+(2.*XJ(J)+1.)*(2.*XJ(J)+1.)*
     $ RHO(IT,J)*GAMMB(W,T(IT),TAU1,TAU2)
!!! avant cette boucle contenait K au lieu de IK en contradiction
!!! avec K Boltzman
      DO 25 I=1,3
      IK=3*(J-1)+I
      XROT(IT)=XROT(IT)+FOX(IK)*(2.*XJ(J)+1.)*(2.*XJ(J+I)+1.)*
     $(RHO(IT,J  )*GAMMB(W-WW(I,J),T(IT),TAU1,TAU2)
     $+RHO(IT,J+I)*GAMMB(W+WW(I,J),T(IT),TAU1,TAU2))
   25 CONTINUE
  125 CONTINUE
c     COEF=CW*C1*IK*ASQR*QUADS1*I10(IT)
      COEF=CW*C1*K*ASQR*QUADS1*I10(IT)
      XROT(IT)=XROT(IT)*COEF
      XTRA(IT)=XTRA(IT)*COEF
  200 CONTINUE
      RETURN
      END
