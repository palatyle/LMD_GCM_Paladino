      FUNCTION GAMFCN(W,T,TAU1,TAU2)
      IMPLICIT REAL (A-H,O-Z)
C
C     COMPUTES THE LINE SHAPE FOR PRESSURE-INDUCED H2-H2 AND H2-HE
C     TRANSITIONS, FROM THE SEMI-EMIRICAL FORMULAE OF BIRNBAUM ET AL.
C     SEE EG., GEORGE BIRNBAUM AND E. RICHARD COHEN, CANADIAN JOURNAL
C     OF PHYSICS, VOL. 54, 593 (1976).
C     NOTE THAT "GAMFCN" IS A POOR NAME FOR THIS ROUTINE; THE GAMMA
C     OF BIRNBAUM AND COHEN IS NOT THE USUAL "GAMMA FUNCTION".
C
      save hk,pi,tau10,tau20
      DATA HK/7.638315E-12/,PI/3.141593/,TAU10/0.0/,TAU20/0.0/
      logical first
      data first/.true./
      save first

C     NOTE: HK = 1.05450E-27 / 1.38054E-16
C
C***********************************************************************
      save tau12,z2,hbh

      if (first)
     s  print*,'WARNING!!! ON rajoute des valeurs a 0.',
     s   'Est-ce bien raisonable... dans GAMFCN'
      first=.false.
      IF (TAU1 .NE. TAU10) GO TO 10
      IF (TAU2 .EQ. TAU20) GO TO 20
 10   TAU12 = TAU1 * TAU1
      TAU22 = TAU2 * TAU2
      HBH = 0.5 * HK / T
      Z2 = SQRT(TAU22 + HBH**2) / TAU1
      TAU10 = TAU1
      TAU20 = TAU2
 20   WSQR = W * W
      Z = SQRT(1.0+WSQR*TAU12) * Z2
      IF (Z .LE. 1.0) GO TO 50
C     COMPUTE K1 BESSEL FUNCTION USING POLYNOMIAL APPROXIMATION
      A = 1.0 / Z
      BK1 = 1.253314 + .4699927*A
      B = A * A
      BK1 = BK1 - .1468583*B
      B = B * A
      BK1 = BK1 + .1280427*B
      B = B * A
      BK1 = BK1 - .1736432*B
      B = B * A
      BK1 = BK1 + .2847618*B
      B = B * A
      BK1 = BK1 - .4594342*B
      B = B * A
      BK1 = BK1 + .6283381*B
      B = B * A
      BK1 = BK1 - .6632295*B
      B = B * A
      BK1 = BK1 + .5050239*B
      B = B * A
      BK1 = BK1 - .2581304*B
      B = B * A
      BK1 = BK1 + .7880001E-01*B
      B = B * A
      BK1 = BK1 - .1082418E-01*B
      BK1 = EXP(-Z) * BK1 * SQRT(A)
      GO TO 100
C     COMPUTE K1 BESSEL FUNCTION USING SERIES EXPANSION
 50   A = 0.5 * Z
      B = .5772157 + LOG(A)
      C = A * A
      BK1 = 1.0/Z + A*(B-0.5)
      A = A * C
      BK1 = BK1 + A*.2500000E+00*(0.5+(B-1.500000)*2.0)
      A = A * C
      BK1 = BK1 + A*.2777777E-01*(0.5+(B-1.833333)*3.0)
      A = A * C
      BK1 = BK1 + A*.1736110E-02*(0.5+(B-2.083333)*4.0)
      A = A * C
      BK1 = BK1 + A*.6944439E-04*(0.5+(B-2.283333)*5.0)
      A = A * C
      BK1 = BK1 + A*.1929009E-05*(0.5+(B-2.449999)*6.0)
      A = A * C
      BK1 = BK1 + A*.3936752E-07*(0.5+(B-2.592855)*7.0)
      A = A * C
      BK1 = BK1 + A*.6151173E-09*(0.5+(B-2.717855)*8.0)
 100  CONTINUE
      GAMFCN = TAU1/PI * EXP(TAU2/TAU1+HBH*W) * Z*BK1 / (1.0+WSQR*TAU12)
      RETURN
      END
