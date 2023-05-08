      SUBROUTINE MXVA (SA,IAC,IAR,SB,IB,SC,IC,NRA,NCA)
      IMPLICIT NONE

      INTEGER IAC,IAR,IB,IC,NRA,NCA
      REAL SA(*), SB(*), SC(*)

      INTEGER I,J

c     initialize product
      DO 100 I = 1, NRA
       SC( 1 + (I-1)*IC ) = 0.
c      ( C(i) := 0. )
100   CONTINUE
c     multiply matrix from sa and vector from sb
      DO 220 J = 1, NCA
        DO 210 I = 1, NRA
          SC( 1 + (I-1)*IC )
     &      = SC( 1 + (I-1)*IC )
     &        + SA( 1 + (I-1)*IAC + (J-1)*IAR )
     &          * SB( 1 + (J-1)*IB )
c         ( C(i) := C(i) + A(i,j)*B(j) )
210      CONTINUE
220   CONTINUE
      RETURN
      END
