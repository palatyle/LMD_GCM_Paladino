      SUBROUTINE MXV(A,NRA,B,NCA,C)
      IMPLICIT NONE

      INTEGER NRA,NCA
      REAL A(NRA,NCA), B(NCA), C(NRA)

      INTEGER I,J

c     initialize product
      DO 100 I = 1, NRA
         C(I) = 0.
100   CONTINUE
c     multiply matrix A and vector B
      DO 220 J = 1, NCA
         DO 210 I = 1, NRA
            C(I) = C(I) + A(I,J)*B(J)
210     CONTINUE
220   CONTINUE
      RETURN
      END
