       SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
c
c
c    Multiplication of matrix A by matrix B ,and result in matrix C
c
c
c      NAR : number of   rows  for matrix  A and  matrix  C
c      NAC : number of columns for matrix  A and number of rows
c                              for matrix  B
c      NBC : number of columns for matrix  B and  matrix  C
c
c
       IMPLICIT NONE
c
       INTEGER  NAR,NAC,NBC, I,J,K
       REAL     A(NAR,NAC), B(NAC,NBC), C(NAR,NBC)
c
       DO 10 J = 1,NBC
       DO 10 I = 1,NAR
       C(I,J)=0.
10     CONTINUE
c
       DO 50 J = 1,NBC
       DO 40 I = 1,NAR

       DO 20 K = 1,NAC
       C(I,J)  = C(I,J) + A(I,K) * B(K,J)
20     CONTINUE

40     CONTINUE
50     CONTINUE

       RETURN
       END
