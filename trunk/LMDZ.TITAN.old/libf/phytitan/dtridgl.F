      SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)
C DOUBLE PRESCISION VERSION OF TRIDGL
c     PARAMETER (NMAX=201)
      PARAMETER (NMAX=401)
      IMPLICIT REAL  (A-H,O-Z)
      DIMENSION AF(L),BF(L),CF(L),DF(L),XK(L)
      DIMENSION AS(NMAX),DS(NMAX)
C* THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
C*  EQUATIONS. THE FORM OF THE EQUATIONS ARE:
C*  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
C* WHERE I=1,L  LESS THAN 103.
C* ..............REVIEWED -CP........
      AS(L) = AF(L)/BF(L)
      DS(L) = DF(L)/BF(L)
      DO 10 I=2,L
           X=1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
           AS(L+1-I)=AF(L+1-I)*X
           DS(L+1-I)=(DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
   10 CONTINUE
      XK(1)=DS(1)
      DO 20 I=2,L
           XKB=XK(I-1)
           XK(I)=DS(I)-AS(I)*XKB
   20 CONTINUE
  910 FORMAT(/,8X,'AF(I)',7X,'BF(I)',7X,'CF(I)',7X,'DF(I)',7X,
     *         'AS(I)',7X,'DS(I)',7X,'XK(I)',/)
  915 FORMAT((3X,7(1X,1PE11.4)))
      RETURN
      END
