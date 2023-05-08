      INTEGER FUNCTION lnblnk (letter)

C--------------------------------------------------------
C Fonction qui determine la longeur d'un string sans les
C blancs qui suivent. Le critere pour determiner la fin du
C string est, trois blancs de suite
C---------------------------------------------------------
C     ARGUMENTS
C     +++++++++
C     letter: CHARACTER*xxx (xxx < imax)
C             le string dont on determine la longuer
C     lnblnk: INTEGER
C             le nombre de characteres
C
C     PARAMETER
C     +++++++++
C     imax : INTEGER
C            le nombre maximale de character que peut contenir le string
C            a traiter

      IMPLICIT NONE
      INTEGER i,imax
      PARAMETER (imax = 256)
      CHARACTER(len=*) letter

      i=0

10    i=i+1
      IF (letter(i:i+3) . EQ . '   ') GOTO 20
      GOTO 10

20    lnblnk=i-1

      RETURN
      END

