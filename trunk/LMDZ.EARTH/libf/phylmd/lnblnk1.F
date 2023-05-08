!
! $Header$
!
      INTEGER FUNCTION lnblnk1 (letter)

C--------------------------------------------------------
C Fonction qui determine la longeur d'un string sans les
C blancs qui suivent. Le critere pour determiner la fin du
C string est, trois blancs de suite
C---------------------------------------------------------
C     ARGUMENTS
C     +++++++++
C     letter: CHARACTER*xxx (xxx < imax)
C             le string dont on determine la longuer
C     lnblnk1: INTEGER
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
c     CHARACTER*256 letter
      CHARACTER*4 letter

      i=0

10    i=i+1
c     IF (letter(i:i+1) . EQ . ' ') THEN
      IF (letter(i:i) . EQ . ' ') THEN
c      print*,'i=',i,'letter(i:i+1)=',letter(i:i+1)
c      print*,'i=',i
       GOTO 20
      ELSE
       GOTO 10
      ENDIF

20    lnblnk1=i-1
c     PRINT*,'lnblnk1=',lnblnk1

      RETURN
      END

