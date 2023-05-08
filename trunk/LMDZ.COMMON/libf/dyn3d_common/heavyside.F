!
! $Header$
!
c
c
       FUNCTION heavyside(a)

c      ...   P. Le Van  ....
c
       IMPLICIT NONE

       REAL(KIND=8) heavyside , a

       IF ( a.LE.0. )  THEN
         heavyside = 0.
       ELSE
         heavyside = 1.
       ENDIF

       RETURN
       END


