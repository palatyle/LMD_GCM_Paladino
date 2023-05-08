c
c $Header: /users/lmdz/cvsroot/LMDZ.3.3/libf/dyn3d/heavyside.F,v 1.2 2001/04/13 15:49:10 lmdz Exp $
c
       FUNCTION heavyside(a)

c      ...   P. Le Van  ....
c
       IMPLICIT NONE

       REAL*8 heavyside , a

       IF ( a.LE.0. )  THEN
         heavyside = 0.
       ELSE
         heavyside = 1.
       ENDIF

       RETURN
       END


