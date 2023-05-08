!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/hgardfou.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE hgardfou (t,tsol,text)

      use dimphy
      IMPLICIT none
c======================================================================
c Verifier la temperature
c======================================================================
#include "YOMCST.h"
      REAL t(klon,klev), tsol(klon)
      CHARACTER*(*) text
C
      INTEGER i, k
      REAL zt(klon)
      INTEGER jadrs(klon), jbad
      LOGICAL ok
c
      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall /.TRUE./
      IF (firstcall) THEN
         PRINT*, 'hgardfou garantit la temperature dans [20,1200] K'
         firstcall = .FALSE.
      ENDIF
c
      ok = .TRUE.
      DO k = 1, klev
         DO i = 1, klon
            zt(i) = t(i,k)
         ENDDO
#ifdef CRAY
         CALL WHENFGT(klon, zt, 1, 1200.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).GT.1200.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,k,temperature =',jadrs(i),k,zt(jadrs(i))
           ENDDO
         ENDIF
#ifdef CRAY
         CALL WHENFLT(klon, zt, 1, 20.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).LT.20.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,k,temperature =',jadrs(i),k,zt(jadrs(i))
           ENDDO
         ENDIF
      ENDDO
c
         DO i = 1, klon
            zt(i) = tsol(i)
         ENDDO
#ifdef CRAY
         CALL WHENFGT(klon, zt, 1, 1200.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).GT.1200.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,temperature =',jadrs(i),zt(jadrs(i))
           ENDDO
         ENDIF
#ifdef CRAY
         CALL WHENFLT(klon, zt, 1, 20.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).LT.20.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,temperature =',jadrs(i),zt(jadrs(i))
           ENDDO
         ENDIF
c
      IF (.NOT. ok) THEN
         PRINT*, 'hgardfou s arrete ', text
         CALL abort
      ENDIF

      RETURN
      END
