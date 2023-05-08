!
! $Id: hgardfou.F 1279 2009-12-10 09:02:56Z fairhead $
      SUBROUTINE hgardfou (t,tsol,text)
      use dimphy
      use phys_state_var_mod
      IMPLICIT none
c======================================================================
c Verifier la temperature
c======================================================================
#include "dimensions.h"
#include "YOMCST.h"
#include "indicesol.h"
      REAL t(klon,klev), tsol(klon,nbsrf)
      CHARACTER*(*) text
      character (len=20) :: modname = 'hgardfou'
      character (len=80) :: abort_message
C
      INTEGER i, k, nsrf
      REAL zt(klon)
      INTEGER jadrs(klon), jbad
      LOGICAL ok
c
      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall /.TRUE./
c$OMP THREADPRIVATE(firstcall)

      IF (firstcall) THEN
         PRINT*, 'hgardfou garantit la temperature dans [100,370] K'
         firstcall = .FALSE.
c        DO i = 1, klon
c         print*,'i=',i,'rlon=',rlon(i),'rlat=',rlat(i)
c        ENDDO
c
      ENDIF
c
      ok = .TRUE.
      DO k = 1, klev
         DO i = 1, klon
            zt(i) = t(i,k)
         ENDDO
#ifdef CRAY
         CALL WHENFGT(klon, zt, 1, 370.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).GT.370.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,k,temperature,lon,lat,pourc ter,oce,lic,sic =',
     $       jadrs(i),k,zt(jadrs(i)),rlon(jadrs(i)),rlat(jadrs(i)),
     $       (pctsrf(jadrs(i),nsrf),nsrf=1,nbsrf)
           ENDDO
         ENDIF
#ifdef CRAY
         CALL WHENFLT(klon, zt, 1, 100.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
!         IF (zt(i).LT.100.0) THEN
         IF (zt(i).LT.50.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
             PRINT *,'i,k,temperature,lon,lat,pourc ter,oce,lic,sic =',
     $       jadrs(i),k,zt(jadrs(i)),rlon(jadrs(i)),rlat(jadrs(i)),
     $       (pctsrf(jadrs(i),nsrf),nsrf=1,nbsrf)
           ENDDO
         ENDIF
      ENDDO
c
      DO nsrf = 1, nbsrf
         DO i = 1, klon
            zt(i) = tsol(i,nsrf)
         ENDDO
#ifdef CRAY
         CALL WHENFGT(klon, zt, 1, 370.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
         IF (zt(i).GT.370.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
            PRINT *,'i,nsrf,temperature,lon,lat,pourc ter,oce,lic,sic ='
     $      ,jadrs(i),nsrf,zt(jadrs(i)),rlon(jadrs(i)),rlat(jadrs(i))
     $      ,pctsrf(jadrs(i),nsrf)
           ENDDO
         ENDIF
#ifdef CRAY
         CALL WHENFLT(klon, zt, 1, 100.0, jadrs, jbad)
#else
         jbad = 0
         DO i = 1, klon
!         IF (zt(i).LT.100.0) THEN
         IF (zt(i).LT.50.0) THEN
            jbad = jbad + 1
            jadrs(jbad) = i
         ENDIF
         ENDDO
#endif
         IF (jbad .GT. 0) THEN
           ok = .FALSE.
           DO i = 1, jbad
            PRINT *,'i,nsrf,temperature,lon,lat,pourc ter,oce,lic,sic ='
     $      ,jadrs(i),nsrf,zt(jadrs(i)),rlon(jadrs(i)),rlat(jadrs(i))
     $      ,pctsrf(jadrs(i),nsrf)
           ENDDO
         ENDIF
      ENDDO
c
      IF (.NOT. ok) THEN
         abort_message= 'hgardfou s arrete '//text
         CALL abort_gcm (modname,abort_message,1)
      ENDIF

      RETURN
      END
