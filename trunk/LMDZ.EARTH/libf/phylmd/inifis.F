!
! $Id: inifis.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE inifis(ngrid,nlayer,
     $           punjours,
     $           pdayref,ptimestep,
     $           plat,plon,parea,
     $           prad,pg,pr,pcpp)
      USE dimphy
      IMPLICIT NONE
c
c=======================================================================
c
c   subject:
c   --------
c
c   Initialisation for the physical parametrisations of the LMD 
c   martian atmospheric general circulation modele.
c
c   author: Frederic Hourdin 15 / 10 /93
c   -------
c
c   arguments:
c   ----------
c
c   input:
c   ------
c
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    pdayref               Day of reference for the simulation
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   declarations:
c   -------------
 
cym#include "dimensions.h"
cym#include "dimphy.h"

      INCLUDE 'iniprint.h'
      REAL prad,pg,pr,pcpp,punjours
 
      INTEGER ngrid,nlayer
      REAL plat(ngrid),plon(ngrid),parea(klon)
      INTEGER pdayref
 
      REAL ptimestep
      CHARACTER (LEN=20) :: modname='inifis'
      CHARACTER (LEN=80) :: abort_message

 
      IF (nlayer.NE.klev) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'nlayer     = ',nlayer
         PRINT*,'klev   = ',klev
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      ENDIF

      IF (ngrid.NE.klon) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'klon   = ',klon
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      ENDIF

      RETURN
9999  continue
      abort_message = 'Cette version demande les fichier rnatur.dat
     & et surf.def'
      CALL abort_gcm (modname,abort_message,1)

      END
