!
! $Id: iniphysiq.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
c
      SUBROUTINE iniphysiq(ngrid,nlayer,
     $           punjours,
     $           pdayref,ptimestep,
     $           plat,plon,parea,pcu,pcv,
     $           prad,pg,pr,pcpp)
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE comgeomphy

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
cym#include "comgeomphy.h"
#include "YOMCST.h"
      REAL prad,pg,pr,pcpp,punjours
 
      INTEGER ngrid,nlayer
      REAL plat(ngrid),plon(ngrid),parea(klon_glo)
      REAL pcu(klon_glo),pcv(klon_glo)
      INTEGER pdayref
      INTEGER :: ibegin,iend,offset
 
      REAL ptimestep
      CHARACTER (LEN=20) :: modname='iniphysiq'
      CHARACTER (LEN=80) :: abort_message
 
      IF (nlayer.NE.klev) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'nlayer     = ',nlayer
         PRINT*,'klev   = ',klev
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      ENDIF

      IF (ngrid.NE.klon_glo) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'klon   = ',klon_glo
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      ENDIF
c$OMP PARALLEL PRIVATE(ibegin,iend) 
c$OMP+         SHARED(parea,pcu,pcv,plon,plat)
      
      offset=klon_mpi_begin-1
      airephy(1:klon_omp)=parea(offset+klon_omp_begin:
     &                          offset+klon_omp_end)
      cuphy(1:klon_omp)=pcu(offset+klon_omp_begin:offset+klon_omp_end)
      cvphy(1:klon_omp)=pcv(offset+klon_omp_begin:offset+klon_omp_end)
      rlond(1:klon_omp)=plon(offset+klon_omp_begin:offset+klon_omp_end)
      rlatd(1:klon_omp)=plat(offset+klon_omp_begin:offset+klon_omp_end)

      call suphel

c$OMP END PARALLEL

      print*,'ATTENTION !!! TRAVAILLER SUR INIPHYSIQ'
      print*,'CONTROLE DES LATITUDES, LONGITUDES, PARAMETRES ...'

      RETURN
9999  CONTINUE
      abort_message ='Cette version demande les fichier rnatur.dat
     & et surf.def'
      CALL abort_gcm (modname,abort_message,1)

      END
