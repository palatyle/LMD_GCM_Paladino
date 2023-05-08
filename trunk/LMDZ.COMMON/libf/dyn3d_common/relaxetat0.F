      SUBROUTINE relaxetat0(fichnom)

      use netcdf, only: nf90_get_var
      use comuforc_h
      use ioipsl_getincom

      IMPLICIT NONE

!=======================================================================
! Author: A. Spiga
! Purpose: Read a simple relaxation field
!=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "netcdf.inc"
#include "iniprint.h"

! INPUTS
      CHARACTER*(*) fichnom

! VARIABLES
      REAL ucov(iip1, jjp1,llm)
      REAL press(llm)
      INTEGER ierr, nid, nvarid
      INTEGER i,j,l,ij

      !! OPEN NETCDF FILE
      ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        write(lunout,*)'relaxetat0: Pb d''ouverture du fichier'
        write(lunout,*)' ierr = ', ierr
        CALL ABORT
      ENDIF

      !! ALLOCATE RELAXATION FIELD
      IF (.NOT.ALLOCATED(uforc)) ALLOCATE(uforc(ip1jmp1,llm))

      !! READ RELAXATION FIELD
      ierr = NF_INQ_VARID (nid, "ucov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <ucov> est absent"
         CALL abort
      ENDIF
      ierr = nf90_get_var(nid, nvarid, ucov) 
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <ucov>"
         CALL abort
      ENDIF

      !! TRANSFER TO SHARED VARIABLE IN comuforc_h MODULE
      !! TBD --> check the relaxation field's dimensions are OK!
      do l=1,llm
        do j=1,jjp1
           do i=1,iip1
              ij=(j-1)*iip1+i
              uforc(ij,l)=ucov(i,j,l)
           enddo
        enddo
      enddo

      !! READ THE RELAXATION RATE IN SETTINGS
      facwind = 1000.
      CALL getin('facwind',facwind)

      !! CALCULATE ATTENUATION FACTOR WITH ALTITUDE
      ierr = NF_INQ_VARID (nid, "presnivs", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Le champ <presnivs> est absent"
         CALL abort
      ENDIF
      ierr = nf90_get_var(nid, nvarid, press)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynetat0: Lecture echouee pour <press>"
         CALL abort
      ENDIF
      phalf = 1.e5 !! at this pressure, attenua is 0.5
      CALL getin('phalf',phalf)
      !! attenua goes from 1 at lowermost levels to 0 at uppermost levels
      IF (.NOT.ALLOCATED(attenua)) ALLOCATE(attenua(llm))
      attenua = 0.5*(1.+tanh(log10(press/phalf)))

      !! END
      ierr = NF_CLOSE(nid)
      RETURN
      END
