      SUBROUTINE surfini(ngrid,nq,qsurf,albedo,albedo_bareground)

      USE surfdat_h, only: albedodat
      use planetwide_mod, only: planetwide_maxval, planetwide_minval
      use radinc_h, only : L_NSPECTV

      IMPLICIT NONE
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccc                                                                 cccccccccccccc
cccccccccccccc   Spectral Albedo Initialisation - Routine modified by MT2015.  cccccccccccccc
cccccccccccccc                                                                 cccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c--------------------
c   Declarations:
c--------------------

      INTEGER,INTENT(IN) :: ngrid
      INTEGER,INTENT(IN) :: nq
      REAL,INTENT(OUT) :: albedo(ngrid,L_NSPECTV)
      REAL,INTENT(OUT) :: albedo_bareground(ngrid)
      REAL,INTENT(IN) :: qsurf(ngrid,nq) ! tracer on surface (kg/m2)

      INTEGER :: ig,nw
      REAL :: min_albedo,max_albedo

c=======================================================================

      ! We get the bare ground albedo from the start files.
      DO ig=1,ngrid
         albedo_bareground(ig)=albedodat(ig)
	 DO nw=1,L_NSPECTV
	    albedo(ig,nw)=albedo_bareground(ig)
	 ENDDO
      ENDDO
      call planetwide_minval(albedo_bareground,min_albedo)
      call planetwide_maxval(albedo_bareground,max_albedo)
      write(*,*) 'surfini: minimum bare ground albedo',min_albedo
      write(*,*) 'surfini: maximum bare ground albedo',max_albedo

    
      call planetwide_minval(albedo,min_albedo)
      call planetwide_maxval(albedo,max_albedo)
      write(*,*) 'surfini: minimum corrected initial albedo',min_albedo
      write(*,*) 'surfini: maximum corrected initial albedo',max_albedo


      END
