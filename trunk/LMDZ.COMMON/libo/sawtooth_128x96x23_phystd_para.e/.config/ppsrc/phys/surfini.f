










      SUBROUTINE surfini(ngrid,nq,qsurf,albedo,albedo_bareground,
     &                   albedo_snow_SPECTV,albedo_co2_ice_SPECTV)

      USE surfdat_h, only: albedodat
      USE tracer_h, only: igcm_co2_ice
      use planetwide_mod, only: planetwide_maxval, planetwide_minval
      use radinc_h, only : L_NSPECTV
      use callkeys_mod, only : albedosnow, albedoco2ice

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
      REAL,INTENT(OUT) :: albedo_snow_SPECTV(L_NSPECTV)
      REAL,INTENT(OUT) :: albedo_co2_ice_SPECTV(L_NSPECTV)
      REAL,INTENT(IN) :: qsurf(ngrid,nq) ! tracer on surface (kg/m2)

      INTEGER :: ig,nw
      REAL :: min_albedo,max_albedo

c=======================================================================

      ! Step 1 : Initialisation of the Spectral Albedos.
      DO nw=1,L_NSPECTV
         albedo_snow_SPECTV(nw)=albedosnow
         albedo_co2_ice_SPECTV(nw)=albedoco2ice
      ENDDO


      ! Step 2 : We get the bare ground albedo from the start files.
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


      ! Step 3 : We modify the albedo considering some CO2 at the surface. We dont take into account water ice (this is made in hydrol after the first timestep) ...
      if (igcm_co2_ice.ne.0) then
         DO ig=1,ngrid
            IF (qsurf(ig,igcm_co2_ice) .GT. 1.) THEN ! This was changed by MT2015. Condition for ~1mm of CO2 ice deposit.
	       DO nw=1,L_NSPECTV
	          albedo(ig,nw)=albedo_co2_ice_SPECTV(nw)
	       ENDDO
            END IF   
         ENDDO   
      else
         write(*,*) "surfini: No CO2 ice tracer on surface  ..."
         write(*,*) "         and therefore no albedo change."
      endif     
      call planetwide_minval(albedo,min_albedo)
      call planetwide_maxval(albedo,max_albedo)
      write(*,*) 'surfini: minimum corrected initial albedo',min_albedo
      write(*,*) 'surfini: maximum corrected initial albedo',max_albedo


      END
