      SUBROUTINE spectral_albedo_calc(albedo_snow_SPECTV)


      use callkeys_mod, only : albedosnow
      use radinc_h, only : L_NSPECTV
      use radcommon_h, only: BWNV

      IMPLICIT NONE
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!
!
!   Spectral Albedo Loading ...
!
!!!!!!!!!
!
! - Snow Albedo. Added by MT (10/2015)
!
! What has to be added next : Spectral albedo of CO2 ice, of oceans, ...
!    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      REAL albedo_snow_SPECTV(L_NSPECTV) ! Snow spectral albedo.
      REAL lambda_temp(L_NSPECTV)        ! Visible band wavelengths.
      INTEGER nw


      DO nw=1,L_NSPECTV
      
         ! Corresponding Middle-Band Wavelength Number Calculation.
         lambda_temp(nw)=0.5*(10000.0/BWNV(nw)+10000.0/BWNV(nw+1))
	 
         ! Spectral Snow Albedo Calculation.
	 call albedo_snow_calc(lambda_temp(nw),        &
	                      albedo_snow_SPECTV(nw),  & 
	                      albedosnow)
         
      ENDDO

      RETURN
      END
