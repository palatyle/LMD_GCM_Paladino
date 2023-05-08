!**********************************************************************

!	diffusion.h
	
!**********************************************************************

      real*8 Pdiff	
      real*8 tdiffmin
      real*8 dzres

      parameter (Pdiff=1.e-1)      ! pressure below which diffusion is computed
      parameter (tdiffmin=5.d0)
      parameter (dzres=0.1d0)	 ! grid resolution (km) for diffusion
      
