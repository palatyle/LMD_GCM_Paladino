!**********************************************************************

!	diffusion.h
	
!**********************************************************************

      real*8 Pdiff	
      real*8 tdiffmin
      real*8 dzres

      parameter (Pdiff=15.)      ! pressure below which diffusion is computed
      parameter (tdiffmin=5d0)
      parameter (dzres=2d0)	 ! grid resolution (km) for diffusion
      
