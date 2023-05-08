SUBROUTINE fill_save_2d( &
tsurf,&
fluxsurf_sw_tot,&
fluxsurf_lw,&
ust,&
hfx,&
ngrid,&
output_tab2d)
   
IMPLICIT NONE
   
include "wrf_output_2d.h"
   
INTEGER :: ngrid
REAL, DIMENSION(ngrid) :: tsurf
REAL, DIMENSION(ngrid) :: fluxsurf_sw_tot
REAL, DIMENSION(ngrid) :: fluxsurf_lw
REAL, DIMENSION(ngrid) :: ust
REAL, DIMENSION(ngrid) :: hfx
REAL, DIMENSION(ngrid,n2d) :: output_tab2d
   
output_tab2d(:,ind_TSURF)=tsurf(:)
output_tab2d(:,ind_SWDOWN)=fluxsurf_sw_tot(:)
output_tab2d(:,ind_LWDOWN)=fluxsurf_lw(:)
output_tab2d(:,ind_UST)=ust(:)
output_tab2d(:,ind_HFX)=hfx(:)
   
END SUBROUTINE fill_save_2d
