SUBROUTINE print_debug_phys (i,debug_lev,text)

use dimphy
use phys_local_var_mod
use phys_state_var_mod
IMPLICIT NONE
integer i,debug_lev
CHARACTER*(*) text


integer k

print*,'PLANTAGE POUR LE POINT i=',i,text
print*,'l    u, v, T, q, ql'
DO k = 1, klev
   write(*,'(i3,2f8.4,3f14.4,2e14.2)') k,rlon(i),rlat(i),u_seri(i,k),v_seri(i,k),t_seri(i,k),q_seri(i,k),ql_seri(i,k)
ENDDO

RETURN
END
