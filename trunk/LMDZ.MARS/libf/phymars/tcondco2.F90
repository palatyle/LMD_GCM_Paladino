SUBROUTINE tcondco2(ngrid,nlay,p,q,tcond)
       USE comcstfi_h
       use conc_mod, only: mmean

IMPLICIT NONE 

!---------------------------------------------------
! Condensation temperature for co2 ice; based on 
! the saturation in co2sat.F JA17
!--------------------------------------------------i

integer, intent(in) :: ngrid,nlay
real, intent(in), dimension(ngrid,nlay):: p,q
double precision, intent(out), dimension(ngrid,nlay):: tcond ! CO2 condensation temperature	(atm)
double precision:: A,B,pco2
integer:: ig,l

A=dlog(1.382d12)
B=-3182.48

DO l=1,nlay
   DO ig=1,ngrid
      pco2 = q(ig,l) * (mmean(ig,l)/44.01) * p(ig,l)
      tcond(ig,l)=B/(dlog(pco2)-A)
    enddo
enddo

end
