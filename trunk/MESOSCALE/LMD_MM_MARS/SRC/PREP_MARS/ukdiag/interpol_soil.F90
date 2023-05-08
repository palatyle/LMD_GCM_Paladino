! ====================================================================
!
! This subroutine converts arrays written on regolith levels to arrays
! on soil levels (sometimes different layers had to be specified for
! the regolith scheme (compared to the soil layers) to ensure numerical
! stability at high obliquities)
!
! Liam Steele March 2015
!
! ====================================================================

subroutine interpol_soil(mlayer,mlayer_r,nsoil,nregmx,id,old_array,new_array)
implicit none

integer :: ik, jk, it, ii, ij
integer :: nsoil, nregmx
real :: lay1, alpha
integer, dimension(4) :: id
real, dimension(nsoil) :: mlayer(0:nsoil-1), layer(1:nsoil)
real, dimension(nreg) :: mlayer_r(0:nregmx-1), layer_r(1:nregmx)
real, dimension(id(1),id(2),id(3),id(4)) :: old_array
real, dimension(id(1),id(2),nsoil,id(4)) :: new_array
integer rloop
real :: lb, ub, lb_r, ub_r
real :: frac(nsoil,nregmx)
integer :: frac_layers(nsoil,nregmx)
integer :: frac_count(nsoil)
integer :: bl(nregmx)

! Calculate soil lower layer boundaries
lay1 = sqrt(mlayer(0)*mlayer(1))
alpha = mlayer(1)/mlayer(0)
do jk=1,nsoil
  layer(jk) = lay1*(alpha**(jk-1))
enddo

! Calculate regolith lower layer boundaries
lay1 = 0.5e-2
if (nregmx.eq.20) then
  alpha = 2.2
else if (nregmx.eq.30) then
  alpha = 1.44
else if (nregmx.eq.40) then
  alpha = 1.26
else if (nregmx.eq.50) then
  alpha = 1.18
else
  print*, 'nregmx is not equal to a specified value'
  stop
endif
layer_r(1) = lay1
do jk=2,10
  layer_r(jk) =  layer_r(jk-1) + lay1
enddo
do jk=11,nregmx
  layer_r(jk) = layer_r(jk-1) + alpha*(layer_r(jk-1)-layer_r(jk-2))
enddo

! Find which soil layer each regolith layer lies within
bl(:) = 0
do ik=0,nregmx-1
  do jk=0,nsoil-1
    if (mlayer_r(ik).ge.mlayer(jk)) bl(ik+1) = jk+1
  enddo
enddo

! Calculate which regolith layers overlap each soil layer, and their contributing fraction
do jk=1,nsoil
  lb = layer(jk)
  if (jk.eq.1) ub = 0.0
  if (jk.gt.1) ub = layer(jk-1)
  frac_count(jk) = 0
  
  do ik=1,nregmx
    lb_r = layer_r(ik)
    if (ik.eq.1) ub_r = 0.0
    if (ik.gt.1) ub_r = layer_r(ik-1)
    
    if (ub.ge.ub_r .and. lb.le.lb_r) then
      ! Whole soil layer is inside regolith layer
      frac_count(jk) = frac_count(jk) + 1
      frac(jk,frac_count(jk)) = 1.0
      frac_layers(jk,frac_count(jk)) = ik
    else if (ub_r.ge.ub .and. lb_r.le.lb) then
      ! Whole regolith layer is inside soil layer
      frac_count(jk) = frac_count(jk) + 1
      frac(jk,frac_count(jk)) = (lb_r-ub_r)/(lb-ub)
      frac_layers(jk,frac_count(jk)) = ik
    else
      if (lb_r.gt.ub .and. ub_r.le.ub) then
        ! Upper overlap
        frac_count(jk) = frac_count(jk) + 1
        frac(jk,frac_count(jk)) = (lb_r-ub)/(lb-ub)
        frac_layers(jk,frac_count(jk)) = ik
      else if (lb_r.gt.lb .and. ub_r.le.lb) then
        ! Lower overlap
        frac_count(jk) = frac_count(jk) + 1
        frac(jk,frac_count(jk)) = (lb-ub_r)/(lb-ub)
        frac_layers(jk,frac_count(jk)) = ik
      endif
    endif
  enddo
enddo

! Convert array from regolith layers to soil layers
new_array(:,:,:,:) = 0.0
do ik=0,nsoil-1
  do rloop=1,frac_count(ik+1)
    new_array(:,:,ik+1,:) = new_array(:,:,ik+1,:) + old_array(:,:,frac_layers(ik+1,rloop),:)*frac(ik+1,rloop)/sum(frac(ik+1,1:frac_count(ik+1)))
  enddo
enddo

return
end
