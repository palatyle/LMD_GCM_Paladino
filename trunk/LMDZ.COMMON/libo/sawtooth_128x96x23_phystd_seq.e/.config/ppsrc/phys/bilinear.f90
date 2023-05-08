










!-------------------------------------------------------------------------
subroutine bilinear(f,f11,f21,f12,f22,x,x1,x2,y,y1,y2)  
! Used for interpolation of continuum data

  implicit none

  real*8 x,y,x1,x2,y1,y2
  real*8 f,f11,f12,f21,f22,fA,fB

  ! 1st in x-direction
  fA=f11*(x2-x)/(x2-x1)+f21*(x-x1)/(x2-x1)
  fB=f12*(x2-x)/(x2-x1)+f22*(x-x1)/(x2-x1)

  ! then in y-direction
  f=fA*(y2-y)/(y2-y1)+fB*(y-y1)/(y2-y1)

  return
end subroutine bilinear
