SUBROUTINE tcondwater(nbpts,p,q,tcond)
IMPLICIT NONE 

!---------------------------------------------------
! Condensation temperature for water ice; based on 
! the saturation in watersat.F AP14
!--------------------------------------------------i

integer, intent(in) :: nbpts
real, intent(in), dimension(nbpts):: p,q
real, intent(out), dimension(nbpts):: tcond

real, dimension(nbpts):: res
real:: x
integer:: i

do i=1,nbpts
   !write(*,*) "i en cours", i, "sur nbpts=", nbpts
   !write(*,*) "q(i)",q(i),"p(i)",p(i)
   x=alog10(MAX(q(i),1e-16)*p(i)/(100.*0.41)) ! max pour erreur q<=0
   ! attention change le 0.41 de place le 10 juin 2014 car priorites
   ! fortran dans watersat.F
   !write(*,*) "x tcondwater AP14! :) :) :) :P", x
   !res(i) = 2.52826991e+02+ 2.39287870e+01*x+ 2.27275932e+00*x**2
   !        + 2.21832905e-01*x**3+ 2.23453930e-02*x**4+2.26075106e-03*x**5
   !        + 2.12411064e-04*x**6+1.64642075e-05*x**7+9.22615632e-07*x**8
   !        + 3.18958825e-08*x**9+5.00656720e-10*x**10 degre 10: trop!
   res(i) = 2.52846556e+02+ 2.39229653e+01*x+ 2.21333897e+00*x**2  &
                + 1.79977992e-01*x**3+ 1.00068175e-02*x**4+2.55145012e-04*x**5
   !write(*,*) "rex(x) tcondwater AP14! :) :) :) :P", res(i)
enddo

tcond=res

return 
end
!polynome de degre 5 pas 0.0001
!polynomial coefs [  2.52846556e+02   2.39229653e+01   2.21333897e+00
!1.79977992e-01
!   1.00068175e-02   2.55145012e-04]
!maximum des abs(difference) 0.0604390646026

! polynome de degre 12
!polynomial coefs [  2.52826992e+02   2.39287716e+01   2.27274870e+00
!2.21863471e-01
 !  2.23765903e-02   2.26605393e-03   2.07841624e-04   1.37374700e-05
 !  2.45106231e-07  -6.16151111e-08  -6.96026651e-09  -3.22690558e-10
!  -5.86804217e-12]   maximum des abs(difference) 2.73827428146e-05
! polynome de degre 6
!polynomial coefs [  2.52831053e+02   2.39333049e+01   2.26006967e+00
!2.06350715e-01
!   1.56882616e-02   7.83034223e-04   1.77637297e-05] 
!maximum des abs(difference) 0.013723768413
! polynome de degre 8
! polynomial coefs [  2.52827042e+02   2.39294477e+01   2.27281607e+00
!2.20637577e-01
!   2.12491384e-02   1.82810482e-03   1.19638702e-04   4.91167244e-06
!   9.12272144e-08]
!maximum des abs(difference) 0.000725035990172
!polynome de degre 9
!polynomial coefs [  2.52826985e+02   2.39289186e+01   2.27286626e+00
!2.21623591e-01
!   2.20538443e-02   2.11335665e-03   1.73307328e-04   1.04920165e-05
!   3.94201385e-07   6.70874574e-09]
!maximum des abs(difference) 0.000168806463876
!polynome de degre 10
!polynomial coefs [  2.52826991e+02   2.39287870e+01   2.27275932e+00
!2.21832905e-01
!   2.23453930e-02   2.26075106e-03   2.12411064e-04   1.64642075e-05
!   9.22615632e-07   3.18958825e-08   5.00656720e-10]
!maximum des abs(difference) 3.96286844477e-05
