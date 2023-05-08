      SUBROUTINE iniorbit
     $     (paphelie,pperiheli,pyear_day,pperi_day,pobliq)
      use planete_h, only: aphelie, periheli, year_day, peri_day,
     &                     obliquit, unitastr, e_elips, p_elips,
     &                     timeperi
      use comcstfi_h, only: pi
      IMPLICIT NONE

!=======================================================================
! Initialisation of orbital parameters (stored in planete_h module)
!=======================================================================

!   Arguments:
!   ----------

      REAL,INTENT(IN) :: paphelie,pperiheli,pyear_day,pperi_day,pobliq

!   Local:
!   ------

      REAL zxref,zanom,zz,zx0,zdx
      INTEGER iter

!-----------------------------------------------------------------------

      pi=2.*asin(1.)

      aphelie =paphelie
      periheli=pperiheli
      year_day=pyear_day
      obliquit=pobliq
      peri_day=pperi_day

      PRINT*,'iniorbit: Perihelion in Mkm  ',periheli
      PRINT*,'iniorbit: Aphelion  in Mkm  ',aphelie 
      PRINT*,'iniorbit: Obliquity in degrees  :',obliquit
      unitastr=149.597927 ! 1 UA, in Mkm
      e_elips=(aphelie-periheli)/(periheli+aphelie)
      p_elips=0.5*(periheli+aphelie)*(1-e_elips*e_elips)/unitastr

      print*,'iniorbit: e_elips',e_elips
      print*,'iniorbit: p_elips',p_elips

!-----------------------------------------------------------------------
! compute polar angle and distance to the Sun:
! -------------------------------------------------------

!  compute mean anomaly zanom

      zz=(year_day-pperi_day)/year_day
      zanom=2.*pi*(zz-nint(zz))
      zxref=abs(zanom)
      PRINT*,'iniorbit: zanom  ',zanom

!  solve equation  zx0 - e * sin (zx0) = zxref for eccentric anomaly zx0
!  using Newton method

      zx0=zxref+e_elips*sin(zxref)
      DO iter=1,100
         zdx=-(zx0-e_elips*sin(zx0)-zxref)/(1.-e_elips*cos(zx0))
         if(abs(zdx).le.(1.e-12)) exit
         zx0=zx0+zdx
      ENDDO

      zx0=zx0+zdx
      if(zanom.lt.0.) zx0=-zx0
      PRINT*,'iniorbit: zx0   ',zx0

      timeperi=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
      PRINT*,'iniorbit: Perihelion solar long. Ls (deg)=',
     &       360.-timeperi*180./pi

      END
