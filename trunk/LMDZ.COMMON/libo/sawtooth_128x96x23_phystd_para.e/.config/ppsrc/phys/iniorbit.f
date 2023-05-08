










      SUBROUTINE iniorbit
     $     (papoastr,pperiastr,pyear_day,pperi_day,pobliq)
      
      USE planete_mod, only: apoastr, periastr, year_day, obliquit,
     &                       peri_day, e_elips, p_elips, timeperi
      use comcstfi_mod, only: pi
      IMPLICIT NONE

!=======================================================================
! Initialisation of orbital parameters (stored in planete_h module)
!=======================================================================

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: papoastr,pperiastr,pyear_day,pperi_day,pobliq

c   Local:
c   ------

      REAL zxref,zanom,zz,zx0,zdx
      INTEGER iter

c-----------------------------------------------------------------------

      pi=2.*asin(1.)

      apoastr =papoastr
      periastr=pperiastr
      year_day=pyear_day
      obliquit=pobliq
      peri_day=pperi_day

      PRINT*,'iniorbit: Periastron in AU  ',periastr
      PRINT*,'iniorbit: Apoastron in AU  ',apoastr 
      PRINT*,'iniorbit: Obliquity in degrees  :',obliquit


      e_elips=(apoastr-periastr)/(periastr+apoastr)
      p_elips=0.5*(periastr+apoastr)*(1-e_elips*e_elips)

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
