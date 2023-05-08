c ********************************************************
      subroutine call_dayperi(Lsperi,e_elips,dayperi,year_day)
      IMPLICIT NONE

c=======================================================================
c
c   Objet:
c   ------
c   Computing the Martian date (number of sols since Ls=0 = spring
c   equinox) at perihelion
c
c   Arguments:
c   ----------
c
c   Input:
c   ------
c   Lsperi      solar longitude (Ls) of perohelion (rad)
c   e_elips       Excentricity
c
c   output
c   ------
c   dayperi       Martian date at perihelion (sol)
c   real year_day      ! number of sols per Mars yar
c-----------------------------------------------------------------------

c arguments:
c ----------
       real Lsperi,e_elips,dayperi
       real year_day      ! number of sols per Mars yar


c      Local
       real x1, x2,pi


       pi=2.*asin(1.)

       x1 = sqrt( (1-e_elips)/(1+e_elips) )
       x2 = e_elips*sqrt(1-e_elips**2)

       dayperi = 0.5*(year_day/pi)*
     &   ( 2*atan(x1*tan(0.5*Lsperi))
     &     -x2*sin(Lsperi)/(1+e_elips*cos(Lsperi)) )
       if(dayperi.lt.0) dayperi=dayperi+year_day
       return
       end
