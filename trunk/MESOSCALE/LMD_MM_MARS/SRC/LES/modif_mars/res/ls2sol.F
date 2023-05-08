cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function ls2sol(ls)

c  Returns solar longitude, Ls (in deg.), from day number (in sol),
c  where sol=0=Ls=0 at the northern hemisphere spring equinox

      implicit none

c  Arguments:
      real ls

c  Local:
      double precision xref,zx0,zteta,zz
c	xref: mean anomaly, zteta: true anomaly, zx0: eccentric anomaly
      double precision year_day 
      double precision peri_day,timeperi,e_elips
      double precision pi,degrad 
      parameter (year_day=668.6d0) ! number of sols in a amartian year
c      data peri_day /485.0/
      parameter (peri_day=485.35d0) ! date (in sols) of perihelion
c  timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
      parameter (timeperi=1.90258341759902d0)
      parameter (e_elips=0.0934d0)  ! eccentricity of orbit
      parameter (pi=3.14159265358979d0)
      parameter (degrad=57.2957795130823d0)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            ls2sol = 0.0
         else
            ls2sol = year_day
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      ls2sol = zz*year_day + peri_day
      if (ls2sol.lt.0.0) ls2sol = ls2sol + year_day
      if (ls2sol.ge.year_day) ls2sol = ls2sol - year_day

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
