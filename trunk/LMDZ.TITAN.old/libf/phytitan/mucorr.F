      SUBROUTINE mucorr(npts,lon_sun, plat, pmu, pfract)
      IMPLICIT NONE

c=======================================================================
c
c   Calcul of equivalent solar angle and and fraction of day whithout 
c   diurnal cycle.
c
c   parmeters :
c   -----------
c
c      Input :
c      -------
c         npts             number of points
c         lon_sun          solar longitude (radians!!)
c         plat(npts)       latitude (en degres)
c
c      Output :
c      --------
c         pmu(npts)          equivalent cosinus of the solar angle
c         pfract(npts)       fractionnal day
c
c=======================================================================

c-----------------------------------------------------------------------
c
c    0. Declarations :
c    -----------------

#include "YOMCST.h"
#include "comorbit.h"

c     Arguments :
c     -----------
      INTEGER npts
      REAL plat(npts),pmu(npts), pfract(npts)
      REAL lon_sun
c
c     Local variables :
c     -----------------
      INTEGER j
      REAL z,cz,sz,tz,phi,cphi,sphi,tphi
      REAL ap,a,t,b,tp
      real pdeclin,incl

c-----------------------------------------------------------------------

c verifs
c     print*,"LATITUDES=",plat(1),plat(npts/2),plat(npts)
c     print*,"zls=",lon_sun*180/RPI

      incl=obliquit * RPI / 180.                ! obliquite en radian
      pdeclin = ASIN (SIN(lon_sun)*SIN(incl) ) ! declin en radian
c     print*,'npts,pdeclin',npts,pdeclin*180./RPI
      z = pdeclin
      cz = cos (z)
      sz = sin (z)
c      print*,'cz,sz',cz,sz

      DO 20 j = 1, npts

         phi = plat(j)*RPI/180.  ! latitude en radian
         cphi = cos(phi)
         if (cphi.le.1.e-9) cphi=1.e-9
         sphi = sin(phi)
         tphi = sphi / cphi
         b = cphi * cz
         t = -tphi * sz / cz
         a = 1.0 - t*t
         ap = a

         IF(t.eq.0.) then
            t=0.5*RPI
         ELSE
            IF (a.lt.0.) a = 0.
            t = sqrt(a) / t
            IF (t.lt.0.) then
               tp = -atan (-t) + RPI
            ELSE
               tp = atan(t)
            ENDIF
            t = tp
         ENDIF
   
         pmu(j) = (sphi*sz*t) / RPI + b*sin(t)/RPI
         pfract(j) = t / RPI
         IF (ap .lt.0.) then
            pmu(j) = sphi * sz
            pfract(j) = 1.0
         ENDIF

         IF (pmu(j).le.0.0) pmu(j) = 0.0
         pmu(j) = pmu(j) / pfract(j)
         IF (pmu(j).eq.0.) pfract(j) = 0.

   20 CONTINUE
c        call dump2d(48,31,pfract(2),'FRACT      ')
c        call dump2d(48,31,pmu(2),'MU0        ')
c        stop
                                  
c-----------------------------------------------------------------------
c   correction de rotondite:
c   ------------------------

c        print*,'dans mucorr avant correction rotondite'
c        print*,'pmu(1)=',pmu(1),' pmu(npts/2)=',pmu(npts/2)
c        print*,'pfract(1)=',pfract(1),' pfract(npts/2)=',pfract(npts/2)
         
      DO 30 j=1,npts
c !!!!!!
         pmu(j)=sqrt(1224.*pmu(j)*pmu(j)+1.)/35.
30    CONTINUE

c        print*,'dans mucorr apres correction rotondite'
c        print*,'pmu(1)=',pmu(1),' pmu(npts/2)=',pmu(npts/2)
         
c     print*,"pmu=",pmu(1),pmu(npts/2),pmu(npts)
c     print*,"pfract=",pfract(1),pfract(npts/2),pfract(npts)
      RETURN
      END
