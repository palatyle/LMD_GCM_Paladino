      SUBROUTINE mufract(jjm,pdecli, plat, pmu, pfract)
      IMPLICIT NONE
c
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
c         jjm              number of points
c         pdecli           solar declinaison
c         plat(jjm)        latitude
c
c      Output :
c      --------
c         pmu(jjm)          equivalent cosinus of the solar angle
c         pfract(jjm)       fractionnal day
c
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0. Declarations :
c    -----------------
c
c     Arguments :
c     -----------

      INTEGER jjm
      REAL plat(jjm),pmu(jjm), pfract(jjm)
      REAL pdecli
c
c     Local variables :
c     -----------------

      INTEGER j
      REAL pi
      REAL z,cz,sz,tz,phi,cphi,sphi,tphi
      REAL ap,a,t,b
c
c=======================================================================
c
      pi = 4. * atan(1.0)
      z = pdecli
      cz = cos (z*pi/180.)
      sz = sin (z*pi/180.)
c
      DO 20 j = 1, jjm
c
         phi = plat(j)
         cphi = cos(phi)
         if (cphi.le.1.e-9) cphi=1.e-9
         sphi = sin(phi)
         tphi = sphi / cphi
         b = cphi * cz
         t = -tphi * sz / cz
         a = 1.0 - t*t
         ap = a
         IF(t.eq.0.) THEN
            t=0.5*pi
         ELSE
            IF (a.lt.0.) a = 0.
            t = sqrt(a) / t
            IF (t.lt.0.) THEN
               t = -atan (-t) + pi
            ELSE
               t = atan(t)
            ENDIF
         ENDIF
         pmu(j) = (sphi*sz*t) / pi + b*sin(t)/pi
         pfract(j) = t / pi
         IF (ap .lt.0.) THEN
            pmu(j) = sphi * sz
            pfract(j) = 1.0
         ENDIF
         IF (pmu(j).le.0.0) pmu(j) = 0.0
         pmu(j) = pmu(j) / pfract(j)
         IF (pmu(j).eq.0.) pfract(j) = 0.
c
   20 CONTINUE
c
      RETURN
      END
c
c* end of mufract
c=======================================================================
c
