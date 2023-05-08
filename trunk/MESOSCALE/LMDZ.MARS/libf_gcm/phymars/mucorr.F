      SUBROUTINE mucorr(npts,pdeclin, plat, pmu, pfract,phaut,prad)
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
c         pdeclin          solar declinaison
c         plat(npts)        latitude
c         phaut            hauteur typique de l'atmosphere
c         prad             rayon planetaire
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

c     Arguments :
c     -----------
      INTEGER npts
      REAL plat(npts),pmu(npts), pfract(npts)
      REAL phaut,prad,pdeclin
c
c     Local variables :
c     -----------------
      INTEGER j
      REAL pi
      REAL z,cz,sz,tz,phi,cphi,sphi,tphi
      REAL ap,a,t,b, tp
      REAL alph

c-----------------------------------------------------------------------

      pi = 4. * atan(1.0)
      z = pdeclin
      cz = cos (z)
      sz = sin (z)

      DO 20 j = 1, npts

         phi = plat(j)
         cphi = cos(phi)
         if (cphi.le.1.e-9) cphi=1.e-9
         sphi = sin(phi)
         tphi = sphi / cphi
         b = cphi * cz
         t = -tphi * sz / cz
         a = 1.0 - t*t
         ap = a
   
         IF(t.eq.0.) then
            tp=0.5*pi
         ELSE
            IF (a.lt.0.) a = 0.
            t = sqrt(a) / t
            IF (t.lt.0.) then
               tp = -atan (-t) + pi
            ELSE
               tp = atan(t)
            ENDIF
         ENDIF
         t = tp
   
         pmu(j) = (sphi*sz*t) / pi + b*sin(t)/pi
         pfract(j) = t / pi
         IF (ap .lt.0.) then
            pmu(j) = sphi * sz
            pfract(j) = 1.0
         ENDIF
         IF (pmu(j).le.0.0) pmu(j) = 0.0
         pmu(j) = pmu(j) / pfract(j)
         IF (pmu(j).eq.0.) pfract(j) = 0.

   20 CONTINUE

c-----------------------------------------------------------------------
c   correction de rotondite:
c   ------------------------

      alph=phaut/prad
      DO 30 j=1,npts
c !!!!!!
         pmu(j)=sqrt(1224.*pmu(j)*pmu(j)+1.)/35.
c    $          (sqrt(alph*alph*pmu(j)*pmu(j)+2.*alph+1.)-alph*pmu(j))
30    CONTINUE

      RETURN
      END
