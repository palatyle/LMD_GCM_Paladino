      SUBROUTINE SULW
C     
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "callkeys.h"
#include "yomlw.h"
C     
C     PURPOSE.
C     --------
C     INITIALIZE YOMLW, THE COMMON THAT CONTAINS COEFFICIENTS
C     NEEDED TO RUN THE LONGWAVE RADIATION ROUTINES
C     
C     AUTHOR.
C     -------
C     JEAN-JACQUES MORCRETTE  *ECMWF*
c     Simplified F. Forget 01/2000
C     
C     ------------------------------------------------------------------
C     ----------------------------------------------------------------
C     
C     *       1.    SET VALUES.
C     -----------
C     
C     
C--   REFERENCE TEMPERATURES FOR THE PLANCK FUNCTIONS ----- TSTAND
C     TEMPERATURE DEPENDENCE OF THE ABSORPTION ----- TREF
      DATA TSTAND,TREF / 200.0, 200.0 /
C     
C--   ROOTS AND WEIGHTS FOR THE 2-POINT GAUSSIAN QUADRATURE
c     DATA (RT1(IG1),IG1=1,2) / -0.577350269, +0.577350269 /
c     DATA (WG1(IG1),IG1=1,2) /  1.0        ,  1.0         /

c-----------------------------------------------------------------------
C--   COEFFICIENTS OF THE POLYNOMIALS GIVING THE PLANCK FUNCTIONS
C     
      DATA xp/

     s     0.63849788E+01, 0.30969419E+02, 0.44790835E+02,
     s     0.52651048E+01,-0.18799237E+02, 0.92836181E+01,

     s     0.26166790E+02, 0.12348011E+03, 0.17868306E+03,
     s     0.33657659E+02,-0.66869343E+02, 0.21017507E+02,
c     
c     s  0.47045285E+02, 0.12234377E+03, 0.61873275E+02,
c     s -0.31971883E+02, 0.59168472E+01, 0.91927407E+01,
c     
c     s  0.11101254E+02, 0.86037325E+02, 0.25892695E+03,
c     s  0.35582991E+03, 0.16958020E+03,-0.41311413E+02/
c     
     s     0.11101254E+02, 0.86037325E+02, 0.25892695E+03,
     s     0.35582991E+03, 0.16958020E+03,-0.41311413E+02,

     s     0.47045285E+02, 0.12234377E+03, 0.61873275E+02,
     s     -0.31971883E+02, 0.59168472E+01, 0.91927407E+01/

c     
c-----------------------------------------------------------------------
c     temperature dependency of absorber amounts:
c     -------------------------------------------

      DATA at
     s     / 0.694E-03, 0.272E-02, 0.275E-02, 0.178E-01/
      DATA bt
     s     / 0.328E-05, 0.298E-05,-0.705E-04,-0.163E-04/

c-----------------------------------------------------------------------
      DATA (ga(ie, 1),ie=1,4)
     s     / 0.288231E-04, 0.170794E-01,-0.339714E-01, 0.000000E+00/
      DATA (gb(ie, 1),ie=1,4)
     s     / 0.288231E-04, 0.145426E-01, 0.543812E+00, 0.100000E+01/

      DATA (ga(ie, 2),ie=1,4)
     s     / 0.289299E-01, 0.190634E+01, 0.384061E+01, 0.000000E+00/
      DATA (gb(ie, 2),ie=1,4)
     s     / 0.289299E-01, 0.189485E+01, 0.600363E+01, 0.100000E+01/

      DATA cst_voigt / 0.500E-02, 0.100E-01, 0.150E-01, 0.100E+00/

      return
      end
