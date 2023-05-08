      SUBROUTINE swrayleigh(kdlon,knu,ppsol,prmu,prayl)
                                                   
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing total rayleigh scat atmospheric optical depth  
c
c   author: F.Forget 
c   ------
c
c   input:
c   ----- 
c   kdlon             Number of gridpoint of horizontal grid
c    knu   :   Solar band # (1 or 2)
c   ppsol             surface pressure (Pa) 
c   prmu         cos of solar zenith angle (=1 when sun at zenith)
c           (CORRECTED for high zenith angle (atmosphere), unlike mu0)
c
c   output:
c   -------
c   prayl       column optical depth in each model column
c
c=======================================================================
#include "comcstfi.h"

c-----------------------------------------------------------------------
c
c    Declarations :
c    --------------
c
c    Input/Output
c    ------------
      INTEGER kdlon, knu

      real ppsol(kdlon),prmu(kdlon),prayl(kdlon)
c
c    Local variables :
c    -----------------
      integer JL, K
c
c   local saved variables
c   ---------------------
c     rayleigh scattering coefficients (from Morcrete et al.EARTH model !)
      real cray(2,6)

      DATA (CRAY(1,K),K=1,6) /
     S     .428937E-01, .890743E+00,-.288555E+01,
     S     .522744E+01,-.469173E+01, .161645E+01/

      DATA (CRAY(2,K),K=1,6) /
     S     .697200E-02, .173297E-01,-.850903E-01,
     S     .248261E+00,-.302031E+00, .129662E+00/
      save cray
c----------------------------------------------------------------------

      DO JL = 1 , KDLON

c        Total Rayleigh Optical thickness on Earth at 101325 Pa level
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        WARNING : the CRAY coefficients are only valid for
c        Spectral interval used in Earth model !!!
c        (i.e. 0.25-0.68 micron and 0.68-4.00 micron
c         ---> should be modified for Mars model !!!

         PRAYL(JL) = CRAY(KNU,1) + PRMU(JL) * (CRAY(KNU,2) + PRMU(JL)
     S        * (CRAY(KNU,3) + PRMU(JL) * (CRAY(KNU,4) + PRMU(JL)
     S        * (CRAY(KNU,5) + PRMU(JL) *   CRAY(KNU,6)       ))))

c        Total local Rayleigh Optical thickness on Mars
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           - PRAYL above is the total Rayleigh Optical thickness
c             at  101325Pa for Earth atmosphere
c           -> Total Rayleigh Optical thickness on Mars
c              = PRAYL*(Psurf/101325)*9.81/g *2.5
c              (Extinction coeff of CO2 = 2.5 * N2)

c     (Comment the following line to get back bugged version before 01/2000)
c        PRAYL(JL) = PRAYL(JL) * PPSOL(JL) * 2.42e-4/g

      END DO


      return
      end 

