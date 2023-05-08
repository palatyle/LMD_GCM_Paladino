      SUBROUTINE co2sat(naersize,t,p,psat)
c       SUBROUTINE co2sat(naersize,t,p,qsat) JA
      IMPLICIT NONE

c=======================================================================
c
c
c  now:  straight psat of CO2 (or qsat of CO2 but need of mmean)
c
c=======================================================================

c   declarations:
c   -------------
c   arguments:
c   ----------

c   INPUT
      integer naersize 
      real t(naersize) , p(naersize)
c   OUTPUT
c      real qsat(naersize) JA
      real psat(naersize)

c   local:
c   ------
      INTEGER i
      REAL r2,r3,r4 , To, es
      SAVE r2,r3,r4
      DATA r2,r3,r4/611.14,21.875,7.66/
      SAVE To
      DATA To/273.16/
          
      do i=1,naersize


c        pression de vapeur saturante (James et al. 1992):

          psat(i)  = 1.382 * 1e12 * exp(-3182.48/t(i)) !; (Pa)

c         OR:

c         qsat(i) = psat/p(i)*44.01/mmean ! Need of updated information on mmean
c         qsat(i) = max(qsat(i), 1.e-30)


      enddo
c      qsat=psat JA
          

      RETURN
      END

