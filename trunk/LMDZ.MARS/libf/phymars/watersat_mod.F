      MODULE watersat_mod
      
      IMPLICIT NONE
      
      CONTAINS 
      
      SUBROUTINE watersat(naersize,t,p,qsat)
      
      IMPLICIT NONE

c=======================================================================
c
c  Water mass mixing ratio at saturation (kg/kg) for a given pressure (Pa)
c   and Temperature (K) array
c
c  Move in module by Deborah BARDET 02/07/18
c=======================================================================

c   declarations:
c   -------------
c   arguments:
c   ----------

c   INPUT
      integer, intent(in) :: naersize 
      real, intent(in) :: t(naersize)
      real, intent(in) :: p(naersize)
c   OUTPUT
      real, intent(out) :: qsat(naersize)

c   local:
c   ------
      INTEGER i
      REAL,PARAMETER :: r2 = 611.14
      REAL,PARAMETER :: r3 = 21.875
      REAL,PARAMETER :: r4 = 7.66
      REAL,PARAMETER :: To = 273.16
      REAL es
      
      
      do i=1,naersize
c        pression de vapeur saturante :
c         es =r2*exp(r3*(t(i)-To)/(t(i)-r4))
c         qsat(i)=es/p(i)*0.41  !rapport de melange
        es=2.07023-0.00320991*t(i)-2484.896/t(i)+3.56654*alog10(t(i))
        es=10.**es 
        qsat(i)=(es*100.)/p(i)*0.41 
        qsat(i) = max(qsat(i), 1.e-30)

c        qss = specific humidity
c        qss=622*10**(2.07023-0.0032991*Tsurf(ig)-2484.896/Tsurf(ig)
c    $                  +3.56654*log10(Tsurf(ig)))

      enddo

      END SUBROUTINE watersat
      
      END MODULE watersat_mod
