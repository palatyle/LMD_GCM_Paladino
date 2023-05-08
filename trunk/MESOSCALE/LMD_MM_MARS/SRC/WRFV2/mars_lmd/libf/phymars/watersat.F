      SUBROUTINE watersat(nsize,t,p,qsat)
      IMPLICIT NONE

c=======================================================================
c
c  Water mass mixing ratio at saturation (kg/kg) for a given pressure (Pa)
c   and Temperature (K) array
c
c=======================================================================

c   declarations:
c   -------------
c   arguments:
c   ----------

c   INPUT
      integer nsize 
      real t(nsize) , p(nsize)
c   OUTPUT
      real qsat(nsize)

c   local:
c   ------
      INTEGER i
      REAL r2,r3,r4 , To, es
      SAVE r2,r3,r4
      DATA r2,r3,r4/611.14,21.875,7.66/
      SAVE To
      DATA To/273.16/
          
      do i=1,nsize
c        pression de vapeur saturante :
c         es =r2*exp(r3*(t(i)-To)/(t(i)-r4))
c         qsat(i)=es/p(i)*0.41  !rapport de melange
        es=2.07023-0.00320991*t(i)-2484.896/t(i)+3.56654*alog10(t(i))
        es=10.**es 
        qsat(i)=(es*100.)/p(i)*0.41 

c        qss = specific humidity
c        qss=622*10**(2.07023-0.0032991*Tsurf(ig)-2484.896/Tsurf(ig)
c    $                  +3.56654*log10(Tsurf(ig)))

      enddo

      RETURN
      END

