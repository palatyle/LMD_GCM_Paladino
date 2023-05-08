C*    Coefficients for the longwave radiation subroutines
c-----------------------------------------------------------------------

      COMMON/yomlw/at(2,2),bt(2,2)
     S ,tref,xp(6,nir),tstand
     S ,ga(npademx,nabsmx),gb(npademx,nabsmx)
     S ,cst_voigt(2,nabsmx), gcp , nlaylte

      REAL at,bt, tref,xp,tstand,ga,gb, cst_voigt
      REAL gcp    ! = g/cpp

c     Number of layers on which LTE calculations (in lw and sw) are performed 
c     (Computed in nlthermeq) :
      INTEGER nlaylte       


      COMMON/ksicom/
     S xi(ngridmx,nuco2,0:nflev+1,0:nflev+1)
     S ,xi_ground(ngridmx,nuco2),xi_emis(ngridmx,nuco2,nflev-1)

      real xi, xi_ground, xi_emis
