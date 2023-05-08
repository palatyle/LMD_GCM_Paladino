      SUBROUTINE vdif_cd( ngrid,nlay,pz0,pg,pz,pu,pv,pts,ph,pcdv,pcdh)
      IMPLICIT NONE
c=======================================================================
c
c   Subject: computation of the surface drag coefficient using the
c   -------  approch developed by Loui for ECMWF.
c
c   Author: Frederic Hourdin  15 /10 /93
c   -------
c
c   Arguments:
c   ----------
c
c   inputs:
c   ------
c     ngrid            size of the horizontal grid
c     pg               gravity (m s -2)
c     pz(ngrid)        height of the first atmospheric layer
c     pu(ngrid)        u component of the wind in that layer
c     pv(ngrid)        v component of the wind in that layer
c     pts(ngrid)       surfacte temperature
c     ph(ngrid)        potential temperature T*(p/ps)^kappa
c
c   outputs:
c   --------
c     pcdv(ngrid)      Cd for the wind
c     pcdh(ngrid)      Cd for potential temperature
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

c   Arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL pz0
      REAL pg,pz(ngrid,nlay)
      REAL pu(ngrid,nlay),pv(ngrid,nlay)
      REAL pts(ngrid,nlay),ph(ngrid,nlay)
      REAL pcdv(ngrid),pcdh(ngrid)

c   Local:
c   ------

      INTEGER ig

      REAL zu2,z1,zri,zcd0,zz

      REAL karman,b,c,d,c2b,c3bc,c3b,umin2
      LOGICAL firstcal
      DATA karman,b,c,d,umin2/.4,5.,5.,5.,1.e-12/
      DATA firstcal/.true./
      SAVE b,c,d,karman,c2b,c3bc,c3b,firstcal,umin2
!$OMP THREADPRIVATE(b,c,d,karman,c2b,c3bc,c3b,firstcal,umin2)

c-----------------------------------------------------------------------
c   couche de surface:
c   ------------------

! simplified calculation

      DO ig=1,ngrid
         z1=1.E+0 + pz(ig,1)/pz0
         zcd0=karman/log(z1)
         zcd0=zcd0*zcd0
         pcdv(ig)=zcd0
         pcdh(ig)=zcd0
      ENDDO

c-----------------------------------------------------------------------

      RETURN
      END
