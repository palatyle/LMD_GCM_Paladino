      SUBROUTINE gwstress
     *         (  nlon  , nlev
     *         , kkcrit, ksect, kkhlim, ktest, kkcrith, kcrit, kkenvh
     *         , kknu
     *         , prho  , pstab , pvph  , pstd, psig
     *         , pmea , ppic , pval  , ptfr  , ptau  
     *         , pgeom1 , pgamma , pd1  , pd2   , pdmod , pnu 
     *         , zeff )
c
c**** *gwstress*
c
c     purpose.
c     --------
c  Compute the surface stress due to Gravity Waves, according
c  to the Phillips (1979) theory of 3-D flow above 
c  anisotropic elliptic ridges.

C  The stress is reduced two account for cut-off flow over
C  hill.  The flow only see that part of the ridge located
c  above the blocked layer (see zeff).
c
c**   interface.
c     ----------
c     call *gwstress*  from *gwdrag*
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c     ==== outputs ===
c
c        implicit arguments :   none
c        --------------------
c
c     method.
c     -------
c
c
c     externals.
c     ----------
c
c
c     reference.
c     ----------
c
c   LOTT and MILLER (1997)  &  LOTT (1999)
c
c     author.
c     -------
c
c     modifications.
c     --------------
c     f. lott put the new gwd on ifs      22/11/93
c
c-----------------------------------------------------------------------
      use dimphy
      implicit none

#include "YOMCST.h"
#include "YOEGWD.h"

c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
      integer nlon,nlev
      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon),ksect(nlon),
     *        kkhlim(nlon),ktest(nlon),kkenvh(nlon),kknu(nlon)
c
      real prho(nlon,nlev+1),pstab(nlon,nlev+1),ptau(nlon,nlev+1),
     *     pvph(nlon,nlev+1),ptfr(nlon),
     *     pgeom1(nlon,nlev),pstd(nlon)
c
      real pd1(nlon),pd2(nlon),pnu(nlon),psig(nlon),pgamma(nlon)
      real pmea(nlon),ppic(nlon),pval(nlon)
      real pdmod(nlon)
      real zeff(nlon) ! effective height seen by the flow when there is blocking
c
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------

      integer jl
c
c-----------------------------------------------------------------------
c
c*       0.3   functions
c              ---------
c     ------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
c      PRINT *,' in gwstress'
 100  continue
c
c*         3.1     gravity wave stress.
c
  300 continue
c
c
      zeff = 0.
      do 301 jl=kidia,kfdia
      if(ktest(jl).eq.1) then
      
c  effective mountain height above the blocked flow
  
         zeff(jl)=ppic(jl)-pval(jl)
         if(kkenvh(jl).lt.klev)then
         zeff(jl)=amin1(GFRCRIT*pvph(jl,klev+1)/sqrt(pstab(jl,klev+1))
     c              ,zeff(jl))
         endif

      
        ptau(jl,klev+1)=gkdrag*prho(jl,klev+1)
     *     *psig(jl)*pdmod(jl)/4./pstd(jl)
     *     *pvph(jl,klev+1)*sqrt(pstab(jl,klev+1))
     *     *zeff(jl)**2


c  too small value of stress or  low level flow include critical level
c  or low level flow:  gravity wave stress nul.
                
c       lo=(ptau(jl,klev+1).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
c    *      .or.(pvph(jl,klev+1).lt.gvcrit)
c       if(lo) ptau(jl,klev+1)=0.0
      
c      print *,jl,ptau(jl,klev+1)

      else
      
          ptau(jl,klev+1)=0.0
          
      endif

  301 continue

c      write(21)(ptau(jl,klev+1),jl=kidia,kfdia)
 
      return
      end


