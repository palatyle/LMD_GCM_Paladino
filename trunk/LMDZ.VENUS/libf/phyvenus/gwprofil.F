      subroutine gwprofil
     *         ( nlon, nlev
     *         , kgwd ,kdx  , ktest
     *         , kkcrit, kkcrith, kcrit ,  kkenvh, kknu,kknu2
     *         , kkbreak
     *         , paphm1, prho   , pstab , ptfr , pvph , pri , ptau
     *         , pdmod   , pnu   , psig ,pgamma, pstd, ppic,pval)

C**** *gwprofil*
C
C     purpose.
C     --------
C
C**   interface.
C     ----------
C          from *gwdrag*
C
C        explicit arguments :
C        --------------------
C     ==== inputs ===
C
C     ==== outputs ===
C
C        implicit arguments :   none
C        --------------------
C
C     method:
C     -------
C     the stress profile for gravity waves is computed as follows:
C     it decreases linearly with heights from the ground 
C     to the low-level indicated by kkcrith,
C     to simulates lee waves or 
C     low-level gravity wave breaking.
C     above it is constant, except when the waves encounter a critical
C     level (kcrit) or when they break.
C     The stress is also uniformly distributed above the level
C     ntop.                                          
C
      use dimphy
      IMPLICIT NONE

#include "YOMCST.h"
#include "YOEGWD.h"

C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      integer nlon,nlev,kgwd
      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon)
     *       ,kdx(nlon),ktest(nlon)
     *       ,kkenvh(nlon),kknu(nlon),kknu2(nlon),kkbreak(nlon)
C
      real paphm1(nlon,nlev+1), pstab(nlon,nlev+1),
     *     prho  (nlon,nlev+1), pvph (nlon,nlev+1),
     *     pri   (nlon,nlev+1), ptfr (nlon), ptau(nlon,nlev+1)
     
      real pdmod (nlon) , pnu (nlon) , psig(nlon),
     *     pgamma(nlon) , pstd(nlon) , ppic(nlon), pval(nlon)
     
C-----------------------------------------------------------------------
C
C*       0.2   local arrays
C              ------------
C
      integer jl,jk
      real zsqr,zalfa,zriw,zdel,zb,zalpha,zdz2n,zdelp,zdelpt

      real zdz2 (klon,klev) , znorm(klon) , zoro(klon)
      real ztau (klon,klev+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
C      print *,' entree gwprofil' 
 100  CONTINUE
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
      do 400 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
      zoro(jl)=psig(jl)*pdmod(jl)/4./pstd(jl)
      ztau(jl,klev+1)=ptau(jl,klev+1)
c     print *,jl,ptau(jl,klev+1)
      ztau(jl,kkcrith(jl))=grahilo*ptau(jl,klev+1)
      endif
  400 continue
  
C
      do 430 jk=klev+1,1,-1
C
C
C*         4.1    constant shear stress until top of the
C                 low-level breaking/trapped layer
  410 CONTINUE
C
      do 411 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
           if(jk.gt.kkcrith(jl)) then
           zdelp=paphm1(jl,jk)-paphm1(jl,klev+1) 
           zdelpt=paphm1(jl,kkcrith(jl))-paphm1(jl,klev+1) 
           ptau(jl,jk)=ztau(jl,klev+1)+zdelp/zdelpt*
     c                 (ztau(jl,kkcrith(jl))-ztau(jl,klev+1))
           else                    
           ptau(jl,jk)=ztau(jl,kkcrith(jl))
           endif
       endif
 411  continue             
C
C*         4.15   constant shear stress until the top of the
C                 low level flow layer.
 415  continue
C        
C
C*         4.2    wave displacement at next level.
C
  420 continue
C
  430 continue

C
C*         4.4    wave richardson number, new wave displacement
C*                and stress:  breaking evaluation and critical 
C                 level
C


      do 440 jk=klev,1,-1

      do 441 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
      znorm(jl)=prho(jl,jk)*sqrt(pstab(jl,jk))*pvph(jl,jk)
      zdz2(jl,jk)=ptau(jl,jk)/amax1(znorm(jl),gssec)/zoro(jl)
      endif
  441 continue

      do 442 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
          if(jk.lt.kkcrith(jl)) then
          if((ptau(jl,jk+1).lt.gtsec).or.(jk.le.kcrit(jl))) then
             ptau(jl,jk)=0.0
          else
               zsqr=sqrt(pri(jl,jk))
               zalfa=sqrt(pstab(jl,jk)*zdz2(jl,jk))/pvph(jl,jk)
               zriw=pri(jl,jk)*(1.-zalfa)/(1+zalfa*zsqr)**2
               if(zriw.lt.grcrit) then
c                 print *,' breaking!!!',ptau(jl,jk),zsqr
                  zdel=4./zsqr/grcrit+1./grcrit**2+4./grcrit
                  zb=1./grcrit+2./zsqr
                  zalpha=0.5*(-zb+sqrt(zdel))
                  zdz2n=(pvph(jl,jk)*zalpha)**2/pstab(jl,jk)
                  ptau(jl,jk)=znorm(jl)*zdz2n*zoro(jl)
               endif
                
               ptau(jl,jk)=amin1(ptau(jl,jk),ptau(jl,jk+1))
                  
          endif
          endif
      endif
  442 continue
  440 continue

C  REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

      do 530 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
         ztau(jl,kkcrith(jl)-1)=ptau(jl,kkcrith(jl)-1)
         ztau(jl,ntop)=ptau(jl,ntop)
      endif
 530  continue      

      do 531 jk=1,klev
      
      do 532 jl=kidia,kfdia
      if(ktest(jl).eq.1)then
                
         if(jk.gt.kkcrith(jl)-1)then

          zdelp=paphm1(jl,jk)-paphm1(jl,klev+1    )
          zdelpt=paphm1(jl,kkcrith(jl)-1)-paphm1(jl,klev+1    )
          ptau(jl,jk)=ztau(jl,klev+1    ) +
     .                (ztau(jl,kkcrith(jl)-1)-ztau(jl,klev+1    ) )*
     .                zdelp/zdelpt
     
        endif
      endif
            
 532  continue    
 
C  REORGANISATION AT THE MODEL TOP....

      do 533 jl=kidia,kfdia
      if(ktest(jl).eq.1)then

         if(jk.lt.ntop)then

          zdelp =paphm1(jl,ntop)
          zdelpt=paphm1(jl,jk)
          ptau(jl,jk)=ztau(jl,ntop)*zdelpt/zdelp 
c         ptau(jl,jk)=ztau(jl,ntop)                

        endif

      endif

 533  continue

 
 531  continue        

c Yo, this is Venus.
      do jl=kidia,kfdia
        do jk=klev,1,-1
          if(ktest(jl).eq.1)then
            if(jk.lt.kkbreak(jl)) ptau(jl,jk)=0.0
          endif
        enddo
      enddo

                          
! Venus: resolve waves
      do jk=klev,1,-1
      do jl=kidia,kfdia
      if(ktest(jl).eq.1)then
      ! if surface stress greater than threshold
      if (ztau(jl,klev+1) .ge. taubs) then 
           ! then enforce same stress in the atmosphere up to the predefined level 
           if  ((jk.gt.levbs)) then
             ptau(jl,jk) = ztau(jl,klev+1)
           ! and zero above
           elseif (jk.le.levbs) then
             ptau(jl,jk) = 0.
           endif
!      else
          !if (jk.le.klev-1) ptau(jl,jk) = 0.
!          ptau(jl,jk) = 0.
      endif
      endif
      enddo
      enddo



 123   format(i4,1x,20(f6.3,1x))


      return
      end

