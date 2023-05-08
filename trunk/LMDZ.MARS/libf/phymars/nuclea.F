*******************************************************
*                                                     *
      subroutine nuclea(ph2o,temp,sat,n_ccn,nucrate)
      USE comcstfi_h
      implicit none
*                                                     *
*   This subroutine computes the nucleation rate      *
*   as given in Pruppacher & Klett (1978) in the      *
*   case of water ice forming on a solid substrate.   *
*     Definition refined by Keese (jgr,1989)          *
*   Authors: F. Montmessin                            *
*     Adapted for the LMD/GCM by J.-B. Madeleine      *
*     (October 2011)                                  *
*     Optimisation by A. Spiga (February 2012)        *  
*******************************************************

#include "microphys.h"

c     Inputs
      DOUBLE PRECISION ph2o,sat
      DOUBLE PRECISION n_ccn(nbin_cld)
      REAL temp

c     Output
   !   DOUBLE PRECISION nucrate(nbin_cld)
      REAL nucrate(nbin_cld)

c     Local variables
      DOUBLE PRECISION nh2o
      DOUBLE PRECISION sig      ! Water-ice/air surface tension  (N.m)
      external sig
      DOUBLE PRECISION rstar    ! Radius of the critical germ (m)
      DOUBLE PRECISION gstar    ! # of molecules forming a critical embryo
      DOUBLE PRECISION fistar   ! Activation energy required to form a critical embryo (J)
!      DOUBLE PRECISION zeldov   ! Zeldovitch factor (no dim)
      DOUBLE PRECISION fshape   ! function defined at the end of the file
      DOUBLE PRECISION deltaf

c     Ratio rstar/radius of the nucleating dust particle
c     double precision xratio
      
      double precision mtetalocal ! local mteta in double precision

      double precision fshapesimple,zefshape


      integer i
      
      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

c     *************************************************

      mtetalocal = mteta  !! use mtetalocal for better performance

cccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc ESSAIS TN MTETA = F (T) cccccccccccccc
c      if (temp .gt. 200) then
c         mtetalocal = mtetalocal
c      else if (temp .lt. 190) then
c         mtetalocal = mtetalocal-0.05
c      else
c         mtetalocal = mtetalocal - (190-temp)*0.005
c      endif
c----------------exp law, see Trainer 2008, J. Phys. Chem. C 2009, 113, 2036\u20132040
       !mtetalocal = max(mtetalocal - 6005*exp(-0.065*temp),0.1)
       !mtetalocal = max(mtetalocal - 6005*exp(-0.068*temp),0.1)
               !print*, mtetalocal, temp
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc 
      IF (firstcall) THEN
          print*, ' '  
          print*, 'dear user, please keep in mind that'
          print*, 'contact parameter IS constant'
          !print*, 'contact parameter IS NOT constant:'
          !print*, 'max(mteta - 6005*exp(-0.065*temp),0.1)'
          !print*, 'max(mteta - 6005*exp(-0.068*temp),0.1)'
          print*, ' '  
         firstcall=.false.
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
   

      if (sat .gt. 1.) then    ! minimum condition to activate nucleation

        nh2o   = ph2o / kbz / temp
        rstar  = 2. * sig(temp) * vo1 / (rgp*temp*dlog(sat))
        gstar  = 4. * nav * pi * (rstar * rstar * rstar) / (3.*vo1)
        
        fshapesimple = (2.+mtetalocal)*(1.-mtetalocal)*(1.-mtetalocal)
     &                   / 4.

c       Loop over size bins
        do 200 i=1,nbin_cld

          if ( n_ccn(i) .lt. 1e-10 ) then
c           no dust, no need to compute nucleation!
            nucrate(i)=0.
            goto 200
          endif

          if (rad_cld(i).gt.3000.*rstar) then
            zefshape = fshapesimple
          else
            zefshape = fshape(mtetalocal,rad_cld(i)/rstar)
          endif

          fistar = (4./3.*pi) * sig(temp) * (rstar * rstar) * 
     &             zefshape
          deltaf = (2.*desorp-surfdif-fistar)/
     &             (kbz*temp)
          deltaf = min( max(deltaf, -100.d0), 100.d0)

          if (deltaf.eq.-100.) then
            nucrate(i) = 0.
          else
            nucrate(i)= real(sqrt ( fistar /
     &               (3.*pi*kbz*temp*(gstar*gstar)) )
     &                  * kbz * temp * rstar
     &                  * rstar * 4. * pi
     &                  * ( nh2o*rad_cld(i) )
     &                  * ( nh2o*rad_cld(i) )
     &                  / ( zefshape * nus * m0 )
     &                  * dexp (deltaf))
          endif

200     continue

      else

        do i=1,nbin_cld
          nucrate(i) = 0.
        enddo

      endif

      return
      end

*********************************************************
      double precision function fshape(cost,rap)
      implicit none
*        function computing the f(m,x) factor           *
* related to energy required to form a critical embryo  *
*********************************************************

      double precision cost,rap
      double precision yeah

          !! PHI
          yeah = sqrt( 1. - 2.*cost*rap + rap*rap )
          !! FSHAPE = TERM A
          fshape = (1.-cost*rap) / yeah
          fshape = fshape * fshape * fshape
          fshape = 1. + fshape
          !! ... + TERM B
          yeah = (rap-cost)/yeah
          fshape = fshape + 
     & rap*rap*rap*(2.-3.*yeah+yeah*yeah*yeah)
          !! ... + TERM C 
          fshape = fshape + 3. * cost * rap * rap * (yeah-1.)
          !! FACTOR 1/2
          fshape = 0.5*fshape

      return 
      end
