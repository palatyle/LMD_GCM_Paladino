*******************************************************
*                                                     *
      subroutine nucleaco2(pco2,temp,sat,n_ccn,nucrate,
     &           n_ccn_h2oice,rad_h2oice,nucrate_h2oice,
     &           vo2co2)
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
*   CO2 nucleation routine dev. by Constantino        *
*     Listowski and Joachim Audouard (2016-2017),     *
*     adapted from the water ice nucleation     
* It computes two different nucleation rates : one 
* on the dust CCN distribution and the other one on
* the water ice particles distribution
*******************************************************
 ! nucrate = output
 ! nucrate_h2o en sortie aussi : 
!nucleation sur dust et h2o separement ici

      include "microphys.h"
      include "callkeys.h"

c     Inputs
      DOUBLE PRECISION pco2,sat,vo2co2
      DOUBLE PRECISION n_ccn(nbinco2_cld), n_ccn_h2oice(nbinco2_cld)
      REAL temp !temperature
c     Output
      DOUBLE PRECISION nucrate(nbinco2_cld)
      DOUBLE PRECISION nucrate_h2oice(nbinco2_cld) ! h2o as substrate
      double precision rad_h2oice(nbinco2_cld) ! h2o ice grid (as substrate)

c     Local variables
      DOUBLE PRECISION nco2
      DOUBLE PRECISION rstar    ! Radius of the critical germ (m)
      DOUBLE PRECISION gstar    ! # of molecules forming a critical embryo
      DOUBLE PRECISION fistar   ! Activation energy required to form a critical embryo (J)
      DOUBLE PRECISION fshapeco2   ! function defined at the end of the file
      DOUBLE PRECISION deltaf      
      double precision mtetalocal,mtetalocalh ! local mteta in double precision
      double precision fshapeco2simple,zefshapeco2
      integer i
c     *************************************************

      mtetalocal = dble(mtetaco2)  !! use mtetalocal for better performance
      mtetalocalh=dble(mteta)


      IF (sat .gt. 1.) THEN    ! minimum condition to activate nucleation

        nco2   = pco2 / kbz / temp
        rstar  = 2. * sigco2 * vo2co2 / (kbz*temp*dlog(sat))
        gstar  = 4. * pi * (rstar * rstar * rstar) / (3.*vo2co2)
        
       fshapeco2simple = (2.+mtetalocal)*(1.-mtetalocal)*(1.-mtetalocal)
     &                   / 4.

c       Loop over size bins
        do i=1,nbinco2_cld
c            write(*,*) "IN NUCLEA, i, RAD_CLDCO2(i) = ",i, rad_cldco2(i),
c     &          n_ccn(i)

          if ( n_ccn(i) .lt. 1e-10 ) then
c           no dust, no need to compute nucleation!
            nucrate(i)=0.
c            goto 210
c          endif
          else
            if (rad_cldco2(i).gt.3000.*rstar) then
              zefshapeco2 = fshapeco2simple
            else
             zefshapeco2 = fshapeco2(mtetalocal,rad_cldco2(i)/rstar)
            endif

            fistar = (4./3.*pi) * sigco2 * (rstar * rstar) * 
     &             zefshapeco2
            deltaf = (2.*desorpco2-surfdifco2-fistar)/
     &             (kbz*temp)
            deltaf = min( max(deltaf, -100.d0), 100.d0)

            if (deltaf.eq.-100.) then
                nucrate(i) = 0.
            else
                nucrate(i)= dble(sqrt ( fistar /
     &               (3.*pi*kbz*temp*(gstar*gstar)) )
     &                  * kbz * temp * rstar
     &                  * rstar * 4. * pi
     &                  * ( nco2*rad_cldco2(i) )
     &                  * ( nco2*rad_cldco2(i) )
     &                  / ( zefshapeco2 * nusco2 * m0co2 )
     &                  * dexp (deltaf))

            
            endif
          endif ! if n_ccn(i) .lt. 1e-10

          if (co2useh2o) then

            if ( n_ccn_h2oice(i) .lt. 1e-10 ) then
c           no H2O ice, no need to compute nucleation!
               nucrate_h2oice(i)=0.
            else   
               if (rad_h2oice(i).gt.3000.*rstar) then
                 zefshapeco2 = (2.+mtetalocalh)*(1.-mtetalocalh)*
     &                (1.-mtetalocalh) / 4.
               else  ! same m for dust/h2o ice
                 zefshapeco2 = fshapeco2(mtetalocalh,
     &                            (rad_h2oice(i)/rstar))
               endif

               fistar = (4./3.*pi) * sigco2 * (rstar * rstar) * 
     &             zefshapeco2
              deltaf = (2.*desorpco2-surfdifco2-fistar)/
     &             (kbz*temp)
              deltaf = min( max(deltaf, -100.d0), 100.d0)

              if (deltaf.eq.-100.) then
                  nucrate_h2oice(i) = 0.
              else
                  nucrate_h2oice(i)= dble(sqrt ( fistar /
     &               (3.*pi*kbz*temp*(gstar*gstar)) )
     &                  * kbz * temp * rstar
     &                  * rstar * 4. * pi
     &                  * ( nco2*rad_h2oice(i) )
     &                  * ( nco2*rad_h2oice(i) )
     &                  / ( zefshapeco2 * nusco2 * m0co2 )
     &                  * dexp (deltaf))
              endif
            endif
          endif
        enddo

      ELSE ! parcelle d'air non saturée

        do i=1,nbinco2_cld
          nucrate(i) = 0.
          nucrate_h2oice(i) = 0.
        enddo

      ENDIF ! if (sat .gt. 1.) 

      end

*********************************************************
      double precision function fshapeco2(cost,rap)
      implicit none
*        function computing the f(m,x) factor           *
* related to energy required to form a critical embryo  *
*********************************************************

      double precision cost,rap
      double precision yeah

          !! PHI
          yeah = sqrt( 1. - 2.*cost*rap + rap*rap )
          !! FSHAPECO2 = TERM A
          fshapeco2 = (1.-cost*rap) / yeah
          fshapeco2 = fshapeco2 * fshapeco2 * fshapeco2
          fshapeco2 = 1. + fshapeco2
          !! ... + TERM B
          yeah = (rap-cost)/yeah
          fshapeco2 = fshapeco2 + 
     & rap*rap*rap*(2.-3.*yeah+yeah*yeah*yeah)
          !! ... + TERM C 
          fshapeco2 = fshapeco2 + 3. * cost * rap * rap * (yeah-1.)
          !! FACTOR 1/2
          fshapeco2 = 0.5*fshapeco2

      end
