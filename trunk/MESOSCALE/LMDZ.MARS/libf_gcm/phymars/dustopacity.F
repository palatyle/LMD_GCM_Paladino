      SUBROUTINE dustopacity(ngrid,nlayer,nq,zday,pplay,pplev,ls,pq,
     $    tauref,tau,aerosol)
                                                   
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing aerosol optical depth (dust opacity) 
c   In each layers
c
c   author: F.Forget 
c   ------
c   update F. Montmessin (water ice scheme) 
c      and S. Lebonnois (12/06/2003) compatibility dust/ice/chemistry
c   
c   input:
c   ----- 
c   ngrid             Number of gridpoint of horizontal grid
c   nlayer            Number of layer
c   nq                Number of tracer
c   ls                Solar longitude (Ls) , radian
c   pplay,pplev       pressure (Pa) in the middle and boundary of each layer
c   pq                Dust mixing ratio (used if tracer =T and active=T).
c
c   output:
c   -------
c   tauref       Prescribed mean column optical depth at 700 Pa 
c   tau          Column total visible dust optical depth at each point
c   aerosol      aerosol(ig,l,1) is the dust optical
c                depth in layer l, grid point ig

c
c=======================================================================
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "comgeomfi.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "tracer.h"
#include "planete.h"

c-----------------------------------------------------------------------
c
c    Declarations :
c    --------------
c
c    Input/Output
c    ------------
      INTEGER ngrid,nlayer,nq

      REAL ls,zday,expfactor    
      REAL pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL pq(ngrid,nlayer,nq)
      REAL tauref(ngrid), tau(ngrid,naerkind)
      REAL aerosol(ngrid,nlayer,naerkind)
c
c    Local variables :
c    -----------------
      INTEGER l,ig,iq
      real topdust(ngridmx)
      real zlsconst, zp
      real taueq,tauS,tauN
      real r0,reff,coefsize
c
c   local saved variables
c   ---------------------

      REAL topdust0(ngridmx) 
      SAVE topdust0

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall


c----------------------------------------------------------------------

c     Initialisation
c     --------------

      IF (firstcall) THEN

c        altitude of the top of the aerosol layer (km) at Ls=2.76rad:
c        in the Viking year scenario
         DO ig=1,ngrid
             topdust0(ig)=60. -22.*SIN(lati(ig))**2
         END DO
         firstcall=.false.
      END IF

c  -------------------------------------------------------------
c     1) Prescribed dust  (if tracer=F or active=F)
c  -------------------------------------------------------------
      IF ((.not.tracer) .or. (.not.active)) THEN 

c       Vertical column optical depth at 700.Pa 
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF(iaervar.eq.1) THEN
           do ig=1, ngridmx
            tauref(ig)=max(tauvis,1.e-9) ! tauvis=cste as read in starfi
          end do
        ELSE IF (iaervar.eq.2) THEN   ! << "Viking" Scenario>>

          tauref(1) = 0.7+.3*cos(ls+80.*pi/180.) ! like seen by VL1
          do ig=2,ngrid
            tauref(ig) = tauref(1)
          end do

        ELSE IF (iaervar.eq.3) THEN  ! << "MGS" scenario >>

           taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls-4.363)))**14
           tauS= 0.1 +(0.5-0.1)  *(cos(0.5*(ls-4.363)))**14
           tauN = 0.1
c	   if (peri_day.eq.150) then
c	     tauS=0.1
c	     tauN=0.1 +(0.5-0.1)  *(cos(0.5*(ls+pi-4.363)))**14
c	     taueq= 0.2 +(0.5-0.2) *(cos(0.5*(ls+pi-4.363)))**14
c           endif
	     
           do ig=1,ngrid/2  ! Northern hemisphere
             tauref(ig)= tauN +
     &       (taueq-tauN)*0.5*(1+tanh((45-lati(ig)*180./pi)*6/60))
           end do
           do ig=ngrid/2+1, ngridmx  ! Southern hemisphere
             tauref(ig)= tauS +
     &       (taueq-tauS)*0.5*(1+tanh((45+lati(ig)*180./pi)*6/60))
           end do

        ELSE IF (iaervar.eq.4) THEN  ! << "TES scenario >>
           call readtesassim(ngrid,nlayer,zday,pplev,tauref)

        ELSE IF (iaervar.eq.5) THEN   ! << Escalier Scenario>>
c         tauref(1) = 0.2
c         if ((ls.ge.210.*pi/180.).and.(ls.le.330.*pi/180.))
c    &                              tauref(1) = 2.5
          tauref(1) = 2.5
          if ((ls.ge.30.*pi/180.).and.(ls.le.150.*pi/180.))
     &                              tauref(1) = .2

          do ig=2,ngrid
            tauref(ig) = tauref(1)
          end do
        ELSE
          stop 'problem with iaervar in dustopacity.F'
        ENDIF

c       Altitude of the top of the dust layer
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        zlsconst=SIN(ls-2.76)
        if (iddist.eq.1) then
          do ig=1,ngrid
             topdust(ig)=topdustref         ! constant dust layer top
          end do

        else if (iddist.eq.2) then          ! "Viking" scenario
          do ig=1,ngrid
            topdust(ig)=topdust0(ig)+18.*zlsconst
          end do

        else if(iddist.eq.3) then         !"MGS" scenario
          do ig=1,ngrid
            topdust(ig)=60.+18.*zlsconst
     &                -(32+18*zlsconst)*sin(lati(ig))**4
     &                 - 8*zlsconst*(sin(lati(ig)))**5
          end do
        endif


c       Optical depth in each layer :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(iddist.ge.1) then
         expfactor=0.
          DO l=1,nlayer
             DO ig=1,ngrid
             if(pplay(ig,l).gt.700.
     $                        /(988.**(topdust(ig)/70.))) then
                zp=(700./pplay(ig,l))**(70./topdust(ig))
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)
              else    
               expfactor=1.e-3
              endif       
                aerosol(ig,l,1)= tauref(ig)/700. *
     s           (pplev(ig,l)-pplev(ig,l+1))
     &           *expfactor
c     s           *max( exp(.007*(1.-max(zp,1.))) , 1.E-3 )
             ENDDO
          ENDDO
c     changement dans le calcul de la distribution verticale          
c     dans le cas des scenarios de poussieres assimiles
c        if (iaervar.eq.4) THEN  ! TES
c        call zerophys(ngrid*naerkind,tau)
c
c        do l=1,nlayer
c           do ig=1,ngrid
c                tau(ig,1)=tau(ig,1)+ aerosol(ig,l,1)
c           end do
c        end do
c        do l=1,nlayer
c           do ig=1,ngrid
c               aerosol(ig,l,1)=aerosol(ig,l,1)*tauref(ig)/tau(ig,1)
c     $     *(pplev(ig,1)/700)
c           end do
c        end do
c        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
        else if(iddist.eq.0) then   
c         old dust vertical distribution function (pollack90)
          DO l=1,nlayer
             DO ig=1,ngrid
                zp=700./pplay(ig,l)
                aerosol(ig,l,1)= tauref(ig)/700. *
     s           (pplev(ig,l)-pplev(ig,l+1))
     s           *max( exp(.03*(1.-max(zp,1.))) , 1.E-3 )
             ENDDO
          ENDDO
        end if

c  ---------------------------------------------------------------------
c     2) Transported radiatively active dust  (if tracer=T and active=T)
c ----------------------------------------------------------------------
      ELSE  IF ((tracer) .and. (active)) THEN 
c     The dust opacity is computed from q

c     a) "doubleq" technique (transport of mass and number mixing ratio)
c        ~~~~~~~~~~~~~~~~~~~
       if(doubleq) then 

         call zerophys(ngrid*nlayer*naerkind,aerosol)

c        Computing effective radius :
         do  l=1,nlayer
           do ig=1, ngrid
              r0=
     &        (r3n_q*pq(ig,l,1)/max(pq(ig,l,2),0.01))**(1./3.)
              r0=min(max(r0,1.e-10),500.e-6)
              reff= ref_r0 * r0
cc           If  reff is small, the transported dust mean Qext
c            is reduced from the reference dust Qext by a factor "coefsize"

             coefsize=min(max(2.52e6*reff-0.043 ,0.)   ,1.)

cc              It is added 1.e-8 to pq to avoid low

                aerosol(ig,l,1)=aerosol(ig,l,1)+ 1.E-8 +
     &         ( 0.75*Qext(1)*coefsize/(rho_dust*reff))
     &          * (pq(ig,l,1))*
     &          (pplev(ig,l)-pplev(ig,l+1))/g
           end do
         end do
         call zerophys(ngrid,tauref)

c     b) Size bin technique (each aerosol can contribute to opacity))
c        ~~~~~~~~~~~~~~~~~~
       else
c        The dust opacity is computed from q
         call zerophys(ngrid*nlayer*naerkind,aerosol)
         do iq=1,dustbin
           do l=1,nlayer
              do ig=1,ngrid
cc               qextrhor(iq) is  (3/4)*Qext/(rho*reff)
cc               It is added 1.e-8 to pq to avoid low
                 aerosol(ig,l,1)=aerosol(ig,l,1)+
     &           qextrhor(iq)* (pq(ig,l,iq) + 1.e-8)*
     &           (pplev(ig,l)-pplev(ig,l+1))/g
              end do
           end do
         end do
         call zerophys(ngrid,tauref)
       end if  ! (doubleq)
      END IF   ! (dust scenario)

c --------------------------------------------------------------------------
c Column integrated visible optical depth in each point (used for diagnostic)
c --------------------------------------------------------------------------
      call zerophys(ngrid*naerkind,tau)
      do l=1,nlayer
         do ig=1,ngrid
               tau(ig,1)=tau(ig,1)+ aerosol(ig,l,1)
         end do
      end do

      return
      end 

