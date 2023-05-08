      SUBROUTINE aeropacity(ngrid,nlayer,nq,zday,pplay,pplev,ls,pq,
     &    tauref,tau,aerosol,reffrad,
     &    QREFvis3d,QREFir3d,omegaREFvis3d,omegaREFir3d)
                                                   
! to use  'getin'
      USE ioipsl_getincom 
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing aerosol optical depth in each gridbox.
c
c   author: F.Forget 
c   ------
c   update F. Montmessin (water ice scheme) 
c      and S. Lebonnois (12/06/2003) compatibility dust/ice/chemistry
c   update J.-B. Madeleine 2008-2009:
c       - added 3D scattering by aerosols;
c       - dustopacity transferred from physiq.F to callradite.F,
c           and renamed into aeropacity.F;
c   
c   input:
c   ----- 
c   ngrid             Number of gridpoint of horizontal grid
c   nlayer            Number of layer
c   nq                Number of tracer
c   zday                  Date (time since Ls=0, in martian days)
c   ls                Solar longitude (Ls) , radian
c   pplay,pplev       pressure (Pa) in the middle and boundary of each layer
c   pq                Dust mixing ratio (used if tracer =T and active=T).
c   reffrad(ngrid,nlayer,naerkind)  Aerosol effective radius
c   QREFvis3d(ngridmx,nlayermx,naerkind) \ 3d extinction coefficients
c   QREFir3d(ngridmx,nlayermx,naerkind)  / at reference wavelengths;
c   omegaREFvis3d(ngridmx,nlayermx,naerkind) \ 3d single scat. albedo
c   omegaREFir3d(ngridmx,nlayermx,naerkind)  / at reference wavelengths;
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
      REAL reffrad(ngrid,nlayer,naerkind)
      REAL QREFvis3d(ngridmx,nlayermx,naerkind)
      REAL QREFir3d(ngridmx,nlayermx,naerkind)
      REAL omegaREFvis3d(ngridmx,nlayermx,naerkind)
      REAL omegaREFir3d(ngridmx,nlayermx,naerkind)
c
c    Local variables :
c    -----------------
      INTEGER l,ig,iq,i,j
      INTEGER iaer           ! Aerosol index
      real topdust(ngridmx)
      real zlsconst, zp
      real taueq,tauS,tauN
      real r0,reff,coefsize
c     Mean Qext(vis)/Qext(ir) profile
      real msolsir(nlayermx,naerkind)
c     Mean Qext(ir)/Qabs(ir) profile
      real mqextsqabs(nlayermx,naerkind)
c     Variables used when multiple particle sizes are used
c       for dust or water ice particles in the radiative transfer
c       (see callradite.F for more information).
      REAL taudusttmp(ngridmx)! Temporary dust opacity
                               !   used before scaling
      REAL taudustvis(ngridmx) ! Dust opacity after scaling
      REAL taudusttes(ngridmx) ! Dust opacity at IR ref. wav. as
                               !   "seen" by the GCM.
      REAL taucloudvis(ngridmx)! Cloud opacity at visible
                               !   reference wavelength
      REAL taucloudtes(ngridmx)! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)
c
c   local saved variables
c   ---------------------

      REAL topdust0(ngridmx) 
      SAVE topdust0

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

! indexes of water ice and dust tracers:
      INTEGER,SAVE :: nqdust(nqmx) ! to store the indexes of dust tracers
      INTEGER,SAVE :: i_ice=0  ! water ice
      CHARACTER(LEN=20) :: tracername ! to temporarly store text

      call zerophys(ngrid*naerkind,tau)

! identify tracers

      IF (firstcall) THEN
        ! identify tracers which are dust
        i=0
        DO iq=1,nq
          tracername=noms(iq)
          IF (tracername(1:4).eq."dust") THEN
          i=i+1
          nqdust(i)=iq
          ENDIF
        ENDDO

        IF (water.AND.activice) THEN
          i_ice=igcm_h2o_ice
          write(*,*) "aeropacity: i_ice=",i_ice
        ENDIF

c       altitude of the top of the aerosol layer (km) at Ls=2.76rad:
c       in the Viking year scenario
        DO ig=1,ngrid
            topdust0(ig)=60. -22.*SIN(lati(ig))**2
        END DO

c       typical profile of solsir and (1-w)^(-1):
        call zerophys(nlayer*naerkind,msolsir)
        call zerophys(nlayer*naerkind,mqextsqabs)
        WRITE(*,*) "Typical profiles of solsir and Qext/Qabs(IR):"
        DO iaer = 1, naerkind ! Loop on aerosol kind
          WRITE(*,*) "Aerosol # ",iaer
          DO l=1,nlayer
            DO ig=1,ngridmx
              msolsir(l,iaer)=msolsir(l,iaer)+
     &              QREFvis3d(ig,l,iaer)/
     &              QREFir3d(ig,l,iaer)
              mqextsqabs(l,iaer)=mqextsqabs(l,iaer)+
     &              (1.E0-omegaREFir3d(ig,l,iaer))**(-1)
            ENDDO
            msolsir(l,iaer)=msolsir(l,iaer)/REAL(ngridmx)
            mqextsqabs(l,iaer)=mqextsqabs(l,iaer)/REAL(ngridmx)
          ENDDO
          WRITE(*,*) "solsir: ",msolsir(:,iaer)
          WRITE(*,*) "Qext/Qabs(IR): ",mqextsqabs(:,iaer)
        ENDDO

!       load value of tauvis from callphys.def (if given there,
!       otherwise default value read from starfi.nc file will be used)
        call getin("tauvis",tauvis)

        firstcall=.false.

      END IF

      DO iaer = 1, naerkind ! Loop on aerosol kind
c     --------------------------------------------
        aerkind: SELECT CASE (iaer)
c==================================================================
        CASE(1) aerkind                             ! Dust (iaer=1)
c==================================================================

c  -------------------------------------------------------------
c     1) Prescribed dust  (if tracer=F or active=F)
c  -------------------------------------------------------------
      IF ((.not.tracer) .or. (.not.active)) THEN 

c       Vertical column optical depth at 700.Pa 
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF(iaervar.eq.1) THEN 
           do ig=1, ngridmx
            tauref(ig)=max(tauvis,1.e-9) ! tauvis=cste (set in callphys.def
                                         ! or read in starfi
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

        ELSE IF ((iaervar.eq.4).or.
     &           ((iaervar.ge.24).and.(iaervar.le.26)))
     &       THEN  ! << "TES assimilated dust scenarios >>
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

        ELSE IF (iaervar.gt.99) THEN  ! << input netcdf file >>
c*************WRF
c
c 2. customized dust opacity field
c   ex: from assimilation
       call meso_readtesassim(ngrid,nlayer,zday,pplev,tauref,
     . iaervar)
c
c*************WRF

        ELSE
          stop 'problem with iaervar in aeropacity.F'
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
          CALL zerophys(ngrid,taudusttmp)
          DO l=1,nlayer
            DO ig=1,ngrid
c             Typical mixing ratio profile 
              if(pplay(ig,l).gt.700.
     $                        /(988.**(topdust(ig)/70.))) then
                zp=(700./pplay(ig,l))**(70./topdust(ig))
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)
              else    
                expfactor=1.e-3
              endif
c             Vertical scaling function
              aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1)) *
     &          expfactor *
     &          QREFvis3d(ig,l,iaer) / QREFvis3d(ig,1,iaer)
c             Scaling factor
              taudusttmp(ig)=taudusttmp(ig)+aerosol(ig,l,iaer)
            ENDDO
          ENDDO

c         Rescaling each layer to reproduce the choosen (or
c           assimilated) dust extinction opacity at visible
c           reference wavelength, which is originally scaled
c           to an equivalent 700Pa pressure surface.
          DO l=1,nlayer
            DO ig=1,ngrid
              aerosol(ig,l,iaer) = tauref(ig) *
     &                     pplev(ig,1) / 700.E0 *
     &                     aerosol(ig,l,iaer) / taudusttmp(ig)
            ENDDO
          ENDDO

          CALL zerophys(ngrid,taudustvis)
          CALL zerophys(ngrid,taudusttes)
          DO l=1,nlayer
            DO ig=1,ngrid
              taudustvis(ig) = taudustvis(ig) + aerosol(ig,l,iaer)
              taudusttes(ig) = taudusttes(ig) + aerosol(ig,l,iaer)*
     &          QREFir3d(ig,l,iaer) / QREFvis3d(ig,l,iaer)*
     &          ( 1.E0 - omegaREFir3d(ig,l,iaer) )
            ENDDO
          ENDDO

c         Outputs
          IF (ngrid.NE.1) THEN
            CALL WRITEDIAGFI(ngridmx,'taudustTES','dust abs IR',
     &        ' ',2,taudusttes)
            CALL wstats(ngridmx,'taudustTES','dust abs IR',
     &        ' ',2,taudusttes)
          ELSE
            CALL writeg1d(ngrid,1,taudusttes,'taudusttes','NU')
          ENDIF

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
     &          * (pq(ig,l,nqdust(1)))*
c               only one dust bin to use with doubleq
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
     &           qextrhor(nqdust(iq))*(pq(ig,l,nqdust(iq))+1.e-8)*
     &           (pplev(ig,l)-pplev(ig,l+1))/g
              end do
           end do
         end do
         call zerophys(ngrid,tauref)
       end if  ! (doubleq)
      END IF   ! (dust scenario)


c==================================================================
        CASE(2) aerkind               ! Water ice crystals (iaer=2)
c==================================================================

      IF (water.AND.activice) THEN
c       1. Initialization
        CALL zerophys(ngrid*nlayer,aerosol(1,1,iaer))
        CALL zerophys(ngrid,taucloudvis)
        CALL zerophys(ngrid,taucloudtes)
c       2. Opacity calculation
        DO ig=1, ngrid
          DO l=1,nlayer
            aerosol(ig,l,iaer) =
     &        (  0.75 * QREFvis3d(ig,l,iaer) /
     &        ( rho_ice * reffrad(ig,l,iaer) )  ) *
     &        ( pq(ig,l,i_ice) + 1.E-8 ) *
     &        ( pplev(ig,l) - pplev(ig,l+1) ) / g
            taucloudvis(ig) = taucloudvis(ig) + aerosol(ig,l,iaer)
            taucloudtes(ig) = taucloudtes(ig) + aerosol(ig,l,iaer)*
     &        QREFir3d(ig,l,iaer) / QREFvis3d(ig,l,iaer) *
     &        ( 1.E0 - omegaREFir3d(ig,l,iaer) )
          ENDDO
        ENDDO
c       3. Outputs
        IF (ngrid.NE.1) THEN
          CALL WRITEDIAGFI(ngridmx,'tauTES','tauabs IR refwvl',
     &      ' ',2,taucloudtes)
          CALL wstats(ngridmx,'tauTES','tauabs IR refwvl',
     &      ' ',2,taucloudtes)
        ELSE
          CALL writeg1d(ngrid,1,taucloudtes,'tautes','NU')
        ENDIF
      ENDIF

c==================================================================
        END SELECT aerkind

c -----------------------------------------------------------------
c Column integrated visible optical depth in each point
c -----------------------------------------------------------------

      do l=1,nlayer
         do ig=1,ngrid
               tau(ig,iaer) = tau(ig,iaer) + aerosol(ig,l,iaer)
         end do
      end do

c     -----------------------------------
      ENDDO ! iaer (loop on aerosol kind)

      return
      end 
