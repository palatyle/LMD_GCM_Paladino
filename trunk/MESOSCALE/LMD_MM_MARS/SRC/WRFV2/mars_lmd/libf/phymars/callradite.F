      SUBROUTINE callradite(icount,ngrid,nlayer,
     $     aerosol,albedo,
     $     emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,
     $     dtlw,dtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,fluxtop_sw)


       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Subroutine designed to call the main canonic
c   radiative transfer subroutine "lwmain" et "swmain"
c   to compute radiative heating and cooling rate and
c   radiative fluxes to the surface.
c
c   These calculations are only valid on the part of the atmosphere
c   where Local Thermal Equilibrium (NLTE) is verified. In practice
c   The calculations are only performed for the first "nlaylte"
c   parameters (nlaylte is calculated by subroutine "nlthermeq"
c   and stored in common "yomlw.h"
c      
c
c
c   The purpose of this subroutine is
c      1) Make some initial calculation at first call
c      2) Split the calculation in several sub-grid
c        ("sub-domain") to save memory and
c        be able run on a workstation at high resolution
c        The sub-grid size is defined in dimradmars.h
c       3) call "lwmain" and "swmain"
c
c
c   author:   
c   ------
c           Francois Forget / Christophe Hourdin
c
c   This version modified to only calculate radiative tendencies over
c   layers 1..NFLEV (set in dimradmars.h).  Returns zero for higher
c   layers, if any.
c   In other routines, nlayermx -> nflev.
c   Routines affected: lwflux, lwi, lwmain, lwxb, lwxd, lwxn.
c   SRL 7/2000
c
c   definition:
c   ----------
c   Here, solar band#1 is spectral interval between "long1vis" and "long2vis"
c   set in dimradmars.h 
c   Here, solar band#2 is spectral interval between "long2vis" and "long3vis"
c   set in dimradmars.h 
c
c   input:
c   ----- 
c   icount                counter of call to subroutine physic by gcm
c   ngrid                 number of gridpoint of horizontal grid
c   nlayer                Number of layer
c
c   aerosol(ngrid,nlayer,naerkind)    aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars.h , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.  
c 
c   albedo (ngrid,2)      hemispheric surface albedo
c                         albedo (i,1) : mean albedo for solar band#1 
c                                        (see below)
c                         albedo (i,2) : mean albedo for solar band#2
c                                        (see below)
c   mu0(ngridmx)           cos of solar zenith angle (=1 when sun at zenith)
c   pplay(ngrid,nlayer)    pressure (Pa) in the middle of each layer
c   pplev(ngrid,nlayer+1)  pressure (Pa) at boundaries of each layer
c   pt(ngrid,nlayer)       atmospheric temperature in each layer (K)
c   tsurf(ngrid)           surface temperature (K)
c   fract(ngridmx)         day fraction of the time interval 
c                          =1 during the full day ; =0 during the night
c   declin                 latitude of subsolar point
c   dist_sol               sun-Mars distance (AU)
c   igout                  coordinate of analysed point for debugging
c
c  output:
c  -------
c dtlw (ngrid,nlayer)       longwave (IR) heating rate (K/s)
c dtsw(ngrid,nlayer)        shortwave (Solar) heating rate (K/s)
c fluxsurf_lw(ngrid)        surface downward flux tota LW (thermal IR) (W.m-2)
c fluxsurf_sw(ngrid,1)      surface downward flux SW for solar band#1 (W.m-2)
c fluxsurf_sw(ngrid,2)      surface downward flux SW for solar band#2 (W.m-2)
c
c fluxtop_lw(ngrid)         outgoing upward flux tota LW (thermal IR) (W.m-2)
c fluxtop_sw(ngrid,1)       outgoing upward flux SW for solar band#1 (W.m-2)
c fluxtop_sw(ngrid,2)       outgoing upward flux SW for solar band#2 (W.m-2)

c=======================================================================
c
c    Declarations :
c    -------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "yomlw.h"


c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER Icount        
      INTEGER ngrid,nlayer  
      INTEGER igout

      REAL aerosol(ngrid,nlayer,naerkind)
      REAL albedo(ngrid,2),emis(ngrid)

      REAL pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL pt(ngrid,nlayer)
      REAL tsurf(ngrid)
      REAL dist_sol,mu0(ngrid),fract(ngrid)
      REAL dtlw(ngridmx,nlayermx),dtsw(ngridmx,nlayermx)
      REAL fluxsurf_lw(ngridmx), fluxtop_lw(ngridmx)
      REAL fluxsurf_sw(ngridmx,2), fluxtop_sw(ngridmx,2)
      REAL flux(ngridmx,6)


c
c    Local variables :
c    -----------------

      INTEGER l,ig, n
      INTEGER jd,j,ig0,nd

      real  cste_mars ! solar constant on Mars (Wm-2)
      REAL ptlev(ngridmx,nlayermx+1)


      INTEGER ndomain
      parameter (ndomain = (ngridmx-1) / ndomainsz + 1)

c     Thermal IR net radiative budget (W m-2)
      real znetrad(ndomainsz,nflev) 

      real zfluxd_sw(ndomainsz,nflev+1,2)
      real zfluxu_sw(ndomainsz,nflev+1,2)


      REAL zplev(ndomainsz,nflev+1)
      REAL zztlev(ndomainsz,nflev+1)
      REAL zplay(ndomainsz,nflev)
      REAL zt(ndomainsz,nflev)
      REAL zaerosol(ndomainsz,nflev,naerkind)
      REAL zalbedo(ndomainsz,2)
      REAL zdp(ndomainsz,nflev)
      REAL zdt0(ndomainsz)

      REAL zzdtlw(ndomainsz,nflev)
      REAL zzdtsw(ndomainsz,nflev)
      REAL zzflux(ndomainsz,6)
      real zrmuz


c   local saved variables
c   ---------------------

      real pview(ngridmx)
      save pview
      
      real zco2   ! volume fraction of CO2 in Mars atmosphere
      DATA zco2/0.95/
      SAVE zco2

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

c----------------------------------------------------------------------

c     Initialisation
c     --------------

      IF (firstcall) THEN
         DO ig=1,ngrid
            pview(ig)=1.66     ! cosecant of viewing angle
         ENDDO
         gcp = g/cpp
         CALL SUAER
         CALL SULW

         write(*,*) 'Splitting radiative calculations: ',
     $              ' ngridmx,ngrid,ndomainsz,ndomain',
     $                ngridmx,ngrid,ndomainsz,ndomain
         if (ngridmx .EQ. 1) then
           if (ndomainsz .NE. 1) then
             print*
             print*,'ATTENTION !!!'
             print*,'pour tourner en 1D, '
             print*,'fixer ndomainsz=1 dans phymars/dimradmars.h'
             print*
             call exit(1)
           endif
         endif
         firstcall=.false.
      END IF

c     Starting loop on sub-domain
c     ----------------------------

      DO jd=1,ndomain
        ig0=(jd-1)*ndomainsz
        if (jd.eq.ndomain) then
         nd=ngridmx-ig0
        else
         nd=ndomainsz
        endif
        

c       Spliting input variable in sub-domain input variables
c       ---------------------------------------------------
        do l=1,nlaylte+1
         do ig = 1,nd
          zplev(ig,l) = pplev(ig0+ig,l)
         enddo
        enddo

        do l=1,nlaylte
         do ig = 1,nd
          zplay(ig,l) = pplay(ig0+ig,l)
          zt(ig,l) = pt(ig0+ig,l)

c         Thickness of each layer (Pa) :
          zdp(ig,l)= pplev(ig0+ig,l) - pplev(ig0+ig,l+1)
         enddo
        enddo

        do n=1,naerkind
          do l=1,nlaylte
            do ig=1,nd
              zaerosol(ig,l,n) = aerosol(ig0+ig,l,n)
            enddo
          enddo
        enddo

        do j=1,2
          do ig = 1,nd
           zalbedo(ig,j) = albedo(ig0+ig,j)
          enddo
        enddo

c       Intermediate  levels: (computing tlev)
c       ---------------------------------------
c       Extrapolation for the air temperature above the surface
        DO ig=1,nd
              zztlev(ig,1)=zt(ig,1)+
     s        (zplev(ig,1)-zplay(ig,1))*
     s        (zt(ig,1)-zt(ig,2))/(zplay(ig,1)-zplay(ig,2))

              zdt0(ig) = tsurf(ig0+ig) - zztlev(ig,1)
        ENDDO

        DO l=2,nlaylte
         DO ig=1,nd
               zztlev(ig,l)=0.5*(zt(ig,l-1)+zt(ig,l))
         ENDDO
        ENDDO

        DO ig=1,nd
           zztlev(ig,nlaylte+1)=zt(ig,nlaylte)
        ENDDO


c       Longwave ("lw") radiative transfer (= thermal infrared)
c       -------------------------------------------------------

        !PRINT*,'ig0',ig0
        !PRINT*,'icount',icount
        !PRINT*,'nd',nd
        !PRINT*,'nflev',nflev
        !PRINT*,'zdp',zdp
        !PRINT*,'zdt0',zdt0
        !PRINT*,'emis',emis(ig0+1)
        !PRINT*,'zplev',zplev
        !PRINT*,'zztlev',zztlev
        !PRINT*,'zt',zt
        !PRINT*,'zaerosol',zaerosol
        !PRINT*,'zzdtlw',zzdtlw
        !PRINT*,'fluxsurf_lw',fluxsurf_lw(ig0+1)
        !PRINT*,'fluxtop_lw',fluxtop_lw(ig0+1)
        !PRINT*,'znetrad',znetrad

        call lwmain (ig0,icount,nd,nflev
     .        ,zdp,zdt0,emis(ig0+1),zplev,zztlev,zt
     .        ,zaerosol,zzdtlw
     .        ,fluxsurf_lw(ig0+1),fluxtop_lw(ig0+1)
     .        ,znetrad)

c       Shortwave ("sw") radiative transfer (= solar radiation)
c       -------------------------------------------------------
c          Mars solar constant (W m-2)
c          1370 W.m-2 is the solar constant at 1 AU.
           cste_mars=1370./(dist_sol*dist_sol)

           call swmain ( nd, nflev, 
     S     cste_mars, zalbedo, 
     S     mu0(ig0+1), zdp, zplev, zaerosol, fract(ig0+1),
     S     zzdtsw, zfluxd_sw, zfluxu_sw)

c       Un-spliting output variable from sub-domain input variables
c       ------------------------------------------------------------

        do l=1,nlaylte
         do ig = 1,nd
          dtlw(ig0+ig,l) = zzdtlw(ig,l)
          dtsw(ig0+ig,l) = zzdtsw(ig,l)
         enddo
        enddo

        do l=1,nlaylte+1
         do ig = 1,nd
          ptlev(ig0+ig,l) = zztlev(ig,l)
         enddo
        enddo

        do ig = 1,nd
          fluxsurf_sw(ig0+ig,1) = zfluxd_sw(ig,1,1)
          fluxsurf_sw(ig0+ig,2) = zfluxd_sw(ig,1,2)
          fluxtop_sw(ig0+ig,1) = zfluxu_sw(ig,nlaylte+1,1)
          fluxtop_sw(ig0+ig,2) = zfluxu_sw(ig,nlaylte+1,2)
        enddo


      ENDDO         !   (boucle jd=1, ndomain)


c     Zero tendencies for any remaining layers between nlaylte and nlayer
      if (nlayer.gt.nlaylte) then
         do l = nlaylte+1, nlayer
            do ig = 1, ngrid
               dtlw(ig, l) = 0.
               dtsw(ig, l) = 0.
            enddo
         enddo
      endif
        

c     Output for debugging if lwrite=T
c     --------------------------------
c     Write all nlayer layers, even though only nlaylte layers may have
c     non-zero tendencies.

         IF(lwrite) THEN
            PRINT*,'Diagnotique for the radiation'
            PRINT*,'albedo, emissiv, mu0,fract,fluxsurf_lw,fluxsurf_sw'
            PRINT*,albedo(igout,1),emis(igout),mu0(igout),
     s           fract(igout), fluxsurf_lw(igout),
     $     fluxsurf_sw(igout,1)+fluxsurf_sw(igout,2)
            PRINT*,'Tlay Tlev Play Plev dT/dt SW dT/dt LW (K/s)'
            PRINT*,'daysec',daysec
            DO l=1,nlayer
               PRINT*,pt(igout,l),ptlev(igout,l),
     s         pplay(igout,l),pplev(igout,l),
     s         dtsw(igout,l),dtlw(igout,l)
            ENDDO
         ENDIF


      return
      end 

