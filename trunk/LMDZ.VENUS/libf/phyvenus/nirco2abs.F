      SUBROUTINE nirco2abs(nlon,nlev,nplay,dist_sol,nq,pq,
     $     mu0,fract,pdtnirco2)

       use dimphy
       use geometry_mod, only: longitude_deg, latitude_deg
       use chemparam_mod, only: i_co2, i_o
c       use compo_hedin83_mod2 


       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing heating rate due to
c   absorption by CO2 in the near-infrared
c   This version includes NLTE effects
c
c   (Scheme to be described in Forget et al., JGR, 2003)
c   (old Scheme described in Forget et al., JGR, 1999)
c
c   This version updated with a new functional fit,
c   see NLTE correction-factor of Lopez-Valverde et al (1998)
c   Stephen Lewis 2000
c
c   apr 2019 d.quirino   Improving NLTE params, SOIR/SPICAV Temp comparison 
c   oct 2014 g.gilli     Coupling with photochemical model
C   jan 2014 g.gilli     Revision (following martian non-lte param)    
C   jun 2013 l.salmi     First adaptation to Venus and NIR NLTE param
c   jul 2011 malv+fgg    New corrections for NLTE implemented
c   08/2002 : correction for bug when running with diurnal=F
c
c   author:  Frederic Hourdin 1996 
c   ------
c            Francois Forget 1999 
c
c   input:
c   ----- 
c   nlon                 number of gridpoint of horizontal grid
c   nlev                Number of layer
c   dist_sol              sun-Venus distance (AU)
c   mu0(nlon)          
c   fract(nlon)        day fraction of the time interval
c   declin                latitude of subslar point
c
c   output:
c   -------
c
c   pdtnirco2(nlon,nlev)      Heating rate (K/sec)
c
c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c

#include "YOMCST.h"
#include "clesphys.h"
c#include "comdiurn.h"
#include "nirdata.h"
c#include "tracer.h"
#include "mmol.h"
c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      integer,intent(in) :: nlon ! number of (horizontal) grid points
      integer,intent(in) :: nlev ! number of atmospheric layers

      real,intent(in) :: nplay(nlon,nlev) ! Pressure
      real,intent(in) :: dist_sol ! Sun-Venus distance (in AU)
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: pq(nlon,nlev,nq) ! mass mixing ratio tracers
      real,intent(in) :: mu0(nlon) ! solar angle
      real,intent(in) :: fract(nlon) ! day fraction of the time interval
c      real,intent(in) :: declin ! latitude of sub-solar point
      real :: co2vmr_gcm(nlon,nlev), o3pvmr_gcm(nlon,nlev)
 
      real,intent(out) :: pdtnirco2(nlon,nlev) ! heating rate (K/sec) 

c
c    Local variables :
c    -----------------
      INTEGER l,ig, n, nstep,i
      REAL co2heat0, zmu(nlon)

c     special diurnal=F 
      real mu0_int(nlon),fract_int(nlon),zday_int
      real ztim1,ztim2,ztim3,step

c
c   local saved variables
c   ---------------------
      logical,save :: firstcall=.true.
      integer,save :: ico2=0 ! index of "co2" tracer
      integer,save :: io=0 ! index of "o" tracer

cccc     parameters for CO2 heating fit
c
c     n_a  =  heating rate for Venusian day at p0, r0, mu =0 [K day-1]
c     Here p0 = p_cloud top [Pa]
c     n_p0 = is a pressure below which non LTE effects are significant [Pa]
c     n_a Solar heating [K/Eday] at the cloud top, taken from Crisps table     
 
      real n_a, n_p0, n_b, p_ctop

    
cc "Nominal" values used in Gilli+2'17
c       parameter (n_a = 18.13/86400.0)     !c     K/Eday  ---> K/sec   
c       parameter (p_ctop=13.2e2)
c       parameter (n_p0=0.008) 

cc "New" values used to improve SPICAV/SOIR Temperature comparision (D.Quirino)
       parameter (n_a = 15.92/86400.0)     !c     K/Eday  ---> K/sec   
       parameter (p_ctop=19.85e2)
       parameter (n_p0=0.1)  
       parameter (n_b=1.362)

c    -- NLTE Param v2  --
C       parameter (n_p0=0.01)  
c       parameter (n_b = 1.3)
   

c     Variables added to implement NLTE correction factor (feb 2011)
      real    pyy(nlev)
      real    cor1(nlev),oldoco2(nlev),alfa2(nlev)
      real    p2011,cociente1,merge
      real    cor0,oco2gcm
!!!!
c      real :: pic27(nlon,nlev), pic27b(nlon,nlev)
c      real :: pic43(nlon,nlev), picnir(nlon,nlev)

c     co2heat is the heating by CO2 at p_ctop=13.2e2 for a zero zenithal angle.

      co2heat0=n_a*(0.72/dist_sol)**2     

CCCCCC   TEST: reduce by X% nir Heating
c      co2heat0  = co2heat0 * 0.8

c----------------------------------------------------------------------
      
c     Initialisation
c     --------------
      if (firstcall) then
        if (nircorr.eq.1) then
c          ! we will need co2 and o tracers
          ico2= i_co2
          if (ico2==0) then
            write(*,*) "nirco2abs error: I need a CO2 tracer"
            write(*,*) "     when running with nircorr==1"
           stop
          endif
          io=i_o
          if (io==0) then
            write(*,*) "nirco2abs error: I need an O tracer"
            write(*,*) "     when running with nircorr==1"
            stop
          endif
        endif
        firstcall=.false.
      endif

      
c      
c     Simple calcul for a given sun incident angle (if cycle_diurne=T)
c     --------------------------------------------

      IF (cycle_diurne) THEN  

         do ig=1,nlon    
            zmu(ig)=sqrt(1224.*mu0(ig)*mu0(ig)+1.)/35.

           
            if(nircorr.eq.1) then
               do l=1,nlev
                  pyy(l)=nplay(ig,l)
               enddo

               call interpnir(cor1,pyy,nlev,corgcm,pres1d,npres)
               call interpnir(oldoco2,pyy,nlev,oco21d,pres1d,npres)
               call interpnir(alfa2,pyy,nlev,alfa,pres1d,npres)
               
            endif

            do l=1,nlev
     
c           Calculations for the O/CO2 correction
               if(nircorr.eq.1) then
                  cor0=1./(1.+n_p0/nplay(ig,l))**n_b
                  if(pq(ig,l,ico2) .gt. 1.e-6) then
                     oco2gcm=pq(ig,l,io)/pq(ig,l,ico2)

                  else
                     oco2gcm=1.e6
                  endif
                  cociente1=oco2gcm/oldoco2(l)
                  
c                  WRITE(*,*) "nirco2abs line 211", l, cociente1

                  merge=alog10(cociente1)*alfa2(l)+alog10(cor0)*
     $                 (1.-alfa2(l))
                  merge=10**merge
                  p2011=sqrt(merge)*cor0

               else if (nircorr.eq.0) then
                  p2011=1.
                  cor1(l)=1.
               endif

              if(fract(ig).gt.0.) pdtnirco2(ig,l)=
     &             co2heat0*sqrt((p_ctop*zmu(ig))/nplay(ig,l))
     &             /(1.+n_p0/nplay(ig,l))**n_b
c           Corrections from tabulation
     $              * cor1(l) * p2011
              
          enddo
         enddo
	 
c     Averaging over diurnal cycle (if diurnal=F)
c     -------------------------------------------
c     NIR CO2 abs is slightly non linear. To remove the diurnal
c     cycle, it is better to average the heating rate over 1 day rather
c     than using the mean mu0 computed by mucorr in physiq.F (FF, 1998)

      ELSE      ! if (.not.diurnal) then
         nstep = 20    ! number of integration step /sol
         do n=1,nstep

            zday_int = (n-1)/float(nstep)

            CALL zenang(0.,zday_int,RDAY/nstep,
     &                  latitude_deg,longitude_deg,
     &                  mu0_int,fract_int)

            do ig=1,nlon
               zmu(ig)=sqrt(1224.*mu0_int(ig)*mu0_int(ig)+1.)/35.

               if(nircorr.eq.1) then
                  do l=1,nlev
                     pyy(l)=nplay(ig,l)
                  enddo
                 call interpnir(cor1,pyy,nlev,corgcm,pres1d,npres)
                 call interpnir(oldoco2,pyy,nlev,oco21d,pres1d,npres)
                 call interpnir(alfa2,pyy,nlev,alfa,pres1d,npres)
               endif
c

               do l=1,nlev
c           Calculations for the O/CO2 correction
               if(nircorr.eq.1) then
                  cor0=1./(1.+n_p0/nplay(ig,l))**n_b
                  oco2gcm=pq(ig,l,io)/pq(ig,l,ico2)
                  cociente1=oco2gcm/oldoco2(l)
                  merge=alog10(cociente1)*alfa2(l)+alog10(cor0)*
     $                 (1.-alfa2(l))
                  merge=10**merge
                  p2011=sqrt(merge)*cor0

               else if (nircorr.eq.0) then
                  p2011=1.
                  cor1(l)=1.
               endif

               if(fract_int(ig).gt.0.) pdtnirco2(ig,l)=
     &              pdtnirco2(ig,l) + (1/float(nstep))*
     &              co2heat0*sqrt((p_ctop*zmu(ig))/nplay(ig,l))
     &              /(1.+n_p0/nplay(ig,l))**n_b
!     Corrections from tabulation
     $              * cor1(l) * p2011

               enddo
            enddo
         end do
      

      END IF  

      return
      end 

      
      subroutine interpnir(escout,p,nlev,escin,pin,nl)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlev) on pressure grid p(nlev).
C
      real escout(nlev),p(nlev)
      real escin(nl),pin(nl),wm,wp
      integer nl,nlev,n1,n,nm,np
      do n1=1,nlev
         if(p(n1) .gt. 1500. .or. p(n1) .lt. 1.0e-13) then
c            escout(n1) = 0.0
            escout(n1) = 1.e-15
         else
            do n = 1,nl-1
               if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
                  wp=1.0 - wm
               endif
            enddo
            escout(n1) = escin(nm)*wm + escin(np)*wp
         endif
      enddo
      return
      end
