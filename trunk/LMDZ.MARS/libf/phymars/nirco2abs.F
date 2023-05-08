      SUBROUTINE nirco2abs(ngrid,nlayer,pplay,dist_sol,nq,pq,
     $     mu0,fract,declin,pdtnirco2)
                                                   
       use tracer_mod, only: igcm_co2, igcm_o
       use comgeomfi_h, only: sinlon, coslon, sinlat, coslat
       USE comcstfi_h, ONLY: pi
       USE time_phylmdz_mod, ONLY: daysec
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
c   jul 2011 malv+fgg    New corrections for NLTE implemented
c   08/2002 : correction for bug when running with diurnal=F
c
c   author:  Frederic Hourdin 1996 
c   ------
c            Francois Forget 1999 
c
c   input:
c   ----- 
c   ngrid                 number of gridpoint of horizontal grid
c   nlayer                Number of layer
c   dist_sol              sun-Mars distance (AU)
c   mu0(ngrid)          
c   fract(ngrid)          day fraction of the time interval
c   declin                latitude of subslar point
c
c   output:
c   -------
c
c   pdtnirco2(ngrid,nlayer)      Heating rate (K/s)
c
c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c
#include "callkeys.h"
#include "nirdata.h"

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      integer,intent(in) :: ngrid ! number of (horizontal) grid points
      integer,intent(in) :: nlayer ! number of atmospheric layers
      real,intent(in) :: pplay(ngrid,nlayer) ! Pressure
      real,intent(in) :: dist_sol ! Sun-Mars distance (in AU)
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: pq(ngrid,nlayer,nq) ! tracers
      real,intent(in) :: mu0(ngrid) ! solar angle
      real,intent(in) :: fract(ngrid) ! day fraction of the time interval
      real,intent(in) :: declin ! latitude of sub-solar point
      
      real,intent(out) :: pdtnirco2(ngrid,nlayer) ! heating rate (K/s)
c
c    Local variables :
c    -----------------
      INTEGER l,ig, n, nstep,i
      REAL co2heat0, zmu(ngrid)

c     special diurnal=F 
      real mu0_int(ngrid),fract_int(ngrid),zday_int
      real ztim1,ztim2,ztim3,step

c
c   local saved variables
c   ---------------------
      logical,save :: firstcall=.true.
      integer,save :: ico2=0 ! index of "co2" tracer
      integer,save :: io=0 ! index of "o" tracer
c     p0noonlte is a pressure below which non LTE effects are significant.
c     REAL p0nonlte
c     DATA p0nonlte/7.5e-3/
c     SAVE p0nonlte

c     parameters for CO2 heating fit
      real n_a, n_p0, n_b
      parameter (n_a=1.1956475)
      parameter (n_b=1.9628251)
      parameter (n_p0=0.0015888279)

c     Variables added to implement NLTE correction factor (feb 2011)
      real    pyy(nlayer)
      real    cor1(nlayer),oldoco2(nlayer),alfa2(nlayer)
      real    p2011,cociente1,merge
      real    cor0,oco2gcm

c----------------------------------------------------------------------

c     Initialisation
c     --------------
      ! AS: OK firstcall absolute
      if (firstcall) then
        if (nircorr.eq.1) then
          ! we will need co2 and o tracers
          ico2=igcm_co2
          if (ico2==0) then
            write(*,*) "nirco2abs error: I need a CO2 tracer"
            write(*,*) "     when running with nircorr==1"
            stop
          endif
          io=igcm_o
          if (io==0) then
            write(*,*) "nirco2abs error: I need an O tracer"
            write(*,*) "     when running with nircorr==1"
            stop
          endif
        endif
        firstcall=.false.
      endif


c     co2heat is the heating by CO2 at 700Pa for a zero zenithal angle.
      co2heat0=n_a*(1.52/dist_sol)**2/daysec

c     Simple calcul for a given sun incident angle (if diurnal=T)
c     --------------------------------------------

      IF (diurnal) THEN  
         do ig=1,ngrid
            zmu(ig)=sqrt(1224.*mu0(ig)*mu0(ig)+1.)/35.

            if(nircorr.eq.1) then
               do l=1,nlayer
                  pyy(l)=pplay(ig,l)
               enddo

               call interpnir(cor1,pyy,nlayer,corgcm,pres1d,npres)
               call interpnir(oldoco2,pyy,nlayer,oco21d,pres1d,npres)
               call interpnir(alfa2,pyy,nlayer,alfa,pres1d,npres)
            endif

            do l=1,nlayer
!           Calculations for the O/CO2 correction
               if(nircorr.eq.1) then
                  cor0=1./(1.+n_p0/pplay(ig,l))**n_b
                  if(pq(ig,l,ico2).gt.1.e-6) then
                     oco2gcm=pq(ig,l,io)/pq(ig,l,ico2)
                  else
                     oco2gcm=1.e6
                  endif
                  cociente1=oco2gcm/oldoco2(l)
                  merge=alog10(cociente1)*alfa2(l)+alog10(cor0)*
     $                 (1.-alfa2(l))
                  merge=10**merge
                  p2011=sqrt(merge)*cor0
               else if (nircorr.eq.0) then
                  p2011=1.
                  cor1(l)=1.
               endif

               if(fract(ig).gt.0.) pdtnirco2(ig,l)=
     &             co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
     &             /(1.+n_p0/pplay(ig,l))**n_b
!           Corrections from tabulation
     $              * cor1(l) * p2011
c           OLD SCHEME (forget et al. 1999)
c    s           co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
c    s          / (1+p0nonlte/pplay(ig,l))
           enddo
         enddo
         

c     Averaging over diurnal cycle (if diurnal=F)
c     -------------------------------------------
c     NIR CO2 abs is slightly non linear. To remove the diurnal
c     cycle, it is better to average the heating rate over 1 day rather
c     than using the mean mu0 computed by mucorr in physiq.F (FF, 1998)

      ELSE      ! if (.not.diurnal) then

         nstep = 20   ! number of integration step /sol
         mu0_int(1:ngrid) = 0.
         ztim1 = 0.
         do n=1,nstep
            zday_int = (n-1)/float(nstep)
            ztim2=COS(declin)*COS(2.*pi*(zday_int-.5))
            ztim3=-COS(declin)*SIN(2.*pi*(zday_int-.5))
            CALL solang(ngrid,sinlon,coslon,sinlat,coslat,
     s             ztim1,ztim2,ztim3,
     s             mu0_int,fract_int)
            do ig=1,ngrid               
               zmu(ig)=sqrt(1224.*mu0_int(ig)*mu0_int(ig)+1.)/35.

               if(nircorr.eq.1) then
                  do l=1,nlayer
                     pyy(l)=pplay(ig,l)
                  enddo
               call interpnir(cor1,pyy,nlayer,corgcm,pres1d,npres)
               call interpnir(oldoco2,pyy,nlayer,oco21d,pres1d,npres)
               call interpnir(alfa2,pyy,nlayer,alfa,pres1d,npres)
               endif

               do l=1,nlayer
                  if(nircorr.eq.1) then
                     cor0=1./(1.+n_p0/pplay(ig,l))**n_b
                     oco2gcm=pq(ig,l,io)/pq(ig,l,ico2)
                     cociente1=oco2gcm/oldoco2(l)
                     merge=alog10(cociente1)*alfa2(l)+alog10(cor0)*
     $                    (1.-alfa2(l))
                     merge=10**merge
                     p2011=sqrt(merge)*cor0
                  else if (nircorr.eq.0) then
                     p2011=1.
                     cor1(l)=1.
                  endif

                  if(fract_int(ig).gt.0.) pdtnirco2(ig,l)=
     &                 pdtnirco2(ig,l) + (1/float(nstep))*
     &                 co2heat0*sqrt((700.*zmu(ig))/pplay(ig,l))
     &                 /(1.+n_p0/pplay(ig,l))**n_b
!                      Corrections from tabulation
     $                 * cor1(l) * p2011
               enddo
            enddo
         end do
      END IF

      return
      end 


      
      subroutine interpnir(escout,p,nlayer,escin,pin,nl)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlayer) on pressure grid p(nlayer).
C
      real escout(nlayer),p(nlayer)
      real escin(nl),pin(nl),wm,wp
      integer nl,nlayer,n1,n,nm,np
      do n1=1,nlayer
         if(p(n1) .gt. 1500. .or. p(n1) .lt. 1.0e-13) then
            escout(n1) = 0.0
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
