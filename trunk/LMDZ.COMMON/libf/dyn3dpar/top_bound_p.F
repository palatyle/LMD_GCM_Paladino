!
! $Id: top_bound_p.F 1793 2013-07-18 07:13:18Z emillour $
!
      SUBROUTINE top_bound_p(vcov,ucov,teta,masse,dt,ducov)
      USE parallel_lmdz
      USE comvert_mod, ONLY: presnivs,preff,scaleheight
      USE comconst_mod, ONLY: iflag_top_bound,tau_top_bound,
     .			mode_top_bound
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"


c ..  DISSIPATION LINEAIRE A HAUT NIVEAU, RUN MESO,
C     F. LOTT DEC. 2006
c                                 (  10/12/06  )

c=======================================================================
c
c   Auteur:  F. LOTT  
c   -------
c
c   Objet:
c   ------
c
c   Dissipation linéaire (ex top_bound de la physique)
c
c=======================================================================

! top_bound sponge layer model:
! Quenching is modeled as: A(t)=Am+A0*exp(-lambda*dt)
! where Am is the zonal average of the field (or zero), and lambda the inverse
! of the characteristic quenching/relaxation time scale
! Thus, assuming Am to be time-independent, field at time t+dt is given by:
! A(t+dt)=A(t)-(A(t)-Am)*(1-exp(-lambda*t))
! Moreover lambda can be a function of model level (see below), and relaxation
! can be toward the average zonal field or just zero (see below).

! NB: top_bound sponge is only called from leapfrog if ok_strato=.true.

! sponge parameters: (loaded/set in conf_gcm.F ; stored in comconst_mod)
!    iflag_top_bound=0 for no sponge
!    iflag_top_bound=1 for sponge over 4 topmost layers
!    iflag_top_bound=2 for sponge from top to ~1% of top layer pressure
!    mode_top_bound=0: no relaxation
!    mode_top_bound=1: u and v relax towards 0
!    mode_top_bound=2: u and v relax towards their zonal mean
!    mode_top_bound=3: u,v and pot. temp. relax towards their zonal mean
!    tau_top_bound : inverse of charactericstic relaxation time scale at
!                       the topmost layer (Hz)


#include "comdissipn.h"
#include "iniprint.h"

c   Arguments:
c   ----------

      real,intent(inout) :: ucov(iip1,jjp1,llm) ! covariant zonal wind
      real,intent(inout) :: vcov(iip1,jjm,llm) ! covariant meridional wind
      real,intent(inout) :: teta(iip1,jjp1,llm) ! potential temperature
      real,intent(in) :: masse(iip1,jjp1,llm) ! mass of atmosphere 
      real,intent(in) :: dt ! time step (s) of sponge model
      real,intent(out) :: ducov(iip1,jjp1,llm) ! increment on ucov due to sponge

c   Local:
c   ------
      REAL massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm),zm
      REAL uzon(jjp1,llm),vzon(jjm,llm),tzon(jjp1,llm)
      
      integer i	
      REAL,SAVE :: rdamp(llm) ! quenching coefficient
      real,save :: lambda(llm) ! inverse or quenching time scale (Hz)
      LOGICAL,SAVE :: first=.true.
      INTEGER j,l,jjb,jje


      if (first) then
c$OMP BARRIER
c$OMP MASTER
         if (iflag_top_bound == 1) then
! sponge quenching over the topmost 4 atmospheric layers
             lambda(:)=0.
             lambda(llm)=tau_top_bound
             lambda(llm-1)=tau_top_bound/2.
             lambda(llm-2)=tau_top_bound/4.
             lambda(llm-3)=tau_top_bound/8.
         else if (iflag_top_bound == 2) then
! sponge quenching over topmost layers down to pressures which are
! higher than 100 times the topmost layer pressure
             lambda(:)=tau_top_bound
     s       *max(presnivs(llm)/presnivs(:)-0.01,0.)
         endif

! quenching coefficient rdamp(:)
!         rdamp(:)=dt*lambda(:) ! Explicit Euler approx.
         rdamp(:)=1.-exp(-lambda(:)*dt)

         write(lunout,*)'TOP_BOUND mode',mode_top_bound
         write(lunout,*)'Sponge layer coefficients'
         write(lunout,*)'p (Pa)  z(km)  tau(s)   1./tau (Hz)'
         do l=1,llm
           if (rdamp(l).ne.0.) then
             write(lunout,'(6(1pe12.4,1x))')
     &        presnivs(l),log(preff/presnivs(l))*scaleheight,
     &           1./lambda(l),lambda(l)
           endif
         enddo
         first=.false.
c$OMP END MASTER
c$OMP BARRIER
      endif ! of if (first)


      CALL massbar_p(masse,massebx,masseby)

      ! compute zonal average of vcov (or set it to zero)
      if (mode_top_bound.ge.2) then
       jjb=jj_begin
       jje=jj_end
       IF (pole_sud) jje=jj_end-1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
        do j=jjb,jje
          zm=0.
          vzon(j,l)=0
          do i=1,iim
! NB: we can work using vcov zonal mean rather than v since the
! cv coefficient (which relates the two) only varies with latitudes 
            vzon(j,l)=vzon(j,l)+vcov(i,j,l)*masseby(i,j,l)
            zm=zm+masseby(i,j,l)
          enddo
          vzon(j,l)=vzon(j,l)/zm
        enddo
       enddo
c$OMP END DO NOWAIT   
      else
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
         vzon(:,l)=0.
       enddo
c$OMP END DO NOWAIT
      endif ! of if (mode_top_bound.ge.2)

      ! compute zonal average of u (or set it to zero)
      if (mode_top_bound.ge.2) then
       jjb=jj_begin
       jje=jj_end
       IF (pole_nord) jjb=jj_begin+1
       IF (pole_sud)  jje=jj_end-1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
        do j=jjb,jje
          uzon(j,l)=0.
          zm=0.
          do i=1,iim
            uzon(j,l)=uzon(j,l)+massebx(i,j,l)*ucov(i,j,l)/cu(i,j)
            zm=zm+massebx(i,j,l)
          enddo
          uzon(j,l)=uzon(j,l)/zm
        enddo
       enddo
c$OMP END DO NOWAIT
      else
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
         uzon(:,l)=0.
       enddo
c$OMP END DO NOWAIT
      endif ! of if (mode_top_bound.ge.2)

      ! compute zonal average of potential temperature, if necessary
      if (mode_top_bound.ge.3) then
       jjb=jj_begin
       jje=jj_end
       IF (pole_nord) jjb=jj_begin+1
       IF (pole_sud)  jje=jj_end-1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)   
       do l=1,llm
        do j=jjb,jje
          zm=0.
          tzon(j,l)=0.
          do i=1,iim
            tzon(j,l)=tzon(j,l)+teta(i,j,l)*masse(i,j,l)
            zm=zm+masse(i,j,l)
          enddo
          tzon(j,l)=tzon(j,l)/zm
        enddo
       enddo
c$OMP END DO NOWAIT
      endif ! of if (mode_top_bound.ge.3)

      if (mode_top_bound.ge.1) then
       ! Apply sponge quenching on vcov:
       jjb=jj_begin
       jje=jj_end
       IF (pole_sud) jje=jj_end-1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
        do j=jjb,jje
          do i=1,iip1
            vcov(i,j,l)=vcov(i,j,l)
     &                  -rdamp(l)*(vcov(i,j,l)-vzon(j,l))
          enddo
        enddo
       enddo
c$OMP END DO NOWAIT

       ! Apply sponge quenching on ucov:
       jjb=jj_begin
       jje=jj_end
       IF (pole_nord) jjb=jj_begin+1
       IF (pole_sud)  jje=jj_end-1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
        do j=jjb,jje
          do i=1,iip1
            ducov(i,j,l)=-rdamp(l)*(ucov(i,j,l)-cu(i,j)*uzon(j,l))
            ucov(i,j,l)=ucov(i,j,l)
     &                  +ducov(i,j,l)
          enddo
       enddo
       enddo
c$OMP END DO NOWAIT
      endif ! of if (mode_top_bound.ge.1)

      if (mode_top_bound.ge.3) then   
       ! Apply sponge quenching on teta:
       jjb=jj_begin
       jje=jj_end
       IF (pole_nord) jjb=jj_begin+1
       IF (pole_sud)  jje=jj_end-1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       do l=1,llm
        do j=jjb,jje
          do i=1,iip1
            teta(i,j,l)=teta(i,j,l)
     &                  -rdamp(l)*(teta(i,j,l)-tzon(j,l))
          enddo
       enddo
       enddo
c$OMP END DO NOWAIT
      endif ! of if (mode_top_bond.ge.3)

      END
