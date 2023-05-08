      SUBROUTINE soil(ngrid,nsoil,firstcall,ptherm_i,
     s          ptimestep,ptsrf,ptsoil,
     s          pcapcal,pfluxgrd)
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  Frederic Hourdin     30/01/92
c   -------
c
c   objet:  computation of : the soil temperature evolution
c   ------                   the surfacic heat capacity "Capcal"
c                            the surface conduction flux pcapcal
c
c
c   Method: implicit time integration
c   -------
c   Consecutive ground temperatures are related by:
c           T(k+1) = C(k) + D(k)*T(k)  (1)
c   the coefficients C and D are computed at the t-dt time-step.
c   Routine structure:
c   1)new temperatures are computed  using (1)
c   2)C and D coefficients are computed from the new temperature
c     profile for the t+dt time-step
c   3)the coefficients A and B are computed where the diffusive
c     fluxes at the t+dt time-step is given by
c            Fdiff = A + B Ts(t+dt)
c     or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
c            with F0 = A + B (Ts(t))
c                 Capcal = B*dt
c           
c   Interface:
c   ----------
c
c   Arguments:
c   ----------
c   ngird               number of grid-points
c   ptimestep              physical timestep (s)
c   pto(ngrid,nsoil)     temperature at time-step t (K)
c   ptn(ngrid,nsoil)     temperature at time step t+dt (K)
c   pcapcal(ngrid)      specific heat (W*m-2*s*K-1)
c   pfluxgrd(ngrid)      surface diffusive flux from ground (Wm-2)
c   
c=======================================================================
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"

c-----------------------------------------------------------------------
c  arguments
c  ---------

      INTEGER ngrid,nsoil
      REAL ptimestep
      REAL ptsrf(ngridmx),ptsoil(ngrid,nsoilmx),ptherm_i(ngridmx)
      REAL pcapcal(ngrid),pfluxgrd(ngrid)
      LOGICAL firstcall

c-----------------------------------------------------------------------
c  local arrays
c  ------------

      INTEGER ig,jk
      REAL za(ngridmx),zb(ngridmx)
      REAL zdz2(nsoilmx),z1(ngridmx)
      REAL min_period,dalph_soil

c   local saved variables:
c   ----------------------
      REAL dz1(nsoilmx),dz2(nsoilmx)
      REAL zc(ngridmx,nsoilmx),zd(ngridmx,nsoilmx)
      REAL lambda
      SAVE dz1,dz2,zc,zd,lambda

c-----------------------------------------------------------------------
c   Depthts:
c   --------

      REAL fz,rk,fz1,rk1,rk2
      fz(rk)=fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)

      IF (firstcall) THEN

c-----------------------------------------------------------------------
c   ground levels 
c   grnd=z/l where l is the skin depth of the diurnal cycle:
c   --------------------------------------------------------

         min_period=887.75
c        min_period=1800.
c        min_period=20000.
         dalph_soil=2.

         OPEN(99,file='soil.def',status='old',form='formatted',err=9999)
         READ(99,*) min_period
         READ(99,*) dalph_soil
         PRINT*,'Discretization for the soil model'
         PRINT*,'First level e-folding depth',min_period,
     s   '   dalph',dalph_soil
         CLOSE(99)
9999     CONTINUE

c   la premiere couche represente un dixieme de cycle diurne
         fz1=sqrt(min_period/3.14)

         DO jk=1,nsoil
            rk1=jk
            rk2=jk-1
            dz2(jk)=fz(rk1)-fz(rk2)
         ENDDO
         DO jk=1,nsoil-1
            rk1=jk+.5
            rk2=jk-.5
            dz1(jk)=1./(fz(rk1)-fz(rk2))
         ENDDO
         lambda=fz(.5)*dz1(1)
         PRINT*,'full layers, intermediate layers (secoonds)'
         DO jk=1,nsoil
            rk=jk
            rk1=jk+.5
            rk2=jk-.5
            PRINT*,fz(rk1)*fz(rk2)*3.14,
     s      fz(rk)*fz(rk)*3.14
         ENDDO

c   Initialisations:
c   ----------------

      ELSE
c-----------------------------------------------------------------------
c   Computation of the soil temperatures using the Cgrd and Dgrd
c  coefficient computed at the previous time-step:
c  -----------------------------------------------

c    surface temperature
         DO ig=1,ngrid
            ptsoil(ig,1)=(lambda*zc(ig,1)+ptsrf(ig))/
     s      (lambda*(1.-zd(ig,1))+1.)
         ENDDO

c   other temperatures
         DO jk=1,nsoil-1
            DO ig=1,ngrid
               ptsoil(ig,jk+1)=zc(ig,jk)+zd(ig,jk)*ptsoil(ig,jk)
            ENDDO
         ENDDO

      ENDIF
c-----------------------------------------------------------------------
c   Computation of the Cgrd and Dgrd coefficient for the next step:
c   ---------------------------------------------------------------

      DO jk=1,nsoil
         zdz2(jk)=dz2(jk)/ptimestep
      ENDDO

      DO ig=1,ngrid
         z1(ig)=zdz2(nsoil)+dz1(nsoil-1)
         zc(ig,nsoil-1)=zdz2(nsoil)*ptsoil(ig,nsoil)/z1(ig)
         zd(ig,nsoil-1)=dz1(nsoil-1)/z1(ig)
      ENDDO

      DO jk=nsoil-1,2,-1
         DO ig=1,ngrid
            z1(ig)=1./(zdz2(jk)+dz1(jk-1)+dz1(jk)*(1.-zd(ig,jk)))
            zc(ig,jk-1)=
     s      (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*zc(ig,jk))*z1(ig)
            zd(ig,jk-1)=dz1(jk-1)*z1(ig)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c   computation of the surface diffusive flux from ground and
c   calorific capacity of the ground:
c   ---------------------------------

      DO ig=1,ngrid
         pfluxgrd(ig)=ptherm_i(ig)*dz1(1)*
     s   (zc(ig,1)+(zd(ig,1)-1.)*ptsoil(ig,1))
         pcapcal(ig)=ptherm_i(ig)*
     s   (dz2(1)+ptimestep*(1.-zd(ig,1))*dz1(1))
         z1(ig)=lambda*(1.-zd(ig,1))+1.
         pcapcal(ig)=pcapcal(ig)/z1(ig)
         pfluxgrd(ig)=pfluxgrd(ig)
     s   +pcapcal(ig)*(ptsoil(ig,1)*z1(ig)-lambda*zc(ig,1)-ptsrf(ig))
     s   /ptimestep
      ENDDO

      RETURN
      END
