!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/soil.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE soil(ptimestep, knon, ptsrf, ptsoil,
     s          pcapcal, pfluxgrd)

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
c   ptimestep            physical timestep (s)
c   ptsrf(klon)          surface temperature at time-step t (K)
c   ptsoil(klon,nsoilmx) temperature inside the ground (K)
c   pcapcal(klon)        surfacic specific heat (W*m-2*s*K-1)
c   pfluxgrd(klon)       surface diffusive flux from ground (Wm-2)
c   
c=======================================================================
c   declarations:
c   -------------

      use dimphy
      IMPLICIT NONE
#include "YOMCST.h"
#include "dimsoil.h"
#include "clesphys.h"

c-----------------------------------------------------------------------
c  arguments
c  ---------

      REAL ptimestep
      INTEGER knon
      REAL ptsrf(klon),ptsoil(klon,nsoilmx)
      REAL pcapcal(klon),pfluxgrd(klon)

c-----------------------------------------------------------------------
c  local arrays
c  ------------

      INTEGER ig,jk
      REAL zdz2(nsoilmx),z1(klon)
      REAL min_period,dalph_soil
      REAL ztherm_i(klon)

c   local saved variables:
c   ----------------------
      REAL dz1(nsoilmx),dz2(nsoilmx)
      REAL,allocatable :: zc(:,:),zd(:,:)
      REAL lambda
      SAVE dz1,dz2,zc,zd,lambda
      LOGICAL firstcall
      SAVE firstcall

      DATA firstcall/.true./

c-----------------------------------------------------------------------
c   Depthts:
c   --------

      REAL fz,rk,fz1,rk1,rk2
      fz(rk)=fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)
      pfluxgrd(:) = 0.
 
      DO ig = 1, knon
          ztherm_i(ig)   = inertie
      ENDDO

      IF (firstcall) THEN

      allocate(zc(klon,nsoilmx),zd(klon,nsoilmx))

c-----------------------------------------------------------------------
c   ground levels 
c   grnd=z/l where l is the skin depth of the diurnal cycle:
c   --------------------------------------------------------

c VENUS : A REVOIR !!!!
         min_period=20000. ! en secondes
         dalph_soil=2.    ! rapport entre les epaisseurs de 2 couches succ.

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

         DO jk=1,nsoilmx
            rk1=jk
            rk2=jk-1
            dz2(jk)=fz(rk1)-fz(rk2)
         ENDDO
         DO jk=1,nsoilmx-1
            rk1=jk+.5
            rk2=jk-.5
            dz1(jk)=1./(fz(rk1)-fz(rk2))
         ENDDO
         lambda=fz(.5)*dz1(1)
         PRINT*,'full layers, intermediate layers (seconds)'
         DO jk=1,nsoilmx
            rk=jk
            rk1=jk+.5
            rk2=jk-.5
            PRINT *,'fz=',
     .               fz(rk1)*fz(rk2)*3.14,fz(rk)*fz(rk)*3.14
         ENDDO
         firstcall =.false.

      ELSE   !--not firstcall
c-----------------------------------------------------------------------
c   Computation of the soil temperatures using the Cgrd and Dgrd
c  coefficient computed at the previous time-step:
c  -----------------------------------------------

c    surface temperature
         DO ig=1,knon
            ptsoil(ig,1)=(lambda*zc(ig,1)+ptsrf(ig))/
     s      (lambda*(1.-zd(ig,1))+1.)
         ENDDO

c   other temperatures
         DO jk=1,nsoilmx-1
            DO ig=1,knon
               ptsoil(ig,jk+1)=zc(ig,jk)+zd(ig,jk)*ptsoil(ig,jk)
            ENDDO
         ENDDO

      ENDIF !--not firstcall
c-----------------------------------------------------------------------
c   Computation of the Cgrd and Dgrd coefficient for the next step:
c   ---------------------------------------------------------------

      DO jk=1,nsoilmx
         zdz2(jk)=dz2(jk)/ptimestep
      ENDDO

      DO ig=1,knon
         z1(ig)=zdz2(nsoilmx)+dz1(nsoilmx-1)
         zc(ig,nsoilmx-1)=
     $       zdz2(nsoilmx)*ptsoil(ig,nsoilmx)/z1(ig)
         zd(ig,nsoilmx-1)=dz1(nsoilmx-1)/z1(ig)
      ENDDO

      DO jk=nsoilmx-1,2,-1
         DO ig=1,knon
            z1(ig)=1./(zdz2(jk)+dz1(jk-1)+dz1(jk)
     $         *(1.-zd(ig,jk)))
            zc(ig,jk-1)=
     s      (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*zc(ig,jk))
     $          *z1(ig)
            zd(ig,jk-1)=dz1(jk-1)*z1(ig)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c   computation of the surface diffusive flux from ground and
c   calorific capacity of the ground:
c   ---------------------------------

      DO ig=1,knon
         pfluxgrd(ig)=ztherm_i(ig)*dz1(1)*
     s   (zc(ig,1)+(zd(ig,1)-1.)*ptsoil(ig,1))
         pcapcal(ig)=ztherm_i(ig)*
     s   (dz2(1)+ptimestep*(1.-zd(ig,1))*dz1(1))
         z1(ig)=lambda*(1.-zd(ig,1))+1.
         pcapcal(ig)=pcapcal(ig)/z1(ig)
         pfluxgrd(ig) = pfluxgrd(ig)
     s   + pcapcal(ig) * (ptsoil(ig,1) * z1(ig)
     $       - lambda * zc(ig,1)
     $       - ptsrf(ig))
     s   /ptimestep
      ENDDO

      RETURN
      END
