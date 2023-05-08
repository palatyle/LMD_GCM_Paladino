      SUBROUTINE vdifc(ngrid,nlay,nq,co2ice,ppopsk,
     $                ptimestep,pcapcal,lecrit,
     $                pplay,pplev,pzlay,pzlev,pz0,
     $                pu,pv,ph,pq,ptsrf,pemis,pqsurf,
     $                pdufi,pdvfi,pdhfi,pdqfi,pfluxsrf,
     $                pdudif,pdvdif,pdhdif,pdtsrf,pq2,
     $                pdqdif,pdqsdif)
      IMPLICIT NONE

c=======================================================================
c
c   subject:
c   --------
c   Turbulent diffusion (mixing) for potential T, U, V and tracer
c
c   Shema implicite
c   On commence par rajouter au variables x la tendance physique
c   et on resoult en fait:
c      x(t+1) =  x(t) + dt * (dx/dt)phys(t)  +  dt * (dx/dt)difv(t+1)
c
c   author:
c   ------
c      Hourdin/Forget/Fournier
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "surfdat.h"
#include "comgeomfi.h"
#include "tracer.h"

#include "watercap.h"
c
c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL ptimestep
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      REAL pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
      REAL pu(ngrid,nlay),pv(ngrid,nlay),ph(ngrid,nlay)
      REAL ptsrf(ngrid),pemis(ngrid)
      REAL pdufi(ngrid,nlay),pdvfi(ngrid,nlay),pdhfi(ngrid,nlay)
      REAL pfluxsrf(ngrid)
      REAL pdudif(ngrid,nlay),pdvdif(ngrid,nlay),pdhdif(ngrid,nlay)
      REAL pdtsrf(ngrid),pcapcal(ngrid)
      REAL pq2(ngrid,nlay+1)

c    Argument added for condensation:
      REAL co2ice (ngrid), ppopsk(ngrid,nlay)
      logical lecrit
      REAL pz0

c    Traceurs :
      integer nq 
      REAL pqsurf(ngrid,nq)
      real pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq) 
      real pdqdif(ngrid,nlay,nq) 
      real pdqsdif(ngrid,nq) 
      
c   local:
c   ------

      INTEGER ilev,ig,ilay,nlev,ierr

      REAL z4st,zdplanck(ngridmx)
      REAL zkv(ngridmx,nlayermx+1),zkh(ngridmx,nlayermx+1)
	  REAL zcdv(ngridmx),zcdh(ngridmx)
	  REAL zcdv_true(ngridmx),zcdh_true(ngridmx)
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx)
      REAL zh(ngridmx,nlayermx)
      REAL ztsrf2(ngridmx)
      REAL z1(ngridmx),z2(ngridmx)
      REAL za(ngridmx,nlayermx),zb(ngridmx,nlayermx)
      REAL zb0(ngridmx,nlayermx)
      REAL zc(ngridmx,nlayermx),zd(ngridmx,nlayermx)
      REAL zcst1
      REAL zu2

      EXTERNAL SSUM,SCOPY
      REAL SSUM
      LOGICAL firstcall
      SAVE firstcall

c     variable added for CO2 condensation:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      REAL hh , zhcond(ngridmx,nlayermx), tconds
       REAL latcond,tcond1mb
      REAL acond,bcond
      SAVE acond,bcond
      DATA latcond,tcond1mb/5.9e5,136.27/

c    Tracers :
c    ~~~~~~~ 
      INTEGER iq
      REAL zq(ngridmx,nlayermx,nqmx)
      REAL zq1temp(ngridmx)
      REAL rho(ngridmx) ! near surface air density
      REAL qsat(ngridmx)
      DATA firstcall/.true./

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans coefdifv'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF
c        To compute: Tcond= 1./(bcond-acond*log(.0095*p)) (p in pascal)
         bcond=1./tcond1mb
         acond=r/latcond
         PRINT*,'In vdifc: Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
         PRINT*,'acond,bcond,ccond',acond,bcond,acond

         firstcall=.false.
      ENDIF





c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

      nlev=nlay+1

c    ** calcul de rho*dz et dt*rho/dz=dt*rho**2 g/dp
c       avec rho=p/RT=p/ (R Theta) (p/ps)**kappa
c       ----------------------------------------

      DO ilay=1,nlay
         DO ig=1,ngrid
            za(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/g
         ENDDO
      ENDDO

      zcst1=4.*g*ptimestep/(r*r)
      DO ilev=2,nlev-1
         DO ig=1,ngrid
            zb0(ig,ilev)=pplev(ig,ilev)*
     s      (pplev(ig,1)/pplev(ig,ilev))**rcp /
     s      (ph(ig,ilev-1)+ph(ig,ilev))
            zb0(ig,ilev)=zcst1*zb0(ig,ilev)*zb0(ig,ilev)/
     s      (pplay(ig,ilev-1)-pplay(ig,ilev))
         ENDDO
      ENDDO
      DO ig=1,ngrid
		 zb0(ig,1)=ptimestep*pplev(ig,1)/(r*ptsrf(ig))
      ENDDO

c    ** diagnostique pour l'initialisation
c       ----------------------------------

      IF(lecrit) THEN
         ig=ngrid/2+1
         PRINT*,'Pression (mbar) ,altitude (km),u,v,theta, rho dz'
         DO ilay=1,nlay
            WRITE(*,'(6f11.5)')
     s      .01*pplay(ig,ilay),.001*pzlay(ig,ilay),
     s      pu(ig,ilay),pv(ig,ilay),ph(ig,ilay),za(ig,ilay)
         ENDDO
         PRINT*,'Pression (mbar) ,altitude (km),zb'
         DO ilev=1,nlay
            WRITE(*,'(3f15.7)')
     s      .01*pplev(ig,ilev),.001*pzlev(ig,ilev),
     s      zb0(ig,ilev)
         ENDDO
      ENDIF

c     Potential Condensation temperature:
c     -----------------------------------

c     if (callcond) then 
c       DO ilev=1,nlay
c         DO ig=1,ngrid
c           zhcond(ig,ilev) =
c    &  (1./(bcond-acond*log(.0095*pplay(ig,ilev))))/ppopsk(ig,ilev)
c         END DO
c       END DO
c     else
        call zerophys(ngrid*nlay,zhcond)
c     end if


c-----------------------------------------------------------------------
c   2. ajout des tendances physiques
c      -----------------------------

      DO ilev=1,nlay
         DO ig=1,ngrid
            zu(ig,ilev)=pu(ig,ilev)+pdufi(ig,ilev)*ptimestep
            zv(ig,ilev)=pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep
            zh(ig,ilev)=ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep
            zh(ig,ilev)=max(zh(ig,ilev),zhcond(ig,ilev))
         ENDDO
      ENDDO
      if(tracer) then
        DO iq =1, nq
         DO ilev=1,nlay
           DO ig=1,ngrid
              zq(ig,ilev,iq)=pq(ig,ilev,iq)+pdqfi(ig,ilev,iq)*ptimestep
           ENDDO
         ENDDO
        ENDDO
      end if

c-----------------------------------------------------------------------
c   3. schema de turbulence
c      --------------------

c    ** source d'energie cinetique turbulente a la surface
c       (condition aux limites du schema de diffusion turbulente
c       dans la couche limite
c       ---------------------

      CALL vdif_cd( ngrid,nlay,pz0,g,pzlay,pu,pv,ptsrf,ph
     &             ,zcdv_true,zcdh_true)
      DO ig=1,ngrid
        zu2=pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1)
        zcdv(ig)=zcdv_true(ig)*sqrt(zu2)
        zcdh(ig)=zcdh_true(ig)*sqrt(zu2)
      ENDDO

c    ** schema de diffusion turbulente dans la couche limite
c       ---------------------------------------------------- 

        CALL vdif_kc(ptimestep,g,pzlev,pzlay
     &              ,pu,pv,ph,zcdv_true
     &              ,pq2,zkv,zkh)



c    ** diagnostique pour le schema de turbulence
c       -----------------------------------------

      IF(lecrit) THEN
         PRINT*
         PRINT*,'Diagnostic for the vertical turbulent mixing'
         PRINT*,'Cd for momentum and potential temperature'

         PRINT*,zcdv(ngrid/2+1),zcdh(ngrid/2+1)
         PRINT*,'Mixing coefficient for momentum and pot.temp.'
         DO ilev=1,nlay
            PRINT*,zkv(ngrid/2+1,ilev),zkh(ngrid/2+1,ilev)
         ENDDO
      ENDIF




c-----------------------------------------------------------------------
c   4. inversion pour l'implicite sur u
c      --------------------------------

c    ** l'equation est 
c       u(t+1) =  u(t) + dt * {(du/dt)phys}(t)  +  dt * {(du/dt)difv}(t+1)
c       avec
c       /zu/ = u(t) + dt * {(du/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(du/dt)difv}(t+1) = dt * {(d/dz)[ Ku (du/dz) ]}(t+1)
c       donc les entrees sont /zcdv/ pour la condition a la limite sol
c       et /zkv/ = Ku
 
      CALL multipl((nlay-1)*ngrid,zkv(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdv,zb0,zb)

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zu(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $         zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zu(ig,ilay)+
     $         zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zu(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zu(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zu(ig,ilay-1)
         ENDDO
      ENDDO





c-----------------------------------------------------------------------
c   5. inversion pour l'implicite sur v
c      --------------------------------

c    ** l'equation est 
c       v(t+1) =  v(t) + dt * {(dv/dt)phys}(t)  +  dt * {(dv/dt)difv}(t+1)
c       avec
c       /zv/ = v(t) + dt * {(dv/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(dv/dt)difv}(t+1) = dt * {(d/dz)[ Kv (dv/dz) ]}(t+1)
c       donc les entrees sont /zcdv/ pour la condition a la limite sol
c       et /zkv/ = Kv

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zv(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $         zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zv(ig,ilay)+
     $         zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zv(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zv(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zv(ig,ilay-1)
         ENDDO
      ENDDO





c-----------------------------------------------------------------------
c   6. inversion pour l'implicite sur h sans oublier le couplage
c      avec le sol (conduction)
c      ------------------------

c    ** l'equation est 
c       h(t+1) =  h(t) + dt * {(dh/dt)phys}(t)  +  dt * {(dh/dt)difv}(t+1)
c       avec
c       /zh/ = h(t) + dt * {(dh/dt)phys}(t)   (voir paragraphe 2.)
c       et
c       dt * {(dh/dt)difv}(t+1) = dt * {(d/dz)[ Kh (dh/dz) ]}(t+1)
c       donc les entrees sont /zcdh/ pour la condition de raccord au sol
c       et /zkh/ = Kh
c       -------------

      CALL multipl((nlay-1)*ngrid,zkh(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdh,zb0,zb)

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zh(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $         zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zh(ig,ilay)+
     $         zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

c    ** calcul de (d Planck / dT) a la temperature d'interface
c       ------------------------------------------------------

      z4st=4.*5.67e-8*ptimestep
      DO ig=1,ngrid
         zdplanck(ig)=z4st*pemis(ig)*ptsrf(ig)*ptsrf(ig)*ptsrf(ig)
      ENDDO

c    ** calcul de la temperature_d'interface et de sa tendance.
c       on ecrit que la somme des flux est nulle a l'interface
c       a t + \delta t,
c       c'est a dire le flux radiatif a {t + \delta t}
c       + le flux turbulent a {t + \delta t} 
c            qui s'ecrit K (T1-Tsurf) avec T1 = d1 Tsurf + c1
c            (notation K dt = /cpp*b/)        
c       + le flux dans le sol a t
c       + l'evolution du flux dans le sol lorsque la temperature d'interface
c       passe de sa valeur a t a sa valeur a {t + \delta t}.
c       ----------------------------------------------------

      DO ig=1,ngrid
         z1(ig)=pcapcal(ig)*ptsrf(ig)+cpp*zb(ig,1)*zc(ig,1)
     s     +zdplanck(ig)*ptsrf(ig)+ pfluxsrf(ig)*ptimestep
         z2(ig)= pcapcal(ig)+cpp*zb(ig,1)*(1.-zd(ig,1))+zdplanck(ig)
         ztsrf2(ig)=z1(ig)/z2(ig)
         pdtsrf(ig)=(ztsrf2(ig)-ptsrf(ig))/ptimestep

c        Modif speciale CO2 condensation:
c        tconds = 1./(bcond-acond*log(.0095*pplev(ig,1)))
c        if ((callcond).and.
c    &      ((co2ice(ig).ne.0).or.(ztsrf2(ig).lt.tconds)))then
c           zh(ig,1)=zc(ig,1)+zd(ig,1)*tconds
c        else
            zh(ig,1)=zc(ig,1)+zd(ig,1)*ztsrf2(ig)
c        end if
      ENDDO

c    ** et a partir de la temperature au sol on remonte 
c       -----------------------------------------------

      DO ilay=2,nlay
         DO ig=1,ngrid
            hh = max( zh(ig,ilay-1) , zhcond(ig,ilay-1) ) ! modif co2cond
            zh(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*hh
         ENDDO
      ENDDO


c-----------------------------------------------------------------------
c   TRACERS
c   -------

      if(tracer) then

c     Using the wind modified by friction for lifting and  sublimation
c     ----------------------------------------------------------------
        DO ig=1,ngrid
          zu2=zu(ig,1)*zu(ig,1)+zv(ig,1)*zv(ig,1)
          zcdv(ig)=zcdv_true(ig)*sqrt(zu2)
          zcdh(ig)=zcdh_true(ig)*sqrt(zu2)
        ENDDO

c       Calcul du flux vertical au bas de la premiere couche (dust) :
c       -----------------------------------------------------------
        do ig=1,ngridmx  
          rho(ig) = zb0(ig,1) /ptimestep
          zb(ig,1) = 0.
        end do
c       Dust lifting:
        if (lifting) then 
           call dustlift(ngrid,nlay,nq,rho,zcdh_true,zcdh,co2ice,
     $                  pdqsdif)
       
           else
           call zerophys(ngrid*nq,pdqsdif)
        end if

c       OU calcul de la valeur de q a la surface (water)  :
c       ----------------------------------------
        if (water) then 
            call watersat(ngridmx,ptsrf,pplev(1,1),qsat)
        end if

c      Inversion pour l'implicite sur q 
c       --------------------------------
        do iq=1,nq
          CALL multipl((nlay-1)*ngrid,zkh(1,2),zb0(1,2),zb(1,2))

          if ((water).and.(iq.eq.nqmx)) then 
c            This line is required to account for turbulent transport 
c            from surface (e.g. ice) to mid-layer of atmosphere:
             CALL multipl(ngrid,zcdv,zb0,zb(1,1))
             CALL multipl(ngrid,dryness,zb(1,1),zb(1,1)) 
          end if

          DO ig=1,ngrid
               z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
               zc(ig,nlay)=za(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
               zd(ig,nlay)=zb(ig,nlay)*z1(ig)
          ENDDO
  
          DO ilay=nlay-1,2,-1
               DO ig=1,ngrid
                z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $           zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
                zc(ig,ilay)=(za(ig,ilay)*zq(ig,ilay,iq)+
     $           zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
                zd(ig,ilay)=zb(ig,ilay)*z1(ig)
               ENDDO
          ENDDO
          DO ig=1,ngrid
                z1(ig)=1./(za(ig,1)+zb(ig,1)+
     $           zb(ig,2)*(1.-zd(ig,2)))
                zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $         zb(ig,2)*zc(ig,2) +
     $        (-pdqsdif(ig,iq)) *ptimestep) *z1(ig)  !tracer flux from surface
          ENDDO
  
          IF ((water).and.(iq.eq.nqmx)) then 
c           Calculation for turbulent exchange with the surface (for ice)
            DO ig=1,ngrid
              zd(ig,1)=zb(ig,1)*z1(ig)
              zq1temp(ig)=zc(ig,1)+ zd(ig,1)*qsat(ig)

              pdqsdif(ig,nq)=rho(ig)*dryness(ig)*zcdv(ig)
     &                       *(zq1temp(ig)-qsat(ig))
c             write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
            END DO
            DO ig=1,ngrid
              if(.not.watercaptag(ig)) then
              if ((-pdqsdif(ig,nq)*ptimestep).gt.pqsurf(ig,nq)) then
c                 write(*,*)'on sublime plus que qsurf!'
                  pdqsdif(ig,nq) = -pqsurf(ig,nq)/ptimestep
c                 write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
                  z1(ig)=1./(za(ig,1)+ zb(ig,2)*(1.-zd(ig,2)))
                  zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $            zb(ig,2)*zc(ig,2) +
     $            (-pdqsdif(ig,iq)) *ptimestep) *z1(ig)
                  zq1temp(ig)=zc(ig,1)
              endif   
              endif   
c             Starting upward calculations for water :
               zq(ig,1,iq)=zq1temp(ig)
            ENDDO
          ELSE
c           Starting upward calculations for simple mixing of tracer (dust)
            DO ig=1,ngrid
               zq(ig,1,iq)=zc(ig,1)
            ENDDO
          END IF

          DO ilay=2,nlay
             DO ig=1,ngrid
                zq(ig,ilay,iq)=zc(ig,ilay)+zd(ig,ilay)*zq(ig,ilay-1,iq)
             ENDDO
          ENDDO
        enddo
      end if

c-----------------------------------------------------------------------
c   8. calcul final des tendances de la diffusion verticale
c      ----------------------------------------------------

      DO ilev = 1, nlay
         DO ig=1,ngrid
            pdudif(ig,ilev)=(    zu(ig,ilev)-
     $      (pu(ig,ilev)+pdufi(ig,ilev)*ptimestep)    )/ptimestep
            pdvdif(ig,ilev)=(    zv(ig,ilev)-
     $      (pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep)    )/ptimestep
            hh = max(ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep ,
     $           zhcond(ig,ilev))        ! modif co2cond
            pdhdif(ig,ilev)=( zh(ig,ilev)- hh )/ptimestep
         ENDDO
      ENDDO

    
      if (tracer) then 
        DO iq = 1, nq
          DO ilev = 1, nlay
            DO ig=1,ngrid
              pdqdif(ig,ilev,iq)=(zq(ig,ilev,iq)-
     $      (pq(ig,ilev,iq) + pdqfi(ig,ilev,iq)*ptimestep))/ptimestep
            ENDDO
          ENDDO
        ENDDO
      end if

c    ** diagnostique final 
c       ------------------

      IF(lecrit) THEN
         PRINT*,'In vdif'
         PRINT*,'Ts (t) and Ts (t+st)'
         WRITE(*,'(a10,3a15)')
     s   'theta(t)','theta(t+dt)','u(t)','u(t+dt)'
         PRINT*,ptsrf(ngrid/2+1),ztsrf2(ngrid/2+1)
         DO ilev=1,nlay
            WRITE(*,'(4f15.7)')
     s      ph(ngrid/2+1,ilev),zh(ngrid/2+1,ilev),
     s      pu(ngrid/2+1,ilev),zu(ngrid/2+1,ilev)

         ENDDO
      ENDIF

      RETURN
      END
