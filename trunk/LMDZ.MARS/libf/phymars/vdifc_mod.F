      MODULE vdifc_mod
      
      IMPLICIT NONE
      
      CONTAINS
      
      SUBROUTINE vdifc(ngrid,nlay,nq,co2ice,ppopsk,
     $                ptimestep,pcapcal,lecrit,
     $                pplay,pplev,pzlay,pzlev,pz0,
     $                pu,pv,ph,pq,ptsrf,pemis,pqsurf,
     $                pdufi,pdvfi,pdhfi,pdqfi,pfluxsrf,
     $                pdudif,pdvdif,pdhdif,pdtsrf,pq2,
     $                pdqdif,pdqsdif,wstar,zcdv_true,zcdh_true,
     $                hfmax,pcondicea_co2microp,sensibFlux,
     $                dustliftday,local_time,watercap, dwatercap_dif)

      use tracer_mod, only: noms, igcm_dust_mass, igcm_dust_number,
     &                      igcm_dust_submicron, igcm_h2o_vap,
     &                      igcm_h2o_ice, alpha_lift,
     &                      igcm_stormdust_mass, igcm_stormdust_number
      use surfdat_h, only: watercaptag, frost_albedo_threshold, dryness
      USE comcstfi_h, ONLY: cpp, r, rcp, g
      use watersat_mod, only: watersat
      use turb_mod, only: turb_resolved, ustar, tstar
      use compute_dtau_mod, only: ti_injection_sol,tf_injection_sol
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

      include "callkeys.h"
      include "microphys.h"

c
c   arguments:
c   ----------

      INTEGER,INTENT(IN) :: ngrid,nlay
      REAL,INTENT(IN) :: ptimestep
      REAL,INTENT(IN) :: pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      REAL,INTENT(IN) :: pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
      REAL,INTENT(IN) :: pu(ngrid,nlay),pv(ngrid,nlay)
      REAL,INTENT(IN) :: ph(ngrid,nlay)
      REAL,INTENT(IN) :: ptsrf(ngrid),pemis(ngrid)
      REAL,INTENT(IN) :: pdufi(ngrid,nlay),pdvfi(ngrid,nlay)
      REAL,INTENT(IN) :: pdhfi(ngrid,nlay)
      REAL,INTENT(IN) :: pfluxsrf(ngrid)
      REAL,INTENT(OUT) :: pdudif(ngrid,nlay),pdvdif(ngrid,nlay)
      REAL,INTENT(OUT) :: pdtsrf(ngrid),pdhdif(ngrid,nlay)
      REAL,INTENT(IN) :: pcapcal(ngrid)
      REAL,INTENT(INOUT) :: pq2(ngrid,nlay+1)

c    Argument added for condensation:
      REAL,INTENT(IN) :: co2ice (ngrid), ppopsk(ngrid,nlay)
      logical,INTENT(IN) :: lecrit
      REAL,INTENT(IN) :: pcondicea_co2microp(ngrid,nlay)! tendency due to CO2 condensation (kg/kg.s-1)
      
      REAL,INTENT(IN) :: pz0(ngrid) ! surface roughness length (m)

c    Argument added to account for subgrid gustiness :

      REAL,INTENT(IN) :: wstar(ngrid), hfmax(ngrid)!, zi(ngrid)

c    Traceurs :
      integer,intent(in) :: nq 
      REAL,INTENT(IN) :: pqsurf(ngrid,nq)
      real,intent(in) :: pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq) 
      real,intent(out) :: pdqdif(ngrid,nlay,nq) 
      real,intent(out) :: pdqsdif(ngrid,nq)
      REAL,INTENT(in) :: dustliftday(ngrid)
      REAL,INTENT(in) :: local_time(ngrid)
      
c   local:
c   ------

      REAL :: pt(ngrid,nlay)
 
      INTEGER ilev,ig,ilay,nlev

      REAL z4st,zdplanck(ngrid)
      REAL zkv(ngrid,nlay+1),zkh(ngrid,nlay+1)
      REAL zkq(ngrid,nlay+1)
      REAL zcdv(ngrid),zcdh(ngrid)
      REAL zcdv_true(ngrid),zcdh_true(ngrid)    ! drag coeff are used by the LES to recompute u* and hfx
      REAL zu(ngrid,nlay),zv(ngrid,nlay)
      REAL zh(ngrid,nlay)
      REAL ztsrf2(ngrid)
      REAL z1(ngrid),z2(ngrid)
      REAL za(ngrid,nlay),zb(ngrid,nlay)
      REAL zb0(ngrid,nlay)
      REAL zc(ngrid,nlay),zd(ngrid,nlay)
      REAL zcst1
      REAL zu2(ngrid)

      EXTERNAL SSUM,SCOPY
      REAL SSUM
      LOGICAL,SAVE :: firstcall=.true.


c     variable added for CO2 condensation:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      REAL hh , zhcond(ngrid,nlay)
      REAL,PARAMETER :: latcond=5.9e5
      REAL,PARAMETER :: tcond1mb=136.27
      REAL,SAVE :: acond,bcond
      
c     For latent heat release from ground water ice sublimation    
      REAL tsrf_lh(ngrid) ! temporary surface temperature
      REAL lh ! latent heat, formulation given in the Technical Document: 
              ! "Modeling water ice sublimation under Phoenix-like conditions", Montmessin et al. 2004  
!      REAL tsrf_lw(ngrid)
!      REAL alpha
!      REAL T1,T2
!      SAVE T1,T2
!      DATA T1,T2/-604.3,1080.7/ ! zeros of latent heat equation for ice

c    Tracers :
c    ~~~~~~~ 
      INTEGER iq
      REAL zq(ngrid,nlay,nq)
      REAL zq1temp(ngrid)
      REAL rho(ngrid) ! near surface air density
      REAL qsat(ngrid)

      REAL kmixmin

c    Argument added for surface water ice budget:
      REAL,INTENT(IN) :: watercap(ngrid)
      REAL,INTENT(OUT) :: dwatercap_dif(ngrid)

c    Mass-variation scheme :
c    ~~~~~~~

      INTEGER j,l
      REAL zcondicea(ngrid,nlay)
      REAL zt(ngrid,nlay),ztcond(ngrid,nlay+1)
      REAL betam(ngrid,nlay),dmice(ngrid,nlay)
      REAL pdtc(ngrid,nlay)
      REAL zhs(ngrid,nlay)
      REAL,SAVE :: ccond

c     Theta_m formulation for mass-variation scheme :
c     ~~~~~~~

      INTEGER,SAVE :: ico2
      INTEGER llnt(ngrid)
      REAL,SAVE :: m_co2, m_noco2, A , B
      REAL vmr_co2(ngrid,nlay)
      REAL qco2,mmean

      REAL,INTENT(OUT) :: sensibFlux(ngrid)

c    ** un petit test de coherence
c       --------------------------

      ! AS: OK firstcall absolute
      IF (firstcall) THEN
c        To compute: Tcond= 1./(bcond-acond*log(.0095*p)) (p in pascal)
         bcond=1./tcond1mb
         acond=r/latcond
         ccond=cpp/(g*latcond)
         PRINT*,'In vdifc: Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
         PRINT*,'          acond,bcond,ccond',acond,bcond,ccond


         ico2=0

         if (tracer) then
c          Prepare Special treatment if one of the tracer is CO2 gas
           do iq=1,nq
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)
c               Compute A and B coefficient use to compute
c               mean molecular mass Mair defined by
c               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
c               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             endif
           enddo
         end if

        firstcall=.false.
      ENDIF


c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

      nlev=nlay+1

      ! initialize output tendencies to zero:
      pdudif(1:ngrid,1:nlay)=0
      pdvdif(1:ngrid,1:nlay)=0
      pdhdif(1:ngrid,1:nlay)=0
      pdtsrf(1:ngrid)=0
      pdqdif(1:ngrid,1:nlay,1:nq)=0
      pdqsdif(1:ngrid,1:nq)=0
      dwatercap_dif(1:ngrid)=0

c    ** calcul de rho*dz et dt*rho/dz=dt*rho**2 g/dp
c       avec rho=p/RT=p/ (R Theta) (p/ps)**kappa
c       ----------------------------------------

      DO ilay=1,nlay
         DO ig=1,ngrid
            za(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/g
! Mass variation scheme:
            betam(ig,ilay)=-za(ig,ilay)*latcond/(cpp*ppopsk(ig,ilay))
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

c     -----------------------------------
c     Potential Condensation temperature:
c     -----------------------------------

c     Compute CO2 Volume mixing ratio
c     -------------------------------
      if (callcond.and.(ico2.ne.0)) then
         DO ilev=1,nlay
            DO ig=1,ngrid
              qco2=MAX(1.E-30
     &             ,pq(ig,ilev,ico2)+pdqfi(ig,ilev,ico2)*ptimestep)
c             Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
              mmean=1/(A*qco2 +B)
              vmr_co2(ig,ilev) = qco2*mmean/m_co2
            ENDDO
         ENDDO
      else
         DO ilev=1,nlay
            DO ig=1,ngrid
              vmr_co2(ig,ilev)=0.95
            ENDDO
         ENDDO
      end if

c     forecast of atmospheric temperature zt and frost temperature ztcond
c     --------------------------------------------------------------------

      if (callcond) then 
        DO ilev=1,nlay
          DO ig=1,ngrid
              ztcond(ig,ilev)=
     &      1./(bcond-acond*log(.01*vmr_co2(ig,ilev)*pplay(ig,ilev)))
            if (pplay(ig,ilev).lt.1e-4) ztcond(ig,ilev)=0.0 !mars Monica
!            zhcond(ig,ilev) =
!     &  (1./(bcond-acond*log(.0095*pplay(ig,ilev))))/ppopsk(ig,ilev)
              zhcond(ig,ilev) = ztcond(ig,ilev)/ppopsk(ig,ilev)
          END DO
        END DO
        ztcond(:,nlay+1)=ztcond(:,nlay)
      else
         zhcond(:,:) = 0
         ztcond(:,:) = 0
      end if


c-----------------------------------------------------------------------
c   2. ajout des tendances physiques
c      -----------------------------

      DO ilev=1,nlay
         DO ig=1,ngrid
            zu(ig,ilev)=pu(ig,ilev)+pdufi(ig,ilev)*ptimestep
            zv(ig,ilev)=pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep
            zh(ig,ilev)=ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep
!            zh(ig,ilev)=max(zh(ig,ilev),zhcond(ig,ilev))
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

      CALL vdif_cd(ngrid,nlay,pz0,g,pzlay,pu,pv,wstar,ptsrf,ph
     &             ,zcdv_true,zcdh_true)

        zu2(:)=pu(:,1)*pu(:,1)+pv(:,1)*pv(:,1)

        IF (callrichsl) THEN
          zcdv(:)=zcdv_true(:)*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)
          zcdh(:)=zcdh_true(:)*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)

           ustar(:)=sqrt(zcdv_true(:))*sqrt(zu2(:)+
     &     (log(1.+0.7*wstar(:) + 2.3*wstar(:)**2))**2)

           tstar(:)=0.
           DO ig=1,ngrid
              IF (zcdh_true(ig) .ne. 0.) THEN      ! When Cd=Ch=0, u*=t*=0
                 tstar(ig)=(ph(ig,1)-ptsrf(ig))*zcdh(ig)/ustar(ig)
              ENDIF
           ENDDO

        ELSE
           zcdv(:)=zcdv_true(:)*sqrt(zu2(:))     ! 1 / bulk aerodynamic momentum conductance 
           zcdh(:)=zcdh_true(:)*sqrt(zu2(:))     ! 1 / bulk aerodynamic heat conductance
           ustar(:)=sqrt(zcdv_true(:))*sqrt(zu2(:))
           tstar(:)=(ph(:,1)-ptsrf(:))*zcdh_true(:)/sqrt(zcdv_true(:))
        ENDIF

! Some usefull diagnostics for the new surface layer parametrization :

!        call WRITEDIAGFI(ngrid,'vdifc_zcdv_true',
!     &              'momentum drag','no units',
!     &                         2,zcdv_true)
!        call WRITEDIAGFI(ngrid,'vdifc_zcdh_true',
!     &              'heat drag','no units',
!     &                         2,zcdh_true)
!        call WRITEDIAGFI(ngrid,'vdifc_ust',
!     &              'friction velocity','m/s',2,ust)
!       call WRITEDIAGFI(ngrid,'vdifc_tst',
!     &              'friction temperature','K',2,tst)
!        call WRITEDIAGFI(ngrid,'vdifc_zcdv',
!     &              'aerodyn momentum conductance','m/s',
!     &                         2,zcdv)
!        call WRITEDIAGFI(ngrid,'vdifc_zcdh',
!     &              'aerodyn heat conductance','m/s',
!     &                         2,zcdh)

c    ** schema de diffusion turbulente dans la couche limite
c       ---------------------------------------------------- 
       IF (.not. callyamada4) THEN

       CALL vdif_kc(ngrid,nlay,nq,ptimestep,g,pzlev,pzlay
     &              ,pu,pv,ph,zcdv_true
     &              ,pq2,zkv,zkh,zq)

      ELSE

      pt(:,:)=ph(:,:)*ppopsk(:,:)
      CALL yamada4(ngrid,nlay,nq,ptimestep,g,r,pplev,pt
     s   ,pzlev,pzlay,pu,pv,ph,pq,zcdv_true,pq2,zkv,zkh,zkq,ustar
     s   ,9)
      ENDIF

      if ((doubleq).and.(ngrid.eq.1)) then
        kmixmin = 80. !80.! minimum eddy mix coeff in 1D
        do ilev=1,nlay
          do ig=1,ngrid
           zkh(ig,ilev) = max(kmixmin,zkh(ig,ilev))
           zkv(ig,ilev) = max(kmixmin,zkv(ig,ilev))
          end do
        end do
      end if

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
 
      zb(1:ngrid,2:nlay)=zkv(1:ngrid,2:nlay)*zb0(1:ngrid,2:nlay)
      zb(1:ngrid,1)=zcdv(1:ngrid)*zb0(1:ngrid,1)

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

c Mass variation scheme:
      zb(1:ngrid,2:nlay)=zkh(1:ngrid,2:nlay)*zb0(1:ngrid,2:nlay)
      zb(1:ngrid,1)=zcdh(1:ngrid)*zb0(1:ngrid,1)

c on initialise dm c
      
      pdtc(:,:)=0.
      zt(:,:)=0.
      dmice(:,:)=0.

c    ** calcul de (d Planck / dT) a la temperature d'interface
c       ------------------------------------------------------

      z4st=4.*5.67e-8*ptimestep
      IF (tke_heat_flux .eq. 0.) THEN
      DO ig=1,ngrid
         zdplanck(ig)=z4st*pemis(ig)*ptsrf(ig)*ptsrf(ig)*ptsrf(ig)
      ENDDO
      ELSE
         zdplanck(:)=0.
      ENDIF

! calcul de zc et zd pour la couche top en prenant en compte le terme
! de variation de masse (on fait une boucle pour que \E7a converge)

! Identification des points de grilles qui ont besoin de la correction

      llnt(:)=1
      IF (.not.turb_resolved) THEN
      IF (callcond) THEN
       DO ig=1,ngrid
         DO l=1,nlay
            if(zh(ig,l) .lt. zhcond(ig,l)) then
               llnt(ig)=300  
! 200 and 100 do not go beyond month 9 with normal dissipation
               goto 5
            endif
         ENDDO
5      continue
       ENDDO
      ENDIF

      ENDIF

      DO ig=1,ngrid

! Initialization of z1 and zd, which do not depend on dmice

      z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
      zd(ig,nlay)=zb(ig,nlay)*z1(ig)

      DO ilay=nlay-1,1,-1
          z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $        zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
          zd(ig,ilay)=zb(ig,ilay)*z1(ig)
      ENDDO

! Convergence loop

      DO j=1,llnt(ig)

            z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
            zc(ig,nlay)=za(ig,nlay)*zh(ig,nlay)
     &      -betam(ig,nlay)*dmice(ig,nlay)
            zc(ig,nlay)=zc(ig,nlay)*z1(ig)
!            zd(ig,nlay)=zb(ig,nlay)*z1(ig)

! calcul de zc et zd pour les couches du haut vers le bas

             DO ilay=nlay-1,1,-1
               z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $            zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
               zc(ig,ilay)=(za(ig,ilay)*zh(ig,ilay)+
     $            zb(ig,ilay+1)*zc(ig,ilay+1)-
     $            betam(ig,ilay)*dmice(ig,ilay))*z1(ig)
!               zd(ig,ilay)=zb(ig,ilay)*z1(ig)
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

         z1(ig)=pcapcal(ig)*ptsrf(ig)+cpp*zb(ig,1)*zc(ig,1)
     s     +zdplanck(ig)*ptsrf(ig)+ pfluxsrf(ig)*ptimestep
         z2(ig)= pcapcal(ig)+cpp*zb(ig,1)*(1.-zd(ig,1))+zdplanck(ig)
         ztsrf2(ig)=z1(ig)/z2(ig)
!         pdtsrf(ig)=(ztsrf2(ig)-ptsrf(ig))/ptimestep  !incremented outside loop
            zhs(ig,1)=zc(ig,1)+zd(ig,1)*ztsrf2(ig)

c    ** et a partir de la temperature au sol on remonte 
c       -----------------------------------------------

            DO ilay=2,nlay
               zhs(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zhs(ig,ilay-1)
            ENDDO
            DO ilay=1,nlay
               zt(ig,ilay)=zhs(ig,ilay)*ppopsk(ig,ilay)
            ENDDO

c      Condensation/sublimation in the atmosphere
c      ------------------------------------------
c      (computation of zcondicea and dmice)

      IF (.NOT. co2clouds) then
       DO l=nlay , 1, -1
           IF(zt(ig,l).LT.ztcond(ig,l)) THEN
               pdtc(ig,l)=(ztcond(ig,l) - zt(ig,l))/ptimestep
               zcondicea(ig,l)=(pplev(ig,l)-pplev(ig,l+1))
     &                        *ccond*pdtc(ig,l)
              dmice(ig,l)= dmice(ig,l) + zcondicea(ig,l)*ptimestep
            END IF
       ENDDO
      ELSE
       DO l=nlay , 1, -1
         zcondicea(ig,l)= 0.!pcondicea_co2microp(ig,l)* 
c     &                        (pplev(ig,l) - pplev(ig,l+1))/g
         dmice(ig,l)= 0.!dmice(ig,l) + zcondicea(ig,l)*ptimestep
         pdtc(ig,l)=0.
       ENDDO    
      ENDIF
      
       ENDDO!of Do j=1,XXX 

      ENDDO   !of Do ig=1,ngrid

      pdtsrf(:)=(ztsrf2(:)-ptsrf(:))/ptimestep
      
      DO ig=1,ngrid  ! computing sensible heat flux (atm => surface)
         sensibFlux(ig)=cpp*zb(ig,1)/ptimestep*(zhs(ig,1)-ztsrf2(ig))
      ENDDO

c-----------------------------------------------------------------------
c   TRACERS
c   -------

      if(tracer) then

c     Using the wind modified by friction for lifting and  sublimation
c     ----------------------------------------------------------------

!     This is computed above and takes into account surface-atmosphere flux 
!     enhancement by subgrid gustiness and atmospheric-stability related
!     variations of transfer coefficients.

!        DO ig=1,ngrid
!          zu2(ig)=zu(ig,1)*zu(ig,1)+zv(ig,1)*zv(ig,1)
!          zcdv(ig)=zcdv_true(ig)*sqrt(zu2(ig))
!          zcdh(ig)=zcdh_true(ig)*sqrt(zu2(ig))
!        ENDDO

c       Calcul du flux vertical au bas de la premiere couche (dust) :
c       -----------------------------------------------------------
        do ig=1,ngrid  
          rho(ig) = zb0(ig,1) /ptimestep
c          zb(ig,1) = 0.
        end do
c       Dust lifting:
        if (lifting) then
#ifndef MESOSCALE
           if (doubleq.AND.submicron) then
             do ig=1,ngrid
c              if(co2ice(ig).lt.1) then
                 pdqsdif(ig,igcm_dust_mass) =
     &             -alpha_lift(igcm_dust_mass)  
                 pdqsdif(ig,igcm_dust_number) = 
     &             -alpha_lift(igcm_dust_number)  
                 pdqsdif(ig,igcm_dust_submicron) =
     &             -alpha_lift(igcm_dust_submicron)
c              end if
             end do
           else if (doubleq) then
             if (dustinjection.eq.0) then !injection scheme 0 (old) 
             				  !or 2 (injection in CL)
              do ig=1,ngrid
               if(co2ice(ig).lt.1) then ! pas de soulevement si glace CO2
                 pdqsdif(ig,igcm_dust_mass) =
     &             -alpha_lift(igcm_dust_mass)  
                 pdqsdif(ig,igcm_dust_number) = 
     &             -alpha_lift(igcm_dust_number) 
               end if 
              end do
             elseif(dustinjection.eq.1)then ! dust injection scheme = 1 injection from surface
              do ig=1,ngrid
               if(co2ice(ig).lt.1) then ! pas de soulevement si glace CO2
                IF((ti_injection_sol.LE.local_time(ig)).and.
     &            	(local_time(ig).LE.tf_injection_sol)) THEN
     		  if (rdstorm) then !Rocket dust storm scheme
             	    	pdqsdif(ig,igcm_stormdust_mass) = 
     &             		-alpha_lift(igcm_stormdust_mass)
     &				*dustliftday(ig)
               		pdqsdif(ig,igcm_stormdust_number) = 
     &             		-alpha_lift(igcm_stormdust_number)
     &				*dustliftday(ig)
                 	pdqsdif(ig,igcm_dust_mass)= 0.
                        pdqsdif(ig,igcm_dust_number)= 0.
                  else
                 	pdqsdif(ig,igcm_dust_mass)=
     &            		-dustliftday(ig)*
     &            		alpha_lift(igcm_dust_mass)               
                 	pdqsdif(ig,igcm_dust_number)= 
     &             		-dustliftday(ig)*
     &            		alpha_lift(igcm_dust_number)
		  endif
             	  if (submicron) then
                   	pdqsdif(ig,igcm_dust_submicron) = 0.
                  endif ! if (submicron)
     		ELSE ! outside dust injection time frame
                  pdqsdif(ig,igcm_dust_mass)= 0.
                  pdqsdif(ig,igcm_dust_number)= 0.
                  if (rdstorm) then      
                 	pdqsdif(ig,igcm_stormdust_mass)= 0.
                        pdqsdif(ig,igcm_stormdust_number)= 0.
                  end if
     		ENDIF
               
                end if ! of if(co2ice(ig).lt.1)
              end do
             endif ! end if dustinjection
           else if (submicron) then
             do ig=1,ngrid
                 pdqsdif(ig,igcm_dust_submicron) =
     &             -alpha_lift(igcm_dust_submicron)
             end do
           else
#endif
            call dustlift(ngrid,nlay,nq,rho,zcdh_true,zcdh,co2ice,
     &                    pdqsdif)
#ifndef MESOSCALE
           endif !doubleq.AND.submicron
#endif
        else
           pdqsdif(1:ngrid,1:nq) = 0.
        end if

c       OU calcul de la valeur de q a la surface (water)  :
c       ----------------------------------------
        if (water) then 
            call watersat(ngrid,ptsrf,pplev(1,1),qsat)
        end if

c      Inversion pour l'implicite sur q 
c       --------------------------------
        do iq=1,nq  !for all tracers including stormdust
          zb(1:ngrid,2:nlay)=zkh(1:ngrid,2:nlay)*zb0(1:ngrid,2:nlay)

          if ((water).and.(iq.eq.igcm_h2o_vap)) then 
c            This line is required to account for turbulent transport 
c            from surface (e.g. ice) to mid-layer of atmosphere:
             zb(1:ngrid,1)=zcdv(1:ngrid)*zb0(1:ngrid,1)
             zb(1:ngrid,1)=dryness(1:ngrid)*zb(1:ngrid,1)
          else ! (re)-initialize zb(:,1)
             zb(1:ngrid,1)=0
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

          if (water.and.(iq.eq.igcm_h2o_ice)) then
            ! special case for water ice tracer: do not include
            ! h2o ice tracer from surface (which is set when handling
           ! h2o vapour case (see further down).
            DO ig=1,ngrid
                z1(ig)=1./(za(ig,1)+zb(ig,1)+
     $           zb(ig,2)*(1.-zd(ig,2)))
                zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $         zb(ig,2)*zc(ig,2)) *z1(ig)
            ENDDO
          else ! general case
            DO ig=1,ngrid
                z1(ig)=1./(za(ig,1)+zb(ig,1)+
     $           zb(ig,2)*(1.-zd(ig,2)))
                zc(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $         zb(ig,2)*zc(ig,2) +
     $        (-pdqsdif(ig,iq)) *ptimestep) *z1(ig)  !tracer flux from surface
            ENDDO
          endif ! of if (water.and.(iq.eq.igcm_h2o_ice))
  
          IF ((water).and.(iq.eq.igcm_h2o_vap)) then 
c           Calculation for turbulent exchange with the surface (for ice)
            DO ig=1,ngrid
              zd(ig,1)=zb(ig,1)*z1(ig)
              zq1temp(ig)=zc(ig,1)+ zd(ig,1)*qsat(ig)

              pdqsdif(ig,igcm_h2o_ice)=rho(ig)*dryness(ig)*zcdv(ig)
     &                       *(zq1temp(ig)-qsat(ig))
c             write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
            END DO

            DO ig=1,ngrid
              if ((-pdqsdif(ig,igcm_h2o_ice)*ptimestep)
     &             .gt.(pqsurf(ig,igcm_h2o_ice))) then
c              write(*,*)'subliming more than available frost:  qsurf!'
                if(watercaptag(ig)) then
c              dwatercap_dif same sign as pdqsdif (pos. toward ground)
                   dwatercap_dif(ig) = -pdqsdif(ig,igcm_h2o_ice)
     &                         - pqsurf(ig,igcm_h2o_ice)/ptimestep
                   pdqsdif(ig,igcm_h2o_ice)=
     &                         -pqsurf(ig,igcm_h2o_ice)/ptimestep
c                  write(*,*) "subliming more than available frost:  qsurf!"
                else  !  (case not.watercaptag(ig)) 
                   pdqsdif(ig,igcm_h2o_ice)=
     &                         -pqsurf(ig,igcm_h2o_ice)/ptimestep
                   dwatercap_dif(ig) = 0.0

                endif   
c                 write(*,*)'flux vers le sol=',pdqsdif(ig,nq)
                  z1(ig)=1./(za(ig,1)+ zb(ig,2)*(1.-zd(ig,2)))
                  zc(ig,1)=(za(ig,1)*zq(ig,1,igcm_h2o_vap)+
     $            zb(ig,2)*zc(ig,2) + (dwatercap_dif(ig)-
c     $            zb(ig,2)*zc(ig,2) + (-
     $            pdqsdif(ig,igcm_h2o_ice)) *ptimestep) *z1(ig)
                  zq1temp(ig)=zc(ig,1)
              endif ! if (.not.watercaptag(ig))
c             Starting upward calculations for water :
              zq(ig,1,igcm_h2o_vap)=zq1temp(ig)
              
!c             Take into account H2O latent heat in surface energy budget
!c             We solve dT/dt = (2834.3-0.28*(T-To)-0.004*(T-To)^2)*1e3*iceflux/cpp
!              tsrf_lw(ig) = ptsrf(ig) + pdtsrf(ig) *ptimestep
!             
!              alpha = exp(-4*abs(T1-T2)*pdqsdif(ig,igcm_h2o_ice)
!     &            *ptimestep/pcapcal(ig))
!
!              tsrf_lw(ig) = (tsrf_lw(ig)*(T2-alpha*T1)+T1*T2*(alpha-1))
!     &                      /(tsrf_lw(ig)*(1-alpha)+alpha*T2-T1)  ! surface temperature at t+1
!
!              pdtsrf(ig) = (tsrf_lw(ig)-ptsrf(ig))/ptimestep
c             Take into account the H2O latent heat impact on the surface temperature
              if (latentheat_surfwater) then
                tsrf_lh(ig) = ptsrf(ig) + pdtsrf(ig) *ptimestep
                lh=(2834.3-0.28*(tsrf_lh(ig)-To)-
     &              0.004*(tsrf_lh(ig)-To)*(tsrf_lh(ig)-To))*1.e+3
                pdtsrf(ig)= pdtsrf(ig) + (dwatercap_dif(ig)+
     &               pdqsdif(ig,igcm_h2o_ice))*lh /pcapcal(ig)
              endif

               if(pqsurf(ig,igcm_h2o_ice)
     &           +(dwatercap_dif(ig)+pdqsdif(ig,igcm_h2o_ice))*ptimestep
     &           .gt.frost_albedo_threshold) ! if there is still ice, T cannot exceed To
     &           pdtsrf(ig) = min(pdtsrf(ig),(To-ptsrf(ig))/ptimestep) ! ice melt case
     
            ENDDO ! of DO ig=1,ngrid
          ELSE
c           Starting upward calculations for simple mixing of tracer (dust)
            DO ig=1,ngrid
               zq(ig,1,iq)=zc(ig,1)
            ENDDO
          END IF ! of IF ((water).and.(iq.eq.igcm_h2o_vap))

          DO ilay=2,nlay
             DO ig=1,ngrid
               zq(ig,ilay,iq)=zc(ig,ilay)+zd(ig,ilay)*zq(ig,ilay-1,iq)
             ENDDO
          ENDDO
        enddo ! of do iq=1,nq
      end if ! of if(tracer)
      

c-----------------------------------------------------------------------
c   8. calcul final des tendances de la diffusion verticale
c      ----------------------------------------------------

      DO ilev = 1, nlay
         DO ig=1,ngrid
            pdudif(ig,ilev)=(    zu(ig,ilev)-
     $      (pu(ig,ilev)+pdufi(ig,ilev)*ptimestep)    )/ptimestep
            pdvdif(ig,ilev)=(    zv(ig,ilev)-
     $      (pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep)    )/ptimestep
            hh = ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep 
     $  + (latcond*dmice(ig,ilev)/cpp)/ppopsk(ig,ilev)
            pdhdif(ig,ilev)=( zhs(ig,ilev)- hh )/ptimestep
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
     s      ph(ngrid/2+1,ilev),zhs(ngrid/2+1,ilev),
     s      pu(ngrid/2+1,ilev),zu(ngrid/2+1,ilev)

         ENDDO
      ENDIF

      END SUBROUTINE vdifc

      END MODULE vdifc_mod
