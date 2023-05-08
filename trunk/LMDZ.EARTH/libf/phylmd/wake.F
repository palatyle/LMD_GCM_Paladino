!
! $Id: wake.F 1403 2010-07-01 09:02:53Z fairhead $
!
      Subroutine WAKE (p,ph,pi,dtime,sigd_con
     :                ,te0,qe0,omgb
     :                ,dtdwn,dqdwn,amdwn,amup,dta,dqa
     :                ,wdtPBL,wdqPBL,udtPBL,udqPBL
     o                ,deltatw,deltaqw,dth,hw,sigmaw,wape,fip,gfl
     o                ,dtls,dqls
     o                ,ktopw,omgbdth,dp_omgb,wdens
     o                ,tu,qu
     o                ,dtKE,dqKE
     o                ,dtPBL,dqPBL
     o                ,omg,dp_deltomg,spread
     o                ,Cstar,d_deltat_gw
     o                ,d_deltatw2,d_deltaqw2)


***************************************************************
*                                                             *
* WAKE                                                        *
*      retour a un Pupper fixe                                *
*                                                             *
* written by   :  GRANDPEIX Jean-Yves   09/03/2000            *
* modified by :   ROEHRIG Romain        01/29/2007            *
***************************************************************
c
      use dimphy
      IMPLICIT none
c============================================================================
C
C
C   But : Decrire le comportement des poches froides apparaissant dans les
C        grands systemes convectifs, et fournir l'energie disponible pour
C        le declenchement de nouvelles colonnes convectives.
C
C   Variables d'etat : deltatw    : ecart de temperature wake-undisturbed area
C                      deltaqw    : ecart d'humidite wake-undisturbed area
C                      sigmaw     : fraction d'aire occupee par la poche.
C
C   Variable de sortie : 
c
c			 wape : WAke Potential Energy
c                        fip  : Front Incident Power (W/m2) - ALP
c                        gfl  : Gust Front Length per unit area (m-1)
C                        dtls : large scale temperature tendency due to wake
C                        dqls : large scale humidity tendency due to wake
C                        hw   : hauteur de la poche
C                     dp_omgb : vertical gradient of large scale omega
C                     wdens   : densite de poches
C                      omgbdth: flux of Delta_Theta transported by LS omega
C                      dtKE   : differential heating (wake - unpertubed)
C                      dqKE   : differential moistening (wake - unpertubed)
C                      omg    : Delta_omg =vertical velocity diff. wake-undist. (Pa/s)
C                 dp_deltomg  : vertical gradient of omg (s-1)
C                     spread  : spreading term in dt_wake and dq_wake
C                 deltatw     : updated temperature difference (T_w-T_u).
C                 deltaqw     : updated humidity difference (q_w-q_u).
C                 sigmaw      : updated wake fractional area.
C                 d_deltat_gw : delta T tendency due to GW
c
C   Variables d'entree : 
c
c		         aire : aire de la maille
c			 te0  : temperature dans l'environnement  (K)
C                        qe0  : humidite dans l'environnement     (kg/kg)
C                        omgb : vitesse verticale moyenne sur la maille (Pa/s)
C                        dtdwn: source de chaleur due aux descentes (K/s)
C                        dqdwn: source d'humidite due aux descentes (kg/kg/s)
C			 dta  : source de chaleur due courants satures et detrain  (K/s)
C			 dqa  : source d'humidite due aux courants satures et detra (kg/kg/s)
C                        amdwn: flux de masse total des descentes, par unite de
C                                surface de la maille (kg/m2/s)
C                        amup : flux de masse total des ascendances, par unite de
C                                surface de la maille (kg/m2/s)
C                        p    : pressions aux milieux des couches (Pa)
C                        ph   : pressions aux interfaces (Pa)
C                        pi  : (p/p_0)**kapa (adim)
C                        dtime: increment temporel (s)
c
C   Variables internes :
c
c			 rhow : masse volumique de la poche froide
C                        rho  : environment density at P levels
C                        rhoh : environment density at Ph levels
C                        te   : environment temperature | may change within
C                        qe   : environment humidity    | sub-time-stepping
C                        the  : environment potential temperature
C                        thu  : potential temperature in undisturbed area
C                        tu   :  temperature  in undisturbed area
C                        qu   : humidity in undisturbed area
C                      dp_omgb: vertical gradient og LS omega
C                      omgbw  : wake average vertical omega
C                     dp_omgbw: vertical gradient of omgbw
C                      omgbdq : flux of Delta_q transported by LS omega
C                        dth  : potential temperature diff. wake-undist.
C                        th1  : first pot. temp. for vertical advection (=thu)
C                        th2  : second pot. temp. for vertical advection (=thw)
C                        q1   : first humidity for vertical advection
C                        q2   : second humidity for vertical advection
C                     d_deltatw   : terme de redistribution pour deltatw
C                     d_deltaqw   : terme de redistribution pour deltaqw
C                      deltatw0   : deltatw initial
C                      deltaqw0   : deltaqw initial
C                      hw0    : hw initial
C                      sigmaw0: sigmaw initial
C                      amflux : horizontal mass flux through wake boundary
C                      wdens_ref: initial number of wakes per unit area (3D) or per
C                               unit length (2D), at the beginning of each time step
C                      Tgw    : 1 sur la période de onde de gravité
c                      Cgw    : vitesse de propagation de onde de gravité
c                      LL     : distance entre 2 poches

c-------------------------------------------------------------------------
c          Déclaration de variables
c-------------------------------------------------------------------------

#include "dimensions.h"
#include "YOMCST.h"
#include "cvthermo.h"
#include "iniprint.h"

c Arguments en entree
c--------------------

      REAL, dimension(klon,klev) :: p, pi
      REAL, dimension(klon,klev+1) :: ph, omgb
      REAL dtime
      REAL, dimension(klon,klev) :: te0,qe0
      REAL, dimension(klon,klev) :: dtdwn, dqdwn
      REAL, dimension(klon,klev) :: wdtPBL,wdqPBL
      REAL, dimension(klon,klev) :: udtPBL,udqPBL
      REAL, dimension(klon,klev) :: amdwn, amup
      REAL, dimension(klon,klev) :: dta, dqa
      REAL, dimension(klon) :: sigd_con

c Sorties
c--------

      REAL, dimension(klon,klev) :: deltatw, deltaqw, dth
      REAL, dimension(klon,klev) :: tu, qu
      REAL, dimension(klon,klev) :: dtls, dqls
      REAL, dimension(klon,klev) :: dtKE, dqKE
      REAL, dimension(klon,klev) :: dtPBL, dqPBL
      REAL, dimension(klon,klev) :: spread
      REAL, dimension(klon,klev) :: d_deltatgw
      REAL, dimension(klon,klev) :: d_deltatw2, d_deltaqw2
      REAL, dimension(klon,klev+1) :: omgbdth, omg
      REAL, dimension(klon,klev) :: dp_omgb, dp_deltomg
      REAL, dimension(klon,klev) :: d_deltat_gw
      REAL, dimension(klon) :: hw, sigmaw, wape, fip, gfl, Cstar
      REAL, dimension(klon) :: wdens
      INTEGER, dimension(klon) :: ktopw

c Variables internes
c-------------------

c Variables à fixer
      REAL ALON
      REAL coefgw
      REAL :: wdens_ref
      REAL stark
      REAL alpk
      REAL delta_t_min
      INTEGER nsub
      REAL dtimesub
      REAL sigmad, hwmin,wapecut
      REAL :: sigmaw_max
      REAL :: dens_rate
      REAL wdens0
cIM 080208
      LOGICAL, dimension(klon) :: gwake

c Variables de sauvegarde
      REAL, dimension(klon,klev) :: deltatw0
      REAL, dimension(klon,klev) :: deltaqw0
      REAL, dimension(klon,klev) :: te, qe
      REAL, dimension(klon) :: sigmaw0, sigmaw1

c Variables pour les GW
      REAL, DIMENSION(klon) :: LL
      REAL, dimension(klon,klev) :: N2
      REAL, dimension(klon,klev) :: Cgw
      REAL, dimension(klon,klev) :: Tgw

c Variables liées au calcul de hw
      REAL, DIMENSION(klon) :: ptop_provis, ptop, ptop_new
      REAL, DIMENSION(klon) :: sum_dth
      REAL, DIMENSION(klon) :: dthmin
      REAL, DIMENSION(klon) :: z, dz, hw0
      INTEGER, DIMENSION(klon) :: ktop, kupper

c Sub-timestep tendencies and related variables
       REAL d_deltatw(klon,klev),d_deltaqw(klon,klev)
       REAL d_te(klon,klev),d_qe(klon,klev)
       REAL d_sigmaw(klon),alpha(klon)
       REAL q0_min(klon),q1_min(klon)
       LOGICAL wk_adv(klon), OK_qx_qw(klon)
       REAL epsilon
       DATA epsilon/1.e-15/

c Autres variables internes
      INTEGER isubstep, k, i

      REAL, DIMENSION(klon) :: sum_thu, sum_tu, sum_qu,sum_thvu
      REAL, DIMENSION(klon) :: sum_dq, sum_rho
      REAL, DIMENSION(klon) :: sum_dtdwn, sum_dqdwn
      REAL, DIMENSION(klon) :: av_thu, av_tu, av_qu, av_thvu
      REAL, DIMENSION(klon) :: av_dth, av_dq, av_rho
      REAL, DIMENSION(klon) :: av_dtdwn, av_dqdwn

      REAL, DIMENSION(klon,klev) :: rho, rhow
      REAL, DIMENSION(klon,klev+1) :: rhoh
      REAL, DIMENSION(klon,klev) :: rhow_moyen
      REAL, DIMENSION(klon,klev) :: zh
      REAL, DIMENSION(klon,klev+1) :: zhh
      REAL, DIMENSION(klon,klev) :: epaisseur1, epaisseur2

      REAL, DIMENSION(klon,klev) :: the, thu

!      REAL, DIMENSION(klon,klev) :: d_deltatw, d_deltaqw

      REAL, DIMENSION(klon,klev+1) :: omgbw
      REAL, DIMENSION(klon) :: pupper
      REAL, DIMENSION(klon) :: omgtop
      REAL, DIMENSION(klon,klev) :: dp_omgbw
      REAL, DIMENSION(klon) :: ztop, dztop
      REAL, DIMENSION(klon,klev) :: alpha_up
      
      REAL, dimension(klon) :: RRe1, RRe2
      REAL :: RRd1, RRd2
      REAL, DIMENSION(klon,klev) :: Th1, Th2, q1, q2
      REAL, DIMENSION(klon,klev) :: D_Th1, D_Th2, D_dth
      REAL, DIMENSION(klon,klev) :: D_q1, D_q2, D_dq
      REAL, DIMENSION(klon,klev) :: omgbdq

      REAL, dimension(klon) :: ff, gg
      REAL, dimension(klon) :: wape2, Cstar2, heff

      REAL, DIMENSION(klon,klev) :: Crep
      REAL Crep_upper, Crep_sol

      REAL, DIMENSION(klon,klev) :: ppi

ccc nrlmd
      real, dimension(klon) :: death_rate,nat_rate
      real, dimension(klon,klev) :: entr
      real, dimension(klon,klev) :: detr

C-------------------------------------------------------------------------
c         Initialisations
c-------------------------------------------------------------------------

c      print*, 'wake initialisations'

c   Essais d'initialisation avec sigmaw = 0.02 et hw = 10.
c-------------------------------------------------------------------------

      DATA wapecut,sigmad, hwmin /5.,.02,10./
ccc nrlmd
      DATA sigmaw_max /0.4/
      DATA dens_rate /0.1/
ccc
C Longueur de maille (en m)
c-------------------------------------------------------------------------

c      ALON = 3.e5
      ALON = 1.e6


C Configuration de coefgw,stark,wdens (22/02/06 by YU Jingmei)
c
c      coefgw : Coefficient pour les ondes de gravité
c       stark : Coefficient k dans Cstar=k*sqrt(2*WAPE)
c       wdens : Densité de poche froide par maille
c-------------------------------------------------------------------------

ccc nrlmd      coefgw=10
c      coefgw=1
c      wdens0 = 1.0/(alon**2)
ccc nrlmd      wdens = 1.0/(alon**2)
ccc nrlmd      stark = 0.50
cCRtest
ccc nrlmd      alpk=0.1
c      alpk = 1.0
c      alpk = 0.5
c      alpk = 0.05
c
       stark  = 0.33
       Alpk   = 0.25
       wdens_ref  = 8.e-12
       coefgw = 4.
      Crep_upper=0.9
      Crep_sol=1.0

ccc nrlmd Lecture du fichier wake_param.data
      OPEN(99,file='wake_param.data',status='old',
     $          form='formatted',err=9999)
      READ(99,*,end=9998) stark
      READ(99,*,end=9998) Alpk
      READ(99,*,end=9998) wdens_ref
      READ(99,*,end=9998) coefgw
9998  Continue
      CLOSE(99)
9999  Continue
c
c   Initialisation de toutes des densites a wdens_ref.
c   Les densites peuvent evoluer si les poches debordent
c   (voir au tout debut de la boucle sur les substeps)
      wdens = wdens_ref
c
c      print*,'stark',stark
c      print*,'alpk',alpk
c      print*,'wdens',wdens
c      print*,'coefgw',coefgw
ccc
C Minimum value for |T_wake - T_undist|. Used for wake top definition
c-------------------------------------------------------------------------

      delta_t_min = 0.2

C 1. - Save initial values and initialize tendencies
C --------------------------------------------------

      DO k=1,klev
      DO i=1, klon
        ppi(i,k)=pi(i,k)
	deltatw0(i,k) = deltatw(i,k)
	deltaqw0(i,k)= deltaqw(i,k)
	te(i,k) = te0(i,k)
	qe(i,k) = qe0(i,k)
	dtls(i,k) = 0.
	dqls(i,k) = 0.
        d_deltat_gw(i,k)=0.
        d_te(i,k) = 0.
        d_qe(i,k) = 0.
        d_deltatw(i,k) = 0.
        d_deltaqw(i,k) = 0.
!IM 060508 beg
        d_deltatw2(i,k)=0.
        d_deltaqw2(i,k)=0.
!IM 060508 end
      ENDDO
      ENDDO
c      sigmaw1=sigmaw
c      IF (sigd_con.GT.sigmaw1) THEN
c      print*, 'sigmaw,sigd_con', sigmaw, sigd_con
c      ENDIF
      DO i=1, klon
cc      sigmaw(i) = amax1(sigmaw(i),sigd_con(i))
      sigmaw(i) = amax1(sigmaw(i),sigmad)
      sigmaw(i) = amin1(sigmaw(i),0.99)
      sigmaw0(i) = sigmaw(i)
      wape(i) = 0.
      wape2(i) = 0.
      d_sigmaw(i) = 0.
      ktopw(i) = 0
      ENDDO
C
C
C 2. - Prognostic part
C --------------------
C
C
C 2.1 - Undisturbed area and Wake integrals
C ---------------------------------------------------------

      DO i=1, klon
      z(i) = 0.
      ktop(i)=0
      kupper(i) = 0
      sum_thu(i) = 0.
      sum_tu(i) = 0.
      sum_qu(i) = 0.
      sum_thvu(i) = 0.
      sum_dth(i) = 0.
      sum_dq(i) = 0.
      sum_rho(i) = 0.
      sum_dtdwn(i) = 0.
      sum_dqdwn(i) = 0.

      av_thu(i) = 0.
      av_tu(i) =0.
      av_qu(i) =0.
      av_thvu(i) = 0.
      av_dth(i) = 0.
      av_dq(i) = 0.
      av_rho(i) =0.
      av_dtdwn(i) =0.
      av_dqdwn(i) = 0.
      ENDDO
c
c Distance between wakes
       DO i = 1,klon
        LL(i) = (1-sqrt(sigmaw(i)))/sqrt(wdens(i))
       ENDDO
C Potential temperatures and humidity
c----------------------------------------------------------
      DO k =1,klev
       DO i=1, klon 
!        write(*,*)'wake 1',i,k,rd,te(i,k)
        rho(i,k) = p(i,k)/(rd*te(i,k))
!        write(*,*)'wake 2',rho(i,k)
        IF(k .eq. 1) THEN
!        write(*,*)'wake 3',i,k,rd,te(i,k)
          rhoh(i,k) = ph(i,k)/(rd*te(i,k))
!        write(*,*)'wake 4',i,k,rd,te(i,k)
          zhh(i,k)=0
        ELSE
!          write(*,*)'wake 5',rd,(te(i,k)+te(i,k-1))
          rhoh(i,k) = ph(i,k)*2./(rd*(te(i,k)+te(i,k-1)))
!          write(*,*)'wake 6',(-rhoh(i,k)*RG)+zhh(i,k-1)
          zhh(i,k)=(ph(i,k)-ph(i,k-1))/(-rhoh(i,k)*RG)+zhh(i,k-1)
        ENDIF
!          write(*,*)'wake 7',ppi(i,k)
        the(i,k) = te(i,k)/ppi(i,k)
        thu(i,k) = (te(i,k) - deltatw(i,k)*sigmaw(i))/ppi(i,k)
        tu(i,k) = te(i,k) - deltatw(i,k)*sigmaw(i)
        qu(i,k)  =  qe(i,k) - deltaqw(i,k)*sigmaw(i)
!          write(*,*)'wake 8',(rd*(te(i,k)+deltatw(i,k)))
        rhow(i,k) = p(i,k)/(rd*(te(i,k)+deltatw(i,k)))
        dth(i,k) = deltatw(i,k)/ppi(i,k)
       ENDDO
      ENDDO
        
      DO k = 1, klev-1
      DO i=1, klon 
        IF(k.eq.1) THEN
          N2(i,k)=0
        ELSE
          N2(i,k)=amax1(0.,-RG**2/the(i,k)*rho(i,k)*(the(i,k+1)-
     $            the(i,k-1))/(p(i,k+1)-p(i,k-1)))
        ENDIF
        ZH(i,k)=(zhh(i,k)+zhh(i,k+1))/2

        Cgw(i,k)=sqrt(N2(i,k))*ZH(i,k)
        Tgw(i,k)=coefgw*Cgw(i,k)/LL(i)
      ENDDO
      ENDDO

      DO i=1, klon
      N2(i,klev)=0
      ZH(i,klev)=0
      Cgw(i,klev)=0
      Tgw(i,klev)=0
      ENDDO

c  Calcul de la masse volumique moyenne de la colonne   (bdlmd)
c-----------------------------------------------------------------

      DO k=1,klev
       DO i=1, klon
        epaisseur1(i,k)=0.
        epaisseur2(i,k)=0.
       ENDDO
      ENDDO

      DO i=1, klon
      epaisseur1(i,1)= -(ph(i,2)-ph(i,1))/(rho(i,1)*rg)+1.
      epaisseur2(i,1)= -(ph(i,2)-ph(i,1))/(rho(i,1)*rg)+1.
      rhow_moyen(i,1) = rhow(i,1)
      ENDDO

      DO k = 2, klev
      DO i=1, klon
        epaisseur1(i,k)= -(ph(i,k+1)-ph(i,k))/(rho(i,k)*rg) +1.
        epaisseur2(i,k)=epaisseur2(i,k-1)+epaisseur1(i,k)
        rhow_moyen(i,k) = (rhow_moyen(i,k-1)*epaisseur2(i,k-1)+
     $                 rhow(i,k)*epaisseur1(i,k))/epaisseur2(i,k)
      ENDDO
      ENDDO

C
C Choose an integration bound well above wake top
c-----------------------------------------------------------------
c
C       Pupper = 50000.  ! melting level
c       Pupper = 60000.
c       Pupper = 80000.  ! essais pour case_e
       DO i = 1,klon
        Pupper(i) = 0.6*ph(i,1)
        Pupper(i) = max(Pupper(i), 45000.)
ccc        Pupper(i) = 60000.
       ENDDO

C
C    Determine Wake top pressure (Ptop) from buoyancy integral
C    --------------------------------------------------------
c
c-1/ Pressure of the level where dth becomes less than delta_t_min.

      DO i=1,klon
      ptop_provis(i)=ph(i,1)
      ENDDO
      DO k= 2,klev
      DO i=1,klon
c
cIM v3JYG; ptop_provis(i).LT. ph(i,1)
c
        IF (dth(i,k) .GT. -delta_t_min .and.
     $      dth(i,k-1).LT. -delta_t_min .and.
     $      ptop_provis(i).EQ. ph(i,1)) THEN
          ptop_provis(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)
     $          - (dth(i,k-1)+delta_t_min)*p(i,k)) /
     $          (dth(i,k) - dth(i,k-1))
        ENDIF
      ENDDO
      ENDDO

c-2/ dth integral

      DO i=1,klon
      sum_dth(i) = 0.
      dthmin(i) = -delta_t_min
      z(i) = 0.
      ENDDO

      DO k = 1,klev
      DO i=1,klon
        dz(i) = -(amax1(ph(i,k+1),ptop_provis(i))-Ph(i,k))/(rho(i,k)*rg)
        IF (dz(i) .gt. 0) THEN
          z(i) = z(i)+dz(i)
          sum_dth(i) = sum_dth(i) + dth(i,k)*dz(i)
          dthmin(i) = amin1(dthmin(i),dth(i,k))
        ENDIF
      ENDDO
      ENDDO

c-3/ height of triangle with area= sum_dth and base = dthmin

      DO i=1,klon
      hw0(i) = 2.*sum_dth(i)/amin1(dthmin(i),-0.5)
      hw0(i) = amax1(hwmin,hw0(i))
      ENDDO

c-4/ now, get Ptop

      DO i=1,klon
      z(i) = 0.
      ptop(i) = ph(i,1)
      ENDDO

      DO k = 1,klev
      DO i=1,klon
        dz(i) = amin1(-(ph(i,k+1)-ph(i,k))/(rho(i,k)*rg),hw0(i)-z(i))
        IF (dz(i) .gt. 0) THEN
         z(i) = z(i)+dz(i)
         ptop(i) = ph(i,k)-rho(i,k)*rg*dz(i)
        ENDIF
      ENDDO
      ENDDO


C-5/ Determination de ktop et kupper

      DO k=klev,1,-1
      DO i=1,klon
        IF (ph(i,k+1) .lt. ptop(i)) ktop(i)=k
        IF (ph(i,k+1) .lt. pupper(i)) kupper(i)=k
      ENDDO
      ENDDO

c      On evite kupper = 1
      DO i=1,klon
        kupper(i) = max(kupper(i),2)
      ENDDO


c-6/ Correct ktop and ptop

      DO i = 1,klon
        ptop_new(i)=ptop(i)
      ENDDO
      DO k= klev,2,-1
      DO i=1,klon
        IF (k .LE. ktop(i) .and.
     $      ptop_new(i) .EQ. ptop(i) .and.
     $      dth(i,k) .GT. -delta_t_min .and.
     $      dth(i,k-1).LT. -delta_t_min) THEN
          ptop_new(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)
     $          - (dth(i,k-1)+delta_t_min)*p(i,k)) /
     $          (dth(i,k) - dth(i,k-1))
        ENDIF
      ENDDO
      ENDDO

      DO i=1,klon
        ptop(i) = ptop_new(i)
      ENDDO

      DO k=klev,1,-1
      DO i=1,klon
        IF (ph(i,k+1) .lt. ptop(i)) ktop(i)=k
      ENDDO
      ENDDO
c
c-5/ Set deltatw & deltaqw to 0 above kupper
c
      DO k = 1,klev
      DO i=1,klon
       IF (k.GE. kupper(i)) THEN
        deltatw(i,k) = 0.
        deltaqw(i,k) = 0.
       ENDIF
      ENDDO
      ENDDO
c
C
C Vertical gradient of LS omega
C
      DO k = 1,klev
      DO i=1,klon
       IF (k.LE. kupper(i)) THEN
        dp_omgb(i,k) = (omgb(i,k+1) - omgb(i,k))/(ph(i,k+1)-ph(i,k))
       ENDIF
      ENDDO
      ENDDO
C
C Integrals (and wake top level number)
C --------------------------------------
C
C Initialize sum_thvu to 1st level virt. pot. temp.

      DO i=1,klon
      z(i) = 1.
      dz(i) = 1.
      sum_thvu(i) =  thu(i,1)*(1.+eps*qu(i,1))*dz(i)
      sum_dth(i) = 0.
      ENDDO

      DO k = 1,klev
      DO i=1,klon
        dz(i) = -(amax1(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
        IF (dz(i) .GT. 0) THEN
         z(i) = z(i)+dz(i)
         sum_thu(i) = sum_thu(i) + thu(i,k)*dz(i)
         sum_tu(i) = sum_tu(i) + tu(i,k)*dz(i)
         sum_qu(i) = sum_qu(i) + qu(i,k)*dz(i)
         sum_thvu(i) = sum_thvu(i) + thu(i,k)*(1.+eps*qu(i,k))*dz(i)
         sum_dth(i) = sum_dth(i) + dth(i,k)*dz(i)
         sum_dq(i) = sum_dq(i) + deltaqw(i,k)*dz(i)
         sum_rho(i) = sum_rho(i) + rhow(i,k)*dz(i)
         sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i,k)*dz(i)
         sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i,k)*dz(i)
        ENDIF
      ENDDO
      ENDDO
c
      DO i=1,klon
        hw0(i) = z(i)
      ENDDO
c
C
C 2.1 - WAPE and mean forcing computation
C ---------------------------------------
C
C ---------------------------------------
C
C Means

      DO i=1,klon
      av_thu(i) = sum_thu(i)/hw0(i)
      av_tu(i) = sum_tu(i)/hw0(i)
      av_qu(i) = sum_qu(i)/hw0(i)
      av_thvu(i) = sum_thvu(i)/hw0(i)
c      av_thve = sum_thve/hw0
      av_dth(i) = sum_dth(i)/hw0(i)
      av_dq(i) = sum_dq(i)/hw0(i)
      av_rho(i) = sum_rho(i)/hw0(i)
      av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
      av_dqdwn(i) = sum_dqdwn(i)/hw0(i)

      wape(i) = - rg*hw0(i)*(av_dth(i)
     $     + eps*(av_thu(i)*av_dq(i)+av_dth(i)*av_qu(i)+av_dth(i)*
     $     av_dq(i) ))/av_thvu(i)
      ENDDO
C
C 2.2 Prognostic variable update
C ------------------------------
C
C Filter out bad wakes

      DO k = 1,klev
       DO i=1,klon
        IF ( wape(i) .LT. 0.) THEN
          deltatw(i,k) = 0.
          deltaqw(i,k) = 0.
          dth(i,k) = 0.
        ENDIF
       ENDDO
      ENDDO
c
      DO i=1,klon
      IF ( wape(i) .LT. 0.) THEN
        wape(i) = 0.
        Cstar(i) = 0.
        hw(i) = hwmin
        sigmaw(i) = amax1(sigmad,sigd_con(i))
        fip(i) = 0.
        gwake(i) = .FALSE.
      ELSE
        Cstar(i) = stark*sqrt(2.*wape(i))
        gwake(i) = .TRUE.
      ENDIF
      ENDDO

c
c Check qx and qw positivity
c --------------------------
      DO i = 1,klon
       q0_min(i)=min(  (qe(i,1)-sigmaw(i)*deltaqw(i,1)),
     $              (qe(i,1)+(1.-sigmaw(i))*deltaqw(i,1))  )
      ENDDO
      DO k = 2,klev
      DO i = 1,klon
        q1_min(i)=min(  (qe(i,k)-sigmaw(i)*deltaqw(i,k)),
     $              (qe(i,k)+(1.-sigmaw(i))*deltaqw(i,k))  )
        IF (q1_min(i).le.q0_min(i)) THEN
          q0_min(i)=q1_min(i)
        ENDIF
      ENDDO
      ENDDO
c
      DO i = 1,klon
       OK_qx_qw(i) = q0_min(i) .GE. 0.
       alpha(i) = 1.
      ENDDO
c
CC -----------------------------------------------------------------
C    Sub-time-stepping
C    -----------------
C
      nsub=10
      dtimesub=dtime/nsub
c
c------------------------------------------------------------
      DO isubstep = 1,nsub
c------------------------------------------------------------
c
c wk_adv is the logical flag enabling wake evolution in the time advance loop
      DO i = 1,klon
       wk_adv(i) = OK_qx_qw(i) .AND. alpha(i) .GE. 1.
      ENDDO
c
ccc nrlmd   Ajout d'un recalcul de wdens dans le cas d'un entrainement négatif de ktop à kupper --------
ccc           On calcule pour cela une densité wdens0 pour laquelle on aurait un entrainement nul ---
      DO i=1,klon
cc       print *,' isubstep,wk_adv(i),cstar(i),wape(i) ',
cc     $           isubstep,wk_adv(i),cstar(i),wape(i)
        IF (wk_adv(i) .AND. cstar(i).GT.0.01) THEN
           omg(i,kupper(i)+1) = - Rg*amdwn(i,kupper(i)+1)/sigmaw(i)
     $                + Rg*amup(i,kupper(i)+1)/(1.-sigmaw(i))
           wdens0 = ( sigmaw(i) / (4.*3.14) ) *
     $     ( (1.-sigmaw(i)) * omg(i,kupper(i)+1) /
     $     ( (ph(i,1)-pupper(i)) * cstar(i) )  ) **(2)
         IF ( wdens(i) .LE. wdens0*1.1 ) THEN
            wdens(i) = wdens0
         ENDIF
cc	   print*,'omg(i,kupper(i)+1),wdens0,wdens(i),cstar(i)
cc     $     ,ph(i,1)-pupper(i)',
cc     $             omg(i,kupper(i)+1),wdens0,wdens(i),cstar(i)
cc     $     ,ph(i,1)-pupper(i)
        ENDIF
      ENDDO

ccc nrlmd

      DO i=1,klon
       IF (wk_adv(i)) THEN
        gfl(i) = 2.*sqrt(3.14*wdens(i)*sigmaw(i))
        sigmaw(i)=amin1(sigmaw(i),sigmaw_max)
       ENDIF
      ENDDO
      DO i=1,klon
        IF (wk_adv(i)) THEN
ccc nrlmd          Introduction du taux de mortalité des poches et test sur sigmaw_max=0.4
ccc         d_sigmaw(i) = gfl(i)*Cstar(i)*dtimesub
           IF (sigmaw(i).ge.sigmaw_max) THEN
           death_rate(i)=gfl(i)*Cstar(i)/sigmaw(i)
           ELSE
             death_rate(i)=0.
           END IF
        d_sigmaw(i) = gfl(i)*Cstar(i)*dtimesub
     $               - death_rate(i)*sigmaw(i)*dtimesub
c     $              - nat_rate(i)*sigmaw(i)*dtimesub
cc        print*, 'd_sigmaw(i),sigmaw(i),gfl(i),Cstar(i),wape(i),
cc     $  death_rate(i),ktop(i),kupper(i)',
cc     $	         d_sigmaw(i),sigmaw(i),gfl(i),Cstar(i),wape(i),
cc     $  death_rate(i),ktop(i),kupper(i)

c        sigmaw(i) =sigmaw(i) + gfl(i)*Cstar(i)*dtimesub
c        sigmaw(i) =min(sigmaw(i),0.99)     !!!!!!!!
c        wdens = wdens0/(10.*sigmaw)
c        sigmaw =max(sigmaw,sigd_con)
c        sigmaw =max(sigmaw,sigmad)
        ENDIF
      ENDDO
C
C
c calcul de la difference de vitesse verticale poche - zone non perturbee
cIM 060208 differences par rapport au code initial; init. a 0 dp_deltomg
cIM 060208 et omg sur les niveaux de 1 a klev+1, alors que avant l'on definit
cIM 060208 au niveau k=1..?
      DO k= 1,klev
      DO i = 1,klon
      if (wk_adv(i)) THEN !!! nrlmd
        dp_deltomg(i,k)=0.
      end if
      ENDDO
      ENDDO
      DO k= 1,klev+1
      DO i = 1,klon
      if (wk_adv(i)) THEN !!! nrlmd
        omg(i,k)=0.
      end if
      ENDDO
      ENDDO
c
      DO i=1,klon
        IF (wk_adv(i)) THEN
        z(i)= 0.
        omg(i,1) = 0.
        dp_deltomg(i,1) = -(gfl(i)*Cstar(i))/(sigmaw(i) * (1-sigmaw(i)))
        ENDIF
      ENDDO
c
      DO k= 2,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. ktop(i)) THEN
          dz(i) = -(ph(i,k)-ph(i,k-1))/(rho(i,k-1)*rg)
          z(i) = z(i)+dz(i)
          dp_deltomg(i,k)= dp_deltomg(i,1)
          omg(i,k)= dp_deltomg(i,1)*z(i)
       ENDIF
      ENDDO
      ENDDO
c
      DO i = 1,klon
        IF (wk_adv(i)) THEN
        dztop(i)=-(ptop(i)-ph(i,ktop(i)))/(rho(i,ktop(i))*rg)
        ztop(i) = z(i)+dztop(i)
        omgtop(i)=dp_deltomg(i,1)*ztop(i)
        ENDIF
      ENDDO
c
c        -----------------
c        From m/s to Pa/s
c        -----------------
c
       DO i=1,klon
        IF (wk_adv(i)) THEN
        omgtop(i) = -rho(i,ktop(i))*rg*omgtop(i)
        dp_deltomg(i,1) = omgtop(i)/(ptop(i)-ph(i,1))
        ENDIF
       ENDDO
c
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. ktop(i)) THEN
          omg(i,k) = - rho(i,k)*rg*omg(i,k)
          dp_deltomg(i,k) = dp_deltomg(i,1)
       ENDIF
      ENDDO
      ENDDO
c
c   raccordement lineaire de omg de ptop a pupper

      DO i=1,klon
      IF ( wk_adv(i) .AND. kupper(i) .GT. ktop(i)) THEN
        omg(i,kupper(i)+1) = - Rg*amdwn(i,kupper(i)+1)/sigmaw(i)
     $                + Rg*amup(i,kupper(i)+1)/(1.-sigmaw(i))
        dp_deltomg(i,kupper(i)) = (omgtop(i)-omg(i,kupper(i)+1))/
     $                     (ptop(i)-pupper(i))
      ENDIF
      ENDDO
c
cc      DO i=1,klon
cc        print*,'Pente entre 0 et kupper (référence)'
cc     $   	,omg(i,kupper(i)+1)/(pupper(i)-ph(i,1))
cc        print*,'Pente entre ktop et kupper'
cc     $  	,(omg(i,kupper(i)+1)-omgtop(i))/(pupper(i)-ptop(i))
cc      ENDDO
cc
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .GT. ktop(i) .AND. k .LE. kupper(i)) THEN
          dp_deltomg(i,k) = dp_deltomg(i,kupper(i))
          omg(i,k) = omgtop(i)+(ph(i,k)-ptop(i))*dp_deltomg(i,kupper(i))
       ENDIF
      ENDDO
      ENDDO
ccc nrlmd
cc      DO i=1,klon
cc      print*,'deltaw_ktop,deltaw_conv',omgtop(i),omg(i,kupper(i)+1)
cc      END DO
ccc
c
c
c--    Compute wake average vertical velocity omgbw
c
c
      DO k = 1,klev+1
      DO i=1,klon
        IF ( wk_adv(i)) THEN
        omgbw(i,k) = omgb(i,k)+(1.-sigmaw(i))*omg(i,k)
        ENDIF
      ENDDO
      ENDDO
c--    and its vertical gradient dp_omgbw
c
      DO k = 1,klev
      DO i=1,klon
        IF ( wk_adv(i)) THEN
        dp_omgbw(i,k) = (omgbw(i,k+1)-omgbw(i,k))/(ph(i,k+1)-ph(i,k))
        ENDIF
      ENDDO
      ENDDO
C
c--    Upstream coefficients for omgb velocity
c--    (alpha_up(k) is the coefficient of the value at level k)
c--    (1-alpha_up(k) is the coefficient of the value at level k-1)
      DO k = 1,klev
      DO i=1,klon
        IF ( wk_adv(i)) THEN
         alpha_up(i,k) = 0.
         IF (omgb(i,k) .GT. 0.) alpha_up(i,k) = 1.
        ENDIF
      ENDDO
      ENDDO

c  Matrix expressing [The,deltatw] from  [Th1,Th2]

      DO i=1,klon
        IF ( wk_adv(i)) THEN
         RRe1(i) = 1.-sigmaw(i)
         RRe2(i) = sigmaw(i)
        ENDIF
      ENDDO
      RRd1 = -1.
      RRd2 = 1.
c
c--    Get [Th1,Th2], dth and [q1,q2]
c
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)+1) THEN
        dth(i,k) = deltatw(i,k)/ppi(i,k)
        Th1(i,k) = the(i,k) - sigmaw(i)     *dth(i,k)   ! undisturbed area
        Th2(i,k) = the(i,k) + (1.-sigmaw(i))*dth(i,k)   ! wake
        q1(i,k) = qe(i,k) - sigmaw(i)     *deltaqw(i,k) ! undisturbed area
        q2(i,k) = qe(i,k) + (1.-sigmaw(i))*deltaqw(i,k) ! wake
       ENDIF
      ENDDO
      ENDDO

      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
       D_Th1(i,1) = 0.
       D_Th2(i,1) = 0.
       D_dth(i,1) = 0.
       D_q1(i,1) = 0.
       D_q2(i,1) = 0.
       D_dq(i,1) = 0.
       end if
      ENDDO

      DO k= 2,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)+1) THEN
        D_Th1(i,k) = Th1(i,k-1)-Th1(i,k)
        D_Th2(i,k) = Th2(i,k-1)-Th2(i,k)
        D_dth(i,k) = dth(i,k-1)-dth(i,k)
        D_q1(i,k) = q1(i,k-1)-q1(i,k)
        D_q2(i,k) = q2(i,k-1)-q2(i,k)
        D_dq(i,k) = deltaqw(i,k-1)-deltaqw(i,k)
       ENDIF
      ENDDO
      ENDDO

      DO i=1,klon
        IF( wk_adv(i)) THEN
         omgbdth(i,1) = 0.
         omgbdq(i,1) = 0.
        ENDIF
      ENDDO

      DO k= 2,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)+1) THEN  !   loop on interfaces
        omgbdth(i,k) = omgb(i,k)*(    dth(i,k-1) -     dth(i,k))
        omgbdq(i,k)  = omgb(i,k)*(deltaqw(i,k-1) - deltaqw(i,k))
       ENDIF
      ENDDO
      ENDDO
c
c-----------------------------------------------------------------
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)-1) THEN
c-----------------------------------------------------------------
c
c   Compute redistribution (advective) term
c
         d_deltatw(i,k) =
     $             dtimesub/(Ph(i,k)-Ph(i,k+1))*(
     $       RRd1*omg(i,k  )*sigmaw(i)     *D_Th1(i,k)
     $      -RRd2*omg(i,k+1)*(1.-sigmaw(i))*D_Th2(i,k+1)
     $      -(1.-alpha_up(i,k))*omgbdth(i,k) - alpha_up(i,k+1)*
     $      omgbdth(i,k+1))*ppi(i,k)
c         print*,'d_deltatw=',d_deltatw(i,k)
c
         d_deltaqw(i,k) =
     $             dtimesub/(Ph(i,k)-Ph(i,k+1))*(
     $       RRd1*omg(i,k  )*sigmaw(i)     *D_q1(i,k)
     $      -RRd2*omg(i,k+1)*(1.-sigmaw(i))*D_q2(i,k+1)
     $      -(1.-alpha_up(i,k))*omgbdq(i,k) - alpha_up(i,k+1)*
     $      omgbdq(i,k+1))
c         print*,'d_deltaqw=',d_deltaqw(i,k)
c
c   and increment large scale tendencies
c

c
C
CC -----------------------------------------------------------------
         d_te(i,k) =  dtimesub*(
     $        ( RRe1(i)*omg(i,k  )*sigmaw(i)     *D_Th1(i,k)
     $         -RRe2(i)*omg(i,k+1)*(1.-sigmaw(i))*D_Th2(i,k+1) )
     $               /(Ph(i,k)-Ph(i,k+1))
ccc nrlmd     $         -sigmaw(i)*(1.-sigmaw(i))*dth(i,k)*dp_deltomg(i,k)
     $         -sigmaw(i)*(1.-sigmaw(i))*dth(i,k)
     $         *(omg(i,k)-omg(i,k+1))/(Ph(i,k)-Ph(i,k+1))
ccc
     $                      )*ppi(i,k)
c
         d_qe(i,k) =  dtimesub*(
     $        ( RRe1(i)*omg(i,k  )*sigmaw(i)     *D_q1(i,k)
     $         -RRe2(i)*omg(i,k+1)*(1.-sigmaw(i))*D_q2(i,k+1) )
     $               /(Ph(i,k)-Ph(i,k+1))
ccc nrlmd     $         -sigmaw(i)*(1.-sigmaw(i))*deltaqw(i,k)*dp_deltomg(i,k)
     $         -sigmaw(i)*(1.-sigmaw(i))*deltaqw(i,k)
     $         *(omg(i,k)-omg(i,k+1))/(Ph(i,k)-Ph(i,k+1))
ccc
     $                      )
ccc nrlmd
       ELSE IF(wk_adv(i) .AND. k .EQ. kupper(i)) THEN
         d_te(i,k) =  dtimesub*(
     $       ( RRe1(i)*omg(i,k  )*sigmaw(i)     *D_Th1(i,k)
     $        /(Ph(i,k)-Ph(i,k+1)))
     $                       )*ppi(i,k)

         d_qe(i,k) =  dtimesub*(
     $       ( RRe1(i)*omg(i,k  )*sigmaw(i)     *D_q1(i,k)       
     $        /(Ph(i,k)-Ph(i,k+1)))
     $                       )

       ENDIF
ccc
      ENDDO
      ENDDO
c------------------------------------------------------------------
C
C   Increment state variables

      DO k= 1,klev
      DO i = 1,klon
ccc nrlmd       IF( wk_adv(i) .AND. k .LE. kupper(i)-1) THEN
        IF( wk_adv(i) .AND. k .LE. kupper(i)) THEN
ccc


c
c Coefficient de répartition

        Crep(i,k)=Crep_sol*(ph(i,kupper(i))-ph(i,k))/(ph(i,kupper(i))
     $          -ph(i,1))
        Crep(i,k)=Crep(i,k)+Crep_upper*(ph(i,1)-ph(i,k))/(p(i,1)-
     $          ph(i,kupper(i)))
        

c Reintroduce compensating subsidence term.

c        dtKE(k)=(dtdwn(k)*Crep(k))/sigmaw
c        dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k))
c     .                   /(1-sigmaw)
c        dqKE(k)=(dqdwn(k)*Crep(k))/sigmaw
c        dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k))
c     .                   /(1-sigmaw)
c
c        dtKE(k)=(dtdwn(k)*Crep(k)+(1-Crep(k))*dta(k))/sigmaw
c        dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k)*Crep(k))
c     .                   /(1-sigmaw)
c        dqKE(k)=(dqdwn(k)*Crep(k)+(1-Crep(k))*dqa(k))/sigmaw
c        dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k)*Crep(k))
c     .                   /(1-sigmaw)

        dtKE(i,k)=(dtdwn(i,k)/sigmaw(i) - dta(i,k)/(1.-sigmaw(i)))
        dqKE(i,k)=(dqdwn(i,k)/sigmaw(i) - dqa(i,k)/(1.-sigmaw(i)))
c        print*,'dtKE= ',dtKE(i,k),' dqKE= ',dqKE(i,k)
c
        dtPBL(i,k)=(wdtPBL(i,k)/sigmaw(i) - udtPBL(i,k)/(1.-sigmaw(i)))
        dqPBL(i,k)=(wdqPBL(i,k)/sigmaw(i) - udqPBL(i,k)/(1.-sigmaw(i)))
c        print*,'dtPBL= ',dtPBL(i,k),' dqPBL= ',dqPBL(i,k)
c
ccc nrlmd          Prise en compte du taux de mortalité
ccc               Définitions de entr, detr
        detr(i,k)=0.

        entr(i,k)=detr(i,k)+gfl(i)*cstar(i)+
     $          sigmaw(i)*(1.-sigmaw(i))*dp_deltomg(i,k)

        spread(i,k) = (entr(i,k)-detr(i,k))/sigmaw(i)
ccc        spread(i,k) = (1.-sigmaw(i))*dp_deltomg(i,k)+gfl(i)*Cstar(i)/
ccc     $  sigmaw(i)


c ajout d'un effet onde de gravité -Tgw(k)*deltatw(k) 03/02/06 YU Jingmei

!      write(lunout,*)'wake.F ',i,k, dtimesub,d_deltat_gw(i,k),
!     &  Tgw(i,k),deltatw(i,k)
        d_deltat_gw(i,k)=d_deltat_gw(i,k)-Tgw(i,k)*deltatw(i,k)*
     $  dtimesub
!      write(lunout,*)'wake.F ',i,k, dtimesub,d_deltatw(i,k)
        ff(i)=d_deltatw(i,k)/dtimesub

c Sans GW
c
c        deltatw(k)=deltatw(k)+dtimesub*(ff+dtKE(k)-spread(k)*deltatw(k))
c
c GW formule 1
c
c        deltatw(k) = deltatw(k)+dtimesub*
c     $         (ff+dtKE(k) - spread(k)*deltatw(k)-Tgw(k)*deltatw(k))
c
c GW formule 2

        IF (dtimesub*Tgw(i,k).lt.1.e-10) THEN
          d_deltatw(i,k) = dtimesub*
     $       (ff(i)+dtKE(i,k)+dtPBL(i,k)
ccc     $       -spread(i,k)*deltatw(i,k)
     $       - entr(i,k)*deltatw(i,k)/sigmaw(i)
     $       - (death_rate(i)*sigmaw(i)+detr(i,k))*deltatw(i,k)
     $       / (1.-sigmaw(i))
ccc
     $       -Tgw(i,k)*deltatw(i,k))
        ELSE
           d_deltatw(i,k) = 1/Tgw(i,k)*(1-exp(-dtimesub*
     $       Tgw(i,k)))*
     $       (ff(i)+dtKE(i,k)+dtPBL(i,k)
ccc     $       -spread(i,k)*deltatw(i,k)
     $       - entr(i,k)*deltatw(i,k)/sigmaw(i)
     $       - (death_rate(i)*sigmaw(i)+detr(i,k))*deltatw(i,k)
     $       / (1.-sigmaw(i))
ccc
     $       -Tgw(i,k)*deltatw(i,k))
        ENDIF

        dth(i,k) = deltatw(i,k)/ppi(i,k)

        gg(i)=d_deltaqw(i,k)/dtimesub

       d_deltaqw(i,k) = dtimesub*(gg(i)+ dqKE(i,k)+dqPBL(i,k)
ccc     $     -spread(i,k)*deltaqw(i,k))
     $        - entr(i,k)*deltaqw(i,k)/sigmaw(i)
     $        - (death_rate(i)*sigmaw(i)+detr(i,k))*deltaqw(i,k)
     $        /(1.-sigmaw(i)))
ccc

ccc nrlmd
ccc       d_deltatw2(i,k)=d_deltatw2(i,k)+d_deltatw(i,k)
ccc       d_deltaqw2(i,k)=d_deltaqw2(i,k)+d_deltaqw(i,k)
ccc
       ENDIF
      ENDDO
      ENDDO

C
C   Scale tendencies so that water vapour remains positive in w and x.
C
      call wake_vec_modulation(klon,klev,wk_adv,epsilon,qe,d_qe,deltaqw,
     $                d_deltaqw,sigmaw,d_sigmaw,alpha)
c
ccc nrlmd
cc      print*,'alpha'
cc      do i=1,klon
cc         print*,alpha(i)
cc      end do
ccc
      DO k = 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)) THEN
        d_te(i,k)=alpha(i)*d_te(i,k)
        d_qe(i,k)=alpha(i)*d_qe(i,k)
        d_deltatw(i,k)=alpha(i)*d_deltatw(i,k)
        d_deltaqw(i,k)=alpha(i)*d_deltaqw(i,k)
        d_deltat_gw(i,k)=alpha(i)*d_deltat_gw(i,k)
       ENDIF
      ENDDO
      ENDDO
      DO i = 1,klon
       IF( wk_adv(i)) THEN
        d_sigmaw(i)=alpha(i)*d_sigmaw(i)
       ENDIF
      ENDDO

C   Update large scale variables and wake variables
cIM 060208 manque DO i + remplace DO k=1,kupper(i)
cIM 060208     DO k = 1,kupper(i)
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)) THEN
        dtls(i,k)=dtls(i,k)+d_te(i,k)
        dqls(i,k)=dqls(i,k)+d_qe(i,k)
ccc nrlmd
        d_deltatw2(i,k)=d_deltatw2(i,k)+d_deltatw(i,k)
        d_deltaqw2(i,k)=d_deltaqw2(i,k)+d_deltaqw(i,k)
ccc
       ENDIF
      ENDDO
      ENDDO
      DO k= 1,klev
      DO i = 1,klon
       IF( wk_adv(i) .AND. k .LE. kupper(i)) THEN
        te(i,k) = te0(i,k) + dtls(i,k)
        qe(i,k) = qe0(i,k) + dqls(i,k)
        the(i,k) = te(i,k)/ppi(i,k)
        deltatw(i,k) = deltatw(i,k)+d_deltatw(i,k)
        deltaqw(i,k) = deltaqw(i,k)+d_deltaqw(i,k)
        dth(i,k) = deltatw(i,k)/ppi(i,k)
cc      print*,'k,qx,qw',k,qe(i,k)-sigmaw(i)*deltaqw(i,k)
cc     $        ,qe(i,k)+(1-sigmaw(i))*deltaqw(i,k)
       ENDIF
      ENDDO
      ENDDO
      DO i = 1,klon
       IF( wk_adv(i)) THEN
        sigmaw(i) = sigmaw(i)+d_sigmaw(i)
       ENDIF
      ENDDO
c
C
c     Determine Ptop from buoyancy integral
c     ---------------------------------------
c
c-     1/ Pressure of the level where dth changes sign.
c
      DO i=1,klon
       IF ( wk_adv(i)) THEN
        Ptop_provis(i)=ph(i,1)
       ENDIF
      ENDDO
c
      DO k= 2,klev
      DO i=1,klon
        IF ( wk_adv(i) .AND.
     $       Ptop_provis(i) .EQ. ph(i,1) .AND.
     $      dth(i,k) .GT. -delta_t_min .and.
     $      dth(i,k-1).LT. -delta_t_min) THEN
          Ptop_provis(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)
     $          - (dth(i,k-1)+delta_t_min)*p(i,k)) /(dth(i,k)
     $          - dth(i,k-1))
        ENDIF
      ENDDO
      ENDDO
c
c-     2/ dth integral
c
      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
       sum_dth(i) = 0.
       dthmin(i) = -delta_t_min
       z(i) = 0.
       end if
      ENDDO

      DO k = 1,klev
      DO i=1,klon
       IF ( wk_adv(i)) THEN
        dz(i) = -(amax1(ph(i,k+1),Ptop_provis(i))-Ph(i,k))/(rho(i,k)*rg)
        IF (dz(i) .gt. 0) THEN
         z(i) = z(i)+dz(i)
         sum_dth(i) = sum_dth(i) + dth(i,k)*dz(i)
         dthmin(i) = amin1(dthmin(i),dth(i,k))
        ENDIF
       ENDIF
      ENDDO
      ENDDO
c
c-     3/ height of triangle with area= sum_dth and base = dthmin

      DO i=1,klon
       IF ( wk_adv(i)) THEN
         hw(i) = 2.*sum_dth(i)/amin1(dthmin(i),-0.5)
         hw(i) = amax1(hwmin,hw(i))
       ENDIF
      ENDDO
c
c-     4/ now, get Ptop
c
      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
       ktop(i) = 0
       z(i)=0.
       end if
      ENDDO
c
      DO k = 1,klev
      DO i=1,klon
       IF ( wk_adv(i)) THEN
        dz(i) = amin1(-(ph(i,k+1)-Ph(i,k))/(rho(i,k)*rg),hw(i)-z(i))
        IF (dz(i) .gt. 0) THEN
         z(i) = z(i)+dz(i)
         Ptop(i) = Ph(i,k)-rho(i,k)*rg*dz(i)
         ktop(i) = k
        ENDIF
       ENDIF
      ENDDO
      ENDDO
c
c      4.5/Correct ktop and ptop
c
      DO i=1,klon
       IF ( wk_adv(i)) THEN
        Ptop_new(i)=ptop(i)
       ENDIF
      ENDDO
c
      DO k= klev,2,-1
      DO i=1,klon
cIM v3JYG; IF (k .GE. ktop(i)
       IF ( wk_adv(i) .AND.
     $      k .LE. ktop(i) .AND.
     $      ptop_new(i) .EQ. ptop(i) .AND.
     $      dth(i,k) .GT. -delta_t_min .and.
     $      dth(i,k-1).LT. -delta_t_min) THEN
          Ptop_new(i) = ((dth(i,k)+delta_t_min)*p(i,k-1)
     $          - (dth(i,k-1)+delta_t_min)*p(i,k)) /(dth(i,k)
     $          - dth(i,k-1))
        ENDIF
      ENDDO
      ENDDO
c
c
      DO i=1,klon
       IF ( wk_adv(i)) THEN
        ptop(i) = ptop_new(i)
       ENDIF
      ENDDO

      DO k=klev,1,-1
      DO i=1,klon
      if (wk_adv(i)) then !!! nrlmd
        IF (ph(i,k+1) .LT. ptop(i)) ktop(i)=k
      end if
      ENDDO
      ENDDO
c
c      5/ Set deltatw & deltaqw to 0 above kupper
c
      DO k = 1,klev
      DO i=1,klon
        IF ( wk_adv(i) .AND. k .GE. kupper(i)) THEN
         deltatw(i,k) = 0.
         deltaqw(i,k) = 0.
        ENDIF
      ENDDO
      ENDDO
c
C
c-------------Cstar computation---------------------------------
      DO i=1, klon
       if (wk_adv(i)) then !!! nrlmd
      sum_thu(i) = 0.
      sum_tu(i) = 0.
      sum_qu(i) = 0.
      sum_thvu(i) = 0.
      sum_dth(i) = 0.
      sum_dq(i) = 0.
      sum_rho(i) = 0.
      sum_dtdwn(i) = 0.
      sum_dqdwn(i) = 0.

      av_thu(i) = 0.
      av_tu(i) =0.
      av_qu(i) =0.
      av_thvu(i) = 0.
      av_dth(i) = 0.
      av_dq(i) = 0.
      av_rho(i) =0.
      av_dtdwn(i) =0.
      av_dqdwn(i) = 0.
       end if
      ENDDO
C
C Integrals (and wake top level number)
C --------------------------------------
C
C Initialize sum_thvu to 1st level virt. pot. temp.

      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
      z(i) = 1.
      dz(i) = 1.
      sum_thvu(i) =  thu(i,1)*(1.+eps*qu(i,1))*dz(i)
      sum_dth(i) = 0.
       end if
      ENDDO

      DO k = 1,klev
      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
        dz(i) = -(max(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
        IF (dz(i) .GT. 0) THEN
         z(i) = z(i)+dz(i)
         sum_thu(i) = sum_thu(i) + thu(i,k)*dz(i)
         sum_tu(i) = sum_tu(i) + tu(i,k)*dz(i)
         sum_qu(i) = sum_qu(i) + qu(i,k)*dz(i)
         sum_thvu(i) = sum_thvu(i) + thu(i,k)*(1.+eps*qu(i,k))*dz(i)
         sum_dth(i) = sum_dth(i) + dth(i,k)*dz(i)
         sum_dq(i) = sum_dq(i) + deltaqw(i,k)*dz(i)
         sum_rho(i) = sum_rho(i) + rhow(i,k)*dz(i)
         sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i,k)*dz(i)
         sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i,k)*dz(i)
        ENDIF
       end if
      ENDDO
      ENDDO
c
      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
        hw0(i) = z(i)
       end if
      ENDDO
c
C
C - WAPE and mean forcing computation
C ---------------------------------------
C
C ---------------------------------------
C
C Means

      DO i=1,klon
       if (wk_adv(i)) then !!! nrlmd
       av_thu(i) = sum_thu(i)/hw0(i)
       av_tu(i) = sum_tu(i)/hw0(i)
       av_qu(i) = sum_qu(i)/hw0(i)
       av_thvu(i) = sum_thvu(i)/hw0(i)
       av_dth(i) = sum_dth(i)/hw0(i)
       av_dq(i) = sum_dq(i)/hw0(i)
       av_rho(i) = sum_rho(i)/hw0(i)
       av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
       av_dqdwn(i) = sum_dqdwn(i)/hw0(i)
c
       wape(i) = - rg*hw0(i)*(av_dth(i)
     $     + eps*(av_thu(i)*av_dq(i)+av_dth(i)*av_qu(i)+av_dth(i)*
     $     av_dq(i) ))/av_thvu(i)
       end if
      ENDDO
C
C Filter out bad wakes

      DO k = 1,klev
       DO i=1,klon
        if (wk_adv(i)) then !!! nrlmd
        IF ( wape(i) .LT. 0.) THEN
          deltatw(i,k) = 0.
          deltaqw(i,k) = 0.
          dth(i,k) = 0.
        ENDIF
	end if
       ENDDO
      ENDDO
c
      DO i=1,klon
      if (wk_adv(i)) then !!! nrlmd
      IF ( wape(i) .LT. 0.) THEN
        wape(i) = 0.
        Cstar(i) = 0.
        hw(i) = hwmin
        sigmaw(i) = max(sigmad,sigd_con(i))
        fip(i) = 0.
        gwake(i) = .FALSE.
      ELSE
        Cstar(i) = stark*sqrt(2.*wape(i))
        gwake(i) = .TRUE.
      ENDIF
      end if
      ENDDO

       ENDDO      ! end sub-timestep loop
C
C -----------------------------------------------------------------
c   Get back to tendencies per second
c
      DO k = 1,klev
      DO i=1,klon

ccc nrlmd        IF ( wk_adv(i) .AND. k .LE. kupper(i)) THEN
        IF ( OK_qx_qw(i) .AND. k .LE. kupper(i)) THEN
ccc
         dtls(i,k) = dtls(i,k)/dtime
         dqls(i,k) = dqls(i,k)/dtime
         d_deltatw2(i,k)=d_deltatw2(i,k)/dtime
         d_deltaqw2(i,k)=d_deltaqw2(i,k)/dtime
         d_deltat_gw(i,k) = d_deltat_gw(i,k)/dtime
cc      print*,'k,dqls,omg,entr,detr',k,dqls(i,k),omg(i,k),entr(i,k)
cc     $         ,death_rate(i)*sigmaw(i)
        ENDIF
      ENDDO
      ENDDO

c
c----------------------------------------------------------
c   Determine wake final state; recompute wape, cstar, ktop;
c   filter out bad wakes.
c----------------------------------------------------------
c
C 2.1 - Undisturbed area and Wake integrals
C ---------------------------------------------------------

      DO i=1,klon
ccc nrlmd       if (wk_adv(i)) then !!! nrlmd
      if (OK_qx_qw(i)) then
ccc
        z(i) = 0.
        sum_thu(i) = 0.
        sum_tu(i) = 0.
        sum_qu(i) = 0.
        sum_thvu(i) = 0.
        sum_dth(i) = 0.
        sum_dq(i) = 0.
        sum_rho(i) = 0.
        sum_dtdwn(i) = 0.
        sum_dqdwn(i) = 0.

        av_thu(i) = 0.
        av_tu(i) =0.
        av_qu(i) =0.
        av_thvu(i) = 0.
        av_dth(i) = 0.
        av_dq(i) = 0.
        av_rho(i) =0.
        av_dtdwn(i) =0.
        av_dqdwn(i) = 0.
       end if	
      ENDDO
C Potential temperatures and humidity
c----------------------------------------------------------

      DO k =1,klev
      DO i=1,klon
ccc nrlmd       IF ( wk_adv(i)) THEN
       if (OK_qx_qw(i)) then
ccc
        rho(i,k) = p(i,k)/(rd*te(i,k))
        IF(k .eq. 1) THEN
          rhoh(i,k) = ph(i,k)/(rd*te(i,k))
          zhh(i,k)=0
        ELSE
          rhoh(i,k) = ph(i,k)*2./(rd*(te(i,k)+te(i,k-1)))
          zhh(i,k)=(ph(i,k)-ph(i,k-1))/(-rhoh(i,k)*RG)+zhh(i,k-1)
        ENDIF
        the(i,k) = te(i,k)/ppi(i,k)
        thu(i,k) = (te(i,k) - deltatw(i,k)*sigmaw(i))/ppi(i,k)
        tu(i,k) = te(i,k) - deltatw(i,k)*sigmaw(i)
        qu(i,k)  =  qe(i,k) - deltaqw(i,k)*sigmaw(i)
        rhow(i,k) = p(i,k)/(rd*(te(i,k)+deltatw(i,k)))
        dth(i,k) = deltatw(i,k)/ppi(i,k)
       ENDIF
      ENDDO
      ENDDO

C Integrals (and wake top level number)
C -----------------------------------------------------------

C Initialize sum_thvu to 1st level virt. pot. temp.

      DO i=1,klon
ccc nrlmd       IF ( wk_adv(i)) THEN
      if (OK_qx_qw(i)) then
ccc
        z(i) = 1.
        dz(i) = 1.
        sum_thvu(i) =  thu(i,1)*(1.+eps*qu(i,1))*dz(i)
        sum_dth(i) = 0.
      ENDIF
      ENDDO

      DO k = 1,klev
      DO i=1,klon
ccc nrlmd       IF ( wk_adv(i)) THEN
       if (OK_qx_qw(i)) then
ccc
        dz(i) = -(amax1(ph(i,k+1),ptop(i))-ph(i,k))/(rho(i,k)*rg)
        IF (dz(i) .GT. 0) THEN
         z(i) = z(i)+dz(i)
         sum_thu(i) = sum_thu(i) + thu(i,k)*dz(i)
         sum_tu(i) = sum_tu(i) + tu(i,k)*dz(i)
         sum_qu(i) = sum_qu(i) + qu(i,k)*dz(i)
         sum_thvu(i) = sum_thvu(i) + thu(i,k)*(1.+eps*qu(i,k))*dz(i)
         sum_dth(i) = sum_dth(i) + dth(i,k)*dz(i)
         sum_dq(i) = sum_dq(i) + deltaqw(i,k)*dz(i)
         sum_rho(i) = sum_rho(i) + rhow(i,k)*dz(i)
         sum_dtdwn(i) = sum_dtdwn(i) + dtdwn(i,k)*dz(i)
         sum_dqdwn(i) = sum_dqdwn(i) + dqdwn(i,k)*dz(i)
        ENDIF
       ENDIF
      ENDDO
      ENDDO
c
      DO i=1,klon
ccc nrlmd       IF ( wk_adv(i)) THEN
       if (OK_qx_qw(i)) then
ccc
        hw0(i) = z(i)
       ENDIF
      ENDDO
c
C - WAPE and mean forcing computation
C-------------------------------------------------------------

C Means

      DO i=1, klon
ccc nrlmd       IF ( wk_adv(i)) THEN
      if (OK_qx_qw(i)) then
ccc
        av_thu(i) = sum_thu(i)/hw0(i)
        av_tu(i) = sum_tu(i)/hw0(i)
        av_qu(i) = sum_qu(i)/hw0(i)
        av_thvu(i) = sum_thvu(i)/hw0(i)
        av_dth(i) = sum_dth(i)/hw0(i)
        av_dq(i) = sum_dq(i)/hw0(i)
        av_rho(i) = sum_rho(i)/hw0(i)
        av_dtdwn(i) = sum_dtdwn(i)/hw0(i)
        av_dqdwn(i) = sum_dqdwn(i)/hw0(i)

        wape2(i) = - rg*hw0(i)*(av_dth(i)
     $     + eps*(av_thu(i)*av_dq(i)+av_dth(i)*av_qu(i)
     $     + av_dth(i)*av_dq(i) ))/av_thvu(i)
       ENDIF
      ENDDO

C Prognostic variable update
C ------------------------------------------------------------

C Filter out bad wakes
c
      DO k = 1,klev
      DO i=1,klon
ccc nrlmd        IF ( wk_adv(i) .AND. wape2(i) .LT. 0.) THEN
      if (OK_qx_qw(i) .AND. wape2(i) .LT. 0.) then
ccc
          deltatw(i,k) = 0.
          deltaqw(i,k) = 0.
          dth(i,k) = 0.
        ENDIF
      ENDDO
      ENDDO
c

      DO i=1, klon
ccc nrlmd       IF ( wk_adv(i)) THEN
      if (OK_qx_qw(i)) then
ccc
       IF ( wape2(i) .LT. 0.) THEN
        wape2(i) = 0.
        Cstar2(i) = 0.
        hw(i) = hwmin
        sigmaw(i) = amax1(sigmad,sigd_con(i))
        fip(i) = 0.
        gwake(i) = .FALSE.
      ELSE
        if(prt_level.ge.10) print*,'wape2>0'
        Cstar2(i) = stark*sqrt(2.*wape2(i))
        gwake(i) = .TRUE.
      ENDIF
      ENDIF
      ENDDO
c
      DO i=1, klon
ccc nrlmd       IF ( wk_adv(i)) THEN
       if (OK_qx_qw(i)) then
ccc
        ktopw(i) = ktop(i)
       ENDIF
      ENDDO
c
      DO i=1, klon
ccc nrlmd       IF ( wk_adv(i)) THEN
       if (OK_qx_qw(i)) then
ccc
       IF (ktopw(i) .gt. 0 .and. gwake(i)) then

Cjyg1     Utilisation d'un h_efficace constant ( ~ feeding layer)
ccc       heff = 600.
C      Utilisation de la hauteur hw
cc       heff = 0.7*hw
       heff(i) = hw(i)

       FIP(i) = 0.5*rho(i,ktopw(i))*Cstar2(i)**3*heff(i)*2*
     $      sqrt(sigmaw(i)*wdens(i)*3.14)
       FIP(i) = alpk * FIP(i)
Cjyg2
       ELSE
         FIP(i) = 0.
       ENDIF
       ENDIF
      ENDDO
c
C   Limitation de sigmaw

ccc nrlmd
c       DO i=1,klon
c         IF (OK_qx_qw(i)) THEN
c	   IF (sigmaw(i).GE.sigmaw_max) sigmaw(i)=sigmaw_max
c	 ENDIF
c       ENDDO
ccc
      DO k = 1,klev
       DO i=1, klon

ccc nrlmd      On maintient désormais constant sigmaw en régime permanent
ccc      IF ((sigmaw(i).GT.sigmaw_max).or.
        IF     ( ((wape(i).ge.wape2(i)).and.(wape2(i).le.1.0)).or.
     $         (ktopw(i).le.2) .OR.
     $         .not. OK_qx_qw(i) ) THEN
ccc
          dtls(i,k) = 0.
          dqls(i,k) = 0.
          deltatw(i,k) = 0.
          deltaqw(i,k) = 0.
        ENDIF
       ENDDO
      ENDDO
c
ccc nrlmd      On maintient désormais constant sigmaw en régime permanent
      DO i=1, klon
        IF  ( ((wape(i).ge.wape2(i)).and.(wape2(i).le.1.0)).or.
     $        (ktopw(i).le.2) .OR.
     $        .not. OK_qx_qw(i)   ) THEN
         wape(i) = 0.
         cstar(i)=0.
         hw(i) = hwmin
         sigmaw(i) = sigmad
         fip(i) = 0.
        ELSE
         wape(i) = wape2(i)
         cstar(i)=cstar2(i)
        ENDIF
cc        print*,'wape wape2 ktopw OK_qx_qw =',
cc     $          wape(i),wape2(i),ktopw(i),OK_qx_qw(i)
      ENDDO
c
c
      RETURN
      END

      SUBROUTINE wake_vec_modulation(nlon,nl,wk_adv,epsilon,qe,d_qe,
     $           deltaqw,d_deltaqw,sigmaw,d_sigmaw,alpha)
c------------------------------------------------------
cDtermination du coefficient alpha tel que les tendances
c corriges alpha*d_G, pour toutes les grandeurs G, correspondent
c a une humidite positive dans la zone (x) et dans la zone (w).
c------------------------------------------------------
c
 
c  Input
      REAL qe(nlon,nl),d_qe(nlon,nl)
      REAL deltaqw(nlon,nl),d_deltaqw(nlon,nl)
      REAL sigmaw(nlon),d_sigmaw(nlon)
      LOGICAL wk_adv(nlon)
      INTEGER nl,nlon
c  Output
      REAL alpha(nlon)
c  Internal variables
      REAL zeta(nlon,nl)
      REAL alpha1(nlon)
      REAL x,a,b,c,discrim
      REAL epsilon
!      DATA epsilon/1.e-15/
c
      DO k=1,nl
      DO i = 1,nlon
       IF (wk_adv(i)) THEN
        IF ((deltaqw(i,k)+d_deltaqw(i,k)).ge.0.) then
         zeta(i,k)=0.
        ELSE
         zeta(i,k)=1.
        END IF
       ENDIF
      ENDDO
      DO i = 1,nlon
       IF (wk_adv(i)) THEN
        x = qe(i,k)+(zeta(i,k)-sigmaw(i))*deltaqw(i,k)
     $    + d_qe(i,k)+(zeta(i,k)-sigmaw(i))*d_deltaqw(i,k)
     $    - d_sigmaw(i)*(deltaqw(i,k)+d_deltaqw(i,k))
        a = -d_sigmaw(i)*d_deltaqw(i,k)
        b = d_qe(i,k)+(zeta(i,k)-sigmaw(i))*d_deltaqw(i,k)
     $    - deltaqw(i,k)*d_sigmaw(i)
        c = qe(i,k)+(zeta(i,k)-sigmaw(i))*deltaqw(i,k)+epsilon
        discrim = b*b-4.*a*c
c      print*, 'x, a, b, c, discrim', x, a, b, c, discrim
        IF (a+b .GE. 0.) THEN !! Condition suffisante pour la positivité de ovap
         alpha1(i)=1.
        ELSE
         IF (x .GE. 0.) THEN
            alpha1(i)=1.
         ELSE
              IF (a .GT. 0.) THEN
                 alpha1(i)=0.9*min(   (2.*c)/(-b+sqrt(discrim)),
     $                        (-b+sqrt(discrim))/(2.*a)   )
              ELSE IF (a .eq. 0.) then
                 alpha1(i)=0.9*(-c/b)
              ELSE
c         print*,'a,b,c discrim',a,b,c discrim
                 alpha1(i)=0.9*max(   (2.*c)/(-b+sqrt(discrim)),
     $                        (-b+sqrt(discrim))/(2.*a)   )
              ENDIF
         ENDIF
        ENDIF
       alpha(i) = min(alpha(i),alpha1(i))
       ENDIF
      ENDDO
      ENDDO
!
      return
      end

      Subroutine WAKE_scal (p,ph,ppi,dtime,sigd_con
     :                ,te0,qe0,omgb
     :                ,dtdwn,dqdwn,amdwn,amup,dta,dqa
     :                ,wdtPBL,wdqPBL,udtPBL,udqPBL
     o                ,deltatw,deltaqw,dth,hw,sigmaw,wape,fip,gfl
     o                ,dtls,dqls
     o                ,ktopw,omgbdth,dp_omgb,wdens
     o                ,tu,qu
     o                ,dtKE,dqKE
     o                ,dtPBL,dqPBL
     o                ,omg,dp_deltomg,spread
     o                ,Cstar,d_deltat_gw
     o                ,d_deltatw2,d_deltaqw2)

***************************************************************
*                                                             *
* WAKE                                                        *
*      retour a un Pupper fixe                                *
*                                                             *
* written by   :  GRANDPEIX Jean-Yves   09/03/2000            *
* modified by :   ROEHRIG Romain        01/29/2007            *
***************************************************************
c
      USE dimphy
      IMPLICIT none
c============================================================================
C
C
C   But : Decrire le comportement des poches froides apparaissant dans les
C        grands systemes convectifs, et fournir l'energie disponible pour
C        le declenchement de nouvelles colonnes convectives.
C
C   Variables d'etat : deltatw    : ecart de temperature wake-undisturbed area
C                      deltaqw    : ecart d'humidite wake-undisturbed area
C                      sigmaw     : fraction d'aire occupee par la poche.
C
C   Variable de sortie : 
c
c			 wape : WAke Potential Energy
c                        fip  : Front Incident Power (W/m2) - ALP
c                        gfl  : Gust Front Length per unit area (m-1)
C                        dtls : large scale temperature tendency due to wake
C                        dqls : large scale humidity tendency due to wake
C                        hw   : hauteur de la poche
C                     dp_omgb : vertical gradient of large scale omega
C                      omgbdth: flux of Delta_Theta transported by LS omega
C                      dtKE   : differential heating (wake - unpertubed)
C                      dqKE   : differential moistening (wake - unpertubed)
C                      omg    : Delta_omg =vertical velocity diff. wake-undist. (Pa/s)
C                 dp_deltomg  : vertical gradient of omg (s-1)
C                     spread  : spreading term in dt_wake and dq_wake
C                 deltatw     : updated temperature difference (T_w-T_u).
C                 deltaqw     : updated humidity difference (q_w-q_u).
C                 sigmaw      : updated wake fractional area.
C                 d_deltat_gw : delta T tendency due to GW
c
C   Variables d'entree : 
c
c		         aire : aire de la maille
c			 te0  : temperature dans l'environnement  (K)
C                        qe0  : humidite dans l'environnement     (kg/kg)
C                        omgb : vitesse verticale moyenne sur la maille (Pa/s)
C                        dtdwn: source de chaleur due aux descentes (K/s)
C                        dqdwn: source d'humidite due aux descentes (kg/kg/s)
C			 dta  : source de chaleur due courants satures et detrain  (K/s)
C			 dqa  : source d'humidite due aux courants satures et detra (kg/kg/s)
C                        amdwn: flux de masse total des descentes, par unite de
C                                surface de la maille (kg/m2/s)
C                        amup : flux de masse total des ascendances, par unite de
C                                surface de la maille (kg/m2/s)
C                        p    : pressions aux milieux des couches (Pa)
C                        ph   : pressions aux interfaces (Pa)
C                        ppi  : (p/p_0)**kapa (adim)
C                        dtime: increment temporel (s)
c
C   Variables internes :
c
c			 rhow : masse volumique de la poche froide
C                        rho  : environment density at P levels
C                        rhoh : environment density at Ph levels
C                        te   : environment temperature | may change within
C                        qe   : environment humidity    | sub-time-stepping
C                        the  : environment potential temperature
C                        thu  : potential temperature in undisturbed area
C                        tu   :  temperature  in undisturbed area
C                        qu   : humidity in undisturbed area
C                      dp_omgb: vertical gradient og LS omega
C                      omgbw  : wake average vertical omega
C                     dp_omgbw: vertical gradient of omgbw
C                      omgbdq : flux of Delta_q transported by LS omega
C                        dth  : potential temperature diff. wake-undist.
C                        th1  : first pot. temp. for vertical advection (=thu)
C                        th2  : second pot. temp. for vertical advection (=thw)
C                        q1   : first humidity for vertical advection
C                        q2   : second humidity for vertical advection
C                     d_deltatw   : terme de redistribution pour deltatw
C                     d_deltaqw   : terme de redistribution pour deltaqw
C                      deltatw0   : deltatw initial
C                      deltaqw0   : deltaqw initial
C                      hw0    : hw initial
C                      sigmaw0: sigmaw initial
C                      amflux : horizontal mass flux through wake boundary
C                      wdens  : number of wakes per unit area (3D) or per
C                               unit length (2D)
C                      Tgw    : 1 sur la période de onde de gravité
c                      Cgw    : vitesse de propagation de onde de gravité
c                      LL     : distance entre 2 poches

c-------------------------------------------------------------------------
c          Déclaration de variables
c-------------------------------------------------------------------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"
#include "cvthermo.h"
#include "iniprint.h"

c Arguments en entree
c--------------------

      REAL p(klev),ph(klev+1),ppi(klev)
      REAL dtime
      REAL te0(klev),qe0(klev)
      REAL omgb(klev+1)
      REAL dtdwn(klev), dqdwn(klev)
      REAL wdtPBL(klev),wdqPBL(klev)
      REAL udtPBL(klev),udqPBL(klev)
      REAL amdwn(klev), amup(klev)
      REAL dta(klev), dqa(klev)
      REAL sigd_con

c Sorties
c--------

      REAL deltatw(klev), deltaqw(klev), dth(klev)
      REAL tu(klev), qu(klev)
      REAL dtls(klev), dqls(klev)
      REAL dtKE(klev), dqKE(klev)
      REAL dtPBL(klev), dqPBL(klev)
      REAL spread(klev)
      REAL d_deltatgw(klev)
      REAL d_deltatw2(klev), d_deltaqw2(klev)
      REAL omgbdth(klev+1), omg(klev+1)
      REAL dp_omgb(klev), dp_deltomg(klev)
      REAL d_deltat_gw(klev)
      REAL hw, sigmaw, wape, fip, gfl, Cstar
      INTEGER ktopw

c Variables internes
c-------------------

c Variables à fixer
      REAL ALON
      REAL coefgw
      REAL wdens0, wdens
      REAL stark
      REAL alpk
      REAL delta_t_min
      REAL Pupper
      INTEGER nsub
      REAL dtimesub
      REAL sigmad, hwmin

c Variables de sauvegarde
      REAL deltatw0(klev)
      REAL deltaqw0(klev)
      REAL te(klev), qe(klev)
      REAL sigmaw0, sigmaw1

c Variables pour les GW
      REAL LL
      REAL N2(klev)
      REAL Cgw(klev)
      REAL Tgw(klev)

c Variables liées au calcul de hw
      REAL ptop_provis, ptop, ptop_new
      REAL sum_dth
      REAL dthmin
      REAL z, dz, hw0
      INTEGER ktop, kupper

c Autres variables internes
      INTEGER isubstep, k

      REAL sum_thu, sum_tu, sum_qu,sum_thvu
      REAL sum_dq, sum_rho
      REAL sum_dtdwn, sum_dqdwn
      REAL av_thu, av_tu, av_qu, av_thvu
      REAL av_dth, av_dq, av_rho
      REAL av_dtdwn, av_dqdwn

      REAL rho(klev), rhoh(klev+1), rhow(klev)
      REAL rhow_moyen(klev)
      REAL zh(klev), zhh(klev+1)
      REAL epaisseur1(klev), epaisseur2(klev)

      REAL the(klev), thu(klev)

      REAL d_deltatw(klev), d_deltaqw(klev)

      REAL omgbw(klev+1), omgtop
      REAL dp_omgbw(klev)
      REAL ztop, dztop
      REAL alpha_up(klev)
      
      REAL RRe1, RRe2, RRd1, RRd2
      REAL Th1(klev), Th2(klev), q1(klev), q2(klev)
      REAL D_Th1(klev), D_Th2(klev), D_dth(klev)
      REAL D_q1(klev), D_q2(klev), D_dq(klev)
      REAL omgbdq(klev)

      REAL ff, gg
      REAL wape2, Cstar2, heff

      REAL Crep(klev)
      REAL Crep_upper, Crep_sol

C-------------------------------------------------------------------------
c         Initialisations
c-------------------------------------------------------------------------

c      print*, 'wake initialisations'

c   Essais d'initialisation avec sigmaw = 0.02 et hw = 10.
c-------------------------------------------------------------------------

      DATA sigmad, hwmin /.02,10./

C Longueur de maille (en m)
c-------------------------------------------------------------------------

c      ALON = 3.e5
      ALON = 1.e6


C Configuration de coefgw,stark,wdens (22/02/06 by YU Jingmei)
c
c      coefgw : Coefficient pour les ondes de gravité
c       stark : Coefficient k dans Cstar=k*sqrt(2*WAPE)
c       wdens : Densité de poche froide par maille
c-------------------------------------------------------------------------

      coefgw=10
c      coefgw=1
c      wdens0 = 1.0/(alon**2)   
      wdens = 1.0/(alon**2)       
      stark = 0.50
cCRtest
      alpk=0.1
c      alpk = 1.0 
c      alpk = 0.5
c      alpk = 0.05
      Crep_upper=0.9
      Crep_sol=1.0


C Minimum value for |T_wake - T_undist|. Used for wake top definition
c-------------------------------------------------------------------------

      delta_t_min = 0.2


C 1. - Save initial values and initialize tendencies
C --------------------------------------------------

      DO k=1,klev
	deltatw0(k) = deltatw(k)
	deltaqw0(k)= deltaqw(k)
	te(k) = te0(k)
	qe(k) = qe0(k)
	dtls(k) = 0.
	dqls(k) = 0.
        d_deltat_gw(k)=0.
        d_deltatw2(k)=0.
        d_deltaqw2(k)=0.
      ENDDO
c      sigmaw1=sigmaw
c      IF (sigd_con.GT.sigmaw1) THEN
c      print*, 'sigmaw,sigd_con', sigmaw, sigd_con
c      ENDIF
      sigmaw = max(sigmaw,sigd_con)
      sigmaw = max(sigmaw,sigmad)
      sigmaw = min(sigmaw,0.99)
      sigmaw0 = sigmaw
c      wdens=wdens0/(10.*sigmaw)
c      IF (sigd_con.GT.sigmaw1) THEN
c      print*, 'sigmaw1,sigd1', sigmaw, sigd_con
c      ENDIF

C 2. - Prognostic part
C =========================================================

c      print *, 'prognostic wake computation'


C 2.1 - Undisturbed area and Wake integrals
C ---------------------------------------------------------

      z = 0.
      ktop=0
      kupper = 0
      sum_thu = 0.
      sum_tu = 0.
      sum_qu = 0.
      sum_thvu = 0.
      sum_dth = 0.
      sum_dq = 0.
      sum_rho = 0.
      sum_dtdwn = 0.
      sum_dqdwn = 0.

      av_thu = 0.
      av_tu =0.
      av_qu =0.
      av_thvu = 0.
      av_dth = 0.
      av_dq = 0.
      av_rho =0.
      av_dtdwn =0.
      av_dqdwn = 0.

C Potential temperatures and humidity
c----------------------------------------------------------

      DO k =1,klev
        rho(k) = p(k)/(rd*te(k))
        IF(k .eq. 1) THEN
          rhoh(k) = ph(k)/(rd*te(k))
          zhh(k)=0
        ELSE
          rhoh(k) = ph(k)*2./(rd*(te(k)+te(k-1)))
          zhh(k)=(ph(k)-ph(k-1))/(-rhoh(k)*RG)+zhh(k-1)
        ENDIF
        the(k) = te(k)/ppi(k)
        thu(k) = (te(k) - deltatw(k)*sigmaw)/ppi(k)
        tu(k) = te(k) - deltatw(k)*sigmaw
        qu(k)  =  qe(k) - deltaqw(k)*sigmaw
        rhow(k) = p(k)/(rd*(te(k)+deltatw(k)))
        dth(k) = deltatw(k)/ppi(k)
        LL = (1-sqrt(sigmaw))/sqrt(wdens)       
      ENDDO
        
      DO k = 1, klev-1
        IF(k.eq.1) THEN
          N2(k)=0
        ELSE
          N2(k)=max(0.,-RG**2/the(k)*rho(k)*(the(k+1)-the(k-1))
     $           /(p(k+1)-p(k-1)))
        ENDIF
        ZH(k)=(zhh(k)+zhh(k+1))/2

        Cgw(k)=sqrt(N2(k))*ZH(k)
        Tgw(k)=coefgw*Cgw(k)/LL
      ENDDO
         
      N2(klev)=0
      ZH(klev)=0
      Cgw(klev)=0
      Tgw(klev)=0

c  Calcul de la masse volumique moyenne de la colonne
c-----------------------------------------------------------------

      DO k=1,klev
        epaisseur1(k)=0.
        epaisseur2(k)=0.
      ENDDO

      epaisseur1(1)= -(Ph(2)-Ph(1))/(rho(1)*rg)+1.
      epaisseur2(1)= -(Ph(2)-Ph(1))/(rho(1)*rg)+1.
      rhow_moyen(1) = rhow(1)

      DO k = 2, klev
        epaisseur1(k)= -(Ph(k+1)-Ph(k))/(rho(k)*rg) +1.
        epaisseur2(k)=epaisseur2(k-1)+epaisseur1(k)
        rhow_moyen(k) = (rhow_moyen(k-1)*epaisseur2(k-1)+
     $                 rhow(k)*epaisseur1(k))/epaisseur2(k)
      ENDDO


C Choose an integration bound well above wake top
c-----------------------------------------------------------------

c       Pupper = 50000.  ! melting level
       Pupper = 60000.
c       Pupper = 70000.


C    Determine Wake top pressure (Ptop) from buoyancy integral
C-----------------------------------------------------------------

c-1/ Pressure of the level where dth becomes less than delta_t_min.

      Ptop_provis=ph(1)
      DO k= 2,klev
        IF (dth(k) .GT. -delta_t_min .and.
     $      dth(k-1).LT. -delta_t_min) THEN
          Ptop_provis = ((dth(k)+delta_t_min)*p(k-1)
     $          - (dth(k-1)+delta_t_min)*p(k)) /(dth(k) - dth(k-1))
          GO TO 25
        ENDIF
      ENDDO
25    CONTINUE

c-2/ dth integral

      sum_dth = 0.
      dthmin = -delta_t_min
      z = 0.

      DO k = 1,klev
        dz = -(max(ph(k+1),Ptop_provis)-Ph(k))/(rho(k)*rg)
        IF (dz .le. 0) GO TO 40
        z = z+dz
        sum_dth = sum_dth + dth(k)*dz
        dthmin = min(dthmin,dth(k))
      ENDDO
40    CONTINUE

c-3/ height of triangle with area= sum_dth and base = dthmin

      hw0 = 2.*sum_dth/min(dthmin,-0.5)
      hw0 = max(hwmin,hw0)

c-4/ now, get Ptop

      z = 0.
      ptop = ph(1)

      DO k = 1,klev
        dz = min(-(ph(k+1)-Ph(k))/(rho(k)*rg),hw0-z)
        IF (dz .le. 0) GO TO 45
        z = z+dz
        Ptop = Ph(k)-rho(k)*rg*dz
      ENDDO
45    CONTINUE


C-5/ Determination de ktop et kupper

      DO k=klev,1,-1
        IF (ph(k+1) .lt. ptop) ktop=k
        IF (ph(k+1) .lt. pupper) kupper=k
      ENDDO

c-6/ Correct ktop and ptop

      Ptop_new=ptop
      DO k= ktop,2,-1
        IF (dth(k) .GT. -delta_t_min .and.
     $      dth(k-1).LT. -delta_t_min) THEN
          Ptop_new = ((dth(k)+delta_t_min)*p(k-1)
     $          - (dth(k-1)+delta_t_min)*p(k)) /(dth(k) - dth(k-1))
          GO TO 225
        ENDIF
      ENDDO
225   CONTINUE

      ptop = ptop_new

      DO k=klev,1,-1
        IF (ph(k+1) .lt. ptop) ktop=k
      ENDDO

c Set deltatw & deltaqw to 0 above kupper
c-----------------------------------------------------------

      DO k = kupper,klev
        deltatw(k) = 0.
        deltaqw(k) = 0.
      ENDDO


C Vertical gradient of LS omega
C------------------------------------------------------------

      DO k = 1,kupper
        dp_omgb(k) = (omgb(k+1) - omgb(k))/(ph(k+1)-ph(k))
      ENDDO


C Integrals (and wake top level number)
C -----------------------------------------------------------

C Initialize sum_thvu to 1st level virt. pot. temp.

      z = 1.
      dz = 1.
      sum_thvu =  thu(1)*(1.+eps*qu(1))*dz
      sum_dth = 0.

      DO k = 1,klev
        dz = -(max(ph(k+1),Ptop)-Ph(k))/(rho(k)*rg)
        IF (dz .LE. 0) GO TO 50
        z = z+dz
        sum_thu = sum_thu + thu(k)*dz
        sum_tu = sum_tu + tu(k)*dz
        sum_qu = sum_qu + qu(k)*dz
        sum_thvu = sum_thvu + thu(k)*(1.+eps*qu(k))*dz
        sum_dth = sum_dth + dth(k)*dz
        sum_dq = sum_dq + deltaqw(k)*dz
        sum_rho = sum_rho + rhow(k)*dz
        sum_dtdwn = sum_dtdwn + dtdwn(k)*dz
        sum_dqdwn = sum_dqdwn + dqdwn(k)*dz
      ENDDO
50    CONTINUE

      hw0 = z

C 2.1 - WAPE and mean forcing computation
C-------------------------------------------------------------

C Means

      av_thu = sum_thu/hw0
      av_tu = sum_tu/hw0
      av_qu = sum_qu/hw0
      av_thvu = sum_thvu/hw0
c      av_thve = sum_thve/hw0
      av_dth = sum_dth/hw0
      av_dq = sum_dq/hw0
      av_rho = sum_rho/hw0
      av_dtdwn = sum_dtdwn/hw0
      av_dqdwn = sum_dqdwn/hw0

      wape = - rg*hw0*(av_dth
     $     + eps*(av_thu*av_dq+av_dth*av_qu+av_dth*av_dq ))/av_thvu

C 2.2 Prognostic variable update
C ------------------------------------------------------------

C Filter out bad wakes

      IF ( wape .LT. 0.) THEN
        if(prt_level.ge.10) print*,'wape<0'
        wape = 0.
        hw = hwmin
        sigmaw = max(sigmad,sigd_con)
        fip = 0.
        DO k = 1,klev
          deltatw(k) = 0.
          deltaqw(k) = 0.
          dth(k) = 0.
        ENDDO
      ELSE
        if(prt_level.ge.10) print*,'wape>0'
        Cstar = stark*sqrt(2.*wape)
      ENDIF

C------------------------------------------------------------------
C    Sub-time-stepping
C------------------------------------------------------------------

c      nsub=36
      nsub=10
      dtimesub=dtime/nsub

c------------------------------------------------------------
      DO isubstep = 1,nsub
c------------------------------------------------------------

c        print*,'---------------','substep=',isubstep,'-------------'

c  Evolution of sigmaw


        gfl = 2.*sqrt(3.14*wdens*sigmaw)            

        sigmaw =sigmaw + gfl*Cstar*dtimesub
        sigmaw =min(sigmaw,0.99)     !!!!!!!!
c        wdens = wdens0/(10.*sigmaw)
c        sigmaw =max(sigmaw,sigd_con)
c        sigmaw =max(sigmaw,sigmad)

c calcul de la difference de vitesse verticale poche - zone non perturbee

        z= 0.
        dp_deltomg(1:klev)=0.
        omg(1:klev+1)=0.

        omg(1) = 0.
        dp_deltomg(1) = -(gfl*Cstar)/(sigmaw * (1-sigmaw))

        DO k=2,ktop
          dz = -(Ph(k)-Ph(k-1))/(rho(k-1)*rg)
          z = z+dz
          dp_deltomg(k)= dp_deltomg(1)
          omg(k)= dp_deltomg(1)*z
        ENDDO

        dztop=-(Ptop-Ph(ktop))/(rho(ktop)*rg)
        ztop = z+dztop
        omgtop=dp_deltomg(1)*ztop


c Conversion de la vitesse verticale de m/s a Pa/s

        omgtop = -rho(ktop)*rg*omgtop
        dp_deltomg(1) = omgtop/(ptop-ph(1))

        DO k = 1,ktop
          omg(k) = - rho(k)*rg*omg(k)
          dp_deltomg(k) = dp_deltomg(1)
        ENDDO

c   raccordement lineaire de omg de ptop a pupper

      IF (kupper .GT. ktop) THEN
        omg(kupper+1) = - Rg*amdwn(kupper+1)/sigmaw
     $                + Rg*amup(kupper+1)/(1.-sigmaw)
        dp_deltomg(kupper) = (omgtop-omg(kupper+1))/(Ptop-Pupper)
        DO k=ktop+1,kupper
          dp_deltomg(k) = dp_deltomg(kupper)
          omg(k) = omgtop+(ph(k)-Ptop)*dp_deltomg(kupper)
        ENDDO
      ENDIF

c   Compute wake average vertical velocity omgbw

      DO k = 1,klev+1
        omgbw(k) = omgb(k)+(1.-sigmaw)*omg(k)
      ENDDO

c  and its vertical gradient dp_omgbw

      DO k = 1,klev
        dp_omgbw(k) = (omgbw(k+1)-omgbw(k))/(ph(k+1)-ph(k))
      ENDDO


c  Upstream coefficients for omgb velocity
c--    (alpha_up(k) is the coefficient of the value at level k)
c--    (1-alpha_up(k) is the coefficient of the value at level k-1)

      DO k = 1,klev
       alpha_up(k) = 0.
       IF (omgb(k) .GT. 0.) alpha_up(k) = 1.
      ENDDO

c  Matrix expressing [The,deltatw] from  [Th1,Th2]

      RRe1 = 1.-sigmaw
      RRe2 = sigmaw
      RRd1 = -1.
      RRd2 = 1.

c Get [Th1,Th2], dth and [q1,q2]

      DO k = 1,kupper+1
        dth(k) = deltatw(k)/ppi(k)
        Th1(k) = the(k) - sigmaw     *dth(k)   ! undisturbed area
        Th2(k) = the(k) + (1.-sigmaw)*dth(k)   ! wake
        q1(k) = qe(k) - sigmaw     *deltaqw(k) ! undisturbed area
        q2(k) = qe(k) + (1.-sigmaw)*deltaqw(k) ! wake
      ENDDO

      D_Th1(1) = 0.
      D_Th2(1) = 0.
      D_dth(1) = 0.
      D_q1(1) = 0.
      D_q2(1) = 0.
      D_dq(1) = 0.

      DO k = 2,kupper+1 !   loop on interfaces
        D_Th1(k) = Th1(k-1)-Th1(k)
        D_Th2(k) = Th2(k-1)-Th2(k)
        D_dth(k) = dth(k-1)-dth(k)
        D_q1(k) = q1(k-1)-q1(k)
        D_q2(k) = q2(k-1)-q2(k)
        D_dq(k) = deltaqw(k-1)-deltaqw(k)
      ENDDO

      omgbdth(1) = 0.
      omgbdq(1) = 0.

      DO k = 2,kupper+1  !   loop on interfaces
        omgbdth(k) = omgb(k)*(    dth(k-1) -     dth(k))
        omgbdq(k)  = omgb(k)*(deltaqw(k-1) - deltaqw(k))
      ENDDO


c-----------------------------------------------------------------
      DO k=1,kupper-1
c-----------------------------------------------------------------
c
c   Compute redistribution (advective) term
c
         d_deltatw(k) =
     $             dtimesub/(Ph(k)-Ph(k+1))*(
     $       RRd1*omg(k  )*sigmaw     *D_Th1(k)
     $      -RRd2*omg(k+1)*(1.-sigmaw)*D_Th2(k+1)
     $      -(1.-alpha_up(k))*omgbdth(k) - alpha_up(k+1)*omgbdth(k+1)
     $                      )*ppi(k)
c         print*,'d_deltatw=',d_deltatw(k)
c
         d_deltaqw(k) =
     $             dtimesub/(Ph(k)-Ph(k+1))*(
     $       RRd1*omg(k  )*sigmaw     *D_q1(k)
     $      -RRd2*omg(k+1)*(1.-sigmaw)*D_q2(k+1)
     $      -(1.-alpha_up(k))*omgbdq(k) - alpha_up(k+1)*omgbdq(k+1)
     $                      )
c         print*,'d_deltaqw=',d_deltaqw(k)
c
c   and increment large scale tendencies
c
         dtls(k) = dtls(k) +
     $               dtimesub*(
     $        ( RRe1*omg(k  )*sigmaw     *D_Th1(k)
     $         -RRe2*omg(k+1)*(1.-sigmaw)*D_Th2(k+1) )
     $               /(Ph(k)-Ph(k+1))
     $         -sigmaw*(1.-sigmaw)*dth(k)*dp_deltomg(k)
     $                      )*ppi(k)
c         print*,'dtls=',dtls(k)
c
         dqls(k) = dqls(k) +
     $               dtimesub*(
     $        ( RRe1*omg(k  )*sigmaw     *D_q1(k)
     $         -RRe2*omg(k+1)*(1.-sigmaw)*D_q2(k+1) )
     $               /(Ph(k)-Ph(k+1))
     $         -sigmaw*(1.-sigmaw)*deltaqw(k)*dp_deltomg(k)
     $                      )
c         print*,'dqls=',dqls(k)

c-------------------------------------------------------------------
      ENDDO
c------------------------------------------------------------------

C   Increment state variables

      DO k = 1,kupper-1

c Coefficient de répartition

        Crep(k)=Crep_sol*(ph(kupper)-ph(k))/(ph(kupper)-ph(1))
        Crep(k)=Crep(k)+Crep_upper*(ph(1)-ph(k))/(p(1)-ph(kupper))
        

c Reintroduce compensating subsidence term.

c        dtKE(k)=(dtdwn(k)*Crep(k))/sigmaw
c        dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k))
c     .                   /(1-sigmaw)
c        dqKE(k)=(dqdwn(k)*Crep(k))/sigmaw
c        dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k))
c     .                   /(1-sigmaw)
c
c        dtKE(k)=(dtdwn(k)*Crep(k)+(1-Crep(k))*dta(k))/sigmaw
c        dtKE(k)=dtKE(k)-(dtdwn(k)*(1-Crep(k))+dta(k)*Crep(k))
c     .                   /(1-sigmaw)
c        dqKE(k)=(dqdwn(k)*Crep(k)+(1-Crep(k))*dqa(k))/sigmaw
c        dqKE(k)=dqKE(k)-(dqdwn(k)*(1-Crep(k))+dqa(k)*Crep(k))
c     .                   /(1-sigmaw)

        dtKE(k)=(dtdwn(k)/sigmaw - dta(k)/(1.-sigmaw))
        dqKE(k)=(dqdwn(k)/sigmaw - dqa(k)/(1.-sigmaw))
c        print*,'dtKE=',dtKE(k)
c        print*,'dqKE=',dqKE(k)
c
        dtPBL(k)=(wdtPBL(k)/sigmaw - udtPBL(k)/(1.-sigmaw))
        dqPBL(k)=(wdqPBL(k)/sigmaw - udqPBL(k)/(1.-sigmaw))
c
        spread(k) = (1.-sigmaw)*dp_deltomg(k)+gfl*Cstar/sigmaw
c        print*,'spread=',spread(k)


c ajout d'un effet onde de gravité -Tgw(k)*deltatw(k) 03/02/06 YU Jingmei

        d_deltat_gw(k)=d_deltat_gw(k)-Tgw(k)*deltatw(k)* dtimesub
c        print*,'d_delta_gw=',d_deltat_gw(k)
        ff=d_deltatw(k)/dtimesub

c Sans GW
c
c        deltatw(k)=deltatw(k)+dtimesub*(ff+dtKE(k)-spread(k)*deltatw(k)) 
c
c GW formule 1
c
c        deltatw(k) = deltatw(k)+dtimesub*
c     $         (ff+dtKE(k) - spread(k)*deltatw(k)-Tgw(k)*deltatw(k))
c
c GW formule 2

        IF (dtimesub*Tgw(k).lt.1.e-10) THEN
          deltatw(k) = deltatw(k)+dtimesub*
     $          (ff+dtKE(k)+dtPBL(k) 
     $          - spread(k)*deltatw(k)-Tgw(k)*deltatw(k))
        ELSE
           deltatw(k) = deltatw(k)+1/Tgw(k)*(1-exp(-dtimesub*Tgw(k)))*
     $          (ff+dtKE(k)+dtPBL(k)
     $          - spread(k)*deltatw(k)-Tgw(k)*deltatw(k))
        ENDIF
   
        dth(k) = deltatw(k)/ppi(k)

        gg=d_deltaqw(k)/dtimesub

       deltaqw(k) = deltaqw(k) +
     $         dtimesub*(gg+ dqKE(k)+dqPBL(k) - spread(k)*deltaqw(k))

       d_deltatw2(k)=d_deltatw2(k)+d_deltatw(k)
       d_deltaqw2(k)=d_deltaqw2(k)+d_deltaqw(k)
      ENDDO

C   And update large scale variables

      DO k = 1,kupper
        te(k) = te0(k) + dtls(k)
        qe(k) = qe0(k) + dqls(k)
        the(k) = te(k)/ppi(k)
      ENDDO

c     Determine Ptop from buoyancy integral
c----------------------------------------------------------------------

c-1/ Pressure of the level where dth changes sign.

      Ptop_provis=ph(1)

      DO k= 2,klev
        IF (dth(k) .GT. -delta_t_min .and.
     $      dth(k-1).LT. -delta_t_min) THEN
          Ptop_provis = ((dth(k)+delta_t_min)*p(k-1)
     $          - (dth(k-1)+delta_t_min)*p(k)) /(dth(k) - dth(k-1))
        GO TO 65
        ENDIF
      ENDDO
65    CONTINUE

c-2/ dth integral

      sum_dth = 0.
      dthmin = -delta_t_min
      z = 0.

      DO k = 1,klev
        dz = -(max(ph(k+1),Ptop_provis)-Ph(k))/(rho(k)*rg)
        IF (dz .le. 0) GO TO 70
        z = z+dz
        sum_dth = sum_dth + dth(k)*dz
        dthmin = min(dthmin,dth(k))
      ENDDO
70    CONTINUE

c-3/ height of triangle with area= sum_dth and base = dthmin

      hw = 2.*sum_dth/min(dthmin,-0.5)
      hw = max(hwmin,hw)

c-4/ now, get Ptop

      ktop = 0
      z=0.

      DO k = 1,klev
        dz = min(-(ph(k+1)-Ph(k))/(rho(k)*rg),hw-z)
        IF (dz .le. 0) GO TO 75
        z = z+dz
        Ptop = Ph(k)-rho(k)*rg*dz
        ktop = k
      ENDDO
75    CONTINUE

c-5/Correct ktop and ptop

      Ptop_new=ptop

      DO k= ktop,2,-1
        IF (dth(k) .GT. -delta_t_min .and.
     $      dth(k-1).LT. -delta_t_min) THEN
          Ptop_new = ((dth(k)+delta_t_min)*p(k-1)
     $          - (dth(k-1)+delta_t_min)*p(k)) /(dth(k) - dth(k-1))
          GO TO 275
        ENDIF
      ENDDO
275   CONTINUE

      ptop = ptop_new

      DO k=klev,1,-1
        IF (ph(k+1) .LT. ptop) ktop=k
      ENDDO

c-6/ Set deltatw & deltaqw to 0 above kupper

      DO k = kupper,klev
        deltatw(k) = 0.
        deltaqw(k) = 0.
      ENDDO

c------------------------------------------------------------------
      ENDDO      ! end sub-timestep loop
C -----------------------------------------------------------------

c   Get back to tendencies per second

      DO k = 1,kupper-1
        dtls(k) = dtls(k)/dtime
        dqls(k) = dqls(k)/dtime
        d_deltatw2(k)=d_deltatw2(k)/dtime
        d_deltaqw2(k)=d_deltaqw2(k)/dtime
        d_deltat_gw(k) = d_deltat_gw(k)/dtime
      ENDDO

C 2.1 - Undisturbed area and Wake integrals
C ---------------------------------------------------------

      z = 0.
      sum_thu = 0.
      sum_tu = 0.
      sum_qu = 0.
      sum_thvu = 0.
      sum_dth = 0.
      sum_dq = 0.
      sum_rho = 0.
      sum_dtdwn = 0.
      sum_dqdwn = 0.

      av_thu = 0.
      av_tu =0.
      av_qu =0.
      av_thvu = 0.
      av_dth = 0.
      av_dq = 0.
      av_rho =0.
      av_dtdwn =0.
      av_dqdwn = 0.

C Potential temperatures and humidity
c----------------------------------------------------------

      DO k =1,klev
        rho(k) = p(k)/(rd*te(k))
        IF(k .eq. 1) THEN
          rhoh(k) = ph(k)/(rd*te(k))
          zhh(k)=0
        ELSE
          rhoh(k) = ph(k)*2./(rd*(te(k)+te(k-1)))
          zhh(k)=(ph(k)-ph(k-1))/(-rhoh(k)*RG)+zhh(k-1)
        ENDIF
        the(k) = te(k)/ppi(k)
        thu(k) = (te(k) - deltatw(k)*sigmaw)/ppi(k)
        tu(k) = te(k) - deltatw(k)*sigmaw
        qu(k)  =  qe(k) - deltaqw(k)*sigmaw
        rhow(k) = p(k)/(rd*(te(k)+deltatw(k)))
        dth(k) = deltatw(k)/ppi(k)
       
      ENDDO

C Integrals (and wake top level number)
C -----------------------------------------------------------

C Initialize sum_thvu to 1st level virt. pot. temp.

      z = 1.
      dz = 1.
      sum_thvu =  thu(1)*(1.+eps*qu(1))*dz
      sum_dth = 0.

      DO k = 1,klev
        dz = -(max(ph(k+1),Ptop)-Ph(k))/(rho(k)*rg)

        IF (dz .LE. 0) GO TO 51
        z = z+dz
        sum_thu = sum_thu + thu(k)*dz
        sum_tu = sum_tu + tu(k)*dz
        sum_qu = sum_qu + qu(k)*dz
        sum_thvu = sum_thvu + thu(k)*(1.+eps*qu(k))*dz
        sum_dth = sum_dth + dth(k)*dz
        sum_dq = sum_dq + deltaqw(k)*dz
        sum_rho = sum_rho + rhow(k)*dz
        sum_dtdwn = sum_dtdwn + dtdwn(k)*dz
        sum_dqdwn = sum_dqdwn + dqdwn(k)*dz
      ENDDO
 51   CONTINUE

      hw0 = z

C 2.1 - WAPE and mean forcing computation
C-------------------------------------------------------------

C Means

      av_thu = sum_thu/hw0
      av_tu = sum_tu/hw0
      av_qu = sum_qu/hw0
      av_thvu = sum_thvu/hw0
      av_dth = sum_dth/hw0
      av_dq = sum_dq/hw0
      av_rho = sum_rho/hw0
      av_dtdwn = sum_dtdwn/hw0
      av_dqdwn = sum_dqdwn/hw0

      wape2 = - rg*hw0*(av_dth
     $     + eps*(av_thu*av_dq+av_dth*av_qu+av_dth*av_dq ))/av_thvu


C 2.2 Prognostic variable update
C ------------------------------------------------------------

C Filter out bad wakes

      IF ( wape2 .LT. 0.) THEN
        if(prt_level.ge.10) print*,'wape2<0'
        wape2 = 0.
        hw = hwmin
        sigmaw = max(sigmad,sigd_con)
        fip = 0.
        DO k = 1,klev
          deltatw(k) = 0.
          deltaqw(k) = 0.
          dth(k) = 0.
        ENDDO
      ELSE
        if(prt_level.ge.10) print*,'wape2>0'
        Cstar2 = stark*sqrt(2.*wape2)

      ENDIF

      ktopw = ktop

      IF (ktopw .gt. 0) then

Cjyg1     Utilisation d'un h_efficace constant ( ~ feeding layer)
ccc       heff = 600.
C      Utilisation de la hauteur hw
cc       heff = 0.7*hw
       heff = hw

       FIP = 0.5*rho(ktopw)*Cstar2**3*heff*2*sqrt(sigmaw*wdens*3.14)
       FIP = alpk * FIP
Cjyg2
       ELSE
         FIP = 0.
       ENDIF


C   Limitation de sigmaw
c
C   sécurité : si le wake occuppe plus de 90 % de la surface de la maille,
C              alors il disparait en se mélangeant à la partie undisturbed

! correction NICOLAS     .     ((wape.ge.wape2).and.(wape2.le.1.0))) THEN
      IF ((sigmaw.GT.0.9).or.
     .     ((wape.ge.wape2).and.(wape2.le.1.0)).or.(ktopw.le.2)) THEN
cIM cf NR/JYG 251108    .     ((wape.ge.wape2).and.(wape2.le.1.0))) THEN
c      IF (sigmaw.GT.0.9) THEN
        DO k = 1,klev
          dtls(k) = 0.
          dqls(k) = 0.
          deltatw(k) = 0.
          deltaqw(k) = 0.
        ENDDO
        wape = 0.
        hw = hwmin
        sigmaw = sigmad
        fip = 0.
      ENDIF

      RETURN
      END



