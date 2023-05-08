!******************************************************************************
!*     venus_SAS_composition SUBROUTINE
!*     modified from
!*     PROGRAM PSC_MODEL_E
!*     by A. Määttänen
!*     subroutine for LMDZ+photochemistry VENUS
!*     by A. Stolzenbach
!*
!*     Input/Output files:
!*     -------------------
!*
!----------------------------------------------------------------------------
      SUBROUTINE new_cloud_venus(
     + nblev, nblon,
     + TT,PP,
     + mrt_wv,mrt_sa,
     + mr_wv,mr_sa)

      USE chemparam_mod
      IMPLICIT NONE

#include "YOMCST.h"

      INTEGER, INTENT(IN) :: nblon  ! nombre de points horizontaux
      INTEGER, INTENT(IN) :: nblev  ! nombre de couches verticales

!----------------------------------------------------------------------------
!     Ambient air state variables:
      REAL, INTENT(IN), DIMENSION(nblon,nblev) :: mrt_wv,mrt_sa,
     +                                            TT,PP
      REAL, INTENT(INOUT), DIMENSION(nblon,nblev) :: mr_wv,mr_sa
!----------------------------------------------------------------------------
      INTEGER :: ilon, ilev, imode
!----------------------------------------------------------------------------
!     Thermodynamic functions:
      REAL :: ROSAS
!----------------------------------------------------------------------------
!     Auxilary variables:
      REAL :: NH2SO4,NH2O
      REAL :: H2SO4_liq,H2O_liq
      REAL :: CONCM
      REAL :: MCONDTOT
      REAL :: RMODE
      REAL :: WSAFLAG
      REAL :: K_SAV
!----------------------------------------------------------------------------
!     Ridder's Method variables:
      REAL :: WVMIN, WVMAX, WVACC

      INTEGER :: NBROOT

      INTEGER :: MAXITE
      PARAMETER(MAXITE=20)

      INTEGER :: NBRAC
      PARAMETER(NBRAC=5)

      INTEGER :: FLAG
!----------------------------------------------------------------------------
c     Ratio radius shell model du mode 3
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	Ce ratio correspond aux mesures effectuées par J. Cimino (1982), Icarus
c     Fixer ce parametre a 0 revient a une gouttelette pure en liquide acide sulfurique
c     ATTENTION ! DOIT ETRE COHERENT AVEC new_cloud_sedim !
      REAL, PARAMETER :: qrad = 0.97
      REAL :: qmass
c	masse volumique du coeur (kg.m-3)
c     ATTENTION ! DOIT ETRE COHERENT AVEC new_cloud_sedim !
      REAL, PARAMETER :: rho_core = 2500.0
!----------------------------------------------------------------------------
!     External functions needed:
      REAL :: IRFRMWV
!----------------------------------------------------------------------------


! >>> Program starts here:

!AM Venus
! These aerosols will then be given an equilibrium composition for the given size distribution

  ! Hanna Vehkamäki and Markku Kulmala and Ismo Napari
  ! and Kari E. J. Lehtinen and Claudia Timmreck and Madis Noppel and Ari Laaksonen, 2002,
  ! An improved parameterization for sulfuric acid/water nucleation rates for tropospheric
  !and stratospheric conditions, () J. Geophys. Res., 107, PP. 4622-4631

!===========================================
!     Debut boucle sur niveau et lat,lon
!===========================================
!     Init, tous les points=0, cela met les niveaux > cloudmax et < cloudmin a 0
      NBRTOT(:,:,:)=0.0E+0
      WH2SO4(:,:)=0.0E+0
      rho_droplet(:,:)=0.0E+0

      DO ilev=cloudmin, cloudmax
      DO ilon=1, nblon

!       Boucle sur les modes
        RMODE=0.0E+0
        K_SAV = 0.0

        DO imode=1, nbr_mode
          IF (K_MASS(ilon,ilev,imode).GT.K_SAV) THEN
!       RMODE est le rayon modal de la distribution en volume du mode le plus
!       representatif pour la Mtot
            RMODE=R_MEDIAN(ilon,ilev,imode)*
     &      EXP(2.*(DLOG(STDDEV(ilon,ilev,imode))**2.))
            K_SAV=K_MASS(ilon,ilev,imode)
          ENDIF
        ENDDO ! FIN boucle imode

!       Initialisation des bornes pour WV
        WVMIN=1.E-90
        WVMAX=mrt_wv(ilon,ilev)

!       Accuracy de WVeq
        WVACC=WVMAX*1.0E-3

!       BRACWV borne la fonction f(WV) - WV = 0
!       de WV=0 à WV=WVtot on cherche l'intervalle où f(WV) - WV = 0
!       avec précisément f(WVliq de WSA<=WVinput) + WVinput - WVtot = 0
!       Elle fait appel à la fct/ssrtine ITERWV()

        CALL BRACWV(TT(ilon,ilev),PP(ilon,ilev),WVMIN,WVMAX,NBRAC,
     &  RMODE,mrt_wv(ilon,ilev),mrt_sa(ilon,ilev),FLAG,WSAFLAG,NBROOT)

        SELECT CASE(FLAG)

        CASE(1)
!         Cas NROOT=1 ou NROOT>1 mais dans un intervalle restreint WVTOT (cas courant)
!       IRFRMWV Ridder's method pour trouver, sur [WVmin,WVmax], WVo tel que f(WVo) - WVo = 0
!       Elle fait appel  la fct/ssrtine ITERWV()

          WH2SO4(ilon,ilev)=IRFRMWV(TT(ilon,ilev),PP(ilon,ilev),
     &    WVMIN,WVMAX,WVACC,MAXITE,RMODE,
     &    mrt_wv(ilon,ilev),mrt_sa(ilon,ilev),NBROOT)

          rho_droplet(ilon,ilev)=ROSAS(TT(ilon,ilev),WH2SO4(ilon,ilev))

!          IF (rho_droplet(ilon,ilev).LT.1100.) THEN
!            PRINT*,'PROBLEM RHO_DROPLET'
!            PRINT*,'rho_droplet',rho_droplet(ilon,ilev)
!            PRINT*,'T',TT(ilon,ilev),'WSA',WH2SO4(ilon,ilev)
!            PRINT*,'RHODROPLET',ROSAS(TT(ilon,ilev),WH2SO4(ilon,ilev))
!            PRINT*,'FLAG',FLAG,'NROOT',NBROOT
!            STOP
!          ENDIF

          CONCM= PP(ilon,ilev)/(1.3806488E-23*TT(ilon,ilev)) !air number density, molec/m3

	      NH2SO4=mrt_sa(ilon,ilev)*CONCM
	      NH2O=mrt_wv(ilon,ilev)*CONCM

          CALL CALCM_SAT(NH2SO4,NH2O,WH2SO4(ilon,ilev),
     &       rho_droplet(ilon,ilev),TT(ilon,ilev),
     &       H2SO4_liq,H2O_liq,MCONDTOT)

!       Boucle sur les modes
! *********************************************
!	AVEC MODE 3 type J. Cimino 1982, Icarus
! *********************************************
!	Mode 1 et 2 avec les Kmass modifies
          DO imode=1, nbr_mode-1
!           calcul qmass
	      qmass = (rho_core*qrad**3)/
     &      (rho_core*qrad**3 + rho_droplet(ilon,ilev)*(1.-qrad**3))

            IF (K_MASS(ilon,ilev,imode).GT.0.) THEN
              NBRTOT(ilon,ilev,imode)= 1.E-6*3./(4.*RPI)*
     &        K_MASS(ilon,ilev,imode)/(1.-qmass*K_MASS(ilon,ilev,3))*
     &        MCONDTOT*EXP(-4.5*DLOG(STDDEV(ilon,ilev,imode))**2.)/
     &        (R_MEDIAN(ilon,ilev,imode)**3.)
            ELSE
              NBRTOT(ilon,ilev,imode)=0.0E+0
            ENDIF
          ENDDO
!	Mode 3 reste identique car on veut le N correspondant aux mesures de Knollenberg
          IF (K_MASS(ilon,ilev,3).GT.0.) THEN
            NBRTOT(ilon,ilev,3)= 1.E-6*3./(4.*RPI)*
     &      K_MASS(ilon,ilev,3)*MCONDTOT*
     &      EXP(-4.5*DLOG(STDDEV(ilon,ilev,3))**2.)/
     &      (R_MEDIAN(ilon,ilev,3)**3.)
          ELSE
              NBRTOT(ilon,ilev,3)=0.0E+0
          ENDIF

!       Passage de #/m3 en VMR
          H2O_liq=H2O_liq/CONCM
          H2SO4_liq=H2SO4_liq/CONCM

          mr_wv(ilon,ilev)=mrt_wv(ilon,ilev)-H2O_liq
          mr_sa(ilon,ilev)=mrt_sa(ilon,ilev)-H2SO4_liq

!          Problemes quand on a condense tout, on peut obtenir des -1e-24
!		aprs la soustraction et conversion de ND  VMR
          IF (mr_wv(ilon,ilev).LE.0.0) mr_wv(ilon,ilev)=1.0E-30
          IF (mr_sa(ilon,ilev).LE.0.0) mr_sa(ilon,ilev)=1.0E-30



        CASE(2)
!       Cas NROOT=0 mais proche de 0

          WH2SO4(ilon,ilev)=WSAFLAG

          rho_droplet(ilon,ilev)=ROSAS(TT(ilon,ilev),WH2SO4(ilon,ilev))

!     ATTENTION ce IF ne sert a rien en fait, juste a retenir une situation
!     ubuesque dans mon code ou sans ce IF les valeurs de rho_droplets sont
!     incoherentes avec TT et WH2SO4 (a priori lorsque NTOT=0)
!     Juste le fait de METTRE un IF fait que rho_droplet a la bonne valeur
!     donne par ROSAS (cf test externe en Python), sinon, la valeur est trop
!     basse (de l'ordre de 1000 kg/m3) et correspond parfois a la valeur avec
!     WSA=0.1 (pas totalement sur)
!     En tous cas, incoherent avec ce qui est attendue pour le WSA et T donnee
!     La version avec le IF (rho<1100 & WSA>0.1) est CORRECTE, rho_droplet a
!     la bonne valeur (tests externes Python confirment)

          IF ((rho_droplet(ilon,ilev).LT.1100.).AND.
     &      (WH2SO4(ilon,ilev).GT.0.1))THEN
            PRINT*,'PROBLEM RHO_DROPLET'
            PRINT*,'rho_droplet',rho_droplet(ilon,ilev)
            PRINT*,'T',TT(ilon,ilev),'WSA',WH2SO4(ilon,ilev)
            PRINT*,'RHODROPLET',ROSAS(TT(ilon,ilev),WH2SO4(ilon,ilev))
            PRINT*,'FLAG',FLAG,'NROOT',NBROOT
            STOP
          ENDIF


          CONCM= PP(ilon,ilev)/(1.3806488E-23*TT(ilon,ilev)) !air number density, molec/m3

	    NH2SO4=mrt_sa(ilon,ilev)*CONCM
	    NH2O=mrt_wv(ilon,ilev)*CONCM

          CALL CALCM_SAT(NH2SO4,NH2O,WH2SO4(ilon,ilev),
     &       rho_droplet(ilon,ilev),TT(ilon,ilev),
     &       H2SO4_liq,H2O_liq,MCONDTOT)

!       Boucle sur les modes
! *********************************************
!	AVEC MODE 3 type J. Cimino 1982, Icarus
! *********************************************
!	Mode 1 et 2 avec alcul coeff*Kmass
          DO imode=1, nbr_mode-1
!           calcul qmass
	      qmass = (rho_core*qrad**3)/
     &      (rho_core*qrad**3 + rho_droplet(ilon,ilev)*(1.-qrad**3))

            IF (K_MASS(ilon,ilev,imode).GT.0.) THEN
              NBRTOT(ilon,ilev,imode)= 1.E-6*3./(4.*RPI)*
     &        K_MASS(ilon,ilev,imode)/(1.-qmass*K_MASS(ilon,ilev,3))*
     &        MCONDTOT*EXP(-4.5*DLOG(STDDEV(ilon,ilev,imode))**2.)/
     &        (R_MEDIAN(ilon,ilev,imode)**3.)
            ELSE
              NBRTOT(ilon,ilev,imode)=0.0E+0
            ENDIF
          ENDDO
!	Mode 3 reste identique car on veut le N correspondant aux mesures de Knollenberg
          IF (K_MASS(ilon,ilev,3).GT.0.) THEN
            NBRTOT(ilon,ilev,3)= 1.E-6*3./(4.*RPI)*
     &      K_MASS(ilon,ilev,3)*MCONDTOT*
     &      EXP(-4.5*DLOG(STDDEV(ilon,ilev,3))**2.)/
     &      (R_MEDIAN(ilon,ilev,3)**3.)
          ELSE
              NBRTOT(ilon,ilev,3)=0.0E+0
          ENDIF


!       Passage de #/m3 en VMR
          H2O_liq=H2O_liq/CONCM
          H2SO4_liq=H2SO4_liq/CONCM

          mr_wv(ilon,ilev)=mrt_wv(ilon,ilev)-H2O_liq
          mr_sa(ilon,ilev)=mrt_sa(ilon,ilev)-H2SO4_liq

!          Problmes quand on a condense tout, on peut obtenir des -1e-24
!		aprs la soustraction et conversion de ND  VMR
          IF (mr_wv(ilon,ilev).LE.0.0) mr_wv(ilon,ilev)=1.0E-30
          IF (mr_sa(ilon,ilev).LE.0.0) mr_sa(ilon,ilev)=1.0E-30

        CASE(3)
!         Cas 0 NROOT
            mr_wv(ilon,ilev)=mrt_wv(ilon,ilev)
            mr_sa(ilon,ilev)=mrt_sa(ilon,ilev)
            rho_droplet(ilon,ilev)=0.0E+0
            WH2SO4(ilon,ilev)=0.0E+0
            DO imode=1, nbr_mode
              NBRTOT(ilon,ilev,imode)=0.0E+0
            ENDDO

        END SELECT
      ENDDO   !FIN boucle ilon
      ENDDO   !FIN boucle ilev

      END SUBROUTINE new_cloud_venus


!******************************************************************
      SUBROUTINE CALCM_SAT(H2SO4,H2O,WSA,DENSO4,
     + T,H2SO4COND,H2OCOND,RMTOT)

!     DERIVE NO (TOTAL NUMBER OF AEROSOL PARTICLES CONCENTRATION)
!     FROM TOTAL H2SO4 AND RMOD/SIGMA OF AEROSOL LOG-NORMAL
!                                       SIZE DISTRIBTUION
!     ASSUMING ALL THE H2SO4 ABOVE MIXTURE SAT PRESSURE modified by H2SO4 activity IS CONDENSED
!    ---------------------------------------------------------------
!     INPUT:
!     H2SO4: #/m3 of total H2SO4
!	H2O  : #/m3 of total H2O
!     WSA: aerosol H2SO4 weight fraction (fraction)
!     DENSO4: aerosol volumic mass (kg/m3 = aerosol mass/aerosol volume)
!       for total mass, almost same result with ro=1.67 gr/cm3
!     RSTDEV: standard deviation of aerosol distribution (no unit)
!     RADIUS: MEDIAN radius (m)
!     T: temperature (K)
!
!     OUTPUT:
!     RMTOT: Total condensed "Mass" (M_tot_distrib / rho_droplet), sans dimension
!            mais rho_droplet et M_tot_distrib doivent etre de meme dimension
!	H2OCOND
!	H2SO4COND



      IMPLICIT NONE

      REAL, INTENT(IN) :: H2SO4, H2O, WSA
      REAL, INTENT(IN) :: DENSO4, T
      REAL, INTENT(OUT) :: H2OCOND, H2SO4COND, RMTOT
!     working variables
      REAL :: RMH2S4
      REAL :: DND2,pstand,lpar,acidps
      REAL :: x1, satpacid
      REAL , DIMENSION(2):: act
!
!     masse of an H2SO4 molecule (kg)
      RMH2S4=98.078/(6.02214129E+26)

      pstand=1.01325E+5 !Pa  1 atm pressure

        x1=(WSA/98.08)/(WSA/98.08 + ((1.-WSA)/18.0153))

        call zeleznik(x1,t,act)

!pure acid satur vapor pressure
        lpar= -11.695+DLOG(pstand) ! Zeleznik
        acidps=1/360.15-1.0/t+0.38/545.
     + *(1.0+DLOG(360.15/t)-360.15/t)
        acidps = 10156.0*acidps +lpar
        acidps = DEXP(acidps)    !Pa

!acid sat.vap.PP over mixture (flat surface):
        satpacid=act(2)*acidps ! Pa

!       Conversion from Pa to N.D #/m3
        DND2=satpacid/(1.3806488E-23*T)
!       Conversion from N.D #/m3 TO #/cm3
!        DND2=DND2*1.d-6

!	H2SO4COND N.D #/m3 condensee ssi H2SO4>H2SO4sat
	IF (H2SO4.GE.DND2) THEN
	H2SO4COND=H2SO4-DND2
!	calcul de H2O cond correspondant a H2SO4 cond
	H2OCOND=H2SO4COND*98.078*(1.0-WSA)/(18.0153*WSA)

!     RMTOT: = Mass of H2SO4 satur per m3 of air/ Mass of sulfuric acid part of droplet solution per m3
!     RMTOT = Mtot/rho_droplet
!	RMTOT=M_distrib/rho_droplet

	RMTOT=H2SO4COND*RMH2S4/(DENSO4*WSA)

!	Si on a H2SO4<H2SO4sat on ne condense rien et NDTOT=0
	ELSE
	H2SO4COND=0.0E+0
	H2OCOND=0.0E+0
	RMTOT=0.0E+0
	END IF

!	Test si H2O en defaut H2Ocond>H2O dispo
	IF ((H2OCOND.GT.H2O).AND.(H2SO4.GE.DND2)) THEN

!     Si H2O en dfaut, on as pas le bon WSA!
!     En effet, normalement, on a exactement le WSA correspondant a
!     WVg + WVl = WVtot
!     Dans les cas o WVtot, SAtot sont trs faibles (Upper Haze) ou
!     quand T est grand (Lower Haze), le modle reprsente mal le WSA
!     cf carte NCL, avec des max erreur absolue de 0.1 sur le WSA

!      PRINT*,'PROBLEM H2O EN DEFAUT'
!      PRINT*,'H2OCOND',H2OCOND,'H2O',H2O
!      PRINT*,'WSA',WSA,'RHO',DENSO4
!      STOP


!	On peut alors condenser tout le H2O dispo
	H2OCOND=H2O
!	On met alors egalement a jour le H2SO4 cond correspondant au H2O cond
	H2SO4COND=H2OCOND*18.0153*WSA/(98.078*(1.0-WSA))

!     RMTOT: = Mass of H2SO4 satur per m3 of air/ Mass of sulfuric acid part of droplet solution per m3
!	RMTOT=Volume of aerosol m3 /m3 of air
!	Volume of aerosol/m3 air

	RMTOT=H2SO4COND*RMH2S4/(DENSO4*WSA)

      END IF

      END SUBROUTINE CALCM_SAT

