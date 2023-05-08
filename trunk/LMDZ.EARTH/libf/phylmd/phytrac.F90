!$Id $

SUBROUTINE phytrac(                            &
     nstep,     julien,   gmtime,   debutphy,  &
     lafin,     pdtphys,  u, v,     t_seri,     &
     paprs,     pplay,    pmfu,     pmfd,      &
     pen_u,     pde_u,    pen_d,    pde_d,     &
     cdragh,    coefh,    fm_therm, entr_therm,&
     yu1,       yv1,      ftsol,    pctsrf,    &
     xlat,      frac_impa,frac_nucl,xlon,      &
     presnivs,  pphis,    pphi,     albsol,    &
     sh,        rh,       cldfra,   rneb,      &
     diafra,    cldliq,   itop_con, ibas_con,  &
     pmflxr,    pmflxs,   prfl,     psfl,      &
     da,        phi,      mp,       upwd,      &
     dnwd,      aerosol_couple,     flxmass_w, &
     tau_aero,  piz_aero,  cg_aero, ccm,       &
     rfname,                                   &
     tr_seri)         
!
!======================================================================
! Auteur(s) FH
! Objet: Moniteur general des tendances traceurs
!======================================================================

  USE ioipsl
  USE dimphy
  USE infotrac
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE comgeomphy
  USE iophy
  USE traclmdz_mod
  USE tracinca_mod
  USE control_mod



  IMPLICIT NONE

  INCLUDE "YOMCST.h"
  INCLUDE "dimensions.h"
  INCLUDE "indicesol.h"
  INCLUDE "clesphys.h"
  INCLUDE "temps.h"
  INCLUDE "paramet.h"
  INCLUDE "thermcell.h"
!==========================================================================
!                   -- ARGUMENT DESCRIPTION --
!==========================================================================

! Input arguments
!----------------
!Configuration grille,temps:
  INTEGER,INTENT(IN) :: nstep      ! Appel physique
  INTEGER,INTENT(IN) :: julien     ! Jour julien
  REAL,INTENT(IN)    :: gmtime
  REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
  LOGICAL,INTENT(IN) :: debutphy   ! le flag de l'initialisation de la physique
  LOGICAL,INTENT(IN) :: lafin      ! le flag de la fin de la physique
  
  REAL,DIMENSION(klon),INTENT(IN) :: xlat    ! latitudes pour chaque point 
  REAL,DIMENSION(klon),INTENT(IN) :: xlon    ! longitudes pour chaque point 
!
!Physique: 
!--------
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: u       ! variable not used
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: v       ! variable not used
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: rh      ! humidite relative
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel
  REAL,DIMENSION(klon),INTENT(IN)        :: pphis
  REAL,DIMENSION(klev),INTENT(IN)        :: presnivs 
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldliq  ! eau liquide nuageuse
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: cldfra  ! fraction nuageuse (tous les nuages)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: diafra  ! fraction nuageuse (convection ou stratus artificiels)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: rneb    ! fraction nuageuse (grande echelle)
  INTEGER,DIMENSION(klon),INTENT(IN)     :: itop_con
  INTEGER,DIMENSION(klon),INTENT(IN)     :: ibas_con
  REAL,DIMENSION(klon),INTENT(IN)        :: albsol  ! albedo surface
!
!Convection:
!----------
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfu  ! flux de masse dans le panache montant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfd  ! flux de masse dans le panache descendant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_u ! flux entraine dans le panache montant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_u ! flux detraine dans le panache montant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_d ! flux entraine dans le panache descendant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_d ! flux detraine dans le panache descendant

!...Tiedke     
  REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: pmflxr, pmflxs ! Flux precipitant de pluie, neige aux interfaces [convection]
  REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: prfl, psfl ! Flux precipitant de pluie, neige aux interfaces [large-scale]

  LOGICAL,INTENT(IN)                       :: aerosol_couple
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: flxmass_w
  REAL,DIMENSION(klon,klev,9,2),INTENT(IN) :: tau_aero
  REAL,DIMENSION(klon,klev,9,2),INTENT(IN) :: piz_aero
  REAL,DIMENSION(klon,klev,9,2),INTENT(IN) :: cg_aero
  CHARACTER(len=4),DIMENSION(9),INTENT(IN) :: rfname 
  REAL,DIMENSION(klon,klev,2),INTENT(IN)   :: ccm 
!... K.Emanuel
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: da
  REAL,DIMENSION(klon,klev,klev),INTENT(IN):: phi
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: mp
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: upwd      ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: dnwd      ! saturated downdraft mass flux
!
!Thermiques:
!----------
  REAL,DIMENSION(klon,klev+1),INTENT(IN)   :: fm_therm
  REAL,DIMENSION(klon,klev),INTENT(IN)     :: entr_therm
!
!Couche limite:
!--------------
!
  REAL,DIMENSION(klon),INTENT(IN)      :: cdragh ! coeff drag pour T et Q
  REAL,DIMENSION(klon,klev),INTENT(IN) :: coefh  ! coeff melange CL (m**2/s)
  REAL,DIMENSION(klon),INTENT(IN)      :: yu1    ! vents au premier niveau
  REAL,DIMENSION(klon),INTENT(IN)      :: yv1    ! vents au premier niveau
!
!Lessivage:
!----------
!
! pour le ON-LINE
!
  REAL,DIMENSION(klon,klev),INTENT(IN) :: frac_impa ! fraction d'aerosols non impactes
  REAL,DIMENSION(klon,klev),INTENT(IN) :: frac_nucl ! fraction d'aerosols non nuclees

! Arguments necessaires pour les sources et puits de traceur:
  REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol  ! Temperature du sol (surf)(Kelvin)
  REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf ! Pourcentage de sol (nature du sol)


! Output argument
!----------------
  REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT) :: tr_seri ! Concentration Traceur [U/KgA]  

!=======================================================================================
!                        -- LOCAL VARIABLES --
!=======================================================================================

  INTEGER :: i, k, it
  INTEGER :: nsplit

!Sources et Reservoirs de traceurs (ex:Radon):
!--------------------------------------------
!
  REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: source  ! a voir lorsque le flux de surface est prescrit 
!$OMP THREADPRIVATE(source)

!
!Entrees/Sorties: (cf ini_histrac.h et write_histrac.h)  
!---------------
  INTEGER                   :: iiq, ierr
  INTEGER                   :: nhori, nvert
  REAL                      :: zsto, zout, zjulian
  INTEGER,SAVE              :: nid_tra     ! pointe vers le fichier histrac.nc         
!$OMP THREADPRIVATE(nid_tra)
  REAL,DIMENSION(klon)      :: zx_tmp_fi2d ! variable temporaire grille physique
  INTEGER                   :: itau_w      ! pas de temps ecriture = nstep + itau_phy
  LOGICAL,PARAMETER :: ok_sync=.TRUE.

!
! Nature du traceur
!------------------
  LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: aerosol  ! aerosol(it) = true  => aerosol => lessivage
!$OMP THREADPRIVATE(aerosol)                        ! aerosol(it) = false => gaz
  REAL,DIMENSION(klon,klev)             :: delp     ! epaisseur de couche (Pa)
!
! Tendances de traceurs (Td):
!------------------------
!
  REAL,DIMENSION(klon,klev)      :: d_tr     ! Td dans l'atmosphere
  REAL,DIMENSION(klon,klev,nbtr) :: d_tr_cl  ! Td couche limite/traceur
  REAL,DIMENSION(klon,klev,nbtr) :: d_tr_cv  ! Td convection/traceur
  REAL,DIMENSION(klon,klev,nbtr) :: d_tr_th  ! Td thermique
  REAL,DIMENSION(klon,klev,nbtr) :: d_tr_lessi_impa ! Td du lessivage par impaction
  REAL,DIMENSION(klon,klev,nbtr) :: d_tr_lessi_nucl ! Td du lessivage par nucleation 
!
! Physique
!----------   
  REAL,DIMENSION(klon,klev,nbtr) :: flestottr ! flux de lessivage dans chaque couche 
  REAL,DIMENSION(klon,klev)      :: zmasse    ! densité atmosphérique Kg/m2
  REAL,DIMENSION(klon,klev)      :: ztra_th
  
!Controles:
!---------
  LOGICAL,SAVE :: couchelimite=.TRUE.
  LOGICAL,SAVE :: convection=.TRUE.
  LOGICAL,SAVE :: lessivage
!$OMP THREADPRIVATE(couchelimite,convection,lessivage)

  CHARACTER(len=8),DIMENSION(nbtr) :: solsym


!######################################################################
!                    -- INITIALIZATION --
!######################################################################
  IF (debutphy) THEN
     WRITE(*,*) 'FIRST TIME IN PHYTRAC : pdtphys(sec) = ',pdtphys,'ecrit_tra (sec) = ',ecrit_tra
     ALLOCATE( source(klon,nbtr), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('phytrac', 'pb in allocation 1',1)
     
     ALLOCATE( aerosol(nbtr), stat=ierr)
     IF (ierr /= 0) CALL abort_gcm('phytrac', 'pb in allocation 2',1)
     

     ! Initialize module for specific tracers
     SELECT CASE(type_trac)
     CASE('lmdz')
        CALL traclmdz_init(pctsrf, ftsol, tr_seri, t_seri, pplay, sh, pdtphys, aerosol, lessivage)
     CASE('inca')
        source(:,:)=0.
        CALL tracinca_init(aerosol,lessivage)
     END SELECT
!
! Initialize diagnostic output
! ----------------------------
#ifdef CPP_IOIPSL
!     INCLUDE "ini_histrac.h"
#endif
  END IF
!############################################ END INITIALIZATION #######

  DO k=1,klev
     DO i=1,klon
        zmasse(i,k)=(paprs(i,k)-paprs(i,k+1))/rg
     END DO
  END DO

!===============================================================================
!    -- Do specific treatment according to chemestry model or local LMDZ tracers
!      
!===============================================================================
  SELECT CASE(type_trac)
  CASE('lmdz')
     !    -- Traitement des traceurs avec traclmdz
     CALL traclmdz(nstep, julien, gmtime, pdtphys, t_seri, paprs, pplay, &
          cdragh,  coefh, yu1, yv1, ftsol, pctsrf, xlat, xlon, couchelimite, &
          sh, tr_seri, source, solsym, d_tr_cl, zmasse)
  CASE('inca')
     !    -- CHIMIE INCA  config_inca = aero or chem --

     CALL tracinca(&
          nstep,    julien,   gmtime,         lafin,     &
          pdtphys,  t_seri,   paprs,          pplay,     &
          pmfu,     ftsol,    pctsrf,         pphis,     &
          pphi,     albsol,   sh,             rh,        &
          cldfra,   rneb,     diafra,         cldliq,    &
          itop_con, ibas_con, pmflxr,         pmflxs,    &
          prfl,     psfl,     aerosol_couple, flxmass_w, &
          tau_aero, piz_aero, cg_aero,        ccm,       &
          rfname,                                        &
          tr_seri,  source,   solsym)      
  END SELECT

!======================================================================
!       -- Calcul de l'effet de la convection --
!======================================================================
  IF (convection) THEN
     DO it=1, nbtr
        IF ( conv_flg(it) == 0 ) CYCLE
        
        IF (iflag_con.LT.2) THEN
           d_tr_cv(:,:,:)=0.
        ELSE IF (iflag_con.EQ.2) THEN
!..Tiedke
           CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
                pplay, paprs, tr_seri(:,:,it), d_tr_cv(:,:,it))
        ELSE
!..K.Emanuel
           CALL cvltr(pdtphys, da, phi, mp, paprs,pplay, tr_seri(:,:,it),&
                upwd,dnwd,d_tr_cv(:,:,it))
        END IF

        DO k = 1, klev
           DO i = 1, klon        
              tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_cv(i,k,it)
           END DO
        END DO

        CALL minmaxqfi(tr_seri(:,:,it),0.,1.e33,'convection it = '//solsym(it))
             
     END DO ! nbtr
  END IF ! convection

!======================================================================
!    -- Calcul de l'effet des thermiques --
!======================================================================

  DO it=1,nbtr
     DO k=1,klev
        DO i=1,klon
           d_tr_th(i,k,it)=0.
           tr_seri(i,k,it)=MAX(tr_seri(i,k,it),0.)
           tr_seri(i,k,it)=MIN(tr_seri(i,k,it),1.e10)
        END DO
     END DO
  END DO
  
  IF (iflag_thermals.GT.0) THEN   
     nsplit=10
     DO it=1, nbtr
        DO isplit=1,nsplit

           CALL dqthermcell(klon,klev,pdtphys/nsplit, &
                fm_therm,entr_therm,zmasse, &
                tr_seri(1:klon,1:klev,it),d_tr,ztra_th)

           DO k=1,klev
              DO i=1,klon
                 d_tr(i,k)=pdtphys*d_tr(i,k)/nsplit
                 d_tr_th(i,k,it)=d_tr_th(i,k,it)+d_tr(i,k)
                 tr_seri(i,k,it)=MAX(tr_seri(i,k,it)+d_tr(i,k),0.)
              END DO
           END DO
        END DO ! nsplit
     END DO ! it
  END IF ! Thermiques

!======================================================================
!     -- Calcul de l'effet de la couche limite --
!======================================================================

  IF (couchelimite) THEN

     DO k = 1, klev
        DO i = 1, klon
           delp(i,k) = paprs(i,k)-paprs(i,k+1)
        END DO
     END DO

     DO it=1, nbtr
        
        IF( pbl_flg(it) /= 0 ) THEN
        
           CALL cltrac(pdtphys, coefh,t_seri,       &
                tr_seri(:,:,it), source(:,it),      &
                paprs, pplay, delp,                 &
                d_tr_cl(:,:,it))
           
           DO k = 1, klev
              DO i = 1, klon
                 tr_seri(i,k,it) = tr_seri(i,k,it) + d_tr_cl(i,k,it)
              END DO
           END DO
        END IF

     END DO
     
  END IF ! couche limite


!======================================================================
!   Calcul de l'effet de la precipitation
!======================================================================

  IF (lessivage) THEN
     
     d_tr_lessi_nucl(:,:,:) = 0. 
     d_tr_lessi_impa(:,:,:) = 0.
     flestottr(:,:,:) = 0. 
!=========================
! LESSIVAGE LARGE SCALE : 
!=========================

! Tendance des aerosols nuclees et impactes 
! -----------------------------------------
     DO it = 1, nbtr
        IF (aerosol(it)) THEN
           DO k = 1, klev
              DO i = 1, klon
                 d_tr_lessi_nucl(i,k,it) = d_tr_lessi_nucl(i,k,it) +    &
                      ( 1 - frac_nucl(i,k) )*tr_seri(i,k,it)
                 d_tr_lessi_impa(i,k,it) = d_tr_lessi_impa(i,k,it) +    &
                      ( 1 - frac_impa(i,k) )*tr_seri(i,k,it)

!
! Flux lessivage total 
! ------------------------------------------------------------
                 flestottr(i,k,it) = flestottr(i,k,it) -   &
                      ( d_tr_lessi_nucl(i,k,it)   +        &
                      d_tr_lessi_impa(i,k,it) ) *          &
                      ( paprs(i,k)-paprs(i,k+1) ) /        &
                      (RG * pdtphys)
!
! Mise a jour des traceurs due a l'impaction,nucleation 
! ----------------------------------------------------------------------
                 tr_seri(i,k,it)=tr_seri(i,k,it)*frac_impa(i,k)*frac_nucl(i,k)
              END DO
           END DO
        END IF
     END DO
     
  END IF ! lessivage

!=============================================================
!   Ecriture des sorties
!=============================================================
#ifdef CPP_IOIPSL
!  INCLUDE "write_histrac.h"
#endif

END SUBROUTINE phytrac
