!$Id $
!
MODULE tracinca_mod
!
! This module prepares and calls the INCA main subroutines. 
!

CONTAINS

  SUBROUTINE tracinca_init(aerosol,lessivage)
    ! This subroutine initialize some control varaibles. 

    USE infotrac
    IMPLICIT NONE
    
    ! Output variables
    LOGICAL,DIMENSION(nbtr), INTENT(OUT) :: aerosol
    LOGICAL,INTENT(OUT) :: lessivage
    
    
    ! Initialization
    lessivage  =.FALSE.
    aerosol(:) = .FALSE.
        
  END SUBROUTINE tracinca_init

  SUBROUTINE tracinca(                                &
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

!========================================================
!    -- CHIMIE INCA --
!========================================================

    USE dimphy
    USE infotrac
    USE vampir
    USE comgeomphy
    USE control_mod

    
    IMPLICIT NONE
    
    INCLUDE "indicesol.h"
    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"

!==========================================================================
!                   -- DESCRIPTION DES ARGUMENTS --
!==========================================================================


! EN ENTREE ...
!
!Configuration grille,temps:
    INTEGER,INTENT(IN) :: nstep      ! Appel physique
    INTEGER,INTENT(IN) :: julien     ! Jour julien
    REAL,INTENT(IN)    :: gmtime
    REAL,INTENT(IN)    :: pdtphys    ! Pas d'integration pour la physique (seconde)
    LOGICAL,INTENT(IN) :: lafin      ! le flag de la fin de la physique
    

!Physique: 
!--------
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: sh      ! humidite specifique
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: rh      ! humidite relative
    REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs   ! pression pour chaque inter-couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)   :: pphi    ! geopotentiel
    REAL,DIMENSION(klon),INTENT(IN)        :: pphis
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

! Arguments necessaires pour les sources et puits de traceur:
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol  ! Temperature du sol (surf)(Kelvin)
    REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: pctsrf ! Pourcentage de sol f(nature du sol)


  ! InOutput argument
    REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT) :: tr_seri ! Concentration Traceur [U/KgA]  

  ! Output arguments
    REAL,DIMENSION(klon,nbtr), INTENT(OUT)        :: source  ! a voir lorsque le flux de surface est prescrit 
    CHARACTER(len=8),DIMENSION(nbtr), INTENT(OUT) :: solsym

!=======================================================================================
!                        -- VARIABLES LOCALES TRACEURS --
!=======================================================================================

    INTEGER :: k
    REAL,DIMENSION(klon,klev) :: pdel
    REAL    :: calday
    INTEGER :: ncsec

    CALL VTe(VTphysiq)
    CALL VTb(VTinca)
    
    calday = REAL(julien) + gmtime
    ncsec  = NINT (86400.*gmtime)
     
    DO k = 1, klev
       pdel(:,k) = paprs(:,k) - paprs (:,k+1)
    END DO
    
    IF (config_inca == 'aero') THEN
#ifdef INCA
       CALL aerosolmain(                    &
            aerosol_couple,tr_seri,pdtphys, &
            pplay,pdel,prfl,pmflxr,psfl,    &
            pmflxs,pmfu,itop_con,ibas_con,  &
            pphi,airephy,nstep,rneb,t_seri, &      
            rh,tau_aero,piz_aero,cg_aero,   &
            rfname,ccm,lafin)
#endif
    END IF

#ifdef INCA
    CALL chemmain (tr_seri, &   !mmr
         nstep,      & !nstep
         calday,     & !calday
         julien,     & !ncdate
         ncsec,      & !ncsec
         1,          & !lat
         pdtphys,    & !delt
         paprs(1,1), & !ps
         pplay,      & !pmid
         pdel,       & !pdel
         airephy,    &
         pctsrf(1,1),& !oro
         ftsol,      & !tsurf
         albsol,     & !albs
         pphi,       & !zma
         pphis,      & !phis
         cldfra,     & !cldfr
         rneb,       & !cldfr_st
         diafra,     & !cldfr_cv
         itop_con,   & !cldtop
         ibas_con,   & !cldbot
         cldliq,     & !cwat
         prfl,       & !flxrst
         pmflxr,     & !flxrcv
         psfl,       & !flxsst
         pmflxs,     & !flxscv
         pmfu,       & !flxupd
         flxmass_w,  & !flxmass_w
         t_seri,     & !tfld
         sh,         & !sh
         rh,         & !rh
         iip1,       & !nx
         jjp1,       & !ny
         source,     &
         solsym)
#endif
    
    CALL VTe(VTinca)
    CALL VTb(VTphysiq)
    
    
  END SUBROUTINE tracinca


END MODULE tracinca_mod
