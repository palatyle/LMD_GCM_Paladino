module radlwsw_m

  IMPLICIT NONE

contains

SUBROUTINE radlwsw( &
   dist, rmu0, fract, &
   paprs, pplay,tsol,alb1, alb2, &
   t,q,wo,&
   cldfra, cldemi, cldtaupd,&
   ok_ade, ok_aie,&
   tau_aero, piz_aero, cg_aero,&
   cldtaupi, new_aod, &
   qsat, flwc, fiwc, &
   heat,heat0,cool,cool0,radsol,albpla,&
   topsw,toplw,solsw,sollw,&
   sollwdown,&
   topsw0,toplw0,solsw0,sollw0,&
   lwdn0, lwdn, lwup0, lwup,&
   swdn0, swdn, swup0, swup,&
   topswad_aero, solswad_aero,&
   topswai_aero, solswai_aero, &
   topswad0_aero, solswad0_aero,&
   topsw_aero, topsw0_aero,&
   solsw_aero, solsw0_aero, &
   topswcf_aero, solswcf_aero)



  USE DIMPHY
  use assert_m, only: assert

  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
  ! Objet: interface entre le modele et les rayonnements
  ! Arguments:
  ! dist-----input-R- distance astronomique terre-soleil
  ! rmu0-----input-R- cosinus de l'angle zenithal
  ! fract----input-R- duree d'ensoleillement normalisee
  ! co2_ppm--input-R- concentration du gaz carbonique (en ppm)
  ! paprs----input-R- pression a inter-couche (Pa)
  ! pplay----input-R- pression au milieu de couche (Pa)
  ! tsol-----input-R- temperature du sol (en K)
  ! alb1-----input-R- albedo du sol(entre 0 et 1) dans l'interval visible 
  ! alb2-----input-R- albedo du sol(entre 0 et 1) dans l'interval proche infra-rouge   
  ! t--------input-R- temperature (K)
  ! q--------input-R- vapeur d'eau (en kg/kg)
  ! cldfra---input-R- fraction nuageuse (entre 0 et 1)
  ! cldtaupd---input-R- epaisseur optique des nuages dans le visible (present-day value)
  ! cldemi---input-R- emissivite des nuages dans l'IR (entre 0 et 1)
  ! ok_ade---input-L- apply the Aerosol Direct Effect or not?
  ! ok_aie---input-L- apply the Aerosol Indirect Effect or not?
  ! tau_ae, piz_ae, cg_ae-input-R- aerosol optical properties (calculated in aeropt.F)
  ! cldtaupi-input-R- epaisseur optique des nuages dans le visible
  !                   calculated for pre-industrial (pi) aerosol concentrations, i.e. with smaller
  !                   droplet concentration, thus larger droplets, thus generally cdltaupi cldtaupd
  !                   it is needed for the diagnostics of the aerosol indirect radiative forcing      
  !
  ! heat-----output-R- echauffement atmospherique (visible) (K/jour)
  ! cool-----output-R- refroidissement dans l'IR (K/jour)
  ! radsol---output-R- bilan radiatif net au sol (W/m**2) (+ vers le bas)
  ! albpla---output-R- albedo planetaire (entre 0 et 1)
  ! topsw----output-R- flux solaire net au sommet de l'atm.
  ! toplw----output-R- ray. IR montant au sommet de l'atmosphere
  ! solsw----output-R- flux solaire net a la surface
  ! sollw----output-R- ray. IR montant a la surface
  ! solswad---output-R- ray. solaire net absorbe a la surface (aerosol dir)
  ! topswad---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol dir)
  ! solswai---output-R- ray. solaire net absorbe a la surface (aerosol ind)
  ! topswai---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol ind)
  !
  ! ATTENTION: swai and swad have to be interpreted in the following manner:
  ! ---------
  ! ok_ade=F & ok_aie=F -both are zero
  ! ok_ade=T & ok_aie=F -aerosol direct forcing is F_{AD} = topsw-topswad
  !                        indirect is zero
  ! ok_ade=F & ok_aie=T -aerosol indirect forcing is F_{AI} = topsw-topswai
  !                        direct is zero
  ! ok_ade=T & ok_aie=T -aerosol indirect forcing is F_{AI} = topsw-topswai
  !                        aerosol direct forcing is F_{AD} = topswai-topswad
  !
  
  !======================================================================
  
  ! ====================================================================
  ! Adapte au modele de chimie INCA par Celine Deandreis & Anne Cozic -- 2009
  ! 1 = ZERO    
  ! 2 = AER total    
  ! 3 = NAT    
  ! 4 = BC    
  ! 5 = SO4    
  ! 6 = POM    
  ! 7 = DUST    
  ! 8 = SS    
  ! 9 = NO3    
  ! 
  ! ====================================================================
  include "YOETHF.h"
  include "YOMCST.h"
  include "clesphys.h"
  include "iniprint.h"

! Input arguments
  REAL,    INTENT(in)  :: dist
  REAL,    INTENT(in)  :: rmu0(KLON), fract(KLON)
  REAL,    INTENT(in)  :: paprs(KLON,KLEV+1), pplay(KLON,KLEV)
  REAL,    INTENT(in)  :: alb1(KLON), alb2(KLON), tsol(KLON)
  REAL,    INTENT(in)  :: t(KLON,KLEV), q(KLON,KLEV)

  REAL, INTENT(in):: wo(:, :, :) ! dimension(KLON,KLEV, 1 or 2)
  ! column-density of ozone in a layer, in kilo-Dobsons
  ! "wo(:, :, 1)" is for the average day-night field, 
  ! "wo(:, :, 2)" is for daylight time.

  LOGICAL, INTENT(in)  :: ok_ade, ok_aie                                 ! switches whether to use aerosol direct (indirect) effects or not
  REAL,    INTENT(in)  :: cldfra(KLON,KLEV), cldemi(KLON,KLEV), cldtaupd(KLON,KLEV)
  REAL,    INTENT(in)  :: tau_aero(KLON,KLEV,9,2)                        ! aerosol optical properties (see aeropt.F)
  REAL,    INTENT(in)  :: piz_aero(KLON,KLEV,9,2)                        ! aerosol optical properties (see aeropt.F)
  REAL,    INTENT(in)  :: cg_aero(KLON,KLEV,9,2)                         ! aerosol optical properties (see aeropt.F)
  REAL,    INTENT(in)  :: cldtaupi(KLON,KLEV)                            ! cloud optical thickness for pre-industrial aerosol concentrations
  LOGICAL, INTENT(in)  :: new_aod                                        ! flag pour retrouver les resultats exacts de l'AR4 dans le cas ou l'on ne travaille qu'avec les sulfates
  REAL,    INTENT(in)  :: qsat(klon,klev) ! Variable pour iflag_rrtm=1
  REAL,    INTENT(in)  :: flwc(klon,klev) ! Variable pour iflag_rrtm=1
  REAL,    INTENT(in)  :: fiwc(klon,klev) ! Variable pour iflag_rrtm=1

! Output arguments
  REAL,    INTENT(out) :: heat(KLON,KLEV), cool(KLON,KLEV)
  REAL,    INTENT(out) :: heat0(KLON,KLEV), cool0(KLON,KLEV)
  REAL,    INTENT(out) :: radsol(KLON), topsw(KLON), toplw(KLON)
  REAL,    INTENT(out) :: solsw(KLON), sollw(KLON), albpla(KLON)
  REAL,    INTENT(out) :: topsw0(KLON), toplw0(KLON), solsw0(KLON), sollw0(KLON)
  REAL,    INTENT(out) :: sollwdown(KLON)
  REAL,    INTENT(out) :: swdn(KLON,kflev+1),swdn0(KLON,kflev+1)
  REAL,    INTENT(out) :: swup(KLON,kflev+1),swup0(KLON,kflev+1)
  REAL,    INTENT(out) :: lwdn(KLON,kflev+1),lwdn0(KLON,kflev+1)
  REAL,    INTENT(out) :: lwup(KLON,kflev+1),lwup0(KLON,kflev+1)
  REAL,    INTENT(out) :: topswad_aero(KLON), solswad_aero(KLON)         ! output: aerosol direct forcing at TOA and surface
  REAL,    INTENT(out) :: topswai_aero(KLON), solswai_aero(KLON)         ! output: aerosol indirect forcing atTOA and surface
  REAL, DIMENSION(klon), INTENT(out)    :: topswad0_aero 
  REAL, DIMENSION(klon), INTENT(out)    :: solswad0_aero
  REAL, DIMENSION(kdlon,9), INTENT(out) :: topsw_aero
  REAL, DIMENSION(kdlon,9), INTENT(out) :: topsw0_aero
  REAL, DIMENSION(kdlon,9), INTENT(out) :: solsw_aero
  REAL, DIMENSION(kdlon,9), INTENT(out) :: solsw0_aero
  REAL, DIMENSION(kdlon,3), INTENT(out) :: topswcf_aero
  REAL, DIMENSION(kdlon,3), INTENT(out) :: solswcf_aero

! Local variables
  REAL(KIND=8) ZFSUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSUP0(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN0(KDLON,KFLEV+1)
  REAL(KIND=8) ZFLUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZFLDN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFLUP0(KDLON,KFLEV+1)
  REAL(KIND=8) ZFLDN0(KDLON,KFLEV+1)
  REAL(KIND=8) zx_alpha1, zx_alpha2
  INTEGER k, kk, i, j, iof, nb_gr
  REAL(KIND=8) PSCT
  REAL(KIND=8) PALBD(kdlon,2), PALBP(kdlon,2)
  REAL(KIND=8) PEMIS(kdlon), PDT0(kdlon), PVIEW(kdlon)
  REAL(KIND=8) PPSOL(kdlon), PDP(kdlon,KLEV)
  REAL(KIND=8) PTL(kdlon,kflev+1), PPMB(kdlon,kflev+1)
  REAL(KIND=8) PTAVE(kdlon,kflev)
  REAL(KIND=8) PWV(kdlon,kflev), PQS(kdlon,kflev)

  real(kind=8) POZON(kdlon, kflev, size(wo, 3)) ! mass fraction of ozone
  ! "POZON(:, :, 1)" is for the average day-night field, 
  ! "POZON(:, :, 2)" is for daylight time.

  REAL(KIND=8) PAER(kdlon,kflev,5)
  REAL(KIND=8) PCLDLD(kdlon,kflev)
  REAL(KIND=8) PCLDLU(kdlon,kflev)
  REAL(KIND=8) PCLDSW(kdlon,kflev)
  REAL(KIND=8) PTAU(kdlon,2,kflev)
  REAL(KIND=8) POMEGA(kdlon,2,kflev)
  REAL(KIND=8) PCG(kdlon,2,kflev)
  REAL(KIND=8) zfract(kdlon), zrmu0(kdlon), zdist
  REAL(KIND=8) zheat(kdlon,kflev), zcool(kdlon,kflev)
  REAL(KIND=8) zheat0(kdlon,kflev), zcool0(kdlon,kflev)
  REAL(KIND=8) ztopsw(kdlon), ztoplw(kdlon)
  REAL(KIND=8) zsolsw(kdlon), zsollw(kdlon), zalbpla(kdlon)
  REAL(KIND=8) zsollwdown(kdlon)
  REAL(KIND=8) ztopsw0(kdlon), ztoplw0(kdlon)
  REAL(KIND=8) zsolsw0(kdlon), zsollw0(kdlon)
  REAL(KIND=8) zznormcp
  REAL(KIND=8) tauaero(kdlon,kflev,9,2)                     ! aer opt properties
  REAL(KIND=8) pizaero(kdlon,kflev,9,2)
  REAL(KIND=8) cgaero(kdlon,kflev,9,2)
  REAL(KIND=8) PTAUA(kdlon,2,kflev)                         ! present-day value of cloud opt thickness (PTAU is pre-industrial value), local use
  REAL(KIND=8) POMEGAA(kdlon,2,kflev)                       ! dito for single scatt albedo
  REAL(KIND=8) ztopswadaero(kdlon), zsolswadaero(kdlon)     ! Aerosol direct forcing at TOAand surface
  REAL(KIND=8) ztopswad0aero(kdlon), zsolswad0aero(kdlon)   ! Aerosol direct forcing at TOAand surface
  REAL(KIND=8) ztopswaiaero(kdlon), zsolswaiaero(kdlon)     ! dito, indirect
  REAL(KIND=8) ztopsw_aero(kdlon,9), ztopsw0_aero(kdlon,9)
  REAL(KIND=8) zsolsw_aero(kdlon,9), zsolsw0_aero(kdlon,9)
  REAL(KIND=8) ztopswcf_aero(kdlon,3), zsolswcf_aero(kdlon,3)     
  real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

  call assert(size(wo, 1) == klon, size(wo, 2) == klev, "radlwsw wo")
  ! initialisation
  tauaero(:,:,:,:)=0.
  pizaero(:,:,:,:)=0.
  cgaero(:,:,:,:)=0.
  
  !
  !-------------------------------------------
  nb_gr = KLON / kdlon
  IF (nb_gr*kdlon .NE. KLON) THEN
      PRINT*, "kdlon mauvais:", KLON, kdlon, nb_gr
      CALL abort
  ENDIF
  IF (kflev .NE. KLEV) THEN
      PRINT*, "kflev differe de KLEV, kflev, KLEV"
      CALL abort
  ENDIF
  !-------------------------------------------
  DO k = 1, KLEV
    DO i = 1, KLON
      heat(i,k)=0.
      cool(i,k)=0.
      heat0(i,k)=0.
      cool0(i,k)=0.
    ENDDO
  ENDDO
  !
  zdist = dist
  !
  PSCT = solaire/zdist/zdist
  DO j = 1, nb_gr
    iof = kdlon*(j-1)
    DO i = 1, kdlon
      zfract(i) = fract(iof+i)
      zrmu0(i) = rmu0(iof+i)
      PALBD(i,1) = alb1(iof+i)
      PALBD(i,2) = alb2(iof+i)
      PALBP(i,1) = alb1(iof+i)
      PALBP(i,2) = alb2(iof+i)
      PEMIS(i) = 1.0 
      PVIEW(i) = 1.66
      PPSOL(i) = paprs(iof+i,1)
      zx_alpha1 = (paprs(iof+i,1)-pplay(iof+i,2))/(pplay(iof+i,1)-pplay(iof+i,2))
      zx_alpha2 = 1.0 - zx_alpha1
      PTL(i,1) = t(iof+i,1) * zx_alpha1 + t(iof+i,2) * zx_alpha2
      PTL(i,KLEV+1) = t(iof+i,KLEV)
      PDT0(i) = tsol(iof+i) - PTL(i,1)
    ENDDO
    DO k = 2, kflev
      DO i = 1, kdlon
        PTL(i,k) = (t(iof+i,k)+t(iof+i,k-1))*0.5
      ENDDO
    ENDDO
    DO k = 1, kflev
      DO i = 1, kdlon
        PDP(i,k) = paprs(iof+i,k)-paprs(iof+i,k+1)
        PTAVE(i,k) = t(iof+i,k)
        PWV(i,k) = MAX (q(iof+i,k), 1.0e-12)
        PQS(i,k) = PWV(i,k)
        POZON(i,k, :) = wo(iof+i, k, :) * RG * dobson_u * 1e3 &
             / (paprs(iof+i, k) - paprs(iof+i, k+1))
        PCLDLD(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
        PCLDLU(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
        PCLDSW(i,k) = cldfra(iof+i,k)
        PTAU(i,1,k) = MAX(cldtaupi(iof+i,k), 1.0e-05)! 1e-12 serait instable
        PTAU(i,2,k) = MAX(cldtaupi(iof+i,k), 1.0e-05)! pour 32-bit machines
        POMEGA(i,1,k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAU(i,1,k))
        POMEGA(i,2,k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i,2,k))
        PCG(i,1,k) = 0.865
        PCG(i,2,k) = 0.910
        !-
        ! Introduced for aerosol indirect forcings.
        ! The following values use the cloud optical thickness calculated from
        ! present-day aerosol concentrations whereas the quantities without the
        ! "A" at the end are for pre-industial (natural-only) aerosol concentrations
        !
        PTAUA(i,1,k) = MAX(cldtaupd(iof+i,k), 1.0e-05)! 1e-12 serait instable
        PTAUA(i,2,k) = MAX(cldtaupd(iof+i,k), 1.0e-05)! pour 32-bit machines
        POMEGAA(i,1,k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAUA(i,1,k))
        POMEGAA(i,2,k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAUA(i,2,k))
      ENDDO
    ENDDO
    !
    DO k = 1, kflev+1
      DO i = 1, kdlon
        PPMB(i,k) = paprs(iof+i,k)/100.0
      ENDDO
    ENDDO
    !
    DO kk = 1, 5
      DO k = 1, kflev
        DO i = 1, kdlon
          PAER(i,k,kk) = 1.0E-15
        ENDDO
      ENDDO
    ENDDO
    DO k = 1, kflev
      DO i = 1, kdlon
        tauaero(i,k,:,1)=tau_aero(iof+i,k,:,1)
        pizaero(i,k,:,1)=piz_aero(iof+i,k,:,1)
        cgaero(i,k,:,1) =cg_aero(iof+i,k,:,1)
        tauaero(i,k,:,2)=tau_aero(iof+i,k,:,2)
        pizaero(i,k,:,2)=piz_aero(iof+i,k,:,2)
        cgaero(i,k,:,2) =cg_aero(iof+i,k,:,2)
      ENDDO
    ENDDO

!
!===== iflag_rrtm ================================================
!      
    IF (iflag_rrtm == 0) THEN
       ! Old radiation scheme, used for AR4 runs
       ! average day-night ozone for longwave
       CALL LW_LMDAR4(&
            PPMB, PDP,&
            PPSOL,PDT0,PEMIS,&
            PTL, PTAVE, PWV, POZON(:, :, 1), PAER,&
            PCLDLD,PCLDLU,&
            PVIEW,&
            zcool, zcool0,&
            ztoplw,zsollw,ztoplw0,zsollw0,&
            zsollwdown,&
            ZFLUP, ZFLDN, ZFLUP0,ZFLDN0)

       ! daylight ozone, if we have it, for short wave
       IF (.NOT. new_aod) THEN 
          ! use old version
          CALL SW_LMDAR4(PSCT, zrmu0, zfract,&
               PPMB, PDP, &
               PPSOL, PALBD, PALBP,&
               PTAVE, PWV, PQS, POZON(:, :, size(wo, 3)), PAER,&
               PCLDSW, PTAU, POMEGA, PCG,&
               zheat, zheat0,&
               zalbpla,ztopsw,zsolsw,ztopsw0,zsolsw0,&
               ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,&
               tau_aero(:,:,5,:), piz_aero(:,:,5,:), cg_aero(:,:,5,:),& 
               PTAUA, POMEGAA,&
               ztopswadaero,zsolswadaero,&
               ztopswaiaero,zsolswaiaero,& 
               ok_ade, ok_aie) 
          
       ELSE ! new_aod=T         
          CALL SW_AEROAR4(PSCT, zrmu0, zfract,&
               PPMB, PDP,&
               PPSOL, PALBD, PALBP,&
               PTAVE, PWV, PQS, POZON(:, :, size(wo, 3)), PAER,&
               PCLDSW, PTAU, POMEGA, PCG,&
               zheat, zheat0,&
               zalbpla,ztopsw,zsolsw,ztopsw0,zsolsw0,&
               ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,&
               tauaero, pizaero, cgaero, &
               PTAUA, POMEGAA,&
               ztopswadaero,zsolswadaero,&
               ztopswad0aero,zsolswad0aero,&
               ztopswaiaero,zsolswaiaero, & 
               ztopsw_aero,ztopsw0_aero,&
               zsolsw_aero,zsolsw0_aero,&
               ztopswcf_aero,zsolswcf_aero, & 
               ok_ade, ok_aie) 
          
       ENDIF

    ELSE  
!===== iflag_rrtm=1, on passe dans SW via RECMWFL ===============
       WRITE(lunout,*) "Option iflag_rrtm=T ne fonctionne pas encore !!!"
       CALL abort_gcm('radlwsw','iflag_rrtm=T not valid',1) 

    ENDIF ! iflag_rrtm
!======================================================================

    DO i = 1, kdlon
      radsol(iof+i) = zsolsw(i) + zsollw(i)
      topsw(iof+i) = ztopsw(i)
      toplw(iof+i) = ztoplw(i)
      solsw(iof+i) = zsolsw(i)
      sollw(iof+i) = zsollw(i)
      sollwdown(iof+i) = zsollwdown(i)
      DO k = 1, kflev+1
        lwdn0 ( iof+i,k)   = ZFLDN0 ( i,k)
        lwdn  ( iof+i,k)   = ZFLDN  ( i,k)
        lwup0 ( iof+i,k)   = ZFLUP0 ( i,k)
        lwup  ( iof+i,k)   = ZFLUP  ( i,k)
      ENDDO
      topsw0(iof+i) = ztopsw0(i)
      toplw0(iof+i) = ztoplw0(i)
      solsw0(iof+i) = zsolsw0(i)
      sollw0(iof+i) = zsollw0(i)
      albpla(iof+i) = zalbpla(i)

      DO k = 1, kflev+1
        swdn0 ( iof+i,k)   = ZFSDN0 ( i,k)
        swdn  ( iof+i,k)   = ZFSDN  ( i,k)
        swup0 ( iof+i,k)   = ZFSUP0 ( i,k)
        swup  ( iof+i,k)   = ZFSUP  ( i,k)
      ENDDO
    ENDDO
    !-transform the aerosol forcings, if they have
    ! to be calculated
    IF (ok_ade) THEN
        DO i = 1, kdlon
          topswad_aero(iof+i) = ztopswadaero(i)
          topswad0_aero(iof+i) = ztopswad0aero(i)
          solswad_aero(iof+i) = zsolswadaero(i)
          solswad0_aero(iof+i) = zsolswad0aero(i)
! MS the following lines seem to be wrong, why is iof on right hand side???
!          topsw_aero(iof+i,:) = ztopsw_aero(iof+i,:)
!          topsw0_aero(iof+i,:) = ztopsw0_aero(iof+i,:)
!          solsw_aero(iof+i,:) = zsolsw_aero(iof+i,:)
!          solsw0_aero(iof+i,:) = zsolsw0_aero(iof+i,:)
          topsw_aero(iof+i,:) = ztopsw_aero(i,:)
          topsw0_aero(iof+i,:) = ztopsw0_aero(i,:)
          solsw_aero(iof+i,:) = zsolsw_aero(i,:)
          solsw0_aero(iof+i,:) = zsolsw0_aero(i,:)
          topswcf_aero(iof+i,:) = ztopswcf_aero(i,:)
          solswcf_aero(iof+i,:) = zsolswcf_aero(i,:)          
        ENDDO
    ELSE
        DO i = 1, kdlon
          topswad_aero(iof+i) = 0.0
          solswad_aero(iof+i) = 0.0
          topswad0_aero(iof+i) = 0.0
          solswad0_aero(iof+i) = 0.0
          topsw_aero(iof+i,:) = 0.
          topsw0_aero(iof+i,:) =0.
          solsw_aero(iof+i,:) = 0.
          solsw0_aero(iof+i,:) = 0.
        ENDDO
    ENDIF
    IF (ok_aie) THEN
        DO i = 1, kdlon
          topswai_aero(iof+i) = ztopswaiaero(i)
          solswai_aero(iof+i) = zsolswaiaero(i)
        ENDDO
    ELSE
        DO i = 1, kdlon
          topswai_aero(iof+i) = 0.0
          solswai_aero(iof+i) = 0.0
        ENDDO
    ENDIF
    DO k = 1, kflev
      DO i = 1, kdlon
        !        scale factor to take into account the difference between
        !        dry air and watter vapour scpecifi! heat capacity
        zznormcp=1.0+RVTMP2*PWV(i,k)
        heat(iof+i,k) = zheat(i,k)/zznormcp
        cool(iof+i,k) = zcool(i,k)/zznormcp
        heat0(iof+i,k) = zheat0(i,k)/zznormcp
        cool0(iof+i,k) = zcool0(i,k)/zznormcp
      ENDDO
    ENDDO

 ENDDO ! j = 1, nb_gr

END SUBROUTINE radlwsw

end module radlwsw_m
