!
! $Id$
!
SUBROUTINE SW_AEROAR4(PSCT, PRMU0, PFRAC, &
     PPMB, PDP, &
     PPSOL, PALBD, PALBP,&
     PTAVE, PWV, PQS, POZON, PAER,&
     PCLDSW, PTAU, POMEGA, PCG,&
     PHEAT, PHEAT0,&
     PALBPLA,PTOPSW,PSOLSW,PTOPSW0,PSOLSW0,&
     ZFSUP,ZFSDN,ZFSUP0,ZFSDN0,&
     tauaero, pizaero, cgaero,&
     PTAUA, POMEGAA,&
     PTOPSWADAERO,PSOLSWADAERO,&
     PTOPSWAD0AERO,PSOLSWAD0AERO,&
     PTOPSWAIAERO,PSOLSWAIAERO,&
     PTOPSWAERO,PTOPSW0AERO,&
     PSOLSWAERO,PSOLSW0AERO,&
     PTOPSWCFAERO,PSOLSWCFAERO,&
     ok_ade, ok_aie )

  USE dimphy

  IMPLICIT NONE

#include "YOMCST.h"
#include "clesphys.h"
  !
  !     ------------------------------------------------------------------
  !
  !     PURPOSE.
  !     --------
  !
  !          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
  !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
  !
  !     METHOD.
  !     -------
  !
  !          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
  !          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
  !          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)
  !
  !     REFERENCE.
  !     ----------
  !
  !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
  !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
  !
  !     AUTHOR.
  !     -------
  !        JEAN-JACQUES MORCRETTE  *ECMWF*
  !
  !     MODIFICATIONS.
  !     --------------
  !        ORIGINAL : 89-07-14
  !        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
  !        03-11-27   J. QUAAS Introduce aerosol forcings (based on BOUCHER)
  !        09-04      A. COZIC - C.DEANDREIS Indroduce NAT/BC/POM/DUST/SS aerosol forcing
  !     ------------------------------------------------------------------
  !
  !* ARGUMENTS:
  !
  REAL(KIND=8) PSCT  ! constante solaire (valeur conseillee: 1370)

  REAL(KIND=8) PPSOL(KDLON)        ! SURFACE PRESSURE (PA)
  REAL(KIND=8) PDP(KDLON,KFLEV)    ! LAYER THICKNESS (PA)
  REAL(KIND=8) PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)

  REAL(KIND=8) PRMU0(KDLON)  ! COSINE OF ZENITHAL ANGLE
  REAL(KIND=8) PFRAC(KDLON)  ! fraction de la journee

  REAL(KIND=8) PTAVE(KDLON,KFLEV)  ! LAYER TEMPERATURE (K)
  REAL(KIND=8) PWV(KDLON,KFLEV)    ! SPECIFI! HUMIDITY (KG/KG)
  REAL(KIND=8) PQS(KDLON,KFLEV)    ! SATURATED WATER VAPOUR (KG/KG)
  REAL(KIND=8) POZON(KDLON,KFLEV)  ! OZONE CONCENTRATION (KG/KG)
  REAL(KIND=8) PAER(KDLON,KFLEV,5) ! AEROSOLS' OPTICAL THICKNESS

  REAL(KIND=8) PALBD(KDLON,2)  ! albedo du sol (lumiere diffuse)
  REAL(KIND=8) PALBP(KDLON,2)  ! albedo du sol (lumiere parallele)

  REAL(KIND=8) PCLDSW(KDLON,KFLEV)    ! CLOUD FRACTION
  REAL(KIND=8) PTAU(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS
  REAL(KIND=8) PCG(KDLON,2,KFLEV)     ! ASYMETRY FACTOR
  REAL(KIND=8) POMEGA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO

  REAL(KIND=8) PHEAT(KDLON,KFLEV) ! SHORTWAVE HEATING (K/DAY)
  REAL(KIND=8) PHEAT0(KDLON,KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
  REAL(KIND=8) PALBPLA(KDLON)     ! PLANETARY ALBEDO
  REAL(KIND=8) PTOPSW(KDLON)      ! SHORTWAVE FLUX AT T.O.A.
  REAL(KIND=8) PSOLSW(KDLON)      ! SHORTWAVE FLUX AT SURFACE
  REAL(KIND=8) PTOPSW0(KDLON)     ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
  REAL(KIND=8) PSOLSW0(KDLON)     ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
  !
  !* LOCAL VARIABLES:
  !
  real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

  REAL(KIND=8) ZOZ(KDLON,KFLEV)
  ! column-density of ozone in layer, in kilo-Dobsons

  REAL(KIND=8) ZAKI(KDLON,2)     
  REAL(KIND=8) ZCLD(KDLON,KFLEV)
  REAL(KIND=8) ZCLEAR(KDLON) 
  REAL(KIND=8) ZDSIG(KDLON,KFLEV)
  REAL(KIND=8) ZFACT(KDLON)
  REAL(KIND=8) ZFD(KDLON,KFLEV+1)
  REAL(KIND=8) ZFDOWN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFU(KDLON,KFLEV+1)
  REAL(KIND=8) ZFUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZRMU(KDLON)
  REAL(KIND=8) ZSEC(KDLON)
  REAL(KIND=8) ZUD(KDLON,5,KFLEV+1)
  REAL(KIND=8) ZCLDSW0(KDLON,KFLEV)

  REAL(KIND=8) ZFSUP(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSUP0(KDLON,KFLEV+1)
  REAL(KIND=8) ZFSDN0(KDLON,KFLEV+1)

  INTEGER inu, jl, jk, i, k, kpl1

  INTEGER swpas  ! Every swpas steps, sw is calculated
  PARAMETER(swpas=1)

  INTEGER, SAVE :: itapsw = 0
  !$OMP THREADPRIVATE(itapsw)
  LOGICAL, SAVE :: appel1er = .TRUE.
  !$OMP THREADPRIVATE(appel1er)
  LOGICAL, SAVE :: initialized = .FALSE.
  !$OMP THREADPRIVATE(initialized)

  !jq-Introduced for aerosol forcings
  REAL(KIND=8), SAVE :: flag_aer
  !$OMP THREADPRIVATE(flag_aer)

  LOGICAL ok_ade, ok_aie    ! use aerosol forcings or not?
  REAL(KIND=8) tauaero(kdlon,kflev,9,2)  ! aerosol optical properties
  REAL(KIND=8) pizaero(kdlon,kflev,9,2)  ! (see aeropt.F)
  REAL(KIND=8) cgaero(kdlon,kflev,9,2)   ! -"-
  REAL(KIND=8) PTAUA(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS (pre-industrial value)
  REAL(KIND=8) POMEGAA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
  REAL(KIND=8) PTOPSWADAERO(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
  REAL(KIND=8) PSOLSWADAERO(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
  REAL(KIND=8) PTOPSWAD0AERO(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL DIR)
  REAL(KIND=8) PSOLSWAD0AERO(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL DIR)
  REAL(KIND=8) PTOPSWAIAERO(KDLON)     ! SHORTWAVE FLUX AT T.O.A.(+AEROSOL IND)
  REAL(KIND=8) PSOLSWAIAERO(KDLON)     ! SHORTWAVE FLUX AT SURFACE(+AEROSOL IND)
  REAL(KIND=8) PTOPSWAERO(KDLON,9)	 ! SW TOA AS DRF nat & ant 
  REAL(KIND=8) PTOPSW0AERO(KDLON,9)	 ! SW SRF AS DRF nat & ant 
  REAL(KIND=8) PSOLSWAERO(KDLON,9)	 ! SW TOA CS DRF nat & ant
  REAL(KIND=8) PSOLSW0AERO(KDLON,9)	 ! SW SRF CS DRF nat & ant
  REAL(KIND=8) PTOPSWCFAERO(KDLON,3)   !  SW TOA AS cloudRF nat & ant 
  REAL(KIND=8) PSOLSWCFAERO(KDLON,3)   !  SW SRF AS cloudRF nat & ant 

  !jq - Fluxes including aerosol effects
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAD_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAD_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAD_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAD_AERO)
  !jq - Fluxes including aerosol effects
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAD0_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAD0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAD0_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAD0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSUPAI_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSUPAI_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE :: ZFSDNAI_AERO(:,:)
  !$OMP THREADPRIVATE(ZFSDNAI_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSUP_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSUP_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSDN_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSDN_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSUP0_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSUP0_AERO)
  REAL(KIND=8),ALLOCATABLE,SAVE ::  ZFSDN0_AERO(:,:,:)
  !$OMP THREADPRIVATE(ZFSDN0_AERO)

! Key to define the aerosol effect acting on climate
! 0: aerosol feedback active according to ok_ade, ok_aie  DEFAULT 
! 1: no feedback , zero aerosol fluxes are used for climate, diagnostics according to ok_ade_ok_aie
! 2: feedback according to total aerosol direct effect used for climate, diagnostics according to ok_ade, ok_aie
! 3: feedback according to natural aerosol direct effect used for climate, diagnostics according to ok_ade_ok_aie

  INTEGER,SAVE :: AEROSOLFEEDBACK_ACTIVE = 0
!$OMP THREADPRIVATE(AEROSOLFEEDBACK_ACTIVE)  

      CHARACTER (LEN=20) :: modname='sw_aeroAR4'
      CHARACTER (LEN=80) :: abort_message

  IF ((.not. ok_ade) .and. (AEROSOLFEEDBACK_ACTIVE .ge. 2)) THEN
     abort_message ='Error: direct effect is not activated but assumed to be active - see sw_aeroAR4.F90'
     CALL abort_gcm (modname,abort_message,1)
  ENDIF
  AEROSOLFEEDBACK_ACTIVE=MIN(MAX(AEROSOLFEEDBACK_ACTIVE,0),3)
  IF  (AEROSOLFEEDBACK_ACTIVE .gt. 3) THEN
     abort_message ='Error: AEROSOLFEEDBACK_ACTIVE options go only until 3'
     CALL abort_gcm (modname,abort_message,1)
  ENDIF

  IF(.NOT.initialized) THEN
     flag_aer=0.
     initialized=.TRUE.
     ALLOCATE(ZFSUPAD_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAD_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSUPAD0_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAD0_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSUPAI_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSDNAI_AERO(KDLON,KFLEV+1))
     ALLOCATE(ZFSUP_AERO (KDLON,KFLEV+1,9))
     ALLOCATE(ZFSDN_AERO (KDLON,KFLEV+1,9))
     ALLOCATE(ZFSUP0_AERO(KDLON,KFLEV+1,9))
     ALLOCATE(ZFSDN0_AERO(KDLON,KFLEV+1,9))
     ZFSUPAD_AERO(:,:)=0.
     ZFSDNAD_AERO(:,:)=0.
     ZFSUPAD0_AERO(:,:)=0.
     ZFSDNAD0_AERO(:,:)=0.
     ZFSUPAI_AERO(:,:)=0.
     ZFSDNAI_AERO(:,:)=0.
     ZFSUP_AERO (:,:,:)=0.
     ZFSDN_AERO (:,:,:)=0.
     ZFSUP0_AERO(:,:,:)=0.
     ZFSDN0_AERO(:,:,:)=0.
  ENDIF

  IF (appel1er) THEN
     PRINT*, 'SW calling frequency : ', swpas
     PRINT*, "   In general, it should be 1"
     appel1er = .FALSE.
  ENDIF
  !     ------------------------------------------------------------------
  IF (MOD(itapsw,swpas).EQ.0) THEN

     DO JK = 1 , KFLEV
        DO JL = 1, KDLON
           ZCLDSW0(JL,JK) = 0.0
           ZOZ(JL,JK) = POZON(JL,JK)*46.6968/RG &
                *PDP(JL,JK)*(101325.0/PPSOL(JL))
        ENDDO
     ENDDO

! clear sky is either computed IF no direct effect is asked for, or for extended diag) 
     IF (( lev_histmth .ge. 4 ) .or. ( .not. ok_ade )) THEN    

     ! clear-sky: zero aerosol effect
     flag_aer=0.0
     CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
          PRMU0,PFRAC,PTAVE,PWV,&
          ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
     INU = 1
     CALL SW1S_LMDAR4(INU,PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          ZFD, ZFU)
     INU = 2
     CALL SW2S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          PWV, PQS,&
          ZFDOWN, ZFUP)
     DO JK = 1 , KFLEV+1
        DO JL = 1, KDLON
           ZFSUP0_AERO(JL,JK,1) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
           ZFSDN0_AERO(JL,JK,1) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
        ENDDO
     ENDDO
     ENDIF

! cloudy sky is either computed IF no indirect effect is asked for, or for extended diag) 
     IF (( lev_histmth .ge. 4 ) .or. ( .not. ok_aie )) THEN    
     ! cloudy-sky: zero aerosol effect
     flag_aer=0.0
     CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
          PRMU0,PFRAC,PTAVE,PWV,&
          ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
     INU = 1
     CALL SW1S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          ZFD, ZFU)
     INU = 2
     CALL SW2S_LMDAR4(INU, PAER, flag_aer, &
          tauaero(:,:,1,:), pizaero(:,:,1,:), cgaero(:,:,1,:),&
          ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
          ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
          PWV, PQS,&
          ZFDOWN, ZFUP)

     DO JK = 1 , KFLEV+1
        DO JL = 1, KDLON
           ZFSUP_AERO(JL,JK,1) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
           ZFSDN_AERO(JL,JK,1) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
        ENDDO
     ENDDO
     ENDIF


     IF (ok_ade) THEN

        ! clear sky (Anne Cozic 03/07/2007) direct effect of total aerosol
        ! CAS AER (2)
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP0_AERO(JL,JK,2) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL) 
              ZFSDN0_AERO(JL,JK,2) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO

! cloudy sky is either computed IF no indirect effect is asked for, or for extended diag) 
        IF (( lev_histmth .ge. 2 ) .or. (.not. ok_aie)) THEN  
        ! cloudy-sky aerosol direct effect of total aerosol
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,2) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL) 
              ZFSDN_AERO(JL,JK,2) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO
        ENDIF

! natural aeroosl clear sky is  computed  for extended diag) 
        IF ( lev_histmth .ge. 4 ) THEN            
        ! clear sky direct effect natural aerosol
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,ZCLDSW0,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP0_AERO(JL,JK,3) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN0_AERO(JL,JK,3) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
           ENDDO
        ENDDO
        ENDIF

! cloud sky natural is for extended diagnostics
        IF ( lev_histmth .ge. 2 ) THEN
        ! cloudy-sky direct effect natural aerosol
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,3,:), pizaero(:,:,3,:), cgaero(:,:,3,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)

        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,3) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN_AERO(JL,JK,3) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
           ENDDO
        ENDDO
        ENDIF

     ENDIF ! ok_ade

! cloudy sky needs to be computed in all cases IF ok_aie is activated
     IF (ok_aie) THEN
        !jq   cloudy-sky + aerosol direct + aerosol indirect of total aerosol
        flag_aer=1.0
        CALL SWU_LMDAR4(PSCT,PCLDSW,PPMB,PPSOL,&
             PRMU0,PFRAC,PTAVE,PWV,&
             ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
        INU = 1
        CALL SW1S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             ZFD, ZFU)
        INU = 2
        CALL SW2S_LMDAR4(INU, PAER, flag_aer,&
             tauaero(:,:,2,:), pizaero(:,:,2,:), cgaero(:,:,2,:),&
             ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,&
             ZDSIG, POMEGAA, ZOZ, ZRMU, ZSEC, PTAUA, ZUD,&
             PWV, PQS,&
             ZFDOWN, ZFUP)
        DO JK = 1 , KFLEV+1
           DO JL = 1, KDLON
              ZFSUP_AERO(JL,JK,4) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
              ZFSDN_AERO(JL,JK,4) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL) 
           ENDDO
        ENDDO
     ENDIF ! ok_aie      

     itapsw = 0
  ENDIF
  itapsw = itapsw + 1

  IF  ( AEROSOLFEEDBACK_ACTIVE .eq. 0) THEN
  IF ( ok_ade .and. ok_aie  ) THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,4)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,4)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,2)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,2)
  ENDIF
  IF ( ok_ade .and. (.not. ok_aie) )  THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,2)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,2)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,2)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,2)
  ENDIF

  IF ( (.not. ok_ade) .and. ok_aie  )  THEN
    print*,'Warning: indirect effect in cloudy regions includes direct aerosol effect'
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,4)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,4)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,1)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,1)
  ENDIF
  IF ((.not. ok_ade) .and. (.not. ok_aie)) THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,1)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,1)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,1)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,1)
  ENDIF

! MS the following allows to compute the forcing diagostics without
! letting the aerosol forcing act on the meteorology
! SEE logic above
  ELSEIF  ( AEROSOLFEEDBACK_ACTIVE .gt. 0) THEN
    ZFSUP(:,:) =    ZFSUP_AERO(:,:,AEROSOLFEEDBACK_ACTIVE)
    ZFSDN(:,:) =    ZFSDN_AERO(:,:,AEROSOLFEEDBACK_ACTIVE)
    ZFSUP0(:,:) =   ZFSUP0_AERO(:,:,AEROSOLFEEDBACK_ACTIVE)
    ZFSDN0(:,:) =   ZFSDN0_AERO(:,:,AEROSOLFEEDBACK_ACTIVE)
  ENDIF
  

  DO k = 1, KFLEV
     kpl1 = k+1
     DO i = 1, KDLON
        PHEAT(i,k) = -(ZFSUP(i,kpl1)-ZFSUP(i,k))-(ZFSDN(i,k)-ZFSDN(i,kpl1))
        PHEAT(i,k) = PHEAT(i,k) * RDAY*RG/RCPD / PDP(i,k)
        PHEAT0(i,k) = -(ZFSUP0(i,kpl1)-ZFSUP0(i,k))-(ZFSDN0(i,k)-ZFSDN0(i,kpl1))
        PHEAT0(i,k) = PHEAT0(i,k) * RDAY*RG/RCPD / PDP(i,k)
     ENDDO
  ENDDO

  DO i = 1, KDLON
! effective SW surface albedo calculation
     PALBPLA(i) = ZFSUP(i,KFLEV+1)/(ZFSDN(i,KFLEV+1)+1.0e-20)
     
! clear sky net fluxes at TOA and SRF
     PSOLSW0(i) = ZFSDN0(i,1) - ZFSUP0(i,1)
     PTOPSW0(i) = ZFSDN0(i,KFLEV+1) - ZFSUP0(i,KFLEV+1)

! cloudy sky net fluxes at TOA and SRF
     PSOLSW(i) = ZFSDN(i,1) - ZFSUP(i,1)
     PTOPSW(i) = ZFSDN(i,KFLEV+1) - ZFSUP(i,KFLEV+1)


! net anthropogenic forcing direct and 1st indirect effect diagnostics
! requires a natural aerosol field read and used 
! Difference of net fluxes from double call to radiation


IF (ok_ade) THEN

! indices 1: natural; 2 anthropogenic 
! TOA/SRF all sky natural forcing
     PSOLSWAERO(i,1) = (ZFSDN_AERO(i,1,3) - ZFSUP_AERO(i,1,3))-(ZFSDN_AERO(i,1,1) - ZFSUP_AERO(i,1,1))
     PTOPSWAERO(i,1) = (ZFSDN_AERO(i,KFLEV+1,3) - ZFSUP_AERO(i,KFLEV+1,3))- (ZFSDN_AERO(i,KFLEV+1,1) - ZFSUP_AERO(i,KFLEV+1,1))

! TOA/SRF all sky anthropogenic forcing
     PSOLSWAERO(i,2) = (ZFSDN_AERO(i,1,2) - ZFSUP_AERO(i,1,2))-(ZFSDN_AERO(i,1,3) - ZFSUP_AERO(i,1,3))
     PTOPSWAERO(i,2) = (ZFSDN_AERO(i,KFLEV+1,2) - ZFSUP_AERO(i,KFLEV+1,2))- (ZFSDN_AERO(i,KFLEV+1,3) - ZFSUP_AERO(i,KFLEV+1,3))

! TOA/SRF clear sky natural forcing
     PSOLSW0AERO(i,1) = (ZFSDN0_AERO(i,1,3) - ZFSUP0_AERO(i,1,3))-(ZFSDN0_AERO(i,1,1) - ZFSUP0_AERO(i,1,1))
     PTOPSW0AERO(i,1) = (ZFSDN0_AERO(i,KFLEV+1,3) - ZFSUP0_AERO(i,KFLEV+1,3))-(ZFSDN0_AERO(i,KFLEV+1,1) - ZFSUP0_AERO(i,KFLEV+1,1))

! TOA/SRF clear sky anthropogenic forcing
     PSOLSW0AERO(i,2) = (ZFSDN0_AERO(i,1,2) - ZFSUP0_AERO(i,1,2))-(ZFSDN0_AERO(i,1,3) - ZFSUP0_AERO(i,1,3))
     PTOPSW0AERO(i,2) = (ZFSDN0_AERO(i,KFLEV+1,2) - ZFSUP0_AERO(i,KFLEV+1,2))-(ZFSDN0_AERO(i,KFLEV+1,3) - ZFSUP0_AERO(i,KFLEV+1,3))

! Cloud forcing indices 1: natural; 2 anthropogenic; 3: zero aerosol direct effect
! Instantaneously computed cloudy sky direct aerosol effect, cloud forcing due to aerosols above clouds
! natural
     PSOLSWCFAERO(i,1) = PSOLSWAERO(i,1) - PSOLSW0AERO(i,1)
     PTOPSWCFAERO(i,1) = PTOPSWAERO(i,1) - PTOPSW0AERO(i,1)

! Instantaneously computed cloudy SKY DIRECT aerosol effect, cloud forcing due to aerosols above clouds
! anthropogenic
     PSOLSWCFAERO(i,2) = PSOLSWAERO(i,2) - PSOLSW0AERO(i,2)
     PTOPSWCFAERO(i,2) = PTOPSWAERO(i,2) - PTOPSW0AERO(i,2)

! Cloudforcing without aerosol
! zero
     PSOLSWCFAERO(i,3) = (ZFSDN_AERO(i,1,1) - ZFSUP_AERO(i,1,1))-(ZFSDN0_AERO(i,1,1) - ZFSUP0_AERO(i,1,1))
     PTOPSWCFAERO(i,3) = (ZFSDN_AERO(i,KFLEV+1,1) - ZFSUP_AERO(i,KFLEV+1,1))- (ZFSDN0_AERO(i,KFLEV+1,1) - ZFSUP0_AERO(i,KFLEV+1,1))

! direct anthropogenic forcing , as in old LMDzT, however differences of net fluxes
     PSOLSWADAERO(i) = PSOLSWAERO(i,2)
     PTOPSWADAERO(i) = PTOPSWAERO(i,2)
     PSOLSWAD0AERO(i) = PSOLSW0AERO(i,2)
     PTOPSWAD0AERO(i) = PTOPSW0AERO(i,2)

ENDIF


IF (ok_aie) THEN
     PSOLSWAIAERO(i) = (ZFSDN_AERO(i,1,4) - ZFSUP_AERO(i,1,4))-(ZFSDN_AERO(i,1,2) - ZFSUP_AERO(i,1,2))
     PTOPSWAIAERO(i) = (ZFSDN_AERO(i,KFLEV+1,4) - ZFSUP_AERO(i,KFLEV+1,4))-(ZFSDN_AERO(i,KFLEV+1,2) - ZFSUP_AERO(i,KFLEV+1,2))
ENDIF

  ENDDO
END SUBROUTINE SW_AEROAR4
