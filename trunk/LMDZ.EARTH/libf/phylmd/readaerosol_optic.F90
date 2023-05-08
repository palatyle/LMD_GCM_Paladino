! $Id$
!
SUBROUTINE readaerosol_optic(debut, new_aod, flag_aerosol, itap, rjourvrai, &
     pdtphys, pplay, paprs, t_seri, rhcl, presnivs, &
     mass_solu_aero, mass_solu_aero_pi, &
     tau_aero, piz_aero, cg_aero, &
     tausum_aero, tau3d_aero )

! This routine will :
! 1) recevie the aerosols(already read and interpolated) corresponding to flag_aerosol
! 2) calculate the optical properties for the aerosols
!
  
  USE dimphy
  USE aero_mod
  USE phys_local_var_mod, only: sconcso4,sconcoa,sconcbc,sconcss,sconcdust, &
      concso4,concoa,concbc,concss,concdust,loadso4,loadoa,loadbc,loadss,loaddust, &
      load_tmp1,load_tmp2,load_tmp3,load_tmp4,load_tmp5,load_tmp6,load_tmp7
  IMPLICIT NONE

! Input arguments
!****************************************************************************************
  LOGICAL, INTENT(IN)                      :: debut
  LOGICAL, INTENT(IN)                      :: new_aod
  INTEGER, INTENT(IN)                      :: flag_aerosol
  INTEGER, INTENT(IN)                      :: itap
  REAL, INTENT(IN)                         :: rjourvrai
  REAL, INTENT(IN)                         :: pdtphys
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay
  REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: t_seri
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: rhcl   ! humidite relative ciel clair
  REAL, DIMENSION(klev), INTENT(IN)        :: presnivs

! Output arguments
!****************************************************************************************
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero    ! Total mass for all soluble aerosols
  REAL, DIMENSION(klon,klev), INTENT(OUT)     :: mass_solu_aero_pi !     -"-     preindustrial values
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(OUT) :: tau_aero    ! Aerosol optical thickness
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(OUT) :: piz_aero    ! Single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(OUT) :: cg_aero     ! asymmetry parameter aerosol
  REAL, DIMENSION(klon,nwave,naero_spc), INTENT(OUT)       :: tausum_aero
  REAL, DIMENSION(klon,klev,nwave,naero_spc), INTENT(OUT)  :: tau3d_aero

! Local variables
!****************************************************************************************
  REAL, DIMENSION(klon)        :: aerindex ! POLDER aerosol index 
  REAL, DIMENSION(klon,klev)   :: sulfate  ! SO4 aerosol concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: bcsol    ! BC soluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: bcins    ! BC insoluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: pomsol   ! POM soluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: pomins   ! POM insoluble concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: cidust    ! DUST aerosol concentration  [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sscoarse  ! SS Coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sssupco   ! SS Super Coarse concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: ssacu     ! SS Acumulation concentration [ug/m3]
  REAL, DIMENSION(klon,klev)   :: sulfate_pi
  REAL, DIMENSION(klon,klev)   :: bcsol_pi
  REAL, DIMENSION(klon,klev)   :: bcins_pi
  REAL, DIMENSION(klon,klev)   :: pomsol_pi
  REAL, DIMENSION(klon,klev)   :: pomins_pi
  REAL, DIMENSION(klon,klev)   :: cidust_pi
  REAL, DIMENSION(klon,klev)   :: sscoarse_pi
  REAL, DIMENSION(klon,klev)   :: sssupco_pi
  REAL, DIMENSION(klon,klev)   :: ssacu_pi
  REAL, DIMENSION(klon,klev)   :: pdel
  REAL, DIMENSION(klon,klev,naero_spc) :: m_allaer
  REAL, DIMENSION(klon,klev,naero_spc) :: m_allaer_pi !RAF  
!  REAL, DIMENSION(klon,naero_tot)      :: fractnat_allaer !RAF delete??

  INTEGER :: k, i
  
!****************************************************************************************
! 1) Get aerosol mass
!    
!****************************************************************************************
! Read and interpolate sulfate
  IF ( flag_aerosol .EQ. 1 .OR. &
       flag_aerosol .EQ. 6 ) THEN 

     CALL readaerosol_interp(id_ASSO4M, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, sulfate, sulfate_pi,loadso4)
  ELSE
     sulfate(:,:) = 0. ; sulfate_pi(:,:) = 0.
     loadso4=0.
  END IF

! Read and interpolate bcsol and bcins
  IF ( flag_aerosol .EQ. 2 .OR. &
       flag_aerosol .EQ. 6 ) THEN 

     ! Get bc aerosol distribution 
     CALL readaerosol_interp(id_ASBCM, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, bcsol, bcsol_pi, load_tmp1 )
     CALL readaerosol_interp(id_AIBCM, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, bcins, bcins_pi, load_tmp2 )
     loadbc(:)=load_tmp1(:)+load_tmp2(:)
  ELSE
     bcsol(:,:) = 0. ; bcsol_pi(:,:) = 0.
     bcins(:,:) = 0. ; bcins_pi(:,:) = 0.
     loadbc=0.
  END IF


! Read and interpolate pomsol and pomins
  IF ( flag_aerosol .EQ. 3 .OR. &
       flag_aerosol .EQ. 6 ) THEN

     CALL readaerosol_interp(id_ASPOMM, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, pomsol, pomsol_pi, load_tmp3)
     CALL readaerosol_interp(id_AIPOMM, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, pomins, pomins_pi, load_tmp4)
     loadoa(:)=load_tmp3(:)+load_tmp4(:)
  ELSE
     pomsol(:,:) = 0. ; pomsol_pi(:,:) = 0.
     pomins(:,:) = 0. ; pomins_pi(:,:) = 0.
     loadoa=0.
  END IF


! Read and interpolate csssm, ssssm, assssm
  IF (flag_aerosol .EQ. 4 .OR. &
      flag_aerosol .EQ. 6 ) THEN 

      CALL readaerosol_interp(id_SSSSM ,itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, sssupco, sssupco_pi, load_tmp5) 
      CALL readaerosol_interp(id_CSSSM ,itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, sscoarse,sscoarse_pi, load_tmp6) 
      CALL readaerosol_interp(id_ASSSM ,itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, ssacu, ssacu_pi, load_tmp7) 
     loadss(:)=load_tmp5(:)+load_tmp6(:)+load_tmp7(:)
  ELSE
     sscoarse(:,:) = 0. ; sscoarse_pi(:,:) = 0. 
     ssacu(:,:)    = 0. ; ssacu_pi(:,:) = 0. 
     sssupco(:,:)  = 0. ; sssupco_pi = 0. 
     loadss=0.
  ENDIF

! Read and interpolate cidustm
  IF (flag_aerosol .EQ. 5 .OR.  &
      flag_aerosol .EQ. 6 ) THEN 

      CALL readaerosol_interp(id_CIDUSTM, itap, pdtphys, rjourvrai, debut, pplay, paprs, t_seri, cidust, cidust_pi, loaddust) 

  ELSE
      cidust(:,:) = 0. ; cidust_pi(:,:) = 0. 
      loaddust=0.
  ENDIF

!
! Store all aerosols in one variable
!
  m_allaer(:,:,id_ASBCM)  = bcsol(:,:)        ! ASBCM
  m_allaer(:,:,id_ASPOMM) = pomsol(:,:)       ! ASPOMM
  m_allaer(:,:,id_ASSO4M) = sulfate(:,:)      ! ASSO4M (= SO4) 
  m_allaer(:,:,id_CSSO4M) = 0.                ! CSSO4M 
  m_allaer(:,:,id_SSSSM)  = sssupco(:,:)      ! SSSSM
  m_allaer(:,:,id_CSSSM)  = sscoarse(:,:)     ! CSSSM
  m_allaer(:,:,id_ASSSM)  = ssacu(:,:)        ! ASSSM
  m_allaer(:,:,id_CIDUSTM)= cidust(:,:)       ! CIDUSTM
  m_allaer(:,:,id_AIBCM)  = bcins(:,:)        ! AIBCM
  m_allaer(:,:,id_AIPOMM) = pomins(:,:)       ! AIPOMM

!RAF
  m_allaer_pi(:,:,1)  = bcsol_pi(:,:)        ! ASBCM pre-ind
  m_allaer_pi(:,:,2)  = pomsol_pi(:,:)       ! ASPOMM pre-ind
  m_allaer_pi(:,:,3)  = sulfate_pi(:,:)      ! ASSO4M (= SO4) pre-ind
  m_allaer_pi(:,:,4)  = 0.                ! CSSO4M pre-ind
  m_allaer_pi(:,:,5)  = sssupco_pi(:,:)      ! SSSSM pre-ind
  m_allaer_pi(:,:,6)  = sscoarse_pi(:,:)     ! CSSSM pre-ind
  m_allaer_pi(:,:,7)  = ssacu_pi(:,:)        ! ASSSM pre-ind
  m_allaer_pi(:,:,8)  = cidust_pi(:,:)       ! CIDUSTM pre-ind
  m_allaer_pi(:,:,9)  = bcins_pi(:,:)        ! AIBCM pre-ind
  m_allaer_pi(:,:,10) = pomins_pi(:,:)       ! AIPOMM pre-ind

!
! Calculate the total mass of all soluble aersosols
!
  mass_solu_aero(:,:)    = sulfate(:,:)    + bcsol(:,:)    + pomsol(:,:) !   + &
!       sscoarse(:,:)    + ssacu(:,:)    + sssupco(:,:) 
  mass_solu_aero_pi(:,:) = sulfate_pi(:,:) + bcsol_pi(:,:) + pomsol_pi(:,:) ! + &
!       sscoarse_pi(:,:) + ssacu_pi(:,:) + sssupco_pi(:,:)

!****************************************************************************************
! 2) Calculate optical properties for the aerosols
!
!****************************************************************************************
  DO k = 1, klev
     DO i = 1, klon
        pdel(i,k) = paprs(i,k) - paprs (i,k+1)
     END DO
  END DO

  IF (new_aod) THEN 

! RAF delete??     fractnat_allaer(:,:) = 0.
! RAF fractnat_allaer -> m_allaer_pi

     CALL aeropt_2bands( &
          pdel, m_allaer, pdtphys, rhcl, & 
          tau_aero, piz_aero, cg_aero,   &
          m_allaer_pi, flag_aerosol, &
          pplay, t_seri, presnivs) 
     
     ! aeropt_5wv only for validation and diagnostics.
     CALL aeropt_5wv(                    &
          pdel, m_allaer,                &
          pdtphys, rhcl, aerindex,       & 
          flag_aerosol, pplay, t_seri,   &
          tausum_aero, tau3d_aero, presnivs)
  ELSE

     CALL aeropt(pplay, paprs, t_seri, sulfate, rhcl, &
          tau_aero(:,:,id_ASSO4M,:), piz_aero(:,:,id_ASSO4M,:), cg_aero(:,:,id_ASSO4M,:), aerindex)
     
  END IF


! Diagnostics calculation for CMIP5 protocol
  sconcso4(:)=m_allaer(:,1,id_ASSO4M)*1.e-9
  sconcoa(:)=(m_allaer(:,1,id_ASPOMM)+m_allaer(:,1,id_AIPOMM))*1.e-9
  sconcbc(:)=(m_allaer(:,1,id_ASBCM)+m_allaer(:,1,id_AIBCM))*1.e-9
  sconcss(:)=(m_allaer(:,1,id_ASSSM)+m_allaer(:,1,id_CSSSM)+m_allaer(:,1,id_SSSSM))*1.e-9
  sconcdust(:)=m_allaer(:,1,id_CIDUSTM)*1.e-9
  concso4(:,:)=m_allaer(:,:,id_ASSO4M)*1.e-9
  concoa(:,:)=(m_allaer(:,:,id_ASPOMM)+m_allaer(:,:,id_AIPOMM))*1.e-9
  concbc(:,:)=(m_allaer(:,:,id_ASBCM)+m_allaer(:,:,id_AIBCM))*1.e-9
  concss(:,:)=(m_allaer(:,:,id_ASSSM)+m_allaer(:,:,id_CSSSM)+m_allaer(:,:,id_SSSSM))*1.e-9
  concdust(:,:)=m_allaer(:,:,id_CIDUSTM)*1.e-9


END SUBROUTINE readaerosol_optic
