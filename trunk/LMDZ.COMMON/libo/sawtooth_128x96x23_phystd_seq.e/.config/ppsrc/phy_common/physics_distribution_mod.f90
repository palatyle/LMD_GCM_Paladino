










!
!$Id: physics_distribution_mod.F90 2351 2015-08-25 15:14:59Z emillour $
!
MODULE physics_distribution_mod


CONTAINS

  SUBROUTINE init_physics_distribution(grid_type, nvertex, &
                                       nbp, nbp_lon, nbp_lat, nbp_lev, &
                                       communicator)
  USE mod_phys_lmdz_para, ONLY: init_phys_lmdz_para, klon_omp
  USE mod_grid_phy_lmdz, ONLY: init_grid_phy_lmdz
  USE dimphy, ONLY : Init_dimphy

  IMPLICIT NONE
    INTEGER,INTENT(IN) :: grid_type 
    INTEGER,INTENT(IN) :: nvertex 
    INTEGER,INTENT(IN) :: nbp           
    INTEGER,INTENT(IN) :: nbp_lon
    INTEGER,INTENT(IN) :: nbp_lat
    INTEGER,INTENT(IN) :: nbp_lev
    INTEGER,INTENT(IN) :: communicator


    CALL init_grid_phy_lmdz(grid_type,nvertex, nbp_lon,nbp_lat,nbp_lev)
    CALL init_phys_lmdz_para(nbp,nbp_lon, nbp_lat, communicator)
!$OMP PARALLEL
    CALL init_dimphy(klon_omp,nbp_lev)
!$OMP END PARALLEL

  END SUBROUTINE init_physics_distribution  

!SUBROUTINE Init_Phys_lmdz(iim,jjp1,llm,nb_proc,distrib)
!  USE mod_phys_lmdz_para, ONLY: Init_phys_lmdz_para!, klon_omp
!  USE mod_grid_phy_lmdz, ONLY: Init_grid_phy_lmdz!, nbp_lev
!  USE dimphy, ONLY : Init_dimphy
!  USE infotrac_phy, ONLY : type_trac
!#ifdef REPROBUS
!  USE CHEM_REP, ONLY : Init_chem_rep_phys
!#endif

!  IMPLICIT NONE
  
!    INTEGER,INTENT(in) :: iim
!    INTEGER,INTENT(in) :: jjp1
!    INTEGER,INTENT(in) :: llm
!    INTEGER,INTENT(in) :: nb_proc
!    INTEGER,INTENT(in) :: distrib(0:nb_proc-1)


!    CALL Init_grid_phy_lmdz(iim,jjp1,llm)
!    CALL Init_phys_lmdz_para(iim,jjp1,nb_proc,distrib)
!!$OMP PARALLEL
!    CALL Init_dimphy(klon_omp,nbp_lev)
!
!! Initialization of Reprobus
!    IF (type_trac == 'repr') THEN
!#ifdef REPROBUS
!       CALL Init_chem_rep_phys(klon_omp,nbp_lev)
!#endif
!    END IF
!
!!$OMP END PARALLEL
 
!END SUBROUTINE Init_Phys_lmdz  








END MODULE physics_distribution_mod

