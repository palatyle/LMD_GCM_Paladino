!
!$Header$
!
SUBROUTINE Init_Phys_lmdz(iim,jjp1,llm,nb_proc,distrib)
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz
  USE dimphy, ONLY : Init_dimphy
  IMPLICIT NONE
  
    INTEGER,INTENT(in) :: iim
    INTEGER,INTENT(in) :: jjp1
    INTEGER,INTENT(in) :: llm
    INTEGER,INTENT(in) :: nb_proc
    INTEGER,INTENT(in) :: distrib(0:nb_proc-1)


    CALL Init_grid_phy_lmdz(iim,jjp1,llm)
    CALL Init_phys_lmdz_para(iim,jjp1,nb_proc,distrib)
!$OMP PARALLEL
    CALL Init_dimphy(klon_omp,nbp_lev)
!$OMP END PARALLEL
 
END SUBROUTINE Init_Phys_lmdz  
