MODULE update_outputs_physiq_mod

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_surf( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            MARS_MODE,&
            M_TSURF,M_CO2ICE,&
            M_H2OICE)

   use phys_state_var_mod, only : ftsol

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER, INTENT(IN) :: MARS_MODE
   INTEGER :: i,j,subs
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &
     M_TSURF,M_CO2ICE,M_H2OICE

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !-------------------------------------------------------!
     ! Save key variables for restart and output and nesting !  
     !-------------------------------------------------------!
     M_TSURF(i,j) = ftsol(subs)

   ENDDO
   ENDDO

END SUBROUTINE update_outputs_physiq_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_soil( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            nsoil,&
            M_TSOIL)

   use phys_state_var_mod, only : ftsoil

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,nsoil
   INTEGER :: i,j,subs
   REAL, DIMENSION( ims:ime, nsoil, jms:jme ), INTENT(INOUT)  :: &
     M_TSOIL

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !-------------------------------------------------------!
     ! Save key variables for restart and output and nesting !  
     !-------------------------------------------------------!
     M_TSOIL(i,:,j) = ftsoil(subs,:)

   ENDDO
   ENDDO

END SUBROUTINE update_outputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            M_FLUXRAD)

   use phys_state_var_mod, only : radsol

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER :: i,j,subs
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: M_FLUXRAD

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !-------------------------------------------------------!
     ! Save key variables for restart and output and nesting !  
     !-------------------------------------------------------!
     M_FLUXRAD(i,j) = radsol(subs)

   ENDDO
   ENDDO

END SUBROUTINE update_outputs_physiq_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_turb( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,kps,kpe,&
            M_Q2,M_WSTAR,&
            HFMAX,ZMAX,USTM,HFX)

   use turb_mod, only: q2,wstar,yustar,sens

   implicit none

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,kps,kpe
   INTEGER :: i,j,subs
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &
     M_WSTAR,HFMAX,ZMAX,USTM,HFX
   REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(INOUT) :: M_Q2

   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !-------------------------------------------------------!
     ! Save key variables for restart and output and nesting !  
     !-------------------------------------------------------!
     M_Q2(i,kps:kpe+1,j) = q2(subs,:)
     M_WSTAR(i,j) = wstar(subs)

     !! output only (arrays already in phys modules)
     USTM(i,j) = yustar(subs)
     HFX(i,j) = sens(subs)

  ENDDO
  ENDDO

END SUBROUTINE update_outputs_physiq_turb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_diag( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,kps,kpe,&
            SWDOWNZ,TAU_DUST,QSURFDUST,&
            MTOT,ICETOT,TAU_ICE,&
            HR_SW,HR_LW,HR_DYN,DT,DTRAD,&
            RDUST,VMR_ICE,RICE,&
            CLOUDFRAC,TOTCLOUDFRAC,&
            RAIN,SNOW,REEVAP,SURFRAIN,ALBEQ,FLUXTOP_DN,FLUXABS_SW,FLUXTOP_LW,FLUXSURF_SW,&
            FLUXSURF_LW,FLXGRD,LSCEZ,H2OICE_REFF,LATENT_HF)
   
   USE comm_wrf !! to get fields to be written from physiq

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,kps,kpe
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
     SWDOWNZ,TAU_DUST,QSURFDUST,&
     MTOT,ICETOT,TAU_ICE,&
     TOTCLOUDFRAC,ALBEQ,FLUXTOP_DN,FLUXABS_SW,FLUXTOP_LW,FLUXSURF_SW,&
     FLUXSURF_LW,FLXGRD,LATENT_HF,REEVAP,SURFRAIN
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) :: &
     HR_SW,HR_LW,RDUST,VMR_ICE,RICE,CLOUDFRAC,HR_DYN,DT,DTRAD,RAIN,SNOW,&
     LSCEZ,H2OICE_REFF
   INTEGER :: i,j,subs


   DO j = jps,jpe
   DO i = ips,ipe

     !-----------------------------------!
     ! 1D subscript for physics "cursor" !
     !-----------------------------------!
     subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

     !! get diagnostics from physics
     SWDOWNZ(i,j) = comm_SWDOWNZ(subs)
     MTOT(i,j) = comm_MTOT(subs) 
     HR_SW(i,kps:kpe,j) = comm_HR_SW(subs,kps:kpe)
     HR_LW(i,kps:kpe,j) = comm_HR_LW(subs,kps:kpe)
     HR_DYN(i,kps:kpe,j) = comm_HR_DYN(subs,kps:kpe)
     DT(i,kps:kpe,j) = comm_DT(subs,kps:kpe)
     DTRAD(i,kps:kpe,j) = comm_DT_RAD(subs,kps:kpe)
     DTVDF(i,kps:kpe,j) = comm_DT_VDF(subs,kps:kpe)
     DTAJS(i,kps:kpe,j) = comm_DT_AJS(subs,kps:kpe)

   ENDDO
   ENDDO

   CALL deallocate_comm_wrf

END SUBROUTINE update_outputs_physiq_diag

END MODULE update_outputs_physiq_mod



