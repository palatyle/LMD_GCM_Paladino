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

   use surfdat_h, only: tsurf, co2ice, qsurf

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
     M_CO2ICE(i,j) = co2ice(subs)
     M_TSURF(i,j) = tsurf(subs)
     SELECT CASE (MARS_MODE)
      CASE (1,11,12)
        M_H2OICE(i,j) = qsurf(subs,2)  !! see above Tracer at surface
     END SELECT

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

   use comsoil_h, only: tsoil

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
     M_TSOIL(i,:,j) = tsoil(subs,:)

   ENDDO
   ENDDO

END SUBROUTINE update_outputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            M_FLUXRAD)

   use dimradmars_mod, only: fluxrad

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
     M_FLUXRAD(i,j) = fluxrad(subs)

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

   use turb_mod, only: q2,wstar,ustar,sensibFlux,&
                        hfmax_th,zmax_th
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
     HFMAX(i,j) = HFMAX_TH(subs)
     ZMAX(i,j) = ZMAX_TH(subs)
     USTM(i,j) = ustar(subs)
     HFX(i,j) = sensibFlux(subs)

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
     TAU_DUST(i,j) = comm_TAU_DUST(subs)
     QSURFDUST(i,j) = comm_QSURFDUST(subs)
     MTOT(i,j) = comm_MTOT(subs)
     ICETOT(i,j) = comm_ICETOT(subs)
     TAU_ICE(i,j) = comm_TAU_ICE(subs)
     HR_SW(i,kps:kpe,j) = comm_HR_SW(subs,kps:kpe)
     HR_LW(i,kps:kpe,j) = comm_HR_LW(subs,kps:kpe)
     RDUST(i,kps:kpe,j) = comm_RDUST(subs,kps:kpe)
     VMR_ICE(i,kps:kpe,j) = comm_VMR_ICE(subs,kps:kpe)
     RICE(i,kps:kpe,j) = comm_RICE(subs,kps:kpe)

   ENDDO
   ENDDO

   CALL deallocate_comm_wrf

END SUBROUTINE update_outputs_physiq_diag

END MODULE update_outputs_physiq_mod



