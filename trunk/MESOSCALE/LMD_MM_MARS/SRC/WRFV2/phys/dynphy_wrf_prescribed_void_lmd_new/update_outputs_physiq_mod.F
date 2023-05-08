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

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   INTEGER, INTENT(IN) :: MARS_MODE
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &
     M_TSURF,M_CO2ICE,M_H2OICE

END SUBROUTINE update_outputs_physiq_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_soil( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            nsoil,&
            M_TSOIL)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,nsoil
   REAL, DIMENSION( ims:ime, nsoil, jms:jme ), INTENT(INOUT)  :: &
     M_TSOIL

END SUBROUTINE update_outputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            M_FLUXRAD)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: M_FLUXRAD

END SUBROUTINE update_outputs_physiq_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_outputs_physiq_turb( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,kps,kpe,&
            M_Q2,M_WSTAR,&
            HFMAX,ZMAX,USTM,HFX)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,kps,kpe
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &
     M_WSTAR,HFMAX,ZMAX,USTM,HFX
   REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(INOUT) :: M_Q2

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
   REAL, DIMENSION(kms:kme) :: phr_sw,phr_lw,phr_dyn
   INTEGER :: i,j,subs

  open(17,file='prescribed_sw.txt',form='formatted',status='old')
  open(18,file='prescribed_lw.txt',form='formatted',status='old')
  open(19,file='prescribed_dyn.txt',form='formatted',status='old')
  rewind(17)
  rewind(18)
  rewind(19)
  DO k=1,kpe-kps+1
    read(17,*) phr_sw(k)
    read(18,*) phr_lw(k)
    read(19,*) phr_dyn(k)
  ENDDO
  close(19)
  close(18)
  close(17)
 !print*,'phr_lw',phr_lw 
  DO j=jps,jpe 
  DO i=ips,ipe
    HR_SW(i,:,j)=phr_sw(:)
    HR_LW(i,:,j)=phr_lw(:)
    HR_DYN(i,:,j)=phr_dyn(:)
  ENDDO
  ENDDO

END SUBROUTINE update_outputs_physiq_diag

END MODULE update_outputs_physiq_mod



