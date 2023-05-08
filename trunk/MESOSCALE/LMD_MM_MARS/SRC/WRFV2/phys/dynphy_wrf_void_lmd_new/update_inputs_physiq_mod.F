MODULE update_inputs_physiq_mod

CONTAINS

!SUBROUTINE update_inputs_physiq_time
!SUBROUTINE update_inputs_physiq_tracers
!SUBROUTINE update_inputs_physiq_constants
!SUBROUTINE update_inputs_physiq_geom
!SUBROUTINE update_inputs_physiq_surf
!SUBROUTINE update_inputs_physiq_soil
!SUBROUTINE update_inputs_physiq_turb
!SUBROUTINE update_inputs_physiq_rad
!SUBROUTINE update_inputs_physiq_slope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_time(&
            JULYR,JULDAY,GMT,&
            elaps,&
            lct_input,lon_input,ls_input,&
            MY)

  INTEGER, INTENT(IN) :: JULDAY, JULYR
  REAL, INTENT(IN) :: GMT,elaps,lon_input,ls_input,lct_input
  REAL,INTENT(OUT) :: MY

END SUBROUTINE update_inputs_physiq_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_tracers(nq,MARS_MODE)

  INTEGER, INTENT(IN) :: nq,MARS_MODE

END SUBROUTINE update_inputs_physiq_tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_constants

END SUBROUTINE update_inputs_physiq_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_geom( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,ngrid,nlayer,&
            DX,DY,MSFT,&
            lat_input, lon_input,&
            XLAT,XLONG)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,ngrid,nlayer
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     MSFT,XLAT,XLONG
   REAL, INTENT(IN) :: dx,dy
   REAL, INTENT(IN) :: lat_input, lon_input

END SUBROUTINE update_inputs_physiq_geom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_surf( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,MARS_MODE,&
            M_ALBEDO,CST_AL,&
            M_TSURF,M_EMISS,M_CO2ICE,&
            M_GW,M_Z0,CST_Z0,&
            M_H2OICE,&
            phisfi_val)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,MARS_MODE
   INTEGER :: i,j,subs,nlast
   REAL, INTENT(IN  ) :: CST_AL, phisfi_val, CST_Z0
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     M_ALBEDO,M_TSURF,M_EMISS,M_CO2ICE,M_H2OICE,M_Z0
   REAL, DIMENSION( ims:ime, 5, jms:jme ), INTENT(IN   )  :: M_GW  

END SUBROUTINE update_inputs_physiq_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_soil( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,nsoil,&
            M_TI,CST_TI,&
            M_ISOIL,M_DSOIL,&
            M_TSOIL,M_TSURF)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR,nsoil
   INTEGER :: i,j,subs,nlast,k
   REAL, INTENT(IN  ) :: CST_TI 
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: &
     M_TI, M_TSURF
   REAL, DIMENSION( ims:ime, nsoil, jms:jme ), INTENT(IN)  :: &
     M_TSOIL, M_ISOIL, M_DSOIL

END SUBROUTINE update_inputs_physiq_soil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_turb( &
            ims,ime,jms,jme,kms,kme,&
            ips,ipe,jps,jpe,&
            RESTART,isles,&
            M_Q2,M_WSTAR)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: M_WSTAR
   REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(IN) :: M_Q2
   LOGICAL, INTENT(IN ) :: RESTART,isles

END SUBROUTINE update_inputs_physiq_turb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_rad( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            RESTART,&
            M_FLUXRAD)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: M_FLUXRAD
   LOGICAL, INTENT(IN ) :: RESTART

END SUBROUTINE update_inputs_physiq_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_inputs_physiq_slope( &
            ims,ime,jms,jme,&
            ips,ipe,jps,jpe,&
            JULYR,&
            SLPX,SLPY)

   INTEGER, INTENT(IN) :: ims,ime,jms,jme
   INTEGER, INTENT(IN) :: ips,ipe,jps,jpe,JULYR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)  :: SLPX,SLPY

END SUBROUTINE update_inputs_physiq_slope

END MODULE update_inputs_physiq_mod
