!
! $Id: $
!
      SUBROUTINE phyredem (fichnom)

      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE iophy
      USE phys_state_var_mod
      USE iostart, only : open_restartphy,close_restartphy, & 
                          put_var,put_field
      use geometry_mod, only: longitude_deg, latitude_deg
      USE time_phylmdz_mod, only: day_end, annee_ref, itau_phy, raz_date

      implicit none
!======================================================================
! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
! Objet: Lecture de l'etat initial pour la physique
!======================================================================
#include "netcdf.inc"
#include "dimsoil.h"
#include "clesphys.h"
#include "tabcontrol.h"
!======================================================================

character(len=*),intent(in) :: fichnom
REAL    :: tab_cntrl(length)
integer :: isoil
CHARACTER(len=2) :: str2


! open file

      CALL open_restartphy(fichnom)

! tab_cntrl() contains run parameters

      tab_cntrl(:)=0.0
 
      tab_cntrl(1) = dtime
      tab_cntrl(2) = radpas
      tab_cntrl(3) = 0.0
      tab_cntrl(4) = solaire
      tab_cntrl(5) = 0
      tab_cntrl(6) = nbapp_rad

      IF( cycle_diurne ) tab_cntrl( 7 ) = 1.
      IF(   soil_model ) tab_cntrl( 8 ) = 1.
      IF(     ok_orodr ) tab_cntrl(10 ) = 1.
      IF(     ok_orolf ) tab_cntrl(11 ) = 1.
      IF( ok_gw_nonoro ) tab_cntrl(12 ) = 1.

      tab_cntrl(13) = day_end
      tab_cntrl(14) = annee_ref
      tab_cntrl(15) = itau_phy

      CALL put_var("controle","Parametres de controle",tab_cntrl)

! coordinates

      CALL put_field("longitude", &
                     "Longitudes de la grille physique",longitude_deg)
     
      CALL put_field("latitude", &
                     "Latitudes de la grille physique",latitude_deg)

! variables

      CALL put_field("TS","Temperature de surface",ftsol)

      DO isoil=1, nsoilmx
        IF (isoil.LE.99) THEN
        WRITE(str2,'(i2.2)') isoil
        CALL put_field("Tsoil"//str2, &
                       "Temperature du sol No."//str2,ftsoil(:,isoil))
        ELSE
        PRINT*, "Trop de couches"
        CALL abort
        ENDIF
      ENDDO

      CALL put_field("ALBE","albedo de surface",falbe)
      CALL put_field("solsw","Rayonnement solaire a la surface",solsw)
      CALL put_field("sollw","Rayonnement IR a la surface",sollw)
      CALL put_field("fder","Derive de flux",fder)
      CALL put_field("dlw","Derivee flux IR",dlw)
      CALL put_field("sollwdown","Flux IR vers le bas a la surface",sollwdown)
      CALL put_field("RADS","Rayonnement net a la surface",radsol)
      CALL put_field("ZMEA","zmea Orographie sous-maille",zmea)
      CALL put_field("ZSTD","zstd Orographie sous-maille",zstd)
      CALL put_field("ZSIG","zsig Orographie sous-maille",zsig)
      CALL put_field("ZGAM","zgam Orographie sous-maille",zgam)
      CALL put_field("ZTHE","zthe Orographie sous-maille",zthe)
      CALL put_field("ZPIC","zpic Orographie sous-maille",zpic)
      CALL put_field("ZVAL","zval Orographie sous-maille",zval)

      CALL put_field("TANCIEN","T Previous iteration",t_ancien)

! close file

      CALL close_restartphy
!$OMP BARRIER

      END SUBROUTINE phyredem
