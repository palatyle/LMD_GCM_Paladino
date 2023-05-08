!
! $Id: phyredem.F 1457 2010-11-22 15:19:25Z musat $
!
c
      SUBROUTINE phyredem (fichnom)

      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE fonte_neige_mod,  ONLY : fonte_neige_final
      USE pbl_surface_mod,  ONLY : pbl_surface_final
      USE phys_state_var_mod
      USE iostart
      USE traclmdz_mod, ONLY : traclmdz_to_restart
      USE infotrac
      USE control_mod
      USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, co2_send

      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Ecriture de l'etat de redemarrage pour la physique
c======================================================================
#include "netcdf.inc"
#include "indicesol.h"
#include "dimsoil.h"
#include "clesphys.h"
#include "temps.h"
#include "thermcell.h"
#include "compbl.h"
c======================================================================
      CHARACTER*(*) fichnom

c les variables globales ecrites dans le fichier restart

      
      REAL tsoil(klon,nsoilmx,nbsrf)
      REAL tslab(klon), seaice(klon)
      REAL qsurf(klon,nbsrf)
      REAL qsol(klon)
      REAL snow(klon,nbsrf)
      REAL evap(klon,nbsrf)
      real fder(klon)
      REAL frugs(klon,nbsrf)
      REAL agesno(klon,nbsrf)
      REAL run_off_lic_0(klon)
      REAL trs(klon,nbtr)
c
      INTEGER nid, nvarid, idim1, idim2, idim3
      INTEGER ierr
      INTEGER length
      PARAMETER (length=100)
      REAL tab_cntrl(length)
c
      INTEGER isoil, nsrf
      CHARACTER (len=7) :: str7
      CHARACTER (len=2) :: str2
      INTEGER           :: it, iiq
      
c======================================================================
c 
c Get variables which will be written to restart file from module 
c pbl_surface_mod
      CALL pbl_surface_final(qsol, fder, snow, qsurf, 
     $     evap, frugs, agesno, tsoil)

c Get a variable calculated in module fonte_neige_mod
      CALL fonte_neige_final(run_off_lic_0)

c======================================================================

      CALL open_restartphy(fichnom)
      
      DO ierr = 1, length
         tab_cntrl(ierr) = 0.0
      ENDDO
      tab_cntrl(1) = dtime
      tab_cntrl(2) = radpas
c co2_ppm : current value of atmospheric CO2
      tab_cntrl(3) = co2_ppm
      tab_cntrl(4) = solaire
      tab_cntrl(5) = iflag_con
      tab_cntrl(6) = nbapp_rad

      IF( cycle_diurne ) tab_cntrl( 7 ) = 1.
      IF(   soil_model ) tab_cntrl( 8 ) = 1.
      IF(     new_oliq ) tab_cntrl( 9 ) = 1.
      IF(     ok_orodr ) tab_cntrl(10 ) = 1.
      IF(     ok_orolf ) tab_cntrl(11 ) = 1.

      tab_cntrl(13) = day_end
      tab_cntrl(14) = annee_ref
      tab_cntrl(15) = itau_phy

c co2_ppm0 : initial value of atmospheric CO2
      tab_cntrl(16) = co2_ppm0
c
      CALL put_var("controle","Parametres de controle",tab_cntrl)
c

      CALL put_field("longitude",
     .               "Longitudes de la grille physique",rlon)
     
      CALL put_field("latitude","Latitudes de la grille physique",rlat)

c
C PB ajout du masque terre/mer
C
      CALL put_field("masque","masque terre mer",zmasq)

c BP ajout des fraction de chaque sous-surface
C
C 1. fraction de terre 
C
      CALL put_field("FTER","fraction de continent",pctsrf(:,is_ter))
C 
C 2. Fraction de glace de terre
C 
      CALL put_field("FLIC","fraction glace de terre",pctsrf(:,is_lic))
C
C 3. fraction ocean
C
      CALL put_field("FOCE","fraction ocean",pctsrf(:,is_oce))
C
C 4. Fraction glace de mer
C
      CALL put_field("FSIC","fraction glace mer",pctsrf(:,is_sic))
C
C
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("TS"//str2,"Temperature de surface No."//str2,
     .                    ftsol(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO
c
      DO nsrf = 1, nbsrf
        DO isoil=1, nsoilmx
          IF (isoil.LE.99 .AND. nsrf.LE.99) THEN
            WRITE(str7,'(i2.2,"srf",i2.2)') isoil,nsrf
            CALL put_field("Tsoil"//str7,"Temperature du sol No."//str7,
     .                     tsoil(:,isoil,nsrf))
          ELSE
            PRINT*, "Trop de couches"
            CALL abort
          ENDIF
        ENDDO
      ENDDO
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("QS"//str2,"Humidite de surface No."//str2,
     .                   qsurf(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      END DO
C
      CALL put_field("QSOL","Eau dans le sol (mm)",qsol)
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("ALBE"//str2,"albedo de surface No."//str2,
     .                   falb1(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO

      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("ALBLW"//str2,"albedo LW de surface No."//str2,
     .                   falb2(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO
c
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("EVAP"//str2,"Evaporation de surface No."//str2
     .                   ,evap(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO

c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
          WRITE(str2,'(i2.2)') nsrf
          CALL put_field("SNOW"//str2,"Neige de surface No."//str2,
     .                   snow(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO

c
      CALL put_field("RADS","Rayonnement net a la surface",radsol)
c
      CALL put_field("solsw","Rayonnement solaire a la surface",solsw)
c
      CALL put_field("sollw","Rayonnement IF a la surface",sollw)
c
      CALL put_field("fder","Derive de flux",fder)
c
      CALL put_field("rain_f","precipitation liquide",rain_fall)
c
      CALL put_field("snow_f", "precipitation solide",snow_fall)
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
        WRITE(str2,'(i2.2)') nsrf
          CALL put_field("RUG"//str2,"rugosite de surface No."//str2,
     .         frugs(:,nsrf))
        ELSE
          PRINT*, "Trop de sous-mailles"
          CALL abort
        ENDIF
      ENDDO
c
      DO nsrf = 1, nbsrf
        IF (nsrf.LE.99) THEN
            WRITE(str2,'(i2.2)') nsrf
            CALL put_field("AGESNO"//str2,
     .                     "Age de la neige surface No."//str2,
     .                     agesno(:,nsrf))
        ELSE
            PRINT*, "Trop de sous-mailles"
            CALL abort
        ENDIF
      ENDDO
c
      CALL put_field("ZMEA","ZMEA",zmea)
c
      CALL put_field("ZSTD","ZSTD",zstd)
      
      CALL put_field("ZSIG","ZSIG",zsig)
      
      CALL put_field("ZGAM","ZGAM",zgam)
      
      CALL put_field("ZTHE","ZTHE",zthe)
      
      CALL put_field("ZPIC","ZPIC",zpic)
      
      CALL put_field("ZVAL","ZVAL",zval)
      
      CALL put_field("RUGSREL","RUGSREL",rugoro)
      
      CALL put_field("TANCIEN","TANCIEN",t_ancien)
      
      CALL put_field("QANCIEN","QANCIEN",q_ancien)
      
      CALL put_field("RUGMER","Longueur de rugosite sur mer",
     .               frugs(:,is_oce))
      
      CALL put_field("CLWCON","Eau liquide convective",clwcon(:,1))
      
      CALL put_field("RNEBCON","Nebulosite convective",rnebcon(:,1))
      
      CALL put_field("RATQS", "Ratqs",ratqs(:,1))
c
c run_off_lic_0
c
      CALL put_field("RUNOFFLIC0","Runofflic0",run_off_lic_0)
c
c
!!!!!!!!!!!!!!!!!!!! DEB TKE PBL !!!!!!!!!!!!!!!!!!!!!!!!!
c
      IF (iflag_pbl>1) then
        DO nsrf = 1, nbsrf
          IF (nsrf.LE.99) THEN
            WRITE(str2,'(i2.2)') nsrf
            CALL put_field("TKE"//str2,"Energ. Cineti. Turb."//str2,
     .                     pbl_tke(:,1:klev+1,nsrf))
          ELSE
            PRINT*, "Trop de sous-mailles"
            CALL abort
          ENDIF
        ENDDO
      ENDIF

!!!!!!!!!!!!!!!!!!!! FIN TKE PBL !!!!!!!!!!!!!!!!!!!!!!!!!
cIM ajout zmax0, f0, ema_work1, ema_work2
cIM wake_deltat, wake_deltaq, wake_s, wake_cstar, wake_pe, wake_fip
      
      CALL put_field("ZMAX0","ZMAX0",zmax0)
      
      CALL put_field("F0","F0",f0)
      
      CALL put_field("EMA_WORK1","EMA_WORK1",ema_work1)
      
      CALL put_field("EMA_WORK2","EMA_WORK2",ema_work2)
      
c wake_deltat
      CALL put_field("WAKE_DELTAT","WAKE_DELTAT",wake_deltat)

      CALL put_field("WAKE_DELTAQ","WAKE_DELTAQ",wake_deltaq)
      
      CALL put_field("WAKE_S","WAKE_S",wake_s)
      
      CALL put_field("WAKE_CSTAR","WAKE_CSTAR",wake_cstar)
      
      CALL put_field("WAKE_PE","WAKE_PE",wake_pe)

      CALL put_field("WAKE_FIP","WAKE_FIP",wake_fip)

c thermiques

      CALL put_field("FM_THERM","FM_THERM",fm_therm)

      CALL put_field("ENTR_THERM","ENTR_THERM",entr_therm)

      CALL put_field("DETR_THERM","DETR_THERM",detr_therm)

! trs from traclmdz_mod
      IF (type_trac == 'lmdz') THEN
         CALL traclmdz_to_restart(trs)
         DO it=1,nbtr
            iiq=niadv(it+2)
            CALL put_field("trs_"//tname(iiq),"",trs(:,it))
         END DO
         IF (carbon_cycle_cpl) THEN
            IF (.NOT. ALLOCATED(co2_send)) THEN
               ! This is the case of create_etat0_limit, ce0l
               ALLOCATE(co2_send(klon))
               co2_send(:) = co2_ppm0
            END IF
            CALL put_field("co2_send","co2_ppm for coupling",co2_send)
         END IF
      END IF

      CALL close_restartphy
!$OMP BARRIER
      RETURN
      END
