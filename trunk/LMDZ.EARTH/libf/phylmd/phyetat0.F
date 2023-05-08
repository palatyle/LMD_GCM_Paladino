!
! $Id: phyetat0.F 1457 2010-11-22 15:19:25Z musat $
!
c
c
      SUBROUTINE phyetat0 (fichnom,
     .           clesphy0,
     .           tabcntr0)

      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE iophy
      USE ocean_cpl_mod,    ONLY : ocean_cpl_init
      USE fonte_neige_mod,  ONLY : fonte_neige_init
      USE pbl_surface_mod,  ONLY : pbl_surface_init
      USE surface_data,     ONLY : type_ocean
      USE phys_state_var_mod
      USE iostart
      USE write_field_phy
      USE infotrac
      USE traclmdz_mod,    ONLY : traclmdz_from_restart
      USE carbon_cycle_mod,ONLY :
     &     carbon_cycle_tr, carbon_cycle_cpl, co2_send

      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Lecture de l'etat initial pour la physique
c======================================================================
#include "dimensions.h"
#include "netcdf.inc"
#include "indicesol.h"
#include "dimsoil.h"
#include "clesphys.h"
#include "temps.h"
#include "thermcell.h"
#include "compbl.h"
c======================================================================
      CHARACTER*(*) fichnom

c les variables globales lues dans le fichier restart

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
      REAL fractint(klon)
      REAL trs(klon,nbtr)

      CHARACTER*6 ocean_in
      LOGICAL ok_veget_in

      INTEGER        longcles
      PARAMETER    ( longcles = 20 )
      REAL clesphy0( longcles )
c
      REAL xmin, xmax
c
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf, isoil ,k
      INTEGER length
      PARAMETER (length=100)
      INTEGER it, iiq
      REAL tab_cntrl(length), tabcntr0(length)
      CHARACTER*7 str7
      CHARACTER*2 str2
      LOGICAL :: found

c FH1D
c     real iolat(jjm+1)
      real iolat(jjm+1-1/iim)
c
c Ouvrir le fichier contenant l'etat initial:
c

     
      CALL open_startphy(fichnom)
      

c
c Lecture des parametres de controle:
c
      CALL get_var("controle",tab_cntrl)
       
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! Les constantes de la physiques sont lues dans la physique seulement.
! Les egalites du type
!             tab_cntrl( 5 )=clesphy0(1)
! sont remplacees par
!             clesphy0(1)=tab_cntrl( 5 )
! On inverse aussi la logique.
! On remplit les tab_cntrl avec les parametres lus dans les .def
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         DO i = 1, length
           tabcntr0( i ) = tab_cntrl( i )
         ENDDO
c
         tab_cntrl(1)=dtime
         tab_cntrl(2)=radpas

c co2_ppm : value from the previous time step
         IF (carbon_cycle_tr .OR. carbon_cycle_cpl) THEN
            co2_ppm = tab_cntrl(3)
            RCO2    = co2_ppm * 1.0e-06  * 44.011/28.97 
c ELSE : keep value from .def
         END IF

c co2_ppm0 : initial value of atmospheric CO2 (from create_etat0_limit.e .def)
         co2_ppm0   = tab_cntrl(16)

         solaire_etat0      = tab_cntrl(4)
         tab_cntrl(5)=iflag_con
         tab_cntrl(6)=nbapp_rad

         if (cycle_diurne) tab_cntrl( 7) =1.
         if (soil_model) tab_cntrl( 8) =1.
         if (new_oliq) tab_cntrl( 9) =1.
         if (ok_orodr) tab_cntrl(10) =1.
         if (ok_orolf) tab_cntrl(11) =1.
         if (ok_limitvrai) tab_cntrl(12) =1.


      itau_phy = tab_cntrl(15)

       
      clesphy0(1)=tab_cntrl( 5 )
      clesphy0(2)=tab_cntrl( 6 )
      clesphy0(3)=tab_cntrl( 7 )
      clesphy0(4)=tab_cntrl( 8 )
      clesphy0(5)=tab_cntrl( 9 )
      clesphy0(6)=tab_cntrl( 10 )
      clesphy0(7)=tab_cntrl( 11 )
      clesphy0(8)=tab_cntrl( 12 )

c
c Lecture des latitudes (coordonnees):
c
      CALL get_field("latitude",rlat)

c
c Lecture des longitudes (coordonnees):
c
      CALL get_field("longitude",rlon)

C
C
C Lecture du masque terre mer
C
      CALL get_field("masque",zmasq,found)
      IF (.NOT. found) THEN
        PRINT*, 'phyetat0: Le champ <masque> est absent'
        PRINT *, 'fichier startphy non compatible avec phyetat0'
      ENDIF

       
C Lecture des fractions pour chaque sous-surface
C
C initialisation des sous-surfaces
C
      pctsrf = 0.
C
C fraction de terre
C

      CALL get_field("FTER",pctsrf(:,is_ter),found)
      IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FTER> est absent'

C
C fraction de glace de terre
C
      CALL get_field("FLIC",pctsrf(:,is_lic),found)
      IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FLIC> est absent'

C
C fraction d'ocean
C
      CALL get_field("FOCE",pctsrf(:,is_oce),found)
      IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FOCE> est absent'

C
C fraction glace de mer
C
      CALL get_field("FSIC",pctsrf(:,is_sic),found)
      IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FSIC> est absent'

C
C  Verification de l'adequation entre le masque et les sous-surfaces
C
      fractint( 1 : klon) = pctsrf(1 : klon, is_ter) 
     $    + pctsrf(1 : klon, is_lic)
      DO i = 1 , klon
        IF ( abs(fractint(i) - zmasq(i) ) .GT. EPSFRA ) THEN
            WRITE(*,*) 'phyetat0: attention fraction terre pas ', 
     $          'coherente ', i, zmasq(i), pctsrf(i, is_ter)
     $          ,pctsrf(i, is_lic)
            WRITE(*,*) 'Je force la coherence zmasq=fractint'
            zmasq(i) = fractint(i)
        ENDIF 
      END DO 
      fractint (1 : klon) =  pctsrf(1 : klon, is_oce) 
     $    + pctsrf(1 : klon, is_sic)
      DO i = 1 , klon
        IF ( abs( fractint(i) - (1. - zmasq(i))) .GT. EPSFRA ) THEN
            WRITE(*,*) 'phyetat0 attention fraction ocean pas ', 
     $          'coherente ', i, zmasq(i) , pctsrf(i, is_oce)
     $          ,pctsrf(i, is_sic)
            WRITE(*,*) 'Je force la coherence zmasq=fractint'
            zmasq(i) = fractint(i)
        ENDIF 
      END DO 

C
c Lecture des temperatures du sol:
c

       CALL get_field("TS",ftsol(:,1),found)
       IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <TS> est absent'
         PRINT*, '          Mais je vais essayer de lire TS**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("TS"//str2,ftsol(:,nsrf))

           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(ftsol(i,nsrf),xmin)
              xmax = MAX(ftsol(i,nsrf),xmax)
           ENDDO
           PRINT*,'Temperature du sol TS**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <TS> est present'
         PRINT*, '          J ignore donc les autres temperatures TS**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(ftsol(i,1),xmin)
            xmax = MAX(ftsol(i,1),xmax)
         ENDDO
         PRINT*,'Temperature du sol <TS>', xmin, xmax
         DO nsrf = 2, nbsrf
         DO i = 1, klon
            ftsol(i,nsrf) = ftsol(i,1)
         ENDDO
         ENDDO
      ENDIF

c
c Lecture des temperatures du sol profond:
c
      DO nsrf = 1, nbsrf
        DO isoil=1, nsoilmx
          IF (isoil.GT.99 .AND. nsrf.GT.99) THEN
            PRINT*, "Trop de couches ou sous-mailles"
            CALL abort
          ENDIF
          WRITE(str7,'(i2.2,"srf",i2.2)') isoil, nsrf
          
          CALL get_field('Tsoil'//str7,tsoil(:,isoil,nsrf),found)
          IF (.NOT. found) THEN
            PRINT*, "phyetat0: Le champ <Tsoil"//str7//"> est absent"
            PRINT*, "          Il prend donc la valeur de surface"
            DO i=1, klon
               tsoil(i,isoil,nsrf)=ftsol(i,nsrf)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c
c Lecture de l'humidite de l'air juste au dessus du sol:
c

      CALL get_field("QS",qsurf(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <QS> est absent'
         PRINT*, '          Mais je vais essayer de lire QS**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("QS"//str2,qsurf(:,nsrf))
           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(qsurf(i,nsrf),xmin)
              xmax = MAX(qsurf(i,nsrf),xmax)
           ENDDO
           PRINT*,'Humidite pres du sol QS**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <QS> est present'
         PRINT*, '          J ignore donc les autres humidites QS**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(qsurf(i,1),xmin)
            xmax = MAX(qsurf(i,1),xmax)
         ENDDO
         PRINT*,'Humidite pres du sol <QS>', xmin, xmax
         DO nsrf = 2, nbsrf
           DO i = 1, klon
             qsurf(i,nsrf) = qsurf(i,1)
           ENDDO
         ENDDO
      ENDIF

C
C Eau dans le sol (pour le modele de sol "bucket")
C
      CALL get_field("QSOL",qsol,found)
      IF (.NOT. found) THEN
        PRINT*, 'phyetat0: Le champ <QSOL> est absent'
        PRINT*, '          Valeur par defaut nulle'
          qsol(:)=0.
      ENDIF

      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
        xmin = MIN(qsol(i),xmin)
        xmax = MAX(qsol(i),xmax)
      ENDDO
      PRINT*,'Eau dans le sol (mm) <QSOL>', xmin, xmax

c
c Lecture de neige au sol:
c

      CALL get_field("SNOW",snow(:,1),found)
      IF (.NOT. found) THEN
        PRINT*, 'phyetat0: Le champ <SNOW> est absent'
        PRINT*, '          Mais je vais essayer de lire SNOW**'
        DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
            PRINT*, "Trop de sous-mailles"
            CALL abort
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          CALL get_field( "SNOW"//str2,snow(:,nsrf))
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO i = 1, klon
            xmin = MIN(snow(i,nsrf),xmin)
            xmax = MAX(snow(i,nsrf),xmax)
          ENDDO
          PRINT*,'Neige du sol SNOW**:', nsrf, xmin, xmax
        ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <SNOW> est present'
         PRINT*, '          J ignore donc les autres neiges SNOW**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(snow(i,1),xmin)
            xmax = MAX(snow(i,1),xmax)
         ENDDO
         PRINT*,'Neige du sol <SNOW>', xmin, xmax
         DO nsrf = 2, nbsrf
         DO i = 1, klon
            snow(i,nsrf) = snow(i,1)
         ENDDO
         ENDDO
      ENDIF
c
c Lecture de albedo de l'interval visible au sol:
c
      CALL get_field("ALBE",falb1(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <ALBE> est absent'
         PRINT*, '          Mais je vais essayer de lire ALBE**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("ALBE"//str2,falb1(:,nsrf))
           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(falb1(i,nsrf),xmin)
              xmax = MAX(falb1(i,nsrf),xmax)
           ENDDO
           PRINT*,'Albedo du sol ALBE**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <ALBE> est present'
         PRINT*, '          J ignore donc les autres ALBE**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(falb1(i,1),xmin)
            xmax = MAX(falb1(i,1),xmax)
         ENDDO
         PRINT*,'Neige du sol <ALBE>', xmin, xmax
         DO nsrf = 2, nbsrf
           DO i = 1, klon
            falb1(i,nsrf) = falb1(i,1)
           ENDDO
         ENDDO
      ENDIF

c
c Lecture de albedo au sol dans l'interval proche infra-rouge:
c
      CALL get_field("ALBLW",falb2(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <ALBLW> est absent'
         PRINT*, '          Mais je vais prendre ALBE**'
         DO nsrf = 1, nbsrf
           DO i = 1, klon
             falb2(i,nsrf) = falb1(i,nsrf)
           ENDDO
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <ALBLW> est present'
         PRINT*, '          J ignore donc les autres ALBLW**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(falb2(i,1),xmin)
            xmax = MAX(falb2(i,1),xmax)
         ENDDO
         PRINT*,'Neige du sol <ALBLW>', xmin, xmax
         DO nsrf = 2, nbsrf
           DO i = 1, klon
             falb2(i,nsrf) = falb2(i,1)
           ENDDO
         ENDDO
      ENDIF
c
c Lecture de evaporation:  
c
      CALL get_field("EVAP",evap(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <EVAP> est absent'
         PRINT*, '          Mais je vais essayer de lire EVAP**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("EVAP"//str2, evap(:,nsrf))
           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(evap(i,nsrf),xmin)
              xmax = MAX(evap(i,nsrf),xmax)
           ENDDO
           PRINT*,'evap du sol EVAP**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <EVAP> est present'
         PRINT*, '          J ignore donc les autres EVAP**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(evap(i,1),xmin)
            xmax = MAX(evap(i,1),xmax)
         ENDDO
         PRINT*,'Evap du sol <EVAP>', xmin, xmax
         DO nsrf = 2, nbsrf
         DO i = 1, klon
            evap(i,nsrf) = evap(i,1)
         ENDDO
         ENDDO
      ENDIF
c
c Lecture precipitation liquide:
c
      CALL get_field("rain_f",rain_fall)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(rain_fall(i),xmin)
         xmax = MAX(rain_fall(i),xmax)
      ENDDO
      PRINT*,'Precipitation liquide rain_f:', xmin, xmax
c
c Lecture precipitation solide:
c
      CALL get_field("snow_f",snow_fall)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(snow_fall(i),xmin)
         xmax = MAX(snow_fall(i),xmax)
      ENDDO
      PRINT*,'Precipitation solide snow_f:', xmin, xmax
c
c Lecture rayonnement solaire au sol:
c
      CALL get_field("solsw",solsw,found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <solsw> est absent'
         PRINT*, 'mis a zero'
         solsw(:) = 0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(solsw(i),xmin)
         xmax = MAX(solsw(i),xmax)
      ENDDO
      PRINT*,'Rayonnement solaire au sol solsw:', xmin, xmax
c
c Lecture rayonnement IF au sol:
c
      CALL get_field("sollw",sollw,found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <sollw> est absent'
         PRINT*, 'mis a zero'
         sollw = 0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(sollw(i),xmin)
         xmax = MAX(sollw(i),xmax)
      ENDDO
      PRINT*,'Rayonnement IF au sol sollw:', xmin, xmax
      
c
c Lecture derive des flux:
c
      CALL get_field("fder",fder,found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <fder> est absent'
         PRINT*, 'mis a zero'
         fder = 0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(fder(i),xmin)
         xmax = MAX(fder(i),xmax)
      ENDDO
      PRINT*,'Derive des flux fder:', xmin, xmax

c
c Lecture du rayonnement net au sol:
c
      CALL get_field("RADS",radsol)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(radsol(i),xmin)
         xmax = MAX(radsol(i),xmax)
      ENDDO
      PRINT*,'Rayonnement net au sol radsol:', xmin, xmax
c
c Lecture de la longueur de rugosite 
c
c
      CALL get_field("RUG",frugs(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <RUG> est absent'
         PRINT*, '          Mais je vais essayer de lire RUG**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("RUG"//str2,frugs(:,nsrf))
           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(frugs(i,nsrf),xmin)
              xmax = MAX(frugs(i,nsrf),xmax)
           ENDDO
           PRINT*,'rugosite du sol RUG**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <RUG> est present'
         PRINT*, '          J ignore donc les autres RUG**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(frugs(i,1),xmin)
            xmax = MAX(frugs(i,1),xmax)
         ENDDO
         PRINT*,'rugosite <RUG>', xmin, xmax
         DO nsrf = 2, nbsrf
         DO i = 1, klon
            frugs(i,nsrf) = frugs(i,1)
         ENDDO
         ENDDO
      ENDIF

c
c Lecture de l'age de la neige:
c
      CALL get_field("AGESNO",agesno(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, 'phyetat0: Le champ <AGESNO> est absent'
         PRINT*, '          Mais je vais essayer de lire AGESNO**'
         DO nsrf = 1, nbsrf
           IF (nsrf.GT.99) THEN
             PRINT*, "Trop de sous-mailles"
             CALL abort
           ENDIF
           WRITE(str2,'(i2.2)') nsrf
           CALL get_field("AGESNO"//str2,agesno(:,nsrf),found)
           IF (.NOT. found) THEN
              PRINT*, "phyetat0: Le champ <AGESNO"//str2//"> est absent"
              agesno = 50.0
           ENDIF
           xmin = 1.0E+20
           xmax = -1.0E+20
           DO i = 1, klon
              xmin = MIN(agesno(i,nsrf),xmin)
              xmax = MAX(agesno(i,nsrf),xmax)
           ENDDO
           PRINT*,'Age de la neige AGESNO**:', nsrf, xmin, xmax
         ENDDO
      ELSE
         PRINT*, 'phyetat0: Le champ <AGESNO> est present'
         PRINT*, '          J ignore donc les autres AGESNO**'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(agesno(i,1),xmin)
            xmax = MAX(agesno(i,1),xmax)
         ENDDO
         PRINT*,'Age de la neige <AGESNO>', xmin, xmax
         DO nsrf = 2, nbsrf
         DO i = 1, klon
            agesno(i,nsrf) = agesno(i,1)
         ENDDO
         ENDDO
      ENDIF

c
      CALL get_field("ZMEA", zmea)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zmea(i),xmin)
         xmax = MAX(zmea(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zmea:', xmin, xmax
c
c
      CALL get_field("ZSTD",zstd)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zstd(i),xmin)
         xmax = MAX(zstd(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zstd:', xmin, xmax
c
c
      CALL get_field("ZSIG",zsig)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zsig(i),xmin)
         xmax = MAX(zsig(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zsig:', xmin, xmax
c
c
      CALL get_field("ZGAM",zgam)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zgam(i),xmin)
         xmax = MAX(zgam(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zgam:', xmin, xmax
c
c
      CALL get_field("ZTHE",zthe)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zthe(i),xmin)
         xmax = MAX(zthe(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zthe:', xmin, xmax
c
c
      CALL get_field("ZPIC",zpic)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zpic(i),xmin)
         xmax = MAX(zpic(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zpic:', xmin, xmax
c
      CALL get_field("ZVAL",zval)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zval(i),xmin)
         xmax = MAX(zval(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zval:', xmin, xmax
c
c
      CALL get_field("RUGSREL",rugoro)
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(rugoro(i),xmin)
         xmax = MAX(rugoro(i),xmax)
      ENDDO
      PRINT*,'Rugosite relief (ecart-type) rugsrel:', xmin, xmax
c
c
     
c
      ancien_ok = .TRUE.

      CALL get_field("TANCIEN",t_ancien,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <TANCIEN> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         ancien_ok = .FALSE.
      ENDIF


      CALL get_field("QANCIEN",q_ancien,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <QANCIEN> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         ancien_ok = .FALSE.
      ENDIF

      u_ancien = 0.0   !AXC: We don't have u_ancien and v_ancien in the start
      v_ancien = 0.0   !AXC: files, therefore they have to be initialized.
c

      clwcon=0.
      CALL get_field("CLWCON",clwcon(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ CLWCON est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(clwcon)
      xmax = MAXval(clwcon)
      PRINT*,'Eau liquide convective (ecart-type) clwcon:', xmin, xmax
c
      rnebcon = 0.
      CALL get_field("RNEBCON",rnebcon(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ RNEBCON est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(rnebcon)
      xmax = MAXval(rnebcon)
      PRINT*,'Nebulosite convective (ecart-type) rnebcon:', xmin, xmax

c
c Lecture ratqs
c
      ratqs=0.
      CALL get_field("RATQS",ratqs(:,1),found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <RATQS> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(ratqs)
      xmax = MAXval(ratqs)
      PRINT*,'(ecart-type) ratqs:', xmin, xmax
c
c Lecture run_off_lic_0
c
      CALL get_field("RUNOFFLIC0",run_off_lic_0,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <RUNOFFLIC0> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         run_off_lic_0 = 0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(run_off_lic_0)
      xmax = MAXval(run_off_lic_0)
      PRINT*,'(ecart-type) run_off_lic_0:', xmin, xmax


c Lecture de l'energie cinetique turbulente
c

      IF (iflag_pbl>1) then
        DO nsrf = 1, nbsrf
          IF (nsrf.GT.99) THEN
            PRINT*, "Trop de sous-mailles"
            CALL abort
          ENDIF
          WRITE(str2,'(i2.2)') nsrf
          CALL get_field("TKE"//str2,pbl_tke(:,1:klev+1,nsrf),found)
          IF (.NOT. found) THEN
            PRINT*, "phyetat0: <TKE"//str2//"> est absent"
            pbl_tke(:,:,nsrf)=1.e-8
          ENDIF
          xmin = 1.0E+20
          xmax = -1.0E+20
          DO k = 1, klev+1
            DO i = 1, klon
              xmin = MIN(pbl_tke(i,k,nsrf),xmin)
              xmax = MAX(pbl_tke(i,k,nsrf),xmax)
            ENDDO
          ENDDO
          PRINT*,'Temperature du sol TKE**:', nsrf, xmin, xmax
        ENDDO
      ENDIF
c
c zmax0
      CALL get_field("ZMAX0",zmax0,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <ZMAX0> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        zmax0=40.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(zmax0)
      xmax = MAXval(zmax0)
      PRINT*,'(ecart-type) zmax0:', xmin, xmax
c
c           f0(ig)=1.e-5
c f0
      CALL get_field("F0",f0,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <f0> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         f0=1.e-5
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(f0)
      xmax = MAXval(f0)
      PRINT*,'(ecart-type) f0:', xmin, xmax
c
c ema_work1
c
      CALL get_field("EMA_WORK1",ema_work1,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <EMA_WORK1> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        ema_work1=0.
      ELSE
        xmin = 1.0E+20
        xmax = -1.0E+20
        DO k = 1, klev
          DO i = 1, klon
            xmin = MIN(ema_work1(i,k),xmin)
            xmax = MAX(ema_work1(i,k),xmax)
          ENDDO
        ENDDO
        PRINT*,'ema_work1:', xmin, xmax
      ENDIF
c
c ema_work2
c
      CALL get_field("EMA_WORK2",ema_work2,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <EMA_WORK2> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        ema_work2=0.
      ELSE
        xmin = 1.0E+20
        xmax = -1.0E+20
        DO k = 1, klev
          DO i = 1, klon
            xmin = MIN(ema_work2(i,k),xmin)
            xmax = MAX(ema_work2(i,k),xmax)
          ENDDO
        ENDDO
        PRINT*,'ema_work2:', xmin, xmax
      ENDIF
c
c wake_deltat
c
      CALL get_field("WAKE_DELTAT",wake_deltat,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <WAKE_DELTAT> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        wake_deltat=0.
      ELSE
        xmin = 1.0E+20
        xmax = -1.0E+20
        DO k = 1, klev
          DO i = 1, klon
            xmin = MIN(wake_deltat(i,k),xmin)
            xmax = MAX(wake_deltat(i,k),xmax)
          ENDDO
        ENDDO
        PRINT*,'wake_deltat:', xmin, xmax
      ENDIF
c
c wake_deltaq
c   
      CALL get_field("WAKE_DELTAQ",wake_deltaq,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <WAKE_DELTAQ> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        wake_deltaq=0.
      ELSE
        xmin = 1.0E+20
        xmax = -1.0E+20
        DO k = 1, klev
          DO i = 1, klon
            xmin = MIN(wake_deltaq(i,k),xmin)
            xmax = MAX(wake_deltaq(i,k),xmax)
          ENDDO
        ENDDO
        PRINT*,'wake_deltaq:', xmin, xmax
      ENDIF
c
c wake_s
c
      CALL get_field("WAKE_S",wake_s,found)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <WAKE_S> est absent"
        PRINT*, "Depart legerement fausse. Mais je continue"
        wake_s=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(wake_s)
      xmax = MAXval(wake_s)
      PRINT*,'(ecart-type) wake_s:', xmin, xmax
c
c wake_cstar
c
      CALL get_field("WAKE_CSTAR",wake_cstar,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <WAKE_CSTAR> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         wake_cstar=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(wake_cstar)
      xmax = MAXval(wake_cstar)
      PRINT*,'(ecart-type) wake_cstar:', xmin, xmax
c
c wake_pe
c
      CALL get_field("WAKE_PE",wake_pe,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <WAKE_PE> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         wake_pe=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(wake_pe)
      xmax = MAXval(wake_pe)
      PRINT*,'(ecart-type) wake_pe:', xmin, xmax
c
c wake_fip
c
      CALL get_field("WAKE_FIP",wake_fip,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <WAKE_FIP> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         wake_fip=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(wake_fip)
      xmax = MAXval(wake_fip)
      PRINT*,'(ecart-type) wake_fip:', xmin, xmax
c
c  thermiques
c

      CALL get_field("FM_THERM",fm_therm,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <fm_therm> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         fm_therm=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(fm_therm)
      xmax = MAXval(fm_therm)
      PRINT*,'(ecart-type) fm_therm:', xmin, xmax

      CALL get_field("ENTR_THERM",entr_therm,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <entr_therm> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         entr_therm=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(entr_therm)
      xmax = MAXval(entr_therm)
      PRINT*,'(ecart-type) entr_therm:', xmin, xmax

      CALL get_field("DETR_THERM",detr_therm,found)
      IF (.NOT. found) THEN
         PRINT*, "phyetat0: Le champ <detr_therm> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         detr_therm=0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINval(detr_therm)
      xmax = MAXval(detr_therm)
      PRINT*,'(ecart-type) detr_therm:', xmin, xmax



c
c Read and send field trs to traclmdz
c
      IF (type_trac == 'lmdz') THEN
         DO it=1,nbtr
            iiq=niadv(it+2)
            CALL get_field("trs_"//tname(iiq),trs(:,it),found)
            IF (.NOT. found) THEN
               PRINT*, 
     $           "phyetat0: Le champ <trs_"//tname(iiq)//"> est absent"
               PRINT*, "Depart legerement fausse. Mais je continue"
               trs(:,it) = 0.
            ENDIF
            xmin = 1.0E+20
            xmax = -1.0E+20
            xmin = MINval(trs(:,it))
            xmax = MAXval(trs(:,it))
            PRINT*,"(ecart-type) trs_"//tname(iiq)//" :", xmin, xmax

         END DO
         CALL traclmdz_from_restart(trs)

         IF (carbon_cycle_cpl) THEN
            ALLOCATE(co2_send(klon), stat=ierr)
            IF (ierr /= 0) CALL abort_gcm
     &           ('phyetat0','pb allocation co2_send',1)
            CALL get_field("co2_send",co2_send,found)
            IF (.NOT. found) THEN
               PRINT*,"phyetat0: Le champ <co2_send> est absent"
               PRINT*,"Initialisation uniforme a co2_ppm=",co2_ppm
               co2_send(:) = co2_ppm
            END IF
         END IF
      END IF


c on ferme le fichier
      CALL close_startphy

      CALL init_iophy_new(rlat,rlon)
      	

c
c Initialize module pbl_surface_mod 
c
      CALL pbl_surface_init(qsol, fder, snow, qsurf,
     $     evap, frugs, agesno, tsoil)

c Initialize module ocean_cpl_mod for the case of coupled ocean
      IF ( type_ocean == 'couple' ) THEN
         CALL ocean_cpl_init(dtime, rlon, rlat)
      ENDIF
c
c Initilialize module fonte_neige_mod      
c
      CALL fonte_neige_init(run_off_lic_0)


      RETURN
      END
