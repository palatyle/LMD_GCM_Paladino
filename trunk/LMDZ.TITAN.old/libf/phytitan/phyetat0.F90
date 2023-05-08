!
! $Id $
!
subroutine phyetat0(fichnom)
! Load initial state for the physics
! and do some resulting initializations

      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE iophy
      USE phys_state_var_mod
      USE iostart
      USE geometry_mod,  only: latitude_deg,longitude_deg
      USE time_phylmdz_mod, only: itau_phy, raz_date

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

character(len=*),intent(in) :: fichnom ! input file name
REAL    :: xmin, xmax
LOGICAL :: found
REAL    :: tab_cntrl(length)
integer :: i,isoil
CHARACTER(len=2) :: str2
REAL :: lon_startphy(klon), lat_startphy(klon)

! les variables globales lues dans le fichier restart

! open physics initial state file:
call open_startphy(fichnom)

!
! Lecture des parametres de controle:
!
      CALL get_var("controle",tab_cntrl,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <controle> est absent'
         CALL abort
      ENDIF
       
      DO i = 1, length
           tabcntr0( i ) = tab_cntrl( i )
      ENDDO


      dtime        = tab_cntrl(1)
      radpas       = tab_cntrl(2)
      chimpas      = tab_cntrl(3)
      lsinit       = tab_cntrl(17)

      itau_phy = tab_cntrl(15)

! Attention si raz_date est active :
! il faut remettre a zero itau_phy apres phyetat0 !
! et verifier que lsinit est proche de 0.
      IF (raz_date.eq.1) THEN
        itau_phy=0
        if ((lsinit.gt.3.).and.(lsinit.lt.357.)) then
          PRINT*, 'phyetat0: raz_date=1 and ls different from 0.'
          PRINT*, 'When raz_date=1, we reset the initial date'
          PRINT*, 'to spring equinox, Ls=0., so the start files'
          PRINT*, 'should be within a couple of degrees from Ls=0.'
          PRINT*, 'or the circulation will be too far from equilibrium'
          CALL abort
        endif
      ENDIF

! read latitudes and make a sanity check (because already known from dyn)
call get_field("latitude",lat_startphy,found)
IF (.not.found) THEN
  PRINT*, 'phyetat0: Le champ <latitude> est absent'
  CALL abort
ENDIF
DO i=1,klon
  IF (ABS(lat_startphy(i)-latitude_deg(i))>=0.01) THEN
    WRITE(*,*) "phyetat0: Warning! Latitude discrepancy wrt startphy file:",&
               " i=",i," lat_startphy(i)=",lat_startphy(i),&
               " latitude_deg(i)=",latitude_deg(i)
    CALL abort
  ENDIF
ENDDO

! read longitudes and make a sanity check (because already known from dyn)
call get_field("longitude",lon_startphy,found)
IF (.not.found) THEN
  PRINT*, 'phyetat0: Le champ <longitude> est absent'
  CALL abort
ENDIF
DO i=1,klon
  IF (ABS(lon_startphy(i)-longitude_deg(i))>=0.01) THEN
    WRITE(*,*) "phyetat0: Warning! Longitude discrepancy wrt startphy file:",&
               " i=",i," lon_startphy(i)=",lon_startphy(i),&
               " longitude_deg(i)=",longitude_deg(i)
    CALL abort
  ENDIF
ENDDO

! read in other variables here ...

! Lecture des temperatures du sol:

       CALL get_field("TS",ftsol(:),found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <TS> est absent'
         PRINT*, "phyetat0: Lecture echouee pour <TS>"
         CALL abort
      ELSE
         PRINT*, 'phyetat0: Le champ <TS> est present'
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(ftsol(i),xmin)
            xmax = MAX(ftsol(i),xmax)
         ENDDO
         PRINT*,'Temperature du sol <TS>', xmin, xmax
      ENDIF


! Lecture des temperatures du sol profond:

      DO isoil=1, nsoilmx
      IF (isoil.GT.99) THEN
         PRINT*, "Trop de couches"
         CALL abort
      ENDIF
      WRITE(str2,'(i2.2)') isoil
      CALL get_field('Tsoil'//str2,ftsoil(:,isoil),found)
      IF (.not.found) THEN
         PRINT*, "phyetat0: Le champ <Tsoil"//str2//"> est absent"
         PRINT*, "          Il prend donc la valeur de surface"
         DO i=1, klon
             ftsoil(i,isoil)=ftsol(i)
         ENDDO
      ENDIF
      ENDDO

! Lecture de albedo au sol:

      CALL get_field("ALBE", falbe,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ALBE> est absent'
         PRINT*, "phyetat0: Lecture echouee pour <ALBE>"
         CALL abort
      ELSE
         xmin = 1.0E+20
         xmax = -1.0E+20
         DO i = 1, klon
            xmin = MIN(falbe(i),xmin)
            xmax = MAX(falbe(i),xmax)
         ENDDO
         PRINT*,'Albedo du sol <ALBE>', xmin, xmax
      ENDIF

! Lecture rayonnement solaire au sol:

      CALL get_field("solsw",solsw,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <solsw> est absent'
         PRINT*, 'mis a zero'
         solsw = 0.
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(solsw(i),xmin)
         xmax = MAX(solsw(i),xmax)
      ENDDO
      PRINT*,'Rayonnement solaire au sol solsw:', xmin, xmax

! Lecture rayonnement IF au sol:

      CALL get_field("sollw",sollw,found)
      IF (.not.found) THEN
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

! Lecture derive des flux:

      CALL get_field("fder",fder,found)
      IF (.not.found) THEN
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

! Lecture du rayonnement net au sol:

      CALL get_field("RADS",radsol,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <RADS> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(radsol(i),xmin)
         xmax = MAX(radsol(i),xmax)
      ENDDO
      PRINT*,'Rayonnement net au sol radsol:', xmin, xmax

! Lecture de l'orographie sous-maille si ok_orodr:

      if(ok_orodr) then
     
      CALL get_field("ZMEA",zmea,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zmea(i),xmin)
         xmax = MAX(zmea(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zmea:', xmin, xmax

      CALL get_field("ZSTD",zstd,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zstd(i),xmin)
         xmax = MAX(zstd(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zstd:', xmin, xmax

      CALL get_field("ZSIG",zsig,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zsig(i),xmin)
         xmax = MAX(zsig(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zsig:', xmin, xmax

      CALL get_field("ZGAM",zgam,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zgam(i),xmin)
         xmax = MAX(zgam(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zgam:', xmin, xmax

      CALL get_field("ZTHE",zthe,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zthe(i),xmin)
         xmax = MAX(zthe(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zthe:', xmin, xmax

      CALL get_field("ZPIC",zpic,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZPIC> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zpic(i),xmin)
         xmax = MAX(zpic(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zpic:', xmin, xmax

      CALL get_field("ZVAL",zval,found)
      IF (.not.found) THEN
         PRINT*, 'phyetat0: Le champ <ZVAL> est absent'
         CALL abort
      ENDIF
      xmin = 1.0E+20
      xmax = -1.0E+20
      DO i = 1, klon
         xmin = MIN(zval(i),xmin)
         xmax = MAX(zval(i),xmax)
      ENDDO
      PRINT*,'OROGRAPHIE SOUS-MAILLE zval:', xmin, xmax

      else 
         zmea = 0.
         zstd = 0.
         zsig = 0.
         zgam = 0.
         zthe = 0.
         zpic = 0.
         zval = 0.

      endif   ! fin test sur ok_orodr

!     Par defaut on cree 2 bandes de methane au pole Nord et au pole Sud
!     (entre 75 et 85 degres de latitude) de 2 metres.
!     Les poles sont sec !
      resch4(1) = 0.    ! pole nord = 1 point
      DO i=2,klon
          if ((latitude_deg(i).ge.75..and.latitude_deg(i).le.85.).or.  &
              (latitude_deg(i).ge.-85.and.latitude_deg(i).le.-75.)) then
            resch4(i) = 2.
          else
            resch4(i) = 0.
          endif
      ENDDO
      resch4(klon) = 0.   ! pole sud = 1 point

      CALL get_field("RESCH4",resch4,found)
      IF (.not.found) THEN
         PRINT*, "phyetat0: Le champ <RESCH4> est absent"
         PRINT*, "Pas de reservoir de methane mais je continue..."
         PRINT*, "Pour info, je met 2 metres de methane sur 2 bandes"
         PRINT*, "comprises entre 75 et 85 degres de latitude dans  "
         PRINT*, "chaque hemisphere."         
      ENDIF

! Lecture de TANCIEN:

      ancien_ok = .TRUE.

      CALL get_field("TANCIEN",t_ancien,found)
      IF (.not.found) THEN
         PRINT*, "phyetat0: Le champ <TANCIEN> est absent"
         PRINT*, "Depart legerement fausse. Mais je continue"
         ancien_ok = .FALSE.
      ENDIF

! close file
call close_startphy

! do some more initializations
call init_iophy_new(latitude_deg,longitude_deg)

end subroutine phyetat0
