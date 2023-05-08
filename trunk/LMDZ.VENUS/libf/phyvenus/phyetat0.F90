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
      use geometry_mod, only: longitude_deg, latitude_deg
      USE time_phylmdz_mod, only: itau_phy, raz_date, pdtphys
      USE ioipsl_getin_p_mod, only: getin_p

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
LOGICAL :: found
REAL    :: tab_cntrl(length)
integer :: i,isoil
CHARACTER(len=2) :: str2
REAL :: lon_startphy(klon), lat_startphy(klon)
REAL :: surface_albedo

! les variables globales lues dans le fichier restart

! open physics initial state file:
if (startphy_file) then
  call open_startphy(fichnom)
endif

!
! Load control parameters:
!
IF (startphy_file) THEN
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

  itau_phy = tab_cntrl(15)

! Attention si raz_date est active :
! il faut remettre a zero itau_phy apres phyetat0 !
  IF (raz_date.eq.1) THEN
    itau_phy=0
  ENDIF

ELSE
  tabcntr0(:)=1 ! dummy initialization
  ! Initialize parameter or get values from def files
  dtime=pdtphys
  radpas=1
  itau_phy=0
ENDIF ! of IF (startphy_file)

IF (startphy_file) THEN
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
ENDIF ! of IF (startphy_file)

! read in other variables here ...

IF (startphy_file) THEN
  ! Load surface temperature:
  CALL get_field("TS",ftsol(:),found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <TS> est absent'
    PRINT*, "phyetat0: Lecture echouee pour <TS>"
    CALL abort
  ELSE
    PRINT*, 'phyetat0: Le champ <TS> est present'
    PRINT*,'Temperature du sol <TS>', minval(ftsol), maxval(ftsol)
  ENDIF
ELSE
  ! Dummy initialization, but in fact this is later handled in physiq
  ftsol(:)=0
ENDIF ! of IF (startphy_file)

IF (startphy_file) THEN
  ! Load sub-surface temperatures:
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
ELSE
  ! Dummy initialization, but in fact this is later handled in physiq
  ftsoil(:,:)=0
ENDIF ! of IF (startphy_file)

IF (startphy_file) THEN
  ! Load surface albedo:
  CALL get_field("ALBE", falbe,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ALBE> est absent'
    PRINT*, "phyetat0: Lecture echouee pour <ALBE>"
    CALL abort
  ENDIF
ELSE
  ! Dummy initialization: read value from def file 
  surface_albedo=0.5 ! default
  CALL getin_p("surface_albedo",surface_albedo)
  falbe(:)=surface_albedo
ENDIF ! of IF (startphy_file)
PRINT*,'Albedo du sol <ALBE>', minval(falbe), maxval(falbe)

IF (startphy_file) THEN
  ! Lecture rayonnement solaire au sol:
  CALL get_field("solsw",solsw,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <solsw> est absent'
    PRINT*, 'mis a zero'
    solsw = 0.
  ENDIF
ELSE
  ! Dummy initialization
  solsw(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Rayonnement solaire au sol solsw:', minval(solsw), maxval(solsw)

IF (startphy_file) THEN
  ! Lecture rayonnement IR au sol:
  CALL get_field("sollw",sollw,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <sollw> est absent'
    PRINT*, 'mis a zero'
    sollw = 0.
  ENDIF
ELSE
  ! Dummy initialization
  sollw(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Rayonnement IR au sol sollw:', minval(sollw), maxval(solsw)

IF (startphy_file) THEN
  ! Lecture derive des flux:
  CALL get_field("fder",fder,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <fder> est absent'
    PRINT*, 'mis a zero'
    fder = 0.
  ENDIF
ELSE
  ! Dummy initialization
  fder(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Derive des flux fder:', minval(fder), maxval(fder)

IF (startphy_file) THEN
  ! Lecture derive flux IR:
  CALL get_field("dlw",dlw,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <dlw> est absent'
    PRINT*, 'mis a zero'
    dlw = 0.
  ENDIF
ELSE
  ! Dummy initialization
  dlw(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Derive flux IR dlw:', minval(dlw), maxval(dlw)

IF (startphy_file) THEN
  ! Lecture rayonnement IR vers le bas au sol:
  CALL get_field("sollwdown",sollwdown,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <sollwdown> est absent'
    PRINT*, 'mis a zero'
    sollwdown = 0.
  ENDIF
ELSE
  ! Dummy initialization
  sollwdown(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Flux IR vers le bas au sol sollwdown:', minval(sollwdown), maxval(sollwdown)

IF (startphy_file) THEN
  ! Lecture du rayonnement net au sol:
  CALL get_field("RADS",radsol,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <RADS> est absent'
    CALL abort
  ENDIF
ELSE
  ! Dummy initialization
  radsol(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'Rayonnement net au sol radsol:', minval(radsol), maxval(radsol)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZMEA",zmea,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZMEA> est absent'
    PRINT*, 'mis a zero'
    zmea=0.
  ENDIF
ELSE
  zmea(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zmea:', minval(zmea), maxval(zmea)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZSTD",zstd,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZSTD> est absent'
    PRINT*, 'mis a zero'
    zstd=0.
  ENDIF
ELSE
  zstd(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zstd:', minval(zstd), maxval(zstd)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZSIG",zsig,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZSIG> est absent'
    PRINT*, 'mis a zero'
    zsig=0.
  ENDIF
ELSE
  zsig(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zsig:', minval(zsig), maxval(zsig)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZGAM",zgam,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZGAM> est absent'
    PRINT*, 'mis a zero'
    zgam=0.
  ENDIF
ELSE
  zgam(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zgam:', minval(zgam), maxval(zgam)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZTHE",zthe,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZTHE> est absent'
    PRINT*, 'mis a zero'
    zthe=0.
  ENDIF
ELSE
  zthe(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zthe:', minval(zthe), maxval(zthe)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZPIC",zpic,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZPIC> est absent'
    PRINT*, 'mis a zero'
    zpic=0.
  ENDIF
ELSE
  zpic(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zpic:', minval(zpic), maxval(zpic)

IF (startphy_file) THEN
  ! Load sub-grid scale orography parameters:
  CALL get_field("ZVAL",zval,found)
  IF (.not.found) THEN
    PRINT*, 'phyetat0: Le champ <ZVAL> est absent'
    PRINT*, 'mis a zero'
    zval=0.
  ENDIF
ELSE
  zval(:)=0
ENDIF ! of IF (startphy_file)
PRINT*,'OROGRAPHIE SOUS-MAILLE zval:', minval(zval), maxval(zval)

IF (startphy_file) THEN
  ! Lecture de TANCIEN:
  ancien_ok = .TRUE.

  CALL get_field("TANCIEN",t_ancien,found)
  IF (.not.found) THEN
    PRINT*, "phyetat0: Le champ <TANCIEN> est absent"
    PRINT*, "Depart legerement fausse. Mais je continue"
    ancien_ok = .FALSE.
  ENDIF
ELSE
  ancien_ok=.false.
ENDIF

! close file
IF (startphy_file) call close_startphy

! do some more initializations
call init_iophy_new(latitude_deg,longitude_deg)

end subroutine phyetat0
