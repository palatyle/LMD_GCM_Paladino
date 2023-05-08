SUBROUTINE suaer
use dimradmars_mod, only: longrefvis, longrefir, nsizemax, long1vis, &
                    long2vis, long3vis, long1ir, long2ir, long1co2, &
                    long2co2, nsun, nir,&
                    naerkind, name_iaer, &
                    iaer_dust_conrath,iaer_dust_doubleq,&
                    iaer_dust_submicron,iaer_h2o_ice,&
                    iaer_stormdust_doubleq,iaer_topdust_doubleq,&
                    file_id,radiustab, gvis, omegavis, &
                    QVISsQREF, gIR, omegaIR, &
                    QIRsQREF, QREFvis, QREFir, &
                    omegaREFvis, omegaREFir, &
                    nsize
use datafile_mod, only: datadir
IMPLICIT NONE
!==================================================================
!     Purpose.
!     --------
!     initialize yomaer, the common that contains the
!     radiative characteristics of the aerosols
!     
!     AUTHORS.
!     -------- 
!     Richard Fournier (1996) Francois Forget (1996)
!     Frederic Hourdin
!     Jean-jacques morcrette *ECMWF*
!     MODIF Francois Forget (2000)
!     MODIF Franck Montmessin (add water ice)
!     MODIF J.-B. Madeleine 2008W27
!       - Optical properties read in ASCII files
!       - Add varying radius of the particules
!
!    Summary.
!    --------
!
!    Read the optical properties -> Mean  -> Variable assignment
!                  (ASCII files)                  (see yomaer.h)
!    wvl(nwvl)                      longsun
!    ep(nwvl)                       epav     QVISsQREF(nsun)
!    omeg(nwvl)                     omegav   omegavis(nsun)
!    gfactor(nwvl)                  gav      gvis(nsun)
!    
!==================================================================

! Includes:

include "callkeys.h"

! Optical properties (read in external ASCII files)
INTEGER          :: nwvl     ! Number of wavelengths in
                             ! the domain (VIS or IR)
REAL, DIMENSION(:),&
  ALLOCATABLE, SAVE :: wvl       ! Wavelength axis
REAL, DIMENSION(:),&
  ALLOCATABLE, SAVE :: radiusdyn ! Particle size axis

REAL, DIMENSION(:,:),&
  ALLOCATABLE, SAVE :: ep,&    ! Extinction coefficient Qext
                       omeg,&  ! Single Scattering Albedo
                       gfactor ! Assymetry Factor

! Local variables:

INTEGER :: iaer                ! Aerosol index
INTEGER :: idomain             ! Domain index (1=VIS,2=IR)
INTEGER :: iir                 ! IR channel index
                               ! iir=1: 15um CO2 bands
                               ! iir=2 : CO2 band wings
                               ! iir=3 : 9 um band
                               ! iir=4 : Far IR
INTEGER :: isun                ! Solar band index
INTEGER :: isize               ! Particle size index
INTEGER :: jfile               ! ASCII file scan index
INTEGER :: file_unit = 60
LOGICAL :: file_ok, endwhile
CHARACTER(LEN=132) :: scanline ! ASCII file scanning line
INTEGER :: read_ok

! I/O  of "aerave" (subroutine averaging spectrally
!   sing.scat.parameters)

REAL tsun            ! Sun brightness temperature (for SW)
REAL tsol            ! Surface reference brightness temp (LW)
REAL longref         ! reference wavelengths
REAL longsun(nsun+1) ! solar band boundaries
REAL longir(nir+1)   ! IR band boundaries
REAL epref           ! reference extinction ep
                     ! at wavelength "longref"
REAL epav(nir)       ! average ep
                     ! (= <Qext>/Qext(longref) if epref=1)
REAL omegav(nir)     ! Average sing.scat.albedo
REAL gav(nir)        ! Average assymetry parameter

!==================================================================
!---- Please indicate the names of the optical property files below
!     Please also choose the reference wavelengths of each aerosol

      DO iaer = 1, naerkind ! Loop on aerosol kind
        aerkind: SELECT CASE (name_iaer(iaer))
!==================================================================
        CASE("dust_conrath") aerkind      ! Typical dust profile
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_dustvis_TM.dat'     !M.Wolff TM
!       file_id(iaer,1) = 'optprop_dustvis_clancy.dat' !Clancy-Lee
!       file_id(iaer,1) = 'optprop_dustvis_ockert.dat' !Ockert-Bell
!       Infrared domain:
        file_id(iaer,2) = 'optprop_dustir_TM.dat'      !M.Wolff TM
!       Toon-Forget + solsir=2 using Clancy-Lee
!       file_id(iaer,2) = 'optprop_dustir_clancy.dat'
!       Toon-Forget + solsir=2 using Ockert-Bell
!       file_id(iaer,2) = 'optprop_dustir_ockert.dat'
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6
!                     For dust: change readtesassim accordingly;
!       Reference wavelength in the infrared:
        longrefir(iaer)=dustrefir
!==================================================================
        CASE("dust_doubleq") aerkind! Two-moment scheme for dust
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_dustvis_TM_n50.dat' !T-Matrix
!       file_id(iaer,1) = 'optprop_dustvis_n50.dat'    !Mie
!       Infrared domain:
        file_id(iaer,2) = 'optprop_dustir_n50.dat'     !Mie
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6
!       If not equal to 0.67e-6 -> change readtesassim accordingly;
!       Reference wavelength in the infrared:
        longrefir(iaer)=dustrefir
!==================================================================
        CASE("dust_submicron") aerkind   ! Small dust population
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_dustvis_01um_TM.dat' !M.Wolff
!       Infrared domain:
        file_id(iaer,2) = 'optprop_dustir_01um_TM.dat'  !M.Wolff
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6
!       If not equal to 0.67e-6 -> change readtesassim accordingly;
!       Reference wavelength in the infrared:
        longrefir(iaer)=dustrefir
!==================================================================
        CASE("h2o_ice") aerkind             ! Water ice crystals
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_icevis_n30.dat' !Warren
!       file_id(iaer,1) = 'optprop_icevis.dat'     !Warren
!       Infrared domain:
        file_id(iaer,2) = 'optprop_iceir_n30.dat'  !Warren
!       file_id(iaer,2) = 'optprop_iceir.dat'      !Warren
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6  ! 1.5um OMEGA/MEx
!       Reference wavelength in the infrared:
        longrefir(iaer)=12.1E-6  ! 825cm-1 TES/MGS
!==================================================================
        CASE("stormdust_doubleq") aerkind   ! Two-moment scheme for stormdust - radiative properties
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_dustvis_TM_n50.dat' !T-Matrix
!       Infrared domain:
        file_id(iaer,2) = 'optprop_dustir_n50.dat'     !Mie
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6
!       If not equal to 0.67e-6 -> change readtesassim accordingly;
!       Reference wavelength in the infrared:
        longrefir(iaer)=dustrefir
!==================================================================
        CASE("topdust_doubleq") aerkind   ! Two-moment scheme for topdust - radiative properties
!==================================================================
!       Visible domain:
        file_id(iaer,1) = 'optprop_dustvis_TM_n50.dat' !T-Matrix
!       Infrared domain:
        file_id(iaer,2) = 'optprop_dustir_n50.dat'     !Mie
!       Reference wavelength in the visible:
        longrefvis(iaer)=0.67E-6
!       If not equal to 0.67e-6 -> change readtesassim accordingly;
!       Reference wavelength in the infrared:
        longrefir(iaer)=dustrefir
!==================================================================
        END SELECT aerkind
!==================================================================
        WRITE(*,*) "Scatterer: ",trim(name_iaer(iaer))
        WRITE(*,*) "  corresponding files: "
        WRITE(*,*) "VIS: ",trim(file_id(iaer,1))
        WRITE(*,*) "IR : ",trim(file_id(iaer,2))
!==================================================================
      ENDDO ! iaer (loop on aerosol kind)

! Initializations:

radiustab(1:naerkind,1:2,1:nsizemax)=0

gvis(1:nsun,1:naerkind,1:nsizemax)=0
omegavis(1:nsun,1:naerkind,1:nsizemax)=0
QVISsQREF(1:nsun,1:naerkind,1:nsizemax)=0

gIR(1:nir,1:naerkind,1:nsizemax)=0
omegaIR(1:nir,1:naerkind,1:nsizemax)=0
QIRsQREF(1:nir,1:naerkind,1:nsizemax)=0

QREFvis(1:naerkind,1:nsizemax)=0
QREFir(1:naerkind,1:nsizemax)=0
omegaREFvis(1:naerkind,1:nsizemax)=0
omegaREFir(1:naerkind,1:nsizemax)=0

DO iaer = 1, naerkind ! Loop on aerosol kind
  DO idomain = 1, 2   ! Loop on radiation domain (VIS or IR)
!==================================================================
! 1. READ OPTICAL PROPERTIES
!==================================================================

!       1.1 Open the ASCII file

INQUIRE(FILE=TRIM(datadir)//&
  '/'//TRIM(file_id(iaer,idomain)),&
  EXIST=file_ok)
IF(.NOT.file_ok) THEN
  write(*,*)'Problem opening ',&
    TRIM(file_id(iaer,idomain))
  write(*,*)'It should be in: ',&
    TRIM(datadir)
  write(*,*)'1) You can change this directory address in callfis.def with'
  write(*,*)'   datadir=/path/to/datafiles'
  write(*,*)'2) If ',&
    TRIM(file_id(iaer,idomain)),&
    ' is a LMD reference datafile, it'
  write(*,*)' can be obtained online on:'
  write(*,*)' http://www.lmd.jussieu.fr/',&
    '~lmdz/planets/mars/datadir'
  write(*,*)'3) If the name of the file is wrong, you can'
  write(*,*)' change it in file phymars/suaer.F90. Just'
  write(*,*)' modify the variable called file_id.'
  CALL ABORT 
ENDIF
OPEN(UNIT=file_unit,&
  FILE=TRIM(datadir)//&
  '/'//TRIM(file_id(iaer,idomain)),&
  FORM='formatted')

!       1.2 Allocate the optical property table

jfile = 1
endwhile = .false.
DO WHILE (.NOT.endwhile)
  READ(file_unit,*,iostat=read_ok) scanline
  if (read_ok.ne.0) then
    write(*,*)' readoptprop: Error reading file',&
    TRIM(datadir)//&
    '/'//TRIM(file_id(iaer,idomain))
    call abort
  endif
  IF ((scanline(1:1) .ne. '#').and.&
    (scanline(1:1) .ne. ' ')) THEN
    BACKSPACE(file_unit)
    reading1_seq: SELECT CASE (jfile) ! ====================
    CASE(1) reading1_seq ! nwvl ----------------------------
        read(file_unit,*,iostat=read_ok) nwvl
        if (read_ok.ne.0) then
          write(*,*)' readoptprop: Error while reading line:',&
          trim(scanline)
          write(*,*)'   of file',&
          TRIM(datadir)//&
          '/'//TRIM(file_id(iaer,idomain))
          call abort
        endif
        jfile = jfile+1
    CASE(2) reading1_seq ! nsize ---------------------------
        read(file_unit,*,iostat=read_ok) nsize(iaer,idomain)
        if (read_ok.ne.0) then
          write(*,*)' readoptprop: Error while reading line:',&
          trim(scanline)
          write(*,*)'   of file',&
          TRIM(datadir)//&
          '/'//TRIM(file_id(iaer,idomain))
          call abort
        endif
        endwhile = .true.
    CASE DEFAULT reading1_seq ! ---------------------------- 
        WRITE(*,*) 'readoptprop: ',& 
          'Error while loading optical properties.' 
        CALL ABORT 
    END SELECT reading1_seq ! ==============================
  ENDIF
ENDDO

ALLOCATE(wvl(nwvl))                            ! wvl
ALLOCATE(radiusdyn(nsize(iaer,idomain)))       ! radiusdyn
ALLOCATE(ep(nwvl,nsize(iaer,idomain)))         ! ep
ALLOCATE(omeg(nwvl,nsize(iaer,idomain)))       ! omeg 
ALLOCATE(gfactor(nwvl,nsize(iaer,idomain)))    ! g

!       1.3 Read the data

jfile = 1
endwhile = .false.
DO WHILE (.NOT.endwhile)
   READ(file_unit,*) scanline
  IF ((scanline(1:1) .ne. '#').and.&
    (scanline(1:1) .ne. ' ')) THEN
    BACKSPACE(file_unit)
    reading2_seq: SELECT CASE (jfile) ! ====================
    CASE(1) reading2_seq ! wvl -----------------------------
        read(file_unit,*) wvl
        jfile = jfile+1
    CASE(2) reading2_seq ! radiusdyn -----------------------
        read(file_unit,*) radiusdyn
        jfile = jfile+1
    CASE(3) reading2_seq ! ep ------------------------------
        isize = 1
        DO WHILE (isize .le. nsize(iaer,idomain))
          READ(file_unit,*) scanline
          IF ((scanline(1:1) .ne. '#').and.&
            (scanline(1:1) .ne. ' ')) THEN
          BACKSPACE(file_unit)
          read(file_unit,*) ep(:,isize)
          isize = isize + 1
          ENDIF
        ENDDO
        jfile = jfile+1
    CASE(4) reading2_seq ! omeg ----------------------------
        isize = 1
        DO WHILE (isize .le. nsize(iaer,idomain))
          READ(file_unit,*) scanline
          IF ((scanline(1:1) .ne. '#').and.&
            (scanline(1:1) .ne. ' ')) THEN
          BACKSPACE(file_unit)
          read(file_unit,*) omeg(:,isize)
          isize = isize + 1
          ENDIF
        ENDDO
        jfile = jfile+1
    CASE(5) reading2_seq ! gfactor -------------------------
        isize = 1
        DO WHILE (isize .le. nsize(iaer,idomain))
          READ(file_unit,*) scanline
          IF ((scanline(1:1) .ne. '#').and.&
            (scanline(1:1) .ne. ' ')) THEN
          BACKSPACE(file_unit)
          read(file_unit,*) gfactor(:,isize)
          isize = isize + 1
          ENDIF
        ENDDO
        endwhile = .true.
    CASE DEFAULT reading2_seq ! ----------------------------
        WRITE(*,*) 'suaer.F90: ',&
          'Error while loading optical properties.'
        CALL ABORT
    END SELECT reading2_seq ! ==============================
  ENDIF
ENDDO

!       1.4 Close the file

CLOSE(file_unit)

!==================================================================
! 2. AVERAGED PROPERTIES AND VARIABLE ASSIGNMENTS
!==================================================================
domain: SELECT CASE (idomain)
!==================================================================
CASE(1) domain !                   VISIBLE DOMAIN (idomain=1)
!==================================================================

! 2.1 Parameters
  tsun=6000.E+0
  longsun(1)=long1vis
  longsun(2)=long2vis
  longsun(3)=long3vis
  longref=longrefvis(iaer)
  epref=1.E+0

DO isize=1,nsize(iaer,idomain)
! test that there is enough room to store the data
 if (isize.gt.nsizemax) then
   write(*,*) "suaer: Error ! nsizemax is too small!"
   write(*,*) "       nsizemax=",nsizemax
   write(*,*) "       you must increase the value of nsizemax"
   write(*,*) "       in dimradmars_mod !"
   stop
 endif
! ------------------------------------------------
! 2.2 Save the particle sizes
  radiustab(iaer,idomain,isize)=radiusdyn(isize)
! 2.3 Averaged optical properties (GCM channels)
! Notice: Aerave also computes the extinction coefficient and
!   single scattering albedo at reference wavelength
!   (called QREFvis and OMEGAREFvis, same in the IR,
!   and not epref, which is a different parameter);
!   Reference wavelengths are defined for each aerosol in
!   dimradmars_mod.

  CALL aerave ( nwvl,&
       wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
       longref,epref,tsun,&
       nsun,longsun, epav,omegav,gav,&
       QREFvis(iaer,isize),omegaREFvis(iaer,isize) )
! 2.4 Variable assignements (declared by yomaer.h)
  DO isun=1,nsun
    QVISsQREF(isun,iaer,isize)=epav(isun)
    gvis(isun,iaer,isize)=gav(isun)
    omegavis(isun,iaer,isize)=omegav(isun)
  END DO
! 2.5 Output display
!  WRITE(*,*) 'Les donnees spectrales :'
!  WRITE(*,*) 'Solaire (SW) ---->'
!  WRITE(*,*) 'Aerosol number: ', iaer
!  WRITE(*,*) 'Rayon aerosol: ', radiustab(iaer,idomain,isize)
!  WRITE(*,*) '<Qext>/Qext(longrefvis) ; omega ; g'
!  DO isun=1,nsun
!    WRITE(*,*) QVISsQREF(isun,iaer,isize),&
!         omegavis(isun,iaer,isize),&
!         gvis(isun,iaer,isize)
!  ENDDO
!  WRITE(*,*) 'QREFvis(',iaer,isize,') = ',QREFvis(iaer,isize)
!  WRITE(*,*) 'omegaREFvis(',iaer,isize,') = ',&
!                                      omegaREFvis(iaer,isize)
! ------------------------------------------------
ENDDO

!==================================================================
CASE(2) domain !                  INFRARED DOMAIN (idomain=2)
!==================================================================

DO isize=1,nsize(iaer,idomain) ! ----------------------------------

! 2.1 solsir is not used anymore; division of Qext(IR) by solsir
!     has to be done in the input ASCII files (if necessary).

! 2.2 Save the particle sizes
  radiustab(iaer,idomain,isize)=radiusdyn(isize)

! 2.3 Parameters

  tsol=215.D+0
  longir(1)=long1ir
  longir(2)=long1co2
  longir(3)=long2co2
  longir(4)=long2ir
  longref=longrefir(iaer)
  epref=1.E+0

! 2.4 Averaged optical properties (GCM channels)
!           epav is <QIR>/Qext(longrefir) since epref=1
! Notice: Aerave also Computes the extinction coefficient at
!   reference wavelength (called QREFvis or QREFir,
!   and not epref, which is a different parameter);
!   Reference wavelengths are defined for each aerosol in
!   dimradmar_mod.

  CALL aerave ( nwvl,&
       wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
       longref,epref,tsol,&
       nir-1,longir,epav,omegav,gav,&
       QREFir(iaer,isize),omegaREFir(iaer,isize) )
!  WRITE(*,*) 'QREFir(',iaer,isize,') = ',QREFir(iaer,isize)
!  WRITE(*,*) 'omegaREFir(',iaer,isize,') = ',&
!                                      omegaREFir(iaer,isize)

! 2.5 Computing  <QIR>/Qext(longrefvis)

  DO iir=1,nir-1
!    WRITE(*,*) 'QIRsQREFir Channel ',iir,': ',epav(iir)
    epav(iir)=  epav(iir) * QREFir(iaer,isize) / &
                            QREFvis(iaer,isize)
  ENDDO
!  WRITE(*,*) 'Aerosol number', iaer
!  WRITE(*,*) 'Particle size: ',radiustab(iaer,idomain,isize)
!  WRITE(*,*) 'Rapport Solaire/IR:',&
!             QREFvis(iaer,isize) / QREFir(iaer,isize)

! 2.6 Variable assignements
!           (variables are declared by yomaer.h)

!         Single scattering properties
!           in each of the "nir" bands
!           (cf. dimradmars_mod)

! iir=1 : central 15um CO2 bands   
  QIRsQREF(1,iaer,isize)=epav(2)
  omegaIR(1,iaer,isize)=omegav(2)
  gIR(1,iaer,isize)=gav(2)

! iir=2 : CO2 band wings
!           (same properties than for central part)
  QIRsQREF(2,iaer,isize)=epav(2)
  omegaIR(2,iaer,isize)=omegav(2)
  gIR(2,iaer,isize)=gav(2)

! iir=3 : 9 um band [long1ir - long1co2]
  QIRsQREF(3,iaer,isize)=epav(1)
  omegaIR(3,iaer,isize)=omegav(1)
  gIR(3,iaer,isize)=gav(1)

! iir=4 : Far IR    [long2co2 - long2ir]
  QIRsQREF(4,iaer,isize)=epav(3)
  omegaIR(4,iaer,isize)=omegav(3)
  gIR(4,iaer,isize)=gav(3)

! 2.7 Output display

!  WRITE(*,*) 'AEROSOL PROPERTIES: Number ',iaer
!  WRITE(*,*) 'Thermal IR (LW) ---->'
!  WRITE(*,*) 'Particle size: ',radiustab(iaer,idomain,isize)
!  WRITE(*,*) '<Qext>/Qext(longrefvis) ; omega ; g'
!  DO iir=1,nir
!    WRITE(*,*) QIRsQREF(iir,iaer,isize),omegaIR(iir,iaer,isize),&
!         gIR(iir,iaer,isize)
!  ENDDO
!  WRITE(*,*) 'CO2: <Qabs>/Qext(longrefvis) = ',&
!       QIRsQREF(1,iaer,isize)*(1-omegaIR(1,iaer,isize))
!  WRITE(*,*) ''

ENDDO ! isize (particle size) -------------------------------------

END SELECT domain
!==================================================================
! 3. Deallocate temporary variables (read in ASCII files)
!==================================================================

DEALLOCATE(wvl)        ! wvl
DEALLOCATE(radiusdyn)  ! radiusdyn
DEALLOCATE(ep)         ! ep 
DEALLOCATE(omeg)       ! omeg 
DEALLOCATE(gfactor)    ! g

  END DO ! Loop on iaer
END DO   ! Loop on idomain
!==================================================================
RETURN
END
