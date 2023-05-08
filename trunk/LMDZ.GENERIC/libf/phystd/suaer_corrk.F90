      subroutine suaer_corrk

      ! inputs
      use radinc_h,    only: L_NSPECTI,L_NSPECTV,nsizemax,naerkind
      use radcommon_h, only: blamv,blami,lamrefir,lamrefvis
      use datafile_mod, only: datadir, aerdir

      ! outputs
      use radcommon_h, only: QVISsQREF,omegavis,gvis,QIRsQREF,omegair,gir
      use radcommon_h, only: radiustab,nsize,tstellar
      use radcommon_h, only: qrefvis,qrefir,omegarefir !,omegarefvis
      use aerosol_mod
      use callkeys_mod, only: tplanet, optprop_back2lay_vis, optprop_back2lay_ir

      implicit none

!==================================================================
!     Purpose
!     -------
!     Initialize all radiative aerosol properties
!
!     Notes
!     -----
!     Reads the optical properties -> Mean  -> Variable assignment
!     (ASCII files)                  (see radcommon_h.F90)
!     wvl(nwvl)                      longsun
!     ep(nwvl)                       epav     QVISsQREF(L_NSPECTV)
!     omeg(nwvl)                     omegav   omegavis(L_NSPECTV)
!     gfactor(nwvl)                  gav      gvis(L_NSPECTV)
!     
!     Authors
!     ------- 
!     Richard Fournier (1996) Francois Forget (1996)
!     Frederic Hourdin
!     Jean-jacques morcrette *ECMWF*
!     MODIF Francois Forget (2000)
!     MODIF Franck Montmessin (add water ice)
!     MODIF J.-B. Madeleine 2008W27
!     - Optical properties read in ASCII files
!     - Add varying radius of the particules
!     - Add water ice clouds
!     MODIF R. Wordsworth (2009)
!     - generalisation to arbitrary spectral bands 
!
!==================================================================

!     Optical properties (read in external ASCII files)
      INTEGER,SAVE      :: nwvl  ! Number of wavelengths in
                                ! the domain (VIS or IR), read by master

!      REAL             :: solsir ! visible to infrared ratio
!                                ! (tau_0.67um/tau_9um)

      REAL, DIMENSION(:),&
      ALLOCATABLE, SAVE :: wvl  ! Wavelength axis, read by master
      REAL, DIMENSION(:),&
      ALLOCATABLE, SAVE :: radiusdyn ! Particle size axis, read by master

      REAL, DIMENSION(:,:),&
      ALLOCATABLE, SAVE :: ep,& ! Extinction coefficient Qext, read by master
      omeg,&                    ! Single Scattering Albedo, read by master
      gfactor                   ! Assymetry Factor, read by master

!     Local variables:

      INTEGER :: iaer           ! Aerosol index
      INTEGER :: idomain        ! Domain index (1=VIS,2=IR)
      INTEGER :: ilw            ! longwave index
      INTEGER :: isw            ! shortwave index
      INTEGER :: isize          ! Particle size index
      INTEGER :: jfile          ! ASCII file scan index
      INTEGER :: file_unit = 60
      LOGICAL :: file_ok, endwhile
      CHARACTER(LEN=132) :: scanline ! ASCII file scanning line

!     I/O  of "aerave" (subroutine that spectrally averages
!     the single scattering parameters)

      REAL lamref                      ! reference wavelengths
      REAL epref                       ! reference extinction ep

!      REAL epav(L_NSPECTI)            ! Average ep (= <Qext>/Qext(lamref) if epref=1)
!      REAL omegav(L_NSPECTI)          ! Average single scattering albedo
!      REAL gav(L_NSPECTI)             ! Average assymetry parameter

      REAL epavVI(L_NSPECTV)            ! Average ep (= <Qext>/Qext(lamref) if epref=1)
      REAL omegavVI(L_NSPECTV)          ! Average single scattering albedo
      REAL gavVI(L_NSPECTV)             ! Average assymetry parameter

      REAL epavIR(L_NSPECTI)            ! Average ep (= <Qext>/Qext(lamref) if epref=1)
      REAL omegavIR(L_NSPECTI)          ! Average single scattering albedo
      REAL gavIR(L_NSPECTI)             ! Average assymetry parameter
      
      logical forgetit                  ! use francois' data?
      integer iwvl

!     Local saved variables:

      CHARACTER(LEN=30), DIMENSION(naerkind,2), SAVE :: file_id
!$OMP THREADPRIVATE(file_id)
!---- Please indicate the names of the optical property files below
!     Please also choose the reference wavelengths of each aerosol      

      if (noaero) then
        print*, 'naerkind= 0'
      endif
      do iaer=1,naerkind
       if (iaer.eq.iaero_co2) then
        forgetit=.true.
          if (.not.noaero) then
              print*, 'naerkind= co2', iaer
          end if
!     visible
        if(forgetit)then
           file_id(iaer,1) = 'optprop_co2_vis_ff.dat' ! Francois' values
        else
           file_id(iaer,1) = 'optprop_co2ice_vis_n50.dat'
        endif
        lamrefvis(iaer)=1.5E-6   ! 1.5um OMEGA/MEx ???

!     infrared
        if(forgetit)then
           file_id(iaer,2) = 'optprop_co2_ir_ff.dat' ! Francois' values
        else
           file_id(iaer,2) = 'optprop_co2ice_ir_n50.dat'
        endif
        lamrefir(iaer)=12.1E-6   ! 825cm-1 TES/MGS ???
       endif ! CO2 aerosols  
! NOTE: these lamref's are currently for dust, not CO2 ice.
! JB thinks this shouldn't matter too much, but it needs to be 
! fixed / decided for the final version.

       if (iaer.eq.iaero_h2o) then
        print*, 'naerkind= h2o', iaer

!     visible
         file_id(iaer,1) = 'optprop_icevis_n50.dat'
         lamrefvis(iaer)=1.5E-6   ! 1.5um OMEGA/MEx
!     infrared
         file_id(iaer,2) = 'optprop_iceir_n50.dat'
         lamrefir(iaer)=12.1E-6   ! 825cm-1 TES/MGS
       endif

       if (iaer.eq.iaero_dust) then
        print*, 'naerkind= dust', iaer

!     visible
         file_id(iaer,1) = 'optprop_dustvis_n50.dat'
         !lamrefvis(iaer)=1.5E-6   ! 1.5um OMEGA/MEx
         lamrefvis(iaer)=0.67e-6
!     infrared
         file_id(iaer,2) = 'optprop_dustir_n50.dat'
         !lamrefir(iaer)=12.1E-6   ! 825cm-1 TES/MGS
         lamrefir(iaer)=9.3E-6
       endif 

       if (iaer.eq.iaero_h2so4) then
         print*, 'naerkind= h2so4', iaer

!     visible
         file_id(iaer,1) = 'optprop_h2so4vis_n20.dat'
         !lamrefir(iaer)=  doesn't exist?
         lamrefvis(iaer)=1.5E-6   ! no idea, must find
!     infrared
         file_id(iaer,2) = 'optprop_h2so4ir_n20.dat'
         !lamrefir(iaer)=12.1E-6   ! 825cm-1 TES/MGS
         lamrefir(iaer)=9.3E-6 ! no idea, must find
! added by LK
       endif

       if (iaer.eq.iaero_back2lay) then
         print*, 'naerkind= back2lay', iaer
         
!     visible
         file_id(iaer,1) = TRIM(optprop_back2lay_vis)
         lamrefvis(iaer)=0.8E-6  ! 
!     infrared
         file_id(iaer,2) = TRIM(optprop_back2lay_ir)
         lamrefir(iaer)=6.E-6  ! 
! added by SG
       endif
      
       if (iaer.eq.iaero_nh3) then
         print*, 'naerkind= nh3', iaer

!     visible
         file_id(iaer,1) = 'optprop_nh3ice_vis.dat'
         lamrefvis(iaer)=0.756E-6  ! 
!     infrared
         file_id(iaer,2) = 'optprop_nh3ice_ir.dat'
         lamrefir(iaer)=6.E-6  ! 
! added by SG
       endif

       if (iaer.eq.iaero_aurora) then
         print*, 'naerkind= aurora', iaer

!     visible
         file_id(iaer,1) = 'optprop_aurora_vis.dat'
         lamrefvis(iaer)=0.3E-6  ! 
!     infrared
         file_id(iaer,2) = 'optprop_aurora_ir.dat'
         lamrefir(iaer)=6.E-6  ! 
! added by SG
       endif
 
       
      enddo

!------------------------------------------------------------------

!     Initializations:

      radiustab(:,:,:) = 0.0
      gvis(:,:,:)      = 0.0
      omegavis(:,:,:)  = 0.0
      QVISsQREF(:,:,:) = 0.0
      gir(:,:,:)       = 0.0
      omegair(:,:,:)   = 0.0
      QIRsQREF(:,:,:)  = 0.0

      DO iaer = 1, naerkind     ! Loop on aerosol kind
         DO idomain = 1, 2      ! Loop on radiation domain (VIS or IR)
!==================================================================
!     1. READ OPTICAL PROPERTIES
!==================================================================

!     1.1 Open the ASCII file

!$OMP MASTER

            INQUIRE(FILE=TRIM(datadir)//'/'//TRIM(aerdir)//&
                    '/'//TRIM(file_id(iaer,idomain)),&
                    EXIST=file_ok)
            IF (file_ok) THEN
              OPEN(UNIT=file_unit,&
                   FILE=TRIM(datadir)//'/'//TRIM(aerdir)//&
                        '/'//TRIM(file_id(iaer,idomain)),&
                   FORM='formatted')
            ELSE
             ! In ye old days these files were stored in datadir;
             ! let's be retro-compatible
              INQUIRE(FILE=TRIM(datadir)//&
                      '/'//TRIM(file_id(iaer,idomain)),&
                      EXIST=file_ok)
              IF (file_ok) THEN
                OPEN(UNIT=file_unit,&
                     FILE=TRIM(datadir)//&
                          '/'//TRIM(file_id(iaer,idomain)),&
                     FORM='formatted')
              ENDIF              
            ENDIF
            IF(.NOT.file_ok) THEN
               write(*,*)'suaer_corrk: Problem opening ',&
               TRIM(file_id(iaer,idomain))
               write(*,*)'It should be in: ',TRIM(datadir)//'/'//TRIM(aerdir)
               write(*,*)'1) You can set this directory address ',&
               'in callphys.def with:'
               write(*,*)' datadir = /absolute/path/to/datagcm'
               write(*,*)'2) If ',&
              TRIM(file_id(iaer,idomain)),&
               ' is a LMD reference datafile, it'
               write(*,*)' can be obtained online at:'
               write(*,*)' http://www.lmd.jussieu.fr/',&
               '~lmdz/planets/LMDZ.GENERIC/datagcm/'
               CALL ABORT 
            ENDIF

!     1.2 Allocate the optical property table

            jfile = 1
            endwhile = .false.
            DO WHILE (.NOT.endwhile)
               READ(file_unit,*) scanline
               IF ((scanline(1:1) .ne. '#').and.&
               (scanline(1:1) .ne. ' ')) THEN
               BACKSPACE(file_unit)
               reading1_seq: SELECT CASE (jfile) ! ====================
            CASE(1) reading1_seq ! nwvl ----------------------------
               read(file_unit,*) nwvl
               jfile = jfile+1
            CASE(2) reading1_seq ! nsize ---------------------------
               read(file_unit,*) nsize(iaer,idomain)
               endwhile = .true.
            CASE DEFAULT reading1_seq ! ---------------------------- 
               WRITE(*,*) 'readoptprop: ',& 
               'Error while loading optical properties.' 
               CALL ABORT 
            END SELECT reading1_seq ! ==============================
         ENDIF
      ENDDO

      ALLOCATE(wvl(nwvl))       ! wvl
      ALLOCATE(radiusdyn(nsize(iaer,idomain))) ! radiusdyn
      ALLOCATE(ep(nwvl,nsize(iaer,idomain))) ! ep
      ALLOCATE(omeg(nwvl,nsize(iaer,idomain))) ! omeg 
      ALLOCATE(gfactor(nwvl,nsize(iaer,idomain))) ! g


!     1.3 Read the data

      jfile = 1
      endwhile = .false.
      DO WHILE (.NOT.endwhile)
         READ(file_unit,*) scanline
         IF ((scanline(1:1) .ne. '#').and.&
         (scanline(1:1) .ne. ' ')) THEN
         BACKSPACE(file_unit)
         reading2_seq: SELECT CASE (jfile) ! ====================
      CASE(1) reading2_seq      ! wvl -----------------------------
         read(file_unit,*) wvl
         jfile = jfile+1
      CASE(2) reading2_seq      ! radiusdyn -----------------------
         read(file_unit,*) radiusdyn
         jfile = jfile+1
      CASE(3) reading2_seq      ! ep ------------------------------
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
      CASE(4) reading2_seq      ! omeg ----------------------------
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
      CASE(5) reading2_seq      ! gfactor -------------------------
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

      jfile = jfile+1
      IF ((idomain.NE.iaero_co2).OR.(iaer.NE.iaero_co2)) THEN
         endwhile = .true.
      ENDIF
      CASE(6) reading2_seq
         endwhile = .true.
      CASE DEFAULT reading2_seq ! ----------------------------
         WRITE(*,*) 'readoptprop: ',&
         'Error while loading optical properties.'
         CALL ABORT
      END SELECT reading2_seq   ! ==============================
      ENDIF
      ENDDO

!     1.4 Close the file

      CLOSE(file_unit)

!     1.5 If Francois' data, convert wvl to metres
       if(iaer.eq.iaero_co2)then
         if(forgetit)then
         !   print*,'please sort out forgetit for naerkind>1'
            do iwvl=1,nwvl
               wvl(iwvl)=wvl(iwvl)*1.e-6
            enddo
         endif
      endif

!$OMP END MASTER
!$OMP BARRIER





!==================================================================
!     2. AVERAGED PROPERTIES AND VARIABLE ASSIGNMENTS
!==================================================================
      domain: SELECT CASE (idomain)
!==================================================================
      CASE(1) domain            !       VISIBLE DOMAIN (idomain=1)
!==================================================================

         lamref=lamrefvis(iaer)
         epref=1.E+0

         DO isize=1,nsize(iaer,idomain)

!     Save the particle sizes
            radiustab(iaer,idomain,isize)=radiusdyn(isize)

!     Averaged optical properties (GCM channels)

!            CALL aerave ( nwvl,&
!            wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
!            lamref,epref,tstellar,&
!            L_NSPECTV,blamv,epavVI,&
!            omegavVI,gavVI,QREFvis(iaer,isize))!,omegaREFir(iaer,isize))
            CALL aerave_new ( nwvl,&
            wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
            lamref,epref,tstellar,&
            L_NSPECTV,blamv,epavVI,&
            omegavVI,gavVI,QREFvis(iaer,isize),omegaREFir(iaer,isize))

!     Variable assignments (declared in radcommon)
            DO isw=1,L_NSPECTV
               QVISsQREF(isw,iaer,isize)=epavVI(isw)
               gvis(isw,iaer,isize)=gavVI(isw)
               omegavis(isw,iaer,isize)=omegavVI(isw)
            END DO

         ENDDO

!==================================================================
      CASE(2) domain            !      INFRARED DOMAIN (idomain=2)
!==================================================================


         DO isize=1,nsize(iaer,idomain) ! ----------------------------------

            lamref=lamrefir(iaer)
            epref=1.E+0

!     Save the particle sizes
            radiustab(iaer,idomain,isize)=radiusdyn(isize)

!     Averaged optical properties (GCM channels)

!     epav is <QIR>/Qext(lamrefir) since epref=1
!     Note: aerave also computes the extinction coefficient at
!     the reference wavelength. This is called QREFvis or QREFir
!     (not epref, which is a different parameter).
!     Reference wavelengths SHOULD BE defined for each aerosol in
!     radcommon_h.F90

!            CALL aerave ( nwvl,&
!            wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
!            lamref,epref,tplanet,&
!            L_NSPECTI,blami,epavIR,&
!            omegavIR,gavIR,QREFir(iaer,isize))!,omegaREFir(iaer,isize))
            CALL aerave_new ( nwvl,&
            wvl(:),ep(:,isize),omeg(:,isize),gfactor(:,isize),&
            lamref,epref,tplanet,&
            L_NSPECTI,blami,epavIR,&
            omegavIR,gavIR,QREFir(iaer,isize),omegaREFir(iaer,isize))


!     Variable assignments (declared in radcommon)
            DO ilw=1,L_NSPECTI
               QIRsQREF(ilw,iaer,isize)=epavIR(ilw)
               gir(ilw,iaer,isize)=gavIR(ilw)
               omegair(ilw,iaer,isize)=omegavIR(ilw)
            END DO


      ENDDO                     ! isize (particle size) -------------------------------------

      END SELECT domain

!========================================================================
!     3. Deallocate temporary variables that were read in the ASCII files
!========================================================================

!$OMP BARRIER
!$OMP MASTER
      IF (ALLOCATED(wvl)) DEALLOCATE(wvl)                 ! wvl
      IF (ALLOCATED(radiusdyn)) DEALLOCATE(radiusdyn)     ! radiusdyn
      IF (ALLOCATED(ep)) DEALLOCATE(ep)                   ! ep 
      IF (ALLOCATED(omeg)) DEALLOCATE(omeg)          	  ! omeg 
      IF (ALLOCATED(gfactor)) DEALLOCATE(gfactor)         ! g
!$OMP END MASTER
!$OMP BARRIER

      END DO                    ! Loop on iaer
      END DO                    ! Loop on idomain
!========================================================================

      RETURN



    END subroutine suaer_corrk
      
