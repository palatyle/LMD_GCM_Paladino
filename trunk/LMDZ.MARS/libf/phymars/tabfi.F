c=======================================================================
      SUBROUTINE tabfi(nid,Lmodif,tab0,day_ini,lmax,p_rad,
     .                 p_omeg,p_g,p_mugaz,p_daysec,time)
c=======================================================================
c
c   C. Hourdin 15/11/96
c
c   Object:        Lecture du tab_cntrl physique dans un fichier 
c   ------            et initialisation des constantes physiques
c
c   Arguments:
c   ----------
c
c     Inputs:
c     ------
c
c      - nid:    unitne logique du fichier ou on va lire le tab_cntrl    
c                      (ouvert dans le programme appellant) 
c
c                 si nid=0:
c                       pas de lecture du tab_cntrl mais
c                       Valeurs par default des constantes physiques
c        
c      - tab0:    Offset de tab_cntrl a partir duquel sont ranges 
c                  les parametres physiques (50 pour start_archive)
c
c      - Lmodif:  si on souhaite modifier les constantes  Lmodif = 1 = TRUE
c
c
c     Outputs:
c     --------
c
c      - day_ini: tab_cntrl(tab0+3) (Dans les cas ou l'on souhaite
c                              comparer avec le day_ini dynamique)
c
c      - lmax:    tab_cntrl(tab0+2) (pour test avec nlayer)
c
c      - p_rad
c      - p_omeg   !
c      - p_g      ! Constantes physiques ayant des 
c      - p_mugaz  ! homonymes dynamiques
c      - p_daysec !
c
c=======================================================================
! to use  'getin_p'
      use ioipsl_getin_p_mod, only: getin_p

      use comsoil_h, only: volcapa ! soil volumetric heat capacity
      use surfdat_h, only: z0_default, emissiv, emisice, albedice,
     &                     iceradius, dtemisice, iceradius
      use dimradmars_mod, only: tauvis
      use iostart, only: get_var
      use mod_phys_lmdz_para, only: is_parallel
      use comcstfi_h, only: g, mugaz, omeg, rad, rcp
      use time_phylmdz_mod, only: daysec, dtphys
      use planete_h, only: aphelie, emin_turb, lmixmin, obliquit,
     &                     peri_day, periheli, year_day
      implicit none
 
      include "netcdf.inc"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c Arguments
c ---------
      INTEGER,INTENT(IN) :: nid,tab0
      INTEGER*4,INTENT(OUT) :: day_ini
      INTEGER,INTENT(IN) :: Lmodif
      INTEGER,INTENT(OUT) :: lmax
      REAL,INTENT(OUT) :: p_rad,p_omeg,p_g,p_mugaz,p_daysec,time

c Variables
c ---------
      INTEGER :: nvarid
      REAL :: peri_ls
      INTEGER length
      parameter (length = 100)
      REAL tab_cntrl(length) ! array in which are stored the run's parameters
      INTEGER  ierr
      INTEGER size
      CHARACTER modif*20
      LOGICAL :: found
      CHARACTER(len=5) :: modname="tabfi"

      write(*,*)"tabfi: nid=",nid," tab0=",tab0," Lmodif=",Lmodif
      
c-----------------------------------------------------------------------
c  Initialization of various physical constants to defaut values (nid = 0 case)
c-----------------------------------------------------------------------
      IF (nid.eq.0) then

      ! to avoid further issues with writing 
      tab_cntrl(:)=0
      lmax=0 
      day_ini=0
      time = 0
 
c Reference pressure
c-------------------------------------
c     pressrf = 670.            ! Pression de reference (Pa) ~650

c Infos about Mars for the dynamics and physics
c----------------------------------------------------------
      rad=3397200.          ! radius of Mars (m)  ~3397200 m
      daysec=88775.         ! length of a sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec)       ! rotation rate  (rad.s-1)
      g=3.72                ! gravity (m.s-2) ~3.72
      mugaz=43.49           ! Molar mass of the atmosphere (g.mol-1) ~43.49
      rcp=.256793         ! = r/cp  ~0.256793

c Informations about Mars, only for physics
c-----------------------------------------------------
      year_day = 669.       !Modif FH: length of year (sols) ~668.6
      periheli = 206.66         ! min. Sun-Mars distance (Mkm) ~206.66
      aphelie = 249.22          ! max. Sun-Mars distance (Mkm) ~249.22
      peri_day =  485.    ! date of perihelion (sols since N. spring)
      obliquit = 25.19         ! Obliquity of the planet (deg) ~25.19

c Boundary layer and turbulence
c----------------------------
      z0_default =  1.e-2       ! surface roughness (m) ~0.01
      emin_turb = 1.e-6         ! minimal energy ~1.e-8
      lmixmin = 30              ! mixing length ~100

c Optical properties of polar caps and ground emissivity
c-----------------------------------------------------
      emissiv=.95               ! Emissivity of martian soil ~.95
      emisice(1)=0.95           ! Emissivity of northern cap
      emisice(2)=0.95           ! Emissivity of southern cap
      albedice(1)=0.65          ! Albedo of northern cap
      albedice(2)=0.65          ! Albedo of southern cap
      iceradius(1) = 100.e-6    ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6    ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 0.4   ! time scale for snow metamorphism (north)
      dtemisice(2) = 0.4   ! time scale for snow metamorphism (south)

c dust aerosol properties
c---------------------------------
      tauvis= 0.2          ! mean visible optical depth

c  Ancien code radiatif (non utilise avec le code d'apres 03/96)
c---------------------------------------------------------------
c     tauir= 0.  ! .2  ratio (mean IR opt.depth)/Visible
c     scatalb=0. ! .86 scaterring albedo visible (~.86)
c     asfact=0.  ! .79 assymetrie factor visible   (~.79)
c     day0 = 0   ! = 0 en general !!!

c soil properties
      volcapa = 1.e6 ! soil volumetric heat capacity (in comsoil.h)
      ELSE
c-----------------------------------------------------------------------
c  Initialization of physical constants by reading array tab_cntrl(:)
c		which contains these parameters	(nid != 0 case)
c-----------------------------------------------------------------------
c Read 'controle' array
c
       call get_var("controle",tab_cntrl,found)
       if (.not.found) then
         call abort_physic(modname, 
     &        "tabfi: Failed reading <controle> array",1)
       else
         write(*,*)'tabfi: tab_cntrl',tab_cntrl
       endif
c
c  Initialization of some physical constants
c informations on physics grid
      lmax = nint(tab_cntrl(tab0+2))
      day_ini = tab_cntrl(tab0+3)
      time = tab_cntrl(tab0+4)
      write (*,*) 'IN tabfi day_ini=',day_ini
c Informations about planet Mars for dynamics and physics
      rad = tab_cntrl(tab0+5)
      omeg = tab_cntrl(tab0+6)
      g = tab_cntrl(tab0+7)
      mugaz = tab_cntrl(tab0+8)
      rcp = tab_cntrl(tab0+9)
      daysec = tab_cntrl(tab0+10)
      dtphys = tab_cntrl(tab0+11)
c Informations about planet Mars for the physics only
      year_day = tab_cntrl(tab0+14)
      periheli = tab_cntrl(tab0+15)
      aphelie = tab_cntrl(tab0+16)
      peri_day = tab_cntrl(tab0+17)
      obliquit = tab_cntrl(tab0+18)
c boundary layer and turbeulence
      z0_default = tab_cntrl(tab0+19)
      lmixmin = tab_cntrl(tab0+20)
      emin_turb = tab_cntrl(tab0+21)
c optical properties of polar caps and ground emissivity
      albedice(1)= tab_cntrl(tab0+22)
      albedice(2)= tab_cntrl(tab0+23)
      emisice(1) = tab_cntrl(tab0+24)
      emisice(2) = tab_cntrl(tab0+25)
      emissiv    = tab_cntrl(tab0+26)
      tauvis     = tab_cntrl(tab0+27)  ! dust opt. depth vis.
      iceradius(1)= tab_cntrl(tab0+31) ! mean scat radius of CO2 snow (north)
      iceradius(2)= tab_cntrl(tab0+32) ! mean scat radius of CO2 snow (south)
      dtemisice(1)= tab_cntrl(tab0+33) !time scale for snow metamorphism (north)
      dtemisice(2)= tab_cntrl(tab0+34) !time scale for snow metamorphism (south)
c soil properties
      volcapa = tab_cntrl(tab0+35) ! volumetric heat capacity
c-----------------------------------------------------------------------
c	Save some constants for later use (as routine arguments)
c-----------------------------------------------------------------------
      p_omeg = omeg
      p_g = g
      p_mugaz = mugaz
      p_daysec = daysec
      p_rad=rad

      ENDIF    ! end of (nid = 0) 

c-----------------------------------------------------------------------
c	Write physical constants to output before modifying them
c-----------------------------------------------------------------------
 
   6  FORMAT(a20,e15.6,e15.6)
   5  FORMAT(a20,f12.2,f12.2)
 
      write(*,*) '*****************************************************'
      write(*,*) 'Reading tab_cntrl when calling tabfi before changes'
      write(*,*) '*****************************************************'
      write(*,5) '(1)        = ngrid?',tab_cntrl(tab0+1)
      write(*,5) '(2)            lmax',tab_cntrl(tab0+2),real(lmax)
      write(*,5) '(3)         day_ini',tab_cntrl(tab0+3),real(day_ini)
      write(*,5) '(5)             rad',tab_cntrl(tab0+5),rad
      write(*,5) '(10)         daysec',tab_cntrl(tab0+10),daysec
      write(*,6) '(6)            omeg',tab_cntrl(tab0+6),omeg
      write(*,5) '(7)               g',tab_cntrl(tab0+7),g
      write(*,5) '(8)           mugaz',tab_cntrl(tab0+8),mugaz
      write(*,5) '(9)             rcp',tab_cntrl(tab0+9),rcp
      write(*,6) '(11)        dtphys?',tab_cntrl(tab0+11),dtphys

      write(*,5) '(14)       year_day',tab_cntrl(tab0+14),year_day
      write(*,5) '(15)       periheli',tab_cntrl(tab0+15),periheli
      write(*,5) '(16)        aphelie',tab_cntrl(tab0+16),aphelie
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit

      write(*,6) '(19)     z0_default',tab_cntrl(tab0+19),z0_default
      write(*,6) '(21)      emin_turb',tab_cntrl(tab0+21),emin_turb
      write(*,5) '(20)        lmixmin',tab_cntrl(tab0+20),lmixmin

      write(*,5) '(26)        emissiv',tab_cntrl(tab0+26),emissiv
      write(*,5) '(24)     emisice(1)',tab_cntrl(tab0+24),emisice(1)
      write(*,5) '(25)     emisice(2)',tab_cntrl(tab0+25),emisice(2)
      write(*,5) '(22)    albedice(1)',tab_cntrl(tab0+22),albedice(1)
      write(*,5) '(23)    albedice(2)',tab_cntrl(tab0+23),albedice(2)
      write(*,6) '(31)   iceradius(1)',tab_cntrl(tab0+31),iceradius(1)
      write(*,6) '(32)   iceradius(2)',tab_cntrl(tab0+32),iceradius(2)
      write(*,5) '(33)   dtemisice(1)',tab_cntrl(tab0+33),dtemisice(1)
      write(*,5) '(34)   dtemisice(2)',tab_cntrl(tab0+34),dtemisice(2)

      write(*,5) '(27)         tauvis',tab_cntrl(tab0+27),tauvis

      write(*,5) '(35)        volcapa',tab_cntrl(tab0+35),volcapa

      write(*,*)
      write(*,*) 'Lmodif in tabfi!!!!!!!',Lmodif

c-----------------------------------------------------------------------
c	 Modifications...
! NB: Modifying controls should only be done by newstart, and in seq mode
      if ((Lmodif.eq.1).and.is_parallel) then
        write(*,*) "tabfi: Error modifying tab_control should",
     &             " only happen in serial mode (eg: by newstart)"
        stop
      endif
c-----------------------------------------------------------------------

      IF(Lmodif.eq.1) then

      write(*,*)
      write(*,*) 'Change values in tab_cntrl ? :'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*) '(Current values given above)'
      write(*,*)
      write(*,*) '(3)          day_ini : Initial day (=0 at Ls=0)'
      write(*,*) '(19)              z0 : default surface roughness (m)'
      write(*,*) '(21)       emin_turb :  minimal energy (PBL)'
      write(*,*) '(20)         lmixmin : mixing length (PBL)'
      write(*,*) '(26)         emissiv : ground emissivity'
      write(*,*) '(24 et 25)   emisice : CO2 ice max emissivity '
      write(*,*) '(22 et 23)  albedice : CO2 ice cap albedos'
      write(*,*) '(31 et 32) iceradius : mean scat radius of CO2 snow'
      write(*,*) '(33 et 34) dtemisice : time scale for snow',
     &           ' metamorphism'
      write(*,*) '(27)        tauvis : mean dust vis. reference ',
     &           'opacity'
      write(*,*) '(35)         volcapa : soil volumetric heat capacity'
      write(*,*) '(18)        obliquit : planet obliquity (deg)'
      write(*,*) '(17)      peri_day : perihelion date (sol since Ls=0)'
      write(*,*) '(  )      peri_ls : perihelion date (Ls since Ls=0)'
      write(*,*) '(15)      periheli : min. sun-mars dist (Mkm)'
      write(*,*) '(16)      aphelie  : max. sun-mars dist (Mkm)'
      write(*,*)
 
 
      do ! neverending loop
        write(*,*)
        write(*,*)
        write(*,*) 'Changes to perform ?'
        write(*,*) '   (enter keyword or return )'
        write(*,*)
        read(*,fmt='(a20)') modif
        if (modif(1:1) .eq. ' ') goto 999
 
        write(*,*)
        write(*,*) trim(modif) , ' : '

        if (trim(modif) .eq. 'day_ini') then
          write(*,*) 'current value:',day_ini
          write(*,*) 'enter new value:'
 101      read(*,*,iostat=ierr) day_ini
          if(ierr.ne.0) goto 101
          write(*,*) ' '
          write(*,*) 'day_ini (new value):',day_ini

        else if (trim(modif) .eq. 'z0') then
          write(*,*) 'current value (m):',z0_default
          write(*,*) 'enter new value (m):'
 102      read(*,*,iostat=ierr) z0_default
          if(ierr.ne.0) goto 102
          write(*,*) ' '
          write(*,*) ' z0 (new value):',z0_default

        else if (trim(modif) .eq. 'emin_turb') then
          write(*,*) 'current value:',emin_turb
          write(*,*) 'enter new value:'
 103      read(*,*,iostat=ierr) emin_turb
          if(ierr.ne.0) goto 103
          write(*,*) ' '
          write(*,*) ' emin_turb (new value):',emin_turb

        else if (trim(modif) .eq. 'lmixmin') then
          write(*,*) 'current value:',lmixmin
          write(*,*) 'enter new value:'
 104      read(*,*,iostat=ierr) lmixmin
          if(ierr.ne.0) goto 104
          write(*,*) ' '
          write(*,*) ' lmixmin (new value):',lmixmin

        else if (trim(modif) .eq. 'emissiv') then
          write(*,*) 'current value:',emissiv
          write(*,*) 'enter new value:'
 105      read(*,*,iostat=ierr) emissiv
          if(ierr.ne.0) goto 105
          write(*,*) ' '
          write(*,*) ' emissiv (new value):',emissiv

        else if (trim(modif) .eq. 'emisice') then
          write(*,*) 'current value emisice(1) North:',emisice(1)
          write(*,*) 'enter new value:'
 106      read(*,*,iostat=ierr) emisice(1)
          if(ierr.ne.0) goto 106
          write(*,*) 
          write(*,*) ' emisice(1) (new value):',emisice(1)
          write(*,*)

          write(*,*) 'current value emisice(2) South:',emisice(2)
          write(*,*) 'enter new value:'
 107      read(*,*,iostat=ierr) emisice(2)
          if(ierr.ne.0) goto 107
          write(*,*) 
          write(*,*) ' emisice(2) (new value):',emisice(2)

        else if (trim(modif) .eq. 'albedice') then
          write(*,*) 'current value albedice(1) North:',albedice(1)
          write(*,*) 'enter new value:'
 108      read(*,*,iostat=ierr) albedice(1)
          if(ierr.ne.0) goto 108
          write(*,*) 
          write(*,*) ' albedice(1) (new value):',albedice(1)
          write(*,*)

          write(*,*) 'current value albedice(2) South:',albedice(2)
          write(*,*) 'enter new value:'
 109      read(*,*,iostat=ierr) albedice(2)
          if(ierr.ne.0) goto 109
          write(*,*) 
          write(*,*) ' albedice(2) (new value):',albedice(2)

        else if (trim(modif) .eq. 'iceradius') then
          write(*,*) 'current value iceradius(1) North:',iceradius(1)
          write(*,*) 'enter new value:'
 110      read(*,*,iostat=ierr) iceradius(1)
          if(ierr.ne.0) goto 110
          write(*,*) 
          write(*,*) ' iceradius(1) (new value):',iceradius(1)
          write(*,*)

          write(*,*) 'current value iceradius(2) South:',iceradius(2)
          write(*,*) 'enter new value:'
 111      read(*,*,iostat=ierr) iceradius(2)
          if(ierr.ne.0) goto 111
          write(*,*) 
          write(*,*) ' iceradius(2) (new value):',iceradius(2)

        else if (trim(modif) .eq. 'dtemisice') then
          write(*,*) 'current value dtemisice(1) North:',dtemisice(1)
          write(*,*) 'enter new value:'
 112      read(*,*,iostat=ierr) dtemisice(1)
          if(ierr.ne.0) goto 112
          write(*,*) 
          write(*,*) ' dtemisice(1) (new value):',dtemisice(1)
          write(*,*)

          write(*,*) 'current value dtemisice(2) South:',dtemisice(2)
          write(*,*) 'enter new value:'
 113      read(*,*,iostat=ierr) dtemisice(2)
          if(ierr.ne.0) goto 113
          write(*,*) 
          write(*,*) ' dtemisice(2) (new value):',dtemisice(2)

        else if (trim(modif) .eq. 'tauvis') then
          write(*,*) 'current value:',tauvis
          write(*,*) 'enter new value:'
 114      read(*,*,iostat=ierr) tauvis
          if(ierr.ne.0) goto 114
          write(*,*) 
          write(*,*) ' tauvis (new value):',tauvis

        else if (trim(modif) .eq. 'obliquit') then
          write(*,*) 'current value:',obliquit
          write(*,*) 'obliquit should be 25.19 on current Mars'
          write(*,*) 'enter new value:'
 115      read(*,*,iostat=ierr) obliquit
          if(ierr.ne.0) goto 115
          write(*,*) 
          write(*,*) ' obliquit (new value):',obliquit

        else if (trim(modif) .eq. 'peri_day') then
          write(*,*) 'current value:',peri_day
          write(*,*) 'peri_day should be 485 sols on current Mars'
          write(*,*) 'enter new value:'
 116      read(*,*,iostat=ierr) peri_day
          if(ierr.ne.0) goto 116
          write(*,*) 
          write(*,*) ' peri_day (new value):',peri_day
          
        else if (trim(modif) .eq. 'peri_ls') then
          write(*,*) 'peri_ls value is not stored in start files,'
          write(*,*) 'but it should be 251 degrees on current Mars'
          write(*,*) '(peri_day should be 485 sols on current Mars)'
          write(*,*) 'enter new value:'
 1160     read(*,*,iostat=ierr) peri_ls
          if(ierr.ne.0) goto 1160
          write(*,*) 
          write(*,*) 'peri_ls asked:',peri_ls
          write(*,*) 'for aphelion =',aphelie
          write(*,*) 'perihelion =',periheli
          write(*,*) 'and',year_day,'sols/year'
          call lsp2solp(peri_ls,peri_day,aphelie,periheli,year_day)
          write(*,*) 'peri_day (new value):',peri_day


        else if (trim(modif) .eq. 'periheli') then
          write(*,*) 'current value:',periheli
          write(*,*) 'perihelion should be 206.66 on current Mars'
          write(*,*) 'enter new value:'
 117      read(*,*,iostat=ierr) periheli
          if(ierr.ne.0) goto 117
          write(*,*) 
          write(*,*) ' periheli (new value):',periheli
 
        else if (trim(modif) .eq. 'aphelie') then
          write(*,*) 'current value:',aphelie
          write(*,*) 'aphelion should be 249.22 on current Mars'
          write(*,*) 'enter new value:'
 118      read(*,*,iostat=ierr) aphelie
          if(ierr.ne.0) goto 118
          write(*,*) 
          write(*,*) ' aphelie (new value):',aphelie
 
        else if (trim(modif) .eq. 'volcapa') then
          write(*,*) 'current value:',volcapa
          write(*,*) 'enter new value:'
 119      read(*,*,iostat=ierr) volcapa
          if(ierr.ne.0) goto 119
          write(*,*) 
          write(*,*) ' volcapa (new value):',volcapa
 
        endif
      enddo ! of do ! neverending loop

 999  continue

c-----------------------------------------------------------------------
c	Write values of physical constants after modifications
c-----------------------------------------------------------------------
 
      write(*,*) '*****************************************************'
      write(*,*) 'Reading tab_cntrl when calling tabfi AFTER changes'
      write(*,*) '*****************************************************'
      write(*,5) '(1)        = ngrid?',tab_cntrl(tab0+1)
      write(*,5) '(2)            lmax',tab_cntrl(tab0+2),real(lmax)
      write(*,5) '(3)         day_ini',tab_cntrl(tab0+3),real(day_ini)
      write(*,5) '(5)             rad',tab_cntrl(tab0+5),rad
      write(*,5) '(10)         daysec',tab_cntrl(tab0+10),daysec
      write(*,6) '(6)            omeg',tab_cntrl(tab0+6),omeg
      write(*,5) '(7)               g',tab_cntrl(tab0+7),g
      write(*,5) '(8)           mugaz',tab_cntrl(tab0+8),mugaz
      write(*,5) '(9)             rcp',tab_cntrl(tab0+9),rcp
      write(*,6) '(11)        dtphys?',tab_cntrl(tab0+11),dtphys
 
      write(*,5) '(14)       year_day',tab_cntrl(tab0+14),year_day
      write(*,5) '(15)       periheli',tab_cntrl(tab0+15),periheli
      write(*,5) '(16)        aphelie',tab_cntrl(tab0+16),aphelie
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit
 
      write(*,6) '(19)     z0_default',tab_cntrl(tab0+19),z0_default
      write(*,6) '(21)      emin_turb',tab_cntrl(tab0+21),emin_turb
      write(*,5) '(20)        lmixmin',tab_cntrl(tab0+20),lmixmin
 
      write(*,5) '(26)        emissiv',tab_cntrl(tab0+26),emissiv
      write(*,5) '(24)     emisice(1)',tab_cntrl(tab0+24),emisice(1)
      write(*,5) '(25)     emisice(2)',tab_cntrl(tab0+25),emisice(2)
      write(*,5) '(22)    albedice(1)',tab_cntrl(tab0+22),albedice(1)
      write(*,5) '(23)    albedice(2)',tab_cntrl(tab0+23),albedice(2)
      write(*,6) '(31)   iceradius(1)',tab_cntrl(tab0+31),iceradius(1)
      write(*,6) '(32)   iceradius(2)',tab_cntrl(tab0+32),iceradius(2)
      write(*,5) '(33)   dtemisice(1)',tab_cntrl(tab0+33),dtemisice(1)
      write(*,5) '(34)   dtemisice(2)',tab_cntrl(tab0+34),dtemisice(2)
 
      write(*,5) '(27)         tauvis',tab_cntrl(tab0+27),tauvis

      write(*,5) '(35)        volcapa',tab_cntrl(tab0+35),volcapa

      write(*,*)  
      write(*,*) 

      ENDIF ! of if (Lmodif == 1)

c-----------------------------------------------------------------------
c	Case when using a start file from before March 1996 (without iceradius...
c-----------------------------------------------------------------------
      if (iceradius(1).eq.0) then 
         iceradius(1) = 100.e-6
         iceradius(2) = 100.e-6
         dtemisice(1) = 0.4
         dtemisice(2) = 0.4
         write (*,*) ' tabfi: WARNING : old initialisation file'
         write (*,*) 'iceradius set to',iceradius(1),iceradius(2)  
         write (*,*) 'dtemisice set to',dtemisice(1),dtemisice(2)  
       end if

c-----------------------------------------------------------------------
      end



      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gives sol at perihelion for ls at perihelion (for precession cycles)
      subroutine lsp2solp(lsp,solp,aphelie,periheli,year_day)

      implicit none
!  Arguments:
      real lsp     ! Input: ls at perihelion
      real solp    ! Output: sol at perihelion
      real aphelie,periheli,year_day ! Input: parameters
 
!  Local:
      double precision zx0 ! eccentric anomaly at Ls=0
      double precision e_elips
      double precision pi,degrad 
      
      parameter (pi=3.14159265358979d0)
      parameter (degrad=57.2957795130823d0)

      e_elips=(aphelie-periheli)/(aphelie+periheli)      
      zx0 = -2.0*datan(dtan(0.5*lsp/degrad)
     .          *dsqrt((1.-e_elips)/(1.+e_elips)))
      if (zx0 .le. 0.) zx0 = zx0 + 2.*pi
      
      solp  = year_day*(1.-(zx0-e_elips*dsin(zx0))/(2.*pi))


      end subroutine lsp2solp



