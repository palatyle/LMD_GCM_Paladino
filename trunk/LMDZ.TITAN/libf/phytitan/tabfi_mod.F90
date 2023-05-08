MODULE tabfi_mod

IMPLICIT NONE

CONTAINS

!=======================================================================
      SUBROUTINE tabfi(ngrid,nid,Lmodif,tab0,day_ini,lmax,p_rad, &
                       p_omeg,p_g,p_cpp,p_mugaz,p_daysec,time)
!=======================================================================
!
!   C. Hourdin 15/11/96
!
!   Object:        Lecture du tab_cntrl physique dans un fichier 
!   ------            et initialisation des constantes physiques
!
!   Arguments:
!   ----------
!
!     Inputs:
!     ------
!
!      - nid:    unitne logique du fichier ou on va lire le tab_cntrl    
!                      (ouvert dans le programme appellant) 
!
!                 si nid=0:
!                       pas de lecture du tab_cntrl mais
!                       Valeurs par default des constantes physiques
!        
!      - tab0:    Offset de tab_cntrl a partir duquel sont ranges 
!                  les parametres physiques (50 pour start_archive)
!
!      - Lmodif:  si on souhaite modifier les constantes  Lmodif = 1 = TRUE
!
!
!     Outputs:
!     --------
!
!      - day_ini: tab_cntrl(tab0+3) (Dans les cas ou l'on souhaite
!                              comparer avec le day_ini dynamique)
!
!      - lmax:    tab_cntrl(tab0+2) (pour test avec nlayer)
!
!      - p_rad
!      - p_omeg   !
!      - p_g      ! Constantes physiques ayant des 
!      - p_mugaz  ! homonymes dynamiques
!      - p_daysec !
!
!=======================================================================
! to use  'getin_p'
      use ioipsl_getin_p_mod, only: getin_p

      use surfdat_h, only: emisice, iceradius, dtemisice, &
                           emissiv
      use comsoil_h, only: volcapa
      use iostart, only: get_var
      use mod_phys_lmdz_para, only: is_parallel
      use planete_mod, only: year_day, periastr, apoastr, peri_day, &
                             obliquit, z0, lmixmin, emin_turb
      use comcstfi_mod, only: rad, omeg, g, mugaz, rcp, cpp, r
      use time_phylmdz_mod, only: dtphys, daysec
      use callkeys_mod, only: check_cpp_match,force_cpp
      implicit none
 
      include "netcdf.inc"

!-----------------------------------------------------------------------
!   Declarations
!-----------------------------------------------------------------------

! Arguments
! ---------
      INTEGER,INTENT(IN) :: ngrid,nid,tab0
      INTEGER*4,INTENT(OUT) :: day_ini
      INTEGER,INTENT(IN) :: Lmodif
      INTEGER,INTENT(OUT) :: lmax
      REAL,INTENT(OUT) :: p_rad,p_omeg,p_g,p_cpp,p_mugaz,p_daysec,time

! Variables
! ---------
      INTEGER,PARAMETER :: length=100
      REAL tab_cntrl(length) ! array in which are stored the run's parameters
      INTEGER  ierr,nvarid
      INTEGER size
      CHARACTER modif*20
      LOGICAL :: found
      CHARACTER(len=5) :: modname="tabfi"
      
      write(*,*)"tabfi: nid=",nid," tab0=",tab0," Lmodif=",Lmodif

      IF (nid.eq.0) then
!-----------------------------------------------------------------------
!  Initialization of various physical constants to defaut values (nid = 0 case)
!-----------------------------------------------------------------------
        tab_cntrl(:)=0
        lmax=0 ! not used anyways
        !day_ini already set via inifis
        time=0 
! Informations about planet for dynamics and physics
        ! rad,cpp,g,r,rcp already initialized by inifis
        omeg=-999.
        call getin_p("omega",omeg)
        if (omeg.eq.-999.) then
          call abort_physic(modname,"Missing value for omega in def files!",1)
        endif
        mugaz=(8.3144621/r)*1.E3
        ! daysec already set by inifis
        ! dtphys alread set by inifis
! Informations about planet for the physics only
        year_day=-999. ! length of year, in standard days
        call getin_p("year_day",year_day)
        if (year_day.eq.-999.) then
          call abort_physic(modname, &
               "Missing value for year_day in def files!",1)
        endif
        periastr=-999.
        call getin_p("periastron",periastr)
        if (periastr.eq.-999.) then
          call abort_physic(modname, &
               "Missing value for periastron in def files!",1)
        endif
        apoastr=-999.
        call getin_p("apoastron",apoastr)
        if (apoastr.eq.-999.) then
          call abort_physic(modname, &
               "Missing value for apoastron in def files!",1)
        endif
        peri_day=-999.
        call getin_p("periastron_day",peri_day)
        if (peri_day.eq.-999.) then
          call abort_physic(modname, &
               "Missing value for periastron date in def files!",1)
        endif
        obliquit=-999.
        call getin_p("obliquity",obliquit)
        if (obliquit.eq.-999.) then
          call abort_physic(modname, &
               "Missing value for obliquity in def files!",1)
        endif
! boundary layer and turbulence
        z0=1.e-2 ! surface roughness length (m)
        lmixmin=30
        emin_turb=1.e-6
! optical properties of polar caps and ground emissivity
        emisice(:)=0
        emissiv=0
        iceradius(:)=1.e-6 ! mean scat radius of CO2 snow
        dtemisice(:)=0 !time scale for snow metamorphism
        volcapa=1000000 ! volumetric heat capacity of subsurface
        
      ELSE
!-----------------------------------------------------------------------
!  Initialization of physical constants by reading array tab_cntrl(:)
!		which contains these parameters	(nid != 0 case)
!-----------------------------------------------------------------------
! Read 'controle' array
!

       call get_var("controle",tab_cntrl,found)
       if (.not.found) then
         call abort_physic(modname,"Failed reading <controle> array",1)
       else
         write(*,*)'tabfi: tab_cntrl',tab_cntrl
       endif
!
!  Initialization of some physical constants
! informations on physics grid
!      if(ngrid.ne.tab_cntrl(tab0+1)) then
!         print*,'tabfi: WARNING !!! tab_cntrl(tab0+1).ne.ngrid'
!         print*,tab_cntrl(tab0+1),ngrid
!      endif
      lmax = nint(tab_cntrl(tab0+2))
      day_ini = tab_cntrl(tab0+3)
      time = tab_cntrl(tab0+4)
      write (*,*) 'IN tabfi day_ini=',day_ini
! Informations about planet for dynamics and physics
      rad = tab_cntrl(tab0+5)
      omeg = tab_cntrl(tab0+6)
      g = tab_cntrl(tab0+7)
      mugaz = tab_cntrl(tab0+8)
      rcp = tab_cntrl(tab0+9)
      cpp=(8.314511/(mugaz/1000.0))/rcp
      daysec = tab_cntrl(tab0+10)
      dtphys = tab_cntrl(tab0+11)
! Informations about planet for the physics only
      year_day = tab_cntrl(tab0+14)
      periastr = tab_cntrl(tab0+15)
      apoastr = tab_cntrl(tab0+16)
      peri_day = tab_cntrl(tab0+17)
      obliquit = tab_cntrl(tab0+18)
! boundary layer and turbulence
      z0 = tab_cntrl(tab0+19)
      lmixmin = tab_cntrl(tab0+20)
      emin_turb = tab_cntrl(tab0+21)
! optical properties of polar caps and ground emissivity
      emisice(1) = tab_cntrl(tab0+24)
      emisice(2) = tab_cntrl(tab0+25)
      emissiv    = tab_cntrl(tab0+26)
      iceradius(1)= tab_cntrl(tab0+31) ! mean scat radius of CO2 snow (north)
      iceradius(2)= tab_cntrl(tab0+32) ! mean scat radius of CO2 snow (south)
      dtemisice(1)= tab_cntrl(tab0+33) !time scale for snow metamorphism (north)
      dtemisice(2)= tab_cntrl(tab0+34) !time scale for snow metamorphism (south)
! soil properties
      volcapa = tab_cntrl(tab0+35) ! volumetric heat capacity
!-----------------------------------------------------------------------
!	Save some constants for later use (as routine arguments)
!-----------------------------------------------------------------------
      p_omeg = omeg
      p_g = g
      p_cpp = cpp
      p_mugaz = mugaz
      p_daysec = daysec
      p_rad=rad

      ENDIF    ! end of (nid = 0) 

!-----------------------------------------------------------------------
!	Write physical constants to output before modifying them
!-----------------------------------------------------------------------
 
   6  FORMAT(a20,e15.6,e15.6)
   5  FORMAT(a20,f12.2,f12.2)
 
      write(*,*) '*****************************************************'
      write(*,*) 'Reading tab_cntrl when calling tabfi before changes'
      write(*,*) '*****************************************************'
      write(*,5) '(1)        = ngrid?',tab_cntrl(tab0+1),float(ngrid)
      write(*,5) '(2)            lmax',tab_cntrl(tab0+2),float(lmax)
      write(*,5) '(3)         day_ini',tab_cntrl(tab0+3),float(day_ini)
      write(*,5) '(5)             rad',tab_cntrl(tab0+5),rad
      write(*,5) '(10)         daysec',tab_cntrl(tab0+10),daysec
      write(*,6) '(6)            omeg',tab_cntrl(tab0+6),omeg
      write(*,5) '(7)               g',tab_cntrl(tab0+7),g
      write(*,5) '(8)           mugaz',tab_cntrl(tab0+8),mugaz
      write(*,5) '(9)             rcp',tab_cntrl(tab0+9),rcp
      write(*,6) '(11)        dtphys?',tab_cntrl(tab0+11),dtphys

      write(*,5) '(14)       year_day',tab_cntrl(tab0+14),year_day
      write(*,5) '(15)       periastr',tab_cntrl(tab0+15),periastr
      write(*,5) '(16)        apoastr',tab_cntrl(tab0+16),apoastr
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit

      write(*,6) '(19)             z0',tab_cntrl(tab0+19),z0
      write(*,6) '(21)      emin_turb',tab_cntrl(tab0+21),emin_turb
      write(*,5) '(20)        lmixmin',tab_cntrl(tab0+20),lmixmin

      write(*,5) '(26)        emissiv',tab_cntrl(tab0+26),emissiv
      write(*,5) '(24)     emisice(1)',tab_cntrl(tab0+24),emisice(1)
      write(*,5) '(25)     emisice(2)',tab_cntrl(tab0+25),emisice(2)
      write(*,6) '(31)   iceradius(1)',tab_cntrl(tab0+31),iceradius(1)
      write(*,6) '(32)   iceradius(2)',tab_cntrl(tab0+32),iceradius(2)
      write(*,5) '(33)   dtemisice(1)',tab_cntrl(tab0+33),dtemisice(1)
      write(*,5) '(34)   dtemisice(2)',tab_cntrl(tab0+34),dtemisice(2)

      write(*,5) '(35)        volcapa',tab_cntrl(tab0+35),volcapa

      write(*,*)
      write(*,*) 'Lmodif in tabfi!!!!!!!',Lmodif

!-----------------------------------------------------------------------
!	 Modifications...
! NB: Modifying controls should only be done by newstart, and in seq mode
      if ((Lmodif.eq.1).and.is_parallel) then
        write(*,*) "tabfi: Error modifying tab_control should", &
                   " only happen in serial mode (eg: by newstart)"
        stop
      endif
!-----------------------------------------------------------------------

      IF(Lmodif.eq.1) then

      write(*,*)
      write(*,*) 'Change values in tab_cntrl ? :'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*) '(Current values given above)'
      write(*,*)
      write(*,*) '(3)          day_ini : Initial day (=0 at Ls=0)'
      write(*,*) '(19)              z0 :  surface roughness (m)'
      write(*,*) '(21)       emin_turb :  minimal energy (PBL)'
      write(*,*) '(20)         lmixmin : mixing length (PBL)'
      write(*,*) '(26)         emissiv : ground emissivity'
      write(*,*) '(24 et 25)   emisice : CO2 ice max emissivity '
      write(*,*) '(31 et 32) iceradius : mean scat radius of CO2 snow'
      write(*,*) '(33 et 34) dtemisice : time scale for snow metamorphism'
      write(*,*) '(35)      volcapa : soil volumetric heat capacity'
      write(*,*) '(18)     obliquit : planet obliquity (deg)'
      write(*,*) '(17)     peri_day : periastron date (sols since Ls=0)'
      write(*,*) '(15)     periastr : min. star-planet dist (UA)'
      write(*,*) '(16)     apoastr  : max. star-planet (UA)'
      write(*,*) '(14)     year_day : length of year (in sols)'
      write(*,*) '(5)      rad      : radius of the planet (m)'
      write(*,*) '(6)      omeg     : planet rotation rate (rad/s)'
      write(*,*) '(7)      g        : gravity (m/s2)'
      write(*,*) '(8)      mugaz    : molecular mass '
      write(*,*) '                       of the atmosphere (g/mol)'
      write(*,*) '(9)      rcp      : r/Cp'
      write(*,*) '(8)+(9)  calc_cpp_mugaz : r/Cp and mugaz '
      write(*,*) '                 computed from gases.def'
      write(*,*) '(10)     daysec   : length of a sol (s)'
      write(*,*)
 
 
      do while(modif(1:1).ne.'hello')
        write(*,*)
        write(*,*)
        write(*,*) 'Changes to perform ?'
        write(*,*) '   (enter keyword or return )'
        write(*,*)
        read(*,fmt='(a20)') modif
        if (modif(1:1) .eq. ' ') goto 999
 
        write(*,*)
        write(*,*) modif(1:len_trim(modif)) , ' : '

        if (modif(1:len_trim(modif)) .eq. 'day_ini') then
          write(*,*) 'current value:',day_ini
          write(*,*) 'enter new value:'
 101      read(*,*,iostat=ierr) day_ini
          if(ierr.ne.0) goto 101
          write(*,*) ' '
          write(*,*) 'day_ini (new value):',day_ini

        else if (modif(1:len_trim(modif)) .eq. 'z0') then
          write(*,*) 'current value:',z0
          write(*,*) 'enter new value:'
 102      read(*,*,iostat=ierr) z0
          if(ierr.ne.0) goto 102
          write(*,*) ' '
          write(*,*) ' z0 (new value):',z0

        else if (modif(1:len_trim(modif)) .eq. 'emin_turb') then
          write(*,*) 'current value:',emin_turb
          write(*,*) 'enter new value:'
 103      read(*,*,iostat=ierr) emin_turb
          if(ierr.ne.0) goto 103
          write(*,*) ' '
          write(*,*) ' emin_turb (new value):',emin_turb

        else if (modif(1:len_trim(modif)) .eq. 'lmixmin') then
          write(*,*) 'current value:',lmixmin
          write(*,*) 'enter new value:'
 104      read(*,*,iostat=ierr) lmixmin
          if(ierr.ne.0) goto 104
          write(*,*) ' '
          write(*,*) ' lmixmin (new value):',lmixmin

        else if (modif(1:len_trim(modif)) .eq. 'emissiv') then
          write(*,*) 'current value:',emissiv
          write(*,*) 'enter new value:'
 105      read(*,*,iostat=ierr) emissiv
          if(ierr.ne.0) goto 105
          write(*,*) ' '
          write(*,*) ' emissiv (new value):',emissiv

        else if (modif(1:len_trim(modif)) .eq. 'emisice') then
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

        else if (modif(1:len_trim(modif)) .eq. 'iceradius') then
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

        else if (modif(1:len_trim(modif)) .eq. 'dtemisice') then
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

        else if (modif(1:len_trim(modif)) .eq. 'obliquit') then
          write(*,*) 'current value:',obliquit
          write(*,*) 'obliquit should be 25.19 on current Mars'
          write(*,*) 'enter new value:'
 115      read(*,*,iostat=ierr) obliquit
          if(ierr.ne.0) goto 115
          write(*,*) 
          write(*,*) ' obliquit (new value):',obliquit

        else if (modif(1:len_trim(modif)) .eq. 'peri_day') then
          write(*,*) 'current value:',peri_day
          write(*,*) 'peri_day should be 485 on current Mars'
          write(*,*) 'enter new value:'
 116      read(*,*,iostat=ierr) peri_day
          if(ierr.ne.0) goto 116
          write(*,*) 
          write(*,*) ' peri_day (new value):',peri_day

        else if (modif(1:len_trim(modif)) .eq. 'periastr') then
          write(*,*) 'current value:',periastr
          write(*,*) 'periastr should be 206.66 on present-day Mars'
          write(*,*) 'enter new value:'
 117      read(*,*,iostat=ierr) periastr
          if(ierr.ne.0) goto 117
          write(*,*) 
          write(*,*) ' periastr (new value):',periastr
 
        else if (modif(1:len_trim(modif)) .eq. 'apoastr') then
          write(*,*) 'current value:',apoastr
          write(*,*) 'apoastr should be 249.22 on present-day Mars'
          write(*,*) 'enter new value:'
 118      read(*,*,iostat=ierr) apoastr
          if(ierr.ne.0) goto 118
          write(*,*) 
          write(*,*) ' apoastr (new value):',apoastr
 
        else if (modif(1:len_trim(modif)) .eq. 'volcapa') then
          write(*,*) 'current value:',volcapa
          write(*,*) 'enter new value:'
 119      read(*,*,iostat=ierr) volcapa
          if(ierr.ne.0) goto 119
          write(*,*) 
          write(*,*) ' volcapa (new value):',volcapa
        
        else if (modif(1:len_trim(modif)).eq.'rad') then
          write(*,*) 'current value:',rad
          write(*,*) 'enter new value:'
 120      read(*,*,iostat=ierr) rad
          if(ierr.ne.0) goto 120
          write(*,*) 
          write(*,*) ' rad (new value):',rad

        else if (modif(1:len_trim(modif)).eq.'omeg') then
          write(*,*) 'current value:',omeg
          write(*,*) 'enter new value:'
 121      read(*,*,iostat=ierr) omeg
          if(ierr.ne.0) goto 121
          write(*,*) 
          write(*,*) ' omeg (new value):',omeg
        
        else if (modif(1:len_trim(modif)).eq.'g') then
          write(*,*) 'current value:',g
          write(*,*) 'enter new value:'
 122      read(*,*,iostat=ierr) g
          if(ierr.ne.0) goto 122
          write(*,*) 
          write(*,*) ' g (new value):',g

        else if (modif(1:len_trim(modif)).eq.'mugaz') then
          write(*,*) 'current value:',mugaz
          write(*,*) 'enter new value:'
 123      read(*,*,iostat=ierr) mugaz
          if(ierr.ne.0) goto 123
          write(*,*) 
          write(*,*) ' mugaz (new value):',mugaz
          r=8.314511/(mugaz/1000.0)
          write(*,*) ' R (new value):',r

        else if (modif(1:len_trim(modif)).eq.'rcp') then
          write(*,*) 'current value:',rcp
          write(*,*) 'enter new value:'
 124      read(*,*,iostat=ierr) rcp
          if(ierr.ne.0) goto 124
          write(*,*) 
          write(*,*) ' rcp (new value):',rcp
          r=8.314511/(mugaz/1000.0)
          cpp=r/rcp
          write(*,*) ' cpp (new value):',cpp

        else if (modif(1:len_trim(modif)).eq.'calc_cpp_mugaz') then
          write(*,*) 'current value rcp, mugaz:',rcp,mugaz
          check_cpp_match=.false.
	  force_cpp=.false.
	  call su_gases
	  call calc_cpp_mugaz
          write(*,*) 
          write(*,*) ' cpp (new value):',cpp
          write(*,*) ' mugaz (new value):',mugaz
          r=8.314511/(mugaz/1000.0)
          rcp=r/cpp
          write(*,*) ' rcp (new value):',rcp
	  
        else if (modif(1:len_trim(modif)).eq.'daysec') then
          write(*,*) 'current value:',daysec
          write(*,*) 'enter new value:'
 125      read(*,*,iostat=ierr) daysec
          if(ierr.ne.0) goto 125
          write(*,*) 
          write(*,*) ' daysec (new value):',daysec

!         added by RW!
        else if (modif(1:len_trim(modif)).eq.'year_day') then
          write(*,*) 'current value:',year_day
          write(*,*) 'enter new value:' 
 126      read(*,*,iostat=ierr) year_day
          if(ierr.ne.0) goto 126
          write(*,*)
          write(*,*) ' year_day (new value):',year_day

        endif
      enddo ! of do while(modif(1:1).ne.'hello')

 999  continue

!----------------------------------------------------------------------
!	Write values of physical constants after modifications
!-----------------------------------------------------------------------
 
      write(*,*) '*****************************************************'
      write(*,*) 'Reading tab_cntrl when calling tabfi AFTER changes'
      write(*,*) '*****************************************************'
      write(*,5) '(1)        = ngrid?',tab_cntrl(tab0+1),float(ngrid)
      write(*,5) '(2)            lmax',tab_cntrl(tab0+2),float(lmax)
      write(*,5) '(3)         day_ini',tab_cntrl(tab0+3),float(day_ini)
      write(*,5) '(5)             rad',tab_cntrl(tab0+5),rad
      write(*,5) '(10)         daysec',tab_cntrl(tab0+10),daysec
      write(*,6) '(6)            omeg',tab_cntrl(tab0+6),omeg
      write(*,5) '(7)               g',tab_cntrl(tab0+7),g
      write(*,5) '(8)           mugaz',tab_cntrl(tab0+8),mugaz
      write(*,5) '(9)             rcp',tab_cntrl(tab0+9),rcp
      write(*,6) '(11)        dtphys?',tab_cntrl(tab0+11),dtphys
 
      write(*,5) '(14)       year_day',tab_cntrl(tab0+14),year_day
      write(*,5) '(15)       periastr',tab_cntrl(tab0+15),periastr
      write(*,5) '(16)        apoastr',tab_cntrl(tab0+16),apoastr
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit
 
      write(*,6) '(19)             z0',tab_cntrl(tab0+19),z0
      write(*,6) '(21)      emin_turb',tab_cntrl(tab0+21),emin_turb
      write(*,5) '(20)        lmixmin',tab_cntrl(tab0+20),lmixmin
 
      write(*,5) '(26)        emissiv',tab_cntrl(tab0+26),emissiv
      write(*,5) '(24)     emisice(1)',tab_cntrl(tab0+24),emisice(1)
      write(*,5) '(25)     emisice(2)',tab_cntrl(tab0+25),emisice(2)
      write(*,6) '(31)   iceradius(1)',tab_cntrl(tab0+31),iceradius(1)
      write(*,6) '(32)   iceradius(2)',tab_cntrl(tab0+32),iceradius(2)
      write(*,5) '(33)   dtemisice(1)',tab_cntrl(tab0+33),dtemisice(1)
      write(*,5) '(34)   dtemisice(2)',tab_cntrl(tab0+34),dtemisice(2)
 
      write(*,5) '(35)        volcapa',tab_cntrl(tab0+35),volcapa

      write(*,*)  
      write(*,*) 

      ENDIF ! of if (Lmodif == 1)

!-----------------------------------------------------------------------
!	Save some constants for later use (as routine arguments)
!-----------------------------------------------------------------------
      p_omeg = omeg
      p_g = g
      p_cpp = cpp
      p_mugaz = mugaz
      p_daysec = daysec
      p_rad=rad


      END SUBROUTINE tabfi

end module tabfi_mod
