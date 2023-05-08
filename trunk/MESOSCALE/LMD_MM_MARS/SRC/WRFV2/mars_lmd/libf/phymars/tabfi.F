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
c      - lmax:    tab_cntrl(tab0+2) (pour test avec nlayermx)
c
c      - p_rad
c      - p_omeg   !
c      - p_g      ! Constantes physiques ayant des 
c      - p_mugaz  ! homonymes dynamiques
c      - p_daysec !
c
c=======================================================================
 
      implicit none
 
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "comgeomfi.h"
#include "planete.h"
#include "surfdat.h"
#include "netcdf.inc"
#include "dimradmars.h"
#include "yomaer.h"

c-----------------------------------------------------------------------
c   Declarations
c-----------------------------------------------------------------------

c Arguments
c ---------
      INTEGER nid,nvarid,tab0
      INTEGER*4 day_ini
      INTEGER Lmodif
      INTEGER lmax
      REAL p_rad,p_omeg,p_g,p_mugaz,p_daysec,time

c Variables
c ---------
      INTEGER length
      parameter (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER  ierr
      INTEGER size
      CHARACTER modif*20

c Fonctions DRS et autres
c -----------------------
      INTEGER setname, cluvdb, getdat

      INTEGER lnblnk

c-----------------------------------------------------------------------
c  Initialisation des constantes physiques par defaut (cas nid = 0)
c-----------------------------------------------------------------------
      IF (nid.eq.0) then
 
c Pression de reference sur la planete
c-------------------------------------
c     pressrf = 670.            ! Pression de reference (Pa) ~650

c Info sur la Planete Mars pour la dynamique et la physique
c----------------------------------------------------------
      rad=3397200.          ! rayon de mars (m)  ~3397200 m
      daysec=88775.         ! duree du sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec)       ! vitesse de rotation (rad.s-1)
      g=3.72                ! gravite (m.s-2) ~3.72
      mugaz=43.49           ! Masse molaire de l''atm (g.mol-1) ~43.49
      rcp=.256793         ! = r/cp  ~0.256793

c Info sur la Planete Mars pour la physique uniquement
c-----------------------------------------------------
      year_day = 669.       !modif FH  duree de l''annee (sols) ~668.6
      periheli = 206.66         ! dist.min. soleil-mars (Mkm) ~206.66
      aphelie = 249.22          ! dist.max. soleil-mars (Mkm) ~249.22
      peri_day =  485.    ! date du perihelie (sols depuis printemps)
      obliquit = 25.19         ! Obliquite de la planete (deg) ~25.19

c Couche limite et Turbulence
c----------------------------
      z0 =  1.e-2               ! surface roughness (m) ~0.01
      emin_turb = 1.e-6         ! energie minimale ~1.e-8
      lmixmin = 30              ! longueur de melange ~100

c propriete optiques des calottes et emissivite du sol
c-----------------------------------------------------
      emissiv=.95               ! Emissivite du sol martien ~.95
      emisice(1)=0.95           ! Emissivite calotte nord
      emisice(2)=0.95           ! Emissivite calotte sud
      albedice(1)=0.65          ! Albedo calotte nord
      albedice(2)=0.65          ! Albedo calotte sud
      iceradius(1) = 100.e-6    ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6    ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 0.4   ! time scale for snow metamorphism (north)
      dtemisice(2) = 0.4   ! time scale for snow metamorphism (south)

c Proprietes des poussiere aerosol
c---------------------------------
      tauvis= 0.2          ! profondeur optique visible moyenne

c  Ancien code radiatif (non utilise avec le code d'apres 03/96)
c---------------------------------------------------------------
c     tauir= 0.  ! .2  ratio (mean IR opt.depth)/Visible
c     scatalb=0. ! .86 scaterring albedo visible (~.86)
c     asfact=0.  ! .79 assymetrie factor visible   (~.79)
c     day0 = 0   ! = 0 en general !!!
 
c-----------------------------------------------------------------------
c  Initialisation des constantes physiques par lecture du tableau 
c		des parametres du run (tab_cntrl) 	(cas nid != 0)
c-----------------------------------------------------------------------
      ELSE

c Lecture du tableau
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Tabfi: Le champ <controle> est absent"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Tabfi: Lecture echoue pour <controle>"
         CALL abort
      ENDIF

       print*,'tab_cntrl',tab_cntrl
c
c  Initialisation des constantes physiques
c infor sur la grille physique
      if(ngridmx.ne.tab_cntrl(tab0+1)) then
         print*,'WARNING !!! tab_cntrl(tab0+1).ne.ngridmx'
         print*,tab_cntrl(tab0+1),ngridmx
      endif
      lmax = nint(tab_cntrl(tab0+2))
      day_ini = tab_cntrl(tab0+3)
      time = tab_cntrl(tab0+4)
      write (*,*) 'IN tabfi day_ini=',day_ini
c Info sur la Planete Mars pour la dynamique et la physique
      rad = tab_cntrl(tab0+5)
      omeg = tab_cntrl(tab0+6)
      g = tab_cntrl(tab0+7)
      mugaz = tab_cntrl(tab0+8)
      rcp = tab_cntrl(tab0+9)
      daysec = tab_cntrl(tab0+10)
      dtphys = tab_cntrl(tab0+11)
c Info sur la Planete Mars pour la physique uniquement
      year_day = tab_cntrl(tab0+14)
      periheli = tab_cntrl(tab0+15)
      aphelie = tab_cntrl(tab0+16)
      peri_day = tab_cntrl(tab0+17)
      obliquit = tab_cntrl(tab0+18)
c couche limite et turbeulence
      z0 = tab_cntrl(tab0+19)
      lmixmin = tab_cntrl(tab0+20)
      emin_turb = tab_cntrl(tab0+21)
c proprietes optiques des calottes et emissivite du sol
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
c-----------------------------------------------------------------------
c	Sauvegarde des constantes pour passages en arguments
c-----------------------------------------------------------------------
      p_omeg = omeg
      p_g = g
      p_mugaz = mugaz
      p_daysec = daysec
      p_rad=rad

      ENDIF    ! (du if nid = 0) 

c-----------------------------------------------------------------------
c	Impression des constantes physiques avant modif
c-----------------------------------------------------------------------
 
   6  FORMAT(a20,e15.6,e15.6)
   5  FORMAT(a20,f12.2,f12.2)
 
      write(*,*) '*****************************************************'
      write(*,*) 'lecture du tab_cntrl a l appel de tabfi avant modif'
      write(*,*) '*****************************************************'
      write(*,5) '(1)      = ngridmx?',tab_cntrl(tab0+1),float(ngridmx)
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
      write(*,5) '(15)       periheli',tab_cntrl(tab0+15),periheli
      write(*,5) '(16)        aphelie',tab_cntrl(tab0+16),aphelie
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit

      write(*,6) '(19)             z0',tab_cntrl(tab0+19),z0
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

      write(*,*)
      write(*,*) 'Lmodif dans tabfi!!!!!!!',Lmodif

c-----------------------------------------------------------------------
c	 Modifications...
c-----------------------------------------------------------------------

      IF(Lmodif.eq.1) then

      write(*,*)
      write(*,*) 'Modifications des valeurs du tab_cntrl:'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*) '(valeurs actuelles indiquees ci-dessus)'
      write(*,*)
      write(*,*) '(3)          day_ini : Jour initial (=0 a Ls=0)'
      write(*,*) '(19)              z0 :  surface roughness (m)'
      write(*,*) '(21)       emin_turb :  energie minimale'
      write(*,*) '(20)         lmixmin : longueur de melange'
      write(*,*) '(26)         emissiv : Emissivite du sol martien'
      write(*,*) '(24 et 25)   emisice : Emissivite des calottes '
      write(*,*) '(22 et 23)  albedice : Albedo des calotte'
      write(*,*) '(31 et 32) iceradius : mean scat radius of CO2 snow'
      write(*,*) '(33 et 34) dtemisice : time scale for snow',
     . ' metamorphism'
      write(*,*) '(27)          tauvis : profondeur optique visible',
     . ' moyenne'
      write(*,*) '(18)        obliquit : planet obliquity (deg)'
      write(*,*) '(17)      peri_day : perihelion date (sol since Ls=0)'
      write(*,*) '(15)      periheli : min. sun-mars dist (Mkm)'
      write(*,*) '(16)      aphelie  : max. sun-mars dist (Mkm)'
      write(*,*)
 
 
      do while(modif(1:1).ne.'hello')
        write(*,*)
        write(*,*)
        write(*,*) 'modification a faire ?'
        write(*,*) '   (entrer le mot cle ou return pour finir)'
        write(*,*)
        read(*,fmt='(a20)') modif
        if (modif(1:1) .eq. ' ') goto 999
 
        write(*,*)
        write(*,*) modif(1:lnblnk(modif)) , ' : '

        if (modif(1:lnblnk(modif)) .eq. 'day_ini') then
          write(*,*) 'valeur actuelle:',day_ini
          write(*,*) 'entrer la nouvelle valeur:'
 101      read(*,*,iostat=ierr) day_ini
          if(ierr.ne.0) goto 101
          write(*,*) ' '
          write(*,*) 'day_ini (nouvelle valeur):',day_ini

        else if (modif(1:lnblnk(modif)) .eq. 'z0') then
          write(*,*) 'valeur actuelle:',z0
          write(*,*) 'entrer la nouvelle valeur:'
 102      read(*,*,iostat=ierr) z0
          if(ierr.ne.0) goto 102
          write(*,*) ' '
          write(*,*) ' z0 (nouvelle valeur):',z0

        else if (modif(1:lnblnk(modif)) .eq. 'emin_turb') then
          write(*,*) 'valeur actuelle:',emin_turb
          write(*,*) 'entrer la nouvelle valeur:'
 103      read(*,*,iostat=ierr) emin_turb
          if(ierr.ne.0) goto 103
          write(*,*) ' '
          write(*,*) ' emin_turb (nouvelle valeur):',emin_turb

        else if (modif(1:lnblnk(modif)) .eq. 'lmixmin') then
          write(*,*) 'valeur actuelle:',lmixmin
          write(*,*) 'entrer la nouvelle valeur:'
 104      read(*,*,iostat=ierr) lmixmin
          if(ierr.ne.0) goto 104
          write(*,*) ' '
          write(*,*) ' lmixmin (nouvelle valeur):',lmixmin

        else if (modif(1:lnblnk(modif)) .eq. 'emissiv') then
          write(*,*) 'valeur actuelle:',emissiv
          write(*,*) 'entrer la nouvelle valeur:'
 105      read(*,*,iostat=ierr) emissiv
          if(ierr.ne.0) goto 105
          write(*,*) ' '
          write(*,*) ' emissiv (nouvelle valeur):',emissiv

        else if (modif(1:lnblnk(modif)) .eq. 'emisice') then
          write(*,*) 'valeur actuelle emisice(1) Nord:',emisice(1)
          write(*,*) 'entrer la nouvelle valeur:'
 106      read(*,*,iostat=ierr) emisice(1)
          if(ierr.ne.0) goto 106
          write(*,*) 
          write(*,*) ' emisice(1) (nouvelle valeur):',emisice(1)
          write(*,*)

          write(*,*) 'valeur actuelle emisice(2) Sud:',emisice(2)
          write(*,*) 'entrer la nouvelle valeur:'
 107      read(*,*,iostat=ierr) emisice(2)
          if(ierr.ne.0) goto 107
          write(*,*) 
          write(*,*) ' emisice(2) (nouvelle valeur):',emisice(2)

        else if (modif(1:lnblnk(modif)) .eq. 'albedice') then
          write(*,*) 'valeur actuelle albedice(1) Nord:',albedice(1)
          write(*,*) 'entrer la nouvelle valeur:'
 108      read(*,*,iostat=ierr) albedice(1)
          if(ierr.ne.0) goto 108
          write(*,*) 
          write(*,*) ' albedice(1) (nouvelle valeur):',albedice(1)
          write(*,*)

          write(*,*) 'valeur actuelle albedice(2) Sud:',albedice(2)
          write(*,*) 'entrer la nouvelle valeur:'
 109      read(*,*,iostat=ierr) albedice(2)
          if(ierr.ne.0) goto 109
          write(*,*) 
          write(*,*) ' albedice(2) (nouvelle valeur):',albedice(2)

        else if (modif(1:lnblnk(modif)) .eq. 'iceradius') then
          write(*,*) 'valeur actuelle iceradius(1) Nord:',iceradius(1)
          write(*,*) 'entrer la nouvelle valeur:'
 110      read(*,*,iostat=ierr) iceradius(1)
          if(ierr.ne.0) goto 110
          write(*,*) 
          write(*,*) ' iceradius(1) (nouvelle valeur):',iceradius(1)
          write(*,*)

          write(*,*) 'valeur actuelle iceradius(2) Sud:',iceradius(2)
          write(*,*) 'entrer la nouvelle valeur:'
 111      read(*,*,iostat=ierr) iceradius(2)
          if(ierr.ne.0) goto 111
          write(*,*) 
          write(*,*) ' iceradius(2) (nouvelle valeur):',iceradius(2)

        else if (modif(1:lnblnk(modif)) .eq. 'dtemisice') then
          write(*,*) 'valeur actuelle dtemisice(1) Nord:',dtemisice(1)
          write(*,*) 'entrer la nouvelle valeur:'
 112      read(*,*,iostat=ierr) dtemisice(1)
          if(ierr.ne.0) goto 112
          write(*,*) 
          write(*,*) ' dtemisice(1) (nouvelle valeur):',dtemisice(1)
          write(*,*)

          write(*,*) 'valeur actuelle dtemisice(2) Sud:',dtemisice(2)
          write(*,*) 'entrer la nouvelle valeur:'
 113      read(*,*,iostat=ierr) dtemisice(2)
          if(ierr.ne.0) goto 113
          write(*,*) 
          write(*,*) ' dtemisice(2) (nouvelle valeur):',dtemisice(2)

        else if (modif(1:lnblnk(modif)) .eq. 'tauvis') then
          write(*,*) 'valeur actuelle:',tauvis
          write(*,*) 'entrer la nouvelle valeur:'
 114      read(*,*,iostat=ierr) tauvis
          if(ierr.ne.0) goto 114
          write(*,*) 
          write(*,*) ' tauvis (nouvelle valeur):',tauvis
        else if (modif(1:lnblnk(modif)) .eq. 'obliquit') then
          write(*,*) 'valeur actuelle:',obliquit
          write(*,*) 'obliquit should be 25.19 on current Mars'
          write(*,*) 'entrer la nouvelle valeur:'
 115      read(*,*,iostat=ierr) obliquit
          if(ierr.ne.0) goto 115
          write(*,*) 
          write(*,*) ' obliquit (nouvelle valeur):',obliquit

        else if (modif(1:lnblnk(modif)) .eq. 'peri_day') then
          write(*,*) 'valeur actuelle:',peri_day
          write(*,*) 'peri_day should be 485 on current Mars'
          write(*,*) 'entrer la nouvelle valeur:'
 116      read(*,*,iostat=ierr) peri_day
          if(ierr.ne.0) goto 116
          write(*,*) 
          write(*,*) ' peri_day (nouvelle valeur):',peri_day

        else if (modif(1:lnblnk(modif)) .eq. 'periheli') then
          write(*,*) 'valeur actuelle:',periheli
          write(*,*) 'periheli should be 206.66 on current Mars'
          write(*,*) 'entrer la nouvelle valeur:'
 117      read(*,*,iostat=ierr) periheli
          if(ierr.ne.0) goto 117
          write(*,*) 
          write(*,*) ' periheli (nouvelle valeur):',periheli
 
        else if (modif(1:lnblnk(modif)) .eq. 'aphelie') then
          write(*,*) 'valeur actuelle:',aphelie
          write(*,*) 'aphelie should be 249.22 on current Mars'
          write(*,*) 'entrer la nouvelle valeur:'
 118      read(*,*,iostat=ierr) aphelie
          if(ierr.ne.0) goto 118
          write(*,*) 
          write(*,*) ' aphelie (nouvelle valeur):',aphelie
 
        endif
      enddo

 999  continue

c-----------------------------------------------------------------------
c	Impression des constantes physiques apres modif si modif
c-----------------------------------------------------------------------
 
      write(*,*) '*****************************************************'
      write(*,*) 'lecture du tab_cntrl a l appel de tabfi apres modif'
      write(*,*) '*****************************************************'
      write(*,5) '(1)      = ngridmx?',tab_cntrl(tab0+1),float(ngridmx)
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
      write(*,5) '(15)       periheli',tab_cntrl(tab0+15),periheli
      write(*,5) '(16)        aphelie',tab_cntrl(tab0+16),aphelie
      write(*,5) '(17)       peri_day',tab_cntrl(tab0+17),peri_day
      write(*,5) '(18)       obliquit',tab_cntrl(tab0+18),obliquit
 
      write(*,6) '(19)             z0',tab_cntrl(tab0+19),z0
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

      write(*,*)  
      write(*,*) 

      ENDIF 			!	(du if Lmodif = 1)

c-----------------------------------------------------------------------
c	Case when using a start file from before March 1996 (without iceradius...
c-----------------------------------------------------------------------
      if (iceradius(1).eq.0) then 
         iceradius(1) = 100.e-6
         iceradius(2) = 100.e-6
         dtemisice(1) = 0.4
         dtemisice(2) = 0.4
         write (*,*) ' WARNING : old initialisation file'
         write (*,*) 'iceradius set to',iceradius(1),iceradius(2)  
         write (*,*) 'dtemisice set to',dtemisice(1),dtemisice(2)  
       end if

c-----------------------------------------------------------------------
      end
