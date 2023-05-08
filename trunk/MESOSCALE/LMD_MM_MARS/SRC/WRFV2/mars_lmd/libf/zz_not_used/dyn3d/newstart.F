C======================================================================
      PROGRAM newstart
c=======================================================================
c
c
c   Auteur:   Christophe Hourdin/Francois Forget/Yann Wanherdrick
c   ------
c             Derniere modif : 12/03
c
c
c   Objet:  Create or modify the initial state for the LMD Mars GCM
c   -----           (fichiers NetCDF start et startfi)
c
c
c=======================================================================

      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "surfdat.h"
#include "dimradmars.h"
#include "yomaer.h"
#include "planete.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom2.h"
#include "control.h"
#include "logic.h"
#include "description.h"
#include "ener.h"
#include "temps.h"
#include "lmdstd.h"
#include "comdissnew.h"
#include "clesph0.h"
#include "serre.h"
#include "netcdf.inc"

c=======================================================================
c   Declarations
c=======================================================================

c Variables dimension du fichier "start_archive"
c------------------------------------
      CHARACTER	relief*3

c et autres:
c----------
      INTEGER lnblnk
      EXTERNAL lnblnk

c Variables pour les lectures NetCDF des fichiers "start_archive" 
c--------------------------------------------------
      INTEGER nid_dyn, nid_fi,nid,nvarid
      INTEGER size
      INTEGER length
      parameter (length = 100)
      INTEGER tab0
      INTEGER   NB_ETATMAX
      parameter (NB_ETATMAX = 100)

      REAL  date
      REAL p_rad,p_omeg,p_g,p_mugaz,p_daysec

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL phis(iip1,jjp1)
      REAL q(iip1,jjp1,llm,nqmx)               ! champs advectes

c autre variables dynamique nouvelle grille
c------------------------------------------
      REAL pks(iip1,jjp1)
      REAL w(iip1,jjp1,llm+1)
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dh(ip1jmp1,llm),dp(ip1jmp1)
      REAL phi(iip1,jjp1,llm)

      integer klatdat,klongdat
      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------
      real zmeaS(iip1,jjp1),zstdS(iip1,jjp1)
      real zsigS(iip1,jjp1),zgamS(iip1,jjp1),ztheS(iip1,jjp1)

c variable physique
c------------------
      REAL tsurf(ngridmx),tsoil(ngridmx,nsoilmx),co2ice(ngridmx)
      REAL emis(ngridmx),qsurf(ngridmx,nqmx)
      REAL q2(ngridmx,nlayermx+1)
      REAL rnaturfi(ngridmx)
      real alb(iip1,jjp1),albfi(ngridmx)
      real ith(iip1,jjp1),ithfi(ngridmx)
      REAL latfi(ngridmx),lonfi(ngridmx),airefi(ngridmx)

      INTEGER i,j,l,idum,ig
      real mugaz ! masse molaire de l''atmosphere

      EXTERNAL iniconst,geopot,inigeom
      integer ierr, nbetat
      integer ismin
      external ismin

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      REAL p(iip1,jjp1)
      REAL t(iip1,jjp1,llm)
      real phisold_newgrid(iip1,jjp1)
      REAL :: teta(iip1, jjp1, llm)
      REAL :: pk(iip1,jjp1,llm)
      REAL :: pkf(iip1,jjp1,llm)
      REAL :: ps(iip1, jjp1)
      REAL :: masse(iip1,jjp1,llm)
      REAL :: xpn,xps,xppn(iim),xpps(iim)
      REAL :: p3d(iip1, jjp1, llm+1)
      REAL :: beta(iip1,jjp1,llm)
      REAL dteta(ip1jmp1,llm)

c Variable de l'ancienne grille 
c------------------------------
      real time
      real tab_cntrl(100)
      real tab_cntrl_bis(100)

c variables diverses
c-------------------
      real choix_1
      character*80      fichnom
      integer Lmodif,iq,thermo
      character modif*20
      real z_reel(iip1,jjp1)
      real tsud,albsud,alb_bb,ith_bb,Tiso
      real ptoto,pcap,patm,airetot,ptotn,patmn
      real ssum
      character*1 yes
      logical :: flagiso=.false.
      real val

      INTEGER :: itau

c sortie visu pour les champs dynamiques
c---------------------------------------
      INTEGER :: visuid
      real :: time_step,t_ops,t_wrt
      CHARACTER*80 :: visu_file

      cpp    = 744.499 ! au lieu de 1004.70885 (Terre)
      preff  = 610.    ! au lieu de 101325. (Terre)
      pa= 20           ! au lieu de 500 (Terre)

c=======================================================================
c   Choix du fichier de demarrage a modifier
c=======================================================================

      write(*,*) 'A partir de quoi souhaitez vous creer vos etats',
     .                  ' initiaux ?'
      write(*,*) '    0 - d un fichier start_archive'
      write(*,*) '    1 - d un fichier start'
 
c-----------------------------------------------------------------------
c   Ouverture des fichiers a modifier (start ou start_archive)
c-----------------------------------------------------------------------

      DO
         read(*,*,iostat=ierr) choix_1
         if ((choix_1 /= 0).OR.(choix_1 /=1)) EXIT
      ENDDO

c     Ouverture de start_archive
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (choix_1.eq.0) then

        write(*,*) 'Creation d un etat initial a partir de'
        write(*,*) './start_archive.dat et ./start_archive.dic'
        write(*,*)
        fichnom = 'start_archive.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Pb d''ouverture du fichier ',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF
        tab0 = 50 
        Lmodif = 1

c     OU Ouverture de start
c     ~~~~~~~~~~~~~~~~~~~~~
      else
        write(*,*) 'Creation d un etat initial a partir de'
        write(*,*) './start.nc et ./startfi.nc'
        write(*,*)
        fichnom = 'start.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid_dyn)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Pb d''ouverture du fichier ',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF
 
        fichnom = 'startfi.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid_fi)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Pb d''ouverture du fichier ',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF

        tab0 = 0 
        Lmodif = 0

      endif

c-----------------------------------------------------------------------
c Lecture du tableau des parametres du run (pour la dynamique)
c-----------------------------------------------------------------------

      if (choix_1.eq.0) then

        write(*,*) 'lecture tab_cntrl START_ARCHIVE'
c
        ierr = NF_INQ_VARID (nid, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
c
      else if (choix_1.eq.1) then

        write(*,*) 'lecture tab_cntrl START'
c
        ierr = NF_INQ_VARID (nid_dyn, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid_dyn, nvarid, tab_cntrl)
#else
        ierr = NF_GET_VAR_REAL(nid_dyn, nvarid, tab_cntrl)
#endif
c
        write(*,*) 'lecture tab_cntrl STARTFI'
c
        ierr = NF_INQ_VARID (nid_fi, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid_fi, nvarid, tab_cntrl_bis)
#else
        ierr = NF_GET_VAR_REAL(nid_fi, nvarid, tab_cntrl_bis)
#endif
c
        do i=1,50
          tab_cntrl(i+50)=tab_cntrl_bis(i)
        enddo
      write(*,*) 'Impression de tab_cntrl', tab_cntrl
      do i=1,100
        write(*,*) i,tab_cntrl(i)
      enddo
      
      endif
c-----------------------------------------------------------------------
c		Initialisation des constantes dynamique
c-----------------------------------------------------------------------

      kappa = tab_cntrl(9) 
      etot0 = tab_cntrl(12)
      ptot0 = tab_cntrl(13)
      ztot0 = tab_cntrl(14)
      stot0 = tab_cntrl(15)
      ang0 = tab_cntrl(16)
      write(*,*) "kappa,etot0,ptot0,ztot0,stot0,ang0"
      write(*,*) kappa,etot0,ptot0,ztot0,stot0,ang0

c-----------------------------------------------------------------------
c   Lecture du tab_cntrl et initialisation des constantes physiques
c  - pour start:  Lmodif = 0 => pas de modifications possibles
c                  (modif dans le tabfi de readfi + loin)
c  - pour start_archive:  Lmodif = 1 => modifications possibles
c-----------------------------------------------------------------------
      if (choix_1.eq.0) then
         call tabfi (nid,Lmodif,tab0,day_ini,lllm,p_rad,
     .            p_omeg,p_g,p_mugaz,p_daysec,time)
      else if (choix_1.eq.1) then
         call tabfi (nid_fi,Lmodif,tab0,day_ini,lllm,p_rad,
     .            p_omeg,p_g,p_mugaz,p_daysec,time)
      endif

      rad = p_rad
      omeg = p_omeg
      g = p_g
      mugaz = p_mugaz
      daysec = p_daysec
      write(*,*) 'aire',aire


c=======================================================================
c  INITIALISATIONS DIVERSES
c=======================================================================

      day_step=180
      CALL defrun_new( 99, .TRUE. )
      CALL iniconst 
      CALL inigeom
      idum=-1
      idum=0

c Initialisation coordonnees /aires
c -------------------------------

      latfi(1)=rlatu(1)
      lonfi(1)=0.
      DO j=2,jjm
         DO i=1,iim
            latfi((j-2)*iim+1+i)=rlatu(j)
            lonfi((j-2)*iim+1+i)=rlonv(i)
         ENDDO
      ENDDO
      latfi(ngridmx)=rlatu(jjp1)
      lonfi(ngridmx)=0.
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,aire,airefi)


c=======================================================================
c  Lecture des fichiers (start ou start_archive)
c=======================================================================

      if (choix_1.eq.0) then

        write(*,*) 'lecture Fichier START_ARCHIVE'
        CALL lect_start_archive(date,tsurf,tsoil,emis,q2,
     .   t,ucov,vcov,ps,co2ice,teta,phisold_newgrid,q,qsurf,nid)
        
        ierr= NF_CLOSE(nid)

      else if (choix_1.eq.1) then !  c'est l'appel a tabfi de phyeta0 qui
                                  !  permet de changer les valeurs du 
                                  !  tab_cntrl Lmodif=1
        tab0=0
        Lmodif=1 ! Lmodif set to 1 to allow modifications in phyeta0                           
        write(*,*) 'lecture Fichier START'
        fichnom = 'start.nc'
        CALL dynetat0(fichnom,nqmx,vcov,ucov,teta,q,masse,
     .       ps,phis,time)

        write(*,*) 'lecture Fichier STARTFI'
        fichnom = 'startfi.nc'
        CALL phyetat0 (fichnom,tab0,Lmodif,nsoilmx,nqmx,
     .        day_ini,time,
     .        tsurf,tsoil,emis,q2,qsurf,co2ice)
        
        do i=1,ngridmx
          albfi(i) = albedodat(i)
          ithfi(i) = inertiedat(i)
        enddo

      else 
        CALL exit(1)
      endif

      dtvr   = daysec/FLOAT(day_step)
      dtphys   = dtvr * FLOAT(iphysiq)

c=======================================================================
c   lecture topographie, albedo, inertie thermique, relief sous-maille
c=======================================================================

      if (choix_1.ne.1) then  ! pour ne pas avoir besoin du fichier 
                              ! surface.dat dans le cas des start

c do while((relief(1:3).ne.'mol').AND.(relief(1:3).ne.'pla'))
c       write(*,*)
c       write(*,*) 'choix du relief (mola,pla)'
c       write(*,*) '(Topographie MGS MOLA, plat)'
c       read(*,fmt='(a3)') relief
        relief="mola"
c     enddo

      CALL datareadnc(relief,phis,alb,ith,zmeaS,zstdS,zsigS,zgamS,
     .          ztheS)

      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,phis,phisfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ith,ithfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zmeaS,zmea)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zstdS,zstd)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zsigS,zsig)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zgamS,zgam)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ztheS,zthe)

      endif

c=======================================================================
c 
c=======================================================================
 
 888  continue

      write(*,*)
      write(*,*)
      write(*,*) 'Autres Modifications possibles :'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*)
      write(*,*) 'flat : no topography ("aquaplanet")'
      write(*,*) 'bilball : albedo, inertie thermique uniforme'
      write(*,*) 'coldspole : sous sol froid et haut albedo au pole sud'
      write(*,*) 'q=0 : ALL tracer =zero'
      write(*,*) 'q=x : give a specific uniform value to one tracer'
      write(*,*) 'ini_q : tracers initialisation for chemistry, water an
     $d ice   '
      write(*,*) 'ini_q-H2O : tracers initialisation for chemistry and 
     $ice '
      write(*,*) 'ini_q-iceH2O : tracers initialisation for chemistry on
     $ly '
      write(*,*) 'watercapn : H20 ice sur la calotte permanente nord '
      write(*,*) 'watercaps : H20 ice sur la calotte permanente sud '
      write(*,*) 'wetstart  : start with a wet atmosphere'
      write(*,*) 'isotherm : Temperatures isothermes et vents nuls '
      write(*,*) 'co2ice=0 : elimination des calottes polaires de CO2'
      write(*,*) 'ptot : pression totale'


      do while(modif(1:1).ne.'hello')
        write(*,*)
        write(*,*) 'modification a faire ?'
        write(*,*) '   (entrer le mot cle ou return pour finir)'
        write(*,*)

        read(*,fmt='(a20)') modif
        if (modif(1:1) .eq. ' ') goto 999
        write(*,*)
        write(*,*) modif(1:lnblnk(modif)) , ' : '

c       'flat : no topography ("aquaplanet")'
c       -------------------------------------
        if (modif(1:lnblnk(modif)) .eq. 'flat') then
          CALL initial0(ip1jmp1,z_reel)
          CALL multscal(ip1jmp1,z_reel,g,phis)
          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,phis,phisfi)
          write(*,*) 'topography set to zero.'
          write(*,*) 'WARNING : the subgrid topography parameters',
     &    ' were not set to zero ! => set calllott to F'                    
          if (choix_1.ne.0) then
              write(*,*)'WARNING : reading start, the pressure is not'
     &                  ,' corrected for topography change'
              write(*,*) '=> you should start with start_archive !'
          end if

c       bilball : albedo, inertie thermique uniforme
c       --------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'bilball') then
          write(*,*) 'Valeur constante pour albedo et iner.therm:'
          write(*,*) 'Nouvelle valeur pour albedo (ex: 0.25) ?'
 101      read(*,*,iostat=ierr) alb_bb
          if(ierr.ne.0) goto 101
          write(*,*)
          write(*,*) ' albedo uniforme (nouvelle valeur):',alb_bb
          write(*,*)

          write(*,*) 'Nouvelle valeur pour iner.therm (ex: 247) ?'
 102      read(*,*,iostat=ierr) ith_bb
          if(ierr.ne.0) goto 102
          write(*,*) 'iner.therm uniforme (nouvelle valeur):',ith_bb
          DO j=1,jjp1
             DO i=1,iip1
                alb(i,j) = alb_bb
                ith(i,j) = ith_bb
             END DO
          END DO
          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ith,ithfi)
          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)

c       coldspole : sous-sol de la calotte sud toujours froid
c       -----------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'coldspole') then
          write(*,*)'nouvelle valeur de la temperature du',
     &   'sous sol de la calotte permanente sud ? (ex: 141 K)'
 103      read(*,*,iostat=ierr) tsud
          if(ierr.ne.0) goto 103
          write(*,*)
          write(*,*) ' nouvelle valeur de la temperature:',tsud
c         nouvelle temperature sous la calotte permanente
          do l=2,nsoilmx
               tsoil(ngridmx,l) =  tsud
          end do


          write(*,*)'nouvelle valeur de l albedo',
     &   'de la calotte permanente sud ? (ex: 0.75 K)'
 104      read(*,*,iostat=ierr) albsud
          if(ierr.ne.0) goto 104
          write(*,*)

c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         Option 1:  Seul l'albedo du pole est modifié :    
          albfi(ngridmx)=albsud
          write(*,*) 'ig=',ngridmx,'   albedo perennial cap ',
     &    albfi(ngridmx)

c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~é 
c          Option 2 A haute resolution : coordonnée de la vrai calotte ~    
c           DO j=1,jjp1
c             DO i=1,iip1
c                ig=1+(j-2)*iim +i
c                if(j.eq.1) ig=1
c                if(j.eq.jjp1) ig=ngridmx
c                if ((rlatu(j)*180./pi.lt.-84.).and.
c     &            (rlatu(j)*180./pi.gt.-91.).and.
c     &            (rlonv(i)*180./pi.gt.-91.).and.
c     &            (rlonv(i)*180./pi.lt.0.))         then
cc    albedo de la calotte permanente fixe a albsud
c                   alb(i,j)=albsud
c                   write(*,*) 'lat=',rlatu(j)*180./pi,
c     &                      ' lon=',rlonv(i)*180./pi
cc     fin de la condition sur les limites de la calotte permanente
c                end if
c             ENDDO
c          ENDDO
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c         CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)


c       ptot : Modification de la pression totale glace + atm actuelle
c       --------------------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'ptot') then

c         calcul de la pression totale glace + atm actuelle
          patm=0.
          airetot=0.
          pcap=0.
          DO j=1,jjp1
             DO i=1,iim
                ig=1+(j-2)*iim +i
                if(j.eq.1) ig=1
                if(j.eq.jjp1) ig=ngridmx
                patm = patm + ps(i,j)*aire(i,j)
                airetot= airetot + aire(i,j)
                pcap = pcap + aire(i,j)*co2ice(ig)*g
             ENDDO
          ENDDO
          ptoto = pcap + patm

          print*,'Pression totale au sol actuelle (co2 ice + atm) ',
     &       ptoto/airetot

          print*,'nouvelle valeur?'
          read(*,*) ptotn
          ptotn=ptotn*airetot
          patmn=ptotn-pcap
          print*,'ptoto,patm,ptotn,patmn'
          print*,ptoto,patm,ptotn,patmn
          print*,'Facteur mult de la pression (atm only)', patmn/patm
          do j=1,jjp1
             do i=1,iip1
                ps(i,j)=ps(i,j)*patmn/patm
             enddo
          enddo

c        Correction pour la conservation des traceurs
         yes=' '
         do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*) 'Do you wish to conserve tracer total mass (y)',
     &              ' or tracer mixing ratio (n) ?'
             read(*,fmt='(a)') yes
         end do

         if (yes.eq.'y') then
           write(*,*) 'OK : conservation of tracer total mass'
           DO iq =1, nqmx
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=q(i,j,l,iq)*patm/patmn
                  ENDDO
               ENDDO
             ENDDO
           ENDDO
          else
            write(*,*) 'OK : conservation of tracer mixing ratio'
          end if

c       q=0 : traceurs a zero
c       ---------------------
        else if (modif(1:lnblnk(modif)) .eq. 'q=0') then
c          mise a 0 des q (traceurs)
          write(*,*) 'Traceurs mis a 0 (1.E-30 en fait)'
           DO iq =1, nqmx
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=1.e-30
                  ENDDO
               ENDDO
             ENDDO
           ENDDO

c          mise a 0 des qsurf (traceurs a la surface)
           DO iq =1, nqmx
             DO ig=1,ngridmx
                 qsurf(ig,iq)=0.
             ENDDO
           ENDDO

c       q=x : initialise tracer manually 
c       --------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'q=x') then
             write(*,*) 'Which tracer do you want to modify ?'
             write(*,*) '(choose between 1 and ',nqmx,')'
             read(*,*) iq 
             write(*,*)'mixing ratio of tracer #',iq, ' ? (kg/kg)'
             read(*,*) val
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=val
                  ENDDO
               ENDDO
             ENDDO
             write(*,*) 'SURFACE value of tracer #', iq , ' ? (kg/m2)'
             read(*,*) val
             DO ig=1,ngridmx
                 qsurf(ig,iq)=val
             ENDDO

c       ini_q : traceurs initialises pour la chimie
c       -----------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'ini_q') then
c         For more than 32 layers, possible to initiate thermosphere only     
          thermo=0
          yes=' '
          if (llm.gt.32) then 
            do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*)'',
     &     'initialisation for thermosphere only? (y/n)'
            read(*,fmt='(a)') yes
            if (yes.eq.'y') then
            thermo=1 
            else
            thermo=0
            endif
            enddo  
          endif
          
              call inichim_newstart(q,ps,sig,nqmx,latfi,lonfi,airefi,
     $   thermo)
          write(*,*) 'Especes chimiques initialisees'

        if (thermo.eq.0) then 
c          mise a 0 des qsurf (traceurs a la surface)
           DO iq =1, nqmx
             DO ig=1,ngridmx
                 qsurf(ig,iq)=0.
             ENDDO
           ENDDO
        endif   

c       ini_q-H2O : idem sauf le dernier traceur (H2O)
c       -----------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'ini_q-H2O') then
          ! for more than 32 layers, possible to initiate thermosphere only     
          thermo=0
          yes=' '
          if(llm.gt.32) then
            do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*)'',
     &      'initialisation for thermosphere only? (y/n)'
            read(*,fmt='(a)') yes
            if (yes.eq.'y') then 
            thermo=1 
            else
            thermo=0
            endif
            enddo
          endif
              call inichim_newstart(q,ps,sig,nqmx-1,latfi,lonfi,airefi,
     $   thermo)
          write(*,*) 'Esp. chim. initialisees sauf derniere (H2O)'

        if (thermo.eq.0) then 
c          mise a 0 des qsurf (traceurs a la surface)
           DO iq =1, nqmx-1
             DO ig=1,ngridmx
                 qsurf(ig,iq)=0.
             ENDDO
           ENDDO
        endif

c       ini_q-iceH2O : idem sauf ice et H2O
c       -----------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'ini_q-iceH2O') then
c         For more than 32 layers, possible to initiate thermosphere only     
          thermo=0
          yes=' '
          if(llm.gt.32) then
            do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*)'',
     &      'initialisation for thermosphere only? (y/n)'
            read(*,fmt='(a)') yes
            if (yes.eq.'y') then 
            thermo=1 
            else
            thermo=0
            endif
            enddo
          endif
     
         call inichim_newstart(q,ps,sig,nqmx-2,latfi,lonfi,airefi,
     $   thermo)
          write(*,*) 'Esp. chim. initialisees sauf ice et H2O'

        if (thermo.eq.0) then 
c          mise a 0 des qsurf (traceurs a la surface)
           DO iq =1, nqmx-2
             DO ig=1,ngridmx
                 qsurf(ig,iq)=0.
             ENDDO 
           ENDDO
        endif
c      wetstart : wet atmosphere with a north to south gradient
c      --------------------------------------------------------
       else if (modif(1:lnblnk(modif)) .eq. 'wetstart') then
          DO l=1,llm
            DO j=1,jjp1
              DO i=1,iip1
                q(i,j,l,nqmx)=150.e-6 * (rlatu(j)+pi/2.) / pi
              ENDDO
            ENDDO
          ENDDO

         write(*,*) 'Water mass mixing ratio at north pole='
     *               ,q(1,1,1,nqmx)
         write(*,*) '---------------------------south pole='
     *               ,q(1,jjp1,1,nqmx)

c      watercapn : H20 ice sur la calotte permanente nord
c      --------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'watercapn') then
           do ig=1,ngridmx
             j=(ig-2)/iim +2
              if(ig.eq.1) j=1
              if (rlatu(j)*180./pi.gt.80.) then
                   qsurf(ig,nqmx)=1.e5
                   write(*,*) 'ig=',ig,'    H2O ice mass (kg/m2)= ',
     &              qsurf(ig,nqmx)
                   write(*,*)'     ==> Ice mesh South boundary (deg)= ',
     &              rlatv(j)*180./pi
                end if
           enddo

c      watercaps : H20 ice sur la calotte permanente sud
c      -------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'watercaps') then
           do ig=1,ngridmx
               j=(ig-2)/iim +2
               if(ig.eq.1) j=1
               if (rlatu(j)*180./pi.lt.-80.) then
                   qsurf(ig,nqmx)=1.e5
                   write(*,*) 'ig=',ig,'   H2O ice mass (kg/m2)= ',
     &              qsurf(ig,nqmx)
                   write(*,*)'     ==> Ice mesh North boundary (deg)= ',
     &              rlatv(j-1)*180./pi
               end if
           enddo

c       isotherm : Temperatures isothermes et vents nuls
c       ------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'isotherm') then

          write(*,*)'Isothermal temperature of the atmosphere, 
     &           surface and subsurface'
          write(*,*) 'Value of THIS temperature ? :'
 203      read(*,*,iostat=ierr) Tiso
          if(ierr.ne.0) goto 203

          do ig=1, ngridmx
            tsurf(ig) = Tiso
          end do
          do l=2,nsoilmx
            do ig=1, ngridmx
              tsoil(ig,l) = Tiso
            end do
          end do
          flagiso=.true.
          call initial0(llm*ip1jmp1,ucov)
          call initial0(llm*ip1jm,vcov)
          call initial0(ngridmx*(llm+1),q2)

c       co2ice=0 : ellimination des calottes polaires de CO2'
c       ------------------------------------------------
        else if (modif(1:lnblnk(modif)) .eq. 'co2ice=0') then
           do ig=1,ngridmx
              co2ice(ig)=0
              emis(ig)=emis(ngridmx/2)
           end do
        else
          goto 888
        end if
      end do

 999  continue

 
c=======================================================================
c   Correction de la pression pour nouvelle grille (menu 0)
c=======================================================================

      if (choix_1.eq.0) then
        r = 1000.*8.31/mugaz

        do j=1,jjp1
          do i=1,iip1
             ps(i,j) = ps(i,j) * 
     .            exp((phisold_newgrid(i,j)-phis(i,j)) /
     .                                  (t(i,j,1) * r))
          end do
        end do
  
c periodicite de ps en longitude
        do j=1,jjp1
          ps(1,j) = ps(iip1,j)
        end do
      end if

c=======================================================================
c=======================================================================

c=======================================================================
c    Initialisation de la physique / ecriture de newstartfi :
c=======================================================================


      CALL inifilr 
      CALL pression(ip1jmp1, ap, bp, ps, p3d)

c-----------------------------------------------------------------------
c   Initialisation  pks:
c-----------------------------------------------------------------------

      CALL exner_hyb(ip1jmp1, ps, p3d, beta, pks, pk, pkf)
! Calcul de la temperature potentielle teta

      if (flagiso) then
          DO l=1,llm
             DO j=1,jjp1
                DO i=1,iim
                   teta(i,j,l) = Tiso * cpp/pk(i,j,l)
                ENDDO
                teta (iip1,j,l)= teta (1,j,l)
             ENDDO
          ENDDO
      else if (choix_1.eq.0) then
         DO l=1,llm
            DO j=1,jjp1
               DO i=1,iim
                  teta(i,j,l) = t(i,j,l) * cpp/pk(i,j,l)
               ENDDO
               teta (iip1,j,l)= teta (1,j,l)
            ENDDO
         ENDDO
      endif

C Calcul intermediaire
c
      if (choix_1.eq.0) then
         CALL massdair( p3d, masse  )
c
         print *,' ALPHAX ',alphax

         DO  l = 1, llm
           DO  i    = 1, iim
             xppn(i) = aire( i, 1   ) * masse(  i     ,  1   , l )
             xpps(i) = aire( i,jjp1 ) * masse(  i     , jjp1 , l )
           ENDDO
             xpn      = SUM(xppn)/apoln
             xps      = SUM(xpps)/apols
           DO i   = 1, iip1
             masse(   i   ,   1     ,  l )   = xpn
             masse(   i   ,   jjp1  ,  l )   = xps
           ENDDO
         ENDDO
      endif
      phis(iip1,:) = phis(1,:)

      CALL inidissip ( lstardis, nitergdiv, nitergrot, niterh,
     *                tetagdiv, tetagrot , tetatemp  )
      itau=0
      if (choix_1.eq.0) then
         day_ini=int(date)
      endif
c
      CALL geopot  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )

      CALL caldyn0( itau,ucov,vcov,teta,ps,masse,pk,phis ,
     *                phi,w, pbaru,pbarv,day_ini+time )
c     CALL caldyn
c    $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
c    $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, day_ini )

      CALL dynredem0("restart.nc",day_ini,anneeref,phis,nqmx)
      CALL dynredem1("restart.nc",0.0,vcov,ucov,teta,q,nqmx,masse,ps) 
C
C Ecriture etat initial physique
C
      call physdem1("restartfi.nc",lonfi,latfi,nsoilmx,nqmx,
     .              dtphys,float(day_ini),
     .              time,tsurf,tsoil,co2ice,emis,q2,qsurf,
     .              airefi,albfi,ithfi,zmea,zstd,zsig,zgam,zthe)

c=======================================================================
c	Formats 
c=======================================================================

   1  FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lllm =',i4,2x,'lue sur le fichier demar
     *rage est differente de la valeur parametree llm =',i4//)

      end
