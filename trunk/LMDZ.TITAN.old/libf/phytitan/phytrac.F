      SUBROUTINE phytrac (firstcall,lastcall,
     .                   nqmax,nmicro,ptimestep,appkim,dtkim,
     .                   pplev,pplay,delp,ptemp,pmu0,pfract,pdecli,
     .                   lonsol,
     .                   pu,pv,pzlev,pzlay,ftsol,
     .                   tr_seri,qaer,d_tr_mph,d_tr_kim,
     .                   fclat,reservoir)

c======================================================================
c    S. Lebonnois, mai 2008
c
c  Arguments:
c
c firstcall----input-L-variable logique indiquant le premier passage
c lastcall-----input-L-variable logique indiquant le dernier passage
c nqmax--------input-I-nombre de traceurs (total)
c nmicro-------input-I-nombre de traceurs microphysiques !! doivent etre toujours en premiers!!
c ptimestep----input-R-pas d integration pour la physique (seconde)
c appkim-------input-I-appel a la chimie
c dtkim--------input-R-pas de temps chimique (seconde)
c pplev--------input-R-pression pour chaque inter-couche (en Pa)
c pplay--------input-R-pression pour chaque couche (en Pa)
c delp---------input-R-epaisseur d'une couche (en Pa)
c ptemp--------input-R-temperature (K)
c pmu0---------input-R-cos angle zenithal
c pfract-------input-R-fractional day
c pdecli-------input-R-declinaison en radian
c lonsol-------input-R-longitude solaire en radian
c pu-----------input-R-vitesse dans la direction X (de O a E) en m/s (1ere couche)
c pv-----------input-R-vitesse Y (de S a N) en m/s                   (1ere couche)
c pzlev--------input-R-altitude pour chaque inter-couche (en m)
c pzlay--------input-R-altitude pour chaque couche (en m)
c ftsol--------input-R-temperature au sol (en K)
c tr_seri------input-R-mass mixing ratio traceurs (kg/kg) 
c d_tr_mph----output-R-tendance microphysique de "qx" (kg/kg/s)
c d_tr_kim----output-R-tendance chimique de "qx" (kg/kg/s)
c fclat--------output-R-flux de chaleur latente d'evaporation du reservoir CH4 (J/m2/s)
c reservoir----outpur-R-un reservoir de surface !!! (m)
c======================================================================
      USE infotrac_phy, ONLY: tname
      use dimphy
      USE common_mod, only: rmcbar,xfbar,ncount,
     &      flxesp_i,tau_drop,tau_aer,solesp,precip,
     &      evapch4,occcld_m,occcld,satch4,satc2h6,satc2h2,rmcloud
      USE moyzon_mod
      USE write_field_phy
      USE logic_mod, ONLY: moyzon_ch,moyzon_mu
      IMPLICIT none
#include "dimensions.h"
#include "clesphys.h"
#include "YOMCST.h"
#include "microtab.h"
#include "varmuphy.h"
#include "itemps.h"

c======================================================================
c Variables argument:
c
      LOGICAL firstcall,lastcall
      INTEGER nqmax,nmicro,nlat,appkim
      REAL ptimestep,dtkim
      REAL pplev(klon,klev+1),pplay(klon,klev+1),delp(klon,klev)
      REAL ptemp(klon,klev)
      REAL pmu0(klon), pfract(klon), pdecli, lonsol
      REAL pu(klon),pv(klon)
      REAL pzlev(klon,klev+1),pzlay(klon,klev)
      REAL ftsol(klon)
      REAL tr_seri(klon,klev,nqmax)
      REAL qaer(klon,klev,nqmax)
      REAL d_tr_mph(klon,klev,nqmax),d_tr_kim(klon,klev,nqmax)
      REAL fclat(klon)
      REAL reservoir(klon)

c======================================================================
c Local variables
      REAL qaer0(klon,klev,nqmax)
      REAL prec(klon,5)

      REAL rcloud(klon,klev,nrad),xfrac(klon,klev,4)

      REAL vcl,nuc,xgsn,xmsn,xesn,xasn


      ReAL gaz1(klon,klev),gaz2(klon,klev),gaz3(klon,klev)

      REAL socccld

c grandeurs en moyennes zonales
      REAL zplev(klon,klev+1),zplay(klon,klev)
      REAL zzlev(klon,klev+1),zzlay(klon,klev)
      REAL ztemp(klon,klev), delpbar(klon,klev)
      real temp_eq(klev),press_eq(klev)
      REAL qaer0bar(klon,klev,nqmax)   ! et non nmicro... Permet nmicro=0.
      REAL zdqmufi(klon,klev,nqmax)
      REAL ychim(klon,klev,nqmax-nmicro) 
c La saturation n est calculee qu une seule fois: sauvegarde qysat
c La chimie n est pas calculee tous les pas, il faut donc
c                      sauvegarder les sorties de la chimie
      REAL,save,allocatable :: qysat(:,:),pdyfi(:,:,:) 
      
      character*10 nomqy(nqmax-nmicro+1)
      integer      i,j,k,l,iq,ig0
      
c    indice des esp chimiques utilisees dans la microfi  
      integer icldch4,icldc2h6,icldc2h2
      save icldch4,icldc2h6,icldc2h2
     
      real fte,ftm,Lvch4

      REAL tmp,ex,kmin,kmax,dqsq
      REAL dqch4

c======================================================================
c======================================================================

      if (firstcall) then
       allocate(qysat(klev,nqmax-nmicro),pdyfi(klon,klev,nqmax-nmicro))

c  -------- Quelques verifications au demarrage sur les tailles des tableaux.
         IF (microfi.ge.1) then
c        Faire de la microphysique sans traceurs... bon courage !
           if (nmicro.le.0) then
             print*,"aLeRtE cRiTiQuE !!!"
             print*,"Vous faites de la microphysique sans traceurs"
             print*,"microphysique..."
             print*,"Je m'arrete et vous laisse reflechir !"
             stop
           endif
c        Nombre de traceurs incompatibles avec la microphysique. 
           if ((nmicro.ne.ntype*nrad).and.(clouds.eq.1)) then
             print*,"aLeRtE cRiTiQuE !!!"
             print*,"Nb trac imcompatible avec la microphysique."
             print*,nmicro,ntype*nrad
             stop
           endif
           if ((nmicro.ne.nrad).and.(clouds.eq.0)) then
             print*,"aLeRtE cRiTiQuE !!!"
             print*,"Nb trac imcompatible avec la microphysique."
             print*,nmicro,nrad
             stop
           endif
         ENDIF

      endif  ! firstcall

c RAZ des sorties : les moyennes se font directement dans IOIPSL :
c
          flxesp_i(:,:,:) = 0.
          tau_drop(:,:)   = 0.
          tau_aer(:,:,:)  = 0.
          solesp(:,:,:)   = 0.
          precip(:,:)     = 0.   ! c'est uniquement une sortie en um/s
c
          prec(:,:)       = 0.   ! c'est la variable temporaire des precipitions de la microfi
                                 ! prec est en m (metre precipitable)

c-----------------------------------------------------------------------
c   convertion moyennes zonales et changement d unites pour microphy
c   ---------------------------------

c     print*,'CONVERSION 2D ET CHANGEMENT UNITES (PHYTRAC)'

c   -------------------
c   Gestion de la temperature et de la pression :
c   Utilisation des moyennes zonales:

c   soit la chimie est active, soit la microphysique se fait en 2D.
      IF (chimi.or.microfi.eq.1) THEN
        zplev(:,:) = zplevbar(:,:)
        zplay(:,:) = zplaybar(:,:)
        zzlev(:,:) = zzlevbar(:,:)
        zzlay(:,:) = zzlaybar(:,:)
        ztemp(:,:) = ztfibar(:,:)
        ychim = 0.0
      ENDIF

c  Si la microphysique est faite en 2D:
      IF (microfi.eq.1) THEN
        DO l=1,llm
          DO i = 1, klon
            delpbar(i,l) = zplevbar(i,l) - zplevbar(i,l+1)
          ENDDO
        ENDDO
c   Traceurs microphysiques: passage en extensif: n/kg --> n/m^2
        DO iq=1,nmicro
         qaer(:,:,iq) = zqfibar(:,:,iq)*delpbar(:,:)/RG
        DO l=1,klev
          DO i = 1, klon
              if (qaer(i,l,iq).lt.0.) then
        print*,"NEGS ICI ICI !!!!",qaer(i,l,iq),i,l,iq
                qaer(i,l,iq)=0.
c          stop
              endif
              if (delpbar(i,l).lt.0.) then
        print*,"NEGS DELP ICI !!!!",i,l,iq,delpbar(i,l)
           stop
              endif
          ENDDO
        ENDDO
         qaer0(:,:,iq)= tr_seri(:,:,iq)*delp(:,:)/RG
         qaer0bar(:,:,iq) = qaer(:,:,iq)
        ENDDO
      ENDIF

c  Si la microphysique est faite en 3D:
      IF (microfi.eq.2) THEN
        zplev(:,:) = pplev(:,:)
        zplay(:,:) = pplay(:,:)
        zzlev(:,:) = pzlev(:,:)
        zzlay(:,:) = pzlay(:,:)
        ztemp(:,:) = ptemp(:,:)    
c   Traceurs microphysiques: passage en extensif: n/kg --> n/m^2
        DO iq=1,nmicro
         qaer(:,:,iq) = tr_seri(:,:,iq)*delp(:,:)/RG
         qaer0(:,:,iq)= tr_seri(:,:,iq)*delp(:,:)/RG
        ENDDO
      ENDIF

      do l=1,llm
	 temp_eq  = tmoy
	 press_eq = playmoy/100. ! en mbar
      enddo

c   -------------------
c    Extraction des gaz pour les nuages
      IF ((microfi.ge.1).and.(clouds.eq.1)) THEN

c     recuperation des indices des gaz qui nous interesse       
      if (firstcall) then
          icldch4=-1
          icldc2h6=-1
          icldc2h2=-1
          do i=1,nqmax
            if (tname(i).eq."CH4") then
              icldch4=i
c              ich4=i
            elseif (tname(i).eq."C2H6") then
              icldc2h6=i
            elseif (tname(i).eq."C2H2") then
              icldc2h2=i
            endif
          enddo
          if (icldch4 .eq.-1 .or. 
     &        icldc2h6.eq.-1 .or.
     &        icldc2h2.eq.-1 ) then
            print*, "Sacrebleu !!!"
            print*, "Vous voulez faire des nuages sans gaz."
            print*, "Mais vous etes inconscient. Je vais m'arreter la"
            print*, "pour vous laisser reflechir au probleme"
            STOP
          endif
      endif      ! firstcall

c     Saturation et fraction molaire CLOUD 
c     Calcul des saturations pour les esp chimique de la muphy des nuages.
c     On le fait ici pour les sortir dans physiq.F sans avoir a surcharger la routine. 
c     Elles passent ensuite dans un common pour passer dans les I/O.

        DO l=1,llm
          DO i = 1, klon
            call ch4sat(ptemp(i,l),pplay(i,l),tmp) !tmp en kg/kg !
            satch4(i,l) = tr_seri(i,l,icldch4)/(tmp*28./16.)

            call c2h6sat(ptemp(i,l),pplay(i,l),tmp)
            satc2h6(i,l) =tr_seri(i,l,icldc2h6)/(tmp*28./30.)

            call c2h2sat(ptemp(i,l),pplay(i,l),tmp)
            satc2h2(i,l) =tr_seri(i,l,icldc2h2)/(tmp*28./26.)
          ENDDO
        ENDDO

c   Copie des gaz (en 3D)  <== UNIQUEMENT SI ON FAIT DES NUAGES
        if (moyzon_mu) then
         gaz1(:,:) = zqfibar(:,:,icldch4)
         gaz2(:,:) = zqfibar(:,:,icldc2h6)
         gaz3(:,:) = zqfibar(:,:,icldc2h2)
        else
         gaz1(:,:) = tr_seri(:,:,icldch4)
         gaz2(:,:) = tr_seri(:,:,icldc2h6)
         gaz3(:,:) = tr_seri(:,:,icldc2h2)
        endif
        
      endif      ! microfi.ge.1 + clouds.eq.1
c   -------------------

c AUTRES TRACEURS
      
      if (nqmax.gt.nmicro) then
       do iq=nmicro+1,nqmax
        if (moyzon_ch) then
          ychim(:,:,iq-nmicro) = zqfibar(:,:,iq)
        else
          ychim(:,:,iq-nmicro) = tr_seri(:,:,iq)
        endif
        nomqy(iq-nmicro) = tname(iq)
c        print*,iq-nmicro,nomqy(iq-nmicro)
       enddo
       nomqy(nqmax-nmicro+1) = "HV"
      endif

c-----------------------------------------------------------------------
c   initialisation des qysat au premier appel:
c   ---------------------------------

c!! ATTENTION, qysat pris uniquement a l'equateur
c!!  justifie puisque dans cette region, les var de t et p sont faibles...

      if (firstcall .and. chimi .and.(nqmax.gt.nmicro)) then
           call inicondens(nqmax-nmicro,press_eq,temp_eq,nomqy,qysat)
      endif

c-----------------------------------------------------------------------
c     Appel de la microphysique   en 2D/3D !!!!!!
c    --------------------------

       IF(firstcall) THEN 
        print*,'MICROPHYSIQUE ',MICROFI
       ENDIF

c       call begintime(tt0)
       IF (MICROFI.eq.0) THEN
c        PAS DE MICROPHYSIQUE :
         IF (firstcall) THEN
          print*,'MICROPHYSIQUE OFF-LINE',MICROFI
         ENDIF
       ELSE
         zdqmufi = 0.  ! ne sert que pour chimi pour condensation
         call muphys(klon,
     &        zplev,zplay,zzlev,zzlay,
     &        ztemp,qaer,gaz1,gaz2,gaz3,
     &        nmicro,ptimestep,
     &        pmu0,pfract,
c -------- sorties diagnostiques
     &        flxesp_i,
     &        tau_drop,tau_aer,
     &        solesp,prec)

c    NOTES :
c    Ici toutes nos sorties sont des champs 3D...(meme les diagnostiques)
c    On a rien a faire mis a part copier les dq dans les d_tr

       ENDIF
c       call endtime(tt0,tt1)
c       ttmuphys=ttmuphys+tt1
       
c-----------------------------------------------------------------------
c     Gestion des sources
c    -------------
c
       IF (clouds.eq.1) THEN
        IF (microfi.eq.1) THEN
c        On repasse les gaz en 3D si on a fait de la microphysique en 2D
         gaz1(:,:)=gaz1(:,:)*tr_seri(:,:,icldch4)/zqfibar(:,:,icldch4)
         gaz2(:,:)=gaz2(:,:)*tr_seri(:,:,icldc2h6)/zqfibar(:,:,icldc2h6)
         gaz3(:,:)=gaz3(:,:)*tr_seri(:,:,icldc2h2)/zqfibar(:,:,icldc2h2)
        ENDIF
c       Mise a jour du reservoir de CH4 (ie : seul le CH4 remplit le reservoir)
        DO i=1,klon
          reservoir(i) = reservoir(i)+prec(i,1)
        ENDDO 

        CALL sources(klon,klev,ptimestep,z0,
     &                pu,pv,pplev,pzlay,pzlev,
     &                gaz1,gaz2,gaz3,
     &                ftsol,evapch4,reservoir) 
 
       ENDIF
c-----------------------------------------------------------------------
c     Condensation
c    -------------

      IF ((chimi).and.(nqmax.gt.nmicro)) then

c   tendance (en /s) passee sur zdqmufi(nmicro+1 a nqmax)
c        print*,'Condensation'

        do iq=1,nqmax-nmicro
           do l=1,llm
              do j=1,klon
                 if (ychim(j,l,iq).gt.qysat(l,iq)) then
           zdqmufi(j,l,nmicro+iq)= (-ychim(j,l,iq)+qysat(l,iq)) !delta y
     .                             / ptimestep                  ! / dt
                 endif
              enddo
           enddo
        enddo

      ENDIF

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c eventuellement, modif initiale de la compo
c
c   tendance (en /s) passee sur zdqmufi(nmicro+1 a nqmax)
c
c     if (firstcall .and. chimi .and.(nqmax.gt.nmicro)) then
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!remise de CH4 a 1.5%!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       do iq=1,nqmax-nmicro
c         if (nomqy(iq).eq."CH4") then
c          do l=1,llm
c             do j=1,klon
c                if (ychim(j,l,iq).le.0.015) then
c          zdqmufi(j,l,nmicro+iq)= (-ychim(j,l,iq)+0.015) !delta y
c    .                            / ptimestep                  ! / dt
c                endif
c             enddo
c          enddo
c         endif
c       enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!remise de C2H2 a 1.e-5 max !!!!!!!!!!!!!
c!!!remise de C2H6 a 3.e-5 max !!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       do iq=1,nqmax-nmicro
c         if (nomqy(iq).eq."C2H2") then
c          do l=1,llm
c             do j=1,klon
c                if (ychim(j,l,iq).gt.1.e-5) then
c          zdqmufi(j,l,nmicro+iq)= (-ychim(j,l,iq)+1.e-5) !delta y
c    .                            / ptimestep                  ! / dt
c                endif
c             enddo
c          enddo
c         endif
c         if (nomqy(iq).eq."C2H6") then
c          do l=1,llm
c             do j=1,klon
c                if (ychim(j,l,iq).gt.3.e-5) then
c          zdqmufi(j,l,nmicro+iq)= (-ychim(j,l,iq)+3.e-5) !delta y
c    .                            / ptimestep                  ! / dt
c                endif
c             enddo
c          enddo
c         endif
c       enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
c ----- commentaire de fin (mise a jour des profil de fraction molaire)
          
c-----------------------------------------------------------------------
c     Appel de la chimie
c    --------------------------

      if((appkim.eq.1).and.(chimi)) then
        print*,'On passe dans la CHIMIE'

c       do iq=1,nqmax-nmicro
c         if (nomqy(iq).eq."C2H2") then
c           print*,"C2H2top=",ychim(:,klev,iq)
c         endif
c       enddo

c Appel Chimie
c ------------
       CALL calchim(klon,nqmax-nmicro,ychim,nomqy,pdecli,lonsol,dtkim,
     .              ztemp,zplay,zplev,zzlay,zzlev,
     .              pdyfi)   
c ychim ne doit pas etre modifie, pdyfi en /s
	 
      endif
      
c-----------------------------------------------------------------------
c   retour des tendances vers 3D
c   ---------------------------------

c     call WriteField_phy('phytrac_qaer01',qaer(:,:,1),klev)

c===============================
c TRACEURS MICROPHYSIQUES
c                                        
c===============================
c ---> pas de microphysique
       IF (microfi.eq.0) THEN
         DO iq=1,nmicro
           d_tr_mph(:,:,iq)=0.
         ENDDO
       ENDIF
c===============================
c ---> microphysique 2D
c     IF (microfi.eq.1) THEN
c        DO iq=1,nmicro
c          DO l=1,llm
c            DO i=1,klon
c  on repasse le champ de traceurs en 3D (pas les tendances)
c qaer est ce qui entre dans muphy, donc la moyenne zonale
c qaer0 est la valeur initiale du champ
c qaer0bar est la moyenne zonale initiale
c la variation relative pour une bande de latitude est donc (qaer/qaer0bar)
c la nouvelle valeur en un point (3D) est donc qaer0*(qaer/qaer0bar)
c et la tendance: qaer0*(qaer/qaer0bar)-qaer0
c    un petit patch : 
c    Si la moyenne zonale au depart est "nulle" :
c    On a quand meme le droit de produire des traceurs dans la cellule.
c    On considere donc que la valeur de sortie 3D correspond a la valeur de sortie 2D.
c    Cela permet aussi entre autre d eviter les NaN pour les traceurs des nuages !
c    (au dessus de la tropo pas de nuages donc qaer(nrad+1:ntype*nrad) = 0 !!!)
c              IF (qaer0bar(i,l,iq).gt.1e-100) THEN
c                  qaer(i,l,iq) = qaer0(i,l,iq) *
c    &             qaer(i,l,iq)/qaer0bar(i,l,iq)
c              ENDIF
c        La tendance correspond a (qaer-qaer0)/ptimestep
c              d_tr_mph(i,l,iq) = (qaer(i,l,iq)-qaer0(i,l,iq))/
c    &                            ptimestep
c            ENDDO
c          ENDDO
c        ENDDO
c ---> microphysique 3D
c      ELSEIF(microfi.gt.1) THEN 
c        DO iq=1,nmicro
c          d_tr_mph(:,:,iq)=(qaer(:,:,iq)-qaer0(:,:,iq))/ptimestep
c        ENDDO
c      ENDIF   ! microfi

c      DO iq=1,nmicro
c  Traceurs microphysiques: passage en intensif: n/m^2 --> n/kg
c        d_tr_mph(:,:,iq) = d_tr_mph(:,:,iq)*RG/delp(:,:)
c      ENDDO

c===============================
c TOUT CE QUI EST AU-DESSUS NE MARCHE PAS: PLEIN DE NEGS...
c CA MARCHE EN 3D, MAIS PAS EN MOY ZON...
c===============================

c ---> microphysique 2D
      IF (microfi.eq.1) THEN
         DO iq=1,nmicro
           DO l=1,llm
             DO i=1,klon
c ici, qaer correspond a la moy zonale modifiee par la microphys.
c  Traceurs microphysiques: passage en intensif: n/m^2 --> n/kg
c en mettant ceci:
               d_tr_mph(i,l,iq) = (qaer(i,l,iq)*RG/delpbar(i,l)
     &                            -qaer0(i,l,iq)*RG/delp(i,l))/
     &                            ptimestep
c on remplace le champ 3D initial par la valeur modifiee de sa moyenne zonale
c => on remet un champ uniforme en zonal...
             ENDDO
           ENDDO
         ENDDO
c ---> microphysique 3D
       ELSEIF(microfi.gt.1) THEN 
         DO iq=1,nmicro
           d_tr_mph(:,:,iq)=(qaer(:,:,iq)-qaer0(:,:,iq))/ptimestep
c  Traceurs microphysiques: passage en intensif: n/m^2 --> n/kg
           d_tr_mph(:,:,iq) = d_tr_mph(:,:,iq)*RG/delp(:,:)
         ENDDO
       ENDIF   ! microfi

c===============================

c AUTRES TRACEURS

      if ((chimi).and.(nqmax.gt.nmicro)) then
c on passe de pdyfi (tendance chimique en /s calculee quand chimie appelee)
c          a  d_tr_kim (tendance chimique 3D en /s, passee a physiq)
c  et de zdqmufi a d_tr_mph (tendance condensation 3D en /s passee a physiq)

      DO iq=nmicro+1,nqmax
         d_tr_kim(:,:,iq) = pdyfi(:,:,iq-nmicro)
     &             *tr_seri(:,:,iq)/ychim(:,:,iq-nmicro)
         d_tr_mph(:,:,iq) = zdqmufi(:,:,iq)
     &             *tr_seri(:,:,iq)/ychim(:,:,iq-nmicro)
      ENDDO

      endif   ! chimi

c--------------------------------------------------
c  CONDENSATION VIA MICROFI
c----------------------
c La microphysique avec nuages doit se faire obligatoirement en 3D.  (FAUX ACTUELLEMENT)
c Rien n empeche de faire la chimie en 2D. Cependant pour prendre en compte la 
c condensation due a la microfi (en 3D) on recalcule la tendance finale pour
c les especes concernees (CH4, C2H6 pour le moment). 
       IF (microfi.ge.1.and.clouds.eq.1) THEN
c     condensation CH4
          d_tr_mph(:,:,icldch4) =(gaz1(:,:)-tr_seri(:,:,icldch4))
     &                            /ptimestep
c     condensation C2H6
          d_tr_mph(:,:,icldc2h6)=(gaz2(:,:)-tr_seri(:,:,icldc2h6))
     &                            /ptimestep
c     condensation C2H2
          d_tr_mph(:,:,icldc2h2)=(gaz3(:,:)-tr_seri(:,:,icldc2h2))
     &                            /ptimestep
       ENDIF

c--------------------------------------------------
c  MISE A JOUR CH4 : (pour refixer la fraction 
c                     molaire)
c--------------------------------------------------
c       IF (firstcall) THEN
c         do i=1,klon
c           do j=1,llm
c             call ch4sat(ptemp(i,j),pplay(i,j),tmp) !tmp en kg/kg !
c             tmp=0.95*0.85*tmp*28./16.
c             if (pplay(i,j).lt.20000.) then
c               dqch4 = 1.4e-2          
c             else
c               dqch4 = tmp
c             endif
c             d_tr_mph(i,j,icldch4)=(-tr_seri(i,j,icldch4)+dqch4)/ 
c     &       ptimestep
c           enddo
c         enddo
c         
c       ENDIF

c--------------------------------------------------
c  CONVERSION PRECIPITATION : 
c  en microns/secondes
c--------------------------------------------------
        precip = prec * 1.e6 / ptimestep


c--------------------------------------------------
c CALCUL DU FLUX DE CHALEUR LATENTE D EVAPORATION 
c DU METHANE 
c--------------------------------------------------
       IF (clouds.eq.1) THEN
         DO i=1,klon
           fte= (1.-ftsol(i)/305.5)
           ftm= (1.-ftsol(i)/190.5)
           if(ftm.le.1.e-3) ftm=1.e-3
           if(fte.le.1.e-3) fte=1.e-3
           Lvch4 =8.314*190.4*
     &     (7.08*ftm**0.354+10.95*1.1e-2*ftm**0.456)
     &     /mch4
           ! evapch4 en m3/m2 {ok}
           ! 425 en kg/m3
           ! Lv en J/kg       {ok}
           ! ptimestep en s   {ok}
           fclat(i)=(evapch4(i)*Lvch4*rhoi_ch4)   ! en J/m2/s
         ENDDO
       ENDIF

c--------------------------------------------------
c  GESTION DES RAYONS DE GOUTTES POUR TR
c--------------------------------------------------
       IF (clouds.eq.1) THEN

c Calcul du rayon des gouttes par bin ...
c----------------------------------------
         DO i=1,klon
           DO j=1,klev 
             DO iq=1,nrad
*      Rayon minimum selon la quantite de noyaux
               IF (qaer(i,j,iq+nrad) .le. 1.e-5) THEN
                 rcloud(i,j,iq) = 1.e-10
               ELSE
                 rcloud(i,j,iq)=
     &           ((qaer(i,j,iq+2*nrad)/qaer(i,j,iq+nrad)+
     &           qaer(i,j,iq+3*nrad)/qaer(i,j,iq+nrad) +
     &           v_e(iq))*0.75/RPI)**(1./3.)
               ENDIF 
             ENDDO
           ENDDO
         ENDDO

c .... et de leur rayon moyen total (tt bins confondu)
c------------------------------------------------------
         DO i=1,klon
           socccld=0.
           DO j=klev,1,-1    !de haut en bas pour le calcul des opacites
             vcl=0.
             nuc=0.
             xgsn=0.
             xmsn=0.
             xesn=0.
             xasn=0.
             DO iq=1,nrad
               vcl=vcl+qaer(i,j,iq+2*nrad)+
     &         qaer(i,j,iq+3*nrad)+ 
     &         qaer(i,j,iq+4*nrad)+ 
     &         v_e(iq)*qaer(i,j,iq+nrad)            ! volume des gouttes
               nuc=nuc+qaer(i,j,iq+nrad)            ! nombre de noyaux
               xgsn=xgsn+qaer(i,j,iq+nrad)*v_e(iq)  ! volume de noyaux
               xmsn=xmsn+qaer(i,j,iq+2*nrad)        ! volume de methane
               xesn=xesn+qaer(i,j,iq+3*nrad)        ! volume d' ethane
               xasn=xasn+qaer(i,j,iq+4*nrad)        ! volume d' acethylene
             ENDDO 
             IF (nuc .le.  1.e-5) THEN
               rmcloud(i,j)=1.0e-10
               xfrac(i,j,:)=0.
             ELSE
               IF(xgsn/vcl.lt.0.  .or. xgsn/vcl.gt.1.001)
     &         print*, 'PB AVEC XFRAC:', i,j,xgsn,vcl 
               rmcloud(i,j)=          ! rayon moyen des gouttes
     &         (vcl/nuc*0.75/RPI)**(1./3.)
               xfrac(i,j,1)=xgsn/vcl         ! fraction volumique noyau/goutte
               xfrac(i,j,2)=xmsn/vcl         ! fraction volumique CH4/goutte
               xfrac(i,j,3)=xesn/vcl         ! fraction volumique C2H6/goutte
               xfrac(i,j,4)=xasn/vcl         ! fraction volumique C2H2/goutte
c              calcul du rayon moyen (moyenne temporelle)
               rmcbar(i,j)=rmcbar(i,j)+rmcloud(i,j)
               xfbar(i,j,:)=xfbar(i,j,:)+xfrac(i,j,:)
               ncount(i,j) = ncount(i,j)+1
             ENDIF
             socccld=socccld+RPI*(rmcloud(i,j)**2.)*nuc
             occcld(i,j)=socccld
           ENDDO
         ENDDO
c
c      OCCCLD
c      Calcul le nombre d occurence d un nuage 
c      d opacite comprise en kmin et kmax
c          k        kmin            kmax
c          1   0.0000000      0.10000000    
c          2  0.10000000      0.17782794    
c          3  0.17782794      0.31622776    
c          4  0.31622776      0.56234139    
c          5  0.56234139       1.0000000    
c          6   1.0000000       1.7782795    
c          7   1.7782795       3.1622777    
c          8   3.1622777       5.6234136    
c          9   5.6234136       10.000000    
c         10   10.000000       17.782795    
c         11   17.782795       31.622778    
c         12   31.622778       100000.00
c
c        mise a zero de occld_m
         occcld_m=0.
         DO i=1,klon
           DO j=1,klev
             DO k=1,12
               ex=10.**(0.25)
               kmin=0.
               kmax=1.e5
               if(k.ne.1)  kmin=0.1*ex**(k-2)
               if(k.ne.12) kmax=0.1*ex**(k-1)
               if(occcld(i,j).ge.kmin .and. occcld(i,j).lt.kmax)
     &         occcld_m(i,j,k)=1.
             ENDDO
           ENDDO
         ENDDO
       ENDIF  ! fin condition clouds => pas besoin de calculer des rayons

      RETURN
      END

