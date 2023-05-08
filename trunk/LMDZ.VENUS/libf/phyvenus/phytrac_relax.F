!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/phytrac.F,v 1.16 2006/03/24 15:06:23 lmdzadmin Exp $
!
c
c
      SUBROUTINE phytrac_relax (debutphy,lafin,
     I                    nqmax,
     I                    nlon,
     I                    nlev, 
     I                    pdtphys,
     I                    pplay,
     O                    tr_seri)

c======================================================================
c Auteur(s) FH
c Objet: Moniteur general des tendances traceurs
c
cAA Remarques en vrac:
cAA--------------------
cAA 1/ le call phytrac se fait avec nqmax 
c
c SL: Janvier 2014
c Version developed by E. Marcq for pseudo-chemistry relaxation
c See Marcq&Lebonnois 2013.
c
c======================================================================
      USE infotrac_phy, ONLY: nqtot, tname
      use dimphy
      USE chemparam_mod,only:M_tr
      IMPLICIT none
#include "YOMCST.h"
#include "clesphys.h"
c======================================================================

c Arguments:

c   EN ENTREE:
c   ==========

      logical debutphy       ! le flag de l'initialisation de la physique
      logical lafin          ! le flag de la fin de la physique
      integer nqmax ! nombre de traceurs auxquels on applique la physique
      integer nlon  ! nombre de points horizontaux
      integer nlev  ! nombre de couches verticales
      real pdtphys  ! pas d'integration pour la physique (seconde)
      real pplay(nlon,nlev)  ! pression pour le mileu de chaque couche (en Pa)

c   EN ENTREE/SORTIE:
c   =================

      real tr_seri(nlon,nlev,nqmax) ! traceur  

cAA ----------------------------
cAA  VARIABLES LOCALES TRACEURS
cAA ----------------------------

C les traceurs

c===================
c it--------indice de traceur
c k,i---------indices long, vert
c===================
c Variables deja declarees dont on a besoin pour traceurs   
c   k,i,it,tr_seri(nlon,nlev,nqmax),pplay(nlon,nlev), 
      integer nqCO_OCS
c      real pzero,gamma
c      parameter (pzero=85000.)
c      parameter (gamma=5000.)
      REAL alpha
      real deltatr(nlon,nlev,nqtot) ! ecart au profil de ref zprof
      real,save,allocatable :: zprof(:,:)
      real,save,allocatable :: tau(:,:) ! temps de relaxation vers le profil (s)
c======================================================================

      INTEGER i, k, it

c Variables liees a l'ecriture de la bande histoire physique

c Variables locales pour effectuer les appels en serie
c----------------------------------------------------

      REAL d_tr(nlon,nlev) ! tendances de traceurs 

      character*20 modname
      character*80 abort_message

c======================================================================

      modname = 'phytrac_relax'
c TRACEURS TYPE CO ET OCS
      nqCO_OCS   = 6

c !!!!!! ATTENTION: A REVOIR / A VERIFIER
c   les traceurs sont en mass mixing ratio dans la dyn
c   et convertis en frac mol avec M_tr dans les sorties...
c   Mettre le profil de rappel en mass mixing ratio !!

      print*,"METTRE A JOUR phytrac_relax"
      stop

c---------
c debutphy
c---------
      if (debutphy) then
         print*,"DEBUT PHYTRAC"
         print*,"PHYTRAC: RELAXATION"
         allocate(zprof(nlev,nqtot),tau(nlev,nqtot))

         ALLOCATE(M_tr(nqtot))
      
c=============================================================
c=============================================================
c=============================================================
c   Initialisation des traceurs
c=============================================================
c=============================================================
c=============================================================

C=========================================================================
C=========================================================================

c II) Declaration d'un profil vertical de traceur OK
c 
c zprof   = profil de rappel 
c
c 1 -> CO ; 2 -> OCS 
c def des profils en log(a) = a * log(P) + b par morceaux, cf. pollack et al
c tr_seri en ppm
c (initialisation seulement si ceux-ci sont nuls)

c ICI, ON UTILISE 3 CONSTANTES DE TEMPS DIFFERENTES POUR CHAQUE,
c DONC TRACEURS 1 A 3 POUR CO ET 4 A 6 POUR OCS
C=========================================================================


c Constantes de rappel:

       print*,"INIT TAU"
       do k=1,nlev
         tau(k,1)=1.e6
         tau(k,2)=1.e7
         tau(k,3)=1.e8
         tau(k,4)=1.e6
         tau(k,5)=1.e7
         tau(k,6)=1.e8
       enddo

c CO

      do it=1,3 
       print*,"INIT ZPROF ",tname(it)
       M_tr(it)=28.       ! CO
       do k=1,nlev
         zprof(k,it)=0.
c pour l'instant, tau fixe, mais possibilite de le faire varier avec z
        if (pplay(nlon/2,k) >= 4.8e6) then
           zprof(k,it)=14.
        endif
        if ((pplay(nlon/2,k)<=4.8e6).and.(pplay(nlon/2,k)>=1.9e6)) then
           alpha=(log(pplay(nlon/2,k))-log(1.9e6))/
     .     (log(4.8e6)-log(1.9e6))
           zprof(k,it)=20.*(14./20.)**alpha
        endif
        if ((pplay(nlon/2,k)<=1.9e6).and.(pplay(nlon/2,k)>=1.5e5)) then
           alpha=(log(pplay(nlon/2,k))-log(1.5e5))/
     .     (log(1.9e6)-log(1.5e5))
           zprof(k,it)=39.*(20./39.)**alpha
        endif
        if ((pplay(nlon/2,k)<=1.5e5).and.(pplay(nlon/2,k)>=1.1e4)) then
           alpha=(log(pplay(nlon/2,k))-log(1.1e4))/
     .     (log(2.73e5)-log(1.1e4))
           zprof(k,it)=50.*(39./50.)**alpha
        endif
        if ((pplay(nlon/2,k)<=1.1e4).and.(pplay(nlon/2,k)>=1.3e3)) then
           alpha=(log(pplay(nlon/2,k))-log(1.3e3))/
     .     (log(1.1e4)-log(1.3e3))
           zprof(k,it)=2.*(50./2.)**alpha
        endif
        if ((pplay(nlon/2,k)<=1.3e3).and.(pplay(nlon/2,k)>=2.4)) then
           alpha=(log(pplay(nlon/2,k))-log(2.4))/
     .     (log(1.3e3)-log(2.4))
           zprof(k,it)=1000.*(2./1000.)**alpha
        endif
        if (pplay(nlon/2,k) <= 2.4) then
           zprof(k,it)=1000.
        endif
       enddo 
       print*,zprof(:,it)
 
c OCS
       print*,"INIT ZPROF ",tname(it+3)
       M_tr(it+3)=60.       ! OCS
       do k=1,nlev
         zprof(k,it+3)=0.
         if (pplay(nlon/2,k) >= 4.8e6) then
           zprof(k,it+3)=30.
         endif
         if ((pplay(nlon/2,k)<=4.8e6).and.(pplay(nlon/2,k)>=9.4e5))
     *   then
           alpha=(log(pplay(nlon/2,k))-log(9.4e5))/
     *     (log(4.8e6)-log(9.4e5))
           zprof(k,it+3)=20.*(30/20.)**alpha
         endif
         if ((pplay(nlon/2,k)<=9.4e5).and.(pplay(nlon/2,k)>=4.724e5))
     *   then
           alpha=(log(pplay(nlon/2,k))-log(4.724e5))/
     *     (log(9.4e5)-log(4.724e5))
           zprof(k,it+3)=0.5*(20/0.5)**alpha
         endif
         if ((pplay(nlon/2,k)<=4.724e5).and.(pplay(nlon/2,k)>=1.1e4))
     *   then
           alpha=(log(pplay(nlon/2,k))-log(1.1e4))/
     *     (log(4.724e5)-log(1.1e4))
           zprof(k,it+3)=0.005*(0.5/0.005)**alpha
         endif
         if (pplay(nlon/2,k)<=1.1e4) then
           zprof(k,it+3)=0.
         endif
       end do
       print*,zprof(:,it+3)
      enddo 

c Initialisation du traceur s'il est nul:
       do it=1,nqCO_OCS
        if ((tr_seri(nlon/2,1,it).eq.0.).and.
     .      (tr_seri(nlon/2,nlev/2,it).eq.0.).and.
     .      (tr_seri(nlon/2,nlev,it).eq.0.)) then
         print*,"INITIALISATION DE ",tname(it)
         do k=1,nlev
           do i=1,nlon
             tr_seri(i,k,it) = zprof(k,it)
           enddo
         enddo
        endif
       enddo

C=========================================================================
C=========================================================================

c-------------
c fin debutphy
c-------------
      ENDIF  ! fin debutphy 

c======================================================================
c Rappel vers un profil
c======================================================================
         do it=1,nqCO_OCS
           do k=1,nlev
             do i=1,nlon
c VERIF
           if (tr_seri(i,k,it).lt.0) then
             print*,"Traceur negatif AVANT rappel:",i,k,it
             stop
           endif
c FIN VERIF

           deltatr(i,k,it) = (-tr_seri(i,k,it)+zprof(k,it))/tau(k,it)
           tr_seri(i,k,it) =  tr_seri(i,k,it) + deltatr(i,k,it)*pdtphys

c VERIF
           if (tr_seri(i,k,it).lt.0) then
             print*,"APRES rappel:",i,k,it,
     .  deltatr(i,k,it),zprof(k,it),tr_seri(i,k,it),pdtphys/tau(k,it)
             stop
           endif
c FIN VERIF
             enddo
           enddo
         enddo

c======================================================================
c======================================================================


      RETURN
      END
      
      
