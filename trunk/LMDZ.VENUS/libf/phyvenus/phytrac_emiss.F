!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/phytrac.F,v 1.16 2006/03/24 15:06:23 lmdzadmin Exp $
!
c
c
      SUBROUTINE phytrac_emiss (timesimu,
     I                    debutphy,
     I                    lafin,
     I                    nqmax,
     I                    nlon,
     I                    nlev, 
     I                    pdtphys,
     I                    paprs,
     I                    xlat,xlon,
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
c Version developed for surface emission
c Maybe could be used just to compute the 'source' variable from physiq
c
c======================================================================
      USE infotrac_phy, ONLY: nqtot
      use dimphy
      USE geometry_mod, only: cell_area
      USE chemparam_mod,only:M_tr
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
      IMPLICIT none
#include "YOMCST.h"
#include "clesphys.h"
c======================================================================

c Arguments:

c   EN ENTREE:
c   ==========

      real timesimu   ! duree depuis debut simu (s)
      logical debutphy       ! le flag de l'initialisation de la physique
      logical lafin          ! le flag de la fin de la physique
      integer nqmax ! nombre de traceurs auxquels on applique la physique
      integer nlon  ! nombre de points horizontaux
      integer nlev  ! nombre de couches verticales
      real pdtphys  ! pas d'integration pour la physique (seconde)
      real paprs(nlon,nlev+1)  ! pression pour chaque inter-couche (en Pa)
      REAL xlat(nlon)       ! latitudes pour chaque point 
      REAL xlon(nlon)       ! longitudes pour chaque point 

c   EN ENTREE/SORTIE:
c   =================

      real tr_seri(nlon,nlev,nqmax) ! traceur  

cAA ----------------------------
cAA  VARIABLES LOCALES TRACEURS
cAA ----------------------------

c pour emission volcan
      real :: deltatr(klon,klev,nqtot)

      integer,parameter :: nblat=5,nblon=4,nbz=3
      integer,parameter :: Nemiss=0     ! duree emission (Ed)
      integer,save :: Nheight(nbz)      ! layer emission
      real,save :: so2_quantity         ! quantity so2 (kg)
      real,save :: lat_volcan(nblat),lon_volcan(nblon)
      real,save :: area_emiss(nblat,nblon)
      integer,save :: ig_volcan(nblat,nblon)

      INTEGER i, k, it
      integer ilat,ilon,iz
      real    deltalat,deltalon
c======================================================================

c EMISSION TRACEURS 

c---------
c debutphy
c---------      
      if (debutphy) then

        print*,"DEBUT PHYTRAC"
        print*,"PHYTRAC: EMISSION"

        ALLOCATE(M_tr(nqtot))
        M_tr(:)=64.                 ! SO2
	
C=========================================================================
c Caracteristiques des traceurs emis:
C=========================================================================

c nombre total de traceur
         if (nbz*nblat*nblon .gt. nqtot) then
            print*, nbz*nblat*nblon, nqtot
            write(*,*) "Attention, pas assez de traceurs"
            write(*,*) "le dernier sera bien le dernier"
         endif

c quantite en kg
         so2_quantity = 20.*10.**9.

c height (in layer index)
         Nheight(1) =  6  ! ~ 1 km
         Nheight(2) = 16  ! ~ 25 km
         Nheight(3) = 24  ! ~ 50 km

c localisation volcan
         lat_volcan(1) =  70.
         lat_volcan(2) =  35.
         lat_volcan(3) =   0.
         lat_volcan(4) = -35.
         lat_volcan(5) = -70.
         lon_volcan(1) = -125.
         lon_volcan(2) =  -35.
         lon_volcan(3) =   55.
         lon_volcan(4) =  145.
	 
	 ig_volcan(ilat,ilon)= 0
         if ((nbp_lon*nbp_lat)==1) then ! running a 1D simulation
           deltalat=180.
           deltalon=360.
         else
           deltalat = 180./(nbp_lat-1)
           deltalon = 360./nbp_lon
         endif

         do i=1,nlon
          do ilat=1,nblat
           do ilon=1,nblon
            if ((xlat(i).ge.lat_volcan(ilat))
     &     .and.((xlat(i)-deltalat).lt.lat_volcan(ilat))
     &     .and.(xlon(i).le.lon_volcan(ilon))
     &     .and.((xlon(i)+deltalon).gt.lon_volcan(ilon)) ) then
             ig_volcan(ilat,ilon)= i
             area_emiss(ilat,ilon) = cell_area(i)
	     print*,"Lat,lon=",ilat,ilon," OK"
            end if
           end do
          end do
         end do

c Reinit des traceurs si necessaire
         if (reinit_trac) then
	   tr_seri(:,:,:)=0.
c CAS N2 TRACEUR PASSIF POUR EVALUER PROFIL SOUS 7 KM
c J'ai mis Nemiss=0 !
           do i=1,klon
            do k=1,klev
              tr_seri(i,k,:)=max(min(0.035,
     &          0.035*(1.-log(paprs(i,k)/6.e6)/log(9.e6/6.e6))),0.)
            enddo
           enddo
c FIN CAS N2 PASSIF
	 endif
	 
C=========================================================================
C=========================================================================
      ENDIF  ! fin debutphy 
c-------------
c fin debutphy
c-------------

c======================================================================
c Emission d'un traceur pendant un certain temps
c necessite raz_date=1 dans run.def
c et reinit_trac=y
c======================================================================
       deltatr(:,:,:) = 0.

c source appliquee pendant Nemiss Ed
       if (timesimu .lt. 86400*Nemiss) then

c emet les traceurs qui sont presents sur la grille
        do ilat  = 1,nblat
        do ilon  = 1,nblon
         if (ig_volcan(ilat,ilon).ne.0) then
	 
          do iz = 1,nbz
           it=min( (iz-1)*nblat*nblon+(ilat-1)*nblon+ilon , nqtot )	    
           i=ig_volcan(ilat,ilon)
	   
c injection dans une seule cellule:
c source en kg/kg/s
c            deltatr(i,Nheight(iz),it) = so2_quantity/(86400.*Nemiss) ! kg/s
c     $ *RG/( area_emiss(ilat,ilon)
c     $      *(paprs(i,Nheight(iz))-paprs(i,Nheight(iz)+1)) )    ! /kg (masse cellule)
     
c            tr_seri(i,Nheight(iz),it) = tr_seri(i,Nheight(iz),it)
c     $      + deltatr(i,Nheight(iz),it)*pdtphys

c injection dans toute la colonne (a faire):
            do k=1,Nheight(iz)
             deltatr(i,k,it) = so2_quantity/(86400.*Nemiss) ! kg/s
     $  *RG/( area_emiss(ilat,ilon)
     $       *(paprs(i,1)-paprs(i,Nheight(iz)+1)) )    ! /kg (masse colonne)
     
             tr_seri(i,k,it) = tr_seri(i,k,it)+deltatr(i,k,it)*pdtphys
	    end do
	    
          end do
	  
         endif  ! ig_volcan!=0
        end do
        end do

       end if  ! duree emission
       
c======================================================================
c======================================================================

      RETURN
      END
