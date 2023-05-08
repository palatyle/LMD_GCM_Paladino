!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/radlwsw.F,v 1.2 2004/10/27 10:14:46 lmdzadmin Exp $
!
      SUBROUTINE radlwsw(dist, rmu0, fract, zzlev,
     .                  paprs, pplay,tsol, t)
c      
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
c Objet: interface entre le modele et les rayonnements
c Arguments:
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c zzlev----input-R- altitude a inter-couche (m)
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c tsol-----input-R- temperature du sol (en K)
c t--------input-R- temperature (K)
      
c MODIFS pour multimatrices ksi SPECIFIQUE VENUS
c   S. Lebonnois    20/12/2006
c   corrections     13/07/2007
c   New ksi matrix: possibility of different cloud model fct of lat   05/2014

c With extension NLTE (G. Gilli, 2014)

c Ksi matrices latitudinaly interpolated  (I. Garate-Lopez, 2016)

c======================================================================
      use dimphy
      USE geometry_mod, ONLY: latitude_deg
      USE phys_state_var_mod, only: heat,cool,radsol,
     .      topsw,toplw,solsw,sollw,sollwdown,lwnet,swnet
      use write_field_phy
      IMPLICIT none
#include "YOMCST.h"
#include "clesphys.h" 
#include "comcstVE.h"
#include "nlteparams.h"

!===========
! Arguments
!===========
      real rmu0(klon), fract(klon), dist

      REAL zzlev(klon,klev+1)
      real paprs(klon,klev+1), pplay(klon,klev)
      real tsol(klon)
      real t(klon,klev)

!===========
! Local
!===========
      INTEGER k, kk, i, j, band

      REAL   PPB(klev+1)
      REAL   PPA(klev)

      REAL   zfract, zrmu0,latdeg

      REAL   zheat(klev), zcool(klev)
      real   temp(klev),znivs(klev+1)
      REAL   ZFSNET(klev+1),ZFLNET(klev+1)
      REAL   ztopsw, ztoplw
      REAL   zsolsw, zsollw
cIM BEG
      REAL   zsollwdown
cIM END
      real,save,allocatable :: ksive(:,:,:,:) ! ksi matrixes in Vincent's file

      real    psi(0:klev+1,0:klev+1)
      real    deltapsi(0:klev+1,0:klev+1)
      real    pt0(0:klev+1)
      real    bplck(0:klev+1,nnuve)    ! Planck luminances in table layers
      real    y(0:klev,nnuve)          ! temporary variable for Planck
      real    zdblay(0:klev+1,nnuve)   ! temperature gradient of planck function
      integer mat0,lat,ips,isza,ips0,isza0
      real    factp,factz,ksi
c ------- for lat-interp ----------------
	integer mat0A, mat0B, latA, latB, kasua
	integer	ipsA, ipsB, iszaA, iszaB, ips0A, ips0B, isza0A, isza0B
	real	lat_deg, latA_deg, latB_deg
	real	factlat, k1, k2, k3, k4
c --------------------------------------
      logical firstcall
      data    firstcall/.true./
      save    firstcall
      
cERROR		! For checking if the file it's being read
c-------------------------------------------
c  Initialisations
c-----------------

      if (firstcall) then

c ---------- ksive --------------
        allocate(ksive(0:klev+1,0:klev+1,nnuve,nbmat))
        call load_ksi(ksive)

      endif ! firstcall
c-------------------------------------------

      DO k = 1, klev
       DO i = 1, klon
         heat(i,k)=0.
         cool(i,k)=0.
       ENDDO
      ENDDO


c+++++++ BOUCLE SUR LA GRILLE +++++++++++++++++++++++++
      DO j = 1, klon

c======================================================================
c  Initialisations
c ---------------

       DO k = 1, klev
        zheat(k) = 0.0
        zcool(k) = 0.0
       ENDDO
c       zheat(1:klev)=0.0	!Explicit loop (no change in performance)
c       zcool(1:klev)=0.0

       DO k = 1, klev+1
        ZFLNET(k) = 0.0
        ZFSNET(k) = 0.0
       ENDDO
c       ZFLNET(1:klev+1)=0.0
c       ZFSNET(1:klev+1)=0.0

       ztopsw = 0.0
       ztoplw = 0.0
       zsolsw = 0.0
       zsollw = 0.0
       zsollwdown = 0.0
      
       zfract = fract(j)
       zrmu0 = rmu0(j)

       DO k = 1, klev+1
         PPB(k) = paprs(j,k)/1.e5
       ENDDO
       DO k = 1,klev
         PPA(k) = pplay(j,k)/1.e5
       ENDDO

       pt0(0)  = tsol(j)
       DO k = 1, klev
         pt0(k) = t(j,k)
       ENDDO
       pt0(klev+1) = 0.
       
       DO k = 0,klev+1
       DO i = 0,klev+1
        psi(i,k) = 0.   ! positif quand nrj de i->k
        deltapsi(i,k) = 0. 
       ENDDO
       ENDDO

c======================================================================
c Getting psi and deltapsi
c ------------------------

c Planck function
c ---------------
      do band=1,nnuve
        do k=0,klev
c B(T,l) = al/(exp(bl/T)-1)
         y(k,band) = exp(bl(band)/pt0(k))-1.
         bplck(k,band) = al(band)/(y(k,band))
         zdblay(k,band)= al(band)*bl(band)*exp(bl(band)/pt0(k))/
     .                  ((pt0(k)*pt0(k))*(y(k,band)*y(k,band)))
        enddo
        bplck(klev+1,band) = 0.0
        zdblay(klev+1,band)= 0.0
      enddo

c finding the right matrixes
c --------------------------

	mat0  = 0
	mat0A = 0
	mat0B = 0

c    Latitude
c    --------
	lat  = 0
	latA = 0
	latB = 0

c	write(*,*) 'nlatve:', nlatve

	lat_deg = abs(latitude_deg(j))

c	if (nlatve.eq.1) then   ! clouds are taken as uniform
	if ((nlatve.eq.1).or.(lat_deg.le.25.)) then
	    lat  = 1
        elseif (lat_deg.le.50.) then
	    lat  = 1
	    latA = 1
	    latB = 2
	    latA_deg = 25.0
	    latB_deg = 55.0
        elseif (lat_deg.le.55.) then
	    lat  = 2
	    latA = 1
	    latB = 2
	    latA_deg = 25.0
	    latB_deg = 55.0
        elseif (lat_deg.le.60.) then
	    lat  = 2
	    latA = 2
	    latB = 3
	    latA_deg = 55.0
	    latB_deg = 65.0
        elseif (lat_deg.le.65.) then
	    lat  = 3
	    latA = 2
	    latB = 3
	    latA_deg = 55.0
	    latB_deg = 65.0
        elseif (lat_deg.le.70.) then
	    lat  = 3
	    latA = 3
	    latB = 4
	    latA_deg = 65.0
	    latB_deg = 75.0
        elseif (lat_deg.le.75.) then
	    lat  = 4
	    latA = 3
	    latB = 4
	    latA_deg = 65.0
	    latB_deg = 75.0
        elseif (lat_deg.le.80.) then
	    lat  = 4
	    latA = 4
	    latB = 5
	    latA_deg = 75.0
	    latB_deg = 85.0
        elseif (lat_deg.le.85.) then
            lat = 5
	    latA = 4
	    latB = 5
	    latA_deg = 75.0
	    latB_deg = 85.0
        else
            lat = 5
	endif

c        write(*,*) 'Lat',lat,'LatA',latA,'LatB',latB

	factlat = 0
	if (latA.gt.0) then
	  factlat = (lat_deg - latA_deg) / (latB_deg - latA_deg)
	endif

c       write (*,*) 'Factor de correccion:', factlat


c    Pressure at Surface
c    -------------------    

       ips0=0
       ips0A=0
       ips0B=0
       if (nbpsve(lat).gt.1) then            ! Interpolation on ps
       do ips=1,nbpsve(lat)-1
         if (  (psurfve(ips,lat).ge.paprs(j,1))
     .  .and.(psurfve(ips+1,lat).lt.paprs(j,1)) ) then
              ips0  = ips
c             print*,'ig=',j,'  ips0=',ips
              factp = (paprs(j,1)         -psurfve(ips0,lat))
     .               /(psurfve(ips0+1,lat)-psurfve(ips0,lat))
              exit
         endif
       enddo
       else            ! Only one ps, no interpolation
        ips0=1
       endif

       if (latA.eq.lat) then
         ips0A=ips0
       else
         if (latA.gt.0) then
           if (nbpsve(latA).gt.1) then
	     do ipsA=1,nbpsve(latA)-1
	       if (  (psurfve(ipsA,latA).ge.paprs(j,1))
     .        .and.(psurfve(ipsA+1,latA).lt.paprs(j,1)) ) then
                ips0A = ipsA
                exit
	       endif
	     enddo
	   else            ! Only one ps, no interpolation
	     ips0A=1
	   endif     ! nbpsve(latA).gt.1
         endif       ! latA.gt.0   (if latA=0 ips0A is not used, so it doesn't matter)
       endif         ! latA.eq.lat

       if (latB.eq.lat) then
         ips0B=ips0
       else
         if (latB.gt.0) then
           if (nbpsve(latB).gt.1) then
	     do ipsB=1,nbpsve(latB)-1
	       if (  (psurfve(ipsB,latB).ge.paprs(j,1))
     .        .and.(psurfve(ipsB+1,latB).lt.paprs(j,1)) ) then
                 ips0B = ipsB
                 exit
	       endif
	     enddo
	   else
	     ips0B=1
	   endif     ! nbpsve(latB).gt.1
         endif       ! latB.gt.0   (if latB=0 ips0B is not used, so it doesn't matter)
       endif         ! latB.eq.lat


c    Solar Zenith Angle
c    ------------------

       isza0=0
       isza0A=0
       isza0B=0
       if (nbszave(lat).gt.1) then
        do isza=1,nbszave(lat)-1
         if (  (szave(isza,lat).ge.zrmu0)
     .  .and.(szave(isza+1,lat).lt.zrmu0) ) then
              isza0  = isza
c             print*,'ig=',j,'  isza0=',isza
              factz = (zrmu0             -szave(isza0,lat))
     .               /(szave(isza0+1,lat)-szave(isza0,lat))
              exit
         endif
        enddo
       else            ! Only one sza, no interpolation
        isza0=-99
       endif


       if (latA.eq.lat) then
         isza0A=isza0
       else
         if (latA.gt.0) then
           if (nbszave(latA).gt.1) then
	     do iszaA=1,nbszave(latA)-1
	       if (  (szave(iszaA,latA).ge.zrmu0)
     .        .and.(szave(iszaA+1,latA).lt.zrmu0) ) then
                 isza0A = iszaA
                 exit
	       endif
	     enddo
	   else            ! Only one sza, no interpolation
	     isza0A=-99
	   endif     ! nbszave(latA).gt.1
         endif       ! latA.gt.0   (if latA=0 isza0A is not used, so it doesn't matter)
       endif         ! latA.eq.lat

       if (latB.eq.lat) then
         isza0B=isza0
       else
         if (latB.gt.0) then
           if (nbszave(latB).gt.1) then
             ! init to avoid outside values (near midnight so similar compo...)
             isza0B = nbszave(latB) 
	     do iszaB=1,nbszave(latB)-1
	       if (  (szave(iszaB,latB).ge.zrmu0)
     .        .and.(szave(iszaB+1,latB).lt.zrmu0) ) then
                 isza0B = iszaB
                 exit
	       endif
	     enddo
	   else            ! Only one sza, no interpolation
	     isza0B=-99
	   endif     ! nbszave(latB).gt.1
         endif       ! latB.gt.0   (if latB=0 isza0B is not used, so it doesn't matter)
       endif         ! latB.eq.lat

c        write(*,*) 'nbszave', nbszave(lat),'nbpsve(lat)',nbpsve(lat)

       
c -------- Probleme aux bords
c surf press lower than the lowest surf pres in matrices
       if ((ips0.eq.0).and.(psurfve(nbpsve(lat),lat).gt.paprs(j,1)))
     . then
              ips0  = nbpsve(lat)-1
              print*,'Extrapolation! ig=',j,'  ips0=',ips0
              factp = (paprs(j,1)         -psurfve(ips0,lat))
     .               /(psurfve(ips0+1,lat)-psurfve(ips0,lat))
       endif
c surf press higher than the highest surf pres in matrices
       if ((ips0.eq.0).and.(psurfve(1,lat).le.paprs(j,1))) then
              ips0  = 1
              print*,'Extrapolation! ig=',j,'  ips0=',ips0
              factp = (paprs(j,1)         -psurfve(ips0,lat))
     .               /(psurfve(ips0+1,lat)-psurfve(ips0,lat))
       endif

c this has to be done for ips0A and ips0B also...
       if (latA.eq.lat) then
	 ips0A = ips0
       else
         if (latA.gt.0) then
          if ((ips0A.eq.0).and.
     .        (psurfve(nbpsve(latA),latA).gt.paprs(j,1))) then
              ips0A  = nbpsve(latA)-1
          endif
          if ((ips0A.eq.0).and.(psurfve(1,latA).le.paprs(j,1))) then
              ips0A  = 1
          endif
         endif
       endif
       if (latB.eq.lat) then
	 ips0B = ips0
       else
         if (latB.gt.0) then
          if ((ips0B.eq.0).and.
     .        (psurfve(nbpsve(latB),latB).gt.paprs(j,1))) then
              ips0B  = nbpsve(latB)-1
          endif
          if ((ips0B.eq.0).and.(psurfve(1,latB).le.paprs(j,1))) then
              ips0B  = 1
          endif
         endif
       endif       

c ---------

       if ((ips0.eq.0).or.(isza0.eq.0)) then
         write(*,*) 'Finding the right matrix in radlwsw'
         print*,'Interpolation problem, grid point ig=',j
         print*,'psurf = ',paprs(j,1),' mu0 = ',zrmu0
         stop
       endif

       if (isza0.eq.-99) then
         mat0  = indexve(lat) +ips0
         if (latA.gt.0) then
           mat0A = indexve(latA)+ips0A
           mat0B = indexve(latB)+ips0B
         endif
       else
         mat0  = indexve(lat) +(isza0 -1)*nbpsve(lat) +ips0
         if (latA.gt.0) then
           mat0A = indexve(latA)+(isza0A-1)*nbpsve(latA)+ips0A
           mat0B = indexve(latB)+(isza0B-1)*nbpsve(latB)+ips0B
         endif
       endif

c        write(*,*) 'Second revision> Lat',lat,'LatA',latA,'LatB',latB
   
c interpolation of ksi and computation of psi,deltapsi
c ----------------------------------------------------

       if (isza0.eq.-99) then
         if (latA.gt.0) then		! Not being in the two extremal bins

           do band=1,nnuve
             do k=0,klev+1
               do i=k+1,klev+1
	         k1 = ksive(i,k,band,mat0A)*(1-factlat)
     .              + ksive(i,k,band,mat0B)*factlat
	         k2 = ksive(i,k,band,mat0A+1)*(1-factlat)
     .              + ksive(i,k,band,mat0B+1)*factlat
                 ksi = k1*(1-factp) + k2*factp
                 psi(i,k) = psi(i,k) + 
     .                      RPI*ksi*(bplck(i,band)-bplck(k,band))
c ONLY NEEDED IF IMPLICIT CHOSEN IN LW_VENUS_VE (not the case right now)
c                 deltapsi(i,k) = deltapsi(i,k) + RPI*ksi*zdblay(i,band)
c                 deltapsi(k,i) = deltapsi(k,i) + RPI*ksi*zdblay(k,band)

	         kasua=1
               enddo
             enddo
           enddo
             do k=0,klev+1
               do i=k+1,klev+1
                 psi(k,i) = -psi(i,k) 
               enddo
             enddo

         else 		! latA=0 --> extremal bins

           do band=1,nnuve
             do k=0,klev+1
               do i=k+1,klev+1
                 ksi = ksive(i,k,band,mat0)*(1-factp)
     .               + ksive(i,k,band,mat0+1)*factp
                 psi(i,k) = psi(i,k) + 
     .                      RPI*ksi*(bplck(i,band)-bplck(k,band))
c ONLY NEEDED IF IMPLICIT CHOSEN IN LW_VENUS_VE (not the case right now)
c                 deltapsi(i,k) = deltapsi(i,k) + RPI*ksi*zdblay(i,band)
c                 deltapsi(k,i) = deltapsi(k,i) + RPI*ksi*zdblay(k,band)

	         kasua=2
               enddo
             enddo
           enddo
             do k=0,klev+1
               do i=k+1,klev+1
                 psi(k,i) = -psi(i,k) 
               enddo
             enddo

         endif 		! latA.gt.0

       else		! isza0=!-99

         if (latA.gt.0) then		! Not being in the two extremal bins

           do band=1,nnuve
             do k=0,klev+1
               do i=k+1,klev+1
	         k1 = ksive(i,k,band,mat0A)*(1-factlat)
     .              + ksive(i,k,band,mat0B)*factlat
	         k2 = ksive(i,k,band,mat0A+1)*(1-factlat)
     .              + ksive(i,k,band,mat0B+1)*factlat
	         k3 = ksive(i,k,band,mat0A+nbpsve(latA))*(1-factlat)
     .              + ksive(i,k,band,mat0B+nbpsve(latB))*factlat
	         k4 = ksive(i,k,band,mat0A+nbpsve(latA)+1)*(1-factlat)
     .              + ksive(i,k,band,mat0B+nbpsve(latB)+1)*factlat
	         ksi = ( k1*(1-factp) + k2*factp )*(1-factz)
     .              + ( k3*(1-factp) + k4*factp )*factz
                 psi(i,k) = psi(i,k) + 
     .                      RPI*ksi*(bplck(i,band)-bplck(k,band))
c ONLY NEEDED IF IMPLICIT CHOSEN IN LW_VENUS_VE (not the case right now)
c                 deltapsi(i,k) = deltapsi(i,k) + RPI*ksi*zdblay(i,band)
c                 deltapsi(k,i) = deltapsi(k,i) + RPI*ksi*zdblay(k,band)

	        kasua=3
               enddo
             enddo
           enddo
             do k=0,klev+1
               do i=k+1,klev+1
                 psi(k,i) = -psi(i,k) 
               enddo
             enddo

         else 		! latA=0 --> extremal bins

           do band=1,nnuve
             do k=0,klev+1
               do i=k+1,klev+1
                 ksi = ksive(i,k,band,mat0)*(1-factp)*(1-factz)
     .               + ksive(i,k,band,mat0+1)*factp  *(1-factz)
     .               + ksive(i,k,band,mat0+nbpsve(lat))*(1-factp)*factz
     .               + ksive(i,k,band,mat0+nbpsve(lat)+1)*factp  *factz
                 psi(i,k) = psi(i,k) + 
     .                      RPI*ksi*(bplck(i,band)-bplck(k,band))
c ONLY NEEDED IF IMPLICIT CHOSEN IN LW_VENUS_VE (not the case right now)
c                 deltapsi(i,k) = deltapsi(i,k) + RPI*ksi*zdblay(i,band)
c                 deltapsi(k,i) = deltapsi(k,i) + RPI*ksi*zdblay(k,band)

	         kasua=4
               enddo
             enddo
           enddo
             do k=0,klev+1
               do i=k+1,klev+1
                 psi(k,i) = -psi(i,k) 
               enddo
             enddo

         endif 		! latA.gt.0
       endif		! isza0.eq.-99

c	write(*,*) 'Kasua:', kasua

c======================================================================
c LW call
c---------
      temp(1:klev)=t(j,1:klev)
      CALL LW_venus_ve(
     .        PPB,temp,psi,deltapsi,
     .        zcool,
     .        ztoplw,zsollw,
     .        zsollwdown,ZFLNET)
c---------
c SW call
c---------
      znivs=zzlev(j,:)
      latdeg=abs(latitude_deg(j))

c      CALL SW_venus_ve_1Dglobave(zrmu0, zfract,   ! pour moy globale
c      CALL SW_venus_ve(zrmu0, zfract,
c     S        PPB,temp,znivs,
c     S        zheat, 
c     S        ztopsw,zsolsw,ZFSNET)

c      CALL SW_venus_cl_1Dglobave(zrmu0,zfract,   ! pour moy globale
c      CALL SW_venus_cl(zrmu0,zfract,
c      CALL SW_venus_dc_1Dglobave(zrmu0,zfract,   ! pour moy globale
c      CALL SW_venus_dc(zrmu0,zfract,
c      CALL SW_venus_rh_1Dglobave(zrmu0,zfract,   ! pour moy globale
c     S        PPB,temp, 
      CALL SW_venus_rh(zrmu0,zfract,latdeg,
     S        PPA,PPB,temp, 
     S        zheat, 
     S        ztopsw,zsolsw,ZFSNET)
c======================================================================
         radsol(j) = zsolsw - zsollw  ! + vers bas
         topsw(j) = ztopsw            ! + vers bas
         toplw(j) = ztoplw            ! + vers haut
         solsw(j) = zsolsw            ! + vers bas
         sollw(j) = -zsollw           ! + vers bas
         sollwdown(j) = zsollwdown    ! + vers bas

         DO k = 1, klev+1
         lwnet  (j,k)   = ZFLNET(k)
         swnet  (j,k)   = ZFSNET(k)
         ENDDO

c
C heat/cool with upper atmosphere 
C
      IF(callnlte) THEN
         DO k = 1,nlaylte
           heat(j,k) = zheat(k)
           cool(j,k) = zcool(k)
         ENDDO
c     Zero tendencies for any remaining layers between nlaylte and klev
       if (klev.gt.nlaylte) then
         do k = nlaylte+1,  klev
           heat(j,k) = 0.
           cool(j,k) = 0.
         enddo
       endif    
      ELSE
         DO k = 1, klev
           heat(j,k) = zheat(k)
           cool(j,k) = zcool(k)
         ENDDO
      ENDIF ! callnlte
      ENDDO ! of DO j = 1, klon
c+++++++ FIN BOUCLE SUR LA GRILLE +++++++++++++++++++++++++

! for tests: write output fields...
!      call writefield_phy('radlwsw_heat',heat,klev)
!      call writefield_phy('radlwsw_cool',cool,klev)
!      call writefield_phy('radlwsw_radsol',radsol,1)
!      call writefield_phy('radlwsw_topsw',topsw,1)
!      call writefield_phy('radlwsw_toplw',toplw,1)
!      call writefield_phy('radlwsw_solsw',solsw,1)
!      call writefield_phy('radlwsw_sollw',sollw,1)
!      call writefield_phy('radlwsw_sollwdown',sollwdown,1)
!      call writefield_phy('radlwsw_swnet',swnet,klev+1)
!      call writefield_phy('radlwsw_lwnet',lwnet,klev+1)

c tests

c     j = klon/2
c     j = 1
c     print*,'mu0=',rmu0(j)
c     print*,'   net flux vis   HEAT(K/Eday)'
c     do k=1,klev
c     print*,k,ZFSNET(k),heat(j,k)*86400.
c     enddo
c     print*,'   net flux IR    COOL(K/Eday)'
c     do k=1,klev
c     print*,k,ZFLNET(k),cool(j,k)*86400.
c     enddo

      firstcall = .false.
      RETURN
      END
