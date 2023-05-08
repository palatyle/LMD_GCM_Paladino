      subroutine ballon (iam,dtphys,rjour,rsec,plat,plon,
     i                   temp, p, u, v, geop)

      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      implicit none

c======================================================================
c Auteur: S. Lebonnois (LMD/CNRS) date: 20091201
c Object: Compute balloon trajectories.
C No outputs, every quantities are written on the iam+ Files. 
c 
c Called by physiq.F if flag ballons activated:
c
c      integer ballons
c (...)
c      ballons  = 1         ! (in initialisations)
c (...)
C  OUVERTURE DES FICHIERS FORMATTES CONTENANT LES POSITIONS ET VITESSES
C  DES BALLONS
c      if (ballons.eq.1) then
c      open(30,file='ballons-lat.out',form='formatted')
c      open(31,file='ballons-lon.out',form='formatted')
c      open(32,file='ballons-u.out',form='formatted')
c      open(33,file='ballons-v.out',form='formatted')
c      open(34,file='ballons-alt.out',form='formatted')
c      write(*,*)'Ouverture des ballons*.out'
c      endif !ballons
c (...)
C        CALL ballon(30,pdtphys,rjourvrai,gmtime,rlatd,rlond,
CC    C               t,pplay,u,v,pphi)   ! alt above surface (smoothed for GCM)
C     C               t,pplay,u,v,zphi)   ! alt above planet average radius
c (...)
C  FERMETURE DES FICHIERS FORMATTES CONTENANT LES POSITIONS ET VITESSES
C  DES BALLONS
c      if (ballons.eq.1) then
c        write(*,*)'Fermeture des ballons*.out'
c        close(30)                                     
c        close(31)                                     
c        close(32)                                     
c        close(33)                                     
c        close(34)                                     
c      endif !ballons
c
C======================================================================
c Explicit Arguments:
c ==================
c iam-----input-I-File number where latitudes are written
c                 It is a formatted file that has been opened
c                 in physiq.F
c   other files: iam+1=longitudes
c                iam+2=zonal speeds
c                iam+3=meridional speeds
c                iam+4=altitudes
c dtphys--input-R-pas de temps physique
c rjour---input-R-Jour compte depuis le debut de la simu (run.def)
c rsec----input-R-Seconde de la journee
c plat ---input-R-Latitude en degres
c plon ---input-R-Longitude en degres
c temp----input-R-Temperature (K) at model levels
c p-------input-R-Pressure (Pa) at model levels
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c geop----input-R-Geopotential !! above surface OR average radius
c
c
c Implicit Arguments:
c ===================
c
c iim--common-I: Number of longitude intervals
c jjm--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                iim*(jjm-1)+2 for instance
c klev-common-I: Number of vertical layers
c RPI,RKBOL--common-R: Pi, KBoltzman
c RDAY,RA,RG-common-R: day length in s, planet radius, gravity
c======================================================================
c Local Variables:
c ================
c
c nb    ---I: number of balloons (parameter)
c phib  ---R: Latitude of balloon in radians
c lamb  ---R: Longitude of balloon in radians
c lognb ---R: log(density) of balloon
c ub    ---R: zonal speed of balloon
c vb    ---R: meridional speed of balloon
c altb  ---R: altitude of balloon
c zlat  ---R: Latitude in radians
c zlon  ---R: Longitude in radians
c logn  ---R: log(density)
c alt   ---R: altitude !! above surface OR average radius
c ull   ---R: zonal wind for one balloon on the lognb surface
c vll   ---R: meridional wind for one balloon on the lognb surface
c aal   ---R: altitude for one balloon on the lognb surface
c======================================================================

#include "YOMCST.h"
c
c ARGUMENTS
c
      INTEGER iam
      REAL dtphys,rjour,rsec,plat(klon),plon(klon)
      REAL temp(klon,klev),p(klon,klev)
      REAL u(klon,klev),v(klon,klev),geop(klon,klev)
c
c Variables locales:
c
      INTEGER i,j,k,l,nb,n
      parameter (nb=20)  !! Adjust the format on line 100 !!
      INTEGER jj,ii,ll

      REAL,SAVE,ALLOCATABLE :: zlon(:),zlat(:)

      REAL time
      REAL logn(klon,klev),ull(klon),vll(klon)
      REAL alt(klon,klev),aal(klon)
      real ub(nb),vb(nb),phib(nb),lamb(nb),lognb(nb),altb(nb)
      save phib,lamb,lognb

      REAL factalt

c RungeKutta order - If not RK, Nrk=1
      integer Nrk,irk
      parameter (Nrk=1)
      real    dtrk

      logical first
      save first
      data first/.true./

      time = rjour*RDAY+rsec
      logn(:,:) = log10(p(:,:)/(RKBOL*temp(:,:)))
      alt(:,:)  = geop(:,:)/RG

c---------------------------------------------
C INITIALIZATIONS
c---------------------------------------------
      if (first) then

      print*,"BALLOONS ACTIVATED"

      allocate(zlon(nbp_lon+1))
      allocate(zlat(nbp_lat))
      
C Latitudes:
      zlat(1)=plat(1)*RPI/180.
      do j = 2,nbp_lat-1
         k=(j-2)*nbp_lon+2
         zlat(j)=plat(k)*RPI/180.
      enddo
      zlat(nbp_lat)=plat(klon)*RPI/180.

C Longitudes:
      do i = 1,nbp_lon
         k=i+1
         zlon(i)=plon(k)*RPI/180.
      enddo
      zlon(nbp_lon+1)=zlon(1)+2.*RPI

c verif init     lat de 90 à -90, lon de -180 à 180
c     print*,"Latitudes:",zlat*180./RPI
c     print*,"Longitudes:",zlon*180./RPI
c     stop

c initial positions of balloons (in degrees for lat/lon)
      do j=1,5
      do i=1,4
        k=(j-1)*4+i
      phib(k)= (j-1)*20.*RPI/180.
      lamb(k)= (i-3)*90.*RPI/180.   ! de -180 à 90
c     lognb(k)= log10(5.e4/(RKBOL*300.)) ! ~55km in VIRA model
      lognb(k)= log10(5.e5/(RKBOL*300.)) ! 5 bars (for Blamont, mai2015)
      enddo
      enddo
      print*,"Balloon density (m^-3)=",10.**(lognb(1))

c     print*,"log(density) profile:"
c     do l=1,klev
c        print*,logn(klon/2,l)
c     enddo
c     stop !verif init

      first=.false.
      endif ! first
c---------------------------------------------

c-------------------------------------------------
c loop over the balloons
c-------------------------------------------------
      do n=1,nb

c Interpolation in altitudes
c-------------------------------------------------
        do k=1,klon
         ll=1 ! en bas
         do l=2,klev
          if (lognb(n).lt.logn(k,l)) ll=l
         enddo
         factalt= (lognb(n)-logn(k,ll))/(logn(k,ll+1)-logn(k,ll))
         ull(k) =   u(k,ll+1)*factalt +   u(k,ll)*(1-factalt)
         vll(k) =   v(k,ll+1)*factalt +   v(k,ll)*(1-factalt)
         aal(k) = alt(k,ll+1)*factalt + alt(k,ll)*(1-factalt)
        enddo

c Interpolation in latitudes and longitudes
c-------------------------------------------
        call wind_interp(ull,vll,aal,zlat,zlon,
     .                   phib(n),lamb(n),ub(n),vb(n),altb(n))
        
      enddo ! over balloons
c-------------------------------------------------

c-------------------------------------------------
c Output of positions and speed at time
c-------------------------------------------------

c Venus regardee a l'envers: lon et lat inversees
      write(iam,  100) time, (-1)*phib*180./RPI
      write(iam+1,100) time, (-1)*lamb*180./RPI
      write(iam+2,100) time, ub
c Venus regardee a l'envers: v inversee
      write(iam+3,100) time, (-1)*vb
      write(iam+4,100) time, altb
c     stop !verif init

c !!!!!!!!!!!!!!!! nb !!!!!!!!!!!!!!!!!
100   format(E14.7,20(1x,E12.5))

c-------------------------------------------------
c Implementation: positions at time+dt 
c RK order Nrk
c-------------------------------------------------

      dtrk = dtphys/Nrk
      time=time+dtrk

      do n=1,nb
        call pos_implem(phib(n),lamb(n),ub(n),vb(n),dtrk)
      enddo

      if (Nrk.gt.1) then
       do irk=2,Nrk
        do n=1,nb
          time=time+dtrk
          call wind_interp(ull,vll,aal,zlat,zlon,
     .                   phib(n),lamb(n),ub(n),vb(n),altb(n))
          call pos_implem(phib(n),lamb(n),ub(n),vb(n),dtrk)
        enddo
       enddo
      endif

      end

c======================================================================
c======================================================================
c======================================================================

      subroutine wind_interp(map_u,map_v,map_a,latit,longit,
     .                       phi,lam,ubal,vbal,abal)

      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      implicit none

c======================================================================
c Auteur: S. Lebonnois (LMD/CNRS) date: 20091201
c Object: interpolate balloon speed from its position.
C======================================================================
c Explicit Arguments:
c ==================
c map_u ---R: zonal wind on the lognb surface
c map_v ---R: meridional wind on the lognb surface
c map_a ---R: altitude on the lognb surface
c latit ---R: Latitude in radians
c longit---R: Longitude in radians
c phi   ---R: Latitude of balloon in radians
c lam   ---R: Longitude of balloon in radians
c ubal  ---R: zonal speed of balloon
c vbal  ---R: meridional speed of balloon
c abal  ---R: altitude of balloon
c======================================================================
c Local Variables:
c ================
c
c ujj   ---R: zonal wind interpolated in latitude
c vjj   ---R: meridional wind interpolated in latitude
c ajj   ---R: altitude interpolated in latitude
c======================================================================

#include "YOMCST.h"
c
c ARGUMENTS
c
      real map_u(klon),map_v(klon),map_a(klon)
      real latit(nbp_lat),longit(nbp_lon)
      real phi,lam,ubal,vbal,abal
c
c Variables locales:
c
      INTEGER i,j,k
      INTEGER jj,ii
      REAL    ujj(nbp_lon+1),vjj(nbp_lon+1),ajj(nbp_lon+1)
      REAL    factlat,factlon

c Interpolation in latitudes
c-------------------------------------------------
        jj=1  ! POLE NORD 
        do j=2,nbp_lat-1
          if (phi.lt.latit(j)) jj=j
        enddo
        factlat  = (phi-latit(jj))/(latit(jj+1)-latit(jj))

c pole nord
        if (jj.eq.1) then
         do i=1,nbp_lon
          ujj(i) = map_u(i+1)*factlat + map_u(1)*(1-factlat)
          vjj(i) = map_v(i+1)*factlat + map_v(1)*(1-factlat)
          ajj(i) = map_a(i+1)*factlat + map_a(1)*(1-factlat)
         enddo
c pole sud
        elseif (jj.eq.nbp_lat-1) then
         do i=1,nbp_lon
          k = (jj-2)*nbp_lon+1+i
          ujj(i) = map_u(klon)*factlat + map_u(k)*(1-factlat)
          vjj(i) = map_v(klon)*factlat + map_v(k)*(1-factlat)
          ajj(i) = map_a(klon)*factlat + map_a(k)*(1-factlat)
         enddo
c autres latitudes
        else
         do i=1,nbp_lon
          k = (jj-2)*nbp_lon+1+i
          ujj(i) = map_u(k+nbp_lon)*factlat + map_u(k)*(1-factlat)
          vjj(i) = map_v(k+nbp_lon)*factlat + map_v(k)*(1-factlat)
          ajj(i) = map_a(k+nbp_lon)*factlat + map_a(k)*(1-factlat)
         enddo
        endif
        ujj(nbp_lon+1)=ujj(1)
        vjj(nbp_lon+1)=vjj(1)
        ajj(nbp_lon+1)=ajj(1)

c Interpolation in longitudes
c-------------------------------------------------
        ii=1  ! lon=-180
        do i=2,nbp_lon
          if (lam.gt.longit(i)) ii=i
        enddo
        factlon = (lam-longit(ii))/(longit(ii+1)-longit(ii))
        ubal    = ujj(ii+1)*factlon + ujj(ii)*(1-factlon)
        vbal    = vjj(ii+1)*factlon + vjj(ii)*(1-factlon)
        abal    = ajj(ii+1)*factlon + ajj(ii)*(1-factlon)

      end

c======================================================================
c======================================================================
c======================================================================

      subroutine pos_implem(phi,lam,ubal,vbal,dt)

      use dimphy
      implicit none

c======================================================================
c Auteur: S. Lebonnois (LMD/CNRS) date: 20091201
c Object: implementation of balloon position.
C======================================================================
c Explicit Arguments:
c ==================
c phi   ---R: Latitude of balloon in radians
c lam   ---R: Longitude of balloon in radians
c ubal  ---R: zonal speed of balloon
c vbal  ---R: meridional speed of balloon
c dt    ---R: time step
c======================================================================

#include "YOMCST.h"
c
c ARGUMENTS
c
      real phi,lam,ubal,vbal,abal,dt

c incrementation longitude
        lam = lam + ubal*dt/(RA*cos(phi))
c maintenue entre -PI et PI:
        do while (lam.ge.RPI)      
              lam=lam-2*RPI
        enddo
        do while (lam.lt.(-1.*RPI))  
              lam=lam+2*RPI
        enddo
c incrementation latitude
        phi = phi + vbal*dt/RA
c maintenue entre -PI/2 et PI/2:
        if (phi.ge.( 0.5*RPI)) phi=    RPI-phi 
        if (phi.le.(-0.5*RPI)) phi=-1.*RPI-phi 

      end
