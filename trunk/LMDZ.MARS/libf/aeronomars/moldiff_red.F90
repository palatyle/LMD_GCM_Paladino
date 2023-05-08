subroutine moldiff_red(ngrid,nlayer,nq,pplay,pplev,pt,pdt,pq,pdq,&
                       ptimestep,zzlay,pdteuv,pdtconduc,pdqdiff)

use tracer_mod, only: noms, mmol
use geometry_mod, only: cell_area

implicit none

include "diffusion.h"

! July 2014 JYC ADD BALISTIC Transport coupling to compute wup for H and H2



!
! Input/Output
!
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of advected tracers
      real ptimestep
      real pplay(ngrid,nlayer)
      real zzlay(ngrid,nlayer)
      real pplev(ngrid,nlayer+1)
      real pq(ngrid,nlayer,nq)
      real pdq(ngrid,nlayer,nq)
      real pt(ngrid,nlayer)
      real pdt(ngrid,nlayer)
      real pdteuv(ngrid,nlayer)
      real pdtconduc(ngrid,nlayer)
      real pdqdiff(ngrid,nlayer,nq)
!
! Local
!
      

!      real hco2(ncompdiff),ho

      integer,dimension(nq) :: indic_diff
      integer ig,iq,nz,l,k,n,nn,p,ij0
      integer istep,il,gcn,ntime,nlraf
      real*8 masse
      real*8 masseU,kBolt,RgazP,Rmars,Grav,Mmars
      real*8 rho0,D0,T0,H0,time0,dZ,time,dZraf,tdiff,Zmin,Zmax
      real*8 FacEsc,invsgmu,PhiEscH,PhiEscH2
      real*8 hp(nlayer)
      real*8 pp(nlayer)
      real*8 pint(nlayer)
      real*8 tt(nlayer),tnew(nlayer),tint(nlayer)
      real*8 zz(nlayer)
      real*8,dimension(:,:),allocatable,save :: qq,qnew,qint,FacMass
      real*8,dimension(:,:),allocatable,save :: rhoK,rhokinit
      real*8 rhoT(nlayer)
      real*8 dmmeandz(nlayer)
      real*8 massemoy(nlayer)
      real*8,dimension(:),allocatable :: Praf,Traf,Rraf,Mraf,Nraf,Draf,Hraf,Wraf
      real*8,dimension(:),allocatable :: Zraf,Tdiffraf
      real*8,dimension(:),allocatable :: Prafold,Mrafold
      real*8,dimension(:,:),allocatable :: Qraf,Rrafk,Nrafk
      real*8,dimension(:,:),allocatable :: Rrafkold
      real*8,dimension(:,:),allocatable :: Drafmol,Hrafmol,Wrafmol,Tdiffrafmol
      real*8,dimension(:),allocatable :: Atri,Btri,Ctri,Dtri,Xtri,Tad,Dad,Zad,rhoad
      real*8,dimension(:),allocatable :: alpha,beta,gama,delta,eps
      real*8,dimension(:),allocatable,save ::  wi,Wad,Uthermal,Lambdaexo,Hspecie
      real*8,dimension(:),allocatable,save :: Mtot1,Mtot2,Mraf1,Mraf2
      integer,parameter :: ListeDiffNb=15
      character(len=20),dimension(ListeDiffNb) :: ListeDiff
      real*8,parameter :: pi=3.141592653
      real*8,parameter :: g=3.72d0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     tracer numbering in the molecular diffusion
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  We need the index of escaping species 

      integer,save :: i_h2  
      integer,save :: i_h
! vertical index limit for the molecular diffusion
      integer,save :: il0  

! Tracer indexes in the GCM:
!      integer,save :: g_co2=0
!      integer,save :: g_co=0
!      integer,save :: g_o=0
!      integer,save :: g_o1d=0
!      integer,save :: g_o2=0
!      integer,save :: g_o3=0
!      integer,save :: g_h=0
!      integer,save :: g_h2=0
!      integer,save :: g_oh=0
!      integer,save :: g_ho2=0
!      integer,save :: g_h2o2=0
!      integer,save :: g_n2=0
!      integer,save :: g_ar=0
!      integer,save :: g_h2o=0
!      integer,save :: g_n=0
!      integer,save :: g_no=0
!      integer,save :: g_no2=0
!      integer,save :: g_n2d=0
!      integer,save :: g_oplus=0
!      integer,save :: g_hplus=0

      integer,save :: ncompdiff
      integer,dimension(:),allocatable,save :: gcmind ! array of GCM indexes
      integer ierr,compteur

      logical,save :: firstcall=.true.
!      real abfac(ncompdiff)
      real,dimension(:,:),allocatable,save :: dij
      real,save :: step


! Initializations at first call
      if (firstcall) then

! Liste des especes qui sont diffuses si elles sont presentes
!	ListeDiff=['co2','o','n2','ar','co','h2','h','d2','d','hd','o2','oh','ho2','h2o_vap','h2o2','o1d','n','no','no2']
	ListeDiff(1)='co2'
	ListeDiff(2)='o'
	ListeDiff(3)='n2'
	ListeDiff(4)='ar'
	ListeDiff(5)='co'
	ListeDiff(6)='h2'
	ListeDiff(7)='h'
	ListeDiff(8)='d2'
	ListeDiff(9)='hd'
	ListeDiff(10)='d'
	ListeDiff(11)='o2'
	ListeDiff(12)='h2o_vap'
	ListeDiff(13)='o3'
	ListeDiff(14)='n'
        ListeDiff(15)='he'
	i_h=1000
	i_h2=1000
! On regarde quelles especes diffusables sont presentes

	ncompdiff=0
	indic_diff(1:nq)=0
	
	do nn=1,nq
	do n=1,ListeDiffNb
	if (ListeDiff(n) .eq. noms(nn)) then
	indic_diff(nn)=1
!	if (noms(nn) .eq. 'h') i_h=n
!	if (noms(nn) .eq. 'h2') i_h2=n
	endif

	enddo
	if (indic_diff(nn) .eq. 1) then
	print*,'specie ', noms(nn), 'diffused in moldiff_red'
	ncompdiff=ncompdiff+1
	endif
	enddo
	print*,'number of diffused species:',ncompdiff
	allocate(gcmind(ncompdiff))

! Store gcm indexes in gcmind
	n=0
	do nn=1,nq
	if (indic_diff(nn) .eq. 1) then
	n=n+1
	gcmind(n)=nn
	if (noms(nn) .eq. 'h') i_h=n
        if (noms(nn) .eq. 'h2') i_h2=n
	endif
	enddo

!	print*,'gcmind',gcmind,i_h,i_h2

! find vertical index above which diffusion is computed

	do l=1,nlayer
	if (pplay(1,l) .gt. Pdiff) then
	il0=l
	endif
	enddo
	il0=il0+1
	print*,'vertical index for diffusion',il0,pplay(1,il0)

	allocate(dij(ncompdiff,ncompdiff))

        call moldiffcoeff_red(dij,indic_diff,gcmind,ncompdiff)
        print*,'MOLDIFF  EXO'

! allocatation des tableaux dependants du nombre d especes diffusees
	allocate(qq(nlayer,ncompdiff))
	allocate(qnew(nlayer,ncompdiff))
	allocate(qint(nlayer,ncompdiff))
	allocate(FacMass(nlayer,ncompdiff))
	allocate(rhok(nlayer,ncompdiff))
	allocate(rhokinit(nlayer,ncompdiff))

	allocate(wi(ncompdiff))
	allocate(wad(ncompdiff))
	allocate(uthermal(ncompdiff))
	allocate(lambdaexo(ncompdiff))
	allocate(Hspecie(ncompdiff))
	allocate(Mtot1(ncompdiff))
	allocate(Mtot2(ncompdiff))
	allocate(Mraf1(ncompdiff))
	allocate(Mraf2(ncompdiff))

        firstcall= .false.
	step=1
      endif ! of if (firstcall)

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      masseU=1.660538782d-27
      kbolt=1.3806504d-23
      RgazP=8.314472 	
!      PI =3.141592653
!      g=3.72d0
      Rmars=3390000d0 ! Used to compute escape parameter no need to be more accurate
	Grav=6.67d-11
	Mmars=6.4d23
	ij0=6000 ! For test
	invsgmu=1d0/g/masseU
	
! Compute the wup(ig) for H and H2 using the balistic code from R Yelle

	PhiEscH=0D0
	PhiEscH2=0D0

      do ig=1,ngrid
	pp=dble(pplay(ig,:))

! Update the temperature

!	CALL TMNEW(pt(ig,:),pdt(ig,:),pdtconduc(ig,:),pdteuv(ig,:)      & 
!     &  ,tt,ptimestep,nlayer,ig)
        do l=1,nlayer
          tt(l)=pt(ig,l)*1D0+(pdt(ig,l)*dble(ptimestep)+                &
                              pdtconduc(ig,l)*dble(ptimestep)+          &
                              pdteuv(ig,l)*dble(ptimestep))
          ! to cach Nans...
          if (tt(l).ne.tt(l)) then
            print*,'Err TMNEW',ig,l,tt(l),pt(ig,l),                     &
              pdt(ig,l),pdtconduc(ig,l),pdteuv(ig,l),dble(ptimestep)
          endif
        enddo ! of do l=1,nlayer

! Update the mass mixing ratios modified by other processes

!	CALL QMNEW(pq(ig,:,:),pdq(ig,:,:),qq,ptimestep,nlayer,        &
!     &  ncompdiff,gcmind,ig)
        do iq=1,ncompdiff
         do l=1,nlayer
          qq(l,iq)=pq(ig,l,gcmind(iq))*1D0+(                            &
                           pdq(ig,l,gcmind(iq))*dble(ptimestep))
          qq(l,iq)=max(qq(l,iq),1d-30)
         enddo ! of do l=1,nlayer
        enddo ! of do iq=1,ncompdiff

! Compute the Pressure scale height

	CALL HSCALE(pp,hp,nlayer)

! Compute the atmospheric mass (in Dalton)

	CALL MMOY(massemoy,mmol,qq,gcmind,nlayer,ncompdiff)

! Compute the vertical gradient of atmospheric mass

	CALL DMMOY(massemoy,hp,dmmeandz,nlayer)

! Compute the altitude of each layer

	CALL ZVERT(pp,tt,massemoy,zz,nlayer,ig)

! Compute the total mass density (kg/m3)
	
	CALL RHOTOT(pp,tt,massemoy,qq,RHOT,RHOK,nlayer,ncompdiff)
	RHOKINIT=RHOK

! Compute total mass of each specie before new grid
! For conservation tests 
! The conservation is not always fulfilled especially 
! for species very far from diffusion equilibrium "photochemical species"
! e.g. OH, O(1D)

	Mtot1(1:ncompdiff)=0d0

	do l=il0,nlayer
	do nn=1,ncompdiff
	Mtot1(nn)=Mtot1(nn)+1d0/g*qq(l,nn)*                             &
     &  (dble(pplev(ig,l))-dble(pplev(ig,l+1)))
	enddo
	enddo

	Zmin=zz(il0)
	Zmax=zz(nlayer)


! If Zmax > 4000 km there is a problem / stop 

	if (Zmax .gt. 4000000.) then
	Print*,'Zmax too high',ig,zmax,zmin
	do l=1,nlayer
	print*,'old',zz(l),pt(ig,l),pdteuv(ig,l),pdq(ig,l,:)
	print*,'l',l,rhot(l),tt(l),pp(l),massemoy(l),qq(l,:)
	enddo
	stop
	endif

! The number of diffusion layers is variable 
! and fixed by the resolution (dzres) specified in diffusion.h
! I fix a minimum  number of layers 40 

	nlraf=int((Zmax-Zmin)/1000./dzres)+1
	if (nlraf .le. 40)  nlraf=40

!	if (nlraf .ge. 200) print*,ig,nlraf,Zmin,Zmax
        
         ! allocate arrays:
	 allocate(Praf(nlraf),Traf(nlraf),Rraf(nlraf),Mraf(nlraf))
         allocate(Nraf(nlraf),Draf(nlraf),Hraf(nlraf),Wraf(nlraf))
         allocate(Zraf(nlraf),Tdiffraf(nlraf))
         allocate(Prafold(nlraf),Mrafold(nlraf))
         allocate(Qraf(nlraf,ncompdiff),Rrafk(nlraf,ncompdiff),Nrafk(nlraf,ncompdiff))
         allocate(Rrafkold(nlraf,ncompdiff))
         allocate(Drafmol(nlraf,ncompdiff),Hrafmol(nlraf,ncompdiff))
         allocate(Wrafmol(nlraf,ncompdiff),Tdiffrafmol(nlraf,ncompdiff))
         allocate(Atri(nlraf),Btri(nlraf),Ctri(nlraf),Dtri(nlraf),Xtri(nlraf))
         allocate(Tad(nlraf),Dad(nlraf),Zad(nlraf),rhoad(nlraf))
         allocate(alpha(nlraf),beta(nlraf),gama(nlraf),delta(nlraf),eps(nlraf))

! before beginning, I use a better vertical resolution above il0, 
! altitude grid reinterpolation 
! The diffusion is solved on an altitude grid with constant step dz

	CALL UPPER_RESOL(pp,tt,zz,massemoy,RHOT,RHOK,                   &
     &  qq,mmol,gcmind,Praf,Traf,Qraf,Mraf,Zraf,                        &
     &  Nraf,Nrafk,Rraf,Rrafk,il0,nlraf,ncompdiff,nlayer,ig)

	Prafold=Praf
	Rrafkold=Rrafk
	Mrafold=Mraf
	
! We reddo interpolation of the gcm grid to estimate mass loss due to interpolation processes.

	CALL GCMGRID_P(Zraf,Praf,Qraf,Traf,Nrafk,Rrafk,qq,qint,tt,tint  &
     &    ,pp,mmol,gcmind,nlraf,ncompdiff,nlayer,ig)

! We compute the mass correction factor of each specie at each pressure level

	CALL CORRMASS(qq,qint,FacMass,nlayer,ncompdiff)

! Altitude step

	Dzraf=Zraf(2)-Zraf(1)

! Total mass computed on the altitude grid

	Mraf1(1:ncompdiff)=0d0
	do nn=1,ncompdiff
	do l=1,nlraf
	Mraf1(nn)=Mraf1(nn)+Rrafk(l,nn)*Dzraf
	enddo
	enddo

! Reupdate values for mass conservation

!	do l=1,nlraf
!	print*,'test',l,Nraf(l),Praf(l)
!	do nn=1,ncompdiff
!	Rrafk(l,nn)=RrafK(l,nn)*Mtot1(nn)/Mraf1(nn)
!	enddo
!	Rraf(l)=sum(Rrafk(l,:))
!	do nn=1,ncompdiff
!	Qraf(l,nn)=Rrafk(l,nn)/Rraf(l)
!	Nrafk(l,nn)=Rrafk(l,nn)/dble(mmol(gcmind(nn)))/masseU
!	enddo
!	Mraf(l)=1d0/sum(Qraf(l,:)/dble(mmol(gcmind(:))))
!	Nraf(l)=Rraf(l)/Mraf(l)/masseU
!	Praf(l)=kbolt*Traf(l)*Nraf(l)
!	print*,'test',l,Nraf(l),Praf(l)
!	enddo

!	do l=1,nlayer
!	print*,'l',l,zz(l),pp(l),tt(l),sum(qq(l,:)),massemoy(l)
!	enddo

! The diffusion is computed above il0 computed from Pdiff in diffusion.h
! No change below il0
	
	do l=1,nlayer
         qnew(l,:)=qq(l,:) ! No effet below il0
       enddo

! all species treated independently

! Upper boundary velocity
! Jeans escape for H and H2 
! Comparison with and without escape still to be done
! No escape for other species


	do nn=1,ncompdiff
	Uthermal(nn)=sqrt(2d0*kbolt*Traf(nlraf)/masseU/                 &
     &  dble(mmol(gcmind(nn))))
	Lambdaexo(nn)=masseU*dble(mmol(gcmind(nn)))*Grav*Mmars/         &
     &  (Rmars+Zraf(nlraf))/kbolt/Traf(nlraf)
	wi(nn)=0D0
	if (nn .eq. i_h .or. nn .eq. i_h2) then 
	wi(nn)=Uthermal(nn)/2d0/sqrt(PI)*exp(-Lambdaexo(nn))*           &
     &  (Lambdaexo(nn)+1d0)
	endif
	enddo

!	print*,'wi',wi(i_h),wi(i_h2),Uthermal,Lambdaexo,mmol(gcmind(:))
!	print*,'wi',wi
!	stop

! Compute time step for diffusion

! Loop on species

	T0=Traf(nlraf)
	rho0=1d0

        do nn=1,ncompdiff
	masse=dble(mmol(gcmind(nn)))
! DIffusion coefficient
	CALL DCOEFF(nn,dij,Praf,Traf,Nraf,Nrafk,Draf,nlraf,ncompdiff)
	Drafmol(:,nn)=Draf(:)
! Scale height of the density of the specie
	CALL HSCALEREAL(nn,Nrafk,Dzraf,Hraf,nlraf,ncompdiff)
	Hrafmol(:,nn)=Hraf(:)
!	Hspecie(nn)=kbolt*T0/masse*invsgmu
! Computation of the diffusion vertical velocity of the specie
	CALL VELVERT(nn,Traf,Hraf,Draf,Dzraf,masse,Wraf,nlraf)
	Wrafmol(:,nn)=Wraf(:)
! Computation of the diffusion time 
	CALL TIMEDIFF(nn,Hraf,Wraf,Tdiffraf,nlraf)
	Tdiffrafmol(:,nn)=Tdiffraf
	enddo
! We use a lower time step
	Tdiff=minval(Tdiffrafmol)
	Tdiff=minval(Tdiffrafmol(nlraf,:))*Mraf(nlraf)
! Some problems when H is dominant 
! The time step is chosen function of atmospheric mass at higher level 
! It is not very nice

!	if (ig .eq. ij0) then 
!	print*,'test',ig,tdiff,tdiffmin,minloc(Tdiffrafmol),minloc(Tdiffrafmol(nlraf,:))
!	endif
	if (tdiff .lt. tdiffmin*Mraf(nlraf)) tdiff=tdiffmin*Mraf(nlraf)

! Number of time step
	ntime=anint(dble(ptimestep)/tdiff)
!	print*,'ptime',ig,ptimestep,tdiff,ntime,tdiffmin,Mraf(nlraf)
! Adimensionned temperature

	do l=1,nlraf
	Tad(l)=Traf(l)/T0
	enddo

	do istep=1,ntime
	do nn=1,ncompdiff

	Draf(1:nlraf)=Drafmol(1:nlraf,nn)

! Parameters to adimension the problem 

	H0=kbolt*T0/dble(mmol(gcmind(nn)))*invsgmu
	D0=Draf(nlraf)
	Time0=H0*H0/D0
	Time=Tdiff/Time0

! Adimensions

	do l=1,nlraf
	Dad(l)=Draf(l)/D0
	Zad(l)=Zraf(l)/H0
	enddo
	Wad(nn)=wi(nn)*Time0/H0
	DZ=Zad(2)-Zad(1)
	FacEsc=exp(-wad(nn)*DZ/Dad(nlraf-1))

	do l=1,nlraf
	RhoAd(l)=Rrafk(l,nn)/rho0
	enddo

! Compute intermediary coefficients

	CALL DIFFPARAM(Tad,Dad,DZ,alpha,beta,gama,delta,eps,nlraf)

! First way to include escape need to be validated
! Compute escape factor (exp(-ueff*dz/D(nz)))

! Compute matrix coefficients

	CALL MATCOEFF(alpha,beta,gama,delta,eps,Dad,rhoAd,FacEsc,       &
     &   dz,time,Atri,Btri,Ctri,Dtri,nlraf)

	Xtri(:)=0D0

! SOLVE SYSTEM

	CALL tridag(Atri,Btri,Ctri,Dtri,Xtri,nlraf)

!	Xtri=rhoAd

!	if (ig .eq. ij0 .and. (nn .eq. 1 .or. nn .eq. 5 .or. nn .eq. 6 .or. nn .eq. 16)) then
!	do l=1,nlraf
!	if (Xtri(l) .lt. 0.) then
!	print*,'l',l,rhoAd(l)*rho0,Xtri(l)*rho0,nn,Tad(l),Zad(l),Dad(l)
!	stop
!	endif
!	enddo
!	endif

! Check mass conservation of speci

!	CALL CheckMass(rhoAd,Xtri,nlraf,nn)

! Update mass density of the species

	do l=1,nlraf
	Rrafk(l,nn)=rho0*Xtri(l)
	if (Rrafk(l,nn) .ne. Rrafk(l,nn) .or.                           &
     &  Rrafk(l,nn) .lt. 0 .and. nn .eq. 16) then 

! Test if n(CO2) < 0 skip diffusion (should never happen)

	print*,'pb moldiff',istep,ig,l,nn,Rrafk(l,nn),tdiff,& 
     &  rho0*Rhoad(l),Zmin,Zmax,ntime
	print*,'Atri',Atri
	print*,'Btri',Btri
	print*,'Ctri',Ctri
	print*,'Dtri',Dtri
	print*,'Xtri',Xtri
	print*,'alpha',alpha
	print*,'beta',beta
	print*,'gama',gama
	print*,'delta',delta
	print*,'eps',eps
	print*,'Dad',Dad
	print*,'Temp',Traf
	print*,'alt',Zraf
	print*,'Mraf',Mraf
	stop
!	pdqdiff(1:ngrid,1:nlayer,1:nq)=0.
!	return
!	Rrafk(l,nn)=1D-30*Rraf(l)
	Rrafk(l,nn)=rho0*Rhoad(l)

	endif

	enddo

	enddo ! END Species Loop

! Update total mass density

	do l=1,nlraf
	Rraf(l)=sum(Rrafk(l,:))
	enddo

! Compute new mass average at each altitude level and new pressure

	do l=1,nlraf
	do nn=1,ncompdiff
	Qraf(l,nn)=Rrafk(l,nn)/Rraf(l)
	Nrafk(l,nn)=Rrafk(l,nn)/dble(mmol(gcmind(nn)))/masseU
	enddo
	Mraf(l)=1d0/sum(Qraf(l,:)/dble(mmol(gcmind(:))))
	Nraf(l)=Rraf(l)/Mraf(l)/masseU
	Praf(l)=Nraf(l)*kbolt*Traf(l)
	enddo

	enddo ! END time Loop

! Compute the total mass of each species to check mass conservation	

	Mraf2(1:ncompdiff)=0d0
        do nn=1,ncompdiff
        do l=1,nlraf
        Mraf2(nn)=Mraf2(nn)+Rrafk(l,nn)*Dzraf
        enddo
        enddo

!        print*,'Mraf',Mraf2

! Reinterpolate values on the GCM pressure levels

	CALL GCMGRID_P2(Zraf,Praf,Qraf,Traf,Nrafk,Rrafk,qq,qnew,tt,tnew,&
     &   pp,mmol,gcmind,nlraf,ncompdiff,nlayer,FacMass,ig)

	CALL RHOTOT(pp,tt,massemoy,qnew,RHOT,RHOK,nlayer,ncompdiff)

! Update total escape flux of H and H2 (if q was really qnew, but not forget we will output
! the trend only at the end

	PhiEscH=PhiEscH+wi(i_h)*Nrafk(nlraf,i_h)*cell_area(ig) ! in s-1
	PhiEscH2=PhiEscH2+wi(i_h2)*Nrafk(nlraf,i_h2)*cell_area(ig) ! in s-1 (U in m/s, aire in m2, Nrafk in m-3)
!	print*,'test',ig,wi(i_h),Nrafk(nlraf,i_h),wi(i_h2),Nrafk(nlraf,i_h2),cell_area(ig),PhiEscH,PhiEscH2,i_h,i_h2
!	stop


	if (ig .eq. ij0) then 
	do l=il0,nlayer
        write(*,'(i2,1x,19(e12.4,1x))') l,zz(l),tt(l),RHOK(l,1)/sum(RHOK(l,:)),RHOKINIT(l,1)/sum(RHOKINIT(l,:)),&
     &  RHOK(l,2)/sum(RHOK(l,:)),RHOKINIT(l,2)/sum(RHOKINIT(l,:)),&
     &  RHOK(l,6)/sum(RHOK(l,:)),RHOKINIT(l,6)/sum(RHOKINIT(l,:)),&
     &  RHOK(l,5)/sum(RHOK(l,:)),RHOKINIT(l,5)/sum(RHOKINIT(l,:)),&
     &  RHOK(l,7)/sum(RHOK(l,:)),RHOKINIT(l,7)/sum(RHOKINIT(l,:))
        enddo
	endif

! Compute total mass of each specie on the GCM vertical grid

	Mtot2(1:ncompdiff)=0d0

        do l=il0,nlayer
        do nn=1,ncompdiff
        Mtot2(nn)=Mtot2(nn)+1d0/g*qnew(l,nn)*                           &
     &  (dble(pplev(ig,l))-dble(pplev(ig,l+1)))
        enddo
        enddo

! Check mass  conservation of each specie on column

!	do nn=1,ncompdiff
!	CALL CheckMass2(qq,qnew,pplev(ig,:),il0,nlayer,nn,ncompdiff)
!	enddo

! Compute the diffusion trends du to diffusion

	do l=1,nlayer
	do nn=1,ncompdiff
	pdqdiff(ig,l,gcmind(nn))=(qnew(l,nn)-qq(l,nn))/ptimestep
	enddo
	enddo

! deallocation des tableaux

	deallocate(Praf,Traf,Rraf,Mraf)
        deallocate(Nraf,Draf,Hraf,Wraf)
        deallocate(Zraf,Tdiffraf)
        deallocate(Prafold,Mrafold)
        deallocate(Qraf,Rrafk,Nrafk)
        deallocate(Rrafkold)
        deallocate(Drafmol,Hrafmol)
        deallocate(Wrafmol,Tdiffrafmol)
        deallocate(Atri,Btri,Ctri,Dtri,Xtri)
        deallocate(Tad,Dad,Zad,rhoad)
        deallocate(alpha,beta,gama,delta,eps)


       enddo  ! ig loop		
!	print*,'Escape flux H, H2 (s-1)',PhiEscH,PhiEscH2

      return
      end

!    ********************************************************************
!    ********************************************************************
!    ********************************************************************
 
!     JYC subtroutine solving MX = Y where M is defined as a block tridiagonal 
!	matrix (Thomas algorithm), tested on a example

      subroutine tridagbloc(M,F,X,n1,n2)
      parameter (nmax=100)
      real*8 M(n1*n2,n1*n2),F(n1*n2),X(n1*n2)
      real*8 A(n1,n1,n2),B(n1,n1,n2),C(n1,n1,n2),D(n1,n2)
      real*8 at(n1,n1),bt(n1,n1),ct(n1,n1),dt(n1),gamt(n1,n1),y(n1,n1)
      real*8 alf(n1,n1),gam(n1,n1,n2),alfinv(n1,n1)
      real*8 uvec(n1,n2),uvect(n1),vvect(n1),xt(n1)
      real*8 indx(n1)
      real*8 rhu
      integer n1,n2
      integer i,p,q    

	X(:)=0.
! Define the bloc matrix A,B,C and the vector D
      A(1:n1,1:n1,1)=M(1:n1,1:n1)
      C(1:n1,1:n1,1)=M(1:n1,n1+1:2*n1)
      D(1:n1,1)=F(1:n1)

      do i=2,n2-1
      A(1:n1,1:n1,i)=M((i-1)*n1+1:i*n1,(i-1)*n1+1:i*n1)
      B(1:n1,1:n1,i)=M((i-1)*n1+1:i*n1,(i-2)*n1+1:(i-1)*n1)
      C(1:n1,1:n1,i)=M((i-1)*n1+1:i*n1,i*n1+1:(i+1)*n1)
      D(1:n1,i)=F((i-1)*n1+1:i*n1)
      enddo
      A(1:n1,1:n1,n2)=M((n2-1)*n1+1:n2*n1,(n2-1)*n1+1:n2*n1)
      B(1:n1,1:n1,n2)=M((n2-1)*n1+1:n2*n1,(n2-2)*n1+1:(n2-1)*n1)
      D(1:n1,n2)=F((n2-1)*n1+1:n2*n1)

! Initialization
	y(:,:)=0. 
	do i=1,n1
	y(i,i)=1.
	enddo
       
	at(:,:)=A(:,:,1)
        ct(:,:)=C(:,:,1)
        dt(:)=D(:,1)
        call ludcmp(at,n1,n1,indx,rhu,ierr)
        do p=1,n1
         call lubksb(at,n1,n1,indx,y(1,p))
         do q=1,n1
         alfinv(q,p)=y(q,p)
         enddo
        enddo
      gamt=matmul(alfinv,ct)
      gam(:,:,1)=gamt(:,:)     
      uvect=matmul(alfinv,dt)
      uvec(:,1)=uvect

      do i=2,n2-1
        y(:,:)=0.
       do j=1,n1
         y(j,j)=1.
       enddo  
      bt(:,:)=B(:,:,i)
      at(:,:)=A(:,:,i)-matmul(bt,gamt)
      ct(:,:)=C(:,:,i)
      dt(:)=D(:,i)
      call ludcmp(at,n1,n1,indx,rhu,ierr)
        do p=1,n1
         call lubksb(at,n1,n1,indx,y(1,p))
         do q=1,n1
         alfinv(q,p)=y(q,p)
         enddo
        enddo
      gamt=matmul(alfinv,ct)
      gam(:,:,i)=gamt
      vvect=dt-matmul(bt,uvect)
      uvect=matmul(alfinv,vvect)
      uvec(:,i)=uvect 
      enddo
      bt=B(:,:,n2)
      dt=D(:,n2)
      at=A(:,:,n2)-matmul(bt,gamt)
      vvect=dt-matmul(bt,uvect)
       y(:,:)=0.
       do j=1,n1
         y(j,j)=1.
       enddo
       call ludcmp(at,n1,n1,indx,rhu,ierr)
        do p=1,n1
         call lubksb(at,n1,n1,indx,y(1,p))
         do q=1,n1
         alfinv(q,p)=y(q,p)
         enddo
        enddo
      xt=matmul(alfinv,vvect)
      X((n2-1)*n1+1 :n1*n2)=xt
      do i=n2-1,1,-1
      gamt=gam(:,:,i)
      xt=X(i*n1+1:n1*n2)
      uvect=uvec(:,i) 
      vvect=matmul(gamt,xt)
      X((i-1)*n1+1:i*n1)=uvect-vvect
      enddo

      end

      subroutine tridag(a,b,c,r,u,n)
!      parameter (nmax=4000)
!      dimension gam(nmax),a(n),b(n),c(n),r(n),u(n)
      real*8 gam(n),a(n),b(n),c(n),r(n),u(n)
      if(b(1).eq.0.)then
        stop 'tridag: error: b(1)=0 !!! '
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          stop 'tridag: error: bet=0 !!! '
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      end

!    ********************************************************************
!    ********************************************************************
!    ********************************************************************

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      implicit none

      integer i,j,n,np,ii,ll
      real*8 sum
      real*8 a(np,np),indx(np),b(np)

!      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

!    ********************************************************************
!    ********************************************************************
!    ********************************************************************

      SUBROUTINE LUDCMP(A,N,NP,INDX,D,ierr)

      implicit none

      integer n,np,nmax,i,j,k,imax
      real*8 d,tiny,aamax
      real*8 a(np,np),indx(np)
      integer ierr  ! error =0 if OK, =1 if problem

      PARAMETER (NMAX=100,TINY=1.0E-20)
!      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      real*8 sum,vv(nmax),dum

      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) then
            write(*,*) 'In moldiff: Problem in LUDCMP with matrix A'
            write(*,*) 'Singular matrix ?'
            write(*,*) 'Matrix A = ', A
!            stop
!           TO DEBUG :
            ierr =1
            return
!           stop
        END IF 

        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      ierr =0
      RETURN
      END

	SUBROUTINE TMNEW(T1,DT1,DT2,DT3,T2,dtime,nl,ig)
	IMPLICIT NONE

	INTEGER,INTENT(IN) :: nl,ig
	REAL,INTENT(IN),DIMENSION(nl) :: T1,DT1,DT2,DT3
	REAL*8,INTENT(OUT),DIMENSION(nl) :: T2
	REAL,INTENT(IN) :: dtime
	INTEGER :: l
	DO l=1,nl
	T2(l)=T1(l)*1D0+(DT1(l)*dble(dtime)+                            &
     &  DT2(l)*dble(dtime)+                                             &
     &  DT3(l)*dble(dtime))*1D0
	if (T2(l) .ne. T2(l)) then
	print*,'Err TMNEW',ig,l,T2(l),T1(l),dT1(l),DT2(l),              &
     &  DT3(l),dtime,dble(dtime)
	endif

	ENDDO
	END

	SUBROUTINE QMNEW(Q1,DQ,Q2,dtime,nl,nq,gc,ig)
        use tracer_mod, only: nqmx
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: nl,nq
        INTEGER,INTENT(IN) :: ig
	INTEGER,INTENT(IN),dimension(nq) :: gc
        REAL,INTENT(IN),DIMENSION(nl,nqmx) :: Q1,DQ
        REAL*8,INTENT(OUT),DIMENSION(nl,nq) :: Q2
        REAL,INTENT(IN) :: dtime
        INTEGER :: l,iq
        DO l=1,nl
	DO iq=1,nq
        Q2(l,iq)=Q1(l,gc(iq))*1D0+(DQ(l,gc(iq))*dble(dtime))*1D0
	Q2(l,iq)=max(Q2(l,iq),1d-30)
        ENDDO
	ENDDO
	END

	SUBROUTINE HSCALE(p,hp,nl)
	IMPLICIT NONE

	INTEGER :: nl
	INTEGER :: l
	REAL*8,dimension(nl) :: P
	REAL*8,DIMENSION(nl) :: Hp
	
	hp(1)=-log(P(2)/P(1))
        hp(nl)=-log(P(nl)/P(nl-1))

	DO l=2,nl-1
	hp(l)=-log(P(l+1)/P(l-1))
	ENDDO
	END

	SUBROUTINE MMOY(massemoy,mmol,qq,gc,nl,nq)
	use tracer_mod, only: nqmx
        IMPLICIT NONE

	INTEGER :: nl,nq,l
	INTEGER,dimension(nq) :: gc
	REAL*8,DIMENSION(nl,nq) :: qq
	REAL*8,DIMENSION(nl) :: massemoy
	REAL,DIMENSION(nqmx) :: MMOL


	do l=1,nl
	massemoy(l)=1D0/sum(qq(l,:)/dble(mmol(gc(:))))
	enddo

	END

	SUBROUTINE DMMOY(M,H,DM,nl)
	IMPLICIT NONE	
	INTEGER :: nl,l
	REAL*8,DIMENSION(nl) :: M,H,DM

	DM(1)=(-3D0*M(1)+4D0*M(2)-M(3))/2D0/H(1)
	DM(nl)=(3D0*M(nl)-4D0*M(nl-1)+M(nl-2))/2D0/H(nl)

	do l=2,nl-1
        DM(l)=(M(l+1)-M(l-1))/H(l)
        enddo

	END

	SUBROUTINE ZVERT(P,T,M,Z,nl,ig)
	IMPLICIT NONE
	INTEGER :: nl,l,ig
	REAL*8,dimension(nl) :: P,T,M,Z,H
	REAL*8 :: z0
	REAL*8 :: kbolt,masseU,Konst,g,Hpm
	masseU=1.660538782d-27
        kbolt=1.3806504d-23
	Konst=kbolt/masseU
	g=3.72D0

	z0=0d0
	Z(1)=z0
	H(1)=Konst*T(1)/M(1)/g

	do l=2,nl
	H(l)=Konst*T(l)/M(l)/g
	Hpm=H(l-1)
	Z(l)=z(l-1)-Hpm*log(P(l)/P(l-1))
	if (Z(l) .ne. Z(l)) then
	print*,'pb',l,ig
	print*,'P',P
	print*,'T',T
	print*,'M',M
	print*,'Z',Z
	print*,'Hpm',Hpm
	endif
	enddo

	END

	SUBROUTINE RHOTOT(P,T,M,qq,rhoN,rhoK,nl,nq)
	IMPLICIT NONE

	REAL*8 :: kbolt,masseU,Konst
	INTEGER :: nl,nq,l,iq
	REAL*8,DIMENSION(nl) :: P,T,M,RHON
	REAL*8,DIMENSION(nl,nq) :: RHOK,qq
	masseU=1.660538782d-27
        kbolt=1.3806504d-23
	Konst=Kbolt/masseU

	do l=1,nl
	RHON(l)=P(l)*M(l)/T(l)/Konst
	do iq=1,nq
	RHOK(l,iq)=qq(l,iq)*RHON(l)
	enddo
	enddo

	END

	SUBROUTINE UPPER_RESOL(P,T,Z,M,R,Rk,                            &
     & qq,mmol,gc,Praf,Traf,Qraf,Mraf,Zraf,                             &
     & Nraf,Nrafk,Rraf,Rrafk,il,nl,nq,nlx,ig)
        use tracer_mod, only: nqmx
	IMPLICIT NONE
	
	INTEGER :: nl,nq,il,l,i,iq,nlx,iz,ig
	INTEGER :: gc(nq)
	INTEGER,DIMENSION(1) :: indz
	REAL*8, DIMENSION(nlx) :: P,T,Z,M,R
	REAL*8, DIMENSION(nlx,nq) :: qq,Rk
	REAL*8, DIMENSION(nl) :: Praf,Traf,Mraf,Zraf,Nraf,Rraf
	REAL*8 :: kbolt,masseU,Konst,g
	REAL*8, DIMENSION(nl,nq) :: Qraf,Rrafk,Nrafk
	REAL*8 :: facZ,dZ,H
	REAL,DIMENSION(nqmx) :: mmol
	masseU=1.660538782d-27
        kbolt=1.3806504d-23
        Konst=Kbolt/masseU
	g=3.72d0


	Zraf(1)=z(il)
	Praf(1)=P(il)
	Traf(1)=T(il)
	Nraf(1)=Praf(1)/kbolt/Traf(1)
	do iq=1,nq
	Qraf(1,iq)=qq(il,iq)
	enddo
	Mraf(1)=1d0/sum(Qraf(1,:)/dble(mmol(gc(:))))
	Rraf(1)=Mraf(1)*masseU*Nraf(1)
	do iq=1,nq
        Rrafk(1,iq)=Rraf(1)*Qraf(1,iq)
	Nrafk(1,iq)=Rrafk(1,iq)/masseU/dble(mmol(gc(iq)))
	enddo
	Zraf(nl)=z(nlx)

	do l=2,nl-1
	Zraf(l)=Zraf(1)+(Zraf(nl)-Zraf(1))/dble(nl-1)*dble((l-1))
	indz=maxloc(z,mask=z <= Zraf(l))
	iz=indz(1)
	if (iz .lt. 1 .or. iz .gt. nlx) then
	print*,'danger',iz,nl,Zraf(l),l,Zraf(1),Zraf(nl),z,P,T,ig
	stop
	endif
	dZ=Zraf(l)-Zraf(l-1)
!	dZ=Zraf(l)-z(iz)
	facz=(Zraf(l)-z(iz))/(z(iz+1)-z(iz))
        Traf(l)=T(iz)+(T(iz+1)-T(iz))*facz
        do iq=1,nq
!        Qraf(l,iq)=qq(iz,iq)+(qq(iz+1,iq)-qq(iz,iq))*facz
	Rrafk(l,iq)=Rk(iz,iq)+(Rk(iz+1,iq)-Rk(iz,iq))*facZ
	Rrafk(l,iq)=Rk(iz,iq)*(Rk(iz+1,iq)/Rk(iz,iq))**facZ
        enddo
!	Mraf(l)=1D0/(sum(qraf(l,:)/dble(mmol(gc(:)))))
	Rraf(l)=sum(Rrafk(l,:))
	do iq=1,nq
	Qraf(l,iq)=Rrafk(l,iq)/Rraf(l)
	enddo
	Mraf(l)=1D0/(sum(qraf(l,:)/dble(mmol(gc(:)))))
!	H=Konst*Traf(l)/Mraf(l)/g
!	H=Konst*T(iz)/M(iz)/g
!	Praf(l)=P(iz)*exp(-dZ/H)
!	Praf(l)=Praf(l-1)*exp(-dZ/H)
!	print*,'iz',l,iz,Praf(il-1)*exp(-dZ/H),z(iz),z(iz+1),H
	Nraf(l)=Rraf(l)/Mraf(l)/masseU
	Praf(l)=Nraf(l)*kbolt*Traf(l)
!	Rraf(l)=Nraf(l)*Mraf(l)*masseU
	do iq=1,nq
!	Rrafk(l,iq)=Rraf(l)*Qraf(l,iq)
	Nrafk(l,iq)=Rrafk(l,iq)/dble(mmol(gc(iq)))/masseU
	if (Nrafk(l,iq) .lt. 0. .or.                                    & 
     &  Nrafk(l,iq) .ne. Nrafk(l,iq)) then
	print*,'pb interpolation',l,iq,Nrafk(l,iq),Rrafk(l,iq),         &
     &  Qraf(l,iq),Rk(iz,iq),Rk(iz+1,iq),facZ,Zraf(l),z(iz)
	stop
	endif
        enddo
	enddo
	Zraf(nl)=z(nlx)
	Traf(nl)=T(nlx)
	do iq=1,nq
	Rrafk(nl,iq)=Rk(nlx,iq)
	Qraf(nl,iq)=Rk(nlx,iq)/R(nlx)
	Nrafk(nl,iq)=Rk(nlx,iq)/dble(mmol(gc(iq)))/masseU
	enddo
	Nraf(nl)=sum(Nrafk(nl,:))
	Praf(nl)=Nraf(nl)*kbolt*Traf(nl)
	Mraf(nl)=1D0/sum(Qraf(nl,:)/dble(mmol(gc(:))))
	END

	SUBROUTINE CORRMASS(qq,qint,FacMass,nl,nq)
	IMPLICIT NONE
	INTEGER :: nl,nq,l,nn
	REAL*8,DIMENSION(nl,nq) :: qq,qint,FacMass

	do nn=1,nq
	do l=1,nl
	FacMass(l,nn)=qq(l,nn)/qint(l,nn)
	enddo
	enddo

	END	


	SUBROUTINE DCOEFF(nn,Dij,P,T,N,Nk,D,nl,nq)
	IMPLICIT NONE
	REAL*8,DIMENSION(nl) :: N,T,D,P
	REAL*8,DIMENSION(nl,nq) :: Nk
	INTEGER :: nn,nl,nq,l,iq
	REAL,DIMENSION(nq,nq)  :: Dij
	REAL*8 :: interm,P0,T0,ptfac,dfac

	P0=1D5
	T0=273d0
	

	do l=1,nl
	ptfac=(P0/P(l))*(T(l)/T0)**1.75d0
	D(l)=0d0
	interm=0d0
	do iq=1,nq
	if (iq .ne. nn) then
	dfac=dble(dij(nn,iq))*ptfac
	interm=interm+Nk(l,iq)/N(l)/dfac
	endif
	enddo
!	D(l)=1d0/interm
        D(l)=(1D0-Nk(l,nn)/N(l))/interm
	enddo
	END
 
	SUBROUTINE HSCALEREAL(nn,Nk,Dz,H,nl,nq)
	IMPLICIT NONE
	INTEGER :: nn,nl,nq,l
	REAL*8,DIMENSION(nl) :: H
	REAL*8,DIMENSION(nl,nq) :: Nk
	REAL*8 :: Dz
	
	H(1)=(-3D0*Nk(1,nn)+4d0*NK(2,nn)-Nk(3,nn))/(2D0*DZ)/            &
     &  NK(1,nn)

	H(1)=-1D0/H(1)

	DO l=2,nl-1
	H(l)=(Nk(l+1,nn)-NK(l-1,nn))/(2D0*DZ)/NK(l,nn)
	H(l)=-1D0/H(l)
	ENDDO

	H(nl)=(3D0*Nk(nl,nn)-4D0*Nk(nl-1,nn)+Nk(nl-2,nn))/(2D0*DZ)/     &
     &  Nk(nl,nn) 
	H(nl)=-1D0/H(nl)

!	do l=1,nl
!	if (abs(H(l)) .lt. 100.) then
!	print*,'H',l,H(l),Nk(l,nn),nn
!	endif
!	enddo

	END
       
	SUBROUTINE VELVERT(nn,T,H,D,Dz,masse,W,nl)
	IMPLICIT NONE
	INTEGER :: l,nl,nn
	REAL*8,DIMENSION(nl) :: T,H,D,W,DT
	REAL*8 :: Dz,Hmol,masse
	REAL*8 :: kbolt,masseU,Konst,g
	masseU=1.660538782d-27
        kbolt=1.3806504d-23
        Konst=Kbolt/masseU
        g=3.72d0

	DT(1)=1D0/T(1)*(-3D0*T(1)+4D0*T(2)-T(3))/(2D0*DZ)
	Hmol=Konst*T(1)/masse/g
	W(1)=-D(1)*(1D0/H(1)-1D0/Hmol-DT(1))

	DO l=2,nl-1
	DT(l)=1D0/T(l)*(T(l+1)-T(l-1))/(2D0*DZ)
	Hmol=Konst*T(l)/masse/g
	W(l)=-D(l)*(1D0/H(l)-1D0/Hmol-DT(l))
	ENDDO

	DT(nl)=1D0/T(nl)*(3d0*T(nl)-4D0*T(nl-1)+T(nl-2))/(2D0*DZ)
	Hmol=Konst*T(nl)/masse/g
	W(nl)=-D(nl)*(1D0/H(nl)-1D0/Hmol-DT(nl))

!	do l=1,nl
!	print*,'W',W(l),D(l),H(l),DT(l)
!	enddo

	END

	SUBROUTINE TIMEDIFF(nn,H,W,TIME,nl)
	IMPLICIT NONE
	INTEGER :: nn,nl,l
	REAL*8,DIMENSION(nl) :: W,H,TIME

	DO l=1,nl
	TIME(l)=abs(H(l)/W(l))
	if (TIME(l) .lt. 1.D-4) then 
!	print*,'low tdiff',H(l),W(l),nn,l
	endif
	ENDDO

	END


	SUBROUTINE DIFFPARAM(T,D,dz,alpha,beta,gama,delta,eps,nl)
	IMPLICIT NONE
	INTEGER :: nl,l
	REAL*8,DIMENSION(nl) :: T,D
	REAL*8 :: DZ,DZinv
	REAL*8,DIMENSION(nl) :: alpha,beta,gama,delta,eps

! Compute alpha,beta and delta
! lower altitude values
	DZinv=1d0/(2D0*DZ)

	beta(1)=1d0/T(1)
	alpha(1)=beta(1)*(-3D0*T(1)+4D0*T(2)-T(3))*Dzinv
	delta(1)=(-3D0*D(1)+4D0*D(2)-D(3))*Dzinv

	beta(2)=1d0/T(2)
	alpha(2)=beta(2)*(T(3)-T(1))*Dzinv
	delta(2)=(D(3)-D(1))*Dzinv

!	do l=2,nl-1
!	beta(l)=1D0/T(l)
!	alpha(l)=beta(l)*(T(l+1)-T(l-1))*Dzinv
!	delta(l)=(D(l+1)-D(l-1))*Dzinv
!	end do   

	do l=3,nl-1
	beta(l)=1D0/T(l)
	alpha(l)=beta(l)*(T(l+1)-T(l-1))*Dzinv
	delta(l)=(D(l+1)-D(l-1))*Dzinv
	gama(l-1)=(beta(l)-beta(l-2))*Dzinv
	eps(l-1)=(alpha(l)-alpha(l-2))*Dzinv
	enddo

! Upper altitude values

	beta(nl)=1D0/T(nl)
	alpha(nl)=beta(nl)*(3D0*T(nl)-4D0*T(nl-1)+T(nl-2))*Dzinv     
	delta(nl)=(3D0*D(nl)-4D0*D(nl-1)+D(nl-2))*Dzinv

! Compute the gama and eps coefficients
! Lower altitude values

	gama(1)=(-3D0*beta(1)+4D0*beta(2)-beta(3))*Dzinv
	eps(1)=(-3D0*alpha(1)+4D0*alpha(2)-alpha(3))*Dzinv

	gama(nl-1)=(beta(nl)-beta(nl-2))*Dzinv
	eps(nl-1)=(alpha(nl)-alpha(nl-2))*Dzinv

!	do l=2,nl-1
!	gama(l)=(beta(l+1)-beta(l-1))*Dzinv
!	eps(l)=(alpha(l+1)-alpha(l-1))*Dzinv
!	end do

	gama(nl)=(3D0*beta(nl)-4D0*beta(nl-1)+beta(nl-2))*Dzinv
	eps(nl)=(3D0*alpha(nl)-4D0*alpha(nl-1)+alpha(nl-2))*Dzinv

!	do l=1,nl
!	print*,'test diffparam',alpha(l),beta(l),delta(l),gama(l),eps(l)
!	enddo
!	stop

	END


	SUBROUTINE MATCOEFF(alpha,beta,gama,delta,eps,Dad,rhoad,        &
     &  FacEsc,dz,dt,A,B,C,D,nl)
	IMPLICIT NONE
	INTEGER :: nl,l
	REAL*8, DIMENSION(nl) :: alpha,beta,gama,delta,eps,Dad,RHoad
	REAL*8 :: dz,dt,del1,del2,del3,FacEsc
	REAL*8, DIMENSION(nl) :: A,B,C,D
	del1=dt/2d0/dz
	del2=dt/dz/dz
	del3=dt

! lower boundary coefficients no change
	A(1)=0d0
	B(1)=1d0
	C(1)=0d0
	D(1)=rhoAd(1)

	do l=2,nl-1
	A(l)=(delta(l)+(alpha(l)+beta(l))*Dad(l))*del1-Dad(l)*del2
	B(l)=-(delta(l)*(alpha(l)+beta(l))+Dad(l)*(gama(l)+eps(l)))*del3
	B(l)=B(l)+1d0+2d0*Dad(l)*del2
	C(l)=-(delta(l)+(alpha(l)+beta(l))*Dad(l))*del1-Dad(l)*del2
	D(l)=rhoAD(l)
	enddo


! Upper boundary coefficients Diffusion profile
	C(nl)=0d0
	B(nl)=-1d0
	A(nl)=exp(-dZ*(alpha(nl)+beta(nl)))*FacEsc
	D(nl)=0D0


	END

	SUBROUTINE Checkmass(X,Y,nl,nn)
	IMPLICIT NONE

	INTEGER :: nl,nn
	REAL*8,DIMENSION(nl) :: X,Y
	REAL*8 Xtot,Ytot

	Xtot=sum(X)
	Ytot=sum(Y)

	if (abs((Xtot-Ytot)/Xtot) .gt. 1d-3) then
	print*,'no conservation for mass',Xtot,Ytot,nn
	endif
	END

	SUBROUTINE Checkmass2(qold,qnew,P,il,nl,nn,nq)
	IMPLICIT NONE
	INTEGER :: nl,nn,l,nq,il
	REAL,DIMENSION(nl+1) :: P
	REAL*8,DIMENSION(nl,nq) :: qold,qnew
	REAL*8 :: DM,Mold,Mnew,g
	g=3.72d0
	DM=0d0
	Mold=0d0
	Mnew=0d0
	DO l=il,nl
	DM=DM+(qnew(l,nn)-qold(l,nn))*(dble(P(l))-dble(P(l+1)))/g
	Mold=Mold+qold(l,nn)*(dble(P(l))-dble(P(l+1)))/g
	Mnew=Mnew+qnew(l,nn)*(dble(P(l))-dble(P(l+1)))/g
!	print*,'l',l,qold(l,nn),qnew(l,nn),Mold,Mnew,DM,P(l),P(l+1)
	ENDDO
	IF (abs(DM/Mold) .gt. 1d-2) THEN
	Print*,'We dont conserve mas',nn,DM,Mold,Mnew,DM/Mold
	ENDIF

	END

	SUBROUTINE GCMGRID_P(Z,P,Q,T,Nk,Rk,qq,qnew,tt,tnew,             &
     &    pp,M,gc,nl,nq,nlx,ig)
	use tracer_mod, only: nqmx
        IMPLICIT NONE
	INTEGER :: nl,nq,nlx,il,nn,iP,ig,compteur
	INTEGER,DIMENSION(1) :: indP
	INTEGER,DIMENSION(nq) :: gc 
	REAL*8,DIMENSION(nl) :: Z,P,T
	REAL*8,DIMENSION(nl,nq) :: Q,Nk,Rk
	REAL,DIMENSION(nqmx) :: M
	REAL*8,DIMENSION(nq) :: nNew
	REAL*8,DIMENSION(nlx) :: pp,tt,tnew
	REAL*8,DIMENSION(nlx) :: rhonew
	REAL*8,DIMENSION(nlx,nq) :: qq,qnew,rhoknew
	REAL*8 :: kbolt,masseU,Konst,g,Dz,facP,Hi
        REAL*8 :: Znew,Znew2,Pnew,Pnew2
        masseU=1.660538782d-27
        kbolt=1.3806504d-23
        Konst=Kbolt/masseU
        g=3.72d0
	Dz=Z(2)-Z(1)
	Znew=Z(nl)
	Znew2=Znew+dz
!	print*,'dz',Znew,Znew2,dz
	nNew(1:nq)=Nk(nl,1:nq)
	Pnew=P(nl)

	do il=1,nlx
!	print*,'il',il,pp(il),P(1),P(nl)
	if (pp(il) .ge. P(1)) then
	qnew(il,:)=qq(il,:)
	tnew(il)=tt(il)
	endif
	if (pp(il) .lt. P(1)) then
	if (pp(il) .gt. P(nl)) then
	indP=maxloc(P,mask=P <  pp(il))
        iP=indP(1)-1
	if (iP .lt. 1 .or. iP .gt. nl) then
        print*,'danger 2',iP,nl,pp(il)
        endif
	facP=(pp(il)-P(ip))/(P(ip+1)-P(ip))
!	print*,'P',P(ip),P(ip+1),facP,indP,iP

!	do nn=1,nq
!	qnew(il,nn)=Q(iP,nn)+
!     &  (Q(ip+1,nn)-Q(ip,nn))*facP
!	enddo

	do nn=1,nq
	rhoknew(il,nn)=Rk(iP,nn)+                                       &
     &  (Rk(ip+1,nn)-Rk(ip,nn))*facP
	enddo
	tnew(il)=T(iP)+(T(iP+1)-T(iP))*facP
	rhonew(il)=sum(rhoknew(il,:))
	do nn=1,nq
	qnew(il,nn)=rhoknew(il,nn)/rhonew(il)
	enddo

	else ! pp < P(nl) need to extrapolate density of each specie
	Pnew2=Pnew

	compteur=0
	do while (pnew2 .ge. pp(il))	
	compteur=compteur+1
	do nn=1,nq
	Hi=Konst*T(nl)/dble(M(gc(nn)))/g
	Nnew(nn)=Nnew(nn)*exp(-dZ/Hi)
	enddo
	Pnew=Pnew2
	Pnew2=kbolt*T(nl)*sum(Nnew(:))
	Znew=Znew2
	Znew2=Znew2+Dz
	if (compteur .ge. 100000) then
	print*,'error moldiff_red infinite loop'
	print*,ig,il,pp(il),tt(nl),pnew2,qnew(il,:),Znew2
	stop
	endif
!	print*,'test',Pnew2,Znew2,Nnew(nq),pp(il)
	enddo
	
	facP=(pp(il)-Pnew)/(Pnew2-Pnew)

!	do nn=1,nq
!	qnew(il,nn)=dble(M(gc(nn)))*Nnew(nn)
!     &  /sum(dble(M(gc(:)))*Nnew(:))
!	enddo

	do nn=1,nq
	rhoknew(il,nn)=dble(M(gc(nn)))*Nnew(nn)
	enddo
	rhonew(il)=sum(rhoknew(il,:))
	do nn=1,nq
	qnew(il,nn)=rhoknew(il,nn)/rhonew(il)
	enddo
	tnew(il)=T(nl)
	endif
	endif
	enddo

	END

	SUBROUTINE GCMGRID_P2(Z,P,Q,T,Nk,Rk,qq,qnew,tt,tnew             &
     &   ,pp,M,gc,nl,nq,nlx,facM,ig)
        use tracer_mod, only: nqmx
        IMPLICIT NONE
        INTEGER :: nl,nq,nlx,il,nn,iP,ig,compteur
        INTEGER,DIMENSION(1) :: indP
        INTEGER,DIMENSION(nq) :: gc
        REAL*8,DIMENSION(nl) :: Z,P,T
        REAL*8,DIMENSION(nl,nq) :: Q,Nk,Rk
        REAL,DIMENSION(nqmx) :: M
        REAL*8,DIMENSION(nq) :: nNew
        REAL*8,DIMENSION(nlx) :: pp,rhonew,tt,tnew
        REAL*8,DIMENSION(nlx,nq) :: qq,qnew,facM,rhoknew
        REAL*8 :: kbolt,masseU,Konst,g,Dz,facP,Hi
        REAL*8 :: Znew,Znew2,Pnew,Pnew2
        masseU=1.660538782d-27
        kbolt=1.3806504d-23
        Konst=Kbolt/masseU
        g=3.72d0
        Dz=Z(2)-Z(1)
        Znew=Z(nl)
        Znew2=Znew+dz
!       print*,'dz',Znew,Znew2,dz
        nNew(1:nq)=Nk(nl,1:nq)
        Pnew=P(nl)

        do il=1,nlx
!       print*,'il',il,pp(il),P(1),P(nl)
        if (pp(il) .ge. P(1)) then
        qnew(il,:)=qq(il,:)
	tnew(il)=tt(il)
        endif
        if (pp(il) .lt. P(1)) then
        if (pp(il) .gt. P(nl)) then
	indP=maxloc(P,mask=P <  pp(il))
        iP=indP(1)-1
	if (iP .lt. 1 .or. iP .gt. nl) then
        print*,'danger 3',iP,nl,pp(il)
        endif
        facP=(pp(il)-P(ip))/(P(ip+1)-P(ip))
!       print*,'P',P(ip),P(ip+1),facP,indP,iP

!        do nn=1,nq
!        qnew(il,nn)=Q(iP,nn)+
!     &  (Q(ip+1,nn)-Q(ip,nn))*facP
!        enddo

	do nn=1,nq
	rhoknew(il,nn)=(RK(iP,nn)+                                       &
     &  (RK(iP+1,nn)-Rk(iP,nn))*facP)*facM(il,nn)
	enddo
	tnew(il)=T(iP)+(T(ip+1)-T(iP))*facP
	rhonew(il)=sum(rhoknew(il,:))
	do nn=1,nq
	qnew(il,nn)=rhoknew(il,nn)/rhonew(il)
	enddo

        else ! pp < P(nl) need to extrapolate density of each specie
        Pnew2=Pnew

	compteur=0
        do while (pnew2 .ge. pp(il))    
	compteur=compteur+1
        do nn=1,nq
        Hi=Konst*T(nl)/dble(M(gc(nn)))/g
        Nnew(nn)=Nnew(nn)*exp(-dZ/Hi)
        enddo
        Pnew=Pnew2
        Pnew2=kbolt*T(nl)*sum(Nnew(:))
        Znew=Znew2
        Znew2=Znew2+Dz
	if (compteur .ge. 100000) then
	print*,'pb moldiff_red infinite loop'
	print*,ig,nl,T(nl),pnew2,qnew(il,:),Znew2
	stop
	endif

!       print*,'test',Pnew2,Znew2,Nnew(nq),pp(il)
        enddo
        
        facP=(pp(il)-Pnew)/(Pnew2-Pnew)

!        do nn=1,nq
!        qnew(il,nn)=dble(M(gc(nn)))*Nnew(nn)
!     &  /sum(dble(M(gc(:)))*Nnew(:))
!        enddo

	do nn=1,nq
        rhoknew(il,nn)=dble(M(gc(nn)))*Nnew(nn)*FacM(il,nn)
        enddo
        rhonew(il)=sum(rhoknew(il,:))
        do nn=1,nq
        qnew(il,nn)=rhoknew(il,nn)/rhonew(il)
        enddo
	tnew(il)=T(nl)

        endif
        endif
        enddo

        END
        
