!
! $Id: iniacademic.F90 1625 2012-05-09 13:14:48Z lguez $
!
SUBROUTINE iniacademic(vcov,ucov,teta,q,masse,ps,phis,time_0)

  USE filtreg_mod, ONLY: inifilr
  USE infotrac, ONLY : nqtot, ok_isotopes, iso_num, zone_num, &
                       iqpere, tnat, alpha_ideal, iso_indnum, &
                       phase_num, iqiso, ok_iso_verif
  USE control_mod, ONLY: day_step,planet_type
#ifdef CPP_IOIPSL
  USE IOIPSL, ONLY: getin
#else
  ! if not using IOIPSL, we still need to use (a local version of) getin
  USE ioipsl_getincom, ONLY: getin
#endif
  USE Write_Field
  use exner_hyb_m, only: exner_hyb
  use exner_milieu_m, only: exner_milieu
  USE comvert_mod, ONLY: ap,bp,preff,presnivs,pressure_exner
  USE comconst_mod, ONLY: im,jm,daysec,dtvr,kappa,cpp,g,pi
  USE logic_mod, ONLY: iflag_phys,read_start
  USE temps_mod, ONLY: annee_ref,day_ref,day_ini
  USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0

  !   Author:    Frederic Hourdin      original: 15/01/93
  ! The forcing defined here is from Held and Suarez, 1994, Bulletin
  ! of the American Meteorological Society, 75, 1825.

  IMPLICIT NONE

  !   Declararations:
  !   ---------------

  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "academic.h"
  include "iniprint.h"

  !   Arguments:
  !   ----------

  REAL,INTENT(OUT) :: time_0

  !   fields
  REAL,INTENT(OUT) :: vcov(ip1jm,llm) ! meridional covariant wind
  REAL,INTENT(OUT) :: ucov(ip1jmp1,llm) ! zonal covariant wind
  REAL,INTENT(OUT) :: teta(ip1jmp1,llm) ! potential temperature (K)
  REAL,INTENT(OUT) :: q(ip1jmp1,llm,nqtot) ! advected tracers (.../kg_of_air)
  REAL,INTENT(OUT) :: ps(ip1jmp1) ! surface pressure (Pa)
  REAL,INTENT(OUT) :: masse(ip1jmp1,llm) ! air mass in grid cell (kg)
  REAL,INTENT(OUT) :: phis(ip1jmp1) ! surface geopotential

  !   Local:
  !   ------

  REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
  REAL pks(ip1jmp1)                      ! exner au  sol
  REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
  REAL phi(ip1jmp1,llm)                  ! geopotentiel
  REAL ddsin,zsig,tetapv,w_pv  ! variables auxiliaires
  real tetastrat ! potential temperature in the stratosphere, in K
  real tetajl(jjp1,llm)
  INTEGER i,j,l,lsup,ij

  REAL teta0,ttp,delt_y,delt_z,eps ! Constantes pour profil de T
  REAL k_f,k_c_a,k_c_s         ! Constantes de rappel
  LOGICAL ok_geost             ! Initialisation vent geost. ou nul
  LOGICAL ok_pv                ! Polar Vortex
  REAL phi_pv,dphi_pv,gam_pv   ! Constantes pour polar vortex 

  real zz,ran1
  integer idum

  REAL alpha(ip1jmp1,llm),beta(ip1jmp1,llm),zdtvr
  
  character(len=*),parameter :: modname="iniacademic"
  character(len=80) :: abort_message


  ! Sanity check: verify that options selected by user are not incompatible
  if ((iflag_phys==1).and. .not. read_start) then
    write(lunout,*) trim(modname)," error: if read_start is set to ", &
    " false then iflag_phys should not be 1"
    write(lunout,*) "You most likely want an aquaplanet initialisation", &
    " (iflag_phys >= 100)"
    call abort_gcm(modname,"incompatible iflag_phys==1 and read_start==.false.",1)
  endif
  
  !-----------------------------------------------------------------------
  ! 1. Initializations for Earth-like case
  ! --------------------------------------
  !
  ! initialize planet radius, rotation rate,...
  call conf_planete

  time_0=0.
  day_ref=1
  annee_ref=0

  im         = iim
  jm         = jjm
  day_ini    = 1
  dtvr    = daysec/REAL(day_step)
  zdtvr=dtvr
  etot0      = 0.
  ptot0      = 0.
  ztot0      = 0.
  stot0      = 0.
  ang0       = 0.

  if (llm == 1) then
     ! specific initializations for the shallow water case
     kappa=1
  endif

  CALL iniconst
  CALL inigeom
  CALL inifilr

  if (llm == 1) then
     ! initialize fields for the shallow water case, if required
     if (.not.read_start) then
        phis(:)=0.
        q(:,:,:)=0
        CALL sw_case_williamson91_6(vcov,ucov,teta,masse,ps)
     endif
  endif

  academic_case: if (iflag_phys >= 2) then
     ! initializations

     ! 1. local parameters
     ! by convention, winter is in the southern hemisphere
     ! Geostrophic wind or no wind?
     ok_geost=.TRUE.
     CALL getin('ok_geost',ok_geost)
     ! Constants for Newtonian relaxation and friction
     k_f=1.                !friction 
     CALL getin('k_j',k_f)
     k_f=1./(daysec*k_f)
     k_c_s=4.  !cooling surface
     CALL getin('k_c_s',k_c_s)
     k_c_s=1./(daysec*k_c_s)
     k_c_a=40. !cooling free atm
     CALL getin('k_c_a',k_c_a)
     k_c_a=1./(daysec*k_c_a)
     ! Constants for Teta equilibrium profile
     teta0=315.     ! mean Teta (S.H. 315K)
     CALL getin('teta0',teta0)
     ttp=200.       ! Tropopause temperature (S.H. 200K)
     CALL getin('ttp',ttp)
     eps=0.         ! Deviation to N-S symmetry(~0-20K)
     CALL getin('eps',eps)
     delt_y=60.     ! Merid Temp. Gradient (S.H. 60K)
     CALL getin('delt_y',delt_y)
     delt_z=10.     ! Vertical Gradient (S.H. 10K)
     CALL getin('delt_z',delt_z)
     ! Polar vortex
     ok_pv=.false.
     CALL getin('ok_pv',ok_pv)
     phi_pv=-50.            ! Latitude of edge of vortex
     CALL getin('phi_pv',phi_pv)
     phi_pv=phi_pv*pi/180.
     dphi_pv=5.             ! Width of the edge
     CALL getin('dphi_pv',dphi_pv)
     dphi_pv=dphi_pv*pi/180.
     gam_pv=4.              ! -dT/dz vortex (in K/km)
     CALL getin('gam_pv',gam_pv)

     ! 2. Initialize fields towards which to relax
     ! Friction
     knewt_g=k_c_a
     DO l=1,llm
        zsig=presnivs(l)/preff
        knewt_t(l)=(k_c_s-k_c_a)*MAX(0.,(zsig-0.7)/0.3)
        kfrict(l)=k_f*MAX(0.,(zsig-0.7)/0.3)
     ENDDO
     DO j=1,jjp1
        clat4((j-1)*iip1+1:j*iip1)=cos(rlatu(j))**4
     ENDDO

     ! Potential temperature 
     DO l=1,llm
        zsig=presnivs(l)/preff
        tetastrat=ttp*zsig**(-kappa)
        tetapv=tetastrat
        IF ((ok_pv).AND.(zsig.LT.0.1)) THEN
           tetapv=tetastrat*(zsig*10.)**(kappa*cpp*gam_pv/1000./g)
        ENDIF
        DO j=1,jjp1
           ! Troposphere
           ddsin=sin(rlatu(j))
           tetajl(j,l)=teta0-delt_y*ddsin*ddsin+eps*ddsin &
                -delt_z*(1.-ddsin*ddsin)*log(zsig)
           if (planet_type=="giant") then
             tetajl(j,l)=teta0+(delt_y*                   &
                ((sin(rlatu(j)*3.14159*eps+0.0001))**2)   &
                / ((rlatu(j)*3.14159*eps+0.0001)**2))     &
                -delt_z*log(zsig)
           endif
           ! Profil stratospherique isotherme (+vortex)
           w_pv=(1.-tanh((rlatu(j)-phi_pv)/dphi_pv))/2.
           tetastrat=tetastrat*(1.-w_pv)+tetapv*w_pv             
           tetajl(j,l)=MAX(tetajl(j,l),tetastrat)  
        ENDDO
     ENDDO

     !          CALL writefield('theta_eq',tetajl)

     do l=1,llm
        do j=1,jjp1
           do i=1,iip1
              ij=(j-1)*iip1+i
              tetarappel(ij,l)=tetajl(j,l)
           enddo
        enddo
     enddo

     ! 3. Initialize fields (if necessary)
     IF (.NOT. read_start) THEN
        ! surface pressure
        if (iflag_phys>2) then
           ! specific value for CMIP5 aqua/terra planets
           ! "Specify the initial dry mass to be equivalent to
           !  a global mean surface pressure (101325 minus 245) Pa."
           ps(:)=101080.  
        else
           ! use reference surface pressure
           ps(:)=preff
        endif
        
        ! ground geopotential
        phis(:)=0.

        CALL pression ( ip1jmp1, ap, bp, ps, p       )
        if (pressure_exner) then
          CALL exner_hyb( ip1jmp1, ps, p, pks, pk)
        else
          call exner_milieu(ip1jmp1,ps,p,pks,pk)
        endif
        CALL massdair(p,masse)

        ! bulk initialization of temperature
        teta(:,:)=tetarappel(:,:)

        ! geopotential
        CALL geopot(ip1jmp1,teta,pk,pks,phis,phi)

        ! winds
        if (ok_geost) then
           call ugeostr(phi,ucov)
        else
           ucov(:,:)=0.
        endif
        vcov(:,:)=0.

        ! bulk initialization of tracers
        if (planet_type=="earth") then
           ! Earth: first two tracers will be water
           do i=1,nqtot
              if (i == 1) q(:,:,i)=1.e-10
              if (i == 2) q(:,:,i)=1.e-15
              if (i.gt.2) q(:,:,i)=0.

              ! CRisi: init des isotopes
              ! distill de Rayleigh très simplifiée
              if (ok_isotopes) then
                if ((iso_num(i).gt.0).and.(zone_num(i).eq.0)) then          
                   q(:,:,i)=q(:,:,iqpere(i))             &
      &                  *tnat(iso_num(i))               &
      &                  *(q(:,:,iqpere(i))/30.e-3)      &
      &                  **(alpha_ideal(iso_num(i))-1)
                endif                
                if ((iso_num(i).gt.0).and.(zone_num(i).eq.1)) then
                  q(:,:,i)=q(:,:,iqiso(iso_indnum(i),phase_num(i)))
                endif
              endif !if (ok_isotopes) then

           enddo
        else
           q(:,:,:)=0
        endif ! of if (planet_type=="earth")

        if (ok_iso_verif) then
! Ehouarn: this will onyly work in serial mode
!           call check_isotopes_seq(q,1,ip1jmp1,'iniacademic_loc')
        endif !if (ok_iso_verif) then

        ! add random perturbation to temperature
        idum  = -1
        zz = ran1(idum)
        idum  = 0
        do l=1,llm
           do ij=iip2,ip1jm
              teta(ij,l)=teta(ij,l)*(1.+0.005*ran1(idum))
           enddo
        enddo

        ! maintain periodicity in longitude
        do l=1,llm
           do ij=1,ip1jmp1,iip1
              teta(ij+iim,l)=teta(ij,l)
           enddo
        enddo

     ENDIF ! of IF (.NOT. read_start)
  endif academic_case

END SUBROUTINE iniacademic
