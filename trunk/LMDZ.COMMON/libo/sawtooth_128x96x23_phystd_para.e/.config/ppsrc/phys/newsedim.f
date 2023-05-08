










      SUBROUTINE newsedim(ngrid,nlay,naersize,ptimestep,
     &  pplev,masse,epaisseur,pt,rd,rho,pqi,wq)
      
      use comcstfi_mod, only: r, g
      IMPLICIT NONE

!==================================================================
!     
!     Purpose
!     -------
!      Calculates sedimentation of 1 tracer 
!      of radius rd (m) and density rho (kg.m-3) 
!     
!==================================================================

!-----------------------------------------------------------------------
!   declarations
!   ------------

!   arguments
!   ---------

      integer,intent(in) :: ngrid  ! number of atmospheric columns
      integer,intent(in) :: nlay  ! number of atmospheric layers
      integer,intent(in) :: naersize   ! number of particle sizes (1 or number
                                       ! of grid boxes)
      real,intent(in) :: ptimestep ! physics time step (s)
      real,intent(in) :: pplev(ngrid,nlay+1)   ! inter-layer pressures (Pa)
      real,intent(in) :: pt(ngrid,nlay)    ! mid-layer temperatures (K)
      real,intent(in) :: masse (ngrid,nlay) ! atmospheric mass (kg)
      real,intent(in) :: epaisseur (ngrid,nlay)  ! thickness of atm. layers (m)
      real,intent(in) :: rd(naersize) ! particle radius (m)
      real,intent(in) :: rho ! particle density (kg.m-3)
      real,intent(inout) :: pqi(ngrid,nlay)  ! tracer   (e.g. ?/kg)
      real,intent(out) :: wq(ngrid,nlay+1)  ! flux of tracer during timestep (?/m-2)
      
c   local:
c   ------

      INTEGER l,ig, k, i
      REAL rfall

      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)

c    Traceurs :
c    ~~~~~~~~ 
      real traversee (ngrid,nlay)
      real vstokes(ngrid,nlay)
      real w(ngrid,nlay)
      real ptop, dztop, Ep, Stra


c    Physical constant
c    ~~~~~~~~~~~~~~~~~
c     Gas molecular viscosity (N.s.m-2)
      real,parameter :: visc=1.e-5       ! CO2
c     Effective gas molecular radius (m)
      real,parameter :: molrad=2.2e-10   ! CO2

c     local and saved variable
      real,save :: a,b
!$OMP THREADPRIVATE(a,b)

c    ** un petit test de coherence
c       --------------------------


      !print*,'temporary fixed particle rad in newsedim!!'

      IF (firstcall) THEN
         firstcall=.false.



!=======================================================================
!     Preliminary calculations for sedimenation velocity

!       - Constant to compute stokes speed simple formulae
!        (Vstokes =  b * rho* r**2   avec   b= (2/9) * rho * g / visc
         b = 2./9. * g / visc

!       - Constant  to compute gas mean free path
!        l= (T/P)*a, avec a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
         a = 0.707*8.31/(4*3.1416* molrad**2  * 6.023e23)

c       - Correction to account for non-spherical shape (Murphy et al.  1990)
c   (correction = 0.85 for irregular particles, 0.5 for disk shaped particles)
c        a = a    * 0.85


      ENDIF

c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c     Sedimentation velocity (m/s)
c     ~~~~~~~~~~~~~~~~~~~~~~
c     (stokes law corrected for low pressure by the Cunningham
c     slip-flow correction  according to Rossow (Icarus 36, 1-50, 1978)
       
        do  l=1,nlay
          do ig=1, ngrid
            if (naersize.eq.1) then 
              rfall=rd(1)
            else
              i=ngrid*(l-1)+ig
              rfall=rd(i) ! how can this be correct!!?
            endif  

            vstokes(ig,l) = b * rho * rfall**2 *
     &      (1 + 1.333* ( a*pt(ig,l)/pplev(ig,l) )/rfall)

c           Layer crossing time (s) :
            traversee(ig,l)= epaisseur(ig,l)/vstokes(ig,l)
          end do
        end do


c     Calcul de la masse d'atmosphere correspondant a q transferee
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     (e.g. on recherche le niveau  en dessous de laquelle le traceur
c      va traverser le niveau intercouche l : "dztop" est sa hauteur
c      au dessus de l (m), "ptop" est sa pression (Pa))
      do  l=1,nlay
        do ig=1, ngrid
             
             dztop = vstokes(ig,l)*  ptimestep 
             Ep=0
             k=0
           w(ig,l) = 0. !! JF+AS ajout initialisation (LK MARS)
c **************************************************************
c            Simple Method
cc             w(ig,l) =
cc     &       (1- exp(-dztop*g/(r*pt(ig,l))))*pplev(ig,l) / g
cc             write(*,*) 'OK simple method l,w =', l, w(ig,l)
cc            write(*,*) 'OK simple method dztop =', dztop
           w(ig,l) = 1. - exp(-dztop*g/(r*pt(ig,l)))
             !!! Diagnostic: JF. Fix: AS. Date: 05/11
             !!! Probleme arrondi avec la quantite ci-dessus
             !!! ---> vaut 0 pour -dztop*g/(r*pt(ig,l)) trop petit
             !!! ---> dans ce cas on utilise le developpement limite !
             !!! ---> exp(-x) = 1 - x lorsque x --> 0 avec une erreur de x^2 / 2  

             IF ( w(ig,l) .eq. 0. ) THEN
                w(ig,l) = ( dztop*g/(r*pt(ig,l)) ) * pplev(ig,l) / g
             ELSE
                w(ig,l) = w(ig,l) * pplev(ig,l) / g
             ENDIF
! LK borrowed simple method from Mars model (AS/JF)

!**************************************************************
cccc         Complex method :
           if (dztop.gt.epaisseur(ig,l)) then
cccc            Cas ou on "epuise" la couche l : On calcule le flux
cccc            Venant de dessus en tenant compte de la variation de Vstokes

               Ep= epaisseur(ig,l)
               Stra= traversee(ig,l)
               do while(dztop.gt.Ep.and.l+k+1.le.nlay)
                 k=k+1
                 dztop= Ep + vstokes(ig,l+k)*(ptimestep -Stra)
                 Ep = Ep + epaisseur(ig,l+k)
                 Stra = Stra + traversee(ig,l+k)
               enddo 
               Ep = Ep - epaisseur(ig,l+k)
!             ptop=pplev(ig,l+k)*exp(-(dztop-Ep)*g/(r*pt(ig,l+k)))
             ptop=exp(-(dztop-Ep)*g/(r*pt(ig,l+k)))
	     IF ( ptop .eq. 1. ) THEN
                !PRINT*, 'newsedim: exposant trop petit ', ig, l
                ptop=pplev(ig,l+k) * ( 1. - (dztop-Ep)*g/(r*pt(ig,l+k)))
             ELSE
                ptop=pplev(ig,l+k) * ptop
             ENDIF

             w(ig,l) = (pplev(ig,l) - ptop)/g

            endif   !!! complex method
c
cc           write(*,*) 'OK new    method l,w =', l, w(ig,l)
cc           write(*,*) 'OK new    method dztop =', dztop
cc       if(l.eq.7)write(*,*)'l=7,k,pplev,Ptop',pplev(ig,l),Ptop
cc       if(l.eq.7)write(*,*)'l=7,dztop,Ep',dztop,Ep
cc            if(l.eq.6)write(*,*)'l=6,k, w',k, w(1,l)
cc            if(l.eq.7)write(*,*)'l=7,k, w',k, w(1,l)
cc            if(l.eq.8)write(*,*)'l=8,k, w',k, w(1,l)
c **************************************************************

        end do
      end do

      call vlz_fi(ngrid,nlay,pqi,2.,masse,w,wq)
c         write(*,*) ' newsed: wq(6), wq(7), q(6)',
c    &                wq(1,6),wq(1,7),pqi(1,6)

      END
