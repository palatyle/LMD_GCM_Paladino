!*****************************************************************
!
!     Photochemical routine 
!
!     Author: Franck Lefevre
!     ------
!
!     Version: 10/11/2015
!
!*****************************************************************

subroutine photochemistry_asis(nlayer, nq, ngrid,                               &
                          ig, lswitch, zycol, sza, fractcol, ptimestep, press,  &
                          temp, dens, zmmean, dist_sol, surfdust1d,             &
                          surfice1d, jo3, tau, iter)

      use callkeys_mod
implicit none

#include "chimiedata.h" 

!===================================================================
!     inputs:
!===================================================================

integer, intent(in) :: nlayer ! number of atmospheric layers
integer, intent(in) :: nq     ! number of tracers
integer,intent(in) :: ngrid   ! number of atmospheric columns
integer :: ig                 ! grid point index
      
real :: sza                   ! solar zenith angle (deg)
real :: fractcol              ! day fraction
real :: ptimestep             ! physics timestep (s)
real :: press(nlayer)         ! pressure (hpa)
real :: temp(nlayer)          ! temperature (k)
real :: dens(nlayer)          ! density (cm-3)
real :: zmmean(nlayer)        ! mean molar mass (g/mole)
real :: dist_sol              ! sun distance (au) 
real :: surfdust1d(nlayer)    ! dust surface area (cm2/cm3)
real :: surfice1d(nlayer)     ! ice surface area (cm2/cm3)
real :: tau                   ! optical depth at 7 hpa

!===================================================================
!     input/output:
!===================================================================
      
real :: zycol(nlayer,nq)       ! chemical species volume mixing ratio

!===================================================================
!     output:
!===================================================================
      
integer :: iter(nlayer)        ! iteration counter
real    :: jo3(nlayer)         ! photodissociation rate o3 -> o1d

!===================================================================
!     local:
!===================================================================

integer :: phychemrat         ! (physical timestep)/(nominal chemical timestep)
integer :: j_o3_o1d, ilev
integer :: iesp, nesp
integer :: lswitch

logical, save :: firstcall = .true.

parameter (nesp = 17)         ! number of species in the chemistry

! tracer indexes in the chemistry:

integer,parameter :: i_co2  =  1
integer,parameter :: i_co   =  2
integer,parameter :: i_o    =  3
integer,parameter :: i_o1d  =  4
integer,parameter :: i_o2   =  5
integer,parameter :: i_o3   =  6
integer,parameter :: i_h    =  7
integer,parameter :: i_h2   =  8
integer,parameter :: i_oh   =  9
integer,parameter :: i_ho2  = 10
integer,parameter :: i_h2o2 = 11
integer,parameter :: i_h2o  = 12
integer,parameter :: i_n    = 13
integer,parameter :: i_n2d  = 14
integer,parameter :: i_no   = 15
integer,parameter :: i_no2  = 16
integer,parameter :: i_n2   = 17

real :: stimestep           ! standard timestep for the chemistry (s) 
real :: ctimestep           ! real timestep for the chemistry (s) 
real :: dt_guess            ! first-guess timestep (s) 
real :: dt_corrected        ! corrected timestep (s) 
real :: dt_min = 1.         ! minimum allowed timestep (s) 
real :: dtg                 ! correction factor for the timestep (s) 
real :: j(nlayer,nd)        ! interpolated photolysis rates (s-1)
real :: time                ! internal time (between 0 and ptimestep, in s)

real, dimension(nlayer,nesp)            :: rm   ! mixing ratios 
real (kind = 8), dimension(nesp)        :: cold ! number densities at previous timestep (molecule.cm-3) 
real (kind = 8), dimension(nlayer,nesp) :: c    ! number densities at current timestep (molecule.cm-3) 
real (kind = 8), dimension(nesp)        :: cnew ! number densities at next timestep (molecule.cm-3) 
 
! reaction rates
 
real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4
logical :: hetero_dust, hetero_ice

! matrix

real (kind = 8), dimension(nesp,nesp) :: mat, mat1
integer, dimension(nesp)              :: indx
integer                               :: code

! production and loss terms (for first-guess solution only)

real (kind = 8), dimension(nesp) :: prod, loss

! curvatures

real :: ratio, curv, e, e1, e2, e3

!===================================================================
!     initialisation of the reaction indexes
!===================================================================

if (firstcall) then
   print*,'photochemistry: initialize indexes'
   call indice(i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
               i_n, i_n2d, i_no, i_no2, i_n2)
   firstcall = .false.
end if

!===================================================================
!     initialisation of mixing ratios and densities       
!===================================================================

call gcmtochim(nlayer, nq, zycol, lswitch, nesp,               &
               i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
               i_n, i_n2d, i_no, i_no2, i_n2, dens, rm, c) 

!===================================================================
!     interpolation of photolysis rates in the lookup table      
!===================================================================

call photolysis_asis(nlayer, ngrid, lswitch, press, temp, sza, fractcol, tau, zmmean, dist_sol, & 
                     rm(:,i_co2), rm(:,i_o3), v_phot)

! save o3 photolysis for output

j_o3_o1d = 5
jo3(:) = v_phot(:,j_o3_o1d)

!===================================================================
!     reaction rates                                     
!===================================================================
!     switches for heterogeneous chemistry
!     hetero_ice  : reactions on ice clouds
!     hetero_dust : reactions on dust    
!===================================================================

hetero_dust = .false.
hetero_ice  = .false.

call reactionrates(nlayer, lswitch, dens, c(:,i_co2), c(:,i_o2), &
                   press, temp, hetero_dust, hetero_ice,         &
                   surfdust1d, surfice1d, v_phot, v_3, v_4)

!===================================================================
!     stimestep : standard chemical timestep (s)                           
!     ctimestep : real chemical timestep (s),
!                 taking into account the physical timestep                           
!===================================================================

stimestep = 600. ! standard value : 10 mn

phychemrat = nint(ptimestep/stimestep)
phychemrat = 1

ctimestep = ptimestep/real(phychemrat)

!print*, "stimestep  = ", stimestep
!print*, "ptimestep  = ", ptimestep
!print*, "phychemrat = ", phychemrat
!print*, "ctimestep  = ", ctimestep
!stop

!===================================================================
!     loop over levels         
!===================================================================

do ilev = 1,lswitch - 1

!  initializations

   time = 0.
   iter(ilev) = 0
   dt_guess = ctimestep
   cold(:) = c(ilev,:)

!  internal loop for the chemistry

   do while (time < ptimestep)

   iter(ilev) = iter(ilev) + 1    ! iteration counter

!  first-guess: fill matrix

   call fill_matrix(ilev, mat1, prod, loss, c, nesp, nlayer, v_phot, v_3, v_4)

!  adaptative evaluation of the sub time step

   call define_dt(nesp, dt_corrected, dt_guess, ctimestep, cold(:), c(ilev,:),  &
                  mat1, prod, loss, dens(ilev))

   if (time + dt_corrected > ptimestep) then
      dt_corrected = ptimestep - time
   end if

!  if (dt_corrected /= dt_guess) then  ! the timestep has been modified

!  form the matrix identity + mat*dt_corrected

   mat(:,:) = mat1(:,:)*dt_corrected
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

!  solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

   cnew(:) = c(ilev,:)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry_asis error, missing LAPACK routine dgesv"
   stop
#endif

!  end if

!  eliminate small values

   where (cnew(:)/dens(ilev) < 1.e-30)
      cnew(:) = 0.
   end where

!  update concentrations

   cold(:)   = c(ilev,:)
   c(ilev,:) = cnew(:)
   cnew(:)   = 0.

!  increment internal time

   time = time + dt_corrected
   dt_guess = dt_corrected     ! first-guess timestep for next iteration

   end do ! while (time < ptimestep)

end do ! ilev

!===================================================================
!     save chemical species for the gcm       
!===================================================================

call chimtogcm(nlayer, nq, zycol, lswitch, nesp,              &
               i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,      &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,              &
               i_n, i_n2d, i_no, i_no2, i_n2, dens, c) 

contains

!================================================================

 subroutine define_dt(nesp, dtnew, dtold, ctimestep, cold, ccur, mat1, &
                      prod, loss, dens)

!================================================================
! iterative evaluation of the appropriate time step dtnew
! according to curvature criterion based on
! e = 2 Rtol [r Cn+1 -(1-r) Cn + Cn-1 ]/[(1+r) Cn]
! with r = (tn - tn-1)/(tn+1 - tn)
!================================================================

implicit none

! input

integer :: nesp  ! number of species in the chemistry

real :: dtold, ctimestep
real (kind = 8), dimension(nesp)      :: cold, ccur
real (kind = 8), dimension(nesp,nesp) :: mat1
real (kind = 8), dimension(nesp)      :: prod, loss
real                                  :: dens

! output

real :: dtnew

! local

real (kind = 8), dimension(nesp)      :: cnew
real (kind = 8), dimension(nesp,nesp) :: mat
real (kind = 8) :: atol, ratio, e, es, coef

integer                  :: code, iesp, iter
integer, dimension(nesp) :: indx

real :: dttest

! parameters

real (kind = 8), parameter :: dtmin   = 10.      ! minimum time step (s)
real (kind = 8), parameter :: vmrtol  = 1.e-11   ! absolute tolerance on vmr
real (kind = 8), parameter :: rtol    = 1./0.05   ! 1/rtol recommended value : 0.1-0.02
integer,         parameter :: niter   = 3        ! number of iterations
real (kind = 8), parameter :: coefmax = 2.
real (kind = 8), parameter :: coefmin = 0.1
logical                    :: fast_guess = .true.

dttest = dtold   ! dttest = dtold = dt_guess

atol = vmrtol*dens ! absolute tolerance in molecule.cm-3

do iter = 1,niter

if (fast_guess) then

! first guess : fast semi-implicit method

   do iesp = 1, nesp
      cnew(iesp) = (ccur(iesp) + prod(iesp)*dttest)/(1. + loss(iesp)*dttest)
   end do

else

! first guess : form the matrix identity + mat*dt_guess

   mat(:,:) = mat1(:,:)*dttest
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

! form right-hand side (RHS) of the system

   cnew(:) = ccur(:)

! solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

#ifdef LAPACK
      call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry_asis error, missing LAPACK routine dgesv"
   stop
#endif

end if

! ratio old/new subtimestep

ratio = dtold/dttest

! e : local error indicatocitr

e = 0.

do iesp = 1,nesp
   es = 2.*abs((ratio*cnew(iesp) - (1. + ratio)*ccur(iesp) + cold(iesp))   &
         /(1. + ratio)/max(ccur(iesp),atol))

   if (es > e) then
      e = es
   end if
end do
e = rtol*e

! timestep correction

coef = max(coefmin, min(coefmax,0.8/sqrt(e)))

dttest = max(dtmin,dttest*coef)
dttest = min(ctimestep,dttest)

end do ! iter

! new timestep

dtnew = dttest

end subroutine define_dt




!======================================================================

 subroutine reactionrates(nlayer,                               &
                          lswitch, dens, co2, o2, press, t,     &
                          hetero_dust, hetero_ice,              &
                          surfdust1d, surfice1d,                &
                          v_phot, v_3, v_4)
 
!================================================================
! compute reaction rates                                        !
!----------------------------------------------------------------
! reaction               type                array              !
!----------------------------------------------------------------
! A + B    --> C + D     bimolecular         v_4                !
! A + A    --> B + C     quadratic           v_3                !
! A + C    --> B + C     quenching           v_phot             !
! A + ice  --> B + C     heterogeneous       v_phot             !
!================================================================

use comcstfi_mod

implicit none

#include "chimiedata.h"

!----------------------------------------------------------------------
!     input
!----------------------------------------------------------------------

integer, intent(in)     :: nlayer            ! number of atmospheric layers
integer                 :: lswitch           ! interface level between lower
                                             ! atmosphere and thermosphere chemistries
real, dimension(nlayer) :: dens              ! total number density (molecule.cm-3)
real, dimension(nlayer) :: press             ! pressure (hPa)
real, dimension(nlayer) :: t                 ! temperature (K)
real, dimension(nlayer) :: surfdust1d        ! dust surface area (cm2.cm-3)
real, dimension(nlayer) :: surfice1d         ! ice surface area (cm2.cm-3)
real (kind = 8), dimension(nlayer) :: co2    ! co2 number density (molecule.cm-3)
real (kind = 8), dimension(nlayer) :: o2     ! o2 number density (molecule.cm-3)
logical :: hetero_dust, hetero_ice           ! switches for heterogeneous chemistry

!----------------------------------------------------------------------
!     output
!----------------------------------------------------------------------

real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4

!----------------------------------------------------------------------
!     local
!----------------------------------------------------------------------

integer          :: ilev
integer          :: nb_phot, nb_reaction_3, nb_reaction_4
real :: ak0, ak1, xpo, rate
real :: k1a0, k1b0, k1ainf, k1a, k1b, fc, fx, x, y, gam
real, dimension(nlayer) :: deq
real, dimension(nlayer) :: a001, a002, a003,                           &
                             b001, b002, b003, b004, b005, b006, b007,   &
                             b008, b009,                                 &
                             c001, c002, c003, c004, c005, c006, c007,   &
                             c008, c009, c010, c011, c012, c013, c014,   &
                             c015, c016, c017, c018,                     &
                             d001, d002, d003, d004, d005, d006, d007,   &
                             d008, d009,                                 &
                             e001, e002,                                 &
                             h001, h002, h003, h004, h005

!----------------------------------------------------------------------
!     initialisation
!----------------------------------------------------------------------

      nb_phot       = 11       ! jmars.20140930 reduit de 13 a 11
      nb_reaction_3 = 0
      nb_reaction_4 = 0

!----------------------------------------------------------------------
!     oxygen reactions 
!----------------------------------------------------------------------

!---  a001: o + o2 + co2 -> o3 + co2

!     jpl 2003
!
!     co2 efficiency as a third body (2.075)
!     from sehested et al., j. geophys. res., 100, 1995.

      a001(:) = 2.075*6.0e-34*(t(:)/300.)**(-2.4)*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a001(:)

!---  a002: o + o + co2 -> o2 + co2

!     Tsang and Hampson, J. Chem. Phys. Ref. Data, 15, 1087, 1986

!     a002(:) = 2.5*5.2e-35*exp(900./t(:))*dens(:)

!     Campbell and Gray, Chem. Phys. Lett., 18, 607, 1973

!     a002(:) = 1.2e-32*(300./t(:))**(2.0)*dens(:)  ! yung expression
      a002(:) = 2.5*9.46e-34*exp(485./t(:))*dens(:) ! nist expression

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = a002(:)

!---  a003: o + o3 -> o2 + o2

!     jpl 2003

      a003(:) = 8.0e-12*exp(-2060./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a003(:)

!----------------------------------------------------------------------
!     o(1d) reactions
!----------------------------------------------------------------------

!---  b001: o(1d) + co2  -> o + co2

!     jpl 2006

      b001(:) = 7.5e-11*exp(115./t(:))
   
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b001(:)*co2(:)

!---  b002: o(1d) + h2o  -> oh + oh

!     jpl 2006
 
      b002(:) = 1.63e-10*exp(60./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b002(:)

!---  b003: o(1d) + h2  -> oh + h

!     jpl 2011

      b003(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b003(:)

!---  b004: o(1d) + o2  -> o + o2

!     jpl 2006

      b004(:) = 3.3e-11*exp(55./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b004(:)*o2(:)
    
!---  b005: o(1d) + o3  -> o2 + o2

!     jpl 2003

      b005(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b005(:)
    
!---  b006: o(1d) + o3  -> o2 + o + o

!     jpl 2003

      b006(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b006(:)
    
!---  b007: o(1d) + ch4 -> ch3 + oh

!     jpl 2003

      b007(:) = 1.5e-10*0.75

!---  b008: o(1d) + ch4 -> ch3o + h

!     jpl 2003

      b008(:) = 1.5e-10*0.20
!
!---  b009: o(1d) + ch4 -> ch2o + h2

!     jpl 2003

      b009(:) = 1.5e-10*0.05

!----------------------------------------------------------------------
!     hydrogen reactions
!----------------------------------------------------------------------

!---  c001: o + ho2 -> oh + o2

!     jpl 2003

      c001(:) = 3.0e-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c001(:)

!---  c002: o + oh -> o2 + h

!     jpl 2011

      c002(:) = 1.8e-11*exp(180./t(:))

!     robertson and smith, j. chem. phys. a 110, 6673, 2006

!     c002(:) = 11.2e-11*t(:)**(-0.32)*exp(177./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c002(:)

!---  c003: h + o3 -> oh + o2

!     jpl 2003

      c003(:) = 1.4e-10*exp(-470./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c003(:)

!---  c004: h + ho2 -> oh + oh

!     jpl 2006

      c004(:) = 7.2e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c004(:)

!---  c005: h + ho2 -> h2 + o2

!     jpl 2006

      c005(:) = 6.9e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c005(:)

!---  c006: h + ho2 -> h2o + o

!     jpl 2006

      c006(:) = 1.6e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c006(:)

!---  c007: oh + ho2 -> h2o + o2

!     jpl 2003

!     canty et al., grl, 2006 suggest to increase this rate
!     by 20%. not done here.

      c007(:) = 4.8e-11*exp(250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c007(:)

!---  c008: ho2 + ho2 -> h2o2 + o2

!     jpl 2006

!     c008(:) = 3.5e-13*exp(430./t(:))

!     christensen et al., grl, 13, 2002

      c008(:) = 1.5e-12*exp(19./t(:))

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c008(:)

!---  c009: oh + h2o2 -> h2o + ho2

!     jpl 2006

      c009(:) = 1.8e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c009(:)

!---  c010: oh + h2 -> h2o + h

!     jpl 2006

      c010(:) = 2.8e-12*exp(-1800./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c010(:)

!---  c011: h + o2 + co2 -> ho2 + co2

!     jpl 2011

      do ilev = 1,lswitch-1
         ak0 = 2.5*4.4e-32*(t(ilev)/300.)**(-1.3)
         ak1 = 7.5e-11*(t(ilev)/300.)**(0.2)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2)
         c011(ilev) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c011(:)

!---  c012: o + h2o2 -> oh + ho2

!     jpl 2003

      c012(:) = 1.4e-12*exp(-2000./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c012(:)

!---  c013: oh + oh -> h2o + o

!     jpl 2006

      c013(:) = 1.8e-12

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c013(:)

!---  c014: oh + o3 -> ho2 + o2

!     jpl 2003

      c014(:) = 1.7e-12*exp(-940./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c014(:)

!---  c015: ho2 + o3 -> oh + o2 + o2

!     jpl 2003

      c015(:) = 1.0e-14*exp(-490./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c015(:)

!---  c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2

!     jpl 2011

      c016(:) = 2.5*2.1e-33*exp(920./t(:))*dens(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c016(:)

!---  c017: oh + oh + co2 -> h2o2 + co2

!     jpl 2003

      do ilev = 1,lswitch-1
         ak0 = 2.5*6.9e-31*(t(ilev)/300.)**(-1.0)
         ak1 = 2.6e-11*(t(ilev)/300.)**(0.0)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2)
         c017(ilev) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c017(:)

!---  c018: h + h + co2 -> h2 + co2

!     baulch et al., 2005

      c018(:) = 2.5*1.8e-30*(t(:)**(-1.0))*dens(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c018(:)

!----------------------------------------------------------------------
!     nitrogen reactions
!----------------------------------------------------------------------

!---  d001: no2 + o -> no + o2

!     jpl 2006

      d001(:) = 5.1e-12*exp(210./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d001(:)

!---  d002: no + o3 -> no2 + o2

!     jpl 2006

      d002(:) = 3.0e-12*exp(-1500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d002(:)

!---  d003: no + ho2 -> no2 + oh

!     jpl 2011

      d003(:) = 3.3e-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d003(:)

!---  d004: n + no -> n2 + o

!     jpl 2011

      d004(:) = 2.1e-11*exp(100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d004(:)

!---  d005: n + o2 -> no + o

!     jpl 2011

      d005(:) = 1.5e-11*exp(-3600./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d005(:)

!---  d006: no2 + h -> no + oh

!     jpl 2011

      d006(:) = 4.0e-10*exp(-340./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d006(:)

!---  d007: n + o -> no
 
      d007(:) = 2.8e-17*(300./t(:))**0.5

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d007(:)

!---  d008: n + ho2 -> no + oh

!     brune et al., j. chem. phys., 87, 1983 

      d008(:) = 2.19e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d008(:)

!---  d009: n + oh -> no + h

!     atkinson et al., j. phys. chem. ref. data, 18, 881, 1989 

      d009(:) = 3.8e-11*exp(85./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d009(:)

!----------------------------------------------------------------------
!     carbon reactions
!----------------------------------------------------------------------

!---  e001: oh + co -> co2 + h

!     jpl 2003

!     e001(:) = 1.5e-13*(1 + 0.6*press(:)/1013.)

!     mccabe et al., grl, 28, 3135, 2001

!     e001(:) = 1.57e-13 + 3.54e-33*dens(:)

!     jpl 2006

!     ak0 = 1.5e-13*(t(:)/300.)**(0.6)
!     ak1 = 2.1e-9*(t(:)/300.)**(6.1)
!     rate1 = ak0/(1. + ak0/(ak1/dens(:)))
!     xpo1 = 1./(1. + alog10(ak0/(ak1/dens(:)))**2)

!     ak0 = 5.9e-33*(t(:)/300.)**(-1.4)
!     ak1 = 1.1e-12*(t(:)/300.)**(1.3)
!     rate2 = (ak0*dens(:))/(1. + ak0*dens(:)/ak1)
!     xpo2 = 1./(1. + alog10((ak0*dens(:))/ak1)**2)

!     e001(:) = rate1*0.6**xpo1 + rate2*0.6**xpo2

!     joshi et al., 2006

      do ilev = 1,lswitch-1
         k1a0 = 1.34*2.5*dens(ilev)                                  &
               *1/(1/(3.62e-26*t(ilev)**(-2.739)*exp(-20./t(ilev)))  &
               + 1/(6.48e-33*t(ilev)**(0.14)*exp(-57./t(ilev))))     ! typo in paper corrected
         k1b0 = 1.17e-19*t(ilev)**(2.053)*exp(139./t(ilev))          &
              + 9.56e-12*t(ilev)**(-0.664)*exp(-167./t(ilev))
         k1ainf = 1.52e-17*t(ilev)**(1.858)*exp(28.8/t(ilev))        &
                + 4.78e-8*t(ilev)**(-1.851)*exp(-318./t(ilev))
         x = k1a0/(k1ainf - k1b0)
         y = k1b0/(k1ainf - k1b0)
         fc = 0.628*exp(-1223./t(ilev)) + (1. - 0.628)*exp(-39./t(ilev))  &
            + exp(-t(ilev)/255.)
         fx = fc**(1./(1. + (alog(x))**2))                           ! typo in paper corrected
         k1a = k1a0*((1. + y)/(1. + x))*fx
         k1b = k1b0*(1./(1.+x))*fx

         e001(ilev) = k1a + k1b
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e001(:)

!---  e002: o + co + m -> co2 + m

!     tsang and hampson, 1986.

      e002(:) = 2.5*6.5e-33*exp(-2184./t(:))*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e002(:)

!----------------------------------------------------------------------
!     heterogeneous chemistry 
!----------------------------------------------------------------------

      if (hetero_ice) then

!        k = (surface*v*gamma)/4 (s-1)
!        v = 100*sqrt(8rt/(pi*m))  (cm s-1)
 
!---     h001: ho2 + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.025
      
         gam = 0.025
         h001(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(33.e-3*pi))*gam/4.
 
!        h002: oh + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.03
 
         gam = 0.03
         h002(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(17.e-3*pi))*gam/4.

!---     h003: h2o2 + ice -> products
 
!        gamma = 0.    test value
 
         gam = 0.
         h003(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(34.e-3*pi))*gam/4.
      else
         h001(:) = 0.
         h002(:) = 0.
         h003(:) = 0.
      end if

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h001(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h002(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h003(:)

      if (hetero_dust) then
 
!---     h004: ho2 + dust -> products
 
!        jacob, 2000: gamma = 0.2
!        see dereus et al., atm. chem. phys., 2005
 
         gam = 0.2
         h004(:) = surfdust1d(:)  &
                   *100.*sqrt(8.*8.31*t(:)/(33.e-3*pi))*gam/4.
 
!---     h005: h2o2 + dust -> products
 
!        gamma = 5.e-4
!        see dereus et al., atm. chem. phys., 2005
 
         gam = 5.e-4
         h005(:) = surfdust1d(:)  &
                   *100.*sqrt(8.*8.31*t(:)/(34.e-3*pi))*gam/4.
      else
         h004(:) = 0.
         h005(:) = 0.
      end if
 
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h004(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h005(:)

end subroutine reactionrates

!======================================================================

 subroutine fill_matrix(ilev, mat, prod, loss, c, nesp, nlayer, v_phot, v_3, v_4)

!======================================================================
! filling of the jacobian matrix
!======================================================================

use types_asis

implicit none

#include "chimiedata.h"

! input

integer             :: ilev    ! level index
integer             :: nesp    ! number of species in the chemistry
integer, intent(in) :: nlayer  ! number of atmospheric layers

real (kind = 8), dimension(nlayer,nesp)              :: c    ! number densities
real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4

! output

real (kind = 8), dimension(nesp,nesp), intent(out) :: mat  ! matrix
real (kind = 8), dimension(nesp), intent(out)      :: prod, loss

! local

integer :: iesp
integer :: ind_phot_2,ind_phot_4,ind_phot_6
integer :: ind_3_2,ind_3_4,ind_3_6
integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
integer :: iphot,i3,i4

real(kind = jprb) :: eps, eps_4  ! implicit/explicit coefficient

! initialisations 

mat(:,:) = 0.
prod(:)  = 0.
loss(:)  = 0.

! photodissociations
! or reactions a + c -> b + c
! or reactions a + ice -> b + c

do iphot = 1,nb_phot_max

  ind_phot_2 = indice_phot(iphot)%z2
  ind_phot_4 = indice_phot(iphot)%z4
  ind_phot_6 = indice_phot(iphot)%z6

  mat(ind_phot_2,ind_phot_2) = mat(ind_phot_2,ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  mat(ind_phot_4,ind_phot_2) = mat(ind_phot_4,ind_phot_2) - indice_phot(iphot)%z3*v_phot(ilev,iphot)
  mat(ind_phot_6,ind_phot_2) = mat(ind_phot_6,ind_phot_2) - indice_phot(iphot)%z5*v_phot(ilev,iphot)

  loss(ind_phot_2) = loss(ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  prod(ind_phot_4) = prod(ind_phot_4) + indice_phot(iphot)%z3*v_phot(ilev,iphot)*c(ilev,ind_phot_2)
  prod(ind_phot_6) = prod(ind_phot_6) + indice_phot(iphot)%z5*v_phot(ilev,iphot)*c(ilev,ind_phot_2)

end do

! reactions a + a -> b + c 

do i3 = 1,nb_reaction_3_max

  ind_3_2 = indice_3(i3)%z2
  ind_3_4 = indice_3(i3)%z4
  ind_3_6 = indice_3(i3)%z6

  mat(ind_3_2,ind_3_2) = mat(ind_3_2,ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_4,ind_3_2) = mat(ind_3_4,ind_3_2) - indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_6,ind_3_2) = mat(ind_3_6,ind_3_2) - indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)

  loss(ind_3_2) = loss(ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  prod(ind_3_4) = prod(ind_3_4) + indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)
  prod(ind_3_6) = prod(ind_3_6) + indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)

end do

! reactions a + b -> c + d 

eps = 1.d-10

do i4 = 1,nb_reaction_4_max

  ind_4_2 = indice_4(i4)%z2
  ind_4_4 = indice_4(i4)%z4
  ind_4_6 = indice_4(i4)%z6
  ind_4_8 = indice_4(i4)%z8

  eps_4 = abs(c(ilev,ind_4_2))/(abs(c(ilev,ind_4_2)) + abs(c(ilev,ind_4_4)) + eps)
  eps_4 = min(eps_4,1.0_jprb)

  mat(ind_4_2,ind_4_2) = mat(ind_4_2,ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4) 
  mat(ind_4_2,ind_4_4) = mat(ind_4_2,ind_4_4) + indice_4(i4)%z1*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_4,ind_4_2) = mat(ind_4_4,ind_4_2) + indice_4(i4)%z3*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_4,ind_4_4) = mat(ind_4_4,ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)   
  mat(ind_4_6,ind_4_2) = mat(ind_4_6,ind_4_2) - indice_4(i4)%z5*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_6,ind_4_4) = mat(ind_4_6,ind_4_4) - indice_4(i4)%z5*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_8,ind_4_2) = mat(ind_4_8,ind_4_2) - indice_4(i4)%z7*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_8,ind_4_4) = mat(ind_4_8,ind_4_4) - indice_4(i4)%z7*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)

  loss(ind_4_2) = loss(ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*c(ilev,ind_4_4)
  loss(ind_4_4) = loss(ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*c(ilev,ind_4_2)
  prod(ind_4_6) = prod(ind_4_6) + indice_4(i4)%z5*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)
  prod(ind_4_8) = prod(ind_4_8) + indice_4(i4)%z7*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)

end do

end subroutine fill_matrix

!================================================================

 subroutine indice(i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
                   i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
                   i_n, i_n2d, i_no, i_no2, i_n2)

!================================================================
! set the "indice" arrays used to fill the jacobian matrix      !
!----------------------------------------------------------------
! reaction                                   type               !
!----------------------------------------------------------------
! A + hv   --> B + C     photolysis          indice_phot        ! 
! A + B    --> C + D     bimolecular         indice_4           !
! A + A    --> B + C     quadratic           indice_3           !
! A + C    --> B + C     quenching           indice_phot        !
! A + ice  --> B + C     heterogeneous       indice_phot        !
!================================================================

use types_asis

implicit none

#include "chimiedata.h"

! input

integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
           i_n, i_n2d, i_no, i_no2, i_n2

! local

integer :: nb_phot, nb_reaction_3, nb_reaction_4
integer :: i_dummy

allocate (indice_phot(nb_phot_max))
allocate (indice_3(nb_reaction_3_max))
allocate (indice_4(nb_reaction_4_max))

i_dummy = 1

nb_phot       = 0
nb_reaction_3 = 0
nb_reaction_4 = 0

!===========================================================
!      O2 + hv -> O + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 2.0, i_o, 0.0, i_dummy) 

!===========================================================
!      O2 + hv -> O + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 1.0, i_o, 1.0, i_o1d) 

!===========================================================
!      CO2 + hv -> CO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o) 

!===========================================================
!      CO2 + hv -> CO + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o1d) 

!===========================================================
!      O3 + hv -> O2 + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2, 1.0, i_o1d) 

!===========================================================
!      O3 + hv -> O2 + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2, 1.0, i_o) 

!===========================================================
!      H2O + hv -> H + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o, 1.0, i_h, 1.0, i_oh) 

!===========================================================
!      H2O2 + hv -> OH + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      HO2 + hv -> OH + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 1.0, i_oh, 1.0, i_o) 

!===========================================================
!      NO + hv -> N + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_n, 1.0, i_o) 

!===========================================================
!      NO2 + hv -> NO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      a001 : O + O2 + CO2 -> O3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o2, 1.0, i_o3, 0.0, i_dummy) 

!===========================================================
!      a002 : O + O + CO2 -> O2 + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_o, 1.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      a003 : O + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      b001 : O(1D) + CO2 -> O + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy) 

!===========================================================
!      b002 : O(1D) + H2O -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2o, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      b003 : O(1D) + H2 -> OH + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2, 1.0, i_oh, 1.0, i_h) 

!===========================================================
!      b004 : O(1D) + O2 -> O + O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy) 

!===========================================================
!      b005 : O(1D) + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      b006 : O(1D) + O3 -> O2 + O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 1.0, i_o2, 2.0, i_o) 

!===========================================================
!      c001 : O + HO2 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_ho2, 1.0, i_oh, 1.0, i_o2) 

!===========================================================
!      c002 : O + OH -> O2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_oh, 1.0, i_o2, 1.0, i_h) 

!===========================================================
!      c003 : H + O3 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o3, 1.0, i_oh, 1.0, i_o2) 

!===========================================================
!      c004 : H + HO2 -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      c005 : H + HO2 -> H2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2, 1.0, i_o2) 

!===========================================================
!      c006 : H + HO2 -> H2O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o) 

!===========================================================
!      c007 : OH + HO2 -> H2O + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o2) 

!===========================================================
!      c008 : HO2 + HO2 -> H2O2 + O2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2) 

!===========================================================
!      c009 : OH + H2O2 -> H2O + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2o2, 1.0, i_h2o, 1.0, i_ho2) 

!===========================================================
!      c010 : OH + H2 -> H2O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2, 1.0, i_h2o, 1.0, i_h) 

!===========================================================
!      c011 : H + O2 + CO2 -> HO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o2, 1.0, i_ho2, 0.0, i_dummy) 

!===========================================================
!      c012 : O + H2O2 -> OH + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_h2o2, 1.0, i_oh, 1.0, i_ho2) 

!===========================================================
!      c013 : OH + OH -> H2O + O
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o, 1.0, i_o) 

!===========================================================
!      c014 : OH + O3 -> HO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_o3, 1.0, i_ho2, 1.0, i_o2) 

!===========================================================
!      c015 : HO2 + O3 -> OH + O2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_o3, 1.0, i_oh, 2.0, i_o2) 

!===========================================================
!      c016 : HO2 + HO2 + CO2 -> H2O2 + O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2) 

!===========================================================
!      c017 : OH + OH + CO2 -> H2O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o2, 0.0, i_dummy) 

!===========================================================
!      c018 : H + H + CO2 -> H2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_h, 1.0, i_h2, 0.0, i_dummy) 

!===========================================================
!      d001 : NO2 + O -> NO + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_o, 1.0, i_no, 1.0, i_o2) 

!===========================================================
!      d002 : NO + O3 -> NO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_o3, 1.0, i_no2, 1.0, i_o2) 

!===========================================================
!      d003 : NO + HO2 -> NO2 + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_ho2, 1.0, i_no2, 1.0, i_oh) 

!===========================================================
!      d004 : N + NO -> N2 + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_no, 1.0, i_n2, 1.0, i_o) 

!===========================================================
!      d005 : N + O2 -> NO + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      d006 : NO2 + H -> NO + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_h, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d007 : N + O -> NO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o, 1.0, i_no, 0.0, i_dummy) 

!===========================================================
!      d008 : N + HO2 -> NO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_ho2, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d009 : N + OH -> NO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_oh, 1.0, i_no, 1.0, i_h) 

!===========================================================
!      e001 : CO + OH -> CO2 + H 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_oh, 1.0, i_co2, 1.0, i_h) 

!===========================================================
!      e002 : CO + O + M -> CO2 + M 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_o, 1.0, i_co2, 0.0, i_dummy) 

!===========================================================
!      h001: HO2 + ice -> products
!            treated as
!            HO2 -> 0.5 H2O + 0.75 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 0.5, i_h2o, 0.75, i_o2) 

!===========================================================
!      h002: OH + ice -> products
!            treated as
!            OH -> 0.5 H2O + 0.25 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_oh, 0.5, i_h2o, 0.25, i_o2) 

!===========================================================
!      h003: H2O2 + ice -> products
!            treated as
!            H2O2 -> H2O + 0.5 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 1.0, i_h2o, 0.5, i_o2) 

!===========================================================
!      h004: HO2 + dust -> products
!            treated as
!            HO2 -> 0.5 H2O + 0.75 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 0.5, i_h2o, 0.75, i_o2) 

!===========================================================
!      h005: H2O2 + dust -> products
!            treated as
!            H2O2 -> H2O + 0.5 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 1.0, i_h2o, 0.5, i_o2) 

!===========================================================
!  check dimensions 
!===========================================================

print*, 'nb_phot       = ', nb_phot
print*, 'nb_reaction_4 = ', nb_reaction_4
print*, 'nb_reaction_3 = ', nb_reaction_3

if ((nb_phot /= nb_phot_max)             .or.  &
    (nb_reaction_3 /= nb_reaction_3_max) .or.  &
    (nb_reaction_4 /= nb_reaction_4_max)) then
   print*, 'wrong dimensions in indice' 
   stop
end if  

end subroutine indice

!*****************************************************************

      subroutine gcmtochim(nlayer, nq, zycol, lswitch, nesp,         &
                           i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
                           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,         &
                           i_n, i_n2d, i_no, i_no2, i_n2,            &
                           dens, rm, c) 

!*****************************************************************

      use tracer_h, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,         &
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,  &
     &                      igcm_ho2, igcm_h2o2, igcm_h2o_vap,           &
     &                      igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2

      use callkeys_mod

      implicit none


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     input:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      integer, intent(in) :: nlayer ! number of atmospheric layers
      integer, intent(in) :: nq     ! number of tracers in the gcm
      integer :: nesp               ! number of species in the chemistry
      integer :: lswitch            ! interface level between chemistries

      integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,           &
                 i_h2, i_oh, i_ho2, i_h2o2, i_h2o,                   &
                 i_n, i_n2d, i_no, i_no2, i_n2

      real :: zycol(nlayer,nq)      ! volume mixing ratios in the gcm
      real :: dens(nlayer)          ! total number density (molecule.cm-3) 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     output:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real, dimension(nlayer,nesp)            :: rm ! volume mixing ratios 
      real (kind = 8), dimension(nlayer,nesp) :: c  ! number densities (molecule.cm-3) 
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     local:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer      :: l, iesp
      logical,save :: firstcall = .true.
      
      
!     first call initializations

      if (firstcall) then

!       identify the indexes of the tracers we need

         if (igcm_co2 == 0) then
            write(*,*) "gcmtochim: Error; no CO2 tracer !!!"
            stop
         endif
         if (igcm_co == 0) then
            write(*,*) "gcmtochim: Error; no CO tracer !!!"
            stop
         end if
         if (igcm_o == 0) then
            write(*,*) "gcmtochim: Error; no O tracer !!!"
            stop
         end if
         if (igcm_o1d == 0) then
            write(*,*) "gcmtochim: Error; no O1D tracer !!!"
            stop
         end if
         if (igcm_o2 == 0) then
            write(*,*) "gcmtochim: Error; no O2 tracer !!!"
            stop
         end if
         if (igcm_o3 == 0) then
            write(*,*) "gcmtochim: Error; no O3 tracer !!!"
            stop
         end if
         if (igcm_h == 0) then
            write(*,*) "gcmtochim: Error; no H tracer !!!"
            stop
         end if
         if (igcm_h2 == 0) then
            write(*,*) "gcmtochim: Error; no H2 tracer !!!"
            stop
         end if
         if (igcm_oh == 0) then
            write(*,*) "gcmtochim: Error; no OH tracer !!!"
            stop
         end if
         if (igcm_ho2 == 0) then
            write(*,*) "gcmtochim: Error; no HO2 tracer !!!"
            stop
         end if
         if (igcm_h2o2 == 0) then
            write(*,*) "gcmtochim: Error; no H2O2 tracer !!!"
            stop
         end if
         if (igcm_n == 0) then
            write(*,*) "gcmtochim: Error; no N tracer !!!"
!            stop
         end if
         if (igcm_n2d == 0) then
            write(*,*) "gcmtochim: Error; no N2D tracer !!!"
!            stop
         end if
         if (igcm_no == 0) then
            write(*,*) "gcmtochim: Error; no NO tracer !!!"
!            stop
         end if
         if (igcm_no2 == 0) then
            write(*,*) "gcmtochim: Error; no NO2 tracer !!!"
!            stop
         end if
         if (igcm_n2 == 0) then
            write(*,*) "gcmtochim: Error; no N2 tracer !!!"
            stop
         end if
         if (igcm_h2o_vap == 0) then
            write(*,*) "gcmtochim: Error; no water vapor tracer !!!"
            stop
         end if
         firstcall = .false.
      end if ! of if (firstcall)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initialise mixing ratios
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1,lswitch-1
         rm(l,i_co2)  = zycol(l, igcm_co2)
         rm(l,i_co)   = zycol(l, igcm_co)
         rm(l,i_o)    = zycol(l, igcm_o)
         rm(l,i_o1d)  = zycol(l, igcm_o1d)
         rm(l,i_o2)   = zycol(l, igcm_o2)
         rm(l,i_o3)   = zycol(l, igcm_o3)
         rm(l,i_h)    = zycol(l, igcm_h)
         rm(l,i_h2)   = zycol(l, igcm_h2)
         rm(l,i_oh)   = zycol(l, igcm_oh)
         rm(l,i_ho2)  = zycol(l, igcm_ho2)
         rm(l,i_h2o2) = zycol(l, igcm_h2o2)
         rm(l,i_h2o)  = zycol(l, igcm_h2o_vap)
         rm(l,i_n)    = zycol(l, igcm_n)
         rm(l,i_n2d)  = zycol(l, igcm_n2d)
         rm(l,i_no)   = zycol(l, igcm_no)
         rm(l,i_no2)  = zycol(l, igcm_no2)
         rm(l,i_n2)   = zycol(l, igcm_n2)
      end do 

      where (rm(:,:) < 1.e-30)
         rm(:,:) = 0.
      end where

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initialise number densities
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iesp = 1,nesp
         do l = 1,lswitch-1
            c(l,iesp) = rm(l,iesp)*dens(l)
         end do
      end do

      end subroutine gcmtochim

!*****************************************************************
 
      subroutine chimtogcm(nlayer, nq, zycol, lswitch, nesp,         &
                           i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
                           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,         &
                           i_n, i_n2d, i_no, i_no2, i_n2, dens, c) 
 
!*****************************************************************
 
      use tracer_h, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,          &
                            igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,   &
                            igcm_ho2, igcm_h2o2, igcm_h2o_vap,            &
                            igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2

      use callkeys_mod

      implicit none


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     inputs:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      integer, intent(in) :: nlayer  ! number of atmospheric layers
      integer, intent(in) :: nq      ! number of tracers in the gcm
      integer :: nesp                ! number of species in the chemistry
      integer :: lswitch             ! interface level between chemistries
      integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
                 i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
                 i_n, i_n2d, i_no, i_no2, i_n2

      real :: dens(nlayer)     ! total number density (molecule.cm-3) 
      real (kind = 8), dimension(nlayer,nesp) :: c  ! number densities (molecule.cm-3) 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     output:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
      real zycol(nlayer,nq)  ! volume mixing ratios in the gcm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     local:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer l
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     save mixing ratios for the gcm
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1,lswitch-1
         zycol(l, igcm_co2)     = c(l,i_co2)/dens(l) 
         zycol(l, igcm_co)      = c(l,i_co)/dens(l) 
         zycol(l, igcm_o)       = c(l,i_o)/dens(l) 
         zycol(l, igcm_o1d)     = c(l,i_o1d)/dens(l)
         zycol(l, igcm_o2)      = c(l,i_o2)/dens(l) 
         zycol(l, igcm_o3)      = c(l,i_o3)/dens(l) 
         zycol(l, igcm_h)       = c(l,i_h)/dens(l)  
         zycol(l, igcm_h2)      = c(l,i_h2)/dens(l) 
         zycol(l, igcm_oh)      = c(l,i_oh)/dens(l) 
         zycol(l, igcm_ho2)     = c(l,i_ho2)/dens(l) 
         zycol(l, igcm_h2o2)    = c(l,i_h2o2)/dens(l)
         zycol(l, igcm_h2o_vap) = c(l,i_h2o)/dens(l)
         zycol(l, igcm_n)       = c(l,i_n)/dens(l)
         zycol(l, igcm_n2d)     = c(l,i_n2d)/dens(l)
         zycol(l, igcm_no)      = c(l,i_no)/dens(l)
         zycol(l, igcm_no2)     = c(l,i_no2)/dens(l)
         zycol(l, igcm_n2)      = c(l,i_n2)/dens(l)
      end do 

      end subroutine chimtogcm

end subroutine photochemistry_asis

