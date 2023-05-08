!*****************************************************************
!
!     Photochemical routine 
!
!     Author: Franck Lefevre
!             Benjamin Charnay
!     ------
!
!     Version: 27/05/2016
!
!*****************************************************************

subroutine photochemistry_asis_early_earth(nlayer, nq, ngrid,                   &
                          ig, lswitch, zycol, sza, fractcol, ptimestep, press,  &
                          temp, dens, zmmean, dist_sol, surfdust1d,          &
                          surfice1d, jo3, jch4, tau, iter)

      use callkeys_mod
implicit none

#include "chimiedata_early_earth.h" 

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
real    :: jch4(nlayer)        ! photodissociation rate ch4 -> ch3 + h

!===================================================================
!     local:
!===================================================================

integer :: phychemrat         ! (physical timestep)/(nominal chemical timestep)
integer :: j_o3_o1d, j_ch4_ch3_h, ilev
integer :: j_ch4_1ch2_h2, j_ch4_3ch2_h_h, j_ch4_ch_h2_h
integer :: iesp, nesp
integer :: lswitch

logical, save :: firstcall = .true.

parameter (nesp = 35)         ! number of species in the chemistry

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

integer,parameter :: i_ch4  = 18
integer,parameter :: i_ch3  = 19
integer,parameter :: i_ch   = 20
integer,parameter :: i_3ch2 = 21
integer,parameter :: i_1ch2 = 22
integer,parameter :: i_cho  = 23
integer,parameter :: i_ch2o = 24
integer,parameter :: i_ch3o = 25
integer,parameter :: i_c    = 26
integer,parameter :: i_c2   = 27
integer,parameter :: i_c2h  = 28
integer,parameter :: i_c2h2 = 29
integer,parameter :: i_c2h3 = 30
integer,parameter :: i_c2h4 = 31
integer,parameter :: i_c2h6 = 32
integer,parameter :: i_ch2co = 33
integer,parameter :: i_ch3co = 34
integer,parameter :: i_hcaer = 35



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
               i_n, i_n2d, i_no, i_no2, i_n2,                  &
               i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,    &
               i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,       &
               i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer)
   firstcall = .false.
end if

!===================================================================
!     initialisation of mixing ratios and densities       
!===================================================================

call gcmtochim(nlayer, nq, zycol, lswitch, nesp,               &
               i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
               i_n, i_n2d, i_no, i_no2, i_n2,                  &
               i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,    &
               i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,       &
               i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer, dens, rm, c) 

!===================================================================
!     interpolation of photolysis rates in the lookup table      
!===================================================================

call photolysis_asis_early_earth(nlayer, ngrid, lswitch, press, temp, sza, fractcol,tau, zmmean, dist_sol, & 
                     rm(:,i_co2), rm(:,i_o3), rm(:,i_ch4), v_phot)

! save o3 photolysis for output

j_o3_o1d = 5
jo3(:) = v_phot(:,j_o3_o1d)

j_ch4_ch3_h = 10
j_ch4_1ch2_h2=11
j_ch4_3ch2_h_h=12
j_ch4_ch_h2_h=13
jch4(:) = v_phot(:,j_ch4_ch3_h)+v_phot(:,j_ch4_1ch2_h2)+v_phot(:,j_ch4_3ch2_h_h)+v_phot(:,j_ch4_ch_h2_h)

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
#   write(*,*) "photochemistry_asis error, missing LAPACK routine dgesv"
#   stop
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
               i_n, i_n2d, i_no, i_no2, i_n2,                 &
               i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,   &
               i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,      &
               i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer, dens, c) 

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

#include "chimiedata_early_earth.h"

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
                             h001, h002, h003, h004, h005,               &       
                             f001, f002, f003, f004, f005, f006, f007,   &
                             f008, f009, f010, f011, f012, f013, f014,   &
                             f015, f016, f017, f018, f019, f020, f021,   &
                             f022, f023, f024, f025, f026, f027, f028,   &
                             f029, f030, f031, f032, f033, f034, f035,   &
                             f036, f037, f038, f039, f040, f041, f042,   &
                             f043, f044, f045, f046, f047, f048, f049,   &
                             f050, f051, f052, f053, f054, f055, f056,   &
                             f057, f058, f059, f060, f061, f062, f063,   &
                             f064, f065, f066, f067, f068, f069, f070,   &
                             f071, f072, f073, f074, f075, f076, f077,   &
                             f078, f079, f080, f081, f082


!----------------------------------------------------------------------
!     initialisation
!----------------------------------------------------------------------

      nb_phot       = 26       ! jearly_earth à 28 ou 26 (sans n/no/no2)
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

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b007(:)
!---  b008: o(1d) + ch4 -> ch3o + h

!     jpl 2003

      b008(:) = 1.5e-10*0.20

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b008(:)
!
!---  b009: o(1d) + ch4 -> ch2o + h2

!     jpl 2003

      b009(:) = 1.5e-10*0.05

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b009(:)
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

if(1.eq.0) then

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

endif

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


if(1.eq.1) then

!----------------------------------------------------------------------
!     methane/ethane reactions
!----------------------------------------------------------------------
!---  f001: 3ch2 + 3ch2 -> c2h2 + h2

!     Braun 1970
      
      f001(:) = 5.3e-11

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f001(:)

!---  f002: c + h2 + m -> 3ch2 

!     Zahnle 1986

      f002(:) = 8.3e-11*8.75e-31*exp(524./t(:))*dens(:)/(dens(:)*8.75e-31*exp(524./t(:))+8.3e-11)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f002(:)

!---  f003: c + o2 -> co + o

!     Donovan and Husain 1970

      f003(:) =3.3e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f003(:)

!---  f004: c + oh -> co + h

!     Giguere and Huebner 1978

      f004(:) =4.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f004(:)

!---  f005: c2 + ch4 -> c2h + ch3

!     Pitts 1982

      f005(:) =5.05e-11*exp(-297./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f005(:)

!---  f006: c2 + h2 -> c2h + h

!     Pitts 1982

      f006(:) =1.77e-10*exp(-1469./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f006(:)

!---  f007: c2 + o -> c + co

!     Prasad and Huntress 1980

      f007(:) =5.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f007(:)


!---  f008: c2 + o2 -> co + co

!     Baughcum and Oldenborg 1984

      f008(:) =1.5e-11*exp(-550./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f008(:)

!---  f009: c2h + c2h2 -> hcaer + h

!     Stephens 1987

      f009(:) =1.5e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f009(:)

!---  f010: c2h + ch4 -> c2h2 + ch3

!     Allen 1992 Lander 1990

      f010(:) =6.94e-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f010(:)


!---  f011: c2h + h + m -> c2h2 + m

!     Tsang 1986

      f011(:) =1.26e-18*exp(-721./t(:))*t(:)**-3.1*3.0e-10*dens(:)/(3.0e-10+dens(:)*1.26e-18*exp(-721./t(:))*t(:)**-3.1)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f011(:)

!---  f012: c2h + h2 -> c2h2 + h

!     Allen 1992 Stephens 1987

      f012(:) =5.58e-11*exp(-1443./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f012(:)

!---  f013: c2h + o -> co + ch

!     Zahnle 1986

      f013(:) =1.0e-10*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f013(:)

!---  f014: c2h + o2 -> co + cho

!     Brown 1981

      f014(:) =2.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f014(:)

!---  f015: c2h2 + h + m -> c2h3 + m

!     Romani 1993

      f015(:) =2.6e-31*8.3e-11*exp(-1374./t(:))*dens(:)/(2.6e-31*dens(:) + 8.3e-11*exp(-1374./t(:)) )

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f015(:)

!---  f016: c2h2 + o -> 3ch2 + co

!     Zahnle 1986

      f016(:) =2.9e-11*exp(-1600./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f016(:)

!---  f017: c2h2 + oh + m -> ch2co + h +  m

!     Perry and Williamson 1982

      f017(:) =5.8e-31*exp(1258./t(:))*1.4e-12*exp(388./t(:))*dens(:)/(dens(:)*5.8e-31*exp(1258./t(:)) + 1.4e-12*exp(388./t(:)))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f017(:)

!---  f018: c2h2 + oh -> co + ch3

!     Hampson and Garvin 1977

      f018(:) = 2.0e-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f018(:)

!---  f019: c2h3 + c2h3 -> c2h4 + c2h2

!     Fahr 1991

      f019(:) = 2.4e-11

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f019(:)

!---  f020: c2h3 + ch3 -> c2h2 + ch4

!     Fahr 1991

      f020(:) = 3.4e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f020(:)

!---  f021: c2h3 + ch4 -> c2h4 + ch3

!     Tsang 1986

      f021(:) = 2.4e-24*exp(-2754./t(:))*t(:)**4.02

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f021(:)

!---  f022: c2h3 + h -> c2h2 + h2

!     Warnatz 1984

      f022(:) = 3.3e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f022(:)

!---  f023: c2h3 + h2 -> c2h4 + h

!     Allen 1992

      f023(:) = 2.6e-13*exp(-2646./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f023(:)

!---  f024: c2h3 + o -> ch2co + h

!     Hoyermann 1981

      f024(:) = 5.5e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f024(:)

!---  f025: c2h3 + oh -> c2h2 + h2o

!     Benson 1967

      f025(:) = 8.3e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f025(:)

!---  f026: c2h4 + o -> cho + ch3

!     Hampson and Garvin 1977

      f026(:) = 5.5e-12*exp(-565./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f026(:)

!---  f027: c2h4 + oh -> ch2o + ch3

!     Hampson and Garvin 1977

      f027(:) = 2.2e-12*exp(385./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f027(:)

!---  f028: ch + ch4 + m -> c2h4 + h + m

!     Romani 1993

      f028(:) = min(2.5e-11*exp(200./t(:)),1.7e-10)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f028(:)

!---  f029: ch + co2 -> cho + co

!     Berman 1982

      f029(:) = 5.9e-12*exp(-350./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f029(:)

!---  f030: ch + h -> c + h2

!     Becker 1989

      f030(:) = 1.4e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f030(:)

!---  f031: ch + h2 -> 3ch2 + h

!     Zabarnick 1986

      f031(:) = 2.38e-10*exp(-1760./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f031(:)

!---  f032: ch + h2 + m -> ch3 + m

!     Romani 1993

      f032(:) = 8.75e-31*exp(524./t(:))*8.3e-11*dens(:)/(dens(:)*8.75e-31*exp(524./t(:)) + 8.3e-11)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f032(:)

!---  f033: ch + o -> co + h

!     Messing 1981

      f033(:) = 9.5e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f033(:)

!---  f034: ch + o2 -> co + oh

!     Butler 1981

      f034(:) = 5.9e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f034(:)

!---  f035: 1ch2 + ch4 -> ch3 + ch3

!     Bohland 1985

      f035(:) = 7.14e-12*exp(-5050./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f035(:)

!---  f036: 1ch2 + co2 -> ch2o + co

!     Zahnle 86

      f036(:) = 1e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f036(:)

!---  f037: 1ch2 + h2 -> 3ch2 + h2

!     Romani 1993

      f037(:) = 1.26e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f037(:)

!---  f038: 1ch2 + h2 -> ch3 + h

!     Tsang 1986

      f038(:) = 5.0e-15

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f038(:)

!---  f039: 1ch2 + m -> 3ch2 + m

!     Ashfold 1981

      f039(:) = 8.8e-12

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = f039(:)

!---  f040: 1ch2 + o2 -> cho + oh

!     Ashfold 1981

      f040(:) = 3.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f040(:)

!---  f041: 3ch2 + c2h3 -> ch3 + c2h2

!     Tsang 1986

      f041(:) = 3.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f041(:)

!---  f042: 3ch2 + ch3 -> c2h4 + h

!     Tsang 1986

      f042(:) = 7.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f042(:)

!---  f043: 3ch2 + co + m -> ch2co + m

!     Yung 1984

      f043(:) = 1.0e-28*1e-15*dens(:)/(dens(:)*1.0e-28 + 1e-15)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f043(:)


!---  f044: 3ch2 + co2 -> ch2o + co

!     Laufer 1981

      f044(:) = 1.0e-14

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f044(:)

!---  f045: 3ch2 + h  -> ch + h2

!     Zabarnick  1986

      f045(:) = 4.7e-10*exp(-370./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f045(:)

!---  f046: 3ch2 + h + m -> ch3 +  m

!     Gladstone  1996

      f046(:) = 1.5e-10*3.1e-30*exp(475./t(:))/(1.5e-10 + dens(:)*3.1e-30*exp(475./t(:)))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f046(:)

!---  f047: 3ch2 + o -> ch + oh

!     Huebner  1980

      f047(:) = 8.0e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f047(:)

!---  f048: 3ch2 + o -> co + h+h

!     Homann  1983

      f048(:) = 8.3e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f048(:)

!---  f049: 3ch2 + o -> cho + h

!     Homann  1983

      f049(:) = 1.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f049(:)

!---  f050: 3ch2 + o2 -> hco + oh

!     Baulch  1994

      f050(:) = 4.1e-11*exp(-750./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f050(:)

!---  f051: ch2co + h -> ch3 + co

!     Michael  1979

      f051(:) = 1.9e-11*exp(-1725./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f051(:)

!---  f052: ch2co + o -> h2co + co

!     Lee  1980

      f052(:) = 3.3e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f052(:)

!---  f053: ch3 + ch3 + m -> c2h6 + m

!     Wagner  1988

      do ilev = 1,lswitch-1
         ak0 = 1.17e-25*exp(-500./t(ilev))*(t(ilev)/300.)**(-3.75)
         ak1 = 3.0e-11*(t(ilev)/300.)**(-1.0)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2.)
         f053(ilev) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f053(:)

!---  f054: ch3 + co + m -> ch3co + m

!     Watkins  1974

      f054(:) = 1.4e-32*exp(-3000./t(:))*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f054(:)

!---  f055: ch3 + h + m -> ch4 + m

!     Baulch  1994


      do ilev = 1,lswitch-1
         ak0 = 1.0e-28*(t(ilev)/300.)**(-1.8)
         ak1 = 2.e-10*(t(ilev)/300.)**(-0.4)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2.)
         f055(ilev) = rate*0.6**xpo
      end do


      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f055(:)

!---  f056: ch3 + ch2o -> ch4 + cho

!     Baulch  1994

      f056(:) = 1.6e-16*exp(899./t(:))*(t(:)/298.)**6.1
!      f056(:) = 4.9e-15*(t(:)/298.)**4.4 *exp(-2450./t(:)) 

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f056(:)

!---  f057: ch3 + cho -> ch4 + co

!     Tsang  1986

      f057(:) = 5.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f057(:)

!---  f058: ch3 + o -> ch2o + h

!     Sander  2006

      f058(:) = 1.1e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f058(:)

!---  f059: ch3 + o2 -> ch2o + oh

!     Sander  2006


      do ilev = 1,lswitch-1
         ak0 = 4.5e-31*(t(ilev)/300.)**(-3.)
         ak1 = 1.8e-12*(t(ilev)/300.)**(-1.7)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2.)
         f059(ilev) = rate*0.6**xpo
      end do


      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f059(:)

!---  f060: ch3 + o3 -> ch2o + ho2

!     Sander  2006

      f060(:) = 5.4e-12*exp(-220./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f060(:)

!---  f061: ch3 + o3 -> ch3o + o2

!     Sander  2006

      f061(:) = 5.4e-12*exp(-220.0/t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f061(:)

!---  f062: 3ch2 + c2h3 -> ch3 + c2h2

!     Tsang  1986

      f062(:) = 3.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f062(:)

!---  f063: ch3 + oh -> ch3o + h

!     Jasper  2007

      f063(:) = 9.3e-11*exp(-1606/t(:))*(t(:)/298.)
!      f063(:) = 1.3e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f063(:)

!---  f064: ch3 + oh -> co + h2 + h2

!     Fenimore [1969]

      f064(:) = 6.7e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f064(:)

!---  f065: ch3co + ch3 -> c2h6 + co

!     Adachi  1981

      f065(:) = 5.4e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f065(:)

!---  f066: ch3co + ch3 -> ch4 + ch2co

!     Adachi  1981

      f066(:) = 8.6e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f066(:)

!---  f067: ch3co + h -> ch4 + co

!     Zahnle  1986

      f067(:) = 1.0e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f067(:)

!---  f068: ch3co + o -> ch2o + cho

!     Zahnle  1986

      f068(:) = 5.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f068(:)

!---  f069: ch3o + co -> ch3 + co2

!     Wen  1989

      f069(:) = 2.6e-11*exp(-5940.0/t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f069(:)

!---  f070: ch4 + o -> ch3 + oh

!     Tsang  1986

      f070(:) = 8.75E-12*exp(-4330.0/t(:))*(t(:)/298.0)**1.5
!      f070(:) = 5.8E-11*exp(-4450.0/t(:))     

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f070(:)

!---  f071: ch4 + oh -> ch3 + h2o

!     Sander  2006

      f071(:) = 2.45e-12*exp(-1775.0/t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f071(:)

!---  f072: h + co + m -> cho + m

!     Baulch  1994

      f072(:) = 1.4e-34*exp(-100.0/t(:))*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f072(:)

!---  f073: h + cho -> h2 + co

!     Baulch  1992

      f073(:) = 1.8e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f073(:)

!---  f074: ch2o + h -> h2 + cho

!     Baulch  1994

      f074(:) = 2.14e-12*exp(-1090/t(:))*(t(:)/298)**1.62

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f074(:)

!---  f075: ch2o + o -> cho + oh

!     Sander  2006

      f075(:) = 3.4e-11*exp(-1600/t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f075(:)

!---  f076: ch2o + oh -> h2o + cho

!     Sander  2006

      f076(:) = 5.5e-12*exp(125/t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f076(:)

!---  f077: cho + ch2o -> ch3o + co

!     Wen  1989

      f077(:) = 3.8e-17

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f077(:)

!---  f078: cho + cho -> ch2o + co

!     Tsang  1986

      f078(:) = 4.5e-11

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f078(:)

!---  f079: cho + o2 -> ho2 + co

!     Sander 2006

      f079(:) = 5.2e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f079(:)

!---  f080: o + cho -> h + co2

!     Tsang 1986

      f080(:) = 5.0e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f080(:)

!---  f081: o + cho -> oh + co

!     Hampson 1977

      f081(:) = 1.0e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f081(:)

!---  f082: oh + cho -> h2o + co

!     Tsang 1986

      f082(:) = 1.0e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f082(:)


endif 




!----------------------------------------------------------------------
!     heterogeneous chemistry 
!----------------------------------------------------------------------

      if (1.eq.0) then
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

      endif
end subroutine reactionrates

!======================================================================

 subroutine fill_matrix(ilev, mat, prod, loss, c, nesp, nlayer, v_phot, v_3, v_4)

!======================================================================
! filling of the jacobian matrix
!======================================================================

use types_asis

implicit none

#include "chimiedata_early_earth.h"

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
                   i_n, i_n2d, i_no, i_no2, i_n2,                  &
                   i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,    &
                   i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,       &
                   i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer)

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

#include "chimiedata_early_earth.h"

! input

integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
           i_n, i_n2d, i_no, i_no2, i_n2,                  &
           i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,    &
           i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,       &
           i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer

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

!nb_phot = nb_phot + 1

!indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_n, 1.0, i_o) 

!===========================================================
!      NO2 + hv -> NO + O
!===========================================================

!nb_phot = nb_phot + 1

!indice_phot(nb_phot) = z3spec(1.0, i_no2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      HNO3 + hv -> NO2 + OH
!===========================================================

!nb_phot = nb_phot + 1

!indice_phot(nb_phot) = z3spec(1.0, i_hno3, 1.0, i_no2, 1.0, i_oh)

!===========================================================
!      HNO4 + hv -> NO2 + HO2
!===========================================================

!nb_phot = nb_phot + 1

!indice_phot(nb_phot) = z3spec(1.0, i_hno4, 1.0, i_no2, 1.0, i_ho2)

if(1.eq.1) then
!===========================================================
!      CH4 + hv -> CH3 + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch4, 1.0, i_ch3, 1.0, i_h)

!===========================================================
!      CH4 + hv -> 1CH2 + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch4, 1.0, i_1ch2, 1.0, i_h2)

!===========================================================
!      CH4 + hv -> 1CH2 + H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch4, 1.0, i_3ch2, 2.0, i_h)

!===========================================================
!      CH4 + hv -> CH + H2 + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(0.5, i_ch4, 1.0, i_ch, 0.0, i_dummy)

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(0.5, i_ch4, 1.0, i_h2, 1.0, i_h)

!if(1.eq.0) then
!===========================================================
!      CH2O + hv -> CHO + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch2o, 1.0, i_cho, 1.0, i_h)

!===========================================================
!      CH2O + hv -> CO + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch2o, 1.0, i_co, 1.0, i_h2)

!===========================================================
!      C2H6 + hv -> CH4 + 1CH2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h6, 1.0, i_ch4, 1.0, i_1ch2)

!===========================================================
!      C2H6 + hv -> C2H2 + H2 + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h6, 1.0, i_c2h2, 2.0, i_h2)

!===========================================================
!      C2H6 + hv -> C2H4 + H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h6, 1.0, i_c2h4, 2.0, i_h)

!===========================================================
!      C2H6 + hv -> C2H4 + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h6, 1.0, i_c2h4, 1.0, i_h2)

!===========================================================
!      C2H6 + hv -> CH3 + CH3
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h6, 2.0, i_ch3, 0.0, i_dummy)

!===========================================================
!      C2H2 + hv -> C2H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h2, 1.0, i_c2h, 1.0, i_h)

!===========================================================
!      C2H2 + hv -> C2 + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h2, 1.0, i_c2, 1.0, i_h2)

!===========================================================
!      C2H4 + hv -> C2H2 + H2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h4, 1.0, i_c2h2, 1.0, i_h2)

!===========================================================
!      C2H4 + hv -> C2H2 + H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_c2h4, 1.0, i_c2h2, 2.0, i_h)

!===========================================================
!      CH2CO + hv -> 3CH2 + CO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ch2co, 1.0, i_3ch2, 1.0, i_co)

endif

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
!      b007 : O(1D) + CH4 -> CH3 + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_ch4, 1.0, i_ch3, 1.0, i_oh)

!===========================================================
!      b008 : O(1D) + CH4 -> CH3O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_ch4, 1.0, i_ch3o, 1.0, i_h)

!===========================================================
!      b009 : O(1D) + CH4 -> CH2O + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_ch4, 1.0, i_ch2o, 1.0, i_h2)


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

if(1.eq.0) then

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

endif

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


if(1.eq.1) then
!===========================================================
!      f001 : 3CH2 + 3CH2 -> C2H2 + H2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_3ch2, 1.0, i_c2h2, 1.0, i_h2)

!===========================================================
!      f002 : C + H2 + m -> 3CH2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c, 1.0, i_h2, 1.0, i_3ch2, 0.0, i_dummy)

!===========================================================
!      f003 : C + O2  -> CO + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c, 1.0, i_o2, 1.0, i_co, 1.0, i_o)

!===========================================================
!      f004 : C + OH -> CO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c, 1.0, i_oh, 1.0, i_co, 1.0, i_h)

!===========================================================
!      f005 : C2 + CH4 -> C2H + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2, 1.0, i_ch4, 1.0, i_c2h, 1.0, i_ch3)

!===========================================================
!      f006 : C2 + H2 -> C2H + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2, 1.0, i_h2, 1.0, i_c2h, 1.0, i_h)

!===========================================================
!      f007 : C2 + O -> C + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2, 1.0, i_o, 1.0, i_c, 1.0, i_co)

!===========================================================
!      f008 : C2 + O2 -> CO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2, 1.0, i_o2, 2.0, i_co, 0.0, i_dummy)

!===========================================================
!      f009 : C2H + C2H2 -> HCAER + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_c2h2, 1.0, i_hcaer, 1.0, i_h)

!===========================================================
!      f010 : C2H + CH4 -> C2H2 + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_ch4, 1.0, i_c2h2, 1.0, i_ch3)

!===========================================================
!      f011 : C2H + H + m -> C2H2 + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_h, 1.0, i_c2h2, 0.0, i_dummy)

!===========================================================
!      f012 : C2H + H2 -> C2H2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_h2, 1.0, i_c2h2, 1.0, i_h)

!===========================================================
!      f013 : C2H + O -> CO + CH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_o, 1.0, i_co, 1.0, i_ch)

!===========================================================
!      f014 : C2H + O2 -> CO + CHO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h, 1.0, i_o2, 1.0, i_co, 1.0, i_cho)

!===========================================================
!      f015 : C2H2 + H + m -> C2H3 + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h2, 1.0, i_h, 1.0, i_c2h3, 0.0, i_dummy)

!===========================================================
!      f016 : C2H2 + O -> 3CH2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h2, 1.0, i_o, 1.0, i_3ch2, 1.0, i_co)

!===========================================================
!      f017 : C2H2 + OH + m -> CH2CO + H + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h2, 1.0, i_oh, 1.0, i_ch2co, 1.0, i_h)

!===========================================================
!      f018 : C2H2 + OH -> CO + CH3 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h2, 1.0, i_oh, 1.0, i_co, 1.0, i_ch3)

!===========================================================
!      f019 : C2H3 + C2H3 -> C2H4 + C2H2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_c2h3, 1.0, i_c2h4, 1.0, i_c2h2)

!===========================================================
!      f020 : C2H3 + CH3 -> C2H2 + CH4
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_ch3, 1.0, i_c2h2, 1.0, i_ch4)

!===========================================================
!      f021 : C2H3 + CH4 -> C2H4 + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_ch4, 1.0, i_c2h4, 1.0, i_ch3)

!===========================================================
!      f022 : C2H3 + H -> C2H2 + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_h, 1.0, i_c2h2, 1.0, i_h2)

!===========================================================
!      f023 : C2H3 + H2 -> C2H4 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_h2, 1.0, i_c2h4, 1.0, i_h)

!===========================================================
!      f024 : C2H3 + O -> CH2CO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_o, 1.0, i_ch2co, 1.0, i_h)

!===========================================================
!      f025 : C2H3 + OH -> C2H2 + H2O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h3, 1.0, i_oh, 1.0, i_c2h2, 1.0, i_h2o)

!===========================================================
!      f026 : C2H4 + O -> CHO + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h4, 1.0, i_o, 1.0, i_cho, 1.0, i_ch3)

!===========================================================
!      f027 : C2H4 + OH -> CH2O + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_c2h4, 1.0, i_oh, 1.0, i_ch2O, 1.0, i_ch3)

!===========================================================
!      f028 : CH + CH4 -> C2H4 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_ch4, 1.0, i_c2h4, 1.0, i_h)

!===========================================================
!      f029 : CH + CO2 -> CHO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_co2, 1.0, i_cho, 1.0, i_co)

!===========================================================
!      f030 : CH + H -> C + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_h, 1.0, i_c, 1.0, i_h2)

!===========================================================
!      f031 : CH + H2 -> 3CH2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_h2, 1.0, i_3ch2, 1.0, i_h)

!===========================================================
!      f032 : CH + H2 + m -> CH3 +m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_h2, 1.0, i_ch3, 0.0, i_dummy)

!===========================================================
!      f033 : CH + O -> CO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_o, 1.0, i_co, 1.0, i_h)

!===========================================================
!      f034 : CH + O2 -> CO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch, 1.0, i_o2, 1.0, i_co, 1.0, i_oh)

!===========================================================
!      f035 : 1CH2 + CH4 -> CH3 + CH3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_1ch2, 1.0, i_ch4, 2.0, i_ch3, 0.0, i_dummy)

!===========================================================
!      f036 : 1CH2 + CO2 -> CH2O + CO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_1ch2, 1.0, i_co2, 1.0, i_ch2o, 1.0, i_co)

!===========================================================
!      f037 : 1CH2 + H2 -> 3CH2 + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_1ch2, 1.0, i_h2, 1.0, i_3ch2, 1.0,i_h2)

!===========================================================
!      f038 : 1CH2 + H2 -> CH3 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_1ch2, 1.0, i_h2, 1.0, i_ch3, 1.0, i_h)

!===========================================================
!      f039 : 1CH2 + m -> 3CH2 + m
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_1ch2, 1.0, i_3ch2, 0.0, i_dummy)


!===========================================================
!      f040 : 1CH2 + O2 -> CHO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_1ch2, 1.0, i_o2, 1.0, i_cho, 1.0, i_oh)

!===========================================================
!      f041 : 3CH2 + C2H3 -> CH3 + C2H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_c2h3, 1.0, i_ch3, 1.0, i_c2h2)

!===========================================================
!      f042 : 3CH2 + CH3 -> C2H4 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_ch3, 1.0, i_c2h4, 1.0, i_h)

!===========================================================
!      f043 : 3CH2 + CO + m -> CH2CO + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_co, 1.0, i_ch2co, 0.0, i_dummy)

!===========================================================
!      f044 : 3CH2 + CO2 -> CH2O + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_co2, 1.0, i_ch2o, 1.0, i_co)

!===========================================================
!      f045 : 3CH2 + H -> CH + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_h, 1.0, i_ch, 1.0, i_h2)

!===========================================================
!      f046 : 3CH2 + H + m -> CH3 + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_h, 1.0, i_ch3, 0.0, i_dummy)

!===========================================================
!      f047 : 3CH2 + O -> CH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_o, 1.0, i_ch, 1.0, i_oh)

!===========================================================
!      f048 : 3CH2 + O -> CO + H + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_o, 1.0, i_co, 2.0, i_h)

!===========================================================
!      f049 : 3CH2 + O -> CHO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_o, 1.0, i_cho, 1.0, i_h)

!===========================================================
!      f050 : 3CH2 + O2 -> CHO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_o2, 1.0, i_cho, 1.0, i_oh)

!===========================================================
!      f051 : CH2CO + H -> CH3 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch2co, 1.0, i_h, 1.0, i_ch3, 1.0, i_co)

!===========================================================
!      f052 : CH2CO + O -> CH2O + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch2co, 1.0, i_o, 1.0, i_ch2o, 1.0, i_co)

!===========================================================
!      f053 : CH3 + CH3 + m -> C2H6 + m
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ch3, 1.0, i_c2h6, 0.0, i_dummy)

!===========================================================
!      f054 : CH3 + CO + m -> CH3CO + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_co, 1.0, i_ch3co, 0.0, i_dummy)

!===========================================================
!      f055 : CH3 + H + m -> CH4 + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_h, 1.0, i_ch4, 0.0, i_dummy)

!===========================================================
!      f056 : CH3 + CH2O -> CH4 + CHO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_ch2o, 1.0, i_ch4, 1.0, i_cho)

!===========================================================
!      f057 : CH3 + CHO -> CH4 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_cho, 1.0, i_ch4, 1.0, i_co)

!===========================================================
!      f058 : CH3 + O -> CH2O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_o, 1.0, i_ch2o, 1.0, i_h)

!===========================================================
!      f059 : CH3 + O2 -> CH2O + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_o2, 1.0, i_ch2o, 1.0, i_oh)

!===========================================================
!      f060 : CH3 + O3 -> CH2O + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_o3, 1.0, i_ch2o, 1.0, i_ho2)

!===========================================================
!      f061 : CH3 + O3 -> CH3O + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_o3, 1.0, i_ch3o, 1.0, i_o2)

!===========================================================
!      f062 : 3CH2 + C2H3 -> CH3 + CH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_3ch2, 1.0, i_c2h3, 1.0, i_ch3, 1.0, i_ch)

!===========================================================
!      f063 : CH3 + OH -> CH3O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_oh, 1.0, i_ch3o, 1.0, i_h)

!===========================================================
!      f064 : CH3 + OH -> CO + H2 + H2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3, 1.0, i_oh, 1.0, i_co, 2.0, i_h2)

!===========================================================
!      f065 : CH3CO + CH3 -> C2H6 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3co, 1.0, i_ch3, 1.0, i_c2h6, 1.0, i_co)

!===========================================================
!      f066 : CH3CO + CH3 -> CH4 + CH2CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3co, 1.0, i_ch3, 1.0, i_ch4, 1.0, i_ch2co)

!===========================================================
!      f067 : CH3CO + H -> CH4 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3co, 1.0, i_h, 1.0, i_ch4, 1.0, i_co)

!===========================================================
!      f068 : CH3CO + O -> CH2O + CHO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3co, 1.0, i_o, 1.0, i_ch2o, 1.0, i_cho)

!===========================================================
!      f069 : CH3O + CO -> CH3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch3o, 1.0, i_co, 1.0, i_ch3, 1.0, i_co2)

!===========================================================
!      f070 : CH4 + O -> CH3 + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch4, 1.0, i_o, 1.0, i_ch3, 1.0, i_oh)

!===========================================================
!      f071 : CH4 + OH -> CH3 + H2O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch4, 1.0, i_oh, 1.0, i_ch3, 1.0, i_h2o)

!===========================================================
!      f072 : H + CO + m -> CHO + m
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_co, 1.0, i_cho, 0.0, i_dummy)

!===========================================================
!      f073 : H + CHO -> H2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_cho, 1.0, i_h2, 1.0, i_co)

!===========================================================
!      f074 : CH2O + H -> H2 + CHO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch2o, 1.0, i_h, 1.0, i_h2, 1.0, i_cho)

!===========================================================
!      f075 : CH2O + O -> CHO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch2o, 1.0, i_o, 1.0, i_cho, 1.0, i_oh)

!===========================================================
!      f076 : CH2O + OH -> H2O + CHO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ch2o, 1.0, i_oh, 1.0, i_h2o, 1.0, i_cho)

!===========================================================
!      f077 : CHO + CH2O -> CH3O + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cho, 1.0, i_ch2o, 1.0, i_ch3o, 1.0, i_co)

!===========================================================
!      f078 : CHO + CHO -> CH2O + CO
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_cho, 1.0, i_ch2o, 1.0, i_co)

!===========================================================
!      f079 : CHO + O2 -> HO2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cho, 1.0, i_o2, 1.0, i_ho2, 1.0, i_co)

!===========================================================
!      f080 : O + CHO -> H + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_cho, 1.0, i_h, 1.0, i_co2)

!===========================================================
!      f081 : CHO + O2 -> OH + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cho, 1.0, i_o2, 1.0, i_oh, 1.0, i_co)

!===========================================================
!      f082 : OH + CHO -> H2O + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_cho, 1.0, i_h2o, 1.0, i_co)

endif

if(1.eq.0) then
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

endif
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
   print*, 'nb_phot_max       = ', nb_phot_max
   print*, 'nb_reaction_4_max = ', nb_reaction_4_max
   print*, 'nb_reaction_3_max = ', nb_reaction_3_max
   stop
end if  

end subroutine indice

!*****************************************************************

      subroutine gcmtochim(nlayer, nq, zycol, lswitch, nesp,         &
                           i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
                           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,         &
                           i_n, i_n2d, i_no, i_no2, i_n2,            &
                           i_ch4, i_ch3, i_3ch2, i_1ch2, i_cho,      &
                           i_ch2o,i_ch3o, i_c, i_c2, i_c2h, i_c2h2,  &
                           i_c2h3, i_c2h4, i_c2h6, i_ch2co, i_ch3co, &
                           i_hcaer, dens, rm, c) 

!*****************************************************************

      use tracer_h, only:  igcm_co2, igcm_co, igcm_o, igcm_o1d,         &
                           igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,  &
                           igcm_ho2, igcm_h2o2, igcm_h2o_vap,           &
                           igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2, &
                           igcm_ch4, igcm_ch3, igcm_ch, igcm_3ch2,       &
                           igcm_1ch2, igcm_cho, igcm_ch2o, igcm_ch3o,    &
                           igcm_c, igcm_c2, igcm_c2h, igcm_c2h2,         &
                           igcm_c2h3, igcm_c2h4, igcm_c2h6, igcm_ch2co,  &
                           igcm_ch3co, igcm_hcaer


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
                 i_n, i_n2d, i_no, i_no2, i_n2, i_ch4,               &
                 i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,               &
                 i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,           &
                 i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer

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
            stop
         end if
         if (igcm_n2d == 0) then
            write(*,*) "gcmtochim: Error; no N2D tracer !!!"
            stop
         end if
         if (igcm_no == 0) then
            write(*,*) "gcmtochim: Error; no NO tracer !!!"
            stop
         end if
         if (igcm_no2 == 0) then
            write(*,*) "gcmtochim: Error; no NO2 tracer !!!"
            stop
         end if
         if (igcm_n2 == 0) then
            write(*,*) "gcmtochim: Error; no N2 tracer !!!"
            stop
         end if
         if (igcm_h2o_vap == 0) then
            write(*,*) "gcmtochim: Error; no water vapor tracer !!!"
            stop
         end if

         if (igcm_ch4 == 0) then
            write(*,*) "gcmtochim: Error; no CH4 tracer !!!"
            stop
         end if
         if (igcm_ch3 == 0) then
            write(*,*) "gcmtochim: Error; no CH3 tracer !!!"
            stop
         end if
         if (igcm_ch == 0) then
            write(*,*) "gcmtochim: Error; no CH tracer !!!"
            stop
         end if
         if (igcm_3ch2 == 0) then
            write(*,*) "gcmtochim: Error; no 3CH2 tracer !!!"
            stop
         end if
         if (igcm_1ch2 == 0) then
            write(*,*) "gcmtochim: Error; no 1CH2 tracer !!!"
            stop
         end if
         if (igcm_cho == 0) then
            write(*,*) "gcmtochim: Error; no CHO tracer !!!"
            stop
         end if
         if (igcm_ch2o == 0) then
            write(*,*) "gcmtochim: Error; no CH2O tracer !!!"
            stop
         end if
         if (igcm_ch3o == 0) then
            write(*,*) "gcmtochim: Error; no CH3O tracer !!!"
            stop
         end if
         if (igcm_c == 0) then
            write(*,*) "gcmtochim: Error; no C tracer !!!"
            stop
         end if
         if (igcm_c2 == 0) then
            write(*,*) "gcmtochim: Error; no C2 tracer !!!"
            stop
         end if
         if (igcm_c2h == 0) then
            write(*,*) "gcmtochim: Error; no C2H tracer !!!"
            stop
         end if
         if (igcm_c2h2 == 0) then
            write(*,*) "gcmtochim: Error; no C2H2 tracer !!!"
            stop
         end if
         if (igcm_c2h3 == 0) then
            write(*,*) "gcmtochim: Error; no C2H3 tracer !!!"
            stop
         end if
         if (igcm_c2h4 == 0) then
            write(*,*) "gcmtochim: Error; no C2H4 tracer !!!"
            stop
         end if
         if (igcm_c2h6 == 0) then
            write(*,*) "gcmtochim: Error; no C2H6 tracer !!!"
            stop
         end if
         if (igcm_ch2co == 0) then
            write(*,*) "gcmtochim: Error; no CH2CO tracer !!!"
            stop
         end if
         if (igcm_ch3co == 0) then
            write(*,*) "gcmtochim: Error; no CH3CO tracer !!!"
            stop
         end if
         if (igcm_hcaer == 0) then
            write(*,*) "gcmtochim: Error; no HCAER tracer !!!"
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
         rm(l,i_ch4)   = zycol(l, igcm_ch4)
         rm(l,i_ch3)   = zycol(l, igcm_ch3)
         rm(l,i_ch)   = zycol(l, igcm_ch)
         rm(l,i_3ch2)   = zycol(l, igcm_3ch2)
         rm(l,i_1ch2)   = zycol(l, igcm_1ch2)
         rm(l,i_cho)   = zycol(l, igcm_cho)
         rm(l,i_ch2o)   = zycol(l, igcm_ch2o)
         rm(l,i_ch3o)   = zycol(l, igcm_ch3o)
         rm(l,i_c)   = zycol(l, igcm_c)
         rm(l,i_c2)   = zycol(l, igcm_c2)
         rm(l,i_c2h)   = zycol(l, igcm_c2h)
         rm(l,i_c2h2)   = zycol(l, igcm_c2h2)
         rm(l,i_c2h3)   = zycol(l, igcm_c2h3)
         rm(l,i_c2h4)   = zycol(l, igcm_c2h4)
         rm(l,i_c2h6)   = zycol(l, igcm_c2h6)
         rm(l,i_ch2co)   = zycol(l, igcm_ch2co)
         rm(l,i_ch3co)   = zycol(l, igcm_ch3co)
         rm(l,i_hcaer)   = zycol(l, igcm_hcaer)
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
                           i_n, i_n2d, i_no, i_no2, i_n2, i_ch4,     &
                           i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,     &
                           i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3, &
                           i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer, dens, c) 
 
!*****************************************************************
 
      use tracer_h, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,            &
                            igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,   &
                            igcm_ho2, igcm_h2o2, igcm_h2o_vap,            &
                            igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2, &
                            igcm_ch4, igcm_ch3, igcm_ch, igcm_3ch2,       &
                            igcm_1ch2, igcm_cho, igcm_ch2o, igcm_ch3o,    &
                            igcm_c, igcm_c2, igcm_c2h, igcm_c2h2,         &
                            igcm_c2h3, igcm_c2h4, igcm_c2h6, igcm_ch2co,  &
                            igcm_ch3co, igcm_hcaer

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
                 i_n, i_n2d, i_no, i_no2, i_n2, i_ch4,           &
                 i_ch3, i_3ch2, i_1ch2, i_cho, i_ch2o,           &
                 i_ch3o, i_c, i_c2, i_c2h, i_c2h2, i_c2h3,       &
                 i_c2h4, i_c2h6, i_ch2co, i_ch3co, i_hcaer



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

         zycol(l, igcm_ch4)     = c(l,i_ch4)/dens(l)
         zycol(l, igcm_ch3)     = c(l,i_ch3)/dens(l)
         zycol(l, igcm_ch)     = c(l,i_ch)/dens(l)
         zycol(l, igcm_3ch2)     = c(l,i_3ch2)/dens(l)
         zycol(l, igcm_1ch2)     = c(l,i_1ch2)/dens(l)
         zycol(l, igcm_cho)     = c(l,i_cho)/dens(l)
         zycol(l, igcm_ch2o)     = c(l,i_ch2o)/dens(l)
         zycol(l, igcm_ch3o)     = c(l,i_ch3o)/dens(l)
         zycol(l, igcm_c)     = c(l,i_c)/dens(l)
         zycol(l, igcm_c2)     = c(l,i_c2)/dens(l)
         zycol(l, igcm_c2h)     = c(l,i_c2h)/dens(l)
         zycol(l, igcm_c2h2)     = c(l,i_c2h2)/dens(l)
         zycol(l, igcm_c2h3)     = c(l,i_c2h3)/dens(l)
         zycol(l, igcm_c2h4)     = c(l,i_c2h4)/dens(l)
         zycol(l, igcm_c2h6)     = c(l,i_c2h6)/dens(l)
         zycol(l, igcm_ch2co)     = c(l,i_ch2co)/dens(l)
         zycol(l, igcm_ch3co)     = c(l,i_ch3co)/dens(l)
         zycol(l, igcm_hcaer)     = c(l,i_hcaer)/dens(l)


      end do 

      end subroutine chimtogcm

end subroutine photochemistry_asis_early_earth

