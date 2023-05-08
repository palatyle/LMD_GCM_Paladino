subroutine photochemistry_venus(nz, n_lon, ptimestep, p, t, tr, mumean, sza_input, nesp, iter)

use chemparam_mod
      
implicit none

!===================================================================
!     input:
!===================================================================

integer, intent(in) :: nz         ! number of atmospheric layers
integer, intent(in) :: nesp       ! number of tracers in traceur.def

real, dimension(nz) :: p          ! pressure (hpa)
real, dimension(nz) :: t          ! temperature (k)
real, dimension(nz) :: mumean     ! mean molecular mass (g/mol)
real :: ptimestep                 ! physics timestep (s)
real :: sza_input                 ! solar zenith angle (degrees)

integer :: n_lon                  ! for 1D test

!===================================================================
!     input/output:
!===================================================================

real, dimension(nz,nesp) :: tr    ! tracer mixing ratio

!===================================================================
!     output:
!===================================================================

integer :: iter(nz)               ! iteration counter

!===================================================================
!     local:
!===================================================================

logical, save :: firstcall = .true.

real, dimension(nz)  :: conc      ! total number density (molecule.cm-3)
real, dimension(nz)  :: surfice1d, surfdust1d

! photolysis lookup table

integer, parameter :: nj = 19, nztable = 201, nsza = 27, nso2 = 13
real, dimension(nso2,nsza,nztable,nj), save :: jphot
real, dimension(nztable), save :: table_colair
real, dimension(nso2,nztable), save :: table_colso2
real, dimension(nsza), save :: table_sza
real :: dist_sol

! number densities

real, dimension(nesp)    :: cold  ! number densities at previous timestep (molecule.cm-3) 
real, dimension(nz,nesp) :: c     ! number densities at current timestep (molecule.cm-3) 
real, dimension(nesp)    :: cnew  ! number densities at next timestep (molecule.cm-3) 
      
! timesteps

real :: ctimestep           ! standard timestep for the chemistry (s) 
real :: dt_guess            ! first-guess timestep (s) 
real :: dt_corrected        ! corrected timestep (s) 
real :: time                ! internal time (between 0 and ptimestep, in s)
integer :: phychemrat

! reaction rates

integer, parameter :: nb_phot_max       = 30
integer, parameter :: nb_reaction_3_max = 12
integer, parameter :: nb_reaction_4_max = 87

real, dimension(nz,nb_phot_max)       :: v_phot
real, dimension(nz,nb_reaction_3_max) :: v_3
real, dimension(nz,nb_reaction_4_max) :: v_4

logical,parameter :: hetero_ice  = .false.
logical,parameter :: hetero_dust = .false.

! matrix

real, dimension(nesp,nesp) :: mat, mat1
integer, dimension(nesp)   :: indx
integer                    :: code

! production and loss terms (for first-guess solution only)

real, dimension(nesp) :: prod, loss

! indexes

integer :: i, iesp, iz

if (firstcall) then
!===================================================================
!     read photolysis lookup table
!===================================================================

   call init_chimie(nj, nztable, nsza, nso2, jphot, table_colair, table_colso2, table_sza)

!===================================================================
!     initialisation of the reaction indexes
!===================================================================

   call indice(nb_phot_max, nb_reaction_3_max, nb_reaction_4_max)

   firstcall = .false.
end if

! cloud and dust surfaces set to zero for the moment

surfice1d(:)  = 0.
surfdust1d(:) = 0.
      
!===================================================================
!   number densities (molecule.cm-3)
!===================================================================

do iz = 1,nz
   conc(iz) = p(iz)/(1.38E-19*t(iz))
   c(iz,:) = tr(iz,:)*conc(iz) 
end do
      
!===================================================================
!    photodissociations         
!===================================================================

! dist_sol : sun-venus distance (au)

dist_sol = 0.72333

call phot(nj, nztable, nsza, nso2, sza_input, dist_sol, mumean, tr(:,i_co2), tr(:,i_so2),         &
          jphot, table_colair, table_colso2, table_sza, nz, nb_phot_max, t, p, v_phot)

!===================================================================
!     reaction rates                                     
!===================================================================
                   
call krates(hetero_ice, hetero_dust, nz, nesp, nj, c, conc, t, p, nb_phot_max, nb_reaction_3_max, &
            nb_reaction_4_max, v_3, v_4, v_phot, sza_input)

!===================================================================
!     ctimestep : standard chemical timestep (s), defined as 
!                 the fraction phychemrat of the physical timestep                           
!===================================================================
       
phychemrat = 1

ctimestep  = ptimestep/real(phychemrat)

!===================================================================
!     loop over levels         
!===================================================================

do iz = 1,nz

!  initializations

   time = 0.
   iter(iz) = 0
   dt_guess = ctimestep
   cold(:) = c(iz,:)     

!  internal loop for the chemistry
          
   do while (time < ptimestep)
      
   iter(iz) = iter(iz) + 1

!  first-guess: fill matrix

   call fill_matrix(iz, mat1, prod, loss, c, nesp, nz,                    &
                    nb_reaction_3_max, nb_reaction_4_max, nb_phot_max,    &
                    v_phot, v_3, v_4)

!  adaptative evaluation of the sub time step

   call define_dt(nesp, dt_corrected, dt_guess, ctimestep, cold(:), c(iz,:), &
                  mat1, prod, loss, conc(iz))

   if (time + dt_corrected > ptimestep) then
      dt_corrected = ptimestep - time
   end if

!  form the matrix identity + mat*dt_corrected

   mat(:,:) = mat1(:,:)*dt_corrected
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

!  solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

   cnew(:) = c(iz,:)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   stop
#endif

!  eliminate small values

   where (cnew(:)/conc(iz) < 1.e-30)
      cnew(:) = 0.
   end where

!  update concentrations

   cold(:) = c(iz,:)
   c(iz,:) = cnew(:)
   cnew(:) = 0.

!  increment internal time

   time = time + dt_corrected
   dt_guess = dt_corrected     ! first-guess timestep for next iteration

   end do ! while (time < ptimestep)

!  save mixing ratios

   tr(iz,:)  = max(c(iz,:)/conc(iz), 1.e-30)
                       
end do  ! end of loop over vertical levels
   
!==================
!!!!! MODEL 1D !!!! ==> n_lon = 1 !!!!
!==================

IF(n_lon .EQ. 1) THEN
PRINT*,'On est en 1D'
!PRINT*,"DEBUT rate_save"
CALL rate_save(nz,p(:),t(:),tr(:,:),nesp,v_phot(:,:),v_3(:,:),v_4(:,:))
!PRINT*,"FIN rate_save"
END IF
      
end subroutine photochemistry_venus

!======================================================================

 subroutine init_chimie(nj, nztable, nsza, nso2, jphot, table_colair, &
                        table_colso2, table_sza)

!======================================================================

implicit none

! photolysis lookup table

integer, INTENT(IN) :: nj, nztable, nsza, nso2
real, INTENT(OUT), dimension(nso2,nsza,nztable,nj) :: jphot
real, INTENT(OUT), dimension(nztable) :: table_colair
real, INTENT(OUT), dimension(nso2,nztable) :: table_colso2
real, INTENT(OUT), dimension(nsza) :: table_sza

integer           :: iz, isza, iozo, iso2, ij
character(len=44) :: jvenus

! lecture de la table des j

jphot(:,:,:,:) = 0.

jvenus = 'jvenus.dat'
open(30, form = 'formatted', status = 'old', file = jvenus)
print*,'lecture de jvenus = ', jvenus

do iso2 = 1,nso2
   do isza = 1,nsza
      do iz = nztable,1,-1
         read(30,*) table_colair(iz), table_colso2(iso2,iz), table_sza(isza)
         read(30,'(7e11.4)') (jphot(iso2,isza,iz,ij), ij = 1,nj)
         do ij = 1,nj
            if (jphot(iso2,isza,iz,ij) == 1.E-30) then 
               jphot(iso2,isza,iz,ij) = 0.
            end if
         end do
      end do
   end do
end do

close(30)
print*,'lecture de la table des j ok.'

end subroutine init_chimie

!======================================================================

 subroutine indice(nb_phot_max, nb_reaction_3_max, nb_reaction_4_max)

!================================================================
! set the "indice" arrays used to fill the jacobian matrix      !
!----------------------------------------------------------------
! reaction               type                array              !
!----------------------------------------------------------------
! A + hv   --> B + C     photolysis          indice_phot        ! 
! A + B    --> C + D     bimolecular         indice_4           !
! A + A    --> B + C     quadratic           indice_3           !
! A + C    --> B + C     quenching           indice_phot        !
! A + ice  --> B + C     heterogeneous       indice_phot        !
!================================================================

use types_asis
use chemparam_mod
 
implicit none

! input

integer :: nb_phot_max, nb_reaction_3_max, nb_reaction_4_max

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
!      O3 + hv -> O2(Dg) + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2dg, 1.0, i_o1d)

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
!      HO2 + hv -> OH + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 1.0, i_oh, 1.0, i_o)

!===========================================================
!      H2O2 + hv -> OH + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 2.0, i_oh, 0.0, i_dummy)

!===========================================================
!      HCl + hv -> H + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hcl, 1.0, i_h, 1.0, i_cl)

!===========================================================
!      Cl2 + hv -> Cl + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_cl2, 2.0, i_cl, 0.0, i_dummy)

!===========================================================
!      HOCl + hv -> OH + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hocl, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      SO2 + hv -> SO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so2, 1.0, i_so, 1.0, i_o)

!===========================================================
!      SO + hv -> S + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so, 1.0, i_s, 1.0, i_o)

!===========================================================
!      SO3 + hv -> SO2 + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so3, 1.0, i_so2, 1.0, i_o)

!===========================================================
!      ClO + hv -> Cl + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_clo, 1.0, i_cl, 1.0, i_o)

!===========================================================
!      OCS + hv -> CO + S
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ocs, 1.0, i_co, 1.0, i_s)

!===========================================================
!      COCl2 + hv -> Cl + Cl + CO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_cocl2, 2.0, i_cl, 1.0, i_co)

!===========================================================
!      H2SO4 + hv -> SO3 + H2O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2so4, 1.0, i_so3, 1.0, i_h2o)

!===========================================================
!      a001 : O + O2 + CO2 -> O3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o2, 1.0, i_o3, 0.0, i_dummy)

!===========================================================
!      a002 : O + O + CO2 -> O2(Dg) + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_o, 1.0, i_o2dg, 0.0, i_dummy)

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
!      f001 : HCl + O(1D) -> OH + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o1d, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      f002 : HCl + O(1D) -> H + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o1d, 1.0, i_h, 1.0, i_clo)

!===========================================================
!      f003 : HCl + O -> OH + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      f004 : HCl + OH -> H2O + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_oh, 1.0, i_h2o, 1.0, i_cl)

!===========================================================
!      f005 : ClO + O -> Cl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_o, 1.0, i_cl, 1.0, i_o2)

!===========================================================
!      f006 : ClO + OH -> Cl + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_oh, 1.0, i_cl, 1.0, i_ho2)

!===========================================================
!      f007 : ClO + OH -> HCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_oh, 1.0, i_hcl, 1.0, i_o2)

!===========================================================
!      f008 : Cl + H2 -> HCl + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h2, 1.0, i_hcl, 1.0, i_h)

!===========================================================
!      f009 : Cl + O3 -> ClO + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_o3, 1.0, i_clo, 1.0, i_o2)

!===========================================================
!      f010 : Cl + HO2 -> ClO + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_ho2, 1.0, i_clo, 1.0, i_oh)

!===========================================================
!      f011 : Cl + HO2 -> HCl + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_ho2, 1.0, i_hcl, 1.0, i_o2)

!===========================================================
!      f012 : Cl + H2O2 -> HCl + HO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h2o2, 1.0, i_hcl, 1.0, i_ho2)

!===========================================================
!      f013 : Cl + CO + CO2 -> ClCO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_co, 1.0, i_clco, 0.0, i_dummy)

!===========================================================
!      f014 : ClCO + CO2 -> Cl + CO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_clco, 1.0, i_cl, 1.0, i_co)

!===========================================================
!      f015 : ClCO + O2 + CO2 -> ClCO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o2, 1.0, i_clco3, 0.0, i_dummy)

!===========================================================
!      f016 : 0.5 ClCO3 + 0.5 Cl -> Cl
!             0.5 ClCO3 + 0.5 Cl -> ClO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_cl, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_cl, 1.0, i_clo, 1.0, i_co2)

!===========================================================
!      f017 : 0.5 ClCO3 + 0.5 O -> Cl
!             0.5 ClCO3 + 0.5 O -> O2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_o, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_o, 1.0, i_o2, 1.0, i_co2)

!===========================================================
!      f018 : ClO + HO2 -> HOCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_ho2, 1.0, i_hocl, 1.0, i_o2)

!===========================================================
!      f019 : OH + HOCl -> H2O + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hocl, 1.0, i_h2o, 1.0, i_clo)

!===========================================================
!      f020 : O + HOCl -> OH + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hocl, 1.0, i_oh, 1.0, i_clo)

!===========================================================
!      f021 : Cl + Cl + CO2 -> Cl2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_cl, 1.0, i_cl2, 0.0, i_dummy)

!===========================================================
!      f022 : ClCO + O -> Cl + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o, 1.0, i_cl, 1.0, i_co2)

!===========================================================
!      f023 : Cl2 + O(1D) -> Cl + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_o1d, 1.0, i_cl, 1.0, i_clo)

!===========================================================
!      f024 : Cl2 + H -> HCl + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_h, 1.0, i_hcl, 1.0, i_cl)

!===========================================================
!      f025 : Cl + ClCO -> Cl2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_clco, 1.0, i_cl2, 1.0, i_co)

!===========================================================
!      f026 : ClCO + ClCO -> COCl2 + CO
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_clco, 1.0, i_cocl2, 1.0, i_co)

!===========================================================
!      f027 : Cl + SO2 + CO2 -> ClSO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_so2, 1.0, i_clso2, 0.0, i_dummy)

!===========================================================
!      f028 : ClSO2 + O -> SO2 + ClO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_o, 1.0, i_so2, 1.0, i_clo)

!===========================================================
!      f029 : ClSO2 + H -> SO2 + HCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_h, 1.0, i_so2, 1.0, i_hcl)

!===========================================================
!      f030 : ClSO2 + ClSO2 -> Cl2 + SO2 + SO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_clso2, 1.0, i_cl2, 2.0, i_so2)

!===========================================================
!      f031 : Cl + O + CO2 -> ClO + CO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_o, 1.0, i_clo, 0.0, i_dummy)

!===========================================================
!      f032 : Cl2 + O -> ClO + Cl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_o, 1.0, i_clo, 1.0, i_cl)

!===========================================================
!      f033 : ClCO + OH -> HOCl + CO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_oh, 1.0, i_hocl, 1.0, i_co)

!===========================================================
!      f034 : Cl2 + OH -> Cl + HOCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_oh, 1.0, i_cl, 1.0, i_hocl)

!===========================================================
!      f035 : ClCO + O -> CO + ClO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o, 1.0, i_co, 1.0, i_clo)

!===========================================================
!      f036 : ClCO + Cl2 -> COCl2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_cl2, 1.0, i_cocl2, 1.0, i_cl)

!===========================================================
!      f037 : HCl + H -> H2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_h, 1.0, i_h2, 1.0, i_cl)

!===========================================================
!      f038 : ClCO + H -> HCl + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_h, 1.0, i_hcl, 1.0, i_co)

!===========================================================
!      f039 : Cl + H + M -> HCl + M
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h, 1.0, i_hcl, 0.0, i_dummy)

!===========================================================
!      g001 : S + O2 -> SO + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o2, 1.0, i_so, 1.0, i_o)

!===========================================================
!      g002 : S + O3 -> SO + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o3, 1.0, i_so, 1.0, i_o2)

!===========================================================
!      g003 : SO + O2 -> SO2 + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o2, 1.0, i_so2, 1.0, i_o)

!===========================================================
!      g004 : SO + O3 -> SO2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o3, 1.0, i_so2, 1.0, i_o2)

!===========================================================
!      g005 : SO + OH -> SO2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_oh, 1.0, i_so2, 1.0, i_h)

!===========================================================
!      g006 : S + OH -> SO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_oh, 1.0, i_so, 1.0, i_h)

!===========================================================
!      g007 : SO + O + CO2 -> SO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o, 1.0, i_so2, 0.0, i_dummy)

!===========================================================
!      g008 : SO + HO2 -> SO2 + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_ho2, 1.0, i_so2, 1.0, i_oh)

!===========================================================
!      g009 : SO2 + O + CO2 -> SO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so2, 1.0, i_o, 1.0, i_so3, 0.0, i_dummy)

!===========================================================
!      g010 : S + O + CO2 -> SO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o, 1.0, i_so, 0.0, i_dummy)

!===========================================================
!      g011 : SO3 + H2O -> H2SO4
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_h2o, 1.0, i_h2so4, 0.0, i_dummy)

!===========================================================
!      g012 : SO + ClO -> SO2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_clo, 1.0, i_so2, 1.0, i_cl)

!===========================================================
!      g013 : SO + SO3 -> SO2 + SO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_so3, 2.0, i_so2, 0.0, i_dummy)

!===========================================================
!      g014 : SO3 + O -> SO2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_o, 1.0, i_so2, 1.0, i_o2)

!===========================================================
!      g015 : SO + SO + CO2 -> S2O2 + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_s2o2, 0.0, i_dummy)

!===========================================================
!      g016 : S2O2 + CO2 -> SO + SO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_s2o2, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!      g017 : 0.5 ClCO3 + 0.5 SO -> Cl  
!             0.5 ClCO3 + 0.5 SO -> SO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_so, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_so, 1.0, i_so2, 1.0, i_co2)

!===========================================================
!      g018 : S + CO + CO2 -> OCS + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_co, 1.0, i_ocs, 0.0, i_dummy)

!===========================================================
!      g019 : ClCO + S -> OCS + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_s, 1.0, i_ocs, 1.0, i_cl)

!===========================================================
!      g020 : SO2 + OH + CO2 -> HSO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so2, 1.0, i_oh, 1.0, i_hso3, 0.0, i_dummy)

!===========================================================
!      g021 : HSO3 + O2 -> HO2 + SO3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hso3, 1.0, i_o2, 1.0, i_ho2, 1.0, i_so3)

!===========================================================
!      g022 : S + S + CO2 -> S2 + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_s, 1.0, i_s2, 0.0, i_dummy)

!===========================================================
!      g023 : S2 + hv -> S + S 
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_s2, 2.0, i_s, 0.0, i_dummy)

!===========================================================
!      g024 : S2 + O -> SO + S
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s2, 1.0, i_o, 1.0, i_so, 1.0, i_s)

!===========================================================
!      g025 : S + OCS -> S2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_ocs, 1.0, i_s2, 1.0, i_co)

!===========================================================
!      g026 : OCS + O -> SO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ocs, 1.0, i_o, 1.0, i_so, 1.0, i_co)

!===========================================================
!      g027 : S + SO3 -> SO2 + SO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_so3, 1.0, i_so2, 1.0, i_so)

!===========================================================
!      g028 : S + HO2 -> SO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_ho2, 1.0, i_so, 1.0, i_oh)

!===========================================================
!      g029 : S + ClO -> SO + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_clo, 1.0, i_so, 1.0, i_cl)

!===========================================================
!      g030: h2so4 + h2o -> so3 + h2o + h2o
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2so4, 1.0, i_so3, 1.0, i_h2o)

!===========================================================
!      g031: so3 + ocs -> s2o2 +  co2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_ocs, 1.0, i_s2o2, 1.0, i_co2)

!===========================================================
!      g032: s2o2 + ocs -> co + so2 + s2
!===========================================================
!	decomposee en
!	0.5 s2o2 + 0.5 ocs -> co
!	0.5 s2o2 + 0.5 ocs -> so2 + s2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_s2o2, 0.5, i_ocs, 1.0, i_co, 0.0, i_dummy) 

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_s2o2, 0.5, i_ocs, 1.0, i_so2, 1.0, i_s2) 

!===========================================================
!      g033: so + so -> so2 + s
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_so2, 1.0, i_s)

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
!      i001: O2(Dg) + CO2 -> O2 + CO2 + hv
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2dg, 1.0, i_o2, 0.0, i_dummy)

!===========================================================
!      i002: O2(Dg) -> O2 + hv
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2dg, 1.0, i_o2, 0.0, i_dummy)

!===========================================================
!  check dimensions 
!===========================================================

print*, 'nb_phot       = ', nb_phot
print*, 'nb_reaction_4 = ', nb_reaction_4
print*, 'nb_reaction_3 = ', nb_reaction_3

!print*, 'check dimension'
if ((nb_phot /= nb_phot_max)             .or.  &
    (nb_reaction_3 /= nb_reaction_3_max) .or.  &
    (nb_reaction_4 /= nb_reaction_4_max)) then
  print*, 'wrong dimensions in indice' 
  stop
end if  

end subroutine indice

!===========================================================

subroutine phot(nj, nztable, nsza, nso2, sza_input, dist_sol, mumean, rmco2, rmso2,                     &
                jphot, table_colair, table_colso2, table_sza, nz, nb_phot_max, t, p, v_phot)

!===========================================================

implicit none

integer, INTENT(IN) :: nz
integer, INTENT(IN) :: nj, nztable, nsza, nso2

real, INTENT(IN), dimension(nz) :: t, p
real, INTENT(IN), dimension(nz) :: mumean        ! [g/mol]
real, INTENT(IN), dimension(nso2,nsza,nztable,nj) :: jphot
real, INTENT(IN), dimension(nso2,nztable) :: table_colso2
real, INTENT(IN), dimension(nztable) :: table_colair
real, INTENT(IN), dimension(nsza) :: table_sza
real, INTENT(IN), dimension(nz) :: rmco2, rmso2
real, INTENT(IN) :: sza_input, dist_sol

real, dimension(nz,nj) :: j
real, dimension(nz) :: coef, col, colso2
real, dimension(nso2) :: colref
real, dimension(2,2,2) :: poids
real :: cicol, cisza, ciso2
real :: avogadro, gvenus, dp

integer :: indcol, indsza, indso2
integer :: isza, iz, i, iso2, ij

integer :: nb_phot_max
real, dimension(nz,nb_phot_max), INTENT(INOUT) :: v_phot

!mugaz    = 43.44E-3
avogadro = 6.022E+23
gvenus   = 8.87

! day/night test

if (sza_input <= 95.) then      ! day

! interpolation in solar zenith angle

indsza = nsza - 1
do isza = 1,nsza
   if (table_sza(isza) >= sza_input) then
      indsza = min(indsza,isza - 1)
      indsza = max(indsza, 1)
   end if
end do

cisza = (sza_input - table_sza(indsza))                       &
       /(table_sza(indsza + 1) - table_sza(indsza))

!    print*, 'indsza    = ', indsza
!    print*, 'table_sza = ', table_sza(indsza)
!    print*, 'cisza     = ', cisza

! co2 and so2 columns

coef(nz)   = avogadro/(gvenus*mumean(nz)*1.E-3)*1.E-4
col(nz)    = coef(nz)*rmco2(nz)*p(nz)*100.
colso2(nz) = coef(nz)*rmso2(nz)*p(nz)*100.

do iz = nz-1, 1, -1
!   print*,"L2490 new_photochemistry", iz,mumean(iz)
   dp = (p(iz) - p(iz+1))*100.
   coef(iz)   = avogadro/(gvenus*mumean(iz)*1.E-3)*1.E-4
   col(iz)    = col(iz+1) + coef(iz)*(rmco2(iz+1) + rmco2(iz))*0.5*dp
   col(iz)    = min(col(iz), table_colair(1))
   colso2(iz) = colso2(iz+1) + coef(iz)*(rmso2(iz+1) + rmso2(iz))*0.5*dp
   colso2(iz) = min(colso2(iz), table_colso2(nso2,1))
end do

! loop over altitude

do iz = 1,nz

! interpolation in co2 column

   do i = 1,nztable
      if (table_colair(i) < col(iz)) then
         cicol = (log(col(iz)) - log(table_colair(i)))           &
                /(log(table_colair(i-1)) - log(table_colair(i)))
         indcol = i - 1
         exit
      end if
   end do

! interpolation in so2 column

! initialize indso2 and ciso2 in case colref is never larger
! than the gcm so2 column.

   indso2 = nso2 - 1
   ciso2 = 1.

! search for the index indso2 between which interpolate

   do iso2 = 1,nso2 
      colref(iso2) = cicol*table_colso2(iso2,indcol)               &
                   + (1.-cicol)*table_colso2(iso2,indcol+1)
      if (colref(iso2) > colso2(iz)) then
         ciso2 = (colso2(iz) - colref(iso2-1))                     &
                /(colref(iso2) - colref(iso2-1))
         indso2 = iso2 - 1
         exit
      end if
   end do

! 4-dimensional interpolation weights

! poids(so2,sza_input,co2)

   poids(1,1,1) = (1.-ciso2)*(1.-cisza)*    cicol 
   poids(1,1,2) = (1.-ciso2)*(1.-cisza)*(1.-cicol)
   poids(1,2,1) = (1.-ciso2)*    cisza *    cicol 
   poids(1,2,2) = (1.-ciso2)*    cisza *(1.-cicol)
   poids(2,1,1) =     ciso2 *(1.-cisza)*    cicol 
   poids(2,1,2) =     ciso2 *(1.-cisza)*(1.-cicol)
   poids(2,2,1) =     ciso2 *    cisza *    cicol 
   poids(2,2,2) =     ciso2 *    cisza *(1.-cicol)

! 4-dimensional interpolation in the lookup table

   do ij = 1,nj
      j(iz,ij) =                                          &
      poids(1,1,1)*jphot(indso2  ,indsza  ,indcol  ,ij)   &
    + poids(1,1,2)*jphot(indso2  ,indsza  ,indcol+1,ij)   &
    + poids(1,2,1)*jphot(indso2  ,indsza+1,indcol  ,ij)   &
    + poids(1,2,2)*jphot(indso2  ,indsza+1,indcol+1,ij)   &
    + poids(2,1,1)*jphot(indso2+1,indsza  ,indcol  ,ij)   &
    + poids(2,1,2)*jphot(indso2+1,indsza  ,indcol+1,ij)   &
    + poids(2,2,1)*jphot(indso2+1,indsza+1,indcol  ,ij)   &
    + poids(2,2,2)*jphot(indso2+1,indsza+1,indcol+1,ij)
   end do

end do           ! end of loop over altitude

else             ! night
   j(:,:) = 0.
end if

! photodissociation rates numbering in the lookup table

!    1     o2 + hv     -> o + o
!    2     o2 + hv     -> o + o(1d)
!    3     co2 + hv    -> co + o
!    4     co2 + hv    -> co + o(1d)
!    5     o3 + hv     -> o2(Dg) + o(1d)
!    6     o3 + hv     -> o2 + o
!    7     h2o + hv    -> h + oh
!    8     ho2 + hv    -> oh + o
!    9     h2o2 + hv   -> oh + oh
!    10    hcl + hv    -> h + cl
!    11    cl2 + hv    -> cl + cl
!    12    hocl + hv   -> oh + cl
!    13    so2 + hv    -> so + o
!    14    so + hv     -> s + o
!    15    so3 + hv    -> so2 + o
!    16    clo + hv    -> cl + o
!    17    ocs + hv    -> co + s
!    18    cocl2 + hv  -> cl + cl + co
!    19    h2so4 + hv  -> so3 + h2o

! fill v_phot array

do ij = 1,nj
   v_phot(:,ij) = j(:,ij)
end do

!PRINT*,'sza_input: ',sza_input
!IF (sza_input.le.40.5 .AND. sza_input.gt.39.5) THEN
!open(200, form = 'formatted')
!100    format(e20.6)
!write(200,100)(v_phot(:,19))
!stop
!END IF

end subroutine phot

!======================================================================

 subroutine krates(hetero_ice,hetero_dust, nz, nesp, nj, c, conc, t, p, nb_phot_max, nb_reaction_3_max, nb_reaction_4_max, v_3, v_4, v_phot,sza_input)
 
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


 USE chemparam_mod
implicit none

real, INTENT(IN), dimension(nz)  :: t, p, conc
integer, INTENT(IN) :: nesp, nj, nz
real, INTENT(IN) :: sza_input
integer :: iz
logical, INTENT(IN) :: hetero_ice, hetero_dust
real    :: ak0, ak1, xpo, rate, pi, gam, bid
real    :: k1a0, k1b0, k1ainf, k1a, k1b, fc, fx, x, y
real, dimension(nz)  :: surfice1d, surfdust1d
real, dimension(nz) :: deq
real, dimension(nz) :: a001, a002, a003,                           &
                       b001, b002, b003, b004, b005, b006, b007,   &
                       b008, b009,                                 &
                       c001, c002, c003, c004, c005, c006, c007,   &
                       c008, c009, c010, c011, c012, c013, c014,   &
                       c015, c016, c017, c018,                     &
                       d001, d002, d003,                           &
                       e001, e002, e003, e004, e005, e006, e007,   &
                       e008, e009, e010, e011, e012, e013, e014,   &
                       e015, e016, e017, e018, e019, e020, e021,   &
                       e022, e023, e024, e025, e026, e027, e028,   &
                       e029, e030, e031, e032, e033, e034, e035,   &
                       e036, e037, e038, e039, e040, e041, e042,   &
                       e043,                                       &
                       f001, f002, f003, f004, f005, f006, f007,   &
                       f008, f009, f010, f011, f012, f013, f014,   &
                       f015, f016, f017, f018, f019, f020, f021,   &
                       f022, f023, f024, f025, f026, f027, f028,   &
                       f029, f030, f031, f032, f033, f034, f035,   &
                       f036, f037, f038, f039,                     &
                       g001, g002, g003, g004, g005, g006, g007,   &
                       g008, g009, g010, g011, g012, g013, g014,   &
                       g015, g016, g017, g018, g019, g020, g021,   &
                       g022, g023, g024, g025, g026, g027, g028,   &
                       g029, g030, g031, g032, g033,               &
                       h001, h002, h003,                           &
                       i001, i002

real, INTENT(IN), dimension(nz,nesp) :: c

integer :: nb_phot_max, nb_reaction_3_max, nb_reaction_4_max, nb_reaction_3, nb_reaction_4, nb_phot

real, dimension(nz,nb_phot_max) :: v_phot
real, dimension(nz,nb_reaction_3_max) :: v_3
real, dimension(nz,nb_reaction_4_max) :: v_4

pi = acos(-1.)

nb_phot       = nj
nb_reaction_3 = 0
nb_reaction_4 = 0

!----------------------------------------------------------------------
!        reactions avec ox
!----------------------------------------------------------------------

!---  a001: o + o2 + co2 -> o3 + co2

!     jpl 2003

      a001(:) = 2.5*6.0E-34*(t(:)/300.)**(-2.4)*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a001(:)

!---  a002: o + o + co2 -> o2 + co2

!     Tsang and Hampson, J. Chem. Phys. Ref. Data, 15, 1087, 1986

!     a002(:) = 2.5*5.2E-35*exp(900./t(:))*conc(:)

!     Campbell and Gray, Chem. Phys. Lett., 18, 607, 1973

      a002(:) = 2.5*9.46E-34*exp(485./t(:))*conc(:) ! nist expression

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = a002(:)

!---  a003: o + o3 -> o2 + o2

!     jpl 2003

      a003(:) = 8.0E-12*exp(-2060./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a003(:)

!----------------------------------------------------------------------
!        reactions avec o(1d)
!----------------------------------------------------------------------

!---  b001: o(1d) + co2  -> o + co2

!     jpl 2006

      b001(:) = 7.5E-11*exp(115./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b001(:)*c(:,i_co2)

!---  b002: o(1d) + h2o  -> oh + oh

!     jpl 2006
 
      b002(:) = 1.63E-10*exp(60./t(:))
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b002(:)

!---  b003: o(1d) + h2  -> oh + h

!     jpl 2003
!      b003(:) = 1.1E-10

!      jpl 2016      
      b003(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b003(:)

!---  b004: o(1d) + o2  -> o + o2

!     jpl 2006

      b004(:) = 3.3E-11*exp(55./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b004(:)*c(:,i_o2)
    
!---  b005: o(1d) + o3  -> o2 + o2

!     jpl 2003

      b005(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b005(:)
    
!---  b006: o(1d) + o3  -> o2 + o + o

!     jpl 2003

      b006(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b006(:)
    
!----------------------------------------------------------------------
!        reactions des hox    
!----------------------------------------------------------------------

!---  c001: o + ho2 -> oh + o2

!     jpl 2003

      c001(:) = 3.0E-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c001(:)

!---  c002: o + oh -> o2 + h

!     jpl 2003
!      c002(:) = 2.2E-11*exp(120./t(:))

!     jpl 2016  
      c002(:) = 1.8E-11*exp(180./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c002(:)

!---  c003: h + o3 -> oh + o2

!     jpl 2003

      c003(:) = 1.4E-10*exp(-470./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c003(:)

!---  c004: h + ho2 -> oh + oh

!     jpl 2006

      c004(:) = 7.2E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c004(:)

!---  c005: h + ho2 -> h2 + o2

!     jpl 2006

      c005(:) = 6.9E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c005(:)

!---  c006: h + ho2 -> h2o + o

!     jpl 2006

      c006(:) = 1.6E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c006(:)

!---  c007: oh + ho2 -> h2o + o2

!     jpl 2003

      c007(:) = 4.8E-11*exp(250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c007(:)

!---  c008: ho2 + ho2 -> h2o2 + o2

!     jpl 2006

!     c008(:) = 3.5E-13*exp(430./t(:))

!     christensen et al., grl, 13, 2002
!      c008(:) = 1.5E-12*exp(19./t(:))

!     jpl 2016
      c008(:) = 3.0E-13*exp(460./t(:))

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c008(:)

!---  c009: oh + h2o2 -> h2o + ho2

!     jpl 2006

      c009(:) = 1.8E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c009(:)

!---  c010: oh + h2 -> h2o + h

!     jpl 2006

      c010(:) = 2.8E-12*exp(-1800./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c010(:)

!---  c011: h + o2 + co2 -> ho2 + co2

!     jpl 2006

      do iz = 1,nz
         ak0 = 2.5*4.4E-32*(t(iz)/300.)**(-1.3)
         ak1 = 4.7E-11*(t(iz)/300.)**(-0.2)

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         c011(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c011(:)

!---  c012: o + h2o2 -> oh + ho2

!     jpl 2003

      c012(:) = 1.4E-12*exp(-2000./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c012(:)

!---  c013: oh + oh -> h2o + o

!     jpl 2006

      c013(:) = 1.8E-12

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c013(:)

!---  c014: oh + o3 -> ho2 + o2

!     jpl 2003

      c014(:) = 1.7E-12*exp(-940./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c014(:)

!---  c015: ho2 + o3 -> oh + o2 + o2

!     jpl 2003

      c015(:) = 1.0E-14*exp(-490./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c015(:)

!---  c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2

!     jpl 2003
!      c016(:) = 2.5*1.7E-33*exp(1000./t(:))*conc(:)

!     jpl 2016
      c016(:) = 2.5*2.1E-33*exp(920./t(:))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c016(:)

!---  c017: oh + oh + co2 -> h2o2 + co2

!     jpl 2003

      do iz = 1,nz
         ak0 = 2.5*6.9E-31*(t(iz)/300.)**(-1.0)
         ak1 = 2.6E-11*(t(iz)/300.)**(0.0)

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         c017(iz) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c017(:)

!---  c018: h + h + co2 -> h2 + co2

!     baulch et al., 2005

      c018(:) = 2.5*1.8E-30*(t(:)**(-1.0))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c018(:)

!----------------------------------------------------------------------
!        reactions des composes azotes
!----------------------------------------------------------------------

!---  d001: no2 + o -> no + o2

!     jpl 2006

      d001(:) = 5.1E-12*exp(210./t(:))

!---  d002: no + o3 -> no2 + o2

!     jpl 2006

      d002(:) = 3.0E-12*exp(-1500./t(:))

!---  d003: no + ho2 -> no2 + oh

!     jpl 2006

      d003(:) = 3.5E-12*exp(250./t(:))

!----------------------------------------------------------------------
!        reactions des composes carbones
!----------------------------------------------------------------------

!---  e001: oh + co -> co2 + h

!     jpl 2003

!     e001(:) = 1.5E-13*(1 + 0.6*p(:)/1013.)

!     mccabe et al., grl, 28, 3135, 2001

!     e001(:) = 1.57E-13 + 3.54E-33*conc(:)

!     jpl 2006

!     ak0 = 1.5E-13*(t(:)/300.)**(0.6)
!     ak1 = 2.1E-9*(t(:)/300.)**(6.1)
!     rate1 = ak0/(1. + ak0/(ak1/conc(:)))
!     xpo1 = 1./(1. + alog10(ak0/(ak1/conc(:)))**2)

!     ak0 = 5.9E-33*(t(:)/300.)**(-1.4)
!     ak1 = 1.1E-12*(t(:)/300.)**(1.3)
!     rate2 = (ak0*conc(:))/(1. + ak0*conc(:)/ak1)
!     xpo2 = 1./(1. + alog10((ak0*conc(:))/ak1)**2)

!     e001(:) = rate1*0.6**xpo1 + rate2*0.6**xpo2

!     joshi et al., 2006

      do iz = 1,nz
         k1a0 = 1.34*2.5*conc(iz)                                &
               *1/(1/(3.62E-26*t(iz)**(-2.739)*exp(-20./t(iz)))  &
               + 1/(6.48E-33*t(iz)**(0.14)*exp(-57./t(iz))))     ! corrige de l'erreur publi
         k1b0 = 1.17E-19*t(iz)**(2.053)*exp(139./t(iz))          &
              + 9.56E-12*t(iz)**(-0.664)*exp(-167./t(iz))
         k1ainf = 1.52E-17*t(iz)**(1.858)*exp(28.8/t(iz))        &
                + 4.78E-8*t(iz)**(-1.851)*exp(-318./t(iz))
         x = k1a0/(k1ainf - k1b0)
         y = k1b0/(k1ainf - k1b0)
         fc = 0.628*exp(-1223./t(iz)) + (1. - 0.628)*exp(-39./t(iz))  &
            + exp(-t(iz)/255.)
         fx = fc**(1./(1. + (alog(x))**2))                   ! corrige de l'erreur publi
         k1a = k1a0*((1. + y)/(1. + x))*fx
         k1b = k1b0*(1./(1.+x))*fx

         e001(iz) = k1a + k1b
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e001(:)

!---  e002: o + co + m -> co2 + m

!     tsang and hampson, 1986.

      e002(:) = 2.5*6.5E-33*exp(-2184./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e002(:)

!----------------------------------------------------------------------
!        reactions des composes chlores
!----------------------------------------------------------------------

!---  f001: hcl + o(1d) -> oh + cl

!     jpl 2011

      f001(:) = 1.0E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f001(:)

!---  f002: hcl + o(1d) -> h + clo

!     jpl 2011

      f002(:) = 3.6E-11
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f002(:)

!---  f003: hcl + o -> oh + cl

!     jpl 2006

      f003(:) = 1.0E-11*exp(-3300./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f003(:)

!---  f004: hcl + oh -> h2o + cl

!     jpl 2006
!      f004(:) = 2.6E-12*exp(-350./t(:))

!     jpl 2016
      f004(:) = 1.8E-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f004(:)

!---  f005: clo + o -> cl + o2

!     jpl 2006

      f005(:) = 2.8E-11*exp(85./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f005(:)

!---  f006: clo + oh -> cl + ho2

!     jpl 2006

      f006(:) = 7.4E-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f006(:)

!---  f007: clo + oh -> hcl + o2

!     jpl 2006

      f007(:) = 6.0E-13*exp(230./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f007(:)

!---  f008: cl + h2 -> hcl + h

!     jpl 2006

      f008(:) = 3.05E-11*exp(-2270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f008(:)

!---  f009: cl + o3 -> clo + o2

!     jpl 2006

      f009(:) = 2.3E-11*exp(-200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f009(:)

!---  f010: cl + ho2 -> clo + oh

!     jpl 2006
!      f010(:) = 4.1E-11*exp(-450./t(:))

!     jpl 2016
      f010(:) = 3.6E-11*exp(-375./t(:)) 
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f010(:)

!---  f011: cl + ho2 -> hcl + o2

!     jpl 2006
!      f011(:) = 1.8E-11*exp(170./t(:))

!     jpl 2016
      f011(:) = 1.4E-11*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f011(:)

!---  f012: cl + h2o2 -> hcl + ho2

!     jpl 2006

      f012(:) = 1.1E-11*exp(-980./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f012(:)

!---  f013: cl + co + co2 -> clco + co2

!     jpl 2011 + nicovich et al., j. phys. chem., 1990

      f013(:) = 3.2*1.3E-33*(t(:)/300.)**(-3.8)*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f013(:)

!---  f014: clco + co2 -> cl + co + co2

!     jpl 2011

!     deq(:) = 3.2*3.5E-25*exp(3730./t(:))

!     mills, 1998

      deq(:) = 1.6E-25*exp(4000./t(:))

      f014(:) = f013(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = f014(:)*conc(:)

!     do iz = 1, nz
!     print*, z(iz), t(iz), f013(iz), f014(iz), v_phot(iz,nb_phot)
!     end do
!     stop

!---  f015: clco + o2 + m -> clco3 + m

!     yung and demore, icarus, 51, 199-247, 1982.

      f015(:) = 5.7E-15*exp(500./t(:))*conc(:)   &
               /(1.e17 + 0.05*conc(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f015(:)

!---  f016: clco3 + cl -> cl + clo + co2

!     yung and demore, icarus, 51, 199-247, 1982.

!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2

      f016(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f016(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f016(:)

!---  f017: clco3 + o -> cl + o2 + co2

!     yung and demore, icarus, 51, 199-247, 1982.

!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2

      f017(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f017(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f017(:)

!---  f018: clo + ho2  -> hocl + o2

      f018(:) = 2.7E-12*exp(220./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f018(:)

!---  f019: oh + hocl -> h2o + clo

      f019(:) = 3.0E-12*exp(-500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f019(:)

!---  f020: o + hocl -> oh + clo

      f020(:) = 1.7E-13

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f020(:)

!---  f021: cl + cl + co2 -> cl2 + co2

!     donohoue et al., j. phys. chem. a, 109, 7732-7741, 2005

!     f021(:) = 2.5*8.4E-33*exp(850.*(1./t(:) - 1./298.))*conc(:)

!     valeur utilisee par Zhang et al., 2011:

      f021(:) = 2.6E-33*exp(900./t(:))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f021(:)

!---  f022: clco + o -> cl + co2

!     yung et al., icarus, 1982 (estimated)

      f022(:) = 3.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f022(:)

!---  f023: cl2 + o(1d) -> cl + clo

!     jpl 2011

      f023(:) = 2.0E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f023(:)

!---  f024: cl2 + h  -> hcl + cl

!     baulch et al., j. phys. chem. ref. data, 1981

      f024(:) = 1.43E-10*exp(-591./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f024(:)

!---  f025: cl + clco  -> cl2 + co

!     baulch et al., j. phys. chem. ref. data, 1981

      f025(:) = 2.16E-9*exp(-1670./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f025(:)

!---  f026: clco + clco  -> cocl2 + co

!     zhang et al., icarus, 2011 (estimated)

      f026(:) = 5.0E-11

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f026(:)

!---  f027: cl + so2 + co2  -> clso2 + co2

!     mills, phd, 1998

      f027(:) = 1.3E-34*exp(940./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f027(:)

!---  f028: clso2 + o  -> so2 + clo

!     mills, phd, 1998

      f028(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f028(:)

!---  f029: clso2 + h  -> so2 + hcl

!     mills, phd, 1998

      f029(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f029(:)

!---  f030: clso2 + clso2  -> cl2 + so2 + so2

!     moses et al. 2002

      f030(:) = 5.0E-13

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f030(:)

!---  f031: cl + o + co2  -> clo + co2

!     yung and demore, 1999 (estimated)

      f031(:) = 5.0E-32*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f031(:)

!---  f032: cl2 + o -> clo + cl

!     mills, phd, 1998

      f032(:) = 7.4E-12*exp(-1650./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f032(:)

!---  f033: clco + oh -> hocl + co

!     mills, phd, 1998

      f033(:) = 1.5E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f033(:)

!---  f034: cl2 + oh -> cl + hocl

!     jpl 2011

      f034(:) = 2.6E-12*exp(-1100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f034(:)

!---  f035: clco + o -> co + clo

!     yung and demore, 1982

      f035(:) = 3.0E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f035(:)

!---  f036: clco + cl2 -> cocl2 + cl

!     ohta, bull. chem. soc. jpn., 1983

      f036(:) = 6.45E-2*f015(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f036(:)

!---  f037: hcl + h -> h2 + cl

!     mills, phd, 1998

      f037(:) = 1.5E-11*exp(-1750./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f037(:)

!---  f038: clco + h -> hcl + co

!     yung and demore, 1982

      f038(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f038(:)

!---  f039: cl + h + m -> hcl + m

!     yung and demore, 1982 (estimate)

      f039(:) = 1.0E-32*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f039(:)

!----------------------------------------------------------------------
!        reactions des composes soufres
!----------------------------------------------------------------------

!---  g001: s + o2 -> so + o

!      g001(:) = 2.3E-12

!     jpl 2016
      g001(:) = 1.6E-12*exp(100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g001(:)

!---  g002: s + o3 -> so + o2

      g002(:) = 1.2E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g002(:)

!---  g003: so + o2 -> so2 + o

!      g003(:) = 1.25E-13*exp(-2190./t(:))

!     jpl 2016
       g003(:) = 1.6E-13*exp(-2280./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g003(:)

!---  g004: so + o3 -> so2 + o2

      g004(:) = 3.4E-12*exp(-1100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g004(:)

!---  g005: so + oh -> so2 + h

      g005(:) = 2.7E-11*exp(335./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g005(:)

!---  g006: s + oh -> so + h

      g006(:) = 6.6E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g006(:)

!---  g007: so + o + co2 -> so2 + co2

!     singleton and cvetanovic, j. phys. chem. ref. data, 1988
!     measured with co2 as third body

      do iz = 1,nz
         ak0 = 4.2E-30
         ak1 = 5.3E-11

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g007(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g007(:)

!---  g008: so + ho2 -> so2 + oh 

      g008(:) = 2.8E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g008(:)

!---  g009: so2 + o + co2 -> so3 + co2

!     jpl 2011
!     Naido 2005

!      do iz = 1,nz
!         ak0 = 2.5*1.8E-33*(t(iz)/300.)**(2.0)
!         ak1 = 4.2E-14*(t(iz)/300.)**(1.8)

!         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
!         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
!         g009(iz) = rate*0.6**xpo
!         g009(iz) = 0.0E+0
!      end do

      do iz = 1,nz
         ak0 = 5.*9.5*1.E-23*(t(iz)**(-3.0))*EXP(-2400./t(iz))
         ak1 = 6.1*1.E-13*EXP(-850./t(iz))
         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         fc = 0.558*EXP(-t(iz)/316.)+0.442*EXP(-t(iz)/7442.)
         g009(iz) = rate*fc**xpo
!         g009(iz) = 0.0E+0
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g009(:)

!---  g010: s + o + co2 -> so + co2

!     zhang et al., icarus, 2011

      g010(:) = 1.5E-34*exp(900./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g010(:)

!---  g011: so3 + h2o + M -> h2so4 + M
!---  avec M = h2o

      DO iz=1,nz
!     jpl 2011
!      g011(:) = 8.5E-21*exp(6540./t(:))*c(:,i_h2o)
      g011(iz) = 2.26E-23*MAX(t(iz),100.)*exp(6540./MAX(t(iz),100.)) &
                *c(iz,i_h2o)
      g011(iz) = g011(iz)*1.0E-20
!      g011(:) = 0. ! SANS H2SO4
      ENDDO

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g011(:)

!---  g012: so + clo -> so2 + cl 

!     jpl 2011

      g012(:) = 2.8E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g012(:)

!---  g013: so + so3 -> so2 + so2 

!     chung et al., int. j. chem. kinet., 1975

      g013(:) = 2.0E-15

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g013(:)

!---  g014: so3 + o -> so2 + o2 

!     jacob and winkler, j. chem. soc. faraday trans. 1, 1972

      g014(:) = 2.32E-16*exp(-487./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g014(:)

!---  g015: so + so + co2 -> s2o2 + co2

!     herron and huie, chem. phys. lett., 1980.

      do iz = 1,nz
         ak0 = 2.5*4.4E-31
         ak1 = 1.0E-11

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g015(iz) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g015(:)

!---  g016: s2o2 + co2 -> so + so + co2

!     mills, phd, 1998

      deq(:) = 2.5*1.0E-28*exp(6000./t(:))

      g016(:) = g015(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g016(:)*conc(:)

!---  g017: clco3 + so -> cl + so2 + co2

!     mills, phd, 1998

!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2

      g017(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g017(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g017(:)

!---  g018: s + co + co2 -> ocs + co2

!     zhang et al., icarus, 2011 (estimate?)

      g018(:) = 2.5*4.0E-33*exp(-1940./t(:))*conc(:)
      
!      g018(:) = 0.0E+0

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g018(:)

!---  g019: clco + s -> ocs + cl 

!     zhang et al., icarus, 2011

      g019(:) = 3.0E-12

!      g019(:) = 0.0E+0

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g019(:)

!---  g020: so2 + oh + co2 -> hso3 + co2

!     jpl 2011

      do iz = 1,nz
         ak0 = 2.5*3.3E-31*(t(iz)/300.)**(-4.3)
         ak1 = 1.6E-12

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g020(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g020(:)

!---  g021: hso3 + o2 -> ho2 + so3 

!     jpl 2011

      g021(:) = 1.3E-12*exp(-330./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g021(:)

!---  g022: s + s + co2 -> s2 + co2

!     nicholas et al., j. chem. soc. faraday trans. 1, 1979

      do iz = 1,nz
         ak0 = 1.19E-29
         ak1 = 1.0E-10

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g022(iz) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g022(:)

!---  g023: s2 + co2 -> s + s + co2

!     chase et al., 1985

!      deq(:) = 2.68E-25*exp(50860./t(:))

!      g023(:) = g022(:)/(deq(:)*conc(:))

!      nb_phot = nb_phot + 1
!      v_phot(:,nb_phot) = g023(:)*conc(:)

!     Changement pour g023 ==> s2 + hv -> s + s
!     Pas encore inclu dans la table jphot

!     Pas de photodissociation sous le nuage photochimique 
      g023(1:28)  = 0.0E+0

!     Dependance en sza pour la photodissociation
!     moins de W.m-2      
      IF (sza_input.LT.90.0) THEN
        g023(29:50) = 6.5E-3*COS(sza_input*pi/180.0)
      ELSE
        g023(29:50) = 0.0E+0
      END IF
           
      nb_phot = nb_phot + 1
      
      v_phot(:,nb_phot) = g023(:)

!---  g024: s2 + o -> so + s 

!     zhang et al., icarus, 2011

      g024(:) = 2.2E-11*exp(-84./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g024(:)

!---  g025: s + ocs -> s2 +  co

!     lu et al., j. chem. phys., 2006

      g025(:) = 6.63E-20*(t(:)**2.57)*exp(-1180./t(:))
      
!      g025(:) = 0.0E+0

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g025(:)

!---  g026: ocs + o -> so + co

!     atkinson et al., 2004

      g026(:) = 1.60E-11*exp(-2150./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g026(:)

!---  g027: s + so3 -> so2 +  so

!     moses et al., 2002

      g027(:) = 1.0E-16

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g027(:)

!---  g028: s + ho2 -> so +  oh

!     yung and demore, 1982

      g028(:) = 3.0E-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g028(:)

!---  g029: s + clo -> so +  cl

!     moses et al., 2002

      g029(:) = 4.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g029(:)

!---  g030: h2so4 + h2o -> so3 + h2o + h2o 

!     krasnopolsky , 2007 

      g030(:) = 7.0E-14*exp(-5170./t(:))

!      g030(:) = g011(:)/(deq(:)*c(:,i_h2o))
!      g030(:) = 0.0E+0

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g030(:)*c(:,i_h2o)
 !     v_phot(:,nb_phot) = 0.0E+0
      
!---  g031: so3 + ocs -> s2o2 +  co2

!     krasnopolsky , 2007 

      g031(:) = 1.0E-11*exp(-10000./t(:))
!      g031(:) = 0.0E+0
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g031(:)
      
!---  g032: s2o2 + ocs -> co + so2 + s2

!	decomposee en
!	0.5 s2o2 + 0.5 ocs -> co
!	0.5 s2o2 + 0.5 ocs -> so2 + s2

!     krasnopolsky , 2007 

      g032(:) = 1.0E-20
!      g032(:) = 0.0E+0
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g032(:)
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g032(:)

!      g033: so + so -> so2 + s
!  Krasnopolsky 2012 from Martinez & Heron 1983 or Moses et al 2002

!      g033(:) = 3.5E-15
      g033(:) =1.0E-12*exp(-1700.0/t(:))
      
      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g033(:)
            
!----------------------------------------------------------------------
!     heterogeneous chemistry 
!----------------------------------------------------------------------

      if (hetero_ice) then

!        k = (surface*v*gamma)/4 (s-1)
!        v = 100*sqrt(8rt/(pi*m))  (cm s-1)
 
!---     h001: ho2 + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.025
      
         gam = 0.025
         h001(:) = surfice1d(:)*1.E-8       &
                   *100.*sqrt(8.*8.31*t(:)/(33.E-3*pi))*gam/4.
 
!        h002: oh + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.03
 
         gam = 0.03
         h002(:) = surfice1d(:)*1.E-8       &
                   *100.*sqrt(8.*8.31*t(:)/(17.E-3*pi))*gam/4.

!---     h003: h2o2 + ice -> products
 
!        gamma = 0.    test value
 
         gam = 0.
         h003(:) = surfice1d(:)*1.E-8        &
                   *100.*sqrt(8.*8.31*t(:)/(34.E-3*pi))*gam/4.
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

!     do iz = 1,nz
!        print*, z(iz), surfice1d(iz), h001(iz), h002(iz)
!     end do
!     stop

!     print*, 'krates : nb_phot       = ', nb_phot
!     print*, 'krates : nb_reaction_4 = ', nb_reaction_4
!     print*, 'krates : nb_reaction_3 = ', nb_reaction_3
!     stop

!----------------------------------------------------------------------
!     reactions avec 02(Dg) 
!----------------------------------------------------------------------

!---     i001: O2(Dg) + CO2 -> O2 + CO2 + hv

!        Krasnopolsky (2010a)

      i001(:) = 1.E-20

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = i001(:)*c(:,i_co2)

!---     i002: O2(Dg) -> O2 + hv

!        Lafferty et al; (1998)

      i002(:) = 2.2E-4

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = i002(:)

return
end subroutine krates

!======================================================================

 subroutine fill_matrix(ilev, mat, prod, loss, c, nesp, nlayer,            &
                        nb_reaction_3_max, nb_reaction_4_max, nb_phot_max, &
                        v_phot, v_3, v_4)

!======================================================================
! filling of the jacobian matrix
!======================================================================

use types_asis

implicit none

! input

integer             :: ilev    ! level index
integer             :: nesp    ! number of species in the chemistry
integer, intent(in) :: nlayer  ! number of atmospheric layers
integer, intent(in) :: nb_reaction_3_max 
                               ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max
                               ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max
                               ! number of processes treated numerically as photodissociations

real, dimension(nlayer,nesp)              :: c    ! number densities
real, dimension(nlayer,      nb_phot_max) :: v_phot
real, dimension(nlayer,nb_reaction_3_max) :: v_3
real, dimension(nlayer,nb_reaction_4_max) :: v_4

! output

real, dimension(nesp,nesp), intent(out) :: mat  ! matrix
real, dimension(nesp), intent(out)      :: prod, loss

! local

integer :: iesp
integer :: ind_phot_2,ind_phot_4,ind_phot_6
integer :: ind_3_2,ind_3_4,ind_3_6
integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
integer :: iphot,i3,i4

real :: eps, eps_4  ! implicit/explicit coefficient

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
  eps_4 = min(eps_4,1.0)

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
real, dimension(nesp)      :: cold, ccur
real, dimension(nesp,nesp) :: mat1
real, dimension(nesp)      :: prod, loss
real                       :: dens

! output

real :: dtnew

! local

real, dimension(nesp)      :: cnew
real, dimension(nesp,nesp) :: mat
real :: atol, ratio, e, es, coef

integer                  :: code, iesp, iter
integer, dimension(nesp) :: indx
integer :: imax

real :: dttest

! parameters

real, parameter    :: dtmin   = 10.      ! minimum time step (s)
real, parameter    :: vmrtol  = 1.e-11   ! absolute tolerance on vmr
real, parameter    :: rtol    = 0.05     ! rtol recommended value : 0.1-0.02
integer, parameter :: niter   = 3        ! number of iterations
real, parameter    :: coefmax = 2.
real, parameter    :: coefmin = 0.1 
logical            :: fast_guess = .true.

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
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   stop
#endif

end if

! ratio old/new subtimestep

ratio = dtold/dttest

! e : local error indicator

e = 0.

do iesp = 1,nesp
   es = 2.*abs((ratio*cnew(iesp) - (1. + ratio)*ccur(iesp) + cold(iesp))   &
         /(1. + ratio)/max(ccur(iesp)*rtol,atol))

   if (es > e) then
      e = es
      imax = iesp
   end if
end do

! timestep correction

coef = max(coefmin, min(coefmax,0.8/sqrt(e)))

dttest = max(dtmin,dttest*coef)
dttest = min(ctimestep,dttest)

end do ! iter

! new timestep

dtnew = dttest

end subroutine define_dt

!======================================================================

      SUBROUTINE  rate_save(            &
                           n_lev,       &
                           pres,        &
                           temperature, &
                           traceur,     &
                           nq_max,      &
                           vphot,       &
                           v3,          &
                           v4)      
!==================
!!!!! MODEL 1D !!!! ==> n_lon = 1 !!!!
!==================
! Ici on a les variables pour le modele 1D, surtout pour la sauvegarde des taux de prod/consom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PENSER a changer les conditions de time_tot
!time_tot=nbr_pdt*(nbr_jour-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE chemparam_mod
      IMPLICIT none
            

!INTEGER, PARAMETER :: time_tot=6000*1

INTEGER :: unit_loc, ierr_loc           ! unite de lecture de "rcm1d.def"
      
INTEGER, SAVE :: time_tot,nbr_pdt,nbr_jour
INTEGER, SAVE :: cpt_time, cpt_time_rate
DOUBLE PRECISION, DIMENSION(n_lev,126) :: rate_day
DOUBLE PRECISION, DIMENSION(n_lev,126) :: rate_night
DOUBLE PRECISION :: rate_local
DOUBLE PRECISION :: concentration(n_lev)
DOUBLE PRECISION :: pres(n_lev)
DOUBLE PRECISION :: temperature(n_lev)
DOUBLE PRECISION :: traceur(n_lev,nq_max)
      
INTEGER :: n_lev, nq_max
INTEGER :: i_lev, i_react, i_v

INTEGER :: i
 
LOGICAL, SAVE :: f_call = .true.

integer, parameter :: nb_phot_max = 31
integer, parameter :: nb_reaction_3_max = 11
integer, parameter :: nb_reaction_4_max = 84
      
real, dimension(n_lev,nb_phot_max) :: vphot
real, dimension(n_lev,nb_reaction_3_max) :: v3
real, dimension(n_lev,nb_reaction_4_max) :: v4

!PRINT*,"DEBUT subroutine rate_save" 


      IF (f_call) THEN
! ------------------------------------------------------
!  Lecture des parametres dans "rcm1d.def" 
! ------------------------------------------------------

!   Opening parameters file "rcm1d.def"
!   ---------------------------------------
      unit_loc =98
      OPEN(unit_loc,file='rcm1d.def',status='old',form='formatted'  &
          ,iostat=ierr_loc)

      IF(ierr_loc.ne.0) THEN
        write(*,*) 'Problem to open "rcm1d.def'
        write(*,*) 'Is it there ?'
        stop
      ELSE
        write(*,*) 'open rcm1d.def success '
      END IF

      do i=1, 2
        read (unit_loc, *)
      end do

      PRINT *,'nombre de pas de temps par jour ?'
      READ(unit_loc,*) nbr_pdt
      print*,nbr_pdt

      PRINT *,'nombre de jours simules ?'
      READ(unit_loc,*) nbr_jour
      print*,nbr_jour
      
 
      
      time_tot = nbr_pdt*(nbr_jour-1)
      PRINT *,'nombre de PdT avant calcul des taux production/consommation ?'
      PRINT*,time_tot
      
      PRINT*,'nlev',n_lev
           
         cpt_time = 1
         cpt_time_rate = 1
         f_call = .false.
         PRINT*,"f_call: ",f_call
         rate_night(:,:)=0.
         rate_day(:,:)=0.
      
      END IF           
      
!      PRINT*,"P	T"
!      PRINT*,pres,temperature
 	     
      IF (cpt_time .GE. time_tot) THEN

! 	PRINT*,'cpt_time',cpt_time
 	      
         DO i_lev=1, n_lev
         concentration(i_lev) = pres(i_lev)/(1.3806488E-19 * temperature(i_lev))     
         END DO
         
         IF (((cpt_time_rate .GE. 1).AND.(cpt_time_rate .LE. (nbr_pdt/4))).OR. &
         (cpt_time_rate .GT. (3*(nbr_pdt/4)))) THEN
         
!===============================
!        !!!! NUIT !!!!
!===============================
!	PRINT*,'NUIT'
	
           DO i_lev=1, n_lev
           i_react=1
           i_v=1
!===============================
!    1     o2 + hv     -> o + o
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    2     o2 + hv     -> o + o(1d)
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    3     co2 + hv    -> co + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    4     co2 + hv    -> co + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    5     o3 + hv     -> o2(Dg) + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    6     o3 + hv     -> o2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    7     h2o + hv    -> h + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    8     ho2 + hv    -> oh + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    9     h2o2 + hv   -> oh + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    10    hcl + hv    -> h + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    11    cl2 + hv    -> cl + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    12    hocl + hv   -> oh + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    13    so2 + hv    -> so + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    14    so + hv     -> s + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    15    so3 + hv    -> so2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    16    clo + hv    -> cl + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    17    ocs + hv    -> co + s
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ocs)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    18    cocl2 + hv  -> cl + cl + co
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cocl2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    19    h2so4 + hv  -> so3 + h2o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2so4)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 20 b001 o(1d) + co2 -> o + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 21 b004 o(1d) + o2 -> o + o2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 22 f014 clco + co2 -> cl + co + co2
!===============================
            rate_local = vphot(i_lev,22)*traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 23 g016 s2o2 + co2 -> 2so + co2
!===============================
            rate_local = vphot(i_lev,23)*traceur(i_lev,i_s2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 24 g023 s2 + co2 -> 2s + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 25 h001 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 26 h002 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!--- 27 h003 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!---  30 i001 o2(Dg) + CO2 -> O2 + CO2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 31 i002 o2(Dg) -> O2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1

! DEBUT DES REACTION V3
		i_v = i_v - nb_phot_max
!===============================
!--- 32 a002: o + o + co2 -> o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 33 c008: ho2 + ho2 -> h2o2 + o2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 34 c013: oh + oh -> h2o + o
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 35 c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 36 c017: oh + oh + co2 -> h2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 37 c018: h + h + co2 -> h2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 38 f021: cl + cl + co2 -> cl2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 39 f026: clco + clco  -> cocl2 + co
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 40 f030: clso2 + clso2  -> cl2 + so2 + so2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_clso2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 41 g015: so + so + co2 -> s2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 42 g022: s + s + co2 -> s2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 

! DEBUT DES REACTION V4

		i_v = i_v - nb_reaction_3_max

!===============================
!--- 43 a001: o + o2 + co2 -> o3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 44 a003: o + o3 -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 45 b002: o(1d) + h2o  -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 46 b003: o(1d) + h2  -> oh + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 47 b005: o(1d) + o3  -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 48 b006: o(1d) + o3  -> o2 + o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 49 c001: o + ho2 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 50 c002: o + oh -> o2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 51 c003: h + o3 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 52 c004: h + ho2 -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 53 c005: h + ho2 -> h2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 54 c006: h + ho2 -> h2o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 55 c007: oh + ho2 -> h2o + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 56 c009: oh + h2o2 -> h2o + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 57 c010: oh + h2 -> h2o + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 58 c011: h + o2 + co2 -> ho2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 59 c012: o + h2o2 -> oh + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 60 c014: oh + o3 -> ho2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 61 c015: ho2 + o3 -> oh + o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 62 e001: oh + co -> co2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 63 e002: o + co + m -> co2 + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 64 f001: hcl + o(1d) -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 65 f002: hcl + o(1d) -> h + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 66 f003: hcl + o -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 67 f004: hcl + oh -> h2o + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 68 f005: clo + o -> cl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 69 f006: clo + oh -> cl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 70 f007: clo + oh -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 71 f008: cl + h2 -> hcl + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 72 f009: cl + o3 -> clo + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 73 f010: cl + ho2 -> clo + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 74 f011: cl + ho2 -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 75 f012: cl + h2o2 -> hcl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 76 f013: cl + co + co2 -> clco + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 77 f015: clco + o2 + m -> clco3 + m
!===============================
		rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 78 & 79 f016: clco3 + cl -> cl + clo + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 80 & 81 f017: clco3 + o -> cl + o2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 82 f018: clo + ho2  -> hocl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 83 f019: oh + hocl -> h2o + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_hocl)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 84 f020: o + hocl -> oh + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 85 f022: clco + o -> cl + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 86 f023: cl2 + o(1d) -> cl + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 87 f024: cl2 + h  -> hcl + cl
!==============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 88 f025: cl + clco  -> cl2 + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 89 f027: cl + so2 + co2  -> clso2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 90 f028: clso2 + o  -> so2 + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 91 f029: clso2 + h  -> so2 + hcl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 92 f031: cl + o + co2  -> clo + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 93 f032: cl2 + o -> clo + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 94 f033: clco + oh -> hocl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 95 f034: cl2 + oh -> cl + hocl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 96 f035: clco + o -> co + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 97 f036: clco + cl2 -> cocl2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 98 f037: hcl + h -> h2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 99 f038: clco + h -> hcl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 100 f039: cl + h + m -> hcl + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 101 g001: s + o2 -> so + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 102 g002: s + o3 -> so + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 103 g003: so + o2 -> so2 + o
!===============================
             rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 104 g004: so + o3 -> so2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 105 g005: so + oh -> so2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 106 g006: s + oh -> so + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 107 g007: so + o + co2 -> so2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 108 g008: so + ho2 -> so2 + oh 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 109 g009: so2 + o + co2 -> so3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 110 g010: s + o + co2 -> so + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 111 g011: so3 + h2o -> h2so4
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 112 g012: so + clo -> so2 + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 113 g013: so + so3 -> so2 + so2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 114 g014: so3 + o -> so2 + o2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 115 & 116 g017: clco3 + so -> cl + so2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 117 g018: s + co + co2 -> ocs + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 118 g019: clco + s -> ocs + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 119 g020: so2 + oh + co2 -> hso3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 120 g021: hso3 + o2 -> ho2 + so3 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_hso3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 121 g024: s2 + o -> so + s 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 122 g025: s + ocs -> s2 +  co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 123 g026: ocs + o -> so + co
!=============================== 
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 124 g027: s + so3 -> so2 +  so
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 125 g028: s + ho2 -> so +  oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 126 g029: s + clo -> so +  cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1   
           END DO       
         ELSE
!===============================         
!        !!!! JOUR !!!!
!===============================
!	PRINT*,'JOUR'
	
           DO i_lev=1, n_lev
           i_react=1
           i_v=1
!===============================
!    1     o2 + hv     -> o + o
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    2     o2 + hv     -> o + o(1d)
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    3     co2 + hv    -> co + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    4     co2 + hv    -> co + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    5     o3 + hv     -> o2(Dg) + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    6     o3 + hv     -> o2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    7     h2o + hv    -> h + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    8     ho2 + hv    -> oh + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    9     h2o2 + hv   -> oh + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    10    hcl + hv    -> h + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    11    cl2 + hv    -> cl + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    12    hocl + hv   -> oh + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    13    so2 + hv    -> so + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    14    so + hv     -> s + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    15    so3 + hv    -> so2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    16    clo + hv    -> cl + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    17    ocs + hv    -> co + s
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ocs)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    18    cocl2 + hv  -> cl + cl + co
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cocl2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    19    h2so4 + hv  -> so3 + h2o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2so4)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    20     o(1d) + co2 -> o + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    21    o(1d) + o2 -> o + o2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    22    clco + co2 -> cl + co + co2
!===============================
            rate_local = vphot(i_lev,22)*traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    23    s2o2 + co2 -> 2so + co2
!===============================
            rate_local = vphot(i_lev,23)*traceur(i_lev,i_s2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    24    s2 + co2 -> 2s + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    25    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    26    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    27    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    28    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    29    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    30    o2(Dg) + CO2 -> O2 + CO2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    31    o2(Dg) -> O2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
		
! DEBUT DES REACTION V3
		i_v = i_v - nb_phot_max
!===============================
!--- 32 a002: o + o + co2 -> o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 33 c008: ho2 + ho2 -> h2o2 + o2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 34 c013: oh + oh -> h2o + o
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 35 c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 36 c017: oh + oh + co2 -> h2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 37 c018: h + h + co2 -> h2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 38 f021: cl + cl + co2 -> cl2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 39 f026: clco + clco  -> cocl2 + co
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 40 f030: clso2 + clso2  -> cl2 + so2 + so2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_clso2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 41 g015: so + so + co2 -> s2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 42 g022: s + s + co2 -> s2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 

! DEBUT DES REACTION V4

		i_v = i_v - nb_reaction_3_max

!===============================
!--- 43 a001: o + o2 + co2 -> o3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 44 a003: o + o3 -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 45 b002: o(1d) + h2o  -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 46 b003: o(1d) + h2  -> oh + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 47 b005: o(1d) + o3  -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 48 b006: o(1d) + o3  -> o2 + o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 49 c001: o + ho2 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 50 c002: o + oh -> o2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 51 c003: h + o3 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 52 c004: h + ho2 -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 53 c005: h + ho2 -> h2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 54 c006: h + ho2 -> h2o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 55 c007: oh + ho2 -> h2o + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 56 c009: oh + h2o2 -> h2o + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 57 c010: oh + h2 -> h2o + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 58 c011: h + o2 + co2 -> ho2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 59 c012: o + h2o2 -> oh + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 60 c014: oh + o3 -> ho2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 61 c015: ho2 + o3 -> oh + o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 62 e001: oh + co -> co2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 63 e002: o + co + m -> co2 + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 64 f001: hcl + o(1d) -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 65 f002: hcl + o(1d) -> h + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 66 f003: hcl + o -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 67 f004: hcl + oh -> h2o + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 68 f005: clo + o -> cl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 69 f006: clo + oh -> cl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 70 f007: clo + oh -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 71 f008: cl + h2 -> hcl + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 72 f009: cl + o3 -> clo + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 73 f010: cl + ho2 -> clo + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 74 f011: cl + ho2 -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 75 f012: cl + h2o2 -> hcl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 76 f013: cl + co + co2 -> clco + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 77 f015: clco + o2 + m -> clco3 + m
!===============================
		rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 78 & 79 f016: clco3 + cl -> cl + clo + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 80 & 81 f017: clco3 + o -> cl + o2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 82 f018: clo + ho2  -> hocl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 83 f019: oh + hocl -> h2o + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_hocl)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 84 f020: o + hocl -> oh + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 85 f022: clco + o -> cl + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 86 f023: cl2 + o(1d) -> cl + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 87 f024: cl2 + h  -> hcl + cl
!==============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 88 f025: cl + clco  -> cl2 + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 89 f027: cl + so2 + co2  -> clso2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 90 f028: clso2 + o  -> so2 + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 91 f029: clso2 + h  -> so2 + hcl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 92 f031: cl + o + co2  -> clo + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 93 f032: cl2 + o -> clo + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 94 f033: clco + oh -> hocl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 95 f034: cl2 + oh -> cl + hocl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 96 f035: clco + o -> co + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 97 f036: clco + cl2 -> cocl2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 98 f037: hcl + h -> h2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 99 f038: clco + h -> hcl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 100 f039: cl + h + m -> hcl + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 101 g001: s + o2 -> so + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 102 g002: s + o3 -> so + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 103 g003: so + o2 -> so2 + o
!===============================
             rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 104 g004: so + o3 -> so2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 105 g005: so + oh -> so2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 106 g006: s + oh -> so + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 107 g007: so + o + co2 -> so2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 108 g008: so + ho2 -> so2 + oh 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 109 g009: so2 + o + co2 -> so3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 110 g010: s + o + co2 -> so + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 111 g011: so3 + h2o -> h2so4
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 112 g012: so + clo -> so2 + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 113 g013: so + so3 -> so2 + so2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 114 g014: so3 + o -> so2 + o2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 115 & 116 g017: clco3 + so -> cl + so2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 117 g018: s + co + co2 -> ocs + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 118 g019: clco + s -> ocs + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 119 g020: so2 + oh + co2 -> hso3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 120 g021: hso3 + o2 -> ho2 + so3 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_hso3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 121 g024: s2 + o -> so + s 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 122 g025: s + ocs -> s2 +  co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 123 g026: ocs + o -> so + co
!=============================== 
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 124 g027: s + so3 -> so2 +  so
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 125 g028: s + ho2 -> so +  oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 126 g029: s + clo -> so +  cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
           END DO
         END IF
         cpt_time_rate = cpt_time_rate + 1
      
      END IF
            
      IF (cpt_time .EQ. (time_tot+nbr_pdt)) THEN
               OPEN(100,file='profile_rate_day.csv')
               DO i_lev=1,n_lev
               write (100,"(128(e15.8,','))")pres(i_lev), temperature(i_lev), (rate_day(i_lev,i_react),i_react=1,126)
               END DO
               
               OPEN(101,file='profile_rate_night.csv')
               DO i_lev=1,n_lev
               write (101,"(128(e15.8,','))") pres(i_lev), temperature(i_lev), (rate_night(i_lev,i_react),i_react=1,126)
               END DO
               
               OPEN(102,file='profile_rate_fullday.csv')
               rate_day=(rate_day+rate_night)/2.
               DO i_lev=1,n_lev
               write (102,"(128(e15.8,','))") pres(i_lev), temperature(i_lev), (rate_day(i_lev,i_react),i_react=1,126)
               END DO
               
               PRINT*,"pression top",pres(n_lev)
               PRINT*,"temp top",temperature(n_lev)
               
      END IF
      
      cpt_time = cpt_time + 1
      
      END     
