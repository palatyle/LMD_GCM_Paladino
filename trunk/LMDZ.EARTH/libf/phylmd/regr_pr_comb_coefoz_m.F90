! $Id$
module regr_pr_comb_coefoz_m

  implicit none

  ! The five module variables declared here are on the partial
  ! "physics" grid.
  ! The value of each variable for index "(i, k)" is at longitude
  ! "rlon(i)", latitude "rlat(i)" and middle of layer "k".

  real, allocatable, save:: c_Mob(:, :)
  ! (sum of Mobidic terms in the net mass production rate of ozone
  ! by chemistry, per unit mass of air, in s-1)

  real, allocatable, save:: a2(:, :)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to ozone mass fraction, in s-1)

  real, allocatable, save:: a4_mass(:, :)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to temperature, in s-1 K-1)

  real, allocatable, save:: a6_mass(:, :)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to mass column-density of ozone above, in m2 s-1 kg-1)

  real, allocatable, save:: r_het_interm(:, :)
  ! (net mass production rate by heterogeneous chemistry, per unit
  ! mass of ozone, corrected for chlorine content and latitude, but
  ! not for temperature and sun direction, in s-1)

  !$omp threadprivate(c_Mob, a2, a4_mass, a6_mass, r_het_interm)

contains

  subroutine alloc_coefoz

    ! This procedure is called once per run.
    ! It allocates module variables.

    use dimphy, only: klon

    ! Variables local to the procedure:
    include "dimensions.h"

    !---------------------------------------

    !$omp master
    print *, "Call sequence information: alloc_coefoz"
    !$omp end master
    allocate(c_Mob(klon, llm), a2(klon, llm), a4_mass(klon, llm))
    allocate(a6_mass(klon, llm), r_het_interm(klon, llm))

  end subroutine alloc_coefoz

  !*******************************************************

  subroutine regr_pr_comb_coefoz(julien, rlat, paprs, pplay)

    ! "regr_pr_comb_coefoz" stands for "regrid pressure combine
    ! coefficients ozone".

    ! In this subroutine:
    ! -- the master thread of the root process reads from a file all
    !    eight coefficients for ozone chemistry, at the current day;
    ! -- the coefficients are packed to the "physics" horizontal grid
    !    and scattered to all threads of all processes;
    ! -- in all the threads of all the processes, the coefficients are
    !    regridded in pressure to the LMDZ vertical grid;
    ! -- in all the threads of all the processes, the eight
    !    coefficients are combined to define the five module variables.

    use netcdf95, only: nf95_open, nf95_close
    use netcdf, only: nf90_nowrite
    use assert_m, only: assert
    use dimphy, only: klon
    use mod_phys_lmdz_mpi_data, only: is_mpi_root
    use regr_pr_av_m, only: regr_pr_av
    use regr_pr_int_m, only: regr_pr_int
    use press_coefoz_m, only: press_in_edg, plev

    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    REAL, intent(in):: rlat(:)
    ! (latitude on the partial "physics" grid, in degrees)

    real, intent(in):: paprs(:, :) ! (klon, llm + 1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(in):: pplay(:, :) ! (klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    ! Variables local to the procedure:

    include "dimensions.h"
    integer ncid ! for NetCDF

    real coefoz(klon, llm, 7)
    ! (temporary storage for 7 ozone coefficients)
    ! (On the partial "physics" grid.
    ! "coefoz(i, k, :)" is at longitude "rlon(i)", latitude "rlat(i)",
    ! middle of layer "k".)

    real a6(klon, llm)
    ! (derivative of "P_net_Mob" with respect to column-density of ozone
    ! above, in cm2 s-1)
    ! (On the partial "physics" grid.
    ! "a6(i, k)" is at longitude "rlon(i)", latitude "rlat(i)",
    ! middle of layer "k".)

    real, parameter:: amu = 1.6605402e-27 ! atomic mass unit, in kg

    real, parameter:: Clx = 3.8e-9
    ! (total chlorine content in the upper stratosphere)

    integer k

    !------------------------------------

    !!print *, "Call sequence information: regr_pr_comb_coefoz"
    call assert((/size(rlat), size(paprs, 1), size(pplay, 1)/) == klon, &
         "regr_pr_comb_coefoz klon")
    call assert((/size(paprs, 2) - 1, size(pplay, 2)/) == llm, &
         "regr_pr_comb_coefoz llm")

    !$omp master
    if (is_mpi_root) call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)
    !$omp end master

    call regr_pr_av(ncid, (/"a2       ", "a4       ", "a6       ", &
         "P_net_Mob", "r_Mob    ", "temp_Mob ", "R_Het    "/), julien, &
         press_in_edg, paprs, coefoz)
    a2 = coefoz(:, :, 1)
    a4_mass = coefoz(:, :, 2) * 48. / 29.

    ! Compute "a6_mass" avoiding underflow, do not divide by 1e4
    ! before dividing by molecular mass:
    a6_mass = coefoz(:, :, 3) / (1e4 * 29. * amu)
    ! (factor 1e4: conversion from cm2 to m2)

    ! We can overwrite "coefoz(:, :, 1)", which was saved to "a2":
    call regr_pr_int(ncid, "Sigma_Mob", julien, plev, pplay, top_value=0., &
         v3=coefoz(:, :, 1))

    ! Combine coefficients to get "c_Mob":
    c_mob = (coefoz(:, :, 4) - a2 * coefoz(:, :, 5) &
         - coefoz(:, :, 3) * coefoz(:, :, 1)) * 48. / 29. &
         - a4_mass * coefoz(:, :, 6)

    r_het_interm = coefoz(:, :, 7)
    ! Heterogeneous chemistry is only at high latitudes:
    forall (k = 1: llm)
       where (abs(rlat) <= 45.) r_het_interm(:, k) = 0.
    end forall
    r_het_interm = r_het_interm * (Clx / 3.8e-9)**2

    !$omp master
    if (is_mpi_root) call nf95_close(ncid)
    !$omp end master

  end subroutine regr_pr_comb_coefoz

end module regr_pr_comb_coefoz_m
