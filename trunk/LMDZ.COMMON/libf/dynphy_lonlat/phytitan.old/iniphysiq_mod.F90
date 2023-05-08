!
! $Id: iniphysiq.F90 2225 2015-03-11 14:55:23Z emillour $
!
MODULE iniphysiq_mod

CONTAINS

SUBROUTINE iniphysiq(ii,jj,nlayer, &
                     nbp, communicator, &
                     punjours, pdayref,ptimestep, &
                     rlatudyn,rlatvdyn,rlonudyn,rlonvdyn, &
                     airedyn,cudyn,cvdyn, &
                     prad,pg,pr,pcpp,iflag_phys)

  USE temps_mod, ONLY: annee_ref, day_ref, day_ini, day_end
  USE comconst_mod, ONLY: cpp
  USE cpdet_phy_mod, ONLY: init_cpdet_phy
  USE infotrac, ONLY: nqtot, tname, ttext
  USE logic_mod, ONLY: iflag_trac
  USE infotrac_phy, ONLY: init_infotrac_phy
  USE time_phylmdz_mod, ONLY: init_time
  USE inigeomphy_mod, ONLY: inigeomphy
  USE control_mod, ONLY: nday
  USE dimphy, ONLY: init_dimphy
  USE mod_phys_lmdz_para, ONLY: klon_omp ! number of columns (on local omp grid)
  IMPLICIT NONE

  ! =======================================================================
  ! Initialisation of the physical constants and some positional and
  ! geometrical arrays for the physics
  ! =======================================================================

  include "YOMCST.h"
  include "iniprint.h"

  REAL, INTENT (IN) :: prad ! radius of the planet (m)
  REAL, INTENT (IN) :: pg ! gravitational acceleration (m/s2)
  REAL, INTENT (IN) :: pr ! ! reduced gas constant R/mu
  REAL, INTENT (IN) :: pcpp ! specific heat Cp
  REAL, INTENT (IN) :: punjours ! length (in s) of a standard day
  INTEGER, INTENT (IN) :: nlayer ! number of atmospheric layers
  INTEGER, INTENT (IN) :: ii ! number of atmospheric columns along longitudes
  INTEGER, INTENT (IN) :: jj ! number of atompsheric columns along latitudes
  INTEGER, INTENT(IN) :: nbp ! number of physics columns for this MPI process
  INTEGER, INTENT(IN) :: communicator ! MPI communicator
  REAL, INTENT (IN) :: rlatudyn(jj+1) ! latitudes of the physics grid
  REAL, INTENT (IN) :: rlatvdyn(jj) ! latitude boundaries of the physics grid
  REAL, INTENT (IN) :: rlonvdyn(ii+1) ! longitudes of the physics grid
  REAL, INTENT (IN) :: rlonudyn(ii+1) ! longitude boundaries of the physics grid
  REAL, INTENT (IN) :: airedyn(ii+1,jj+1) ! area of the dynamics grid (m2)
  REAL, INTENT (IN) :: cudyn((ii+1)*(jj+1)) ! cu coeff. (u_covariant = cu * u)
  REAL, INTENT (IN) :: cvdyn((ii+1)*jj) ! cv coeff. (v_covariant = cv * v)
  INTEGER, INTENT (IN) :: pdayref ! reference day of for the simulation
  REAL, INTENT (IN) :: ptimestep !physics time step (s)
  INTEGER, INTENT (IN) :: iflag_phys ! type of physics to be called

  CHARACTER (LEN=20) :: modname = 'iniphysiq'
  CHARACTER (LEN=80) :: abort_message

  ! the common part for all planetary physics
  !------------------------------------------
  ! --> initialize physics distribution, global fields and geometry
  ! (i.e. things in phy_common or dynphy_lonlat)
  CALL inigeomphy(ii,jj,nlayer, &
               nbp, communicator, &
               rlatudyn,rlatvdyn, &
               rlonudyn,rlonvdyn, &
               airedyn,cudyn,cvdyn)

  ! the distinct part for all planetary physics  (ie. things in phytitan)
  !------------------------------------------

!$OMP PARALLEL

  ! Initialize dimphy module
  call init_dimphy(klon_omp,nlayer)

  ! Initialize some physical constants
  call suphec
  
  ! Initialize cpdet_phy module
  call init_cpdet_phy(cpp)

  ! Initialize some "temporal and calendar" related variables
  CALL init_time(annee_ref,day_ref,day_ini,day_end,nday,ptimestep)

  ! Initialize tracers in physics
  CALL init_infotrac_phy(iflag_trac,nqtot,tname,ttext)

!$OMP END PARALLEL


  ! check that physical constants set in 'suphec' are coherent
  ! with values set in the dynamics:
  IF (rday/=punjours) THEN
    WRITE (lunout, *) 'iniphysiq: length of day discrepancy!!!'
    WRITE (lunout, *) '  in the dynamics punjours=', punjours
    WRITE (lunout, *) '   but in the physics RDAY=', rday
    IF (abs(rday-punjours)>0.01*punjours) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'length of day discrepancy'
      CALL abort_gcm(modname, abort_message, 1)
    END IF
  END IF

  IF (rg/=pg) THEN
    WRITE (lunout, *) 'iniphysiq: gravity discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pg=', pg
    WRITE (lunout, *) '  but in the physics RG=', rg
    IF (abs(rg-pg)>0.01*pg) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'gravity discrepancy'
      CALL abort_gcm(modname, abort_message, 1)
    END IF
  END IF
  IF (ra/=prad) THEN
    WRITE (lunout, *) 'iniphysiq: planet radius discrepancy !!!'
    WRITE (lunout, *) '   in the dynamics prad=', prad
    WRITE (lunout, *) '  but in the physics RA=', ra
    IF (abs(ra-prad)>0.01*prad) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'planet radius discrepancy'
      CALL abort_gcm(modname, abort_message, 1)
    END IF
  END IF
  IF (rd/=pr) THEN
    WRITE (lunout, *) 'iniphysiq: reduced gas constant discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pr=', pr
    WRITE (lunout, *) '  but in the physics RD=', rd
    IF (abs(rd-pr)>0.01*pr) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'reduced gas constant discrepancy'
      CALL abort_gcm(modname, abort_message, 1)
    END IF
  END IF
  IF (rcpd/=pcpp) THEN
    WRITE (lunout, *) 'iniphysiq: specific heat discrepancy !!!'
    WRITE (lunout, *) '     in the dynamics pcpp=', pcpp
    WRITE (lunout, *) '  but in the physics RCPD=', rcpd
    IF (abs(rcpd-pcpp)>0.01*pcpp) THEN
        ! stop here if the relative difference is more than 1%
      abort_message = 'specific heat discrepancy'
      CALL abort_gcm(modname, abort_message, 1)
    END IF
  END IF

END SUBROUTINE iniphysiq

END MODULE iniphysiq_mod
