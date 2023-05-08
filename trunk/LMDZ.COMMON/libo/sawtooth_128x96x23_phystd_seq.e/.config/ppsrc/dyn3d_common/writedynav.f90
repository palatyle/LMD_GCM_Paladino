










! $Id: writedynav.F90 1612 2012-01-31 10:11:48Z lguez $

subroutine writedynav(time, vcov, ucov, teta, ppk, phi, q, masse, ps, phis)

  USE infotrac, ONLY : nqtot, ttext
  use com_io_dyn_mod, only : histaveid, histvaveid, histuaveid
  USE comconst_mod, ONLY: cpp
  USE temps_mod, ONLY: itau_dyn

  implicit none

  !   Ecriture du fichier histoire au format IOIPSL

  !   Appels succesifs des routines: histwrite

  !   Entree:
  !      time: temps de l'ecriture
  !      vcov: vents v covariants
  !      ucov: vents u covariants
  !      teta: temperature potentielle
  !      phi : geopotentiel instantane
  !      q   : traceurs
  !      masse: masse
  !      ps   :pression au sol
  !      phis : geopotentiel au sol

  !   L. Fairhead, LMD, 03/99

  !   Declarations
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "iniprint.h"

  !   Arguments

  REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm) 
  REAL teta(ip1jmp1*llm), phi(ip1jmp1, llm), ppk(ip1jmp1*llm)     
  REAL ps(ip1jmp1), masse(ip1jmp1, llm)                   
  REAL phis(ip1jmp1)                  
  REAL q(ip1jmp1, llm, nqtot)
  integer time

  write(lunout, *) "writedynav: Warning this routine should not be", &
       " used without ioipsl"
  ! of #ifdef CPP_IOIPSL

end subroutine writedynav
