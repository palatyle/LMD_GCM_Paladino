! $Id: writedynav.F90 1612 2012-01-31 10:11:48Z lguez $

subroutine writedynav(time, vcov, ucov, teta, ppk, phi, q, masse, ps, phis)

#ifdef CPP_IOIPSL
  USE ioipsl
#endif
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

#ifdef CPP_IOIPSL
  ! This routine needs IOIPSL to work
  !   Variables locales

  integer ndex2d(ip1jmp1), ndexu(ip1jmp1*llm), ndexv(ip1jm*llm)
  INTEGER iq, ii, ll
  real tm(ip1jmp1*llm)
  REAL vnat(ip1jm, llm), unat(ip1jmp1, llm) 
  logical ok_sync
  integer itau_w

  !-----------------------------------------------------------------

  !  Initialisations

  ndexu = 0
  ndexv = 0
  ndex2d = 0
  ok_sync = .TRUE.
  tm = 999.999
  vnat = 999.999
  unat = 999.999
  itau_w = itau_dyn + time

  ! Passage aux composantes naturelles du vent
  call covnat(llm, ucov, vcov, unat, vnat)

  !  Appels a histwrite pour l'ecriture des variables a sauvegarder

  !  Vents U 

  call histwrite(histuaveid, 'u', itau_w, unat,  &
       iip1*jjp1*llm, ndexu)

  !  Vents V

  call histwrite(histvaveid, 'v', itau_w, vnat,  &
       iip1*jjm*llm, ndexv)

  !  Temperature potentielle moyennee

  call histwrite(histaveid, 'theta', itau_w, teta,  &
       iip1*jjp1*llm, ndexu)

  !  Temperature moyennee

  do ii = 1, ijp1llm
     tm(ii) = teta(ii) * ppk(ii)/cpp
  enddo
  call histwrite(histaveid, 'temp', itau_w, tm,  &
       iip1*jjp1*llm, ndexu)

  !  Geopotentiel

  call histwrite(histaveid, 'phi', itau_w, phi,  &
       iip1*jjp1*llm, ndexu)

  !  Traceurs

  !  DO iq=1, nqtot
  !       call histwrite(histaveid, ttext(iq), itau_w, q(:, :, iq), &
  !                   iip1*jjp1*llm, ndexu)
  ! enddo

  !  Masse

  call histwrite(histaveid, 'masse', itau_w, masse,  &
       iip1*jjp1*llm, ndexu)

  !  Pression au sol

  call histwrite(histaveid, 'ps', itau_w, ps, iip1*jjp1, ndex2d)

  ! Geopotentiel au sol

  ! call histwrite(histaveid, 'phis', itau_w, phis, iip1*jjp1, ndex2d)

  if (ok_sync) then
     call histsync(histaveid)
     call histsync(histvaveid)
     call histsync(histuaveid)
  ENDIF

#else
  write(lunout, *) "writedynav: Warning this routine should not be", &
       " used without ioipsl"
#endif
  ! of #ifdef CPP_IOIPSL

end subroutine writedynav
