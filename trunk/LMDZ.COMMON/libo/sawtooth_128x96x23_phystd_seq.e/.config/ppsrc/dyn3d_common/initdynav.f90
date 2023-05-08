










! $Id: initdynav.F90 1611 2012-01-25 14:31:54Z lguez $

subroutine initdynav(day0,anne0,tstep,t_ops,t_wrt)

  USE infotrac, ONLY : nqtot, ttext
  use com_io_dyn_mod, only : histaveid,histvaveid,histuaveid, &
       dynhistave_file,dynhistvave_file,dynhistuave_file
  USE comvert_mod, ONLY: presnivs
  USE comconst_mod, ONLY: pi
  USE temps_mod, ONLY: itau_dyn
  implicit none


  !   Routine d'initialisation des ecritures des fichiers histoires LMDZ
  !   au format IOIPSL. Initialisation du fichier histoire moyenne.

  !   Appels succesifs des routines: histbeg
  !                                  histhori
  !                                  histver
  !                                  histdef
  !                                  histend

  !   Entree:

  !      infile: nom du fichier histoire a creer
  !      day0,anne0: date de reference
  !      tstep : frequence d'ecriture
  !      t_ops: frequence de l'operation pour IOIPSL
  !      t_wrt: frequence d'ecriture sur le fichier


  !   L. Fairhead, LMD, 03/99

  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "iniprint.h"

  !   Arguments

  integer day0, anne0
  real tstep, t_ops, t_wrt

  write(lunout,*)"initdynav: Warning this routine should not be", &
       " used without ioipsl"
  ! of #ifdef CPP_IOIPSL

end subroutine initdynav
