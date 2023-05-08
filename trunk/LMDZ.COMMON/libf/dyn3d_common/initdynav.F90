! $Id: initdynav.F90 1611 2012-01-25 14:31:54Z lguez $

subroutine initdynav(day0,anne0,tstep,t_ops,t_wrt)

#ifdef CPP_IOIPSL
  USE IOIPSL
#endif
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

#ifdef CPP_IOIPSL
  ! This routine needs IOIPSL to work
  !   Variables locales

  integer tau0
  real zjulian
  integer iq
  real rlong(iip1,jjp1), rlat(iip1,jjp1)
  integer uhoriid, vhoriid, thoriid, zvertiid
  integer ii,jj
  integer zan, dayref

  !--------------------------------------------------------------------

  !  Initialisations

  pi = 4. * atan (1.)

  !  Appel a histbeg: creation du fichier netcdf et initialisations diverses


  zan = anne0
  dayref = day0
  CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
  tau0 = itau_dyn

  do jj = 1, jjp1
     do ii = 1, iip1
        rlong(ii,jj) = rlonv(ii) * 180. / pi
        rlat(ii,jj)  = rlatu(jj) * 180. / pi
     enddo
  enddo

  ! Creation de 3 fichiers pour les differentes grilles horizontales
  ! Restriction de IOIPSL: seulement 2 coordonnees dans le meme fichier
  ! Grille Scalaire       
  call histbeg(dynhistave_file, iip1, rlong(:,1), jjp1, rlat(1,:), &
       1, iip1, 1, jjp1, &
       tau0, zjulian, tstep, thoriid,histaveid)

  ! Creation du fichier histoire pour les grilles en V et U (oblige
  ! pour l'instant, IOIPSL ne permet pas de grilles avec des nombres
  ! de point differents dans  un meme fichier)
  ! Grille V
  do jj = 1, jjm
     do ii = 1, iip1
        rlong(ii,jj) = rlonv(ii) * 180. / pi
        rlat(ii,jj) = rlatv(jj) * 180. / pi
     enddo
  enddo

  call histbeg(dynhistvave_file, iip1, rlong(:,1), jjm, rlat(1,:), &
       1, iip1, 1, jjm, &
       tau0, zjulian, tstep, vhoriid,histvaveid)
  ! Grille U
  do jj = 1, jjp1
     do ii = 1, iip1
        rlong(ii,jj) = rlonu(ii) * 180. / pi
        rlat(ii,jj) = rlatu(jj) * 180. / pi
     enddo
  enddo

  call histbeg(dynhistuave_file, iip1, rlong(:,1),jjp1, rlat(1,:), &
       1, iip1, 1, jjp1, &
       tau0, zjulian, tstep, uhoriid,histuaveid)

  !  Appel a histvert pour la grille verticale

  call histvert(histaveid,'presnivs','Niveaux Pression approximatifs','mb', &
       llm, presnivs/100., zvertiid,'down')
  call histvert(histuaveid,'presnivs','Niveaux Pression approximatifs','mb', &
       llm, presnivs/100., zvertiid,'down')
  call histvert(histvaveid,'presnivs','Niveaux Pression approximatifs','mb', &
       llm, presnivs/100., zvertiid,'down')

  !  Appels a histdef pour la definition des variables a sauvegarder

  !  Vents U

  call histdef(histuaveid, 'u', 'vent u moyen ', &
       'm/s', iip1, jjp1, uhoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Vents V

  call histdef(histvaveid, 'v', 'vent v moyen', &
       'm/s', iip1, jjm, vhoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)


  !  Temperature

  call histdef(histaveid, 'temp', 'temperature moyenne', 'K', &
       iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Temperature potentielle

  call histdef(histaveid, 'theta', 'temperature potentielle', 'K', &
       iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Geopotentiel

  call histdef(histaveid, 'phi', 'geopotentiel moyen', '-', &
       iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Traceurs

  !        DO iq=1,nqtot
  !          call histdef(histaveid, ttext(iq), ttext(iq), '-', &
  !                  iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
  !                  32, 'ave(X)', t_ops, t_wrt)
  !        enddo

  !  Masse

  call histdef(histaveid, 'masse', 'masse', 'kg', &
       iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Pression au sol

  call histdef(histaveid, 'ps', 'pression naturelle au sol', 'Pa', &
       iip1, jjp1, thoriid, 1, 1, 1, -99, &
       32, 'ave(X)', t_ops, t_wrt)

  !  Geopotentiel au sol

  !      call histdef(histaveid, 'phis', 'geopotentiel au sol', '-', &
  !                  iip1, jjp1, thoriid, 1, 1, 1, -99, &
  !                  32, 'ave(X)', t_ops, t_wrt)

  call histend(histaveid)
  call histend(histuaveid)
  call histend(histvaveid)
#else
  write(lunout,*)"initdynav: Warning this routine should not be", &
       " used without ioipsl"
#endif
  ! of #ifdef CPP_IOIPSL

end subroutine initdynav
