!
! $Id: inithist.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine inithist(day0,anne0,tstep,t_ops,t_wrt)

#ifdef CPP_IOIPSL
       USE IOIPSL
#endif
       USE infotrac, ONLY : nqtot, ttext
       use com_io_dyn_mod, only : histid,histvid,histuid,               &
     &                        dynhist_file,dynhistv_file,dynhistu_file
       USE comvert_mod, ONLY: presnivs
       USE comconst_mod, ONLY: pi
       USE temps_mod, ONLY: itau_dyn

      implicit none

C
C   Routine d'initialisation des ecritures des fichiers histoires LMDZ
C   au format IOIPSL
C
C   Appels succesifs des routines: histbeg
C                                  histhori
C                                  histver
C                                  histdef
C                                  histend
C
C   Entree:
C
C      infile: nom du fichier histoire a creer
C      day0,anne0: date de reference
C      tstep: duree du pas de temps en seconde
C      t_ops: frequence de l'operation pour IOIPSL
C      t_wrt: frequence d'ecriture sur le fichier
C      nq: nombre de traceurs
C
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "iniprint.h"

C   Arguments
C
      integer day0, anne0
      real tstep, t_ops, t_wrt

#ifdef CPP_IOIPSL
! This routine needs IOIPSL to work
C   Variables locales
C
      integer tau0
      real zjulian
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1)
      integer uhoriid, vhoriid, thoriid, zvertiid
      integer ii,jj
      integer zan, dayref
C
C  Initialisations
C
      pi = 4. * atan (1.)
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = anne0
      dayref = day0
      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
      tau0 = itau_dyn
      
! -------------------------------------------------------------
! Creation des 3 fichiers pour les grilles horizontales U,V,Scal
! -------------------------------------------------------------
!Grille U      
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonu(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo
       
      call histbeg(dynhistu_file, iip1, rlong(:,1), jjp1, rlat(1,:),
     .             1, iip1, 1, jjp1,
     .             tau0, zjulian, tstep, uhoriid, histuid)

! Grille V
      do jj = 1, jjm
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatv(jj) * 180. / pi
        enddo
      enddo

      call histbeg(dynhistv_file, iip1, rlong(:,1), jjm, rlat(1,:),
     .             1, iip1, 1, jjm,
     .             tau0, zjulian, tstep, vhoriid, histvid)

!Grille Scalaire
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo

      call histbeg(dynhist_file, iip1, rlong(:,1), jjp1, rlat(1,:),
     .             1, iip1, 1, jjp1,
     .             tau0, zjulian, tstep, thoriid, histid)
! -------------------------------------------------------------
C  Appel a histvert pour la grille verticale
! -------------------------------------------------------------
      call histvert(histid, 'presnivs', 'Niveaux pression','mb',
     .              llm, presnivs/100., zvertiid,'down')
      call histvert(histvid, 'presnivs', 'Niveaux pression','mb',
     .              llm, presnivs/100., zvertiid,'down')
      call histvert(histuid, 'presnivs', 'Niveaux pression','mb',
     .              llm, presnivs/100., zvertiid,'down')
C
! -------------------------------------------------------------
C  Appels a histdef pour la definition des variables a sauvegarder
! -------------------------------------------------------------
C
C  Vents U
C
      call histdef(histuid, 'u', 'vent u', 'm/s',
     .             iip1, jjp1, uhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Vents V
C
      call histdef(histvid, 'v', 'vent v', 'm/s',
     .             iip1, jjm, vhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Temperature potentielle
C
      call histdef(histid, 'teta', 'temperature potentielle', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Geopotentiel
C
      call histdef(histid, 'phi', 'geopotentiel', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Traceurs
C
!
!        DO iq=1,nqtot
!          call histdef(histid, ttext(iq),  ttext(iq), '-',
!     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
!     .             32, 'inst(X)', t_ops, t_wrt)
!        enddo
!C
C  Masse
C
      call histdef(histid, 'masse', 'masse', 'kg',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(histid, 'ps', 'pression naturelle au sol', 'Pa',
     .             iip1, jjp1, thoriid, 1, 1, 1, -99,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Geopotentiel au sol
!C
!      call histdef(histid, 'phis', 'geopotentiel au sol', '-',
!     .             iip1, jjp1, thoriid, 1, 1, 1, -99,
!     .             32, 'inst(X)', t_ops, t_wrt)
!C
C  Fin
C
      call histend(histid)
      call histend(histuid)
      call histend(histvid)
#else
! tell the user this routine should be run with ioipsl
      write(lunout,*)"inithist: Warning this routine should not be",
     &               " used without ioipsl"
#endif
! of #ifdef CPP_IOIPSL
      return
      end
