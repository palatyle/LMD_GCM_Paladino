!
! $Id: initdynav_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine initdynav_p(infile,day0,anne0,tstep,t_ops,t_wrt,fileid)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
       USE IOIPSL
#endif
       use parallel_lmdz
       use Write_field
       use misc_mod
       USE infotrac
       USE comvert_mod, ONLY: nivsigs
       USE comconst_mod, ONLY: pi
       USE temps_mod, ONLY: itau_dyn

      implicit none

C
C   Routine d'initialisation des ecritures des fichiers histoires LMDZ
C   au format IOIPSL. Initialisation du fichier histoire moyenne.
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
C      tstep : frequence d'ecriture
C      t_ops: frequence de l'operation pour IOIPSL
C      t_wrt: frequence d'ecriture sur le fichier
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
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
      character*(*) infile
      integer*4 day0, anne0
      real tstep, t_ops, t_wrt
      integer fileid

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      integer thoriid, zvertiid
      integer tau0
      real zjulian
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1)
      integer ii,jj
      integer zan, dayref
      integer :: jjb,jje,jjn

! definition du domaine d'ecriture pour le rebuild

      INTEGER,DIMENSION(2) :: ddid
      INTEGER,DIMENSION(2) :: dsg
      INTEGER,DIMENSION(2) :: dsl
      INTEGER,DIMENSION(2) :: dpf
      INTEGER,DIMENSION(2) :: dpl
      INTEGER,DIMENSION(2) :: dhs
      INTEGER,DIMENSION(2) :: dhe 
      
      INTEGER :: dynave_domain_id
      
      if (adjust) return
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
      
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj)  = rlatu(jj) * 180. / pi
        enddo
      enddo

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb

      ddid=(/ 1,2 /)
      dsg=(/ iip1,jjp1 /)
      dsl=(/ iip1,jjn /)
      dpf=(/ 1,jjb /)
      dpl=(/ iip1,jje /)
      dhs=(/ 0,0 /)
      dhe=(/ 0,0 /)

      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',dynave_domain_id)
             
      call histbeg(trim(infile),iip1, rlong(:,1), jjn, rlat(1,jjb:jje),
     .             1, iip1, 1, jjn,tau0, zjulian, tstep, thoriid,
     .             fileid,dynave_domain_id)

C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'sigss', 'Niveaux sigma','Pa',
     .              llm, nivsigs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder
C
C  Vents U
C
      write(6,*)'inithistave',tstep
      call histdef(fileid, 'u', 'vents u scalaires moyennes',
     .             'm/s', iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)

C
C  Vents V
C
      call histdef(fileid, 'v', 'vents v scalaires moyennes',
     .             'm/s', iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)

C
C  Temperature
C
      call histdef(fileid, 'temp', 'temperature moyennee', 'K',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Temperature potentielle
C
      call histdef(fileid, 'theta', 'temperature potentielle', 'K',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)


C
C  Geopotentiel
C
      call histdef(fileid, 'phi', 'geopotentiel moyenne', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Traceurs
C
        DO iq=1,nqtot
          call histdef(fileid, ttext(iq), ttext(iq), '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', t_ops, t_wrt)
        enddo
C
C  Masse
C
      call histdef(fileid, 'masse', 'masse', 'kg',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(fileid, 'ps', 'pression naturelle au sol', 'Pa',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(fileid, 'phis', 'geopotentiel au sol', '-',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'ave(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(fileid)
#else
      write(lunout,*)'initdynav_p: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      return
      end
