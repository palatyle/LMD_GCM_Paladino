!
! $Id: inithist_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine inithist_p(infile,day0,anne0,tstep,t_ops,t_wrt,
     .                      fileid,filevid)

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
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C      filevid:ID du fichier netcdf pour la grille v
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
      integer fileid, filevid

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      integer tau0
      real zjulian
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1)
      integer uhoriid, vhoriid, thoriid, zvertiid
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
      
      INTEGER :: dynu_domain_id
      INTEGER :: dynv_domain_id

C
C  Initialisations
C
      if (adjust) return
       
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
          rlong(ii,jj) = rlonu(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
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
     .                 'box',dynu_domain_id)
      
       call histbeg(trim(infile),iip1, rlong(:,1), jjn, 
     .              rlat(1,jjb:jje), 1, iip1, 1, jjn, tau0,
     .              zjulian, tstep, uhoriid, fileid,dynu_domain_id)
C
C  Creation du fichier histoire pour la grille en V (oblige pour l'instant,
C  IOIPSL ne permet pas de grilles avec des nombres de point differents dans 
C  un meme fichier)

      do jj = 1, jjm
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatv(jj) * 180. / pi
        enddo
      enddo

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb
      if (pole_sud) jje=jj_end-1
      if (pole_sud) jjn=jj_nb-1

      ddid=(/ 1,2 /)
      dsg=(/ iip1,jjm /)
      dsl=(/ iip1,jjn /)
      dpf=(/ 1,jjb /)
      dpl=(/ iip1,jje /)
      dhs=(/ 0,0 /)
      dhe=(/ 0,0 /)

      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',dynv_domain_id)
      
      call histbeg('dyn_histv', iip1, rlong(:,1), jjn, rlat(1,jjb:jje),
     .             1, iip1, 1, jjn, tau0, zjulian, tstep, vhoriid, 
     .             filevid,dynv_domain_id)
C
C  Appel a histhori pour rajouter les autres grilles horizontales
C
      
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb

      call histhori(fileid, iip1, rlong(:,jjb:jje),jjn,rlat(:,jjb:jje),
     .              'scalar','Grille points scalaires', thoriid)
C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'sig_s', 'Niveaux sigma','-',
     .              llm, nivsigs, zvertiid)
C Pour le fichier V
      call histvert(filevid, 'sig_s', 'Niveaux sigma','-',
     .              llm, nivsigs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder
C
C  Vents U
C
      jjn=jj_nb

      call histdef(fileid, 'ucov', 'vents u covariants', 'm/s',
     .             iip1, jjn, uhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Vents V
C
      if (pole_sud) jjn=jj_nb-1
      
      call histdef(filevid, 'vcov', 'vents v covariants', 'm/s',
     .             iip1, jjn, vhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Temperature potentielle
C
      jjn=jj_nb
      
      call histdef(fileid, 'teta', 'temperature potentielle', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Geopotentiel
C
      call histdef(fileid, 'phi', 'geopotentiel instantane', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Traceurs
C
        DO iq=1,nqtot
          call histdef(fileid, ttext(iq),  ttext(iq), '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
        enddo
C
C  Masse
C
      call histdef(fileid, 'masse', 'masse', 'kg',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(fileid, 'ps', 'pression naturelle au sol', 'Pa',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pression au sol
C
      call histdef(fileid, 'phis', 'geopotentiel au sol', '-',
     .             iip1, jjn, thoriid, 1, 1, 1, -99,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(fileid)
      call histend(filevid)
#else
      write(lunout,*)'inithist_p: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      return
      end
