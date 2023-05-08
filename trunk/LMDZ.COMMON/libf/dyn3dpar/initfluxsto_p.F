!
! $Id: initfluxsto_p.F 1438 2010-10-08 10:19:34Z jghattas $
!
      subroutine initfluxsto_p
     .  (infile,tstep,t_ops,t_wrt,
     .                    fileid,filevid,filedid)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
       USE IOIPSL
#endif
       use parallel_lmdz
       use Write_field
       use misc_mod
       USE comvert_mod, ONLY: nivsigs
       USE comconst_mod, ONLY: pi
       USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn
       
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
      real tstep, t_ops, t_wrt
      integer fileid, filevid,filedid

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      real nivd(1)
      integer tau0
      real zjulian
      character*3 str
      character*10 ctrac
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1),rl(1,1)
      integer uhoriid, vhoriid, thoriid, zvertiid,dhoriid,dvertiid
      integer ii,jj
      integer zan, idayref
      logical ok_sync
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
      pi = 4. * atan (1.)
      str='q  '
      ctrac = 'traceur   '
      ok_sync = .true.
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = annee_ref
      idayref = day_ref
      CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
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
       
      call histbeg(trim(infile),iip1, rlong(:,1), jjn, rlat(1,jjb:jje),
     .             1, iip1, 1, jjn, tau0, zjulian, tstep, uhoriid,
     .             fileid,dynu_domain_id)
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
     
      call histbeg('fluxstokev',iip1, rlong(:,1), jjn, rlat(1,jjb:jje),
     .             1, iip1, 1, jjn,tau0, zjulian, tstep, vhoriid,
     .             filevid,dynv_domain_id)
	
      rl(1,1) = 1.	
      
      if (mpi_rank==0) then
          
        call histbeg('defstoke.nc', 1, rl, 1, rl,
     .               1, 1, 1, 1,
     .               tau0, zjulian, tstep, dhoriid, filedid)
     
      endif
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
     .             'scalar','Grille points scalaires', thoriid)
	
C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'sig_s', 'Niveaux sigma',
     . 'sigma_level',
     .              llm, nivsigs, zvertiid)
C Pour le fichier V
      call histvert(filevid, 'sig_s', 'Niveaux sigma',
     .  'sigma_level',
     .              llm, nivsigs, zvertiid)
c pour le fichier def
      if (mpi_rank==0) then
         nivd(1) = 1
         call histvert(filedid, 'sig_s', 'Niveaux sigma',
     .        'sigma_level',
     .        1, nivd, dvertiid)
      endif
C
C  Appels a histdef pour la definition des variables a sauvegarder
	
	CALL histdef(fileid, "phis", "Surface geop. height", "-",
     .                iip1,jjn,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)

         CALL histdef(fileid, "aire", "Grid area", "-",
     .                iip1,jjn,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
	
        if (mpi_rank==0) then
	
	CALL histdef(filedid, "dtvr", "tps dyn", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
        
         CALL histdef(filedid, "istdyn", "tps stock", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
         
         CALL histdef(filedid, "istphy", "tps stock phy", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)

        endif
C
C Masse 
C
      call histdef(fileid, 'masse', 'Masse', 'kg',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pbaru 
C
      call histdef(fileid, 'pbaru', 'flx de masse zonal', 'kg m/s',
     .             iip1, jjn, uhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Pbarv 
C
      if (pole_sud) jjn=jj_nb-1
      
      call histdef(filevid, 'pbarv', 'flx de masse mer', 'kg m/s',
     .             iip1, jjn, vhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  w 
C
      if (pole_sud) jjn=jj_nb
      call histdef(fileid, 'w', 'flx de masse vert', 'kg m/s',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Temperature potentielle
C
      call histdef(fileid, 'teta', 'temperature potentielle', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C

C
C Geopotentiel 
C
      call histdef(fileid, 'phi', 'geopotentiel instantane', '-',
     .             iip1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(fileid)
      call histend(filevid)
      if (mpi_rank==0) call histend(filedid)
      if (ok_sync) then
        call histsync(fileid)
        call histsync(filevid)
        if (mpi_rank==0) call histsync(filedid)
      endif
	
#else
      write(lunout,*)'initfluxsto_p: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      return
      end
