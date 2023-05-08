!
! $Id: initphysto.F90 1454 2010-11-18 12:01:24Z fairhead $
!
SUBROUTINE initphysto(infile,tstep,t_ops,t_wrt,fileid)
  
  USE dimphy
  USE mod_phys_lmdz_para
  USE IOIPSL
  USE iophy
  USE control_mod
  
  IMPLICIT NONE

!
!   Routine d'initialisation des ecritures des fichiers histoires LMDZ
!   au format IOIPSL
!
!   Appels succesifs des routines: histbeg
!                                  histhori
!                                  histver
!                                  histdef
!                                  histend
!
!   Entree:
!
!      infile: nom du fichier histoire a creer
!      day0,anne0: date de reference
!      tstep: duree du pas de temps en seconde
!      t_ops: frequence de l'operation pour IOIPSL
!      t_wrt: frequence d'ecriture sur le fichier
!
!   Sortie:
!      fileid: ID du fichier netcdf cree
!
!   L. Fairhead, LMD, 03/99
!
! =====================================================================
!
!   Declarations
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comconst.h"
  INCLUDE "comgeom.h"
  INCLUDE "temps.h"
  INCLUDE "logic.h"
  INCLUDE "description.h"
  INCLUDE "serre.h"
  INCLUDE "indicesol.h"

!   Arguments
  CHARACTER(len=*), INTENT(IN) :: infile
  REAL, INTENT(IN)             :: tstep
  REAL, INTENT(IN)             :: t_ops
  REAL, INTENT(IN)             :: t_wrt
  INTEGER, INTENT(OUT)         :: fileid

! Variables locales
  INTEGER nhoriid, i
  INTEGER l,k
  REAL nivsigs(llm)
  INTEGER tau0
  REAL zjulian
  INTEGER iq
  INTEGER uhoriid, vhoriid, thoriid, zvertiid
  INTEGER ii,jj
  INTEGER zan, idayref
  LOGICAL ok_sync
  REAL zx_lon(iim,jjm+1), zx_lat(iim,jjm+1)
  CHARACTER(len=12) :: nvar

!  Initialisations
!
  pi = 4. * ATAN (1.)
  ok_sync= .TRUE.
!
!  Appel a histbeg: creation du fichier netcdf et initialisations diverses
!         

  zan = annee_ref
  idayref = day_ref
  CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
  tau0 = 0
  
  CALL histbeg_phy(infile,tau0, zjulian, tstep, &
       nhoriid, fileid)

!$OMP MASTER	
!  Appel a histvert pour la grille verticale
!
  DO l=1,llm
     nivsigs(l)=REAL(l)
  ENDDO
  
  CALL histvert(fileid, 'sig_s', 'Niveaux sigma', &
       'sigma_level', &
       llm, nivsigs, zvertiid)
!
!  Appels a histdef pour la definition des variables a sauvegarder
!
  CALL histdef(fileid, "phis", "Surface geop. height", "-", &
       iim,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)
  
  CALL histdef(fileid, "aire", "Grid area", "-", &
       iim,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)

  CALL histdef(fileid, "longitudes", "longitudes", "-", &
       iim,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)

  CALL histdef(fileid, "latitudes", "latitudes", "-", &
       iim,jj_nb,nhoriid, 1,1,1, -99, 32, &
       "once", t_ops, t_wrt)
! T 
  CALL histdef(fileid, 't', 'Temperature', 'K', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! mfu 
  CALL histdef(fileid, 'mfu', 'flx m. pan. mt', 'kg m/s',iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! mfd 
  CALL histdef(fileid, 'mfd', 'flx m. pan. des', 'kg m/s',iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! en_u 
  CALL histdef(fileid, 'en_u', 'flx ent pan mt', 'kg m/s', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! de_u 
  CALL histdef(fileid, 'de_u', 'flx det pan mt', 'kg m/s',iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! en_d 
  CALL histdef(fileid, 'en_d', 'flx ent pan dt', 'kg m/s', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! de_d 
  CALL histdef(fileid, 'de_d', 'flx det pan dt', 'kg m/s', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! coefh
  CALL histdef(fileid, "coefh", " ", " ", iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, "inst(X)", t_ops, t_wrt)
! fm_th
  CALL histdef(fileid, "fm_th", " ", " ",iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, "inst(X)", t_ops, t_wrt)
! en_th
  CALL histdef(fileid, "en_th", " ", " ",iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, "inst(X)", t_ops, t_wrt)
! frac_impa
  CALL histdef(fileid, 'frac_impa', ' ', ' ',iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! frac_nucl
  CALL histdef(fileid, 'frac_nucl', ' ', ' ',iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid,32, 'inst(X)', t_ops, t_wrt)
! pyu1
  CALL histdef(fileid, "pyu1", " ", " ", iim,jj_nb,nhoriid, &
       1,1,1, -99, 32, "inst(X)", t_ops, t_wrt)
! pyv1
  CALL histdef(fileid, "pyv1", " ", " ", iim,jj_nb,nhoriid, &
       1,1,1, -99, 32,"inst(X)", t_ops, t_wrt)    
! ftsol1
  CALL histdef(fileid, "ftsol1", " ", " ",iim, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol2
  CALL histdef(fileid, "ftsol2", " ", " ",iim, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol3
  CALL histdef(fileid, "ftsol3", " ", " ", iim, jj_nb, nhoriid, &
       1, 1,1, -99,32, "inst(X)", t_ops, t_wrt)
! ftsol4
  CALL histdef(fileid, "ftsol4", " ", " ",iim, jj_nb, nhoriid, &
       1, 1,1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf1
  CALL histdef(fileid, "psrf1", " ", " ",iim, jj_nb, nhoriid, &
       1, 1, 1, -99,32, "inst(X)", t_ops, t_wrt)
! psrf2
  CALL histdef(fileid, "psrf2", " ", " ",iim, jj_nb, nhoriid, &
       1, 1, 1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf3
  CALL histdef(fileid, "psrf3", " ", " ",iim, jj_nb, nhoriid, &
       1, 1, 1, -99, 32, "inst(X)", t_ops, t_wrt)
! psrf4
  CALL histdef(fileid, "psrf4", " ", " ", iim, jj_nb, nhoriid, &
       1, 1, 1, -99,32, "inst(X)", t_ops, t_wrt)
! sh
  CALL histdef(fileid, 'sh', '', '', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! da
  CALL histdef(fileid, 'da', '', '', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! mp
  CALL histdef(fileid, 'mp', '', '', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! upwd
  CALL histdef(fileid, 'upwd', '', '', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
! dnwd
  CALL histdef(fileid, 'dnwd', '', '', iim, jj_nb, nhoriid, &
       llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)

! phi
  DO k=1,llm
     IF (k<10) THEN
        WRITE(nvar,'(i1)') k
     ELSE IF (k<100) THEN
        WRITE(nvar,'(i2)') k
     ELSE
        WRITE(nvar,'(i3)') k
     END IF
     nvar='phi_lev'//trim(nvar)
     
     CALL histdef(fileid, nvar, '', '', iim, jj_nb, nhoriid, &
          llm, 1, llm, zvertiid, 32, 'inst(X)', t_ops, t_wrt)
  END DO

  CALL histend(fileid)
  IF (ok_sync) CALL histsync
!$OMP END MASTER
	
END SUBROUTINE initphysto
