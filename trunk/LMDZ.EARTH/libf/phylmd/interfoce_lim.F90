!
! $Header$
!
SUBROUTINE interfoce_lim(itime, dtime, jour, &
     knon, knindex, &
     debut,  &
     lmt_sst_p, pctsrf_new_p)
  
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  
  IMPLICIT NONE
  
  INCLUDE "indicesol.h"
  INCLUDE "netcdf.inc"

! Cette routine sert d'interface entre le modele atmospherique et un fichier
! de conditions aux limites
!
! L. Fairhead 02/2000
!
! input:
!   itime        numero du pas de temps courant
!   dtime        pas de temps de la physique (en s)
!   jour         jour a lire dans l'annee
!   nisurf       index de la surface a traiter (1 = sol continental)
!   knon         nombre de points dans le domaine a traiter
!   knindex      index des points de la surface a traiter
!   klon         taille de la grille
!   debut        logical: 1er appel a la physique (initialisation)
!
! output:
!   lmt_sst_p      SST lues dans le fichier de CL
!   pctsrf_new-p   sous-maille fractionnelle
!


! Parametres d'entree
!****************************************************************************************
  INTEGER, INTENT(IN)                       :: itime
  INTEGER, INTENT(IN)                       :: jour
  INTEGER, INTENT(IN)                       :: knon
  INTEGER, DIMENSION(klon_loc), INTENT(IN)  :: knindex
  REAL   , INTENT(IN)                       :: dtime
  LOGICAL, INTENT(IN)                       :: debut
  
! Parametres de sortie
!****************************************************************************************
  REAL, INTENT(OUT), DIMENSION(klon_loc)       :: lmt_sst_p
  REAL, INTENT(OUT), DIMENSION(klon_loc,nbsrf) :: pctsrf_new_p


! Variables locales avec l'attribut SAVE
!****************************************************************************************
! frequence de lecture des conditions limites (en pas de physique) 
  INTEGER,SAVE                              :: lmt_pas   
  !$OMP THREADPRIVATE(lmt_pas)
! pour indiquer que le jour a lire est deja lu pour une surface precedente
  LOGICAL,SAVE                              :: deja_lu   
  !$OMP THREADPRIVATE(deja_lu)
  INTEGER,SAVE                              :: jour_lu 
  !$OMP THREADPRIVATE(jour_lu)
  CHARACTER (len = 20),SAVE                 :: fich ='limit.nc'
  !$OMP THREADPRIVATE(fich)
  LOGICAL, SAVE                             :: newlmt = .TRUE.
  !$OMP THREADPRIVATE(newlmt)
  LOGICAL, SAVE                             :: check = .FALSE.
  !$OMP THREADPRIVATE(check)
  REAL, ALLOCATABLE , SAVE, DIMENSION(:)    :: sst_lu_p
  !$OMP THREADPRIVATE(sst_lu_p)
  REAL, ALLOCATABLE , SAVE, DIMENSION(:,:)  :: pct_tmp_p
  !$OMP THREADPRIVATE(pct_tmp_p)

! Variables locales 
!****************************************************************************************
  INTEGER                                   :: nid, nvarid
  INTEGER                                   :: ii
  INTEGER                                   :: ierr
  INTEGER, DIMENSION(2)                     :: start, epais
  CHARACTER (len = 20)                      :: modname = 'interfoce_lim'
  CHARACTER (len = 80)                      :: abort_message
  REAL, DIMENSION(klon_glo,nbsrf)           :: pctsrf_new
  REAL, DIMENSION(klon_glo,nbsrf)           :: pct_tmp
  REAL, DIMENSION(klon_glo)                 :: sst_lu
  REAL, DIMENSION(klon_glo)                 :: nat_lu
!
! Fin declaration
!****************************************************************************************

!****************************************************************************************
! Start calculation
!
!****************************************************************************************
  IF (debut .AND. .NOT. ALLOCATED(sst_lu_p)) THEN
     lmt_pas = NINT(86400./dtime * 1.0) ! pour une lecture une fois par jour
     jour_lu = jour - 1
     ALLOCATE(sst_lu_p(klon_loc))
     ALLOCATE(pct_tmp_p(klon_loc,nbsrf))
  ENDIF
  
  IF ((jour - jour_lu) /= 0) deja_lu = .FALSE.
  
  IF (check) WRITE(*,*) modname, ' :: jour, jour_lu, deja_lu', jour, jour_lu, deja_lu 
  IF (check) WRITE(*,*) modname, ' :: itime, lmt_pas ', itime, lmt_pas,dtime

!****************************************************************************************
! Ouverture et lecture du fichier pour le master process si c'est le bon moment
!
!****************************************************************************************
! Tester d'abord si c'est le moment de lire le fichier
  IF (MOD(itime-1, lmt_pas) == 0 .AND. .NOT. deja_lu) THEN

!$OMP MASTER
     IF (is_mpi_root) THEN

        fich = TRIM(fich)
        ierr = NF_OPEN (fich, NF_NOWRITE,nid)
        IF (ierr.NE.NF_NOERR) THEN
           abort_message = 'Pb d''ouverture du fichier de conditions aux limites'
           CALL abort_gcm(modname,abort_message,1)
        ENDIF

        ! La tranche de donnees a lire:

        start(1) = 1
        start(2) = jour
        epais(1) = klon_glo
        epais(2) = 1

        IF (newlmt) THEN
           !
           ! Fraction "ocean" 
           !
           ierr = NF_INQ_VARID(nid, 'FOCE', nvarid)
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Le champ <FOCE> est absent'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais,pct_tmp(1,is_oce))
#else
           ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_oce))
#endif
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Lecture echouee pour <FOCE>'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
           !
           ! Fraction "glace de mer" 
           !
           ierr = NF_INQ_VARID(nid, 'FSIC', nvarid)
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Le champ <FSIC> est absent'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais,pct_tmp(1,is_sic))
#else
           ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_sic))
#endif
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Lecture echouee pour <FSIC>'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
           !
           ! Fraction "terre" 
           !
           ierr = NF_INQ_VARID(nid, 'FTER', nvarid)
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Le champ <FTER> est absent'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais,pct_tmp(1,is_ter))
#else
           ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_ter))
#endif
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Lecture echouee pour <FTER>'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
           !
           ! Fraction "glacier terre" 
           !
           ierr = NF_INQ_VARID(nid, 'FLIC', nvarid)
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Le champ <FLIC> est absent'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais,pct_tmp(1,is_lic))
#else
           ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais,pct_tmp(1,is_lic))
#endif
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Lecture echouee pour <FLIC>'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
           !
        ELSE  ! on en est toujours a rnatur
           ! 
           ierr = NF_INQ_VARID(nid, 'NAT', nvarid)
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Le champ <NAT> est absent'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
#ifdef NC_DOUBLE
           ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais, nat_lu)
#else
           ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, nat_lu)
#endif
           IF (ierr /= NF_NOERR) THEN
              abort_message = 'Lecture echouee pour <NAT>'
              CALL abort_gcm(modname,abort_message,1)
           ENDIF
!
! Remplissage des fractions de surface
! nat = 0, 1, 2, 3 pour ocean, terre, glacier, seaice
! 
           pct_tmp = 0.0
           DO ii = 1, klon_glo
              pct_tmp(ii,NINT(nat_lu(ii)) + 1) = 1.
           ENDDO

!
!  On se retrouve avec ocean en 1 et terre en 2 alors qu'on veut le contraire
!
           pctsrf_new = pct_tmp
           pctsrf_new (:,2)= pct_tmp (:,1)
           pctsrf_new (:,1)= pct_tmp (:,2)
           pct_tmp = pctsrf_new 
        ENDIF ! fin test sur newlmt
!
! Lecture SST
!
        ierr = NF_INQ_VARID(nid, 'SST', nvarid)
        IF (ierr /= NF_NOERR) THEN
           abort_message = 'Le champ <SST> est absent'
           CALL abort_gcm(modname,abort_message,1)
        ENDIF
#ifdef NC_DOUBLE
        ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start,epais, sst_lu)
#else
        ierr = NF_GET_VARA_REAL(nid,nvarid,start,epais, sst_lu)
#endif
        IF (ierr /= NF_NOERR) THEN
           abort_message = 'Lecture echouee pour <SST>'
           CALL abort_gcm(modname,abort_message,1)
        ENDIF
          
!****************************************************************************************
! Fin de lecture, fermeture de fichier
!
!****************************************************************************************
        ierr = NF_CLOSE(nid)
     ENDIF ! is_mpi_root

!$OMP END MASTER
!$OMP BARRIER


!****************************************************************************************
! Distribue les variables sur tous les processus
!
!****************************************************************************************
     CALL Scatter(sst_lu,sst_lu_p)
     CALL Scatter(pct_tmp(:,is_oce),pct_tmp_p(:,is_oce))
     CALL Scatter(pct_tmp(:,is_sic),pct_tmp_p(:,is_sic))
     deja_lu = .TRUE.
     jour_lu = jour
  ENDIF

!****************************************************************************************
! Recopie des variables dans les champs de sortie
!
!****************************************************************************************
  lmt_sst_p = 999999999.
  
  DO ii = 1, knon
     lmt_sst_p(ii) = sst_lu_p(knindex(ii))
  ENDDO
  
  DO ii=1,klon_loc
     pctsrf_new_p(ii,is_oce)=pct_tmp_p(ii,is_oce)
     pctsrf_new_p(ii,is_sic)=pct_tmp_p(ii,is_sic)
  ENDDO
  
  
END SUBROUTINE interfoce_lim
