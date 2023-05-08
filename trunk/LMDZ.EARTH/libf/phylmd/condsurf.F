c $Header$
c
      SUBROUTINE condsurf( jour, jourvrai, lmt_bils )
      USE dimphy
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      IMPLICIT none
c
c I. Musat 05.2005
c
c Lire chaque jour le bilan de chaleur au sol issu 
c d'un run atmospherique afin de l'utiliser dans
c dans un run "slab" ocean 
c -----------------------------------------
c jour     : input  , numero du jour a lire
c jourvrai : input  , vrai jour de la simulation  
c
c lmt_bils: bilan chaleur au sol (a utiliser pour "slab-ocean")
c
#include "netcdf.inc"
      INTEGER nid, nvarid
      INTEGER debut(2)
      INTEGER epais(2)
c
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "indicesol.h"
#include "temps.h"
#include "clesphys.h"
c
      INTEGER     nannemax
      PARAMETER ( nannemax = 60 )
c
      INTEGER jour, jourvrai
      REAL lmt_bils(klon) !bilan chaleur au sol
c
c Variables locales:
      INTEGER ig, i, kt, ierr
      LOGICAL ok
      INTEGER anneelim,anneemax
      CHARACTER*20 fich
      
      REAL :: lmt_bils_glo(klon_glo)
      
cc
cc   .....................................................................
cc
cc    Pour lire le fichier limit correspondant vraiment  a l'annee de la
cc     simulation en cours , il suffit de mettre  ok_limitvrai = .TRUE.
cc
cc   ......................................................................
c
c
      
      IF (jour.LT.0 .OR. jour.GT.(360-1)) THEN
         PRINT*,'Le jour demande n est pas correct: ', jour
         CALL ABORT
      ENDIF
c
       anneelim  = annee_ref
       anneemax  = annee_ref + nannemax
c
c
       IF( ok_limitvrai )       THEN
          DO  kt = 1, nannemax
           IF(jourvrai.LE. (kt-1)*360 + 359  )  THEN
              WRITE(fich,'("limit",i4,".nc")') anneelim
c             PRINT *,' Fichier  Limite ',fich
              GO TO 100
             ENDIF
           anneelim = anneelim + 1
          ENDDO

         PRINT *,' PBS ! Le jour a lire sur le fichier limit ne se '
         PRINT *,' trouve pas sur les ',nannemax,' annees a partir de '
         PRINT *,' l annee de debut', annee_ref
         CALL EXIT(1)
c
100     CONTINUE
c
       ELSE
     
            WRITE(fich,'("limitNEW.nc")') 
c           PRINT *,' Fichier  Limite ',fich
       ENDIF
c
c Ouvrir le fichier en format NetCDF:
c
c$OMP MASTER
      IF (is_mpi_root) THEN
      ierr = NF_OPEN (fich, NF_NOWRITE,nid)
      IF (ierr.NE.NF_NOERR) THEN
        WRITE(6,*)' Pb d''ouverture du fichier ', fich
        WRITE(6,*)' Le fichier limit ',fich,' (avec 4 chiffres , pour' 
        WRITE(6,*)'       l an 2000 )  ,  n existe  pas !  ' 
        WRITE(6,*)' ierr = ', ierr
        CALL EXIT(1)
      ENDIF
c     DO k = 1, jour
c La tranche de donnees a lire:
c
      debut(1) = 1
      debut(2) = jourvrai
      epais(1) = klon_glo
      epais(2) = 1
c Bilan flux de chaleur au sol:
c
      ierr = NF_INQ_VARID (nid, "BILS", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "condsurf: Le champ <BILS> est absent"
         CALL abort
      ENDIF
      PRINT*,'debut,epais',debut,epais,'jour,jourvrai',jour,jourvrai
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,debut,epais,lmt_bils_glo)
#else
      ierr = NF_GET_VARA_REAL(nid, nvarid,debut,epais,lmt_bils_glo)
#endif
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "condsurf: Lecture echouee pour <BILS>"
         CALL abort
      ENDIF
c     ENDDO !k = 1, jour
c
c Fermer le fichier:
c
      ierr = NF_CLOSE(nid)
      
      ENDIF ! is_mpi_root==0

c$OMP END MASTER
      CALL scatter(lmt_bils_glo,lmt_bils)
            
c
c
c     PRINT*, 'lmt_bils est lu pour jour: ', jour
c
      RETURN
      END
