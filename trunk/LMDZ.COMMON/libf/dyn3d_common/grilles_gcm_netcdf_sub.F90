!
! $Header$
!
! This subroutine creates the file grilles_gcm.nc containg longitudes and
! latitudes in degrees for grid u and v. This subroutine is called from
! ce0l if grilles_gcm_netcdf=TRUE. This subroutine corresponds to the first 
! part in the program create_fausse_var.
!
SUBROUTINE grilles_gcm_netcdf_sub(masque,phis)

  USE comvert_mod, ONLY: pa,preff,presnivs
  USE comconst_mod, ONLY: daysec,rad,omeg,g,kappa,cpp,pi

  IMPLICIT NONE

  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comgeom.h"
  INCLUDE "netcdf.inc"


  REAL,DIMENSION(iip1,jjp1),INTENT(IN)  :: masque ! masque terre/mer
  REAL,DIMENSION(iip1,jjp1),INTENT(IN)  :: phis   ! geopotentiel au sol

  REAL temp(iim+1,jjm+1)
  ! Attributs netcdf sortie
  INTEGER ncid_out,rcode_out
  INTEGER out_lonuid,out_lonvid,out_latuid,out_latvid,out_levid
  INTEGER out_varid
  INTEGER out_lonudim,out_lonvdim
  INTEGER out_latudim,out_latvdim,out_dim(3)
  INTEGER out_levdim

  INTEGER start(4),COUNT(4)

  INTEGER status,i,j
  REAL rlatudeg(jjp1),rlatvdeg(jjm),rlevdeg(llm)
  REAL rlonudeg(iip1),rlonvdeg(iip1)

  REAL dlon1(iip1),dlon2(iip1),dlat1(jjp1),dlat2(jjp1)
  REAL acoslat,dxkm,dykm,resol(iip1,jjp1)
  REAL,DIMENSION(iip1,jjp1)  :: phis_loc
  INTEGER masque_int(iip1,jjp1)
  INTEGER :: phis_id
  INTEGER :: area_id
  INTEGER :: mask_id
  
  rad = 6400000
  omeg = 7.272205e-05
  g = 9.8
  kappa = 0.285716
  daysec = 86400
  cpp = 1004.70885

  preff = 101325.
  pa= 50000.

  CALL conf_gcm( 99, .TRUE. )
  CALL iniconst
  CALL inigeom

  DO j=1,jjp1
     rlatudeg(j)=rlatu(j)*180./pi
  ENDDO
  DO j=1,jjm
     rlatvdeg(j)=rlatv(j)*180./pi
  ENDDO

  DO i=1,iip1
     rlonudeg(i)=rlonu(i)*180./pi + 360.
     rlonvdeg(i)=rlonv(i)*180./pi + 360.
  ENDDO


  !  2 ----- OUVERTURE DE LA SORTIE NETCDF
  ! ---------------------------------------------------
  ! CREATION OUTPUT
  ! ouverture fichier netcdf de sortie out
  status=NF_CREATE('grilles_gcm.nc',NF_NOCLOBBER,ncid_out)
  status=NF_DEF_DIM(ncid_out,'lonu',iim+1,out_lonudim)
  status=NF_DEF_DIM(ncid_out,'lonv',iim+1,out_lonvdim)
  status=NF_DEF_DIM(ncid_out,'latu',jjm+1,out_latudim)
  status=NF_DEF_DIM(ncid_out,'latv',jjm,out_latvdim)


  !   Longitudes en u
  status=NF_DEF_VAR(ncid_out,'lonu',NF_FLOAT,1,out_lonudim, out_lonuid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_lonuid,'units', 12,'degrees_east')
  status=NF_PUT_ATT_TEXT(ncid_out,out_lonuid,'long_name',9,'Longitude en u')

  !   Longitudes en v
  status=NF_DEF_VAR(ncid_out,'lonv',NF_FLOAT,1,out_lonvdim, out_lonvid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_lonvid,'units', 12,'degrees_east')
  status=NF_PUT_ATT_TEXT(ncid_out,out_lonvid,'long_name', 9,'Longitude en v')

  !   Latitude en u
  status=NF_DEF_VAR(ncid_out,'latu',NF_FLOAT,1,out_latudim, out_latuid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_latuid,'units', 13,'degrees_north')
  status=NF_PUT_ATT_TEXT(ncid_out,out_latuid,'long_name', 8,'Latitude en u')

  !  Latitude en v
  status=NF_DEF_VAR(ncid_out,'latv',NF_FLOAT,1,out_latvdim, out_latvid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_latvid,'units', 13,'degrees_north')
  status=NF_PUT_ATT_TEXT(ncid_out,out_latvid,'long_name', 8,'Latitude en v')

  !   ecriture de la grille u
  out_dim(1)=out_lonudim
  out_dim(2)=out_latudim
  status=NF_DEF_VAR(ncid_out,'grille_u',NF_FLOAT,2,out_dim, out_varid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'units', 6,'Kelvin')
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'long_name', 16,'Grille aux point u')

  !   ecriture de la grille v
  out_dim(1)=out_lonvdim
  out_dim(2)=out_latvdim
  status=NF_DEF_VAR(ncid_out,'grille_v',NF_FLOAT,2,out_dim, out_varid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'units', 6,'Kelvin')
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'long_name', 16,'Grille aux point v')

  !   ecriture de la grille u
  out_dim(1)=out_lonvdim
  out_dim(2)=out_latudim
  status=NF_DEF_VAR(ncid_out,'grille_s',NF_FLOAT,2,out_dim, out_varid)
  CALL handle_err(status)
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'units', 6,'Kelvin')
  status=NF_PUT_ATT_TEXT(ncid_out,out_varid,'long_name',16,'Grille aux point u')

  status=NF_ENDDEF(ncid_out)
  ! 5) ----- FERMETURE DES FICHIERS NETCDF------------------
  ! --------------------------------------------------------
  ! 3-b- Ecriture de la grille pour la sortie
  ! rajoute l'ecriture de la grille

#ifdef NC_DOUBLE
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_lonuid,1,iim+1,rlonudeg)
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_lonvid,1,iim+1,rlonvdeg)
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_latuid,1,jjm+1,rlatudeg)
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_latvid,1,jjm,rlatvdeg)
#else
  status=NF_PUT_VARA_REAL(ncid_out,out_lonuid,1,iim+1,rlonudeg)
  status=NF_PUT_VARA_REAL(ncid_out,out_lonvid,1,iim+1,rlonvdeg)
  status=NF_PUT_VARA_REAL(ncid_out,out_latuid,1,jjm+1,rlatudeg)
  status=NF_PUT_VARA_REAL(ncid_out,out_latvid,1,jjm,rlatvdeg)
#endif

  start(1)=1
  start(2)=1
  start(3)=1
  start(4)=1

  COUNT(1)=iim+1
  COUNT(2)=jjm+1
  COUNT(3)=1
  COUNT(4)=1

  DO j=1,jjm+1
     DO i=1,iim+1
        temp(i,j)=MOD(i,2)+MOD(j,2)
     ENDDO
  ENDDO

#ifdef NC_DOUBLE
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_varid,start, count,temp)
#else
  status=NF_PUT_VARA_REAL(ncid_out,out_varid,start, count,temp)
#endif

  ! On re-ouvre le fichier pour rajouter 4 nouvelles variables necessaire pour INCA
! lev - phis - aire - mask
  rlevdeg(:) = presnivs
  phis_loc(:,:) = phis(:,:)/g

! niveaux de pression verticaux
  status = NF_REDEF (ncid_out)
  status=NF_DEF_DIM(ncid_out,'lev',llm,out_levdim)
  
! fields
  out_dim(1)=out_lonvdim
  out_dim(2)=out_latudim

  status = nf_def_var(ncid_out,'phis',NF_FLOAT,2,out_dim,phis_id)
  CALL handle_err(status)
  status = nf_def_var(ncid_out,'aire',NF_FLOAT,2,out_dim,area_id)
  CALL handle_err(status)
  status = nf_def_var(ncid_out,'mask',NF_INT  ,2,out_dim,mask_id)
  CALL handle_err(status)

  status=NF_ENDDEF(ncid_out)

  ! ecriture des variables
#ifdef NC_DOUBLE
  status=NF_PUT_VARA_DOUBLE(ncid_out,out_levid,1,llm,rlevdeg)
#else
  status=NF_PUT_VARA_REAL(ncid_out,out_levid,1,llm,rlevdeg)
#endif

  start(1)=1
  start(2)=1
  start(3)=1
  start(4)=0
  COUNT(1)=iip1
  COUNT(2)=jjp1
  COUNT(3)=1
  COUNT(4)=0

  status = nf_put_vara_double(ncid_out, phis_id,start,count, phis_loc)
  status = nf_put_vara_double(ncid_out, area_id,start,count, aire)
  masque_int(:,:) = nINT(masque(:,:))
  status = nf_put_vara_int(ncid_out, mask_id,start,count,masque_int)
  CALL handle_err(status)
  
  ! fermeture du fichier netcdf
  CALL ncclos(ncid_out,rcode_out)

END SUBROUTINE grilles_gcm_netcdf_sub



SUBROUTINE handle_err(status)
  INCLUDE "netcdf.inc"

  INTEGER status
  IF (status.NE.nf_noerr) THEN
     PRINT *,NF_STRERROR(status)
     CALL abort_gcm('grilles_gcm_netcdf','netcdf error',1)
  ENDIF
END SUBROUTINE handle_err

