!
! $Id: startvar.F90 1425 2010-09-02 13:45:23Z lguez $
!
!*******************************************************************************
!
MODULE startvar
!
!*******************************************************************************
!
!-------------------------------------------------------------------------------
! Purpose: Access data from the database of atmospheric to initialize the model.
!-------------------------------------------------------------------------------
! Comments:
!
!    *  This module is designed to work for Earth (and with ioipsl)
!
!    *  There are three ways to acces data, depending on the type of field
!  which needs to be extracted. In any case the call should come after a restget
!  and should be of the type :                     CALL startget(...)
!
!  - A 1D variable on the physical grid :
!    CALL startget_phys1d((varname, iml, jml,  lon_in,  lat_in,  nbindex,              &
!           champ, val_exp,      jml2, lon_in2, lat_in2, ibar )
!
!  - A 2D variable on the dynamical grid :
!    CALL startget_phys2d(varname, iml, jml,  lon_in,  lat_in,                        &
!           champ, val_exp,      jml2, lon_in2, lat_in2, ibar )             
!
!  - A 3D variable on the dynamical grid :
!    CALL startget_dyn((varname, iml, jml,  lon_in,  lat_in,  lml, pls, workvar,    &
!           champ, val_exp,      jml2, lon_in2, lat_in2, ibar )
!
!    *  Data needs to be in NetCDF format
!
!    *  Variables should have the following names in the files:
!            'RELIEF' : High resolution orography 
!            'ST'     : Surface temperature
!            'CDSW'   : Soil moisture
!            'Z'      : Surface geopotential
!            'SP'     : Surface pressure
!            'U'      : East ward wind
!            'V'      : Northward wind
!            'TEMP'   : Temperature
!            'R'      : Relative humidity
!
!   *   There is a big mess with the longitude size. Should it be iml or iml+1 ?
!  I have chosen to use the iml+1 as an argument to this routine and we declare
!  internaly smaller fields when needed. This needs to be cleared once and for
!  all in LMDZ. A convention is required.
!-------------------------------------------------------------------------------
  USE ioipsl
  IMPLICIT NONE

  PRIVATE
  PUBLIC startget_phys2d, startget_phys1d, startget_dyn
!  INTERFACE startget
!    MODULE PROCEDURE startget_phys1d, startget_phys2d, startget_dyn
!  END INTERFACE

  REAL,    SAVE :: deg2rad,  pi
  INTEGER, SAVE ::           iml_rel,  jml_rel
  INTEGER, SAVE :: fid_phys, iml_phys, jml_phys
  INTEGER, SAVE :: fid_dyn,  iml_dyn,  jml_dyn,  llm_dyn,  ttm_dyn
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: lon_phys, lon_dyn
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: lat_phys, lat_dyn
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: lon_rug, lon_alb, lon_rel
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: lat_rug, lat_alb, lat_rel
  REAL, DIMENSION(:),     ALLOCATABLE, TARGET, SAVE :: levdyn_ini
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: relief, zstd, zsig, zgam
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: zthe, zpic, zval
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: rugo, phis, tsol, qsol
  REAL, DIMENSION(:,:),   ALLOCATABLE, TARGET, SAVE :: psol_dyn
  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET, SAVE :: var_ana3d

   CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE startget_phys1d(varname, iml, jml, lon_in, lat_in, nbindex, champ,  &
                           val_exp ,jml2, lon_in2, lat_in2, ibar)
!
!-------------------------------------------------------------------------------
! Comment:
!   This routine only works if the variable does not exist or is constant.
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*),         INTENT(IN)    :: varname
  INTEGER,                  INTENT(IN)    :: iml, jml
  REAL, DIMENSION(iml),     INTENT(IN)    :: lon_in
  REAL, DIMENSION(jml),     INTENT(IN)    :: lat_in
  INTEGER,                  INTENT(IN)    :: nbindex
  REAL, DIMENSION(nbindex), INTENT(INOUT) :: champ
  REAL,                     INTENT(IN)    :: val_exp
  INTEGER,                  INTENT(IN)    :: jml2
  REAL, DIMENSION(iml),     INTENT(IN)    :: lon_in2
  REAL, DIMENSION(jml2),    INTENT(IN)    :: lat_in2
  LOGICAL,                  INTENT(IN)    :: ibar
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  REAL, DIMENSION(:,:), POINTER :: v2d
!-------------------------------------------------------------------------------
  v2d=>NULL()
  IF(MINVAL(champ)==MAXVAL(champ).AND.MINVAL(champ)==val_exp) THEN

!--- CHECKING IF THE FIELD IS KNOWN ; READING UNALLOCATED FILES
    SELECT CASE(varname)
      CASE('tsol')
        IF(.NOT.ALLOCATED(tsol))                                               &
         CALL start_init_phys(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
      CASE('qsol')
        IF(.NOT.ALLOCATED(qsol))                                               &
         CALL start_init_phys(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
      CASE('psol')
        IF(.NOT.ALLOCATED(psol_dyn))                                           &
         CALL start_init_dyn (iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
      CASE('zmea','zstd','zsig','zgam','zthe','zpic','zval')
        IF(.NOT.ALLOCATED(relief))                                             &
         CALL start_init_orog(iml,jml,lon_in,lat_in)
      CASE('rads','snow','tslab','seaice','rugmer','agsno')
      CASE DEFAULT
        WRITE(lunout,*)'startget_phys1d'
        WRITE(lunout,*)'No rule is present to extract variable '//TRIM(varname)&
                     //' from any data set'; STOP
    END SELECT

!--- SELECTING v2d FOR WANTED VARIABLE AND CHEKING ITS SIZE
    SELECT CASE(varname)
      CASE('rads','snow','tslab','seaice');  champ=0.0
      CASE('rugmer');                        champ(:)=0.001
      CASE('agsno');                         champ(:)=50.0
      CASE DEFAULT
        SELECT CASE(varname)
          CASE('tsol'); v2d=>tsol
          CASE('qsol'); v2d=>qsol
          CASE('psol'); v2d=>psol_dyn
          CASE('zmea'); v2d=>relief
          CASE('zstd'); v2d=>zstd
          CASE('zsig'); v2d=>zsig
          CASE('zgam'); v2d=>zgam
          CASE('zthe'); v2d=>zthe
          CASE('zpic'); v2d=>zpic
          CASE('zval'); v2d=>zval
        END SELECT
        IF(SIZE(v2d)/=SIZE(lon_in)*SIZE(lat_in)) THEN
         WRITE(lunout,*)'STARTVAR module has been initialized to the wrong size'
         STOP
        END IF
        CALL gr_dyn_fi(1,iml,jml,nbindex,v2d,champ)
    END SELECT

  ELSE

!--- SOME FIELDS ARE CAUGHT: MAY BE NEEDED FOR A 3D INTEPROLATION
    SELECT CASE(varname)
      CASE('tsol')
        IF(.NOT.ALLOCATED(tsol)) ALLOCATE(tsol(iml,jml))
        CALL gr_fi_dyn(1,iml,jml,nbindex,champ,tsol)
    END SELECT

  END IF

END SUBROUTINE  startget_phys1d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE startget_phys2d(varname, iml, jml, lon_in, lat_in, champ, val_exp,  &
                           jml2, lon_in2, lat_in2 , ibar)
!
!-------------------------------------------------------------------------------
! Comment:
!   This routine only works if the variable does not exist or is constant.
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*),         INTENT(IN)           :: varname
  INTEGER,                  INTENT(IN)           :: iml, jml
  REAL, DIMENSION(iml),     INTENT(IN)           :: lon_in
  REAL, DIMENSION(jml),     INTENT(IN)           :: lat_in
  REAL, DIMENSION(iml,jml), INTENT(INOUT)        :: champ
  REAL,                     INTENT(IN)           :: val_exp
  INTEGER,                  INTENT(IN)           :: jml2
  REAL, DIMENSION(iml),     INTENT(IN)           :: lon_in2
  REAL, DIMENSION(jml2),    INTENT(IN)           :: lat_in2
  LOGICAL,                  INTENT(IN)           :: ibar
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  REAL, DIMENSION(:,:), POINTER :: v2d=>NULL()
!-------------------------------------------------------------------------------
  v2d=>NULL()
  IF(MINVAL(champ)==MAXVAL(champ).AND.MINVAL(champ)==val_exp) THEN

!--- CHECKING IF THE FIELD IS KNOWN ; READING UNALLOCATED FILES
    SELECT CASE(varname)
      CASE('psol')
        IF(.NOT.ALLOCATED(psol_dyn))                                           &
          CALL start_init_dyn (iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
      CASE('relief')
        IF(.NOT.ALLOCATED(relief)) CALL start_init_orog(iml,jml,lon_in,lat_in)
      CASE('surfgeo')
        IF(.NOT.ALLOCATED(phis)) CALL start_init_orog(iml,jml,lon_in,lat_in)
      CASE('rugosite')
        IF(.NOT.ALLOCATED(rugo)) CALL start_init_orog(iml,jml,lon_in,lat_in)
      CASE DEFAULT
        WRITE(lunout,*)'startget_phys2d'
        WRITE(lunout,*)'No rule is present to extract variable '//TRIM(varname)&
                     //' from any data set'; STOP
    END SELECT

!--- SELECTING v2d FOR WANTED VARIABLE AND CHEKING ITS SIZE
    SELECT CASE(varname)
      CASE('psol');     v2d=>psol_dyn
      CASE('relief');   v2d=>relief
      CASE('rugosite'); v2d=>rugo
      CASE('surfgeo');  v2d=>phis
    END SELECT
    IF(SIZE(champ)/=SIZE(v2d)) THEN
      WRITE(lunout,*) 'STARTVAR module has been initialized to the wrong size'
      STOP
    END IF

    champ(:,:)=v2d(:,:)

  ELSE

!--- SOME FIELDS ARE CAUGHT: MAY BE NEEDED FOR A 3D INTEPROLATION
    SELECT CASE(varname)
      CASE ('surfgeo')
        IF(.NOT.ALLOCATED(phis)) ALLOCATE(phis(iml,jml))
        IF(SIZE(phis)/=SIZE(champ)) THEN
         WRITE(lunout,*)'STARTVAR module has been initialized to the wrong size'
         STOP
        END IF
        phis(:,:)=champ(:,:)
    END SELECT

  END IF

END SUBROUTINE startget_phys2d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE startget_dyn(varname,  lon_in,  lat_in, pls,workvar,&
                     champ, val_exp, lon_in2, lat_in2, ibar)

      use assert_eq_m, only: assert_eq
      USE comconst_mod

!-------------------------------------------------------------------------------
! Comment:
!   This routine only works if the variable does not exist or is constant.
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN)    :: varname
  REAL, INTENT(IN)    :: lon_in(:) ! dim(iml)
  REAL, INTENT(IN)    :: lat_in(:) ! dim(jml)
  REAL, INTENT(IN)    :: pls(:, :, :) ! dim(iml, jml, lml)
  REAL, INTENT(IN)    :: workvar(:, :, :) ! dim(iml, jml, lml)
  REAL, INTENT(INOUT) :: champ(:, :, :) ! dim(iml, jml, lml)
  REAL, INTENT(IN)    :: val_exp
  REAL, INTENT(IN)    :: lon_in2(:) ! dim(iml)
  REAL, INTENT(IN)    :: lat_in2(:) ! dim(jml2)
  LOGICAL,                      INTENT(IN)    :: ibar
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
  INTEGER    :: iml, jml
  INTEGER    :: lml
  INTEGER    :: jml2
  REAL, DIMENSION(:,:,:), POINTER :: v3d=>NULL()
  CHARACTER(LEN=10) :: vname
  INTEGER :: il
  REAL    :: xppn, xpps
!-------------------------------------------------------------------------------
  NULLIFY(v3d)
  IF(MINVAL(champ)==MAXVAL(champ).AND.MINVAL(champ)==val_exp) THEN

      iml = assert_eq((/size(lon_in), size(pls, 1), size(workvar, 1), &
     &     size(champ, 1), size(lon_in2)/), "startget_dyn iml")
      jml = assert_eq(size(lat_in), size(pls, 2), size(workvar, 2),   &
     &     size(champ, 2), "startget_dyn jml")
      lml = assert_eq(size(pls, 3), size(workvar, 3), size(champ, 3), &
     &     "startget_dyn lml")
      jml2 = size(lat_in2)

!--- READING UNALLOCATED FILES
    IF(.NOT.ALLOCATED(psol_dyn))                                               &
      CALL start_init_dyn(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)

!--- CHECKING IF THE FIELD IS KNOWN AND INTERPOLATING 3D FIELDS
    SELECT CASE(varname)
      CASE('u');        vname='U'
      CASE('v');        vname='V'
      CASE('t','tpot'); vname='TEMP'
      CASE('q');        vname='R'
      CASE DEFAULT
        WRITE(lunout,*)'startget_dyn'
        WRITE(lunout,*)'No rule is present to extract variable '//TRIM(varname)&
                //' from any data set'; STOP
    END SELECT
    CALL start_inter_3d(TRIM(vname), iml, jml, lml, lon_in, lat_in, jml2,      &
                        lon_in2, lat_in2,  pls, champ,ibar )

!--- COMPUTING THE REQUIRED FILED
    SELECT CASE(varname)
      CASE('u')                                            !--- Eastward wind
        DO il=1,lml; champ(:,:,il)=champ(:,:,il)*cu(:,1:jml); END DO
        champ(iml,:,:)=champ(1,:,:)

      CASE('v')                                            !--- Northward wind
        DO il=1,lml; champ(:,:,il)=champ(:,:,il)*cv(:,1:jml); END DO
        champ(iml,:,:)=champ(1,:,:)

      CASE('tpot')                                         !--- Temperature
        IF(MINVAL(workvar)/=MAXVAL(workvar)) THEN
          champ=champ*cpp/workvar
          DO il=1,lml
            xppn = SUM(aire(:,1  )*champ(:,1  ,il))/apoln
            xpps = SUM(aire(:,jml)*champ(:,jml,il))/apols
            champ(:,1  ,il) = xppn
            champ(:,jml,il) = xpps
          END DO
        ELSE
          WRITE(lunout,*)'Could not compute potential temperature as the'
          WRITE(lunout,*)'Exner function is missing or constant.'; STOP
        END IF

      CASE('q')                                            !--- Relat. humidity
        IF(MINVAL(workvar)/=MAXVAL(workvar)) THEN
          champ=0.01*champ*workvar
          WHERE(champ<0.) champ=1.0E-10
          DO il=1,lml
            xppn = SUM(aire(:,1  )*champ(:,1  ,il))/apoln
            xpps = SUM(aire(:,jml)*champ(:,jml,il))/apols
            champ(:,1  ,il) = xppn
            champ(:,jml,il) = xpps
          END DO
        ELSE
          WRITE(lunout,*)'Could not compute specific humidity as the'
          WRITE(lunout,*)'saturated humidity is missing or constant.'; STOP
        END IF

    END SELECT

  END IF

END SUBROUTINE startget_dyn
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_orog(iml,jml,lon_in,lat_in)
      USE comconst_mod
      USE grid_atob_m, ONLY: rugsoro

!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,                  INTENT(IN)           :: iml, jml
  REAL, DIMENSION(iml),     INTENT(IN)           :: lon_in
  REAL, DIMENSION(jml),     INTENT(IN)           :: lat_in
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  CHARACTER(LEN=25)     :: title
  CHARACTER(LEN=120)    :: orofname
  LOGICAL               :: check=.TRUE.
  REAL,    DIMENSION(1) :: lev
  REAL                  :: date, dt
  INTEGER, DIMENSION(1) :: itau
  INTEGER               :: fid, llm_tmp, ttm_tmp
  REAL,    DIMENSION(:,:), ALLOCATABLE :: relief_hi, tmp_var
  REAL,    DIMENSION(:),   ALLOCATABLE :: lon_rad, lat_rad, lon_ini, lat_ini
!-------------------------------------------------------------------------------
  pi=2.0*ASIN(1.0); deg2rad=pi/180.0

  orofname = 'Relief.nc'; title='RELIEF'
  IF(check) WRITE(lunout,*)'Reading the high resolution orography'
  CALL flininfo(orofname, iml_rel, jml_rel, llm_tmp, ttm_tmp, fid)

  ALLOCATE(lat_rel(iml_rel,jml_rel),lon_rel(iml_rel,jml_rel))
  CALL flinopen(orofname, .FALSE., iml_rel, jml_rel, llm_tmp, lon_rel, lat_rel,&
                lev, ttm_tmp, itau, date, dt, fid)
  ALLOCATE(relief_hi(iml_rel,jml_rel))
  CALL flinget(fid, title, iml_rel, jml_rel, llm_tmp, ttm_tmp, 1, 1, relief_hi)
  CALL flinclo(fid)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_rel),lat_ini(jml_rel))
  lon_ini(:)=lon_rel(:,1); IF(MAXVAL(lon_rel)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_rel(1,:); IF(MAXVAL(lat_rel)>pi) lat_ini=lat_ini*deg2rad

!--- FIELDS ARE PROCESSED TO BE ON STANDARD ANGULAR DOMAINS
  ALLOCATE(lon_rad(iml_rel),lat_rad(jml_rel))
  CALL conf_dat2d(title, iml_rel, jml_rel, lon_ini, lat_ini, lon_rad, lat_rad, &
                  relief_hi, .FALSE.)
  DEALLOCATE(lon_ini,lat_ini)

!--- COMPUTING THE REQUIRED FIELDS USING ROUTINE grid_noro
  IF(check) WRITE(lunout,*)'Computes all parameters needed for gravity wave dra&
     &g code'

  ALLOCATE(phis(iml,jml))      ! Geopotentiel au sol
  ALLOCATE(zstd(iml,jml))      ! Deviation standard de l'orographie sous-maille
  ALLOCATE(zsig(iml,jml))      ! Pente de l'orographie sous-maille 
  ALLOCATE(zgam(iml,jml))      ! Anisotropie de l'orographie sous maille
  ALLOCATE(zthe(iml,jml))      ! Orientation axe +grande pente d'oro sous maille
  ALLOCATE(zpic(iml,jml))      ! Hauteur pics de la SSO
  ALLOCATE(zval(iml,jml))      ! Hauteur vallees de la SSO
  ALLOCATE(relief(iml,jml))    ! Orographie moyenne

  CALL grid_noro(iml_rel, jml_rel, lon_rad, lat_rad, relief_hi, iml-1, jml,    &
       lon_in, lat_in, phis, relief, zstd, zsig, zgam, zthe, zpic, zval)
  phis = phis * g

!--- SURFACE ROUGHNESS COMPUTATION (UNUSED FOR THE MOMENT !!! )
  IF(check) WRITE(lunout,*)'Compute surface roughness induced by the orography'
  ALLOCATE(rugo   (iml  ,jml))
  ALLOCATE(tmp_var(iml-1,jml))
  CALL rugsoro(lon_rad, lat_rad, relief_hi,      &
       lon_in, lat_in, tmp_var)
  rugo(1:iml-1,:)=tmp_var; rugo(iml,:)=tmp_var(1,:)
  DEALLOCATE(relief_hi,tmp_var,lon_rad,lat_rad)
  RETURN

END SUBROUTINE start_init_orog
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_phys(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,               INTENT(IN) :: iml, jml
  REAL, DIMENSION(iml),  INTENT(IN) :: lon_in
  REAL, DIMENSION(jml),  INTENT(IN) :: lat_in
  INTEGER,               INTENT(IN) :: jml2
  REAL, DIMENSION(iml),  INTENT(IN) :: lon_in2
  REAL, DIMENSION(jml2), INTENT(IN) :: lat_in2
  LOGICAL,               INTENT(IN) :: ibar
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  CHARACTER(LEN=25)     :: title
  CHARACTER(LEN=120)    :: physfname
  LOGICAL               :: check=.TRUE.
  REAL                  :: date, dt
  INTEGER, DIMENSION(1) :: itau
  INTEGER               :: llm_tmp, ttm_tmp
  REAL,    DIMENSION(:,:), ALLOCATABLE :: var_ana
  REAL,    DIMENSION(:),   ALLOCATABLE :: lon_rad, lat_rad, lon_ini, lat_ini
  REAL,    DIMENSION(:),   ALLOCATABLE :: levphys_ini
!-------------------------------------------------------------------------------
  physfname = 'ECPHY.nc'; pi=2.0*ASIN(1.0); deg2rad=pi/180.0
  IF(check) WRITE(lunout,*)'Opening the surface analysis'
  CALL flininfo(physfname, iml_phys, jml_phys, llm_tmp, ttm_tmp, fid_phys)

  ALLOCATE(lat_phys(iml_phys,jml_phys))
  ALLOCATE(lon_phys(iml_phys,jml_phys))
  ALLOCATE(levphys_ini(llm_tmp))
  CALL flinopen(physfname, .FALSE., iml_phys, jml_phys, llm_tmp, lon_phys,     &
                lat_phys, levphys_ini, ttm_tmp, itau, date, dt, fid_phys)
  DEALLOCATE(levphys_ini)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_phys),lat_ini(jml_phys))
  lon_ini(:)=lon_phys(:,1); IF(MAXVAL(lon_phys)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_phys(1,:); IF(MAXVAL(lat_phys)>pi) lat_ini=lat_ini*deg2rad

  ALLOCATE(var_ana(iml_phys,jml_phys),lon_rad(iml_phys),lat_rad(jml_phys))

!--- SURFACE TEMPERATURE
  title='ST'
  ALLOCATE(tsol(iml,jml))
  CALL flinget(fid_phys,title,iml_phys,jml_phys,llm_tmp,ttm_tmp,1,1,var_ana)
  CALL conf_dat2d(title,iml_phys, jml_phys, lon_ini, lat_ini, lon_rad, lat_rad,&
                  var_ana , ibar  )
  CALL interp_startvar(title, ibar, .TRUE.,                                    &
      iml_phys, jml_phys, lon_rad, lat_rad, var_ana, iml, jml, jml-1,          &
      lon_in,   lat_in,   lon_in2, lat_in2, tsol)

!--- SOIL MOISTURE
  title='CDSW'
  ALLOCATE(qsol(iml,jml))
  CALL flinget(fid_phys,title,iml_phys,jml_phys,llm_tmp,ttm_tmp,1,1,var_ana)
  CALL conf_dat2d(title,iml_phys, jml_phys, lon_ini, lat_ini, lon_rad, lat_rad,&
                  var_ana, ibar  )
  CALL interp_startvar(title, ibar, .TRUE.,                                    &
      iml_phys, jml_phys, lon_rad, lat_rad, var_ana, iml, jml, jml-1,          &
      lon_in,   lat_in,   lon_in2, lat_in2, qsol)

  CALL flinclo(fid_phys)

  DEALLOCATE(var_ana,lon_rad,lat_rad,lon_ini,lat_ini)

END SUBROUTINE start_init_phys
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_init_dyn(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
!
!-------------------------------------------------------------------------------
! Arguments:
  INTEGER,               INTENT(IN) :: iml, jml
  REAL, DIMENSION(iml),  INTENT(IN) :: lon_in
  REAL, DIMENSION(jml),  INTENT(IN) :: lat_in
  INTEGER,               INTENT(IN) :: jml2
  REAL, DIMENSION(iml),  INTENT(IN) :: lon_in2
  REAL, DIMENSION(jml2), INTENT(IN) :: lat_in2
  LOGICAL,               INTENT(IN) :: ibar
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
  CHARACTER(LEN=25)     :: title
  CHARACTER(LEN=120)    :: physfname
  LOGICAL               :: check=.TRUE.
  REAL                  :: date, dt
  INTEGER, DIMENSION(1) :: itau
  INTEGER               :: i, j
  REAL,    DIMENSION(:,:), ALLOCATABLE :: var_ana, z
  REAL,    DIMENSION(:),   ALLOCATABLE :: lon_rad, lat_rad, lon_ini, lat_ini
  REAL,    DIMENSION(:),   ALLOCATABLE :: xppn, xpps
!-------------------------------------------------------------------------------

!--- KINETIC ENERGY
  physfname = 'ECDYN.nc'; pi=2.0*ASIN(1.0); deg2rad=pi/180.0
  IF(check) WRITE(lunout,*) 'Opening the surface analysis'
  CALL flininfo(physfname, iml_dyn, jml_dyn, llm_dyn, ttm_dyn, fid_dyn)
  IF(check) WRITE(lunout,*) 'Values read: ', iml_dyn, jml_dyn, llm_dyn, ttm_dyn

  ALLOCATE(lat_dyn(iml_dyn,jml_dyn))
  ALLOCATE(lon_dyn(iml_dyn,jml_dyn))
  ALLOCATE(levdyn_ini(llm_dyn))
  CALL flinopen(physfname, .FALSE., iml_dyn, jml_dyn, llm_dyn, lon_dyn,lat_dyn,&
                levdyn_ini, ttm_dyn, itau, date, dt, fid_dyn)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_dyn),lat_ini(jml_dyn))
  lon_ini(:)=lon_dyn(:,1); IF(MAXVAL(lon_dyn)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_dyn(1,:); IF(MAXVAL(lat_dyn)>pi) lat_ini=lat_ini*deg2rad

  ALLOCATE(var_ana(iml_dyn,jml_dyn),lon_rad(iml_dyn),lat_rad(jml_dyn))

!--- SURFACE GEOPOTENTIAL
  title='Z'
  ALLOCATE(z(iml,jml))
  CALL flinget(fid_dyn, title, iml_dyn, jml_dyn, 0, ttm_dyn, 1, 1, var_ana)
  CALL conf_dat2d(title, iml_dyn, jml_dyn, lon_ini, lat_ini, lon_rad, lat_rad, &
                  var_ana, ibar)
  CALL interp_startvar(title, ibar, .TRUE.,                                    &
      iml_dyn, jml_dyn, lon_rad, lat_rad, var_ana, iml, jml, jml-1,            &
      lon_in,  lat_in,  lon_in2, lat_in2, z)

!--- SURFACE PRESSURE
  title='SP'
  ALLOCATE(psol_dyn(iml,jml))
  CALL flinget(fid_dyn, title, iml_dyn, jml_dyn, 0, ttm_dyn, 1, 1, var_ana)
  CALL conf_dat2d(title, iml_dyn, jml_dyn, lon_ini, lat_ini, lon_rad, lat_rad, &
                  var_ana, ibar)
  CALL interp_startvar(title, ibar, .TRUE.,                                    &
      iml_dyn, jml_dyn, lon_rad, lat_rad, var_ana, iml, jml, jml-1,            &
      lon_in,  lat_in,  lon_in2, lat_in2, psol_dyn)

  DEALLOCATE(var_ana,lon_rad,lat_rad,lon_ini,lat_ini)

!--- ALLOCATION OF VARIABLES CREATED IN OR COMING FROM RESTART FILE
  IF(.NOT.ALLOCATED(tsol)) THEN
    CALL start_init_phys(iml,jml,lon_in,lat_in,jml2,lon_in2,lat_in2,ibar)
  ELSE
    IF(SIZE(tsol)/=SIZE(psol_dyn)) THEN
      WRITE(lunout,*)'start_init_dyn :'
      WRITE(lunout,*)'The temperature field we have does not have the right size'
      STOP
    END IF
  END IF

  IF(.NOT.ALLOCATED(phis)) THEN
    CALL start_init_orog(iml,jml,lon_in,lat_in)
  ELSE
    IF(SIZE(phis)/=SIZE(psol_dyn)) THEN
      WRITE(lunout,*)'start_init_dyn :'
      WRITE(lunout,*)'The orography field we have does not have the right size'
      STOP
    END IF
  END IF

!--- PSOL IS COMPUTED IN PASCALS
  DO j = 1, jml
    DO i = 1, iml-1
      psol_dyn(i,j) = psol_dyn(i,j)*(1.0+(z(i,j)-phis(i,j))/287.0/tsol(i,j))
    END DO
    psol_dyn(iml,j) = psol_dyn(1,j)
  END DO
  DEALLOCATE(z)

  ALLOCATE(xppn(iml-1),xpps(iml-1)) 
  DO i = 1, iml-1
    xppn(i) = aire( i,1) * psol_dyn( i,1)
    xpps(i) = aire( i,jml) * psol_dyn( i,jml)
  END DO
  DO i = 1, iml
    psol_dyn(i,1  ) = SUM(xppn)/apoln
    psol_dyn(i,jml) = SUM(xpps)/apols
  END DO
  DEALLOCATE(xppn,xpps) 

  RETURN

END SUBROUTINE start_init_dyn
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE start_inter_3d(varname, iml, jml, lml, lon_in, lat_in, jml2, &
     lon_in2, lat_in2, pls_in, var3d, ibar)

  use pchsp_95_m, only: pchsp_95
  use pchfe_95_m, only: pchfe_95

! Arguments:
  CHARACTER(LEN=*),             INTENT(IN)    :: varname
  INTEGER,                      INTENT(IN)    :: iml, jml, lml
  REAL, DIMENSION(iml),         INTENT(IN)    :: lon_in
  REAL, DIMENSION(jml),         INTENT(IN)    :: lat_in
  INTEGER,                      INTENT(IN)    :: jml2
  REAL, DIMENSION(iml),         INTENT(IN)    :: lon_in2
  REAL, DIMENSION(jml2),        INTENT(IN)    :: lat_in2
  REAL, DIMENSION(iml, jml, lml), INTENT(IN)    :: pls_in
  REAL, DIMENSION(iml, jml, lml), INTENT(OUT)   :: var3d
  LOGICAL,                      INTENT(IN)    :: ibar
!----------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  LOGICAL:: check=.TRUE., skip
  REAL                  chmin, chmax
  INTEGER ii, ij, il, ierr
  integer n_extrap ! number of extrapolated points
  REAL, DIMENSION(iml, jml, llm_dyn):: var_tmp3d
  REAL,    DIMENSION(:),     ALLOCATABLE :: lon_rad, lat_rad, lon_ini, lat_ini
  REAL, DIMENSION(llm_dyn):: lev_dyn, ax, ay, yder

!---------------------------------------------------------------------------
  IF(check) WRITE(lunout, *)'Going into flinget to extract the 3D  field.'
  IF(check) WRITE(lunout, *) fid_dyn, varname, iml_dyn, jml_dyn, llm_dyn, &
       ttm_dyn
  IF(check) WRITE(lunout, *) 'Allocating space for interpolation', iml, jml, &
       llm_dyn

  IF(.NOT.ALLOCATED(var_ana3d)) ALLOCATE(var_ana3d(iml_dyn, jml_dyn, llm_dyn))
  CALL flinget(fid_dyn, varname, iml_dyn, jml_dyn, llm_dyn, ttm_dyn, 1, 1, &
       var_ana3d)

!--- IF ANGLES ARE IN DEGREES, THEY ARE CONVERTED INTO RADIANS
  ALLOCATE(lon_ini(iml_dyn), lat_ini(jml_dyn))
  lon_ini(:)=lon_dyn(:, 1); IF(MAXVAL(lon_dyn)>pi) lon_ini=lon_ini*deg2rad
  lat_ini(:)=lat_dyn(1, :); IF(MAXVAL(lat_dyn)>pi) lat_ini=lat_ini*deg2rad

!--- FIELDS ARE PROCESSED TO BE ON STANDARD ANGULAR DOMAINS
  ALLOCATE(lon_rad(iml_dyn), lat_rad(jml_dyn))
  CALL conf_dat3d (varname, iml_dyn, jml_dyn, llm_dyn, lon_ini, lat_ini,      &
                   levdyn_ini, lon_rad, lat_rad, lev_dyn, var_ana3d, ibar)
  DEALLOCATE(lon_ini, lat_ini)

!--- COMPUTING THE REQUIRED FIELDS USING ROUTINE grid_noro
  DO il=1, llm_dyn
    CALL interp_startvar(varname, ibar, il==1, iml_dyn, jml_dyn, lon_rad, &
         lat_rad, var_ana3d(:, :, il), iml, jml, jml2, lon_in, lat_in, &
         lon_in2, lat_in2, var_tmp3d(:, :, il))
  END DO
  DEALLOCATE(lon_rad, lat_rad)

!--- VERTICAL INTERPOLATION IS PERFORMED FROM TOP OF ATMOSPHERE TO GROUND
  ax = lev_dyn(llm_dyn:1:-1) 
  skip = .false.
  n_extrap = 0
  DO ij=1, jml
    DO ii=1, iml-1
      ay = var_tmp3d(ii, ij, llm_dyn:1:-1)
      yder = pchsp_95(ax, ay, ibeg=2, iend=2, vc_beg=0., vc_end=0.)
      CALL pchfe_95(ax, ay, yder, skip, pls_in(ii, ij, lml:1:-1), &
           var3d(ii, ij, lml:1:-1), ierr)
      if (ierr < 0) stop 1
      n_extrap = n_extrap + ierr
    END DO
  END DO
  if (n_extrap /= 0) then
     print *, "start_inter_3d pchfe_95: n_extrap = ", n_extrap
  end if
  var3d(iml, :, :) = var3d(1, :, :) 

  DO il=1, lml
    CALL minmax(iml*jml, var3d(1, 1, il), chmin, chmax)
    WRITE(lunout, *)' '//TRIM(varname)//'  min max l ', il, chmin, chmax
  END DO

END SUBROUTINE start_inter_3d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE interp_startvar(vname, ibar, ibeg, ii, jj,    lon,  lat,  vari,     &
                                 i1, j1, j2, lon1, lat1, lon2, lat2, varo)
!
!-------------------------------------------------------------------------------

  USE inter_barxy_m, only: inter_barxy
  USE grid_atob_m, only: grille_m

! Arguments:
  CHARACTER(LEN=*),       INTENT(IN)  :: vname
  LOGICAL,                INTENT(IN)  :: ibar, ibeg
  INTEGER,                INTENT(IN)  :: ii, jj
  REAL, DIMENSION(ii),    INTENT(IN)  :: lon
  REAL, DIMENSION(jj),    INTENT(IN)  :: lat
  REAL, DIMENSION(ii,jj), INTENT(IN)  :: vari
  INTEGER,                INTENT(IN)  :: i1, j1, j2
  REAL, DIMENSION(i1),    INTENT(IN)  :: lon1
  REAL, DIMENSION(j1),    INTENT(IN)  :: lat1
  REAL, DIMENSION(i1),    INTENT(IN)  :: lon2
  REAL, DIMENSION(j2),    INTENT(IN)  :: lat2
  REAL, DIMENSION(i1,j1), INTENT(OUT) :: varo
!-------------------------------------------------------------------------------
! Local variables:
#include "iniprint.h"
  REAL, DIMENSION(i1-1,j1) :: vtmp
!-------------------------------------------------------------------------------
  IF(ibar) THEN
    IF(ibeg) THEN
      WRITE(lunout,*)                                                          &
               '---------------------------------------------------------------'
      WRITE(lunout,*)                                                          &
 '$$$ Utilisation de l interpolation barycentrique  pour  '//TRIM(vname)//' $$$'
      WRITE(lunout,*)                                                          &
               '---------------------------------------------------------------'
    END IF
    CALL inter_barxy(lon, lat(:jj-1), vari, lon2(:i1-1), lat2(:j2), vtmp)
  ELSE
    CALL grille_m   (lon, lat, vari, lon1, lat1, vtmp)
  END IF
  CALL gr_int_dyn(vtmp, varo, i1-1, j1)

END SUBROUTINE interp_startvar
!
!-------------------------------------------------------------------------------

END MODULE startvar
!
!*******************************************************************************
