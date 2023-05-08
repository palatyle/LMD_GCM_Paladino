










MODULE dynredem_mod

  USE netcdf, ONLY: NF90_NOERR, NF90_DOUBLE, NF90_STRERROR, &
                    NF90_REDEF, NF90_ENDDEF, NF90_PUT_VAR, &
                    NF90_PUT_ATT, NF90_GET_VAR, NF90_INQ_VARID, &
                    NF90_DEF_VAR
  PRIVATE
  PUBLIC :: dynredem_write_u, dynredem_write_v, dynredem_read_u, err
  PUBLIC :: cre_var, get_var1, put_var1, put_var2, fil, modname, msg
  include "dimensions.h"
  include "paramet.h"
  CHARACTER(LEN=256), SAVE :: fil, modname
  INTEGER,            SAVE :: nvarid


CONTAINS


!===============================================================================
!
SUBROUTINE dynredem_write_u(ncid,id,var,ll,nb)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN) :: ncid
  CHARACTER(LEN=*), INTENT(IN) :: id
  REAL,             INTENT(IN) :: var(iip1,jjp1,ll)
  INTEGER,          INTENT(IN) :: ll
  INTEGER,          INTENT(IN) :: nb
!===============================================================================
! Local variables:
  INTEGER :: start(4), count(4)
!===============================================================================
  IF (ll.eq.1) THEN
    start(:)=[1,1,nb,1]
    count(:)=[iip1,jjp1,1,1]
  ELSE
    start(:)=[1,1,1,nb]
    count(:)=[iip1,jjp1,ll,1]
  ENDIF
  CALL err(NF90_INQ_VARID(ncid,id,nvarid),"inq",id)
  CALL err(NF90_PUT_VAR(ncid,nvarid,var,start,count),"put",id)
  
END SUBROUTINE dynredem_write_u
!
!===============================================================================


!===============================================================================
!
SUBROUTINE dynredem_write_v(ncid,id,var,ll,nb)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN) :: ncid
  CHARACTER(LEN=*), INTENT(IN) :: id
  REAL,             INTENT(IN) :: var(iip1,jjm,ll)
  INTEGER,          INTENT(IN) :: ll
  INTEGER,          INTENT(IN) :: nb
!===============================================================================
! Local variables:
  INTEGER :: start(4), count(4)
!===============================================================================
  IF (ll.eq.1) THEN
    start(:)=[1,1,nb,1]
    count(:)=[iip1,jjm,1,1]
  ELSE
    start(:)=[1,1,1,nb]
    count(:)=[iip1,jjm,ll,1]
  ENDIF
  CALL err(NF90_INQ_VARID(ncid,id,nvarid),"inq",id)
  CALL err(NF90_PUT_VAR(ncid,nvarid,var,start,count),"put",id)
  
END SUBROUTINE dynredem_write_v
!
!===============================================================================


!===============================================================================
!
SUBROUTINE dynredem_read_u(ncid,id,var,ll)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,          INTENT(IN)  :: ncid
  CHARACTER(LEN=*), INTENT(IN)  :: id
  REAL,             INTENT(OUT) :: var(iip1,jjp1,ll)
  INTEGER,          INTENT(IN)  :: ll
!===============================================================================
! Local variables:
  INTEGER :: start(4), count(4)
!===============================================================================
  start(:)=[1,1,1,1]; count(:)=[iip1,jjp1,ll,1]
  CALL err(NF90_INQ_VARID(ncid,id,nvarid),"inq",id)
  CALL err(NF90_GET_VAR(ncid,nvarid,var,start,count),"get",id)
  
END SUBROUTINE dynredem_read_u    
!
!===============================================================================


!===============================================================================
!
SUBROUTINE cre_var(ncid,var,title,did,units)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ncid
  CHARACTER(LEN=*),           INTENT(IN) :: var, title
  INTEGER,                    INTENT(IN) :: did(:)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
!===============================================================================
  CALL err(NF90_DEF_VAR(ncid,var,NF90_DOUBLE,did,nvarid),"inq",var)
  IF(title/="")      CALL err(NF90_PUT_ATT(ncid,nvarid,"title",title),var)
  IF(PRESENT(units)) CALL err(NF90_PUT_ATT(ncid,nvarid,"units",units),var)

END SUBROUTINE cre_var
!
!===============================================================================


!===============================================================================
!
SUBROUTINE put_var1(ncid,var,title,did,v,units)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ncid
  CHARACTER(LEN=*),           INTENT(IN) :: var, title
  INTEGER,                    INTENT(IN) :: did(1)
  REAL,                       INTENT(IN) :: v(:)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
!===============================================================================
  IF(     PRESENT(units)) CALL cre_var(ncid,var,title,did,units)
  IF(.NOT.PRESENT(units)) CALL cre_var(ncid,var,title,did)
  CALL err(NF90_ENDDEF(ncid))
  CALL err(NF90_PUT_VAR(ncid,nvarid,v),"put",var)
  CALL err(NF90_REDEF(ncid))

END SUBROUTINE put_var1
!
!===============================================================================


!===============================================================================
!
SUBROUTINE put_var2(ncid,var,title,did,v,units)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ncid
  CHARACTER(LEN=*),           INTENT(IN) :: var, title
  INTEGER,                    INTENT(IN) :: did(2)
  REAL,                       INTENT(IN) :: v(:,:)
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
!===============================================================================
  IF(     PRESENT(units)) CALL cre_var(ncid,var,title,did,units)
  IF(.NOT.PRESENT(units)) CALL cre_var(ncid,var,title,did)
  CALL err(NF90_ENDDEF(ncid))
  CALL err(NF90_PUT_VAR(ncid,nvarid,v),"put",var)
  CALL err(NF90_REDEF(ncid))

END SUBROUTINE put_var2
!
!===============================================================================


!===============================================================================
!
FUNCTION msg(typ,nam)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  CHARACTER(LEN=256)                     :: msg    !--- STANDARDIZED MESSAGE
  CHARACTER(LEN=*),           INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nam    !--- FIELD NAME
!===============================================================================
  SELECT CASE(typ)
    CASE('open');  msg="Opening failed for <"//TRIM(fil)//">"
    CASE('close'); msg="Closing failed for <"//TRIM(fil)//">"
    CASE('get');   msg="Reading failed for <"//TRIM(nam)//">"
    CASE('put');   msg="Writting failed for <"//TRIM(nam)//">"
    CASE('inq');   msg="Missing field <"//TRIM(nam)//">"
    CASE('fnd');   msg="Found field <"//TRIM(nam)//">"
  END SELECT
  msg=TRIM(msg)//" in file <"//TRIM(fil)//">"

END FUNCTION msg
!
!===============================================================================


!===============================================================================
!
SUBROUTINE err(ierr,typ,nam)
!
!===============================================================================
  IMPLICIT NONE
!===============================================================================
! Arguments:
  INTEGER,                    INTENT(IN) :: ierr   !--- NetCDF ERROR CODE
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nam    !--- FIELD NAME
!===============================================================================
  IF(ierr==NF90_NoERR) RETURN
  IF(.NOT.PRESENT(typ)) THEN
    CALL ABORT_gcm(modname,NF90_STRERROR(ierr),ierr)
  ELSE
    CALL ABORT_gcm(modname,msg(typ,nam),ierr)
  END IF

END SUBROUTINE err
!
!===============================================================================

END MODULE dynredem_mod   
