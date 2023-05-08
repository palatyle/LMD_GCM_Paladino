










!
! $Id: dynredem.F 1635 2012-07-12 11:37:16Z lguez $
!
!
SUBROUTINE dynredem0(fichnom,iday_end,phis)
  USE infotrac, ONLY: nqtot, tname, ttext
  USE netcdf, ONLY: NF90_CREATE, NF90_DEF_DIM, NF90_INQ_VARID, NF90_GLOBAL,    &
                    NF90_CLOSE,  NF90_PUT_ATT, NF90_UNLIMITED, NF90_CLOBBER
  USE dynredem_mod, ONLY: cre_var, put_var1, put_var2, err, modname, fil
  use netcdf95, only: NF95_PUT_VAR
  use control_mod, only : planet_type
  USE comvert_mod, ONLY: ap,bp,aps,bps,presnivs,pseudoalt,pa,preff, &
      			nivsig,nivsigs
  USE comconst_mod, ONLY: daysec,dtvr,rad,omeg,g,cpp,kappa,pi
  USE logic_mod, ONLY: fxyhypb,ysinus
  USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy, &
      			taux,tauy
  USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn,itaufin, &
      			start_time,hour_ini
  USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0

  IMPLICIT NONE
!=======================================================================
! Writting the NetCDF restart file (initialisation)
!=======================================================================
!   Declarations:
!   -------------
  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "netcdf.inc"
  include "iniprint.h"

!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom          !--- FILE NAME
  INTEGER,          INTENT(IN) :: iday_end         !--- 
  REAL,             INTENT(IN) :: phis(iip1, jjp1) !--- GROUND GEOPOTENTIAL
!===============================================================================
! Local variables:
  INTEGER :: iq,l
  INTEGER, PARAMETER :: length=100
  REAL :: tab_cntrl(length) ! run parameters
  INTEGER :: ierr
  CHARACTER(LEN=80) :: abort_message

!   For NetCDF:
  CHARACTER(LEN=30) :: unites
  INTEGER :: indexID
  INTEGER :: rlonuID, rlonvID, rlatuID, rlatvID
  INTEGER :: sID, sigID, nID, vID, timID
  INTEGER :: yyears0, jjour0, mmois0
  REAL    :: zan0, zjulian, hours

  CHARACTER(len=12) :: start_file_type="earth" ! default start file type
  INTEGER :: idecal

!===============================================================================
  ! fill dynredem_mod module variables
  modname='dynredem0'; fil=fichnom

! set yyears0, mmois0, jjour0 to 0,1,1 (hours is not used)
  yyears0=0
  mmois0=1
  jjour0=1

  !!! AS: idecal is a hack to be able to read planeto starts...
  !!!     .... while keeping everything OK for LMDZ EARTH
  if ((planet_type.eq."generic").or.(planet_type.eq."mars")) then
    write(lunout,*) trim(modname),' : Planeto-like start file'
    start_file_type="planeto"
    idecal = 4
  else
    if ( planet_type.eq."titan" ) then
      ! Titan inherited Earth-like start files with idecal=5
      write(lunout,*) trim(modname),' : Titan start file'
    else
      write(lunout,*) trim(modname),' : Earth-like start file'
    endif
    idecal = 5
  endif

  tab_cntrl(:)  = 0.
  tab_cntrl(1)  = REAL(iim)
  tab_cntrl(2)  = REAL(jjm)
  tab_cntrl(3)  = REAL(llm)
  if (start_file_type.eq."earth") then
    tab_cntrl(4)=REAL(day_ref)
  else
    !tab_cntrl(4)=REAL(day_end)
    tab_cntrl(4)=REAL(iday_end)
  endif
  tab_cntrl(5)  = REAL(annee_ref)
  tab_cntrl(idecal+1)  = rad
  tab_cntrl(idecal+2)  = omeg
  tab_cntrl(idecal+3)  = g
  tab_cntrl(idecal+4)  = cpp
  tab_cntrl(idecal+5) = kappa
  tab_cntrl(idecal+6) = daysec
  tab_cntrl(idecal+7) = dtvr
  tab_cntrl(idecal+8) = etot0
  tab_cntrl(idecal+9) = ptot0
  tab_cntrl(idecal+10) = ztot0
  tab_cntrl(idecal+11) = stot0
  tab_cntrl(idecal+12) = ang0
  tab_cntrl(idecal+13) = pa
  tab_cntrl(idecal+14) = preff

!    .....    parameters for the zoom      ......   
  tab_cntrl(idecal+15)  = clon
  tab_cntrl(idecal+16)  = clat
  tab_cntrl(idecal+17)  = grossismx
  tab_cntrl(idecal+18)  = grossismy
!
  IF ( fxyhypb )   THEN
    tab_cntrl(idecal+19) = 1.
    tab_cntrl(idecal+20) = dzoomx
    tab_cntrl(idecal+21) = dzoomy
    tab_cntrl(idecal+22) = 0.
    tab_cntrl(idecal+23) = taux
    tab_cntrl(idecal+24) = tauy
  ELSE
    tab_cntrl(idecal+19) = 0.
    tab_cntrl(idecal+20) = dzoomx
    tab_cntrl(idecal+21) = dzoomy
    tab_cntrl(idecal+22) = 0.
    tab_cntrl(idecal+23) = 0.
    tab_cntrl(idecal+24) = 0.
    IF( ysinus )  tab_cntrl(idecal+22) = 1.
  ENDIF

  if (start_file_type.eq."earth") then
    tab_cntrl(idecal+25) = REAL(iday_end)
    tab_cntrl(idecal+26) = REAL(itau_dyn + itaufin)
! start_time: start_time of simulation (not necessarily 0.)
    tab_cntrl(idecal+27) = start_time
  endif

  if (planet_type=="mars") then ! For Mars only
    tab_cntrl(29)=hour_ini
  endif

!--- File creation
  CALL err(NF90_CREATE(fichnom,NF90_CLOBBER,nid))

!--- Some global attributes
  CALL err(NF90_PUT_ATT(nid,NF90_GLOBAL,"title","Fichier demarrage dynamique"))

!--- Dimensions
  if (start_file_type.eq."earth") then
    CALL err(NF90_DEF_DIM(nid,"index", length, indexID))
    CALL err(NF90_DEF_DIM(nid,"rlonu", iip1,   rlonuID))
    CALL err(NF90_DEF_DIM(nid,"rlatu", jjp1,   rlatuID))
    CALL err(NF90_DEF_DIM(nid,"rlonv", iip1,   rlonvID))
    CALL err(NF90_DEF_DIM(nid,"rlatv", jjm,    rlatvID))
    CALL err(NF90_DEF_DIM(nid,"sigs",  llm,        sID))
    CALL err(NF90_DEF_DIM(nid,"sig",   llmp1,    sigID))
    CALL err(NF90_DEF_DIM(nid,"temps", NF90_UNLIMITED, timID))
  else
    CALL err(NF90_DEF_DIM(nid,"index", length, indexID))
    CALL err(NF90_DEF_DIM(nid,"rlonu", iip1,   rlonuID))
    CALL err(NF90_DEF_DIM(nid,"latitude", jjp1,   rlatuID))
    CALL err(NF90_DEF_DIM(nid,"longitude", iip1,   rlonvID))
    CALL err(NF90_DEF_DIM(nid,"rlatv", jjm,    rlatvID))
    CALL err(NF90_DEF_DIM(nid,"altitude",  llm,        sID))
    CALL err(NF90_DEF_DIM(nid,"interlayer",   llmp1,    sigID))
    CALL err(NF90_DEF_DIM(nid,"Time", NF90_UNLIMITED, timID))
  endif

!--- Define and save invariant fields
  CALL put_var1(nid,"controle","Parametres de controle" ,[indexID],tab_cntrl)
  CALL put_var1(nid,"rlonu"   ,"Longitudes des points U",[rlonuID],rlonu)
  CALL put_var1(nid,"rlatu"   ,"Latitudes des points U" ,[rlatuID],rlatu)
  CALL put_var1(nid,"rlonv"   ,"Longitudes des points V",[rlonvID],rlonv)
  CALL put_var1(nid,"rlatv"   ,"Latitudes des points V" ,[rlatvID],rlatv)
  if (start_file_type.eq."earth") then
    CALL put_var1(nid,"nivsigs" ,"Numero naturel des couches s"    ,[sID]  ,nivsigs)
    CALL put_var1(nid,"nivsig"  ,"Numero naturel des couches sigma",[sigID],nivsig)
  endif ! of if (start_file_type.eq."earth")
  CALL put_var1(nid,"ap"      ,"Coefficient A pour hybride"      ,[sigID],ap)
  CALL put_var1(nid,"bp"      ,"Coefficient B pour hybride"      ,[sigID],bp)
  if (start_file_type.ne."earth") then
    CALL put_var1(nid,"aps","Coef AS: hybrid pressure at midlayers",[sID],aps)
    CALL put_var1(nid,"bps","Coef BS: hybrid sigma at midlayers",[sID],bps)
  endif ! of if (start_file_type.eq."earth")
  CALL put_var1(nid,"presnivs",""                                ,[sID]  ,presnivs)
  if (start_file_type.ne."earth") then
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR(nid,"latitude",NF_DOUBLE,1,rlatuID,vID)
        ierr =NF_PUT_ATT_TEXT(nid,vID,'units',13,"degrees_north")
        ierr = NF_PUT_ATT_TEXT (nid,vID,"long_name", 14, &
              "North latitude")
        ierr = NF_ENDDEF(nid)
        call NF95_PUT_VAR(nid,vID,rlatu*180/pi)
!
        ierr = NF_REDEF (nid)
        ierr=NF_DEF_VAR(nid,"longitude",NF_DOUBLE,1,rlonvID,vID)
        ierr = NF_PUT_ATT_TEXT (nid,vID,"long_name", 14, &
              "East longitude")
        ierr = NF_PUT_ATT_TEXT(nid,vID,'units',12,"degrees_east")
        ierr = NF_ENDDEF(nid)
        call NF95_PUT_VAR(nid,vID,rlonv*180/pi)
!
        ierr = NF_REDEF (nid)
        ierr = NF_DEF_VAR (nid, "altitude", NF_DOUBLE, 1, &
             sID,vID)
        ierr = NF_PUT_ATT_TEXT(nid,vID,"long_name",10,"pseudo-alt")
        ierr = NF_PUT_ATT_TEXT (nid,vID,'units',2,"km")
        ierr = NF_PUT_ATT_TEXT (nid,vID,'positive',2,"up")
        ierr = NF_ENDDEF(nid)
        call NF95_PUT_VAR(nid,vID,pseudoalt)
        CALL err(NF_REDEF(nid))
  endif ! of if (start_file_type.ne."earth")

! covariant <-> contravariant <-> natural conversion coefficients
  CALL put_var2(nid,"cu","Coefficient de passage pour U",[rlonuID,rlatuID],cu)
  CALL put_var2(nid,"cv","Coefficient de passage pour V",[rlonvID,rlatvID],cv)
  CALL put_var2(nid,"aire","Aires de chaque maille"     ,[rlonvID,rlatuID],aire)
  CALL put_var2(nid,"phisinit","Geopotentiel au sol"    ,[rlonvID,rlatuID],phis)


! Define variables that will be stored later:
  WRITE(unites,"('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')"),&
               yyears0,mmois0,jjour0
  IF (start_file_type.eq."earth") THEN
    CALL cre_var(nid,"temps","Temps de simulation",[timID],unites)
  ELSE
    CALL cre_var(nid,"Time","Temps de simulation",[timID],unites)
  ENDIF

  CALL cre_var(nid,"ucov" ,"Vitesse U"  ,[rlonuID,rlatuID,sID,timID])
  CALL cre_var(nid,"vcov" ,"Vitesse V"  ,[rlonvID,rlatvID,sID,timID])
  CALL cre_var(nid,"teta" ,"Temperature",[rlonvID,rlatuID,sID,timID])

  IF(nqtot.GE.1) THEN
    DO iq=1,nqtot
      CALL cre_var(nid,tname(iq),ttext(iq),[rlonvID,rlatuID,sID,timID])
    END DO
  ENDIF

  CALL cre_var(nid,"masse","Masse d air"    ,[rlonvID,rlatuID,sID,timID])
  CALL cre_var(nid,"ps"   ,"Pression au sol",[rlonvID,rlatuID    ,timID])

  CALL err(NF90_CLOSE (nid)) ! close file

  WRITE(lunout,*)TRIM(modname)//': iim,jjm,llm,iday_end',iim,jjm,llm,iday_end
  WRITE(lunout,*)TRIM(modname)//': rad,omeg,g,cpp,kappa',rad,omeg,g,cpp,kappa

END SUBROUTINE dynredem0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dynredem1(fichnom,time,vcov,ucov,teta,q,masse,ps)
!
!-------------------------------------------------------------------------------
! Purpose: Write the NetCDF restart file (append).
!-------------------------------------------------------------------------------
  USE infotrac, ONLY: nqtot, tname, type_trac
  USE control_mod, only : planet_type
  USE netcdf,   ONLY: NF90_OPEN,  NF90_NOWRITE, NF90_GET_VAR, NF90_INQ_VARID,  &
                      NF90_CLOSE, NF90_WRITE,   NF90_PUT_VAR, NF90_NoErr
  use netcdf95, only: NF95_PUT_VAR
  USE temps_mod, ONLY: itaufin,itau_dyn
  USE dynredem_mod, ONLY: dynredem_write_u, dynredem_write_v, dynredem_read_u, &
                          err, modname, fil, msg
 
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "netcdf.inc"
  include "comgeom.h"
  include "iniprint.h"
!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom              !-- FILE NAME
  REAL, INTENT(IN)    ::  time                         !-- TIME
  REAL, INTENT(IN)    ::  vcov(iip1,jjm, llm)          !-- V COVARIANT WIND
  REAL, INTENT(IN)    ::  ucov(iip1,jjp1,llm)          !-- U COVARIANT WIND
  REAL, INTENT(IN)    ::  teta(iip1,jjp1,llm)          !-- POTENTIAL TEMPERATURE
  REAL, INTENT(INOUT) ::     q(iip1,jjp1,llm,nqtot)    !-- TRACERS
  REAL, INTENT(IN)    :: masse(iip1,jjp1,llm)          !-- MASS PER CELL
  REAL, INTENT(IN)    ::    ps(iip1,jjp1)              !-- GROUND PRESSURE
!===============================================================================
! Local variables:
  INTEGER :: l, iq, nid, vID, ierr, nid_trac, vID_trac
  INTEGER,SAVE :: nb=0
  INTEGER, PARAMETER :: length=100
  REAL               :: tab_cntrl(length) ! tableau des parametres du run
  CHARACTER(LEN=256) :: var, dum
  LOGICAL            :: lread_inca
  CHARACTER(LEN=80) :: abort_message
  CHARACTER(len=12) :: start_file_type="earth" ! default start file type

  ! fill dynredem_mod module variables
  modname='dynredem1'; fil=fichnom

  if ((planet_type.eq."generic").or.(planet_type.eq."mars")) then
      write(lunout,*) trim(modname),' : Planeto-like start file'
      start_file_type="planeto"
  else
      write(lunout,*) trim(modname),' : Earth-like start file'
  endif 

  CALL err(NF90_OPEN(fil,NF90_WRITE,nid),"open",fil)

!--- Write/extend time coordinate
  nb = nb + 1
  if (start_file_type.eq."earth") then
        ierr = NF_INQ_VARID(nid, "temps", vID)
        IF (ierr .NE. NF_NOERR) THEN
          write(lunout,*) NF_STRERROR(ierr)
          abort_message='Variable temps n est pas definie'
          CALL abort_gcm(modname,abort_message,ierr)
        ENDIF
 else
        ierr = NF_INQ_VARID(nid,"Time", vID)
        IF (ierr .NE. NF_NOERR) THEN
          write(lunout,*) NF_STRERROR(ierr)
          abort_message='Variable Time not found'
          CALL abort_gcm(modname,abort_message,ierr)
        ENDIF
  endif ! of if (start_file_type.eq."earth") 
  call NF95_PUT_VAR(nid,vID,time,start=(/nb/))
  WRITE(lunout,*)TRIM(modname)//": Saving for ", nb, time

!--- Rewrite control table (itaufin undefined in dynredem0)
  var="controle"
  CALL err(NF90_INQ_VARID(nid,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(nid,vID,tab_cntrl),"get",var)
  if (start_file_type=="earth") then
    tab_cntrl(31) = REAL(itau_dyn + itaufin)
  else
    tab_cntrl(31) = 0
  endif
  CALL err(NF90_INQ_VARID(nid,var,vID),"inq",var)
  CALL err(NF90_PUT_VAR(nid,vID,tab_cntrl),"put",var)

!--- Save fields
  CALL dynredem_write_u(nid,"ucov" ,ucov ,llm, nb)
  CALL dynredem_write_v(nid,"vcov" ,vcov ,llm, nb)
  CALL dynredem_write_u(nid,"teta" ,teta ,llm, nb)
  CALL dynredem_write_u(nid,"masse",masse,llm, nb)
  CALL dynredem_write_u(nid,"ps"   ,ps   ,1, nb)

!--- Tracers in file "start_trac.nc" (added by Anne)
  lread_inca=.FALSE.; fil="start_trac.nc"
  IF(type_trac=='inca') INQUIRE(FILE=fil,EXIST=lread_inca)
  IF(lread_inca) CALL err(NF90_OPEN(fil,NF90_NOWRITE,nid_trac),"open")

!--- Save tracers
  IF(nqtot.GE.1) THEN
    DO iq=1,nqtot
      var=tname(iq); ierr=-1
      IF(lread_inca) THEN                  !--- Possibly read from "start_trac.nc"
        fil="start_trac.nc"
        ierr=NF90_INQ_VARID(nid_trac,var,vID_trac)
        dum='inq'; IF(ierr==NF90_NoErr) dum='fnd'
        WRITE(lunout,*)msg(dum,var)


        IF(ierr==NF90_NoErr) CALL dynredem_read_u(nid_trac,var,q(:,:,:,iq),llm)
      END IF ! of IF(lread_inca)
      fil=fichnom
      CALL dynredem_write_u(nid,var,q(:,:,:,iq),llm,nb)
    END DO ! of DO iq=1,nqtot
  ENDIF ! of IF(nqtot.GE.1)

  CALL err(NF90_CLOSE(nid),"close")
  fil="start_trac.nc"
  IF(lread_inca) CALL err(NF90_CLOSE(nid_trac),"close")

END SUBROUTINE dynredem1

