










!
! $Id $
!
SUBROUTINE dynetat0(fichnom,vcov,ucov,teta,q,masse,ps,phis,time0)

      USE infotrac, only: tname, nqtot, zone_num, iso_indnum,&
                          iso_num, phase_num, alpha_ideal, iqiso, &
                          ok_isotopes, iqpere, tnat
      use netcdf, only: nf90_open,NF90_NOWRITE,nf90_noerr,nf90_strerror, &
                        nf90_get_var, nf90_inq_varid, nf90_inq_dimid, &
                        nf90_inquire_dimension,nf90_close

      use control_mod, only : planet_type, timestart
      USE comvert_mod, ONLY: pa,preff
      USE comconst_mod, ONLY: im,jm,lllm,daysec,dtvr, &
     			rad,omeg,g,cpp,kappa,pi
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy
      USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn, &
      			start_time,day_ini,hour_ini
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0

      IMPLICIT NONE

!=======================================================================
!
! Read initial confitions file
!
!=======================================================================

  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"
  include "iniprint.h"

!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom          !--- FILE NAME
  REAL, INTENT(OUT) ::  vcov(iip1,jjm, llm)        !--- V COVARIANT WIND
  REAL, INTENT(OUT) ::  ucov(iip1,jjp1,llm)        !--- U COVARIANT WIND
  REAL, INTENT(OUT) ::  teta(iip1,jjp1,llm)        !--- POTENTIAL TEMP.
  REAL, INTENT(OUT) ::     q(iip1,jjp1,llm,nqtot)  !--- TRACERS
  REAL, INTENT(OUT) :: masse(iip1,jjp1,llm)        !--- MASS PER CELL
  REAL, INTENT(OUT) ::    ps(iip1,jjp1)            !--- GROUND PRESSURE
  REAL, INTENT(OUT) ::  phis(iip1,jjp1)            !--- GEOPOTENTIAL
  REAL,INTENT(OUT) :: time0
!===============================================================================
!   Local Variables 
  CHARACTER(LEN=256) :: msg, var, modname
  INTEGER,PARAMETER :: length=100
  INTEGER :: iq, fID, vID, idecal
  REAL :: tab_cntrl(length) ! array containing run parameters
  INTEGER :: ierr
  CHARACTER(len=12) :: start_file_type="earth" ! default start file type

  REAL,ALLOCATABLE :: time(:) ! times stored in start
  INTEGER :: timelen ! number of times stored in the file
  INTEGER :: indextime ! index of selected time
  !REAL  hour_ini ! fraction of day of stored date. Equivalent of day_ini, but 0=<hour_ini<1

  INTEGER :: edges(4),corner(4)
  INTEGER :: i

!-----------------------------------------------------------------------
  modname="dynetat0"

!  Open initial state NetCDF file
  var=fichnom
  CALL err(NF90_OPEN(var,NF90_NOWRITE,fID),"open",var)
!
  CALL get_var1("controle",tab_cntrl)

      !!! AS: idecal is a hack to be able to read planeto starts...
      !!!     .... while keeping everything OK for LMDZ EARTH
      if ((planet_type.eq."generic").or.(planet_type.eq."mars")) then
          write(lunout,*)'dynetat0 : Planeto-like start file'
          start_file_type="planeto"
          idecal = 4
          annee_ref  = 2000
      else
          if (planet_type.eq."titan") then
             ! Titan inherited Earth-like start files with idecal=5
             write(lunout,*)'dynetat0 : Titan start file'
          else
             write(lunout,*)'dynetat0 : Earth-like start file'
          endif
          idecal = 5
          annee_ref  = tab_cntrl(5)
      endif


      im         = tab_cntrl(1)
      jm         = tab_cntrl(2)
      lllm       = tab_cntrl(3)
      if (start_file_type.eq."earth") then
        day_ref    = tab_cntrl(4)
      else
        day_ini    = tab_cntrl(4)
        day_ref=0
      endif
      rad        = tab_cntrl(idecal+1)
      omeg       = tab_cntrl(idecal+2)
      g          = tab_cntrl(idecal+3)
      cpp        = tab_cntrl(idecal+4)
      kappa      = tab_cntrl(idecal+5)
      daysec     = tab_cntrl(idecal+6)
      dtvr       = tab_cntrl(idecal+7)
      etot0      = tab_cntrl(idecal+8)
      ptot0      = tab_cntrl(idecal+9)
      ztot0      = tab_cntrl(idecal+10)
      stot0      = tab_cntrl(idecal+11)
      ang0       = tab_cntrl(idecal+12)
      pa         = tab_cntrl(idecal+13)
      preff      = tab_cntrl(idecal+14)
!
      clon       = tab_cntrl(idecal+15)
      clat       = tab_cntrl(idecal+16)
      grossismx  = tab_cntrl(idecal+17)
      grossismy  = tab_cntrl(idecal+18)
!
      IF ( tab_cntrl(idecal+19).EQ.1. )  THEN
        fxyhypb  = .TRUE.
!        dzoomx   = tab_cntrl(25)
!        dzoomy   = tab_cntrl(26)
!        taux     = tab_cntrl(28)
!        tauy     = tab_cntrl(29)
      ELSE
        fxyhypb = .FALSE.
        ysinus  = .FALSE.
        IF( tab_cntrl(idecal+22).EQ.1. ) ysinus = .TRUE. 
      ENDIF

      if (planet_type=="mars") then ! so far this is only for Mars
        hour_ini = tab_cntrl(29)
      else
        hour_ini=0
      endif

      if (start_file_type.eq."earth") then
        day_ini = tab_cntrl(30)
        itau_dyn = tab_cntrl(31)
        start_time = tab_cntrl(32)
      else
        day_ini=tab_cntrl(4)
        itau_dyn=0
        start_time=0
      endif
!   .................................................................
!
!
  WRITE(lunout,*)trim(modname)//': rad,omeg,g,cpp,kappa ', &
                     rad,omeg,g,cpp,kappa

  CALL check_dim(im,iim,'im','iim')
  CALL check_dim(jm,jjm,'jm','jjm')
  CALL check_dim(llm,lllm,'llm','lllm')

  CALL get_var1("rlonu",rlonu)
  CALL get_var1("rlatu",rlatu)
  CALL get_var1("rlonv",rlonv)
  CALL get_var1("rlatv",rlatv)

  CALL get_var2("cu"   ,cu)
  CALL get_var2("cv"   ,cv)

  CALL get_var2("aire" ,aire)
  CALL get_var2("phisinit",phis)

! read time axis
      ierr = nf90_inq_varid (fID, "temps", vID)
      IF (ierr .NE. nf90_noerr) THEN
        write(lunout,*)"dynetat0: Le champ <temps> est absent"
        write(lunout,*)"dynetat0: J essaie <Time>"
        ierr = nf90_inq_varid (fID, "Time", vID)
        IF (ierr .NE. nf90_noerr) THEN
           write(lunout,*)"dynetat0: Le champ <Time> est absent"
           write(lunout,*)trim(nf90_strerror(ierr))
           CALL ABORT_gcm("dynetat0", "", 1)
        ENDIF
        ! Get the length of the "Time" dimension
        ierr = nf90_inq_dimid(fID,"Time",vID)
        ierr = nf90_inquire_dimension(fID,vID,len=timelen)
        allocate(time(timelen))
        ! Then look for the "Time" variable
        ierr  =nf90_inq_varid(fID,"Time",vID)
        ierr = nf90_get_var(fID, vID, time)
        IF (ierr .NE. nf90_noerr) THEN
           write(lunout,*)"dynetat0: Lecture echouee <Time>"
           write(lunout,*)trim(nf90_strerror(ierr))
           CALL ABORT_gcm("dynetat0", "", 1)
        ENDIF
      ELSE   
        ! Get the length of the "temps" dimension
        ierr = nf90_inq_dimid(fID,"temps",vID)
        ierr = nf90_inquire_dimension(fID,vID,len=timelen)
        allocate(time(timelen))
        ! Then look for the "temps" variable
        ierr = nf90_inq_varid (fID, "temps", vID)
        ierr = nf90_get_var(fID, vID, time)
        IF (ierr .NE. nf90_noerr) THEN
           write(lunout,*)"dynetat0: Lecture echouee <temps>"
           write(lunout,*)trim(nf90_strerror(ierr))
           CALL ABORT_gcm("dynetat0", "", 1)
        ENDIF
      ENDIF

! select the desired time
      IF (timestart .lt. 0) THEN  ! default: we use the last time value
        indextime = timelen
      ELSE  ! else we look for the desired value in the time axis
       indextime = 0
        DO i=1,timelen
          IF (abs(time(i) - timestart) .lt. 0.01) THEN
             indextime = i
             EXIT
          ENDIF
        ENDDO
        IF (indextime .eq. 0) THEN
          write(lunout,*)"Time", timestart," is not in " &
                                            //trim(fichnom)//"!!"
          write(lunout,*)"Stored times are:"
          DO i=1,timelen
             PRINT*, time(i)
          ENDDO
          CALL ABORT_gcm("dynetat0", "", 1)
        ENDIF
      ENDIF

      if (planet_type=="mars") then
        ! In start the absolute date is day_ini + hour_ini + time
        ! For now on, in the GCM dynamics, it is day_ini + time0
        time0 = time(indextime) + hour_ini
        day_ini = day_ini + INT(time0)
        time0 = time0 - INT(time0) ! time0 devient le nouveau hour_ini
        hour_ini = time0
      else
        time0 = time(indextime) 
      endif
      
      PRINT*, "dynetat0: Selected time ",time(indextime), &
              " at index ",indextime
      
      DEALLOCATE(time)

! read vcov
  CALL get_var3v_t("vcov",vcov,indextime)

! read ucov
  CALL get_var3u_t("ucov",ucov,indextime)
 
! read teta (same corner/edges as ucov)
  CALL get_var3u_t("teta",teta,indextime)

! read tracers (same corner/edges as ucov)
  corner(1)=1
  corner(2)=1
  corner(3)=1
  corner(4)=indextime
  edges(1)=iip1
  edges(2)=jjp1
  edges(3)=llm
  edges(4)=1
  IF(nqtot.GE.1) THEN
      DO iq=1,nqtot
        ierr= nf90_inq_varid(fID,tname(iq),vID)
        IF (ierr .NE. nf90_noerr) THEN
           write(lunout,*)TRIM(modname)//": Tracer <"//TRIM(var)//"> is missing"
           write(lunout,*)"         It is hence initialized to zero"
           q(:,:,:,iq)=0.
           IF (planet_type=="earth") THEN
            !--- CRisi: for isotops, theoretical initialization using very simplified
            !           Rayleigh distillation las.
            IF(ok_isotopes.AND.iso_num(iq)>0) THEN
             IF(zone_num(iq)==0) q(:,:,:,iq)=q(:,:,:,iqpere(iq))*tnat(iso_num(iq))    &
             &             *(q(:,:,:,iqpere(iq))/30.e-3)**(alpha_ideal(iso_num(iq))-1)
             IF(zone_num(iq)==1) q(:,:,:,iq)=q(:,:,:,iqiso(iso_indnum(iq),phase_num(iq)))
            END IF
           ENDIF
        ELSE
           ierr=nf90_get_var(fID,vID,q(:,:,:,iq),corner,edges)
          IF (ierr .NE. nf90_noerr) THEN
            write(lunout,*)"dynetat0: Lecture echouee pour " &
                                      //trim(tname(iq))
            write(lunout,*)trim(nf90_strerror(ierr))
            CALL ABORT_gcm("dynetat0", "", 1)
          ENDIF
        ENDIF
      ENDDO
  ENDIF

!read masse (same corner/edges as ucov)
  CALL get_var3u_t("masse",masse,indextime)

!read ps
  CALL get_var2_t("ps",ps,indextime)

  CALL err(NF90_CLOSE(fID),"close",fichnom)

  if (planet_type/="mars") then
    day_ini=day_ini+INT(time0) ! obsolete stuff ; 0<time<1 anyways
    time0=time0-INT(time0)
  endif


  CONTAINS

SUBROUTINE check_dim(n1,n2,str1,str2)
  INTEGER,          INTENT(IN) :: n1, n2
  CHARACTER(LEN=*), INTENT(IN) :: str1, str2
  CHARACTER(LEN=256) :: s1, s2
  IF(n1/=n2) THEN
    s1='value of '//TRIM(str1)//' ='
    s2=' read in starting file differs from parametrized '//TRIM(str2)//' ='
    WRITE(msg,'(10x,a,i4,2x,a,i4)')TRIM(s1),n1,TRIM(s2),n2
    CALL ABORT_gcm(TRIM(modname),TRIM(msg),1)
  END IF
END SUBROUTINE check_dim


SUBROUTINE get_var1(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var1


SUBROUTINE get_var2(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var2

SUBROUTINE get_var2_t(var,v,indextime)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:)
  INTEGER, INTENT(IN) :: indextime
  corner(1)=1
  corner(2)=1
  corner(3)=indextime
  edges(1)=iip1
  edges(2)=jjp1
  edges(3)=1
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v,corner,edges),"get",var)
END SUBROUTINE get_var2_t


SUBROUTINE get_var3(var,v) ! on U grid
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var3

SUBROUTINE get_var3u_t(var,v,indextime) ! on U grid
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:)
  INTEGER, INTENT(IN) :: indextime
  corner(1)=1
  corner(2)=1
  corner(3)=1
  corner(4)=indextime
  edges(1)=iip1
  edges(2)=jjp1
  edges(3)=llm
  edges(4)=1
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v,corner,edges),"get",var)
END SUBROUTINE get_var3u_t

SUBROUTINE get_var3v_t(var,v,indextime) ! on V grid
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:)
  INTEGER, INTENT(IN) :: indextime
  corner(1)=1
  corner(2)=1
  corner(3)=1
  corner(4)=indextime
  edges(1)=iip1
  edges(2)=jjm
  edges(3)=llm
  edges(4)=1
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v,corner,edges),"get",var)
END SUBROUTINE get_var3v_t

SUBROUTINE err(ierr,typ,nam)
  INTEGER,          INTENT(IN) :: ierr   !--- NetCDF ERROR CODE
  CHARACTER(LEN=*), INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), INTENT(IN) :: nam    !--- FIELD/FILE NAME
  IF(ierr==NF90_NoERR) RETURN
  SELECT CASE(typ)
    CASE('inq');   msg="Field <"//TRIM(nam)//"> is missing"
    CASE('get');   msg="Reading failed for <"//TRIM(nam)//">"
    CASE('open');  msg="File opening failed for <"//TRIM(nam)//">"
    CASE('close'); msg="File closing failed for <"//TRIM(nam)//">"
  END SELECT
  CALL ABORT_gcm(TRIM(modname),TRIM(msg),ierr)
END SUBROUTINE err

END SUBROUTINE dynetat0
