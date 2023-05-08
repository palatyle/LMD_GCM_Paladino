      SUBROUTINE dynredem0(fichnom,idayref,anneeref,phis,nq)
      IMPLICIT NONE
c=======================================================================
c Ecriture du fichier de redemarrage sous format NetCDF (initialisation)
c=======================================================================
c   Declarations:
c   -------------
#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom.h"
#include "temps.h"
#include "ener.h"
#include "logic.h"
#include "netcdf.inc"
#include "description.h"
#include "serre.h"

c   Arguments:
c   ----------
      INTEGER*4 idayref,anneeref
      REAL phis(ip1jmp1)
      CHARACTER*(*) fichnom
      INTEGER nq

c   Local:
c   ------
      INTEGER iq,l
      CHARACTER str3*3
      INTEGER length
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr
      character*20 modname
      character*80 abort_message

c   Variables locales pour NetCDF:
c
      INTEGER dims2(2), dims3(3), dims4(4)
      INTEGER idim_index
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_llm, idim_llmp1
      INTEGER idim_tim
      INTEGER nid,nvarid

      REAL zan0,zjulian,hours
      REAL sigs(llm)
      INTEGER yyears0,jjour0, mmois0
      character*30 unites


c-----------------------------------------------------------------------
      modname='dynredem'
      do l=1,llm
         sigs(l)=float(l)
      enddo

      DO l=1,length
       tab_cntrl(l) = 0.
      ENDDO
       tab_cntrl(1)  = FLOAT(iim)
       tab_cntrl(2)  = FLOAT(jjm)
       tab_cntrl(3)  = FLOAT(llm)
       tab_cntrl(4)  = FLOAT(idayref)
       tab_cntrl(5)  = rad
       tab_cntrl(6)  = omeg
       tab_cntrl(7)  = g
       tab_cntrl(8)  = cpp
       tab_cntrl(9) = kappa
       tab_cntrl(10) = daysec
       tab_cntrl(11) = dtvr
       tab_cntrl(12) = etot0
       tab_cntrl(13) = ptot0
       tab_cntrl(14) = ztot0
       tab_cntrl(15) = stot0
       tab_cntrl(16) = ang0
       tab_cntrl(17) = pa
       tab_cntrl(18) = preff
c
c    .....    parametres  pour le zoom      ......   

       tab_cntrl(19)  = clon
       tab_cntrl(20)  = clat
       tab_cntrl(21)  = grossismx
       tab_cntrl(22)  = grossismy
c
      IF ( fxyhypb )   THEN
       tab_cntrl(23) = 1.
       tab_cntrl(24) = dzoomx
       tab_cntrl(25) = dzoomy
       tab_cntrl(26) = 0.
       tab_cntrl(27) = taux
       tab_cntrl(28) = tauy
      ELSE
       tab_cntrl(23) = 0.
       tab_cntrl(24) = dzoomx
       tab_cntrl(25) = dzoomy
       tab_cntrl(26) = 0.
       tab_cntrl(27) = 0.
       tab_cntrl(28) = 0.
       IF( ysinus )  tab_cntrl(26) = 1. 
      ENDIF
c
c    .........................................................
c
c Creation du fichier:
c
      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
         WRITE(6,*)" Pb d ouverture du fichier "//fichnom
         WRITE(6,*)' ierr = ', ierr
         CALL ABORT
      ENDIF
c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid,NF_GLOBAL,"title",18,
     .                       "Dynamic start file")
c
c Definir les dimensions du fichiers:
c
      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
      ierr = NF_DEF_DIM (nid, "interlayer", llmp1, idim_llmp1)
      ierr = NF_DEF_DIM (nid, "Time", NF_UNLIMITED, idim_tim)

c     CHAMPS AJOUTES POUR LA VISUALISATION T,ps, etc... avec Grads ou ferret:
c     ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
c     ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
c     ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
c
c Definir et enregistrer certains champs invariants:
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"controle",NF_DOUBLE,1,idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Parametres de controle")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 30,
     .          "Coef A: hybrid pressure levels"  )
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 27,
     .      "Coef B: hybrid sigma levels")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 37,
     .      "Coef AS: hybrid pressure in midlayers")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 31,
     .      "Coef BS: hybrid sigma midlayers")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_llm,nvarid)
#endif
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
c ------------------------------------------------------------------
c ------------------------------------------------------------------
c  Variable uniquement pour visualisation avec Grads ou Ferret
c ------------------------------------------------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"latitude",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"latitude",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
      ierr =NF_PUT_ATT_TEXT(nid,nvarid,'units',13,"degrees_north")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "North latitude")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu/pi*180)
#endif
c ----------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid,"longitude", NF_DOUBLE, 1, idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR(nid,"longitude", NF_FLOAT, 1, idim_rlonv,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "East longitude")
      ierr = NF_PUT_ATT_TEXT(nid,nvarid,'units',12,"degrees_east")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv/pi*180)
#endif
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "altitude", NF_DOUBLE, 1,
     .       idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "altitude", NF_FLOAT, 1,
     .       idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name",10,"pseudo-alt")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'units',2,"km")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'positive',2,"up")
  
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,pseudoalt)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pseudoalt)
#endif
  
  
c ----------------------
c ----------------------
c
c Coefficients de passage cov. <-> contra. <--> naturel
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonu
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cu",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cu",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
c
c ----------------------
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatv
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cv",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cv",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
c
c ----------------------
c Aire de chaque maille:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aire",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aire",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Aires de chaque maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif
c
c ----------------------
c Geopentiel au sol:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"phisinit",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Geopotentiel au sol")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
#endif
c
c ----------------------
c Definir les variables pour pouvoir les enregistrer plus tard:
c
      ierr = NF_REDEF (nid) ! entrer dans le mode de definition
c
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"Time",NF_DOUBLE,1,idim_tim,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"Time",NF_FLOAT,1,idim_tim,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Temps de simulation")
      write(unites,200)yyears0,mmois0,jjour0
200   format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30,
     .                         unites)

c
      dims4(1) = idim_rlonu
      dims4(2) = idim_rlatu
      dims4(3) = idim_llm
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ucov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ucov",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse U")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatv
      dims4(3) = idim_llm
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"vcov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"vcov",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse V")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_llm
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"teta",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"teta",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 11,
     .                       "Temperature")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_llm
      dims4(4) = idim_tim
      IF(nq.GE.1) THEN
         DO iq=1,nq
            IF (iq.GT.99) THEN
               PRINT*, "Trop de traceurs"
               CALL abort
            ELSE
               str3(1:1)='q'
               WRITE(str3(2:3),'(i2.2)') iq
#ifdef NC_DOUBLE
               ierr = NF_DEF_VAR (nid,str3,NF_DOUBLE,4,dims4,nvarid)
#else
               ierr = NF_DEF_VAR (nid,str3,NF_FLOAT,4,dims4,nvarid)
#endif
               ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
     .                          "Traceurs "//str3)
            ENDIF
         ENDDO
      ENDIF
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_llm
      dims4(4) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"masse",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"masse",NF_FLOAT,4,dims4,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
     .                       "C est quoi ?")
c
      dims3(1) = idim_rlonv
      dims3(2) = idim_rlatu
      dims3(3) = idim_tim
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ps",NF_DOUBLE,3,dims3,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 15,
     .                       "Pression au sol")
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
      ierr = NF_CLOSE(nid) ! fermer le fichier

      PRINT*,'iim,jjm,llm,idayref',iim,jjm,llm,idayref
      PRINT*,'rad,omeg,g,cpp,kappa',
     ,        rad,omeg,g,cpp,kappa

      RETURN
      END

c ================================================================
      SUBROUTINE dynredem1(fichnom,time,
     .                     vcov,ucov,teta,q,nq,masse,ps)
      IMPLICIT NONE
c=================================================================
c  Ecriture du fichier de redemarrage sous format NetCDF
c=================================================================
#include "dimensions.h"
#include "paramet.h"
#include "description.h"
#include "netcdf.inc"
#include "comvert.h"
#include "comgeom.h"

      INTEGER nq, l
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL q(iip1,jjp1,llm,nqmx)
      REAL q3d(iip1,jjp1,llm) !temporary variable
      CHARACTER*(*) fichnom
     
      REAL time
      INTEGER nid, nvarid
      INTEGER ierr
      INTEGER iq
      CHARACTER str3*3
      character*20 modname
      character*80 abort_message
c
      INTEGER nb,i,j
      SAVE nb
      DATA nb / 0 /

      modname = 'dynredem1'
      ierr = NF_OPEN(fichnom, NF_WRITE, nid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Pb. d ouverture "//fichnom
         CALL abort
      ENDIF

c  Ecriture/extension de la coordonnee temps

      nb = nb + 1
      ierr = NF_INQ_VARID(nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         print *, NF_STRERROR(ierr)
         abort_message='Variable Time n est pas definie'
         CALL abort_gcm(modname,abort_message,ierr)
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR1_DOUBLE (nid,nvarid,nb,time)
#else
      ierr = NF_PUT_VAR1_REAL (nid,nvarid,nb,time)
#endif
      PRINT*, "Enregistrement pour ", nb, time

c  Ecriture des champs
c
      ierr = NF_INQ_VARID(nid, "ucov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable ucov n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ucov)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ucov)
#endif

      ierr = NF_INQ_VARID(nid, "vcov", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable vcov n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,vcov)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,vcov)
#endif

      ierr = NF_INQ_VARID(nid, "teta", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable teta n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,teta)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,teta)
#endif

      IF (nq.GT.99) THEN
         PRINT*, "Trop de traceurs"
         CALL abort
      ENDIF
      IF(nq.GE.1) THEN
         DO iq=1,nq
            str3(1:1)='q'
            WRITE(str3(2:3),'(i2.2)') iq
            ierr = NF_INQ_VARID(nid, str3, nvarid)
            IF (ierr .NE. NF_NOERR) THEN
               PRINT*, "Variable "//str3//" n est pas definie"
               CALL abort
            ENDIF
            do l=1,llm
               do j=1,jjp1
                  do i=1,iip1
                     q3d(i,j,l)=q(i,j,l,iq)
                  enddo
               enddo
            enddo
#ifdef NC_DOUBLE
            ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,q3d)
#else
            ierr = NF_PUT_VAR_REAL (nid,nvarid,q3d)
#endif
            IF (ierr .NE. NF_NOERR) THEN
               PRINT*, "Error: ", NF_STRERROR(ierr)
               CALL abort
            ENDIF
         ENDDO
      ENDIF
c
      ierr = NF_INQ_VARID(nid, "masse", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable masse n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,masse)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,masse)
#endif
c
      ierr = NF_INQ_VARID(nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Variable ps n est pas definie"
         CALL abort
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ps)
#endif

      ierr = NF_CLOSE(nid)
c
      RETURN
      END

