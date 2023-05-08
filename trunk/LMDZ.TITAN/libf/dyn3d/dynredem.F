      SUBROUTINE dynredem0(fichnom,idayref,phis)
      use infotrac, only: tname,nqtot
      USE comvert_mod, ONLY: ap,bp,aps,bps,pa,preff,presnivs,pseudoalt
      USE comconst_mod, ONLY: daysec,dtvr,rad,omeg,g,cpp,kappa,pi
      USE logic_mod, ONLY: fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,
     .			dzoomx,dzoomy,taux,tauy
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      IMPLICIT NONE
c=======================================================================
c Ecriture du fichier de redemarrage sous format NetCDF (initialisation)
c=======================================================================
c   Declarations:
c   -------------
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "netcdf.inc"
!#include "advtrac.h"
c   Arguments:
c   ----------
      INTEGER*4 idayref
      REAL phis(ip1jmp1)
      CHARACTER*(*) fichnom

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
      character(len=80) :: txt ! to store some text

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
      data yyears0 /1/
      data jjour0 /1/
      data mmois0 /1/
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
         WRITE(6,*)" Failed creating file "//fichnom
         WRITE(6,*)' ierr = ', ierr
         CALL ABORT
      ENDIF
c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid,NF_GLOBAL,"title",18,
     .                       "Dynamic start file")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title in file "//fichnom
        call abort
      endif
c
c Definir les dimensions du fichiers:
c
      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <index> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <rlonu> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <latitude> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <longitude> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <rlatv> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <altitude> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "interlayer", llmp1, idim_llmp1)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <interlayer> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_DEF_DIM (nid, "Time", NF_UNLIMITED, idim_tim)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining dimension <Time> ",
     &             "in file "//fichnom
        call abort
      endif


c     CHAMPS AJOUTES POUR LA VISUALISATION T,ps, etc... avec Grads ou ferret:
c     ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_rlatu)
c     ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_rlonv)
c     ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode",
     &             "in file "//fichnom
        call abort
      endif


c
c Definir et enregistrer certains champs invariants:
c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"controle",NF_DOUBLE,1,idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <controle> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Parametres de controle")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <controle> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <controle> ",
     &             "in file "//fichnom
        call abort
!      else
!       write(*,*) "dynredem0: controle(1)=",tab_cntrl(1)
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <rlonu> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points U")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <rlonu> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <rlonu> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <rlatu> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points U")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <rlatu> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <rlatu> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <rlonv> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points V")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <rlonv> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <rlonv> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <rlatv> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points V")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <rlatv> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <rlatv> ",
     &             "in file "//fichnom
        call abort
      endif
c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <ap> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 30,
     .          "Coef A: hybrid pressure levels"  )
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <ap> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <ap> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_llmp1,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_llmp1,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <bp> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 27,
     .      "Coef B: hybrid sigma levels")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <bp> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <bp> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
!      ierr = NF_REDEF (nid)
!      if (ierr.ne.NF_NOERR) then
!        write(*,*) "dynredem0: Failed to switch back to define mode"
!        call abort
!      endif
!#ifdef NC_DOUBLE
!      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_llmp1,nvarid)
!#else
!      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_llmp1,nvarid)
!#endif
!      if (ierr.ne.NF_NOERR) then
!        write(*,*) "dynredem0: Failed defining <ap> ",
!     &             "in file "//fichnom
!        call abort
!      endif
!
!      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 30,
!     .          "Coef A: hybrid pressure levels"  )
!      if (ierr.ne.NF_NOERR) then
!        write(*,*) "dynredem0: Failed writing title attribute ",
!     &             "for <ap> in file "//fichnom
!
!      ierr = NF_ENDDEF(nid)
!      if (ierr.ne.NF_NOERR) then
!        write(*,*) "dynredem0: Failed to switch out of define mode"
!        call abort
!      endif
!#ifdef NC_DOUBLE
!      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
!#else
!      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
!#endif
!      if (ierr.ne.NF_NOERR) then
!        write(*,*) "dynredem0: Failed writing <bp> ",
!     &             "in file "//fichnom
!        call abort
!      endif
c
c ----------------------
c
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <aps> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 37,
     .      "Coef AS: hybrid pressure at midlayers")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <aps> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <aps> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bps",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bps",NF_FLOAT,1,idim_llm,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <bps> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 34,
     .      "Coef BS: hybrid sigma at midlayers")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute ",
     &             "for <bps> in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <bps> ",
     &             "in file "//fichnom
        call abort
      endif

c
c ----------------------
c
      ierr = NF_REDEF (nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_llm,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <presniv> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_ENDDEF(nid)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch out of define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing <presniv> ",
     &             "in file "//fichnom
        call abort
      endif

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
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed to switch back to define mode"
        call abort
      endif

#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"Time",NF_DOUBLE,1,idim_tim,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"Time",NF_FLOAT,1,idim_tim,nvarid)
#endif
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed defining <Time> ",
     &             "in file "//fichnom
        call abort
      endif

      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     &                       "Temps de simulation")
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing title attribute",
     &             "for <Time> in file "//fichnom
        call abort
      endif


      write(unites,200)yyears0,mmois0,jjour0
200   format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30,
     .                         unites)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "dynredem0: Failed writing units attribute",
     &             "for <Time> in file "//fichnom
        call abort
      endif


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
      IF(nqtot.GE.1) THEN
         DO iq=1,nqtot
            IF (iq.GT.99) THEN
               PRINT*, "Trop de traceurs"
               CALL abort
            ELSE
!               str3(1:1)='q'
!               WRITE(str3(2:3),'(i2.2)') iq
!#ifdef NC_DOUBLE
!               ierr = NF_DEF_VAR (nid,str3,NF_DOUBLE,4,dims4,nvarid)
!#else
!               ierr = NF_DEF_VAR (nid,str3,NF_FLOAT,4,dims4,nvarid)
!#endif
!               ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
!     .                          "Traceurs "//str3)
             txt="Traceur "//trim(tname(iq))
#ifdef NC_DOUBLE
               ierr=NF_DEF_VAR(nid,tname(iq),NF_DOUBLE,4,dims4,nvarid)
#else
               ierr=NF_DEF_VAR(nid,tname(iq),NF_FLOAT,4,dims4,nvarid)
#endif
               ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",
     .                  len_trim(txt),trim(txt))
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

      write(*,*)'dynredem0: iim,jjm,llm,idayref',iim,jjm,llm,idayref
      write(*,*)'dynredem0: rad,omeg,g,cpp,kappa',
     &        rad,omeg,g,cpp,kappa
      
!      stop "dynredem0 halt"
      
      RETURN
      END

c ================================================================
c ================================================================

      SUBROUTINE dynredem1(fichnom,time,
     .                     vcov,ucov,teta,q,masse,ps)
      use infotrac, only: nqtot, tname
      IMPLICIT NONE
c=================================================================
c  Ecriture du fichier de redemarrage sous format NetCDF
c=================================================================
#include "dimensions.h"
#include "paramet.h"
#include "netcdf.inc"
#include "comgeom.h"
!#include"advtrac.h"

      INTEGER l
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL q(iip1,jjp1,llm,nqtot)
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

      IF (nqtot.GT.99) THEN
         PRINT*, "Trop de traceurs"
         CALL abort
      ENDIF
      IF(nqtot.GE.1) THEN
         DO iq=1,nqtot
!            str3(1:1)='q'
!            WRITE(str3(2:3),'(i2.2)') iq
!            ierr = NF_INQ_VARID(nid, str3, nvarid)
            ierr=NF_INQ_VARID(nid,tname(iq),nvarid)
            IF (ierr .NE. NF_NOERR) THEN
!               PRINT*, "Variable "//str3//" n est pas definie"
              PRINT*,"Variable "//trim(tname(iq))//" n est pas definie"
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

