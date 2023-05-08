      PROGRAM anldiag_nc
      IMPLICIT NONE
c======================================================================
c
c
c Program designed to read output files from the MARS gcm
c or the Mars Climate database.
c
c    Francois Forget 1999
c    Version used by monica to interpolate in zaeroid and p coordinate

c
c=======================================================================
c   declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comdissip.h"
#include "comvert.h"
#include "comgeom2.h"
#include "logic.h"
#include "temps.h"
#include "control.h"
#include "ener.h"
#include "description.h"
#include "netcdf.inc"

      INTEGER itau,nbpas,nbpasmx, coor,sor,varoutsize
      integer varspecsize
      integer varout(6),varspec(nqmx)
      integer itau0
      PARAMETER(nbpasmx=1000000)
      REAL temps(nbpasmx),y
      INTEGER unitlec
      INTEGER i,j,l,irec,itautot,ierr,iq,n,ff,k
      real, dimension(:),allocatable :: lat,lon,alt
      integer:: nout,latdimout,londimout,altdimout,timedimout,timevarout
      integer:: latdim,londim,altdim,dimid
      integer:: latvar,lonvar,altvar,timevar
      integer:: latlen,lonlen,altlen,timelen
      integer, dimension(4) :: edges,corner,id
      integer, dimension(3) :: edges2,corner2

      INTEGER   unit,nvarid,nid,varid
      REAL mugaz

      REAL missing
      PARAMETER(missing=1E+20)
      REAL valid_range(2)
      DATA valid_range /0., 300.0/
c   variables meteo dans la grille verticale GCM
c   --------------------------------------------
      REAL v(iip1,jjp1,llm),u(iip1,jjp1,llm),w(iip1,jjp1,llm)
      REAL t(iip1,jjp1,llm),rho(iip1,jjp1,llm),phi(iip1,jjp1,llm)
      REAL ps(iip1,jjp1) , tsurf(iip1,jjp1)
      REAL Rnew(iip1,jjp1,llm)
      real euv(iip1,jjp1,llm),con(iip1,jjp1,llm),nir(iip1,jjp1,llm)
      real nspec(iip1,jjp1,llm,nqmx),spec(iip1,jjp1,llm)

      REAL v_sd(iip1,jjp1,llm),u_sd(iip1,jjp1,llm)
      real w_sd(iip1,jjp1,llm),t_sd(iip1,jjp1,llm)
      real rho_sd(iip1,jjp1,llm)
      real euv_sd(iip1,jjp1,llm),con_sd(iip1,jjp1,llm)
      real nir_sd(iip1,jjp1,llm)
      real nspec_sd(iip1,jjp1,llm,nqmx)

c   ALtitude of the GCM levels
c   -------------------------
      REAL zlocal(iip1,jjp1,llm) ! altitude above the local surface (m)
      REAL zareoid(iip1,jjp1,llm) ! altitude above the mola areoid (m)
      real tmean   ! use to integrate hydrostatic equation
      real Rmean   ! use to integrate hydrostatic equation

c    Also used to integrate hydrostatic equation:
      real g0 ! exact mean gravity at radius=3396.km (LEMOINE ET AL 2001) 
      data    g0 /3.7257964/
      data    a0 /3396.0/  ! radius=3396.km
      real gz  ! used to estimate gravity in altitude...


c   variables meteo en coordonnee de pression
c   --------------------------------------------
      REAL vp(iip1,jjp1,llm),up(iip1,jjp1,llm),wp(iip1,jjp1,llm)
      REAL tp(iip1,jjp1,llm),rhop(iip1,jjp1,llm)
      REAL nspecp(iip1,jjp1,llm,nqmx)
      REAL phip(iip1,jjp1,llm)
      
      REAL vp_sd(iip1,jjp1,llm),up_sd(iip1,jjp1,llm)
      real wp_sd(iip1,jjp1,llm),tp_sd(iip1,jjp1,llm)
      real rhop_sd(iip1,jjp1,llm),nspecp_sd(iip1,jjp1,llm,nqmx)
 
c     Niveaux de pression (Pa)
      real pref(llm)

c    Version 32 layers : ************
c     data pref/
c     &1000.,884., 783. , 610., 475, 370, 288, 224, 174, 140.,
c     &115 , 80.6481, 48.0445, 26.6881, 14.1462, 7.29102, 3.70004,
c     &1.86241, 0.933516, 0.466924, 0.233296, 0.116503, 5.81631E-02,
c     & 2.90335E-02,1.44919E-02, 7.23328E-03, 3.61027E-03, 1.80193E-03,
c     & 8.99364E-04,4.48881E-04, 2.15769E-04, 1.3E-04/

c    Version 50 layers : ************
      data pref/
     &1400., 1100., 800., 565., 400., 280., 200., 140., 100.,
     &70., 50., 35., 25., 18., 13., 9., 6.,
     &3.63227,1.803712,0.8956969,0.4447900,0.2208762,0.1096839,
     &5.4467395E-02,2.7047679E-02,1.3431421E-02,
     &6.6698161E-03,3.3121328E-03,
     &1.6447564E-03,8.1676187E-04,4.0559194E-04,2.0141099E-04,
     &1.0001774E-04,4.9667346E-05,2.4664072E-05,1.2247816E-05,
     &6.0820857E-06,3.0202741E-06,1.4998238E-06,7.4479044E-07,
     &3.6985199E-07,1.8366305E-07,9.1204370E-08,4.5290751E-08,
     &2.2490720E-08,1.1168562E-08,5.5461431E-09,2.7541311E-09,
     &1.3676602E-09,6.7915956E-10/

c   variables meteo en coordonnee de zareoid
c   --------------------------------------------
      REAL va(iip1,jjp1,llm),ua(iip1,jjp1,llm),wa(iip1,jjp1,llm)
      REAL ta(iip1,jjp1,llm),rhoa(iip1,jjp1,llm),ra(iip1,jjp1,llm)
      real euva(iip1,jjp1,llm),cona(iip1,jjp1,llm),nira(iip1,jjp1,llm)
      real nspeca(iip1,jjp1,llm,nqmx)

      REAL va_sd(iip1,jjp1,llm),ua_sd(iip1,jjp1,llm)
      real wa_sd(iip1,jjp1,llm),ta_sd(iip1,jjp1,llm)
      real rhoa_sd(iip1,jjp1,llm)
      real euva_sd(iip1,jjp1,llm),cona_sd(iip1,jjp1,llm)
      real nira_sd(iip1,jjp1,llm)
      real nspeca_sd(iip1,jjp1,llm,nqmx)

      REAL vij(llm),uij(llm),zaij(llm),wij(llm),rij(llm)
      REAL tij(llm),rhoij(llm+1),sigij(llm)
      real euvij(llm),conij(llm),nirij(llm)
      real nspecij(llm+1,nqmx)

      REAL vij_sd(llm),uij_sd(llm)
      real wij_sd(llm),tij_sd(llm),rhoij_sd(llm+1)
      real euvij_sd(llm),conij_sd(llm),nirij_sd(llm)
      real nspecij_sd(llm+1,nqmx)

c     Niveaux de zareoid (m)
      real zaref(llm)
      real p_zaref,p_zaij(llm+1) ! used to interpolate rho

       data zaref/    
      & -6.E+3,-4.E+3,-2.E+3,0.E+3,2.E+3,4.E+3,6.E+3,8.E+3,10.E+3,12.E+3
      & ,14.E+3,16.E+3,18.E+3,22.E+3,26.E+3,30.E+3,34.E+3
      & ,38.E+3,44.E+3,50.E+3,56.E+3,62.E+3,68.E+3,74.E+3
      & ,80.E+3,86.E+3,92.E+3,98.E+3,104.E+3,110.E+3,116.E+3
      & ,122.E+3,128.E+3,134.E+3,140.E+3,146.E+3,152.E+3,158.E+3
      & ,164.E+3,170.E+3,176.E+3,182.E+3,188.E+3,194.E+3,200.E+3
      & ,206.E+3,212.E+3,218.E+3,224.E+3,230.E+3/


c   Surface data (sometime needed for some analysis):
c   ------------------------------------------------
      character relief*3
      real phis(iip1,jjp1)
      real alb(iip1,jjp1),ith(iip1,jjp1)
      real zmeaS(iip1,jjp1),zstdS(iip1,jjp1)
      real zsigS(iip1,jjp1),zgamS(iip1,jjp1),ztheS(iip1,jjp1)


      real utime
      real localtime(iip1)
      common/temporaire/localtime

      INTEGER*4 day0
      integer nmoy(jjp1,llm),tp1,tp2,pass

      CHARACTER*1 yes,yescomp

      integer iformat
       
c   declarations de l'interface avec mywrite:
c   -----------------------------------------

      CHARACTER file*80
      CHARACTER nomfich*60
      CHARACTER nomfile(2)*60,nomfileout(2)*60
      integer fichsize,fstats
      CHARACTER fdiag*60

c   externe:
c   --------

      EXTERNAL iniconst,inigeom,covcont,mywrite
      EXTERNAL inifilr,exner
      EXTERNAL solarlong,coordij,moy2
      EXTERNAL SSUM
      REAL SSUM
      EXTERNAL lnblnk
      INTEGER lnblnk

c   Dust
c   ----

#include "fxyprim.h"

      character*10  noms(nqmx),varnoms(6)
c      data          noms/"co2","co","o","o1d","o2","o3","h","h2",
c     $                   "oh","ho2","h2o2", "n2", "h2o"/
c      data          varnoms/"T","U","V","W","RHO","Species"/


c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

      unitlec=11
      itautot=0


c     Lecture du fichier a lire
c     -------------------------

      print*,' ' 
      PRINT*,'File to process: 1)DIAGFI 2)STATS 3)BOTH'
      READ(5,'(i)') fichsize
      if(fichsize.ne.2)then
        i=1
        PRINT*,'   enter name of diagfi file (w/o .nc)'
        READ(5,'(a)') nomfile(i)
        fdiag=nomfile(i)
        PRINT*,'   do you want to rename the output file?  
     &  1) Yes  2) No'
        READ(5,'(i)') sor
        if (sor.eq.1) then
          PRINT*,'   enter label for output file'
          READ(5,'(a)') nomfileout(i)
        elseif(sor.eq.2) then
          nomfileout(i)=nomfile(i)
        endif
      endif
      if(fichsize.ne.1)then
        i=i+1
        PRINT*,'   enter name of stats file (w/o .nc)'
        READ(5,'(a)') nomfile(i)
        PRINT*,'   do you want to rename the output file?  
     &  1) Yes  2) No'
        READ(5,'(i)') sor
        if (sor.eq.1) then
          PRINT*,'   enter label for output file'
          READ(5,'(a)') nomfileout(i)
        elseif(sor.eq.2) then
          nomfileout(i)=nomfile(i)
        endif
        if(fichsize.eq.2) then
          PRINT*,'   enter name of diagfi file for controle (w/o .nc)'
          read(5,'(a)') fdiag
        endif
        fstats=1
      endif
      if (fichsize.ge.2) fichsize=fichsize-1
      print*,' '
      print*,'Output in: 1) netcdf 2) IDL 3) both ?'
      READ(5,'(i)') sor
      print*,' '
      print*,'Coordinates in: 1) pressure 2) Zareoid 3) both  ?'
      READ(5,'(i)') coor

c ************* Partie a moderniser !!! FF 2004 ************
      print*,' '
      print*,'Fields:  T  U  V  W  Rho  Species  All   ?'
      print*,'Type total number of desired outputs or 6 for All   '
      READ(5,'(i)') varoutsize
      if(varoutsize.lt.6) then
        print*,' '
        print*,'Fields:  T  U  V  W  Rho  Species '
        print*,'  #   :  1  2  3  4   5      6'
        print*,' '
        do i=1,varoutsize
          print*,'Type number of variable ',i
          READ(5,'(i)') varout(i)
        enddo
      endif
      if(varoutsize.eq.6) then
        do i=1,6
          varout(i)=i
        enddo
      endif
      do i=1,varoutsize
        if (varout(i).eq.6) then
           print*,'Type total number of desired species or 13 for All'
           READ(5,'(i2)') varspecsize
           if(varspecsize.ne.13) then
             print*,' '
             print*,'Select Species: CO2 CO O O(1D) O2 O3' 
             print*,'                 1  2  3   4   5  6'
             print*,'                 H H2 OH HO2 H2O2 N2 H2O'
             print*,'                 7 8  9  10   11  12 13'
             do n=1,varspecsize
               print*,'Type number of species ',i
               read(5,'(i)') varspec(n)
             enddo
           else
             do n=1,varspecsize
               varspec(n)=n
             enddo
           endif
         endif
       enddo
      if(varoutsize.gt.6) print*,'Value must be <= 6'
     
c ***************** Fin partie a moderniser ***************

! LOOP ON FILES
! =================================================================
      do ff=1,fichsize
 
      file=nomfile(ff)  
      nomfich=nomfile(ff)  
      file=file(1:lnblnk(file))
      PRINT*,'file',file

c     Ouverture fichiers :
c     ------------------
      write(*,*) "opening "//trim(nomfich)//"..."
      ierr = NF_OPEN(trim(adjustl(nomfich))//".nc",NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening file '//nomfich
         stop ""
      endif

! Control & lecture on dimensions
! =================================================================
      ierr=NF_INQ_DIMID(nid,"latitude",latdim)
      ierr=NF_INQ_VARID(nid,"latitude",latvar)
      if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <latitude> is missing'
      stop ""
      endif
      ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
!  write(*,*) "latlen: ",latlen

      ierr=NF_INQ_DIMID(nid,"longitude",londim)
      ierr=NF_INQ_VARID(nid,"longitude",lonvar)
      if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <longitude> is lacking'
      stop ""
      endif
      ierr=NF_INQ_DIMLEN(nid,londim,lonlen)
!  write(*,*) "lonlen: ",lonlen

      ierr=NF_INQ_DIMID(nid,"altitude",altdim)
      ierr=NF_INQ_VARID(nid,"altitude",altvar)
      if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <altitude> is lacking'
      stop ""
      endif

      ierr=NF_INQ_DIMLEN(nid,altdim,altlen)
!  write(*,*) "altlen: ",altlen

      if (ff.eq.1) then
        allocate(lat(latlen))
        allocate(lon(lonlen))
        allocate(alt(altlen))
      endif

#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid,latvar,lat)
      ierr = NF_GET_VAR_DOUBLE(nid,lonvar,lon)
      ierr = NF_GET_VAR_DOUBLE(nid,altvar,alt)
#else
      ierr = NF_GET_VAR_REAL(nid,latvar,lat)
      ierr = NF_GET_VAR_REAL(nid,lonvar,lon)
      ierr = NF_GET_VAR_REAL(nid,altvar,alt)
#endif

c   Lecture Time :

      ierr= NF_INQ_DIMID (nid,"Time",dimid)
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'anl_NC: Le champ <Time> est absent'
          CALL abort
        ENDIF

      ierr= NF_INQ_DIMLEN (nid,dimid,nbpas)
      ierr = NF_INQ_VARID (nid, "Time", nvarid)
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, temps)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, temps)
#endif
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'anl_NC: Lecture echouee pour <Time>'
          CALL abort
        ENDIF
        PRINT*,'temps',(temps(itau),itau=1,10)
! =================================================================

! SETTING OUTPUT FILES

      file=nomfileout(ff)  
      nomfich=nomfileout(ff)  
      file=file(1:lnblnk(file))

c    layers in zareoid : ************
          do l=1,llm
            alt(l)=zaref(l)/1000.
          enddo

! =================================================================
! Output in netcdf
! =================================================================
      if (sor.ne.2) then

        if(coor.ne.2) then 					! PRESSURE COORDINATES
 
          ierr = NF_CREATE(trim(adjustl(nomfich))//"_P.nc",
     &                     NF_CLOBBER, nout)
          if (ierr.NE.NF_NOERR) THEN
            write(*,*)' Pb d''ouverture du fichier '
            write(*,*)' ierr = ', ierr
            STOP
          ENDIF
          ierr = NF_PUT_ATT_TEXT(nout, NF_GLOBAL, "title", 24,
     &                          "Pressure coordinates ")
          ierr = NF_DEF_DIM(nout,"latitude", size(lat), latdimout)
          ierr = NF_DEF_DIM(nout,"longitude", size(lon), londimout)
          ierr = NF_DEF_DIM(nout,"altitude", size(pref), altdimout)
          ierr = NF_DEF_DIM(nout,"Time", NF_UNLIMITED, timedimout)
          ierr = NF_ENDDEF(nout)

        endif

        if(coor.ne.1) then 					! ZAEROID COORDINATES 

          ierr = NF_CREATE(trim(adjustl(nomfich))//"_Z.nc",
     &                     NF_CLOBBER, nout)
          IF(ierr.NE.NF_NOERR) THEN
            write(*,*)' Pb d''ouverture du fichier '
            write(*,*)' ierr = ', ierr
            STOP
          ENDIF
          ierr = NF_PUT_ATT_TEXT(nout, NF_GLOBAL, "title", 24,
     &                           "zaeroid coordinates")
          ierr = NF_DEF_DIM(nout,"latitude", size(lat), latdimout)
          ierr = NF_DEF_DIM(nout,"longitude", size(lon), londimout)
          ierr = NF_DEF_DIM(nout,"altitude ", size(pref), altdimout)
          ierr = NF_DEF_DIM(nout,"Time", NF_UNLIMITED, timedimout)
          ierr = NF_ENDDEF(nout)

        endif

        call def_var(nout,"Time","Time","days since 0000-00-0
     &  00:00:00",1,(/timedimout/),timevarout,ierr)

        ierr = NF_REDEF (nout)       
        call def_var(nout,"latitude","latitude","degrees_north",1,
     &               (/latdimout/),nvarid,ierr)

#ifdef NC_DOUBLE
        ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,lat)
#else
        ierr = NF_PUT_VAR_REAL (nout,nvarid,lat)
#endif
        ierr = NF_REDEF (nout) 
        call def_var(nout,"longitude","East longitude","degrees_east",1,
     &               (/londimout/),nvarid,ierr)
#ifdef NC_DOUBLE
        ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,lon)
#else
        ierr = NF_PUT_VAR_REAL (nout,nvarid,lon)
#endif
        ierr = NF_REDEF (nout)

        if(coor.ne.2) then 					! PRESSURE CCOORDINATES
#ifdef NC_DOUBLE
          ierr = NF_DEF_VAR (nout,"altitude",NF_DOUBLE,1,
     &                       (/altdimout/), nvarid)
#else
          ierr = NF_DEF_VAR (nout,"altitude",NF_FLOAT,1,
     &                       (/altdimout/),nvarid)
#endif

          ierr = NF_PUT_ATT_TEXT (nout,nvarid,"long_name",22,"altitude
     &                            pression")
          ierr = NF_PUT_ATT_TEXT (nout,nvarid,'units',2,"Pa")
          ierr = NF_PUT_ATT_TEXT (nout,nvarid,'positive',4,"down")
          ierr = NF_ENDDEF(nout)
#ifdef NC_DOUBLE
          ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,pref)
#else
          ierr = NF_PUT_VAR_REAL (nout,nvarid,pref)
#endif
        endif

        if(coor.ne.1) then 					! ZAEROID COORDINATES 

#ifdef NC_DOUBLE
          ierr = NF_DEF_VAR (nout,"altitude",NF_DOUBLE,1,
     &                       (/altdimout/), nvarid)
#else
          ierr = NF_DEF_VAR (nout,"altitude",NF_FLOAT,1,
     &                       (/altdimout/),nvarid)
#endif
          ierr = NF_PUT_ATT_TEXT (nout,nvarid,"long_name",20,"altitude
     &                           zaeroid")
          ierr = NF_PUT_ATT_TEXT (nout,nvarid,'units',2,"km")
          ierr = NF_PUT_ATT_TEXT (nout,nvarid,'positive',2,"up")
          ierr = NF_ENDDEF(nout)
#ifdef NC_DOUBLE
          ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,alt)
#else
          ierr = NF_PUT_VAR_REAL (nout,nvarid,alt)
#endif
        endif

      endif !end if netcdf  

! Output in idl
! =================================================================
      if (sor.ne.1) then

        if(coor.ne.2) then 					! PRESSURE COORDINATES
          do i=1,varoutsize
            j=121+i
            if(varout(i).eq.6) then
              do iq=1,varspecsize
                file=trim(nomfileout(ff))//'_P_'
     &               //trim(noms(iq))//'.dat'  
                j=121+i+(iq-1)
                open(j,file=trim(file), status='unknown')
                write(j,117) noms(iq) 
                write(j,120) iim,jjp1,llm
                write(j,119) lon
                write(j,119) lat
                write(j,119) pref
              enddo
            else
              file=trim(nomfileout(ff))//'_P_'
     &            //trim(varnoms(varout(i)))//'.dat'  
              open(j,file=trim(file), status='unknown')
              write(j,117) varnoms(varout(i)) 
              write(j,120) iim,jjp1,llm
              write(j,119) lon
              write(j,119) lat
              write(j,119) pref
            endif
          enddo 
        endif

        if(coor.ne.1) then 					! ZAEROID COORDINATES
          do i=1,varoutsize
            j=1221+i
            if(varout(i).eq.6) then
              do iq=1,varspecsize
                file=trim(nomfileout(ff))//'_Z_'
     &               //trim(noms(iq))//'.dat'  
                j=1221+i+(iq-1)
                open(j,file=trim(file), status='unknown')
                write(j,117) noms(iq)
                write(j,120) iim,jjp1,llm
                write(j,119) lon
                write(j,119) lat
                write(j,119) alt
              enddo
            else
              file=trim(nomfileout(ff))//'_Z_'
     &            //trim(varnoms(varout(i)))//'.dat'  
              open(j,file=trim(file), status='unknown')
              write(j,117) varnoms(varout(i))
              write(j,120) iim,jjp1,llm
              write(j,119) lon
              write(j,119) lat
              write(j,119) alt
            endif
          enddo
              
c              file=trim(nomfileout(ff))
c              j=j+1
c              open(j,file=file(1:lnblnk(file))//'_Z_EUV.dat',
c     &               status='unknown')
c              write(j,117) 'EUV'
c              write(j,120) iim,jjp1,llm
c              write(j,119) lon
c              write(j,119) lat
c              write(j,119) alt
c              j=j+1
c              open(j,file=file(1:lnblnk(file))//'_Z_CON.dat',
c     &               status='unknown')
c              write(j,117) 'CONDUCTION'
c              write(j,120) iim,jjp1,llm
c              write(j,119) lon
c              write(j,119) lat
c              write(j,119) alt
c              j=j+1
c              open(j,file=file(1:lnblnk(file))//'_Z_NIR.dat',
c     &               status='unknown')
c              write(j,117) 'NIR'
c              write(j,120) iim,jjp1,llm
c              write(j,119) lon
c              write(j,119) lat
c              write(j,119) alt
        endif

      endif
120    format(3i3)
119    format(12(e10.4,1x))
118    format(2i3)
117    format(6(a7))
116    format(13(a6))
! =================================================================
 

800   continue ! LOOP SUR LES FICHIERS


! =================================================================
! INITIALIZATION OF PARAMETERS
! =================================================================

      rad=3397200.              ! rayon de mars (m)  ~3397200 m
      daysec=88775.             ! duree du sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec) ! vitesse de rotation (rad.s-1)
      g=3.72                    ! gravite (m.s-2) ~3.72
      mugaz=43.49               ! Masse molaire de l'atm (g.mol-1) ~43.49
      kappa=.256793             ! = r/cp  ~0.256793
      r = 191.18213
      pi=2.*asin(1.)
      ecritphy =1
      iphysiq=1
      day_ini=0.
      day_step=1

      CALL readhead_NC(trim(adjustl(fdiag))//".nc",day0,phis,R)
      WRITE (*,*) 'day0 = ' , day0

      CALL iniconst
      CALL inigeom
      CALL inifilr

c     If needed : getting topography, albedo, thermal inertia...
c     -------------------------------------------------------
c     "phis" is the surface geopotential = topography (m) * g 
       relief='mola'    ! Topography MOLA par defaut

        CALL datareadnc(relief,phis,alb,ith,zmeaS,zstdS,zsigS,zgamS,
      .          ztheS)


! =================================================================
! LOOP IN TIME
! =================================================================

      DO itau=1,nbpas

! FILE READING : 
! =================================================================
        print*,'Timestep = ',itau,'/',nbpas

c Cas special stats pour la base de donnee : FF 2004
c        call lect_var(nid,'tsurf',3,ierr,itau,	tsurf)
        call lect_var(nid,'ps',3,ierr,itau,ps)

        call lect_var(nid,'temp',4,ierr,itau,t)
        call lect_var(nid,'u',4,ierr,itau,u)
        call lect_var(nid,'v',4,ierr,itau,v)
        call lect_var(nid,'w',4,ierr,itau,		w)
        call lect_var(nid,'rho',4,ierr,itau,	rho)
        


        do n=1,varspecsize
          i=varspec(n)
        call lect_var(nid,'n_'//trim(noms(i)),4,ierr,itau,spec)
        do l=1,llm
          do j=1,jjp1
            do i=1,iip1
              nspec(i,j,l,n)=spec(i,j,l)
            enddo
          enddo
        enddo
        enddo

        if(fstats*ff-fichsize .eq. 0) then
          call lect_var(nid,'temp_sd',4,ierr,itau,	t_sd)
          call lect_var(nid,'u_sd',4,ierr,itau,		u_sd)
          call lect_var(nid,'v_sd',4,ierr,itau,		v_sd)
c          call lect_var(nid,'w_sd',4,ierr,itau,		w_sd)
c          call lect_var(nid,'rho_sd',4,ierr,itau,	rho_sd)
          do n=1,varspecsize
            i=varspec(n)
            call lect_var(nid,'n_'//trim(noms(i))//'_sd',4,ierr,
     &                    itau,	spec)
            do l=1,llm
              do j=1,jjp1
                do i=1,iip1
                  nspec_sd(i,j,l,n)=spec(i,j,l)
                enddo
              enddo
            enddo
          enddo
        endif

! =================================================================
! COMPUTING ALTITUDE OF GCM LEVELS IN GCM GRID
! =================================================================

c  ---------------------------------------------------------------
c   Computing Molecular Gas Constant at every grid point (J kg-1 K-1)
c   later used to integrate the hydrostatic equation
c  ---------------------------------------------------------------
        do l=1,llm
          do j=1,jjp1
            do i=1,iip1
c              Rnew = P/(rho*T)
              Rnew(i,j,l) =  (aps(l) + bps(l)*ps(i,j))
     &                      / ( rho(i,j,l) * T(i,j,l) )
            enddo
          enddo
        enddo


c      ALTITUDE OF GCM LEVELS 
c      ----------------------
c      zlocal 	! altitude above the local surface (m)
c      zareoid 	! altitude above the mola areoid (m)
c      phis(i,j)/g ! altitude of the local surface above the mola areoid (m)

        do j=1,jjp1
          do i=1,iip1
            sigij(1)=aps(1)/ps(i,j)+bps(1)
            zlocal(i,j,1)=-log(sigij(1))* Rnew(i,j,1)*t(i,j,1)/g
            zareoid(i,j,1) = zlocal(i,j,1) + phis(i,j)/g0

            do l=2,llm
              sigij(l)=aps(l)/ps(i,j)+bps(l)
              tmean=t(i,j,l)
              if(t(i,j,l).ne.t(i,j,l-1))
     &          tmean=(t(i,j,l)-t(i,j,l-1))/log(t(i,j,l)/t(i,j,l-1))
                Rmean = 0.5*(Rnew(i,j,l) + Rnew(i,j,l-1))
                gz= g0*a0**2/(a0+zareoid(i,j,l-1))**2

              zlocal(i,j,l)=zlocal(i,j,l-1) 
     &                -log(sigij(l)/sigij(l-1))*Rmean*tmean/gz
              zareoid(i,j,l) = zlocal(i,j,l) + phis(i,j)/g
            enddo
  
          enddo
        enddo


c Local time
c ----------
c universal time (in stats oe MCD files)
        utime = (itau-1)*2.

c local time (in stats file)
        DO i = 1,iim
          localtime(i)=utime + 12.*rlonv(i)/pi
          if(localtime(i).gt.24) localtime(i)=localtime(i)-24.
          if(localtime(i).lt.0) localtime(i)=localtime(i)+24.
        ENDDO

c c Passage en niveaux de pression
c ------------------------------
        if(coor.ne.2) then 					! PRESSURE COORDINATES
          call xsig2p (ps,T,pref,llm,Tp)
          call xsig2p (ps,u,pref,llm,up)
          call xsig2p (ps,v,pref,llm,vp)
          call xsig2p (ps,w,pref,llm,wp)
          call xsig2p (ps,rho,pref,llm,rhop)
          do iq=1,varspecsize
            call xsig2p (ps,nspec(1,1,1,iq),pref,llm,nspecp(1,1,1,iq))
          enddo
           if(fstats*ff-fichsize .eq. 0) then
             call xsig2p (ps,u_sd,pref,llm,up_sd)
             call xsig2p (ps,v_sd,pref,llm,vp_sd)
             call xsig2p(ps,t_sd,pref,llm,Tp_sd)
             call xsig2p (ps,w_sd,pref,llm,wp_sd)
             call xsig2p (ps,rho_sd,pref,llm,rhop_sd)
             do iq=1,varspecsize
               call xsig2p (ps,nspec_sd(1,1,1,iq),pref,llm,
     $                      nspecp_sd(1,1,1,iq))
             enddo
           endif
         endif

c c Passage en niveaux de zareoid
c ------------------------------
        if(coor.ne.1) then 					! ZAEROID COORDINATES
        do j=1,jjp1
          do i=1,iip1
            do l=1,llm
              Tij(l)=T(i,j,l)
              uij(l)=u(i,j,l)
              vij(l)=v(i,j,l)
              wij(l)=w(i,j,l)
              rhoij(l)=rho(i,j,l)
              do iq=1,varspecsize
                nspecij(l,iq)=nspec(i,j,l,iq)
                if(l.eq.llm) nspecij(l+1,iq)=nspecij(l,iq)*1.e-3
              enddo
c              euvij(l)=euv(i,j,l)
c              conij(l)=con(i,j,l)
c              nirij(l)=nir(i,j,l)
              zaij(l)=zareoid(i,j,l)
              p_zaij(l)=exp(-zareoid(i,j,l)/1.e4) ! used to interpolate rho
              if(fstats*ff-fichsize .eq. 0) then
                Tij_sd(l)=T_sd(i,j,l)
                uij_sd(l)=u_sd(i,j,l)
                vij_sd(l)=v_sd(i,j,l)
                wij_sd(l)=w_sd(i,j,l)
                rhoij_sd(l)=rho_sd(i,j,l)
                do iq=1,varspecsize
                  nspecij_sd(l,iq)=nspec_sd(i,j,l,iq)
                  if(l.eq.llm) nspecij_sd(l+1,iq)=nspecij_sd(l,iq)*1.e-3
                enddo
              endif
            enddo

            do l=1,llm
              call interpolf(zaref(l),y,zaij,Tij,llm)
              Ta(i,j,l)=y
              call interpolf(zaref(l),y,zaij,uij,llm)
              ua(i,j,l)=y
              call interpolf(zaref(l),y,zaij,vij,llm)
              va(i,j,l)=y
              call interpolf(zaref(l),y,zaij,wij,llm)
              wa(i,j,l)=y
c              call interpolf(zaref(l),y,zaij,euvij,llm)
c              euva(i,j,l)=y
c              call interpolf(zaref(l),y,zaij,conij,llm)
c              cona(i,j,l)=y
c              call interpolf(zaref(l),y,zaij,nirij,llm)
c              nira(i,j,l)=y

              p_zaref=exp(-zaref(l)/1.e4)
              call interpolf(p_zaref,y,p_zaij,rhoij,llm+1)
              rhoa(i,j,l)=y
               
              do iq=1,varspecsize
                call interpolf(p_zaref,y,p_zaij,nspecij(1,iq),llm+1)
c                call interpolf(zaref(l),y,zaij,nspecij(1,iq),llm)
                nspeca(i,j,l,iq)=y
              enddo

              if(fstats*ff-fichsize .eq. 0) then
                call interpolf(zaref(l),y,zaij,Tij_sd,llm)
                Ta_sd(i,j,l)=y
                call interpolf(zaref(l),y,zaij,uij_sd,llm)
                ua_sd(i,j,l)=y
                call interpolf(zaref(l),y,zaij,vij_sd,llm)
                va_sd(i,j,l)=y
                call interpolf(zaref(l),y,zaij,wij_sd,llm)
                wa_sd(i,j,l)=y
                call interpolf(p_zaref,y,p_zaij,rhoij_sd,llm+1)
                rhoa_sd(i,j,l)=y
                do iq=1,varspecsize
                  call interpolf(p_zaref,y,p_zaij,nspecij_sd(1,iq),
     &                           llm+1)
c                  call interpolf(zaref(l),y,zaij,nspecij_sd(1,iq),llm)
                  nspeca_sd(i,j,l,iq)=y
                enddo
              endif

            enddo
          enddo
        enddo

c            do l=1,llm
c        print*,T(32,24,l),Ta(32,24,l)
c        print*,nspec(32,24,l,1),nspeca(32,24,l,1)
c            enddo

        endif
 
c Output at each time-step, in netcdf format
! =================================================================
        if (sor .ne. 2) then

          ierr= NF_INQ_VARID(nout,"Time",nvarid)
#ifdef NC_DOUBLE
          ierr= NF_PUT_VARA_DOUBLE(nout,nvarid,itau,1,temps(itau))
#else
          ierr= NF_PUT_VARA_REAL(nout,nvarid,itau,1,temps(itau))
#endif
          if (ierr.ne.NF_NOERR) then
            write(*,*) "time problem ",NF_STRERROR(ierr)
            stop
          endif

          if (varoutsize.eq.6) n=1
          if (varoutsize.ne.6) n=varoutsize
          do i=1,n
            if(varout(i).eq.1 .or. varoutsize.eq.6) then
             if(coor.ne.1) 
     &          call put_var(nout,'temp','temp','K',4,ierr,itau,Ta)
             if(coor.ne.2) 
     &          call put_var(nout,'temp','temp','K',4,ierr,itau,Tp)
              if(fstats*ff-fichsize.eq.0) then
                if(coor.ne.1) call put_var(nout,'temp_sd','temp_sd','K',
     &                                      4,ierr,itau,Ta_sd)
                if(coor.ne.2) call put_var(nout,'temp_sd','temp_sd','K',
     &                                      4,ierr,itau,Tp_sd)
              endif
            endif
            if(varout(i).eq.2 .or. varoutsize.eq.6)  then
              if(coor.ne.1) 
     &          call put_var(nout,'u','u','m s-1',4,ierr,itau,ua)
              if(coor.ne.2) 
     &          call put_var(nout,'u','u','m s-1',4,ierr,itau,up)
              if(fstats*ff-fichsize.eq.0) then
                if(coor.ne.1) call put_var(nout,'u_sd','u_sd','m s-1',
     &                                      4,ierr,itau,ua_sd)
                if(coor.ne.2) call put_var(nout,'u_sd','u_sd','m s-1',
     &                                      4,ierr,itau,up_sd)
              endif
            endif
            if(varout(i).eq.3 .or. varoutsize.eq.6) then
              if(coor.ne.1) 
     &          call put_var(nout,'v','v','m s-1',4,ierr,itau,va)
              if(coor.ne.2) 
     &          call put_var(nout,'v','v','m s-1',4,ierr,itau,vp)
              if(fstats*ff-fichsize.eq.0) then
                if(coor.ne.1) call put_var(nout,'v_sd','v_sd','m s-1',
     &                                      4,ierr,itau,va_sd)
                if(coor.ne.2) call put_var(nout,'v_sd','v_sd','m s-1',
     &                                      4,ierr,itau,vp_sd)
              endif
            endif
            if(varout(i).eq.4 .or. varoutsize.eq.6) then
              if(coor.ne.1) 
     &          call put_var(nout,'w','w','m s-1',4,ierr,itau,wa)
              if(coor.ne.2) 
     &          call put_var(nout,'w','w','m s-1',4,ierr,itau,wp)
              if(fstats*ff-fichsize.eq.0) then
                if(coor.ne.1) call put_var(nout,'w_sd','w_sd','m s-1',
     &                                      4,ierr,itau,wa_sd)
                if(coor.ne.2) call put_var(nout,'w_sd','w_sd','m s-1',
     &                                      4,ierr,itau,wp_sd)
              endif
            endif
            if(varout(i).eq.5 .or. varoutsize.eq.6)  then
              if(coor.ne.1) call put_var(nout,'rho','rho','',4,ierr,
     &                                    itau,rhoa)
              if(coor.ne.2) call put_var(nout,'rho','rho','',4,ierr,
     &                                    itau,rhop)
              if(fstats*ff-fichsize.eq.0) then
                if(coor.ne.1) call put_var(nout,'rho_sd','rho_sd','',
     &                                      4,ierr,itau,rhoa_sd)
                if(coor.ne.2) call put_var(nout,'rho_sd','rho_sd','',
     &                                      4,ierr,itau,rhop_sd)
              endif
            endif
            if(varout(i).eq.6 .or. varoutsize.eq.6) then
              do iq=1,varspecsize
                j=varspec(iq)
                if(coor.ne.1) call put_var(nout,'n_'//trim(noms(j)),
     &       'n_'//trim(noms(j)),'cm-3',4,ierr,itau,nspeca(1,1,1,iq))
                if(coor.ne.2) call put_var(nout,'n_'//trim(noms(j)),
     &         'n_'//trim(noms(j)),'cm-3',4,ierr,itau,nspecp(1,1,1,iq))
                if(fstats*ff-fichsize.eq.0) then
                  if(coor.ne.1) call put_var(nout,
     &                        'n_'//trim(noms(j))//'_sd',
     &                        'n_'//trim(noms(j))//'_sd','cm-3',4,ierr,
     &                           itau,nspeca_sd(1,1,1,iq))
                  if(coor.ne.2) call put_var(nout,
     &                        'n_'//trim(noms(j))//'_sd',
     &                        'n_'//trim(noms(j))//'_sd','cm-3',4,ierr,
     &                           itau,nspecp_sd(1,1,1,iq))
                endif
              enddo
            endif     
          enddo
        endif !endif writing data on Netcdf file 

c Output at each time-step, ASCII file for IDL
! =================================================================
c In pressure and zareoid coordinates:

      if (sor .ne. 1) then
       do n=1,varoutsize
         if(varout(n).eq.1) then
           do l=1,llm
             if(coor.ne.2) then
               k=121+n
               write(k,124) ((Tp(i,j,l),i=1,iim),j=1,jjp1)
             elseif(coor.ne.1) then
               k=1221+n
               write(k,124) ((Ta(i,j,l),i=1,iim),j=1,jjp1)
             endif
           enddo
         endif
         if(varout(n).eq.2) then
           do l=1,llm
             if(coor.ne.2) then
               k=121+n
               write(k,124) ((up(i,j,l),i=1,iim),j=1,jjp1)
             elseif(coor.ne.1)then
               k=1221+n
               write(k,124) ((ua(i,j,l),i=1,iim),j=1,jjp1)
             endif
           enddo
         endif 
         if(varout(n).eq.3) then
           do l=1,llm
             if(coor.ne.2) then
               k=121+n
               write(k,124) ((vp(i,j,l),i=1,iim),j=1,jjp1)
             elseif(coor.ne.1) then
               k=1221+n
               write(k,124) ((va(i,j,l),i=1,iim),j=1,jjp1)
             endif
           enddo
         endif
         if(varout(n).eq.4) then
           do l=1,llm
             if(coor.ne.2) then
               k=121+n
               write(122,124) ((wp(i,j,l),i=1,iim),j=1,jjp1)
             elseif(coor.ne.1) then
               k=1221+n
               write(k,124) ((wa(i,j,l),i=1,iim),j=1,jjp1)
             endif
           enddo
         endif
         if(varout(n).eq.5) then
           do l=1,llm
            if(coor.ne.2) then
              k=121+n
              write(k,124) ((rhop(i,j,l),i=1,iim),j=1,jjp1)
            elseif(coor.ne.1) then
               k=1221+n
              write(k,124) ((rhoa(i,j,l),i=1,iim),j=1,jjp1)
             endif
           enddo
         endif
         if(varout(n).eq.6) then
           do iq=1,varspecsize
             do l=1,llm
               if(coor.ne.1) then
                 k=1221+n+(iq-1)
                 write(k,124) ((nspeca(i,j,l,iq),i=1,iim),j=1,jjp1)
               endif
             enddo
           enddo
         endif     
       enddo
c       k=k+1
c       do l=1,llm
c         if(coor.ne.1)write(k,124) ((euva(i,j,l),i=1,iim),j=1,jjp1)
c       enddo
c       k=k+1
c       do l=1,llm
c         if(coor.ne.1)write(k,124) ((cona(i,j,l),i=1,iim),j=1,jjp1)
c       enddo
c       k=k+1
c       do l=1,llm
c         if(coor.ne.1)write(k,124) ((nira(i,j,l),i=1,iim),j=1,jjp1)
c       enddo
      endif
124   format(12(e10.4,1x))


       ENDDO		!LOOP in TIME
! =================================================================

 900   continue


c  ECRITURE FINALE si besoin
c  *************************
      ierr= NF_CLOSE(nid)
      ierr=NF_CLOSE(nout)
      close (33)
      close (122)
      close (1222)
       
      ENDDO		!LOOP in FILES
! =================================================================

9999  PRINT*,'Fin '
1000  format(a5,3x,i4,i3,x,a39)
7777  FORMAT ('latitude/longitude',4f7.1)

      END

c*******************************************************

      subroutine  xsig2p (ps,qsig,pref,npref,qp)
      IMPLICIT NONE
c=======================================================================
c
c   Francois Forget    Avril 1996
c
c Cette subroutine calcule les champs interpole en niveaux de pression
c different niveaux de pression pref
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "comconst.h"

c
c  ARGUMENTS
c  ---------

c Inputs:
      
      REAL qsig(iip1,jjp1,llm)
      REAL ps(iip1,jjp1)
      integer npref
      real pref(npref)

c outputs
      REAL qp(iip1,jjp1,npref)


c Variables du modele
c -------------------

      REAL lnp(llm)  , q(llm)
      REAL x,y
      INTEGER i,j,l ,n


      logical firstcall
      save firstcall
      data firstcall/.true./
 
#include "fxyprim.h"
c**********************************************************************
         if (firstcall) then
           write(*,*) 'kappa = ', kappa
           write(*,*) 'cpp = ', cpp
           firstcall = .false.
         end if

        DO j=1,jjp1
          DO i=1,iip1
            do l=1,llm
              lnp(l) = -log(aps(l)+bps(l)*ps(i,j))
              q(l) = qsig(i,j,l)
            end do
            do n =1,npref
              if ( (pref(n).le.ps(i,j)) .and.
     &          (pref(n).ge.(ps(i,j)*bps(llm)+aps(llm))) ) then
                    x = -log(pref(n))
                    call interpolf(x,y,lnp,q,llm)
                    qp(i,j,n) = y
              else
                    qp(i,j,n) = 1E+20 
              end if
            end do
          ENDDO
        ENDDO

      return
      end 

      subroutine  missing_value(nout,nvarid,valid_range,missing)
      IMPLICIT NONE
c=======================================================================
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------
      include "netcdf.inc"

      INTEGER nout,nvarid,ierr
      REAL missing
      REAL valid_range(2)
c
c  ARGUMENTS
c  ---------

c Inputs:
      
             ierr= NF_PUT_ATT_REAL(nout,nvarid,'valid_range',
     $  NF_FLOAT,2,valid_range)
        IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'anl_NC: valid_range attribution failed'
            WRITE(*,*) 'NF_NOERR', NF_NOERR
            CALL abort
        ENDIF

#ifdef NC_DOUBLE
             ierr= NF_PUT_ATT_DOUBLE(nout,nvarid,'missing_value',
     $  NF_DOUBLE,1,missing)
#else
              ierr= NF_PUT_ATT_REAL(nout,nvarid,'missing_value',
     $  NF_FLOAT,1,missing)

#endif

        IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'anl_NC: missing value attribution failed'
            WRITE(*,*) 'NF_NOERR', NF_NOERR
            CALL abort
        ENDIF
      return
      end
**********************************************************
      subroutine lect_var(nid,name,nbdim,ierr,itau,
     $ var)
      IMPLICIT NONE
c=======================================================================
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------
#include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "comconst.h"

c inputs
      character (len=*) :: name
      INTEGER nid,nbdim,nvarid,ierr,itau
      integer, dimension(nbdim) :: edges,corner
c output var
      integer, dimension(nbdim) :: var
c  ---------
        if (nbdim.eq.4) then
         corner=(/1,1,1,itau/)
         edges=(/iip1,jjp1,llm,1/)
        endif

        if (nbdim.eq.3) then
            corner=(/1,1,itau/)
            edges=(/iip1,jjp1,1/)
        endif


          ierr = NF_INQ_VARID (nid,adjustl(name), nvarid)
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid, nvarid,corner,edges, var)
#else
          ierr = NF_GET_VARA_REAL(nid, nvarid,corner,edges, var)
#endif
          IF (ierr.NE.NF_NOERR) THEN
            write(*,*) 'anl_NC: reading failed for ', name
            CALL abort
          ENDIF

      return
      end


      subroutine put_var(nout,name,title,units,nbdim,ierr,itau,
     $ var)
      IMPLICIT NONE
c=======================================================================
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------
#include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "comconst.h"

      character (len=*) :: title,units,name
      INTEGER nout,nbdim,nvarid,ierr,itau
      integer, dimension(nbdim) :: edges,corner,id,var
      REAL valid_range(2)
      DATA valid_range /0., 300.0/
c  ---------
        if (nbdim.eq.4) then
         ierr= NF_INQ_DIMID(nout,"longitude",id(1))
         ierr= NF_INQ_DIMID(nout,"latitude",id(2))
         ierr= NF_INQ_DIMID(nout,"altitude",id(3))
         ierr= NF_INQ_DIMID(nout,"Time",id(4))
         corner=(/1,1,1,itau/)
         edges=(/iip1,jjp1,llm,1/)
        endif

        if (nbdim.eq.3) then
             ierr= NF_INQ_DIMID(nout,"longitude",id(1))
             ierr= NF_INQ_DIMID(nout,"latitude",id(2))
             ierr= NF_INQ_DIMID(nout,"Time",id(3))
            corner=(/1,1,itau/)
            edges=(/iip1,jjp1,1/)
        endif


        ierr = NF_INQ_VARID (nout,adjustl(name), nvarid)
          if (ierr /= NF_NOERR) then
!  choix du nom des coordonnees
! Creation de la variable si elle n'existait pas
             write (*,*) "====================="
             write (*,*) "creation de  ",name
         call def_var(nout,adjustl(title),adjustl(name),adjustl(units),
     $   nbdim,id,nvarid,ierr)
          if (name.eq.'temp') then
             ierr = NF_REDEF (nout)
             call missing_value(nout,nvarid,valid_range,1E+20)
             ierr = NF_ENDDEF(nout)
          endif
          endif
#ifdef NC_DOUBLE
        ierr = NF_PUT_VARA_DOUBLE(nout,nvarid,corner,edges,var)
#else
        ierr = NF_PUT_VARA_REAL(nout,nvarid,corner,edges,var)
#endif
        IF (ierr.NE.NF_NOERR) THEN
            write(*,*) 'anl_NC: writing failed for ',name
            CALL abort
        ENDIF

                
      return
      end


      Subroutine interpolf(x,y,xd,yd,nd)

c******************************************************
c   SUBROUTINE   (interpol)
c interpolation, give y = f(x) with array xd,yd known, size nd

c  Version with CONSTANT values oustide limits
c**********************************************************


c Variable declaration
c --------------------
c  Arguments :
      real x,y
      real xd(*),yd(*)
      integer nd
c  internal
      integer i,j
      real y_undefined

c run
c ---
c      y_undefined=-999.999
      y_undefined=1.e20

      y=0.
      if ((x.le.xd(1)).and.(x.le.xd(nd))) then
        if (xd(1).lt.xd(nd)) y = y_undefined ! yd(1)
        if (xd(1).ge.xd(nd)) y = y_undefined ! yd(nd)
      else if ((x.ge.xd(1)).and.(x.ge.xd(nd))) then
        if (xd(1).lt.xd(nd)) y = y_undefined ! yd(nd)
        if (xd(1).ge.xd(nd)) y = y_undefined ! yd(1)
c       y = yd(nd)
      else
        do i=1,nd-1
         if ( ( (x.ge.xd(i)).and.(x.lt.xd(i+1)) )
     &     .or. ( (x.le.xd(i)).and.(x.gt.xd(i+1)) ) ) then
           y=yd(i)+(x-xd(i))*(yd(i+1)-yd(i))/(xd(i+1)-xd(i))
           goto 99
         end if
        end do
      end if

c     write (*,*)' x , y' , x,y
c     do i=1,nd
c       write (*,*)' i, xd , yd' , xd(i),yd(i)
c     end do
c     stop

 99   continue

      end
