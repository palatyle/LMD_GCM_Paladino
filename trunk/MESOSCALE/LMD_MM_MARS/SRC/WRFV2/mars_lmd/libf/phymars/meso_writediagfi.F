      subroutine meso_writediagfi(ngrid,nom,titre,unite,dim,px)

c=======================================================================
c
c	CAREFUL: THIS IS A VERSION TO BE USED WITH WRF !!!
c
c	... CHECK THE ****WRF lines
c
c=======================================================================
!  Ecriture de variables diagnostiques au choix dans la physique 
!  dans un fichier NetCDF nomme  'diagfi'. Ces variables peuvent etre
!  3d (ex : temperature), 2d (ex : temperature de surface), ou
!  0d (pour un scalaire qui ne depend que du temps : ex : la longitude
!  solaire)
!  Dans la version 2000, la periode d'ecriture est celle de 
!  "ecritphy " regle dans le fichier de controle de run :  run.def
!
!    writediagfi peut etre appele de n'importe quelle subroutine
!    de la physique, plusieurs fois. L'initialisation et la creation du
!    fichier se fait au tout premier appel.
!
! WARNING : les variables dynamique (u,v,t,q,ps)
!  sauvees par writediagfi avec une
! date donnee sont legerement differentes que dans le fichier histoire car 
! on ne leur a pas encore ajoute de la dissipation et de la physique !!!
! IL est  RECOMMANDE d'ajouter les tendance physique a ces variables
! avant l'ecriture dans diagfi (cf. physiq.F)
!  
!
!  parametres (input) :
!  ----------
!      ngrid : nombres de point ou est calcule la physique
!                (ngrid = 2+(jjm-1)*iim - 1/jjm)
!                 (= nlon ou klon dans la physique terrestre)
!      
!      unit : unite logique du fichier de sortie (toujours la meme)
!      nom  : nom de la variable a sortir (chaine de caracteres)
!      titre: titre de la variable (chaine de caracteres)
!      unite : unite de la variable (chaine de caracteres)
!      px : variable a sortir (real 0, 2, ou 3d)
!      dim : dimension de px : 0, 2, ou 3 dimensions
!
!=================================================================
 
      implicit none

! Arguments variable

! Local variables
#include "dimensions.h"
#include "dimphys.h"
#include "paramet.h"
#include "control.h"
#include "comvert.h"
#include "comgeom.h"
#include "description.h"
#include "netcdf.inc"
#include "temps.h"
#include "surfdat.h"

      integer ngrid
	
!****WRF
        INTEGER idim_alt,idim_lon,idim_lat
        REAL px(ngrid,nlayermx)
        REAL dx3(wiim,wjjm,nlayermx),dx2(wiim,wjjm)
!****WRF
!!      REAL px(ngrid,llm)
!!      REAL dx3(iip1,jjp1,llm),dx2(iip1,jjp1)

      REAL dx0

      real date
      character (len=*) :: nom,titre,unite

      REAL phis(ip1jmp1)

      integer dim
      integer irythme
      integer ierr
      integer iq
      integer i,j,l,zmax , ig0

      integer zitau
      character firstnom*10
      SAVE firstnom
      SAVE zitau
      SAVE date
      data firstnom /'1234567890'/
      data zitau /0/

! Ajouts
      integer, save :: ntime=0
      integer :: idim,varid
      integer :: nid
      character (len =50):: fichnom
      integer, dimension(4) :: id
      integer, dimension(4) :: edges,corner
      

!***************************************************************
!Sortie des variables au rythme voulu

      irythme = int(ecritphy) ! sortie au rythme de ecritphy
        !print*,irythme 
!     irythme = iconser  ! sortie au rythme des variables de controle
!     irythme = iphysiq  ! sortie a tous les pas physique
!     irythme = iecri*day_step ! sortie au rythme des fichiers histoires
!     irythme = periodav*day_step ! sortie au rythme des fichiers histmoy

!***************************************************************
 
c     nom=trim((nom))
c     unite=trim((unite))
c     titre=trim((titre))

c
c lat / lon / alt dimensions from WRF (in dimphys.h)
c
        idim_lat=wjjm
c        write (*,*) "lat",idim_lat
        idim_lon=wiim
c        write (*,*) "lon",idim_lon
        idim_alt=nlayermx
c        write (*,*) "alt",idim_alt     



! initialisation de 'firstnom' / ouverture du fichier NetCDF
! ----------------------------------------------------------
! (Au tout premier appel de la subroutine durant le run.)

      if(ngrid.eq.1) return ! don't use writediagfi with 1D version

      fichnom="diagfi.nc"
      if (firstnom.eq.'1234567890') then
         firstnom = nom

!****WRF
!         call gr_fi_dyn(1,ngrid,iip1,jjp1,phisfi,phis)
!****WRF

! assign Logical Unit: ouverture du fichier NetCDF
         ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
! Creation de la dimension (temps)
         ierr = nf_def_dim(nid,"Time",NF_UNLIMITED,idim)
        ierr = NF_DEF_DIM (nid, "latitude", wjjm, idim_lat)
        ierr = NF_DEF_DIM (nid, "longitude", wiim, idim_lon)
        ierr = NF_DEF_DIM (nid, "altitude", nlayermx,idim_alt)


#ifdef NC_DOUBLE
         ierr = NF_DEF_VAR (nid, "Time", NF_DOUBLE, 1, idim,varid)
#else
         ierr = NF_DEF_VAR (nid, "Time", NF_FLOAT, 1, idim,varid)
#endif
         ierr = NF_PUT_ATT_TEXT (nid, varid, "long_name",
     .          4,"Time")
         ierr = NF_PUT_ATT_TEXT(nid, varid,'units',29,
     .          "days since 0000-00-0 00:00:00")
         ierr = NF_ENDDEF(nid)
! ecriture de l'entete du fichier (longitudes, latitudes ... relief)

!!****caca caca
!                phis(:)=0.   
!         call iniwrite(nid,day_ini , phis)
!!
!        ierr = NF_DEF_DIM (nid, "latitude", wjjm, idim_lat)
!        ierr = NF_DEF_DIM (nid, "longitude", wiim, idim_lon)
!        ierr = NF_DEF_DIM (nid, "altitude", nlayermx, idim_alt)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid, "latitude", NF_DOUBLE, 1, idim_lat,varid)
#else
      ierr =NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim_lat,varid)
#endif
      ierr =NF_PUT_ATT_TEXT(nid,varid,'units',13,"degrees_north")
      ierr = NF_PUT_ATT_TEXT (nid,varid,"long_name", 14,
     .      "North latitude")
      ierr = NF_ENDDEF(nid)
!#ifdef NC_DOUBLE
!      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu/pi*180)
!#else
!      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu/pi*180)
!#endif
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid,"longitude", NF_DOUBLE, 1, idim_lon,varid)
#else
      ierr = NF_DEF_VAR(nid,"longitude", NF_FLOAT, 1,idim_lon,varid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,varid,"long_name",14,
     .      "East longitude")
      ierr =NF_PUT_ATT_TEXT(nid,varid,'units',12,"degrees_east")
      ierr = NF_ENDDEF(nid)
!#ifdef NC_DOUBLE
!      ierr = NF_PUT_VAR_DOUBLE(nid,nvarid,rlonv/pi*180)
!#else
!      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv/pi*180)
!#endif
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "altitude", NF_DOUBLE, 1,
     .       idim_alt,varid)
#else
      ierr = NF_DEF_VAR (nid, "altitude", NF_FLOAT, 1,
     .       idim_alt,varid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,varid,"long_name",10,"pseudo-alt")
      ierr = NF_PUT_ATT_TEXT (nid,varid,'units',2,"km")
      ierr = NF_PUT_ATT_TEXT (nid,varid,'positive',2,"up")

      ierr = NF_ENDDEF(nid)
!#ifdef NC_DOUBLE
!      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,pseudoalt)
!#else
!      ierr = NF_PUT_VAR_REAL (nid,nvarid,pseudoalt)
!#endif



         zitau = -1
      else
         ierr = NF_OPEN(fichnom,NF_WRITE,nid)
      endif

!incrementation de temps a chaque premier appel de WRITEDIAGFI ds la physique
!------------------------------------------------------------------------
      if (nom.eq.firstnom) then
          !! -- elapsed seconds --   
          zitau = zitau + iphysiq
          !! iphysiq is ptimestep          
      end if

!--------------------------------------------------------
!Sortie des variables au rythme voulu
!--------------------------------------------------------

      if ( MOD(zitau+1,irythme) .eq.0.) then

!Calcul/ecriture/extension de la coordonnee temps (date en jours)
!--------------------------------------------------------
!(effectuee a chaque premier appel de writediagfi ds la physique)
!On date les sorties comme histoire, c.a.d un pas de temps plus loin !

!Detection du changement de temps
        if (nom.eq.firstnom) then
           ntime=ntime+1
!           date= float (zitau +1)/float (day_step)
           date=float(zitau+1)/3700.
!
           ierr= NF_INQ_VARID(nid,"Time",varid)
#ifdef NC_DOUBLE
           ierr= NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with time"
              write(*,*) 'ierr=', ierr   
c             call abort
           endif

           write(6,*)'WRITEDIAGFI: date= ', date
        end if

!Cas Variable 3D
!---------------
        if (dim.eq.3) then

!         Passage variable physique -->  variable dynamique

!****WRF

           DO l=1,nlayermx
             DO j=1,wjjm
                ig0=(j-1)*wiim 
                DO i=1,wiim
                   dx3(i,j,l)=px(ig0+i,l)
                ENDDO
             ENDDO
           ENDDO


!!           DO l=1,llm
!!             DO i=1,iip1
!!                dx3(i,1,l)=px(1,l)
!!                dx3(i,jjp1,l)=px(ngrid,l)
!!             ENDDO
!!             DO j=2,jjm
!!                ig0= 1+(j-2)*iim
!!                DO i=1,iim
!!                   dx3(i,j,l)=px(ig0+i,l)
!!                ENDDO
!!                dx3(iip1,j,l)=dx3(1,j,l)
!!             ENDDO
!!           ENDDO

!****WRF

!         Ecriture du champs

!         write (*,*) 'In  writediagfi, on sauve:  ' , nom
!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
! nom de la variable
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! choix du nom des coordonnees
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"altitude",id(3))
              ierr= NF_INQ_DIMID(nid,"Time",id(4))

! Creation de la variable si elle n'existait pas

              write (*,*) "====================="
              write (*,*) "DIAGFI: creation de ",nom
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           endif 

           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime
!****WRF
!!           edges(1)=iip1
!!           edges(2)=jjp1
           edges(1)=wiim
           edges(2)=wjjm
           edges(3)=nlayermx
!****WRF

           edges(4)=1
#ifdef NC_DOUBLE
           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx3)
#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

!Cas Variable 2D
!---------------

        else if (dim.eq.2) then

!         Passage variable physique -->  physique dynamique

!****WRF

!           DO l=1,nlayermx
             DO j=1,wjjm
                ig0=(j-1)*wiim 
                DO i=1,wiim
                   dx2(i,j)=px(ig0+i,1)
                ENDDO
             ENDDO
!           ENDDO

!!             DO i=1,iip1
!!                dx2(i,1)=px(1,1)
!!                dx2(i,jjp1)=px(ngrid,1)
!!             ENDDO
!!             DO j=2,jjm
!!                ig0= 1+(j-2)*iim
!!                DO i=1,iim
!!                   dx2(i,j)=px(ig0+i,1)
!!                ENDDO
!!                dx2(iip1,j)=dx2(1,j)
!!             ENDDO

!****WRF

!         Ecriture du champs

!         write (*,*) 'In  writediagfi, on sauve:  ' , nom
!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
!  choix du nom des coordonnees
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"Time",id(3))

! Creation de la variable si elle n'existait pas

              write (*,*) "====================="
              write (*,*) "DIAGFI: creation de ",nom

              call def_var(nid,nom,titre,unite,3,id,varid,ierr)

           endif

           corner(1)=1
           corner(2)=1
           corner(3)=ntime
!****WRF
!!           edges(1)=iip1
!!           edges(2)=jjp1
           edges(1)=wiim
           edges(2)=wjjm
!****WRF
           edges(3)=1


#ifdef NC_DOUBLE
           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx2) 
#else         
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx2)
#endif     

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

!Cas Variable 0D (scalaire dependant du temps)
!---------------------------------------------

        else if (dim.eq.0) then
           dx0 = px (1,1)

!         Ecriture du champs

!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
!  choix du nom des coordonnees
              ierr= NF_INQ_DIMID(nid,"Time",id(1))

! Creation de la variable si elle n'existait pas

              write (*,*) "====================="
              write (*,*) "DIAGFI: creation de ",nom

              call def_var(nid,nom,titre,unite,1,id,varid,ierr)

           endif

           corner(1)=ntime
           edges(1)=1

#ifdef NC_DOUBLE
           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx0)  
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx0)
#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

        endif

      endif

      ierr= NF_CLOSE(nid)

      end
