      subroutine writediagfi(ngrid,nom,titre,unite,dim,px)

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

! Commons
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

! Arguments on input:
      integer ngrid
      character (len=*) :: nom,titre,unite
      integer dim
      real px(ngrid,llm)

! Local variables:

      real dx3(iip1,jjp1,llm) ! to store a 3D data set
      real dx2(iip1,jjp1)     ! to store a 2D (surface) data set
      real dx0

      real date

      REAL phis(ip1jmp1)

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
!     irythme = iconser  ! sortie au rythme des variables de controle
!     irythme = iphysiq  ! sortie a tous les pas physique
!     irythme = iecri*day_step ! sortie au rythme des fichiers histoires
!     irythme = periodav*day_step ! sortie au rythme des fichiers histmoy

!***************************************************************

! The following test is here to enforce that writediagfi is not used with the
! 1D version of the GCM
      if (ngrid.eq.1) return
      
c     nom=trim((nom))
c     unite=trim((unite))
c     titre=trim((titre))

! Initialisation of 'firstnom' and create/open the "diagfi.nc" NetCDF file
! ------------------------------------------------------------------------
! (Au tout premier appel de la subroutine durant le run.)

      fichnom="diagfi.nc"

      if (firstnom.eq.'1234567890') then ! .true. for the very first call
      !  to this subroutine; now set 'firstnom'
         firstnom = nom

         call gr_fi_dyn(1,ngrid,iip1,jjp1,phisfi,phis)
         ! Create the NetCDF file
         ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
         ! Define the 'Time' dimension
         ierr = nf_def_dim(nid,"Time",NF_UNLIMITED,idim)
         ! Define the 'Time' variable
#ifdef NC_DOUBLE
         ierr = NF_DEF_VAR (nid, "Time", NF_DOUBLE, 1, idim,varid)
#else
         ierr = NF_DEF_VAR (nid, "Time", NF_FLOAT, 1, idim,varid)
#endif
         ! Add a long_name attribute
         ierr = NF_PUT_ATT_TEXT (nid, varid, "long_name",
     .          4,"Time")
         ! Add a units attribute
         ierr = NF_PUT_ATT_TEXT(nid, varid,'units',29,
     .          "days since 0000-00-0 00:00:00")
         ! Switch out of NetCDF Define mode
         ierr = NF_ENDDEF(nid)

         ! write "header" of file (longitudes, latitudes, geopotential, ...)
         call iniwrite(nid,day_ini,phis)
         
         zitau = -1 ! initialize zitau
      else
         ! Open the NetCDF file
         ierr = NF_OPEN(fichnom,NF_WRITE,nid)
      endif ! if (firstnom.eq.'1234567890')

! Increment time index 'zitau' if it is the "fist call" (at given time level)
! to writediagfi
!------------------------------------------------------------------------
      if (nom.eq.firstnom) then
          zitau = zitau + iphysiq
      end if

!--------------------------------------------------------
! Write the variables to output file if it's time to do so
!--------------------------------------------------------

      if ( MOD(zitau+1,irythme) .eq.0.) then

! Compute/write/extend 'Time' coordinate (date given in days)
! (done every "first call" (at given time level) to writediagfi)
! Note: date is incremented as 1 step ahead of physics time
!       (like the 'histoire' outputs)
!--------------------------------------------------------

        if (nom.eq.firstnom) then
        ! We have identified a "first call" (at given date)
           ntime=ntime+1 ! increment # of stored time steps
           ! compute corresponding date (in days and fractions thereof)
           date= float (zitau +1)/float (day_step)
           ! Get NetCDF ID of 'Time' variable
           ierr= NF_INQ_VARID(nid,"Time",varid)
           ! Write (append) the new date to the 'Time' array
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
        end if ! of if (nom.eq.firstnom)

!Case of a 3D variable
!---------------------
        if (dim.eq.3) then

!         Passage variable physique -->  variable dynamique
!         recast (copy) variable from physics grid to dynamics grid
           DO l=1,llm
             DO i=1,iip1
                dx3(i,1,l)=px(1,l)
                dx3(i,jjp1,l)=px(ngrid,l)
             ENDDO
             DO j=2,jjm
                ig0= 1+(j-2)*iim
                DO i=1,iim
                   dx3(i,j,l)=px(ig0+i,l)
                ENDDO
                dx3(iip1,j,l)=dx3(1,j,l)
             ENDDO
           ENDDO

!         Ecriture du champs

!         write (*,*) 'In  writediagfi, on sauve:  ' , nom
!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
! name of the variable
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"altitude",id(3))
              ierr= NF_INQ_DIMID(nid,"Time",id(4))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           endif 

           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime

           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=llm
           edges(4)=1
#ifdef NC_DOUBLE
           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx3)
#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

!Case of a 2D variable
!---------------------

        else if (dim.eq.2) then

!         Passage variable physique -->  physique dynamique
!         recast (copy) variable from physics grid to dynamics grid

             DO i=1,iip1
                dx2(i,1)=px(1,1)
                dx2(i,jjp1)=px(ngrid,1)
             ENDDO
             DO j=2,jjm
                ig0= 1+(j-2)*iim
                DO i=1,iim
                   dx2(i,j)=px(ig0+i,1)
                ENDDO
                dx2(iip1,j)=dx2(1,j)
             ENDDO

!         write (*,*) 'In  writediagfi, on sauve:  ' , nom
!         write (*,*) 'In  writediagfi. Estimated date = ' ,date
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"Time",id(3))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom

              call def_var(nid,nom,titre,unite,3,id,varid,ierr)

           endif

           corner(1)=1
           corner(2)=1
           corner(3)=ntime
           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=1


#ifdef NC_DOUBLE
           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx2) 
#else         
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx2)
#endif     

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

!Case of a 0D variable (ie: a time-dependent scalar)
!---------------------------------------------------

        else if (dim.eq.0) then
           dx0 = px (1,1)

           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"Time",id(1))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",nom

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
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
c             call abort
           endif 

        endif ! of if (dim.eq.3) elseif(dim.eq.2)...

      endif ! of if ( MOD(zitau+1,irythme) .eq.0.)

      ierr= NF_CLOSE(nid)

      end
