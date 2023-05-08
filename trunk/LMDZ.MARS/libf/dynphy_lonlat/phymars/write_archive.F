c=======================================================================
      subroutine write_archive(nid,ntime,nom,titre,unite,dim,px)
c=======================================================================
c
c
c   Date:    01/1997
c   ----
c
c   Objet:   Ecriture de champs sur grille scalaire (iip1*jjp1)
c   -----    dans un fichier DRS nomme "start_archive"
c
c    Il faut au prealable avoir cree un entete avec un "call ini_archive".
c    Ces variables peuvent etre 3d (ex: temperature), 2d (ex: temperature
c    de surface), ou 0d (pour un scalaire qui ne depend que du temps)
c    (ex: la longitude solaire)
c
c
c   Arguments: 
c   ----------
c
c     Inputs:
c     ------
c
c		  nid      Unite logique du fichier "start_archive"
c         nom      nom du champ a ecrire dans le fichier "start_archive"
c         titre    titre de la variable dans le fichier DRS "start_archive"
c         unite    unite de la variable ....
c         dim      dimension de la variable a ecrire
c         px       tableau contenant la variable a ecrire
c
c
c=======================================================================

      use comsoil_h, only: nsoilmx
      implicit none

#include "dimensions.h"
#include "paramet.h"
!#include "control.h"
#include "comgeom.h"
#include "netcdf.inc"

c-----------------------------------------------------------------------
c	Declarations   
c-----------------------------------------------------------------------

c Arguments:

      INTEGER nid
      integer ntime ! time index
      integer dim 
      REAL px(iip1,jjp1,llm) 

      CHARACTER*(*) nom, titre, unite

      integer ierr


c local
      integer, dimension(4) :: edges,corner,id
      integer :: varid,i,j,l
c-----------------------------------------------------------------------
c      Ecriture du champs dans le fichier            (3 cas)      
c-----------------------------------------------------------------------

! For an atmospheric 3D Variable
!--------------------------------
        if (dim.eq.3) then

!         Ecriture du champs

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
              write (*,*) "creation de ",nom
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           endif

! mars s'arranger pour qu'il n'y ai plus besoin de ca

c          do l=1,llm
c             do j=1,jjp1
c                do i=1,iip1
c                   pxbis(i,j,l)=px(i,j,llm-l+1)
c                enddo
c             enddo
c          enddo
           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime

           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=llm
           edges(4)=1
#ifdef NC_DOUBLE
           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,px)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,px)
#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in write_archive"
              write(*,*) "***** with ",nom," ",nf_STRERROR(ierr)
              call abort
           endif


! For a subterranean 3D Variable
!-------------------------------

        else if (dim.eq.-3) then
	! get variables' ID, if it exists
	ierr=NF_INQ_VARID(nid,nom,varid)
	
	 if (ierr.ne.NF_NOERR) then ! variable not defined yet
	  ! build related coordinates
	  ierr=NF_INQ_DIMID(nid,"longitude",id(1))
	  ierr=NF_INQ_DIMID(nid,"latitude",id(2))
	  ierr=NF_INQ_DIMID(nid,"subsurface_layers",id(3))
	  if (ierr.ne.NF_NOERR) then
	   write(*,*)"write_archive: dimension <subsurface_layers>",
     &               " is missing !!!"
	   call abort
	  endif
          ierr=NF_INQ_DIMID(nid,"Time",id(4))
	  
	  ! define the variable
	  write(*,*)"====================="
	  write(*,*)"defining ",nom
	  call def_var(nid,nom,titre,unite,4,id,varid,ierr)
	  
	 endif

        ! build cedges and corners
        corner(1)=1
        corner(2)=1
        corner(3)=1
        corner(4)=ntime

        edges(1)=iip1
        edges(2)=jjp1
        edges(3)=nsoilmx
        edges(4)=1
#ifdef NC_DOUBLE
           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,px)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,px)
#endif


! For a surface 2D Variable
!--------------------------

        else if (dim.eq.2) then

!         Ecriture du champs

           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
!  choix du nom des coordonnees
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"Time",id(3))

! Creation de la variable si elle n'existait pas

              write (*,*) "====================="
              write (*,*) "creation de ",nom

              call def_var(nid,nom,titre,unite,3,id,varid,ierr)

           endif

           corner(1)=1
           corner(2)=1
           corner(3)=ntime
           edges(1)=iip1
           edges(2)=jjp1
           edges(3)=1


#ifdef NC_DOUBLE
           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,px)
#else         
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,px)
#endif     

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in write_archive"
              write(*,*) "***** with ",nom,nf_STRERROR(ierr)
              call abort
           endif


!Cas Variable 0D (scalaire dependant du temps)
!---------------------------------------------

        else if (dim.eq.0) then

!         Ecriture du champs

           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
!  choix du nom des coordonnees
              ierr= NF_INQ_DIMID(nid,"Time",id(1))

! Creation de la variable si elle n'existait pas

              write (*,*) "====================="
              write (*,*) "creation de ",nom

              call def_var(nid,nom,titre,unite,1,id,varid,ierr)

           endif

           corner(1)=ntime
           edges(1)=1

#ifdef NC_DOUBLE
           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,px)
#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,px)
#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in write_archive"
              write(*,*) "***** with ",nom,nf_STRERROR(ierr)
              call abort
           endif

        else
	  write(*,*) "write_archive: dim=",dim," ?!?"
	  call abort
        endif ! of if (dim.eq.3) else if (dim.eq.-3) ....

      return
      end

