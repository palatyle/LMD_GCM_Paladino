










      subroutine writediagspecIR(ngrid,nom,titre,unite,dimpx,px)

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
!      dimpx : dimension de px : 0, 2, ou 3 dimensions
!
!=================================================================
!
!      This is a modified version that accepts spectrally varying input
!      RW (2010)
!
!=================================================================
 
! Addition by RW (2010) to allow OLR to be saved in .nc format
      use radinc_h, only : L_NSPECTI
      use geometry_mod, only: cell_area
      use mod_phys_lmdz_para, only : is_mpi_root, is_master, gather
      use mod_grid_phy_lmdz, only : klon_glo, Grid1Dto2D_glo,
     &                              nbp_lon, nbp_lat
      use time_phylmdz_mod, only: ecritphy, iphysiq, day_step, day_ini
      use callkeys_mod, only: iradia

      implicit none

      include "netcdf.inc"

! Arguments on input:
      integer,intent(in) :: ngrid
      character (len=*),intent(in) :: nom,titre,unite
      integer,intent(in) :: dimpx
      real,intent(in) :: px(ngrid,L_NSPECTI)

! Local variables:

!      real dx3(iip1,jjp1,llm) ! to store a 3D data set
!      real dx2(iip1,jjp1)     ! to store a 2D (surface) data set
!      real dx0

      integer irythme
      integer ierr
      integer iq
      integer i,j,l,zmax , ig0

      integer,save :: zitau=0
      character(len=20),save :: firstnom='1234567890'
      real,save :: date
!$OMP THREADPRIVATE(firstnom,zitau,date)

! Ajouts
      integer, save :: ntime=0
!$OMP THREADPRIVATE(ntime)
      integer :: idim,varid
      integer :: nid
      character (len =50):: fichnom
      integer, dimension(4) :: id
      integer, dimension(4) :: edges,corner

      real area((nbp_lon+1),nbp_lat)
! added by RDW for OLR output
      real dx3(nbp_lon+1,nbp_lat,L_NSPECTI) ! to store the data set
      real dx3_1d(1,L_NSPECTI) ! to store the data with 1D model

      real areafi_glo(ngrid) ! mesh area on global physics grid

!***************************************************************
!Sortie des variables au rythme voulu

      irythme = ecritphy*iradia ! sortie au rythme de ecritphy*iradia
!EM+JL if the spetra need to be output more frequently, need to define a ecritSpec...
!     irythme = iphysiq  ! sortie a tous les pas physique


!***************************************************************

! Initialisation of 'firstnom' and create/open the "diagfi.nc" NetCDF file
! ------------------------------------------------------------------------
! (Au tout premier appel de la subroutine durant le run.)

      fichnom="diagspecIR.nc"

      if (firstnom.eq.'1234567890') then ! .true. for the very first call
      !  to this subroutine; now set 'firstnom'
         firstnom = nom
         ! just to be sure, check that firstnom is large enough to hold nom
         if (len_trim(firstnom).lt.len_trim(nom)) then
           write(*,*) "writediagspecIR: Error !!!"
           write(*,*) "   firstnom string not long enough!!"
           write(*,*) "   increase its size to at least ",len_trim(nom)
           stop
         endif

         areafi_glo(:)=cell_area(:)
         ! Create the NetCDF file
         if (is_master) then
         ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
         ! Define the 'Time' dimension
         ierr = nf_def_dim(nid,"Time",NF_UNLIMITED,idim)
         ! Define the 'Time' variable
         ierr = NF_DEF_VAR (nid, "Time", NF_DOUBLE, 1, idim,varid)
         ! Add a long_name attribute
         ierr = NF_PUT_ATT_TEXT (nid, varid, "long_name",
     .          4,"Time")
         ! Add a units attribute
         ierr = NF_PUT_ATT_TEXT(nid, varid,'units',29,
     .          "days since 0000-00-0 00:00:00")
         ! Switch out of NetCDF Define mode
         ierr = NF_ENDDEF(nid)

         ! Build area()
         IF (klon_glo>1) THEN
          do i=1,nbp_lon+1 ! poles
           ! divide at the poles by nbp_lon
           area(i,1)=areafi_glo(1)/nbp_lon
           area(i,nbp_lat)=areafi_glo(klon_glo)/nbp_lon
          enddo
          do j=2,nbp_lat-1
           ig0= 1+(j-2)*nbp_lon
           do i=1,nbp_lon
              area(i,j)=areafi_glo(ig0+i)
           enddo
           ! handle redundant point in longitude
           area(nbp_lon+1,j)=area(1,j)
          enddo
         ENDIF

         ! write "header" of file (longitudes, latitudes, area, ...)
         IF (klon_glo>1) THEN ! general 3D case
           call iniwrite_specIR(nid,day_ini,area,nbp_lon+1,nbp_lat)
         ELSE
           call iniwrite_specIR(nid,day_ini,areafi_glo(1),1,1)
         ENDIF
         endif ! of if (is_master)

         zitau = -1 ! initialize zitau
      else
         if (is_master) then
           ! Open the NetCDF file
           ierr = NF_OPEN(fichnom,NF_WRITE,nid)
         endif
      endif ! if (firstnom.eq.'1234567890')

! Increment time index 'zitau' if it is the "firstcall" (at given time level)
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

           if (is_master) then
             ! Get NetCDF ID of 'Time' variable
             ierr= NF_INQ_VARID(nid,"Time",varid)

             ! Write (append) the new date to the 'Time' array
             ierr= NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
             if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagspec_nc"
              write(*,*) "***** with time"
              write(*,*) 'ierr=', ierr   
c             call abort
             endif

             write(6,*)'WRITEDIAGSPECIR: date= ', date
           endif ! of if (is_master)
        end if ! of if (nom.eq.firstnom)


 
!Case of a 3D variable
!---------------------
        if (dimpx.eq.3) then

!         A. Recast (copy) variable from physics grid to dynamics grid
          IF (klon_glo>1) THEN ! General case
           DO l=1,L_NSPECTI
             DO i=1,nbp_lon+1
                dx3(i,1,l)=px(1,l)
                dx3(i,nbp_lat,l)=px(ngrid,l)
             ENDDO
             DO j=2,nbp_lat-1
                ig0= 1+(j-2)*nbp_lon
                DO i=1,nbp_lon
                   dx3(i,j,l)=px(ig0+i,l)
                ENDDO
                dx3(nbp_lon+1,j,l)=dx3(1,j,l)
             ENDDO
           ENDDO
          ELSE ! 1D model case
            dx3_1d(1,1:L_NSPECTI)=px(1,1:L_NSPECTI)
          ENDIF

!         B. Write (append) the variable to the NetCDF file
          if (is_master) then

! name of the variable
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"longitude",id(1))
              ierr= NF_INQ_DIMID(nid,"latitude",id(2))
              ierr= NF_INQ_DIMID(nid,"IR_Wavenumber",id(3))
              ierr= NF_INQ_DIMID(nid,"Time",id(4))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGSPECIR: creating variable ",nom
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           endif

           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime

           IF (klon_glo==1) THEN
             edges(1)=1
           ELSE
             edges(1)=nbp_lon+1
           ENDIF
           edges(2)=nbp_lat
           edges(3)=L_NSPECTI
           edges(4)=1
           IF (klon_glo>1) THEN ! General case
             ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3)
           ELSE
             ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3_1d)
           ENDIF

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagspec"
              write(*,*) "***** with ",nom
              write(*,*) 'ierr=', ierr
             call abort
           endif 

          endif ! of if (is_master)

        endif ! of if (dimpx.eq.3)

      endif ! of if ( MOD(zitau+1,irythme) .eq.0.)

      ! Close the NetCDF file
      if (is_master) then
        ierr= NF_CLOSE(nid)
      endif

      end
