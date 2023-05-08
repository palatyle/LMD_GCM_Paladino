      subroutine writediagfi(ngrid,nom,titre,unite,dim,px)

!  Ecriture de variables diagnostiques au choix dans la physique 
!  dans un fichier NetCDF nomme  'diagfi'. Ces variables peuvent etre
!  3d (ex : temperature), 2d (ex : temperature de surface), ou
!  0d (pour un scalaire qui ne depend que du temps : ex : la longitude
!  solaire)
!  (ou encore 1d, dans le cas de testphys1d, pour sortir une colonne)
!  La periode d'ecriture est donnee par 
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
! Modifs: Aug.2010 Ehouarn: enforce outputs to be real*4
!         Oct 2011 Francois: enable having a 'diagfi.def' file to select
!                            at runtime, which variables to put in file
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
!      px : variable a sortir (real 0, 1, 2, ou 3d)
!      dim : dimension de px : 0, 1, 2, ou 3 dimensions
!
!=================================================================
      use surfdat_h, only: phisfi
      use geometry_mod, only: cell_area
      use time_phylmdz_mod, only: ecritphy, day_step, iphysiq, day_ini
      USE mod_phys_lmdz_para, only : is_parallel, is_mpi_root,
     &                               is_master, gather
      USE mod_grid_phy_lmdz, only : klon_glo, Grid1Dto2D_glo,
     &                              nbp_lon, nbp_lat, nbp_lev
      implicit none

! Commons
      include "netcdf.inc"

! Arguments on input:
      integer,intent(in) :: ngrid
      character (len=*),intent(in) :: nom,titre,unite
      integer,intent(in) :: dim
      real,intent(in) :: px(ngrid,nbp_lev)

! Local variables:

      real*4 dx3(nbp_lon+1,nbp_lat,nbp_lev) ! to store a 3D data set
      real*4 dx2(nbp_lon+1,nbp_lat)     ! to store a 2D (surface) data set
      real*4 dx1(nbp_lev)           ! to store a 1D (column) data set
      real*4 dx0
      real*4 dx3_1d(1,nbp_lev) ! to store a profile with 1D model
      real*4 dx2_1d ! to store a surface value with 1D model

      real*4,save :: date
!$OMP THREADPRIVATE(date)

      REAL phis((nbp_lon+1),nbp_lat)
      REAL area((nbp_lon+1),nbp_lat)

      integer irythme
      integer ierr,ierr2
      integer i,j,l, ig0

      integer,save :: zitau=0
      character(len=20),save :: firstnom='1234567890'
!$OMP THREADPRIVATE(zitau,firstnom)

! Ajouts
      integer, save :: ntime=0
!$OMP THREADPRIVATE(ntime)
      integer :: idim,varid
      integer :: nid
      character(len=*),parameter :: fichnom="diagfi.nc"
      integer, dimension(4) :: id
      integer, dimension(4) :: edges,corner

! Added to use diagfi.def to select output variable
      logical,save :: diagfi_def
      logical :: getout
      integer,save :: n_nom_def
      integer :: n
      integer,parameter :: n_nom_def_max=199
      character(len=120),save :: nom_def(n_nom_def_max)
      logical,save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall) 	!diagfi_def,n_nom_def,nom_def read in diagfi.def
      
#ifdef CPP_PARA
! Added to work in parallel mode
      real dx3_glop(klon_glo,nbp_lev)
      real dx3_glo(nbp_lon,nbp_lat,nbp_lev) ! to store a global 3D data set
      real dx2_glop(klon_glo)
      real dx2_glo(nbp_lon,nbp_lat)     ! to store a global 2D (surface) data set
      real px2(ngrid)
!      real dx1_glo(nbp_lev)          ! to store a 1D (column) data set
!      real dx0_glo
      real phisfi_glo(klon_glo) ! surface geopotential on global physics grid
      real areafi_glo(klon_glo) ! mesh area on global physics grid
#else
      real phisfi_glo(ngrid) ! surface geopotential on global physics grid
      real areafi_glo(ngrid) ! mesh area on global physics grid
#endif

!***************************************************************
!Sortie des variables au rythme voulu

      irythme = int(ecritphy) ! output rate set by ecritphy

!***************************************************************

! At very first call, check if there is a "diagfi.def" to use and read it
! -----------------------------------------------------------------------
      IF (firstcall) THEN
         firstcall=.false.

!$OMP MASTER
  !      Open diagfi.def definition file if there is one:
         open(99,file="diagfi.def",status='old',form='formatted',
     s   iostat=ierr2)

         if(ierr2.eq.0) then
            diagfi_def=.true.
            write(*,*) "******************"
            write(*,*) "Reading diagfi.def"
            write(*,*) "******************"
            do n=1,n_nom_def_max
              read(99,fmt='(a)',end=88) nom_def(n)
              write(*,*) 'Output in diagfi: ', trim(nom_def(n))
            end do 
 88         continue
            if (n.ge.n_nom_def_max) then
               write(*,*)"n_nom_def_max too small in writediagfi.F:",n
               stop
            end if 
            n_nom_def=n-1
            close(99)
         else
            diagfi_def=.false.
         endif
!$OMP END MASTER
!$OMP BARRIER
      END IF ! of IF (firstcall)

! Get out of write_diagfi if there is diagfi.def AND variable not listed
!  ---------------------------------------------------------------------
      if (diagfi_def) then
          getout=.true.
          do n=1,n_nom_def
             if(trim(nom_def(n)).eq.nom) getout=.false.
          end do
          if (getout) return
      end if

! Initialisation of 'firstnom' and create/open the "diagfi.nc" NetCDF file
! ------------------------------------------------------------------------
! (at very first call to the subroutine, in accordance with diagfi.def)

      if (firstnom.eq.'1234567890') then ! .true. for the very first valid
      !   call to this subroutine; now set 'firstnom'
         firstnom = nom
         ! just to be sure, check that firstnom is large enough to hold nom
         if (len_trim(firstnom).lt.len_trim(nom)) then
           write(*,*) "writediagfi: Error !!!"
           write(*,*) "   firstnom string not long enough!!"
           write(*,*) "   increase its size to at least ",len_trim(nom)
           stop
         endif
         
         zitau = -1 ! initialize zitau

#ifdef CPP_PARA
          ! Gather phisfi() geopotential on physics grid
          call Gather(phisfi,phisfi_glo)
          ! Gather cell_area() mesh area on physics grid
          call Gather(cell_area,areafi_glo)
#else
         phisfi_glo(:)=phisfi(:)
         areafi_glo(:)=cell_area(:)
#endif

         !! parallel: we cannot use the usual writediagfi method
!!         call iophys_ini
         if (is_master) then
         ! only the master is required to do this

         ! Create the NetCDF file
         ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
         ! Define the 'Time' dimension
         ierr = nf_def_dim(nid,"Time",NF_UNLIMITED,idim)
         ! Define the 'Time' variable
!#ifdef NC_DOUBLE
!         ierr = NF_DEF_VAR (nid, "Time", NF_DOUBLE, 1, idim,varid)
!#else
         ierr = NF_DEF_VAR (nid, "Time", NF_FLOAT, 1, idim,varid)
!#endif
         ! Add a long_name attribute
         ierr = NF_PUT_ATT_TEXT (nid, varid, "long_name",
     .          4,"Time")
         ! Add a units attribute
         ierr = NF_PUT_ATT_TEXT(nid, varid,'units',29,
     .          "days since 0000-00-0 00:00:00")
         ! Switch out of NetCDF Define mode
         ierr = NF_ENDDEF(nid)

         ! Build phis() and area()
         IF (klon_glo>1) THEN
          do i=1,nbp_lon+1 ! poles
           phis(i,1)=phisfi_glo(1)
           phis(i,nbp_lat)=phisfi_glo(klon_glo)
           ! for area, divide at the poles by nbp_lon
           area(i,1)=areafi_glo(1)/nbp_lon
           area(i,nbp_lat)=areafi_glo(klon_glo)/nbp_lon
          enddo
          do j=2,nbp_lat-1
           ig0= 1+(j-2)*nbp_lon
           do i=1,nbp_lon
              phis(i,j)=phisfi_glo(ig0+i)
              area(i,j)=areafi_glo(ig0+i)
           enddo
           ! handle redundant point in longitude
           phis(nbp_lon+1,j)=phis(1,j)
           area(nbp_lon+1,j)=area(1,j)
          enddo
         ENDIF
         
         ! write "header" of file (longitudes, latitudes, geopotential, ...)
         IF (klon_glo>1) THEN ! general 3D case
           call iniwrite(nid,day_ini,phis,area,nbp_lon+1,nbp_lat)
         ELSE ! 1D model
           call iniwrite(nid,day_ini,phisfi_glo(1),areafi_glo(1),1,1)
         ENDIF

         endif ! of if (is_master)

      else

         if (is_master) then
           ! only the master is required to do this

           ! Open the NetCDF file
           ierr = NF_OPEN(fichnom,NF_WRITE,nid)
         endif ! of if (is_master)

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
!--------------------------------------------------------

        if (is_master) then
           ! only the master is required to do this
        if (nom.eq.firstnom) then
        ! We have identified a "first call" (at given date)
           ntime=ntime+1 ! increment # of stored time steps
           ! compute corresponding date (in days and fractions thereof)
           date=(zitau +1.)/day_step
           ! Get NetCDF ID of 'Time' variable
           ierr= NF_INQ_VARID(nid,"Time",varid)
           ! Write (append) the new date to the 'Time' array
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
!#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with time"
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)  
c             call abort
           endif

           write(6,*)'WRITEDIAGFI: date= ', date
        end if ! of if (nom.eq.firstnom)

        endif ! of if (is_master)

!Case of a 3D variable
!---------------------
        if (dim.eq.3) then

#ifdef CPP_PARA
          ! Gather field on a "global" (without redundant longitude) array
          call Gather(px,dx3_glop)
!$OMP MASTER
          if (is_mpi_root) then
            call Grid1Dto2D_glo(dx3_glop,dx3_glo)
            ! copy dx3_glo() to dx3(:) and add redundant longitude
            dx3(1:nbp_lon,:,:)=dx3_glo(1:nbp_lon,:,:)
            dx3(nbp_lon+1,:,:)=dx3(1,:,:)
          endif
!$OMP END MASTER
!$OMP BARRIER
#else
!         Passage variable physique -->  variable dynamique
!         recast (copy) variable from physics grid to dynamics grid
          IF (klon_glo>1) THEN ! General case
           DO l=1,nbp_lev
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
           dx3_1d(1,1:nbp_lev)=px(1,1:nbp_lev)
          ENDIF
#endif
!         Ecriture du champs

          if (is_master) then
           ! only the master writes to output
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
              write (*,*) "DIAGFI: creating variable ",trim(nom)
              call def_var(nid,nom,titre,unite,4,id,varid,ierr)

           else
             if (ntime==0) then
              write(*,*) "DIAGFI Error: failed creating variable ",
     &                   trim(nom)
              write(*,*) "it seems it already exists!"
              stop
             endif
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
           edges(3)=nbp_lev
           edges(4)=1
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx3)
!#else
!           write(*,*)"test:  nid=",nid," varid=",varid
!           write(*,*)"       corner()=",corner
!           write(*,*)"       edges()=",edges
!           write(*,*)"       dx3()=",dx3
           IF (klon_glo>1) THEN ! General case
             ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx3)
           ELSE
             ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx3_1d)
           ENDIF
!#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagfi"
              write(*,*) "***** with dx3: ",trim(nom)
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
              stop
           endif 

          endif !of if (is_master)

!Case of a 2D variable
!---------------------

        else if (dim.eq.2) then

#ifdef CPP_PARA
          ! Gather field on a "global" (without redundant longitude) array
          px2(:)=px(:,1)
          call Gather(px2,dx2_glop)
!$OMP MASTER
          if (is_mpi_root) then
            call Grid1Dto2D_glo(dx2_glop,dx2_glo)
            ! copy dx2_glo() to dx2(:) and add redundant longitude
            dx2(1:nbp_lon,:)=dx2_glo(1:nbp_lon,:)
            dx2(nbp_lon+1,:)=dx2(1,:)
          endif
!$OMP END MASTER
!$OMP BARRIER
#else

!         Passage variable physique -->  physique dynamique
!         recast (copy) variable from physics grid to dynamics grid
          IF (klon_glo>1) THEN ! General case
             DO i=1,nbp_lon+1
                dx2(i,1)=px(1,1)
                dx2(i,nbp_lat)=px(ngrid,1)
             ENDDO
             DO j=2,nbp_lat-1
                ig0= 1+(j-2)*nbp_lon
                DO i=1,nbp_lon
                   dx2(i,j)=px(ig0+i,1)
                ENDDO
                dx2(nbp_lon+1,j)=dx2(1,j)
             ENDDO
          ELSE ! 1D model case
            dx2_1d=px(1,1)
          ENDIF
#endif

          if (is_master) then
           ! only the master writes to output
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
              write (*,*) "DIAGFI: creating variable ",trim(nom)

              call def_var(nid,nom,titre,unite,3,id,varid,ierr)

           else
             if (ntime==0) then
              write(*,*) "DIAGFI Error: failed creating variable ",
     &                   trim(nom)
              write(*,*) "it seems it already exists!"
              stop
             endif
           endif

           corner(1)=1
           corner(2)=1
           corner(3)=ntime
           IF (klon_glo==1) THEN
             edges(1)=1
           ELSE
             edges(1)=nbp_lon+1
           ENDIF
           edges(2)=nbp_lat
           edges(3)=1


!#ifdef NC_DOUBLE
!           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx2) 
!#else         
           IF (klon_glo>1) THEN ! General case
             ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx2)
           ELSE
             ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx2_1d)
           ENDIF
!#endif     

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with dx2: ",trim(nom)
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
              stop
           endif 

          endif !of if (is_master)

!Case of a 1D variable (ie: a column)
!---------------------------------------------------

       else if (dim.eq.1) then
         if (is_parallel) then
           write(*,*) "writediagfi error: dim=1 not implemented ",
     &                 "in parallel mode. Problem for ",trim(nom)
           stop
         endif
!         Passage variable physique -->  physique dynamique
!         recast (copy) variable from physics grid to dynamics grid
          do l=1,nbp_lev
            dx1(l)=px(1,l)
          enddo
          
          ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"altitude",id(1))
              ierr= NF_INQ_DIMID(nid,"Time",id(2))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",trim(nom)

              call def_var(nid,nom,titre,unite,2,id,varid,ierr)
              
           else
             if (ntime==0) then
              write(*,*) "DIAGFI Error: failed creating variable ",
     &                   trim(nom)
              write(*,*) "it seems it already exists!"
              stop
             endif
           endif
           
           corner(1)=1
           corner(2)=ntime
           
           edges(1)=nbp_lev
           edges(2)=1
!#ifdef NC_DOUBLE
!           ierr= NF_PUT_VARA_DOUBLE(nid,varid,corner,edges,dx1)
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx1)
!#endif

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR problem in writediagfi"
              write(*,*) "***** with dx1: ",trim(nom)
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
              stop
           endif 

!Case of a 0D variable (ie: a time-dependent scalar)
!---------------------------------------------------

        else if (dim.eq.0) then

           dx0 = px (1,1)

          if (is_master) then
           ! only the master writes to output
           ierr= NF_INQ_VARID(nid,nom,varid)
           if (ierr /= NF_NOERR) then
! corresponding dimensions
              ierr= NF_INQ_DIMID(nid,"Time",id(1))

! Create the variable if it doesn't exist yet

              write (*,*) "=========================="
              write (*,*) "DIAGFI: creating variable ",trim(nom)

              call def_var(nid,nom,titre,unite,1,id,varid,ierr)

           else
             if (ntime==0) then
              write(*,*) "DIAGFI Error: failed creating variable ",
     &                   trim(nom)
              write(*,*) "it seems it already exists!"
              stop
             endif
           endif

           corner(1)=ntime
           edges(1)=1

!#ifdef NC_DOUBLE
!           ierr = NF_PUT_VARA_DOUBLE (nid,varid,corner,edges,dx0)  
!#else
           ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,dx0)
!#endif
           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi"
              write(*,*) "***** with dx0: ",trim(nom)
              write(*,*) 'ierr=', ierr,": ",NF_STRERROR(ierr)
              stop
           endif 

          endif !of if (is_master)

        endif ! of if (dim.eq.3) elseif(dim.eq.2)...

      endif ! of if ( MOD(zitau+1,irythme) .eq.0.)

      if (is_master) then
        ierr= NF_CLOSE(nid)
      endif

      end
