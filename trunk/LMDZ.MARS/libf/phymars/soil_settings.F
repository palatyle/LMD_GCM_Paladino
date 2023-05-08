      subroutine soil_settings(nid,ngrid,nsoil,tsurf,tsoil,indextime)

!      use netcdf
      use comsoil_h, only: layer, mlayer, inertiedat, volcapa
      use iostart, only: inquire_field_ndims, get_var, get_field,
     &                   inquire_field, inquire_dimension_length
      implicit none

!======================================================================
!  Author: Ehouarn Millour (07/2006)
!
!  Purpose: Read and/or initialise soil depths and properties
!
! Modifications: Aug.2010 EM : use NetCDF90 to load variables (enables using
!                      r4 or r8 restarts independently of having compiled
!                      the GCM in r4 or r8)
!                June 2013 TN : Possibility to read files with a time axis
!
!
!  This subroutine reads from a NetCDF file (opened by the caller)
!  of "startfi.nc" format.
!  The various actions and variable read/initialized are:
!  1. Check out the number of soil layers (in datafile); if it isn't equal
!     to nsoil, then some interpolation will be required
!     Also check if data in file "startfi.nc" is in older format (ie:
!     thermal inertia was depth-independent; and there was no "depth"
!     coordinate.
!     Read/build layer (and midlayer) depths
!  2. Read volumetric specific heat (or initialise it to default value)
!  3. Read Thermal inertia
!  4. Read soil temperatures
!  5. Interpolate thermal inertia and temperature on the new grid, if
!     necessary
!======================================================================

!======================================================================
!  arguments
!  ---------
!  inputs:
      integer,intent(in) :: nid	! Input Netcdf file ID 
      integer,intent(in) :: ngrid	! # of horizontal grid points
      integer,intent(in) :: nsoil	! # of soil layers
      real,intent(in) :: tsurf(ngrid)   ! surface temperature
      integer,intent(in) :: indextime	! position on time axis
!  output:
      real,intent(out) :: tsoil(ngrid,nsoil)	! soil temperature

!======================================================================
! local variables:
      integer ierr	! status (returned by NetCDF functions)
      integer nvarid	! ID of NetCDF variable
      integer dimid	! ID of NetCDF dimension
      integer dimlen	! length along the "depth" dimension
      integer ndims	! # of dimensions of read <inertiedat> data
      integer ig,iloop	! loop counters
      
      integer edges(3),corner(3) ! to read a specific time

      logical :: olddepthdef=.false. ! flag
      logical :: interpol=.false. ! flag: true if interpolation will be requiered

      ! to store "old" values
      real,dimension(:),allocatable :: surfinertia !surface thermal inertia
      real,dimension(:),allocatable :: oldmlayer
      real,dimension(:,:),allocatable :: oldinertiedat
      real,dimension(:,:),allocatable :: oldtsoil
      
      ! for interpolation
      real,dimension(:),allocatable :: oldgrid
      real,dimension(:),allocatable :: oldval
      real,dimension(:),allocatable :: newval

      real alpha,lay1 ! coefficients for building layers
      real xmin,xmax ! to display min and max of a field

      real,parameter :: default_volcapa=1.e6

      logical :: found,ok
      
!======================================================================

! 1. Depth coordinate
! -------------------

! 1.1 Start by reading how many layers of soil there are

!	ierr=nf90_inq_dimid(nid,"subsurface_layers",dimid)
!	if (ierr.ne.nf90_noerr) then
!	 write(*,*)'soil_settings: Problem reading <subsurface_layers>'
!         write(*,*)trim(nf90_strerror(ierr))
!	 call abort
!	endif

!	ierr=nf90_inquire_dimension(nid,dimid,len=dimlen)
!	if (ierr.ne.nf90_noerr) then
!	 write(*,*)'soil_settings: Problem getting <subsurface_layers>',
!     &             'length'
!	 write(*,*)trim(nf90_strerror(ierr))
!         call abort
!	endif
        dimlen=inquire_dimension_length("subsurface_layers")

        if (dimlen.ne.nsoil) then
          write(*,*)'soil_settings: Interpolation of soil temperature ',
     &              'and thermal inertia will be required!'
	! if dimlen doesn't match nsoil, then interpolation of
	! soil temperatures and thermal inertia will be requiered
	  interpol=.true.
	endif

! 1.2 Find out the # of dimensions <inertiedat> was defined as using
!     (in ye old days, thermal inertia was only given at the "surface")
      ! Look for thermal inertia data
!      ierr=nf90_inq_varid(nid,"inertiedat",nvarid)
!      if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Field <inertiedat> not found!'
!         write(*,*)trim(nf90_strerror(ierr))
!         call abort
!      endif
!
!      ! Read the # of dimensions <inertidat> was defined as using
!      ierr=nf90_inquire_variable(nid,nvarid,ndims=ndims)
!      ! if (ndims.eq.1) then we have the "old 2D-surface" format
      ndims=inquire_field_ndims("inertiedat")

! 1.3 Read depths values or set olddepthdef flag and values
      if (ndims.eq.1) then ! we know that there is none
        write(*,*)'soil_settings: no <soildepth> field expected'
	write(*,*)'building one mimicking old definitions'
        olddepthdef=.true.
	interpol=.true.
        ! allocate oldmlayer
        if (.not.allocated(oldmlayer)) then
          allocate(oldmlayer(dimlen),stat=ierr)
          if (ierr.ne.0) then
            write(*,*) 'soil_settings: failed allocation of oldmlayer!'
            stop
          endif
        endif
	do iloop=1,dimlen
	  oldmlayer(iloop)=sqrt(887.75/3.14)*((2.**(iloop-0.5))-1.)
	enddo
      else ! Look for depth 
!        ierr=nf90_inq_varid(nid,"soildepth",nvarid)
!        if (ierr.ne.nf90_noerr) then
!          write(*,*)'soil_settings: Field <soildepth> not found!'
!          write(*,*)trim(nf90_strerror(ierr))
!	  call abort
!        endif
        ! read <depth> coordinate
        if (interpol) then !put values in oldmlayer
!          ierr=nf90_get_var(nid,nvarid,oldmlayer)
!          if (ierr.ne.nf90_noerr) then
!           write(*,*)'soil_settings: Problem while reading <soildepth>'
!           write(*,*)trim(nf90_strerror(ierr))
!           call abort
!          endif
          call get_var("soildepth",oldmlayer,found)
          if (.not.found) then
            write(*,*)'soil_settings: Problem while reading <soildepth>'
          endif
        else ! put values in mlayer
!          ierr=nf90_get_var(nid,nvarid,mlayer)
!          if (ierr.ne.nf90_noerr) then
!           write(*,*)'soil_settings: Problem while reading <soildepth>'
!           write(*,*)trim(nf90_strerror(ierr))
!           call abort
!          endif
          call get_var("soildepth",mlayer,found)
          if (.not.found) then
            write(*,*)'soil_settings: Problem while reading <soildepth>'
          endif
        endif !of if (interpol)
      endif !of if (ndims.eq.1)

! 1.4 Build mlayer(), if necessary
      if (interpol) then
      ! default mlayer distribution, following a power law:
      !  mlayer(k)=lay1*alpha**(k-1/2)
        lay1=2.e-4
	alpha=2
        do iloop=0,nsoil-1
	  mlayer(iloop)=lay1*(alpha**(iloop-0.5))
	enddo
      endif

! 1.5 Build layer(); following the same law as mlayer()
      ! Assuming layer distribution follows mid-layer law:
      ! layer(k)=lay1*alpha**(k-1)
      lay1=sqrt(mlayer(0)*mlayer(1))
      alpha=mlayer(1)/mlayer(0)
      do iloop=1,nsoil
        layer(iloop)=lay1*(alpha**(iloop-1))
      enddo

! 2. Volumetric heat capacity (note: it is declared in comsoil_h)
! ---------------------------
! "volcapa" is (so far) 0D and written in "controle" table of startfi file 
! volcapa is read or set when "controle" is read (see tabfi.F)
! Just in case, we check here that it is not zero. If it is, we
! set it to "default_volcapa"

      if (volcapa.le.0.0) then
        write(*,*)'soil_settings: Warning, volcapa = ',volcapa
	write(*,*)'               That doesn t seem right'
        write(*,*)'        Initializing Volumetric heat capacity to ',
     &             default_volcapa
	volcapa=default_volcapa
      endif
! Look for it
!      ierr=NF_INQ_VARID(nid,"volcapa",nvarid)
!      if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Field <volcapa> not found!'
!         write(*,*)'Initializing Volumetric heat capacity to ',
!     &             default_volcapa
!         volcapa=default_volcapa
!      else
!#ifdef NC_DOUBLE
!       ierr = NF_GET_VAR_DOUBLE(nid,nvarid,volcapa)
!#else
!       ierr = NF_GET_VAR_REAL(nid,nvarid,volcapa)
!#endif
!        if (ierr.ne.nf90_noerr) then
!	  write(*,*)'soil_settings: Problem while reading <volcapa>'
!         call abort
!	endif
!      endif

! 3. Thermal inertia (note: it is declared in comsoil_h)
! ------------------

! 3.1 Look (again) for thermal inertia data (to reset nvarid)
!      ierr=nf90_inq_varid(nid,"inertiedat",nvarid)
!      if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Field <inertiedat> not found!'
!         write(*,*)trim(nf90_strerror(ierr))
!         call abort
!      endif

! 3.2 Knowing the # of dimensions <inertidat> was defined as using,
!     read/build thermal inertia

      if (ndims.eq.1) then ! "old 2D-surface" format
       write(*,*)'soil_settings: Thermal inertia is only given as surfac
     &e data!'
       ! Read Surface thermal inertia
       allocate(surfinertia(ngrid))
!       ierr=nf90_get_var(nid,nvarid,surfinertia)
!        if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Problem while reading <inertiedat>'
!         write(*,*)trim(nf90_strerror(ierr))
!         call abort
!        endif
       call get_field("inertiedat",surfinertia,found)
       if (.not.found) then
         write(*,*) "soil_settings: Failed loading <inertiedat>"
         call abort
       endif
       
       write(*,*)' => Building soil thermal inertia (using reference sur
     &face thermal inertia)'
       do iloop=1,nsoil
         inertiedat(:,iloop)=surfinertia(:)
       enddo
       deallocate(surfinertia)

      else ! "3D surface+depth" format
       if (interpol) then ! put values in oldinertiedat
         if (.not.allocated(oldinertiedat)) then
           allocate(oldinertiedat(ngrid,dimlen),stat=ierr)
           if (ierr.ne.0) then
            write(*,*) 'soil_settings: failed allocation of ',
     &                 'oldinertiedat!'
            stop
           endif
         endif ! of if (.not.allocated(oldinertiedat))
!         ierr=nf90_get_var(nid,nvarid,oldinertiedat)
!        if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Problem while reading <inertiedat>'
!         write(*,*)trim(nf90_strerror(ierr))
!         call abort
!        endif
        call get_field("inertiedat",oldinertiedat,found)
        if (.not.found) then
          write(*,*) "soil_settings: Failed loading <inertiedat>"
          call abort
        endif
       else ! put values in therm_i
!        ierr=nf90_get_var(nid,nvarid,inertiedat)
!        if (ierr.NE.nf90_noerr) then
!         write(*,*)'soil_settings: Problem while reading <inertiedat>'
!         write(*,*)trim(nf90_strerror(ierr))
!         call abort
         call get_field("inertiedat",inertiedat,found)
         if (.not.found) then
           write(*,*) "soil_settings: Failed loading <inertiedat>"
           call abort
         endif
!        endif
       endif ! of if (interpol)
      endif ! of if (ndims.eq.1)
      
! 4. Read soil temperatures
! -------------------------

!      ierr=nf90_inq_varid(nid,"tsoil",nvarid)
      ok=inquire_field("tsoil")
!      if (ierr.ne.nf90_noerr) then
      if (.not.ok) then
        write(*,*)'soil_settings: Field <tsoil> not found!'
	write(*,*)' => Building <tsoil> from surface values <tsurf>'
	do iloop=1,nsoil
	  tsoil(:,iloop)=tsurf(:)
	enddo
      else ! <tsoil> found
       if (interpol) then ! put values in oldtsoil
         if (.not.allocated(oldtsoil)) then
           allocate(oldtsoil(ngrid,dimlen),stat=ierr)
           if (ierr.ne.0) then
             write(*,*) 'soil_settings: failed allocation of ',
     &                  'oldtsoil!'
             stop
           endif
         endif
!        ierr=nf90_get_var(nid,nvarid,oldtsoil)
!        if (ierr.ne.nf90_noerr) then
!	 write(*,*)'soil_settings: Problem while reading <tsoil>'
!         write(*,*)trim(nf90_strerror(ierr))
!	 call abort
!	endif
         call get_field("tsoil",oldtsoil,found)
         if (.not.found) then
           write(*,*) "soil_settings: Failed loading <tsoil>"
           call abort
         endif
       else ! put values in tsoil
!        corner(1)=1
!        corner(2)=1
!        corner(3)=indextime
!        edges(1)=ngrid
!        edges(2)=nsoil
!        edges(3)=1
!        !ierr=nf90_get_var(nid,nvarid,tsoil,corner,edges)
!        ierr=nf90_get_var(nid,nvarid,tsoil)
!        if (ierr.ne.nf90_noerr) then
!	 write(*,*)'soil_settings: Problem while reading <tsoil>'
!         write(*,*)trim(nf90_strerror(ierr))
!	 call abort
!	endif
         call get_field("tsoil",tsoil,found,timeindex=indextime)
         if (.not.found) then
           write(*,*) "soil_settings: Failed loading <tsoil>"
           call abort
         endif
       endif ! of if (interpol) 
      endif! of if (.not.ok)

! 5. If necessary, interpolate soil temperatures and thermal inertias
! -------------------------------------------------------------------

      if (olddepthdef) then
      ! tsoil needs to be interpolated, but not therm_i
        allocate(oldgrid(dimlen+1))
        allocate(oldval(dimlen+1))
	allocate(newval(nsoil))

        do ig=1,ngrid
	  ! copy values
	  oldval(1)=tsurf(ig)
	  oldval(2:dimlen+1)=oldtsoil(ig,1:dimlen)
	  ! build vertical coordinate
	  oldgrid(1)=0. ! ground
	  oldgrid(2:dimlen+1)=oldmlayer(1:dimlen)*
     &                (inertiedat(ig,1)/volcapa)
	  ! interpolate
	  call interp_line(oldgrid,oldval,dimlen+1,mlayer,newval,nsoil)
	  ! copy result in tsoil
	  tsoil(ig,:)=newval(:)
	enddo

        ! cleanup
	deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)
	interpol=.false. ! no need for interpolation any more      
      endif !of if (olddepthdef)

      if (interpol) then
      write(*,*)'soil_settings: Vertical interpolation along new grid'
      ! interpolate soil temperatures and thermal inertias
        if (.not.allocated(oldgrid)) then
          allocate(oldgrid(dimlen+1))
          allocate(oldval(dimlen+1))
	  allocate(newval(nsoil))
        endif

      ! thermal inertia
        do ig=1,ngrid
	  ! copy data
	  oldval(1:dimlen)=oldinertiedat(ig,dimlen)
	  ! interpolate
	  call interp_line(oldmlayer,oldval,dimlen,mlayer,newval,nsoil)
	  !copy result in inertiedat
	  inertiedat(ig,:)=newval(:)
	enddo
        
      ! soil temperature
        ! vertical coordinate
	oldgrid(1)=0.0
	oldgrid(2:dimlen+1)=oldmlayer(1:dimlen)
        do ig=1,ngrid
	  ! data
	  oldval(1)=tsurf(ig)
	  oldval(2:dimlen+1)=oldtsoil(ig,1:dimlen)
	  ! interpolate
	  call interp_line(oldgrid,oldval,dimlen+1,mlayer,newval,nsoil)
	  ! copy result in inertiedat
	  tsoil(ig,:)=newval(:)
	enddo
	
	!cleanup
        deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)
	deallocate(oldinertiedat)
	deallocate(oldtsoil)
      endif ! of if (interpol)
      
! 6. Report min and max values of soil temperatures and thermal inertias
! ----------------------------------------------------------------------

      write(*,*)
      write(*,*)'Soil volumetric heat capacity:',volcapa

      xmin = MINVAL(inertiedat)
      xmax = MAXVAL(inertiedat)
      write(*,*)'Soil thermal inertia <inertiedat>:',xmin,xmax

      xmin = 1.0E+20
      xmax = -1.0E+20
      xmin = MINVAL(tsoil)
      xmax = MAXVAL(tsoil)
      write(*,*)'Soil temperature <tsoil>:',xmin,xmax

      end
