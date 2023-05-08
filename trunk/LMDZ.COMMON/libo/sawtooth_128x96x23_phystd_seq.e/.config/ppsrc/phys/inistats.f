










      subroutine inistats(ierr)

      use statto_mod, only: istats,istime
      use mod_phys_lmdz_para, only : is_master
      USE vertical_layers_mod, ONLY: ap,bp,aps,bps,preff,
     &                               pseudoalt,presnivs
      USE nrtype, ONLY: pi
      USE time_phylmdz_mod, ONLY: daysec,dtphys
      USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev
      implicit none

      include "netcdf.inc"

      integer,intent(out) :: ierr
      integer :: nid
      integer :: l,nsteppd
      real, dimension(nbp_lev) ::  sig_s
      real,allocatable :: lon_reg_ext(:) ! extended longitudes
      integer :: idim_lat,idim_lon,idim_llm,idim_llmp1,idim_time
      real, dimension(istime) :: lt
      integer :: nvarid


      IF (nbp_lon*nbp_lat==1) THEN
        ! 1D model
        ALLOCATE(lon_reg_ext(1))
      ELSE
        ! 3D model
        ALLOCATE(lon_reg_ext(nbp_lon+1))
      ENDIF
      
      write (*,*) 
      write (*,*) '                        || STATS ||'
      write (*,*) 
      write (*,*) 'daysec',daysec
      write (*,*) 'dtphys',dtphys
      nsteppd=nint(daysec/dtphys)
      write (*,*) 'nsteppd=',nsteppd
      if (abs(float(nsteppd)-daysec/dtphys).gt.1.e-8*daysec)
     &   stop'Dans Instat:  1jour .ne. n pas physiques'

      if(mod(nsteppd,istime).ne.0)
     &   stop'Dans Instat:  1jour .ne. n*istime pas physiques'

      istats=nsteppd/istime
      write (*,*) 'istats=',istats
      write (*,*) 'Storing ',istime,'times per day'
      write (*,*) 'thus every ',istats,'physical timestep '
      write (*,*) 

      do l= 1, nbp_lev
         sig_s(l)=((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
         pseudoalt(l)=-10.*log(presnivs(l)/preff)   
      enddo
      
      lon_reg_ext(1:nbp_lon)=lon_reg(1:nbp_lon)
      IF (nbp_lon*nbp_lat/=1) THEN
        ! In 3D, add extra redundant point (180 degrees,
        ! since lon_reg starts at -180)
        lon_reg_ext(nbp_lon+1)=-lon_reg_ext(1)
      ENDIF

      if (is_master) then
      ! only the master needs do this

      ierr = NF_CREATE("stats.nc",IOR(NF_CLOBBER,NF_64BIT_OFFSET),nid)
      if (ierr.ne.NF_NOERR) then
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif

      ierr = NF_DEF_DIM (nid, "latitude", nbp_lat, idim_lat)
      IF (nbp_lon*nbp_lat==1) THEN
        ierr = NF_DEF_DIM (nid, "longitude", 1, idim_lon)
      ELSE
        ierr = NF_DEF_DIM (nid, "longitude", nbp_lon+1, idim_lon)
      ENDIF
      ierr = NF_DEF_DIM (nid, "altitude", nbp_lev, idim_llm)
      ierr = NF_DEF_DIM (nid, "llmp1", nbp_lev+1, idim_llmp1)
      ierr = NF_DEF_DIM (nid, "Time", NF_UNLIMITED, idim_time)

      ierr = NF_ENDDEF(nid)
      call def_var_stats(nid,"Time","Time",
     &            "days since 0000-00-0 00:00:00",1,
     &            idim_time,nvarid,ierr)
! Time is initialised later by mkstats subroutine

      call def_var_stats(nid,"latitude","latitude",
     &            "degrees_north",1,idim_lat,nvarid,ierr)
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lat_reg/pi*180)
      call def_var_stats(nid,"longitude","East longitude",
     &            "degrees_east",1,idim_lon,nvarid,ierr)
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lon_reg_ext/pi*180)

! Niveaux verticaux, aps et bps
      ierr = NF_REDEF (nid)
      ierr = NF_DEF_VAR (nid,"altitude", NF_DOUBLE, 1,idim_llm,nvarid)
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name",8,"altitude")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'units',2,"km")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'positive',2,"up")
      ierr = NF_ENDDEF(nid)
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,pseudoalt)
      call def_var_stats(nid,"aps","hybrid pressure at midlayers"
     & ," ",1,idim_llm,nvarid,ierr)
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)

      call def_var_stats(nid,"bps","hybrid sigma at midlayers"
     & ," ",1,idim_llm,nvarid,ierr)
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)

      ierr=NF_CLOSE(nid)

      endif ! of if (is_master)
      end subroutine inistats

