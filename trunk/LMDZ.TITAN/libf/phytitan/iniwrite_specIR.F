      SUBROUTINE iniwrite_specIR(nid,idayref,area,nbplon,nbplat)

      use radinc_h, only: L_NSPECTI
      use radcommon_h, only: WNOI,DWNI
      use comcstfi_mod, only: rad, omeg, g, mugaz, rcp, pi
      use time_phylmdz_mod, only: daysec, dtphys
!      USE logic_mod, ONLY: fxyhypb,ysinus
!      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy
!      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev

      implicit none

c=======================================================================
c
c   Auteur:  L. Fairhead  ,  P. Le Van, Y. Wanherdrick, F. Forget
c   -------
c
c   Objet:
c   ------
c
c   'Initialize' the diagfi_spec.nc file: write down dimensions as well
c   as time-independent fields (e.g: geopotential, mesh area, ...)
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      include "netcdf.inc"

c   Arguments:
c   ----------

      integer,intent(in) :: nid        ! NetCDF file ID
      INTEGER*4,intent(in) :: idayref  ! date (initial date for this run)
      real,intent(in) :: area(nbplon,nbplat) ! mesh area (m2)
      integer,intent(in) :: nbplon,nbplat ! sizes of area

c   Local:
c   ------
      INTEGER length,l
      parameter (length = 100)
      REAL tab_cntrl(length) ! run parameters are stored in this array
      INTEGER ierr
      REAl,ALLOCATABLE :: lon_reg_ext(:) ! extended longitudes

      integer :: nvarid,idim_index,idim_rlonu,idim_rlonv
      integer :: idim_rlatu,idim_rlatv,idim_llmp1,idim_llm
 !     integer :: idim_nsoilmx ! "subsurface_layers" dimension ID #
      integer :: idim_bandsIR ! "IR Wavenumber" dimension ID #
      integer, dimension(2) :: id  

c-----------------------------------------------------------------------

      IF (nbp_lon*nbp_lat==1) THEN
        ! 1D model
        ALLOCATE(lon_reg_ext(1))
      ELSE
        ! 3D model
        ALLOCATE(lon_reg_ext(nbp_lon+1))
      ENDIF

      DO l=1,length
         tab_cntrl(l)=0.
      ENDDO
      tab_cntrl(1)  = FLOAT(nbp_lon)
      tab_cntrl(2)  = FLOAT(nbp_lat-1)
      tab_cntrl(3)  = FLOAT(nbp_lev)
      tab_cntrl(4)  = FLOAT(idayref)
      tab_cntrl(5)  = rad
      tab_cntrl(6)  = omeg
      tab_cntrl(7)  = g
      tab_cntrl(8)  = mugaz
      tab_cntrl(9)  = rcp
      tab_cntrl(10) = daysec
      tab_cntrl(11) = dtphys
!      tab_cntrl(12) = etot0
!      tab_cntrl(13) = ptot0
!      tab_cntrl(14) = ztot0
!      tab_cntrl(15) = stot0
!      tab_cntrl(16) = ang0
c
c    ..........    P.Le Van  ( ajout le 8/04/96 )    .........
c         .....        parametres  pour le zoom          ......   
!      tab_cntrl(17)  = clon
!      tab_cntrl(18)  = clat
!      tab_cntrl(19)  = grossismx
!      tab_cntrl(20)  = grossismy
c
c     .....   ajout  le 6/05/97 et le 15/10/97  .......
c
!      IF ( fxyhypb )   THEN
!        tab_cntrl(21) = 1.
!        tab_cntrl(22) = dzoomx
!        tab_cntrl(23) = dzoomy
!      ELSE
!        tab_cntrl(21) = 0.
!        tab_cntrl(22) = dzoomx
!        tab_cntrl(23) = dzoomy
!        tab_cntrl(24) = 0.
!        IF( ysinus )  tab_cntrl(24) = 1.
!      ENDIF

c    .........................................................

! Define dimensions
    
      ierr = NF_REDEF (nid)

      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "latitude", nbp_lat, idim_rlatu)
      IF (nbp_lon*nbp_lat==1) THEN
        ierr = NF_DEF_DIM (nid, "longitude", 1, idim_rlonv)
      ELSE
        ierr = NF_DEF_DIM (nid, "longitude", nbp_lon+1, idim_rlonv)
      ENDIF
      ierr = NF_DEF_DIM (nid, "IR_Wavenumber",L_NSPECTI,idim_bandsIR)

      ierr = NF_ENDDEF(nid)

c  Contol parameters for this run
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "controle", NF_DOUBLE, 1, 
     .       idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "controle", NF_FLOAT, 1, 
     .       idim_index,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"title", 18,
     .                       "Control parameters")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif

c --------------------------
c  longitudes and latitudes


      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr =NF_DEF_VAR(nid, "latitude", NF_DOUBLE, 1, idim_rlatu,nvarid)
#else
      ierr =NF_DEF_VAR (nid, "latitude", NF_FLOAT, 1, idim_rlatu,nvarid)
#endif
      ierr =NF_PUT_ATT_TEXT(nid,nvarid,'units',13,"degrees_north")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 14,
     .      "North latitude")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lat_reg/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,lat_reg/pi*180)
#endif
c
c --------------------------
      lon_reg_ext(1:nbp_lon)=lon_reg(1:nbp_lon)
      !add extra redundant point (180 degrees, since lon_reg starts at -180
      IF (nbp_lon*nbp_lat/=1) THEN
        ! In 3D, add extra redundant point (180 degrees,
        ! since lon_reg starts at -180)
        lon_reg_ext(nbp_lon+1)=-lon_reg_ext(1)
      ENDIF
      
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
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,lon_reg_ext/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,lon_reg_ext/pi*180)
#endif
c

!-------------------------------
! Number of bands in the IR
!-------------------------------

      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"IR_Wavenumber",NF_DOUBLE,1,
     .                                idim_bandsIR,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"IR_Wavenumber",NF_FLOAT,1,
     .                                idim_bandsIR,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 34,
     .                        "Band mid frequency in the infrared")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",5,"cm^-1")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,dble(WNOI))
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,real(WNOI))
#endif

!-------------------------------
! Width of bands in the IR
!-------------------------------

      ierr=NF_REDEF (nid) ! Enter NetCDF (re-)define mode
      ! define variable
#ifdef NC_DOUBLE
      ierr=NF_DEF_VAR(nid,"IR_Bandwidth",NF_DOUBLE,1,
     .                                idim_bandsIR,nvarid)
#else
      ierr=NF_DEF_VAR(nid,"IR_Bandwidth",NF_FLOAT,1,
     .                                idim_bandsIR,nvarid)
#endif
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"long_name", 25,
     .                        "Bandwidth in the infrared")
      ierr=NF_PUT_ATT_TEXT (nid,nvarid,"units",5,"cm^-1")
      ierr=NF_ENDDEF(nid) ! Leave NetCDF define mode
      ! write variable
#ifdef NC_DOUBLE
      ierr=NF_PUT_VAR_DOUBLE (nid,nvarid,dble(DWNI))
#else
      ierr=NF_PUT_VAR_REAL (nid,nvarid,real(DWNI))
#endif
c
c --------------------------
c  Mesh area and conversion coefficients cov. <-> contra. <--> natural


      id(1)=idim_rlonv
      id(2)=idim_rlatu
c
c --------------------------
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid, "aire", NF_DOUBLE, 2, id,nvarid)
#else
      ierr = NF_DEF_VAR (nid, "aire", NF_FLOAT, 2, id,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Mesh area")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,area)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,area)
#endif


      END
