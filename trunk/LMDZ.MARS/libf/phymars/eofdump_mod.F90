module eofdump_mod
! this module controls the production of data for EOFs
! it won't work if run in parallel (but it's OK, we don't use it anymore...)
! Mainly kept for reference.
implicit none
! Dump profiles for EOFs every ieofs physics timesteps,
! starting at first call;
integer :: ieofs
! Dump profiles every eofskip points in each direction
! on the model grid.
integer, parameter :: eofskip = 4
! Units for writing EOF header and data:
integer, parameter :: uehead = 82, uedata = 83

contains

      subroutine eofdump(ngrid,nlayer,u,v,t,rho,ps)

      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      implicit none
!
!     Dumps profiles for calculation of variability EOFs 
!     Modified to include rho, FF 09/2004
!     Corrected small bug in sampling rate/count, EM 11/2007
!
!

      integer,intent(in) :: ngrid ! total number of physics grid points
      integer,intent(in) :: nlayer ! number of atmospheric layers
      real,intent(in) :: u(ngrid,nlayer)
      real,intent(in) :: v(ngrid,nlayer)
      real,intent(in) :: t(ngrid,nlayer)
      real,intent(in) :: rho(ngrid,nlayer)
      real,intent(in) :: ps(ngrid)
      integer,save :: count=0
      integer i,j,l, ig

      LOGICAL,SAVE :: firstcall=.true.

!-------------------------------------------------------
!     Initialization at first call:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (firstcall) THEN
        write(*,*) 'CALL ineofdump'
        CALL ineofdump(ngrid,nlayer)
        firstcall=.false.
      END IF

!-------------------------------------------------------
!     Dumps every ieofs physics timesteps
!
!      write(*,*)'eofdump:count=',count,' ps(1)=',ps(1)
!      if ((ieofs.gt.0).and.(mod(count,ieofs).eq.0)) then
      if (mod(count+1,ieofs).eq.0) then
!        write(*,*)'eofdump: dump --> ps(1)=',ps(1)
        do i=1,nbp_lon,eofskip
          do j=1+eofskip/2,nbp_lat,eofskip
            ig = 1+ (j-2)*nbp_lon +i
#ifdef NC_DOUBLE
            write(uedata) (real(u(ig,l)),l=1,nlayer)
            write(uedata) (real(v(ig,l)),l=1,nlayer)
            write(uedata) (real(t(ig,l)),l=1,nlayer)
            write(uedata) (real(rho(ig,l)),l=1,nlayer)
            write(uedata) real(ps(ig))
#else
            write(uedata) (u(ig,l),l=1,nlayer)
            write(uedata) (v(ig,l),l=1,nlayer)
            write(uedata) (t(ig,l),l=1,nlayer)
            write(uedata) (rho(ig,l),l=1,nlayer)
            write(uedata) ps(ig)
#endif
          enddo
        enddo
      endif
      count=count+1
 
      end subroutine eofdump


      subroutine ineofdump(ngrid,nlayer)

      use geometry_mod, only: longitude, latitude
      use nrtype, only: pi
      use time_phylmdz_mod, only: daysec, dtphys
      USE vertical_layers_mod, ONLY: aps,bps
      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      implicit none
!
!     Initialise dumping of profiles for EOF calculations
!

      integer,intent(in) :: ngrid ! total number of physics grid points
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer ig,i,j,l    
      logical,save :: firstcall=.true.
      integer,save :: npgrid


      if (firstcall) then
         npgrid=ngrid+2*(nbp_lon-1)
         firstcall=.false.
      endif

!
!     Set frequency for dumping at once per day
!
      ieofs=nint(daysec/dtphys)
      if (abs(float(ieofs)-daysec/dtphys).gt.1.e-8*daysec) &
         stop 'In ineofdump:  1 day .ne. n physics timesteps'
!
!     Header
!
      open(uehead,file='profiles.hdr',form='formatted')
      write(uehead,*) 0.E+0,0,0,ieofs,1,0
      write(uehead,*) nbp_lon,npgrid/nbp_lon,npgrid,nlayer

      do i=1,nbp_lon,eofskip
        do j=1+eofskip/2,nbp_lat,eofskip    
          ig = 1+ (j-2)*nbp_lon +i
          if(j.eq.1) stop 'Problem in ineofdump.F'
          if(j.eq.nbp_lat) stop 'Problem in ineofdump.F'
#ifdef NC_DOUBLE
          write(uehead,*) real(longitude(ig)*180./pi),real(latitude(ig)*180./pi)
#else
          write(uehead,*) longitude(ig)*180./pi, latitude(ig)*180./pi
#endif
!         write(*,*) 'eof grid j=',j,' lat= ', latitude(ig)*180./pi
        enddo
      enddo

#ifdef NC_DOUBLE
      write(uehead,*) real(aps)
      write(uehead,*) real(bps)
#else
      write(uehead,*) aps
      write(uehead,*) bps
#endif
      close(uehead)
!
!     Main profile file
!
      open(uedata,file='profiles.dat',form='unformatted')
      end subroutine ineofdump

end module eofdump_mod
