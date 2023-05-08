
!  statto:
!     This include file controls the production of statistics.
!     Some variables could be set in a namelist, but it is easier to
!     do it here since arrays can then be dimensioned using parameters
!     and values shouldn't have to change too often.   SRL

!     Calculate stats every istats physics timesteps, starting at first
!     call.  If istats=0 then don't do statistics at all.  Check value
!     if number of physics timesteps changes.
	integer istats

!     Calculate itime independent sums and sums of squares,
!     example, istat=1,istime=1 gives a single time mean
	integer, parameter :: istime=12

!     Number of 2D and 3D variables on which to do statistics.
	integer n2dvar, n3dvar
	parameter (n2dvar = 8, n3dvar = 5)

!     Units for writing stats header and data
	integer usdata

!     count tab to know the variable record
        integer count(istime)

!     Record of the number of stores made for each time.
	integer nstore(istime)

! Size of the "controle" array
        integer, parameter :: cntrlsize=15

!       common /sttcom/ dummy,nstore,istats,usdata
        common /sttcom/ nstore,istats,usdata,count
