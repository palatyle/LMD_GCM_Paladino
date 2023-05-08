c
c  eofdump:
c	This include file controls the production of data for EOFs.
c
c	Dump profiles for EOFs every ieofs physics timesteps,
c       starting at first call.
	integer ieofs
c
c	Only dump profiles every eofskip points in each direction
c       on the model grid.
	integer eofskip
	parameter (eofskip = 4)
c
c	Units for writing EOF header and data
	integer uehead, uedata
	parameter (uehead = 82, uedata = 83)
c
	common /eofcom/ ieofs
c
      INTEGER npgrid
      PARAMETER(npgrid=ngridmx+2*(iim-1))

