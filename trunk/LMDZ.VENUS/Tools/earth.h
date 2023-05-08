! Parameters needed to integrate hydrostatic equation:

real,parameter :: g0=9.80616
!g0: exact mean gravity at radius=6371.22km

real,parameter :: a0=6371.22E3
!a0: 'mean' radius=6371.22km

real,parameter :: R0=287.1 ! molecular gas constant

real,parameter :: psref=1.0e5 ! reference pressure at surface (Pa)

real,parameter :: omega=7.29e-5 ! angular rotation speed (s-1)

real,parameter :: localday=86400. ! local day (s)

character (len=5),parameter :: planet="Earth"

real,parameter :: cp0=1004.64 !doit etre egal a cpp (dyn) et RCPD (phy)
real,parameter :: t0=0.
real,parameter :: nu=0.

