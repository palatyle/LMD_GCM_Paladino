! Parameters needed to integrate hydrostatic equation:

real,parameter :: g0=8.87
!g0: exact mean gravity at radius=6051.km

real,parameter :: a0=6051.E3
!a0: 'mean' radius=6051.km

real,parameter :: R0=191.4 ! molecular gas constant

real,parameter :: psref=9.2e6 ! reference pressure at surface (Pa)

real,parameter :: omega=2.992677e-7 ! angular rotation speed (s-1)

real,parameter :: localday=1.0087e7 ! local day (s)

character (len=5),parameter :: planet="Venus"

real,parameter :: cp0=1000. !doit etre egal a cpp (dyn) et RCPD (phy)
real,parameter :: t0=460.
real,parameter :: nu=0.35

