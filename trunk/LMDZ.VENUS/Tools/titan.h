! Parameters needed to integrate hydrostatic equation:

real,parameter :: g0=1.35
!g0: exact mean gravity at radius=2575.km

real,parameter :: a0=2575.E3
!a0: 'mean' radius=2575.km

real,parameter :: R0=296.9 ! molecular gas constant

real,parameter :: psref=1.4e5 ! reference pressure at surface (Pa)

real,parameter :: omega=4.5238899E-06 ! angular rotation speed (s-1)

real,parameter :: localday=1.37889e6 ! local day (s)

character (len=5),parameter :: planet="Titan"

real,parameter :: cp0=1039. !doit etre egal a cpp (dyn) et RCPD (phy)
real,parameter :: t0=0.
real,parameter :: nu=0.
