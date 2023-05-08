module types_asis

implicit none

integer, parameter :: jpim = 4
integer, parameter :: jprb = 8

type z3spec
	real(kind=jprb)    :: z1
        integer(kind=jpim) :: z2
	real(kind=jprb)    :: z3
        integer(kind=jpim) :: z4
	real(kind=jprb)    :: z5
        integer(kind=jpim) :: z6
end type z3spec
type z4spec
	real(kind=jprb)    :: z1
        integer(kind=jpim) :: z2
	real(kind=jprb)    :: z3
        integer(kind=jpim) :: z4
	real(kind=jprb)    :: z5
        integer(kind=jpim) :: z6
	real(kind=jprb)    :: z7
        integer(kind=jpim) :: z8
end type z4spec

! indexes for the jacobian matrix

type(z3spec), allocatable, save :: indice_phot(:)
type(z3spec), allocatable, save :: indice_3(:)
type(z4spec), allocatable, save :: indice_4(:)

end module types_asis
