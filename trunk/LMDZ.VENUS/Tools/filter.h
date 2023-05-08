! Tuning parameters for fft.F90

! Low  cutting frequency, in Hz    : fcoup1
 real, parameter :: fcoup1=1.e-6

! High cutting frequency, in Hz    : fcoup2
 real, parameter :: fcoup2=3.5e-6

! Half-width of the filters, in Hz : width
 real, parameter :: width=4.e-7

! Choice of output files:
!                                  (U,     V,      W,     T)
 logical,dimension(4) :: ok_out=(/.true.,.true.,.false.,.true./)

