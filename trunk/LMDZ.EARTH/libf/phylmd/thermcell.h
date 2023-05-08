      integer            :: iflag_thermals,nsplit_thermals
      real,parameter     :: r_aspect_thermals=2.,l_mix_thermals=30.
      real               :: tau_thermals
      integer,parameter  :: w2di_thermals=1
      integer            :: isplit

      integer            :: iflag_coupl,iflag_clos,iflag_wake
      integer            :: iflag_thermals_ed,iflag_thermals_optflux

      common/ctherm1/iflag_thermals,nsplit_thermals
      common/ctherm2/tau_thermals
      common/ctherm4/iflag_coupl,iflag_clos,iflag_wake
      common/ctherm5/iflag_thermals_ed,iflag_thermals_optflux

!$OMP THREADPRIVATE(/ctherm1/,/ctherm2/,/ctherm4/,/ctherm5/)
