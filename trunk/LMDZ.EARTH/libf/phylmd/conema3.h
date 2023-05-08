!
! $Header$
!-- Modified by : Filiberti M-A 06/2005
!
      real epmax             ! 0.993
      logical ok_adj_ema      ! F
      integer iflag_clw      ! 0
	  integer iflag_cvl_sigd
      real sig1feed      ! 1.
      real sig2feed      ! 0.95

      common/comconema1/epmax,ok_adj_ema,iflag_clw,sig1feed,sig2feed
      common/comconema2/iflag_cvl_sigd

!      common/comconema/epmax,ok_adj_ema,iflag_clw
!$OMP THREADPRIVATE(/comconema1/)
!$OMP THREADPRIVATE(/comconema2/)
