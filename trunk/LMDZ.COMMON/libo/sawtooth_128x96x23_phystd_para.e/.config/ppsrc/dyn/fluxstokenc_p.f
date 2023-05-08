










!
! $Id: fluxstokenc_p.F 1454 2010-11-18 12:01:24Z fairhead $
!
      SUBROUTINE fluxstokenc_p(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )
      write(lunout,*)
     & 'fluxstokenc: Needs IOIPSL to function'
! of #ifdef CPP_IOIPSL
      RETURN
      END
