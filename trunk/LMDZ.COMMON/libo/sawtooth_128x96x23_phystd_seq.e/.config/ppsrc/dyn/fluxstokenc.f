










!
! $Id: fluxstokenc.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )
      write(lunout,*)
     & 'fluxstokenc: Needs IOIPSL to function'
! of #ifdef CPP_IOIPSL
      RETURN
      END
