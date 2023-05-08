      SUBROUTINE iniconst

      use control_mod, only: iphysiq, dissip_period
      USE comconst_mod, ONLY: im,jm,lllm,imp1,jmp1,lllmp1,lllmm1,
     .			dtdiss,dtvr,dtphys,pi,r,cpp,kappa,unsim
      IMPLICIT NONE
c
c      P. Le Van
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------
c
#include "dimensions.h"
#include "paramet.h"
!#include "control.h"


      INTEGER :: l
c
c
      EXTERNAL disvert
c
c-----------------------------------------------------------------------
c   dimension des boucles:
c   ----------------------

      im      = iim
      jm      = jjm
      lllm    = llm
      imp1    = iim 
      jmp1    = jjm + 1
      lllmm1  = llm - 1
      lllmp1  = llm + 1

c-----------------------------------------------------------------------

      dtdiss  = dissip_period * dtvr
      dtphys  = iphysiq * dtvr
      unsim   = 1./iim
      pi      = 2.*ASIN( 1. )

!      print*,'iphysiq=',iphysiq
!      print*,'dtvr=',dtvr
!      print*,'dtphys=',dtphys
!      stop

c-----------------------------------------------------------------------
c

      r       = cpp * kappa

      PRINT*,' R  CPP  Kappa ',  r , cpp,  kappa
c
c-----------------------------------------------------------------------
      CALL disvert
c
c
      RETURN
      END
