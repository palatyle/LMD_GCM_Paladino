      SUBROUTINE iniconst

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
#include "comconst.h"
#include "temps.h"
#include "control.h"
#include "comvert.h"


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

      dtdiss  = idissip * dtvr
      dtphys  = iphysiq * dtvr
      unsim   = 1./iim
      pi      = 2.*ASIN( 1. )

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
