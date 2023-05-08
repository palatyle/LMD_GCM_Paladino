










!
! $Id: iniconst.F90 1625 2012-05-09 13:14:48Z lguez $
!
SUBROUTINE iniconst

  USE control_mod
  ! if not using IOIPSL, we still need to use (a local version of) getin
  use ioipsl_getincom
  USE comvert_mod, ONLY: disvert_type,pressure_exner
  USE comconst_mod, ONLY: im,jm,lllm,imp1,jmp1,lllmm1,lllmp1,	&
		dtphys,dtvr,unsim,r,cpp,kappa,pi

  IMPLICIT NONE
  !
  !      P. Le Van
  !
  !   Declarations:
  !   -------------
  !
  include "dimensions.h"
  include "paramet.h"
  include "iniprint.h"

  character(len=*),parameter :: modname="iniconst"
  character(len=80) :: abort_message
  !
  !
  !
  !-----------------------------------------------------------------------
  !   dimension des boucles:
  !   ----------------------

  im      = iim
  jm      = jjm
  lllm    = llm
  imp1    = iim 
  jmp1    = jjm + 1
  lllmm1  = llm - 1
  lllmp1  = llm + 1

  !-----------------------------------------------------------------------

  dtphys  = iphysiq * dtvr
  unsim   = 1./iim
  pi      = 2.*ASIN( 1. )

  !-----------------------------------------------------------------------
  !

  r       = cpp * kappa

  write(lunout,*) trim(modname),': R  CP  Kappa ',r,cpp,kappa
  !
  !-----------------------------------------------------------------------

  ! vertical discretization: default behavior depends on planet_type flag
  if (planet_type=="earth") then
     disvert_type=1
  else
     disvert_type=2
  endif
  ! but user can also specify using one or the other in run.def:
  call getin('disvert_type',disvert_type)
  write(lunout,*) trim(modname),': disvert_type=',disvert_type

  pressure_exner = disvert_type == 1 ! default value
  call getin('pressure_exner', pressure_exner)

  if (disvert_type==1) then
     ! standard case for Earth (automatic generation of levels)
     call disvert()
  else if (disvert_type==2) then
     ! standard case for planets (levels generated using z2sig.def file)
     call disvert_noterre
  else
     write(abort_message,*) "Wrong value for disvert_type: ", disvert_type
     call abort_gcm(modname,abort_message,0)
  endif

END SUBROUTINE iniconst
