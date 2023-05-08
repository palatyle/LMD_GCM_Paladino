subroutine su_gases(nlayer,tracer)

  use gases_h

  implicit none

  integer,intent(in) :: nlayer
  logical,intent(in) :: tracer
  
  integer igas, ierr
  
  !==================================================================
  !     
  !     Purpose
  !     -------
  !     Load atmospheric composition info
  !     
  !     Authors
  !     -------
  !     R. Wordsworth (2011)
  !     Allocatable arrays by A. Spiga (2011)
  !     Titan's version : J. Vatant d'Ollone (2017)
  !              * force igas_X labels and def gnom for N2,CH4,H2
  !              * get rid of variable species ( enrichment dimension will be set in corrk routines )
  !
  !==================================================================

!$OMP MASTER


  ngasmx   = 3 ! N2, CH4, H2


  ! load reference level and reference molar fractions from file 'gases.def'
  open(90,file='gases.def',status='old',form='formatted',iostat=ierr)
  
  if (ierr.eq.0) then
     write(*,*) "sugases.F90: reading file gases.def"
     read(90,*) ! header

     ! We allocate gfrac and we set gas molar fractions for the reference level only.
     ! This will be useful for the cpp_mu and rayleigh routines
     ! Other gas molar fractions are now set in routine callprofilgases
     
     write(*,*) 'sugases.F90: allocating and reading gas molar fractions from reference level in gases.def...'
     
     if(.not.allocated(gfrac)) allocate(gfrac(ngasmx,nlayer))
     
     read(90,*,iostat=ierr) nivref     
     if (ierr.ne.0) then
        write(*,*) "sugases.F90: error reading reference level"
        write(*,*) "   (first line of gases.def) "
        call abort
     endif
     
     print*, "layer", nivref, "is reference level found in gases.def ..."
     
     do igas=1,ngasmx
        read(90,*,iostat=ierr) gfrac(igas,nivref)
        if (ierr.ne.0) then
           write(*,*) 'sugases.F90: error reading reference gas molar fractions in gases.def... aborting'
           call abort
        endif
     enddo                  !of do igas=1,ngasmx


     ! We force gnom = (N2, CH4, H2) and igas_X for Titan
     
     igas_N2  = 1
     igas_CH4 = 2
     igas_H2  = 3


     if (.not.allocated(gnom)) allocate(gnom(ngasmx))
     gnom(igas_N2)  = "N2_"
     gnom(igas_CH4) = "CH4"
     gnom(igas_H2)  = "H2_"
     
     
  else
     write(*,*) 'Cannot find required file "gases.def"'
     call abort
  endif
  
  close(90)
!$OMP END MASTER
!$OMP BARRIER

end subroutine su_gases
