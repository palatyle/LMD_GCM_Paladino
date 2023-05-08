subroutine call_profilgases(nlayer)

  use gases_h

  implicit none


  !============================================================================
  !
  !     Purpose
  !     -------
  !     Set atmospheric composition for Titan main gases (N2, CH4, H2) from :
  !       - Interactive tracers if we run with chemistry
  !       or
  !       - CH4 vertical profile file (profile.def) at firstcall only
  !
  !     Author
  !     ------
  !     J. Vatant d'Ollone (2017)
  !
  !============================================================================


  !--------------------------
  ! 0. Declarations
  ! -------------------------

  integer, intent(in) :: nlayer
  
  
  logical, save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)


  integer ilay, ierr

  
  ! -----------------------------------------------
  ! 1. First case : no interactive CH4 tracer
  ! -----------------------------------------------

  if (firstcall) then

!$OMP MASTER
     
     ! Load CH4 vertical profile from file 'profile.def'
     open(90,file='profile.def',status='old',form='formatted',iostat=ierr)

     if (ierr.eq.0) then
        write(*,*) "call_profilgases.F90: reading file profile.def"
        read(90,*) ! header

        write(*,*) "call_profilgases.F90: reading CH4 vertical profile..."

        do ilay=1,nlayer
           read(90,*,iostat=ierr) gfrac(igas_CH4,ilay)
           if (ierr.ne.0) then
              write(*,*) 'call_profilgases.F90: error reading CH4 vertical profile... aborting'
              call abort
           endif          
        enddo        
        
        ! Then set H2 (fixed) and N2 (what remains)

        gfrac(igas_H2,:) = 1.0E-3
        gfrac(igas_N2,:) = 1.0 - ( gfrac(igas_CH4,:) + gfrac(igas_H2,:) )
        
     else
        write(*,*) 'Cannot find required file "profile.def"'
        call abort
     endif

     close(90)
!$OMP END MASTER
!$OMP BARRIER
     
     firstcall=.false.
  endif ! if (firstcall)


  ! -----------------------------------------------
  ! 2. Second case : interactive CH4 tracer
  ! -----------------------------------------------
  
end subroutine call_profilgases
