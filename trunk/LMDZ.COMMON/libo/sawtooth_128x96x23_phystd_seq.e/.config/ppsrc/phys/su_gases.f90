










subroutine su_gases

  use gases_h

  implicit none

  integer igas, ierr, count

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
  !     
  !==================================================================

!$OMP MASTER
  ! load gas names from file 'gases.def'
  open(90,file='gases.def',status='old',form='formatted',iostat=ierr)
  if (ierr.eq.0) then
     write(*,*) "sugases.F90: reading file gases.def"
     read(90,*)
     read(90,*,iostat=ierr) ngasmx
     if (ierr.ne.0) then
        write(*,*) "sugases.F90: error reading number of gases"
        write(*,*) "   (first line of gases.def) "
        call abort
     endif

     print*,ngasmx, " gases found in gases.def. Allocating names and molar fractions..."

     if (.not.allocated(gnom)) allocate(gnom(ngasmx))
     do igas=1,ngasmx
        read(90,*,iostat=ierr) gnom(igas)
        if (ierr.ne.0) then
           write(*,*) 'sugases.F90: error reading gas names in gases.def...'
           call abort
        endif
     enddo                  !of do igas=1,ngasmx

     vgas=0
     if(.not.allocated(gfrac)) allocate(gfrac(ngasmx))
     do igas=1,ngasmx
        read(90,*,iostat=ierr) gfrac(igas)
        if (ierr.ne.0) then
           write(*,*) 'sugases.F90: error reading gas molar fractions in gases.def...'
           call abort
        endif

        ! find variable gas (if any)
        if(gfrac(igas).eq.-1.0)then 
           if(vgas.eq.0)then
              vgas=igas
           else
              print*,'You seem to be choosing two variable gases'
              print*,'Check that gases.def is correct'
              call abort
           endif
        endif

     enddo                  !of do igas=1,ngasmx


     ! assign the 'igas_X' labels
     count=0
     do igas=1,ngasmx
        if (trim(gnom(igas)).eq."H2_" .or. trim(gnom(igas)).eq."H2") then
           igas_H2=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."He_" .or. trim(gnom(igas)).eq."He") then
           igas_He=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."H2O") then
           igas_H2O=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."CO2") then
           igas_CO2=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."CO_" .or. trim(gnom(igas)).eq."CO") then
           igas_CO=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."N2_" .or. trim(gnom(igas)).eq."N2") then
           igas_N2=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."O2_" .or. trim(gnom(igas)).eq."O2") then
           igas_O2=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."SO2") then
           igas_SO2=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."H2S") then
           igas_H2S=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."CH4") then
           igas_CH4=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."NH3") then
           igas_NH3=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."C2H6") then
           igas_C2H6=igas
           count=count+1
        elseif (trim(gnom(igas)).eq."C2H2") then
           igas_C2H2=igas
           count=count+1
        endif
     enddo

     if(count.ne.ngasmx)then
        print*,'Mismatch between ngas and number of recognised gases in sugas_corrk.F90.'
        print*,'Either we haven`t managed to assign all the gases, or there are duplicates.'
        print*,'Please try again.'
     endif

  else
     write(*,*) 'Cannot find required file "gases.def"'
     call abort
  endif
  close(90)
!$OMP END MASTER
!$OMP BARRIER

end subroutine su_gases
