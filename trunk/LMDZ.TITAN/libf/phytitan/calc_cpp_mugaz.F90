      subroutine calc_cpp_mugaz

!==================================================================
!     Purpose
!     -------
!     Check to see if the atmospheric specific heat capacity and
!     mean molar mass for the gas mixture defined in gases.def
!     corresponds to what we're using. If it doesn't, abort run
!     unless option 'check_cpp_match' is set to false in 
!     callphys.def.
!
!     Authors
!     ------- 
!     Robin Wordsworth (2009)
!     A. Spiga: make the routine OK with latest changes in rcm1d
!
!==================================================================

      use gases_h
      use comcstfi_mod, only: cpp, mugaz
      use callkeys_mod, only: check_cpp_match,force_cpp
      implicit none

      real cpp_c   
      real mugaz_c

      integer igas

      cpp_c   = 0.0
      mugaz_c = 0.0


! compute mugaz
      do igas=1,ngasmx

         ! all values at 300 K from Engineering Toolbox
         if(igas.eq.igas_N2)then
            mugaz_c = mugaz_c + 28.01*gfrac(igas,nivref)
         elseif(igas.eq.igas_H2)then
            mugaz_c = mugaz_c + 2.01*gfrac(igas,nivref)
         elseif(igas.eq.igas_CH4)then
            mugaz_c = mugaz_c + 16.04*gfrac(igas,nivref)
         else
            print*,'Error in calc_cpp_mugaz: Gas species not recognised!'
            call abort
         endif
      enddo

!compute cpp
      do igas=1,ngasmx

         ! all values at 300 K from Engineering Toolbox
         if(igas.eq.igas_N2)then
            cpp_c   = cpp_c   + 1.040*gfrac(igas,nivref)*28.01/mugaz_c
         elseif(igas.eq.igas_H2)then
            cpp_c   = cpp_c   + 14.31*gfrac(igas,nivref)*2.01/mugaz_c
         elseif(igas.eq.igas_CH4)then
            cpp_c   = cpp_c   + 2.226*gfrac(igas,nivref)*16.04/mugaz_c
         else
            print*,'Error in calc_cpp_mugaz: Gas species not recognised!'
            call abort
         endif
      enddo

      cpp_c = 1000.0*cpp_c

      print*,'Cp in calc_cpp_mugaz is ',cpp_c,'J kg^-1 K^-1'
      print*,'Mg in calc_cpp_mugaz is ',mugaz_c,'amu'
      print*,'Predefined Cp in physics is ',cpp,'J kg^-1 K^-1'
      print*,'Predefined Mg in physics is ',mugaz,'amu'

      if (check_cpp_match) then
         print*,'REQUEST TO CHECK cpp_match :'
         if((abs(1.-cpp/cpp_c).gt.1.e-6) .or.  &
              (abs(1.-mugaz/mugaz_c).gt.1.e-6)) then
            ! Ehouarn: tolerate a small mismatch between computed/stored values
            print*,'--> Values do not match!'
            print*,'    Either adjust cpp / mugaz via newstart to calculated values,'
            print*,'    or set check_cpp_match to .false. in callphys.def.'
            stop
         else
            print*,'--> OK. Settings match composition.'
         endif
      endif

      if (.not.force_cpp) then
          print*,'*** Setting cpp & mugaz to computations in calc_cpp_mugaz.'
          mugaz = mugaz_c
          cpp = cpp_c
      else
          print*,'*** Setting cpp & mugaz to predefined values.'
      endif


      return
    end subroutine calc_cpp_mugaz
