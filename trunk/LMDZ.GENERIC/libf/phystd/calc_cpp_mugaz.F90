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

         if(igas.eq.vgas)then
            ! ignore variable gas in cpp calculation
         else
            ! all values at 300 K from Engineering Toolbox
            if(igas.eq.igas_CO2)then
               mugaz_c = mugaz_c + 44.01*gfrac(igas)
            elseif(igas.eq.igas_N2)then
               mugaz_c = mugaz_c + 28.01*gfrac(igas)
            elseif(igas.eq.igas_H2)then
               mugaz_c = mugaz_c + 2.01*gfrac(igas)
            elseif(igas.eq.igas_He)then
               mugaz_c = mugaz_c + 4.003*gfrac(igas)
            elseif(igas.eq.igas_H2O)then
               mugaz_c = mugaz_c + 18.02*gfrac(igas)
            elseif(igas.eq.igas_SO2)then
               mugaz_c = mugaz_c + 64.066*gfrac(igas)
            elseif(igas.eq.igas_H2S)then
               mugaz_c = mugaz_c + 34.08*gfrac(igas)
            elseif(igas.eq.igas_CH4)then
               mugaz_c = mugaz_c + 16.04*gfrac(igas)
            elseif(igas.eq.igas_NH3)then
               mugaz_c = mugaz_c + 17.03*gfrac(igas)
            elseif(igas.eq.igas_C2H6)then 
               ! C2H6 http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=28
               mugaz_c = mugaz_c + 30.07*gfrac(igas)
            elseif(igas.eq.igas_C2H2)then
               ! C2H2 http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=1
               mugaz_c = mugaz_c + 26.04*gfrac(igas)
            else
               print*,'Error in calc_cpp_mugaz: Gas species not recognised!'
               call abort
            endif
         endif

      enddo

!compute cpp
      do igas=1,ngasmx

         if(igas.eq.vgas)then
            ! ignore variable gas in cpp calculation
         else
            ! all values at 300 K from Engineering Toolbox
            if(igas.eq.igas_CO2)then
               !cpp_c   = cpp_c   + 0.744*gfrac(igas) ! @ ~210 K (better for
               !Mars conditions) 
               cpp_c   = cpp_c   + 0.846*gfrac(igas)*44.01/mugaz_c
            elseif(igas.eq.igas_N2)then
               cpp_c   = cpp_c   + 1.040*gfrac(igas)*28.01/mugaz_c
            elseif(igas.eq.igas_H2)then
               cpp_c   = cpp_c   + 14.31*gfrac(igas)*2.01/mugaz_c
            elseif(igas.eq.igas_He)then
               cpp_c   = cpp_c   + 5.19*gfrac(igas)*4.003/mugaz_c
            elseif(igas.eq.igas_H2O)then
               cpp_c   = cpp_c   + 1.864*gfrac(igas)*18.02/mugaz_c
            elseif(igas.eq.igas_SO2)then
               cpp_c   = cpp_c   + 0.64*gfrac(igas)*64.066/mugaz_c
            elseif(igas.eq.igas_H2S)then
               cpp_c   = cpp_c   + 1.003*gfrac(igas)*34.08/mugaz_c ! from wikipedia...
            elseif(igas.eq.igas_CH4)then
               cpp_c   = cpp_c   + 2.226*gfrac(igas)*16.04/mugaz_c
            elseif(igas.eq.igas_NH3)then
               cpp_c   = cpp_c   + 2.175*gfrac(igas)*17.03/mugaz_c
               print*,'WARNING, cpp for NH3 may be for liquid'
            elseif(igas.eq.igas_C2H6)then
               ! C2H6
               ! http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=28
               cpp_c   = cpp_c   + 1.763*gfrac(igas)*30.07/mugaz_c
            elseif(igas.eq.igas_C2H2)then
               ! C2H2
               ! http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=1
               cpp_c   = cpp_c   + 1.575*gfrac(igas)*26.04/mugaz_c
            else
               print*,'Error in calc_cpp_mugaz: Gas species not recognised!'
               call abort
            endif
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
