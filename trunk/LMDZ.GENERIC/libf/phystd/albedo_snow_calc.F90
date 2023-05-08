      subroutine albedo_snow_calc(lambda,albedo,A_sw)

      real lambda_sw
      real lambda_lw
      real A_sw
      real A_lw
      real albedo
      real lambda
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! Routine written by Martin TURBET (July 2015)
      !!!!!
      !!!!! Albedo expression derived from Yoshi et al. and Warren et al. articles.
      !!!!!
      !!
      !! Short Wave Snow Albedo 'A_sw' (at 0.5 microns) must be chosen in callphys.def file.
      !! 
      !!!!!
      !!!!!
      !! Some Calibration Informations :
      !!!!!
      !!!!!
      !!
      !! 1. For Pure Snow, A_sw (at 0.5 microns) should be equal to ~0.95
      !!    -> This gives an equivalent integrated snow albedo of ~0.73 for the Sun.
      !!
      !! 2. For Dusty Snow, A_sw (at 0.5 microns) should be equal to ~0.50.
      !!    -> This gives an equivalent integrated snow albedo of ~0.39 for the Sun.
      !!
      !! 3. For 'Realistic' Snow, A_sw (at 0.5 microns) should be equal to ~0.645.
      !!    -> This gives an equivalent integrated snow albedo of ~0.50 for the Sun.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      A_lw=5.0D-2
      lambda_sw=9.0D-1
      lambda_lw=1.4
      
      if (lambda .le. lambda_sw) then
         albedo=A_sw
      else if (lambda .ge. lambda_lw) then
         albedo=A_lw
      else if ( (lambda .gt. lambda_sw) .and. &
             (lambda .lt. lambda_lw) ) then
         albedo=A_sw-(lambda-lambda_sw)*(A_sw-A_lw) &
           /(lambda_lw-lambda_sw)  
      end if

      return
      end
