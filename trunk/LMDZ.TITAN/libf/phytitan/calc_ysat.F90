SUBROUTINE calc_ysat(nlon,nlay,press,temp,ysat)

  !     ==============================================================================
  !     Purpose
  !     -------
  !     Compute saturation profiles (mol/mol) for chemical tracers
  !     
  !     Authors
  !     -------
  !     S. Lebonnois 
  !     -> inicondens.F in old Titan, with T,P in planetary average
  !     J. Vatant d'Ollone (2017)
  !     -> Adapt to new physics, move to F90 and compute every grid point 
  !     ==============================================================================


  !-----------------------------------------------------------------------
  !     Declarations:
  !     -------------

  USE comchem_h, only: nkim, cnames

  IMPLICIT NONE

  !     Arguments :
  !     -----------
  INTEGER, INTENT(IN)                           :: nlon  ! # of grid points in the chunk
  INTEGER, INTENT(IN)                           :: nlay  ! # of vertical layes 
  REAL, DIMENSION(nlon,nlay),      INTENT(IN)   :: press ! Mid-layers pressure (Pa)
  REAL, DIMENSION(nlon,nlay),      INTENT(IN)   :: temp  ! Mid-layers temperature (K)
  REAL, DIMENSION(nlon,nlay,nkim), INTENT(OUT)  :: ysat  ! Saturation profiles (mol/mol)

  !     Local variables :
  !     -----------------
  INTEGER                           ::  ic
  REAL,DIMENSION(:,:), ALLOCATABLE  ::  x
  !     -------------------------------------------------------------------------

  ALLOCATE(x(nlon,nlay))

  !     By default, ysat=1, meaning we do not condense
  ysat(:,:,:) = 1.0

  !     Main loop

  do ic=1,nkim

     if(trim(cnames(ic)).eq."CH4") then

        where (temp(:,:).lt.90.65)                                      
           ysat(:,:,ic) = 10.0**(4.42507e0 - ( ( ( 1165560.7e0 / temp(:,:)                  &
                -  115352.19e0 ) / temp(:,:) + 4055.6016e0 ) / temp(:,:)                    &
                + 453.92414e0 ) / temp(:,:) ) / press(:,:) * 1013.25e0
        elsewhere
           ysat(:,:,ic) = 10.0**(3.901408e0 - ( ( 154567.02e0 / temp(:,:) - 1598.8512e0 )   &
                / temp(:,:) + 437.54809e0 ) / temp(:,:)) / press(:,:) * 1013.25e0
        endwhere

        ! Forcing CH4 to 1.4% minimum               
        where (ysat(:,:,ic).lt.0.014) ysat(:,:,ic) = 0.014

     else if(trim(cnames(ic)).eq."C2H2") then

        ysat(:,:,ic) = 10.0**(6.09748e0-1644.1e0/temp(:,:)+7.42346e0                        &
             * alog10(1.0e3/temp(:,:)) ) / press(:,:)*1013.25e0/760.0

     else if(trim(cnames(ic)).eq."C2H4") then

        where (temp(:,:).lt.89.0)
           ysat(:,:,ic) = 10.0**(1.5477e0 + (1.0e0/temp(:,:) - 0.011e0)                     &
                * (16537.0e0*(1.0e0/temp(:,:) - 0.011e0) - 1038.1e0))                       &
                / press(:,:) * 1.01325e0 / 760.0
        elsewhere (temp(:,:).lt.104.0)
           ysat(:,:,ic) = 10.0**(8.724e0 - 901.6e0/(temp(:,:) - 2.555e0) )                  &
                / press(:,:) * 1013.25e0 / 760.0
        elsewhere (temp(:,:).lt.120.0)
           ysat(:,:,ic) = 10.0**(50.79e0 - 1703.0e0/temp(:,:) - 17.141e0                    &
                * alog10(temp(:,:)) ) / press(:,:) * 1013.25e0 / 760.0
        elsewhere (temp(:,:).lt.155.0)
           ysat(:,:,ic) = 10.0**(6.74756e0 - 585.0e0/(temp(:,:) - 18.16e0) )                &
                / press(:,:) * 1013.25e0 / 760.0
        endwhere

     else if(trim(cnames(ic)).eq."C2H6") then

        where (temp(:,:).lt.90.)
           ysat(:,:,ic) = 10.0**(10.01e0-1085.0e0/(temp(:,:)-0.561e0) )                     &
                / press(:,:) * 1013.25e0 / 760.0e0
        elsewhere
           ysat(:,:,ic) = 10.0**(5.9366e0 - 1086.17e0/temp(:,:) + 3.83464e0                 &
                * alog10(1.0e3/temp(:,:)) ) / press(:,:)*1013.25e0/760.0
        endwhere

     else if((trim(cnames(ic)).eq."CH3CCH") .or. (trim(cnames(ic)).eq."CH2CCH2")) then

        ysat(:,:,ic) = 10.0**(2.8808e0 - 4.5e0*(249.9e0 - temp(:,:))                        &
             / (1.15e0*temp(:,:) - 37.485e0) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C3H6")  then

        ysat(:,:,ic) = 10.0**(7.4463e0 - 1028.5654e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C3H8")  then

        ysat(:,:,ic) = 10.0**(7.217e0 - 994.30251e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if((trim(cnames(ic)).eq."C4H2").or.(trim(cnames(ic)).eq."C4H2s")) then

        ysat(:,:,ic) = 10.0**(96.26781e0 - 4651.872e0/temp(:,:) - 31.68595e0                &
             * alog10(temp(:,:)) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C4H4")  then

        ysat(:,:,ic)= 1.0e3 * exp(9.3898e0 - 2203.57/(temp(:,:)-43.15e0) ) / press(:,:)

     else if(trim(cnames(ic)).eq."C4H6")  then

        ysat(:,:,ic)= 10.0**(2.8808e0 - 4.6e0*(262.3e0 - temp(:,:))                         &
             /(1.15e0*temp(:,:) - 39.345e0) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C4H10")  then

        ysat(:,:,ic)= 10.0**(8.446e0 - 1461.2e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C6H2")  then

        ysat(:,:,ic)= 10.0**(4.666e0 - 4956e0/temp(:,:) + 25.845e0 *                        &
             alog10(1.0e3/temp(:,:)) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C8H2")  then

        ysat(:,:,ic)= 10.0**(3.95e0 - 6613e0/temp(:,:) + 35.055e0 *                         &
             alog10(1.0e3/temp(:,:)) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."AC6H6")  then

        x = 1.0e0 - temp(:,:) / 562.2e0
        ysat(:,:,ic)= 48.9e3 * exp( ( 1.33213 * x**1.5 - 6.98273 * x                        &
             - x**3 * (2.62863 + 3.33399 * x**3) ) * 562.2e0/temp(:,:) ) / press(:,:)

     else if(trim(cnames(ic)).eq."HCN")  then

        ysat(:,:,ic)= 10.0**(8.6165e0 - 1516.5e0/(temp(:,:) - 26.2e0) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."CH3CN")  then

        ysat(:,:,ic)= 10.0**(8.458e0 - 1911.7e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C2H3CN")  then

        ysat(:,:,ic)= 10.0**(9.3051e0 - 2782.21/(temp(:,:) - 51.15e0) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."NCCN")  then

        ysat(:,:,ic)=  10.0**(7.454e0 - 1832e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."HC3N")  then

        ysat(:,:,ic)= 10.0**(7.7446e0 - 1453.5609e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     else if(trim(cnames(ic)).eq."C4N2")  then

        ysat(:,:,ic)= 10.0**(8.269e0 - 2155.0e0/temp(:,:) ) / press(:,:) * 1013.25e0 / 760.0e0

     endif
  enddo

END SUBROUTINE calc_ysat
