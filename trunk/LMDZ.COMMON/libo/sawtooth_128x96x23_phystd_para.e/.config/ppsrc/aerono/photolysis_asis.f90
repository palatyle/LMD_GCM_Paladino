










!==========================================================================

      subroutine photolysis_asis(nlayer, ngrid,                                 &
                                 lswitch, press, temp, sza, fractcol, tauref,   &
                                 zmmean, dist_sol, rmco2, rmo3, v_phot)

!==========================================================================

      use comcstfi_mod
      use callkeys_mod

      implicit none

!--------------------------------------------
!     data for photochemistry
!--------------------------------------------

!--------------------------------------------
!     dimensions of photolysis lookup table
!--------------------------------------------

      integer, parameter :: nd    = 13  ! species
      integer, parameter :: nz    = 143 ! altitude
      integer, parameter :: nozo  = 7   ! ozone
      integer, parameter :: nsza  = 27  ! solar zenith angle
      integer, parameter :: ntemp = 4   ! temperature
      integer, parameter :: ntau  = 8   ! dust

!--------------------------------------------

      real, parameter :: kb = 1.3806e-23

      common/chimiedata/jphot,colairtab,table_ozo

      real jphot(ntemp,nsza,nz,nozo,ntau,nd)
      real colairtab(nz)
      real szatab(nsza)
      real table_ozo(nozo)
      real tautab(ntau)

      data szatab/0.,  5., 10., 15., 20., 25.,                          &
     &            30., 35., 40., 45., 50., 55.,                         &
     &            60., 65., 70., 75., 80., 82.,                         &
     &            84., 86., 88., 90., 91., 92.,                         &
     &            93., 94., 95./

      data tautab/0., 0.2, 0.4, 0.6, 0.8, 1., 2., 4./

!--------------------------------------------
!     number of reactions in ASIS solver
!--------------------------------------------

      integer, parameter :: nb_phot_max       = 18
      integer, parameter :: nb_reaction_3_max = 6
      integer, parameter :: nb_reaction_4_max = 30

!==========================================================================
!     input:
!==========================================================================
       
      integer, intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: ngrid   ! number of atmospheric columns
      integer :: lswitch            ! interface level between chemistries
      real :: press(nlayer)         ! pressure (hPa)
      real :: temp(nlayer)          ! temperature (K)
      real :: sza                   ! solar zenith angle (deg)
      real :: fractcol              ! day fraction
      real :: tauref                ! optical depth at 7 hpa
      real :: zmmean(nlayer)        ! mean molecular mass (g)
      real :: dist_sol              ! sun distance (AU) 
      real :: rmco2(nlayer)         ! co2 mixing ratio
      real :: rmo3(nlayer)          ! ozone mixing ratio

!==========================================================================
!     output: interpolated photodissociation rates (s-1)
!==========================================================================

      real (kind = 8), dimension(nlayer,nb_phot_max) :: v_phot

!==========================================================================
!     local:
!==========================================================================

      integer :: icol, ij, indsza, indtau, indcol, indozo, indtemp,     &
                 iozo, isza, itau, it, l

      integer :: j_o2_o, j_o2_o1d, j_co2_o, j_co2_o1d, j_o3_o1d,        &
                 j_o3_o, j_h2o, j_hdo, j_h2o2, j_ho2, j_no, j_no2,      &
                 j_hno3, j_hno4,                                        &
                 j_ch4_ch3_h, j_ch4_1ch2_h2, j_ch4_3ch2_h_h,            &
                 j_ch4_ch_h2_h, j_ch3o2h, j_ch2o_hco, j_ch2o_co,        &
                 j_ch3oh, j_c2h6, j_hcl, j_hocl, j_clo, j_so2, j_so,    &
                 j_h2s, j_so3

      real :: col(nlayer)                 ! overhead air column   (molecule cm-2)
      real :: colo3(nlayer)               ! overhead ozone column (molecule cm-2)
      real :: poids(2,2,2,2,2)            ! 5D interpolation weights 
      real :: tref                        ! temperature  at 1.9 hPa in the gcm (K)
      real :: table_temp(ntemp)           ! temperatures at 1.9 hPa in jmars   (K)
      real :: cinf, csup, cicol, ciozo, cisza, citemp, citau
      real :: colo3min, dp, coef
      real :: ratio_o3(nlayer)
      real :: tau
      real :: j(nlayer,nd)

!==========================================================================
!     day/night criterion
!==========================================================================

      if (sza <= 95.) then 

!==========================================================================
!     temperatures at 1.9 hPa in lookup table
!==========================================================================
      
      table_temp(1) = 226.2
      table_temp(2) = 206.2
      table_temp(3) = 186.2
      table_temp(4) = 169.8

!==========================================================================
!     interpolation in solar zenith angle
!==========================================================================
 
      indsza = nsza - 1
      do isza = 1,nsza
         if (szatab(isza) >= sza) then
            indsza = isza - 1
            indsza = max(indsza, 1)
            exit
         end if
      end do
      cisza = (sza - szatab(indsza))  &
             /(szatab(indsza + 1) - szatab(indsza))

!==========================================================================
!     interpolation in dust (tau)
!==========================================================================

      tau = min(tauref, tautab(ntau))
      tau = max(tau, tautab(1))

      indtau = ntau - 1
      do itau = 1,ntau
         if (tautab(itau) >= tau) then
            indtau = itau - 1
            indtau = max(indtau, 1)
            exit
         end if
      end do
      citau = (tau - tautab(indtau))     &
             /(tautab(indtau + 1) - tautab(indtau))

!==========================================================================
!     co2 and ozone columns
!==========================================================================

!     co2 column at model top (molecule.cm-2)

      col(lswitch-1) = 6.022e22*rmco2(lswitch-1)*press(lswitch-1)*100.  &
                       /(zmmean(lswitch-1)*g)

!     ozone column at model top

      colo3(lswitch-1) = 0.

!     co2 and ozone columns for other levels (molecule.cm-2)

      do l = lswitch-2,1,-1
         dp = (press(l) - press(l+1))*100.
         col(l) = col(l+1) + (rmco2(l+1) + rmco2(l))*0.5   &
                             *6.022e22*dp/(zmmean(l)*g)
         col(l) = min(col(l), colairtab(1))
         colo3(l) = colo3(l+1) + (rmo3(l+1) + rmo3(l))*0.5 &
                                 *6.022e22*dp/(zmmean(l)*g)
      end do

!     ratio ozone column/minimal theoretical column (0.1 micron-atm) 

!     ro3 = 7.171e-10 is the o3 mixing ratio for a uniform
!     profile giving a column 0.1 micron-atmosphere at
!     a surface pressure of 10 hpa.

      do l = 1,lswitch-1
         colo3min    = col(l)*7.171e-10
         ratio_o3(l) = colo3(l)/colo3min
         ratio_o3(l) = min(ratio_o3(l), table_ozo(nozo)*10.)
         ratio_o3(l) = max(ratio_o3(l), 1.)
      end do

!==========================================================================
!     temperature dependence
!==========================================================================

!     1) search for temperature at 1.9 hPa (tref): vertical interpolation

      tref = temp(1)
      do l = (lswitch-1)-1,1,-1
         if (press(l) > 1.9) then
            cinf = (press(l) - 1.9)/(press(l) - press(l+1))
            csup = 1. - cinf
            tref = cinf*temp(l+1) + csup*temp(l)
            exit
         end if
      end do 

!     2) interpolation in temperature

      tref = min(tref,table_temp(1))
      tref = max(tref,table_temp(ntemp))

      do it = 2, ntemp
         if (table_temp(it) <= tref) then
            citemp = (log(tref) - log(table_temp(it)))              &
                    /(log(table_temp(it-1)) - log(table_temp(it)))
            indtemp = it - 1
            exit
         end if
      end do

!==========================================================================
!     loop over vertical levels
!==========================================================================

      do l = 1,lswitch-1

!     interpolation in air column

         indcol = nz - 1
         do icol = 1,nz
            if (colairtab(icol) < col(l)) then
               indcol = icol - 1
               exit
            end if
         end do
         cicol = (log(col(l)) - log(colairtab(indcol + 1)))              &
                /(log(colairtab(indcol)) - log(colairtab(indcol + 1)))

!     interpolation in ozone column

         indozo = nozo - 1
         do iozo = 1,nozo
            if (table_ozo(iozo)*10. >= ratio_o3(l)) then
               indozo = iozo - 1
               indozo = max(indozo, 1)
               exit
            end if
         end do
         ciozo = (ratio_o3(l) - table_ozo(indozo)*10.)             &
                /(table_ozo(indozo + 1)*10. - table_ozo(indozo)*10.)

!     4-dimensional interpolation weights

!     poids(temp,sza,co2,o3,tau)

         poids(1,1,1,1,1) = citemp*(1.-cisza)*cicol*(1.-ciozo)*(1.-citau)
         poids(1,1,1,2,1) = citemp*(1.-cisza)*cicol*ciozo*(1.-citau)
         poids(1,1,2,1,1) = citemp*(1.-cisza)*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(1,1,2,2,1) = citemp*(1.-cisza)*(1.-cicol)*ciozo*(1.-citau)
         poids(1,2,1,1,1) = citemp*cisza*cicol*(1.-ciozo)*(1.-citau)
         poids(1,2,1,2,1) = citemp*cisza*cicol*ciozo*(1.-citau)
         poids(1,2,2,1,1) = citemp*cisza*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(1,2,2,2,1) = citemp*cisza*(1.-cicol)*ciozo*(1.-citau)
         poids(2,1,1,1,1) = (1.-citemp)*(1.-cisza)*cicol*(1.-ciozo)*(1.-citau)
         poids(2,1,1,2,1) = (1.-citemp)*(1.-cisza)*cicol*ciozo*(1.-citau)
         poids(2,1,2,1,1) = (1.-citemp)*(1.-cisza)*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(2,1,2,2,1) = (1.-citemp)*(1.-cisza)*(1.-cicol)*ciozo*(1.-citau)
         poids(2,2,1,1,1) = (1.-citemp)*cisza*cicol*(1.-ciozo)*(1.-citau)
         poids(2,2,1,2,1) = (1.-citemp)*cisza*cicol*ciozo*(1.-citau)
         poids(2,2,2,1,1) = (1.-citemp)*cisza*(1.-cicol)*(1.-ciozo)*(1.-citau)
         poids(2,2,2,2,1) = (1.-citemp)*cisza*(1.-cicol)*ciozo*(1.-citau)
!
         poids(1,1,1,1,2) = citemp*(1.-cisza)*cicol*(1.-ciozo)*citau
         poids(1,1,1,2,2) = citemp*(1.-cisza)*cicol*ciozo*citau
         poids(1,1,2,1,2) = citemp*(1.-cisza)*(1.-cicol)*(1.-ciozo)*citau
         poids(1,1,2,2,2) = citemp*(1.-cisza)*(1.-cicol)*ciozo*citau
         poids(1,2,1,1,2) = citemp*cisza*cicol*(1.-ciozo)*citau
         poids(1,2,1,2,2) = citemp*cisza*cicol*ciozo*citau
         poids(1,2,2,1,2) = citemp*cisza*(1.-cicol)*(1.-ciozo)*citau
         poids(1,2,2,2,2) = citemp*cisza*(1.-cicol)*ciozo*citau
         poids(2,1,1,1,2) = (1.-citemp)*(1.-cisza)*cicol*(1.-ciozo)*citau
         poids(2,1,1,2,2) = (1.-citemp)*(1.-cisza)*cicol*ciozo*citau
         poids(2,1,2,1,2) = (1.-citemp)*(1.-cisza)*(1.-cicol)*(1.-ciozo)*citau
         poids(2,1,2,2,2) = (1.-citemp)*(1.-cisza)*(1.-cicol)*ciozo*citau
         poids(2,2,1,1,2) = (1.-citemp)*cisza*cicol*(1.-ciozo)*citau
         poids(2,2,1,2,2) = (1.-citemp)*cisza*cicol*ciozo*citau
         poids(2,2,2,1,2) = (1.-citemp)*cisza*(1.-cicol)*(1.-ciozo)*citau
         poids(2,2,2,2,2) = (1.-citemp)*cisza*(1.-cicol)*ciozo*citau

!     4-dimensional interpolation in the lookup table

         do ij = 1,nd
            j(l,ij) =                                                                & 
            poids(1,1,1,1,1)*jphot(indtemp,indsza,indcol,indozo,indtau,ij)           &
          + poids(1,1,1,2,1)*jphot(indtemp,indsza,indcol,indozo+1,indtau,ij)         &
          + poids(1,1,2,1,1)*jphot(indtemp,indsza,indcol+1,indozo,indtau,ij)         &
          + poids(1,1,2,2,1)*jphot(indtemp,indsza,indcol+1,indozo+1,indtau,ij)       &
          + poids(1,2,1,1,1)*jphot(indtemp,indsza+1,indcol,indozo,indtau,ij)         &
          + poids(1,2,1,2,1)*jphot(indtemp,indsza+1,indcol,indozo+1,indtau,ij)       &
          + poids(1,2,2,1,1)*jphot(indtemp,indsza+1,indcol+1,indozo,indtau,ij)       &
          + poids(1,2,2,2,1)*jphot(indtemp,indsza+1,indcol+1,indozo+1,indtau,ij)     &
          + poids(2,1,1,1,1)*jphot(indtemp+1,indsza,indcol,indozo,indtau,ij)         &
          + poids(2,1,1,2,1)*jphot(indtemp+1,indsza,indcol,indozo+1,indtau,ij)       &
          + poids(2,1,2,1,1)*jphot(indtemp+1,indsza,indcol+1,indozo,indtau,ij)       &
          + poids(2,1,2,2,1)*jphot(indtemp+1,indsza,indcol+1,indozo+1,indtau,ij)     &
          + poids(2,2,1,1,1)*jphot(indtemp+1,indsza+1,indcol,indozo,indtau,ij)       &
          + poids(2,2,1,2,1)*jphot(indtemp+1,indsza+1,indcol,indozo+1,indtau,ij)     &
          + poids(2,2,2,1,1)*jphot(indtemp+1,indsza+1,indcol+1,indozo,indtau,ij)     &
          + poids(2,2,2,2,1)*jphot(indtemp+1,indsza+1,indcol+1,indozo+1,indtau,ij)   &
!
          + poids(1,1,1,1,2)*jphot(indtemp,indsza,indcol,indozo,indtau+1,ij)         &
          + poids(1,1,1,2,2)*jphot(indtemp,indsza,indcol,indozo+1,indtau+1,ij)       &
          + poids(1,1,2,1,2)*jphot(indtemp,indsza,indcol+1,indozo,indtau+1,ij)       &
          + poids(1,1,2,2,2)*jphot(indtemp,indsza,indcol+1,indozo+1,indtau+1,ij)     &
          + poids(1,2,1,1,2)*jphot(indtemp,indsza+1,indcol,indozo,indtau+1,ij)       &
          + poids(1,2,1,2,2)*jphot(indtemp,indsza+1,indcol,indozo+1,indtau+1,ij)     &
          + poids(1,2,2,1,2)*jphot(indtemp,indsza+1,indcol+1,indozo,indtau+1,ij)     &
          + poids(1,2,2,2,2)*jphot(indtemp,indsza+1,indcol+1,indozo+1,indtau+1,ij)   &
          + poids(2,1,1,1,2)*jphot(indtemp+1,indsza,indcol,indozo,indtau+1,ij)       &
          + poids(2,1,1,2,2)*jphot(indtemp+1,indsza,indcol,indozo+1,indtau+1,ij)     &
          + poids(2,1,2,1,2)*jphot(indtemp+1,indsza,indcol+1,indozo,indtau+1,ij)     &
          + poids(2,1,2,2,2)*jphot(indtemp+1,indsza,indcol+1,indozo+1,indtau+1,ij)   &
          + poids(2,2,1,1,2)*jphot(indtemp+1,indsza+1,indcol,indozo,indtau+1,ij)     &
          + poids(2,2,1,2,2)*jphot(indtemp+1,indsza+1,indcol,indozo+1,indtau+1,ij)   &
          + poids(2,2,2,1,2)*jphot(indtemp+1,indsza+1,indcol+1,indozo,indtau+1,ij)   &
          + poids(2,2,2,2,2)*jphot(indtemp+1,indsza+1,indcol+1,indozo+1,indtau+1,ij)
         end do

!     correction for sun distance

         do ij = 1,nd
!            j(l,ij) = j(l,ij)*(1.52/dist_sol)**2.
            j(l,ij) = j(l,ij)*(1.0/dist_sol)**2.

            ! Only during daylight.
            if((ngrid.eq.1))then
                  j(l,ij)= j(l,ij)* 0.25 ! globally averaged = divide by 4
            elseif(diurnal .eqv. .false.) then
                  j(l,ij)= j(l,ij) * fractcol
            endif
         end do 

!==========================================================================
!     end of loop over vertical levels
!==========================================================================

      end do

      else

!==========================================================================
!     night
!==========================================================================

         j(:,:) = 0.

      end if

! photodissociation rates numbering in the lookup table

! jmars.20140930


      j_o2_o         =  1      ! o2 + hv     -> o + o
      j_o2_o1d       =  2      ! o2 + hv     -> o + o(1d)
      j_co2_o        =  3      ! co2 + hv    -> co + o
      j_co2_o1d      =  4      ! co2 + hv    -> co + o(1d)
      j_o3_o1d       =  5      ! o3 + hv     -> o2 + o(1d)
      j_o3_o         =  6      ! o3 + hv     -> o2 + o
      j_h2o          =  7      ! h2o + hv    -> h + oh
      j_h2o2         =  8      ! h2o2 + hv   -> oh + oh
      j_ho2          =  9      ! ho2 + hv    -> oh + o
      j_no           =  10     ! no + hv     -> n + o
      j_no2          =  11     ! no2 + hv    -> no + o
      j_hno3         =  12     ! hno3 + hv   -> no2 + oh
      j_hno4         =  13     ! hno4 + hv   -> no2 + ho2

! jmars.20111014

!     j_o2_o         =  1      ! o2 + hv     -> o + o
!     j_o2_o1d       =  2      ! o2 + hv     -> o + o(1d)
!     j_co2_o        =  3      ! co2 + hv    -> co + o
!     j_co2_o1d      =  4      ! co2 + hv    -> co + o(1d)
!     j_o3_o1d       =  5      ! o3 + hv     -> o2 + o(1d)
!     j_o3_o         =  6      ! o3 + hv     -> o2 + o
!     j_h2o          =  7      ! h2o + hv    -> h + oh
!     j_hdo          =  8      ! hdo + hv    -> d + oh
!     j_h2o2         =  9      ! h2o2 + hv   -> oh + oh
!     j_ho2          =  10     ! ho2 + hv    -> oh + o
!     j_no2          =  11     ! no2 + hv    -> no + o
!     j_ch4_ch3_h    =  12     ! ch4 + hv    -> ch3 + h
!     j_ch4_1ch2_h2  =  13     ! ch4 + hv    -> 1ch2 + h2
!     j_ch4_3ch2_h_h =  14     ! ch4 + hv    -> 3ch2 + h + h
!     j_ch4_ch_h2_h  =  15     ! ch4 + hv    -> ch + h2 + h
!     j_ch3o2h       =  16     ! ch3o2h + hv -> ch3o + oh
!     j_ch2o_hco     =  17     ! ch2o + hv   -> h + hco
!     j_ch2o_co      =  18     ! ch2o + hv   -> h2 + co
!     j_ch3oh        =  19     ! ch3oh + hv  -> ch3o + h
!     j_c2h6         =  20     ! c2h6 + hv   -> products
!     j_hcl          =  21     ! hcl + hv    -> h + cl
!     j_hocl         =  22     ! hocl + hv   -> oh + cl
!     j_clo          =  23     ! clo + hv    -> o + cl
!     j_so2          =  24     ! so2 + hv    -> so + o
!     j_so           =  25     ! so + hv     -> s + o
!     j_h2s          =  26     ! h2s + hv    -> hs + s
!     j_so3          =  27     ! so2 + hv    -> so2 + o
!     j_hno3         =  28     ! hno3 + hv   -> oh + no2
!     j_hno4         =  29     ! hno4 + hv   -> ho2 + no2

! fill v_phot array

      v_phot(:,:) = 0.

      do l = 1,lswitch-1
         v_phot(l, 1) = j(l,j_o2_o)
         v_phot(l, 2) = j(l,j_o2_o1d)
         v_phot(l, 3) = j(l,j_co2_o)
         v_phot(l, 4) = j(l,j_co2_o1d)
         v_phot(l, 5) = j(l,j_o3_o1d)
         v_phot(l, 6) = j(l,j_o3_o)
         v_phot(l, 7) = j(l,j_h2o)
         v_phot(l, 8) = j(l,j_h2o2)
         v_phot(l, 9) = j(l,j_ho2)
         v_phot(l,10) = j(l,j_no)
         v_phot(l,11) = j(l,j_no2)
      end do

      return
      end

!*****************************************************************
