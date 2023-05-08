     subroutine interpolateH2Ocont_CKD(wn,temp,presS,presF,abcoef,firstcall,ind)

!==================================================================
!     
!     Purpose
!     -------
!     Calculates the H2O continuum opacity, using a lookup table from
!     Clough (2005).
!
!     Authors
!     -------
!     R. Wordsworth (2011)
!     
!==================================================================

      use datafile_mod, only: datadir
      implicit none

      ! input
      double precision wn                 ! wavenumber             (cm^-1)
      double precision temp               ! temperature            (Kelvin)
      double precision presS              ! self-pressure          (Pascals)
      double precision presF              ! foreign (air) pressure (Pascals)

      ! output
      double precision abcoef             ! absorption coefficient (m^-1)

      integer nS,nT
      parameter(nS=1001)
      parameter(nT=11)

      double precision kB
      parameter(kB=1.3806488e-23)

      double precision amagatS, amagatF, abcoefS, abcoefF, Nmolec
      double precision wn_arr(nS)
      double precision temp_arr(nT)
      double precision abs_arrS(nS,nT)
      double precision abs_arrF(nS,nT)
      double precision data_tmp(nT)

      integer k,ind
      logical firstcall

      save wn_arr, temp_arr, abs_arrS, abs_arrF !read by master

      character*100 dt_file
      integer strlen,ios

      amagatS=(273.15/temp)*(presS/101325.0)
      amagatF=(273.15/temp)*(presF/101325.0)

      if(firstcall)then ! called by sugas_corrk only
         print*,'----------------------------------------------------'
         print*,'Initialising H2O continuum from MT_CKD data...'

!     1.1 Open the ASCII files

!$OMP MASTER
         ! nu array
         dt_file=TRIM(datadir)//'/continuum_data/H2O_CONT_NU.dat'
         open(33,file=dt_file,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from interpolateH2O_cont'
           write(*,*) 'Data file ',trim(dt_file),' not found.'
           write(*,*)'Check that your path to datagcm:',trim(datadir)
           write(*,*)' is correct. You can change it in callphys.def with:'
           write(*,*)' datadir = /absolute/path/to/datagcm'
           write(*,*)' Also check that there is a continuum_data/H2O_CONT_NU.dat there.'
           call abort
         else
            do k=1,nS
               read(33,*) wn_arr(k)
            enddo
         endif
         close(33)

         ! self broadening
         dt_file=TRIM(datadir)//'/continuum_data/H2O_CONT_SELF.dat'
         open(34,file=dt_file,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from interpolateH2O_cont'
           write(*,*) 'Data file ',trim(dt_file),' not found.'
           write(*,*)'Check that your path to datagcm:',trim(datadir)
           write(*,*)' is correct. You can change it in callphys.def with:'
           write(*,*)' datadir = /absolute/path/to/datagcm'
           write(*,*)' Also check that there is a continuum_data/H2O_CONT_SELF.dat there.'
           call abort
         else
            do k=1,nS
               read(34,*) data_tmp
               abs_arrS(k,1:nT)=data_tmp(1:nT)
            end do
         endif
         close(34)

         ! foreign (N2+O2+Ar) broadening
         dt_file=TRIM(datadir)//'/continuum_data/H2O_CONT_FOREIGN.dat'
         open(35,file=dt_file,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from interpolateH2O_cont'
           write(*,*) 'Data file ',trim(dt_file),' not found.'
           write(*,*)'Check that your path to datagcm:',trim(datadir)
           write(*,*)' is correct. You can change it in callphys.def with:'
           write(*,*)' datadir = /absolute/path/to/datagcm'
           write(*,*)' Also check that there is a continuum_data/H2O_CONT_FOREIGN.dat there.'
           call abort
         else
            do k=1,nS
               read(35,*) data_tmp
               abs_arrF(k,1:nT)=data_tmp(1:nT)
            end do
         endif
         close(35)

         temp_arr(1)  = 200.
         temp_arr(2)  = 250.
         temp_arr(3)  = 300.
         temp_arr(4)  = 350.
         temp_arr(5)  = 400.
         temp_arr(6)  = 450.
         temp_arr(7)  = 500.
         temp_arr(8)  = 550.
         temp_arr(9)  = 600.
         temp_arr(10) = 650.
         temp_arr(11) = 700.

         print*,'interpolateH2Ocont: At wavenumber ',wn,' cm^-1'
         print*,'   temperature ',temp,' K'
         print*,'   H2O pressure ',presS,' Pa'
         print*,'   air pressure ',presF,' Pa'
!$OMP END MASTER
!$OMP BARRIER
	 
      endif

      call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arrS,wn,temp,abcoefS,ind)
!      print*,'the self absorption is ',abcoefS,' cm^2 molecule^-1'

      call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arrF,wn,temp,abcoefF,ind)
!      print*,'the foreign absorption is ',abcoefF,' cm^2 molecule^-1'

!      print*,'We have ',amagatS,' amagats of H2O vapour'
!      print*,'and ',amagatF,' amagats of air'

      abcoef = abcoefS*amagatS + abcoefF*amagatF ! Eq. (15) in Clough (1989)
      abcoef = abcoef*(presS/(presF+presS))      ! take H2O mixing ratio into account
                                                    ! abs coeffs are given per molecule of H2O

      Nmolec = (presS+presF)/(kB*temp)           ! assume ideal gas
!      print*,'Total number of molecules per m^3 is',Nmolec

      abcoef = abcoef*Nmolec/(100.0**2)          ! convert to m^-1
!      print*,'So the total absorption is ',abcoef,' m^-1'
!      print*,'And optical depth / km : ',1000.0*abcoef


      if(wn.gt.500 .and. wn.lt.1400)then
      elseif(wn.gt.2100 .and. wn.lt.3000)then
      else
         abcoef = 0.0
      endif

      ! unlike for Rayleigh scattering, we do not currently weight by the BB function
      ! however our bands are normally thin, so this is no big deal.


      return
    end subroutine interpolateH2Ocont_CKD

