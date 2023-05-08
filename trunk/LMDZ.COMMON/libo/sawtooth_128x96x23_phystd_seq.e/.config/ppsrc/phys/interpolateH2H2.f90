










     subroutine interpolateH2H2(wn,temp,pres,abcoef,firstcall,ind)

!==================================================================
!     
!     Purpose
!     -------
!     Calculates the H2-H2 CIA opacity, using a lookup table from
!     HITRAN (2011 or later)
!
!     Authors
!     -------
!     R. Wordsworth (2011)
!
!     + J.Vatant d'Ollone (2019) 
!        - Enable updated HITRAN file (Karman2019,Fletcher2018)
!           (2018 one should be default for giant planets)
!==================================================================

      use callkeys_mod, only: versH2H2cia
      use datafile_mod, only: datadir

      implicit none

      ! input
      double precision wn                 ! wavenumber             (cm^-1)
      double precision temp               ! temperature            (Kelvin)
      double precision pres               ! pressure               (Pascals)

      ! output
      double precision abcoef             ! absorption coefficient (m^-1)

      integer nS,nT
      parameter(nT=10)

      double precision, parameter :: losch = 2.6867774e19
      ! Loschmit's number (molecule cm^-3 at STP) 
      ! converts cm^5 molecule^-2 --> cm^-1 amagat^-2
      ! see Richard et al. 2011, JQSRT for details

      double precision amagat
      double precision temp_arr(nT)
      
      double precision, dimension(:),   allocatable :: wn_arr
      double precision, dimension(:,:), allocatable :: abs_arr

      integer k,iT
      logical firstcall

      save nS, wn_arr, temp_arr, abs_arr !read by master

      character*100 dt_file
      integer ios

      character(LEN=*), parameter :: fmat11 = "(A20,F10.3,F10.3,I7,F7.1,E10.3,F5.3)"
      character(LEN=*), parameter :: fmat18 = "(A12,A3,A5,F10.6,F10.4,I7,F7.3,E10.3,F5.3)"

      character*20 bleh
      double precision blah, Ttemp
      integer nres

      integer ind

      if(temp.gt.400)then
         print*,'Your temperatures are too high for this H2-H2 CIA dataset. If you '
         print*,'really want to run simulations with hydrogen at T > 400 K, contact'
         print*,'Robin Wordsworth [rwordsworth@uchicago.edu].'
         stop
      endif

      amagat = (273.15/temp)*(pres/101325.0)

      if(firstcall)then ! called by sugas_corrk only
         print*,'----------------------------------------------------'
         print*,'Initialising H2-H2 continuum from HITRAN database...'

!     1.1 Open the ASCII files and set nS according to version
         ! Only two possible versions for now : 2011 or 2018 (sanity check in inifis_mod)
         if (versH2H2cia.eq.2011) then
           dt_file=TRIM(datadir)//'/continuum_data/H2-H2_norm_2011.cia'
           nS = 2428
         else if (versH2H2cia.eq.2018) then
           dt_file=TRIM(datadir)//'/continuum_data/H2-H2_norm_2018.cia'
           nS = 9600
         endif

         if(.not.allocated(wn_arr))  allocate(wn_arr(nS)) 
         if(.not.allocated(abs_arr)) allocate(abs_arr(nS,nT)) 

!$OMP MASTER
         open(33,file=dt_file,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from interpolateH2H2'
           write(*,*) 'Data file ',trim(dt_file),' not found.'
           write(*,*) 'Check that your path to datagcm:',trim(datadir)
           write(*,*) 'is correct. You can change it in callphys.def with:'
           write(*,*) 'datadir = /absolute/path/to/datagcm'
           write(*,*) 'Also check that the continuum data continuum_data/H2-H2_norm_2011.cia or H2-H2_norm_2018.cia is there.'
           call abort
         else
         
         if(versH2H2cia.eq.2011) then
           write(*,*) '... You are using H2-H2 CIA from 2011 but you should use more recent data available on HITRAN !'
           write(*,*) '... (Especially if you are running a giant planet atmosphere)'
           write(*,*) '... Just find out the H2-H2_norm_2018.cia, put it in your datadir and have a look at interpolateH2H2.F90 ! .'
         endif

            do iT=1,nT
               
               ! Only two possibles values for now : 2011 or 2018 (sanity check in inifis_mod)
               if(versH2H2cia.eq.2011) then
                 read(33,fmat11) bleh,blah,blah,nres,Ttemp
               else if (versH2H2cia.eq.2018) then
                 read(33,fmat18) bleh,bleh,bleh,blah,blah,nres,Ttemp
               endif

               if(nS.ne.nres)then
                  print*,'Resolution given in file: ',trim(dt_file)
                  print*,'is ',nres,', which does not match nS.'
                  print*,'Please adjust nS value in interpolateH2H2.F90'
                  stop
               endif
               temp_arr(iT)=Ttemp

               do k=1,nS
                  read(33,*) wn_arr(k),abs_arr(k,it)
               end do

            end do
      
         endif
         close(33)
!$OMP END MASTER
!$OMP BARRIER

         print*,'interpolateH2H2: At wavenumber ',wn,' cm^-1'
         print*,'   temperature ',temp,' K'
         print*,'   pressure ',pres,' Pa'

      endif

         call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arr,wn,temp,abcoef,ind)

         !print*,'the absorption is ',abcoef,' cm^5 molecule^-2'
         !print*,'or ',abcoef*losch**2,' cm^-1 amagat^-2'

         abcoef=abcoef*losch**2*100.0*amagat**2 ! convert to m^-1

         !print*,'We have ',amagat,' amagats of H2'
         !print*,'So the absorption is ',abcoef,' m^-1'

         ! unlike for Rayleigh scattering, we do not currently weight by the BB function
         ! however our bands are normally thin, so this is no big deal.

      return
    end subroutine interpolateH2H2
