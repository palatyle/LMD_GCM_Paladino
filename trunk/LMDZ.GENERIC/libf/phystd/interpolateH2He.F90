     subroutine interpolateH2He(wn,temp,presH2,presHe,abcoef,firstcall,ind)

!==================================================================
!     
!     Purpose
!     -------
!     Calculates the H2-He CIA opacity, using a lookup table from
!     HITRAN (2011)
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
      double precision presH2             ! H2 partial pressure    (Pascals)
      double precision presHe             ! He partial pressure    (Pascals)

      ! output
      double precision abcoef             ! absorption coefficient (m^-1)

      integer nS,nT
      parameter(nS=2428)
      parameter(nT=10)

      double precision, parameter :: losch = 2.6867774e19
      ! Loschmit's number (molecule cm^-3 at STP) 
      ! converts cm^5 molecule^-2 --> cm^-1 amagat^-2
      ! see Richard et al. 2011, JQSRT for details

      double precision amagatH2
      double precision amagatHe
      double precision wn_arr(nS)
      double precision temp_arr(nT)
      double precision abs_arr(nS,nT)

      integer k,iT
      logical firstcall

      save wn_arr, temp_arr, abs_arr !read by master

      character*100 dt_file
      integer strlen,ios

      character(LEN=*), parameter :: fmat1 = "(A20,F10.3,F10.3,I7,F7.1,E10.3,F5.3)"

      character*20 bleh
      double precision blah, Ttemp
      integer nres

      integer ind
 
      if(temp.gt.400)then
         print*,'Your temperatures are too high for this H2-He CIA dataset.'
         print*,'Please run mixed H2-He atmospheres below T = 400 K.'      
         stop
      endif

      amagatH2 = (273.15/temp)*(presH2/101325.0)
      amagatHe = (273.15/temp)*(presHe/101325.0)

      if(firstcall)then ! called by sugas_corrk only
         print*,'----------------------------------------------------'
         print*,'Initialising H2-He continuum from HITRAN database...'

!     1.1 Open the ASCII files
         dt_file=TRIM(datadir)//'/continuum_data/H2-He_norm_2011.cia'
	 
!$OMP MASTER
         open(33,file=dt_file,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from interpolateH2He'
           write(*,*) 'Data file ',trim(dt_file),' not found.'
           write(*,*) 'Check that your path to datagcm:',trim(datadir)
           write(*,*) 'is correct. You can change it in callphys.def with:'
           write(*,*) 'datadir = /absolute/path/to/datagcm'
           write(*,*) 'Also check that the continuum data continuum_data/H2-He_norm_2011.cia is there.'
           call abort
         else

            do iT=1,nT

               read(33,fmat1) bleh,blah,blah,nres,Ttemp
               if(nS.ne.nres)then
                  print*,'Resolution given in file: ',trim(dt_file)
                  print*,'is ',nres,', which does not match nS.'
                  print*,'Please adjust nS value in interpolateH2He.F90'
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

         print*,'interpolateH2He: At wavenumber ',wn,' cm^-1'
         print*,'   temperature                 ',temp,' K'
         print*,'   H2 partial pressure         ',presH2,' Pa'
         print*,'   and He partial pressure     ',presHe,' Pa'

      endif

         call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arr,wn,temp,abcoef,ind)

         !print*,'the absorption is ',abcoef,' cm^5 molecule^-2'
         !print*,'or ',abcoef*losch**2,' cm^-1 amagat^-2'

         abcoef=abcoef*losch**2*100.0*amagatH2*amagatHe ! convert to m^-1

         !print*,'We have ',amagatH2,' amagats of H2'
         !print*,'and     ',amagatHe,' amagats of He'
         !print*,'So the absorption is ',abcoef,' m^-1'

         ! unlike for Rayleigh scattering, we do not currently weight by the BB function
         ! however our bands are normally thin, so this is no big deal.

      return
    end subroutine interpolateH2He
