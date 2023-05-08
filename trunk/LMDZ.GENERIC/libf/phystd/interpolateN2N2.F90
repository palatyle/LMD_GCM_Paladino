subroutine interpolateN2N2(wn,temp,pres,abcoef,firstcall,ind)

  !==================================================================
  !     
  !     Purpose
  !     -------
  !     Calculates the N2-N2 CIA opacity, using a lookup table from
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
  double precision pres               ! pressure               (Pascals)
  integer :: ind

  ! output
  double precision abcoef             ! absorption coefficient (m^-1)

  integer nS,nT
  parameter(nS=582)
  parameter(nT=10)

  double precision, parameter :: losch = 2.6867774e19
  ! Loschmit's number (molecule cm^-3 at STP) 
  ! converts cm^5 molecule^-2 --> cm^-1 amagat^-2
  ! see Richard et al. 2011, JQSRT for details

  double precision amagat
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


  if(temp.gt.400)then
     print*,'Your temperatures are too high for this N2-N2 CIA dataset.'
     print*,'Currently, HITRAN provides data for this pair in the range 40-400 K.'
     stop
  endif

  amagat = (273.15/temp)*(pres/101325.0)

  if(firstcall)then ! called by sugas_corrk only
     print*,'----------------------------------------------------'
     print*,'Initialising N2-N2 continuum from HITRAN database...'

     !     1.1 Open the ASCII files
     dt_file=TRIM(datadir)//'/continuum_data/N2-N2_2011.cia'

!$OMP MASTER
     open(33,file=dt_file,form='formatted',status='old',iostat=ios)
     if (ios.ne.0) then        ! file not found
        write(*,*) 'Error from interpolateN2N2'
        write(*,*) 'Data file ',trim(dt_file),' not found.'
        write(*,*) 'Check that your path to datagcm:',trim(datadir)
        write(*,*) 'is correct. You can change it in callphys.def with:'
        write(*,*) 'datadir = /absolute/path/to/datagcm'
        write(*,*) 'Also check that the continuum data continuum_data/N2-N2_norm_2011.cia is there.'
        call abort
     else

        do iT=1,nT

           read(33,fmat1) bleh,blah,blah,nres,Ttemp
           if(nS.ne.nres)then
              print*,'Resolution given in file: ',trim(dt_file)
              print*,'is ',nres,', which does not match nS.'
              print*,'Please adjust nS value in interpolateN2N2.F90'
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

     print*,'interpolateN2N2: At wavenumber ',wn,' cm^-1'
     print*,'   temperature ',temp,' K'
     print*,'   pressure ',pres,' Pa'

  endif
     call bilinearbig(nS,nT,wn_arr,temp_arr,abs_arr,wn,temp,abcoef,ind)

!     print*,'the absorption is ',abcoef,' cm^5 molecule^-2'
!     print*,'or ',abcoef*losch**2,' cm^-1 amagat^-2'

     abcoef=abcoef*losch**2*100.0*amagat**2 ! convert to m^-1
!     abcoef=0.

!     print*,'We have ',amagat,' amagats of N2'
!     print*,'So the absorption is ',abcoef,' m^-1'

!     unlike for Rayleigh scattering, we do not currently weight by the BB function
!     however our bands are normally thin, so this is no big deal.


  return
end subroutine interpolateN2N2

