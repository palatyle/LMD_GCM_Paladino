










      subroutine setspv

!==================================================================
!     
!     Purpose
!     -------
!     Set up spectral intervals, stellar spectrum and Rayleigh 
!     opacity in the shortwave. 
!     
!     Authors
!     ------- 
!     Adapted from setspv in the NASA Ames radiative code by
!     Robin Wordsworth (2009).
!
!     Called by
!     ---------
!     callcorrk.F
!     
!     Calls
!     -----
!     ave_stelspec.F
!     
!==================================================================

      use radinc_h,    only: L_NSPECTV, corrkdir, banddir
      use radcommon_h, only: BWNV,BLAMV,WNOV,DWNV,WAVEV, &
                             STELLARF,TAURAY
      use datafile_mod, only: datadir
      use callkeys_mod, only: Fat1AU,rayleigh

      implicit none

      logical file_ok

      integer N, M, file_entries

      character(len=30)  :: temp1
      character(len=200) :: file_id
      character(len=200) :: file_path

      real*8 :: lastband(2)

      real*8 STELLAR(L_NSPECTV)
      real*8 sum, dummy

      !! used to count lines
      integer :: nb
      integer :: ierr

!=======================================================================
!     Set up spectral bands - wavenumber [cm^(-1)]. Go from smaller to
!     larger wavenumbers, the same as in the IR.

      write(temp1,'(i2.2)') L_NSPECTV
      file_id='/corrk_data/'//trim(adjustl(banddir))//'/narrowbands_VI.in' 
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by setspv.F90, exiting.'
         write(*,*)'Check that your path to datagcm:',trim(datadir)
         write(*,*)' is correct. You can change it in callphys.def with:'
         write(*,*)' datadir = /absolute/path/to/datagcm'
         write(*,*)'Also check that the corrkdir you chose in callphys.def exists.'
         call abort
      endif
	
!$OMP MASTER        
      nb=0
      ierr=0
      ! check that the file contains the right number of bands 
      open(131,file=file_path,form='formatted')
      read(131,*,iostat=ierr) file_entries
      do while (ierr==0)
        read(131,*,iostat=ierr) dummy
        if (ierr==0) nb=nb+1
      enddo
      close(131)

      write(*,*) 'setspv: L_NSPECTV = ',L_NSPECTV, 'in the model '
      write(*,*) '        there are   ',nb, 'entries in ',TRIM(file_path)
      if(nb.ne.L_NSPECTV) then
         write(*,*) 'MISMATCH !! I stop here'
         call abort
      endif

      ! load and display the data
      open(111,file=file_path,form='formatted')
      read(111,*) 
       do M=1,L_NSPECTV-1
         read(111,*) BWNV(M)
      end do
      read(111,*) lastband
      close(111)
      BWNV(L_NSPECTV)  =lastband(1)
      BWNV(L_NSPECTV+1)=lastband(2)
!$OMP END MASTER
!$OMP BARRIER

      print*,'setspv: VI band limits:'
      do M=1,L_NSPECTV+1
         print*,m,'-->',BWNV(M),' cm^-1'
      end do
      print*,' '

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTV
         WNOV(M)  = 0.5*(BWNV(M+1)+BWNV(M))
         DWNV(M)  = BWNV(M+1)-BWNV(M)
         WAVEV(M) = 1.0E+4/WNOV(M)
         BLAMV(M) = 0.01/BWNV(M)
      end do
      BLAMV(M) = 0.01/BWNV(M) ! wavelength in METERS for aerosol stuff
!     note M=L_NSPECTV+1 after loop due to Fortran bizarreness

!=======================================================================
!     Set up stellar spectrum

      write(*,*)'setspv: Interpolating stellar spectrum from the hires data...'
      call ave_stelspec(STELLAR)

!     Sum the stellar flux, and write out the result.  
      sum = 0.0  
      do N=1,L_NSPECTV
         STELLARF(N) = STELLAR(N) * Fat1AU
         sum         = sum+STELLARF(N)
      end do
      write(6,'("setspv: Stellar flux at 1 AU = ",f7.2," W m-2")') sum
      print*,' '


!=======================================================================
!     Set up the wavelength independent part of the Rayleigh scattering.
!     The pressure dependent part will be computed elsewhere (OPTCV).
!     WAVEV is in microns.  There is no Rayleigh scattering in the IR.

      if(rayleigh) then
         call calc_rayleigh
      else
         print*,'setspv: No Rayleigh scattering, check for NaN in output!'
         do N=1,L_NSPECTV
            TAURAY(N) = 1E-16
         end do
      endif

      RETURN
    END subroutine setspv
