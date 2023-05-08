










      subroutine setspi

!==================================================================
!     
!     Purpose
!     -------
!     Set up spectral intervals and Planck function in the longwave.
!     
!     Authors
!     ------- 
!     Adapted from setspi in the NASA Ames radiative code by
!     Robin Wordsworth (2009).
!     
!     Called by
!     ---------
!     callcorrk.F
!     
!     Calls
!     -----
!     none
!     
!==================================================================

      use radinc_h,    only: L_NSPECTI,corrkdir,banddir,NTstart,NTstop,NTfac
      use radcommon_h, only: BWNI,BLAMI,WNOI,DWNI,WAVEI,planckir,sigma
      use datafile_mod, only: datadir
      use comcstfi_mod, only: pi

      implicit none

      logical file_ok
      integer nw, nt, m, mm, file_entries
      real*8 a, b, ans, y, bpa, bma, T, dummy

      character(len=30)  :: temp1
      character(len=200) :: file_id
      character(len=200) :: file_path

!     C1 and C2 values from Goody and Yung (2nd edition)  MKS units
!     These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4

      real*8 :: c1 = 3.741832D-16 ! W m^-2
      real*8 :: c2 = 1.438786D-2  ! m K
      
      real*8 :: lastband(2), plancksum

      !! used to count lines
      integer :: nb
      integer :: ierr

      logical forceEC, planckcheck

      real*8 :: x(12) = [ -0.981560634246719D0,  -0.904117256370475D0, &
      -0.769902674194305D0,  -0.587317954286617D0,                     &
      -0.367831498998180D0,  -0.125233408511469D0,                     &
       0.125233408511469D0,   0.367831498998180D0,                     &
       0.587317954286617D0,   0.769902674194305D0,                     &
       0.904117256370475D0,   0.981560634246719D0  ]

      real*8 :: w(12) = [  0.047175336386512D0,   0.106939325995318D0, &
           0.160078328543346D0,   0.203167426723066D0,                 &
           0.233492536538355D0,   0.249147045813403D0,                 &
           0.249147045813403D0,   0.233492536538355D0,                 &
           0.203167426723066D0,   0.160078328543346D0,                 &
           0.106939325995318D0,   0.047175336386512D0  ]
      mm=0

      forceEC=.true.
      planckcheck=.true.

!=======================================================================
!     Set up spectral bands - wavenumber [cm^(-1)]. Go from smaller to
!     larger wavenumbers.

      write(temp1,'(i2.2)') L_NSPECTI
      !file_id='/corrk_data/' // corrkdir(1:LEN_TRIM(corrkdir)) // '/narrowbands_IR.in'
      file_id='/corrk_data/'//trim(adjustl(banddir))//'/narrowbands_IR.in' 
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by setspi.F90, exiting.'
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
!        write(*,*) 'setspi: file_entries:',dummy,'ierr=',ierr
        if (ierr==0) nb=nb+1
      enddo
      close(131)

      write(*,*) 'setspi: L_NSPECTI = ',L_NSPECTI, 'in the model '
      write(*,*) '        there are   ',nb, 'entries in ',TRIM(file_path)
      if(nb.ne.L_NSPECTI) then
         write(*,*) 'MISMATCH !! I stop here'
         call abort
      endif

      ! load and display the data
      open(111,file=file_path,form='formatted')
      read(111,*) 
      do M=1,L_NSPECTI-1
         read(111,*) BWNI(M)
      end do
      read(111,*) lastband
      close(111)
      BWNI(L_NSPECTI)  =lastband(1)
      BWNI(L_NSPECTI+1)=lastband(2)
!$OMP END MASTER
!$OMP BARRIER

      print*,''
      print*,'setspi: IR band limits:'
      do M=1,L_NSPECTI+1
         print*,m,'-->',BWNI(M),' cm^-1'
      end do

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTI
         WNOI(M)  = 0.5D0*(BWNI(M+1)+BWNI(M))
         DWNI(M)  = BWNI(M+1)-BWNI(M)
         WAVEI(M) = 1.0D+4/WNOI(M)
         BLAMI(M) = 0.01D0/BWNI(M)         
      end do
      BLAMI(M) = 0.01D0/BWNI(M)
!     note M=L_NSPECTI+1 after loop due to Fortran bizarreness

!=======================================================================
!     For each IR wavelength interval, compute the integral of B(T), the
!     Planck function, divided by the wavelength interval, in cm-1.  The
!     integration is in MKS units, the final answer is the same as the
!     original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.

      print*,''
      print*,'setspi: Current Planck integration range:'
      print*,'T = ',dble(NTstart)/NTfac, ' to ',dble(NTstop)/NTfac,' K.'

      IF(.NOT.ALLOCATED(planckir)) ALLOCATE(planckir(L_NSPECTI,NTstop-NTstart+1))

      do NW=1,L_NSPECTI
         a = 1.0D-2/BWNI(NW+1)
         b = 1.0D-2/BWNI(NW)
         bpa = (b+a)/2.0D0
         bma = (b-a)/2.0D0
         do nt=NTstart,NTstop
            T   = dble(NT)/NTfac
            ans = 0.0D0

            do mm=1,12
               y    = bma*x(mm)+bpa
               ans  = ans + w(mm)*c1/(y**5*(exp(c2/(y*T))-1.0D0))
            end do

            planckir(NW,nt-NTstart+1) = ans*bma/(PI*DWNI(NW))
         end do
      end do
         
      ! force planck=sigma*eps*T^4 for each temperature in array
      if(forceEC)then
         print*,'setspi: Force F=sigma*eps*T^4 for all values of T!'
         do nt=NTstart,NTstop
            plancksum=0.0D0
            T=dble(NT)/NTfac
       
            do NW=1,L_NSPECTI
               plancksum=plancksum+  &
                  planckir(NW,nt-NTstart+1)*DWNI(NW)*pi
            end do

            do NW=1,L_NSPECTI
               planckir(NW,nt-NTstart+1)=     &
                  planckir(NW,nt-NTstart+1)*  &
                          sigma*(dble(nt)/NTfac)**4/plancksum
            end do
         end do
      endif

      if(planckcheck)then
         ! check energy conservation at lower temperature boundary
         plancksum=0.0D0
         nt=NTstart
         do NW=1,L_NSPECTI
            plancksum=plancksum+planckir(NW,nt-NTstart+1)*DWNI(NW)*pi
         end do
         print*,'setspi: At lower limit:'
         print*,'in model sig*T^4 = ',plancksum,' W m^-2'
         print*,'actual sig*T^4   = ',sigma*(dble(nt)/NTfac)**4,' W m^-2'
         
         ! check energy conservation at upper temperature boundary
         plancksum=0.0D0
         nt=NTstop
         do NW=1,L_NSPECTI
            plancksum=plancksum+planckir(NW,nt-NTstart+1)*DWNI(NW)*pi
         end do
         print*,'setspi: At upper limit:'
         print*,'in model sig*T^4 = ',plancksum,' W m^-2'
         print*,'actual sig*T^4   = ',sigma*(dble(nt)/NTfac)**4,' W m^-2'
         print*,''
      endif

      return
    end subroutine setspi
