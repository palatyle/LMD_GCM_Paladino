










      subroutine ave_stelspec(STELLAR)

!==================================================================
!     
!     Purpose
!     -------
!     Average the chosen high resolution stellar spectrum over the
!     visible bands in the model.
!     
!     Authors
!     ------- 
!     Robin Wordsworth (2010).
!     Generalized to very late spectral types (and Brown dwarfs) Jeremy Leconte (2012)
!
!     Called by
!     ---------
!     setspv.F
!     
!     Calls
!     -----
!     none
!     
!==================================================================

      use radinc_h, only: L_NSPECTV
      use radcommon_h, only: BWNV, DWNV, tstellar
      use datafile_mod, only: datadir
      use callkeys_mod, only: stelbbody,stelTbb,startype

      implicit none

      real*8 STELLAR(L_NSPECTV)
!      integer startype 

      integer Nfine
      integer,parameter :: Nfineband=200
      integer ifine,band

      real,allocatable,save :: lam(:),stel_f(:) 	!read by master
      real lamm,lamp
      real dl

      character(len=100)  :: file_id,file_id_lam
      character(len=200) :: file_path,file_path_lam

      real lam_temp
      double precision stel_temp
      
      integer :: ios ! file opening/reading status

      STELLAR(:)=0.0

      print*,'enter ave_stellspec'
      if(stelbbody)then
         tstellar=stelTbb
	 Nfine=L_NSPECTV*Nfineband
	 do band=1,L_NSPECTV
	    lamm=10000.0/BWNV(band+1)
	    lamp=10000.0/BWNV(band)
	    dl=(lamp-lamm)/(Nfineband)
	    do ifine=1,Nfineband
	       lam_temp=lamm+(lamp-lamm)*(ifine-1.)/(Nfineband)
               call blackl(dble(lam_temp*1e-6),dble(tstellar),stel_temp)
	       STELLAR(band)=STELLAR(band)+stel_temp*dl
	    enddo	    
	 end do
         STELLAR(1:L_NSPECTV)=STELLAR(1:L_NSPECTV)/sum(STELLAR(1:L_NSPECTV))
      else
         ! load high resolution stellar data
         Select Case(startype)
	   Case(1)
            file_id='/stellar_spectra/sol.txt'
            tstellar=5800.
            file_id_lam='/stellar_spectra/lam.txt'
            Nfine=5000
	   Case(2)
            file_id='/stellar_spectra/gl581.txt'
            tstellar=3200.
            file_id_lam='/stellar_spectra/lam.txt'
            Nfine=5000
	   Case(3)
            file_id='/stellar_spectra/adleo.txt'
            tstellar=3200.
            file_id_lam='/stellar_spectra/lam.txt'
            Nfine=5000
	   Case(4)
            file_id='/stellar_spectra/gj644.txt'
            print*,'Find out tstellar before using this star!'
            call abort
            file_id_lam='/stellar_spectra/lam.txt'
            Nfine=5000
	   Case(5)
            file_id='/stellar_spectra/hd128167.txt'
            tstellar=6700. ! Segura et al. (2003)
            file_id_lam='/stellar_spectra/lam.txt'
            Nfine=5000
	   Case(6)
            file_id='/stellar_spectra/BD_Teff-1600K.txt'
            tstellar=1600. 
            file_id_lam='/stellar_spectra/lamBD.txt'
            Nfine=5000
	   Case(7)
            file_id='/stellar_spectra/BD_Teff-1000K.txt'
            tstellar=1000. 
            file_id_lam='/stellar_spectra/lamBD.txt'
            Nfine=5000
	   Case(8)
            file_id='/stellar_spectra/Flux_K5_Teff4700_logg4.5_Met-0.5_BTsettle.dat'
            tstellar=4700. 
            file_id_lam='/stellar_spectra/lambda_K5_Teff4700_logg4.5_Met-0.5_BTsettle.dat'
            Nfine=3986
           Case(9)
            file_id='/stellar_spectra/Flux_TRAPPIST1.dat'
            tstellar=2550. 
            file_id_lam='/stellar_spectra/lambda_TRAPPIST1.dat'
            Nfine=5000
           Case(10)
            file_id='/stellar_spectra/Flux_Proxima.dat'
            tstellar=3050. 
            file_id_lam='/stellar_spectra/lambda_Proxima.dat'
            Nfine=5000
	   Case Default
            print*,'Error: unknown star type chosen'
            call abort
         End Select

!$OMP MASTER
	 allocate(lam(Nfine),stel_f(Nfine))

         file_path_lam=TRIM(datadir)//TRIM(file_id_lam)
         open(110,file=file_path_lam,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from ave_stelspec'
           write(*,*) 'Data file ',trim(file_id_lam),' not found.'
           write(*,*)'Check that your path to datagcm:',trim(datadir)
           write(*,*)' is correct. You can change it in callphys.def with:'
           write(*,*)' datadir = /absolute/path/to/datagcm'
           write(*,*)' Also check that there is a ',trim(file_id_lam),' there.'
           call abort
         else
           do ifine=1,Nfine
             read(110,*) lam(ifine)
           enddo
           close(110)
         endif


         ! load high resolution wavenumber data
         file_path=TRIM(datadir)//TRIM(file_id)
         open(111,file=file_path,form='formatted',status='old',iostat=ios)
         if (ios.ne.0) then        ! file not found
           write(*,*) 'Error from ave_stelspec'
           write(*,*) 'Data file ',trim(file_id),' not found.'
           write(*,*)'Check that your path to datagcm:',trim(datadir)
           write(*,*)' is correct. You can change it in callphys.def with:'
           write(*,*)' datadir = /absolute/path/to/datagcm'
           write(*,*)' Also check that there is a ',trim(file_id),' there.'
           call abort
         else
           do ifine=1,Nfine
             read(111,*) stel_f(ifine)
           enddo
           close(111)
         endif
!$OMP END MASTER
!$OMP BARRIER
	 
         ! sum data by band
         band=1
	 Do while(lam(1).lt. real(10000.0/BWNV(band+1)))
	    if (band.gt.L_NSPECTV-1) exit
            band=band+1
	 enddo
	 dl=lam(2)-lam(1)
         STELLAR(band)=STELLAR(band)+stel_f(1)*dl
         do ifine = 2,Nfine
            if(lam(ifine) .gt. real(10000.0/BWNV(band)))then
               band=band-1
            endif
            if(band .lt. 1) exit
	    dl=lam(ifine)-lam(ifine-1)
            STELLAR(band)=STELLAR(band)+stel_f(ifine)*dl
         end do
	       
	 
         STELLAR(1:L_NSPECTV)=STELLAR(1:L_NSPECTV)/sum(STELLAR(1:L_NSPECTV))
!$OMP BARRIER
!$OMP MASTER
	 if (allocated(lam)) deallocate(lam)
	 if (allocated(stel_f)) deallocate(stel_f)
!$OMP END MASTER
!$OMP BARRIER         
      endif

      end subroutine ave_stelspec
