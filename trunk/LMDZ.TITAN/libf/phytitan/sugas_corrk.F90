      subroutine sugas_corrk

!==================================================================
!
!     Purpose
!     -------
!     Set up gaseous absorption parameters used by the radiation code.
!     This subroutine is a replacement for the old 'setrad', which contained
!     both absorption and scattering data.
!
!     Authors
!     -------
!     Adapted and generalised from the NASA Ames code by Robin Wordsworth (2009)
!     Added double gray case by Jeremy Leconte (2012)
!     New HITRAN continuum data section by RW (2012)
!
!     Summary
!     -------
!
!==================================================================

      use radinc_h
      use radcommon_h, only : pgasref,pfgasref,pgasmin,pgasmax
      use radcommon_h, only : tgasref,tgasmin,tgasmax
      use radcommon_h, only : gasv,gasi,FZEROI,FZEROV,gweight, w_cum
      use radcommon_h, only : WNOI,WNOV
      use radcommon_h, only : gasi_recomb, gasv_recomb, pqrold, useptold
      use radcommon_h, only : radvar_mask, radvar_indx, permut_idx
      use datafile_mod, only: datadir, corrkdir, banddir
      use comcstfi_mod, only: mugaz
      use gases_h
      use ioipsl_getin_p_mod, only: getin_p
      use callkeys_mod, only: graybody,callgasvis, continuum, tracer, corrk_recombin
      use tracer_h, only: noms, nqtot_p
      implicit none

!==================================================================

      logical file_ok

      integer n, nt, np, nh, ng, nw, m, i

      character(len=200) :: file_id
      character(len=500) :: file_path

      ! ALLOCATABLE ARRAYS -- AS 12/2011
      REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE,SAVE :: gasi8, gasv8 	!read by master
      character*20,allocatable,DIMENSION(:),SAVE :: gastype ! for name check, read by master

      real*8 x, xi(4), yi(4), ans, p
!     For gray case (JL12)
      real kappa_IR, kappa_VI, IR_VI_wnlimit
      integer nVI_limit,nIR_limit

      integer ngas, igas, jgas

      double precision testcont ! for continuum absorption initialisation

      integer :: dummy
      
      integer ngasvar
      logical found

!=======================================================================
!     Load variable species data, exit if we have wrong database
      file_id='/corrk_data/' // TRIM(corrkdir) // '/Q.dat'
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by sugas_corrk.F90, exiting.'
         write(*,*)'Check that your path to datagcm:',trim(datadir)
         write(*,*)' is correct. You can change it in callphys.def with:'
         write(*,*)' datadir = /absolute/path/to/datagcm'
         write(*,*)'Also check that the corrkdir you chose in callphys.def exists.'
         call abort
      endif

!$OMP MASTER
      ! check that database matches varactive toggle
      open(111,file=TRIM(file_path),form='formatted')
      read(111,*) ngas

      if(ngas.lt.1)then
         print*,ngas,' species in database [',               &
                     corrkdir(1:LEN_TRIM(corrkdir)),           &
                '], radiative code cannot handle this.'
         call abort
      endif
      
      if ( corrk_recombin ) then
         L_REFVAR =  ngas
         
         ! dynamical allocations and read from Q.dat
         IF ( .NOT. ALLOCATED( gastype ) ) ALLOCATE( gastype( L_REFVAR ) )
         IF ( .NOT. ALLOCATED( radvar_mask ) ) ALLOCATE( radvar_mask( L_REFVAR ) )         
         IF ( .NOT. ALLOCATED( radvar_indx ) ) ALLOCATE( radvar_indx( L_REFVAR ) )
         
         ngasvar = 0
         
         write(*,*) ' '
         do igas=1,L_REFVAR
           read(111,*) gastype(igas)
           read(111,*) radvar_mask(igas) ! boolean
           print *, 'Found radiatively active specie : ',trim(gastype(igas)), ' : Variable =', radvar_mask(igas)
           if ( radvar_mask(igas) ) ngasvar = ngasvar+1
         enddo
         
         print *, ' '
         print *, 'TOTAL number of species found for radiative transfer = ', L_REFVAR
         print *, '      including : ', ngasvar, 'variable species.'
         print *, ' '
        
         if (.not. tracer .and. ngasvar .gt. 0) then
            write(*,*) 'sugas_corrk:error: You try to run variable species for corrk recombining'
            write(*,*) 'in the radiative transfer but without tracers : I will stop here ! '
            call abort
         endif
         
         ! NB (JVO 18) We should enable here to force composition even if no associated tracers
      
         ! Check that we have matching tracers for the variable radiative active species
         if ( ngasvar .GT. 0 ) then

           do igas=1,L_REFVAR
            if ( radvar_mask(igas) ) then 
            
               found = .false.
               do jgas=1,nqtot_p
                 if ( trim(noms(jgas)) == trim(gastype(igas)) .OR. &
                      trim(noms(jgas)) == trim(gastype(igas))//'_') then
                    radvar_indx(igas) = jgas
                    found=.true.
                    exit
                 endif
               enddo
               
               if (.not. found) then
                 write(*,*) "sugas_corrk:error: "//trim(gastype(igas))//" is missing from tracers."
                 call abort
               endif
               
             endif
           enddo ! igas=1,L_REFVAR
           
         endif ! if ngasvar .gt. 0
      
      else
        L_REFVAR = 1 ! JVO 2017 : those will not be used, just set to 1 to keep the code running
      endif ! if corrk_recombin

!=======================================================================
!     Set the weighting in g-space

      file_id='/corrk_data/' // TRIM(corrkdir) // '/g.dat'
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by sugas_corrk.F90, exiting.'
         write(*,*)'Check that your path to datagcm:',trim(datadir)
         write(*,*)' is correct. You can change it in callphys.def with:'
         write(*,*)' datadir = /absolute/path/to/datagcm'
         write(*,*)'Also check that the corrkdir you chose in callphys.def exists.'
         call abort
      endif
      
      ! check the array size is correct, load the coefficients
      open(111,file=TRIM(file_path),form='formatted')
      read(111,*) L_NGAUSS
      IF( .NOT. ALLOCATED( gweight ) ) ALLOCATE( GWEIGHT(L_NGAUSS) )
      read(111,*) gweight
      close(111)
 
      ! display the values
      print*,'Correlated-k g-space grid:'
      do n=1,L_NGAUSS
         print*,n,'.',gweight(n)
      end do
      print*,''

      IF (.NOT. ALLOCATED(permut_idx)) ALLOCATE(permut_idx(L_NGAUSS*L_NGAUSS)) ! for corr-k recombining
      permut_idx = (/(i, i=1,L_NGAUSS*L_NGAUSS)/) ! for the recombin_corrk firstcall

      IF (.NOT. ALLOCATED(w_cum)) ALLOCATE(w_cum(L_NGAUSS)) ! for corr-k recombining      
      w_cum(1)= gweight(1)
      DO n=2,L_NGAUSS
        w_cum(n) = w_cum(n-1)+gweight(n)
      ENDDO

!=======================================================================
!     Set the reference pressure and temperature arrays.  These are
!     the pressures and temperatures at which we have k-coefficients.

!-----------------------------------------------------------------------
! pressure

      file_id='/corrk_data/' // TRIM(corrkdir) // '/p.dat'
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by sugas_corrk.F90, exiting.'
         write(*,*)'Check that your path to datagcm:',trim(datadir)
         write(*,*)' is correct. You can change it in callphys.def with:'
         write(*,*)' datadir = /absolute/path/to/datagcm'
         write(*,*)'Also check that the corrkdir you chose in callphys.def exists.'
         call abort
      endif
     
      ! get array size, load the coefficients
      open(111,file=TRIM(file_path),form='formatted')
      read(111,*) L_NPREF
      IF( .NOT. ALLOCATED( pgasref ) ) ALLOCATE( PGASREF(L_NPREF) )
      read(111,*) pgasref
      close(111)
      L_PINT = (L_NPREF-1)*5+1
      IF( .NOT. ALLOCATED( pfgasref ) ) ALLOCATE( PFGASREF(L_PINT) )

      ! display the values
      print*,'Correlated-k pressure grid (mBar):'
      do n=1,L_NPREF
         print*,n,'. 1 x 10^',pgasref(n),' mBar'
      end do
      print*,''

      ! save the min / max matrix values
      pgasmin = 10.0**pgasref(1)
      pgasmax = 10.0**pgasref(L_NPREF)

      ! interpolate to finer grid, adapted to uneven grids
      do n=1,L_NPREF-1
         do m=1,5
            pfgasref((n-1)*5+m) = pgasref(n)+(m-1)*(pgasref(n+1) - pgasref(n))/5.
         end do
      end do
      pfgasref(L_PINT) = pgasref(L_NPREF)

!-----------------------------------------------------------------------
! temperature

      file_id='/corrk_data/' // TRIM(corrkdir) // '/T.dat'
      file_path=TRIM(datadir)//TRIM(file_id)

      ! check that the file exists
      inquire(FILE=file_path,EXIST=file_ok)
      if(.not.file_ok) then
         write(*,*)'The file ',TRIM(file_path)
         write(*,*)'was not found by sugas_corrk.F90, exiting.'
         write(*,*)'Check that your path to datagcm:',trim(datadir)
         write(*,*)' is correct. You can change it in callphys.def with:'
         write(*,*)' datadir = /absolute/path/to/datagcm'
         write(*,*)'Also check that the corrkdir you chose in callphys.def exists.'
         call abort
      endif

      ! get array size, load the coefficients
      open(111,file=TRIM(file_path),form='formatted')
      read(111,*) L_NTREF
      IF( .NOT. ALLOCATED( tgasref ) ) ALLOCATE( TGASREF(L_NTREF) )
      read(111,*) tgasref
      close(111)

      ! display the values
      print*,'Correlated-k temperature grid:'
      do n=1,L_NTREF
         print*,n,'.',tgasref(n),' K'
      end do

      ! save the min / max matrix values
      tgasmin = tgasref(1)
      tgasmax = tgasref(L_NTREF)

      IF( .NOT. ALLOCATED( gasi8 ) ) ALLOCATE( gasi8(L_NTREF,L_NPREF,L_REFVAR,L_NSPECTI,L_NGAUSS) )
      IF( .NOT. ALLOCATED( gasv8 ) ) ALLOCATE( gasv8(L_NTREF,L_NPREF,L_REFVAR,L_NSPECTV,L_NGAUSS) )
!$OMP END MASTER
!$OMP BARRIER

!-----------------------------------------------------------------------
! allocate the multidimensional arrays in radcommon_h
      IF( .NOT. ALLOCATED( gasi ) ) ALLOCATE( gasi(L_NTREF,L_PINT,L_REFVAR,L_NSPECTI,L_NGAUSS) )
      IF( .NOT. ALLOCATED( gasv ) ) ALLOCATE( gasv(L_NTREF,L_PINT,L_REFVAR,L_NSPECTV,L_NGAUSS) )

      if (corrk_recombin) then
        IF( .NOT. ALLOCATED( gasi_recomb ) ) ALLOCATE( gasi_recomb(L_NTREF,L_PINT,L_NSPECTI,L_NGAUSS) )
        IF( .NOT. ALLOCATED( gasv_recomb ) ) ALLOCATE( gasv_recomb(L_NTREF,L_PINT,L_NSPECTV,L_NGAUSS) )
        IF( .NOT. ALLOCATED( pqrold ) ) ALLOCATE( pqrold(L_PINT,L_REFVAR) )
        IF( .NOT. ALLOCATED( useptold ) ) ALLOCATE( useptold(L_PINT,L_NTREF) )
        pqrold(:,:) = 0.0
        useptold(:,:) = .false.
      endif
      
      ! display the values
      print*,''
      print*,'Correlated-k matrix size:' 
      print*,'[',L_NTREF,',',L_NPREF,',',L_REFVAR,',',L_NGAUSS,']' 

!=======================================================================
!     Get gaseous k-coefficients and interpolate onto finer pressure grid


!        wavelength used to separate IR from VI in graybody. We will need that anyway
         IR_VI_wnlimit=3000.
         write(*,*)"graybody: Visible / Infrared separation set at",10000./IR_VI_wnlimit,"um"
	 
	 nVI_limit=0
	 Do nw=1,L_NSPECTV
	    if ((WNOV(nw).gt.IR_VI_wnlimit).and.(L_NSPECTV.gt.1)) then
	       nVI_limit=nw-1
	       exit
	    endif
	 End do
	 nIR_limit=L_NSPECTI
	 Do nw=1,L_NSPECTI
	    if ((WNOI(nw).gt.IR_VI_wnlimit).and.(L_NSPECTI.gt.1)) then
	       nIR_limit=nw-1
	       exit
	    endif
	 End do

      if (graybody) then
!        constant absorption coefficient in visible
         write(*,*)"graybody: constant absorption coefficient in visible:"
         kappa_VI=-100000.
         call getin_p("kappa_VI",kappa_VI)
         write(*,*)" kappa_VI = ",kappa_VI
	 kappa_VI=kappa_VI*1.e4* mugaz * 1.672621e-27	 ! conversion from m^2/kg to cm^2/molecule         
      
!        constant absorption coefficient in IR
         write(*,*)"graybody: constant absorption coefficient in InfraRed:"
         kappa_IR=-100000.
         call getin_p("kappa_IR",kappa_IR)
         write(*,*)" kappa_IR = ",kappa_IR	 
         kappa_IR=kappa_IR*1.e4* mugaz * 1.672621e-27	 ! conversion from m^2/kg to cm^2/molecule 

         write(*,*)"graybody: Visible / Infrared separation set at band: IR=",nIR_limit,", VI=",nVI_limit
	       
      Else
         kappa_VI=1.e-30      
         kappa_IR=1.e-30        
      End if

!$OMP MASTER         
!      print*,corrkdir(1:4)
      ! VISIBLE
      if (callgasvis) then
         if ((corrkdir(1:4).eq.'null'))then   !(TRIM(corrkdir).eq.'null_LowTeffStar')) then
            gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)=0.0
            print*,'using no corrk data'
            print*,'Visible corrk gaseous absorption is set to zero if graybody=F'
	 else
             
            if (corrk_recombin .eqv. .false. ) then
              file_id='/corrk_data/'//trim(adjustl(banddir))//'/corrk_gcm_VI.dat' 
              file_path=TRIM(datadir)//TRIM(file_id)
              
              ! check that the file exists
              inquire(FILE=file_path,EXIST=file_ok)
              if(.not.file_ok) then
                 write(*,*)'The file ',TRIM(file_path)
                 write(*,*)'was not found by sugas_corrk.F90.'
                 write(*,*)'Are you sure you have absorption data for these bands?'
                 call abort
              endif
           
              open(111,file=TRIM(file_path),form='formatted')
              read(111,*) gasv8
              close(111)
              
            else
              
              do igas=1,L_REFVAR
                ! JVO 18 : To comply with Q.dat watch out for 2-letters gas trailing underscore :
                ! e.g. if ever you have N2_ (symmetric so shouldn't happen) file should be ..._N2_.dat
                file_id='/corrk_data/'//trim(adjustl(banddir))//'/corrk_gcm_VI_'//trim(adjustl(gastype(igas)))//'.dat' 
                file_path=TRIM(datadir)//TRIM(file_id)
                
                ! check that the file exists
                inquire(FILE=file_path,EXIST=file_ok)
                if(.not.file_ok) then
                   write(*,*)'The file ',TRIM(file_path)
                   write(*,*)'was not found by sugas_corrk.F90.'
                   write(*,*)'Are you sure you have absorption data for these bands and this specie?'
                   call abort
                 endif
           
                open(111,file=TRIM(file_path),form='formatted')
                read(111,*) gasv8(:,:,igas,:,:)
                close(111)
              enddo              
              
            endif !  if corrk_recombin
            
	 end if

         if(nVI_limit.eq.0) then
	    gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)=   &
	          gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)+kappa_VI
         else if (nVI_limit.eq.L_NSPECTV) then
	    gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)=   &
	          gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)+kappa_IR
	 else
	    gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:nVI_limit,1:L_NGAUSS)=   &
	          gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:nVI_limit,1:L_NGAUSS)+kappa_IR
	    gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,nVI_limit+1:L_NSPECTV,1:L_NGAUSS)=   &
	          gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,nVI_limit+1:L_NSPECTV,1:L_NGAUSS)+kappa_VI
	 end if
      else
         print*,'Visible corrk gaseous absorption is set to zero.'
         gasv8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTV,1:L_NGAUSS)=0.0
      endif
!$OMP END MASTER
!$OMP BARRIER

      ! INFRA-RED
      if ((corrkdir(1:4).eq.'null'))then       !.or.(TRIM(corrkdir).eq.'null_LowTeffStar')) then
         print*,'Infrared corrk gaseous absorption is set to zero if graybody=F'
!$OMP MASTER         
         gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTI,1:L_NGAUSS)=0.0
!$OMP END MASTER
!$OMP BARRIER
      else 
      
         if ( corrk_recombin .eqv. .false. ) then
           file_id='/corrk_data/'//trim(adjustl(banddir))//'/corrk_gcm_IR.dat' 
           file_path=TRIM(datadir)//TRIM(file_id)
        
           ! check that the file exists
           inquire(FILE=file_path,EXIST=file_ok)
           if(.not.file_ok) then
              write(*,*)'The file ',TRIM(file_path)
              write(*,*)'was not found by sugas_corrk.F90.'
              write(*,*)'Are you sure you have absorption data for these bands?'
              call abort
           endif
           
!$OMP MASTER          	        
           open(111,file=TRIM(file_path),form='formatted')
           read(111,*) gasi8
           close(111)
!$OMP END MASTER
!$OMP BARRIER

         else
         
           do igas=1,L_REFVAR
             ! JVO 18 : To comply with Q.dat watch out for 2-letters gas trailing underscore :
             ! e.g. if ever tou have N2_ (symmetric so shouldn't happen) file should be ..._N2_.dat
             file_id='/corrk_data/'//trim(adjustl(banddir))//'/corrk_gcm_IR_'//trim(adjustl(gastype(igas)))//'.dat' 
             file_path=TRIM(datadir)//TRIM(file_id)
          
             ! check that the file exists
             inquire(FILE=file_path,EXIST=file_ok)
             if(.not.file_ok) then
                write(*,*)'The file ',TRIM(file_path)
                write(*,*)'was not found by sugas_corrk.F90.'
                write(*,*)'Are you sure you have absorption data for these bands and this specie?'
                call abort
             endif
             
!$OMP MASTER          	        
             open(111,file=TRIM(file_path),form='formatted')
             read(111,*) gasi8(:,:,igas,:,:)
             close(111)
!$OMP END MASTER
!$OMP BARRIER
           enddo ! do igas=1,L_REFVAR
           
         endif! if corrk_recombin
         
     
         ! 'fzero' is a currently unused feature that allows optimisation
         ! of the radiative transfer by neglecting bands where absorption
         ! is close to zero. As it could be useful in the future, this 
         ! section of the code has been kept commented and not erased.
         ! RW 7/3/12.

         do nw=1,L_NSPECTI
            fzeroi(nw) = 0.d0
!            do nt=1,L_NTREF
!               do np=1,L_NPREF
!                  do nh=1,L_REFVAR
!                     do ng = 1,L_NGAUSS
!                        if(gasi8(nt,np,nh,nw,ng).lt.1.0e-25)then
!                           fzeroi(nw)=fzeroi(nw)+1.d0
!                        endif
!                     end do
!                  end do
!               end do
!            end do
!            fzeroi(nw)=fzeroi(nw)/dble(L_NTREF*L_NPREF*L_REFVAR*L_NGAUSS)
         end do

         do nw=1,L_NSPECTV
            fzerov(nw) = 0.d0
!            do nt=1,L_NTREF
!               do np=1,L_NPREF
!                  do nh=1,L_REFVAR
!                     do ng = 1,L_NGAUSS
!                        if(gasv8(nt,np,nh,nw,ng).lt.1.0e-25)then
!                           fzerov(nw)=fzerov(nw)+1.d0
!                        endif
!                     end do
!                  end do
!               end do
!            end do
!            fzerov(nw)=fzerov(nw)/dble(L_NTREF*L_NPREF*L_REFVAR*L_NGAUSS)
         end do

      endif

!$OMP MASTER         	     
      if(nIR_limit.eq.0) then
         gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTI,1:L_NGAUSS)=   &
	          gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTI,1:L_NGAUSS)+kappa_VI
      else if (nIR_limit.eq.L_NSPECTI) then
	 gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTI,1:L_NGAUSS)=   &
	          gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:L_NSPECTI,1:L_NGAUSS)+kappa_IR
      else
	 gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:nIR_limit,1:L_NGAUSS)=   &
	          gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,1:nIR_limit,1:L_NGAUSS)+kappa_IR
	 gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,nIR_limit+1:L_NSPECTI,1:L_NGAUSS)=   &
	          gasi8(1:L_NTREF,1:L_NPREF,1:L_REFVAR,nIR_limit+1:L_NSPECTI,1:L_NGAUSS)+kappa_VI
      end if


!     Take log10 of the values - this is what we will interpolate.
!     Smallest value is 1.0E-200.

      do nt=1,L_NTREF
         do np=1,L_NPREF
            do nh=1,L_REFVAR
               do ng = 1,L_NGAUSS

                  do nw=1,L_NSPECTV
                     if(gasv8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        gasv8(nt,np,nh,nw,ng) = log10(gasv8(nt,np,nh,nw,ng))
                     else
                        gasv8(nt,np,nh,nw,ng) = -200.0
                     end if
                  end do

                  do nw=1,L_NSPECTI
                     if(gasi8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        gasi8(nt,np,nh,nw,ng) = log10(gasi8(nt,np,nh,nw,ng))
                     else
                        gasi8(nt,np,nh,nw,ng) = -200.0
                     end if
                  end do
                  
               end do
            end do
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER

!     Interpolate the values:  first the longwave

      do nt=1,L_NTREF
         do nh=1,L_REFVAR
            do nw=1,L_NSPECTI
               do ng=1,L_NGAUSS

!     First, the initial interval

                  n = 1 
                  do m=1,5
                     x     = pfgasref(m)
                     xi(1) = pgasref(n)
                     xi(2) = pgasref(n+1)
                     xi(3) = pgasref(n+2)
                     xi(4) = pgasref(n+3)
                     yi(1) = gasi8(nt,n,nh,nw,ng)
                     yi(2) = gasi8(nt,n+1,nh,nw,ng)
                     yi(3) = gasi8(nt,n+2,nh,nw,ng)
                     yi(4) = gasi8(nt,n+3,nh,nw,ng)
                     call lagrange(x,xi,yi,ans)
                     gasi(nt,m,nh,nw,ng) = 10.0**ans
                  end do 
                  
                  do n=2,L_NPREF-2
                     do m=1,5
                        i     = (n-1)*5+m
                        x     = pfgasref(i)
                        xi(1) = pgasref(n-1)
                        xi(2) = pgasref(n)
                        xi(3) = pgasref(n+1)
                        xi(4) = pgasref(n+2)
                        yi(1) = gasi8(nt,n-1,nh,nw,ng)
                        yi(2) = gasi8(nt,n,nh,nw,ng)
                        yi(3) = gasi8(nt,n+1,nh,nw,ng)
                        yi(4) = gasi8(nt,n+2,nh,nw,ng)
                        call lagrange(x,xi,yi,ans)
                        gasi(nt,i,nh,nw,ng) = 10.0**ans
                     end do 
                  end do

!     Now, get the last interval

                  n = L_NPREF-1                 
                  do m=1,5
                     i     = (n-1)*5+m
                     x     = pfgasref(i)
                     xi(1) = pgasref(n-2)
                     xi(2) = pgasref(n-1)
                     xi(3) = pgasref(n)
                     xi(4) = pgasref(n+1)
                     yi(1) = gasi8(nt,n-2,nh,nw,ng)
                     yi(2) = gasi8(nt,n-1,nh,nw,ng)
                     yi(3) = gasi8(nt,n,nh,nw,ng)
                     yi(4) = gasi8(nt,n+1,nh,nw,ng)
                     call lagrange(x,xi,yi,ans)
                     gasi(nt,i,nh,nw,ng) = 10.0**ans
                  end do  

!     Fill the last pressure point

                  gasi(nt,L_PINT,nh,nw,ng) = &
                       10.0**gasi8(nt,L_NPREF,nh,nw,ng)

               end do
            end do
         end do
      end do

!     Interpolate the values:  now the shortwave

      do nt=1,L_NTREF
         do nh=1,L_REFVAR
            do nw=1,L_NSPECTV
               do ng=1,L_NGAUSS

!     First, the initial interval

                  n = 1 
                  do m=1,5
                     x     = pfgasref(m)
                     xi(1) = pgasref(n)
                     xi(2) = pgasref(n+1)
                     xi(3) = pgasref(n+2)
                     xi(4) = pgasref(n+3)
                     yi(1) = gasv8(nt,n,nh,nw,ng)
                     yi(2) = gasv8(nt,n+1,nh,nw,ng)
                     yi(3) = gasv8(nt,n+2,nh,nw,ng)
                     yi(4) = gasv8(nt,n+3,nh,nw,ng)
                     call lagrange(x,xi,yi,ans)
                     gasv(nt,m,nh,nw,ng) = 10.0**ans
                  end do 
                  
                  do n=2,L_NPREF-2
                     do m=1,5
                        i     = (n-1)*5+m
                        x     = pfgasref(i)
                        xi(1) = pgasref(n-1)
                        xi(2) = pgasref(n)
                        xi(3) = pgasref(n+1)
                        xi(4) = pgasref(n+2)
                        yi(1) = gasv8(nt,n-1,nh,nw,ng)
                        yi(2) = gasv8(nt,n,nh,nw,ng)
                        yi(3) = gasv8(nt,n+1,nh,nw,ng)
                        yi(4) = gasv8(nt,n+2,nh,nw,ng)
                        call lagrange(x,xi,yi,ans)
                        gasv(nt,i,nh,nw,ng) = 10.0**ans
                     end do 
                  end do

!     Now, get the last interval

                  n = L_NPREF-1
                  do m=1,5
                     i     = (n-1)*5+m
                     x     = pfgasref(i)
                     xi(1) = pgasref(n-2)
                     xi(2) = pgasref(n-1)
                     xi(3) = pgasref(n)
                     xi(4) = pgasref(n+1)
                     yi(1) = gasv8(nt,n-2,nh,nw,ng)
                     yi(2) = gasv8(nt,n-1,nh,nw,ng)
                     yi(3) = gasv8(nt,n,nh,nw,ng)
                     yi(4) = gasv8(nt,n+1,nh,nw,ng)
                     call lagrange(x,xi,yi,ans)
                     gasv(nt,i,nh,nw,ng) = 10.0**ans
                  end do  

!     Fill the last pressure point

                  gasv(nt,L_PINT,nh,nw,ng) = &
                      10.0**gasv8(nt,L_NPREF,nh,nw,ng)
                  
               end do
            end do
         end do
      end do
      
      if (corrk_recombin) then
        gasi_recomb(:,:,:,:) = kappa_IR ! non-zero init
        gasv_recomb(:,:,:,:) = kappa_VI ! non-zero init
      endif

!=======================================================================
!     Initialise the continuum absorption data
      if(continuum)then
      do igas=1,ngasmx

         if (igas .eq. igas_N2) then

            dummy = -9999
            call interpolateN2N2(100.D+0,250.D+0,17500.D+0,testcont,.true.,dummy)

         elseif (igas .eq. igas_H2) then

            ! first do self-induced absorption
            dummy = -9999
            call interpolateH2H2(500.D+0,250.D+0,17500.D+0,testcont,.true.,dummy)
            ! then cross-interactions with other gases
            do jgas=1,ngasmx
               if (jgas .eq. igas_N2) then
                  dummy = -9999
                  call interpolateN2H2(592.D+0,278.15D+0,200000.D+0,10000.D+0,testcont,.true.,dummy)
               endif
            enddo

         elseif (igas .eq. igas_CH4) then

            ! first do self-induced absorption
            dummy = -9999
            call interpolateCH4CH4(200.D+0,200.D+0,7500.D+0,testcont,.true.,dummy)
            ! then cross-interactions with other gases
            do jgas=1,ngasmx
               if (jgas .eq. igas_N2) then
                  dummy = -9999
                  call interpolateN2CH4(200.D+0,250.0D+0,100000.D+0,5000.D+0,testcont,.true.,dummy)
               endif
            enddo

         endif  
         
      enddo
      endif

      print*,'----------------------------------------------------'
      print*,'And that`s all we have. It`s possible that other'
      print*,'continuum absorption may be present, but if it is we'
      print*,'don`t yet have data for it...'
      print*,''

!     Deallocate local arrays
!$OMP BARRIER
!$OMP MASTER
      IF( ALLOCATED( gasi8 ) ) DEALLOCATE( gasi8 )
      IF( ALLOCATED( gasv8 ) ) DEALLOCATE( gasv8 )
      IF( ALLOCATED( pgasref ) ) DEALLOCATE( pgasref )
      IF( ALLOCATED( gastype ) ) DEALLOCATE( gastype )
!$OMP END MASTER
!$OMP BARRIER

      return
    end subroutine sugas_corrk
