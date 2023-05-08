      program generate_kmatrix
      implicit none

!     -------------------------------------------------------------
!     Purpose: to read hires spectra and produce kdistribution data
!     Authour: Robin Wordsworth (2010)
!     -------------------------------------------------------------

      integer nmolec,mol
      integer i,j,k,im
      integer Nb

      character*4 fname
      character*32 fnum

      integer iT, iP, iQ, Nlevs

      double precision nuT, sigT, kT

      integer Nq
      parameter (Nq=16) ! we choose this here

      integer Nmax
      parameter (Nmax=1000) ! maximum number of levels

      integer Nmol_max
      parameter (Nmol_max=10) ! maximum number of radiatively active species

      integer Nkmax
      ! large array = slow calculation... this makes a huge difference
      !parameter (Nkmax=50000) ! works for IR, 27 bands
      !parameter (Nkmax=100000) ! works for IR, 27 bands
      !parameter (Nkmax=140000) ! works for IR, 32 bands
      !parameter (Nkmax=250000) ! works for VI, 36 bands
      parameter (Nkmax=400000) ! works for IR, 8 bands
      !parameter (Nkmax=400000) ! works for VI, 12 bands
      !parameter (Nkmax=200000) ! works for IR, 30 bands

      integer NbMAX
      parameter (NbMAX=100)

      integer Ngres
      parameter (Ngres=200)

      character*10 molec_names(1:Nmol_max)
      double precision T(1:Nmax), P(1:Nmax), Q(1:Nmax)


      double precision pi
      parameter (pi=3.14159265)
      double precision x
      double precision y(1:Nq+2)

      double precision nu_min(1:NbMAX) ! make allocatable?
      double precision nu_max(1:NbMAX)

      integer quad, band

      double precision ktemp(1:Nq)
      double precision w(1:Nq), gw(1:Nq)

      double precision, allocatable :: kdist_temp(:)
      double precision, allocatable :: kdist(:,:,:,:,:)

      double precision nuband(1:Nkmax), kband(1:Nkmax), nuwband(1:Nkmax) ! better static I think
      double precision nuwbtot

      double precision gfine(1:Ngres)      
      double precision klevs(1:Ngres)
      double precision addline(1:Nkmax)

      double precision kmax, kmin, lkmax, lkmin, dlogK, klev
      integer ic, ib, ig, ik1, ik2

      integer filerr
      logical filexi

      logical addCIA
      integer iCO2
      double precision k1,k2,kCIA,amg,Pp

      real truemean, distmean, etot
      

      addCIA=.false.

      print*,'Creating the GCM correlated-k data...'
      print*,' '

      !--------------------------------------------------------
      ! load temperature, pressure, mixing ratio data
      open(9,file='T.dat')
      read(9,*) iT
      do i=1,iT
         read(9,*) T(i)
      enddo
      close(9)

      open(10,file='p.dat')
      read(10,*) iP
      do i=1,iP
         read(10,*) P(i)
         P(i)=10.**P(i)  ! note we keep P in millibars for the conversion later
      enddo
      close(10)

      open(11,file='Q.dat')
      read(11,*) nmolec
      do mol=1,nmolec
         read(11,*) molec_names(mol)
      enddo
      read(11,*) iQ
      do i=1,iQ
         read(11,*) Q(i)
      enddo
      close(11)

      Nlevs=iT*iP*iQ ! total number of layers

      print*,'Temperature layers:  ',iT
      print*,'Pressure layers:     ',iP
      print*,'Mixing ratio layers: ',iQ
      print*,'Total:               ',Nlevs


      if(addCIA)then

         if(nmolec.gt.2)then
            print*,'I cant deal with this situation yet'
            call abort
         endif

         iCO2=0
         do mol=1,nmolec
            if(molec_names(mol).eq.'CO2')then
               iCO2=mol
            endif
         enddo

         !if(iCO2.eq.1)then
         print*,' '
         print*,'--- Adding CO2 CIA ---'
         !else
         !   print*,'Warning, not CO2+H2O mixture, need CO2 mixing ratio'
         !   call abort
         !endif
      endif

      !--------------------------------------------------------
      ! load spectral data
      open(12,file='narrowbands.in')
      read(12,*) Nb

      do band=1,Nb
         read(12,*) nu_min(band),nu_max(band)
      enddo
      close(12)

      print*,' '
      print*,'Number of spectral bands: ',Nb

      print*,'Band limits:'
      do band=1,Nb
         print*,band,'-->',nu_min(band),nu_max(band),' cm^-1'
      end do

      !--------------------------------------------------------
      ! create g-space data
      do quad=0,Nq+1
         x=dble(quad)/dble(Nq+2)
         !y(quad+1)=sin(pi*x)*exp(-5*x)
         !y(quad+1)=sin(pi*x)*exp(-9*x)  ! CO2_H2Ovar 38x36
         y(quad+1)=sin(pi*x)*exp(-12*x)
         !y(quad+1)=sin(pi*x)*exp(-15*x) ! CO2_CH4
      enddo

      do quad=1,Nq
         w(quad)=y(quad+1)
      enddo
      w=w/sum(w)

      print*,' '
      print*,'Number of g-space intervals: ',Nq

      print*,'g-space weighting:'
      do quad=1,Nq
         print*,quad,'-->',w(quad)
      end do

      ! write g-space data
      open(14,file='g.dat')
      write(14,*) Nq+1

      do quad=1,Nq
         write(14,*) w(quad)
      enddo
      write(14,*) 0.0
      close(14)

      ! cumulative sum for g-vector
      gw(1)=w(1)
      do quad=2,Nq
         gw(quad)=gw(quad-1)+w(quad)
      enddo

      ! allocate 5D matrices
      allocate(kdist_temp(Nq))
      allocate(kdist(iT,iP,iQ,Nb,Nq+1))


      nuT=0.0
      do i=1,iQ
         do j=1,iP


            if(addCIA.and.(iCO2.eq.1))then 
               Pp=P(j)*(1.0-Q(i))
            else
               Pp=P(j)*Q(i)
            endif
            !if(addCIA) Pp=P(j)*(1.0-Q(i))


            do k=1,iT             
               im = (i-1)*iP*iT + (j-1)*iT + k

               !-----------------------------------------------
               ! open kspectrum file
               write(fnum,*) im
               if(im<10)then
                   fname='k00'//trim(adjustl(fnum))
               elseif(im<100)then
                   fname='k0'//trim(adjustl(fnum))
               else
                   fname='k'//trim(adjustl(fnum))
               endif

               ! check that the file exists
               inquire(FILE='hires_spectrum/'//fname,EXIST=filexi)
               if(.not.filexi) then
                  write(*,*)'The file ',fname(1:LEN_TRIM(fname))
                  write(*,*)'was not found in the folder hires_spectrum, exiting.'
                  call abort
               endif

               print*,'Opening file ',fname(1:LEN_TRIM(fname))
               open(17,file='hires_spectrum/'//fname)

               nuT=0.0
               do band=1,Nb
               print*,'Band ---> ',band

                  !--------------------------------------------
                  ! write kdata to temporary band arrays
                  ! and find max / min k-values in band
                  do while (nuT.lt.nu_min(band))
                      read(17,*,IOSTAT=filerr) nuT, sigT, kT
                   ! ideally need to catch all errors here (ask Ehouarn)
                      if(filerr.lt.0)then
                         print*,    &
                             'nu_min greater than highest value in file'
                         print*,'IOSTAT=',filerr
                         print*,'nu_min=',nu_min(band)
                         print*,'nu_endoffile=',nuT
                         call abort
                      endif
                  enddo

                  ic=1
                  kmin=kT
                  kmax=kT
                  nuband(:)=0
                  kband(:)=0

                  do while (nuT.lt.nu_max(band))

                      nuband(ic)=nuT
                      kband(ic)=kT

                      ! insert CIA option here
                      if(addCIA)then

                         k1=0.0
                         k2=0.0
                         if(nuT.lt.500.0) call getspc(T(k),nuT,k1)
                         if((nuT.gt.1000.0).and.(nuT.lt.1800.0)) call baranov(T(k),nuT,k2)
                         ! cm^-1.amagat^2

                         amg=273.15D+0/T(k)*(Pp/1013.25) ! amagats of CO2

                         kCIA=(k1+k2)*1.0D+2*amg**2.0D+0 ! m^-1
                         kband(ic)=kband(ic)+kCIA
                   	 kT=kband(ic)

                         !kband(ic)=kCIA
			 !kT=kCIA
			 !if(nuband(ic).gt.1250.0)then
			 !	print*,'P=',P(i)
			 !	print*,'T=',T(k)
		         !	print*,'amg=',amg
			 !	print*,'nu=',nuband(ic)
                         !       print*,'k1=',k1
                         !       print*,'k2=',k2
                         !       print*,'kCIA=',kCIA
		  	 !	call abort
                         !endif

                      endif

                      if(kT.lt.kmin)then
                        kmin=kT
                      endif
                      if(kT.gt.kmax)then ! bug corrected!
                        kmax=kT
                      endif

                      read(17,*,IOSTAT=filerr) nuT, sigT, kT
                      if(filerr.lt.0)then
                         print*,    &
                             'nu_max greater than highest value in file'
                         print*,'IOSTAT=',filerr
                         print*,'nu_max=',nu_max(band)
                         print*,'nu_endoffile=',nuT
                         call abort
                      endif


                      ic=ic+1
                      
                      if(ic.gt.Nkmax)then
                         print*,'Spectral resolution too low for these bands, exiting!'
                         print*,'Increase Nkmax in generate_kmatrix.F90.'
                         call abort
                      endif

                  enddo
                  nuband(ic)=nuT
                  kband(ic)=kT

                  if(kmax.le.0.0)then
                  !if(kmax.le.1.0e-15)then
                     print*,'kmax=',kmax
                     print*,'Setting values to zero in this band.'
                     kdist(k,j,i,band,1:(Nq+1))=0.0
                  else

                  if(kmin.lt.1e-30)then
                     kmin=1.0e-30
                  endif

                  !--------------------------------------------
                  ! calculate delta nu values
                  nuwband(:)=0.0
                  do ib=2,(ic-1)
                      nuwband(ib)=(nuband(ib+1)-nuband(ib-1))/2
                  enddo
                  nuwband(1)=nuwband(2)
                  nuwband(ic)=nuwband(ic-1) ! causes tiny boundary error; whatever
                  nuwbtot=sum(nuwband)

                  lkmax=log10(kmax)
                  lkmin=log10(kmin)
                  
                  dlogK=(lkmax-lkmin)/(Ngres-1)                 

                  !--------------------------------------------
                  ! calculate the g function
                  gfine(:)=0.0
                  do ig=1,Ngres
                    klev=lkmin + dble(ig)*dlogK
                    klev=10**klev
                    klevs(ig)=klev

                    addline(:)=0
                    do ib=1,ic
                        if(kband(ib).lt.klev)then
                           addline(ib)=1.0
                        endif
                    enddo
                    gfine(ig)=sum(addline(:)*nuwband(:))/nuwbtot

                  enddo
                  print*,'gnorm    = ',gfine(Ngres)

                  !--------------------------------------------
                  ! invert the g function
                  ik1=1
                  ik2=1
                  kdist_temp(:)=0.0
                  do quad=1,Nq
                    
                    ik2=ik1

                    gfind: do while (gfine(ik2).lt.gw(quad))
                       ik2=ik2+1
                       if(ik2.gt.Ngres)then
                          print*,'g-inversion overflow...'
                          ik2=Ngres
                          exit gfind
                       endif
                    enddo gfind

                    !if(ik2.gt.Ngres)then ! avoid out-of-bounds on final step
                    !   ik2=Ngres
                    !endif


                     if(.not.(ik1.ge.ik2))then ! this could probably be improved
                        !kdist_temp(quad)=sum(klevs(ik1:ik2))/(ik2-ik1) ! a bug was here!
                        kdist_temp(quad)=sum(klevs(ik1:ik2))/(ik2+1-ik1)


                     elseif(quad.ne.1)then
                        kdist_temp(quad)=kdist_temp(quad-1)
                     else
                        kdist_temp(quad)=kmin
                     endif

                     !-----------------------------------------
                     ! assign to 5D matrix and convert to cm^2 / molecule
                     kdist(k,j,i,band,quad)=1.3807e-21 *               &
                                kdist_temp(quad) * T(k)/P(j)

                     ik1=ik2
                  enddo ! g-space
                  kdist(k,j,i,band,Nq+1)=0.0

                  truemean = sum(kband(1:ic)*nuwband(1:ic))/sum(nuwband(1:ic))
                  distmean = sum(kdist_temp(:)*w(:))/sum(w)
                  etot     = 100*abs(truemean-distmean)/(distmean+truemean)
                  print*,'truemean = ',truemean
                  print*,'distmean = ',distmean
                  !print*,'kdisttemp = ',kdist_temp(Nq)
                  !print*,'kdistmax  = ',kmax
                  print*,'etot (%) = ',etot
                  !print*,'w  = ',w

                  endif ! if kmax .le. 0.0 

               enddo ! bands

            enddo
         enddo
      enddo

      !--------------------------------------------------------
      ! save data and clean up
      close(17)
      open(14,file='corrk_gcm.dat',form='formatted')
      write(14,*) kdist
      close(14)
      print*,' '
      print*,'Correlated-k data saved in corrk_gcm.dat'
      print*,' '
      print*,'Now you must change the values of L_NTREF etc. in radinc_h.F90 (if necessary)'
      print*,'and the variable "corrkdata" in callphys.def'

      deallocate(kdist_temp)
      deallocate(kdist)

    end program generate_kmatrix
