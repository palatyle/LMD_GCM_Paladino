      program make_composition_gen
      implicit none

!     -------------------------------------------------------------   
!     Purpose: to create "composition.in" file readable by kspectrum
!     Authour: Robin Wordsworth (2010)
!     -------------------------------------------------------------

      include 'max.inc'
      include 'formats.inc'
      integer strlen,nmolec,nlev_mf,mol,imf,ind,indf
      integer i,j,k,im,m
      integer nb,na,ne
      double precision alt(1:Nmax)
      double precision pres(1:Nmax)
      double precision temp(1:Nmax)

      double precision P_mf(1:Nmax)
      double precision x(1:Nmax,1:Nmol_max)
      double precision xlev(1:Nmol_max)
      double precision xconst(1:Nmol_max)
      character*100 composition_file,model_file
      character*100 planet_descriptor
      character*200 string
      character*20 s
      character*10 molec_names(1:Nmol_max),molname
      character*1 spch

      integer iP, iT, iQ
      double precision P(1:Nmax), T(1:Nmax), Q(1:Nmax)

      spch=' '

      print*,'Name of atmosphere / planet:'
      read*,planet_descriptor

      composition_file='./composition.in'

!     ! nmolec=2 ! total number of molecules

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
         P(i)=10.**P(i)/1013.25
      enddo
      close(10)

      open(11,file='Q.dat')
      read(11,*) nmolec

      print*,'nmolec=',nmolec

      do mol=1,nmolec
         read(11,*) molec_names(mol)
      enddo
      read(11,*) iQ
      do i=1,iQ
         read(11,*) Q(i)
      enddo
      close(11)

      m=iP*iT*iQ ! total number of layers

      print*,'Temperature layers:  ',iT
      print*,'Pressure layers:     ',iP
      print*,'Mixing ratio layers: ',iQ
      print*,'Total:               ',m

      do mol=1,nmolec-1
         print*,'Please enter vmr of ',molec_names(mol)
         read(*,*) xconst(mol)
      end do


!      if(planet_descriptor(1:5).eq.'Earth')then
!         print*,'Assuming CO2 mixing ratio = 5.88x10^-4.' ! this is MASS mixing ratio. ERROR!
!      endif
     
      do i=1,iQ
         do j=1,iP              
            do k=1,iT

               im = (i-1)*iP*iT + (j-1)*iT + k

               alt(im)=0.D0
               pres(im)=P(j)
               temp(im)=T(k)

             !  if(nmolec.eq.1)then
             !     x(im,1)=Q(i)
             !  elseif(nmolec.eq.2)then
               do mol=1,nmolec-1
                  x(im,mol)=xconst(mol)*(1.0-Q(i))
               enddo
               x(im,nmolec)=Q(i)

! not finished!!!!

!                  if(planet_descriptor(1:5).eq.'Earth')then
!                     x(im,1)=5.88e-4
!                  endif
!
!                  x(im,nmolec)=Q(i)
!               else
!                  print*,'nmolec=',nmolec
!                  print*,'Dont know how this works...?'
!                  stop
!               endif

            enddo
         enddo
      enddo

      open(12,file=composition_file(1:strlen(composition_file)))
      string='Atmospheric composition input '
      string=string(1:strlen(string))//' data file for planet:'
      string=string(1:strlen(string))//' '//

     &     planet_descriptor(1:strlen(planet_descriptor))
      write(12,'(a)') string(1:strlen(string))
      write(12,39) 'Number of atmospheric levels: ',m
      write(12,34) 'Number of molecules: ',nmolec
      write(12,*) 
      string='      z (km)     /    P (atm)     /     T (K) *'


      do mol=1,nmolec

         molname=molec_names(mol)

         nb=5
         if (strlen(molname).eq.1) then
            na=7
         else if (strlen(molname).eq.2) then
            na=6
         else if (strlen(molname).eq.3) then
            na=5
         endif
         s='*'
         do ne=1,nb
            s=s(1:strlen(s)-1)//spch//'*'
         enddo
         s=s(1:strlen(s)-1)//'/*'
         do ne=1,na
            s=s(1:strlen(s)-1)//spch//'*'
         enddo
         
         string=string(1:strlen(string)-1)
     &        //s(1:strlen(s)-1)
     &        //'x['
     &        //molname(1:strlen(molname))
     &        //']*'
      enddo
      write(12,*) string(1:strlen(string)-1)

      do im=1,m

         do mol=1,nmolec
            xlev(mol)=x(im,mol)
         enddo                  ! mol

         write(12,50) alt(im),pres(im),temp(im),(xlev(mol),mol=1,nmolec)

      enddo

      write(12,*) 
      close(12)

      write(*,*) 'Output file successfully generated:'
      write(*,*) composition_file(1:strlen(composition_file))
      
      end
