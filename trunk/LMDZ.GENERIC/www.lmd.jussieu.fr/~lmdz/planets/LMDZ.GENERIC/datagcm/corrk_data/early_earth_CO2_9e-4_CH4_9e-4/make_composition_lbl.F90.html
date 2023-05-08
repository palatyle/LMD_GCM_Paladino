      program make_composition_gen
      implicit none

!     -------------------------------------------------------------   
!     Purpose: to create "composition.in" file readable by kspectrum
!     Authour: Robin Wordsworth (2010)
!     -------------------------------------------------------------

      !include 'max.inc'
      include 'formats.inc'
      integer strlen,nmolec, mol!,nlev_mf,mol,imf,ind,indf
      integer i,j,k,im,m
      integer Nlevs
      integer nb,na,ne
      

      integer Nmax
      parameter(Nmax=25)

      integer Nmol_max
      parameter(Nmol_max=5)

      double precision alt(1:Nmax)
      double precision pres(1:Nmax)
      double precision temp(1:Nmax)

      double precision P_mf(1:Nmax)
      double precision x(1:Nmax,1:Nmol_max)
      double precision xlev(1:Nmol_max)
      character*100 composition_file,model_file
      character*100 planet_descriptor
      character*200 string
      character*20 s
      character*10 molec_names(1:Nmol_max),molname
      character*1 spch

      integer iP, iT, iQ
      double precision P(1:Nmax), T(1:Nmax), Q(1:Nmax)

      double precision Hscale

      spch=' '

      print*,'Name of atmosphere / planet:'
      read*,planet_descriptor

      composition_file='./composition.in'

!     ! nmolec=2 ! total number of molecules

      T(:)=0.0
      P(:)=0.0
      Q(:)=0.0

      ! load general data     
      open(8,file='general.dat')
      read(8,*) Hscale
      close(8)

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
         P(i)=P(i)/101325.0
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


      if(iT.ne.iP .or. iT.ne.iQ)then
         print*,'Number of vertical levels inconsistent, exiting.'
         call abort
      else
         Nlevs=iT
      endif

      print*,'Vertical layers:  ',Nlevs
   
      if(planet_descriptor(1:5).eq.'Earth')then
         print*,'Assuming CO2 mixing ratio = 5.88x10^-4.'
      endif
      



      alt(:)=0.0
      pres(:)=0.0
      temp(:)=0.0
      x(:,:)=0.0

      do im=1,Nlevs

         if(im.gt.1)then
            alt(im) = -(Hscale/1000.0)*log(P(im)/P(1)) 
         endif

         pres(im)=P(im)
         temp(im)=T(im)

         if(nmolec.eq.1)then
            x(im,1)=Q(im)
         elseif(nmolec.eq.2)then
            x(im,1)=1-Q(im)
            
            if(planet_descriptor(1:5).eq.'Earth')then
               x(im,1)=5.88e-4
            endif
            
            x(im,2)=Q(im)
         else
            print*,'Dont know how this works...?'
            stop
         endif

      enddo

      
      open(12,file=composition_file(1:strlen(composition_file)))
      string='Atmospheric composition input '
      string=string(1:strlen(string))//' data file for planet:'
      string=string(1:strlen(string))//' '//                       &
          planet_descriptor(1:strlen(planet_descriptor))
      write(12,'(a)') string(1:strlen(string))
      write(12,39) 'Number of atmospheric levels: ',Nlevs
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
         
         string=string(1:strlen(string)-1)     &
             //s(1:strlen(s)-1)                & 
              //'x['                           & 
              //molname(1:strlen(molname))     &
              //']*'
      enddo
      write(12,*) string(1:strlen(string)-1)

      ! write data to file, level by level
      do im=1,Nlevs
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






