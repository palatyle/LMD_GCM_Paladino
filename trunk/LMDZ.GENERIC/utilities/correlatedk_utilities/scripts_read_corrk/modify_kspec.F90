      program mod_kmatrix
      implicit none

      integer i,j,k
      integer Nb_IR, Nb_VI, N_g

      integer iT, iP, iQ

      double precision pi
      parameter (pi=3.14159265)

      double precision, allocatable :: kdist_IR(:,:,:,:,:)
      double precision, allocatable :: kdist_VI(:,:,:,:,:)      

      iT=10 # Number of temperature points
      iP=9 # Number of pressure points
      iQ=11 # Number of mixing ratio points

      Nb_IR=38 # Number of IR band
      Nb_VI=36 # Number of VI band
      
      N_g=17 # Number of gauss points

      allocate(kdist_IR(iT,iP,iQ,Nb_IR,N_g))
      allocate(kdist_VI(iT,iP,iQ,Nb_VI,N_g))
      
      open(1,file='corrk_gcm_IR.dat',form='formatted')
      read(1,*) kdist_IR
      close(1)

      open(2,file='corrk_gcm_VI.dat',form='formatted')
      read(2,*) kdist_VI
      close(2)

      do i=1,iQ
         do j=1,iP
            do k=1,N_g
               kdist_IR(1,j,i,1:Nb_IR,k)=kdist_IR(3,j,i,1:Nb_IR,k) # Here we manipulate the corrk IR file
               kdist_VI(1,j,i,1:Nb_VI,k)=kdist_VI(3,j,i,1:Nb_VI,k) # Here we manipulate the corrk VI file
            enddo
         enddo
      enddo

      open(11,file='corrk_gcm_IR_mod.dat',form='formatted')
      write(11,*) kdist_IR
      close(11)

      open(12,file='corrk_gcm_VI_mod.dat',form='formatted')
      write(12,*) kdist_VI
      close(12)

      deallocate(kdist_IR)
      deallocate(kdist_VI)

    end program mod_kmatrix
