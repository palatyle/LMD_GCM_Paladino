module write_field_p
implicit none
  
  interface WriteField_p
    module procedure Write_field3d_p,Write_Field2d_p,Write_Field1d_p
  end interface WriteField_p
  
  contains
  
  subroutine write_field1D_p(name,Field)
    USE parallel_lmdz
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=1
    character(len=*)   :: name
    real, dimension(:) :: Field
    real, dimension(:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1)))
    New_Field(:)=Field(:)
    call Gather_Field(New_Field,dim(1),1,0)
    
    if (MPI_Rank==0) call WriteField(name,New_Field)
    
    end subroutine write_field1D_p

  subroutine write_field2D_p(name,Field)
    USE parallel_lmdz
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=2
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    real, dimension(:,:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1),Dim(2)))
    New_Field(:,:)=Field(:,:)
    call Gather_Field(New_Field(1,1),dim(1)*dim(2),1,0)
    
    if (MPI_Rank==0) call WriteField(name,New_Field)
    
     
  end subroutine write_field2D_p
  
  subroutine write_field3D_p(name,Field)
    USE parallel_lmdz
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=3
    character(len=*)   :: name
    real, dimension(:,:,:) :: Field
    real, dimension(:,:,:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1),Dim(2),Dim(3)))
    New_Field(:,:,:)=Field(:,:,:)
    call Gather_Field(New_Field(1,1,1),dim(1)*dim(2),dim(3),0)
    
   if (MPI_Rank==0) call WriteField(name,New_Field)
    
  end subroutine write_field3D_p  

end module write_field_p
  
