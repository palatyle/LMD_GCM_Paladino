










!
! $Id: write_field.F90 1279 2009-12-10 09:02:56Z fairhead $
!
module write_field
implicit none

  integer, parameter :: MaxWriteField = 100
  integer, dimension(MaxWriteField),save :: FieldId
  integer, dimension(MaxWriteField),save :: FieldVarId
  integer, dimension(MaxWriteField),save :: FieldIndex
  character(len=255), dimension(MaxWriteField) ::  FieldName 
   
  integer,save :: NbField = 0
  
  interface WriteField
    module procedure WriteField3d,WriteField2d,WriteField1d
  end interface WriteField
  contains
  
    function GetFieldIndex(name)
    implicit none
      integer          :: GetFieldindex
      character(len=*) :: name
    
      character(len=255) :: TrueName
      integer            :: i
       
      
      TrueName=TRIM(ADJUSTL(name))
    
      GetFieldIndex=-1
      do i=1,NbField
        if (TrueName==FieldName(i)) then
          GetFieldIndex=i
          exit
        endif
      enddo
    end function GetFieldIndex
 
    subroutine WriteField3d(name,Field)
    implicit none
      character(len=*) :: name
      real, dimension(:,:,:) :: Field 
      integer, dimension(3) :: Dim
      
      Dim=shape(Field)
      call WriteField_gen(name,Field,Dim(1),Dim(2),Dim(3))  
  
    end subroutine WriteField3d
    
    subroutine WriteField2d(name,Field)
    implicit none
      character(len=*) :: name
      real, dimension(:,:) :: Field 
      integer, dimension(2) :: Dim
      
      Dim=shape(Field)
      call WriteField_gen(name,Field,Dim(1),Dim(2),1)  
  
    end subroutine WriteField2d
    
    subroutine WriteField1d(name,Field)
    implicit none
      character(len=*) :: name
      real, dimension(:) :: Field 
      integer, dimension(1) :: Dim
      
      Dim=shape(Field)
      call WriteField_gen(name,Field,Dim(1),1,1)  
  
    end subroutine WriteField1d
        
    subroutine WriteField_gen(name,Field,dimx,dimy,dimz)
    implicit none
    include 'netcdf.inc'
      character(len=*) :: name
      integer :: dimx,dimy,dimz
      real,dimension(dimx,dimy,dimz) :: Field
      integer,dimension(dimx*dimy*dimz) :: ndex
      integer :: status
      integer :: index
      integer :: start(4)
      integer :: count(4)
      
           
      Index=GetFieldIndex(name)
      if (Index==-1) then
        call CreateNewField(name,dimx,dimy,dimz)
	Index=GetFieldIndex(name)
      else
        FieldIndex(Index)=FieldIndex(Index)+1.
      endif
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=FieldIndex(Index)

      count(1)=dimx
      count(2)=dimy
      count(3)=dimz
      count(4)=1

      status = NF_PUT_VARA_DOUBLE(FieldId(Index),FieldVarId(Index),start,count,Field)
      status = NF_SYNC(FieldId(Index))
      
    end subroutine WriteField_gen
       
    subroutine CreateNewField(name,dimx,dimy,dimz)
    implicit none
    include 'netcdf.inc'  
      character(len=*) :: name
      integer :: dimx,dimy,dimz
      integer :: TabDim(4)
      integer :: status
      
      
      NbField=NbField+1
      FieldName(NbField)=TRIM(ADJUSTL(name))
      FieldIndex(NbField)=1
      
      
      status = NF_CREATE(TRIM(ADJUSTL(name))//'.nc', NF_CLOBBER, FieldId(NbField))
      status = NF_DEF_DIM(FieldId(NbField),'X',dimx,TabDim(1))
      status = NF_DEF_DIM(FieldId(NbField),'Y',dimy,TabDim(2))
      status = NF_DEF_DIM(FieldId(NbField),'Z',dimz,TabDim(3))
      status = NF_DEF_DIM(FieldId(NbField),'iter',NF_UNLIMITED,TabDim(4))
      status = NF_DEF_VAR(FieldId(NbField),FieldName(NbField),NF_DOUBLE,4,TabDim,FieldVarId(NbField))
      status = NF_ENDDEF(FieldId(NbField))

    end subroutine CreateNewField
    
    
    
  subroutine write_field1D(name,Field)
    implicit none
  
    integer, parameter :: MaxDim=1
    character(len=*)   :: name
    real, dimension(:) :: Field
    real, dimension(:),allocatable :: New_Field
    character(len=20) :: str
    integer, dimension(MaxDim) :: Dim
    integer :: i,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize 
    integer :: pos
    character(len=255) :: form
    character(len=255) :: MaxLen
    
    
    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'",//)')
    Dim=shape(Field)
    MaxLen=int2str(len(trim(int2str(Dim(1)))))
    ColumnSize=20+6+3+len(trim(int2str(Dim(1))))
    Nb=0
    Pos=2
    do i=1,Dim(1)
      nb=nb+1
      
      if (MOD(nb,NbCol)==0) then
        form='(t'//trim(int2str(pos))// ',i'//trim(MaxLen) //'," ---> ",g22.16,/)'
        Pos=2
      else
        form='(t'//trim(int2str(pos))// ',i'//trim(MaxLen) //'," ---> ",g22.16," | ",)'
        Pos=Pos+ColumnSize
      endif
      write (id,form,advance='no') i,Field(i)
    enddo
     
    close(id)

  end subroutine write_field1D

  subroutine write_field2D(name,Field)
    implicit none
  
    integer, parameter :: MaxDim=2
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    real, dimension(:,:),allocatable :: New_Field
    character(len=20) :: str
    integer, dimension(MaxDim) :: Dim
    integer :: i,j,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize 
    integer :: pos,offset
    character(len=255) :: form
    character(len=255) :: spacing
    
    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'",//)')
    
    Dim=shape(Field)
    offset=len(trim(int2str(Dim(1))))+len(trim(int2str(Dim(2))))+3
    ColumnSize=20+6+3+offset

    spacing='(t2,"'//repeat('-',ColumnSize*NbCol)//'")'
    
    do i=1,Dim(2)
      nb=0
      Pos=2
      do j=1,Dim(1)
        nb=nb+1
      
        if (MOD(nb,NbCol)==0) then
          form='(t'//trim(int2str(pos))//            &
               ',"('//trim(int2str(j))//','          &
                    //trim(int2str(i))//')",t'       & 
                    //trim(int2str(pos+offset))     &    
                    //'," ---> ",g22.16,/)'
          Pos=2
        else
          form='(t'//trim(int2str(pos))//            &
               ',"('//trim(int2str(j))//','          &
                    //trim(int2str(i))//')",t'       & 
                    //trim(int2str(pos+offset))     &    
                    //'," ---> ",g22.16," | ")'
          Pos=Pos+ColumnSize
        endif
        write (id,form,advance='no') Field(j,i)
      enddo
      if (MOD(nb,NbCol)==0) then
        write (id,spacing)
      else
        write (id,'("")')
        write (id,spacing)
      endif
    enddo
     
  end subroutine write_field2D
  
  subroutine write_field3D(name,Field)
    implicit none
  
    integer, parameter :: MaxDim=3
    character(len=*)   :: name
    real, dimension(:,:,:) :: Field
    real, dimension(:,:,:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    integer :: i,j,k,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize 
    integer :: pos,offset
    character(len=255) :: form
    character(len=255) :: spacing

    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'"//)')
    
    Dim=shape(Field)
    offset=len(trim(int2str(Dim(1))))+len(trim(int2str(Dim(2))))+len(trim(int2str(Dim(3))))+4
    ColumnSize=22+6+3+offset

!    open(unit=id,file=name,form=formatted
   
    spacing='(t2,"'//repeat('-',ColumnSize*NbCol)//'")'
    
    do i=1,Dim(3)
    
      do j=1,Dim(2)
        nb=0
        Pos=2
        
        do k=1,Dim(1)
        nb=nb+1
      
          if (MOD(nb,NbCol)==0) then
            form='(t'//trim(int2str(pos))//            &
                 ',"('//trim(int2str(k))//','          &
                      //trim(int2str(j))//','          &
                      //trim(int2str(i))//')",t'       & 
                      //trim(int2str(pos+offset))      &    
                      //'," ---> ",g22.16,/)'
           Pos=2
          else
            form='(t'//trim(int2str(pos))//            &
                 ',"('//trim(int2str(k))//','          &
                      //trim(int2str(j))//','          &
                      //trim(int2str(i))//')",t'       & 
                      //trim(int2str(pos+offset))      &    
                      //'," ---> ",g22.16," | ")'
! dépent de l'implémention, sur compaq, c'est necessaire
!            Pos=Pos+ColumnSize
          endif
          write (id,form,advance='no') Field(k,j,i)
        enddo
        if (MOD(nb,NbCol)==0) then
          write (id,spacing)
        else
          write (id,'("")')
          write (id,spacing)
        endif
      enddo
      write (id,spacing)
    enddo
    
    close(id)
  
  end subroutine write_field3D  
  
  function int2str(int)
    implicit none
    integer, parameter :: MaxLen=10
    integer,intent(in) :: int
    character(len=MaxLen) :: int2str
    logical :: flag
    integer :: i
    flag=.true.
    
    i=int
    
    int2str=''
    do while (flag)
      int2str=CHAR(MOD(i,10)+48)//int2str
      i=i/10
      if (i==0) flag=.false.
    enddo
  end function int2str

end module write_field
  
