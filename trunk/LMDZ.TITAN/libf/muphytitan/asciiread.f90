! Copyright Jérémie Burgalat (2010-2015,2017)
! 
! jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to provide configuration 
! file and command line arguments parsing features to Fortran programs.
! 
! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: asciiread.f90
!! summary: ASCII data file reader source file
!! author: burgalat
!! date: 2013-2015,2017
MODULE ASCIIREAD
  !! ASCII data file reader module
  !!
  !! This module provides a single generic method that can be used to read 2D/3D
  !! data array from ASCII file.
  !!
  !! ``` 
  !! FUNCTION read_data(path,data) RESULT(err)
  !! ```
  !!
  !! Where:
  !!
  !! - __path__ is a string with the path of the data file.
  !! - __data__ is an output __allocatable__ 2D/3D array of real(kind=8) values.
  !!
  !! ## Expected Format of the data file
  !!
  !! The input file:
  !!   - must use blank space(s) as value delimiter.
  !!   - must have a regular number of columns, that is each data line must
  !!     have the same number of columns. 
  !!   - can contains any number of empty lines and/or comment line (i.e. line
  !!     where first non-blank character is "#"). All other lines are assumed
  !!     to be data.
  !! Moreover, in the case of 3D data, it must use a SINGLE empty line for "depth" block separator.  
  !!
  !! Error occured when:
  !!
  !! - path does not refer to a existing file (1)
  !! - No free logical unit available (1)
  !! - the file does not have regular data-columns number (5)
  !! - at least a value cannot be cast in double precision (5)
  !!
  !! On error, __data__ array is not allocated.
  !!
  !! ## 3D data structure
  !!
  !! In the case of 3D data the method assumes that the input files consists in _D_ blocks
  !! of _R_ lines with _C_ columns. Each block must be separated by a single empty line
  !! and each columns must separated by one or more blank spaces (no tabulation ALLOWED).
  !!  
  !! On success, the shape of the 3D output array will be _data(R,C,D)_.
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END, IOSTAT_EOR
  !USE STRING_OP, ONLY: tokenize,from_string, st_slen
  USE STRING_OP
  USE ERRORS
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: noerror,error, error_to_string,aborting
  PUBLIC :: readline,read_data, OPERATOR(/=), OPERATOR(==)

  !! Global interface to reading methods
  INTERFACE read_data
    MODULE PROCEDURE read_data_2d, read_data_3d
  END INTERFACE

  CONTAINS 


  FUNCTION read_data_3d(path,data3d,delimiter) RESULT(err)
    !! Read an input ASCII data file and stores its content in a 3D array
    !!
    !! The function reads an ASCII file and saves its values in a real(kind=8) 3D-array.
    !! 
    !! The input file:
    !!
    !! - must have a regular number of columns, that is each data line must have the same number
    !!   of columns (according to the delimiter used).  
    !! - must use a SINGLE empty line for "depth" block separator.
    !! - can contains any number of comment lines (i.e. line where first non-blank character is "#").
    !!   All other lines (except empty lines) are assumed to be data.
    !! 
    !! Error occured when:
    !!    - Path does not refer to a existing file (-11)
    !!    - No free logical unit available (-1)
    !!    - The file does not have regular data-columns number (-5)
    !!    - At least a value cannot be cast in double precision (-10)
    !!
    !! The method assumes the input files consists in _D_ block of _R_ lines
    !! with _C_ columns. Each block must be separated by a single empty line and
    !! each columns must separated by one or more blank spaces (no tabulation ALLOWED).
    !! 
    !! On success, the shape of the 3D output array will be _output(R,C,D)_.
    !! On error, the 3D output array is __not allocated__.
    CHARACTER(len=*), INTENT(in)                             :: path      !! Path of the input data file 
    REAL(kind=8), INTENT(out), DIMENSION(:,:,:), ALLOCATABLE :: data3d    !! 3D-array with the output values (double precision)
    CHARACTER(len=*), INTENT(in), OPTIONAL                   :: delimiter !! Optional column delimiter(s)
    TYPE(error) :: err                                                    !! Error status of the function
    LOGICAL                                           :: ok
    INTEGER                                           :: i,lc,tlc
    INTEGER                                           :: ndr,ndc,ndd
    INTEGER                                           :: ir,jc,kd,lu
    REAL(kind=8), DIMENSION(:), ALLOCATABLE           :: tmp
    CHARACTER(len=5)                                  :: slc
    CHARACTER(len=15)                                 :: i2s
    CHARACTER(len=:), ALLOCATABLE                     :: line,lm1,zdelim
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: wrds

    zdelim = CHAR(9)//CHAR(32)
    IF (PRESENT(delimiter)) zdelim = delimiter
    ! Check file status
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (.NOT.ok) THEN
      err = error(trim(path)//": no such file",-1) ; RETURN
    ENDIF
    lu = free_lun() 
    IF (lu == -1) THEN ; err = error("Cannot find available logical unit...",-1) ; RETURN ; ENDIF
    ! Open file
    OPEN(lu,FILE=TRIM(path),STATUS='OLD',ACTION='READ')

    ! First pass : 
    ! ------------
    !   - get size (rows, columns, depth) 
    !   - check size consistendcy
    !   - check value type
    lc = 0 ; tlc = 0 
    ndr = -1 ; ndc = -1 ; ndd = 1 
    DO WHILE(readline(lu,line)) 
      lm1 = line
      ! Read the line
      lc = lc + 1 ; WRITE(slc,'(I5)') lc ; slc = ADJUSTL(slc)
      ! skip comment line
      IF (INDEX(TRIM(ADJUSTL(line)),"#") == 1) CYCLE
      ! An empty line: new 2D block 
      IF (LEN_TRIM(line) == 0) THEN
        ndd = ndd + 1
        IF (ndr < 0) THEN
          ndr = tlc
        ELSEIF (tlc /= ndr) THEN
          WRITE(i2s,'(I15)') ndd ; i2s = ADJUSTL(i2s)
          err = error(trim(path)//":Invalid number of lines in block #"//i2s//"(line "//TRIM(slc)//")",-5)
          RETURN
        ENDIF
        tlc = 0
        CYCLE
      ENDIF
      tlc = tlc + 1
      ! Splits line in words
      IF (.NOT.tokenize(line,wrds,zdelim,.true.)) THEN 
        ! cannot tokenize
        err = error(trim(path)//": Cannot parse line "//TRIM(slc),-5)
        RETURN
      ELSEIF (.NOT.from_string(wrds,tmp)) THEN
        ! cannot cast values
        err = error(trim(path)//": Cannot cast values at line "//TRIM(slc),-5)
        RETURN
      ELSEIF (ndc > 0 .AND. ndc /= SIZE(wrds)) THEN 
        ! current number of columns not equal to last one
        err = error(trim(path)//": Invalid number of columns (line "//TRIM(slc)//")",-5)
        RETURN
      ENDIF
      IF (ndc == -1) ndc = SIZE(wrds)
      IF (ALLOCATED(wrds)) DEALLOCATE(wrds)
      IF (ALLOCATED(tmp)) DEALLOCATE(tmp)
    ENDDO

    ! NOTE:
    ! there is a possible bug if data file ends by an empty line:
    ! we will have an extra empty bloc !
    !   current patch: we save the last line of the file and check it:
    !     - if we have empty line, we reduce ndd by one.
    IF (LEN_TRIM(lm1) == 0) ndd = ndd-1

    ! Rewind input data file
    REWIND(lu)
    ! Allocate memory
    ALLOCATE(data3d(ndr,ndc,ndd))
    ir = 0 ; kd = 1 ; 
    DO WHILE(readline(lu,line)) 
      IF (INDEX(TRIM(ADJUSTL(line)),"#") == 1) CYCLE
      ir = ir + 1
      ! empty line update block subscripts 
      IF (LEN_TRIM(line) == 0) THEN
        kd = kd + 1 ; ir = 0 ; CYCLE
      ENDIF
      ok = tokenize(line,wrds,zdelim,.true.)
      DO jc = 1,ndc ; ok = from_string(wrds(jc),data3d(ir,jc,kd)) ; ENDDO
      IF (ALLOCATED(wrds)) DEALLOCATE(wrds)
    ENDDO
    CLOSE(lu)
  END FUNCTION read_data_3d

  FUNCTION read_data_2d(path,data2d,delimiter) RESULT(err)
    !! Read an input ASCII data file and stores its content in a 2D array
    !!
    !! The function reads an ASCII file and saves its values in a real(kind=8) 2D-array.
    !! 
    !! The input file:
    !!
    !! - can contains any number of empty lines and/or comment line (i.e. line where first 
    !!   non-blank character is "#"). All other lines are assumed to be data.
    !! - must have a regular number of columns, that is each data line must have the same 
    !!   number of columns. 
    !! - must use blank space(s) as value delimiter.
    !! 
    !! Error occured when:
    !!
    !! - Path does not refer to a existing file (-1)
    !! - No free logical unit available (-1)
    !! - The file does not have regular data-columns number (-5)
    !! - At least a value cannot be cast in double precision (-5)
    !!
    !! On error, the 2D output array is __not allocated__.
    USE FSYSTEM
    CHARACTER(len=*), INTENT(in)                           :: path      !! Path of the input data file 
    REAL(kind=8), INTENT(out), DIMENSION(:,:), ALLOCATABLE :: data2d    !! 2D-array with the output values (double precision)
    CHARACTER(len=*), INTENT(in), OPTIONAL                 :: delimiter !! Optional column delimiter(s)
    TYPE(error) :: err                                                  !! Error status of function.
    LOGICAL                                           :: ok
    INTEGER                                           :: i,e,vc,lc
    INTEGER                                           :: nl,nc,lu
    REAL(kind=8), DIMENSION(:), ALLOCATABLE           :: tmp
    CHARACTER(len=5)                                  :: slc
    CHARACTER(len=:), ALLOCATABLE                     :: line,zdelim
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: wrds
    zdelim = CHAR(9)//CHAR(32)
    IF (PRESENT(delimiter)) zdelim = delimiter
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (.NOT.ok) THEN
      err = error(trim(path)//": no such file",-1) ; RETURN
    ENDIF
    lu = free_lun() 
    IF (lu == -1) THEN ; err = error("Cannot find available logical unit...",-1) ; RETURN ; ENDIF
    OPEN(lu,FILE=TRIM(path),STATUS='OLD',ACTION='READ')
    vc = 0 ; lc=0 ; ok = .true.
    ! Read the file twice :)
    lc=0 ; vc = 0 ; nc=-1
    ! First pass : get number of row values and checks everything !
    DO 
      ! Read the line
      IF (.NOT.readline(lu,line)) EXIT
      lc = lc + 1 ; WRITE(slc,'(I5)') lc ; slc = ADJUSTL(slc)
      ! skip empty/comment line
      IF (INDEX(TRIM(ADJUSTL(line)),"#") == 1.OR.LEN_TRIM(line) == 0) CYCLE
      ! update row counter
      vc = vc + 1
      IF (.NOT.tokenize(line,wrds,zdelim,.true.)) THEN 
        ! cannot tokenize
        err = error(trim(path)//": Cannot parse line "//TRIM(slc),-5)
        RETURN
      ELSEIF (.NOT.from_string(wrds,tmp)) THEN
        ! cannot cast values
        err = error(trim(path)//": Cannot cast values at line "//TRIM(slc),-5)
        RETURN
      ELSEIF (nc > 0 .AND. nc /= SIZE(tmp)) THEN 
        ! current number of columns not equal to last one
        err = error(trim(path)//": Invalid number of columns (line "//TRIM(slc)//")",-5)
        RETURN
      ENDIF
      IF (nc == -1) nc = SIZE(wrds)
      IF (ALLOCATED(wrds)) DEALLOCATE(wrds)
      IF (ALLOCATED(tmp)) DEALLOCATE(tmp)
    ENDDO
    ! Rewind input data file
    REWIND(lu)
    nl = vc
    ! allocate memory
    ALLOCATE(data2d(nl,nc))
    ! Second pass : saves values :)
    vc = 0 
    DO WHILE(vc <= nl)
      ! Reads the line
      IF (.NOT.readline(lu,line)) EXIT
      ! Check if we have comment or null string
      IF (INDEX(TRIM(ADJUSTL(line)),"#") == 1.OR.LEN_TRIM(line) == 0) CYCLE
      vc = vc + 1
      ok = tokenize(line,wrds,zdelim,.true.)
      DO i = 1,nc
        ok = from_string(wrds(i),data2d(vc,i))
      ENDDO
      IF (ALLOCATED(wrds)) DEALLOCATE(wrds)
    ENDDO
    CLOSE(lu)
    RETURN
  END FUNCTION read_data_2d

  FUNCTION readline(lun,line) RESULT(not_eof)
    !! Read a complete line
    !! 
    !! Each time, it is called, the function reads a complete of the file opened in __lun__ 
    !! logical unit and returns .false. if EOF has been reached, .true. otherwise.
    !!
    !! The function is intended to read a file line by line:
    !!
    !! ```
    !! lu = 1
    !! open(lu,file="path/to/the/file/to/read")
    !! l=0   ! A line counter
    !! DO WHILE(readline(lu,line))
    !!   l = l + 1
    !!   WRITE(*,'("L",I2.2,": ",(a))') il,line
    !! ENDDO
    !! CLOSE(1)
    !! ```
    INTEGER, INTENT(in)                        :: lun  !! Logical unit with the opened file to read. 
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: line !! Processed line 
    LOGICAL  :: not_eof                                !! .true. if EOF has NOT been reached yet, .false. otherwise
    CHARACTER(len=50) :: buf
    INTEGER           :: e,sz
    not_eof = .true. ; line = ''
    DO
      READ(lun,'(a)',ADVANCE="no",SIZE=sz,IOSTAT=e) buf
      IF (e == IOSTAT_END) THEN
        not_eof = .false.
        IF (sz > 0) line=line//buf(1:sz)
        EXIT
      ELSE IF (e == IOSTAT_EOR) THEN
        line = line//buf(1:sz)
        EXIT
      ELSE
        line = line//buf
      ENDIF
    ENDDO
  END FUNCTION readline

END MODULE ASCIIREAD
