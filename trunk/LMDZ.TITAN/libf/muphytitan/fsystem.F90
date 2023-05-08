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

!! file: fsystem.F90
!! summary: File system methods source file.
!! author: J. Burgalat
!! date: 2013-2015,2017


#include "defined.h"

MODULE FSYSTEM
  !! File system methods module
  USE, INTRINSIC :: ISO_C_BINDING
  USE ERRORS
  IMPLICIT NONE

  PUBLIC 

  PRIVATE :: get_umask
  PRIVATE :: c2t

  INTEGER, PARAMETER :: MAX_PATH = 512 !! Maximum length of a path

  TYPE, PUBLIC :: chrono
    !! Define a simple chronometer
    !!
    !! This object can be used to get an approximation of the execution of some piece of code.
    REAL(kind=8), PRIVATE    :: cpu_start = 0d0   
      !! Starting CPU time
    INTEGER(kind=8), PRIVATE :: clock_start = 0d0 
      !! Starting clock time
    LOGICAL, PRIVATE         :: on_run = .false.
      !! Chrono running state.
#if HAVE_FTNPROC
    CONTAINS 
      PROCEDURE :: is_running => chrono_is_running 
      PROCEDURE :: start      => chrono_start
      PROCEDURE :: stop       => chrono_stop
      PROCEDURE :: reset      => chrono_reset
      PROCEDURE :: get        => chrono_get
#endif
  END TYPE chrono

#ifndef FORD_DOC
  ! C interfaces
  INTERFACE
    FUNCTION strlen_c(s) RESULT(length) bind(C,name="strlen")
      !! Get length of C-string up to (but not including) the terminator
      IMPORT C_PTR, C_SIZE_T
      TYPE(C_PTR), INTENT(in), VALUE :: s !! C string (a C_PTR type)
      INTEGER(kind=C_SIZE_T) :: length    !! An integer with the size of the string.
    END FUNCTION strlen_c

    SUBROUTINE free_c(ptr) bind(C,name="free")
      !! Free memory used by a C pointer
      IMPORT C_PTR
      TYPE(C_PTR), INTENT(in), VALUE :: ptr !! TYPE(C_PTR) object with the underlying C pointer to free
    END SUBROUTINE free_c

    FUNCTION errno_c() BIND(C,name="c_get_errno")
      !! Get last error numero
      IMPORT C_INT
      INTEGER(kind=C_INT) :: errno_c !! Last errno
    END FUNCTION errno_c

    FUNCTION usleep_c(usec) BIND(C,name="usleep") 
      !! (attemps to) Sleep for a given number of microseconds
      IMPORT C_INT
      INTEGER(kind=C_INT), INTENT(in), VALUE :: usec !! Number of microseconds to sleep
      INTEGER(kind=C_INT) :: usleep_c !! An integer with 0 on success, last errno otherwise
    END FUNCTION usleep_c

    FUNCTION getgid_c() BIND(C, name="getgid")
      !! Get Group ID
      IMPORT C_INT
      INTEGER(kind=C_INT) :: getgid_c !! Group identifier
    END FUNCTION getgid_c

    FUNCTION getpid_c() BIND(C, name="getpid")
      !! Get Process ID
      IMPORT C_INT
      INTEGER(kind=C_INT) :: getpid_c !! Current process identifier
    END FUNCTION getpid_c

    FUNCTION getuid_c() BIND(C, name="getuid")
      !! Get User ID
      IMPORT C_INT
      INTEGER(kind=C_INT) :: getuid_c !! User identifier
    END FUNCTION getuid_c

    FUNCTION umask_c() BIND(C,name="c_umask")
      !! Get the current umask of the session
      IMPORT C_INT
      INTEGER(kind=C_INT) :: umask_c !! Current umask value in decimal system
    END FUNCTION umask_c

    FUNCTION access_c(path,perm) BIND(C,name="c_access")
      !! Check if path is accessible for current user 
      IMPORT c_char, C_INT
      CHARACTER(len=c_char), INTENT(in)      :: path(*)  !! Path to check
      INTEGER(kind=C_INT), INTENT(in), VALUE :: perm     !! User's permission to check
      INTEGER(kind=C_INT)                    :: access_c !! 0 on success, last errno on failure
    END FUNCTION access_c

    FUNCTION create_c(path,mode,asfile,forced) BIND(C,name="c_create") 
      !! Create a directory or a file in given path
      IMPORT c_char, C_INT
      CHARACTER(len=c_char), INTENT(in)      :: path(*)   !! Path to create
      INTEGER(kind=C_INT), INTENT(in), VALUE :: mode,   & !! Decimal permission of the path
                                                asfile, & !! 0 to create a directory, any other value to create file
                                                forced    !! non-zero value to force the creation of intermediate directories
      INTEGER(kind=C_INT)                    :: create_c  !! 0 on success, last errno otherwise
    END FUNCTION create_c

    FUNCTION uname_c(uid) BIND(C, name="c_uname")
      !! Get the name of the given user id
      IMPORT C_INT, c_ptr
      INTEGER(kind=C_INT), INTENT(in), VALUE :: uid     !! User id
      TYPE(C_PTR)                            :: uname_c !! C_PTR to the underlying char* pointer storing user name
    END FUNCTION uname_c 

    FUNCTION gname_c(gid) BIND(C, name="c_gname")
      !! Get the name of the given group id
      IMPORT C_INT, c_ptr
      INTEGER(kind=C_INT), INTENT(in), VALUE :: gid     !! Group id
      TYPE(C_PTR)                            :: gname_c !! C_PTR to the underlying char* pointer storing group name
    END FUNCTION gname_c 

    FUNCTION dirname_c(path) BIND(C,name="c_dirname") 
      !! Get the directory name of the path
      IMPORT c_char, c_ptr
      CHARACTER(kind=c_char), INTENT(in)  :: path(*)   !! Input path
      TYPE(C_PTR)                         :: dirname_c !! C_PTR to the underlying char* pointer storing dirname
    END FUNCTION dirname_c

    FUNCTION basename_c(path) BIND(C,name="c_basename")
      !! Get the basename of the path
      IMPORT c_char, c_ptr
      CHARACTER(kind=c_char), INTENT(in)  :: path(*)    !! Input path
      TYPE(C_PTR)                         :: basename_c !! C_PTR to the underlying char* pointer sotring basename
    END FUNCTION basename_c

    FUNCTION getcwd_c() BIND(C,name="c_getcwd") 
      !! Get the current working directory
      IMPORT c_ptr
      TYPE(C_PTR) :: getcwd_c !! C_PTR to the underlying char* pointer storing current working directory
    END FUNCTION getcwd_c

    FUNCTION realpath_c(path) BIND(C, name="c_realpath")
      !! Get the real path from given path
      IMPORT c_char, c_ptr
      CHARACTER(kind=c_char), INTENT(in)  :: path(*)    !! The path to expand
      TYPE(C_PTR)                         :: realpath_c !! C_PTR to the underlying char* pointer storing realpath
    END FUNCTION realpath_c

    FUNCTION relpath_c(fname,reldir) BIND(C, name="c_relpath")
      !! Get the relative path of path from another
      IMPORT c_char, c_ptr
      CHARACTER(kind=c_char), INTENT(in) :: fname(*), & !! Path to process
                                            reldir(*)   !! New base path
      TYPE(C_PTR)                        :: relpath_c   !! C_PTR to the underlying char* pointer storing relative path
    END FUNCTION

    FUNCTION rename_c(input,output) BIND(C,name="c_rename")
      !! Rename a path
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in) :: input(*)  !! Path to rename
      CHARACTER(kind=c_char), INTENT(in) :: output(*) !! New name of the path
      INTEGER(kind=C_INT)                :: rename_c  !! 0 on success, last errno on failure 
    END FUNCTION rename_c

    FUNCTION chmod_c(path,mode) BIND(C,name="c_chmod")
      !! Change file/directory permissions
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)     :: path(*) !! Path to modify
      INTEGER(kind=C_INT), INTENT(in), VALUE :: mode    !! New decimal permissions of the path to set
      INTEGER(kind=C_INT)                    :: chmod_c !! 0 on success, last errno on failure 
    END FUNCTION chmod_c

    FUNCTION chdir_c(new) BIND(C,name="c_chdir")
      !! Change current directory
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)  :: new(*)  !! Path of the new working directory
      INTEGER(kind=C_INT)                 :: chdir_c !! 0 on success, last errno on failure 
    END FUNCTION chdir_c

    FUNCTION mkdir_c(dirname,mode) BIND(C,name="c_mkdir")
      !! Create directory
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)     :: dirname(*) !! Path of the directory to create
      INTEGER(kind=C_INT), INTENT(in), VALUE :: mode       !! Decimal permission to set 
      INTEGER(kind=C_INT)                    :: mkdir_c    !! 0 on success, last errno on failure 
    END FUNCTION mkdir_c

    FUNCTION mkdirp_c(dirname,mode) BIND(C,name="c_mkdirp")
      !! Create directory recursively
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)     :: dirname(*) !! Path of the directory to create
      INTEGER(kind=C_INT), INTENT(in), VALUE :: mode       !! Decimal permission to set 
      INTEGER(kind=C_INT)                    :: mkdirp_c   !! 0 on success, last errno on failure 
    END FUNCTION mkdirp_c

    FUNCTION copy_c(to,from) BIND(C,name="c_copy") 
      !! Copy a file.
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)  :: to(*)    !! Destination path.
      CHARACTER(kind=c_char), INTENT(in)  :: from(*)  !! Input file path to copy.
      INTEGER(kind=C_INT)                 :: copy_c !! 0 on success, 1 on failure.
    END FUNCTION copy_c

    FUNCTION remove_c(path) BIND(C,name="c_remove") 
      !! Remove a file (or a directory) from the filesystem
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in)  :: path(*)  !! Path to delete
      INTEGER(kind=C_INT)                 :: remove_c !! 0 on success, last errno on failure 
    END FUNCTION remove_c

    FUNCTION rmdir_c(dirpath) BIND(C,name="c_rmdir")
      !! Remove empty directory
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in) :: dirpath(*) !! Directory to delete
      INTEGER(kind=C_INT)                :: rmdir_c    !! 0 on success, last errno on failure 
    END FUNCTION rmdir_c

    FUNCTION rmdirf_c(dirpath) BIND(C,name="c_rmdir_f")
      !! Remove directory (forced)
      IMPORT c_char, C_INT
      CHARACTER(kind=c_char), INTENT(in) :: dirpath(*) !! Directory to delete
      INTEGER(kind=C_INT)                :: rmdirf_c   !! 0 on success, last errno on failure 
    END FUNCTION rmdirf_c

    FUNCTION fstat_c(pa, pe, ln, ty, ui, gi, si, at, mt, ct) BIND(C, name='c_fstat')
      !! Get file informations
      IMPORT c_char, c_long, C_INT
      CHARACTER(kind=c_char), INTENT(in)  :: pa(*)   !! Path of the file
      INTEGER(kind=C_INT), INTENT(out)    :: pe      !! Decimal permission of the path
      INTEGER(kind=C_INT), INTENT(out)    :: ln      !! Number of links
      INTEGER(kind=C_INT), INTENT(out)    :: ty      !! Type of the path
      INTEGER(kind=C_INT), INTENT(out)    :: ui      !! User ID of the path
      INTEGER(kind=C_INT), INTENT(out)    :: gi      !! Group ID of the path
      INTEGER(kind=c_long), INTENT(out)   :: si      !! Size of the path
      CHARACTER(kind=c_char), INTENT(out) :: at(20)  !! Last access date
      CHARACTER(kind=c_char), INTENT(out) :: mt(20)  !! Last modification date
      CHARACTER(kind=c_char), INTENT(out) :: ct(20)  !! Creation date
      INTEGER(kind=C_INT)                 :: fstat_c !! 0 on success, last errno on failure
    END FUNCTION fstat_c

    FUNCTION termsize_c(r,c) BIND(C, name='c_termsize')
      !! Get terminal window size
      IMPORT C_INT
      INTEGER(kind=C_INT), INTENT(out) :: r, &       !! Number of rows
                                          c          !! Number of columns
      INTEGER(kind=C_INT)              :: termsize_c !! 0 on success, last errno on failure 
    END FUNCTION termsize_c

    FUNCTION getCurrentRSS_c() BIND(C, name="c_getCurrentRSS")
      !! Get the current resident set size memory in bytes.
      IMPORT  C_SIZE_T
      INTEGER(kind=C_SIZE_T) :: getCurrentRSS_c !! Current resident set size in bytes (0 if not available).
    END FUNCTION getCurrentRSS_c 

    FUNCTION getPeakRSS_c() BIND(C, name="c_getPeakRSS")
      !! Get the peak resident set size memory in bytes.
      IMPORT  C_SIZE_T
      INTEGER(kind=C_SIZE_T) :: getPeakRSS_c !! Peak resident set size in bytes (0 if not available).
    END FUNCTION getPeakRSS_c 

    FUNCTION getSystemMemory_c(total,avail,free) BIND(C, name='c_getSystemMemory')
      !! Get global memory informations.
      IMPORT C_LONG_LONG,C_INT
      INTEGER(kind=C_LONG_LONG), INTENT(out) :: total             !! Total available memory.
      INTEGER(kind=C_LONG_LONG), INTENT(out) :: avail             !! Current available memory.
      INTEGER(kind=C_LONG_LONG), INTENT(out) :: free              !! Current free memory.
      INTEGER(kind=C_INT)                    :: getSystemMemory_c !! status, 0 on success, 1 otherwise. 
    END FUNCTION getSystemMemory_c
  END INTERFACE
#endif

  CONTAINS

  FUNCTION fstring(string) RESULT(str)
    !! Convert C string to  Fortran string 
    !!
    !! The method copies the input C string up to the last C_NULL_CHAR found (not including it),
    !! and returns the converted Fortran string.
    !! All other C_NULL_CHAR found in the C string are removed.
    CHARACTER(len=*), INTENT(in) :: string !! A string from C 
    CHARACTER(len=:), ALLOCATABLE :: str   !! Converted fortran string
    INTEGER :: i,idx 
    str = ""
    idx = INDEX(string,C_NULL_CHAR,.true.)
    IF (idx == 0) THEN
      str = string
    ELSE
      DO i=1,idx-1
        IF (string(i:i) /= C_NULL_CHAR) str = str//string(i:i)
      ENDDO
    ENDIF
    str = TRIM(str)
  END FUNCTION fstring

  FUNCTION cstr2fstr(cstr) RESULT(fstr)
    !! Get a Fortran (allocatable) string from a C string
    !!
    !! The method build the fortran string from a TYPE(C_PTR) object that represent a
    !! C char\* pointer string. 
    !! @note
    !! If __cstr__ is not allocated (i.e. the C_PTR is not associated) or if it is set
    !! to a C empty string (i.e. '\0') then the method returns an empty string.
    !! @attention
    !! The method does not free the underlying C string and it should be free using 
    !! the subroutine free_c(_cstr_).
    TYPE(C_PTR), INTENT(in) :: cstr
      !! A TYPE(C_PTR) that represent the pointer to the C char array.
    CHARACTER(len=:), ALLOCATABLE :: fstr
      !! An allocatable Fortran string with the content of the input char array.
    CHARACTER(len=1,kind=C_CHAR), DIMENSION(:), POINTER :: pchars
    INTEGER                                             :: i,length
    IF (.NOT.C_ASSOCIATED(cstr)) THEN
      fstr = ""
      RETURN
    ENDIF
    length = INT(strlen_c(cstr), kind=4)
    IF (length ==0) THEN
      fstr = ""
      RETURN
    ENDIF
    CALL C_F_POINTER(cstr,pchars,(/length/))
    ALLOCATE(CHARACTER(len=length) :: fstr)
    DO i=1,length
      fstr(i:i) = pchars(i)
    ENDDO
  END FUNCTION cstr2fstr


  FUNCTION cstring(string) RESULT(str)
    !> convert Fortran string to cstring 
    !!
    !! The method returns a copy of the input string suitable for C functions argument.
    !! @note 
    !! Input string is trimmed during computations
    CHARACTER(len=*), INTENT(in) :: string
      !! A fortran string
    CHARACTER(len=:,kind=C_CHAR), ALLOCATABLE :: str
      !! Same string as __string__ except that C_NULL_CHAR is appended at the end
    INTEGER :: slen
    slen = LEN_TRIM(string)
    ALLOCATE(CHARACTER(len=slen+1,kind=C_CHAR) :: str)
    str(:slen) = TRIM(string) ; str(slen+1:slen+1) = C_NULL_CHAR
  END FUNCTION cstring

!===============================================================================
! C WRAPPER FUNCTIONS/SUBROUTINES
!===============================================================================

  FUNCTION fs_getgid() RESULT(ret) 
    !! Get Group ID
    INTEGER(kind=4) :: ret !! An integer with the group identifier
    ret = INT(getgid_c(),kind=4) 
    RETURN
  END FUNCTION fs_getgid

  FUNCTION fs_getpid() RESULT(ret)
    !! Get Process ID
    INTEGER(kind=4) :: ret !! An integer with the current process identifier
    ret = INT(getpid_c(),kind=4)
    RETURN
  END FUNCTION fs_getpid

  FUNCTION fs_getuid() RESULT(ret) 
    !! Get User ID
    INTEGER(kind=4) :: ret !! An integer with the user identifier
    ret = INT(getuid_c(),kind=4)
    RETURN
  END FUNCTION fs_getuid

  FUNCTION fs_gname(gid) RESULT(gname)
    !! Get a group name from a group id
    INTEGER, INTENT(in) :: gid             !! User id to check
    CHARACTER(len=:), ALLOCATABLE :: gname !! A string with the name of the group id
    TYPE(C_PTR) :: zname
    zname = gname_c(gid)
    IF (.NOT.C_ASSOCIATED(zname)) THEN
      gname = "" 
    ELSE
      gname = cstr2fstr(zname)
    ENDIF
    CALL free_c(zname)
  END FUNCTION fs_gname

  FUNCTION fs_uname(uid) RESULT(uname)
    !! Get a user name from a user id
    INTEGER, INTENT(in) :: uid             !! User id to check
    CHARACTER(len=:), ALLOCATABLE :: uname !! A string with the name of the user id
    TYPE(C_PTR) :: zname
    zname = gname_c(uid)
    IF (.NOT.C_ASSOCIATED(zname)) THEN
      uname = "" 
    ELSE
      uname = cstr2fstr(zname)
    ENDIF
    CALL free_c(zname)
  END FUNCTION fs_uname

  FUNCTION fs_dirname(path)  RESULT(opath)
    !! Get the parent directory path of the given path
    CHARACTER(len=*), INTENT(in)  :: path
      !! A string with a (valid) path
    CHARACTER(len=:), ALLOCATABLE :: opath 
      !! A Fortran allocated string with the parent directory path or an empty string if method fails
    TYPE(C_PTR) :: zpath
    IF (LEN_TRIM(path) == 0) THEN
      opath = ""
      RETURN
    ENDIF
    zpath = dirname_c(cstring(ADJUSTL(path)))
    IF (.NOT.C_ASSOCIATED(zpath)) THEN
      opath = ""
    ELSE
      opath = cstr2fstr(zpath)
    ENDIF
    CALL free_c(zpath)
  END FUNCTION fs_dirname

  FUNCTION fs_basename(path) RESULT(opath)
    !! Get the base name of the path
    CHARACTER(len=*), INTENT(in)  :: path
      !! A string with a (valid) path
    CHARACTER(len=:), ALLOCATABLE :: opath 
      !! The basename of the path or an empty string if method fails
    TYPE(C_PTR) :: zpath
    IF (LEN_TRIM(path) == 0) THEN
      opath = ""
      RETURN
    ENDIF
    zpath = basename_c(cstring(ADJUSTL(path)))
    IF (.NOT.C_ASSOCIATED(zpath)) THEN
      opath = ""
    ELSE
      opath = cstr2fstr(zpath)
    ENDIF
    CALL free_c(zpath)
  END FUNCTION fs_basename

  FUNCTION fs_realpath(path) RESULT(opath)
    !! Get the real path of the path
    !!
    !! The method computes the absolute path of the given path using C realpath function.
    !! @note 
    !! If the input path is empty then current working directory is returned.
    CHARACTER(len=*), INTENT(in)  :: path
      !! A string with a (valid) path
    CHARACTER(len=:), ALLOCATABLE :: opath 
      !! The absolute of the path or an empty string if method fails
    TYPE(C_PTR) :: zpath
    zpath = realpath_c(cstring(ADJUSTL(path)))
    IF (.NOT.C_ASSOCIATED(zpath)) THEN
      opath = ""
    ELSE
      opath = cstr2fstr(zpath)
    ENDIF
    CALL free_c(zpath)
  END FUNCTION fs_realpath

  FUNCTION fs_relpath(path,reldir) RESULT(res)
    !! Get the relative representation of two paths
    !!
    !! The method computes the relative representation of __path__ from __reldir__ if possible. 
    !! If no common prefix is found, the method returns __path__.
    CHARACTER(len=*), INTENT(in) :: path, & !! Path to be computed relative to reldir
                                    reldir  !! A directory path from which output should be relative to
    CHARACTER(len=:), ALLOCATABLE :: res    !! An allocated string with the resulting path
    TYPE(C_PTR) :: zpath
    zpath = relpath_c(cstring(ADJUSTL(path)),cstring(ADJUSTL(reldir)))
    IF (.NOT.C_ASSOCIATED(zpath)) THEN
      res = TRIM(ADJUSTL(path))
    ELSE
      res = cstr2fstr(zpath)
    ENDIF 
    CALL free_c(zpath)
  END FUNCTION fs_relpath

  FUNCTION fs_getcwd() RESULT(path) 
    !! Get the current working directory
    CHARACTER(len=:), ALLOCATABLE :: path
      !! The current working directory or an empty string if method fails
    TYPE(C_PTR) :: zpath
    zpath = getcwd_c()
    IF (C_ASSOCIATED(zpath)) THEN
      path = cstr2fstr(zpath)
    ELSE
      path = ""
    ENDIF
    CALL free_c(zpath)
    RETURN
  END FUNCTION fs_getcwd

  FUNCTION fs_copy(input,output) RESULT(ret)
    !! Copy input file into output file.
    CHARACTER(len=*), INTENT(in)  :: input  !! Input file path to copy.
    CHARACTER(len=*), INTENT(in)  :: output !! Output file path destination.
    LOGICAL :: ret                          !! True on success, false otherwise.
    IF (LEN_TRIM(input) == 0 .OR. LEN_TRIM(output) == 0 .OR. input == output) THEN
      ret = .false.
    ELSE
      ret = INT(copy_c(cstring(ADJUSTL(output)),cstring(ADJUSTL(input)))) == 0
    ENDIF
    RETURN
  END FUNCTION fs_copy

  FUNCTION fs_remove(path) RESULT(ret)
    !! Delete the file/directory pointed by the given path
    CHARACTER(len=*), INTENT(in)  :: path !! A string with the (valid) file path to delete
    LOGICAL :: ret                        !! True on success, false otherwise.
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false.
    ELSE
      ret = INT(remove_c(cstring(ADJUSTL(path)))) == 0
    ENDIF
    RETURN
  END FUNCTION fs_remove

  FUNCTION fs_rename(old, new) RESULT(ret)
    !! Rename a path
    CHARACTER(len=*), INTENT(in) :: old, & !! A string with the (valid) path to rename
                                    new    !! A string with the new name of the path
    LOGICAL :: ret                         !! True on success, false otherwise.
    IF (LEN_TRIM(old) == 0.OR.LEN_TRIM(new) == 0) THEN
      ret = .false. 
    ELSE
      ret = INT(rename_c(cstring(ADJUSTL(old)),cstring(ADJUSTL(new)))) == 0
    ENDIF
    RETURN
  END FUNCTION fs_rename

  FUNCTION fs_chmod(path, mode) RESULT(ret)
    !! Change file/directory permissions
    CHARACTER(len=*), INTENT(in) :: path !! Path to modify
    INTEGER, INTENT(in)          :: mode !! New octal permissions of the file
    LOGICAL  :: ret                      !! True on success, false otherwise.
    INTEGER(kind=C_INT) :: zmode
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false. 
    ELSE
      zmode = INT(oct_2_dec(mode),kind=C_INT)
      ret = INT(chmod_c(cstring(ADJUSTL(path)), zmode)) == 0
    ENDIF
    RETURN
  END FUNCTION fs_chmod

  FUNCTION fs_chdir(path) RESULT(ret)
    !! Change current working directory
    CHARACTER(len=*), INTENT(in) :: path !! Path of the new working directory
    LOGICAL :: ret                       !! True on success, false otherwise.
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false. 
    ELSE
      ret = INT(chdir_c(cstring(ADJUSTL(path)))) == 0
    ENDIF
    RETURN
  END FUNCTION fs_chdir

  FUNCTION fs_mkdir(path, mode, permissive) RESULT(ret)
    !! Create directory
    !!
    !! The method attempts to create a new directory pointed by __path__ with the permission
    !! given by mode.
    CHARACTER(len=*), INTENT(in)  :: path 
      !! The path to modify
    INTEGER, INTENT(in), OPTIONAL :: mode
      !! Optional octal permission to set for the new directory
    LOGICAL, INTENT(in), OPTIONAL :: permissive
      !! Optional boolean with .true. to create intermediate directories in the path
    LOGICAL :: ret
      !! True on success, false otherwise.
    INTEGER :: zmode
    LOGICAL :: zperm
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false. 
    ELSE
      zmode = oct_2_dec(744) 
      IF (PRESENT(mode)) THEN
        IF (.NOT.chk_pm(mode)) THEN 
          ret = .false. ; RETURN
        ENDIF
        zmode = oct_2_dec(mode)
      ENDIF
      zperm = .false. ; IF (PRESENT(permissive)) zperm = permissive 
      IF (zperm) THEN
        ret = INT(mkdirp_c(cstring(ADJUSTL(path)),INT(zmode,kind=C_INT))) == 0
      ELSE
        ret = INT(mkdir_c(cstring(ADJUSTL(path)),INT(zmode,kind=C_INT))) == 0
      ENDIF
    ENDIF
    RETURN
  END FUNCTION fs_mkdir

  FUNCTION fs_rmdir(path,forced) RESULT(ret)
    !! Remove directory
    !!
    !! By default, the function removes an __empty__ directory. If __forced__ is given and set 
    !! to .true. then the function recursively deletes the directory and __ALL__ its content.
    CHARACTER(len=*), INTENT(in)  :: path
      !! The path of the directory to delete
    LOGICAL, INTENT(in), OPTIONAL :: forced
      !! Optional boolean with @ti{.true.} to remove all contents of the directory.
    LOGICAL :: ret
      !! True on success, false otherwise.
    LOGICAL :: zforce 
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false. 
    ELSE
      zforce = .false. ; IF (PRESENT(forced)) zforce = forced
      IF (.NOT.zforce) THEN
        ret = INT(rmdir_c(cstring(ADJUSTL(path)))) == 0
      ELSE
        ret = INT(rmdirf_c(cstring(ADJUSTL(path)))) == 0
      ENDIF
    ENDIF
    RETURN
  END FUNCTION fs_rmdir

  FUNCTION fs_stat(path,type,perm,nlnks,uid,gid,fsize,atime,mtime,ctime) RESULT(ret)
    !! Get some informations about a path
    !!
    !! The method retrieves various informations about the input path using fstat C function. 
    !! The type of path as returned in __type__ argument is can take the following values:
    !!
    !! - 0, a file
    !! - 1, a link to a file
    !! - 2, a directory
    !! - 3, a link to a directory
    !! - 4, other (fifo, socket, block special, char special...)
    CHARACTER(len=*), INTENT(in)             :: path     !! Input path
    INTEGER, INTENT(out), OPTIONAL           :: type,  & !! Optional type of path (see function documentation).
                                                perm,  & !! Optional permission of the path
                                                nlnks, & !! Optional number of links to the path 
                                                uid,   & !! Optional user id
                                                gid      !! Optional group id
    INTEGER(kind=8), INTENT(out), OPTIONAL   :: fsize    !! Optional file size
    CHARACTER(len=19), INTENT(out), OPTIONAL :: atime, & !! Optional last access time
                                                mtime, & !! Optional last modification time
                                                ctime    !! Optional creation time
    LOGICAL :: ret                                       !! True on success, false otherwise.
    INTEGER                       :: ty,pe,ln,ud,gd 
    INTEGER(kind=8)               :: fs
    CHARACTER(len=:), ALLOCATABLE :: at,mt,ct
    INTEGER(kind=C_INT)           :: p,l,t,u,g
    INTEGER(kind=c_long)          :: f
    CHARACTER(len=20,kind=C_CHAR) :: ta,tm,tc
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false.; RETURN
    ELSE IF (.NOT.(PRESENT(type)  .OR. PRESENT(perm)  .OR. PRESENT(nlnks) .OR. &
                   PRESENT(uid)   .OR. PRESENT(gid)   .OR. PRESENT(fsize) .OR. &
                   PRESENT(atime) .OR. PRESENT(mtime) .OR. PRESENT(ctime))) THEN
      ret = .true.
    ELSE
      ! set default values
      pe=-1 ; ty=-1 ; ud=-1 ; gd=-1 ; fs=-1 ; at="" ; mt="" ; ct=""
      ret = INT(fstat_c(cstring(ADJUSTL(path)),p,l,t,u,g,f,ta,tm,tc)) == 0
      IF (ret) THEN
        pe=INT(p) ; ln=INT(l) ; ty=INT(t) ; ud=INT(u) ; gd=INT(g) 
        fs=INT(f,kind=8) 
        at = fstring(ta)
        mt = fstring(tm)
        ct = fstring(tc)
      ENDIF
      IF (PRESENT(type))  type  = ty 
      IF (PRESENT(perm))  perm  = pe
      IF (PRESENT(nlnks)) nlnks = ln
      IF (PRESENT(uid))   uid   = ud
      IF (PRESENT(gid))   gid   = gd
      IF (PRESENT(fsize)) fsize = fs
      IF (PRESENT(atime)) atime = at
      IF (PRESENT(mtime)) mtime = mt
      IF (PRESENT(ctime)) ctime = ct
    ENDIF
    RETURN
  END FUNCTION fs_stat

  FUNCTION fs_isdir(path) RESULT (ret)
    !! Check if a path is a directory
    !!
    !! The method is just a wrapper of [[fsystem(module):fs_stat(function)]] to get only specific 
    !! information about __path__ type.
    CHARACTER(len=*), INTENT(in) :: path !! The path to check
    LOGICAL :: ret                       !! .true. if the path is a directory, .false. otherwise. 
    INTEGER :: ty
    ret = fs_stat(path,type=ty)
    ret = ret.AND.(ty==2.or.ty==3) 
    RETURN
  END FUNCTION fs_isdir

  FUNCTION fs_isfile(path) RESULT (ret)
    !! Check if a path is a file 
    !!
    !! The method is just a wrapper of [[fsystem(module):fs_stat(function)]] to get only specific 
    !! information about __path__ type.
    CHARACTER(len=*), INTENT(in) :: path !! The path to check
    LOGICAL :: ret                       !! .true. if the path is a file, .false. otherwise. 
    INTEGER :: ty
    ret=fs_stat(path,type=ty)
    ret = ret.and.(ty==0.or.ty==1)
    RETURN
  END FUNCTION fs_isfile

  FUNCTION fs_islink(path) RESULT (ret)
    !! Check if a path is a link 
    !!
    !! The method is just a wrapper of [[fsystem(module):fs_stat(function)]] to get only specific 
    !! information about __path__ type.
    CHARACTER(len=*), INTENT(in) :: path !! The path to check
    LOGICAL :: ret                       !! .true. if the path is a link, .false. otherwise. 
    INTEGER :: ty 
    ret=fs_stat(path,type=ty)
    ret = ret.and.(ty==1.or.ty==3)
    RETURN
  END FUNCTION fs_islink

  FUNCTION fs_access(path,permission) RESULT(ret)
    !! Check if a path is accessible for current user
    !!
    !! The method checks if the given path is accessible for the current user. By default,
    !! it does not check for specific permissions. If __permission__ is given it should be
    !! an integer between 0 and 7 resulting from the possible combinations:
    !!
    !! - 0 : Checks for path existence (default)
    !! - 1 : Checks for EXECUTE permission
    !! - 2 : Checks for WRITE permission
    !! - 4 : Checks for READ permission 
    CHARACTER(len=*), INTENT(in)  :: path       !! Path to check
    INTEGER, INTENT(in), OPTIONAL :: permission !! Optional permission to check
    LOGICAL :: ret                              !! True on success, false otherwise.
    INTEGER(kind=C_INT) :: zp
    IF (LEN_TRIM(path) == 0) THEN
      ret = .false. 
    ELSE
      zp = 0 ; IF (PRESENT(permission)) zp = INT(permission,kind=C_INT)
      ! Defaults are set in the C function.
      ret = INT(access_c(cstring(ADJUSTL(path)),zp)) == 0
    ENDIF
    RETURN
  END FUNCTION fs_access

  FUNCTION fs_split_ext(path, base, ext, absolute) RESULT(ret)
    !! Split given path into base,extension
    !!
    !! The __base__ of a path is conventionnally defined as all characters before the last dot of the path. 
    !! The extension (__ext__) of the path gathers consequently all characters from the last dot to the end 
    !! of the string.
    !! @note
    !! If the basename of the path begins by a dot then the path is assumed to be an hidden file (directory). 
    !! __ext__ will then be empty.
    CHARACTER(len=*), INTENT(in)               :: path     !! Path to split 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: base, &  !! Output base of the path
                                                  ext      !! Output extension of the path
    LOGICAL, INTENT(in), OPTIONAL              :: absolute !! .true. to return absolute path 
    LOGICAL                       :: ret                   !! .true. on success, .false. otherwise. 
    LOGICAL                       :: zabs
    INTEGER                       :: p
    CHARACTER(len=:), ALLOCATABLE :: d,b,apath
    base = "" ; ext = "" 
    ret = .false.
    IF (LEN_TRIM(path) == 0) THEN
      RETURN
    ENDIF
    zabs = .false. ; IF (PRESENT(absolute)) zabs = absolute
    apath = TRIM(path)
    IF (zabs) THEN
      apath = fs_realpath(path) ; IF (LEN_TRIM(apath) == 0) RETURN
    ENDIF 
    d = fs_dirname(apath) ; IF (LEN_TRIM(d) == 0) RETURN
    b = fs_basename(apath) ; IF (LEN_TRIM(b) == 0) RETURN
    p = INDEX(b,".",.true.)
    ! If dot is set as first char of basename : it's an hidden file
    IF (p > 1) THEN
      ext = b(p:) ; base = TRIM(d)//"/"//b(:p-1) 
    ELSE
      base = TRIM(apath) 
    ENDIF
    ret = .true.
    RETURN
  END FUNCTION fs_split_ext

  FUNCTION fs_create(path, mode, type, permissive) RESULT(ret)
    !! Create a directory/file 
    !!
    !! The method creates the file/directory pointed by given __path__.
    !! If __type__ is not given, the method builds the path as :
    !!
    !! -# A file if the basename of the path contains an extension
    !! -# A directory in any other cases.
    !!
    !! Otherwise __type__ should be set to "f" for file or "d" for directory.
    !!
    !! Unless __permissive__ is set to .true., the method will fails if intermediate
    !! directories in the path do not exist.
    CHARACTER(len=*), INTENT(in)           :: path        !! Path to create 
    INTEGER, INTENT(in), OPTIONAL          :: mode        !! Optional octal permisions to set
    CHARACTER(len=1), INTENT(in), OPTIONAL :: type        !! Optional type of path to create
    LOGICAL, INTENT(in), OPTIONAL          :: permissive  !! .true. to create intermediate directories if not existing
    LOGICAL :: ret                                        !! True on success, false otherwise.
    INTEGER                       :: zmd,zt,zp
    CHARACTER(len=:), ALLOCATABLE :: b,e 
    ret = .false.
    ! Checking for existence
    IF (LEN_TRIM(path) == 0) THEN
      RETURN
    ELSE IF (fs_access(path)) THEN
      RETURN 
    ENDIF
    ! Set type of path
    IF (PRESENT(type)) THEN
      IF (.NOT.(type(1:1)=="f".OR.type(1:1)=="d")) THEN
        RETURN 
      ELSE
        zt=0 ; IF (type(1:1)=="f") zt = 1
      ENDIF
    ELSE
      IF(.NOT.fs_split_ext(path,b,e)) RETURN
      zt = 0 ; IF (LEN_TRIM(e) /= 0) zt=1
    ENDIF
    ! set permissions according to type
    IF (zt == 0) THEN
      zmd = oct_2_dec(777)-get_umask() 
    ELSE
      zmd = oct_2_dec(666) -get_umask()
    ENDIF
    ! Check mode argument if present
    IF (PRESENT(mode)) THEN
      IF(.NOT.chk_pm(mode)) THEN
        ! not a valid permission : We raise an error and abort
        RETURN
      ELSE
        zmd = oct_2_dec(mode)
      ENDIF
    ENDIF
    zp = 0 ; IF(PRESENT(permissive)) THEN ; IF(permissive) zp=1 ; ENDIF
    ret = INT(create_c(cstring(ADJUSTL(path)),INT(zmd,kind=C_INT),INT(zt,kind=C_INT),INT(zp,kind=C_INT))) == 0
    RETURN
  END FUNCTION fs_create

  FUNCTION fs_get_parent(path, n) RESULT(opath)
    !! Get the nth parent of the given path
    !! 
    !! The method first resolves the given path using [[fsystem(module):fs_realpath(function)]] 
    !! to get an absolute path.
    !! @note 
    !! If __n__ is greater than the maximum parent level of the path, "/" is returned.
    CHARACTER(len=*), INTENT(in)  :: path
      !! Input path
    INTEGER, INTENT(in), OPTIONAL :: n 
      !! The level of the parent to get
    CHARACTER(len=:), ALLOCATABLE :: opath
      !! The nth parent of the given path, or an empty string if the parent can not be computed 
    CHARACTER(len=:), ALLOCATABLE :: zp
    INTEGER                       :: i,mx,zl,mzl
    opath = "" 
    zl = 1 ; IF (PRESENT(n)) zl = MAX(n,1)
    IF (LEN_TRIM(path) == 0) THEN
      RETURN
    ENDIF
    ! Gets the absolute path
    zp = fs_realpath(TRIM(ADJUSTL(path)))
    IF (LEN_TRIM(zp) == 0) RETURN
    ! removing trailing / (only if it's not the first ^^)
    mx = LEN_TRIM(zp) ; IF (zp(mx:mx)=="/".AND.mx/=1) zp(mx:mx) = ""
    ! compute maximum level
    mzl = 1 ; DO i=1,mx ; IF(zp(i:i) == '/') mzl=mzl+1 ; ENDDO
    i=0
    DO 
      mx = INDEX(zp(1:mx),'/',.true.) ; i=i+1
      IF (mx==0.OR.i>=zl.OR.i>=mzl) EXIT 
      mx = mx - 1
    ENDDO
    IF (mx >= 1) THEN
      opath = zp(1:MAX(1,mx-1))
    ELSE 
      opath = "/" 
    ENDIF
    RETURN
  END FUNCTION fs_get_parent

  SUBROUTINE fs_termsize(row, column)
    !! Get the current terminal window size
    !! @attention
    !! If the program is redirected to a file (and maybe some other device), the C
    !! function can raise an error. In that case, the default values (20,80) are
    !! returned by the C function and thus the subroutine !
    INTEGER, INTENT(out) :: row,   &  !! Number of rows of the window
                            column    !! Number of columns of the window
    INTEGER(kind=C_INT) :: r, c, ret
    ret = termsize_c(r,c)
    row = INT(r) ; column = INT(c)
    RETURN
  END SUBROUTINE fs_termsize

  SUBROUTINE fs_usleep(usec)
    !! Sleep for a given number of microseconds
    !! @note 
    !! Currently if C usleep function failed, the system... does not sleep ! 
    INTEGER, INTENT(in) :: usec !! The number of microseconds to sleep for
    INTEGER(kind=C_INT) :: ret 
    ! usleep expects useconds_t (unsigned int) which is given here as a 4-bytes int
    ret = usleep_c(INT(usec,kind=C_INT))
  END SUBROUTINE fs_usleep

  SUBROUTINE fs_msleep(msec)
    !! Sleep for a given number of milliseconds
    INTEGER, INTENT(in) :: msec !! The number of milliseconds to sleep for
    CALL fs_usleep(msec*1000)
  END SUBROUTINE fs_msleep

  FUNCTION fs_get_memory(peak,units) RESULT(mem)
    !! Get the memory usage of the current process.
    LOGICAL, INTENT(in), OPTIONAL          :: peak  !! True to retrieve the peak RSS memory, otherwise retrieve the current RSS memory. Default to False.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: units !! Output units: either 'B' (Bytes),'KB' (Kilo-),'MB' (Mega-),'GB' (Giga-). Default to 'B'.
    REAL(kind=8)                           :: mem   !! Memory usage.
    LOGICAL          :: zpeak
    CHARACTER(len=2) :: zunits
    zpeak = .false. ; IF (PRESENT(peak)) zpeak = peak
    zunits = 'B '   ; IF (PRESENT(units)) zunits = units
    IF (zunits /= 'B' .AND. zunits /= 'KB' .AND. zunits /= 'MB' .AND. zunits /= 'GB') zunits = 'B '
    IF (zpeak) THEN
      mem = REAL(getPeakRSS_c(),kind=8)
    ELSE
      mem = REAL(getCurrentRSS_c(),kind=8)
    ENDIF
    IF (zunits == 'KB') THEN
      mem = mem / 1024d0
    ELSE IF (zunits == 'MB') THEN
      mem = mem / 1048576d0
    ELSE IF (zunits == 'GB') THEN
      mem = mem / 1073741824d0
    ENDIF
    RETURN
  END FUNCTION fs_get_memory

  FUNCTION fs_get_system_memory(total,available,free,units) RESULT(ret)
    !! Get informations about system memory.
    !!
    !! If no informations is available, output arguments are set to 0 and the method returns false.
    REAL(kind=8), INTENT(out), OPTIONAL    :: total      !! Total available memory.
    REAL(kind=8), INTENT(out), OPTIONAL    :: available  !! Current available memory.
    REAL(kind=8), INTENT(out), OPTIONAL    :: free       !! Current free memory.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: units      !! Output units: either 'B' (Bytes),'KB' (Kilo-),'MB' (Mega-),'GB' (Giga-). Default to 'B'.
    LOGICAL                                :: ret        !! True on success, false otherwise.
    LOGICAL          :: zpeak
    CHARACTER(len=2) :: zunits
    INTEGER(kind=8)  :: ztot,zava,zfre    

    zunits = 'B '   ; IF (PRESENT(units)) zunits = units
    IF (zunits /= 'B' .AND. zunits /= 'KB' .AND. zunits /= 'MB' .AND. zunits /= 'GB') zunits = 'B '
    ret = INT(getSystemMemory_c(ztot,zava,zfre),kind=4) == 0
    ztot = ztot * 1024 ; zava = zava * 1024 ; zfre = zfre * 1024 

    IF (PRESENT(total))     total     = ztot
    IF (PRESENT(available)) available = zava
    IF (PRESENT(free))      free      = zfre
    IF (.NOT.ret) RETURN

    IF (zunits == 'KB') THEN
      IF (PRESENT(total))     total     = ztot / 1024d0
      IF (PRESENT(available)) available = zava / 1024d0
      IF (PRESENT(free))      free      = zfre / 1024d0
    ELSE IF (zunits == 'MB') THEN
      IF (PRESENT(total))     total     = ztot / 1048576d0
      IF (PRESENT(available)) available = zava / 1048576d0
      IF (PRESENT(free))      free      = zfre / 1048576d0
    ELSE IF (zunits == 'GB') THEN
      IF (PRESENT(total))     total     = ztot / 1073741824d0
      IF (PRESENT(available)) available = zava / 1073741824d0
      IF (PRESENT(free))      free      = zfre / 1073741824d0
    ENDIF
    RETURN
  END FUNCTION fs_get_system_memory


!===============================================================================
! MODULE MISCELLANEOUS METHODS
!===============================================================================

  FUNCTION oct_2_dec(octal) RESULT(res)
    !> Octal to decimal conversion
    !! 
    !! The method converts the octal number ranging from 0 to 777 in the decimal system.
    !! @attention
    !! If the __octal__ number is out of range then the method returns 384 (600 in octal).
    INTEGER, INTENT(in) :: octal !! The octal value to convert
    INTEGER :: res               !! The converted decimal value
    INTEGER :: o,d,i
    IF (octal < 0 .OR. octal > 777) THEN
      res = 384 ; RETURN ! --> 600 in octal : rw-------
    ENDIF
    d = 0 ; i = 0 ; o =  octal
    DO WHILE(o/=0)
      d=d+mod(o,10)*8**i ; i=i+1 ; o=o/10
    ENDDO
    res=d
    RETURN 
  END FUNCTION oct_2_dec

  FUNCTION dec_2_oct(decimal) RESULT(res)
    !! Decimal to octal conversion
    !! The method converts the decimal number ranging from 0 to 511 in the octal system.
    !! @attention
    !! If the __decimal__ number is out of range, then it the method returns 600 (384 in decimal).
    INTEGER, INTENT(in) :: decimal !! The decimal value to convert
    INTEGER :: res                 !! The converted octal value
    ! - LOCAL
    INTEGER :: o,d,i,m
    IF (decimal < 0 .OR. decimal > 511) THEN
      res = 600 ;  RETURN ! --> 384 in decimal : rw-------
    ENDIF
    o=0 ; d = decimal ; i=0 ; m=0
    DO WHILE(d/=0)
      d=d/8 ; m=m+1
    ENDDO
    m=m-1 ; d=decimal
    DO i=0,m
      o=o+mod(d,8)*10**i ; d=d/8
    ENDDO
    res = o
    RETURN
  END FUNCTION dec_2_oct

  FUNCTION sp_2_op(str) RESULT(oct)
    !! Get octal number of string representation's permission
    CHARACTER(len=3),INTENT(in) :: str !! The permission to convert
    INTEGER :: oct                     !! Octal value of the string permission on succes, -1 otherwise. 
    oct = -1
    IF (LEN_TRIM(str) /= 3) RETURN
    SELECT CASE(str)
      CASE("---")  ; oct = 0 
      CASE("--x")  ; oct = 1
      CASE("-w-")  ; oct = 2
      CASE("-wx")  ; oct = 3
      CASE("r--")  ; oct = 4
      CASE("r-x")  ; oct = 5
      CASE("rw-")  ; oct = 6
      CASE("rwx")  ; oct = 7
      CASE DEFAULT 
        oct = -1 ; RETURN
    END SELECT 
    RETURN
  END FUNCTION sp_2_op

  FUNCTION op_2_sp(oct) RESULT(str)
    !! Get string representation of the octal number's permission
    INTEGER, INTENT(in) :: oct !! Octal number to convert
    CHARACTER(len=3) :: str    !! String representation of the octal number on succes, 'ukn' otherwise
    SELECT CASE(oct)
      CASE(0) ; str="---"
      CASE(1) ; str="--x"
      CASE(2) ; str="-w-"
      CASE(3) ; str="-wx"
      CASE(4) ; str="r--"
      CASE(5) ; str="r-x"
      CASE(6) ; str="rw-"
      CASE(7) ; str="rwx"
      CASE DEFAULT 
        str='ukn' ; RETURN
    END SELECT 
    RETURN
  END FUNCTION op_2_sp

  FUNCTION str_perm(oct_perm) RESULT(ret)
    !! Get the string representation of the given permission mask
    INTEGER, INTENT(in) :: oct_perm !! The octal representation of the permission 
    CHARACTER(len=9) :: ret      !! String representation of the octal number on succes, 'ukn' otherwise
    INTEGER :: u,g,o
    IF (.NOT.chk_pm(oct_perm)) THEN 
      ret = "ukn" ; RETURN
    ENDIF
    u=int(oct_perm/100) ; g=int((oct_perm-u*100)/10) ; o=int(oct_perm-u*100-g*10)
    ret(1:3) = op_2_sp(u) ; ret(4:6) = op_2_sp(g) ; ret(7:9) = op_2_sp(o)
    RETURN
  END FUNCTION str_perm

  FUNCTION oct_perm(str) RESULT(ret)
    !! Get the string representation of the given permission mask
    CHARACTER(len=9), INTENT(in) :: str !! The string representation of the permission
    INTEGER :: ret                      !! Octal permission on success, -1 otherwise
    ! - LOCAL
    INTEGER :: u,g,o
    u = sp_2_op(str(1:3)) ; g = sp_2_op(str(4:6)) ; o = sp_2_op(str(7:9))
    IF (u==-1.OR.g==-1.OR.o==-1) THEN
      ret = -1 ; RETURN
    ELSE
      ret = u*100 + g*10 + o
    ENDIF
    RETURN
  END FUNCTION oct_perm

  FUNCTION chk_pm(perm) RESULT(valid)
    !! Check if the given permission is valid
    INTEGER, INTENT(in) :: perm !! Octal permission mask
    LOGICAL :: valid            !! .true. if the permission mask is valid, .false. otherwise
    INTEGER :: u,g,o
    u=int(perm/100) ; g=int((perm-u*100)/10) ; o=int(perm-u*100-g*10)
    valid = (u>=0.AND.u<=7).AND.(g>=0.AND.g<=7).AND.(o>=0.AND.o<=7)
    RETURN
  END FUNCTION chk_pm

  FUNCTION get_umask() RESULT(mask)
    !! Get the umask value of the current session
    INTEGER :: mask !! Current umask value in decimal system
    mask = INT(umask_c())
    RETURN
  END FUNCTION get_umask

  FUNCTION sz2str(file_size) RESULT(fstr)
    !! Get a human readable file size
    INTEGER(kind=8), INTENT(in) :: file_size !! File size (assumed to be bytes)
    CHARACTER(len=50) :: fstr                !! Size in a human readable format
    ! - LOCAL
    INTEGER                                   :: cc
    REAL(kind=8)                              :: zfs
    CHARACTER(len=2), DIMENSION(6), PARAMETER :: sn =  &
                       (/'B ','KB','MB','GB','TB','PB'/)
    zfs=DBLE(file_size)
    DO cc=1,size(sn)-1 ; IF (zfs<1024.) EXIT ; zfs=zfs/1024. ; ENDDO
    IF (MOD(zfs,1.0) == 0) THEN
      WRITE(fstr,'(I50)') INT(zfs) ; fstr = TRIM(ADJUSTL(fstr))//sn(cc)
    ELSE
      WRITE(fstr,'(F50.2)') zfs ; fstr = TRIM(ADJUSTL(fstr))//sn(cc)
    ENDIF
    RETURN
  END FUNCTION sz2str

  FUNCTION chrono_is_running(this) RESULT (ret)
    !! Get chrono's state.
    OBJECT(chrono), INTENT(in) :: this !! Chrono object reference.
    LOGICAL :: ret                    !! Running state.
    ret = this%on_run
    RETURN
  END FUNCTION chrono_is_running

  SUBROUTINE chrono_start(this)
    !! Start the chrono. 
    !!
    !! @note
    !! Calling the method multiple times without explicitly stopping the chrono
    !! [[chrono(type):stop(bound)]] does nothing (except for the first called).
    OBJECT(chrono), INTENT(inout) :: this  !! Chrono object reference.
    IF (.NOT.this%on_run) THEN
      CALL CPU_TIME(this%cpu_start)
      CALL SYSTEM_CLOCK(this%clock_start)
    ENDIF
    this%on_run = .true.
  END SUBROUTINE chrono_start 

  SUBROUTINE chrono_stop(this)
    !! Stop the chrono.
    OBJECT(chrono), INTENT(inout) :: this !! Chrono object reference.
    REAL(kind=8)    :: ecpu
    INTEGER(kind=8) :: eclk,nbm,nbr
    this%on_run = .false.
  END SUBROUTINE chrono_stop

  SUBROUTINE chrono_reset(this)
    !! Reset the chrono's internal elapsed times.
    OBJECT(chrono), INTENT(inout) :: this !! Chrono object reference.
    CALL CPU_TIME(this%cpu_start)
    CALL SYSTEM_CLOCK(this%clock_start)
  END SUBROUTINE chrono_reset

  SUBROUTINE chrono_get(this,cpu,clock,units) 
    !! Get elapsed time since last call of start or reset methods.
    !! 
    !! The method computes the time elapsed in two ways :
    !!
    !! - If the [[fsystem(module):chrono(type)]] is not running, the method retruns 0.
    !! - Otherwise, elapsed time since the last call of 
    !!   [[chrono(type):start(bound)]] (or [[chrono(type):reset(bound)]]).
    OBJECT(chrono), INTENT(in)             :: this
      !! Chrono object reference.
    REAL(kind=8), INTENT(out), OPTIONAL    :: cpu
      !! Elapsed cpu time in seconds by default (see units argument).
    REAL(kind=8), INTENT(out), OPTIONAL    :: clock 
      !! Elapsed system clock time in seconds by default (see units argument).
    CHARACTER(len=2), INTENT(in), OPTIONAL :: units
      !! A two characters wide string with the units to convert in. Units should
      !! be one of the following : 'ms', 's' (default), 'm', 'h' or 'd'. 
    CHARACTER(len=2) :: zu
    REAL(kind=8)     :: cu, fact
    INTEGER(kind=8)  :: ck, r, m
    IF (this%on_run) THEN
      IF (PRESENT(cpu)) THEN
        CALL CPU_TIME(cu)
        cpu = (cu - this%cpu_start)
      ENDIF
      IF (PRESENT(clock)) THEN
        CALL SYSTEM_CLOCK(ck,r,m) 
        clock = c2t(ck,this%clock_start,r,m)
      ENDIF
    ELSE
      IF (PRESENT(cpu))   cpu = 0d0
      IF (PRESENT(clock)) clock = 0d0
    ENDIF
    fact = 1d0
    zu = 's' 
    IF (PRESENT(units))  THEN
      zu = units
      SELECT CASE(zu)
        CASE ('d') ; fact = 3600d0*24.
        CASE ('h') ; fact = 3600d0
        CASE ('m') ; fact = 60d0
        CASE ('ms') ; fact = 1d-3
        CASE DEFAULT ; fact = 1d0
      END SELECT
    ENDIF
    IF (PRESENT(cpu)) cpu = cpu / fact 
    IF (PRESENT(clock)) clock = clock / fact
  END SUBROUTINE chrono_get

  FUNCTION c2t(e,i,r,m) RESULT(time)
    !! Get the real-time between two clock counts from system_clock.
    INTEGER(kind=8), INTENT(in) :: e !! Final clock count
    INTEGER(kind=8), INTENT(in) :: i !! Initial clock count 
    INTEGER(kind=8), INTENT(in) :: r !! Clock count rate
    INTEGER(kind=8), INTENT(in) :: m !! Maximum Clock count value
    REAL(kind=8)    :: time          !! Time in seconds
    INTEGER(kind=8) :: nc
    nc = e-i ; IF (e < i) nc = nc+m
    time = REAL(nc,kind=8)/r
    RETURN
  END FUNCTION c2t
END MODULE FSYSTEM

