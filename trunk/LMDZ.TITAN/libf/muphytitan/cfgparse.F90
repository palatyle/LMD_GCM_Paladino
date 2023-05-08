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

!! file: cfgparse.F90
!! summary: Configuration file parser source file.
!! author: J. Burgalat
!! date: 2013-2015,2017

#include "defined.h"

MODULE CFGPARSE
  !! Configuration file parsing module
  !!
  !! This module defines a set of derived types as well as methods to parse configuration files.
  !!
  !! If you only wish to have an overview of cfgparse usage, you'd better go
  !! [here](|url|/page/swift/p02_cfgparse.html).
  !! @todo
  !! Add interpolation from environment and/or parser options.
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE ERRORS
  USE STRING_OP
  USE FSYSTEM
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cfg_clear, cfg_read_config, cfg_write_config, &
            cfg_get_value, cfg_set_value, cfg_count, cfg_check_name, &
            cfg_has_option, cfg_has_section, &
            cfg_option_names,cfg_section_names, &
            cfg_remove_option, cfg_remove_section, &
            cfg_sort_options

  PUBLIC :: noerror,error, error_to_string, aborting

  ! some public definitions from other modules
  ! from strings
  PUBLIC :: to_lower,st_slen, st_llen

  PUBLIC :: OPERATOR(==), OPERATOR(/=), ASSIGNMENT(=)

  TYPE, PUBLIC :: option
    !! Define an option
    CHARACTER(len=st_slen), PRIVATE :: name = ""       !! Name of the option
    CHARACTER(len=st_slen), PRIVATE :: section = ""    !! Associated section name.
    TYPE(words), PRIVATE            :: values          !! Values of the option
  END TYPE option

  TYPE, PUBLIC :: cfgparser
    !! Define a parser of options
    !!
    !! A [[cfgparser(type)]] stores [[option(type)]] objects.
    TYPE(option), DIMENSION(:), ALLOCATABLE  :: options !! list of options.
#if HAVE_FTNPROC
    CONTAINS
    PROCEDURE, PRIVATE :: cp_get_rv_sc
    PROCEDURE, PRIVATE :: cp_get_dv_sc
    PROCEDURE, PRIVATE :: cp_get_iv_sc
    PROCEDURE, PRIVATE :: cp_get_lv_sc
    PROCEDURE, PRIVATE :: cp_get_cv_sc
    PROCEDURE, PRIVATE :: cp_get_sv_sc
    PROCEDURE, PRIVATE :: cp_get_rv_ve
    PROCEDURE, PRIVATE :: cp_get_dv_ve
    PROCEDURE, PRIVATE :: cp_get_iv_ve
    PROCEDURE, PRIVATE :: cp_get_lv_ve
    PROCEDURE, PRIVATE :: cp_get_cv_ve
    PROCEDURE, PRIVATE :: cp_get_sv_ve
    PROCEDURE, PRIVATE :: cp_set_rv_sc
    PROCEDURE, PRIVATE :: cp_set_dv_sc
    PROCEDURE, PRIVATE :: cp_set_iv_sc
    PROCEDURE, PRIVATE :: cp_set_lv_sc
    PROCEDURE, PRIVATE :: cp_set_cv_sc
    PROCEDURE, PRIVATE :: cp_set_sv_sc
    PROCEDURE, PRIVATE :: cp_set_rv_ve
    PROCEDURE, PRIVATE :: cp_set_dv_ve
    PROCEDURE, PRIVATE :: cp_set_iv_ve
    PROCEDURE, PRIVATE :: cp_set_lv_ve
    PROCEDURE, PRIVATE :: cp_set_cv_ve
    PROCEDURE, PRIVATE :: cp_set_sv_ve
    !> Read configuration file
    PROCEDURE, PUBLIC  :: read_config => cfg_read_config
    !> Write configuration file
    PROCEDURE, PUBLIC  :: write_config => cfg_write_config
    !> Get the number of sections in the parser
    PROCEDURE, PUBLIC  :: count => cfg_count
    !> Clean the parser (delete all options, free memory)
    PROCEDURE, PUBLIC  :: clear => cfg_clear
    !> Get the names of the user-defined sections in the parser
    PROCEDURE, PUBLIC  :: section_names => cfg_section_names
    !> Get the list of options names
    PROCEDURE, PUBLIC  :: option_names => cfg_option_names
    !> Check if parser has option by name
    PROCEDURE, PUBLIC  :: has_option => cfg_has_option
    !> Check if parser has section by name
    PROCEDURE, PUBLIC  :: has_section => cfg_has_section
    !> Remove an option from the parser.
    PROCEDURE, PUBLIC :: remove_option  => cfg_remove_option
    !> Remove a section (and all the associated options) from the parser.
    PROCEDURE, PUBLIC :: remove_section => cfg_remove_section
    !> Get value(s) of an option in the parser by name
    !!
    !! ```
    !! FUNCTION cfg_get_value(this,name,output) RESULT(error)
    !! ```
    !!
    !! The method attempts to get the value(s) of an option that matches _name_ in _this_ parser.
    !!
    !! On error, __output__ argument is undefined (that is, left unchanged
    !! for scalar versions, **unallocated** for vector version).
    !!
    !! Errors occur in the following situations:
    !! - The option has no value (-6)
    !! - The option does not exist (-7)
    !! - The option's value cannot be cast in the desired type (-10)
    GENERIC, PUBLIC :: get_value     => cp_get_rv_sc,cp_get_dv_sc,cp_get_iv_sc, &
                                        cp_get_lv_sc,cp_get_cv_sc,cp_get_sv_sc, &
                                        cp_get_rv_ve,cp_get_dv_ve,cp_get_iv_ve, &
                                        cp_get_lv_ve,cp_get_cv_ve,cp_get_sv_ve
    !> Set value(s) of an option in the parser by name
    !!
    !! ```
    !! FUNCTION cfg_set_value(this,name,input,create) RESULT(error)
    !! ```
    !!
    !! The method searches for the option matching the given _name_ in _this_ parser and sets new
    !! values.
    !!
    !! If _create_ is set to .true. (default to .false.) the method creates the option if does not
    !! exist in _this_ parser.
    !! @warning
    !! In such case, if the given is not valid, an assertion is raised !
    !!
    !! On error (i.e. no option matches the given _name_), no values are set.
    GENERIC, PUBLIC :: set_value     => cp_set_rv_sc,cp_set_dv_sc,cp_set_iv_sc, &
                                        cp_set_lv_sc,cp_set_cv_sc,cp_set_sv_sc, &
                                        cp_set_rv_ve,cp_set_dv_ve,cp_set_iv_ve, &
                                        cp_set_lv_ve,cp_set_cv_ve,cp_set_sv_ve
#endif

  END TYPE cfgparser

  !> Get value(s) of an option in the parser by name.
  !!
  !! ```
  !! FUNCTION cfg_get_value(parser,name,output) RESULT(error)
  !! ```
  !!
  !! The method attempts to get the value(s) of an option that matches _name_ in _this_ parser.
  !!
  !! On error, __output__ argument is undefined (that is, left unchanged
  !! for scalar versions, **unallocated** for vector version).
  !!
  !! Errors occur in the following situations:
  !! - The option has no value (-6)
  !! - The option does not exist (-7)
  !! - The option's value cannot be cast in the desired type (-10)
  INTERFACE cfg_get_value
    MODULE PROCEDURE cp_get_rv_sc,cp_get_dv_sc,cp_get_iv_sc, &
                     cp_get_lv_sc,cp_get_cv_sc,cp_get_sv_sc, &
                     cp_get_rv_ve,cp_get_dv_ve,cp_get_iv_ve, &
                     cp_get_lv_ve,cp_get_cv_ve,cp_get_sv_ve
  END INTERFACE

    !> Set value(s) of an option in the parser by name
    !!
    !! ```
    !! FUNCTION set_value(this,name,input,create) RESULT(error)
    !! ```
    !!
    !! The method searches for the option matching the given _name_ in _this_ parser  and sets new
    !! values.
    !!
    !! If _create_ is set to .true. (default to .false.) the method quietly create the option if does not
    !! exist in _this_ parser.
    !! @warning
    !! In such case, if the given __name__ is not valid, an error is raised !
    !!
    !! On error (i.e. no option matches the given _name_), no values are set.
    INTERFACE cfg_set_value
      MODULE PROCEDURE :: cp_set_rv_sc,cp_set_dv_sc,cp_set_iv_sc, &
                          cp_set_lv_sc,cp_set_cv_sc,cp_set_sv_sc, &
                          cp_set_rv_ve,cp_set_dv_ve,cp_set_iv_ve, &
                          cp_set_lv_ve,cp_set_cv_ve,cp_set_sv_ve
    END INTERFACE

    !> Derived type assignment operator
    !!
    !! This interface defines the assignment operator for the containers defined in the module.
    INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE cp_affect_sc, op_affect_sc
    END INTERFACE

  CONTAINS

  SUBROUTINE op_affect_sc(this,other)
    !! Option object assignment operator subroutine
    !!
    !! The subroutine assigns __other__ to __this__.
    TYPE(option), INTENT(inout) :: this  !! An option object to be assigned
    TYPE(option), INTENT(in)    :: other !! An option object to assign
    CALL words_clear(this%values) ! normally not needed
    this%name = other%name
    this%section = other%section
    this%values = other%values
  END SUBROUTINE op_affect_sc

  FUNCTION op_valid(opt) RESULT(ok)
    !! Check whether or not the option is valid (i.e. has name)
    TYPE(option), INTENT(in)      :: opt  !! An option object
    LOGICAL :: ok !! True if the option is valid, false otherwise.
    ok = LEN_TRIM(opt%name) > 0
  END FUNCTION op_valid

  SUBROUTINE op_clear(opt)
    !! Clear and invalid the given option.
    TYPE(option), INTENT(inout)      :: opt  !! An option object
    opt%name = ''
    opt%section = ''
    CALL words_clear(opt%values)
  END SUBROUTINE op_clear

  FUNCTION op_full_name(opt) RESULT(name)
    !! Get the full name of an option.
    !!
    !! @note
    !! If no section is defined in the option (that should not happen), "__default__" is used
    !! as the section part of the full name.
    TYPE(option), INTENT(in)      :: opt  !! An option object
    CHARACTER(len=:), ALLOCATABLE :: name !! The fullname of the option
    IF (LEN_TRIM(opt%section) == 0) THEN
      name = "__default__/"//TRIM(opt%name)
    ELSE
      name = TRIM(opt%section)//"/"//TRIM(opt%name)
    ENDIF
  END FUNCTION op_full_name

  FUNCTION op_split_name(fname,sname,pname) RESULT(err)
    !> Split a full name in section and option names
    !!
    !! The method splits a full name into (section,option) names:
    !!
    !! - Option basename is always set in lower-case.
    !! - If any, section name case is left unmodified.
    !!
    !! A full name simply consists in a section name and an option name separated by a single "/".
    !!
    !! The method never checks the validity of the output names. Consider using [[cfg_check_name(function)]]
    !! to do so.
    !! @note
    !! If _fname_ does not contains any "/", the method sets the special name "\_\_default\_\_" for the output
    !! section name.
    !! @note
    !! On failure, output arguments are set to empty strings.
    !! @warning
    !! If _fname_ ends with a "/", an error (-9, invalid name) is raised: the method always assumes it can
    !! find an option part in the name.
    CHARACTER(len=*), INTENT(in)               :: fname    !! A name to split
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: sname, & !! Section part of the name
                                                  pname    !! Option part of the name
    TYPE(error)                                :: err      !! Error status of the method
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: tfname
    err = noerror ; pname = "" ; sname = ""
    tfname = op_format(fname,sname,pname)
    IF (LEN_TRIM(tfname) == 0) err = error("Invalid option name ("//TRIM(fname)//")",-9)
  END FUNCTION op_split_name

  FUNCTION op_format(name,sname,pname) RESULT(oname)
    !! Format the input name to be consistent with character case requirements.
    !!
    !! Given a **name**, the method tries to split in section/option names.
    !! Then it converts the option part in lower-case.
    !!
    !! If no section part is found (not '/' or set as first character of **name**), the
    !! special section name `__default__` is set.
    !!
    !! If **name** ends with a '/', it is an error and the method returns an empty string.
    CHARACTER(len=*), INTENT(in)  :: name  !! Name to format.
    CHARACTER(len=:), ALLOCATABLE, INTENT(out), OPTIONAL :: sname !! Section part of the name (optional output)
    CHARACTER(len=:), ALLOCATABLE, INTENT(out), OPTIONAL :: pname !! Option part of the name (optional output)
    CHARACTER(len=:), ALLOCATABLE :: oname                        !! Formatted full option name.
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: zsname,zpname
    zpname = "" ; zsname = ""
    ! splits input name in sname, pname
    idx = INDEX(name,'/')
    IF (idx == LEN_TRIM(name)) THEN
      oname = ''
      IF (PRESENT(sname)) sname = ''
      IF (PRESENT(pname)) pname = ''
      RETURN
    ELSE IF (idx <= 1) THEN
      zsname = "__default__" ; zpname = to_lower(TRIM(name))
      IF (idx == 1) zpname=zpname(2:)
    ELSE
      zsname = name(:idx-1)
      zpname = to_lower(name(idx+1:LEN_TRIM(name)))
    ENDIF
    oname = zsname//"/"//zpname
    IF (PRESENT(sname)) sname = zsname
    IF (PRESENT(pname)) pname = zpname 
  END FUNCTION op_format

  FUNCTION op_greater_than(left,right) RESULT(ret)
    !> greater than operator for option derived type.
    !!
    !! the comparison is made on section and option name (in alphabetical order).
    TYPE(option), INTENT(in) :: left  !! LHS option.
    TYPE(option), INTENT(in) :: right !! RHS option.
    LOGICAL :: ret
      !! .true. if LHS is _greater_ than RHS (based on section and option name)
    ret = LGT(op_full_name(left),op_full_name(right))
  END FUNCTION op_greater_than

  FUNCTION op_less_than(left,right) RESULT(ret)
    !> Less than operator for option derived type.
    !!
    !! the comparison is made on section and option name (in alphabetical order).
    TYPE(option), INTENT(in) :: left  !! LHS option.
    TYPE(option), INTENT(in) :: right !! RHS option.
    LOGICAL :: ret
      !! .true. if LHS is _less_ than RHS (based on section and option name)
    ret = LLT(op_full_name(left),op_full_name(right))
  END FUNCTION op_less_than

  FUNCTION op_to_str(opt,num_values) RESULT(str)
    !! Get the string representation of a option object
    !! @note
    !! If the object is not valid an empty string is returned.
    TYPE(option), INTENT(in) :: opt
      !! A option object
    INTEGER, INTENT(in), OPTIONAL :: num_values
      !! Optional integer with the number of values to print per line
    CHARACTER(len=:), ALLOCATABLE :: str
      !! An allocated string with the representation of the option object
    LOGICAL                                           :: ret
    INTEGER                                           :: nv,np,i
    CHARACTER(len=:), ALLOCATABLE                     :: nspcs
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: vec
    nv = 8 ; IF (PRESENT(num_values)) nv = MAX(1,num_values)
    str = ""
    str = TRIM(opt%name)//" = " ; np = LEN(str)
    ALLOCATE(CHARACTER(len=np) :: nspcs) ; nspcs(1:) = " "
    ! stores the error but do not check...
    ret = words_to_vector(opt%values,vec)
    IF (.NOT.ALLOCATED(vec)) RETURN
    DO i=1,SIZE(vec)
      IF (string_is(TRIM(vec(i))) == st_string.AND.TRIM(vec(i))/="NULL") THEN
        str = str//'"'//remove_quotes(TRIM(vec(i)))//'",'
      ELSE
        str = str//TRIM(vec(i))//','
      ENDIF
      IF (MOD(i,nv) == 0) THEN
        str = str//NEW_LINE('A')//nspcs
      ELSE
        str = str//CHAR(32)
      ENDIF
    ENDDO
    str = TRIM(str)
    IF (str(LEN(str):) == ",") str=str(:LEN(str)-1)
  END FUNCTION op_to_str

  !-------------------------------
  ! DERIVED TYPE cfgparser METHODS
  !-------------------------------

  SUBROUTINE cfg_clear(parser)
    !! Clear the cfgparser object ("destructor")
    !!
    !! This subroutine clears the given parser (deallocates memory).
    OBJECT(cfgparser), INTENT(inout) :: parser !! A cfgparser object to clear
    INTEGER :: i
    IF (ALLOCATED(parser%options)) THEN
      DO i = 1, SIZE(parser%options)
        CALL op_clear(parser%options(i))
      ENDDO
      DEALLOCATE(parser%options)
    ENDIF
  END SUBROUTINE cfg_clear


  FUNCTION cfg_check_name(name) RESULT(valid)
    !! Check if a name is valid.
    !!
    !! If **name** contains a '/' it is assumed to be a full option name. In such case
    !! both parts of the name are checked against section/option names requirements (see below).
    !!
    !! Otherwise it is assumed to be the basename of the option.
    !! 
    !! A valid option (base) name is an alphanumeric sequence in lower-case that always begin by
    !! a letter.
    !!
    !! A valid section name is and alphanumeric sequence (in any case) that always begins by
    !! by a letter.
    CHARACTER(len=*), INTENT(in) :: name !! A string with the name to check.
    LOGICAL :: valid                     !! .true. if the name is valid, .false. otherwise
    INTEGER                       :: i
    CHARACTER(len=26), PARAMETER  :: alpha  = "abcdefghijklmnopqrstuvwxyz"
    CHARACTER(len=26), PARAMETER  :: ualpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    CHARACTER(len=12), PARAMETER  :: num    = "0123456789_"
    CHARACTER(len=:), ALLOCATABLE :: pname,sname
    TYPE(error)                   :: err
    valid = .false.
    i = INDEX(TRIM(name),"/")
    IF (i /= 0) THEN
      err = op_split_name(name,sname,pname)
      IF (err /= 0) THEN
        RETURN
      ENDIF
    ELSE
      pname = to_lower(TRIM(name))
      sname = "__default__"
    ENDIF
    ! check option:
    i = INDEX(pname,CHAR(32))
    IF (i /= 0.OR.LEN_TRIM(pname) <= 0) RETURN
    valid = (VERIFY(pname(1:1),alpha) == 0 .AND.VERIFY(TRIM(pname),alpha//num) == 0)
    IF (.NOT.valid) RETURN
    ! check section
    IF (sname == "__default__") THEN
      valid = .true.
      RETURN
    ELSE
      i = INDEX(sname,CHAR(32))
      IF (i /= 0.OR.LEN_TRIM(sname) <= 0) RETURN
      valid = (VERIFY(sname(1:1),ualpha//alpha) == 0 .AND.VERIFY(TRIM(sname),ualpha//alpha//num) == 0)
    ENDIF
  END FUNCTION cfg_check_name

  FUNCTION cfg_count(this,section) RESULT(num)
    !! Get the total number of option in the parser.
    !!
    !! If a section name is given in argument, the method returns the count for the given section only.
    !!
    !! To get the number of top-level options (i.e. that belongs to the default section) the keyword \_\_default\_\_
    !! should be set for the section argument.
    !!
    !! If __section__ is not defined in the parser, the method returns 0.
    !!
    !! @internal
    !! If no options are defined, then it implies that the internal vector of options is
    !! not allocated.
    OBJECT(cfgparser), INTENT(in)          :: this    !! A cfgparser object to search in
    CHARACTER(len=*), INTENT(in), OPTIONAL :: section !! Optional section name to search in.
    INTEGER :: num                                    !! Number of current options registered in the parser.
    INTEGER :: i
    num = 0
    IF(.NOT.ALLOCATED(this%options)) RETURN
    IF (.NOT.PRESENT(section)) THEN
      num = SIZE(this%options)  
    ELSE
      DO i=1, SIZE(this%options)
        IF (this%options(i)%section == section) num = num+1
      ENDDO
    ENDIF
  END FUNCTION cfg_count

  FUNCTION cfg_section_names(this) RESULT(list)
    !! Get the list of user-defined section names
    !! @note
    !! If the parser does not have user-defined sections, the vector is still
    !! allocated but with 0 elements.
    OBJECT(cfgparser), INTENT(in)                     :: this !! A cfgparser object to process.
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: list !! List of section names.
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER :: i,j,k,no
    LOGICAL :: found
    no = cfg_count(this)
    IF (no == 0) THEN
      ALLOCATE(list(0))
      RETURN
    ENDIF
    ALLOCATE(tmp(no))
    tmp(1) = ''
    ! get the first non-default section
    DO i=1,no
      IF (this%options(i)%section /= "__default__") THEN
        tmp(1) = this%options(i)%section
        EXIT
      ENDIF
    ENDDO
    ! no user defined section
    IF (LEN_TRIM(tmp(1)) == 0) THEN
      DEALLOCATE(tmp)
      ALLOCATE(list(0))
      RETURN
    ENDIF
    k = 1
    DO i=1,no
      found = .false.
      DO j=1,k
        ! Found a match so start looking again
        found = (tmp(j) == this%options(i)%section .OR. &
                 this%options(i)%section == "__default__")
        IF (found) EXIT
      ENDDO
      IF (.NOT.found) THEN
        k = k + 1
        tmp(k) = this%options(i)%section
      ENDIF
    ENDDO
    ALLOCATE(list(k))
    list(1:k) = tmp(1:k)
    DEALLOCATE(tmp)
  END FUNCTION cfg_section_names

  FUNCTION cfg_option_names(this,secname) RESULT(list)
    !! Get the list of option names.
    !!
    !! @note
    !! If the parser does not have options, the vector is still allocated but with 0 elements.
    OBJECT(cfgparser), INTENT(in)                     :: this    !! A cfgparser object to process.
    CHARACTER(len=*), INTENT(in), OPTIONAL            :: secname !! Optional section name to search in.
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: list    !! List of option names.
    INTEGER               :: j,i,no,nso
    no = cfg_count(this)
    IF (no == 0) THEN
       ALLOCATE(list(no)) ; RETURN
    ENDIF
    IF (PRESENT(secname)) THEN
      IF (.NOT.cfg_has_section(this,TRIM(secname))) THEN
        ALLOCATE(list(no)) ; RETURN
      ELSE
        nso = 0
        DO i=1,no ; IF (this%options(i)%section == TRIM(secname)) nso = nso + 1 ; ENDDO
        ALLOCATE(list(nso))
        IF (nso == 0) RETURN
        j = 1
        DO i=1,no
          IF (this%options(i)%section == TRIM(secname)) THEN
            list(j) = TRIM(this%options(i)%section)//"/"//TRIM(this%options(i)%name) ; j=j+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      ALLOCATE(list(no))
      DO i=1,no
        IF (this%options(i)%section == "__default__") THEN
          list(i) = TRIM(this%options(i)%name)
        ELSE
          list(i) = TRIM(this%options(i)%section)//"/"//TRIM(this%options(i)%name)
        ENDIF
      ENDDO
    ENDIF
  END FUNCTION cfg_option_names

  FUNCTION cfg_has_section(this,name) RESULT(yes)
    !! Check if parser has section by name
    !!
    !! @note
    !! Keep in mind that section name in the configuration are case-sensitive.
    OBJECT(cfgparser), INTENT(in) :: this !! cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name !! Name of the section to search
    LOGICAL :: yes                        !! .true. if the section exists, .false. otherwise
    INTEGER                       :: i,no
    yes = .false.
    no = cfg_count(this)
    IF (no == 0) RETURN
    DO i = 1,no
      IF (this%options(i)%section == name) THEN
        yes = .true.
        RETURN
      ENDIF
    ENDDO
  END FUNCTION cfg_has_section

  FUNCTION cfg_has_option(this,name) RESULT(yes)
    !! Check if parser has option by name
    OBJECT(cfgparser), INTENT(in) :: this !! A cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name !! (extended) Name of the option to search
    LOGICAL :: yes                        !! .true. if the option is found, .false. otherwise
    CHARACTER(len=:), ALLOCATABLE :: pname,zname
    INTEGER                       :: i,no,iscan
    yes = .false.
    no = cfg_count(this)
    IF (no == 0) RETURN
    zname = op_format(name)
    IF (LEN_TRIM(zname) == 0) RETURN
    DO i = 1,no
      pname = op_full_name(this%options(i))
      IF (pname == zname) THEN
        yes = .true.
        RETURN
      ENDIF
    ENDDO
  END FUNCTION cfg_has_option

  SUBROUTINE cfg_sort_options(this)
    !! Sort the options in the parser (alphabetiCALLy).
    OBJECT(cfgparser), INTENT(inout) :: this !! A cfgparser object
    INTEGER :: no
    no = cfg_count(this)
    IF (no == 0) RETURN
    CALL insertionSort(this%options)
  END SUBROUTINE cfg_sort_options

  SUBROUTINE cfg_remove_option(this,name)
    !! Remove an option from parser by name.
    OBJECT(cfgparser), INTENT(inout) :: this !! A cfgparser object to search in
    CHARACTER(len=*), INTENT(in)     :: name !! The name of the option to remove
    INTEGER                                 :: no,idx,i,j
    TYPE(option), DIMENSION(:), ALLOCATABLE :: tmp
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) RETURN
    no = cfg_count(this)
    ! only one opt
    IF (no == 1) THEN
      CALL op_clear(this%options(1))
      DEALLOCATE(this%options)
      RETURN
    ENDIF
    ALLOCATE(tmp(no-1))
    j = 1
    DO i=1,no
      IF (i /= idx) THEN
        tmp(j) = this%options(i)
        j= j+1
      ENDIF
      CALL op_clear(this%options(i))
    ENDDO
    DEALLOCATE(this%options)
    ALLOCATE(this%options(no-1))
    DO i=1,no-1
      this%options(i) = tmp(i)
      CALL op_clear(tmp(i))
    ENDDO
    DEALLOCATE(tmp)
  END SUBROUTINE cfg_remove_option

  SUBROUTINE cfg_remove_section(this,name)
    !! Remove a section from parser by name.
    !!
    !! The method removes all the options that belong to the given section name.
    OBJECT(cfgparser), INTENT(inout) :: this
      !! A cfgparser object to search in
    CHARACTER(len=*), INTENT(in)     :: name
      !! The name of the section to remove
    INTEGER                                 :: no,i,j,icount
    INTEGER, DIMENSION(:), ALLOCATABLE      :: idxs,itmp
    TYPE(option), DIMENSION(:), ALLOCATABLE :: tmp
    no = cfg_count(this)
    IF (no == 0) RETURN
    ALLOCATE(itmp(no))
    itmp(:) = -1
    icount = 0
    DO i=1,no
      IF (TRIM(this%options(i)%section) == TRIM(name)) THEN
        itmp(icount+1) = i
        icount = icount + 1
      ENDIF
    ENDDO
    IF (icount == 0) RETURN
    ALLOCATE(idxs(icount))
    idxs(:) = itmp(1:icount)
    !DEALLOCATE(itmp)
    ! all options matches (should be rarely the case): remove all
    IF (icount == no) THEN
       DO i=1,no ; CALL op_clear(this%options(i)); ENDDO
       DEALLOCATE(this%options)
       RETURN
    ENDIF
    ALLOCATE(tmp(icount))
    j = 1
    DO i=1,no
      IF (ANY(idxs == i)) THEN
        tmp(j) = this%options(i)
        j = j + 1
      ENDIF
      CALL op_clear(this%options(i))
    ENDDO
    DEALLOCATE(idxs)
    DEALLOCATE(this%options)
    ALLOCATE(this%options(icount))
    DO i=1,icount
      this%options(i) = tmp(i)
      CALL op_clear(tmp(i))
    ENDDO
    DEALLOCATE(tmp)
  END SUBROUTINE cfg_remove_section

  FUNCTION cfg_read_config(parser,path,override) RESULT(err)
    !! Read configuration file
    !!
    !! @note
    !! If the library support C bindings, the method can read included files which are defined
    !! by the __#include <...>__ directive (see [p_cfgparse](here) from more details).
    OBJECT(cfgparser), INTENT(inout) :: parser
      !! A cfgparser object that will store the configuration
    CHARACTER(len=*), INTENT(in)     :: path
      !! Path of the configuration file
    LOGICAL, INTENT(in), OPTIONAL    :: override
      !! An optional boolean flag with .true. to override previous value of duplicated options instead of raising an error.
    TYPE(error) :: err
      !! An error with the first error encountered
    INTEGER                       :: i
    LOGICAL                       :: zoverride,ok
    TYPE(words)                   :: incfiles
    CHARACTER(len=:), ALLOCATABLE :: name
    CHARACTER(len=st_slen)        :: isec
    err = noerror
    zoverride = .false. ; IF (PRESENT(override)) zoverride = override
    isec = "__default__"
    name = TRIM(path)
    i = INDEX(name,'/',.true.) ; IF (i /= 0) name = name(i+1:)
    IF (i == 0) THEN
      name = fs_realpath("./"//path)
    ELSE
      name = fs_realpath(path)
    ENDIF
    CALL words_append(incfiles,name)
    INQUIRE(FILE=TRIM(path), EXIST=ok)
    IF (.NOT.ok) THEN
      err = error(TRIM(path)//": no such file",-11)
    ELSE
      err = read_include(parser,TRIM(path),isec,incfiles,zoverride)
    ENDIF
    call words_clear(incfiles)
  END FUNCTION cfg_read_config

  FUNCTION cfg_write_config(this,lu,num_values) RESULT(err)
    !> Write the content of the parser in the logical unit.
    !!
    !! the method expects the logical unit to be opened and does not close it
    !! at the end of the process.
    OBJECT(cfgparser), INTENT(inout) :: this
      !! Parser to write.
    INTEGER, INTENT(in) :: lu
      !! Logical unit. It should be already opened or an error is raised.
    INTEGER, INTENT(in), OPTIONAL :: num_values
      !! Optional integer with the number of values to print per line for options (default to 8).
    TYPE(error) :: err
      !! Error status
    CHARACTER(len=st_slen) :: sname
    LOGICAL                :: ok
    INTEGER                :: nv,no,i
    err = noerror
    INQUIRE(UNIT=lu,OPENED=ok)
    IF (.NOT.ok) THEN
      err = error("Logical unit not opened",-15)
      RETURN
    ENDIF
    no = cfg_count(this)
    IF (no == 0) THEN
      err = error("No options to write",-7)
      RETURN
    ENDIF
    ! sort options.
    CALL cfg_sort_options(this)
    nv = 8 ; IF (PRESENT(num_values)) nv = MAX(1,num_values)
    sname = this%options(1)%section
    IF (sname /= "__default__") &
      WRITE(lu,'(a)') "[ "//TRIM(sname)//" ]"
    DO i=1,no
      IF (this%options(i)%section /= sname) THEN
        sname = this%options(i)%section
        ! write section first
        WRITE(lu,*)
        WRITE(lu,'(a)') "[ "//TRIM(sname)//" ]"
      ENDIF
      WRITE(lu,'(a)') op_to_str(this%options(i),nv)
    ENDDO
  END FUNCTION cfg_write_config

  ! internal (private methods)

  SUBROUTINE cp_affect_sc(this,other)
    !! cfgparser object assignment operator subroutine
    !!
    !! The subroutine assigns __other__ to __this__.
    TYPE(cfgparser), INTENT(inout) :: this  !! A cfgparser object to be assigned
    TYPE(cfgparser), INTENT(in)    :: other !! A cfgparser object to assign
    INTEGER :: i,ono
    CALL cfg_clear(this)
    ono = cfg_count(other)
    IF (ono == 0) RETURN
    ALLOCATE(this%options(ono))
    DO i=1,ono
      this%options(i) = other%options(i)
    ENDDO
    RETURN
  END SUBROUTINE cp_affect_sc

  FUNCTION cp_get_opt_idx(this,name) RESULT(idx)
    !! Get the index of an option by name in the parser.
    !!
    !! The method searches in the parser for the option with the given (full) __name__.
    !! If found, it returns the index of the option in the internal vector of options. Otherwise
    !! -1 is returned.
    OBJECT(cfgparser), INTENT(in) :: this !! A cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name !! A string with the name of the option
    INTEGER :: idx                        !! Index of the option (-1 if not found).
    CHARACTER(len=:), ALLOCATABLE :: zname,pname
    INTEGER                       :: no,i
    idx = -1
    no = cfg_count(this)
    IF (no == 0) RETURN
    zname = op_format(name) ! prepare the name to search.
    IF (LEN_TRIM(zname) == 0) RETURN
    DO i=1,no
      pname = op_full_name(this%options(i))
      IF (pname == zname) THEN
        idx = i
        RETURN
      ENDIF
    ENDDO
  END FUNCTION cp_get_opt_idx

  FUNCTION cp_update_opt(this,sname,pname,values) RESULT(err)
    !! Update an option in the parser.
    !!
    !! The method attempts to update the option in the parser.
    !!
    !! If __sname__ is set to empty string, the method searches for the option
    !! in the default section.
    !!
    !! If no option is found, The the option is appended in the parser. Otherwise it is updated
    !! with the content of __values__.
    !!
    !! If the option name is not valid, the method does nothing and -9 error status is returned.
    !!
    !! @internal
    !! The method performs the same kind of operations than the setters except that it
    !! expects raw data ([[string_op(module):words(type)]]).
    OBJECT(cfgparser), INTENT(inout) :: this   !! cfgparser object to process.
    CHARACTER(len=*), INTENT(in)     :: sname  !! Name of the section.
    CHARACTER(len=*), INTENT(in)     :: pname  !! Basename of the option.
    TYPE(words), INTENT(in)          :: values !! Raw values.
    TYPE(error)                      :: err    !! Error status.
    CHARACTER(len=:), ALLOCATABLE :: zsname,fname
    INTEGER                       :: i
    err = noerror
    zsname = TRIM(sname)
    IF (LEN_TRIM(sname) == 0) zsname = "__default__"
    fname = zsname//"/"//to_lower(TRIM(pname))
    IF (.NOT.cfg_check_name(fname)) THEN
       err = error("Invalid option (no name)",-9)
       RETURN
    ENDIF
    i = cp_get_opt_idx(this,fname)
    IF (i /= -1) THEN
      CALL words_clear(this%options(i)%values)
      this%options(i)%values = values
    ELSE
      err = cp_add_opt(this,zsname,pname,values)
    ENDIF
  END FUNCTION cp_update_opt

  FUNCTION cp_add_opt(this,sname,pname,values) RESULT(err)
    !! Add an option to the parser.
    !!
    !! In order to add an option to the default section, _sname_ should be left empty or set to "\_\_default\_\_".
    !!
    !! The following error code can be returned:
    !!  - 0, no error.
    !!  - -8, the option already exists.
    !!  - -9, option name is not valid.
    OBJECT(cfgparser), INTENT(inout) :: this
      !! A cfgparser object to process.
    CHARACTER(len=*), INTENT(in)     :: sname
      !! Section name.
    CHARACTER(len=*), INTENT(in)     :: pname
      !! Option basename.
    TYPE(words), INTENT(in)          :: values
      !! Values to set.
    TYPE(error) :: err
      !! Return error status.
    CHARACTER(len=:), ALLOCATABLE           :: zsname,fname
    TYPE(option), DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER                                 :: no,i
    TYPE(option)                            :: sca

    err = noerror
    zsname = sname
    no = cfg_count(this)
    IF (LEN_TRIM(zsname) == 0) zsname = "__default__"
    fname = TRIM(zsname)//"/"//to_lower(TRIM(pname))
    ! check name
    IF (.NOT.cfg_check_name(fname)) THEN
      err = error("Invalid option name '"//fname//"'",-9)
      RETURN
    ENDIF
    ! check if opt exists in the parser
    IF (no > 0) THEN
      IF (cp_get_opt_idx(this,fname) /= -1) THEN
        err = error("Duplicate option '"//TRIM(pname)//"' in "//TRIM(zsname),-8)
        RETURN
      ENDIF
    ENDIF

    ! build option
    CALL op_clear(sca)
    sca%name = to_lower(TRIM(pname))
    sca%section = zsname
    sca%values = values

    IF (no == 0) THEN
      ! no options yet -> allocate
      ALLOCATE(this%options(1))
    ELSE
      ! parser has options: increase this%options size (ugly copy).
      ALLOCATE(tmp(no))
      DO i =1,no
        tmp(i) = this%options(i)
        CALL op_clear(this%options(i))
      ENDDO
      DEALLOCATE(this%options)
      ALLOCATE(this%options(no+1))
      DO i =1,no
        this%options(i) = tmp(i)
        CALL op_clear(tmp(i))
      ENDDO
      DEALLOCATE(tmp)
    ENDIF
    ! always add the option at the end.
    this%options(no+1) = sca
    CALL op_clear(sca)
  END FUNCTION cp_add_opt

  FUNCTION cp_get_rv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (real/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    !!  - -10, the value cannot be converted in the output type.
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    REAL(kind=4), INTENT(out)     :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      tmp = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      IF (LEN(tmp) == 0) THEN
        err = error("Option "//TRIM(name)//" has no value",-6)
      ELSE
        IF(.NOT.from_string(tmp,output)) &
        err = error(TRIM(name)//": Cannot convert "//tmp//" to real.",-10)
      ENDIF
    ENDIF
  END FUNCTION cp_get_rv_sc

  FUNCTION cp_get_dv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (double/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    !!  - -10, the value cannot be converted in the output type.
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    REAL(kind=8), INTENT(out)     :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      tmp = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      IF (LEN(tmp) == 0) THEN
        err = error("Option "//TRIM(name)//" has no value",-6)
      ELSE
        IF(.NOT.from_string(tmp,output)) &
        err = error(TRIM(name)//": Cannot convert "//tmp//" to double.",-10)
      ENDIF
    ENDIF
  END FUNCTION cp_get_dv_sc

  FUNCTION cp_get_iv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (integer/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    !!  - -10, the value cannot be converted in the output type.
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    INTEGER, INTENT(out)          :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      tmp = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      IF (LEN(tmp) == 0) THEN
        err = error("Option "//TRIM(name)//" has no value",-6)
      ELSE
        IF(.NOT.from_string(tmp,output)) &
        err = error(TRIM(name)//": Cannot convert "//tmp//" to integer.",-10)
      ENDIF
    ENDIF
  END FUNCTION cp_get_iv_sc

  FUNCTION cp_get_lv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (logical/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    !!  - -10, the value cannot be converted in the output type.
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    LOGICAL, INTENT(out)          :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      tmp = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      IF (LEN(tmp) == 0) THEN
        err = error("Option "//TRIM(name)//" has no value",-6)
      ELSE
        IF(.NOT.from_string(tmp,output)) &
        err = error(TRIM(name)//": Cannot convert "//tmp//" to logical.",-10)
      ENDIF
    ENDIF
  END FUNCTION cp_get_lv_sc

  FUNCTION cp_get_cv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (complex/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    !!  - -10, the value cannot be converted in the output type.
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    COMPLEX, INTENT(out)          :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      tmp = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      IF (LEN(tmp) == 0) THEN
        err = error("Option "//TRIM(name)//" has no value",-6)
      ELSE
        IF(.NOT.from_string(tmp,output)) &
        err = error(TRIM(name)//": Cannot convert "//tmp//" to complex.",-10)
      ENDIF
    ENDIF
  END FUNCTION cp_get_cv_sc

  FUNCTION cp_get_sv_sc(this,name,output) RESULT(err)
    !! Get the first value of an option in the parser by name (string/scalar)
    !!
    !! The following error status can be returned by the method:
    !!  - -7, no option matches the given name.
    !!  - -6, the option does not have value(s).
    OBJECT(cfgparser), INTENT(in) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)  :: name   !! (Full) Name of the option to get
    CHARACTER(len=*), INTENT(out) :: output !! Output value
    TYPE(error) :: err
      !! Error status
    INTEGER :: idx
    !CHARACTER(len=:), ALLOCATABLE :: tmp
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values)== 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      output = TRIM(ADJUSTL(words_get(this%options(idx)%values,1)))
      err = noerror
    ENDIF
  END FUNCTION cp_get_sv_sc

  FUNCTION cp_get_rv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (real/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                        :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                         :: name   !! (Full) Name of the option to get
    REAL(kind=4), INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err
      !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    LOGICAL                                           :: ok
    INTEGER                                           :: i,idx
    CHARACTER(len=15)                                 :: i2s
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ALLOCATE(output(words_length(this%options(idx)%values)))
      ok = words_to_vector(this%options(idx)%values,tmp)
      DO i=1, SIZE(tmp)
        WRITE(i2s,'(I15)') i ; i2s=ADJUSTL(i2s)
        IF (LEN_TRIM(ADJUSTL(tmp(i))) == 0) THEN
          err = error("Cannot get value #"//TRIM(i2s)//" from option "//TRIM(name),-6)
          DEALLOCATE(output) ; EXIT
        ELSE IF (.NOT.from_string(tmp(i), output(i))) THEN
          err = error("Cannot convert value #"//TRIM(i2s)//" from option "//TRIM(name),-10)
          DEALLOCATE(output) ; EXIT
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE(tmp)
    RETURN
  END FUNCTION cp_get_rv_ve

  FUNCTION cp_get_dv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (double/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                        :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                         :: name   !! (Full) Name of the option to get
    REAL(kind=8), INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err
      !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    LOGICAL                                           :: ok
    INTEGER                                           :: i,idx
    CHARACTER(len=15)                                 :: i2s
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ALLOCATE(output(words_length(this%options(idx)%values)))
      ok = words_to_vector(this%options(idx)%values,tmp)
      DO i=1, SIZE(tmp)
        WRITE(i2s,'(I15)') i ; i2s=ADJUSTL(i2s)
        IF (LEN_TRIM(ADJUSTL(tmp(i))) == 0) THEN
          err = error("Cannot get value #"//TRIM(i2s)//" from option "//TRIM(name),-6)
          DEALLOCATE(output) ; EXIT
        ELSE IF (.NOT.from_string(tmp(i), output(i))) THEN
          err = error("Cannot convert value #"//TRIM(i2s)//" from option "//TRIM(name),-10)
          DEALLOCATE(output) ; EXIT
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE(tmp)
    RETURN
  END FUNCTION cp_get_dv_ve

  FUNCTION cp_get_iv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (integer/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                   :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                    :: name   !! (Full) Name of the option to get
    INTEGER, INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err
      !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    LOGICAL                                           :: ok
    INTEGER                                           :: i,idx
    CHARACTER(len=15)                                 :: i2s
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ALLOCATE(output(words_length(this%options(idx)%values)))
      ok = words_to_vector(this%options(idx)%values,tmp)
      DO i=1, SIZE(tmp)
        WRITE(i2s,'(I15)') i ; i2s=ADJUSTL(i2s)
        IF (LEN_TRIM(ADJUSTL(tmp(i))) == 0) THEN
          err = error("Cannot get value #"//TRIM(i2s)//" from option "//TRIM(name),-6)
          DEALLOCATE(output) ; EXIT
        ELSE IF (.NOT.from_string(tmp(i), output(i))) THEN
          err = error("Cannot convert value #"//TRIM(i2s)//" from option "//TRIM(name),-10)
          DEALLOCATE(output) ; EXIT
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE(tmp)
    RETURN
  END FUNCTION cp_get_iv_ve

  FUNCTION cp_get_lv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (logical/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                   :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                    :: name   !! (Full) Name of the option to get
    LOGICAL, INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err
      !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    LOGICAL                                           :: ok
    INTEGER                                           :: i,idx
    CHARACTER(len=15)                                 :: i2s
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ALLOCATE(output(words_length(this%options(idx)%values)))
      ok = words_to_vector(this%options(idx)%values,tmp)
      DO i=1, SIZE(tmp)
        WRITE(i2s,'(I15)') i ; i2s=ADJUSTL(i2s)
        IF (LEN_TRIM(ADJUSTL(tmp(i))) == 0) THEN
          err = error("Cannot get value #"//TRIM(i2s)//" from option "//TRIM(name),-6)
          DEALLOCATE(output) ; EXIT
        ELSE IF (.NOT.from_string(tmp(i), output(i))) THEN
          err = error("Cannot convert value #"//TRIM(i2s)//" from option "//TRIM(name),-10)
          DEALLOCATE(output) ; EXIT
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE(tmp)
    RETURN
  END FUNCTION cp_get_lv_ve

  FUNCTION cp_get_cv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (complex/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                   :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                    :: name   !! (Full) Name of the option to get
    COMPLEX, INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err
      !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tmp
    LOGICAL                                           :: ok
    INTEGER                                           :: i,idx
    CHARACTER(len=15)                                 :: i2s
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ALLOCATE(output(words_length(this%options(idx)%values)))
      ok = words_to_vector(this%options(idx)%values,tmp)
      DO i=1, SIZE(tmp)
        WRITE(i2s,'(I15)') i ; i2s=ADJUSTL(i2s)
        IF (LEN_TRIM(ADJUSTL(tmp(i))) == 0) THEN
          err = error("Cannot get value #"//TRIM(i2s)//" from option "//TRIM(name),-6)
          DEALLOCATE(output) ; EXIT
        ELSE IF (.NOT.from_string(tmp(i), output(i))) THEN
          err = error("Cannot convert value #"//TRIM(i2s)//" from option "//TRIM(name),-10)
          DEALLOCATE(output) ; EXIT
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE(tmp)
    RETURN
  END FUNCTION cp_get_cv_ve

  FUNCTION cp_get_sv_ve(this,name,output) RESULT(err)
    !! Get the value(s) of an option in the parser by name (string/vector)
    !!
    !! On error, the output vector is not allocated.
    OBJECT(cfgparser), INTENT(in)                            :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)                             :: name   !! (Full) Name of the option to get
    CHARACTER(len=*), INTENT(out), DIMENSION(:), ALLOCATABLE :: output !! Output values
    TYPE(error) :: err                                                 !! Error status of the method (see [[cfgparser(type):get_value(bound)]] documentation)
    LOGICAL :: ok
    INTEGER :: idx
    err = noerror
    idx = cp_get_opt_idx(this,name)
    IF (idx == -1) THEN
      err = error("Option "//TRIM(name)//" does not exist",-7)
      RETURN
    ENDIF
    IF (words_length(this%options(idx)%values) == 0) THEN
      err = error("Option "//TRIM(name)//" has no value",-6)
    ELSE
      ok = words_to_vector(this%options(idx)%values,output)
    ENDIF
    RETURN
  END FUNCTION cp_get_sv_ve

  FUNCTION cp_set_rv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (real/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    REAL(kind=4), INTENT(in)         :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,to_string(input))
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_rv_sc

  FUNCTION cp_set_dv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (double/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    REAL(kind=8), INTENT(in)         :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,to_string(input))
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_dv_sc

  FUNCTION cp_set_iv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (double/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    INTEGER, INTENT(in)              :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,to_string(input))
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_iv_sc

  FUNCTION cp_set_lv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (logical/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    LOGICAL, INTENT(in)              :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,to_string(input))
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_lv_sc

  FUNCTION cp_set_cv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (complex/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    COMPLEX, INTENT(in)              :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,to_string(input))
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_cv_sc

  FUNCTION cp_set_sv_sc(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (string/scalar)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout) :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)     :: name   !! (Full) Name of the option to set
    CHARACTER(len=*), INTENT(in)     :: input  !! Input value
    LOGICAL, INTENT(in), OPTIONAL    :: create !! .true. to create option if it does not exist (default to false).
    TYPE(error)                      :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    CALL words_append(values,input)
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_sv_sc

  FUNCTION cp_set_rv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (real/vector)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)       :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)           :: name   !! (Full) Name of the option to get
    REAL(kind=4), INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL          :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                            :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,to_string(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_rv_ve

  FUNCTION cp_set_dv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (double/vector))
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)       :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)           :: name   !! (Full) Name of the option to get
    REAL(kind=8), INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL          :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                            :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,to_string(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_dv_ve

  FUNCTION cp_set_iv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (integer/vector)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)  :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)      :: name   !! (Full) Name of the option to get
    INTEGER, INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL     :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                       :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,to_string(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_iv_ve

  FUNCTION cp_set_lv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (logical/vector)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)  :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)      :: name   !! (Full) Name of the option to get
    LOGICAL, INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL     :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                       :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,to_string(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_lv_ve

  FUNCTION cp_set_cv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (complex/vector)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)  :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)      :: name   !! (Full) Name of the option to get
    COMPLEX, INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL     :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                       :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,to_string(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_cv_ve

  FUNCTION cp_set_sv_ve(this,name,input,create) RESULT(err)
    !! Set new value for the given option by name (complex/vector)
    !!
    !! If _create_ is given to .true., the method will add a new option if it does not exist in
    !! the parser.
    !! In such case, an error ((-9, invalid name)) is raised if the option name is not valid.
    !!
    !! In other case, if the option is not defined in the parser the error status is set to -7.
    OBJECT(cfgparser), INTENT(inout)           :: this   !! Cfgparser object
    CHARACTER(len=*), INTENT(in)               :: name   !! (Full) Name of the option to get
    CHARACTER(len=*), INTENT(in), DIMENSION(:) :: input  !! Input values
    LOGICAL, INTENT(in), OPTIONAL              :: create !! .true. to create option if it does not exist (default to false)
    TYPE(error)                                :: err    !! Error status
    LOGICAL                       :: zcreate
    INTEGER                       :: i,idx
    CHARACTER(len=:), ALLOCATABLE :: sname,pname
    TYPE(words) :: values
    zcreate = .false. ; IF (PRESENT(create)) zcreate = create
    err = noerror
    idx = cp_get_opt_idx(this,name)
    DO i=1,SIZE(input) ; CALL words_append(values,trim(input(i))); ENDDO
    IF (idx == -1) THEN
      IF (zcreate) THEN
        err = op_split_name(name,sname,pname)
        IF (err == 0) err = cp_add_opt(this,sname,pname,values)
      ELSE
        err = error("Option "//TRIM(name)//" does not exist",-7)
      ENDIF
    ELSE
      this%options(idx)%values = values
    ENDIF
    CALL words_clear(values)
  END FUNCTION cp_set_sv_ve

  ! i/o functions
  !--------------

  RECURSIVE FUNCTION read_include(parser,path,isec,ipaths,override) RESULT(err)
    !! Read and parse an included configuration file (internal)
    !! @note
    !! On error, the cfgparser object is left unchanged.
    TYPE(cfgparser), INTENT(inout)        :: parser
      !! A cfgparser object that will store the configuration
    CHARACTER(len=*), INTENT(in)          :: path
      !! A string with the path of the input file to read
    CHARACTER(len=st_slen), INTENT(inout) :: isec
      !! Current section name
    TYPE(words), INTENT(inout)            :: ipaths
      !! List of paths of already included files
    LOGICAL, INTENT(in), OPTIONAL         :: override
      !! An optional boolean flag with .true. to override previous value of duplicated options instead of raising an error.
    TYPE(error) :: err
      !! An error with the first error encountered
    TYPE(option)                  :: curopt
    LOGICAL                       :: zoverride,ok,has_opt
    INTEGER                       :: lineno,lu,i
    CHARACTER(len=2), PARAMETER   :: space = CHAR(32)//","    ! check if , is really wanted... A: YES space are the delimiter of the words internal object !
    CHARACTER(len=2), PARAMETER   :: blanks = CHAR(9)//CHAR(32) ! currently not used because blanks truncate.
    CHARACTER(len=15)             :: sln
    CHARACTER(len=:), ALLOCATABLE :: fulp,dirp,basp
    CHARACTER(len=:), ALLOCATABLE :: curval,ipath
    CHARACTER(len=:), ALLOCATABLE :: line,name,value
    INTEGER, PARAMETER            :: cfg_UNKNOWN = -1, &
                                     cfg_SECTION =  0, &
                                     cfg_OPTION  =  1

    zoverride = .false. ; IF (PRESENT(override)) zoverride = override
    ! initialize local variables
    curval = '' ; line   = '' ; name  = '' ; value = ''
    lineno = 0  ; lu = free_lun()
    IF (LEN_TRIM(isec) == 0) isec = "__default__"
    i = INDEX(TRIM(path),"/",.true.)
    IF (i == 0) THEN
      fulp = fs_realpath("./"//TRIM(ADJUSTL(path)))
    ELSE
      fulp = fs_realpath(TRIM(ADJUSTL(path)))
    ENDIF
    basp = fs_basename(fulp)
    dirp = fs_dirname(fulp)
    ! check for file
    INQUIRE(FILE=TRIM(path),EXIST=ok)
    IF (.NOT.ok) THEN
      err = error(TRIM(path)//": no such file",-11)
      RETURN
    ENDIF
    ! check for lun
    IF (lu == -1) THEN ; err = error("No available logical unit",-12) ; RETURN ; ENDIF
    OPEN(lu,FILE=TRIM(path),STATUS='old',ACTION='READ')
    DO WHILE(readline(lu,line))
      lineno = lineno + 1
      WRITE(sln,'(I15)') lineno ; sln=ADJUSTL(sln)
      ! comment or blank line ?
      IF (is_comment(line,ipath)) THEN
        ! check for includes
        IF (LEN(ipath) > 0) THEN
          ! 1) get relative path
          ipath = fs_relpath(ipath,dirp)
          ! 2) compute asbolute path
          ipath = TRIM(dirp)//"/"//TRIM(ipath)
          ipath = fs_realpath(ipath)
          IF (.NOT.check_include(ipaths,ipath)) THEN
            ipath = fs_basename(ipath)
            err = error(basp//'(L'//TRIM(sln)//"): Circular include &
                        &reference to "//ipath,-14)
            EXIT
          ENDIF
          IF (op_valid(curopt) .AND. LEN(curval) > 0) THEN
            CALL words_extend(curopt%values,TRIM(ADJUSTL(curval)),space,.true.,.true.)
            IF (zoverride) THEN
              err = cp_update_opt(parser,curopt%section,curopt%name,curopt%values)
            ELSE
              err = cp_add_opt(parser,curopt%section,curopt%name,curopt%values)
            ENDIF
            CALL op_clear(curopt); curval = ''
          ENDIF
          err = read_include(parser,ipath,isec,ipaths,zoverride)
          IF (err /= 0) EXIT
        ENDIF
        CYCLE
      ENDIF
      ! continuation line ?
      IF (SCAN(line(1:1),blanks) /= 0 .AND. op_valid(curopt)) THEN
          IF (LEN(curval) == 0) THEN
            curval = strip_comment(line)
          ELSE
            curval = curval//CHAR(32)//strip_comment(line)
          ENDIF
      ELSE
       ! 1. Remove comment part and left adjust line
       line = strip_comment(line)
       ! a section header or option header?
       SELECT CASE (get_kind(line,name,value))
         CASE(cfg_SECTION)
           ! 1. add current value to current option (if any)
           ! 2. update "isec" variable
           IF (op_valid(curopt)) THEN
              IF (LEN(curval) > 0) &
              CALL words_extend(curopt%values,TRIM(ADJUSTL(curval)),space,.true.,.true.)
              IF (zoverride) THEN
                err = cp_update_opt(parser,curopt%section,curopt%name,curopt%values)
              ELSE
                err = cp_add_opt(parser,curopt%section,curopt%name,curopt%values)
              ENDIF
           ENDIF
           CALL op_clear(curopt) ; curval = ''
           IF (cfg_has_section(parser,name) .AND. &
               TRIM(name)/="__default__"    .AND. &
               .NOT.zoverride) THEN
             err = error(basp//'(L'//TRIM(sln)//"): Duplicate section '"//name,-8)
             EXIT
           ENDIF
           isec = TRIM(name)
         CASE(cfg_OPTION)
           ! 1. add current value to current option (if any)
           ! 2. search for option in cursect:
           !    --> duplicate option error if it exists
           !    --> create new option if it does not exist (using curopt)
           ! 3. curval is set to value
           ! 4. update curval
           IF (op_valid(curopt)) THEN
              IF (LEN(curval) > 0) &
              CALL words_extend(curopt%values,TRIM(ADJUSTL(curval)),space,.true.,.true.)
              IF (zoverride) THEN
                err = cp_update_opt(parser,curopt%section,curopt%name,curopt%values)
              ELSE
                err = cp_add_opt(parser,curopt%section,curopt%name,curopt%values)
              ENDIF
           ENDIF
           CALL op_clear(curopt) ; curval = ''
           has_opt = cfg_has_option(parser,TRIM(isec)//"/"//TRIM(name))

           IF (has_opt.AND..NOT.zoverride) THEN
             ! it is an error: no duplicate allowed
             err = error(basp//'(L'//TRIM(sln)//"): Duplicate option '"//TRIM(name)//"' in "//isec,-8)
             EXIT
           ENDIF
           curopt%name = TRIM(name)
           curopt%section = TRIM(isec)
           CALL words_clear(curopt%values)
           curval = value
         CASE(cfg_UNKNOWN)
           ! unknown handles also invalid name: it is a critical error
           IF (err == -9) EXIT
       END SELECT
      ENDIF
    ENDDO
    IF (op_valid(curopt)) THEN
      IF (LEN(curval) > 0) &
      CALL words_extend(curopt%values,TRIM(ADJUSTL(curval)),space,.true.,.true.)
      IF (zoverride) THEN
        err = cp_update_opt(parser,curopt%section,curopt%name,curopt%values)
      ELSE
        err = cp_add_opt(parser,curopt%section,curopt%name,curopt%values)
      ENDIF
    ENDIF
    CALL op_clear(curopt) ; curval = ''

    CLOSE(lu)

  CONTAINS
    FUNCTION get_kind(string,name,value) RESULT(kind)
      !! Split input line and attempt to guess its relevant kind of statement
      !!
      !! The input line is searched for section header format or option assignment.
      !!
      !! - If line begins with '[', has ']' and no '=#' between '[' and ']'
      !!   it is a section header.
      !! - Otherwise, if line has '=', without '#' before '=', it is an option.
      !!
      !! Then the method returns an integer with the kind flag of the statement which is one of
      !! -1 (cfg_UNKNOWN), 0 (cfg_SECTION) or 1 (cfg_OPTION).
      CHARACTER(len=*), INTENT(in)               :: string  !! Input string to process
      CHARACTER(len=:), INTENT(out), ALLOCATABLE :: name, & !! Name of the relevant option/section if any, otherwise an empty string.
                                                    value   !! Value of the relevant option (if any), otherwise an empty string
      INTEGER :: kind                                       !! An integer with the kind of statement.
      CHARACTER(len=:), ALLOCATABLE :: copy
      INTEGER                       :: bi,ei
      CHARACTER(len=2), PARAMETER   :: quotes=CHAR(34)//CHAR(39)
      kind = cfg_UNKNOWN
      ! get a trimmed (and left adjusted) copy
      copy = TRIM(string)
      ! Is it a section ?
      !   ---> search for subscripts of '[' and ']'
      !   ---> check that '[' is 1st char and ']' is last char
      bi = INDEX(copy,'[') ; ei = INDEX(copy,']')
      IF (bi == 1 .AND. ei == LEN(copy) .AND. bi < ei) THEN
        ! it is a section header
        kind = cfg_SECTION
        ! get section name: adjust and trim to remove extra blank spaces
        name = TRIM(ADJUSTL(copy(bi+1:ei-1)))
        ! hack cfg_check_name: append '/a' to get a valid option part to test
        IF (TRIM(name) /= "__default__" .AND. .NOT.cfg_check_name(name//"/a")) THEN
          kind = cfg_UNKNOWN
          err = error("Invalid section name ("//name//")",-9)
          RETURN
        ENDIF
        value = ''
      ELSE
        ! Is it an option ?
        !   --> search for '=' and check if it is set before
        !       1st quote (if any)
        bi = INDEX(copy,"=")
        ! search for quotes
        ei = SCAN(copy,quotes) ; IF (ei==0) ei = LEN(copy)+1
        IF (bi /= 0 .AND. bi < ei) THEN
          kind = cfg_OPTION
          name = to_lower(TRIM(copy(1:bi-1)))
          IF (.NOT.cfg_check_name(name)) THEN
            kind = cfg_UNKNOWN
            err = error("Invalid option name ("//TRIM(name)//")",-9)
            RETURN
          ENDIF
          IF (bi == LEN(copy)) THEN
            value = ''
          ELSE
            value = TRIM(copy(bi+1:))
          ENDIF
        ELSE
          ! not an option and not a section: malformed statement !
          err = error('Malformed statement at line '//TRIM(sln),-13)
        ENDIF
      ENDIF
      RETURN
    END FUNCTION get_kind

    FUNCTION strip_comment(line) RESULT(stripped)
      !! Replace comments part of a string by blank spaces
      !! The method replaces every characters after '#' (included) by spaces.
      !! @note
      !! Output string is also left adjusted, thus only trailing blank can be present.
      CHARACTER(len=*), INTENT(in) :: line !! A string to process
      CHARACTER(len=LEN(line)) :: stripped !! A string of same length than 'line' but without comment(s)

      INTEGER :: idx
      stripped = ADJUSTL(line)
      idx = INDEX(stripped,"#")
      IF (idx > 0) stripped(idx:) = REPEAT(CHAR(32),LEN(line)-idx+1)
      RETURN
    END FUNCTION strip_comment

    FUNCTION readline(lun,string) RESULT(not_eof)
      !! Read a complete line
      !!
      !! Each time it is CALLed, the function reads a complete of the file opened in 'lun' logical
      !! unit and returns .false. if EOF has been reached, .true. otherwise.
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
      INTEGER, INTENT(in)                        :: lun     !! Logical unit with the opened file to read.
      CHARACTER(len=:), INTENT(out), ALLOCATABLE :: string  !! Output processed line
      LOGICAL                                    :: not_eof !! .true. if EOF has NOT been reached yet, .false. otherwise
      CHARACTER(len=50) :: buf
      INTEGER           :: e,sz
      not_eof = .true. ; string = ''
      DO
        READ(lun,'(a)',ADVANCE="no",SIZE=sz,IOSTAT=e) buf
        IF (e == IOSTAT_END) THEN
          not_eof = .false.
          IF (sz > 0) THEN
            string=string//buf(1:sz)
          ENDIF
          EXIT
        ELSE IF (e == IOSTAT_EOR) THEN
          string = string//buf(1:sz)
          EXIT
        ELSE
          string = string//TRIM(buf)
        ENDIF
      ENDDO
    END FUNCTION readline

    FUNCTION is_comment(str,incpath) RESULT(res)
      !! Check if line is a comment or an empty string
      !! @warning
      !! Currently, if an '#include' statement is found, the method assumes a single path is set after the directive.
      CHARACTER(len=*), INTENT(in)               :: str
        !! The string to check
      CHARACTER(len=:), INTENT(out), ALLOCATABLE :: incpath
        !! A string with the filepath to be included if '#include' statement is found, empty string otherwise
      LOGICAL :: res
        !! .true. if line is a comment or an empty string, .false. otherwise
      CHARACTER(len=:), ALLOCATABLE :: copy
      res = .false. ; incpath = ''
      copy = TRIM(ADJUSTL(str))
      IF (LEN(copy) == 0) THEN
        res = .true.
      ELSE IF (INDEX(copy,"#") == 1) THEN
        res = .true.
        ! search for include statement
        ! IMPORTANT: assume that there is only a path after include statement
        IF (INDEX(copy,"#include ") == 1) incpath = remove_quotes(TRIM(ADJUSTL(copy(10:))))
      ENDIF
      RETURN
    END FUNCTION is_comment

    FUNCTION check_include(list,incpath) RESULT(ok)
      !! Check if path is not in list
      !! @note
      !! If path is not in list it is added to the list.
      TYPE(words), INTENT(inout)   :: list    !! A list of paths
      CHARACTER(len=*), INTENT(in) :: incpath !! Path to check in list
      LOGICAL :: ok                           !! .true. if 'path' is __not__ in list, .false. otherwise
      CALL words_reset(list)
      ok = .true.
      DO WHILE(words_valid(list))
        IF (TRIM(incpath) == TRIM(words_current(list))) THEN
          ok = .false. ; EXIT
        ENDIF
        CALL words_next(list)
      ENDDO
      IF (ok) CALL words_append(list,TRIM(incpath))
    END FUNCTION check_include

  END FUNCTION read_include

  ! insertion sort... internal

  SUBROUTINE insertionSort(opts)
    !! Sort an array of Options using insertion sort algorithm
    TYPE(option), INTENT(inout), DIMENSION(:) :: opts !! Array to sort.
    TYPE(option) :: temp
    INTEGER :: i, j
    DO i = 2, SIZE(opts)
      j = i - 1
      temp = opts(i)
      DO WHILE (j>=1) ! .AND. op_greater_than(opts(j),temp))
        IF (op_greater_than(opts(j),temp)) THEN
        opts(j+1) = opts(j)
        j = j - 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      opts(j+1) = temp
      CALL op_clear(temp)
    ENDDO
  END SUBROUTINE insertionSort

END MODULE CFGPARSE

