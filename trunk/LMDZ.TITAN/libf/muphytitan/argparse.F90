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

!! file: argparse.F90
!! summary: Command-line parser source file.
!! author: J. Burgalat
!! date: 2013-2015,2017

#include "defined.h"

MODULE ARGPARSE
  !> Command-line parsing module
  !!
  !! Here are described all the public members of argparse module.
  !! For your own sanity, private methods that call Ancient Gods powers through
  !! evil black magic rituals are not described here.
  !! 
  !! If you only wish to have an overview of argparse usage, you'd better go
  !! [here](|url|/page/swift/p01_argparse.html).

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT
  USE ERRORS
  USE FSYSTEM, ONLY : fs_termsize
  USE STRING_OP, getpar => format_paragraph, splitstr  => format_string  
  IMPLICIT NONE

  PRIVATE

  ! Public members and imported features
  ! from errors: export everything (operator are exported latter)
  PUBLIC :: noerror,error, error_to_string,aborting
  ! from strings
  PUBLIC :: st_slen, st_llen
  PUBLIC :: stderr, stdout
  ! argparse module
  PUBLIC :: new_argparser,              &
            argparser_clear,            &
            argparser_add_option,       &
            argparser_add_positionals,  &
            argparser_throw_error,      & 
            argparser_parse,            &
            argparser_help,             &
            argparser_get_positional,   &
            argparser_get_value,        &
            argparser_reset_values,     &
            argparser_found,            &
            argparser_get_num_values,   &
            argparser_found_positional, &
            argparser_get_num_positional

  PUBLIC :: OPERATOR(==), OPERATOR(/=), ASSIGNMENT(=)

  ! ===========================
  ! PARAMETERS (INTRISIC TYPES)
  ! ===========================

  INTEGER, PARAMETER, PUBLIC :: ap_string  = st_string 
    !! String value type identifier.
  INTEGER, PARAMETER, PUBLIC :: ap_complex = st_complex
    !! Complex value type identifier.
  INTEGER, PARAMETER, PUBLIC :: ap_logical = st_logical
    !! Logical value type identifier.
  INTEGER, PARAMETER, PUBLIC :: ap_integer = st_integer
    !! Integer value type identifier.
  INTEGER, PARAMETER, PUBLIC :: ap_real = st_real
    !! Real value type identifier.

  !> List of all available actions

  INTEGER, PARAMETER, PUBLIC :: ap_store  = 1 
    !! store action ID : Each time the option is seen, values are replaced.
  INTEGER, PARAMETER, PUBLIC :: ap_append = 2 
    !! append action ID : Each time the option is seen, values are appended. 
  INTEGER, PARAMETER, PUBLIC :: ap_count  = 3 
    !! count action ID : increase a counter each time the option is seen.
  INTEGER, PARAMETER, PUBLIC :: ap_help   = 4 
    !! help action ID : help is requested !

  !> List of all available actions
  INTEGER, DIMENSION(4), PARAMETER, PRIVATE :: ap_actions = (/ap_store,  &
                                                              ap_append, &
                                                              ap_count,  &
                                                              ap_help/)
  !> List of all recognized types by the parser
  INTEGER, DIMENSION(5), PARAMETER, PRIVATE :: ap_types = (/ap_string,  & 
                                                            ap_logical, &
                                                            ap_complex, &
                                                            ap_integer, &
                                                            ap_real/) 
  !> The unknown flag
  !!
  !! This flag is only intended to initialize flags. It is set by default during initialization 
  !! and quielty replaced by default flags, if user does not provide the relevant feature.
  INTEGER, PARAMETER :: ap_undef = -1
    
  !> Add an option to the parser
  !!
  !! ```
  !! FUNCTION argparser_add_option(this,dest,sflag,lflag,type,action,default,nrec,help,meta) RESULT(err)
  !!          argparser_add_option(this,dest,flag,type,action,default,nrec,help,meta) RESULT(err)
  !! ```
  !!
  !! The function defines a new argument based on input parameters, checks it and finally sets it 
  !! in the parser. 
  !! 
  !! In its first version both short (`sflag`) and long (`lflag`) options flags are mandatory. In its second
  !! form, a single flag (`flag`) is expected: the method will automatically deduce if it belongs to short or
  !! a long option flag based on the number of hyphens given.
  !! 
  !! `type` value should be one of the following module constants (which are aliases from [[string_op(module)]]):
  !!
  !! - `ap_string` ([[string_op(module):st_string(variable)]])
  !! - `ap_complex` ([[string_op(module):st_complex(variable)]])
  !! - `ap_logical` ([[string_op(module):st_logical(variable)]])
  !! - `ap_integer` ([[string_op(module):st_integer(variable)]])
  !! - `ap_real` ([[string_op(module):st_real(variable)]])
  !! 
  !! `action` value should be one of the following module constants:
  !!
  !! - [[argparse(module):ap_store(variable)]]
  !! - [[argparse(module):ap_append(variable)]]
  !! - [[argparse(module):ap_count(variable)]]
  !! - [[argparse(module):ap_help(variable)]]
  !!
  !! `nrec` string can take the following forms:
  !!
  !! tag   | description
  !! :---: | : -------------
  !!  "?"  | zero or one argument's value
  !!  "*"  | any number of arguments
  !!  "+"  | one or more argument's value(s)
  !!  "X"  | Exactly X values. Where X is the string representation of an integer (0 is accepted).
  !!
  !! On success, `noerror` is returned. Otherwise -9 error code is returned. Errors are only 
  !! produced by misuse of the function arguments. In such case, the program should be
  !! stopped: note that such error should not occur in _released_ programs.
  INTERFACE argparser_add_option
    MODULE PROCEDURE ap_add_option_1, ap_add_option_2
  END INTERFACE
  
  !> Get positional argument value(s)
  INTERFACE argparser_get_positional
    MODULE PROCEDURE ap_get_positional_sc, ap_get_positional_ve
  END INTERFACE 

  !> Get optional argument value(s)
  !!
  !! ```
  !! FUNCTION argparser_get_value(this,name,output) RESULT(err)
  !! ```
  !! 
  !! This is the generic method that can be used to retrieve any kind of argument value(s) from the parser for a given
  !! argument name (as defined by the `dest` argument of [[argparse(module):argparser_add_option(interface)]].
  !! All the methods have the same dummy arguments only `output` dummy argument differs in type and shape.
  !! 
  !! @note
  !! For string vector, `output` is expected to be an allocatable vector of **assumed length**
  !! strings (thus string length is left to user responsability).
  !! A good compromise for strings length is to use the [[string_op(module):st_slen(variable)]] 
  !! parameter.
  INTERFACE argparser_get_value
    MODULE PROCEDURE ap_get_dv_sc, ap_get_rv_sc, ap_get_iv_sc, ap_get_lv_sc, &
                     ap_get_cv_sc, ap_get_sv_sc, ap_get_dv_ve, ap_get_rv_ve, &
                     ap_get_iv_ve, ap_get_lv_ve, ap_get_cv_ve, ap_get_sv_ve
  END INTERFACE

  !> Interface to [[argparse(module):argc(type)]] getters
  !!
  !! All the functions have the same prototype, only kind and type of arguments changed 
  !! from a function to the other.
  INTERFACE argc_get_value                    
    MODULE PROCEDURE ac_get_dv_sc, ac_get_rv_sc, ac_get_iv_sc, ac_get_lv_sc, &
                     ac_get_cv_sc, ac_get_sv_sc, ac_get_dv_ve, ac_get_rv_ve, &
                     ac_get_iv_ve, ac_get_lv_ve, ac_get_cv_ve, ac_get_sv_ve
  END INTERFACE

  !> argc equality operator
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE ac_equals_arg
  END INTERFACE

  !> argc inequality operator
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE ac_differs_arg
  END INTERFACE

  !> argc assignment statement
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE ac_affect_arg
    MODULE PROCEDURE ap_affect_parser
  END INTERFACE

  !> argc *destructors* interface
  INTERFACE clear_argc
    MODULE PROCEDURE ac_clear_arg_sc, ac_clear_arg_ve
  END INTERFACE


  TYPE, PRIVATE :: argc
    !! Defines a command-line argument.
    !!
    !! An [[argparse(module):argc(type)]] object stores all information about a command-line 
    !! argument, that is:
    !!
    !! - its name
    !! - its optional flags 
    !! - its type
    !! - its action
    !! - and finally its values
    PRIVATE
    INTEGER          :: ptype = ap_logical
      !! Type flag (an integer from enum argparse::ap_types)
    INTEGER          :: paction = ap_store
      !! Action flag (an integer from enum argparse::ap_actions)
    INTEGER          :: nrec = 0
      !! Number of values
    LOGICAL          :: fnd = .false.
      !! A boolean flag set to `.true.` if argument has been found in the command-line.
    TYPE(words)      :: values
      !! Values of the argument
    CHARACTER(len=2) :: sflag = "  "
      !! Short flag option
    TYPE(words)      :: meta
      !! Meta variable name(s) of the argument
#if HAVE_FTNDTSTR 
    CHARACTER(len=:), ALLOCATABLE :: default
      !! Default flag 
    CHARACTER(len=:), ALLOCATABLE :: name
      !! Name of the argument (needed to check and retrieve its value(s))
    CHARACTER(len=:), ALLOCATABLE :: lflag 
      !! Long flag option (st_short_len max chars !)
    CHARACTER(len=:), ALLOCATABLE :: help
      !! Help about the argument
#else
    CHARACTER(len=st_slen) :: default = ""
      !! Default flag 
    CHARACTER(len=st_slen) :: name
      !! Name of the argument (needed to check and retrieve its value(s))
    CHARACTER(len=st_slen) :: lflag = ""
      !! Long flag option (st_short_len max chars !)
    CHARACTER(len = 200)   :: help = ""
      !! Help about the argument
#endif
  END TYPE argc


  TYPE, PUBLIC :: argparser
    !! Command-line parser
    !!
    !! This is the main object of the module. It stores definitions of CLI arguments and 
    !! their value(s) once the command-line have been parsed.
    TYPE(argc), PRIVATE, ALLOCATABLE, DIMENSION(:) :: args
      !! List of defined arguments
    INTEGER, PRIVATE :: nargs = 0
      !! Size of args
    INTEGER, PRIVATE :: parsed = -1
      !! Parsing control flag
      !!
      !! The parsing flag determines if the command line have been parsed :
      !!
      !!   - -1 : not parsed yet
      !!   - 0  : unsuccessfully parsed
      !!   - 1  : successfully parsed
    TYPE(argc), PRIVATE :: posals
      !! Positionals arguments (defined as a single argc object)
    LOGICAL, PRIVATE :: have_posal = .false.
      !! Positional control flag
#if HAVE_FTNDTSTR
    CHARACTER(len=:), PRIVATE, ALLOCATABLE :: usg
      !! Program command usage 
    CHARACTER(len=:), PRIVATE, ALLOCATABLE :: descr
      !! Program help description
    CHARACTER(len=:), PRIVATE, ALLOCATABLE :: eplg 
      !! Program help epilog
#else
    CHARACTER(len=st_llen), PRIVATE :: usg
      !! Program command usage 
    CHARACTER(len=st_llen), PRIVATE :: descr
      !! Program help description
    CHARACTER(len=st_llen), PRIVATE :: eplg 
      !! Program help epilog
#endif
    INTEGER, PRIVATE :: mxhlpos = 20 
      !! Position of the short help for options
    INTEGER, PRIVATE :: width = 0
      !! Maximum width of the help 
    LOGICAL, PRIVATE :: init = .false.
      !! Initialization control flag
#if HAVE_FTNPROC    

    CONTAINS
    PROCEDURE, PRIVATE :: ap_add_option_1
    PROCEDURE, PRIVATE :: ap_add_option_2
    PROCEDURE, PRIVATE :: ap_get_positional_sc
    PROCEDURE, PRIVATE :: ap_get_positional_ve
    PROCEDURE, PRIVATE :: ap_get_dv_sc
    PROCEDURE, PRIVATE :: ap_get_rv_sc
    PROCEDURE, PRIVATE :: ap_get_iv_sc
    PROCEDURE, PRIVATE :: ap_get_lv_sc
    PROCEDURE, PRIVATE :: ap_get_cv_sc
    PROCEDURE, PRIVATE :: ap_get_sv_sc
    PROCEDURE, PRIVATE :: ap_get_dv_ve
    PROCEDURE, PRIVATE :: ap_get_rv_ve
    PROCEDURE, PRIVATE :: ap_get_iv_ve
    PROCEDURE, PRIVATE :: ap_get_lv_ve
    PROCEDURE, PRIVATE :: ap_get_cv_ve
    PROCEDURE, PRIVATE :: ap_get_sv_ve
    PROCEDURE, PRIVATE :: ap_check_state
    PROCEDURE, PUBLIC  :: throw_error        => argparser_throw_error
      !! Throw an error and exit the program
    PROCEDURE, PUBLIC  :: parse              => argparser_parse
      !! Parse the command-line (or the given input string).
    PROCEDURE, PUBLIC  :: help               => argparser_help  
      !! Compute and print help
    PROCEDURE, PUBLIC  :: found              => argparser_found
      !! Check if an optional argument has been found on the command-line
    PROCEDURE, PUBLIC  :: get_num_values     => argparser_get_num_values
      !! Get the actual number of values stored in an argument.
    PROCEDURE, PUBLIC  :: found_positional   => argparser_found_positional
      !! Check if positional argument(s) has been found on the command-line
    PROCEDURE, PUBLIC  :: get_num_positional => argparser_get_num_positional
      !! Get the actual number of values stored as positionals.
    PROCEDURE, PUBLIC  :: add_positionals    => argparser_add_positionals
      !! Add positionals definitions in the parser.
    !> Add optional argument definition in the parser.
    GENERIC, PUBLIC    :: add_option         => ap_add_option_1, &
                                                ap_add_option_2
    !> Get the values of the positionals stored in the parser.                                                
    GENERIC, PUBLIC    :: get_positional     => ap_get_positional_sc, &
                                                ap_get_positional_ve
    !> Get the value(s) of the given argument stored in the parser.
    GENERIC, PUBLIC    :: get_value          => ap_get_dv_sc, &
                                                ap_get_rv_sc, &
                                                ap_get_iv_sc, &
                                                ap_get_lv_sc, &
                                                ap_get_cv_sc, &
                                                ap_get_sv_sc, &
                                                ap_get_dv_ve, &
                                                ap_get_rv_ve, &
                                                ap_get_iv_ve, &
                                                ap_get_lv_ve, &
                                                ap_get_cv_ve, &
                                                ap_get_sv_ve
#endif       
  END TYPE argparser

  CONTAINS

  ! argparser main methods (public)
  ! -------------------------------

  FUNCTION new_argparser(usg, dsc, epg, add_help, width, max_help_pos) RESULT(this) 
    !! Initialize an argparser object.
    !! 
    !! The method initializes (properly) an [[argparse(module):argparser(type)]] object.
    !! Even if all the arguments are optional, it is mandatory to **call** the method 
    !! before using an argparser object.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: usg
      !! An optional string with the command line usage of the program. If it is not given
      !! command-line usage is automatically built from informations set in the parser.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: dsc
      !! An optional string with the short description of the program.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: epg
      !! An optional string with the epilog of the program's help
    LOGICAL, INTENT(in), OPTIONAL          :: add_help
      !! An optional boolean flag with `.true.` to automatically set an option for program's help. 
      !! Note that, the option flags `-h` and `--help` are no more available in that case.
    INTEGER, INTENT(in), OPTIONAL          :: width
      !! An optional integer with the maximum width the help text.
    INTEGER, INTENT(in), OPTIONAL          :: max_help_pos
      !! An optional integer with the maximum position of the help string for each option of 
      !! the program when help is requested. Note that this value is just an indicator. The 
      !! helper computes the minimum position between this value and the maximum length of the 
      !! options flags.
    TYPE(argparser) :: this 
      !! An initialized argparse object.
    INTEGER     :: zh
    TYPE(error) :: err
    ! We always clear the parser
    CALL argparser_clear(this)
    ! Set keywords
    IF (PRESENT(usg)) THEN ; this%usg=usg   ; ELSE ; this%usg=''   ; ENDIF
    IF (PRESENT(dsc)) THEN ; this%descr=dsc ; ELSE ; this%descr='' ; ENDIF
    IF (PRESENT(epg)) THEN ; this%eplg=epg  ; ELSE ; this%eplg=''  ; ENDIF
    CALL fs_termsize(zh,this%width) 
    IF (PRESENT(width)) this%width = MAX(width,50)
    IF(PRESENT(max_help_pos)) this%mxhlpos = MAX(5,max_help_pos)
    this%init = .true. ! before adding option !!!
    IF (PRESENT(add_help)) THEN
      IF (add_help) &
        err = argparser_add_option(this,'help',sflag='-h',lflag='--help', &
                    action=ap_help, help="Print this help and quit")
    ENDIF
    RETURN
  END FUNCTION new_argparser

  SUBROUTINE argparser_clear(this)
    !! Clear the parser
    !! The subroutine is used as finalization subroutine of argparser type.
    !! Once the method is called, the object is no more usable until it is
    !! (re)initialized by calling argparse::new_argparser.
    !! @note If **fccp** has not been built with support for finalization subroutine,
    !! it should be called whenever the argparser object is no more used.
    TYPE(argparser), INTENT(inout) :: this 
      !! An argparser object
    IF (ALLOCATED(this%args)) THEN
      CALL clear_argc(this%args)
      DEALLOCATE(this%args)
    ENDIF
    ! maybe we should should set this%posals as allocatable
    CALL clear_argc(this%posals)
    this%nargs      = 0
    this%parsed     = -1
    this%have_posal = .false.
    this%usg        = ''
    this%descr      = ''
    this%eplg       = ''
    this%mxhlpos    = 20
    this%width      = 0
    this%init       = .false.
  END SUBROUTINE argparser_clear

  SUBROUTINE argparser_reset_values(this)
    !! Reset all arguments values in the parser
    !! The method only deletes arguments value(s) currently stored in the parser.
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    ! Check if list is empty
    IF (this%nargs == 0) RETURN
    CALL words_clear(this%args(:)%values)
    CALL words_clear(this%posals%values)
  END SUBROUTINE argparser_reset_values

  FUNCTION argparser_add_positionals(this,nargs,meta,help) RESULT(err)
    !! Add positional arguments definition to the parser.
    !!
    !! The method initializes the entry for positional arguments in the parser.
    !! Positional arguments are always seen by the parser as strings and the 
    !! default associated action is 'store'.
    OBJECT(argparser), INTENT(inout)                     :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                         :: nargs
      !! A string with the expected number of specified values for the option
    CHARACTER(len=*), INTENT(in), DIMENSION(:), OPTIONAL :: meta
      !! A vector of strings with the the displayed value name(s) of the positionals in the help command
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: help
      !! An optional string with a short description of the positional argument(s) 
    TYPE(error) :: err
      !! Error object with the first error encountered in the process.
    INTEGER                       :: ty,ac
    CHARACTER(len=:), ALLOCATABLE :: sf,lf,de
    err = noerror
    IF (.NOT.this%init) THEN 
      err = error("argparse: parser not initialized yet",-1) 
      RETURN
    ENDIF
    IF (this%have_posal) THEN
      err = error('argparse: positionals arguments already defined',-8)
      RETURN
    ELSE
      this%posals%name = 'positional'
      sf = '  ' ; lf = '' ; ty = ap_string ; ac = ap_store ; de = ''
      IF (PRESENT(help)) THEN
        this%posals%help = TRIM(help)
      ELSE
        this%posals%help = ''
      ENDIF
      IF (PRESENT(meta)) THEN
        err = ac_check_and_set(this%posals,sf,lf,ty,ac,de,nargs,meta,.false.)
      ELSE
        err = ac_check_and_set(this%posals,sf,lf,ty,ac,de,nargs,check_flag=.false.)
      ENDIF
      IF (err /= noerror) THEN
        RETURN
      ENDIF
      this%have_posal = this%posals%nrec /= 0 
    ENDIF
    RETURN
  END FUNCTION argparser_add_positionals

  FUNCTION argparser_parse(this,cmd_line,auto) RESULT(err)
    !! Parse the command line
    !! The method parses either the command-line or the string `cmd_line`
    !! given as optional argument and fills the parser's arguments.
    !! @note
    !! If `cmd_line` is provided it should not contains the name of the program or 
    !! the parsing process will certainly failed: program name will be seen as the 
    !! first positional argument and all tokens of the string will then be seen as 
    !! positional.
    OBJECT(argparser), INTENT(inout)       :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cmd_line
      !! An optional string to parse that substitute for the actual command-line.
    LOGICAL, INTENT(in), OPTIONAL          :: auto
      !! An optional boolean flag with `.true.` to instruct the parser wether to perform 
      !! automatic actions or not when error occur during parsing. If `auto` is enabled,
      !! then the parser dumps program's usage and stops the program on error.
    TYPE(error) :: err
      !! Error object with the first error encountered in the process.
    CHARACTER(len=:), ALLOCATABLE :: cline,z
    LOGICAL                       :: zauto
    LOGICAL                       :: rhelp 
    INTEGER                       :: l 
    TYPE(words)                   :: cmd_tokens
    err = noerror
    IF (.NOT.this%init)  THEN
      err = error("parser not initialized yet",-1) ; RETURN
    ENDIF
    rhelp = .false.
    zauto = .false. ; IF (PRESENT(auto)) zauto = auto
    IF (PRESENT(cmd_line)) THEN
      ALLOCATE(cline,source=cmd_line)
    ELSE
      CALL GET_COMMAND(length=l) 
      ALLOCATE(CHARACTER(len=l) :: z) ; CALL GET_COMMAND(z)
      CALL GET_COMMAND_ARGUMENT(0,length=l)
      IF (l >= LEN_TRIM(z)) THEN ; cline='' ; ELSE ; cline=z(l+1:LEN(z)) ; ENDIF
      DEALLOCATE(z)
    ENDIF
    ! reset parsing status
    this%parsed = -1
    CALL argparser_reset_values(this)
    DO
      ! Do we have to launch the turbine ?
      IF (LEN_TRIM(cline) == 0) THEN
        IF (this%have_posal) err=error("Wrong number of arguments",-16)
        EXIT ! ... No :)
      ELSE
        err = ap_split_cmd(this,cline,cmd_tokens,rhelp) 
        ! we only stops processing if :
        !   - the internal error (string length) is raised
        !   - help flag found AND auto has been set
        IF (err /= noerror .OR. (rhelp.AND.zauto)) EXIT
      ENDIF
      CALL words_reset(cmd_tokens) ! not mandatory... at least theoretically
      ! Parses the options
      err = ap_parse_options(this,cmd_tokens,rhelp) 
      IF (err /= noerror) EXIT
      ! exit loop if help is requested. Parser is not completely filled but we
      ! expect someone to use the help action..
      IF (rhelp) EXIT
      ! Parses positionals
      err = ap_parse_positionals(this,cmd_tokens) 
      EXIT ! A one iterated loop :)
    ENDDO
    IF (err /= 0) THEN
      CALL argparser_reset_values(this)
    ELSE
      this%parsed = 1
    ENDIF
    IF (zauto) THEN
      IF (rhelp) CALL argparser_help(this)
      IF (err /= 0) CALL argparser_throw_error(this,err,2)
    ENDIF
    RETURN 
  END FUNCTION argparser_parse

  SUBROUTINE argparser_help(this)
    !! Print help and exit program
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    CHARACTER(len=:), ALLOCATABLE :: helpme
    !!!! WARNING we set no indication here !!!
    IF (.NOT.this%init) RETURN
    helpme = ap_gen_help(this)

    WRITE(stdout,'(a)') helpme

    CALL argparser_clear(this)
    ! Finally we exit the program
    CALL EXIT(0)
  END SUBROUTINE argparser_help

  SUBROUTINE argparser_throw_error(this,err,exit_id)
    !! Dump error on standard error and exit
    !!
    !! The method performs the following actions:
    !!  
    !! - Print the usage command of the program
    !! - Dump the provided @p error message
    !! - Call parser's clean-up subroutine (if a *cleaner* callback has been given during the
    !!   parser's initialization, see [[argparse(module):new_argparser(function)]] documentation
    !! - Stop the program
    !!
    !! The error message is always printed in standard error output.
    !! @note 
    !! If errors::error::id is 0 the method does nothing.
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    TYPE(error), INTENT(in)          :: err
      !! An error object with the error to print
    INTEGER, INTENT(in), OPTIONAL    :: exit_id
      !! An optional integer with the exit code (default to 2)
    CHARACTER(len=:), ALLOCATABLE    :: pgn
    TYPE(error) :: zerr
    IF (err == 0) RETURN
    zerr = error(err%msg,2)
    IF (PRESENT(exit_id)) zerr%id=exit_id
    CALL ap_format_usage(this)
    WRITE(stderr,'(a)') TRIM(this%usg)//NEW_LINE('A')
    ! clean the parser
    CALL argparser_clear(this)
    pgn = get_progname()
    WRITE(stderr,'(a)') pgn//": error: "//TRIM(err%msg)
    CALL EXIT(err%id)
  END SUBROUTINE argparser_throw_error

  FUNCTION argparser_found(this,argname) RESULT(found)
    !! Check wether an argument has been found in the command-line.
    !! @note 
    !! Keep in mind that arguments in the parser always have a default
    !! value. This method is not intended to check if an argument has a value but
    !! only if it has been seen on the command-line !
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: argname
      !! Name of the argument to check
    LOGICAL :: found
      !! `.true.` if the option has been parsed, `.false.` otherwise
    INTEGER  :: idx
    idx = ap_get_arg_index(this,argname)
    IF (idx == -1) THEN
      found = .false.
    ELSE
      found = ac_found(this%args(idx))
    ENDIF
  END FUNCTION argparser_found

  FUNCTION argparser_get_num_values(this,argname) RESULT(num)
    !! Get the actual number of values stored in an argument.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: argname
      !! Name of the argument to check.
    INTEGER :: num 
      !! The number of actual values stored in the argument
    INTEGER  :: idx
    idx = ap_get_arg_index(this,argname)
    IF (idx == -1) THEN
      num = 0
    ELSE
      num = words_length(this%args(idx)%values)
    ENDIF
  END FUNCTION argparser_get_num_values

  FUNCTION argparser_found_positional(this) RESULT(ret)
    !! Check if positional(s) has been found in the command line.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    LOGICAL :: ret
    !! `.true.` if found, `.false.` otherwise
    TYPE(error) :: err
    ret = .false.
    IF (this%have_posal) THEN
      ret = this%posals%fnd
    ENDIF
  END FUNCTION argparser_found_positional

  FUNCTION argparser_get_num_positional(this) RESULT(ret)
    !! Get the actual number of positional argument values stored in the parser .
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    INTEGER :: ret
      !! The number of actual positionals arguments
    ret = 0 
    IF (this%have_posal) THEN
      ret = words_length(this%posals%values)
    ENDIF
  END FUNCTION argparser_get_num_positional 

  ! argparser private methods
  ! -------------------------

  FUNCTION ap_check_state(this) RESULT(err)
    !! Check current parser state 
    !! The method returns an error based on the current parser's state:
    !! - Parser is ready (0)
    !! - parsing not done yet (-19)
    !! - parsing (already) failed (-20)
    !! - parser is **NOT** initialized (-1)
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    TYPE(error) :: err
      !! Error object with the *status* of the parser 
    err = noerror
    IF (this%parsed == -1) THEN
      err = error("argparse: Command-line not parsed (yet)",-19)
    ELSE IF (this%parsed == 0) THEN
      err = error("argparse: command-line parsing failed",-20)
    ELSE IF (.NOT.this%init) THEN
      err = error("argparse: parser not initialized yet",-1) 
    ENDIF
    RETURN
  END FUNCTION ap_check_state

  SUBROUTINE ap_append_arg(this,arg)
    !! Append an argument to the list of arguments.
    TYPE(argparser), INTENT(inout) :: this
      !! An argparser object reference
    TYPE(argc), INTENT(in)         :: arg
      !! An [[argparse(module):argc(type)]] object with the argument to add
    INTEGER                               :: i
    TYPE(argc), ALLOCATABLE, DIMENSION(:) :: tmp
    TYPE(error) :: err 
    IF (.NOT.this%init)  THEN
      err = error("parser not initialized yet",-1) ; RETURN
    ENDIF
    ! Empty list : we create it
    IF (this%nargs == 0) THEN
      ALLOCATE(this%args(1))
      this%args(1) = arg 
      this%nargs = 1 
      RETURN
    ENDIF
    ! Adds a new argument to the vector of arguments
    ! we will make a test with move_alloc but i'm not sure it does everything
    ! the way we want (keep in mind that there is some pointer in argc members
    ! !!!)
    ALLOCATE(tmp(this%nargs))
    DO i=1,this%nargs ; tmp(i) = this%args(i) ; ENDDO
    CALL clear_argc(this%args) 
    DEALLOCATE(this%args) 
    this%nargs = this%nargs+1 ; ALLOCATE(this%args(this%nargs))
    DO i=1,this%nargs-1 ; this%args(i) = tmp(i) ; ENDDO
    CALL clear_argc(tmp)
    DEALLOCATE(tmp)
    this%args(i) = arg
    RETURN
  END SUBROUTINE ap_append_arg

  FUNCTION ap_get_arg_index(this, name, sflag, lflag)  RESULT(idx)
    !! Get an argument by name or option flag in the parser.
    TYPE(argparser), INTENT(in), TARGET    :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name
      !! An optional string with the name of the argument
    CHARACTER(len=*), INTENT(in), OPTIONAL :: lflag
      !! An optional string with the long option flag of the argument
    CHARACTER(len=2), INTENT(in), OPTIONAL :: sflag
      !! An optional string with the short option flag of the argument
    INTEGER :: idx
      !! Index of the argument in the internal vector of argument or -1 if no argument is found.
    INTEGER                       :: i,nn,nl,ns
    CHARACTER(LEN=:), ALLOCATABLE :: lna, llf
    CHARACTER(LEN=2)              :: lsf
    idx = -1
    IF (.NOT.this%init) RETURN
    lna="" ; IF(PRESENT(name)) lna = name   ; nn = LEN_TRIM(lna)
    lsf="" ; IF(PRESENT(sflag)) lsf = sflag ; ns = LEN_TRIM(lsf)
    llf="" ; IF(PRESENT(lflag)) llf = lflag ; nl = LEN_TRIM(llf)
    IF (nn == 0 .AND. ns == 0 .AND. nl == 0) RETURN
    ! empty parser
    IF (this%nargs == 0) RETURN
    DO i=1, this%nargs 
      IF ((nn /= 0 .AND. TRIM(this%args(i)%name)  == TRIM(lna)) .OR. &
          (ns /= 0 .AND. TRIM(this%args(i)%sflag) == TRIM(lsf)) .OR. &
          (nl /= 0 .AND. TRIM(this%args(i)%lflag) == TRIM(llf))) THEN
          idx = i ; RETURN
      ENDIF
    ENDDO
    RETURN
  END FUNCTION ap_get_arg_index


  FUNCTION ap_add_option_1(this,dest,sflag,lflag,type,action,default,nrec,help,meta) RESULT(err)
    !! Add an argument to the parser (interface #1)
    !!
    !! The function defines a new argument based on input parameters, checks it and finally sets it 
    !! in the parser. Both **short and long options flags** are mandatory input arguments of the function.
    !! 
    !!  `type` value should be one of the following module constants (which are aliases from [[string_op(module)]]):
    !!  - ap_string ([[string_op(module):st_string(variable)]])
    !!  - ap_complex ([[string_op(module):st_complex(variable)]])
    !!  - ap_logical ([[string_op(module):st_logical(variable)]])
    !!  - ap_integer ([[string_op(module):st_integer(variable)]])
    !!  - ap_real ([[string_op(module):st_real(variable)]])
    !! 
    !!  `action` value should be one of the following module constants:
    !!  - [[argparse(module):ap_store(variable)]]
    !!  - [[argparse(module):ap_append(variable)]]
    !!  - [[argparse(module):ap_count(variable)]]
    !!  - [[argparse(module):ap_help(variable)]]
    !!
    !!  `nrec` string can take the following forms:
    !!  tag | description
    !!  :-: | : -------------
    !!   ?  | zero or one argument's value
    !!   *  | any number of arguments
    !!   +  | one or more argument's value(s)
    !!   X  | Exactly X values. Where X is the string representation of an integer (0 is accepted).
    !!
    !! See also ap_add_option_2 documentation.
    OBJECT(argparser), INTENT(inout)                     :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                         :: dest
      !! A string with the name of the argument
    CHARACTER(len=2), INTENT(in)                         :: sflag
      !! A two characters string with the short option flag of the argument
    CHARACTER(len=*), INTENT(in)                         :: lflag
      !! A string (3 characters minimum) with the long option flag of the argument
    INTEGER, INTENT(in), OPTIONAL                        :: type
      !! An integer with the type of the argument 
    INTEGER, INTENT(in), OPTIONAL                        :: action
      !! An integer with the action of the argument 
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: default
      !! A string with the default value of the argument if not provided in the CLI
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: nrec
      !! A string with the expected number of specified values for the argument in the CLI. 
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: help
      !! A string with a short description of the argument
    CHARACTER(len=*), DIMENSION(:), INTENT(in), OPTIONAL :: meta
      !! A vector of strings with the displayed value's name(s) of the argument in the help command
    TYPE(error) :: err
      !! Error object with the first error encountered in the process.
    CHARACTER(len=:), ALLOCATABLE :: de,na,he
    INTEGER                       :: ty,ac
    TYPE(argc)                    :: arg
    err = noerror
    he = ""       ; IF (PRESENT(help))    he = TRIM(help)
    na = ""       ; IF (PRESENT(nrec))    na = TRIM(nrec)
    ty = ap_undef ; IF (PRESENT(type))    ty = type
    ac = ap_undef ; IF (PRESENT(action))  ac = action
    de =''        ; IF (PRESENT(default)) de = TRIM(default)
    IF (.NOT.this%init)  THEN
      err = error("argparse: parser not initialized yet",-1) 
      RETURN
    ENDIF
    IF (LEN_TRIM(dest) == 0) THEN
      err = error("argparse: Invalid argument (empty dest)",-2)
      RETURN
    ENDIF
    arg%name = TRIM(dest) ; arg%help = TRIM(he)
    IF (PRESENT(meta)) THEN
      err = ac_check_and_set(arg,sflag,lflag,ty,ac,de,na,meta)
    ELSE
      err = ac_check_and_set(arg,sflag,lflag,ty,ac,de,na)
    ENDIF
    IF (err /= noerror) RETURN
    err = ap_check_in_parser(this,arg)
    IF (err /= noerror) RETURN
    CALL ap_append_arg(this,arg)
    CALL clear_argc(arg)
  END FUNCTION ap_add_option_1

  FUNCTION ap_add_option_2(this,dest,flag,type,action,default,nrec,help,meta) RESULT(err)
    !! Add an argument to the parser (interface #2)
    !!
    !! The function is a wrapper to ap_add_option_1. In this version, 
    !! only one option flag is required. The method only checks for the (trimmed) length of **flag** in
    !! order to choose wether it is a **short** or **long** option flag. Then the function simply calls 
    !! ap_add_option_1 to set the argument.
    !! 
    !! Other dummy arguments have the same meaning as in ap_add_option_1.
    OBJECT(argparser), INTENT(inout)                     :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                         :: dest
      !! A string with the name of the argument
    CHARACTER(len=*), INTENT(in)                         :: flag
      !! A string with either the short or long option flag
    INTEGER, INTENT(in), OPTIONAL                        :: type
      !! A string with the type of the argument
    INTEGER, INTENT(in), OPTIONAL                        :: action
      !! A string with the action of the argument 
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: default
      !! A string with the default value of the argument if not provided in the CLI
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: nrec
      !! A string with the expected number of specified values for the argument in the CLI. 
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: help
      !! A string with a short description of the argument
    CHARACTER(len=*), DIMENSION(:), INTENT(in), OPTIONAL :: meta
      !! A vector of strings with the displayed value's name(s) of the argument in the help command
    TYPE(error) :: err
      !! Error object with the first error encountered in the process.
    CHARACTER(len=:), ALLOCATABLE :: sf,lf,de,na,he
    INTEGER                       :: ty,ac
    err = noerror
    sf = '  ' ; lf = ''
    IF (LEN_TRIM(flag) == 2) THEN
      sf = TRIM(flag)
    ELSE
      lf = TRIM(flag)
    ENDIF
    he = ""       ; IF (PRESENT(help))    he = TRIM(help)
    na = ""       ; IF (PRESENT(nrec))    na = TRIM(nrec)
    ty = ap_undef ; IF (PRESENT(type))    ty = type
    ac = ap_undef ; IF (PRESENT(action))  ac = action
    de = ''       ; IF (PRESENT(default)) de = TRIM(default)
    ! Calling "true" function
    IF (PRESENT(meta)) THEN
      err = ap_add_option_1(this,dest,sf,lf,ty,ac,de,na,he,meta)
    ELSE
      err = ap_add_option_1(this,dest,sf,lf,ty,ac,de,na,he)
    ENDIF
    RETURN
  END FUNCTION ap_add_option_2

  FUNCTION ap_check_in_parser(this,arg) RESULT(err)
    !! Check if an argument is already set in the parser
    !! An argument is assumed to be already set in the parser if either its name,
    !! short or long option flag is already defined (in the parser).
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    TYPE(argc), INTENT(in)        :: arg
      !! An argc object to check
    TYPE(error) :: err
      !! Error object with -8 id if object is already set, no error otherwise.
    INTEGER :: i
    err = noerror
    IF (this%nargs == 0 ) RETURN
    DO i=1,this%nargs
    ! note : we could have use the == operator  but we want to get an specific
    ! error message)
      IF (TRIM(this%args(i)%name) == TRIM(arg%name)) THEN
        err = error("argparse: argument '"//TRIM(arg%name)//"' already defined",-8) ; RETURN
      ENDIF
      IF (LEN_TRIM(this%args(i)%sflag) > 0 .AND. &
          TRIM(this%args(i)%sflag) == TRIM(arg%sflag)) THEN
        err = error("argparse: flag '"//TRIM(arg%sflag)//"' already defined",-8) ; RETURN
      ENDIF
      IF (LEN_TRIM(this%args(i)%lflag) > 0 .AND. &
          TRIM(this%args(i)%lflag) == TRIM(arg%lflag)) THEN
        err = error("argparse: flag '"//TRIM(arg%lflag)//"' already defined",-8) ; RETURN
      ENDIF
    ENDDO
    RETURN
  END FUNCTION ap_check_in_parser

  FUNCTION ap_parse_options(this,cmd,help_req) RESULT(err)
    !! Parse options of the internal command line
    !! This (internal) function manages the parsing of the command line options. 
    OBJECT(argparser), INTENT(inout), TARGET :: this
      !! An argparser object reference
    TYPE(words), INTENT(inout)               :: cmd
      !! A [[string_op(module):words(type)]] object with the command-line to parse
    LOGICAL, INTENT(out)                     :: help_req
      !! An output logical flag with `.true.` if help option has been found, `.false.` otherwise
    TYPE(error) :: err
      !! Error object with the first error encountered in the process
    CHARACTER(len=1), PARAMETER   :: sq = CHAR(39)
    CHARACTER(len=:), ALLOCATABLE :: elt
    INTEGER                       :: arg_idx
    INTEGER                       :: i,nv,ic
    err = noerror ; arg_idx = -1
    help_req = .false.
    DO WHILE(words_valid(cmd))
      ! get current element
      elt = words_current(cmd)
      ! check element kind: is it an option flag (-1/0)or a value (1)?
      ic = ap_check_string(this,elt,arg_idx)
      IF (ic <= 0) THEN
        IF (arg_idx /= -1) THEN
          err = ap_fill_argument(this,cmd,this%args(arg_idx))
          IF (err == 0 .AND. this%args(arg_idx)%paction == ap_help) THEN
            this%parsed = 1
            err = argparser_get_value(this,'help',help_req)
            EXIT
          ENDIF
          IF (err /= 0) EXIT
        ELSE
          err = error("Unknown argument ('"//elt//"')",-9)
          EXIT
        ENDIF
      ELSE
        ! We are in the positionals   !!!
        IF (TRIM(elt) == '--') CALL words_next(cmd)
        EXIT
      ENDIF
      ! iterates to next value
      CALL words_next(cmd)
    ENDDO 

    ! Do we need to check for error here ?
    IF (err /= 0) THEN
      this%parsed = 0
      arg_idx = -1
      RETURN
    ELSE IF (help_req) THEN
      arg_idx = -1
      RETURN
    ENDIF

    ! Must check here for argument with append action if they have the correct
    ! number of argument
    DO i=1,this%nargs
      nv = words_length(this%args(i)%values)
      IF (this%args(i)%fnd                  .AND. &
          this%args(i)%paction == ap_append .AND. &
          this%args(i)%nrec > 0             .AND. &
          nv /= this%args(i)%nrec) THEN
        IF (this%args(i)%nrec < nv) THEN
          err = ac_fmt_val_err(this%args(i),-18) ! extra values
        ELSE 
          err = ac_fmt_val_err(this%args(i),-17) ! missing values
        ENDIF
      ENDIF 
    ENDDO
    IF (err /= 0) this%parsed = 0
    RETURN
  END FUNCTION ap_parse_options 

  FUNCTION ap_parse_positionals(this,cmd) RESULT(err)
    !! Parse positional arguments of the internal command line
    !! This (internal) function manages the parsing of the command line positional arguments. 
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    TYPE(words), INTENT(inout)       :: cmd
      !! A [[string_op(module):words(type)]] object with the command-line to parse
    TYPE(error) :: err
      !! Error object with the first error encountered in the process
    INTEGER                       :: na
    CHARACTER(len=:), ALLOCATABLE :: elt
    err = noerror
    ! cmd iterator is located at the first positional argument

    ! no positional required but current word is valid
    ! Either : no positional required but valid element is present 
    !     Or : positional required but no valid element is present
    IF ((this%have_posal.AND..NOT.words_valid(cmd)) .OR. &
        (.NOT.this%have_posal.AND.words_valid(cmd))) THEN
      err = error("Wrong number of arguments",-16) ; RETURN
    ENDIF
    ! ugly patch : we must clear this%posals%values because of the automatic
    ! setting of default value
    CALL words_clear(this%posals%values)
    this%posals%fnd = .true.
    DO 
      na = words_length(this%posals%values)
      IF (words_valid(cmd)) THEN
        ! Gets the element value
        elt = words_current(cmd)
        ! and add it
        CALL words_append(this%posals%values,elt)
      ELSE
        ! no more elements: check against the number of expected records
        ! 1 or + elements expected but nothing has been saved
        IF (this%posals%nrec == -3 .AND. na > 1) THEN
          err = error("Extra value(s) found (positionals arguments)",-18)
        ELSE IF (this%posals%nrec == -2 .AND. na == 0) THEN
          err = error("Missing values (positionals arguments)",-17)
        ELSE IF (this%posals%nrec > 0 .AND. na /= this%posals%nrec) THEN
          IF (na > this%posals%nrec) THEN
            err = error("Extra value(s) found (positionals arguments)",-18)
          ELSE
            err = error("Missing values (positionals arguments)",-17)
          ENDIF
        ENDIF
        EXIT
      ENDIF
      ! get to the next element
      CALL words_next(cmd)
    ENDDO
    IF (err /= noerror) THEN
      this%posals%fnd = .false.
      this%parsed = 0
    ENDIF
    RETURN
  END FUNCTION ap_parse_positionals

  FUNCTION ap_fill_argument(this,cmd,arg) RESULT(err)
    !! Fill an argument with values
    !!
    !! The function parses remaining parts of the command line from the position of 
    !! the given argument and attempts to retrieve its value(s) (if any). 
    !! Several tests are performed and may raise errors. The function always stops
    !! at the first error encountered which can be one of the following :
    !! - Invalid value type (conversion check failed)
    !! - Missing values (can not get as much values as defined in the argument)
    !! - Extra values found (more values can be retrieved from the command line
    !!   than defined in the argument)
    OBJECT(argparser), INTENT(inout)  :: this
      !! An argparser object reference
    TYPE(words), INTENT(inout)        :: cmd
      !! The command line
    TYPE(argc), INTENT(inout), TARGET :: arg
      !! An argc object reference with the argument currently processed. 
    TYPE(error) :: err
      !! Error object with the first error encountered in the process
    INTEGER                       :: ca, isopt, itmp 
    LOGICAL                       :: ltmp
    CHARACTER(len=:), ALLOCATABLE :: elt
    err = noerror
    ! We always reset arguments value if we are in store mode.
    ! If one wants to set an option several times and keeps track of the
    ! previous values, she must use 'append'
    IF (arg%paction == ap_store) CALL words_clear(arg%values)
    ! Gets values as a function of the expected number of records
    IF (arg%nrec == 0) THEN
      ! argparser_parse (main parsing method) reset all values: for 0 nrec
      ! we should have at least one value saved which is the default.
      ! if it is not the case we set default as values of the argument
      IF (words_length(arg%values) == 0) CALL words_append(arg%values,arg%default)
      ! NOW HANDLING 0 records case:
      ! we do not have to consume any argument but :
      !  - trigger a boolean (set it to .not.default)
      !  - increase a counter (if action == 'count')
      IF (arg%paction == ap_count) THEN
        elt = words_pop(arg%values)
        READ(elt,*) itmp ; itmp = itmp + 1 
        CALL words_append(arg%values,TO_STRING(itmp))
      ELSE IF (arg%ptype == ap_logical) THEN 
        elt = words_pop(arg%values)
        READ(elt,*) ltmp 
        CALL words_append(arg%values,TO_STRING(.NOT.ltmp))
      ENDIF
    ELSE
      ! For any other case, the algorithm is quite simple :
      ! We consume tokens of the command-line until either the end or 
      ! the next option (or '--' separator)
      ! When the exit condition is met we perform some tests based on the
      ! expected argument and sets the returned error object...
      ca = 0 ; IF(arg%paction == ap_append) ca = words_length(arg%values)
      DO
        CALL words_next(cmd) ; elt = words_current(cmd) ; isopt = ap_check_string(this,elt)
        ! We have a "valid" value
        IF (((.NOT.words_valid(cmd)).EQV.(isopt<=0)).AND.TRIM(elt) /= '--') THEN
          ! we have a value that we should check for conversion !
          err = ac_check_value(arg,elt)
          IF (err /= 0) THEN
            err = ac_fmt_val_err(arg,-2)
          ELSE
            CALL words_append(arg%values,elt) ; ca = ca + 1
          ENDIF
        ELSE
          ! 3 cases are possible :
          !    1) we have consumed all argument of command line
          !    2) current argument is not a value !
          !    3) current argument the separator '--' 
          IF (isopt <= 0 .OR. TRIM(elt)=='--') CALL words_previous(cmd)
          IF (arg%nrec == -2 .AND. words_length(arg%values) == 0) & 
            err = ac_fmt_val_err(arg,-17)
          IF (arg%paction /= ap_append .AND. arg%nrec > 0 .AND. ca /= arg%nrec) &
          THEN
            IF (ca > arg%nrec) THEN
              err = ac_fmt_val_err(arg,-18)
            ELSE
              err = ac_fmt_val_err(arg,-17)
            ENDIF
          ENDIF
          EXIT
        ENDIF
        IF (err /= noerror) EXIT
        IF (arg%nrec == -3) EXIT
        IF (ca == arg%nrec) EXIT
      ENDDO
    ENDIF
    arg%fnd = (err == noerror)
    RETURN
  END FUNCTION ap_fill_argument

  FUNCTION ap_split_cmd(this,string,new_cmd,rhelp) RESULT(err) 
    !! Preprocess the command line
    !! The function reads and splits the given string so merged options/values 
    !! are splitted and saves the resulting string elements in a list of words.
    !! @warning
    !! For compilers that does not support allocatable strings in derived types,
    !! computation are highly dependent of [[string_op(module):st_slen(variable):
    !! tokens length are limited by this parameter. 
    IMPLICIT NONE
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: string
      !! A string to process
    TYPE(words), INTENT(out)      :: new_cmd
      !! An output [[string_op(module):words(type)]] object with the processed command line
    LOGICAL, INTENT(out)          :: rhelp   
      !! An output boolean flag with `.true.` if help is requested, `.false.` otherwise
    TYPE(error) :: err
      !! Error object with the first error encountered in the process
    INTEGER                       :: isopt,j,tl,res
    CHARACTER(len=:), ALLOCATABLE :: elt
    TYPE(words)                   :: splitted 
    INTEGER                       :: arg_idx
    err = noerror ; rhelp = .false.
    IF (LEN_TRIM(string) == 0) THEN
      err = error('internal error (empty string)',-255) 
      RETURN
    ENDIF
    ! split input command line in words !
    call words_extend(splitted,string," ",.true.)
    ! reset iterator
    CALL words_reset(splitted)
    DO WHILE(words_valid(splitted))
      elt = words_current(splitted) ; tl = LEN_TRIM(elt)
      isopt = ap_check_string(this,TRIM(elt),arg_idx)
      ! we have a short option : maybe we need to split it
      IF (isopt == -1 ) THEN
        DO j=2,tl
          res = ap_check_string(this,"-"//elt(j:j),arg_idx)
          ! great we have another short option flag
          IF (res == -1 .AND. arg_idx /= -1) THEN ! another short option ! 
            rhelp = (this%args(arg_idx)%paction == ap_help.OR.rhelp)
            ! we must not set some argument's values here  for help argument
            ! if auto is not set the parse method will disable the option 
            ! during next! parsing process
            CALL words_append(new_cmd,"-"//elt(j:j))
          ELSE
            IF (j == 2) THEN ! no more option flag !
              CALL words_append(new_cmd,TRIM(elt(j-1:)))
            ELSE
              CALL words_append(new_cmd,TRIM(elt(j:)))
            ENDIF
            EXIT ! WE MUST EXIT !!!!!!!
          ENDIF
        ENDDO
      ELSE
        IF (isopt == 0.AND.arg_idx /= -1) &
        rhelp = (this%args(arg_idx)%paction == ap_help.OR.rhelp)
        ! we must not set some argument's values here for help argument
        ! if auto is not set the parse method will disable the option during 
        ! next parsing process
        CALL words_append(new_cmd,TRIM(elt))
      ENDIF
      ! Iterates to next word
      CALL words_next(splitted)
    ENDDO
    CALL words_clear(splitted)
    RETURN
  END FUNCTION ap_split_cmd
   
  FUNCTION ap_check_string(this,string,idx) RESULT(ret)
    !! Check if a string is an option flag
    !! The function checks if the given string is either a short, long option flag or a value.
    !! @warning
    !! An empty string is considered as a value !
    OBJECT(argparser), INTENT(in)              :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)               :: string
      !! A string to check
    INTEGER, INTENT(out), OPTIONAL             :: idx
      !! An optional output intger with the index of the afferent argument in the parser (-1 if not found)
    INTEGER :: ret 
      !! Return code with the following possible values:
      !! - -1 if the string is a SHORT option flag
      !! - 0 if the string is a LONG option flag
      !! - 1 if it considered as a value
    CHARACTER(len=52), PARAMETER :: alpha = "abcdefghijklmnopqrstuvwxyz&
                                            &ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    ! return code can be either :
    !  -1 : it's a short option flag
    !   0 : it's a long option flag
    !   1 : it's a value
    ! The combination of the (optional) output index and the return code
    ! allows to check if an option is known or not
    ret = 1 
    ! '--' is special : it is a separator that is seen as a value.
    IF (TRIM(string) == '--') RETURN
    IF (PRESENT(idx)) idx = -1
    ! First we will check the two first characters of the string
    IF (LEN_TRIM(string) >= 2) THEN
      IF (string(1:1) == ACHAR(45)) THEN
        ! It's a long option flag !
        IF (string(2:2) == ACHAR(45)) THEN
          ret = 0
          IF (PRESENT(idx)) idx = ap_get_arg_index(this,lflag=TRIM(string))
        ELSE
          ! It's a short option flag !
          IF (INDEX(alpha, string(2:2)) /= 0) THEN
            ret = -1
            IF (PRESENT(idx)) idx = ap_get_arg_index(this,sflag=string(1:2))
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_check_string

  FUNCTION ap_gen_help(this) RESULT(hlp)
    !! Generate the help string
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    CHARACTER(len=:), ALLOCATABLE :: hlp
      !! The formatted help message
    CHARACTER(len=:), ALLOCATABLE :: copt,spc,opts,text
    INTEGER                       :: i,j,ia,optmw,i1,io,n,zw,zh
    IF (this%width == 0) CALL fs_termsize(zh,this%width) 
    zw = this%width
    ! Sets usage
    CALL ap_format_usage(this,optmw)
    ! computes the help position
    optmw = optmw +2 + 2 ! +2 for indentation + 2 for spaces after option
    optmw = MIN(MAX(5,this%mxhlpos-1),optmw)
    hlp = TRIM(this%usg)//NEW_LINE('A')//NEW_LINE('A')
    ! Sets description
    IF (LEN_TRIM(this%descr) /= 0) & 
      hlp=hlp//getpar(this%descr,zw,2)//NEW_LINE('A')//NEW_LINE('A')
    ! Sets positionals
    IF (this%have_posal) THEN
      IF (LEN_TRIM(this%posals%help) /= 0) THEN
        hlp=hlp//'positionals:'//NEW_LINE('A')
        copt = ac_get_opt_str(this%posals) ; n = LEN_TRIM(copt)
        hlp=hlp//getpar(copt,zw,idt1=2)
        ! we must put help on a new line : no need to extends opt with spaces
        IF (n > optmw-2) THEN
          hlp=hlp//NEW_LINE('A')//REPEAT(CHAR(32),optmw)
        ELSE
          hlp=hlp//REPEAT(CHAR(32),optmw-n-2)
        ENDIF
        spc = REPEAT(CHAR(32),optmw)
        text = getpar(this%posals%help,zw-optmw,0,0)
        j=1
        DO WHILE(j <= LEN(text))
          i = INDEX(text(j:),NEW_LINE('A'))
          IF (i /= 0) THEN
            hlp=hlp//text(j:j+i-1)//spc
            j=j+i
          ELSE
            hlp=hlp//text(j:) ; EXIT
          ENDIF
        ENDDO
        hlp=hlp//NEW_LINE('A')//NEW_LINE('A')
      ENDIF
    ENDIF
    ! Sets options
    IF (this%nargs > 0) THEN
      opts=''
      DO ia=1, this%nargs
        IF (LEN_TRIM(this%args(ia)%help) /= 0) THEN
          copt = ac_get_opt_str(this%args(ia)) ; n = LEN_TRIM(copt)
          opts=opts//getpar(copt,zw,idt1=2)
          ! we must put help on a new line : no need to extends opt with spaces
          IF (n > optmw-2) THEN
            opts=opts//NEW_LINE('A')//REPEAT(CHAR(32),optmw)
          ELSE
            opts=opts//REPEAT(CHAR(32),optmw-n-2)
          ENDIF
          spc = REPEAT(CHAR(32),optmw)
          text = getpar(this%args(ia)%help,zw-optmw,0,0)
          j=1
          DO WHILE(j <= LEN(text))
            i = INDEX(text(j:),NEW_LINE('A'))
            IF (i /= 0) THEN
              opts=opts//text(j:j+i-1)//spc
              j=j+i
            ELSE
              opts=opts//text(j:) ; EXIT
            ENDIF
          ENDDO
          opts=opts//NEW_LINE('A')
        ENDIF
      ENDDO
      IF (LEN_TRIM(opts) > 0) hlp=hlp//'options:'//NEW_LINE('A')//opts
    ENDIF
    IF (LEN_TRIM(this%eplg) /= 0) THEN
      hlp=hlp//NEW_LINE('A')//getpar(this%eplg,zw,2)//NEW_LINE('A')
    ENDIF
    RETURN
  END FUNCTION ap_gen_help

  SUBROUTINE ap_format_usage(this,optmw)
    !! Format command line usage.
    !!
    !! The subroutine creates and formats the command line usage of the 
    !! given argparser object. If [[argparser(type):usg(variable)]] is already set (i.e. not empty) 
    !! then the method only computes the maximum width of the option part of the usage command 
    !! (see `optmw` argument description) if needed. In the other case, the method builds the usage 
    !! line based on the arguments stored in the parser.
    OBJECT(argparser), INTENT(inout) :: this
      !! An argparser object reference
    INTEGER, INTENT(out), OPTIONAL   :: optmw
      !! An optional integer with the maximum length of the option part of the usage command. 
      !! This variable is intended to set a fancy indentation while printing option in the helper.
    CHARACTER(len=:), ALLOCATABLE :: usage, idts, copt,pgn
    INTEGER                       :: i,omw,idt,cl,wopt
    ! Get/Checks/Sets maximum width
    ! get program name
    pgn = get_progname()
    omw = 0
    IF (LEN_TRIM(this%usg) == 0) THEN
      idt = 8 + LEN_TRIM(pgn)
      ALLOCATE(CHARACTER(LEN=idt+1) :: idts)
      idts(1:1) = NEW_LINE('A') ; idts(2:idt+1) = CHAR(32)
      ! Builds the usage string
      usage = "usage: "//TRIM(pgn)//CHAR(32)
      ! Looping over arguments
      cl = idt
      DO i=1,this%nargs
        ! gets max length for help options
        copt = ac_get_opt_str(this%args(i)) ; wopt=LEN_TRIM(copt)
        IF (LEN_TRIM(this%args(i)%help) /= 0.AND. wopt > omw) omw = wopt
        !IF (LEN_TRIM(copt) > omw) omw = LEN_TRIM(copt) ; copt=""
        ! Retrieve current argument usage string
        copt = ac_get_usg_opt_str(this%args(i))
        ! Set a new line if it does not hold width
        IF (cl+LEN_TRIM(copt) > this%width) THEN
          cl = idt ; usage = usage(:)//idts(:)
        ENDIF
        ! Write the argument usage + a blank space !
        usage = usage(:)//TRIM(copt)//CHAR(32)
        cl = cl + LEN_TRIM(copt)+1
      ENDDO
      IF (PRESENT(optmw)) optmw = omw
      ! handling positional
      IF (this%have_posal) THEN
        copt = ac_get_usg_opt_str(this%posals)
        IF (LEN_TRIM(copt) > 0) THEN
          ! print newline if it cannot fit
          IF (cl+LEN_TRIM(copt) > this%width) usage = usage(:)//idts(:)
          usage = usage(:)//TRIM(copt)
        ENDIF
      ENDIF 
      this%usg = usage
    ELSE
      IF (PRESENT(optmw)) THEN
        DO i=1,this%nargs
          copt = ac_get_opt_str(this%args(i)) ; wopt=LEN_TRIM(copt)
          IF (LEN_TRIM(this%args(i)%help) /= 0.AND.wopt > omw) omw = wopt
          copt = ""
        ENDDO
        optmw = omw
      ENDIF
    ENDIF 
  END SUBROUTINE ap_format_usage

  SUBROUTINE ap_affect_parser(this,other)
    !! Argparser assignment operator subroutine
    !! The subroutine assigns `other` to `this`
    TYPE(argparser), INTENT(out) :: this
      !! An argparser object to be assigned
    TYPE(argparser), INTENT(in)  :: other
      !! An argparser object to assign
    INTEGER :: i
    IF (other%nargs > 0) THEN
      this%nargs = other%nargs
      ALLOCATE(this%args(this%nargs))
      DO i=1,this%nargs ; this%args(i) = other%args(i) ; ENDDO
    ENDIF
    this%have_posal = other%have_posal
    this%posals     = other%posals
    this%parsed     = other%parsed
#if HAVE_FTNDTSTR
    IF (ALLOCATED(other%usg))   this%usg   = other%usg
    IF (ALLOCATED(other%descr)) this%descr = other%descr
    IF (ALLOCATED(other%eplg))  this%eplg  = other%eplg
#else
    this%usg   = other%usg
    this%descr = other%descr
    this%eplg  = other%eplg
#endif
    this%mxhlpos = other%mxhlpos
    this%width   = other%width
    this%init    = other%init
  END SUBROUTINE ap_affect_parser

  ! argc methods
  ! ------------

  SUBROUTINE ac_affect_arg(this,other)
    !! Argc object assignment operator subroutine
    !! The subroutine assigns `other` to `this`
    TYPE(argc), INTENT(out) :: this
      !! An argc object to be assigned
    TYPE(argc), INTENT(in)  :: other
      !! An argc object to assign
    this%nrec   = other%nrec
    this%paction = other%paction
    this%ptype   = other%ptype
    this%fnd  = other%fnd
    this%sflag  = other%sflag
#if HAVE_FTNDTSTR
    IF (ALLOCATED(other%name))    this%name    = other%name
    IF (ALLOCATED(other%lflag))   this%lflag   = other%lflag
    IF (ALLOCATED(other%help))    this%help    = other%help
    IF (ALLOCATED(other%default)) this%default = other%default
#else
    this%name    = other%name
    this%lflag   = other%lflag
    this%default = other%default
    this%help    = other%help
#endif
    CALL words_clear(this%values)
    CALL words_clear(this%meta)
    this%values = other%values
    this%meta   = other%meta
    RETURN
  END SUBROUTINE ac_affect_arg

  FUNCTION ac_equals_arg(this,other) RESULT(ret)
    !! Check if two arguments are identical
    !! The method checks if two arguments are equal based on their name, short option flag and long 
    !! option flag.Two arguments are considered equals if at least one of these three members is equal 
    !! and not empty.
    TYPE(argc), INTENT(in) :: this
      !! First argc object to compare
    TYPE(argc), INTENT(in) :: other
      !! Second argc object to compare
    LOGICAL :: ret
      !! `.true.` if the two objects are equal, `.false.` otherwise.
    CHARACTER(len=:), ALLOCATABLE :: tna,ona,tsf,osf,tlf,olf
    INTEGER                       :: tn,ts,tl,on,os,ol
    ret = .false.
    tsf = TRIM(this%sflag) ; osf = TRIM(this%sflag)
#if HAVE_FTNDTSTR
    tna="" ; IF (ALLOCATED(this%name))   tna=TRIM(this%name)
    ona="" ; IF (ALLOCATED(other%name))  ona=TRIM(other%name)
    tlf="" ; IF (ALLOCATED(this%lflag))  tlf=TRIM(this%lflag)
    olf="" ; IF (ALLOCATED(other%lflag)) olf=TRIM(other%lflag)
#else
    tna=TRIM(this%name)  ; ona=TRIM(other%name)
    tlf=TRIM(this%lflag) ; olf=TRIM(other%lflag)
#endif
    tn=LEN_TRIM(tna) ; on=LEN_TRIM(ona) 
    tl=LEN_TRIM(tlf) ; ol=LEN_TRIM(olf)
    ts=LEN_TRIM(tsf) ; os=LEN_TRIM(osf) 
    ! check on name :
    ! Returns True if at least of name, sflag and lflag is set to non-empty
    ! string and is indentical for the two argument !
    ret = ((tn/=0 .AND. on==tn) .AND. tna == ona) .OR. &
          ((tl/=0 .AND. ol==tl) .AND. tlf == olf) .OR. &
          ((ts/=0 .AND. os==ol) .AND. tsf == osf) 
    DEALLOCATE(tna,ona,tsf,osf,tlf,olf)
  END FUNCTION ac_equals_arg

  FUNCTION ac_differs_arg(this,other) RESULT(ret)
    !! Check if two arguments are different
    !! The method checks if two arguments are different based on their names, short option flag 
    !! and long option flag.
    !! @note 
    !! This function is the extact contrary of [[argparse(module):ac_equals_arg(function)]] !
    TYPE(argc), INTENT(in) :: this
      !! First argc object to compare
    TYPE(argc), INTENT(in) :: other
      !! Second argc object to compare
    LOGICAL :: ret
      !! `.true.` if the two objects are different, `.false.` otherwise.
    ret = .NOT. ac_equals_arg(this,other)
  END FUNCTION ac_differs_arg

  SUBROUTINE ac_clear_arg_sc(arg)
    !! argc destructor (scalar)
    !! The subroutine frees all memory used by an argc object and resets its member to 
    !! default values.
    TYPE(argc), INTENT(inout) :: arg
      !! An argc object to free
    CALL words_clear(arg%values)
    CALL words_clear(arg%meta)
    arg%ptype   = ap_logical
    arg%paction = ap_store
    arg%nrec    = 0
    arg%fnd     = .false.
    arg%sflag   = "  "
#if HAVE_FTNDTSTR
    IF (ALLOCATED(arg%default)) DEALLOCATE(arg%default)
    IF (ALLOCATED(arg%name))    DEALLOCATE(arg%name)
    IF (ALLOCATED(arg%lflag))   DEALLOCATE(arg%lflag)
    IF (ALLOCATED(arg%help))    DEALLOCATE(arg%help)
#else
    arg%default = ""
    arg%name    = ""
    arg%help    = ""
    arg%lflag   = ""
#endif
  END SUBROUTINE ac_clear_arg_sc

  SUBROUTINE ac_clear_arg_ve(args)
    !! argc destructor (vector)
    TYPE(argc), INTENT(inout), DIMENSION(:), TARGET :: args
      !! A vector of argc objects to free
    INTEGER :: i
    CALL words_clear(args%values)
    CALL words_clear(args%meta)
    DO i=1, SIZE(args)
      args(i)%ptype   = ap_logical
      args(i)%paction = ap_store
      args(i)%nrec    = 0
      args(i)%fnd     = .false.
      args(i)%sflag    = "  "
#if HAVE_FTNDTSTR
      IF (ALLOCATED(args(i)%default)) DEALLOCATE(args(i)%default)
      IF (ALLOCATED(args(i)%name))    DEALLOCATE(args(i)%name)
      IF (ALLOCATED(args(i)%lflag))   DEALLOCATE(args(i)%lflag)
      IF (ALLOCATED(args(i)%help))    DEALLOCATE(args(i)%help)
#else
      args(i)%default = ""
      args(i)%name    = ""
      args(i)%help    = ""
      args(i)%lflag   = ""
#endif
    ENDDO 
  END SUBROUTINE ac_clear_arg_ve

  FUNCTION ac_found(this) RESULT(yes)
    !! Check if argument has been found in the command line
    TYPE(argc), INTENT(in) :: this
      !! An argc object
    LOGICAL :: yes
      !! `.true.` if the option has been parsed, `.false.` otherwise
    yes = this%fnd
  END FUNCTION ac_found

  FUNCTION ac_get_num_values(this) RESULT(num)
    !! Get the actual number of values stored in the argument
    TYPE(argc), INTENT(in) :: this
      !! An argc object
    INTEGER :: num 
      !! The number of values stored in the argument
    num = words_length(this%values)
  END FUNCTION ac_get_num_values

  FUNCTION ac_get_usg_opt_str(arg) RESULT(line)
    !! Build and format the option string for the usage part of the help message
    !! The function is private part of the help builder. It creates the 
    !! option string part of a given argument.
    TYPE(argc), INTENT(in), TARGET :: arg
      !! An argc object
    CHARACTER(len=:), ALLOCATABLE :: line,meta
      !! Allocated string with the option flags string
    line=""
    IF (LEN_TRIM(arg%sflag) > 0) THEN
      line="["//TRIM(arg%sflag)
    ELSE IF (LEN_TRIM(arg%lflag) > 0) THEN
      line="["//TRIM(arg%lflag)
    ENDIF
    meta = TRIM(words_get(arg%meta,1))
    SELECT CASE (arg%nrec)
      CASE (-3)
        line=line(:)//" ["//meta//"]"
      CASE (-2)
        line=line(:)//CHAR(32)//meta//" [...]"
      CASE (-1)
        line=line(:)//" ["//meta//" [...]]"
      CASE (0)
        ! DO NOTHING BUT MUST BE EXPLICITELY SET
      CASE DEFAULT
        meta = words_to_string(arg%meta," ")
        line = line(:)//CHAR(32)//meta
    END SELECT
    IF (line(1:1) == "[") line=line(:)//"]"
    RETURN
  END FUNCTION ac_get_usg_opt_str

  FUNCTION ac_get_opt_str(arg) RESULT(line)
    !! Build and format the option flag string for the option part of the help message
    !! The function is private part of the help builder. It creates the 
    !! option string part of a given argument for the options part of the help.
    TYPE(argc), INTENT(in), TARGET :: arg
      !! An argc object
    CHARACTER(len=:), ALLOCATABLE :: line
      !! Allocated string with the option flags string
    CHARACTER(len=:), ALLOCATABLE :: values,m
    ! creates values string
    m = TRIM(words_get(arg%meta,1))
    SELECT CASE (arg%nrec)
      CASE(-3)
        values="["//m//"]"
      CASE(-2)
        values=m//" [...]"
      CASE(-1)
        values="["//m//" [...]]"
      CASE (0)
        values=""
      CASE DEFAULT
        values = words_to_string(arg%meta," ")
     END SELECT
     ! build the final string
     ! -s values, --long values
     line=""
     IF (LEN_TRIM(arg%sflag) > 0) THEN
       line=TRIM(arg%sflag)//CHAR(32)//TRIM(values)
       IF (LEN_TRIM(arg%lflag) > 0) line = line(:)//", "
     ENDIF
     IF (LEN_TRIM(arg%lflag) > 0) THEN
         line = line(:)//TRIM(arg%lflag)//CHAR(32)//TRIM(values)
     ENDIF
     ! this line handles positionals :)
     IF (arg%name == 'positional') line = TRIM(values)
    RETURN
  END FUNCTION ac_get_opt_str

  FUNCTION ac_fmt_val_err(arg,id) RESULT(ret)
    !! Format specific error for argc values
    !! The function formats argparse specific errors when extra (missing) values
    !! are (not) set or when given values are not consistent with argument's type.
    !! For each of these errors, the basic error is updated with a more precise 
    !! message. 
    TYPE(argc), INTENT(in) :: arg
      !! An argc object
    INTEGER, INTENT(in)    :: id
      !! Integer with the error id (-3, -15 or -16 in fact)
    TYPE(error) :: ret
      !! Error object with the updated error message.
    CHARACTER(len=:), ALLOCATABLE :: msg,nv
    IF (LEN_TRIM(arg%sflag) /= 0) THEN
      msg=' ('//arg%sflag
    ELSE ! note that normally if no sflag are present lflag should be set
      msg=' ('//TRIM(arg%lflag)
    ENDIF
    nv = to_string(arg%nrec) ; IF (arg%nrec < 0) nv='1'
    ! we only handle ids : -2 (invalid arg), -17 (missing values), -18 (extra values)
    SELECT CASE(id)
      CASE (-2)  ! invalid value: cannot cast ->
        msg = msg//" option expects '"//apt2str(arg%ptype)//"' values)"
      CASE (-17) ! missing values -> -17
        IF (arg%nrec == -2) THEN
          msg = msg//' takes at least '//nv//' value(s))'
        ELSE 
          msg = msg//' takes exactly '//nv//' value(s))'
        ENDIF
      CASE (-18) ! extra values -> -18
        IF (arg%nrec == -3) THEN 
          msg = msg//' takes at most '//nv//' value(s))'
        ELSE
          msg = msg//' takes exactly '//nv//' value(s))'
        ENDIF
      CASE DEFAULT
        msg=''
    END SELECT
    ret = error(msg,id)
  END FUNCTION ac_fmt_val_err

  FUNCTION ac_check_and_set(this,sf,lf,ty,ac,de,na,meta,check_flag) RESULT(ret)
    !! Interface to all argc's member tests
    !! The function calls all the tests to perform on argc members. Some of these tests can 
    !! alter argc's member values to fit argparser's requirements. 
    !!
    !! On success, `noerror` is returned. Otherwise -9 error code is returned. Errors are only 
    !! produced by misuse of the function arguments. In such case, the program should be
    !! stopped: note that such ezrror should not occur in _released_ programs.
    TYPE(argc), INTENT(inout)                            :: this
      !! An argc object
    CHARACTER(len=2), INTENT(in)                         :: sf
      !! A string (2 chars wide) with the short option flag
    CHARACTER(len=*), INTENT(in)                         :: lf
      !! A string with the long option flag
    INTEGER, INTENT(in)                                  :: ty
      !! An integer with the type of the argument
    INTEGER, INTENT(in)                                  :: ac
      !! An integer with the action of the argument
    CHARACTER(len=*), INTENT(in)                         :: de
      !! A string with the default value of the argument
    CHARACTER(len=*), INTENT(in)                         :: na
      !! A string pattern with the expected number of values for the argument
    CHARACTER(len=*), INTENT(in), DIMENSION(:), OPTIONAL :: meta
      !! An optional vector of strings with the Meta-name of the values
    LOGICAL, INTENT(in), OPTIONAL                        :: check_flag
      !! An optional boolean flag hat instructs the method wether to check for option flag 
      !! or not. By default this test is enabled but it should be disabled if one wants to 
      !! check for POSITIONAL arguments as they do not have option flags.
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    LOGICAL :: zflag
    zflag = .true. ; IF (PRESENT(check_flag)) zflag = check_flag
    ! 1) CHECKING FLAG SYNTAX
    IF (zflag) THEN
      ret = ac_check_flags(this,sf,lf) ; IF (ret /= 0) RETURN
    ELSE
      this%sflag = sf ; this%lflag = lf
    ENDIF
    ! 2) CHECKING AND SETTING action, type, default value and number of record
    ret = ac_check_ac_ty_de_na(this,ac,ty,de,na) ; IF (ret /= 0) RETURN
    ! 3) Sets/updates meta name
    IF (PRESENT(meta)) THEN
      CALL ac_set_meta(this,meta)
    ELSE
      CALL ac_set_meta(this)
    ENDIF
    RETURN
  END FUNCTION ac_check_and_set

  FUNCTION ac_check_ac_ty_de_na(this,ac,ty,de,na) RESULT(ret)
    !! Check and set argc's action, type and default value
    !! The method checks if input argument's options are valid and update the argc object 
    !! consequently. 
    !!
    !! On success, `noerror` is returned. Otherwise -9 error code is returned. Errors are only 
    !! produced by misuse of the function arguments.
    TYPE(argc), INTENT(inout)  :: this
      !! An argc object to update
    INTEGER, INTENT(in)          :: ac
      !! An integer with the action to set and check 
    INTEGER, INTENT(in)          :: ty
      !! An integer with the type to set and check
    CHARACTER(len=*), INTENT(in) :: de
      !! A string with the default value to set and check
    CHARACTER(len=*), INTENT(in) :: na 
      !! A string with the expected number of value to set and check
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=:), ALLOCATABLE :: eprf,zna
    TYPE(error)                   :: uty,ina
    ret = noerror
    eprf = 'argparse: '//"Invalid argument `"//TRIM(this%name)//"'"
    uty = error(eprf//" (type)",-9) 
    ina = error(eprf//" (inconsistent nargs)",-9)
    zna = TRIM(na)
    ! Checks action
    IF(ANY(ap_actions == ac).OR.ac == ap_undef) THEN
      this%paction = ac 
    ELSE
      ret = error(eprf//" (action)",-9) ; RETURN 
    ENDIF
    ! Checks and sets type and default as a function of the action
    SELECT CASE(this%paction)
      ! HELP: fixed in any case: 
      CASE(ap_help)
        this%default = 'F'
        this%ptype = ap_logical
        this%nrec = 0
      ! COUNT:
      ! we always use "hard-coded" stuff and do not warn if dev has made
      ! mistakes...
      CASE(ap_count)
        ! default settings of the count action
        this%ptype = ap_integer ; this%default = '0' ; this%nrec = 0
        ! check and set default value
        ret = set_def_val()
      ! STORE, APPEND actions
      CASE (ap_store, ap_append) 
        ! set type 
        IF (ty == ap_undef) THEN
          this%ptype= ap_integer 
        ELSEIF (ANY(ap_types == ty)) THEN
          this%ptype = ty
        ELSE
          ret = uty ; RETURN
        ENDIF
        ! check and set default value
        ret = set_def_val() ; IF (ret /= 0) RETURN
        ! check for nargs (if na is empty then we set "*")
        ret = set_nrec("*") 
      ! UNDEFINED:
      !   -> 1) set to action to store
      !   ->
      ! try to define a logical trigger and modifies it with user-defined
      ! characteristics
      CASE (ap_undef)
        ! 1) always define store action
        this%paction = ap_store 
        ! 2) set type and nrec:
        !    2.1) If type is undef:
        !         - to default value type if default is given
        !         - ap_logical otherwiset type 
        !    2.2) set nrec
        !        - if final type is ap_logical set nrec to 0
        !        - otherwise to *
        If (ty == ap_undef) THEN
          ! no explicit type : define logical trigger first
          this%ptype = ap_logical ; this%nrec = 0 
          ! icheck the default value given
          IF (LEN_TRIM(de) > 0) THEN
            this%ptype = ap_types(string_is(de)+1)
          ELSE
            this%ptype = ap_logical
          ENDIF
          IF (this%ptype == ap_logical) THEN
            ret = set_nrec("0") 
          ELSE
            ret = set_nrec("1") 
          ENDIF
          ret = set_def_val() 
          IF (ret /= 0) RETURN
        ! type is given
        ELSE IF (ANY(ty == ap_types)) THEN
          ! known type given :
          !  check default value and nrec: -> if na not given set "*"
          this%ptype = ty 
          ret = set_def_val() ; IF (ret /= 0) RETURN
          IF (this%ptype == ap_logical) THEN
            ret = set_nrec("0") 
          ELSE
            ret = set_nrec("1")
          ENDIF
        ELSE
          ! unknown type => bad end !
          ret = uty ; RETURN 
        ENDIF
    END SELECT 
    ! set default value as first value if ret is noerror ....
    IF (ret == 0) CALL words_append(this%values,this%default)
    RETURN

    CONTAINS

      FUNCTION set_nrec(base) RESULT(terr)
        !! Check and set argument's expected number of records
        !! The method compares `na` value with the expected and known flags and decides
        !! wether to raise an error or defines nrec member of the argument object. If `na` 
        !! is empty then `base` is used.
        CHARACTER(len=1),INTENT(in) :: base
          !! Base value of nrec
        TYPE(error) :: terr
          !! Error object with the return status of the function
        terr = noerror
        ! check for nargs:
        IF (LEN(zna) == 0) zna=base
        SELECT CASE(zna)
          CASE("*") ; this%nrec = -1
          CASE("+") ; this%nrec = -2
          CASE("?") ; this%nrec = -3
            IF (this%paction == ap_append) terr = ina
          CASE DEFAULT ; this%nrec = 0
            ! check numeric characters
            IF (VERIFY(zna,"0123456789")==0) READ(zna,*) this%nrec
            IF (this%nrec == 0) terr = ina 
        END SELECT
      END FUNCTION set_nrec

      FUNCTION set_def_val() RESULT(terr)
        !! Check and set argument's default value
        !! The method compares `de` value with the type already stored in the argument and 
        !! decides wether to raise an error or to save `de` as argument's default value.
        !! If `de` is empty then it sets a default value according to argument's type.
        TYPE(error) :: terr
          !! Error object with the return status of the function
        INTEGER :: t
        terr = noerror
        IF (LEN_TRIM(de) /= 0) THEN
          this%default = de ; t = string_is(de) 
          IF (t /= this%ptype) THEN
            terr = error(eprf//" (inconsistent default value: expected '"// &
                           TRIM(st_type_names(this%ptype))//"', found '"// &
                           TRIM(st_type_names(t))//"')",-9)
            RETURN
          ENDIF
        ELSE
          SELECT CASE (this%ptype)
            CASE(ap_string)  ; this%default = ''
            CASE(ap_logical) ; this%default = 'F'
            CASE(ap_complex) ; this%default = '(0d0,0d0)'
            CASE(ap_integer) ; this%default = '0'
            CASE(ap_real)    ; this%default = '0d0'
          END SELECT
        ENDIF
        RETURN
      END FUNCTION set_def_val 
  END FUNCTION ac_check_ac_ty_de_na

  FUNCTION ac_check_flags(this,sflag,lflag) RESULT(ret)
    !! Set and check if argument's option flags are valid
    !! The method first sets the given flags in the argc object and then checks if
    !! its option flags are valid :
    !!    - A valid short option flag (`sflag`) starts with `-' followed by the
    !!      regex pattern : @p [a-zA-Z]
    !!    - A valid long option flag (`lflag`) starts with `--' followed by the
    !!      regex pattern : @p [-a-zA-Z_0-9]+
    !! @note
    !! Although the two arguments are mandatory one can be an empty string
    TYPE(argc), INTENT(inout)    :: this
      !! An argc object
    CHARACTER(len=2), INTENT(in) :: sflag
      !! A 2-characters wide string with the short option flag
    CHARACTER(len=*), INTENT(in) :: lflag 
      !! A string (at least 3 characters wide) with the long option flag
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: ic
    ret = noerror
    this%sflag = sflag ; this%lflag = lflag
    ! 1) Check null strings !
    IF (LEN_TRIM(this%sflag) == 0 .AND. LEN_TRIM(this%lflag) == 0) THEN
      ret = error("argparse: Invalid argument (empty option flags)",-2) ; RETURN
    ENDIF
    ! 2) Check short option flag:
    IF (LEN_TRIM(this%sflag) > 0) THEN
      IF (this%sflag(1:1) /= achar(45)) THEN
        ret = error("argparse: Invalid argument (short option)",-2) ; RETURN
      ELSE
        SELECT CASE(iachar(this%sflag(2:2)))
          CASE(97:122,65:90) ; ret = noerror
          CASE DEFAULT
            ret = error("argparse: Invalid argument (short option)",-2) ; RETURN
        END SELECT
      ENDIF
    ENDIF
    ! 3) Check long option flag
    IF (LEN_TRIM(this%lflag) > 0) THEN
      ! We must have at least 3 chars !
      IF (LEN_TRIM(this%lflag) < 3) THEN
        ret = error("argparse: Invalid argument (long option)",-2) ; RETURN
      ELSE IF (this%lflag(1:2) /= achar(45)//achar(45)) THEN
        ret = error("argparse: Invalid argument (long option)",-2) ; RETURN
      ELSE
        DO ic=3,LEN_TRIM(this%lflag)
          SELECT CASE(iachar(this%lflag(ic:ic)))
            !  corresponds to [-a-zA-Z0-9_]
            CASE(45,48:57,65:90,95,97:122) ; ret = noerror
            CASE DEFAULT
              ret = error("argparse: Invalid argument (long option)",-2) ; RETURN
              EXIT
          END SELECT
        ENDDO
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ac_check_flags

  SUBROUTINE ac_set_meta(this, meta)
    !! Set meta-variable of the given argc object
    !! The method set meta-variable in the argc object. If no `meta` are given, the method
    !! uses argument's name to set the values. 
    !! @warning
    !! To be effective, this subroutine must be called after argparse::chk_opt_nargs
    TYPE(argc), INTENT(inout)                            :: this
      !! An argc object reference
    CHARACTER(len=*), INTENT(in), DIMENSION(:), OPTIONAL :: meta
      !! An optional vector of strings with the meta-variable(s) to set
    INTEGER                       :: i,j,ms,blk
    CHARACTER(len=:), ALLOCATABLE :: zmeta
    ! clear meta values (not needed normally)
    CALL words_clear(this%meta)
    IF (PRESENT(meta)) THEN
      SELECT CASE(this%nrec)
        CASE(-3,-2,-1)
          zmeta = to_upper(meta(1))
          blk = INDEX(TRIM(zmeta),CHAR(32)) - 1
          IF (blk <= 0) blk=LEN_TRIM(zmeta)
          CALL words_append(this%meta,zmeta(1:blk))
        CASE(0)
          CALL words_append(this%meta,"")
        CASE DEFAULT
          ms = SIZE(meta) ; j = 0
          DO i=1,this%nrec
            j=j+1 ; IF (j>ms) j=1
            zmeta = to_upper(meta(j))
            blk=INDEX(TRIM(zmeta),CHAR(32))-1 
            IF (blk <= 0) blk=LEN_TRIM(zmeta)
            CALL words_append(this%meta,zmeta(1:blk))
          ENDDO
      END SELECT
    ELSE
      zmeta = to_upper(TRIM(this%name))
      SELECT CASE(this%nrec)
        CASE(-3,-2,-1)
          CALL words_append(this%meta,zmeta)
        CASE DEFAULT
          DO i=1,this%nrec
            CALL words_append(this%meta,zmeta)
          ENDDO
      END SELECT
    ENDIF
    RETURN
  END SUBROUTINE ac_set_meta

  FUNCTION ac_check_value(this,str) RESULT(err)
    !! Check if given string is a valid value for the argument
    TYPE(argc), INTENT(in)       :: this 
      !! An argc object reference
    CHARACTER(len=*), INTENT(in) :: str
      !! A string with the value to check
    TYPE(error) :: err
      !! Error object with the first error encountered in the process
    INTEGER :: ty
    err = noerror
    ! Special conditions for strings: any kind of value is ok if we asked for
    ! strings
    IF (this%ptype == ap_string) RETURN
    ty = string_is(str)
    ! special handling for integer: they can be seen as real
    IF (this%ptype == ap_real) THEN
      IF (ty < 3)  err = error("Cannot cast value",-10)
    ELSE
      IF (ty /= this%ptype) err = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_check_value

  !===========
  !  GETTERS
  !===========

  ! argparser getters
  ! -----------------

  FUNCTION ap_get_positional_sc(this,idx,value) RESULT(ret)
    !! Get positional arguments value at given index
    !! @warning 
    !! On error, the status of `value` is undefined.
    OBJECT(argparser), INTENT(in)              :: this
      !! An argparser object reference
    INTEGER, INTENT(in)                        :: idx
      !! Subscript of the positional argument value to get
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: value
      !! Output raw value of the positional. If `idx` is out of range, `value` is set to an 
      !! empty string.
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    IF (.NOT.this%have_posal) THEN
      ret = error("argparse: No positional argument(s) defined",-7)
    ELSE IF (idx <= 0 .OR. idx > words_length(this%posals%values)) THEN
      ret = error("argparse: index out of range", -3)
    ELSE
      value = words_get(this%posals%values,idx)
      ret = noerror
    ENDIF
    RETURN
  END FUNCTION ap_get_positional_sc

  FUNCTION ap_get_positional_ve(this,values) RESULT(ret)
    !! Get all positional arguments value
    !! @warning 
    !! On error, the status of `values` is undefined.
    OBJECT(argparser), INTENT(in)                                  :: this
      !! An argparser object reference
    CHARACTER(len=st_slen), INTENT(out), ALLOCATABLE, DIMENSION(:) :: values
      !! An allocatable vector of **assumed length** strings with the value(s) of all 
      !! positionals arguments found. 
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    LOGICAL :: ok
    ret = noerror
    IF (.NOT.this%have_posal) THEN
      ret = error("argparse: No positional argument(s) defined",-7)
    ELSE IF (words_length(this%posals%values) == 0) THEN
      ret = error("argparse: No positional argument(s) values found",-6)
    ELSE
      ok = words_to_vector(this%posals%values,values)
    ENDIF
    RETURN
  END FUNCTION ap_get_positional_ve

  FUNCTION ap_get_dv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `REAL(kind=8)` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` value is undefined.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: name
      !! Name of the argument
    REAL(kind=8), INTENT(out)     :: output
      !! A scalar with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_dv_sc

  FUNCTION ap_get_rv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `REAL(kind=4)` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` value is undefined.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: name
      !! Name of the argument
    REAL(kind=4), INTENT(out)     :: output
      !! A scalar with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_rv_sc

  FUNCTION ap_get_iv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `INTEGER` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` value is undefined.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: name
      !! Name of the argument
    INTEGER, INTENT(out)          :: output
      !! A scalar with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_iv_sc

  FUNCTION ap_get_lv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `LOGICAL` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` value is undefined.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: name
      !! Name of the argument
    LOGICAL, INTENT(out)          :: output
      !! A scalar with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_lv_sc

  FUNCTION ap_get_cv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `COMPLEX` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` value is undefined.
    OBJECT(argparser), INTENT(in) :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)  :: name
      !! Name of the argument
    COMPLEX, INTENT(out)          :: output
      !! A scalar with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_cv_sc

  FUNCTION ap_get_sv_sc(this, name, output) RESULT(ret)
    !! Get a scalar `STRING` value from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)              :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)               :: name
      !! Name of the argument
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: output
      !! An allocatable string with the first value of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
    RETURN
  END FUNCTION ap_get_sv_sc

  FUNCTION ap_get_dv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `REAL(kind=8)` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                        :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                         :: name
      !! Name of the argument
    REAL(kind=8), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_dv_ve

  FUNCTION ap_get_rv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `REAL(kind=4)` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                        :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                         :: name
      !! Name of the argument
    REAL(kind=4), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_rv_ve

  FUNCTION ap_get_iv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `INTEGER` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                   :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                    :: name
      !! Name of the argument
    INTEGER, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_iv_ve

  FUNCTION ap_get_lv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `LOGICAL` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                   :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                    :: name
      !! Name of the argument
    LOGICAL, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_lv_ve

  FUNCTION ap_get_cv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `COMPLEX` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                   :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                    :: name
      !! Name of the argument
    COMPLEX, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_cv_ve

  FUNCTION ap_get_sv_ve(this, name, output) RESULT(ret)
    !! Get a vector of `STRING` values from given argument
    !! The error returned by the method can be either:
    !! -  0  : No error
    !! - -1  : parser has not been initialized
    !! - -7  : argument not found (i.e. does not set in the parser)
    !! - -19 : parsing not done yet
    !! - -20 : (previous) parsing failed
    !! - -21 : inconsistent destination type
    !! @note
    !! If no error occured, the function always set a value which is the default value
    !! if the argument has not been parsed. Otherwise, `output` status is undefined.
    OBJECT(argparser), INTENT(in)                            :: this
      !! An argparser object reference
    CHARACTER(len=*), INTENT(in)                             :: name
      !! Name of the argument
    CHARACTER(len=*), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! An allocatable vector of **assumed length** strings with the values of the argument
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    INTEGER :: idx
    ret = ap_check_state(this)
    IF (ret == 0) THEN
      idx = ap_get_arg_index(this,name)
      IF (idx == -1) THEN
        ret = error("argparse: Argument not found ("//TRIM(name)//")",-7)
      ELSE
        ret = argc_get_value(this%args(idx),output)
      ENDIF
    ENDIF
  END FUNCTION ap_get_sv_ve

  ! argc getters
  ! ------------

  !> Gets a scalar @p REAL(kind=8) value from given argument
  !! @param[in,out] this An argc object 
  !! @param[out] output A scalar with the first value of the argument
  !! @return An errors::error object with -21 if the destination variable's type
  !! is inconsistent, 0 otherwise.
  FUNCTION ac_get_dv_sc(this, output) RESULT(ret)
    !! Get a scalar `REAL(kind=8)` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` value is undefined.
    TYPE(argc), INTENT(in)  :: this
      !! An argc object
    REAL(kind=8), INTENT(out) :: output
      !! Output value
    TYPE(error) :: ret
      !! Error object with the -21 if the destination variable's type is inconsistent, 0 otherwise
    ret = noerror
    IF (this%ptype /= ap_real) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got REAL(kind=8))",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      IF (.NOT.from_string(this%default,output)) &
        ret = error("Cannot cast value",-10)
    ELSE
      IF (.NOT.from_string(words_get(this%values,1), output)) &
        ret = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_get_dv_sc

  FUNCTION ac_get_rv_sc(this, output) RESULT(ret)
    !! Get a scalar `REAL(kind=4)` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` value is undefined.
    TYPE(argc), INTENT(in)    :: this
      !! An argc object
    REAL(kind=4), INTENT(out) :: output
      !! Output value
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_real) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got REAL(kind=4))",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      IF (.NOT.from_string(this%default,output)) &
        ret = error("Cannot cast value",-10)
    ELSE
      IF (.NOT.from_string(words_get(this%values,1), output)) &
        ret = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_get_rv_sc

  FUNCTION ac_get_iv_sc(this, output) RESULT(ret)
    !! Get a scalar `INTEGER` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` value is undefined.
    TYPE(argc), INTENT(in) :: this
      !! An argc object
    INTEGER, INTENT(out)   :: output
      !! Output value
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_integer) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got INTEGER)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      IF (.NOT.from_string(this%default,output)) &
        ret = error("Cannot cast value",-10)
    ELSE
      IF (.NOT.from_string(words_get(this%values,1), output)) &
        ret = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_get_iv_sc

  FUNCTION ac_get_lv_sc(this, output) RESULT(ret)
    !! Get a scalar `INTEGER` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` value is undefined.
    TYPE(argc), INTENT(in) :: this
      !! An argc object
    LOGICAL, INTENT(out)   :: output
      !! Output value
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_logical) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got LOGICAL)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      IF (.NOT.from_string(this%default,output)) &
        ret = error("Cannot cast value",-10)
    ELSE
      IF (.NOT.from_string(words_get(this%values,1), output)) &
        ret = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_get_lv_sc

  FUNCTION ac_get_cv_sc(this, output) RESULT(ret)
    !! Get a scalar `COMPLEX` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` value is undefined.
    TYPE(argc), INTENT(in) :: this
      !! An argc object
    COMPLEX, INTENT(out)   :: output
      !! Ouput value
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_complex) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got COMPLEX)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      IF (.NOT.from_string(this%default,output)) &
        ret = error("Cannot cast value",-10)
    ELSE
      IF (.NOT.from_string(words_get(this%values,1), output)) &
        ret = error("Cannot cast value",-10)
    ENDIF
    RETURN
  END FUNCTION ac_get_cv_sc

  FUNCTION ac_get_sv_sc(this, output) RESULT(ret)
    !! Get a scalar `STRING` value from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                     :: this
      !! An argc object
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: output
      !! Output value
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_string) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got STRING)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      output = this%default
    ELSE
      output = words_get(this%values,1)
    ENDIF
  END FUNCTION ac_get_sv_sc

  FUNCTION ac_get_dv_ve(this, output) RESULT(ret)
    !! Get a vector of `REAL(kind=8)` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                               :: this
      !! Argc object
    REAL(kind=8), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: tmp
    LOGICAL                                           :: ok
    ret = noerror
    IF (this%ptype /= ap_real) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got REAL(kind=8))",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      ok = from_string(this%default,output(1))
      output(1:SIZE(output)) = output(1)
    ELSE
      IF (ALLOCATED(output)) DEALLOCATE(output)
      ALLOCATE(output(words_length(this%values)))
      ok = words_to_vector(this%values,tmp)
      ok = from_string(tmp, output)
    ENDIF
  END FUNCTION ac_get_dv_ve

  FUNCTION ac_get_rv_ve(this, output) RESULT(ret)
    !! Get a vector of `REAL(kind=4)` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                               :: this
      !! Argc object
    REAL(kind=4), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: tmp
    LOGICAL                                           :: ok
    ret = noerror
    IF (this%ptype /= ap_real) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got REAL(kind=4))",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      ok = from_string(this%default,output(1))
      output(1:SIZE(output)) = output(1)
    ELSE
      IF (ALLOCATED(output)) DEALLOCATE(output)
      ALLOCATE(output(words_length(this%values)))
      ok = words_to_vector(this%values,tmp)
      ok = from_string(tmp, output)
    ENDIF
  END FUNCTION ac_get_rv_ve

  FUNCTION ac_get_iv_ve(this, output) RESULT(ret)
    !! Get a vector of `INTEGER` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                          :: this
      !! Argc object
    INTEGER, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: tmp
    LOGICAL                                           :: ok
    ret = noerror
    IF (this%ptype /= ap_integer) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got INTEGER)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      ok = from_string(this%default,output(1))
      output(1:SIZE(output)) = output(1)
    ELSE
      IF (ALLOCATED(output)) DEALLOCATE(output)
      ALLOCATE(output(words_length(this%values)))
      ok = words_to_vector(this%values,tmp)
      ok = from_string(tmp, output)
    ENDIF
  END FUNCTION ac_get_iv_ve

  FUNCTION ac_get_lv_ve(this, output) RESULT(ret)
    !! Get a vector of `LOGICAL` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                          :: this
      !! Argc object
    LOGICAL, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: tmp
    LOGICAL                                           :: ok
    ret = noerror
    IF (this%ptype /= ap_logical) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got LOGICAL)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      ok = from_string(this%default,output(1))
      output(1:SIZE(output)) = output(1)
    ELSE
      IF (ALLOCATED(output)) DEALLOCATE(output)
      ALLOCATE(output(words_length(this%values)))
      ok = words_to_vector(this%values,tmp)
      ok = from_string(tmp, output)
    ENDIF
  END FUNCTION ac_get_lv_ve

  FUNCTION ac_get_cv_ve(this, output) RESULT(ret)
    !! Get a vector of `COMPLEX` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                          :: this
      !! Argc object
    COMPLEX, INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    CHARACTER(len=st_slen), ALLOCATABLE, DIMENSION(:) :: tmp
    LOGICAL                                           :: ok
    ret = noerror
    IF (this%ptype /= ap_complex) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got COMPLEX)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      ok = from_string(this%default,output(1))
      output(1:SIZE(output)) = output(1)
    ELSE
      IF (ALLOCATED(output)) DEALLOCATE(output)
      ALLOCATE(output(words_length(this%values)))
      ok = words_to_vector(this%values,tmp)
      ok = from_string(tmp, output)
    ENDIF
  END FUNCTION ac_get_cv_ve

  FUNCTION ac_get_sv_ve(this, output) RESULT(ret)
    !! Get a vector of `STRING` values from given argument
    !! If no error occured, the function always returns at least a value (whatever the parser's 
    !! state is) which is the default value if no specific values are set in the argument.
    !! Otherwise, `output` status is undefined.
    TYPE(argc), INTENT(in)                                         :: this
      !! Argc object
    CHARACTER(len=st_slen), INTENT(out), ALLOCATABLE, DIMENSION(:) :: output
      !! Output values
    TYPE(error) :: ret
      !! Error object with the first error encountered in the process
    ret = noerror
    IF (this%ptype /= ap_string) THEN
      ret = error("argparse: invalid type for output (expected `"// &
                  apt2str(this%ptype)//"' got STRING)",-21)
    ELSEIF (words_length(this%values) == 0) THEN
      ALLOCATE(output(MAX(this%nrec,1)))
      output(1:SIZE(output)) = TRIM(this%default)
    ELSE
      IF (.NOT.words_to_vector(this%values,output)) DEALLOCATE(output)
    ENDIF
  END FUNCTION ac_get_sv_ve

  ! miscellaneous methods
  ! ---------------------

  FUNCTION apt2str(ap_t) RESULT(str)
    !! Get the string representation of argparse types constants
    INTEGER, INTENT(in) :: ap_t
      !! One of ap_logical, ap_integer, ap_real, ap_complex or ap_string module constants.
    CHARACTER(len=:), ALLOCATABLE :: str
      !! String representation of the type.
    SELECT CASE(ap_t)
      CASE(ap_logical) ; str = 'logical'
      CASE(ap_integer) ; str = 'integer'
      CASE(ap_real)    ; str = 'real'
      CASE(ap_complex) ; str = 'complex'
      CASE(ap_string)  ; str = 'string'
      CASE DEFAULT     ; str = 'unknown'
    END SELECT
    RETURN
  END FUNCTION apt2str

  FUNCTION apa2str(ap_a) RESULT(str)
    !! Get the string representation of argparse actions constants 
    INTEGER, INTENT(in) :: ap_a
     !! One of ap_store, ap_append,cap_count or ap_help module constants
    CHARACTER(len=:), ALLOCATABLE :: str
      !! String representation of the action.
    SELECT CASE(ap_a)
      CASE(ap_store)  ; str = 'store'
      CASE(ap_append) ; str = 'append'
      CASE(ap_count)  ; str = 'count'
      CASE(ap_help)   ; str = 'help'
      CASE DEFAULT    ; str = 'unknown'
    END SELECT
    RETURN
  END FUNCTION apa2str

  FUNCTION get_progname() RESULT(name)
    !> Get the name of the program
    CHARACTER(len=:), ALLOCATABLE :: name !! The name of the program
    INTEGER :: c
    CALL GET_COMMAND_ARGUMENT(0,length=c)
    ALLOCATE(CHARACTER(len=c) :: name)
    CALL GET_COMMAND_ARGUMENT(0,value=name)
    c = INDEX(name, "/", back=.true.)
    IF (c /= 0.AND.c /= LEN_TRIM(name)) name = TRIM(name(c+1:))
  END FUNCTION get_progname

END MODULE ARGPARSE

