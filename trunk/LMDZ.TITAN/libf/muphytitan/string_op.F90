! Copyright Jérémie Burgalat (2013-2015,2017)
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

!! file: strings.F90
!! summary: Strings manipulation source file
!! author: J. Burgalat
!! date: 2013-2015,2017

#include "defined.h"

MODULE STRING_OP
  !! Fortran strings manipulation module
  !! 
  !! This module provides methods and objects to manipulate Fortran (allocatable) strings. It defines
  !! a doubly linked-list of strings, [[string_op(module):words(type)]] and several methods to format
  !! strings or convert them in other intrinsic types. 
  USE ERRORS
  IMPLICIT NONE
  
  PRIVATE



  PUBLIC ::  str2dble_sc,str2dble_ve,str2real_sc,str2real_ve
  ! errors module (not used but propagated)
  PUBLIC :: stdout,stderr,noerror,error, error_to_string,aborting
  
  ! misc module methods
  PUBLIC :: to_string, from_string, string_is, remove_quotes, format_string,     &
            format_paragraph, strip_newline, tokenize, str_length, endswith, &
            startswith, to_lower, to_upper, add_csi,      &
            del_csi, reset_csi, str_remove, str_replace,&
            fancy

  ! words object related methods
  PUBLIC :: words_length, words_insert, words_append, words_prepend, words_get, &
            words_set, words_get_max_width, words_get_total_width, words_pop,   &
            words_remove, words_next, words_previous, words_reset,              &
            words_valid, words_current, words_extend, words_reverse,            &
            words_reversed, words_dump, words_to_string, words_to_vector,       &
            words_clear

  PRIVATE :: fis_affect_int, fis_affect_bool, fis_affect_real,        &
             fis_affect_double, fis_affect_cplx, fis_affect_dcplx,    &
             fis_cat_int, fis_cat_bool, fis_cat_real, fis_cat_double, &
             fis_cat_cplx, fis_cat_dcplx, fis_cat_int_inv,            &
             fis_cat_bool_inv, fis_cat_real_inv, fis_cat_double_inv,  &
             fis_cat_cplx_inv, fis_cat_dcplx_inv

  ! Operators
  PUBLIC :: ASSIGNMENT(=), OPERATOR(/=), OPERATOR(==)

  INTEGER, PUBLIC, PARAMETER :: st_string  = 1 !! String type ID
  INTEGER, PUBLIC, PARAMETER :: st_logical = 2 !! Logical type ID
  INTEGER, PUBLIC, PARAMETER :: st_complex = 3 !! Complex type ID
  INTEGER, PUBLIC, PARAMETER :: st_integer = 4 !! Integer type ID
  INTEGER, PUBLIC, PARAMETER :: st_real    = 5 !! Real type ID
  
  !> List of types names
  CHARACTER(len=*), DIMENSION(5), PARAMETER, PUBLIC :: st_type_names = &
  (/ 'string ', 'logical', 'complex', 'integer', 'real   '/)

  INTEGER, PUBLIC, PARAMETER :: st_slen = SSLEN !! Maximum short string length
  INTEGER, PUBLIC, PARAMETER :: st_llen = SLLEN !! Maximum long string length


  
  INTEGER, PUBLIC, PARAMETER :: FC_BLACK     = 30 !! Black foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_RED       = 31 !! Red foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_GREEN     = 32 !! Green foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_YELLOW    = 33 !! Yellow foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_BLUE      = 34 !! Blue foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_MAGENTA   = 35 !! Magenta foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_CYAN      = 36 !! Cyan foreground csi code
  INTEGER, PUBLIC, PARAMETER :: FC_WHITE     = 37 !! White foreground csi code
  INTEGER, PUBLIC, PARAMETER :: BG_BLACK     = 40 !! Black foreground csi code
  INTEGER, PUBLIC, PARAMETER :: BG_RED       = 41 !! Black background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_GREEN     = 42 !! Green background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_YELLOW    = 43 !! Yellow background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_BLUE      = 44 !! Blue background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_MAGENTA   = 45 !! Magenta background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_CYAN      = 46 !! Cyan background csi code
  INTEGER, PUBLIC, PARAMETER :: BG_WHITE     = 47 !! White background csi code
  INTEGER, PUBLIC, PARAMETER :: ST_NORMAL    =  0 !! Normal (regular) attribute
  INTEGER, PUBLIC, PARAMETER :: ST_BOLD      =  1 !! Bold (brighter) attribute
  INTEGER, PUBLIC, PARAMETER :: ST_ITALIC    =  3 !! Italic attribute (sometimes reverse video or underline)
  INTEGER, PUBLIC, PARAMETER :: ST_UNDERLINE =  4 !! Underline attribute
  INTEGER, PUBLIC, PARAMETER :: ST_BLINK     =  5 !! Slow blink mode
  !> List of all attributes in a vector
  INTEGER, PUBLIC, PARAMETER, DIMENSION(21) :: attributes = [FC_BLACK,     &
                                                             FC_RED,       &
                                                             FC_GREEN,     & 
                                                             FC_YELLOW,    &
                                                             FC_BLUE,      &
                                                             FC_MAGENTA,   &
                                                             FC_CYAN,      &
                                                             FC_WHITE,     &
                                                             BG_BLACK,     &
                                                             BG_RED,       &
                                                             BG_GREEN,     &
                                                             BG_YELLOW,    &
                                                             BG_BLUE,      &
                                                             BG_MAGENTA,   &
                                                             BG_CYAN,      &
                                                             BG_WHITE,     & 
                                                             ST_NORMAL,    &
                                                             ST_BOLD,      &
                                                             ST_ITALIC,    &
                                                             ST_UNDERLINE, &
                                                             ST_BLINK      &
                                                            ]     
  
  !> Aliases for CSI codes.
  CHARACTER(len=2), DIMENSION(21), PARAMETER, PUBLIC :: csis =(/ &
                "fk", "fr", "fg", "fy", "fb", "fm", "fc", "fw", &
                "bk", "br", "bg", "by", "bb", "bm", "bc", "bw", &
                 "sn", "sb", "si", "su", "sk"/)

  !> [[words(type)]] object assignement interface
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE ws_affect
  END INTERFACE
  
  !> Clear either a scalar or a vector of list of [[words(type)]]
  !!
  !! The interface encapsulates words _destructors_, that deallocate memory used 
  !! by the given list(s) of words. This method should be called anytime words 
  !! object(s) is no longer used to avoid memory leaks.
  !! @note
  !! If the library support Derived type finalization, calling destructor is not
  !! mandatory.
  INTERFACE words_clear
    MODULE PROCEDURE ws_clear_sc, ws_clear_ve
  END INTERFACE

  !> Extend a given [[words(type)]] object either by another or by a string
  !!
  !! The interface encapsulates two subroutines:
  !!
  !! - [[ws_extend_ws(subroutine)]](this,other) which extends __this__ by __other__ 
  !!   (both are words objects).
  !! - [[ws_extend_str(subroutine)]](this,str,delimiter,merge) which splits __str__
  !!   according to __delimiter__ (and optionally __merge__) and then extends 
  !!   __this__ with the resulting tokens.
  INTERFACE words_extend
    MODULE PROCEDURE ws_extend_ws,ws_extend_str
  END INTERFACE

  !> Convert an intrinsic type value to a string
  !!
  !! This (very) generic interface provides conversion functions from
  !! intrinsic types to ALLOCATED string.
  !!
  !! ```
  !! (1)  FUNCTION to_string(value)     RESULT(str)
  !! (2)  FUNCTION to_string(value,fmt) RESULT(str)
  !! ```
  !! Where :
  !!
  !! - __value__ is the value to convert
  !! - __fmt__ is a string the format descriptor of the output string. Surrounding
  !!   parenthesis can be omitted.
  !! - __str__ is an allocatable string with the converted value in output, or an empty
  !!   string if the conversion failed. 
  INTERFACE to_string
    MODULE PROCEDURE int2str_as,int2str_fs
    MODULE PROCEDURE log2str_as,log2str_fs
    MODULE PROCEDURE real2str_as,real2str_fs
    MODULE PROCEDURE dble2str_as,dble2str_fs
    MODULE PROCEDURE cplx2str_as,cplx2str_fs
    MODULE PROCEDURE dcplx2str_as,dcplx2str_fs
  END INTERFACE
  
  !> Convert a string into an intrisinc type
  !!
  !! All methods defined in the interface are functions which take in arguments,
  !! a string (input) and an output variable with the relevant type (or vectors of both). 
  !! They always return an error object which is set to -5 error code (i.e. cannot cast value)
  !! on error, otherwise [[errors(module):noerror(variable)]].
  INTERFACE from_string
    MODULE PROCEDURE str2int_sc,str2log_sc,str2real_sc,str2dble_sc,str2cplx_sc
    MODULE PROCEDURE str2int_ve,str2log_ve,str2real_ve,str2dble_ve,str2cplx_ve
  END INTERFACE

  !> Compute a fancy string
  !!
  !! The generic interface adds CSI codes to the given value and returns a fortran intrinsic string.
  !!
  !! This is convinient wrapper to [[string_op(module):to_string(interface)]] and 
  !! [[string_op(module):add_csi(function)]].
  !!
  !! ```fortran
  !! FUNCTION fancy(value, flags) RESULT(str)
  !! ```
  !!
  !! **value** can be either a INTEGER,  REAL, REAL(kind=8), COMPLEX, COMPLEX(kind=8) or STRING.
  !!
  !! **flags is a vector of (string) attributes that can take the values as defined in the second
  !! column of the [list of csi attributes](|url|/page/swift/p02_strings.html#fancy-strings).
  INTERFACE fancy
    MODULE PROCEDURE fancy_fstr, fancy_int, fancy_real, fancy_double, fancy_cplx, fancy_dcplx
  END INTERFACE

  !> Overloaded string assignment operator interface
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE fis_affect_int, fis_affect_bool, fis_affect_real,    &
                     fis_affect_double, fis_affect_cplx, fis_affect_dcplx
  END INTERFACE

  !> Overloaded string concatentation operator interface
  INTERFACE OPERATOR(//)
    MODULE PROCEDURE fis_cat_int, fis_cat_bool, fis_cat_real, fis_cat_double, &
                     fis_cat_cplx, fis_cat_dcplx
    MODULE PROCEDURE fis_cat_int_inv, fis_cat_bool_inv, fis_cat_real_inv,     &
                     fis_cat_double_inv, fis_cat_cplx_inv, fis_cat_dcplx_inv
  END INTERFACE

  !> Define a linked word
  !!
  !! Linked words are only intended to be used within a words type.
  !! It's part of the doubly linked list words.
  TYPE, PUBLIC :: word
#if HAVE_FTNDTSTR    
    CHARACTER(len=:), ALLOCATABLE :: value !! Value of the word
#else    
    !> Value of the word
    !!
    !! @warning 
    !! It is always limited to [[string_op(module):st_slen(variable)]] characters.
    CHARACTER(len=st_slen)        :: value = ''
#endif
    TYPE(word), PRIVATE, POINTER  :: next => null() !! Next word in the list of words
    TYPE(word), PRIVATE, POINTER  :: prev => null() !! Previous word in the list of words
  END TYPE word
  
  !> Define a list of words
  TYPE, PUBLIC :: words
    INTEGER, PRIVATE             :: nw = 0         !! Number of word in the list
    TYPE(word), PRIVATE, POINTER :: head => null() !! First word in the list
    TYPE(word), PRIVATE, POINTER :: tail => null() !! Last word in the list
    TYPE(word), PRIVATE, POINTER :: iter => null() !! Current word (iterator)
#if HAVE_FTNPROC 
    CONTAINS
    PROCEDURE, PRIVATE :: ws_extend_ws
    PROCEDURE, PRIVATE :: ws_extend_str
    PROCEDURE, PUBLIC :: length      => words_length
      !! Get the number of words in the list
    PROCEDURE, PUBLIC :: insert      => words_insert
      !! Insert a word at given index
    PROCEDURE, PUBLIC :: append      => words_append
      !! Append a word at the end of the list 
    PROCEDURE, PUBLIC :: prepend     => words_prepend
      !! Prepend a word at the beginning of the list 
    PROCEDURE, PUBLIC :: get         => words_get
      !! Get the word at given index
    PROCEDURE, PUBLIC :: set         => words_set
      !! Set a word at given index
    PROCEDURE, PUBLIC :: max_width   => words_get_max_width
      !! Get the width of the biggest word in the list
    PROCEDURE, PUBLIC :: total_width => words_get_total_width
      !! Get the total width of the words stored in the list
    PROCEDURE, PUBLIC :: reverse     => words_reverse
      !! Reverse the list in place
    PROCEDURE, PUBLIC :: reversed    => words_reversed
      !! Get a reversed copy of the list 
    PROCEDURE, PUBLIC :: dump        => words_dump
      !! Dump words of the list (on per line)
    PROCEDURE, PUBLIC :: tostring    => words_to_string
      !! Convert the list in a single string
    PROCEDURE, PUBLIC :: to_vector   => words_to_vector
      !! Convert the list in a vector
    PROCEDURE, PUBLIC :: pop         => words_pop
      !! Pop a word from the list and returns it 
    PROCEDURE, PUBLIC :: remove      => words_remove
      !! Remove a word from the list
    PROCEDURE, PUBLIC :: next        => words_next
      !! Go to the next word in the list
    PROCEDURE, PUBLIC :: previous    => words_previous
      !! Go to the previous word in the list
    PROCEDURE, PUBLIC :: reset       => words_reset
      !! Reset the list's iterator
    PROCEDURE, PUBLIC :: valid       => words_valid
      !! Check if iterator position is valid
    PROCEDURE, PUBLIC :: current     => words_current
      !! Get the current word in the list
    GENERIC, PUBLIC :: extend => ws_extend_ws,ws_extend_str
      !! Extend a list with either a string or another list of words
#endif
  END TYPE words
  
  CONTAINS
  
  FUNCTION word_length(this) RESULT(lgth)
    !! Get the trimmed length of the word object
    TYPE(word), INTENT(in) :: this
      !! A word object
    INTEGER :: lgth
      !! The length of the word's value (without trailing spaces)
#if HAVE_FTNDTSTR    
    IF (.NOT.ALLOCATED(this%value)) THEN
      lgth = 0 ; RETURN
    ENDIF
#endif    
    lgth = LEN_TRIM(this%value)
    RETURN
  END FUNCTION word_length
  
  SUBROUTINE disconnect_word(this)
    !! Disconnect a word object
    !!
    !! The object is no more connected to its neighbours which are connected together.
    !! @note 
    !! After this method is called the object is no longer connected to its parent words 
    !! object and should be deallocated in order to avoid memory leaks.
    TYPE(word), INTENT(inout) :: this
      !! A word object to disconnect
    TYPE(word), POINTER :: pw,nw
    nw => this%next ; pw => this%prev
    IF (ASSOCIATED(nw)) nw%prev => pw
    IF (ASSOCIATED(pw)) pw%next => nw
    RETURN
  END SUBROUTINE disconnect_word
  
  SUBROUTINE ws_affect(this,other)
    !! words object assignment operator subroutine
    TYPE(words), INTENT(out) :: this
      !! A words object to be assigned
    TYPE(words), INTENT(in)  :: other
      !! A words object to assign
    TYPE(word), POINTER :: cur
    CALL ws_clear_sc(this)
    IF (other%nw == 0) THEN
      RETURN
    ELSE
      cur => other%head
      DO WHILE(associated(cur))
#if HAVE_FTNDTSTR      
        IF (.NOT.ALLOCATED(cur%value)) THEN
          CALL words_append(this,"")
        ELSE
          CALL words_append(this,cur%value)
        ENDIF
#else
        CALL words_append(this,cur%value)
#endif
        IF (ASSOCIATED(cur,other%iter)) this%iter => this%tail
        cur => cur%next
      ENDDO
    ENDIF
    RETURN
  END SUBROUTINE ws_affect
  
  SUBROUTINE ini_word(this,value)
    !! Initialize the first word of a list of words
    !!
    !! This subroutine is not a constructor. It is only intended to set the first word 
    !! object in a words object.
    TYPE(words), INTENT(inout)   :: this
      !! A words object reference
    CHARACTER(len=*), INTENT(in) :: value
      !! A string with the word used to initialize the list 
    ALLOCATE(this%head)
    this%tail => this%head
    ASSIGN_DTSTR(value,this%tail%value)
    this%nw = 1
    RETURN
  END SUBROUTINE ini_word

  SUBROUTINE ws_clear_sc(obj)
    !! Clear a list of words
    !!
    !! This subroutine deallocates all memory used by the given words object.
    !! @warning 
    !! The subroutine should be called whenever a words is no more used (e.g. at 
    !! the end of the current scope), otherwise memory leaks could occur.
    TYPE(words),INTENT(inout), TARGET :: obj
      !! A words object to clear
    TYPE(word), POINTER :: cur,next
    IF (obj%nw == 0) RETURN
    cur => obj%head 
    DO WHILE(ASSOCIATED(cur))
      next => cur%next
      CALL disconnect_word(cur)
#if HAVE_FTNDTSTR
      IF (ALLOCATED(cur%value)) DEALLOCATE(cur%value)
#endif
      DEALLOCATE(cur)
      cur => next
    ENDDO
    obj%nw = 0
    obj%head => null() ; obj%tail => null() ; obj%iter => null()
  END SUBROUTINE ws_clear_sc

  SUBROUTINE ws_clear_ve(objs)
    !! Clear a vector of list of words
    !!
    !! This subroutine deallocates all memory used by the given vector of words objects.
    !! @warning 
    !! The subroutine should be called whenever a words is no more used (e.g. at the end 
    !! of the current scope), otherwise memory leaks could occur.
    TYPE(words),INTENT(inout), DIMENSION(:) :: objs
      !! A vector of words objects to clear
    TYPE(word), POINTER :: cur,next
    INTEGER             :: i
    DO i=1,SIZE(objs)
      call ws_clear_sc(objs(i))
    ENDDO
  END SUBROUTINE ws_clear_ve

  SUBROUTINE ws_extend_ws(this, other)
    !! Extend a list of words with another one
    OBJECT(words), INTENT(inout) :: this
      !! A words object to extend
    TYPE(words), INTENT(in)     :: other
      !! A words object to extend with
    TYPE(word), POINTER :: cw
    IF (other%nw == 0) RETURN
    cw => other%head
    DO WHILE(ASSOCIATED(cw))
      CALL words_append(this,cw%value) ; cw => cw%next
    ENDDO
    RETURN
  END SUBROUTINE ws_extend_ws

  SUBROUTINE ws_extend_str(this,str,delimiter,merge,protect) 
    !> Extend a list of word with a given string
    !! @details The method adds a new list of words to the current list by 
    !! splitting a string using a set of delimiters.
    !! 
    !!   - If __delimiter__ is not given, THEN blank space is used.
    !!   - __delimiter__ can be a string of any length, but each character of 
    !!     the sequence is seen as a single delimiter. Each time one of these 
    !!     special character is seen on the string, it is splitted.
    !!   - If __protect__ is set to .true. THEN delimiter enclosed by
    !!     either single or double quotes are protected.
    !!   - The optional argument __merge__ instructs the method wether to merge
    !!     or not successive delimiters in the string.
    !! 
    !! For example, considering the following string:
    !! <center>@verbatim "I like coffee and bananas." @endverbatim</center>
    !!   - Used with only __delimiter__ = " e", the method returns the list:
    !!     <center>"I","lik","","coff","","","and","bananas"</center>
    !!   - Used with both __delimiter__ = " e" and __merge__ = .true. :
    !!     <center>"I","lik","coff","and","bananas"</center>
    !! @warning
    !! The method does not trim or adjust the input string. Consequently, it can
    !! add several empty words at the end of the list if the string is not well
    !! defined.
    !! @warning To avoid such problems, consider using TRIM() and ADJUSTL()
    !! function on __str__ actual argument when calling this subroutine.
    OBJECT(words), INTENT(inout), TARGET   :: this
      !! A words object to extend
    CHARACTER(len=*), INTENT(in)           :: str
      !! A string to split in words
    CHARACTER(len=*), INTENT(in), OPTIONAL :: delimiter
      !! An optional string with the words delimiters (default to blank space). 
    LOGICAL, INTENT(in), OPTIONAL          :: merge
      !! An optional boolean control flag that instructs the method
      !! wether to merge or not successive delimiters (default to .false.)
    LOGICAL, INTENT(in), OPTIONAL          :: protect
      !! An optional boolean flag with .true. to indicate that 
      !! delimiter characters between quotes are protected
    ! - LOCAL
    INTEGER                       :: sl,p,i,j,stat
    LOGICAL                       :: zmerge,zprotect,indq,insq,outer
    CHARACTER(len=:), ALLOCATABLE :: seps
    CHARACTER(len=:), ALLOCATABLE :: curw
    CHARACTER(len=1), PARAMETER   :: sq = CHAR(39) ! single quote ascii code
    CHARACTER(len=1), PARAMETER   :: dq = CHAR(34) ! double quotes ascii code
    stat=0 ; p=1 ; indq = .false. ; insq = .false. 
    seps = ' '
    zmerge = .false. ; IF (PRESENT(merge)) zmerge = merge
    zprotect = .true. ; IF (PRESENT(protect)) zprotect = protect
    IF (PRESENT(delimiter)) THEN
      IF (LEN(delimiter) > 0) seps = delimiter
    ENDIF
    sl = LEN(str) ; IF (sl == 0) RETURN
    outer =     (INDEX(str,sq) == 1 .AND. INDEX(str,sq,.true.) == LEN(str)) & 
            .OR.(INDEX(str,dq) == 1 .AND. INDEX(str,dq,.true.) == LEN(str))
    ! no delimiter found or (have outer quotes and should protect)
    IF (SCAN(str,seps) == 0.OR.(outer.AND.zprotect)) THEN
      CALL words_append(this,remove_quotes(str)) 
      RETURN
    ENDIF
    ! We have to loop... 
    i = 1 ; curw=''
    DO
      IF (i > sl) EXIT
      p = SCAN(str(i:),seps) ! position du delimiteur
      IF (p == 0) THEN
        ! a gerer
        curw = curw//TRIM(str(i:))
        CALL words_append(this,TRIM(str(i:))) ; EXIT
        curw=''
      ELSE
        IF (zprotect) THEN
          j=i
          ! starting state
          DO WHILE(j<i+p)
            IF (str(j:j) == sq.AND. .NOT.indq) insq = .NOT.insq
            IF (str(j:j) == dq.AND. .NOT.insq) indq = .NOT.indq
            j = j+1
          ENDDO
          IF ((insq.AND.INDEX(str(j:),"'")/=0) .OR. &
              (indq.AND.INDEX(str(j:),'"')/=0)) THEN
            curw=curw//str(i:i+p-1)
            i=i+p ; CYCLE
          ENDIF
        ENDIF
        IF (p == 1) THEN
          IF (.NOT.zmerge) THEN
            curw=''
            CALL words_append(this,curw)
          ENDIF
          i = i + 1 ; CYCLE
        ELSE
          curw=curw//str(i:i+p-2)
          CALL words_append(this,curw)
          curw = ''
          i = i + p
        ENDIF
      ENDIF
    ENDDO
    DEALLOCATE(curw) ; DEALLOCATE(seps)
    !IF (zprotect) THEN
    !  ! catching unbalanced quotes
    !  IF (insq .OR. indq) &
    !  WRITE(*,'(a)') "extends:warning: unbalanced quotes"
    !ENDIF
    RETURN
  END SUBROUTINE ws_extend_str

  FUNCTION ws_get_ptr(this,idx) RESULT(pted)
    !! Get the pointer of the word object at given index
    !!
    !! The method returns the pointer of the word object at the given index. 
    !! If index is out of range a null poitner is returned.
    OBJECT(words), INTENT(in) :: this
      !! A words object 
    INTEGER, INTENT(in)       :: idx
      !! An integer with the index of the desired object in __this__
    TYPE(word), POINTER :: pted
      !! A pointer to the selected word object.
    INTEGER :: i
    pted => null()
    IF (idx < 1 .OR. idx > words_length(this)) THEN
      RETURN
    ENDIF
    IF (idx > (this%nw+1)/2) THEN
      pted => this%tail
      DO i=1,this%nw - idx ; pted => pted%prev ; ENDDO
    ELSE
      pted => this%head
      DO i=1,idx-1 ; pted => pted%next ; ENDDO
    ENDIF
    RETURN
  END FUNCTION ws_get_ptr

  FUNCTION words_length(this) RESULT(res)
    !! Get the size of the words object.
    !!
    !! The method returns the number of words stored in the given list of words.
    OBJECT(words), INTENT(in) :: this !! A words object.
    INTEGER :: res                    !! The number of words in the object.
    res = this%nw
    RETURN
  END FUNCTION words_length

  SUBROUTINE words_insert(this, idx, value)
    !! Insert a word before given index in a list of words.
    !!
    !! The method inserts a new word before the given index in the list of words. If the given index is out
    !! of range, the method prepends/appends the object based on the index value.
    OBJECT(words), INTENT(inout)  :: this
      !! A words object.
    INTEGER, INTENT(in)          :: idx
      !! An integer with the index of an object in the list. The new object will be inserted before that index.
    CHARACTER(len=*), INTENT(in) :: value
      !! A string with the word to insert in the list.
    TYPE(word), POINTER :: welt,nx,pv
    INTEGER             :: i
    welt => null() ; nx => null() ; pv => null()
    IF (this%nw == 0) THEN
      CALL ini_word(this,value)
    ELSE IF (idx > this%nw) THEN
      this%nw = this%nw + 1
      welt => this%tail
      allocate(this%tail)
      ASSIGN_DTSTR(value,this%tail%value)
      this%tail%prev => welt
      this%tail%prev%next => this%tail
    ELSE IF (idx <= 1) THEN
      this%nw = this%nw + 1
      welt => this%head
      allocate(this%head)
      ASSIGN_DTSTR(value,this%head%value)
      this%head%next => welt
      this%head%next%prev => this%head
    ELSE
      IF (idx > (this%nw+1)/2) THEN
        nx => this%tail 
        DO i=1, this%nw - idx ; nx => nx%prev ; ENDDO
      ELSE
        nx => this%head 
        DO i=1, idx-1 ; nx => nx%next ; ENDDO
      ENDIF
      pv => nx%prev
      allocate(welt)
      ASSIGN_DTSTR(value,welt%value)
      welt%prev => pv ; welt%next => nx
      pv%next => welt ; nx%prev => welt
      this%nw = this%nw + 1
    ENDIF
    RETURN
  END SUBROUTINE words_insert

  SUBROUTINE words_append(this,value)
    !! Append a word to the list of word
    !! 
    !! The method appends a word to the list of word. This is a convinient wrapper to 
    !! [[string_op(module)::words_insert(subroutine)]] to add a new word at the beginning of the list.
    OBJECT(words), INTENT(inout) :: this  !! A words object 
    CHARACTER(len=*), INTENT(in) :: value !! A string to append
    !CALL words_insert(this,this%nw+1,value)
    type(word), pointer :: np

    call words_insert(this,this%nw+1,value)

    !! If the list is empty
    !if (this%nw == 0) then
    !  call ini_word(this, value)
    !  return
    !end if
    !! Add new element to the end
    !this%nw = this%nw + 1
    !np => this%tail
    !allocate(this%tail)
    !this%tail%value = TRIM(value)
    !this%tail%prev => np
    !this%tail%prev%next => this%tail 
    RETURN
  END SUBROUTINE words_append

  SUBROUTINE words_prepend(this,value)
    !! Prepend a word to the list of word
    !!
    !! The method prepends a word to the list of word. This is a convinient wrapper to
    !! [[string_op(module)::words_insert(subroutine)]] to add a new word at the end of the list.
    OBJECT(words), INTENT(inout) :: this  !! A words object 
    CHARACTER(len=*), INTENT(in) :: value !! A string to prepend
    CALL words_insert(this,0,value)
    RETURN
  END SUBROUTINE words_prepend

  FUNCTION words_get(this,idx,case) RESULT (res)
    !! Get the word's value at given index
    !! 
    !! The method attempts to get the word's value at the given index. If index is out of range
    !! an empty string is returned.
    !! @note 
    !! The returned string is always trimmed.
    OBJECT(words), INTENT(in)              :: this
      !! A words object reference
    INTEGER, INTENT(in)                    :: idx
      !! An integer with the index of a word in the list
    CHARACTER(len=5), INTENT(in), OPTIONAL :: case 
      !! An optional string with either 'upper' or 'lower' to get the value converted in the relevant case
    CHARACTER(len=:), ALLOCATABLE :: res
      !! The value of the word stored at given index in the list of words
    TYPE(word), POINTER :: cur
    cur => ws_get_ptr(this,idx)
    IF (.not.associated(cur)) THEN
      res = '' ; RETURN
    ENDIF
    IF (PRESENT(case)) THEN
      IF (case == "upper") res = to_upper(cur%value) 
      IF (case == "lower") res = to_lower(cur%value)
    ELSE
      res = TRIM(cur%value)
    ENDIF
    RETURN
  END FUNCTION words_get

  SUBROUTINE words_set(this,idx,value)
    !! Set a new value to a word object in the list of words at given index
    !!
    !! The method sets a new word at given index. If index is out of range, the method simply does nothing.
    OBJECT(words), INTENT(inout) :: this  !! A words object
    INTEGER, INTENT(in)          :: idx   !! An integer with the index of the word object to modify in the list
    CHARACTER(len=*), INTENT(in) :: value !! A string with the new value to set
    TYPE(word), POINTER :: cur
    cur => ws_get_ptr(this,idx)
    IF (.NOT.ASSOCIATED(cur)) RETURN
    cur%value = value
  END SUBROUTINE words_set

  FUNCTION words_get_max_width(this) RESULT(res)
    !! Get the longest word's width in the words object
    !!
    !! The method computes and returns the longest (trimmed) word's width in the words object.
    OBJECT(words), INTENT(in) :: this !! A words object 
    INTEGER :: res                    !! An integer with the maximum width (0 if the list is empty)
    TYPE(word), POINTER :: cur
    res = 0
    IF (this%nw == 0) RETURN
    cur => this%head ; res = word_length(cur)
    DO WHILE(ASSOCIATED(cur%next))
      cur => cur%next
      IF (word_length(cur) > res) res = word_length(cur)
    ENDDO
    RETURN
  END FUNCTION words_get_max_width

  FUNCTION words_get_total_width(this) RESULT(width)
    !! Get the total width of all words stored in the list of words
    !! 
    !! The method computes and returns the total width of all words stored in 
    !! the list of words.
    !! @note 
    !! Total width is computed using strings::word_length so it only takes
    !! into account trimmed words (without trailing blanks)
    !! @note 
    !! If csi codes have been added to words elements they are counted in the width.
    OBJECT(words), INTENT(in) :: this !! A words object
    INTEGER :: width                  !! Total length of the list of words
    TYPE(word), POINTER :: cur
    width = 0
    IF (this%nw == 0) RETURN
    cur => this%head ; width = word_length(cur)
    DO WHILE(ASSOCIATED(cur%next))
      cur => cur%next
      width = width + word_length(cur)
    ENDDO
    cur => null()
    RETURN
  END FUNCTION words_get_total_width

  SUBROUTINE words_reverse(this)
    !! Reverse the list of words in-place
    OBJECT(words), INTENT(inout) :: this
      !! A words object to reverse
    TYPE(word), POINTER :: loop,iwc,iwp
    IF (this%nw <= 1) RETURN
    loop => this%head ; iwc=> this%head ; iwp=> null()
    DO WHILE(ASSOCIATED(loop%next))
      loop => loop%next
      iwp => iwc%prev ; iwc%prev => iwc%next ; iwc%next => iwp
      iwc => loop
    ENDDO
    iwp=>this%tail%prev ; this%tail%prev=>this%tail%next ; this%tail%next=>iwp
    iwc => this%head ; this%head => this%tail ; this%tail => iwc
    loop => null() ; iwc => null() ; iwp => null()
    RETURN
  END SUBROUTINE words_reverse

  FUNCTION words_reversed(this) RESULT(res)
    !! Get a reversed copy of the list of words
    OBJECT(words), INTENT(in) :: this
      !! A words object to reverse
    TYPE(words) :: res
      !! A reversed copy of the input list of words
    TYPE(word),POINTER  :: cur
    IF(this%nw == 0) RETURN
    cur => this%tail
    DO WHILE(ASSOCIATED(cur))
      CALL words_append(res,cur%value)
      IF (ASSOCIATED(cur,this%iter)) res%iter => res%tail 
      cur => cur%prev
    ENDDO
    cur => null()
    RETURN
  END FUNCTION words_reversed

  SUBROUTINE words_dump(this,lun)
    !! Dump the list of words
    !! 
    !! The method dumps on the given logical unit the elements of the list one by line.
    OBJECT(words), INTENT(in)     :: this
      !! A words object to dump
    INTEGER, INTENT(in), OPTIONAL :: lun
      !! An optional integer with the printing logical unit. If not given, the list is dumped on 
      !! standard output stream.
    TYPE(word), POINTER :: cur
    INTEGER             :: lu
    IF (this%nw == 0) RETURN
    lu=6 ; IF (PRESENT(lun)) lu = lun
    cur => this%head
    DO WHILE(ASSOCIATED(cur))
      WRITE(lu,'(a)') TRIM(cur%value)
      cur => cur%next
    ENDDO
    cur => null()
    RETURN
  END SUBROUTINE words_dump

  FUNCTION words_to_string(this, delimiter,protect) RESULT(str)
    !! Convert the list of words into a string
    !!
    !! The method converts the list of words into a string. In output, string is always
    !! allocated even if the list is empty.
    !!
    !! If `protect` is set to .true. (default to .false.), each word is enclose between
    !! double quotes. This option can be used to perform operation on the string before
    !! setting it back as values of the list of words.
    OBJECT(words), INTENT(in)              :: this
      !! A words object 
    CHARACTER(len=*), INTENT(in), OPTIONAL :: delimiter
      !! An optional string used as delimiter between each words
    LOGICAL, INTENT(in), OPTIONAL :: protect
      !! Optional control flag with .true. to protect each word during concatentation.
    CHARACTER(len=:), ALLOCATABLE :: str
      !! An allocatable string with the list of words joined by the given delimiter (if any)
    TYPE(word), POINTER :: cur
    LOGICAL :: zprotect
    zprotect = .false. ; IF (PRESENT(protect)) zprotect = protect
    str = ''
    IF (this%nw == 0) RETURN
    cur => this%head
    DO WHILE(ASSOCIATED(cur))
      IF (zprotect) THEN
        str=str//'"'//TRIM(cur%value)//'"'
      ELSE
        str=str//TRIM(cur%value)
      ENDIF
      IF (PRESENT(delimiter).AND..NOT.ASSOCIATED(cur,this%tail)) &
        str=str//delimiter
      cur => cur%next
    ENDDO
    RETURN
  END FUNCTION words_to_string

  FUNCTION words_to_vector(this,ret) RESULT(ok)
    !! Convert the list of words into a vector of strings
    !!
    !! The method attempts to convert the list of words in a vector of strings.
    !! If _this_ list of words is empty, the output vector is allocated with 0 elements and the method returns
    !! .false., otherwise it returns .true.
    !! @note
    !! If elements in __this__ words object are wider than [[string_op(module):st_slen(variable)]], output 
    !! values will be truncated.
    OBJECT(words), INTENT(in)                                      :: this
      !! A words object reference
    CHARACTER(len=st_slen), INTENT(out), ALLOCATABLE, DIMENSION(:) :: ret
      !! An allocatable vector of assumed length string with the words of __this__
    LOGICAL             :: ok
      !! Return status.
    INTEGER             :: l,mw
    TYPE(word), POINTER :: iw
    ok = .true.
    l = words_length(this)
    IF (l == 0) THEN
      ALLOCATE(ret(0))
      ok = .false.
      RETURN
    ENDIF
    ALLOCATE(ret(l)) ; mw = LEN(ret(l))
    ret(1:l) = ' ' ! really needed ?
    iw => this%head ; l=1
    DO WHILE(ASSOCIATED(iw))
       ret(l) = TRIM(iw%value) ; l=l+1 ; iw => iw%next
    ENDDO
  END FUNCTION words_to_vector

  FUNCTION words_pop(this,idx,move_forward) RESULT(value)
    !! Pop a word in the list of words
    !!
    !! The method removes the word of the list at given index and returns it. If no index is given, 
    !! last word of the list is removed.
    !!
    !! If the index is out of range, the method does nothing and returns an empty string.
    !!
    !! By default, if the iterator is located on the item to be removed, it is moved backward before
    !! deletion occurs. If __move\_forward__ is set to .true., the iterator is moved forward.
    OBJECT(words), INTENT(inout)  :: this
     !! A words object
    INTEGER, INTENT(in), OPTIONAL :: idx
      !! Optional index of the word to delete
    LOGICAL, INTENT(in), OPTIONAL :: move_forward 
      !! Move the iterator forward if needed. By default the iterator is moved backward. 
    CHARACTER(len=:), ALLOCATABLE :: value
      !! The word's value at given index 
    LOGICAL             :: zforward
    INTEGER             :: zidx
    TYPE(word), POINTER :: cur
    zidx=words_length(this) ; IF (PRESENT(idx)) zidx = idx
    zforward = .false. ; IF (PRESENT(move_forward)) zforward = move_forward
    cur => ws_get_ptr(this,zidx)
    IF (.NOT.ASSOCIATED(cur)) THEN
      value = '' ; RETURN
    ELSE IF (ASSOCIATED(cur,this%iter)) THEN
      IF (zforward) THEN
        CALL words_next(this)
      ELSE
        CALL words_previous(this)
      ENDIF
    ENDIF
    value = TRIM(cur%value)
    CALL disconnect_word(cur)
    DEALLOCATE(cur)
    this%nw = this%nw - 1
    RETURN
  END FUNCTION words_pop

  SUBROUTINE words_remove(this,idx,move_forward)
    !! Remove the word of the list at given index
    !!
    !! The method removes the word of the list at given index. If no index is given, last word 
    !! of the list is removed.
    !!
    !! If the index is out of range, the method does nothing.
    !!
    !! By default, if the iterator is located on the item to be removed, it is moved backward before
    !! deletion occurs. If __move\_forward__ is set to .true., the iterator is moved forward.
    OBJECT(words), INTENT(inout)  :: this
      !! A words object
    INTEGER, INTENT(in), OPTIONAL :: idx
      !! Index of the word to delete
    LOGICAL, INTENT(in), OPTIONAL :: move_forward
      !! Move the iterator forward if needed. By default the iterator is moved backward. 
    LOGICAL             :: zforward 
    INTEGER             :: zidx
    TYPE(word), POINTER :: cur
    zidx=words_length(this) ; IF(PRESENT(idx)) zidx = idx
    zforward = .false. ; IF (PRESENT(move_forward)) zforward = move_forward
    cur => ws_get_ptr(this,idx)
    IF (.NOT.ASSOCIATED(cur)) THEN
      RETURN
    ELSE IF (ASSOCIATED(cur,this%iter)) THEN
      IF (zforward) THEN
        CALL words_next(this)
      ELSE
        CALL words_previous(this)
      ENDIF
    ENDIF
    CALL disconnect_word(cur)
    DEALLOCATE(cur)
    this%nw = this%nw - 1
    RETURN
  END SUBROUTINE words_remove

  SUBROUTINE words_next(this)
    !! Go to the next word in the list
    OBJECT(words), INTENT(inout) :: this !! A words object
    IF (ASSOCIATED(this%iter)) this%iter => this%iter%next
  END SUBROUTINE words_next

  SUBROUTINE words_previous(this)
    !! Go to the previous word in the list
    OBJECT(words), INTENT(inout) :: this !! A words object
    IF (ASSOCIATED(this%iter)) this%iter => this%iter%prev
  END SUBROUTINE words_previous

  FUNCTION words_valid(this) RESULT(ret)
    !! Check if the current iterated word is valid 
    OBJECT(words), INTENT(in) :: this !! A words object
    LOGICAL :: ret                    !! A logical flag with .true. if the current iterated word is valid
    ret = associated(this%iter)
  END FUNCTION words_valid

  FUNCTION words_current(this) RESULT(wrd)
    !! Get current word value
    OBJECT(words), INTENT(in) :: this
      !! A words object
    CHARACTER(len=:), ALLOCATABLE :: wrd
      !! A string with the value of the current word or __an unallocated string__ if current word 
      !! is not valid (see [[string_op(module):words_valid(function)]]).
    IF (ASSOCIATED(this%iter)) THEN
      wrd = this%iter%value
    ENDIF
  END FUNCTION words_current

  SUBROUTINE words_reset(this,to_end)
    !! Reset the iterator 
    !!
    !! The method resets the iterator either at the beginning or at the end of the list of words 
    !! (if __to_end__ is set to .true.).
    OBJECT(words), INTENT(inout)  :: this   !! A words object 
    LOGICAL, INTENT(in), OPTIONAL :: to_end !! An optional logical flag with .true. to reset the iterator at the end of the list
    this%iter => this%head
    IF (PRESENT(to_end)) THEN
      IF (to_end) this%iter => this%tail
    ENDIF
  END SUBROUTINE words_reset

  ! Fancy string methods
  ! --------------------

  FUNCTION tokenize(str,vector,delimiter,merge,protect) RESULT(ok)
    !! Tokenize a string.
    CHARACTER(len=*), INTENT(in)                             :: str
      !! A string to tokenize
    CHARACTER(len=*), INTENT(out), DIMENSION(:), ALLOCATABLE :: vector
      !! An allocatable vector of strings with the tokens found. If string cannot be tokenized, 
      !! the vector is __allocated to 0 elements__ and the method returns .false..
    CHARACTER(len=*), INTENT(in), OPTIONAL                   :: delimiter
      !! An optional string with the words delimiters. It is set to blank space by default. 
      !! Note that each character is seen as a single delimiter.
    LOGICAL, INTENT(in), OPTIONAL                            :: merge
      !! An optional boolean control flag with .true. that instructs the method whether to 
      !! merge or not successive delimiters. Default to .false.
    LOGICAL, INTENT(in), OPTIONAL                            :: protect
      !! An optional boolean flag with .true. to indicate that delimiter characters between 
      !! quotes are protected. Default to .true.
    LOGICAL :: ok
      !! Return status (.true. on success).
    CHARACTER(len=:), ALLOCATABLE :: seps
    TYPE(words)                   :: tmp
    LOGICAL                       :: zmerge,zprotect
    TYPE(word), POINTER           :: wrd
    integer                       :: i,nw
    ok = .true.
    zmerge = .false. ; zprotect = .true. ; seps = ' ' 
    IF (PRESENT(merge)) zmerge = merge
    IF (PRESENT(protect)) zprotect = protect
    IF (PRESENT(delimiter)) THEN
      IF (LEN(delimiter) > 0 ) seps = delimiter
    ENDIF
    call words_extend(tmp,str,seps,zmerge,zprotect)
    nw = tmp%nw
    i =  1
    ALLOCATE(vector(tmp%nw))
    IF (nw > 0) THEN
      wrd => tmp%head
      DO WHILE(ASSOCIATED(wrd))
        vector(i) = TRIM(wrd%value)
        wrd => wrd%next
        i = i+1
      ENDDO
    ELSE
      ok = .false.
    ENDIF
    call words_clear(tmp)
    RETURN
  END FUNCTION tokenize

  FUNCTION remove_quotes(str) RESULT(ostr)
    !! Strips outer quotes from string
    !!
    !! The function removes only external quotes from the input string
    !! and returns the result in an allocatable string.
    !! Quotes are removed only if they are the first and last non blank
    !! characters. Either double and single quotes are stripped without distinction.
    !! The output string is trimmed from leading and trailing blank spaces (after quotes removal !)
    CHARACTER(len=*), INTENT(in)  :: str  !! A string to check
    CHARACTER(len=:), ALLOCATABLE :: ostr !! A string without external quotes (if any).  
    CHARACTER(len=1), PARAMETER   :: sq=CHAR(39), dq=CHAR(34)
    CHARACTER(len=2), PARAMETER   :: dsq=CHAR(39)//CHAR(34)
    INTEGER                       :: i, j
    IF (LEN_TRIM(str) == 0) RETURN
    ostr = TRIM(ADJUSTL(str))
    i = SCAN(ostr,sq//dq) ; j = SCAN(ostr,sq//dq,.true.)
    IF (i == j) RETURN
    IF (i /= 1) i = 0
    IF (j /= LEN(ostr)) j = LEN(ostr)+1
    ostr = TRIM(ostr(i+1:j-1))
    RETURN
  END FUNCTION remove_quotes

  FUNCTION string_is(str) RESULT(ret)
    !! Check if string represents an intrinsic type
    !!
    !! The method checks if the given string represents an intrinsic type. Both logical and complex type 
    !! are checked in a strict way :
    !!
    !! - A string is a logical if it is one of the following value: __.false.__, __.true.__, __F__, __T__.
    !! - A string is potentially a complex if it has the following format: __(\*\*\*,\*\*\*)__ where 
    !!   __\*\*\*__ is checked to see wether it is numerical or not.
    !!
    !! Valid numerical values can take the following forms:
    !! ```
    !!   [0-9]
    !!   [0-9]*.?[0-9]*?([ed][+-]?[0-9]+)?
    !! ```
    !! Obviously if returned value is greater than 3, the string can be converted in 
    !! floating point value.
    !!
    !! Empty input string is simply considered to be of string type !  
    CHARACTER(len=*), INTENT(in) :: str
      !! A string to check
    INTEGER :: ret
      !! An integer with the intrinsic type related to the string.
      !!
      !! Types are one of the following parameters
      !!
      !! - [[string_op(module):st_string(variable)]] (1) for string
      !! - [[string_op(module):st_logical(variable)]] (2) for logical
      !! - [[string_op(module):st_complex(variable)]] (3) for complex
      !! - [[string_op(module):st_integer(variable)]] (4) for integer
      !! - [[string_op(module):st_real(variable)]] (5) for floating point value
    CHARACTER(len=:), ALLOCATABLE :: zs,zzs
    INTEGER :: j,l
    ret = 1 ; IF (LEN_TRIM(str) == 0) RETURN
    zs = to_lower(TRIM(ADJUSTL(str))) ; j = INDEX(zs,',') ; l = len(zs)
    IF (zs(1:1)=='('.AND.zs(l:l) == ')'.AND.j==INDEX(zs,',')) THEN
      IF (j == 2 .OR. j == l-1) RETURN
      zzs = TRIM(ADJUSTL(zs(2:j-1))) ; IF (what_(zzs) < 3) RETURN
      zzs = TRIM(ADJUSTL(zs(j+1:l-1))) ; ret = what_(zzs)
      IF (ret > 3) THEN ; ret = 3 ; ELSE ; ret = 1 ; ENDIF
    ELSE
      ret = what_(zs)
    ENDIF
    CONTAINS
      FUNCTION what_(s) RESULT(is)
        !! Check if the given string is numerical, logical or a simple string
        !! @note
        !! Input string should be in lower case, otherwise, the method will give a a wrong result.
        !! @warning
        !! The test performed for logical checking is quite strict : A string is considered as logical
        !! if and only if it is one of the following values : __.false.__, __.true.__, __F__, __T__.
        CHARACTER(len=*), INTENT(in) :: s
          !! A string to check
        INTEGER :: is
          !! An integer with : __1__ for string, __2__ for logical, __4__ for integer and __5__ for real
        LOGICAL                      :: dec,fdot,fexp
        INTEGER                      :: i
        CHARACTER(len=24), PARAMETER :: aset='abcfghijklmnopqrstuvwxyz'
        CHARACTER(len=10), PARAMETER :: iset='1234567890'
        CHARACTER(len=2),  PARAMETER :: dset='ed'
        CHARACTER(len=2),  PARAMETER :: sset='+-'
        CHARACTER(len=7),  PARAMETER :: slog(4) = (/'.true. ','.false.',&
                                                    't      ','f      '/)
        is = -1 ; dec = .false. ; fdot = dec ; fexp = fdot
        DO i = 1,LEN(s)
          IF (i == 1) THEN
            ! string does not start by [+-\.\d]
            IF (VERIFY(s(i:i),'.'//iset//sset) /= 0) THEN
              is = 1 ; EXIT
            ENDIF
            ! update control flag for decimal part
            dec = s(i:i) == '.' ; fdot = dec
          ELSE
            ! check if char is in [a-z]
            IF(VERIFY(s(i:i),aset) == 0) THEN
              dec=.false. ; is = 1 ; EXIT
            ELSE IF (s(i:i) == '.') THEN
              ! check for dot in decimal/exponent part (==> not a number
              IF (fdot.OR.fexp) THEN
                dec = .false. ; is = 1 ; EXIT
              ENDIF
            ELSE IF (VERIFY(s(i:i),dset)==0) THEN
              IF (fexp) THEN
                dec = .false. ; is = 1 ; EXIT
              ENDIF
            ELSE IF (VERIFY(s(i:i),sset) == 0) THEN
              IF (VERIFY(s(i-1:i-1),dset) /= 0) THEN
                dec = .false. ; is = 1 ; EXIT
              ENDIF
            ENDIF
            fdot = (fdot .OR. s(i:i) == '.')
            fexp = (fexp .OR. VERIFY(s(i:i), dset) == 0)
          ENDIF
        ENDDO
        ! it is a string
        IF (is == 1) THEN
          ! but have the format of a logical
          IF (any(slog == s)) is = 2
        ELSE
          IF ((fexp.AND.SCAN(s(LEN(s):LEN(s)),dset) /= 0)) THEN
            is = 1
          ELSE
            is = 4
            IF (fdot.OR.fexp) is = 5
          ENDIF
        ENDIF
      END FUNCTION what_
  END FUNCTION string_is

  FUNCTION format_string(str,idt1,idto) RESULT(output)
    !! Format the given string
    !!
    !! This function only replaces all '\\n' escape sequence in the given string by NEW_LINE() character.
    !! The output string is eventually indented if optional arguments are set.
    !! @warning
    !! __idto__ is relative to __idt1__ !
    CHARACTER(len=*), INTENT(in)  :: str     !! The string to format
    INTEGER, INTENT(in), OPTIONAL :: idt1, & !! An optional integer with the indentation level of the first output line (default to 0)
                                     idto    !! An optional integer with the indentation level of all other output lines (default to 0)
    CHARACTER(len=:), ALLOCATABLE :: output  !! An allocatable string with the output formatted string.
    ! - LOCAL
    INTEGER :: i,c,ti,mx
    CHARACTER(len=:), ALLOCATABLE :: idts
    IF (LEN_TRIM(str) == 0) THEN
      ALLOCATE(output,source='') ; RETURN
    ENDIF
    i=0 ; IF (PRESENT(idt1)) i = MAX(i,idt1) 
    ALLOCATE(CHARACTER(len=i) :: output) 
    IF (i > 0) output(1:i) = CHAR(32) 
    ! i0 is relative to i1 and must be >= 0
    IF (PRESENT(idto)) i = MAX(i+idto,0)
    ALLOCATE(CHARACTER(len=i+1) :: idts)
    idts(1:1) = NEW_LINE('A') ; IF (i>1) idts(2:) = CHAR(32) 
    ! Builds output string 
    c=1 ; mx = LEN_TRIM(str)
    i = INDEX(str(c:),'\n') ; ti = c+i-1
    IF (i == 0) THEN
      output=output//TRIM(str(ti+1:mx)) 
    ELSE
      output=output//TRIM(str(c:ti-1)) ; c=ti+2 
      DO
        i = INDEX(str(c:),"\n") ; ti = c+i-1
        IF (i == 0) THEN
          output=output//TRIM(str(ti+1:mx)) ; c = mx+1 
        ELSE
          output=output//idts//str(c:ti-1) ; c = ti+2
        ENDIF
        IF (c > mx) EXIT
      ENDDO
    ENDIF
    ! print a newline if we have \n at the end of the string
    IF (INDEX(TRIM(str),'\n',.true.) == mx-1.AND.TRIM(str) /= '\n') &
    output=output//idts(1:1)
  END FUNCTION format_string

  FUNCTION format_paragraph(str,width,idt1,idto) RESULT(output)
    !! Split and format a string over several lines 
    !! 
    !! The function splits an input string in words so output lines fit (almost) in __width__ characters. 
    !! The method handles indentation level (defined as leading blank spaces). It also accounts for known 
    !! csi (see [[string_op(module):attributes(variable)]]).
    !! @note
    !! Words are considered indivisible and thus output lines can sometimes exceed the maximum width if 
    !! there is not enough space to put a word (with the associated indentation if given). The default
    !! behavior in that case is to print the word in a new line (with the correct leading blank spaces).
    !! @warning
    !! If __width__, __idt1__ and/or __idto__ have inconsistent values (e.g. __width__ <= __idt1__), the
    !! method still computes the paragraph, but each words will be set on a new line with the appropriate
    !! indentation.
    CHARACTER(len=*), INTENT(in)  :: str    !! string with the content to split 
    INTEGER, INTENT(in)           :: width  !! An positive integer with the maximum width of a line
    INTEGER, INTENT(in), OPTIONAL :: idt1   !! An optional integer with the indentation level of the first output line
    INTEGER, INTENT(in), OPTIONAL :: idto   !! An optional integer with the indentation level of the other output lines
    CHARACTER(len=:), ALLOCATABLE :: output !! An allocatable string with the output content
    CHARACTER(len=:), ALLOCATABLE :: idts,zs
    INTEGER                       :: l1,lo,zmx,zw,cc,j,jj,l
    zw = abs(width) ; zs = strip_newline(str)
    zmx = LEN_TRIM(zs)
    IF (zmx == 0) THEN
      ALLOCATE(output,source='') ; RETURN
    ENDIF
    l1=0 ; IF (PRESENT(idt1)) l1 = MAX(l1,idt1)
    ALLOCATE(CHARACTER(len=l1) :: output)
    IF (l1 > 0) output(1:l1) = CHAR(32)
    lo=l1 ; IF (PRESENT(idto)) lo = MAX(l1+idto,0)
    ALLOCATE(CHARACTER(len=lo+1) :: idts)
    idts(1:1) = NEW_LINE('A') ; IF (lo>=1) idts(2:len(idts)) = CHAR(32)
    ! Prints a message if user is just stupid...
    IF (lo+1 > zw .OR. l1+1 > zw) THEN
      output = str ; RETURN
    ENDIF
    ! check if can just return the string as is 
    IF (zmx + l1 <= zw) THEN
      output=output//TRIM(zs) ; RETURN
    ENDIF
    j=1 ; jj=1+l1 
    DO 
      ! Gets next blank in input string
      cc = INDEX(TRIM(zs(j:)),CHAR(32))
      ! no more blank
      ! Gets total length of csi between zs(j:j+cc-1)
      ! this value will be substracted to each length test
      IF (cc == 0) THEN
        l = csis_length(zs(j:)) 
        IF (jj-1+LEN_TRIM(zs(j:))-l > zw) THEN
          output = output//idts
        ENDIF
        output=output//TRIM(zs(j:))
        EXIT ! we are at the last word : we must exit the infinite loop !
      ELSE
        l = csis_length(zs(j:j+cc-1)) 
        IF (cc+jj-1-l > zw) THEN
          output=output//idts//zs(j:j+cc-1) ; jj = lo+1+cc+1 - l
        ELSE
          output=output//zs(j:j+cc-1) ; jj = jj + cc - l
        ENDIF
      ENDIF
      j = j + cc
    ENDDO
    CONTAINS
    FUNCTION csis_length(str) RESULT(value)
      ! - DUMMY
      CHARACTER(len=*), INTENT(in) :: str
      ! - RESULT
      INTEGER :: value
      ! - LOCAL
      INTEGER :: jc,iesc,im
      LOGICAL :: tcsi
      value = 0 
      jc=1
      DO 
        IF (jc>LEN(str)) EXIT
        ! search for escape
        iesc = INDEX(str(jc:),CHAR(27))
        IF (iesc == 0) EXIT 
        ! search for m
        im = INDEX(str(jc+iesc:),"m")
        ! no m in the string after ESC --> this could not be a csi
        IF (im == 0) EXIT
        ! check if this is really a csi and updates length
        tcsi = is_csi(str(jc+iesc-1:jc+iesc+im-1))
        jc = jc + iesc 
        IF (tcsi) THEN
          value=value+im+1
          jc=jc+im
        ENDIF
      ENDDO
    END FUNCTION csis_length
  END FUNCTION format_paragraph

  FUNCTION strip_newline(str,rpl) RESULT(stripped)
    !! Replace newline escape sequences by spaces
    !!
    !! The function replaces newline (both '\\n' escape sequence and Fortran NEW_LINE() character) in the 
    !! given string and returns the resulting string.
    CHARACTER(len=*), INTENT(in)           :: str !! A string to process
    CHARACTER(len=1), INTENT(in), OPTIONAL :: rpl !! A optional single character used as substitution of escape sequences (blank space by default)
    CHARACTER(len=:), ALLOCATABLE :: stripped     !! An allocatable string with all newline sequences replaced by blank space or __rpl__ if given 
    CHARACTER(len=1) :: zrp
    INTEGER          :: i, j, ns 
    zrp = CHAR(32) ; IF(PRESENT(rpl)) zrp = rpl
    IF (str == NEW_LINE('A')) THEN
      stripped = zrp ; RETURN 
    ENDIF
    ns = LEN_TRIM(str)
    IF (ns == 0) THEN
      ALLOCATE(stripped,source='') ; RETURN
    ENDIF
    ALLOCATE(CHARACTER(len=ns) :: stripped) ; stripped(1:ns) = CHAR(32)
    i=1 ; j=1
    DO 
      IF (str(i:i) == NEW_LINE('A')) THEN
        stripped(j:j) = zrp 
      ELSE IF (i < ns) THEN
          IF (str(i:i+1) == "\n") THEN
            stripped(j:j) = zrp ; i=i+1
          ELSE
            stripped(j:j) = str(i:i) 
          ENDIF
      ELSE
        stripped(j:j) = str(i:i) 
      ENDIF
      j=j+1 ; i=i+1
      IF (i > ns .OR. j > ns) EXIT
    ENDDO
    IF (j < ns) stripped = stripped(1:j)
    RETURN
  END FUNCTION strip_newline

  FUNCTION str_length(str) RESULT(res)
    !! Get the length of the string object
    !! 
    !! The method computes the length of the string. It differs from LEN intrinsic function as
    !! it does not account for extra-characters of csi codes.
    CHARACTER(len=*), INTENT(in) :: str !! String to process
    INTEGER :: res                      !! The actual length of string (i.e. does not account for csi codes)
    CHARACTER(len=:), ALLOCATABLE :: tmp
    res = 0 
    IF (LEN(str) /= 0) THEN
      tmp = reset_csi(str)
      res = LEN(tmp)
      DEALLOCATE(tmp)
    ENDIF
    RETURN
  END FUNCTION str_length

  FUNCTION to_lower(str1) RESULT(str)
    !! Convert the string in lower case
    !!
    !! The method converts the input string in lower case and accounts for
    !! possible csi codes in the string.
    CHARACTER(len=*), INTENT(in) :: str1 !! Input string to convert
    CHARACTER(len=:), ALLOCATABLE :: str !! A copy of the string in lower case
    INTEGER :: i,ic
    IF (LEN(str1) /= 0) THEN
      str = str1
      DO i = 1, len(str1)
        ic = ichar(str1(i:i))
        IF (ic >= 65 .AND. ic < 90) str(i:i) = char(ic + 32)
      ENDDO
    ELSE
      str=''
    ENDIF
  END FUNCTION to_lower

  FUNCTION to_upper(str1) RESULT(str)
    !! Convert the string in upper case
    !!
    !! The method converts the input string in upper case and accounts for
    !! possible csi codes in the string.
    CHARACTER(len=*), INTENT(in) :: str1 !! Input string to convert
    CHARACTER(len=:), ALLOCATABLE :: str !! A copy of the string in upper case
    INTEGER :: j,i,ic,icsi,lcsi
    IF (LEN(str1) > 0) THEN
      str = str1 
      i = 1
      DO
        IF (i > LEN(str)) EXIT
        icsi = str_index_of_csi(str(i:),lcsi)
        IF (icsi == 0) THEN
          ! no more csi the end of string is upper case converted
          DO j=i,LEN(str)
            ic = ichar(str(j:j))
            IF (ic >= 97 .AND. ic < 122) str(j:j) = char(ic-32)
          ENDDO
          RETURN
        ELSE IF (icsi == 1) THEN
          i = i + lcsi
        ELSE IF (icsi > 1) THEN
          ! csi is not the first word: we convert in upper case until its
          ! position THEN copy the csi and get back in the loop
          DO j=i,i+icsi-2
            ic = ichar(str(j:j))
            IF (ic >= 97 .AND. ic < 122) str(j:j) = char(ic-32)
          ENDDO
          i = i + icsi + lcsi-1 
        ENDIF
      ENDDO
    ELSE
      str=''
    ENDIF
  END FUNCTION to_upper

 FUNCTION str_remove(string,substring,back,all) RESULT(str)
   !! Remove substring from current string 
   !! 
   !! The function removes the first occurence of __substring__ in __string__ or all 
   !! its occurences if __all__ is explicitly set to .true..
    CHARACTER(len=*), INTENT(in)  :: string    !! A string to search in
    CHARACTER(len=*), INTENT(in)  :: substring !! A string to search and removes from __string__
    LOGICAL, INTENT(in), OPTIONAL :: back, &   !! An optional boolean flag with .true. to begin search at the end of the string
                                     all       !! An optional boolean flag with .true. to remove all occurences of __substring__
    CHARACTER(len=:), ALLOCATABLE :: str       !! An allocatable string with __substring__ occurence(s) removed
    LOGICAL :: zb,za
    INTEGER :: is,j,zboff
    str=''
    zb = .false. ; za = .false.
    IF (PRESENT(back)) zb = back
    IF (PRESENT(all)) za = all
    IF (za) zb=.false.
    zboff = 0 ; IF (zb) zboff = 1
    IF (LEN(string) == 0) RETURN
    j=1 
    DO 
      IF (j>LEN(string)) EXIT
      ! search for substring
      is = INDEX(string(j:),substring,back=zb)
      IF (is == 0) THEN
        ! substring is not found : we get the last part of the string and return
        str = str//string(j:) ; RETURN 
      ELSE IF (is == 1) THEN
        j = j + LEN(substring)
      ELSE
        ! substring is not at the begin of the string : saves the string
        str = str//string(j:j+is-2)
        j = j + is+LEN(substring)-1
      ENDIF
      ! if we only want to str_remove ONE occurence we exit if substring
      ! has been found
      IF (.NOT.(is==0.OR.za)) EXIT 
    ENDDO
    IF (j <= LEN(string).AND..NOT.zb) str=str//string(j:)
    RETURN 
  END FUNCTION str_remove

 FUNCTION str_replace(string,old,new,back,all) RESULT(str)
    !! Replace substring from current string 
    !!
    !! The function replaces the first occurence of __old__ in __string__ by 
    !! __new__ or all its occurence(s) if __all__ is explicitly set to .true..
    CHARACTER(len=*), INTENT(in)  :: string  !! A string to search in
    CHARACTER(len=*), INTENT(in)  :: old,  & !! A string to search and replace
                                     new     !! A string to substitute to __old__
    LOGICAL, INTENT(in), OPTIONAL :: back, & !! An optional boolean flag with .true. to begin search at the end of the string
                                     all     !! An optional boolean flag with .true. to replace all occurences of __old__
    CHARACTER(len=:), ALLOCATABLE :: str     !! An allocatable string with occurence(s) of __old__ replaced by __new__
    LOGICAL :: zb,za
    INTEGER :: is,j
    str=''
    zb = .false. ; za = .false.
    IF (PRESENT(back)) zb = back
    IF (PRESENT(all)) za = all
    IF (za) zb = .NOT.za
    IF (LEN(string) == 0) RETURN 
    j=1 
    DO 
      IF (j>LEN(string)) EXIT
      ! search for "old"
      is = INDEX(string(j:),old,back=zb)
      IF (is == 0) THEN
        ! "old" is not found : we get the last part of the string and return
        str = str//string(j:) ; RETURN 
      ELSE IF (is == 1) THEN
        str = str//new
        j = j + LEN(old)
      ELSE
        ! "old" is not at the begin of the string : saves the string
        str = str//string(j:j+is-2)//new
        j = j + is + LEN(old) - 1 
      ENDIF
      IF (.NOT.(is==0.OR.za)) EXIT 
    ENDDO
    IF (j <= LEN(str)) str=str//string(j:)
    RETURN 
  END FUNCTION str_replace

  FUNCTION endswith(string,substring,icase) RESULT(ret)
    !! Check if string ends by substring 
    CHARACTER(len=*), INTENT(in)  :: string
      !! @param[in] string A string to check
    CHARACTER(len=*), INTENT(in)  :: substring
      !! A string to search in __string__
    LOGICAL, INTENT(in), OPTIONAL :: icase 
      !! An optional boolean flag with .true. to perform insensitive case search
    LOGICAL :: ret
      !! .true. if __string__ ends by __substring__, .false. otherwise.
    CHARACTER(len=:), ALLOCATABLE :: zthis,zstr
    INTEGER                       :: idx 
    LOGICAL                       :: noc 
    ret = .false.
    noc = .false. ; IF (PRESENT(icase)) noc = icase
    IF (LEN(string) == 0 .OR. LEN(substring) == 0) RETURN
    zthis = reset_csi(string) ; zstr=reset_csi(substring)
    IF (noc) THEN
      idx = INDEX(to_lower(zthis),to_lower(zstr),.true.)
    ELSE
      idx = INDEX(zthis,zstr,.true.)
    ENDIF
    IF (idx == 0.OR.idx+str_length(zstr)-1 /= str_length(zthis)) RETURN
    ret=.true.
  END FUNCTION endswith

  FUNCTION startswith(string,substring,icase) RESULT(ret)
    !! Check if string starts by substring 
    CHARACTER(len=*), INTENT(in)  :: string
      !! A string to check
    CHARACTER(len=*), INTENT(in)  :: substring
      !! A string to search in __string__
    LOGICAL, INTENT(in), OPTIONAL :: icase 
      !! An optional boolean flag with .true. to perform insensitive case search
    LOGICAL :: ret
      !! .true. if __string__ starts by __substring__, .false. otherwise.
    CHARACTER(len=:), ALLOCATABLE :: zthis,zstr
    INTEGER                       :: idx 
    LOGICAL                       :: noc 
    ret = .false.
    noc = .false. ; IF (PRESENT(icase)) noc = icase
    IF (LEN(string) == 0 .OR. LEN(substring) == 0) RETURN
    zthis = reset_csi(string) ; zstr=reset_csi(substring)
    IF (noc) THEN
      idx = INDEX(to_lower(zthis),to_lower(zstr))
    ELSE
      idx = INDEX(zthis,zstr)
    ENDIF
    IF (idx /= 1) RETURN
    ret=.true.
  END FUNCTION startswith

  ! CSI related functions 
  ! ---------------------

  FUNCTION add_csi(string,attrs) RESULT(str)
    !! Set csi attributes to the given string object
    !!
    !! The function adds csi (ANSI escape sequences) to the given string and
    !! returns a copy of it.
    CHARACTER(len=*), INTENT(in)      :: string
      !! @param[in] string A string object reference 
    INTEGER, INTENT(in), DIMENSION(:) :: attrs
      !! A vector of integers with the code to add. Each __attrs__ value should refers to one i
      !! of [[string_op(module):attributes(variable)]] values.
    CHARACTER(len=:), ALLOCATABLE :: str
      !! An allocatable string with new csi codes added.
    INTEGER                       :: j,iesc,im
    CHARACTER(len=:), ALLOCATABLE :: tmp,csi
    CHARACTER(len=4), PARAMETER   :: rcsi = CHAR(27)//"[0m"
    str=''
    ! 1) Check for input string
    IF (LEN(string) == 0) RETURN
    ! 2) Removes last <ESC>[0m if any and initializes output string
    ! we must remove only the last <ESC>[0m if any
    IF (INDEX(string,rcsi,.true.) == LEN(string)-3) THEN
      tmp = str_remove(string,rcsi,back=.true.)
    ELSE
      tmp = string
    ENDIF
    ! 3) Add all the given csi preceded by <ESC>[0m at the beginning of the string 
    !    if it does not start by an ANSI sequence
    IF (INDEX(tmp,CHAR(27)//"[") /= 1) &
    tmp = str_add_to_csi(rcsi,attrs)//tmp
    ! Loops on new string and updates csi codes
    j=1 
    DO 
      IF (j>LEN(tmp)) EXIT
      ! search for escape
      iesc = INDEX(tmp(j:),CHAR(27))
      IF (iesc == 0) THEN
        ! no more ESC : cat until end of input string and exit
        str = str//tmp(j:) ; EXIT
      ELSE IF (iesc > 1) THEN
        ! ESC is not first char: copy until ESC
        str = str//tmp(j:j+iesc-2)
      ENDIF
      ! search for m
      im = INDEX(tmp(j+iesc:),"m")
      ! no m in the string after ESC --> copy string (INCLUDING ESC) and leave
      IF (im == 0) THEN
        str = str//tmp(j+iesc-1:)
        RETURN
      ENDIF
      csi = tmp(j+iesc-1:j+iesc+im-1)
      ! we have a csi: we add new codes to it
      IF (is_csi(csi)) THEN
        csi = str_add_to_csi(csi,attrs)
      ENDIF
      str = str//csi
      j = j + iesc + im
    ENDDO
    IF (INDEX(str,rcsi,.true.) /= LEN(str)-3) str = str//rcsi
    RETURN 
  END FUNCTION add_csi

  FUNCTION del_csi(string,attrs) RESULT(str)
    !! Remove attributes to the given string
    !!
    !! The function removes list of csi (ANSI escape sequences) from the given 
    !! string and returns a copy of it.
    CHARACTER(len=*), INTENT(in)      :: string
      !! Input string 
    INTEGER, INTENT(in), DIMENSION(:) :: attrs 
      !! A vector of integers with the code to remove. Each __attrs__ value should 
      !! refers to one of [[string_op(module):attributes(variable)]] values.
    CHARACTER(len=:), ALLOCATABLE :: str
      !! An allocatable string with csi codes from __list__ removed
    LOGICAL                                           :: ok
    INTEGER                                           :: j,iesc,im
    CHARACTER(len=:), ALLOCATABLE                     :: tmp,csi,csis
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: tks
    CHARACTER(len=4), PARAMETER                       :: rcsi = CHAR(27)//"[0m"
    str=''
    IF (LEN(string) == 0) RETURN
    ! remove last <ESC>[0m if found at the end of the string
    IF (INDEX(string,rcsi,.true.) == LEN(string)-3) THEN
      tmp = str_remove(string,rcsi,back=.true.)
    ELSE
      tmp = string
    ENDIF
    ! Loops on new string and updates csi codes
    j=1 ; csis=""
    DO 
      IF (j>LEN(tmp)) EXIT
      ! search for escape
      iesc = INDEX(tmp(j:),CHAR(27))
      IF (iesc == 0) THEN
        ! no more ESC : cat until end of input string and exit
        str = str//tmp(j:) ; EXIT
      ELSE IF (iesc > 1) THEN
        ! ESC is not first char: copy until ESC
        str = str//tmp(j:j+iesc-2)
      ENDIF
      ! search for m
      im = INDEX(tmp(j+iesc:),"m")
      ! no m in the string after ESC --> copy string (INCLUDING ESC) and leave
      IF (im == 0) THEN
        str = str//tmp(j+iesc-1:)
        RETURN
      ENDIF
      csi = tmp(j+iesc-1:j+iesc+im-1)
      ! we have a csi: we add new codes to it
      IF (is_csi(csi)) THEN
        csi = str_del_from_csi(csi,attrs)
      ENDIF
      csis=csis//csi//"|"
      str = str//csi
      j = j + iesc + im
    ENDDO
    ! Add <ESC>[0m at the end of string if not found
    IF (INDEX(str,rcsi,.true.) /= LEN(str)-3) str = str//rcsi
    ! resets all attributes if we only have <ESC>[0m in final list 
    ok = tokenize(csis(1:LEN(csis)-1),tks,"|") 
    IF (ALL(tks == rcsi)) str = reset_csi(str)
    DEALLOCATE(tks)
    RETURN 
  END FUNCTION del_csi

  FUNCTION reset_csi(string) RESULT(str)
    !! Reset all csi codes of the string
    !! 
    !! The method removes __all__ the known escape sequences from the input string.
    CHARACTER(len=*), INTENT(in) :: string
      !! Input string
    CHARACTER(len=:), ALLOCATABLE :: str 
      !! An allocatable string with the copy of input string stripped off csi codes.
    INTEGER :: j,iesc,im
    LOGICAL :: tcsi
    str = ""
    IF (LEN(string) == 0) RETURN 
    j=1 
    DO 
      IF (j>LEN(string)) EXIT
      ! search for escape
      iesc = INDEX(string(j:),CHAR(27))
      IF (iesc == 0) THEN
        str = str//string(j:) ; EXIT
      ENDIF
      ! search for m
      im = INDEX(string(j+iesc:),"m")
      ! no m in the string after ESC --> copy string (INCLUDING ESC) and leave
      IF (im == 0) THEN
        str = str//string(j+iesc-1:)
        RETURN
      ENDIF
      ! csi includes everything between ESC and m (excluding them):
      ! to check for csi it should begin by [ and then be a list of integers
      ! separated by ;
      tcsi = is_csi(string(j+iesc-1:j+iesc+im-1))
      IF (iesc > 1) THEN
        str = str//string(j:j+iesc-2)
      ENDIF
      j = j + iesc ; IF (tcsi) j=j+im
    ENDDO
    RETURN 
  END FUNCTION reset_csi

  FUNCTION is_csi(value) RESULT(yes)
    !! Check if string is a known csi
    !! 
    !! The function only check for known csi code which are defined in [[string_op(module):attributes(variable)]].
    CHARACTER(len=*), INTENT(in) :: value
      !! A Fortran intrinsic string to check
    LOGICAL :: yes
      !! .true. if it is a known csi, .false. otherwise
    LOGICAL                                           :: ok
    CHARACTER(len=:), ALLOCATABLE                     :: tmp 
    TYPE(words)                                       :: wtks
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: stks
    INTEGER, DIMENSION(:), ALLOCATABLE                :: nums
    INTEGER                                           :: i
    yes = .false.
    IF (LEN(value) < 4) RETURN
    tmp = value(3:len(value)-1)
    call words_extend(wtks,tmp,';')
    ok = words_to_vector(wtks,stks)
    CALL ws_clear_sc(wtks)
    IF (.NOT.ok) RETURN
    ! if we cannot convert strings to integers : it is not a csi
    IF (.NOT.from_string(stks,nums)) RETURN
    DEALLOCATE(stks)
    DO i=1, SIZE(nums)
      IF (.NOT.ANY(attributes == nums(i))) RETURN
    ENDDO
    yes = .true.
  END FUNCTION is_csi

  FUNCTION str_add_to_csi(csi,list) RESULT(ncsi)
    !! Add a new list of codes to the input csi string
    !! 
    !! The method adds all the csi codes given in __list__ that are known by the module and not
    !! already present in the input csi.
    CHARACTER(len=*), INTENT(in)      :: csi
      !! A string with the input csi. It __must__ begin with "<ESC>[" and ends with "m".
    INTEGER, INTENT(in), DIMENSION(:) :: list
      !! A vector of integers with the csi code to add. Each value of __list__ should be one of 
      !! [[string_op(module):attributes(variable)]] values. All unknown values are filtered out as well 
      !! as csi code already present in input __csi__.
    CHARACTER(len=:), ALLOCATABLE :: ncsi 
      !! A new csi string or the input __csi__ if some "errors" occured (the input csi could not 
      !! be tokenized or none of __list__ values are left after filtering).
    LOGICAL                                            :: ok 
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE  :: tks
    CHARACTER(len=:), ALLOCATABLE                      :: tmp
    INTEGER, DIMENSION(:), ALLOCATABLE                 :: zlist,nums
    INTEGER                                            :: i,j,ni,no
    ! 1) Filter input list :
    ! 1.1) Gets the list of current csi codes 
    ncsi = csi(3:len(csi)-1) 
    ok = tokenize(ncsi,tks,"; ",merge=.true.)
    IF (.NOT.from_string(tks,nums)) THEN
      ncsi = csi
      RETURN
    ENDIF
    DEALLOCATE(tks)
    ! 1.2) Filter input list of new flags to add 
    ! counts number of valid flags
    j=0 
    DO i=1,SIZE(list) 
      ! new flags must be in attributes but NOT in nums
      IF (ANY(attributes==list(i).AND..NOT.ANY(nums == list(i)))) j=j+1
    ENDDO
    ! No "valid" flags -> returns old csi
    IF (j == 0) THEN ; ncsi = csi ; RETURN ; ENDIF
    ni = SIZE(nums) ; no = j + ni
    ALLOCATE(zlist(no)) ; zlist(1:ni) = nums(:) ; j = ni
    DO i=1,SIZE(list) 
      ! new flags must be in attributes but NOT in nums
      IF (ANY(attributes==list(i).AND..NOT.ANY(nums == list(i)))) THEN
        j=j+1 ; zlist(j) = list(i)
      ENDIF 
    ENDDO
    DEALLOCATE(nums)
    ! 2) Builds new csi
    !    Here we explictly set the first flag to 0 (i.e. reset attributes)...
    ncsi = CHAR(27)//"[0;"
    DO i=1,no
      ! ... So we get rid of all "0" flag in the list 
      IF (zlist(i) /= 0) THEN
        tmp = to_string(zlist(i))
        IF (LEN_TRIM(tmp) == 0) THEN
          ncsi = csi ; RETURN
        ENDIF
        ncsi = ncsi//tmp
        IF (i /= no) ncsi = ncsi//";"
      ENDIF
    ENDDO
    ncsi = ncsi//"m" 
  END FUNCTION str_add_to_csi

  FUNCTION str_del_from_csi(csi,list) RESULT(ncsi)
    !! Remove a list of codes from the input csi string
    !!
    !! The method removes all the csi codes given in __list__ that are known by the
    !! module and already present in the input csi.
    CHARACTER(len=*), INTENT(in)      :: csi
      !! An intrinsic Fortran string with the input csi. It __must__ begin with "<ESC>[" and ends with "m".
    INTEGER, INTENT(in), DIMENSION(:) :: list
      !! A vector of integers with the csi code to remove. Each value of __list__ should be one of 
      !! [[string_op(module):attributes(variable)]] values. All unknown values are filtered out.
    CHARACTER(len=:), ALLOCATABLE :: ncsi 
      !! A new csi string or the input __csi__ if some "errors" occured (the input csi could not 
      !! be tokenized or none of __list__ values are left after filtering).
    LOGICAL                                            :: ok
    CHARACTER(len=LEN(csi)), DIMENSION(:), ALLOCATABLE :: tks
    CHARACTER(len=:), ALLOCATABLE                      :: tmp
    INTEGER, DIMENSION(:), ALLOCATABLE                 :: nums
    INTEGER                                            :: i
    ncsi = csi(3:len(csi)-1) 
    ok = tokenize(ncsi,tks,"; ",merge=.true.) 
    IF (.NOT.from_string(tks,nums)) THEN 
      ncsi = csi 
      RETURN 
    ENDIF
    DEALLOCATE(tks)
    tmp=""
    DO i=1, SIZE(nums)
      IF (ALL(nums(i) /= list).AND.nums(i) /= 0) THEN
        ! no need to check for to_string status : it is always ok !
        tmp=tmp//to_string(nums(i))//";"
      ENDIF
    ENDDO
    IF (LEN_TRIM(tmp) /= 0) THEN
      ncsi=CHAR(27)//"[0;"//tmp(1:LEN(tmp)-1)//"m"
    ELSE
      ncsi=CHAR(27)//"[0m"
    ENDIF
  END FUNCTION str_del_from_csi

  FUNCTION str_index_of_csi(str,length) RESULT(pos)
    !! Get the position of the first known csi in string
    !!
    !! The method searches for the first known csi in string. The csi must contain known codes 
    !! (i.e. values of [[string_op(module):attributes(variable)]]).
    CHARACTER(len=*), INTENT(in) :: str    !! A string to search in
    INTEGER, INTENT(out)         :: length !! Length of the csi in the string
    INTEGER                      :: pos    !! Position of the first csi found. It is set to 0 if no csi has been found.
    INTEGER :: iesc,im
    pos = 0 ; length = 0
    ! we need at least 4 chars to create a csi
    IF (LEN_TRIM(str) < 4) RETURN 
    iesc = INDEX(str,CHAR(27))
    IF (iesc == 0) RETURN
    ! search for m
    im = INDEX(str(iesc:),"m")
    ! no m in the string after ESC --> copy string (INCLUDING ESC) and leave
    IF (im == 0) RETURN
    IF (.NOT.is_csi(str(iesc:iesc+im-1))) RETURN
    pos = iesc ; length = im
  END FUNCTION str_index_of_csi

  ! String conversion functions
  ! ---------------------------

  FUNCTION str2int_sc(str, value) RESULT(ret)
    !! Convert string value to integer value (scalar)
    CHARACTER(len=*), INTENT(in) :: str   !! String to convert
    INTEGER, INTENT(out)         :: value !! Output value
    LOGICAL :: ret                        !! Return status (.true. on success)
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; zs = remove_quotes(str)
    IF (string_is(zs) /= st_integer) THEN
      ret = .false.
    ELSE
      READ(zs, *) value
    ENDIF
    RETURN
  END FUNCTION str2int_sc

  FUNCTION str2log_sc(str, value) RESULT(ret)
    !! Convert string value to logical value (scalar)
    CHARACTER(len=*), INTENT(in) :: str   !! String to convert
    LOGICAL, INTENT(out)         :: value !! Output value
    LOGICAL :: ret                        !! Return status (.true. on success)
    CHARACTER(len=:), ALLOCATABLE :: zs
    integer :: r
    ret = .true. ; zs = remove_quotes(str)
    r = string_is(zs)
    IF (string_is(zs) /= st_logical) THEN
      ret = .false.
    ELSE
      READ(zs, *) value
    ENDIF
    RETURN
  END FUNCTION str2log_sc

  FUNCTION str2real_sc(str, value) RESULT(ret)
    !! Convert string value to simple precision floating precision value (scalar)
    CHARACTER(len=*), INTENT(in) :: str   !! String to convert
    REAL(kind=4), INTENT(out)    :: value !! Output value
    LOGICAL :: ret                        !! Return status (.true. on success)
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true.; zs = remove_quotes(str)
    IF (string_is(zs) < st_integer) THEN
      ret = .false.
    ELSE
      READ(zs, *) value
    ENDIF
    RETURN
  END FUNCTION str2real_sc

  FUNCTION str2dble_sc(str, value) RESULT(ret)
    !! Convert string value to double precision floating precision value (scalar)
    CHARACTER(len=*), INTENT(in) :: str   !! String to convert
    REAL(kind=8), INTENT(out)    :: value !! Output value
    LOGICAL :: ret                        !! Return status (.true. on success)
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; zs = remove_quotes(str)
    IF (string_is(zs) < st_integer) THEN
      ret = .false.
    ELSE
      READ(zs, *) value
    ENDIF
    RETURN
  END FUNCTION str2dble_sc

  FUNCTION str2cplx_sc(str, value) RESULT(ret)
    !! Convert string value to complex value (scalar)
    CHARACTER(len=*), INTENT(in) :: str   !! String to convert
    COMPLEX(kind=4), INTENT(out) :: value !! Output value
    LOGICAL :: ret                        !! Return status (.true. on success)
    ! - LOCAL 
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; zs = remove_quotes(str)
    IF (string_is(zs) /= st_complex) THEN
      ret = .false.
    ELSE
      READ(zs, *) value
    ENDIF
    RETURN
  END FUNCTION str2cplx_sc

  FUNCTION str2int_ve(str, value) RESULT(ret)
    !! Convert strings values to integer values (vector)
    CHARACTER(len=*), INTENT(in), DIMENSION(:)      :: str   !! Vector of strings to convert
    INTEGER, INTENT(out), DIMENSION(:), ALLOCATABLE :: value !! Vector of output values
    LOGICAL :: ret                                           !! Return status (.true. on success)
    INTEGER                       :: i,ns
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; ns = SIZE(str) ; ALLOCATE(value(ns))
    DO i=1,ns
      zs = remove_quotes(str(i))
      IF (string_is(zs) /= st_integer) THEN
        ret = .false. ; DEALLOCATE(value) ; RETURN
      ELSE
        READ(zs, *) value(i)
      ENDIF
    ENDDO
    RETURN
  END FUNCTION str2int_ve

  FUNCTION str2log_ve(str, value) RESULT(ret)
    !! Convert strings values to logical values (vector)
    CHARACTER(len=*), INTENT(in), DIMENSION(:)      :: str   !! Vector of strings to convert
    LOGICAL, INTENT(out), DIMENSION(:), ALLOCATABLE :: value !! Vector of output values
    LOGICAL :: ret                                           !! Return status (.true. on success)
    INTEGER                       :: i,ns
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; ns = SIZE(str) ; ALLOCATE(value(ns))
    DO i=1,ns
      zs = remove_quotes(str(i))
      IF (string_is(zs) /= st_logical) THEN
        ret = .false. ; DEALLOCATE(value) ; RETURN
      ELSE
        READ(zs, *) value(i)
      ENDIF
    ENDDO
    RETURN
  END FUNCTION str2log_ve

  FUNCTION str2real_ve(str, value) RESULT(ret)
    !! Convert strings values to simple precision floating point values (vector)
    CHARACTER(len=*), INTENT(in), DIMENSION(:)           :: str   !! Vector of strings to convert
    REAL(kind=4), INTENT(out), DIMENSION(:), ALLOCATABLE :: value !! Vector of output values
    LOGICAL :: ret                                                !! Return status (.true. on success) 
    INTEGER                       :: i,ns
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; ns = SIZE(str) ; ALLOCATE(value(ns))
    DO i=1,ns
      IF (string_is(zs) < st_integer) THEN
        ret = .false. ; DEALLOCATE(value) ; RETURN
      ELSE
        READ(zs, *) value(i)
      ENDIF
    ENDDO
    RETURN
  END FUNCTION str2real_ve

  FUNCTION str2dble_ve(str, value) RESULT(ret)
    !! Convert strings values to double precision floating point values (vector)
    CHARACTER(len=*), INTENT(in), DIMENSION(:)           :: str   !! Vector of strings to convert
    REAL(kind=8), INTENT(out), DIMENSION(:), ALLOCATABLE :: value !! Vector of output values
    LOGICAL :: ret                                                !! Return status (.true. on success)
    INTEGER                       :: i,ns
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; ns = SIZE(str) ; ALLOCATE(value(ns))
    DO i=1,ns
      zs = remove_quotes(str(i))
      IF (string_is(zs) < st_integer) THEN
        ret = .false. ; DEALLOCATE(value) ; RETURN
      ELSE
        READ(zs, *) value(i)
      ENDIF
    ENDDO
    RETURN
  END FUNCTION str2dble_ve

  FUNCTION str2cplx_ve(str, value) RESULT(ret)
    !! Convert strings values to complex values (vector)
    CHARACTER(len=*), INTENT(in), DIMENSION(:)              :: str   !! Vector of strings to convert
    COMPLEX(kind=4), INTENT(out), DIMENSION(:), ALLOCATABLE :: value !! Vector of output values
    LOGICAL :: ret                                                   !! Return status (.true. on success)
    INTEGER                       :: i,ns
    CHARACTER(len=:), ALLOCATABLE :: zs
    ret = .true. ; ns = SIZE(str) ; ALLOCATE(value(ns))
    DO i=1,ns
      zs = remove_quotes(str(i))
      IF (string_is(zs) /= st_complex) THEN
        ret = .false. ; DEALLOCATE(value) ; RETURN
      ELSE
        READ(zs, *) value(i)
      ENDIF
    ENDDO
    RETURN
  END FUNCTION str2cplx_ve

  FUNCTION int2str_as(value) RESULT(str)
    !! Convert an integer value to string (auto format / string result) 
    INTEGER, INTENT(in)           :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=DIGITS(value)) :: str)
    WRITE(str,*,iostat=err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION int2str_as

  FUNCTION log2str_as(value) RESULT(str)
    !! Convert a logical value to string (auto format / string result) 
    LOGICAL, INTENT(in)           :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=2) :: str)
    WRITE(str, *, IOSTAT = err) value
    str=TRIM(ADJUSTL(str))
    IF (err /= 0) str = ''
    RETURN
  END FUNCTION log2str_as

  FUNCTION real2str_as(value) RESULT(str)
    !! Convert a simple precision floating point value to string (auto format / string result) 
    REAL(kind=4), INTENT(in)      :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=DIGITS(value)) ::str)
    WRITE(str,*, IOSTAT = err) value
    str=TRIM(ADJUSTL(str))
    IF (err /= 0)  str = '' 
    RETURN
  END FUNCTION real2str_as

  FUNCTION dble2str_as(value) RESULT(str)
    !! Convert a double precision floating point value to string (auto format / string result) 
    REAL(kind=8), INTENT(in)      :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=DIGITS(value)) ::str)
    WRITE(str,*, IOSTAT = err) value
    str=TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION dble2str_as

  FUNCTION cplx2str_as(value) RESULT(str)
    !! Convert a complex value to string (auto format / string result) 
    COMPLEX(kind=4), INTENT(in)   :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err,sl
    sl = DIGITS(REAL(value))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str, *, IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION cplx2str_as

  FUNCTION dcplx2str_as(value) RESULT(str)
    !! Convert a complex value to string (auto format / string result) 
    COMPLEX(kind=8), INTENT(in)   :: value !! Value to convert
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err,sl
    sl = DIGITS(REAL(value))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str, *, IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION dcplx2str_as

  FUNCTION int2str_fs(value, fmt) RESULT(str)
    !! Convert an integer value to string (user format / string result) 
    INTEGER, INTENT(in)           :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION int2str_fs

  FUNCTION log2str_fs(value, fmt) RESULT(str)
    !! Convert a logical value to string (user format / string result) 
    LOGICAL, INTENT(in)           :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str=TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION log2str_fs

  FUNCTION real2str_fs(value, fmt) RESULT(str)
    !! Convert a simple precision floating point value to string (user format / string result) 
    REAL(kind=4), INTENT(in)      :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION real2str_fs

  FUNCTION dble2str_fs(value, fmt) RESULT(str)
    !! Convert a double precision floating point value to string (user format / string result) 
    REAL(kind=8), INTENT(in)      :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION dble2str_fs

  FUNCTION cplx2str_fs(value, fmt) RESULT(str)
    !! Convert a complex value to string (user format / string result) 
    COMPLEX(kind=4), INTENT(in)   :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION cplx2str_fs

  FUNCTION dcplx2str_fs(value, fmt) RESULT(str)
    !! Convert a complex value to string (user format / string result) 
    COMPLEX(kind=8), INTENT(in)   :: value !! Value to convert
    CHARACTER(len=*), INTENT(in)  :: fmt   !! String format
    CHARACTER(len=:), ALLOCATABLE :: str   !! String with the converted value in output
    INTEGER :: err
    ALLOCATE(CHARACTER(len=st_slen) :: str)
    WRITE(str, '('//fmt//')', IOSTAT = err) value
    str = TRIM(ADJUSTL(str))
    IF (err /= 0) str = '' 
    RETURN
  END FUNCTION dcplx2str_fs

  ! Extended strings features
  ! ---------------------------

  FUNCTION fis_cat_int(str1,int2) RESULT(str)
    !! Concatenate a fortran intrinsic string with a integer.
    CHARACTER(len=*), INTENT(in)  :: str1 !! String to concatenate
    INTEGER, INTENT(in)           :: int2 !! Integer to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str  !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(int2)) :: str) 
    WRITE(str,*) int2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str1//str
    RETURN
  END FUNCTION fis_cat_int

    !! @param[in] int2 An integer to concatenate
    !! @param[in] str1 A string to concatenate
    !! @return An allocatable string with the concatenation of input values.
  FUNCTION fis_cat_int_inv(int2,str1) RESULT(str)
    !! Concatenate a fortran intrinsic string with a integer (reversed).
    INTEGER, INTENT(in)           :: int2 !! Integer to concatenate 
    CHARACTER(len=*), INTENT(in)  :: str1 !! String to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str  !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(int2)) :: str) 
    WRITE(str,*) int2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_int_inv

  FUNCTION fis_cat_bool(str1,bool2) RESULT(str)
    !! Concatenate a string with a logical
    CHARACTER(len=*), INTENT(in)  :: str1  !! String to concatenate
    LOGICAL, INTENT(in)           :: bool2 !! Logical to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str   !! Output string
    CHARACTER(len=2) ::tmp
    WRITE(tmp,*) bool2 
    str=TRIM(ADJUSTL(tmp))
    IF (LEN(str1) /= 0) str = str1//str
    RETURN
  END FUNCTION fis_cat_bool

  FUNCTION fis_cat_bool_inv(bool2,str1) RESULT(str)
    !! Concatenate a string with a logical (reversed)
    LOGICAL, INTENT(in)           :: bool2 !! Logical to concatenate
    CHARACTER(len=*), INTENT(in)  :: str1    !! String to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str     !! Output string
    CHARACTER(len=2) ::tmp 
    WRITE(tmp,*) bool2 
    str = TRIM(ADJUSTL(tmp))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_bool_inv

  FUNCTION fis_cat_real(str1,real2) RESULT(str)
    !! Concatenate a string with a real simple precision
    CHARACTER(len=*), INTENT(in)  :: str1    !! String to concatenate
    REAL(kind=4), INTENT(in)      :: real2 !! Simple precision real to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str     !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(real2)) :: str) 
    WRITE(str,*) real2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str=str1//str 
    RETURN
  END FUNCTION fis_cat_real

  FUNCTION fis_cat_real_inv(real2,str1) RESULT(str)
    !! Concatenate a string with a real simple precision (reversed)
    REAL(kind=4), INTENT(in)      :: real2 !! Simple precision real to concatenate
    CHARACTER(len=*), INTENT(in)  :: str1  !! String to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str   !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(real2)) :: str) 
    WRITE(str,*) real2  ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_real_inv

  FUNCTION fis_cat_double(str1,double2) RESULT(str)
    !! Concatenate a string with a real double precision
    CHARACTER(len=*), INTENT(in)  :: str1    !! String to concatenate
    REAL(kind=8), INTENT(in)      :: double2 !! Double precision real to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str     !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(double2)) :: str) 
    WRITE(str,*) double2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str=str1//str 
    RETURN
  END FUNCTION fis_cat_double

  FUNCTION fis_cat_double_inv(double2,str1) RESULT(str)
    !! Concatenate a string with a real double precision (reversed)
    REAL(kind=8), INTENT(in)      :: double2 !! Double precision real to concatenate
    CHARACTER(len=*), INTENT(in)  :: str1    !! String to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str     !! Output string
    ALLOCATE(CHARACTER(len=DIGITS(double2)) :: str) 
    WRITE(str,*) double2  ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_double_inv

  FUNCTION fis_cat_cplx(str1,cplx2) RESULT(str)
    !! Concatenate a string with a complex 
    CHARACTER(len=*), INTENT(in)  :: str1  !! String to concatenate 
    COMPLEX(kind=4), INTENT(in)   :: cplx2 !! Complex value to concatenate 
    CHARACTER(len=:), ALLOCATABLE :: str   !! Output string
    INTEGER :: sl
    sl = DIGITS(REAL(cplx2))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str,*) cplx2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_cplx

  FUNCTION fis_cat_cplx_inv(cplx2,str1) RESULT(str)
    !! Concatenate a string with a complex (reversed)
    COMPLEX(kind=4), INTENT(in)   :: cplx2 !! Complex value to concatenate 
    CHARACTER(len=*), INTENT(in)  :: str1  !! String to concatenate 
    CHARACTER(len=:), ALLOCATABLE :: str   !! Output string
    INTEGER :: sl
    sl = DIGITS(REAL(cplx2))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str,*) cplx2
    str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_cplx_inv

  FUNCTION fis_cat_dcplx(str1,dcplx2) RESULT(str)
    !! Concatenate a string with a double precision complex 
    CHARACTER(len=*), INTENT(in)  :: str1   !! String to concatenate
    COMPLEX(kind=8), INTENT(in)   :: dcplx2 !! Complex value to concatenate 
    CHARACTER(len=:), ALLOCATABLE :: str    !! Output string
    INTEGER :: sl
    sl = DIGITS(REAL(dcplx2))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str,*) dcplx2 ; str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_dcplx

  FUNCTION fis_cat_dcplx_inv(dcplx2,str1) RESULT(str)
    !! Concatenate a string with a double precision complex (reversed)
    COMPLEX(kind=8), INTENT(in)   :: dcplx2 !! Complex value to concatenate 
    CHARACTER(len=*), INTENT(in)  :: str1   !! string to concatenate
    CHARACTER(len=:), ALLOCATABLE :: str    !! Output string
    INTEGER :: sl
    sl = DIGITS(REAL(dcplx2))*2+3
    ALLOCATE(CHARACTER(len=sl) :: str)
    WRITE(str,*) dcplx2
    str = TRIM(ADJUSTL(str))
    IF (LEN(str1) /= 0) str = str//str1
    RETURN
  END FUNCTION fis_cat_dcplx_inv

  SUBROUTINE fis_affect_int(str,int)
    !! Assignment subroutine (using intrinsic integer) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str !! Output string to be assigned
    INTEGER, INTENT(in)                        :: int !! Input value to assign 
    str = fis_cat_int('',int)
  END SUBROUTINE fis_affect_int

  SUBROUTINE fis_affect_bool(str,bool)
    !! Assignment subroutine (using intrinsic logical) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str  !! Output string to be assigned
    LOGICAL, INTENT(in)                        :: bool !! Input value to assign
    str = fis_cat_bool('',bool)
  END SUBROUTINE fis_affect_bool

  SUBROUTINE fis_affect_real(str,float)
    !! Assignment subroutine (using intrinsic real) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str   !! Output string to be assigned
    REAL(kind=4), INTENT(in)                   :: float !! Input value to assign
    str = fis_cat_real('',float)
  END SUBROUTINE fis_affect_real

  SUBROUTINE fis_affect_double(str,double)
    !! Assignment subroutine (using intrinsic real(kind=8)) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str    !! Output string to be assigned
    REAL(kind=8), INTENT(in)                   :: double !! Input value to assign 
    str = fis_cat_double('',double)
  END SUBROUTINE fis_affect_double

  SUBROUTINE fis_affect_cplx(str,cplx)
    !! Assignment subroutine (using intrinsic complex) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str  !! Output string to be assigned
    COMPLEX(kind=4), INTENT(in)                :: cplx !! Input value to assign
    str = fis_cat_cplx('',cplx)
  END SUBROUTINE fis_affect_cplx

  SUBROUTINE fis_affect_dcplx(str,dcplx)
    !! Assignment subroutine (using intrinsic complex(kind=8)) 
    CHARACTER(len=:), INTENT(out), ALLOCATABLE :: str   !! Output string to be assigned
    COMPLEX(kind=8), INTENT(in)                :: dcplx !! Input value to assign
    str = fis_cat_dcplx('',dcplx)
  END SUBROUTINE fis_affect_dcplx

  FUNCTION get_attrs_indexes(flags) RESULT(codes)
    !! Convert a list of csi flags into a csi codes.
    !!
    !! Only know CSI codes are returned. If no known CSI are found the outputput vector is
    !! allocated with 0 elements.
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags !! CSI attributes flags 
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes         !! CSI codes.
    INTEGER :: i,j,n
    n = 0
    ALLOCATE(codes(SIZE(flags)))
    codes(:) = -1
    DO i = 1, SIZE(flags)
      DO j = 1, SIZE(attributes)
        IF (to_lower(flags(i)) == csis(j)) THEN
          n = n + 1
          codes(n) = attributes(j)
          EXIT
        ENDIF
      ENDDO
    ENDDO
    IF (n > 0) THEN
      codes = codes(1:n)
    ELSE
      DEALLOCATE(codes)
      ALLOCATE(codes(0))
    ENDIF
  END FUNCTION get_attrs_indexes

  FUNCTION fancy_fstr(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given (fortran intrinsic) string.
    CHARACTER(len=*), INTENT(in)               :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format (unused for this overload)
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    codes = get_attrs_indexes(flags) 
    IF (SIZE(codes) == 0) THEN
      output = value ; RETURN
    ELSE
      output = add_csi(value,codes)
    ENDIF
  END FUNCTION fancy_fstr

  FUNCTION fancy_int(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given integer value.
    INTEGER, INTENT(in)                        :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_int

  FUNCTION fancy_bool(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given logical value.
    LOGICAL, INTENT(in)                        :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_bool

  FUNCTION fancy_real(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given real value (simple precision).
    REAL(kind=4), INTENT(in)                   :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_real 

  FUNCTION fancy_double(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given real value (double precision).
    REAL(kind=8), INTENT(in)                   :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_double

  FUNCTION fancy_cplx(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given complex value (simple precision).
    COMPLEX(kind=4), INTENT(in)                :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_cplx
  
  FUNCTION fancy_dcplx(value,flags,fmt) RESULT(output)
    !! Compute a fancy string from the given complex value (double precision).
    COMPLEX(kind=8), INTENT(in)                :: value  !! String object reference 
    CHARACTER(len=2), DIMENSION(:), INTENT(in) :: flags  !! CSI attributes flags
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: fmt    !! Optional format. If given it must be a valid Fortran format.
    CHARACTER(len=:), ALLOCATABLE              :: output !! Output fortran instrinsic string
    INTEGER, DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=:), ALLOCATABLE      :: tmp 
    codes = get_attrs_indexes(flags) 
    IF (PRESENT(fmt)) THEN ; tmp = to_string(value,fmt) ; ELSE ; tmp = to_string(value) ; ENDIF
    IF (SIZE(codes) /= 0) THEN ; output = add_csi(tmp,codes) ; ELSE ; output = tmp ; ENDIF
  END FUNCTION fancy_dcplx

END MODULE STRING_OP
