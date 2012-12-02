      subroutine ch_cap ( ch )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if

      return
      end
      function ch_eqi ( c1, c2 )

c*********************************************************************72
c
cc CH_EQI is a case insensitive comparison of two characters for equality.
c
c  Example:
c
c    CH_EQI ( 'A', 'a' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C1, C2, the characters to compare.
c
c    Output, logical CH_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap == c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
      end if

      return
      end
      function ch_index_last ( s, ch )

c*********************************************************************72
c
cc CH_INDEX_LAST is the last occurrence of a character in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be searched.
c
c    Input, character CH, the character to be searched for.
c
c    Output, integer CH_INDEX_LAST, the location of the last 
c    occurrence of the character in the string, or -1 if it does not occur.
c
      implicit none

      character ch
      integer ch_index_last
      integer i
      character * ( * )  s
      integer s_len_trim
      integer s_length

      ch_index_last = -1
      s_length = s_len_trim ( s )

      do i = s_length, 1, -1

        if ( s(i:i) .eq. ch ) then
          ch_index_last = i
          return
        end if

      end do
 
      return
      end
      function ch_is_digit ( c )

c*********************************************************************72
c
cc CH_IS_DIGIT returns TRUE if a character is a decimal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the character to be analyzed.
c
c    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
c
      implicit none

      character c
      logical ch_is_digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
      else
        ch_is_digit = .false.
      end if

      return
      end
      subroutine ch_low ( ch )

!*********************************************************************72
!
!! CH_LOW lowercases a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to be lowercased.
!
      implicit none

      character ch
      integer i

      i = ichar ( ch )

      if ( 65 .le. i .and. i .le. 90 ) then
        ch = char ( i + 32 )
      end if
 
      return
      end
      subroutine ch_to_digit ( c, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT returns the integer value of a base 10 digit.
c
c  Example:
c
c     C   DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    ...  ...
c    '9'    9
c    ' '    0
c    'X'   -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the decimal digit, '0' through '9' or blank
c    are legal.
c
c    Output, integer DIGIT, the corresponding integer value.  If C was
c    'illegal', then DIGIT is -1.
c
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      function ch_to_rot13 ( ch )

c*********************************************************************72
c
cc CH_TO_ROT13 converts a character to its ROT13 equivalent.
c
c  Discussion:
c
c    Two applications of CH_TO_ROT13 to a character will return the original.
c
c    As a further scrambling, digits are similarly rotated using
c    a "ROT5" scheme.
c
c  Example:
c
c    Input:  Output:
c
c    a       n
c    C       P
c    J       W
c    1       6
c    5       0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, the character to be converted.
c
c    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
c
      implicit none

      character ch
      character ch_to_rot13
      integer itemp

      itemp = ichar ( ch )
c
c  [0:4] -> [5:9]
c
      if ( 48 .le. itemp .and. itemp .le. 52 ) then
        itemp = itemp + 5
c
c  [5:9] -> [0:4]
c
      else if ( 53 .le. itemp .and. itemp .le. 57 ) then
        itemp = itemp - 5
c
c  [A:M] -> [N:Z]
c
      else if ( 65 .le. itemp .and. itemp .le. 77 ) then
        itemp = itemp + 13
c
c  [N:Z] -> [A:M]
c
      else if ( 78 .le. itemp .and. itemp .le. 90 ) then
        itemp = itemp - 13
c
c  [a:m] -> [n:z]
c
      else if ( 97 .le. itemp .and. itemp .le. 109 ) then
        itemp = itemp + 13
c
c  [n:z] -> [a:m]
c
      else if ( 110 .le. itemp .and. itemp .le. 122 ) then
        itemp = itemp - 13
      end if
 
      ch_to_rot13 = char ( itemp )
 
      return
      end
      subroutine digit_inc ( c )

c*********************************************************************72
c
cc DIGIT_INC increments a decimal digit.
c
c  Example:
c
c    Input  Output
c    -----  ------
c    '0'    '1'
c    '1'    '2'
c    ...
c    '8'    '9'
c    '9'    '0'
c    'A'    'A'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character C, a digit to be incremented.
c
      implicit none

      character c
      integer digit

      call ch_to_digit ( c, digit )

      if ( digit .eq. -1 ) then
        return
      end if

      digit = digit + 1

      if ( digit .eq. 10 ) then
        digit = 0
      end if

      call digit_to_ch ( digit, c )

      return
      end
      subroutine digit_to_ch ( digit, c )

c*********************************************************************72
c
cc DIGIT_TO_CH returns the character representation of a decimal digit.
c
c  Example:
c
c    DIGIT   C
c    -----  ---
c      0    '0'
c      1    '1'
c    ...    ...
c      9    '9'
c     17    '*'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit value between 0 and 9.
c
c    Output, character C, the corresponding character, or '*' if DIGIT
c    was illegal.
c
      implicit none

      character c
      integer digit

      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine file_column_count ( input_file_name, column_num )

c*********************************************************************72
c
cc FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
c
c  Discussion:
c
c    The file is assumed to be a simple text file.
c
c    Most lines of the file is presumed to consist of COLUMN_NUM words,
c    separated by spaces.  There may also be some blank lines, and some
c    comment lines,
c    which have a "#" in column 1.
c
c    The routine tries to find the first non-comment non-blank line and
c    counts the number of words in that line.
c
c    If all lines are blanks or comments, it goes back and tries to analyze
c    a comment line.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the file.
c
c    Output, integer COLUMN_NUM, the number of columns in the file.
c
      implicit none

      integer column_num
      logical got_one
      character * ( * ) input_file_name
      integer input_unit
      character * ( 255 ) line
      integer s_len_trim
c
c  Open the file.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
     &  status = 'old',  form = 'formatted', access = 'sequential' )
c
c  Read one line, but skip blank lines and comment lines.
c
      got_one = .false.

10    continue

        read ( input_unit, '(a)', err = 20 ) line

        if ( s_len_trim ( line ) .eq. 0 ) then
          go to 10
        end if

        if ( line(1:1) .eq. '#' ) then
          go to 10
        end if

        got_one = .true.
        go to 20

      go to 10

20    continue

      if ( .not. got_one ) then

        rewind ( input_unit )

30      continue

          read ( input_unit, '(a)', err = 40 ) line

          if ( s_len_trim ( line ) .eq. 0 ) then
            go to 30
          end if

          got_one = .true.
          go to 40

        go to 30

40    continue

      end if

      close ( unit = input_unit )

      if ( .not. got_one ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning.'
        write ( *, '(a)' ) '  The file does not contain any data.'
        column_num = -1
        return
      end if

      call s_word_count ( line, column_num )

      return
      end
      subroutine file_delete ( file_name )

c*********************************************************************72
c
cc FILE_DELETE deletes a file if it exists.
c
c  Discussion:
c
c    You might want to call this routine to get rid of any old copy
c    of a file, before trying to open a new copy with the OPEN argument:
c
c      status = 'new'.
c
c    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
c
c    For instance, on the SGI, the most recent version of the FORTRAN
c    compiler seems to go crazy when I open an unformatted direct
c    access file this way.  It creates an enormous file (of somewhat
c    random size).  The problem goes away if I delete any old copy
c    using this routine, and then open a fresh copy with
c    " STATUS = 'NEW' ".  It's a scary world.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) FILE_NAME, the name of the file to be deleted.
c
      implicit none

      character*(*) file_name
      logical lfile
      integer s_len_trim
      integer unit
c
c  Does the file exist?
c
      inquire (
     &  file = file_name,
     &  exist = lfile )

      if ( .not. lfile ) then
        return
      end if
c
c  Can we get a FORTRAN unit number?
c
      call get_unit ( unit )

      if ( unit .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_DELETE - Error!'
        write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
        return
      end if
c
c  Can we open the file?
c  
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE:'
      write ( *, '(a)' ) '  Open "' // trim ( file_name ) // '".'

      open (
     &  unit = unit,
     &  file = file_name,
     &  status = 'old',
     &  err = 10 )
c
c  Can we close the file with "Delete" status?
c
      write ( *, '(a)' ) '  Delete "' // trim ( file_name ) // '".'

      close (
     &  unit = unit,
     &  status = 'delete',
     &  err = 20 )

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE - Error!'
      write ( *, '(a)' ) 
     &  '  Could not open "' // trim ( file_name ) // '".'
      return

20    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE - Error!'
      write ( *, '(a)' ) 
     &  '  Could not delete "' // trim ( file_name ) // '".'

      return
      end
      function file_exist ( file_name )

c*****************************************************************************80
c
cc FILE_EXIST reports whether a file exists.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) FILE_NAME, the name of the file.
c
c    Output, logical FILE_EXIST, is TRUE if the file exists.
c
      implicit none

      logical file_exist
      character*(*) file_name

      inquire ( file = file_name, exist = file_exist )

      return
      end
      subroutine file_name_ext_get ( file_name, i, j )

c*********************************************************************72
c
cc FILE_NAME_EXT_GET determines the "extension" of a file name.
c
c  Discussion:
c
c    The "extension" of a filename is the string of characters
c    that appears after the LAST period in the name.  A file
c    with no period, or with a period as the last character
c    in the name, has a "null" extension.
c
c    Blanks are unusual in filenames.  This routine ignores all
c    trailing blanks, but will treat initial or internal blanks
c    as regular characters acceptable in a file name.
c
c  Example:
c
c    FILE_NAME   I  J
c
c    bob.for     4  7
c    N.B.C.D     6  7
c    Naomi.      6  6
c    Arthur     -1 -1 
c    .com        1  1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_NAME, a file name to be examined.
c
c    Output, integer I, J, the indices of the first and 
c    last characters in the file extension.  
c    If no period occurs in FILE_NAME, then
c      I = J = -1;
c    Otherwise,
c      I is the position of the LAST period in FILE_NAME, and J is the
c      position of the last nonblank character following the period.
c
      implicit none

      character * ( * ) file_name
      integer i
      integer j
      integer ch_index_last

      i = ch_index_last ( file_name, '.' )

      if ( i .eq. -1 ) then

        j = -1

      else

        j = len_trim ( file_name )

      end if

      return
      end
      subroutine file_name_ext_swap ( file_name, ext )

c*********************************************************************72
c
cc FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
c
c  Discussion:
c
c    The "extension" of a filename is the string of characters
c    that appears after the LAST period in the name.  A file
c    with no period, or with a period as the last character
c    in the name, has a "null" extension.
c
c  Example:
c
c          Input           Output
c    ================     =========
c    FILE_NAME    EXT     FILE_NAME
c
c    bob.for      obj     bob.obj
c    bob.bob.bob  txt     bob.bob.txt
c    bob          yak     bob.yak
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) FILE_NAME, a file name.
c    On output, the extension of the file has been changed.
c
c    Input, character * ( * ) EXT, the extension to be used on the output
c    copy of FILE_NAME, replacing the current extension if any.
c
      implicit none

      character * ( * ) ext
      character * ( * ) file_name
      integer i
      integer j
      integer len_max
      integer len_name

      len_max = len ( file_name )
      len_name = len_trim ( file_name )

      call file_name_ext_get ( file_name, i, j )

      if ( i .eq. 0 ) then

        if ( len_max .lt. len_name + 1 ) then
          return
        end if

        len_name = len_name + 1
        file_name(len_name:len_name) = '.'
        i = len_name + 1

      else

        i = i + 1
        file_name(i:j) = ' '

      end if

      file_name(i:) = ext

      return
      end
      subroutine filename_inc ( filename )

c*********************************************************************72
c
cc FILENAME_INC increments a partially numeric filename.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  Non-numeric letters of the name are unaffected.
c
c    If the name is empty, then the routine stops.
c
c    If the name contains no digits, the empty string is returned.
c
c  Example:
c
c      Input          Output
c      -----          ------
c      a7to11.txt     a7to12.txt
c      a7to99.txt     a8to00.txt
c      a9to99.txt     a0to00.txt
c      cat.txt        ' '
c      ' '            STOP!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) FILENAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      integer change
      integer digit
      character*(*) filename
      integer i
      integer lens

      lens = len_trim ( filename )

      if ( lens .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
        write ( *, '(a)' ) '  The input string is empty.'
        stop
      end if

      change = 0

      do i = lens, 1, -1

        c = filename(i:i)

        if ( ch_is_digit ( c ) ) then

          change = change + 1

          call digit_inc ( c )

          filename(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do
c
c  No digits were found.  Return blank.
c
      if ( change .eq. 0 ) then
        filename = ' '
        return
      end if

      return
      end
      subroutine filename_inc_nowrap ( filename )

c*********************************************************************72
c
cc FILENAME_INC_NOWRAP increments a partially numeric filename.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  Non-numeric letters of the name are unaffected.
c
c    If the (nonempty) name contains no digits, or all the digits are
c    9, then the empty string is returned.
c
c    If the empty string is input, the routine stops.
c
c  Example:
c
c      Input            Output
c      -----            ------
c      'a7to11.txt'     'a7to12.txt'
c      'a7to99.txt'     'a8to00.txt'
c      'a8to99.txt'     'a9to00.txt'
c      'a9to99.txt'     ' '
c      'cat.txt'        ' '
c      ' '              STOP!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) FILENAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      integer carry
      logical ch_is_digit
      integer change
      integer digit
      character*(*) filename
      integer i
      integer lens

      lens = len_trim ( filename )

      if ( lens .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
        write ( *, '(a)' ) '  The input string is empty.'
        stop
      end if

      change = 0
      carry = 0

      do i = lens, 1, -1

        c = filename(i:i)

        if ( ch_is_digit ( c ) ) then

          change = change + 1
          carry = 0

          call digit_inc ( c )

          filename(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

          carry = 1

        end if

      end do
c
c  Unsatisfied carry.  The input digits were all 9.  Return blank.
c
      if ( carry == 1 ) then
        filename = ' '
        return
      end if
c
c  No digits were found.  Return blank.
c
      if ( change .eq. 0 ) then
        filename = ' '
        return
      end if

      return
      end
      subroutine file_row_count ( input_file_name, row_num )

c*********************************************************************72
c
cc FILE_ROW_COUNT counts the number of row records in a file.
c
c  Discussion:
c
c    It does not count lines that are blank, or that begin with a
c    comment symbol '#'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Output, integer ROW_NUM, the number of rows found.
c
      implicit none

      integer bad_num
      integer comment_num
      integer ierror
      character * ( * ) input_file_name
      integer input_status
      integer input_unit
      character * ( 255 ) line
      integer record_num
      integer row_num
      integer s_len_trim

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
     &  status = 'old' )

      comment_num = 0
      row_num = 0
      record_num = 0
      bad_num = 0

10    continue

        read ( input_unit, '(a)', err = 20, end = 20 ) line

        record_num = record_num + 1

        if ( line(1:1) .eq. '#' ) then
          comment_num = comment_num + 1
          go to 10
        end if

        if ( s_len_trim ( line ) .eq. 0 ) then
          comment_num = comment_num + 1
          go to 10
        end if

        row_num = row_num + 1

      go to 10

20    continue

      close ( unit = input_unit )

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

      return
      end
      function s_index_last ( s, sub )

c*********************************************************************72
c
cc S_INDEX_LAST finds the LAST occurrence of a given substring.
c
c  Discussion:
c
c    It returns the location in the string at which the substring SUB is
c    first found, or 0 if the substring does not occur at all.
c
c    The routine is also trailing blank insensitive.  This is very
c    important for those cases where you have stored information in
c    larger variables.  If S is of length 80, and SUB is of
c    length 80, then if S = 'FRED' and SUB = 'RED', a match would
c    not be reported by the standard FORTRAN INDEX, because it treats
c    both variables as being 80 characters long!  This routine assumes that
c    trailing blanks represent garbage!
c
c    This means that this routine cannot be used to find, say, the last
c    occurrence of a substring 'A ', since it assumes the blank space
c    was not specified by the user, but is, rather, padding by the
c    system.  However, as a special case, this routine can properly handle
c    the case where either S or SUB is all blanks.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be searched.
c
c    Input, character * ( * ) SUB, the substring to search for.
c
c    Output, integer S_INDEX_LAST.  0 if SUB does not occur in
c    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
c    where LENS is the length of SUB, and is the last place
c    this happens.
c
      implicit none

      integer i
      integer j
      integer llen2
      character * ( * ) s
      integer s_index_last
      integer s_length
      character * ( * ) sub

      s_index_last = 0

      s_length = len_trim ( s )
      llen2 = len_trim ( sub )
c
c  In case S or SUB is blanks, use LEN.
c
      if ( s_length .eq. 0 ) then
        s_length = len ( s )
      end if

      if ( llen2 .eq. 0 ) then
        llen2 = len ( sub )
      end if

      if ( s_length .lt. llen2 ) then
        return
      end if

      do j = 1, s_length + 1 - llen2

        i = s_length + 2 - llen2 - j

        if ( s(i:i+llen2-1) .eq. sub ) then
          s_index_last = i
          return
        end if

      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine s_low ( s )

c*********************************************************************72
c
cc S_LOW replaces all uppercase letters by lowercase ones.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to be
c    transformed.  On output, the string is all lowercase.
c
      implicit none

      integer i
      character * ( * ) s
      integer s_length

      s_length = len_trim ( s )
 
      do i = 1, s_length
        call ch_low ( s(i:i) )
      end do
 
      return
      end
      subroutine s_reverse ( s )

c*********************************************************************72
c
cc S_REVERSE reverses the characters in a string.
c
c  Example:
c
c    Input        Output
c
c    ' Cat'       'taC '
c    'Goo gol  '  'log ooG  '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to reverse.
c    Trailing blanks are ignored.
c
      implicit none

      character ch
      integer i
      integer j
      character * ( * ) s
      integer s_len_trim
      integer s_length

      s_length = s_len_trim ( s )
 
      do i = 1, s_length / 2
        j = s_length + 1 - i
        ch     = s(i:i)
        s(i:i) = s(j:j)
        s(j:j) = ch
      end do
 
      return
      end
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

      return
      end
      subroutine s_to_r8 ( s, dval, ierror, length )

c*********************************************************************72
c
cc S_TO_R8 reads an R8 from a string.
c
c  Discussion:
c
c    The routine will read as many characters as possible until it reaches
c    the end of the string, or encounters a character which cannot be
c    part of the number.
c
c    Legal input is:
c
c       1 blanks,
c       2 '+' or '-' sign,
c       2.5 blanks
c       3 integer part,
c       4 decimal point,
c       5 fraction part,
c       6 'E' or 'e' or 'D' or 'd', exponent marker,
c       7 exponent sign,
c       8 exponent integer part,
c       9 exponent decimal point,
c      10 exponent fraction part,
c      11 blanks,
c      12 final comma or semicolon,
c
c    with most quantities optional.
c
c  Example:
c
c    S                 DVAL
c
c    '1'               1.0
c    '     1   '       1.0
c    '1A'              1.0
c    '12,34,56'        12.0
c    '  34 7'          34.0
c    '-1E2ABCD'        -100.0
c    '-1X2ABCD'        -1.0
c    ' 2E-1'           0.2
c    '23.45'           23.45
c    '-4.2E+2'         -420.0
c    '17d2'            1700.0
c    '-14e-2'         -0.14
c    'e2'              100.0
c    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate at the end of the string, or when no more
c    characters can be read to form a legal real.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, double precision DVAL, the value read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    1, 2, 6 or 7, the input number was garbled.  The
c    value of IERROR is the last type of input successfully
c    read.  For instance, 1 means initial blanks, 2 means
c    a plus or minus sign, and so on.
c
c    Output, integer LENGTH, the number of characters read
c    to form the number, including any terminating
c    characters such as a trailing comma or blanks.
c
      implicit none

      logical ch_eqi
      character c
      double precision dval
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer length
      integer nchar
      integer ndig
      double precision rbot
      double precision rexp
      double precision rtop
      character * ( * ) s
      integer s_len_trim

      nchar = s_len_trim ( s )

      ierror = 0
      dval = 0.0D+00
      length = -1
      isgn = 1
      rtop = 0
      rbot = 1
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

10    continue

        length = length + 1

        if ( nchar .lt. length+1 ) then
          go to 20
        end if

        c = s(length+1:length+1)
c
c  Blank character.
c
        if ( c .eq. ' ' ) then

          if ( ihave .eq. 2 ) then

          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            iterm = 1
          else if ( 1 .lt. ihave ) then
            ihave = 11
          end if
c
c  Comma.
c
        else if ( c .eq. ',' .or. c .eq. ';' ) then

          if ( ihave .ne. 1 ) then
            iterm = 1
            ihave = 12
            length = length + 1
          end if
c
c  Minus sign.
c
        else if ( c .eq. '-' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
            isgn = -1
          else if ( ihave .eq. 6 ) then
            ihave = 7
            jsgn = -1
          else
            iterm = 1
          end if
c
c  Plus sign.
c
        else if ( c .eq. '+' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
          else if ( ihave .eq. 6 ) then
            ihave = 7
          else
            iterm = 1
          end if
c
c  Decimal point.
c
        else if ( c .eq. '.' ) then

          if ( ihave .lt. 4 ) then
            ihave = 4
          else if ( 6 .le. ihave .and. ihave .le. 8 ) then
            ihave = 9
          else
            iterm = 1
          end if
c
c  Scientific notation exponent marker.
c
        else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

          if ( ihave .lt. 6 ) then
            ihave = 6
          else
            iterm = 1
          end if
c
c  Digit.
c
        else if ( ihave .lt. 11 .and. lle ( '0', c ) 
     &    .and. lle ( c, '9' ) ) then

          if ( ihave .le. 2 ) then
            ihave = 3
          else if ( ihave .eq. 4 ) then
            ihave = 5
          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            ihave = 8
          else if ( ihave .eq. 9 ) then
            ihave = 10
          end if

          call ch_to_digit ( c, ndig )

          if ( ihave .eq. 3 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
          else if ( ihave .eq. 5 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
            rbot = 10.0D+00 * rbot
          else if ( ihave .eq. 8 ) then
            jtop = 10 * jtop + ndig
          else if ( ihave .eq. 10 ) then
            jtop = 10 * jtop + ndig
            jbot = 10 * jbot
          end if
c
c  Anything else is regarded as a terminator.
c
        else
          iterm = 1
        end if
c
c  If we haven't seen a terminator, and we haven't examined the
c  entire string, go get the next character.
c
        if ( iterm .eq. 1 ) then
          go to 20
        end if

        go to 10

20    continue
c
c  If we haven't seen a terminator, and we have examined the
c  entire string, then we're done, and LENGTH is equal to NCHAR.
c
      if ( iterm .ne. 1 .and. length+1 .eq. nchar ) then
        length = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c  Not if we terminated in states 1, 2, 6 or 7.
c
      if ( ihave .eq. 1 .or. ihave .eq. 2 .or. 
     &     ihave .eq. 6 .or. ihave .eq. 7 ) then
        ierror = ihave
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a,a)' ) '    ', s
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0D+00
      else
        if ( jbot .eq. 1 ) then
          rexp = 10.0D+00 ** ( jsgn * jtop )
        else
          rexp = 10.0D+00 ** ( dble ( jsgn * jtop ) / dble ( jbot ) )
        end if
      end if

      dval = dble ( isgn ) * rexp * rtop / rbot

      return
      end
      subroutine s_to_r8vec ( s, n, rvec, ierror )

c*********************************************************************72
c
cc S_TO_R8VEC reads an R8VEC from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be read.
c
c    Input, integer N, the number of values expected.
c
c    Output, double precision RVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer  n

      integer i
      integer ierror
      integer ilo
      integer lchar
      double precision rvec(n)
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + lchar

        go to 10

      end if

20    continue

      return
      end
      subroutine s_to_rot13 ( s )

c*********************************************************************72
c
cc S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
c
c  Discussion:
c
c    Two applications of the routine will return the original string.
c
c  Example:
c
c    Input:                      Output:
c
c    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
c    Cher                        Pure
c    James Thurston Howell       Wnzrf Guhefgba Ubjryy
c    0123456789                  5678901234
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, a string to be "rotated".
c
      implicit none

      character ch_to_rot13
      integer i
      character * ( * ) s
      integer s_length

      s_length = len_trim ( s )
 
      do i = 1, s_length
        s(i:i) = ch_to_rot13 ( s(i:i) )
      end do
 
      return
      end
      subroutine s_word_count ( s, nword )

c*********************************************************************72
c
cc S_WORD_COUNT counts the number of "words" in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be examined.
c
c    Output, integer NWORD, the number of "words" in the string.
c    Words are presumed to be separated by one or more blanks.
c
      implicit none

      logical blank
      integer i
      integer lens
      integer nword
      character * ( * ) s

      nword = 0
      lens = len ( s )

      if ( lens .le. 0 ) then
        return
      end if

      blank = .true.

      do i = 1, lens

        if ( s(i:i) .eq. ' ' ) then
          blank = .true.
        else if ( blank ) then
          nword = nword + 1
          blank = .false.
        end if

      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
