!  spacer_data_convert.f  03 May 2000
!
      program spacer_data_convert
!
!***********************************************************************
!
!! SPACER_DATA_CONVERT converts a data file to SPACER format.
!
!
!  Discussion:
!
!    The input data file has a fairly simple form.  
!
      integer, parameter :: maxn = 100
!
      integer ierror
      character ( len = 80 ) input_filename
      character ( len = 10 ) name(maxn)
      integer ncol
      integer nrow
      character ( len = 80 ) output_filename
      real xx(maxn,maxn)
!
      write ( *, * ) ' '
      write ( *, * ) 'SPACER_DATA_CONVERT'
      write ( *, * ) '  Read Hugh''s data, '
      write ( *, * ) '  Print it out,'
      write ( *, * ) '  Write it to a file for SPACER.'
      write ( *, * ) ' '
!
!  Read Hugh's data.
!
      input_filename = 'input_data.txt'

      call input_read ( input_filename, maxn, nrow, ncol, xx, name, 
     &  ierror )

      if ( ierror /= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SPACER_DATA_CONVERT'
        write ( *, * ) '  Failure while reading the input data.'
        stop
      end if
!
!  Print the data.
!
      call input_print ( maxn, nrow, ncol, xx, name )
!
!  Write the data to a file.
!
      output_filename = 'spacer_data.txt'

      call output_write ( output_filename, maxn, nrow, ncol, xx,
     &  ierror )

      if ( ierror /= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SPACER_DATA_CONVERT'
        write ( *, * ) '  Failure while writing the output data.'
        stop
      end if

      write ( *, * ) ' '
      write ( *, * ) 'SPACER_DATA_CONVERT'
      write ( *, * ) '  Normal end of execution.'

      stop
      end
      subroutine c_cap ( c )
!
!***********************************************************************
!
!! C_CAP capitalizes a single character.
!
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
!    Input/output, character*1 C, the character to capitalize.
!
      character*1 c
      integer itemp
!
      itemp = ichar ( c )

      if ( 97 <= itemp .and. itemp <= 122 ) then
        c = char ( itemp - 32 )
      end if

      return
      end
      function c_eqi ( c1, c2 )
!
!***********************************************************************
!
!! C_EQI is a case insensitive comparison of two characters for equality.
!
!
!  Examples:
!
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*1 C1, C2, the characters to compare.
!
!    Output, logical C_EQI, the result of the comparison.
!
      logical c_eqi
      character*1 c1
      character*1 c2
      character*1 cc1
      character*1 cc2
!
      cc1 = c1
      cc2 = c2

      call c_cap ( cc1 )
      call c_cap ( cc2 )

      if ( cc1 == cc2 ) then
        c_eqi = .true.
      else
        c_eqi = .false.
      end if

      return
      end
      subroutine c_to_digit ( c, digit )
!
!***********************************************************************
!
!! C_TO_DIGIT returns the integer value of a base 10 digit.
!
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*1 C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
      character*1 c
      integer digit
!
      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c == ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      subroutine file_delete ( file_name )
!
!***********************************************************************
!
!! FILE_DELETE deletes a named file if it exists.
!
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) FILE_NAME, the name of the file to be deleted.
!
      character*80 ctemp
      character*(*) file_name
      integer iunit
      integer lens
      logical lfile
!
!  Does the file exist?
!
      inquire (
     &  file = file_name,
     &  exist = lfile )

      if ( .not. lfile ) then
        return
      end if
!
!  Get a free unit number.
!
      call get_unit ( iunit )

      if ( iunit == 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'FILE_DELETE: Warning!'
        write ( *, * ) '  A free FORTRAN unit could not be found.'
        return
      end if

      write ( *, * ) ' '
      ctemp = 'FILE_DELETE: deleting old version of ' // file_name
      lens = len_trim ( ctemp )
      write ( *, '(a)' ) ctemp(1:lens)

      open (
     &  unit = iunit,
     &  file = file_name,
     &  status = 'old',
     &  err = 10 )

      close (
     &  unit = iunit,
     &  status = 'delete' )

      return

10    continue

      write ( *, * ) ' '
      write ( *, * ) 'FILE_DELETE: Warning!'
      write ( *, * ) '  Could not open the file.'

      return
      end
      subroutine get_unit ( iunit )
!
!***********************************************************************
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    If IUNIT < 0, then an I/O error occurred while trying to inquire
!    on the status of unit abs(IUNIT).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
      integer i
      integer iunit
      logical lopen

      iunit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 ) then

          iunit = -i

          inquire (
     &      unit = i,
     &      opened = lopen,
     &      err = 10 )

          if ( .not. lopen ) then
            iunit = i
            return
          end if

        end if

      end do
!
!  No free unit was found.
!
      iunit = 0

      return
!
!  An I/O error occurred during an INQUIRE.
!
10    continue

      return
      end
      subroutine i_extract ( s, i, ierror )
!
!***********************************************************************
!
!! I_EXTRACT "extracts" an integer from the beginning of a string.
!
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character*(*) S; on input, a string from
!    whose beginning an integer is to be extracted.  On output,
!    the integer, if found, has been removed.
!
!    Output, integer I.  If IERROR is 0, then I contains the
!    "next" integer read from S; otherwise I is 0.
!
!    Output, integer IERROR.
!    0, no error.
!    nonzero, an integer could not be extracted from the beginning of the
!    string.  I is 0 and S is unchanged.
!
      integer i
      integer ierror
      integer lchar
      character*(*) s
!
      i = 0

      call s_to_i ( s, i, ierror, lchar )

      if ( ierror /= 0 .or. lchar == 0 ) then
        ierror = 1
        i = 0
      else
        call s_shift_left ( s, lchar )
      end if

      return
      end
      subroutine input_print ( maxn, nrow, ncol, xx, name )
!
!***********************************************************************
!
!! INPUT_PRINT prints the input data.
!
!
!  Modified:
!
!    03 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer NROW, the number of rows of data.
!
!    Input, integer NCOL, the number of columns of data.
!
!    Input, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Input, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
      integer maxn
!
      integer i
      integer j
      integer jhi
      integer jlo
      character ( len = 10 ) name(maxn)
      integer ncol
      integer nrow
      real xx(maxn,maxn)
!
      write ( *, * ) ' '
      write ( *, * ) 'Object names:'
      write ( *, * ) ' '
      do i = 1, nrow
        write ( *, '(i6,2x,a10)' ) i, name(i)
      end do

      write ( *, * ) ' '
      write ( *, * ) 'Object distances:'
      write ( *, * ) ' '
      do i = 1, nrow
        do jlo = 1, ncol, 10
          jhi = min ( jlo + 9, ncol )
          if ( jlo == 1 ) then
            write ( *, '(i6,2x,10f8.3)' ) 
     &        i, ( xx(i,j), j = jlo, jhi )
          else
            write ( *, '(6x,2x,10f8.3)' )
     &           ( xx(i,j), j = jlo, jhi )
          end if
        end do
      end do

      return
      end
      subroutine input_read ( input_filename, maxn, nrow, ncol, xx, 
     &  name, ierror )
!
!***********************************************************************
!
!! INPUT_READ reads the distance data from a file.
!
!
!  Modified:
!
!    03 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Output, integer NROW, the number of rows of data.
!
!    Output, integer NCOL, the number of columns of data.
!
!    Output, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, character ( len = 10 ) NAME(MAXN), the names of the objects.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, an error occurred.
!
      integer maxn
!
      integer i
      integer ierror
      character ( len = 80 ) input_filename
      integer input_unit
      integer j
      character ( len = 10 ) name(maxn)
      integer ncol
      integer nrow
      character ( len = 200 ) s
      real xx(maxn,maxn)
!
      ierror = 0

      write ( *, * ) ' '
      write ( *, * ) 'INPUT_READ'
      write ( *, * ) '  Reading data from ' // trim ( input_filename )

      call get_unit ( input_unit )

      open (
     &  unit = input_unit,
     &  file = input_filename,
     &  status = 'old',
     &  err = 10 )

      read ( input_unit, '(a)', end = 10 ) s

      call i_extract ( s, nrow, ierror )
      ncol = nrow

      do i = 1, nrow

        read ( input_unit, '(a)', end = 10 ) s

        call word_extract ( s, name(i) )

        do j = 1, ncol

20        continue

          call r_extract ( s, xx(i,j), ierror )

          if ( ierror /= 0 ) then
            read ( input_unit, '(a)', end = 10 ) s
            go to 20
          end if

        end do

      end do

      close (
     &  unit = input_unit )

      return

10    continue

      ierror = 1

      write ( *, * ) ' '
      write ( *, * ) 'INPUT_READ'
      write ( *, * ) '  An error occurred!'

      return
      end
      subroutine output_write ( output_filename, maxn, nrow, ncol, xx,
     &  ierror )
!
!***********************************************************************
!
!! OUTPUT_WRITE writes the distance data to a file in SPACER format.
!
!
!  Modified:
!
!    03 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 80 ) OUTPUT_FILENAME, the name of the output file.
!
!    Input, integer MAXN, the maximum values for NROW and NCOL.
!
!    Input, integer NROW, the number of rows of data.
!
!    Input, integer NCOL, the number of columns of data.
!
!    Input, real XX(MAXN,MAXN), the NROW by NCOL distance data.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, an error occurred.
!
      integer maxn
!
      integer i
      integer ierror
      integer j
      integer ncol
      integer nrow
      character ( len = 80 ) output_filename
      integer output_unit
      real xx(maxn,maxn)
!
      ierror = 0

      write ( *, * ) ' '
      write ( *, * ) 'OUTPUT_WRITE'
      write ( *, * ) '  Writing data to ' // trim ( output_filename )

      call file_delete ( output_filename )

      call get_unit ( output_unit )

      open (
     &  unit = output_unit,
     &  file = output_filename,
     &  status = 'new',
     &  err = 10 )

      write ( output_unit, '(i6)' ) nrow

      do i = 1, nrow
        write ( output_unit, '(10f8.3)' ) ( xx(i,j), j = 1, i )
      end do

      close ( unit = output_unit )

      return

10    continue

      ierror = 1

      write ( *, * ) ' '
      write ( *, * ) 'OUTPUT_WRITE - Fatal error!'
      write ( *, * ) '  Could not open the output file.'
      return
      end
      subroutine r_extract ( s, r, ierror )
!
!***********************************************************************
!
!! R_EXTRACT "extracts" a real from the beginning of a string.
!
!
!  Modified:
!
!    02 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character*(*) S; on input, a string from
!    whose beginning a real is to be extracted.  On output,
!    the real, if found, has been removed.
!
!    Output, real R.  If IERROR is 0, then R contains the
!    next real read from the string; otherwise R is 0.
!
!    Output, integer IERROR.
!    0, no error.
!    nonzero, a real could not be extracted from the beginning of the
!    string.  R is 0.0 and S is unchanged.
!
      integer ierror
      integer lchar
      real r
      character*(*) s
!
      r = 0.0

      call s_to_r ( s, r, ierror, lchar )

      if ( ierror /= 0 .or. lchar == 0 ) then
        ierror = 1
        r = 0.0
      else
        call s_shift_left ( s, lchar )
      end if

      return
      end
      subroutine s_shift_left ( s, ishft )
!
!***********************************************************************
!
!! S_SHIFT_LEFT shifts the characters in a string to the left and blank pads.
!
!
!  Discussion:
!
!    A shift of 2 would change "Violin" to "olin  ".
!    A shift of -2 would change "Violin" to "  Violin".
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character*(*) S, the string to be shifted.
!
!    Input, integer ISHFT, the number of positions to the
!    left to shift the characters.
!
      integer i
      integer ishft
      integer nchar
      character*(*) s
!
      nchar = len ( s )

      if ( ishft > 0 ) then

        do i = 1, nchar - ishft
          s(i:i) = s(i+ishft:i+ishft)
        end do

        do i = nchar - ishft + 1, nchar
          s(i:i) = ' '
        end do

      else if ( ishft < 0 ) then

        do i = nchar, - ishft + 1, - 1
          s(i:i) = s(i+ishft:i+ishft)
        end do

        do i = - ishft, 1, -1
          s(i:i) = ' '
        end do

      end if

      return
      end
      subroutine s_to_i ( s, ival, ierror, last )
!
!***********************************************************************
!
!! S_TO_I reads an integer value from a string.
!
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character that was
!    part of the representation of IVAL.
!
      character*1 c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer last
      integer lens
      character*(*) s
!
      ierror = 0
      istate = 0

      isgn = 1
      ival = 0

      lens = len ( s )

      i = 0

10    continue

      i = i + 1

      c = s(i:i)

      if ( istate == 0 ) then

        if ( c == ' ' ) then

        else if ( c == '-' ) then
          istate = 1
          isgn = -1
        else if ( c == '+' ) then
          istate = 1
          isgn = + 1
        else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
          istate = 2
          ival = ichar ( c ) - ichar ( '0' )
        else
          ierror = 1
          return
        end if

      else if ( istate == 1 ) then

        if ( c == ' ' ) then

        else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
          istate = 2
          ival = ichar ( c ) - ichar ( '0' )
        else
          ierror = 1
          return
        end if

      else if ( istate == 2 ) then

        if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
          ival = 10 * ival + ichar ( c ) - ichar ( '0' )
        else
          istate = 3
        end if

      end if
!
!  Continue or exit?
!
      if ( istate == 3 ) then
        ival = isgn * ival
        last = i - 1
        return
      else if ( i >= lens ) then
        if ( istate == 2 ) then
          ival = isgn * ival
          last = lens
        else
          ierror = 1
          last = 0
        end if
        return
      end if

      go to 10
      end
      subroutine s_to_r ( s, r, ierror, lchar )
!
!***********************************************************************
!
!! S_TO_R reads a real number from a string.
!
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
      logical c_eqi
      character*1 c
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer lchar
      integer nchar
      integer ndig
      real r
      real rbot
      real rexp
      real rtop
      character*(*) s
      character ( len = 1 ), parameter :: TAB = char ( 9 )
!
      nchar = len ( s )
      ierror = 0
      r = 0.0
      lchar = - 1
      isgn = 1
      rtop = 0.0
      rbot = 1.0
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

10    continue

      lchar = lchar + 1
      c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
      if ( c == ' ' .or. c == TAB ) then
!
!  20 November 1993
!
!  I would like to allow input like "+ 2", where there is a space
!  between the plus and the number.  So I am going to comment out
!  this line, because I think that's all that's keeping me from
!  doing this.
!
!       if ( ihave == 2 .or.
!    &     ihave == 6 .or.
!    &     ihave == 7 ) then

        if ( ihave == 2 ) then

        else if ( ihave == 6 .or. ihave == 7 ) then
          iterm = 1
        else if ( ihave > 1 ) then
          ihave = 11
        end if
!
!  Comma.
!
      else if ( c == ',' .or. c == ';' ) then

        if ( ihave /= 1 ) then
          iterm = 1
          ihave = 12
          lchar = lchar + 1
        end if
!
!  Minus sign.
!
      else if ( c == '-' ) then

        if ( ihave == 1 ) then
          ihave = 2
          isgn = - 1
        else if ( ihave == 6 ) then
          ihave = 7
          jsgn = - 1
        else
          iterm = 1
        end if
!
!  Plus sign.
!
      else if ( c == '+' ) then

        if ( ihave == 1 ) then
          ihave = 2
        else if ( ihave == 6 ) then
          ihave = 7
        else
          iterm = 1
        end if
!
!  Decimal point.
!
      else if ( c == '.' ) then

        if ( ihave < 4 ) then
          ihave = 4
        else if ( ihave >= 6 .and. ihave <= 8 ) then
          ihave = 9
        else
          iterm = 1
        end if
!
!  Exponent marker.
!
      else if (
     &  c_eqi ( c, 'E' ) .or.
     &  c_eqi ( c, 'D' ) ) then

        if ( ihave < 6 ) then
          ihave = 6
        else
          iterm = 1
        end if
!
!  Digit.
!
      else if ( ihave < 11 .and.
     &  lge ( c, '0' ) .and.
     &  lle ( c, '9' ) ) then

        if ( ihave <= 2 ) then
          ihave = 3
        else if ( ihave == 4 ) then
          ihave = 5
        else if ( ihave == 6 .or. ihave == 7 ) then
          ihave = 8
        else if ( ihave == 9 ) then
          ihave = 10
        end if

        call c_to_digit ( c, ndig )

        if ( ihave == 3 ) then
          rtop = 10.0 * rtop + real ( ndig )
        else if ( ihave == 5 ) then
          rtop = 10.0 * rtop + real ( ndig )
          rbot = 10.0 * rbot
        else if ( ihave == 8 ) then
          jtop = 10 * jtop + ndig
        else if ( ihave == 10 ) then
          jtop = 10 * jtop + ndig
          jbot = 10 * jbot
        end if
!
!  Anything else is regarded as a terminator.
!
      else
        iterm = 1
      end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
      if ( iterm /= 1 .and. lchar+1 < nchar ) then
        go to 10
      end if
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
      if ( iterm /= 1 .and. lchar+1 == nchar ) then
        lchar = nchar
      end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
      if (
     &  ihave == 1 .or.
     &  ihave == 2 .or.
     &  ihave == 6 .or.
     &  ihave == 7 ) then

        ierror = ihave

        return
      end if
!
!  Number seems OK.  Form it.
!
      if ( jtop == 0 ) then
        rexp = 1.0
      else

        if ( jbot == 1 ) then
          rexp = 10.0**( jsgn * jtop )
        else
          rexp = jsgn * jtop
          rexp = rexp / jbot
          rexp = 10.0**rexp
        end if

      end if

      r = isgn * rexp * rtop / rbot

      return
      end
      subroutine word_extract ( s, word )
!
!***********************************************************************
!
!! WORD_EXTRACT extracts the next word from a string.
!
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character*(*) S, the string.  On output, the first
!    word has been removed, and the string has been shifted left.
!
!    Output, character*(*) WORD, the leading word of the string.
!
      integer iget1
      integer iget2
      integer lchar
      character*(*) s
      character*(*) word
!
      word = ' '

      lchar = len ( s )
!
!  Find the first nonblank.
!
      iget1 = 0

10    continue

      iget1 = iget1 + 1

      if ( iget1 > lchar ) then
        return
      end if

      if ( s(iget1:iget1) == ' ' ) then
        go to 10
      end if
!
!  Now look for the last contiguous nonblank
!
      iget2 = iget1

20    continue

      if ( iget2 < lchar ) then
        if ( s(iget2+1:iget2+1) /= ' ' ) then
          iget2 = iget2 + 1
          go to 20
        end if
      end if
!
!  Copy the word, and shift the string.
!
      word = s(iget1:iget2)

      call s_shift_left ( s, iget2 )

      return
      end
