      function a_to_i4 ( ch )

c*********************************************************************72
c
cc A_TO_I4 returns the index of an alphabetic character.
c
c  Discussion:
c
c    Instead of ICHAR, we now use the IACHAR function, which
c    guarantees the ASCII collating sequence.
c
c  Example:
c
c    CH  A_TO_I4
c
c    'A'   1
c    'B'   2
c    ...
c    'Z'  26
c    'a'  27
c    'b'  28
c    ...
c    'z'  52
c    '$'   0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, a character.
c
c    Output, integer A_TO_I4, is the alphabetic index of the character,
c    between 1 and 26 if the character is a capital letter,
c    between 27 and 52 if it is lower case, and 0 otherwise.
c
      implicit none

      integer a_to_i4
      integer cap_shift
      parameter ( cap_shift = 64 )
      character ch
      integer low_shift
      parameter ( low_shift = 96 )

      if ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) then
        a_to_i4 = ichar ( ch ) - cap_shift
      else if ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) then
        a_to_i4 = ichar ( ch ) - low_shift + 26
      else
        a_to_i4 = 0
      end if

      return
      end
      subroutine binary_to_i4 ( s, i )

c*********************************************************************72
c
cc BINARY_TO_I4 converts a binary representation into an integer value.
c
c  Example:
c
c        S        I
c
c      '101'      5
c    '-1000'     -8
c        '1'      1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, the binary representation.
c
c    Output, integer I, the integer whose representation was input.
c
      implicit none

      character c
      integer i
      integer ichr
      integer isgn
      character*(*) s
      integer s_len_trim
      integer s_length
      integer state

      s_length = s_len_trim ( s )

      i = 0
      ichr = 1
      state = 0
      isgn = 1

10    continue

      if ( ichr .le. s_length ) then

        c = s(ichr:ichr)
c
c  Blank.
c
        if ( c .eq. ' ' ) then

          if ( state .eq. 2 ) then
            state = 3
          end if
c
c  Sign, + or -.
c
        else if ( c .eq. '-' ) then

          if ( state .eq. 0 ) then
            state = 1
            isgn = -1
          else
            state = -1
          end if

        else if ( c .eq. '+' ) then

          if ( state .eq. 0 ) then
            state = 1
          else
            state = -1
          end if
c
c  Digit, 0 or 1.
c
        else if ( c .eq. '1' ) then

          i = 2 * i
          i = i + 1
          state = 2

        else if ( c .eq. '0' ) then

          i = 2 * i
          state = 2
c
c  Illegal or unknown sign.
c
        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BINARY_TO_I4 - Serious error!'
          write ( *, '(a)' ) '  Illegal digit = "' // c // '"'
          write ( *, '(a)' ) '  Conversion halted prematurely!'
          return

        end if

        if ( state .eq. -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BINARY_TO_I4 - Serious error!'
          write ( *, '(a)' ) '  Unable to decipher input!'
          return
        end if

        if ( 3 .le. state ) then
          go to 20
        end if

        ichr = ichr + 1

        end if

      go to 10
c
c  Apply the sign.
c
20    continue

      i = isgn * i

      return
      end
      subroutine binary_to_r4 ( s, r )

c*********************************************************************72
c
cc BINARY_TO_R4 converts a binary representation into an R4.
c
c  Discussion:
c
c    An "R4" value is simply a real number to be stored as a
c    variable of type "real".
c
c  Example:
c
c        S                         R
c
c    -1010.11                   -10.75
c    0.011011                     0.4218750
c    0.01010101010101010101010    0.3333333
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the binary representation.
c
c    Output, real R, the real number.
c
      implicit none

      character c
      integer ichr
      integer intval
      integer isgn
      integer power
      real r
      character * ( * ) s
      integer s_length
      integer state

      s_length = len_trim ( s )

      intval = 0
      ichr = 1
      state = 0
      isgn = 1
      r = 0.0E+00
      power = 0

10    continue

      if ( ichr .le. s_length ) then

        c = s(ichr:ichr)
c
c  Blank.
c
        if ( c .eq. ' ' ) then

          if ( state .eq. 4 ) then
            state = 5
          end if
c
c  Sign, + or -.
c
        else if ( c .eq. '-' ) then

          if ( state .eq. 0 ) then
            state = 1
            isgn = -1
          else
            state = -1
          end if

        else if ( c .eq. '+' ) then

          if ( state .eq. 0 ) then
            state = 1
          else
            state = -1
          end if
c
c  Digit, 0 or 1.
c
        else if ( c .eq. '1' ) then

          intval = 2 * intval + 1

          if ( state .eq. 0 .or. state .eq. 1 ) then
            state = 2
          else if ( state .eq. 3 ) then
            state = 4
          end if

          if ( state .eq. 4 ) then
            power = power + 1
          end if

        else if ( c .eq. '0' ) then

          intval = 2 * intval

          if ( state .eq. 0 .or. state .eq. 1 ) then
            state = 2
          else if ( state .eq. 3 ) then
            state = 4
          end if

          if ( state .eq. 4 ) then
            power = power + 1
          end if
c
c  Decimal point.
c
        else if ( c .eq. '.' ) then

          if ( state .le. 2 ) then
            state = 3
          else
            state = -1
          end if
c
c  Illegal or unknown sign.
c
        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BINARY_TO_R4 - Serious error!'
          write ( *, '(a)' ) '  Illegal character = "' // c // '"'
          write ( *, '(a)' ) '  Conversion halted prematurely!'
          stop

        end if

        if ( state .eq. -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BINARY_TO_R4 - Serious error!'
          write ( *, '(a)' ) '  Unable to decipher input!'
          stop
        end if

        if ( 5 .le. state ) then
          go to 20
        end if

        ichr = ichr + 1

        go to 10

      end if

20    continue
c
c  Apply the sign and the scale factor.
c
      r = real ( isgn * intval ) / 2.0E+00**power

      return
      end
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
      function ch_index_first ( s, ch )

c*********************************************************************72
c
cc CH_INDEX_FIRST is the first occurrence of a character in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, the string to be searched.
c
c    Input, character CH, the character to be searched for.
c
c    Output, integer CH_INDEX_FIRST, the location of the first
c    occurrence of the character in the string, or -1 if it does not occur.
c
      implicit none

      character ch
      integer ch_index_first
      integer i
      character*(*) s
      integer s_len_trim
      integer s_length

      ch_index_first = -1
      s_length = s_len_trim ( s )

      do i = 1, s_length

        if ( s(i:i) .eq. ch ) then
          ch_index_first = i
          return
        end if

      end do

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
      function ch_is_alpha ( ch )

c*********************************************************************72
c
cc CH_IS_ALPHA is TRUE if CH is an alphabetic character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, a character to check.
c
c    Output, logical CH_IS_ALPHA is TRUE if CH is an alphabetic character.
c
      implicit none

      character ch
      logical ch_is_alpha

      if ( ( lle ( 'a', ch ) .and. lle ( ch, 'z' ) ) .or.
     &     ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) ) then
        ch_is_alpha = .true.
      else
        ch_is_alpha = .false.
      end if

      return
      end
      function ch_is_alphanumeric ( ch )

c*********************************************************************72
c
cc CH_IS_ALPHANUMERIC is TRUE if CH is alphanumeric.
c
c  Discussion:
c
c    Alphanumeric characters are 'A' through 'Z', 'a' through 'z' and
c    '0' through '9'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, the character to be checked.
c
c    Output, logical CH_IS_ALPHANUMERIC, is TRUE if the character is
c    alphabetic or numeric.
c
      implicit none

      character ch
      logical ch_is_alphanumeric
      integer i

      i = ichar ( ch )

      if ( ( 65 .le. i .and. i .le. 90 ) .or.
     &     ( 97 .le. i .and. i .le. 122 ) .or.
     &     ( 48 .le. i .and. i .le. 57 ) ) then

        ch_is_alphanumeric = .true.

      else

        ch_is_alphanumeric = .false.

      end if

      return
      end
      function ch_is_control ( ch )

c*********************************************************************72
c
cc CH_IS_CONTROL is TRUE if a character is a control character.
c
c  Discussion:
c
c    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, the character to be tested.
c
c    Output, logical CH_IS_CONTROL, TRUE if the character is a control
c    character, and FALSE otherwise.
c
      implicit none

      character ch
      logical ch_is_control

      if ( ichar ( ch ) .le. 31 .or. 127 .le. ichar ( ch ) ) then
        ch_is_control = .true.
      else
        ch_is_control = .false.
      end if

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
      function ch_is_space ( ch )

c*********************************************************************72
c
cc CH_IS_SPACE is TRUE if a character is a whitespace character.
c
c  Discussion:
c
c    A whitespace character is a space, a form feed, a newline,
c    a carriage return, a tab, or a vertical tab.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, a character to check.
c
c    Output, logical CH_IS_SPACE is TRUE if the character is a whitespace
c    character.
c
      implicit none

      character ch
      logical ch_is_space

      if ( ch .eq. ' ' ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 12 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 10 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 13 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 9 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 11 ) ) then
        ch_is_space = .true.
      else
        ch_is_space = .false.
      end if

      return
      end
      function ch_is_upper ( ch )

c*********************************************************************72
c
cc CH_IS_UPPER is TRUE if CH is an upper case letter.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, the character to be analyzed.
c
c    Output, logical CH_IS_UPPER, is TRUE if CH is an upper case letter.
c
      implicit none

      character ch
      logical ch_is_upper

      if ( lle ( 'A', ch ) .and. lle ( ch, 'Z' ) ) then
        ch_is_upper = .true.
      else
        ch_is_upper = .false.
      end if

      return
      end
      subroutine ch_low ( ch )

c*********************************************************************72
c
cc CH_LOW lowercases a single character.
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
c    Input/output, character CH, the character to be lowercased.
c
      implicit none

      character ch
      integer i

      i = ichar ( ch )

      if ( 65 .le. i .and. i .le. 90 ) then
        ch = char ( i + 32 )
      end if

      return
      end
      function ch_roman_to_i4 ( ch )

c*********************************************************************72
c
cc CH_ROMAN_TO_I4 returns the integer value of a single Roman digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, a Roman digit.
c
c    Output, integer CH_ROMAN_TO_I4, the value of the Roman
c    numeral.  If the Roman numeral was not recognized, 0 is returned.
c
      implicit none

      character ch
      integer ch_roman_to_i4
      integer i

      if ( ch .eq. 'M' .or. ch .eq. 'm' ) then
        i = 1000
      else if ( ch .eq. 'D' .or. ch .eq. 'd' ) then
        i = 500
      else if ( ch .eq. 'C' .or. ch .eq. 'c' ) then
        i = 100
      else if ( ch .eq. 'L' .or. ch .eq. 'l' ) then
        i = 50
      else if ( ch .eq. 'X' .or. ch .eq. 'x' ) then
        i = 10
      else if ( ch .eq. 'V' .or. ch .eq. 'v' ) then
        i = 5
      else if ( ch .eq. 'I' .or. ch .eq. 'i' .or.
     &          ch .eq. 'J' .or. ch .eq. 'j' ) then
        i = 1
      else
        i = 0
      end if

      ch_roman_to_i4 = i

      return
      end
      function ch_scrabble ( tile )

c*********************************************************************72
c
cc CH_SCRABBLE returns the character on a given Scrabble tile.
c
c  Discussion:
c
c    The tiles are numbered 1 to 100, and are labeled 'A' through 'Z',
c    plus two blanks.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer TILE, the index of the desired Scrabble tile.
c
c    Output, character CH_SCRABBLE, the character on the given tile.
c
      implicit none

      character ch_scrabble
      character scrabble(100)
      integer tile

      save scrabble

      data scrabble /
     &  'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B',
     &  'B', 'C', 'C', 'D', 'D', 'D', 'D', 'E', 'E', 'E',
     &  'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'F',
     &  'F', 'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'I',
     &  'I', 'I', 'I', 'I', 'I', 'J', 'K', 'L', 'L', 'L',
     &  'L', 'M', 'M', 'N', 'N', 'N', 'N', 'N', 'N', 'O',
     &  'O', 'O', 'O', 'O', 'O', 'O', 'O', 'P', 'P', 'Q',
     &  'R', 'R', 'R', 'R', 'R', 'R', 'S', 'S', 'S', 'S',
     &  'T', 'T', 'T', 'T', 'T', 'T', 'U', 'U', 'U', 'U',
     &  'V', 'V', 'W', 'W', 'X', 'X', 'Y', 'Z', ' ', ' ' /

      if ( 1 .le. tile .and. tile .le. 100 ) then
        ch_scrabble = scrabble(tile)
      else
        ch_scrabble = '?'
      end if

      return
      end
      subroutine ch_swap ( ch1, ch2 )

c*********************************************************************72
c
cc CH_SWAP swaps two characters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 July 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH1, CH2.  On output, the values
c    have been interchanged.
c
      implicit none

      character ch1
      character ch2
      character ch3

      ch3 = ch1
      ch1 = ch2
      ch2 = ch3

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
      subroutine ch_to_digit_bin ( ch, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT_BIN returns the integer value of a binary digit.
c
c  Discussion:
c
c    This routine handles other traditional binary pairs of "digits"
c    besides '0' and '1'.
c
c  Example:
c
c     CH  DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    'T'    1
c    'F'    0
c    'Y'    1
c    'N'    0
c    '+'    1
c    '-'    0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, the binary digit.
c
c    Output, integer DIGIT, the corresponding integer value.
c    If CH was 'illegal', then DIGIT is -1.
c
      implicit none

      character ch
      integer digit

      if ( ch .eq. '0' .or.
     &     ch .eq. 'F' .or.
     &     ch .eq. 'f' .or.
     &     ch .eq. '-' .or.
     &     ch .eq. 'N' .or.
     &     ch .eq. 'n' ) then

        digit = 0

      else if ( ch .eq. '1' .or.
     &          ch .eq. 'T' .or.
     &          ch .eq. 't' .or.
     &          ch .eq. '+' .or.
     &          ch .eq. 'Y' .or.
     &          ch .eq. 'y' ) then

        digit = 1

      else

        digit = -1

      end if

      return
      end
      subroutine ch_to_digit_hex ( ch, i )

c*********************************************************************72
c
cc CH_TO_DIGIT_HEX returns the integer value of a hexadecimal digit.
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
c    Input, character CH, the hexadecimal digit, '0'
c    through '9', or 'A' through 'F', or also 'a' through 'f'
c    are allowed.
c
c    Output, integer I, the corresponding integer, or -1 if CH was illegal.
c
      implicit none

      character ch
      integer i

      i = ichar ( ch )

      if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

        i = i - 48

      else if ( 65 .le. i .and. i .le. 70 ) then

        i = i - 55

      else if ( 97 .le. i .and. i .le. 102 ) then

        i = i - 87

      else if ( ch .eq. ' ' ) then

        i = 0

      else

        i = -1

      end if

      return
      end
      subroutine ch_to_digit_oct ( ch, i )

c*********************************************************************72
c
cc CH_TO_DIGIT_OCT returns the integer value of an octal digit.
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
c    Input, character CH, the octal digit, '0' through '7'.
c
c    Output, integer I, the corresponding integer value, or
c    -1 if CH was illegal.
c
      implicit none

      character ch
      integer i

      i = ichar ( ch )

      if ( lle ( '0', ch ) .and. lle ( ch, '7' ) ) then

        i = i - 48

      else if ( ch .eq. ' ' ) then

        i = 0

      else

        i = -1

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
      function ch_to_scrabble ( ch )

c*********************************************************************72
c
cc CH_TO_SCRABBLE returns the Scrabble index of a character.
c
c  Discussion:
c
c    'A' through 'Z' have indices 1 through 26, and blank is index 27.
c    Case is ignored.  All other characters return index -1.
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
c    Input, character CH, the character.
c
c    Output, integer CH_TO_SCRABBLE, the Scrabble index of the character.
c
      implicit none

      integer a_to_i4
      character ch
      character ch_copy
      integer ch_to_scrabble
      integer ic

      if ( ch .eq. ' ' ) then
        ch_to_scrabble = 27
        return
      end if

      ch_copy = ch
      call ch_cap ( ch_copy )
      ic = a_to_i4 ( ch_copy )

      if ( 1 .le. ic .and. ic .le. 26 ) then
        ch_to_scrabble = ic
      else
        ch_to_scrabble = -1
      end if

      return
      end
      function ch_uniform ( clo, chi, seed )

c*********************************************************************72
c
cc CH_UNIFORM returns a scaled pseudorandom CH.
c
c  Discussion:
c
c    A CH is an alphabetic character value.
c
c    The value is scaled to lie between characters CLO and CHI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CLO, CHI, the minimum and maximum acceptable characters.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, character CH_UNIFORM, the randomly chosen character.
c
      implicit none

      character ch_uniform
      character chi
      character clo
      integer i
      integer ihi
      integer ilo
      real r4_uniform_01
      integer seed

      if ( seed == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CH_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      ilo = ichar ( clo )
      ihi = ichar ( chi )

      i = ilo + int ( r4_uniform_01 ( seed ) * real ( ihi + 1 - ilo ) )

      i = max ( i, ilo )
      i = min ( i, ihi )

      ch_uniform = char ( i )

      return
      end
      subroutine digit_bin_to_ch ( i, ch )

c*********************************************************************72
c
cc DIGIT_BIN_TO_CH returns the character representation of a binary digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer, between 0 and 1.
c
c    Output, character CH, the character representation of the integer.
c
      implicit none

      character ch
      integer i

      if ( i .eq. 0 ) then
        ch = '0'
      else if ( i .eq. 1 ) then
        ch = '1'
      else
        ch = '*'
      end if

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
      subroutine digit_oct_to_ch ( i, ch )

c*********************************************************************72
c
cc DIGIT_OCT_TO_CH returns the character representation of an octal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer, between 0 and 7.
c
c    Output, character CH, the character representation of the integer.
c
      character ch
      integer i

      if ( 0 .le. i .and. i .le. 7 ) then
        ch = char ( i + 48 )
      else
        ch = '*'
      end if

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
      subroutine file_name_inc ( file_name )

c*********************************************************************72
c
cc FILE_NAME_INC generates the next filename in a series.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  If this number is all 9's on input, the output number
c    is all 0's.  Non-numeric letters of the name are unaffected, and
c    if the name contains no digits, then nothing is done.
c
c  Example:
c
c      Input          Output
c      -----          ------
c      a7to11.txt     a7to12.txt
c      a7to99.txt     a8to00.txt
c      a9to99.txt     a0to00.txt
c      cat.txt        cat.txt
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) FILE_NAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      character*(*) file_name
      integer i
      integer lens

      lens = len ( file_name )

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do

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
      subroutine hex_to_i4 ( s, intval )

c***********************************************************************72
c
cc HEX_TO_I4 converts a hexadecimal string to its integer value.
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
c    Input, character ( len = * ) S, the string of hexadecimal digits.
c
c    Output, integer INTVAL, the corresponding integer value.
c
      implicit none

      integer first
      integer idig
      integer intval
      integer isgn
      integer j
      character ( len = * ) s
      integer s_len_trim
      integer s_length

      s_length = s_len_trim ( s )
c
c  Determine if there is a plus or minus sign.
c
      isgn = 1

      first = s_length + 1

      do j = 1, s_length

        if ( s(j:j) .eq. '-' ) then
          isgn = -1
        else if ( s(j:j) .eq. '+' ) then
          isgn = + 1
        else if ( s(j:j) .ne. ' ' ) then
          first = j
          exit
        end if

      end do
c
c  Read the numeric portion of the string.
c
      intval = 0

      do j = first, s_length
        call ch_to_digit_hex ( s(j:j), idig )
        intval = intval * 16 + idig
      end do

      intval = isgn * intval

      return
      end
      subroutine hex_to_s ( hex, s )

c*********************************************************************72
c
cc HEX_TO_S converts a hexadecimal string into characters.
c
c  Example:
c
c    Input:
c
c      '414243'
c
c    Output:
c
c      'ABC'.
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
c    Input, character ( len = * ) HEX, a string of pairs of hexadecimal values.
c
c    Output, character ( len = * ) S, the corresponding character string.
c
      implicit none

      character ( len = * ) hex
      integer i
      integer intval
      integer j
      integer ndo
      integer nhex
      character ( len = * ) s
      integer s_len_trim
      integer s_length

      s_length = len ( s )
      nhex = s_len_trim ( hex )
      ndo = min ( nhex / 2, s_length )

      s = ' '

      do i = 1, ndo
        j = 2 * i - 1
        call hex_to_i4 ( hex(j:j+1), intval )
        s(i:i) = char ( intval )
      end do

      return
      end
      function i4_gcd ( i, j )

c*********************************************************************72
c
cc I4_GCD finds the greatest common divisor of I and J.
c
c  Discussion:
c
c    Only the absolute values of I and J are
c    considered, so that the result is always nonnegative.
c
c    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
c
c    If I and J have no common factor, I4_GCD is returned as 1.
c
c    Otherwise, using the Euclidean algorithm, I4_GCD is the
c    largest common factor of I and J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common divisor
c    is desired.
c
c    Output, integer I4_GCD, the greatest common divisor of I and J.
c
      implicit none

      integer i
      integer i4_gcd
      integer ip
      integer iq
      integer ir
      integer j

      i4_gcd = 1
c
c  Return immediately if either I or J is zero.
c
      if ( i .eq. 0 ) then
        i4_gcd = max ( 1, abs ( j ) )
        return
      else if ( j .eq. 0 ) then
        i4_gcd = max ( 1, abs ( i ) )
        return
      end if
c
c  Set IP to the larger of I and J, IQ to the smaller.
c  This way, we can alter IP and IQ as we go.
c
      ip = max ( abs ( i ), abs ( j ) )
      iq = min ( abs ( i ), abs ( j ) )
c
c  Carry out the Euclidean algorithm.
c
10    continue

        ir = mod ( ip, iq )

        if ( ir .eq. 0 ) then
          go to 20
        end if

        ip = iq
        iq = ir

      go to 10

20    continue

      i4_gcd = iq

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function i4_length ( i )

c*********************************************************************72
c
cc I4_LENGTH computes the number of characters needed to print an integer.
c
c  Example:
c
c         I    I4_LENGTH
c
c         0       1
c         1       1
c        -1       2
c      1952       4
c    123456       6
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
c    Input, integer I, the integer whose length is desired.
c
c    Output, integer I4_LENGTH, the number of characters required
c    to print the integer.
c
      implicit none

      integer i
      integer i_copy
      integer i4_length

      if ( i .lt. 0 ) then
        i4_length = 1
      else if ( i .eq. 0 ) then
        i4_length = 1
        return
      else if ( 0 .lt. i ) then
        i4_length = 0
      end if

      i_copy = abs ( i )

 10   continue

      if ( i_copy == 0 ) go to 20

        i4_length = i4_length + 1

        i_copy = i_copy / 10

      go to 10

20    continue

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      function i4_to_a ( i )

c*********************************************************************72
c
cc I4_TO_A returns the I-th alphabetic character.
c
c  Example:
c
c    I  I4_TO_A
c
c   -8  ' '
c    0  ' '
c    1  'A'
c    2  'B'
c   ..
c   26  'Z'
c   27  'a'
c   52  'z'
c   53  ' '
c   99  ' '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the letter to be returned.
c    0 is a space;
c    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
c    27 through 52 requests 'a' through 'z', (ASCII 97:122);
c
c    Output, character I4_TO_A, the requested alphabetic letter.
c
      implicit none

      integer cap_shift
      parameter ( cap_shift = 64 )
      integer i
      character i4_to_a
      integer low_shift
      parameter ( low_shift = 96 )

      if ( i .le. 0 ) then
        i4_to_a = ' '
      else if ( 1 .le. i .and. i .le. 26 ) then
        i4_to_a = char ( cap_shift + i )
      else if ( 27 .le. i .and. i .le. 52 ) then
        i4_to_a = char ( low_shift + i - 26 )
      else if ( 53 .le. i ) then
        i4_to_a = ' '
      end if

      return
      end
      subroutine i4_to_amino_code ( i, ch )

c*********************************************************************72
c
cc I4_TO_AMINO_CODE converts an integer to an amino code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl Branden, John Tooze,
c    Introduction to Protein Structure,
c    Garland Publishing, 1991.
c
c  Parameters:
c
c    Input, integer I, the index of an amino acid, between 1
c    and 23.
c
c    Output, character CH, the one letter code for an amino acid.
c
      implicit none

      integer n
      parameter ( n = 23 )

      character ch
      character ch_table(n)
      integer i

      save ch_table

      data ch_table /
     &  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
     &  'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
     &  'X', 'Y', 'Z' /

      if ( 1 .le. i .and. i .le. 23 ) then
        ch = ch_table(i)
      else
        ch = '?'
      end if

      return
      end
      subroutine i4_to_hex_digit ( i, ch )

c*********************************************************************72
c
cc I4_TO_HEX_DIGIT converts a (small) I4 to a hexadecimal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer, between 0 and 15.
c
c    Output, character CH, the hexadecimal digit corresponding to the integer.
c
      implicit none

      character ch
      integer i

      if ( 0 .le. i .and. i .le. 9 ) then
        ch = char ( i + 48 )
      else if ( 10 .le. i .and. i .le. 15 ) then
        ch = char ( i + 55 )
      else
        ch = '*'
      end if

      return
      end
      function i4_to_isbn ( i )

c*********************************************************************72
c
cc I4_TO_ISBN converts an I4 to an ISBN digit.
c
c  Discussion:
c
c    Only the integers 0 through 10 can be input.  The representation
c    of 10 is 'X'.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, integer I, an integer between 0 and 10.
c
c    Output, character I4_TO_ISBN, the ISBN character code of the integer.
c    If I is illegal, then I4_TO_ISBN is set to '?'.
c
      implicit none

      integer i
      character i4_to_isbn

           if ( i .eq. 0 ) then
        i4_to_isbn = '0'
      else if ( i .eq. 1 ) then
        i4_to_isbn = '1'
      else if ( i .eq. 2 ) then
        i4_to_isbn = '2'
      else if ( i .eq. 3 ) then
        i4_to_isbn = '3'
      else if ( i .eq. 4 ) then
        i4_to_isbn = '4'
      else if ( i .eq. 5 ) then
        i4_to_isbn = '5'
      else if ( i .eq. 6 ) then
        i4_to_isbn = '6'
      else if ( i .eq. 7 ) then
        i4_to_isbn = '7'
      else if ( i .eq. 8 ) then
        i4_to_isbn = '8'
      else if ( i .eq. 9 ) then
        i4_to_isbn = '9'
      else if ( i .eq. 10 ) then
        i4_to_isbn = 'X'
      else
        i4_to_isbn = '?'
      end if

      return
      end
      subroutine i4_to_month_abb ( m, month_abb )

c*********************************************************************72
c
cc I4_TO_MONTH_ABB returns the 3 character abbreviation of a given month.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the index of the month, which should
c    be between 1 and 12.
c
c    Output, character * 3 MONTH_ABB, the month abbreviation
c
      implicit none

      character * 3 abb(12)
      integer m
      character * 3 month_abb

      save abb

      data abb /
     &  'Jan', 'Feb', 'Mar', 'Apr',
     &  'May', 'Jun', 'Jul', 'Aug',
     &  'Sep', 'Oct', 'Nov', 'Dec' /

      if ( m .lt. 1 .or. 12 .lt. m ) then

        month_abb = '???'

      else

        month_abb = abb(m)

      end if

      return
      end
      subroutine i4_to_s_left ( i4, s )

c*********************************************************************72
c
cc I4_TO_S_LEFT converts an I4 to a left-justified string.
c
c  Discussion:
c
c    An I4 is an integer.
c
c  Example:
c
c    Assume that S is 6 characters long:
c
c        I4  S
c
c         1  1
c        -1  -1
c         0  0
c      1952  1952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, an integer to be converted.
c
c    Output, character * ( * ) S, the representation of the integer.
c    The integer will be left-justified.  If there is not enough space,
c    the string will be filled with stars.
c
      implicit none

      character c
      integer i
      integer i4
      integer idig
      integer ihi
      integer ilo
      integer ipos
      integer ival
      character * ( * ) s

      s = ' '

      ilo = 1
      ihi = len ( s )

      if ( ihi .le. 0 ) then
        return
      end if
c
c  Make a copy of the integer.
c
      ival = i4
c
c  Handle the negative sign.
c
      if ( ival .lt. 0 ) then

        if ( ihi .le. 1 ) then
          s(1:1) = '*'
          return
        end if

        ival = -ival
        s(1:1) = '-'
        ilo = 2

      end if
c
c  The absolute value of the integer goes into S(ILO:IHI).
c
      ipos = ihi
c
c  Find the last digit of IVAL, strip it off, and stick it into the string.
c
10    continue

        idig = mod ( ival, 10 )
        ival = ival / 10

        if ( ipos .lt. ilo ) then
          do i = 1, ihi
            s(i:i) = '*'
          end do
          return
        end if

        call digit_to_ch ( idig, c )

        s(ipos:ipos) = c
        ipos = ipos - 1

        if ( ival .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue
c
c  Shift the string to the left.
c
      s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
      s(ilo+ihi-ipos:ihi) = ' '

      return
      end
      subroutine i4_to_s_zero ( intval, s )

c*********************************************************************72
c
cc I4_TO_S_ZERO converts an integer to a string, with zero padding.
c
c  Example:
c
c    Assume that S is 6 characters long:
c
c    INTVAL  S
c
c         1  000001
c        -1  -00001
c         0  000000
c      1952  001952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, an integer to be converted.
c
c    Output, character * ( * ) S, the representation of the integer.
c    The integer will be right justified, and zero padded.
c    If there is not enough space, the string will be filled with stars.
c
      implicit none

      character c
      integer i
      integer idig
      integer ihi
      integer ilo
      integer intval
      integer ipos
      integer ival
      character * ( * ) s

      s = ' '

      ilo = 1
      ihi = len ( s )

      if ( ihi .le. 0 ) then
        return
      end if
c
c  Make a copy of the integer.
c
      ival = intval
c
c  Handle the negative sign.
c
      if ( ival .lt. 0 ) then

        if ( ihi .le. 1 ) then
          s(1:1) = '*'
          return
        end if

        ival = - ival
        s(1:1) = '-'
        ilo = 2

      end if
c
c  Working from right to left, strip off the digits of the integer
c  and place them into S(ILO:IHI).
c
      ipos = ihi

10    continue

      if ( ival .ne. 0 .or. ipos .eq. ihi ) then

        idig = mod ( ival, 10 )
        ival = ival / 10

        if ( ipos < ilo ) then
          do i = 1, ihi
            s(i:i) = '*'
          end do
          return
        end if

        call digit_to_ch ( idig, c )

        s(ipos:ipos) = c
        ipos = ipos - 1

        go to 10

      end if 
c
c  Fill the empties with zeroes.
c
      do i = ilo, ipos
        s(i:i) = '0'
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
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer s_len_trim
      character*(*) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,i12)' ) i, a(i)
      end do

      return
      end
      function isbn_to_i4 ( c )

c*********************************************************************72
c
cc ISBN_TO_I4 converts an ISBN character into an integer.
c
c  Discussion:
c
c    The characters '0' through '9' stand for themselves, but
c    the character 'X' or 'x' stands for 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, character C, the ISBN character code to be converted.
c
c    Output, integer ISBN_TO_I4, the numeric value of the character
c    code, between 0 and 10.  This value is returned as -1 if C is
c    not a valid character code.
c
      implicit none

      character c
      integer isbn_to_i4

           if ( c .eq. '0' ) then
        isbn_to_i4 = 0
      else if ( c .eq. '1' ) then
        isbn_to_i4 = 1
      else if ( c .eq. '2' ) then
        isbn_to_i4 = 2
      else if ( c .eq. '3' ) then
        isbn_to_i4 = 3
      else if ( c .eq. '4' ) then
        isbn_to_i4 = 4
      else if ( c .eq. '5' ) then
        isbn_to_i4 = 5
      else if ( c .eq. '6' ) then
        isbn_to_i4 = 6
      else if ( c .eq. '7' ) then
        isbn_to_i4 = 7
      else if ( c .eq. '8' ) then
        isbn_to_i4 = 8
      else if ( c .eq. '9' ) then
        isbn_to_i4 = 9
      else if ( c .eq. 'X' .or. c == 'x' ) then
        isbn_to_i4 = 10
      else
        isbn_to_i4 = -1
      end if

      return
      end
      subroutine number_inc ( s )

c*********************************************************************72
c
cc NUMBER_INC increments the integer represented by a string.
c
c  Example:
c
c    Input      Output
c    -----      ------
c    '17'       '18'
c    'cat3'     'cat4'
c    '2for9'    '3for0'
c    '99thump'  '00thump'
c
c  Discussion:
c
c    If the string contains characters that are not digits, they will
c    simply be ignored.  If the integer is all 9's on input, then
c    the output will be all 0's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, a string representing an integer.
c
      implicit none

      logical ch_is_digit
      integer i
      character*(*) s
      integer s_len
      integer s_len_trim

      s_len = s_len_trim ( s )

      do i = s_len, 1, -1

        if ( ch_is_digit ( s(i:i) ) ) then

          call digit_inc ( s(i:i) )

          if ( s(i:i) .ne. '0' ) then
            return
          end if

        end if

      end do

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
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
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      if ( seed == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
      subroutine s_adjustl ( s )

c*********************************************************************72
c
cc S_ADJUSTL flushes a string left.
c
c  Discussion:
c
c    Both blanks and tabs are treated as "white space".
c
c    This routine is similar to the FORTRAN90 ADJUSTL routine.
c
c  Example:
c
c    Input             Output
c
c    '     Hello'      'Hello     '
c    ' Hi therec  '    'Hi therec   '
c    'Fred  '          'Fred  '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 Jun3 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S.
c    On input, S is a string of characters.
c    On output, any initial blank or tab characters have been cut.
c
      implicit none

      integer i
      integer nonb
      character * ( * ) s
      integer s_length
      character tab

      tab = char ( 9 )
c
c  Check the length of the string to the last nonblank.
c  If nonpositive, return.
c
      s_length = len_trim ( s )

      if ( s_length .le. 0 ) then
        return
      end if
c
c  Find NONB, the location of the first nonblank, nontab.
c
      nonb = 0

      do i = 1, s_length

        if ( s(i:i) .ne. ' ' .and. s(i:i) .ne. tab ) then
          nonb = i
          go to 10
        end if

      end do

10    continue

      if ( nonb .eq. 0 ) then
        s = ' '
        return
      end if
c
c  Shift the string left.
c
      if ( 1 .lt. nonb ) then
        do i = 1, s_length + 1 - nonb
          s(i:i) = s(i+nonb-1:i+nonb-1)
        end do
      end if
c
c  Blank out the end of the string.
c
      s(s_length+2-nonb:s_length) = ' '

      return
      end
      subroutine s_adjustr ( s )

c*********************************************************************72
c
cc S_ADJUSTR flushes a string right.
c
c  Example:
c
c    Input             Output
c    'Hello     '      '     Hello'
c    ' Hi there!  '    '   Hi there!'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, on output, trailing blank
c    characters have been cut, and pasted back onto the front.
c
      implicit none

      integer i
      integer nonb
      character * ( * ) s
      integer s_length
c
c  Check the full length of the string.
c
      s_length = len ( s )
c
c  Find the occurrence of the last nonblank.
c
      nonb = len_trim ( s )
c
c  Shift the string right.
c
      do i = s_length, s_length + 1 - nonb, -1
        s(i:i) = s(i-s_length+nonb:i-s_length+nonb)
      end do
c
c  Blank out the beginning of the string.
c
      s(1:s_length-nonb) = ' '

      return
      end
      subroutine s_after_ss_copy ( s1, ss, s2 )

c*********************************************************************72
c
cc S_AFTER_SS_COPY copies a string after a given substring.
c
c  Discussion:
c
c    S1 and S2 can be the same object, in which case the string is
c    overwritten by a copy of itself after the substring.
c
c  Example:
c
c    Input:
c
c      S1 = 'ABCDEFGH'
c      SS = 'EF'
c
c    Output:
c
c      S2 = 'GH'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S1, the string to be copied.
c
c    Input, character*(*) SS, the substring after which the copy begins.
c
c    Output, character*(*) S2, the copied portion of S.
c
      implicit none

      integer first
      integer last
      integer last_s2
      character*(*)  s1
      integer s1_length
      character*(*)  s2
      character*(*)  ss
c
c  Find the first occurrence of the substring.
c
      first = index ( s1, ss )
c
c  If the substring doesn't occur at all, then S2 is blank.
c
      if ( first .eq. 0 ) then
        s2 = ' '
        return
      end if
c
c  Redefine FIRST to point to the first character to copy after
c  the substring.
c
      first = first + len ( ss )
c
c  Measure the two strings.
c
      s1_length = len ( s1 )
      last_s2 = len ( s2 )
c
c  Adjust effective length of S if S2 is short.
c
      last = min ( s1_length, last_s2 + first - 1 )
c
c  Copy the string.
c
      s2(1:s1_length+1-first) = s1(first:s1_length)
c
c  Clear out the rest of the copy.
c
      s2(s1_length+2-first:last_s2) = ' '

      return
      end
      function s_begin ( s1, s2 )

c*********************************************************************72
c
cc S_BEGIN is TRUE if one string matches the beginning of the other.
c
c  Discussion:
c
c    The strings are compared, ignoring blanks, spaces and capitalization.
c
c  Example:
c
c     S1              S2      S_BEGIN
c
c    'Bob'          'BOB'     TRUE
c    '  B  o b '    ' bo b'   TRUE
c    'Bob'          'Bobby'   TRUE
c    'Bobo'         'Bobb'    FALSE
c    ' '            'Bob'     FALSE    (Do not allow a blank to match
c                                       anything but another blank string.)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( ) S1, S2, the strings to be compared.
c
c    Output, logical S_BEGIN, is TRUE if the strings match up to
c    the end of the shorter string, ignoring case.
c
      implicit none

      logical ch_eqi
      integer i1
      integer i2
      logical s_begin
      character * ( * )  s1
      integer s1_length
      character * ( * )  s2
      integer s2_length

      s1_length = len_trim ( s1 )
      s2_length = len_trim ( s2 )
c
c  If either string is blank, then both must be blank to match.
c  Otherwise, a blank string matches anything, which is not
c  what most people want.
c
      if ( s1_length .eq. 0 .or. s2_length .eq. 0 ) then

        if ( s1_length .eq. 0 .and. s2_length .eq. 0 ) then
          s_begin = .true.
        else
          s_begin = .false.
        end if

        return

      end if

      i1 = 0
      i2 = 0
c
c  Find the next nonblank in S1.
c
10    continue

20      continue

          i1 = i1 + 1

          if ( s1_length .lt. i1 ) then
            s_begin = .true.
            return
          end if

          if ( s1(i1:i1) .ne. ' ' ) then
            go to 30
          end if

        go to 20

30      continue
c
c  Find the next nonblank in S2.
c
40      continue

          i2 = i2 + 1

          if ( s2_length .lt. i2 ) then
            s_begin = .true.
            return
          end if

          if ( s2(i2:i2) .ne. ' ' ) then
            go to 50
          end if

        go to 40

50      continue
c
c  If the characters match, get the next pair.
c
        if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
          go to 60
        end if

      go to 10

60    continue

      s_begin = .false.

      return
      end
      subroutine s_blank_delete ( s )

c*********************************************************************72
c
cc S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
c
c  Discussion:
c
c    All TAB characters are also removed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      character ch
      integer get
      integer put
      character*(*) s
      integer s_len_trim
      integer s_length
      character tab

      tab = char ( 9 )

      put = 0
      s_length = s_len_trim ( s )

      do get = 1, s_length

        ch = s(get:get)

        if ( ch .ne. ' ' .and. ch .ne. tab ) then
          put = put + 1
          s(put:put) = ch
        end if

      end do

      s(put+1:s_length) = ' '

      return
      end
      subroutine s_blanks_delete ( s )

c*********************************************************************72
c
cc S_BLANKS_DELETE replaces consecutive blanks by one blank.
c
c  Discussion:
c
c    Thanks to Bill Richmond for pointing out a programming flaw which
c    meant that, as characters were slid to the left through multiple
c    blanks, their original images were not blanked out.  This problem
c    is easiest resolved by using a copy of the string.
c
c    The remaining characters are left justified and right padded with blanks.
c    TAB characters are converted to spaces.
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
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      integer i
      integer j
      character newchr
      character oldchr
      character*(*) s
      character*255 s_copy
      integer s_length
      character tab

      tab = char ( 9 )

      s_length = len ( s )

      j = 0
      s_copy(1:s_length) = s(1:s_length)
      s(1:s_length) = ' '

      newchr = ' '

      do i = 1, s_length

        oldchr = newchr
        newchr = s_copy(i:i)

        if ( newchr == tab ) then
          newchr = ' '
        end if

        if ( oldchr .ne. ' ' .or. newchr .ne. ' ' ) then
          j = j + 1
          s(j:j) = newchr
        end if

      end do

      return
      end
      subroutine s_cap ( s )

c*********************************************************************72
c
cc S_CAP replaces any lowercase letters by uppercase ones in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to be transformed.
c
      implicit none

      character ch
      integer i
      character * ( * )  s
      integer s_length

      s_length = len_trim ( s )

      do i = 1, s_length

        ch = s(i:i)
        call ch_cap ( ch )
        s(i:i) = ch

      end do

      return
      end
      subroutine s_cat ( s1, s2, s3 )

c*********************************************************************72
c
cc S_CAT concatenates two strings to make a third string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S1, the "prefix" string.
c
c    Input, character * ( * ) S2, the "postfix" string.
c
c    Output, character * ( * ) S3, the string made by
c    concatenating S1 and S2, ignoring any trailing blanks.
c
      implicit none

      character * ( * ) s1
      character * ( * ) s2
      character * ( * ) s3

      s3 = trim ( s1 ) // trim ( s2 )

      return
      end
      subroutine s_cat1 ( s1, s2, s3 )

c*********************************************************************72
c
cc S_CAT1 concatenates two strings, with a single blank separator.
c
c  Example:
c
c    S1       S2       S
c
c    'cat'    'dog'    'cat dog'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S1, the "prefix" string.
c
c    Input, character * ( * ) S2, the "postfix" string.
c
c    Output, character * ( * ) S3, the string made by concatenating
c    S1, a blank, and S2, ignoring any trailing blanks.
c
      implicit none

      character * ( * ) s1
      character * ( * ) s2
      character * ( * ) s3

      if ( len_trim ( s1 ) .eq. 0 .and. len_trim ( s2 ) .eq. 0 ) then
        s3 = ' '
      else if ( len_trim ( s1 ) .eq. 0 ) then
        s3 = s2
      else if ( len_trim ( s2 ) .eq. 0 ) then
        s3 = s1
      else
        s3 = trim ( s1 ) // ' ' // trim ( s2 )
      end if

      return
      end
      subroutine s_ch_blank ( s, ch )

c*********************************************************************72
c
cc S_CH_BLANK replaces each occurrence of a particular character by a blank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to be transformed.
c
c    Input, character CH, the character to be removed.
c
      implicit none

      character ch
      integer i
      character * ( * ) s
      integer s_length

      s_length = len_trim ( s )

      do i = 1, s_length

        if ( s(i:i) .eq. ch ) then
          s(i:i) = ' '
        end if

      end do

      return
      end
      subroutine s_ch_count ( s, ch, ch_count )

c*********************************************************************72
c
cc S_CH_COUNT counts occurrences of a particular character in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string.
c
c    Input, character CH, the character to be counted.
c
c    Output, integer CH_COUNT, the number of occurrences.
c 
      implicit none

      character ch
      integer ch_count
      integer i
      character * ( * ) s
      integer s_length

      ch_count = 0

      s_length = len ( s )

      do i = 1, s_length

        if ( s(i:i) .eq. ch ) then
          ch_count = ch_count + 1
        end if

      end do

      return
      end
      subroutine s_ch_delete ( s, ch )

c*********************************************************************72
c
cc S_CH_DELETE removes all occurrences of a character from a string.
c
c  Discussion:
c
c    Each time the given character is found in the string, the characters
c    to the right of the string are shifted over one position.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to be transformed.
c
c    Input, character CH, the character to be removed.
c
      implicit none

      character ch
      integer get
      integer put
      character * ( * ) s
      integer s_length

      s_length = len_trim ( s )

      put = 1

      do get = 1, s_length

        if ( s(get:get) .eq. ch ) then

        else if ( put .eq. get ) then
          put = put + 1
        else
          s(put:put) = s(get:get)
          put = put + 1
        end if

      end do

      s(put:s_length) = ' '

      return
      end
      subroutine s_chop ( s, ilo, ihi )

c*********************************************************************72
c
cc S_CHOP "chops out" a portion of a string, and closes up the hole.
c
c  Example:
c
c    S = 'Fred is not a jerk!'
c
c    call s_chop ( S, 9, 12 )
c
c    S = 'Fred is a jerk!    '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string to be transformed.
c
c    Input, integer ILO, IHI, the locations of the first and last
c    characters to be removed.
c
      implicit none

      integer ihi
      integer ihi2
      integer ilo
      integer ilo2
      character * ( * )  s
      integer s_length

      s_length = len ( s )

      ilo2 = max ( ilo, 1 )
      ihi2 = min ( ihi, s_length )

      if ( ihi2 .lt. ilo2 ) then
        return
      end if

      s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
      s(s_length+ilo2-ihi2:s_length) = ' '

      return
      end
      subroutine s_control_blank ( s )

c*********************************************************************72
c
cc S_CONTROL_BLANK replaces control characters with blanks.
c
c  Discussion:
c
c    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      logical ch_is_control
      integer i
      character*(*) s
      integer s_len_trim
      integer s_length

      s_length = s_len_trim ( s )

      do i = 1, s_length
        if ( ch_is_control ( s(i:i) ) ) then
          s(i:i) = ' '
        end if
      end do

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
      function s_first_nonblank ( s )

c*********************************************************************72
c
cc S_FIRST_NONBLANK returns the location of the first nonblank.
c
c  Discussion:
c
c    If all characters are blanks, a 0 is returned.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) S, the string to be examined.
c
c    Output, integer S_FIRST_NONBLANK, the location of the first
c    nonblank character in the string, or 0 if all are blank.
c
      implicit none

      integer i
      character * ( * ) s
      integer s_first_nonblank
      integer s_length

      s_length = len ( s )

      do i = 1, s_length

        if ( s(i:i) .ne. ' ' ) then
          s_first_nonblank = i
          return
        end if

      end do

      s_first_nonblank = 0

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
      function s_indexi ( s, sub )

c*********************************************************************72
c
cc S_INDEXI is a case-insensitive INDEX function.
c
c  Discussion:
c
c    The function returns the location in the string at which the
c    substring SUB is first found, or 0 if the substring does not
c    occur at all.
c
c    The routine is also trailing blank insensitive.  This is very
c    important for those cases where you have stored information in
c    larger variables.  If S is of length 80, and SUB is of
c    length 80, then if S = 'FRED' and SUB = 'RED', a match would
c    not be reported by the standard FORTRAN INDEX, because it treats
c    both variables as being 80 characters long!  This routine assumes that
c    trailing blanks represent garbage!
c
c    Because of the suppression of trailing blanks, this routine cannot be
c    used to find, say, the first occurrence of the two-character
c    string 'A '.  However, this routine treats as a special case the
c    occurrence where S or SUB is entirely blank.  Thus you can
c    use this routine to search for occurrences of double or triple blanks
c    in a string, for example, although INDEX itself would be just as
c    suitable for that problem.
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
c    Output, integer S_INDEXI.  0 if SUB does not occur in
c    the string.  Otherwise S(S_INDEXI:S_INDEXI+LENS-1) = SUB,
c    where LENS is the length of SUB, and is the first place
c    this happens.  However, note that this routine ignores case,
c    unlike the standard FORTRAN INDEX function.
c
      implicit none

      integer i
      integer llen2
      character * ( * ) s
      logical s_eqi
      integer s_indexi
      integer s_length
      character * ( * ) sub

      s_indexi = 0

      s_length = len_trim ( s )
      llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
      if ( s_length .eq. 0 ) then
        s_length = len ( s )
      end if

      if ( llen2 .eq. 0 ) then
        llen2 = len ( sub )
      end if

      if ( s_length .lt. llen2 ) then
        return
      end if

      do i = 1, s_length + 1 - llen2

        if ( s_eqi ( s(i:i+llen2-1), sub ) ) then
          s_indexi = i
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
      subroutine s_to_i4vec ( s, n, ivec, ierror )

c*********************************************************************72
c
cc S_TO_I4VEC reads an I4VEC from a string.
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
c    Output, integer IVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer ilo
      integer ivec(n)
      integer length
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + length

      go to 10

      end if

20    continue

      return
      end
      function s_to_l ( s )

c*********************************************************************72
c
cc S_TO_L reads a logical value from a string.
c
c  Discussion:
c
c    There are several ways of representing logical data that this routine
c    recognizes:
c
c      False   True
c      -----   ----
c
c      0       1
c      F       T
c      f       t
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be read.
c
c    Output, logical S_TO_L, the logical value read from the string.
c
      implicit none

      integer i
      character * ( * ) s
      integer s_length
      logical s_to_l

      s_length = len_trim ( s )

      do i = 1, s_length

        if ( s(i:i) .eq. '0' .or.
     &       s(i:i) .eq. 'F' .or.
     &       s(i:i) .eq. 'f' ) then
          s_to_l = .false.
          return
        else if ( s(i:i) .eq. '1' .or.
     &            s(i:i) .eq. 'T' .or.
     &            s(i:i) .eq. 't' ) then
          s_to_l = .true.
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_TO_L - Fatal error!'
      write ( *, '(a)' ) '  Input string did not contain logical data.'

      stop
      end
      subroutine s_to_r4 ( s, dval, ierror, length )

c*********************************************************************72
c
cc S_TO_R4 reads an R4 from a string.
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
c    02 January 2009
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
c    Output, real DVAL, the value read from the string.
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
      real dval
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
      real rbot
      real rexp
      real rtop
      character * ( * ) s
      integer s_len_trim

      nchar = s_len_trim ( s )

      ierror = 0
      dval = 0.0E+00
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
            rtop = 10.0E+00 * rtop + real ( ndig )
          else if ( ihave .eq. 5 ) then
            rtop = 10.0E+00 * rtop + real ( ndig )
            rbot = 10.0E+00 * rbot
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
        write ( *, '(a)' ) 'S_TO_R4 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a,a)' ) '    ', s
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0E+00
      else
        if ( jbot .eq. 1 ) then
          rexp = 10.0E+00 ** ( jsgn * jtop )
        else
          rexp = 10.0E+00 ** ( real ( jsgn * jtop ) / real ( jbot ) )
        end if
      end if

      dval = real ( isgn ) * rexp * rtop / rbot

      return
      end
      subroutine s_to_r4vec ( s, n, rvec, ierror )

c*********************************************************************72
c
cc S_TO_R4VEC reads an R4VEC from a string.
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
c    Input, character * ( * ) S, the string to be read.
c
c    Input, integer N, the number of values expected.
c
c    Output, real RVEC(N), the values read from the string.
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
      real rvec(n)
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_r4 ( s(ilo:), rvec(i), ierror, lchar )

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
      subroutine s_word_extract_first ( s, w )

c*********************************************************************72
c
cc S_WORD_EXTRACT_FIRST extracts the first word from a string.
c
c  Discussion:
c
c    A "word" is a string of characters terminated by a blank or
c    the end of the string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string.  On output, the first
c    word has been removed, and the remaining string has been shifted left.
c
c    Output, character * ( * ) W, the leading word of the string.
c
      implicit none

      integer get1
      integer get2
      character * ( * ) s
      integer s_len_trim
      integer s_length
      character * ( * ) w

      w = ' '

      s_length = s_len_trim ( s )

      if ( s_length .lt. 1 ) then
        return
      end if
c
c  Find the first nonblank.
c
      get1 = 0

10    continue

        get1 = get1 + 1

        if ( s_length .lt. get1 ) then
          return
        end if

        if ( s(get1:get1) .ne. ' ' ) then
          go to 20
        end if

      go to 10

20    continue
c
c  Look for the last contiguous nonblank.
c
      get2 = get1

30    continue

        if ( s_length .le. get2 ) then
          go to 40
        end if

        if ( s(get2+1:get2+1) .eq. ' ' ) then
          go to 40
        end if

        get2 = get2 + 1

      go to 30

40    continue
c
c  Copy the word.
c
      w = s(get1:get2)
c
c  Shift the string.
c
      s(1:get2) = ' '
      call s_adjustl ( s )

      return
      end
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf
c    This FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
