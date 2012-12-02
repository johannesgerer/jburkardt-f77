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
      subroutine file_name_inc ( file_name )

c*********************************************************************72
c
cc FILE_NAME_INC increments a partially numeric filename.
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
c    Input/output, character*(*) FILE_NAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      integer change
      integer digit
      character*(*) file_name
      integer i
      integer lens

      lens = len_trim ( file_name )

      if ( lens .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
        write ( *, '(a)' ) '  The input string is empty.'
        stop
      end if

      change = 0

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          change = change + 1

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do
c
c  No digits were found.  Return blank.
c
      if ( change .eq. 0 ) then
        file_name = ' '
        return
      end if

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
