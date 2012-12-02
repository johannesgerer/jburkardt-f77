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
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

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
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      subroutine jed_to_weekday ( jed, w, f )

c*********************************************************************72
c
cc JED_TO_WEEKDAY computes the day of the week from a JED.
c
c  Discussion:
c
c    BC 4713/01/01 => JED = 0.0 was noon on a Monday.
c
c    jedmod = mod ( 0.0D+00, 7.0D+00 ) = 0.0D+00
c    j = mod ( nint ( 0 ), 7 ) = 0
c    f = ( 0.0D+00 + 0.5D+00 ) - real ( j ) = 0.5D+00
c    w = i4_wrap ( 0 + 2, 1, 7 ) = 2 = MONDAY
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
c  Reference:
c
c    Edward Richards,
c    Mapping Time, The Calendar and Its History,
c    Oxford, 1999.
c
c  Parameters:
c
c    Input, double precision JED, the Julian Ephemeris Date.
c
c    Output, integer W, the day of the week of the date.
c    The days are numbered from Sunday through Saturday, 1 through 7.
c
c    Output, double precision F, the fractional part of the day.
c
      implicit none

      double precision f
      integer i4_wrap
      integer j
      double precision jed
      double precision jedmod
      integer w

      jedmod = mod ( jed, 7.0D+00 )

      j = mod ( nint ( jedmod ), 7 )

      f = ( jedmod + 0.5D+00 ) - dble ( j )

      w = i4_wrap ( j + 2, 1, 7 )

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
      subroutine weekday_to_name_common ( w, s )

c*********************************************************************72
c
cc WEEKDAY_TO_NAME_COMMON returns the name of a Common weekday.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer W, the weekday index.
c
c    Output, character * ( * ) S, the weekday name.
c
      implicit none

      character * 9 name(7)
      character * ( * ) s
      integer w
      integer w2

      save name

      data name /
     &  'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday',
     &  'Thursday ', 'Friday   ', 'Saturday ' /
c
c  Make a local copy of the weekday number.
c
      w2 = w
c
c  Return the weekday name.
c
      s = name ( w2 )

      return
      end
      subroutine weekday_values ( n_data, y, m, d, w )

c*********************************************************************72
c
cc WEEKDAY_VALUES returns the day of the week for various dates.
c
c  Discussion:
c
c    The CE or Common Era calendar is used, under the
c    hybrid Julian/Gregorian Calendar, with a transition from Julian
c    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
c
c    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
c    years BC/BCE are indicated by a negative year value.
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
c  Reference:
c
c    Edward Reingold, Nachum Dershowitz,
c    Calendrical Calculations: The Millennium Edition,
c    Cambridge University Press, 2001,
c    ISBN: 0 521 77752 6
c    LC: CE12.R45.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0
c    before the first call.  On each call, the routine increments N_DATA by 1,
c    and returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer Y, M, D, the Common Era date.
c
c    Output, integer W, the day of the week.  Sunday = 1.
c
      implicit none

      integern_max
      parameter ( n_max = 34 )

      integer d
      integer d_vec(n_max)
      integer m
      integer m_vec(n_max)
      integer n_data
      integer w
      integer w_vec(n_max)
      integer y
      integer y_vec(n_max)

      save d_vec
      save m_vec
      save w_vec
      save y_vec

      data d_vec /
     &  30,
     &   8,
     &  26,
     &   3,
     &   7,
     &  18,
     &   7,
     &  19,
     &  14,
     &  18,
     &  16,
     &   3,
     &  26,
     &  20,
     &   4,
     &  25,
     &  31,
     &   9,
     &  24,
     &  10,
     &  30,
     &  24,
     &  19,
     &   2,
     &  27,
     &  19,
     &  25,
     &  29,
     &  19,
     &   7,
     &  17,
     &  25,
     &  10,
     &  18 /
      data m_vec /
     &   7,
     &  12,
     &   9,
     &  10,
     &   1,
     &   5,
     &  11,
     &   4,
     &  10,
     &   5,
     &   3,
     &   3,
     &   3,
     &   4,
     &   6,
     &   1,
     &   3,
     &   9,
     &   2,
     &   6,
     &   6,
     &   7,
     &   6,
     &   8,
     &   3,
     &   4,
     &   8,
     &   9,
     &   4,
     &  10,
     &   3,
     &   2,
     &  11,
     &   7 /
      data w_vec /
     &  1,
     &  4,
     &  4,
     &  1,
     &  4,
     &  2,
     &  7,
     &  1,
     &  7,
     &  1,
     &  6,
     &  7,
     &  6,
     &  1,
     &  1,
     &  4,
     &  7,
     &  7,
     &  7,
     &  4,
     &  1,
     &  6,
     &  1,
     &  2,
     &  4,
     &  1,
     &  1,
     &  2,
     &  2,
     &  5,
     &  3,
     &  1,
     &  4,
     &  1 /
      data y_vec /
     &  - 587,
     &  - 169,
     &     70,
     &    135,
     &    470,
     &    576,
     &    694,
     &   1013,
     &   1066,
     &   1096,
     &   1190,
     &   1240,
     &   1288,
     &   1298,
     &   1391,
     &   1436,
     &   1492,
     &   1553,
     &   1560,
     &   1648,
     &   1680,
     &   1716,
     &   1768,
     &   1819,
     &   1839,
     &   1903,
     &   1929,
     &   1941,
     &   1943,
     &   1943,
     &   1992,
     &   1996,
     &   2038,
     &   2094 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        y = 0
        m = 0
        d = 0
        w = 0
      else
        y = y_vec(n_data)
        m = m_vec(n_data)
        d = d_vec(n_data)
        w = w_vec(n_data)
      end if

      return
      end
      subroutine y_common_to_astronomical ( y, y2 )

c*********************************************************************72
c
cc Y_COMMON_TO_ASTRONOMICAL converts a Common year to an Astronomical year.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y, the Common year.
c
c    Output, integer Y2, the Astronomical year.
c
      implicit none

      integer y
      integer y2

      if ( y .lt. 0 ) then
        y2 = y + 1
      else if ( y .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Y_COMMON_TO_ASTRONOMICAL - Fatal error!'
        write ( *, '(a)' ) '  COMMON calendar does not have a year 0.'
        stop
      else
        y2 = y
      end if

      return
      end
      subroutine ymd_to_s_common ( y, m, d, s )

c*********************************************************************72
c
cc YMD_TO_S_COMMON writes a Common YMD date into a string.
c
c  Format:
c
c    CE YYYY/MM/DD
c    BCE YYYY/MM/DD
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y, M, D, the YMD date.
c
c    Output, character * ( * ) S, a representation of the date.
c
      implicit none

      integer d
      integer d2
      integer ierror
      integer m
      integer m2
      character * ( 20 ) s1
      character * ( 2 ) s2
      character * ( * )  s
      integer y
      integer y2
c
c  Copy the input.
c
      y2 = y
      m2 = m
      d2 = d

      if ( 0 .le. y2 ) then
        s1 = 'CE '
        call i4_to_s_left ( y2, s1(4:) )
      else
        s1 = 'BCE '
        call i4_to_s_left (  - y2, s1(5:) )
      end if

      call s_cat ( s1, '/', s1 )

      call i4_to_s_zero ( m2, s2 )

      call s_cat ( s1, s2, s1 )

      call s_cat ( s1, '/', s1 )

      call i4_to_s_zero ( d2, s2 )

      call s_cat ( s1, s2, s )

      return
      end
      subroutine ymd_to_weekday_common ( y, m, d, w )

c*********************************************************************72
c
cc YMD_TO_WEEKDAY_COMMON returns the weekday of a Common YMD date.
c
c  Discussion:
c
c    The "common" calendar is meant to be the calendar which is Julian up to
c    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y, M, D, the YMD date.
c
c    Output, integer W, is the week day number of the date, with
c    1 for Sunday, through 7 for Saturday.
c
      implicit none

      integer d
      double precision f
      double precision f2
      double precision jed
      integer m
      integer w
      integer y

      f = 0.5D+00

      call ymdf_to_jed_common ( y, m, d, f, jed )

      call jed_to_weekday ( jed, w, f2 )

      return
      end
      subroutine ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

c*********************************************************************72
c
cc YMDF_COMPARE compares two YMDF dates.
c
c  Discussion:
c
c    The comparison should work for a pair of dates in any calendar.
c
c    No check is made that the dates are actually legitimate.  It is
c    assumed that the calling routine has already ensured that the
c    dates are properly "normalized".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y1, M1, D1, double precision F1, the
c    first YMDF date.
c
c    Input, integer Y2, M2, D2, double precision F2, the
c    second YMDF date.
c
c    Output, character CMP:
c    '<' if date 1 precedes date 2;
c    '=' if date 1 equals date 2;
c    '>' if date 1 follows date 2;
c
      implicit none

      character cmp
      integer d1
      integer d2
      double precision f1
      double precision f2
      integer m1
      integer m2
      integer y1
      integer y2

      cmp = '?'
c
c  Compare years...
c
      if ( y1 .lt. y2 ) then
        cmp = '<'
      else if ( y1 .gt. y2 ) then
        cmp = '>'
      else
c
c  ...if necessary, compare months in equal years...
c
        if ( m1 .lt. m2 ) then
          cmp = '<'
        else if ( m1 .gt. m2 ) then
          cmp = '>'
        else
c
c  ...if necessary, compare days in equal months...
c
          if ( d1 .lt. d2 ) then
            cmp = '<'
          else if ( d1 .gt. d2 ) then
            cmp = '>'
          else
c
c  ...if necessary, compare fractional parts.
c
            if ( f1 .lt. f2 ) then
              cmp = '<'
            else if ( f1 .gt. f2 ) then
              cmp = '>'
            else
              cmp = '='
            end if

          end if

        end if

      end if

      return
      end
      subroutine ymdf_to_jed_common ( y, m, d, f, jed )

c*********************************************************************72
c
cc YMDF_TO_JED_COMMON converts a Common YMDF date to a JED.
c
c  Discussion:
c
c    The "common" calendar is meant to be the calendar which is Julian up to
c    day JED = 2299160, and Gregorian from day JED = 2299161 and after.
c
c    The Julian Ephemeris Date is essentially a count of the number
c    of days that have elapsed since noon, 1 January 4713 BC, at
c    Greenwich, England.  Strictly speaking, the Julian Ephemeris Date
c    is counted from noon, and thus day "0" began at noon on 1 January 4713 BC,
c    and ended at noon on 2 January 4713 BC.
c
c    The Julian Ephemeris Date was devised by Joseph Scaliger in 1583.
c
c    The Julian Ephemeris Date has been adopted by astronomers as
c    a convenient reference for dates.
c
c  Example:
c
c       Y   M     D         JED
c    --------------     -------
c    BC 4713 Jan  1           0
c    AD 1968 May 23     2440000
c    AD 1984 Dec 31     2446065
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y, M, D, double precision F, the YMDF date.
c
c    Output, double precision JED, the Julian Ephemeris Date.
c
      implicit none

      character cmp
      integer d
      integer d1
      integer d2
      double precision f
      double precision f1
      double precision f2
      integer ierror
      double precision jed
      integer m
      integer m1
      integer m2
      integer y
      integer y1
      integer y2
c
c  Copy the month and year.
c
      y1 = y
      m1 = m
      d1 = d
      f1 = f

      y2 = 1582
      m2 = 10
      d2 = 4+1
      f2 = 0.0D+00

      call ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

      if ( cmp .eq. '<' ) then
        call ymdf_to_jed_julian ( y1, m1, d1, f1, jed )
        return
      end if
c
c  Use the Gregorian calendar for dates strictly after 1752/9/13.
c
      y2 = 1582
      m2 = 10
      d2 = 15-1
      f2 = 0.0D+00

      call ymdf_compare ( y1, m1, d1, f1, y2, m2, d2, f2, cmp )

      if ( cmp .eq. '>' ) then
        call ymdf_to_jed_gregorian ( y1, m1, d1, f1, jed )
        return
      end if

      jed = -1.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'YMDF_TO_JED_COMMON - Error!'
      write ( *, '(a)' ) '  Illegal date!'

      return
      end
      subroutine ymdf_to_jed_gregorian ( y, m, d, f, jed )

c*********************************************************************72
c
cc YMDF_TO_JED_GREGORIAN converts a Gregorian YMDF date to a JED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Richards,
c    Algorithm E,
c    Mapping Time, The Calendar and Its History,
c    Oxford, 1999, pages 323-324.
c
c  Parameters:
c
c    Input, integer Y, M, D, double precision F, the YMDF date.
c
c    Output, double precision JED, the corresponding JED.
c
      implicit none

      integer d
      integer d_prime
      double precision f
      integer g
      integer ierror
      double precision jed
      integer j1
      integer j2
      integer m
      integer m_prime
      integer y
      integer y2
      integer y_prime
c
c  Account for the missing year 0 by moving negative years up one.
c
      call y_common_to_astronomical ( y, y2 )
c
c  Convert the calendar date to a computational date.
c
      y_prime = y2 + 4716 - ( 14 - m ) / 12
      m_prime = mod ( m + 9, 12 )
      d_prime = d - 1
c
c  Convert the computational date to a JED.
c
      j1 = ( 1461 * y_prime ) / 4

      j2 = ( 153 * m_prime + 2 ) / 5

      g = ( 3 * ( ( y_prime + 184 ) / 100 ) / 4 ) - 38

      jed = dble ( j1 + j2 + d_prime - 1401 - g ) - 0.5D+00
      jed = jed + f

      return
      end
      subroutine ymdf_to_jed_julian ( y, m, d, f, jed )

c*********************************************************************72
c
cc YMDF_TO_JED_JULIAN converts a Julian YMDF date to a JED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Richards,
c    Algorithm E,
c    Mapping Time, The Calendar and Its History,
c    Oxford, 1999, pages 323-324.
c
c  Parameters:
c
c    Input, integer Y, M, D, double precision F, the YMDF date.
c
c    Output, double precision JED, the Julian Ephemeris Date.
c
      implicit none

      integer d
      integer d_prime
      double precision f
      integer ierror
      double precision jed
      integer j1
      integer j2
      integer m
      integer m_prime
      integer y
      integer y2
      integer y_prime
c
c  Account for the missing year 0 by moving negative years up one.
c
      call y_common_to_astronomical ( y, y2 )
c
c  Convert the calendar date to a computational date.
c
      y_prime = y2 + 4716 - ( 14 - m ) / 12
      m_prime = mod ( m + 9, 12 )
      d_prime = d - 1
c
c  Convert the computational date to a JED.
c
      j1 = ( 1461 * y_prime ) / 4

      j2 = ( 153 * m_prime + 2 ) / 5

      jed = dble ( j1 + j2 + d_prime - 1401 ) - 0.5D+00
      jed = jed + f

      return
      end
