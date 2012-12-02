      subroutine doomsday_gregorian ( y, m, d, w )

c*****************************************************************************80
c
cc DOOMSDAY_GREGORIAN: weekday given any date in Gregorian calendar.
c
c  Discussion:
c
c    This procedure does not include any procedure to switch to the Julian
c    calendar for dates early enough that that calendar was used instead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Conway,
c    Tomorrow is the Day After Doomsday,
c    Eureka,
c    Volume 36, October 1973, pages 28-31.
c
c  Parameters:
c
c    Input, integer Y, M, D, the year, month and day of the date.
c    Note that the year must be positive.
c
c    Output, integer W, the weekday of the date.
c
      implicit none

      integer anchor(4)
      integer c
      integer d
      integer drd
      integer drdr
      integer i4_wrap
      integer l
      integer m
      integer mdoom(12)
      integer w
      integer y
      integer ydoom
      logical year_is_leap_gregorian
      integer yy
      integer yy12d
      integer yy12r
      integer yy12r4d

      save anchor
      save mdoom

      data anchor / 1, 6, 4, 3 /
      data mdoom / 3, 28, 0, 4, 9, 6, 11, 8, 5, 10, 7, 12 /
c
c  Refuse to handle Y <= 0.
c
      if ( y .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DOOMSDAY_GREGORIAN - Fatal error!'
        write ( *, '(a)' ) '  Y <= 0.'
        stop
      end if
c
c  Determine the century C.
c
      c = y / 100
c
c  Determine the last two digits of the year, YY
c
      yy = mod ( y, 100 )
c
c  Divide the last two digits of the year by 12.
c
      yy12d = yy / 12
      yy12r = mod ( yy, 12 ) 
      yy12r4d = yy12r / 4
      drd = yy12d + yy12r + yy12r4d

      drdr = mod ( drd, 7 )
      ydoom = anchor( mod ( c-1, 4 ) + 1 ) + drdr
      ydoom = i4_wrap ( ydoom, 1, 7 )
c
c  If M = 1 or 2, and Y is a leap year, add 1.
c
      if ( ( m .eq. 1 .or. m .eq. 2 ) .and. 
     &  year_is_leap_gregorian ( y ) ) then
        l = 1
      else
        l = 0
      end if

      w = ydoom + ( d -  mdoom(m) - l )
      w = i4_wrap ( w, 1, 7 )


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
      function year_is_leap_gregorian ( y )

c*********************************************************************72
c
cc YEAR_IS_LEAP_GREGORIAN returns TRUE if the Gregorian year was a leap year.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer Y, the year to be checked.
c
c    Output, logical YEAR_IS_LEAP_GREGORIAN, TRUE if the year was a leap year,
c    FALSE otherwise.
c
      implicit none

      integer y
      logical year_is_leap_gregorian

      if ( y .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'YEAR_IS_LEAP_GREGORIAN - Fatal error'
        write ( *, '(a)' ) 
     &    '  This function will not accept nonpositive years.'
        stop
      end if

      if ( mod ( y, 400 ) .eq. 0 ) then
        year_is_leap_gregorian = .true.
      else if ( mod ( y, 100 ) .eq. 0 ) then
        year_is_leap_gregorian = .false.
      else if ( mod ( y, 4 ) .eq. 0 ) then
        year_is_leap_gregorian = .true.
      else
        year_is_leap_gregorian = .false.
      end if

      return
      end
