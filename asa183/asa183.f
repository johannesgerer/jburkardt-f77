      function r4_random ( s1, s2, s3 )

c*********************************************************************72
c
cc R4_RANDOM returns a pseudorandom number between 0 and 1.
c
c  Discussion:
c
c    This function returns a pseudo-random number rectangularly distributed
c    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
c    of Applied Statistics (1984) volume 33), not as claimed in the
c    original article.
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Brian Wichman, David Hill.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Brian Wichman, David Hill,
c    Algorithm AS 183: An Efficient and Portable Pseudo-Random
c    Number Generator,
c    Applied Statistics,
c    Volume 31, Number 2, 1982, pages 188-190.
c
c  Parameters:
c
c    Input/output, integer S1, S2, S3, three values used as the
c    seed for the sequence.  These values should be positive
c    integers between 1 and 30,000.
c
c    Output, real R4_RANDOM, the next value in the sequence.
c
      implicit none

      integer s1
      integer s2
      integer s3
      real r4_random

      s1 = mod ( 171 * s1, 30269 )
      s2 = mod ( 172 * s2, 30307 )
      s3 = mod ( 170 * s3, 30323 )
 
      r4_random = mod ( real ( s1 ) / 30269.0E+00 
     &                + real ( s2 ) / 30307.0E+00 
     &                + real ( s3 ) / 30323.0E+00, 1.0E+00 )

      return
      end
      function r4_uni ( s1, s2 )

c*********************************************************************72
c
cc R4_UNI returns a pseudorandom number between 0 and 1.
c
c  Discussion:
c
c    This function generates uniformly distributed pseudorandom numbers
c    between 0 and 1, using the 32-bit generator from figure 3 of
c    the article by L'Ecuyer.
c
c    The cycle length is claimed to be 2.30584E+18.
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Pascal original version by Pierre L'Ecuyer
c    Modifications by John Burkardt
c
c  Reference:
c
c    Pierre LEcuyer,
c    Efficient and Portable Combined Random Number Generators,
c    Communications of the ACM,
c    Volume 31, Number 6, June 1988, pages 742-751.
c
c  Parameters:
c
c    Input/output, integer S1, S2, two values used as the
c    seed for the sequence.  On first call, the user should initialize
c    S1 to a value between 1 and 2147483562;  S2 should be initialized
c    to a value between 1 and 2147483398.
c
c    Output, real R4_UNI, the next value in the sequence.
c
      implicit none

      integer k
      real r4_uni
      integer s1
      integer s2
      integer z
    
      k = s1 / 53668
      s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
      if ( s1 .lt. 0 ) then
        s1 = s1 + 2147483563
      end if

      k = s2 / 52774
      s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
      if ( s2 .lt. 0 ) then
        s2 = s2 + 2147483399
      end if

      z = s1 - s2
      if ( z .lt. 1 ) then
        z = z + 2147483562
      end if

      r4_uni = real ( z ) / 2147483563.0E+00

      return
      end
      function r8_random ( s1, s2, s3 )

c*********************************************************************72
c
cc R8_RANDOM returns a pseudorandom number between 0 and 1.
c
c  Discussion:
c
c    This function returns a pseudo-random number rectangularly distributed
c    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
c    of Applied Statistics (1984) volume 33), not as claimed in the
c    original article.
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Brian Wichman, David Hill.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Brian Wichman, David Hill,
c    Algorithm AS 183: An Efficient and Portable Pseudo-Random
c    Number Generator,
c    Applied Statistics,
c    Volume 31, Number 2, 1982, pages 188-190.
c
c  Parameters:
c
c    Input/output, integer S1, S2, S3, three values used as the
c    seed for the sequence.  These values should be positive
c    integers between 1 and 30,000.
c
c    Output, double precision R8_RANDOM, the next value in the sequence.
c
      implicit none

      integer s1
      integer s2
      integer s3
      double precision r8_random

      s1 = mod ( 171 * s1, 30269 )
      s2 = mod ( 172 * s2, 30307 )
      s3 = mod ( 170 * s3, 30323 )
 
      r8_random = mod ( dble ( s1 ) / 30269.0D+00 
     &                + dble ( s2 ) / 30307.0D+00 
     &                + dble ( s3 ) / 30323.0D+00, 1.0D+00 )

      return
      end
      function r8_uni ( s1, s2 )

c*********************************************************************72
c
cc R8_UNI returns a pseudorandom number between 0 and 1.
c
c  Discussion:
c
c    This function generates uniformly distributed pseudorandom numbers
c    between 0 and 1, using the 32-bit generator from figure 3 of
c    the article by L'Ecuyer.
c
c    The cycle length is claimed to be 2.30584E+18.
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Pascal original version by Pierre L'Ecuyer
c    Modifications by John Burkardt
c
c  Reference:
c
c    Pierre LEcuyer,
c    Efficient and Portable Combined Random Number Generators,
c    Communications of the ACM,
c    Volume 31, Number 6, June 1988, pages 742-751.
c
c  Parameters:
c
c    Input/output, integer S1, S2, two values used as the
c    seed for the sequence.  On first call, the user should initialize
c    S1 to a value between 1 and 2147483562;  S2 should be initialized
c    to a value between 1 and 2147483398.
c
c    Output, double precision R8_UNI, the next value in the sequence.
c
      implicit none

      integer k
      double precision r8_uni
      integer s1
      integer s2
      integer z
    
      k = s1 / 53668
      s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
      if ( s1 .lt. 0 ) then
        s1 = s1 + 2147483563
      end if

      k = s2 / 52774
      s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
      if ( s2 .lt. 0 ) then
        s2 = s2 + 2147483399
      end if

      z = s1 - s2
      if ( z .lt. 1 ) then
        z = z + 2147483562
      end if

      r8_uni = dble ( z ) / 2147483563.0D+00

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
