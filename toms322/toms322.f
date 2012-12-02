      function prob ( lm, ln, x )

c*********************************************************************72
c
cc PROB evaluates several probability density functions.
c
c  Discussion:
c
c    Probabilities are returned as the integral of P(q)dq for q
c    in the range zero to X for the F-ratio and chi-square distributions
c    and from -X to +X for the Student's T and Normal distributions.
c
c    The returned probability will be in the range 0 to 1, unless an error
c    occurred, in which case -1.0 is returned.
c
c    Recoded by LW, 28-Mar-94.
c
c    The interpretation of the meaning of the input parameters for this 
c    routine depends on the quantity being computed.
c
c    For F-ratios:
c
c      LM = numerator degrees of freedom, with 1 < LM < 300.
c      LN = denominator degrees of freedom, with 1 < LN < 5000. 
c      X = F-ratio
c
c    For Student's T (two tailed):
c
c      LM = 1
c      LN = Degrees of freedom, with 1 < LN < 5000.
c      X = The square of t
c
c    For normal deviates (two tailed):
c
c      LM = 1
c      LN = 5000
c      X = The square of the deviate
c
c    For chi-square:
c
c      LM = Degrees of freedom, with 1 < LM < 300.
c      LN = 5000
c      X = Chi-square/LM
c
c  Modified:
c
c    28 December 2007
c
c  Author:
c
c    Egon Dorrer
c    Modifications by John Burkardt
c
c  Reference:
c
c    Egon Dorrer, 
c    Algorithm 322: F-Distribution,
c    Communications of the ACM, 
c    Volume 11, Number 2, 1968, pages 116-117.
c
c    JBF Field, 
c    Certification of Algorithm 322,
c    Communications of the ACM, 
c    Volume 12, Number 1, 1969, page 39.
c
c    Hubert Tolman, 
c    Remark on Algorithm 322,
c    Communications of the ACM,
c    Volume 14, Number 2, 1979, page 117.
c
c  Parameters:
c
c    Input, integer LM, integer LN, double precision X, specify the 
c    distribution and its parameters.
c
c    Output, double precision PROB, the value of the probability density function at X.
c
      implicit none

      integer a
      integer b
      double precision d
      integer i
      integer j
      integer lm
      integer ln
      integer m
      integer n
      double precision p
      double precision prob
      double precision w
      double precision x
      double precision y
      double precision z
      double precision zk

      m  = min ( lm, 300 )
      n  = min ( ln, 5000 )

      if ( min ( m, n ) .lt. 1 ) then
        prob = -1.0D+00
        return
      end if

      a = 2 * ( m / 2 ) - m + 2
      b = 2 * ( n / 2 ) - n + 2
      w = x * dble ( m ) / dble ( n )
      z = 1.0D+00 / ( 1.0D+00 + w )

      if ( a .eq. 1 ) then

        if ( b .eq. 1 ) then
          p = sqrt ( w )
          y = 0.3183098862D+00
          d = y * z / p
          p = 2.0D+00 * y * atan ( p )
        else
          p = sqrt ( w * z )
          d = p * z / ( 2.0D+00 * w )
        end if

      else

        if ( b .eq. 1 ) then
          p = sqrt ( z )
          d = z * p / 2.0D+00
          p = 1.0D+00 - p
        else
          d = z * z
          p = w * z
        end if

      end if

      y = 2.0D+00 * w / z

      if ( a .eq. 1 ) then

        do j = b + 2, n, 2
          d = ( 1.0D+00 + dble ( a ) / dble ( j - 2 ) ) * d * z
          p = p + d * y / dble ( j - 1 )
        end do

      else

        zk = z**( ( n - 1 ) / 2 )
        d  = d * zk * dble ( n ) / dble ( b )
        p  = p * zk + w * z * ( zk - 1.0D+00 ) / ( z - 1.0D+00 )

      end if

      y = w * z
      z = 2.0D+00 / z
      b = n - 2

      do i = a + 2, m, 2
        j = i + b
        d = ( y * d * dble ( j ) ) / dble ( i - 2 )
        p = p - z * d / dble ( j )
      end do

      prob = max ( 0.0D+00, min ( 1.0D+00, p ) )

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
