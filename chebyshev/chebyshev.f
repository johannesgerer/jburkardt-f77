      subroutine chebyshev_coefficients ( a, b, n, f, c )

c*********************************************************************72
c
cc CHEBYSHEV_COEFFICIENTS determines Chebyshev interpolation coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Second Edition,
c    Cambridge University Press, 1992,
c    ISBN: 0-521-43064-X,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, double precision A, B, the domain of definition.
c
c    Input, integer N, the order of the interpolant.
c
c    Input, double precision, external :: F ( X ), an external function.
c
c    Output, double precision C(N), the Chebyshev coefficients.
c
      implicit none

      integer n

      double precision a
      double precision angle
      double precision b
      double precision c(n)
      double precision f
      external f
      double precision fx(n)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      do i = 1, n
        angle = dble ( 2 * i - 1 ) * pi / dble ( 2 * n )
        x = cos ( angle )
        x = 0.5D+00 * ( a + b ) + x * 0.5D+00 * ( b - a )
        fx(i) = f ( x );
      end do

      do i = 1, n
        c(i) = 0.0D+00
        do j = 1, n
          angle = dble ( ( i - 1 ) * ( 2 * j - 1 ) ) * pi 
     &      / dble ( 2 * n )
          c(i) = c(i) + fx(j) * cos ( angle )
        end do
      end do

      do i = 1, n
        c(i) = 2.0D+00 * c(i) / dble ( n )
      end do

      return
      end
      subroutine chebyshev_interpolant ( a, b, n, c, m, x, cf )

c*********************************************************************72
c
cc CHEBYSHEV_INTERPOLANT evaluates a Chebyshev interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Second Edition,
c    Cambridge University Press, 1992,
c    ISBN: 0-521-43064-X,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, double precision A, B, the domain of definition.
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision C(N), the Chebyshev coefficients.
c
c    Input, integer M, the number of points.
c
c    Input, double precision X(M), the point at which the polynomial is
c    to be evaluated.
c
c    Output, double precision CF(M), the value of the Chebyshev
c    polynomial at X.
c
      implicit none

      integer m
      integer n

      double precision a
      double precision b
      double precision c(n)
      double precision cf(m)
      double precision di
      double precision dip1
      double precision dip2
      integer i
      integer j
      double precision x(m)
      double precision y

      do j = 1, m

        dip1 = 0.0D+00
        di = 0.0D+00
        y = ( 2.0D+00 * x(j) - a  - b ) / ( b - a )

        do i = n, 2, -1
          dip2 = dip1
          dip1 = di
          di = 2.0D+00 * y * dip1 - dip2 + c(i)
        end do

        cf(j) = y * di - dip1 + 0.5D+00 * c(1)

      end do

      return
      end
      subroutine chebyshev_zeros ( n, x )

c*********************************************************************72
c
cc CHEBYSHEV_ZEROS returns zeroes of the Chebyshev polynomial T(N)(X).
c
c  Discussion:
c
c    We produce the Chebyshev zeros in ascending order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision X(N), the zeroes of T(N)(X).
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)

      do i = 1, n
        angle = dble ( 2 * ( n - i ) + 1 ) * pi / dble ( 2 * n )
        x(i) = cos ( angle )
      end do

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
