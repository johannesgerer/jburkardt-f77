      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2006
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine sine_transform_data ( n, d, s )

c*********************************************************************72
c
cc SINE_TRANSFORM_DATA does a sine transform on a vector of data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, double precision D(N), the vector of data.
c
c    Output, double precision S(N), the sine transform coefficients.
c
      implicit none

      integer n

      double precision angle
      double precision d(n)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision s(n)

      do i = 1, n
        s(i) = 0.0D+00
        do j = 1, n
          angle = pi * dble ( i * j ) / dble ( n + 1 )
          s(i) = s(i) + sin ( angle ) * d(j)
        end do
        s(i) = s(i) * sqrt ( 2.0D+00 / dble ( n + 1 ) )
      end do

      return
      end
      subroutine sine_transform_function ( n, a, b, f, s )

c*********************************************************************72
c
cc SINE_TRANSFORM_FUNCTION does a sine transform on functional data.
c
c  Discussion:
c
c    The interval [A,B] is divided into N+1 intervals using N+2 points,
c    which are indexed by 0 through N+1.
c
c    The original function F(X) is regarded as the sum of a linear function 
c    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
c    which is 0 at A and B.
c
c    The sine transform coefficients for F2 are then computed.
c
c    To recover the interpolant of F(X), it is necessary to combine the
c    linear part F1 with the sine transform interpolant:
c
c      Interp(F)(X) = F1(X) + F2(X)
c
c    This can be done by calling SINE_TRANSFORM_INTERPOLANT().
c    
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, double precision A, B, the interval endpoints.
c
c    Input, external, double precision F, a pointer to the function.
c
c    Output, double precision S(N), the sine transform coefficients.
c
      implicit none

      integer n

      double precision a
      double precision angle
      double precision b
      double precision f
      external f
      double precision f2(n)
      double precision fa
      double precision fb
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision s(n)
      double precision x(n)
c
c  Evenly spaced points between A and B, but omitting
c  A and B themselves.
c
      do i = 1, n
        x(i) = ( dble ( n - i + 1 ) * a   
     &         + dble (     i     ) * b ) 
     &         / dble ( n     + 1 )
      end do
c
c  Subtract F1(X) from F(X) to get F2(X).
c
      fa = f ( a )
      fb = f ( b )

      do i = 1, n
        f2(i) = f ( x(i) )                
     &        - ( ( b - x(i)     ) * fa   
     &        +   (     x(i) - a ) * fb ) 
     &        /   ( b        - a )
      end do
c
c  Compute the sine transform of F2(X).
c
      do i = 1, n
        s(i) = 0.0D+00
        do j = 1, n
          angle = pi * dble ( i * j ) / dble ( n + 1 )
          s(i) = s(i) + sin ( angle ) * f2(j)
        end do
        s(i) = s(i) * sqrt ( 2.0D+00 / dble ( n + 1 ) )
      end do

      return
      end
      subroutine sine_transform_interpolant ( n, a, b, fa, fb, s, nx, 
     &  x, value )

c*********************************************************************72
c
cc SINE_TRANSFORM_INTERPOLANT evaluates the sine transform interpolant.
c
c  Discussion:
c
c    The interval [A,B] is divided into N+1 intervals using N+2 points,
c    which are indexed by 0 through N+1.
c
c    The original function F(X) is regarded as the sum of a linear function 
c    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
c    which is 0 at A and B.
c
c    The function F2 has been approximated using the sine transform,
c    and the interpolant is then evaluated as:
c
c      Interp(F)(X) = F1(X) + F2(X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of terms in the approximation.
c
c    Input, double precision A, B, the interval over which the approximant 
c    was defined.
c
c    Input, double precision FA, FB, the function values at A and B.
c
c    Input, double precision S(N), the approximant coefficients.
c
c    Input, integer NX, the number of evaluation points.
c
c    Input, double precision X(NX), the evaluation points.
c
c    Output, double precision VALUE(NX), the value of the interpolant.
c
      implicit none

      integer n
      integer nx

      double precision a
      double precision angle
      double precision b
      double precision f1(nx)
      double precision f2(nx)
      double precision fa
      double precision fb
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision s(n)
      double precision value(nx)
      double precision x(nx)
c
c  Compute linear function F1(X).
c
      do i = 1, nx
        f1(i) = ( ( b - x(i)     ) * fa   
     &          + (     x(i) - a ) * fb ) 
     &          / ( b        - a )
      end do
c
c  Compute sine interpolant F2(X).
c
      do i = 1, nx
        f2(i) = 0.0D+00
        do j = 1, n
          angle = dble ( j ) * ( x(i) - a ) * pi / ( b - a )
          f2(i) = f2(i) + s(j) * sin ( angle )
        end do
      end do

      do i = 1, nx
        f2(i) = f2(i) * sqrt ( 2.0D+00 / dble ( n + 1 ) )
      end do
c
c  Interpolant = F1 + F2.
c
      do i = 1, nx
        value(i) = f1(i) + f2(i)
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
