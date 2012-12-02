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
      subroutine i4vec_max ( n, a, amax )

c*********************************************************************72
c
cc I4VEC_MAX computes the maximum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMAX, the value of the largest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amax
      integer i

      amax = a(1)

      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      subroutine i4vec_mean ( n, a, mean )

c*********************************************************************72
c
cc I4VEC_MEAN returns the mean of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector whose mean is desired.
c
c    Output, double precision MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      double precision mean

      mean = 0.0D+00
      do i = 1, n
        mean = mean + dble ( a(i) )
      end do
      mean = mean / dble ( n )

      return
      end
      subroutine i4vec_min ( n, a, amin )

c*********************************************************************72
c
cc I4VEC_MIN computes the minimum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amin
      integer i

      amin = a(1)

      do i = 2, n
        amin = min ( amin, a(i) )
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
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      subroutine i4vec_variance ( n, a, variance )

c*********************************************************************72
c
cc I4VEC_VARIANCE returns the variance of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector whose variance is desired.
c
c    Output, double precision VARIANCE, the variance of the vector entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      double precision mean
      double precision variance

      if ( n .lt. 2 ) then

        variance = 0.0D+00

      else

        mean = 0.0D+00
        do i = 1, n
          mean = mean + dble ( a(i) )
        end do
        mean = mean / dble ( n )

        variance = 0.0D+00
        do i = 1, n
          variance = variance + ( dble ( a(i) ) - mean )**2
        end do

        variance = variance / dble ( n - 1 )

      end if

      return
      end
      subroutine poisson_fixed_events ( lambda, event_num, seed, t, w )

c*********************************************************************72
c
cc POISSON_FIXED_EVENTS waits for a given number of Poisson events.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision LAMBDA, the average number of events per 
c    unit time.
c
c    Input, integer EVENT_NUM, the number of events to wait for.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision T(0:EVENT_NUM), the time at which a total 
c    of 0, 1, 2, ... and EVENT_NUM events were observed.
c
c    Output, double precision W(0:EVENT_NUM), the waiting time until the
c    I-th event occurred.
c
      implicit none

      integer event_num

      integer i
      double precision lambda
      integer seed
      double precision t(0:event_num)
      double precision w(0:event_num)
c
c  Poisson waiting times follow an exponential distribution.
c
      w(0) = 0.0D+00
      call r8vec_uniform_01 ( event_num, seed, w(1:event_num) )
      do i = 1, event_num
        w(i) = - log ( w(i) ) / lambda
      end do
c
c  The time til event I is the sum of the waiting times 0 through I.
c
      call r8vec_cum ( event_num + 1, w, t )

      return
      end
      subroutine poisson_fixed_time ( lambda, time, seed, n )

c*********************************************************************72
c
cc POISSON_FIXED_TIME counts the Poisson events in a fied time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision LAMBDA, the average number of events 
c    per unit time.
c
c    Input, double precision TIME, the amount of time to observe.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer N, the number of Poisson events observed.
c
      implicit none

      double precision dt
      double precision lambda
      integer n
      double precision r8_uniform_01
      integer seed
      double precision t
      double precision time
      double precision u

      n = 0
      t = 0.0D+00

10    continue

      if ( t .lt. time ) then
        u = r8_uniform_01 ( seed )
        dt = - log ( u ) / lambda
        n = n + 1
        t = t + dt
        go to 10
      end if

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
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
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
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
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
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8vec_cum ( n, a, a_cum )

c*********************************************************************72
c
cc R8VEC_CUM computes the cumulutive sums of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    Input:
c
c      A = (/ 1.0, 2.0, 3.0, 4.0 /)
c
c    Output:
c
c      A_CUM = (/ 1.0, 3.0, 6.0, 10.0 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector to be summed.
c
c    Output, double precision A_CUM(N), the cumulative sums.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_cum(n)
      integer i

      a_cum(1) = a(1)

      do i = 2, n
        a_cum(i) = a_cum(i-1) + a(i)
      end do

      return
      end
      subroutine r8vec_max ( n, a, amax )

c*********************************************************************72
c
cc R8VEC_MAX returns the maximum value in an R8VEC.
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
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMAX, the value of the largest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amax
      integer i

      amax = a(1)
      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      subroutine r8vec_mean ( n, a, mean )

c*********************************************************************72
c
cc R8VEC_MEAN returns the mean of an R8VEC.
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
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector whose mean is desired.
c
c    Output, double precision MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean

      mean = 0.0D+00
      do i = 1, n
        mean = mean + a(i)
      end do
      mean = mean / dble ( n )

      return
      end
      subroutine r8vec_midspace ( n, a, b, x )

c*********************************************************************72
c
cc R8VEC_MIDSPACE creates a vector of linearly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    This function divides the interval [a,b] into n subintervals, and then
c    returns the midpoints of those subintervals.
c
c  Example:
c
c    N = 5, A = 10, B = 20
c    X = [ 11, 13, 15, 17, 19 ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A, B, the endpoints of the interval.
c
c    Output, double precision X(N), a vector of linearly spaced data.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision x(n)

      do i = 1, n
        x(i) = ( dble ( 2 * n - 2 * i + 1 ) * a 
     &         + dble (         2 * i - 1 ) * b ) 
     &         / dble ( 2 * n )
      end do

      return
      end
      subroutine r8vec_min ( n, a, amin )

c*********************************************************************72
c
cc R8VEC_MIN returns the minimum value in an R8VEC.
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
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amin
      integer i

      amin = a(1)
      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
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
      subroutine r8vec_variance ( n, a, variance )

c*********************************************************************72
c
cc R8VEC_VARIANCE returns the variance of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The variance of a vector X of length N is defined as
c
c      mean ( X(1:n) ) = sum ( X(1:n) ) / n
c
c      var ( X(1:n) ) = sum ( ( X(1:n) - mean )^2 ) / ( n - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c    N should be at least 2.
c
c    Input, double precision A(N), the vector.
c
c    Output, double precision VARIANCE, the variance of the vector.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean
      double precision variance

      if ( n .lt. 2 ) then

        variance = 0.0D+00

      else

        mean = 0.0D+00
        do i = 1, n
          mean = mean + a(i)
        end do
        mean = mean / dble ( n )

        variance = 0.0D+00
        do i = 1, n
          variance = variance + ( a(i) - mean )**2
        end do
        variance = variance / dble ( n - 1 )

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
