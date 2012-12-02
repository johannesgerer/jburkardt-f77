      subroutine cdelay2 ( m, q )

c*********************************************************************72
c
cc CDELAY2 is a circular buffer implementation of M-fold delay.
c
c  Example:
c
c    Suppose we call CDELAY2 12 times, always with M = 3, and with
c    Q having the input value 3 on the first call.  Q will go through
c    the following sequence of values over the 12 calls:
c
c    I   M  Qin  Qout
c
c    1   3   3   2
c    2   3   2   1
c    3   3   1   0
c    4   3   0   3
c    5   3   3   2
c    6   3   2   1
c    7   3   1   0
c    8   3   0   3
c    9   3   3   2
c   10   3   2   1
c   11   3   1   0
c   12   3   0   3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    Original C version by Sophocles Orfanidis.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer M, the maximum value that Q can have.
c
c    Input/output, integer Q, a counter which is decremented 
c    on every call.
c    However, the value "after" 0 is M.  
c
      implicit none

      integer m
      integer q
c
c  Decrement the offset.
c
      q = q - 1
c
c  Q = - 1 wraps to Q = M.
c
      call wrap2 ( m, q )

      return
      end
      subroutine corr ( n, x, m, r )

c*********************************************************************72
c
cc CORR computes the sample correlation of a signal sample.
c
c  Discussion:
c
c    The sample correlation is defined, for 0 <= i < N, as
c
c      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * X(j)
c
c    The sample correlation is an estimate of the correlation function.
c
c    It is usually the case that the signal X is assumed to
c    have zero mean.  Here, we compute the mean and adjust the
c    calculation accordingly:
c
c      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) 
c        ( X(i+j) - Xbar ) * ( X(j) - Xbar )
c
c    Experience suggests that only the first 5 or 10 percent of
c    the lags are statistically reliable, so that one might choose
c    M = N / 20 or M = N / 10, for instance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer N, the number of equally spaced signal
c    samples.
c
c    Input, double precision X(0:N-1), the signal samples.
c
c    Input, integer M, the maximum lag to consider.
c    0 <= M < N.
c
c    Output, double precision R(0:M), the sample correlations.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      double precision r(0:m)
      double precision x(n)
      double precision xbar

      do i = 0, m
        r(i) = 0.0D+00
      end do

      xbar = 0.0D+00
      do j = 0, n - 1
        xbar = xbar + x(j)
      end do
      xbar = xbar / dble ( n )

      do i = 0, m
        do j = 0, n - i - 1
          r(i) = r(i) + ( x(i+j) - xbar ) * ( x(j) - xbar )
        end do
      end do

      do i = 0, m
        r(i) = r(i) / dble ( n )
      end do

      return
      end
      subroutine cross_corr ( n, x, y, m, r )

c*********************************************************************72
c
cc CROSS_CORR computes the sample cross correlation between two signal samples.
c
c  Discussion:
c
c    The sample cross correlation is defined, for 0 <= i < N, as
c
c      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * Y(j)
c
c    The sample cross correlation is an estimate of the cross 
c    correlation function.
c
c    It is usually the case that the signals X and Y are assumed to
c    have zero mean.  Here, we compute the means and adjust the
c    calculation accordingly:
c
c      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) 
c        ( X(i+j) - Xbar ) * ( Y(j) - Ybar )
c
c    Experience suggests that only the first 5 or 10 percent of
c    the lags are statistically reliable, so that one might choose
c    M = N / 20 or M = N / 10, for instance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer N, the number of equally spaced signal
c    samples.
c
c    Input, double precision X(0:N-1), Y(0:N-1), the signal samples.
c
c    Input, integer M, the maximum lag to consider.
c    0 <= M < N.
c
c    Output, double precision R(0:M), the sample correlations.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      double precision r(0:m)
      double precision x(0:n-1)
      double precision xbar
      double precision y(0:n-1)
      double precision ybar

      do i = 0, m
        r(i) = 0.0D+00
      end do

      xbar = 0.0D+00
      do j = 0, n - 1
        xbar = xbar + x(j)
      end do
      xbar = xbar / dble ( n )

      ybar = 0.0D+00
      do j = 0, n - 1
        ybar = ybar + y(j)
      end do
      ybar = ybar / dble ( n )

      do i = 0, m
        do j = 0, n - i - 1
          r(i) = r(i) + ( x(i+j) - xbar ) * ( y(j) - ybar )
        end do
      end do

      do i = 0, m
        r(i) = r(i) / dble ( n )
      end do

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
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
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
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
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

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
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
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      function ran1f ( b, u, q, seed )

c*********************************************************************72
c
cc RAN1F is a 1/F random number generator.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    Original C version by Sophocles Orfanidis.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer B, the number of signals to combine.
c    For this algorithm, B cannot be more than 31c
c
c    Input/output, double precision U(B), the signals to combine.  It is 
c    expected that each of the initial values of U will be drawn from a 
c    distribution with zero mean.
c
c    Input/output, integer Q(B), a set of counters that determine 
c    when each entry of U is to be updated.
c
c    Input/output, integr SEED, a seed for the random number generator.
c
c    Output, double precision RAN1F, the value.
c
      implicit none

      integer b

      integer i
      integer j
      integer q(b)
      double precision ran1f
      double precision ranh
      integer seed
      double precision u(b)
      double precision y

      if ( 31 .lt. b ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAN1F - Fatal error!'
        write ( *, '(a)' ) '  32 <= B, too many signals.'
        stop
      end if

      y = 0.0D+00

      j = 1
      do i = 1, b
        y = y + ranh ( j, u(i), q(i), seed )
        j = j * 2
      end do

      if ( 0 .lt. b ) then
        y = y / dble ( b )
      end if

      ran1f = y

      return
      end
      function ranh ( d, u, q, seed )

c*********************************************************************72
c
c  Purpose:
c
c    RANH is a hold random number generator of period D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    Original C version by Sophocles Orfanidis.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer D, the hold period.  D must be at least 1.
c
c    Input/output, double precision U, a value to be held until Q has 
c    decremented to 0, when Q will be reset to D, and U will be randomly reset.
c
c    Input/output, integer Q, a counter which is decremented by 1 
c    on each call until reaching 0.
c
c    Input/output, integr SEED, a seed for the random number generator.
c
c    Output, double precision RANH, the input value of U.
c
      implicit none

      integer d
      integer q
      double precision r8_uniform_01
      double precision ranh
      integer seed
      double precision u
      double precision y

      if ( d .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANH - Fatal error!'
        write ( *, '(a)' ) '  D < 1.'
        stop
      end if
c
c  Hold this sample for D calls.
c
      y = u
c
c  Decrement Q and wrap mod D.
c
      call cdelay2 ( d - 1, q )
c
c  Every D calls, get a new U with zero mean.
c
      if ( q .eq. 0 ) then
        u = r8_uniform_01 ( seed )
        u = 2.0D+00 * u - 1.0D+00
      end if

      ranh = y

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
      subroutine wrap2 ( m, q )

c*****************************************************************************80
c
cc WRAP2 is a circular wrap of the pointer offset Q.
c
c  Discussion:
c
c    Input values of Q between 0 and M are "legal".
c    Values of Q below 0 are incremented by M + 1 until they are legal.
c    Values of Q above M are decremented by M + 1 until they become legal.
c    The legal value is the output value of the function.
c
c  Example:
c
c    M  Qin  Qout
c
c    3  -5   3
c    3  -4   0
c    3  -3   1
c    3  -2   2
c    3  -1   3
c    3   0   0
c    3   1   1
c    3   2   2
c    3   3   3
c    3   4   0
c    3   5   1
c    3   6   2
c    3   7   3
c    3   8   0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    Original C version by Sophocles Orfanidis.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Sophocles Orfanidis,
c    Introduction to Signal Processing,
c    Prentice-Hall, 1995,
c    ISBN: 0-13-209172-0,
c    LC: TK5102.5.O246.
c
c  Parameters:
c
c    Input, integer M, the maximum acceptable value for outputs.
c    M must be at least 0.
c
c    Input/output, integer Q, the value to be wrapped.
c
      implicit none

      integer m
      integer q

      if ( m .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRAP2 - Fatal error!'
        write ( *, '(a)' ) '  M < 0.'
        stop
      end if
c
c  When Q = M + 1, it wraps to Q = 0.
c
10    continue

      if ( m .lt. q ) then
        q = q - m - 1
        go to 10
      end if
c
c  When Q = - 1, it wraps to Q = M.
c
20    continue

      if ( q .lt. 0 ) then
        q = q + m + 1
        go to 20
      end if

      return
      end
