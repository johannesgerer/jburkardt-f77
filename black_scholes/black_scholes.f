      subroutine asset_path ( s0, mu, sigma, t1, n, seed, s )

c*********************************************************************72
c
cc ASSET_PATH simulates the behavior of an asset price over time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    Black-Scholes for Scientific Computing Students,
c    Computing in Science and Engineering,
c    November/December 2004, Volume 6, Number 6, pages 72-79.
c
c  Parameters:
c
c    Input, double precision S0, the asset price at time 0.
c
c    Input, double precision MU, the expected growth rate.
c
c    Input, double precision SIGMA, the volatility of the asset.
c
c    Input, double precision T1, the expiry date.
c
c    Input, integer N, the number of steps to take between 0 and T1.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision S(0:N), the option values from time 0 to T1 
c    in equal steps.
c
      implicit none

      integer n

      double precision dt
      integer i
      double precision mu
      double precision p
      double precision r(n)
      double precision s(0:n)
      double precision s0
      integer seed
      double precision sigma
      double precision t1

      dt = t1 / dble ( n )

      call r8vec_normal_01 ( n, seed, r )

      s(0) = s0
      p = s0
      do i = 1, n
        p = p * exp ( ( mu - sigma * sigma ) * dt 
     &    + sigma * sqrt ( dt ) * r(i) )
        s(i) = p
      end do

      return
      end
      subroutine binomial ( s0, e, r, sigma, t1, m, c )

c*********************************************************************72
c
cc BINOMIAL uses the binomial method for a European call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    Black-Scholes for Scientific Computing Students,
c    Computing in Science and Engineering,
c    November/December 2004, Volume 6, Number 6, pages 72-79.
c
c  Parameters:
c
c    Input, double precision S0, the asset price at time 0.
c
c    Input, double precision E, the exercise price.
c
c    Input, double precision R, the interest rate.
c
c    Input, double precision SIGMA, the volatility of the asset.
c
c    Input, double precision T1, the expiry date.
c
c    Input, integer M, the number of steps to take 
c    between 0 and T1.
c
c    Output, double precision C, the option value at time 0.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision dt
      double precision e
      integer i
      integer m
      integer n
      double precision p
      double precision r
      double precision s0
      double precision sigma
      double precision t1
      double precision u
      double precision w(1:m+1)
c
c  Time stepsize.
c
      dt = t1 / dble ( m )

      a = 0.5D+00 * ( exp ( - r * dt ) + exp ( ( r + sigma**2 ) * dt ) )

      d = a - sqrt ( a * a - 1.0D+00 )
      u = a + sqrt ( a * a - 1.0D+00 )

      p = ( exp ( r * dt ) - d ) / ( u - d )

      do i = 1, m + 1
        w(i) = max ( s0 * d**(m+1-i) * u**(i-1) - e, 0.0D+00 )
      end do
c
c  Trace backwards to get the option value at time 0.
c
      do n = m, 1, -1
        do i = 1, n
          w(i) = ( 1.0D+00 - p ) * w(i) + p * w(i+1)
        end do
      end do

      do i = 1, m + 1
        w(i) = exp ( - r * t1 ) * w(i)
      end do

      c = w(1)

      return
      end
      subroutine bsf ( s0, t0, e, r, sigma, t1, c )

c*********************************************************************72
c
cc BSF evaluates the Black-Scholes formula for a European call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    Black-Scholes for Scientific Computing Students,
c    Computing in Science and Engineering,
c    November/December 2004, Volume 6, Number 6, pages 72-79.
c
c  Parameters:
c
c    Input, double precision S0, the asset price at time T0.
c
c    Input, double precision T0, the time at which the asset price is known.
c
c    Input, double precision E, the exercise price.
c
c    Input, double precision R, the interest rate.
c
c    Input, double precision SIGMA, the volatility of the asset.
c
c    Input, double precision T1, the expiry date.
c
c    Output, double precision C, the value of the call option.
c
      implicit none

      double precision c
      double precision d1
      double precision d2
      double precision e
      double precision n1
      double precision n2
      double precision r
      double precision s0
      double precision sigma
      double precision t0
      double precision t1
      double precision tau

      tau = t1 - t0

      if ( 0.0D+00 .lt. tau ) then

        d1 = ( log ( s0 / e ) + ( r + 0.5D+00 * sigma * sigma ) * tau )
     &    / ( sigma * sqrt ( tau ) )

        d2 = d1 - sigma * sqrt ( tau )

        n1 = 0.5D+00 * ( 1.0D+00 + erf ( d1 / sqrt ( 2.0D+00 ) ) )
        n2 = 0.5D+00 * ( 1.0D+00 + erf ( d2 / sqrt ( 2.0D+00 ) ) )

        c = s0 * n1 - e * exp ( - r * tau ) * n2

      else

        c = max ( s0 - e, 0.0D+00 )

      end if

      return
      end
      subroutine forward ( e, r, sigma, t1, nx, nt, smax, u )

c*********************************************************************72
c
cc FORWARD uses the forward difference method to value a European call option.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    Black-Scholes for Scientific Computing Students,
c    Computing in Science and Engineering,
c    November/December 2004, Volume 6, Number 6, pages 72-79.
c
c  Parameters:
c
c    Input, double precision E, the exercise price.
c
c    Input, double precision R, the interest rate.
c
c    Input, double precision SIGMA, the volatility of the asset.
c
c    Input, double precision T1, the expiry date.
c
c    Input, integer NX, the number of "space" steps used to 
c    divide the interval [0,L].
c
c    Input, integer NT, the number of time steps.
c
c    Input, double precision SMAX, the maximum value of S to consider.
c
c    Output, double precision U(NX-1,NT+1), the value of the European 
c    call option.
c
      implicit none

      integer nt
      integer nx

      double precision a(2:nx-1)
      double precision b(1:nx-1)
      double precision c(1:nx-2)
      double precision dt
      double precision dx
      double precision e
      integer i
      integer j
      double precision p
      double precision r
      double precision sigma
      double precision smax
      double precision t
      double precision t1
      double precision u(nx-1,nt+1)
      double precision u0

      dt = t1 / dble ( nt )
      dx = smax / dble ( nx )

      do i = 1, nx - 1
        b(i) = 1.0D+00 - r * dt - dt * ( sigma * i )**2
      end do

      do i = 1, nx - 2
        c(i) = 0.5D+00 * dt * ( sigma * i )**2 + 0.5D+00 * dt * r * i
      end do

      do i = 2, nx - 1
        a(i) = 0.5D+00 * dt * ( sigma * i )**2 - 0.5D+00 * dt * r * i
      end do

      u0 = 0.0D+00
      do i = 1, nx - 1
        u0 = u0 + dx
        u(i,1) = max ( u0 - e, 0.0D+00 )
      end do
      
      do j = 1, nt

        t = dble ( j - 1 ) * t1 / dble ( nt )

        p = 0.5D+00 * dt * ( nx - 1 ) 
     &    * ( sigma * sigma * ( nx - 1 ) + r ) 
     &    * ( smax - e * exp ( - r * t ) )

        do i = 1, nx - 1
          u(i,j+1) = b(i) * u(i,j)
        end do

        do i = 1, nx - 2
          u(i,j+1) = u(i,j+1) + c(i) * u(i+1,j)
        end do

        do i = 2, nx - 1
          u(i,j+1) = u(i,j+1) + a(i) * u(i-1,j)
        end do

        u(nx-1,j+1) = u(nx-1,j+1) + p

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
      subroutine mc ( s0, e, r, sigma, t1, m, seed, conf )

c*********************************************************************72
c
cc MC uses Monte Carlo valuation on a European call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    Black-Scholes for Scientific Computing Students,
c    Computing in Science and Engineering,
c    November/December 2004, Volume 6, Number 6, pages 72-79.
c
c  Parameters:
c
c    Input, double precision S0, the asset price at time 0.
c
c    Input, double precision E, the exercise price.
c
c    Input, double precision R, the interest rate.
c
c    Input, double precision SIGMA, the volatility of the asset.
c
c    Input, double precision T1, the expiry date.
c
c    Input, integer M, the number of simulations.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision CONF(2), the estimated range of the valuation.
c
      implicit none

      integer m

      double precision conf(2)
      double precision e
      integer i
      double precision pmean
      double precision pvals(m)
      double precision r
      double precision s0
      integer seed
      double precision sigma
      double precision std
      double precision svals(m)
      double precision t1
      double precision u(m)
      double precision width

      call r8vec_normal_01 ( m, seed, u )

      do i = 1, m
        svals(i) = s0 * exp ( ( r - 0.5D+00 * sigma * sigma ) * t1 
     &    + sigma * sqrt ( t1 ) * u(i) )
      end do

      do i = 1, m
        pvals(i) = exp ( - r * t1 ) * max ( svals(i) - e, 0.0D+00 )
      end do

      pmean = 0.0D+00
      do i = 1, m
        pmean = pmean + pvals(i)
      end do
      pmean = pmean / dble ( m )

      std = 0.0D+00
      do i = 1, m
        std = std + ( pvals(i) - pmean ) ** 2
      end do
      std = sqrt ( std / dble ( m - 1 ) )

      width = 1.96D+00 * std / sqrt ( dble ( m ) )

      conf(1) = pmean - width
      conf(2) = pmean + width

      return
      end
      function r8_normal_01 ( seed )

c*********************************************************************72
c
cc R8_NORMAL_01 returns a unit pseudonormal R8.
c
c  Discussion:
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, the code can use the second
c    value that it calculated.
c
c    However, if the user has changed the SEED value between calls,
c    the routine automatically resets itself and discards the saved data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_uniform_01
      integer seed
      integer seed1
      integer seed2
      integer seed3
      integer used
      double precision v1
      double precision v2

      save seed1
      save seed2
      save seed3
      save used
      save v2

      data seed2 / 0 /
      data used / 0 /
      data v2 / 0.0D+00 /
c
c  If USED is odd, but the input SEED does not match
c  the output SEED on the previous call, then the user has changed
c  the seed.  Wipe out internal memory.
c
      if ( mod ( used, 2 ) .eq. 1 ) then

        if ( seed .ne. seed2 ) then
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        end if

      end if
c
c  If USED is even, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        seed1 = seed

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed

        r2 = r8_uniform_01 ( seed )

        seed3 = seed

        v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

        r8_normal_01 = v1
        seed = seed2
c
c  If USED is odd (and the input SEED matched the output value from
c  the previous call), return the second normal and its corresponding seed.
c
      else

        r8_normal_01 = v2
        seed = seed3

      end if

      used = used + 1

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
      subroutine r8vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
c
c    The Box-Muller method is used, which is efficient, but
c    generates an even number of values each time.  On any call
c    to this routine, an even number of new values are generated.
c    Depending on the situation, one value may be left over.
c    In that case, it is saved for the next call.
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
c  Parameters:
c
c    Input, integer N, the number of values desired.  If N is negative,
c    then the code will flush its internal memory; in particular,
c    if there is a saved value to be used on the next call, it is
c    instead discarded.  This is useful if the user has reset the
c    random number seed, for instance.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision X(N), a sample of the standard normal PDF.
c
c  Local parameters:
c
c    Local, integer MADE, records the number of values that have
c    been computed.  On input with negative N, this value overwrites
c    the return value of N, so the user can get an accounting of
c    how much work has been done.
c
c    Local, integer SAVED, is 0 or 1 depending on whether there is a
c    single saved value left over from the previous call.
c
c    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
c    X that we need to compute.  This starts off as 1:N, but is adjusted
c    if we have a saved value that can be immediately stored in X(1),
c    and so on.
c
c    Local, double precision Y, the value saved from the previous call, if
c    SAVED is 1.
c
      implicit none

      integer n

      integer i
      integer m
      integer made
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(2)
      double precision r8_uniform_01
      integer saved
      integer seed
      double precision x(n)
      integer x_hi_index
      integer x_lo_index
      double precision y

      save made
      save saved
      save y

      data made / 0 /
      data saved / 0 /
      data y / 0.0D+00 /
c
c  I'd like to allow the user to reset the internal data.
c  But this won't work properly if we have a saved value Y.
c  I'm making a crock option that allows the user to signal
c  explicitly that any internal memory should be flushed,
c  by passing in a negative value for N.
c
      if ( n .lt. 0 ) then
        n = made
        made = 0
        saved = 0
        y = 0.0D+00
        return
      else if ( n .eq. 0 ) then
        return
      end if
c
c  Record the range of X we need to fill in.
c
      x_lo_index = 1
      x_hi_index = n
c
c  Use up the old value, if we have it.
c
      if ( saved .eq. 1 ) then
        x(1) = y
        saved = 0
        x_lo_index = 2
      end if
c
c  Maybe we don't need any more values.
c
      if ( x_hi_index - x_lo_index + 1 .eq. 0 ) then
c
c  If we need just one new value, do that here to avoid null arrays.
c
      else if ( x_hi_index - x_lo_index + 1 .eq. 1 ) then

        r(1) = r8_uniform_01 ( seed )

        if ( r(1) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal errorc'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        r(2) = r8_uniform_01 ( seed )

        x(x_hi_index) =
     &           sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * cos ( 2.0D+00 * pi * r(2) )
        y =      sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + 2
c
c  If we require an even number of values, that's easy.
c
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) .eq. 0 ) then

        do i = x_lo_index, x_hi_index, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        made = made + x_hi_index - x_lo_index + 1
c
c  If we require an odd number of values, we generate an even number,
c  and handle the last pair specially, storing one in X(N), and
c  saving the other for later.
c
      else

        do i = x_lo_index, x_hi_index - 1, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        call r8vec_uniform_01 ( 2, seed, r )

        x(n) = sqrt ( -2.0D+00 * log ( r(1) ) )
     &    * cos ( 2.0D+00 * pi * r(1) )

        y = sqrt ( -2.0D+00 * log ( r(2) ) )
     &    * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + x_hi_index - x_lo_index + 2

      end if

      return
      end
      subroutine r8vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_PART prints "part" of an R8VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer max_print
      character*(*) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '  ........  ..............'
        i = n

        write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(2x,i8,a1,1x,g14.6,a)' )
     &    i, ':', a(i), '...more entries...'

      end if

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
      subroutine r8vec_write ( output_filename, n, x )

c*********************************************************************72
c
cc R8VEC_WRITE writes an R8VEC file.
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
c    10 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      double precision x(n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. n ) then
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, '(2x,g24.16)' ) x(j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

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
