      subroutine brownian_displacement_display ( k, n, d, t, dsq, 
     &  header )

c*********************************************************************72
c
cc BROWNIAN_DISPLACEMENT_DISPLAY displays average Brownian motion displacement.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer K, the number of repetitions.
c
c    Input, integer N, the number of time steps.  
c
c    Input, double precision D, the diffusion coefficient.
c
c    Input, double precision T, the total time.
c
c    Input, double precision DSQ(K,N), the displacements over time for 
c    each repetition.
c
c    Input, character * ( * ) HEADER, an identifier for the output files.
c
      implicit none

      integer k
      integer n

      character * ( 80 ) command_filename
      integer command_unit
      double precision d
      character * ( 80 ) data_filename
      integer data_unit
      double precision dsq(k,n)
      double precision dsq_ave
      double precision dsq_ideal
      character * ( * ) header
      integer i
      integer i4_uniform
      integer ii(5)
      integer j
      integer seed
      double precision t
      double precision ti

      seed = 123456789
c
c  Choose 5 paths at random.
c
      do j = 1, 5
        ii(j) = i4_uniform ( 1, k, seed )
      end do
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = trim ( header ) // '_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do j = 1, n
        ti = dble ( j - 1 ) * t / dble ( n - 1 )
        dsq_ave = 0.0D+00
        do i = 1, k
          dsq_ave = dsq_ave + dsq(i,j)
        end do
        dsq_ave = dsq_ave / dble ( k )
        dsq_ideal = d * ti
        write ( data_unit, '(8(2x,g14.6))' ) 
     &    ti, dsq(ii(1:5),j), dsq_ave, dsq_ideal
      end do

      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BROWNIAN_DISPLACEMENT data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = trim ( header ) // '_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' // 
     &  trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "' 
     &  // trim ( header ) // '.png"'
      write ( command_unit, '(a)' ) 'set xlabel "T"'
      write ( command_unit, '(a)' ) 'set ylabel "D^2"'
      write ( command_unit, '(a)' ) 'set title "Squared ' //
     &  'displacement (Red), Predicted (Black), Samples (Blue)"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'plot "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:2 title "sample 1" linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:3 title "sample 2" linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:4 title "sample 3" linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:5 title "sample 4" linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:6 title "sample 5" linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:7 title "Averaged" lw 3 linecolor rgb "red", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:8 title "Ideal" lw 3 linecolor rgb "black"'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) 
     &  '  BROWNIAN_DISPLACEMENT plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine brownian_displacement_simulation ( k, n, m, d, t, 
     &  seed, dsq )

c*********************************************************************72
c
cc BROWNIAN_DISPLACEMENT_SIMULATION simulates Brownian displacement.
c
c  Discussion:
c
c    This function computes the square of the distance of the Brownian
c    particle from the starting point, repeating this calculation 
c    several times.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer K, the number of repetitions.
c
c    Input, integer N, the number of time steps to take, plus 1.
c
c    Input, integer M, the spatial dimension.
c
c    Input, double precision D, the diffusion coefficient.  
c    Computationally, this is simply a scale factor between time and space.
c
c    Input, double precision T, the total time.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision DSQ(K,N), the displacements over time for each 
c    repetition.  DSQ(:,1) is 0.0, because we include the displacement at the 
c    initial time. 
c
      implicit none

      integer k
      integer m
      integer n

      double precision d
      double precision dt
      double precision dsq(k,n)
      integer i
      integer i2
      integer j
      double precision s
      integer seed
      double precision t
      double precision x(m,n)

      do i = 1, k

        call brownian_motion_simulation ( m, n, d, t, seed, x )

        do j = 1, n
          s = 0.0D+00
          do i2 = 1, m
            s = s + x(i2,j)**2
          end do
          dsq(i,j) = s
        end do

      end do

      return
      end
      subroutine brownian_motion_display ( m, n, x, header )

c*********************************************************************72
c
cc BROWNIAN_MOTION_DISPLAY displays successive Brownian motion positions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c    M should be 1 or 2.
c
c    Input, integer N, the number of time steps. 
c
c    Input, double precision X(M,N), the particle positions.
c
c    Input, character ( len = * ) HEADER, an identifier for the output files.
c
      implicit none

      integer m
      integer n

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      character * ( * ) header
      integer i
      double precision t
      double precision x(m,n)

      if ( m .ne. 1 .and. m .ne. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BROWNIAN_MOTION_DISPLAY - Fatal error!'
        write ( *, '(a)' ) '  This routine can only handle M = 1 or 2.'
        stop
      end if
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = trim ( header ) // '_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      if ( m .eq. 1 ) then
        do i = 1, n
          t = dble ( i - 1 ) / dble ( n - 1 )
          write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, x(1,i)
        end do
      else if ( m .eq. 2 ) then
        do i = 1, n
          write ( data_unit, '(2x,g14.6,2x,g14.6)' ) x(1,i), x(2,i)
        end do
      end if

      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BROWNIAN_MOTION data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = trim ( header ) // '_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' 
     &  // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "' 
     &  // trim ( header ) // '.png"'
      write ( command_unit, '(a)' ) 'set xlabel "X"'
      write ( command_unit, '(a)' ) 'set ylabel "T"'
      write ( command_unit, '(a)' ) 'set title "Brownian motion in 1D"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) 
     &  // '" using 1:2'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  BROWNIAN_MOTION commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine brownian_motion_simulation ( m, n, d, t, seed, x )

c*********************************************************************72
c
cc BROWNIAN_MOTION_SIMULATION simulates Brownian motion.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of time steps to take, plus 1. 
c
c    Input, double precision D, the diffusion coefficient.  
c
c    Input, double precision T, the total time.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision X(M,N), the initial position at time 0.0, and 
c    the N-1 successive locations of the particle.
c
      implicit none

      integer m
      integer n

      double precision d
      double precision dt
      double precision dx(m)
      integer i
      integer j
      double precision norm_dx
      double precision r8_normal_01
      double precision s
      integer seed
      double precision t
      double precision x(m,n)
c
c  Set the time step.
c
      dt = t / dble ( n - 1 )
c
c  Start at the origin.
c
      do i = 1, m
        x(i,1) = 0.0D+00
      end do
c
c  Take N - 1 steps.
c
      do j = 2, n
c
c  S is the stepsize.
c
        s = sqrt ( d * dt ) * r8_normal_01 ( seed )
c
c  Direction DX is random, unit norm.
c
        if ( m .eq. 1 ) then
          dx(1) = s
        else
          call r8vec_normal_01 ( m, seed, dx )
          norm_dx = 0.0D+00
          do i = 1, m
            norm_dx = norm_dx + dx(i)**2
          end do
          norm_dx = sqrt ( norm_dx )
          do i = 1, m
            dx(i) = s * dx(i) / norm_dx
          end do
        end if
c
c  Add the step to the current position.
c
        do i = 1, m
          x(i,j) = x(i,j-1) + dx(i) 
        end do

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
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
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
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

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
