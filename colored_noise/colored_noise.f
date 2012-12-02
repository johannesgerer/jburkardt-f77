      subroutine f_alpha ( n, q_d, alpha, seed, x )

c*********************************************************************72
c
cc F_ALPHA generates a 1/F^ALPHA noise sequence.
c
c  Discussion:
c
c    Thanks to Miro Stoyanov for pointing out that the second half of
c    the data returned by the inverse Fourier transform should be
c    discarded, 24 August 2010.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    Original C version by Todd Walter.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, integer N, the number of samples to generate.
c
c    Input, double precision Q_D, the variance of the noise.
c
c    Input, double precision ALPHA, the exponent for the noise.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision X(N), a sequence sampled with the given
c    power law.
c
      implicit none

      integer n

      double precision alpha
      double precision h_a(n)
      double precision h_azero
      double precision h_b(n)
      double precision hfa(2*n)
      integer i
      double precision q_d
      double precision r8_normal_01
      integer seed
      double precision w_a(n)
      double precision w_azero
      double precision w_b(n)
      double precision wfa(2*n)
      double precision wi
      double precision wr
      double precision x(n)
      double precision x2(2*n)
c
c  Set the deviation of the noise.
c
      q_d = sqrt ( q_d )
c
c  Generate the coefficients Hk.
c
      hfa(1) = 1.0D+00
      do i = 2, n
        hfa(i) = hfa(i-1) * ( 0.5D+00 * alpha + dble ( i - 2 ) ) 
     &    / ( dble ( i - 1 ) )
      end do
      do i = n + 1, 2 * n
        hfa(i) = 0.0D+00
      end do
c
c  Fill Wk with white noise.
c
      do i = 1, n
        wfa(i) = q_d * r8_normal_01 ( seed )
      end do
      do i = n + 1, 2 * n
        wfa(i) = 0.0D+00
      end do
c
c  Perform the discrete Fourier transforms of Hk and Wk.
c
      call r8vec_sftf ( 2 * n, hfa, h_azero, h_a, h_b )

      call r8vec_sftf ( 2 * n, wfa, w_azero, w_a, w_b )
c
c  Multiply the two complex vectors.
c
      w_azero = w_azero * h_azero

      do i = 1, n
        wr = w_a(i)
        wi = w_b(i)
        w_a(i) = wr * h_a(i) - wi * h_b(i)
        w_b(i) = wi * h_a(i) + wr * h_b(i)
      end do
c
c  This scaling is introduced only to match the behavior
c  of the Numerical Recipes code...
c
      w_azero = w_azero * dble ( 2 * n )

      do i = 1, n - 1
        w_a(i) = w_a(i) * dble ( n )
        w_b(i) = w_b(i) * dble ( n )
      end do

      w_a(n) = w_a(n) * dble ( 2 * n )
      w_b(n) = w_b(n) * dble ( 2 * n )
c
c  Take the inverse Fourier transform of the result.
c
      call r8vec_sftb ( 2 * n, w_azero, w_a, w_b, x2 )
c
c  Only return the first N inverse Fourier transform values.
c
      do i = 1, n
        x(i) = x2(i)
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
      if ( mod ( used, 2 ) == 1 ) then

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
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
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
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

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
c    16 September 2003
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
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '......  ..............'
        i = n

        write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(i6,a,1x,g14.6,a)' ) 
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      subroutine r8vec_sftb ( n, azero, a, b, r )

c*********************************************************************72
c
cc R8VEC_SFTB computes a "slow" backward Fourier transform of real data.
c
c  Discussion:
c
c    SFTB and SFTF are inverses of each other.  If we begin with data
c    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
c    resulting R vector, we should get back the original AZERO, A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, double precision AZERO, the constant Fourier coefficient.
c
c    Input, double precision A(N/2), B(N/2), the Fourier coefficients.
c
c    Output, double precision R(N), the reconstructed data sequence.
c
      implicit none

      integer n

      double precision a(n/2)
      double precision azero
      double precision b(n/2)
      integer i
      integer k
      double precision r(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta

      do i = 1, n
        r(i) = azero
      end do

      do i = 1, n
        do k = 1, n / 2
          theta = dble ( k * ( i - 1 ) * 2 ) * pi / dble ( n )
          r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
        end do
      end do

      return
      end
      subroutine r8vec_sftf ( n, r, azero, a, b )

c*********************************************************************72
c
cc R8VEC_SFTF computes a "slow" forward Fourier transform of real data.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
c    A, and B, we should get back the original R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision R(N), the data to be transformed.
c
c    Output, double precision AZERO, = sum ( 1 <= I <= N ) R(I) / N.
c
c    Output, double precision A(N/2), B(N/2), the Fourier coefficients.
c
      implicit none

      integer n

      double precision a(1:n/2)
      double precision azero
      double precision b(1:n/2)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(n)
      double precision theta

      azero = 0.0D+00
      do i = 1, n
        azero = azero + r(i)
      end do
      azero = azero / dble ( n )

      do i = 1, n / 2

        a(i) = 0.0D+00
        b(i) = 0.0D+00

        do j = 1, n
          theta = dble ( 2 * i * ( j - 1 ) ) * pi 
     &      / dble ( n )
          a(i) = a(i) + r(j) * cos ( theta )
          b(i) = b(i) + r(j) * sin ( theta )
        end do

        a(i) = a(i) / dble ( n )
        b(i) = b(i) / dble ( n )

        if ( i .ne. ( n / 2 ) ) then
          a(i) = 2.0D+00 * a(i)
          b(i) = 2.0D+00 * b(i)
        end if

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
