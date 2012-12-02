      function i8_bit_hi1 ( n )

c*********************************************************************72
c
cc I8_BIT_HI1 returns the position of the high 1 bit base 2 in an I8.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c  Example:
c
c       N    Binary    Hi 1
c    ----    --------  ----
c       0           0     0
c       1           1     1
c       2          10     2
c       3          11     2 
c       4         100     3
c       5         101     3
c       6         110     3
c       7         111     3
c       8        1000     4
c       9        1001     4
c      10        1010     4
c      11        1011     4
c      12        1100     4
c      13        1101     4
c      14        1110     4
c      15        1111     4
c      16       10000     5
c      17       10001     5
c    1023  1111111111    10
c    1024 10000000000    11
c    1025 10000000001    11
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer*8 N, the integer to be measured.
c    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1
c    will always be 0.
c
c    Output, integer*8 I8_BIT_HI1, the number of bits base 2.
c
      implicit none

      integer*8 bit
      integer*8 i8_bit_hi1
      integer*8 i
      integer*8 n

      i = n
      bit = 0

10    continue

        if ( i .le. 0 ) then
          go to 20
        end if

        bit = bit + 1
        i = i / 2

      go to 10

20    continue

      i8_bit_hi1 = bit

      return
      end
      function i8_bit_lo0 ( n )

c*********************************************************************72
c
cc I8_BIT_LO0 returns the position of the low 0 bit base 2 in an I8.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c  Example:
c
c       N    Binary    Lo 0
c    ----    --------  ----
c       0           0     1
c       1           1     2
c       2          10     1
c       3          11     3 
c       4         100     1
c       5         101     2
c       6         110     1
c       7         111     4
c       8        1000     1
c       9        1001     2
c      10        1010     1
c      11        1011     3
c      12        1100     1
c      13        1101     2
c      14        1110     1
c      15        1111     5
c      16       10000     1
c      17       10001     2
c    1023  1111111111     1
c    1024 10000000000     1
c    1025 10000000001     1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer*8 N, the integer to be measured.
c    N should be nonnegative.
c
c    Output, integer*8 I8_BIT_LO0, the position of the low 1 bit.
c
      implicit none

      integer*8 bit
      integer*8 i
      integer*8 i2
      integer*8 i8_bit_lo0
      integer*8 n

      bit = 0
      i = n

10    continue

        bit = bit + 1
        i2 = i / 2

        if ( i .eq. 2 * i2 ) then
          go to 20
        end if

        i = i2

      go to 10

20    continue

      i8_bit_lo0 = bit

      return
      end
      function i8_choose ( n, k )

c*********************************************************************72
c
cc I8_CHOOSE computes the binomial coefficient C(N,K) as an I8.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in integer arithmetic.
c
c    The formula used is:
c
c      C(N,K) = Nc / ( Kc * (N-K)c )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer*8 N, K, are the values of N and K.
c
c    Output, integer*8 I8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer*8 i
      integer*8 i8_choose
      integer*8 k
      integer*8 mn
      integer*8 mx
      integer*8 n
      integer*8 value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0

      else if ( mn .eq. 0 ) then

        value = 1

      else

        mx = max ( k, n - k )
        value = mx + 1

        do i = 2, mn
          value = ( value * ( mx + i ) ) / i
        end do

      end if

      i8_choose = value

      return
      end
      function i8_huge ( )

c*********************************************************************72
c
cc I8_HUGE returns a "huge" I8.
c
c  Discussion:
c
c    On an IEEE 32 bit machine, I8_HUGE should be 2**63 - 1, and its
c    bit pattern should be
c
c     0111111111111111111111111111111111111111111111111111111111111111
c
c    In this case, its numerical value is 9223372036854775807.
c
c    Integer*8 variables might not be available with your compiler.
c
c    The method of defining the literal value with an "_8" suffix
c    might not be acceptable to your compiler.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer*8 I8_HUGE, a "huge" I8.
c
      implicit none

      integer*8 i8_huge

      i8_huge = 9223372036854775807_8

      return
      end
      function i8_huge_normalizer ( )

c*********************************************************************72
c
cc I8_HUGE_NORMALIZER returns the "normalizer" for I8_HUGE.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c    The value returned is 1 / ( I8_HUGE + 1 ).
c
c    For any I8, it should be the case that
c
c     -1 < I8 * I8_HUGE_NORMALIZER < 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision I8_HUGE_NORMALIZER, the "normalizer"
c    for I8_HUGE.
c
      implicit none

      double precision i8_huge_normalizer

      i8_huge_normalizer = 1.084202172485504434007D-19

      return
      end
      function i8_power ( i, j )

c*********************************************************************72
c
cc I8_POWER returns the integer power of an I8.
c
c  Discussion:
c
c    An I8 is a double precision integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer*8 I, J, the base and the power.  
c    J should be nonnegative.
c
c    Output, integer*8 I8_POWER, the value of I^J.
c
      implicit none

      integer*8 i
      integer*8 i8_power
      integer*8 j
      integer*8 k

      if ( j .lt. 0 ) then

        if ( i .eq. 1 ) then
          i8_power = 1
        else if ( i .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I8_POWER - Fatal error!'
          write ( *, '(a)' ) '  I^J requested, with I = 0, J negative.'
          stop
        else
          i8_power = 0
        end if

      else if ( j .eq. 0 ) then

        if ( i .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I8_POWER - Fatal error!'
          write ( *, '(a)' ) '  I^J requested, with I = 0 and J = 0.'
          stop
        else
          i8_power = 1
        end if

      else if ( j .eq. 1 ) then

        i8_power = i

      else

        i8_power = 1
        do k = 1, j
          i8_power = i8_power * i
        end do

      end if

      return
      end
      function i8_uniform ( a, b, seed )

c*********************************************************************72
c
cc I8_UNIFORM returns a scaled pseudorandom I8.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c    Note that ALL integer variables in this routine are
c    of type integer*8!
c
c    Such "double precision integers" might not be available with
c    your compiler.
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
c    22 May 2008
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
c    Input, integer*8 A, B, the limits of the interval.
c
c    Input/output, integer*8 SEED, the "seed" value, which 
c    should NOT be 0.  On output, SEED has been updated.
c
c    Output, integer*8 I8_UNIFORM, a number between A and B.
c
      implicit none

      integer*8 a
      integer*8 b
      integer*8 i8_uniform
      double precision r
      double precision r8i8_uniform_01
      integer*8 seed
      integer*8 value

      if ( seed .eq. 0 ) then
       write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      r = r8i8_uniform_01 ( seed )
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0D+00 - r ) * ( dble ( min ( a, b ) ) - 0.5D+00 ) 
     &  +             r   * ( dble ( max ( a, b ) ) + 0.5D+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i8_uniform = value

      return
      end
      function i8_xor ( i, j )

c*********************************************************************72
c
cc I8_XOR calculates the exclusive OR of two I8's.
c
c  Discussion:
c
c    An I8 is an integer*8 value.
c
c    FORTRAN offers the function IEOR ( I, J ) which should be used instead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c   John Burkardt
c
c  Parameters:
c
c    Input, integer*8 I, J, two values whose exclusive OR is needed.
c
c    Output, integer*8 I8_XOR, the exclusive OR of I and J.
c
      implicit none

      integer*8 i
      integer*8 i1
      integer*8 i2
      integer*8 i8_xor
      integer*8 j
      integer*8 j1
      integer*8 j2
      integer*8 k
      integer*8 l

      i1 = i
      j1 = j
      k = 0
      l = 1

10    continue

        if ( i1 .eq. 0 .and. j1 .eq. 0 ) then
          go to 20
        end if

        i2 = i1 / 2
        j2 = j1 / 2
 
        if ( 
     &    ( ( i1 .eq. 2 * i2 ) .and. ( j1 .ne. 2 * j2 ) ) .or. 
     &    ( ( i1 .ne. 2 * i2 ) .and. ( j1 .eq. 2 * j2 ) ) ) then
          k = k + l
        end if

        i1 = i2
        j1 = j2
        l = 2 * l

      go to 10

20    continue

      i8_xor = k

      return
      end
      function r8i8_uniform ( a, b, seed )

c*********************************************************************72
c
cc R8I8_UNIFORM returns a scaled pseudorandom R8 using an I8 seed.
c
c  Discussion:
c
c    An R8 is a real ( kind = 8 ) value.
c
c    An I8 is an integer*8 value.
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
c    22 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer*8 SEED, the "seed" value, which should
c    NOT be 0.  On output, SEED has been updated.
c
c    Output, double precision R8I8_UNIFORM, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer*8 i8_huge
      integer*8 k
      double precision r8i8_uniform
      integer*8 seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8I8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i8_huge ( )
      end if

      r8i8_uniform = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      return
      end
      function r8i8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8I8_UNIFORM_01 returns a unit pseudorandom R8 using an I8 seed.
c
c  Discussion:
c
c    An R8 is a real ( kind = 8 ) value.
c
c    An I8 is an integer*8 value.
c
c    This routine implements the recursion
c
c      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8I8_UNIFORM_01
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
c    31 May 2007
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
c    Input/output, integer*8 SEED, the "seed" value, which should
c    NOT be 0. On output, SEED has been updated.
c
c    Output, real ( kind = 8 ) R8I8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer*8 i8_huge
      integer*8 k
      double precision r8i8_uniform_01
      integer*8 seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8I8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i8_huge ( )
      end if

      r8i8_uniform_01 = dble ( seed ) * 4.656612875D-10

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
