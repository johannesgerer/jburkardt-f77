      function i4_bit_hi1 ( n )

c*********************************************************************72
c
cc I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
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
c    22 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be measured.
c    N should be nonnegative.  If N is nonpositive, the function
c    will always be 0.
c
c    Output, integer I4_BIT_HI1, the position of the highest bit.
c
      implicit none

      integer bit
      integer i
      integer i4_bit_hi1
      integer n

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

      i4_bit_hi1 = bit

      return
      end
      function i4_bit_lo0 ( n )

c*********************************************************************72
c
cc I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
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
c    16 February 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be measured.
c    N should be nonnegative.
c
c    Output, integer I4_BIT_LO0, the position of the low 1 bit.
c
      implicit none

      integer bit
      integer i
      integer i2
      integer i4_bit_lo0
      integer n

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

      i4_bit_lo0 = bit

      return
      end
      function i4_bit_lo1 ( n )

c*********************************************************************72
c
cc I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Example:
c
c       N    Binary    Lo 1
c    ----    --------  ----
c       0           0     0
c       1           1     1
c       2          10     2
c       3          11     1
c       4         100     3
c       5         101     1
c       6         110     2
c       7         111     1
c       8        1000     4
c       9        1001     1
c      10        1010     2
c      11        1011     1
c      12        1100     3
c      13        1101     1
c      14        1110     2
c      15        1111     1
c      16       10000     5
c      17       10001     1
c    1023  1111111111     1
c    1024 10000000000    11
c    1025 10000000001     1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 February 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be measured.
c    N should be nonnegative.
c
c    Output, integer I4_BIT_LO1, the position of the low 1 bit.
c
      implicit none

      integer bit
      integer i
      integer i2
      integer i4_bit_lo1
      integer n

      bit = 0
      i = n

10    continue

        bit = bit + 1
        i2 = i / 2

        if ( i /= 2 * i2 ) then
          go to 20
        end if

        i = i2

      go to 10

20    continue

      i4_bit_lo1 = bit

      return
      end
      function i4_bit_reverse ( i, n )

c*********************************************************************72
c
cc I4_BIT_REVERSE reverses the bits in an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Example:
c
c       I      N  2^N     I4_BIT_REVERSE ( I, N )
c    ----    --------  -----------------------
c       0      0    1     0
c       1      0    1     1
c
c       0      3    8     0
c       1      3    8     4
c       2      3    8     2
c       3      3    8     6
c       4      3    8     1
c       5      3    8     5
c       6      3    8     3
c       7      3    8     7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be bit reversed.
c    I should be nonnegative.  Normally I < 2^N.
c
c    Input, integer N, indicates the number of bits to
c    be reverse (N+1) or the base with respect to which the integer is to
c    be reversed (2^N).  N should be nonnegative.
c
c    Output, integer I4_BIT_REVERSE, the bit reversed value.
c
      implicit none

      integer b
      integer i
      integer i4_bit_reverse
      integer j
      integer n
      integer value

      if ( i .lt. 0 ) then

        value = -1

      else if ( n .lt. 0 ) then

        value = -1

      else

        b = 2**n
        j = mod ( i, b )

        value = 0

10      continue

          if ( b .eq. 1 ) then

            value = value + j
            j = 0
            go to 20

          else

            if ( mod ( j, 2 ) .eq. 1 ) then
              value = value + b / 2
              j = j - 1
            end if

            j = j / 2
            b = b / 2

          end if

        go to 10

20      continue

      end if

      i4_bit_reverse = value

      return
      end
      function i4_ceiling ( r )

c*********************************************************************72
c
cc I4_CEILING rounds an R8 "up" to the nearest I4.
c
c  Example:
c
c     R     Value
c
c    -1.1  -1
c    -1.0  -1
c    -0.9   0
c     0.0   0
c     5.0   5
c     5.1   6
c     5.9   6
c     6.0   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the value to be rounded up.
c
c    Output, integer I4_CEILING, the rounded value.
c
      implicit none

      double precision r
      integer i4_ceiling
      integer value

      value = int ( r )
      if ( dble ( value ) .lt. r ) then
        value = value + 1
      end if

      i4_ceiling = value

      return
      end
      function i4_characteristic ( q )

c*********************************************************************72
c
cc I4_CHARACTERISTIC gives the characteristic for an I4.
c
c  Discussion:
c
c    For any positive integer Q, the characteristic is:
c
c    Q, if Q is a prime;
c    P, if Q = P**N for some prime P and some integer N;
c    0, otherwise, that is, if Q is negative, 0, 1, or the product
c       of more than one distinct prime.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 December 2004
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Harald Niederreiter,
c    Algorithm 738:
c    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
c    ACM Transactions on Mathematical Software,
c    Volume 20, Number 4, 1994, pages 494-495.
c
c  Parameters:
c
c    Input, integer Q, the value to be tested.
c
c    Output, integer I4_CHARACTERISTIC, the characteristic of Q.
c
      implicit none

      integer i
      integer i4_characteristic
      integer i_max
      integer q
      integer q_copy

      if ( q .le. 1 ) then
        i4_characteristic = 0
        return
      end if
c
c  If Q is not prime, then there is at least one prime factor
c  of Q no greater than SQRT(Q)+1.
c
c  A faster code would only consider prime values of I,
c  but that entails storing a table of primes and limiting the
c  size of Q.  Simplicity and flexibility for nowc
c
      i_max = int ( sqrt ( real ( q ) ) ) + 1
      q_copy = q

      do i = 2, i_max

        if ( mod ( q_copy, i ) .eq. 0 ) then

10        continue

          if ( mod ( q_copy, i ) .eq. 0 ) then
            q_copy = q_copy / i
            go to 10
          end if

          if ( q_copy .eq. 1 ) then
            i4_characteristic = i
          else
            i4_characteristic = 0
          end if

          return

        end if

      end do
c
c  If no factor was found, then Q is prime.
c
      i4_characteristic = q

      return
      end
      function i4_choose ( n, k )

c*********************************************************************72
c
cc I4_CHOOSE computes the binomial coefficient C(N,K).
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in integer arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2007
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
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer I4_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer i4_choose
      integer k
      integer mn
      integer mx
      integer n
      integer value

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

      i4_choose = value

      return
      end
      function i4_div_rounded ( a, b )

c*********************************************************************72
c
cc I4_DIV_ROUNDED computes the rounded result of I4 division.
c
c  Discussion:
c
c    This routine computes C = A / B, where A, B and C are integers
c    and C is the closest integer value to the exact real result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, B, the number to be divided,
c    and the divisor.
c
c    Output, integer I4_DIV_ROUNDED, the rounded result
c    of the division.
c
      implicit none

      integer a
      integer a_abs
      integer b
      integer b_abs
      integer i4_div_rounded
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer value

      if ( a .eq. 0 .and. b .eq. 0 ) then

        value = i4_huge
 
      else if ( a .eq. 0 ) then

        value = 0

      else if ( b .eq. 0 ) then

        if ( a .lt. 0 ) then
          value = - i4_huge
        else
          value = + i4_huge
        end if

      else

        a_abs = abs ( a )
        b_abs = abs ( b )

        value = a_abs / b_abs
c
c  Round the value.
c
        if ( ( 2 * value + 1 ) * b_abs .lt. 2 * a_abs ) then
          value = value + 1
        end if
c
c  Set the sign.
c
        if ( ( a .lt. 0 .and. 0 .lt. b ) .or. 
     &       ( 0 .lt. a .and. b .lt. 0 ) ) then
          value = - value
        end if

      end if

      i4_div_rounded = value

      return
      end
      function i4_divp ( i, j )

c*********************************************************************72
c
cc I4_DIVP returns the smallest multiple of J greater than or equal to I.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Example:
c
c    I  J  I4_DIVP(I,J)
c
c    0  4    0
c    1  4    1
c    2  4    1
c    3  4    1
c    4  4    1
c    5  4    2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be analyzed.
c
c    Input, integer J, the number, multiples of which will
c    be compared against I.  J may not be zero.
c
c    Output, integer I4_DIVP, the smallest multiple of J that
c    is greater than or equal to I.
c
      implicit none

      integer i
      integer i4_divp
      integer j

      if ( j .ne. 0 ) then
        i4_divp = 1 + ( i - 1 ) / j
      else
        i4_divp = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_DIVP - Fatal error!'
        write ( *, '(a)' ) '  The input value of J was zero.'
        stop
      end if

      return
      end
      function i4_even ( i )

c*********************************************************************72
c
cc I4_EVEN returns TRUE if an I4 is even.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be tested.
c
c    Output, logical I4_EVEN, is TRUE if I is even.
c
      implicit none

      integer i
      logical i4_even
      integer :: i4_two = 2

      i4_even = ( mod ( i, i4_two ) .eq. 0 )

      return
      end
      function i4_factorial ( n )

c*********************************************************************72
c
cc I4_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c    0 <= N <= 13 is required.
c
c    Output, integer I4_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer i4_factorial
      integer n

      i4_factorial = 1

      if ( 13 .lt. n ) then
        i4_factorial = - 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
        write ( *, '(a)' )
     &  '  I4_FACTORIAL(N) cannot be computed as an integer'
        write ( *, '(a)' ) '  for 13 < N.'
        write ( *, '(a,i8)' ) '  Input value N = ', n
        stop
      end if

      do i = 1, n
        i4_factorial = i4_factorial * i
      end do

      return
      end
      function i4_factorial2 ( n )

c*********************************************************************72
c
cc I4_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial 
c    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
c
c    Output, integer I4_FACTORIAL2, the value of N!!.
c
      implicit none

      integer i4_factorial2
      integer n
      integer n_copy

      if ( n .lt. 1 ) then
        i4_factorial2 = 1
        return
      end if

      n_copy = n
      i4_factorial2 = 1

10    continue

      if ( 1 .lt. n_copy ) then
        i4_factorial2 = i4_factorial2 * n_copy
        n_copy = n_copy - 2
        go to 10
      end if

      return
      end
      function i4_floor ( r )

c*********************************************************************72
c
cc I4_FLOOR rounds an R8 "down" (towards -infinity) to the nearest I4.
c
c  Example:
c
c    R     Value
c
c   -1.1  -2
c   -1.0  -1
c   -0.9  -1
c    0.0   0
c    5.0   5
c    5.1   5
c    5.9   5
c    6.0   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the value to be rounded down.
c
c    Output, integer I4_FLOOR, the rounded value.
c
      implicit none

      integer i4_floor
      double precision r
      integer value

      value = int ( r )
      if ( r .lt. dble ( value ) ) then
        value = value - 1
      end if

      i4_floor = value

      return
      end
      subroutine i4_fraction ( i, j, k )

c*********************************************************************72
c
cc I4_FRACTION computes a ratio and returns an integer result.
c
c  Discussion:
c
c    Given integer variables I and J, FORTRAN will evaluate the expression
c    "I/J" using integer arithmetic.  This routine, which carries out the
c    same operation, is thus not needed in FORTRAN.  It is provided simply
c    to match the corresponding function in MATLAB, where the default
c    result of "I/J" is a real number.
c
c  Example:
c
c       I     J     Real     K = I4_FRACTION ( I, J)
c
c       1     2     0.5      0
c       8     4     2.00     2
c       9     4     2.25     2
c       7     4     1.75     1
c      -7     4    -1.75    -1
c       7    -4    -1.75    -1
c      -7    -4     1.75     1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the arguments.
c
c    Output, integer K, the value of the ratio.
c
      implicit none

      integer i
      integer j
      integer k

      k = i / j

      return
      end
      function i4_gcd ( i, j )

c*********************************************************************72
c
cc I4_GCD finds the greatest common divisor of I and J.
c
c  Discussion:
c
c    Only the absolute values of I and J are
c    considered, so that the result is always nonnegative.
c
c    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
c
c    If I and J have no common factor, I4_GCD is returned as 1.
c
c    Otherwise, using the Euclidean algorithm, I4_GCD is the
c    largest common factor of I and J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common divisor
c    is desired.
c
c    Output, integer I4_GCD, the greatest common divisor of I and J.
c
      implicit none

      integer i
      integer i4_gcd
      integer ip
      integer iq
      integer ir
      integer j

      i4_gcd = 1
c
c  Return immediately if either I or J is zero.
c
      if ( i .eq. 0 ) then
        i4_gcd = max ( 1, abs ( j ) )
        return
      else if ( j .eq. 0 ) then
        i4_gcd = max ( 1, abs ( i ) )
        return
      end if
c
c  Set IP to the larger of I and J, IQ to the smaller.
c  This way, we can alter IP and IQ as we go.
c
      ip = max ( abs ( i ), abs ( j ) )
      iq = min ( abs ( i ), abs ( j ) )
c
c  Carry out the Euclidean algorithm.
c
10    continue

        ir = mod ( ip, iq )

        if ( ir .eq. 0 ) then
          go to 20
        end if

        ip = iq
        iq = ir

      go to 10

20    continue

      i4_gcd = iq

      return
      end
      function i4_gcdb ( i, j, k )

c*********************************************************************72
c
cc I4_GCDB finds the greatest common divisor of the form K^N of two I4's.
c
c  Discussion:
c
c    Note that if J is negative, I4_GCDB will also be negative.
c    This is because it is likely that the caller is forming
c    the fraction I/J, and so any minus sign should be
c    factored out of J.
c
c    If I and J are both zero, I4_GCDB is returned as 1.
c
c    If I is zero and J is not, I4_GCDB is returned as J,
c    and vice versa.
c
c    If I and J are nonzero, and have no common divisor of the
c    form K^N, I4_GCDB is returned as 1.
c
c    Otherwise, I4_GCDB is returned as the largest common divisor
c    of the form K^N shared by I and J.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common
c    divisor K^N is desired.
c
c    Input, integer K, the possible divisor of I and J.
c
c    Output, integer I4_GCDB, the greatest common divisor of
c    the form K^N shared by I and J.
c
      implicit none

      integer i
      integer icopy
      integer i4_gcdb
      integer j
      integer jcopy
      integer k

      i4_gcdb = 1
c
c  If both I and J are zero, I4_GCDB is 1.
c
      if ( i .eq. 0 .and. j .eq. 0 ) then
        i4_gcdb = 1
        return
      end if
c
c  If just one of I and J is zero, I4_GCDB is the other one.
c
      if ( i .eq. 0 ) then
        i4_gcdb = j
        return
      else if ( j .eq. 0 ) then
        i4_gcdb = i
        return
      end if
c
c  Divide out K as long as you can.
c
      if ( 0 .lt. j ) then
        i4_gcdb = 1
      else
        i4_gcdb = -1
      end if

      icopy = i
      jcopy = j

10    continue

        if ( mod ( icopy, k ) .ne. 0 .or.
     &       mod ( jcopy, k ) .ne. 0 ) then
          go to 20
        end if

        i4_gcdb = i4_gcdb * k
        icopy = icopy / k
        jcopy = jcopy / k

      go to 10

20    continue

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function i4_huge_normalizer ( )

c*********************************************************************72
c
cc I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
c
c  Discussion:
c
c    The value returned is 1 / ( I4_HUGE + 1 ).
c
c    For any I4, it should be the case that
c
c     -1 < I4 * I4_HUGE_NORMALIZER < 1.
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
c    Output, double precision I4_HUGE_NORMALIZER, the "normalizer"
c    for I4_HUGE.
c
      implicit none

      double precision i4_huge_normalizer

      i4_huge_normalizer = 4.656612873077392578125D-10

      return
      end
      function i4_is_power_of_2 ( n )

c*********************************************************************72
c
cc I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
c
c  Discussion:
c
c    The powers of 2 are 1, 2, 4, 8, 16, and so on.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be tested.
c
c    Output, logical I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
c
      implicit none

      logical i4_is_power_of_2
      integer n
      integer n_copy

      n_copy = n
      i4_is_power_of_2 = .false.

      if ( n_copy .le. 0 ) then
        return
      end if

10    continue

      if ( n_copy .ne. 1 ) then

        if ( mod ( n_copy, 2 ) .eq. 1 ) then
          return
        end if

        n_copy = n_copy / 2

        go to 10

      end if

      i4_is_power_of_2 = .true.

      return
      end
      function i4_is_prime ( n )

c*********************************************************************72
c
cc I4_IS_PRIME reports whether an I4 is prime.
c
c  Discussion:
c
c    A simple, unoptimized sieve of Erasthosthenes is used to
c    check whether N can be divided by any integer between 2
c    and SQRT(N).
c
c    Note that negative numbers, 0 and 1 are not considered prime.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be tested.
c
c    Output, logical I4_IS_PRIME, is TRUE if N is prime, and FALSE
c    otherwise.
c
      implicit none

      integer i
      logical i4_is_prime
      integer n
      integer nhi

      if ( n .le. 0 ) then
        i4_is_prime = .false.
        return
      end if

      if ( n .eq. 1 ) then
        i4_is_prime = .false.
        return
      end if

      if ( n .le. 3 ) then
        i4_is_prime = .true.
        return
      end if

      nhi = int ( sqrt ( real ( n ) ) )

      do i = 2, nhi
        if ( mod ( n, i ) .eq. 0 ) then
          i4_is_prime = .false.
          return
        end if
      end do

      i4_is_prime = .true.

      return
      end
      function i4_lcm ( i, j )

c*********************************************************************72
c
cc I4_LCM computes the least common multiple of two I4's.
c
c  Discussion:
c
c    The least common multiple may be defined as
c
c      LCM(I,J) = ABS( I * J ) / GCD(I,J)
c
c    where GCD(I,J) is the greatest common divisor of I and J.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the integers whose I4_LCM is desired.
c
c    Output, integer I4_LCM, the least common multiple of I and J.
c    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
c
      implicit none

      integer i
      integer i4_gcd
      integer j
      integer i4_lcm

      i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

      return
      end
      function i4_log_10 ( i )

c*********************************************************************72
c
cc I4_LOG_10 returns the integer part of the logarithm base 10 of ABS(X).
c
c  Discussion:
c
c    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
c
c  Example:
c
c        I  I4_LOG_10
c    -----  --------
c        0    0
c        1    0
c        2    0
c        9    0
c       10    1
c       11    1
c       99    1
c      100    2
c      101    2
c      999    2
c     1000    3
c     1001    3
c     9999    3
c    10000    4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 10 is desired.
c
c    Output, integer I4_LOG_10, the integer part of the logarithm base 10 of
c    the absolute value of X.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_10
      integer ten_pow

      if ( i .eq. 0 ) then

        i4_log_10 = 0

      else

        i4_log_10 = 0
        ten_pow = 10

        i_abs = abs ( i )

10      continue

        if ( ten_pow .le. i_abs ) then
          i4_log_10 = i4_log_10 + 1
          ten_pow = ten_pow * 10
          go to 10
        end if

      end if

      return
      end
      function i4_log_2 ( i )

c*********************************************************************72
c
cc I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
c
c  Discussion:
c
c    For positive I4_LOG_2(I), it should be true that
c      2^I4_LOG_2(X) .le. |I| < 2^(I4_LOG_2(I)+1).
c    The special case of I4_LOG_2(0) returns -HUGE().
c
c    An I4 is an integer value.
c
c  Example:
c
c     I  I4_LOG_2
c
c     0  -1
c     1,  0
c     2,  1
c     3,  1
c     4,  2
c     5,  2
c     6,  2
c     7,  2
c     8,  3
c     9,  3
c    10,  3
c   127,  6
c   128,  7
c   129,  7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 2
c    is desired.
c
c    Output, integer I4_LOG_2, the integer part of the
c    logarithm base 2 of the absolute value of I.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_2
      integer i4_huge
      parameter ( i4_huge = 2147483647 )

      if ( i .eq. 0 ) then

        i4_log_2 = - i4_huge

      else

        i4_log_2 = 0

        i_abs = abs ( i )

10      continue

        if ( 2 .le. i_abs ) then
          i_abs = i_abs / 2
          i4_log_2 = i4_log_2 + 1
          go to 10
        end if

      end if

      return
      end
      function i4_log_i4 ( i4, j4 )

c*********************************************************************72
c
cc I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
c
c  Discussion:
c
c    Only the integer part of the logarithm is returned.
c
c    If
c
c      K4 = I4_LOG_J4 ( I4, J4 ),
c
c    then we ordinarily have
c
c      J4^(K4-1) < I4 .le. J4^K4.
c
c    The base J4 should be positive, and at least 2.  If J4 is negative,
c    a computation is made using the absolute value of J4.  If J4 is
c    -1, 0, or 1, the logarithm is returned as 0.
c
c    The number I4 should be positive and at least 2.  If I4 is negative,
c    a computation is made using the absolute value of I4.  If I4 is
c    -1, 0, or 1, then the logarithm is returned as 0.
c
c    An I4 is an integer value.
c
c  Example:
c
c    I4  J4  K4
c
c     0   3   0
c     1   3   0
c     2   3   0
c     3   3   1
c     4   3   1
c     8   3   1
c     9   3   2
c    10   3   2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, the number whose logarithm is desired.
c
c    Input, integer J4, the base of the logarithms.
c
c    Output, integer I4_LOG_I4, the integer part of the logarithm
c    base abs(J4) of abs(I4).
c
      implicit none

      integer i4
      integer i4_abs
      integer i4_log_i4
      integer j4
      integer j4_abs
      integer value

      value = 0

      i4_abs = abs ( i4 )

      if ( 2 .le. i4_abs ) then

        j4_abs = abs ( j4 )

        if ( 2 .le. j4_abs ) then

10        continue

          if ( j4_abs .le. i4_abs ) then
            i4_abs = i4_abs / j4_abs
            value = value + 1
            go to 10
          end if

        end if

      end if

      i4_log_i4 = value

      return
      end
      function i4_log_r8 ( x, b )

c*********************************************************************72
c
cc I4_LOG_R8 returns the logarithm of an I4 to an R8 base.
c
c  Discussion:
c
c    The base B should be positive, but in any case only the absolute
c    value of B is considered.
c
c    The number X whose logarithm is desired should be positive, but
c    in any case only the absolute value of X is considered.
c
c    An I4 is an integer value.
c
c    An R8 is a double precision value.
c
c  Example:
c
c    If B is greater than 1, and X is positive:
c
c    if 1/B^2  <  X .le. 1/B   I4_LOG_R8(X) = -1,
c    if 1/B    <  X .le. 1     I4_LOG_R8(X) = 0,
c    if 1      .le. X <  B,    I4_LOG_R8(X) = 0,
c    if B      .le. X <  B^2   I4_LOG_R8(X) = 1,
c    if B^2    .le. X <  B^3   I4_LOG_R8(X) = 2.
c
c    For positive I4_LOG_R8(X), it should be true that
c
c      ABS(B)^I4_LOG_R8(X) .le. ABS(X) < ABS(B)^(I4_LOG_R8(X)+1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the number whose logarithm base B is
c    desired.  If X is 0, then I4_LOG_B is returned as -I4_HUGE().
c
c    Input, double precision B, the absolute value of the base of the
c    logarithms.  B must not be -1, 0, or 1.
c
c    Output, integer I4_LOG_R8, the integer part of the logarithm
c    base abs(B) of abs(X).
c
      implicit none

      double precision b
      double precision b_abs
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_log_r8
      integer value_sign
      integer x
      double precision x_abs

      if ( x .eq. 0 ) then
        i4_log_r8 = - i4_huge
        return
      end if

      b_abs = abs ( b )
      i4_log_r8 = 0

      if ( b_abs .eq. 1.0D+00 ) then
        return
      end if

      if ( b .eq. 0.0D+00 ) then
        return
      end if

      x_abs = abs ( dble ( x ) )

      if ( b_abs .lt. 1.0D+00 ) then
        value_sign = -1
        b_abs = 1.0D+00 / b_abs
      else
        value_sign = +1
      end if

      if ( 1.0D+00 .le. x_abs .and. x_abs .lt. b_abs ) then
        i4_log_r8 = value_sign * i4_log_r8
        return
      end if

10    continue

      if ( b_abs .lt. x_abs ) then
        x_abs = x_abs / b_abs
        i4_log_r8 = i4_log_r8 + 1
        go to 10
      end if

20    continue

      if ( x_abs * b_abs .le. 1.0D+00 ) then
        x_abs = x_abs * b_abs
        i4_log_r8 = i4_log_r8 - 1
        go to 20
      end if
c
c  If the absolute value of the base was less than 1, we inverted
c  earlier.  Now negate the logarithm to account for that.
c
      i4_log_r8 = value_sign * i4_log_r8

      return
      end
      subroutine i4_mant ( x, s, j, k, l )

c*********************************************************************72
c
cc I4_MANT computes the "mantissa" of a double precision number.
c
c  Discussion:
c
c    I4_MANT computes the "mantissa" or "fraction part" of a real
c    number X, which it stores as a pair of integers, (J/K).
c
c    It also computes the sign, and the integer part of the logarithm
c    (base 2) of X.
c
c    On return:
c
c      X = S * (J/K) * 2^L
c
c    where
c
c      S is +1 or -1,
c      K is a power of 2,
c      1 .le. (J/K) .lt. 2,
c      L is an integer.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number to be decomposed.
c
c    Output, integer S, the "sign" of the number.
c    S will be -1 if X is less than 0, and +1 if X is greater
c    than or equal to zero.
c
c    Output, integer J, the top part of the mantissa fraction.
c
c    Output, integer K, the bottom part of the mantissa
c    fraction.  K is a power of 2.
c
c    Output, integer L, the integer part of the logarithm
c    (base 2) of X.
c
      implicit none

      integer j
      integer k
      integer l
      integer s
      double precision x
      double precision xtemp
c
c  1: Handle the special case of 0.
c
      if ( x .eq. 0.0D+00 ) then
        s = 1
        j = 0
        k = 1
        l = 0
        return
      end if
c
c  2: Determine the sign S.
c
      if ( 0.0D+00 .lt. x ) then
        s = 1
        xtemp = x
      else
        s = -1
        xtemp = -x
      end if
c
c  3: Force XTEMP to lie between 1 and 2, and compute the logarithm L.
c
      l = 0

10    continue

      if ( 2.0D+00 .le. xtemp ) then
        xtemp = xtemp / 2.0D+00
        l = l + 1
        go to 10
      end if

20    continue

      if ( xtemp .lt. 1.0D+00 ) then
        xtemp = xtemp * 2.0D+00
        l = l - 1
        go to 20
      end if
c
c  4: Now strip out the mantissa as J/K.
c
      j = 0
      k = 1

30    continue

        j = 2 * j

        if ( 1.0D+00 .le. xtemp ) then
          j = j + 1
          xtemp = xtemp - 1.0D+00
        end if

        if ( xtemp .eq. 0.0D+00 ) then
          go to 40
        end if

        k = 2 * k
        xtemp = xtemp * 2.0D+00

      go to 30

40    continue

      return
      end
      subroutine i4_mod_inv ( b, n, y )

c*********************************************************************72
c
cc I4_MOD_INV calculates the inverse of B mod N.
c
c  Discussion:
c
c    This function uses the extended Euclidean algorithm.
c
c    Unless the algorithm fails, the output value Y will satisfy
c
c      ( B * Y ) mod N = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2011
c
c  Author:
c
c    Original MATLAB version by Wade Trappe, Lawrence Washington.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wade Trappe, Lawrence Washington,
c    Introduction to Cryptography with Coding Theory,
c    Prentice Hall, 2005,
c    ISBN13: 978-0131862395,
c    LC: QA268.T73.
c
c  Parameters:
c
c    Input, integer B, the value whose inverse is desired.
c    B must not be 0, or a multiple of N.  However, B can be negative.
c
c    Input, integer N, the value with respect to which the inverse
c    is desired.  N must be 2 or greater.
c
c    Output, integer Y, the inverse of B mod N.  However, if the
c    inverse does not exist, Y is returned as 0.
c
      implicit none

      integer b
      integer b0
      integer n
      integer n0
      integer q
      integer r
      integer t
      integer t0
      integer temp
      integer y

      n0 = n
      b0 = abs ( b )
      t0 = 0
      t = 1

      q = ( n0 / b0 )
      r = n0 - q * b0

10    continue

      if ( 0 .lt. r ) then

        temp = t0 - q * t

        if ( 0 .le. temp ) then
          temp =     mod (   temp, n )
        else
          temp = n - mod ( - temp, n )
        end if

        n0 = b0
        b0 = r
        t0 = t
        t = temp

        q = ( n0 / b0 )
        r = n0 - q * b0

        go to 10

      end if

      if ( b0 .ne. 1 ) then
        y = 0
      else
        y = mod ( t, n )
        if ( b .lt. 0 ) then
          y = - y
        end if
      end if

      return
      end
      subroutine i4_moddiv ( n, d, m, r )

c*********************************************************************72
c
cc I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
c
c  Discussion:
c
c    The formula used is:
c
c      N = M * D + R
c
c      0 <= || R || < || D ||
c
c    and R has the sign of N.
c
c    An I4 is an integer value.
c
c  Example:
c
c      N         D       M      R
c
c     107       50      2      7
c     107      -50     -2      7
c    -107       50     -2     -7
c    -107      -50      2     -7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number to be decomposed.
c
c    Input, integer D, the divisor.  D may not be zero.
c
c    Output, integer M, the number of times N
c    is evenly divided by D.
c
c    Output, integer R, a remainder.
c
      implicit none

      integer d
      integer m
      integer n
      integer r

      if ( d .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODDIV - Fatal error!'
        write ( *, '(a)' ) '  Input divisor D = 0'
        stop
      end if

      m = n / d
      r = n - d * m

      return
      end
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

      return
      end
      function i4_mop ( i )

c*********************************************************************72
c
cc I4_MOP returns the I-th power of -1 as an I4 value.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, integer I4_MOP, the I-th power of -1.
c
      implicit none

      integer i
      integer i4_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        i4_mop = 1
      else
        i4_mop = -1
      end if

      return
      end
      function i4_odd ( i )

c*********************************************************************72
c
cc I4_ODD returns TRUE if an I4 is odd.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be tested.
c
c    Output, logical I4_ODD, is TRUE if I is odd.
c
      implicit none

      integer i
      logical i4_odd

      i4_odd = ( mod ( i + 1, 2 ) .eq. 0 )

      return
      end
      function i4_power ( i, j )

c*********************************************************************72
c
cc I4_POWER returns the integer power of an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the base and the power.
c    J should be nonnegative.
c
c    Output, integer I4_POWER, the value of I^J.
c
      implicit none

      integer i
      integer i4_power
      integer j
      integer k

      if ( j .lt. 0 ) then

        if ( i .eq. 1 ) then
          i4_power = 1
        else if ( i .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_POWER - Fatal error!'
          write ( *, '(a)' ) '  I^J requested, with I = 0, J negative.'
          stop
        else
          i4_power = 0
        end if

      else if ( j .eq. 0 ) then

        if ( i .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_POWER - Fatal error!'
          write ( *, '(a)' ) '  I^J requested, with I = 0 and J = 0.'
          stop
        else
          i4_power = 1
        end if

      else if ( j .eq. 1 ) then

        i4_power = i

      else

        i4_power = 1
        do k = 1, j
          i4_power = i4_power * i
        end do

      end if

      return
      end
      function i4_sign ( x )

c*********************************************************************72
c
cc I4_SIGN evaluates the sign of an I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the number whose sign is desired.
c
c    Output, integer I4_SIGN, the sign of the number.
c
      implicit none

      integer i4_sign
      integer x

      if ( x .lt. 0 ) then
        i4_sign = -1
      else
        i4_sign = +1
      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      subroutine i4_swap3 ( i, j, k )

c*********************************************************************72
c
cc I4_SWAP3 swaps three I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J, K.  On output, the values of I, J, and K
c    have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k
      integer l

      l = i
      i = j
      j = k
      k = l

      return
      end
      subroutine i4_to_angle ( i, angle )

c*********************************************************************72
c
cc I4_TO_ANGLE maps I4's to points on a circle.
c
c  Discussion:
c
c    The angles are intended to be used to select colors on a color
c    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
c    magenta.
c
c    An I4 is an integer value.
c
c  Example:
c
c     I   X      ANGLE
c
c     0   0/3      0
c     1   1/3    120
c     2   2/3    240
c
c     3   1/6     60
c     4   3/6    180
c     5   5/6    300
c
c     6   1/12    30
c     7   3/12    90
c     8   5/12   150
c     9   7/12   210
c    10   9/12   270
c    11  11/12   330
c
c    12   1/24    15
c    13   3/24    45
c    14   5/24    75
c    etc
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the desired color.
c
c    Output, double precision ANGLE, an angle, measured in degrees,
c    between 0 and 360.
c
      implicit none

      double precision angle
      integer i
      integer i4_log_2
      integer i1
      integer i2
      integer i3
      integer i4

      if ( 0 .le. abs ( i ) .and. abs ( i ) .le. 2 ) then

        angle = 120.0D+00 * dble ( abs ( i ) )

      else

        i1 = i4_log_2 ( abs ( i ) / 3 )
        i2 = abs ( i ) + 1 - 3 * 2**i1
        i3 = 2 * ( i2 - 1 ) + 1
        i4 = 3 * 2**( i1 + 1 )

        angle = 360.0D+00 * dble ( i3 ) / dble ( i4 )

      end if

      return
      end
      subroutine i4_to_digits_binary ( i, n, c )

c*********************************************************************72
c
cc I4_TO_DIGITS_BINARY produces the binary digits of an I4.
c
c  Discussion:
c
c    An I4 is an integer.
c
c  Example:
c
c     I    N     C               Binary
c    --  ---   ---         ------------
c     0    1   0                      0
c     0    2   0, 0                  00
c     1    3   1, 0, 0              100
c     2    3   0, 1, 0              010
c     3    3   1, 1, 0              011
c     4    3   0, 0, 1              100
c     8    3   0, 0, 0           (1)000
c     8    5   0, 0, 0, 1, 0      01000
c    -8    5   0, 0, 0, 1, 0  (-) 01000
c
c     0    3   0, 0, 0
c     1    3   1, 0, 0
c     2    3   0, 1, 0
c     3    3   1, 1, 0
c     4    3   0, 0, 1
c     5    3   1, 0, 1
c     6    3   0, 1, 1
c     7    3   1, 1, 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, an integer to be represented.
c
c    Input, integer N, the number of binary digits to produce.
c
c    Output, integer C(N), the first N binary digits of I,
c    with C(1) being the units digit.
c
      implicit none

      integer n

      integer c(n)
      integer i
      integer i_copy
      integer j

      i_copy = abs ( i )

      do j = 1, n

        c(j) = mod ( i_copy, 2 )
        i_copy = i_copy / 2

      end do

      return
      end
      subroutine i4_to_digits_decimal ( i, n, digit )

c*********************************************************************72
c
cc I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be analyzed.
c
c    Input, integer N, the number of digits to determine.
c
c    Output, integer DIGIT(N), the last N decimal digits of I.
c    DIGIT(I) is the "coefficient" of 10**(I-1).
c
      implicit none

      integer n

      integer digit(n)
      integer i
      integer i_copy
      integer i4_ten
      parameter ( i4_ten = 10 )
      integer j

      i_copy = i

      do j = 1, n
        digit(j) = mod ( i_copy, i4_ten )
        i_copy = ( i_copy - digit(j) ) / 10
      end do

      return
      end
      subroutine i4_to_fac ( intval, prime_num, npower )

c*********************************************************************72
c
cc I4_TO_FAC converts an I4 into a product of prime factors.
c
c  Discussion:
c
c    This routine will fail if the input integer is not positive,
c    or if PRIME_NUM is too small to account for the factors of the integer.
c
c    An I4 is an integer value.
c
c    The formula is:
c
c      INTVAL = Product ( 1 <= I <= PRIME_NUM ) PRIME(I)**NPOWER(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, the integer to be factored.
c
c    Input, integer PRIME_NUM, the number of prime factors for
c    which storage has been allocated.
c
c    Output, integer NPOWER(PRIME_NUM), the powers of the primes.
c
      implicit none

      integer prime_num

      integer i
      integer intcopy
      integer intval
      integer npower(prime_num)
      integer p
      integer prime

      if ( intval .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_FAC - Fatal error!'
        write ( *, '(a)' ) '  Input integer is not positive.'
        stop
      end if
c
c  Try dividing the remainder by each prime.
c
      intcopy = intval

      do i = 1, prime_num

        npower(i) = 0

        p = prime ( i )

10      continue

        if ( mod ( intcopy, p ) .eq. 0 ) then
          npower(i) = npower(i) + 1
          intcopy = intcopy / p
          go to 10
        end if

      end do

      return
      end
      subroutine i4_to_halton ( dim_num, step, seed, leap, base, r )

c*********************************************************************72
c
cc I4_TO_HALTON computes one element of a leaped Halton subsequence.
c
c  Discussion:
c
c    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
c    sequences, each generated by a particular base.
c
c    This routine selects elements of a "leaped" subsequence of the
c    Halton sequence.  The subsequence elements are indexed by a
c    quantity called STEP, which starts at 0.  The STEP-th subsequence
c    element is simply element
c
c      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
c
c    of the original Halton sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Halton,
c    On the efficiency of certain quasi-random sequences of points
c    in evaluating multi-dimensional integrals,
c    Numerische Mathematik,
c    Volume 2, 1960, pages 84-90.
c
c    John Halton, GB Smith,
c    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
c    Communications of the ACM,
c    Volume 7, 1964, pages 701-702.
c
c    Ladislav Kocis, William Whiten,
c    Computational Investigations of Low-Discrepancy Sequences,
c    ACM Transactions on Mathematical Software,
c    Volume 23, Number 2, 1997, pages 266-294.
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c    1 <= DIM_NUM is required.
c
c    Input, integer STEP, the index of the subsequence element.
c    0 <= STEP is required.
c
c    Input, integer SEED(DIM_NUM), the Halton sequence index corresponding
c    to STEP = 0.
c    0 <= SEED(1:DIM_NUM) is required.
c
c    Input, integer LEAP(DIM_NUM), the successive jumps in the Halton sequence.
c    1 <= LEAP(1:DIM_NUM) is required.
c
c    Input, integer BASE(DIM_NUM), the Halton bases.
c    1 < BASE(1:DIM_NUM) is required.
c
c    Output, double precision R(DIM_NUM), the STEP-th element of the leaped
c    Halton subsequence.
c
      implicit none

      integer dim_num

      integer base(dim_num)
      double precision base_inv
      integer digit
      integer i
      integer leap(dim_num)
      double precision r(dim_num)
      integer seed(dim_num)
      integer seed2
      integer step
c
c  Check the input.
c
      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM < 1.'
        stop
      end if

      if ( step .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
        write ( *, '(a)' ) ' STEP < 0.'
        stop
      end if

      do i = 1, dim_num

        if ( seed(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
          write ( *, '(a)' ) '  Some SEED(*) < 0.'
          stop
        end if

      end do

      do i = 1, dim_num

        if ( leap(i) .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
          write ( *, '(a)' ) '  Some LEAP < 1.'
          stop
        end if

      end do

      do i = 1, dim_num

        if ( base(i) .le. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
          write ( *, '(a)' ) '  Some BASE <= 1.'
          stop
        end if

      end do
c
c  Calculate the data.
c
      do i = 1, dim_num

        seed2 = seed(i) + step * leap(i)

        r(i) = 0.0D+00

        base_inv = 1.0D+00 / dble ( base(i) )

10      continue

        if ( seed2 .ne. 0 ) then
          digit = mod ( seed2, base(i) )
          r(i) = r(i) + dble ( digit ) * base_inv
          base_inv = base_inv / dble ( base(i) )
          seed2 = seed2 / base(i)
          go to 10
        end if

      end do

      return
      end
      function i4_to_isbn ( i )

c*********************************************************************72
c
cc I4_TO_ISBN converts an I4 to an ISBN digit.
c
c  Discussion:
c
c    Only the integers 0 through 10 can be input.  The representation
c    of 10 is 'X'.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, integer I, an integer between 0 and 10.
c
c    Output, character I4_TO_ISBN, the ISBN character code of the integer.
c    If I is illegal, then I4_TO_ISBN is set to '?'.
c
      implicit none

      integer i
      character i4_to_isbn

           if ( i .eq. 0 ) then
        i4_to_isbn = '0'
      else if ( i .eq. 1 ) then
        i4_to_isbn = '1'
      else if ( i .eq. 2 ) then
        i4_to_isbn = '2'
      else if ( i .eq. 3 ) then
        i4_to_isbn = '3'
      else if ( i .eq. 4 ) then
        i4_to_isbn = '4'
      else if ( i .eq. 5 ) then
        i4_to_isbn = '5'
      else if ( i .eq. 6 ) then
        i4_to_isbn = '6'
      else if ( i .eq. 7 ) then
        i4_to_isbn = '7'
      else if ( i .eq. 8 ) then
        i4_to_isbn = '8'
      else if ( i .eq. 9 ) then
        i4_to_isbn = '9'
      else if ( i .eq. 10 ) then
        i4_to_isbn = 'X'
      else
        i4_to_isbn = '?'
      end if

      return
      end
      function i4_to_l ( i4 )

c*********************************************************************72
c
cc I4_TO_L converts an I4 to a logical value.
c
c  Discussion:
c
c    0 is FALSE, and anything else if TRUE.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, an integer.
c
c    Output, logical I4_TO_L, the logical value of I4.
c
      implicit none

      integer i4
      logical i4_to_l
      logical value

      value = ( i4 .ne. 0 )

      i4_to_l = value

      return
      end
      function i4_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM_AB, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_uniform_ab
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
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

      i4_uniform_ab = value

      return
      end
      subroutine i4_unswap3 ( i, j, k )

c*********************************************************************72
c
cc I4_UNSWAP3 unswaps three I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J, K.  On output, the values of I, J, and K
c    have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k
      integer l

      l = k
      k = j
      j = i
      i = l

      return
      end
      function i4_walsh_1d ( x, digit )

c*********************************************************************72
c
cc I4_WALSH_1D evaluates the Walsh function.
c
c  Discussion:
c
c    Consider the binary representation of X, and number the digits
c    in descending order, from leading to lowest, with the units digit
c    being numbered 0.
c
c    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the Walsh function.
c
c    Input, integer DIGIT, the index of the Walsh function.
c
c    Output, integer I4_WALSH_1D, the value of the Walsh function.
c
      implicit none

      integer digit
      integer i4_walsh_1d
      integer n
      double precision x
      double precision x_copy
c
c  Hide the effect of the sign of X.
c
      x_copy = abs ( x )
c
c  If DIGIT is positive, divide by 2 DIGIT times.
c  If DIGIT is negative, multiply by 2 (-DIGIT) times.
c
      x_copy = x_copy / 2.0D+00**digit
c
c  Make it an integer.
c  Because it's positive, and we're using INT, we don't change the
c  units digit.
c
      n = int ( x_copy )
c
c  Is the units digit odd or even?
c
      if ( mod ( n, 2 ) .eq. 0 ) then
        i4_walsh_1d = 0
      else
        i4_walsh_1d = 1
      end if

      return
      end
      function i4_width ( i )

c*********************************************************************72
c
cc I4_WIDTH returns the "width" of an I4.
c
c  Discussion:
c
c    The width of an integer is the number of characters necessary to print it.
c
c    The width of an integer can be useful when setting the appropriate output
c    format for a vector or array of values.
c
c    An I4 is an integer value.
c
c  Example:
c
c        I  I4_WIDTH
c    -----  -------
c    -1234    5
c     -123    4
c      -12    3
c       -1    2
c        0    1
c        1    1
c       12    2
c      123    3
c     1234    4
c    12345    5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose width is desired.
c
c    Output, integer I4_WIDTH, the number of characters
c    necessary to represent the integer in base 10, including a negative
c    sign if necessary.
c
      implicit none

      integer i
      integer i4_log_10
      integer i4_width

      if ( 0 .le. i ) then
        i4_width = i4_log_10 ( i ) + 1
      else
        i4_width = i4_log_10 ( i ) + 2
      end if

      return
      end
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      function i4_xor ( i, j )

c*********************************************************************72
c
cc I4_XOR calculates the exclusive OR of two I4's.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The FORTRAN intrinsinc IEOR should be used instead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c   John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two values whose exclusive OR is needed.
c
c    Output, integer I4_XOR, the exclusive OR of I and J.
c
      implicit none

      integer i
      integer i1
      integer i2
      integer i4_xor
      integer j
      integer j1
      integer j2
      integer k
      integer l

      i1 = i
      j1 = j
      k = 0
      l = 1

10    continue

      if ( i1 .ne. 0 .or. j1 .ne. 0 ) then

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

      end if

      i4_xor = k

      return
      end
      subroutine i43mat_flip_cols ( m, n, a )

c*********************************************************************72
c
cc I43MAT_FLIP_COLS swaps the columns of an I43MAT.
c
c  Discussion:
c
c    An I43MAT is a matrix, each of whose entries is an I43,
c    a triple of I4's.
c
c    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
c    and N counts the "rows".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(3,M,N), the matrix whose columns
c    are to be flipped.
c
      implicit none

      integer m
      integer n

      integer a(3,m,n)
      integer b
      integer i
      integer j
      integer k

      do k = 1, n / 2
        do j = 1, m
          do i = 1, 3
            b            = a(i,j,    k)
            a(i,j,    k) = a(i,j,n+1-k)
            a(i,j,n+1-k) = b
          end do
        end do
      end do

      return
      end
      subroutine i43mat_flip_rows ( m, n, a )

c*********************************************************************72
c
cc I43MAT_FLIP_ROWS swaps the rows of an I43MAT.
c
c  Discussion:
c
c    An I43MAT is a matrix, each of whose entries is an I43,
c    a triple of I4's.
c
c    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
c    and N counts the "rows".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(3,M,N), the matrix whose rows
c    are to be flipped.
c
      implicit none

      integer m
      integer n

      integer a(3,m,n)
      integer b
      integer i
      integer j
      integer k

      do k = 1, n
        do j = 1, m / 2
          do i = 1, 3
            b            = a(i,    j,k)
            a(i,    j,k) = a(i,m+1-j,k)
            a(i,m+1-j,k) = b
          end do
        end do
      end do

      return
      end
      subroutine i4block_print ( l, m, n, a, title )

c*********************************************************************72
c
cc I4BLOCK_PRINT prints an I4BLOCK.
c
c  Discussion:
c
c    An I4BLOCK is a 3D array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L, M, N, the dimensions of the block.
c
c    Input, integer A(L,M,N), the matrix to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer l
      integer m
      integer n

      integer a(l,m,n)
      integer i
      integer j
      integer jhi
      integer jlo
      integer k
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do k = 1, n

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  K = ', k

        do jlo = 1, m, 10
          jhi = min ( jlo + 10 - 1, m )
          write ( *, '(a)' ) ' '
          write ( *, '(8x,a2,10(2x,i6))' ) 'J:', ( j, j = jlo, jhi )
          write ( *, '(7x,a2)' ) 'I:'
          do i = 1, l
            write ( *, '(2x,i6,a1,1x,10(2x,i6))' ) 
     &        i, ':', a(i,jlo:jhi,k)
          end do
        end do

      end do

      return
      end
      subroutine i4col_compare ( m, n, a, i, j, isgn )

c*********************************************************************72
c
cc I4COL_COMPARE compares columns I and J of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, I = 2, J = 4
c
c      A = (
c        1  2  3  4
c        5  6  7  8
c        9 10 11 12 )
c
c    Output:
c
c      ISGN = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of N columns of
c    vectors of length M.
c
c    Input, integer I, J, the columns to be compared.
c    I and J must be between 1 and N.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, column I < column J,
c     0, column I = column J,
c    +1, column J < column I.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer isgn
      integer j
      integer k
c
c  Check.
c
      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' )
     &    '  Column index I = ', i, ' is less than 1.'
        stop
      end if

      if ( n .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' )
     &    '  N = ', n, ' is less than column index I = ', i
        stop
      end if

      if ( j .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' )
     &    '  Column index J = ', j, ' is less than 1.'
        stop
      end if

      if ( n .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' )
     &    '  N = ', n, ' is less than column index J = ', j
        stop
      end if

      isgn = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. m ) then

        if ( a(k,i) .lt. a(k,j) ) then
          isgn = -1
          return
        else if ( a(k,j) .lt. a(k,i) ) then
          isgn = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4col_find ( m, n, a, ivec, col )

c*********************************************************************72
c
cc I4COL_FIND searches an I4COL for a particular column value.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Example:
c
c    M = 3, N = 4,
c
c    A = (
c      1  2  3  4
c      5  6  7  8
c      9 10 11 12 )
c
c    IVEC = ( 3, 7, 11 )
c
c    COL = 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in
c    the table.  M is also the length of IVEC.
c
c    Input, integer A(M,N), an array of N columns of vectors
c    of length M.
c
c    Input, integer IVEC(M), a vector to be matched with the data
c    in the array.
c
c    Output, integer COL, the index of the first column of
c    the table which exactly matches every entry of IVEC, or -1 if no match
c    could be found.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer col
      integer ivec(m)
      integer j

      if ( m .le. 0 ) then
        col = -1
        return
      end if

      do j = 1, n

        i = 1

10      continue

        if ( ivec(i) .eq. a(i,j) ) then

          if ( i .eq. m ) then
            col = j
            return
          end if

          i = i + 1
          go to 10

        end if

      end do

      col = -1

      return
      end
      subroutine i4col_find_item ( m, n, a, item, row, col )

c*********************************************************************72
c
cc I4COL_FIND_ITEM searches an I4COL for a given scalar value.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in
c    the table.
c
c    Input, integer A(M,N), an array of N columns of vectors
c    of length M.
c
c    Input, integer ITEM, the value to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM.  The search
c    is conducted by columns.  If the item is not found, then
c    ROW = COL = -1.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer item
      integer j
      integer row

      do j = 1, n
        do i = 1, m
          if ( a(i,j) .eq. item ) then
            row = i
            col = j
            return
          end if
        end do
      end do

      row = -1
      col = -1

      return
      end
      subroutine i4col_find_pair_wrap ( m, n, a, item1, item2, row,
     &  col )

c*********************************************************************72
c
cc I4COL_FIND_PAIR_WRAP searches an I4COL for a pair of items.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    The items (ITEM1, ITEM2) must occur consecutively.
c    However, wrapping is allowed, that is, if ITEM1 occurs
c    in the last row, and ITEM2 "follows" it in the first row
c    of the same column, a match is declared.
c
c    If the pair of items is not found, then ROW = COL = -1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the array.
c
c    Input, integer A(M,N), the array to search.
c
c    Input, integer ITEM1, ITEM2, the values to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM1 followed immediately
c    by ITEM2.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer i2
      integer item1
      integer item2
      integer j
      integer row

      do j = 1, n
        do i = 1, m

          if ( a(i,j) .eq. item1 ) then

            i2 = i + 1

            if ( m .lt. i2 ) then
              i2 = 1
            end if

            if ( a(i2,j) .eq. item2 ) then
              row = i
              col = j
              return
            end if

          end if

        end do
      end do

      row = -1
      col = -1

      return
      end
      subroutine i4col_first_index ( m, n, a, first_index )

c*********************************************************************72
c
cc I4COL_FIRST_INDEX indexes the first occurrence of values in an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values.
c    It is regarded as an array of N columns of length M.
c
c    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
c    the first column whose entries are equal to A(1:M,J).
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
c    Input, integer M, N, the number of rows and columns of A.
c    The length of an "element" of A, and the number of "elements".
c
c    Input, integer A(M,N), the array.
c
c    Output, integer FIRST_INDEX(N), the first occurrence index.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical all_equal
      integer first_index(n)
      integer i
      integer j1
      integer j2

      do i = 1, n
        first_index(i) = -1
      end do

      do j1 = 1, n

        if ( first_index(j1) .eq. -1 ) then

          first_index(j1) = j1

          do j2 = j1 + 1, n
            all_equal = .true.
            do i = 1, m
              if ( a(i,j1) .ne. a(i,j2) ) then
                all_equal = .false.
              end if
            end do
            if ( all_equal ) then
              first_index(j2) = j1
            end if
          end do

        end if

      end do

      return
      end
      subroutine i4col_sort_a ( m, n, a )

c*********************************************************************72
c
cc I4COL_SORT_A ascending sorts an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    In lexicographic order, the statement "X < Y", applied to two real
c    vectors X and Y of length M, means that there is some index I, with
c    1 <= I <= M, with the property that
c
c      X(J) = Y(J) for J < I,
c    and
c      X(I) < Y(I).
c
c    In other words, the first time they differ, X is smaller.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A, and the length of
c    a vector of data.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of N columns of M-vectors.
c    On output, the columns of A have been sorted in ascending
c    lexicographic order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer indx
      integer isgn
      integer j

      if ( m .le. 0 ) then
        return
      end if

      if ( n .le. 1 ) then
        return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4col_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4col_compare ( m, n, a, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4col_sort_d ( m, n, a )

c*********************************************************************72
c
cc I4COL_SORT_D descending sorts an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    In lexicographic order, the statement "X < Y", applied to two real
c    vectors X and Y of length M, means that there is some index I, with
c    1 <= I <= M, with the property that
c
c      X(J) = Y(J) for J < I,
c    and
c      X(I) < Y(I).
c
c    In other words, the first time they differ, X is smaller.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A, and the length of
c    a vector of data.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of N columns of M-vectors.
c    On output, the columns of A have been sorted in descending
c    lexicographic order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer indx
      integer isgn
      integer j

      if ( m .le. 0 ) then
        return
      end if

      if ( n .le. 1 ) then
        return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4col_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4col_compare ( m, n, a, i, j, isgn )
          isgn = -isgn

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4col_sort2_a ( m, n, a )

c*********************************************************************72
c
cc I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A, and the length
c    of a vector of data.
c
c    Input/output, integer A(M,N).
c    On input, the array of N columns of M vectors.
c    On output, the elements of each column of A have been sorted in ascending
c    order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer indx
      integer isgn
      integer j
      integer t

      if ( m .le. 1 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if
c
c  Initialize.
c
      do col = 1, n

        i = 0
        indx = 0
        isgn = 0
        j = 0
c
c  Call the external heap sorter.
c
10      continue

          call sort_heap_external ( m, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
          if ( 0 .lt. indx ) then

            t        = a(i,col)
            a(i,col) = a(j,col)
            a(j,col) = t
c
c  Compare the I and J objects.
c
          else if ( indx .lt. 0 ) then

            if ( a(j,col) .lt. a(i,col) ) then
              isgn = +1
            else
              isgn = -1
            end if

          else if ( indx .eq. 0 ) then

            go to 20

          end if

        go to 10

20    continue

      end do

      return
      end
      subroutine i4col_sort2_d ( m, n, a )

c*********************************************************************72
c
cc I4COL_SORT2_D descending sorts elements of each column of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A, and the length
c    of a vector of data.
c
c    Input/output, integer A(M,N).
c    On input, the array of N columns of M vectors.
c    On output, the elements of each column of A have been sorted in descending
c    order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer indx
      integer isgn
      integer  j
      integer t

      if ( m .le. 1 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if
c
c  Initialize.
c
      do col = 1, n

        i = 0
        indx = 0
        isgn = 0
        j = 0
c
c  Call the external heap sorter.
c
10      continue

          call sort_heap_external ( m, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
          if ( 0 .lt. indx ) then

            t        = a(i,col)
            a(i,col) = a(j,col)
            a(j,col) = t
c
c  Compare the I and J objects.
c
          else if ( indx .lt. 0 ) then

            if ( a(i,col) .lt. a(j,col) ) then
              isgn = +1
            else
              isgn = -1
            end if

          else if ( indx .eq. 0 ) then

            go to 20

          end if

        go to 10

20      continue

      end do

      return
      end
      subroutine i4col_sorted_singleton_count ( m, n, a, singleton_num )

c*********************************************************************72
c
cc I4COL_SORTED_SINGLETON_COUNT counts singletons in an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    The columns of the array may be ascending or descending sorted.
c
c    A "singleton" is an item that occurs exactly once.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), a sorted array, containing
c    N columns of data.
c
c    Output, integer SINGLETON_NUM, the number of singletons.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical differ_from_next
      logical differ_from_previous
      integer i
      integer j
      integer singleton_num

      singleton_num = 0

      if ( n .le. 0 ) then
        return
      end if

      differ_from_next = .true.

      do j = 1, n

        differ_from_previous = differ_from_next

        if ( j .lt. n ) then

          differ_from_next = .false.

          do i = 1, m
            if ( a(i,j) .ne. a(i,j+1) ) then
              differ_from_next = .true.
              go to 10
            end if
          end do

10        continue

        else

          differ_from_next = .true.

        end if

        if ( differ_from_previous .and. differ_from_next ) then
          singleton_num = singleton_num + 1
        end if

      end do

      return
      end
      subroutine i4col_sorted_unique ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4COL_SORTED_UNIQUE keeps unique elements in a sorted I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    The array can be sorted into ascending or descending order.
c    The important point is that identical elements must be stored
c    in adjacent positions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c   John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A, and the length of
c    a vector of data.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the sorted array of N columns of M-vectors.
c    On output, a sorted array of columns of M-vectors.
c
c    Output, integer UNIQUE_NUM, the number of unique columns of A.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i1
      integer i2
      integer j1
      integer j2
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      j1 = 1

      do j2 = 2, n

        do i = 1, m

          if ( a(i,j1) .ne. a(i,j2) ) then
            j1 = j1 + 1
            do i2 = 1, m
              a(i2,j1) = a(i2,j2)
            end do
            go to 10
          end if

        end do

10      continue

      end do

      unique_num = j1

      return
      end
      subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
c
c  Discussion:
c
c    The columns of the array may be ascending or descending sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), a sorted array, containing
c    N columns of data.
c
c    Output, integer UNIQUE_NUM, the number of unique columns.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j1
      integer j2
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1
      j1 = 1

      do j2 = 2, n

        do i = 1, m
          if ( a(i,j1) .ne. a(i,j2) ) then
            unique_num = unique_num + 1
            j1 = j2
            go to 10
          end if
        end do

10      continue

      end do

      return
      end
      subroutine i4col_swap ( m, n, a, j1, j2 )

c*********************************************************************72
c
cc I4COL_SWAP swaps columns J1 and J2 of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, J1 = 2, J2 = 4
c
c      A = (
c        1  2  3  4
c        5  6  7  8
c        9 10 11 12 )
c
c    Output:
c
c      A = (
c        1  4  3  2
c        5  8  7  6
c        9 12 11 10 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the array.
c
c    Input/output, integer A(M,N), an array of N columns
c    of length M.
c
c    Input, integer J1, J2, the columns to be swapped.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j1
      integer j2
      integer t

      if ( j1 .lt. 1 .or. n .lt. j1 .or.
     &     j2 .lt. 1 .or. n .lt. j2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
        write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
        write ( *, '(a,i8)' ) '  J1 =    ', j1
        write ( *, '(a,i8)' ) '  J2 =    ', j2
        write ( *, '(a,i8)' ) '  N =     ', n
        stop

      end if

      if ( j1 .eq. j2 ) then
        return
      end if

      do i = 1, m
        t       = a(i,j1)
        a(i,j1) = a(i,j2)
        a(i,j2) = t
      end do

      return
      end
      subroutine i4col_unique_index ( m, n, a, unique_index )

c*********************************************************************72
c
cc I4COL_UNIQUE_INDEX indexes the first occurrence of values in an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values.
c    It is regarded as an array of N columns of length M.
c
c    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
c    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
c    gathered in order, then
c
c      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
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
c    Input, integer M, N, the number of rows and columns of A.
c    The length of an "element" of A, and the number of "elements".
c
c    Input, integer A(M,N), the array.
c
c    Output, integer UNIQUE_INDEX(N), the unique index.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical all_equal
      integer i
      integer j1
      integer j2
      integer unique_index(n)
      integer unique_num

      do i = 1, n
        unique_index(i) = -1
      end do
      unique_num = 0

      do j1 = 1, n

        if ( unique_index(j1) .eq. -1 ) then

          unique_num = unique_num + 1
          unique_index(j1) = unique_num

          do j2 = j1 + 1, n
            all_equal = .true.
            do i = 1, m
              if ( a(i,j1) .ne. a(i,j2) ) then
                all_equal = .false.
              end if
            end do
            if ( all_equal ) then
              unique_index(j2) = unique_num
            end if
          end do

        end if

      end do

      return
      end
      subroutine i4i4_sort_a ( i1, i2, j1, j2 )

c*********************************************************************72
c
cc I4I4_SORT_A ascending sorts a pair of integers.
c
c  Discussion:
c
c    An I4I4 is a pair of integers, regarded as a single data item.
c
c    The program allows the reasonable call:
c
c      call i4i4_sort_a ( i1, i2, i1, i2 )
c
c    and this will return the reasonable result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I1, I2, the values to sort.
c
c    Output, integer J1, J2, the sorted values.
c
      implicit none

      integer i1
      integer i2
      integer j1
      integer j2
      integer k1
      integer k2
c
c  Copy arguments, so that the user can make "reasonable" calls like:
c
c    call i4i4_sort_a ( i1, i2, i1, i2 )
c
      k1 = i1
      k2 = i2

      j1 = min ( k1, k2 )
      j2 = max ( k1, k2 )

      return
      end
      subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

c*********************************************************************72
c
cc I4I4I4_SORT_A ascending sorts a triple of integers.
c
c  Discussion:
c
c    An I4I4I4 is a triple of integers, regarded as a single data item.
c
c    The program allows the reasonable call:
c
c      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
c
c    and this will return the reasonable result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I1, I2, I3, the values to sort.
c
c    Output, integer J1, J2, J3, the sorted values.
c
      implicit none

      integer i1
      integer i2
      integer i3
      integer j1
      integer j2
      integer j3
      integer k1
      integer k2
      integer k3
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
      k1 = i1
      k2 = i2
      k3 = i3

      j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
      j2 = min ( max ( k1, k2 ),
     &    min ( max ( k2, k3 ), max ( k3, k1 ) ) )
      j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )

      return
      end
      subroutine i4list_print ( n, first, list_num, list, title )

c*********************************************************************72
c
cc I4LIST_PRINT prints an I4LIST.
c
c  Discussion:
c
c    An I4LIST is a list of integers grouped into N segments.
c    An index vector locates the first entry of each segment.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of segments.
c
c    Input, integer FIRST(N+1), indexes the first entry
c    of each segment.
c
c    Input, integer LIST_NUM, the number of entries.
c
c    Input, integer LIST(LIST_NUM), the data.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer list_num
      integer n

      integer first(n+1)
      integer i
      integer j
      integer jhi
      integer jlo
      integer list(list_num)
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      do i = 1, n

        do jlo = first(i), first(i+1) - 1, 5
          jhi = min ( jlo + 4, first(i+1) - 1 )
          if ( jlo .eq. first(i) ) then
            write ( *, '(i5,a,5(2x,i8))' )
     &        i, ':', ( list(j), j = jlo, jhi )
          else
            write ( *, '(6x,5(2x,i8))' )
     &                ( list(j), j = jlo, jhi )
          end if
        end do

      end do

      return
      end
      subroutine i4mat_border_add ( m, n, table, table2 )

c*********************************************************************72
c
cc I4MAT_BORDER_ADD adds a "border" to an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    We suppose the input data gives values of a quantity on nodes
c    in the interior of a 2D grid, and we wish to create a new table
c    with additional positions for the nodes that would be on the
c    border of the 2D grid.
c
c                  0 0 0 0 0 0
c      * * * *     0 * * * * 0
c      * * * * --> 0 * * * * 0
c      * * * *     0 * * * * 0
c                  0 0 0 0 0 0
c
c    The illustration suggests the situation in which a 3 by 4 array
c    is input, and a 5 by 6 array is to be output.
c
c    The old data is shifted to its correct positions in the new array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input,integer TABLE(M,N), the table data.
c
c    Output, integer TABLE2(M+2,N+2), the augmented table data.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer table(m,n)
      integer table2(m+2,n+2)

      do j = 1, n + 2
        table2(1,j) = 0
      end do

      do j = 1, n + 2
        table2(m+2,j) = 0
      end do

      do i = 2, m + 1
        table2(i,1) = 0
      end do

      do i = 2, m + 1
        table2(i,n+2) = 0
      end do

      do j = 2, n + 1
        do i = 2, m + 1
          table2(i,j) = table(i-1,j-1)
        end do
      end do

      return
      end
      subroutine i4mat_border_cut ( m, n, table, table2 )

c*********************************************************************72
c
cc I4MAT_BORDER_CUT cuts the "border" of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    We suppose the input data gives values of a quantity on nodes
c    on a 2D grid, and we wish to create a new table corresponding only
c    to those nodes in the interior of the 2D grid.
c
c      0 0 0 0 0 0
c      0 * * * * 0    * * * *
c      0 * * * * 0 -> * * * *
c      0 * * * * 0    * * * *
c      0 0 0 0 0 0
c
c    The illustration suggests the situation in which a 5 by 6 array
c    is input, and a 3 by 4 array is to be output.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, integer TABLE(M,N), the table data.
c
c    Output, integer TABLE2(M-2,N-2), the new table data.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer table(m,n)
      integer table2(m-2,n-2)

      if ( m .le. 2 .or. n .le. 2 ) then
        return
      end if

      do j = 1, n - 2
        do i = 1, m - 2
          table2(i,j) = table(i+1,j+1)
        end do
      end do

      return
      end
      subroutine i4mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc I4MAT_COPY copies an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, integer A1(M,N), the matrix to copy.
c
c    Output, integer A2(M,N), the copy.
c
      implicit none

      integer m
      integer n

      integer a1(m,n)
      integer a2(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
        end do
      end do

      return
      end
      subroutine i4mat_elim ( m, n, a )

c*********************************************************************72
c
cc I4MAT_ELIM carries out exact Gauss elimination on an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input/output, integer A(M,N).  On input, the M by N matrix to
c    be Gauss eliminated.  On output, the Gauss-eliminated matrix.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer amax
      integer i
      integer icol(n)
      integer ifact
      integer i4_gcd
      integer imax
      integer imult
      integer irow(m)
      integer iswap
      integer j
      integer jcol
      integer jmult
      integer temp
c
c  Initialize the swap parity counter.
c
      iswap = 1
c
c  For each column JCOL...
c
      do jcol = 1, min ( m, n )
c
c  Find the maximum element in rows JCOL through M.
c
        amax = abs ( a(jcol,jcol) )
        imax = jcol

        do i = jcol + 1, m
          if ( amax .lt. abs ( a(i,jcol) ) ) then
            amax = abs ( a(i,jcol) )
            imax = i
          end if
        end do
c
c  If the maximum entry is nonzero, then...
c
        if ( amax .ne. 0 ) then
c
c  If the maximum entry does not occur in row JCOL, then swap rows.
c
          if ( imax .ne. jcol ) then
            iswap = - iswap
            do j = 1, n
              temp      = a(jcol,j)
              a(jcol,j) = a(imax,j)
              a(imax,j) = temp
            end do
          end if
c
c  Eliminate all nonzero entries in column JCOL, below the diagonal entry.
c
          do i = jcol + 1, m

            if ( a(i,jcol) .ne. 0 ) then

              jmult = a(i,jcol)
              imult = a(jcol,jcol)
              ifact = i4_gcd ( imult, jmult )
              imult = imult / ifact
              jmult = jmult / ifact

              do j = jcol, n
                a(i,j) = jmult * a(jcol,j) - imult * a(i,j)
              end do

            end if

          end do
c
c  Remove any row or column factors.
c
          call i4mat_red ( m, n, a, irow, icol )

        end if

      end do

      return
      end
      subroutine i4mat_flip_cols ( m, n, a )

c*********************************************************************72
c
cc I4MAT_FLIP_COLS swaps the columns of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an integer matrix.
c
c    To "flip" the columns of an I4MAT is to start with something like
c
c      11 12 13 14 15
c      21 22 23 24 25
c      31 32 33 34 35
c      41 42 43 44 45
c      51 52 53 54 55
c
c    and return
c
c      15 14 13 12 11
c      25 24 23 22 21
c      35 34 33 32 31
c      45 44 43 42 41
c      55 54 53 52 51
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), the matrix whose columns
c    are to be flipped.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b
      integer i
      integer j

      do j = 1, n / 2
        do i = 1, m
          b          = a(i,    j)
          a(i,    j) = a(i,n+1-j)
          a(i,n+1-j) = b
        end do
      end do

      return
      end
      subroutine i4mat_flip_rows ( m, n, a )

c*********************************************************************72
c
cc I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an integer matrix.
c
c    To "flip" the rows of an I4MAT is to start with something like
c
c      11 12 13 14 15
c      21 22 23 24 25
c      31 32 33 34 35
c      41 42 43 44 45
c      51 52 53 54 55
c
c    and return
c
c      51 52 53 54 55
c      41 42 43 44 45
c      31 32 33 34 35
c      21 22 23 24 25
c      11 12 13 14 15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), the matrix whose rows
c    are to be flipped.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b
      integer i
      integer j

      do i = 1, m / 2
        do j = 1, n
          b          = a(    i,j)
          a(    i,j) = a(m+1-i,j)
          a(m+1-i,j) = b
        end do
      end do

      return
      end
      subroutine i4mat_histogram ( m, n, a, histo_num, histo_gram )

c*********************************************************************72
c
cc I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    It is assumed that the entries in the vector A are nonnegative.
c    Only values between 0 and HISTO_NUM will be histogrammed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of A.
c
c    Input, integer A(M,N), the array to examine.
c
c    Input, integer HISTO_NUM, the maximum value for which a
c    histogram entry will be computed.
c
c    Output, integer HISTO_GRAM(0:HISTO_NUM), contains the
c    number of entries of A with the values of 0 through HISTO_NUM.
c
      implicit none

      integer histo_num
      integer m
      integer n

      integer a(m,n)
      integer histo_gram(0:histo_num)
      integer i
      integer j

      do i = 0, histo_num
        histo_gram(i) = 0
      end do

      do j = 1, n
        do i = 1, m

          if ( 0 .le. a(i,j) .and. a(i,j) .le. histo_num ) then
            histo_gram(a(i,j)) = histo_gram(a(i,j)) + 1
          end if

        end do
      end do

      return
      end
      subroutine i4mat_indicator ( m, n, table )

c*********************************************************************72
c
cc I4MAT_INDICATOR sets up an "indicator" I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    The value of each entry suggests its location, as in:
c
c      11  12  13  14
c      21  22  23  24
c      31  32  33  34
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Output, integer TABLE(M,N), the table.
c
      implicit none

      integer m
      integer n

      integer fac
      integer i
      integer i4_log_10
      integer j
      integer table(m,n)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, m
        do j = 1, n
          table(i,j) = fac * i + j
        end do
      end do

      return
      end
      subroutine i4mat_l1_inverse ( n, a, b )

c*********************************************************************72
c
cc I4MAT_L1_INVERSE inverts a unit lower triangular I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c    A unit lower triangular matrix is a matrix with only 1's on the main
c    diagonal, and only 0's above the main diagonal.
c
c    The inverse of an integer unit lower triangular matrix is also
c    an integer unit lower triangular matrix.
c
c    This routine can invert a matrix in place, that is, with no extra
c    storage.  If the matrix is stored in A, then the call
c
c      call i4mat_l1_inverse ( n, a, a )
c
c    will result in A being overwritten by its inverse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, number of rows and columns in the matrix.
c
c    Input, integer  A(N,N), the unit lower triangular matrix.
c
c    Output, integer B(N,N), the inverse matrix.
c
      implicit none

      integer n

      integer a(n,n)
      integer b(n,n)
      integer i
      integer j
      integer k

      do i = 1, n

        do j = 1, i - 1
          b(i,j) = 0
          do k = 1, i - 1
            b(i,j) = b(i,j) - a(i,k) * b(k,j)
          end do
        end do

        b(i,i) = 1
        do j = i + 1, n
          b(i,j) = 0
        end do

      end do

      return
      end
      function i4mat_max ( m, n, a )

c*********************************************************************72
c
cc I4MAT_MAX returns the maximum of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the M by N matrix.
c
c    Output, integer I4MAT_MAX, the maximum entry of A.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i4_huge
      integer i4mat_max
      integer j

      i4mat_max = - i4_huge ( 1 )

      do j = 1, n
        do i = 1, m
          i4mat_max = max ( i4mat_max, a(i,j) )
        end do
      end do

      return
      end
      subroutine i4mat_max_index ( m, n, a, i_max, j_max )

c*********************************************************************72
c
cc I4MAT_MAX_INDEX returns the location of the maximum of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the M by N matrix.
c
c    Output, integer I_MAX, J_MAX, the indices of the
c    maximum entry of A.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i_max
      integer j
      integer j_max

      i_max = -1;
      j_max = -1;

      do j = 1, n
        do i = 1, m
          if ( i .eq. 1 .and. j .eq. 1 ) then
            i_max = i
            j_max = j
          else if ( a(i_max,j_max) .lt. a(i,j) ) then
            i_max = i
            j_max = j
          end if
        end do
      end do

      return
      end
      function i4mat_min ( m, n, a )

c*********************************************************************72
c
cc I4MAT_MIN returns the minimum of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the M by N matrix.
c
c    Output, integer I4MAT_MIN, the minimum entry of A.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i4_huge
      integer i4mat_min
      integer j

      i4mat_min = i4_huge ( 1 )

      do j = 1, n
        do i = 1, m
          i4mat_min = min ( i4mat_min, a(i,j) )
        end do
      end do

      return
      end
      subroutine i4mat_min_index ( m, n, a, i_min, j_min )

c*********************************************************************72
c
cc I4MAT_MIN_INDEX returns the location of the minimum of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the M by N matrix.
c
c    Output, integer I_MIN, J_MIN, the indices of the
c    minimum entry of A.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i_min
      integer j
      integer j_min

      i_min = -1
      j_min = -1

      do j = 1, n
        do i = 1, m
          if ( i .eq. 1 .and. j .eq. 1 ) then
            i_min = i
            j_min = j
          else if ( a(i,j) .lt. a(i_min,j_min) ) then
            i_min = i
            j_min = j
          end if
        end do
      end do

      return
      end
      subroutine i4mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc I4MAT_MM multiplies two I4MAT's.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, integer A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, integer C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      integer a(n1,n2)
      integer b(n2,n3)
      integer c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      return
      end
      subroutine i4mat_perm ( n, a, p )

c*********************************************************************72
c
cc I4MAT_PERM permutes the rows and columns of a square I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c   07 June 2010
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, integer A(N,N).
c    On input, the matrix to be permuted.
c    On output, the permuted matrix.
c
c    Input, integer P(N), the permutation.  P(I) is the new
c    number of row and column I.
c
      implicit none

      integer n

      integer a(n,n)
      integer base
      parameter ( base = 1 )
      integer i
      integer i1
      integer ierror
      integer is
      integer it
      integer j
      integer j1
      integer j2
      integer k
      integer lc
      integer nc
      integer p(n)
      integer t

      call perm_check ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4MAT_PERM - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if

      call perm_cycle ( n, p, is, nc, 1 )

      do i = 1, n

        i1 = - p(i)

        if ( 0 .lt. i1 ) then

          lc = 0

10        continue

            i1 = p(i1)
            lc = lc + 1

            if ( i1 .le. 0 ) then
              go to 20
            end if

          go to 10

20        continue

          i1 = i

          do j = 1, n

            if ( p(j) .le. 0 ) then

              j2 = j
              k = lc

30            continue

                j1 = j2
                it = a(i1,j1)

40              continue

                  i1 = abs ( p(i1) )
                  j1 = abs ( p(j1) )

                  t        = a(i1,j1)
                  a(i1,j1) = it
                  it       = t

                  if ( j1 .ne. j2 ) then
                    go to 40
                  end if

                  k = k - 1

                  if ( i1 .eq. i ) then
                    go to 50
                  end if

                go to 40

50              continue

                j2 = abs ( p(j2) )

                if ( k .eq. 0 ) then
                  go to 60
                end if

              go to 30

60            continue

            end if

          end do

        end if

      end do
c
c  Restore the positive signs of the data.
c
      do i = 1, n
        p(i) = abs ( p(i) )
      end do

      return
      end
      subroutine i4mat_perm_uniform ( n, a, seed )

c*********************************************************************72
c
cc I4MAT_PERM_UNIFORM selects a random permutation of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c    The matrix is assumed to be square.  A single permutation is
c    applied to both rows and columns.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns
c    in the array.
c
c    Input/output, integer A(N,N), the array to be permuted.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
      implicit none

      integer n

      integer a(n,n)
      integer i
      integer i4_uniform_ab
      integer j
      integer k1
      integer k2
      integer seed
      integer t
c
c  Permute the rows and columns together.
c
      do k1 = 1, n

        k2 = i4_uniform_ab ( k1, n, seed )

        do j = 1, n
          t       = a(k2,j)
          a(k2,j) = a(k1,j)
          a(k1,j) = t
        end do

        do i = 1, n
          t       = a(i,k2)
          a(i,k2) = a(i,k1)
          a(i,k1) = t
        end do

      end do

      return
      end
      subroutine i4mat_perm2_uniform ( m, n, a, seed )

c*********************************************************************72
c
cc I4MAT_PERM2_UNIFORM selects a random permutation of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c    The matrix may be rectangular.  Separate permutations are
c    applied to the rows and columns.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), the array to be permuted.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i4_uniform_ab
      integer i2
      integer j
      integer j2
      integer seed
      integer t
c
c  Permute the rows.
c
      do i = 1, m
        i2 = i4_uniform_ab ( i, m, seed )
        do j = 1, n
          t       = a(i2,j)
          a(i2,j) = a(i,j)
          a(i,j)  = t
        end do
      end do
c
c  Permute the columns.
c
      do j = 1, n
        j2 = i4_uniform_ab ( j, n, seed )
        do i = 1, m
          t = a(i,j2)
          a(i,j2) = a(i,j)
          a(i,j)  = t
        end do
      end do

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the matrix to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      character*(*) title

      ilo = 1
      ihi = m
      jlo = 1
      jhi = n

      call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

      return
      end
      subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc I4MAT_PRINT_SOME prints some of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
      character*(8) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8)' ) j
        end do

        write ( *, '(''  Col '',10a8)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine i4mat_red ( m, n, a, row, col )

c*********************************************************************72
c
cc I4MAT_RED divides out common factors in a row or column of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in the matrix.
c
c    Input, integer N, the number of columns in the matrix.
c
c    Input/output, integer A(M,N), on input, the M by N matrix
c    to be reduced.  On output, A has been reduced.  The greatest common
c    factor in any row or column is 1.
c
!    Output, integer ROW(M), the row factors that were divided out.
!
!    Output, integer COL(N), the column factors that were divided
!    out.
!
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer arow(n)
      integer col(n)
      integer factor
      integer i
      integer j
      integer row(m)

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMAT_RED - Fatal error!'
        write ( *, '(a)' ) '  M must be greater than 0.'
        write ( *, '(a,i8)' ) '  Input M = ', m
        stop
      end if

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMAT_RED - Fatal error!'
        write ( *, '(a)' ) '  N must be greater than 0.'
        write ( *, '(a,i8)' ) '  Input N = ', n
        stop
      end if
c
c  Remove factors common to a column.
c
      do j = 1, n
        call i4vec_red ( m, a(1,j), factor )
        col(j) = factor
      end do
c
c  Remove factors common to a row.
c
      do i = 1, m
        do j = 1, n
          arow(j) = a(i,j)
        end do
        call i4vec_red ( n, arow, factor )
        row(i) = factor
      end do

      return
      end
      subroutine i4mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    39 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      character * ( * ) title

      call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi,
     &   jhi, title )

c*********************************************************************72
c
cc I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
      character*8 ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer  j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )  title

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(i8)' ) i
        end do

        write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
        write ( *, '(a)' ) '  Col'
        write ( *, '(a)' ) ' '

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

          do i2 = 1, inc

            i = i2lo - 1 + i2

            write ( ctemp(i2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      subroutine i4mat_u1_inverse ( n, a, b )

c*********************************************************************72
c
cc I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    A unit upper triangular matrix is a matrix with only 1's on the main
c    diagonal, and only 0's below the main diagonal.
c
c    The inverse of an integer unit upper triangular matrix is also
c    an integer unit upper triangular matrix.
c
c    This routine can invert a matrix in place, that is, with no extra
c    storage.  If the matrix is stored in A, then the call
c
c      call i4mat_u1_inverse ( n, a, a )
c
c    will result in A being overwritten by its inverse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, number of rows and columns in the matrix.
c
c    Input, integer A(N,N), the unit upper triangular matrix.
c
c    Output, integer B(N,N), the inverse matrix.
c
      implicit none

      integer n

      integer a(n,n)
      integer b(n,n)
      integer dot
      integer i
      integer j
      integer k

      do j = n, 1, -1

        b(j,j) = 1
        do i = j + 1, n
          b(i,j) = 0
        end do

        do i = j - 1, 1, -1
          dot = 0
          do k = i + 1, j
            dot = dot + a(i,k) * b(k,j)
          end do
          b(i,j) = - dot
        end do

      end do

      return
      end
      subroutine i4mat_uniform_ab ( m, n, a, b, seed, x )

c*********************************************************************72
c
cc I4MAT_UNIFORM_AB returns a scaled pseudorandom I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    The pseudorandom numbers should be uniformly distributed
c    between A and B.
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
c    Input, integer M, N, the row and column dimensions of the matrix.
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer X(M,N), a matrix of values between A and B.
c
      implicit none

      integer m
      integer n

      integer a
      integer b
      integer i
      integer j
      integer k
      real r
      integer seed
      integer value
      integer x(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4MAT_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

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
     &      +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
          value = nint ( r )

          value = max ( value, min ( a, b ) )
          value = min ( value, max ( a, b ) )

          x(i,j) = value

        end do
      end do

      return
      end
      subroutine i4mat_zero ( m, n, a )

c*********************************************************************72
c
cc I4MAT_ZERO zeroes out an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the row and column dimensions of the matrix.
c
c    Output, integer A(M,N), a matrix of zeroes.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = 0
        end do
      end do

      return
      end
      subroutine i4row_compare ( m, n, a, i, j, isgn )

c*********************************************************************72
c
cc I4ROW_COMPARE compares two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c  Example:
c
c    Input:
c
c    M = 3, N = 4, I = 2, J = 3
c
c    A = (
c      1  2  3  4
c      5  6  7  8
c      9 10 11 12 )
c
c    Output:
c
c      ISGN = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of M rows of vectors
c    of length N.
c
c    Input, integer I, J, the rows to be compared.
c    I and J must be between 1 and M.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, row I .lt. row J,
c     0, row I = row J,
c    +1, row J .lt. row I.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer isgn
      integer j
      integer k
c
c  Check that I and J are legal.
c
      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is less than 1.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      else if ( m .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is out of bounds.'
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      if ( j .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is less than 1.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      else if ( m .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is out of bounds.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      isgn = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. n ) then

        if ( a(i,k) .lt. a(j,k) ) then
          isgn = -1
          return
        else if ( a(j,k) .lt. a(i,k) ) then
          isgn = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4row_find_item ( m, n, a, item, row, col )

c*********************************************************************72
c
cc I4ROW_FIND_ITEM searches the rows of an I4ROW for a given value.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), the table to search.
c
c    Input, integer ITEM, the value to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM.  The search
c    is conducted by rows.  If the item is not found, then
c    ROW = COL = -1.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer item
      integer j
      integer row

      row = -1
      col = -1

      do i = 1, m
        do j = 1, n
          if ( a(i,j) .eq. item ) then
            row = i
            col = j
            return
          end if
        end do
      end do

      return
      end
      subroutine i4row_find_pair_wrap ( m, n, a, item1, item2, row,
     &  col )

c*********************************************************************72
c
cc I4ROW_FIND_PAIR_WRAP searches rows of an I4ROW for a pair of items.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c    The items must occur consecutively, with ITEM1 occurring
c    first.  However, wrapping is allowed.  That is, if ITEM1
c    occurs in the last column, and ITEM2 in the first, this
c    is also regarded as a match.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), the table to search.
c
c    Input, integer ITEM1, ITEM2, the values to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM1 followed immediately
c    by ITEM2.  The search is conducted by rows.  If the pair of
c    items is not found, then ROW = COL = -1.  If COL = N,
c    the ITEM1 occurs in column N and ITEM2 occurs in column 1.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer item1
      integer item2
      integer j
      integer jp1
      integer row

      row = -1
      col = -1

      do i = 1, m
        do j = 1, n

          if ( a(i,j) .eq. item1 ) then

            if ( j .lt. n ) then
              jp1 = j + 1
            else
              jp1 = 1
            end if

            if ( a(i,jp1) .eq. item2 ) then
              row = i
              col = j
              return
            end if

          end if

        end do
      end do

      return
      end
      subroutine i4row_max ( m, n, a, amax )

c*********************************************************************72
c
cc I4ROW_MAX returns the maximums of the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), the array to be examined.
c
c    Output, integer AMAX(M), the maximums of the rows
c    of the array.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer amax(m)
      integer i
      integer j

      do i = 1, m

        amax(i) = a(i,1)
        do j = 2, n
          if ( amax(i) .lt. a(i,j) ) then
            amax(i) = a(i,j)
          end if
        end do

      end do

      return
      end
      subroutine i4row_mean ( m, n, a, mean )

c*********************************************************************72
c
cc I4ROW_MEAN returns the means of the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of data.
c
c    Output, double precision MEAN(M), the mean of each row.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j
      double precision mean(m)

      do i = 1, m
        mean(i) = 0.0D+00
        do j = 1, n
          mean(i) = mean(i) + dble ( a(i,j) )
        end do
        mean(i) = mean(i) / dble ( n )
      end do

      return
      end
      subroutine i4row_min ( m, n, a, amin )

c*********************************************************************72
c
cc I4ROW_MIN returns the minimums of the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), the array to be examined.
c
c    Output, integer AMIN(M), the minimums of the rows.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer amin(m)
      integer i
      integer j

      do i = 1, m

        amin(i) = a(i,1)
        do j = 2, n
          if ( a(i,j) .lt. amin(i) ) then
            amin(i) = a(i,j)
          end if
        end do

      end do

      return
      end
      subroutine i4row_sort_a ( m, n, a )

c*********************************************************************72
c
cc I4ROW_SORT_A ascending sorts the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c    In lexicographic order, the statement "X .lt. Y", applied to two
c    vectors X and Y of length M, means that there is some index I, with
c    1 .le. I .le. M, with the property that
c
c      X(J) = Y(J) for J .lt. I,
c    and
c      X(I) .lt. Y(I).
c
c    In other words, X is less than Y if, at the first index where they
c    differ, the X value is less than the Y value.
c
c  Example:
c
c    Input:
c
c      M = 5, N = 3
c
c      A =
c        3  2  1
c        2  4  3
c        3  1  8
c        2  4  2
c        1  9  9
c
c    Output:
c
c      A =
c        1  9  9
c        2  4  2
c        2  4  3
c        3  1  8
c        3  2  1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of M rows of N-vectors.
c    On output, the rows of A have been sorted in ascending
c    lexicographic order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer indx
      integer isgn
      integer j

      if ( m .le. 1 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( m, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4row_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4row_compare ( m, n, a, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4row_sorted_unique ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4ROW_SORTED_UNIQUE keeps unique elements in an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c    The rows of the array may be ascending or descending sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), a sorted array, containing
c    M rows of data.  On output, the first UNIQUE_NUM rows
c    contain the unique rows.
c
c    Output, integer UNIQUE_NUM, the number of unique rows.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical equal
      integer i1
      integer i2
      integer j
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      i1 = 1

      do i2 = 2, m

        equal = .true.

        do j = 1, n
          if ( a(i1,j) .ne. a(i2,j) ) then
            equal = .false.
            go to 10
          end if
        end do

10      continue

        if ( .not. equal ) then
          i1 = i1 + 1
          do j = 1, n
            a(i1,j) = a(i2,j)
          end do
        end if

      end do

      unique_num = i1

      return
      end
      subroutine i4row_sorted_unique_count ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c    The rows of the array may be ascending or descending sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), a sorted array, containing
c    M rows of data.
c
c    Output, integer UNIQUE_NUM, the number of unique rows.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical equal
      integer i1
      integer i2
      integer j
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1
      i1 = 1

      do i2 = 2, m

        equal = .true.

        do j = 1, n
          if ( a(i1,j) .ne. a(i2,j) ) then
            equal = .false.
            go to 10
          end if
        end do

10      continue

        if ( .not. equal ) then
          unique_num = unique_num + 1
          i1 = i2
        end if

      end do

      return
      end
      subroutine i4row_swap ( m, n, a, i1, i2 )

c*********************************************************************72
c
cc I4ROW_SWAP swaps two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), an array of data.
c
c    Input, integer I1, I2, the two rows to swap.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i1
      integer i2
      integer row(n)
c
c  Check.
c
      if ( i1 .lt. 1 .or. m .lt. i1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I1 is out of range.'
        stop
      end if

      if ( i2 .lt. 1 .or. m .lt. i2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I2 is out of range.'
        stop
      end if

      if ( i1 .eq. i2 ) then
        return
      end if

      row(1:n)  = a(i1,1:n)
      a(i1,1:n) = a(i2,1:n)
      a(i2,1:n) = row(1:n)

      return
      end
      subroutine i4row_variance ( m, n, a, variance )

c*********************************************************************72
c
cc I4ROW_VARIANCE returns the variances of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), the array of data.
c
c    Output, double precision VARIANCE(M), the variance of each row.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j
      double precision mean
      double precision variance(m)

      if ( n .lt. 2 ) then

        do i = 1, m
          variance(i) = 0.0D+00
        end do

      else

        do i = 1, m

          mean = 0.0D+00
          do j = 1, n
            mean = mean + dble ( a(i,j) )
          end do
          mean = mean / dble ( n )

          variance(i) = 0.0D+00
          do j = 1, n
            variance(i) = variance(i) + ( dble ( a(i,j) ) - mean )**2
          end do

          variance(i) = variance(i) / dble ( n - 1 )

        end do

      end if

      return
      end
      subroutine i4vec_add ( n, a, b, c )

c*********************************************************************72
c
cc I4VEC_ADD computes C = A + B for I4VEC's.
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
c    28 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer A(N), the first vector.
c
c    Input, integer B(N), the second vector.
c
c    Output, integer C(N), the sum of the vectors.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer c(n)
      integer i

      do i = 1, n
        c(i) = a(i) + b(i)
      end do

      return
      end
      function i4vec_all_nonpositive ( n, a )

c*********************************************************************72
c
cc I4VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for I4VEC's.
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
c    08 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_ALL_NONPOSITIVE is TRUE if all entries
c    of A are less than or equal to 0.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_all_nonpositive

      do i = 1, n
        if ( 0 .lt. a(i) ) then
          i4vec_all_nonpositive = .false.
          return
        end if
      end do

      i4vec_all_nonpositive = .true.

      return
      end
      subroutine i4vec_amax ( n, a, aamax )

c*********************************************************************72
c
cc I4VEC_AMAX returns the largest magnitude in an I4VEC.
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
c    01 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be searched.
c
c    Output, integer AAMAX, the value of the entry of
c    largest magnitude.
c
      implicit none

      integer n

      integer a(n)
      integer aamax
      integer i

      if ( n .le. 0 ) then

        aamax = 0

      else

        aamax = abs ( a(1) )

        do i = 2, n
          aamax = max ( aamax, abs ( a(i) ) )
        end do

      end if

      return
      end
      subroutine i4vec_amax_index ( n, a, amax_index )

c*********************************************************************72
c
cc I4VEC_AMAX_INDEX returns the index of the largest magnitude in an I4VEC.
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
c    Input, integer A(N), the vector to be searched.
c
c    Output, integer AMAX_INDEX, the index of the entry
c    of largest magnitude.
c
      implicit none

      integer n

      integer a(n)
      integer aamax
      integer i
      integer amax_index

      if ( n .le. 0 ) then

        amax_index = 0

      else

        aamax = abs ( a(1) )
        amax_index = 1

        do i = 2, n

          if ( aamax .lt. abs ( a(i) ) ) then
            aamax = abs ( a(i) )
            amax_index = i
          end if

        end do

      end if

      return
      end
      subroutine i4vec_amin ( n, a, aamin )

c*********************************************************************72
c
cc I4VEC_AMIN returns the smallest magnitude in an I4VEC.
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
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries to be checked.
c
c    Input, integer A(N), the vector to be checked.
c
c    Output, integer AAMIN, the value of the smallest magnitude.
c
      implicit none

      integer n

      integer a(n)
      integer aamin
      integer i

      if ( n .le. 0 ) then

        aamin = 0

      else

        aamin = abs ( a(1) )

        do i = 2, n
          aamin = min ( aamin, abs ( a(i) ) )
        end do

      end if

      return
      end
      subroutine i4vec_amin_index ( n, a, amin_index )

c*********************************************************************72
c
cc I4VEC_AMIN_INDEX returns the index of the smallest magnitude in an I4VEC.
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
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries to be checked.
c
c    Input, integer A(N), the vector to be checked.
c
c    Output, integer AMIN_INDEX, the entry of the smallest
c    magnitude.
c
      implicit none

      integer n

      integer a(n)
      integer aamin
      integer i
      integer amin_index

      if ( n .le. 0 ) then

        amin_index = 0

      else

        aamin = a(1)
        amin_index = 1

        do i = 2, n

          if ( abs ( a(i) ) .lt. aamin ) then
            aamin = abs ( a(i) )
            amin_index = i
          end if

        end do

      end if

      return
      end
      subroutine i4vec_aminz ( n, a, aminz )

c*********************************************************************72
c
cc I4VEC_AMINZ returns the smallest nonzero magnitude in an I4VEC.
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
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries to be checked.
c
c    Input, integer A(N), the vector to be checked.
c
c    Output, integer AMINZ, the value of the smallest nonzero
c    magnitude.  If all entries are zero, AMINZ is 0.
c
      implicit none

      integer n

      integer a(n)
      integer aminz
      integer i
      integer iset

      aminz = 0
      iset = 0

      do i = 1, n

        if ( a(i) .ne. 0 ) then

          if ( iset .eq. 0 ) then
            aminz = abs ( a(i) )
            iset = 1
          else
            aminz = min ( aminz, abs ( a(i) ) )
          end if

        end if

      end do

      return
      end
      subroutine i4vec_aminz_index ( n, a, aminz_index )

c*********************************************************************72
c
cc I4VEC_AMINZ_INDEX returns the smallest nonzero magnitude in an I4VEC.
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
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries to be checked.
c
c    Input, integer A(N), the vector to be checked.
c
c    Output, integer AMINZ_INDEX, the entry of the smallest
c    nonzero magnitude.  If all entries are zero, AMINZ_INDEX is 0.
c
      implicit none

      integer n

      integer a(n)
      integer aminz
      integer i
      integer aminz_index

      aminz = 0
      aminz_index = 0

      do i = 1, n

        if ( a(i) .ne. 0 ) then

          if ( aminz_index .eq. 0 .or. abs ( a(i) ) .lt. aminz ) then
            aminz = abs ( a(i) )
            aminz_index = i
          end if

        end if

      end do

      return
      end
      function i4vec_any_lt ( n, a, b )

c*********************************************************************72
c
cc I4VEC_ANY_LT: ( any ( A < B ) ) for I4VEC's.
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
c    28 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer A(N), the first vector.
c
c    Input, integer B(N), the second vector.
c
c    Output, logical I4VEC_ANY_LT is TRUE if any entry
c    of A is less than the corresponding entry of B.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer i
      logical i4vec_any_lt

      do i = 1, n
        if ( a(i) .lt. b(i) ) then
          i4vec_any_lt = .true.
          return
        end if
      end do

      i4vec_any_lt = .false.

      return
      end
      function i4vec_any_negative ( n, a )

c*********************************************************************72
c
cc I4VEC_ANY_NEGATIVE: ( any A < 0 ) for I4VEC's.
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
c    08 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_ANY_NEGATIVE is TRUE if any entry is negative.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_any_negative

      do i = 1, n
        if ( a(i) .lt. 0 ) then
          i4vec_any_negative = .true.
          return
        end if
      end do

      i4vec_any_negative = .false.

      return
      end
      function i4vec_any_nonzero ( n, a )

c*********************************************************************72
c
cc I4VEC_ANY_NONZERO: ( any A nonzero ) for I4VEC's.
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
c    25 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_ANY_NONZERO is TRUE if any entry is nonzero.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_any_nonzero

      do i = 1, n
        if ( a(i) .ne. 0 ) then
          i4vec_any_nonzero = .true.
          return
        end if
      end do

      i4vec_any_nonzero = .false.

      return
      end
      subroutine i4vec_ascend_sub ( n, a, length, sub )

c*********************************************************************72
c
cc I4VEC_ASCEND_SUB computes the longest ascending subsequence of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The subsequence is required to be strictly increasing.
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
c    Input, integer N, the length of the vector.
c
c    Input, integer A(N), the vector to be examined.
c
c    Output, integer LENGTH, the length of the longest
c    increasing subsequence.
c
c    Output, integer SUB(N), contains in entries 1 through LENGTH
c    a longest increasing subsequence of A.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer k
      integer length
      integer sub(n)
      integer top(n)
      integer top_prev(n)

      do i = 1, n
        top(i) = 0
        top_prev(i) = 0
        sub(i) = 0
      end do

      if ( n .le. 0 ) then
        length = 0
        return
      end if

      length = 0

      do i = 1, n

        k = -1

        do j = 1, length
          if ( a(i) .le. a(top(j)) ) then
            k = j
            go to 10
          end if
        end do

10      continue

        if ( k .eq. -1 ) then
          length = length + 1
          k = length
        end if

        top(k) = i

        if ( 1 .lt. k ) then
          top_prev(i) = top(k-1)
        else
          top_prev(i) = 0
        end if

      end do
c
c  Extract the subsequence.
c
      j = top(length)
      sub(length) = a(j)

      do i = length - 1, 1, -1
        j = top_prev(j)
        sub(i) = a(j)
      end do

      return
      end
      function i4vec_ascends ( n, x )

c*********************************************************************72
c
cc I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    X = ( -8, 1, 2, 3, 7, 7, 9 )
c
c    I4VEC_ASCENDS = TRUE
c
c    The sequence is not required to be strictly ascending.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, logical I4VEC_ASCENDS, is TRUE if the entries of X ascend.
c
      implicit none

      integer n

      integer i
      logical i4vec_ascends
      integer x(n)

      i4vec_ascends = .false.

      do i = 1, n - 1
        if ( x(i+1) .lt. x(i) ) then
          return
        end if
      end do

      i4vec_ascends = .true.

      return
      end
      subroutine i4vec_axpy ( n, ia, x, incx, y, incy )

c*********************************************************************72
c
cc I4VEC_AXPY:  Y(I) := Y(I) + A * X(I).
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    If X and Y are simple vectors, then IAXPY is equivalent to:
c
c      DO I = 1, N
c        Y(I) = Y(I) + IA * X(I)
c      END DO
c
c    However, by using the increments correctly, IAXPY can also be used
c    to manipulate rows or columns of matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of X and Y.
c
c    Input, integer IA, the scalar value by which each entry
c    of X is multiplied before being added to Y.
c
c    Input, integer X(*), the vector, a multiple of which is to be
c    added to Y.
c
c    Input, integer INCX, the increment between successive
c    entries of X.
c
c    Input/output, integer Y(*).
c    On output, each entry of Y has been increased by
c    IA times the corresponding entry of X.
c
c    Input, integer INCY, the increment between successive
c    entries of Y.
c
      implicit none

      integer i
      integer ia
      integer incx
      integer incy
      integer indx
      integer indy
      integer n
      integer x(*)
      integer y(*)

      indx = 1
      indy = 1

      do i = 1, n

        y(indy) = y(indy) + ia * x(indx)

        indx = indx + incx
        indy = indy + incy

      end do

      return
      end
      subroutine i4vec_bracket ( n, a, xval, left, right )

c*********************************************************************72
c
cc I4VEC_BRACKET searches a sorted I4VEC for successive brackets of a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    If the values in the vector are thought of as defining intervals
c    on the number line, then this routine searches for the interval
c    containing the given value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of input array.
c
c    Input, integer A(N), an array that has been sorted
c    into ascending order.
c
c    Input, integer XVAL, a value to be bracketed.
c
c    Output, integer LEFT, RIGHT, the results of the search.
c    In the most common case, 1 <= LEFT .lt. LEFT + 1 = RIGHT <= N,
c    and A(LEFT) <= XVAL <= A(RIGHT).
c
c    Special cases:
c      Value is less than all data values:
c        LEFT = -1, RIGHT = 1, and XVAL .lt. A(RIGHT).
c      Value is greater than all data values:
c        LEFT = N, RIGHT = -1, and A(LEFT) .lt. XVAL.
c      Value is equal to a data value:
c        LEFT = RIGHT, and A(LEFT) = A(RIGHT) = XVAL.
c
      implicit none

      integer n

      integer a(n)
      integer high
      integer left
      integer low
      integer mid
      integer right
      integer xval
c
c  XVAL .lt. A(1).
c
      if ( xval .lt. a(1) ) then
        left = -1
        right = 1
c
c  A(N) .lt. XVAL.
c
      else if ( a(n) .lt. xval ) then
        left = n
        right = -1
c
c  N = 1
c
      else if ( n .eq. 1 ) then
        left = 1
        right = 1
c
c  A(1) <= XVAL <= A(N).
c
      else

        low = 1
        high = n - 1

10      continue

          mid = ( low + high ) / 2

          if ( high .lt. low ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'I4VEC_BRACKET - Fatal error!'
            write ( *, '(a)' ) '  Algorithm or data failure.'
            stop
          end if

          if ( a(mid) .eq. xval ) then
            left = mid
            right = mid
            go to 20
          else if ( a(mid+1) .eq. xval ) then
            left = mid + 1
            right = mid + 1
            go to 20
          else if ( a(mid) .lt. xval .and. xval .lt. a(mid+1) ) then
            left = mid
            right = mid + 1
            go to 20
          else if ( a(mid+1) .lt. xval ) then
            low = mid + 1
          else if ( xval .lt. a(mid) ) then
            high = mid - 1
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine i4vec_compare ( n, a1, a2, isgn )

c*********************************************************************72
c
cc I4VEC_COMPARE compares two I4VEC's.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c    The lexicographic ordering is used.
c
c  Example:
c
c    Input:
c
c      A1 = ( 2, 6, 2 )
c      A2 = ( 2, 8, 12 )
c
c    Output:
c
c      ISGN = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, integer A1(N), A2(N), the vectors to be compared.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, A1 < A2,
c     0, A1 = A2,
c    +1, A2 < A1.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer isgn
      integer k

      isgn = 0

      k = 1

10    continue

      if ( k .le. n ) then

        if ( a1(k) .lt. a2(k) ) then
          isgn = -1
          return
        else if ( a2(k) .lt. a1(k) ) then
          isgn = + 1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_COPY copies an I4VEC.
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
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer A1(N), the vector to be copied.
c
c    Output, integer A2(N), a copy of A1.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      subroutine i4vec_cum ( n, a, a_cum )

c*********************************************************************72
c
cc I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    Input:
c
c      A = (/ 1, 2, 3, 4 /)
c
c    Output:
c
c      A_CUM = (/ 1, 3, 6, 10 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be summed.
c
c    Output, integer A_CUM(N), the cumulative sum of the
c    entries of A.
c
      implicit none

      integer n

      integer a(n)
      integer a_cum(n)
      integer i

      a_cum(1) = a(1)

      do i = 2, n
        a_cum(i) = a_cum(i-1) + a(i)
      end do

      return
      end
      subroutine i4vec_cum0 ( n, a, a_cum )

c*********************************************************************72
c
cc I4VEC_CUM0 computes the cumulutive sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    This routine returns a vector of length N+1, with the first value
c    being 0.
c
c  Example:
c
c    Input:
c
c      A = (/ 1, 2, 3, 4 /)
c
c    Output:
c
c      A_CUM = (/ 0, 1, 3, 6, 10 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be summed.
c
c    Output, integer A_CUM(0:N), the cumulative sum of the
c    entries of A.
c
      implicit none

      integer n

      integer a(n)
      integer a_cum(0:n)
      integer i

      a_cum(0) = 0

      do i = 1, n
        a_cum(i) = a_cum(i-1) + a(i)
      end do

      return
      end
      function i4vec_descends ( n, x )

c*********************************************************************72
c
cc I4VEC_DESCENDS determines if an I4VEC is decreasing.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    X = ( 9, 7, 7, 3, 2, 1, -8 )
c
c    I4VEC_DESCENDS = TRUE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, logical I4VEC_DESCENDS, is TRUE if the entries of X descend.
c
      implicit none

      integer n

      integer i
      logical i4vec_descends
      integer x(n)

      i4vec_descends = .false.

      do i = 1, n - 1
        if ( x(i) .lt. x(i+1) ) then
          return
        end if
      end do

      i4vec_descends = .true.

      return
      end
      subroutine i4vec_direct_product ( factor_index, factor_order,
     &  factor_value, factor_num, point_num, x )

c*********************************************************************72
c
cc I4VEC_DIRECT_PRODUCT creates a direct product of I4VEC's.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    To explain what is going on here, suppose we had to construct
c    a multidimensional quadrature rule as the product of K rules
c    for 1D quadrature.
c
c    The product rule will be represented as a list of points and weights.
c
c    The J-th item in the product rule will be associated with
c      item J1 of 1D rule 1,
c      item J2 of 1D rule 2,
c      ...,
c      item JK of 1D rule K.
c
c    In particular,
c      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
c    and
c      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
c
c    So we can construct the quadrature rule if we can properly
c    distribute the information in the 1D quadrature rules.
c
c    This routine carries out that task for the abscissas X.
c
c    Another way to do this would be to compute, one by one, the
c    set of all possible indices (J1,J2,...,JK), and then index
c    the appropriate information.  An advantage of the method shown
c    here is that you can process the K-th set of information and
c    then discard it.
c
c  Example:
c
c    Rule 1:
c      Order = 4
c      X(1:4) = ( 1, 2, 3, 4 )
c
c    Rule 2:
c      Order = 3
c      X(1:3) = ( 10, 20, 30 )
c
c    Rule 3:
c      Order = 2
c      X(1:2) = ( 100, 200 )
c
c    Product Rule:
c      Order = 24
c      X(1:24) =
c        ( 1, 10, 100 )
c        ( 2, 10, 100 )
c        ( 3, 10, 100 )
c        ( 4, 10, 100 )
c        ( 1, 20, 100 )
c        ( 2, 20, 100 )
c        ( 3, 20, 100 )
c        ( 4, 20, 100 )
c        ( 1, 30, 100 )
c        ( 2, 30, 100 )
c        ( 3, 30, 100 )
c        ( 4, 30, 100 )
c        ( 1, 10, 200 )
c        ( 2, 10, 200 )
c        ( 3, 10, 200 )
c        ( 4, 10, 200 )
c        ( 1, 20, 200 )
c        ( 2, 20, 200 )
c        ( 3, 20, 200 )
c        ( 4, 20, 200 )
c        ( 1, 30, 200 )
c        ( 2, 30, 200 )
c        ( 3, 30, 200 )
c        ( 4, 30, 200 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR_INDEX, the index of the factor being
c    processed.  The first factor processed must be factor 1c
c
c    Input, integer FACTOR_ORDER, the order of the factor.
c
c    Input, integer FACTOR_VALUE(FACTOR_ORDER), the factor values
c    for factor FACTOR_INDEX.
c
c    Input, integer FACTOR_NUM, the number of factors.
c
c    Input, integer POINT_NUM, the number of elements in the
c    direct product.
c
c    Input/output, integer X(FACTOR_NUM,POINT_NUM), the elements of the
c    direct product, which are built up gradually.
c
c  Local Parameters:
c
c    Local, integer START, the first location of a block of
c    values to set.
c
c    Local, integer CONTIG, the number of consecutive values
c    to set.
c
c    Local, integer SKIP, the distance from the current value
c    of START to the next location of a block of values to set.
c
c    Local, integer REP, the number of blocks of values to set.
c
      implicit none

      integer factor_num
      integer factor_order
      integer point_num

      integer contig
      integer factor_index
      integer factor_value(factor_order)
      integer i
      integer j
      integer k
      integer rep
      integer skip
      integer start
      integer x(factor_num,point_num)

      save contig
      save rep
      save skip

      data contig / -1 /
      data rep / -1 /
      data skip / -1 /

      if ( factor_index .eq. 1 ) then
        contig = 1
        skip = 1
        rep = point_num
        do j = 1, point_num
          do i = 1, factor_num
            x(i,j) = 0
          end do
        end do
      end if

      rep = rep / factor_order
      skip = skip * factor_order

      do j = 1, factor_order

        start = 1 + ( j - 1 ) * contig

        do k = 1, rep
          do i = start, start + contig - 1
            x(factor_index,i) = factor_value(j)
          end do
          start = start + skip
        end do

      end do

      contig = contig * factor_order

      return
      end
      subroutine i4vec_direct_product2 ( factor_index, factor_order,
     &  factor_value, factor_num, point_num, w )

c*********************************************************************72
c
cc I4VEC_DIRECT_PRODUCT2 creates a direct product of I4VEC's.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    To explain what is going on here, suppose we had to construct
c    a multidimensional quadrature rule as the product of K rules
c    for 1D quadrature.
c
c    The product rule will be represented as a list of points and weights.
c
c    The J-th item in the product rule will be associated with
c      item J1 of 1D rule 1,
c      item J2 of 1D rule 2,
c      ...,
c      item JK of 1D rule K.
c
c    In particular,
c      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
c    and
c      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
c
c    So we can construct the quadrature rule if we can properly
c    distribute the information in the 1D quadrature rules.
c
c    This routine carries out the task involving the weights W.
c
c    Another way to do this would be to compute, one by one, the
c    set of all possible indices (J1,J2,...,JK), and then index
c    the appropriate information.  An advantage of the method shown
c    here is that you can process the K-th set of information and
c    then discard it.
c
c  Example:
c
c    Rule 1:
c      Order = 4
c      W(1:4) = ( 2, 3, 5, 7 )
c
c    Rule 2:
c      Order = 3
c      W(1:3) = ( 11, 13, 17 )
c
c    Rule 3:
c      Order = 2
c      W(1:2) = ( 19, 23 )
c
c    Product Rule:
c      Order = 24
c      W(1:24) =
c        ( 2 * 11 * 19 )
c        ( 3 * 11 * 19 )
c        ( 4 * 11 * 19 )
c        ( 7 * 11 * 19 )
c        ( 2 * 13 * 19 )
c        ( 3 * 13 * 19 )
c        ( 5 * 13 * 19 )
c        ( 7 * 13 * 19 )
c        ( 2 * 17 * 19 )
c        ( 3 * 17 * 19 )
c        ( 5 * 17 * 19 )
c        ( 7 * 17 * 19 )
c        ( 2 * 11 * 23 )
c        ( 3 * 11 * 23 )
c        ( 5 * 11 * 23 )
c        ( 7 * 11 * 23 )
c        ( 2 * 13 * 23 )
c        ( 3 * 13 * 23 )
c        ( 5 * 13 * 23 )
c        ( 7 * 13 * 23 )
c        ( 2 * 17 * 23 )
c        ( 3 * 17 * 23 )
c        ( 5 * 17 * 23 )
c        ( 7 * 17 * 23 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR_INDEX, the index of the factor being
c    processed.  The first factor processed must be factor 1c
c
c    Input, integer FACTOR_ORDER, the order of the factor.
c
c    Input, integer FACTOR_VALUE(FACTOR_ORDER), the factor values
c    for factor FACTOR_INDEX.
c
c    Input, integer FACTOR_NUM, the number of factors.
c
c    Input, integer POINT_NUM, the number of elements in the
c    direct product.
c
c    Input/output, integer W(POINT_NUM), the elements of the
c    direct product, which are built up gradually.
c
c  Local Parameters:
c
c    Local, integer START, the first location of a block of
c    values to set.
c
c    Local, integer CONTIG, the number of consecutive values to
c    set.
c
c    Local, integer SKIP, the distance from the current value
c    of START to the next location of a block of values to set.
c
c    Local, integer REP, the number of blocks of values to set.
c
      implicit none

      integer factor_num
      integer factor_order
      integer point_num

      integer contig
      integer factor_index
      integer factor_value(factor_order)
      integer i
      integer j
      integer k
      integer rep
      integer skip
      integer start
      integer w(point_num)

      save contig
      save rep
      save skip

      data contig / -1 /
      data rep / -1 /
      data skip / -1 /

      if ( factor_index .eq. 1 ) then
        contig = 1
        skip = 1
        rep = point_num
        do i = 1, point_num
          w(i) = 1
        end do
      end if

      rep = rep / factor_order
      skip = skip * factor_order

      do j = 1, factor_order

        start = 1 + ( j - 1 ) * contig

        do k = 1, rep
          do i = start, start + contig - 1
            w(i) = w(i) * factor_value(j)
          end do
          start = start + skip
        end do

      end do

      contig = contig * factor_order

      return
      end
      function i4vec_dot_product ( n, x, y )

c*********************************************************************72
c
cc I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
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
c    19 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), Y(N), the arrays.
c
c    Output, integer I4VEC_DOT_PRODUCT, the dot product of X and Y.
c
      implicit none

      integer n

      integer i
      integer i4vec_dot_product
      integer value
      integer x(n)
      integer y(n)

      value = 0
      do i = 1, n
        value = value + x(i) * y(i)
      end do

      i4vec_dot_product = value

      return
      end
      function i4vec_eq ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_EQ is true if every pair of entries in two I4VECs is equal.
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
c    13 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, integer A1(N), A2(N), two vectors to compare.
c
c    Output, logical I4VEC_EQ.
c    I4VEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
c    and .FALSE. otherwise.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      logical i4vec_eq

      i4vec_eq = .false.

      do i = 1, n
        if ( a1(i) .ne. a2(i) ) then
          return
        end if
      end do

      i4vec_eq = .true.

      return
      end
      function i4vec_even_all ( n, a )

c*********************************************************************72
c
cc I4VEC_EVEN_ALL is TRUE if all entries of an I4VEC are even.
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
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_EVEN_ALL, TRUE if all entries are even.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_even_all

      do i = 1, n
        if ( mod ( a(i), 2 ) .ne. 0 ) then
          i4vec_even_all = .false.
          return
        end if
      end do

      i4vec_even_all = .true.

      return
      end
      function i4vec_even_any ( n, a )

c*********************************************************************72
c
cc I4VEC_EVEN_ANY is TRUE if any entry of an I4VEC is even.
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
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_EVEN_ANY, TRUE if any entry is even.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_even_any

      do i = 1, n
        if ( mod ( a(i), 2 ) .eq. 0 ) then
          i4vec_even_any = .true.
          return
        end if
      end do

      i4vec_even_any = .false.

      return
      end
      subroutine i4vec_find ( n, a, value, location )

!*********************************************************************72
!
!! I4VEC_FIND finds the first occurrence of a value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input, integer A(N), the array.
!
!    Input, integer VALUE, the value being sought.
!
!    Output, integer LOCATION, the first location in A where VALUE occurs,
!    or -1 if VALUE never occurs.
!
      implicit none

      integer n

      integer a(n)
      integer i
      integer location
      integer value

      location = -1

      do i = 1, n

        if ( a(i) .eq. value ) then
          location = i
          return
        end if

      end do

      return
      end
      subroutine i4vec_first_index ( n, a, first_index )

c*********************************************************************72
c
cc I4VEC_FIRST_INDEX indexes the first occurrence of values in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
c    the first occurrence of the value A(I).
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
c    Input, integer N, the number of elements of A.
c
c    Input, integer A(N), the array.
c
c    Output, integer FIRST_INDEX(N), the first occurrence index.
c
      implicit none

      integer n

      integer a(n)
      integer first_index(n)
      integer i
      integer j

      do i = 1, n
        first_index(i) = -1
      end do

      do i = 1, n

        if ( first_index(i) .eq. -1 ) then

          first_index(i) = i

          do j = i + 1, n
            if ( a(i) .eq. a(j) ) then
             first_index(j) = i
            end if
          end do

        end if

      end do

      return
      end
      subroutine i4vec_frac ( n, a, k, frac )

c*********************************************************************72
c
cc I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Hoare's algorithm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), array to search.  On output,
c    the elements of A have been somewhat rearranged.
c
c    Input, integer K, the fractile to be sought.  If K = 1, the
c    minimum entry is sought.  If K = N, the maximum is sought.
c    Other values of K search for the entry which is K-th in size.
c    K must be at least 1, and no greater than N.
c
c    Output, integer FRAC, the value of the K-th fractile of A.
c
      implicit none

      integer n

      integer a(n)
      integer frac
      integer i
      integer iryt
      integer ix
      integer j
      integer k
      integer left
      integer t

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
        stop
      end if

      if ( k .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
        stop
      end if

      if ( n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
        stop
      end if

      left = 1
      iryt = n

10    continue

        if ( iryt .le. left ) then
          frac = a(k)
          go to 60
        end if

        ix = a(k)
        i = left
        j = iryt

20      continue

          if ( j .lt. i ) then

            if ( j .lt. k ) then
              left = i
            end if

            if ( k .lt. i ) then
              iryt = j
            end if

            go to 50

          end if
c
c  Find I so that IX <= A(I).
c
30        continue

          if ( a(i) .lt. ix ) then
            i = i + 1
            go to 30
          end if
c
c  Find J so that A(J) <= IX.
c
40        continue

          if ( ix .lt. a(j) ) then
            j = j - 1
            go to 40
          end if

          if ( i .le. j ) then

            t    = a(i)
            a(i) = a(j)
            a(j) = t

            i = i + 1
            j = j - 1

          end if

        go to 20

50      continue

      go to 10

60    continue

      return
      end
      subroutine i4vec_gcd ( n, v, gcd )

c*********************************************************************72
c
cc I4VEC_GCD returns the greatest common divisor of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The value GCD returned has the property that it is the greatest integer
c    which evenly divides every entry of V.
c
c    The entries in V may be negative.
c
c    Any zero entries in V are ignored.  If all entries of V are zero,
c    GCD is returned as 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of V.
c
c    Input, integer V(N), the vector.
c
c    Output, integer GCD, the greatest common divisor of V.
c
      implicit none

      integer n

      integer gcd
      integer i
      integer i4_gcd
      integer v(n)

      gcd = 0

      do i = 1, n

        if ( v(i) .ne. 0 ) then
          if ( gcd .eq. 0 ) then
            gcd = abs ( v(i) )
          else
            gcd = i4_gcd ( gcd, v(i) )
          end if
        end if

      end do
c
c  If GCD is 0, that can only happen because all entries of V are zero.
c
      if ( gcd .eq. 0 ) then
        gcd = 1
      end if

      return
      end
      subroutine i4vec_heap_a ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_A reorders an I4VEC into an ascending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An ascending heap is an array A with the property that, for every index J,
c    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
c    2*J and 2*J+1 are legal).
c
c                  A(1)
c                /      \
c            A(2)         A(3)
c          /     \        /  \
c      A(4)       A(5)  A(6) A(7)
c      /  \       /   \
c    A(8) A(9) A(10) A(11)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input/output, integer A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ifree
      integer key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n / 2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = a(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the smaller of the two values,
c  and update M if necessary.
c
            if ( a(m+1) .lt. a(m) ) then
              m = m + 1
            end if

          end if
c
c  If the small descendant is smaller than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( key .le. a(m) ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot.
c
20      continue

        a(ifree) = key

      end do

      return
      end
      subroutine i4vec_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_D reorders an I4VEC into an descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    A descending heap is an array A with the property that, for every index J,
c    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
c    2*J and 2*J+1 are legal).
c
c                  A(1)
c                /      \
c            A(2)         A(3)
c          /     \        /  \
c      A(4)       A(5)  A(6) A(7)
c      /  \       /   \
c    A(8) A(9) A(10) A(11)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input/output, integer A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ifree
      integer key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n/2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = a(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the larger of the two values,
c  and update M if necessary.
c
            if ( a(m) .lt. a(m+1) ) then
              m = m + 1
            end if

          end if
c
c  If the large descendant is larger than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( a(m) .le. key ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot IFREE.
c
20      continue

        a(ifree) = key

      end do

      return
      end
      subroutine i4vec_heap_d_extract ( n, a, value )

c*********************************************************************72
c
cc I4VEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In other words, the routine finds the maximum value in the
c    heap, returns that value to the user, deletes that value from
c    the heap, and restores the heap to its proper form.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, 2001,
c    ISBN: 0262032937.
c
c  Parameters:
c
c    Input/output, integer N, the number of items in the heap.
c
c    Input/output, integer A(N), the heap.
c
c    Output, integer VALUE, the item of maximum value, which has
c    been removed from the heap.
c
      implicit none

      integer a(*)
      integer n
      integer value

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_HEAP_D_EXTRACT - Fatal error!'
        write ( *, '(a)' ) '  The heap is empty.'
        stop
      end if
c
c  Get the maximum value.
c
      value = a(1)

      if ( n .eq. 1 ) then
        n = 0
        return
      end if
c
c  Shift the last value down.
c
      a(1) = a(n)
c
c  Restore the heap structure.
c
      n = n - 1
      call i4vec_sort_heap_d ( n, a )

      return
      end
      subroutine i4vec_heap_d_insert ( n, a, value )

c*********************************************************************72
c
cc I4VEC_HEAP_D_INSERT inserts a new I4 into a descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 150.
c
c  Parameters:
c
c    Input/output, integer N, the number of items in the heap.
c
c    Input/output, integer A(N), the heap.
c
c    Input, integer VALUE, the value to be inserted.
c
      implicit none

      integer a(*)
      integer i
      integer n
      integer parent
      integer value

      n = n + 1
      i = n

10    continue

      if ( 1 .lt. i ) then

        parent = i / 2

        if ( value .le. a(parent) ) then
          go to 20
        end if

        a(i) = a(parent)
        i = parent

        go to 10

      end if

20    continue

      a(i) = value

      return
      end
      subroutine i4vec_heap_d_max ( n, a, val_max )

c*********************************************************************72
c
cc I4VEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 150.
c
c  Parameters:
c
c    Input, integer N, the number of items in the heap.
c
c    Input, integer A(N), the heap.
c
c    Output, integer VAL_MAX, the maximum value in the heap.
c
      implicit none

      integer n

      integer a(n)
      integer val_max

      val_max = a(1)

      return
      end
      function i4vec_index ( n, a, aval )

c*********************************************************************72
c
cc I4VEC_INDEX returns the location of the first occurrence of a given value.
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
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be searched.
c
c    Input, integer AVAL, the value to be indexed.
c
c    Output, integer I4VEC_INDEX, the first location in A which
c    has the value AVAL, or 0 if no such index exists.
c
      implicit none

      integer n

      integer a(n)
      integer aval
      integer i
      integer i4vec_index

      do i = 1, n
        if ( a(i) .eq. aval ) then
          i4vec_index = i
          return
        end if
      end do

      i4vec_index = 0

      return
      end
      subroutine i4vec_index_delete_all ( n, x, indx, xval )

c*********************************************************************72
c
cc I4VEC_INDEX_DELETE_ALL deletes a value in an indexed sorted I4VEC.
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
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer N, the size of the current list.
c
c    Input/output, integer X(N), the list.
c
c    Input/output, integer INDX(N), the sort index of the list.
c
c    Input, integer XVAL, the value to be sought.
c
      implicit none

      integer n

      integer equal
      integer equal1
      integer equal2
      integer get
      integer i
      integer indx(*)
      integer less
      integer more
      integer put
      integer x(*)
      integer xval

      if ( n .lt. 1 ) then
        n = 0
        return
      end if

      call i4vec_index_search ( n, x, indx, xval, less, equal, more )

      if ( equal .eq. 0 ) then
        return
      end if

      equal1 = equal

10    continue

        if ( equal1 .le. 1 ) then
          go to 20
        end if

        if ( x(indx(equal1-1)) .ne. xval ) then
          go to 20
        end if

        equal1 = equal1 - 1

      go to 10

20    continue

      equal2 = equal

30    continue

        if ( n .le. equal2 ) then
          go to 40
        end if

        if ( x(indx(equal2+1)) .ne. xval ) then
          go to 40
        end if

        equal2 = equal2 + 1

      go to 30

40    continue
c
c  Discard certain X values.
c
      put = 0

      do get = 1, n

        if ( x(get) .ne. xval ) then
          put = put + 1
          x(put) = x(get)
        end if

      end do

      do i = put + 1, n
        x(i) = 0
      end do
c
c  Adjust the INDX values.
c
      do equal = equal1, equal2
        do i = 1, n
          if ( indx(equal) .lt. indx(i) ) then
            indx(i) = indx(i) - 1
          end if
        end do
      end do
c
c  Discard certain INDX values.
c
      do i = equal1, equal1 + n - equal2 - 1
        indx(i) = indx(i-equal1+equal2+1)
      end do

      do i = equal1 + n - equal2, n
        indx(i) = 0
      end do
c
c  Adjust N.
c
      n = put

      return
      end
      subroutine i4vec_index_delete_dupes ( n, x, indx, n2, x2, indx2 )

c*********************************************************************72
c
cc I4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The output quantities N2, X2, and INDX2 are computed from the
c    input quantities by sorting, and eliminating duplicates.
c
c    The output arrays should be dimensioned of size N, unless the user
c    knows in advance what the value of N2 will be.
c
c    The output arrays may be identified with the input arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the input list.
c
c    Input, integer X(N), the list.
c
c    Input, integer INDX(N), the sort index of the list.
c
c    Output, integer N2, the number of unique entries in X.
c
c    Output, integer X2(N2), a copy of the list which has
c    been sorted, and made unique.
c
c    Output, integer INDX2(N2), the sort index of the new list.
c
      implicit none

      integer n

      integer i
      integer indx(n)
      integer indx2(n)
      integer n2
      integer n3
      integer x(n)
      integer x2(n)
      integer x3(n)

      i = 0
      n3 = 0

10    continue

        i = i + 1

        if ( n .lt. i ) then
          go to 20
        end if

        if ( 1 .lt. i ) then
          if ( x(indx(i)) .eq. x3(n3) ) then
            go to 10
          end if
        end if

        n3 = n3 + 1
        x3(n3) = x(indx(i))

      go to 10

20    continue
!
!  Copy data into output arrays.
!
      n2 = n3
      do i = 1, n2
        x2(i) = x3(i)
      end do

      call i4vec_indicator ( n2, indx2 )

      return
      end
      subroutine i4vec_index_insert_unique ( n, x, indx, xval )

c*********************************************************************72
c
cc I4VEC_INDEX_INSERT_UNIQUE inserts a unique I4 into an indexed sorted I4VEC.
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
c    25 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer N, the size of the current list.
c    If the input value XVAL does not already occur in X, then N is increased.
c
c    Input/output, integer X(N), the list.
c    If the input value XVAL does not already occur in X, then it is added
c    to X.
c
c    Input/output, integer INDX(N), the sort index of the list.
c    If the input value XVAL does not already occur in X, then INDX is updated.
c
c    Input, integer XVAL, the value which will be inserted into
c    the X vector if it is not there already.
c
      implicit none

      integer n

      integer equal
      integer i
      integer indx(*)
      integer less
      integer more
      integer x(*)
      integer xval

      if ( n .le. 0 ) then
        n = 1
        x(1) = xval
        indx(1) = 1
        return
      end if
c
c  Does XVAL already occur in X?
c
      call i4vec_index_search ( n, x, indx, xval, less, equal, more )

      if ( equal .eq. 0 ) then
        x(n+1) = xval
        do i = n, more, - 1
          indx(i+1) = indx(i)
        end do
        indx(more) = n + 1
        n = n + 1
      end if

      return
      end
      subroutine i4vec_index_search ( n, x, indx, xval, less, equal, 
     &  more )

c*********************************************************************72
c
cc I4VEC_INDEX_SEARCH searches for an I4 in an indexed sorted I4VEC.
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
c    25 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the current list.
c
c    Input, integer X(N), the list.
c
c    Input, integer INDX(N), the sort index of the list.
c
c    Input, integer XVAL, the value to be sought.
c
c    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
c    entries of X that are just less than, equal to, and just greater
c    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
c    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
c    is the greatest entry of X, then MORE is N+1.
c
      implicit none

      integer n

      integer equal
      integer hi
      integer indx(n)
      integer less
      integer lo
      integer mid
      integer more
      integer x(n)
      integer xhi
      integer xlo
      integer xmid
      integer xval

      if ( n .le. 0 ) then
        less = 0
        equal = 0
        more = 0
        return
      end if

      lo = 1
      hi = n
      xlo = x(indx(lo))
      xhi = x(indx(hi))

      if ( xval .lt. xlo ) then
        less = 0
        equal = 0
        more = 1
        return
      else if ( xval .eq. xlo ) then
        less = 0
        equal = 1
        more = 2
        return
      end if

      if ( xhi .lt. xval ) then
        less = n
        equal = 0
        more = n + 1
        return
      else if ( xval .eq. xhi ) then
        less = n - 1
        equal = n
        more = n + 1
        return
      end if

10    continue

        if ( lo + 1 .eq. hi ) then
          less = lo
          equal = 0
          more = hi
          go to 20
        end if

        mid = ( lo + hi ) / 2
        xmid = x(indx(mid))

        if ( xval .eq. xmid ) then
          equal = mid
          less = mid - 1
          more = mid + 1
          go to 20
        else if ( xval .lt. xmid ) then
          hi = mid
        else if ( xmid .lt. xval ) then
          lo = mid
        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec_index_sort_unique ( n, x, n2, x2, indx2 )

c*********************************************************************72
c
cc I4VEC_INDEX_SORT_UNIQUE creates a sorted unique index for an I4VEC.
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
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the current list.
c
c    Input, integer X(N), the list.
c
c    Output, integer N2, the number of unique elements in X.
c
c    Output, integer X2(N2), a list of the unique elements of X.
c
c    Output, integer INDX2(N2), the sort index of the list.
c
      implicit none

      integer n

      integer i
      integer indx2(n)
      integer n2
      integer x(n)
      integer x2(n)

      n2 = 0

      do i = 1, n
        call i4vec_index_insert_unique ( n2, x2, indx2, x(i) )
      end do

      do i = n2 + 1, n
        x2(i) = -1
        indx2(i) = -1
      end do

      return
      end
      subroutine i4vec_indexed_heap_d ( n, a, indx )

c*********************************************************************72
c
cc I4VEC_INDEXED_HEAP_D creates a descending heap from an indexed I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
c    each referencing an entry of the data vector.
c
c    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
c    we have:
c      A(INDX(2*J))   <= A(INDX(J))
c    and
c      A(INDX(2*J+1)) <= A(INDX(J))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the index array.
c
c    Input, integer A(*), the array that is indexed.
c
c    Input/output, integer INDX(N), the index array.
c    Each entry of INDX must be a valid index for the array A.
c    On output, the indices have been reordered into a descending heap.
c
      implicit none

      integer n

      integer a(*)
      integer i
      integer ifree
      integer indx(n)
      integer key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n / 2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = indx(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the larger of the two values,
c  and update M if necessary.
c
            if ( a(indx(m)) .lt. a(indx(m+1)) ) then
              m = m + 1
            end if

          end if
c
c  If the large descendant is larger than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( a(indx(m)) .le. a(key) ) then
            go to 20
          end if

          indx(ifree) = indx(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot IFREE.
c
        indx(ifree) = key

      end do

20    continue

      return
      end
      subroutine i4vec_indexed_heap_d_extract ( n, a, indx,
     &  indx_extract )

c*********************************************************************72
c
cc I4VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
c    each referencing an entry of the data vector.
c
c    The routine finds the maximum value in the heap, returns that value to the
c    user, deletes that value from the heap, and restores the heap to its
c    proper form.
c
c    Note that the argument N must be a variable, which will be decremented
c    before return, and that INDX will hold one less value on output than it
c    held on input.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, 2001,
c    ISBN: 0262032937,
c    LC: QA76.C662.
c
c  Parameters:
c
c    Input/output, integer N, the number of items in the
c    index vector.
c
c    Input, integer A(*), the data vector.
c
c    Input/output, integer INDX(N), the index vector.
c
c    Output, integer INDX_EXTRACT, the index in A of the item of
c    maximum value, which has now been removed from the heap.
c
      implicit none

      integer a(*)
      integer indx(*)
      integer indx_extract
      integer n

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!'
        write ( *, '(a)' ) '  The heap is empty.'
        stop
      end if
c
c  Get the index of the maximum value.
c
      indx_extract = indx(1)

      if ( n .eq. 1 ) then
        n = 0
        return
      end if
c
c  Shift the last index down.
c
      indx(1) = indx(n)
c
c  Restore the heap structure.
c
      n = n - 1
      call i4vec_indexed_heap_d ( n, a, indx )

      return
      end
      subroutine i4vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

c*********************************************************************72
c
cc I4VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
c    each referencing an entry of the data vector.
c
c    Note that the argument N must be a variable, and will be incremented before
c    return, and that INDX must be able to hold one more entry on output than
c    it held on input.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, 2001,
c    ISBN: 0262032937,
c    LC: QA76.C662.
c
c  Parameters:
c
c    Input/output, integer N, the number of items in the
c    index vector.
c
c    Input, integer A(*), the data vector.
c
c    Input/output, integer INDX(N), the index vector.
c
c    Input, integer INDX_INSERT, the index in A of the value
c    to be inserted into the heap.
c
      implicit none

      integer a(*)
      integer i
      integer indx(*)
      integer indx_insert
      integer n
      integer parent

      n = n + 1
      i = n

10    continue

      if ( 1 .lt. i ) then

        parent = i / 2

        if ( a(indx_insert) .le. a(indx(parent)) ) then
          go to 20
        end if

        indx(i) = indx(parent)
        i = parent

        go to 10

      end if

20    continue

      indx(i) = indx_insert

      return
      end
      subroutine i4vec_indexed_heap_d_max ( n, a, indx, indx_max )

c*********************************************************************72
c
cc I4VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
c    each referencing an entry of the data vector.
c
c    This is one of three functions needed to model a priority queue.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, 2001,
c    ISBN: 0262032937,
c    LC: QA76.C662.
c
c  Parameters:
c
c    Input, integer N, the number of items in the index vector.
c
c    Input, integer A(*), the data vector.
c
c    Input, integer INDX(N), the index vector.
c
c    Output, integer INDX_MAX, the index in A of the maximum value
c    in the heap.
c
      implicit none

      integer n

      integer a(*)
      integer indx(n)
      integer indx_max

      indx_max = indx(1)

      return
      end
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
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
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      subroutine i4vec_insert ( n, a, pos, value )

c*********************************************************************72
c
cc I4VEC_INSERT inserts a value into an I4VEC.
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
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the array on input.
c
c    Input/output, integer A(N+1), the array.  On input, A is
c    assumed to contain N entries.  On output, A actually contains N+1 entries.
c
c    Input, integer POS, the position to be assigned the new entry.
c    1 <= POS <= N+1.
c
c    Input, integer VALUE, the value to be inserted at the given
c    position.
c
      implicit none

      integer n

      integer a(n+1)
      integer i
      integer pos
      integer value

      if ( pos .lt. 1 .or. n + 1 .lt. pos ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_INSERT - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal insertion position = ', pos
        stop

      else

        do i = n + 1, pos + 1, -1
          a(i) = a(i-1)
        end do

        a(pos) = value

      end if

      return
      end
      subroutine i4vec_lcm ( n, v, lcm )

c*********************************************************************72
c
cc I4VEC_LCM returns the least common multiple of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The value LCM returned has the property that it is the smallest integer
c    which is evenly divisible by every element of V.
c
c    The entries in V may be negative.
c
c    If any entry of V is 0, then LCM is 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of V.
c
c    Input, integer V(N), the vector.
c
c    Output, integer LCM, the least common multiple of V.
c
      implicit none

      integer n

      integer i
      integer i4_lcm
      integer lcm
      integer v(n)

      lcm = 1

      do i = 1, n

        if ( v(i) .eq. 0 ) then
          lcm = 0
          return
        end if

        lcm = i4_lcm ( lcm, v(i) )

      end do

      return
      end
      subroutine i4vec_mask_print ( n, a, mask_num, mask, title )

c*********************************************************************72
c
cc I4VEC_MASK_PRINT prints a masked I4VEC.
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
c    05 July 2009
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
c    Input, integer MASK_NUM, the number of masked elements.
c
c    Input, integer MASK(MASK_NUM), the indices of the vector
c    to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer mask_num
      integer n

      integer a(n)
      integer i
      integer mask(mask_num)
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Masked vector printout:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      do i = 1, mask_num
        write ( *, '(2x,i8,a,1x,i8,2x,i10)' )
     &    i, ':', mask(i), a(mask(i))
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
      subroutine i4vec_max_index ( n, a, max_index )

c*********************************************************************72
c
cc I4VEC_MAX_INDEX computes the index of a maximum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    If more than one element has the maximum value, this routine returns
c    the index of the first such element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
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
c    Output, integer MAX_INDEX, the index of the largest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amax
      integer i
      integer max_index

      if ( n .le. 0 ) then

        max_index = 0

      else

        amax = a(1)
        max_index = 1

        do i = 2, n

          if ( amax .lt. a(i) ) then
            amax = a(i)
            max_index = i
          end if

        end do

      end if

      return
      end
      function i4vec_max_index_last ( n, x )

c*********************************************************************72
c
cc I4VEC_MAX_INDEX_LAST returns the last maximal element location in an I4VEC
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    X = ( 5, 1, 2, 5, 0, 5, 3 )
c
c    I4VEC_MAX_INDEX_LAST = 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, integer I4VEC_MAX_INDEX_LAST, the index of the
c    last element of X of maximal value.
c
      implicit none

      integer n

      integer i
      integer i4vec_max_index_last
      integer max_last
      integer x(n)

      i4vec_max_index_last = 0

      do i = 1, n
        if ( i .eq. 1 ) then
          i4vec_max_index_last = 1
          max_last = x(1)
        else if ( max_last .le. x(i) ) then
          i4vec_max_index_last = i
          max_last = x(i)
        end if
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
      subroutine i4vec_median ( n, a, median )

c*********************************************************************72
c
cc I4VEC_MEDIAN returns the median of an unsorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Hoare's algorithm is used.  The values of the vector are
c    rearranged by this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), the array to search.  On output,
c    the order of the elements of A has been somewhat changed.
c
c    Output, integer MEDIAN, the value of the median of A.
c
      implicit none

      integer n

      integer a(n)
      integer k
      integer median

      k = ( n + 1 ) / 2

      call i4vec_frac ( n, a, k, median )

      return
      end
      subroutine i4vec_merge_a ( na, a, nb, b, nc, c )

c*********************************************************************72
c
cc I4VEC_MERGE_A merges two ascending sorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The elements of A and B should be sorted in ascending order.
c
c    The elements in the output array C will also be in ascending order,
c    and unique.
c
c    The output vector C may share storage with A or B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, integer A(NA), the first sorted array.
c
c    Input, integer NB, the dimension of B.
c
c    Input, integer B(NB), the second sorted array.
c
c    Output, integer NC, the number of elements in the output
c    array.  Note that C should usually be dimensioned at least NA+NB in the
c    calling routine.
c
c    Output, integer C(NC), the merged unique sorted array.
c
      implicit none

      integer na
      integer nb

      integer a(na)
      integer b(nb)
      integer c(na+nb)
      integer d(na+nb)
      integer j
      integer ja
      integer jb
      integer na2
      integer nb2
      integer nc
      integer order

      na2 = na
      nb2 = nb

      ja = 0
      jb = 0
      nc = 0

      call i4vec_order_type ( na2, a, order )

      if ( order .lt. 0 .or. 2 .lt. order ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
        write ( *, '(a)') '  Input array A is not ascending sorted.'
        stop
      end if

      call i4vec_order_type ( nb2, b, order )

      if ( order .lt. 0 .or. 2 .lt. order ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
        write ( *, '(a)' ) '  Input array B is not ascending sorted.'
        stop
      end if

10    continue
c
c  If we've used up all the entries of A, stick the rest of B on the end.
c
        if ( na2 .le. ja ) then

          do j = 1, nb2 - jb
            jb = jb + 1
            if ( nc .eq. 0 .or. d(nc) .lt. b(jb) ) then
              nc = nc + 1
              d(nc) = b(jb)
            end if
          end do

          do j = 1, nc
            c(j) = d(j)
          end do

          go to 20
c
c  If we've used up all the entries of B, stick the rest of A on the end.
c
        else if ( nb2 .le. jb ) then

          do j = 1, na2 - ja
            ja = ja + 1
            if ( nc .eq. 0 .or. d(nc) .lt. a(ja) ) then
              nc = nc + 1
              d(nc) = a(ja)
            end if
          end do

          do j = 1, nc
            c(j) = d(j)
          end do

          go to 20
c
c  Otherwise, if the next entry of A is smaller, that's our candidate.
c
        else if ( a(ja+1) .le. b(jb+1) ) then

          ja = ja + 1
          if ( nc .eq. 0 .or. d(nc) .lt. a(ja) ) then
            nc = nc + 1
            d(nc) = a(ja)
          end if
c
c  ...or if the next entry of B is the smaller, consider that.
c
        else

          jb = jb + 1
          if ( nc .eq. 0 .or. d(nc) .lt. b(jb) ) then
            nc = nc + 1
            d(nc) = b(jb)
          end if

        end if

      go to 10

20    continue

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
      subroutine i4vec_min_index ( n, a, imin )

c*********************************************************************72
c
cc I4VEC_MIN_INDEX computes the index of the minimum element of an I4VEC.
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
c    06 July 2009
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
c    Output, integer IMIN, the index of the smallest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amin
      integer i
      integer imin

      if ( n .le. 0 ) then

        imin = 0

      else

        amin = a(1)
        imin = 1

        do i = 2, n

          if ( a(i) .lt. amin ) then
            amin = a(i)
            imin = i
          end if

        end do

      end if

      return
      end
      subroutine i4vec_min_mv ( m, n, u, v, w )

c*********************************************************************72
c
cc I4VEC_MIN_MV determines U(1:N) /\ V for vectors U and a single vector V.
c
c  Discussion:
c
c    For two vectors U and V, each of length M, we define
c
c      ( U /\ V ) (I) = min ( U(I), V(I) ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the dimension of the vectors.
c
c    Input, integer N, the number of vectors in U.
c
c    Input, integer U(M,N), N vectors, each of length M.
c
c    Input, integer V(M), a vector of length M.
c
c    Output, integer W(M,N), the value of U /\ W.
c
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      integer ( kind = 4 ) u(m,n)
      integer ( kind = 4 ) v(m)
      integer ( kind = 4 ) w(m,n)

      do j = 1, n
        do i = 1, m
          w(i,j) = min ( u(i,j), v(i) )
        end do
      end do

      return
      end
      function i4vec_nonzero_count ( n, a )

c*********************************************************************72
c
cc I4VEC_NONZERO_COUNT counts the nonzero entries in an I4VEC.
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
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input, integer A(N), an array.
c
c    Output, integer I4VEC_NONZERO_COUNT, the number of
c    nonzero entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_nonzero_count

      i4vec_nonzero_count = 0

      do i = 1, n
        if ( a(i) .ne. 0 ) then
          i4vec_nonzero_count = i4vec_nonzero_count + 1
        end if
      end do

      return
      end
      subroutine i4vec_nonzero_first ( n, x, nz, indx )

c*********************************************************************72
c
cc I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The routine preserves the ordering of the nonzero entries.  It counts
c    the nonzeros, and returns an index vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer X(N), the vector to be shifted.
c
c    Output, integer NZ, the number of nonzero entries in
c    the vector.
c
c    Output, integer INDX(N), contains the original location
c    of each entry.
c
      implicit none

      integer n

      integer indx(n)
      integer j
      integer k
      integer nz
      integer x(n)

      nz = 0

      do j = 1, n
        indx(j) = j
      end do

      j = 0

10    continue

      if ( j .lt. n ) then

        j = j + 1

        if ( x(j) .ne. 0 ) then

          nz = nz + 1

          if ( nz .ne. j ) then

            x(nz) = x(j)
            x(j) = 0

            k = indx(nz)
            indx(nz) = j
            indx(j) = k

          end if
        end if

        go to 10

      end if

      return
      end
      function i4vec_norm_l0 ( n, a )

c*********************************************************************72
c
cc I4VEC_NORM_L0 returns the l0 "norm" of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The l0 "norm" simply counts the number of nonzero entries in the vector.
c    It is not a true norm, but has some similarities to one.  It is useful
c    in the study of compressive sensing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, integer I4VEC_NORM_L0, the value of the norm.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_norm_l0
      integer value

      value = 0
      do i = 1, n
        if ( a(i) .ne. 0 ) then
          value = value + 1
        end if
      end do

      i4vec_norm_l0 = value

      return
      end
      function i4vec_odd_all ( n, a )

c*********************************************************************72
c
cc I4VEC_ODD_ALL is TRUE if all entries of an I4VEC are odd.
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
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_ODD_ALL, TRUE if all entries are odd.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_odd_all

      do i = 1, n
        if ( mod ( a(i), 2 ) .ne. 1 ) then
          i4vec_odd_all = .false.
          return
        end if
      end do

      i4vec_odd_all = .true.

      return
      end
      function i4vec_odd_any ( n, a )

c*********************************************************************72
c
cc I4VEC_ODD_ANY is TRUE if any entry of an I4VEC is odd.
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
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, logical I4VEC_ODD_ANY, TRUE if any entry is odd.
c
      implicit none

      integer n

      integer a(n)
      integer i
      logical i4vec_odd_any

      do i = 1, n
        if ( mod ( a(i), 2 ) .eq. 1 ) then
          i4vec_odd_any = .true.
          return
        end if
      end do

      i4vec_odd_any = .false.

      return
      end
      subroutine i4vec_order_type ( n, a, order )

c*********************************************************************72
c
cc I4VEC_ORDER_TYPE determines if I4VEC is (non)strictly ascending/descending.
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
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the array.
c
c    Input, integer A(N), the array to be checked.
c
c    Output, integer ORDER, order indicator:
c    -1, no discernable order;
c    0, all entries are equal;
c    1, ascending order;
c    2, strictly ascending order;
c    3, descending order;
c    4, strictly descending order.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer order
c
c  Search for the first value not equal to A(1).
c
      i = 1

10    continue

        i = i + 1

        if ( n .lt. i ) then
          order = 0
          return
        end if

        if ( a(1) .lt. a(i) ) then

          if ( i .eq. 2 ) then
            order = 2
          else
            order = 1
          end if

          go to 20

        else if ( a(i) .lt. a(1) ) then

          if ( i .eq. 2 ) then
            order = 4
          else
            order = 3
          end if

          go to 20

        end if

      go to 10
c
c  Now we have a "direction".  Examine subsequent entries.
c
20    continue

      if ( i .lt. n ) then

        i = i + 1

        if ( order .eq. 1 ) then

          if ( a(i) .lt. a(i-1) ) then
            order = -1
            go to 30
          end if

        else if ( order .eq. 2 ) then

          if ( a(i) .lt. a(i-1) ) then
            order = -1
            go to 30
          else if ( a(i) .eq. a(i-1) ) then
            order = 1
          end if

        else if ( order .eq. 3 ) then

          if ( a(i-1) .lt. a(i) ) then
            order = -1
            go to 30
          end if

        else if ( order .eq. 4 ) then

          if ( a(i-1) .lt. a(i) ) then
            order = -1
            go to 30
          else if ( a(i) .eq. a(i-1) ) then
            order = 3
          end if

        end if

        go to 20

      end if

30    continue

      return
      end
      function i4vec_pairwise_prime ( n, a )

c*********************************************************************72
c
cc I4VEC_PAIRWISE_PRIME checks whether an I4VEC's entries are pairwise prime.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Two positive integers I and J are pairwise prime if they have no common
c    factor greater than 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to check.
c
c    Input, integer A(N), the vector of integers.
c
c    Output, logical I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
c    is pairwise prime.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_gcd
      logical i4vec_pairwise_prime
      integer j

      i4vec_pairwise_prime = .false.

      do i = 1, n
        do j = i + 1, n
          if ( i4_gcd ( a(i), a(j) ) /= 1 ) then
            return
          end if
        end do
      end do

      i4vec_pairwise_prime = .true.

      return
      end
      subroutine i4vec_part ( n, nval, a )

c*********************************************************************72
c
cc I4VEC_PART partitions an integer NVAL into N nearly equal parts.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    Input:
c
c      N = 5, NVAL = 17
c
c    Output:
c
c      A = ( 4, 4, 3, 3, 3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer NVAL, the integer to be partitioned.
c    NVAL may be positive, zero, or negative.
c
c    Output, integer A(N), the partition of NVAL.  The entries of
c    A add up to NVAL.  The entries of A are either all equal, or
c    differ by at most 1.  The entries of A all have the same sign
c    as NVAL, and the "largest" entries occur first.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer nval

      do i = 1, n
        a(i) = 0
      end do

      if ( 0 .lt. nval ) then

        j = 1
        do i = 1, nval
          a(j) = a(j) + 1
          j = j + 1
          if ( n .lt. j ) then
            j = 1
          end if
        end do

      else if ( nval .lt. 0 ) then

        j = 1
        do i = nval, -1
          a(j) = a(j) - 1
          j = j + 1
          if ( n .lt. j ) then
            j = 1
          end if
        end do

      end if

      return
      end
      subroutine i4vec_part_quick_a ( n, a, l, r )

c*********************************************************************72
c
cc I4VEC_PART_QUICK_A reorders an I4VEC as part of a quick sort.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The routine reorders the entries of A.  Using A(1) as a key,
c    all entries of A that are less than or equal to the key will
c    precede the key which precedes all entries that are greater than the key.
c
c  Example:
c
c    Input:
c
c      N = 8
c
c      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
c
c    Output:
c
c      L = 3, R = 6
c
c      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
c            -------        -------
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of A.
c
c    Input/output, integer A(N).  On input, the array to be
c    checked.  On output, A has been reordered as described above.
c
c    Output, integer L, R, the indices of A that define the
c    three segments.
c    Let KEY = the input value of A(1).  Then
c    I <= L                 A(I) < KEY;
c         L < I < R         A(I) = KEY;
c                 R <= I    KEY < A(I).
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer key
      integer l
      integer m
      integer r
      integer t

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_PART_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        l = 0
        r = 2
        return
      end if

      key = a(1)
      m = 1
c
c  The elements of unknown size have indices between L+1 and R-1.
c
      l = 1
      r = n + 1

      do i = 2, n

        if ( key .lt. a(l+1) ) then
          r = r - 1
          t      = a(r)
          a(r)   = a(l+1)
          a(l+2) = t
        else if ( a(l+1) .eq. key ) then
          m = m + 1
          t      = a(m)
          a(m)   = a(l+1)
          a(l+1) = t
          l = l + 1
        else if ( a(l+1) .lt. key ) then
          l = l + 1
        end if

      end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
      do i = 1, l - m
        a(i) = a(i+m)
      end do
!
!  Out of bounds here, occasionally.
!
      l = l - m

      do i = 1, m
        a(l+i) = key
      end do

      return
      end
      subroutine i4vec_permute ( n, p, a )

c*********************************************************************72
c
cc I4VEC_PERMUTE permutes an I4VEC in place.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    This routine permutes an array of integer "objects", but the same
c    logic can be used to permute an array of objects of any arithmetic
c    type, or an array of objects of any complexity.  The only temporary
c    storage required is enough to store a single object.  The number
c    of data movements made is N + the number of cycles of order 2 or more,
c    which is never more than N + N/2.
c
c  Example:
c
c    Input:
c
c      N = 5
c      P = (   2,   4,   5,   1,   3 )
c      A = (   1,   2,   3,   4,   5 )
c
c    Output:
c
c      A    = (   2,   4,   5,   1,   3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer P(N), the permutation.  P(I) = J means
c    that the I-th element of the output array should be the J-th
c    element of the input array.
c
c    Input/output, integer A(N), the array to be permuted.
c
      implicit none

      integer n

      integer a(n)
      integer a_temp
      integer base
      parameter ( base = 1 )
      integer i
      integer ierror
      integer iget
      integer iput
      integer istart
      integer p(n)

      call perm_check ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if
c
c  Search for the next element of the permutation that has not been used.
c
      do istart = 1, n

        if ( p(istart) .lt. 0 ) then

          go to 20

        else if ( p(istart) .eq. istart ) then

          p(istart) = - p(istart)
          go to 20

        else

          a_temp = a(istart)
          iget = istart
c
c  Copy the new value into the vacated entry.
c
10        continue

            iput = iget
            iget = p(iget)

            p(iput) = - p(iput)

            if ( iget .lt. 1 .or. n .lt. iget ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
              write ( *, '(a)' ) '  An index is out of range.'
              write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
              stop
            end if

            if ( iget .eq. istart ) then
              a(iput) = a_temp
              go to 20
            end if

            a(iput) = a(iget)

          go to 10

        end if

20      continue

      end do
c
c  Restore the signs of the entries.
c
      do i = 1, n
        p(1:n) = - p(1:n)
      end do

      return
      end
      subroutine i4vec_permute_uniform ( n, a, seed )

c*********************************************************************72
c
cc I4VEC_PERMUTE_UNIFORM randomly permutes an I4VEC.
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
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input/output, integer A(N), the array to be permuted.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
      implicit none

      integer n

      integer a(n)
      integer p(n)
      integer seed

      call perm_uniform ( n, seed, p )

      call i4vec_permute ( n, a, p )

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
      subroutine i4vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc I4VEC_PRINT_PART prints "part" of an I4VEC.
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
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
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
          write ( *, '(2x,i8,a1,1x,i8)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a1,1x,i8)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '  ........  ........'
        i = n

        write ( *, '(2x,i8,a1,1x,i8)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a1,1x,i8)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(2x,i8,a1,1x,i8,a)' )
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      subroutine i4vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc I4VEC_PRINT_SOME prints "some" of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
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
c    17 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
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
          write ( *, '(i6,a,1x,2x,i10)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(i6,a,1x,2x,i10)' ) i, ':', a(i)
        end do
        write ( *, '(a)' ) '......  ..............'
        i = n
        write ( *, '(i6,a,1x,2x,i10)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i6,a,1x,2x,i10)' ) i, ':', a(i)
        end do
        i = max_print
        write ( *, '(i6,a,1x,i10,2x,a)' )
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      function i4vec_product ( n, a )

c*********************************************************************72
c
cc I4VEC_PRODUCT returns the product of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    PRODUCT function:
c
c      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
c
c    In MATLAB, this facility is offered by the built in
c    PROD function:
c
c      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2007
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
c    Output, integer I4VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_product

      i4vec_product = 1
      do i = 1, n
        i4vec_product = i4vec_product * a(i)
      end do

      return
      end
      subroutine i4vec_red ( n, a, factor )

c*********************************************************************72
c
cc I4VEC_RED divides out common factors in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    On output, the entries of A have no common factor
c    greater than 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer A(N), the vector to be reduced.
c
c    Output, integer FACTOR, the common factor that was divided
c    out.
c
      implicit none

      integer n

      integer a(n)
      integer factor
      integer i
      integer i4_gcd
c
c  Find the smallest nonzero value.
c
      factor = 0

      do i = 1, n

        if ( a(i) .ne. 0 ) then

          if ( factor .eq. 0 ) then
            factor = abs ( a(i) )
          else
            factor = min ( factor, abs ( a(i) ) )
          end if

        end if

      end do

      if ( factor .eq. 0 ) then
        return
      end if
c
c  Find the greatest common factor of the entire vector.
c
      do i = 1, n
        factor = i4_gcd ( a(i), factor )
      end do

      if ( factor .eq. 1 ) then
        return
      end if
c
c  Divide out the common factor.
c
      do i = 1, n
        a(i) = a(i) / factor
      end do

      return
      end
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11, 12, 13, 14, 15 ).
c
c    Output:
c
c      A = ( 15, 14, 13, 12, 11 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N), the array to be reversed.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer t

      do i = 1, n / 2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
      end do

      return
      end
      subroutine i4vec_rotate ( n, m, a )

c*********************************************************************72
c
cc I4VEC_ROTATE rotates an I4VEC in place.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    Input:
c
c      N = 5, M = 2
c      A = ( 1, 2, 3, 4, 5 )
c
c    Output:
c
c      A = ( 4, 5, 1, 2, 3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer M, the number of positions to the right that
c    each element should be moved.  Elements that shift pass position
c    N "wrap around" to the beginning of the array.
c
c    Input/output, integer A(N), the array to be rotated.
c
      implicit none

      integer n

      integer a(n)
      integer i4_modp
      integer iget
      integer iput
      integer istart
      integer m
      integer mcopy
      integer nset
      integer temp
c
c  Force M to be positive, between 0 and N-1.
c
      mcopy = i4_modp ( m, n )

      if ( mcopy .eq. 0 ) then
        return
      end if

      istart = 0
      nset = 0

10    continue

        istart = istart + 1

        if ( n .lt. istart ) then
          go to 40
        end if

        temp = a(istart)
        iget = istart
c
c  Copy the new value into the vacated entry.
c
20      continue

          iput = iget

          iget = iget - mcopy

          if ( iget .lt. 1 ) then
            iget = iget + n
          end if

          if ( iget .eq. istart ) then
            go to 30
          end if

          a(iput) = a(iget)
          nset = nset + 1

        go to 20

30      continue

        a(iput) = temp
        nset = nset + 1

        if ( n .le. nset ) then
          go to 40
        end if

      go to 10

40    continue

      return
      end
      subroutine i4vec_run_count ( n, a, run_count )

c*********************************************************************72
c
cc I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    A run is a sequence of equal values.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be examined.
c
c    Output, integer RUN_COUNT, the number of runs.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer run_count
      integer test

      run_count = 0

      if ( n .lt. 1 ) then
        return
      end if

      test = 0

      do i = 1, n

        if ( i .eq. 1 .or. a(i) .ne. test ) then
          run_count = run_count + 1
          test = a(i)
        end if

      end do

      return
      end
      subroutine i4vec_search_binary_a ( n, a, b, indx )

c*********************************************************************72
c
cc I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Binary search is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.9,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 26.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, integer A(N), the array to be searched.  A must
c    be sorted in ascending order.
c
c    Input, integer B, the value to be searched for.
c
c    Output, integer INDX, the result of the search.
c    -1, B does not occur in A.
c    I, A(I) = B.
c
      implicit none

      integer n

      integer a(n)
      integer b
      integer high
      integer indx
      integer low
      integer mid

      indx = - 1

      low = 1
      high = n

10    continue

      if ( low .le. high ) then

        mid = ( low + high ) / 2

        if ( a(mid) .eq. b ) then
          indx = mid
          go to 20
        else if ( a(mid) .lt. b ) then
          low = mid + 1
        else if ( b .lt. a(mid) ) then
          high = mid - 1
        end if

        go to 10

      end if

20    continue

      return
      end
      subroutine i4vec_search_binary_d ( n, a, b, indx )

c*********************************************************************72
c
cc I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Binary search is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.9,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 26.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, integer A(N), the array to be searched.  A must
c    be sorted in descending order.
c
c    Input, integer B, the value to be searched for.
c
c    Output, integer INDX, the result of the search.
c    -1, B does not occur in A.
c    I, A(I) = B.
c
      implicit none

      integer n

      integer a(n)
      integer b
      integer high
      integer indx
      integer low
      integer mid

      indx = - 1

      low = 1
      high = n

10    continue

      if ( low .le. high ) then

        mid = ( low + high ) / 2

        if ( a(mid) .eq. b ) then
          indx = mid
          go to 20
        else if ( b .lt. a(mid) ) then
          low = mid + 1
        else if ( a(mid) .lt. b ) then
          high = mid - 1
        end if

        go to 10

      end if

20    continue

      return
      end
      subroutine i4vec_sort_bubble_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
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
c    02 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j

      do i = 1, n - 1
        do j = i + 1, n
          if ( a(j) .lt. a(i) ) then
            call i4_swap ( a(i), a(j) )
          end if
        end do
      end do

      return
      end
      subroutine i4vec_sort_bubble_d ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_BUBBLE_D descending sorts an I4VEC using bubble sort.
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
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer k

      do i = 1, n - 1
        do j = i + 1, n
          if ( a(i) .lt. a(j) ) then
            k    = a(i)
            a(i) = a(j)
            a(j) = k
          end if
        end do
      end do

      return
      end
      subroutine i4vec_sort_heap_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
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
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if
c
c  1: Put A into descending heap form.
c
      call i4vec_heap_d ( n, a )
c
c  2: Sort A.
c
c  The largest object in the heap is in A(1).
c  Move it to position A(N).
c
      call i4_swap ( a(1), a(n) )
c
c  Consider the diminished heap of size N1.
c
      do n1 = n - 1, 2, -1
c
c  Restore the heap structure of A(1) through A(N1).
c
        call i4vec_heap_d ( n1, a )
c
c  Take the largest object from A(1) and move it to A(N1).
c
        call i4_swap ( a(1), a(n1) )

      end do

      return
      end
      subroutine i4vec_sort_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
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
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if
c
c  1: Put A into ascending heap form.
c
      call i4vec_heap_a ( n, a )
c
c  2: Sort A.
c
c  The smallest object in the heap is in A(1).
c  Move it to position A(N).
c
      call i4_swap ( a(1), a(n) )
c
c  Consider the diminished heap of size N1.
c
      do n1 = n - 1, 2, -1
c
c  Restore the heap structure of A(1) through A(N1).
c
        call i4vec_heap_a ( n1, a )
c
c  Take the smallest object from A(1) and move it to A(N1).
c
        call i4_swap ( a(1), a(n1) )

      end do

      return
      end
      subroutine i4vec_sort_heap_index_a ( n, a, indx )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(INDX(1:N)) is sorted,
c
c    or explicitly, by the call
c
c      call i4vec_permute ( n, indx, a )
c
c    after which A(1:N) is sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), an array to be index-sorted.
c
c    Output, integer INDX(N), the sort index.  The
c    I-th element of the sorted array is A(INDX(I)).
c
      implicit none

      integer n

      integer a(n)
      integer aval
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer j
      integer l

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, n
        indx(i) = i
      end do

      if ( n .eq. 1 ) then
        return
      end if

      l = n / 2 + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          aval = a(indxt)

        else

          indxt = indx(ir)
          aval = a(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then
            if ( a(indx(j)) .lt. a(indx(j+1)) ) then
              j = j + 1
            end if
          end if

          if ( aval .lt. a(indx(j)) ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

      return
      end
      subroutine i4vec_sort_heap_index_d ( n, a, indx )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(INDX(1:N)) is sorted,
c
c    or explicitly, by the call
c
c      call i4vec_permute ( n, indx, a )
c
c    after which A(1:N) is sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), an array to be index-sorted.
c
c    Output, integer INDX(N), the sort index.  The
c    I-th element of the sorted array is A(INDX(I)).
c
      implicit none

      integer n

      integer a(n)
      integer aval
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer j
      integer l

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, n
        indx(i) = i
      end do

      if ( n .eq. 1 ) then
        return
      end if

      l = n / 2 + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          aval = a(indxt)

        else

          indxt = indx(ir)
          aval = a(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then
            if ( a(indx(j+1)) .lt. a(indx(j)) ) then
              j = j + 1
            end if
          end if

          if ( a(indx(j)) .lt. aval ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

      return
      end
      subroutine i4vec_sort_insert_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
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
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.1,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 11.
c
c  Parameters:
c
c    Input, integer N, the number of items in the vector.
c    N must be positive.
c
c    Input/output, integer A(N).
c    On input, A contains data to be sorted.
c    On output, the entries of A have been sorted in ascending order.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer x

      do i = 2, n

        x = a(i)

        j = i - 1

10      continue

        if ( 1 .le. j ) then

          if ( a(j) .le. x ) then
            go to 20
          end if

          a(j+1) = a(j)
          j = j - 1

          go to 10

        end if

20      continue

        a(j+1) = x

      end do

      return
      end
      subroutine i4vec_sort_insert_d ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
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
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.1,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 11.
c
c  Parameters:
c
c    Input, integer N, the number of items in the vector.
c    N must be positive.
c
c    Input/output, integer A(N).
c    On input, A contains data to be sorted.
c    On output, the entries of A have been sorted in descending order.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer x

      do i = 2, n

        x = a(i)

        j = i - 1

10      continue

        if ( 1 .le. j ) then

          if ( x .le. a(j) ) then
            go to 20
          end if

          a(j+1) = a(j)
          j = j - 1
          go to 10

        end if

20      continue

        a(j+1) = x

      end do

      return
      end
      subroutine i4vec_sort_quick_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_QUICK_A ascending sorts an I4VEC using quick sort.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    Input:
c
c      N = 7
c
c      A = (/ 6, 7, 3, 2, 9, 1, 8 /)
c
c    Output:
c
c      A = (/ 1, 2, 3, 6, 7, 8, 9 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted.
c    On output, the array has been sorted.
c
      implicit none

      integer level_max
      parameter ( level_max = 30 )
      integer n

      integer a(n)
      integer base
      integer l_segment
      integer level
      integer n_segment
      integer rsave(level_max)
      integer r_segment

      if ( n .le. 1 ) then
        return
      end if

      level = 1
      rsave(level) = n + 1
      base = 1
      n_segment = n

10    continue
c
c  Partition the segment.
c
        call i4vec_part_quick_a ( n_segment, a(base), l_segment,
     &    r_segment )
c
c  If the left segment has more than one element, we need to partition it.
c
        if ( 1 .lt. l_segment ) then

          if ( level_max .lt. level ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'I4VEC_SORT_QUICK_A - Fatal error!'
            write ( *, '(a,i8)' )
     &        '  Exceeding recursion maximum of ', level_max
            stop
          end if

          level = level + 1
          n_segment = l_segment
          rsave(level) = r_segment + base - 1
c
c  The left segment and the middle segment are sorted.
c  Must the right segment be partitioned?
c
        else if ( r_segment .lt. n_segment ) then

          n_segment = n_segment + 1 - r_segment
          base = base + r_segment - 1
c
c  Otherwise, we back up a level if there is an earlier one.
c
        else

20        continue

            if ( level .le. 1 ) then
              return
            end if

            base = rsave(level)
            n_segment = rsave(level-1) - rsave(level)
            level = level - 1

            if ( 0 .lt. n_segment ) then
              go to 10
            end if

          go to 20

        end if

      go to 10

      return
      end
      subroutine i4vec_sort_shell_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_SHELL_A ascending sorts an I4VEC using Shell's sort.
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
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, an array to be sorted.
c    On output, the sorted array.
c
      implicit none

      integer n

      integer a(n)
      integer asave
      integer i
      integer ifree
      integer inc
      integer ipow
      integer j
      integer k
      integer maxpow

      if ( n .le. 1 ) then
        return
      end if
c
c  Determine the smallest MAXPOW so that
c    N <= ( 3**MAXPOW - 1 ) / 2
c
      maxpow = 1

10    continue

      if ( 3**maxpow .lt. 2 * n + 1 ) then
        maxpow = maxpow + 1
        go to 10
      end if

      if ( 1 .lt. maxpow ) then
        maxpow = maxpow - 1
      end if
c
c  Now sort groups of size ( 3**IPOW - 1 ) / 2.
c
      do ipow = maxpow, 1, -1

        inc = ( 3**ipow - 1 ) / 2
c
c  Sort the values with indices equal to K mod INC.
c
        do k = 1, inc
c
c  Insertion sort of the items with index
c  INC+K, 2*INC+K, 3*INC+K, ...
c
          do i = inc + k, n, inc

            asave = a(i)
            ifree = i
            j = i - inc

20          continue

              if ( j .lt. 1 ) then
                go to 30
              end if

              if ( a(j) .le. asave ) then
                go to 30
              end if

              ifree = j
              a(j+inc) = a(j)
              j = j - inc

            go to 20

30          continue

            a(ifree) = asave

          end do

        end do

      end do

      return
      end
      subroutine i4vec_sorted_undex ( x_num, x_val, x_unique_num,
     &  undx, xdnu )

c*********************************************************************72
c
cc I4VEC_SORTED_UNDEX returns unique sorted indexes for a sorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c    The goal of this routine is to determine a vector UNDX,
c    which points, to the unique elements of X, in sorted order,
c    and a vector XDNU, which identifies, for each entry of X, the index of
c    the unique sorted element of X.
c
c    This is all done with index vectors, so that the elements of
c    X are never moved.
c
c    Assuming X is already sorted, we examine the entries of X in order,
c    noting the unique entries, creating the entries of XDNU and
c    UNDX as we go.
c
c    Once this process has been completed, the vector X could be
c    replaced by a compressed vector XU, containing the unique entries
c    of X in sorted order, using the formula
c
c      XU(I) = X(UNDX(I)).
c
c    We could then, if we wished, reconstruct the entire vector X, or
c    any element of it, by index, as follows:
c
c      X(I) = XU(XDNU(I)).
c
c    We could then replace X by the combination of XU and XDNU.
c
c    Later, when we need the I-th entry of X, we can locate it as
c    the XDNU(I)-th entry of XU.
c
c    Here is an example of a vector X, the unique sort and
c    inverse unique sort vectors and the compressed unique sorted vector.
c
c      I    X    XU  Undx  Xdnu
c    ----+----+----+-----+-----+
c      1 | 11 |  11    1     1
c      2 | 11 |  22    5     1
c      3 | 11 |  33    8     1
c      4 | 11 |  55    9     1
c      5 | 22 |              2
c      6 | 22 |              2
c      7 | 22 |              2
c      8 | 33 |              3
c      9 | 55 |              4
c
c    UNDX(3) = 8 means that unique sorted item(3) is at X(8).
c    XDNU(6) = 2 means that X(6) is at unique sorted item(2).
c
c    XU(XDNU(I))) = X(I).
c    XU(I)        = X(UNDX(I)).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of data values.
c
c    Input, integer X_VAL(X_NUM), the data values.
c
c    Input, integer X_UNIQUE_NUM, the number of unique values i
c    n X_VAL.  This value is only required for languages in which the size of
c    UNDX must be known in advance.
c
c    Output, integer UNDX(X_UNIQUE_NUM), the UNDX vector.
c
c    Output, integer XDNU(X_NUM), the XDNU vector.
c
      implicit none

      integer x_num
      integer x_unique_num

      integer i
      integer j
      integer undx(x_unique_num)
      integer x_val(x_num)
      integer xdnu(x_num)
!
!  Walk through the sorted array.
!
      i = 1

      j = 1
      undx(j) = i

      xdnu(i) = j

      do i = 2, x_num

        if ( x_val(i) .ne. x_val(undx(j)) ) then
          j = j + 1
          undx(j) = i
        end if

        xdnu(i) = j

      end do

      return
      end
      subroutine i4vec_sorted_unique ( n, a, unique_num )

c*********************************************************************72
c
cc I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
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
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements in A.
c
c    Input/output, integer A(N).  On input, the sorted
c    integer array.  On output, the unique elements in A.
c
c    Output, integer UNIQUE_NUM, the number of unique elements in A.
c
      implicit none

      integer n

      integer a(n)
      integer itest
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1

      do itest = 2, n

        if ( a(itest) .ne. a(unique_num) ) then
          unique_num = unique_num + 1
          a(unique_num) = a(itest)
        end if

      end do

      return
      end
      subroutine i4vec_sorted_unique_count ( n, a, unique_num )

c*********************************************************************72
c
cc I4VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Because the array is sorted, this algorithm is O(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input, integer A(N), the sorted array to examine.
c
c    Output, integer UNIQUE_NUM, the number of unique elements of A.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer unique_num

      if ( n .lt. 1 ) then
        unique_num = 0
        return
      end if

      unique_num = 1

      do i = 2, n

        if ( a(i-1) .ne. a(i) ) then
          unique_num = unique_num + 1
        end if

      end do

      return
      end
      subroutine i4vec_sorted_unique_hist ( n, a, maxuniq, unique_num,
     &  auniq, acount )

c*********************************************************************72
c
cc I4VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted I4VEC.
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
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input, integer A(N), the array to examine.  The elements of A
c    should have been sorted.
c
c    Input, integer MAXUNIQ, the maximum number of unique elements
c    that can be handled.  If there are more than MAXUNIQ unique
c    elements in A, the excess will be ignored.
c
c    Output, integer UNIQUE_NUM, the number of unique elements.
c
c    Output, integer AUNIQ(UNIQUE_NUM), the unique elements of A.
c
c    Output, integer ACOUNT(UNIQUE_NUM), the number of times
c    each element of AUNIQ occurs in A.
c
      implicit none

      integer maxuniq
      integer n

      integer a(n)
      integer acount(maxuniq)
      integer auniq(maxuniq)
      integer i
      integer unique_num
!
!  Start taking statistics.
!
      unique_num = 0

      do i = 1, n

        if ( i .eq. 1 ) then

          unique_num = 1
          auniq(unique_num) = a(1)
          acount(unique_num) = 1

        else if ( a(i) .eq. auniq(unique_num) ) then

          acount(unique_num) = acount(unique_num) + 1

        else if ( unique_num .lt. maxuniq ) then

          unique_num = unique_num + 1
          auniq(unique_num) = a(i)
          acount(unique_num) = 1

        end if

      end do

      return
      end
      subroutine i4vec_split ( n, a, split, split_index )

c*********************************************************************72
c
cc I4VEC_SPLIT "splits" an unsorted I4VEC based on a splitting value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    If the vector is already sorted, it is simpler to do a binary search
c    on the data than to call this routine.
c
c    The vector is not assumed to be sorted before input, and is not
c    sorted during processing.  If sorting is not needed, then it is
c    more efficient to use this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), the array to split.  On output,
c    all the entries of A that are less than or equal to SPLIT
c    are in A(1:SPLIT_INDEX).
c
c    Input, integer SPLIT, the value used to split the vector.
c    It is not necessary that any value of A actually equal SPLIT.
c
c    Output, integer SPLIT_INDEX, indicates the position of the
c    last entry of the split vector that is less than or equal to SPLIT.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i1
      integer i2
      integer i3
      integer j1
      integer j2
      integer j3
      integer split
      integer split_index
c
c  Partition the vector into A1, A2, A3, where
c    A1 = A(I1:J1) holds values <= SPLIT,
c    A2 = A(I2:J2) holds untested values,
c    A3 = A(I3:J3) holds values > SPLIT.
c
      i1 = 1
      j1 = 0

      i2 = 1
      j2 = n

      i3 = n + 1
      j3 = n
c
c  Pick the next item from A2, and move it into A1 or A3.
c  Adjust indices appropriately.
c
      do i = 1, n

        if ( a(i2) .le. split ) then
          i2 = i2 + 1
          j1 = j1 + 1
        else
          call i4_swap ( a(i2), a(i3-1) )
          i3 = i3 - 1
          j2 = j2 - 1
        end if

      end do

      split_index = j1

      return
      end
      subroutine i4vec_std ( n, a, std )

c*********************************************************************72
c
cc I4VEC_STD returns the standard deviation of an I4VEC.
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
c    Output, double precision STD, the standard deviation of the vector entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      double precision mean
      double precision std

      if ( n .lt. 2 ) then

        std = 0.0D+00

      else

        mean = 0.0D+00
        do i = 1, n
          mean = mean + dble ( a(i) )
        end do
        mean = mean / dble ( n )

        std = 0.0D+00
        do i = 1, n
          std = std + ( dble ( a(i) ) - mean )**2
        end do

        std = sqrt ( std / dble ( n - 1 ) )

      end if

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    SUM function:
c
c      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
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
c    Output, integer I4VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_sum

      i4vec_sum = 0

      do i = 1, n
        i4vec_sum = i4vec_sum + a(i)
      end do

      return
      end
      subroutine i4vec_swap ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_SWAP swaps the entries of two I4VEC's.
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
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the arrays.
c
c    Input/output, integer A1(N), A2(N), the vectors to swap.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer a3
      integer i

      do i = 1, n
        a3    = a1(i)
        a1(i) = a2(i)
        a2(i) = a3
      end do

      return
      end
      subroutine i4vec_transpose_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
c    TITLE = 'My vector:  '
c
c    My vector:
c
c        1    2    3    4    5
c        6    7    8    9   10
c       11
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
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
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ihi
      integer ilo
      character ( len = 11 ) string
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do ilo = 1, n, 5
        ihi = min ( ilo + 5 - 1, n )
        write ( *, '(5i12)' ) ( a(i), i = ilo, ihi)
      end do

      return
      end
      subroutine i4vec_undex ( x_num, x_val, x_unique_num, undx, xdnu )

c*********************************************************************72
c
cc I4VEC_UNDEX returns unique sorted indexes for an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The goal of this routine is to determine a vector UNDX,
c    which points, to the unique elements of X, in sorted order,
c    and a vector XDNU, which identifies, for each entry of X, the index of
c    the unique sorted element of X.
c
c    This is all done with index vectors, so that the elements of
c    X are never moved.
c
c    The first step of the algorithm requires the indexed sorting
c    of X, which creates arrays INDX and XDNI.  (If all the entries
c    of X are unique, then these arrays are the same as UNDX and XDNU.)
c
c    We then use INDX to examine the entries of X in sorted order,
c    noting the unique entries, creating the entries of XDNU and
c    UNDX as we go.
c
c    Once this process has been completed, the vector X could be
c    replaced by a compressed vector XU, containing the unique entries
c    of X in sorted order, using the formula
c
c      XU(1:X_UNIQUE_NUM) = X(UNDX(1:X_UNIQUE_NUM)).
c
c    We could then, if we wished, reconstruct the entire vector X, or
c    any element of it, by index, as follows:
c
c      X(I) = XU(XDNU(I)).
c
c    We could then replace X by the combination of XU and XDNU.
c
c    Later, when we need the I-th entry of X, we can locate it as
c    the XDNU(I)-th entry of XU.
c
c    Here is an example of a vector X, the sort and inverse sort
c    index vectors, and the unique sort and inverse unique sort vectors
c    and the compressed unique sorted vector.
c
c      I    X  Indx  Xdni      XU  Undx  Xdnu
c    ----+----+-----+-----+-------+-----+-----+
c      1 | 11     1     1 |    11     1     1
c      2 | 22     3     5 |    22     2     2
c      3 | 11     6     2 |    33     4     1
c      4 | 33     9     8 |    55     5     3
c      5 | 55     2     9 |                 4
c      6 | 11     7     3 |                 1
c      7 | 22     8     6 |                 2
c      8 | 22     4     7 |                 2
c      9 | 11     5     4 |                 1
c
c    INDX(2) = 3 means that sorted item(2) is X(3).
c    XDNI(2) = 5 means that X(2) is sorted item(5).
c
c    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
c    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
c
c    XU(XDNU(I))) = X(I).
c    XU(I)        = X(UNDX(I)).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of data values.
c
c    Input, integer X_VAL(X_NUM), the data values.
c
c    Input, integer X_UNIQUE_NUM, the number of unique values
c    in X_VAL.  This value is only required for languages in which the size of
c    UNDX must be known in advance.
c
c    Output, integer UNDX(X_UNIQUE_NUM), the UNDX vector.
c
c    Output, integer XDNU(X_NUM), the XDNU vector.
c
      implicit none

      integer x_num
      integer x_unique_num

      integer i
      integer indx(x_num)
      integer j
      integer undx(x_unique_num)
      integer x_val(x_num)
      integer xdnu(x_num)
c
c  Implicitly sort the array.
c
      call i4vec_sort_heap_index_a ( x_num, x_val, indx )
c
c  Walk through the implicitly sorted array.
c
      i = 1

      j = 1
      undx(j) = indx(i)

      xdnu(indx(i)) = j

      do i = 2, x_num

        if ( x_val(indx(i)) .ne. x_val(undx(j)) ) then
          j = j + 1
          undx(j) = indx(i)
        end if

        xdnu(indx(i)) = j

      end do

      return
      end
      subroutine i4vec_uniform_ab ( n, a, b, seed, x )

c*********************************************************************72
c
cc I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The pseudorandom numbers should be uniformly distributed
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
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer X(N), a vector of numbers between A and B.
c
      implicit none

      integer n

      integer a
      integer b
      integer i
      integer k
      real r
      integer seed
      integer value
      integer x(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

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
     &    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
        value = nint ( r )

        value = max ( value, min ( a, b ) )
        value = min ( value, max ( a, b ) )

        x(i) = value

      end do

      return
      end
      subroutine i4vec_unique_count ( n, a, unique_num )

c*********************************************************************72
c
cc I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Because the array is unsorted, this algorithm is O(N^2).
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
c    Input, integer N, the number of elements of A.
c
c    Input, integer A(N), the unsorted array to examine.
c
c    Output, integer UNIQUE_NUM, the number of unique elements
c    of A.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer unique_num

      unique_num = 0

      do i = 1, n

        unique_num = unique_num + 1

        do j = 1, i - 1

          if ( a(i) .eq. a(j) ) then
            unique_num = unique_num - 1
            go to 10
          end if

        end do

10      continue

      end do

      return
      end
      subroutine i4vec_unique_index ( n, a, unique_index )

c*********************************************************************72
c
cc I4VEC_UNIQUE_INDEX indexes the first occurrence of values in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    For element A(I) of the vector, FIRST_UNIQUE(I) is the uniqueness index
c    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
c    gathered in order, then
c
c      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
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
c    Input, integer N, the number of elements of A.
c
c    Input, integer A(N), the array.
c
c    Output, integer UNIQUE_INDEX(N), the unique index.
c
      implicit none

      integer  n

      integer a(n)
      integer i
      integer j
      integer unique_index(n)
      integer unique_num

      do i = 1, n
        unique_index(i) = -1
      end do
      unique_num = 0

      do i = 1, n

        if ( unique_index(i) .eq. -1 ) then

          unique_num = unique_num + 1
          unique_index(i) = unique_num

          do j = i + 1, n
            if ( a(i) .eq. a(j) ) then
              unique_index(j) = unique_num
            end if
          end do

        end if

      end do

      return
      end
      subroutine i4vec_value_index ( n, a, value, max_index, n_index,
     &  value_index )

c*********************************************************************72
c
cc I4VEC_VALUE_INDEX indexes entries equal to a given value in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Example:
c
c    Input:
c
c      N = 10
c      A = (  2, 3, 1, 3, 2, 4, 2, 3, 5, 3 )
c      X_VALUE = 3
c
c    Output:
c
c      N_INDEX = 4
c      VALUE_INDEX = ( 2, 4, 8, 10 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer A(N), the array to be indexed.
c
c    Input, integer VALUE, a value to be searched for.
c
c    Input, integer MAX_INDEX, the maximum number of indices
c    to find.
c
c    Output, integer N_INDEX, the number of entries equal to VALUE.
c
c    Output, integer VALUE_INDEX(MAX_INDEX), the indices of entries
c    equal to VALUE.
c
      implicit none

      integer max_index
      integer n

      integer a(n)
      integer i
      integer n_index
      integer value
      integer value_index(max_index)

      n_index = 0

      do i = 1, n

        if ( a(i) .eq. value ) then

          if ( max_index .le. n_index ) then
            return
          end if

          n_index = n_index + 1
          value_index(n_index) = i

        end if

      end do

      return
      end
      subroutine i4vec_value_num ( n, a, value, value_num )

c*********************************************************************72
c
cc I4VEC_VALUE_NUM counts entries equal to a given value in an I4VEC.
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
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer A(N), the array to be indexed.
c
c    Input, integer VALUE, a value to be searched for.
c
c    Input, integer VALUE_NUM, the number of times the
c    value occurs.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer value
      integer value_num

      value_num = 0

      do i = 1, n

        if ( a(i) .eq. value ) then
          value_num = value_num + 1
        end if

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
      function i4vec_width ( n, a )

c*********************************************************************72
c
cc I4VEC_WIDTH returns the "width" of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The width of an integer vector is simply the maximum of the widths of
c    its entries.
c
c    The width of a single integer is the number of characters
c    necessary to print it.
c
c    The width of an integer vector can be useful when the vector is
c    to be printed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector.
c
c    Output, integer I4VEC_WIDTH, the width of the vector.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_width
      integer i4vec_width

      i4vec_width = -1

      do i = 1, n
        i4vec_width = max ( i4vec_width, i4_width ( a(i) ) )
      end do

      return
      end
      subroutine i4vec_zero ( n, a )

c*********************************************************************72
c
cc I4VEC_ZERO sets the entries of an I4VEC to 0.
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
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Output, integer A(N), the vector, which has been set to zero.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = 0
      end do

      return
      end
      subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

c*********************************************************************72
c
cc I4VEC2_COMPARE compares pairs of integers stored in two vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data items.
c
c    Input, integer A1(N), A2(N), contain the two components
c    of each item.
c
c    Input, integer I, J, the items to be compared.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, item I .lt. item J,
c     0, item I = item J,
c    +1, item J .lt. item I.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer isgn
      integer j

      isgn = 0

           if ( a1(i) .lt. a1(j) ) then

        isgn = -1

      else if ( a1(i) .eq. a1(j) ) then

             if ( a2(i) .lt. a2(j) ) then
          isgn = -1
        else if ( a2(i) .lt. a2(j) ) then
          isgn = 0
        else if ( a2(j) .lt. a2(i) ) then
          isgn = +1
        end if

      else if ( a1(j) .lt. a1(i) ) then

        isgn = +1

      end if

      return
      end
      subroutine i4vec2_print ( n, a, b, title )

c*********************************************************************72
c
cc I4VEC2_PRINT prints a pair of integer vectors.
c
c  Discussion:
c
c    An I4VEC2 is a pair of I4VEC's.
c
c    An I4VEC is a vector of I4's.
c
c    Entry K of an I4VEC2 is the pair of values located
c    at the K-th entries of the two I4VEC's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), B(N), the vectors to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer i
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i10,2x,i10)' ) i, ':', a(i), b(i)
      end do

      return
      end
      subroutine i4vec2_sort_a ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
c
c  Discussion:
c
c    Each item to be sorted is a pair of integers (I,J), with the I
c    and J values stored in separate vectors A1 and A2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items of data.
c
c    Input/output, integer A1(N), A2(N), the data to be sorted.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer indx
      integer isgn
      integer j
      integer temp

      if ( n .le. 1 ) then
        return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          temp  = a1(i)
          a1(i) = a1(j)
          a1(j) = temp

          temp  = a2(i)
          a2(i) = a2(j)
          a2(j) = temp
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4vec2_compare ( n, a1, a2, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec2_sort_d ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC2_SORT_D descending sorts a vector of pairs of integers.
c
c  Discussion:
c
c    An I4VEC2 is a pair of I4VEC's.
c
c    An I4VEC is a vector of I4's.
c
c    Entry K of an I4VEC2 is the pair of values located
c    at the K-th entries of the two I4VEC's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items of data.
c
c    Input/output, integer A1(N), A2(N), the data to be sorted.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer indx
      integer isgn
      integer j
      integer t

      if ( n .le. 1 ) then
       return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          t     = a1(i)
          a1(i) = a1(j)
          a1(j) = t

          t     = a2(i)
          a2(i) = a2(j)
          a2(j) = t
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4vec2_compare ( n, a1, a2, i, j, isgn )
          isgn = -isgn

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

c*********************************************************************72
c
cc I4VEC2_SORTED_UNIQUE gets the unique elements in a sorted I4VEC2.
c
c  Discussion:
c
c    Item I is stored as the pair A1(I), A2(I).
c
c    The items must have been sorted, or at least it must be the
c    case that equal items are stored in adjacent vector locations.
c
c    If the items were not sorted, then this routine will only
c    replace a string of equal values by a single representative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items.
c
c    Input/output, integer A1(N), A2(N).
c    On input, the array of N items.
c    On output, an array of unique items.
c
c    Output, integer UNIQUE_NUM, the number of unique items.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer itest
      integer unique_num

      unique_num = 0

      if ( n .le. 0 ) then
        return
      end if

      unique_num = 1

      do itest = 2, n

        if ( a1(itest) .ne. a1(unique_num) .or.
     &       a2(itest) .ne. a2(unique_num) ) then

          unique_num = unique_num + 1

          a1(unique_num) = a1(itest)
          a2(unique_num) = a2(itest)

        end if

      end do

      return
      end
      function l_to_i4 ( l )

c*********************************************************************72
c
cc L_TO_I4 converts a logical value to an I4.
c
c  Discussion:
c
c    0 is FALSE, and anything else if TRUE.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, logical L, a logical value.
c
c    Output, integer L_TO_I4, the integer value of L.
c
      implicit none

      logical l
      integer l_to_i4
      integer value

      if ( l ) then
        value = 1
      else
        value = 0
      end if

      l_to_i4 = value

      return
      end
      subroutine perm_check ( n, p, base, ierror )

c*********************************************************************72
c
cc PERM_CHECK checks that a vector represents a permutation.
c
c  Discussion:
c
c    The routine verifies that each of the integers from BASE to
c    to BASE+N-1 occurs among the N entries of the permutation.
c
c    Set the input quantity BASE to 0, if P is a 0-based permutation,
c    or to 1 if P is a 1-based permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the array to check.
c
c    Input, integer BASE, the index base.
c
c    Output, integer IERROR, error flag.
c    0, the array represents a permutation.
c    nonzero, the array does not represent a permutation.  The smallest
c    missing value is equal to IERROR.
c
      implicit none

      integer n

      integer base
      integer find
      integer ierror
      integer p(n)
      integer seek

      ierror = 0

      do seek = base, base + n - 1

        ierror = 1

        do find = 1, n
          if ( p(find) .eq. seek ) then
            ierror = 0
            go to 10
          end if
        end do

10      continue

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a proper permutation.'
          stop
        end if

      end do

      return
      end
      subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

c*********************************************************************72
c
cc PERM_CYCLE analyzes a permutation.
c
c  Discussion:
c
c    The routine will count cycles, find the sign of a permutation,
c    and tag a permutation.
c
c  Example:
c
c    Input:
c
c      N = 9
c      IOPT = 1
c      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      NCYCLE = 3
c      ISGN = +1
c      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer IOPT, requests tagging.
c    0, the permutation will not be tagged.
c    1, the permutation will be tagged.
c
c    Input/output, integer P(N).  On input, P describes a
c    permutation, in the sense that entry I is to be moved to P(I).
c    If IOPT = 0, then P will not be changed by this routine.
c    If IOPT = 1, then on output, P will be "tagged".  That is,
c    one element of every cycle in P will be negated.  In this way,
c    a user can traverse a cycle by starting at any entry I1 of P
c    which is negative, moving to I2 = ABS(P(I1)), then to
c    P(I2), and so on, until returning to I1.
c
c    Output, integer ISGN, the "sign" of the permutation, which is
c    +1 if the permutation is even, -1 if odd.  Every permutation
c    may be produced by a certain number of pairwise switches.
c    If the number of switches is even, the permutation itself is
c    called even.
c
c    Output, integer NCYCLE, the number of cycles in the permutation.
c
      implicit none

      integer n

      integer base
      integer i
      integer i1
      integer i2
      integer ierror
      integer iopt
      integer is
      integer isgn
      integer ncycle
      integer p(n)

      base = 1
      call perm_check ( n, p, base, ierror )

      is = 1
      ncycle = n

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          ncycle = ncycle - 1
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        if ( iopt .ne. 0 ) then
          is = -sign ( 1, p(i) )
        end if

        p(i) = sign ( p(i), is )

      end do

      isgn = 1 - 2 * mod ( n - ncycle, 2 )

      return
      end
      subroutine perm_uniform ( n, seed, p )

c*********************************************************************72
c
cc PERM_UNIFORM selects a random permutation of N objects.
c
c  Discussion:
c
c    The routine assumes the objects are labeled 1, 2, ... N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer P(N), a permutation of ( 1, 2, ..., N ), in standard
c    index form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform_ab
      integer j
      integer p(n)
      integer pk
      integer seed

      do i = 1, n
        p(i) = i
      end do

      do i = 1, n
        j = i4_uniform_ab ( i, n, seed )
        pk = p(i)
        p(i) = p(j)
        p(j) = pk
      end do

      return
      end
      function prime ( n )

c*********************************************************************72
c
cc PRIME returns any of the first PRIME_MAX prime numbers.
c
c  Discussion:
c
c    PRIME_MAX is 1600, and the largest prime stored is 13499.
c
c    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Daniel Zwillinger,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996, pages 95-98.
c
c  Parameters:
c
c    Input, integer N, the index of the desired prime number.
c    In general, is should be true that 0 <= N <= PRIME_MAX.
c    N = -1 returns PRIME_MAX, the index of the largest prime available.
c    N = 0 is legal, returning PRIME = 1.
c
c    Output, integer PRIME, the N-th prime.  If N is out of range,
c    PRIME is returned as -1.
c
      implicit none

      integer prime_max
      parameter ( prime_max = 1600 )

      integer i
      integer n
      integer npvec(prime_max)
      integer prime

      save npvec

      data ( npvec(i), i = 1, 100 ) /
     &      2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     &     31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     &     73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
     &    127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
     &    179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
     &    233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
     &    283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
     &    353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
     &    419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
     &    467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /

      data ( npvec(i), i = 101, 200 ) /
     &    547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
     &    607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
     &    661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
     &    739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
     &    811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
     &    877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
     &    947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     &   1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     &   1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     &   1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /

      data ( npvec(i), i = 201, 300 ) /
     &   1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
     &   1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
     &   1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
     &   1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
     &   1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
     &   1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
     &   1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
     &   1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
     &   1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     &   1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /

      data ( npvec(i), i = 301, 400 ) /
     &   1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     &   2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     &   2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     &   2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     &   2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
     &   2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
     &   2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
     &   2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
     &   2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
     &   2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /

      data ( npvec(i), i = 401, 500 ) /
     &   2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
     &   2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
     &   2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
     &   3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
     &   3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
     &   3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
     &   3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
     &   3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
     &   3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
     &   3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /

      data ( npvec(i), i = 501, 600 ) /
     &   3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     &   3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     &   3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     &   3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     &   3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     &   4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     &   4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     &   4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     &   4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     &   4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /

      data ( npvec(i), i = 601, 700 ) /
     &   4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
     &   4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
     &   4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
     &   4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
     &   4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
     &   4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
     &   4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
     &   5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
     &   5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
     &   5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /

      data ( npvec(i), i = 701, 800 ) /
     &   5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     &   5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     &   5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     &   5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     &   5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     &   5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     &   5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     &   5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     &   5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     &   6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /

      data ( npvec(i), i = 801, 900 ) /
     &   6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     &   6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     &   6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
     &   6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
     &   6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
     &   6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
     &   6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
     &   6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
     &   6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
     &   6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /

      data ( npvec(i), i = 901, 1000 ) /
     &   7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
     &   7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
     &   7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
     &   7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
     &   7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
     &   7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
     &   7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
     &   7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
     &   7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
     &   7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /

      data ( npvec(i), i = 1001, 1100 ) /
     &   7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     &   8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     &   8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     &   8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     &   8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     &   8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     &   8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     &   8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     &   8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
     &   8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /

      data ( npvec(i), i = 1101, 1200 ) /
     &   8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     &   8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     &   9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     &   9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     &   9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     &   9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     &   9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     &   9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     &   9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     &   9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /

      data ( npvec(i), i = 1201, 1300 ) /
     &   9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     &   9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     &   9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
     &  10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
     &  10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
     &  10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
     &  10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
     &  10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
     &  10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
     &  10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /

      data ( npvec(i), i = 1301, 1400 ) /
     &  10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
     &  10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
     &  10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
     &  10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
     &  11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
     &  11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
     &  11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
     &  11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
     &  11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
     &  11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /

      data ( npvec(i), i = 1401, 1500 ) /
     &  11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
     &  11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
     &  11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
     &  11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
     &  12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
     &  12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
     &  12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
     &  12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
     &  12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
     &  12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /

      data ( npvec(i), i = 1501, 1600 ) /
     &  12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
     &  12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
     &  12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
     &  12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
     &  12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
     &  13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
     &  13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
     &  13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
     &  13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
     &  13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /

      if ( n .eq. -1 ) then
        prime = prime_max
      else if ( n .eq. 0 ) then
        prime = 1
      else if ( n .le. prime_max ) then
        prime = npvec(n)
      else
        prime = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRIME - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
        write ( *, '(a,i8)' )
     &    '  N should be between 1 and PRIME_MAX =', prime_max
        stop
      end if

      return
      end
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
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
