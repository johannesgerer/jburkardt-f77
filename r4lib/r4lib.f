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
      subroutine i4int_to_r4int ( imin, imax, i, rmin, rmax, r )

c*********************************************************************72
c
cc I4INT_TO_R4INT maps an I4INT to an R4INT.
c
c  Discussion:
c
c    The formula used is:
c
c      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 September 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IMIN, IMAX, the range.
c
c    Input, integer I, the integer to be converted.
c
c    Input, real RMIN, RMAX, the range.
c
c    Output, real R, the corresponding value in [RMIN,RMAX].
c
      implicit none

      integer i
      integer imax
      integer imin
      real r
      real rmax
      real rmin

      if ( imax .eq. imin ) then

        r = 0.5E+00 * ( rmin + rmax )

      else

        r = ( real ( imax - i        ) * rmin   
     &      + real (        i - imin ) * rmax ) 
     &      / real ( imax     - imin )

      end if

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
      subroutine legendre_zeros ( n, x )

c*********************************************************************72
c
cc LEGENDRE_ZEROS computes the zeros of the Legendre polynomial of degree N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2011
c
c  Author:
c
c    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the order.
c    0 .lt. N.
c
c    Output, real X(N), the abscissas.
c
      implicit none

      integer n

      real d1
      real d2pn
      real d3pn
      real d4pn
      real dp
      real dpn
      real e1
      real fx
      real h
      integer i
      integer iback
      integer k
      integer m
      integer mp1mi
      integer ncopy
      integer nmove
      real p
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real pk
      real pkm1
      real pkp1
      real t
      real u
      real v
      real x(n)
      real x0
      real xtemp

      e1 = real ( n * ( n + 1 ) )

      m = ( n + 1 ) / 2

      do i = 1, m

        mp1mi = m + 1 - i

        t = real ( 4 * i - 1 ) * pi 
     &    / real ( 4 * n + 2 )

        x0 = cos ( t ) * ( 1.0E+00 - ( 1.0E+00 - 1.0E+00 
     &    / real ( n ) ) 
     &    / real ( 8 * n * n ) )

        pkm1 = 1.0E+00
        pk = x0

        do k = 2, n
          pkp1 = 2.0E+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) 
     &      / real ( k )
          pkm1 = pk
          pk = pkp1
        end do

        d1 = real ( n ) * ( pkm1 - x0 * pk )

        dpn = d1 / ( 1.0E+00 - x0 ) / ( 1.0E+00 + x0 )

        d2pn = ( 2.0E+00 * x0 * dpn - e1 * pk ) / ( 1.0E+00 - x0 ) 
     &    / ( 1.0E+00 + x0 )

        d3pn = ( 4.0E+00 * x0 * d2pn + ( 2.0E+00 - e1 ) * dpn ) 
     &    / ( 1.0E+00 - x0 ) / ( 1.0E+00 + x0 )

        d4pn = ( 6.0E+00 * x0 * d3pn + ( 6.0E+00 - e1 ) * d2pn ) 
     &    / ( 1.0E+00 - x0 ) / ( 1.0E+00 + x0 )

        u = pk / dpn
        v = d2pn / dpn
c
c  Initial approximation H:
c
        h = - u * ( 1.0E+00 + 0.5E+00 * u * ( v + u * ( v * v - d3pn / 
     &    ( 3.0E+00 * dpn ) ) ) )
c
c  Refine H using one step of Newton's method:
c
        p = pk + h * ( dpn + 0.5E+00 * h * ( d2pn + h / 3.0E+00 
     &    * ( d3pn + 0.25E+00 * h * d4pn ) ) )

        dp = dpn + h * ( d2pn + 0.5E+00 * h * 
     &    ( d3pn + h * d4pn / 3.0E+00 ) )

        h = h - p / dp

        xtemp = x0 + h

        x(mp1mi) = xtemp

        fx = d1 - h * e1 * ( pk + 0.5E+00 * h * ( dpn + h / 3.0E+00 
     &    * ( d2pn + 0.25E+00 * h 
     &    * ( d3pn + 0.2E+00 * h * d4pn ) ) ) )

      end do

      if ( mod ( n, 2 ) .eq. 1 ) then
        x(1) = 0.0E+00
      end if
c
c  Shift the data up.
c
      nmove = ( n + 1 ) / 2
      ncopy = n - nmove

      do i = 1, nmove
        iback = n + 1 - i
        x(iback) = x(iback-ncopy)
      end do
c
c  Reflect values for the negative abscissas.
c
      do i = 1, n - nmove
        x(i) = - x(n+1-i)
      end do

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
c    Second Edition,
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
      integer i4_uniform
      integer j
      integer p(n)
      integer pk
      integer seed

      do i = 1, n
        p(i) = i
      end do

      do i = 1, n
        j = i4_uniform ( i, n, seed )
        pk = p(i)
        p(i) = p(j)
        p(j) = pk
      end do

      return
      end
      function r4_abs ( x )

c*********************************************************************72
c
cc R4_ABS returns the absolute value of an R4.
c
c  Discussion:
c
c    FORTRAN90 supplies the ABS function, which should be used instead
c    of this function!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose absolute value is desired.
c
c    Output, real R4_ABS, the absolute value of X.
c
      implicit none

      real r4_abs
      real x

      if ( 0.0 .le. x ) then
        r4_abs = + x
      else
        r4_abs = - x
      end if

      return
      end
      function r4_add ( x, y )

c*********************************************************************72
c
cc R4_ADD adds two R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, Y, the numbers to be added.
c
c    Output, real R4_ADD, the sum of X and Y.
c
      implicit none

      real r4_add
      real x
      real y

      r4_add = x + y

      return
      end
      function r4_aint ( x )

c********************************************************************72
c
cc R4_AINT truncates an R4 argument to an integer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2011
c
c  Author:
c
c    John Burkardt.
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real VALUE, the truncated version of X.
c
      implicit none

      real r4_aint
      real value
      real x

      if ( x .lt. 0.0E+00 ) then
        value = - int ( abs ( x ) )
      else
        value =   int ( abs ( x ) )
      end if

      r4_aint = value

      return
      end
      function r4_atan ( y, x )

c*********************************************************************72
c
cc R4_ATAN computes the inverse tangent of the ratio Y / X.
c
c  Discussion:
c
c    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
c    the built in functions ATAN and ATAN2 already do.
c
c    However:
c
c    * R4_ATAN always returns a positive angle, between 0 and 2 PI,
c      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
c      and [-PI,+PI] respectively;
c
c    * R4_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
c     function by contrast always returns an angle in the first or fourth
c     quadrants.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real Y, X, two quantities which represent the
c    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
c
c    Output, real R4_ATAN, an angle between 0 and 2 * PI, whose
c    tangent is (Y/X), and which lies in the appropriate quadrant so that
c    the signs of its cosine and sine match those of X and Y.
c
      implicit none

      real abs_x
      real abs_y
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r4_atan
      real theta
      real theta_0
      real x
      real y
c
c  Special cases:
c
      if ( x .eq. 0.0E+00 ) then

        if ( 0.0E+00 .lt. y ) then
          theta = pi / 2.0E+00
        else if ( y .lt. 0.0E+00 ) then
          theta = 3.0E+00 * pi / 2.0E+00
        else if ( y .eq. 0.0E+00 ) then
          theta = 0.0E+00
        end if

      else if ( y .eq. 0.0E+00 ) then

        if ( 0.0E+00 .lt. x ) then
          theta = 0.0E+00
        else if ( x .lt. 0.0E+00 ) then
          theta = pi
        end if
c
c  We assume that ATAN2 is correct when both arguments are positive.
c
      else

        abs_y = abs ( y )
        abs_x = abs ( x )

        theta_0 = atan2 ( abs_y, abs_x )

        if ( 0.0E+00 .lt. x .and. 0.0E+00 .lt. y ) then
          theta = theta_0
        else if ( x .lt. 0.0E+00 .and. 0.0E+00 .lt. y ) then
          theta = pi - theta_0
        else if ( x .lt. 0.0E+00 .and. y .lt. 0.0E+00 ) then
          theta = pi + theta_0
        else if ( 0.0E+00 .lt. x .and. y .lt. 0.0E+00 ) then
          theta = 2.0E+00 * pi - theta_0
        end if

      end if

      r4_atan = theta

      return
      end
      function r4_cas ( x )

c*********************************************************************72
c
cc R4_CAS returns the "casine" of an R4 value.
c
c  Discussion:
c
c    The "casine", used in the discrete Hartley transform, is abbreviated
c    CAS(X), and defined by:
c
c      CAS(X) = cos ( X ) + sin( X )
c             = sqrt ( 2 ) * sin ( X + pi/4 )
c             = sqrt ( 2 ) * cos ( X - pi/4 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ralph Hartley,
c    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
c    Proceedings of the Institute of Radio Engineers,
c    Volume 30, pages 144-150, 1942.
c
c  Parameters:
c
c    Input, real X, the number whose casine is desired.
c
c    Output, real R4_CAS, the casine of X, which will be between
c    plus or minus the square root of 2.
c
      implicit none

      real r4_cas
      real x

      r4_cas = cos ( x ) + sin ( x )

      return
      end
      function r4_ceiling ( r )

c*********************************************************************72
c
cc R4_CEILING rounds an R4 "up" to the nearest integer.
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
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R, the value to be rounded up.
c
c    Output, integer R4_CEILING, the rounded value.
c
      implicit none

      real r
      integer r4_ceiling
      integer value

      value = int ( r )
      if ( real ( value ) .lt. r ) then
        value = value + 1
      end if

      r4_ceiling = value

      return
      end
      function r4_choose ( n, k )

c*********************************************************************72
c
cc R4_CHOOSE computes the binomial coefficient C(N,K) as an R4.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R4 arithmetic.
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
c    11 August 2010
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
c    Output, real R4_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      real r4_choose
      real value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0E+00

      else if ( mn .eq. 0 ) then

        value = 1.0E+00

      else

        mx = max ( k, n - k )
        value = real ( mx + 1 )

        do i = 2, mn
          value = ( value * real ( mx + i ) ) / real ( i )
        end do

      end if

      r4_choose = value

      return
      end
      function r4_chop ( place, x )

c*********************************************************************72
c
cc R4_CHOP chops an R4 to a given number of binary places.
c
c  Example:
c
c    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
c
c    The following values would be returned for the 'chopped' value of
c    3.875:
c
c    PLACE  Value
c
c       1      2
c       2      3     = 2 + 1
c       3      3.5   = 2 + 1 + 1/2
c       4      3.75  = 2 + 1 + 1/2 + 1/4
c       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PLACE, the number of binary places to preserve.
c    PLACE = 0 means return the integer part of X.
c    PLACE = 1 means return the value of X, correct to 1/2.
c    PLACE = 2 means return the value of X, correct to 1/4.
c    PLACE = -1 means return the value of X, correct to 2.
c
c    Input, real X, the number to be chopped.
c
c    Output, real R4_CHOP, the chopped number.
c
      implicit none

      real fac
      integer place
      real r4_chop
      real r4_log_2
      real r4_sign
      real s
      integer temp
      real x

      s = r4_sign ( x )
      temp = int ( r4_log_2 ( abs ( x ) ) )
      fac = 2.0E+00**( temp - place + 1 )
      r4_chop = s * real ( int ( abs ( x ) / fac ) ) * fac

      return
      end
      function r4_csqrt ( x )

c*********************************************************************72
c
cc R4_CSQRT returns the complex square root of an R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose square root is desired.
c
c    Output, double complex R4_CSQRT, the square root of X:
c
      implicit none

      real argument
      real magnitude
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      double complex r4_csqrt
      real x

      if ( 0.0E+00 .lt. x ) then
        magnitude = x
        argument = 0.0E+00
      else if ( 0.0E+00 .eq. x ) then
        magnitude = 0.0E+00
        argument = 0.0E+00
      else if ( x .lt. 0.0E+00 ) then
        magnitude = -x
        argument = pi
      end if

      magnitude = sqrt ( magnitude )
      argument = argument / 2.0E+00

      r4_csqrt = magnitude
     &  * dcmplx ( cos ( argument ), sin ( argument ) )

      return
      end
      function r4_cube_root ( x )

c*********************************************************************72
c
cc R4_CUBE_ROOT returns the cube root of an R4.
c
c  Discussion:
c
c    This routine is designed to avoid the possible problems that can occur
c    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose cube root is desired.
c
c    Output, real R4_CUBE_ROOT, the cube root of X.
c
      implicit none

      real r4_cube_root
      real value
      real x

      if ( 0.0E+00 .lt. x ) then
        value = x**(1.0E+00/3.0E+00)
      else if ( x .eq. 0.0E+00 ) then
        value = 0.0E+00
      else
        value = - ( abs ( x ) )**(1.0E+00/3.0E+00)
      end if

      r4_cube_root = value

      return
      end
      function r4_diff ( x, y, n )

c*********************************************************************72
c
cc R4_DIFF computes the difference of two R4's to a specified accuracy.
c
c  Discussion:
c
c    The user controls how many binary digits of accuracy
c    are to be used.
c
c    N determines the accuracy of the value of the result.  If N = 10,
c    for example, only 11 binary places will be used in the arithmetic.
c    In general, only N+1 binary places will be used.
c
c    N may be zero.  However, a negative value of N should
c    not be used, since this will cause both X and Y to look like 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, Y, the two values whose difference is desired.
c
c    Input, integer N, the number of binary digits to use.
c
c    Output, real R4_DIFF, the value of X-Y.
c
      implicit none

      real cx
      real cy
      integer n
      real pow2
      real r4_diff
      real size
      real x
      real y

      if ( x .eq. y ) then
        r4_diff = 0.0E+00
        return
      end if

      pow2 = 2.0E+00**n
c
c  Compute the magnitude of X and Y, and take the larger of the
c  two.  At least one of the two values is not zero!
c
      size = max ( abs ( x ), abs ( y ) )
c
c  Make normalized copies of X and Y.  One of the two values will
c  actually be equal to 1.
c
      cx = x / size
      cy = y / size
c
c  Here's where rounding comes in.  We know that the larger of the
c  the two values equals 1.  We multiply both values by 2**N,
c  where N+1 is the number of binary digits of accuracy we want
c  to use, truncate the values, and divide back by 2**N.
c
      cx = real ( int ( cx * pow2 + sign ( 0.5E+00, cx ) ) ) / pow2
      cy = real ( int ( cy * pow2 + sign ( 0.5E+00, cy ) ) ) / pow2
c
c  Take the difference now.
c
      r4_diff = cx - cy
c
c  Undo the scaling.
c
      r4_diff = r4_diff * size

      return
      end
      subroutine r4_digit ( x, idigit, digit )

c*********************************************************************72
c
cc R4_DIGIT returns a particular decimal digit of an R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real real X, the number whose NDIG-th decimal digit
c    is desired.  If X is zero, all digits will be returned as 0.
c
c    Input, integer IDIGIT, the position of the desired decimal
c    digit.  A value of 1 means the leading digit, a value of 2 the second digit
c    and so on.
c
c    Output, integer DIGIT, the value of the IDIGIT-th decimal
c    digit of X.
c
      implicit none

      integer digit
      integer i
      integer idigit
      integer ival
      real x
      real xcopy

      if ( x .eq. 0.0E+00 ) then
        digit = 0
        return
      end if

      if ( idigit .le. 0 ) then
        digit = 0
        return
      end if
c
c  Set XCOPY = X, and then force XCOPY to lie between 1 and 10.
c
      xcopy = abs ( x )

10    continue

      if ( xcopy .lt. 1.0E+00 ) then
        xcopy = xcopy * 10.0E+00
        go to 10
      end if

20    continue

      if ( 10.0E+00 .le. xcopy ) then
        xcopy = xcopy / 10.0E+00
        go to 20
      end if

      do i = 1, idigit
        ival = int ( xcopy )
        xcopy = ( xcopy - ival ) * 10.0E+00
      end do

      digit = ival

      return
      end
      function r4_epsilon ( )

c*********************************************************************72
c
cc R4_EPSILON returns the R4 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 December 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_EPSILON, the R4 roundoff unit.
c
      implicit none

      real r4_epsilon

      r4_epsilon = 1.19209290E-07

      return
      end
      function r4_epsilon_compute ( )

c*********************************************************************72
c
cc R4_EPSILON_COMPUTE returns the R4 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the property
c    that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_EPSILON_COMPUTE, the R4 roundoff unit.
c
      implicit none

      real one
      real r4_add
      real r4_epsilon_compute
      real temp
      real test
      real value

      one = real ( 1 )

      value = one
      temp = value / 2.0E+00
      test = r4_add ( one, temp )

10    continue

      if ( one .lt. test ) then
        value = temp
        temp = value / 2.0E+00
        test = r4_add ( one, temp )
        go to 10
      end if

      r4_epsilon_compute = value

      return
      end
      function r4_exp ( x )

c*********************************************************************72
c
cc R4_EXP computes the exponential function, avoiding overflow and underflow.
c
c  Discussion:
c
c    My experience with the G95 compiler has included many unpleasant
c    floating point exceptions when very small arguments are given to
c    the exponential function.
c
c    This routine is designed to avoid such problems.
c
c    Ideally, the rule would be:
c
c                    X <= log ( TINY ) => R4_EXP ( X ) = 0
c    log ( HUGE ) <= X                 => R4_EXP ( X ) = HUGE
c
c    However, the G95 math library seems to produce infinity for
c    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
c    included a fudge factor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument of the exponential function.
c
c    Output, real R4_EXP, the value of exp ( X ).
c
      implicit none

      real log_max
      parameter ( log_max = 88.71397E+00 )
      real log_min
      parameter ( log_min = -87.34528E+00 )
      real r4_exp
      real r4_huge
      real x

      if ( x .le. log_min ) then
        r4_exp = 0.0E+00
      else if ( x .lt. log_max ) then
        r4_exp = exp ( x )
      else
        r4_exp = r4_huge ( )
      end if

      return
      end
      function r4_factorial ( n )

c*********************************************************************72
c
cc R4_FACTORIAL computes the factorial of N.
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
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, real R4_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      real r4_factorial

      r4_factorial = 1.0E+00

      do i = 1, n
        r4_factorial = r4_factorial * real ( i )
      end do

      return
      end
      function r4_factorial2 ( n )

c*********************************************************************72
c
cc R4_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N   Value
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
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
c    Input, integer N, the argument of the double factorial
c    function.  If N is less than 1, R4_FACTORIAL2 is returned as 1.0.
c
c    Output, real R4_FACTORIAL2, the value.
c
      implicit none

      integer n
      real r4_factorial2
      real r4_n

      if ( n .lt. 1 ) then
        r4_factorial2 = 1.0E+00
        return
      end if

      r4_n = real ( n )
      r4_factorial2 = 1.0E+00

10    continue

      if ( 1.0E+00 .lt. r4_n ) then
        r4_factorial2 = r4_factorial2 * r4_n
        r4_n = r4_n - 2.0E+00
        go to 10
      end if

      return
      end
      function r4_floor ( r )

c*********************************************************************72
c
cc R4_FLOOR rounds an R4 "down" (towards -infinity) to the next integer.
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
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R, the value to be rounded down.
c
c    Output, integer R4_FLOOR, the rounded value.
c
      implicit none

      real r
      integer r4_floor
      integer value

      value = int ( r )
      if ( r .lt. real ( value ) ) then
        value = value - 1
      end if

      r4_floor = value

      return
      end
      function r4_fraction ( i, j )

c*********************************************************************72
c
cc R4_FRACTION uses real arithmetic on an integer ratio.
c
c  Discussion:
c
c    Given integer variables I and J, both FORTRAN and C will evaluate
c    an expression such as "I/J" using what is called "integer division",
c    with the result being an integer.  It is often convenient to express
c    the parts of a fraction as integers but expect the result to be computed
c    using real arithmetic.  This function carries out that operation.
c
c  Example:
c
c       I     J   I/J  R4_FRACTION
c
c       1     2     0  0.5
c       7     4     1  1.75
c       8     4     2  2.00
c       9     4     2  2.25
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the arguments.
c
c    Output, real R4_FRACTION, the value of the ratio.
c
      implicit none

      integer i
      integer j
      real r4_fraction

      r4_fraction = real ( i ) / real ( j )

      return
      end
      function r4_fractional ( x )

c*********************************************************************72
c
cc R4_FRACTIONAL returns the fraction part of an R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real R4_FRACTIONAL, the fraction part of X.
c
      implicit none

      real r4_fractional
      real x

      r4_fractional = abs ( x ) - real ( int ( abs ( x ) ) )

      return
      end
      function r4_huge ( )

c*********************************************************************72
c
cc R4_HUGE returns a "huge" R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_HUGE, a huge number.
c
      implicit none

      real r4_huge

      r4_huge = 1.0E+30

      return
      end
      function r4_in_01 ( a )

c*********************************************************************72
c
cc R4_IN_01 is TRUE if an R4 is in the range [0,1].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A, the value.
c
c    Output, logical R4_IN_01, is TRUE if 0 <= A <= 1.
c
      implicit none

      real a
      logical r4_in_01
      logical value

      if ( a .lt. 0.0E+00 .or. 1.0E+00 .lt. a ) then
        value = .false.
      else
        value = .true.
      end if

      r4_in_01 = value

      return
      end
      function r4_is_int ( r )

c*********************************************************************72
c
cc R4_IS_INT determines if an R4 represents an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R, the number to be checked.
c
c    Output, logical R4_IS_INT, is TRUE if R is an integer value.
c
      implicit none

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      real r
      logical r4_is_int

      if ( real ( i4_huge ) .lt. r ) then
        r4_is_int = .false.
      else if ( r .lt. - real ( i4_huge ) ) then
        r4_is_int = .false.
      else if ( r .eq. real ( int ( r ) ) ) then
        r4_is_int = .true.
      else
        r4_is_int = .false.
      end if

      return
      end
      function r4_log_2 ( x )

c*********************************************************************72
c
cc R4_LOG_2 returns the logarithm base 2 of an R4.
c
c  Discussion:
c
c    value = Log ( |X| ) / Log ( 2.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose base 2 logarithm is desired.
c    X should not be 0.
c
c    Output, real R4_LOG_2, the logarithm base 2 of the absolute
c    value of X.  It should be true that |X| = 2**D_LOG_2.
c
      implicit none

      real r4_huge
      real r4_log_2
      real x

      if ( x .eq. 0.0E+00 ) then
        r4_log_2 = - r4_huge ( )
      else
        r4_log_2 = log ( abs ( x ) ) / log ( 2.0E+00 )
      end if

      return
      end
      function r4_log_10 ( x )

c*********************************************************************72
c
cc R4_LOG_10 returns the logarithm base 10 of an R4.
c
c  Discussion:
c
c    value = Log10 ( |X| )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose base 2 logarithm is desired.
c    X should not be 0.
c
c    Output, real R4_LOG_10, the logarithm base 10 of the absolute
c    value of X.  It should be true that |X| = 10**R_LOG_10.
c
      implicit none

      real r4_huge
      real r4_log_10
      real x

      if ( x .eq. 0.0E+00 ) then
        r4_log_10 = - r4_huge ( x )
      else
        r4_log_10 = log10 ( abs ( x ) )
      end if

      return
      end
      function r4_log_b ( x, b )

c*********************************************************************72
c
cc R4_LOG_B returns the logarithm base B of an R4.
c
c  Discussion:
c
c    value = log ( |X| ) / log ( |B| )
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
c    Input, real X, the number whose base B logarithm is desired.
c    X should not be 0.
c
c    Input, real B, the base, which should not be 0, 1 or -1.
c
c    Output, real R4_LOG_B, the logarithm base B of the absolute
c    value of X.  It should be true that |X| = |B|**R4_LOG_B.
c
      implicit none

      real b
      real r4_huge
      real r4_log_b
      real x

      if ( b .eq. 0.0E+00 .or.
     &     b .eq. 1.0E+00 .or.
     &     b .eq. - 1.0E+00 ) then
        r4_log_b = - r4_huge ( x )
      else if ( abs ( x ) .eq. 0.0E+00 ) then
        r4_log_b = - r4_huge ( x )
      else
        r4_log_b = log ( abs ( x ) ) / log ( abs ( b ) )
      end if

      return
      end
      subroutine r4_mant ( x, s, r, l )

c*********************************************************************72
c
cc R4_MANT computes the "mantissa" or "fraction part" of an R4.
c
c  Discussion:
c
c    X = S * R * 2**L
c
c    S is +1 or -1,
c    R is an R4 value between 1.0 and 2.0,
c    L is an integer.
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
c    Input, real X, the number to be decomposed.
c
c    Output, integer S, the "sign" of the number.
c    S will be -1 if X is less than 0, and +1 if X is greater
c    than or equal to zero.
c
c    Output, real R, the mantissa of X.  R will be greater
c    than or equal to 1, and strictly less than 2.  The one
c    exception occurs if X is zero, in which case R will also
c    be zero.
c
c    Output, integer L, the integer part of the logarithm
c    (base 2) of X.
c
      implicit none

      integer l
      real r
      integer s
      real x
c
c  Determine the sign.
c
      if ( x .lt. 0.0E+00 ) then
        s = -1
      else
        s = 1
      end if
c
c  Set R to the absolute value of X, and L to zero.
c  Then force R to lie between 1 and 2.
c
      if ( x .lt. 0.0E+00 ) then
        r = -x
      else
        r = x
      end if

      l = 0
c
c  Time to bail out if X is zero.
c
      if ( x .eq. 0.0E+00 ) then
        return
      end if

10    continue

      if ( 2.0E+00 .le. r ) then
        r = r / 2.0E+00
        l = l + 1
        go to 10
      end if

20    continue

      if ( r .lt. 1.0E+00 ) then
        r = r * 2.0E+00
        l = l - 1
        go to 20
      end if

      return
      end
      function r4_mod ( x, y )

c*********************************************************************72
c
cc R4_MOD returns the remainder of R4 division.
c
c  Discussion:
c
c    If
c      REM = R4_MOD ( X, Y )
c      RMULT = ( X - REM ) / Y
c    then
c      X = Y * RMULT + REM
c    where REM has the same sign as X, and abs ( REM ) < Y.
c
c  Example:
c
c        X         Y     R4_MOD  R4_MOD Factorization
c
c      107        50       7      107 =  2 *  50 + 7
c      107       -50       7      107 = -2 * -50 + 7
c     -107        50      -7     -107 = -2 *  50 - 7
c     -107       -50      -7     -107 =  2 * -50 - 7
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
c    Input, real X, the number to be divided.
c
c    Input, real Y, the number that divides X.
c
c    Output, real R4_MOD, the remainder when X is divided by Y.
c
      implicit none

      real r4_mod
      real x
      real y

      if ( y .eq. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_MOD - Fatal error!'
        write ( *, '(a,g14.6)' ) '  R4_MOD ( X, Y ) called with Y = ', y
        stop
      end if

      r4_mod = x - real ( int ( x / y ) ) * y

      if ( x .lt. 0.0E+00 .and. 0.0E+00 .lt. r4_mod ) then
        r4_mod = r4_mod - abs ( y )
      else if ( 0.0E+00 .lt. x .and. r4_mod .lt. 0.0E+00 ) then
        r4_mod = r4_mod + abs ( y )
      end if

      return
      end
      function r4_modp ( x, y )

c*********************************************************************72
c
cc R4_MODP returns the nonnegative remainder of R4 division.
c
c  Discussion:
c
c    If
c      REM = R4_MODP ( X, Y )
c      RMULT = ( X - REM ) / Y
c    then
c      X = Y * RMULT + REM
c    where REM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360.0) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, R4_MODP(A,360.0) is between 0 and 360, always.
c
c  Example:
c
c        X         Y     MOD R4_MODP  R4_MODP Factorization
c
c      107        50       7       7    107 =  2 *  50 + 7
c      107       -50       7       7    107 = -2 * -50 + 7
c     -107        50      -7      43   -107 = -3 *  50 + 43
c     -107       -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number to be divided.
c
c    Input, real Y, the number that divides X.
c
c    Output, real R4_MODP, the nonnegative remainder
c    when X is divided by Y.
c
      implicit none

      real r4_modp
      real x
      real y

      if ( y .eq. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_MODP - Fatal error!'
        write ( *, '(a,g14.6)' )
     &    '  R4_MODP ( X, Y ) called with Y = ', y
        stop
      end if

      r4_modp = mod ( x, y )

      if ( r4_modp .lt. 0.0E+00 ) then
        r4_modp = r4_modp + abs ( y )
      end if

      return
      end
      function r4_mop ( i )

c*********************************************************************72
c
cc R4_MOP returns the I-th power of -1 as an R4.
c
c  Discussion:
c
c    An R4 is a real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, real R4_MOP, the I-th power of -1.
c
      implicit none

      integer i
      real r4_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        r4_mop = + 1.0E+00
      else
        r4_mop = - 1.0E+00
      end if

      return
      end
      function r4_nint ( x )

c*********************************************************************72
c
cc R4_NINT returns the nearest integer to an R4.
c
c  Example:
c
c        X        R4_NINT
c
c      1.3         1
c      1.4         1
c      1.5         1 or 2
c      1.6         2
c      0.0         0
c     -0.7        -1
c     -1.1        -1
c     -1.6        -2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the value.
c
c    Output, integer R4_NINT, the nearest integer to X.
c
      implicit none

      integer r4_nint
      integer s
      real x

      if ( x .lt. 0.0E+00 ) then
        s = - 1
      else
        s = + 1
      end if

      r4_nint = s * int ( abs ( x ) + 0.5E+00 )

      return
      end
      function r4_normal ( a, b, seed )

c*********************************************************************72
c
cc R4_NORMAL returns a scaled pseudonormal R4.
c
c  Discussion:
c
c    The normal probability distribution function (PDF) is sampled,
c    with mean A and standard deviation B.
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
c    Input, real A, the mean of the PDF.
c
c    Input, real B, the standard deviation of the PDF.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, real R4_NORMAL, a sample of the normal PDF.
c
      implicit none

      real a
      real b
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r1
      real r2
      real r4_normal
      real r4_uniform_01
      integer seed
      integer seed2
      integer used
      real x
      real y

      save seed2
      save used
      save y

      data seed2 / 0 /
      data used / 0 /
      data y / 0.0E+00 /
c
c  On odd numbered calls, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        r1 = r4_uniform_01 ( seed )

        if ( r1 .eq. 0.0E+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R4_NORMAL - Fatal error!'
          write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed
        r2 = r4_uniform_01 ( seed2 )

        x = sqrt ( -2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * pi * r2 )
        y = sqrt ( -2.0E+00 * log ( r1 ) ) * sin ( 2.0E+00 * pi * r2 )
c
c  On odd calls, return the second normal and its corresponding seed.
c
      else

        seed = seed2
        x = y

      end if

      used = used + 1

      r4_normal = a + b * x

      return
      end
      function r4_normal_01 ( seed )

c*********************************************************************72
c
cc R4_NORMAL_01 returns a unit pseudonormal real R4.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, essentially, the input value of
c    SEED is ignored, since the code saves the second normal random value.
c
c    If you didn't know this, you might be confused since, usually, the
c    output of a random number generator can be completely controlled by
c    the input value of the SEED.  If I were more careful, I could rewrite
c    this routine so that it would distinguish between cases where the input
c    value of SEED is the output value from the previous call (all is well)
c    and those cases where it is not (the user has decided to do something
c    new.  Restart the uniform random number sequence.)  But I'll leave
c    that for later.
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
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, real R4_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r1
      real r2
      real r4_normal_01
      real r4_uniform_01
      integer seed
      integer seed2
      integer used
      real x
      real y

      save seed2
      save used
      save y

      data seed2 / 0 /
      data used / 0 /
      data y / 0.0E+00 /
c
c  On odd numbered calls, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        r1 = r4_uniform_01 ( seed )

        if ( r1 .eq. 0.0E+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R4_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed
        r2 = r4_uniform_01 ( seed2 )

        x = sqrt ( -2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * pi * r2 )
        y = sqrt ( -2.0E+00 * log ( r1 ) ) * sin ( 2.0E+00 * pi * r2 )
c
c  On odd calls, return the second normal and its corresponding seed.
c
      else

        seed = seed2
        x = y

      end if

      used = used + 1

      r4_normal_01 = x

      return
      end
      function r4_pi ( )

c*********************************************************************72
c
cc R4_PI returns the value of pi as an R4.
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
c    Output, real R4_PI, the value of pi.
c
      implicit none

      real r4_pi

      r4_pi = 3.1415926E+00

      return
      end
      function r4_power ( r, p )

c*********************************************************************72
c
cc R4_POWER computes the P-th power of an R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R, the base.
c
c    Input, integer P, the power, which may be negative.
c
c    Output, real R4_POWER, the value of the P-th power of R.
c
      implicit none

      integer p
      real r
      real r4_power
      real value
c
c  Special case.  R^0 = 1.
c
      if ( p .eq. 0 ) then

        value = 1.0E+00
c
c  Special case.  Positive powers of 0 are 0.
c  For negative powers of 0, we go ahead and compute R**P,
c  relying on the software to complain.
c
      else if ( r .eq. 0.0E+00 ) then

        if ( 0 .lt. p ) then
          value = 0.0E+00
        else
          value = r**p
        end if

      else if ( 1 .le. p ) then
        value = r**p
      else
        value = 1.0E+00 / r**(-p)
      end if

      r4_power = value

      return
      end
      subroutine r4_power_fast ( r, p, rp, mults )

c*********************************************************************72
c
cc R4_POWER_FAST computes an integer power of an R4.
c
c  Discussion:
c
c    Obviously, R**P can be computed using P-1 multiplications.
c
c    However, R**P can also be computed using at most 2*LOG2(P) multiplications.
c    To do the calculation this way, let N = LOG2(P).
c    Compute A, A**2, A**4, ..., A**N by N-1 successive squarings.
c    Start the value of R**P at A, and each time that there is a 1 in
c    the binary expansion of P, multiply by the current result of the squarings.
c
c    This algorithm is not optimal.  For small exponents, and for special
c    cases, the result can be computed even more quickly.
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
c    Input, real R, the base.
c
c    Input, integer P, the power, which may be negative.
c
c    Output, real RP, the value of R**P.
c
c    Output, integer MULTS, the number of multiplications
c    and divisions.
c
      implicit none

      integer mults
      integer p
      integer p_mag
      integer p_sign
      real r
      real r2
      real rp

      mults = 0
c
c  Special bases.
c
      if ( r .eq. 1.0E+00 ) then
        rp = 1.0E+00
        return
      end if

      if ( r .eq. -1.0E+00 ) then

        if ( mod ( p, 2 ) .eq. 1 ) then
          rp = -1.0E+00
        else
          rp = 1.0E+00
        end if

        return

      end if

      if ( r .eq. 0.0E+00 ) then

        if ( p .le. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R4_POWER_FAST - Fatal error!'
          write ( *, '(a)' )
     &      '  Base R is zero, and exponent is negative.'
          write ( *, '(a,i8)' ) '  Exponent P = ', p
          stop
        end if

        rp = 0.0E+00
        return

      end if
c
c  Special powers.
c
      if ( p .eq. -1 ) then
        rp = 1.0E+00 / r
        mults = mults + 1
        return
      else if ( p .eq. 0 ) then
        rp = 1.0E+00
        return
      else if ( p .eq. 1 ) then
        rp = r
        return
      end if
c
c  Some work to do.
c
      p_mag = abs ( p )
      p_sign = sign ( 1, p )

      rp = 1.0E+00
      r2 = r

10    continue

      if ( 0 .lt. p_mag ) then

        if ( mod ( p_mag, 2 ) .eq. 1 ) then
          rp = rp * r2
          mults = mults + 1
        end if

        p_mag = p_mag / 2
        r2 = r2 * r2
        mults = mults + 1

        go to 10

      end if

      if ( p_sign .eq. -1 ) then
        rp = 1.0E+00 / rp
        mults = mults + 1
      end if

      return
      end
      function r4_pythag ( a, b )

c*********************************************************************72
c
cc R4_PYTHAG computes sqrt ( A * A + B * B ), avoiding overflow and underflow.
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
c    Input, real A, B, the values for which sqrt ( A * A + B * B )
c    is desired.
c
c    Output, real R4_PYTHAG, the value of sqrt ( A * A + B * B ).
c
      implicit none

      real a
      real a_abs
      real b
      real b_abs
      real r4_pythag

      a_abs = abs ( a )
      b_abs = abs ( b )

      if ( b_abs .lt. a_abs ) then
        r4_pythag = a_abs *
     &    sqrt ( 1.0E+00 + ( b_abs / a_abs ) * ( b_abs / a_abs ) )
      else if ( b_abs .eq. 0.0E+00 ) then
        r4_pythag = 0.0E+00
      else if ( a_abs .le. b_abs ) then
        r4_pythag = b_abs *
     &    sqrt ( 1.0E+00 + ( a_abs / b_abs ) * ( a_abs / b_abs ) )
      end if

      return
      end
      subroutine r4_round2 ( nplace, x, xround )

c*********************************************************************72
c
cc R4_ROUND2 rounds an R4 to a specified number of binary digits.
c
c  Discussion:
c
c    Assume that the input quantity X has the form
c
c      X = S * J * 2^L
c
c    where S is plus or minus 1, L is an integer, and J is a binary
c    mantissa which is either exactly zero, or greater than or equal
c    to 0.5 and strictly less than 1.0.
c
c    Then on return, XROUND will satisfy
c
c      XROUND = S * K * 2^L
c
c    where S and L are unchanged, and K is a binary mantissa which
c    agrees with J in the first NPLACE binary digits and is zero
c    thereafter.
c
c    If NPLACE is 0, XROUND will always be zero.
c
c    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
c
c    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
c    or 0.75.
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
c    Input, integer NPLACE, the number of binary digits to
c    preserve.  NPLACE should be 0 or positive.
c
c    Input, real X, the number to be decomposed.
c
c    Output, real XROUND, the rounded value of X.
c
      implicit none

      integer iplace
      integer l
      integer nplace
      integer s
      real x
      real xmant
      real xround
      real xtemp

      xround = 0.0E+00
c
c  1: Handle the special case of 0.
c
      if ( x .eq. 0.0E+00 ) then
        return
      end if

      if ( nplace .le. 0 ) then
        return
      end if
c
c  2: Determine the sign S.
c
      if ( 0.0E+00 .lt. x ) then
        s = 1
        xtemp = x
      else
        s = -1
        xtemp = -x
      end if
c
c  3: Force XTEMP to lie between 1 and 2, and compute the
c  logarithm L.
c
      l = 0

10    continue

      if ( 2.0E+00 .le. xtemp ) then
        xtemp = xtemp / 2.0E+00
        l = l + 1
        go to 10
      end if

20    continue

      if ( xtemp .lt. 1.0E+00 ) then
        xtemp = xtemp * 2.0E+00
        l = l - 1
        go to 20
      end if
c
c  4: Strip out the digits of the mantissa as XMANT, and decrease L.
c
      xmant = 0.0E+00
      iplace = 0

30    continue

        xmant = 2.0E+00 * xmant

        if ( 1.0E+00 .le. xtemp ) then
          xmant = xmant + 1.0E+00
          xtemp = xtemp - 1.0E+00
        end if

        iplace = iplace + 1

        if ( xtemp .eq. 0.0E+00 .or. nplace .le. iplace ) then
          xround = s * xmant * 2.0E+00**l
          go to 40
        end if

        l = l - 1
        xtemp = xtemp * 2.0E+00

      go to 30

40    continue

      return
      end
      subroutine r4_roundb ( base, nplace, x, xround )

c*********************************************************************72
c
cc R4_ROUNDB rounds an R4 to a given number of digits in a given base.
c
c  Discussion:
c
c    The code does not seem to do a good job of rounding when
c    the base is negativec
c
c    Assume that the input quantity X has the form
c
c      X = S * J * BASE^L
c
c    where S is plus or minus 1, L is an integer, and J is a
c    mantissa base BASE which is either exactly zero, or greater
c    than or equal to (1/BASE) and less than 1.0.
c
c    Then on return, XROUND will satisfy
c
c      XROUND = S * K * BASE^L
c
c    where S and L are unchanged, and K is a mantissa base BASE
c    which agrees with J in the first NPLACE digits and is zero
c    thereafter.
c
c    Note that because of rounding, for most bases, most numbers
c    with a fractional quantities cannot be stored exactly in the
c    computer, and hence will have trailing "bogus" digits.
c
c    If NPLACE is 0, XROUND will always be zero.
c
c    If NPLACE is 1, the mantissa of XROUND will be 0,
c    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
c
c    If NPLACE is 2, the mantissa of XROUND will be 0,
c    BASE/BASE^2, (BASE+1)/BASE^2, ...,
c    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
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
c    Input, integer BASE, the base of the arithmetic.
c    BASE must not be zero.  Theoretically, BASE may be negative.
c
c    Input, integer NPLACE, the number of digits base BASE to
c    preserve.  NPLACE should be 0 or positive.
c
c    Input, real X, the number to be decomposed.
c
c    Output, real XROUND, the rounded value of X.
c
      implicit none

      integer base
      integer iplace
      integer is
      integer js
      integer l
      integer nplace
      real x
      real xmant
      real xround
      real xtemp

      xround = 0.0E+00
c
c  0: Error checks.
c
      if ( base .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_ROUNDB - Fatal error!'
        write ( *, '(a)' ) '  The base BASE cannot be zero.'
        stop
      end if
c
c  1: Handle the special case of 0.
c
      if ( x .eq. 0.0E+00 ) then
        return
      end if

      if ( nplace .le. 0 ) then
        return
      end if
c
c  2: Determine the sign IS.
c
      if ( 0.0E+00 .lt. x ) then
        is = 1
        xtemp = x
      else
        is = -1
        xtemp = -x
      end if
c
c  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
c  logarithm L.
c
      l = 0

10    continue

      if ( abs ( base ) .le. abs ( xtemp ) ) then

        xtemp = xtemp / real ( base )

        if ( xtemp .lt. 0.0E+00 ) then
          is = -is
          xtemp = -xtemp
        end if

        l = l + 1

        go to 10

      end if

20    continue

      if ( abs ( xtemp ) .lt. 1.0E+00 ) then

        xtemp = xtemp * base

        if ( xtemp .lt. 0.0E+00 ) then
          is = -is
          xtemp = -xtemp
        end if

        l = l - 1

        go to 20

      end if
c
c  4: Now strip out the digits of the mantissa as XMANT, and
c  decrease L.
c
      xmant = 0.0E+00
      iplace = 0
      js = is

30    continue

        xmant = base * xmant

        if ( xmant .lt. 0.0E+00 ) then
          js = -js
          xmant = -xmant
        end if

        if ( 1.0E+00 .le. xtemp ) then
          xmant = xmant + int ( xtemp )
          xtemp = xtemp - int ( xtemp )
        end if

        iplace = iplace + 1

        if ( xtemp .eq. 0.0E+00 .or. nplace .le. iplace ) then
          xround = js * xmant * real ( base )**l
          go to 40
        end if

        l = l - 1
        xtemp = xtemp * base

        if ( xtemp .lt. 0.0E+00 ) then
          is = -is
          xtemp = -xtemp
        end if

      go to 30

40    continue

      return
      end
      subroutine r4_roundx ( nplace, x, xround )

c*********************************************************************72
c
cc R4_ROUNDX rounds an R4.
c
c  Discussion:
c
c    Assume that the input quantity X has the form
c
c      X = S * J * 10^L
c
c    where S is plus or minus 1, L is an integer, and J is a decimal
c    mantissa which is either exactly zero, or greater than or equal
c    to 0.1 and less than 1.0.
c
c    Then on return, XROUND will satisfy
c
c      XROUND = S * K * 10^L
c
c    where S and L are unchanged, and K is a decimal mantissa which
c    agrees with J in the first NPLACE decimal digits and is zero
c    thereafter.
c
c    Note that because of rounding, most decimal fraction quantities
c    cannot be stored exactly in the computer, and hence will have
c    trailing "bogus" digits.
c
c    If NPLACE is 0, XROUND will always be zero.
c
c    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
c    0.2, ..., or 0.9.
c
c    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
c    0.03, ..., 0.98, 0.99.
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
c    Input, integer NPLACE, the number of decimal digits to
c    preserve.  NPLACE should be 0 or positive.
c
c    Input, real X, the number to be decomposed.
c
c    Output, real XROUND, the rounded value of X.
c
      implicit none

      integer iplace
      integer is
      integer l
      integer nplace
      real x
      real xmant
      real xround
      real xtemp

      xround = 0.0E+00
c
c  1: Handle the special case of 0.
c
      if ( x .eq. 0.0E+00 ) then
        return
      end if

      if ( nplace .le. 0 ) then
        return
      end if
c
c  2: Determine the sign IS.
c
      if ( 0.0E+00 .lt. x ) then
        is = 1
        xtemp = x
      else
        is = -1
        xtemp = -x
      end if
c
c  3: Force XTEMP to lie between 1 and 10, and compute the
c  logarithm L.
c
      l = 0

10    continue

      if ( 10.0E+00 .le. x ) then
        xtemp = xtemp / 10.0E+00
        l = l + 1
        go to 10
      end if

20    continue

      if ( xtemp .lt. 1.0E+00 ) then
        xtemp = xtemp * 10.0E+00
        l = l - 1
        go to 20
      end if
c
c  4: Now strip out the digits of the mantissa as XMANT, and
c  decrease L.
c
      xmant = 0.0E+00
      iplace = 0

30    continue

        xmant = 10.0E+00 * xmant

        if ( 1.0E+00 .le. xtemp ) then
          xmant = xmant + int ( xtemp )
          xtemp = xtemp - int ( xtemp )
        end if

        iplace = iplace + 1

        if ( xtemp .eq. 0.0E+00 .or. nplace .le. iplace ) then
          xround = is * xmant * ( 10.0E+00**l )
          go to 40
        end if

        l = l - 1
        xtemp = xtemp * 10.0E+00

      go to 30

40    continue

      return
      end
      function r4_sign ( x )

c*********************************************************************72
c
cc R4_SIGN returns the sign of an R4.
c
c  Discussion:
c
c    value = -1 if X < 0;
c    value =  0 if X => 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose sign is desired.
c
c    Output, real R4_SIGN, the sign of X:
c
      implicit none

      real r4_sign
      real x

      if ( x .lt. 0.0E+00 ) then
        r4_sign = -1.0E+00
      else
        r4_sign = +1.0E+00
      end if

      return
      end
      function r4_sign_opposite ( r1, r2 )

c*********************************************************************72
c
cc R4_SIGN_OPPOSITE is TRUE if two R4's are not of the same sign.
c
c  Discussion:
c
c    This test could be coded numerically as
c
c      if ( r1 * r2 <= 0.0 ) then ...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R1, R2, the values to check.
c
c    Output, logical R4_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
c    or ( R2 <= 0 and 0 <= R1 ).
c
      implicit none

      real r1
      real r2
      logical r4_sign_opposite

      r4_sign_opposite = ( r1 <= 0.0E+00 .and. 0.0E+00 <= r2 ) .or.
     &                   ( r2 <= 0.0E+00 .and. 0.0E+00 <= r1 )

      return
      end
      function r4_sign_opposite_strict ( r1, r2 )

c*********************************************************************72
c
cc R4_SIGN_OPPOSITE_STRICT is TRUE if two R4's are strictly of opposite sign.
c
c  Discussion:
c
c    This test could be coded numerically as
c
c      if ( r1 * r2 < 0.0 ) then ...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R1, R2, the values to check.
c
c    Output, logical R4_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
c    or ( R2 < 0 and 0 < R1 ).
c
      implicit none

      real r1
      real r2
      logical r4_sign_opposite_strict

      r4_sign_opposite_strict = ( r1 < 0.0E+00 .and. 0.0E+00 < r2 ) .or.
     &                          ( r2 < 0.0E+00 .and. 0.0E+00 < r1 )

      return
      end
      subroutine r4_swap ( x, y )

c*********************************************************************72
c
cc R4_SWAP swaps two R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      real x
      real y
      real z

      z = x
      x = y
      y = z

      return
      end
      subroutine r4_swap3 ( x, y, z )

c*********************************************************************72
c
cc R4_SWAP3 swaps three R4's.
c
c  Example:
c
c    Input:
c
c      X = 1, Y = 2, Z = 3
c
c    Output:
c
c      X = 2, Y = 3, Z = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real X, Y, Z, three values to be swapped.
c
      implicit none

      real w
      real x
      real y
      real z

      w = x
      x = y
      y = z
      z = w

      return
      end
      function r4_tiny ( )

c*********************************************************************72
c
cc R4_TINY returns the smallest positive R4.
c
c  Discussion:
c
c    FORTRAN90 provides a built-in routine TINY ( X ) that
c    is more suitable for this purpose.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_TINY, a "tiny" value.
c
      implicit none

      real r4_tiny

      r4_tiny = 0.1175494350822E-37

      return
      end
      subroutine r4_to_r4_discrete ( r, rmin, rmax, nr, rd )

c*********************************************************************72
c
cc R4_TO_R4_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
c
c  Formula:
c
c    if ( R < RMIN ) then
c      RD = RMIN
c    else if ( RMAX < R ) then
c      RD = RMAX
c    else
c      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
c      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
c
c    In the special case where NR = 1, when
c
c      XD = 0.5 * ( RMAX + RMIN )
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
c    Input, real R, the number to be converted.
c
c    Input, real RMAX, RMIN, the maximum and minimum
c    values for RD.
c
c    Input, integer NR, the number of allowed values for XD.
c    NR should be at least 1.
c
c    Output, real RD, the corresponding discrete value.
c
      implicit none

      integer f
      integer nr
      real r
      real rd
      real rmax
      real rmin
c
c  Check for errors.
c
      if ( nr .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_TO_R4_DISCRETE - Fatal error!'
        write ( *, '(a,i8)' ) '  NR = ', nr
        write ( *, '(a)' ) '  but NR must be at least 1.'
        stop
      end if

      if ( nr .eq. 1 ) then
        rd = 0.5E+00 * ( rmin + rmax )
        return
      end if

      if ( rmax .eq. rmin ) then
        rd = rmax
        return
      end if

      f = nint ( real ( nr ) * ( rmax - r ) / ( rmax - rmin ) )
      f = max ( f, 0 )
      f = min ( f, nr )

      rd = ( real (      f ) * rmin
     &     + real ( nr - f ) * rmax )
     &     / real ( nr     )

      return
      end
      subroutine r4_to_dhms ( r, d, h, m, s )

c*********************************************************************72
c
cc R4_TO_DHMS converts decimal days into days, hours, minutes, seconds.
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
c    Input, real R, a decimal number representing a time
c    period measured in days.
c
c    Output, integer D, H, M, S, the equivalent number of days,
c    hours, minutes and seconds.
c
      implicit none

      integer d
      integer h
      integer m
      real r
      real r_copy
      integer s

      r_copy = abs ( r )

      d = int ( r_copy )

      r_copy = r_copy - d
      r_copy = 24.0E+00 * r_copy
      h = int ( r_copy )

      r_copy = r_copy - h
      r_copy = 60.0E+00 * r_copy
      m = int ( r_copy )

      r_copy = r_copy - m
      r_copy = 60.0E+00 * r_copy
      s = int ( r_copy )

      if ( r .lt. 0.0E+00 ) then
        d = -d
        h = -h
        m = -m
        s = -s
      end if

      return
      end
      subroutine r4_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

c*********************************************************************72
c
cc R4_TO_I4 maps X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
c
c  Formula:
c
c    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
c    IX := min ( IX, max ( IXMIN, IXMAX ) )
c    IX := max ( IX, min ( IXMIN, IXMAX ) )
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
c    Input, real X, the number to be converted.
c
c    Input, real XMIN, XMAX, the range.  XMAX and
c    XMIN must not be equal.  It is not necessary that XMIN be less than XMAX.
c
c    Input, integer IXMIN, IXMAX, the allowed range of the output
c    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
c    It is not necessary that IXMIN be less than IXMAX.
c
c    Output, integer IX, the value in the range [IXMIN,IXMAX] that
c    corresponds to X.
c
      implicit none

      integer ix
      integer ixmax
      integer ixmin
      real temp
      real x
      real xmax
      real xmin

      if ( xmax .eq. xmin ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  XMAX = XMIN, making a zero divisor.'
        write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
        write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
        stop
      end if

      temp =
     &    ( ( xmax - x        ) * real ( ixmin )
     &    + (        x - xmin ) * real ( ixmax ) )
     &    / ( xmax     - xmin )

      if ( 0.0E+00 .le. temp ) then
        temp = temp + 0.5E+00
      else
        temp = temp - 0.5E+00
      end if

      ix = int ( temp )

      return
      end
      function r4_uniform ( a, b, seed )

c*********************************************************************72
c
cc R4_UNIFORM returns a scaled pseudorandom R4.
c
c  Discussion:
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
c    29 January 2005
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
c    Input, real A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM, a number strictly between A and B.
c
      implicit none

      real a
      real b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r4_uniform

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM - Fatal error!'
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
      r4_uniform = a + ( b - a )
     &  * real ( dble ( seed ) * 4.656612875E-10 )

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
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
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875E-10 )

      return
      end
      subroutine r4_unswap3 ( x, y, z )

c*********************************************************************72
c
cc R4_UNSWAP3 unswaps three R4's.
c
c  Example:
c
c    Input:
c
c      X = 2, Y = 3, Z = 1
c
c    Output:
c
c      X = 1, Y = 2, Z = 3
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
c    Input/output, real X, Y, Z, three values to be swapped.
c
      implicit none

      real w
      real x
      real y
      real z

      w = z
      z = y
      y = x
      x = w

      return
      end
      function r4_walsh_1d ( x, digit )

c*********************************************************************72
c
cc R4_WALSH_1D evaluates the Walsh function.
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
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument of the Walsh function.
c
c    Input, integer DIGIT, the index of the Walsh function.
c
c    Output, real R4_WALSH_1D, the value of the Walsh function.
c
      implicit none

      integer digit
      integer n
      real r4_walsh_1d
      real x
      real x_copy
c
c  Hide the effect of the sign of X.
c
      x_copy = abs ( x )
c
c  If DIGIT is positive, divide by 2 DIGIT times.
c  If DIGIT is negative, multiply by 2 (-DIGIT) times.
c
      x_copy = x_copy / 2.0E+00**digit
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
        r4_walsh_1d = 0.0E+00
      else
        r4_walsh_1d = 1.0E+00
      end if

      return
      end
      function r4_wrap ( r, rlo, rhi )

c*********************************************************************72
c
cc R4_WRAP forces an R4 to lie between given limits by wrapping.
c
c  Discussion:
c
c    An R4 is a real value.
c
c  Example:
c
c    RLO = 4.0, RHI = 8.0
c
c     R  Value
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
c    04 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real R, a value.
c
c    Input, real RLO, RHI, the desired bounds.
c
c    Output, real R4_WRAP, a "wrapped" version of the value.
c
      implicit none

      integer n
      real r
      real r4_wrap
      real rhi
      real rhi2
      real rlo
      real rlo2
      real rwide
      real value
c
c  Guarantee RLO2 .lt. RHI2.
c
      rlo2 = min ( rlo, rhi )
      rhi2 = max ( rlo, rhi )
c
c  Find the width.
c
      rwide = rhi2 - rlo2
c
c  Add enough copies of (RHI2-RLO2) to R so that the
c  result ends up in the interval RLO2 - RHI2.
c
      if ( rwide .eq. 0.0E+00 ) then
        value = rlo
      else if ( r .lt. rlo2 ) then
        n = int ( ( rlo2 - r ) / rwide ) + 1
        value = r + n * rwide
        if ( value .eq. rhi ) then
          value = rlo
        end if
      else
        n = int ( ( r - rlo2 ) / rwide )
        value = r - n * rwide
        if ( value .eq. rlo ) then
          value = rhi
        end if
      end if

      r4_wrap = value

      return
      end
      subroutine r42_cheby ( n, alo, ahi, a )

c*********************************************************************72
c
cc R42_CHEBY sets up the Chebyshev abscissas in an R4 interval.
c
c  Discussion:
c
c    The routine sets up a vector of X values spaced between the values
c    XLO and XHI in a similar way to the spacing of the Chebyshev
c    points of the same order in the interval [-1,1].
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
c    Input, integer N, the number of points to compute.
c
c    Input, real ALO, AHI, the range.
c
c    Output, real A(N), the computed X values.
c
      implicit none

      integer n

      real a(n)
      real ahi
      real alo
      real arg
      integer i
      real pi
      parameter ( pi = 3.141592653589793E+00 )

      if ( n .eq. 1 ) then

        a(1) = 0.5E+00 * ( alo + ahi )

      else if ( 1 .lt. n ) then

        do i = 1, n

          arg = real ( 2 * i - 1 ) * pi / real ( 2 * n )

          a(i) = 0.5E+00 * ( ( 1.0E+00 + cos ( arg ) ) * alo
     &                     + ( 1.0E+00 - cos ( arg ) ) * ahi )

        end do

      end if

      return
      end
      function r42_dist_l2 ( a1, a2 )

c*********************************************************************72
c
cc R42_DIST_L2 returns the L2 distance between a pair of R42's.
c
c  Discussion:
c
c    An R42 is a vector of type R4, with two entries.
c
c    The vector L2 norm is defined as:
c
c      sqrt ( sum ( 1 <= I <= N ) A(I) * A(I) ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 September 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A1(2), A2(2), the vectors.
c
c    Output, real R42_DIST_L2, the L2 norm of the distance
c    between A1 and A2.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      real a1(dim_num)
      real a2(dim_num)
      real r42_dist_l2

      r42_dist_l2 = sqrt ( 
     &  ( a1(1) - a2(1) )**2 + ( a1(2) - a2(2) )**2 )

      return
      end
      subroutine r4col_compare ( m, n, a, i, j, value )

c*********************************************************************72
c
cc R4COL_COMPARE compares columns in an R4COL.
c
c  Discussion:
c
c    An R4COL is an M by N array of R4 values.
c    It is regarded as an array of N columns of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, I = 2, J = 4
c
c      A = (
c        1.  2.  3.  4.
c        5.  6.  7.  8.
c        9. 10. 11. 12. )
c
c    Output:
c
c      VALUE = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, real A(M,N), the M by N array.
c
c    Input, integer I, J, the columns to be compared.
c    I and J must be between 1 and N.
c
c    Output, integer VALUE, the results of the comparison:
c    -1, column I < column J,
c     0, column I = column J,
c    +1, column J < column I.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer isgn
      integer j
      integer k
      integer value
c
c  Check.
c
      if ( i .lt. 1 .or. n .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4COL_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Column index I is out of bounds.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      end if

      if ( j .lt. 1 .or. n .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4COL_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Column index J is out of bounds.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      end if

      value = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. m ) then

        if ( a(k,i) .lt. a(k,j) ) then
          value = -1
          return
        else if ( a(k,j) .lt. a(k,i) ) then
          value = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine r4mat_border_add ( m, n, table, table2 )

c*********************************************************************72
c
cc R4MAT_BORDER_ADD adds a "border" to an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4's.
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
c    Input, real TABLE(M,N), the table data.
c
c    Output, real TABLE2(M+2,N+2), the augmented table data.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      real table(m,n)
      real table2(m+2,n+2)

      do j = 1, n + 2
        table2(1,j) = 0.0E+00
      end do

      do j = 1, n + 2
        table2(m+2,j) = 0.0E+00
      end do

      do i = 2, m + 1
        table2(i,1) = 0.0E+00
      end do

      do i = 2, m + 1
        table2(i,n+2) = 0.0E+00
      end do

      do j = 2, n + 1
        do i = 2, m + 1
          table2(i,j) = table(i-1,j-1)
        end do
      end do

      return
      end
      subroutine r4mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc R4MAT_COPY copies an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Input, real A1(M,N), the matrix to be copied.
c
c    Output, real A2(M,N), a copy of the matrix.
c
      implicit none

      integer m
      integer n

      real a1(m,n)
      real a2(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
        end do
      end do

      return
      end
      subroutine r4mat_identity ( n, a )

c*********************************************************************72
c
cc R4MAT_IDENTITY stores the identity matrix in an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of A.
c
c    Output, real A(N,N), the N by N identity matrix.
c
      implicit none

      integer n

      real a(n,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, n
          if ( i .eq. j ) then
            a(i,i) = 1.0E+00
          else
            a(i,j) = 0.0E+00
          end if
        end do
      end do

      return
      end
      subroutine r4mat_indicator ( m, n, table )

c*********************************************************************72
c
cc R4MAT_INDICATOR sets up an "indicator" R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4's.
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
c    22 January 2011
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
c    Output, real TABLE(M,N), the table.
c
      implicit none

      integer m
      integer n

      integer fac
      integer i
      integer i4_log_10
      integer j
      real table(m,n)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, m
        do j = 1, n
          table(i,j) = real ( fac * i + j )
        end do
      end do

      return
      end
      function r4mat_max ( m, n, a )

c*********************************************************************72
c
cc R4MAT_MAX returns the maximum entry of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2008
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
c    Input, real A(M,N), the matrix.
c
c    Output, real R4MAT_MAX, the maximum entry of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_max
      real value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = max ( value, a(i,j) )
        end do
      end do

      r4mat_max = value

      return
      end
      function r4mat_min ( m, n, a )

c*********************************************************************72
c
cc R4MAT_MIN returns the minimum entry of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2008
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
c    Input, real A(M,N), the matrix.
c
c    Output, real R4MAT_MIN, the minimum entry of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_min
      real value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = min ( value, a(i,j) )
        end do
      end do

      r4mat_min = value

      return
      end
      subroutine r4mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R4MAT_MM multiplies two R4MAT's.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
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
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, real A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, real C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      real a(n1,n2)
      real b(n2,n3)
      real c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0E+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      return
      end
      function r4mat_norm_l1 ( m, n, a )

c*********************************************************************72
c
cc R4MAT_NORM_L1 returns the matrix L1 norm of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c    The matrix L1 norm is defined as:
c
c      R4MAT_NORM_L1 = max ( 1 <= J <= N )
c        sum ( 1 <= I <= M ) abs ( A(I,J) ).
c
c    The matrix L1 norm is derived from the vector L1 norm, and
c    satisifies:
c
c      r4vec_norm_l1 ( A * x ) <= r4mat_norm_l1 ( A ) * r4vec_norm_l1 ( x ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 November 2011
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
c    Input, real A(M,N), the matrix whose L1 norm is desired.
c
c    Output, real R4MAT_NORM_L1, the L1 norm of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_norm_l1
      real temp

      r4mat_norm_l1 = 0.0E+00

      do j = 1, n
        temp = 0.0E+00
        do i = 1, m
          temp = temp + abs ( a(i,j) )
        end do
        r4mat_norm_l1 = max ( r4mat_norm_l1, temp )
      end do

      return
      end
      subroutine r4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R4MAT_PRINT prints an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
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
c    Input, real A(M,N), the matrix.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      character*(*) title

      call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R4MAT_PRINT_SOME prints some of an R4MAT.
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
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, real A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      real a(m,n)
      character * ( 14 ) ctemp(incx)
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
      character * ( * ) title

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
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r4vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc R4VEC_COPY copies an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, real A1(N), the vector to be copied.
c
c    Output, real A2(N), a copy of A1.
c
      implicit none

      integer n

      real a1(n)
      real a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      function r4vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R4VEC_DOT_PRODUCT finds the dot product of a pair of R4VEC's.
c
c  Discussion:
c
c    An R4VEC is a vector of R4 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, real V1(N), V2(N), the vectors.
c
c    Output, real R4VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      real r4vec_dot_product
      real v1(n)
      real v2(n)
      real value

      value = 0.0E+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r4vec_dot_product = value

      return
      end
      subroutine r4vec_heap_a ( n, a )

c*********************************************************************72
c
cc R4VEC_HEAP_A reorders an R4VEC into an ascending heap.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
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
c    23 January 2011
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
c    Input/output, real A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer ifree
      real key
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
      subroutine r4vec_heap_d ( n, a )

c*********************************************************************72
c
cc R4VEC_HEAP_D reorders an R4VEC into an descending heap.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
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
c    23 January 2011
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
c    Input/output, real A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer ifree
      real key
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

20      continue
c
c  Once there is no more shifting to do, KEY moves into the free spot IFREE.
c
        a(ifree) = key

      end do

      return
      end
      subroutine r4vec_max ( n, a, amax )

c*********************************************************************72
c
cc R4VEC_MAX returns the maximum value in an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, real A(N), the array.
c
c    Output, real AMAX, the value of the largest entry.
c
      implicit none

      integer n

      real a(n)
      real amax
      integer i

      amax = a(1)
      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      subroutine r4vec_mean ( n, a, mean )

c*********************************************************************72
c
cc R4VEC_MEAN returns the mean of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real A(N), the vector whose mean is desired.
c
c    Output, real MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      real a(n)
      integer i
      real mean

      mean = 0.0E+00
      do i = 1, n
        mean = mean + a(i)
      end do
      mean = mean / real ( n )

      return
      end
      subroutine r4vec_min ( n, a, amin )

c*********************************************************************72
c
cc R4VEC_MIN returns the minimum value in an R4VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, real A(N), the array.
c
c    Output, real AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      real a(n)
      real amin
      integer i

      amin = a(1)
      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      function r4vec_norm ( n, a )

c*********************************************************************72
c
cc R4VEC_NORM returns the L2 norm of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4 values.
c
c    The vector L2 norm is defined as:
c
c      R4VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, real A(N), the vector whose L2 norm is desired.
c
c    Output, real R4VEC_NORM, the L2 norm of A.
c
      implicit none

      integer n

      real a(n)
      integer i
      real r4vec_norm
      real value

      value = 0.0E+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r4vec_norm = value

      return
      end
      function r4vec_norm_l1 ( n, a )

c*********************************************************************72
c
cc R4VEC_NORM_L1 returns the L1 norm of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c    The vector L1 norm is defined as:
c
c      R4VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, real A(N), the vector whose L1 norm is desired.
c
c    Output, real R4VEC_NORM_L1, the L1 norm of A.
c
      implicit none

      integer n

      real a(n)
      integer i
      real r4vec_norm_l1

      r4vec_norm_l1 = 0.0E+00
      do i = 1, n
        r4vec_norm_l1 = r4vec_norm_l1 + abs ( a(i) )
      end do

      return
      end
      function r4vec_norm_l2 ( n, a )

c*********************************************************************72
c
cc R4VEC_NORM_L2 returns the L2 norm of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4 values.
c
c    The vector L2 norm is defined as:
c
c      R4VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, real A(N), the vector whose L2 norm is desired.
c
c    Output, real R4VEC_NORM_L2, the L2 norm of A.
c
      implicit none

      integer n

      real a(n)
      integer i
      real r4vec_norm_l2
      real value

      value = 0.0E+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r4vec_norm_l2 = value

      return
      end
      subroutine r4vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R4VEC_NORMAL_01 returns a unit pseudonormal R4VEC.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
c
c    Before calling this routine, the user may call RANDOM_SEED
c    in order to set the seed of the random number generator.
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
c    23 January 2011
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
c    Output, real X(N), a sample of the standard normal PDF.
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
c    Local, real Y, the value saved from the previous call, if
c    SAVED is 1.
c
      implicit none

      integer n

      integer i
      integer m
      integer made
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r(2)
      real r4_uniform_01
      integer saved
      integer seed
      real x(n)
      integer x_hi_index
      integer x_lo_index
      real y

      save made
      save saved
      save y

      data made / 0 /
      data saved / 0 /
      data y / 0.0E+00 /
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
        y = 0.0E+00
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

        r(1) = r4_uniform_01 ( seed )

        if ( r(1) .eq. 0.0E+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R4VEC_NORMAL_01 - Fatal errorc'
          write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
          stop
        end if

        r(2) = r4_uniform_01 ( seed )

        x(x_hi_index) =
     &           sqrt ( -2.0E+00 * log ( r(1) ) )
     &           * cos ( 2.0E+00 * pi * r(2) )
        y =      sqrt ( -2.0E+00 * log ( r(1) ) )
     &           * sin ( 2.0E+00 * pi * r(2) )

        saved = 1

        made = made + 2
c
c  If we require an even number of values, that's easy.
c
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) .eq. 0 ) then

        do i = x_lo_index, x_hi_index, 2

          call r4vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0E+00 * log ( r(1) ) )
     &      * cos ( 2.0E+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0E+00 * log ( r(1) ) )
     &      * sin ( 2.0E+00 * pi * r(2) )

        end do

        made = made + x_hi_index - x_lo_index + 1
c
c  If we require an odd number of values, we generate an even number,
c  and handle the last pair specially, storing one in X(N), and
c  saving the other for later.
c
      else

        do i = x_lo_index, x_hi_index - 1, 2

          call r4vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0E+00 * log ( r(1) ) )
     &      * cos ( 2.0E+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0E+00 * log ( r(1) ) )
     &      * sin ( 2.0E+00 * pi * r(2) )

        end do

        call r4vec_uniform_01 ( 2, seed, r )

        x(n) = sqrt ( -2.0E+00 * log ( r(1) ) )
     &    * cos ( 2.0E+00 * pi * r(1) )

        y = sqrt ( -2.0E+00 * log ( r(2) ) )
     &    * sin ( 2.0E+00 * pi * r(2) )

        saved = 1

        made = made + x_hi_index - x_lo_index + 2

      end if

      return
      end
      subroutine r4vec_permute_cyclic ( n, k, a )

c*********************************************************************72
c
cc R4VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c    For 0 <= K < N, this function cyclically permutes the input vector
c    to have the form
c
c     ( A(K+1), A(K+2), ..., A(N), A(1), ..., A(K) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer K, the increment used.
c
c    Input/output, real A(N), the array to be permuted.
c
      implicit none

      integer n

      real a(n)
      real b(n)
      integer i
      integer i4_modp
      integer i4_wrap
      integer ipk
      integer k

      do i = 1, n
        ipk = i4_wrap ( i + k, 1, n )
        b(i) = a(ipk)
      end do

      do i = 1, n
        a(i) = b(i)
      end  do

      return
      end
      subroutine r4vec_print ( n, a, title )

c*********************************************************************72
c
cc R4VEC_PRINT prints an R4VEC.
c
c  Discussion:
c
c    An R4VEC is an array of R4's.
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
c    Input, real A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
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
      subroutine r4vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc R4VEC_PRINT_PART prints "part" of an R4VEC.
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
c    Input, real A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
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
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '  ........  ..............'
        i = n

        write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(2x,i8,a,1x,g14.6,a)' )
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      subroutine r4vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R4VEC_PRINT_SOME prints "some" of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is an array of R4's.
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
c    Input, real A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
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
      subroutine r4vec_sort_bubble_a ( n, a )

c*********************************************************************72
c
cc R4VEC_SORT_BUBBLE_A ascending sorts an R4VEC using bubble sort.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c    Bubble sort is simple to program, but inefficient.  It should not
c    be used for large arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, real A(N).
c    On input, an unsorted array.
c    On output, the array has been sorted.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer j
      real t

      do i = 1, n - 1
        do j = i + 1, n
          if ( a(j) .lt. a(i) ) then
            t    = a(i)
            a(i) = a(j)
            a(j) = t
          end if
        end do
      end do

      return
      end
      subroutine r4vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM returns a scaled pseudorandom R4VEC.
c
c  Discussion:
c
c    An R4VEC is an array of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2005
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
c    Input, integer M, the number of entries in the vector.
c
c    Input, real A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      real a
      real b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = a + ( b - a ) * real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2006
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine r4vec_variance ( n, a, variance )

c*********************************************************************72
c
cc R4VEC_VARIANCE returns the variance of an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c    The variance of a vector X of length N is defined as
c
c      mean ( X(1:n) ) = sum ( X(1:n) ) / n
c
c      var ( X(1:n) ) = sum ( ( X(1:n) - mean )**2 ) / ( n - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
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
c    Input, real A(N), the vector.
c
c    Output, real VARIANCE, the variance of the vector.
c
      implicit none

      integer n

      real a(n)
      integer i
      real mean
      real variance

      if ( n .lt. 2 ) then

        variance = 0.0E+00

      else

        mean = 0.0E+00
        do i = 1, n
          mean = mean + a(i)
        end do
        mean = mean / real ( n )

        variance = 0.0E+00
        do i = 1, n
          variance = variance + ( a(i) - mean )**2
        end do
        variance = variance / real ( n - 1 )

      end if

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
