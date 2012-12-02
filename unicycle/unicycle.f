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
      subroutine perm_check ( n, p )

c*********************************************************************72
c
cc PERM_CHECK checks that a vector represents a permutation.
c
c  Discussion:
c
c    The routine verifies that each of the integers from 1
c    to N occurs among the N entries of the permutation.
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
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the permutation, in standard index form.
c
      implicit none

      integer n

      integer ierror
      integer ifind
      integer iseek
      integer p(n)

      ierror = 0

      do iseek = 1, n

        ierror = iseek

        do ifind = 1, n
          if ( p(ifind) .eq. iseek ) then
            ierror = 0
            exit
          end if
        end do

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a proper permutation.  In particular,'
          write ( *, '(a,i8)' ) '  it is missing the value ', ierror
          stop
        end if

      end do

      return
      end
      subroutine perm_inverse ( n, p, pinv  )

c*********************************************************************72
c
cc PERM_INVERSE computes the inverse of a permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer PINV(N), the inverse permutation.
c
      implicit none

      integer n

      integer i
      integer p(n)
      integer pinv(n)
c
c  Check.
c
      call perm_check ( n, p )

      do i = 1, n
        pinv(p(i)) = i
      end do

      return
      end
      function perm_is_unicycle ( n, p )

c*********************************************************************72
c
cc PERM_IS_UNICYCLE is TRUE if a permutation is a unicycle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects in the permutation.
c
c    Input, integer P(N), the permutation.
c
c    Output, logical PERM_IS_UNICYCLE, is TRUE if the permutation is a unicycle.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(n)
      logical perm_is_unicycle

      perm_is_unicycle = .false.

      call perm_check ( n, p )
c
c  From 1, you must be able to take N-1 steps without reaching 1...
c
      i = 1
      do j = 1, n - 1
        i = p(i)
        if ( i .eq. 1 ) then
          return
        end if
      end do
c
c  ...and the N-th step must reach 1.
c
      i = p(i)
      if ( i .eq. 1 ) then
        perm_is_unicycle = .true.
      end if

      return
      end
      subroutine perm_lex_next ( n, p, rank )

c*********************************************************************72
c
cc PERM_LEX_NEXT computes the lexicographic permutation successor.
c
c  Example:
c
c    RANK  Permutation
c
c       0  1 2 3 4
c       1  1 2 4 3
c       2  1 3 2 4
c       3  1 3 4 2
c       4  1 4 2 3
c       5  1 4 3 2
c       6  2 1 3 4
c       ...
c      23  4 3 2 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input/output, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is -1.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(n)
      integer rank
      integer temp
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        call i4vec_indicator ( n, p )
        rank = 0
        return
      end if
c
c  Check.
c
      call perm_check ( n, p )
c
c  Seek I, the highest index for which the next element is bigger.
c
      i = n - 1

10    continue

        if ( i .le. 0 ) then
          go to 20
        end if

        if ( p(i) .le. p(i+1) ) then
          go to 20
        end if

        i = i - 1

      go to 10

20    continue
c
c  If no I could be found, then we have reach the final permutation,
c  N, N-1, ..., 2, 1.  Time to start over again.
c
      if ( i .eq. 0 ) then
        call i4vec_indicator ( n, p )
        rank = -1
      else
c
c  Otherwise, look for the the highest index after I whose element
c  is bigger than I's.  We know that I+1 is one such value, so the
c  loop will never fail.
c
        j = n

30      continue

        if ( p(j) .lt. p(i) ) then
          j = j - 1
          go to 30
        end if
c
c  Interchange elements I and J.
c
        temp = p(i)
        p(i) = p(j)
        p(j) = temp
c
c  Reverse the elements from I+1 to N.
c
        call i4vec_reverse ( n - i, p(i+1:n) )
        rank = rank + 1

      end if

      return
      end
      subroutine perm_lex_rank ( n, p, rank )

c*********************************************************************72
c
cc PERM_LEX_RANK computes the lexicographic rank of a permutation.
c
c  Discussion:
c
c    The original code altered the input permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer RANK, the rank of the permutation.
c
      implicit none

      integer n

      integer i
      integer i4_factorial
      integer j
      integer p(n)
      integer pcopy(n)
      integer rank
c
c  Check.
c
      call perm_check ( n, p )

      rank = 0
      do i = 1, n
        pcopy(i) = p(i)
      end do

      do j = 1, n

        rank = rank + ( pcopy(j) - 1 ) * i4_factorial ( n - j )

        do i = j + 1, n
          if ( pcopy(j) .lt. pcopy(i) ) then
            pcopy(i) = pcopy(i) - 1
          end if
        end do

      end do

      return
      end
      subroutine perm_lex_unrank ( n, rank, p )

c*********************************************************************72
c
cc PERM_LEX_UNRANK computes the permutation of given lexicographic rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer RANK, the rank of the permutation.
c
c    Output, integer P(N), describes the permutation.
c
      implicit none

      integer n

      integer d
      integer i
      integer i4_factorial
      integer j
      integer nperm
      integer p(n)
      integer rank
      integer rank_copy
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call perm_enum ( n, nperm )

      if ( rank .lt. 0 .or. nperm .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      p(n) = 1

      do j = 1, n - 1

        d = mod ( rank_copy, i4_factorial ( j + 1 ) ) 
     &    / i4_factorial ( j )
        rank_copy = rank_copy - d * i4_factorial ( j )
        p(n-j) = d + 1

        do i = n - j + 1, n

          if ( d .lt. p(i) ) then
            p(i) = p(i) + 1
          end if

        end do

      end do

      return
      end
      subroutine perm_print ( n, p, title )

c*********************************************************************72
c
cc PERM_PRINT prints a permutation.
c
c  Example:
c
c    Input:
c
c      P = 7 2 4 1 5 3 6
c
c    Printed output:
c
c      "This is the permutation:"
c
c      1 2 3 4 5 6 7
c      7 2 4 1 5 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P(N), the permutation, in standard index form.
c
c    Input, character ( len = * ) TITLE, a title.
c    If no title is supplied, then only the permutation is printed.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer p(n)
      character * ( * ) title

      if ( len_trim ( title ) .ne. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(a)' ) ' '
          write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      else

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      end if

      return
      end
      subroutine perm_random ( n, seed, p )

c*********************************************************************72
c
cc PERM_RANDOM selects a random permutation of N objects.
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
c    12 May 2002
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
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, integer P(N), a permutation of ( 1, 2, ..., N ),
c    in standard index form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer p(n)
      integer seed
      integer t

      call i4vec_indicator ( n, p )

      do i = 1, n - 1

        j = i4_uniform ( i, n, seed )

        t    = p(i)
        p(i) = p(j)
        p(j) = t

      end do

      return
      end
      subroutine perm_enum ( n, nperm )

c*********************************************************************72
c
cc PERM_ENUM enumerates the permutations on N digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be nonnegative.
c
c    Output, integer NPERM, the number of distinct elements.
c
      implicit none

      integer i4_factorial
      integer n
      integer nperm

      nperm = i4_factorial ( n )

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
      subroutine unicycle_check ( n, p )

c*********************************************************************72
c
cc UNICYCLE_CHECK checks that a vector represents a unicycle.
c
c  Discussion:
c
c    A unicycle is a permutation with a single cycle.  This might be called
c    a cyclic permutation, except that that name is used with at least two
c    other meanings.  So, to be clear, a unicycle is a permutation of N
c    objects in which each object is returned to itself precisely after 
c    N applications of the permutation.
c
c    This routine verifies that each of the integers from 1
c    to N occurs among the N entries of the permutation.
c
c    Any permutation of the integers 1 to N describes a unicycle.
c    The permutation ( 3, 4, 2, 1 ) indicates that the unicycle
c    sends 3 to 4, 4 to 2, 2 to 1 and 1 to 3.  This is the sequential
c    description of a unicycle.
c
c    The standard sequence "rotates" the permutation so that it begins
c    with 1.  The above sequence becomes a standard sequence when
c    written as ( 1, 3, 4, 2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the unicycle sequence vector
c
      implicit none

      integer n

      integer ierror
      integer ifind
      integer iseek
      integer p(n)

      ierror = 0

      do iseek = 1, n

        ierror = iseek

        do ifind = 1, n
          if ( p(ifind) .eq. iseek ) then
            ierror = 0
            exit
          end if
        end do

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'UNICYCLE_CHECK - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a unicycle.  In particular, the'
          write ( *, '(a,i8)' ) '  array is missing the value ', ierror
          stop
        end if

      end do

      return
      end
      subroutine unicycle_enum ( n, num )

c*********************************************************************72
c
cc UNICYCLE_ENUM enumerates the unicycles.
c
c  Discussion:
c
c    Each standard sequence corresponds to a unique unicycle.  Since the
c    first element of a standard sequence is always 1, the number of standard
c    sequences, and hence the number of unicycles, is (n-1)!.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicyle.
c
c    Output, integer NUM, the number of unicycles.
c
      implicit none

      integer i4_factorial
      integer n
      integer num

      num = i4_factorial ( n - 1 )

      return
      end
      subroutine unicycle_index ( n, u, u_index )

c*********************************************************************72
c
cc UNICYCLE_INDEX returns the index form of a unicycle.
c
c  Example:
c
c    N = 4
c
c    U       = 1 3 4 2
c    U_INDEX = 3 1 4 2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycles.
c
c    Input, integer U(N), the unicycle sequence vector.
c
c    Output, integer U_INDEX(N), the unicycle index vector.
c
      implicit none

      integer n

      integer i
      integer i4_wrap
      integer ip1
      integer u(n)
      integer u_index(n)

      do i = 1, n
        ip1 = i4_wrap ( i + 1, 1, n )
        u_index(u(i)) = u(ip1)
      end do
      
      return
      end
      subroutine unicycle_index_print ( n, u_index, title )

c*********************************************************************72
c
cc UNICYCLE_INDEX_PRINT prints a unicycle given in index form.
c
c  Example:
c
c    Input:
c
c      U_INDEX = 7 1 4 5 2 3 6
c
c    Printed output:
c
c      1 2 3 4 5 6 7
c      7 1 4 5 2 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycle.
c
c    Input, integer U_INDEX(N), the unicycle index vector.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer u_index(n)
      character * ( * ) title

      if ( 0 .lt. len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
      end if

      do ilo = 1, n, inc
        ihi = min ( n, ilo + inc - 1 )
        write ( *, '(a)' ) ' '
        write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
        write ( *, '(2x,20i4)' ) u_index(ilo:ihi)
      end do

      return
      end
      subroutine unicycle_index_to_sequence ( n, u_index, u )

c*********************************************************************72
c
cc UNICYCLE_INDEX_TO_SEQUENCE converts a unicycle from index to sequence form.
c
c  Example:
c
c    N = 4
c
c    U_INDEX = 3 1 4 2
c    U       = 1 3 4 2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycles.
c
c    Output, integer U_INDEX(N), the unicycle index vector.
c
c    Input, integer U(N), the unicycle sequence vector.
c
      implicit none

      integer n

      integer i
      integer j
      integer u(n)
      integer u_index(n)

      u(1) = 1
      i = 1

      do j = 2, n

        i = u_index(i)
        u(j) = i

        if ( i .eq. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'UNICYCLE_INDEX_TO_SEQUENCE - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The index vector does not represent a unicycle.'
          write ( *, '(a,i4,a,i4,a)' ) 
     &      '  On step ', j, ' u_index(', i, ') = 1.'
          stop
        end if

      end do
      
      return
      end
      subroutine unicycle_inverse ( n, u, u_inverse )

c*********************************************************************72
c
cc UNICYCLE_INVERSE returns the inverse of a unicycle.
c
c  Example:
c
c    N = 4
c
c    U         = 1 3 4 2
c    U_INVERSE = 1 2 4 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycles.
c
c    Input, integer U(N), the unicycle sequence vector.
c
c    Output, integer U_INVERSE(N), the inverse unicycle.
c
      implicit none

      integer n

      integer i
      integer u(n)
      integer u_inverse(n)

      u_inverse(1) = 1
      do i = 2, n
        u_inverse(i) = u(n+2-i)
      end do

      return
      end
      subroutine unicycle_next ( n, u, rank )

c*********************************************************************72
c
cc UNICYCLE_NEXT generates unicycles in lexical order, one at a time.
c
c  Example:
c
c    N = 4
c
c    1   1 2 3 4
c    2   1 2 4 3
c    3   1 3 2 4
c    4   1 3 4 2
c    5   1 4 2 3
c    6   1 4 3 2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycles.
c
c    Input/output, integer U(N); on first call with MORE = FALSE,
c    this value is not used.  Otherwise, the input value is the previous
c    unicycle.  The output value is the next unicycle.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is -1.
c
      implicit none

      integer n

      integer i
      integer p(n-1)
      integer rank
      integer u(n)

      if ( rank .eq. -1 ) then
        u(1) = 1
      else
        do i = 1, n - 1
          p(i) = u(i+1) - 1
        end do
      end if

      call perm_lex_next ( n - 1, p, rank )
     
      do i = 2, n
        u(i) = p(i-1) + 1
      end do

      return
      end
      subroutine unicycle_print ( n, u, title )

c*********************************************************************72
c
cc UNICYCLE_PRINT prints a unicycle given in sequence form.
c
c  Example:
c
c    Input:
c
c      U = 7 1 4 5 2 3 6
c
c    Printed output:
c
c      7 1 4 5 2 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the unicycle.
c
c    Input, integer U(N), the unicycle sequence vector.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer u(n)
      character * ( * ) title

      if ( 0 .lt. len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
        write ( *, '(a)' ) ' '
      end if

      do ilo = 1, n, inc
        ihi = min ( n, ilo + inc - 1 )
        write ( *, '(2x,20i4)' ) ( u(i), i = ilo, ihi )
      end do

      return
      end
      subroutine unicycle_random ( n, seed, u )

c*********************************************************************72
c
cc UNICYCLE_RANDOM selects a random unicycle of N objects.
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
c    22 August 2012
c
c  Author:
c
c    John Burkardt.
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
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, integer U(N), a unicycle in sequence form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer u(n)
      integer seed
      integer t

      call i4vec_indicator ( n, u )

      do i = 2, n
        j = i4_uniform ( i, n, seed )
        t = u(i)
        u(i) = u(j)
        u(j) = t
      end do

      return
      end
      subroutine unicycle_rank ( n, u, rank )

c*********************************************************************72
c
cc UNICYCLE_RANK computes the rank of a unicycle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt.
c
c  Parameters:
c
c    Input, integer N, the order of the unicycle.
c
c    Input, integer U(N), a unicycle in sequence form.
c
c    Output, integer RANK, the rank of the unicycle.
c
      implicit none

      integer n

      integer i
      integer p(n-1)
      integer rank
      integer u(n)

      do i = 1, n - 1
        p(i) = u(i+1) - 1
      end do

      call perm_lex_rank ( n - 1, p, rank )

      return
      end
      subroutine unicycle_unrank ( n, rank, u )

c*********************************************************************72
c
cc UNICYCLE_UNRANK "unranks" a unicycle.
c
c  Discussion:
c
c    That is, given a rank, it computes the corresponding unicycle.
c
c    The value of the rank should be between 0 and (N-1)!-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472,
c    LC: QA164.S79.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set.
c
c    Input, integer RANK, the desired rank of the permutation.
c
c    Output, integer U(N), the unicycle.
c
      implicit none

      integer n

      integer i
      integer p(n-1)
      integer rank
      integer u(n)

      call perm_lex_unrank ( n - 1, rank, p )

      u(1) = 1
      do i = 2, n
        u(i) = p(i-1) + 1
      end do

      return
      end
