      subroutine bicycle_lock ( c, step )

c*********************************************************************72
c
cc BICYCLE_LOCK finds the combination on a typical bicycle lock.
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
c    Input, integer C, the combination, a value between 0 and 999.
c
c    Output, integer STEP, the step on which the combination 
c    was found.  A value of -1 means the combination was not found.
c
      implicit none

      integer a
      integer c
      integer step

      step = -1

      do a = 0, 999

        if ( a .eq. c ) then
          step = a + 1
          go to 10
        end if
      
      end do

10    continue

      return
      end
      subroutine combination_lock ( m, n, c, step )

c*********************************************************************72
c
cc COMBINATION_LOCK determines the combination of a lock.
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
c    Input, integer M, the number of dials.
c
c    Input, integer N, the number of symbols on each dial.
c    We assume the symbols are the integers 0 to N-1.
c
c    Input, integer C(M), the combination.
c
c    Output, integer STEP, the step on which the combination 
c    was found.  A value of -1 means the combination was not found.
c
      implicit none

      integer m

      integer a(m)
      integer c(m)
      integer i
      logical i4vec_eq
      logical more
      integer n
      integer step
c
c  Starting with the guess (0, 0, ... 0),
c  generate every possible combination, in order, and try it.
c
      more = .false.
      do i = 1, m
        a(i) = 0
      end do
      step = 0

10    continue
      
        call combination_next ( m, n, a, more )

        if ( .not. more ) then
          step = -1
          go to 20
        end if

        step = step + 1

        if ( i4vec_eq ( m, a, c ) ) then
          go to 20
        end if
      
      go to 10

20    continue

      return
      end
      subroutine combination_next ( m, base, a, more )

c*********************************************************************72
c
cc COMBINATION_NEXT generates lock combinations in lex order.
c
c  Discussion:
c
c    The vectors are produced in lexical order, starting with
c    (0,0,...,0),
c    (0,0,...,1),
c    ...
c    (BASE-1,BASE-1,...,BASE-1).
c
c  Example:
c
c    M = 2,
c    BASE = 3
c
c    0   0
c    0   1
c    0   2
c    1   0
c    1   1
c    1   2
c    2   0
c    2   1
c    2   2
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
c    Input, integer M, the size of the vectors to be used.
c
c    Input, integer BASE, the base to be used.  BASE = 2 will
c    give vectors of 0's and 1's, for instance.
c
c    Input/output, integer A(M).  The input value of A is
c    not important on the first call.  Thereafter, it should simply be the 
c    output value from the previous call.  The output value is the next vector
c    in the sequence.
c
c    Input/output, logical MORE.  The input value should be FALSE on the first 
c    call, and TRUE on subsequent calls.  The output value will be TRUE as long 
c    as the next vector could be computed, and FALSE once there are no more.
c
      implicit none

      integer m

      integer a(m)
      integer base
      integer i
      logical more

      if ( .not. more ) then

        do i = 1, m
          a(i) = 0
        end do
        more = .true.

      else
          
        do i = m, 1, -1

          a(i) = a(i) + 1

          if ( a(i) .lt. base ) then
            return
          end if

          a(i) = 0

        end do

        more = .false.

      end if

      return
      end
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ) then
        seed = seed - 1
      end if

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
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
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

      i4_uniform = value

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
