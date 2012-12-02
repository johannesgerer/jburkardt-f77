      program main

c*********************************************************************72
c
cc MAIN is the main program for MGS_PRB.
c
c  Discussion:
c
c    MGS_TEST gives some test data to the MGS function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4_uniform
      integer m
      integer n
      integer seed
      integer test

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MGS_PRB:'
      write ( *, '(a)' ) '  Test cases for MGS.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, 4

        m = i4_uniform ( 1, 20, seed )
        n = i4_uniform ( 1, 20, seed )

        call test_case ( m, n, seed )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MGS_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test_case ( m, n, seed )

c*********************************************************************72
c
cc TEST_CASE sets up a test case.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      real a(m,n)
      real a_qr(m,n)
      real a_save(m,n)
      real a_hi
      real a_lo
      real diff
      integer i
      integer j
      integer k
      real q(m,n)
      real r(n,n)
      integer seed
      integer test
 
      a_lo = -10.0
      a_hi = +10.0

      call r4mat_uniform ( m, n, a_lo, a_hi, seed, a )

      do j = 1, n
        do i = 1, m
          a_save(i,j) = a(i,j)
        end do
      end do

      do j = 1, n
        do i = 1, m
          q(i,j) = 0.0
        end do
      end do

      do j = 1, n
        do i = 1, n
          r(i,j) = 0.0
        end do
      end do

      call mgs ( m, n, a, r, q )

      do j = 1, n
        do i = 1, m
          a_qr(i,j) = 0.0
          do k = 1, n
            a_qr(i,j) = a_qr(i,j) + q(i,k) * r(k,j)
          end do
        end do
      end do

      diff = 0.0
      do j = 1, n
        do i = 1, m
          diff = diff + ( a_save(i,j) - a_qr(i,j) )**2
        end do
      end do

      diff = sqrt ( diff ) / sqrt ( real ( m * n ) )

      write ( *, * ) m, n, diff

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
      subroutine r4mat_uniform ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R4MAT_UNIFORM returns a scaled pseudorandom R4MAT.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input, real A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      real a
      real b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      real r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4MAT_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = a + ( b - a ) * real ( seed ) * 4.656612875E-10

        end do
      end do

      return
      end
