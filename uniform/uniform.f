      function c4_uniform_01 ( seed )

c*********************************************************************72
c
cc C4_UNIFORM_01 returns a unit pseudorandom C4.
c
c  Discussion:
c
c    A C4 is a complex single precision value.
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
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
c    Output, complex C4_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      complex c4_uniform_01
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      real pi
      parameter ( pi = 3.1415926E+00 )
      real r
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0E+00 * pi 
     &  * real ( dble ( seed ) * 4.656612875D-10 )

      c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

      return
      end
      subroutine c4mat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
c
c  Discussion:
c
c    A C4MAT is a matrix of complex single precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
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
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      complex c(m,n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      real r
      integer k
      real pi
      parameter ( pi = 3.1415926E+00 )
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4MAT_UNIFORM_01 - Fatal error!'
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

          r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          theta = 2.0D+00 * pi 
     &      * real ( dble ( seed ) * 4.656612875D-10 )

          c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ) )

        end do

      end do

      return
      end
      subroutine c4vec_uniform_01 ( n, seed, c )

c*********************************************************************72
c
cc C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
c
c  Discussion:
c
c    A C4VEC is a vector of complex single precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
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
c    Input, integer N, the number of values to compute.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(N), the pseudorandom complex vector.
c
      implicit none

      integer n

      complex c(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      real r
      integer k
      real pi
      parameter ( pi = 3.141526E+00 )
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        theta = 2.0E+00 * pi * real ( dble ( seed ) * 4.656612875D-10 )

        c(i) = r * cmplx ( cos ( theta ), sin ( theta ) )

      end do

      return
      end
      function c8_uniform_01 ( seed )

c*********************************************************************72
c
cc C8_UNIFORM_01 returns a unit pseudorandom C8.
c
c  Discussion:
c
c    A C8 is a complex double precision value.
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2006
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
c    Output, double complex Z_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      integer seed
      double precision theta
      double complex c8_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( dble ( seed ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

      c8_uniform_01 = r * dcmplx ( dcos ( theta ), dsin ( theta ) )

      return
      end
      subroutine c8mat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of complex double precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
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
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      double complex c(m,n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      double precision r
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seed
      double precision theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C8MAT_UNIFORM_01 - Fatal error!'
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

          r = sqrt ( dble ( seed ) * 4.656612875D-10 )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

          c(i,j) = r * dcmplx ( dcos ( theta ), dsin ( theta ) )

        end do

      end do

      return
      end
      subroutine c8vec_uniform_01 ( n, seed, c )

c*********************************************************************72
c
cc C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
c
c  Discussion:
c
c    A C8VEC is a vector of complex double precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
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
c    Input, integer N, the number of values to compute.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C(N), the pseudorandom complex vector.
c
      implicit none

      integer n

      double complex c(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      integer seed
      double precision theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r = sqrt ( dble ( seed ) * 4.656612875D-10 )

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

        c(i) = r * dcmplx ( dcos ( theta ), dsin ( theta ) )

      end do

      return
      end
      function ch_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc CH_UNIFORM_AB returns a scaled pseudorandom CH.
c
c  Discussion:
c
c    A CH is an alphabetic character value.
c
c    The value is scaled to lie between characters A and B.
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
c    Input, character CLO, CHI, the minimum and maximum acceptable characters.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, character CH_UNIFORM, the randomly chosen character.
c
      implicit none

      character a
      character b
      character ch_uniform_ab
      integer i
      integer ihi
      integer ilo
      real r4_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CH_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      ilo = ichar ( a )
      ihi = ichar ( b )

      i = ilo + int ( r4_uniform_01 ( seed ) * real ( ihi + 1 - ilo ) )
 
      i = max ( i, ilo )
      i = min ( i, ihi )

      ch_uniform_ab = char ( i )

      return
      end
      subroutine congruence ( a, b, c, ierror, x )

c*********************************************************************72
c
cc CONGRUENCE solves a congruence of the form A * X = C ( mod B ).
c
c  Discussion:
c
c    A, B and C are given integers.  The equation is solvable if and only
c    if the greatest common divisor of A and B also divides C.
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
c  Reference:
c
c    Eric Weisstein, editor,
c    CRC Concise Encylopedia of Mathematics,
c    CRC Press, 1998, page 446.
c
c  Parameters:
c
c    Input, integer A, B, C, the coefficients of the Diophantine equation.
c
c    Output, integer IERROR, error flag.
c    0, no error, X was computed.
c    1, A = B = 0, C is nonzero.
c    2, A = 0, B and C nonzero, but C is not a multiple of B.
c    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
c    4, A, B, C nonzero, but GCD of A and B does not divide C.
c    5, algorithm ran out of internal space.
c
c    Output, integer X, the solution of the Diophantine equation.
c    X will be between 0 and B-1.
c
      implicit none

      integer nmax
      parameter ( nmax = 100 )

      integer a
      integer a_copy
      integer a_mag
      integer a_sign
      integer b
      integer b_copy
      integer b_mag
      integer b_sign
      integer c
      integer c_copy
      integer g
      integer i4_gcd
      integer ierror
      integer k
      integer n
      integer q(nmax)
      logical swap
      integer x
      integer y
      integer z
c
c  Defaults for output parameters.
c
      ierror = 0
      x = 0
      y = 0
c
c  Special cases.
c
      if ( a .eq. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .eq. 0 .and. b .eq. 0 .and. c .ne. 0 ) then
        ierror = 1
        x = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .ne. 0 ) then
        x = 0
        if ( mod ( c, b ) .ne. 0 ) then
          ierror = 2
        end if
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c /= 0 ) then
        x = c / a
        if ( mod ( c, a ) .ne. 0 ) then
          ierror = 3
        end if
        return
      else if ( a .ne. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
c       g = i4_gcd ( a, b )
c       x = b / g
        x = 0
        return
      end if
c
c  Handle the "general" case: A, B and C are nonzero.
c
c  Step 1: Compute the GCD of A and B, which must also divide C.
c
      g = i4_gcd ( a, b )

      if ( mod ( c, g ) .ne. 0 ) then
        ierror = 4
        return
      end if

      a_copy = a / g
      b_copy = b / g
      c_copy = c / g
c
c  Step 2: Split A and B into sign and magnitude.
c
      a_mag = abs ( a_copy )
      a_sign = sign ( 1, a_copy )
      b_mag = abs ( b_copy )
      b_sign = sign ( 1, b_copy )
c
c  Another special case, A_MAG = 1 or B_MAG = 1.
c
      if ( a_mag .eq. 1 ) then
        x = a_sign * c_copy
        return
      else if ( b_mag .eq. 1 ) then
        x = 0
        return
      end if
c
c  Step 3: Produce the Euclidean remainder sequence.
c
      if ( b_mag .le. a_mag ) then

        swap = .false.
        q(1) = a_mag
        q(2) = b_mag

      else

        swap = .true.
        q(1) = b_mag
        q(2) = a_mag

      end if

      n = 3

10    continue

        q(n) = mod ( q(n-2), q(n-1) )

        if ( q(n) .eq. 1 ) then
          go to 20
        end if

        n = n + 1

        if ( nmax .lt. n ) then
          ierror = 5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CONGRUENCE - Fatal error!'
          write ( *, '(a)' ) '  Exceeded number of iterations.'
          stop
        end if

      go to 10

20    continue
c
c  Step 4: Go backwards to solve X * A_MAG + Y * B_MAG = 1.
c
      y = 0
      do k = n, 2, -1
        x = y
        y = ( 1 - x * q(k-1) ) / q(k)
      end do
c
c  Step 5: Undo the swapping.
c
      if ( swap ) then
        z = x
        x = y
        y = z
      end if
c
c  Step 6: Apply signs to X and Y so that X * A + Y * B = 1.
c
      x = x * a_sign
c
c  Step 7: Multiply by C, so that X * A + Y * B = C.
c
      x = x * c_copy
c
c  Step 8: Force 0 <= X < B.
c
      x = mod ( x, b )
c
c  Step 9: Force positivity.
c
      if ( x .lt. 0 ) then
        x = x + b
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
      function i4_seed_advance ( seed )

c*********************************************************************72
c
cc I4_SEED_ADVANCE "advances" the seed.
c
c  Discussion:
c
c    This routine implements one step of the recursion
c
c      SEED = ( 16807 * SEED ) mod ( 2^31 - 1 )
c
c    This version of the routine does not check whether the input value of
c    SEED is zero.  If the input value is zero, the output value will be zero.
c
c    If we repeatedly use the output of SEED_ADVANCE as the next input, 
c    and we start with SEED = 12345, then the first few iterates are:
c
c         Input      Output
c          SEED        SEED
c
c         12345   207482415
c     207482415  1790989824
c    1790989824  2035175616
c    2035175616    77048696
c      77048696    24794531
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
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
c    Input, integer SEED, the seed value.
c
c    Output, integer I4_SEED_ADVANCE, the "next" seed.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_seed_advance
      integer k
      integer seed
      integer seed_new

      k = seed / 127773

      seed_new = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed_new < 0 ) then
        seed_new = seed_new + i4_huge
      end if

      i4_seed_advance = seed_new 

      return
      end
      function i4_uniform_0i ( seed )

c*********************************************************************72
c
cc I4_UNIFORM_0I returns a pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    SEED = SEED * (7^5) mod (2^31 - 1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 September 2005
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
c    Input/output, integer SEED, the integer "seed" used to generate
c    the output value.  SEED should not be 0.
c
c    Output, integer I4_UNIFORM_0I, a uniform random value between
c    1 and 2^31-1.
c
c  Local parameters:
c
c    IA = 7^5
c    IB = 2^15
c    IB16 = 2^16
c    IP = 2^31-1
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_uniform_0i
      integer ia
      parameter ( ia = 16807 )
      integer ib15
      parameter ( ib15 = 32768 )
      integer ib16
      parameter ( ib16 = 65536 )
      integer ip
      parameter ( ip = 2147483647 )
      integer iprhi
      integer ixhi
      integer k
      integer leftlo
      integer loxa
      integer seed
c
c  Don't let SEED be 0.
c
      if ( seed .eq. 0 ) then
        seed = i4_huge
      end if
c
c  Get the 15 high order bits of SEED.
c
      ixhi = seed / ib16
c
c  Get the 16 low bits of SEED and form the low product.
c
      loxa = ( seed - ixhi * ib16 ) * ia
c
c  Get the 15 high order bits of the low product.
c
      leftlo = loxa / ib16
c
c  Form the 31 highest bits of the full product.
c
      iprhi = ixhi * ia + leftlo
c
c  Get overflow past the 31st bit of full product.
c
      k = iprhi / ib15
c
c  Assemble all the parts and presubtract IP.  The parentheses are
c  essential.
c
      seed = ( ( ( loxa - leftlo * ib16 ) - ip ) 
     &        + ( iprhi - k * ib15 ) * ib16 ) + k
c
c  Add IP back in if necessary.
c
      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      i4_uniform_0i = seed

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
      subroutine i4mat_uniform_ab ( m, n, a, b, seed, x )

c*********************************************************************72
c
cc I4MAT_UNIFORM_AB returns a scaled pseudorandom I4MAT.
c
c  Discussion:
c
c    An I4MAT is a matrix of integer values.
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
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
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
            seed = seed + i4_huge
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
      subroutine i4vec_uniform_ab ( n, a, b, seed, x )

c*********************************************************************72
c
cc I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
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
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
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
          seed = seed + i4_huge
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
      function i8_huge ( )

c*********************************************************************72
c
cc I8_HUGE returns a "huge" I8.
c
c  Discussion:
c
c    On an IEEE 32 bit machine, I4_HUGE should be 2**63 - 1, and its
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
c    Output, integer ( kind = 8 ) I8_HUGE, a "huge" I8.
c
      implicit none

      integer*8 i8_huge

      i8_huge = 9223372036854775807_8

      return
      end
      function i8_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc I8_UNIFORM_AB returns a scaled pseudorandom I8.
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
c    Output, integer*8 I8_UNIFORM_AB, a number between A and B.
c
      implicit none

      integer*8 a
      integer*8 b
      integer*8 i8_uniform_ab
      double precision r
      double precision r8i8_uniform_01
      integer*8 seed
      integer*8 value

      if ( seed .eq. 0 ) then
       write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I8_UNIFORM_AB - Fatal error!'
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

      i8_uniform_ab = value

      return
      end
      function l_uniform ( seed )

c*********************************************************************72
c
cc L_UNIFORM returns a pseudorandom L.
c
c  Discussion:
c
c    An L is a LOGICAL value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
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
c    Input/output, integer SEED, the "seed" value, which should
c    NOT be 0.  On output, SEED has been updated.
c
c    Output, logical L_UNIFORM, a pseudorandom logical value.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_huge_half
      parameter ( i4_huge_half = 1073741823 )
      integer k
      logical l_uniform
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'L_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      l_uniform = ( i4_huge_half .lt. seed )

      return
      end
      subroutine lcrg_anbn ( a, b, c, n, an, bn )

c*********************************************************************72
c
cc LCRG_ANBN computes the "N-th power" of a linear congruential generator.
c
c  Discussion:
c
c    We are considering a linear congruential random number generator.
c    The LCRG takes as input an integer value called SEED, and returns
c    an updated value of SEED, 
c
c      SEED(out) = ( a * SEED(in) + b ) mod c.
c
c    and an associated pseudorandom real value
c
c      U = SEED(out) / c.
c
c    In most cases, a user is content to call the LCRG repeatedly, with
c    the updating of SEED being taken care of automatically.
c
c    The purpose of this routine is to determine the values of AN and BN
c    that describe the LCRG that is equivalent to N applications of the
c    original LCRG.
c
c    One use for such a facility would be to do random number computations
c    in parallel.  If each of N processors is to compute many random values, 
c    you can guarantee that they work with distinct random values 
c    by starting with a single value of SEED, using the original LCRG to generate
c    the first N-1 "iterates" of SEED, so that you now have N "seed" values,
c    and from now on, applying the N-th power of the LCRG to the seeds.
c
c    If the K-th processor starts from the K-th seed, it will essentially
c    be computing every N-th entry of the original random number sequence,
c    offset by K.  Thus the individual processors will be using a random
c    number stream as good as the original one, and without repeating, and
c    without having to communicate.
c 
c    To evaluate the N-th value of SEED directly, we start by ignoring 
c    the modular arithmetic, and working out the sequence of calculations
c    as follows:
c
c      SEED(0)   =     SEED.
c      SEED(1)   = a * SEED      + b
c      SEED(2)   = a * SEED(1)   + b = a^2 * SEED           + a * b + b
c      SEED(3)   = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
c      ...
c      SEED(N-1) = a * SEED(N-2) + b 
c
c      SEED(N) = a * SEED(N-1) + b = a^N * SEED 
c                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b
c
c    or, using the geometric series,
c
c      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
c              = AN * SEED + BN
c
c    Thus, from any SEED, we can determine the result of N applications of the
c    original LCRG directly if we can solve
c
c      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,
c
c    and evaluate:
c
c      AN = a^N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Barry Wilkinson, Michael Allen,
c    Parallel Programming:
c    Techniques and Applications Using Networked Workstations and Parallel Computers,
c    Prentice Hall, 
c    ISBN: 0-13-140563-2,
c    LC: QA76.642.W54.
c
c  Parameters:
c
c    Input, integer A, the multiplier for the LCRG.
c
c    Input, integer B, the added value for the LCRG.
c
c    Input, integer C, the base for the modular arithmetic.  
c    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
c    required that 0 < C.
c
c    Input, integer N, the "index", or number of times that the
c    LCRG is to be applied.  It is required that 0 <= N.
c
c    Output, integer AN, BN, the multiplier and added value for
c    the LCRG that represent N applications of the original LCRG.
c
      implicit none

      integer a
      integer am1
      integer an
      integer anm1tb
      integer b
      integer bn
      integer c
      integer ierror
      integer n

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of N = ', n
        stop
      end if

      if ( c .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of C = ', c
        stop
      end if

      if ( n .eq. 0 ) then
        an = 1
        bn = 0
      else if ( n .eq. 1 ) then
        an = a
        bn = b
      else
c
c  Compute A^N.
c
        call power_mod ( a, n, c, an )
c
c  Solve 
c    ( a - 1 ) * BN = ( a^N - 1 ) mod B
c  for BN.
c
        am1 = a - 1
        anm1tb = ( an - 1 ) * b

        call congruence ( am1, c, anm1tb, ierror, bn )

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
          write ( *, '(a)' ) '  An error occurred in CONGRUENCE.'
          write ( *, '(a,i8)' ) '  The error code was IERROR = ', ierror
          stop
        end if

      end if

      return
      end
      subroutine lcrg_evaluate ( a, b, c, x, y )

c*********************************************************************72
c
cc LCRG_EVALUATE evaluates an LCRG, y = ( A * x + B ) mod C.
c
c  Discussion:
c
c    This routine cannot be recommended for production use.  Because we want
c    to do modular arithmetic, but the base is not a power of 2, we need to
c    use "double precision" integers to keep accuracy.
c
c    If we knew the base C, we could try to avoid overflow while not changing
c    precision.
c
c    If the base C was a power of 2, we could rely on the usual properties of 
c    integer arithmetic on computers, in which overflow bits, which are always 
c    ignored, don't actually matter.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, the multiplier for the LCRG.
c
c    Input, integer B, the added value for the LCRG.
c
c    Input, integer C, the base for the modular arithmetic.  
c    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
c    required that 0 < C.
c
c    Input, integer X, the value to be processed.
c
c    Output, integer Y, the processed value.
c
      implicit none

      integer a
      integer*8 a8
      integer b
      integer*8 b8
      integer c
      integer*8 c8
      integer x
      integer*8 x8
      integer y
      integer*8 y8
c
c  To avoid roundoff issues, we need to go to "double precision" integers.
c  (Not available on all planets.)
c
      a8 = a
      b8 = b
      c8 = c
      x8 = x

      y8 = mod ( a8 * x8 + b8, c8 )

      y = int ( y8 )

      if ( y < 0 ) then
        y = y + c
      end if

      return
      end
      subroutine lcrg_seed ( a, b, c, n, seed, seed_lcrg )

c*********************************************************************72
c
cc LCRG_SEED computes the N-th seed of a linear congruential generator.
c
c  Discussion:
c
c    We are considering a linear congruential random number generator.
c    The LCRG takes as input an integer value called SEED, and returns
c    an updated value of SEED, 
c
c      SEED(out) = a * SEED(in) + b, mod c.
c
c   and an associated pseudorandom real value
c
c      U = SEED(out) / c.
c
c    In most cases, a user is content to call the LCRG repeatedly, with
c    the updating of SEED being taken care of automatically.
c
c    The purpose of this routine is to determine the value of SEED that
c    would be output after N successive applications of the LCRG.  This
c    allows the user to know, in advance, what the 1000-th value of
c    SEED would be, for instance.  Obviously, one way to do this is to
c    apply the LCRG formula 1,000 times.  However, it is possible to
c    do this in a more direct and efficient way.
c
c    One use for such a facility would be to do random number computations
c    in parallel.  If each processor is to compute 1,000 values, you can
c    guarantee that they work with distinct random values by starting the
c    first processor with SEED, the second with the value of SEED after 
c    1,000 applications of the LCRG, and so on.
c
c    To evaluate the N-th value of SEED directly, we start by ignoring 
c    the modular arithmetic, and working out the sequence of calculations
c    as follows:
c
c      SEED(0) =     SEED.
c      SEED(1) = a * SEED      + b
c      SEED(2) = a * SEED(1)   + b = a^2 * SEED + a * b + b
c      SEED(3) = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
c      ...
c      SEED(N) = a * SEED(N-1) + b = a^N * SEED 
c                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b
c
c    or, using the geometric series,
c
c      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
c
c    Therefore, we can determine SEED(N) directly if we can solve
c
c      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,
c
c    and evaluated:
c
c      AN = a^N
c
c    Using the formula:
c
c      SEED(N) = AN * SEED + BN, mod c
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, the multiplier for the LCRG.
c
c    Input, integer B, the added value for the LCRG.
c
c    Input, integer C, the base for the modular arithmetic.  For 32 bit
c    arithmetic, this is often 2^31 - 1, or 2147483647.  It is required
c    that 0 < C.
c
c    Input, integer N, the "index", or number of times that the LCRG
c    is to be applied.  It is required that 0 <= N.
c
c    Input, integer SEED, the starting value of SEED.  It is customary
c    that 0 < SEED.
c
c    Output, integer SEED_LCRG, the value of SEED that would be output
c    if the LCRG were applied to the starting value N times.
c
      implicit none

      integer a
      integer an
      integer b
      integer bn
      integer c
      integer ierror
      integer n
      integer seed
      integer seed_lcrg
      integer*8 value2

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_SEED - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of N = ', n
        stop
      end if

      if ( c .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_SEED - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of C = ', c
        stop
      end if

      if ( n .eq. 0 ) then
        seed_lcrg = mod ( seed, c )
        if ( seed_lcrg .lt. 0 ) then
          seed_lcrg = seed_lcrg + c
        end if
        return
      end if
c
c  Get A^N.
c
      call power_mod ( a, n, c, an )
c
c  Solve ( a - 1 ) * BN = ( a^N - 1 ) for BN.
c
c  The LCRG I have been investigating uses B = 0, so this code
c  has not been properly tested yet.
c
      call congruence ( a-1, c, ( an - 1 ) * b, ierror, bn )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_SEED - Fatal error!'
        write ( *, '(a)' ) '  An error occurred in CONGRUENCE.'
        write ( *, '(a,i8)' ) '  The error code was IERROR = ', ierror
        stop
      end if
c
c  Set the new SEED.
c
      value2 = an * seed + bn

      value2 = mod ( value2, c )
c
c  Guarantee that the value is positive.
c
      if ( value2 .lt. 0 ) then
        value2 = value2 + c 
      end if

      seed_lcrg = value2

      return
      end
      subroutine lmat_uniform ( m, n, seed, lmat )

c*********************************************************************72
c
cc LMAT_UNIFORM returns a pseudorandom LMAT.
c
c  Discussion:
c
c    An LMAT is a two dimensional array of LOGICAL values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
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
c    Input, integer M, N, the order of the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should
c    NOT be 0.  On output, SEED has been updated.
c
c    Output, logical LMAT(M,N), a pseudorandom logical matrix.
c
      implicit none

      integer m
      integer n

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_huge_half
      parameter ( i4_huge_half = 1073741823 )
      integer i
      integer j
      integer k
      logical lmat(m,n)
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LMAT_UNIFORM - Fatal error!'
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

          lmat(i,j) = ( i4_huge_half .lt. seed )

        end do

      end do

      return
      end
      subroutine lvec_uniform ( n, seed, lvec )

c*********************************************************************72
c
cc LVEC_UNIFORM returns a pseudorandom LVEC.
c
c  Discussion:
c
c    An LVEC is a vector of LOGICAL values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
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
c    Input, integer N, the order of the vector.
c
c    Input/output, integer SEED, the "seed" value, which should
c    NOT be 0.  On output, SEED has been updated.
c
c    Output, logical LVEC(N), a pseudorandom logical vector.
c
      implicit none

      integer n

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_huge_half
      parameter ( i4_huge_half = 1073741823 )
      integer i
      integer k
      logical lvec(n)
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LVEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        lvec(i) = ( i4_huge_half .lt. seed )

      end do
 
      return
      end
      subroutine power_mod ( a, n, m, x )

c*********************************************************************72
c
cc POWER_MOD computes mod ( A^N, M ).
c
c  Discussion:
c
c    Some programming tricks are used to speed up the computation, and to
c    allow computations in which the value A**N is much too large to
c    store in an integer word.
c
c    First, for efficiency, the power A**N is computed by determining
c    the binary expansion of N, then computing A, A^2, A^4, and so on
c    by repeated squaring, and multiplying only those factors that
c    contribute to A**N.
c
c    Secondly, the intermediate products are immediately "mod'ed", which
c    keeps them small.
c
c    For instance, to compute mod ( A^13, 11 ), we essentially compute
c
c       13 = 1 + 4 + 8
c
c       A^13 = A * A^4 * A^8
c
c       mod ( A^13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).
c
c    Fermat's little theorem says that if P is prime, and A is not divisible
c    by P, then ( A^(P-1) - 1 ) is divisible by P.
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
c    Input, integer A, the base of the expression to be tested.
c    0 <= A is required.
c
c    Input, integer N, the power to which the base is raised.
c    0 <= N is required.
c
c    Input, integer M, the divisor against which the expression is tested.
c    0 < M is required.
c
c    Output, integer X, the remainder when A^N is divided by M.
c    If any input quantity is unacceptable, then the nonsensical value
c    X = -1 is returned.
c
      implicit none

      integer a
      integer*8 a_square2
      integer d
      integer m
      integer*8 m2
      integer n
      integer ncopy
      integer x
      integer*8 x2

      if ( a .lt. 0 ) then
        x = -1
        return
      end if

      if ( m .le. 0 ) then
        x = -1
        return
      end if

      if ( n .lt. 0 ) then
        x = -1
        return
      end if
c
c  A_SQUARE contains the successive squares of A.
c
      a_square2 = a
      x2 = 1
      m2 = m

      ncopy = n

10    continue

      if ( 0 .lt. ncopy ) then

        d = mod ( ncopy, 2 )

        if ( d .eq. 1 ) then
          x2 = mod ( x2 * a_square2, m2 )
        end if

        a_square2 = mod ( a_square2 * a_square2, m2 )
        ncopy = ( ncopy - d ) / 2
        go to 10

      end if
c
c  Fix up X so that it is nonnegative.
c
20    continue

      if ( x2 .lt. 0 ) then
        x2 = x2 + m2
        go to 20
      end if

      x = int ( x2 )

      return
      end
      function r4_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc R4_UNIFORM_AB returns a scaled pseudorandom R4.
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
c    Output, real R4_UNIFORM_AB, a number strictly between A and B.
c
      implicit none

      real a
      real b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r4_uniform_ab

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r4_uniform_ab = a + ( b - a ) 
     &  * real ( dble ( seed ) * 4.656612875D-10 )

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
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r4_uniform_01 = seed / ( 2^31 - 1 )
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
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
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
        seed = seed + i4_huge
      end if

      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
      subroutine r4mat_uniform_ab ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R4MAT_UNIFORM_AB returns a scaled pseudorandom R4MAT.
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
        write ( *, '(a)' ) 'R4MAT_UNIFORM_AB - Fatal error!'
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
      subroutine r4mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      real r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4MAT_UNIFORM_01 - Fatal error!'
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

          r(i,j) = real ( seed ) * 4.656612875E-10

        end do
      end do

      return
      end
      subroutine r4vec_uniform_ab ( n, a, b, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_AB returns a scaled pseudorandom R4VEC.
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
        write ( *, '(a)' ) 'R4VEC_UNIFORM_AB - Fatal error!'
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
      function r8_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM_AB returns a scaled pseudorandom R8.
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
c    06 January 2006
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
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_AB, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_ab
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r8_uniform_ab = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

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

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8col_uniform_ab ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8COL_UNIFORM_AB fills an R8COL with scaled pseudorandom numbers.
c
c  Discussion:
c
c    An R8COL is an array of R8 values, regarded as a set of column vectors.
c
c    The user specifies a minimum and maximum value for each row.
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
c    Input, integer M, N, the number of rows and columns in
c    the array.
c
c    Input, double precision A(M), B(M), the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which
c    should NOT be 0.  On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      double precision a(m)
      double precision b(m)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = a(i) 
     &      + ( b(i) - a(i) ) * dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      function r8i8_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc R8I8_UNIFORM_AB returns a scaled pseudorandom R8 using an I8 seed.
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
c    Output, double precision R8I8_UNIFORM_AB, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer*8 i8_huge
      integer*8 k
      double precision r8i8_uniform_ab
      integer*8 seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8I8_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i8_huge ( )
      end if

      r8i8_uniform_ab = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

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
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
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

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8mat_uniform_ab ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
c
c  Discussion:
c
c    A <= R(I,J) <= B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2005
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
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      double precision a
      double precision b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_AB - Fatal error!'
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

          r(i,j) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8mat_uniform_abvec ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.
c
c  Discussion:
c
c    A(I) <= R(I,J) <= B(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 October 2012
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
c    Input, double precision A(M), B(M), the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      double precision a(m)
      double precision b(m)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_ABVEC - Fatal error!'
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

          r(i,j) = a(i) + ( b(i) - a(i) ) * dble ( seed ) 
     &      * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8row_uniform_ab ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8ROW_UNIFORM_AB fills an R8ROW with scaled pseudorandom numbers.
c
c  Discussion:
c
c    An R8ROWL is an array of R8 values, regarded as a set of row vectors.
c
c    The user specifies a minimum and maximum value for each column.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 2012
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
c    Input, integer M, N, the number of rows and columns in
c    the array.
c
c    Input, double precision A(N), B(N), the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which
c    should NOT be 0.  On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      do i = 1, m

        do j = 1, n

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = a(j) 
     &      + ( b(j) - a(j) ) * dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2004
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
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_uniform_ab ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
c
c  Discussion:
c
c    Each dimension ranges from A to B.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_uniform_abvec ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.
c
c  Discussion:
c
c    Dimension I ranges from A(I) to B(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2012
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
c    Input, double precision A(N), B(N), the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_ABVEC - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = a(i) + ( b(i) - a(i) ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_uniform_unit ( m, seed, w )

c*********************************************************************72
c
cc R8VEC_UNIFORM_UNIT generates a uniformly random unit vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision W(M), a random direction vector,
c    with unit norm.
c
      implicit none

      integer m

      integer i
      double precision norm
      double precision r8vec_norm_l2
      integer seed
      double precision w(m)
c
c  Get N values from a standard normal distribution.
c
      call r8vec_normal_01 ( m, seed, w )
c
c  Compute the length of the vector.
c
      norm = r8vec_norm_l2 ( m, w )
c
c  Normalize the vector.
c
      do i = 1, m
        w(i) = w(i) / norm
      end do

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
