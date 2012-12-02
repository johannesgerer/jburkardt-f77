      function c4_normal_01 ( seed )

c*********************************************************************72
c
cc C4_NORMAL_01 returns a unit pseudonormal C4.
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
c    Output, complex C4_NORMAL_01, a unit pseudonormal value.
c
      implicit none

      complex c4_normal_01
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r4_uniform_01
      integer seed
      real v1
      real v2
      real x_c
      real x_r

      v1 = r4_uniform_01 ( seed )
      v2 = r4_uniform_01 ( seed )

      x_r = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * pi * v2 )
      x_c = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * pi * v2 )

      c4_normal_01 = cmplx ( x_r, x_c )

      return
      end
      function c8_normal_01 ( seed )

c*********************************************************************72
c
cc C8_NORMAL_01 returns a unit pseudonormal C8.
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
c    Output, double complex C8_NORMAL_01, a sample of the PDF.
c
      implicit none

      double complex c8_normal_01
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision v1
      double precision v2
      double precision x_c
      double precision x_r

      v1 = r8_uniform_01 ( seed )
      v2 = r8_uniform_01 ( seed )

      x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * pi * v2 )
      x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * pi * v2 )

      c8_normal_01 = cmplx ( x_r, x_c )

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
      function i4_normal ( a, b, seed )

c*********************************************************************72
c
cc I4_NORMAL returns a scaled pseudonormal I4.
c
c  Discussion:
c
c    The normal probability distribution function (PDF) is sampled,
c    with mean A and standard deviation B.
c
c    The result is then rounded to the nearest integer.
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
c    Output, integer I4_NORMAL, a sample of the normal PDF.
c
      implicit none

      real a
      real b
      integer i4_normal
      real r4_uniform_01
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r1
      real r2
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
          write ( *, '(a)' ) 'I4_NORMAL - Fatal error!'
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

      i4_normal = nint ( a + b * x )

      return
      end
      function i8_normal ( a, b, seed )

c*********************************************************************72
c
cc I8_NORMAL returns a scaled pseudonormal I8.
c
c  Discussion:
c
c    The normal probability distribution function (PDF) is sampled,
c    with mean A and standard deviation B.
c
c    The result is then rounded to the nearest integer.
c
c    I changed my mind, and backed down from using an integer*8 as
c    the seed to r8_uniform_01, so this routine won't work until
c    I decide how to redo it.
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
c    Input, double precision A, the mean of the PDF.
c
c    Input, double precision B, the standard deviation of the PDF.
c
c    Input/output, integer*8 SEED, a seed for the
c    random number generator.
c
c    Output, integer*8 I8_NORMAL, a sample of the normal PDF.
c
      implicit none

      double precision a
      double precision b
      integer*8 i8_normal
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_uniform_01
      integer*8 seed
      integer*8 seed2
      integer*8 used
      double precision x
      double precision y

      save seed2
      save used
      save y

      data seed2 / 0 /
      data used / 0 /
      data y / 0.0D+00 /
c
c  On odd numbered calls, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) == 0 ) then

        r1 = r8_uniform_01 ( seed )

        if ( r1 == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I8_NORMAL - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed
        r2 = r8_uniform_01 ( seed2 )

        x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
c
c  On odd calls, return the second normal and its corresponding seed.
c
      else

        seed = seed2
        x = y

      end if

      used = used + 1

      i8_normal = nint ( a + b * x )

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

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

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

        if ( r(1) .eq. 0.0D+00 ) then
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
      function r8_normal ( a, b, seed )

c*********************************************************************72
c
cc R8_NORMAL returns a scaled pseudonormal R8.
c
c  Discussion:
c
c    The normal probability distribution function (PDF) is sampled,
c    with mean A and standard deviation B.
c
c    I changed my mind, and backed down from using an integer*8 as
c    the seed to r8_uniform_01, so this routine won't work until
c    I decide how to redo it.
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
c    Input, double precision A, the mean of the PDF.
c
c    Input, double precision B, the standard deviation of the PDF.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL, a sample of the normal PDF.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal
      double precision r8_uniform_01
      integer seed
      integer seed2
      integer used
      double precision x
      double precision y

      save seed2
      save used
      save y

      data seed2 / 0 /
      data used / 0 /
      data y / 0.0D+00 /
c
c  On odd numbered calls, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
         end if

        seed2 = seed
        r2 = r8_uniform_01 ( seed2 )

        x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
c
c  On odd calls, return the second normal and its corresponding seed.
c
      else

        seed = seed2
        x = y

      end if

      used = used + 1

      r8_normal = a + b * x

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

      integer k
      double precision r8_uniform_01
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_normal ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8MAT_NORMAL returns a scaled pseudonormal R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input, double precision A, B, the mean and standard deviation.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudonormal values.
c
      implicit none

      integer m
      integer n

      double precision a
      double precision b
      integer seed
      double precision r(m,n)

      call r8vec_normal ( m * n, a, b, seed, r )

      return
      end
      subroutine r8mat_normal_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudonormal values.
c
      implicit none

      integer m
      integer n

      integer seed
      double precision r(m,n)

      call r8vec_normal_01 ( m * n, seed, r )

      return
      end
      subroutine r8vec_normal ( n, a, b, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL returns a scaled pseudonormal R8VEC.
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
c    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
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

      double precision a
      double precision b
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
          write ( *, '(a)' ) 'R8VEC_NORMAL - Fatal errorc'
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

      do i = 1, n
        x(i) = a + b * x(i)
      end do

      return
      end
      subroutine r8vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
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
