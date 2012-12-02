      subroutine comp_next ( n, k, a, more, h, t )

c*********************************************************************72
c
cc COMP_NEXT computes the compositions of the integer N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K parts is an ordered sequence
c    of K nonnegative integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The routine computes one composition on each call until there are no more.
c    For instance, one composition of 6 into 3 parts is
c    3+2+1, another would be 6+0+0.
c
c    On the first call to this routine, set MORE = FALSE.  The routine
c    will compute the first element in the sequence of compositions, and
c    return it, as well as setting MORE = TRUE.  If more compositions
c    are desired, call again, and again.  Each time, the routine will
c    return with a new composition.
c
c    However, when the LAST composition in the sequence is computed
c    and returned, the routine will reset MORE to FALSE, signaling that
c    the end of the sequence has been reached.
c
c  Example:
c
c    The 28 compositions of 6 into three parts are:
c
c      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
c      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
c      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
c      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
c      0 3 3,  0 2 4,  0 1 5,  0 0 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2007
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
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.
c
c    Input/output, integer A(K), the parts of the composition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
c    Input/output, integer H, T, values used by the program.
c    The user should NOT set or alter these quantities.
c
      implicit none

      integer k

      integer a(k)
      integer h
      integer i
      logical more
      integer n
      integer t
c
c  The first computation.
c
      if ( .not. more ) then

        t = n
        h = 0
        a(1) = n
        do i = 2, k
          a(i) = 0
        end do
c
c  The next computation.
c
      else

        if ( 1 < t ) then
          h = 0
        end if

        h = h + 1
        t = a(h)
        a(h) = 0
        a(1) = t - 1
        a(h+1) = a(h+1) + 1

      end if
c
c  This is the last element of the sequence if all the
c  items are in the last slot.
c
      more = ( a(k) .ne. n )

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
      subroutine gm_rule_set ( rule, dim_num, point_num, w, x )

c*********************************************************************72
c
cc GM_RULE_SET sets a Grundmann-Moeller rule.
c
c  Discussion:
c
c    This is a revised version of the calculation which seeks to compute
c    the value of the weight in a cautious way that avoids intermediate
c    overflow.  Thanks to John Peterson for pointing out the problem on
c    26 June 2008.
c
c    This rule returns weights and abscissas of a Grundmann-Moeller 
c    quadrature rule for the DIM_NUM-dimensional unit simplex.
c
c    The dimension POINT_NUM can be determined by calling GM_RULE_SIZE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex 
c    by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c    0 .le. RULE.
c
c    Input, integer DIM_NUM, the spatial dimension.
c    1 .le. DIM_NUM.
c
c    Input, integer POINT_NUM, the number of points in the rule.
c
c    Output, double precision W(POINT_NUM), the weights.
c
c    Output, double precision X(DIM_NUM,POINT_NUM), the abscissas.
c
      implicit none

      integer dim_num
      integer point_num

      integer beta(dim_num+1)
      integer beta_sum
      integer d
      integer h
      integer i
      integer i1
      integer j
      integer k
      logical more
      integer n
      integer one_pm
      integer rule
      integer s
      integer t
      double precision w(point_num)
      double precision weight
      double precision x(dim_num,point_num)

      s = rule
      d = 2 * s + 1
      k = 0
      n = dim_num
      one_pm = 1

      do i = 0, s

        weight = dble ( one_pm )

        do j = 1, max ( n, d, d + n - i )

          if ( j .le. n ) then
            weight = weight * dble ( j )
          end if
          if ( j .le. d ) then
            weight = weight * dble ( d + n - 2 * i )
          end if
          if ( j .le. 2 * s ) then
            weight = weight / 2.0D+00
          end if
          if ( j .le. i ) then
            weight = weight / dble ( j )
          end if
          if ( j .le. d + n - i ) then
            weight = weight / dble ( j )
          end if

        end do

        one_pm = - one_pm

        beta_sum = s - i
        more = .false.
        h = 0
        t = 0

10      continue

          call comp_next ( beta_sum, dim_num + 1, beta, more, h, t )

          k = k + 1

          w(k) = weight

          do i1 = 1, dim_num
            x(i1,k) = dble ( 2 * beta(i1+1) + 1 ) 
     &              / dble ( d + n - 2 * i )
          end do

          if ( .not. more ) then
            go to 20
          end if

        go to 10

20      continue
      
      end do

      return
      end
      subroutine gm_rule_size ( rule, dim_num, point_num )

c*********************************************************************72
c
cc GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
c
c  Discussion:
c
c    This rule returns the value of POINT_NUM, the number of points associated
c    with a GM rule of given index.
c
c    After calling this rule, the user can use the value of POINT_NUM to
c    allocate space for the weight vector as W(POINT_NUM) and the abscissa 
c    vector as X(DIM_NUM,POINT_NUM), and then call GM_RULE_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex 
c    by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c    0 .le. RULE.
c
c    Input, integer DIM_NUM, the spatial dimension.
c    1 .le. DIM_NUM.
c
c    Output, integer POINT_NUM, the number of points in the rule.
c
      implicit none

      integer arg1
      integer dim_num
      integer i4_choose
      integer point_num
      integer rule

      arg1 = dim_num + rule + 1

      point_num = i4_choose ( arg1, rule )

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
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Discussion:
c
c    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
c    bit pattern should be
c
c     01111111111111111111111111111111
c
c    In this case, its numerical value is 2147483647.
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
c  Parameters:
c
c    Output, integer I4_HUGE, a "huge" I4.
c
      implicit none

      integer i4
      integer i4_huge

      i4_huge = 2147483647

      return
      end
      subroutine monomial_value ( dim_num, point_num, x, expon, value )

c*********************************************************************72
c
cc MONOMIAL_VALUE evaluates a monomial.
c
c  Discussion:
c
c    This routine evaluates a monomial of the form
c
c      product ( 1 .le. dim .le. dim_num ) x(dim)^expon(dim)
c
c    where the exponents are nonnegative integers.  Note that
c    if the combination 0^0 is encountered, it should be treated
c    as 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 May 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer POINT_NUM, the number of points at which the
c    monomial is to be evaluated.
c
c    Input, double precision X(DIM_NUM,POINT_NUM), the point coordinates.
c
c    Input, integer EXPON(DIM_NUM), the exponents.
c
c    Output, double precision VALUE(POINT_NUM), the value of the monomial.
c
      implicit none

      integer dim_num
      integer point_num

      integer dim
      integer expon(dim_num)
      integer j
      double precision value(point_num)
      double precision x(dim_num,point_num)

      value(1:point_num) = 1.0D+00

      do dim = 1, dim_num
        if ( 0 .ne. expon(dim) ) then
          do j = 1, point_num
            value(j) = value(j) * x(dim,j) ** expon(dim)
          end do
        end if
      end do

      return
      end
      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial.
c
c  Discussion:
c
c    The formula used is:
c
c      FACTORIAL ( N ) = PRODUCT ( 1 .le. I .le. N ) I
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
c    If N is less than 1, R8_FACTORIAL is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
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
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of double precision values.
c
c    For now, the input quantity SEED is an integer variable.
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
c    Volume 8, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which 
c    should NOT be 0.  On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
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
          seed = seed + i4_huge ( )
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine simplex_unit_monomial_int ( dim_num, expon, value )

c*********************************************************************72
c
cc SIMPLEX_UNIT_MONOMIAL_INT integrates a monomial over a simplex.
c
c  Discussion:
c
c    This routine evaluates a monomial of the form
c
c      product ( 1 .le. dim .le. dim_num ) x(dim)^expon(dim)
c
c    where the exponents are nonnegative integers.  Note that
c    if the combination 0^0 is encountered, it should be treated
c    as 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer EXPON(DIM_NUM), the exponents.
c
c    Output, double precision VALUE, the value of the integral of the
c    monomial.
c
      implicit none

      integer dim_num

      integer dim
      integer expon(dim_num)
      integer i
      integer k
      double precision value
c
c  The first computation ends with VALUE = 1.0;
c
      value = 1.0D+00

      k = 0

      do dim = 1, dim_num

        do i = 1, expon(dim)
          k = k + 1
          value = value * dble ( i ) / dble ( k )
        end do

      end do

      do dim = 1, dim_num

        k = k + 1
        value = value / dble ( k )

      end do

      return
      end
      subroutine simplex_unit_monomial_quadrature ( dim_num, expon, 
     &  point_num, x, w, quad_error )

c*********************************************************************72
c
cc SIMPLEX_UNIT_MONOMIAL_QUADRATURE: quadrature of monomials in a unit simplex.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer EXPON(DIM_NUM), the exponents.
c
c    Input, integer POINT_NUM, the number of points in the rule.
c
c    Input, double precision X(DIM_NUM,POINT_NUM), the quadrature points.
c
c    Input, double precision W(POINT_NUM), the quadrature weights.
c
c    Output, double precision QUAD_ERROR, the quadrature error.
c
      implicit none

      integer dim_num

      double precision exact
      integer expon(dim_num)
      integer point_num
      double precision quad
      double precision quad_error
      double precision r8vec_dot_product
      double precision scale
      double precision value(point_num)
      double precision volume
      double precision w(point_num)
      double precision x(dim_num,point_num)
c
c  Get the exact value of the integral of the unscaled monomial.
c
      call simplex_unit_monomial_int ( dim_num, expon, scale )
c
c  Evaluate the monomial at the quadrature points.
c
      call monomial_value ( dim_num, point_num, x, expon, value )
c
c  Compute the weighted sum and divide by the exact value.
c
      call simplex_unit_volume ( dim_num, volume )
      quad = volume * r8vec_dot_product ( point_num, w, value ) / scale
c
c  Error:
c
      exact = 1.0D+00
      quad_error = abs ( quad - exact )

      return
      end
      subroutine simplex_unit_sample ( dim_num, point_num, seed, x )

c*********************************************************************72
c
cc SIMPLEX_UNIT_SAMPLE returns uniformly random points from a general simplex.
c
c  Discussion:
c
c    The interior of the unit DIM_NUM-dimensional simplex is the set of 
c    points X(1:DIM_NUM) such that each X(I) is nonnegative, and 
c    sum(X(1:DIM_NUM)) .le. 1.
c
c    This routine is valid for any spatial dimension DIM_NUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Reuven Rubinstein,
c    Monte Carlo Optimization, Simulation, and Sensitivity 
c    of Queueing Networks,
c    Krieger, 1992,
c    ISBN: 0894647644,
c    LC: QA298.R79.
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision X(DIM_NUM,POINT_NUM), the points.
c
      implicit none

      integer dim_num
      integer point_num

      double precision e(dim_num+1)
      double precision e_sum
      integer i
      integer j
      double precision r8vec_sum
      integer seed
      double precision x(dim_num,point_num)
c
c  The construction begins by sampling DIM_NUM+1 points from the
c  exponential distribution with parameter 1.
c
      do j = 1, point_num

        call r8vec_uniform_01 ( dim_num + 1, seed, e )

        do i = 1, dim_num + 1
          e(i) = - log ( e(i) )
        end do

        e_sum = r8vec_sum ( dim_num + 1, e )

        do i = 1, dim_num
          x(i,j) = e(i) / e_sum
        end do

      end do

      return
      end
      subroutine simplex_unit_to_general ( dim_num, point_num, t, 
     &  ref, phy )

c*********************************************************************72
c
cc SIMPLEX_UNIT_TO_GENERAL maps the unit simplex to a general simplex.
c
c  Discussion:
c
c    Given that the unit simplex has been mapped to a general simplex
c    with vertices T, compute the images in T, under the same linear
c    mapping, of points whose coordinates in the unit simplex are REF.
c
c    The vertices of the unit simplex are listed as suggested in the
c    following:
c
c      (0,0,0,...,0)
c      (1,0,0,...,0)
c      (0,1,0,...,0)
c      (0,0,1,...,0)
c      (...........)
c      (0,0,0,...,1)
c
c    Thanks to Andrei ("spiritualworlds") for pointing out a mistake in the
c    previous implementation of this routine, 02 March 2008.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer POINT_NUM, the number of points to transform.
c
c    Input, double precision T(DIM_NUM,DIM_NUM+1), the vertices of the
c    general simplex.  
c
c    Input, double precision REF(DIM_NUM,POINT_NUM), points in the 
c    reference triangle.
c
c    Output, double precision PHY(DIM_NUM,POINT_NUM), corresponding points 
c    in the physical triangle.
c
      implicit none

      integer dim_num
      integer point_num

      integer dim
      integer i
      integer j
      double precision phy(dim_num,point_num)
      integer point
      double precision ref(dim_num,point_num)
      double precision t(dim_num,dim_num+1)
      integer vertex
c
c  The image of each point is initially the image of the origin.
c
c  Insofar as the pre-image differs from the origin in a given vertex
c  direction, add that proportion of the difference between the images
c  of the origin and the vertex.
c
      do dim = 1, dim_num 

        do j = 1, point_num
          phy(dim,j) = t(dim,1)
        end do

        do vertex = 2, dim_num + 1

          do j = 1, point_num
            phy(dim,j) = phy(dim,j) 
     &        + ( t(dim,vertex) - t(dim,1) ) * ref(vertex-1,j)
          end do

        end do

      end do

      return
      end
      subroutine simplex_unit_volume ( dim_num, volume )

c*********************************************************************72
c
cc SIMPLEX_UNIT_VOLUME computes the volume of the unit simplex.
c
c  Discussion:
c
c    The formula is simple: volume = 1/N!.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Output, double precision VOLUME, the volume of the cone.
c
      implicit none

      integer i
      integer dim_num
      double precision volume

      volume = 1.0D+00
      do i = 1, dim_num
        volume = volume / dble ( i )
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
