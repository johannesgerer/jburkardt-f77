      subroutine box_nd ( func, dim_num, order, xtab, weight, result, 
     &  eval_num )

c*********************************************************************72
c
cc BOX_ND estimates a multidimensional integral using a product rule.
c
c  Discussion:
c
c    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
c    supplied by the user.  The routine is fairly inflexible.  If
c    you supply a rule for integration from -1 to 1, then your product
c    box must be a product of DIM_NUM copies of the interval [-1,1].
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
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
c    Input, double precision, external FUNC, a routine which evaluates
c    the function to be integrated, of the form:
c      function func ( dim_num, x )
c      integer dim_num
c      double precision func
c      double precision x(dim_num)
c      func = ...
c      return
c      end
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer ORDER, the number of points used in the 1D rule.
c
c    Input, double precision XTAB(ORDER), the abscissas of the 1D rule.
c
c    Input, double precision WEIGHT(ORDER), the weights of the 1D rule.
c
c    Output, double precision RESULT, the approximate value of the integral.
c
c    Output, integer EVAL_NUM, the number of function evaluations.
c
      implicit none

      integer dim_num
      integer order

      integer eval_num
      double precision func
      external func
      integer i
      integer indx(dim_num)
      integer k
      double precision result
      double precision w
      double precision weight(order)
      double precision x(dim_num)
      double precision xtab(order)

      eval_num = 0

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BOX_ND - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM < 1'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( order .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BOX_ND - Fatal error!'
        write ( *, '(a)' ) '  ORDER < 1'
        write ( *, '(a,i8)' ) '  ORDER = ', order
        stop
      end if

      k = 0
      result = 0.0D+00

10    continue

        call tuple_next ( 1, order, dim_num, k, indx )

        if ( k .eq. 0  ) then
          go to 20
        end if

        w = 1.0D+00
        do i = 1, dim_num
          w = w * weight(indx(i))
        end do

        do i = 1, dim_num
          x(i) = xtab(indx(i))
        end do

        result = result + w * func ( dim_num, x )
        eval_num = eval_num + 1

      go to 10

20    continue

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
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
      subroutine monte_carlo_nd ( func, dim_num, a, b, eval_num, seed,
     &  result )

c*********************************************************************72
c
cc MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
c
c  Discussion:
c
c    Unlike the other routines, this routine requires the user to specify
c    the number of function evaluations as an INPUT quantity.
c    No attempt at error estimation is made.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
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
c    Input, double precision, external FUNC, a routine which evaluates
c    the function to be integrated, of the form:
c      function func ( dim_num, x )
c      integer dim_num
c      real ( kind = 8 ) func
c      real ( kind = 8 ) x(dim_num)
c      func = ...
c      return
c      end
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision A(DIM_NUM), B(DIM_NUM), the integration limits.
c
c    Input, integer EVAL_NUM, the number of function evaluations.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision RESULT, the approximate value of the integral.
c
      implicit none

      integer dim_num

      double precision a(dim_num)
      double precision b(dim_num)
      integer dim
      integer eval_num
      double precision func
      external func
      integer i
      double precision result
      integer seed
      double precision volume
      double precision x(dim_num)

      result = 0.0D+00

      do i = 1, eval_num

        call r8vec_uniform_01 ( dim_num, seed, x )

        result = result + func ( dim_num, x )

      end do

      volume = 1.0
      do dim = 1, dim_num
        volume = volume * ( b(dim) - a(dim) )
      end do

      result = result * volume / dble ( eval_num )

      return
      end
      subroutine p5_nd ( func, dim_num, a, b, result, eval_num )

c*********************************************************************72
c
cc P5_ND estimates a multidimensional integral using a formula of exactness 5.
c
c  Discussion:
c
c    The routine uses a method which is exact for polynomials of total
c    degree 5 or less.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    Philip Davis, Philip Rabinowitz.
c
c    FORTRAN90 version by John Burkardt
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
c    Input, double precision, external FUNC, a routine which evaluates
c    the function to be integrated, of the form:
c      function func ( dim_num, x )
c      integer dim_num
c      double precision func
c      double precision x(dim_num)
c      func = ...
c      return
c      end
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision A(DIM_NUM), B(DIM_NUM), the integration limits.
c
c    Output, double precision RESULT, the approximate value of the integral.
c
c    Output, integer EVAL_NUM, the number of function evaluations.
c
      implicit none

      integer dim_num

      double precision a(dim_num)
      double precision a0
      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision b(dim_num)
      integer dim
      double precision en
      integer eval_num
      double precision func
      external func
      integer i
      integer j
      double precision result
      double precision sum1
      double precision sum2
      double precision sum3
      double precision volume
      double precision work(dim_num)

      eval_num = 0

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P5_ND - Fatal error!'
        write ( *, '(a,i8)' ) '  DIM_NUM < 1, DIM_NUM = ', dim_num
        stop
      end if

      a2 = 25.0D+00 / 324.0D+00
      a3 = sqrt ( 0.6D+00 )
      en = dble ( dim_num )
      a0 = ( 25.0D+00 * en * en - 115.0D+00 * en + 162.0D+00 ) 
     &  / 162.0D+00
      a1 = ( 70.0D+00 - 25.0D+00 * en ) / 162.0D+00

      volume = 1.0D+00
      do dim = 1, dim_num
        volume = volume * ( b(dim) - a(dim) )
      end do

      do dim = 1, dim_num
        work(dim) = 0.5D+00 * ( a(dim) + b(dim) )
      end do

      result = 0.0D+00
      if ( volume .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P5_ND - Warning!'
        write ( *, '(a)' ) '  Volume = 0, integral = 0.'
        return
      end if

      sum1 = a0 * func ( dim_num, work )
      eval_num = eval_num + 1

      sum2 = 0.0D+00
      sum3 = 0.0D+00

      do i = 1, dim_num

        work(i) = 0.5D+00 * ( ( a(i) + b(i) ) + a3 * ( b(i) - a(i) ) )
        sum2 = sum2 + func ( dim_num, work )
        eval_num = eval_num + 1

        work(i) = 0.5D+00 * ( ( a(i) + b(i) ) - a3 * ( b(i) - a(i) ) )
        sum2 = sum2 + func ( dim_num, work )
        eval_num = eval_num + 1

        work(i) = 0.5D+00 * ( a(i) + b(i) )

      end do

      if ( 1 .lt. dim_num ) then

        a4 = a3

10      continue

          do i = 1, dim_num - 1

            work(i) = 0.5D+00 * ( ( a(i) + b(i) ) 
     &        + a4 * ( b(i) - a(i) ) )
            a5 = a3

20          continue

              do j = i + 1, dim_num
                work(j) = 0.5D+00 * ( ( a(j) + b(j) ) 
     &            + a5 * ( b(j) - a(j) ) )
                sum3 = sum3 + func ( dim_num, work )
                eval_num = eval_num + 1
                work(j) = 0.5D+00 * ( a(j) + b(j) )
              end do

              a5 = -a5

              if ( 0.0D+00 <= a5 ) then
                go to 30
              end if

            go to 20

30          continue

            work(i) = 0.5D+00 * ( a(i) + b(i) )

          end do

          a4 = -a4

          if ( 0.0D+00 .le. a4 ) then
            go to 40
          end if

        go to 10

40      continue

      end if

      result = volume * ( sum1 + a1 * sum2 + a2 * sum3 )

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
      subroutine romberg_nd ( n, a, b, h, al, m, e, f, aint, key, 
     &  eval_num )

c*********************************************************************72
c
cc ROMBERG_ND approximates an integral over an N-dimensional box.
c
c  Discussion:
c
c    In the reference, this routine was named "NDIMRI".
c
c    The routine uses a Romberg-type method based on the midpoint rule.
c    Initially, each spatial interval [A(I),B(I)] is subdivided into
c    subintervals of width H(I).  On the second iteration, the number
c    of points in each direction is doubled, and thereafter, it is
c    increased by a factor of AL on each iteration.
c
c    This routine, however ragged looking, is close to the form of
c    the original as published in the reference.  It is intended to
c    be available as a benchmark against which more-modernized versions
c    can be tested.  Therefore, great effort has been expended to
c    exert as little effort as possible in cleaning up or improving
c    this version of the algorithm!
c
c  Modified:
c
c    11 September 2006
c
c  Author:
c
c    Philip Davis, Philip Rabinowitz
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
c    Input, integer N, the dimension of the region which is
c    to be integrated.  N must be less than 10.
c
c    Input, double precision A(N), contains the lower limits of integration.
c
c    Input, double precision B(N), contains the upper limits of integration.
c
c    Input, double precision H(N), an array of fractions of the form 1/K,
c    representing the initial spacing in each spatial dimension.
c
c    Input, double precision AL, the growth rate for the number of points in
c    each dimension.  1.5 <= AL <= 2.0.
c
c    Input, integer M, the maximum number of iterations to
c    be performed.  The number of function evaluations on
c    iteration J is at least J**N, which grows very rapidly.
c    M should be small!
c
c    Input, double precision E, an error tolerance for the approximation
c    of the integral.
c
c    Input, double precision, external FUNC, the name of a routine 
c    which evaluates the function to be integrated, of the form:
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c      func = ...
c      return
c      end
c
c    Output, double precision AINT, the approximate value of the integral.
c
c    Output, integer KEY, result key.
c    1, for satisfactory convergence.
c    0, the input was not satisfactory.
c    -1, if the maximum number of iterations was reached without convergence.
c
c    Output, integer EVAL_NUM, the number of function evaluations.
c
      implicit none

      integer n

      integer maxit
      parameter ( maxit = 15 )

      double precision a(n)
      double precision aa(maxit)
      double precision aint
      double precision al
      double precision b(n)
      double precision c(9)
      double precision d(9)
      double precision e
      double precision ee
      double precision en
      integer eval_num
      external f
      double precision f
      double precision g(9)
      double precision h(n)
      integer i
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer i6
      integer i7
      integer i8
      integer i9
      integer k(maxit)
      integer key
      integer l
      integer ll
      integer m
      integer nn(9)
      double precision p(9)
      double precision tol
      double precision u
      double precision v(maxit)
      double precision x(9)

      data k(1) / 1 /
      data k(2) / 2 /
      data tol / 1.0D-07 /

      eval_num = 0
c
c  Make sure the input is satisfactory.
c
      if ( n .lt. 1 ) then
        key = 0
        return
      end if

      if ( n .gt. 9 ) then
        key = 0
        return
      end if
c
c  Original code pretends to allow 15 < M, but silently
c  ignores values greater than 15.
c
      if ( m .lt. 1 .or. 15 .lt. m ) then
        key = 0
        return
      end if

      if ( al .lt. 1.5D+00 .or. al .gt. 2.0D+00 ) then
        key = 0
        return
      end if

      do i = 1, n
        if ( h(i) .le. 0.0D+00 .or. h(i) .gt. 1.0D+00 ) then
          key = 0
          return
        end if
      end do
c
c  Initialization.
c
      do i = 1, n
        d(i) = a(i)
      end do

      ee = dmax1 ( e, tol )
c
c  The original code loops from n to 8 and indexes by i+1,
c  a needless complication.
c
      do i = n+1, 9
        p(i) = 0.0D+00
        nn(i) = 1
        d(i) = 0.0D+00
      end do

      l = 1
c
c  Original code has b(i) - d(i), but this is pointless obfuscation.
c
      do i = 1, n
        c(i) = b(i) - a(i)
      end do
c
c  The Romberg loop.
c
21    continue

      u = 0.0D+00

      do i = 1, n
        g(i) = h(i) / k(l)
        nn(i) = 1.0D+00 / g(i) + 0.5D+00
        p(i) = c(i) * g(i)
      end do

      do i9 = 1, nn(9)
        x(9) = d(9) + p(9) * ( i9 - 0.5D+00 )

        do i8 = 1, nn(8)
          x(8) = d(8) + p(8) * ( i8 - 0.5D+00 )

          do i7 = 1, nn(7)
            x(7) = d(7) + p(7) * ( i7 - 0.5D+00 )

            do i6 = 1, nn(6)
              x(6) = d(6) + p(6) * ( i6 - 0.5D+00 )

              do i5 = 1, nn(5)
                x(5) = d(5) + p(5) * ( i5 - 0.5D+00 )

                do i4 = 1, nn(4)
                  x(4) = d(4) + p(4) * ( i4 - 0.5D+00 )

                  do i3 = 1, nn(3)
                    x(3) = d(3) + p(3) * ( i3 - 0.5D+00 )

                    do i2 = 1, nn(2)
                      x(2) = d(2) + p(2) * ( i2 - 0.5D+00 )

                      do i1 = 1, nn(1)
                        x(1) = d(1) + p(1) * ( i1 - 0.5D+00 )

                        u = u + f ( n, x )
                        eval_num = eval_num + 1

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      do i = 1, n
        u = u * p(i)
      end do

      v(l) = u

      if ( l <= 1 ) then
        aa(1) = v(1)
        l = l + 1
        go to 21
      end if

      en = k(l)

      do ll = 2, l
        i = l + 1 - ll
        v(i) = v(i+1) + ( v(i+1) - v(i) ) 
     &    / ( ( en / k(i) )**2 - 1.0D+00 )
      end do

      aint = v(1)

      if ( abs ( aint - aa(l-1) ) .lt. abs ( aint * ee ) ) then
        key = 1
        return
      end if

      if ( m .le. l ) then
        key = -1
        return
      end if

      aa(l) = aint
      l = l + 1
      k(l) = al * k(l-1)
      go to 21

      end
      subroutine sample_nd ( func, k1, k2, dim_num, est1, err1, dev1, 
     &  est2, err2, dev2, eval_num )

c*********************************************************************72
c
cc SAMPLE_ND estimates a multidimensional integral using sampling.
c
c  Discussion:
c
c    This routine computes two sequences of integral estimates, EST1 
c    and EST2, for indices K going from K1 to K2.  These estimates are 
c    produced by the generation of 'random' abscissas in the region.  
c    The process can become very expensive if high accuracy is needed.
c
c    The total number of function evaluations is
c    4*(K1**DIM_NUM+(K1+1)**DIM_NUM+...+(K2-1)**DIM_NUM+K2**DIM_NUM), and K2
c    should be chosen so as to make this quantity reasonable.
c    In most situations, EST2(K) are much better estimates than EST1(K).
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    Philip Davis, Philip Rabinowitz.
c
c    This FORTRAN77 version by John Burkardt
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
c    Input, double precision, external FUNC, a routine which evaluates
c    the function to be integrated, of the form:
c      function func ( dim_num, x )
c      integer dim_num
c      double precision func
c      double precision x(dim_num)
c      func = ...
c      return
c      end
c
c    Input, integer K1, the beginning index for the iteration.
c    1 <= K1 <= K2.
c
c    Input, integer K2, the final index for the iteration.  K1 <= K2.
c    Increasing K2 increases the accuracy of the calculation,
c    but vastly increases the work and running time of the code.
c
c    Input, integer DIM_NUM, the spatial dimension.  1 <= DIM_NUM <= 10.
c
c    Output, double precision EST1(K2).  Entries K1 through K2 contain
c    successively better estimates of the integral.
c
c    Output, double precision ERR1(K2).  Entries K1 through K2 contain
c    the corresponding estimates of the integration errors.
c
c    Output, double precision DEV1(K2).  Entries K1 through K2 contain
c    estimates of the reliability of the the integration.
c    If consecutive values DEV1(K) and DEV1(K+1) do not differ
c    by more than 10 percent, then ERR1(K) can be taken as
c    a reliable upper bound on the difference between EST1(K)
c    and the true value of the integral.
c
c    Output, double precision EST2(K2).  Entries K2 through K2 contain
c    successively better estimates of the integral.
c
c    Output, double precision ERR2(K2).  Entries K2 through K2 contain
c    the corresponding estimates of the integration errors.
c
c    Output, double precision DEV2(K2).  Entries K2 through K2 contain
c    estimates of the reliability of the the integration.
c    If consecutive values DEV2(K) and DEV2(K+2) do not differ
c    by more than 10 percent, then ERR2(K) can be taken as
c    a reliable upper bound on the difference between EST2(K)
c    and the true value of the integral.
c
c    Output, integer EVAL_NUM, the number of function evaluations.
c
      implicit none

      integer k2
      integer dim_max
      parameter ( dim_max = 10 )
      integer dim_num

      double precision ak
      double precision ak1
      double precision akn
      double precision al(dim_max)
      double precision b
      double precision be(dim_max)
      double precision bk
      double precision d1
      double precision d2
      double precision dev1(k2)
      double precision dev2(k2)
      double precision dex(dim_max)
      integer dim
      double precision err1(k2)
      double precision err2(k2)
      double precision est1(k2)
      double precision est2(k2)
      integer eval_num
      double precision func
      external func
      double precision g
      double precision ga(dim_max)
      integer i
      integer j
      integer k
      integer k1
      integer key
      double precision p1(dim_max)
      double precision p2(dim_max)
      double precision p3(dim_max)
      double precision p4(dim_max)
      double precision s1
      double precision s2
      double precision t
      double precision y1
      double precision y2
      double precision y3
      double precision y4

      save al

      data al /
     &  0.4142135623730950D+00, 
     &  0.7320508075688773D+00, 
     &  0.2360679774997897D+00, 
     &  0.6457513110645906D+00, 
     &  0.3166247903553998D+00, 
     &  0.6055512754639893D+00, 
     &  0.1231056256176605D+00, 
     &  0.3589989435406736D+00, 
     &  0.7958315233127195D+00, 
     &  0.3851648071345040D+00 /

      eval_num = 0
c
c  Check input
c
      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_ND - Fatal errorc'
        write ( *, '(a)' ) '  DIM_NUM must be at least 1,'
        write ( *, '(a,i8)' ) '  but DIM_NUM = ', dim_num
        stop
      end if
 
      if ( dim_max .lt. dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  DIM_NUM must be no more than DIM_MAX = ', dim_max
        write ( *, '(a,i8)' ) '  but DIM_NUM = ', dim_num
        stop
      end if
 
      if ( k1 .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
        write ( *, '(a,i8)' ) '  K1 must be at least 1, but K1 = ', k1
        stop
      end if
 
      if ( k2 .lt. k1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_ND - Fatal error!'
        write ( *, '(a)' ) '  K1 may not be greater than K2, but '
        write ( *, '(a,i8)' ) '  K1 = ', k1
        write ( *, '(a,i8)' ) '  K2 = ', k2
        stop
      end if
 
      do dim = 1, dim_num
        be(dim) = al(dim)
        ga(dim) = al(dim)
        dex(dim) = 0.0D+00
      end do
 
      do k = k1, k2
 
        ak = dble ( k )
        key = 0
        ak1 = ak - 1.1D+00
        s1 = 0.0D+00
        d1 = 0.0D+00
        s2 = 0.0D+00
        d2 = 0.0D+00
        akn = ak**dim_num
        t = sqrt ( ak**dim_num ) * ak
        bk = 1.0D+00 / ak

10      continue 
 
          key = key + 1

          if ( key .eq. 1 ) then
            go to 50
          end if

          key = key - 1
          j = 1
 
30        continue
     
          if ( dex(j) .le. ak1 ) then
            dex(j) = dex(j) + 1.0D+00
            go to 50
          end if
 
          dex(j) = 0.0D+00
          j = j + 1
          if ( j .le. dim_num ) then
            go to 30
          end if

          go to 70
 
50        continue
 
          do i = 1, dim_num

            b = be(i) + al(i)
            if ( 1.0D+00 .lt. b ) then
              b = b - 1.0D+00
            end if

            g = ga(i) + b
            if ( 1.0D+00 .lt. g ) then
              g = g - 1.0D+00
            end if

            be(i) = b + al(i)
            if ( 1.0D+00 .lt. be(i) ) then
              be(i) = be(i) - 1.0D+00
            end if

            ga(i) = be(i) + g
            if ( 1.0D+00 .lt. ga(i) ) then
              ga(i) = ga(i) - 1.0D+00
            end if

            p1(i) = ( dex(i) + g ) * bk
            p2(i) = ( dex(i) + 1.0D+00 - g ) * bk
            p3(i) = ( dex(i) + ga(i) ) * bk
            p4(i) = ( dex(i) + 1.0D+00 - ga(i) ) * bk

          end do
 
          y1 = func ( dim_num, p1 )
          eval_num = eval_num + 1
c
c  There may be an error in the next two lines,
c  but oddly enough, that is how the original reads
c
          y3 = func ( dim_num, p2 )
          eval_num = eval_num + 1
          y2 = func ( dim_num, p3 )
          eval_num = eval_num + 1
          y4 = func ( dim_num, p4 )
          eval_num = eval_num + 1

          s1 = s1 + y1 + y2
          d1 = d1 + ( y1 - y2 )**2
          s2 = s2 + y3 + y4
          d2 = d2 + ( y1 + y3 - y2 - y4 )**2

        go to 10

70      continue
  
        est1(k) = 0.5D+00 * s1 / akn
        err1(k) = 1.5D+00 * sqrt ( d1 ) / akn
        dev1(k) = err1(k) * t
        est2(k) = 0.25D+00 * ( s1 + s2 ) / akn
        err2(k) = 0.75D+00 * sqrt ( d2 ) / akn
        dev2(k) = err2(k) * t * ak
 
      end do
 
      return
      end
      subroutine sum2_nd ( func, xtab, weight, order, dim_num, result, 
     &  eval_num )

c*********************************************************************72
c
cc SUM2_ND estimates a multidimensional integral using a product rule.
c
c  Discussion:
c
c    The routine uses a product rule supplied by the user.
c
c    The region may be a product of any combination of finite,
c    semi-infinite, or infinite intervals.
c
c    For each factor in the region, it is assumed that an integration
c    rule is given, and hence, the region is defined implicitly by
c    the integration rule chosen.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    Philip Davis, Philip Rabinowitz.
c
c    FORTRAN90 version by John Burkardt
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
c    Input, double precision, external FUNC, the name of a routine which evaluates
c    the function to be integrated, of the form:
c      function func ( dim_num, x )
c      integer dim_num
c      double precision func
c      double precision x(dim_num)
c      func = ...
c      return
c      end
c
c    Input, double precision XTAB(DIM_NUM,ORDER_MAX).  XTAB(I,J) is the 
c    I-th abscissa of the J-th rule.
c
c    Input, double precision WEIGHT(DIM_NUM,ORDER_MAX).  WEIGHT(I,J) is the 
c    I-th weight for the J-th rule.
c
c    Input, integer ORDER(DIM_NUM).  ORDER(I) is the number of
c    abscissas to be used in the J-th rule.  ORDER(I) must be
c    greater than 0 and less than or equal to ORDER_MAX.
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Output, double precision RESULT, the approximate value of the integral.
c
c    Output, integer EVAL_NUM, the number of function evaluations.
c
      implicit none

      integer dim_max
      parameter ( dim_max = 10 )
      integer dim_num

      integer dim
      integer eval_num
      double precision func
      external func
      integer i
      integer iwork(dim_max)
      integer k
      integer m1
      integer order(dim_num)
      double precision result
      double precision w1
      double precision weight(dim_num,*)
      double precision work(dim_max)
      double precision xtab(dim_num,*)
c
c  Default values.
c
      result = 0.0D+00
      eval_num = 0

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUM2_ND - Fatal errorc'
        write ( *, '(a)' ) '  DIM_NUM < 1'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if
     
      if ( dim_max .lt. dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUM2_ND - Fatal errorc'
        write ( *, '(a)' ) '  DIM_MAX < DIM_NUM'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      do i = 1, dim_num

        if ( order(i) .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SUM2_ND - Fatal errorc'
          write ( *, '(a)' ) '  ORDER(I) < 1.'
          write ( *, '(a,i8)' ) '  For I = ', i
          write ( *, '(a,i8)' ) '  ORDER(I) = ', order(i)
          stop
        end if

      end do
 
      do dim = 1, dim_num
        iwork(dim) = 1
      end do

10    continue 
 
        k = 1
 
        w1 = 1.0D+00
        do i = 1, dim_num
          m1 = iwork(i)
          work(i) = xtab(i,m1)
          w1 = w1 * weight(i,m1)
        end do

        result = result + w1 * func ( dim_num, work )
        eval_num = eval_num + 1

20      continue

          if ( iwork(k) .ne. order(k) ) then
            go to 30
          end if

          iwork(k) = 1
          k = k + 1

          if ( dim_num < k ) then
            go to 40
          end if

        go to 20
 
30      continue

        iwork(k) = iwork(k) + 1

      go to 10

40    continue

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
      subroutine tuple_next ( m1, m2, n, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT computes the next element of a tuple space.
c
c  Discussion:
c
c    The elements are N vectors.  Each entry is constrained to lie
c    between M1 and M2.  The elements are produced one at a time.
c    The first element is
c      (M1,M1,...,M1),
c    the second element is
c      (M1,M1,...,M1+1),
c    and the last element is
c      (M2,M2,...,M2)
c    Intermediate elements are produced in lexicographic order.
c
c  Example:
c
c    N = 2, M1 = 1, M2 = 3
c
c    INPUT        OUTPUT
c    -------      -------
c    Rank  X      Rank   X
c    ----  ---    -----  ---
c    0     * *    1      1 1
c    1     1 1    2      1 2
c    2     1 2    3      1 3
c    3     1 3    4      2 1
c    4     2 1    5      2 2
c    5     2 2    6      2 3
c    6     2 3    7      3 1
c    7     3 1    8      3 2
c    8     3 2    9      3 3
c    9     3 3    0      0 0
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
c    Input, integer M1, M2, the minimum and maximum entries.
c
c    Input, integer N, the number of components.
c
c    Input/output, integer RANK, counts the elements.
c    On first call, set RANK to 0.  Thereafter, the output value of RANK
c    will indicate the order of the element returned.  When there are no
c    more elements, RANK will be returned as 0.
c
c    Input/output, integer X(N), on input the previous tuple.
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer m1
      integer m2
      integer rank
      integer x(n)

      if ( m2 .lt. m1 ) then
        rank = 0
        return
      end if

      if ( rank .le. 0 ) then

        do i = 1, n
          x(i) = m1
        end do
        rank = 1

      else

        rank = rank + 1
        i = n

10      continue

          if ( x(i) .lt. m2 ) then
            x(i) = x(i) + 1
            go to 20
          end if

          x(i) = m1

          if ( i .eq. 1 ) then
            rank = 0
            do i = 1, n
              x(i) = m1
            end do
            go to 20
          end if

          i = i - 1

        go to 10

20      continue

      end if

      return
      end
