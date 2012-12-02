      subroutine hermite_compute ( order, xtab, weight )

c*********************************************************************72
c
cc HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
c
c  Discussion:
c
c    The abscissas are the zeros of the N-th order Hermite polynomial.
c
c    The integration interval is ( -oo, +oo ).
c
c    The weight function is w(x) = exp ( - x * x ).
c
c    The integral to approximate:
c
c      integral ( -oo .lt. X .lt. +oo ) exp ( - X * X ) * F(X) dX
c
c    The quadrature rule:
c
c      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Input, integer ORDER, the order of the formula.
c
c    Output, double precision XTAB(ORDER), the abscissas.
c
c    Output, double precision WEIGHT(ORDER), the weights.
c
      implicit none

      integer order

      double precision cc
      double precision dp2
      integer i
      double precision p1
      double precision r8_gamma
      double precision s
      double precision temp
      double precision weight(order)
      double precision x
      double precision xtab(order)

      cc = 1.7724538509D+00 * r8_gamma ( dble ( order ) ) 
     &  / ( 2.0D+00**( order - 1 ) )

      s = ( 2.0D+00 * dble ( order ) + 1.0D+00 )**( 1.0D+00 / 6.0D+00 )

      do i = 1, ( order + 1 ) / 2

        if ( i .eq. 1 ) then

          x = s**3 - 1.85575D+00 / s

        else if ( i .eq. 2 ) then

          x = x - 1.14D+00 * ( ( dble ( order ) )**0.426D+00 ) / x

        else if ( i .eq. 3 ) then

          x = 1.86D+00 * x - 0.86D+00 * xtab(1)

        else if ( i .eq. 4 ) then

          x = 1.91D+00 * x - 0.91D+00 * xtab(2)

        else

          x = 2.0D+00 * x - xtab(i-2)

        end if

        call hermite_root ( x, order, dp2, p1 )

        xtab(i) = x
        weight(i) = ( cc / dp2 ) / p1

        xtab(order-i+1) = - x
        weight(order-i+1) = weight(i)

      end do
c
c  Reverse the order.
c
      do i = 1, order / 2
        temp            = xtab(i)
        xtab(i)         = xtab(order+1-i)
        xtab(order+1-i) = temp
      end do

      return
      end
      subroutine hermite_integral ( n, value )

c*********************************************************************72
c
cc HERMITE_INTEGRAL returns the value of a Hermite polynomial integral.
c
c  Discussion:
c
c    H(n) = Integral ( -oo .lt. x .lt. +oo ) x^n exp(-x^2) dx
c
c    H(n) is 0 for n odd.
c
c    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the integral.  
c    0 <= N.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer i4_factorial2
      integer n
      double precision pi_sqrt
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )
      double precision r8_huge
      double precision value

      if ( n .lt. 0 ) then

        value = - r8_huge ( )

      else if ( mod ( n, 2 ) .eq. 1 ) then

        value = 0.0D+00

      else

        value = dble ( i4_factorial2 ( n - 1 ) ) * pi_sqrt 
     &    / 2.0D+00**( n / 2 )

      end if

      return
      end
      subroutine hermite_recur ( p2, dp2, p1, x, order )

c*********************************************************************72
c
cc HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Output, double precision P2, the value of H(ORDER)(X).
c
c    Output, double precision DP2, the value of H'(ORDER)(X).
c
c    Output, double precision P1, the value of H(ORDER-1)(X).
c
c    Input, double precision X, the point at which polynomials are evaluated.
c
c    Input, integer ORDER, the order of the polynomial 
c    to be computed.
c
      implicit none

      integer i
      double precision dp0
      double precision dp1
      double precision dp2
      integer order
      double precision p0
      double precision p1
      double precision p2
      double precision x

      p1 = 1.0D+00
      dp1 = 0.0D+00

      p2 = x
      dp2 = 1.0D+00

      do i = 2, order

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2  = x * p1 - 0.5D+00 * ( dble ( i ) - 1.0D+00 ) * p0
        dp2 = x * dp1 + p1 - 0.5D+00 * ( dble ( i ) - 1.0D+00 ) * dp0

      end do

      return
      end
      subroutine hermite_root ( x, order, dp2, p1 )

c*********************************************************************72
c
cc HERMITE_ROOT improves an approximate root of a Hermite polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Input/output, double precision X, the approximate root, which
c    should be improved on output.
c
c    Input, integer ORDER, the order of the Hermite polynomial.
c
c    Output, double precision DP2, the value of H'(ORDER)(X).
c
c    Output, double precision P1, the value of H(ORDER-1)(X).
c
      implicit none

      double precision d
      double precision dp2
      double precision eps
      parameter ( eps = 1.0D-12 )
      integer order
      double precision p1
      double precision p2
      integer step
      integer step_max
      parameter ( step_max = 10 )
      double precision x

      do step = 1, step_max

        call hermite_recur ( p2, dp2, p1, x, order )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0D+00 ) ) then
          return
        end if

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
      subroutine p00_exact ( problem, exact )

c*********************************************************************72
c
cc P00_EXACT returns the exact integral for any problem.
c
c  Discussion:
c
c    This routine provides a "generic" interface to the exact integral
c    routines for the various problems, and allows a problem to be called
c    by index (PROBLEM) rather than by name.
c
c    In most cases, the "exact" value of the integral is not given;
c    instead a "respectable" approximation is available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Output, double precision EXACT, the exact value of the integral.
c
      implicit none

      double precision exact
      integer problem

      if ( problem .eq. 1 ) then
        call p01_exact ( exact )
      else if ( problem .eq. 2 ) then
        call p02_exact ( exact )
      else if ( problem .eq. 3 ) then
        call p03_exact ( exact )
      else if ( problem .eq. 4 ) then
        call p04_exact ( exact )
      else if ( problem .eq. 5 ) then
        call p05_exact ( exact )
      else if ( problem .eq. 6 ) then
        call p06_exact ( exact )
      else if ( problem .eq. 7 ) then
        call p07_exact ( exact )
      else if ( problem .eq. 8 ) then
        call p08_exact ( exact )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_fun ( problem, option, n, x, f )

c*********************************************************************72
c
cc P00_FUN evaluates the integrand for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer option
      integer problem
      double precision x(n)

      if ( problem .eq. 1 ) then
        call p01_fun ( option, n, x, f )
      else if ( problem .eq. 2 ) then
        call p02_fun ( option, n, x, f )
      else if ( problem .eq. 3 ) then
        call p03_fun ( option, n, x, f )
      else if ( problem .eq. 4 ) then
        call p04_fun ( option, n, x, f )
      else if ( problem .eq. 5 ) then
        call p05_fun ( option, n, x, f )
      else if ( problem .eq. 6 ) then
        call p06_fun ( option, n, x, f )
      else if ( problem .eq. 7 ) then
        call p07_fun ( option, n, x, f )
      else if ( problem .eq. 8 ) then
        call p08_fun ( option, n, x, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FUN - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_gauss_hermite ( problem, order, result )

c*********************************************************************72
c
cc P00_GAUSS_HERMITE applies a Gauss-Hermite quadrature rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, integer ORDER, the order of the Gauss-Laguerre rule
c    to apply.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer order

      double precision f_vec(order)
      integer option
      integer problem
      double precision result
      double precision weight(order)
      double precision xtab(order)

      call hermite_compute ( order, xtab, weight )

      option = 1
      call p00_fun ( problem, option, order, xtab, f_vec )

      result = dot_product ( weight(1:order), f_vec(1:order) )

      return
      end
      subroutine p00_monte_carlo ( problem, order, result )

c*********************************************************************72
c
cc P00_MONTE_CARLO applies a Monte Carlo procedure to Hermite integrals.
c
c  Discussion:
c
c    We wish to estimate the integral:
c
c      I(f) = integral ( -oo .lt. x .lt. +oo ) f(x) exp ( - x * x ) dx
c
c    We do this by a Monte Carlo sampling procedure, in which 
c    we select N points X(1:N) from a standard normal distribution,
c    and estimate
c
c      Q(f) = sum ( 1 <= I <= N ) f(x(i)) / sqrt ( pi )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, integer ORDER, the order of the Gauss-Laguerre rule
c    to apply.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision f_vec(order)
      integer option
      integer order
      double precision pi_sqrt
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )
      integer problem
      double precision r8vec_sum
      double precision result
      integer seed
      double precision weight
      double precision x_vec(order)

      seed = 123456789
      call r8vec_normal_01 ( order, seed, x_vec )

      option = 2
      call p00_fun ( problem, option, order, x_vec, f_vec )

      weight = dble ( order ) / pi_sqrt / sqrt ( 2.0D+00 )

      result = r8vec_sum ( order, f_vec ) / weight

      return
      end
      subroutine p00_problem_num ( problem_num )

c*********************************************************************72
c
cc P00_PROBLEM_NUM returns the number of test integration problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROBLEM_NUM, the number of test problems.
c
      implicit none

      integer problem_num

      problem_num = 8

      return
      end
      subroutine p00_title ( problem, title )

c*********************************************************************72
c
cc P00_TITLE returns the title for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      integer problem
      character ( len = * ) title

      if ( problem .eq. 1 ) then
        call p01_title ( title )
      else if ( problem .eq. 2 ) then
        call p02_title ( title )
      else if ( problem .eq. 3 ) then
        call p03_title ( title )
      else if ( problem .eq. 4 ) then
        call p04_title ( title )
      else if ( problem .eq. 5 ) then
        call p05_title ( title )
      else if ( problem .eq. 6 ) then
        call p06_title ( title )
      else if ( problem .eq. 7 ) then
        call p07_title ( title )
      else if ( problem .eq. 8 ) then
        call p08_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_turing ( problem, h, tol, n, result )

c*********************************************************************72
c
cc P00_TURING applies the Turing quadrature rule.
c
c  Discussion:
c
c    We consider the approximation:
c
c      Integral ( -oo .lt. x .lt. +oo ) f(x) dx
c
c      = h * Sum ( -oo .lt. i .lt. +oo ) f(nh) + error term
c
c    Given H and a tolerance TOL, we start summing at I = 0, and
c    adding one more term in the positive and negative I directions,
c    until the absolute value of the next two terms being added 
c    is less than TOL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Turing,
c    A Method for the Calculation of the Zeta Function,
c    Proceedings of the London Mathematical Society,
c    Volume 48, 1943, pages 180-197.
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, double precision H, the spacing to use.
c
c    Input, double precision TOL, the tolerance.  
c
c    Output, integer N, the number of pairs of steps taken.
c    The actual number of function evaluations is 2*N+1.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision f_vec(2)
      double precision h
      integer n
      integer n_too_many
      parameter ( n_too_many = 100000 )
      integer option
      integer order
      integer problem
      double precision result
      double precision tol
      double precision xtab(2)

      option = 0
      n = 0

      result = 0.0D+00
      order = 1
      xtab(1) = 0.0D+00
      call p00_fun ( problem, option, order, xtab, f_vec )
      result = result + h * f_vec(1)

10    continue

        n = n + 1

        xtab(1) =   dble ( n ) * h
        xtab(2) = - dble ( n ) * h

        order = 2
        call p00_fun ( problem, option, order, xtab, f_vec )

        result = result + h * ( f_vec(1) + f_vec(2) )
c
c  Just do a simple-minded absolute error tolerance check to start with.
c
        if ( abs ( f_vec(1) ) + abs ( f_vec(2) ) <= tol ) then
          go to 20
        end if
c
c  Just in case things go crazy.
c
        if ( n_too_many <= n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P00_TURING - Warning!'
          write ( *, '(a,i8)' ) 
     &      '  Number of steps exceeded N_TOO_MANY = ', n_too_many
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine p01_exact ( exact )

c*********************************************************************72
c
cc P01_EXACT returns the exact integral for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision omega
      parameter ( omega = 1.0D+00 )
      double precision pi_sqrt 
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )

      exact = pi_sqrt * exp ( - omega * omega )

      return
      end
      subroutine p01_fun ( option, n, x, f )

c*********************************************************************72
c
cc P01_FUN evaluates the integrand for problem 1.
c
c  Discussion:
c
c    Squire gives exact value as sqrt(pi) * exp(-w*w).
c
c    Integral ( -oo .lt. x .lt. +oo ) exp(-x*x) cos(2*w*x) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Squire,
c    Comparison of Gauss-Hermite and Midpoint Quadrature with Application
c    to the Voigt Function,
c    in Numerical Integration: 
c    Recent Developments, Software and Applications,
c    edited by Patrick Keast, Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision omega
      parameter ( omega = 1.0D+00 )
      integer option
      double precision x(n)

      do i = 1, n 
        f(i) = cos ( 2.0D+00 * omega * x(i) )
      end do

      if ( option .eq. 0 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - x(i) * x(i) )
        end do
      else if ( option .eq. 1 ) then

      else if ( option .eq. 2 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - 0.5D+00 * x(i) * x(i) )
        end do
      end if

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the title for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'cos(2*omega*x)'

      return
      end
      subroutine p02_exact ( exact )

c*********************************************************************72
c
cc P02_EXACT returns the exact integral for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi_sqrt
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )

      exact = pi_sqrt

      return
      end
      subroutine p02_fun ( option, n, x, f )

c*********************************************************************72
c
cc P02_FUN evaluates the integrand for problem 2.
c
c  Discussion:
c
c    The exact value is sqrt(pi).
c
c    Integral ( -oo .lt. x .lt. +oo ) exp(-x*x) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision x(n)

      do i = 1, n
        f(i) = 1.0D+00
      end do

      if ( option .eq. 0 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - x(i) * x(i) )
        end do
      else if ( option .eq. 1 ) then

      else if ( option .eq. 2 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - 0.5D+00 * x(i) * x(i) )
        end do
      end if

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'exp(-x*x)'

      return
      end
      subroutine p03_exact ( exact )

c*********************************************************************72
c
cc P03_EXACT returns the exact integral for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision p
      parameter ( p = 1.0D+00 )
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )
      double precision q
      parameter ( q = 3.0D+00 )

      exact = pi / ( q * sin ( pi * p / q ) )

      return
      end
      subroutine p03_fun ( option, n, x, f )

c*********************************************************************72
c
cc P03_FUN evaluates the integrand for problem 3.
c
c  Discussion:
c
c    The exact value is pi / (q*sin(pi*p/q) ), assuming 0 .lt. p .lt. q.
c
c    Integral ( -oo .lt. x .lt. +oo ) exp(-px) / ( 1 + exp ( -qx) ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision p
      parameter ( p = 1.0D+00 )
      double precision q
      parameter ( q = 3.0D+00 )
      double precision x(n)

      do i = 1, n
        f(i) = exp ( - p * x(i) ) / ( 1.0D+00 + exp ( -q * x(i) ) )
      end do

      if ( option .eq. 0 ) then

      else if ( option .eq. 1 ) then
        do i = 1, n
          f(i) = f(i) * exp ( + x(i) * x(i) )
        end do
      else if ( option .eq. 2 ) then
        do i = 1, n
          f(i) = f(i) * exp ( + 0.5D+00 * x(i) * x(i) )
        end do
      end if

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'exp(-px) / ( 1 + exp(-qx) )'

      return
      end
      subroutine p04_exact ( exact )

c*********************************************************************72
c
cc P04_EXACT returns the exact integral for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = sqrt ( pi / 2.0D+00 )

      return
      end
      subroutine p04_fun ( option, n, x, f )

c*********************************************************************72
c
cc P04_FUN evaluates the integrand for problem 4.
c
c  Discussion:
c
c    The exact value is sqrt ( pi / 2 )
c
c    Integral ( -oo .lt. x .lt. +oo ) sin ( x*x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision x(n)

      do i = 1, n
        f(i) = sin ( x(i)**2 )
      end do

      if ( option .eq. 0 ) then

      else if ( option .eq. 1 ) then
        do i = 1, n
          f(i) = f(i) * exp ( + x(i) * x(i) )
        end do
      else if ( option .eq. 2 ) then
        do i = 1, n
          f(i) = f(i) * exp ( + 0.5D+00 * x(i) * x(i) )
        end do
      end if

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'sin(x^2)'

      return
      end
      subroutine p05_exact ( exact )

c*********************************************************************72
c
cc P05_EXACT returns the exact integral for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = pi / 3.0D+00

      return
      end
      subroutine p05_fun ( option, n, x, f )

c*********************************************************************72
c
cc P05_FUN evaluates the integrand for problem 5.
c
c  Discussion:
c
c    The exact value is pi / 3.
c
c    Integral ( -oo .lt. x .lt. +oo ) dx / ( (1+x^2) sqrt(4+3x^2) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision x(n)

      do i = 1, n
        f(i) = 1.0D+00 / ( ( 1.0D+00 + x(i)**2 ) 
     &    * sqrt ( 4.0D+00 + 3.0D+00 * x(i)**2 ) )
      end do

      if ( option .eq. 0 ) then

      else if ( option .eq. 1 ) then
        do i = 1, n
          f(i) = f(i) * exp ( x(i)**2 )
        end do
      else if ( option .eq. 2 ) then
        do i = 1, n
          f(i) = f(i) * exp ( 0.5D+00 * x(i)**2 )
        end do
      end if

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = '1/( (1+x^2) sqrt(4+3x^2) )'

      return
      end
      subroutine p06_exact ( exact )

c*********************************************************************72
c
cc P06_EXACT returns the exact integral for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      integer i4_factorial2
      integer m
      double precision pi_sqrt
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )
      double precision r8_huge

      call p06_param ( 'G', 'M', m )

      if ( m .le. -1 ) then

        exact = - r8_huge ( )

      else if ( mod ( m, 2 ) .eq. 1 ) then

        exact = 0.0D+00

      else

        exact = dble ( i4_factorial2 ( m - 1 ) ) * pi_sqrt 
     &    / 2.0D+00**( m / 2 )

      end if

      return
      end
      subroutine p06_fun ( option, n, x, f )

c*********************************************************************72
c
cc P06_FUN evaluates the integrand for problem 6.
c
c  Discussion:
c
c    The exact value is (m-1)cc * sqrt ( pi ) / sqrt ( 2^m ).
c
c    Integral ( -oo .lt. x .lt. +oo ) x^m exp (-x*x) dx
c
c    The parameter M is set by calling P06_PARAM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer m
      integer option
      double precision x(n)

      call p06_param ( 'G', 'M', m )

      do i = 1, n
        f(i) = x(i)**m
      end do

      if ( option .eq. 0 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - x(i) * x(i) )
        end do
      else if ( option .eq. 1 ) then

      else if ( option .eq. 2 ) then
        do i = 1, n 
          f(i) = f(i) * exp ( - 0.5D+00 * x(i) * x(i) )
        end do
      end if

      return
      end
      subroutine p06_param ( action, name, value )

c*********************************************************************72
c
cc P06_PARAM gets or sets parameters for problem 6.
c
c  Discussion:
c
c    The parameter is named "M", and it represents the value of the exponent
c    in the integrand function:
c
c    Integral ( -oo .lt. x .lt. +oo ) x^m exp (-x*x) dx
c
c    M must be greater than -1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ACTION, the action.
c    'S' to set the value,
c    'G' to get the value.
c
c    Input, character NAME, the parameter name.
c    'M', the exponent.
c
c    Input/output, integer VALUE, the parameter value.
c    If ACTION = 'S', then VALUE is an input quantity, and M is set to VALUE.
c    If ACTION = 'G', then VALUE is an output quantity, and VALUE is set to M.
c
      implicit none

      character action
      character name
      integer m
      integer value

      save m

      data m / 0 /

      if ( action .eq. 'S' .or. action .eq. 's' ) then

        if ( value <= -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Parameter M must be greater than -1.'
          stop
        end if

        m = value

      else if ( action .eq. 'G' .or. action .eq. 'g' ) then

        value = m

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized value of ACTION = "' 
     &    // action // '".'
        stop

      end if

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'x^m exp(-x*x)'

      return
      end
      subroutine p07_exact ( exact )

c*********************************************************************72
c
cc P07_EXACT returns the exact integral for problem 7.
c
c  Discussion:
c
c    The 20 digit values of pi^(1/2) and e^(1/4) were computed by Mathematica.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision e_sqrt_sqrt
      parameter ( e_sqrt_sqrt = 1.2840254166877414841D+00 )
      double precision exact
      double precision pi_sqrt 
      parameter ( pi_sqrt = 1.7724538509055160273D+00 )

      exact = 0.25D+00 * pi_sqrt / e_sqrt_sqrt

      return
      end
      subroutine p07_fun ( option, n, x, f )

c*********************************************************************72
c
cc P07_FUN evaluates the integrand for problem 7.
c
c  Discussion:
c
c    The exact value is (1/4) sqrt(pi) / sqrt(sqrt(e)).
c
c    Integral ( -oo .lt. x .lt. +oo ) x^2 cos(x) e^(-x^2) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Prem Kythe, Michael Schaeferkotter,
c    Handbook of Computational Methods for Integration,
c    Chapman and Hall, 2004,
c    ISBN: 1-58488-428-2,
c    LC: QA299.3.K98.
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision x(n)

      if ( option .eq. 0 ) then
        do i = 1, n
          f(i) = x(i)**2 * cos ( x(i) ) * exp ( - x(i)**2 )
        end do
      else if ( option .eq. 1 ) then
        do i = 1, n
          f(i) = x(i)**2 * cos ( x(i) )
        end do
      else if ( option .eq. 2 ) then
        do i = 1, n
          f(i) = x(i)**2 * cos ( x(i) ) 
     &      * exp ( - x(i)**2 / 2.0D+00 )
        end do
      end if

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns the title for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'x^2 cos ( x ) exp(-x*x)'

      return
      end
      subroutine p08_exact ( exact )

c*********************************************************************72
c
cc P08_EXACT returns the exact integral for problem 8.
c
c  Discussion:
c
c    The 20 digit value of the answer was computed by Mathematica.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 3.0088235661136433510D+00

      return
      end
      subroutine p08_fun ( option, n, x, f )

c*********************************************************************72
c
cc P08_FUN evaluates the integrand for problem 8.
c
c  Discussion:
c
c    The exact value is sqrt ( 2 pi ) * HypergeometricU ( -1/2, 0, 1 ).
c
c    Integral ( -oo .lt. x .lt. +oo ) sqrt(1+x*x/2) * exp(-x*x/2) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Prem Kythe, Michael Schaeferkotter,
c    Handbook of Computational Methods for Integration,
c    Chapman and Hall, 2004,
c    ISBN: 1-58488-428-2,
c    LC: QA299.3.K98.
c
c  Parameters:
c
c    Input, integer OPTION:
c    0, integrand is f(x).
c    1, integrand is exp(-x*x) * f(x);
c    2, integrand is exp(-x*x/2) * f(x);
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer option
      double precision x(n)

      do i = 1, n
        f(i) = sqrt ( 1.0D+00 + 0.5D+00 * x(i)**2 )
      end do

      if ( option .eq. 0 ) then
        do i = 1, n
          f(i) = f(i) * exp ( - 0.5D+00 * x(i)**2 )
        end do
      else if ( option .eq. 1 ) then
        do i = 1, n
          f(i) = f(i) * exp ( + 0.5D+00 * x(i)**2 )
        end do
      else if ( option .eq. 2 ) then

      end if

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns the title for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'sqrt(1+x*x/2) * exp(-x*x/2)'

      return
      end
      function r8_csevl ( x, a, n )

c*********************************************************************72
c
cc R8_CSEVL evaluates a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Roger Broucke.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, double precision CS(N), the Chebyshev coefficients.
c
c    Input, integer N, the number of Chebyshev coefficients.
c
c    Output, double precision R8_CSEVL, the Chebyshev series evaluated at X.
c
      implicit none

      integer n

      double precision a(n)
      double precision b0
      double precision b1
      double precision b2
      integer i
      double precision r8_csevl
      double precision twox
      double precision x

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms <= 0.'
        stop
      end if

      if ( 1000 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms > 1000.'
        stop
      end if

      if ( x .lt. -1.1D+00 .or. 1.1D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  X outside (-1,+1)'
        write ( *, '(a,g14.6)' ) '  X = ', x
        stop
      end if

      twox = 2.0D+00 * x
      b1 = 0.0D+00
      b0 = 0.0D+00

      do i = n, 1, -1
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + a(i)
      end do

      r8_csevl = 0.5D+00 * ( b0 - b2 )

      return
      end
      subroutine r8_gaml ( xmin, xmax )

c*********************************************************************72
c
cc R8_GAML evaluates bounds for an R8 argument of the gamma function.
c
c  Discussion:
c
c    This function calculates the minimum and maximum legal bounds 
c    for X in the evaluation of GAMMA ( X ).
c
c    XMIN and XMAX are not the only bounds, but they are the only 
c    non-trivial ones to calculate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Output, double precision XMIN, XMAX, the bounds.
c
      implicit none

      double precision alnbig
      double precision alnsml
      integer i
      integer j
      double precision r8_mach
      double precision xln
      double precision xmax
      double precision xmin
      double precision xold

      alnsml = dlog ( r8_mach ( 1 ) )
      xmin = - alnsml

      do i = 1, 10

        xold = xmin
        xln = dlog ( xmin )
        xmin = xmin - xmin * ( ( xmin + 0.5D+00 ) * xln - xmin 
     &    - 0.2258D+00 + alnsml ) / ( xmin * xln + 0.5D+00 )

        if ( dabs ( xmin - xold ) .lt. 0.005D+00 ) then

          xmin = - xmin + 0.01D+00

          alnbig = dlog ( r8_mach ( 2 ) )
          xmax = alnbig

          do j = 1, 10

            xold = xmax
            xln = dlog ( xmax )
            xmax = xmax - xmax * ( ( xmax - 0.5D+00 ) * xln - xmax 
     &        + 0.9189D+00 - alnbig ) / ( xmax * xln - 0.5D+00 )

            if ( dabs ( xmax - xold ) .lt. 0.005D+00 ) then
              xmax = xmax - 0.01D+00
              xmin = dmax1 ( xmin, - xmax + 1.0D+00 )
              return
            end if

          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAML - Fatal error!'
          write ( *, '(a)' ) '  Unable to find XMAX.'
          stop

        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMIN.'

      stop
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates the gamma function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_GAMMA, the gamma function of X.
c
      implicit none

      double precision dxrel
      double precision gcs(42)
      integer i
      integer n
      integer ngcs
      double precision pi
      double precision r8_csevl
      double precision r8_gamma
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision sinpiy
      double precision sq2pil
      double precision x
      double precision xmax
      double precision xmin
      double precision xsml
      double precision y

      save dxrel
      save gcs
      save ngcs
      save pi
      save sq2pil
      save xmax
      save xmin
      save xsml

      data gcs(  1) / +0.8571195590989331421920062399942D-02 /
      data gcs(  2) / +0.4415381324841006757191315771652D-02 /
      data gcs(  3) / +0.5685043681599363378632664588789D-01 /
      data gcs(  4) / -0.4219835396418560501012500186624D-02 /
      data gcs(  5) / +0.1326808181212460220584006796352D-02 /
      data gcs(  6) / -0.1893024529798880432523947023886D-03 /
      data gcs(  7) / +0.3606925327441245256578082217225D-04 /
      data gcs(  8) / -0.6056761904460864218485548290365D-05 /
      data gcs(  9) / +0.1055829546302283344731823509093D-05 /
      data gcs( 10) / -0.1811967365542384048291855891166D-06 /
      data gcs( 11) / +0.3117724964715322277790254593169D-07 /
      data gcs( 12) / -0.5354219639019687140874081024347D-08 /
      data gcs( 13) / +0.9193275519859588946887786825940D-09 /
      data gcs( 14) / -0.1577941280288339761767423273953D-09 /
      data gcs( 15) / +0.2707980622934954543266540433089D-10 /
      data gcs( 16) / -0.4646818653825730144081661058933D-11 /
      data gcs( 17) / +0.7973350192007419656460767175359D-12 /
      data gcs( 18) / -0.1368078209830916025799499172309D-12 /
      data gcs( 19) / +0.2347319486563800657233471771688D-13 /
      data gcs( 20) / -0.4027432614949066932766570534699D-14 /
      data gcs( 21) / +0.6910051747372100912138336975257D-15 /
      data gcs( 22) / -0.1185584500221992907052387126192D-15 /
      data gcs( 23) / +0.2034148542496373955201026051932D-16 /
      data gcs( 24) / -0.3490054341717405849274012949108D-17 /
      data gcs( 25) / +0.5987993856485305567135051066026D-18 /
      data gcs( 26) / -0.1027378057872228074490069778431D-18 /
      data gcs( 27) / +0.1762702816060529824942759660748D-19 /
      data gcs( 28) / -0.3024320653735306260958772112042D-20 /
      data gcs( 29) / +0.5188914660218397839717833550506D-21 /
      data gcs( 30) / -0.8902770842456576692449251601066D-22 /
      data gcs( 31) / +0.1527474068493342602274596891306D-22 /
      data gcs( 32) / -0.2620731256187362900257328332799D-23 /
      data gcs( 33) / +0.4496464047830538670331046570666D-24 /
      data gcs( 34) / -0.7714712731336877911703901525333D-25 /
      data gcs( 35) / +0.1323635453126044036486572714666D-25 /
      data gcs( 36) / -0.2270999412942928816702313813333D-26 /
      data gcs( 37) / +0.3896418998003991449320816639999D-27 /
      data gcs( 38) / -0.6685198115125953327792127999999D-28 /
      data gcs( 39) / +0.1146998663140024384347613866666D-28 /
      data gcs( 40) / -0.1967938586345134677295103999999D-29 /
      data gcs( 41) / +0.3376448816585338090334890666666D-30 /
      data gcs( 42) / -0.5793070335782135784625493333333D-31 /

      data dxrel / 0.0D+00 /
      data ngcs / 0 /
      data pi / 3.14159265358979323846264338327950D+00 /
      data sq2pil / 0.91893853320467274178032973640562D+00 /
      data xmax / 0.0D+00 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ngcs .eq. 0 ) then
        ngcs = r8_inits ( gcs, 42, 0.1D+00 * r8_mach ( 3 ) )
        call r8_gaml ( xmin, xmax )
        xsml = dexp ( dmax1 ( dlog ( r8_mach ( 1 ) ),
     &    - dlog ( r8_mach ( 2 ) ) ) + 0.01D+00 )
        dxrel = dsqrt ( r8_mach ( 4 ) )
      end if

      y = dabs ( x )

      if ( y .le. 10.0D+00 ) then

        n = int ( x )
        if ( x .lt. 0.0D+00 ) then
          n = n - 1
        end if
        y = x - dble ( n )
        n = n - 1
        r8_gamma = 0.9375D+00 + r8_csevl ( 2.0D+00 * y - 1.0D+00, 
     &    gcs, ngcs )

        if ( n .eq. 0 ) then

          return

        else if ( n .lt. 0 ) then

          n = - n

          if ( x .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is 0.'
            stop
          end if

          if ( x .lt. 0.0D+00 .and. 
     &      x + dble ( n - 2 ) .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is a negative integer.'
            stop
          end if

          if ( x .lt. - 0.5D+00 .and. 
     &      dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Warning!'
            write ( *, '(a)' ) '  X too near a negative integer,'
            write ( *, '(a)' ) '  answer is half precision.'
          end if

          if ( y .lt. xsml ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) 
     &        '  X is so close to zero that Gamma overflows.'
            stop
          end if

          do i = 1, n
            r8_gamma = r8_gamma / ( x + dble ( i - 1 ) )
          end do

        else if ( n .eq. 0 ) then

        else

          do i = 1, n
            r8_gamma = ( y + dble ( i ) ) * r8_gamma
          end do

        end if

      else

        if ( xmax .lt. x ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X so big that Gamma overflows.'
          stop
        end if
c
c  Underflow.
c
        if ( x .lt. xmin ) then
          r8_gamma = 0.0D+00
          return
        end if

        r8_gamma = dexp ( ( y - 0.5D+00 ) * dlog ( y ) - y + sq2pil 
     &    + r8_lgmc ( y ) )

        if ( 0.0D+00 .lt. x ) then
          return
        end if

        if ( dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Warning!'
          write ( *, '(a)' ) '  X too near a negative integer,'
          write ( *, '(a)' ) '  answer is half precision.'
        end if

        sinpiy = dsin ( pi * y )

        if ( sinpiy .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X is a negative integer.'
          stop
        end if

        r8_gamma = - pi / ( y * sinpiy * r8_gamma )

      end if

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
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
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_inits ( dos, nos, eta )

c*********************************************************************72
c
cc R8_INITS initializes a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision DOS(NOS), the Chebyshev coefficients.
c
c    Input, integer NOS, the number of coefficients.
c
c    Input, double precision ETA, the desired accuracy.
c
c    Output, integer R8_INITS, the number of terms of the series needed
c    to ensure the requested accuracy.
c
      implicit none

      integer nos

      double precision dos(nos)
      double precision err 
      double precision eta
      integer i
      integer r8_inits

      if ( nos .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_INITS - Fatal error!'
        write ( *, '(a)' ) '  Number of coefficients < 1.'
        stop
      end if

      err = 0.0D+00

      do i = nos, 1, -1
        err = err + dabs ( dos(i) )
        if ( eta .lt. err ) then
          r8_inits = i
          return
        end if
      end do

      r8_inits = nos
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_INITS - Warning!'
      write ( *, '(a)' ) '  ETA may be too small.'

      return
      end
      function r8_lgmc ( x )

c*********************************************************************72
c
cc R8_LGMC evaluates the log gamma correction factor for an R8 argument.
c
c  Discussion:
c
c    For 10 <= X, compute the log gamma correction factor so that
c
c      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
c                          + ( x - 0.5 ) * log ( x ) - x 
c                          + r8_lgmc ( x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_LGMC, the correction factor.
c
      implicit none

      double precision algmcs(15)
      integer nalgm
      double precision r8_csevl
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision x
      double precision xbig
      double precision xmax

      save algmcs
      save nalgm
      save xbig
      save xmax

      data algmcs(  1) / +0.1666389480451863247205729650822D+00 /
      data algmcs(  2) / -0.1384948176067563840732986059135D-04 /
      data algmcs(  3) / +0.9810825646924729426157171547487D-08 /
      data algmcs(  4) / -0.1809129475572494194263306266719D-10 /
      data algmcs(  5) / +0.6221098041892605227126015543416D-13 /
      data algmcs(  6) / -0.3399615005417721944303330599666D-15 /
      data algmcs(  7) / +0.2683181998482698748957538846666D-17 /
      data algmcs(  8) / -0.2868042435334643284144622399999D-19 /
      data algmcs(  9) / +0.3962837061046434803679306666666D-21 /
      data algmcs( 10) / -0.6831888753985766870111999999999D-23 /
      data algmcs( 11) / +0.1429227355942498147573333333333D-24 /
      data algmcs( 12) / -0.3547598158101070547199999999999D-26 /
      data algmcs( 13) / +0.1025680058010470912000000000000D-27 /
      data algmcs( 14) / -0.3401102254316748799999999999999D-29 /
      data algmcs( 15) / +0.1276642195630062933333333333333D-30 /

      data nalgm / 0 /
      data xbig / 0.0D+00 /
      data xmax / 0.0D+00 /

      if ( nalgm .eq. 0 ) then
        nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) )
        xbig = 1.0D+00 / dsqrt ( r8_mach ( 3 ) )
        xmax = dexp ( dmin1 ( dlog ( r8_mach ( 2 ) / 12.0D+00 ), 
     &    - dlog ( 12.0D+00 * r8_mach ( 1 ) ) ) )
      end if

      if ( x .lt. 10.0D+00 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_LGMC - Fatal error!'
        write ( *, '(a)' ) '  X must be at least 10.'
        stop

      else if ( x .lt. xbig ) then

        r8_lgmc = r8_csevl ( 2.0D+00 * ( 10.0D+00 / x ) 
     &    * ( 10.0D+00 / x ) - 1.0D+00, algmcs, nalgm ) / x

      else if ( x .lt. xmax ) then

        r8_lgmc = 1.0D+00 / ( 12.0D+00 * x )

      else

        r8_lgmc = 0.0D+00

      end if

      return
      end
      function r8_mach ( i )

c*********************************************************************72
c
cc R8_MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    R8_MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = R8_MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    R8_MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
c    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
c    R8_MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision R8_MACH, the value of the constant.
c
      implicit none

      double precision r8_mach
      integer i

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      else if ( i .eq. 1 ) then
        r8_mach = 4.450147717014403D-308
      else if ( i .eq. 2 ) then
        r8_mach = 8.988465674311579D+307
      else if ( i .eq. 3 ) then
        r8_mach = 1.110223024625157D-016
      else if ( i .eq. 4 ) then
        r8_mach = 2.220446049250313D-016
      else if ( i .eq. 5 ) then
        r8_mach = 0.301029995663981D+000
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      end if

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

      double precision r8_uniform_01
      integer k
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
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
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
c    An R8VEC is a vector of R8's.
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
