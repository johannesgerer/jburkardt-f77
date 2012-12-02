      subroutine laguerre_compute ( order, xtab, weight, alpha )

c*********************************************************************72
c
cc LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
c
c  Discussion:
c
c    In the simplest case, ALPHA is 0, and we are approximating the
c    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
c    it is easy to modify the rule to approximate the integral from
c    A to +oo as well.
c
c    If ALPHA is nonzero, then there is no simple way to extend the
c    rule to approximate the integral from A to +oo.  The simplest
c    procedures would be to approximate the integral from 0 to A.
c
c    The integration interval is [ A, +oo ) or [ 0, +oo ).
c
c    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
c
c
c    If the integral to approximate is:
c
c        Integral ( A .le. X < +oo ) exp ( - X ) * F(X) dX
c      or
c        Integral ( 0 .le. X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX
c
c    then the quadrature rule is:
c
c      exp ( - A ) * Sum ( 1 .le. I .le. ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
c    or
c      Sum ( 1 .le. I .le. ORDER ) WEIGHT(I) * F ( XTAB(I) )
c
c
c    If the integral to approximate is:
c
c        Integral ( A .le. X < +oo ) F(X) dX
c      or
c        Integral ( 0 .le. X < +oo ) X^ALPHA * F(X) dX
c
c    then the quadrature rule is:
c
c      exp ( - A ) * Sum ( 1 .le. I .le. ORDER ) 
c        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
c    or
c      Sum ( 1 .le. I .le. ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 March 2000
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
c    Input, integer ORDER, the order of the quadrature rule 
c    to be computed.  ORDER must be at least 1.
c
c    Output, double precision XTAB(ORDER), the abscissas.
c
c    Output, double precision WEIGHT(ORDER), the weights.
c
c    Input, double precision ALPHA, the exponent of the X factor.
c    Set ALPHA = 0.0D+00 for the simplest rule.
c    ALPHA must be nonnegative.
c
      implicit none

      integer order

      double precision alpha
      double precision b(order)
      double precision c(order)
      double precision cc
      double precision dp2
      integer i
      double precision p1
      double precision r1
      double precision r2
      double precision r8_gamma
      double precision ratio
      double precision weight(order)
      double precision x
      double precision xtab(order)
c
c  Set the recursion coefficients.
c
      do i = 1, order
        b(i) = ( alpha + dble ( 2 * i - 1 ) )
      end do

      do i = 1, order
        c(i) = dble ( i - 1 ) * ( alpha + dble ( i - 1 ) )
      end do

      cc = 1.0D+00
      do i = 2, order
        cc = cc * c(i)
      end do
      cc = r8_gamma ( alpha + 1.0D+00 ) * cc

      do i = 1, order
c
c  Compute an estimate for the root.
c
        if ( i .eq. 1 ) then

          x = ( 1.0D+00 + alpha ) * ( 3.0D+00+ 0.92 * alpha ) / 
     &      ( 1.0D+00 + 2.4D+00 * dble ( order ) + 1.8D+00 * alpha )

        else if ( i .eq. 2 ) then

          x = x + ( 15.0D+00 + 6.25D+00 * alpha ) / 
     &      ( 1.0D+00 + 0.9D+00 * alpha + 2.5D+00 * dble ( order ) )

        else

          r1 = ( 1.0D+00 + 2.55D+00 * dble ( i - 2 ) ) 
     &      / ( 1.9D+00 * dble ( i - 2 ) )

          r2 = 1.26D+00 * dble ( i - 2 ) * alpha / 
     &      ( 1.0D+00 + 3.5D+00 * dble ( i - 2 ) )

          ratio = ( r1 + r2 ) / ( 1.0D+00 + 0.3D+00 * alpha )

          x = x + ratio * ( x - xtab(i-2) )

        end if
c
c  Use iteration to find the root.
c
        call laguerre_root ( x, order, alpha, dp2, p1, b, c )
c
c  Set the abscissa and weight.
c
        xtab(i) = x
        weight(i) = ( cc / dp2 ) / p1

      end do

      return
      end
      subroutine laguerre_recur ( p2, dp2, p1, x, order, alpha, b, c )

c*********************************************************************72
c
cc LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 1998
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
c    Output, double precision P2, the value of L(ORDER)(X).
c
c    Output, double precision DP2, the value of L'(ORDER)(X).
c
c    Output, double precision P1, the value of L(ORDER-1)(X).
c
c    Input, double precision X, the point at which polynomials are evaluated.
c
c    Input, integer ORDER, the order of the polynomial 
c    to be computed.
c
c    Input, double precision ALPHA, the exponent of the X factor in the
c    integrand.
c
c    Input, double precision B(ORDER), C(ORDER), the recursion
c    coefficients.
c
      implicit none

      integer order

      double precision alpha
      double precision b(order)
      double precision c(order)
      double precision dp0
      double precision dp1
      double precision dp2
      integer i
      double precision p0
      double precision p1
      double precision p2
      double precision x

      p1 = 1.0D+00
      dp1 = 0.0D+00

      p2 = x - alpha - 1.0D+00
      dp2 = 1.0D+00

      do i = 2, order

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2 = ( x - b(i) ) * p1 - c(i) * p0
        dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

      end do

      return
      end
      subroutine laguerre_root ( x, order, alpha, dp2, p1, b, c )

c*********************************************************************72
c
cc LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2000
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    FORTRAN90 version by John Burkardt.
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
c    Input, integer ORDER, the order of the polynomial 
c    to be computed.
c
c    Input, double precision ALPHA, the exponent of the X factor.
c
c    Output, double precision DP2, the value of L'(ORDER)(X).
c
c    Output, double precision P1, the value of L(ORDER-1)(X).
c
c    Input, double precision B(ORDER), C(ORDER), the recursion coefficients.
c
      implicit none

      integer order

      double precision alpha
      double precision b(order)
      double precision c(order)
      double precision d
      double precision dp2
      double precision eps
      double precision p1
      double precision p2
      double precision r8_epsilon
      integer step
      integer step_max
      parameter ( step_max = 10 )
      double precision x

      eps = r8_epsilon ( )

      do step = 1, step_max

        call laguerre_recur ( p2, dp2, p1, x, order, alpha, b, c )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0D+00 ) ) then
          return
        end if

      end do

      return
      end
      subroutine legendre_compute ( n, x, w )

c*********************************************************************72
c
cc LEGENDRE_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.
c
c  Discussion:
c
c    The integral:
c
c      integral ( -1 .le. x .le. 1 ) f(x) dx
c
c    The quadrature rule:
c
c      sum ( 1 .le. i .le. n ) w(i) * f ( x(i) )
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
c    FORTRAN90 version by John Burkardt.
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
c    Output, double precision X(N), the abscissas.
c
c    Output, double precision W(N), the weights.
c
      implicit none

      integer n

      double precision d1
      double precision d2pn
      double precision d3pn
      double precision d4pn
      double precision dp
      double precision dpn
      double precision e1
      double precision fx
      double precision h
      integer i
      integer iback
      integer k
      integer m
      integer mp1mi
      integer ncopy
      integer nmove
      double precision p
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision pk
      double precision pkm1
      double precision pkp1
      double precision t
      double precision u
      double precision v
      double precision w(n)
      double precision x(n)
      double precision x0
      double precision xtemp

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_DR_COMPUTE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of N = ', n
        stop
      end if

      e1 = dble ( n * ( n + 1 ) )

      m = ( n + 1 ) / 2

      do i = 1, m

        mp1mi = m + 1 - i

        t = dble ( 4 * i - 1 ) * pi 
     &    / dble ( 4 * n + 2 )

        x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 
     &    / dble ( n ) ) 
     &    / dble ( 8 * n * n ) )

        pkm1 = 1.0D+00
        pk = x0

        do k = 2, n
          pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) 
     &      / dble ( k )
          pkm1 = pk
          pk = pkp1
        end do

        d1 = dble ( n ) * ( pkm1 - x0 * pk )

        dpn = d1 / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

        d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 ) 
     &    / ( 1.0D+00 + x0 )

        d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) 
     &    / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

        d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) 
     &    / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

        u = pk / dpn
        v = d2pn / dpn
c
c  Initial approximation H:
c
        h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / 
     &    ( 3.0D+00 * dpn ) ) ) )
c
c  Refine H using one step of Newton's method:
c
        p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 
     &    * ( d3pn + 0.25D+00 * h * d4pn ) ) )

        dp = dpn + h * ( d2pn + 0.5D+00 * h * 
     &    ( d3pn + h * d4pn / 3.0D+00 ) )

        h = h - p / dp

        xtemp = x0 + h

        x(mp1mi) = xtemp

        fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 
     &    * ( d2pn + 0.25D+00 * h 
     &    * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

        w(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp ) 
     &    * ( 1.0D+00 + xtemp ) / ( fx * fx )

      end do

      if ( mod ( n, 2 ) .eq. 1 ) then
        x(1) = 0.0D+00
      end if
c
c  Shift the data up.
c
      nmove = ( n + 1 ) / 2
      ncopy = n - nmove

      do i = 1, nmove
        iback = n + 1 - i
        x(iback) = x(iback-ncopy)
        w(iback) = w(iback-ncopy)
      end do
c
c  Reflect values for the negative abscissas.
c
      do i = 1, n - nmove
        x(i) = - x(n+1-i)
        w(i) = w(n+1-i)
      end do

      return
      end
      subroutine p00_alpha ( problem, alpha )

c*********************************************************************72
c
cc P00_ALPHA returns the value of ALPHA for any problem.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c    The typical or default value is 0.0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha
      integer problem

      if ( problem .eq. 1 ) then
        call p01_alpha ( alpha )
      else if ( problem .eq. 2 ) then
        call p02_alpha ( alpha )
      else if ( problem .eq. 3 ) then
        call p03_alpha ( alpha )
      else if ( problem .eq. 4 ) then
        call p04_alpha ( alpha )
      else if ( problem .eq. 5 ) then
        call p05_alpha ( alpha )
      else if ( problem .eq. 6 ) then
        call p06_alpha ( alpha )
      else if ( problem .eq. 7 ) then
        call p07_alpha ( alpha )
      else if ( problem .eq. 8 ) then
        call p08_alpha ( alpha )
      else if ( problem .eq. 9 ) then
        call p09_alpha ( alpha )
      else if ( problem .eq. 10 ) then
        call p10_alpha ( alpha )
      else if ( problem .eq. 11 ) then
        call p11_alpha ( alpha )
      else if ( problem .eq. 12 ) then
        call p12_alpha ( alpha )
      else if ( problem .eq. 13 ) then
        call p13_alpha ( alpha )
      else if ( problem .eq. 14 ) then
        call p14_alpha ( alpha )
      else if ( problem .eq. 15 ) then
        call p15_alpha ( alpha )
      else if ( problem .eq. 16 ) then
        call p16_alpha ( alpha )
      else if ( problem .eq. 17 ) then
        call p17_alpha ( alpha )
      else if ( problem .eq. 18 ) then
        call p18_alpha ( alpha )
      else if ( problem .eq. 19 ) then
        call p19_alpha ( alpha )
      else if ( problem .eq. 20 ) then
        call p20_alpha ( alpha )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_ALPHA - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
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
c    26 December 2011
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
      else if ( problem .eq. 9 ) then
        call p09_exact ( exact )
      else if ( problem .eq. 10 ) then
        call p10_exact ( exact )
      else if ( problem .eq. 11 ) then
        call p11_exact ( exact )
      else if ( problem .eq. 12 ) then
        call p12_exact ( exact )
      else if ( problem .eq. 13 ) then
        call p13_exact ( exact )
      else if ( problem .eq. 14 ) then
        call p14_exact ( exact )
      else if ( problem .eq. 15 ) then
        call p15_exact ( exact )
      else if ( problem .eq. 16 ) then
        call p16_exact ( exact )
      else if ( problem .eq. 17 ) then
        call p17_exact ( exact )
      else if ( problem .eq. 18 ) then
        call p18_exact ( exact )
      else if ( problem .eq. 19 ) then
        call p19_exact ( exact )
      else if ( problem .eq. 20 ) then
        call p20_exact ( exact )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_exp_transform ( problem, order, result )

c*********************************************************************72
c
cc P00_EXP_TRANSFORM applies an exponential transform and Gauss-Legendre rule.
c
c  Discussion:
c
c    To approximate:
c
c      Integral ( alpha .le. x < +oo ) f(x) dx
c
c    Transform:
c
c      u = exp ( -x )
c      du = - exp ( -x ) dx
c
c      x = - log ( u )
c      dx = - du / u
c
c      x = alpha    => u = exp ( -alpha )
c      x = +oo      => u = 0
c
c    Transformed integral:
c
c      Integral ( 0 < u .le. exp ( -alpha ) ) f ( -log(u) ) du / u
c
c    We apply a Gauss-Legendre rule here, but we could easily use any rule
c    that avoids evaluation at U = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2007
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
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, integer ORDER, the order of the Gauss-Legendre rule 
c    to apply.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision alpha
      double precision f_vec(order)
      integer i
      integer order
      integer problem
      double precision r8vec_dot_product
      double precision result
      double precision u(order)
      double precision u_log(order)
      double precision weight(order)

      call p00_alpha ( problem, alpha )
c
c  Get the abscissas and weights for Gauss-Legendre quadrature.
c
      call legendre_compute ( order, u, weight )
c
c  Modify the weights from [-1,1] to [0,exp(-alpha)].
c
      do i = 1, order
        weight(i) = exp ( -alpha ) * weight(i) / 2.0D+00
      end do
c
c  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
c
      do i = 1, order
        u(i) = ( ( 1.0D+00 + u(i) ) * exp ( - alpha ) 
     &         + ( 1.0D+00 - u(i) ) * 0.0D+00 )       
     &         / ( 2.0D+00              )
      end do
c
c  Define U_LOG = - log ( U )
c
      do i = 1, order
        u_log(i) = - log ( u(i) )
      end do
c
c  Evaluate F ( -LOG(U) ).
c
      call p00_fun ( problem, order, u_log, f_vec )
c
c  The integrand is F ( -LOG(U) ) / U
c
      do i = 1, order
        f_vec(i) = f_vec(i) / u(i)
      end do
c
c  Sum.
c
      result = r8vec_dot_product ( order, weight, f_vec )

      return
      end
      subroutine p00_fun ( problem, n, x, f )

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
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the index of the problem.
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
      integer problem
      double precision x(n)

      if ( problem .eq. 1 ) then
        call p01_fun ( n, x, f )
      else if ( problem .eq. 2 ) then
        call p02_fun ( n, x, f )
      else if ( problem .eq. 3 ) then
        call p03_fun ( n, x, f )
      else if ( problem .eq. 4 ) then
        call p04_fun ( n, x, f )
      else if ( problem .eq. 5 ) then
        call p05_fun ( n, x, f )
      else if ( problem .eq. 6 ) then
        call p06_fun ( n, x, f )
      else if ( problem .eq. 7 ) then
        call p07_fun ( n, x, f )
      else if ( problem .eq. 8 ) then
        call p08_fun ( n, x, f )
      else if ( problem .eq. 9 ) then
        call p09_fun ( n, x, f )
      else if ( problem .eq. 10 ) then
        call p10_fun ( n, x, f )
      else if ( problem .eq. 11 ) then
        call p11_fun ( n, x, f )
      else if ( problem .eq. 12 ) then
        call p12_fun ( n, x, f )
      else if ( problem .eq. 13 ) then
        call p13_fun ( n, x, f )
      else if ( problem .eq. 14 ) then
        call p14_fun ( n, x, f )
      else if ( problem .eq. 15 ) then
        call p15_fun ( n, x, f )
      else if ( problem .eq. 16 ) then
        call p16_fun ( n, x, f )
      else if ( problem .eq. 17 ) then
        call p17_fun ( n, x, f )
      else if ( problem .eq. 18 ) then
        call p18_fun ( n, x, f )
      else if ( problem .eq. 19 ) then
        call p19_fun ( n, x, f )
      else if ( problem .eq. 20 ) then
        call p20_fun ( n, x, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FUN - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_gauss_laguerre ( problem, order, result )

c*********************************************************************72
c
cc P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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

      double precision alpha
      double precision alpha2
      double precision f_vec(order)
      integer i
      integer order
      integer problem
      double precision r8vec_dot_product
      double precision result
      double precision weight(order)
      double precision xtab(order)

      call p00_alpha ( problem, alpha )

      alpha2 = 0.0D+00
      call laguerre_compute ( order, xtab, weight, alpha2 )

      do i = 1, order
        xtab(i) = xtab(i) + alpha
      end do

      call p00_fun ( problem, order, xtab, f_vec )
c
c  The Gauss-Laguerre rule assumes a weight of EXP(-X).
c
c  We need to multiply each F(X) by EXP(X) to implicitly 
c  adjust for this weight.
c
      do i = 1, order
        f_vec(i) = f_vec(i) * exp ( xtab(i) )
      end do

      result = exp ( -alpha ) 
     &  * r8vec_dot_product ( order, weight, f_vec )

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
c    26 December 2011
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

      problem_num = 20

      return
      end
      subroutine p00_rat_transform ( problem, order, result )

c*********************************************************************72
c
cc P00_RAT_TRANSFORM applies a rational transform and Gauss-Legendre rule.
c
c  Discussion:
c
c    To approximate:
c
c      Integral ( alpha .le. x < +oo ) f(x) dx
c
c    Transform:
c
c      u = 1 / ( 1 + x )
c      du = - dx / ( 1 + x )^2
c
c      x = ( 1 - u ) / u
c      dx = - du / u^2
c
c      x = alpha    => u = 1 / ( 1 + alpha )
c      x = +oo      => u = 0
c
c    Transformed integral:
c
c      Integral ( 0 < u .le. 1 / ( 1 + alpha ) ) f ( ( 1 - u ) / u ) du / u^2
c
c    We apply a Gauss-Legendre rule here, but we could easily use any rule
c    that avoids evaluation at U = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2007
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
c    Input, integer PROBLEM, the index of the problem.
c
c    Input, integer ORDER, the order of the Gauss-Legendre rule 
c    to apply.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision alpha
      double precision f_vec(order)
      integer i
      integer order
      integer problem
      double precision r8vec_dot_product
      double precision result
      double precision u(order)
      double precision u_rat(order)
      double precision weight(order)

      call p00_alpha ( problem, alpha )
c
c  Get the abscissas and weights for Gauss-Legendre quadrature.
c
      call legendre_compute ( order, u, weight )
c
c  Modify the weights from [-1,1] to [0,1/(1+alpha)].
c
      do i = 1, order
        weight(i) = weight(i) / 2.0D+00 / ( 1.0D+00 + alpha )
      end do
c
c  Linear transform of abscissas from [-1,1] to [0,1/(1+alpha)].
c
      do i = 1, order
        u(i) = ( ( 1.0D+00 + u(i) ) / ( 1.0D+00 + alpha ) 
     &         + ( 1.0D+00 - u(i) ) * 0.0D+00 )       
     &         / ( 2.0D+00              )
      end do
c
c  Define U_RAT = ( 1 - U ) / U.
c
      do i = 1, order
        u_rat(i) = ( 1.0D+00 - u(i) ) / u(i)
      end do
c
c  Evaluate F ( ( 1 - U ) / U ).
c
      call p00_fun ( problem, order, u_rat, f_vec )
c
c  The integrand is F ( ( 1 - U ) / U ) / U^2
c
      do i = 1, order
        f_vec(i) = f_vec(i) / u(i)**2
      end do
c
c  Sum.
c
      result = r8vec_dot_product ( order, weight, f_vec )

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
c    26 December 2011
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
      else if ( problem .eq. 9 ) then
        call p09_title ( title )
      else if ( problem .eq. 10 ) then
        call p10_title ( title )
      else if ( problem .eq. 11 ) then
        call p11_title ( title )
      else if ( problem .eq. 12 ) then
        call p12_title ( title )
      else if ( problem .eq. 13 ) then
        call p13_title ( title )
      else if ( problem .eq. 14 ) then
        call p14_title ( title )
      else if ( problem .eq. 15 ) then
        call p15_title ( title )
      else if ( problem .eq. 16 ) then
        call p16_title ( title )
      else if ( problem .eq. 17 ) then
        call p17_title ( title )
      else if ( problem .eq. 18 ) then
        call p18_title ( title )
      else if ( problem .eq. 19 ) then
        call p19_title ( title )
      else if ( problem .eq. 20 ) then
        call p20_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p01_alpha ( alpha )

c*********************************************************************72
c
cc P01_ALPHA returns ALPHA for problem 1.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

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
c    27 July 2007
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

      exact = 0.19524754198276439152D+00

      return
      end
      subroutine p01_fun ( n, x, f )

c*********************************************************************72
c
cc P01_FUN evaluates the integrand for problem 1.
c
c  Discussion:
c
c    D&R gives "exact" value as 0.19524753...
c    Mathematica returns        0.19524754198276439152...
c    D&R gives Laguerre(16) as  0.16623627...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) dx / ( x * log(x)^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision x(n)

      do i = 1, n
        f(i) = exp ( -2.0D+00 ) / ( x(i) * log ( x(i) )**2 )
      end do

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
c    27 July 2007
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

      title = '1 / ( x * log(x)^2 )'

      return
      end
      subroutine p02_alpha ( alpha )

c*********************************************************************72
c
cc P02_ALPHA returns ALPHA for problem 2.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

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
c    27 July 2007
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

      exact = 0.32510848278991335198D+00

      return
      end
      subroutine p02_fun ( n, x, f )

c*********************************************************************72
c
cc P02_FUN evaluates the integrand for problem 2.
c
c  Discussion:
c
c    D&R gives "exact" value as 0.32510855...
c    Mathematica returns        0.32510848278991335198...
c    D&R gives Laguerre(16) as  0.19142399...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) dx / ( x * log(x)^(3/2) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision x(n)

      do i = 1, n
        f(i) = exp ( -2.0D+00 ) / ( x(i) * sqrt ( log ( x(i) )**3 ) )
      end do

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
c    27 July 2007
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

      title = '1 / ( x * log(x)^(3/2) )'

      return
      end
      subroutine p03_alpha ( alpha )

c*********************************************************************72
c
cc P03_ALPHA returns ALPHA for problem 3.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

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
c    27 July 2007
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

      exact = 13.628D+00

      return
      end
      subroutine p03_fun ( n, x, f )

c*********************************************************************72
c
cc P03_FUN evaluates the integrand for problem 3.
c
c  Discussion:
c
c    D&R gives "exact" value as 13.628...
c    Mathematica returns        13.440045415012575106...
c    D&R gives Laguerre(16) as   0.44996932...
c
c    This integral is "something of a numerical joke, as it is
c    scarcely distinguishable from the divergent integrand 1/x."
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) dx / ( x^1.01 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision x(n)

      do i = 1, n
        f(i) = exp ( -2.0D+00 ) * 1.0D+00 / x(i)**1.01D+00
      end do

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
c    27 July 2007
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

      title = '1 / ( x^1.01 )'

      return
      end
      subroutine p04_alpha ( alpha )

c*********************************************************************72
c
cc P04_ALPHA returns ALPHA for problem 4.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

      return
      end
      subroutine p04_exact ( exact )

c*********************************************************************72
c
cc P04_EXACT returns the estimated integral for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = -0.0046848541335080643181D+00

      return
      end
      subroutine p04_fun ( n, x, f )

c*********************************************************************72
c
cc P04_FUN evaluates the integrand for problem 4.
c
c  Discussion:
c
c    D&R gives "exact" value as -0.0046984...
c    Mathematica returns        -0.0046848541335080643181...
c    D&R gives Laguerre(16) as  -0.039258696...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) ( sin ( x ) / x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision x(n)

      do i = 1, n
        if ( x(i) .eq. 0.0D+00 ) then
          f(i) = exp ( -2.0D+00 )
        else
          f(i) = exp ( -2.0D+00 ) * sin ( x(i) ) / x(i)
        end if
      end do

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
c    27 July 2007
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

      title = 'Sine integral'

      return
      end
      subroutine p05_alpha ( alpha )

c*********************************************************************72
c
cc P05_ALPHA returns ALPHA for problem 5.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

      return
      end
      subroutine p05_exact ( exact )

c*********************************************************************72
c
cc P05_EXACT returns the estimated integral for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.0015897286158592328774D+00

      return
      end
      subroutine p05_fun ( n, x, f )

c*********************************************************************72
c
cc P05_FUN evaluates the integrand for problem 5.
c
c  Discussion:
c
c    D&R gives "exact" value as  0.00158973...
c    Mathematica returns         0.0015897286158592328774...
c    D&R gives Laguerre(16) as  -0.067859545...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) cos ( pi * x^2 / 2 ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )
      double precision x(n)

      do i = 1, n
        f(i) = exp ( -2.0D+00 ) * cos ( 0.5D+00 * pi * x(i)**2 )
      end do

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
c    27 July 2007
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

      title = 'Fresnel integral'

      return
      end
      subroutine p06_alpha ( alpha )

c*********************************************************************72
c
cc P06_ALPHA returns ALPHA for problem 6.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

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
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.00056103711148387120640D+00

      return
      end
      subroutine p06_fun ( n, x, f )

c*********************************************************************72
c
cc P06_FUN evaluates the integrand for problem 6.
c
c  Discussion:
c
c    D&R gives "exact" value as 0.0005610371...
c    Mathematica returns        0.00056103711148387120640...
c    D&R gives Laguerre(16) as  0.00056100775...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) exp ( -x^2 ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision exponent_min
      parameter ( exponent_min = -80.0D+00 )
      double precision f(n)
      integer i
      double precision x(n)

      do i = 1, n
        if ( - x(i)**2 .lt. exponent_min ) then
          f(i) = 0.0D+00
        else
          f(i) = exp ( -2.0D+00 ) * exp ( - x(i)**2 )
        end if
      end do

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
c    27 July 2007
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

      title = 'Complementary error function'

      return
      end
      subroutine p07_alpha ( alpha )

c*********************************************************************72
c
cc P07_ALPHA returns ALPHA for problem 7.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

      return
      end
      subroutine p07_exact ( exact )

c*********************************************************************72
c
cc P07_EXACT returns the exact integral for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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

      exact = 0.16266891D+00

      return
      end
      subroutine p07_fun ( n, x, f )

c*********************************************************************72
c
cc P07_FUN evaluates the integrand for problem 7.
c
c  Discussion:
c
c    D&R gives "exact" value as 0.16266891...
c    Mathematica does not return a value.
c    D&R gives Laguerre(16) as  0.097083064...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) sin ( x - 1 ) dx
c      / sqrt ( x * ( x - 2 ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 2.0D+00 ) then
          f(i) = 0.0D+00
        else
          f(i) = exp ( -2.0D+00 ) 
     &      * sin ( x(i) - 1.0D+00 ) 
     &      / sqrt ( x(i) * ( x(i) - 2.0D+00 ) )
        end if

      end do

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
c    27 July 2007
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

      title = 'Bessel function'

      return
      end
      subroutine p08_alpha ( alpha )

c*********************************************************************72
c
cc P08_ALPHA returns ALPHA for problem 8.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 2.0D+00

      return
      end
      subroutine p08_exact ( exact )

c*********************************************************************72
c
cc P08_EXACT returns the estimated integral for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.058334852497734677320D+00

      return
      end
      subroutine p08_fun ( n, x, f )

c*********************************************************************72
c
cc P08_FUN evaluates the integrand for problem 8.
c
c  Discussion:
c
c    D&R gives "exact" value as 0.0583349...
c    Mathematica returns        0.058334852497734677320...
c    D&R gives Laguerre(16) as  0.05834841...
c
c  Integral
c
c    exp ( -2 ) Integral ( 2 .le. x < +oo ) x / ( exp ( x ) - 1 ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
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
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision exponent_max
      parameter ( exponent_max = 80.0D+00 )
      double precision f(n)
      integer i
      double precision x(n)

      do i = 1, n
        if ( x(i) .lt. exponent_max ) then
          f(i) = exp ( -2.0D+00 ) * x(i) / ( exp ( x(i) ) - 1.0D+00 )
        else
          f(i) = 0.0D+00
        end if
      end do

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
c    27 July 2007
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

      title = 'Debye function'

      return
      end
      subroutine p09_alpha ( alpha )

c*********************************************************************72
c
cc P09_ALPHA returns ALPHA for problem 9.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p09_exact ( exact )

c*********************************************************************72
c
cc P09_EXACT returns the estimated integral for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 24.0D+00

      return
      end
      subroutine p09_fun ( n, x, f )

c*********************************************************************72
c
cc P09_FUN evaluates the integrand for problem 9.
c
c  Discussion:
c
c    The integral is the definition of the Gamma function for
c    Z = 5, with exact value (Z-1)! = 24.
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) x^4 exp ( -x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision exponent_min
      parameter ( exponent_min = -80.0D+00 )
      double precision f(n)
      integer i
      double precision x(n)

      do i = 1, n
        if ( -x(i) .lt. exponent_min ) then
          f(i) = 0.0D+00
        else
          f(i) = x(i)**4 * exp ( -x(i) )
        end if
      end do

      return
      end
      subroutine p09_title ( title )

c*********************************************************************72
c
cc P09_TITLE returns the title for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = 'Gamma(Z=5) function'

      return
      end
      subroutine p10_alpha ( alpha )

c*********************************************************************72
c
cc P10_ALPHA returns ALPHA for problem 10.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p10_exact ( exact )

c*********************************************************************72
c
cc P10_EXACT returns the estimated integral for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = pi / 2.0D+00

      return
      end
      subroutine p10_fun ( n, x, f )

c*********************************************************************72
c
cc P10_FUN evaluates the integrand for problem 10.
c
c  Discussion:
c
c    S&S gives exact value as pi/2 = 1.5707963267948966192...
c    S&S gives Laguerre(16) as       1.5537377347...
c    S&S gives EXP_TRANSFORM(16) as  1.4293043007...
c    S&S gives RAT_TRANSFORM(16) as  1.5707963267...
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) 1/(1+x*x) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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
      double precision x(n)

      do i = 1, n
        f(i) = 1.0D+00 / ( 1.0D+00 + x(i) * x(i) )
      end do

      return
      end
      subroutine p10_title ( title )

c*********************************************************************72
c
cc P10_TITLE returns the title for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = '1/(1+x*x)'

      return
      end
      subroutine p11_alpha ( alpha )

c*********************************************************************72
c
cc P11_ALPHA returns ALPHA for problem 11.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p11_exact ( exact )

c*********************************************************************72
c
cc P11_EXACT returns the estimated integral for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = pi

      return
      end
      subroutine p11_fun ( n, x, f )

c*********************************************************************72
c
cc P11_FUN evaluates the integrand for problem 11.
c
c  Discussion:
c
c    S&S gives exact value as pi =  3.1415926535897932385...
c    S&S gives Laguerre(16) as      2.6652685196...
c    S&S gives EXP_TRANSFORM(16) as 2.3629036166...
c    S&S gives RAT_TRANSFORM(16) as 3.0360705907... 
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) 1/((1+x)*sqrt(x)) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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
      double precision x(n)

      do i = 1, n
        if ( x(i) .eq. 0.0D+00 ) then
          f(i) = 0.0D+00
        else
          f(i) = 1.0D+00 / ( ( 1.0D+00 + x(i) ) * sqrt ( x(i) ) )
        end if
      end do

      return
      end
      subroutine p11_title ( title )

c*********************************************************************72
c
cc P11_TITLE returns the title for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = '1 / ( (1+x) * sqrt(x) )'

      return
      end
      subroutine p12_alpha ( alpha )

c*********************************************************************72
c
cc P12_ALPHA returns ALPHA for problem 12.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p12_exact ( exact )

c*********************************************************************72
c
cc P12_EXACT returns the estimated integral for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.5D+00

      return
      end
      subroutine p12_fun ( n, x, f )

c*********************************************************************72
c
cc P12_FUN evaluates the integrand for problem 12.
c
c  Discussion:
c
c    S&S gives exact value as pi =  0.5
c    S&S gives Laguerre(16) as      0.5000000000...
c    S&S gives EXP_TRANSFORM(16) as 0.5019065783... 
c    S&S gives RAT_TRANSFORM(16) as 0.4988027685...
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) exp ( -x ) * cos ( x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision exponent_min
      parameter ( exponent_min = -80.0D+00 )
      double precision f(n)
      integer i
      double precision x(n)

      do i = 1, n
        if ( -x(i) .lt. exponent_min ) then
          f(i) = 0.0D+00
        else
          f(i) = exp ( -x(i) ) * cos ( x(i) )
        end if
      end do

      return
      end
      subroutine p12_title ( title )

c*********************************************************************72
c
cc P12_TITLE returns the title for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = 'exp ( - x ) * cos ( x )'

      return
      end
      subroutine p13_alpha ( alpha )

c*********************************************************************72
c
cc P13_ALPHA returns ALPHA for problem 13.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p13_exact ( exact )

c*********************************************************************72
c
cc P13_EXACT returns the estimated integral for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = pi / 2.0D+00

      return
      end
      subroutine p13_fun ( n, x, f )

c*********************************************************************72
c
cc P13_FUN evaluates the integrand for problem 13.
c
c  Discussion:
c
c    S&S gives exact value as pi/2 = 1.5707963267948966192...
c    S&S gives Laguerre(16) as       1.4399523793...
c    S&S gives EXP_TRANSFORM(16) as  1.3045186595...
c    S&S gives RAT_TRANSFORM(16) as  0.2046437026...
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) sin ( x ) / x dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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
      double precision x(n)

      do i = 1, n
        if ( x(i) .eq. 0.0D+00 ) then
          f(i) = 1.0D+00
        else
          f(i) = sin ( x(i) ) / x(i)
        end if
      end do

      return
      end
      subroutine p13_title ( title )

c*********************************************************************72
c
cc P13_TITLE returns the title for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = 'sin(x) / x'

      return
      end
      subroutine p14_alpha ( alpha )

c*********************************************************************72
c
cc P14_ALPHA returns ALPHA for problem 14.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p14_exact ( exact )

c*********************************************************************72
c
cc P14_EXACT returns the estimated integral for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.0634618101722400407D+00

      return
      end
      subroutine p14_fun ( n, x, f )

c*********************************************************************72
c
cc P14_FUN evaluates the integrand for problem 14.
c
c  Discussion:
c
c    S&S gives "exact" value as     1.0634618101...
c    Mathematica returns            1.0634618101722400407...
c    S&S gives Laguerre(16) as      1.0634713425...
c    S&S gives EXP_TRANSFORM(16) as 1.0634618101...
c    S&S gives RAT_TRANSFORM(16) as 1.0634574249...
c
c    The FORTRAN version of this routine, compiled with G95, was getting 
c    a floating point exception when evaluating the integrand
c    and using a Laguerre rule of order 64.  So I have had to truncate
c    the evaluation of the exponential.
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) sin ( exp ( - x ) + exp ( - 4 x ) ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2007
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
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision exponent_min
      parameter ( exponent_min = -80.0D+00 )
      double precision f(n)
      integer i
      double precision x(n)

      do i = 1, n
        if ( - x(i) .lt. exponent_min ) then
          f(i) = 0.0D+00
        else if ( -4.0D+00 * x(i) .lt. exponent_min ) then
          f(i) = sin ( exp ( -x(i) ) )
        else
          f(i) = sin ( exp ( -x(i) ) + exp ( -4.0D+00 * x(i) ) )
        end if

      end do

      return
      end
      subroutine p14_title ( title )

c*********************************************************************72
c
cc P14_TITLE returns the title for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2007
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

      title = 'sin ( exp(-x) + exp(-4x) )'

      return
      end
      subroutine p15_alpha ( alpha )

c*********************************************************************72
c
cc P15_ALPHA returns ALPHA for problem 15.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p15_exact ( exact )

c*********************************************************************72
c
cc P15_EXACT returns the estimated integral for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )

      exact = - pi * log ( 10.0D+00 ) / 20.0D+00

      return
      end
      subroutine p15_fun ( n, x, f )

c*********************************************************************72
c
cc P15_FUN evaluates the integrand for problem 15.
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) log(x) / (1+100*x*x) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise deDoncker-Kapenga, 
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983,
c    ISBN: 3540125531,
c    LC: QA299.3.Q36.
c
c  Parameters:
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
      double precision x(n)

      do i = 1, n
        if ( x(i) .eq. 0.0D+00 ) then
          f(i) = - huge ( x(i) )
        else
          f(i) = log ( x(i) ) / ( 1.0D+00 + 100.0D+00 * x(i) * x(i) )
        end if

      end do

      return
      end
      subroutine p15_title ( title )

c*********************************************************************72
c
cc P15_TITLE returns the title for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
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

      title = 'log(x) / ( 1 + 100 x^2 )'

      return
      end
      subroutine p16_alpha ( alpha )

c*********************************************************************72
c
cc P16_ALPHA returns ALPHA for problem 16.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p16_exact ( exact )

c*********************************************************************72
c
cc P16_EXACT returns the estimated integral for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.0D+00

      return
      end
      subroutine p16_fun ( n, x, f )

c*********************************************************************72
c
cc P16_FUN evaluates the integrand for problem 16.
c
c  Integral
c
c    Integral ( 0 .le. x < +oo ) cos ( pi * x / 2 ) / sqrt ( x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise deDoncker-Kapenga, 
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983,
c    ISBN: 3540125531,
c    LC: QA299.3.Q36.
c
c  Parameters:
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
      double precision pi
      parameter ( pi = 3.1415926535897932385D+00 )
      double precision r8_huge
      double precision x(n)

      do i = 1, n
        if ( x(i) .eq. 0.0D+00 ) then
          f(i) = r8_huge ( x(i) )
        else
          f(i) = cos ( pi * x(i) / 2.0D+00 ) / sqrt ( x(i) )
        end if

      end do

      return
      end
      subroutine p16_title ( title )

c*********************************************************************72
c
cc P16_TITLE returns the title for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2007
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

      title = 'cos(pi x / 2 ) / sqrt(x)'

      return
      end
      subroutine p17_alpha ( alpha )

c*********************************************************************72
c
cc P17_ALPHA returns ALPHA for problem 17.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p17_exact ( exact )

c*********************************************************************72
c
cc P17_EXACT returns the exact integral for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      double precision beta
      parameter ( beta = 2.0D+00 )
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = sqrt ( pi ) * cos ( 0.5D+00 * atan ( 2.0D+00**beta ) ) 
     &  / sqrt ( sqrt ( 1.0D+00 + 0.25D+00**beta) )

      return
      end
      subroutine p17_fun ( n, x, fx )

c*********************************************************************72
c
cc P17_FUN evaluates the integrand for problem 17.
c
c  Integral:
c
c    Integral ( 0 .le. x < +oo ) exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision beta
      parameter ( beta = 2.0D+00 )
      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = exp ( - x(i) / 2.0D+00**beta ) * cos ( x(i) )
     &      / sqrt ( x(i) )
        end if

      end do

      return
      end
      subroutine p17_title ( title )

c*********************************************************************72
c
cc P17_TITLE returns the title for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      title = 'exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x )'

      return
      end
      subroutine p18_alpha ( alpha )

c*********************************************************************72
c
cc P18_ALPHA returns ALPHA for problem 18.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p18_exact ( exact )

c*********************************************************************72
c
cc P18_EXACT returns the exact integral for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      double precision beta
      parameter ( beta = 1.0D+00 )
      double precision exact

      exact = 2.0D+00**( 3.0D+00 * beta + 1.0D+00 )

      return
      end
      subroutine p18_fun ( n, x, fx )

c*********************************************************************72
c
cc P18_FUN evaluates the integrand for problem 18.
c
c  Integral:
c
c    Integral ( 0 .le. x < +oo ) x^2 * exp ( - x / 2^beta ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision beta
      parameter ( beta = 1.0D+00 )
      double precision fx(n)
      double precision x(n)

      fx(1:n) = x(1:n)**2 * exp ( - x(1:n) / 2**beta )

      return
      end
      subroutine p18_title ( title )

c*********************************************************************72
c
cc P18_TITLE returns the title for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      title = 'x^2 * exp ( - x / 2^beta )'

      return
      end
      subroutine p19_alpha ( alpha )

c*********************************************************************72
c
cc P19_ALPHA returns ALPHA for problem 19.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p19_exact ( exact )

c*********************************************************************72
c
cc P19_EXACT returns the exact integral for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      double precision beta
      parameter ( beta = 0.5D+00 )
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      if ( beta .eq. 1.0D+00 ) then
        exact = 1.0D+00 / 10.0D+00
      else
        exact = ( 1.0D+00 - beta ) * pi 
     &    / ( 10.0D+00**beta * sin ( pi * beta ) )
      end if

      return
      end
      subroutine p19_fun ( n, x, fx )

c*********************************************************************72
c
cc P19_FUN evaluates the integrand for problem 19.
c
c  Integral:
c
c    Integral ( 0 .le. x < +oo ) x^(beta-1) / ( 1 + 10 x )^2 dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision beta
      parameter ( beta = 0.5D+00 )
      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( beta .eq. 1.0D+00 ) then
          fx(i) = 1.0D+00 / ( 1.0D+00 + 10.0D+00 * x(i) )**2
        else if ( beta .lt. 1.0D+00 .and. x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = x(i)**( beta - 1.0D+00 ) 
     &      / ( 1.0D+00 + 10.0D+00 * x(i) )**2
        end if

      end do

      return
      end
      subroutine p19_title ( title )

c*********************************************************************72
c
cc P19_TITLE returns the title for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      title = 'x^(beta-1) / ( 1 + 10 x )^2'

      return
      end
      subroutine p20_alpha ( alpha )

c*********************************************************************72
c
cc P20_ALPHA returns ALPHA for problem 10.
c
c  Discussion:
c
c    ALPHA is the lower, finite limit of integration in the integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.0D+00

      return
      end
      subroutine p20_exact ( exact )

c*********************************************************************72
c
cc P20_EXACT returns the exact integral for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      double precision beta
      parameter ( beta = 1.0D+00 )
      double precision exact

      exact = 
     &  ( 
     &    log ( 1.5D+00 ) / 2.0D+00**beta 
     &    - 1.0D+00 / 2.0D+00**( beta + 1.0D+00 ) * 
     &    log ( ( 16.0D+00 + 0.25D+00**beta ) 
     &    / ( 1.0D+00 + 0.25D+00**beta ) ) 
     &    - atan ( 2.0D+00**( beta + 2.0D+00 ) ) 
     &    - atan ( 2.0D+00**beta ) 
     &  ) / ( 1.0D+00 + 0.25D+00**beta )

      return
      end
      subroutine p20_fun ( n, x, fx )

c*********************************************************************72
c
cc P20_FUN evaluates the integrand for problem 20.
c
c  Integral:
c
c    Integral ( 0 .le. x < +oo ) 
c      1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision beta
      parameter ( beta = 1.0D+00 )
      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( ( x(i) - 1.0D+00 )**2 + 0.25D+00**beta .eq. 0.0D+00 .or. 
     &       x(i) .eq. 2.0D+00 ) then

          fx(i) = 0.0D+00

        else

          fx(i) = 1.0D+00 / 
     &      ( 2.0D+00**beta 
     &      * ( ( x(i) - 1.0D+00 )**2 + 0.25D+00**beta ) 
     &      * ( x(i) - 2.0D+00 ) )

        end if

      end do

      return
      end
      subroutine p20_title ( title )

c*********************************************************************72
c
cc P20_TITLE returns the title for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
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

      title = 
     &  '1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) )'

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
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
c    01 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates Gamma(X) for a real argument.
c
c  Discussion:
c
c    This routine calculates the gamma function for a real argument X.
c    Computation is based on an algorithm outlined in reference 1.
c    The program uses rational functions that approximate the gamma
c    function to at least 20 significant decimal digits.  Coefficients
c    for the approximation over the interval (1,2) are unpublished.
c    Those for the approximation for 12 .le. X are from reference 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    An Overview of Software Development for Special Functions,
c    in Numerical Analysis Dundee, 1975,
c    edited by GA Watson,
c    Lecture Notes in Mathematics 506,
c    Springer, 1976.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      integer i
      integer n
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
c
c  Mathematical constants
c
      data sqrtpi /0.9189385332046727417803297D+00/
      data pi /3.1415926535897932384626434D+00/
c
c  Machine dependent parameters
c
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (1,2).
c
      data p/
     & -1.71618513886549492533811d+00,
     &  2.47656508055759199108314d+01,
     & -3.79804256470945635097577d+02,
     &  6.29331155312818442661052d+02,
     &  8.66966202790413211295064d+02,
     & -3.14512729688483675254357d+04,
     & -3.61444134186911729807069d+04,
     &  6.64561438202405440627855d+04/

      data q/
     & -3.08402300119738975254353d+01,
     &  3.15350626979604161529144d+02,
     & -1.01515636749021914166146d+03,
     & -3.10777167157231109440444d+03,
     &  2.25381184209801510330112d+04,
     &  4.75584627752788110767815d+03,
     & -1.34659959864969306392456d+05,
     & -1.15132259675553483497211d+05/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     & -1.910444077728D-03,
     &  8.4171387781295D-04,
     & -5.952379913043012D-04,
     &  7.93650793500350248D-04,
     & -2.777777777777681622553D-03,
     &  8.333333333333333331554247D-02,
     &  5.7083835261D-03/

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
c
c  Argument is negative.
c
      if ( y .le. 0.0D+00 ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. 0.0D+00 ) then

          if ( y1 .ne. aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + 1.0D+00

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Argument is positive.
c
      if ( y .lt. eps ) then
c
c  Argument < EPS.
c
        if ( xminin .le. y ) then
          res = 1.0D+00 / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. 12.0D+00 ) then

        y1 = y
c
c  0.0 < argument < 1.0.
c
        if ( y .lt. 1.0D+00 ) then

          z = y
          y = y + 1.0D+00
c
c  1.0 < argument < 12.0.
c  Reduce argument if necessary.
c
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - 1.0D+00

        end if
c
c  Evaluate approximation for 1.0 < argument < 2.0.
c
        xnum = 0.0D+00
        xden = 1.0D+00
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + 1.0D+00
c
c  Adjust result for case  0.0 < argument < 1.0.
c
        if ( y1 .lt. y ) then

          res = res / y1
c
c  Adjust result for case 2.0 < argument < 12.0.
c
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + 1.0D+00
          end do

        end if

      else
c
c  Evaluate for 12.0 .le. argument.
c
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - 0.5D+00 ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Final adjustments and return.
c
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. 1.0D+00 ) then
        res = fact / res
      end if

      r8_gamma = res

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
