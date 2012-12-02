      subroutine hermite_cubic_integral ( x1, f1, d1, x2, f2, d2, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_INTEGRAL returns the integral of a Hermite cubic polynomial.
c
c  Discussion:
c
c    The integral is taken over the definition interval [X1,X2].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, F1, D1, the left endpoint, function value
c    and derivative.
c
c    Input, double precision X2, F2, D2, the right endpoint, function value
c    and derivative.
c
c    Output, double precision Q, the integral of the Hermite cubic polynomial
c    over the interval X1 <= X <= X2.
c
      implicit none

      double precision d1
      double precision d2
      double precision f1
      double precision f2
      double precision h
      double precision q
      double precision x1
      double precision x2

      h = x2 - x1

      q = 0.5D+00 * h * ( f1 + f2 + h * ( d1 - d2 ) / 6.0D+00 )

      return
      end
      subroutine hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, 
     &  b, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic polynomial from A to B.
c
c  Discussion:
c
c    A and B may be scalars, or one may be a vector, or both
c    may be vectors of the same size.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, F1, D1, the left endpoint, function value
c    and derivative.
c
c    Input, double precision X2, F2, D2, the right endpoint, function value
c    and derivative.
c
c    Input, double precision A, B, the left and right endpoints of the interval
c    of integration.
c
c    Output, double precision Q, the integral of the Hermite cubic polynomial
c    over the interval A <= X <= B.
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision dterm
      double precision f1
      double precision f2
      double precision fterm
      double precision h
      double precision phia1
      double precision phia2
      double precision phib1
      double precision phib2
      double precision psia1
      double precision psia2
      double precision psib1
      double precision psib2
      double precision q
      double precision ta1
      double precision ta2
      double precision tb1
      double precision tb2
      double precision ua1
      double precision ua2
      double precision ub1
      double precision ub2
      double precision x1
      double precision x2

      h = x2 - x1

      ta1 = ( a - x1 ) / h
      ta2 = ( x2 - a ) / h
      tb1 = ( b - x1 ) / h
      tb2 = ( x2 - b ) / h

      ua1 = ta1 * ta1 * ta1
      phia1 = ua1 * ( 2.0D+00 - ta1 )
      psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )

      ua2 = ta2 * ta2 * ta2
      phia2 =  ua2 * ( 2.0D+00 - ta2 )
      psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

      ub1 = tb1 * tb1 * tb1
      phib1 = ub1 * ( 2.0D+00 - tb1 )
      psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )

      ub2 = tb2 * tb2 * tb2
      phib2 =  ub2 * ( 2.0D+00 - tb2 )
      psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

      fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 )
      dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) 
     &  * ( h / 6.0D+00 )

      q = 0.5D+00 * h * ( fterm + dterm )

      return
      end
      subroutine hermite_cubic_lagrange_integral ( x1, x2, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_LAGRANGE_INTEGRAL: Hermite cubic Lagrange integrals.
c
c  Discussion:
c
c    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
c    (F1,D1,F2,D2) satisfies:
c
c      P(X1) = F1,
c      P'(X1) = D1,
c      P(X2) = F2,
c      P'(X2) = D2.
c
c    We can determine four Lagrange polynomials L1(X) through L4(X) so that
c
c      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
c
c    This function returns the integrals of these four polynomials over
c    the domain of definition [X1,X2].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, X2, the endpoints.
c
c    Output, double precision Q(4), the integrals of the Hermite cubic
c    Lagrange polynomials from X1 to X2.
c
      implicit none

      double precision h
      double precision q(4)
      double precision x1
      double precision x2

      h = x2 - x1

      q(1) =   h     / 2.0D+00
      q(2) =   h * h / 12.0D+00
      q(3) =   h     / 2.0D+00
      q(4) = - h * h / 12.0D+00

      return
      end
      subroutine hermite_cubic_lagrange_integrate ( x1, x2, a, b, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_LAGRANGE_INTEGRATE: integrate Hermite cubic Lagrange polys.
c
c  Discussion:
c
c    A and B may be scalars, or one may be a vector, or both
c    may be vectors of the same size.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, real X1, X2, the endpoints of the interval of definition.
c
c    Input, real A, B, the left and right endpoints of the interval
c    of integration.
c
c    Output, real Q(4), the integrals of the Hermite cubic Lagrange polynomials
c    over the interval A <= X <= B.
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision f1
      double precision f2
      double precision h
      double precision phia1
      double precision phia2
      double precision phib1
      double precision phib2
      double precision psia1
      double precision psia2
      double precision psib1
      double precision psib2
      double precision q(4)
      double precision ta1
      double precision ta2
      double precision tb1
      double precision tb2
      double precision ua1
      double precision ua2
      double precision ub1
      double precision ub2
      double precision x1
      double precision x2

      h = x2 - x1
      ta1 = ( a - x1 ) / h
      ta2 = ( x2 - a ) / h
      tb1 = ( b - x1 ) / h
      tb2 = ( x2 - b ) / h

      ua1 = ta1 * ta1 * ta1
      phia1 = ua1 * ( 2.0D+00 - ta1 )
      psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )

      ua2 = ta2 * ta2 * ta2
      phia2 =  ua2 * ( 2.0D+00 - ta2 )
      psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

      ub1 = tb1 * tb1 * tb1
      phib1 = ub1 * ( 2.0D+00 - tb1 )
      psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )

      ub2 = tb2 * tb2 * tb2
      phib2 =  ub2 * ( 2.0D+00 - tb2 )
      psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

      q(1) = 0.5D+00 * h * ( phia2 - phib2 )
      q(2) = 0.5D+00 * h * ( psia2 - psib2 ) * ( h / 6.0D+00 )
      q(3) = 0.5D+00 * h * ( phib1 - phia1 )
      q(4) = 0.5D+00 * h * ( psib1 - psia1 ) * ( h / 6.0D+00 )

      return
      end
      subroutine hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, 
     &  s, t )

c*********************************************************************72
c
cc HERMITE_CUBIC_LAGRANGE_VALUE: evaluate Hermite cubic Lagrange polynomials.
c
c  Discussion:
c
c    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
c    (F1,D1,F2,D2) satisfies:
c
c      P(X1) = F1,
c      P'(X1) = D1,
c      P(X2) = F2,
c      P'(X2) = D2.
c
c    We can determine four Lagrange polynomials L1(X) through L4(X) so that
c
c      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
c
c    This function returns the values and derivatives of these four
c    polynomials.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, X2, the endpoints.
c
c    Input, integer N, the number of sample points.
c
c    Input, double precision X(N), the sample points.
c
c    Output, double precision F(4,N), D(4,N), S(4,N), T(4,N), the value
c    and first three derivatives of the Hermite cubic Lagrange polynomials at X.
c
      implicit none

      integer n

      double precision d(4,n)
      double precision dx(n)
      double precision f(4,n)
      double precision h
      integer j
      double precision s(4,n)
      double precision t(4,n)
      double precision x(n)
      double precision x1
      double precision x2

      h = x2 - x1
      do j = 1, n
        dx(j) = x(j) - x1
      end do
c
c  F1.
c
      do j = 1, n

        f(1,j) = 1.0D+00 + ( dx(j)**2   / h**2 )
     &    * ( - 3.0D+00 + ( dx(j) / h ) *  2.0D+00 )
        d(1,j) = ( dx(j)   / h**2 ) 
     &    * ( - 6.0D+00 + ( dx(j) / h ) *  6.0D+00 )
        s(1,j) = ( 1.0D+00 / h**2 ) 
     &    * ( - 6.0D+00 + ( dx(j) / h ) * 12.0D+00 )
        t(1,j) = ( 1.0D+00 / h**3 ) * 12.0D+00
c
c  D1
c
        f(2,j) = dx(j)   + ( dx(j)**2   / h    ) 
     &    * ( - 2.0D+00 + ( dx(j) / h )           )
        d(2,j) = 1.0D+00 + ( dx(j)      / h    ) 
     &    * ( - 4.0D+00 + ( dx(j) / h ) * 3.0D+00 )
        s(2,j) =           ( 1.0D+00 / h    ) 
     &    * ( - 4.0D+00 + ( dx(j) / h ) * 6.0D+00 )
        t(2,j) =           ( 1.0D+00 / h**2 ) * 6.0D+00
c
c  F2
c
        f(3,j) = ( dx(j)**2   / h**2 ) 
     &    * ( 3.0D+00 -  2.0D+00 * ( dx(j) / h ) )
        d(3,j) = ( dx(j)      / h**2 ) 
     &    * ( 6.0D+00 -  6.0D+00 * ( dx(j) / h ) )
        s(3,j) = ( 1.0D+00 / h**2 ) 
     &    * ( 6.0D+00 - 12.0D+00 * ( dx(j) / h ) )
        t(3,j) = ( 1.0D+00 / h**3 ) 
     &    * (         - 12.0D+00              )
c
c  D2
c
        f(4,j) = ( dx(j)**2   / h ) 
     &    * ( - 1.0D+00 + ( dx(j) / h )           )
        d(4,j) = ( dx(j)      / h ) 
     &    * ( - 2.0D+00 + ( dx(j) / h ) * 3.0D+00 )
        s(4,j) = ( 1.0D+00 / h ) 
     &    * ( - 2.0D+00 + ( dx(j) / h ) * 6.0D+00 )
        t(4,j) = ( 1.0D+00 / h )                            * 6.0D+00

      end do

      return
      end
      subroutine hermite_cubic_spline_integral ( nn, xn, fn, dn, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_SPLINE_INTEGRAL: Hermite cubic spline integral.
c
c  Discussion:
c
c    The integral is taken over the definition interval [X(1),X(NN)].
c
c    Note that if the intervals are equal in size, then the derivative
c    information in DN has no effect on the integral value,
c    except for the first and last entries.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, integer NN, the number of data points.
c
c    Input, double precision XN(NN), the coordinates of the data points.
c    The entries in XN must be in strictly ascending order.
c
c    Input, double precision FN(NN), the function values.
c
c    Input, double precision DN(NN), the derivative values.
c
c    Output, double precision Q, the integral of the Hermite cubic spline
c    over the interval X(1) <= X <= X(NN).
c
      implicit none

      integer nn

      double precision dn(nn)
      double precision fn(nn)
      integer i
      double precision q
      double precision xn(nn)

      q = 0.0D+00
      do i = 1, nn - 1
        q = q + 0.5D+00 * ( xn(i+1) - xn(i) ) * ( fn(i) + fn(i+1) 
     &    + ( xn(i+1) - xn(i) ) * ( dn(i) - dn(i+1) ) / 6.0D+00 )
      end do

      return
      end
      subroutine hermite_cubic_spline_quad_rule ( nn, xn, w )

c*********************************************************************72
c
cc HERMITE_CUBIC_SPLINE_QUAD_RULE: Hermite cubic spline quadrature rule.
c
c  Discussion:
c
c    The integral is taken over the definition interval [X(1),X(NN)].
c
c    Note that if the intervals are equal in size, then the derivative
c    information in DN has no effect on the integral value,
c    except for the first and last entries.
c
c    The quadrature rule is
c
c      Integral ( XN(1) <= X <= XN(NN) ) F(X) dX is approximately
c
c      Sum ( 1 <= I <= NN ) W(1,I) * F(X(I)) + W(2,I) * F'(X(I))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, integer NN, the number of data points.
c
c    Input, double precision XN(NN), the coordinates of the data points.
c    The entries in XN must be in strictly ascending order.
c
c    Output, double precision W(2,NN), the quadrature weights for F(1:NN)
c    and DN(1:NN).
c
      implicit none

      integer nn

      integer j
      double precision w(2,nn)
      double precision xn(nn)

      w(1,1)      = 0.5D+00 * ( xn(2)    - xn(1)      )
      do j = 2, nn - 1
        w(1,j) = 0.5D+00 * ( xn(j+1) - xn(j-1) )
      end do
      w(1,nn)     = 0.5D+00 * ( xn(nn)   - xn(nn-1)   )

      w(2,1)      =   ( xn(2) - xn(1) )**2 / 12.0D+00
      do j = 2, nn - 1
        w(2,j) =   ( xn(j+1) - xn(j-1) ) 
     &    * ( xn(j+1) - 2.0D+00 * xn(j) + xn(j-1) ) / 12.0D+00
      end do
      w(2,nn)     = - ( xn(nn-1) - xn(nn) )**2 / 12.0D+00

      return
      end
      subroutine hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, 
     &  b, q )

c*********************************************************************72
c
cc HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite cubic spline over [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, integer NN, the number of data points.
c
c    Input, double precision XN(NN), the coordinates of the data points.
c    The entries in XN must be in strictly ascending order.
c
c    Input, double precision FN(NN), the function values.
c
c    Input, double precision DN(NN), the derivative values.
c
c    Input, integer N, the number of integration intervals.
c
c    Input, double precision A(N), B(N), the integration endpoints.
c
c    Output, double precision Q(N), the integral over the interval [A,B].
c
      implicit none

      integer n
      integer nn

      double precision a(n)
      double precision aa
      double precision b(n)
      double precision bb
      double precision dn(nn)
      double precision fn(nn)
      integer i
      integer ii
      integer j
      integer k
      double precision q(n)
      double precision qq
      double precision s
      double precision xn(nn)

      do ii = 1, n
        q(ii) = 0.0D+00
      end do

      do ii = 1, n

        if ( a(ii) .le. b(ii) ) then
          aa = a(ii)
          bb = b(ii)
          s = + 1.0D+00
        else
          aa = b(ii)
          bb = a(ii)
          s = - 1.0D+00
        end if

        call r8vec_bracket3 ( nn, xn, aa, i )
        call r8vec_bracket3 ( nn, xn, bb, j )
c
c  Evaluate the polynomial with the appropriate data.
c
        if ( i .eq. j ) then

          call hermite_cubic_integrate ( xn(i), fn(i), dn(i), 
     &      xn(i+1), fn(i+1), dn(i+1), aa, bb, q(ii) )

        else

          call hermite_cubic_integrate ( xn(i), fn(i), dn(i), 
     &      xn(i+1), fn(i+1), dn(i+1), aa, xn(i+1), qq )

          q(ii) = qq

          do k = i + 1, j - 1

            call hermite_cubic_integral ( xn(k), fn(k), dn(k), 
     &        xn(k+1), fn(k+1), dn(k+1), qq )

            q(ii) = q(ii) + qq

          end do

          call hermite_cubic_integrate ( xn(j), fn(j), dn(j), 
     &      xn(j+1), fn(j+1), dn(j+1), xn(j), bb, qq )

          q(ii) = q(ii) + qq

        end if

        q(ii) = s * q(ii)

      end do

      return
      end
      subroutine hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, 
     &  f, d, s, t )

c*********************************************************************72
c
cc HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, integer NN, the number of data points.
c
c    Input, double precision XN(NN), the coordinates of the data points.
c    The entries in XN must be in strictly ascending order.
c
c    Input, double precision FN(NN), the function values.
c
c    Input, double precision DN(NN), the derivative values.
c
c    Input, integer N, the number of sample points.
c
c    Input, double precision X(N), the coordinates of the sample points.
c
c    Output, double precision F(N), the function value at the sample points.
c
c    Output, double precision D(N), the derivative value at the sample points.
c
c    Output, double precision S(N), the second derivative value at the
c    sample points.
c
c    Output, double precision T(N), the third derivative value at the
c    sample points.
c
      implicit none

      integer n
      integer nn

      double precision d(n)
      double precision dn(nn)
      double precision f(n)
      double precision fn(nn)
      integer i
      integer left
      double precision s(n)
      double precision t(n)
      double precision x(n)
      double precision xn(nn)

      do i = 1, n

        call r8vec_bracket3 ( nn, xn, x(i), left )

        call hermite_cubic_value ( xn(left), fn(left), dn(left), 
     &    xn(left+1), fn(left+1), dn(left+1), 1, x(i), f(i), d(i), 
     &    s(i), t(i) )

      end do

      return
      end
      subroutine hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, 
     &  d2, c0, c1, c2, c3 )

c*********************************************************************72
c
cc HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, F1, D1, the left endpoint, function value
c    and derivative.
c
c    Input, double precision X2, F2, D2, the right endpoint, function value
c    and derivative.
c
c    Output, double precision C0, C1, C2, C3, the power form of the polynomial.
c
      implicit none

      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision d
      double precision d1
      double precision d2
      double precision df
      double precision f1
      double precision f2
      double precision h
      double precision x1
      double precision x2

      h =    x2 - x1
      df = ( f2 - f1 ) / h
c
c  Polynomial in terms of X - X1:
c
      c0 = f1
      c1 = d1
      c2 = - ( 2.0D+00 * d1 - 3.0D+00 * df + d2 ) / h
      c3 =   (           d1 - 2.0D+00 * df + d2 ) / h / h
c
c  Shift polynomial to X.
c
      c2 = c2 - x1 * c3
      c1 = c1 - x1 * c2
      c0 = c0 - x1 * c1
      c2 = c2 - x1 * c3
      c1 = c1 - x1 * c2
      c2 = c2 - x1 * c3

      return
      end
      subroutine hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, 
     &  f, d, s, t )

c*********************************************************************72
c
cc HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
c
c  Discussion:
c
c    The input arguments can be vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision X1, F1, D1, the left endpoint, function value
c    and derivative.
c
c    Input, double precision X2, F2, D2, the right endpoint, function value
c    and derivative.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the points at which the Hermite cubic
c    is to be evaluated.
c
c    Output, double precision F(N), D(N), S(N), T(N), the value and first
c    three derivatives of the Hermite cubic at X.
c
      implicit none

      integer n

      double precision c2
      double precision c3
      double precision d(n)
      double precision d1
      double precision d2
      double precision df
      double precision dx(n)
      double precision f(n)
      double precision f1
      double precision f2
      double precision h
      integer i
      double precision s(n)
      double precision t(n)
      double precision x(n)
      double precision x1
      double precision x2

      h =    x2 - x1
      df = ( f2 - f1 ) / h

      c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h
      c3 =   (       d1 - 2.0 * df + d2 ) / h / h

      do i = 1, n
        dx(i) = x(i) - x1
        f(i) = f1 + dx(i) 
     &           * ( d1 + dx(i) * ( c2 + dx(i) * c3 ) )
        d(i) = d1 + dx(i) * ( 2.0D+00 * c2 + dx(i) * 3.0D+00 * c3 )
        s(i) = 2.0D+00 * c2 + dx(i) * 6.0D+00 * c3
        t(i) = 6.0D+00 * c3
      end do

      return
      end
      subroutine power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, 
     &  f1, d1, f2, d2 )

c*********************************************************************72
c
cc POWER_CUBIC_TO_HERMITE_CUBIC converts a power cubic to Hermite form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision C0, C1, C2, C3, the power form of the
c    polynomial.
c
c    Input, double precision X1, X2, the left and right endpoints of
c    the Hermite form.
c
c    Output, double precision F1, D1, the function and derivative values at X1.
c
c    Output, double precision F2, D2, the function and derivative values at X2.
c
      implicit none

      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision d1
      double precision d2
      double precision f1
      double precision f2
      double precision x1
      double precision x2

      f1 = c0 + x1 * ( c1 + x1 * (           c2 + x1           * c3 ) )
      d1 =             c1 + x1 * ( 2.0D+00 * c2 + x1 * 3.0D+00 * c3 )

      f2 = c0 + x2 * ( c1 + x2 * (           c2 + x2           * c3 ) )
      d2 =             c1 + x2 * ( 2.0D+00 * c2 + x2 * 3.0D+00 * c3 )

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
      subroutine r8vec_bracket3 ( n, t, tval, left )

c*********************************************************************72
c
cc R8VEC_BRACKET3 finds the interval containing or nearest a given value.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The routine always returns the index LEFT of the sorted array
c    T with the property that either
c    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
c    *  T .lt. T(LEFT) = T(1), or
c    *  T > T(LEFT+1) = T(N).
c
c    The routine is useful for interpolation problems, where
c    the abscissa must be located within an interval of data
c    abscissas for interpolation, or the "nearest" interval
c    to the (extreme) abscissa must be found so that extrapolation
c    can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of the input array.
c
c    Input, double precision T(N), an array that has been sorted
c    into ascending order.
c
c    Input, double precision TVAL, a value to be bracketed by entries of T.
c
c    Input/output, integer LEFT.
c    On input, if 1 .le. LEFT .le. N-1, LEFT is taken as a suggestion for the
c    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
c    is searched first, followed by the appropriate interval to the left
c    or right.  After that, a binary search is used.
c    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
c    is the closest to TVAL; it either contains TVAL, or else TVAL
c    lies outside the interval [ T(1), T(N) ].
c
      implicit none

      integer n

      integer high
      integer left
      integer low
      integer mid
      double precision t(n)
      double precision tval
c
c  Check the input data.
c
      if ( n .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
        write ( *, '(a)' ) '  N must be at least 2.'
        stop
      end if
c
c  If LEFT is not between 1 and N-1, set it to the middle value.
c
      if ( left .lt. 1 .or. n - 1 .lt. left ) then
        left = ( n + 1 ) / 2
      end if
c
c  CASE 1: TVAL .lt. T(LEFT):
c  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
c
      if ( tval .lt. t(left) ) then

        if ( left .eq. 1 ) then
          return
        else if ( left .eq. 2 ) then
          left = 1
          return
        else if ( t(left-1) .le. tval ) then
          left = left - 1
          return
        else if ( tval .le. t(2) ) then
          left = 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
c
        low = 2
        high = left - 2

10      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 10
c
c  CASE2: T(LEFT+1) .lt. TVAL:
c  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
c
      else if ( t(left+1) .lt. tval ) then

        if ( left .eq. n - 1 ) then
          return
        else if ( left .eq. n - 2 ) then
          left = left + 1
          return
        else if ( tval .le. t(left+2) ) then
          left = left + 1
          return
        else if ( t(n-1) .le. tval ) then
          left = n - 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
c
        low = left + 2
        high = n - 2

20      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 20
c
c  CASE3: T(LEFT) .le. TVAL .le. T(LEFT+1):
c  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
c
      else

      end if

      return
      end
      subroutine r8vec_even ( n, alo, ahi, a )

c*********************************************************************72
c
cc R8VEC_EVEN returns an R8VEC of evenly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    If N is 1, then the midpoint is returned.
c
c    Otherwise, the two endpoints are returned, and N-2 evenly
c    spaced points between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision ALO, AHI, the low and high values.
c
c    Output, double precision A(N), N evenly spaced values.
c    Normally, A(1) = ALO and A(N) = AHI.
c    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
c
      implicit none

      integer n

      double precision a(n)
      double precision ahi
      double precision alo
      integer i

      if ( n .eq. 1 ) then

        a(1) = 0.5D+00 * ( alo + ahi )

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * alo
     &           + dble (     i - 1 ) * ahi )
     &           / dble ( n     - 1 )
        end do

      end if

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
