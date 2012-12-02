      subroutine vandermonde_approx_1d_coef ( n, m, x, y, c )

c*********************************************************************72
c
cc VANDERMONDE_APPROX_1D_COEF computes a 1D polynomial approximant.
c
c  Discussion:
c
c    We assume the approximating function has the form
c
c      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
c
c    We have n data values (x(i),y(i)) which must be approximated:
c
c      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
c
c    This can be cast as an Nx(M+1) linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
c      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
c      [ .................. ] [ ... ] = [ ... ]
c      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
c
c    In the typical case, N is greater than M+1 (we have more data and equations
c    than degrees of freedom) and so a least squares solution is appropriate,
c    in which case the computed polynomial will be a least squares approximant
c    to the data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, integer M, the degree of the polynomial.
c
c    Input, double precision X(N), Y(N), the data values.
c
c    Output, double precision C(0:M), the coefficients of the approximating
c    polynomial.  C(0) is the constant term, and C(M) multiplies X^M.
c
      implicit none

      integer m
      integer n

      double precision a(1:n,0:m)
      double precision c(0:m)
      double precision x(n)
      double precision y(n)

      call vandermonde_approx_1d_matrix ( n, m, x, a )

      call qrsolve ( n, m + 1, a, y, c )

      return
      end
      subroutine vandermonde_approx_1d_matrix ( n, m, x, a )

c*********************************************************************72
c
cc VANDERMONDE_APPROX_1D_MATRIX computes a Vandermonde 1D approximation matrix.
c
c  Discussion:
c
c    We assume the approximant has the form
c
c      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
c
c    We have n data values (x(i),y(i)) which must be approximated:
c
c      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
c
c    This can be cast as an Nx(M+1) linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
c      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
c      [ .................. ] [ ... ] = [ ... ]
c      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, integer M, the degree of the polynomial.
c
c    Input, double precision X(N), the data values.
c
c    Output, double precision A(N,0:M), the Vandermonde matrix for X.
c
      implicit none

      integer m
      integer n

      double precision a(n,0:m)
      integer i
      integer j
      double precision x(n)

      do i = 1, n
        a(1:n,0) = 1.0D+00
      end do

      do j = 1, m
        do i = 1, n
          a(i,j) = a(i,j-1) * x(i)
        end do
      end do

      return
      end
