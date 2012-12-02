      subroutine vandermonde_interp_1d_coef ( n, x, y, c )

c*********************************************************************72
c
cc VANDERMONDE_INTERP_1D_COEF computes a 1D polynomial interpolant.
c
c  Discussion:
c
c    We assume the interpolant has the form
c
c      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
c
c    We have n data values (x(i),y(i)) which must be interpolated:
c
c      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
c
c    This can be cast as an NxN linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
c      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
c      [ ...................... ] [ ... ] = [ ... ]
c      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
c
c    and if the x values are distinct, the system is theoretically
c    invertible, so we can retrieve the coefficient vector c and
c    evaluate the interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, double precision X(N), Y(N), the data values.
c
c    Output, double precision C(N), the coefficients of the interpolating
c    polynomial.  C(1) is the constant term, and C(N) multiplies X^(N-1).
c
      implicit none

      integer n

      double precision a(n,n)
      double precision c(n)
      double precision x(n)
      double precision y(n)

      call vandermonde_interp_1d_matrix ( n, x, a )

      call qr_solve ( n, n, a, y, c )

      return
      end
      subroutine vandermonde_interp_1d_matrix ( n, x, a )

c*********************************************************************72
c
cc VANDERMONDE_INTERP_1D_MATRIX computes a Vandermonde 1D interpolation matrix.
c
c  Discussion:
c
c    We assume the interpolant has the form
c
c      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
c
c    We have n data values (x(i),y(i)) which must be interpolated:
c
c      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
c
c    This can be cast as an NxN linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
c      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
c      [ ...................... ] [ ... ] = [ ... ]
c      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
c
c    and if the x values are distinct, the matrix A is theoretically
c    invertible (though in fact, generally badly conditioned).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, real X(N), the data values.
c
c    Output, real A(N,N), the Vandermonde matrix for X.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
      double precision x(n)

      do i = 1, n
        a(i,1) = 1.0D+00
      end do

      do j = 2, n
        do i = 1, n
          a(i,j) = x(i) * a(i,j-1)
        end do
      end do

      return
      end
