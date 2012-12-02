      function triangle_num ( n )

c*********************************************************************72
c
cc TRIANGLE_NUM returns the N-th triangular number.
c
c  Definition:
c
c    The N-th triangular number T(N) is formed by the sum of the first
c    N integers:
c
c      T(N) = sum ( 1 <= I <= N ) I
c
c    By convention, T(0) = 0.
c
c  Formula:
c
c    T(N) = ( N * ( N + 1 ) ) / 2
c
c  First Values:
c
c     0
c     1
c     3
c     6
c    10
c    15
c    21
c    28
c    36
c    45
c    55
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the desired number, 
c    which must be at least 0.
c
c    Output, integer TRIANGLE_NUM, the N-th triangular number.
c
      implicit none

      integer n
      integer triangle_num

      triangle_num = ( n * ( n + 1 ) ) / 2

      return
      end
      subroutine vandermonde_approx_2d_coef ( n, m, x, y, z, c )

c*********************************************************************72
c
cc VANDERMONDE_APPROX_2D_COEF computes a 2D polynomial approximant.
c
c  Discussion:
c
c    We assume the approximating function has the form of a polynomial
c    in X and Y of total degree M.
c
c      p(x,y) = c00 
c             + c10 * x                + c01 *  y
c             + c20 * x^2   + c11 * xy + c02 * y^2
c             + ...
c             + cm0 * x^(m) + ...      + c0m * y^m.
c
c    If we let T(K) = the K-th triangular number 
c            = sum ( 1 <= I <= K ) I
c    then the number of coefficients in the above polynomial is T(M+1).
c
c    We have n data locations (x(i),y(i)) and values z(i) to approximate:
c
c      p(x(i),y(i)) = z(i)
c
c    This can be cast as an NxT(M+1) linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
c      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
c      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
c      [ ...................... ] [ ... ] = [ ... ]
c      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
c
c    In the typical case, N is greater than T(M+1) (we have more data and 
c    equations than degrees of freedom) and so a least squares solution is 
c    appropriate, in which case the computed polynomial will be a least squares
c    approximant to the data.
c
c    The polynomial defined by the T(M+1) coefficients C could be evaluated 
c    at the Nx2-vector x by the command
c
c      pval = r8poly_value_2d ( m, c, n, x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, integer M, the maximum degree of the polynomial.
c
c    Input, double precision X(N), Y(N) the data locations.
c
c    Input, double precision Z(N), the data values.
c
c    Output, double precision C(T(M+1)), the coefficients of the approximating
c    polynomial.  C(1) is the constant term, and C(T(M+1)) multiplies Y^M.
c
      implicit none

      integer m
      integer n

      double precision a(n,(m+1)*(m+2)/2)
      double precision c(*)
      integer tm
      integer triangle_num
      double precision x(n)
      double precision y(n)
      double precision z(n)

      tm = triangle_num ( m + 1 )

      call vandermonde_approx_2d_matrix ( n, m, tm, x, y, a )

      call qr_solve ( n, tm, a, z, c )

      return
      end
      subroutine vandermonde_approx_2d_matrix ( n, m, tm, x, y, a )

c*********************************************************************72
c
cc VANDERMONDE_APPROX_2D_MATRIX computes a Vandermonde 2D approximation matrix.
c
c  Discussion:
c
c    We assume the approximating function has the form of a polynomial
c    in X and Y of total degree M.
c
c      p(x,y) = c00 
c             + c10 * x                + c01 * y
c             + c20 * x^2   + c11 * xy + c02 * y^2
c             + ...
c             + cm0 * x^(m) + ...      + c0m * y^m.
c
c    If we let T(K) = the K-th triangular number 
c            = sum ( 1 <= I <= K ) I
c    then the number of coefficients in the above polynomial is T(M+1).
c
c    We have n data locations (x(i),y(i)) and values z(i) to approximate:
c
c      p(x(i),y(i)) = z(i)
c
c    This can be cast as an NxT(M+1) linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
c      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
c      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
c      [ ...................... ] [ ... ] = [ ... ]
c      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
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
c    Input, integer TM, the M+1st triangular number.
c
c    Input, double precision X(N), Y(N), the data locations.
c
c    Output, double precision A(N,TM), the Vandermonde matrix for X.
c
      implicit none

      integer n
      integer tm

      double precision a(n,tm)
      integer ex
      integer ey
      integer j
      integer m
      integer s
      double precision x(n)
      double precision y(n)

      j = 0

      do s = 0, m
        do ex = s, 0, -1
          ey = s - ex
          j = j + 1
          a(1:n,j) = x(1:n) ** ex * y(1:n) ** ey 
        end do
      end do
     
      return
      end
