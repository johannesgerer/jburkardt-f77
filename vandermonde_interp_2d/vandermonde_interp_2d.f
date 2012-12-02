      function triangle_num ( n )

c*********************************************************************72
c
cc TRIANGLE_NUM returns the N-th triangular number.
c
c  Discussion:
c
c    The N-th triangular number T(N) is formed by the sum of the first
c    N integers:
c
c      T(N) = sum ( 1 <= I <= N ) I
c
c    By convention, T(0) = 0.
c
c    T(N) can be computed quickly by the formula:
c
c      T(N) = ( N * ( N + 1 ) ) / 2
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
c    11 August 1998
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
      subroutine vandermonde_interp_2d_matrix ( n, m, x, y, a )

c*********************************************************************72
c
cc VANDERMONDE_INTERP_2D_MATRIX computes a Vandermonde 2D interpolation matrix.
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
c    and we assume that N = T(M+1).
c
c    This can be cast as an NxN linear system for the polynomial
c    coefficients:
c
c      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
c      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
c      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
c      [ ...................... ] [ ... ] = [ ... ]
c      [ 1 xn yn  xn^2 ... yn^m ] [ c0n ] = [  zn ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.  It is necessary 
c    that N = T(M+1), where T(K) is the K-th triangular number.
c
c    Input, integer M, the degree of the polynomial.
c
c    Input, double precision X(N), Y(N), the data locations.
c
c    Output, double precision A(N,N), the Vandermonde matrix for X.
c
      implicit none

      integer n

      double precision a(n,n)
      integer ex
      integer ey
      integer i
      integer j
      integer m
      integer s
      integer tmp1
      integer triangle_num
      double precision x(n)
      double precision y(n)

      tmp1 = triangle_num ( m + 1 )

      if ( n .ne. tmp1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_MATRIX - Fatal error!'
        write ( *, '(a)' ) '  For interpolation, we need N = T(M+1).'
        write ( *, '(a,i6)' ) '  But we have N = ', n
        write ( *, '(a,i6)' ) '  M = ', m
        write ( *, '(a,i6)' ) '  and T(M+1) = ', tmp1
        stop
      end if

      j = 0

      do s = 0, m
        do ex = s, 0, -1
          ey = s - ex
          j = j + 1
          do i = 1, n
            a(i,j) = x(i) ** ex * y(i) ** ey
          end do
        end do
      end do
     
      return
      end
