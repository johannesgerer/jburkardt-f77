      subroutine lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc LAGRANGE_APPROX_1D evaluates the Lagrange approximant of degree M.
c
c  Discussion:
c
c    The Lagrange approximant L(M,ND,XD,YD)(X) is a polynomial of
c    degree M which approximates the data (XD(I),YD(I)) for I = 1 to ND.
c
c    We can represent any polynomial of degree M+1 as the sum of the Lagrange 
c    basis functions at the M+1 Chebyshev points.
c
c      L(M)(X) = sum ( 1 <= I <= M+1 ) C(I) LB(M,XC)(X)
c
c    Given our data, we can seek the M+1 unknown coefficients C which minimize
c    the norm of || L(M)(XD(1:ND)) - YD(1:ND) ||.
c
c    Given the coefficients, we can then evaluate the polynomial at the
c    points XI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the polynomial degree.
c
c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, double precision YD(ND), the data values.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision YI(NI), the interpolated values.
c
      implicit none

      integer m
      integer nd
      integer ni

      double precision a
      double precision b
      double precision ld(nd,m+1)
      double precision li(ni,m+1)
      integer nc
      double precision xc(m+1)
      double precision xd(nd)
      double precision xi(ni)
      double precision yc(m+1)
      double precision yd(nd)
      double precision yi(ni)

      nc = m + 1
c
c  Evaluate the Chebyshev points.
c
      a = -1.0D+00
      b = +1.0D+00
      call r8vec_chebyspace ( nc, a, b, xc )
c
c  Evaluate the Lagrange basis functions for the Chebyshev points 
c  at the data points.
c
      call lagrange_basis_1d ( nc, xc, nd, xd, ld )
c
c  The value of the Lagrange approximant at each data point should
c  approximate the data value: LD * YC = YD, where YC are the unknown
c  coefficients.
c
      call qr_solve ( nd, nc, ld, yd, yc )
c
c  Now we want to evaluate the Lagrange approximant at the "interpolant
c  points": LI * YC = YI
c
      call lagrange_basis_1d ( nc, xc, ni, xi, li )

      call r8mat_mv ( ni, nc, li, yc, yi )

      return
      end
      subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb )

c*********************************************************************72
c
cc LAGRANGE_BASIS_1D evaluates the Lagrange basis polynomials.
c
c  Discussion:
c
c    Given ND distinct abscissas, XD(1:ND),
c    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
c    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
c    other abscissas.
c
c    A formal representation is:
c
c      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
c       ( T - T(J) ) / ( T(I) - T(J) )
c
c    This routine accepts a set of NI values at which all the Lagrange
c    basis polynomials should be evaluated.
c
c    Given data values YD at each of the abscissas, the value of the
c    Lagrange interpolating polynomial at each of the interpolation points
c    is then simple to compute by matrix multiplication:
c
c      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision LB(NI,ND), the values
c    of the Lagrange basis polynomials at the interpolation points.
c
      implicit none

      integer nd
      integer ni

      integer i
      integer j
      integer k
      double precision lb(ni,nd)
      double precision xd(nd)
      double precision xi(nd)
c
c  Evaluate the polynomial.
c
      do j = 1, nd
        do k = 1, ni
          lb(k,j) = 1.0D+00
        end do
      end do

      do i = 1, nd

        do j = 1, nd

          if ( j .ne. i ) then

            do k = 1, ni
              lb(k,i) = lb(k,i) * ( xi(k) - xd(j) ) / ( xd(i) - xd(j) )
            end do

          end if

        end do

      end do

      return
      end
