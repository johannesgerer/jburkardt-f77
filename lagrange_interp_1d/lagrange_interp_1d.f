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
c    12 September 2012
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
      double precision xi(ni)
c
c  Evaluate the polynomial.
c
      do j = 1, nd
        do i = 1, ni
          lb(i,j) = 1.0D+00
        end do
      end do

      do k = 1, nd

        do j = 1, nd

          if ( j .ne. k ) then

            do i = 1, ni
              lb(i,k) = lb(i,k) * ( xi(i) - xd(j) ) / ( xd(k) - xd(j) )
            end do

          end if

        end do

      end do

      return
      end
      subroutine lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
c
c  Discussion:
c
c    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
c    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
c    to ND.
c
c    The Lagrange interpolant can be constructed from the Lagrange basis
c    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange 
c    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree 
c    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
c
c    Given data values YD at each of the abscissas, the value of the
c    Lagrange interpolant may be written as
c
c      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 September 2012
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
c    Input, double precision YD(ND), the data values.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision YI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni

      integer i
      integer j
      double precision lb(ni,nd)
      double precision xd(nd)
      double precision yd(nd)
      double precision xi(ni)
      double precision yi(ni)

      call lagrange_basis_1d ( nd, xd, ni, xi, lb )

      call r8mat_mv ( ni, nd, lb, yd, yi )

      return
      end
