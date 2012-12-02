      subroutine nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
c
c  Discussion:
c
c    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
c    constant function which interpolates the data (XD(I),YD(I)) for I = 1
c    to ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 September 2012
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

      double precision d
      double precision d2
      integer i
      integer j
      integer k
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)

      do i = 1, ni

        k = 1
        d = abs ( xi(i) - xd(k) )

        do j = 2, nd

          d2 = abs ( xi(i) - xd(j) )

          if ( d2 .lt. d ) then
            k = j
            d = d2
          end if

        end do

        yi(i) = yd(k)

      end do

      return
      end
