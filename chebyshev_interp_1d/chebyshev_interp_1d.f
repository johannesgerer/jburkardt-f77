      subroutine chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

c*********************************************************************72
c
cc CHEBYSHEV_COEF_1D determines the Chebyshev interpolant coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 September 2012
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
c    Input, double precision XD(ND), the data locations.
c
c    Input, double precision YD(ND), the data values.
c
c    Output, double precision C(ND), the Chebyshev coefficients.
c
c    Output, double precision XMIN, XMAX, the interpolation interval.
c
      implicit none

      integer nd

      double precision a(nd,nd)
      double precision c(nd)
      integer i
      integer j
      double precision x(nd)
      double precision xd(nd)
      double precision xmax
      double precision xmin
      double precision yd(nd)

      if ( nd .eq. 1 ) then
        xmin = xd(1)
        xmax = xd(1)
        c(1) = 1.0D+00
        return
      end if

      call r8vec_min ( nd, xd, xmin )
      call r8vec_max ( nd, xd, xmax )
c
c  Map XD to [-1,+1].
c
      do i = 1, nd
        x(i) = ( 2.0D+00 * xd(i) - xmin - xmax ) / ( xmax - xmin )
      end do
c
c  Form the Chebyshev Vandermonde matrix.
c
      do j = 1, nd
        do i = 1, nd
          a(i,j) = cos ( acos ( x(i) ) * dble ( j - 1 ) )
        end do
      end do 
c
c  Solve for the expansion coefficients.
c
      call qr_solve ( nd, nd, a, yd, c )

      return
      end
      subroutine chebyshev_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc CHEBYSHEV_INTERP_1D determines and evaluates the Chebyshev interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 September 2012
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
c    Input, double precision XD(ND), the data locations.
c
c    Input, double precision YD(ND), the data values.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points, which
c    must be each be in the interval [ min(XD), max(XD)].
c
c    Output, double precision YI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni

      double precision c(nd)
      double precision xd(nd)
      double precision xi(ni)
      double precision xmax
      double precision xmin
      double precision yd(nd)
      double precision yi(ni)

      call chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

      call chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

      return
      end
      subroutine chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

c*********************************************************************72
c
cc CHEBYSHEV_VALUE_1D evaluates a Chebyshev interpolant, given its coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 September 2012
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
c    Input, double precision C(ND), the Chebyshev coefficients.
c
c    Input, double precision XMIN, XMAX, the interpolation interval.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points, which
c    must be each be in the interval [XMIN,XMAX].
c
c    Output, double precision YI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni

      double precision a(ni,nd)
      double precision c(nd)
      integer i
      integer j
      double precision x(ni)
      double precision xi(ni)
      double precision xmax
      double precision xmin
      double precision yi(ni)

      if ( nd .eq. 1 ) then
        yi(1) = c(1)
        return
      end if
c
c  Map XI to [-1,+1].
c
      do i = 1, ni
        x(i) = ( 2.0D+00 * xi(i) - xmin - xmax ) / ( xmax - xmin )
      end do

      do j = 1, nd
        do i = 1, ni
          a(i,j) = cos ( acos ( x(i) ) * dble ( j - 1 ) )
        end do
      end do 

      call r8mat_mv ( ni, nd, a, c, yi )

      return
      end
