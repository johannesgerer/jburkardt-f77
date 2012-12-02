      subroutine lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc LAGCHEBY1_INTERP_1D evaluates the Lagrange Chebyshev 1 interpolant.
c
c  Discussion:
c
c    The weight vector WD computed below is only valid if the data points
c    XD are, as expected, the Chebyshev Type 1 points for [-1,+1], or a linearly 
c    mapped version for [A,B].  The XD values may be computed by:
c
c      xd = r8vec_cheby1space ( nd, a, b )
c
c    for instance.
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
c  Reference:
c
c    Jean-Paul Berrut, Lloyd Trefethen,
c    Barycentric Lagrange Interpolation,
c    SIAM Review,
c    Volume 46, Number 3, September 2004, pages 501-517.
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

      double precision denom(ni)
      integer exact(ni)
      integer i
      integer j
      double precision numer(ni)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_mop
      double precision t
      double precision theta
      double precision wd
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)

      do i = 1, ni
        exact(i) = 0
      end do

      do j = 1, nd

        theta = dble ( 2 * j - 1 ) * pi / dble ( 2 * nd )
        wd = r8_mop ( j + 1 ) * sin ( theta )

        do i = 1, ni

          if ( xi(i) .eq. xd(j) ) then
            exact(i) = j
            numer(i) = yd(j)
            denom(i) = 1.0D+00
          end if

          if ( exact(i) .eq. 0 ) then
            t = wd / ( xi(i) - xd(j) )
            numer(i) = numer(i) + t * yd(j)
            denom(i) = denom(i) + t
          end if

        end do
      end do

      do i = 1, ni
        yi(i) = numer(i) / denom(i)
      end do

      return
      end
      subroutine lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc LAGCHEBY2_INTERP_1D evaluates the Lagrange Chebyshev 2 interpolant.
c
c  Discussion:
c
c    The weight vector WD computed below is only valid if the data points
c    XD are, as expected, the Chebyshev Type 2 points for [-1,+1], or a linearly 
c    mapped version for [A,B].  The XD values may be computed by:
c
c      xd = r8vec_cheby2space ( nd, a, b )
c
c    for instance.
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
c  Reference:
c
c    Jean-Paul Berrut, Lloyd Trefethen,
c    Barycentric Lagrange Interpolation,
c    SIAM Review,
c    Volume 46, Number 3, September 2004, pages 501-517.
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

      double precision denom(ni)
      integer exact(ni)
      integer i
      integer j
      double precision numer(ni)
      double precision r8_mop
      double precision t
      double precision wd
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)

      do i = 1, ni
        exact(i) = 0
      end do

      do j = 1, nd

        wd = r8_mop ( j + 1 )
        if ( j .eq. 1 .or. j .eq. nd ) then
          wd = 0.5D+00 * wd
        end if

        do i = 1, ni

          if ( xi(i) .eq. xd(j) ) then
            exact(i) = j
            numer(i) = yd(j)
            denom(i) = 1.0D+00
          end if

          if ( exact(i) .eq. 0 ) then
            t = wd / ( xi(i) - xd(j) )
            numer(i) = numer(i) + t * yd(j)
            denom(i) = denom(i) + t
          end if

        end do
      end do

      do i = 1, ni
        yi(i) = numer(i) / denom(i)
      end do

      return
      end
      subroutine lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

c*********************************************************************72
c
cc LAGEVEN_VALUE_1D evaluates the Lagrange evenly-spaced interpolant.
c
c  Discussion:
c
c    The weight vector WD computed below is only valid if the data points
c    XD are, as expected, evenly spaced in an interval [A,B] with
c    spacing (B-A)/N.  The XD values might be computed by:
c
c      xd(i) = ( ( 2 * nd - 2 * i + 1 ) * a 
c              + (          2 * i - 1 ) * b ) 
c              / ( 2 * nd             )
c
c    for instance.
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
c  Reference:
c
c    Jean-Paul Berrut, Lloyd Trefethen,
c    Barycentric Lagrange Interpolation,
c    SIAM Review,
c    Volume 46, Number 3, September 2004, pages 501-517.
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

      double precision denom(ni)
      integer exact(ni)
      integer i
      integer j
      double precision numer(ni)
      double precision r8_choose
      double precision r8_mop
      double precision t
      double precision wd
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)

      do i = 1, ni
        exact(i) = 0
      end do

      do j = 1, nd

        wd = r8_mop ( j ) * r8_choose ( nd, j )

        do i = 1, ni

          if ( xi(i) .eq. xd(j) ) then
            exact(i) = j
            numer(i) = yd(j)
            denom(i) = 1.0D+00
          end if

          if ( exact(i) .eq. 0 ) then
            t = wd / ( xi(i) - xd(j) )
            numer(i) = numer(i) + t * yd(j)
            denom(i) = denom(i) + t
          end if

        end do
      end do

      do i = 1, ni
        yi(i) = numer(i) / denom(i)
      end do

      return
      end
