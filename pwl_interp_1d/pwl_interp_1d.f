      subroutine pwl_basis_1d ( nd, xd, k, ni, xi, bk )

c*****************************************************************************80
c
cc PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, integer K, the index of the desired basis function,
c    1 .le. K .le. ND.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision BK(NI), the basis function at the 
c    interpolation points.
c
      implicit none

      integer nd
      integer ni

      double precision bk(ni)
      integer i
      integer k
      double precision t
      double precision xd(nd)
      double precision xi(ni)

      do i = 1, ni
        bk(i) = 0.0D+00
      end do

      if ( nd .eq. 1 ) then
        do i = 1, ni
          bk(i) = 1.0D+00
        end do
        return
      end if

      do i = 1, ni

        if ( k .eq. 1 .and. xi(i) .le. xd(k) ) then

          t = ( xi(i) - xd(k) ) / ( xd(k+1) - xd(k) )
          bk(i) = 1.0D+00 - t

        else if ( k .eq. nd .and. xd(k) .le. xi(i) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          bk(i) = t

        else if ( xd(k-1) .le. xi(i) .and. xi(i) .le. xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          bk(i) = t

        else if ( xd(k) .le. xi(i) .and. xi(i) .le. xd(k+1) ) then

          t = ( xi(i) - xd(k) ) / ( xd(k+1) - xd(k) )
          bk(i) = 1.0D+00 - t

        end if

      end do

      return
      end
      subroutine pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

c*****************************************************************************80
c
cc PWL_INTERP_1D evaluates the piecewise linear interpolant.
c
c  Discussion:
c
c    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
c    linear function which interpolates the data (XD(I),YD(I)) for I = 1
c    to ND.
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
      integer k
      double precision t
      double precision xd(nd)
      double precision yd(nd)
      double precision xi(ni)
      double precision yi(ni)

      do i = 1, ni
        yi(i) = 0.0D+00
      end do

      if ( nd .eq. 1 ) then
        do i = 1, ni
          yi(i) = yd(1)
        end do
        return
      end if

      do i = 1, ni

        if ( xi(i) .le. xd(1) ) then

          t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
          yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

        else if ( xd(nd) .le. xi(i) ) then

          t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

        else

          do k = 2, nd

            if ( xd(k-1) .le. xi(i) .and. xi(i) .le. xd(k) ) then

              t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
              yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
              exit

            end if

          end do

        end if

      end do
      
      return
      end
