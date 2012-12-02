      subroutine pwl_approx_1d ( nd, xd, yd, nc, xc, yc )

c*****************************************************************************80
c
cc PWL_APPROX_1D determines the control values for a PWL approximant.
c
c  Discussion:
c
c    The piecewise linear approximant is defined by NC control pairs 
c    (XC(I),YC(I)) and approximates ND data pairs (XD(I),YD(I)).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 October 2012
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
c    Input, integer NC, the number of control points.
c    NC must be at least 1.
c
c    Input, double precision XC(NC), the control points.
c
c    Output, double precision YC(NC), the control values.
c
      implicit none

      integer nc
      integer nd

      double precision a(nd,nc)
      double precision xc(nc)
      double precision xd(nd)
      double precision yc(nc)
      double precision yd(nd)
c
c  Define the NDxNC linear system that determines the control values.
c
      call pwl_approx_1d_matrix ( nd, xd, yd, nc, xc, a )
c
c  Solve the system.
c
      call qr_solve ( nd, nc, a, yd, yc )

      return
      end
      subroutine pwl_approx_1d_matrix ( nd, xd, yd, nc, xc, a )

c*****************************************************************************80
c
cc PWL_APPROX_1D_MATRIX returns the matrix for the PWL approximant controls.
c
c  Discussion:
c
c    The value of the piecewise linear approximant, using control points XC
c    and control values YC, evaluated at the point XD, can be represented by
c
c      YD = A * YC
c
c    where A is a matrix whose values depend on XC and XD.
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
c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, double precision YD(ND), the data values.
c
c    Input, integer NC, the number of control points.
c    NC must be at least 1.
c
c    Input, double precision XC(NC), the control points.
c
c    Output, double precision A(ND,NC), the matrix.
c
      implicit none

      integer nc
      integer nd

      double precision a(nd,nc)
      integer i
      integer j
      integer k
      double precision t
      double precision xc(nc)
      double precision xd(nd)
      double precision yd(nd)

      do j = 1, nc
        do i = 1, nd
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, nd
        k = nc - 1
        do j = 2, nc - 1
          if ( xd(i) .lt. xc(j) ) then
            k = j - 1
            go to 10
          end if
        end do
10      continue
        t = ( xd(i) - xc(k) ) / ( xc(k+1) - xc(k) )
        a(i,k)   = 1.0D+00 - t
        a(i,k+1) =           t
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
c    23 September 2012
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
              go to 10

            end if

          end do

10        continue

        end if

      end do
      
      return
      end
