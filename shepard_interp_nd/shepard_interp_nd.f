      subroutine shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

c*********************************************************************72
c
cc SHEPARD_INTERP_ND evaluates a Shepard interpolant in M dimensions.
c
c  Discussion:
c
c    This code should be vectorized.
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
c  Reference:
c
c    Donald Shepard,
c    A two-dimensional interpolation function for irregularly spaced data,
c    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
c    ACM, pages 517-524, 1969.
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(M,ND), the data points.
c
c    Input, double precision ZD(ND), the data values.
c
c    Input, double precision P, the power.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(M,NI), the interpolation points.
c
c    Output, double precision ZI(NI), the interpolated values.
c
      implicit none

      integer m
      integer nd
      integer ni

      integer i
      integer i2
      integer j
      double precision p
      double precision r8vec_dot_product
      double precision r8vec_norm_affine
      double precision r8vec_sum
      double precision s
      double precision w(nd)
      double precision xd(m,nd)
      double precision xi(m,nd)
      integer z
      double precision zd(nd)
      double precision zi(ni)

      do i = 1, ni

        if ( p .eq. 0.0D+00 ) then

          do i2 = 1, nd
            w(i2) = 1.0D+00 / dble ( nd )
          end do

        else

          z = -1
          do j = 1, nd
            w(j) = r8vec_norm_affine ( m, xi(1,i), xd(1,j) )
            if ( w(j) .eq. 0.0D+00 ) then
              z = j
              go to 10
            end if
          end do

10    continue

          if ( z .ne. -1 ) then
            do i2 = 1, nd
              w(i2) = 0.0D+00
            end do
            w(z) = 1.0D+00
          else
            do i2 = 1, nd
              w(i2) = 1.0D+00 / w(i2) ** p
            end do
            s = r8vec_sum ( nd, w )
            do i2 = 1, nd
              w(i2) = w(i2) / s
            end do
          end if

        end if

        zi(i) = r8vec_dot_product ( nd, w, zd )

      end do

      return
      end
