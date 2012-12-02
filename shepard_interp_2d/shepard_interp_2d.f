      subroutine shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, 
     &  zi )

c*********************************************************************72
c
cc SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
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
c    26 September 2012
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
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(ND), YD(ND), the data points.
c
c    Input, double precision ZD(ND), the data values.
c
c    Input, double precision P, the power.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), YI(NI), the interpolation points.
c
c    Output, double precision ZI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni

      integer i
      integer j
      integer k
      double precision p
      double precision r8vec_dot_product
      double precision r8vec_sum
      double precision s
      double precision w(nd)
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)
      integer z
      double precision zd(nd)
      double precision zi(ni)

      do i = 1, ni

        if ( p .eq. 0.0D+00 ) then

          w(1:nd) = 1.0D+00 / dble ( nd )

        else

          z = -1
          do j = 1, nd
            w(j) = sqrt ( ( xi(i) - xd(j) ) ** 2 
     &                  + ( yi(i) - yd(j) ) ** 2 )
            if ( w(j) .eq. 0.0D+00 ) then
              z = j
              exit
            end if
          end do

          if ( z .ne. -1 ) then
            do k = 1, nd
              w(k) = 0.0D+00
            end do
            w(z) = 1.0D+00
          else
            do k = 1, nd
              w(k) = 1.0D+00 / w(k) ** p
            end do
            s = r8vec_sum ( nd, w )
            do k = 1, nd
              w(k) = w(k) / s
            end do
          end if

        end if

        zi(i) = r8vec_dot_product ( nd, w, zd )

      end do

      return
      end
