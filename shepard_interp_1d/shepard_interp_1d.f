      subroutine shepard_basis_1d ( nd, xd, k, p, ni, xi, bk )

c*********************************************************************72
c
cc SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 September 2012
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
c    Input, double precision XD(ND), the data points.
c
c    Input, integer K, the index of the desired basis function,
c    1 <= K <= ND.
c
c    Input, double precision P, the power.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision BK(NI), the basis function at the interpolation 
c    points.
c 
      implicit none

      integer nd
      integer ni

      double precision bk(ni)
      integer i
      integer j
      integer k
      double precision p
      double precision s
      double precision w(nd)
      double precision xd(nd)
      double precision xi(ni)
      integer z

      do i = 1, ni

        if ( p .eq. 0.0D+00 ) then

          do j = 1, nd
            w(j) = 1.0D+00 / dble ( nd )
          end do

        else

          z = -1
          do j = 1, nd
            w(j) = abs ( xi(i) - xd(j) )
            if ( w(j) .eq. 0.0D+00 ) then
              z = j
              go to 10
            end if
          end do

10        continue

          if ( z .ne. -1 ) then
            do j = 1, nd
              w(j) = 0.0D+00
            end do
            w(z) = 1.0D+00
          else
            do j = 1, nd
              w(j) = 1.0D+00 / w(j) ** p
            end do
            s = 0.0D+00
            do j = 1, nd
              s = s + w(j)
            end do
            do j = 1, nd
              w(j) = w(j) / s
            end do
          end if

        end if

        bk(i) = w(k)

      end do

      return
      end
      subroutine shepard_interp_1d ( nd, xd, yd, p, ni, xi, yi )

c*********************************************************************72
c
cc SHEPARD_INTERP_1D evaluates a 1D Shepard interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 September 2012
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
c    Input, double precision XD(ND), the data points.
c
c    Input, double precision YD(ND), the data values.
c
c    Input, double precision P, the power.
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
      integer k
      double precision p
      double precision r8vec_dot_product
      double precision s
      double precision w(nd)
      double precision xd(nd)
      double precision xi(ni)
      double precision yd(nd)
      double precision yi(ni)
      integer z

      do i = 1, ni

        if ( p .eq. 0.0D+00 ) then

          do j = 1, nd
            w(j) = 1.0D+00 / dble ( nd )
          end do

        else

          z = 0
          do j = 1, nd
            w(j) = abs ( xi(i) - xd(j) )
            if ( w(j) .eq. 0.0D+00 ) then
              z = j
              go to 10
            end if
          end do

10        continue

          if ( z .ne. 0 ) then
            do j = 1, nd
              w(j) = 0.0D+00
            end do
            w(z) = 1.0D+00
          else
            do j = 1, nd
              w(j) = 1.0D+00 / w(j) ** p
            end do
            s = 0.0D+00
            do j = 1, nd
              s = s + w(j)
            end do
            do j = 1, nd
              w(j) = w(j) / s
            end do
          end if

        end if

        yi(i) = r8vec_dot_product ( nd, w, yd )

      end do

      return
      end
