      subroutine lagrange_basis_function_1d ( mx, xd, i, xi, yi ) 

c*********************************************************************72
c
cc LAGRANGE_BASIS_FUNCTION_1D evaluates a 1D Lagrange basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MX, the degree of the basis function.
c
c    Input, double precision XD(MX+1), the interpolation nodes.
c
c    Input, integer I, the index of the basis function.
c    1 <= I <= MX+1.
c
c    Input, double precision XI, the evaluation point.
c
c    Output, double precision YI, the value of the I-th Lagrange 1D basis 
c    function for the nodes XD, evaluated at XI.
c
      implicit none

      integer mx

      integer i
      integer j
      double precision xd(mx+1)
      double precision xi
      double precision yi

      yi = 1.0D+00

      if ( xi .ne. xd(i) ) then
        do j = 1, mx + 1
          if ( j .ne. i ) then
            yi = yi * ( xi - xd(j) ) / ( xd(i) - xd(j) )
          end if
        end do
      end if

      return
      end
      subroutine lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, 
     &  xi, yi, zi )

c*********************************************************************72
c
cc LAGRANGE_INTERP_2D evaluates the Lagrange interpolant for a product grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MX, MY, the polynomial degree in X and Y.
c
c    Input, double precision XD_1D(MX+1), YD_1D(MY+1), the 1D data locations.
c
c    Input, double precision ZD((MX+1),(MY+1)), the 2D array of data values.
c
c    Input, integer NI, the number of 2D interpolation points.
c
c    Input, double precision XI(NI), YI(NI), the 2D interpolation points.
c
c    Output, double precision ZI(NI), the interpolated values.
c
      implicit none

      integer mx
      integer my
      integer ni

      integer i
      integer j
      integer k
      integer l
      double precision lx
      double precision ly
      double precision xd_1d(mx+1)
      double precision xi(ni)
      double precision yd_1d(my+1)
      double precision yi(ni)
      double precision zd(mx+1,my+1)
      double precision zi(ni)

      do k = 1, ni
        l = 0
        zi(k) = 0.0D+00
        do i = 1, mx + 1
          do j = 1, my + 1
            l = l + 1
            call lagrange_basis_function_1d ( mx, xd_1d, i, xi(k), lx )
            call lagrange_basis_function_1d ( my, yd_1d, j, yi(k), ly )
            zi(k) = zi(k) + zd(i,j) * lx * ly
          end do
        end do
      end do

      return
      end
