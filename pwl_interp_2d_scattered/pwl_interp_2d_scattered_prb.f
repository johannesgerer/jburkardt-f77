      program main

c*********************************************************************72
c
cc PWL_INTERP_2D_SCATTERED_PRB tests PWL_INTERP_2D_SCATTERED_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prob
      integer prob_num


      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PWL_INTERP_2D_SCATTERED_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test PWL_INTERP_2D_SCATTERED.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call test01 ( )
      call test02 ( )
c
c  Numerical tests.
c
      call f00_num ( prob_num )

      do prob = 1, prob_num
        call test03 ( prob )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PWL_INTERP_2D_SCATTERED_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R8TRIS2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 9 )
      integer element_order
      parameter ( element_order = 3 )

      integer element_neighbor(3,2*node_num)
      integer element_num
      double precision node_xy(dim_num,node_num)
      integer triangle(element_order,2*node_num)

      save node_xy

      data node_xy /
     &     0.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     0.2D+00, 0.5D+00, 
     &     0.3D+00, 0.6D+00, 
     &     0.4D+00, 0.5D+00, 
     &     0.6D+00, 0.4D+00, 
     &     0.6D+00, 0.5D+00, 
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  R8TRIS2 computes the Delaunay triangulation'
      write ( *, '(a)' ) '  of a set of nodes in 2D.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, element_num, triangle, 
     &  element_neighbor )

      call triangulation_order3_print ( node_num, element_num, node_xy, 
     &  triangle, element_neighbor )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8TRIS2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 August 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer ni
      parameter ( ni = 25 )
      integer node_num
      parameter ( node_num = 9 )
      integer element_order
      parameter ( element_order = 3 )

      integer element_neighbor(3,2*node_num)
      integer element_num
      integer i
      integer j
      integer k
      double precision node_xy(dim_num,node_num)
      integer triangle(element_order,2*node_num)
      double precision x
      double precision xyi(2,ni)
      double precision y
      double precision zd(node_num)
      double precision zi(ni)

      save node_xy

      data node_xy /
     &     0.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     0.2D+00, 0.5D+00, 
     &     0.3D+00, 0.6D+00, 
     &     0.4D+00, 0.5D+00, 
     &     0.6D+00, 0.4D+00, 
     &     0.6D+00, 0.5D+00, 
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  PWL_INTERP_2D_SCATTERED_VALUE evaluates a'
      write ( *, '(a)' ) '  piecewise linear interpolant of scattered'
      write ( *, '(a)' ) '  data.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, element_num, triangle, 
     &  element_neighbor )
c
c  Define the Z data.
c
      do i = 1, node_num
        x = node_xy(1,i)
        y = node_xy(2,i)
        zd(i) = x + 2.0D+00 * y
      end do
c
c  Define the interpolation points.
c
      k = 0
      do i = 1, 5
        do j = 1, 5
          k = k + 1
          xyi(1,k) = dble ( i - 1 ) / 4.0D+00
          xyi(2,k) = dble ( j - 1 ) / 4.0D+00
        end do
      end do
c
c  Evaluate the interpolant.
c
      call pwl_interp_2d_scattered_value ( node_num, node_xy, zd, 
     &  element_num, triangle, element_neighbor, ni, xyi, zi )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'K Xi(K)  Yi(K)  Zi(K)   Z(X,Y)'
      write ( *, '(a)' ) ' '
      do k = 1, ni
        write ( *, '(2x,i4,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2)' ) 
     &    k, xyi(1,k), xyi(2,k), zi(k), xyi(1,k) + 2.0 * xyi(2,k)
      end do

      return
      end
      subroutine test03 ( prob )

c*********************************************************************72
c
cc TEST03 tests PWL_INTERP_2D_SCATTERED_VALUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd_max
      parameter ( nd_max = 33 )
      integer ni
      parameter ( ni = 25 )

      integer element_neighbor(3,2*nd_max)
      integer element_num
      integer g
      integer i
      integer j
      integer k
      integer nd
      integer prob
      double precision rms
      integer triangle(3,2*nd_max)
      double precision x
      double precision xd(nd_max)
      double precision xi(ni)
      double precision xyd(2,nd_max)
      double precision xyi(2,ni)
      double precision y
      double precision yd(nd_max)
      double precision yi(ni)
      double precision zd(nd_max)
      double precision ze(ni)
      double precision zi(ni)

      g = 2
      call g00_size ( g, nd )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  PWL_INTERP_2D_SCATTERED_VALUE evaluates a'
      write ( *, '(a)' ) 
     &  '  piecewise linear interpolant to scattered data.'
      write ( *, '(a,i2)' ) '  Here, we use grid number ', g
      write ( *, '(a,i4,a)' ) 
     &  '  with ', nd, ' scattered points in the unit square'
      write ( *, '(a,i2)' ) '  on problem ', prob
c
c  Get the data points and evaluate the function there.
c
      call g00_xy ( g, nd, xd, yd )

      call f00_f0 ( prob, nd, xd, yd, zd )

      do j = 1, nd
        xyd(1,j) = xd(j)
        xyd(2,j) = yd(j)
      end do
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( nd, xyd, element_num, triangle, 
     &  element_neighbor )
c
c  Define the interpolation points.
c
      k = 0
      do i = 1, 5
        do j = 1, 5
          k = k + 1
          xyi(1,k) = dble ( 2 * i - 1 ) / 10.0D+00
          xyi(2,k) = dble ( 2 * j - 1 ) / 10.0D+00
        end do
      end do

      do j = 1, ni
        xi(j) = xyi(1,j)
        yi(j) = xyi(2,j)
      end do

      call f00_f0 ( prob, ni, xi, yi, ze )
c
c  Evaluate the interpolant.
c
      call pwl_interp_2d_scattered_value ( nd, xyd, zd, 
     &  element_num, triangle, element_neighbor, ni, xyi, zi )


      rms = 0.0D+00
      do k = 1, ni
        rms = rms + ( zi(k) - ze(k) )**2
      end do
      rms = sqrt ( rms / dble ( ni ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,e10.2)' ) '  RMS error is ', rms

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)'
      write ( *, '(a)' ) ' '

      do k = 1, ni
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    k, xyi(1,k), xyi(2,k), zi(k), ze(k)
      end do

      return
      end
