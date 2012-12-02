      program main

c*********************************************************************72
c
cc RBF_INTERP_2D_TEST tests RBF_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer g
      external phi1
      external phi2
      external phi3
      external phi4
      integer prob
      integer prob_num

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RBF_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the RBF_INTERP_2D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )
      g = 1

      do prob = 1, prob_num
        call test01 ( prob, g, phi1, 'phi1' )
        call test01 ( prob, g, phi2, 'phi2' )
        call test01 ( prob, g, phi3, 'phi3' )
        call test01 ( prob, g, phi4, 'phi4' )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RBF_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, g, phi, phi_name )

c*********************************************************************72
c
cc RBF_INTERP_2D_TEST01 tests RBF_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the problem.
c
c    Input, integer G, the index of the grid.
c
c    Input, external PHI, the radial basis function.
c
c    Input, character * ( * ) PHI_NAME, the name of the radial basis function.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 100 )
      integer ni_max
      parameter ( ni_max = 100 )

      double precision e
      integer g
      integer i
      double precision int_error
      integer j
      integer m
      integer nd
      integer ni
      external phi
      integer prob
      character * ( * ) phi_name
      double precision r0
      double precision r8vec_norm_affine
      double precision volume
      double precision w(nd_max)
      double precision xd(nd_max)
      double precision xmax
      double precision xmin
      double precision xyd(2,nd_max)
      double precision xyi(2,ni_max)
      double precision yd(nd_max)
      double precision ymax
      double precision ymin
      double precision zd(nd_max)
      double precision zi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RBF_INTERP_2D_TEST01:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP_2D problem #', prob
      write ( *, '(a,i6)' ) '  using grid #', g
      write ( *, '(a)' ) 
     &  '  using radial basis function "' // trim ( phi_name ) // '".'

      call g00_size ( g, nd )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call g00_xy ( g, nd, xd, yd )

      call f00_f0 ( prob, nd, xd, yd, zd )

      if ( phi_name .eq. 'phi1' ) then
        call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
      end if

      m = 2

      do i = 1, nd
        xyd(1,i) = xd(i)
        xyd(2,i) = yd(i)
      end do

      call r8vec_max ( nd, xd, xmax )
      call r8vec_min ( nd, xd, xmin )
      call r8vec_max ( nd, yd, ymax )
      call r8vec_min ( nd, yd, ymin )
      volume = ( xmax - xmin ) * ( ymax - ymin )

      e = 1.0D+00 / dble ( m )
      r0 = ( volume / nd ) ** e

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Setting R0 = ', r0

      call rbf_weight ( m, nd, xyd, r0, phi, zd, w )
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do j = 1, ni
        do i = 1, 2
          xyi(i,j) = xyd(i,j)
        end do
      end do

      call rbf_interp ( m, nd, xyd, r0, phi, w, ni, xyi, zi )

      int_error = r8vec_norm_affine ( ni, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
