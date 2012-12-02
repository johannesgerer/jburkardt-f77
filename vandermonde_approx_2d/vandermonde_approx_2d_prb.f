      program main

c*********************************************************************72
c
cc VANDERMONDE_APPROX_2D_TEST tests VANDERMONDE_APPROX_2D.
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

      integer m_test_num
      parameter ( m_test_num = 5 )

      integer j
      integer m
      integer m_test(m_test_num)
      integer grid
      integer prob
      integer prob_num

      save m_test

      data m_test / 0, 1, 2, 4, 8 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_APPROX_2D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the VANDERMONDE_APPROX_2D library.'
      write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )

      do prob = 1, prob_num
        grid = 1
        do j = 1, m_test_num
          m = m_test(j)
          call test01 ( prob, grid, m )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_APPROX_2D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, grd, m )

c*********************************************************************72
c
cc VANDERMONDE_APPROX_2D_TEST01 tests VANDERMONDE_APPROX_2D_MATRIX.
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
c    Input, integer PROB, the problem number.
c
c    Input, integer GRD, the grid number.
c
c    Input, integer M, the total polynomial degree.
c
      implicit none

      integer m
      integer nd_max
      parameter ( nd_max = 100 )
      integer ni_max
      parameter ( ni_max = 100 )

      double precision a(nd_max,(m+1)*(m+2)/2)
      double precision app_error
      double precision c((m+1)*(m+2)/2)
      integer grd
      integer i 
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      integer tm
      integer triangle_num
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision yd(nd_max)
      double precision yi(ni_max)
      double precision zd(nd_max)
      double precision zi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i4)' ) 
     &  '  Approximate data from TEST_INTERP_2D problem #', prob
      write ( *, '(a,i4)' ) 
     &  '  Use grid from TEST_INTERP_2D with index #', grd
      write ( *, '(a,i4)' ) 
     &  '  Using polynomial approximant of total degree ', m

      call g00_size ( grd, nd )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call g00_xy ( grd, nd, xd, yd )

      call f00_f0 ( prob, nd, xd, yd, zd )

      if ( nd .lt. 10 ) then
        call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
      end if
c
c  Compute the Vandermonde matrix.
c
      tm = triangle_num ( m + 1 );
      call vandermonde_approx_2d_matrix ( nd, m, tm, xd, yd, a )
c
c  Solve linear system.
c
      call qr_solve ( nd, tm, a, zd, c )
c
c  #1:  Does approximant match function at data points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
        yi(i) = yd(i)
      end do

      call r8poly_value_2d ( m, c, ni, xi, yi, zi )

      app_error = r8vec_norm_affine ( ni, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 data approximation error = ', app_error

      return
      end
