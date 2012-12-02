      program main

c*********************************************************************72
c
cc VANDERMONDE_INTERP_2D_TPRB tests VANDERMONDE_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2012
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
      integer prob
      integer prob_num

      save m_test

      data m_test / 1, 2, 3, 4, 8 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_PRB:'
      write ( *, '(a)' ) '  MATLAB version'
      write ( *, '(a)' ) '  Test the VANDERMONDE_INTERP_2D library.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )
      do prob = 1, prob_num
        do j = 1, m_test_num
          m = m_test(j)
          call test01 ( prob, m )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_INTERP_2D_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, m )

c*********************************************************************72
c
cc TEST01 tests VANDERMONDE_INTERP_2D_MATRIX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer M, the degree of interpolation.
c
      implicit none

      integer m_max
      parameter ( m_max = 8 )
      integer nd_max
      parameter ( nd_max = ( m_max + 1 ) * ( m_max + 2 ) / 2 )
      integer ni_max
      parameter ( ni_max = ( m_max + 1 ) * ( m_max + 2 ) / 2 )

      double precision a(nd_max,nd_max)
      double precision app_error
      double precision c(nd_max)
      integer i
      integer m
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      integer seed
      integer tmp1
      integer triangle_num
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision yd(nd_max)
      double precision yi(ni_max)
      double precision zd(nd_max)
      double precision zi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP_2D problem #', prob
      write ( *, '(a,i6)' ) 
     &  '  Create an interpolant of total degree ', m
      tmp1 = triangle_num ( m + 1 )
      write ( *, '(a,i6)' ) '  Number of data values needed is', tmp1

      nd = tmp1

      seed = 123456789

      call r8vec_uniform_01 ( nd, seed, xd )
      call r8vec_uniform_01 ( nd, seed, yd )

      call f00_f0 ( prob, nd, xd, yd, zd )

      call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
c
c  Compute the Vandermonde matrix.
c
      call vandermonde_interp_2d_matrix ( nd, m, xd, yd, a )
c
c  Solve linear system.
c
      call qr_solve ( nd, nd, a, zd, c )
c
c  #1:  Does interpolant match function at data points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
        yi(i) = yd(i)
      end do

      call r8poly_value_2d ( m, c, ni, xi, yi, zi )

      app_error = r8vec_norm_affine ( zi - zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 data interpolation error = ', app_error

      return
      end
