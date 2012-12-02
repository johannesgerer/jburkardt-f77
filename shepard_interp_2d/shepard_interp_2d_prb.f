      program main

c*********************************************************************72
c
cc SHEPARD_INTERP_2D_PRB tests SHEPARD_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer p_test_num
      parameter ( p_test_num = 4 )

      integer g
      integer j
      double precision p
      double precision p_test(p_test_num)
      integer prob
      integer prob_num

      save p_test

      data p_test / 1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_2D_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SHEPARD_INTERP_2D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )
      g = 1

      do prob = 1, prob_num

        do j = 1, p_test_num
          p = p_test(j)
          call test01 ( prob, g, p )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_2D_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, g, p )

c*****************************************************************************80
c
cc TEST01 tests SHEPARD_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer G, the grid number.
c
c    Input, double precision P, the power used in the distance weighting.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 100 )
      integer ni_max
      parameter ( ni_max = 100 )

      logical debug
      parameter ( debug = .false. )
      integer g
      integer i
      double precision int_error
      integer nd
      integer ni
      double precision p
      integer prob
      double precision r8vec_norm_affine
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
      write ( *, '(a,i6)' ) '  using grid #', g
      write ( *, '(a,g14.6)' ) 
     &  '  using Shepard interpolation with P = ', p

      call g00_size ( g, nd )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call g00_xy ( g, nd, xd, yd )
      

      call f00_f0 ( prob, nd, xd, yd, zd )

      if ( debug ) then
        call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
      end if
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
        yi(i) = yd(i)
      end do

      call shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi )

      int_error = r8vec_norm_affine ( nd, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
