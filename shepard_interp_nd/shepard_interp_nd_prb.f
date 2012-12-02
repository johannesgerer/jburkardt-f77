      program main

c*********************************************************************72
c
cc SHEPARD_INTERP_ND_TEST tests SHEPARD_INTERP_ND.
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
      implicit none

      integer p_test_num
      parameter ( p_test_num = 4 )

      integer j
      integer m
      integer n1d
      integer nd
      double precision p
      double precision p_test(p_test_num)
      integer prob
      integer prob_num

      save p_test

      data p_test / 1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_ND_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SHEPARD_INTERP_ND library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_ND library.'
c
c  Look at Shepard interpolant on an irregular grid.
c
      nd = 25

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        do m = 2, 5, 3

          do j = 1, p_test_num
            p = p_test(j)
            call test01 ( prob, p, m, nd )
          end do

        end do
      end do
c
c  Look at Shepard interpolant on a regular N1D^M grid.
c
      n1d = 5

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        do m = 2, 5, 3

          do j = 1, p_test_num
            p = p_test(j)
            call test02 ( prob, p, m, n1d )
          end do

        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_ND_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, p, m, nd )

c*********************************************************************72
c
cc TEST01 tests SHEPARD_INTERP on an irregular grid.
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
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, double precision P, the power used in the distance weighting.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer ND, the number of data points.
c
      implicit none

      integer m
      integer nd
      integer ni_max
      parameter ( ni_max = 1000 )

      double precision app_error
      double precision c(m)
      integer i
      double precision int_error
      integer j
      integer ni
      double precision p
      integer prob
      double precision r8vec_norm_affine
      integer seed
      double precision w(m)
      double precision xd(m,nd)
      double precision xi(m,ni_max)
      double precision zd(nd)
      double precision ze(ni_max)
      double precision zi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i4)' ) 
     &  '  Interpolate data from TEST_INTERP_ND problem #', prob
      write ( *, '(a,g14.6)' ) '  using Shepard interpolation with P = ', p
      write ( *, '(a,i4)' ) '  spatial dimension M = ', m
      write ( *, '(a,i4,a)' ) 
     &  '  and an irregular grid of ND = ', nd, ' data points.'
c
c  Set problem parameters:
c
      seed = 123456789
      call r8vec_uniform_01 ( m, seed, c )
      call r8vec_uniform_01 ( m, seed, w )

      call r8mat_uniform_01 ( m, nd, seed, xd )

      call p00_f ( prob, m, c, w, nd, xd, zd )
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do j = 1, ni
        do i = 1, m
          xi(i,j) = xd(i,j)
        end do
      end do

      call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

      int_error = r8vec_norm_affine ( ni, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      ni = 1000
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

      call p00_f ( prob, m, c, w, ni, xi, ze )

      app_error = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples =     ', 
     &  app_error

      return
      end
      subroutine test02 ( prob, p, m, n1d )

c*********************************************************************72
c
cc TEST02 tests SHEPARD_INTERP_ND on a regular N1D^M grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, double precision P, the power used in the distance weighting.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N1D, the number of points in 1D.
c
      implicit none

      integer m
      integer n1d
      integer nd_max
      parameter ( nd_max = 3125 )
      integer ni_max
      parameter ( ni_max = 3125 )

      double precision a
      double precision app_error
      double precision b
      double precision c(m)
      integer i
      double precision int_error
      integer j
      integer nd
      integer ni
      double precision p
      integer prob
      double precision r8vec_norm_affine
      integer seed
      double precision w(m)
      double precision x1d(n1d)
      double precision xd(m,nd_max)
      double precision xi(m,ni_max)
      double precision zd(nd_max)
      double precision ze(ni_max)
      double precision zi(ni_max)
c
c  Set problem parameters:
c
      seed = 123456789
      call r8vec_uniform_01 ( m, seed, c )
      call r8vec_uniform_01 ( m, seed, w )

      nd = n1d ** m

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a,i4)' ) 
     &  '  Interpolate data from TEST_INTERP_ND problem #', prob
      write ( *, '(a,g14.6)' ) 
     &  '  using Shepard interpolation with P = ', p
      write ( *, '(a,i4)' ) '  spatial dimension M = ', m
      write ( *, '(a,i6,a)' ) 
     &  '  and a regular grid of N1D^M = ', nd, ' data points.'

      a = 0.0D+00
      b = 1.0D+00

      call r8vec_linspace ( n1d, a, b, x1d )

      do i = 1, m
        call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
      end do

      call p00_f ( prob, m, c, w, nd, xd, zd )
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do j = 1, nd
        do i = 1, m
          xi(i,j) = xd(i,j)
        end do
      end do

      call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

      int_error = r8vec_norm_affine ( ni, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      ni = 1000
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

      call p00_f ( prob, m, c, w, ni, xi, ze )

      app_error = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples =     ', 
     &  app_error

      return
      end
