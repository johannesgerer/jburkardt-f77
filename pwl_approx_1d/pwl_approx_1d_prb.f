      program main

c*********************************************************************72
c
cc PWL_APPROX_1D_TEST tests PWL_APPROX_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nc_test_num
      parameter ( nc_test_num = 4 )
      integer nd_test_num
      parameter ( nd_test_num = 2 )

      integer j
      integer k
      integer nc
      integer nc_test(nc_test_num)
      integer nd
      integer nd_test(nd_test_num)
      integer prob
      integer prob_num

      save nc_test
      save nd_test

      data nc_test / 2, 4, 8, 16 /
      data nd_test / 16, 64 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PWL_APPROX_1D library'
      write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) '  The test needs the TEST_INTERP_1D library.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num
        do j = 1, nc_test_num
          nc = nc_test(j)
          do k = 1, nd_test_num
            nd = nd_test(k)
            call test01 ( prob, nc, nd )
          end do
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, nc, nd )

c*********************************************************************72
c
cc TEST01 tests PWL_APPROX_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer NC, the number of control points.
c
c    Input, integer ND, the number of data points.
c
      implicit none

      integer nc
      integer nd

      double precision app_error
      integer i
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xc(nc)
      double precision xd(nd)
      double precision xi(nd)
      double precision yc(nc)
      double precision yd(nd)
      double precision yi(nd)
      double precision xmax
      double precision xmin

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i4)' )  
     &  '  Approximate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a,i4)' ) '  Number of control points = ', nc
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      xmin = 0.0D+00
      xmax = 1.0D+00

      call r8vec_linspace ( nd, xmin, xmax, xd )

      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  Determine control values.
c
      call r8vec_linspace ( nc, xmin, xmax, xc )

      call pwl_approx_1d ( nd, xd, yd, nc, xc, yc )
c
c  #1:  Does approximant come close to function at data points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do

      call pwl_interp_1d ( nc, xc, yc, ni, xi, yi )

      app_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per data node = ', app_error

      return
      end
