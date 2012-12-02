      program main

c*********************************************************************72
c
cc LAGRANGE_INTERP_1D_TEST tests LAGRANGE_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd_test_num
      parameter ( nd_test_num = 6 )

      integer j
      integer nd
      integer nd_test(nd_test_num)
      integer prob
      integer prob_num

      save nd_test

      data nd_test / 4, 8, 16, 32, 64, 256 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  These tests need the TEST_INTERP_1D library.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num
        do j = 1, nd_test_num
          nd = nd_test(j)
          call test02 ( prob, nd )
        end do
      end do

      do prob = 1, prob_num
        do j = 1, nd_test_num
          nd = nd_test(j)
          call test03 ( prob, nd )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test02 ( prob, nd )

c*********************************************************************72
c
cc TEST02 tests LAGRANGE_VALUE_1D with evenly spaced data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer ND, the number of data points to use.
c
      implicit none

      integer nd
      integer ni_max
      parameter ( ni_max = 501 )

      double precision a
      double precision b
      integer i
      double precision int_error
      double precision ld
      double precision li
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd)
      double precision xi(ni_max)
      double precision yd(nd)
      double precision yi(ni_max)
      double precision ymax
      double precision ymin

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a,i4)' ) 
     &  '  Interpolate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a)' ) '  Use even spacing for data points.'
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      a = 0.0D+00
      b = 1.0D+00
      
      call r8vec_linspace ( nd, a, b, xd )

      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do

      call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( nd, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error 
c
c  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
c  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
c  (YMAX-YMIN).
c
      call r8vec_min ( nd, yd, ymin )
      call r8vec_max ( nd, yd, ymax )

      ni = 501
      call r8vec_linspace ( ni, a, b, xi )
      call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

      ld = 0.0D+00
      do i = 1, nd - 1
        ld = ld + sqrt ( ( ( xd(i+1) - xd(i) ) / ( b - a ) )**2 
     &                 + ( ( yd(i+1) - yd(i) ) / ( ymax - ymin ) )**2 ) 
      end do

      li = 0.0D+00
      do i = 1, ni - 1
        li = li + sqrt ( ( ( xi(i+1) - xi(i) ) / ( b - a ) )**2 
     &                 + ( ( yi(i+1) - yi(i) ) / ( ymax - ymin ) )**2 )
      end do

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of piecewise linear interpolant = ', ld
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of polynomial interpolant       = ', li

      return
      end
      subroutine test03 ( prob, nd )

c*********************************************************************72
c
cc TEST03 tests LAGRANGE_VALUE_1D with Chebyshev spaced data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer ND, the number of data points to use.
c
      implicit none

      integer nd
      integer ni_max
      parameter ( ni_max = 501 )

      double precision a
      double precision b
      integer i
      double precision int_error
      double precision ld
      double precision li
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd)
      double precision xi(ni_max)
      double precision yd(nd)
      double precision yi(ni_max)
      double precision ymax
      double precision ymin

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a,i4)' ) 
     &  '  Interpolate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a)' ) '  Use Chebyshev spacing for data points.'
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      a = 0.0D+00
      b = 1.0D+00
      call r8vec_chebyspace ( nd, a, b, xd )

      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do

      call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( nd, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error 
c
c  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
c  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
c  (YMAX-YMIN).
c
      call r8vec_min ( nd, yd, ymin )
      call r8vec_max ( nd, yd, ymax )

      ni = 501
      call r8vec_linspace ( ni, a, b, xi )
      call lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

      ld = 0.0D+00
      do i = 1, nd - 1
        ld = ld + sqrt ( ( ( xd(i+1) - xd(i) ) / ( b - a ) )**2 
     &                 + ( ( yd(i+1) - yd(i) ) / ( ymax - ymin ) )**2 ) 
      end do

      li = 0.0D+00
      do i = 1, ni - 1
        li = li + sqrt ( ( ( xi(i+1) - xi(i) ) / ( b - a ) )**2 
     &                 + ( ( yi(i+1) - yi(i) ) / ( ymax - ymin ) )**2 )
      end do

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of piecewise linear interpolant = ', ld
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of polynomial interpolant       = ', li

      return
      end
