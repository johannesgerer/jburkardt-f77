      program main

c*********************************************************************72
c
cc SHEPARD_INTERP_1D_TEST tests SHEPARD_INTERP_1D.
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
      implicit none

      integer p_num
      parameter ( p_num = 5 )

      integer j
      double precision p
      double precision p_test(p_num)
      integer prob
      integer prob_num

      save p_test

      data p_test / 0.0D+00, 1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SHEPARD_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test needs the TEST_INTERP library as well.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        do j = 1, p_num
          p = p_test(j)
          call test01 ( prob, p )
        end do

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SHEPARD_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, p )

c*********************************************************************72
c
cc TEST01 tests SHEPARD_INTERP_1D.
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
      implicit none

      integer nd_max
      parameter ( nd_max = 49 )
      integer ni_max
      parameter ( ni_max = 501 )

      integer dim_num
      integer i
      double precision int_error
      double precision ld
      double precision li
      integer nd
      integer ni
      double precision p
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision xmax
      double precision xmin
      double precision xy(2,nd_max)
      double precision yd(nd_max)
      double precision yi(ni_max)
      double precision ymax
      double precision ymin

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i4)' ) 
     &  '  Interpolate data from TEST_INTERP problem #', prob
      write ( *, '(a,g14.6)' ) 
     &  '  using Shepard interpolation with P = ', p

      call p00_dim_num ( prob, dim_num )

      call p00_data_num ( prob, nd )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call p00_data ( prob, dim_num, nd, xy )
      
      if ( p .eq. 0.0D+00 ) then
        call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
      end if

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do

      call shepard_interp_1d ( nd, xd, yd, p, ni, xi, yi )

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
      call r8vec_min ( nd, xd, xmin )
      call r8vec_max ( nd, xd, xmax )
      call r8vec_min ( nd, yd, ymin )
      call r8vec_max ( nd, yd, ymax )

      ni = 501
      call r8vec_linspace ( ni, xmin, xmax, xi )
      call shepard_interp_1d ( nd, xd, yd, p, ni, xi, yi )

      ld = 0.0D+00
      do i = 1, nd - 1
        ld = ld + sqrt 
     &    ( ( ( xd(i+1) - xd(i) ) / ( xmax - xmin ) ) ** 2 
     &    + ( ( yd(i+1) - yd(i) ) / ( ymax - ymin ) ) ** 2 ) 
      end do

      li = 0.0D+00
      do i = 1, ni - 1
        li = li + sqrt 
     &    ( ( ( xi(i+1) - xi(i) ) / ( xmax - xmin ) ) ** 2 
     &    + ( ( yi(i+1) - yi(i) ) / ( ymax - ymin ) ) ** 2 ) 
      end do

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of piecewise linear interpolant = ', ld
      write ( *, '(a,g14.6)' ) 
     &  '  Normalized length of Shepard interpolant          = ', li

      return
      end
