      program main

c*********************************************************************72
c
cc VANDERMONDE_APPROX_1D_TEST tests VANDERMONDE_APPROX_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_test_num
      parameter ( m_test_num = 8 )

      integer j
      integer m
      integer m_test(m_test_num)
      integer prob
      integer prob_num

      save m_test

      data m_test / 0, 1, 2, 3, 4, 5, 9, 12 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the VANDERMONDE_APPROX_1D library.'
      write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test needs the TEST_INTERP libary as well.'

      call p00_prob_num ( prob_num )
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
      write ( *, '(a)' ) 'VANDERMONDE_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, m )

c*********************************************************************72
c
cc TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer M, the polynomial degree.
c
      implicit none

      integer m
      integer nd_max
      parameter ( nd_max = 49 )
      integer ni_max
      parameter ( ni_max = 501 )

      double precision a(nd_max,0:m)
      double precision app_error
      double precision c(0:m)
      integer i
      double precision ld
      double precision li
      integer nd
      integer ni
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
     &  '  Approximate data from TEST_INTERP problem #', prob

      call p00_data_num ( prob, nd )
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      call p00_data ( prob, 2, nd, xy )
      
      if ( m .eq. 0 ) then
        call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
      end if

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do
c
c  Compute the Vandermonde matrix.
c
      write ( *, '(a,i4)' ) 
     &  '  Using polynomial approximant of degree ', m

      call vandermonde_approx_1d_matrix ( nd, m, xd, a )
c
c  Solve linear system.
c
      call qr_solve ( nd, m + 1, a, yd, c )
c
c  #1:  Does approximant match function at data points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do

      call r8poly_value ( m, c, ni, xi, yi )

      app_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 data approximation error = ', app_error
c
c  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
c  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
c  (YMAX-YMIN).
c
      call r8vec_max ( nd, xd, xmax )
      call r8vec_min ( nd, xd, xmin )
      call r8vec_max ( nd, yd, ymax )
      call r8vec_min ( nd, yd, ymin )

      ni = 501
      call r8vec_linspace ( ni, xmin, xmax, xi )

      call r8poly_value ( m, c, ni, xi, yi )
  
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
     &  '  Normalized length of polynomial approximant       = ', li

      return
      end
