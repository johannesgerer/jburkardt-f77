      program main

c*********************************************************************72
c
cc VANDERMONDE_INTERP_1D_TEST tests VANDERMONDE_INTERP_1D.
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

      integer prob
      integer prob_num

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the VANDERMONDE_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
      write ( *, '(a)' ) '  This test needs the TEST_INTERP library.'
      write ( *, '(a)' ) '  This test needs the CONDITION library.'

      call p00_prob_num ( prob_num )
      do prob = 1, prob_num
        call test01 ( prob )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'VANDERMONDE_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob )

c*********************************************************************72
c
cc TEST01 tests VANDERMONDE_INTERP_1D_MATRIX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 October 2012
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

      double precision a(nd_max,nd_max)
      double precision c(nd_max)
      double precision condition
      logical debug
      parameter ( debug = .false. )
      integer i
      double precision int_error
      double precision ld
      double precision li
      integer m
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
      write ( *, '(a,i2)' ) 
     &  '  Interpolate data from TEST_INTERP problem #', prob

      call p00_data_num ( prob, nd )
      write ( *, '(a,i2)' ) '  Number of data points = ', nd

      call p00_data ( prob, 2, nd, xy )
      
      if ( debug ) then
        call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
      end if

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do
c
c  Choose the degree of the polynomial to be ND - 1.
c
      m = nd - 1
c
c  Compute Vandermonde matrix and get condition number.
c
      call vandermonde_interp_1d_matrix ( nd, xd, a )

      call condition_hager ( nd, a, condition )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  Condition of Vandermonde matrix is ', condition
c
c  Solve linear system.
c
      call qr_solve ( nd, nd, a, yd, c )
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do
      call r8poly_value ( m, c, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

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
     &  '  Normalized length of polynomial interpolant       = ', li

      return
      end
