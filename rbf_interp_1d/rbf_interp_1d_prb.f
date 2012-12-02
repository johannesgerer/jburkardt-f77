      program main

c*********************************************************************72
c
cc RBF_INTERP_1D_TEST tests RBF_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd
      external phi1
      external phi2
      external phi3
      external phi4
      integer prob
      integer prob_num
      double precision r0
      double precision, allocatable :: xd(:)
      double precision xmax
      double precision xmin
      double precision, allocatable :: xy(:,:)

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RBF_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the RBF_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test needs the TEST_INTERP library.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num
c
c  Determine an appropriate value of R0, the spacing parameter.
c
        call p00_data_num ( prob, nd )
        allocate ( xy(1:2,1:nd) )
        call p00_data ( prob, 2, nd, xy )
        allocate ( xd(1:nd) )
        xd(1:nd) = xy(1,1:nd)
        xmax = maxval ( xd )
        xmin = minval ( xd )
        r0 = ( xmax - xmin ) / dble ( nd - 1 )
        deallocate ( xd )
        deallocate ( xy )

        call test01 ( prob, phi1, 'phi1', r0 )
        call test01 ( prob, phi2, 'phi2', r0 )
        call test01 ( prob, phi3, 'phi3', r0 )
        call test01 ( prob, phi4, 'phi4', r0 )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'RBF_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, phi, phi_name, r0 )

c*********************************************************************72
c
cc TEST01 tests RBF_INTERP_1D.
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
c  Parameters:
c
c    Input, integer PROB, the index of the problem.
c
c    Input, external PHI, the name of the radial basis function.
c
c    Input, character * ( * ) PHI_NAME, the name of the radial basis function.
c
c    Input, double precision R0, the scale factor.  Typically, this might be
c    a small multiple of the average distance between points.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 49 )
      integer ni_max
      parameter ( ni_max = 501 )

      integer i
      double precision int_error
      double precision ld
      double precision li
      integer m
      integer nd
      integer ni
      external phi
      character * ( * ) phi_name
      integer prob
      double precision r0
      double precision r8vec_norm_affine
      double precision w(nd_max)
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
      write ( *, '(a)' ) 
     &  '  using radial basis function "' // trim ( phi_name ) // '".'
      write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

      call p00_data_num ( prob, nd )
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      call p00_data ( prob, 2, nd, xy )
      
      call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do

      m = 1
      call rbf_weight ( m, nd, xd, r0, phi, yd, w )
c
c  #1:  Does interpolant match function at interpolation points?
c
      ni = nd
      do i = 1, nd
        xi(i) = xd(i)
      end do
      call rbf_interp ( m, nd, xd, r0, phi, w, ni, xi, yi )

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
      call rbf_interp ( m, nd, xd, r0, phi, w, ni, xi, yi )

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
