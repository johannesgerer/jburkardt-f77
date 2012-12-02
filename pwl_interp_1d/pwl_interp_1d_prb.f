      program main

c*********************************************************************72
c
cc PWL_INTERP_1D_TEST tests PWL_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 September 2012
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
      write ( *, '(a)' ) 'PWL_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PWL_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) '  The test needs the TEST_INTERP library.'

      call p00_prob_num ( prob_num )
      do prob = 1, prob_num
        call pwl_interp_1d_test01 ( prob )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine pwl_interp_1d_test01 ( prob )

c*********************************************************************72
c
cc PWL_INTERP_1D_TEST01 tests PWL_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 49 )
      integer ni_max
      parameter ( ni_max = nd_max )

      integer i
      double precision int_error
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision xy(2,nd_max)
      double precision yd(nd_max)
      double precision yi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_INTERP_1D_TEST01:'
      write ( *, '(a)' ) 
     &  '  PWL_INTERP_1D evaluates the piecewise linear interpolant.'
      write ( *, '(a,i2)' ) 
     &  '  Interpolate data from TEST_INTERP problem #', prob

      call p00_data_num ( prob, nd )
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      call p00_data ( prob, 2, nd, xy )
      
      call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do
c
c  Does interpolant match function at interpolation points?
c
      ni = nd

      do i = 1, ni
        xi(i) = xd(i)
      end do

      call pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
