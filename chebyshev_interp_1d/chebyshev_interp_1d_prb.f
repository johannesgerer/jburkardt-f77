      program main

!*********************************************************************72
!
!! CHEBYSHEV_INTERP_1D_TEST tests CHEBYSHEV_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer prob
      integer prob_num

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CHEBYSHEV_INTERP_1D library.'
      write ( *, '(a)' ) 
     &  '  The QR_SOLVE and R8LIB libraries are needed.'
      write ( *, '(a)' ) 
     &  '  The test needs the TEST_INTERP library as well.'

      call p00_prob_num ( prob_num )
      do prob = 1, prob_num
        call test01 ( prob )
      end do
!
!  Terminate.
!
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob )

!*********************************************************************72
!
!! TEST01 tests CHEBYSHEV_VALUE_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2012
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer nd_max
      parameter ( nd_max = 49 )

      integer i
      double precision int_error
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd_max)
      double precision xi(nd_max)
      double precision xy(2,nd_max)
      double precision yd(nd_max)
      double precision yi(nd_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST01:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP problem #', prob

      call p00_data_num ( prob, nd )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call p00_data ( prob, 2, nd, xy )
      
      call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

      do i = 1, nd
        xd(i) = xy(1,i)
        yd(i) = xy(2,i)
      end do
!
!  #1:  Does interpolant match function at interpolation points?
!
      ni = nd

      do i = 1, ni
        xi(i) = xd(i)
      end do

      call chebyshev_interp_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
