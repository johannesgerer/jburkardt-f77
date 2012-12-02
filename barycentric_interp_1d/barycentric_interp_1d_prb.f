      program main

c*********************************************************************72
c
cc BARYCENTRIC_INTERP_1D_TEST tests BARYCENTRIC_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd_test_num
      parameter ( nd_test_num = 6 )

      integer i
      integer nd
      integer nd_test(nd_test_num)
      integer prob
      integer prob_num

      save nd_test

      data nd_test / 4, 8, 16, 32, 64, 1000 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BARYCENTRIC_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  The tests need the TEST_INTERP_1D libraries.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num
        do i = 1, nd_test_num
          nd = nd_test(i)
          call test01 ( prob, nd )
        end do
      end do

      do prob = 1, prob_num
        do i = 1, nd_test_num
          nd = nd_test(i)
          call test02 ( prob, nd )
        end do
      end do

      do prob = 1, prob_num
        do i = 1, nd_test_num
          nd = nd_test(i)
          call test03 ( prob, nd )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, nd )

c*********************************************************************72
c
cc BARYCENTRIC_INTERP_1D_TEST01 tests LAGCHEBY1_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 September 2012
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

      double precision a
      double precision b
      integer i
      double precision int_error
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd)
      double precision xi(nd)
      double precision yd(nd)
      double precision yi(nd)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST01:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a)' ) 
     &  '  Use Chebyshev Type 1 spacing for data points.'
      write ( *, '(a,i6)' ) '  Number of data points = ', nd
c
c  Define the data.
c
      a =  0.0D+00
      b = +1.0D+00
      call r8vec_cheby1space ( nd, a, b, xd )
      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does the interpolant match the function at the interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do
      call lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
      subroutine test02 ( prob, nd )

c*********************************************************************72
c
cc BARYCENTRIC_INTERP_1D_TEST02 tests LAGCHEBY2_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 September 2012
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

      double precision a
      double precision b
      integer i
      double precision int_error
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd)
      double precision xi(nd)
      double precision yd(nd)
      double precision yi(nd)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST02:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a)' ) 
     &  '  Use Chebyshev Type 2 spacing for data points.'
      write ( *, '(a,i6)' ) '  Number of data points = ', nd
c
c  Define the data.
c
      a =  0.0D+00
      b = +1.0D+00
      call r8vec_cheby2space ( nd, a, b, xd )
      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does the interpolant match the function at the interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do
      call lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
      subroutine test03 ( prob, nd )

c*********************************************************************72
c
cc BARYCENTRIC_INTERP_1D_TEST03 tests LAGEVEN_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 September 2012
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

      double precision a
      double precision b
      integer i
      double precision int_error
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd)
      double precision xi(nd)
      double precision yd(nd)
      double precision yi(nd)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BARYCENTRIC_INTERP_1D_TEST03:'
      write ( *, '(a,i6)' ) 
     &  '  Interpolate data from TEST_INTERP_1D problem #', prob
      write ( *, '(a)' ) '  Use even spacing for data points.'
      write ( *, '(a,i6)' ) '  Number of data points = ', nd
c
c  Define the data.
c
      a =  0.0D+00
      b = +1.0D+00
      call r8vec_midspace ( nd, a, b, xd )
      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does the interpolant match the function at the interpolation points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do
      call lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( ni, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error

      return
      end
