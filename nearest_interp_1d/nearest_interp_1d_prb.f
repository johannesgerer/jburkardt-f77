      program main

c*********************************************************************72
c
cc NEAREST_INTERP_1D_TEST tests NEAREST_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ni
      integer prob
      integer prob_num

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the NEAREST_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  The test needs the TEST_INTERP library.'

      call p00_prob_num ( prob_num )

      ni = 11
      do prob = 1, prob_num
        call nearest_interp_1d_test01 ( prob, ni )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine nearest_interp_1d_test01 ( prob, ni )

c*********************************************************************72
c
cc NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the problem.
c
c    Input, integer NI, the number of interpolation points.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 49 )

      double precision d(2,nd_max)
      integer j
      integer ni
      integer nd
      integer prob
      character * ( 80 ) title
      double precision xd(nd_max)
      double precision xi(ni)
      double precision xd_max
      double precision xd_min
      double precision yd(nd_max)
      double precision yi(ni)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST01'
      write ( *, '(a,i2)' ) 
     &  '  Sample the nearest neighbor interpolant for problem # ', prob

      call p00_data_num ( prob, nd )

      call p00_data ( prob, 2, nd, d )

      do j = 1, nd
        xd(j) = d(1,j)
        yd(j) = d(2,j)
      end do

      call r8vec_min ( nd, xd, xd_min )
      call r8vec_max ( nd, xd, xd_max )

      call r8vec_linspace ( ni, xd_min, xd_max, xi )
      call nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

      write ( title, '(a,i2)' ) 'X, Y for problem ', prob

      call r8vec2_print ( ni, xi, yi, title )

      return
      end

