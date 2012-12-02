      program main

c*********************************************************************72
c
cc LAGRANGE_APPROX_1D_TEST tests LAGRANGE_APPROX_1D.
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
      parameter ( m_test_num = 7 )
      integer nd_test_num
      parameter ( nd_test_num = 3 )

      integer j
      integer k
      integer m
      integer m_test(m_test_num)
      integer nd
      integer nd_test(nd_test_num)
      integer prob
      integer prob_num

      save m_test
      save nd_test

      data m_test / 0, 1, 2, 3, 4, 8, 16 /
      data nd_test / 16, 64, 1000 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAGRANGE_APPROX_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
      write ( *, '(a)' ) 
     &  '  These tests need the TEST_INTERP_1D library.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num
        do j = 1, m_test_num
          m = m_test(j)
          do k = 1, nd_test_num
            nd = nd_test(k)
            call test02 ( prob, m, nd )
          end do
        end do
      end do

      do prob = 1, prob_num
        do j = 1, m_test_num
          m = m_test(j)
          do k = 1, nd_test_num
            nd = nd_test(k)
            call test03 ( prob, m, nd )
          end do
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_APPROX_1D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test02 ( prob, m, nd )

c*********************************************************************72
c
cc TEST02 tests LAGRANGE_APPROX_1D with evenly spaced data
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer M, the polynomial approximant degree.
c
c    Input, integer ND, the number of data points.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 1000 )
      integer ni_max
      parameter ( ni_max = 1000 )

      double precision a
      double precision b
      integer i
      double precision int_error
      integer m
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision yd(nd_max)
      double precision yi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a,i4)' ) 
     &  '  Approximate evenly spaced data from problem #', prob
      write ( *, '(a,i4)' ) '  Use polynomial approximant of degree ', m
      write ( *, '(a,i4)' ) '  Number of data points = ', nd

      a = 0.0D+00
      b = 1.0D+00
      call r8vec_linspace ( nd, a, b, xd )

      call p00_f ( prob, nd, xd, yd )

      if ( nd .lt. 10 ) then
        call r8vec2_print ( nd, xd, yd, '  Data array:' )
      end if
c
c  #1:  Does approximant come close to function at data points?
c
      ni = nd
      do i = 1, ni
        xi(i) = xd(i)
      end do
      call lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( nd, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per data node = ', int_error

      return
      end
      subroutine test03 ( prob, m, nd )

c*********************************************************************72
c
cc  TEST03 tests LAGRANGE_APPROX_1D with Chebyshev spaced data.
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
c    Input, integer PROB, the problem index.
c
c    Input, integer M, the polynomial approximant degree.
c
c    Input, integer ND, the number of data points.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 1000 )
      integer ni_max
      parameter ( ni_max = 1000 )

      double precision a
      double precision b
      integer i
      double precision int_error
      integer m
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd(nd_max)
      double precision xi(ni_max)
      double precision yd(nd_max)
      double precision yi(ni_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a,i4)' ) 
     & '  Approximate Chebyshev-spaced data from problem #', prob
      write ( *, '(a,i4)' ) '  Use polynomial approximant of degree ', m
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
      call lagrange_approx_1d ( m, nd, xd, yd, ni, xi, yi )

      int_error = r8vec_norm_affine ( nd, yi, yd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per data node = ', int_error


      return
      end
