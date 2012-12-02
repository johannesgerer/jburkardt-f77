      program main

c*********************************************************************72
c
cc LAGRANGE_INTERP_2D_TEST tests LAGRANGE_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_test_num
      parameter ( m_test_num = 5 )

      integer i
      integer m
      integer m_test(m_test_num)
      integer prob
      integer prob_num

      save m_test

      data m_test / 1, 2, 3, 4, 8 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_2D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test also needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )
c
c  Numerical tests.
c
      do prob = 1, prob_num
        do i = 1, m_test_num
          m = m_test(i)
          call test01 ( prob, m )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, m )

c*****************************************************************************80
c
cc LAGRANGE_INTERP_2D_TEST01 tests LAGRANGE_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer M, the polynomial degree in each dimension.
c
      implicit none

      integer m_max
      parameter ( m_max = 8 )

      double precision app_error
      integer i
      integer ij
      double precision int_error
      integer j
      integer m
      integer mx
      integer my
      integer nd
      integer ni
      integer prob
      double precision r8vec_norm_affine
      double precision xd((m_max+1)*(m_max+1))
      double precision xd_1d(m_max+1)
      double precision xi((m_max+1)*(m_max+1))
      double precision xi_1d(m_max)
      double precision yd((m_max+1)*(m_max+1))
      double precision yd_1d(m_max+1)
      double precision yi((m_max+1)*(m_max+1))
      double precision yi_1d(m_max)
      double precision zd((m_max+1)*(m_max+1))
      double precision zdm(m_max*m_max)
      double precision zi((m_max+1)*(m_max+1))

      mx = m
      my = m

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'LAGRANGE_INTERP_2D_TEST01:'
      write ( *, '(a,i2)' ) 
     &  '  Interpolate data from TEST_INTERP_2D problem #', prob
      write ( *, '(a,i2,a,i2)' ) 
     &  '  Using polynomial interpolant of product degree ', 
     &  mx, ' x ', my

      nd = ( mx + 1 ) * ( my + 1 )
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call r8vec_chebyspace ( mx + 1, 0.0D+00, 1.0D+00, xd_1d )
      call r8vec_chebyspace ( my + 1, 0.0D+00, 1.0D+00, yd_1d )

      ij = 0
      do j = 1, my + 1
        do i = 1, mx + 1
          ij = ij + 1
          xd(ij) = xd_1d(i)
          yd(ij) = yd_1d(j)
        end do
      end do

      call f00_f0 ( prob, nd, xd, yd, zd )

      if ( nd .le. 20 ) then
        call r8vec3_print ( nd, xd, yd, zd, '  X, Y, Z data:' )
      end if
c
c  #1:  Does interpolant match function at data points?
c
      ni = nd

      do i = 1, ni
        xi(i) = xd(i)
        yi(i) = yd(i)
      end do

      call lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, 
     &  yi, zi )

      if ( ni .le. 20 ) then
        call r8vec3_print ( ni, xi, yi, zi, '  X, Y, Z interpolation:' )
      end if

      int_error = r8vec_norm_affine ( ni, zi, zd ) / dble ( ni )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  RMS data interpolation error = ', int_error
c
c  #2:  Does interpolant approximate data at midpoints?
c
      if ( 1 .lt. nd ) then

        do i = 1, mx
          xi_1d(i) = 0.5D+00 * ( xd_1d(i) + xd_1d(i+1) )
        end do

        do i = 1, my
          yi_1d(i) = 0.5D+00 * ( yd_1d(i) + yd_1d(i+1) )
        end do

        ni = mx * my
       
        ij = 0
        do j = 1, my
          do i = 1, mx
            ij = ij + 1
            xi(ij) = xi_1d(i)
            yi(ij) = yi_1d(j)
          end do
        end do

        call f00_f0 ( prob, ni, xi, yi, zdm )

        call lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, 
     &    yi, zi )

        app_error = r8vec_norm_affine ( ni, zi, zdm ) / dble ( ni )

        write ( *, '(a)' ) ''
        write ( *, '(a,g14.6)' ) '  RMS data approximation error = ', 
     &    app_error

      end if

      return
      end
