      program main

c*********************************************************************72
c
cc PWL_INTERP_2D_TEST tests PWL_INTERP_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_test_num
      parameter ( n_test_num = 5 )

      integer i
      integer n
      integer n_test(n_test_num)
      integer prob
      integer prob_num

      save n_test

      data n_test / 2, 3, 4, 5, 9 /

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PWL_INTERP_2D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) '  The test needs the TEST_INTERP_2D library.'

      call f00_num ( prob_num )
c
c  Numerical tests.
c
      do prob = 1, prob_num
        do i = 1, n_test_num
          n = n_test(i)
          call test01 ( prob, n )
        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_INTERP_2D_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( prob, n )

c*********************************************************************72
c
cc PWL_INTERP_2D_TEST01 tests PWL_INTERP_2D.
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
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer N, the grid size in each dimension.
c
      implicit none

      integer n

      double precision app_error
      integer i
      integer ij
      double precision int_error
      integer j
      integer nd
      integer ni
      integer nxd
      integer nyd
      integer prob
      double precision r8vec_norm_affine
      double precision xd(n*n)
      double precision xd_1d(n)
      double precision xi(n*n)
      double precision xi_1d(n-1)
      double precision yd(n*n)
      double precision yd_1d(n)
      double precision yi(n*n)
      double precision yi_1d(n-1)
      double precision zd(n*n)
      double precision zdm((n-1)*(n-1))
      double precision zi(n*n)

      nxd = n
      nyd = n

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'PWL_INTERP_2D_TEST01:'
      write ( *, '(a,i2)' ) 
     &  '  Interpolate data from TEST_INTERP_2D problem #', prob
      write ( *, '(a,i2,a,i2)' ) 
     &  '  Using polynomial interpolant of product degree ', 
     &  nxd, ' x ', nyd

      nd = nxd * nyd
      write ( *, '(a,i6)' ) '  Number of data points = ', nd

      call r8vec_linspace ( nxd, 0.0D+00, 1.0D+00, xd_1d )
      call r8vec_linspace ( nyd, 0.0D+00, 1.0D+00, yd_1d )

      ij = 0
      do j = 1, nyd
        do i = 1, nxd
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

      call pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi, yi, zi )

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

        do i = 1, n - 1
          xi_1d(i) = 0.5D+00 * ( xd_1d(i) + xd_1d(i+1) )
          yi_1d(i) = 0.5D+00 * ( yd_1d(i) + yd_1d(i+1) )
        end do

        ni = ( nxd - 1 ) * ( nyd - 1 )
        
        ij = 0
        do j = 1, nyd - 1
          do i = 1, nxd - 1
            ij = ij + 1
            xi(ij) = xi_1d(i)
            yi(ij) = yi_1d(j)
          end do
        end do

        call f00_f0 ( prob, ni, xi, yi, zdm )

        call pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi,
     &    yi, zi )

        app_error = r8vec_norm_affine ( ni, zi, zdm ) / dble ( ni )

        write ( *, '(a,g14.6)' ) 
     &    '  RMS data approximation error = ', app_error

      end if

      return
      end
