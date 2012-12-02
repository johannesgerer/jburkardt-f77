      program main

c*********************************************************************72
c
cc RBF_INTERP_ND_TEST tests RBF_INTERP_ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the RBF_INTERP_ND library.'
      write ( *, '(a)' ) '  The R8LIB library is also needed.'

      call rbf_interp_nd_test01 ( )
      call rbf_interp_nd_test02 ( )
      call rbf_interp_nd_test03 ( )
      call rbf_interp_nd_test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine rbf_interp_nd_test01 ( )

c*********************************************************************72
c
cc RBF_INTERP_ND_TEST01 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n1d 
      parameter ( n1d = 5 )
      integer ni
      parameter ( ni = 1000 )

      integer nd
      parameter ( nd = n1d**m )

      double precision a
      double precision app_error
      double precision b
      double precision fd(nd)
      double precision fe(ni)
      double precision fi(ni)
      integer i
      double precision int_error
      integer j
      external phi1
      double precision r0
      double precision r8vec_norm_affine
      integer seed
      double precision w(nd)
      double precision x1d(n1d)
      double precision xd(m,nd)
      double precision xi(m,ni)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST01:'
      write ( *, '(a)' ) 
     &  '  RBF_WEIGHT computes weights for RBF interpolation.'
      write ( *, '(a)' ) 
     &  '  RBF_INTERP_ND evaluates the RBF interpolant.'
      write ( *, '(a)' ) 
     &  '  Use the multiquadratic basis function PHI1(R).'

      a = 0.0D+00
      b = 2.0D+00

      call r8vec_linspace ( n1d, a, b, x1d )

      do i = 1, m
        call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
      end do

      call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

      r0 = ( b - a ) / dble ( n1d )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

      do j = 1, nd
        fd(j) = xd(1,j) * xd(2,j) * exp ( - xd(1,j) * xd(2,j) )
      end do

      call r8vec_print ( nd, fd, '  Function data:' )

      call rbf_weight ( m, nd, xd, r0, phi1, fd, w )

      call r8vec_print ( nd, w, '  Weight vector:' )
c
c  #1: Interpolation test.  Does interpolant match function at interpolation points?
c
      call r8mat_copy ( m, nd, xd, xi )

      call rbf_interp_nd ( m, nd, xd, r0, phi1, w, nd, xi, fi )

      int_error = r8vec_norm_affine ( nd, fd, fi ) / dble ( nd )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      seed = 123456789
      call r8mat_uniform_ab ( m, ni, a, b, seed, xi )
      call rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi, fi )

      do j = 1, ni
        fe(j) = xi(1,j) * xi(2,j) * exp ( - xi(1,j) * xi(2,j) )
      end do

      app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) 
     &  / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples = ', 
     &  app_error

      return
      end
      subroutine rbf_interp_nd_test02 ( )

c*********************************************************************72
c
cc RBF_INTERP_ND_TEST02 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n1d 
      parameter ( n1d = 5 )
      integer ni
      parameter ( ni = 1000 )

      integer nd
      parameter ( nd = n1d**m )

      double precision a
      double precision app_error
      double precision b
      double precision fd(nd)
      double precision fe(ni)
      double precision fi(ni)
      integer i
      double precision int_error
      integer j
      external phi2
      double precision r0
      double precision r8vec_norm_affine
      integer seed
      double precision w(nd)
      double precision x1d(n1d)
      double precision xd(m,nd)
      double precision xi(m,ni)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST02:'
      write ( *, '(a)' ) 
     &  '  RBF_WEIGHT computes weights for RBF interpolation.'
      write ( *, '(a)' ) 
     &  '  RBF_INTERP_ND evaluates the RBF interpolant.'
      write ( *, '(a)' ) 
     &  '  Use the inverse multiquadratic basis function PHI2(R).'

      a = 0.0D+00
      b = 2.0D+00

      call r8vec_linspace ( n1d, a, b, x1d )

      do i = 1, m
        call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
      end do

      call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

      r0 = ( b - a ) / dble ( n1d )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

      do j = 1, nd
        fd(j) = xd(1,j) * xd(2,j) * exp ( - xd(1,j) * xd(2,j) )
      end do

      call r8vec_print ( nd, fd, '  Function data:' )

      call rbf_weight ( m, nd, xd, r0, phi2, fd, w )

      call r8vec_print ( nd, w, '  Weight vector:' )
c
c  #1: Interpolation test.  Does interpolant match function at interpolation points?
c
      call r8mat_copy ( m, nd, xd, xi )

      call rbf_interp_nd ( m, nd, xd, r0, phi2, w, nd, xi, fi )

      int_error = r8vec_norm_affine ( nd, fd, fi ) / dble ( nd )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      seed = 123456789
      call r8mat_uniform_ab ( m, ni, a, b, seed, xi )
      call rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi, fi )

      do j = 1, ni
        fe(j) = xi(1,j) * xi(2,j) * exp ( - xi(1,j) * xi(2,j) )
      end do

      app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) 
     &  / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples = ', 
     &  app_error

      return
      end
      subroutine rbf_interp_nd_test03 ( )

c*********************************************************************72
c
cc RBF_INTERP_ND_TEST03 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n1d 
      parameter ( n1d = 5 )
      integer ni
      parameter ( ni = 1000 )

      integer nd
      parameter ( nd = n1d**m )

      double precision a
      double precision app_error
      double precision b
      double precision fd(nd)
      double precision fe(ni)
      double precision fi(ni)
      integer i
      double precision int_error
      integer j
      external phi3
      double precision r0
      double precision r8vec_norm_affine
      integer seed
      double precision w(nd)
      double precision x1d(n1d)
      double precision xd(m,nd)
      double precision xi(m,ni)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST03:'
      write ( *, '(a)' ) 
     &  '  RBF_WEIGHT computes weights for RBF interpolation.'
      write ( *, '(a)' ) 
     &  '  RBF_INTERP_ND evaluates the RBF interpolant.'
      write ( *, '(a)' ) 
     &  '  Use the thin-plate spline basis function PHI3(R).'

      a = 0.0D+00
      b = 2.0D+00

      call r8vec_linspace ( n1d, a, b, x1d )

      do i = 1, m
        call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
      end do

      call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

      r0 = ( b - a ) / dble ( n1d )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

      do j = 1, nd
        fd(j) = xd(1,j) * xd(2,j) * exp ( - xd(1,j) * xd(2,j) )
      end do

      call r8vec_print ( nd, fd, '  Function data:' )

      call rbf_weight ( m, nd, xd, r0, phi3, fd, w )

      call r8vec_print ( nd, w, '  Weight vector:' )
c
c  #1: Interpolation test.  Does interpolant match function at interpolation points?
c
      call r8mat_copy ( m, nd, xd, xi )

      call rbf_interp_nd ( m, nd, xd, r0, phi3, w, nd, xi, fi )

      int_error = r8vec_norm_affine ( nd, fd, fi ) / dble ( nd )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      seed = 123456789
      call r8mat_uniform_ab ( m, ni, a, b, seed, xi )
      call rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi, fi )

      do j = 1, ni
        fe(j) = xi(1,j) * xi(2,j) * exp ( - xi(1,j) * xi(2,j) )
      end do

      app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) 
     &  / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples = ', 
     &  app_error

      return
      end
      subroutine rbf_interp_nd_test04 ( )

c*********************************************************************72
c
cc RBF_INTERP_ND_TEST04 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n1d 
      parameter ( n1d = 5 )
      integer ni
      parameter ( ni = 1000 )

      integer nd
      parameter ( nd = n1d**m )

      double precision a
      double precision app_error
      double precision b
      double precision fd(nd)
      double precision fe(ni)
      double precision fi(ni)
      integer i
      double precision int_error
      integer j
      external phi4
      double precision r0
      double precision r8vec_norm_affine
      integer seed
      double precision w(nd)
      double precision x1d(n1d)
      double precision xd(m,nd)
      double precision xi(m,ni)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RBF_INTERP_ND_TEST04:'
      write ( *, '(a)' ) 
     &  '  RBF_WEIGHT computes weights for RBF interpolation.'
      write ( *, '(a)' ) 
     &  '  RBF_INTERP_ND evaluates the RBF interpolant.'
      write ( *, '(a)' ) 
     &  '  Use the gaussian basis function PHI4(R).'

      a = 0.0D+00
      b = 2.0D+00

      call r8vec_linspace ( n1d, a, b, x1d )

      do i = 1, m
        call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
      end do

      call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

      r0 = ( b - a ) / dble ( n1d )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

      do j = 1, nd
        fd(j) = xd(1,j) * xd(2,j) * exp ( - xd(1,j) * xd(2,j) )
      end do

      call r8vec_print ( nd, fd, '  Function data:' )

      call rbf_weight ( m, nd, xd, r0, phi4, fd, w )

      call r8vec_print ( nd, w, '  Weight vector:' )
c
c  #1: Interpolation test.  Does interpolant match function at interpolation points?
c
      call r8mat_copy ( m, nd, xd, xi )

      call rbf_interp_nd ( m, nd, xd, r0, phi4, w, nd, xi, fi )

      int_error = r8vec_norm_affine ( nd, fd, fi ) / dble ( nd )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 interpolation error averaged per interpolant node = ', 
     &  int_error
c
c  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
c
      seed = 123456789
      call r8mat_uniform_ab ( m, ni, a, b, seed, xi )
      call rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi, fi )

      do j = 1, ni
        fe(j) = xi(1,j) * xi(2,j) * exp ( - xi(1,j) * xi(2,j) )
      end do

      app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) 
     &  / dble ( ni )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L2 approximation error averaged per 1000 samples = ', 
     &  app_error

      return
      end
