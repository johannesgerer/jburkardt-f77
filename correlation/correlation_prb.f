      program main

c*********************************************************************72
c
cc CORRELATION_PRB tests the CORRELATION library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the CORRELATION library.'

      call correlation_test01 ( )
      call correlation_test02 ( )
      call correlation_test03 ( )
      call correlation_test04 ( )
      call correlation_test05 ( )
      call correlation_test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine correlation_test01 ( )

c*********************************************************************72
c
cc CORRELATION_TEST01 plots the correlation functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )

      double precision c(n)
      double precision e
      double precision nu
      double precision rho(n)
      double precision rho0

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TEST01'
      write ( *, '(a)' ) '  Make plots of correlation functions.'
      write ( *, '(a)' ) ' '
c
c  besselj
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -8.0D+00, 8.0D+00, rho )
      call correlation_besselj ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'besselj', 
     &  'Bessel J correlation' )
c
c  besselk
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -4.0D+00, 4.0D+00, rho )
      call correlation_besselk ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'besselk', 
     &  'Bessel K correlation' )
c
c  circular
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_circular ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'circular', 
     &  'Circular correlation' )
c
c  constant
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_constant ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'constant', 
     &  'Constant correlation' )
c
c  cubic
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_cubic ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'cubic', 
     &  'Cubic correlation' )
c
c  damped_cosine
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -6.0D+00, 6.0D+00, rho )
      call correlation_damped_cosine ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'damped_cosine', 
     &  'Damped cosine correlation' )
c
c  damped_sine
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -12.0D+00, 12.0D+00, rho )
      call correlation_damped_sine ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'damped_sine', 
     &  'Damped sine correlation' )
c
c  exponential
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_exponential ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'exponential', 
     &  'Exponential correlation' )
c
c  gaussian
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_gaussian ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'gaussian', 
     &  'Gaussian correlation' )
c
c  hole
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -6.0D+00, 6.0D+00, rho )
      call correlation_hole ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'hole', 
     &  'Hole correlation' )
c
c  linear
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_linear ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'linear', 
     &  'Linear correlation' )
c
c  matern, nu = 2.5
c
      rho0 = 1.0D+00
      nu = 2.5D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_matern ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'matern', 
     &  'Matern correlation (NU = 2.5)' )
c
c  pentaspherical
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_pentaspherical ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'pentaspherical', 
     &  'Pentaspherical correlation' )
c
c  power, e = 2.0
c
      rho0 = 1.0D+00
      e = 2.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_power ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'power', 
     &  'Power correlation' )
c
c  rational_quadratic
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -4.0D+00, 4.0D+00, rho )
      call correlation_rational_quadratic ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'rational_quadratic', 
     &  'Rational quadratic correlation' )
c
c  spherical
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_spherical ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'spherical', 
     &  'Spherical correlation' )
c
c  white_noise
c
      rho0 = 1.0D+00
      call r8vec_linspace ( n, -2.0D+00, 2.0D+00, rho )
      call correlation_white_noise ( n, rho, rho0, c )
      call correlation_plot ( n, rho, c, 'white_noise', 
     &  'White noise correlation' )

      return
      end
      subroutine correlation_test02 ( )

c*********************************************************************72
c
cc CORRELATION_TEST02 plots sample paths with SAMPLE_PATHS_CHOLESKY.
c
c  Discussion:
c
c    Most paths will be blue, but make the LAST one red so that there will
c    always be one distinguished path that is easy to follow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )
      integer n2
      parameter ( n2 = 3 )

      external correlation_besselj
      external correlation_besselk
      external correlation_circular
      external correlation_constant
      external correlation_cubic
      external correlation_damped_cosine
      external correlation_damped_sine
      external correlation_exponential
      external correlation_gaussian
      external correlation_hole
      external correlation_linear
      external correlation_matern
      external correlation_pentaspherical
      external correlation_power
      external correlation_rational_quadratic
      external correlation_spherical
      external correlation_white_noise
      double precision rho(n)
      double precision rho0
      parameter ( rho0 = 1.0D+00 )
      double precision rhomax
      parameter ( rhomax = 10.0D+00 )
      double precision rhomin
      parameter ( rhomin = 0.0D+00 )
      integer seed
      double precision x(n,n2)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TEST02'
      write ( *, '(a)' ) 
     &  '  SAMPLE_PATHS_CHOLESKY generates sample paths from the'
      write ( *, '(a)' ) 
     &  '  correlation matrix, factored using the Cholesky factor.'
      write ( *, '(a)' ) 
     &  '  It requires that the correlation matrix is nonnegative'
      write ( *, '(a)' ) '  definite.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  SAMPLE_PATHS_EIGEN generates sample paths from the'
      write ( *, '(a)' ) 
     &  '  correlation matrix, factored using the eigen factorization.'
      write ( *, '(a)' ) 
     &  '  If the correlation matrix is not nonnegative definite,'
      write ( *, '(a)' ) 
     &  '  we simply suppress negative eigenvalues.'
      write ( *, '(a)' ) ''

      call r8vec_linspace ( n, rhomin, rhomax, rho )
c
c  besselj
c  Use EIGEN, because CHOLESKY fails.
c
      seed = 123456789
      call sample_paths_eigen ( n, n2, rhomax, rho0, 
     &  correlation_besselj, seed, x )
      call paths_plot ( n, n2, rho, x, 'besselj', 
     &  'Bessel J correlation' )
c
c  besselk
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_besselk, seed, x )
      call paths_plot ( n, n2, rho, x, 'besselk', 
     &  'Bessel K correlation' )
c
c  circular
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_circular, seed, x )
      call paths_plot ( n, n2, rho, x, 'circular', 
     &  'Circular correlation' )
c
c  constant
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_constant, seed, x )
      call paths_plot ( n, n2, rho, x, 'constant', 
     &  'Constant correlation' )
c
c  cubic
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_cubic, seed, x )
      call paths_plot ( n, n2, rho, x, 'cubic', 
     &  'Cubic correlation' )
c
c  damped_cosine
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_damped_cosine, seed, x )
      call paths_plot ( n, n2, rho, x, 'damped_cosine', 
     &  'Damped cosine correlation' )
c
c  damped_sine
c  Use EIGEN, because CHOLESKY fails.
c
      seed = 123456789
      call sample_paths_eigen ( n, n2, rhomax, rho0, 
     &  correlation_damped_sine, seed, x )
      call paths_plot ( n, n2, rho, x, 'damped_sine', 
     &  'Damped sine correlation' )
c
c  exponential
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_exponential, seed, x )
      call paths_plot ( n, n2, rho, x, 'exponential', 
     &  'Exponential correlation' )
c
c  gaussian
c  Use EIGEN, because CHOLESKY fails.
c
      seed = 123456789
      call sample_paths_eigen ( n, n2, rhomax, rho0, 
     &  correlation_gaussian, seed, x )
      call paths_plot ( n, n2, rho, x, 'gaussian', 
     &  'Gaussian correlation' )
c
c  hole
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_hole, seed, x )
      call paths_plot ( n, n2, rho, x, 'hole', 
     &  'Hole correlation' )
c
c  linear
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_linear, seed, x )
      call paths_plot ( n, n2, rho, x, 'linear', 
     &  'Linear correlation' )
c
c  matern ( nu = 2.5 )
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_matern, seed, x )
      call paths_plot ( n, n2, rho, x, 'matern', 
     &  'Matern correlation (nu=2.5)' )
c
c  pentaspherical
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_pentaspherical, seed, x )
      call paths_plot ( n, n2, rho, x, 'pentaspherical', 
     &  'Pentaspherical correlation' )
c
c  power ( e = 2.0 )
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_power, seed, x )
      call paths_plot ( n, n2, rho, x, 'power', 
     &  'Power correlation (e=2.0)' )
c
c  rational_quadratic
c  Use EIGEN, because CHOLESKY fails.
c
      seed = 123456789
      call sample_paths_eigen ( n, n2, rhomax, rho0, 
     &  correlation_rational_quadratic, seed, x )
      call paths_plot ( n, n2, rho, x, 'rational_quadratic', 
     &  'Rational quadratic correlation' )
c
c  spherical
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_spherical, seed, x )
      call paths_plot ( n, n2, rho, x, 'spherical', 
     &  'Spherical correlation' )
c
c  white_noise
c
      seed = 123456789
      call sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation_white_noise, seed, x )
      call paths_plot ( n, n2, rho, x, 'white_noise', 
     &  'White noise correlation' )

      return
      end
      subroutine correlation_test03 ( )

c*********************************************************************72
c
cc CORRELATION_TEST03 plots a correlation function for several values of RH00.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )
      integer n2_max
      parameter ( n2_max = 7 )

      double precision c(n,n2_max)
      integer j
      integer n2
      double precision rho(n)
      double precision rho0(n2_max)
      double precision rhomax
      double precision rhomin

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TEST03'
      write ( *, '(a)' ) '  Make plots of correlation functions for'
      write ( *, '(a)' ) '  a range of correlation lengths.'
      write ( *, '(a)' ) ''
c
c  besselj
c
      n2 = 5
      rho0(1) = 1.0
      rho0(2) = 1.5
      rho0(3) = 2.0
      rho0(4) = 4.0
      rho0(5) = 8.0
      rhomin = -8.0D+00
      rhomax = +8.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_besselj ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'besselj', 
     &  'Bessel J correlation' )
c
c  besselk
c
      n2 = 5
      rho0(1) = 1.0
      rho0(2) = 1.5
      rho0(3) = 2.0
      rho0(4) = 4.0
      rho0(5) = 8.0
      rhomin = -4.0D+00
      rhomax = +4.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_besselk ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'besselk', 
     &  'Bessel K correlation' )
c
c  circular
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_circular ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'circular', 
     &  'Circular correlation' )
c
c  constant
c  1 plot is enough
c
      n2 = 1
      rho0(1) = 1.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_constant ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'constant', 
     &  'Constant correlation' )
c
c  cubic
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -8.0D+00
      rhomax = +8.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_cubic ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'cubic', 
     &  'Cubic correlation' )
c
c  damped_cosine
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -6.0D+00
      rhomax = +6.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_damped_cosine ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'damped_cosine', 
     &  'Damped cosine correlation' )
c
c  damped_sine
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -8.0D+00
      rhomax = +8.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_damped_sine ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'damped_sine', 
     &  'Damped sine correlation' )
c
c  exponential
c
      n2 = 7
      rho0(1) = 0.25
      rho0(2) = 0.5
      rho0(3) = 1.0
      rho0(4) = 1.5
      rho0(5) = 2.0
      rho0(6) = 4.0
      rho0(7) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_exponential ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'exponential', 
     &  'Exponential correlation' )
c
c  gaussian
c
      n2 = 7
      rho0(1) = 0.25
      rho0(2) = 0.5
      rho0(3) = 1.0
      rho0(4) = 1.5
      rho0(5) = 2.0
      rho0(6) = 4.0
      rho0(7) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_gaussian ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'gaussian', 
     &  'Gaussian correlation' )
c
c  hole
c
      n2 = 7
      rho0(1) = 0.25
      rho0(2) = 0.5
      rho0(3) = 1.0
      rho0(4) = 1.5
      rho0(5) = 2.0
      rho0(6) = 4.0
      rho0(7) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_hole ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'hole', 
     &  'Hole correlation' )
c
c  linear
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_linear ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'linear', 
     &  'Linear correlation' )
c
c  matern, nu = 2.5
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_matern ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'matern', 
     &  'Matern correlation (NU = 2.5)' )
c
c  pentaspherical
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_pentaspherical ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'pentaspherical', 
     &  'Pentaspherical correlation' )
c
c  power, e = 2.0
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_power ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'power', 
     &  'Power correlation (E = 2.0)' )
c
c  rational_quadratic
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -4.0D+00
      rhomax = +4.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_rational_quadratic ( n, rho, rho0(j), 
     &    c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 
     &  'rational_quadratic', 'Rational quadratic correlation' )
c
c  spherical
c
      n2 = 6
      rho0(1) = 0.5
      rho0(2) = 1.0
      rho0(3) = 1.5
      rho0(4) = 2.0
      rho0(5) = 4.0
      rho0(6) = 8.0
      rhomin = -8.0D+00
      rhomax = +8.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_spherical ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'spherical', 
     &  'Spherical correlation' )
c
c  white_noise
c  1 plot is enough
c
      n2 = 1
      rho0(1) = 1.0
      rhomin = -2.0D+00
      rhomax = +2.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho )
      do j = 1, n2
        call correlation_white_noise ( n, rho, rho0(j), c(1:n,j) )
      end do
      call correlation_plots ( n, n2, rho, rho0, c, 'white_noise', 
     &  'White noise correlation' )

      return
      end
      subroutine correlation_test04 ( )

c*********************************************************************72
c
cc CORRELATION_TEST04 converts between covariance and correlation matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(n,n)
      double precision k(n,n)
      double precision k2(n,n)
      double precision sigma(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TEST04'
      write ( *, '(a)' ) 
     &  '  Convert between a correlation and a covariance matrix.'

      call minij ( n, n, k )

      call r8mat_print ( n, n, k, '  Covariance matrix K:' )

      call covariance_to_correlation ( n, k, c, sigma )

      call r8mat_print ( n, n, c, '  Correlation matrix C:' )

      call r8vec_print ( n, sigma, '  Variances:' )

      call correlation_to_covariance ( n, c, sigma, k2 )

      call r8mat_print ( n, n, k2, '  Recovered covariance matrix K2:' )

      return
      end
      subroutine correlation_test05 ( )

c*********************************************************************72
c
cc CORRELATION_TEST05 calls CORRELATION_BROWNIAN_DISPLAY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'CORRELATION_TEST05'
      write ( *, '(a)' ) 
     &  '  CORRELATION_BROWNIAN_DISPLAY displays 4 slices of'
      write ( *, '(a)' ) '  the Brownian correlation function.'

      call correlation_brownian_display ( )

      return
      end
      subroutine correlation_test06 ( )

c*********************************************************************72
c
cc CORRELATION_TEST06 plots sample paths with SAMPLE_PATHS2_CHOLESKY/EIGEN/FFT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )
      integer n2
      parameter ( n2 = 3 )

      external correlation_brownian
      double precision rho(n)
      double precision rho0
      double precision rhomax
      double precision rhomin
      integer seed
      double precision x(n,n2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CORRELATION_TEST06'
      write ( *, '(a)' ) '  For non-stationary correlation functions:'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 
     &  '  SAMPLE_PATHS2_CHOLESKY generates sample paths from the'
      write ( *, '(a)' ) 
     &  '  correlation matrix, factored using the Cholesky factor.  It'
      write ( *, '(a)' ) 
     &  '  requires that the correlation matrix' //
     &  ' is nonnegative definite.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 
     &  '  SAMPLE_PATHS2_EIGEN generates sample paths from the'
      write ( *, '(a)' ) 
     &  '  correlation matrix, factored using the eigen factorization.';
      write ( *, '(a)' ) 
     &  '  If the correlation matrix is not nonnegative definite,'
      write ( *, '(a)' ) '  we simply suppress negative eigenvalues.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 
     &  '  SAMPLE_PATHS2_FFT generates sample paths from the'
      write ( *, '(a)' ) 
     &  '  correlation matrix, factored using the FFT factorization';
      write ( *, '(a)' ) 
     &  '  of the correlation matrix after embedding in a circulant.'
      write ( *, '(a)' ) ''
c
c  brownian
c
      rhomin = 0.0D+00
      rhomax = 10.0D+00
      rho0 = 1.0D+00
      seed = 123456789
 
      call sample_paths2_cholesky ( n, n2, rhomin, rhomax, rho0, 
     &  correlation_brownian, seed, x )

      call r8vec_linspace ( n, rhomin, rhomax, rho )

      call paths_plot ( n, n2, rho, x, 'brownian', 
     &  'Brownian correlation' )
 
      return
      end


