      program main

c*********************************************************************72
c
cc COMPASS_SEARCH_TEST tests COMPASS_SEARCH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPASS_SEARCH_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the COMPASS_SEARCH library.'

      call beale_test ( )
      call bohach1_test ( )
      call bohach2_test ( )
      call broyden_test ( )
      call extended_rosenbrock_test ( )
      call goldstein_price_test ( )
      call himmelblau_test ( )
      call local_test ( )
      call mckinnon_test ( )
      call powell_test ( )
      call rosenbrock_test ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPASS_SEARCH_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine beale_test ( )

c*********************************************************************72
c
cc BEALE_TEST works with the Beale function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision, external :: beale
      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BEALE_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Beale function.'
      delta_tol = 0.00001D+00
      delta = 0.1D+00
      k_max = 20000

      x0(1) = 1.0D+00
      x0(2) = 1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', beale ( m, x0 )

      call compass_search ( beale, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i10)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Repeat with more difficult start.
c
      x0(1) = 1.0D+00
      x0(2) = 4.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', beale ( m, x0 )

      call compass_search ( beale, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 3.0D+00
      x(2) = 0.5D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', beale ( m, x )

      return
      end
      subroutine bohach1_test ( )

c*********************************************************************72
c
cc BOHACH1_TEST works with the Bohachevsky function #1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision, external :: bohach1
      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOHACH1_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Bohachevsky function #1.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = 0.5D+00
      x0(2) = 1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', bohach1 ( m, x0 )

      call compass_search ( bohach1, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 0.0D+00
      x(2) = 0.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', bohach1 ( m, x )

      return
      end
      subroutine bohach2_test ( )

c*********************************************************************72
c
cc BOHACH2_TEST works with the Bohachevsky function #2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision, external :: bohach2
      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOHACH2_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Bohachevsky function #2.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = 0.6D+00
      x0(2) = 1.3D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', bohach2 ( m, x0 )

      call compass_search ( bohach2, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 0.0D+00
      x(2) = 0.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', bohach2 ( m, x )

      return
      end
      subroutine broyden_test ( )

c*********************************************************************72
c
cc BROYDEN_TEST works with the Broyden function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision, external :: broyden
      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BROYDEN_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Broyden function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = -0.9D+00
      x0(2) = -1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', broyden ( m, x0 )

      call compass_search ( broyden, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = -0.511547141775014D+00
      x(2) = -0.398160951772036D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', broyden ( m, x )

      return
      end
      subroutine extended_rosenbrock_test ( )

c*********************************************************************72
c
cc EXTENDED_ROSENBROCK_TEST works with the extended Rosenbrock function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 4

      double precision delta
      double precision delta_tol
      double precision, external :: extended_rosenbrock
      double precision fx
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTENDED_ROSENBROCK_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the extended Rosenbrock function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = - 1.2D+00
      x0(2) =   1.0D+00
      x0(3) = - 1.5D+00
      x0(4) =   1.2D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  F(X0) = ', extended_rosenbrock ( m, x0 )

      call compass_search ( extended_rosenbrock, m, x0, delta_tol, 
     &  delta, k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 1.0D+00
      x(2) = 1.0D+00
      x(3) = 1.0D+00
      x(4) = 1.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  F(X*) = ', extended_rosenbrock ( m, x )

      return
      end
      subroutine goldstein_price_test ( )

c*********************************************************************72
c
cc GOLDSTEIN_PRICE_TEST works with the Goldstein-Price function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision delta
      double precision delta_tol
      double precision fx
      double precision, external :: goldstein_price
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GOLDSTEIN_PRICE_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Goldstein-Price function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = -0.5D+00
      x0(2) =  0.25D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', goldstein_price ( m, x0 )

      call compass_search ( goldstein_price, m, x0, delta_tol, 
     &  delta, k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Repeat with more difficult start.
c
      x0(1) = -4.0D+00
      x0(2) =  5.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', goldstein_price ( m, x0 )

      call compass_search ( goldstein_price, m, x0, delta_tol, 
     &  delta, k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) =  0.0D+00
      x(2) = -1.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', goldstein_price ( m, x )

      return
      end
      subroutine himmelblau_test ( )

c*********************************************************************72
c
cc HIMMELBLAU_TEST works with the Himmelblau function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision delta
      double precision delta_tol
      double precision fx
      double precision, external :: himmelblau
      integer k
      integer k_max
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HIMMELBLAU_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Himmelblau function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = 1.0D+00
      x0(2) = 1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

      call compass_search ( himmelblau, m, x0, delta_tol, delta, 
     &  k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x0(1) = -1.0D+00
      x0(2) =  1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

      call compass_search ( himmelblau, m, x0, delta_tol, delta, 
     &  k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x0(1) = -1.0D+00
      x0(2) = -1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

      call compass_search ( himmelblau, m, x0, delta_tol, delta, 
     &  k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x0(1) =  1.0D+00
      x0(2) = -1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', himmelblau ( m, x0 )

      call compass_search ( himmelblau, m, x0, delta_tol, delta, 
     &  k_max, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate Himmelblau minimizers.
c
      x(1) = 3.0D+00
      x(2) = 2.0D+00
      call r8vec_print ( m, x, '  Correct Himmelblau minimizer X1*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

      x(1) =  3.58439D+00
      x(2) = -1.84813D+00
      call r8vec_print ( m, x, '  Correct Himmelblau minimizer X2*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

      x(1) = -3.77934D+00
      x(2) = -3.28317D+00
      call r8vec_print ( m, x, '  Correct Himmelblau minimizer X3*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

      x(1) = -2.80512D+00
      x(2) =  3.13134D+00
      call r8vec_print ( m, x, '  Correct Himmelblau minimizer X4*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', himmelblau ( m, x )

      return
      end
      subroutine local_test ( )

c*********************************************************************72
c
cc LOCAL_TEST works with the Local function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision, external :: local
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOCAL_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Local function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = 1.0D+00
      x0(2) = 1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', local ( m, x0 )

      call compass_search ( local, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate local minimizer.
c
      x(1) = 0.2858054412D+00
      x(2) = 0.2793263206D+00
      call r8vec_print ( m, x, '  Correct local minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', local ( m, x )
c
c  Try for global minimizer.
c
      x0(1) = -15.0D+00
      x0(2) = -35.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', local ( m, x0 )

      call compass_search ( local, m, x0, delta_tol, delta, k_max, x, 
     &  fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate global minimizer.
c
      x(1) = -21.02667179D+00
      x(2) = -36.75997872D+00
      call r8vec_print ( m, x, '  Correct global minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', local ( m, x )

      return
      end
      subroutine mckinnon_test ( )

c*********************************************************************72
c
cc MCKINNON_TEST works with the McKinnon function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision a
      double precision b
      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision, external :: mckinnon
      double precision phi
      double precision tau
      double precision theta
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCKINNON_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the McKinnon function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000
c
c  Test 1
c
      a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
      b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

      phi = 10.0D+00
      tau = 1.0D+00
      theta = 15.0D+00

      call mckinnon_parameters ( 'set', phi, tau, theta )

      x0(1) = a
      x0(2) = b
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) 
     &  '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

      call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x(1) =  0.0D+00
      x(2) = -0.5D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )
c
c  Test 2
c
      a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
      b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

      phi = 60.0D+00
      tau = 2.0D+00
      theta = 6.0D+00

      call mckinnon_parameters ( 'set', phi, tau, theta )

      x0(1) = a
      x0(2) = b
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) 
     &  '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

      call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x(1) =  0.0D+00
      x(2) = -0.5D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )
c
c  Test 3
c
      a = ( 1.0D+00 + sqrt ( 33.0D+00 ) ) / 8.0D+00
      b = ( 1.0D+00 - sqrt ( 33.0D+00 ) ) / 8.0D+00

      phi = 4000.0D+00
      tau = 3.0D+00
      theta = 6.0D+00

      call mckinnon_parameters ( 'set', phi, tau, theta )

      x0(1) = a
      x0(2) = b
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) 
     &  '  PHI = ', phi, ' TAU = ', tau, ' THETA = ', theta
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', mckinnon ( m, x0 )

      call compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k

      x(1) =  0.0D+00
      x(2) = -0.5D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', mckinnon ( m, x )

      return
      end
      subroutine powell_test ( )

c*********************************************************************72
c
cc POWELL_TEST works with the Powell function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 4

      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision, external :: powell
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POWELL_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Powell function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) =  3.0D+00
      x0(2) = -1.0D+00
      x0(3) =  0.0D+00
      x0(4) =  1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', powell ( m, x0 )

      call compass_search ( powell, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 0.0D+00
      x(2) = 0.0D+00
      x(3) = 0.0D+00
      x(4) = 0.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', powell ( m, x )

      return
      end
      subroutine rosenbrock_test ( )

c*********************************************************************72
c
cc ROSENBROCK_TEST works with the Rosenbrock function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: m = 2

      double precision delta
      double precision delta_tol
      double precision fx
      integer k
      integer k_max
      double precision, external :: rosenbrock
      double precision x(m)
      double precision x0(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROSENBROCK_TEST:'
      write ( *, '(a)' ) 
     &  '  Test COMPASS_SEARCH with the Rosenbrock function.'
      delta_tol = 0.00001D+00
      delta = 0.3D+00
      k_max = 20000

      x0(1) = -1.2D+00
      x0(2) =  1.0D+00
      call r8vec_print ( m, x0, '  Initial point X0:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X0) = ', rosenbrock ( m, x0 )

      call compass_search ( rosenbrock, m, x0, delta_tol, delta, k_max, 
     &  x, fx, k, x, fx, k )
      call r8vec_print ( m, x, '  Estimated minimizer X1:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i12)' ) 
     &  '  F(X1) = ', fx, ' number of steps = ', k
c
c  Demonstrate correct minimizer.
c
      x(1) = 1.0D+00
      x(2) = 1.0D+00
      call r8vec_print ( m, x, '  Correct minimizer X*:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', rosenbrock ( m, x )

      return
      end
      function beale ( m, x )

c*********************************************************************72
c
cc BEALE computes the Beale function.
c
c  Discussion:
c
c    This function has a global minimizer:
c
c      X = ( 3.0, 0.5 )
c
c    for which
c
c      F(X) = 0.
c
c    For a relatively easy computation, start the search at
c
c      X = ( 1.0, 1.0 )
c
c    A harder computation starts at
c
c      X = ( 1.0, 4.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Evelyn Beale,
c    On an Iterative Method for Finding a Local Minimum of a Function
c    of More than One Variable,
c    Technical Report 25, 
c    Statistical Techniques Research Group,
c    Princeton University, 1958.
c
c    Richard Brent,
c    Algorithms for Minimization with Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision BEALE, the value of the function at X.
c
      implicit none

      integer m

      double precision beale
      double precision f
      double precision f1
      double precision f2
      double precision f3
      double precision x(m)

      f1 = 1.5D+00   - x(1) * ( 1.0D+00 - x(2)    )
      f2 = 2.25D+00  - x(1) * ( 1.0D+00 - x(2)**2 )
      f3 = 2.625D+00 - x(1) * ( 1.0D+00 - x(2)**3 )

      f = f1**2 + f2**2 + f3**2
     
      beale = f

      return
      end
      function bohach1 ( m, x )

c*********************************************************************72
c
cc BOHACH1 evaluates the Bohachevsky function #1.
c
c  Discussion:
c
c    The minimizer is
c
c      X* = ( 0.0, 0.0 )
c      F(X*) = 0.0
c
c    Suggested starting point:
c
c      X = ( 0.5, 1.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Zbigniew Michalewicz,
c    Genetic Algorithms + Data Structures = Evolution Programs,
c    Third Edition,
c    Springer Verlag, 1996,
c    ISBN: 3-540-60676-9,
c    LC: QA76.618.M53.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision BOHACH1, the value of the function at X.
c
      implicit none

      integer m

      double precision bohach1
      double precision f
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision x(m)

      f =           x(1) * x(1) - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) 
     &  + 2.0D+00 * x(2) * x(2) - 0.4D+00 * cos ( 4.0D+00 * pi * x(2) ) 
     &  + 0.7D+00

      bohach1 = f

      return
      end
      function bohach2 ( m, x )

c*********************************************************************72
c
cc BOHACH2 evaluates the Bohachevsky function #2.
c
c  Discussion:
c
c    The minimizer is
c
c      X* = ( 0.0, 0.0 )
c      F(X*) = 0.0
c
c    Suggested starting point:
c
c      X = ( 0.6, 1.3 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Zbigniew Michalewicz,
c    Genetic Algorithms + Data Structures = Evolution Programs,
c    Third Edition,
c    Springer Verlag, 1996,
c    ISBN: 3-540-60676-9,
c    LC: QA76.618.M53.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision BOHACH2, the value of the function at X.
c
      implicit none

      integer m

      double precision bohach2
      double precision f
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(m)

      f =           x(1) * x(1) 
     &  + 2.0D+00 * x(2) * x(2)
     &  - 0.3D+00 * cos ( 3.0D+00 * pi * x(1) ) 
     &  * cos ( 4.0D+00 * pi * x(2) ) 
     &  + 0.3D+00

      bohach2 = f

      return
      end
      function broyden ( m, x )

c*********************************************************************72
c
cc BROYDEN computes the two dimensional modified Broyden function.
c
c  Discussion:
c
c    This function has a global minimizer:
c
c      X = ( -0.511547141775014, -0.398160951772036 )
c
c    for which
c
c      F(X) = 1.44E-04
c
c    Start the search at
c
c      X = ( -0.9, -1.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Broyden,
c    A class of methods for solving nonlinear simultaneous equations,
c    Mathematics of Computation,
c    Volume 19, 1965, pages 577-593.
c
c    Jorge More, Burton Garbow, Kenneth Hillstrom,
c    Testing unconstrained optimization software,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 1, March 1981, pages 17-41. 
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision BROYDEN, the value of the function at X.
c
      implicit none

      integer m

      double precision broyden
      double precision f
      double precision f1
      double precision f2
      double precision p
      double precision x(m)

      f1 = abs ( ( 3.0D+00 -           x(1) ) * x(1) 
     &  - 2.0D+00 * x(2) + 1.0D+00 )
      f2 = abs ( ( 3.0D+00 - 2.0D+00 * x(2) ) * x(2) 
     &  -           x(1) + 1.0D+00 )

      p = 3.0D+00 / 7.0D+00

      f = f1**p + f2**p
     
      broyden = f

      return
      end
      function extended_rosenbrock ( m, x )

c*********************************************************************72
c
cc EXTENDED_ROSENBROCK computes the extended Rosenbrock function.
c
c  Discussion:
c
c    The number of dimensions is arbitrary, except that it must be even.
c
c    There is a global minimum at X* = (1,1,&), F(X*) = 0.
c
c    The contours are sharply twisted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Howard Rosenbrock,
c    An Automatic Method for Finding the Greatest or Least Value of a Function,
c    Computer Journal,
c    Volume 3, 1960, pages 175-184.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision EXTENDED_ROSENBROCK, the value of the function at X.
c
      implicit none

      integer m

      double precision extended_rosenbrock
      double precision f
      double precision f1
      double precision f2
      integer i
      double precision x(m)

      f = 0.0D+00

      do i = 1, m - 1, 2
        f1 = 1.0D+00 - x(i)
        f2 = 10.0D+00 * ( x(i+1) - x(i)**2 )
        f = f + f1**2 + f2**2
      end do

      extended_rosenbrock = f

      return
      end
      function goldstein_price ( m, x )

c*********************************************************************72
c
cc GOLDSTEIN_PRICE evaluates the Goldstein-Price polynomial.
c
c  Discussion:
c
c    The minimizer is
c
c      X* = ( 0.0, -1.0 )
c      F(X*) = 3.0
c
c    Suggested starting point:
c
c      X = ( -0.5, 0.25 ) (easy convergence)
c      X = ( -4.0, 5.00 ) (harder convergence)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Zbigniew Michalewicz,
c    Genetic Algorithms + Data Structures = Evolution Programs,
c    Third Edition,
c    Springer Verlag, 1996,
c    ISBN: 3-540-60676-9,
c    LC: QA76.618.M53.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision GOLDSTEIN_PRICE, the value of the function at X.
c
      implicit none

      integer m

      double precision a
      double precision b
      double precision c
      double precision d
      double precision f
      double precision goldstein_price
      double precision x(m)

      a = x(1) + x(2) + 1.0D+00

      b = 19.0D+00 - 14.0D+00 * x(1) + 3.0D+00 * x(1) * x(1) 
     &  - 14.0D+00 * x(2) + 6.0D+00 * x(1) * x(2) 
     &  + 3.0D+00 * x(2) * x(2)

      c = 2.0D+00 * x(1) - 3.0D+00 * x(2)

      d = 18.0D+00 - 32.0D+00 * x(1) + 12.0D+00 * x(1) * x(1) 
     &  + 48.0D+00 * x(2) - 36.0D+00 * x(1) * x(2) 
     &  + 27.0D+00 * x(2) * x(2)

      f = ( 1.0D+00 + a * a * b ) * ( 30.0D+00 + c * c * d )

      goldstein_price = f

      return
      end
      function himmelblau ( m, x )

c*********************************************************************72
c
cc HIMMELBLAU computes the Himmelblau function.
c
c  Discussion:
c
c    This function has 4 global minima:
c
c      X* = (  3,        2       ), F(X*) = 0.
c      X* = (  3.58439, -1.84813 ), F(X*) = 0.
c      X* = ( -3.77934, -3.28317 ), F(X*) = 0.
c      X* = ( -2.80512,  3.13134 ), F(X*) = 0.
c
c    Suggested starting points are
c
c      (+1,+1),
c      (-1,+1),
c      (-1,-1),
c      (+1,-1),
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Himmelblau,
c    Applied Nonlinear Programming,
c    McGraw Hill, 1972,
c    ISBN13: 978-0070289215,
c    LC: T57.8.H55.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision HIMMELBLAU, the value of the function at X.
c
      implicit none

      integer m

      double precision f
      double precision himmelblau
      double precision x(m)

      f = ( x(1)**2 + x(2) - 11.0D+00 )**2 
     &  + ( x(1) + x(2)**2 - 7.0D+00 )**2

      himmelblau = f
     
      return
      end
      function local ( m, x )

c*********************************************************************72
c
cc LOCAL computes the local function.
c
c  Discussion:
c
c    This function has a local minimizer:
c
c      X* = ( 0.2858054412&, 0.2793263206&), F(X*) = 5.9225&
c
c    and a global minimizer:
c
c      X* = ( -21.02667179&, -36.75997872&), F(X*) = 0.
c
c    Suggested starting point for local minimizer:
c
c      X = ( 1, 1 ), F(X) = 3.33 * 10^6.
c
c    Suggested starting point for global minimizer:
c
c      X = ( -15, -35), F(X) = 1.49 * 10^8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Himmelblau,
c    Applied Nonlinear Programming,
c    McGraw Hill, 1972,
c    ISBN13: 978-0070289215,
c    LC: T57.8.H55.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision LOCAL, the value of the function at X.
c
      implicit none

      integer m

      double precision f
      double precision local
      double precision x(m)

      f = ( x(1)**2 + 12.0D+00 * x(2) - 1.0D+00 )**2 
     &  + ( 49.0D+00 * x(1)**2 + 49.0D+00 * x(2)**2 + 84.0D+00 * x(1) 
     &  + 2324.0D+00 * x(2) - 681.0D+00 )**2
     
      local = f

      return
      end
      function mckinnon ( m, x )

c*********************************************************************72
c
cc MCKINNON computes the McKinnon function.
c
c  Discussion:
c
c    This function has a global minimizer:
c
c      X* = ( 0.0, -0.5 ), F(X*) = -0.25.
c
c    There are three parameters, TAU, THETA and PHI.
c
c    1 < TAU, then F is strictly convex.
c             and F has continuous first derivatives.
c    2 < TAU, then F has continuous second derivatives.
c    3 < TAU, then F has continuous third derivatives.
c
c    However, this function can cause the Nelder-Mead optimization
c    algorithm to "converge" to a point which is not the minimizer
c    of the function F.
c
c    Sample parameter values which cause problems for Nelder-Mead 
c    include:
c
c      PHI = 10,  TAU = 1, THETA = 15
c      PHI = 60,  TAU = 2, THETA =  6
c      PHI = 400, TAU = 3, THETA =  6
c
c    To get the bad behavior, we also assume the initial simplex has the form
c
c      X1 = (0,0),
c      X2 = (1,1),
c      X3 = (A,B), 
c
c    where 
c
c      A = (1+sqrt(33))/8 =  0.84307&
c      B = (1-sqrt(33))/8 = -0.59307&
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ken McKinnon,
c    Convergence of the Nelder-Mead simplex method to a nonstationary point,
c    SIAM Journal on Optimization,
c    Volume 9, Number 1, 1998, pages 148-158.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision MCKINNON, the value of the function at X.
c
      implicit none

      integer m

      double precision f
      double precision mckinnon
      double precision phi
      double precision tau
      double precision theta
      double precision x(m)

      call mckinnon_parameters ( 'get', phi, tau, theta )

      if ( x(1) .le. 0.0D+00 ) then
        f = theta * phi * abs ( x(1) )**tau + x(2) * ( 1.0D+00 + x(2) )
      else
        f = theta       *       x(1)**tau   + x(2) * ( 1.0D+00 + x(2) )
      end if

      mckinnon = f

      return
      end
      subroutine mckinnon_parameters ( action, phi, tau, theta )

c*********************************************************************72
c
cc MCKINNON_PARAMETERS sets or gets McKinnon function parameters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.  
c    'S', set the parameters.
c    'G', get the parameters.
c
c    Input/output, double precision PHI, TAU, THETA, the parameter values.
c
      implicit none

      character ( len = * ) action
      double precision phi
      double precision, save :: phi_save = 60.0D+00
      double precision tau
      double precision, save :: tau_save = 2.0D+00
      double precision theta
      double precision, save :: theta_save = 6.0D+00

      if ( action(1:1) .eq. 'S' .or. action(1:1) .eq. 's' ) then
        phi_save = phi
        tau_save = tau
        theta_save = theta
      else if ( action(1:1) .eq. 'G' .or. action(1:1) .eq. 'g' ) then
        phi = phi_save
        tau = tau_save
        theta = theta_save
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MCKINNON_PARAMETERS - Fatal error!'
        write ( *, '(a)' ) '  Unexpected value of ACTION.'
        stop
      end if

      return
      end
      function powell ( m, x )

c*********************************************************************72
c
cc POWELL computes the Powell singular quartic function.
c
c  Discussion:
c
c    This function has a global minimizer:
c
c      X* = ( 0.0, 0.0, 0.0, 0.0 ), F(X*) = 0.
c
c    Start the search at
c
c      X = ( 3.0, -1.0, 0.0, 1.0 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Powell,
c    An Iterative Method for Finding Stationary Values of a Function
c    of Several Variables,
c    Computer Journal,
c    Volume 5, 1962, pages 147-151.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision POWELL, the value of the function at X.
c
      implicit none

      integer m

      double precision f
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision powell
      double precision x(m)

      f1 = x(1) + 10.0D+00 * x(2)
      f2 =                                x(3) - x(4)
      f3 =               x(2) - 2.0D+00 * x(3)
      f4 = x(1)                                - x(4)

      f = f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4
     
      powell = f

      return
      end
      function rosenbrock ( m, x )

c*********************************************************************72
c
cc ROSENBROCK computes the Rosenbrock function.
c
c  Discussion:
c
c    There is a global minimum at X* = (1,1), F(X*) = 0.
c
c    The starting point X = ( -1.2, 1.0 ) is suggested.
c
c    The contours are sharply twisted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Howard Rosenbrock,
c    An Automatic Method for Finding the Greatest or Least Value of a Function,
c    Computer Journal,
c    Volume 3, 1960, pages 175-184.
c
c  Parameters:
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X(M), the argument of the function.
c
c    Output, double precision ROSENBROCK, the value of the function at X.
c
      implicit none

      integer m

      double precision f
      double precision rosenbrock
      double precision x(m)

      f = ( 1.0D+00 - x(1) )**2 + 100.0D+00 * ( x(2) - x(1) * x(1) )**2

      rosenbrock = f

      return
      end
