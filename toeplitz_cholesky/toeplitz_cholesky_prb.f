      program main

c*********************************************************************72
c
cc MAIN tests TOEPLITZ_CHOLESKY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOEPLITZ_CHOLESKY library.'

      call toeplitz_cholesky_prb01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      stop
      end
      subroutine toeplitz_cholesky_prb01 ( )

c*********************************************************************72
c
cc TOEPLITZ_CHOLESKY_PRB01 tests TOEPLITZ_CHOLESKY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision a(n,n)
      double precision b(n,n)
      double precision g(2,n)
      double precision l(n,n)
      double precision r(n,n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB01'
      write ( *, '(a)' ) '  Test the factorization of a simple example.'
c
c  TOEPLITZ_CHOLESKY_UPPER.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_UPPER:'

      a(1,1) =  1.0
      a(2,1) =  0.5
      a(3,1) = -0.375

      a(1,2) =  0.5
      a(2,2) =  1.0
      a(3,2) =  0.5

      a(1,3) = -0.375
      a(2,3) =  0.5
      a(3,3) =  1.0

      call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

      call toeplitz_cholesky_upper ( n, a, r )
      call r8mat_print ( n, n, r, 
     &  '  Computed upper Cholesky factor R:' )

      call r8mat_mtm ( n, n, n, r, r, b )
      call r8mat_print ( n, n, b, '  Product R''R:' )
c
c  TOEP_CHOLESKY_UPPER.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEP_CHOLESKY_UPPER:'
  
      g(1,1) =  1.0
      g(2,1) =  0.0

      g(1,2) =  0.5
      g(2,2) =  0.5

      g(1,3) = -0.375
      g(2,3) = -0.375

      call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

      call toep_cholesky_upper ( n, g, r )
      call r8mat_print ( n, n, r, 
     &  '  Computed upper Cholesky factor R:' )

      call r8mat_mtm ( n, n, n, r, r, b )
      call r8mat_print ( n, n, b, '  Product R''R:' )
c
c  TOEPLITZ_CHOLESKY_LOWER.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_LOWER:'

      a(1,1) =  1.0
      a(2,1) =  0.5
      a(3,1) = -0.375

      a(1,2) =  0.5
      a(2,2) =  1.0
      a(3,2) =  0.5

      a(1,3) = -0.375
      a(2,3) =  0.5
      a(3,3) =  1.0

      call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

      call toeplitz_cholesky_lower ( n, a, l )
      call r8mat_print ( n, n, l, 
     &  '  Computed lower Cholesky factor L:' )

      call r8mat_mmt ( n, n, n, l, l, b )
      call r8mat_print ( n, n, b, '  Product LL'':' )
c
c  TOEP_CHOLESKY_LOWER.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TOEP_CHOLESKY_LOWER:'

      g(1,1) =  1.0
      g(2,1) =  0.0

      g(1,2) =  0.5
      g(2,2) =  0.5

      g(1,3) = -0.375
      g(2,3) = -0.375

      call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

      call toep_cholesky_lower ( n, g, l )
      call r8mat_print ( n, n, l, 
     &  '  Computed lower Cholesky factor L:' )

      call r8mat_mmt ( n, n, n, l, l, b )
      call r8mat_print ( n, n, b, '  Product LL'':' )

      return
      end



