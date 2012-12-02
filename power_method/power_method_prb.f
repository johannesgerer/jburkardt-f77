      program main

c*********************************************************************72
c
cc MAIN is the main program for POWER_METHOD_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POWER_METHOD_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version:'
      write ( *, '(a)' ) '  Test the POWER_METHOD library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POWER_METHOD_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 uses POWER_METHOD on the Fibonacci2 matrix.
c
c  Discussion:
c
c    This matrix, despite having a single dominant eigenvalue, will generally
c    converge only very slowly under the power method.  This has to do with
c    the fact that the matrix has only 3 eigenvectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 50 )

      double precision a(n,n)
      double precision cos_x1x2
      double precision ctime
      double precision ctime1
      double precision ctime2
      integer i
      integer it_max
      integer it_num
      double precision lambda
      double precision norm
      double precision phi
      integer seed
      double precision sin_x1x2
      double precision tol
      double precision x(n)
      double precision x2(n)

      call fibonacci2 ( n, a )

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, x )

      it_max = 300
      tol = 0.000001D+0

      phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' )
     &  '  Use the power method on the Fibonacci2 matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Matrix order N       = ', n
      write ( *, '(a,i8)'      ) '  Maximum iterations   = ', it_max
      write ( *, '(a,g14.6)'   ) '  Error tolerance      = ', tol

      call cpu_time ( ctime1 )

      call power_method ( n, a, x, it_max, tol, lambda, it_num )

      call cpu_time ( ctime2 )
      ctime = ctime2 - ctime1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Number of iterations = ', it_num
      write ( *, '(a,g14.6)'   ) '  CPU time             = ', ctime
      write ( *, '(a,f14.10)'  ) '  Estimated eigenvalue = ', lambda
      write ( *, '(a,f14.10)'  ) '  Correct value        = ', phi
      write ( *, '(a,g14.6)'   )
     &  '  Error                = ', abs ( lambda - phi )
c
c  X2 is the exact eigenvector.
c
      x2(1) = 1.0
      do i = 2, n
        x2(i) = phi * x2(i-1)
      end do

      norm = 0.0
      do i = 1, n
        norm = norm + x2(i) * x2(i)
      end do
      norm = sqrt ( norm )

      do i = 1, n
        x2(i) = x2(i) / norm
      end do
c
c  The sine of the angle between X and X2 is a measure of error.
c
      cos_x1x2 = 0.0
      do i = 1, n
        cos_x1x2 = cos_x1x2 + x(i) * x2(i)
      end do

      sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f14.10)' )
     &  '  Sine of angle between true and estimated vectors = ',
     &  sin_x1x2

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 uses POWER_METHOD2 on the Fibonacci2 matrix.
c
c  Discussion:
c
c    This matrix, despite having a single dominant eigenvalue, will generally
c    converge only very slowly under the power method.  This has to do with
c    the fact that the matrix has only 3 eigenvectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 50 )

      double precision a(n,n)
      double precision alpha
      double precision beta
      double precision cos_y1y2
      double precision ctime
      double precision ctime1
      double precision ctime2
      double precision gamma
      integer i
      integer it_max
      integer it_num
      double complex lambda
      double precision phi
      integer seed
      double precision sin_y1y2
      double precision tol
      double complex v(n)
      double precision x(n)

      call fibonacci2 ( n, a )

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, x )

      it_max = 300
      tol = 0.000001D+0

      phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  Use the power method2 on the Fibonacci2 matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Matrix order N       = ', n
      write ( *, '(a,i8)'      ) '  Maximum iterations   = ', it_max
      write ( *, '(a,g14.6)'   ) '  Error tolerance      = ', tol

      call cpu_time ( ctime1 )

      call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )

      call cpu_time ( ctime2 )
      ctime = ctime2 - ctime1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Number of iterations = ', it_num
      write ( *, '(a,g14.6)'   ) '  CPU time             = ', ctime
      write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue = ', lambda
      write ( *, '(a,f14.10)'  ) '  Correct value        = ', phi
      write ( *, '(a,g14.6)'   )
     &  '  Error                = ', abs ( lambda - phi )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 uses POWER_METHOD2 on the TRIS matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 50 )

      double precision a(n,n)
      double precision alpha
      double precision beta
      double precision cos_y1y2
      double precision ctime
      double precision ctime1
      double precision ctime2
      double precision gamma
      integer i
      integer it_max
      integer it_num
      double complex lambda
      double complex lambda_max
      double complex lambda_vec(n)
      double precision phi
      integer seed
      double precision sin_y1y2
      double precision tol
      double complex v(n)
      double precision x(n)
c
c  If ALPHA * GAMMA is negative, this matrix will have complex eigenvalues.
c
      alpha = -1.0D+00
      beta = 10.0D+00
      gamma = +8.0D+0

      call tris ( n, n, alpha, beta, gamma, a )

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, x )

      it_max = 4000
      tol = 0.000001D+0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  Use the power method2 on the TRIS matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Matrix order N       = ', n
      write ( *, '(a,i8)'      ) '  Maximum iterations   = ', it_max
      write ( *, '(a,g14.6)'   ) '  Error tolerance      = ', tol

      call cpu_time ( ctime1 )

      call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )

      call cpu_time ( ctime2 )
      ctime = ctime2 - ctime1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'      ) '  Number of iterations = ', it_num
      write ( *, '(a,g14.6)'   ) '  CPU time             = ', ctime
      write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue = ', lambda

      call tris_eigenvalues ( n, alpha, beta, gamma, lambda_vec )

      lambda_max = lambda_vec(1)
      do i = 2, n
        if ( abs ( lambda_max ) .lt. abs ( lambda_vec(i) ) ) then
          lambda_max = lambda_vec(i)
        end if
      end do

      write ( *, '(a,2f14.10)' ) '  Correct value        = ', lambda_max
      write ( *, '(a,g14.6)'   )
     &  '  Error                = ', abs ( lambda - lambda_max )

      return
      end
