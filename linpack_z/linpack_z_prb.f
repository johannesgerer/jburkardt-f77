      program main

c*********************************************************************72
c
cc MAIN is the main program for LINPACK_Z_PRB.
c
c  Discussion:
c
c    LINPACK_Z_PRB tests the double precision complex LINPACK test routines.
c
c  Modified:
c
c    03 May 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINPACK_Z_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LINPACK_Z library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
      call test22 ( )
      call test23 ( )
      call test24 ( )
      call test25 ( )
      call test26 ( )
      call test27 ( )
      call test28 ( )
      call test29 ( )

      call test30 ( )
      call test31 ( )
      call test32 ( )
      call test33 ( )
      call test34 ( )
      call test345 ( )
      call test35 ( )
      call test36 ( )
      call test37 ( )
!
!  Terminate.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINPACK_Z_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ZCHDC.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(lda,n)
      double complex b(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      double complex work(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a double complex Hermitian '
      write ( *, '(a)' ) '  positive definite matrix,'
      write ( *, '(a)' ) '  ZCHDC computes the Cholesky decomposition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00 )
      a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00 )
      a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00 )

      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00 )
      a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00 )

      a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Decompose the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Decompose the matrix.'

      job = 0
      do i = 1, n
        ipvt(i) = 0
      end do

      call zchdc ( a, lda, n, work, ipvt, job, info )

      if ( info .ne. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  ZCHDC returned INFO = ', info
        write ( *, '(a)' )
     &    '  The matrix is not Hermitian positive definite.'
      end if
c
c  Zero out the lower diagonal.
c
      do i = 2, n
        do j = 1, i-1
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Compute the Cholesky product.
c
      do i = 1, n
        do j = 1, n
          b(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            b(i,j) = b(i,j) + conjg ( a(k,i) ) * a(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product U^H * U: '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests ZCHEX.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer nz

      parameter ( n = 3 )
      parameter ( lda = n )
      parameter ( nz = 1 )

      double complex a(lda,n)
      double complex b(lda,n)
      double complex c(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer l
      double complex s(n)
      integer seed
      double complex work(n)
      double complex z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For a double complex Hermitian'
      write ( *, '(a)' ) '  positive definite matrix,'
      write ( *, '(a)' )
     &  '  ZCHEX can shift columns in a Cholesky factorization.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00 )
      a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00 )
      a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00 )

      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00 )
      a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00 )

      a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do

      do i = 1, n
        z(i) = cmplx ( i, 0.0D+00 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector Z:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,2g14.6)' ) z(i)
      end do
c
c  Decompose the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Decompose the matrix.'

      job = 0
      do i = 1, n
        ipvt(i) = 0
      end do

      call zchdc ( a, lda, n, work, ipvt, job, info )

      if ( info .ne. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  ZCHDC returned INFO = ', info
        write ( *, '(a)' ) '  The matrix is not positive definite.'
      end if
c
c  Zero out the lower diagonal.
c
      do i = 2, n
        do j = 1, i-1
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Right circular shift columns L through K.
c
      k = 1
      l = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' )
     &  '  Right circular shift columns K  = ', k,
     &  ' through L = ', l

      job = 1
      call zchex ( a, lda, n, k, l, z, n, nz, c, s, job )
c
c  Left circular shift columns K+1 through L.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' )
     &  '  Left circular shift columns K+1 = ', k+1,
     &  ' through L = ', l

      job = 2
      call zchex ( a, lda, n, k+1, l, z, n, nz, c, s, job )
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted vector Z:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,2g14.6)' ) z(i)
      end do
c
c  Compute the Cholesky product.
c
      do i = 1, n
        do j = 1, n
          b(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            b(i,j) = b(i,j) + conjg ( a(k,i) ) * a(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted product U'' * U: '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests ZCHUD and ZTRSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer p
      integer ldr
      integer nz

      parameter ( p = 20 )
      parameter ( ldr = p )
      parameter ( nz = 1 )

      double complex b(p)
      double complex beta(p)
      double complex c(p)
      integer i
      integer info
      integer j
      integer job
      double complex r(ldr,p)
      double complex rho(nz)
      double complex row(p)
      double complex s(p)
      integer seed
      double complex x(p)
      double complex y(nz)
      double complex z(p,nz)
      double complex zdotu

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a double complex Hermitian matrix'
      write ( *, '(a)' ) '  ZCHUD updates a Cholesky decomposition.'
      write ( *, '(a)' ) '  ZTRSL solves a triangular linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this example, we use ZCHUD to solve a'
      write ( *, '(a)' ) '  least squares problem R * b = z.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of equations is P = ', p
c
c  Initialize.
c
      do j = 1, p
        do i = 1, p
          r(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      do j = 1, nz
        do i = 1, p
          z(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      do i = 1, p
        x(i) = cmplx ( i, mod ( i, 2 ) )
      end do
c
c  Use ZCHUD to form R, Z and RHO by adding X and Y a row at a time.
c  X is a row of the least squares matrix and Y the right hand side.
c
      seed = 123456789

      do i = 1, p
        call zmat_uniform_01 ( 1, p, seed, row )
        y(1) = zdotu ( p, row, 1, x, 1 )
        rho(1) = cmplx ( 0.0D+00, 0.0D+00 )
        call zchud ( r, ldr, p, row, z, p, nz, y, rho, c, s )
      end do
c
c  Generate the least squares solution, b = inverse ( R ) * Z.
c
      do j = 1, nz

        do i = 1, p
          b(i) = z(i,j)
        end do

        job = 01

        call ztrsl ( r, ldr, p, b, job, info )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Solution vector # ', j
        write ( *, '(a)' ) '  (Should be (1,1) (2,0), (3,1) (4,0) ...)'
        write ( *, '(a)' ) ' '

        do i = 1, p
          if ( i .le. 5 .or. p-5 .lt. i ) then
            write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
          end if
          if ( i .eq. 5 ) then
            write ( *, '(a)' ) '  ......  ..............'
          end if
        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests ZGBCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ml
      integer mu
      integer n
      integer lda

      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( n = 3 )

      parameter ( lda = 2*ml+mu+1 )

      double complex a(lda,n)
      double complex a_save(n,n)
      integer i
      integer i1
      integer i2
      integer ipvt(n)
      integer j
      integer k
      integer m
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For a double complex '
      write ( *, '(a)' ) '  general band storage matrix:'
      write ( *, '(a)' ) '  ZGBCO factors the matrix and estimates the'
      write ( *, '(a)' ) '  reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
      write ( *, '(a,i6)' ) '  The lower band is ML =  ', ml
      write ( *, '(a,i6)' ) '  The upper band is MU =  ', mu
c
c  Set the values of the matrix A.
c
      do j = 1, n
        do i = 1, n
          a_save(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      m = ml + mu + 1

      seed = 123456789

      do j = 1, n
        i1 = max ( 1, j - mu )
        i2 = min ( n, j + ml )
        do i = i1, i2
          k = i - j + m
          a(k,j) = z_uniform_01 ( seed )
          a_save(i,j) = a(k,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zgbco ( a, lda, n, ml, mu, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests ZGBFA and ZGBSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ml
      integer mu
      integer n
      integer lda

      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( n = 3 )

      parameter ( lda = 2*ml+mu+1 )

      double complex a(lda,n)
      double complex a_save(n,n)
      double complex b(n)
      integer i
      integer i1
      integer i2
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer m
      integer seed
      double complex x(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  general band storage matrix:'
      write ( *, '(a)' ) '  ZGBFA factors the matrix;'
      write ( *, '(a)' ) '  ZGBSL solves a factored linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
      write ( *, '(a,i6)' ) '  The lower band is ML =  ', ml
      write ( *, '(a,i6)' ) '  The upper band is MU =  ', mu
c
c  Set the values of the matrix A.
c
      do j = 1, n
        do i = 1, n
          a_save(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      m = ml + mu + 1
      seed = 123456789

      do j = 1, n
        i1 = max ( 1, j - mu )
        i2 = min ( n, j + ml )
        do i = i1, i2
          k = i - j + m
          a(k,j) = z_uniform_01 ( seed )
          a_save(i,j) = a(k,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      do i = 1, n
        x(i) = z_uniform_01 ( seed )
      end do

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a_save(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zgbfa ( a, lda, n, ml, mu, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  ZGBFA returned INFO = ', info
        return
      end if
c
c  Solve the system.
c
      job = 0
      call zgbsl ( a, lda, n, ml, mu, ipvt, b, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests ZGBFA and ZGBDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ml
      integer mu
      integer n
      integer lda

      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( n = 3 )

      parameter ( lda = 2*ml+mu+1 )

      double complex a(lda,n)
      double complex a_save(n,n)
      double complex det(2)
      integer i
      integer i1
      integer i2
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer m
      integer seed
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  general band storage matrix:'
      write ( *, '(a)' ) '  ZGBFA factors the matrix.'
      write ( *, '(a)' ) '  ZGBDI computes the determinant.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
      write ( *, '(a,i6)' ) '  The lower band is ML =  ', ml
      write ( *, '(a,i6)' ) '  The upper band is MU =  ', mu
c
c  Set the values of the matrix A.
c
      do j = 1, n
        do i = 1, n
          a_save(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      m = ml + mu + 1
      seed = 123456789

      do j = 1, n
        i1 = max ( 1, j - mu )
        i2 = min ( n, j + ml )
        do i = i1, i2
          k = i - j + m
          a(k,j) = z_uniform_01 ( seed )
          a_save(i,j) = a(k,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zgbfa ( a, lda, n, ml, mu, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  ZGBFA returned INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      call zgbdi ( a, lda, n, ml, mu, ipvt, det )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', dble ( det(2) )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests ZGECO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      integer i
      integer ipvt(n)
      integer j
      integer lda
      double precision rcond
      integer seed
      double complex z(n)

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  general storage matrix:'
      write ( *, '(a)' ) '  ZGECO factors the matrix and estimates the'
      write ( *, '(a)' ) '  reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      call zmat_uniform_01 ( n, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zgeco ( a, lda, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests ZGEFA and ZGESL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex b(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer lda
      integer seed
      double complex x(n)

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  general storage matrix:'
      write ( *, '(a)' ) '  ZGEFA factors the matrix.'
      write ( *, '(a)' ) '  ZGESL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      call zmat_uniform_01 ( n, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      call zmat_uniform_01 ( n, 1, seed, x )

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zgefa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZGEFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the system.
c
      job = 0
      call zgesl ( a, lda, n, ipvt, b, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests ZGEFA and ZGEDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex a_save(n,n)
      double complex c(n,n)
      double complex det(2)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer lda
      integer seed
      double complex work(n)

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  general storage matrix:'
      write ( *, '(a)' ) '  ZGEFA factors the matrix.'
      write ( *, '(a)' ) '  ZGEDI computes the determinant or inverse.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      call zmat_uniform_01 ( n, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do

      do j = 1, n
        do i = 1, n
          a_save(i,j) = a(i,j)
        end do
      end do
c
c  Factor the matrix A.
c
      call zgefa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZGEFA returned an error flag INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      job = 10
      call zgedi ( a, lda, n, ipvt, det, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6,a,g14.6)' )
     &   '  Determinant = ', det(1), ' * 10** ', dble ( det(2) )
c
c  Get the inverse.
c
      job = 01
      call zgedi ( a, lda, n, ipvt, det, work, job )

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + a(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests ZGTSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double complex b(n)
      double complex c(n)
      double complex d(n)
      double complex e(n)
      integer i
      integer info
      integer seed
      double complex x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  tridiagonal matrix:'
      write ( *, '(a)' ) '  ZGTSL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      seed = 123456789

      c(1) = cmplx ( 0.0D+00, 0.0D+00 )
      call zmat_uniform_01 ( n-1, 1, seed, c(2) )

      call zmat_uniform_01 ( n-1, 1, seed, e )
      e(n) = cmplx ( 0.0D+00, 0.0D+00 )

      do i = 1, n
        d(i) = cmplx ( 0.0D+00, 0.0D+00 )
      end do

      do i = 1, n-1
        d(i) = d(i) - 2.0D+00 * e(i)
      end do

      do i = 2, n
        d(i) = d(i) - 2.0D+00 * c(i)
      end do
c
c  Set the desired solution
c
      do i = 1, n
        x(i) = cmplx ( i, 10 * i )
      end do
c
c  Compute the corresponding right hand side.
c
      b(1) = d(1) * x(1) + e(1) * x(2)
      do i = 2, n-1
        b(i) = c(i) * x(i-1) + d(i) * x(i) + e(i) * x(i+1)
      end do
      b(n) = c(n) * x(n-1) + d(n) * x(n)
c
c  Solve the tridiagonal system.
c
      call zgtsl ( n, c, d, e, b, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests ZHICO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double precision d_uniform_01
      integer i
      integer ipvt(n)
      integer j
      integer lda
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  For a double complex Hermitian matrix:'
      write ( *, '(a)' ) '  ZHICO factors the matrix and estimates'
      write ( *, '(a)' ) '  the reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = conjg ( a(i,j) )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zhico ( a, lda, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests ZHIFA and ZHISL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex b(n)
      double precision d_uniform_01
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer lda
      integer seed
      double complex x(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  For a double complex Hermitian matrix:'
      write ( *, '(a)' ) '  ZHIFA factors the matrix.'
      write ( *, '(a)' ) '  ZHISL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = conjg ( a(i,j) )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      call zmat_uniform_01 ( n, 1, seed, x )

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zhifa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZHIFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the system.
c
      call zhisl ( a, lda, n, ipvt, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests ZHIFA and ZHIDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex a_save(n,n)
      double complex c(n,n)
      double precision d_uniform_01
      double precision det(2)
      integer i
      integer inert(3)
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer lda
      integer seed
      double complex work(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  For a double complex hermitian matrix:'
      write ( *, '(a)' ) '  ZHIFA factors the matrix.'
      write ( *, '(a)' ) '  ZHIDI computes the determinant, inverse,'
      write ( *, '(a)' ) '  or inertia.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = conjg ( a(i,j) )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do

      do j = 1, n
        do i = 1, n
          a_save(i,j) = a(i,j)
        end do
      end do
c
c  Factor the matrix A.
c
      call zhifa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZHIFA returned an error flag INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      job = 010
      call zhidi ( a, lda, n, ipvt, det, inert, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', det(2)
c
c  Get the inertia.
c
      job = 100
      call zhidi ( a, lda, n, ipvt, det, inert, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The inertia:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2x,i6)' ) inert(i)
      end do
c
c  Get the inverse.
c
      job = 001
      call zhidi ( a, lda, n, ipvt, det, inert, work, job )
c
c  Only the upper triangle is set, so the user must set up the
c  lower triangle:
c
      do i = 1, n
        do j = 1, i-1
          a(i,j) = conjg ( a(j,i) )
        end do
      end do

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + a(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests ZHPCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double precision d_uniform_01
      integer i
      integer ipvt(n)
      integer j
      integer k
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  For a double complex Hermitian matrix'
      write ( *, '(a)' ) '  using packed storage,'
      write ( *, '(a)' ) '  ZHPCO factors the matrix and estimates'
      write ( *, '(a)' ) '  the reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n

      seed = 123456789
c
c  Set the values of the matrix A.
c
      k = 0

      do j = 1, n

        do i = 1, j-1
          k = k + 1
          a(k) = z_uniform_01 ( seed )
          a_save(i,j) = a(k)
          a_save(j,i) = conjg ( a(k) )
        end do

        k = k + 1
        a(k) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        a_save(j,j) = a(k)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zhpco ( a, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests ZHPFA and ZHPSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double complex b(n)
      double precision d_uniform_01
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer seed
      double complex x(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  For a double complex Hermitian matrix,'
      write ( *, '(a)' ) '  using packed storage,'
      write ( *, '(a)' ) '  ZHPFA factors the matrix.'
      write ( *, '(a)' ) '  ZHPSL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n

      seed = 123456789
c
c  Set the values of the matrix A.
c
      k = 0

      do j = 1, n

        do i = 1, j-1
          k = k + 1
          a(k) = z_uniform_01 ( seed )
          a_save(i,j) = a(k)
          a_save(j,i) = conjg ( a(k) )
        end do

        k = k + 1
        a(k) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        a_save(j,j) = a(k)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      do i = 1, n
        x(i) = z_uniform_01 ( seed )
      end do

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a_save(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zhpfa ( a, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZHPFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the system.
c
      call zhpsl ( a, n, ipvt, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests ZHPFA and ZHPDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double complex b(n,n)
      double complex c(n,n)
      double precision d_uniform_01
      double precision det(2)
      integer i
      integer inert(3)
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer seed
      double complex work(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  For a double complex hermitian matrix,'
      write ( *, '(a)' ) '  using packed storage,'
      write ( *, '(a)' ) '  ZHPFA factors the matrix.'
      write ( *, '(a)' ) '  ZHPDI computes the determinant, inverse,'
      write ( *, '(a)' ) '  or inertia.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n

      seed = 123456789
c
c  Set the values of the matrix A.
c
      k = 0

      do j = 1, n

        do i = 1, j-1
          k = k + 1
          a(k) = z_uniform_01 ( seed )
          a_save(i,j) = a(k)
          a_save(j,i) = conjg ( a(k) )
        end do

        k = k + 1
        a(k) = cmplx ( d_uniform_01 ( seed ), 0.0D+00 )
        a_save(j,j) = a(k)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zhpfa ( a, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZHPFA returned an error flag INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      job = 010
      call zhpdi ( a, n, ipvt, det, inert, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', det(2)
c
c  Get the inertia.
c
      job = 100
      call zhpdi ( a, n, ipvt, det, inert, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The inertia:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2x,i6)' ) inert(i)
      end do
c
c  Get the inverse.
c
      job = 001
      call zhpdi ( a, n, ipvt, det, inert, work, job )
c
c  Only the upper triangle is set, so the user must set up the
c  lower triangle:
c
      k = 0
      do j = 1, n
        do i = 1, j-1
          k = k + 1
          b(i,j) = a(k)
          b(j,i) = conjg ( a(k) )
        end do
        k = k + 1
        b(j,j) = a(k)
      end do

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + b(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests ZPBCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer m
      integer lda

      parameter ( n = 3 )
      parameter ( m = 1 )

      parameter ( lda = m+1 )

      double complex a(lda,n)
      integer i
      integer info
      double precision rcond
      double complex z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  For a double complex positive definite '
      write ( *, '(a)' ) '  hermitian band matrix,'
      write ( *, '(a)' )
     &  '  ZPBCO estimates the reciprocal condition number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the value of the superdiagonal and diagonal.
c
      a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00 )
      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00 )

      a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00 )
      a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00 )
      a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00 )
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call zpbco ( a, lda, n, m, rcond, z, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  ZPBCO returned INFO = ', info
        write ( *, '(a)' ) '  The factorization was not completed.'
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests ZPBDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer m
      integer lda

      parameter ( n = 3 )
      parameter ( m = 1 )

      parameter ( lda = m+1 )

      double complex a(lda,n)
      double precision det(2)
      integer i
      integer info

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  For a double complex positive definite '
      write ( *, '(a)' ) '  hermitian band matrix,'
      write ( *, '(a)' ) '  ZPBDI computes the determinant as'
      write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
      write ( *, '(a)' ) ' '
c
c  Set the value of the superdiagonal and diagonal.
c
      a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00 )
      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00 )

      a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00 )
      a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00 )
      a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00 )

      call zpbfa ( a, lda, n, m, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' ) '  Error!  ZPBFA returns INFO = ', info
        return
      end if

      call zpbdi ( a, lda, n, m, det )

      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', det(2)

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests ZPBFA and ZPBSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer m
      integer lda

      parameter ( n = 3 )
      parameter ( m = 1 )

      parameter ( lda = m+1 )

      double complex a(lda,n)
      double complex b(n)
      integer i
      integer info

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  For a double complex positive definite '
      write ( *, '(a)' ) '  hermitian band matrix,'
      write ( *, '(a)' ) '  ZPBFA computes the LU factors.'
      write ( *, '(a)' ) '  ZPBSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the value of the superdiagonal and diagonal.
c
      a(1,1) = cmplx ( 0.0000D+00,  0.0000D+00 )
      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(1,3) = cmplx ( 2.0905D+00,  1.1505D+00 )

      a(2,1) = cmplx ( 4.5281D+00,  0.0000D+00 )
      a(2,2) = cmplx ( 5.0371D+00,  0.0000D+00 )
      a(2,3) = cmplx ( 4.7638D+00,  0.0000D+00 )
c
c  Set the right hand side.
c
      b(1) = cmplx (  8.7963D+00, -0.4294D+00 )
      b(2) = cmplx ( 18.4798D+00,  3.6662D+00 )
      b(3) = cmplx ( 18.4724D+00, -2.3010D+00 )
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call zpbfa ( a, lda, n, m, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' ) '  Error!  ZPBFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call zpbsl ( a, lda, n, m, b )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution:'
      write ( *, '(a)' ) '  (Should be roughly (1,2,3)):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
      end do

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests ZPOCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(lda,n)
      integer i
      integer info
      double precision rcond
      double complex z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite matrix,'
      write ( *, '(a)' )
     &  '  ZPOCO estimates the reciprocal condition number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00 )
      a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00 )
      a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00 )

      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00 )
      a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00 )

      a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00 )
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call zpoco ( a, lda, n, rcond, z, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests ZPOFA and ZPODI.
c
c  Discussion:
c
c    ZPOFA factors a positive definite symmetric matrix,
c    and ZPODI can compute the determinant or the inverse.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(lda,n)
      double precision det(2)
      integer i
      integer info
      integer j
      integer job

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite matrix,'
      write ( *, '(a)' ) '  ZPOFA computes the LU factors,'
      write ( *, '(a)' ) '  ZPODI computes the inverse or determinant.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00 )
      a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00 )
      a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00 )

      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00 )
      a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00 )

      a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00 )
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call zpofa ( a, lda, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, ZPOFA returns INFO = ', info
        return
      end if
c
c  Get the determinant and inverse.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Get the determinant and inverse.'

      job = 11
      call zpodi ( a, lda, n, det, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant  = ', det(1), ' * 10 ** ', det(2)
c
c  ZPODI produces only the 'upper half triangle' of the inverse,
c  which is actually symmetric.  Thus, the lower half could be
c  produced by copying from the upper half.  However, the first row
c  of A, as returned, is exactly the first row of the inverse.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First row of inverse:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,6f10.4)' ) ( a(1,j), j = 1, n )

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests ZPOFA and ZPOSL.
c
c  Discussion:
c
c    ZPOFA factors a Hermitian positive definite matrix,
c    and ZPOSL can solve a factored linear system.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(lda,n)
      double complex b(n)
      integer i
      integer info
      integer j
      double complex x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite matrix,'
      write ( *, '(a)' ) '  ZPOFA computes the LU factors.'
      write ( *, '(a)' ) '  ZPOSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = cmplx ( 2.5281D+00,  0.0000D+00 )
      a(2,1) = cmplx ( 2.1341D+00,  0.2147D+00 )
      a(3,1) = cmplx ( 2.4187D+00, -0.2932D+00 )

      a(1,2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(2,2) = cmplx ( 3.0371D+00,  0.0000D+00 )
      a(3,2) = cmplx ( 2.0905D+00, -1.1505D+00 )

      a(1,3) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(2,3) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(3,3) = cmplx ( 2.7638D+00,  0.0000D+00 )
c
c  Set the right hand side.
c
      do i = 1, n
        x(i) = cmplx ( 2 * i - 1, 2 * i  )
      end do

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call zpofa ( a, lda, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, ZPOFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call zposl ( a, lda, n, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution:'
      write ( *, '(a)' ) '  (Should be (1+2i),(3+4i),(5+6i):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests ZPPCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      integer info
      double precision rcond
      double complex z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite packed matrix,'
      write ( *, '(a)' )
     &  '  ZPPCO estimates the reciprocal condition number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1) = cmplx ( 2.5281D+00,  0.0000D+00 )

      a(2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(3) = cmplx ( 3.0371D+00,  0.0000D+00 )

      a(4) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(5) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(6) = cmplx ( 2.7638D+00,  0.0000D+00 )
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition number.'

      call zppco ( a, n, rcond, z, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition number = ', rcond

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests ZPPFA and ZPPDI.
c
c  Discussion:
c
c    ZPPFA factors a Hermitian positive definite packed matrix,
c    and ZPPDI can compute the determinant or the inverse.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex b(n,n)
      double precision det(2)
      integer i
      integer info
      integer j
      integer job
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite packed matrix,'
      write ( *, '(a)' ) '  ZPPFA factors the matrix.'
      write ( *, '(a)' ) '  ZPPDI computes the inverse or determinant.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1) = cmplx ( 2.5281D+00,  0.0000D+00 )

      a(2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(3) = cmplx ( 3.0371D+00,  0.0000D+00 )

      a(4) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(5) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(6) = cmplx ( 2.7638D+00,  0.0000D+00 )
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call zppfa ( a, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, ZPPFA returns INFO = ', info
        return
      end if
c
c  Invert the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Get the determinant and inverse.'

      job = 11
      call zppdi ( a, n, det, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant  = ', det(1), ' * 10 ** ', det(2)
c
c  ZPPDI produces only the 'upper half triangle' of the inverse,
c  which is actually symmetric.  Thus, the lower half could be
c  produced by copying from the upper half.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          b(i,j) = a(k)
          b(j,i) = conjg ( a(k) )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Inverse:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,6f10.4)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests ZPPFA and ZPPSL.
c
c  Discussion:
c
c    ZPOFA factors a Hermitian positive definite packed matrix,
c    and ZPOSL can solve a factored linear system.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex b(n)
      integer i
      integer info
      integer j
      integer k
      double complex x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite packed matrix,'
      write ( *, '(a)' ) '  ZPPFA factors the matrix.'
      write ( *, '(a)' ) '  ZPPSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1) = cmplx ( 2.5281D+00,  0.0000D+00 )

      a(2) = cmplx ( 2.1341D+00, -0.2147D+00 )
      a(3) = cmplx ( 3.0371D+00,  0.0000D+00 )

      a(4) = cmplx ( 2.4187D+00,  0.2932D+00 )
      a(5) = cmplx ( 2.0905D+00,  1.1505D+00 )
      a(6) = cmplx ( 2.7638D+00,  0.0000D+00 )

      b(1) = cmplx ( 20.12350D+00, 28.92670D+00 )
      b(2) = cmplx ( 14.36550D+00, 34.92680D+00 )
      b(3) = cmplx ( 27.69760D+00, 26.03750D+00 )
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call zppfa ( a, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, ZPPFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call zppsl ( a, n, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution:'
      write ( *, '(a)' ) '  (Should be (1+2i),(3+4i),(5+6i):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests ZPTSL.
c
c  Discussion:
c
c    ZPTSL factors and solves a Hermitian positive definite
c    tridiagonal system.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex b(n)
      double complex d(n)
      double complex e(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  For a double complex'
      write ( *, '(a)' ) '  Hermitian positive definite '
      write ( *, '(a)' ) '  tridiagonal matrix,'
      write ( *, '(a)' ) '  ZPTSL factors and solves a linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the value of the superdiagonal and diagonal.
c
      e(1) = cmplx ( 2.1341D+00, -0.2147D+00 )
      e(2) = cmplx ( 2.0905D+00,  1.1505D+00 )
      e(3) = cmplx ( 0.0000D+00,  0.0000D+00 )

      d(1) = cmplx ( 4.5281D+00,  0.0000D+00 )
      d(2) = cmplx ( 5.0371D+00,  0.0000D+00 )
      d(3) = cmplx ( 4.7638D+00,  0.0000D+00 )
c
c  Set the right hand side.
c
      b(1) = cmplx (  8.7963D+00, -0.4294D+00 )
      b(2) = cmplx ( 18.4798D+00,  3.6662D+00 )
      b(3) = cmplx ( 18.4724D+00, -2.3010D+00 )
c
c  Factor and solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix and solve the system.'

      call zptsl ( n, d, e, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution:'
      write ( *, '(a)' ) '  (Should be roughly (1,2,3)):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,2g14.6)' ) i, b(i)
      end do

      return
      end
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests ZQRDC and ZQRSL.
c
c  Discussion:
c
c    ZQRDC and ZQRSL compute the QR factorization, and use it
c    to solve linear systems.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer p
      integer lda

      parameter ( n = 3 )
      parameter ( p = 3 )
      parameter ( lda = n )

      double complex a(lda,p)
      double complex b(lda,p)
      integer i
      integer info
      integer ipvt(p)
      integer j
      integer job
      integer k
      double complex q(n,n)
      double complex qraux(p)
      double complex qty(n)
      double complex qy(n)
      double complex r(n,p)
      double complex rsd(n)
      integer seed
      double complex work(p)
      double complex xb(n)
      double complex y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  For a double complex general matrix,'
      write ( *, '(a)' ) '  ZQRDC computes the QR decomposition of a '
      write ( *, '(a)' )
     &  '  matrix, but does not return Q and R explicitly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Show how Q and R can be recovered using ZQRSL.'
c
c  Set the values of the matrix A.
c
      seed = 123456789

      call zmat_uniform_01 ( n, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(6f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Decompose the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Decompose the matrix.'

      job = 0
      do i = 1, p
        ipvt(i) = 0
      end do

      call zqrdc ( a, lda, n, p, qraux, ipvt, work, job )
c
c  Print out what ZQRDC has stored in A...
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  The packed matrix A which describes Q and R:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(6f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  ...and in QRAUX.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The QRAUX vector, containing some'
      write ( *, '(a)' ) '  additional information defining Q:'
      write ( *, '(a)' ) ' '

      write ( *, '(2x,6f8.4)' ) ( qraux(i), i = 1, n )
c
c  Print out the resulting R factor.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The R factor:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        do j = 1, p

          if ( j .lt. i ) then
            r(i,j) = 0.0D+00
          else
            r(i,j) = a(i,j)
          end if

        end do

        write ( *, '(2x,6f8.4)' ) ( r(i,j), j = 1, p )

      end do
c
c  Call ZQRSL to extract the information about the Q matrix.
c  We do this, essentially, by asking ZQRSL to tell us the
c  value of Q*Y, where Y is a column of the identity matrix.
c
      job = 10000

      do i = 1, n
c
c  Set the vector Y.
c
        do j = 1, n
          y(j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do

        y(i) = cmplx ( 1.0D+00, 0.0D+00 )
c
c  Ask ZQRSL to tell us what Q*Y is.
c
        call zqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd, xb,
     &      job, info )

        if ( info .ne. 0 ) then
          write ( *, '(a,i8)' )
     &      '  Error!  ZQRSL returns INFO = ', info
          return
        end if
c
c  Copy QY into the appropriate column of Q.
c
        do j = 1, n
          q(j,i) = qy(j)
        end do

      end do
c
c  Now print out the Q matrix we have extracted.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Q factor:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f8.4)' ) ( q(i,j), j = 1, n )
      end do
c
c  Compute Q*R to verify that it equals A.
c
      do i = 1, n
        do j = 1, p
          b(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            b(i,j) = b(i,j) + q(i,k) * r(k,j)
          end do
        end do
      end do
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product Q * R:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,6f8.4)' ) ( b(i,j), j = 1, p )
      end do

      return
      end
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests ZSICO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double precision d_uniform_01
      integer i
      integer ipvt(n)
      integer j
      integer lda
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) '  For a double complex symmetric matrix:'
      write ( *, '(a)' ) '  ZSICO factors the matrix and estimates'
      write ( *, '(a)' ) '  the reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = z_uniform_01 ( seed )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zsico ( a, lda, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests ZSIFA and ZSISL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex b(n)
      double precision d_uniform_01
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer lda
      integer seed
      double complex x(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) '  For a double complex symmetric matrix:'
      write ( *, '(a)' ) '  ZSIFA factors the matrix.'
      write ( *, '(a)' ) '  ZSISL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = z_uniform_01 ( seed )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      do i = 1, n
        x(i) = z_uniform_01 ( seed )
      end do

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zsifa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZSIFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the system.
c
      call zsisl ( a, lda, n, ipvt, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests ZSIFA and ZSIDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a(n,n)
      double complex a_save(n,n)
      double complex c(n,n)
      double precision d_uniform_01
      double complex det(2)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer lda
      integer seed
      double complex work(n)
      double complex z_uniform_01

      lda = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) '  For a double complex symmetric matrix:'
      write ( *, '(a)' ) '  ZSIFA factors the matrix.'
      write ( *, '(a)' ) '  ZSIDI computes the determinant or inverse.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the matrix A.
c
      seed = 123456789

      do i = 1, n
        a(i,i) = z_uniform_01 ( seed )
        do j = i+1, n
          a(i,j) = z_uniform_01 ( seed )
          a(j,i) = a(i,j)
        end do
      end do

      do j = 1, n
        do i = 1, n
          a_save(i,j) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zsifa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZSIFA returned an error flag INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      job = 10
      call zsidi ( a, lda, n, ipvt, det, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', dble ( det(2) )
c
c  Get the inverse.
c
      job = 01
      call zsidi ( a, lda, n, ipvt, det, work, job )
c
c  Only the upper triangle is set, so the user must set up the
c  lower triangle:
c
      do i = 1, n
        do j = 1, i-1
          a(i,j) = a(j,i)
        end do
      end do

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + a(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests ZSPCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double precision d_uniform_01
      integer i
      integer ipvt(n)
      integer j
      integer k
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  For a double complex symmetric matrix'
      write ( *, '(a)' ) '  in packed storage,'
      write ( *, '(a)' ) '  ZSPCO factors the matrix and estimates'
      write ( *, '(a)' ) '  the reciprocal condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the packed matrix A.
c
      seed = 123456789

      k = 0

      do j = 1, n

        do i = 1, j-1
          k = k + 1
          a(k) = z_uniform_01 ( seed )
        end do

        k = k + 1
        a(k) = z_uniform_01 ( seed )

      end do
c
c  Copy the packed matrix into a "normal" matrix.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          a_save(i,j) = a(k)
        end do
      end do

      do j = 1, n
        do i = j+1, n
          a_save(i,j) = a_save(j,i)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zspco ( a, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests ZSPFA and ZSPSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double complex b(n)
      double precision d_uniform_01
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer seed
      double complex x(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' ) '  For a double complex symmetric matrix'
      write ( *, '(a)' ) '  in packed storage,'
      write ( *, '(a)' ) '  ZSPFA factors the matrix.'
      write ( *, '(a)' ) '  ZSPSL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the packed matrix A.
c
      seed = 123456789

      k = 0

      do j = 1, n

        do i = 1, j-1
          k = k + 1
          a(k) = z_uniform_01 ( seed )
        end do

        k = k + 1
        a(k) = z_uniform_01 ( seed )

      end do
c
c  Copy the packed matrix into a "normal" matrix.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          a_save(i,j) = a(k)
        end do
      end do

      do j = 1, n
        do i = j+1, n
          a_save(i,j) = a_save(j,i)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      do i = 1, n
        x(i) = z_uniform_01 ( seed )
      end do

      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a_save(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2f8.4)' ) b(i)
      end do
c
c  Factor the matrix A.
c
      call zspfa ( a, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZSPFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the system.
c
      call zspsl ( a, n, ipvt, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      subroutine test33 ( )

c*********************************************************************72
c
cc TEST33 tests ZSPFA and ZSPDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double complex a((n*(n+1))/2)
      double complex a_save(n,n)
      double complex b_save(n,n)
      double complex c(n,n)
      double precision d_uniform_01
      double complex det(2)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer seed
      double complex work(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST33'
      write ( *, '(a)' ) '  For a double complex symmetric matrix'
      write ( *, '(a)' ) '  in packed storage,'
      write ( *, '(a)' ) '  ZSPFA factors the matrix.'
      write ( *, '(a)' ) '  ZSPDI computes the determinant or inverse.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Set the values of the packed matrix A.
c
      seed = 123456789
      k = 0

      do j = 1, n

        do i = 1, j-1
          k= k + 1
          a(k) = z_uniform_01 ( seed )
        end do

        k = k + 1
        a(k) = z_uniform_01 ( seed )

      end do
c
c  Copy the packed matrix into a "normal" matrix.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          a_save(i,j) = a(k)
        end do
      end do

      do j = 1, n
        do i = j+1, n
          a_save(i,j) = a_save(j,i)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( a_save(i,j), j = 1, n )
      end do
c
c  Factor the matrix A.
c
      call zspfa ( a, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' )
     &    '  ZSPFA returned an error flag INFO = ', info
        return
      end if
c
c  Get the determinant.
c
      job = 10
      call zspdi ( a, n, ipvt, det, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6,a,g14.6)' )
     & '  Determinant = ', det(1), ' * 10** ', dble ( det(2) )
c
c  Get the inverse.
c
      job = 01
      call zspdi ( a, n, ipvt, det, work, job )
c
c  Copy the packed matrix into a "normal" matrix.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          b_save(i,j) = a(k)
        end do
      end do

      do j = 1, n
        do i = j+1, n
          b_save(i,j) = b_save(j,i)
        end do
      end do

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + b_save(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test34 ( )

c*********************************************************************72
c
cc TEST34 tests ZSVDC.
c
c  Discussion:
c
c    ZSVDC computes the singular value decomposition:
c
c      A = U * S * conjg-transpose ( V )
c
c  Modified:
c
c    03 May 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 4 )
      parameter ( n = 3 )

      double complex a(m,n)
      double complex b(m,n)
c     double complex e(max(m+1,n))
      double complex e(m+n)
      integer i
      integer info
      integer j
      integer k
      integer lda
      integer ldu
      integer ldv
      integer job
c     complex s(max(m+1,n))
      double complex s(m+n)
      integer seed
      double complex sigma(m,n)
      double complex u(m,m)
      double complex v(n,n)
      double complex work(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST34'
      write ( *, '(a)' )
     &  '  For an MxN matrix A in double complex general storage,'
      write ( *, '(a)' )
     &  '  ZSVDC computes the singular value decomposition:'
      write ( *, '(a)' ) '    A = U * S * V^H'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n
c
c  Set A.
c
      seed = 123456789

      call zmat_uniform_01 ( m, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Decompose the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Decompose the matrix.'

      job = 11
      lda = m
      ldu = m
      ldv = n

      call zsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning:'
        write ( *, '(a,i8)' ) '  ZSVDC returned nonzero INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Singular values:'
      write ( *, '(a)' ) ' '

      do i = 1, min ( m, n )
        write ( *, '(2x,i4,2x,2g14.6)' ) i, s(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,8f10.4)' ) ( u(i,j), j = 1, m )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
      write ( *, '(a)' ) ' '

      do i = 1,  n
        write ( *, '(2x,6f10.4)' ) ( v(i,j), j = 1, n )
      end do

      do j = 1, n
        do i = 1, m
          if ( i .eq. j ) then
            sigma(i,j) = s(i)
          else
            sigma(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          end if
        end do
      end do

      do i = 1, m
        do j = 1, n
          b(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            b(i,j) = b(i,j) + sigma(i,k) * conjg ( v(j,k) )
          end do
        end do
      end do

      do i = 1, m
        do j = 1, n
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, m
            a(i,j) = a(i,j) + u(i,k) * b(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product U * S * V^H (should equal A):'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,6f10.4)' ) ( a(i,j), j = 1, n )
      end do

      return
      end
      subroutine test345 ( )

c*********************************************************************72
c
cc TEST345 tests ZSVDC.
c
c  Discussion:
c
c    ZSVDC computes the singular value decomposition:
c
c      A = U * S * conjg-transpose ( V )
c
c  Modified:
c
c    03 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 4 )
      parameter ( n = 4 )

      double complex a(m,n)
      double complex b(m,n)
      double complex e(m+n)
      double complex eye
      integer i
      integer info
      integer j
      integer job
      integer k
      integer lda
      integer ldu
      integer ldv
      double complex one
      double complex s(m+n)
      double complex sigma(m,n)
      double complex u(m,m)
      double complex v(n,n)
      double complex work(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST345'
      write ( *, '(a)' )
     &  '  For an MxN matrix A in double complex general storage,'
      write ( *, '(a)' )
     &  '  ZSVDC computes the singular value decomposition:'
      write ( *, '(a)' ) '    A = U * S * V^H'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n
c
c  Set A.
c
      one = cmplx ( 1.0D+00, 0.0D+00 )
      eye = cmplx ( 0.0D+00, 1.0D+00 )

      a(1,1) =   one
      a(2,1) = - eye
      a(3,1) = - one
      a(4,1) =   eye

      a(1,2) =   one
      a(2,2) = - one
      a(3,2) = - one
      a(4,2) =   one

      a(1,3) =   one
      a(2,3) =   one
      a(3,3) =   one
      a(4,3) =   one

      a(1,4) =   one
      a(2,4) =   eye
      a(3,4) = - one
      a(4,4) = - eye

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,8f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Decompose the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Decompose the matrix.'

      job = 11
      lda = m
      ldu = m
      ldv = n

      call zsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning:'
        write ( *, '(a,i8)' ) '  ZSVDC returned nonzero INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Singular values:'
      write ( *, '(a)' ) ' '

      do i = 1, min ( m, n )
        write ( *, '(2x,i4,2x,2g14.6)' ) i, s(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,8f10.4)' ) ( u(i,j), j = 1, m )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
      write ( *, '(a)' ) ' '

      do i = 1,  n
        write ( *, '(2x,8f10.4)' ) ( v(i,j), j = 1, n )
      end do

      do j = 1, n
        do i = 1, m
          if ( i .eq. j ) then
            sigma(i,j) = s(i)
          else
            sigma(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          end if
        end do
      end do

      do i = 1, m
        do j = 1, n
          b(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            b(i,j) = b(i,j) + sigma(i,k) * conjg ( v(j,k) )
          end do
        end do
      end do

      do i = 1, m
        do j = 1, n
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, m
            a(i,j) = a(i,j) + u(i,k) * b(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product U * S * V^H (should equal A):'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,8f10.4)' ) ( a(i,j), j = 1, n )
      end do

      return
      end
      subroutine test35 ( )

c*********************************************************************72
c
cc TEST35 tests ZTRCO.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(n,n)
      integer i
      integer j
      integer job
      integer k
      double precision rcond
      integer seed
      double complex z(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST35'
      write ( *, '(a)' ) '  For a double complex triangular matrix,'
      write ( *, '(a)' ) '  ZTRCO estimates the condition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      seed = 123456789

      k = 0
      do i = 1, n
        do j = 1, i
          k = k + 1
          a(i,j) = z_uniform_01 ( seed )
        end do
        do j = i+1, n
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do
c
c  Get the condition of the lower triangular matrix.
c
      job = 0
      call ztrco ( a, lda, n, rcond, z, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition RCOND = ', rcond

      return
      end
      subroutine test36 ( )

c*********************************************************************72
c
cc TEST36 tests ZTRDI.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 3 )
      parameter ( lda = n )

      double complex a(n,n)
      double complex a_save(n,n)
      double complex c(n,n)
      double complex det(2)
      integer i
      integer info
      integer j
      integer job
      integer k
      integer seed
      double complex x(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST36'
      write ( *, '(a)' ) '  For a double complex triangular matrix,'
      write ( *, '(a)' ) '  ZTRDI computes the determinant or inverse.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      seed = 123456789

      k = 0
      do i = 1, n
        do j = 1, i
          k = k + 1
          a(i,j) = z_uniform_01 ( seed )
        end do
        do j = i+1, n
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      do i = 1, n
        do j = 1, n
          a_save(i,j) = a(i,j)
        end do
      end do
c
c  Get the determinant of the lower triangular matrix.
c
      job = 100
      call ztrdi ( a, lda, n, det, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6,a,g14.6)' )
     &  '  Determinant = ', det(1), ' * 10** ', dble ( det(2) )
c
c  Get the inverse of the lower triangular matrix.
c
      job = 010
      call ztrdi ( a, lda, n, det, job, info )

      do i = 1, n
        do j = 1, n
          c(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
          do k = 1, n
            c(i,j) = c(i,j) + a(i,k) * a_save(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product inv(A) * A is '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(10f8.4)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test37 ( )

c*********************************************************************72
c
cc TEST37 tests ZTRSL.
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 10 )
      parameter ( lda = n )

      double complex a(n,n)
      double complex b(n)
      integer i
      integer info
      integer j
      integer job
      integer k
      integer seed
      double complex x(n)
      double complex z_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST37'
      write ( *, '(a)' ) '  For a double complex triangular matrix,'
      write ( *, '(a)' ) '  ZTRSL solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      seed = 123456789

      k = 0
      do i = 1, n
        do j = 1, i
          k = k + 1
          a(i,j) = z_uniform_01 ( seed )
        end do
        do j = i+1, n
          a(i,j) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do
c
c  Set the desired solution
c
      do i = 1, n
        x(i) = cmplx ( i, 10 * i )
      end do
c
c  Compute the corresponding right hand side.
c
      do i = 1, n
        b(i) = cmplx ( 0.0D+00, 0.0D+00 )
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do
c
c  Solve the lower triangular system.
c
      job = 0
      call ztrsl ( a, lda, n, b, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed                     Exact'
      write ( *, '(a)' ) '  Solution                     Solution'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(4g14.6)' ) b(i), x(i)
      end do

      return
      end
      function d_uniform_01 ( seed )

c*********************************************************************72
c
cc D_UNIFORM_01 returns a unit double precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      d_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      D_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision D_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision d_uniform_01
      integer k
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      d_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
      function z_uniform_01 ( seed )

c*********************************************************************72
c
cc Z_UNIFORM_01 returns a unit double precision complex pseudorandom number.
c
c  Discussion:
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Modified:
c
c    27 February 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex Z_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      double precision r
      integer k
      integer seed
      double precision theta
      double complex z_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( dble ( seed ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

      z_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

      return
      end
      subroutine zmat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc ZMAT_UNIFORM_01 returns a unit double complex pseudorandom matrix.
c
c  Discussion:
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      double complex c(m,n)
      integer i
      integer j
      double precision r
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seed
      double precision theta

      do j = 1, n
        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if

          r = sqrt ( dble ( seed ) * 4.656612875D-10 )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if

          theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

          c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ) )

        end do

      end do

      return
      end
