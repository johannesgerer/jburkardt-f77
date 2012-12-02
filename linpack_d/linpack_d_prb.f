      program main

c*********************************************************************72
c
cc MAIN is the main program for LINPACK_D_PRB.
c
c  Discussion:
c
c    LINPACK_D_PRB calls the double precision real LINPACK test routines.
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
      write ( *, '(a)' ) 'LINPACK_D_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LINPACK_D library.'

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
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINPACK_D_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DCHDC.
c
c  Modified:
c
c    01 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 4 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      double precision work(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For double precision, general storage,'
      write ( *, '(a)' ) '  DCHDC computes the Cholesky decomposition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
        a(j,j) = 2.0D+00
      end do

      do i = 1, n-1
        a(i,i+1) = -1.0D+00
      end do

      do i = 2, n
        a(i-1,i) = -1.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
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

      call dchdc ( a, lda, n, work, ipvt, job, info )

      if ( info .ne. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  DCHDC returned INFO = ', info
        write ( *, '(a)' )
     &    '  This means the matrix is not positive definite.'
      end if
c
c  Zero out the lower diagonal.
c
      do i = 2, n
        do j = 1, i-1
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
      end do
c
c  Compute the Cholesky product.
c
      do j = 1, n
        do i = 1, n
          b(i,j) = 0.0D+00
          do k = 1, n
            b(i,j) = b(i,j) + a(k,i) * a(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product U'' * U: '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DCHEX.
c
c  Modified:
c
c    02 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer nz

      parameter ( n = 5 )
      parameter ( lda = n )
      parameter ( nz = 1 )

      double precision a(lda,n)
      double precision b(lda,n)
      double precision c(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer l
      double precision s(n)
      integer seed
      double precision work(n)
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For double precision, general storage,'
      write ( *, '(a)' ) '  DCHEX can shift columns in a Cholesky '
      write ( *, '(a)' ) '  factorization.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
        a(j,j) = 2.0D+00
      end do

      do i = 1, n-1
        a(i,i+1) = -1.0D+00
      end do

      do i = 2, n
        a(i-1,i) = -1.0D+00
      end do

      do i = 1, n
        z(i) = dble ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector Z:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6)' ) z(i)
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

      call dchdc ( a, lda, n, work, ipvt, job, info )

      if ( info .ne. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  DCHDC returned INFO = ', info
        write ( *, '(a)' ) '  The matrix is not positive definite.'
      end if
c
c  Zero out the lower diagonal.
c
      do i = 2, n
        do j = 1, i-1
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
      end do
c
c  Right circular shift columns L through K.
c
      k = 1
      l = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a,i6)' )
     &  '  Right circular shift columns K  = ', k,
     &  ' through L = ', l

      job = 1
      call dchex ( a, lda, n, k, l, z, n, nz, c, s, job )
c
c  Left circular shift columns K+1 through L.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a,i6)' )
     &  '  Left circular shift columns K+1 = ', k+1,
     &  ' through L = ', l

      job = 2
      call dchex ( a, lda, n, k+1, l, z, n, nz, c, s, job )
c
c  Print the factorization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted Cholesky factor U:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted vector Z:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6)' ) z(i)
      end do
c
c  Compute the Cholesky product.
c
      do j = 1, n
        do i = 1, n
          b(i,j) = 0.0D+00
          do k = 1, n
            b(i,j) = b(i,j) + a(k,i) * a(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The shifted product U'' * U: '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DCHUD.
c
c  Modified:
c
c    03 March 2006
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

      double precision b(p)
      double precision beta(p)
      double precision c(p)
      double precision ddot
      integer i
      integer info
      integer j
      integer job
      double precision r(ldr,p)
      double precision rho(nz)
      double precision row(p)
      double precision s(p)
      integer seed
      double precision x(p)
      double precision y(nz)
      double precision z(p,nz)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For double precision, general storage,'
      write ( *, '(a)' ) '  DCHUD updates a Cholesky decomposition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this example, we use DCHUD to solve a'
      write ( *, '(a)' ) '  least squares problem R * b = z.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of equations is P = ', p
c
c  Initialize.
c
      do j = 1, p
        do i = 1, p
          r(i,j) = 0.0D+00
        end do
      end do

      do j = 1, nz
        do i = 1, p
          z(i,j) = 0.0D+00
        end do
      end do

      do i = 1, p
        x(i) = dble ( i )
      end do
c
c  Use DCHUD to form R, Z and RHO by adding X and Y a row at a time.
c  X is a row of the least squares matrix and Y the right hand side.
c
      seed = 123456789

      do i = 1, p
        call dmat_uniform_01 ( 1, p, seed, row )
        y(1) = ddot ( p, row, 1, x, 1 )
        rho(1) = 0.0D+00
        call dchud ( r, ldr, p, row, z, p, nz, y, rho, c, s )
      end do
c
c  Generate the least squares solution, b = inverse ( R ) * Z.
c
      do j = 1, nz

        do i = 1, p
          b(i) = z(i,j)
        end do

        job = 01

        call dtrsl ( r, ldr, p, b, job, info )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Solution vector # ', j
        write ( *, '(a)' ) '  (Should be (1,2,3...,n))'
        write ( *, '(a)' ) ' '

        do i = 1, p
          if ( i .le. 5 .or. p-5 .lt. i ) then
            write ( *, '(2x,i6,2x,g14.6)' ) i, b(i)
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
cc TEST04 tests DGBCO.
c
c  Modified:
c
c    03 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer ml
      integer mu
      integer lda

      parameter ( n = 10 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n)
      integer i
      integer ipivot(n)
      integer j
      integer m
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For a banded matrix in general format,'
      write ( *, '(a)' )
     &  '  DGBCO estimates the reciprocal condition number.'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      m = ml + mu + 1
      write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', m

      do j = 1, n
        a(m-1,j) = -1.0D+00
        a(m,j) = 2.0D+00
        a(m+1,j) = -1.0D+00
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call dgbco ( a, lda, n, ml, mu, ipivot, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition = ', rcond

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests DGBFA and DGBSL.
c
c  Discussion:
c
c    The problem solved here is a larger version of this one:
c
c    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
c                (-1  2 -1  0  0)                          (0)
c                ( 0 -1  2 -1  0)                          (0)
c                ( 0  0 -1  2 -1)                          (0)
c                ( 0  0  0 -1  2)                          (1)
c
c
c    Solution is   (1)
c                  (1)
c                  (1)
c                  (1)
c                  (1)
c
c  Modified:
c
c    04 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Local Parameters:
c
c    N is the number of equations.
c
c    ML is the number of subdiagonals,
c    MU the number of superdiagonals.
c
c    LDA is the leading dimension of the array used to store the
c    matrix, which must be at least 2*ML+MU+1.
c
      implicit none

      integer n
      integer ml
      integer mu
      integer lda

      parameter ( n = 10 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipivot(n)
      integer j
      integer job
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For a banded matrix in general format,'
      write ( *, '(a)' ) '  DGBFA factors the matrix,'
      write ( *, '(a)' ) '  DGBSL solves a factored linear system.'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Set the right hand side B.
c
      b(1) = 1.0D+00
      do i = 2, n-1
        b(i) = 0.0D+00
      end do
      b(n) = 1.0D+00
c
c  Set the matrix A.
c
      m = ml + mu + 1
      write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', m

      do j = 1, n
        a(m-1,j) = -1.0D+00
        a(m,j) = 2.0D+00
        a(m+1,j) = -1.0D+00
      end do
c
c  Call DGBFA to factor the matrix A.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dgbfa ( a, lda, n, ml, mu, ipivot, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i6)' ) '  Error!  DGBFA returns INFO = ', info
        return
      end if
c
c  Call DGBSL to solve the linear system.  The solution
c  is returned in B, that is, it overwrites the right hand side.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      job = 0
      call dgbsl ( a, lda, n, ml, mu, ipivot, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution values:'
      write ( *, '(a)' ) '  (All should be 1):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i6,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests DGBFA and DGBDI.
c
c  Discussion:
c
c    Matrix A is ( 2 -1  0  0  0)
c                (-1  2 -1  0  0)
c                ( 0 -1  2 -1  0)
c                ( 0  0 -1  2 -1)
c                ( 0  0  0 -1  2)
c
c  Modified:
c
c    04 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Local Parameters:
c
c    N is the number of equations.
c
c    ML is the number of subdiagonals,
c    MU the number of superdiagonals.
c
c    LDA is the leading dimension of the array used to store the
c    matrix, which must be at least 2*ML+MU+1.
c
      implicit none

      integer n_max
      integer ml
      integer mu
      integer lda

      parameter ( n_max = 128 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n_max)
      double precision det(2)
      integer i
      integer info
      integer ipivot(n_max)
      integer j
      integer m
      integer n
      integer n_log

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For a banded matrix in general format,'
      write ( *, '(a)' ) '  DGBFA factors the matrix,'
      write ( *, '(a)' ) '  DGBDI computes the determinant as'
      write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Find the determinant of the -1,2,-1 matrix'
      write ( *, '(a)' ) '  for N = 2, 4, 8, 16, 32, 64, 128.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  (For this matrix, det ( A ) = N + 1.)'
c
c  Set the matrix A.
c
      m = ml + mu + 1
      write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', m
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N    Mantissa       Exponent'
      write ( *, '(a)' ) ' '

      n = 1

      do n_log = 1, 7

        n = 2 * n

        do j = 1, n
          do i = 1, lda
            a(i,j) = 0.0D+00
          end do
        end do

        do j = 1, n
          a(m-1,j) = -1.0D+00
          a(m,j) = 2.0D+00
          a(m+1,j) = -1.0D+00
        end do

        call dgbfa ( a, lda, n, ml, mu, ipivot, info )

        if ( info .ne. 0 ) then
          write ( *, '(a,i6)' ) '  Error!  DGBFA returns INFO = ', info
          return
        end if

        call dgbdi ( a, lda, n, ml, mu, ipivot, det )

        write ( *, '(2x,i6,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests DGBFA and DGBSL.
c
c  Discussion:
c
c    DGBFA and DGBSL are for general banded matrices.
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer ml
      integer mu
      integer lda

      parameter ( n = 100 )
      parameter ( ml = 25 )
      parameter ( mu = 25 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer info
      integer ipivot(n)
      integer j
      integer job
      integer m
      double precision temp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  For a banded matrix in general format,'
      write ( *, '(a)' ) '  DGBFA factors the matrix,'
      write ( *, '(a)' ) '  DGBSL solves a factored linear system.'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Assign values to matrix A and right hand side B.
c
c  We want to try a problem with a significant bandwidth.
c
      m = ml + mu + 1
      write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', m

      do j = 1, n

        ilo = max ( 1, j - mu )
        ihi = min ( n, j + ml )

        temp = 0.0D+00
        do i = ilo, ihi
          a(i-j+m,j) = -1.0D+00
          temp = temp - 1.0D+00
        end do

        temp = temp + 1.0D+00
        a(m,j) = 4.0D+00 - temp
        b(j) = 4.0D+00

      end do
c
c  Factor the matrix A.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dgbfa ( a, lda, n, ml, mu, ipivot, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i6)' ) '  Error!  DGBFA returns INFO = ', info
        return
      end if
c
c  Call DGBSL to solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      job = 0
      call dgbsl ( a, lda, n, ml, mu, ipivot, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (All should be 1):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i6,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 calls DGECO and DGESL.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Local Parameters:
c
c    LDA defines the maximum matrix size we will use.
c
      implicit none

      integer lda
      parameter ( lda = 10 )

      double precision a(lda,lda)
      double precision b(lda)
      integer i
      integer ipvt(lda)
      integer job
      integer n
      double precision rcond
      double precision z(lda)

      n = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  DGECO factors a general matrix and computes'
      write ( *, '(a)' ) '  its reciprocal condition number;'
      write ( *, '(a)' ) '  DGESL solves a factored linear system.'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = 1.0D+00
      a(1,2) = 2.0D+00
      a(1,3) = 3.0D+00

      a(2,1) = 4.0D+00
      a(2,2) = 5.0D+00
      a(2,3) = 6.0D+00

      a(3,1) = 7.0D+00
      a(3,2) = 8.0D+00
      a(3,3) = 0.0D+00
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dgeco ( a, lda, n, ipvt, rcond, z )

      write ( *, '(a,g14.6)' )
     &  '  The reciprocal matrix condition number = ', rcond

      if ( rcond + 1.0D+00 .eq. 1.0D+00 ) then
        write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
        return
      end if
c
c  Set a right hand side.
c
      b(1) = 14.0D+00
      b(2) = 32.0D+00
      b(3) = 23.0D+00
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      job = 0
      call dgesl ( a, lda, n, ipvt, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution returned by DGESL'
      write ( *, '(a)' ) '  (Should be (1,2,3))'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do
c
c  A second right hand side can be solved without refactoring a.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call DGESL for a new right hand '
      write ( *, '(a)' ) '  side for the same, factored matrix.'
c
c  Set the right hand side.
c
      b(1) = 1.0D+00
      b(2) = 4.0D+00
      b(3) = 7.0D+00
c
c  Solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve a linear system.'

      job = 0
      call dgesl ( a, lda, n, ipvt, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution returned by DGESL'
      write ( *, '(a)' ) '  (should be (1,0,0))'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do
c
c  The transposed problem A'*X = B can be solved by DGESL
c  as well, without any refactoring.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call DGESL for transposed problem.'
c
c  Set the right hand side.
c
      b(1) =  6.0D+00
      b(2) =  6.0D+00
      b(3) = -3.0D+00
c
c  Solve the transposed system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Call DGESL to solve a transposed linear system.'

      job = 1
      call dgesl ( a, lda, n, ipvt, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution returned by DGESL'
      write ( *, '(a)' ) '  (should be (-1,0,1))'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests DGEFA and DGEDI.
c
c  Modified:
c
c    06 March 2006
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

      double precision a(lda,n)
      double precision det(2)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job
      double precision work(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  DGEFA factors a general matrix;'
      write ( *, '(a)' ) '  DGEDI computes the inverse and determinant'
      write ( *, '(a)' ) '  of a factored matrix.'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = 1.0D+00
      a(1,2) = 2.0D+00
      a(1,3) = 3.0D+00

      a(2,1) = 4.0D+00
      a(2,2) = 5.0D+00
      a(2,3) = 6.0D+00

      a(3,1) = 7.0D+00
      a(3,2) = 8.0D+00
      a(3,3) = 0.0D+00
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix'

      call dgefa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) '  Error!  The matrix is nearly singular!'
        return
      end if
c
c  Get the inverse and determinant.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Get the inverse and determinant'

      job = 11
      call dgedi ( a, lda, n, ipvt, det, work, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  The determinant = ', det(1), ' * 10 ** ', det(2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The inverse matrix:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests DGEFA and DGESL.
c
c  Discussion:
c
c    Solve A*x = b where A is a given matrix, and B a right hand side.
c
c    We will also assume that A is stored in the simplest
c    possible way.
c
c  Modified:
c
c    07 March 2006
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

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  DGEFA factors a general matrix;'
      write ( *, '(a)' ) '  DGESL solves a factored linear system;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = 1.0D+00
      a(1,2) = 2.0D+00
      a(1,3) = 3.0D+00

      a(2,1) = 4.0D+00
      a(2,2) = 5.0D+00
      a(2,3) = 6.0D+00

      a(3,1) = 7.0D+00
      a(3,2) = 8.0D+00
      a(3,3) = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      b(1) = 14.0D+00
      b(2) = 32.0D+00
      b(3) = 23.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The right hand side B is '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix'

      call dgefa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i6)' )
     &    '  DGEFA returned an error flag INFO = ', info
        return
      end if
c
c  If no error occurred, now use DGESL to solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      job = 0
      call dgesl ( a, lda, n, ipvt, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGESL returns the solution:'
      write ( *, '(a)' ) '  (Should be (1,2,3))'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests DGEFA and DGESL.
c
c  Discussion:
c
c    In this example, we solve a relatively large linear system.
c
c  Modified:
c
c    07 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 100 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer job

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  DGEFA factors a general matrix;'
      write ( *, '(a)' ) '  DGESL solves a factored linear system;'
      write ( *, '(a,i6)' ) '  The matrix size is N = ', n
c
c  Assign values to the matrix A and the right hand side B.
c
c  The problem is just an enlarged version of the
c  problem for N = 5, which is:
c
c  Matrix A is ( n -1 -1 -1 -1)    Right hand side B is  (1)
c              (-1  n -1 -1 -1)                          (1)
c              (-1 -1  n -1 -1)                          (1)
c              (-1 -1 -1  n -1)                          (1)
c              (-1 -1 -1 -1  n)                          (1)
c
c  Solution is   (1)
c                (1)
c                (1)
c                (1)
c                (1)
c
      do i = 1, n
        b(i) = 1.0D+00
      end do

      do j = 1, n
        do i = 1, n
          a(i,j) = -1.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = dble ( n )
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dgefa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i6)' )
     &    '  DGEFA returned an error flag INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      job = 0
      call dgesl ( a, lda, n, ipvt, b, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last five solution entries:'
      write ( *, '(a)' ) '  (All of them should be 1.)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i6,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests DGTSL.
c
c  Modified:
c
c    08 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 100 )

      double precision b(n)
      double precision c(n)
      double precision d(n)
      double precision e(n)
      integer i
      integer info

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  For a general tridiagonal matrix,'
      write ( *, '(a)' ) '  DGTSL factors and solves a linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
      write ( *, '(a)' ) ' '
c
c  Set up the linear system, by storing the values of the
c  subdiagonal, diagonal, and superdiagonal in C, D, and E,
c  and the right hand side in B.
c
      c(1) = 0.0D+00
      do i = 2, n
        c(i) = -1.0D+00
      end do

      do i = 1, n
        d(i) = 2.0D+00
      end do

      do i = 1, n-1
        e(i) = -1.0D+00
      end do
      e(n) = 0.0D+00

      do i = 1, n-1
        b(i) = 0.0D+00
      end do
      b(n) = dble ( n + 1 )
c
c  Factor and solve the system in one step.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix and solve the system.'

      call dgtsl ( n, c, d, e, b, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGTSL returns nonzero INFO = ', info
        return
      end if
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 entries of solution:'
      write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests DPBCO.
c
c  Modified:
c
c    08 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 10 )
      parameter ( lda = 2 )

      double precision a(lda,n)
      integer i
      integer info
      integer j
      integer m
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  For a positive definite symmetric '
      write ( *, '(a)' ) '  band matrix, DPBCO estimates the '
      write ( *, '(a)' ) '  reciprocal condition number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the number of nonzero diagonals.
c
      m = 1
c
c  Set the value of the superdiagonal and diagonal.
c
      a(1,1)   =  0.0D+00
      do j = 2, n
        a(1,j) = -1.0D+00
      end do
      do j = 1, n
        a(2,j) =  2.0D+00
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call dpbco ( a, lda, n, m, rcond, z, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests DPBDI.
c
c  Modified:
c
c    09 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      integer lda

      parameter ( n_max = 128 )
      parameter ( lda = 2 )

      double precision a(lda,n_max)
      double precision det(2)
      integer i
      integer info
      integer j
      integer m
      integer n
      integer n_log

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  For a positive definite symmetric band'
      write ( *, '(a)' ) '  matrix, DPBDI computes the determinant as'
      write ( *, '(a)' ) '    det = MANTISSA * 10**EXPONENT'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Find the determinant of the -1,2,-1 matrix'
      write ( *, '(a)' ) '  for N = 2, 4, 8, 16, 32, 64, 128.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  (For this matrix, det ( A ) = N + 1.)'
      write ( *, '(a)' ) ' '
c
c  Set the number of nonzero diagonals.
c
      m = 1

      write ( *, '(a,i8)' )
     &  '  The bandwidth of the matrix is ', 2 * m + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N    Mantissa       Exponent'
      write ( *, '(a)' ) ' '

      n = 1

      do n_log = 1, 7

        n = 2 * n

        do j = 1, n
          do i = 1, lda
            a(i,j) = 0.0D+00
          end do
        end do

        a(1,1) =  0.0D+00

        do j = 2, n
          a(1,j) = -1.0D+00
        end do

        do j = 1, n
          a(2,j) =  2.0D+00
        end do

        call dpbfa ( a, lda, n, m, info )

        if ( info .ne. 0 ) then
          write ( *, '(a,i8)' )
     &      '  Error!  DGBFA returns INFO = ', info
          return
        end if

        call dpbdi ( a, lda, n, m, det )

        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, det(1), det(2)

      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests DPBFA and DPBSL.
c
c  Discussion:
c
c    DPBFA and DPBSL are for a positive definite symmetric band matrix.
c
c    The problem is just an enlarged version of the
c    problem for N = 5, which is:
c
c    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
c                (-1  2 -1  0  0)                          (0)
c                ( 0 -1  2 -1  0)                          (0)
c                ( 0  0 -1  2 -1)                          (0)
c                ( 0  0  0 -1  2)                          (1)
c
c    solution is   (1)
c                  (1)
c                  (1)
c                  (1)
c                  (1)
c  Modified:
c
c    09 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 10 )
      parameter ( lda = 2 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer j
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  For a positive definite symmetric band'
      write ( *, '(a)' ) '  matrix,'
      write ( *, '(a)' ) '  DPBFA computes the LU factors.'
      write ( *, '(a)' ) '  DPBSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Assign values to matrix A and right hand side B.
c
c  Set the right hand side.
c
      b(1) = 1.0D+00
      do i = 2, n-1
        b(i) = 0.0D+00
      end do
      b(n) = 1.0D+00
c
c  Set the number of nonzero diagonals.
c
      m = 1
c
c  Set the value of the superdiagonal and diagonal.
c
      a(1,1) =  0.0D+00

      do j = 2, n
        a(1,j) = -1.0D+00
      end do

      do j = 1, n
        a(2,j) =  2.0D+00
       end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dpbfa ( a, lda, n, m, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' )
     &    '  Error!  DPBFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call dpbsl ( a, lda, n, m, b )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last solution entries:'
      write ( *, '(a)' ) '  (All should be 1):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests DPOCO.
c
c  Modified:
c
c    10 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      integer i
      integer info
      integer j
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
      write ( *, '(a)' ) '  DPOCO estimates the reciprocal condition'
      write ( *, '(a)' ) '  number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = 2.0D+00
        if ( 1 .lt. i ) then
          a(i,i-1) = -1.0D+00
        end if
        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call dpoco ( a, lda, n, rcond, z, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition  = ', rcond

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests DPOFA and DPODI.
c
c  Discussion:
c
c    DPOFA factors a positive definite symmetric matrix,
c    and DPODI can compute the determinant or the inverse.
c
c  Modified:
c
c    10 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision det(2)
      integer i
      integer info
      integer j
      integer job

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
      write ( *, '(a)' ) '  DPOFA computes the LU factors,'
      write ( *, '(a)' ) '  DPODI computes the inverse or determinant.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = 2.0D+00
        if ( 1 .lt. i ) then
          a(i,i-1) = -1.0D+00
        end if
        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dpofa ( a, lda, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, DPOFA returns INFO = ', info
        return
      end if
c
c  Get the determinant and inverse.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Get the determinant and inverse.'

      job = 11
      call dpodi ( a, lda, n, det, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant  = ', det(1), ' * 10 ** ', det(2)
c
c  DPODI produces only the 'upper half triangle' of the inverse,
c  which is actually symmetric.  Thus, the lower half could be
c  produced by copying from the upper half.  However, the first row
c  of A, as returned, is exactly the first row of the inverse.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First row of inverse:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( a(1,j), j = 1, n )

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests DPOFA and DPOSL.
c
c  Discussion:
c
c    DPOFA factors a positive definite symmetric matrix,
c    and DPOSL can solve a factored linear system.
c
c  Modified:
c
c    11 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 20 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer j
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
      write ( *, '(a)' ) '  DPOFA computes the LU factors.'
      write ( *, '(a)' ) '  DPOSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = 2.0D+00
        if ( 1 .lt. i ) then
          a(i,i-1) = -1.0D+00
        end if
        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if
      end do
c
c  Set the right hand side.
c
      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dpofa ( a, lda, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, DPOFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call dposl ( a, lda, n, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .le. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests DPPCO.
c
c  Modified:
c
c    11 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a((n*(n+1))/2)
      integer i
      integer info
      integer j
      integer k
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  For a positive definite symmetric '
      write ( *, '(a)' ) '  packed matrix,'
      write ( *, '(a)' ) '  DPPCO estimates the reciprocal condition '
      write ( *, '(a)' ) '  number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          if ( i .eq. j - 1 ) then
            a(k) = -1.0D+00
          else if ( i .eq. j ) then
            a(k) = 2.0D+00
          else
            a(k) = 0.0D+00
          end if
        end do
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition number.'

      call dppco ( a, n, rcond, z, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Reciprocal condition number = ', rcond

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests DPPFA and DPPDI.
c
c  Discussion:
c
c    DPPFA factors a packed positive definite symmetric matrix,
c    and DPPDI can compute the determinant or the inverse.
c
c  Modified:
c
c    11 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a((n*(n+1))/2)
      double precision b(n,n)
      double precision det(2)
      integer i
      integer info
      integer j
      integer job
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  For a positive definite symmetric '
      write ( *, '(a)' ) '  packed matrix,'
      write ( *, '(a)' ) '  DPPFA factors the matrix.'
      write ( *, '(a)' ) '  DPPDI computes the inverse or determinant.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          if ( i .eq. j - 1 ) then
            a(k) = -1.0D+00
          else if ( i .eq. j ) then
            a(k) = 2.0D+00
          else
            a(k) = 0.0D+00
          end if
        end do
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dppfa ( a, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, DPPFA returns INFO = ', info
        return
      end if
c
c  Invert the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Get the determinant and inverse.'

      job = 11
      call dppdi ( a, n, det, job )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  Determinant  = ', det(1), ' * 10 ** ', det(2)
c
c  DPPDI produces only the 'upper half triangle' of the inverse,
c  which is actually symmetric.  Thus, the lower half could be
c  produced by copying from the upper half.  However, the first row
c  of A, as returned, is exactly the first row of the inverse.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          b(i,j) = a(k)
          b(j,i) = a(k)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Inverse:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( b(i,j), j = 1, n )
      end do

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests DPPFA and DPPSL.
c
c  Discussion:
c
c    DPOFA factors a positive definite symmetric matrix,
c    and DPOSL can solve a factored linear system.
c
c  Modified:
c
c    11 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      double precision a((n*(n+1))/2)
      double precision b(n)
      integer i
      integer info
      integer j
      integer k
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  For a positive definite symmetric '
      write ( *, '(a)' ) '  packed matrix,'
      write ( *, '(a)' ) '  DPPFA factors the matrix.'
      write ( *, '(a)' ) '  DPPSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        b(i) = 0.0D+00
      end do

      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          if ( i .eq. j - 1 ) then
            a(k) = -1.0D+00
            b(i) = b(i) + a(k) * x(j)
            b(j) = b(j) + a(k) * x(i)
          else if ( i .eq. j ) then
            a(k) = 2.0D+00
            b(i) = b(i) + a(k) * x(i)
          else
            a(k) = 0.0D+00
          end if
        end do
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dppfa ( a, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Error, DPPFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call dppsl ( a, n, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests DPTSL.
c
c  Discussion:
c
c    DPTSL factors and solves a positive definite symmetric tridiagonal system.
c
c  Modified:
c
c    11 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      double precision b(n)
      double precision d(n)
      double precision e(n)
      integer i
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  For a positive definite symmetric'
      write ( *, '(a)' ) '  tridiagonal matrix,'
      write ( *, '(a)' ) '  DPTSL factors and solves a linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Set the matrix A.
c
      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        b(i) = 0.0D+00
      end do
      do i = 1, n
        d(i) = 2.0D+00
      end do
      do i = 1, n-1
        e(i) = -1.0D+00
      end do
      e(n) = 0.0D+00

      do i = 1, n

        if ( 1 .lt. i ) then
          b(i) = b(i) + e(i-1) * x(i-1)
        end if
        b(i) = b(i) + d(i) * x(i)
        if ( i .lt. n ) then
          b(i) = b(i) + e(i) * x(i+1)
        end if

      end do
c
c  Factor and solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix and solve the system.'

      call dptsl ( n, d, e, b )
c
c  Print the result.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (Should be 1,2,3,4,5,...,n-1,n):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests DQRDC and DQRSL.
c
c  Discussion:
c
c    DQRDC and DQRSL compute the QR factorization, and use it
c    to solve linear systems.
c
c  Modified:
c
c    12 March 2006
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

      double precision a(lda,p)
      double precision b(lda,p)
      integer i
      integer info
      integer ipvt(p)
      integer j
      integer job
      integer k
      double precision q(n,n)
      double precision qraux(p)
      double precision qty(n)
      double precision qy(n)
      double precision r(n,p)
      double precision rsd(n)
      double precision work(p)
      double precision xb(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  For a general matrix,'
      write ( *, '(a)' ) '  DQRDC computes the QR decomposition of a '
      write ( *, '(a)' ) '  matrix, but does not return Q and R'
      write ( *, '(a)' ) '  explicitly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recover Q and R using DQRSL.'
c
c  Set the matrix A.
c
      a(1,1) = 1.0D+00
      a(2,1) = 1.0D+00
      a(3,1) = 0.0D+00

      a(1,2) = 1.0D+00
      a(2,2) = 0.0D+00
      a(3,2) = 1.0D+00

      a(1,3) = 0.0D+00
      a(2,3) = 1.0D+00
      a(3,3) = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The original matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, p )
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

      call dqrdc ( a, lda, n, p, qraux, ipvt, work, job )
c
c  Print out what DQRDC has stored in A...
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The packed matrix A,'
      write ( *, '(a)' ) '  describing Q and R:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, p )
      end do
c
c  ...and in QRAUX.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The QRAUX vector, containing additional'
      write ( *, '(a)' ) '  information defining Q:'
      write ( *, '(a)' ) ' '

      write ( *, '(2x,5g14.6)' ) ( qraux(i), i = 1, n )
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

        write ( *, '(2x,5g14.6)' ) ( r(i,j), j = 1, p )

      end do
c
c  Call DQRSL to extract the information about the Q matrix.
c  We do this, essentially, by asking DQRSL to tell us the
c  value of Q*Y, where Y is a column of the identity matrix.
c
      job = 10000

      do j = 1, n
c
c  Set the vector Y.
c
        do i = 1, n
          y(i) = 0.0D+00
        end do

        y(j) = 1.0D+00
c
c  Ask DQRSL to tell us what Q*Y is.
c
        call dqrsl ( a, lda, n, p, qraux, y, qy, qty, b, rsd,
     &    xb, job, info )

        if ( info .ne. 0 ) then
          write ( *, '(a,i8)' ) '  Error!  DQRSL returns INFO = ', info
          return
        end if
c
c  Copy QY into the appropriate column of Q.
c
        do i = 1, n
          q(i,j) = qy(i)
        end do

      end do
c
c  Now print out the Q matrix we have extracted.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Q factor:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( q(i,j), j = 1, n )
      end do
c
c  Compute Q*R to verify that it equals A.
c
      do i = 1, n
       do j = 1, p
          b(i,j) = 0.0D+00
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
        write ( *, '(2x,5g14.6)' ) ( b(i,j), j = 1, p )
      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests DSICO.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 100 )
      parameter ( lda = n )

      double precision a(lda,n)
      integer i
      integer ipvt(n)
      integer j
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  For a symmetric indefinite matrix,'
      write ( *, '(a)' ) '  DSICO estimates the reciprocal condition'
      write ( *, '(a)' ) '  number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Assign values to the matrix A.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = 2.0D+00
        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call dsico ( a, lda, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition = ', rcond

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests DSIFA and DSISL.
c
c  Discussion:
c
c    DSIFA and DSISL are for symmetric indefinite matrices.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 100 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipvt(n)
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  For a symmetric indefinite matrix,'
      write ( *, '(a)' ) '  DSIFA factors the matrix,'
      write ( *, '(a)' ) '  DSISL solves a factored linear system,'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Assign values to the matrix A and the right hand side B.
c
      do i = 1, n-1
        b(i) = 0.0D+00
      end do
      b(n)= dble ( n + 1 )

      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,i) = 2.0D+00
        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dsifa ( a, lda, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' ) '  Error!  DSIFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call dsisl ( a, lda, n, ipvt, b )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests DSPCO.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 100 )

      double precision a((n*(n+1))/2)
      integer i
      integer ipvt(n)
      integer j
      integer k
      double precision rcond
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  For a symmetric indefinite packed matrix,'
      write ( *, '(a)' ) '  DSPCO estimates the reciprocal condition'
      write ( *, '(a)' ) '  number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Assign values to the matrix A.
c
      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          if ( i .eq. j ) then
            a(k) = 2.0D+00
          else if ( j .eq. i+1 ) then
            a(k) = -1.0D+00
          else
            a(k) = 0.0D+00
          end if
        end do
      end do
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition.'

      call dspco ( a, n, ipvt, rcond, z )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Estimated reciprocal condition = ', rcond

      return
      end
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests DSPFA and DSPSL.
c
c  Discussion:
c
c    DSPFA and DSPSL are for packed symmetric indefinite matrices.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 100 )

      double precision a((n*(n+1))/2)
      double precision b(n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  For a symmetric indefinite packed matrix,'
      write ( *, '(a)' ) '  DSPFA factors the matrix,'
      write ( *, '(a)' ) '  DSPSL solves a factored linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Assign values to the matrix A and the right hand side B.
c
      do i = 1, n-1
        b(i) = 0.0D+00
      end do
      b(n)= dble ( n + 1 )

      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          if ( i .eq. j ) then
            a(k) = 2.0D+00
          else if ( j .eq. i+1 ) then
            a(k) = -1.0D+00
          else
            a(k) = 0.0D+00
          end if
        end do
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Factor the matrix.'

      call dspfa ( a, n, ipvt, info )

      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' ) '  Error!  DSPFA returns INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call dspsl ( a, n, ipvt, b )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The first and last 5 solution entries:'
      write ( *, '(a)' ) '  (Should be (1,2,3,4,5,...,n-1,n))'
      write ( *, '(a)' ) ' '

      do i = 1, n
        if ( i .le. 5 .or. n-5 .lt. i ) then
          write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
        end if
        if ( i .eq. 5 ) then
          write ( *, '(a)' ) '  ......  ..............'
        end if
      end do

      return
      end
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests DSVDC.
c
c  Discussion:
c
c    SSVDC computes the singular value decomposition:
c
c      A = U * S * V'
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

      parameter ( m = 6 )
      parameter ( n = 4 )

      double precision a(m,n)
c     double precision e(max(m+1,n))
      double precision e(m+n)
      integer i
      integer info
      integer j
      integer job
      integer k
      integer lda
      integer ldu
      integer ldv
c     double precision s(max(m+1,n))
      double precision s(m+n)
      integer seed
      double precision sigma(m,n)
      double precision sv(m,n)
      double precision u(m,m)
      double precision v(n,n)
      double precision work(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) '  For an MxN matrix A in general storage,'
      write ( *, '(a)' ) '  DSVDC computes the singular value '
      write ( *, '(a)' ) '  decomposition:'
      write ( *, '(a)' ) '    A = U * S * V'''
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n
c
c  Set A.
c
      seed = 123456789

      call dmat_uniform_01 ( m, n, seed, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,7f10.4)' ) ( a(i,j), j = 1, n )
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

      call dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning:'
        write ( *, '(a,i8)' ) '  DSVDC returned nonzero INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Singular values:'
      write ( *, '(a)' ) ' '

      do i = 1, min ( m, n )
        write ( *, '(2x,i4,2x,g14.6)' ) i, s(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Left Singular Vector Matrix U:'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,7f10.4)' ) ( u(i,j), j = 1, m )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Right Singular Vector Matrix V:'
      write ( *, '(a)' ) ' '

      do i = 1,  n
        write ( *, '(2x,7f10.4)' ) ( v(i,j), j = 1, n )
      end do

      do i = 1, m
        do j = 1, n
          sigma(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        sigma(i,i) = s(i)
      end do

      do i = 1, m
        do j = 1, n
          sv(i,j) = 0.0D+00
          do k = 1, n
            sv(i,j) = sv(i,j) + sigma(i,k) * v(j,k)
          end do
        end do
      end do

      do i = 1, m
        do j = 1, n
          a(i,j) = 0.0D+00
          do k = 1, n
            a(i,j) = a(i,j) + u(i,k) * sv(k,j)
          end do
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The product U * S * V'' (should equal A):'
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,7f10.4)' ) ( a(i,j), j = 1, n )
      end do

      return
      end
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests DTRCO.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      integer i
      integer j
      integer job
      double precision rcond
      integer seed
      double precision z(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) '  For a triangular matrix,'
      write ( *, '(a)' ) '  DTRCO computes the LU factors and'
      write ( *, '(a)' ) '  computes its reciprocal condition number.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Lower triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = i+1, n
          a(i,j) = 0.0D+00
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lower triangular matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do

      job = 0
c
c  Estimate the condition.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition:'

      call dtrco ( a, lda, n, rcond, z, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  The reciprocal condition number = ', rcond
c
c  Upper triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = 1, i - 1
          a(i,j) = 0.0D+00
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Upper triangular matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do
c
c  Estimate the condition.
c
      job = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the condition:'

      call dtrco ( a, lda, n, rcond, z, job )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  The reciprocal condition number = ', rcond

      return
      end
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests DTRDI.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision det(2)
      integer i
      integer info
      integer j
      integer job
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) '  For a triangular matrix,'
      write ( *, '(a)' ) '  DTRDI computes the determinant or inverse.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Lower triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = i+1, n
          a(i,j) = 0.0D+00
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lower triangular matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do

      job = 110

      call dtrdi ( a, lda, n, det, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  The determinant = ', det(1), ' * 10 ** ', det(2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The inverse matrix:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do
c
c  Upper triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = 1, i - 1
          a(i,j) = 0.0D+00
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Upper triangular matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do

      job = 111

      call dtrdi ( a, lda, n, det, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' )
     &  '  The determinant = ', det(1), ' * 10 ** ', det(2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The inverse matrix:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5g14.6)') ( a(i,j), j = 1, n )
      end do

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests DTRSL.
c
c  Discussion:
c
c    DTRSL solves triangular linear systems.
c
c  Modified:
c
c    12 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer j
      integer job
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  For a triangular matrix,'
      write ( *, '(a)' ) '  DTRSL solves a linear system.'
      write ( *, '(a,i8)' ) '  The matrix size is N = ', n
c
c  Lower triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = i+1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For a lower triangular matrix A,'
      write ( *, '(a)' ) '  solve A * x = b'

      job = 00

      call dtrsl ( a, lda, n, b, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
      end do

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(j,i) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For a lower triangular matrix A,'
      write ( *, '(a)' ) '  solve A'' * x = b'

      job = 10

      call dtrsl ( a, lda, n, b, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
      end do
c
c  Upper triangular matrix A.
c
      call dmat_uniform_01 ( n, n, seed, a )

      do i = 1, n
        do j = 1, i - 1
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For an upper triangular matrix A,'
      write ( *, '(a)' ) '  solve A * x = b'

      job = 01

      call dtrsl ( a, lda, n, b, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
      end do

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(j,i) * x(j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For an upper triangular matrix A,'
      write ( *, '(a)' ) '  solve A'' * x = b'

      job = 11

      call dtrsl ( a, lda, n, b, job, info )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The solution (should be 1,2,3,4,5):'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,g14.6)' ) i, b(i)
      end do

      return
      end
      subroutine dmat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc DMAT_UNIFORM_01 fills an array with unit pseudorandom numbers.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

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
