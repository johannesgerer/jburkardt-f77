      program main

c*********************************************************************72
c
cc MAIN is the main program for LAPACK_PRB.
c
c  Discussion:
c
c    LAPACK_PRB calls some of the LAPACK routines.
c
c  Modified:
c
c    04 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAPACK library.'

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
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAPACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DGBTRF and DGBTRS.
c
c  Discussion:
c
c    The problem is just an enlarged version of the
c    problem for n = 5, which is:
c
c    Matrix A is ( 2 -1  0  0  0)    right hand side b is  (1)
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
c    09 April 2007
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

      parameter ( n = 25 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipiv(n)
      integer j
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general band storage mode (GB):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGBTRF factors a general band matrix.'
      write ( *, '(a)' ) '  DGBTRS solves a factored system.'
      write ( *, '(a)' ) ' '
c
c  Assign values to matrix A and right hand side b.
c
      b(1) = 1.0D+00
      do i = 2, n-1
        b(i) = 0.0D+00
      end do
      b(n) = 1.0D+00
c
c  Zero out the matrix.
c
      do i = 1, lda
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      m = ml + mu + 1
c
c  Superdiagonal,
c  Diagonal,
c  Subdiagonal.
c
      do j = 2, n
        a(m-1,j) = -1.0D+00
      end do

      do j = 1, n
        a(m,j) = 2.0D+00
      end do

      do j = 1, n-1
        a(m+1,j) = -1.0D+00
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Bandwidth is ', m
      write ( *, '(a)' ) ' '

      call dgbtrf ( n, n, ml, mu, a, lda, ipiv, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      call dgbtrs ( 'n', n, ml, mu, 1, a, lda, ipiv, b, n, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a,i8)' ) '  Solution failed, INFO = ', info
        return
      end if

      call r8vec_print_some ( n, b, 5, 
     &  '  Partial solution (all should be 1)' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DGBTRF and DGBTRS.
c
c  Modified:
c
c    09 April 2007
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

      parameter ( n = 25 )
      parameter ( ml = 3 )
      parameter ( mu = 3 )
      parameter ( lda = 2 * ml + mu + 1 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer info
      integer ipiv(n)
      integer j
      integer m
      double precision temp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general band storage mode (GB):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGBTRF factors a general band matrix.'
      write ( *, '(a)' ) '  DGBTRS solves a factored system.'
      write ( *, '(a)' ) ' '
c
c  Assign values to matrix A and right hand side b.
c
c  We want to try a problem with a significant bandwidth.
c
      m = ml + mu + 1

      do j = 1, n
        ilo = max ( 1, j-mu )
        ihi = min ( n, j+ml )
 
        temp = 0.0D+00
        do i = ilo, ihi
          a(i-j+m,j) = -1.0D+00
          temp = temp - 1.0D+00
        end do
        temp = temp + 1.0D+00
 
        a(m,j) = 4.0D+00 - temp
      end do
c
c  Right hand side.
c
      do j = 1, n
        b(j) = 4.0D+00
      end do
c
c  Factor the matrix.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Bandwidth is ', m
 
      call dgbtrf ( n, n, ml, mu, a, lda, ipiv, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      call dgbtrs ( 'n', n, ml, mu, 1, a, lda, ipiv, b, n, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Solution failed, INFO = ', info
      end if

      call r8vec_print_some ( n, b, 5, 
     &  '  Partial solution (all should be 1)' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DGECON and DGETRF.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer lwork

      parameter ( n = 3 )
      parameter ( lda = n )
      parameter ( lwork = 4 * n )

      double precision a(lda,n)
      double precision anorm
      integer info
      integer ipiv(n)
      integer iwork(n)
      double precision rcond
      double precision r8mat_norm_li
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGETRF computes the LU factorization;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGECON computes the condition number '
      write ( *, '(a)' ) '  of a factored matrix'
c
c  Set the matrix.
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

      call r8mat_print ( n, n, a, '  The matrix A:' )
c
c  Compute the L-Infinity norm of the matrix.
c
      anorm = r8mat_norm_li ( n, n, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Matrix L-infinity norm is ', anorm
c
c  Factor the matrix.
c
      call dgetrf ( n, n, a, lda, ipiv, info )
c
c  Get the condition number.
c
      call dgecon ( 'I', n, a, lda, anorm, rcond, work, iwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Condition number calculation failedc'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Matrix reciprocal condition number = ', rcond
 
      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DGESVD.
c
c  Discussion:
c
c    DGESVD computes the singular value decomposition:
c
c      A = U * S * V'
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n
      integer lwork

      parameter ( m = 6 )
      parameter ( n = 4 )
c     parameter ( lwork = 3*min(m,n) + max ( max(m,n), 2*min(m,n) ) )
      parameter ( lwork = 3*4 + 8 )

      double precision a(m,n)
      double precision b(m,n)
      integer i
      integer info
      integer j
      integer k
      integer lda
      integer ldu
      integer ldvt
      character jobu
      character jobvt
      integer mn
      double precision s(m+n)
      integer seed
      double precision sigma(m,n)
      double precision u(m,m)
      double precision vt(n,n)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGESVD computes the singular value '
      write ( *, '(a)' ) '  decomposition:'
      write ( *, '(a)' ) '    A = U * S * V'''
c
c  Set A.
c
      seed = 123456789

      call r8mat_uniform_01 ( m, n, seed, a )

      call r8mat_print ( m, n, a, '  The matrix A:' )
c
c  Compute the eigenvalues and eigenvectors.
c
      jobu = 'A'
      jobvt = 'A'
      lda = m
      ldu = m
      ldvt = n
      
      call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, 
     &  work, lwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
        return
      end if

      mn = min ( m, n )

      call r8vec_print ( mn, s, '  Singular values' )

      call r8mat_print ( m, m, u, '  Left singular vectors U:' )
      call r8mat_print ( n, n, vt, '  Right singular vectors V'':' )

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
          b(i,j) = 0.0D+00
          do k = 1, n
            b(i,j) = b(i,j) + sigma(i,k) * vt(k,j)
          end do
        end do
      end do

      do i = 1, m
        do j = 1, n
          a(i,j) = 0.0D+00
          do k = 1, m
            a(i,j) = a(i,j) + u(i,k) * b(k,j)
          end do
        end do
      end do

      call r8mat_print ( m, n, a, '  The product U * S * V'':' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests DGETRF and DGETRI.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer lwork

      parameter ( n = 3 )
      parameter ( lda = n )
      parameter ( lwork = n )

      double precision a(lda,n)
      integer info
      integer ipiv(n)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGETRF factors a general matrix;'
      write ( *, '(a)' ) '  DGETRI computes the inverse.'
c
c  Set the matrix.
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

      call r8mat_print ( n, n, a, '  The matrix A:' )
c
c  Factor the matrix.
c
      call dgetrf ( n, n, a, lda, ipiv, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
        write ( *, '(a)' ) '  The matrix is numerically singular.'
        return
      end if
c
c  Compute the inverse matrix.
c
      call dgetri ( n, a, lda, ipiv, work, lwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The inversion procedure failedc'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8mat_print ( n, n, a, '  The inverse matrix:' )
 
      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests DGETRF and DGETRS.
c
c  Modified:
c
c    09 April 2007
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
      integer info
      integer ipiv(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGETRF computes the LU factorization;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGETRS solves linear systems using '
      write ( *, '(a)' ) '  the LU factors;'
c
c  Set the matrix.
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

      call r8mat_print ( n, n, a, '  The matrix A:' )
c
c  Factor the matrix.
c
      call dgetrf ( n, n, a, lda, ipiv, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
        write ( *, '(a)' ) '  The matrix is numerically singular.'
        return
      end if
c
c  Set the right hand side.
c
      b(1) = 14.0D+00
      b(2) = 32.0D+00
      b(3) = 23.0D+00

      call r8vec_print ( n, b, '  Right hand side B' )
c
c  Solve the linear system.
c
      call dgetrs ( 'n', n, 1, a, lda, ipiv, b, n, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Solution procedure failedc'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8vec_print ( n, b, '  The solution X' )
 
      return
      end
      subroutine test07 ( )
c
c*********************************************************************72
c
cc TEST07 tests DGETRF and DGETRS.
c
c  Discussion:
c
c    The problem is just an enlarged version of the
c    problem for n = 5, which is:
c
c    Matrix A is ( N -1 -1 -1 -1)    right hand side b is  (1)
c                (-1  N -1 -1 -1)                          (1)
c                (-1 -1  N -1 -1)                          (1)
c                (-1 -1 -1  N -1)                          (1)
c                (-1 -1 -1 -1  N)                          (1)
c
c    Solution is   (1)
c                  (1)
c                  (1)
c                  (1)
c                  (1)
c
c    For this problem, no pivoting is required.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 25 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer ipiv(n)
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGETRF factors a general matrix;'
      write ( *, '(a)' ) '  DGETRS solves a linear system;'
      write ( *, '(a)' ) ' '
c
c  Assign values to matrix A and right hand side b.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = dble ( n )
          else
            a(i,j) = -1.0D+00
          end if
        end do
      end do

      do i = 1, n
        b(i) = 1.0D+00
      end do
c
c  Factor the matrix.
c
      call dgetrf ( n, n, a, lda, ipiv, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix is singular, INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      call dgetrs ( 'n', n, 1, a, lda, ipiv, b, n, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  Solution procedure failed, INFO = ', info
        return
      end if

      call r8vec_print_some ( n, b, 5, 
     &  '  Partial solution (all should be 1)' )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests DGTSV.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 100 )

      double precision b(n)
      double precision c(n-1)
      double precision d(n)
      double precision e(n-1)
      integer i
      integer info
      integer ldb
      integer nrhs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general tridiagonal storage mode (GT):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGTSV factors and solves a linear system'
      write ( *, '(a)' ) '  with a general tridiagonal matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The system is of order N = ', n
      write ( *, '(a)' ) ' '
c
c  Right hand side.
c 
      do i = 1, n-1
        b(i) = 0.0D+00
      end do
      b(n) = n + 1
c
c  Subdiagonal.
c  Diagonal.
c  Superdiagonal.
c 
      do i = 1, n-1
        c(i) = -1.0D+00
      end do

      do i = 1, n
        d(i) = 2.0D+00
      end do

      do i = 1, n-1
        e(i) = -1.0D+00
      end do
 
      nrhs = 1
      ldb = n
c
c  Factor and solve the linear system.
c
      call dgtsv ( n, nrhs, c, d, e, b, ldb, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Solution procedure failed.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8vec_print_some ( n, b, 5, 
     &  '  Partial solution (Should be 1,2,3...)' )

      return
      end
      subroutine test09 ( )
c
c*********************************************************************72
c
cc TEST09 tests DPBTRF.
c
c  Discussion:
c
c    We want to compute the lower triangular Cholesky factor L
c    of a positive definite symmetric band matrix.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer nband

      parameter ( n = 5 )
      parameter ( nband = 1 )

      double precision a(nband+1,n)
      integer i
      integer info
      integer j
      double precision l(nband+1,n)
      double precision l_row(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) 
     &  '  in positive definite band storage mode (PB):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DPBTRF computes'
      write ( *, '(a)' ) '    the lower Cholesky factor A = L*L'' or'
      write ( *, '(a)' ) '    the upper Cholesky factor A = U''*U;'
      write ( *, '(a)' ) ' '
c
c  Zero out the matrix.
c
      do i = 1, nband+1
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Store the diagonal of a symmetric band matrix.
c
      do j = 1, n
        a(1,j) = 2.0D+00
      end do
c
c  Store the subdiagonal of a symmetric band matrix.
c
      do j = 1, n-1
        a(2,j) = -1.0D+00
      end do
c
c  Get the lower triangular Cholesky factor L:
c
      do i = 1, nband+1
        do j = 1, n
          l(i,j) = a(i,j)
        end do
      end do

      call dpbtrf ( 'L', n, nband, l, nband+1, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        return
      end if
c
c  Print the relevant entries of L:
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The lower Cholesky factor L:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        do j = 1, n

          if ( 0 .le. i - j .and. i-j .le. nband ) then
            l_row(j) = l(i-j+1,j)
          else
            l_row(j) = 0.0D+00
          end if

        end do

        write ( *, '(5f10.6)' ) ( l_row(j), j = 1, n )

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests DPBTRF and DPBTRS.
c
c  Discussion:
c
c    The problem is just an enlarged version of the
c    problem for n = 5, which is:
c
c    Matrix A is ( 2 -1  0  0  0)    right hand side b is  (1)
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
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer nband
      integer lda

      parameter ( n = 25 )
      parameter ( nband = 1 )
      parameter ( lda = nband + 1 )

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer info
      integer j 
      integer nrhs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  For a double precision real matrix (D) in'
      write ( *, '(a)' ) '  positive definite band storage mode (PB):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For a positive definite symmetric band'
      write ( *, '(a)' ) '  matrix:'
      write ( *, '(a)' ) '  DPBTRF factors;'
      write ( *, '(a)' ) '  DPBTRS solves linear systems.'
      write ( *, '(a)' ) ' '
c
c  Zero out the matrix.
c
      do i = 1, lda
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Super (and sub) diagonal.
c
      do j = 2, n
        a(1,j) = -1.0D+00
      end do
c
c  Diagonal.
c
      do j = 1, n
        a(2,j) = 2.0D+00
      end do
c
c  Set the right hand side.
c
      b(1) = 1.0D+00
      do i = 2, n-1
        b(i) = 0.0D+00
      end do
      b(n) = 1.0D+00
c
c  Factor the matrix.
c
      call dpbtrf ( 'u', n, nband, a, lda, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        return
      end if
c
c  Solve the linear system.
c
      nrhs = 1
      call dpbtrs ( 'u', n, nband, nrhs, a, lda, b, n, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Solution failed, INFO = ', info
      end if

      call r8vec_print_some ( n, b, 5, 
     &  '  Partial solution (all should be 1)' )

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests DPOTRF and DPOTRI.
c
c  Modified:
c
c    09 April 2007
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

      double precision a(n,n)
      double precision a_inv(n,n)
      integer i
      integer info
      integer j
      integer k
      double precision r(n,n)
      double precision temp(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in positive definite storage mode (PO):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DPOTRF computes the Cholesky factor '
      write ( *, '(a)' ) '    A = R''*R;'
      write ( *, '(a)' ) '  DPOTRI computes the inverse.'
      write ( *, '(a)' ) ' '
c
c  Zero out the matrix.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Subdiagonal.
c
      do i = 2, n
        a(i,i-1) = -1.0D+00
      end do
c
c  Diagonal.
c
      do i = 1, n
        a(i,i) = 2.0D+00
      end do
c
c  Superdiagonal.
c
      do i = 1, n - 1
        a(i,i+1) = -1.0D+00
      end do

      call r8mat_print ( n, n, a, '  The matrix A:' )
c
c  Factor the matrix.
c
      do i = 1, n
        do j = 1, n
          r(i,j) = a(i,j)
        end do
      end do

      call dpotrf ( 'u', n, r, lda, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a,i8)' ) '  DPOTRF returns INFO = ', info
        return
      end if

      do i = 1, n
        do j = 1, i-1
          r(i,j) = 0.0D+00
        end do
      end do

      call r8mat_print ( n, n, r, '  The Cholesky factor R:' )

      do i = 1, n
        do j = 1, n
          temp(i,j) = 0.0D+00
          do k = 1, n
            temp(i,j) = temp(i,j) + r(i,k) * r(k,j)
          end do
        end do
      end do

      call r8mat_print ( n, n, temp, '  The product R'' * R' )
c
c  Compute the inverse matrix.
c
      do i = 1, n
        do j = 1, n
          a_inv(i,j) = r(i,j)
        end do
      end do

      call dpotri ( 'u', n, a_inv, lda, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The inversion procedure failed, INFO = ', info
        return
      end if

      do i = 1, n
        do j = 1, i-1
          a_inv(i,j) = a_inv(j,i)
        end do
      end do

      call r8mat_print ( n, n, a_inv, '  The inverse matrix B:' )

      do i = 1, n
        do j = 1, n
          temp(i,j) = 0.0D+00
          do k = 1, n
            temp(i,j) = temp(i,j) + a_inv(i,k) * a(k,j)
          end do
        end do
      end do

      call r8mat_print ( n, n, temp, '  The product B * A' )

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests DGEQRF and DORGQR.
c
c  Discussion:
c
c    DGEQRF computes the QR factorization of an M by N matrix A:
c
c      A(MxN) = Q(MxK) * R(KxN)
c
c    where K = min ( M, N ).
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n
      integer k
      integer lwork

      parameter ( m = 8 )
      parameter ( n = 6 )
c     parameter ( k = min ( m, n ) )
      parameter ( k = 6 )
      parameter ( lwork = n )

      double precision a(m,n)
      integer i
      integer info
      integer j
      integer l
      integer lda
      double precision q(m,k)
      double precision r(k,n)
      integer seed
      double precision tau(k)
      double precision work(lwork)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGEQRF computes the QR factorization:'
      write ( *, '(a)' ) '    A = Q * R'
      write ( *, '(a)' ) '  DORGQR computes the explicit form of '
      write ( *, '(a)' ) '  the Q factor.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this case, our M x N matrix A has more'
      write ( *, '(a)' ) '  rows than columns:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
c
c  Set A.
c
      call r8mat_uniform_01 ( m, n, seed, a )

      call r8mat_print ( m, n, a, '  The matrix A:' )
c
c  Compute the QR factorization.
c
      lda = m

      call dgeqrf ( m, n, a, lda, tau, work, lwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
        return
      end if

      do i = 1, k
        do j = 1, n
          r(i,j) = 0.0D+00
        end do
      end do

      do i = 1, k
        do j = i, n
          r(i,j) = a(i,j)
        end do
      end do
c
c  Construct Q explicitly.
c
      do i = 1, m
        do j = 1, k
          q(i,j) = a(i,j)
        end do
      end do

      call dorgqr ( m, k, k, q, lda, tau, work, lwork, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
        return
      end if

      call r8mat_print ( m, k, q, '  The Q factor:' )

      call r8mat_print ( k, n, r, '  The R factor:' )

      do i = 1, m
        do j = 1, n
          a(i,j) = 0.0D+00
          do l = 1, k
            a(i,j) = a(i,j) + q(i,l) * r(l,j)
          end do
        end do
      end do

      call r8mat_print ( m, n, a, '  The product Q * R:' )

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests DGEQRF and DORGQR.
c
c  Discussion:
c
c    DGEQRF computes the QR factorization of an M by N matrix A:
c
c      A(MxN) = Q(MxK) * R(KxN)
c
c    where K = min ( M, N ).
c
c    In order to get an M x M matrix Q, we are going to pad our
c    original matrix A with extra columns of zeroesc
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n
      integer n2
      integer k
      integer lwork

      parameter ( m = 8 )
      parameter ( n = 6 )
      parameter ( n2 = 8 )
!     parameter ( k = min ( m, n2 ) )
      parameter ( k = 8 )
      parameter ( lwork = n2 )

      double precision a(m,n)
      double precision a2(m,n2)
      integer i
      integer info
      integer j
      integer l
      integer lda
      double precision q(m,k)
      double precision r(k,n2)
      integer seed
      double precision tau(k)
      double precision work(lwork)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGEQRF computes the QR factorization:'
      write ( *, '(a)' ) '    A = Q * R'
      write ( *, '(a)' ) '  DORGQR computes the explicit form of '
      write ( *, '(a)' ) '  the Q factor.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this case, our M x N matrix A has more'
      write ( *, '(a)' ) '  rows than columns:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Normally, LAPACK will only return an'
      write ( *, '(a)' ) '  M x min(M,N) portion of Q.  When N .lt. M, '
      write ( *, '(a)' ) '  we lose information.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we force the computation of the '
      write ( *, '(a)' ) '  full Q by making a copy of A with N-M '
      write ( *, '(a)' ) '  extra zero columns.'
      write ( *, '(a)' ) ' '
c
c  Set A.
c
      call r8mat_uniform_01 ( m, n, seed, a )

      call r8mat_print ( m, n, a, '  The matrix A:' )
c
c  Copy A, and pad with extra columns.
c
      do i = 1, m
        do j = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      do i = 1, m
        do j = n+1, m
          a2(i,j) = 0.0D+00
        end do
      end do
c
c  Compute the QR factorization.
c
      lda = m

      call dgeqrf ( m, n2, a2, lda, tau, work, lwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
        return
      end if

      do i = 1, k
        do j = 1, n2
          r(i,j) = 0.0D+00
        end do
      end do

      do i = 1, k
        do j = i, n2
          r(i,j) = a2(i,j)
        end do
      end do
c
c  Construct Q explicitly.
c
      do i = 1, m
        do j = 1, k
          q(i,j) = a2(i,j)
        end do
      end do

      call dorgqr ( m, k, k, q, lda, tau, work, lwork, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
        return
      end if

      call r8mat_print ( m, k, q, '  The Q factor:' )

      call r8mat_print ( k, n, r, '  The R factor:' )

      do i = 1, m
        do j = 1, n2
          a2(i,j) = 0.0D+00
          do l = 1, k
            a2(i,j) = a2(i,j) + q(i,l) * r(l,j)
          end do
        end do
      end do

      call r8mat_print ( m, n2, a2, '  The product Q * R:' )

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests DGEQRF and DORGQR.
c
c  Discussion:
c
c    DGEQRF computes the QR factorization of an M by N matrix A:
c
c      A(MxN) = Q(MxK) * R(KxN)
c
c    where K = min ( M, N ).
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n
      integer k
      integer lwork

      parameter ( m = 6 )
      parameter ( n = 8 )
c     parameter ( k = min ( m, n ) )
      parameter ( k = 6 )
      parameter ( lwork = n )

      double precision a(m,n)
      integer i
      integer info
      integer j
      integer l
      integer lda
      double precision q(m,k)
      double precision r(k,n)
      integer seed
      double precision tau(k)
      double precision work(lwork)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in general storage mode (GE):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DGEQRF computes the QR factorization:'
      write ( *, '(a)' ) '    A = Q * R'
      write ( *, '(a)' ) '  DORGQR computes the explicit form of'
      write ( *, '(a)' ) '  the Q factor.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this case, our M x N matrix A has more'
      write ( *, '(a)' ) '  columns than rows:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
c
c  Set A.
c
      call r8mat_uniform_01 ( m, n, seed, a )

      call r8mat_print ( m, n, a, '  The matrix A:' )
c
c  Compute the QR factorization.
c
      lda = m

      call dgeqrf ( m, n, a, lda, tau, work, lwork, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
        return
      end if

      do i = 1, k
        do j = 1, n
          r(i,j) = 0.0D+00
        end do
      end do

      do i = 1, k
        do j = i, n
          r(i,j) = a(i,j)
        end do
      end do
c
c  Construct Q explicitly.
c
      do i = 1, m
        do j = 1, k
          q(i,j) = a(i,j)
        end do
      end do

      call dorgqr ( m, k, k, q, lda, tau, work, lwork, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
        return
      end if

      call r8mat_print ( m, k, q, '  The Q factor:' )

      call r8mat_print ( k, n, r, '  The R factor:' )

      do i = 1, m
        do j = 1, n
          a(i,j) = 0.0D+00
          do l = 1, k
            a(i,j) = a(i,j) + q(i,l) * r(l,j)
          end do
        end do
      end do

      call r8mat_print ( m, n, a, '  The product Q*R' )

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests DSBGVX.
c
c  Discussion:
c
c    DSBGVX deals with the generalized eigenvalue problem:
c
c      A * x = lambda * B * x
c
c    where A and B are symmetric and banded (and stored in LAPACK symmetric
c    band storage mode).  B is additionally assumed to be positive definite.
c
c    This is an "expert" interface, and the user is requesting
c    only some of the eigenvalues and eigenvectors.  In this example,
c    only the largest and smallest (in magnitude) eigenvalues will
c    be requested.
c
c  Modified:
c
c    17 March 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    double precision AB(LDAB,N), contains, on input, the upper or lower
c    triangle of the symmetric band matrix A, stored in the first KA+1 rows
c    of the array AB.  
c    If UPLO = 'U', then
c      AB(KA+1+I-J,J) = A(I,J) for max(1,J-KA) .le. I .le. J;
c    If UPLO = 'L', then
c      AB(1+I-J,J) = A(I,J) for J .le. I .le. min(N,J+KA).
c
c    double precision ABSTOL, the absolute error tolerance for the eigenvalues.
c    If the input value of ABSTOL is not positive, then an appropriate
c    value will be determined internally and used instead.
c
c    double precision BB(LDBB,N), contains, on input, the upper or lower
c    triangle of the positive definite symmetric band matrix B, stored in 
c    the first KB+1 rows of the array BB.  
c    If UPLO = 'U', then
c      BB(KB+1+I-J,J) = B(I,J) for max(1,J-KB) .le. I .le. J;
c    If UPLO = 'L', then
c      BB(1+I-J,J) = B(I,J) for J .le. I .le. min(N,J+KB).
c
c    integer IFAIL(N), if JOBZ = 'V', then if INFO is 0, the first
c    M elements of IFAIL have been set to zero by DSBGVX, but if INFO
c    is nonzero, IFAIL contains the indices of the eigenvalues
c    for which the eigenvectors failed to converge.  If JOBZ = 'N',
c    then IFAIL is not referenced.
c
c    integer IL, IU, the indices of the first (smallest) and last
c    (largest) eigenvalues to be returned.  These values are only
c    used if RANGE = 'I'.  It must be the case that 1 .le. IL .le. IU .le. N.
c
c    Integer INFO, is 0 for a successful computation, 
c    negative if an input argument was illegal (the index of that
c    argument is the value of -INFO), or positive, in which case,
c    if 0 .lt. INFO .le. N, then INFO off diagonal elements of an
c    intermediate tridiagonal form did not converge to zero, or
c    if N .lt. INFO, B is not positive definite and the computation
c    could not be completed.
c
c    integer IWORK(5*N), workspace.
c
c    character JOBZ, is 'N' if only eigenvalues are desired, or 'V'
c    if eigenvectors will also be required.
c
c    Integer KA, the number of superdiagonals (if UPLO = 'U') or
c    subdiagonals (if UPLO = 'L') of A that are nonzero.
c
c    integer KB, the number of superdiagonals (if UPLO = 'U') or
c    subdiagonals (if UPLO = 'L') of B that are nonzero.
c
c    integer LDAB, the leading dimension of the array AB, which
c    must be at least KA+1.
c
c    integer LDBB, the leading dimension of the array BB, which
c    must be at least KB+1.
c
c    integer LDQ, the leading dimension of the array Q.
c    If JOBZ = 'N', then Q is not used, and LDQ should be any
c    positive value such as 1.  If JOBZ = 'V', then LDQ must be
c    at least N.
c
c    integer LDZ, the leading dimension of the array Z.
c    If JOBZ = 'N', then Z is not used, and LDZ should be any
c    positive value such as 1.  If JOBZ = 'V', then LDZ must be
c    at least N.
c
c    integer M, the number of eigenvalues found by DSBGVX.
c
c    integer N, the order of the matrices A and B.
c
c    double precision Q(LDQ,N), if JOBZ = 'V', the N by N matrix used to 
c    reduce the problem to standard form: "C * x = lambda * x"
c    and then to reduce the matrix C to tridiagonal form.  But
c    if JOBZ is not 'V', Q is not referenced.
c
c    character RANGE, specifies which eigenvalues are desired.
c    'A' means all, 'V' means a real interval will be specified in which
c    eigenvalues are to be sought, 'I' means a range of indices will
c    be specified.
c
c    character UPLO, is 'U' if the upper triangles of A and B are stored,
c    'L' if the lower triangles are stored.
c
c    double precision VL, VU, the lower and upper bounds of an interval to be
c    searched for eigenvalues.  In this case, VL must be less than VU.  
c    These values are used only if RANGE = 'V'.
c
c    double precision W(N), the requested eigenvalues, in ascending order.
c
c    double precision WORK(7*N), workspace.
c
c    double precision Z(LDZ,N), if JOBZ = 'V', the I-th column of Z contains
c    the eigenvector associated with the I-th eigenvalue W(I).
c
      implicit none

      integer n
      integer ka
      integer kb
      integer ldab
      integer ldbb
      integer ldq
      integer ldz

      parameter ( n = 11 )
      parameter ( ka = 2 )
      parameter ( kb = 1 )
      parameter ( ldab = ka+1 )
      parameter ( ldbb = kb+1 )
      parameter ( ldq = 1 )
      parameter ( ldz = 1 )

      double precision ab(ldab,n)
      double precision abstol
      double precision bb(ldbb,n)
      integer i
      integer ifail(n)
      integer il
      integer ilo
      integer info
      integer iu
      integer iwork(5*n)
      integer j
      character jobz
      integer m
      double precision q(ldq,n)
      character range
      integer test
      character uplo
      double precision value
      double precision vl
      double precision vu
      double precision w(n)
      double precision work(7*n)
      double precision z(ldz,n)

      abstol = 0.0D+00
      jobz = 'N'
      range = 'I'
      uplo = 'U'
      vl = 0.0D+00  
      vu = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in symmetric band storage mode (SB):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For a symmetric banded NxN matrix A, '
      write ( *, '(a)' ) '  and a symmetric banded positive definite '
      write ( *, '(a)' ) '  NxN matrix B, DSBGVX solves the '
      write ( *, '(a)' ) '  generalized eigenvalue problem'
      write ( *, '(a)' ) '    A * X = LAMBDA * B * X'
      write ( *, '(a)' ) ' '

      do test = 1, 2
c
c  Set A.
c
        do j = 1, n
          ilo = max ( j - ka, 1 )
          do i = ilo, j

            if ( j .eq. i-2 ) then
              value = -1.0D+00
            else if ( j .eq. i-1 ) then
              value = -1.0D+00
            else if ( j .eq. i ) then
              value = +4.0D+00
            else if ( j .eq. i+1 ) then
              value = -1.0D+00
            else if ( j .eq. i+2 ) then
              value = -1.0D+00
            else
              value = 0.0D+00
            end if

            ab(ka+1+i-j,j) = value

          end do
        end do
c
c  Set B.
c
        do j = 1, n
          ilo = max ( j - kb, 1 )
          do i = ilo, j

            if ( j .eq. i-1 ) then
              value = -1.0D+00
            else if ( j .eq. i ) then
              value = +2.0D+00
            else if ( j .eq. i+1 ) then
              value = -1.0D+00
            else
              value = 0.0D+00
            end if

            bb(kb+1+i-j,j) = value

          end do
        end do
c
c  Request the value of the SMALLEST or LARGEST eigenvalue:
c
        if ( test .eq. 1 ) then
          il = 1
          iu = 1
        else if ( test .eq. 2 ) then
          il = n
          iu = n
        end if

        call dsbgvx ( jobz, range, uplo, n, ka, kb, ab, ldab, bb, 
     &    ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, 
     &    iwork, ifail, info )
 
        if ( info .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) 
     &      '  Illegal value for input argument ', -info
          return
        else if ( 0 .lt. info ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  The eigenvalue or eigenvector'
          write ( *, '(a)' ) '  iterations did not converge.'
          cycle
        end if

        call r8vec_print ( m, w, '  Computed eigenvalues' )

      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests DSYEV.
c
c  Modified:
c
c    09 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lwork

      parameter ( n = 7 )
      parameter ( lwork = 3 * n - 1 )

      double precision a(n,n)
      integer info
      character jobz
      double precision lambda(n)
      character uplo
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  For a double precision real matrix (D)'
      write ( *, '(a)' ) '  in symmetric storage mode (SY):'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For a symmetric matrix in general storage,'
      write ( *, '(a)' ) '  DSYEV computes eigenvalues and '
      write ( *, '(a)' ) '  eigenvectors;'
      write ( *, '(a)' ) ' '
c
c  Set A.
c
      call clement2 ( n, a )

      call r8mat_print ( n, n, a, '  The matrix A:' )
c
c  Compute the eigenvalues and eigenvectors.
c
      jobz = 'V'
      uplo = 'U'

      call dsyev ( jobz, uplo, n, a, n, lambda, work, lwork, info )
 
      if ( info .ne. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DSYEV returned nonzero INFO = ', info

      else

        call r8vec_print ( n, lambda, '  The eigenvalues:' )

        if ( jobz .eq. 'V' ) then
          call r8mat_print ( n, n, a, '  The eigenvector matrix:' )
        end if
      
      end if

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests DGEQRF to QR factor a matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer m
      integer n
      integer nrhs
      integer tau_size

      parameter ( m = 5 )
      parameter ( n = 5 )
      parameter ( nrhs = 1 )

      parameter ( lda = m )
      parameter ( ldb = max ( m, n ) )
      parameter ( lwork = n )
      parameter ( tau_size = n )

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer i
      integer info
      integer j
      double precision tau(tau_size)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  Call the standard LAPACK routine'
      write ( *, '(a)' ) '  DGEQRF to get the QR factorization of'
      write ( *, '(a)' ) '  a matrix stored in GE format.'

      do j = 1, n
        do i = 1, m
          if ( j == i - 1 ) then
            a(i,j) = - 1.0D+00
          else if ( j == i ) then
            a(i,j) = 2.0D+00
          else if ( j == i + 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      call r8mat_print ( m, n, a, '  Input matrix:' )

      call dgeqrf ( m, n, a, lda, tau, work, lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRF called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRF returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8mat_print ( m, n, a, '  Factored matrix:' )

      call r8vec_print ( tau_size, tau, '  Tau:' );
c
c  Set up and solve a linear system using DQEQRS.
c
      do i = 1, n
        b(i,1) = 0.0D+00
      end do
      b(n,1) = dble ( n + 1 )

      call dgeqrs ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution from DGEQRS:'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(2x,g14.6)' ) b(i,1)
      end do

      return
      end
      subroutine clement2 ( n, a )

c*********************************************************************72
c
cc CLEMENT2 returns the Clement2 matrix.
c
c  Formula:
c
c    if ( J = I+1 )
c      A(I,J) = sqrt(I*(N-I))
c    else if ( I = J+1 )
c      A(I,J) = sqrt(J*(N-J))
c    else
c      A(I,J) = 0
c
c  Example:
c
c    N = 5
c
c       .    sqrt(4)    .       .       .
c    sqrt(4)    .    sqrt(6)    .       .
c       .    sqrt(6)    .    sqrt(6)    .
c       .       .    sqrt(6)    .    sqrt(4)
c       .       .       .    sqrt(4)    .
c
c  Properties:
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    The diagonal of A is zero.
c
c    A is singular if N is odd.
c
c    About 64 percent of the entries of the inverse of A are zero.
c
c    The eigenvalues are plus and minus the numbers
c
c      N-1, N-3, N-5, ..., (1 or 0).
c
c    If N is even,
c
c      det ( A ) = (-1)**(N/2) * (N-1) * (N+1)**(N/2)
c
c    and if N is odd,
c
c      det ( A ) = 0
c
c  Reference:
c
c    P A Clement,
c    A class of triple-diagonal matrices for test purposes,
c    SIAM Review,
c    Volume 1, 1959, pages 50-52.
c
c  Modified:
c
c    15 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of A.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n

          if ( j .eq. i + 1 ) then
            a(i,j) = sqrt ( dble ( i * ( n - i ) ) )
          else if ( i .eq. j + 1 ) then
            a(i,j) = sqrt ( dble ( j * ( n - j ) ) )
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine dgeqrs ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork,
     &  info )

c*********************************************************************72
c
cc DGEQRS solves a linear system factored by DGEQRF.
c
c  Discussion:
c
c    This routine solves the least squares problem
c      min || A*X - B ||
c    using the QR factorization
c      A = Q*R
c    computed by DGEQRF.
c
c  Parameters:
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  M >= N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of columns of B.  NRHS >= 0.
c
c  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
c          Details of the QR factorization of the original matrix A as
c          returned by DGEQRF.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= M.
c
c  TAU     (input) DOUBLE PRECISION array, dimension (N)
c          Details of the orthogonal matrix Q.
c
c  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
c          On entry, the m-by-nrhs right hand side matrix B.
c          On exit, the n-by-nrhs solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B. LDB >= M.
c
c  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
c
c  LWORK   (input) INTEGER
c          The length of the array WORK.  LWORK must be at least NRHS,
c          and should be at least NRHS*NB, where NB is the block size
c          for this environment.
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer n
      integer nrhs

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer info
      integer m
      double precision one
      parameter ( one = 1.0D+00 )
      double precision tau(n)
      double precision work(lwork)

      info = 0

      if ( m .lt. 0 ) then
        info = - 1
      else if ( n .lt. 0 .or. n .gt. m ) then
        info = - 2
      else if ( nrhs .lt. 0 ) then
        info = - 3
      else if ( lda .lt. max( 1, m ) ) then
        info = - 5
      else if ( ldb .lt. max( 1, m ) ) then
        info = - 8
      else if ( lwork .lt. 1 .or. lwork .lt. nrhs .and. 
     &  m .gt. 0 .and. n .gt. 0 ) then
         info = - 10
      end if

      if ( info .ne. 0 ) then
         call xerbla ( 'DGEQRS', - info )
         return
      end if
c
c  Quick return if possible.
c
      if ( n .eq. 0 .or. nrhs .eq. 0 .or. m .eq. 0 ) then
        return
      end if
c
c  Compute B := Q' * B
c
      call dormqr( 'Left', 'Transpose', m, nrhs, n, a, lda, tau, b, ldb,
     &  work, lwork, info )
c
c  Solve R * X = B(1:n,:)
c
      call dtrsm ( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs,
     &  one, a, lda, b, ldb )

      return
      end
      function r8mat_norm_li ( m, n, a )

c*********************************************************************72
c
cc R8MAT_NORM_LI returns the matrix L-infinity norm of an R8MAT.
c
c  Definition:
c
c    The matrix L-infinity norm is defined as:
c
c      R8MAT_NORM_LI =  Max ( 1 .le. I .le. M ) Sum ( 1 .le. J .le. N ) abs ( A(I,J) ).
c
c    The matrix L-infinity norm is derived from the vector L-infinity norm,
c    and satisifies:
c
c      R8VEC_NORM_LI ( A*x ) .le. R8MAT_NORM_LI ( A ) * R8VEC_NORM_LI ( x ).
c
c  Modified:
c
c    24 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix whose L-infinity norm is
c    desired.
c
c    Output, double precision R8MAT_NORM_LI, the L-infinity norm of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8mat_norm_li
      double precision temp

      r8mat_norm_li = 0.0D+00

      do i = 1, m
        temp = 0.0D+00
        do j = 1, n
          temp = temp + abs ( a(i,j) )
        end do
        r8mat_norm_li = max ( r8mat_norm_li, temp )
      end do

      return
      end
      subroutine r8mat_print ( m, n, a, title )
c
c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c
c  Modified:
c
c    12 September 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title to be printed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character * ( * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Modified:
c
c    12 September 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer incx
      integer m
      integer n

      parameter ( incx = 5 )

      double precision a(m,n)
      character ( len = 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc)
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) .eq. dble ( int ( a(i,j) ) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 fills an array with unit pseudorandom numbers.
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
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Modified:
c
c    22 August 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,g16.8)' ) i, a(i)
      end do

      return
      end
      subroutine r8vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_SOME prints "some" of an R8VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Modified:
c
c    19 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer max_print
      character * ( * ) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print-2
          write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
        end do
        write ( *, '(a)' ) '  ......  ..............'
        i = n
        write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
        end do
        i = max_print
        write ( *, '(2x,i8,2x,g14.6,2x,a)' ) 
     &    i, a(i), '...more entries...'

      end if

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
