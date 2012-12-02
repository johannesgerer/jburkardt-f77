      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS3_PRB.
c
c  Discussion:
c
c    BLAS3_PRB tests the BLAS3 routines.
c
c  Modified:
c
c    15 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS3_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS3 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS3_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests a Cholesky solver.
c
c  Discussion:
c
c    Solve a positive definite symmetric linear system A*X = B.
c
      implicit none

      integer n
      parameter ( n = 20 )

      double precision a(n,n)
      integer i
      integer info
      integer j
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Demonstrate Cholesky factor/solve routines'
      write ( *, '(a)' ) '  built out of level 1, 2 and 3 BLAS.'
      write ( *, '(a)' ) ' '
c
c  Set the entries of the matrix.
c
      do i = 1, n

        do j = 1, n
          a(i,j) = 0.0D+00
        end do

        a(i,i) = 2.0D+00

        if ( 1 .lt. i ) then
          a(i,i-1) = -1.0D+00
        end if

        if ( i .lt. n ) then
          a(i,i+1) = -1.0D+00
        end if

      end do
c
c  Set the right hand side vector.
c
      do i = 1, n-1
        x(i) = 0.0D+00
      end do
      x(n) = dble ( n + 1 )
c
c  Factor the matrix.
c
      call dlltb ( n, a, n, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The matrix is singular!'
        return
      end if
c
c  Solve the system.
c
      call dllts ( n, a, n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 entries of solution:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) ( x(i), i = 1, 10 )

      return
      end
      subroutine dlltb ( n, a, lda, info )

c*********************************************************************72
c
cc DLLTB Cholesky factors a positive definite symmetric matrix.
c
c  Discussion:
c
c    The factorization has the form A = L * L'.
c
c    The parameter NB determines the 'blocking factor' used
c    by the level 3 BLAS routines.
c
c    DLLT is an equivalent routine, but it solves the problem
c    at a lower level.
c
c    DLLTB solves the problem in chunks,
c    which is hoped to allow for greater optimization.
c
      implicit none

      integer lda
      integer n
      integer nb

      parameter ( nb = 64 )

      double precision a(lda,*)
      integer info
      integer j
      integer jb

      info = 0

      do j = 1, n, nb

        jb = min ( nb, n-j+1 )
c
c  Update diagonal block
c
        call dsyrk ( 'lower', 'no transpose', jb, j-1, -1.0D+00,
     &    a(j,1), lda, 1.0D+00, a(j,j), lda )
c
c  Factorize diagonal block and test for non-positive-definiteness.
c
        call dllt ( jb, a(j,j), lda, info )

        if ( info .ne. 0 )then
          info = info + j - 1
          return
        end if

        if ( j + jb .le. n ) then
c
c  Update subdiagonal block.
c
          call dgemm ( 'no transpose', 'transpose', n-j-jb+1, jb,
     &        j-1, -1.0D+00, a(j+jb,1), lda, a(j,1), lda, 1.0D+00,
     &        a(j+jb,j), lda )
c
c  Compute the subdiagonal block of L.
c
          call dtrsm ( 'right', 'lower', 'transpose', 'non-unit',
     &      n-j-jb+1, jb, 1.0D+00, a(j,j), lda, a(j+jb,j), lda )

        end if

      end do

      return
      end
      subroutine dllt ( n, a, lda, info )

c*********************************************************************72
c
cc DLLT computes an L*L' factorization of an SPD matrix.
c
c  Discussion:
c
c    DLLT uses level 2 and level 1 BLAS routines.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      integer info
      integer j
      double precision ddot

      info = 0

      do j = 1, n
c
c  Update A(J,J).
c
        a(j,j) = a(j,j) - ddot ( j-1, a(j,1), lda, a(j,1), lda )
c
c  Compute L(J,J) and test for non-positive-definiteness.
c
        if ( a(j,j) .le. 0.0D+00 ) then
          info = j
          return
        end if

        a(j,j) = sqrt ( a(j,j) )
c
c  Update elements J+1 to N of J-th column.
c
        if ( j .lt. n ) then
          call dgemv ( 'no transpose', n-j, j-1, -1.0D+00, a(j+1,1),
     &      lda, a(j,1), lda, 1.0D+00, a(j+1,j), 1 )
c
c  Compute elements J+1 to N of J-th column of L.
c
          call dscal ( n-j, 1.0D+00/a(j,j), a(j+1,j), 1 )

        end if

      end do

      return
      end
      subroutine dllts ( n, a, lda, b )

c*********************************************************************72
c
cc DLLTS solves A*X=B, for a positive definite symmetric matrix.
c
c  Discussion:
c
c    The matrix A should have been factored by DLLT or DLLTB.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      double precision b(n)
      integer k
      double precision ddot
      double precision t

      do k = 1, n
        t = ddot ( k-1, a(k,1), lda, b(1), 1 )
        b(k) = ( b(k) - t ) / a(k,k)
      end do

      do k = n, 1, -1
        b(k) = b(k) / a(k,k)
        t = - b(k)
        call daxpy ( k-1, t, a(k,1), lda, b(1), 1 )
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
