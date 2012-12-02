      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA006_PRB.
c
c  Discussion:
c
c    ASA006_PRB calls the ASA006 routines.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA006_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA006 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA006_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of CHOLESKY.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )
      integer nn_max
      parameter ( nn_max = (n_max*(n_max+1))/2 )

      double precision a(nn_max)
      double precision diff
      integer i
      integer ifault
      integer j
      integer k
      integer l
      integer n
      integer nn
      integer nullty
      double precision u(nn_max)
      double precision ufull(n_max,n_max)
      double precision utu

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  CHOLESKY computes the Cholesky factor'
      write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
      write ( *, '(a)' ) '  A compressed storage format is used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we look at the matrix A which is'
      write ( *, '(a)' ) '  N+1 on the diagonal and'
      write ( *, '(a)' ) '  N   on the off diagonals.'

      do n = 1, n_max

        nn = ( n * ( n + 1 ) ) / 2
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
        k = 0
        do i = 1, n
          do j = 1, i - 1
            k = k + 1
            a(k) = dble ( n )
          end do
          k = k + 1
          a(k) = dble ( n + 1 )
        end do

        call cholesky ( a, n, nn, u, nullty, ifault )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

        k = 0
        do j = 1, n
          do i = 1, j
            k = k + 1
            ufull(i,j) = u(k)
          end do
          do i = j + 1, n
            ufull(i,j) = 0.0D+00
          end do
        end do
!
!  Compute U' * U and compare to A.
!
        k = 0
        diff = 0.0D+00
        do i = 1, n
          do j = 1, i
            k = k + 1
            utu = 0.0D+00
            do l = 1, n
              utu = utu + ufull(l,i) * ufull(l,j)
            end do
            diff = diff + ( a(k) - utu )**2
          end do
        end do

        diff = sqrt ( diff )

        write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of CHOLESKY.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )
      integer nn_max
      parameter ( nn_max = (n_max*(n_max+1))/2 )

      double precision a(nn_max)
      double precision diff
      integer i
      integer ifault
      integer j
      integer k
      integer l
      integer n
      integer nn
      integer nullty
      double precision u(nn_max)
      double precision ufull(n_max,n_max)
      double precision utu

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  CHOLESKY computes the Cholesky factor'
      write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
      write ( *, '(a)' ) '  A compressed storage format is used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
      write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this particular matrix, we expect the'
      write ( *, '(a)' ) '  errors to grow rapidly.'

      do n = 1, n_max

        nn = ( n * ( n + 1 ) ) / 2
!
!  Set A to the Hilbert matrix.
!
        k = 0
        do i = 1, n
          do j = 1, i
            k = k + 1
            a(k) = 1.0D+00 / dble ( i + j - 1 )
          end do
        end do

        call cholesky ( a, n, nn, u, nullty, ifault )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

        k = 0
        do j = 1, n
          do i = 1, j
            k = k + 1
            ufull(i,j) = u(k)
          end do
          do i = j + 1, n
            ufull(i,j) = 0.0D+00
          end do
        end do
!
!  Compute U' * U and compare to A.
!
        k = 0
        diff = 0.0D+00
        do i = 1, n
          do j = 1, i
            k = k + 1
            utu = 0.0D+00
            do l = 1, n
              utu = utu + ufull(l,i) * ufull(l,j)
            end do
            diff = diff + ( a(k) - utu )**2
          end do
        end do

        diff = sqrt ( diff )

        write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 demonstrates the use of SUBCHL.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )
      integer nn_max
      parameter ( nn_max = (n_max*(n_max+1))/2 )

      double precision a(nn_max)
      integer b(n_max)
      double precision det
      double precision diff
      integer i
      integer ifault
      integer j
      integer k
      integer l
      integer n
      integer nullty
      double precision u(nn_max)
      double precision ufull(n_max,n_max)
      double precision utu

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  SUBCHL computes the Cholesky factor'
      write ( *, '(a)' ) '  of a submatrix '
      write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
      write ( *, '(a)' ) '  A compressed storage format is used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
      write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this particular matrix, we expect the'
      write ( *, '(a)' ) '  errors to grow rapidly.'
!
!  Set A to the N_MAX order Hilbert matrix.
!
      k = 0
      do i = 1, n_max
        do j = 1, i
          k = k + 1
          a(k) = 1.0D+00 / dble ( i + j - 1 )
        end do
      end do

      do n = 1, n_max

        do i = 1, n
          b(i) = i
        end do

        call subchl ( a, b, n, u, nullty, ifault, nn_max, det )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty
        write ( *, '(a,g14.6)' ) '  Matrix determinant DET = ', det

        k = 0
        do j = 1, n
          do i = 1, j
            k = k + 1
            ufull(i,j) = u(k)
          end do
          do i = j + 1, n
            ufull(i,j) = 0.0D+00
          end do
        end do
!
!  Compute U' * U and compare to A.
!
        k = 0
        diff = 0.0D+00
        do i = 1, n
          do j = 1, i
            k = k + 1
            utu = 0.0D+00
            do l = 1, n
              utu = utu + ufull(l,i) * ufull(l,j)
            end do
            diff = diff + ( a(k) - utu )**2
          end do
        end do

        diff = sqrt ( diff )

        write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
      end do

      return
      end

