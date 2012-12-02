      program main

c*********************************************************************72
c
cc MAIN is the main program for NAPACK_PRB.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NAPACK_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the NAPACK library.'

      call test_addchg ( )
      call test_bal ( )
      call test_con ( )
      call test_kcon ( )
      call test_kdet ( )
      call test_kfact ( )
      call test_ksolve ( )
      call test_ktrans ( )
      call test_kvert ( )
      call test_mgrid ( )
      call test_pack ( )
      call test_rpack ( )
      call test_scon ( )
      call test_sdet ( )
      call test_smult ( )
      call test_sort ( )
      call test_sort2 ( )

      call test02 ( )
      call test03 ( )
      call test08 ( )
      call test09 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NAPACK:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test_addchg ( )

c*********************************************************************72
c
cc TEST_ADDCHG tests ADDCHG.
c
c  Modified:
c
c    14 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      real chg(n)
      real dif
      integer i
      integer it
      real size
      real xnew(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_ADDCHG'
      write ( *, '(a)' ) '  ADDCHG monitors the change in a vector'
      write ( *, '(a)' ) '  being updated in an iteration.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        IT    SIZE(X)         CHANGE(X)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        xnew(i) = real ( i )
      end do

      it = 0

10    continue

        it = it + 1

        do i = 1, n
          chg(i) = 0.5 * ( ( real ( i ) / xnew(i) ) - xnew(i) )
        end do

        call addchg ( dif, size, xnew, chg, n )

        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) it, size, dif

      if ( 0.001E+00 * size .lt. dif ) then
        go to 10
      end if

      call r4vec_print ( n, xnew, '  Accepted iterate:' )

      return
      end
      subroutine test_bal ( )

c*********************************************************************72
c
cc TEST_BAL tests BAL.
c
c  Modified:
c
c    19 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      real a(n,n)
      real d(n)
      integer i
      integer j
      real w(n,2)

      do i = 1, n
        do j = 1, n
          a(i,j) = real ( i + j - 1 )
        end do
      end do

      a(1,3) = 100.0 * a(1,3)
      a(2,3) = 100.0 * a(2,3)
      a(4,3) = 100.0 * a(4,3)
      a(5,3) = 100.0 * a(5,3)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BAL'
      write ( *, '(a)' ) '  BAL balances a matrix.'

      call r4mat_print ( n, n, a, '  Matrix A:' );

      call bal ( a, n, n, d, w )

      call r4mat_print ( n, n, a, '  After balancing:' )

      return
      end
      subroutine test_con ( )

c*********************************************************************72
c
cc TEST_CON tests CON.
c
c  Modified:
c
c    15 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real a_con_l1
      real a_inv_norm_l1
      real a_norm_l1
      real b(n,n)
      real b_hager(3+n*(n+1))
      real c(n,n)
      real con
      integer i
      integer j
      integer k
      real r4mat_norm_l1
      real w(n)

      save a

c     data a /
c    &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
c    &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
c    &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
c    &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /
c
c  This matrix needs no pivoting.
c
c  The factors are
c
c    1    0     0   0     7.0 10.00  8.00  7.00
c    0.85 1     0   0     0   -0.57  3.14  3.00
c    0.71 0.25  1   0     0    0     2.50  4.25
c    0.71 0.25 -0.2 1     0    0     0     0.10
c
      data a /
     &  7.0E+00,  6.0E+00,  5.0E+00,  5.0E+00,
     & 10.0E+00,  8.0E+00,  7.0E+00,  7.0E+00,
     &  8.0E+00, 10.0E+00,  9.0E+00,  6.0E+00,
     &  7.0E+00,  9.0E+00, 10.0E+00,  5.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_CON'
      write ( *, '(a)' ) 
     &  '  CON computes the condition number of a matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )

      a_norm_l1 = r4mat_norm_l1 ( n, n, a )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          b_hager(k) = a(i,j)
        end do
      end do

      call fact ( b_hager, n, n )
c
c  Request condition number of A.
c
      a_con_l1 = con ( b_hager, w )
c
c  Retrieve the L/U factors which are packed inside B_HAGER.
c
      k = 4
      do j = 1, n
        do i = 1, n
          k = k + 1
          b(i,j) = b_hager(k)
        end do
        k = k + 1
      end do

      call r4mat_print ( n, n, b, '  LU factor information B: ' )
c
c  Inverse(A) = -41, -17, 10, 68
c                25, 10, -6, -41
c                10,  5, -3, -17,
c                -6, -3, 2, 10
c
      call vert ( a, n, n, w )

      do j = 1, n
        do i = 1, n
          b(i,j) = a(i,j)
        end do
      end do

      call r4mat_print ( n, n, b, '  Inverse matrix' )

      a_inv_norm_l1 = r4mat_norm_l1 ( n, n, b )
 
      write ( *, '(a)') ' '
      write ( *, '(a,g14.6)' ) '  L1 norm of A = ', a_norm_l1
      write ( *, '(a,g14.6)' ) '  L1 norm of inv(A) = ', a_inv_norm_l1
      write ( *, '(a,g14.6)' ) 
     &  '  L1 condition number estimate is ', a_con_l1
      write ( *, '(a,g14.6)' ) 
     &  '  ||A|| * ||inv(A)|| = ', a_norm_l1 * a_inv_norm_l1

      return
      end
      subroutine test_kcon ( )

c*********************************************************************72
c
cc TEST_KCON tests KCON.
c
c  Modified:
c
c    15 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real a_con_l1
      real b(n*n+2*n+2)
      real kcon
      integer i
      integer j
      integer k
      real w(n)

      save a

      data a /
     &  7.0E+00,  6.0E+00,  5.0E+00,  5.0E+00,
     & 10.0E+00,  8.0E+00,  7.0E+00,  7.0E+00,
     &  8.0E+00, 10.0E+00,  9.0E+00,  6.0E+00,
     &  7.0E+00,  9.0E+00, 10.0E+00,  5.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KCON'
      write ( *, '(a)' ) 
     &  '  KCON computes the condition number of a matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do

      call kfact ( b, n, n )
c
c  Request condition number of A.
c
      a_con_l1 = kcon ( b, w )
 
      write ( *, '(a)') ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L1 condition number estimate is ', a_con_l1

      return
      end
      subroutine test_kdet ( )

c*********************************************************************72
c
cc TEST_KDET tests KDET.
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real b(n*n+2*n+2)
      real det
      integer i
      integer iexp
      integer j
      integer k
      real kdet

      save a

      data a /
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
     &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
     &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KDET'
      write ( *, '(a)' ) 
     &  '  KDET computes the determinant of a general matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do

      call kfact ( b, n, n )

      det = kdet ( iexp, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i4)' ) 
     &  '  Determinant = ', det, ' * 10 ^ ', iexp
      write ( *, '(a)' ) '  Correct answer is 1.'

      return
      end
      subroutine test_kfact ( )

c*********************************************************************72
c
cc TEST_KFACT tests KFACT.
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real b(n*n+2*n+2)
      integer i
      integer j
      integer k

      save a

      data a /
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
     &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
     &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KFACT'
      write ( *, '(a)' ) '  KFACT factors a general matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do

      call kfact ( b, n, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L1 norm of matrix A is ', b(3)

      call r4vec_print ( n * n + 2 * n + 2, b, 
     &  '  The factor information:' )

      return
      end
      subroutine test_ksolve ( )

c*********************************************************************72
c
cc TEST_KSOLVE tests KSOLVE.
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real amat(n,n)
      real a(n*n+2*n+2)
      real b(n)
      integer i
      integer j
      integer k
      real x(n)

      save amat
      save b

      data amat /
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
     &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
     &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /

      data b /
     &  57.0, 79.0, 88.0, 86.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KSOLVE'
      write ( *, '(a)' ) '  KSOLVE solves A*x=b.'

      call r4mat_print ( n, n, amat, '  Matrix A: ' )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          a(k) = amat(i,j)
        end do
      end do

      call kfact ( a, n, n )

      call ksolve ( x, a, b )

      call r4vec_print ( n, x, '  Computed solution:' )

      return
      end
      subroutine test_ktrans ( )

c*********************************************************************72
c
cc TEST_KTRANS tests KTRANS.
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real amat(n,n)
      real a(n*n+2*n+2)
      real b(n)
      integer i
      integer j
      integer k
      real x(n)

      save amat
      save b

      data amat /
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
     &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
     &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /

      data b /
     &  57.0, 79.0, 88.0, 86.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KTRANS'
      write ( *, '(a)' ) '  KTRANS solves A''*x=b.'

      call r4mat_print ( n, n, amat, '  Matrix A: ' )

      k = 0
      do j = 1, n
        do i = 1, n
          k = k + 1
          a(k) = amat(i,j)
        end do
      end do

      call kfact ( a, n, n )

      call ktrans ( x, a, b )

      call r4vec_print ( n, x, '  Computed solution:' )

      return
      end
      subroutine test_kvert ( )

c*********************************************************************72
c
cc TEST_KVERT tests KVERT.
c
c  Modified:
c
c    15 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real b(n,n)
      real c(n,n)
      integer i
      integer j
      integer k
      real w(n,2)

      save a

      data a /
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00,
     &  7.0E+00, 10.0E+00,  8.0E+00,  7.0E+00,
     &  6.0E+00,  8.0E+00, 10.0E+00,  9.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00, 10.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_KVERT'
      write ( *, '(a)' ) '  KVERT inverts a general matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )

      do j = 1, n
        do i = 1, n
          b(i,j) = a(i,j)
        end do
      end do

      call kvert ( b, n, n, w )

      call r4mat_print ( n, n, b, '  Computed inverse B: ' )

      do j = 1, n
        do i = 1, n
          c(i,j) = 0.0E+00
          do k = 1, n
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      call r4mat_print ( n, n, c, '  Product A*B: ' )

      return
      end
      subroutine test_mgrid ( )

c*********************************************************************72
c
cc TEST_MGRID tests MGRID.
c
c  Modified:
c
c    21 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MGRID:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MGRID multigrid library.'

      do k = 2, 8
        call mgrid ( k )
      end do

      return
      end
      subroutine test_pack ( )

c*********************************************************************72
c
cc TEST_PACK tests PACK and UNPACK
c
c  Modified:
c
c    15 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer la
      parameter ( la = n + 2 )

      real a(la,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = real ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_PACK'
      write ( *, '(a)' ) 
     &  '  PACK packs a square matrix into contiguous storage.'
      write ( *, '(a)' ) '  UNPACK restores it.'

      call r4mat_print ( la, n, a, '  5X5 matrix in 7x5 storage' )

      call pack ( a, la, n )

      call r4mat_print ( la, n, a, '  After packing:' )

      call unpack ( a, la, n )

      call r4mat_print ( la, n, a, '  After unpacking:' )

      return
      end
      subroutine test_rpack ( )

c*********************************************************************72
c
cc TEST_RPACK tests RPACK
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 5 )
      integer la
      parameter ( la = m + 3 )

      real a(la,n)
      integer i
      integer j

      do i = 1, m
        do j = 1, n
          a(i,j) = real ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_RPACK'
      write ( *, '(a)' ) 
     &  '  RPACK packs a rectangular matrix into contiguous storage.'
      write ( *, '(a)' ) '  RUNPACK reverses the process.'

      call r4mat_print ( la, n, a, '  3X5 matrix in 6x5 storage' )

      call rpack ( a, la, m, n )

      call r4mat_print ( la, n, a, '  After packing:' )

      call runpack ( a, la, m, n )

      call r4mat_print ( la, n, a, '  After unpacking:' )

      return
      end
      subroutine test_scon ( )

c*********************************************************************72
c
cc TEST_SCON tests SCON.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real a_con_l1
      real b(3+(n*(n+1))/2)
      integer i
      integer j
      integer k
      real scon
      real w(n)

      save a

      data a /
     &  7.0E+00,  6.0E+00,  5.0E+00,  5.0E+00,
     &  6.0E+00,  8.0E+00,  7.0E+00,  7.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00,  6.0E+00,
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SCON'
      write ( *, '(a)' ) 
     &  '  SCON computes the condition number of a symmetric matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )
c
c  Copy the rows starting at the diagaonl.
c
      k = 0
      do j = 1, n
        do i = j, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do
c
c  Factor the matrix.
c
      call sfact ( b, n, w )
c
c  Request condition number.
c
      a_con_l1 = scon ( b, w )
 
      write ( *, '(a)') ' '
      write ( *, '(a,g14.6)' ) 
     &  '  L1 condition number estimate is ', a_con_l1

      return
      end
      subroutine test_sdet ( )

c*********************************************************************72
c
cc TEST_SDET tests SDET.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real b(3+(n*(n+1))/2)
      real det
      integer i
      integer iexp
      integer j
      integer k
      real sdet
      real w(n)

      save a

      data a /
     &  7.0E+00,  6.0E+00,  5.0E+00,  5.0E+00,
     &  6.0E+00,  8.0E+00,  7.0E+00,  7.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00,  6.0E+00,
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SDET'
      write ( *, '(a)' ) 
     &  '  SDET computes the determinant of a symmetric matrix.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )
c
c  Copy the rows starting at the diagaonl.
c
      k = 0
      do j = 1, n
        do i = j, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do

      call sfact ( b, n, w )

      det = sdet ( iexp, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,i4)' ) 
     &  '  Determinant = ', det, ' * 10 ^ ', iexp

      return
      end
      subroutine test_smult ( )

c*********************************************************************72
c
cc TEST_SMULT tests SMULT.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      real a(n,n)
      real b(3+(n*(n+1))/2)
      integer i
      integer j
      integer k
      real x(n)
      real y(n)

      save a

      data a /
     &  7.0E+00,  6.0E+00,  5.0E+00,  5.0E+00,
     &  6.0E+00,  8.0E+00,  7.0E+00,  7.0E+00,
     &  5.0E+00,  7.0E+00,  9.0E+00,  6.0E+00,
     &  5.0E+00,  7.0E+00,  6.0E+00,  5.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SMULT'
      write ( *, '(a)' ) 
     &  '  SMULT computes y=A*x for a symmetric matrix A.'

      call r4mat_print ( n, n, a, '  Matrix A: ' )
c
c  Copy the rows starting at the diagaonl.
c
      k = 0
      do j = 1, n
        do i = j, n
          k = k + 1
          b(k) = a(i,j)
        end do
      end do
c
c  Compute the product.
c
      do i = 1, n
        x(i) = real ( i )
      end do

      call r4vec_print ( n, x, '  X:' )

      call smult ( y, x, b, n )

      call r4vec_print ( n, y, '  Y = A*X:' )

      return
      end
      subroutine test_sort ( )

c*********************************************************************72
c
cc TEST_SORT tests SORT.
c
c  Modified:
c
c    18 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer seed
      real x(n)
      real y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SORT'
      write ( *, '(a)' ) '  SORT sorts a real array.'

      seed = 123456789

      call r4vec_uniform_01 ( n, seed, x )

      call r4vec_print ( n, x, '  X before sorting:' )

      call sort ( x, y, n )

      call r4vec_print ( n, x, '  X after sorting:' )

      return
      end
      subroutine test_sort2 ( )

c*********************************************************************72
c
cc TEST_SORT2 tests SORT2.
c
c  Modified:
c
c    19 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer seed
      integer w(n)
      real x(n)
      integer y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SORT2'
      write ( *, '(a)' ) '  SORT2 index sorts a real array.'

      seed = 123456789

      call r4vec_uniform_01 ( n, seed, x )

      call r4vec_print ( n, x, '  X before indexing:' )

      call sort2 ( x, y, w, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I  Y(I)     X(Y(I))'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, y(i), x(y(i))
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc test02 approximates the semiconductor equation at evenly spaced points               
c
      parameter (np1=11)
      real new(np1),old(np1),d,dif,s,size,t
      integer i,n,nm1

      write ( *, * ) ' '
      write ( *, * ) 'TEST02'
      write ( *, * ) '  Approximate the semiconductor equation.'
      write ( *, * ) ' '
c
c  set number of intervals
c
      n=np1-1
c
c  preprocess and initialize variables
c
20    nm1 = n - 1
      d = 1./float(n)
      s = d*d
      do 30 i = 1,n
           old(i) = 0.
30    continue
c
c  compute next iteration
c
40    new(1) = .5*(old(2) + s*exp(-old(1)))
      if ( n .le. 2 ) goto 60
      do 50 i = 2,nm1
           new(i) = .5*(old(i-1) + old(i+1) + s*exp(-old(i)))
50    continue
c
c  update and test error
c
60    call update(dif,size,new,old,nm1)
      call stopit(dif,size,4,5000)
      if ( dif .gt. 0 ) goto 40
      call whatis(dif,size)
c
c  print the solution
c
      write(6,70) n
70    format(i10,'  intervals',/'      t             x(t)'/)
      do i = 1,nm1
           t = i*d
           write(6,80) t,new(i)
80         format(f10.4,f15.4)
      end do

      write(6,100)
100   format(//)
      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc test03 approximates the semiconductor equation at evenly spaced points
c  ( chapter 2 algorithm )                  
c
      real  new(50),old(50),d(50),u(50)
      real  s,t,dif,size

      write ( *, * ) ' '
      write ( *, * ) 'Test03'
      write ( *, * ) '  Semiconductor equation, another algorithm.'
      write ( *, * ) ' '
c
c  set number of intervals
c
      n=10
c            
c  preprocess and initialize variables
c            
30    nm1 = n - 1
      s = n*n
      do i = 1,nm1
        old(i) = 0.
        d(i) = s + s
        u(i) = -s
      end do

      call tfact(u,d,u,nm1)
c
c  compute next iteration
c
50    do 60 i = 1,nm1
           new(i) = exp(-old(i))
60    continue
      call tsolve(new,u,d,u,new)
c
c                   *** update and test error ***
c                   -
      call update(dif,size,new,old,nm1)
      call stopit(dif,size,4,100)
      if ( dif .gt. 0 ) goto 50
      call whatis(dif,size)
c                     --
c                     *** print the answer ***
c
      write(6,70) n
70    format(i11,'  intervals',/'      t               x(t)'/)
      do 90 i = 1,nm1
           t = i/float(n)
           write(6,80) t,new(i)
80         format(f10.5,f20.8)
90    continue
      write(6,100)
100   format(//)
      return
      end
      subroutine primat(a,lda,m,n)

c*********************************************************************72
c
cc primat prints a matrix.
c
      dimension a(lda,n)

      do 10 i=1,m
        write(*,'(1x,6g11.4)')(a(i,j),j=1,n)
10      continue
      return
      end
      subroutine mulmat(a,lda,ma,na,b,ldb,mb,nb,c,ldc,mc,nc)

c*********************************************************************72
c
cc mulmat multiplies two matrices.
c
      dimension a(lda,na)
      dimension b(ldb,nb)
      dimension c(ldc,nc)

      do 30 i=1,ma
        do 20 j=1,nb
          c(i,j)=0.0
          do 10 k=1,na
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
10          continue
20        continue
30      continue
      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests BASIS.
c
      parameter (la=8)
      parameter (lb=la)
      parameter (m=4)
      parameter (n=3)

      real a(la,n)
      real b(lb,n)
      real cutoff
      integer nvec

      a(1,1)=1.0
      a(2,1)=2.0
      a(3,1)=3.0
      a(4,1)=-3.0

      a(1,2)=2.0
      a(2,2)=4.0
      a(3,2)=3.0
      a(4,2)=0.0

      a(1,3)=4.0
      a(2,3)=8.0
      a(3,3)=2.0
      a(4,3)=8.0

      write ( *, * ) ' '
      write ( *, * ) 'TEST08'
      write ( *, * ) '  Test BASIS'
      write ( *, * ) ' '

      call qr ( a, la, m, n )
      cutoff = 0.001E+00
      call basis ( b, lb, nvec, a, cutoff )

      write ( *, * ) 'nvec=', nvec
      call primat ( b, lb, m, n )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc test09 tests QR and under.
c
      parameter (la=8)
      parameter (m=4)
      parameter (n=3)

      dimension a(la,n)
      dimension x(n)
      dimension b(m)
c
c  in order to solve the underdetermined system a*x=b, we have to
c  use the transpose of the matrix we're interested in, and we
c  have to call qr, then under.
c
      a(1,1)=1.0
      a(2,1)=2.0
      a(3,1)=3.0
      a(4,1)=-3.0

      a(1,2)=2.0
      a(2,2)=4.0
      a(3,2)=3.0
      a(4,2)=0.0

      a(1,3)=4.0
      a(2,3)=8.0
      a(3,3)=2.0
      a(4,3)=8.0

      b(1)=1.0
      b(2)=2.0
      b(3)=1.0
      b(4)=1.0

      write ( *, * ) ' '
      write ( *, * ) 'Test09'
      write ( *, * ) '  test qr and under'
      write ( *, * ) ' '

      call qr(a,la,m,n)
      call under(x,a,b)

      call r4vec_print ( n, x, '  Computed solution:' )

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc test11 solves a large system with kfact/kdet/ksolve.
c
c  the problem for n=100 is just an enlarged version of the
c  problem for n=5, which is:
c
c  matrix a is ( 2 -1  0  0  0)    right hand side b is  (1)     
c              (-1  2 -1  0  0)                          (0)
c              ( 0 -1  2 -1  0)                          (0)
c              ( 0  0 -1  2 -1)                          (0)
c              ( 0  0  0 -1  2)                          (1)
c
c  solution is   (1)
c                (1)
c                (1)
c                (1)
c                (1)
c
      parameter (n=100)
c
c  in order to guarantee the space required by kfact, we set
c  lda to n+3.
c
      parameter (lda=n+3)

      dimension a(lda,n)
      dimension b(n)
      real      kcon
      real      kdet
      dimension work(n)

      write ( *, * ) ' '
      write ( *, * ) 'Test11'
      write ( *, '(a)' ) '  For a general matrix,'
      write ( *, '(a)' ) '  KFACT computes the LU factorization;'
      write ( *, '(a)' ) '  KDET computes the determinant.'
      write ( *, '(a)' ) '  KSOLVE solves a linear system.'
      write ( *, * ) ' '
c
c  assign values to matrix a and right hand side b.
c
      b(1) = 1.0E+00
      do i = 2, n - 1
        b(i) = 0.0E+00
      end do
      b(n) = 1.0E+00

      do i=1,n
        do j=1,n
          a(i,j)=0.0
        end do
        a(i,i)=2.0
        if(i.gt.1)a(i,i-1)=-1.0
        if(i.lt.n)a(i,i+1)=-1.0
      end do
c
c  call kfact to factor matrix a.
c
      call kfact(a,lda,n)

      det=kdet(iexp,a)
      write ( *, * ) 'determinant=',det,' * 10 ** ',iexp
      conk = kcon ( a, work )
      write ( *, * ) '  Condition number=',conk
c
c  call ksolve to solve linear system a*x=b once a has been factored.
c
      call ksolve(b,a,b)
      write ( *, * ) ' '
      write ( *, * ) 'first 5 entries of solution (all should be 1)'
      write(*,'(1x,5g14.6)')(b(j),j=1,5)
      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc test12 solves a large tridiagonal system with pfact/pdet/psolve.
c
c  the problem for n=100 is just an enlarged version of the
c  problem for n=5, which is:
c
c  matrix a is ( 2 -1  0  0  0)    right hand side b is  (1)     
c              (-1  2 -1  0  0)                          (0)
c              ( 0 -1  2 -1  0)                          (0)
c              ( 0  0 -1  2 -1)                          (0)
c              ( 0  0  0 -1  2)                          (1)
c
c  solution is   (1)
c                (1)
c                (1)
c                (1)
c                (1)
c
      parameter (n=100)
c
c  in order to guarantee the space required by pfact, we set
c  lda to 5.
c
      parameter (lda=5)

      real a(lda,n)
      real b(n)
      real pcon
      real pdet
      real work(n)

      write ( *, * ) ' '
      write ( *, * ) 'TEST12'
      write ( *, * ) '  PFACT factors a triadiagonal matrix;'
      write ( *, * ) '  PSOLVE solves linear systems;'
      write ( *, * ) '  PDET computes the determinant.'
      write(*,1110)
      write ( *, * ) ' '
c
c  assign values to matrix a and right hand side b.
c
      do i=1,n
        b(i)=0.0
      end do

      do 6 i=1,lda
        do 7 j=1,n
          a(i,j)=0.0
7         continue
6       continue

      do 10 j=1,n
        if(j.lt.n)a(1,j)=-1.0
        a(2,j)=2.0
        if(j.lt.n)a(3,j)=-1.0
10      continue
      b(1)=1.0
      b(n)=1.0
c
c  call pfact to factor matrix a.
c
      call pfact(a,lda,n)

      det=pdet(iexp,a)

      conp = pcon ( a, work )

      write ( *, * ) ' '
      write ( *, * ) '  determinant=',det,' * 10 ** ',iexp
      write ( *, * ) '  condition number=',conp
c
c  call psolve to solve linear system a*x=b once a has been factored.
c
      call psolve(b,a,b)
      write ( *, * ) ' '
      write(*,1020)
      write(*,'(1x,5g14.6)')(b(j),j=1,5)
      return
1020  format(' first 5 entries of solution: (all should be 1)')
1110  format(' demonstrate use of napack routines pfact/pdet/psolve.')
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc test13 tests sfact, sdet, ssolve
c
c  the actual matrix is
c
c  2 4 6
c  4 7 8
c  6 8 5
c
      parameter (na=9)
      parameter (n=3)

      dimension a(na)
      dimension b(n)
      dimension work(n)

      write ( *, * ) ' '
      write ( *, * ) 'Test13'
      write ( *, '(a)' ) '  For a symmetric matrix,'
      write ( *, '(a)' ) '  SFACT computes the LU factorization;'
      write ( *, '(a)' ) '  SDET computes the determinant.'
      write ( *, '(a)' ) '  SSOLVE solves a linear system.'
      write ( *, * ) ' '
c
c  symmetric matrix is stored by columns, from the diagonal down.
c

      a(1)=2.0
      a(2)=4.0
      a(3)=6.0

      a(4)=7.0
      a(5)=8.0

      a(6)=5.0

      b(1)=14.0
      b(2)=21.0
      b(3)=20.0

      call sfact(a,n,b)
      det=sdet(iexp,a)
      write(*,*)'determinant is ',det,' *10**',iexp
      cons=scon(a,work)
      write(*,*)'condition number=',cons
      call ssolve(b,a,b)

      write(*,*)'solution is ',b
      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc test14 tests ROOT.
c
      real a
      real b
      external fsin
      external fexp
      real tol

      a=0.5
      b=1.0
      tol=0.00001

      write ( *, * ) ' '
      write ( *, * ) 'Test14'
      write ( *, * ) ' '
      write ( *, * ) 'demonstrate root, finds x so that f(x)=0'
      write ( *, * ) 'between points a and b where f changes sign.'
      write ( *, * ) ' '
      x=root(a,b,tol,fsin)
      write ( *, * ) 'root of x**3-sin(x)=',x
      fx=fsin(x)
      write ( *, * ) 'function value=',fx
      x=root(a,b,tol,fexp)
      write ( *, * ) 'root of x-exp(-x)=',x
      fx=fexp(x)
      write ( *, * ) 'function value=',fx
      return
      end
      real function fsin(x)

c*********************************************************************72
c
cc fsin defines f(x) = x * x * x - sin ( x )
c
      fsin=x*x*x-sin(x)
      return
      end
      real function fexp(x)

c*********************************************************************72
c
cc fexp defines f(x) = x - exp ( - x ).
c
      fexp=x-exp(-x)
      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc test15 tests TFACT and tsolve.
c
c  demonstrate use of tfact, tsolve, power and whatis to carry out
c  the inverse power method on a tridiagonal matrix.
c
      external prod
      common /comm15/ d(12),u(9)
      complex ex,x(9),ey,y(9)

      n=10
      m=n-1
      t=n*n
      do i=1,m
        d(i)=t+t
        u(i)=-t
        x(i)=cmplx(1.0,0.0)
      end do

      write ( *, * ) ' '
      write ( *, * ) 'Test15'
      write ( *, '(a)' ) '  For a tridiagonal matrix,'
      write ( *, '(a)' ) '  TFACT computes the LU factorization;'
      write ( *, '(a)' ) '  TDET computes the determinant.'
      write ( *, '(a)' ) '  TSOLVE solves a linear system.'
      write ( *, * ) ' '
      write ( *, * ) '  We carry out an inverse iteration scheme'
      write ( *, * ) '  for the smallest eigenvalue.'
      write ( *, * ) ' '

      call tfact(u,d,u,m)
      call power(ex,x,ey,y,m,dif,size,5,100,prod)
      call whatis(dif,size)
      e=1.0/ex
      write ( *, * ) ' '
      write ( *, * ) 'eigenvalue=',e
      write ( *, * ) 'eigenvector='
      do i=1,m
        write ( *, * ) i,real(x(i))
      end do

      return
      end
      subroutine prod(p,v)

c*********************************************************************72
c
cc prod
c
      common /comm15/ d(12),u(9)

      real p(*),v(*)

      call tsolve(p,u,d,u,v)

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc test16 tests czero, which finds zeroes of a polynomial
c
c  Discussion:
c
c    The polynomial is 24 - 50 * x + 35 * x^2 - 10 * x^3 + x^4
c    which has roots 1, 2, 3, and 4.
c
      integer nd
      parameter ( nd = 4 )

      integer nw
      parameter ( nw = 6 * nd + 32 )

      integer i
      complex p(nd)
      complex poly
      complex w(nw)
      complex z(nd)

      p(1) = cmplx (  24.0, 0.0 )
      p(2) = cmplx ( -50.0, 0.0 )
      p(3) = cmplx (  35.0, 0.0 )
      p(4) = cmplx ( -10.0, 0.0 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test16'
      write ( *, '(a)' ) '  CZERO finds polynomial roots.'
      write ( *, '(a)' ) '  Test it on:'
      write ( *, '(a)' ) 
     &  '  p(x) = 24 - 50 * x + 35 * x^2 - 10 * x^3 + x^4'
      write ( *, '(a)' ) ' '
      call czero(z,p,nd,w)

      itemp = int ( w(1) )
      write ( *, '(a)' ) 
     &  'czero found ',itemp,' zeroes of the polynomial:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' value of root        value of polynomial'
      write ( *, '(a)' ) ' '

      do i = 1, itemp
        poly=(((z(i)+p(4))*z(i)+p(3))*z(i)+p(2))*z(i)+p(1)
        write ( *, * ) z(i),poly
      end do

      return
      end
