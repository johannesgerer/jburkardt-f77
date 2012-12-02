      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS1_D_PRB.
c
c  Discussion:
c
c    BLAS1_D_PRB tests the BLAS1 routines.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_D_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS1_D library.'

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
!
!  Terminate.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_D_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DASUM.
c
c  Modified:
c
c    28 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ma
      integer na
      integer nx

      parameter ( lda = 6 )
      parameter ( ma = 5 )
      parameter ( na = 4 )
      parameter ( nx = 10 )

      double precision a(lda,na)
      double precision dasum
      integer i
      integer j
      double precision x(nx)

      do i = 1, nx
        x(i) = (-1.0D+00)**i * dble ( 2 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  DASUM adds the absolute values of '
      write ( *, '(a)' ) '  elements of a double precision vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, nx
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  DASUM ( NX, X, 1 ) =   ',
     &  dasum ( nx, x, 1 )
      write ( *, '(a,g14.6)' ) '  DASUM ( NX/2, X, 2 ) = ',
     &  dasum ( nx/2, x, 2 )
      write ( *, '(a,g14.6)' ) '  DASUM ( 2, X, NX/2 ) = ',
     &  dasum ( 2, x, nx/2 )

      do i = 1, lda
        do j = 1, na
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, ma
        do j = 1, na
          a(i,j) = (-1.0D+00)**(i+j) * dble ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate with a matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, ma
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, na )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  DASUM(MA,A(1,2),1) =   ',
     &  dasum ( ma, a(1,2), 1 )
      write ( *, '(a,g14.6)' ) '  DASUM(NA,A(2,1),LDA) = ',
     &  dasum ( na, a(2,1), lda )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DAXPY.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision da
      integer i
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  DAXPY adds a double precision multiple of '
      write ( *, '(a)' ) '  vector X to vector Y.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      da = 1.0D+00
      call daxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = -2.0D+00
      call daxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = +3.0D+00
      call daxpy ( 3, da, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = -4.0D+00
      call daxpy ( 3, da, x, 1, y, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 1, Y, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DCOPY.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(5,5)
      integer i
      integer j
      double precision x(10)
      double precision y(10)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  DCOPY copies one double precision vector'
      write ( *, '(a)' ) '  into another.'

      do i = 1, 10
        x(i) = dble ( i )
      end do

      do i = 1, 10
        y(i) = dble ( 10 * i )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = dble ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      call dcopy ( 5, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, 10
        y(i) = dble ( 10 * i )
      end do

      call dcopy ( 3, x, 2, y, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 3, X, 2, Y, 3 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      call dcopy ( 5, x, 1, a, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, A, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = dble ( 10 * i + j )
        end do
      end do

      call dcopy ( 5, x, 2, a, 5 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 2, A, 5 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DDOT.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer ldb
      integer ldc

      parameter ( n = 5 )
      parameter ( lda = 10 )
      parameter ( ldb = 7 )
      parameter ( ldc = 6 )

      double precision a(lda,lda)
      double precision b(ldb,ldb)
      double precision c(ldc,ldc)
      integer i
      integer j
      double precision ddot
      double precision sum1
      double precision x(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  DDOT computes the dot product of '
      write ( *, '(a)' ) '  double precision vectors.'

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = - dble ( i )
      end do

      do i = 1, n
        do j = 1, n
          a(i,j) = dble ( i + j )
        end do
      end do

      do i = 1, n
        do j = 1, n
          b(i,j) = dble ( i - j )
        end do
      end do
c
c  To compute a simple dot product of two vectors, use a
c  call like this:
c
      sum1 = ddot ( n, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Dot product of X and Y is ', sum1
c
c  To multiply a ROW of a matrix A times a vector X, we need to
c  specify the increment between successive entries of the row of A:
c
      sum1 = ddot ( n, a(2,1), lda, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Product of row 2 of A and X is ', sum1
c
c  Product of a column of A and a vector is simpler:
c
      sum1 = ddot ( n, a(1,2), 1, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Product of column 2 of A and X is ', sum1
c
c  Here's how matrix multiplication, c = a*b, could be done
c  with DDOT:
c
      do i = 1, n
        do j = 1, n
          c(i,j) = ddot ( n, a(i,1), lda, b(1,j), 1 )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix product computed with DDOT:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests DMACH.
c
c  Modified:
c
c    21 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision dmach

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  DMACH computes several machine-dependent'
      write ( *, '(a)' ) '  double precision arithmetic parameters.'

      write ( *, '(a)' ) ' '
      write ( *, * ) '  DMACH(1)  = machine epsilon = ', dmach ( 1 )
      write ( *, * ) '  DMACH(2)  = a tiny value    = ', dmach ( 2 )
      write ( *, * ) '  DMACH(3)  = a huge value    = ', dmach ( 3 )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests DNRM2.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n + 5 )
c
c  These parameters illustrate the fact that matrices are typically
c  dimensioned with more space than the user requires.
c
      double precision a(lda,lda)
      integer i
      integer incx
      integer j
      double precision dnrm2
      double precision sum1
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  DNRM2 computes the Euclidean norm of '
      write ( *, '(a)' ) '  a double precision vector.'
c
c  Compute the euclidean norm of a vector:
c
      do i = 1, n
        x(i) = dble ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do
      incx = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of X is ',
     &  dnrm2 ( n, x, incx )
c
c  Compute the euclidean norm of a row or column of a matrix:
c
      do i = 1, n
        do j = 1, n
          a(i,j) = dble ( i + j )
        end do
      end do

      incx = lda
      sum1 = dnrm2 ( n, a(2,1), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of row 2 of A is ', sum1

      incx = 1
      sum1 = dnrm2 ( n, a(1,2), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of column 2 of A is ', sum1

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests DROT.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      double precision c
      integer i
      double precision s
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( i * i - 12 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  DROT carries out a double precision '
      write ( *, '(a)' ) '  Givens rotation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      c = 0.5D+00
      s = sqrt ( 1.0D+00 - c * c )
      call drot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( i * i - 12 )
      end do

      c = x(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      s = y(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      call drot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests DROTG.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision a
      double precision b
      double precision c
      double precision r8_uniform_01
      double precision r
      double precision s
      double precision sa
      double precision sb
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 5 )
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  DROTG generates a real Givens rotation'
      write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
      write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        a = r8_uniform_01 ( seed )
        b = r8_uniform_01 ( seed )

        sa = a
        sb = b

        call drotg ( sa, sb, c, s )

        r = sa
        z = sb

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a,g14.6)' ) '  A =  ', a,  '  B =  ', b
        write ( *, '(a,g14.6,a,g14.6)' ) '  C =  ', c,  '  S =  ', s
        write ( *, '(a,g14.6,a,g14.6)' ) '  R =  ', r,  '  Z =  ', z
        write ( *, '(a,g14.6)' ) '   C*A+S*B = ',  c * a + s * b
        write ( *, '(a,g14.6)' ) '  -S*A+C*B = ', -s * a + c * b

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests DSCAL.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      double precision da
      integer i
      double precision x(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  DSCAL multiplies a double precision scalar'
      write ( *, '(a)' ) '  times a vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      da = 5.0D+00
      call dscal ( n, da, x, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSCAL ( N, ', da, ', X, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      da = -2.0D+00
      call dscal ( 3, da, x, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSCAL ( 3, ', da, ', X, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests DSWAP.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      integer i
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  DSWAP swaps two vectors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      call dswap ( n, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DSWAP ( N, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      call dswap ( 3, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSWAP ( 3, X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests IDAMAX.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 11 )

      integer i
      integer i1
      integer incx
      integer idamax
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  IDAMAX returns the index of maximum '
      write ( *, '(a)' ) '  magnitude;'

      do i = 1, n
        x(i) = dble ( mod ( 7 * i, 11 ) ) - dble ( n / 2 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do

      incx = 1

      i1 = idamax ( n, x, incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i1

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests IDAMAX, DAXPY and DSCAL.
c
c  Modified:
c
c    15 May 2006
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

      double precision a(lda,n)
      double precision b(n)
      integer i
      integer idamax
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer l
      double precision t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  Use IDAMAX, DAXPY and DSCAL'
      write ( *, '(a)' ) '  in a Gauss elimination routine.'
c
c  Set the matrix.
c
      do i = 1, n
        do j = 1, n

          if ( i == j ) then
            a(i,j) = 2.0D+00
          else if ( i == j + 1 ) then
            a(i,j) = - 1.0D+00
          else if ( i == j - 1 ) then
            a(i,j) = - 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do
c
c  Set the right hand side.
c
      do i = 1, n-1
        b(i) = 0.0D+00
      end do
      b(n) = dble ( n ) + 1.0D+00

      info = 0

      do k = 1, n - 1

        l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
        ipvt(k) = l

        if ( a(l,k) == 0.0D+00 ) then

          info = k

        else

          if ( l /= k ) then
            t = a(l,k)
            a(l,k) = a(k,k)
            a(k,k) = t
          end if

          t = -1.0D+00 / a(k,k)
          call dscal ( n-k, t, a(k+1,k), 1 )

          do j = k+1, n
            t = a(l,j)
            if ( l /= k ) then
              a(l,j) = a(k,j)
              a(k,j) = t
            end if
            call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
          end do

        end if

      end do

      ipvt(n) = n
      if ( a(n,n) == 0.0D+00 ) then
        info = n
      end if

      if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The matrix is singular.'
        return
      end if

      do k = 1, n-1
        l = ipvt(k)
        t = b(l)
        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if
        call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )
      end do

      do k = n, 1, -1
        b(k) = b(k) / a(k,k)
        t = - b(k)
        call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
      end do

      write ( *, '(a,g14.6)' ) ' '
      write ( *, '(a,g14.6)' ) '  First five entries of solution:'
      write ( *, '(a,g14.6)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(j), j = 1, 5 )

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit double precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      r8_uniform_01
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
c    Paul Bratley, Bennett Fox, Linus Schrage,
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
c    Philip Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
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
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

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
