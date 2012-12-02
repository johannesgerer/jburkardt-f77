      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS1_S_PRB.
c
c  Discussion:
c
c    BLAS1_S_PRB tests the BLAS1 routines.
c
c  Modified:
c
c    14 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_S_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS1_S library.'

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
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_S_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ISAMAX.
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
      integer isamax
      real x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  ISAMAX returns the index of maximum '
      write ( *, '(a)' ) '  magnitude;'

      do i = 1, n
        x(i) = real ( mod ( 7 * i, 11 ) ) - real ( n / 2 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do

      incx = 1

      i1 = isamax ( n, x, incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i1

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests ISAMAX, SAXPY and SSCAL.
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

      real a(lda,n)
      real b(n)
      integer i
      integer info
      integer ipvt(n)
      integer isamax
      integer j
      integer k
      integer l
      real t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Use ISAMAX, SAXPY and SSCAL'
      write ( *, '(a)' ) '  in a Gauss elimination routine.'
c
c  Set the matrix.
c
      do i = 1, n
        do j = 1, n

          if ( i == j ) then
            a(i,j) = 2.0E+00
          else if ( i == j + 1 ) then
            a(i,j) = - 1.0E+00
          else if ( i == j - 1 ) then
            a(i,j) = - 1.0E+00
          else
            a(i,j) = 0.0E+00
          end if

        end do
      end do
c
c  Set the right hand side.
c
      do i = 1, n-1
        b(i) = 0.0E+00
      end do
      b(n) = real ( n ) + 1.0E+00

      info = 0

      do k = 1, n - 1

        l = isamax ( n-k+1, a(k,k), 1 ) + k - 1
        ipvt(k) = l

        if ( a(l,k) == 0.0E+00 ) then

          info = k

        else

          if ( l /= k ) then
            t = a(l,k)
            a(l,k) = a(k,k)
            a(k,k) = t
          end if

          t = -1.0E+00 / a(k,k)
          call sscal ( n-k, t, a(k+1,k), 1 )

          do j = k+1, n
            t = a(l,j)
            if ( l /= k ) then
              a(l,j) = a(k,j)
              a(k,j) = t
            end if
            call saxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
          end do

        end if

      end do

      ipvt(n) = n
      if ( a(n,n) == 0.0E+00 ) then
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
        call saxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )
      end do

      do k = n, 1, -1
        b(k) = b(k) / a(k,k)
        t = - b(k)
        call saxpy ( k-1, t, a(1,k), 1, b(1), 1 )
      end do

      write ( *, '(a,g14.6)' ) ' '
      write ( *, '(a,g14.6)' ) '  First five entries of solution:'
      write ( *, '(a,g14.6)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(j), j = 1, 5 )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests SASUM.
c
c  Modified:
c
c    22 February 2006
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

      real a(lda,na)
      integer i
      integer j
      real sasum
      real x(nx)

      do i = 1, nx
        x(i) = ( -1.0E+00 )**i * real ( 2 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  SASUM adds the absolute values '
      write ( *, '(a)' ) '  of elements of a vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, nx
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  SASUM ( NX, X, 1 ) =   ',
     &   sasum ( nx, x, 1 )
      write ( *, '(a,g14.6)' ) '  SASUM ( NX/2, X, 2 ) = ',
     &   sasum ( nx/2, x, 2 )
      write ( *, '(a,g14.6)' ) '  SASUM ( 2, X, NX/2 ) = ',
     &   sasum ( 2, x, nx/2 )

      do i = 1, ma
        do j = 1, na
          a(i,j) = ( -1.0E+00 )**( i + j ) * real ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate with a matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, ma
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, na )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  SASUM(MA,A(1,2),1) =   ',
     &  sasum ( ma, a(1,2), 1 )
      write ( *, '(a,g14.6)' ) '  SASUM(NA,A(2,1),LDA) = ',
     &  sasum ( na, a(2,1), lda )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests SAXPY.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      real da
      integer i
      real x(n)
      real y(n)

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  SAXPY adds a multiple of vector X to vector Y.'
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

      da = 1.0E+00
      call saxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      da = -2.0E+00
      call saxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      da = +3.0E+00
      call saxpy ( 3, da, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SAXPY ( 3, ', da, ', X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      da = -4.0E+00
      call saxpy ( 3, da, x, 1, y, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SAXPY ( 3, ', da, ', X, 1, Y, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests SCOPY.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a(5,5)
      integer i
      integer j
      real x(10)
      real y(10)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  SCOPY copies one vector into another.'

      do i = 1, 10
        x(i) = real ( i )
      end do

      do i = 1, 10
        y(i) = real ( 10 * i )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = real ( 10 * i + j )
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
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5)
      end do

      call scopy ( 5, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, 10
        y(i) = real ( 10 * i )
      end do

      call scopy ( 3, x, 2, y, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SCOPY ( 3, X, 2, Y, 3 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      call scopy ( 5, x, 1, a, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 1, A, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = real ( 10 * i + j )
        end do
      end do

      call scopy ( 5, x, 2, a, 5 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SCOPY ( 5, X, 2, A, 5 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests SDOT.
c
c  Modified:
c
c    22 February 2006
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

      real a(lda,lda)
      real b(ldb,ldb)
      real c(ldc,ldc)
      integer i
      integer j
      real sdot
      real sum1
      real x(n)
      real y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  SDOT computes the dot product of vectors.'

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = - real ( i )
      end do

      do i = 1, n
        do j = 1, n
          a(i,j) = real ( i + j )
        end do
      end do

      do i = 1, n
        do j = 1, n
          b(i,j) = real ( i - j )
        end do
      end do
c
c  To compute a simple dot product of two vectors, use a
c  call like this:
c
      sum1 = sdot ( n, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Dot product of X and Y is ', sum1
c
c  To multiply a ROW of a matrix A times a vector X, we need to
c  specify the increment between successive entries of the row of A:
c
      sum1 = sdot ( n, a(2,1), lda, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Product of row 2 of A and X is ', sum1
c
c  Product of a column of A and a vector is simpler:
c
      sum1 = sdot ( n, a(1,2), 1, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Product of column 2 of A and X is ', sum1
c
c  Here's how matrix multiplication, c = a*b, could be done
c  with SDOT:
c
      do i = 1, n
        do j = 1, n
          c(i,j) = sdot ( n, a(i,1), lda, b(1,j), 1 )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix product computed with SDOT:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests SMACH.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real smach

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  SMACH computes several machine arithmetic parameters.'

      write ( *, '(a)' ) ' '
      write ( *, * ) '  SMACH(1)  = machine epsilon = ', smach ( 1 )
      write ( *, * ) '  SMACH(2)  = a tiny value    = ', smach ( 2 )
      write ( *, * ) '  SMACH(3)  = a huge value    = ', smach ( 3 )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests SNRM2.
c
c  Modified:
c
c    22 February 2006
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
      real a(lda,lda)
      integer i
      integer incx
      integer j
      real snrm2
      real sum1
      real x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  SNRM2 computes the Euclidean norm '
      write ( *, '(a)' ) '  of a vector.'
c
c  Compute the euclidean norm of a vector:
c
      do i = 1, n
        x(i) = real ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do
      incx = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  The 2-norm of X is ', snrm2 ( n, x, incx )
c
c  Compute the euclidean norm of a row or column of a matrix:
c
      do i = 1, n
        do j = 1, n
          a(i,j) = real ( i + j )
        end do
      end do

      incx = lda
      sum1 = snrm2 ( n, a(2,1), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  The 2-norm of row 2 of A is ', sum1

      incx = 1
      sum1 = snrm2 ( n, a(1,2), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  The 2-norm of column 2 of A is ', sum1

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests SROT.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      real c
      integer i
      real s
      real x(n)
      real y(n)

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = real ( i * i - 12 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  SROT carries out a Givens rotation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      c = 0.5E+00
      s = sqrt ( 1.0E+00 - c * c )
      call srot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  SROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = real ( i * i - 12 )
      end do

      c = x(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      s = y(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      call srot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  SROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests SROTG.
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

      real a
      real b
      real c
      real r
      real r4_uniform_01
      real s
      real sa
      real sb
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 5 )
      real z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  SROTG generates a real Givens rotation'
      write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
      write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        a = r4_uniform_01 ( seed )
        b = r4_uniform_01 ( seed )

        sa = a
        sb = b

        call srotg ( sa, sb, c, s )

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
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests SSCAL.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      real da
      integer i
      real x(n)

      do i = 1, n
        x(i) = real ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  SSCAL multiplies a vector by a scalar.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      da = 5.0E+00
      call sscal ( n, da, x, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SSCAL ( N, ', da, ', X, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      do i = 1, n
        x(i) = real ( i )
      end do

      da = -2.0E+00
      call sscal ( 3, da, x, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SSCAL ( 3, ', da, ', X, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests SSWAP.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      integer i
      real x(n)
      real y(n)

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  SSWAP swaps two vectors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      call sswap ( n, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SSWAP ( N, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = real ( i )
      end do

      do i = 1, n
        y(i) = real ( 100 * i )
      end do

      call sswap ( 3, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  SSWAP ( 3, X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit single precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
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
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

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
