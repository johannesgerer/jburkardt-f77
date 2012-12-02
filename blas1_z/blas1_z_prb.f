      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS1_Z_PRB.
c
c  Discussion:
c
c    BLAS1_Z_PRB tests the BLAS1 routines.
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
      write ( *, '(a)' ) 'BLAS1_Z_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS1_Z library.'

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
      write ( *, '(a)' ) 'BLAS1_Z_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DZASUM.
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

      integer ma
      integer na
      integer nx

      parameter ( ma = 5 )
      parameter ( na = 4 )
      parameter ( nx = 8 )

      double complex a(ma,na)
      double precision dzasum
      integer i
      integer j
      double complex x(nx)

      a(1,1) = ( -3.0D+00,  4.0D+00 )
      a(2,1) = (  2.0D+00,  0.0D+00 )
      a(3,1) = (  3.0D+00, -4.0D+00 )
      a(4,1) = (  2.0D+00,  0.0D+00 )
      a(5,1) = (  2.0D+00, -1.0D+00 )
      a(1,2) = ( -1.0D+00,  1.0D+00 )
      a(2,2) = (  0.0D+00,  5.0D+00 )
      a(3,2) = ( -4.0D+00, -2.0D+00 )
      a(4,2) = ( -4.0D+00,  1.0D+00 )
      a(5,2) = ( -4.0D+00, -3.0D+00 )
      a(1,3) = (  0.0D+00, -2.0D+00 )
      a(2,3) = (  1.0D+00,  3.0D+00 )
      a(3,3) = ( -3.0D+00,  3.0D+00 )
      a(4,3) = ( -3.0D+00,  3.0D+00 )
      a(5,3) = ( -1.0D+00, -2.0D+00 )
      a(1,4) = ( -1.0D+00,  2.0D+00 )
      a(2,4) = (  2.0D+00, -4.0D+00 )
      a(3,4) = (  0.0D+00, -1.0D+00 )
      a(4,4) = (  0.0D+00, -1.0D+00 )
      a(5,4) = ( -2.0D+00,  4.0D+00 )

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )
      x(6) = ( -1.0D+00,  0.0D+00 )
      x(7) = (  0.0D+00, -3.0D+00 )
      x(8) = (  4.0D+00,  0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  DZASUM adds the absolute values of'
      write ( *, '(a)' ) '  elements of a complex vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, nx
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  DZASUM ( NX,   X, 1    ) = ', dzasum ( nx,   x, 1 )
      write ( *, '(a,g14.6)' )
     &  '  DZASUM ( NX/2, X, 2    ) = ', dzasum ( nx/2, x, 2 )
      write ( *, '(a,g14.6)' )
     &  '  DZASUM ( 2,    X, NX/2 ) = ', dzasum ( 2,    x, nx/2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate with a matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, ma
        write ( *, '(4(2x,f6.1,2x,f6.1))' ) ( a(i,j), j = 1, na )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  DZASUM ( MA, A(1,2), 1 )   = ', dzasum ( ma, a(1,2), 1 )
      write ( *, '(a,g14.6)' )
     &  '  DZASUM ( NA, A(2,1), MA ) = ', dzasum ( na, a(2,1), ma )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DZNRM2.
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

      integer n
      parameter ( n = 5 )

      double precision dznrm2
      integer i
      integer incx
      double precision norm
      double complex x(n)

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  DZNRM2 returns the Euclidean norm'
      write ( *, '(a)' ) '  of a complex vector.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '( 2x, i6, 2x, f6.1, 2x, f6.1 )' ) i, x(i)
      end do

      incx = 1
      norm = dznrm2 ( n, x, incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)'    ) '  The L2 norm of X is ', norm

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests IZAMAX.
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

      parameter ( n = 5 )

      integer i
      integer incx
      integer izamax
      double complex x(n)
      double precision zabs1

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  IZAMAX returns the index of maximum '
      write ( *, '(a)' ) '  magnitude;'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The entries and ZABS1 magnitudes:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2f8.4,2x,f8.4)' ) i, x(i), zabs1 ( x(i) )
      end do

      incx = 1

      i = izamax ( n, x, incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Note that this is a 1-based index.'
      write ( *, '(a)' ) '  Note that the L1 norm is used.'

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests ZABS1.
c
c  Modified:
c
c    27 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double complex c
      double precision c_norm
      integer i
      integer seed
      double complex c8_uniform_01
      double precision zabs1

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  ZABS1 returns the L1 norm '
      write ( *, '(a)' ) '  of a double complex number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Real      Imaginary              '
      write ( *, '(a)' ) '      Part      Part           ZABS1(Z)'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        c = 5.0D+00 * c8_uniform_01 ( seed )

        c_norm = zabs1 ( c )

        write ( *, '(2x,2f10.4,5x,f10.4)' ) c, c_norm

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests ZABS2.
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

      double complex c
      double precision c_norm
      integer i
      integer seed
      double complex c8_uniform_01
      double precision zabs2

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  ZABS2 returns the L2 norm '
      write ( *, '(a)' ) '  of a double complex number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Real      Imaginary              '
      write ( *, '(a)' ) '      Part      Part           ZABS2(Z)'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        c = 5.0D+00 * c8_uniform_01 ( seed )

        c_norm = zabs2 ( c )

        write ( *, '(2x,2f10.4,5x,f10.4)' ) c, c_norm

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests ZAXPY.
c
c  Modified:
c
c    27 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      double complex s
      double complex x(n)
      double complex y(n)

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )

      y(1) = ( -1.0D+00,  0.0D+00 )
      y(2) = (  0.0D+00, -3.0D+00 )
      y(3) = (  4.0D+00,  0.0D+00 )
      y(4) = ( -3.0D+00,  4.0D+00 )
      y(5) = ( -2.0D+00,  0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  ZAXPY adds a multiple of one '
      write ( *, '(a)' ) '  double complex vector to another.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
      end do

      s = ( 0.50D+00, -1.00D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  The scalar multiplier is: ', s

      call zaxpy ( n, s, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A * X + Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2f10.6)' ) i, y(i)
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests ZCOPY.
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

      double complex a(5,5)
      integer i
      integer j
      double complex x(10)
      double complex y(10)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  ZCOPY copies one double complex vector '
      write ( *, '(a)' ) '  into another.'

      do i = 1, 10
        x(i) = cmplx ( 10 * i, i )
      end do

      do i = 1, 10
        y(i) = cmplx ( 20 * i, 2 * i )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = cmplx ( 10 * i,  j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,2g14.6)' ) i, x(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,2g14.6)' ) i, y(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,10f7.1)' ) ( a(i,j), j = 1, 5 )
      end do

      call zcopy ( 5, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZCOPY ( 5, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,2g14.6)' ) i, y(i)
      end do

      do i = 1, 10
        y(i) = cmplx ( 20 * i, 2 * i )
      end do

      call zcopy ( 3, x, 2, y, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZCOPY ( 3, X, 2, Y, 3 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,2g14.6)' ) i, y(i)
      end do

      call zcopy ( 5, x, 1, a, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZCOPY ( 5, X, 1, A, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,10f7.1)' ) ( a(i,j), j = 1, 5 )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = cmplx ( 10 * i,  j )
        end do
      end do

      call zcopy ( 5, x, 2, a, 5 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZCOPY ( 5, X, 2, A, 5 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,10f7.1)' ) ( a(i,j), j = 1, 5 )
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests ZDOTC.
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
      parameter ( n = 5 )

      integer i
      double complex x(n)
      double complex x_norm
      double complex xy_dot
      double complex y(n)
      double complex zdotc

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )

      y(1) = ( -1.0D+00,  0.0D+00 )
      y(2) = (  0.0D+00, -3.0D+00 )
      y(3) = (  4.0D+00,  0.0D+00 )
      y(4) = ( -3.0D+00,  4.0D+00 )
      y(5) = ( -2.0D+00,  0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  ZDOTC computes the conjugated dot product'
      write ( *, '(a)' ) '  of two double complex vectors.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      x_norm = zdotc ( n, x, 1, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The square of the norm of X, computed as'
      write ( *, '(a,f10.4,2x,f10.4)' ) '  ZDOTC(X,X) = ', x_norm

      xy_dot = zdotc ( n, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.4,2x,f10.4)' )
     &  '  The dot product X.Y* is ', xy_dot

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests ZDOTU.
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
      parameter ( n = 5 )

      integer i
      double complex x_norm
      double complex xy_dot
      double complex x(n)
      double complex y(n)
      double complex zdotu

      x(1) = (  2.0D+00, -1.0D+00 )
      x(2) = ( -4.0D+00, -2.0D+00 )
      x(3) = (  3.0D+00,  1.0D+00 )
      x(4) = (  2.0D+00,  2.0D+00 )
      x(5) = ( -1.0D+00, -1.0D+00 )

      y(1) = ( -1.0D+00,  0.0D+00 )
      y(2) = (  0.0D+00, -3.0D+00 )
      y(3) = (  4.0D+00,  0.0D+00 )
      y(4) = ( -3.0D+00,  4.0D+00 )
      y(5) = ( -2.0D+00,  0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  ZDOTU computes the unconjugated dot product'
      write ( *, '(a)' ) '  of two double complex vectors.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      x_norm = zdotu ( n, x, 1, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The unconjugated dot product ( X dot X )'
      write ( *, '(a)' )
     &  '  (which is NOT the square of the norm of X!):'
      write ( *, '(a,f10.4,2x,f10.4)' ) '  ZDOTU(X,X) = ', x_norm

      xy_dot = zdotu ( n, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, y(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.4,2x,f10.4)' )
     &  '  The dot product ( X dot Y ) is ', xy_dot

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests ZDROT.
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
      parameter ( n = 6 )

      double precision c
      integer i
      double precision s
      doublecomplex x(n)
      doublecomplex y(n)

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      do i = 1, n
        y(i) = cmplx ( 20 * i, 2 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  ZDROT carries out a Givens rotation'
      write ( *, '(a)' ) '  on a double complex vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,2f10.1,2x,2f10.1)' ) i, x(i), y(i)
      end do

      c = 0.5D+00
      s = sqrt ( 1.0D+00 - c * c )
      call zdrot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  ZDROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,2f10.1,2x,2f10.1)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests ZDSCAL.
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
      parameter ( n = 6 )

      double precision da
      integer i
      double complex x(n)

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  ZDSCAL multiplies a double precision'
      write ( *, '(a)' ) '  scalar times a double complex vector.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      da = 5.0D+00
      call zdscal ( n, da, x, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZDSCAL ( N, ', da, ', X, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      da = -2.0D+00
      call zdscal ( 3, da, x, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZDSCAL ( 3, ', da, ', X, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests ZMACH.
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

      double precision zmach

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  ZMACH computes several machine-dependent'
      write ( *, '(a)' ) '  double complex arithmetic parameters.'

      write ( *, '(a)' ) ' '
      write ( *, * ) '  ZMACH(1)  = machine epsilon = ', zmach ( 1 )
      write ( *, * ) '  ZMACH(2)  = a tiny value    = ', zmach ( 2 )
      write ( *, * ) '  ZMACH(3)  = a huge value    = ', zmach ( 3 )

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests ZROTG.
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

      double complex a
      double complex b
      double precision c
      double complex r
      double complex s
      double complex sa
      double complex sb
      integer seed
      integer test
      integer test_num
      double complex c8_uniform_01

      test_num = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  ZROTG generates a complex Givens rotation'
      write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
      write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        a = c8_uniform_01 ( seed )
        b = c8_uniform_01 ( seed )

        sa = a
        sb = b

        call zrotg ( sa, sb, c, s )

        r = sa

        write ( *, '(a)' ) ' '
        write ( *, '(a,2g14.6)' ) '  A =  ', a
        write ( *, '(a,2g14.6)' ) '  B =  ', b
        write ( *, '(a, g14.6)' ) '  C =  ', c
        write ( *, '(a,2g14.6)' ) '  S =  ', s
        write ( *, '(a,2g14.6)' ) '  R =  ', r
        write ( *, '(a,2g14.6)' )
     &    '         C *A+S*B = ',          c   * a + s * b
        write ( *, '(a,2g14.6)' )
     &    '  -conjg(S)*A+C*B = ', -conjg ( s ) * a + c * b

      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests ZSCAL.
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
      parameter ( n = 6 )

      double complex da
      integer i
      double complex x(n)

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  ZSCAL multiplies a double complex scalar '
      write ( *, '(a)' ) '  times a vector.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      da = ( 5.0D+00, 0.0D+00 )
      call zscal ( n, da, x, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,2f8.4,a)' ) '  ZSCAL ( N, (', da, '), X, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      da = ( -2.0D+00, 1.0D+00 )
      call zscal ( 3, da, x, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,2f8.4,a)' ) '  ZSCAL ( 3, (', da, '), X, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f6.1,2x,f6.1)' ) i, x(i)
      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests ZSIGN1.
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

      double complex c1
      double complex c2
      double complex c3
      integer i
      integer seed
      double complex c8_uniform_01
      double complex zsign1

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  ZSIGN1 ( C1, C2 ) transfers the sign of'
      write ( *, '(a)' ) '  double complex C2 to the ZABS1 magnitude '
      write ( *, '(a)' ) '  of C1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '           C1                    C2                    C3'
      write ( *, '(a,a)' )
     &  '  --------------------  --------------------  ',
     &  '--------------------'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        c1 = 5.0D+00 * c8_uniform_01 ( seed )
        c2 = 5.0D+00 * c8_uniform_01 ( seed )
        c3 = zsign1 ( c1, c2 )

        write ( *, '(2x,2f10.4,2x,2f10.4,2x,2f10.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests ZSIGN2.
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

      double complex c1
      double complex c2
      double complex c3
      integer i
      integer seed
      double complex c8_uniform_01
      double complex zsign2

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  ZSIGN2 ( C1, C2 ) transfers the sign of'
      write ( *, '(a)' ) '  double complex C2 to the ZABS2 magnitude '
      write ( *, '(a)' ) '  of C1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '           C1                    C2                    C3'
      write ( *, '(a,a)' )
     &  '  --------------------  --------------------  ',
     &  '--------------------'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        c1 = 5.0D+00 * c8_uniform_01 ( seed )
        c2 = 5.0D+00 * c8_uniform_01 ( seed )
        c3 = zsign2 ( c1, c2 )

        write ( *, '(2x,2f10.4,2x,2f10.4,2x,2f10.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests ZSWAP.
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
      parameter ( n = 5 )

      integer i
      double complex x(n)
      double complex y(n)

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      do i = 1, n
        y(i) = cmplx ( 20 * i, 2 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  ZSWAP swaps two double complex vectors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
      end do

      call zswap ( n, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ZSWAP ( N, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = cmplx ( 10 * i, i )
      end do

      do i = 1, n
        y(i) = cmplx ( 20 * i, 2 * i )
      end do

      call zswap ( 3, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  ZSWAP ( 3, X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,2f7.1,2x,2f7.1)' ) i, x(i), y(i)
      end do

      return
      end
      function c8_uniform_01 ( seed )

c*********************************************************************72
c
cc C8_UNIFORM_01 returns a unit double precision complex pseudorandom number.
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
c    Output, double complex C8_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      double precision r
      integer k
      integer seed
      double precision theta
      double complex c8_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( dble ( seed ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

      c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

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
