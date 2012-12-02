      program main

c*********************************************************************72
c
cc TOMS423_PRB tests the TOMS423 routines DECOMP and SOLVE.
c
c  Discussion:
c
c    Solve A*x = b where A is a given matrix, and B a right hand side.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer ndim

      parameter ( n = 3 )
      parameter ( ndim = 4 )

      real a(ndim,ndim)
      real b(n)
      real determ
      integer i
      integer ip(n)
      integer j

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS423_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS423 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DECOMP factors a general matrix;'
      write ( *, '(a)' ) '  SOLVE solves a factored linear system;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix dimension NDIM = ', ndim
      write ( *, '(a,i6)' ) '  The number of equations is N = ', n
c
c  Set the values of the matrix A.
c
      a(1,1) = 1.0E+00
      a(1,2) = 2.0E+00
      a(1,3) = 3.0E+00

      a(2,1) = 4.0E+00
      a(2,2) = 5.0E+00
      a(2,3) = 6.0E+00

      a(3,1) = 7.0E+00
      a(3,2) = 8.0E+00
      a(3,3) = 0.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,3g14.6)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      b(1) = 14.0E+00
      b(2) = 32.0E+00
      b(3) = 23.0E+00

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

      call decomp ( n, ndim, a, ip )
c
c  Get the determinant.
c
      determ = ip(n)
      do i = 1, n
        determ = determ * a(i,i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Matrix determinant = ', determ

      if ( ip(n) .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DECOMP reports the matrix is singular.'
        return
      end if
c
c  If no error occurred, now use DGESL to solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solve the linear system.'

      call solve ( n, ndim, a, b, ip )
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SOLVE returns the solution:'
      write ( *, '(a)' ) '  (Should be (1,2,3))'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6)' ) b(i)
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS423_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
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
