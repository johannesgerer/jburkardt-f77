      program main

c*********************************************************************72
c
cc TOMS384_PRB tests the TOMS384 routine SYMQR.
c
c  Modified:
c
c    09 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS384_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS384 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS384_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SYMQR on a full matrix.
c
c  Modified:
c
c    09 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      real a(n,n)
      logical abscnv
      real angle
      real d(n)
      real e(n)
      real eps
      real exact
      integer fail
      integer i
      integer j
      real k0
      integer na
      real pi
      logical trd
      logical vec

      pi = 3.14159265E+00

      k0 = 0.0E+00
      na = n
      eps = 0.00001E+00
      abscnv = .false.
      vec = .true.
      trd = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
c
c  Set the values of the matrix A.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = min ( i, j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix A:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5f10.4)' ) ( a(i,j), j = 1, n )
      end do
c
c  Set the values of the right hand side vector B.
c
      call symqr ( a, d, e, k0, n, na, eps, abscnv, vec, trd,
     &  fail )

      if ( fail .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  FAIL was returned nonzero.'
        write ( *, '(a,i6)' ) '  FAIL = ', fail
        write ( *, '(a)' ) '  Only eigendata FAIL+1, ..., N can be'
        write ( *, '(a)' ) '  relied on.'
      end if
c
c  Print the results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Eigenvalues:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Computed          Exact'
      write ( *, '(a)' ) ' '

      do i = 1, n
        angle = real ( 2 * i - 1 ) * pi / real ( 2 * n + 1 )
        exact = 0.5E+00 / ( 1.0E+00 - cos ( angle ) )
        write ( *, '(2x,g14.6,2x,g14.6)' ) d(i), exact
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Eigenvector matrix:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,5f10.4)' ) ( a(i,j), j = 1, n )
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
