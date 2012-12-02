      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA159_PRB.
c
c  Discussion:
c
c    ASA159_PRB tests the routines in ASA159.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA159_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA159 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA159_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests RCONT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 5 )

      integer a(m,n)
      integer c(n)
      integer i
      integer ierror
      integer jwork(n)
      logical key
      integer ntest
      parameter ( ntest = 10 )
      integer r(m)
      integer seed

      save c
      save r

      data c / 2, 2, 2, 2, 1 /
      data r / 3, 2, 2, 1, 1 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  RCONT2 constructs a random matrix with'
      write ( *, '(a)' ) '  given row and column sums.'

      call i4vec_print ( m, r, '  The rowsum vector:' )
      call i4vec_print ( n, c, '  The columnsum vector: ' )

      key = .false.

      do i = 1, ntest

        call rcont2 ( m, n, r, c, jwork, key, seed, a, ierror )

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' )
     &    '  RCONT2 returned error flag IERROR = ', ierror
          return
        end if

        call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

      end do

      return
      end
