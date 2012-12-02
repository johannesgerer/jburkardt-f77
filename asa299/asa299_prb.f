      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA299_PRB.
c
c  Discussion:
c
c    ASA299_PRB calls the ASA299 test routines.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA299_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Tests for the ASA299 library.'

      call test01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA299_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 tests SIMPLEX_LATTICE_POINT_NEXT.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      integer j
      logical more
      integer t
      parameter ( t = 4 )
      integer x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  SIMPLEX_LATTICE_POINT_NEXT generates'
      write ( *, '(a)' ) '  generates lattice points in the simplex'
      write ( *, '(a)' ) '    0 <= X'
      write ( *, '(a)' ) '    sum ( X(1:N) ) <= T'
      write ( *, '(a,i8)' ) '  Here N = ', n
      write ( *, '(a,i8)' ) '  and T =  ', t
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     Index        X(1)      X(2)      X(3)      X(4)'
      write ( *, '(a)' ) ' '

      more = .false.

      i = 0

10    continue

        call simplex_lattice_point_next ( n, t, more, x )

        i = i + 1

        write ( *, '(2x,i8,2x,4(2x,i8))' ) i, ( x(j), j = 1, n )

        if ( .not. more )  then
          go to 20
        end if

      go to 10

20    continue
 
      return
      end
