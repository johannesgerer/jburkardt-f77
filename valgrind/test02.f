      program main

c*********************************************************************72
c
c  Purpose:
c
c    MAIN is the main program for TEST01.
c
c  Discussion:
c
c    TEST02 has some uninitialized data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2011
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A sample code for analysis by VALGRIND.'

      call junk_data ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine junk_data ( )

c*********************************************************************72
c
c  Purpose:
c
c    JUNK_DATA has some uninitialized variables.
c
c  Discussion:
c
c    VALGRIND's MEMCHECK program monitors uninitialized variables, but does
c    not complain unless such a variable is used in a way that means its
c    value affects the program's results, that is, the value is printed,
c    or computed with.  Simply copying the unitialized data to another variable
c    is of no concern.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 May 2011
c
      implicit none

      integer i
      integer x(10)
c
c  X = { 0, 1, 2, 3, 4, ?a, ?b, ?c, ?d, ?e }.
c
      do i = 1, 5
        x(i) = i - 1
      end do
c
c  Copy some values.
c  X = { 0, 1, ?c, 3, 4, ?b, ?b, ?c, ?d, ?e }.
c
      x(3) = x(8)
      x(6) = x(7)
c
c  Modify some uninitialized entries.
c  Memcheck doesn't seem to care about this.
c
      do i = 1, 10
        x(i) = 2 * x(i)
      end do
c
c  Print X.
c
      do i = 1, 10
        write ( *, '(2x,i2,2x,i2)' ) i, x(i)
      end do

      return
      end
