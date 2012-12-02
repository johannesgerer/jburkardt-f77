      program main

c*********************************************************************72
c
c  Purpose:
c
c    MAIN is the main program for TEST01.
c
c  Discussion:
c
c    TEST01 calls F, which has a memory "leak".  This memory leak can be
c    detected by VALGRID.
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

      integer ( kind = 4 ) n
      parameter ( n = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  A sample code for analysis by VALGRIND.'

      call f ( n )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine f ( n )

c*********************************************************************72
c
c  Purpose:
c
c    F computes N+1 entries of the Fibonacci sequence.
c
c  Discussion:
c
c    Unfortunately, F only allocates space for N entries.  Hence, the
c    assignment of a value to the N+1 entry causes a memory leak.
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

      integer n

      integer i
      integer x(n)

      x(1) = 1
      write ( *, '(2x,i2,2x,i2)' ) 1, x(1)

      x(2) = 1
      write ( *, '(2x,i2,2x,i2)' ) 2, x(2)

      do i = 3, n + 1
        x(i) = x(i-1) + x(i-2)
        write ( *, '(2x,i2,2x,i2)' ) i, x(i)
      end do

      return
      end

