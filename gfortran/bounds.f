      program main

c*********************************************************************72
c
cc MAIN is the main program for BOUNDS.
c
c  Discussion:
c
c    BOUNDS is a FORTRAN77 program in which an illegal array reference is made.
c
c    The GFORTRAN compiler switch "-fbounds-check" will catch this.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOUNDS:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This program uses an illegal memory reference.'
      write ( *, '(a)' ) 
     &  '  In this case, an array element A(11) is read,'
      write ( *, '(a)' ) 
     &  '  although the array is only dimensioned for size 10.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Compilation with the GFORTRAN switch -fbounds-check'
      write ( *, '(a)' ) '  will generate a run-time warning.'
c
c  Initialize.
c
      do i = 1, 10
        a(i) = i
      end do
c
c  Add neighbor.
c  (And accidentally invoke nonexistent A(11)c)
c
      do i = 1, 10
        a(i) = a(i) + a(i+1)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOUNDS:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
