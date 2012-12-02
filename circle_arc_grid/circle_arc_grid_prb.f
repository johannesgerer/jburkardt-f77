      program main

c*********************************************************************72
c
cc MAIN is the main program for CIRCLE_ARC_GRID_PRB.
c
c  Discussion:
c
c    CIRCLE_ARC_GRID_PRB tests CIRCLE_ARC_GRID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CIRCLE_ARC_GRID_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test CIRCLE_ARC_GRID.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CIRCLE_ARC_GRID_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of CIRCLE_ARC_GRID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(2)
      double precision c(2)
      character ( len = 80 ) filename
      double precision r
      double precision xy(2,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Compute points on a 90 degree arc.'

      r = 2.0D+00
      c(1) = 5.0D+00
      c(2) = 5.0D+00
      a(1) = 0.0D+00
      a(2) = 90.0D+00
c
c  Echo the input.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Radius =           ', r
      write ( *, '(a,g14.6,2x,g14.6)' ) 
     &  '  Center =           ', c(1), c(2)
      write ( *, '(a,g14.6)' ) '  Angle 1 =          ', a(1)
      write ( *, '(a,g14.6)' ) '  Angle 2 =          ', a(2)
      write ( *, '(a,i8)' ) '  Number of points = ', n
c
c  Compute the data.
c
      call circle_arc_grid ( r, c, a, n, xy )
c
c  Print a little of the data.
c
      call r82vec_print_part ( n, xy, 5, '  A few of the points:' )
c
c  Write the data.
c
      filename = 'arc.txt'
      call r8mat_write ( filename, 2, n, xy )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Data written to "' // trim ( filename ) // '".'

      return
      end
