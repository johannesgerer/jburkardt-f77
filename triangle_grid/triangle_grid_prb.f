      program main

c*********************************************************************72
c
cc MAIN is the main program for TRIANGLE_GRID_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGLE_GRID_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TRIANGLE_GRID library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGLE_GRID_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests TRIANGLE_GRID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer ng
      parameter ( ng = ((n+1)*(n+2))/2 )

      character * ( 255 ) filename
      integer j
      double precision t(2,3)
      double precision tg(2,ng)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  TRIANGLE_GRID can define a triangular grid of points'
      write ( *, '(a)' ) 
     &  '  with N+1 points on a side, based on any triangle.'

      t(1,1) = 0.0D+00
      t(2,1) = 0.0D+00
      t(1,2) = 1.0D+00
      t(2,2) = 0.0D+00
      t(1,3) = 0.5D+00
      t(2,3) = 0.86602540378443860D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Defining triangle:'
      write ( *, '(a)' ) '     J      X      Y'
      write ( *, '(a)' ) ' '
      do j = 1, 3
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) j, t(1,j), t(2,j)
      end do
      call triangle_grid ( n, t, tg )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     J      X      Y'
      write ( *, '(a)' ) ' '
      do j = 1, ng
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) j, tg(1,j), tg(2,j)
      end do

      filename = 'triangle_grid_test01.xy'

      open ( unit = 1, file = filename, status = 'replace' )
      do j = 1, ng
        write ( 1, '(2x,g14.6,2x,g14.6)' ) tg(1,j), tg(2,j)
      end do
      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data written to "' // trim ( filename ) // '".'


      return
      end

