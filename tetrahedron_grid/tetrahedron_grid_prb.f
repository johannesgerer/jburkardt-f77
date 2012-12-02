      program main

c*****************************************************************************80
c
cc TETRAHEDRON_GRID_TEST tests TETRAHEDRON_GRID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TETRAHEDRON_GRID_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TETRAHEDRON_GRID library.'

      call tetrahedron_grid_test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TETRAHEDRON_GRID_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine tetrahedron_grid_test01 ( )

c*********************************************************************72
c
cc TETRAHEDRON_GRID_TEST01 tests TETRAHEDRON_GRID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 80 ) filename
      integer n
      integer ng
      double precision t(3,4)
      double precision tg(3,500)

      save t

      data t /
     &  0.0D+00, 0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 0.0D+00, 
     &  0.0D+00, 1.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) 
     &  '  TETRAHEDRON_GRID can define a tetrahedral grid'
      write ( *, '(a)' ) 
     &  '  with N+1 points on a side, based on any tetrahedron.'

      n = 10
      write ( *, '(a,i8)' ) '  N = ', n

      call tetrahedron_grid_count ( n, ng )

      call r8mat_print ( 3, 4, t, '  Tetrahedron vertices:' )

      call tetrahedron_grid ( n, t, ng, tg )

      call r83vec_print_part ( ng, tg, 20, 
     &  '  Part of the grid point array:' )

      filename = 'tetrahedron_grid_test01.xyz'

      call r8mat_write ( filename, 3, ng, tg )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Data written to the file "' // trim ( filename ) // '".'

      return
      end
