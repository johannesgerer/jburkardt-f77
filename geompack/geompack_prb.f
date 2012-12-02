      program main

c*********************************************************************72
c
cc MAIN is the main program for GEOMPACK_PRB.
c
c  Discussion:
c
c    GEOMPACK_PRB calls a set of problems for GEOMPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GEOMPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the GEOMPACK library.'

      call test005 ( )
      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GEOMPACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests DIAEDG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 4 )
      integer triangle_num
      parameter ( triangle_num = 2 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision alpha_area
      double precision alpha_ave
      double precision alpha_min_swapped
      double precision alpha_min_unswapped
      integer diaedg
      double precision node_xy(2,node_num)
      integer seed
      logical swap
      integer test
      integer test_num
      parameter ( test_num = 10 )
      integer triangle_node(triangle_order,triangle_num)
      integer value

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) '  DIAEDG determines whether two triangles'
      write ( *, '(a)' ) 
     &  '  with a common edge need to "swap" diagonals.'
      write ( *, '(a)' ) 
     &  '  If swapping is indicated, then ALPHA_MIN should decrease.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Swap   ALPHA_MIN   ALPHA_MIN'
      write ( *, '(a)' ) '         Unswapped   Swapped'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
c
c  Generate a random quadrilateral (1,2,3,4).
c
        call quad_convex_random ( seed, node_xy )
c
c  Does it need swapping?
c
        value = diaedg ( 
     &    node_xy(1,1), node_xy(2,1), 
     &    node_xy(1,2), node_xy(2,2), 
     &    node_xy(1,3), node_xy(2,3), 
     &    node_xy(1,4), node_xy(2,4) )

        if ( value .eq. 1 ) then
          swap = .false.
        else
          swap = .true.
        end if
c
c  Compute ALPHA_MIN unswapped.
c
        triangle_node(1,1) = 1
        triangle_node(2,1) = 2
        triangle_node(3,1) = 3
        triangle_node(1,2) = 1
        triangle_node(2,2) = 3
        triangle_node(3,2) = 4

        call alpha_measure ( node_num, node_xy, triangle_order, 
     &    triangle_num, triangle_node, alpha_min_unswapped, alpha_ave, 
     &    alpha_area )
c
c  Compute ALPHA_MIN swapped.
c
        triangle_node(1,1) = 1
        triangle_node(2,1) = 2
        triangle_node(3,1) = 4
        triangle_node(1,2) = 2
        triangle_node(2,2) = 3
        triangle_node(3,2) = 4

        call alpha_measure ( node_num, node_xy, triangle_order, 
     &    triangle_num, triangle_node, alpha_min_swapped, alpha_ave, 
     &    alpha_area )

        if ( .false. ) then
          call r8mat_transpose_print ( 2, node_num, node_xy, 
     &      '  Quadrilateral' )
        end if

        write ( *, '(2x,3x,l1,2x,f10.6,2x,f10.6)' ) 
     &    swap, alpha_min_unswapped, alpha_min_swapped

      end do

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests POINTS_DELAUNAY_NAIVE_2D.
c
c  Diagram:
c
c   c....3&11....
c   c............
c   c............
c   X..9.........
c   c.....5......
c   c...........6
c   c.4.2...10...
c   c.....8...12.
c   V............
c   c..7.........
c   c......1.....
c   c............
c   c............
c   c----V----X--
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 12 )
      integer dim_num
      parameter ( dim_num = 2 )

      integer maxtri
      parameter ( maxtri = 2 * node_num - 3 )

      integer triangle_node(3,maxtri)
      integer triangle_num
      double precision node_xy(dim_num,node_num )

      save node_xy

      data node_xy /
     &   7.0D+00,  3.0D+00, 
     &   4.0D+00,  7.0D+00, 
     &   5.0D+00, 13.0D+00, 
     &   2.0D+00,  7.0D+00, 
     &   6.0D+00,  9.0D+00, 
     &  12.0D+00, 10.0D+00, 
     &   3.0D+00,  4.0D+00, 
     &   6.0D+00,  6.0D+00, 
     &   3.0D+00, 10.0D+00, 
     &   8.0D+00,  7.0D+00, 
     &   5.0D+00, 13.0D+00, 
     &  10.0D+00,  6.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay'
      write ( *, '(a)' ) '  triangulation of a set of points.'

      call r8mat_transpose_print ( dim_num, node_num, node_xy, 
     &  '  The points:' )

      call points_delaunay_naive_2d ( node_num, node_xy, maxtri, 
     &  triangle_num, triangle_node )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of triangles is TRIANGLE_NUM = ', 
     &  triangle_num

      call i4mat_transpose_print ( 3, triangle_num, triangle_node, 
     &  '  The triangles:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R82VEC_PART_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 12 )
      integer dim_num
      parameter ( dim_num = 2 )

      double precision a(dim_num,node_num)
      integer i
      integer j
      integer l
      integer r
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  R82VEC_PART_QUICK_A reorders a vector'
      write ( *, '(a)' ) '  as part of a quick sort.'
      write ( *, '(a,i12)' ) 
     &  '  Using initial random number seed = ', seed

      call r8mat_uniform_01 ( dim_num, node_num, seed, a )

      do j = 1, node_num
        do i = 1, dim_num
          a(i,j) = 10.0D+00 * a(i,j)
        end do
      end do

      call r8mat_transpose_print ( dim_num, node_num, a, 
     &  '  Before rearrangment:' )

      call r82vec_part_quick_a ( node_num, a, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rearranged array'
      write ( *, '(a,i6)' ) '  Left index =  ', l
      write ( *, '(a,i6)' ) '  Key index =   ', l+1
      write ( *, '(a,i6)' ) '  Right index = ', r

      call r8mat_transpose_print ( dim_num, l,     a(1:dim_num,1:l),   
     &  '  Left half:' )

      call r8mat_transpose_print ( dim_num, 1,     a(1:dim_num,l+1),   
     &  '  Key:' )

      call r8mat_transpose_print ( dim_num, node_num-l-1, 
     &  a(1:dim_num,l+2:node_num), 
     &  '  Right half:' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R82VEC_SORT_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 12 )
      integer dim_num
      parameter ( dim_num = 2 )

      double precision a(dim_num,node_num)
      integer i
      integer j
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  R82VEC_SORT_QUICK_A sorts a vector'
      write ( *, '(a)' ) '  using quick sort.'
      write ( *, '(a,i12)' ) 
     &  '  Using initial random number seed = ', seed

      call r8mat_uniform_01 ( dim_num, node_num, seed, a )

      do j = 1, node_num
        do i = 1, dim_num
          a(i,j) = 10.0D+00 * a(i,j)
        end do
      end do
c
c  Give a few elements the same first component.
c
      a(1,3) = a(1,5)
      a(1,4) = a(1,12)
c
c  Give a few elements the same second component.
c
      a(2,6) = a(2,1)
      a(2,2) = a(2,9)
c
c  Make two entries equal.
c
      a(1,7) = a(1,11)
      a(2,7) = a(2,11)

      call r8mat_transpose_print ( dim_num, node_num, a, 
     &  '  Before rearrangement:' )

      call r82vec_sort_quick_a ( node_num, a )

      call r8mat_transpose_print ( dim_num, node_num, a, 
     &  '  Sorted array:' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R8TRIS2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 9 )
      integer dim_num
      parameter ( dim_num = 2 )

      integer triangle_neighbor(3,2*node_num)
      integer triangle_node(3,2*node_num)
      integer triangle_num
      double precision node_xy(dim_num,node_num)

      save node_xy

      data node_xy /
     &     0.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     0.2D+00, 0.5D+00, 
     &     0.3D+00, 0.6D+00, 
     &     0.4D+00, 0.5D+00, 
     &     0.6D+00, 0.4D+00, 
     &     0.6D+00, 0.5D+00, 
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  R8TRIS2 computes the Delaunay triangulation'
      write ( *, '(a)' ) '  of a pointset in 2D.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, triangle_num, triangle_node, 
     &  triangle_neighbor )

      call triangulation_order3_print ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests TRIANGLE_CIRCUMCENTER_2D;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer test_num
      parameter ( test_num = 4 )

      double precision center(dim_num)
      integer i
      double precision t(dim_num,3)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For a triangle in 2D:'
      write ( *, '(a)' ) 
     &  '  TRIANGLE_CIRCUMCENTER_2D computes the circumcenter.'

      do i = 1, test_num

        if ( i .eq. 1 ) then
          t(1,1) = 0.0D+00
          t(2,1) = 0.0D+00
          t(1,2) = 1.0D+00
          t(2,2) = 0.0D+00
          t(1,3) = 0.0D+00
          t(2,3) = 1.0D+00
        else if ( i .eq. 2 ) then
          t(1,1) = 0.0D+00
          t(2,1) = 0.0D+00
          t(1,2) = 1.0D+00
          t(2,2) = 0.0D+00
          t(1,3) = 0.5D+00
          t(2,3) = 0.86602539D+00
        else if ( i .eq. 3 ) then
          t(1,1) = 0.0D+00
          t(2,1) = 0.0D+00
          t(1,2) = 1.0D+00
          t(2,2) = 0.0D+00
          t(1,3) = 0.5D+00
          t(2,3) = 10.0D+00
        else if ( i .eq. 4 ) then
          t(1,1) = 0.0D+00
          t(2,1) = 0.0D+00
          t(1,2) = 1.0D+00
          t(2,2) = 0.0D+00
          t(1,3) = 10.0D+00
          t(2,3) = 2.0D+00
        end if

        call r8mat_transpose_print ( dim_num, 3, t, 
     &    '  Triangle vertices' )

        call triangle_circumcenter_2d ( t, center )

        call r8vec_print ( dim_num, center, '  Circumcenter :' )

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests TRIANGULATION_ORDER3_PLOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 9 )
      integer triangle_num
      parameter ( triangle_num = 12 )

      character * ( 80 ) file_name
      integer node_show
      double precision node_xy(dim_num,node_num)
      integer triangle_node(3,triangle_num)
      integer triangle_show

      save node_xy
      save triangle_node

      data node_xy /
     &     0.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     0.2D+00, 0.5D+00, 
     &     0.3D+00, 0.6D+00, 
     &     0.4D+00, 0.5D+00, 
     &     0.6D+00, 0.4D+00, 
     &     0.6D+00, 0.5D+00, 
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /
      data triangle_node /
     &     2, 1, 3, 
     &     3, 1, 6,
     &     2, 3, 4, 
     &     4, 3, 5, 
     &     7, 4, 5, 
     &     5, 3, 6, 
     &     7, 5, 6, 
     &     9, 4, 7, 
     &     6, 1, 8, 
     &     7, 6, 8, 
     &     7, 8, 9, 
     &     2, 4, 9 /

      file_name = 'triangulation_plot.eps'
      node_show = 2
      triangle_show = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  TRIANGULATION_ORDER3_PLOT can plot a triangulation.'

      call triangulation_order3_plot ( file_name, node_num, node_xy, 
     &  triangle_num, triangle_node, node_show, triangle_show )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT has created an'
      write ( *, '(a)' ) 
     &  '  Encapsulated PostScript file (EPS) containing'
      write ( *, '(a)' ) '  an image of the triangulation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This file is called "' // trim ( file_name ) //'".'

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests TRIANGULATION_ORDER3_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 9 )
      integer triangle_num
      parameter ( triangle_num = 12 )

      double precision node_xy(2,node_num)
      integer triangle_node(3,triangle_num)
      integer triangle_neighbor(3,triangle_num)

      save node_xy
      save triangle_node
      save triangle_neighbor

      data node_xy /
     &     0.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     0.2D+00, 0.5D+00, 
     &     0.3D+00, 0.6D+00, 
     &     0.4D+00, 0.5D+00, 
     &     0.6D+00, 0.4D+00, 
     &     0.6D+00, 0.5D+00, 
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /
      data triangle_node /
     &     2, 1, 3, 
     &     3, 1, 6, 
     &     2, 3, 4, 
     &     4, 3, 5, 
     &     7, 4, 5, 
     &     5, 3, 6, 
     &     7, 5, 6, 
     &     9, 4, 7, 
     &     6, 1, 8, 
     &     7, 6, 8, 
     &     7, 8, 9, 
     &     2, 4, 9 /
      data triangle_neighbor /
     &     -28,   2,  3, 
     &       1,   9,  6, 
     &       1,   4, 12, 
     &       3,   6,  5, 
     &       8,   4,  7, 
     &       4,   2,  7, 
     &       5,   6, 10, 
     &      12,   5, 11, 
     &       2, -34, 10, 
     &       7,   9, 11, 
     &      10, -38,  8, 
     &       3,   8, -3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) 
     &  '  TRIANGULATION_ORDER3_PRINT prints out a triangulation.'

      call triangulation_order3_print ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

      return
      end
