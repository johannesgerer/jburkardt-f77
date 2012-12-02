      program main

c*********************************************************************72
c
cc MAIN is the main program for TRIANGULATION_PRB.
c
c  Discussion:
c
c    TRIANGULATION_PRB tests routines from the TRIANGULATION library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TRIANGULATION library.'

      call test01 ( )
      call test02 ( )
      call test025 ( )
      call test026 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test125 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
      call test213 ( )
      call test215 ( )
      call test217 ( )
      call test219 ( )
      call test22 ( )
      call test23 ( )
      call test24 ( )
      call test25 ( )
      call test26 ( )
      call test265 ( )
      call test27 ( )

      call test31 ( )
      call test32 ( )
      call test33 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ALPHA_MEASURE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 13 )
      integer triangle_num
      parameter ( triangle_num = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision alpha_area
      double precision alpha_ave
      double precision alpha_min
      double precision node_xy(2,node_num)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  ALPHA_MEASURE returns the ALPHA measure of'
      write ( *, '(a)' ) '  quality of a triangulation.'
c
c  Get the triangulation data.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Compute the triangulation quality.
c
      call alpha_measure ( node_num, node_xy, triangle_order,
     &  triangle_num, triangle_node, alpha_min, alpha_ave, alpha_area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f12.6)' ) '  ALPHA_MIN  = ', alpha_min
      write ( *, '(a,f12.6)' ) '  ALPHA_AVE  = ', alpha_ave
      write ( *, '(a,f12.6)' ) '  ALPHA_AREA = ', alpha_area

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests AREA_MEASURE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 13 )
      integer triangle_num
      parameter ( triangle_num = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision area_ave
      double precision area_max
      double precision area_min
      double precision area_ratio
      double precision area_std
      double precision node_xy(2,node_num)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  AREA_MEASURE returns the AREA measure of'
      write ( *, '(a)' ) '  quality of a triangulation.'
c
c  Get the triangulation data.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Compute the triangulation quality.
c
      call area_measure ( node_num, node_xy, triangle_order,
     &  triangle_num, triangle_node, area_min, area_max, area_ratio,
     &  area_ave, area_std )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f12.6)' ) '  AREA_MIN   = ', area_min
      write ( *, '(a,f12.6)' ) '  AREA_MAX   = ', area_max
      write ( *, '(a,f12.6)' ) '  AREA_RATIO = ', area_ratio
      write ( *, '(a,f12.6)' ) '  AREA_AVE   = ', area_ave
      write ( *, '(a,f12.6)' ) '  AREA_STD   = ', area_std

      return
      end
      subroutine test025 ( )

c*****************************************************************************80
c
cc TEST025 tests DELAUNAY_SWAP_TEST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
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
      integer triangle_node(triangle_order,triangle_num)

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025'
      write ( *, '(a)' )
     &  '  DELAUNAY_SWAP_TEST determines whether two triangles'
      write ( *, '(a)' )
     &  '  with a common edge need to "swap" diagonals.'
      write ( *, '(a)' )
     &  '  If swapping is indicated, then ALPHA_MIN should increase.'
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
        call delaunay_swap_test ( node_xy, swap )
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
      subroutine test026 ( )

c*****************************************************************************80
c
cc TEST026 tests DIAEDG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
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
      integer triangle_node(triangle_order,triangle_num)
      integer value

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026'
      write ( *, '(a)' ) '  DIAEDG determines whether two triangles'
      write ( *, '(a)' )
     &  '  with a common edge need to "swap" diagonals.'
      write ( *, '(a)' )
     &  '  If swapping is indicated, then ALPHA_MIN should increase.'
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
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests NODE_MERGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 15 )
      integer test_num
      parameter ( test_num = 4 )

      integer node
      integer node_rep(node_num)
      double precision node_xy(2,node_num)
      integer rep
      integer rep_num
      integer test
      double precision tolerance
      double precision tolerance_test(test_num)

      save node_xy
      save tolerance_test

      data node_xy /
     &     0.0D+00, 0.0D+00,
     &     1.0D+00, 0.0D+00,
     &     3.0D+00, 0.0D+00,
     &     4.0D+00, 0.0D+00,
     &     1.0D+00, 1.0D+00,
     &     4.0D+00, 1.0D+00,
     &     2.0D+00, 2.0D+00,
     &     3.0D+00, 3.0D+00,
     &     2.0D+00, 3.5D+00,
     &     0.5D+00, 4.0D+00,
     &     1.0D+00, 4.0D+00,
     &     1.5D+00, 4.0D+00,
     &     4.0D+00, 4.0D+00,
     &     1.0D+00, 4.5D+00,
     &     1.0D+00, 4.5D+00 /

      data tolerance_test /
     &  0.01D+00, 0.75D+00, 1.2D+00, 1.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  NODE_MERGE identifies groups of nodes'
      write ( *, '(a)' ) '  that can be merged, with a given tolerance.'

      call r8mat_transpose_print ( 2, node_num, node_xy,
     &  '  Node coordinates:' )

      do test = 1, test_num

        tolerance = tolerance_test(test)

        call node_merge ( dim_num, node_num, node_xy, tolerance,
     &    node_rep )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  TOLERANCE = ', tolerance
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Node  Representatives:'
        write ( *, '(a)' ) ' '

        do node = 1, node_num
          write ( *, '(2x,i8,2x,i8)' ) node, node_rep(node)
        end do
c
c  Make a list of the node representatives.
c
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Rep   Coordinates:'
        write ( *, '(a)' ) ' '

        call i4vec_sort_heap_a ( node_num, node_rep )

        rep_num = 0

        do node = 1, node_num

          if ( 2 .le. node ) then
            if ( node_rep(node-1) .eq. node_rep(node) ) then
              go to 10
            end if
          end if

          rep = node_rep(node)
          rep_num = rep_num + 1

          write ( *, '(2x,i8,2x,2g14.6)' )
     &      rep_num, node_xy(1,rep), node_xy(2,rep)

10        continue

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests NS_ADJ_COL_SET, NS_ADJ_COUNT and NS_ADJ_ROW_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer adj_max
      parameter ( adj_max = 708 )
      integer node_num
      parameter ( node_num = 15 )
      integer triangle_num
      parameter ( triangle_num = 4 )
      integer triangle_order
      parameter ( triangle_order = 6 )
      integer variable_num
      parameter ( variable_num = 36 )

      integer adj_num
      integer adj_col(variable_num+1)
      integer adj_row(adj_max)
      character * ( 80 ) file_name
      integer node
      integer node_show
      integer node_p_variable(node_num)
      integer node_u_variable(node_num)
      integer node_v_variable(node_num)
      double precision node_xy(2,node_num)
      integer num
      integer r
      integer rhi
      integer rlo
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_show
      integer variable

      save node_p_variable
      save node_u_variable
      save node_v_variable
      save node_xy
      save triangle_neighbor
      save triangle_node

      data node_p_variable /
     &  3, -1,  8, -1, 13,
     & -1, -1, -1, -1,
     & 24, -1, 29,
     & -1, -1,
     & 36 /
      data node_u_variable /
     &  1,  4,  6,  9, 11,
     & 14, 16, 18, 20,
     & 22, 25, 27,
     & 30, 32,
     & 34 /
      data node_v_variable /
     &  2,  5,  7, 10, 12,
     & 15, 17, 19, 21,
     & 23, 26, 28,
     & 31, 33,
     &  35 /
      data node_xy /
     & 0.0D+00, 0.0D+00,
     & 0.0D+00, 1.0D+00,
     & 0.0D+00, 2.0D+00,
     & 0.0D+00, 3.0D+00,
     & 0.0D+00, 4.0D+00,
     & 1.0D+00, 0.0D+00,
     & 1.0D+00, 1.0D+00,
     & 1.0D+00, 2.0D+00,
     & 1.0D+00, 3.0D+00,
     & 2.0D+00, 0.0D+00,
     & 2.0D+00, 1.0D+00,
     & 2.0D+00, 2.0D+00,
     & 3.0D+00, 0.0D+00,
     & 3.0D+00, 1.0D+00,
     & 4.0D+00, 0.0D+00 /
      data triangle_neighbor /
     &  -1,  2, -1,
     &   3,  1,  4,
     &   2, -1, -1,
     &  -1, -1,  2 /
      data triangle_node /
     &   1, 10,  3,  6,  7,  2,
     &  12,  3, 10,  8,  7, 11,
     &   3, 12,  5,  8,  9,  4,
     &  10, 15, 12, 13, 14, 11 /

      file_name = 'ns_triangulation.eps'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  For an order 3/order 6 Taylor Hood triangulation'
      write ( *, '(a)' ) '  for Navier Stokes velocity and pressure,'
      write ( *, '(a)' ) '  NS_ADJ_COUNT counts variable adjacencies'
      write ( *, '(a)' ) '    and sets up the sparse compressed column'
      write ( *, '(a)' ) '    column pointer array.'
      write ( *, '(a)' )
     &  '  NS_ADJ_COL_SET sets up the sparse compressed column'
      write ( *, '(a)' ) '    COL vector.'
      write ( *, '(a)' )
     &  '  NS_ADJ_ROW_SET sets up the sparse compressed column'
      write ( *, '(a)' ) '    ROW vector.'
c
c  Plot the example.
c
      node_show = 2
      triangle_show = 2

      call triangulation_order6_plot ( file_name, node_num, node_xy,
     &  triangle_num, triangle_node, node_show, triangle_show )
c
c  Get the count of the variable adjacencies.
c  We don't really need to make this call, since the next
c  call does the calculation as part of getting ADJ_COL.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of variables is ', variable_num

      call ns_adj_count ( node_num, triangle_num, variable_num,
     &  triangle_node, triangle_neighbor, node_u_variable,
     &  node_v_variable, node_p_variable, adj_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)'    ) '  As computed by NS_ADJ_COUNT,'
      write ( *, '(a,i8)' )
     &  '  Number of variable adjacency entries is ', adj_num
c
c  Get the count of the variable adjacencies and the COL vector.
c
      call ns_adj_col_set ( node_num, triangle_num, variable_num,
     &  triangle_node, triangle_neighbor, node_u_variable,
     &  node_v_variable, node_p_variable, adj_num, adj_col )

      write ( *, '(a)' ) ' '
      write ( *, '(a)'    ) '  As computed by NS_ADJ_COL_SET,'
      write ( *, '(a,i8)' )
     &  '  Number of variable adjacency entries is ', adj_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Variable adjacency column pointers:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Variable     First      Last    Number'
      write ( *, '(a)' ) ' '

      do variable = 1, variable_num

        num = adj_col(variable+1) - adj_col(variable)

        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' )
     &    variable, adj_col(variable), adj_col(variable+1)-1, num

      end do
c
c  Get the ROW vector.
c
      call ns_adj_row_set ( node_num, triangle_num, variable_num,
     &  triangle_node, triangle_neighbor, node_u_variable,
     &  node_v_variable, node_p_variable, adj_num, adj_col, adj_row )
c
c  This is a huge array.  We only print out the beginning and end.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Variable adjacency row entries:'
      write ( *, '(a)' ) '  (Partial printout only)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Entry     Row       Col'
      write ( *, '(a)' ) ' '

      do variable = 1, variable_num

        rlo = adj_col(variable)
        rhi = adj_col(variable+1)-1

        if ( variable .le. 3 .or. variable_num - 3 .le. variable ) then

          write ( *, '(a)' ) ' '

          do r = rlo, rhi
            write ( *, '(2x,i8,2x,i8,2x,i8)' )
     &        r, adj_row(r), variable
          end do

        end if

        if ( variable .eq. 3 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  (SKIPPING MANY MANY ENTRIES...)'
          write ( *, '(a)' ) ' '
        end if

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests POINTS_DELAUNAY_NAIVE_2D.
c
c  Diagram:
c
c    !....3&11....
c    !............
c    !............
c    X..9.........
c    !.....5......
c    !...........6
c    !.4.2...10...
c    !.....8...12.
c    V............
c    !..7.........
c    !......1.....
c    !............
c    !............
c    !----V----X--
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtri
      parameter ( maxtri = 20 )
      integer node_num
      parameter ( node_num = 12 )
      integer dim_num
      parameter ( dim_num = 2 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer tri(triangle_order,maxtri)
      integer triangle_num
      double precision node_xy(dim_num,node_num)

      save node_xy

      data node_xy /
     &   7.0D+00,  3.0D+00,
     &   4.0D+00,  7.0D+00,
     &   5.0D+00, 13.0D+00,
     &   2.0D+00,  7.0D+00,
     &   6.0D+00,  9.0D+00,
     &  12.0D+00,  8.0D+00,
     &   3.0D+00,  4.0D+00,
     &   6.0D+00,  6.0D+00,
     &   3.0D+00, 10.0D+00,
     &   8.0D+00,  7.0D+00,
     &   5.0D+00, 13.0D+00,
     &  10.0D+00,  6.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' )
     &  '  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay'
      write ( *, '(a)' ) '    triangulation of a set of nodes.'

      call r8mat_transpose_print ( dim_num, node_num, node_xy,
     &  '  The nodes:' )

      call points_delaunay_naive_2d ( node_num, node_xy, maxtri,
     &  triangle_num, tri )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of triangles is ', triangle_num

      call i4mat_transpose_print ( triangle_order, triangle_num, tri,
     &  '  The triangles:' )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests POINTS_HULL_2D.
c
c  Diagram:
c
c    c....3.......
c    c............
c    c..9.........
c    c.....5......
c    c...........6
c    c.4.2...10...
c    c.....8......
c    c.........12.
c    c..7.........
c    c......1.....
c    c............
c    c............
c    c-----------
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 12 )

      integer i
      integer ival(node_num)
      integer nval
      double precision node_xy(dim_num,node_num)

      save node_xy

      data node_xy /
     &     7.0D+00,  3.0D+00,
     &     4.0D+00,  7.0D+00,
     &     5.0D+00, 13.0D+00,
     &     2.0D+00,  7.0D+00,
     &     6.0D+00,  9.0D+00,
     &    12.0D+00,  8.0D+00,
     &     3.0D+00,  4.0D+00,
     &     6.0D+00,  6.0D+00,
     &     3.0D+00, 10.0D+00,
     &     8.0D+00,  7.0D+00,
     &     5.0D+00, 13.0D+00,
     &    10.0D+00,  6.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  POINTS_HULL_2D computes the convex hull'
      write ( *, '(a)' ) '    of a set of nodes.'

      call r8mat_transpose_print ( dim_num, node_num, node_xy,
     &  '  The nodes:' )

      call points_hull_2d ( node_num, node_xy, nval, ival )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The convex hull is formed by connecting:'
      write ( *, '(a)' ) ' '
      do i = 1, nval
        write ( *, '(2x,i4,2x,i4,2x,2g14.6)' )
     &    i, ival(i), node_xy(1:dim_num,ival(i))
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The correct sequence of nodes is:'
      write ( *, '(a)' ) '  4, 9, 3, 6, 12, 1, 7, (4).'
      write ( *, '(a)' ) ' '

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests Q_MEASURE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 13 )
      integer triangle_num
      parameter ( triangle_num = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision node_xy(2,node_num)
      double precision q_area
      double precision q_ave
      double precision q_max
      double precision q_min
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  Q_MEASURE returns the Q measure of'
      write ( *, '(a)' ) '  quality of a triangulation.'
c
c  Get the triangulation data.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Compute the triangulation quality.
c
      call q_measure ( node_num, node_xy, triangle_order,
     &  triangle_num, triangle_node, q_min, q_max, q_ave, q_area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f12.6)' ) '  Q_MIN  = ', q_min
      write ( *, '(a,f12.6)' ) '  Q_MAX  = ', q_max
      write ( *, '(a,f12.6)' ) '  Q_AVE  = ', q_ave
      write ( *, '(a,f12.6)' ) '  Q_AREA = ', q_area

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests R8TRIS2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 June 2009
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
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision node_xy(dim_num,node_num)
      integer triangle_num
      integer triangle(triangle_order,2*node_num)
      integer triangle_neighbor(3,2*node_num)

      save node_xy

      data node_xy /
     &  0.0D+00, 0.0D+00,
     &  0.0D+00, 1.0D+00,
     &  0.2D+00, 0.5D+00,
     &  0.3D+00, 0.6D+00,
     &  0.4D+00, 0.5D+00,
     &  0.6D+00, 0.4D+00,
     &  0.6D+00, 0.5D+00,
     &  1.0D+00, 0.0D+00,
     &  1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  R8TRIS2 computes the Delaunay'
      write ( *, '(a)' ) '    triangulation of a set of nodes in 2D.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, triangle_num, triangle,
     &  triangle_neighbor )

      call triangulation_order3_print ( node_num, triangle_num,
     &  node_xy, triangle, triangle_neighbor )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer j
      double precision phy(2,n)
      double precision r8_uniform_01
      double precision ref(2,n)
      double precision ref2(2,n)
      integer seed
      double precision t(2,3)

      save t

      data t /
     &  1.0D+00, 1.0D+00,
     &  3.0D+00, 1.0D+00,
     &  2.0D+00, 5.0D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  For an order 3 triangle,'
      write ( *, '(a)' ) '  TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE '
      write ( *, '(a)' )
     &  '    maps a physical point to a reference point.'
      write ( *, '(a)' ) '  TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL '
      write ( *, '(a)' )
     &  '    maps a reference point to a physical point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )'
      write ( *, '(a)' ) ' '

      call triangle_reference_sample ( n, seed, ref )

      call triangle_order3_reference_to_physical ( t, n, ref, phy )

      call triangle_order3_physical_to_reference ( t, n, phy, ref2 )

      do j = 1, n

        write ( *, '(2x,2f8.4,2x,2f8.4,2x,2f8.4)' )
     &    ref(1:2,j), phy(1:2,j), ref2(1:2,j)

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision ref2(2,n)
      integer seed
      double precision t(2,6)

      save t

      data t /
     &  7.0D+00, 2.0D+00,
     &  9.0D+00, 2.0D+00,
     &  7.0D+00, 3.0D+00,
     &  8.0D+00, 2.0D+00,
     &  8.0D+00, 2.5D+00,
     &  7.0D+00, 2.5D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  For an order 6 triangle,'
      write ( *, '(a)' ) '  TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE '
      write ( *, '(a)' )
     &  '    maps a physical point to a reference point.'
      write ( *, '(a)' ) '  TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL '
      write ( *, '(a)' )
     &  '    maps a reference point to a physical point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )'
      write ( *, '(a)' ) ' '

      call triangle_reference_sample ( n, seed, ref )

      call triangle_order6_reference_to_physical ( t, n, ref, phy )

      call triangle_order6_physical_to_reference ( t, n, phy, ref2 )

      do j = 1, n

        write ( *, '(2x,2f8.4,2x,2f8.4,2x,2f8.4)' )
     &    ref(1:2,j), phy(1:2,j), ref2(1:2,j)

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests TRIANGULATION_NODE_ORDER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 36 )
      integer triangle_num
      parameter ( triangle_num = 41 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer node_order(node_num)
      integer triangle_node(triangle_order,triangle_num)

      save triangle_node

      data triangle_node /
     &   1,  8,  7,
     &   1,  2,  8,
     &   2,  9,  8,
     &   2,  3,  9,
     &   3, 10,  9,
     &   3,  4, 10,
     &   4, 11, 10,
     &   4,  5, 11,
     &   5, 12, 11,
     &   5,  6, 12,
     &   7, 14, 13,
     &   7,  8, 14,
     &   8, 15, 14,
     &   8,  9, 15,
     &  11, 18, 17,
     &  11, 12, 18,
     &  13, 20, 19,
     &  13, 14, 20,
     &  14, 21, 20,
     &  14, 15, 21,
     &  15, 22, 21,
     &  15, 16, 22,
     &  16, 23, 22,
     &  16, 17, 23,
     &  17, 24, 23,
     &  17, 18, 24,
     &  19, 26, 25,
     &  19, 20, 26,
     &  21, 28, 27,
     &  21, 22, 28,
     &  25, 30, 29,
     &  25, 26, 30,
     &  26, 31, 30,
     &  27, 32, 31,
     &  27, 28, 32,
     &  29, 34, 33,
     &  29, 30, 34,
     &  30, 35, 34,
     &  30, 31, 35,
     &  31, 36, 35,
     &  31, 32, 36 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  TRIANGULATION_NODE_ORDER computes the order'
      write ( *, '(a)' ) '  of the nodes in a triangulation.'

      call triangulation_node_order ( triangle_order, triangle_num,
     &  triangle_node, node_num, node_order )

      call i4vec_print ( node_num, node_order, '  NODE ORDER:' )

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests TRIANGULATION_ORDER3_ADJ_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 Jne 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer adj_max
      parameter ( adj_max = 69 )
      integer node_max
      parameter ( node_max = 13 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj(adj_max)
      integer adj_num
      integer adj_col(node_max+1)
      integer hole_num
      integer k
      integer node
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_ADJ_COUNT counts node adjacencies'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_ADJ_SET sets node adjacencies.'
c
c  Get the sizes of the example.
c
      call triangulation_order3_example1_size ( node_num, triangle_num,
     &  hole_num )
c
c  Get the data of the example.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Get the count of the node adjacencies.
c
      call triangulation_order3_adj_count ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Adjacency pointers:'
      write ( *, '(a)' ) ' '
      do node = 1, node_num
        write ( *, '(2x,i8,2x,i8,2x,i8)' )
     &    node, adj_col(node), adj_col(node+1)-1
      end do
c
c  Get the node adjacencies.
c
      call triangulation_order3_adj_set ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col, adj )
c
c  Print the node adjacencies.
c
      do node = 1, node_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Nodes adjacent to node ', node
        write ( *, '(a)' ) ' '

        do k = adj_col(node), adj_col(node+1)-1
          write ( *, '(2x,i8)' ) adj(k)
        end do

      end do

      return
      end
      subroutine test125 ( )

c*********************************************************************72
c
cc TEST125 tests TRIANGULATION_ORDER3_ADJ_SET2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer adj_max
      parameter ( adj_max = 69 )
      integer node_max
      parameter ( node_max = 13 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj
      integer adj_num
      integer adj_col(node_max+1)
      integer hole_num
      integer ia(adj_max)
      integer ja(adj_max)
      integer k
      integer node
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST125'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_ADJ_COUNT counts node adjacencies'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_ADJ_SET2 sets node adjacencies'
      write ( *, '(a)' ) '  as a pair of vectors IA(*), JA(*).'
c
c  Get the sizes of the example.
c
      call triangulation_order3_example1_size ( node_num, triangle_num,
     &  hole_num )
c
c  Get the data of the example.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Get the count of the node adjacencies.
c
      call triangulation_order3_adj_count ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Adjacency pointers:'
      write ( *, '(a)' ) ' '
      do node = 1, node_num
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) node, adj_col(node),
     &    adj_col(node+1)-1
      end do
c
c  Get the node adjacencies.
c
      call triangulation_order3_adj_set2 ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col, ia, ja )
c
c  Print the node adjacencies.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Node adjacencies stored in IA(*) and JA(*):'
      write ( *, '(a)' ) ' '

      do adj = 1, adj_num

        write ( *, '(2x,i8,2x,a,i2,a,i2,a)' )
     &    adj, '(', ia(adj), ',', ja(adj), ')'

      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 36 )
      integer triangle_num
      parameter ( triangle_num = 41 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer boundary_edge_num
      character * ( 80 ) file_name
      integer node_show
      double precision node_xy(dim_num,node_num)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_show

      save node_xy
      save triangle_node

      data node_xy /
     &  0.0D+00, 0.0D+00,
     &  1.0D+00, 0.0D+00,
     &  2.0D+00, 0.0D+00,
     &  3.0D+00, 0.0D+00,
     &  4.0D+00, 0.0D+00,
     &  5.0D+00, 0.0D+00,
     &  0.0D+00, 1.0D+00,
     &  1.0D+00, 1.0D+00,
     &  2.0D+00, 1.0D+00,
     &  3.0D+00, 1.0D+00,
     &  4.0D+00, 1.0D+00,
     &  5.0D+00, 1.0D+00,
     &  0.0D+00, 2.0D+00,
     &  1.0D+00, 2.0D+00,
     &  2.0D+00, 2.0D+00,
     &  3.0D+00, 2.0D+00,
     &  4.0D+00, 2.0D+00,
     &  5.0D+00, 2.0D+00,
     &  0.0D+00, 3.0D+00,
     &  1.0D+00, 3.0D+00,
     &  2.0D+00, 3.0D+00,
     &  3.0D+00, 3.0D+00,
     &  4.0D+00, 3.0D+00,
     &  5.0D+00, 3.0D+00,
     &  0.0D+00, 4.0D+00,
     &  1.0D+00, 4.0D+00,
     &  2.0D+00, 4.0D+00,
     &  3.0D+00, 4.0D+00,
     &  0.0D+00, 5.0D+00,
     &  1.0D+00, 5.0D+00,
     &  2.0D+00, 5.0D+00,
     &  3.0D+00, 5.0D+00,
     &  0.0D+00, 6.0D+00,
     &  1.0D+00, 6.0D+00,
     &  2.0D+00, 6.0D+00,
     &  3.0D+00, 6.0D+00 /
      data triangle_node /
     &   1,  8,  7,
     &   1,  2,  8,
     &   2,  9,  8,
     &   2,  3,  9,
     &   3, 10,  9,
     &   3,  4, 10,
     &   4, 11, 10,
     &   4,  5, 11,
     &   5, 12, 11,
     &   5,  6, 12,
     &   7, 14, 13,
     &   7,  8, 14,
     &   8, 15, 14,
     &   8,  9, 15,
     &  11, 18, 17,
     &  11, 12, 18,
     &  13, 20, 19,
     &  13, 14, 20,
     &  14, 21, 20,
     &  14, 15, 21,
     &  15, 22, 21,
     &  15, 16, 22,
     &  16, 23, 22,
     &  16, 17, 23,
     &  17, 24, 23,
     &  17, 18, 24,
     &  19, 26, 25,
     &  19, 20, 26,
     &  21, 28, 27,
     &  21, 22, 28,
     &  25, 30, 29,
     &  25, 26, 30,
     &  26, 31, 30,
     &  27, 32, 31,
     &  27, 28, 32,
     &  29, 34, 33,
     &  29, 30, 34,
     &  30, 35, 34,
     &  30, 31, 35,
     &  31, 36, 35,
     &  31, 32, 36 /

      node_show = 2
      triangle_show = 2
      file_name = 'triangulation_order3_plot2.eps'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the'
      write ( *, '(a)' ) '    boundary edges.'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_PLOT plots the triangulation.'

      call triangulation_order3_plot ( file_name, node_num, node_xy,
     &  triangle_num, triangle_node, node_show, triangle_show )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  An Encapsulated PostScript image of this'
      write ( *, '(a)' )
     &  '  triangulation is in "' // trim ( file_name ) // '".'

      call triangulation_order3_boundary_edge_count ( triangle_num,
     &  triangle_node, boundary_edge_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  Number of boundary edges = ', boundary_edge_num
      write ( *, '(a,i8)' ) '  Correct number =           ', 33

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer boundary_num
      integer hole_num
      integer node_num
      integer triangle_num

      hole_num = 2
      node_num = 36
      triangle_num = 41

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER'
      write ( *, '(a)' )
     &  '  determines the number of edges that lie on the'
      write ( *, '(a)' )
     &  '  boundary of a region that has been triangulated.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num
      write ( *, '(a,i8)' )
     &  '  Number of triangles =      ', triangle_num
      write ( *, '(a,i8)' ) '  Number of holes =          ', hole_num

      call triangulation_order3_boundary_edge_count_euler ( node_num,
     &  triangle_num, hole_num, boundary_num )

      write ( *, '(a,i8)' )
     &  '  Number of boundary edges = ', boundary_num

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests TRIANGULATION_ORDER3_BOUNDARY_NODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 36 )
      integer triangle_num
      parameter ( triangle_num = 41 )
      integer triangle_order
      parameter (triangle_order = 3 )

      integer i
      logical node_boundary(node_num)
      integer triangle_node(triangle_order,triangle_num)

      save triangle_node

      data triangle_node /
     &   1,  8,  7,
     &   1,  2,  8,
     &   2,  9,  8,
     &   2,  3,  9,
     &   3, 10,  9,
     &   3,  4, 10,
     &   4, 11, 10,
     &   4,  5, 11,
     &   5, 12, 11,
     &   5,  6, 12,
     &   7, 14, 13,
     &   7,  8, 14,
     &   8, 15, 14,
     &   8,  9, 15,
     &  11, 18, 17,
     &  11, 12, 18,
     &  13, 20, 19,
     &  13, 14, 20,
     &  14, 21, 20,
     &  14, 15, 21,
     &  15, 22, 21,
     &  15, 16, 22,
     &  16, 23, 22,
     &  16, 17, 23,
     &  17, 24, 23,
     &  17, 18, 24,
     &  19, 26, 25,
     &  19, 20, 26,
     &  21, 28, 27,
     &  21, 22, 28,
     &  25, 30, 29,
     &  25, 26, 30,
     &  26, 31, 30,
     &  27, 32, 31,
     &  27, 28, 32,
     &  29, 34, 33,
     &  29, 30, 34,
     &  30, 35, 34,
     &  30, 31, 35,
     &  31, 36, 35,
     &  31, 32, 36 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_BOUNDARY_NODE determines which'
      write ( *, '(a)' ) '  nodes lie on the boundary.'

      call triangulation_order3_boundary_node ( node_num, triangle_num,
     &  triangle_node, node_boundary )

      call lvec_print ( node_num, node_boundary, '    Node  BN?' )

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests TRIANGULATION_ORDER3_CHECK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer triangle_num
      parameter ( triangle_num = 16 )
      integer node_num
      parameter ( node_num = 13 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i
      integer ierror
      integer isave
      integer node_num2
      integer triangle_num2
      integer triangle_node(triangle_order,triangle_num)

      save triangle_node

      data triangle_node /
     &   3,   4,   1,
     &   3,   1,   2,
     &   3,   2,   8,
     &   2,   1,   5,
     &   8,   2,  13,
     &   8,  13,   9,
     &   3,   8,   9,
     &  13,   2,   5,
     &   9,  13,   7,
     &   7,  13,   5,
     &   6,   7,   5,
     &   9,   7,   6,
     &  10,   9,   6,
     &   6,   5,  12,
     &  11,   6,  12,
     &  10,   6,  11 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_CHECK checks the triangulation.'

      call i4mat_transpose_print ( triangle_order, triangle_num,
     &  triangle_node, '  Triangles:' );
c
c  Pass all tests.
c
      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )

      write ( *, '(a,i8)' ) '  Error code = ', ierror
c
c  Fail test 1.
c
      node_num2 = 2

      call triangulation_order3_check ( node_num2, triangle_num,
     &  triangle_node, ierror )

      write ( *, '(a,i8)' ) '  Error code = ', ierror
c
c  Fail test 2.
c
      triangle_num2 = 0

      call triangulation_order3_check ( node_num, triangle_num2,
     &  triangle_node, ierror )

      write ( *, '(a,i8)' ) '  Error code = ', ierror
c
c  Fail test 3.
c
      isave = triangle_node(2,5)
      triangle_node(2,5) = 0

      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )

      write ( *, '(a,i8)' ) '  Error code = ', ierror
      triangle_node(2,5) = isave
c
c  Fail test 4.
c
      isave = triangle_node(3,10)
      triangle_node(3,10) = 2 * node_num + 1

      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )

      write ( *, '(a,i8)' ) '  Error code = ', ierror
      triangle_node(3,10) = isave
c
c  Fail test 5.
c
      triangle_node(3,4) = 3
      triangle_node(3,8) = 3
      triangle_node(3,10) = 3
      triangle_node(3,11) = 3
      triangle_node(2,14) = 3

      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )
      write ( *, '(a,i8)' ) '  Error code = ', ierror

      triangle_node(3,4) = 5
      triangle_node(3,8) = 5
      triangle_node(3,10) = 5
      triangle_node(3,11) = 5
      triangle_node(2,14) = 5
c
c  Fail test 6.
c
      triangle_node(1,9) = 7
      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )
      write ( *, '(a,i8)' ) '  Error code = ', ierror
      triangle_node(1,9) = 9
c
c  Fail test 7.
c
      triangle_node(3,7) = 2
      call triangulation_order3_check ( node_num, triangle_num,
     &  triangle_node, ierror )
      write ( *, '(a,i8)' ) '  Error code = ', ierror
      triangle_node(3,7) = 9

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests TRIANGULATION_ORDER3_EXAMPLE1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_max
      parameter ( node_max = 13 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer hole_num
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_EXAMPLE1_SIZE gives the sizes'
      write ( *, '(a)' ) '    for an example triangulation;'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_EXAMPLE1 returns the information'
      write ( *, '(a)' ) '    for an example triangulation;'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_PRINT prints a triangulation.'
c
c  Get the sizes.
c
      call triangulation_order3_example1_size ( node_num, triangle_num,
     &  hole_num )

      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      call triangulation_order3_print ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests TRIANGULATION_ORDER3_NEIGHBOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer triangle_num
      parameter ( triangle_num = 16 )
      integer node_num
      parameter ( node_num = 13 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer s1
      integer s2
      integer t1
      integer t2
      integer triangle_node(triangle_order,triangle_num)

      save triangle_node

      data triangle_node /
     &   3,   4,   1,
     &   3,   1,   2,
     &   3,   2,   8,
     &   2,   1,   5,
     &   8,   2,  13,
     &   8,  13,   9,
     &   3,   8,   9,
     &  13,   2,   5,
     &   9,  13,   7,
     &   7,  13,   5,
     &   6,   7,   5,
     &   9,   7,   6,
     &  10,   9,   6,
     &   6,   5,  12,
     &  11,   6,  12,
     &  10,   6,  11 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' ) '  TRIANGULATION_ORDER3_NEIGHBOR determines'
      write ( *, '(a)' ) '  the triangle neighbors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T1  S1  T2  S2'
      write ( *, '(a)' ) ' '
      do t1 = 1, triangle_num
        do s1 = 1, 3
          call triangulation_order3_neighbor ( triangle_num,
     &      triangle_node, t1, s1, t2, s2 )
          write ( *, '(2x,4i4)' ) t1, s1, t2, s2
        end do
      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer triangle_num
      parameter ( triangle_num = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      save triangle_node

      data triangle_node /
     &   3,   4,   1,
     &   3,   1,   2,
     &   3,   2,   8,
     &   2,   1,   5,
     &   8,   2,  13,
     &   8,  13,   9,
     &   3,   8,   9,
     &  13,   2,   5,
     &   9,  13,   7,
     &   7,  13,   5,
     &   6,   7,   5,
     &   9,   7,   6,
     &  10,   9,   6,
     &   6,   5,  12,
     &  11,   6,  12,
     &  10,   6,  11 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines the'
      write ( *, '(a)' )
     &  '    adjacency relationships between triangles.'

      call i4mat_transpose_print ( triangle_order, triangle_num,
     &  triangle_node, '  Triangles:' )

      call triangulation_order3_neighbor_triangles ( triangle_num,
     &  triangle_node, triangle_neighbor )

      call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
     &  '  Triangle neighbors:' )

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests TRIANGULATION_ORDER3_PLOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_num
      parameter ( node_num = 13 )
      integer triangle_num
      parameter ( triangle_num = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      character * ( 80 ) file_name
      integer hole_num
      integer node_show
      double precision node_xy(2,node_num)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_show

      file_name = 'triangulation_order3_plot.eps'
      node_show = 2
      triangle_show = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_PLOT can plot a triangulation.'
c
c  Get the example data.
c
      call triangulation_order3_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Make the plot.
c
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
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests TRIANGULATION_ORDER3_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
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
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision node_xy(2,node_num)
      integer triangle_node(triangle_order,triangle_num)
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
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_PRINT prints out the data.'

      call triangulation_order3_print ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      return
      end
      subroutine test213 ( )

c*********************************************************************72
c
cc TEST213 tests TRIANGULATION_ORDER3_QUAD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_max
      parameter ( node_max = 81 )
      integer quad_num
      parameter ( quad_num = 6 )
      integer triangle_max
      parameter ( triangle_max = 128 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i
      integer j
      integer k
      integer n
      integer n11
      integer n12
      integer n21
      integer n22
      integer node_num
      double precision node_xy(2,node_max)
      external quad_fun
      double precision quad_value
      double precision quad_w(quad_num)
      double precision quad_xy(2,quad_num)
      double precision region_area
      integer test
      integer triangle_node(triangle_order,triangle_max)
      integer test_num
      integer triangle_num
      double precision x
      double precision y

      save quad_w
      save quad_xy

      data quad_w /
     &  0.1666666666666666D+00,
     &  0.1666666666666666D+00,
     &  0.1666666666666666D+00,
     &  0.1666666666666666D+00,
     &  0.1666666666666666D+00,
     &  0.1666666666666666D+00 /
      data quad_xy /
     &  0.659027622374092D+00,  0.231933368553031D+00,
     &  0.659027622374092D+00,  0.109039009072877D+00,
     &  0.231933368553031D+00,  0.659027622374092D+00,
     &  0.231933368553031D+00,  0.109039009072877D+00,
     &  0.109039009072877D+00,  0.659027622374092D+00,
     &  0.109039009072877D+00,  0.231933368553031D+00  /

      test_num = 4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST213'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER3_QUAD can apply a quadrature rule'
      write ( *, '(a)' ) '  to every triangle in a triangulated region,'
      write ( *, '(a)' ) '  and estimate the integral of a function'
      write ( *, '(a)' ) '  over that region.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  NODE_NUM   TRI_NUM  Integral estim  Area of Region'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
c
c  Set up the grid.
c
        n = 2**( test - 1 )
        node_num = ( n + 1 ) * ( n + 1 )

        k = 0
        do j = 1, n + 1
          y = dble ( j - 1 ) / dble ( n + 1 - 1 )
          do i = 1, n + 1
            x = dble ( i - 1 ) / dble ( n + 1 - 1 )
            k = k + 1
            node_xy(1,k) = x
            node_xy(2,k) = y
          end do
        end do
c
c  Set up the triangulation.
c
        triangle_num = 2 * n * n

        k = 0
        do j = 1, n
          do i = 1, n

            n11 = i     + ( j     - 1 ) * ( n + 1 )
            n12 = i     + ( j + 1 - 1 ) * ( n + 1 )
            n21 = i + 1 + ( j     - 1 ) * ( n + 1 )
            n22 = i + 1 + ( j + 1 - 1 ) * ( n + 1 )

            k = k + 1
            triangle_node(1,k) = n11
            triangle_node(2,k) = n21
            triangle_node(3,k) = n12
            k = k + 1
            triangle_node(1,k) = n22
            triangle_node(2,k) = n12
            triangle_node(3,k) = n21

          end do
        end do
c
c  Estimate the integral.
c
        call triangulation_order3_quad ( node_num, node_xy,
     &    triangle_order, triangle_num, triangle_node, quad_fun,
     &    quad_num, quad_xy, quad_w, quad_value, region_area )

        write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' )
     &    node_num, triangle_num, quad_value, region_area

      end do

      return
      end
      subroutine quad_fun ( n, xy_vec, f_vec )

c*********************************************************************72
c
cc QUAD_FUN is a sample integrand function for TRIANGULATION_QUAD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision XY_VEC(2,N), the evaluation points.
c
c    Output, double precision F_VEC(N), the value of the integrand
c    function at the evaluation points.
c
      implicit none

      integer n

      double precision f_vec(n)
      integer i
      double precision xy_vec(2,n)

      do i = 1, n
        f_vec(i) = exp ( xy_vec(1,i)**2 + xy_vec(2,i)**2 )
      end do

      return
      end
      subroutine test215 ( )

c*********************************************************************72
c
cc TEST215 tests TRIANGULATION_ORDER3_REFINE_COMPUTE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_max2
      parameter ( node_max2 = 12 )
      integer node_num1
      parameter ( node_num1 = 5 )
      integer triangle_max2
      parameter ( triangle_max2 = 12 )
      integer triangle_num1
      parameter ( triangle_num1 = 3 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer edge_data(5,3*triangle_num1)
      integer node_num2
      double precision node_xy1(dim_num,node_num1)
      double precision node_xy2(dim_num,node_max2)
      integer triangle_node1(triangle_order,triangle_num1)
      integer triangle_node2(triangle_order,triangle_max2)
      integer triangle_num2

      save node_xy1
      save triangle_node1

      data node_xy1 /
     &     0.0D+00, 0.0D+00,
     &     1.0D+00, 0.0D+00,
     &     0.0D+00, 1.0D+00,
     &     1.0D+00, 1.0D+00,
     &     0.5D+00, 1.5D+00  /
      data triangle_node1 /
     &     1, 2, 3,
     &     4, 3, 2,
     &     3, 4, 5 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST215'
      write ( *, '(a)' ) '  For an order3 triangulation:'
      write ( *, '(a)' ) '  TRIANGULATION_ORDER3_REFINE_SIZE determines'
      write ( *, '(a)' ) '  the size of a refined triangulation.'
      write ( *, '(a)' ) '  TRIANGULATION_ORDER3_REFINE_COMPUTES'
      write ( *, '(a)' ) '  computes the refined triangulation.'

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num1
      write ( *, '(a,i8)' )
     &  '  The number of triangles is ', triangle_num1

      call r8mat_transpose_print ( dim_num, node_num1, node_xy1,
     &  '  The nodes' )

      call i4mat_transpose_print ( triangle_order, triangle_num1,
     &  triangle_node1, '  The triangles:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sizing the refined mesh:'

      call triangulation_order3_refine_size ( node_num1, triangle_num1,
     &  triangle_node1, node_num2, triangle_num2, edge_data )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Information about the refined mesh:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num2
      write ( *, '(a,i8)' )
     &  '  The number of triangles is ', triangle_num2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computing the refined mesh:'

      call triangulation_order3_refine_compute ( node_num1,
     &  triangle_num1, node_xy1, triangle_node1, node_num2,
     &  triangle_num2, edge_data, node_xy2, triangle_node2 )

      call r8mat_transpose_print ( dim_num, node_num2, node_xy2,
     &  '  The refined nodes' )

      call i4mat_transpose_print ( triangle_order, triangle_num2,
     &  triangle_node2, '  The refined triangles:' )

      return
      end
      subroutine test217 ( )

c*********************************************************************72
c
cc TEST217 tests TRIANGULATION_ORDER3_SEARCH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 13 )
      integer test_num
      parameter ( test_num = 10 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision alpha
      double precision beta
      double precision d1
      double precision d2
      double precision d3
      double precision dist
      double precision dnear
      integer edge
      double precision gamma
      integer i1
      integer i2
      integer i3
      integer nnear
      double precision node_xy(dim_num,node_num)
      double precision p(dim_num)
      integer seed
      integer step_num
      integer td(test_num)
      integer test
      integer triangle_num
      integer triangle_node(triangle_order,2*node_num)
      integer triangle_index
      integer triangle_neighbor(3,2*node_num)
      double precision xd(dim_num,test_num)

      save node_xy

      data node_xy /
     &     0.0D+00, 0.0D+00,
     &     2.0D+00, 2.0D+00,
     &    -1.0D+00, 3.0D+00,
     &    -2.0D+00, 2.0D+00,
     &     8.0D+00, 2.0D+00,
     &     9.0D+00, 5.0D+00,
     &     7.0D+00, 4.0D+00,
     &     5.0D+00, 6.0D+00,
     &     6.0D+00, 7.0D+00,
     &     8.0D+00, 8.0D+00,
     &    11.0D+00, 7.0D+00,
     &    10.0D+00, 4.0D+00,
     &     6.0D+00, 4.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST217'
      write ( *, '(a)' )
     &  '  Given a set of nodes NODE_XY, and a single point XD,'
      write ( *, '(a)' ) '  find the nearest node in NODE_XY to XD.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  POINTS_POINT_NEAR_NAIVE_ND uses a naive method.'
      write ( *, '(a)' )
     &  '  TRIANGULATION_SEARCH_DELAUNAY finds a triangle'
      write ( *, '(a)' )
     &  '    containing the point.  Often, one of these vertices'
      write ( *, '(a)' ) '    is the closest point.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, triangle_num, triangle_node,
     &  triangle_neighbor )
c
c  Get the test points.
c
      seed = 123456789

      call triangulation_order3_sample ( node_num, node_xy,
     &  triangle_num, triangle_node, test_num, seed, xd, td )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        Y     Distance  Index'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        p(1) = xd(1,test)
        p(2) = xd(2,test)

        call points_point_near_naive_nd ( dim_num, node_num, node_xy,
     &    p, nnear, dnear )

        write ( *, '(a)' ) ' '
        write ( *, '(a,2f8.4   )' ) '  XD       ', p(1), p(2)
        write ( *, '(a,3f8.4,i8)' )
     &    '  Naive    ', node_xy(1,nnear), node_xy(2,nnear), dnear,
     &    nnear

        call triangulation_search_delaunay ( node_num, node_xy,
     &    triangle_order, triangle_num, triangle_node,
     &    triangle_neighbor, p, triangle_index, alpha, beta, 
     &    gamma, edge, step_num )

        if ( triangle_index .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Error: the search failed.'
          cycle
        end if

        i1 = triangle_node(1,triangle_index)
        d1 = sqrt ( ( p(1) - node_xy(1,i1) )**2
     &            + ( p(2) - node_xy(2,i1) )**2 )

        dist = d1
        nnear = i1

        i2 = triangle_node(2,triangle_index)
        d2 = sqrt ( ( p(1) - node_xy(1,i2) )**2
     &            + ( p(2) - node_xy(2,i2) )**2 )

        if ( d2 .lt. dist ) then
          dnear = d2
          nnear = i2
        end if

        i3 = triangle_node(3,triangle_index)
        d3 = sqrt ( ( p(1) - node_xy(1,i3) )**2
     &            + ( p(2) - node_xy(2,i3) )**2 )

        if ( d3 .lt. dist ) then
          dnear = d3
          nnear = i3
        end if

        write ( *, '(a,3f8.4,i8)' ) '  Delaunay ', node_xy(1,nnear),
     &    node_xy(2,nnear), dnear, nnear

      end do

      return
      end
      subroutine test219 ( )

c*********************************************************************72
c
cc TEST219 tests TRIANGULATION_SEARCH_DELAUNAY, TRIANGULATION_SEARCH_NAIVE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      parameter ( node_num = 13 )
      integer test_num
      parameter ( test_num = 10 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision alpha
      double precision beta
      integer edge
      double precision gamma
      double precision node_xy(dim_num,node_num)
      double precision p_test(dim_num,test_num)
      integer seed
      integer step_num
      integer t_test(test_num)
      integer test
      integer triangle_num
      integer triangle_node(triangle_order,2*node_num)
      integer triangle_index1
      integer triangle_index2
      integer triangle_neighbor(3,2*node_num)

      save node_xy

      data node_xy /
     &     0.0D+00, 0.0D+00,
     &     2.0D+00, 2.0D+00,
     &    -1.0D+00, 3.0D+00,
     &    -2.0D+00, 2.0D+00,
     &     8.0D+00, 2.0D+00,
     &     9.0D+00, 5.0D+00,
     &     7.0D+00, 4.0D+00,
     &     5.0D+00, 6.0D+00,
     &     6.0D+00, 7.0D+00,
     &     8.0D+00, 8.0D+00,
     &    11.0D+00, 7.0D+00,
     &    10.0D+00, 4.0D+00,
     &     6.0D+00, 4.0D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST219'
      write ( *, '(a)' ) '  Given a triangulation, and a point P,'
      write ( *, '(a)' ) '  find the triangle T containing to P.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  TRIANGULATION_SEARCH_NAIVE uses a naive method.'
      write ( *, '(a)' )
     &  '  TRIANGULATION_SEARCH_DELAUNAY uses a method that will'
      write ( *, '(a)' )
     &  '    work fast if the triangulation is Delaunay.'
c
c  Set up the Delaunay triangulation.
c
      call r8tris2 ( node_num, node_xy, triangle_num, triangle_node,
     &  triangle_neighbor )
c
c  Get the test points.
c
      call triangulation_order3_sample ( node_num, node_xy,
     &  triangle_num, triangle_node, test_num, seed, p_test, t_test )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        Y     T Naive   T Fast'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        call triangulation_search_naive ( node_num, node_xy,
     &    triangle_order, triangle_num, triangle_node, p_test(1,test),
     &    triangle_index1 )

        call triangulation_search_delaunay ( node_num, node_xy,
     &    triangle_order, triangle_num, triangle_node,
     &    triangle_neighbor, p_test(1,test), triangle_index2,
     &    alpha, beta, gamma, edge, step_num )

        write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,i8)' )
     &    p_test(1,test), p_test(2,test), triangle_index1,
     &    triangle_index2

      end do

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests TRIANGULATION_ORDER6_ADJ_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer adj_max
      parameter ( adj_max = 432 )
      integer node_max
      parameter ( node_max = 48 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer adj(adj_max)
      integer adj_num
      integer adj_col(node_max+1)
      integer hole_num
      integer k
      integer node
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)

      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_ADJ_COUNT counts node adjacencies'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_ADJ_SET sets node adjacencies.'
c
c  Get the sizes of the example.
c
      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )
c
c  Get the data of the example.
c
      call triangulation_order6_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Get the count of the node adjacencies.
c
      call triangulation_order6_adj_count ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of adjacency entries is ', adj_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Adjacency pointers:'
      write ( *, '(a)' ) ' '
      do node = 1, node_num
        write ( *, '(2x,i8,2x,i8,2x,i8)' )
     &    node, adj_col(node), adj_col(node+1)-1
      end do
c
c  Get the node adjacencies.
c
      call triangulation_order6_adj_set ( node_num, triangle_num,
     &  triangle_node, triangle_neighbor, adj_num, adj_col, adj )
c
c  Print the node adjacencies.
c
      do node = 1, node_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Nodes adjacent to node ', node
        write ( *, '(a)' ) ' '

        do k = adj_col(node), adj_col(node+1)-1
          write ( *, '(2x,i8)' ) adj(k)
        end do

      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_max
      parameter ( node_max = 48 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer boundary_edge_num
      integer hole_num
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the'
      write ( *, '(a)' ) '    boundary edges.'

      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )

      call triangulation_order6_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      call triangulation_order6_boundary_edge_count ( triangle_num,
     &  triangle_node, boundary_edge_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  Number of boundary edges = ', boundary_edge_num
      write ( *, '(a,i8)' ) '  Correct number =           ', 16

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer boundary_num
      integer hole_num
      integer node_num
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER'
      write ( *, '(a)' )
     &  '  determines the number of edges that lie on the'
      write ( *, '(a)' )
     &  '  boundary of a region that has been triangulated.'

      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num
      write ( *, '(a,i8)' )
     &  '  Number of triangles =      ', triangle_num
      write ( *, '(a,i8)' ) '  Number of holes =          ', hole_num

      call triangulation_order6_boundary_edge_count_euler ( node_num,
     &  triangle_num, hole_num, boundary_num )

      write ( *, '(a,i8)' )
     &  '  Number of boundary edges = ', boundary_num
      write ( *, '(a,i8)' ) '  Correct number =           ', 16

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests TRIANGULATION_ORDER6_BOUNDARY_NODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_max
      parameter ( node_max = 48 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      character * ( 80 ) file_name
      integer i
      integer hole_num
      logical node_boundary(node_max)
      integer node_num
      integer node_show
      double precision node_xy(2,node_max)
      integer triangle_num
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_show

      file_name = 'triangulation_order6_plot.eps'
      node_show = 2
      triangle_show = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_BOUNDARY_COUNT counts the boundary'
      write ( *, '(a)' ) '    edges.'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_PLOT plots the triangulation.'

      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )

      call triangulation_order6_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )
c
c  Make the plot.
c
      call triangulation_order6_plot ( file_name, node_num, node_xy,
     &  triangle_num, triangle_node, node_show, triangle_show )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  An Encapsulated PostScript image of this'
      write ( *, '(a)' )
     &  '  triangulation is in "' // trim ( file_name ) // '".'

      call triangulation_order6_boundary_node ( node_num, triangle_num,
     &  triangle_node, node_boundary )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Node  BN?'
      write ( *, '(a)' ) ' '

      do i = 1, node_num
        write ( *, '(2x,i8,2x,l1)' ) i, node_boundary(i)
      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests TRIANGULATION_ORDER6_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer node_max
      parameter ( node_max = 48 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer hole_num
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' ) '  TRIANGULATION_ORDER6_PRINT prints the data.'

      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )

      call triangulation_order6_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      call triangulation_order6_print ( node_num, triangle_num, node_xy,
     &  triangle_node, triangle_neighbor )

      return
      end
      subroutine test265 ( )

c*********************************************************************72
c
cc TEST265 tests TRIANGULATION_ORDER6_REFINE_COMPUTE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_max2
      parameter ( node_max2 = 35 )
      integer node_num1
      parameter ( node_num1 = 12 )
      integer triangle_num1
      parameter ( triangle_num1 = 3 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer triangle_max2
      parameter ( triangle_max2 = triangle_num1 * 4 )

      integer edge_data(5,3*triangle_num1)
      integer node_num2
      double precision node_xy1(dim_num,node_num1)
      double precision node_xy2(dim_num,node_max2)
      integer triangle_node1(triangle_order,triangle_num1)
      integer triangle_node2(triangle_order,triangle_max2)
      integer triangle_num2

      save node_xy1
      save triangle_node1

      data node_xy1 /
     &     0.0D+00, 0.0D+00,
     &     2.0D+00, 0.0D+00,
     &     0.0D+00, 2.0D+00,
     &     2.0D+00, 2.0D+00,
     &     1.0D+00, 3.0D+00,
     &     1.0D+00, 0.0D+00,
     &     0.0D+00, 1.0D+00,
     &     1.0D+00, 1.0D+00,
     &     2.0D+00, 1.0D+00,
     &     1.0D+00, 2.0D+00,
     &     0.5D+00, 2.5D+00,
     &     1.5D+00, 2.5D+00  /
      data triangle_node1 /
     &     1,  2,  3,  6,  8,  7,
     &     4,  3,  2,  9, 10,  8,
     &     3,  4,  5, 10, 12, 11 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST265'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_REFINE_SIZE determines the'
      write ( *, '(a)' ) '  size of a refined triangulation.'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_REFINE_COMPUTES computes the'
      write ( *, '(a)' ) '  refined triangulation.'

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num1
      write ( *, '(a,i8)' )
     &  '  The number of triangles is ', triangle_num1

      call r8mat_transpose_print ( dim_num, node_num1, node_xy1,
     &  '  The nodes' )

      call i4mat_transpose_print ( triangle_order, triangle_num1,
     &  triangle_node1, '  The triangles:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sizing the refined mesh:'

      call triangulation_order6_refine_size ( node_num1, triangle_num1,
     &  triangle_node1, node_num2, triangle_num2, edge_data )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Information about the refined mesh:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num2
      write ( *, '(a,i8)' )
     &  '  The number of triangles is ', triangle_num2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computing the refined mesh:'

      call triangulation_order6_refine_compute ( node_num1,
     &  triangle_num1, node_xy1, triangle_node1, node_num2,
     &  triangle_num2, edge_data, node_xy2, triangle_node2 )

      call r8mat_transpose_print ( dim_num, node_num2, node_xy2,
     &  '  The refined nodes' )

      call i4mat_transpose_print ( triangle_order, triangle_num2,
     &  triangle_node2, '  The refined triangles:' )

      return
      end
      subroutine test27 ( )

!*********************************************************************72
!
!! TEST27 tests TRIANGULATION_ORDER6_VERTEX_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2009
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_max
      parameter ( node_max = 48 )
      integer triangle_max
      parameter ( triangle_max = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer hole_num
      integer midside_num
      integer node_num
      double precision node_xy(2,node_max)
      integer triangle_neighbor(3,triangle_max)
      integer triangle_node(triangle_order,triangle_max)
      integer triangle_num
      integer vertex_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  For an order6 triangulation:'
      write ( *, '(a)' )
     &  '  TRIANGULATION_ORDER6_VERTEX_COUNT counts the '
      write ( *, '(a)' ) '  vertex nodes and midside nodes.'

      call triangulation_order6_example1_size ( node_num, triangle_num,
     &  hole_num )

      call triangulation_order6_example1 ( node_num, triangle_num,
     &  node_xy, triangle_node, triangle_neighbor )

      call triangulation_order6_vertex_count ( node_num, triangle_num,
     &  triangle_node, vertex_num, midside_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes =         ', node_num
      write ( *, '(a,i8)' ) '  Number of vertex nodes =  ', vertex_num
      write ( *, '(a,i8)' ) '  Number of midside nodes = ', midside_num

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests VORONOI_POLYGON_AREA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      parameter ( neighbor_num = 4 )
      integer node_num
      parameter ( node_num = 5 )

      double precision area
      double precision area_correct
      integer center
      integer neighbor_index(neighbor_num)
      double precision node_xy(dim_num,node_num)

      save neighbor_index
      save node_xy

      data neighbor_index / 1, 2, 3, 4 /

      data node_xy /
     &  0.0D+00, 0.0D+00,
     &  1.0D+00, 0.0D+00,
     &  1.0D+00, 1.0D+00,
     &  0.0D+00, 1.0D+00,
     &  0.5D+00, 0.5D+00 /

      area_correct = 0.5D+00
      center = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  VORONOI_POLYGON_AREA computes the area of'
      write ( *, '(a)' ) '  a finite Voronoi polygon.'

      call voronoi_polygon_area ( center, neighbor_num, neighbor_index,
     &  node_num, node_xy, area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The computed area is ', area
      write ( *, '(a,g14.6)' ) '  The correct area is  ', area_correct

      return
      end
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests VORONOI_POLYGON_CENTROID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      parameter ( neighbor_num = 4 )
      integer node_num
      parameter ( node_num = 5 )

      integer center
      double precision centroid(dim_num)
      double precision centroid_exact(dim_num)
      integer neighbor_index(neighbor_num)
      double precision node_xy(dim_num,node_num)

      save centroid_exact
      save neighbor_index
      save node_xy

      data centroid_exact / 0.5D+00, 0.5D+00 /

      data neighbor_index / 1, 2, 3, 4 /

      data node_xy /
     &  0.0D+00, 0.0D+00,
     &  1.0D+00, 0.0D+00,
     &  1.0D+00, 1.0D+00,
     &  0.0D+00, 1.0D+00,
     &  0.5D+00, 0.5D+00 /

      center = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' )
     &  '  VORONOI_POLYGON_CENTROID computes the centroid of'
      write ( *, '(a)' ) '  a finite Voronoi polygon.'

      call voronoi_polygon_centroid ( center, neighbor_num,
     &  neighbor_index, node_num, node_xy, centroid )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' )
     &  '  Computed centroid ', centroid(1), centroid(2)
      write ( *, '(a,2g14.6)' )
     &  '  Correct centroid  ', centroid_exact(1), centroid_exact(2)

      return
      end
      subroutine test33 ( )

c*********************************************************************72
c
cc TEST33 tests VORONOI_POLYGON_VERTICES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      parameter ( neighbor_num = 4 )
      integer node_num
      parameter ( node_num = 5 )

      integer center
      integer i
      integer neighbor_index(neighbor_num)
      double precision node_xy(dim_num,node_num)
      double precision v(2,neighbor_num)

      save neighbor_index
      save node_xy

      data neighbor_index / 1, 2, 3, 4 /

      data node_xy /
     &  0.0D+00, 0.0D+00,
     &  1.0D+00, 0.0D+00,
     &  1.0D+00, 1.0D+00,
     &  0.0D+00, 1.0D+00,
     &  0.5D+00, 0.5D+00 /

      center = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST33'
      write ( *, '(a)' )
     &  '  VORONOI_POLYGON_VERTICES computes the vertices of'
      write ( *, '(a)' ) '  a finite Voronoi polygon.'

      call voronoi_polygon_vertices ( center, neighbor_num,
     &  neighbor_index, node_num, node_xy, v )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Voronoi Polygon Vertex coordinates:'
      write ( *, '(a)' ) ' '
      do i = 1, neighbor_num
        write ( *, '(2x,i8,4x,2g14.6)' ) i, v(1,i), v(2,i)
      end do

      return
      end
