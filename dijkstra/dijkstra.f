      program main

c*********************************************************************72
c
cc MAIN runs an example of Dijkstra's minimum distance algorithm.
c
c  Discussion:
c
c    Given the distance matrix that defines a graph, we seek a list
c    of the minimum distances between node 0 and all other nodes.
c
c    This program sets up a small example problem and solves it.
c
c    The correct minimum distances are:
c
c      0   35   15   45   49   41
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 June 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
      implicit none

      integer nv
      parameter ( nv = 6 )

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer mind(nv)
      integer ohd(nv,nv)

      call timestamp ( )
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  'DIJKSTRA:'
      write ( *, '(a)' )  '  FORTRAN77 version'
      write ( *, '(a)' )  
     &  '  Use Dijkstra''s algorithm to determine the minimum'
      write ( *, '(a)' )  
     &  '  distance from node 1 to each node in a graph,'
      write ( *, '(a)' )  
     &  '  given the distances between each pair of nodes.'
c
c  Initialize the problem data.
c
      call init ( nv, ohd )
c
c  Print the distance matrix.
c
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  '  Distance matrix:'
      write ( *, '(a)' )  ' '
      do i = 1, nv
        write ( *, '(6(2x,i3))' ) ( ohd(i,j), j = 1, nv )
      end do
c
c  Carry out the algorithm.
c
      call dijkstra_distance ( nv, ohd, mind )
c
c  Print the results.
c
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  '  Minimum distances from node 1:'
      write ( *, '(a)' )  ' '
      do i = 1, nv
        write ( *, '(2x,i2,2x,i2)' ) i, mind(i)
      end do
c
c  Terminate.
c
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  'DIJKSTRA:'
      write ( *, '(a)' )  '  Normal end of execution.'

      write ( *, '(a)' )  ' '
      call timestamp ( )

      stop
      end
      subroutine dijkstra_distance ( nv, ohd, mind )

c*********************************************************************72
c
cc DIJKSTRA_DISTANCE uses Dijkstra's minimum distance algorithm.
c
c  Discussion:
c
c    We essentially build a tree.  We start with only node 0 connected
c    to the tree, and this is indicated by setting CONNECTED(0) = TRUE.
c
c    We initialize MIND(I) to the one step distance from node 0 to node I.
c    
c    Now we search among the unconnected nodes for the node MV whose minimum
c    distance is smallest, and connect it to the tree.  For each remaining
c    unconnected node I, we check to see whether the distance from 0 to MV
c    to I is less than that recorded in MIND(I), and if so, we can reduce
c    the distance.
c
c    After NV-1 steps, we have connected all the nodes to 0, and computed
c    the correct minimum distances.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NV, the number of nodes.
c
c    Input, integer OHD(NV,NV), the distance of the direct
c    link between nodes I and J.
c
c    Output, integer MIND(NV), the minimum 
c    distance from node 1 to each node.
c
      implicit none

      integer nv

      logical connected(nv)
      integer i
      integer md
      integer mind(nv)
      integer mv
      integer ohd(nv,nv)
      integer step 
c
c  Start out with only node 1 connected to the tree.
c
      connected(1) = .true.
      do i = 2, nv
        connected(i) = .false.
      end do
c
c  Initialize the minimum distance to the one-step distance.
c
      do i = 1, nv
        mind(i) = ohd(1,i)
      end do
c
c  Attach one more node on each iteration.
c
      do step = 2, nv
c
c  Find the nearest unconnected node.
c
        call find_nearest ( nv, mind, connected, md, mv )

        if ( mv .eq. - 1 ) then
          write ( *, '(a)' )  ' '
          write ( *, '(a)' )  'DIJKSTRA_DISTANCE - Warning!'
          write ( *, '(a)' )  '  Search terminated early.'
          write ( *, '(a)' )  '  Graph might not be connected.'
          return
        end if
c
c  Mark this node as connected.
c
        connected(mv) = .true.
c
c  Having determined the minimum distance to node MV, see if
c  that reduces the minimum distance to other nodes.
c
        call update_mind ( nv, connected, ohd, mv, mind )

      end do

      return
      end
      subroutine find_nearest ( nv, mind, connected, d, v )

c*********************************************************************72
c
cc FIND_NEAREST finds the nearest unconnected node.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NV, the number of nodes.
c
c    Input, integer MIND(NV), the currently computed minimum 
c    distance from node 1 to each node.
c
c    Input, logical CONNECTED(NV), is true for each connected node, whose 
c    minimum distance to node 1 has been determined.
c
c    Output, integer D, the distance from node 1 to the nearest 
c    unconnected node.
c
c    Output, integer V, the index of the nearest unconnected node.
c
      implicit none

      integer nv

      logical connected(nv)
      integer d
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer mind(nv)
      integer v

      d = i4_huge
      v = -1

      do i = 1, nv
        if ( .not. connected(i) .and. mind(i) .lt. d ) then
          d = mind(i)
          v = i
        end if
      end do

      return
      end
      subroutine init ( nv, ohd )

c*********************************************************************72
c
cc INIT initializes the problem data.
c
c  Discussion:
c
c    The graph uses 6 nodes, and has the following diagram and
c    distance matrix:
c
c    N0--15--N2-100--N3           0   40   15  Inf  Inf  Inf
c      \      |     /            40    0   20   10   25    6
c       \     |    /             15   20    0  100  Inf  Inf
c        40  20  10             Inf   10  100    0  Inf  Inf
c          \  |  /              Inf   25  Inf  Inf    0    8
c           \ | /               Inf    6  Inf  Inf    8    0
c            N1
c            / \
c           /   \
c          6    25
c         /       \
c        /         \
c      N5----8-----N4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NV, the number of nodes.
c
c    Output, integer OHD(NV,NV), the distance of the direct
c    link between nodes I and J.
c
      implicit none

      integer nv

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer ohd(nv,nv)

      do j = 1, nv
        do i = 1, nv
          ohd(i,j) = i4_huge
        end do
      end do

      do i = 1, nv
        ohd(i,i) = 0
      end do

      ohd(1,2) = 40
      ohd(1,3) = 15
      ohd(2,3) = 20
      ohd(2,4) = 10
      ohd(2,5) = 25
      ohd(3,4) = 100
      ohd(2,6) = 6
      ohd(5,6) = 8

      ohd(2,1) = 40
      ohd(3,1) = 15
      ohd(3,2) = 20
      ohd(4,2) = 10
      ohd(5,2) = 25
      ohd(4,3) = 100
      ohd(6,2) = 6
      ohd(6,5) = 8

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine update_mind ( nv, connected, ohd, mv, mind )

c*********************************************************************72
c
cc UPDATE_MIND updates the minimum distance vector.
c
c  Discussion:
c
c    We've just determined the minimum distance to node MV.
c
c    For each node I which is not connected yet,
c    check whether the route from node 0 to MV to I is shorter
c    than the currently known minimum distance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NV, the number of nodes.
c
c    Input, logical CONNECTED(NV), is true for each connected node, whose 
c    minimum distance to node 0 has been determined.
c
c    Input, integer OHD(NV,NV), the distance of the direct link 
c    between nodes I and J.
c
c    Input, integer MV, the node whose minimum distance to node 20
c    has just been determined.
c
c    Input/output, integer MIND(NV), the currently computed
c    minimum distances from node 1 to each node.
c
      implicit none

      integer nv

      logical connected(nv)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer mind(nv)
      integer mv
      integer ohd(nv,nv)

      do i = 1, nv
        if ( .not. connected(i) ) then
c
c  If we really use the maximum integer (or something close) to indicate
c  no link, then we'll get burned if we add it to another value
c  Integer arithmetic can 'wrap around', so that 17 + i4_huge becomes
c  a very negative numberc  So first we eliminate the possiblity that
c  the link is infinite.
c
          if ( ohd(mv,i) .lt. i4_huge ) then
            if ( mind(mv) + ohd(mv,i) .lt. mind(i) ) then
              mind(i) = mind(mv) + ohd(mv,i)
            end if
          end if
        end if
      end do

      return
      end
