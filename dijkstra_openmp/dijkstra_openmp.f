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
c    02 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
      implicit none

      include 'omp_lib.h'

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
      write ( *, '(a)' )  'DIJKSTRA_OPENMP:'
      write ( *, '(a)' )  '  FORTRAN77 version'
      write ( *, '(a)' )  
     &  '  Use Dijkstra''s algorithm to determine the minimum'
      write ( *, '(a)' )  
     &  '  distance from node 1 to each node in a graph,'
      write ( *, '(a)' )  
     &  '  given the distances between each pair of nodes.'
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  
     &  '  Although a very small example is considered, we'
      write ( *, '(a)' )  '  demonstrate the use of OpenMP directives'
      write ( *, '(a)' )  '  for parallel execution.'
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
      write ( *, '(a)' )  'DIJKSTRA_OPENMP:'
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
c    02 July 2010
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

      include 'omp_lib.h'

      integer nv

      logical connected(nv)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer md
      integer mind(nv)
      integer mv
      integer my_first
      integer my_id
      integer my_last
      integer my_md
      integer my_mv
      integer my_step
      integer nth
      integer ohd(nv,nv)
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
c  Begin the parallel region.
c
c$omp parallel private ( my_first, my_id, my_last, my_md, my_mv, my_step )
c$omp&         shared ( connected, md, mind, mv, nth, ohd )

      my_id = omp_get_thread_num ( )
      nth = omp_get_num_threads ( ) 
      my_first =   (   my_id       * nv ) / nth + 1
      my_last  =   ( ( my_id + 1 ) * nv ) / nth
c
c  The SINGLE directive means that the block is to be executed by only
c  one thread, and that thread will be whichever one gets here first.
c
c$omp single
      write ( *, '(a)' ) ' '
      write ( *, * ) '  P', my_id, 
     &  ': Parallel region begins with ', nth, ' threads.'
      write ( *, '(a)' ) ' '
c$omp end single
      write ( *, * ) '  P', my_id, ':  First=', my_first, 
     &  '  Last=', my_last
c
c  Attach one more node on each iteration.
c
      do my_step = 2, nv
c
c  Before we compare the results of each thread, set the shared variable 
c  MD to a big value.  Only one thread needs to do this.
c
c$omp single 
        md = i4_huge
        mv = -1
c$omp end single
c
c  Each thread finds the nearest unconnected node in its part of the graph.
c  Some threads might have no unconnected nodes left.
c
        call find_nearest ( my_first, my_last, nv, mind, connected, 
     &    my_md, my_mv )
c
c  In order to determine the minimum of all the MY_MD's, we must insist
c  that only one thread at a time execute this block!
c
c$omp critical
        if ( my_md .lt. md ) then
          md = my_md
          mv = my_mv
        end if
c$omp end critical
c
c  This barrier means that ALL threads have executed the critical
c  block, and therefore MD and MV have the correct value.  Only then
c  can we proceed.
c
c$omp barrier
c
c  If MV is -1, then NO thread found an unconnected node, so we're done early. 
c  OpenMP does not like to BREAK out of a parallel region, so we'll just have 
c  to let the iteration run to the end, while we avoid doing any more updates.
c
c  Otherwise, we connect the nearest node.
c
c$omp single
        if ( mv .ne. - 1 ) then
          connected(mv) = .true.
          write ( *, * ) '  P', my_id, ': Connecting node ', mv
        else
          write ( *, * ) 
     &      '  P', my_id, ': No connecting node on step ', my_step
        end if
c$omp end single
c
c  Again, we don't want any thread to proceed until the value of
c  CONNECTED is updated.
c
c$omp barrier
c
c  Now each thread should update its portion of the MIND vector,
c  by checking to see whether the trip from 0 to MV plus the step
c  from MV to a node is closer than the current record.
c
        if ( mv .ne. -1 ) then
          call update_mind ( my_first, my_last, nv, connected, ohd, 
     &      mv, mind )
        end if
c
c  Before starting the next step of the iteration, we need all threads 
c  to complete the updating, so we set a BARRIER here.
c
c$omp barrier

      end do
c
c  Once all the nodes have been connected, we can exit.
c
c$omp single
      write ( *, * ) ' '
      write ( *, * ) '  P', my_id, ': Exiting parallel region.'
c$omp end single

c$omp end parallel

      return
      end
      subroutine find_nearest ( s, e, nv, mind, connected, d, v )

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
c    02 July 2010
c
c  Author:
c
c    Original C version by Norm Matloff, CS Dept, UC Davis.
c    This FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer S, E, the first and last nodes that 
c    are to be checked.
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
c    unconnected node in the range S to E.
c
c    Output, integer V, the index of the nearest unconnected node
c    in the range S to E.
c
      implicit none

      integer nv

      logical connected(nv)
      integer d
      integer e
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer mind(nv)
      integer s
      integer v

      d = i4_huge
      v = -1

      do i = s, e
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
      subroutine update_mind ( s, e, nv, connected, ohd, mv, mind )

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
c    Input, integer S, E, the first and last nodes that are 
c    to be checked.
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
c    minimum distances from node 1 to each node.  On output, the values for
c    nodes S through E have been updated.
c
      implicit none

      integer nv

      logical connected(nv)
      integer e
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer mind(nv)
      integer mv
      integer ohd(nv,nv)
      integer s

      do i = s, e
        if ( .not. connected(i) ) then
          if ( ohd(mv,i) .lt. i4_huge ) then
            if ( mind(mv) + ohd(mv,i) .lt. mind(i) ) then
              mind(i) = mind(mv) + ohd(mv,i)
            end if
          end if
        end if
      end do

      return
      end
