      function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

c*****************************************************************************80
c
cc DIAEDG chooses a diagonal edge.
c
c  Discussion:
c
c    The routine determines whether 0--2 or 1--3 is the diagonal edge
c    that should be chosen, based on the circumcircle criterion, where
c    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
c    quadrilateral in counterclockwise order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, double precision X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
c    coordinates of the vertices of a quadrilateral, given in
c    counter clockwise order.
c
c    Output, integer DIAEDG, chooses a diagonal:
c    +1, if diagonal edge 02 is chosen;
c    -1, if diagonal edge 13 is chosen;
c     0, if the four vertices are cocircular.
c
      implicit none

      double precision ca
      double precision cb
      integer diaedg
      double precision dx10
      double precision dx12
      double precision dx30
      double precision dx32
      double precision dy10
      double precision dy12
      double precision dy30
      double precision dy32
      double precision r8_epsilon
      double precision s
      double precision tol
      double precision tola
      double precision tolb
      double precision x0
      double precision x1
      double precision x2
      double precision x3
      double precision y0
      double precision y1
      double precision y2
      double precision y3

      tol = 100.0D+00 * r8_epsilon ( )

      dx10 = x1 - x0
      dy10 = y1 - y0
      dx12 = x1 - x2
      dy12 = y1 - y2
      dx30 = x3 - x0
      dy30 = y3 - y0
      dx32 = x3 - x2
      dy32 = y3 - y2

      tola = tol * max ( abs ( dx10 ), abs ( dy10 ), 
     &  abs ( dx30 ), abs ( dy30 ) )

      tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), 
     &  abs ( dx32 ), abs ( dy32 ) )

      ca = dx10 * dx30 + dy10 * dy30
      cb = dx12 * dx32 + dy12 * dy32

      if ( tola .lt. ca .and. tolb .lt. cb ) then

        diaedg = - 1

      else if ( ca .lt. - tola .and. cb .lt. - tolb ) then

        diaedg = 1

      else

        tola = max ( tola, tolb )
        s = ( dx10 * dy30 - dx30 * dy10 ) * cb 
     &    + ( dx32 * dy12 - dx12 * dy32 ) * ca

        if ( tola .lt. s ) then
          diaedg = - 1
        else if ( s .lt. - tola ) then
          diaedg = 1
        else
          diaedg = 0
        end if

      end if

      return
      end
      function i4_sign ( x )

c*********************************************************************72
c
cc I4_SIGN evaluates the sign of an I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the number whose sign is desired.
c
c    Output, integer I4_SIGN, the sign of the number.
c
      implicit none

      integer i4_sign
      integer x

      if ( x .lt. 0 ) then
        i4_sign = -1
      else
        i4_sign = +1
      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      subroutine i4mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    39 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      character * ( * ) title

      call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi,
     &   jhi, title )

c*********************************************************************72
c
cc I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
      character*8 ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer  j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      if ( 0 .lt. len ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )  title
      end if

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(i8)' ) i
        end do

        write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
        write ( *, '(a)' ) '  Col'
        write ( *, '(a)' ) ' '

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

          do i2 = 1, inc

            i = i2lo - 1 + i2

            write ( ctemp(i2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      subroutine i4vec_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_D reorders an I4VEC into an descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c    A descending heap is an array A with the property that, for every index J,
c    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
c    2*J and 2*J+1 are legal).
c
c                  A(1)
c                /      \
c            A(2)         A(3)
c          /     \        /  \
c      A(4)       A(5)  A(6) A(7)
c      /  \       /   \
c    A(8) A(9) A(10) A(11)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms,
c    Academic Press, 1978, second edition,
c    ISBN 0-12-519260-6.
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input/output, integer A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ifree
      integer key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n/2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = a(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the larger of the two values,
c  and update M if necessary.
c
            if ( a(m) .lt. a(m+1) ) then
              m = m + 1
            end if

          end if
c
c  If the large descendant is larger than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( a(m) .le. key ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot IFREE.
c
20      continue

        a(ifree) = key

      end do

      return
      end
      subroutine i4vec_min ( n, a, amin )

c*********************************************************************72
c
cc I4VEC_MIN computes the minimum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amin
      integer i

      amin = a(1)

      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      subroutine i4vec_sort_heap_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms,
c    Academic Press, 1978, second edition,
c    ISBN 0-12-519260-6.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if
c
c  1: Put A into descending heap form.
c
      call i4vec_heap_d ( n, a )
c
c  2: Sort A.
c
c  The largest object in the heap is in A(1).
c  Move it to position A(N).
c
      call i4_swap ( a(1), a(n) )
c
c  Consider the diminished heap of size N1.
c
      do n1 = n - 1, 2, -1
c
c  Restore the heap structure of A(1) through A(N1).
c
        call i4vec_heap_d ( n1, a )
c
c  Take the largest object from A(1) and move it to A(N1).
c
        call i4_swap ( a(1), a(n1) )

      end do

      return
      end
      subroutine i4vec_sorted_unique ( n, a, unique_num )

c*********************************************************************72
c
cc I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements in A.
c
c    Input/output, integer A(N).  On input, the sorted
c    integer array.  On output, the unique elements in A.
c
c    Output, integer UNIQUE_NUM, the number of unique elements in A.
c
      implicit none

      integer n

      integer a(n)
      integer itest
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1

      do itest = 2, n

        if ( a(itest) .ne. a(unique_num) ) then
          unique_num = unique_num + 1
          a(unique_num) = a(itest)
        end if

      end do

      return
      end
      function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

c*********************************************************************72
c
cc LRLINE determines if a point is left of, right or, or on a directed line.
c
c  Discussion:
c
c    The directed line is parallel to, and at a signed distance DV from
c    a directed base line from (XV1,YV1) to (XV2,YV2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, double precision XU, YU, the coordinates of the point whose
c    position relative to the directed line is to be determined.
c
c    Input, double precision XV1, YV1, XV2, YV2, the coordinates of two points
c    that determine the directed base line.
c
c    Input, double precision DV, the signed distance of the directed line
c    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
c    DV is positive for a line to the left of the base line.
c
c    Output, integer ( kind = 4 ) LRLINE, the result:
c    +1, the point is to the right of the directed line;
c     0, the point is on the directed line;
c    -1, the point is to the left of the directed line.
c
      implicit none

      double precision dv
      double precision dx
      double precision dxu
      double precision dy
      double precision dyu
      integer ( kind = 4 ) lrline
      double precision r8_epsilon
      double precision t
      double precision tol
      double precision tolabs
      double precision xu
      double precision xv1
      double precision xv2
      double precision yu
      double precision yv1
      double precision yv2

      tol = 100.0D+00 * r8_epsilon ( tol )

      dx = xv2 - xv1
      dy = yv2 - yv1
      dxu = xu - xv1
      dyu = yu - yv1

      tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), 
     &  abs ( dyu ), abs ( dv ) )

      t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

      if ( tolabs .lt. t ) then
        lrline = 1
      else if ( - tolabs .le. t ) then
        lrline = 0
      else
        lrline = -1
      end if

      return
      end
      subroutine perm_check2 ( n, p, base, ierror )

c*********************************************************************72
c
cc PERM_CHECK2 checks that a vector represents a permutation.
c
c  Discussion:
c
c    The routine verifies that each of the integers from BASE to
c    to BASE+N-1 occurs among the N entries of the permutation.
c
c    Set the input quantity BASE to 0, if P is a 0-based permutation,
c    or to 1 if P is a 1-based permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the array to check.
c
c    Input, integer BASE, the index base.
c
c    Output, integer IERROR, error flag.
c    0, the array represents a permutation.
c    nonzero, the array does not represent a permutation.  The smallest
c    missing value is equal to IERROR.
c
      implicit none

      integer n

      integer base
      integer find
      integer ierror
      integer p(n)
      integer seek

      ierror = 0

      do seek = base, base + n - 1

        ierror = 1

        do find = 1, n
          if ( p(find) .eq. seek ) then
            ierror = 0
            go to 10
          end if
        end do

10      continue

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK2 - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a proper permutation.'
          stop
        end if

      end do

      return
      end
      subroutine perm_inverse ( n, p )

c*********************************************************************72
c
cc PERM_INVERSE inverts a permutation "in place".
c
c  Discussion:
c
c    This algorithm assumes that the entries in the permutation vector are
c    strictly positive.  In particular, the value 0 must not occur.
c
c    When necessary, this function shifts the data temporarily so that
c    this requirement is satisfied.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On output, P describes the inverse permutation
c
      implicit none

      integer n

      integer base
      integer i
      integer i0
      integer i1
      integer i2
      integer ierror
      integer is
      integer p(n)
      integer p_min

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        stop
      end if
c
c  Find the least value, and shift data so it begins at 1.
c
      call i4vec_min ( n, p, p_min )
      base = 1
      do i = 1, n
        p(i) = p(i) - p_min + base
      end do
c
c  Check the permutation.
c
      call perm_check2 ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if
c
c  Invert the permutation.
c
      is = 1

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        is = -sign ( 1, p(i) )
        p(i) = sign ( p(i), is )

      end do

      do i = 1, n

        i1 = -p(i)

        if ( 0 .le. i1 ) then

          i0 = i

20        continue

            i2 = p(i1)
            p(i1) = i0

            if ( i2 .lt. 0 ) then
              go to 30
            end if

            i0 = i1
            i1 = i2

          go to 20

30        continue

        end if

      end do
c
c  Undo the shift.
c
      do i = 1, n
        p(i) = p(i) + p_min - base
      end do

      return
      end
      subroutine pwl_interp_2d_scattered_value ( nd, xyd, zd, t_num, 
     &  t, t_neighbor, ni, xyi, zi )

c*********************************************************************72
c
cc PWL_INTERP_2D_SCATTERED_VALUE evaluates a 2d interpolant of scattered data
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XYD(2,ND), the data point coordinates.
c
c    Input, double precision ZD(ND), the data values.
c
c    Input, integer T_NUM, the number of triangles.
c
c    Input, integer T(3,T_NUM), the triangle information.
c
c    Input, integer T_NEIGHBOR(3,T_NUM), the triangle neighbors.
c
c    Input, integer  NI, the number of interpolation points.
c
c    Input, double precision XYI(2,NI), the interpolation point coordinates.
c
c    Output, double precision ZI(NI), the interpolated values.
c
      implicit none

      integer nd
      integer ni
      integer t_num

      double precision alpha
      double precision beta
      double precision gamma
      integer edge
      integer i
      integer j
      integer step_num
      integer t(3,t_num)
      integer t_neighbor(3,t_num)
      double precision xyd(2,nd)
      double precision xyi(2,ni)
      double precision zd(nd)
      double precision zi(ni)

      do i = 1, ni

        call triangulation_search_delaunay ( nd, xyd, 3, t_num, t, 
     &    t_neighbor, xyi(1:2,i), j, alpha, beta, gamma, edge, 
     &    step_num )

        if ( j .eq. -1 ) then
          zi(i) = -1.0D+00
        end if

        zi(i) = alpha * zd(t(1,j)) 
     &        + beta  * zd(t(2,j)) 
     &        + gamma * zd(t(3,j))

      end do

      return
      end
      subroutine r8tris2 ( node_num, node_xy, triangle_num, 
     &  triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
c
c  Discussion:
c
c    The routine constructs the Delaunay triangulation of a set of 2D vertices
c    using an incremental approach and diagonal edge swaps.  Vertices are
c    first sorted in lexicographically increasing (X,Y) order, and
c    then are inserted one at a time from outside the convex hull.
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
c    Original FORTRAN77 version by Barry Joe.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input/output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
c    of the nodes.  On output, the vertices have been sorted into 
c    dictionary order.
c
c    Output, integer TRIANGLE_NUM, the number of triangles in the 
c    triangulation;  TRIANGLE_NUM is equal to 2*NODE_NUM - NB - 2, where NB is 
c    the number of boundary vertices.
c
c    Output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes that 
c    make up each triangle.  The elements are indices of P.  The vertices of
c    the triangles are in counter clockwise order.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbor list.  Positive elements are indices of TIL; negative 
c    elements are used for links of a counter clockwise linked list of boundary
c    edges;  LINK = -(3*I + J-1) where I, J = triangle, edge index;
c    TRIANGLE_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J 
c    to J+1 (mod 3).
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num

      double precision cmax
      integer dim
      integer e
      integer i
      integer ierr
      integer indx(node_num)
      integer j
      integer k
      integer l
      integer ledg
      integer lr
      integer lrline
      integer ltri
      integer m
      integer m1
      integer m2
      integer n
      double precision node_xy(dim_num,node_num)
      double precision r8_epsilon
      integer redg
      integer rtri
      integer stack(node_num)
      integer t
      double precision tol
      integer top
      integer triangle_neighbor(3,node_num*2)
      integer triangle_num
      integer triangle_node(3,node_num*2)

      tol = 100.0D+00 * r8_epsilon ( )

      ierr = 0
c
c  Sort the vertices by increasing (x,y).
c
      call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

      call r82vec_permute ( node_num, indx, node_xy )
c
c  Make sure that the data nodes are "reasonably" distinct.
c
      m1 = 1

      do i = 2, node_num

        m = m1
        m1 = i

        k = 0

        do j = 1, dim_num

          cmax = max ( abs ( node_xy(j,m) ), abs ( node_xy(j,m1) ) )

          if ( tol * ( cmax + 1.0D+00 ) 
     &         .lt. abs ( node_xy(j,m) - node_xy(j,m1) ) ) then
            k = j
            go to 10
          end if

        end do

10      continue

        if ( k .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
          write ( *, '(a,i8)' ) '  Fails for point number I = ', i
          write ( *, '(a,i8)' ) '  M = ', m
          write ( *, '(a,i8)' ) '  M1 = ', m1
          write ( *, '(a,2g14.6)' ) 
     &      '  NODE_XY(M)  = ', ( node_xy(dim,m), dim = 1, dim_num )
          write ( *, '(a,2g14.6)' ) 
     &      '  NODE_XY(M1) = ', ( node_xy(dim,m1), dim = 1, dim_num )
          ierr = 224
          stop
        end if

      end do
c
c  Starting from nodes M1 and M2, search for a third point M that
c  makes a "healthy" triangle (M1,M2,M)
c
      m1 = 1
      m2 = 2
      j = 3

20    continue

        if ( node_num .lt. j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
          ierr = 225
          stop
        end if

        m = j

        lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), 
     &    node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

        if ( lr .ne. 0 ) then
          go to 30
        end if

        j = j + 1

      go to 20

30    continue
c
c  Set up the triangle information for (M1,M2,M), and for any other
c  triangles you created because points were collinear with M1, M2.
c
      triangle_num = j - 2

      if ( lr .eq. -1 ) then

        triangle_node(1,1) = m1
        triangle_node(2,1) = m2
        triangle_node(3,1) = m
        triangle_neighbor(3,1) = -3

        do i = 2, triangle_num

          m1 = m2
          m2 = i+1

          triangle_node(1,i) = m1
          triangle_node(2,i) = m2
          triangle_node(3,i) = m

          triangle_neighbor(1,i-1) = -3 * i
          triangle_neighbor(2,i-1) = i
          triangle_neighbor(3,i) = i - 1

        end do

        triangle_neighbor(1,triangle_num) = -3 * triangle_num - 1
        triangle_neighbor(2,triangle_num) = -5
        ledg = 2
        ltri = triangle_num

      else

        triangle_node(1,1) = m2
        triangle_node(2,1) = m1
        triangle_node(3,1) = m

        triangle_neighbor(1,1) = -4

        do i = 2, triangle_num

          m1 = m2
          m2 = i+1

          triangle_node(1,i) = m2
          triangle_node(2,i) = m1
          triangle_node(3,i) = m

          triangle_neighbor(3,i-1) = i
          triangle_neighbor(1,i) = -3 * i - 3
          triangle_neighbor(2,i) = i - 1

        end do

        triangle_neighbor(3,triangle_num) = -3 * triangle_num
        triangle_neighbor(2,1) = -3 * triangle_num - 2
        ledg = 2
        ltri = 1

      end if
c
c  Insert the vertices one at a time from outside the convex hull,
c  determine visible boundary edges, and apply diagonal edge swaps until
c  Delaunay triangulation of vertices (so far) is obtained.
c
      top = 0

      do i = j+1, node_num

        m = i
        m1 = triangle_node(ledg,ltri)

        if ( ledg <= 2 ) then
          m2 = triangle_node(ledg+1,ltri)
        else
          m2 = triangle_node(1,ltri)
        end if

        lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), 
     &    node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

        if ( 0 .lt. lr ) then
          rtri = ltri
          redg = ledg
          ltri = 0
        else
          l = -triangle_neighbor(ledg,ltri)
          rtri = l / 3
          redg = mod ( l, 3 ) + 1
        end if

        call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, 
     &    triangle_num, triangle_node, triangle_neighbor, ltri, ledg, 
     &    rtri, redg )

        n = triangle_num + 1
        l = -triangle_neighbor(ledg,ltri)

40      continue

          t = l / 3
          e = mod ( l, 3 ) + 1
          l = -triangle_neighbor(e,t)
          m2 = triangle_node(e,t)

          if ( e <= 2 ) then
            m1 = triangle_node(e+1,t)
          else
            m1 = triangle_node(1,t)
          end if

          triangle_num = triangle_num + 1
          triangle_neighbor(e,t) = triangle_num

          triangle_node(1,triangle_num) = m1
          triangle_node(2,triangle_num) = m2
          triangle_node(3,triangle_num) = m

          triangle_neighbor(1,triangle_num) = t
          triangle_neighbor(2,triangle_num) = triangle_num - 1
          triangle_neighbor(3,triangle_num) = triangle_num + 1

          top = top + 1

          if ( node_num .lt. top ) then
            ierr = 8
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
            write ( *, '(a)' ) '  Stack overflow.'
            stop
          end if

          stack(top) = triangle_num

          if ( t .eq. rtri .and. e .eq. redg ) then
            go to 50
          end if

        go to 40

50      continue

        triangle_neighbor(ledg,ltri) = -3 * n - 1
        triangle_neighbor(2,n) = -3 * triangle_num - 2
        triangle_neighbor(3,triangle_num) = -l

        ltri = n
        ledg = 2

        call swapec ( m, top, ltri, ledg, node_num, node_xy, 
     &    triangle_num, triangle_node, triangle_neighbor, stack, 
     &    ierr )

        if ( ierr .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
          write ( *, '(a)' ) '  Error return from SWAPEC.'
          stop
        end if

      end do
c
c  Now account for the sorting that we did.
c
      do i = 1, 3
        do j = 1, triangle_num
          triangle_node(i,j) = indx ( triangle_node(i,j) )
        end do
      end do

      call perm_inverse ( node_num, indx )

      call r82vec_permute ( node_num, indx, node_xy )

      return
      end
      subroutine swapec ( i, top, btri, bedg, node_num, node_xy, 
     &  triangle_num, triangle_node, triangle_neighbor, stack, ierr )

c*********************************************************************72
c
cc SWAPEC swaps diagonal edges until all triangles are Delaunay.
c
c  Discussion:
c
c    The routine swaps diagonal edges in a 2D triangulation, based on
c    the empty circumcircle criterion, until all triangles are Delaunay,
c    given that I is the index of the new vertex added to the triangulation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, integer I, the index of the new vertex.
c
c    Input/output, integer TOP, the index of the top of the stack.
c    On output, TOP is zero.
c
c    Input/output, integer BTRI, BEDG; on input, if positive, are 
c    the triangle and edge indices of a boundary edge whose updated indices
c    must be recorded.  On output, these may be updated because of swaps.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input/output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the
c    triangle incidence list.  May be updated on output because of swaps.
c
c    Input/output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the
c    triangle neighbor list; negative values are used for links of the 
c    counter-clockwise linked list of boundary edges;  May be updated on output 
c    because of swaps.
c      LINK = -(3*I + J-1) where I, J = triangle, edge index.
c
c    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
c    contain the indices of initial triangles (involving vertex I)
c    put in stack; the edges opposite I should be in interior;  entries
c    TOP+1 through MAXST are used as a stack.
c
c    Output, integer IERR is set to 8 for abnormal return.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer triangle_num

      integer a
      integer b
      integer bedg
      integer btri
      integer c
      integer diaedg
      integer e
      integer ee
      integer em1
      integer ep1
      integer f
      integer fm1
      integer fp1
      integer i
      integer ierr
      integer i4_wrap
      integer l
      double precision node_xy(dim_num,node_num)
      integer r
      integer s
      integer stack(node_num)
      integer swap
      integer t
      integer top
      integer triangle_node(3,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      integer tt
      integer u
      double precision x
      double precision y
c
c  Determine whether triangles in stack are Delaunay, and swap
c  diagonal edge of convex quadrilateral if not.
c
      x = node_xy(1,i)
      y = node_xy(2,i)

10    continue

        if ( top .le. 0 ) then
          go to 40
        end if

        t = stack(top)
        top = top - 1

        if ( triangle_node(1,t) .eq. i ) then
          e = 2
          b = triangle_node(3,t)
        else if ( triangle_node(2,t) .eq. i ) then
          e = 3
          b = triangle_node(1,t)
        else
          e = 1
          b = triangle_node(2,t)
        end if

        a = triangle_node(e,t)
        u = triangle_neighbor(e,t)

        if ( triangle_neighbor(1,u) .eq. t ) then
          f = 1
          c = triangle_node(3,u)
        else if ( triangle_neighbor(2,u) .eq. t ) then
          f = 2
          c = triangle_node(1,u)
        else
          f = 3
          c = triangle_node(2,u)
        end if

        swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), 
     &    node_xy(2,c), node_xy(1,b), node_xy(2,b) )

        if ( swap .eq. 1 ) then

          em1 = i4_wrap ( e - 1, 1, 3 )
          ep1 = i4_wrap ( e + 1, 1, 3 )
          fm1 = i4_wrap ( f - 1, 1, 3 )
          fp1 = i4_wrap ( f + 1, 1, 3 )

          triangle_node(ep1,t) = c
          triangle_node(fp1,u) = i

          r = triangle_neighbor(ep1,t)
          s = triangle_neighbor(fp1,u)

          triangle_neighbor(ep1,t) = u
          triangle_neighbor(fp1,u) = t
          triangle_neighbor(e,t) = s
          triangle_neighbor(f,u) = r

          if ( 0 .lt. triangle_neighbor(fm1,u) ) then
            top = top + 1
            stack(top) = u
          end if

          if ( 0 .lt. s ) then

            if ( triangle_neighbor(1,s) .eq. u ) then
              triangle_neighbor(1,s) = t
            else if ( triangle_neighbor(2,s) .eq. u ) then
              triangle_neighbor(2,s) = t
            else
              triangle_neighbor(3,s) = t
            end if

            top = top + 1

            if ( node_num .lt. top ) then
              ierr = 8
              return
            end if

            stack(top) = t

          else

            if ( u .eq. btri .and. fp1 .eq. bedg ) then
              btri = t
              bedg = e
            end if

            l = - ( 3 * t + e - 1 )
            tt = t
            ee = em1

20          continue

            if ( 0 < triangle_neighbor(ee,tt) ) then

              tt = triangle_neighbor(ee,tt)

              if ( triangle_node(1,tt) .eq. a ) then
                ee = 3
              else if ( triangle_node(2,tt) .eq. a ) then
                ee = 1
              else
                ee = 2
              end if

              go to 20

            end if

            triangle_neighbor(ee,tt) = l

          end if

          if ( 0 .lt. r ) then

            if ( triangle_neighbor(1,r) .eq. t ) then
              triangle_neighbor(1,r) = u
            else if ( triangle_neighbor(2,r) .eq. t ) then
              triangle_neighbor(2,r) = u
            else
              triangle_neighbor(3,r) = u
            end if

          else

            if ( t .eq. btri .and. ep1 .eq. bedg ) then
              btri = u
              bedg = f
            end if

            l = - ( 3 * u + f - 1 )
            tt = u
            ee = fm1

30          continue

            if ( 0 .lt. triangle_neighbor(ee,tt) ) then

              tt = triangle_neighbor(ee,tt)

              if ( triangle_node(1,tt) .eq. b ) then
                ee = 3
              else if ( triangle_node(2,tt) .eq. b ) then
                ee = 1
              else
                ee = 2
              end if

              go to 30

            end if

            triangle_neighbor(ee,tt) = l

          end if

        end if

      go to 10

40    continue

      return
      end
      subroutine triangulation_order3_print ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_PRINT prints information about a triangulation.
c
c  Discussion:
c
c    Triangulations created by R8TRIS2 include extra information encoded
c    in the negative values of TRIANGLE_NEIGHBOR.
c
c    Because some of the nodes counted in NODE_NUM may not actually be
c    used in the triangulation, I needed to compute the true number
c    of vertices.  I added this calculation on 13 October 2001.
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up the triangles.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbors on each side.  If there is no triangle neighbor on 
c    a particular side, the value of TRIANGLE_NEIGHBOR should be negative.  
c    If the triangulation data was created by R8TRIS22, then there is more 
c    information encoded in the negative values.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer boundary_num
      integer i
      integer i4_wrap
      integer j
      integer k
      integer n1
      integer n2
      double precision node_xy(dim_num,node_num)
      integer s
      logical skip
      integer t
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      integer vertex_list(triangle_order*triangle_num)
      integer vertex_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PRINT'
      write ( *, '(a)' ) 
     &  '  Information defining an order3 triangulation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num

      call r8mat_transpose_print ( dim_num, node_num, node_xy, 
     &  '  Node coordinates' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of triangles is ', triangle_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Sets of three nodes are used as vertices of'
      write ( *, '(a)' ) 
     &  '  the triangles.  For each triangle, the nodes'
      write ( *, '(a)' ) '  are listed in counterclockwise order.'

      call i4mat_transpose_print ( 3, triangle_num, triangle_node, 
     &  '  Triangle nodes:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  On each side of a given triangle, there is either'
      write ( *, '(a)' ) 
     &  '  another triangle, or a piece of the convex hull.'
      write ( *, '(a)' ) 
     &  '  For each triangle, we list the indices of the three'
      write ( *, '(a)' ) 
     &  '  neighbors, or (if negative) the codes of the'
      write ( *, '(a)' ) '  segments of the convex hull.'

      call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor, 
     &  '  Triangle neighbors' )
c
c  Determine the number of vertices.  
c
      k = 0
      do j = 1, triangle_num
        do i = 1, 3
          k = k + 1
          vertex_list(k) = triangle_node(i,j)
        end do
      end do

      call i4vec_sort_heap_a ( 3*triangle_num, vertex_list )

      call i4vec_sorted_unique ( 3*triangle_num, vertex_list, 
     &  vertex_num )
c
c  Determine the number of boundary points.
c
      boundary_num = 2 * vertex_num - triangle_num - 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of boundary points is ', boundary_num
       
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The segments that make up the convex hull can be'
      write ( *, '(a)' ) 
     &  '  determined from the negative entries of the triangle'
      write ( *, '(a)' ) '  neighbor list.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #   Tri  Side    N1    N2'
      write ( *, '(a)' ) ' '

      skip = .false.

      k = 0

      do i = 1, triangle_num

        do j = 1, 3

          if ( triangle_neighbor(j,i) .lt. 0 ) then
            s = - triangle_neighbor(j,i)
            t = s / 3

            if ( t .lt. 1 .or. triangle_num .lt. t ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 
     &          '  Sorry, this data does not use the R8TRIS2'
              write ( *, '(a)' ) 
     &          '  convention for convex hull segments.'
              skip = .true.
              exit
            end if

            s = mod ( s, 3 ) + 1
            k = k + 1
            n1 = triangle_node(s,t)
            n2 = triangle_node(i4_wrap(s+1,1,3),t)
            write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) 
     &        k, t, s, n1, n2
          end if

        end do

        if ( skip ) then
          exit
        end if

      end do

      return
      end
      subroutine triangulation_search_delaunay ( node_num, node_xy, 
     &  element_order, element_num, element_node, element_neighbor, 
     &  p, triangle_index, alpha, beta, gamma, edge, step_num )

c*********************************************************************72
c
cc TRIANGULATION_SEARCH_DELAUNAY searches a Delaunay triangulation for a point.
c
c  Discussion:
c
c    The algorithm "walks" from one triangle to its neighboring triangle,
c    and so on, until a triangle is found containing point P, or P is found
c    to be outside the convex hull.
c
c    The algorithm computes the barycentric coordinates of the point with
c    respect to the current triangle.  If all three quantities are positive,
c    the point is contained in the triangle.  If the I-th coordinate is
c    negative, then P lies on the far side of edge I, which is opposite
c    from vertex I.  This gives a hint as to where to search next.
c
c    For a Delaunay triangulation, the search is guaranteed to terminate.
c    For other triangulations, a cycle may occur.
c
c    Note the surprising fact that, even for a Delaunay triangulation of
c    a set of nodes, the nearest node to P need not be one of the
c    vertices of the triangle containing P.
c
c    The code can be called for triangulations of any order, but only
c    the first three nodes in each triangle are considered.  Thus, if
c    higher order triangles are used, and the extra nodes are intended
c    to give the triangle a polygonal shape, these will have no effect,
c    and the results obtained here might be misleading.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2012
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer ELEMENT_ORDER, the order of the triangles.
c
c    Input, integer ELEMENT_NUM, the number of triangles.
c
c    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
c    the nodes that make up each triangle.
c
c    Input, integer ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
c    triangle neighbor list.
c
c    Input, double precision P(2), the coordinates of a point.
c
c    Output, integer TRIANGLE_INDEX, the index of the triangle
c    where the search ended.  If a cycle occurred, then TRIANGLE_INDEX = -1.
c
c    Output, double precision ALPHA, BETA, GAMMA, the barycentric
c    coordinates of the point relative to triangle TRIANGLE_INDEX.
c
c    Output, integer EDGE, indicates the position of the point P in
c    triangle TRIANGLE_INDEX:
c    0, the interior or boundary of the triangle;
c    -1, outside the convex hull of the triangulation, past edge 1;
c    -2, outside the convex hull of the triangulation, past edge 2;
c    -3, outside the convex hull of the triangulation, past edge 3.
c
c    Output, integer STEP_NUM, the number of steps.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer element_num
      integer element_order

      integer a
      double precision alpha
      integer b
      double precision beta
      integer c
      double precision det
      double precision dxp
      double precision dxa
      double precision dxb
      double precision dyp
      double precision dya
      double precision dyb
      integer edge
      double precision gamma
      double precision node_xy(dim_num,node_num)
      double precision p(dim_num)
      integer step_num
      integer element_node(element_order,element_num)
      integer triangle_index
      integer triangle_index_save
      integer element_neighbor(3,element_num)

      save triangle_index_save

      data triangle_index_save / -1 /
c
c  If possible, start with the previous successful value of TRIANGLE_INDEX.
c
      if ( triangle_index_save .lt. 1 .or. 
     &  element_num .lt. triangle_index_save ) then
        triangle_index = ( element_num + 1 ) / 2
      else
        triangle_index = triangle_index_save
      end if

      step_num = - 1
      edge = 0

10    continue

        step_num = step_num + 1

        if ( element_num .lt. step_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      'TRIANGULATION_SEARCH_DELAUNAY - Fatal error!'
          write ( *, '(a)' ) '  The algorithm seems to be cycling.'
          triangle_index = -1
          edge = -1
          stop
        end if
c
c  Get the nodes of triangle TRIANGLE_INDEX.
c
        a = element_node(1,triangle_index)
        b = element_node(2,triangle_index)
        c = element_node(3,triangle_index)
c
c  Using vertex C as a base, compute the distances to vertices A and B,
c  and the point P.
c
        dxa = node_xy(1,a) - node_xy(1,c)
        dya = node_xy(2,a) - node_xy(2,c)

        dxb = node_xy(1,b) - node_xy(1,c)
        dyb = node_xy(2,b) - node_xy(2,c)

        dxp = p(1)         - node_xy(1,c)
        dyp = p(2)         - node_xy(2,c)

        det = dxa * dyb - dya * dxb
c
c  Compute the barycentric coordinates of the point P with respect
c  to this triangle.
c
        alpha = ( dxp * dyb - dyp * dxb ) / det
        beta =  ( dxa * dyp - dya * dxp ) / det
        gamma = 1.0D+00 - alpha - beta
c
c  If the barycentric coordinates are all positive, then the point
c  is inside the triangle and we're done.
c
        if ( 0.0D+00 .le. alpha .and. 
     &       0.0D+00 .le. beta  .and. 
     &       0.0D+00 .le. gamma ) then
          go to 20
        end if
c
c  At least one barycentric coordinate is negative.
c
c  If there is a negative barycentric coordinate for which there exists
c  an opposing triangle neighbor closer to the point, move to that triangle.
c
c  (Two coordinates could be negative, in which case we could go for the
c  most negative one, or the most negative one normalized by the actual
c  distance it represents).
c
        if ( alpha .lt. 0.0D+00 .and. 
     &    0 .lt. element_neighbor(2,triangle_index) ) then
          triangle_index = element_neighbor(2,triangle_index)
          go to 10
        else if ( beta .lt. 0.0D+00 .and. 
     &    0 .lt. element_neighbor(3,triangle_index) ) then
          triangle_index = element_neighbor(3,triangle_index)
          go to 10
        else if ( gamma .lt. 0.0D+00 .and. 
     &    0 .lt. element_neighbor(1,triangle_index) ) then
          triangle_index = element_neighbor(1,triangle_index)
          go to 10
        end if
c
c  All negative barycentric coordinates correspond to vertices opposite
c  sides on the convex hull.
c
c  Note the edge and exit.
c
        if ( alpha .lt. 0.0D+00 ) then
          edge = -2
          go to 20
        else if ( beta .lt. 0.0D+00 ) then
          edge = -3
          go to 20
        else if ( gamma .lt. 0.0D+00 ) then
          edge = -1
          go to 20
        end if

      go to 10

20    continue

      triangle_index_save = triangle_index

      return
      end
      subroutine vbedg ( x, y, node_num, node_xy, triangle_num, 
     &  triangle_node, triangle_neighbor, ltri, ledg, rtri, redg )

c*********************************************************************72
c
cc VBEDG determines which boundary edges are visible to a point.
c
c  Discussion:
c
c    The point (X,Y) is assumed to be outside the convex hull of the
c    region covered by the 2D triangulation.
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
c    Original FORTRAN77 version by Barry Joe.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point outside the
c    convex hull of the current triangulation.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the triangle 
c    incidence list.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbor list; negative values are used for links of a 
c    counter clockwise linked list of boundary edges;
c      LINK = -(3*I + J-1) where I, J = triangle, edge index.
c
c    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these 
c    values are assumed to be already computed and are not changed, else they 
c    are updated.  On output, LTRI is the index of boundary triangle to the 
c    left of the leftmost boundary triangle visible from (X,Y), and LEDG is 
c    the boundary edge of triangle LTRI to the left of the leftmost boundary 
c    edge visible from (X,Y).  1 <= LEDG <= 3.
c
c    Input/output, integer RTRI.  On input, the index of the 
c    boundary triangle to begin the search at.  On output, the index of the 
c    rightmost boundary triangle visible from (X,Y).
c
c    Input/output, integer REDG, the edge of triangle RTRI that 
c    is visible from (X,Y).  1 <= REDG <= 3.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer triangle_num

      integer a
      integer b
      integer e
      integer i4_wrap
      integer l
      logical ldone
      integer ledg
      integer lr
      integer lrline
      integer ltri
      double precision node_xy(2,node_num)
      integer redg
      integer rtri
      integer t
      integer triangle_node(3,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      double precision x
      double precision y
c
c  Find the rightmost visible boundary edge using links, then possibly
c  leftmost visible boundary edge using triangle neighbor information.
c
      if ( ltri .eq. 0 ) then
        ldone = .false.
        ltri = rtri
        ledg = redg
      else
        ldone = .true.
      end if

10    continue

        l = - triangle_neighbor(redg,rtri)
        t = l / 3
        e = mod ( l, 3 ) + 1
        a = triangle_node(e,t)

        if ( e .le. 2 ) then
          b = triangle_node(e+1,t)
        else
          b = triangle_node(1,t)
        end if

        lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), 
     &    node_xy(2,b), 0.0D+00 )

        if ( lr .le. 0 ) then
          go to 20
        end if

        rtri = t
        redg = e

      go to 10

20    continue

      if ( ldone ) then
        return
      end if

      t = ltri
      e = ledg

30    continue

        b = triangle_node(e,t)
        e = e - 1
        e = i4_wrap ( e, 1, 3 )

40      continue

        if ( 0 .lt. triangle_neighbor(e,t) ) then

          t = triangle_neighbor(e,t)

          if ( triangle_node(1,t) .eq. b ) then
            e = 3
          else if ( triangle_node(2,t) .eq. b ) then
            e = 1
          else
            e = 2
          end if

          go to 40

        end if

        a = triangle_node(e,t)

        lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), 
     &    node_xy(2,b), 0.0D+00 )

        if ( lr .le. 0 ) then
          go to 50
        end if

      go to 30

50    continue

      ltri = t
      ledg = e

      return
      end
