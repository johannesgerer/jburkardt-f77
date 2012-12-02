      subroutine bin_search_one_2d ( bin, nset, pset, nbin, bin_start,
     &  bin_next, ptest, found_a_neighbor, i_min, dist_min_sq,
     &  compares )

c*********************************************************************72
c
cc BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
c
c  Discussion:
c
c    The bins are presumed to have been set up by successive calls to:
c
c      R82VEC_BIN_EVEN2,
c      R82VEC_BINNED_REORDER, and
c      R82VEC_BINNED_SORT_A.
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
c    Input, integer BIN(2), the indices of the cell to be examined.
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the coordinates of the points
c    in the set.
c
c    Input, integer NBIN(2), the number of cells in the horizontal and
c    vertical directions.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    indicates the index of the first and last element in the bin, or -1
c    if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(2), the coordinates of the test point.
c
c    Input/output, logical FOUND_A_NEIGHBOR, is set to TRUE if at least
c    one point of PSET is found in the current bin.  Otherwise, it retains its
c    input value.
c
c    Input/output, integer I_MIN, the index of the nearest neighbor in
c    PSET to PTEST, if at least one neighbor has been found.
c
c    Input/output, double precision DIST_MIN_SQ, the square of the distance
c    from the nearest neighbor in PSET to PTEST, if at least one neighbor
c    has been found.
c
c    Input/output, integer COMPARES, the number of elements of PSET whose
c    distance to PTEST has been computed.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer nbin(ndim)
      integer nset

      integer bin(ndim)
      integer bin_next(nset)
      integer bin_start(nbin(1),nbin(2))
      integer compares
      double precision dist_min_sq
      double precision dist_sq
      logical found_a_neighbor
      integer i
      integer i_min
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)

      node = bin_start(bin(1),bin(2))

10    continue

      if ( 0 .lt. node ) then

        found_a_neighbor = .true.

        dist_sq = 0.0
        do i = 1, ndim
          dist_sq = dist_sq + ( ptest(i) - pset(i,node) )**2
        end do

        compares = compares + 1

        if ( dist_sq .lt. dist_min_sq ) then
          dist_min_sq = dist_sq
          i_min = node
        end if

        node = bin_next(node)

        go to 10

      end if

      return
      end
      subroutine bin_to_r8_even ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R8_EVEN returns the limits for a given "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
c    An initial bin takes everything less than A, and a final bin takes
c    everything greater than B.
c
c  Example:
c
c    NBIN = 7, A = 10, B = 20
c
c    BIN      CMIN  CMAX
c
c    1         -HUGE 10
c    2         10    12
c    3         12    14
c    4         14    16
c    5         16    18
c    6         18    20
c    7         20    HUGE
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
c    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
c    If NBIN is 1 or 2, then everything is assigned to bin 1.
c
c    Input, integer BIN, the index of the bin to be considered.
c    If BIN is less than 1, or greater than NBIN, the user will get what
c    the user deserves.
c
c    Input, double precision A, B, the lower and upper limits of the
c    bin interval.  While A is expected to be less than B, the code
c    should return useful results if A is actually greater than B.
c
c    Output, double precision CMIN, CMAX, the minimum and maximum
c    limits on the bin.
c
      implicit none

      double precision a
      double precision b
      integer bin
      double precision cmax
      double precision cmin
      integer nbin
      double precision r8_huge
c
c  Take care of special cases.
c
      if ( nbin .le. 2 ) then
        cmin = - r8_huge ( )
        cmax =   r8_huge ( )
        return
      end if

      if ( b .eq. a ) then
        cmin = - r8_huge ( )
        cmax =   r8_huge ( )
        return
      end if
c
c  Compute the bin limits.
c
      if ( bin .eq. 1 ) then
        cmin = - r8_huge ( )
        cmax = a
      else if ( bin .lt. nbin ) then

        cmin = ( dble ( nbin - bin     ) * a
     &         + dble (        bin - 2 ) * b )
     &         / dble ( nbin       - 2 )

        cmax = ( dble ( nbin - bin - 1 ) * a
     &         + dble (        bin - 1 ) * b )
     &         / dble ( nbin       - 2 )

      else if ( bin .eq. nbin ) then
        cmin = b
        cmax = r8_huge ( )
      else
        cmin = - r8_huge ( )
        cmax =   r8_huge ( )
      end if

      return
      end
      subroutine bin_to_r8_even2 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R8_EVEN2 returns the limits for a given "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
c
c  Example:
c
c    NBIN = 5, A = 10, B = 20
c
c    BIN      CMIN  CMAX
c
c    1         10    12
c    2         12    14
c    3         14    16
c    4         16    18
c    5         18    20
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins.
c
c    Input, integer BIN, the index of the bin to be considered.
c    If BIN is less than 1, or greater than NBIN, the user will get what
c    the user deserves.
c
c    Input, double precision A, B, the lower and upper limits of the bin
c    interval.  While A is expected to be less than B, the code should
c    return useful results if A is actually greater than B.
c
c    Output, double precision CMIN, CMAX, the minimum and maximum limits
c    on the bin.
c
      implicit none

      double precision a
      double precision b
      integer bin
      double precision cmax
      double precision cmin
      integer nbin
      double precision r8_huge
c
c  Compute the bin limits.
c
      if ( bin .lt. 1 ) then
        cmin = - r8_huge ( )
        cmax = a
      else if ( bin .le. nbin ) then
        cmin = ( dble ( nbin - bin + 1 ) * a
     &         + dble (        bin - 1 ) * b )
     &         / dble ( nbin           )
        cmax = ( dble ( nbin - bin     ) * a
     &         + dble (        bin     ) * b )
     &         / dble ( nbin           )
      else if ( nbin .le. bin ) then
        cmin = b
        cmax = r8_huge ( )
      end if

      return
      end
      subroutine bin_to_r82_even ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R82_EVEN returns the limits for a given R82 "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
c    An initial bin takes everything less than A, and a final bin takes
c    everything greater than B.
c
c  Example:
c
c    NBIN = 7, A(1) = 5, B(1) = 15
c              A(2) = 0, B(2) = 20
c
c     BIN         CMIN      CMAX
c    ------   -----------  --------
c    1, 1     -HUGE -HUGE   5     0
c    2, 2       5     0     7     4
c    3, 3       7     4     9     8
c    4, 4       9     8    11    12
c    5, 5      11    12    13    16
c    6, 6      13    16    15    20
c    7, 7      15    20    HUGE HUGE
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
c    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
c    If NBIN is 1 or 2, then everything is assigned to bin 1.
c
c    Input, integer BIN(2), the index of the bin to be considered.
c    If BIN(I) is less than 1, or greater than NBIN, the user will get what
c    the user deserves.
c
c    Input, double precision A(2), B(2), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(2), CMAX(2), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call bin_to_r8_even ( nbin, bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      subroutine bin_to_r82_even2 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R82_EVEN2 returns the limits for a given D82 "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
c
c  Example:
c
c    NBIN = 5, A(1) = 5, B(1) = 15
c              A(2) = 0, B(2) = 20
c
c     BIN         CMIN      CMAX
c    ------   -----------  --------
c    1, 1       5     0     7     4
c    2, 2       7     4     9     8
c    3, 3       9     8    11    12
c    4, 4      11    12    13    16
c    5, 5      13    16    15    20
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
c    Input, integer NBIN, the number of bins.
c
c    Input, integer BIN(2), the index of the bin to be considered.
c
c    Input, double precision A(2), B(2), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(2), CMAX(2), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin, bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      subroutine bin_to_r82_even3 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R82_EVEN3 returns the limits for a given R82 "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A(I) to B(I) is divided into NBIN(I) equal
c    subintervals or bins.
c
c  Example:
c
c    NBIN = (/ 4, 5, /)
c
c    A(1) = 5, B(1) = 15
c    A(2) = 0, B(2) = 20
c
c     BIN         CMIN      CMAX
c    ------   -----------  --------
c    1, 1       5     0     7     4
c    2, 2       7     4     9     8
c    3, 3       9     8    11    12
c    4, 4      11    12    13    16
c    5, 5      13    16    15    20
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
c    Input, integer NBIN(2), the number of bins in each dimension.
c
c    Input, integer BIN(2), the index of the bin to be considered.
c
c    Input, double precision A(2), B(2), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(2), CMAX(2), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      subroutine bin_to_r83_even2 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R83_EVEN2 returns the limits for a given R83 "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
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
c    Input, integer NBIN, the number of bins.
c
c    Input, integer BIN(), the index of the bin to be considered.
c
c    Input, double precision A(3), B(3), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(3), CMAX(3), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin, bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      subroutine bin_to_r83_even3 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R83_EVEN3 returns the limits for a given R83 "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A(I) to B(I) is divided into NBIN(I) equal
c    subintervals or bins.
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
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input, integer BIN(3), the index of the bin to be considered.
c
c    Input, double precision A(3), B(3), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(3), CMAX(3), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      subroutine bin_to_r8vec_even3 ( ndim, nbin, bin, a, b, cmin,
     &  cmax )

c*********************************************************************72
c
cc BIN_TO_R8VEC_EVEN3 returns the limits for a given R8VEC "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A(I) to B(I) is divided into NBIN(I) equal
c    subintervals or bins.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDIM, the dimension of the space.
c
c    Input, integer NBIN(NDIM), the number of bins in each dimension.
c
c    Input, integer BIN(NDIM), the index of the bin to be considered.
c
c    Input, double precision A(NDIM), B(NDIM), the lower and upper limits
c    of the bin interval.  While A(I) is expected to be less than B(I), the
c    code should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(NDIM), CMAX(NDIM), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i),
     &    cmax(i) )
      end do

      return
      end
      function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

c*********************************************************************72
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
c    25 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This version by John Burkardt.
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

      tola = tol * max (
     &  abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )

      tolb = tol * max (
     &  abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

      ca = dx10 * dx30 + dy10 * dy30
      cb = dx12 * dx32 + dy12 * dy32

      if ( tola .lt. ca .and. tolb .lt. cb ) then

        diaedg = - 1

      else if ( ca .lt. -tola .and. cb .lt. -tolb ) then

        diaedg = 1

      else

        tola = max ( tola, tolb )
        s = ( dx10 * dy30 - dx30 * dy10 ) * cb
     &    + ( dx32 * dy12 - dx12 * dy32 ) * ca

        if ( tola .lt. s ) then
          diaedg = -1
        else if ( s .lt. - tola ) then
          diaedg = 1
        else
          diaedg = 0
        end if

      end if

      return
      end
      subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert,
     &  tri_nabe )

c*********************************************************************72
c
cc DTRIS2 constructs a Delaunay triangulation of 2D vertices.
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
c    29 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This version by John Burkardt.
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
c    Input, integer POINT_NUM, the number of vertices.
c
c    Input/output, double precision POINT_XY(2,POINT_NUM), the coordinates
c    of the vertices.  On output, the vertices have been sorted into
c    dictionary order.
c
c    Output, integer TRI_NUM, the number of triangles in the triangulation;
c    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
c    of boundary vertices.
c
c    Output, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
c    The elements are indices of POINT_XY.  The vertices of the triangles are
c    in counter clockwise order.
c
c    Output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list.
c    Positive elements are indices of TIL; negative elements are used for links
c    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
c    where I, J = triangle, edge index; TRI_NABE(J,I) refers to
c    the neighbor along edge from vertex J to J+1 (mod 3).
c
      implicit none

      integer point_num

      double precision cmax
      integer e
      integer i
      integer ierr
      integer indx(point_num)
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
      double precision point_xy(2,point_num)
      double precision r8_epsilon
      integer redg
      integer rtri
      integer stack(point_num)
      integer t
      double precision tol
      integer top
      integer tri_nabe(3,point_num*2)
      integer tri_num
      integer tri_vert(3,point_num*2)

      tol = 100.0D+00 * r8_epsilon ( )

      ierr = 0
c
c  Sort the vertices by increasing (x,y).
c
      call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

      call r82vec_permute ( point_num, point_xy, indx )
c
c  Make sure that the data points are "reasonably" distinct.
c
      m1 = 1

      do i = 2, point_num

        m = m1
        m1 = i

        k = 0

        do j = 1, 2

          cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

          if ( tol * ( cmax + 1.0D+00 )
     &         .lt. abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
            k = j
            go to 10
          end if

        end do

10      continue

        if ( k .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a,i6)' ) '  Fails for point number I = ', i
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  M1 = ', m1
          write ( *, '(a,2g14.6)' )
     &      '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
          write ( *, '(a,2g14.6)' )
     &      '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
          ierr = 224
          return
        end if

      end do
c
c  Starting from points M1 and M2, search for a third point M that
c  makes a "healthy" triangle (M1,M2,M)
c
      m1 = 1
      m2 = 2
      j = 3

20    continue

        if ( point_num .lt. j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          ierr = 225
          return
        end if

        m = j

        lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1),
     &    point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

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
      tri_num = j - 2

      if ( lr .eq. -1 ) then

        tri_vert(1,1) = m1
        tri_vert(2,1) = m2
        tri_vert(3,1) = m
        tri_nabe(3,1) = -3

        do i = 2, tri_num

          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m1
          tri_vert(2,i) = m2
          tri_vert(3,i) = m
          tri_nabe(1,i-1) = -3 * i
          tri_nabe(2,i-1) = i
          tri_nabe(3,i) = i - 1

        end do

        tri_nabe(1,tri_num) = -3 * tri_num - 1
        tri_nabe(2,tri_num) = -5
        ledg = 2
        ltri = tri_num

      else

        tri_vert(1,1) = m2
        tri_vert(2,1) = m1
        tri_vert(3,1) = m
        tri_nabe(1,1) = -4

        do i = 2, tri_num
          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m2
          tri_vert(2,i) = m1
          tri_vert(3,i) = m
          tri_nabe(3,i-1) = i
          tri_nabe(1,i) = -3 * i - 3
          tri_nabe(2,i) = i - 1
        end do

        tri_nabe(3,tri_num) = -3 * tri_num
        tri_nabe(2,1) = -3 * tri_num - 2
        ledg = 2
        ltri = 1

      end if
c
c  Insert the vertices one at a time from outside the convex hull,
c  determine visible boundary edges, and apply diagonal edge swaps until
c  Delaunay triangulation of vertices (so far) is obtained.
c
      top = 0

      do i = j+1, point_num

        m = i
        m1 = tri_vert(ledg,ltri)

        if ( ledg .le. 2 ) then
          m2 = tri_vert(ledg+1,ltri)
        else
          m2 = tri_vert(1,ltri)
        end if

        lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1),
     &    point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

        if ( 0 .lt. lr ) then
          rtri = ltri
          redg = ledg
          ltri = 0
        else
          l = -tri_nabe(ledg,ltri)
          rtri = l / 3
          redg = mod(l,3) + 1
        end if

        call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy,
     &    tri_num, tri_vert, tri_nabe, ltri, ledg, rtri, redg )

        n = tri_num + 1
        l = -tri_nabe(ledg,ltri)

50      continue

          t = l / 3
          e = mod ( l, 3 ) + 1
          l = -tri_nabe(e,t)
          m2 = tri_vert(e,t)

          if ( e .le. 2 ) then
            m1 = tri_vert(e+1,t)
          else
            m1 = tri_vert(1,t)
          end if

          tri_num = tri_num + 1
          tri_nabe(e,t) = tri_num
          tri_vert(1,tri_num) = m1
          tri_vert(2,tri_num) = m2
          tri_vert(3,tri_num) = m
          tri_nabe(1,tri_num) = t
          tri_nabe(2,tri_num) = tri_num - 1
          tri_nabe(3,tri_num) = tri_num + 1
          top = top + 1

          if ( point_num .lt. top ) then
            ierr = 8
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
            write ( *, '(a)' ) '  Stack overflow.'
            return
          end if

          stack(top) = tri_num

          if ( t .eq. rtri .and. e .eq. redg ) then
            go to 60
          end if

        go to 50

60      continue

        tri_nabe(ledg,ltri) = -3 * n - 1
        tri_nabe(2,n) = -3 * tri_num - 2
        tri_nabe(3,tri_num) = -l
        ltri = n
        ledg = 2

        call swapec ( m, top, ltri, ledg, point_num, point_xy,
     &    tri_num, tri_vert, tri_nabe, stack, ierr )

        if ( ierr .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a)' ) '  Error return from SWAPEC.'
          return
        end if

      end do
c
c  Now account for the sorting that we did.
c
      do i = 1, 3
        do j = 1, tri_num
          tri_vert(i,j) = indx ( tri_vert(i,j) )
        end do
      end do

      call perm_inv ( point_num, indx )

      call r82vec_permute ( point_num, point_xy, indx )

      return
      end
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ) then
        seed = seed - 1
      end if

      return
      end
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

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
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
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
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
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
c    30 June 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the matrix to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      character*(*) title

      ilo = 1
      ihi = m
      jlo = 1
      jhi = n

      call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

      return
      end
      subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc I4MAT_PRINT_SOME prints some of an I4MAT.
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
c    04 November 2003
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
c    Input, character*(*) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
      character*(8) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      integer s_len_trim
      character*(*) title

      if ( 0 .lt. s_len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8)' ) j
        end do

        write ( *, '(''  Col '',10a8)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

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
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
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
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
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
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer s_len_trim
      character*(*) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,i12)' ) i, a(i)
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
      subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

c*********************************************************************72
c
cc I4VEC2_COMPARE compares pairs of integers stored in two vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data items.
c
c    Input, integer A1(N), A2(N), contain the two components of each item.
c
c    Input, integer I, J, the items to be compared.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, item I < item J,
c     0, item I = item J,
c    +1, item J < item I.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer isgn
      integer j

      isgn = 0

           if ( a1(i) .lt. a1(j) ) then

        isgn = -1

      else if ( a1(i) .eq. a1(j) ) then

             if ( a2(i) .lt. a2(j) ) then
          isgn = -1
        else if ( a2(i) .lt. a2(j) ) then
          isgn = 0
        else if ( a2(j) .lt. a2(i) ) then
         isgn = +1
        end if

      else if ( a1(j) .lt. a1(i) ) then

        isgn = +1

      end if

      return
      end
      subroutine i4vec2_sort_a ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
c
c  Discussion:
c
c    Each item to be sorted is a pair of integers (I,J), with the I
c    and J values stored in separate vectors A1 and A2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items of data.
c
c    Input/output, integer A1(N), A2(N), the data to be sorted..
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer indx
      integer isgn
      integer j
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4_swap ( a1(i), a1(j) )
          call i4_swap ( a2(i), a2(j) )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4vec2_compare ( n, a1, a2, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec2_sorted_unique ( n, a1, a2, nuniq )

c*********************************************************************72
c
cc I4VEC2_SORTED_UNIQUE keeps the unique elements in a array of pairs of integers.
c
c  Discussion:
c
c    Item I is stored as the pair A1(I), A2(I).
c
c    The items must have been sorted, or at least it must be the
c    case that equal items are stored in adjacent vector locations.
c
c    If the items were not sorted, then this routine will only
c    replace a string of equal values by a single representative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items.
c
c    Input/output, integer A1(N), A2(N).
c    On input, the array of N items.
c    On output, an array of NUNIQ unique items.
c
c    Output, integer NUNIQ, the number of unique items.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer itest
      integer nuniq

      nuniq = 0

      if ( n .le. 0 ) then
        return
      end if

      nuniq = 1

      do itest = 2, n

        if ( a1(itest) .ne. a1(nuniq) .or.
     &       a2(itest) .ne. a2(nuniq) ) then

          nuniq = nuniq + 1

          a1(nuniq) = a1(itest)
          a2(nuniq) = a2(itest)

        end if

      end do

      return
      end
      subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
c
c  Discussion:
c
c    The box has center at (IC,JC), and has half-widths N1 and N2.
c    The indices are exactly those which are between (IC-N1,JC-N2) and
c    (IC+N1,JC+N2) with the property that at least one of I and J
c    is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the half-widths of the box, that is, the
c    maximum distance allowed between (IC,JC) and (I,J).
c
c    Input, integer IC, JC, the central cell of the box.
c
c    Input/output, integer I, J.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I and J are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      logical more
      integer n1
      integer n2

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        return
      end if

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 ) then
        more = .false.
        return
      end if
c
c  Increment J.
c
      j = j + 1
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and.
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k,
     &  more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
c
c  Discussion:
c
c    The box has a central cell of (IC,JC,KC), with a half widths of
c    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3)
c    and (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J,
c    and K is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
c    maximum distances from the central cell allowed for I, J and K.
c
c    Input, integer IC, JC, KC, the central cell of the box.
c
c    Input/output, integer I, J, K.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I, J, and K are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      integer k
      integer kc
      logical more
      integer n1
      integer n2
      integer n3

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        k = kc - n3
        return
      end if

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 .and.
     &  k .eq. kc + n3 ) then
        more = .false.
        return
      end if
c
c  Increment K.
c
      k = k + 1
c
c  Check K.
c
      if ( kc + n3 .lt. k ) then
        k = kc - n3
        j = j + 1
      else if ( k .lt. kc + n3 .and.
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or.
     &    j .eq. jc - n2 .or. j .eq. jc + n2 ) ) then
        return
      else
        k = kc + n3
        return
      end if
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and.
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or.
     &    k .eq. kc - n3 .or. k .eq. kc + n3 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

c*********************************************************************72
c
cc LRLINE determines where a point lies in relation to a directed line.
c
c  Discussion:
c
c    LRLINE determines whether a point is to the left of, right of,
c    or on a directed line parallel to a line through given points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This version by John Burkardt.
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
c    Input, double precision XU, YU, XV1, YV1, XV2, YV2, are vertex
c    coordinates; the directed line is parallel to and at signed distance
c    DV to the left of the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU)
c    is the vertex for which the position relative to the directed line is
c    to be determined.
c
c    Input, double precision DV, the signed distance, positive for left.
c
c    Output, integer LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
c    to the right of, on, or left of the directed line.  LRLINE is 0 if
c    the line degenerates to a point.
c
      implicit none

      double precision dv
      double precision dx
      double precision dxu
      double precision dy
      double precision dyu
      integer lrline
      double precision t
      double precision tol
      parameter ( tol = 0.0000001D+00 )
      double precision tolabs
      double precision xu
      double precision xv1
      double precision xv2
      double precision yu
      double precision yv1
      double precision yv2

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
      else if ( t .lt. -tolabs ) then
        lrline = -1
      end if

      return
      end
      subroutine perm_inv ( n, p )

c*********************************************************************72
c
cc PERM_INV inverts a permutation "in place".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
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

      integer i
      integer i0
      integer i1
      integer i2
      integer is
      integer p(n)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INV - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        stop
      end if

      is = 1

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          i2 = p(i1)
          p(i1) = - i2
          i1 = i2
          go to 10
        end if

        is = - sign ( 1, p(i) )
        p(i) = sign ( p(i), is )

      end do

      do i = 1, n

        i1 = - p(i)

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

      return
      end
      subroutine points_nearest_point_bins_2d ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, p, i_min,
     &  dist_min, compares )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS_2D finds the nearest point to a given point in 2D.
c
c  Discussion:
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = 5:
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 < P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X < P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point P, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if P lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing P, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to P.  We now know that
c       we don't need to search any cell whose points will all be further
c       from P than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       P than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, integer NBIN, the number of cells.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
c    the index of the first and last element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision P(2), the point to be tested.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to P.
c
c    Output, double precision DIST_MIN, the distance between P and PSET(*,I_MIN).
c
c    Output, integer COMPARES, the number of point-to-point comparisons.
c
      implicit none

      integer nbin
      integer ndim
      parameter ( ndim = 2 )
      integer nset

      integer bin(ndim)
      integer bin_last(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin,nbin)
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      integer compares
      integer dim
      double precision dist_min
      double precision dist_min_sq
      double precision dist_sq
      integer i
      integer i_min
      integer ic
      integer il
      integer j
      integer jc
      integer jl
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision p(ndim)
      double precision pset(ndim,nset)
      double precision r8_huge
      double precision search_radius

      compares = 0
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        dist_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        dist_min = 0.0D+00
        do dim = 1, ndim
          dist_min = dist_min + ( p(dim) - pset(dim,1) )**2
        end do
        dist_min = sqrt ( dist_min )
        compares = 1
        i_min = 1
        return
      end if
c
c  Initialize.
c
      dist_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00

      layer_width = r8_huge ( )
      do dim = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(dim) - bin_min(dim) ) / dble ( nbin ) )
      end do
c
c  Determine the bin coordinates of the point P.
c
      call r82_to_bin_even2 ( nbin, bin_min, bin_max, p, bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min,
     &  c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)
c
c  Set
c  * the current layer,
c  * the starting bin of the current layer,
c  * the current bin
c
      layer = 0
      il = ic
      jl = jc
      i = il
      j = jl

10    continue
c
c  Search all legal bins in layer LAYER.
c
20      continue
c
c  Search BIN I, J.
c
          if ( 1 .le. i .and. i .le. nbin .and.
     &         1 .le. j .and. j .le. nbin ) then

            node = bin_start(i,j)

30          continue

            if ( 0 .lt. node ) then

              dist_sq = 0.0D+00
              do dim = 1, ndim
                dist_sq = dist_sq + ( p(dim) - pset(dim,node) )**2
              end do
              compares = compares + 1

              if ( dist_sq .lt. dist_min_sq ) then
                dist_min_sq = dist_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you're done the layer.
c
          more_bins = .true.

40        continue

            if ( i .lt. ic + layer .and. j .eq. jc - layer ) then
              i = i + 1
            else if ( i .eq. ic + layer .and.
     &                j .lt. jc + layer ) then
              j = j + 1
            else if ( ic - layer .lt. i .and. j .eq. jc + layer ) then
              i = i - 1
            else if ( i .eq. ic - layer .and.
     &                jc - layer + 1 .lt. j ) then
              j = j - 1
            else
              more_bins = .false.
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin .and.
     &           1 .le. j .and. j .le. nbin ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20

60      continue
c
c  We've completed this layer.
c  Update the radius of the searched area.
c
        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do dim = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( p(dim) - c_min(dim) ) )
            search_radius = min ( search_radius,
     &        abs ( p(dim) - c_max(dim) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We can stop if:
c
c    * We've found at least one neighbor;
c    AND
c    * We've searched all legal boxes that contain the circle
c      with P at the center and the nearest N on the circumference.
c
       if ( i_min .ne. 0 ) then
         dist_min = sqrt ( dist_min_sq )
         if ( dist_min .le. search_radius ) then
           go to 70
         end if
       end if
c
c  Prepare to search the next layer.
c
        layer = layer + 1

        il = ic - layer
        jl = jc - layer

        i = il
        j = jl

      go to 10

70    continue

      return
      end
      subroutine points_nearest_point_bins_2d_2 ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ptest, i_min,
     &  dist_min, compares )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS_2D_2 finds the nearest point to a given point in 2D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINT_BINS_2D by calling
c    a subroutine to compute the next bin index.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = 5:
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 < P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X < P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we don't need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, integer NBIN, the number of cells.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
c    the index of the first and last element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(2), the coordinates of the test points.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
c
c    Output, double precision DIST_MIN, the distance between PTEST and
c    PSET(*,I_MIN).
c
c    Output, integer COMPARES, the number of point-to-point comparisons.
c
      implicit none

      integer nbin
      integer ndim
      parameter ( ndim = 2 )
      integer nset

      integer bin(ndim)
      integer bin_last(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin,nbin)
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      integer compares
      integer dim
      double precision dist_min
      double precision dist_min_sq
      double precision dist_sq
      integer i
      integer i_min
      integer ic
      integer j
      integer jc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge
      double precision search_radius

      compares = 0
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        dist_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        dist_min = 0.0D+00
        do dim = 1, ndim
          dist_min = dist_min + ( ptest(dim) - pset(dim,1) )**2
        end do
        dist_min = sqrt ( dist_min )
        compares = 1
        i_min = 1
        return
      end if

      layer_width = r8_huge ( )
      do dim = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(dim) - bin_min(dim) ) / dble ( nbin ) )
      end do
c
c  Search for each of the NTEST points.
c
      dist_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
      call r82_to_bin_even2 ( nbin, bin_min, bin_max, ptest, bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min,
     &  c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)

      layer = 0
c
c  Search all legal bins in layer LAYER.
c
10    continue

        more_bins = .false.
        call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &    more_bins )
c
c  In layer LAYER, search each BIN I, J.
c
20      continue

          if ( 1 .le. i .and. i .le. nbin .and.
     &         1 .le. j .and. j .le. nbin ) then

            node = bin_start(i,j)

30          continue
            if ( 0 .lt. node ) then

              dist_sq = 0.0D+00
              do dim = 1, ndim
                dist_sq = dist_sq + ( ptest(dim) - pset(dim,node) )**2
              end do
              compares = compares + 1

              if ( dist_sq .lt. dist_min_sq ) then
                dist_min_sq = dist_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you're done the layer.
c
40        continue

            call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &        more_bins )

            if ( .not. more_bins ) then
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin .and.
     &           1 .le. j .and. j .le. nbin ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20

60      continue
c
c  We've completed layer LAYER.
c  Update the radius of the searched area.
c
        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do dim = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( ptest(dim) - c_min(dim) ) )
            search_radius = min ( search_radius,
     &        abs ( ptest(dim) - c_max(dim) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We are done with PTEST if:
c
c    * We've found at least one neighbor;
c    AND
c    * We've searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
        if ( i_min .ne. 0 ) then
          dist_min = sqrt ( dist_min_sq )
          if ( dist_min .le. search_radius ) then
            go to 70
          end if
        end if

        layer = layer + 1

      go to 10

70    continue
c
c  We are now done with all the layers.
c
      return
      end
      subroutine points_nearest_point_bins_2d_3 ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ptest,
     &  i_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS_2D_3 finds the nearest point to a given point in 2D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
c    user to specify the number of bins in each dimension.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
c    cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = (/ 5, 4 /)
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 < P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X < P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we do not need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the coordinates of the points
c    in the set.
c
c    Input, integer NBIN(2), the number of cells in the horizontal and
c    vertical directions.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and maximum
c    bin values.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    indicates the index of the first and last element in the bin, or -1
c    if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(2), the coordinates of the test point.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer nbin(ndim)
      integer nset

      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin(1),nbin(2))
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      double precision d_min
      double precision d_min_sq
      double precision d_sq
      integer i
      integer i_min
      integer ic
      integer j
      integer jc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge
      double precision search_radius
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        d_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        d_min = 0.0
        do i = 1, ndim
          d_min = d_min + ( ptest(i) - pset(i,1) )**2
        end do
        d_min = sqrt ( d_min )
        i_min = 1
        return
      end if
c
c  The efficiency of the code will suffer if the data in the vector
c
c    ( bin_max(1:ndim) - bin_min(1:ndim) ) / dble ( nbin(1:ndim) )
c
c  varies significantly.
c
      layer_width = r8_huge ( )
      do i = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(i) - bin_min(i) ) / dble ( nbin(i) ) )
      end do

      d_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
      call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max,
     &  ptest, bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max,
     &  c_min, c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)

      layer = 0
c
c  Search all legal bins in layer LAYER.
c
10    continue

        more_bins = .false.

        call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &    more_bins )
c
c  In layer LAYER, search each BIN I, J.
c
20      continue

          if ( 1 .le. i .and. i .le. nbin(1) .and.
     &         1 .le. j .and. j .le. nbin(2) ) then

            node = bin_start(i,j)

30          continue

            if ( 0 .lt. node ) then

              d_sq = 0.0
              do i = 1, ndim
                d_sq = d_sq + ( ptest(i) - pset(i,node) )**2
              end do

              if ( d_sq .lt. d_min_sq ) then
                d_min_sq = d_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you are done the layer.
c
40        continue

            call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &        more_bins )

            if ( .not. more_bins ) then
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin(1) .and.
     &           1 .le. j .and. j .le. nbin(2) ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20
c
c  We have completed layer LAYER.
c  Update the radius of the searched area.
c
60      continue

        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do i = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( ptest(i) - c_min(i) ) )
          end do
          do i = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( ptest(i) - c_max(i) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We are done if:
c
c    * We have found at least one neighbor;
c    AND
c    * We have searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
        if ( i_min .ne. 0 ) then
          d_min = sqrt ( d_min_sq )
          if ( d_min .le. search_radius ) then
            go to 70
          end if
        end if

        layer = layer + 1

      go to 10

70    continue

      return
      end
      subroutine points_nearest_point_bins_3d_3 ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ptest,
     &  i_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS_3D_3 finds the nearest point to a given point in 3D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINTS_BINS_3D by allowing the
c    user to specify the number of bins in each dimension.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    box.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN(1) by NBIN(2) by NBIN(3)
c    regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = (/ 5, 4, 2 /)
c
c             Z LAYER 1                       Z LAYER 2
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |     | 36 | 37 | 38 | 39 | 40 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |     | 31 | 32 | 33 | 34 | 35 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |     | 26 | 27 | 28 | 29 | 30 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 < P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X < P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
c      * P1 is in the same cell as P2, P1.X = P2.X, P1.Y = P2.Y,
c        but P1.Z < P2.Z
c
c    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
c    I, J, K of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J, K cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we do not need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(3,NSET), the coordinates of the points
c    in the set.
c
c    Input, integer NBIN(3), the number of cells in the X, Y and Z directions.
c
c    Input, double precision BIN_MIN(3), BIN_MAX(3), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the first and last
c    element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(3), the coordinates of the test points.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer nbin(ndim)
      integer nset

      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      double precision d_min
      double precision d_min_sq
      double precision d_sq
      integer i
      integer i_min
      integer ic
      integer j
      integer jc
      integer k
      integer kc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge
      double precision search_radius
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        d_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        d_min = 0.0
        do i = 1, ndim
          d_min = d_min + ( ptest(i) - pset(i,1) )**2
        end do
        d_min = sqrt ( d_min )
        i_min = 1
        return
      end if
c
c  The efficiency of the code will suffer if the data in the vector
c
c    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
c
c  varies significantly.
c
      layer_width = r8_huge ( )
      do i = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(i) - bin_min(i) ) / dble ( nbin(i) ) )
      end do

      d_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
      call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max,
     &  ptest(1), bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max,
     &  c_min, c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)
      kc = bin(3)

      layer = 0
c
c  Search all legal bins in layer LAYER.
c
10    continue

        more_bins = .false.

        call index_box2_next_3d ( layer, layer, layer, ic, jc, kc,
     &    i, j, k, more_bins )
c
c  In layer LAYER, search each BIN I, J, K.
c
20      continue

          if ( 1 .le. i .and. i .le. nbin(1) .and.
     &         1 .le. j .and. j .le. nbin(2) .and.
     &         1 .le. k .and. k .le. nbin(3) ) then

            node = bin_start(i,j,k)

30          continue

            if ( 0 .lt. node ) then

              d_sq = 0.0
              do i = 1, ndim
                d_sq = d_sq + ( ptest(i) - pset(i,node) )**2
              end do

              if ( d_sq .lt. d_min_sq ) then
                d_min_sq = d_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J, K.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you are done the layer.
c
40        continue

            call index_box2_next_3d ( layer, layer, layer, ic, jc,
     &        kc, i, j, k, more_bins )

            if ( .not. more_bins ) then
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin(1) .and.
     &           1 .le. j .and. j .le. nbin(2) .and.
     &           1 .le. k .and. k .le. nbin(3) ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20
c
c  We have completed layer LAYER.
c  Update the radius of the searched area.
c
60      continue

        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do i = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( ptest(i) - c_min(i) ) )
          end do
          do i = 1, ndim
            search_radius = min ( search_radius,
     &        abs ( ptest(i) - c_max(i) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We are done with PTEST if:
c
c    * We have found at least one neighbor;
c    AND
c    * We have searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
        if ( i_min .ne. 0 ) then
          d_min = sqrt ( d_min_sq )
          if ( d_min .le. search_radius ) then
            go to 70
          end if
        end if

        layer = layer + 1

      go to 10

70    continue
c
c  We are now done with all the layers.
c
      return
      end
      subroutine points_nearest_point_del_2d ( point_num, xc, xd,
     &  nabes_first, nabes_num, nabes_dim, nabes, nnear, dnear )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_DEL_2D: nearest neighbor in a Delaunay triangulation.
c
c  Discussion:
c
c    A set of points XC is given, along with its Delaunay triangulation.
c    Now a new point XD is given, and we need to know the closest point in XC.
c
c    This algorithm starts at a random point in XC, and then repeatedly moves
c    to a neighboring point that is closer to XD.  This is guaranteed to be
c    possible because the triangulation is Delaunay.  Otherwise, it
c    would be possible to reach a vertex which was not the closest,
c    but for which all neighbors were further away.
c
c    This algorithm is able to handle the case where the point XD lies
c    outside the convex hull.
c
c    The algorithm is very simple to code.  In the most likely
c    case, it should have an expected time complexity of O(sqrt(N)).
c
c    Overhead occurs in the development of the vertex adjacency data
c    structure.  The size of this array should be roughly 6*N on average.
c    Given the list of nodes that make up triangles, the vertex adjacency
c    data can be constructed by storing every pair of nodes (I,J) and (J,I),
c    and sorting the data into dictionary order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, double precision XC(2,POINT_NUM), the set of points.
c
c    Input, double precision XD(2), a point whose nearest neighbor is
c    to be found.
c
c    Input, integer NABES_FIRST(POINT_NUM), the index in NABES of the first
c    neighbor in the list for each node.
c
c    Input, integer NABES_NUM(POINT_NUM), the number of neighbors of each node.
c
c    Input, integer NABES_DIM, the dimension of NABES.
c
c    Input, integer NABES(NABES_DIM), a list of the neighbors of all the nodes.
c    Neighbors of node 1 are listed first, and so on.
c
c    Output, integer NNEAR, the nearest node to XD.
c
c    Output, double precision DNEAR, the distance of the nearest node to XD.
c
      implicit none

      integer nabes_dim
      integer point_num

      double precision dist
      double precision dnear
      integer i
      integer i1
      integer j
      integer nabes(nabes_dim)
      integer nabes_first(point_num)
      integer nabes_num(point_num)
      integer nnear
      integer nnear_old
      double precision x
      double precision x1
      double precision xc(2,point_num)
      double precision xd(2)
      double precision y
      double precision y1

      x = xd(1)
      y = xd(2)
c
c  Select a random vertex.
c
      nnear = 1
      x1 = xc(1,nnear)
      y1 = xc(2,nnear)
      dnear = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )
c
c  From the current vertex, consider all neighboring vertices.
c  For each neighbor, compute the distance to the point.
c  If no neighbor is closer, then the current vertex is the closest.
c  Otherwise, set the current vertex to the neighbor that was closest,
c  and repeat.
c
10    continue

        nnear_old = nnear

        j = nabes_first(nnear_old)

        do i = 1, nabes_num(nnear_old)

          i1 = nabes(j)
          x1 = xc(1,i1)
          y1 = xc(2,i1)
          dist = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )

          if ( dist .lt. dnear ) then
            dnear = dist
            nnear = i1
          end if

          j = j + 1

        end do
c
c  If no neighbor was closer, we're done.
c
        if ( nnear .eq. nnear_old ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine points_nearest_point_naive_2d ( nset, pset, ptest,
     &  j_min, dist_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
c
c  Discussion:
c
c    A naive algorithm is used.  The distance to every point is calculated,
c    in order to determine the smallest.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, double precision PTEST(2), the point whose nearest neighbor
c    is sought.
c
c    Output, integer J_MIN, the index of the nearest point in PSET to P.
c
c    Output, double precision DIST_MIN, the distance between P and PSET(*,J_MIN).
c
      implicit none

      integer nset
      integer ndim
      parameter ( ndim = 2 )

      double precision d
      double precision dist_min
      integer i
      integer j
      integer j_min
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge

      dist_min = r8_huge ( )
      j_min = -1

      do j = 1, nset
        d = 0.0D+00
        do i = 1, ndim
          d = d + ( ptest(i) - pset(i,j) )**2
        end do
        if ( d .lt. dist_min ) then
          dist_min = d
          j_min = j
        end if
      end do

      dist_min = sqrt ( dist_min )

      return
      end
      subroutine points_nearest_point_naive_3d ( nset, pset, ptest,
     &  j_min, dist_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
c
c  Discussion:
c
c    A naive algorithm is used.  The distance to every point is calculated,
c    in order to determine the smallest.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(3,NSET), the points in the set.
c
c    Input, double precision PTEST(3), the point whose nearest neighbor
c    is sought.
c
c    Output, integer J_MIN, the index of the nearest point in PSET to P.
c
c    Output, double precision DIST_MIN, the distance between P and PSET(*,J_MIN).
c
      implicit none

      integer nset
      integer ndim
      parameter ( ndim = 3 )

      double precision d
      double precision dist_min
      integer i
      integer j
      integer j_min
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge

      dist_min = r8_huge ( )
      j_min = -1

      do j = 1, nset
        d = 0.0D+00
        do i = 1, ndim
          d = d + ( ptest(i) - pset(i,j) )**2
        end do
        if ( d .lt. dist_min ) then
          dist_min = d
          j_min = j
        end if
      end do

      dist_min = sqrt ( dist_min )

      return
      end
      subroutine points_nearest_point_naive_nd ( ndim, n, pset, p,
     &  i_min, dist_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
c
c  Discussion:
c
c    A naive algorithm is used.  No attempt is made to optimize the
c    calculation, so there will be N distance calculations done.
c
c    For a large dataset, it would be better to group the points into
c    clusters, so that far fewer distance calculations need to be done.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDIM, the dimension of the points.
c
c    Input, integer N, the number of points in the set.
c
c    Input, double precision PSET(NDIM,N), the coordinates of the points
c    in the set.
c
c    Input, double precision P(NDIM), the point to be tested.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to P.
c
c    Output, double precision DIST_MIN, the distance between P and PSET(*,I_MIN).
c
      implicit none

      integer n
      integer ndim

      double precision d
      integer dim
      double precision dist_min
      integer i
      integer i_min
      double precision p(ndim)
      double precision pset(ndim,n)
      double precision r8_huge

      dist_min = r8_huge ( )
      i_min = -1

      do i = 1, n

        d = 0.0D+00
        do dim = 1, ndim
          d = d + ( p(dim) - pset(ndim,i) )**2
        end do

        if ( d .lt. dist_min ) then
          dist_min = d
          i_min = i
        end if

      end do
c
c  We save a little work by waiting til the end to take the square root.
c
      dist_min = sqrt ( dist_min )

      return
      end
      subroutine points_nearest_points_bins_2d ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ntest,
     &  ptest, i_min, dist_min, compares )

c*********************************************************************72
c
cc POINTS_NEAREST_POINTS_BINS_2D finds the nearest point to given points in 2D.
c
c  Discussion:
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = 5:
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 .lt. P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X .lt. P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y .lt. P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we don't need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, integer NBIN, the number of cells.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
c    the index of the first and last element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, integer NTEST, the number of test points.
c
c    Input, double precision PTEST(2,NTEST), the test points.
c
c    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
c    to PTEST(ITEST).
c
c    Output, double precision DIST_MIN(NTEST), the distance between
c    PTEST(*,ITEST) and PSET(*,I_MIN).
c
c    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
c
      implicit none

      integer nbin
      integer ndim
      parameter ( ndim = 2 )
      integer nset
      integer ntest

      integer bin(ndim)
      integer bin_last(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin,nbin)
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      integer compares(ntest)
      integer dim
      double precision dist_min(ntest)
      double precision dist_min_sq
      double precision dist_sq
      integer i
      integer i_min(ntest)
      integer ic
      integer il
      integer itest
      integer j
      integer jc
      integer jl
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim,ntest)
      double precision r8_huge
      double precision search_radius
      double precision t

      do itest = 1, ntest
        compares(itest) = 0
      end do
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        do itest = 1, ntest
          dist_min(itest) = r8_huge ( )
          i_min(itest) = 0
        end do
        return
      end if

      if ( nset .eq. 1 ) then

        do itest = 1, ntest
          t = 0.0D+00
          do dim = 1, ndim
            t = t + ( ptest(dim,itest) - pset(dim,1) )**2
          end do
          dist_min(itest) = sqrt ( t )
        end do

        do itest = 1, ntest
          compares(itest) = 1
          i_min(itest) = 1
        end do

        return

      end if

      layer_width = r8_huge ( )
      do dim = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(dim) - bin_min(dim) ) / dble ( nbin ) )
      end do
c
c  Search for each of the NTEST points.
c
      do itest = 1, ntest

        dist_min_sq = r8_huge ( )
        i_min(itest) = 0
        search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
        call r82_to_bin_even2 ( nbin, bin_min, bin_max,
     &    ptest(1,itest), bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
        call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max,
     &    c_min, c_max )
c
c  Set the central bin of the layers.
c
        ic = bin(1)
        jc = bin(2)
c
c  Set
c  * the current layer,
c  * the starting bin of the current layer,
c  * the current bin
c
        layer = 0
        il = ic
        jl = jc
        i = il
        j = jl

10      continue
c
c  Search all legal bins in layer LAYER.
c
20        continue
c
c  Search BIN I, J.
c
            if ( 1 .le. i .and. i .le. nbin .and.
     &           1 .le. j .and. j .le. nbin ) then

              node = bin_start(i,j)

30            continue

              if ( 0 .lt. node ) then

                dist_sq = 0.0D+00
                do dim = 1, ndim
                  dist_sq = dist_sq
     &              + ( ptest(dim,itest) - pset(dim,node) )**2
                end do

                compares(itest) = compares(itest) + 1

                if ( dist_sq .lt. dist_min_sq ) then
                  dist_min_sq = dist_sq
                  i_min(itest) = node
                end if

                node = bin_next(node)

                go to 30

              end if

            end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you're done the layer.
c
            more_bins = .true.

40          continue

              if ( i .lt. ic + layer .and. j .eq. jc - layer ) then
                i = i + 1
              else if ( i .eq. ic + layer .and. j .lt. jc + layer ) then
                j = j + 1
              else if ( ic - layer .lt. i .and. j .eq. jc + layer ) then
                i = i - 1
              else if ( i .eq. ic - layer .and.
     &                  jc - layer + 1 .lt. j ) then
                j = j - 1
              else
                more_bins = .false.
                exit
              end if

              if ( 1 .le. i .and. i .le. nbin .and.
     &             1 .le. j .and. j .le. nbin ) then
                go to 50
              end if

            go to 40

50          continue

            if ( .not. more_bins ) then
              go to 60
            end if

          go to 20

60        continue
c
c  We've completed this layer.
c  Update the radius of the searched area.
c
          if ( layer .eq. 0 ) then
            t = r8_huge ( )
            do dim = 1, ndim
              t = min ( t, abs ( ptest(dim,itest) - c_min(dim) ) )
              t = min ( t, abs ( ptest(dim,itest) - c_max(dim) ) )
            end do
            search_radius = t
          else
            search_radius = search_radius + layer_width
          end if
c
c  We are done with PTEST(*,ITEST) if:
c
c    * We've found at least one neighbor;
c    AND
c    * We've searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
         if ( i_min(itest) .ne. 0 ) then
           dist_min(itest) = sqrt ( dist_min_sq )
           if ( dist_min(itest) .le. search_radius ) then
             go to 70
           end if
         end if
c
c  Prepare to search the next layer.
c
          layer = layer + 1

          il = ic - layer
          jl = jc - layer

          i = il
          j = jl

        go to 10

70      continue

      end do

      return
      end
      subroutine points_nearest_points_bins_2d_2 ( nset, pset, nbin,
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ntest,
     &  ptest, i_min, dist_min, compares )

c*********************************************************************72
c
cc POINTS_NEAREST_POINTS_BINS_2D_2 finds the nearest point to given points in 2D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINTS_BINS_2D by calling
c    a subroutine to compute the next bin index.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = 5:
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 < P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X < P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we don't need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, integer NBIN, the number of cells.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
c    the index of the first and last element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, integer NTEST, the number of test points.
c
c    Input, double precision PTEST(2,NTEST), the test points.
c
c    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
c    to PTEST(ITEST).
c
c    Output, double precision DIST_MIN(NTEST), the distance between
c    PTEST(*,ITEST) and PSET(*,I_MIN).
c
c    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
c
      implicit none

      integer nbin
      integer ndim
      parameter ( ndim = 2 )
      integer nset
      integer ntest

      integer bin(ndim)
      integer bin_last(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin,nbin)
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      integer compares(ntest)
      integer dim
      double precision dist_min(ntest)
      double precision dist_min_sq
      double precision dist_sq
      integer i
      integer i_min(ntest)
      integer ic
      integer itest
      integer j
      integer jc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim,ntest)
      double precision r8_huge
      double precision search_radius
      double precision temp

      do itest = 1, ntest
        compares(itest) = 0
      end do
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        do itest = 1, ntest
          dist_min(itest) = r8_huge ( )
          i_min(itest) = 0
        end do
        return
      end if

      if ( nset .eq. 1 ) then
        do itest = 1, ntest
          temp = 0.0D+00
          do dim = 1, ndim
            temp = temp + ( ptest(dim,itest) - pset(dim,1) )**2
          end do
          dist_min(itest) = sqrt ( temp )
        end do
        do itest = 1, ntest
          compares(itest) = 1
          i_min(itest) = 1
        end do
        return
      end if

      layer_width = r8_huge ( )
      do dim = 1, ndim
        layer_width = min ( layer_width,
     &    ( bin_max(dim) - bin_min(dim) ) / dble ( nbin ) )
      end do
c
c  Search for each of the NTEST points.
c
      do itest = 1, ntest

        dist_min_sq = r8_huge ( )
        i_min(itest) = 0
        search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
        call r82_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest),
     &    bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
        call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min,
     &    c_max )
c
c  Set the central bin of the layers.
c
        ic = bin(1)
        jc = bin(2)

        layer = 0
c
c  Search all legal bins in layer LAYER.
c
10      continue

          more_bins = .false.
          call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &      more_bins )
c
c  In layer LAYER, search each BIN I, J.
c
20        continue

            if ( 1 .le. i .and. i .le. nbin .and.
     &           1 .le. j .and. j .le. nbin ) then

              node = bin_start(i,j)

30            continue

              if ( 0 .lt. node ) then

                dist_sq = 0.0D+00
                do dim = 1, ndim
                  dist_sq = dist_sq
     &              + ( ptest(dim,itest) - pset(dim,node) )**2
                end do

                compares(itest) = compares(itest) + 1

                if ( dist_sq .lt. dist_min_sq ) then
                  dist_min_sq = dist_sq
                  i_min(itest) = node
                end if

                node = bin_next(node)
                go to 30

              end if

            end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you're done the layer.
c
40          continue

              call index_box2_next_2d ( layer, layer, ic, jc, i, j,
     &          more_bins )

              if ( .not. more_bins ) then
                go to 50
              end if

              if ( 1 .le. i .and. i .le. nbin .and.
     &             1 .le. j .and. j .le. nbin ) then
                go to 50
              end if

            go to 40

50          continue

            if ( .not. more_bins ) then
              go to 60
            end if

          go to 20

60        continue
c
c  We've completed layer LAYER.
c  Update the radius of the searched area.
c
          if ( layer .eq. 0 ) then
            temp = r8_huge ( )
            do dim = 1, ndim
              temp = min ( temp, abs ( ptest(dim,itest) - c_min(dim) ) )
              temp = min ( temp, abs ( ptest(dim,itest) - c_max(dim) ) )
            end do
            search_radius = temp
          else
            search_radius = search_radius + layer_width
          end if
c
c  We are done with PTEST(*,ITEST) if:
c
c    * We've found at least one neighbor;
c    AND
c    * We've searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
          if ( i_min(itest) .ne. 0 ) then
            dist_min(itest) = sqrt ( dist_min_sq )
            if ( dist_min(itest) .le. search_radius ) then
              go to 70
            end if
          end if

          layer = layer + 1

        go to 10

70      continue
c
c  We are now done with all the layers.
c
      end do
c
c  We are now done with all the test points.
c
      return
      end
      subroutine points_nearest_points_naive_2d ( nset, pset,
     &  ntest, ptest, i_min, dist_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
c
c  Discussion:
c
c    A naive algorithm is used.  The distance to every point is calculated,
c    in order to determine the smallest.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the points in the set.
c
c    Input, integer NTEST, the number of test points.
c
c    Input, double precision PTEST(2,NTEST), the test points.
c
c    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
c    to PTEST(ITEST).
c
c    Output, double precision DIST_MIN(NTEST), the distance between PTEST(*,ITEST)
c    and PSET(*,I_MIN).
c
      implicit none

      integer nset
      integer ndim
      parameter ( ndim = 2 )
      integer ntest

      double precision d
      integer dim
      double precision dist_min(ntest)
      integer i
      integer i_min(ntest)
      integer itest
      double precision pset(ndim,nset)
      double precision ptest(ndim,ntest)
      double precision r8_huge

      do itest = 1, ntest

        dist_min(itest) = r8_huge ( )
        i_min(itest) = -1

        do i = 1, nset

          d = 0.0D+00
          do dim = 1, ndim
            d = d + ( ptest(dim,itest) - pset(dim,i) )**2
          end do

          if ( d .lt. dist_min(itest) ) then
            dist_min(itest) = d
            i_min(itest) = i
          end if
        end do

        dist_min(itest) = sqrt ( dist_min(itest) )

      end do

      return
      end
      subroutine points_nearest_points_naive_3d ( nset, pset,
     &  ntest, ptest, i_min, dist_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
c
c  Discussion:
c
c    A naive algorithm is used.  The distance to every point is calculated,
c    in order to determine the smallest.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(3,NSET), the points in the set.
c
c    Input, integer NTEST, the number of test points.
c
c    Input, double precision PTEST(3,NTEST), the test points.
c
c    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
c    to PTEST(ITEST).
c
c    Output, double precision DIST_MIN(NTEST), the distance between PTEST(*,ITEST)
c    and PSET(*,I_MIN).
c
      implicit none

      integer nset
      integer ndim
      parameter ( ndim = 3 )
      integer ntest

      double precision d
      integer dim
      double precision dist_min(ntest)
      integer i
      integer i_min(ntest)
      integer itest
      double precision pset(ndim,nset)
      double precision ptest(ndim,ntest)
      double precision r8_huge

      do itest = 1, ntest

        dist_min(itest) = r8_huge ( )
        i_min(itest) = -1

        do i = 1, nset

          d = 0.0D+00
          do dim = 1, ndim
            d = d + ( ptest(dim,itest) - pset(dim,i) )**2
          end do

          if ( d .lt. dist_min(itest) ) then
            dist_min(itest) = d
            i_min(itest) = i
          end if
        end do

        dist_min(itest) = sqrt ( dist_min(itest) )

      end do

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      subroutine r8_to_bin_even ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R8_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
c    An initial bin takes everything less than A, and a final bin takes
c    everything greater than B.
c
c  Example:
c
c    NBIN = 7, A = 5, B = 15
c
c    C   BIN
c
c    1    1
c    3    1
c    4.9  1
c    5    2
c    6    2
c    7    3
c    8    3
c    9.5  4
c   13    6
c   14    6
c   15    6
c   15.1  7
c   99    7
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
c    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
c    If NBIN is 1 or 2, then everything is assigned to bin 1.
c
c    Input, double precision A, B, the lower and upper limits of the bin
c    interval.  While A is expected to be less than B, the code should
c    return useful results if A is actually greater than B.
c
c    Input, double precision C, a value to be placed in a bin.
c
c    Output, integer BIN, the index of the bin to which C is assigned.
c
      implicit none

      double precision a
      double precision a2
      double precision b
      double precision b2
      integer bin
      double precision c
      integer nbin
      logical switch
c
c  Take care of special cases.
c
      if ( nbin .lt. 1 ) then
        bin = 0
        return
      else if ( nbin .eq. 1 .or. nbin .eq. 2 ) then
        bin = 1
        return
      end if

      if ( b .eq. a ) then
        bin = 0
        return
      end if
c
c  If the limits are descending, then we switch them now, and
c  unswitch the results at the end.
c
      if ( a .lt. b ) then
        switch = .false.
        a2 = a
        b2 = b
      else
        switch = .true.
        a2 = b
        b2 = a
      end if
c
c  Compute the bin.
c
      if ( c .lt. a2 ) then
        bin = 1
      else if ( c .eq. a2 ) then
        bin = 2
      else if ( c .eq. b2 ) then
        bin = nbin - 1
      else if ( b2 .lt. c ) then
        bin = nbin
      else
        bin = 2 + int ( dble ( nbin - 2 ) * ( c - a2 ) / ( b2 - a2 ) )
        bin = max ( bin, 2 )
        bin = min ( bin, nbin - 1 )
      end if
c
c  Reverse the switching.
c
      if ( switch ) then
        bin = nbin + 1 - bin
      end if

      return
      end
      subroutine r8_to_bin_even2 ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R8_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
c
c  Example:
c
c    NBIN = 5, A = 5, B = 15
c
c    <-1-+-2-+-3-+-4-+-5->
c    5   7   9  11  13  15
c
c
c    C   BIN
c
c    1    1
c    3    1
c    4.9  1
c    5    1
c    6    1
c    7.1  2
c    8    2
c    9.5  3
c   12    4
c   14    5
c   15    5
c   15.1  5
c   99    5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins.
c
c    Input, double precision A, B, the lower and upper limits of the bin
c    interval.  While A is expected to be less than B, the code should
c    return useful results if A is actually greater than B.
c
c    Input, double precision C, a value to be placed in a bin.
c
c    Output, integer BIN, the index of the bin to which C is assigned.
c
      implicit none

      double precision a
      double precision a2
      double precision b
      double precision b2
      integer bin
      double precision c
      integer nbin
      logical switch
c
c  Take care of special cases.
c
      if ( nbin .lt. 1 ) then
        bin = 0
        return
      end if

      if ( nbin .eq. 1 ) then
        bin = 1
        return
      end if

      if ( b .eq. a ) then
        bin = 1
        return
      end if
c
c  If the limits are descending, then we switch them now, and
c  unswitch the results at the end.
c
      if ( a .lt. b ) then
        switch = .false.
        a2 = a
        b2 = b
      else
        switch = .true.
        a2 = b
        b2 = a
      end if
c
c  Compute the bin.
c
      if ( c .le. a2 ) then
        bin = 1
      else if ( b2 .le. c ) then
        bin = nbin
      else
        bin = 1 + int ( dble ( nbin ) * ( c - a2 ) / ( b2 - a2 ) )
        bin = max ( bin, 1 )
        bin = min ( bin, nbin )
      end if
c
c  Reverse the switching.
c
      if ( switch ) then
        bin = nbin + 1 - bin
      end if

      return
      end
      subroutine r8_to_bin_uneven ( nbin, xbin, x, bin )

c*********************************************************************72
c
cc R8_TO_BIN_UNEVEN places X in one of several unevenly spaced bins.
c
c  Discussion:
c
c    The XBIN array is assumed to be sorted.
c
c  Example:
c
c    NBIN = 5
c    XBIN(1:4) = (/ 0.0, 2.0, 8.0, 9.0 /)
c
c    so bins are
c
c    1  ( -Inf,   0 )
c    2  (    0,   2 )
c    3  (    2,   8 )
c    4  (    8,   9 )
c    5  (    9, Inf )
c
c    X   BIN
c
c   -7    1
c   -3    1
c    0    1
c    0.1  2
c    1    2
c    3    3
c    8    3
c    9.5  5
c   13    5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins.
c
c    Input, double precision XBIN(NBIN-1), the dividing values for the bins.
c
c    Input, double precision X, a value to be placed in a bin.
c
c    Output, integer BIN, the index of the bin to which X is assigned.
c
      implicit none

      integer nbin

      integer bin
      double precision x
      double precision xbin(nbin-1)

      bin = 1

10    continue

      if ( bin .lt. nbin ) then

        if ( x .le. xbin(bin) ) then
          go to 20
        end if

        bin = bin + 1

        go to 10

      end if

20    continue

      return
      end
      function r8_uniform ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM returns a scaled pseudorandom R8.
c
c  Discussion:
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
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
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number.
c
      r8_uniform = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number.
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r82_to_bin_even ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R82_TO_BIN_EVEN determines the appropriate "bin" for an R82 value.
c
c  Discussion:
c
c    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN-2
c    equal subintervals or bins.  Boundary bins take care of extreme values.
c
c  Example:
c
c    NBIN = 7, A(1) = 5,  A(2) = 0,
c              B(1) = 15, B(2) = 20.
c
c      C      BIN
c   ------  ------
c    8 -2    3  1
c    0  1    1  2
c    6  9    2  4
c   10 11    4  4
c   14 23    6  7
c   25 13    7  5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins in each dimension.
c    NBIN is normally at least 3.  If NBIN is 1 or 2, then everything
c    is assigned to bin 1.
c
c    Input, double precision A(2), B(2), the lower and upper limits of
c    the bin interval.  While A(I) is expected to be less than B(I), the
c    code should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(2), a value to be placed in a bin.
c
c    Output, integer BIN(2), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call r8_to_bin_even ( nbin, a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r82_to_bin_even2 ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R82_TO_BIN_EVEN2 determines the appropriate "bin" for an R82 value.
c
c  Discussion:
c
c    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
c    equal subintervals or bins.  Boundary bins take care of extreme values.
c
c  Example:
c
c    NBIN = 5, A(1) = 5,  A(2) = 0,
c              B(1) = 15, B(2) = 20.
c
c   20 +    +    +    +    +    +
c        15 | 25 | 35 | 45 | 55
c   16 +----+----+----+----+----+
c        14 | 24 | 34 | 44 | 54
c   12 +----+----+----+----+----+
c        13 | 23 | 33 | 43 | 53
c    8 +----+----+----+----+----+
c        12 | 22 | 32 | 42 | 52
c    4 +----+----+----+----+----+
c        11 | 21 | 31 | 41 | 51
c    0 +    +    +    +    +    +
c      5    7    9   11   13   15
c
c      C      BIN
c   ------  ------
c    8 -2    2  1
c    0  1    1  1
c    6  9    1  3
c   10 11    3  3
c   14 23    5  5
c   25 13    5  4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision A(2), B(2), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(2), a value to be placed in a bin.
c
c    Output, integer BIN(2), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r82_to_bin_even3 ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R82_TO_BIN_EVEN3 determines the appropriate "bin" for an R82 value.
c
c  Discussion:
c
c    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
c    or bins.
c
c  Example:
c
c    NBIN = (/ 4, 5 /),
c
c      A(1) = 1,  A(2) = 0,
c      B(1) = 17, B(2) = 20.
c
c   20 +    +    +    +    +
c        15 | 25 | 35 | 45
c   16 +----+----+----+----+
c        14 | 24 | 34 | 44
c   12 +----+----+----+----+
c        13 | 23 | 33 | 43
c    8 +----+----+----+----+
c        12 | 22 | 32 | 42
c    4 +----+----+----+----+
c        11 | 21 | 31 | 41
c    0 +    +    +    +    +
c      1    5    9   13   17
c
c      C      BIN
c   ------  ------
c    8 -2    2  1
c    0  1    1  1
c    6  9    2  3
c   10 11    3  3
c   14 23    4  5
c   25 13    4  4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN(2), the number of bins in each dimension.
c
c    Input, double precision A(2), B(2), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(2), a value to be placed in a bin.
c
c    Output, integer BIN(2), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r82_uniform ( rlo, rhi, seed, r )

c*********************************************************************72
c
cc R82_UNIFORM returns a random R82 value in a given range.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RLO(2), RHI(2), the minimum and maximum values.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R(2), the randomly chosen value.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer i
      double precision r(ndim)
      double precision r8_uniform
      double precision rhi(ndim)
      double precision rlo(ndim)
      integer seed

      do i = 1, ndim
        r(i) = r8_uniform ( rlo(i), rhi(i), seed )
      end do

      return
      end
      subroutine r82vec_bin_even ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BIN_EVEN bins an R82VEC into evenly spaced bins.
c
c  Discussion:
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN 1D bins in both X and Y directions, making a total
c    of NBIN**2 2D bins.  Each set of 1D bins begins and ends with an
c    "open-ended" bin that catches extreme values.
c
c    The 2D bins are indexed by the X and Y bins that construct them,
c    and ordered lexicographically by these indices:
c
c      1,4 | 2,4 | 3,4 | 4,4
c      ----+-----+-----+----
c      1,3 | 2,3 | 3,3 | 4,3
c      ----+-----+-----+----
c      1,2 | 2,2 | 3,2 | 4,2
c      ----+-----+-----+----
c      1,1 | 2,1 | 3,1 | 4,1
c
c    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(2,N), the R82 data to be binned.
c
c    Input, integer NBIN, the (square root of) the number of bins.
c    NBIN must be at least 3.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the bin limits.
c
c    Output, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the index
c    of the first and last elements of A that went into each bin, or
c    -1 if there were no entries in the bin.
c
c    Output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do j = 1, nbin
        do i = 1, nbin
          bin_last(i,j) = -1
        end do
      end do
      do j = 1, nbin
        do i = 1, nbin
          bin_start(i,j) = -1
        end do
      end do
      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r82_to_bin_even ( nbin, bin_min, bin_max, a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)

        if ( bin_start(i1,i2) .eq. -1 ) then
          bin_start(i1,i2) = j
        else
          k = bin_last(i1,i2)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2) = j

      end do

      return
      end
      subroutine r82vec_bin_even2 ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BIN_EVEN2 bins an R82VEC into evenly spaced bins.
c
c  Discussion:
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN 1D bins in both X and Y directions, making a total
c    of NBIN**2 2D bins.  Each set of 1D bins begins and ends at user
c    specified mininum and maximum values.
c
c    The 2D bins are indexed by the X and Y bins that construct them,
c    and ordered lexicographically by these indices:
c
c      1,4 | 2,4 | 3,4 | 4,4
c      ----+-----+-----+----
c      1,3 | 2,3 | 3,3 | 4,3
c      ----+-----+-----+----
c      1,2 | 2,2 | 3,2 | 4,2
c      ----+-----+-----+----
c      1,1 | 2,1 | 3,1 | 4,1
c
c    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(2,N), the R82 data to be binned.
c
c    Input, integer NBIN, the (square root of) the number of bins.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the bin limits.
c
c    Output, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the
c    index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do j = 1, nbin
        do i = 1, nbin
          bin_last(i,j) = -1
        end do
      end do
      do j = 1, nbin
        do i = 1, nbin
          bin_start(i,j) = -1
        end do
      end do
      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r82_to_bin_even2 ( nbin, bin_min, bin_max, a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)

        if ( bin_start(i1,i2) .eq. -1 ) then
          bin_start(i1,i2) = j
        else
          k = bin_last(i1,i2)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2) = j

      end do

      return
      end
      subroutine r82vec_bin_even3 ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BIN_EVEN3 bins an R82VEC into evenly spaced bins.
c
c  Discussion:
c
c    A different number of bins may be used in each dimension.
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, making a
c    total of NBIN(1) * NBIN(2) 2D bins.  Each set of 1D bins begins and
c    ends at user specified mininum and maximum values.
c
c    The 2D bins are indexed by the X and Y bins that construct them,
c    and ordered lexicographically by these indices:
c
c      1,4 | 2,4 | 3,4 | 4,4 | 5,4
c      ----+-----+-----+-----+-----
c      1,3 | 2,3 | 3,3 | 4,3 | 5,3
c      ----+-----+-----+-----+-----
c      1,2 | 2,2 | 3,2 | 4,2 | 5,2
c      ----+-----+-----+-----+-----
c      1,1 | 2,1 | 3,1 | 4,1 | 5,1
c
c    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1),
c    ..., (5,4).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(2,N), the D2 data to be binned.
c
c    Input, integer NBIN(2), the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the bin limits.
c
c    Output, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    the index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), contains the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do j = 1, nbin(2)
        do i = 1, nbin(1)
          bin_last(i,j) = -1
          bin_start(i,j) = -1
        end do
      end do

      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max,
     &    a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)

        if ( bin_start(i1,i2) .eq. -1 ) then
          bin_start(i1,i2) = j
        else
          k = bin_last(i1,i2)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2) = j

      end do

      return
      end
      subroutine r82vec_binned_reorder ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BINNED_REORDER reorders a binned R82 data vector.
c
c  Discussion:
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the R82 data to be sorted.
c
c    Input, integer NBIN, the (square root of the) number of bins.
c
c    Input/output, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the
c    index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Input/output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin)
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin

        do i2 = 1, nbin

          j = bin_start(i1,i2)

          if ( 0 .lt. j ) then
            bin_start(i1,i2) = k + 1
          end if

10        continue

          if ( 0 .lt. j ) then
            k = k + 1
            bin_last(i1,i2) = k
            a2(1:ndim,k) = a(1:ndim,j)
            j = bin_next(j)
            go to 10
          end if

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin
        do i2 = 1, nbin

          k = bin_last(i1,i2)

          if ( 0 .lt. k ) then
            bin_next(k) = 0
          end if

        end do
      end do

      return
      end
      subroutine r82vec_binned_reorder2 ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BINNED_REORDER2 reorders a binned R82VEC.
c
c  Discussion:
c
c    This routine allows there to be a different number of bins in
c    each dimension.
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the D2 data to be sorted.
c
c    Input, integer NBIN(2), the number of bins in each direction.
c
c    Input/output, integer BIN_START(NBIN(1),NBIN(2)),
c    BIN_LAST(NBIN(1),NBIN(2)), contains the index of the first and last
c    element of A that went into each bin, or -1 if there are no entries
c    in the bin.
c
c    Input/output, integer BIN_NEXT(N), contains the index of the next
c    element of A that follows this element in the same bin.  A value of
c    0 means this is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin(1),nbin(2))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2))
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          j = bin_start(i1,i2)

          if ( 0 .lt. j ) then
            bin_start(i1,i2) = k + 1
          end if

10        continue

          if ( 0 .lt. j ) then
            k = k + 1
            bin_last(i1,i2) = k
            do i = 1, ndim
              a2(i,k) = a(i,j)
            end do
            j = bin_next(j)
            go to 10
          end if

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin(1)
        do i2 = 1, nbin(2)

          k = bin_last(i1,i2)

          if ( 0 .lt. k ) then
            bin_next(k) = 0
          end if

        end do
      end do

      return
      end
      subroutine r82vec_binned_sort_a ( n, a, nbin, bin_start,
     &  bin_last )

c*********************************************************************72
c
cc R82VEC_BINNED_SORT_A sorts each bin of a binned R82VEC.
c
c  Discussion:
c
c    Presumably, the data vector was first binned by R82VEC_BIN_EVEN,
c    then reordered by R82VEC_BINNED_REORDER.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted R82VEC.
c
c    What is true is that if I < J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the R82 data to be sorted.
c
c    Input, integer NBIN, the (square root of the) number of bins.
c
c    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the index
c    of the first and last element of A that went into each bin, or -1
c    if there are no entries in this bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      integer bin_last(nbin,nbin)
      integer bin_start(nbin,nbin)
      integer i1
      integer i2
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin

        do i2 = 1, nbin

          j1 = bin_start(i1,i2)

          if ( 0 .lt. j1 ) then

            j2 = bin_last(i1,i2)

            n1 = j2 + 1 - j1

            if ( 1 .lt. n1 ) then
              call r82vec_sort_quick_a ( n1, a(1,j1) )
            end if

          end if

        end do

      end do

      return
      end
      subroutine r82vec_binned_sort_a2 ( n, a, nbin, bin_start,
     &  bin_last )

c*********************************************************************72
c
cc R82VEC_BINNED_SORT_A2 sorts each bin of a binned R82VEC.
c
c  Discussion:
c
c    This routine allows a different number of bins in each dimension.
c
c    Presumably, the data vector was first binned by R82VEC_BIN_EVEN3,
c    then reordered by R82VEC_BINNED_REORDER2.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted D2 vector.
c
c    What is true is that if I < J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the R2 data to be sorted.
c
c    Input, integer NBIN(2), the number of bins in each dimension.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    the index of the first and last element of A that went into each bin, or -1
c    if there are no entries in this bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin_last(nbin(1),nbin(2))
      integer bin_start(nbin(1),nbin(2))
      integer i1
      integer i2
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          j1 = bin_start(i1,i2)

          if ( 0 .lt. j1 ) then

            j2 = bin_last(i1,i2)

            n1 = j2 + 1 - j1

            if ( 1 .lt. n1 ) then
              call r82vec_sort_quick_a ( n1, a(1,j1) )
            end if

          end if

        end do

      end do

      return
      end
      subroutine r82vec_part_quick_a ( n, a, l, r )

c*********************************************************************72
c
cc R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.
c
c  Discussion:
c
c    The routine reorders the entries of A.  Using A(1:2,1) as a
c    key, all entries of A that are less than or equal to the key will
c    precede the key, which precedes all entries that are greater than the key.
c
c  Example:
c
c    Input:
c
c      N = 8
c
c      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
c
c    Output:
c
c      L = 2, R = 4
c
c      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
c             -----------          ----------------------------------
c             LEFT          KEY    RIGHT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of A.
c
c    Input/output, double precision A(2,N).  On input, the array to be checked.
c    On output, A has been reordered as described above.
c
c    Output, integer L, R, the indices of A that define the three segments.
c    Let KEY = the input value of A(1:2,1).  Then
c    I <= L                 A(1:2,I) < KEY;
c         L < I < R         A(1:2,I) = KEY;
c                 R <= I    A(1:2,I) > KEY.
c
      implicit none

      integer n
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      logical r8vec_eq
      logical r8vec_gt
      logical r8vec_lt
      integer i
      integer j
      double precision key(ndim)
      integer l
      integer m
      integer r

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_PART_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        l = 0
        r = 2
        return
      end if

      do i = 1, ndim
        key(i) = a(i,1)
      end do
      m = 1
c
c  The elements of unknown size have indices between L+1 and R-1.
c
      l = 1
      r = n + 1

      do i = 2, n

        if ( r8vec_gt ( ndim, a(1,l+1), key ) ) then
          r = r - 1
          call r8vec_swap ( ndim, a(1,r), a(1,l+1) )
        else if ( r8vec_eq ( ndim, a(1,l+1), key ) ) then
          m = m + 1
          call r8vec_swap ( ndim, a(1,m), a(1,l+1) )
          l = l + 1
        else if ( r8vec_lt ( ndim, a(1,l+1), key ) ) then
          l = l + 1
        end if

      end do
c
c  Now shift small elements to the left, and KEY elements to center.
c
      do j = 1, l - m
        do i = 1, ndim
          a(i,j) = a(i,j+m)
        end do
      end do

      l = l - m

      do j = l + 1, l + m
        do i = 1, ndim
          a(i,j) = key(i)
        end do
      end do

      return
      end
      subroutine r82vec_permute ( n, a, p )

c*****************************************************************************80
c
cc R82VEC_PERMUTE permutes an R82VEC in place.
c
c  Discussion:
c
c    This routine permutes an array of real "objects", but the same
c    logic can be used to permute an array of objects of any arithmetic
c    type, or an array of objects of any complexity.  The only temporary
c    storage required is enough to store a single object.  The number
c    of data movements made is N + the number of cycles of order 2 or more,
c    which is never more than N + N/2.
c
c  Example:
c
c    Input:
c
c      N = 5
c      P = (   2,    4,    5,    1,    3 )
c      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
c          (11.0, 22.0, 33.0, 44.0, 55.0 )
c
c    Output:
c
c      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
c             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input/output, double precision A(2,N), the array to be permuted.
c
c    Input, integer P(N), the permutation.  P(I) = J means
c    that the I-th element of the output array should be the J-th
c    element of the input array.  P must be a legal permutation
c    of the integers from 1 to N, otherwise the algorithm will
c    fail catastrophically.
c
      implicit none

      integer n

      double precision a(2,n)
      double precision a_temp(2)
      integer i
      integer iget
      integer iput
      integer istart
      integer p(n)
c
c  Search for the next element of the permutation that has not been used.
c
      do istart = 1, n

        if ( p(istart) .lt. 0 ) then

        else if ( p(istart) .eq. istart ) then

          p(istart) = - p(istart)

        else

          do i = 1, 2
            a_temp(i) = a(i,istart)
          end do

          iget = istart
c
c  Copy the new value into the vacated entry.
c
10        continue

            iput = iget
            iget = p(iget)

            p(iput) = - p(iput)

            if ( iget .lt. 1 .or. n .lt. iget ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
              stop
            end if

            if ( iget .eq. istart ) then
              do i = 1, 2
                a(i,iput) = a_temp(i)
              end do
              go to 20
            end if

            do i = 1, 2
              a(i,iput) = a(i,iget)
            end do

          go to 10

20        continue

        end if

      end do
c
c  Restore the signs of the entries.
c
      do i = 1, n
        p(i) = - p(i)
      end do

      return
      end
      subroutine r82vec_print ( n, a, title )

c*********************************************************************72
c
cc R82VEC_PRINT prints an R82VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(2,N), the R82 vector to be printed.
c
c    Input, character ( len = * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      integer i
      integer j
      character * ( * ) title

      if ( title .ne. ' ' ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
      end if

      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(i6,(5g14.6))' ) j, ( a(i,j), i = 1, 2 )
      end do

      return
      end
      subroutine r82vec_sort_quick_a ( n, a )

c*********************************************************************72
c
cc R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, real ( kind = 8 ) A(2,N).
c    On input, the array to be sorted.
c    On output, the array has been sorted.
c
      implicit none

      integer maxlevel
      parameter ( maxlevel = 25 )
      integer n
      integer ndim
      parameter ( ndim = 2 )

      double precision a(ndim,n)
      integer base
      integer l_segment
      integer level
      integer n_segment
      integer rsave(maxlevel)
      integer r_segment

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        return
      end if

      level = 1
      rsave(level) = n + 1
      base = 1
      n_segment = n

10    continue
c
c  Partition the segment.
c
        call r82vec_part_quick_a ( n_segment, a(1,base), l_segment,
     &    r_segment )
c
c  If the left segment has more than one element, we need to partition it.
c
        if ( 1 .lt. l_segment ) then

          if ( maxlevel .lt. level ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal errorc'
            write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ',
     &        maxlevel
            stop
          end if

          level = level + 1
          n_segment = l_segment
          rsave(level) = r_segment + base - 1
c
c  The left segment and the middle segment are sorted.
c  Must the right segment be partitioned?
c
        else if ( r_segment .lt. n_segment ) then

          n_segment = n_segment + 1 - r_segment
          base = base + r_segment - 1
c
c  Otherwise, we back up a level if there is an earlier one.
c
        else

20        continue

            if ( level .le. 1 ) then
              go to 40
            end if

            base = rsave(level)
            n_segment = rsave(level-1) - rsave(level)
            level = level - 1

            if ( 0 .lt. n_segment ) then
              go to 30
            end if

          go to 20

30        continue

        end if

      go to 10

40    continue

      return
      end
      subroutine r82vec_uniform ( n, alo, ahi, seed, a )

c*********************************************************************72
c
cc R82VEC_UNIFORM returns a random R82VEC in a given range.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision ALO(2), AHI(2), the minimum and maximum
c    values allowed for A(1,1:N) and A(2,1:N).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision A(2,N), the vector of randomly chosen values.
c
      implicit none

      integer n

      double precision a(2,n)
      double precision ahi(2)
      double precision alo(2)
      integer i
      integer j
      double precision r8_uniform
      integer seed

      do i = 1, 2
        do j = 1, n
          a(i,j) = r8_uniform ( alo(i), ahi(i), seed )
        end do
      end do

      return
      end
      subroutine r83_to_bin_even2 ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R83_TO_BIN_EVEN2 determines the appropriate "bin" for an R83 value.
c
c  Discussion:
c
c    The intervals [A(I),B(I)] are each divided into NBIN
c    equal subintervals or bins.  Boundary bins take care of extreme values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision A(3), B(3), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(3), a value to be placed in a bin.
c
c    Output, integer BIN(3), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r83_to_bin_even3 ( nbin, a, b, c, bin )

c*****************************************************************************80
c
cc R83_TO_BIN_EVEN3 determines the appropriate "bin" for an R83 value.
c
c  Discussion:
c
c    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
c    or bins.
c
c  Example:
c
c    NBIN = (/ 4, 5, 2 /),
c
c      A(1) = 1,  A(2) = 0,  A(3) = 8
c      B(1) = 17, B(2) = 20, B(3) = 10
c
c
c            8 < Z < 9                    9 < Z < 10
c
c   20 +     +     +     +     +     20 +     +     +     +     +
c        151 | 251 | 351 | 451            152 | 252 | 352 | 452
c   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
c        141 | 241 | 341 | 441            142 | 242 | 342 | 442
c   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
c        131 | 231 | 331 | 431            132 | 232 | 332 | 432
c    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
c        121 | 221 | 321 | 421            122 | 222 | 322 | 422
c    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
c        111 | 211 | 311 | 411            112 | 212 | 312 | 412
c    0 +     +     +     +     +      0 +     +     +     +     +
c      1     5     9    13    17        1     5     9    13    17
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input, double precision A(3), B(3), the lower and upper limits of the
c    bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(3), a value to be placed in a bin.
c
c    Output, integer BIN(3), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r83vec_bin_even2 ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BIN_EVEN2 bins an R83VEC into evenly spaced bins.
c
c  Discussion:
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN 1D bins in each coordinate, making a total
c    of NBIN**NDIM bins.  Each set of 1D bins begins and ends at user
c   specified mininum and maximum values.
c
c    The bins are indexed by the 1D bins that construct them,
c    and ordered lexicographically by these indices:
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(3,N), the data to be binned.
c
c    Input, integer NBIN, the (cube root of) the number of bins.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(3), BIN_MAX(3), the bin limits.
c
c    Output, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
c    the index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin,nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin,nbin)
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do k = 1, nbin
        do j = 1, n
          do i = 1, n
            bin_last(i,j,k) = -1
          end do
        end do
      end do

      do k = 1, nbin
        do j = 1, n
          do i = 1, n
            bin_start(i,j,k) = -1
          end do
        end do
      end do

      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r83_to_bin_even2 ( nbin, bin_min, bin_max, a(1,j),
     &    bin )

        i1 = bin(1)
        i2 = bin(2)
        i3 = bin(3)

        if ( bin_start(i1,i2,i3) .eq. -1 ) then
          bin_start(i1,i2,i3) = j
        else
          k = bin_last(i1,i2,i3)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2,i3) = j

      end do

      return
      end
      subroutine r83vec_bin_even3 ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BIN_EVEN3 bins an R83VEC into evenly spaced bins.
c
c  Discussion:
c
c    A different number of bins may be used in each dimension.
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, and
c    NBIN(3) for Z, making a total of NBIN(1) * NBIN(2) * NBIN(3) 3D bins.
c    Each set of 1D bins begins and ends at user specified mininum and
c    maximum values.
c
c    The 3D bins are indexed by the X, Y and Z bins that construct them,
c    and ordered lexicographically by these indices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(3,N), the R3 data to be binned.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(3), BIN_MAX(3), the bin limits.
c
c    Output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
c    index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), contains the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do k = 1, nbin(3)
        do j = 1, nbin(2)
          do i = 1, nbin(1)
            bin_last(i,j,k) = -1
            bin_start(i,j,k) = -1
          end do
        end do
      end do

      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max,
     &    a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)
        i3 = bin(3)

        if ( bin_start(i1,i2,i3) .eq. -1 ) then
          bin_start(i1,i2,i3) = j
        else
          k = bin_last(i1,i2,i3)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2,i3) = j

      end do

      return
      end
      subroutine r83vec_binned_reorder ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BINNED_REORDER reorders a binned R83 data vector.
c
c  Discussion:
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN, the (cube root of the) number of bins.
c
c    Input/output, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
c    the index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Input/output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin,nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin,nbin)
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin

        do i2 = 1, nbin

          do i3 = 1, nbin

            j = bin_start(i1,i2,i3)

            if ( 0 .lt. j ) then
              bin_start(i1,i2,i3) = k + 1
            end if

10          continue

            if ( 0 .lt. j ) then
              k = k + 1
              bin_last(i1,i2,i3) = k
              do i = 1, ndim
                a2(i,k) = a(i,j)
              end do
              j = bin_next(j)
              go to 10
            end if

          end do

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin
        do i2 = 1, nbin
          do i3 = 1, nbin
            k = bin_last(i1,i2,i3)

            if ( 0 .lt. k ) then
              bin_next(k) = 0
            end if

          end do
        end do
      end do

      return
      end
      subroutine r83vec_binned_reorder2 ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BINNED_REORDER2 reorders a binned R83VEC.
c
c  Discussion:
c
c    This routine allows there to be a different number of bins in
c    each dimension.
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input/output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
c    index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Input/output, integer BIN_NEXT(N), contains the index of the next
c    element of A that follows this element in the same bin.  A value of 0
c    means this is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          do i3 = 1, nbin(3)

            j = bin_start(i1,i2,i3)

            if ( 0 .lt. j ) then
              bin_start(i1,i2,i3) = k + 1
            end if

10          continue

            if ( 0 .lt. j ) then
              k = k + 1
              bin_last(i1,i2,i3) = k
              do i = 1, ndim
                a2(i,k) = a(i,j)
              end do
              j = bin_next(j)
              go to 10
            end if

          end do

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin(1)
        do i2 = 1, nbin(2)
          do i3 = 1, nbin(3)

            k = bin_last(i1,i2,i3)

            if ( 0 .lt. k ) then
              bin_next(k) = 0
            end if

          end do
        end do
      end do

      return
      end
      subroutine r83vec_binned_sort_a ( n, a, nbin, bin_start,
     &  bin_last )

c*********************************************************************72
c
cc R83VEC_BINNED_SORT_A sorts each bin of a binned R83VEC.
c
c  Discussion:
c
c    Presumably, the data vector was first binned by R83VEC_BIN_EVEN,
c    then reordered by R83VEC_BINNED_REORDER.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted R83 vector.
c
c    What is true is that if I < J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN, the (cube root of the) number of bins.
c
c    Input, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
c    the index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in this bin.
c
      implicit none

      integer n
      integer nbin
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      integer bin_last(nbin,nbin,nbin)
      integer bin_start(nbin,nbin,nbin)
      integer i1
      integer i2
      integer i3
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin

        do i2 = 1, nbin

          do i3 = 1, nbin

            j1 = bin_start(i1,i2,i3)

            if ( 0 .lt. j1 ) then

              j2 = bin_last(i1,i2,i3)

              n1 = j2 + 1 - j1

              if ( 1 .lt. n1 ) then
                call r83vec_sort_quick_a ( n1, a(1,j1) )
              end if

            end if

          end do

        end do

      end do

      return
      end
      subroutine r83vec_binned_sort_a2 ( n, a, nbin, bin_start,
     &  bin_last )

c*********************************************************************72
c
cc R83VEC_BINNED_SORT_A2 sorts each bin of a binned R83VEC.
c
c  Discussion:
c
c    This routine allows a different number of bins in each dimension.
c
c    Presumably, the data vector was first binned by R83VEC_BIN_EVEN3,
c    then reordered by R83VEC_BINNED_REORDER2.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted D3 vector.
c
c    What is true is that if I < J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains
c    the index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in this bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer i1
      integer i2
      integer i3
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          do i3 = 1, nbin(3)

            j1 = bin_start(i1,i2,i3)

            if ( 0 .lt. j1 ) then

              j2 = bin_last(i1,i2,i3)

              n1 = j2 + 1 - j1

              if ( 1 .lt. n1 ) then
                call r82vec_sort_quick_a ( n1, a(1,j1) )
              end if

            end if

          end do

        end do

      end do

      return
      end
      subroutine r83vec_part_quick_a ( n, a, l, r )

c*********************************************************************72
c
cc R83VEC_PART_QUICK_A reorders an R83VEC as part of a quick sort.
c
c  Discussion:
c
c    The routine reorders the entries of A.  Using A(1:3,1) as a
c    key, all entries of A that are less than or equal to the key will
c    precede the key, which precedes all entries that are greater than the key.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of A.
c
c    Input/output, double precision A(3,N).  On input, the array to be checked.
c    On output, A has been reordered as described above.
c
c    Output, integer L, R, the indices of A that define the three segments.
c    Let KEY = the input value of A(1:3,1).  Then
c    I <= L                 A(1:3,I) < KEY;
c         L < I < R         A(1:3,I) = KEY;
c                 R <= I    A(1:3,I) > KEY.
c
      implicit none

      integer n
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      integer i
      integer j
      double precision key(ndim)
      integer l
      integer m
      integer r
      logical r8vec_eq
      logical r8vec_gt
      logical r8vec_lt

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83VEC_PART_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        l = 0
        r = 2
        return
      end if

      do i = 1, ndim
        key(i) = a(i,1)
      end do
      m = 1
c
c  The elements of unknown size have indices between L+1 and R-1.
c
      l = 1
      r = n + 1

      do i = 2, n

        if ( r8vec_gt ( ndim, a(1,l+1), key ) ) then
          r = r - 1
          call r8vec_swap ( ndim, a(1,r), a(1,l+1) )
        else if ( r8vec_eq ( ndim, a(1,l+1), key ) ) then
          m = m + 1
          call r8vec_swap ( ndim, a(1,m), a(1:ndim,l+1) )
          l = l + 1
        else if ( r8vec_lt ( ndim, a(1,l+1), key ) ) then
          l = l + 1
        end if

      end do
c
c  Now shift small elements to the left, and KEY elements to center.
c
      do j = 1, l - m
        do i = 1, ndim
          a(i,j) = a(i,j+m)
        end do
      end do

      l = l - m

      do j = l + 1, l + m
        do i = 1, ndim
          a(i,j) = key(i)
        end do
      end do

      return
      end
      subroutine r83vec_print ( n, a, title )

c*****************************************************************************80
c
cc R83VEC_PRINT prints an R83VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(3,N), the R83 vector to be printed.
c
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      integer i
      integer j
      character * ( * ) title

      if ( title .ne. ' ' ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
      end if

      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(i6,(5g14.6))' ) j, ( a(i,j), i = 1, ndim )
      end do

      return
      end
      subroutine r83vec_sort_quick_a ( n, a )

c*********************************************************************72
c
cc R83VEC_SORT_QUICK_A ascending sorts an R83VEC using quick sort.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, double precision A(3,N).
c    On input, the array to be sorted.
c    On output, the array has been sorted.
c
      implicit none

      integer maxlevel
      parameter ( maxlevel = 25 )
      integer n
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      integer base
      integer l_segment
      integer level
      integer n_segment
      integer rsave(maxlevel)
      integer r_segment

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      end if

      if ( n .eq. 1 ) then
        return
      end if

      level = 1
      rsave(level) = n + 1
      base = 1
      n_segment = n

10    continue
c
c  Partition the segment.
c
        call r83vec_part_quick_a ( n_segment, a(1,base), l_segment,
     &    r_segment )
c
c  If the left segment has more than one element, we need to partition it.
c
        if ( 1 .lt. l_segment ) then

          if ( maxlevel < level ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R83VEC_SORT_QUICK_A - Fatal error!'
            write ( *, '(a,i8)' )
     &        '  Exceeding recursion maximum of ', maxlevel
            stop
          end if

          level = level + 1
          n_segment = l_segment
          rsave(level) = r_segment + base - 1
c
c  The left segment and the middle segment are sorted.
c  Must the right segment be partitioned?
c
        else if ( r_segment .lt. n_segment ) then

          n_segment = n_segment + 1 - r_segment
          base = base + r_segment - 1
c
c  Otherwise, we back up a level if there is an earlier one.
c
        else

20        continue

            if ( level .le. 1 ) then
              go to 40
            end if

            base = rsave(level)
            n_segment = rsave(level-1) - rsave(level)
            level = level - 1

            if ( 0 .lt. n_segment ) then
              go to 30
            end if

          go to 20

30        continue

        end if

      go to 10

40    continue

      return
      end
      subroutine r83vec_uniform ( n, alo, ahi, seed, a )

c*********************************************************************72
c
cc R83VEC_UNIFORM returns a random R83VEC in a given range.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision ALO(3), AHI(3), the minimum and maximum values.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision A(3,N), the vector of randomly chosen values.
c
      implicit none

      integer n
      integer ndim
      parameter ( ndim = 3 )

      double precision a(ndim,n)
      double precision ahi(ndim)
      double precision alo(ndim)
      integer i
      integer j
      double precision r8_uniform
      integer seed

      do i = 1, ndim
        do j = 1, n
          a(i,j) = r8_uniform ( alo(i), ahi(i), seed )
        end do
      end do

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title to be printed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine r8mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, character*(*) TITLE, an optional title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character*(*) title

      call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi,
     &  jhi, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT transposed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      integer s_len_trim
      character * ( * ) title
      integer title_len

      title_len = s_len_trim ( title )

      if ( 0 .lt. title_len ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_len)
      end if

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(i8,6x)') i
        end do

        write ( *, '(''       Row'',5a14)' ) ctemp(1:inc)
        write ( *, '(a)' ) '       Col'

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

          do i2 = 1, inc
            i = i2lo - 1 + i2
            write ( ctemp(i2), '(g14.6)' ) a(i,j)
          end do

          write ( *, '(2x,i8,5a14)' ) j, ( ctemp(i), i = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine r8vec_bin ( n, a, nbin, bin_min, bin_max, bin,
     &  bin_limit )

c*********************************************************************72
c
cc R8VEC_BIN bins an R8VEC, returning the population of each bin.
c
c  Discussion:
c
c    The user specifies minimum and maximum bin values, BIN_MIN and
c    BIN_MAX, and the number of bins, NBIN.  This determines a
c    "bin width":
c
c      H = ( BIN_MAX - BIN_MIN ) / NBIN
c
c    so that bin I will count all entries X(J) such that
c
c      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
c
c    The array does NOT have to be sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of A.
c
c    Input, double precision A(N), an (unsorted) array to be binned.
c
c    Input, integer NBIN, the number of bins.  Two extra bins, #0 and
c    #NBIN+1, count extreme values.
c
c    Input, double precision BIN_MIN, BIN_MAX, define the range and size
c    of the bins.  BIN_MIN and BIN_MAX must be distinct.
c    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
c    this, but proper results will be computed if BIN_MAX < BIN_MIN.
c
c    Output, integer BIN(0:NBIN+1).
c    BIN(0) counts entries of A less than BIN_MIN.
c    BIN(NBIN+1) counts entries greater than or equal to BIN_MAX.
c    For 1 <= I <= NBIN, BIN(I) counts the entries X(J) such that
c      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
c    where H is the bin spacing.
c
c    Output, double precision BIN_LIMIT(0:NBIN), the "limits" of the bins.
c    BIN(I) counts the number of entries X(J) such that
c      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
c
      implicit none

      integer n
      integer nbin

      double precision a(n)
      integer bin(0:nbin+1)
      double precision bin_limit(0:nbin)
      double precision bin_max
      double precision bin_min
      integer i
      integer j
      double precision t

      if ( bin_max .eq. bin_min ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BIN - Fatal error!'
        write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
        stop
      end if

      do i = 0, nbin + 1
        bin(i) = 0
      end do

      do i = 1, n

        t = ( a(i) - bin_min ) / ( bin_max - bin_min )

        if ( t .lt. 0.0D+00 ) then
          j = 0
        else if ( 1.0D+00 .le. t ) then
          j = nbin + 1
        else
          j = 1 + int ( dble ( nbin ) * t )
        end if

        bin(j) = bin(j) + 1

      end do
c
c  Compute the bin limits.
c
      do i = 0, nbin
        bin_limit(i) = ( dble ( nbin - i ) * bin_min
     &                 + dble (        i ) * bin_max )
     &                 / dble ( nbin     )
      end do

      return
      end
      subroutine r8vec_bin_even ( n, a, nbin, bin_min, bin_max,
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R8VEC_BIN_EVEN bins an R8VEC into evenly spaced bins.
c
c  Discussion:
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the set.
c
c    Input, double precision A(N), the data to be binned.
c
c    Input, integer NBIN, the number of bins.
c
c    Input, double precision BIN_MIN, BIN_MAX, the bin limits.
c
c    Output, integer BIN_START(NBIN), BIN_LAST(NBIN), the index of the
c    first and last element of A that went into each bin, or -1 if there
c    are no entries in this bin.
c
c    Output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin

      double precision a(n)
      integer bin_last(nbin)
      integer bin_next(n)
      integer bin_start(nbin)
      double precision bin_max
      double precision bin_min
      integer i
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do i = 1, nbin
        bin_last(i) = -1
      end do
      do i = 1, nbin
        bin_start(i) = -1
      end do
      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do i = 1, n

        call r8_to_bin_even ( nbin, bin_min, bin_max, a(i), j )

        if ( bin_start(j) .eq. -1 ) then
          bin_start(j) = i
        else
          k = bin_last(j)
          bin_next(k) = i
        end if

        bin_next(i) = 0

        bin_last(j) = i

      end do

      return
      end
      subroutine r8vec_binned_reorder ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R8VEC_BINNED_REORDER reorders a binned R8VEC.
c
c  Discussion:
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START, BIN_LAST and BIN_NEXT arrays have also been
c    updated so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(N), the data to be sorted.
c
c    Input, integer NBIN, the number of bins.
c
c    Input/output, integer BIN_START(NBIN), BIN_LAST(NBIN), the index of
c    the first and last element of A that went into each bin, or -1 if
c    there are no entries in the bin.
c
c    Input/output, integer BIN_NEXT(N), the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer n
      integer nbin

      double precision a(n)
      double precision a2(n)
      integer bin_last(nbin)
      integer bin_next(n)
      integer bin_start(nbin)
      integer i
      integer j
      integer k

      k = 0

      do i = 1, nbin

        j = bin_start(i)

        if ( 0 .lt. j ) then
          bin_start(i) = k + 1
        end if

10      continue

        if ( 0 .lt. j ) then
          k = k + 1
          bin_last(i) = k
          a2(k) = a(j)
          j = bin_next(j)
          go to 10
        end if

      end do

      do i = 1, n
        a(i) = a2(i)
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i+1
      end do

      do i = 1, nbin
        k = bin_last(i)
        if ( 0 .lt. k ) then
          bin_next(k) = 0
        end if
      end do

      return
      end
      subroutine r8vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

c*********************************************************************72
c
cc R8VEC_BINNED_SORT_A ascending sorts a binned reordered R8VEC.
c
c  Discussion:
c
c    Presumably, the data vector was first binned by R8VEC_BIN_EVEN,
c    then reordered by R8VEC_BINNED_REORDER.  Now, each of the
c    bins of data is sorted one at a time, which results in sorting
c    the entire vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(N), the data to be sorted.
c
c    Input, integer NBIN, the number of bins.
c
c    Input, integer BIN_START(NBIN), BIN_LAST(NBIN), the index of the first
c    and last element of A that went into each bin, or -1 if there were no
c    elements in the bin.
c
      implicit none

      integer n
      integer nbin

      double precision a(n)
      integer bin_last(nbin)
      integer bin_start(nbin)
      integer i
      integer i1
      integer i2
      integer n1

      if ( n .le. 1 ) then
        return
      end if

      do i = 1, nbin

        i1 = bin_start(i)

        if ( 0 .lt. i1 ) then

          i2 = bin_last(i)

          n1 = i2 + 1 - i1

          call r8vec_sort_quick_a ( n1, a(i1:i2) )

        end if

      end do

      return
      end
      subroutine r8vec_bracket ( n, x, xval, left, right )

c*********************************************************************72
c
cc R8VEC_BRACKET searches a sorted array for successive brackets of a value.
c
c  Discussion:
c
c    If the values in the vector are thought of as defining intervals
c    on the real line, then this routine searches for the interval
c    nearest to or containing the given value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of input array.
c
c    Input, double precision X(N), an array that has been sorted into
c    ascending order.
c
c    Input, double precision XVAL, a value to be bracketed.
c
c    Output, integer LEFT, RIGHT, the results of the search.
c    Either:
c      XVAL < X(1), when LEFT = 1, RIGHT = 2;
c      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
c    or
c      X(LEFT) <= XVAL <= X(RIGHT).
c
      implicit none

      integer n

      integer i
      integer left
      integer right
      double precision x(n)
      double precision xval

      do i = 2, n - 1

        if ( xval .lt. x(i) ) then
          left = i - 1
          right = i
          return
        end if

       end do

      left = n - 1
      right = n

      return
      end
      function r8vec_eq ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_EQ is true if every pair of entries in two vectors is equal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision A1(N), A2(N), two vectors to compare.
c
c    Output, logical R8VEC_EQ.
c    R8VEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
c    and .FALSE. otherwise.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_eq

      r8vec_eq = .false.

      do i = 1, n
        if ( a1(i) .ne. a2(i) ) then
          return
        end if
      end do

      r8vec_eq = .true.

      return
      end
      function r8vec_gt ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_GT == ( A1 > A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 > A2  <=>                              A1(1) > A2(1) or
c                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
c                 ...
c                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision A1(N), A2(N), the vectors to be compared.
c
c    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_gt

      r8vec_gt = .false.

      do i = 1, n

        if ( a2(i) .lt. a1(i) ) then
          r8vec_gt = .true.
          return
        else if ( a1(i) .lt. a2(i) ) then
          return
        end if

      end do

      return
      end
      function r8vec_lt ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_LT == ( A1 < A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 < A2  <=>                              A1(1) < A2(1) or
c                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
c                 ...
c                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision A1(N), A2(N), the vectors to be compared.
c
c    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_lt

      r8vec_lt = .false.

      do i = 1, n

        if ( a1(i) .lt. a2(i) ) then
          r8vec_lt = .true.
          return
        else if ( a2(i) .lt. a1(i) ) then
          return
        end if

      end do

      return
      end
      subroutine r8vec_part_quick_a ( n, a, l, r )

c*********************************************************************72
c
cc R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
c
c  Discussion:
c
c    The routine reorders the entries of A.  Using A(1) as the key,
c    all entries of A that are less than or equal to the key will
c    precede the key which precedes all entries that are greater than the key.
c
c  Example:
c
c    Input:
c
c      N = 8
c
c      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
c
c    Output:
c
c      L = 3, R = 6
c
c      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
c            -------        -------
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of A.
c
c    Input/output, double precision A(N).  On input, the array to be checked.
c    On output, A has been reordered as described above.
c
c    Output, integer L, R, the indices of A that define the three segments.
c    Let KEY = the input value of A(1).  Then
c    I <= L                 A(I) < KEY;
c         L < I < R         A(I) = KEY;
c                 R <= I    KEY < A(I).
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision key
      integer l
      integer m
      integer r
      double precision temp

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_PART_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        l = 0
        r = 2
        return
      end if

      key = a(1)
      m = 1
c
c  The elements of unknown size have indices between L+1 and R-1.
c
      l = 1
      r = n + 1

      do i = 2, n

        if ( key .lt. a(l+1) ) then
          r = r - 1
          temp = a(r)
          a(r) = a(l+1)
          a(l+1) = temp
        else if ( a(l+1) .eq. key ) then
          m = m + 1
          temp = a(m)
          a(m) = a(l+1)
          a(l+1) = temp
          l = l + 1
        else if ( a(l+1) .lt. key ) then
          l = l + 1
        end if

      end do
c
c  Now shift small elements to the left, and KEY elements to center.
c
      do i = 1, l - m
        a(i) = a(i+m)
      end do
c
c  Out of bounds here, occasionally
c
      l = l - m

      do i = l + 1, l + m
        a(i) = key
      end do

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is an array of double precision real values.
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
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

      return
      end
      subroutine r8vec_sort_quick_a ( n, a )

c*********************************************************************72
c
cc R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
c
c  Example:
c
c    Input:
c
c      N = 7
c      A = ( 6, 7, 3, 2, 9, 1, 8 )
c
c    Output:
c
c      A = ( 1, 2, 3, 6, 7, 8, 9 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, double precision A(N).
c    On input, the array to be sorted.
c    On output, the array has been sorted.
c
      implicit none

      integer level_max
      parameter ( level_max = 25 )
      integer n

      double precision a(n)
      integer base
      integer l_segment
      integer level
      integer n_segment
      integer rsave(level_max)
      integer r_segment

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      else if ( n .eq. 1 ) then
        return
      end if

      level = 1
      rsave(level) = n + 1
      base = 1
      n_segment = n

10    continue
c
c  Partition the segment.
c
        call r8vec_part_quick_a ( n_segment, a(base), l_segment,
     &    r_segment )
c
c  If the left segment has more than one element, we need to partition it.
c
        if ( 1 .lt. l_segment ) then

          if ( level_max < level ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
            write ( *, '(a,i6)' )
     &        '  Exceeding recursion maximum of ', level_max
            stop
          end if

          level = level + 1
          n_segment = l_segment
          rsave(level) = r_segment + base - 1
c
c  The left segment and the middle segment are sorted.
c  Must the right segment be partitioned?
c
        else if ( r_segment .lt. n_segment ) then

          n_segment = n_segment + 1 - r_segment
          base = base + r_segment - 1
c
c  Otherwise, we back up a level if there is an earlier one.
c
        else

20        continue

            if ( level .le. 1 ) then
              go to 40
            end if

            base = rsave(level)
            n_segment = rsave(level-1) - rsave(level)
            level = level - 1

            if ( 0 .lt. n_segment ) then
              go to 30
            end if

          go to 20

30        continue

        end if

      go to 10

40    continue

      return
      end
      subroutine r8vec_swap ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_SWAP swaps two R8VEC's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the arrays.
c
c    Input/output, double precision A1(N), A2(N), the vectors to swap.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      double precision t

      do i = 1, n
        t     = a1(i)
        a1(i) = a2(i)
        a2(i) = t
      end do

      return
      end
      subroutine r8vec_to_bin_even3 ( ndim, nbin, a, b, c, bin )

c*********************************************************************72
c
cc R8VEC_TO_BIN_EVEN3 determines the appropriate "bin" for a R8VEC value.
c
c  Discussion:
c
c    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
c    or bins.
c
c  Example:
c
c    NDIM = 3
c    NBIN = (/ 4, 5, 2 /),
c
c      A(1) = 1,  A(2) = 0,  A(3) = 8
c      B(1) = 17, B(2) = 20, B(3) = 10
c
c
c            8 < Z < 9                    9 < Z < 10
c
c   20 +     +     +     +     +     20 +     +     +     +     +
c        151 | 251 | 351 | 451            152 | 252 | 352 | 452
c   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
c        141 | 241 | 341 | 441            142 | 242 | 342 | 442
c   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
c        131 | 231 | 331 | 431            132 | 232 | 332 | 432
c    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
c        121 | 221 | 321 | 421            122 | 222 | 322 | 422
c    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
c        111 | 211 | 311 | 411            112 | 212 | 312 | 412
c    0 +     +     +     +     +      0 +     +     +     +     +
c      1     5     9    13    17        1     5     9    13    17
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDIM, the dimension of the space.
c
c    Input, integer NBIN(NDIM), the number of bins in each dimension.
c
c    Input, double precision A(NDIM), B(NDIM), the lower and upper limits of
c    the bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(NDIM), a value to be placed in a bin.
c
c    Output, integer BIN(NDIM), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
      end do

      return
      end
      subroutine r8vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer M, the number of entries in the vector.
c
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original version by Albert Nijenhuis, Herbert Wilf.
c    This version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
      end if

      return
      end
      subroutine swapec ( i, top, btri, bedg, point_num, point_xy,
     &  tri_num, tri_vert, tri_nabe, stack, ierr )

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
c    29 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This version by John Burkardt.
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
c    Input/output, integer BTRI, BEDG; on input, if positive, are the
c    triangle and edge indices of a boundary edge whose updated indices
c    must be recorded.  On output, these may be updated because of swaps.
c
c    Input, intger POINT_NUM, the number of points.
c
c    Input, double precision POINT_XY(2,POINT_NUM), the coordinates of
c    the points.
c
c    Input, integer TRI_NUM, the number of triangles.
c
c    Input/output, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.
c    May be updated on output because of swaps.
c
c    Input/output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list;
c    negative values are used for links of the counter-clockwise linked
c    list of boundary edges;  May be updated on output because of swaps.
c
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

      integer point_num
      integer tri_num

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
      integer r
      integer s
      integer stack(point_num)
      integer swap
      integer t
      integer top
      integer tri_nabe(3,tri_num)
      integer tri_vert(3,tri_num)
      integer tt
      integer u
      double precision point_xy(2,point_num)
      double precision x
      double precision y
c
c    Determine whether triangles in stack are Delaunay, and swap
c    diagonal edge of convex quadrilateral if not.
c
      x = point_xy(1,i)
      y = point_xy(2,i)

10    continue

        if ( top .le. 0 ) then
          go to 40
        end if

        t = stack(top)
        top = top - 1

        if ( tri_vert(1,t) .eq. i ) then
          e = 2
          b = tri_vert(3,t)
        else if ( tri_vert(2,t) .eq. i ) then
          e = 3
          b = tri_vert(1,t)
        else
          e = 1
          b = tri_vert(2,t)
        end if

        a = tri_vert(e,t)
        u = tri_nabe(e,t)

        if ( tri_nabe(1,u) .eq. t ) then
          f = 1
          c = tri_vert(3,u)
        else if ( tri_nabe(2,u) .eq. t ) then
          f = 2
          c = tri_vert(1,u)
        else
          f = 3
          c = tri_vert(2,u)
        end if

        swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a),
     &    point_xy(1,c), point_xy(2,c), point_xy(1,b), point_xy(2,b) )

        if ( swap .eq. 1 ) then

          em1 = i4_wrap ( e - 1, 1, 3 )
          ep1 = i4_wrap ( e + 1, 1, 3 )
          fm1 = i4_wrap ( f - 1, 1, 3 )
          fp1 = i4_wrap ( f + 1, 1, 3 )

          tri_vert(ep1,t) = c
          tri_vert(fp1,u) = i
          r = tri_nabe(ep1,t)
          s = tri_nabe(fp1,u)
          tri_nabe(ep1,t) = u
          tri_nabe(fp1,u) = t
          tri_nabe(e,t) = s
          tri_nabe(f,u) = r

          if ( 0 .lt. tri_nabe(fm1,u) ) then
            top = top + 1
            stack(top) = u
          end if

          if ( 0 .lt. s ) then

            if ( tri_nabe(1,s) .eq. u ) then
              tri_nabe(1,s) = t
            else if ( tri_nabe(2,s) .eq. u ) then
              tri_nabe(2,s) = t
            else
              tri_nabe(3,s) = t
            end if

            top = top + 1

            if ( point_num .lt. top ) then
              ierr = 8
              go to 40
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

            if ( 0 .lt. tri_nabe(ee,tt) ) then

              tt = tri_nabe(ee,tt)

              if ( tri_vert(1,tt) .eq. a ) then
                ee = 3
              else if ( tri_vert(2,tt) .eq. a ) then
                ee = 1
              else
                ee = 2
              end if

              go to 20

            end if

            tri_nabe(ee,tt) = l

          end if

          if ( 0 .lt. r ) then

            if ( tri_nabe(1,r) .eq. t ) then
              tri_nabe(1,r) = u
            else if ( tri_nabe(2,r) .eq. t ) then
              tri_nabe(2,r) = u
            else
              tri_nabe(3,r) = u
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

            if ( 0 .lt. tri_nabe(ee,tt) ) then

              tt = tri_nabe(ee,tt)

              if ( tri_vert(1,tt) .eq. b ) then
                ee = 3
              else if ( tri_vert(2,tt) .eq. b ) then
                ee = 1
              else
                ee = 2
              end if

              go to 30

            end if

            tri_nabe(ee,tt) = l

          end if

        end if

      go to 10

40    continue

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
      subroutine triangle_area_2d ( t, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_2D computes the area of a triangle in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision AREA, the absolute area of the triangle.
c
      implicit none

      double precision area
      double precision t(2,3)

      area = 0.5D+00 * abs (
     &    t(1,1) * ( t(2,2) - t(2,3) )
     &  + t(1,2) * ( t(2,3) - t(2,1) )
     &  + t(1,3) * ( t(2,1) - t(2,2) ) )

      return
      end
      subroutine triangle_sample_2d ( t, seed, p )

c*********************************************************************72
c
cc TRIANGLE_SAMPLE_2D returns a random point in a triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision P(2), a random point in the triangle.
c
      implicit none

      double precision alpha
      double precision beta
      double precision p(2)
      double precision r
      double precision r8_uniform_01
      integer seed
      double precision t(2,3)
      double precision x12
      double precision x13
      double precision y12
      double precision y13

      r = r8_uniform_01 ( seed )
c
c  Interpret R as a percentage of the triangle's area.
c
c  Imagine a line L, parallel to side 1, so that the area between
c  vertex 1 and line L is R percent of the full triangle's area.
c
c  The line L will intersect sides 2 and 3 at a fraction
c  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
c
      alpha = sqrt ( r )
c
c  Determine the coordinates of the points on sides 2 and 3 intersected
c  by line L.
c
      x12 = alpha * t(1,1) + ( 1.0D+00 - alpha ) * t(1,2)
      y12 = alpha * t(2,1) + ( 1.0D+00 - alpha ) * t(2,2)

      x13 = alpha * t(1,1) + ( 1.0D+00 - alpha ) * t(1,3)
      y13 = alpha * t(2,1) + ( 1.0D+00 - alpha ) * t(2,3)
c
c  Now choose, uniformly at random, a point on the line L.
c
      beta = r8_uniform_01 ( seed )

      p(1) = beta * x12 + ( 1.0D+00 - beta ) * x13
      p(2) = beta * y12 + ( 1.0D+00 - beta ) * y13

      return
      end
      subroutine triangulation_nabe_nodes ( point_num, tri_num,
     &  tri_vert, nabes_first, nabes_num, nabes_max, nabes_dim,
     &  nabes )

c*********************************************************************72
c
cc TRIANGULATION_NABE_NODES determines the neighbors of triangulation nodes.
c
c  Example:
c
c    On input, the triangle data structure is:
c
c    Triangle  Nodes
c    --------  ----------
c     1        3,   4,   1
c     2        3,   1,   2
c     3        3,   2,   6
c     4        2,   1,   5
c     5        6,   2,   5
c
c  On output, the auxilliary neighbor arrays are:
c
c    Node  Num  First
c    ----  ---  -----
c     1     4     1
c     2     4     5
c     3     4     9
c     4     2    13
c     5     3    15
c     6     3    18
c
c  and the neighbor array is:
c
c    Position  Node
c    --------  ----
c
c     1        2
c     2        3
c     3        4
c     4        5
c    -----------
c     5        1
c     6        3
c     7        5
c     8        6
c    -----------
c     9        1
c    10        2
c    11        4
c    12        6
c    -----------
c    13        1
c    14        3
c    -----------
c    15        1
c    16        2
c    17        6
c    -----------
c    18        2
c    19        3
c    20        5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer TRI_NUM, the number of triangles.
c
c    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
c
c    Output, integer NABES_FIRST(POINT_NUM), the index in NABES of the first
c    neighbor in the list for each node.
c
c    Output, integer NABES_NUM(POINT_NUM), the number of neighbors of each node.
c
c    Input, integer NABES_MAX, the maximum dimension of NABES.
c
c    Output, integer NABES_DIM, the dimension of NABES.
c
c    Output, integer NABES(NABES_DIM), a list of the neighbors of all the nodes.
c    Neighbors of node 1 are listed first, and so on.
c
      implicit none

      integer nabes_max
      integer point_num
      integer tri_num

      integer i
      integer i_current
      integer j
      integer k
      integer nabe
      integer nabes(nabes_max)
      integer nabes1(nabes_max)
      integer nabes_dim
      integer nabes_first(point_num)
      integer nabes_num(point_num)
      integer nuniq
      integer tri
      integer tri_vert(3,tri_num)
c
c  Step 1.  From the triangle list (I,J,K)
c  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
c
      nabes_dim = 0

      do tri = 1, tri_num

        i = tri_vert(1,tri)
        j = tri_vert(2,tri)
        k = tri_vert(3,tri)

        nabes1(nabes_dim+1) = i
        nabes1(nabes_dim+2) = i
        nabes1(nabes_dim+3) = j
        nabes1(nabes_dim+4) = j
        nabes1(nabes_dim+5) = k
        nabes1(nabes_dim+6) = k

        nabes(nabes_dim+1) = j
        nabes(nabes_dim+2) = k
        nabes(nabes_dim+3) = i
        nabes(nabes_dim+4) = k
        nabes(nabes_dim+5) = i
        nabes(nabes_dim+6) = j

        nabes_dim = nabes_dim + 6
      end do
c
c  Step 2. Dictionary sort the neighbor relations.
c
      call i4vec2_sort_a ( nabes_dim, nabes1, nabes )
c
c  Step 3. Remove duplicate entries.
c
      call i4vec2_sorted_unique ( nabes_dim, nabes1, nabes, nuniq )

      nabes_dim = nuniq
c
c  Step 4. Construct the NABES_NUM and NABES_FIRST data.
c
      do i = 1, point_num
        nabes_num(i) = 0
      end do

      do i = 1, point_num
        nabes_first(i) = 0
      end do

      i_current = 0
      do nabe = 1, nabes_dim
        i = nabes1(nabe)
        if ( i .eq. i_current ) then
          nabes_num(i) = nabes_num(i) + 1
        else
          i_current = i
          nabes_first(i) = nabe
          nabes_num(i) = 1
        end if
      end do

      return
      end
      subroutine triangulation_print ( point_num, tri_num, xc,
     &  tri_vert, tri_nabe, vertex_list )

c*********************************************************************72
c
cc TRIANGULATION_PRINT prints out information defining a Delaunay triangulation.
c
c  Discussion:
c
c    Triangulations created by RTRIS include extra information encoded
c    in the negative values of TRI_NABE.
c
c    Because some of the nodes counted in POINT_NUM may not actually be
c    used in the triangulation, I needed to compute the true number
c    of vertices.  I added this calculation on 13 October 2001.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer TRI_NUM, the number of triangles.
c
c    Input, double precision XC(2,POINT_NUM), the point coordinates.
c
c    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up the triangles.
c
c    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbors on each side.
c    If there is no triangle neighbor on a particular side, the value of
c    TRI_NABE should be negative.  If the triangulation data was created by
c    DTRIS2, then there is more information encoded in the negative values.
c
c    Workspace, integer VERTEX_LIST(3*TRI_NUM).
c
      implicit none

      integer point_num
      integer tri_num

      integer boundary_num
      integer i
      integer i4_wrap
      integer j
      integer k
      integer n1
      integer n2
      integer s
      integer t
      integer tri
      integer tri_nabe(3,tri_num)
      integer tri_vert(3,tri_num)
      integer vertex_list(3*tri_num)
      integer vertex_num
      double precision xc(2,point_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_PRINT'
      write ( *, '(a)' ) '  Information defining a triangulation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of points is ', point_num

      call r8mat_transpose_print ( 2, point_num, xc,
     &  '  Point coordinates' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of triangles is ', tri_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Sets of three points are used as vertices of'
      write ( *, '(a)' )
     &  '  the triangles.  For each triangle, the points'
      write ( *, '(a)' )
     &  '  are listed in counterclockwise order.'

      call i4mat_transpose_print ( 3, tri_num, tri_vert,
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

      call i4mat_transpose_print ( 3, tri_num, tri_nabe,
     &  '  Triangle neighbors' )
c
c  Determine the number of vertices.  This is not the same as the
c  number of pointsc
c
      k = 0
      do tri = 1, tri_num
        do i = 1, 3
          k = k + 1
          vertex_list(k) = tri_vert(i,tri)
        end do
      end do

      call i4vec_sort_heap_a ( 3*tri_num, vertex_list )

      call i4vec_sorted_unique ( 3*tri_num, vertex_list, vertex_num )
c
c  Determine the number of boundary points.
c
      boundary_num = 2 * vertex_num - tri_num - 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The number of boundary points is ',
     &  boundary_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  The segments that make up the convex hull can be'
      write ( *, '(a)' )
     &  '  determined from the negative entries of the triangle'
      write ( *, '(a)' ) '  neighbor list.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  # Tri Side  N1  N2'
      write ( *, '(a)' ) ' '

      k = 0

      do i = 1, tri_num

        do j = 1, 3

          if ( tri_nabe(j,i) .lt. 0 ) then
            s = - tri_nabe(j,i)
            t = s / 3

            if ( t .lt. 1 .or. tri_num .lt. t ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' )
     &          '  Sorry, this data does not use the DTRIS2'
              write ( *, '(a)' )
     &          '  convention for convex hull segments.'
              go to 10
            end if

            s = mod ( s, 3 ) + 1
            k = k + 1
            n1 = tri_vert(s,t)
            n2 = tri_vert(i4_wrap(s+1,1,3),t)
            write ( *, '(5i4)' ) k, t, s, n1, n2
          end if

        end do

      end do

10    continue

      return
      end
      subroutine triangulation_sample_2d ( point_num, xc, tri_num,
     &  tri_vert, num_ran, seed, xd, td )

c*********************************************************************72
c
cc TRIANGULATION_SAMPLE_2D returns random points in a triangulation.
c
c  Discussion:
c
c    It is assumed that the triangulation consists of a set of non-overlapping
c    triangles.
c
c    The point is chosen uniformly in the area covered by the triangulation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points used in the triangulation.
c
c    Input, double precision XC(2,POINT_NUM), the point coordinates.
c
c    Input, integer TRI_NUM, the number of triangles.
c
c    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up the triangles.
c
c    Input, integer NUM_RAN, the number of points to sample.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision XD(2,NUM_RAN), the sample points.
c
c    Output, integer TD(NUM_RAN), the triangle to which each sample point
c    belongs.
c
      implicit none

      integer point_num
      integer num_ran
      integer tri_num

      double precision area
      double precision area_cum(0:tri_num)
      double precision area_total
      integer i
      integer i1
      integer i2
      integer i3
      integer left
      double precision r
      double precision r8_uniform_01
      integer right
      integer seed
      double precision t(2,3)
      integer td(num_ran)
      integer tri_vert(3,tri_num)
      double precision xc(2,point_num)
      double precision xd(2,num_ran)
c
c  Compute the areas of the triangles.
c  Build a cumulative area vector.
c  Convert it to a relative cumulative area vector.
c
      area_cum(0) = 0.0D+00

      do i = 1, tri_num

        i1 = tri_vert(1,i)
        i2 = tri_vert(2,i)
        i3 = tri_vert(3,i)

        t(1,1) = xc(1,i1)
        t(2,1) = xc(2,i1)
        t(1,2) = xc(1,i2)
        t(2,2) = xc(2,i2)
        t(1,3) = xc(1,i3)
        t(2,3) = xc(2,i3)

        call triangle_area_2d ( t, area )

        area_cum(i) = area_cum(i-1) + area

      end do

      area_total = area_cum(tri_num)

      do i = 0, tri_num
        area_cum(i) = area_cum(i) / area_total
      end do
c
c  Pick random values.  A random value R indicates the corresponding triangle
c  whose cumulative relative area contains R.
c
c  Bracket the random value in the cumulative relative areas,
c  indicating a triangle.
c
c  Pick a random point in the triangle.
c
      do i = 1, num_ran

        r = r8_uniform_01 ( seed )

        call r8vec_bracket ( tri_num+1, area_cum, r, left, right )

        td(i) = right - 1

        i1 = tri_vert(1,td(i))
        i2 = tri_vert(2,td(i))
        i3 = tri_vert(3,td(i))

        t(1,1) = xc(1,i1)
        t(2,1) = xc(2,i1)
        t(1,2) = xc(1,i2)
        t(2,2) = xc(2,i2)
        t(1,3) = xc(1,i3)
        t(2,3) = xc(2,i3)

        call triangle_sample_2d ( t, seed, xd(1,i) )

      end do

      return
      end
      subroutine tuple_next2 ( n, xmin, xmax, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT2 computes the next element of an integer tuple space.
c
c  Discussion:
c
c    The elements X are N vectors.
c
c    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
c
c    The elements are produced one at a time.
c
c    The first element is
c      (XMIN(1), XMIN(2), ..., XMIN(N)),
c    the second is (probably)
c      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
c    and the last element is
c      (XMAX(1), XMAX(2), ..., XMAX(N))
c
c    Intermediate elements are produced in a lexicographic order, with
c    the first index more important than the last, and the ordering of
c    values at a fixed index implicitly defined by the sign of
c    XMAX(I) - XMIN(I).
c
c  Example:
c
c    N = 2,
c    XMIN = (/ 1, 10 /)
c    XMAX = (/ 3,  8 /)
c
c    RANK    X
c    ----  -----
c      1   1 10
c      2   1  9
c      3   1  8
c      4   2 10
c      5   2  9
c      6   2  8
c      7   3 10
c      8   3  9
c      9   3  8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components.
c
c    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
c    These values are minimum and maximum only in the sense of the lexicographic
c    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
c    than XMAX(I).
c
c    Input/output, integer RANK, the rank of the item.  On first call,
c    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
c    there are no more items in the sequence.
c
c    Input/output, integer X(N), on input the previous tuple.
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer prod
      integer rank
      integer x(n)
      integer xmin(n)
      integer xmax(n)

      if ( rank .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      prod = 1
      do i = 1, n
        prod = prod * ( 1 + abs ( xmax(i) - xmin(i) ) )
      end do

      if ( prod .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      if ( rank .eq. 0 ) then
        do i = 1, n
          x(i) = xmin(i)
        end do
        rank = 1
        return
      end if

      rank = rank + 1
      i = n

10    continue

        if ( x(i) .ne. xmax(i) ) then
          x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
          go to 20
        end if

        x(i) = xmin(i)

        if ( i .eq. 1 ) then
          rank = 0
          go to 20
        end if

        i = i - 1

      go to 10

20    continue

      return
      end
      subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert,
     &  tri_nabe, ltri, ledg, rtri, redg )

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
c    31 December 2008
c
c  Author:
c
c    Original FORTRAN77 version by Barry Joe.
c    This version by John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Modified:
c
c    29 July 2008
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point outside the
c    convex hull of the current triangulation.
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, double precision POINT_XY(2,POINT_NUM), the coordinates of the
c    vertices.
c
c    Input, integer TRI_NUM, the number of triangles.
c
c    Input, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.
c
c    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list; negative
c    values are used for links of a counter clockwise linked list of boundary
c    edges;
c      LINK = -(3*I + J-1) where I, J = triangle, edge index.
c
c    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these values are
c    assumed to be already computed and are not changed, else they are updated.
c    On output, LTRI is the index of boundary triangle to the left of the
c    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
c    edge of triangle LTRI to the left of the leftmost boundary edge visible
c    from (X,Y).  1 <= LEDG <= 3.
c
c    Input/output, integer RTRI.  On input, the index of the boundary triangle
c    to begin the search at.  On output, the index of the rightmost boundary
c    triangle visible from (X,Y).
c
c    Input/output, integer REDG, the edge of triangle RTRI that is visible
c    from (X,Y).  1 <= REDG <= 3.
c
      implicit none

      integer point_num
      integer tri_num

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
      double precision point_xy(2,point_num)
      integer redg
      integer rtri
      integer t
      integer tri_nabe(3,tri_num)
      integer tri_vert(3,tri_num)
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

        l = -tri_nabe(redg,rtri)
        t = l / 3
        e = mod ( l, 3 ) + 1
        a = tri_vert(e,t)

        if ( e .le. 2 ) then
          b = tri_vert(e+1,t)
        else
          b = tri_vert(1,t)
        end if

        lr = lrline ( x, y, point_xy(1,a), point_xy(2,a),
     &    point_xy(1,b), point_xy(2,b), 0.0D+00 )

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

        b = tri_vert(e,t)
        e = i4_wrap ( e-1, 1, 3 )

40      continue

        if ( 0 .lt. tri_nabe(e,t) ) then

          t = tri_nabe(e,t)

          if ( tri_vert(1,t) .eq. b ) then
            e = 3
          else if ( tri_vert(2,t) .eq. b ) then
            e = 1
          else
            e = 2
          end if

          go to 40

        end if

        a = tri_vert(e,t)

        lr = lrline ( x, y, point_xy(1,a), point_xy(2,a),
     &    point_xy(1,b), point_xy(2,b), 0.0D+00 )

        if ( lr .le. 0 ) then
          go to 50
        end if

      go to 30

50    continue

      ltri = t
      ledg = e

      return
      end
