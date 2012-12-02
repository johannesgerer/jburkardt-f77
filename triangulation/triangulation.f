      subroutine alpha_measure ( n, z, triangle_order, triangle_num, 
     &  triangle_node, alpha_min, alpha_ave, alpha_area )

c*********************************************************************72
c
cc ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
c
c  Discusion:
c
c    The ALPHA measure evaluates the uniformity of the shapes of the triangles
c    defined by a triangulated pointset.
c
c    We compute the minimum angle among all the triangles in the triangulated
c    dataset and divide by the maximum possible value (which, in degrees,
c    is 60).  The best possible value is 1, and the worst 0.  A good
c    triangulation should have an ALPHA score close to 1.
c
c    The code has been modified to 'allow' 6-node triangulations.
c    However, no effort is made to actually process the midside nodes.
c    Only information from the vertices is used.
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
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision Z(2,N), the points.
c
c    Input, integer TRIANGLE_ORDER, the order of the triangles.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
c    the triangulation.
c
c    Output, double precision ALPHA_MIN, the minimum value of ALPHA over all
c    triangles.
c
c    Output, double precision ALPHA_AVE, the value of ALPHA averaged over
c    all triangles.
c
c    Output, double precision ALPHA_AREA, the value of ALPHA averaged over
c    all triangles and weighted by area.
c
      implicit none

      integer n
      integer triangle_num
      integer triangle_order

      double precision a_angle
      integer a_index
      double precision a_x
      double precision a_y
      double precision ab_len
      double precision alpha
      double precision alpha_area
      double precision alpha_ave
      double precision alpha_min
      double precision area
      double precision area_total
      double precision b_angle
      integer b_index
      double precision b_x
      double precision b_y
      double precision bc_len
      double precision c_angle
      integer c_index
      double precision c_x
      double precision c_y
      double precision ca_len
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_acos
      double precision r8_huge
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      double precision z(2,n)

      alpha_min = r8_huge ( )
      alpha_ave = 0.0D+00
      alpha_area = 0.0D+00
      area_total = 0.0D+00

      do triangle = 1, triangle_num

        a_index = triangle_node(1,triangle)
        b_index = triangle_node(2,triangle)
        c_index = triangle_node(3,triangle)

        a_x = z(1,a_index)
        a_y = z(2,a_index)
        b_x = z(1,b_index)
        b_y = z(2,b_index)
        c_x = z(1,c_index)
        c_y = z(2,c_index)

        area = 0.5D+00 * abs ( a_x * ( b_y - c_y ) 
     &                       + b_x * ( c_y - a_y ) 
     &                       + c_x * ( a_y - b_y ) )

        ab_len = sqrt ( ( a_x - b_x )**2 + ( a_y - b_y )**2 )
        bc_len = sqrt ( ( b_x - c_x )**2 + ( b_y - c_y )**2 )
        ca_len = sqrt ( ( c_x - a_x )**2 + ( c_y - a_y )**2 )
c
c  Take care of a ridiculous special case.
c
        if ( ab_len .eq. 0.0D+00 .and. 
     &       bc_len .eq. 0.0D+00 .and. 
     &       ca_len .eq. 0.0D+00 ) then

          a_angle = 2.0D+00 * pi / 3.0D+00
          b_angle = 2.0D+00 * pi / 3.0D+00
          c_angle = 2.0D+00 * pi / 3.0D+00

        else

          if ( ca_len .eq. 0.0D+00 .or. ab_len .eq. 0.0D+00 ) then
            a_angle = pi
          else
            a_angle = 
     &        r8_acos ( ( ca_len**2 + ab_len**2 - bc_len**2 ) 
     &        / ( 2.0D+00 * ca_len * ab_len ) )
          end if

          if ( ab_len .eq. 0.0D+00 .or. bc_len .eq. 0.0D+00 ) then
            b_angle = pi
          else
            b_angle = 
     &        r8_acos ( ( ab_len**2 + bc_len**2 - ca_len**2 ) 
     &        / ( 2.0D+00 * ab_len * bc_len ) )
          end if

          if ( bc_len .eq. 0.0D+00 .or. ca_len .eq. 0.0D+00 ) then
            c_angle = pi
          else
            c_angle = 
     &        r8_acos ( ( bc_len**2 + ca_len**2 - ab_len**2 ) 
     &        / ( 2.0D+00 * bc_len * ca_len ) )
          end if

        end if

        alpha_min = min ( alpha_min, a_angle )
        alpha_min = min ( alpha_min, b_angle )
        alpha_min = min ( alpha_min, c_angle )

        alpha_ave = alpha_ave + alpha_min

        alpha_area = alpha_area + area * alpha_min

        area_total = area_total + area

      end do

      alpha_ave = alpha_ave / dble ( triangle_num )
      alpha_area = alpha_area / area_total
c
c  Normalize angles from [0,pi/3] degrees into qualities in [0,1].
c
      alpha_min = alpha_min * 3.0D+00 / pi
      alpha_ave = alpha_ave * 3.0D+00 / pi
      alpha_area = alpha_area * 3.0D+00 / pi

      return
      end
      function angle_rad_2d ( p1, p2, p3 )

c*********************************************************************72
c
cc ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
c
c  Discussion:
c
c    Except for the zero angle case, it should be true that
c
c      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
c
c        P1
c        /
c       /
c      /
c     /
c    P2--------->P3
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
c    Input, double precision P1(2), P2(2), P3(2), define the rays
c    P1 - P2 and P3 - P2 which define the angle.
c
c    Output, double precision ANGLE_RAD_2D, the angle swept out by the rays,
c    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
c    length, then ANGLE_RAD_2D is set to 0.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle_rad_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)

      p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) 
     &     + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )


      p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) 
     &     - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

      if ( p(1) .eq. 0.0D+00 .and. p(2) .eq. 0.0D+00 ) then
        angle_rad_2d = 0.0D+00
        return
      end if

      angle_rad_2d = atan2 ( p(2), p(1) )

      if ( angle_rad_2d .lt. 0.0D+00 ) then
        angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
      end if

      return
      end
      subroutine area_measure ( n, z, triangle_order, triangle_num, 
     &  triangle_node, area_min, area_max, area_ratio, area_ave, 
     &  area_std )

c*********************************************************************72
c
cc AREA_MEASURE determines the area ratio quality measure.
c
c  Discusion:
c
c    This measure computes the area of every triangle, and returns
c    the ratio of the minimum to the maximum triangle.  A value of
c    1 is "perfect", indicating that all triangles have the same area.
c    A value of 0 is the worst possible result.
c
c    The code has been modified to 'allow' 6-node triangulations.
c    However, no effort is made to actually process the midside nodes.
c    Only information from the vertices is used.
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
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision Z(2,N), the points.
c
c    Input, integer TRIANGLE_ORDER, the order of the triangles.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
c    the triangulation.
c
!    Output, double precision AREA_MIN, AREA_MAX, the minimum and maximum 
!    areas.
!
!    Output, double precision AREA_RATIO, the ratio of the minimum to the
!    maximum area.
!
!    Output, double precision AREA_AVE, the average area.
!
!    Output, double precision AREA_STD, the standard deviation of the areas.
!
      implicit none

      integer n
      integer triangle_num
      integer triangle_order

      double precision area
      double precision area_ave
      double precision area_max
      double precision area_min
      double precision area_ratio
      double precision area_std
      double precision r8_huge
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      double precision x1
      double precision x2
      double precision x3
      double precision y1
      double precision y2
      double precision y3
      double precision z(2,n)

      area_max = 0.0D+00
      area_min = r8_huge ( )
      area_ave = 0.0D+00

      do triangle = 1, triangle_num

        x1 = z(1,triangle_node(1,triangle))
        y1 = z(2,triangle_node(1,triangle))
        x2 = z(1,triangle_node(2,triangle))
        y2 = z(2,triangle_node(2,triangle))
        x3 = z(1,triangle_node(3,triangle))
        y3 = z(2,triangle_node(3,triangle))

        area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) 
     &                       + x2 * ( y3 - y1 ) 
     &                       + x3 * ( y1 - y2 ) )

        area_min = min ( area_min, area )
        area_max = max ( area_max, area )

        area_ave = area_ave + area

      end do

      area_ave = area_ave / dble ( triangle_num )

      area_std = 0.0D+00
      do triangle = 1, triangle_num

        x1 = z(1,triangle_node(1,triangle))
        y1 = z(2,triangle_node(1,triangle))
        x2 = z(1,triangle_node(2,triangle))
        y2 = z(2,triangle_node(2,triangle))
        x3 = z(1,triangle_node(3,triangle))
        y3 = z(2,triangle_node(3,triangle))

        area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) 
     &                       + x2 * ( y3 - y1 ) 
     &                       + x3 * ( y1 - y2 ) )

        area_std = area_std + ( area - area_ave )**2
      end do
      area_std = sqrt ( area_std / dble ( triangle_num ) )

      area_ratio = area_min / area_max

      return
      end
      subroutine bandwidth ( element_order, element_num, element_node, 
     &  ml, mu, m )

c*********************************************************************72
c
cc BANDWIDTH determines the bandwidth associated with a finite element mesh.
c
c  Discussion:
c
c    The quantity computed here is the "geometric" bandwidth determined
c    by the finite element mesh alone.
c
c    If a single finite element variable is associated with each node
c    of the mesh, and if the nodes and variables are numbered in the
c    same way, then the geometric bandwidth is the same as the bandwidth
c    of a typical finite element matrix.
c
c    The bandwidth M is defined in terms of the lower and upper bandwidths:
c
c      M = ML + 1 + MU
c
c    where
c
c      ML = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but earlier column,
c
c      MU = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but later column.
c
c    Because the finite element node adjacency relationship is symmetric,
c    we are guaranteed that ML = MU.
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
c    Input, integer ELEMENT_ORDER, the order of the elements.
c
c    Input, integer ELEMENT_NUM, the number of elements.
c
c    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
c    ELEMENT_NODE(I,J) is the global index of local node I in element J.
c
c    Output, integer ML, MU, the lower and upper bandwidths 
c    of the matrix.
c
c    Output, integer M, the bandwidth of the matrix.
c
      implicit none

      integer element_num
      integer element_order

      integer element
      integer element_node(element_order,element_num)
      integer global_i
      integer global_j
      integer local_i
      integer local_j
      integer m
      integer ml
      integer mu

      ml = 0
      mu = 0

      do element = 1, element_num

        do local_i = 1, element_order
          global_i = element_node(local_i,element)

          do local_j = 1, element_order
            global_j = element_node(local_j,element)

            mu = max ( mu, global_j - global_i )
            ml = max ( ml, global_i - global_j )

          end do
        end do
      end do

      m = ml + 1 + mu

      return
      end
      subroutine delaunay_swap_test ( xy, swap )

c*********************************************************************72
c
cc DELAUNAY_SWAP_TEST performs the Delaunay swap test.
c
c  Discussion:
c
c    The current triangles are formed by nodes (1,2,3) and (1,3,4).
c    if a swap is recommended, the new triangles will be (1,2,4) and (2,3,4).
c
c      4     ?     4
c     / \         /|\
c    1---3  ==>  1 | 3
c     \ /         \|/
c      2           2
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
c  Reference:
c
c    Graham Carey,
c    Computational Grids:
c    Generation, Adaptation and Solution Strategies,
c    Taylor and Francis, 1997,
c    ISBN13: 978-1560326359,
c    LC: QA377.C32.
c
c  Parameters:
c
c    Input, double precision XY(2,4), the coordinates of four points.
c
c    Output, logical SWAP, is TRUE if the triangles (1,2,4) and (2,3,4)
c    are to replace triangles (1,2,3) and (1,3,4).
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      logical swap
      double precision x13
      double precision x14
      double precision x23
      double precision x24
      double precision xy(2,4)
      double precision y13
      double precision y14
      double precision y23
      double precision y24

      x13 = xy(1,1) - xy(1,3)
      x14 = xy(1,1) - xy(1,4)
      x23 = xy(1,2) - xy(1,3)
      x24 = xy(1,2) - xy(1,4)

      y13 = xy(2,1) - xy(2,3)
      y14 = xy(2,1) - xy(2,4)
      y23 = xy(2,2) - xy(2,3)
      y24 = xy(2,2) - xy(2,4)

      a = x13 * x23 + y13 * y23
      b = x24 * y14 - x14 * y24
      c = x23 * y13 - x13 * y23
      d = x24 * x14 + y14 * y24
c
c  The reference gives two initial tests before the
c  main one.  However, there seems to be an error
c  in at least one of these tests.  Since they are
c  intended to avoid error in borderline cases, but
c  instead cause real error in common cases, they are
c  omitted for now.
c
c if ( 0.0D+00 .le. a .and. 0.0D+00 .le. d ) then
c   swap = .true.
c else if ( a .lt. d .and. d .lt. 0.0D+00 ) then
c   swap = .true.
c  else if...

      if ( a * b .lt. c * d ) then
        swap = .true.
      else
        swap = .false.
      end if

      return
      end
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
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

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
      subroutine i4col_compare ( m, n, a, i, j, isgn )

c*********************************************************************72
c
cc I4COL_COMPARE compares columns I and J of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, I = 2, J = 4
c
c      A = (
c        1  2  3  4
c        5  6  7  8
c        9 10 11 12 )
c
c    Output:
c
c      ISGN = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of N columns of 
c    vectors of length M.
c
c    Input, integer I, J, the columns to be compared.
c    I and J must be between 1 and N.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, column I < column J,
c     0, column I = column J,
c    +1, column J < column I.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer isgn
      integer j
      integer k
c
c  Check.
c
      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  Column index I = ', i, ' is less than 1.'
        stop
      end if

      if ( n .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  N = ', n, ' is less than column index I = ', i
        stop
      end if

      if ( j .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  Column index J = ', j, ' is less than 1.'
        stop
      end if

      if ( n .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  N = ', n, ' is less than column index J = ', j
        stop
      end if

      isgn = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. m ) then

        if ( a(k,i) .lt. a(k,j) ) then
          isgn = -1
          return
        else if ( a(k,j) .lt. a(k,i) ) then
          isgn = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4col_sort_a ( m, n, a )

c*********************************************************************72
c
cc I4COL_SORT_A ascending sorts an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    In lexicographic order, the statement "X < Y", applied to two real
c    vectors X and Y of length M, means that there is some index I, with
c    1 <= I <= M, with the property that
c
c      X(J) = Y(J) for J < I,
c    and
c      X(I) < Y(I).
c
c    In other words, the first time they differ, X is smaller.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A, and the length of
c    a vector of data.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of N columns of M-vectors.
c    On output, the columns of A have been sorted in ascending
c    lexicographic order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer indx
      integer isgn
      integer j

      if ( m .le. 0 ) then
        return
      end if

      if ( n .le. 1 ) then
        return
      end if
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

          call i4col_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4col_compare ( m, n, a, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
c
c  Discussion:
c
c    The columns of the array may be ascending or descending sorted.
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
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), a sorted array, containing
c    N columns of data.
c
c    Output, integer UNIQUE_NUM, the number of unique columns.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j1
      integer j2
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1
      j1 = 1

      do j2 = 2, n

        do i = 1, m
          if ( a(i,j1) .ne. a(i,j2) ) then
            unique_num = unique_num + 1
            j1 = j2
            go to 10
          end if
        end do

10      continue

      end do

      return
      end
      subroutine i4col_swap ( m, n, a, j1, j2 )

c*********************************************************************72
c
cc I4COL_SWAP swaps columns J1 and J2 of an I4COL.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, J1 = 2, J2 = 4
c
c      A = (
c        1  2  3  4
c        5  6  7  8
c        9 10 11 12 )
c
c    Output:
c
c      A = (
c        1  4  3  2
c        5  8  7  6
c        9 12 11 10 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the array.
c
c    Input/output, integer A(M,N), an array of N columns 
c    of length M.
c
c    Input, integer J1, J2, the columns to be swapped.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j1
      integer j2
      integer t

      if ( j1 .lt. 1 .or. n .lt. j1 .or. 
     &     j2 .lt. 1 .or. n .lt. j2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
        write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
        write ( *, '(a,i8)' ) '  J1 =    ', j1
        write ( *, '(a,i8)' ) '  J2 =    ', j2
        write ( *, '(a,i8)' ) '  N =     ', n
        stop

      end if

      if ( j1 .eq. j2 ) then
        return
      end if

      do i = 1, m
        t       = a(i,j1)
        a(i,j1) = a(i,j2)
        a(i,j2) = t
      end do

      return
      end
      subroutine i4i4_sort_a ( i1, i2, j1, j2 )

c*********************************************************************72
c
cc I4I4_SORT_A ascending sorts a pair of integers.
c
c  Discussion:
c
c    An I4I4 is a pair of integers, regarded as a single data item.
c
c    The program allows the reasonable call:
c
c      call i4i4_sort_a ( i1, i2, i1, i2 )
c 
c    and this will return the reasonable result.
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
c    Input, integer I1, I2, the values to sort.
c
c    Output, integer J1, J2, the sorted values.
c
      implicit none

      integer i1
      integer i2
      integer j1
      integer j2
      integer k1
      integer k2
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4_sort_a ( i1, i2, i1, i2 )
!
      k1 = i1
      k2 = i2

      j1 = min ( k1, k2 )
      j2 = max ( k1, k2 )

      return
      end
      subroutine i4mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc I4MAT_COPY copies an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of I4 values.
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
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, integer A1(M,N), the matrix to copy.
c
c    Output, integer A2(M,N), the copy.
c
      implicit none

      integer m
      integer n

      integer a1(m,n)
      integer a2(m,n)
      integer i
      integer j
 
      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
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
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11, 12, 13, 14, 15 ).
c
c    Output:
c
c      A = ( 15, 14, 13, 12, 11 ).
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
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N), the array to be reversed.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer t

      do i = 1, n/2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
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
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data items.
c
c    Input, integer A1(N), A2(N), contain the two components 
c    of each item.
c
c    Input, integer I, J, the items to be compared.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, item I .lt. item J,
c     0, item I = item J,
c    +1, item J .lt. item I.
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
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items of data.
c
c    Input/output, integer A1(N), A2(N), the data to be sorted.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i
      integer indx
      integer isgn
      integer j
      integer temp

      if ( n .le. 1 ) then
        return
      end if
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

          temp  = a1(i)
          a1(i) = a1(j)
          a1(j) = temp

          temp  = a2(i)
          a2(i) = a2(j)
          a2(j) = temp
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
      subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

c*********************************************************************72
c
cc I4VEC2_SORTED_UNIQUE gets the unique elements in a sorted I4VEC2.
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
c    03 June 2009
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
c    On output, an array of unique items.
c
c    Output, integer UNIQUE_NUM, the number of unique items.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer itest
      integer unique_num

      unique_num = 0

      if ( n .le. 0 ) then
        return
      end if

      unique_num = 1

      do itest = 2, n

        if ( a1(itest) .ne. a1(unique_num) .or. 
     &       a2(itest) .ne. a2(unique_num) ) then

          unique_num = unique_num + 1

          a1(unique_num) = a1(itest)
          a2(unique_num) = a2(itest)

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
      subroutine lvec_print ( n, a, title )

c*********************************************************************72
c
cc LVEC_PRINT prints an LVEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, logical A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      logical a(n)
      character * ( * ) title

      call lvec_print_some ( n, a, 1, n, title )

      return
      end
      subroutine lvec_print_some ( n, a, i_lo, i_hi, title )

c*********************************************************************72
c
cc LVEC_PRINT_SOME prints "some" of an LVEC.
c
c  Discussion:
c
c    An LVEC is a vector of logical values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, logical A(N), the vector to be printed.
c
c    Input, integer I_LO, I_HI, the first and last indices to 
c    print. The routine expects 1 <= I_LO <= I_HI <= N.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      logical a(n)
      integer i
      integer i_hi
      integer i_lo
      character * ( * ) title

      if ( 0 .lt. len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
      end if

      write ( *, '(a)' ) ' '
      do i = max ( i_lo, 1 ), min ( i_hi, n )
        write ( *, '(2x,i8,2x,l1)' ) i, a(i)
      end do

      return
      end
      subroutine node_merge ( dim_num, node_num, node_xy, tolerance, 
     &  node_rep )

c*********************************************************************72
c
cc NODE_MERGE detects nodes that should be merged.
c
c  Discussion:
c
c    Two nodes "should" be merged if they are within TOLERANCE distance
c    of each other.
c
c    With a tolerance of 0, only exactly equal nodes are counted.
c
c    With a positive tolerance, a pair of nodes inside a circle of
c    radius TOLERANCE result in a count of 1 duplicate.
c
c    However, what do we do if nodes A, B and C are arranged in a line,c
c    with A and B just within TOLERANCE of each other, and B and C just
c    within tolerance of each other?  What we do here is make a choice
c    that can be defended consistently.  A and B define an equivalence
c    class because they are closer than TOLERANCE.  C is then added to
c    this equivalence class, because it is within TOLERANCE of at least
c    on thing in that equivalence class.
c
c    Thus, if 100 nodes are separated pairwise by slightly less
c    than TOLERANCE, a total of 99 duplicates will be counted.
c
c    The program starts out by giving each node its own label.
c    If it finds that two nodes should be merged, then the index of
c    one node is used as the label for both.  This process continues
c    until all nodes have been considered.  The number of unique nodes
c    is the number of unique values in the output quantity NODE_REP.
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
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(DIM_NUM,NODE_NUM), the nodes.
c
c    Input, double precision TOLERANCE, the maximum distance between
c    two nodes regarded as duplicate.
c
c    Output, integer NODE_REP(NODE_NUM), the "representative" of 
c    each node.  NODE_REP(NODE) is the index of a node which is within 
c    TOLERANCE of node NODE, or for which a chain of nodes can be found, all
c    having the same representative, and all of which are pairwise closer 
c    than TOLERANCE.
c
      implicit none

      integer dim_num
      integer node_num

      integer dim
      double precision dist
      integer node_rep(node_num)
      double precision node_xy(dim_num,node_num)
      integer node1
      integer node2
      double precision r8_huge
      integer rep
      double precision rep_dist(node_num)
      double precision tolerance

      do node1 = 1, node_num
        node_rep(node1) = node1
      end do

      do node1 = 1, node_num

        do node2 = 1, node_num
          rep_dist(node2) = r8_huge ( )
        end do

        do node2 = 1, node_num

          dist = 0.0D+00
          do dim = 1, dim_num
            dist = dist + ( node_xy(dim,node1) - node_xy(dim,node2) )**2
          end do
          dist = sqrt ( dist )

          rep = node_rep(node2)

          if ( dist .lt. rep_dist(rep) ) then
            rep_dist(rep) = dist
          end if

        end do

        do node2 = 1, node_num
          rep = node_rep(node2)
          if ( rep_dist(rep) .le. tolerance ) then
            node_rep(node2) = node1
          end if
        end do

      end do

      return
      end
      subroutine ns_adj_col_set ( node_num, triangle_num, variable_num, 
     &  triangle_node, triangle_neighbor, node_u_variable, 
     &  node_v_variable, node_p_variable, adj_num, adj_col )

c*********************************************************************72
c
cc NS_ADJ_COL_SET sets the COL array in a Navier Stokes triangulation.
c
c  Discussion:
c
c    This routine also returns the value of ADJ_NUM, the number of
c    Navier-Stokes variable adjacencies, which should be identical to the 
c    value that would have been computed by calling NS_ADJ_COUNT.
c
c    This routine is called to set up the ADJ_COL array, which indicates
c    the number of entries needed to store each column in the sparse
c    compressed column storage to be used for the adjacency matrix.
c
c    The triangulation is assumed to involve 6-node triangles.
c
c    Variables for the horizontal and vertical velocities are associated
c    with every node.  Variables for the pressure are associated only with
c    the vertex nodes.
c
c    We are interested in determining the number of nonzero entries in the
c    stiffness matrix of the Stokes equations, or the jacobian matrix of
c    the Navier Stokes equations.  To this end, we will say, somewhat
c    too broadly, that two variables are "adjacent" if their associated 
c    nodes both occur in some common element.  This adjacency of variables
c    I and J is taken to be equivalent to the possible nonzeroness of
c    matrix entries A(I,J) and A(J,I).
c
c    A sparse compressed column format is used to store the counts for
c    the nonzeroes.  In other words, while the value ADJ_NUM reports the
c    number of adjacencies, the vector ADJ_COL is sufficient to allow us
c    to properly set up a sparse compressed matrix for the actual storage
c    of the sparse matrix, if we desire to proceed.
c
c  Local Node Numbering:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  6   5  side 2
c       |    \
c    3  |     \
c       |      \
c       1---4---2
c
c         side 1
c
c  Variable Diagram:
c
c      UVP
c       |\
c       | \
c       |  \
c      UV   UV
c       |    \
c       |     \
c       |      \
c      UVP--UV--UVP
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer VARIABLE_NUM, the number of variables.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer NODE_U_VARIABLE(NODE_NUM), 
c    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the 
c    horizontal velocity, vertical velocity and pressure variables associated 
c    with a node, or -1 if no such variable is associated with the node.
c
c    Output, integer ADJ_NUM, the number of 
c    Navier Stokes variable adjacencies.
c
c    Output, integer ADJ_COL(VARIABLE_NUM+1).  Information about 
c    variable J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )
      integer variable_num

      integer adj_num
      integer adj_col(variable_num+1)
      integer i
      integer n1
      integer n2
      integer n3
      integer n4
      integer n5
      integer n6
      integer node
      integer node_p_variable(node_num)
      integer node_u_variable(node_num)
      integer node_v_variable(node_num)
      integer p1
      integer p2
      integer p3
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)
      integer u1
      integer u2
      integer u3
      integer u4
      integer u5
      integer u6
      integer v1
      integer v2
      integer v3
      integer v4
      integer v5
      integer v6
      integer variable

      adj_num = 0
c
c  Set every variable to be adjacent to itself.
c
      do i = 1, variable_num
        adj_col(i) = 1
      end do
c
c  Set every variable to be adjacent to the other variables associated with
c  that node. 
c
c  U <=> V
c  U <=> P (if there is a P variable)
c  V <=> P (if there is a P variable)
c
      do node = 1, node_num

        u1 = node_u_variable(node)
        v1 = node_v_variable(node)
        p1 = node_p_variable(node)

        adj_col(u1) = adj_col(u1) + 1
        adj_col(v1) = adj_col(v1) + 1  

        if ( 0 .lt. p1 ) then
          adj_col(u1) = adj_col(u1) + 1
          adj_col(v1) = adj_col(v1) + 1
          adj_col(p1) = adj_col(p1) + 2
        end if

      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
        n4 = triangle_node(4,triangle)
        n5 = triangle_node(5,triangle)
        n6 = triangle_node(6,triangle)

        u1 = node_u_variable(n1)
        v1 = node_v_variable(n1)
        p1 = node_p_variable(n1)

        u2 = node_u_variable(n2)
        v2 = node_v_variable(n2)
        p2 = node_p_variable(n2)

        u3 = node_u_variable(n3)
        v3 = node_v_variable(n3)
        p3 = node_p_variable(n3)

        u4 = node_u_variable(n4)
        v4 = node_v_variable(n4)

        u5 = node_u_variable(n5)
        v5 = node_v_variable(n5)

        u6 = node_u_variable(n6)
        v6 = node_v_variable(n6)
c
c  For sure, we add the new adjacencies:
c
c    U5 V5 <=> U1 V1 P1
c    U6 V6 <=> U2 V2 P2
c    U4 V4 <=> U3 V3 P3
c    U5 V5 <=> U4 V4
c    U6 V6 <=> U4 V4
c    U6 V6 <=> U5 V5
c
        adj_col(u1) = adj_col(u1) + 2
        adj_col(v1) = adj_col(v1) + 2
        adj_col(p1) = adj_col(p1) + 2

        adj_col(u2) = adj_col(u2) + 2
        adj_col(v2) = adj_col(v2) + 2
        adj_col(p2) = adj_col(p2) + 2

        adj_col(u3) = adj_col(u3) + 2
        adj_col(v3) = adj_col(v3) + 2
        adj_col(p3) = adj_col(p3) + 2

        adj_col(u4) = adj_col(u4) + 7
        adj_col(v4) = adj_col(v4) + 7

        adj_col(u5) = adj_col(u5) + 7
        adj_col(v5) = adj_col(v5) + 7

        adj_col(u6) = adj_col(u6) + 7
        adj_col(v6) = adj_col(v6) + 7
c
c  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
c  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
c  Maybe add
c
c    U1 V1 P1 <=> U2 V2 P2
c    U1 V1 P1 <=> U4 V4
c    U2 V2 P2 <=> U4 V4
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          adj_col(u1) = adj_col(u1) + 5
          adj_col(v1) = adj_col(v1) + 5
          adj_col(p1) = adj_col(p1) + 5

          adj_col(u2) = adj_col(u2) + 5
          adj_col(v2) = adj_col(v2) + 5
          adj_col(p2) = adj_col(p2) + 5

          adj_col(u4) = adj_col(u4) + 6
          adj_col(v4) = adj_col(v4) + 6

        end if
c
c  Maybe add
c
c    U2 V2 P2 <=> U3 V3 P3
c    U2 V2 P2 <=> U5 V5
c    U3 V3 P3 <=> U5 V5
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          adj_col(u2) = adj_col(u2) + 5
          adj_col(v2) = adj_col(v2) + 5
          adj_col(p2) = adj_col(p2) + 5

          adj_col(u3) = adj_col(u3) + 5
          adj_col(v3) = adj_col(v3) + 5
          adj_col(p3) = adj_col(p3) + 5

          adj_col(u5) = adj_col(u5) + 6
          adj_col(v5) = adj_col(v5) + 6

        end if
c
c  Maybe add
c
c    U1 V1 P1 <=> U3 V3 P3
c    U1 V1 P1 <=> U6 V6
c    U3 V3 P3 <=> U6 V6
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          adj_col(u1) = adj_col(u1) + 5
          adj_col(v1) = adj_col(v1) + 5
          adj_col(p1) = adj_col(p1) + 5

          adj_col(u3) = adj_col(u3) + 5
          adj_col(v3) = adj_col(v3) + 5
          adj_col(p3) = adj_col(p3) + 5

          adj_col(u6) = adj_col(u6) + 6
          adj_col(v6) = adj_col(v6) + 6

        end if
          
      end do
c
c  We used ADJ_COL to count the number of entries in each column.
c  Convert it to pointers into the ADJ array.
c
      do i = 1, variable_num
        adj_col(i+1) = adj_col(i)
      end do
      adj_col(1) = 1
      do variable = 2, variable_num + 1
        adj_col(variable) = adj_col(variable-1) + adj_col(variable)
      end do

      adj_num = adj_col(variable_num+1) - 1

      return
      end
      subroutine ns_adj_count ( node_num, triangle_num, variable_num, 
     &  triangle_node, triangle_neighbor, node_u_variable, 
     &  node_v_variable, node_p_variable, adj_num )

c*********************************************************************72
c
cc NS_ADJ_COUNT counts adjacencies in a Navier Stokes triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The value of ADJ_NUM computed and returned by this routine should
c    be identical to the value computed by NS_ADJ_COL_SET.
c
c    The triangulation is assumed to involve 6-node triangles.
c
c    Variables for the horizontal and vertical velocities are associated
c    with every node.  Variables for the pressure are associated only with
c    the vertex nodes.
c
c    We are interested in determining the number of nonzero entries in the
c    stiffness matrix of the Stokes equations, or the jacobian matrix of
c    the Navier Stokes equations.  To this end, we will say, somewhat
c    too broadly, that two variables are "adjacent" if their associated 
c    nodes both occur in some common element.  This adjacency of variables
c    I and J is taken to be equivalent to the possible nonzeroness of
c    matrix entries A(I,J) and A(J,I).
c
c    A sparse compressed column format is used to store the counts for
c    the nonzeroes.  In other words, while the value ADJ_NUM reports the
c    number of adjacencies, the vector ADJ_COL is sufficient to allow us
c    to properly set up a sparse compressed matrix for the actual storage
c    of the sparse matrix, if we desire to proceed.
c
c  Local Node Numbering:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  6   5  side 2
c       |    \
c    3  |     \
c       |      \
c       1---4---2
c
c         side 1
c
c  Variable Diagram:
c
c      UVP
c       |\
c       | \
c       |  \
c      UV   UV
c       |    \
c       |     \
c       |      \
c      UVP--UV--UVP
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer VARIABLE_NUM, the number of variables.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer NODE_U_VARIABLE(NODE_NUM), 
c    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the 
c    horizontal velocity, vertical velocity and pressure variables associated 
c    with a node, or -1 if no such variable is associated with the node.
c
c    Output, integer ADJ_NUM, the value of ADJ_NUM, the number of 
c    Navier Stokes variable adjacencies.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )
      integer variable_num

      integer adj_num
      integer node
      integer node_p_variable(node_num)
      integer node_u_variable(node_num)
      integer node_v_variable(node_num)
      integer p1
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)
      integer variable

      adj_num = 0
c
c  Set every variable to be adjacent to itself.
c
      adj_num = variable_num
c
c  Set every variable to be adjacent to the other variables associated with
c  that node. 
c
c  U <=> V
c  U <=> P (if there is a P variable)
c  V <=> P (if there is a P variable)
c
      do node = 1, node_num

        adj_num = adj_num + 2

        p1 = node_p_variable(node)

        if ( 0 .lt. p1 ) then
          adj_num = adj_num + 4
        end if

      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num
c
c  For sure, we add the new adjacencies:
c
c    U5 V5 <=> U1 V1 P1
c    U6 V6 <=> U2 V2 P2
c    U4 V4 <=> U3 V3 P3
c    U5 V5 <=> U4 V4
c    U6 V6 <=> U4 V4
c    U6 V6 <=> U5 V5
c
        adj_num = adj_num + 60
c
c  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
c  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
c  Maybe add
c
c    U1 V1 P1 <=> U2 V2 P2
c    U1 V1 P1 <=> U4 V4
c    U2 V2 P2 <=> U4 V4
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_num = adj_num + 42
        end if
c
c  Maybe add
c
c    U2 V2 P2 <=> U3 V3 P3
c    U2 V2 P2 <=> U5 V5
c    U3 V3 P3 <=> U5 V5
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_num = adj_num + 42
        end if
c
c  Maybe add
c
c    U1 V1 P1 <=> U3 V3 P3
c    U1 V1 P1 <=> U6 V6
c    U3 V3 P3 <=> U6 V6
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_num = adj_num + 42
        end if
          
      end do

      return
      end
      subroutine ns_adj_insert ( v1, v2, variable_num, adj_num, 
     &  adj_col_free, adj_row )

c*********************************************************************72
c
cc NS_ADJ_INSERT inserts an adjacency into a compressed column adjacency matrix.
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
c  Parameters:
c
c    Input, integer V1, V2, the indices of two items which are 
c    adjacent.
c
c    Input, integer VARIABLE_NUM, the number of items.
c
c    Input, integer ADJ_NUM, the number of entries available 
c    in ADJ_ROW.
c
c    Input/output, integer ADJ_COL_FREE(VARIABLE_NUM), the next 
c    free location in which an entry for a given column can be stored.  On 
c    output, two pointers have been updated.
c
c    Input/output, integer ADJ_ROW(ADJ_NUM), the row indices of 
c    the Navier Stokes variable adjacency matrix.  On output, two new entries 
c    have been added.
c
      integer adj_num
      integer variable_num
        
      integer adj_col_free(variable_num)
      integer adj_row(adj_num)
      integer v1
      integer v2

      adj_row(adj_col_free(v1)) = v2
      adj_col_free(v1) = adj_col_free(v1) + 1

      if ( v1 .eq. v2 ) then
        return
      end if

      adj_row(adj_col_free(v2)) = v1
      adj_col_free(v2) = adj_col_free(v2) + 1

      return
      end
      subroutine ns_adj_row_set ( node_num, triangle_num, variable_num, 
     &  triangle_node, triangle_neighbor, node_u_variable, 
     &  node_v_variable, node_p_variable, adj_num, adj_col, adj_row )

c*********************************************************************72
c
cc NS_ADJ_ROW_SET sets the Navier Stokes sparse compressed column row indices.
c
c  Discussion:
c
c    After NS_ADJ_COUNT has been called to count ADJ_NUM, the number of
c    variable adjacencies and to set up ADJ_COL, the compressed column pointer, 
c    this routine can be called to assign values to ADJ_ROW, the row
c    indices of the sparse compressed column adjacency matrix.
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer VARIABLE_NUM, the number of variables.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer NODE_U_VARIABLE(NODE_NUM), 
c    NODE_V_VARIABLE(NODE_NUM), NODE_P_VARIABLE(NODE_NUM), the index of the 
c    horizontal velocity, vertical velocity and pressure variables associated 
c    with a node, or -1 if no such variable is associated with the node.
c
c    Input, integer ADJ_NUM, the number of Navier Stokes variable
c    adjacencies.
c
c    Input, integer ADJ_COL(VARIABLE_NUM+1).  Information about
c    variable J  is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
c    Output, integer ADJ_ROW(ADJ_NUM), the row indices of the 
c    Navier Stokes variable adjacency matrix.
c
c  Local Parameters:
c
c    Local, integer ADJ_COL_FREE(VARIABLE_NUM), for each column,
c    the location in ADJ_ROW which can store the next row index.
c    
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )
      integer variable_num

      integer adj_num
      integer adj_col(variable_num+1)
      integer adj_col_free(variable_num)
      integer adj_row(adj_num)
      integer i
      integer k1
      integer k2
      integer n1
      integer n2
      integer n3
      integer n4
      integer n5
      integer n6
      integer node
      integer node_p_variable(node_num)
      integer node_u_variable(node_num)
      integer node_v_variable(node_num)
      integer number
      integer p1
      integer p2
      integer p3
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)
      integer u1
      integer u2
      integer u3
      integer u4
      integer u5
      integer u6
      integer v
      integer v1
      integer v2
      integer v3
      integer v4
      integer v5
      integer v6

      do i = 1, adj_num
        adj_row(i) = -1
      end do

      do v = 1, variable_num
        adj_col_free(v) = adj_col(v)
      end do
c
c  Set every variable to be adjacent to itself.
c
      do v = 1, variable_num
        call ns_adj_insert ( v, v, variable_num, adj_num, adj_col_free, 
     &    adj_row )
      end do
c
c  Set every variable to be adjacent to the other variables associated with
c  that node. 
c
c  U <=> V
c  U <=> P (if there is a P variable)
c  V <=> P (if there is a P variable)
c
      do node = 1, node_num

        u1 = node_u_variable(node)
        v1 = node_v_variable(node)
        p1 = node_p_variable(node)

        call ns_adj_insert ( u1, v1, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        if ( 0 .lt. p1 ) then

          call ns_adj_insert ( u1, p1, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( v1, p1, variable_num, adj_num, 
     &      adj_col_free, adj_row )

        end if

      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
        n4 = triangle_node(4,triangle)
        n5 = triangle_node(5,triangle)
        n6 = triangle_node(6,triangle)

        u1 = node_u_variable(n1)
        v1 = node_v_variable(n1)
        p1 = node_p_variable(n1)

        u2 = node_u_variable(n2)
        v2 = node_v_variable(n2)
        p2 = node_p_variable(n2)

        u3 = node_u_variable(n3)
        v3 = node_v_variable(n3)
        p3 = node_p_variable(n3)

        u4 = node_u_variable(n4)
        v4 = node_v_variable(n4)

        u5 = node_u_variable(n5)
        v5 = node_v_variable(n5)

        u6 = node_u_variable(n6)
        v6 = node_v_variable(n6)
c
c  For sure, we add the new adjacencies:
c
c    U5 V5 <=> U1 V1 P1
c    U6 V6 <=> U2 V2 P2
c    U4 V4 <=> U3 V3 P3
c    U5 V5 <=> U4 V4
c    U6 V6 <=> U4 V4
c    U6 V6 <=> U5 V5
c
        call ns_adj_insert ( u5, u1, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u5, v1, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u5, p1, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v5, u1, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v5, v1, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v5, p1, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        call ns_adj_insert ( u6, u2, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u6, v2, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u6, p2, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, u2, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, v2, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, p2, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        call ns_adj_insert ( u4, u3, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u4, v3, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u4, p3, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v4, u3, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v4, v3, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v4, p3, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        call ns_adj_insert ( u5, u4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u5, v4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v5, u4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v5, v4, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        call ns_adj_insert ( u6, u4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u6, v4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, u4, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, v4, variable_num, adj_num, 
     &    adj_col_free, adj_row )

        call ns_adj_insert ( u6, u5, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( u6, v5, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, u5, variable_num, adj_num, 
     &    adj_col_free, adj_row )
        call ns_adj_insert ( v6, v5, variable_num, adj_num, 
     &    adj_col_free, adj_row )
c
c  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
c  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
c  Maybe add
c
c    U1 V1 P1 <=> U2 V2 P2
c    U1 V1 P1 <=> U4 V4
c    U2 V2 P2 <=> U4 V4
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          call ns_adj_insert ( u1, u2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, v2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, p2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, u2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, v2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, p2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, u2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, v2, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, p2, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u1, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u2, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u2, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, u4, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, v4, variable_num, adj_num, 
     &      adj_col_free, adj_row )

        end if
c
c  Maybe add
c
c    U2 V2 P2 <=> U3 V3 P3
c    U2 V2 P2 <=> U5 V5
c    U3 V3 P3 <=> U5 V5
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          call ns_adj_insert ( u2, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u2, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u2, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u2, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u2, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v2, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p2, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u3, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u3, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v3, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v3, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p3, u5, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p3, v5, variable_num, adj_num, 
     &      adj_col_free, adj_row )

        end if
c
c  Maybe add
c
c    U1 V1 P1 <=> U3 V3 P3
c    U1 V1 P1 <=> U6 V6
c    U3 V3 P3 <=> U6 V6
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          call ns_adj_insert ( u1, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, u3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, v3, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, p3, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u1, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u1, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v1, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p1, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )

          call ns_adj_insert ( u3, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( u3, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v3, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( v3, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p3, u6, variable_num, adj_num, 
     &      adj_col_free, adj_row )
          call ns_adj_insert ( p3, v6, variable_num, adj_num, 
     &      adj_col_free, adj_row )

        end if
          
      end do
c
c  Ascending sort the entries for each variable.
c
      do v = 1, variable_num
        k1 = adj_col(v)
        k2 = adj_col(v+1)-1
        number = k2 + 1 - k1
        call i4vec_sort_heap_a ( number, adj_row(k1:k2) )
      end do

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
      subroutine points_delaunay_naive_2d ( node_num, node_xy, maxtri, 
     &  triangle_num, triangle_node )

!*********************************************************************72
!
!! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
!
!  Discussion:
!
!    This routine is only suitable as a demonstration code for small
!    problems.  Its running time is of order NODE_NUM**4.  Much faster
!    algorithms are available.
!
!    Given a set of nodes in the plane, a triangulation is set of
!    triples of distinct nodes, forming triangles, so that every
!    point within the convex hull of the set of nodes is either
!    one of the nodes, or lies on an edge of one or more triangles,
!    or lies within exactly one triangle.
!
!    A Delaunay triangulation is a triangulation with additional
!    properties.
!
!    NODE_NUM must be at least 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.  
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer MAXTRI, the maximum number of triangles.
!
!    Output, integer TRIANGLE_NUM, the number of triangles in 
!    the triangulation.
!
!    Output, integer TRIANGLE_NODE(3,MAXTRI), the indices of 
!    the triangle nodes.
!
      implicit none

      integer, parameter :: dim_num = 2
      integer maxtri
      integer node_num

      logical flag
      integer i
      integer j
      integer k
      integer m
      double precision node_xy(dim_num,node_num)
      integer triangle_node(3,maxtri)
      integer triangle_num
      double precision xn
      double precision yn
      double precision z(node_num)
      double precision zn

      triangle_num = 0

      if ( node_num .lt. 3 ) then
        return
      end if
!
!  Compute Z = X*X + Y*Y.
!
      do i = 1, node_num
        z(i) = node_xy(1,i)**2 + node_xy(2,i)**2
      end do
!  
!  For each triple (I,J,K):
!
      do i = 1, node_num - 2
        do j = i + 1, node_num
          do k = i + 1, node_num

            if ( j .ne. k ) then

              xn = ( node_xy(2,j) - node_xy(2,i) ) * ( z(k) - z(i) ) 
     &           - ( node_xy(2,k) - node_xy(2,i) ) * ( z(j) - z(i) )

              yn = ( node_xy(1,k) - node_xy(1,i) ) * ( z(j) - z(i) ) 
     &           - ( node_xy(1,j) - node_xy(1,i) ) * ( z(k) - z(i) )

              zn = ( node_xy(1,j) - node_xy(1,i) ) 
     &           * ( node_xy(2,k) - node_xy(2,i) ) 
     &           - ( node_xy(1,k) - node_xy(1,i) ) 
     &           * ( node_xy(2,j) - node_xy(2,i) )

              flag = ( zn .lt. 0.0D+00 )

              if ( flag ) then
                do m = 1, node_num
                  flag = flag .and. 
     &              ( ( node_xy(1,m) - node_xy(1,i) ) * xn 
     &              + ( node_xy(2,m) - node_xy(2,i) ) * yn 
     &              + ( z(m)   - z(i) )   * zn .le. 0.0D+00 )
                end do
              end if

              if ( flag ) then
                if ( triangle_num .lt. maxtri ) then
                  triangle_num = triangle_num + 1
                  triangle_node(1,triangle_num) = i
                  triangle_node(2,triangle_num) = j
                  triangle_node(3,triangle_num) = k
                end if
              end if

            end if

          end do
        end do
      end do

      return
      end
      subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

c*********************************************************************72
c
cc POINTS_HULL_2D computes the convex hull of 2D points.
c
c  Discussion:
c
c    The work involved is N*log(H), where N is the number of points, and H is
c    the number of points that are on the hull.
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Output, integer HULL_NUM, the number of nodes that lie on 
c    the convex hull.
c
c    Output, integer HULL(NODE_NUM).  Entries 1 through HULL_NUM 
c    contain the indices of the nodes that form the convex hull, in order.
c
      implicit none

      integer node_num
      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle
      double precision angle_max
      double precision angle_rad_2d
      double precision di
      double precision dr
      integer first
      integer hull(node_num)
      integer hull_num
      integer i
      double precision node_xy(dim_num,node_num)
      double precision p_xy(dim_num)
      integer q
      double precision q_xy(2)
      integer r
      double precision r_xy(2)

      if ( node_num .lt. 1 ) then
        hull_num = 0
        return
      end if
c
c  If NODE_NUM = 1, the hull is the point.
c
      if ( node_num .eq. 1 ) then
        hull_num = 1
        hull(1) = 1
        return
      end if
c
c  If NODE_NUM = 2, then the convex hull is either the two distinct points,
c  or possibly a single (repeated) point.
c
      if ( node_num .eq. 2 ) then

        if ( node_xy(1,1) .ne. node_xy(1,2) .or. 
     &       node_xy(2,1) .ne. node_xy(2,2) ) then
          hull_num = 2
          hull(1) = 1
          hull(2) = 2
        else
          hull_num = 1
          hull(1) = 1
        end if

        return

      end if
c
c  Find the leftmost point, and take the bottom-most in a tie.
c  Call it "Q".
c
      q = 1
      do i = 2, node_num
        if ( node_xy(1,i) .lt. node_xy(1,q) .or. 
     &    ( node_xy(1,i) .eq. node_xy(1,q) .and. 
     &      node_xy(2,i) .lt. node_xy(2,q) ) ) then
          q = i
        end if
      end do

      q_xy(1) = node_xy(1,q)
      q_xy(2) = node_xy(2,q)
c
c  Remember the starting point, so we know when to stop!
c
      first = q
      hull_num = 1
      hull(1) = q
c
c  For the first point, make a dummy previous point, 1 unit south,
c  and call it "P".
c
      p_xy(1) = q_xy(1)
      p_xy(2) = q_xy(2) - 1.0D+00
c
c  Now, having old point P, and current point Q, find the new point R
c  so the angle PQR is maximal.
c
c  Watch out for the possibility that the two nodes are identical.
c
10    continue

        r = 0
        angle_max = 0.0D+00

        do i = 1, node_num

          if ( i .ne. q .and.
     &      ( node_xy(1,i) .ne. q_xy(1) .or.
     &        node_xy(2,i) .ne. q_xy(2) ) ) then

            angle = angle_rad_2d ( p_xy, q_xy, node_xy(1,i) )

            if ( r .eq. 0 .or. angle_max .lt. angle ) then

              r = i
              r_xy(1) = node_xy(1,r)
              r_xy(2) = node_xy(2,r)
              angle_max = angle
c
c  In case of ties, choose the nearer point.
c
            else if ( r .ne. 0 .and. angle .eq. angle_max ) then

              di = ( node_xy(1,i) - q_xy(1) )**2 
     &           + ( node_xy(2,i) - q_xy(2) )**2

              dr = ( r_xy(1)      - q_xy(1) )**2 
     &           + ( r_xy(2)      - q_xy(2) )**2

              if ( di .lt. dr ) then
                r = i
                r_xy(1) = node_xy(1,r)
                r_xy(2) = node_xy(2,r)
                angle_max = angle
              end if

            end if

          end if

        end do
c
c  If we've returned to our starting point, exit.
c
        if ( r .eq. first ) then
          go to 20
        end if

        hull_num = hull_num + 1

        if ( node_num .lt. hull_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
          write ( *, '(a)' ) '  The algorithm failed.'
          stop
        end if
c
c  Add point R to the convex hull.
c
        hull(hull_num) = r
c
c  Set P = Q, Q = R, and prepare to search for the next point R.
c
        q = r

        p_xy(1) = q_xy(1)
        p_xy(2) = q_xy(2)

        q_xy(1) = r_xy(1)
        q_xy(2) = r_xy(2)

      go to 10

20    continue

      return
      end
      subroutine points_point_near_naive_nd ( dim_num, nset, pset, p, 
     &  i_min, d_min )

!*********************************************************************72
!
!! POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, double precision PSET(DIM_NUM,NSET), the points in the set.
!
!    Input, double precision P(DIM_NUM), the point whose nearest neighbor
!    is sought.
!
!    Output, integer I_MIN, the index of the nearest point in 
!    PSET to P.
!
!    Output, double precision D_MIN, the distance between P(*) 
!    and PSET(*,I_MIN).
!
      implicit none

      integer dim_num
      integer nset

      double precision d
      double precision d_min
      integer dim
      integer i
      integer i_min
      double precision p(dim_num)
      double precision pset(dim_num,nset)
      double precision r8_huge

      d_min = r8_huge ( )
      i_min = -1

      do i = 1, nset
        d = 0.0D+00
        do dim = 1, dim_num
          d = d + ( p(dim) - pset(dim,i) )**2
        end do
        if ( d .lt. d_min ) then
          d_min = d
          i_min = i
        end if
      end do

      d_min = sqrt ( d_min )

      return
      end
      subroutine q_measure ( n, z, triangle_order, triangle_num, 
     &  triangle_node, q_min, q_max, q_ave, q_area )

c*********************************************************************72
c
cc Q_MEASURE determines the triangulated pointset quality measure Q.
c
c  Discussion:
c
c    The Q measure evaluates the uniformity of the shapes of the triangles
c    defined by a triangulated pointset.
c
c    For a single triangle T, the value of Q(T) is defined as follows:
c
c      TAU_IN = radius of the inscribed circle,
c      TAU_OUT = radius of the circumscribed circle,
c
c      Q(T) = 2 * TAU_IN / TAU_OUT 
c        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
c
c    where A, B and C are the lengths of the sides of the triangle T.
c
c    The Q measure computes the value of Q(T) for every triangle T in the
c    triangulation, and then computes the minimum of this
c    set of values:
c
c      Q_MEASURE = min ( all T in triangulation ) Q(T)
c
c    In an ideally regular mesh, all triangles would have the same
c    equilateral shape, for which Q = 1.  A good mesh would have
c    0.5 < Q.
c
c    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
c    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
c    triangles.  Generally, a maximal triangulation is expected, namely,
c    a triangulation whose image is a planar graph, but for which the
c    addition of any new triangle would mean the graph was no longer planar.
c    A Delaunay triangulation is a maximal triangulation which maximizes
c    the minimum angle that occurs in any triangle.
c
c    The code has been modified to 'allow' 6-node triangulations.
c    However, no effort is made to actually process the midside nodes.
c    Only information from the vertices is used.
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
c  Reference:
c
c    Max Gunzburger, John Burkardt,
c    Uniformity Measures for Point Samples in Hypercubes.
c
c    Per-Olof Persson, Gilbert Strang,
c    A Simple Mesh Generator in MATLAB,
c    SIAM Review,
c    Volume 46, Number 2, June 2004, pages 329-345.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision Z(2,N), the points.
c
c    Input, integer TRIANGLE_ORDER, the order of the triangles.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
c    the triangulation.
c
c    Output, double precision Q_MIN, Q_MAX, the minimum and maximum values
c    of Q over all triangles.
c
c    Output, double precision Q_AVE, the average value of Q.
c
c    Output, double precision Q_AREA, the average value of Q, weighted by
c    the area of each triangle.
c
      implicit none

      integer n
      integer triangle_num
      integer triangle_order

      integer a_index
      double precision ab_length
      double precision area
      double precision area_total
      integer b_index
      double precision bc_length
      integer c_index
      double precision ca_length
      double precision q
      double precision q_area
      double precision q_ave
      double precision q_max
      double precision q_min
      double precision r8_huge
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      double precision x1
      double precision x2
      double precision x3
      double precision y1
      double precision y2
      double precision y3
      double precision z(2,n)

      q_min = r8_huge ( )
      q_max = - huge ( q_max )
      q_ave = 0.0D+00
      q_area = 0.0D+00
      area_total = 0.0D+00

      do triangle = 1, triangle_num

        a_index = triangle_node(1,triangle)
        b_index = triangle_node(2,triangle)
        c_index = triangle_node(3,triangle)

        ab_length = sqrt ( 
     &      ( z(1,a_index) - z(1,b_index) )**2 
     &    + ( z(2,a_index) - z(2,b_index) )**2 )

        bc_length = sqrt ( 
     &      ( z(1,b_index) - z(1,c_index) )**2 
     &    + ( z(2,b_index) - z(2,c_index) )**2 )

        ca_length = sqrt ( 
     &      ( z(1,c_index) - z(1,a_index) )**2 
     &    + ( z(2,c_index) - z(2,a_index) )**2 )

        q = ( bc_length + ca_length - ab_length ) 
     &    * ( ca_length + ab_length - bc_length ) 
     &    * ( ab_length + bc_length - ca_length ) 
     &    / ( ab_length * bc_length * ca_length )

        x1 = z(1,triangle_node(1,triangle))
        y1 = z(2,triangle_node(1,triangle))
        x2 = z(1,triangle_node(2,triangle))
        y2 = z(2,triangle_node(2,triangle))
        x3 = z(1,triangle_node(3,triangle))
        y3 = z(2,triangle_node(3,triangle))

        area = 0.5D+00 * abs ( x1 * ( y2 - y3 ) 
     &                       + x2 * ( y3 - y1 ) 
     &                       + x3 * ( y1 - y2 ) )
    
        q_min = min ( q_min, q )
        q_max = max ( q_max, q )
        q_ave = q_ave + q
        q_area = q_area + q * area

        area_total = area_total + area

      end do

      q_ave = q_ave / dble ( triangle_num )
      q_area = q_area / area_total

      return
      end
      subroutine quad_convex_random ( seed, xy )

c*********************************************************************72
c
cc QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
c
c  Description:
c
c    The quadrilateral is constrained in that the vertices must all lie
c    with the unit square.
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
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision XY(2,NODE_NUM), the coordinates of the 
c    nodes of the quadrilateral, given in counterclockwise order.
c
      implicit none

      integer node_num
      parameter ( node_num = 4 )

      integer i
      integer ival(node_num)
      integer j
      integer nval
      integer seed
      double precision xy(2,node_num)
      double precision xy_random(2,node_num)

10    continue
c
c  Generate 4 random points.
c
        call r8mat_uniform_01 ( 2, node_num, seed, xy_random )
c
c  Determine the convex hull.
c
        call points_hull_2d ( node_num, xy_random, nval, ival )
c
c  If NVAL .lt. NODE_NUM, then our convex hull is a triangle.
c  Try again.
c
        if ( nval .eq. node_num ) then
          go to 20
        end if

      go to 10
c
c  Make an ordered copy of the random points.
c
20    continue

      do j = 1, nval
        do i = 1, 2
          xy(i,j) = xy_random(i,ival(j))
        end do
      end do

      return
      end
      function r8_acos ( c )

c*********************************************************************72
c
cc R8_ACOS computes the arc cosine function, with argument truncation.
c
c  Discussion:
c
c    If you call your system ACOS routine with an input argument that is
c    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
c    surprise (I did).
c
c    This routine simply truncates arguments outside the range.
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
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision C, the argument.
c
c    Output, double precision R8_ACOS, an angle whose cosine is C.
c
      implicit none

      double precision c
      double precision c2
      double precision r8_acos

      c2 = c
      c2 = max ( c2, -1.0D+00 )
      c2 = min ( c2, +1.0D+00 )

      r8_acos = acos ( c2 )

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
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
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
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r82vec_permute ( n, p, a )

c*********************************************************************72
c
cc R82VEC_PERMUTE permutes an R82VEC in place.
c
c  Discussion:
c
c    An R82VEC is an array of pairs of R8 values.
c
c    The same logic can be used to permute an array of objects of any 
c    arithmetic type, or an array of objects of any complexity.  The only
c    temporary storage required is enough to store a single object.  The number
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
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer P(N), the permutation.  P(I) = J means
c    that the I-th element of the output array should be the J-th
c    element of the input array.  
c
c    Input/output, double precision A(2,N), the array to be permuted.
c
      implicit none

      integer n
      integer, parameter :: dim_num = 2

      double precision a(dim_num,n)
      double precision a_temp(dim_num)
      integer base
      parameter ( base = 1 )
      integer dim
      integer ierror
      integer iget
      integer iput
      integer istart
      integer p(n)

      call perm_check2 ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if
c
c  Search for the next element of the permutation that has not been used.
c
      do istart = 1, n

        if ( p(istart) .lt. 0 ) then

        else if ( p(istart) .eq. istart ) then

          p(istart) = - p(istart)

        else

          do dim = 1, dim_num
            a_temp(dim) = a(dim,istart)
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
              write ( *, '(a)' ) 
     &          '  A permutation index is out of range.'
              write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
              stop
            end if

            if ( iget .eq. istart ) then
              do dim = 1, dim_num
                a(dim,iput) = a_temp(dim)
              end do
              go to 20
            end if

            do dim = 1, dim_num
              a(dim,iput) = a(dim,iget)
            end do

          go to 10

        end if

20      continue

      end do
c
c  Restore the signs of the entries.
c
      do istart = 1, n
        p(istart) = - p(istart)
      end do

      return
      end
      subroutine r82vec_sort_heap_index_a ( n, a, indx )

c*********************************************************************72
c
cc R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
c
c  Discussion:
c
c    An R82VEC is an array of R82's.
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(1:2,INDX(1:N)) is sorted,
c
c    or explicitly, by the call
c
c      call r82vec_permute ( n, indx, a )
c
c    after which A(1:2,I), I = 1 to N is sorted.
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
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(2,N), an array to be index-sorted.
c
c    Output, integer INDX(N), the sort index.  The
c    I-th element of the sorted array is A(1:2,INDX(I)).
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer n

      double precision a(dim_num,n)
      double precision aval(dim_num)
      integer dim
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer j
      integer l

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, n
        indx(i) = i
      end do

      if ( n .eq. 1 ) then
        return
      end if

      l = n / 2 + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          do dim = 1, dim_num
            aval(dim) = a(dim,indxt)
          end do

        else

          indxt = indx(ir)
          do dim = 1, dim_num
            aval(dim) = a(dim,indxt)
          end do
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then
            if (   a(1,indx(j)) .lt.  a(1,indx(j+1)) .or. 
     &           ( a(1,indx(j)) .eq. a(1,indx(j+1)) .and. 
     &             a(2,indx(j)) .lt.  a(2,indx(j+1)) ) ) then
              j = j + 1
            end if
          end if

          if (   aval(1) .lt.  a(1,indx(j)) .or. 
     &         ( aval(1) .eq. a(1,indx(j)) .and. 
     &           aval(2) .lt.  a(2,indx(j)) ) ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

      return
      end
      subroutine r8mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc R8MAT_COPY copies an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
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
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision A1(M,N), the matrix to be copied.
c
c    Output, double precision A2(M,N), a copy of the matrix.
c
      implicit none

      integer m
      integer n

      double precision a1(m,n)
      double precision a2(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
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

      return
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
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
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2004
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
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
c    03 June 2009
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

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision t(dim_num,3)

      area = 0.5D+00 * abs ( 
     &    t(1,1) * ( t(2,2) - t(2,3) ) 
     &  + t(1,2) * ( t(2,3) - t(2,1) ) 
     &  + t(1,3) * ( t(2,1) - t(2,2) ) )

      return
      end
      subroutine triangle_circumcenter_2d ( t, center )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
c
c  Discussion:
c
c    The circumcenter of a triangle is the center of the circumcircle, the
c    circle that passes through the three vertices of the triangle.
c
c    The circumcircle contains the triangle, but it is not necessarily the
c    smallest triangle to do so.
c
c    If all angles of the triangle are no greater than 90 degrees, then
c    the center of the circumscribed circle will lie inside the triangle.
c    Otherwise, the center will lie outside the triangle.
c
c    The circumcenter is the intersection of the perpendicular bisectors
c    of the sides of the triangle.
c
c    In geometry, the circumcenter of a triangle is often symbolized by "O".
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
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision CENTER(2), the circumcenter of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision asq
      double precision bot
      double precision center(dim_num)
      double precision csq
      double precision t(dim_num,3)
      double precision top(dim_num)

      asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
      csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
      
      top(1) =    ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
      top(2) =  - ( t(1,2) - t(1,1) ) * csq + ( t(1,3) - t(1,1) ) * asq

      bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) 
     &      - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

      center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot

      return
      end
      subroutine triangle_order3_physical_to_reference ( t, n, phy,
     &  ref )

c*********************************************************************72
c
cc TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE maps physical points to reference points.
c
c  Discussion:
c
c    Given the vertices of an order 3 physical triangle and a point
c    (X,Y) in the physical triangle, the routine computes the value
c    of the corresponding image point (XSI,ETA) in reference space.
c
c    This routine is also appropriate for an order 4 triangle, assuming
c    that the fourth node is always the centroid of the triangle.
c
c    This routine may be appropriate for an order 6
c    triangle, if the mapping between reference and physical space
c    is linear.  This implies, in particular, that the sides of the
c    image triangle are straight and that the "midside" nodes in the
c    physical triangle are halfway along the sides of
c    the physical triangle.
c
c  Reference Element T3:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  |  \
c    |  |   \
c    |  |    \
c    0  1-----2
c    |
c    +--0--R--1-->
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
c    Input, double precision T(2,3), the X and Y coordinates
c    of the vertices.  The vertices are assumed to be the images of
c    (0,0), (1,0) and (0,1) respectively.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision PHY(2,N), the coordinates of physical points
c    to be transformed.
c
c    Output, double precision REF(2,N), the coordinates of the corresponding
c    points in the reference space.
c
      implicit none

      integer n

      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,3)

      do j = 1, n

        ref(1,j) = ( ( t(2,3) - t(2,1) ) * ( phy(1,j) - t(1,1) )   
     &             - ( t(1,3) - t(1,1) ) * ( phy(2,j) - t(2,1) ) ) 
     &           / ( ( t(2,3) - t(2,1) ) * ( t(1,2)   - t(1,1) )   
     &             - ( t(1,3) - t(1,1) ) * ( t(2,2)   - t(2,1) ) )

        ref(2,j) = ( ( t(1,2) - t(1,1) ) * ( phy(2,j) - t(2,1) )   
     &             - ( t(2,2) - t(2,1) ) * ( phy(1,j) - t(1,1) ) ) 
     &           / ( ( t(2,3) - t(2,1) ) * ( t(1,2)   - t(1,1) )   
     &             - ( t(1,3) - t(1,1) ) * ( t(2,2)   - t(2,1) ) )

      end do

      return
      end
      subroutine triangle_order3_reference_to_physical ( t, n, ref, 
     &  phy )

c*********************************************************************72
c
cc TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL maps T3 reference points to physical points.
c
c  Discussion:
c
c    Given the vertices of an order 3 physical triangle and a point
c    (XSI,ETA) in the reference triangle, the routine computes the value
c    of the corresponding image point (X,Y) in physical space.
c
c    This routine is also appropriate for an order 4 triangle,
c    as long as the fourth node is the centroid of the triangle.
c
c    This routine may also be appropriate for an order 6
c    triangle, if the mapping between reference and physical space
c    is linear.  This implies, in particular, that the sides of the
c    image triangle are straight and that the "midside" nodes in the
c    physical triangle are halfway along the sides of
c    the physical triangle.
c
c  Reference Element T3:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  |  \
c    |  |   \
c    |  |    \
c    0  1-----2
c    |
c    +--0--R--1-->
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
c    Input, double precision T(2,3), the coordinates of the vertices.
c    The vertices are assumed to be the images of (0,0), (1,0) and
c    (0,1) respectively.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision REF(2,N), points in the reference triangle.
c
c    Output, double precision PHY(2,N), corresponding points in the
c    physical triangle.
c
      implicit none

      integer n

      integer i
      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,3)

      do i = 1, 2
        do j = 1, n
          phy(i,j) = t(i,1) * ( 1.0D+00 - ref(1,j) - ref(2,j) )
     &             + t(i,2) *             ref(1,j)
     &             + t(i,3) *                        ref(2,j)
        end do
      end do

      return
      end
      subroutine triangle_order6_physical_to_reference ( t, n, phy, 
     &  ref )

c*********************************************************************72
c
cc TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE maps a physical point to a reference point.
c
c  Discussion:
c
c    Given the vertices of an order 6 physical triangle and a point
c    (X,Y) in the physical triangle, the routine computes the value
c    of the corresponding image point (R,S) in reference space.
c
c    The mapping from (R,S) to (X,Y) has the form:
c
c      X(R,S) = A1 * R * R + B1 * R * S + C1 * S * S
c             + D1 * R     + E1 * S     + F1
c
c      Y(R,S) = A2 * R * R + B2 * R * S + C2 * S * S
c             + D2 * R     + E2 * S     + F2
c
c  Reference Element T3:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  6  5
c    |  |   \
c    |  |    \
c    0  1--4--2
c    |
c    +--0--R--1-->
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
c  Parameters:
c
c    Input, double precision T(2,6), the coordinates of the vertices.  
c    The vertices are assumed to be the images of (0,0), (1,0), (0,1), 
c    (1/2,0), (1/2,1/2) and (0,1/2), in that order.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision PHY(2,N), the coordinates of points in the 
c    physical space.
c
c    Output, double precision REF(2,N), the coordinates of the corresponding
c    points in the reference space.
c
      implicit none

      integer n

      double precision a(2)
      double precision b(2)
      double precision c(2)
      double precision d(2)
      double precision det
      integer dim
      double precision dx(2)
      double precision e(2)
      double precision f(2)
      double precision fun(2)
      double precision fun_norm
      integer i
      integer it
      integer j
      double precision jac(2,2)
      integer it_max
      parameter ( it_max = 10 )
      double precision it_tol
      parameter ( it_tol = 0.000001D+00 )
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,6)
c
c  Set iteration parameters.
c
      do i = 1, 2

        a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2) 
     &         - 4.0D+00 * t(i,4)

        b(i) =   4.0D+00 * t(i,1)
     &         - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

        c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3)
     &                                               - 4.0D+00 * t(i,6)

        d(i) = - 3.0D+00 * t(i,1) -           t(i,2)
     &         + 4.0D+00 * t(i,4)

        e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3)
     &                                               + 4.0D+00 * t(i,6)
        f(i) =             t(i,1)

      end do
c
c  Initialize the points by inverting the linear map.
c
      call triangle_order3_physical_to_reference ( t(1:2,1:3), n, phy, 
     &  ref )
c
c  Carry out the Newton iteration.
c
      do j = 1, n

        do it = 1, it_max

          do dim = 1, 2
            fun(dim) = a(dim) * ref(1,j) * ref(1,j) 
     &             + b(dim) * ref(1,j) * ref(2,j) 
     &             + c(dim) * ref(2,j) * ref(2,j) 
     &             + d(dim) * ref(1,j) 
     &             + e(dim) * ref(2,j) 
     &             + f(dim) 
     &             - phy(dim,j)

          end do

          fun_norm = sqrt ( fun(1) * fun(1) + fun(2) * fun(2) )

          if ( fun_norm .le. it_tol ) then
            go to 10
          end if

          do dim = 1, 2
            jac(dim,1) = 2.0D+00 * a(dim) * ref(1,j) 
     &                 +           b(dim) * ref(2,j) + d(dim)

            jac(dim,2) =           b(dim) * ref(1,j) 
     &                 + 2.0D+00 * c(dim) * ref(2,j) + e(dim)

          end do

          det = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)

          if ( det .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 
     &        'TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE - Fatal error!'
            write ( *, '(a)' ) 
     &        '  The jacobian of the mapping is singular.'
          end if

          dx(1) = (  jac(2,2) * fun(1) - jac(1,2) * fun(2) ) / det
          dx(2) = ( -jac(2,1) * fun(1) + jac(1,1) * fun(2) ) / det

          do dim = 1, 2
            ref(dim,j) = ref(dim,j) - dx(dim)
          end do

        end do

10      continue

      end do

      return
      end
      subroutine triangle_order6_reference_to_physical ( t, n, ref, 
     &  phy )

c*********************************************************************72
c
cc TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL maps T6 reference points to physical points.
c
c  Discussion:
c
c    Given the vertices of an order 6 physical triangle and a point
c    (XSI,ETA) in the reference triangle, the routine computes the value
c    of the corresponding image point (X,Y) in physical space.
c
c    The mapping from (XSI,ETA) to (X,Y) has the form:
c
c      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
c                 + D1 * XSI    + E1 * ETA     + F1
c
c      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
c                 + D2 * XSI    + E2 * ETA     + F2
c
c  Reference Element T6:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  6  5
c    |  |   \
c    |  |    \
c    0  1--4--2
c    |
c    +--0--R--1-->
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
c  Parameters:
c
c    Input, double precision T(2,6), the coordinates of the vertices.
c    The vertices are assumed to be the images of (0,0), (1,0),
c    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision REF(2,N), points in the reference triangle.
c
c    Output, double precision PHY(2,N), corresponding points in the
c    physical triangle.
c
      implicit none

      integer n

      double precision a(2)
      double precision b(2)
      double precision c(2)
      double precision d(2)
      double precision e(2)
      double precision f(2)
      integer i
      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,6)

      do i = 1, 2

        a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    
     &         - 4.0D+00 * t(i,4)

        b(i) =   4.0D+00 * t(i,1)                                       
     &         - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

        c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) 
     &                                               - 4.0D+00 * t(i,6)

        d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    
     &         + 4.0D+00 * t(i,4)

        e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) 
     &                                               + 4.0D+00 * t(i,6)
        f(i) =             t(i,1)

      end do

      do i = 1, 2
        do j = 1, n
          phy(i,j) = a(i) * ref(1,j) * ref(1,j) 
     &             + b(i) * ref(1,j) * ref(2,j) 
     &             + c(i) * ref(2,j) * ref(2,j) 
     &             + d(i) * ref(1,j)
     &             + e(i) * ref(2,j)
     &             + f(i)
        end do
      end do

      return
      end
      subroutine triangle_reference_sample ( n, seed, p )

c*********************************************************************72
c
cc TRIANGLE_REFERENCE_SAMPLE returns random points in the reference triangle.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  |   \  side 2
c       |    \
c    3  |     \
c       |      \
c       1-------2
c
c         side 1
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
c  Parameters:
c
c    Input, integer N, the number of points to generate.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision P(2,N), random points in the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer n

      double precision alpha
      double precision beta
      integer j
      double precision p(dim_num,n)
      double precision r
      double precision r8_uniform_01
      integer seed

      do j = 1, n

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
c  Now choose, uniformly at random, a point on the line L.
c
        beta = r8_uniform_01 ( seed )

        p(1,j) = ( 1.0D+00 - beta ) * alpha
        p(2,j) =             beta   * alpha

      end do

      return
      end
      subroutine triangle_sample ( t, n, seed, p )

c*********************************************************************72
c
cc TRIANGLE_SAMPLE returns random points in a triangle.
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
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, integer N, the number of points to generate.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision P(2,N), random points in the triangle.
c
      implicit none

      integer dim_num 
      parameter ( dim_num = 2 )
      integer n

      double precision alpha(n)
      integer dim
      integer i
      integer j
      double precision p(dim_num,n)
      double precision p12(dim_num,n)
      double precision p13(dim_num,n)
      integer seed
      double precision t(dim_num,3)

      call r8vec_uniform_01 ( n, seed, alpha )
c
c  Interpret R as a percentage of the triangle's area.
c
c  Imagine a line L, parallel to side 1, so that the area between
c  vertex 1 and line L is R percent of the full triangle's area.
c
c  The line L will intersect sides 2 and 3 at a fraction
c  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
c
      do i = 1, n
        alpha(i) = sqrt ( alpha(i) )
      end do
c
c  Determine the coordinates of the points on sides 2 and 3 intersected
c  by line L.
c
      do dim = 1, dim_num

        do j = 1, n
          p12(dim,j) = ( 1.0D+00 - alpha(j) ) * t(dim,1) 
     &                           + alpha(j)   * t(dim,2)

          p13(dim,j) = ( 1.0D+00 - alpha(j) ) * t(dim,1) 
     &                           + alpha(j)   * t(dim,3)
        end do
      end do
c
c  Now choose, uniformly at random, a point on the line L.
c
      call r8vec_uniform_01 ( n, seed, alpha )

      do dim = 1, dim_num
        do j = 1, n
          p(dim,j) = ( 1.0D+00 - alpha(j) ) * p12(dim,j) 
     &                         + alpha(j)   * p13(dim,j)
        end do
      end do

      return
      end
      subroutine triangulation_node_order ( triangle_order, 
     &  triangle_num, triangle_node, node_num, node_order )

c*********************************************************************72
c
cc TRIANGULATION_NODE_ORDER determines the order of nodes in a triangulation.
c
c  Discussion:
c
c    The order of a node is the number of triangles that use that node
c    as a vertex.
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
c  Parameters:
c
c    Input, integer TRIANGLE_ORDER, the order of the triangles.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
c    the nodes that make up the triangles.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Output, integer NODE_ORDER(NODE_NUM), the order of each node.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order

      integer i
      integer node
      integer node_order(node_num)
      integer triangle
      integer triangle_node(triangle_order,triangle_num)

      node_order(1:node_num) = 0

      do triangle = 1, triangle_num
        do i = 1, triangle_order
          node = triangle_node(i,triangle)
          if ( node < 1 .or. node_num < node ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TRIANGULATION_NODE_ORDER - Fatal error!'
            write ( *, '(a)' ) '  Illegal entry in TRIANGLE_NODE.'
            stop
          else
            node_order(node) = node_order(node) + 1
          end if
        end do
      end do

      return
      end
      subroutine triangulation_order3_adj_count ( node_num, 
     &  triangle_num, triangle_node, triangle_neighbor, adj_num, 
     &  adj_col )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The triangulation is assumed to involve 3-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  |   \  side 2
c       |    \
c    3  |     \
c       |      \
c       1-------2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
c
c    A sample grid.
c
c
c    Below, we have a chart that summarizes the adjacency relationships
c    in the sample grid.  On the left, we list the node, and its neighbors,
c    with an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).  On the right, we list the number of adjacencies to
c    lower-indexed nodes, to the node itself, to higher-indexed nodes,
c    the total number of adjacencies for this node, and the location
c    of the first and last entries required to list this set of adjacencies
c    in a single list of all the adjacencies.
c
c    N   Adjacencies                Below  Self   Above   Total First  Last
c
c   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
c    1:  *  2  6                        0     1       2       3     1     3
c    2:  1  *  3  6  7                  1     1       3       5     4     8
c    3:  2  *  4  7  8                  1     1       3       5     9    13
c    4:  3  *  5  8  9                  1     1       3       5    14    18
c    5:  4  *  9 10                     1     1       2       4    19    22
c    6:  1  2  *  7 11                  2     1       2       5    23    27
c    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
c    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
c    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
c   10:  5  9  * 14 15                  2     1       2       5    49    53
c   11:  6  7  * 12 16                  2     1       2       5    54    58
c   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
c   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
c   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
c   15: 10 14  * 19 20                  2     1       2       5    80    84
c   16: 11 12  * 17 21                  2     1       2       5    85    89
c   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
c   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
c   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
c   20: 15 19  * 24 25                  2     1       2       5   111   115
c   21: 16 17  * 22                     2     1       1       4   116   119
c   22: 17 18 21  * 23                  3     1       1       5   120   124
c   23: 18 19 22  * 24                  3     1       1       5   125   129
c   24: 19 20 23  * 25                  3     1       1       5   130   134
c   25: 20 24  *                        2     1       0       3   135   137
c   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle, in counterclockwise order. 
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Output, integer ADJ_NUM, the number of adjacencies.
c
c    Output, integer ADJ_COL(NODE_NUM+1).  Information about 
c    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj_num
      integer adj_col(node_num+1)
      integer i
      integer n1
      integer n2
      integer n3
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      adj_num = 0
c
c  Set every node to be adjacent to itself.
c
      do i = 1, node_num
        adj_col(i) = 1
      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
c
c  Add edge (1,2) if this is the first occurrence,
c  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE < TRIANGLE2).
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n2) = adj_col(n2) + 1
        end if
c
c  Add edge (2,3).
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n3) = adj_col(n3) + 1
        end if
c
c  Add edge (3,1).
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n3) = adj_col(n3) + 1
        end if
          
      end do
c
c  We used ADJ_COL to count the number of entries in each column.
c  Convert it to pointers into the ADJ array.
c
      do i = node_num, 1, -1
        adj_col(i+1) = adj_col(i)
      end do

      adj_col(1) = 1
      do i = 2, node_num+1
        adj_col(i) = adj_col(i-1) + adj_col(i)
      end do

      adj_num = adj_col(node_num+1) - 1

      return
      end
      subroutine triangulation_order3_adj_set ( node_num, triangle_num, 
     &  triangle_node, triangle_neighbor, adj_num, adj_col, adj )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The triangulation is assumed to involve 3-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c    This routine can be used to create the compressed column storage
c    for a linear triangle finite element discretization of 
c    Poisson's equation in two dimensions.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  |   \  side 2
c       |    \
c    3  |     \
c       |      \
c       1-------2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
c
c    A sample grid
c
c
c    Below, we have a chart that summarizes the adjacency relationships
c    in the sample grid.  On the left, we list the node, and its neighbors,
c    with an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).  On the right, we list the number of adjacencies to
c    lower-indexed nodes, to the node itself, to higher-indexed nodes,
c    the total number of adjacencies for this node, and the location
c    of the first and last entries required to list this set of adjacencies
c    in a single list of all the adjacencies.
c
c    N   Adjacencies                Below  Self    Above  Total First  Last
c
c   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
c    1:  *  2  6                        0     1       2       3     1     3
c    2:  1  *  3  6  7                  1     1       3       5     4     8
c    3:  2  *  4  7  8                  1     1       3       5     9    13
c    4:  3  *  5  8  9                  1     1       3       5    14    18
c    5:  4  *  9 10                     1     1       2       4    19    22
c    6:  1  2  *  7 11                  2     1       2       5    23    27
c    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
c    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
c    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
c   10:  5  9  * 14 15                  2     1       2       5    49    53
c   11:  6  7  * 12 16                  2     1       2       5    54    58
c   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
c   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
c   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
c   15: 10 14  * 19 20                  2     1       2       5    80    84
c   16: 11 12  * 17 21                  2     1       2       5    85    89
c   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
c   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
c   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
c   20: 15 19  * 24 25                  2     1       2       5   111   115
c   21: 16 17  * 22                     2     1       1       4   116   119
c   22: 17 18 21  * 23                  3     1       1       5   120   124
c   23: 18 19 22  * 24                  3     1       1       5   125   129
c   24: 19 20 23  * 25                  3     1       1       5   130   134
c   25: 20 24  *                        2     1       0       3   135   137
c   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists the nodes
c    that make up each triangle in counterclockwise order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer ADJ_NUM, the number of adjacencies.
c
c    Input, integer ADJ_COL(NODE_NUM+1).  Information about 
c    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
c    Output, integer ADJ(ADJ_NUM), the adjacency information.
c
      implicit none

      integer adj_num
      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj(adj_num)
      integer adj_col(node_num+1)
      integer adj_copy(node_num)
      integer i
      integer k1
      integer k2
      integer n1
      integer n2
      integer n3
      integer node
      integer number
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      do i = 1, adj_num
        adj(1:adj_num) = -1
      end do

      do node = 1, node_num
        adj_copy(node) = adj_col(node)
      end do
c
c  Set every node to be adjacent to itself.
c
      do node = 1, node_num
        adj(adj_copy(node)) = node
        adj_copy(node) = adj_copy(node) + 1
      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
c
c  Add edge (1,2) if this is the first occurrence,
c  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n1)) = n2
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n2)) = n1
          adj_copy(n2) = adj_copy(n2) + 1
        end if
c
c  Add edge (2,3).
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n2)) = n3
          adj_copy(n2) = adj_copy(n2) + 1
          adj(adj_copy(n3)) = n2
          adj_copy(n3) = adj_copy(n3) + 1
        end if
c
c  Add edge (3,1).
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n1)) = n3
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n3)) = n1
          adj_copy(n3) = adj_copy(n3) + 1
        end if
          
      end do
c
c  Ascending sort the entries for each node.
c
      do node = 1, node_num
        k1 = adj_col(node)
        k2 = adj_col(node+1)-1
        number = k2 + 1 - k1
        call i4vec_sort_heap_a ( number, adj(k1:k2) )
      end do

      return
      end
      subroutine triangulation_order3_adj_set2 ( node_num, triangle_num, 
     &  triangle_node, triangle_neighbor, adj_num, adj_col, ia, ja )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to set up the arrays IA and JA that
c    record which nodes are adjacent in a triangulation.
c
c    The triangulation is assumed to involve 3-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c    This routine can be used to set up the sparse triplet storage
c    for a linear triangle finite element discretization of Poisson's 
c    equation in two dimensions.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  |   \  side 2
c       |    \
c    3  |     \
c       |      \
c       1-------2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
c
c    A sample grid
c
c
c    Below, we have a chart that summarizes the adjacency relationships
c    in the sample grid.  On the left, we list the node, and its neighbors,
c    with an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).  On the right, we list the number of adjacencies to
c    lower-indexed nodes, to the node itself, to higher-indexed nodes,
c    the total number of adjacencies for this node, and the location
c    of the first and last entries required to list this set of adjacencies
c    in a single list of all the adjacencies.
c
c    N   Adjacencies                Below  Self    Above  Total First  Last
c
c   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
c    1:  *  2  6                        0     1       2       3     1     3
c    2:  1  *  3  6  7                  1     1       3       5     4     8
c    3:  2  *  4  7  8                  1     1       3       5     9    13
c    4:  3  *  5  8  9                  1     1       3       5    14    18
c    5:  4  *  9 10                     1     1       2       4    19    22
c    6:  1  2  *  7 11                  2     1       2       5    23    27
c    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
c    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
c    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
c   10:  5  9  * 14 15                  2     1       2       5    49    53
c   11:  6  7  * 12 16                  2     1       2       5    54    58
c   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
c   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
c   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
c   15: 10 14  * 19 20                  2     1       2       5    80    84
c   16: 11 12  * 17 21                  2     1       2       5    85    89
c   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
c   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
c   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
c   20: 15 19  * 24 25                  2     1       2       5   111   115
c   21: 16 17  * 22                     2     1       1       4   116   119
c   22: 17 18 21  * 23                  3     1       1       5   120   124
c   23: 18 19 22  * 24                  3     1       1       5   125   129
c   24: 19 20 23  * 25                  3     1       1       5   130   134
c   25: 20 24  *                        2     1       0       3   135   137
c   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
c
c    For this example, the initial portion of the IA and JA arrays will be:
c
c      (1,1), (1,2), (1,6),
c      (2,1), (2,2), (2,3), (2,6), (2,7),
c      (3,2), (3,3), (3,4), (3,7), (3,8),
c      ...
c      (25,20), (25,24), (25,25)
c
c    for a total of 137 pairs of values.
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists the nodes
c    that make up each triangle in counterclockwise order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer ADJ_NUM, the number of adjacencies.
c
c    Input, integer ADJ_COL(NODE_NUM+1).  Information about 
c    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
c    Output, integer IA(ADJ_NUM), JA(ADJ_NUM), the adjacency
c    information.
c
      implicit none

      integer adj_num
      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj_col(node_num+1)
      integer adj_copy(node_num)
      integer i
      integer ia(adj_num)
      integer ja(adj_num)
      integer k1
      integer k2
      integer n1
      integer n2
      integer n3
      integer node
      integer number
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      do i = 1, adj_num
        ia(i) = -1
      end do

      do i = 1, adj_num
        ja(i) = -1
      end do

      do node = 1, node_num
        adj_copy(node) = adj_col(node)
      end do
c
c  Set every node to be adjacent to itself.
c
      do node = 1, node_num
        ia(adj_copy(node)) = node
        ja(adj_copy(node)) = node
        adj_copy(node) = adj_copy(node) + 1
      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
c
c  Add edge (1,2) if this is the first occurrence,
c  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          ia(adj_copy(n1)) = n1
          ja(adj_copy(n1)) = n2
          adj_copy(n1) = adj_copy(n1) + 1

          ia(adj_copy(n2)) = n2
          ja(adj_copy(n2)) = n1
          adj_copy(n2) = adj_copy(n2) + 1

        end if
c
c  Add edge (2,3).
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          ia(adj_copy(n2)) = n2
          ja(adj_copy(n2)) = n3
          adj_copy(n2) = adj_copy(n2) + 1

          ia(adj_copy(n3)) = n3
          ja(adj_copy(n3)) = n2
          adj_copy(n3) = adj_copy(n3) + 1

        end if
c
c  Add edge (3,1).
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then

          ia(adj_copy(n1)) = n1
          ja(adj_copy(n1)) = n3
          adj_copy(n1) = adj_copy(n1) + 1

          ia(adj_copy(n3)) = n3
          ja(adj_copy(n3)) = n1
          adj_copy(n3) = adj_copy(n3) + 1

        end if
          
      end do
c
c  Lexically sort the IA, JA values.
c
      call i4vec2_sort_a ( adj_num, ia, ja )

      return
      end
      subroutine triangulation_order3_boundary_edge_count ( 
     &  triangle_num, triangle_node, boundary_edge_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the boundary edges.
c
c  Discussion:
c
c    This routine is given a triangulation, an abstract list of triples
c    of nodes.  It is assumed that the nodes in each triangle are listed
c    in a counterclockwise order, although the routine should work 
c    if the nodes are consistently listed in a clockwise order as well.
c
c    It is assumed that each edge of the triangulation is either 
c    * an INTERIOR edge, which is listed twice, once with positive
c      orientation and once with negative orientation, or;
c    * a BOUNDARY edge, which will occur only once.
c
c    This routine should work even if the region has holes - as long
c    as the boundary of the hole comprises more than 3 edgesc
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise 
c    order.
c
c    Output, integer BOUNDARY_EDGE_NUM, the number of boundary 
c    edges.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer boundary_edge_num
      integer edge(2,3*triangle_num)
      integer i
      integer interior_edge_num
      integer j
      integer k
      integer m
      integer n
      integer triangle_node(triangle_order,triangle_num)
      integer unique_num

      m = 2
      n = 3 * triangle_num
c
c  Set up the edge array.
c
      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j) = triangle_node(i,j)
        end do
      end do

      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j+triangle_num) = triangle_node(i+1,j)
        end do
      end do

      do j = 1, triangle_num
        edge(1,j+2*triangle_num) = triangle_node(3,j)
        edge(2,j+2*triangle_num) = triangle_node(1,j)
      end do
c
c  In each column, force the smaller entry to appear first.
c
      do j = 1, n
        i = min ( edge(1,j), edge(2,j) )
        k = max ( edge(1,j), edge(2,j) )
        edge(1,j) = i
        edge(2,j) = k
      end do
c
c  Ascending sort the column array.
c
      call i4col_sort_a ( m, n, edge )
c
c  Get the number of unique columns in EDGE.
c
      call i4col_sorted_unique_count ( m, n, edge, unique_num )

      interior_edge_num = 3 * triangle_num - unique_num

      boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num

      return
      end
      subroutine triangulation_order3_boundary_edge_count_euler ( 
     &  node_num, triangle_num, hole_num, boundary_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
c
c  Discussion:
c
c    We assume we are given information about a triangulation
c    of a set of nodes in the plane.
c
c    Given the number of nodes, triangles and holes, we are going to apply
c    Euler's formula to determine the number of edges that lie on the
c    boundary of the set of nodes.
c
c    The number of faces, including the infinite face and internal holes, 
c    is TRIANGLE_NUM + HOLE_NUM + 1.
c
c    Let BOUNDARY_NUM denote the number of edges on the boundary.
c    Each of the TRIANGLE_NUM triangles uses three edges.  Every edge
c    occurs in two different faces, so the number of edges must be
c    ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2.
c
c    The number of nodes used in the triangulation is NODE_NUM.
c
c    Euler's formula asserts that, for a simple connected figure in the
c    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
c    FACE_NUM faces:
c
c      NODE_NUM - EDGE_NUM + FACE_NUM = 2
c
c    In our context, this becomes
c
c      NODE_NUM - ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2 
c      + TRIANGLE_NUM + HOLE_NUM + 1 = 2
c
c    or
c
c      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - TRIANGLE_NUM - 2
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
c  Reference:
c
c    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
c    Computational Geometry, Section 9.1,
c    Springer, 2000.
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer HOLE_NUM, the number of internal holes.
c
c    Output, integer BOUNDARY_NUM, the number of edges that 
c    lie on the boundary of the triangulation.
c
      implicit none

      integer boundary_num
      integer hole_num
      integer node_num
      integer triangle_num

      boundary_num = 2 * node_num + 2 * hole_num - triangle_num - 2

      return
      end
      subroutine triangulation_order3_boundary_node ( node_num, 
     &  triangle_num, triangle_node, node_boundary )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_BOUNDARY_NODE indicates which nodes are on the boundary.
c
c  Discussion:
c
c    This routine is given a triangulation, an abstract list of triples
c    of nodes.  It is assumed that the nodes in each triangle are listed
c    in a counterclockwise order, although the routine should work 
c    if the nodes are consistently listed in a clockwise order as well.
c
c    It is assumed that each edge of the triangulation is either 
c    * an INTERIOR edge, which is listed twice, once with positive
c      orientation and once with negative orientation, or;
c    * a BOUNDARY edge, which will occur only once.
c
c    This routine should work even if the region has holes - as long
c    as the boundary of the hole comprises more than 3 edges!
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise 
c    order.
c
c    Output, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node
c    is on a boundary edge.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer e1(3*triangle_num)
      integer e2(3*triangle_num)
      integer edge(2,3*triangle_num)
      integer i
      integer j
      integer k
      integer m
      integer n
      integer node
      logical node_boundary(node_num)
      integer triangle_node(triangle_order,triangle_num)

      m = 2
      n = 3 * triangle_num
c
c  Set up the edge array.
c
      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j) = triangle_node(i,j)
        end do
      end do
      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j+triangle_num) = triangle_node(i+1,j)
        end do
      end do
      do j = 1, triangle_num
        edge(1,j+2*triangle_num) = triangle_node(3,j)
        edge(2,j+2*triangle_num) = triangle_node(1,j)
      end do
c
c  In each column, force the smaller entry to appear first.
c
      do j = 1, n
        i = min ( edge(1,j), edge(2,j) )
        k = max ( edge(1,j), edge(2,j) )
        edge(1,j) = i
        edge(2,j) = k
      end do
c
c  Ascending sort the column array.
c
      call i4col_sort_a ( m, n, edge )
c
c  Records which appear twice are internal edges and can be ignored.
c
      do node = 1, node_num
        node_boundary(node) = .false.
      end do

      j = 0

10    continue

      if ( j .lt. 3 * triangle_num ) then

        j = j + 1

        if ( j .eq. 3 * triangle_num ) then
          node_boundary(edge(1,j)) = .true.
          node_boundary(edge(2,j)) = .true.
        else if ( edge(1,j) .eq. edge(1,j+1) .and.
     &            edge(2,j) .eq. edge(2,j+1) ) then
          j = j + 1
        else
          node_boundary(edge(1,j)) = .true.
          node_boundary(edge(2,j)) = .true.
        end if

        go to 10

      end if

      return
      end
      subroutine triangulation_order3_check ( node_num, triangle_num, 
     &  triangle_node, ierror )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_CHECK makes some simple checks on a triangulation.
c
c  Discussion:
c
c    Because this routine does not receive the physical coordinates of
c    the nodes, it cannot ensure that the triangulation is maximal,
c    that is, that no more triangles can be created.
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise 
c    order.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    nonzero, an error occurred, the triangulation is not valid.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer boundary_num
      integer euler
      integer i
      integer ierror
      integer j
      integer triangle_node(triangle_order,triangle_num)
      integer used(node_num)

      ierror = 0
c
c  Checks 1 and 2:
c  NODE_NUM must be at least 3.
c  TRIANGLE_NUM must be at least 1.
c
      if ( node_num .lt. 3 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
        write ( *, '(a)' ) '  The number of nodes is less than 3!'
        return
      end if

      if ( triangle_num .lt. 1 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
        write ( *, '(a)' ) '  The number of triangles is less than 1!'
        return
      end if
c
c  Checks 3 and 4:
c  Verify that all node values are greater than or equal to 1
c  and less than or equal to NODE_NUM.
c
      do i = 1, 3
        do j = 1, triangle_num
          if ( triangle_node(i,j) .lt. 1 ) then
            ierror = 3
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
            write ( *, '(a)' ) '  Some nodes are less than 1!'
            return
          end if
        end do
      end do

      do i = 1, 3
        do j = 1, triangle_num
          if ( node_num .lt. triangle_node(i,j) ) then
            ierror = 4
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
            write ( *, '(a)' ) '  Some nodes are greater than NODE_NUM!'
            return
          end if
        end do
      end do
c
c  Check 5:
c  Verify that every node is used at least once.
c
      do i = 1, node_num
        used(i) = 0
      end do

      do i = 1, 3
        do j = 1, triangle_num
          used(triangle_node(i,j)) = used(triangle_node(i,j)) + 1
        end do
      end do

      do i = 1, node_num
        if ( used(i) .eq. 0 ) then
          ierror = 5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
          write ( *, '(a)' ) 
     &      '  Some nodes are never used as triangle vertices!'
          return
        end if
      end do
c
c  Check 6:
c  Verify that no node is repeated in a triangle.
c
      do i = 1, triangle_num
        if ( triangle_node(1,i) .eq. triangle_node(2,i) .or. 
     &       triangle_node(2,i) .eq. triangle_node(3,i) .or. 
     &       triangle_node(3,i) .eq. triangle_node(1,i) ) then
          ierror = 6
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
          write ( *, '(a)' ) '  A triangle contains a null edge!'
          return
        end if
      end do
c
c  Check 7:
c  Verify that no edge is repeated, and that repeated edges occur in
c  negated pairs.
c
      call triangulation_order3_edge_check ( triangle_num, 
     &  triangle_node, boundary_num, ierror )

      if ( ierror .ne. 0 ) then
        ierror = 7
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
        write ( *, '(a)' ) '  Some edges are repeated,'
        write ( *, '(a)' ) '  or given in the wrong direction!'
        return
      end if
c
c  Check 8:
c  Does the triangulation satisfy Euler's criterion?
c  If not, then the triangulation is not proper.  (For instance, there
c  might be a hole in the interior.)
c
      euler = boundary_num + triangle_num + 2 - 2 * node_num

      if ( euler .ne. 0 ) then
        ierror = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_ORDER3_CHECK - Warning!'
        write ( *, '(a)' ) 
     &    '  The triangulation fails Euler''s criterion!'
        return
      end if

      return
      end
      subroutine triangulation_order3_edge_check ( triangle_num, 
     &  triangle_node, boundary_num, ierror )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_EDGE_CHECK checks the edges of a triangulation.
c
c  Discussion:
c
c    This routine was revised to store the edge data as columns,
c    rather than rows.
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer BOUNDARY_NUM, the number of edges that 
c    lie on the boundary.
c
c    Output, integer IERROR, an error flag.
c    0, no errors were detected.
c    nonzero, an error occurred.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer boundary_num
      integer col(3,3*triangle_num)
      integer i
      integer ierror
      integer j
      integer k
      integer tri
      integer triangle_node(triangle_order,triangle_num)

      ierror = 0
c
c  Step 1.
c  From the list of nodes for triangle T, of the form: (I,J,K)
c  construct the three neighbor relations:
c
c    (I,J,+1) or (J,I,-1),
c    (J,K,+1) or (K,J,-1),
c    (K,I,+1) or (I,K,-1)
c
c  where we choose (I,J,+1) if I .lt. J, or else (J,I,-1) and so on.
c
      do tri = 1, triangle_num

        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)

        if ( i .lt. j ) then
          col(1,3*(tri-1)+1) = i
          col(2,3*(tri-1)+1) = j
          col(3,3*(tri-1)+1) = + 1
        else
          col(1,3*(tri-1)+1) = j
          col(2,3*(tri-1)+1) = i
          col(3,3*(tri-1)+1) = - 1
        end if

        if ( j .lt. k ) then
          col(1,3*(tri-1)+2) = j
          col(2,3*(tri-1)+2) = k
          col(3,3*(tri-1)+2) = + 1
        else
          col(1,3*(tri-1)+2) = k
          col(2,3*(tri-1)+2) = j
          col(3,3*(tri-1)+2) = - 1
        end if

        if ( k .lt. i ) then
          col(1,3*(tri-1)+3) = k
          col(2,3*(tri-1)+3) = i
          col(3,3*(tri-1)+3) = + 1
        else
          col(1,3*(tri-1)+3) = i
          col(2,3*(tri-1)+3) = k
          col(3,3*(tri-1)+3) = - 1
        end if

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c
      call i4col_sort_a ( 3, 3*triangle_num, col )
c
c  Step 3.
c
c  Most records (/ A, B, C /) occur twice, with C being -1, then +1.
c  These records represent internal edges.
c
c  Unpaired records represent edges on the boundary.
c
      i = 0
      boundary_num = 0

10    continue

      if ( i .lt. 3 * triangle_num ) then

        i = i + 1

        if ( i .eq. 3 * triangle_num ) then

          boundary_num = boundary_num + 1

        else

          if ( col(1,i) .eq. col(1,i+1) .and. 
     &         col(2,i) .eq. col(2,i+1) ) then

            if ( col(3,i) .eq. col(3,i+1) ) then
              ierror = 1
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 
     &          'TRIANGULATION_ORDER3_EDGE_CHECK - Warning!'
              write ( *, '(a)' ) '  An edge occurs twice.'
              return
            end if

            i = i + 1

          else

            boundary_num = boundary_num + 1

          end if

        end if

        go to 10

      end if

      return
      end
      subroutine triangulation_order3_example1 ( node_num, 
     &  triangle_num, node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_EXAMPLE1 sets up a sample triangulation.
c
c  Discussion:
c
c    This triangulation is actually a Delaunay triangulation.
c
c    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
c    determined by calling TRIANGULATION_ORDER3_EXAMPLE1_SIZE first.
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.  
c
c    Input, integer TRIANGLE_NUM, the number of triangles.  
c
c    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of the
c    nodes.
c
c    Output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes
c    that make up the triangles.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbors on each side.  Negative values indicate edges that 
c    lie on the exterior.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer node_num_save
      parameter ( node_num_save = 13 )
      integer triangle_num
      integer triangle_num_save
      parameter ( triangle_num_save = 16 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision node_xy(dim_num,node_num)
      double precision node_xy_save(dim_num,node_num_save)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_neighbor_save(3,triangle_num_save)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_node_save(triangle_order,triangle_num_save)

      save node_xy_save
      save triangle_neighbor_save
      save triangle_node_save

      data node_xy_save /
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

      data triangle_neighbor_save /
     &     -4,  -13,    2, 
     &      1,    4,    3, 
     &      2,    5,    7, 
     &      2,  -43,    8, 
     &      3,    8,    6, 
     &      5,    9,    7, 
     &      3,    6,   -3, 
     &      5,    4,   10, 
     &      6,   10,   12, 
     &      9,    8,   11, 
     &     12,   10,   14, 
     &      9,   11,   13, 
     &    -23,   12,   16, 
     &     11,  -47,   15, 
     &     16,   14,  -50, 
     &     13,   15,  -39 /

      data triangle_node_save /
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

      call r8mat_copy ( dim_num, node_num, node_xy_save, node_xy )

      call i4mat_copy ( 3, triangle_num, triangle_neighbor_save,
     &  triangle_neighbor )

      call i4mat_copy ( triangle_order, triangle_num, 
     &  triangle_node_save, triangle_node )

      return
      end
      subroutine triangulation_order3_example1_size ( node_num, 
     &  triangle_num, hole_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_EXAMPLE1_SIZE sets sizes for a sample triangulation.
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
c  Parameters:
c
c    Output, integer NODE_NUM, the number of nodes.  
c
c    Output, integer TRIANGLE_NUM, the number of triangles. 
c
c    Output, integer HOLE_NUM, the number of holes.
c
      implicit none

      integer hole_num
      integer node_num
      integer triangle_num

      hole_num = 0
      node_num = 13
      triangle_num = 16

      return
      end
      subroutine triangulation_order3_example2 ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_EXAMPLE2 returns an example triangulation.
c
c  Discussion:
c
c    This triangulation is actually a Delaunay triangulation.
c
c    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
c    determined by calling TRIANGULATION_ORDER3_EXAMPLE2_SIZE first.
c
c  Diagram:
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of the
c    nodes.
c
c    Output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle, in counterclockwise order.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for 
c    each side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer node_num_save
      parameter ( node_num_save = 25 )
      integer triangle_num
      integer triangle_num_save
      parameter ( triangle_num_save = 32 )
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision node_xy(dim_num,node_num)
      double precision node_xy_save(dim_num,node_num_save)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_neighbor_save(3,triangle_num_save)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_node_save(triangle_order,triangle_num_save)

      save node_xy_save
      save triangle_neighbor_save
      save triangle_node_save

      data node_xy_save /
     &  0.0D+00, 0.0D+00,
     &  1.0D+00, 0.0D+00,
     &  2.0D+00, 0.0D+00,
     &  3.0D+00, 0.0D+00,
     &  4.0D+00, 0.0D+00,
     &  0.0D+00, 1.0D+00,
     &  1.0D+00, 1.0D+00,
     &  2.0D+00, 1.0D+00,
     &  3.0D+00, 1.0D+00,
     &  4.0D+00, 1.0D+00,
     &  0.0D+00, 2.0D+00,
     &  1.0D+00, 2.0D+00,
     &  2.0D+00, 2.0D+00,
     &  3.0D+00, 2.0D+00,
     &  4.0D+00, 2.0D+00,
     &  0.0D+00, 3.0D+00,
     &  1.0D+00, 3.0D+00,
     &  2.0D+00, 3.0D+00,
     &  3.0D+00, 3.0D+00,
     &  4.0D+00, 3.0D+00,
     &  0.0D+00, 4.0D+00,
     &  1.0D+00, 4.0D+00,
     &  2.0D+00, 4.0D+00,
     &  3.0D+00, 4.0D+00,
     &  4.0D+00, 4.0D+00 /

      data triangle_node_save /
     &   1,  2,  6,
     &   7,  6,  2,
     &   2,  3,  7,
     &   8,  7,  3,
     &   3,  4,  8,
     &   9,  8,  4,
     &   4,  5,  9,
     &  10,  9,  5,
     &   6,  7, 11,
     &  12, 11,  7,
     &   7,  8, 12,
     &  13, 12,  8,
     &   8,  9, 13,
     &  14, 13,  9,
     &   9, 10, 14,
     &  15, 14, 10,
     &  11, 12, 16,
     &  17, 16, 12,
     &  12, 13, 17,
     &  18, 17, 13,
     &  13, 14, 18,
     &  19, 18, 14,
     &  14, 15, 19,
     &  20, 19, 15,
     &  16, 17, 21,
     &  22, 21, 17,
     &  17, 18, 22,
     &  23, 22, 18,
     &  18, 19, 23,
     &  24, 23, 19,
     &  19, 20, 24,
     &  25, 24, 20 /

      data triangle_neighbor_save /
     &  -1,  2, -1,
     &   9,  1,  3,
     &  -1,  4,  2,
     &  11,  3,  5,
     &  -1,  6,  4,
     &  13,  5,  7,
     &  -1,  8,  6,
     &  15,  7, -1, 
     &   2, 10, -1,
     &  17,  9, 11,
     &   4, 12, 10,
     &  19, 11, 13,
     &   6, 14, 12,
     &  21, 13, 15,
     &   8, 16, 14,
     &  23, 15, -1, 
     &  10, 18, -1,
     &  25, 17, 19,
     &  12, 20, 18,
     &  27, 19, 21,
     &  14, 22, 20,
     &  29, 21, 23,
     &  16, 24, 22,
     &  31, 23, -1, 
     &  18, 26, -1,
     &  -1, 25, 27,
     &  20, 28, 26,
     &  -1, 27, 29,
     &  22, 30, 28,
     &  -1, 29, 31,
     &  24, 32, 30,
     &  -1, 31, -1 /

      call r8mat_copy ( dim_num, node_num, node_xy_save, node_xy )

      call i4mat_copy ( 3, triangle_num, triangle_neighbor_save,
     &  triangle_neighbor )

      call i4mat_copy ( triangle_order, triangle_num, 
     &  triangle_node_save, triangle_node )
      return
      end
      subroutine triangulation_order3_example2_size ( node_num, 
     &  triangle_num, hole_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_EXAMPLE2_SIZE returns the size of an example.
c
c  Diagram:
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
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
c  Parameters
c
c    Output, integer NODE_NUM, the number of nodes.
c
c    Output, integer TRIANGLE_NUM, the number of triangles.
c
c    Output, integer HOLE_NUM, the number of holes.
c
      implicit none

      integer hole_num
      integer node_num
      integer triangle_num

      node_num = 25
      triangle_num = 32
      hole_num = 0

      return
      end
      subroutine triangulation_order3_neighbor ( triangle_num, 
     &  triangle_node, t1, s1, t2, s2 )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_NEIGHBOR determines a neighbor of a given triangle.
c
c  Discussion:
c
c    A set of nodes is given.  A triangulation of the nodes has been
c    defined and recorded in TRIANGLE.  The TRIANGLE data structure records
c    triangles as sets of three nodes, N1, N2, N3, that implicitly define three
c    sides, being the line segments N1-N2, N2-N3, and N3-N1.
c
c    The nodes of the triangle are listed in counterclockwise order.
c    This means that if two triangles share a side, then the nodes
c    defining that side occur in the order (N1,N2) for one triangle,
c    and (N2,N1) for the other.
c
c    The routine is given a triangle and a side, and asked to find
c    another triangle (if any) that shares that side.  The routine
c    simply searches the TRIANGLE_NODE structure for an occurrence of the
c    nodes in the opposite order.
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that define each triangle.
c
c    Input, integer T1, the index of the triangle.
c
c    Input, integer S1, the index of the triangle side.
c
c    Output, integer T2, the index of the triangle which is 
c    the neighbor to T1 on side S1, or -1 if there is no such neighbor.
c
c    Output, integer S2, the index of the side of triangle T2 
c    which is shared with triangle T1, or -1 if there is no such neighbor.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i4_wrap
      integer n1
      integer n2
      integer s
      integer s1
      integer s2
      integer ss
      integer t
      integer t1
      integer t2
      integer triangle_node(triangle_order,triangle_num)

      t2 = - 1
      s2 = - 1

      n1 = triangle_node(s1,t1)
      ss = s1 + 1
      ss = i4_wrap ( ss, 1, 3 )
      n2 = triangle_node(ss,t1)

      do t = 1, triangle_num
        do s = 1, 3
          if ( triangle_node(s,t) .eq. n1 ) then
            ss = s - 1
            ss = i4_wrap ( ss, 1, 3 )
            if ( triangle_node(ss,t) .eq. n2 ) then
              t2 = t
              s2 = ss
              return
            end if
          end if
        end do
      end do

      return
      end
      subroutine triangulation_order3_neighbor_nodes ( node_num, 
     &  triangle_num, nabes_max, triangle_node, nabes_first, 
     &  nabes_num, nabes_dim, nabes )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_NEIGHBOR_NODES determines triangulation neighbor nodes.
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
c    05 June 2009
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
c    Input, integer NABES_MAX, the maximum dimension of NABES.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer NABES_FIRST(NODE_NUM), the index in NABES 
c    of the first neighbor in the list for each node.
c
c    Output, integer NABES_NUM(NODE_NUM), the number of neighbors
c    of each node.
c
c    Output, integer NABES_DIM, the dimension of NABES.
c
c    Output, integer NABES(NABES_DIM), a list of the neighbors 
c    of all the nodes.  Neighbors of node 1 are listed first, and so on.
c
      implicit none

      integer nabes_max
      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i
      integer i_current
      integer j
      integer k
      integer nabe
      integer nabes(nabes_max)
      integer nabes1(nabes_max)
      integer nabes_dim
      integer nabes_first(node_num)
      integer nabes_num(node_num)
      integer node
      integer tri
      integer triangle_node(triangle_order,triangle_num)
      integer unique_num
c
c  Step 1.  From the triangle list (I,J,K)
c  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
c
      nabes_dim = 0
      do tri = 1, triangle_num
        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)
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
      call i4vec2_sorted_unique ( nabes_dim, nabes1, nabes, unique_num )

      nabes_dim = unique_num
c
c  Step 4. Construct the NABES_NUM and NABES_FIRST data.
c
      do node = 1, node_num
        nabes_num(node) = 0
      end do

      do node = 1, node_num
        nabes_first(node) = 0
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
      subroutine triangulation_order3_neighbor_nodes_print ( node_num, 
     &  nabes_first, nabes_num, nabes_dim, nabes )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_NEIGHBOR_NODES_PRINT prints a node neighbor array.
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer NABES_FIRST(NODE_NUM), the index in NABES 
c    of the first neighbor in the list for each node.
c
c    Input, integer NABES_NUM(NODE_NUM), the number of neighbors 
c    of each node.
c
c    Input, integer NABES_DIM, the dimension of NABES.
c
c    Input, integer NABES(NABES_DIM), a list of the neighbors 
c    of all the nodes.  Neighbors of node 1 are listed first, and so on.
c
      implicit none

      integer nabes_dim
      integer node_num

      integer i
      integer nabes(nabes_dim)
      integer nabes_first(node_num)
      integer nabes_num(node_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Node based arrays:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Node  Neighbors  Index #1'
      write ( *, '(a)' ) ' '
      do i = 1, node_num
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) 
     &    i, nabes_num(i), nabes_first(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The raw neighbor array.'
      write ( *, '(a)' ) ' '
      do i = 1, nabes_dim
        write ( *, '(2x,i8,2x,i8)' ) i, nabes(i)
      end do

      return
      end
      subroutine triangulation_order3_neighbor_triangles ( 
     &  triangle_num, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
c
c  Discussion:
c
c    A triangulation of a set of nodes can be completely described by
c    the coordinates of the nodes, and the list of nodes that make up
c    each triangle.  However, in some cases, it is necessary to know
c    triangle adjacency information, that is, which triangle, if any,
c    is adjacent to a given triangle on a particular side.
c
c    This routine creates a data structure recording this information.
c
c    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
c    data items.
c
c    Note that COL is a work array allocated dynamically inside this
c    routine.  It is possible, for very large values of TRIANGLE_NUM,
c    that the necessary amount of memory will not be accessible, and the
c    routine will fail.  This is a limitation of the implementation of
c    dynamic arrays in FORTRAN.  One way to get around this would be
c    to require the user to declare COL in the calling routine
c    as an allocatable array, get the necessary memory explicitly with
c    an ALLOCATE statement, and then pass COL into this routine.
c
c    Of course, the point of dynamic arrays was to make it easy to
c    hide these sorts of temporary work arrays from the poor user!
c
c  Example:
c
c    The input information from TRIANGLE_NODE:
c
c    Triangle   Nodes
c    --------   ---------------
c     1         3      4      1
c     2         3      1      2
c     3         3      2      8
c     4         2      1      5
c     5         8      2     13
c     6         8     13      9
c     7         3      8      9
c     8        13      2      5
c     9         9     13      7
c    10         7     13      5
c    11         6      7      5
c    12         9      7      6
c    13        10      9      6
c    14         6      5     12
c    15        11      6     12
c    16        10      6     11
c
c    The output information in TRIANGLE_NEIGHBOR:
c
c    Triangle  Neighboring Triangles
c    --------  ---------------------
c
c     1        -1     -1      2
c     2         1      4      3
c     3         2      5      7
c     4         2     -1      8
c     5         3      8      6
c     6         5      9      7
c     7         3      6     -1
c     8         5      4     10
c     9         6     10     12
c    10         9      8     11
c    11        12     10     14
c    12         9     11     13
c    13        -1     12     16
c    14        11     -1     15
c    15        16     14     -1
c    16        13     15     -1
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangles that are direct neighbors of a given triangle.  
c    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
c    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative
c    if there is no neighbor on that side.  In this case, that side of the 
c    triangle lies on the boundary of the triangulation.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer col(4,3*triangle_num)
      integer i
      integer icol
      integer j
      integer k
      integer side1
      integer side2
      integer triangle_neighbor(3,triangle_num)
      integer tri
      integer triangle_node(triangle_order,triangle_num)
      integer tri1
      integer tri2
c
c  Step 1.
c  From the list of nodes for triangle T, of the form: (I,J,K)
c  construct the three neighbor relations:
c
c    (I,J,1,T) or (J,I,1,T),
c    (J,K,2,T) or (K,J,2,T),
c    (K,I,3,T) or (I,K,3,T)
c
c  where we choose (I,J,1,T) if I .lt. J, or else (J,I,1,T)
c
      do tri = 1, triangle_num

        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)

        if ( i .lt. j ) then
          col(1,3*(tri-1)+1) = i
          col(2,3*(tri-1)+1) = j
          col(3,3*(tri-1)+1) = 1
          col(4,3*(tri-1)+1) = tri
        else
          col(1,3*(tri-1)+1) = j
          col(2,3*(tri-1)+1) = i
          col(3,3*(tri-1)+1) = 1
          col(4,3*(tri-1)+1) = tri
        end if

        if ( j .lt. k ) then
          col(1,3*(tri-1)+2) = j
          col(2,3*(tri-1)+2) = k
          col(3,3*(tri-1)+2) = 2
          col(4,3*(tri-1)+2) = tri
        else
          col(1,3*(tri-1)+2) = k
          col(2,3*(tri-1)+2) = j
          col(3,3*(tri-1)+2) = 2
          col(4,3*(tri-1)+2) = tri
        end if

        if ( k .lt. i ) then
          col(1,3*(tri-1)+3) = k
          col(2,3*(tri-1)+3) = i
          col(3,3*(tri-1)+3) = 3
          col(4,3*(tri-1)+3) = tri
        else
          col(1,3*(tri-1)+3) = i
          col(2,3*(tri-1)+3) = k
          col(3,3*(tri-1)+3) = 3
          col(4,3*(tri-1)+3) = tri
        end if

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c  We only intend to sort on columns 1 and 2; the routine we call here
c  sorts on columns 1 through 4 but that won't hurt us.
c
c  What we need is to find cases where two triangles share an edge.
c  Say they share an edge defined by the nodes I and J.  Then there are
c  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
c  we make sure that these two columns occur consecutively.  That will
c  make it easy to notice that the triangles are neighbors.
c
      call i4col_sort_a ( 4, 3*triangle_num, col )
c
c  Step 3. Neighboring triangles show up as consecutive columns with
c  identical first two entries.  Whenever you spot this happening,
c  make the appropriate entries in TRIANGLE_NEIGHBOR.
c
      do i = 1, 3
        do j = 1, triangle_num
          triangle_neighbor(i,j) = -1
        end do
      end do

      icol = 1

10    continue

        if ( 3 * triangle_num .le. icol ) then
          go to 20
        end if

        if ( col(1,icol) .ne. col(1,icol+1) .or. 
     &       col(2,icol) .ne. col(2,icol+1) ) then
          icol = icol + 1
          go to 10
        end if

        side1 = col(3,icol)
        tri1 = col(4,icol)
        side2 = col(3,icol+1)
        tri2 = col(4,icol+1)

        triangle_neighbor(side1,tri1) = tri2
        triangle_neighbor(side2,tri2) = tri1

        icol = icol + 2

      go to 10

20    continue

      return
      end
      subroutine triangulation_order3_plot ( file_name, node_num, 
     &  node_xy, triangle_num, triangle_node, node_show, 
     &  triangle_show )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_PLOT plots a 3-node triangulation of a set of nodes.
c
c  Discussion:
c
c    The triangulation is most usually a Delaunay triangulation,
c    but this is not necessary.
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
c  Parameters:
c
c    Input, character * ( * ) FILE_NAME, the name of the output file.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists, for 
c    each triangle, the indices of the nodes that form the vertices of the 
c    triangle.
c
c    Input, integer NODE_SHOW,
c    0, do not show nodes;
c    1, show nodes;
c    2, show nodes and label them.
c
c    Input, integer TRIANGLE_SHOW,
c    0, do not show triangles;
c    1, show triangles;
c    2, show triangles and label them.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision ave_x
      double precision ave_y
      integer circle_size
      integer delta
      integer e
      character * ( * )  file_name
      integer file_unit
      integer i
      integer i4_wrap
      integer ios
      integer node
      integer node_show
      double precision node_xy(2,node_num)
      character * ( 40 ) string
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_show
      double precision x_max
      double precision x_min
      integer x_ps
      integer x_ps_max
      integer x_ps_max_clip
      integer x_ps_min
      integer x_ps_min_clip
      double precision x_scale
      double precision y_max
      double precision y_min
      integer y_ps
      integer y_ps_max
      integer y_ps_max_clip
      integer y_ps_min
      integer y_ps_min_clip
      double precision y_scale
c
c  We need to do some figuring here, so that we can determine
c  the range of the data, and hence the height and width
c  of the piece of paper.
c
      x_max = node_xy(1,1)
      x_min = node_xy(1,1)
      do node = 2, node_num
        x_max = max ( x_max, node_xy(1,node) )
        x_min = min ( x_min, node_xy(1,node) )
      end do
      x_scale = x_max - x_min

      x_max = x_max + 0.05D+00 * x_scale
      x_min = x_min - 0.05D+00 * x_scale
      x_scale = x_max - x_min

      y_max = node_xy(2,1)
      y_min = node_xy(2,1)
      do node = 2, node_num
        y_max = max ( y_max, node_xy(2,node) )
        y_min = min ( y_min, node_xy(2,node) )
      end do
      y_scale = y_max - y_min

      y_max = y_max + 0.05D+00 * y_scale
      y_min = y_min - 0.05D+00 * y_scale
      y_scale = y_max - y_min

      x_ps_max = 576
      x_ps_max_clip = 594
      x_ps_min = 36
      x_ps_min_clip = 18
      y_ps_max = 666
      y_ps_max_clip = 684
      y_ps_min = 126
      y_ps_min_clip = 108
      
      if ( x_scale .lt. y_scale ) then

        delta = nint ( dble ( x_ps_max - x_ps_min ) 
     &    * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

        x_ps_max = x_ps_max - delta
        x_ps_min = x_ps_min + delta

        x_ps_max_clip = x_ps_max_clip - delta
        x_ps_min_clip = x_ps_min_clip + delta

        x_scale = y_scale

      else if ( y_scale .lt. x_scale ) then

        delta = nint ( dble ( y_ps_max - y_ps_min ) 
     &    * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

        y_ps_max      = y_ps_max - delta
        y_ps_min      = y_ps_min + delta

        y_ps_max_clip = y_ps_max_clip - delta
        y_ps_min_clip = y_ps_min_clip + delta

        y_scale = x_scale

      end if

      call get_unit ( file_unit )

      open ( unit = file_unit, file = file_name, status = 'replace' )

      write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
      write ( file_unit, '(a)' ) 
     &  '%%Creator: triangulation_order3_plot.f'
      write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
      write ( file_unit, '(a)' ) '%%Pages: 1'
      write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ',
     &  x_ps_min, y_ps_min, x_ps_max, y_ps_max
      write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
      write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
      write ( file_unit, '(a)' ) '%%EndComments'
      write ( file_unit, '(a)' ) '%%BeginProlog'
      write ( file_unit, '(a)' ) '/inch {72 mul} def'
      write ( file_unit, '(a)' ) '%%EndProlog'
      write ( file_unit, '(a)' ) '%%Page: 1 1'
      write ( file_unit, '(a)' ) 'save'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Increase line width from default 0.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '2 setlinewidth'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Set the RGB line color to very light gray.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Draw a gray border around the page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_min, ' moveto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_max, y_ps_min, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_max, y_ps_max, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_max, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_min, ' lineto'
      write ( file_unit, '(a)' ) 'stroke'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the font and its size.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '/Times-Roman findfont'
      write ( file_unit, '(a)' ) '0.50 inch scalefont'
      write ( file_unit, '(a)' ) 'setfont'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Print a title.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  210  702  moveto'
      write ( file_unit, '(a)' ) '%  (Triangulation)  show'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' moveto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a)' ) 'clip newpath'
c
c  Draw the nodes.
c
      if ( node_num .le. 200 ) then
        circle_size = 5
      else if ( node_num .le. 500 ) then
        circle_size = 4
      else if ( node_num .le. 1000 ) then
        circle_size = 3
      else if ( node_num .le. 5000 ) then
        circle_size = 2
      else
        circle_size = 1
      end if

      if ( 1 .le. node_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (         node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 
     &      'newpath ', x_ps, y_ps, 
     &      circle_size, '0 360 arc closepath fill'

        end do

      end if
c
c  Label the nodes.
c
      if ( 2 .le. node_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the nodes:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) 
     &    '%  Set the RGB color to darker blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (       + node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( string, '(i4)' ) node
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, 
     &      ' moveto (' // trim ( string ) // ') show'

        end do

      end if
c
c  Draw the triangles.
c
      if ( 1 .le. triangle_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw the triangles.'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          write ( file_unit, '(a)' ) 'newpath'

          do i = 1, 4

            e = i4_wrap ( i, 1, 3 )

            node = triangle_node(e,triangle)

            x_ps = int ( 
     &        ( ( x_max - node_xy(1,node)         ) 
     &        * dble ( x_ps_min )   
     &        + (         node_xy(1,node) - x_min ) 
     &        * dble ( x_ps_max ) )
     &        / ( x_max                   - x_min ) )

            y_ps = int ( 
     &        ( ( y_max - node_xy(2,node)         ) 
     &        * dble ( y_ps_min )   
     &        + (         node_xy(2,node) - y_min ) 
     &        * dble ( y_ps_max ) )
     &        / ( y_max                   - y_min ) )

            if ( i .eq. 1 ) then
              write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &          x_ps, y_ps, ' moveto'
            else
              write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &          x_ps, y_ps, ' lineto'
            end if

          end do

          write ( file_unit, '(a)' ) 'stroke'

        end do

      end if
c
c  Label the triangles.
c
      if ( 2 .le. triangle_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the triangles:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          ave_x = 0.0D+00
          ave_y = 0.0D+00

          do i = 1, 3

            node = triangle_node(i,triangle)

            ave_x = ave_x + node_xy(1,node)
            ave_y = ave_y + node_xy(2,node)

          end do

          ave_x = ave_x / 3.0D+00
          ave_y = ave_y / 3.0D+00

          x_ps = int ( 
     &      ( ( x_max - ave_x         ) * dble ( x_ps_min )   
     &      + (       + ave_x - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max         - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - ave_y         ) * dble ( y_ps_min )   
     &      + (         ave_y - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max         - y_min ) )

          write ( string, '(i4)' ) triangle
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' 
     &      // trim ( string ) // ') show'

        end do

      end if

      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'restore  showpage'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  End of page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%%Trailer'
      write ( file_unit, '(a)' ) '%%EOF'
      close ( unit = file_unit )

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
      subroutine triangulation_order3_quad ( node_num, node_xy, 
     &  triangle_order, triangle_num, triangle_node, quad_fun, 
     &  quad_num, quad_xy, quad_w, quad_value, region_area )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_QUAD approximates an integral over a triangulation.
c
c  Discussion:
c
c    The routine will accept triangulations of order higher than 3.
c    However, only the first three nodes (the vertices) of each
c    triangle will be used.  This will still produce correct results
c    for higher order triangulations, as long as the sides of the
c    triangle are straight.
c
c    We assume that the vertices of each triangle are listed first
c    in the description of higher order triangles, and we assume that
c    the vertices are listed in counterclockwise order.
c
c    The approximation of the integral is made using a quadrature rule 
c    defined on the unit triangle, and supplied by the user.  
c
c    The user also supplies the name of a subroutine, here called "QUAD_FUN", 
c    which evaluates the integrand at a set of points.  The form of
c    this routine is:
c
c      subroutine quad_fun ( n, xy_vec, f_vec )
c      integer n
c      double precision f_vec(n)
c      double precision xy_vec(2,n)
c
c    and it returns in each entry F_VEC(1:N), the value of the integrand
c    at XY_VEC(1:2,1:N).
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes in the 
c    triangulation.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_ORDER, the order of triangles in 
c    the triangulation.
c
c    Input, integer TRIANGLE_NUM, the number of triangles in 
c    the triangulation.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
c    the nodes making up each triangle.
c
c    Input, external QUAD_FUN, the name of the integrand routine.
c
c    Input, integer QUAD_NUM, the order of the quadrature rule.
c
c    Input, double precision QUAD_XY(2,QUAD_NUM), the abscissas of the 
c    quadrature rule, in the unit triangle.
c
c    Input, double precision QUAD_W(QUAD_NUM), the weights of the 
c    quadrature rule.
c
c    Output, double precision QUAD_VALUE, the estimate of the integral
c    of F(X,Y) over the region covered by the triangulation.
c
c    Output, double precision REGION_AREA, the area of the region.
c
      implicit none

      integer node_num
      integer quad_num
      integer triangle_num
      integer triangle_order

      integer i
      integer j
      double precision node_xy(2,node_num)
      double precision quad_f(quad_num)
      external quad_fun
      double precision quad_value
      double precision quad_w(quad_num)
      double precision quad_xy(2,quad_num)
      double precision quad2_xy(2,quad_num)
      double precision r8vec_dot_product
      double precision region_area
      integer triangle
      double precision triangle_area
      integer triangle_node(triangle_order,triangle_num)
      double precision triangle_xy(2,3)

      quad_value = 0.0D+00
      region_area = 0.0D+00

      do triangle = 1, triangle_num

        do j = 1, 3
          do i = 1, 2
            triangle_xy(i,j) = node_xy(i,triangle_node(j,triangle))
          end do
        end do

        call triangle_area_2d ( triangle_xy, triangle_area )

        call triangle_order3_reference_to_physical ( triangle_xy, 
     &    quad_num, quad_xy, quad2_xy )

        call quad_fun ( quad_num, quad2_xy, quad_f )

        quad_value = quad_value + triangle_area 
     &    * r8vec_dot_product ( quad_num, quad_w, quad_f )

        region_area = region_area + triangle_area

      end do

      return
      end
      subroutine triangulation_order3_refine_compute ( node_num1, 
     &  triangle_num1, node_xy1, triangle_node1, node_num2, 
     &  triangle_num2, edge_data, node_xy2, triangle_node2 )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
c
c  Discussion:
c
c    Given a triangle defined by nodes 1, 2, 3, we need to generate
c    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
c    and T4.
c
c    The task is more complicated by the fact that we are working with
c    a mesh of triangles, so that we want to create a node only once,
c    even though it may be shared by other triangles.
c
c          3
c         / \
c        /T3 \
c      13----23
c      / \T4 / \
c     /T1 \ /T2 \
c    1----12-----2
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
c    Input, integer NODE_NUM1, the number of nodes.
c
c    Input, integer TRIANGLE_NUM1, the number of triangles.
c
c    Input, double precision NODE_XY1(2,NODE_NUM1), the nodes.
c
c    Input, integer TRIANGLE_NODE1(3,TRIANGLE_NUM1), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise 
c    order.
c
c    Input, integer NODE_NUM2, the number of nodes in the 
c    refined mesh.
c
c    Input, integer TRIANGLE_NUM2, the number of triangles in 
c    the refined mesh.
c
c    Input, integer EDGE_DATA(5,3*TRIANGLE_NUM1), edge information
c    computed by TRIANGULATION_ORDER3_REFINE_SIZE.
c
c    Output, double precision NODE_XY2(2,NODE_NUM2), the refined nodes.
c
c    Output, integer TRIANGLE_NODE2(3,TRIANGLE_NUM2), the nodes 
c    that make up the triangles in the refined mesh.
c
      implicit none

      integer node_num1
      integer node_num2
      integer triangle_num1
      integer triangle_num2

      integer edge
      integer edge_data(5,3*triangle_num1)
      integer i
      integer j
      integer n1
      integer n1_old
      integer n2
      integer n2_old
      integer node
      double precision node_xy1(2,node_num1)
      double precision node_xy2(2,node_num2)
      integer triangle_node1(3,triangle_num1)
      integer triangle_node2(3,triangle_num2)
      integer triangle1
      integer v1
      integer v2
c
c  Copy the old nodes.
c
      do j = 1, node_num1
        do i = 1, 2
          node_xy2(i,j) = node_xy1(i,j)
        end do
      end do

      do j = 1, triangle_num2
        do i = 1, 3
          triangle_node2(i,j) = -1
        end do
      end do
c
c  We can assign the existing nodes to the new triangles.
c
      do triangle1 = 1, triangle_num1
        triangle_node2(1,(triangle1-1)*4+1) 
     &    = triangle_node1(1,triangle1)
        triangle_node2(2,(triangle1-1)*4+2) 
     &    = triangle_node1(2,triangle1)
        triangle_node2(3,(triangle1-1)*4+3) 
     &    = triangle_node1(3,triangle1)
      end do

      node = node_num1

      n1_old = -1
      n2_old = -1

      do edge = 1, 3 * triangle_num1

        n1 = edge_data(1,edge)
        n2 = edge_data(2,edge)
c
c  If this edge is new, create the coordinates and index for this node.
c
        if ( n1 .ne. n1_old .or. n2 .ne. n2_old ) then
          node = node + 1

          if ( node_num2 .lt. node ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 
     &        'TRIANGLE_MESH_ORDER3_REFINE - Fatal error!'
            write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
            stop
          end if

          node_xy2(1,node) = 
     &      ( node_xy2(1,n1) + node_xy2(1,n2) ) / 2.0D+00
          node_xy2(2,node) = 
     &      ( node_xy2(2,n1) + node_xy2(2,n2) ) / 2.0D+00

          n1_old = n1
          n2_old = n2

        end if
c
c  Assign the node to triangles.
c
        v1 = edge_data(3,edge)
        v2 = edge_data(4,edge)
        triangle1 = edge_data(5,edge)

        if ( v1 .eq. 1 .and. v2 .eq. 2 ) then

          triangle_node2(1,(triangle1-1)*4+2) = node
          triangle_node2(2,(triangle1-1)*4+1) = node
          triangle_node2(3,(triangle1-1)*4+4) = node

        else if ( v1 .eq. 1 .and. v2 .eq. 3 ) then

          triangle_node2(1,(triangle1-1)*4+3) = node
          triangle_node2(2,(triangle1-1)*4+4) = node
          triangle_node2(3,(triangle1-1)*4+1) = node

        else if ( v1 .eq. 2 .and. v2 .eq. 3 ) then

          triangle_node2(1,(triangle1-1)*4+4) = node
          triangle_node2(2,(triangle1-1)*4+3) = node
          triangle_node2(3,(triangle1-1)*4+2) = node

        end if

      end do

      return
      end
      subroutine triangulation_order3_refine_size ( node_num1, 
     &  triangle_num1, triangle_node1, node_num2, triangle_num2, 
     &  edge_data )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
c
c  Discussion:
c
c    Given a triangle defined by nodes 1, 2, 3, we need to generate
c    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
c    and T4.
c
c    The task is more complicated by the fact that we are working with
c    a mesh of triangles, so that we want to create a node only once,
c    even though it may be shared by other triangles.
c
c          3
c         / \
c        /T3 \
c      13----23
c      / \T4 / \
c     /T1 \ /T2 \
c    1----12-----2
c
c    This routine simply determines the sizes of the resulting node 
c    and triangle arrays.
c
c    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
c    data items, one item for every edge of every triangle.  Each
c    data item records, for a given edge, the global indices
c    of the two endpoints, the local indices of the two endpoints,
c    and the index of the triangle.
c
c    Through careful sorting, it is possible to arrange this data in
c    a way that allows the proper generation of the interpolated nodes.
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
c    Input, integer NODE_NUM1, the number of nodes in the 
c    original mesh.
c
c    Input, integer TRIANGLE_NUM1, the number of triangles in the
c    original mesh.
c
c    Input, integer TRIANGLE_NODE1(3,TRIANGLE_NUM1), the indices 
c    of the nodes that form the triangles in the input mesh.
c
c    Output, integer NODE_NUM2, the number of nodes in the refined
c    mesh.
c
c    Output, integer TRIANGLE_NUM2, the number of triangles in 
c    the refined mesh.
c
c    Output, integer EDGE_DATA(5,3*TRIANGLE_NUM1), edge data that 
c    will be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
c
      implicit none

      integer node_num1
      integer triangle_num1

      integer a
      integer b
      integer edge
      integer edge_data(5,3*triangle_num1)
      integer i
      integer j
      integer k
      integer n1
      integer n1_old
      integer n2
      integer n2_old
      integer node_num2
      integer triangle
      integer triangle_node1(3,triangle_num1)
      integer triangle_num2
c
c  Step 1.
c  From the list of nodes for triangle T, of the form: (I,J,K)
c  construct the edge relations:
c
c    (I,J,1,2,T)
c    (I,K,1,3,T)
c    (J,K,2,3,T)
c
c  In order to make matching easier, we reorder each pair of nodes
c  into ascending order.
c
      do triangle = 1, triangle_num1

        i = triangle_node1(1,triangle)
        j = triangle_node1(2,triangle)
        k = triangle_node1(3,triangle)

        call i4i4_sort_a ( i, j, a, b )

        edge_data(1,3*(triangle-1)+1) = a
        edge_data(2,3*(triangle-1)+1) = b
        edge_data(3,3*(triangle-1)+1) = 1
        edge_data(4,3*(triangle-1)+1) = 2
        edge_data(5,3*(triangle-1)+1) = triangle

        call i4i4_sort_a ( i, k, a, b )

        edge_data(1,3*(triangle-1)+2) = a
        edge_data(2,3*(triangle-1)+2) = b
        edge_data(3,3*(triangle-1)+2) = 1
        edge_data(4,3*(triangle-1)+2) = 3
        edge_data(5,3*(triangle-1)+2) = triangle

        call i4i4_sort_a ( j, k, a, b )

        edge_data(1,3*(triangle-1)+3) = a
        edge_data(2,3*(triangle-1)+3) = b
        edge_data(3,3*(triangle-1)+3) = 2
        edge_data(4,3*(triangle-1)+3) = 3
        edge_data(5,3*(triangle-1)+3) = triangle

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c  We only intend to sort on rows 1:2; the routine we call here
c  sorts on the full column but that won't hurt us.
c
c  What we need is to find all cases where triangles share an edge.
c  By sorting the columns of the EDGE_DATA array, we will put shared edges
c  next to each other.
c
      call i4col_sort_a ( 5, 3*triangle_num1, edge_data )
c
c  Step 3. All the triangles which share an edge show up as consecutive
c  columns with identical first two entries.  Figure out how many new
c  nodes there are, and allocate space for their coordinates.
c
      node_num2 = node_num1

      n1_old = -1
      n2_old = -1

      do edge = 1, 3 * triangle_num1
        n1 = edge_data(1,edge)
        n2 = edge_data(2,edge)
        if ( n1 .ne. n1_old .or. n2 .ne. n2_old ) then
          node_num2 = node_num2 + 1
          n1_old = n1
          n2_old = n2
        end if
      end do

      triangle_num2 = 4 * triangle_num1

      return
      end
      subroutine triangulation_order3_sample ( node_num, node_xy, 
     &  triangle_num, triangle_node, num_ran, seed, xd, td )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_SAMPLE returns random points in a triangulation.
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
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes that 
c    make up the triangles.
c
c    Input, integer NUM_RAN, the number of points to sample.
c
c    Input/output, integer SEED, a seed for the random 
c     number generator.
c
c    Output, double precision XD(2,NUM_RAN), the sample points.
c
c    Output, integer TD(NUM_RAN), the triangle to which each 
c    sample point belongs.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer num_ran
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      double precision area
      double precision area_cum(0:triangle_num)
      double precision area_total
      integer i
      integer i1
      integer i2
      integer i3
      integer left
      double precision node_xy(dim_num,node_num)
      double precision r
      double precision r8_uniform_01
      integer right
      integer seed
      double precision t(dim_num,3)
      integer td(num_ran)
      integer triangle_node(triangle_order,triangle_num)
      double precision xd(dim_num,num_ran)
c
c  Compute the areas of the triangles.
c  Build a cumulative area vector.
c  Convert it to a relative cumulative area vector.
c
      area_cum(0) = 0.0D+00

      do i = 1, triangle_num

        i1 = triangle_node(1,i)
        i2 = triangle_node(2,i)
        i3 = triangle_node(3,i)

        t(1,1) = node_xy(1,i1)
        t(2,1) = node_xy(2,i1)
        t(1,2) = node_xy(1,i2)
        t(2,2) = node_xy(2,i2)
        t(1,3) = node_xy(1,i3)
        t(2,3) = node_xy(2,i3)

        call triangle_area_2d ( t, area )

        area_cum(i) = area_cum(i-1) + area

      end do

      area_total = area_cum(triangle_num)

      do i = 0, triangle_num
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

        call r8vec_bracket ( triangle_num+1, area_cum, r, left, right )

        td(i) = right - 1

        i1 = triangle_node(1,td(i))
        i2 = triangle_node(2,td(i))
        i3 = triangle_node(3,td(i))

        t(1,1) = node_xy(1,i1)
        t(2,1) = node_xy(2,i1)
        t(1,2) = node_xy(1,i2)
        t(2,2) = node_xy(2,i2)
        t(1,3) = node_xy(1,i3)
        t(2,3) = node_xy(2,i3)

        call triangle_sample ( t, 1, seed, xd(1,i) )

      end do

      return
      end
      subroutine triangulation_order4_plot ( file_name, node_num, 
     &  node_xy, triangle_num, triangle_node, node_show, 
     &  triangle_show )

c*********************************************************************72
c
cc TRIANGULATION_ORDER4_PLOT plots a 4-node triangulation of a pointset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_NAME, the name of the output file.
c
c    Input, integer NODE_NUM, the number of points.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(4,TRIANGLE_NUM), lists, 
c    for each triangle, the indices of the points that form the vertices 
c    and the centroid of the triangle.
c
c    Input, integer NODE_SHOW,
c    0, do not show nodes;
c    1, show nodes;
c    2, show nodes and label them.
c
c    Input, integer TRIANGLE_SHOW,
c    0, do not show triangles;
c    1, show triangles;
c    2, show triangles and label them.
c
      implicit none

      integer node_num
      integer triangle_num

      double precision ave_x
      double precision ave_y
      integer circle_size
      integer delta
      integer e
      character * ( * )  file_name
      integer file_unit
      integer i
      integer i4_wrap
      integer ios
      integer node
      integer node_show
      double precision node_xy(2,node_num)
      character ( len = 40 ) string
      integer triangle
      integer triangle_node(4,triangle_num)
      integer triangle_show
      double precision x_max
      double precision x_min
      integer x_ps
      integer x_ps_max
      integer x_ps_max_clip
      integer x_ps_min
      integer x_ps_min_clip
      double precision x_scale
      double precision y_max
      double precision y_min
      integer y_ps
      integer y_ps_max
      integer y_ps_max_clip
      integer y_ps_min
      integer y_ps_min_clip
      double precision y_scale

      x_ps_max = 576
      x_ps_max_clip = 594
      x_ps_min = 36
      x_ps_min_clip = 18
      y_ps_max = 666
      y_ps_max_clip = 684
      y_ps_min = 126
      y_ps_min_clip = 108
c
c  We need to do some figuring here, so that we can determine
c  the range of the data, and hence the height and width
c  of the piece of paper.
c
      x_max = node_xy(1,1)
      x_min = node_xy(1,1)
      do node = 2, node_num
        x_max = max ( x_max, node_xy(1,node) )
        x_min = min ( x_min, node_xy(1,node) )
      end do
      x_scale = x_max - x_min

      x_max = x_max + 0.05D+00 * x_scale
      x_min = x_min - 0.05D+00 * x_scale
      x_scale = x_max - x_min

      y_max = node_xy(2,1)
      x_min = node_xy(2,1)
      do node = 2, node_num
        y_max = max ( y_max, node_xy(2,node) )
        y_min = min ( y_min, node_xy(2,node) )
      end do
      y_scale = y_max - y_min

      y_max = y_max + 0.05D+00 * y_scale
      y_min = y_min - 0.05D+00 * y_scale
      y_scale = y_max - y_min

      if ( x_scale .lt. y_scale ) then

        delta = nint ( dble ( x_ps_max - x_ps_min ) 
     &    * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

        x_ps_max = x_ps_max - delta
        x_ps_min = x_ps_min + delta

        x_ps_max_clip = x_ps_max_clip - delta
        x_ps_min_clip = x_ps_min_clip + delta

        x_scale = y_scale

      else if ( y_scale .lt. x_scale ) then

        delta = nint ( dble ( y_ps_max - y_ps_min ) 
     &    * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

        y_ps_max      = y_ps_max - delta
        y_ps_min      = y_ps_min + delta

        y_ps_max_clip = y_ps_max_clip - delta
        y_ps_min_clip = y_ps_min_clip + delta

        y_scale = x_scale

      end if

      call get_unit ( file_unit )

      open ( unit = file_unit, file = file_name, status = 'replace' )

      write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
      write ( file_unit, '(a)' ) 
     &  '%%Creator: triangulation_order4_plot.f'
      write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
      write ( file_unit, '(a)' ) '%%Pages: 1'
      write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) 
     &  '%%BoundingBox: ', 
     &  x_ps_min, y_ps_min, x_ps_max, y_ps_max
      write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
      write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
      write ( file_unit, '(a)' ) '%%EndComments'
      write ( file_unit, '(a)' ) '%%BeginProlog'
      write ( file_unit, '(a)' ) '/inch {72 mul} def'
      write ( file_unit, '(a)' ) '%%EndProlog'
      write ( file_unit, '(a)' ) '%%Page: 1 1'
      write ( file_unit, '(a)' ) 'save'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Increase line width from default 0.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '2 setlinewidth'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Set the RGB line color to very light gray.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Draw a gray border around the page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(2x,i3,2x,i3,2x,a)' ) 
     &  x_ps_min, y_ps_min, ' moveto'
      write ( file_unit, '(2x,i3,2x,i3,2x,a)' ) 
     &  x_ps_max, y_ps_min, ' lineto'
      write ( file_unit, '(2x,i3,2x,i3,2x,a)' ) 
     &  x_ps_max, y_ps_max, ' lineto'
      write ( file_unit, '(2x,i3,2x,i3,2x,a)' ) 
     &  x_ps_min, y_ps_max, ' lineto'
      write ( file_unit, '(2x,i3,2x,i3,2x,a)' ) 
     &  x_ps_min, y_ps_min, ' lineto'
      write ( file_unit, '(a)' ) 'stroke'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the font and its size.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '/Times-Roman findfont'
      write ( file_unit, '(a)' ) '0.50 inch scalefont'
      write ( file_unit, '(a)' ) 'setfont'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Print a title.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  210  702  moveto'
      write ( file_unit, '(a)' ) '%  (Triangulation)  show'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' moveto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a)' ) 'clip newpath'
c
c  Draw the nodes.
c
      if ( node_num .le. 200 ) then
        circle_size = 5
      else if ( node_num .le. 500 ) then
        circle_size = 4
      else if ( node_num .le. 1000 ) then
        circle_size = 3
      else if ( node_num .le. 5000 ) then
        circle_size = 2
      else
        circle_size = 1
      end if

      if ( 1 .le. node_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (         node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 
     &      'newpath ', x_ps, y_ps, 
     &      circle_size, '0 360 arc closepath fill'

        end do

      end if
c
c  Label the nodes.
c
      if ( 2 .le. node_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the nodes:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) 
     &    '%  Set the RGB color to darker blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (       + node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( string, '(i4)' ) node
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, 
     &      ' moveto (' // trim ( string ) // ') show'

        end do

      end if
c
c  Draw the triangles.
c
      if ( 1 .le. triangle_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw the triangles.'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          write ( file_unit, '(a)' ) 'newpath'

          do i = 1, 4

            e = i4_wrap ( i, 1, 3 )

            node = triangle_node(e,triangle)

            x_ps = int ( 
     &        ( ( x_max - node_xy(1,node)         ) 
     &        * dble ( x_ps_min )   
     &        + (         node_xy(1,node) - x_min ) 
     &        * dble ( x_ps_max ) ) 
     &        / ( x_max                   - x_min ) )

            y_ps = int ( 
     &        ( ( y_max - node_xy(2,node)         ) 
     &        * dble ( y_ps_min )   
     &        + (         node_xy(2,node) - y_min ) 
     &        * dble ( y_ps_max ) ) 
     &        / ( y_max                   - y_min ) )

            if ( i .eq. 1 ) then
              write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &          x_ps, y_ps, ' moveto'
            else
              write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &          x_ps, y_ps, ' lineto'
            end if

          end do

          write ( file_unit, '(a)' ) 'stroke'

        end do

      end if
c
c  Label the triangles.
c
      if ( 2 .le. triangle_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the triangles:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          ave_x = 0.0D+00
          ave_y = 0.0D+00

          do i = 1, 3

            node = triangle_node(i,triangle)

            ave_x = ave_x + node_xy(1,node)
            ave_y = ave_y + node_xy(2,node)

          end do

          ave_x = ave_x / 3.0D+00
          ave_y = ave_y / 3.0D+00

          x_ps = int ( 
     &      ( ( x_max - ave_x         ) * dble ( x_ps_min )   
     &      + (       + ave_x - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max         - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - ave_y         ) * dble ( y_ps_min )   
     &      + (         ave_y - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max         - y_min ) )

          write ( string, '(i4)' ) triangle
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' 
     &      // trim ( string ) // ') show'

        end do

      end if

      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'restore  showpage'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  End of page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%%Trailer'
      write ( file_unit, '(a)' ) '%%EOF'
      close ( unit = file_unit )

      return
      end
      subroutine triangulation_order6_adj_count ( node_num, 
     &  triangle_num, triangle_node, triangle_neighbor, adj_num, 
     &  adj_col )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The triangulation is assumed to involve 6-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  6   5  side 2
c       |    \
c    3  |     \
c       |      \
c       1---4---2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\    |\    |
c    | \   | \   |
c   16 17 18 19 20
c    |   \ |   \ |
c    |    \|    \|
c   11-12-13-14-15
c    |\    |\    |
c    | \   | \   |
c    6  7  8  9 10
c    |   \ |   \ |
c    |    \|    \|
c    1--2--3--4--5
c
c    A sample grid.
c
c
c    Below, we have a chart that lists the nodes adjacent to each node, with 
c    an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).
c
c    N   Adjacencies
c
c    1:  *  2  3  6  7 11
c    2:  1  *  3  6  7 11
c    3:  1  2  *  4  5  6  7  8  9 11 12 13
c    4:  3  *  5  8  9 13
c    5:  3  4  *  8  9 10 13 14 15
c    6:  1  2  3  *  7 11
c    7:  1  2  3  6  *  8 11 12 13
c    8:  3  4  5  7  *  9 11 12 13
c    9:  3  4  5  8  * 10 13 14 15
c   10:  5  9  * 13 14 15
c   11:  1  2  3  6  7  8  * 12 13 16 17 21
c   12:  3  7  8 11  * 13 16 17 21
c   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
c   14:  5  9 10 13  * 15 18 19 23
c   15:  5  9 10 13 14  * 18 19 20 23 24 25
c   16: 11 12 13  * 17 21
c   17: 11 12 13 16  * 18 21 22 23
c   18: 13 14 15 17  * 19 21 22 23
c   19: 13 14 15 18  * 20 23 24 25
c   20: 15 19  * 23 24 25
c   21: 11 12 13 16 17 18  * 22 23
c   22: 13 17 18 21  * 23
c   23: 13 14 15 17 18 19 20 21 22  * 24 25
c   24: 15 19 20 23  * 25
c   25: 15 19 20 23 24  *    
c
c    Below, we list the number of adjacencies to lower-indexed nodes, to 
c    the node itself, to higher-indexed nodes, the total number of 
c    adjacencies for this node, and the location of the first and last 
c    entries required to list this set of adjacencies in a single list 
c    of all the adjacencies.
c
c    N   Below  Self   Above   Total First  Last
c
c   --      --    --      --      --   ---     0
c    1:      0     1       5       6     1     6
c    2:      1     1       4       6     7    12
c    3:      2     1       9      12    13    24
c    4:      1     1       4       6    25    30
c    5:      2     1       6       9    31    39
c    6:      3     1       2       6    40    45
c    7:      4     1       4       9    46    54
c    8:      4     1       4       9    55    63
c    9:      4     1       4       9    62    72
c   10:      2     1       3       6    73    78
c   11:      6     1       5      12    79    90
c   12:      4     1       4       9    91    99
c   13:      9     1       9      19   100   118
c   14:      4     1       4       9   119   127
c   15:      5     1       6      12   128   139
c   16:      3     1       2       6   140   145
c   17:      4     1       4       9   146   154
c   18:      4     1       4       9   155   163
c   19:      4     1       4       9   164   172
c   20:      2     1       3       6   173   178
c   21:      6     1       2       9   179   187
c   22:      4     1       1       6   188   193
c   23:      9     1       2      12   194   205
c   24:      4     1       1       6   206   211
c   25:      5     1       0       6   212   217
c   --      --    --      --      --   218   ---
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Output, integer ADJ_NUM, the number of adjacencies.
c
c    Output, integer ADJ_COL(NODE_NUM+1).  Information about 
c    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer adj_num
      integer adj_col(node_num+1)
      integer i
      integer n1
      integer n2
      integer n3
      integer n4
      integer n5
      integer n6
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      adj_num = 0
c
c  Set every node to be adjacent to itself.
c
      adj_col(1:node_num) = 1
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
        n4 = triangle_node(4,triangle)
        n5 = triangle_node(5,triangle)
        n6 = triangle_node(6,triangle)
c
c  For sure, we add the adjacencies:
c    43 / (34)
c    51 / (15)
c    54 / (45)
c    62 / (26)
c    64 / (46)
c    65 / (56)
c
        adj_col(n3) = adj_col(n3) + 1
        adj_col(n4) = adj_col(n4) + 1
        adj_col(n1) = adj_col(n1) + 1
        adj_col(n5) = adj_col(n5) + 1
        adj_col(n4) = adj_col(n4) + 1
        adj_col(n5) = adj_col(n5) + 1
        adj_col(n2) = adj_col(n2) + 1
        adj_col(n6) = adj_col(n6) + 1
        adj_col(n4) = adj_col(n4) + 1
        adj_col(n6) = adj_col(n6) + 1
        adj_col(n5) = adj_col(n5) + 1
        adj_col(n6) = adj_col(n6) + 1
c
c  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
c  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
c  Maybe add
c    21 / 12
c    41 / 14
c    42 / 24
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n4) = adj_col(n4) + 1
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n4) = adj_col(n4) + 1
        end if
c
c  Maybe add
c    32 / 23
c    52 / 25
c    53 / 35
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n3) = adj_col(n3) + 1
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n5) = adj_col(n5) + 1
          adj_col(n3) = adj_col(n3) + 1
          adj_col(n5) = adj_col(n5) + 1
        end if
c
c  Maybe add
c    31 / 13
c    61 / 16
c    63 / 36
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n3) = adj_col(n3) + 1
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n6) = adj_col(n6) + 1
          adj_col(n3) = adj_col(n3) + 1
          adj_col(n6) = adj_col(n6) + 1
        end if
          
      end do
c
c  We used ADJ_COL to count the number of entries in each column.
c  Convert it to pointers into the ADJ array.
c
      do i = node_num, 1, -1
        adj_col(i+1) = adj_col(i)
      end do

      adj_col(1) = 1
      do i = 2, node_num+1
        adj_col(i) = adj_col(i-1) + adj_col(i)
      end do

      adj_num = adj_col(node_num+1) - 1

      return
      end
      subroutine triangulation_order6_adj_set ( node_num, 
     &  triangle_num, triangle_node, triangle_neighbor, adj_num, 
     &  adj_col, adj )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The triangulation is assumed to involve 6-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c    This routine can be used to create the compressed column storage
c    for a quadratic triangle finite element discretization of 
c    Poisson's equation in two dimensions.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  6   5  side 2
c       |    \
c    3  |     \
c       |      \
c       1---4---2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\    |\    |
c    | \   | \   |
c   16 17 18 19 20
c    |   \ |   \ |
c    |    \|    \|
c   11-12-13-14-15
c    |\    |\    |
c    | \   | \   |
c    6  7  8  9 10
c    |   \ |   \ |
c    |    \|    \|
c    1--2--3--4--5
c
c    A sample grid.
c
c
c    Below, we have a chart that lists the nodes adjacent to each node, with 
c    an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).
c
c    N   Adjacencies
c
c    1:  *  2  3  6  7 11
c    2:  1  *  3  6  7 11
c    3:  1  2  *  4  5  6  7  8  9 11 12 13
c    4:  3  *  5  8  9 13
c    5:  3  4  *  8  9 10 13 14 15
c    6:  1  2  3  *  7 11
c    7:  1  2  3  6  *  8 11 12 13
c    8:  3  4  5  7  *  9 11 12 13
c    9:  3  4  5  8  * 10 13 14 15
c   10:  5  9  * 13 14 15
c   11:  1  2  3  6  7  8  * 12 13 16 17 21
c   12:  3  7  8 11  * 13 16 17 21
c   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
c   14:  5  9 10 13  * 15 18 19 23
c   15:  5  9 10 13 14  * 18 19 20 23 24 25
c   16: 11 12 13  * 17 21
c   17: 11 12 13 16  * 18 21 22 23
c   18: 13 14 15 17  * 19 21 22 23
c   19: 13 14 15 18  * 20 23 24 25
c   20: 15 19  * 23 24 25
c   21: 11 12 13 16 17 18  * 22 23
c   22: 13 17 18 21  * 23
c   23: 13 14 15 17 18 19 20 21 22  * 24 25
c   24: 15 19 20 23  * 25
c   25: 15 19 20 23 24  *    
c
c    Below, we list the number of adjacencies to lower-indexed nodes, to 
c    the node itself, to higher-indexed nodes, the total number of 
c    adjacencies for this node, and the location of the first and last 
c    entries required to list this set of adjacencies in a single list 
c    of all the adjacencies.
c
c    N   Below  Self   Above   Total First  Last
c
c   --      --    --      --      --   ---     0
c    1:      0     1       5       6     1     6
c    2:      1     1       4       6     7    12
c    3:      2     1       9      12    13    24
c    4:      1     1       4       6    25    30
c    5:      2     1       6       9    31    39
c    6:      3     1       2       6    40    45
c    7:      4     1       4       9    46    54
c    8:      4     1       4       9    55    63
c    9:      4     1       4       9    62    72
c   10:      2     1       3       6    73    78
c   11:      6     1       5      12    79    90
c   12:      4     1       4       9    91    99
c   13:      9     1       9      19   100   118
c   14:      4     1       4       9   119   127
c   15:      5     1       6      12   128   139
c   16:      3     1       2       6   140   145
c   17:      4     1       4       9   146   154
c   18:      4     1       4       9   155   163
c   19:      4     1       4       9   164   172
c   20:      2     1       3       6   173   178
c   21:      6     1       2       9   179   187
c   22:      4     1       1       6   188   193
c   23:      9     1       2      12   194   205
c   24:      4     1       1       6   206   211
c   25:      5     1       0       6   212   217
c   --      --    --      --      --   218   ---
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the nodes 
c    that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Input, integer ADJ_NUM, the number of adjacencies.
c
c    Input, integer ADJ_COL(NODE_NUM+1).  Information about column
c    J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
c    Output, integer ADJ(ADJ_NUM), the adjacency information.
c
      implicit none

      integer adj_num
      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer adj(adj_num)
      integer adj_copy(node_num)
      integer adj_col(node_num+1)
      integer i
      integer k1
      integer k2
      integer n1
      integer n2
      integer n3
      integer n4
      integer n5
      integer n6
      integer node
      integer number
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      do i = 1, adj_num
        adj(i) = -1
      end do

      do node = 1, node_num
        adj_copy(node) = adj_col(node)
      end do
c
c  Set every node to be adjacent to itself.
c
      do node = 1, node_num
        adj(adj_copy(node)) = node
        adj_copy(node) = adj_copy(node) + 1
      end do
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
        n4 = triangle_node(4,triangle)
        n5 = triangle_node(5,triangle)
        n6 = triangle_node(6,triangle)
c
c  For sure, we add the adjacencies:
c    43 / (34)
c    51 / (15)
c    54 / (45)
c    62 / (26)
c    64 / (46)
c    65 / (56)
c
        adj(adj_copy(n3)) = n4
        adj_copy(n3) = adj_copy(n3) + 1
        adj(adj_copy(n4)) = n3
        adj_copy(n4) = adj_copy(n4) + 1

        adj(adj_copy(n1)) = n5
        adj_copy(n1) = adj_copy(n1) + 1
        adj(adj_copy(n5)) = n1
        adj_copy(n5) = adj_copy(n5) + 1

        adj(adj_copy(n4)) = n5
        adj_copy(n4) = adj_copy(n4) + 1
        adj(adj_copy(n5)) = n4
        adj_copy(n5) = adj_copy(n5) + 1

        adj(adj_copy(n2)) = n6
        adj_copy(n2) = adj_copy(n2) + 1
        adj(adj_copy(n6)) = n2
        adj_copy(n6) = adj_copy(n6) + 1

        adj(adj_copy(n4)) = n6
        adj_copy(n4) = adj_copy(n4) + 1
        adj(adj_copy(n6)) = n4
        adj_copy(n6) = adj_copy(n6) + 1

        adj(adj_copy(n5)) = n6
        adj_copy(n5) = adj_copy(n5) + 1
        adj(adj_copy(n6)) = n5
        adj_copy(n6) = adj_copy(n6) + 1
c
c  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
c  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE .lt. TRIANGLE2).
c
c  Maybe add
c    21 / 12
c    41 / 14
c    42 / 24
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n1)) = n2
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n2)) = n1
          adj_copy(n2) = adj_copy(n2) + 1
          adj(adj_copy(n1)) = n4
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n4)) = n1
          adj_copy(n4) = adj_copy(n4) + 1
          adj(adj_copy(n2)) = n4
          adj_copy(n2) = adj_copy(n2) + 1
          adj(adj_copy(n4)) = n2
          adj_copy(n4) = adj_copy(n4) + 1
        end if
c
c  Maybe add
c    32 / 23
c    52 / 25
c    53 / 35
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n2)) = n3
          adj_copy(n2) = adj_copy(n2) + 1
          adj(adj_copy(n3)) = n2
          adj_copy(n3) = adj_copy(n3) + 1
          adj(adj_copy(n2)) = n5
          adj_copy(n2) = adj_copy(n2) + 1
          adj(adj_copy(n5)) = n2
          adj_copy(n5) = adj_copy(n5) + 1
          adj(adj_copy(n3)) = n5
          adj_copy(n3) = adj_copy(n3) + 1
          adj(adj_copy(n5)) = n3
          adj_copy(n5) = adj_copy(n5) + 1
        end if
c
c  Maybe add
c    31 / 13
c    61 / 16
c    63 / 36
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj(adj_copy(n1)) = n3
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n3)) = n1
          adj_copy(n3) = adj_copy(n3) + 1
          adj(adj_copy(n1)) = n6
          adj_copy(n1) = adj_copy(n1) + 1
          adj(adj_copy(n6)) = n1
          adj_copy(n6) = adj_copy(n6) + 1
          adj(adj_copy(n3)) = n6
          adj_copy(n3) = adj_copy(n3) + 1
          adj(adj_copy(n6)) = n3
          adj_copy(n6) = adj_copy(n6) + 1
        end if
          
      end do
c
c  Ascending sort the entries for each node.
c
      do node = 1, node_num
        k1 = adj_col(node)
        k2 = adj_col(node+1)-1
        number = k2 + 1 - k1
        call i4vec_sort_heap_a ( number, adj(k1) )
      end do

      return
      end
      subroutine triangulation_order6_boundary_edge_count ( 
     &  triangle_num, triangle_node, boundary_edge_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the boundary edges.
c
c  Discussion:
c
c    This routine is given a triangulation, a set of 6-node triangles.
c    It is assumed that, in each list of 6 nodes, the vertices are listed
c    first, in counterclockwise order, followed by the three midside nodes,
c    in counterclockwise order, starting with the node between vertices
c    1 and 2.
c
c    It is assumed that each edge of the triangulation is either 
c    * an INTERIOR edge, which is listed twice, once with positive
c      orientation and once with negative orientation, or;
c    * a BOUNDARY edge, which will occur only once.
c
c    This routine should work even if the region has holes - as long
c    as the boundary of the hole comprises more than 3 edges!
c
c    Except for the dimension of TRIANGLE, this routine is identical
c    to the routine for the order 3 case.
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), the nodes that
c    make up the triangles.  These should be listed in counterclockwise order.
c
c    Output, integer BOUNDARY_EDGE_NUM, the number of boundary 
c    edges.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer boundary_edge_num
      integer e1(3*triangle_num)
      integer e2(3*triangle_num)
      integer edge(2,3*triangle_num)
      integer i
      integer interior_edge_num
      integer j
      integer k
      integer m
      integer n
      integer triangle_node(triangle_order,triangle_num)
      integer unique_num

      m = 2
      n = 3 * triangle_num
c
c  Set up the edge array.
c
      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j) = triangle_node(i,j)
        end do
      end do

      do i = 1, 2
        do j = 1, triangle_num
          edge(i,j+triangle_num) = triangle_node(i+1,j)
        end do
      end do

      do j = 1, triangle_num
        edge(1,j+2*triangle_num) = triangle_node(3,j)
        edge(2,j+2*triangle_num) = triangle_node(1,j)
      end do
c
c  In each column, force the smaller entry to appear first.
c
      do j = 1, n
        i = min ( edge(1,j), edge(2,j) )
        k = max ( edge(1,j), edge(2,j) )
        edge(1,j) = i
        edge(2,j) = k
      end do
c
c  Ascending sort the column array.
c
      call i4col_sort_a ( m, n, edge )
c
c  Get the number of unique columns in EDGE.
c
      call i4col_sorted_unique_count ( m, n, edge, unique_num )

      interior_edge_num = 3 * triangle_num - unique_num

      boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num

      return
      end
      subroutine triangulation_order6_boundary_edge_count_euler ( 
     &  node_num, triangle_num, hole_num, boundary_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
c
c  Discussion:
c
c    We assume we are given information about an order 6 triangulation
c    of a set of nodes in the plane.  
c
c    By ignoring the midside nodes, we can determine the corresponding
c    information for an order 3 triangulation, and apply Euler's formula 
c    to determine the number of edges that lie on the boundary of the 
c    set of nodes.
c
c    Thus, if we have TRIANGLE_NUM triangles, and NODE_NUM nodes, we
c    imagine that each triangle is replaced by 4 triangles, created
c    by adding the edges created by joining the midside nodes.
c
c    Thus, for 4 * TRIANGLE_NUM triangles, we can apply Euler's formula
c    to compute the number of boundary edges.  
c
c    Now, to adjust the data to our order 6 triangles, we divide the
c    number of boundary edges by 2.
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
c  Reference:
c
c    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
c    Computational Geometry, Section 9.1,
c    Springer, 2000.
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer HOLE_NUM, the number of internal nodes.
c
c    Output, integer BOUNDARY_NUM, the number of edges that lie 
c    on the boundary of the triangulation.
c
      implicit none

      integer boundary_num
      integer hole_num
      integer node_num
      integer triangle_num

      boundary_num = 
     &  ( 2 * node_num + 2 * hole_num - 4 * triangle_num - 2 ) / 2

      return
      end
      subroutine triangulation_order6_boundary_node ( node_num, 
     &  triangle_num, triangle_node, node_boundary )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_BOUNDARY_NODE indicates which nodes are on the boundary.
c
c  Discussion:
c
c    This routine is given a triangulation, an abstract list of triples
c    of nodes.  It is assumed that the nodes in each triangle are listed
c    in a counterclockwise order, although the routine should work
c    if the nodes are consistently listed in a clockwise order as well.
c
c    It is assumed that each edge of the triangulation is either
c    * an INTERIOR edge, which is listed twice, once with positive
c      orientation and once with negative orientation, or;
c    * a BOUNDARY edge, which will occur only once.
c
c    This routine should work even if the region has holes - as long
c    as the boundary of the hole comprises more than 3 edgesc
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
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), the nodes
c    that make up the triangles.
c
c    Output, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node
c    is on a boundary edge.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer e1(3*triangle_num)
      integer e2(3*triangle_num)
      integer edge(3,3*triangle_num)
      integer i
      integer j
      integer k
      integer m
      integer n
      integer node
      logical node_boundary(node_num)
      integer triangle_node(triangle_order,triangle_num)

      m = 3
      n = 3 * triangle_num
c
c  Set up the edge array.  The midside node is listed last, as
c  it is not needed for the sorting process.
c
      do j = 1, triangle_num
        edge(1,j) = triangle_node(1,j)
        edge(2,j) = triangle_node(4,j)
        edge(3,j) = triangle_node(2,j)
      end do

      do j = 1, triangle_num
        edge(1,j+triangle_num) = triangle_node(2,j)
        edge(2,j+triangle_num) = triangle_node(5,j)
        edge(3,j+triangle_num) = triangle_node(2,j)
      end do

      do j = 1, triangle_num
        edge(1,j+2*triangle_num) = triangle_node(3,j)
        edge(2,j+2*triangle_num) = triangle_node(6,j)
        edge(3,j+2*triangle_num) = triangle_node(1,j)
      end do
c
c  In each column, force the smaller of the two vertices to appear first.
c
      do j = 1, n
        i = min ( edge(1,j), edge(2,j) )
        k = max ( edge(1,j), edge(2,j) )
        edge(1,j) = i
        edge(2,j) = k
      end do
c
c  Ascending sort the column array.
c
      call i4col_sort_a ( m, n, edge )
c
c  Records which appear twice are internal edges and can be ignored.
c
      do node = 1, node_num
        node_boundary(node) = .false.
      end do

      i = 0

10    continue

      if ( i .lt. 3 * triangle_num ) then

        i = i + 1

        if ( i .eq. 3 * triangle_num ) then
          node_boundary(edge(1:m,i)) = .true.
        else if ( edge(1,i) == edge(1,i+1) .and.
     &            edge(2,i) .eq. edge(2,i+1) .and.
     &            edge(3,i) .eq. edge(3,i+1) ) then
          i = i + 1
        else
          node_boundary(edge(1,i)) = .true.
          node_boundary(edge(2,i)) = .true.
          node_boundary(edge(3,i)) = .true.
        end if

        go to 10

      end if

      return
      end
      subroutine triangulation_order6_example1 ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_EXAMPLE1 sets up a sample triangulation.
c
c  Discussion:
c
c    This triangulation is actually a Delaunay triangulation.
c
c    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
c    determined by calling TRIANGULATION_ORDER6_EXAMPLE1_SIZE first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
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
c    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of 
c    the nodes.
c
c    Output, integer TRIANGLE_ORDER(6,TRIANGLE_NUM), the nodes 
c    that make up the triangles.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbors on each side.  Negative values indicate edges that 
c    lie on the exterior.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer node_num_save
      parameter ( node_num_save = 48 )
      integer triangle_num
      integer triangle_num_save
      parameter ( triangle_num_save = 16 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      double precision node_xy(dim_num,node_num)
      double precision node_xy_save(dim_num,node_num_save)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_neighbor_save(3,triangle_num_save)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_node_save(triangle_order,triangle_num_save)

      save node_xy_save
      save triangle_neighbor_save
      save triangle_node_save

      data node_xy_save /
     &     0.0D+00, 0.0D+00, 
     &     1.0D+00, 0.0D+00, 
     &     2.0D+00, 0.0D+00, 
     &     3.0D+00, 0.0D+00, 
     &     4.0D+00, 0.0D+00, 
     &     5.0D+00, 0.0D+00, 
     &     6.0D+00, 0.0D+00, 
     &     0.0D+00, 1.0D+00, 
     &     1.0D+00, 1.0D+00, 
     &     2.0D+00, 1.0D+00, 
     &     3.0D+00, 1.0D+00, 
     &     4.0D+00, 1.0D+00, 
     &     5.0D+00, 1.0D+00, 
     &     6.0D+00, 1.0D+00, 
     &     0.0D+00, 2.0D+00, 
     &     1.0D+00, 2.0D+00, 
     &     2.0D+00, 2.0D+00, 
     &     3.0D+00, 2.0D+00, 
     &     4.0D+00, 2.0D+00, 
     &     5.0D+00, 2.0D+00, 
     &     6.0D+00, 2.0D+00, 
     &     0.0D+00, 3.0D+00, 
     &     1.0D+00, 3.0D+00, 
     &     2.0D+00, 3.0D+00, 
     &     3.0D+00, 3.0D+00, 
     &     5.0D+00, 3.0D+00, 
     &     6.0D+00, 3.0D+00, 
     &     0.0D+00, 4.0D+00, 
     &     1.0D+00, 4.0D+00, 
     &     2.0D+00, 4.0D+00, 
     &     3.0D+00, 4.0D+00, 
     &     4.0D+00, 4.0D+00, 
     &     5.0D+00, 4.0D+00, 
     &     6.0D+00, 4.0D+00, 
     &     0.0D+00, 5.0D+00, 
     &     1.0D+00, 5.0D+00, 
     &     2.0D+00, 5.0D+00, 
     &     3.0D+00, 5.0D+00, 
     &     4.0D+00, 5.0D+00, 
     &     5.0D+00, 5.0D+00, 
     &     6.0D+00, 5.0D+00, 
     &     0.0D+00, 6.0D+00, 
     &     1.0D+00, 6.0D+00, 
     &     2.0D+00, 6.0D+00, 
     &     3.0D+00, 6.0D+00, 
     &     4.0D+00, 6.0D+00, 
     &     5.0D+00, 6.0D+00, 
     &     6.0D+00, 6.0D+00 /

      data triangle_node_save /
     &   1,  3, 15,  2,  9,  8, 
     &  17, 15,  3, 16,  9, 10, 
     &   5, 19,  3, 12, 11,  4, 
     &  17,  3, 19, 10, 11, 18, 
     &   7, 21,  5, 14, 13,  6, 
     &  19,  5, 21, 12, 13, 20, 
     &  17, 30, 15, 24, 23, 16, 
     &  28, 15, 30, 22, 23, 29, 
     &  30, 17, 32, 24, 25, 31, 
     &  21, 34, 19, 27, 26, 20, 
     &  30, 44, 28, 37, 36, 29, 
     &  42, 28, 44, 35, 36, 43, 
     &  32, 46, 30, 39, 38, 31, 
     &  44, 30, 46, 37, 38, 45, 
     &  32, 34, 46, 33, 40, 39, 
     &  48, 46, 34, 47, 40, 41 /

      data triangle_neighbor_save /
     &  -3,   2,  -5, 
     &   7,   1,   4, 
     &   6,   4, -11, 
     &   2,   3, -14, 
     & -15,   6, -17, 
     &   3,   5,  10, 
     &   9,   8,   2, 
     & -24,   7,  11, 
     &   7, -28,  13, 
     &  27, -31,   6, 
     &   8,  14,  12, 
     & -36,  11, -38, 
     &  15,  14,   9, 
     &  11,  13, -44, 
     & -45,  16,  13, 
     & -48,  15, -50 /

      call r8mat_copy ( dim_num, node_num, node_xy_save, node_xy )

      call i4mat_copy ( 3, triangle_num, triangle_neighbor_save,
     &  triangle_neighbor )

      call i4mat_copy ( triangle_order, triangle_num, 
     &  triangle_node_save, triangle_node )

      return
      end
      subroutine triangulation_order6_example1_size ( node_num, 
     &  triangle_num, hole_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_EXAMPLE1_SIZE sets sizes for a sample triangulation.
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
c  Parameters:
c
c    Output, integer NODE_NUM, the number of nodes.  
c
c    Output, integer TRIANGLE_NUM, the number of triangles. 
c
c    Output, integer HOLE_NUM, the number of holes.
c
      implicit none

      integer hole_num
      integer node_num
      integer triangle_num

      node_num = 48
      triangle_num = 16
      hole_num = 1

      return
      end
      subroutine triangulation_order6_example2 ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_EXAMPLE2 returns an example triangulation.
c
c  Discussion:
c
c    This triangulation is actually a Delaunay triangulation.
c
c    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
c    determined by calling TRIANGULATION_ORDER6_EXAMPLE2_SIZE first.
c
c  Diagram:
c
c   21-22-23-24-25
c    |\  6 |\  8 |
c    | \   | \   |
c   16 17 18 19 20
c    |   \ |   \ |
c    | 5  \| 7  \|
c   11-12-13-14-15
c    |\  2 |\  4 |
c    | \   | \   |
c    6  7  8  9 10
c    | 1 \ | 3 \ |
c    |    \|    \|
c    1--2--3--4--5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of 
c    the nodes.
c
c    Output, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer node_num_save
      parameter ( node_num_save = 25 )
      integer triangle_num
      integer triangle_num_save
      parameter ( triangle_num_save = 8 )
      integer triangle_order
      parameter ( triangle_order = 6 )

      double precision node_xy(dim_num,node_num)
      double precision node_xy_save(dim_num,node_num_save)
      integer triangle_neighbor(3,triangle_num)
      integer triangle_neighbor_save(3,triangle_num_save)
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_node_save(triangle_order,triangle_num_save)

      save node_xy_save
      save triangle_neighbor_save
      save triangle_node_save

      data node_xy_save /
     &  0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 
     &  2.0D+00, 0.0D+00, 
     &  3.0D+00, 0.0D+00, 
     &  4.0D+00, 0.0D+00, 
     &  0.0D+00, 1.0D+00, 
     &  1.0D+00, 1.0D+00, 
     &  2.0D+00, 1.0D+00, 
     &  3.0D+00, 1.0D+00, 
     &  4.0D+00, 1.0D+00, 
     &  0.0D+00, 2.0D+00, 
     &  1.0D+00, 2.0D+00, 
     &  2.0D+00, 2.0D+00, 
     &  3.0D+00, 2.0D+00, 
     &  4.0D+00, 2.0D+00, 
     &  0.0D+00, 3.0D+00, 
     &  1.0D+00, 3.0D+00, 
     &  2.0D+00, 3.0D+00, 
     &  3.0D+00, 3.0D+00, 
     &  4.0D+00, 3.0D+00, 
     &  0.0D+00, 4.0D+00, 
     &  1.0D+00, 4.0D+00, 
     &  2.0D+00, 4.0D+00, 
     &  3.0D+00, 4.0D+00, 
     &  4.0D+00, 4.0D+00  /

      data triangle_node_save /
     &   1,  3, 11,  2,  7,  6, 
     &  13, 11,  3, 12,  7,  8, 
     &   3,  5, 13,  4,  9,  8, 
     &  15, 13,  5, 14,  9, 10, 
     &  11, 13, 21, 12, 17, 16, 
     &  23, 21, 13, 22, 17, 18, 
     &  13, 15, 23, 14, 19, 18, 
     &  25, 23, 15, 24, 19, 20 /

      data triangle_neighbor_save /
     &  -1,  2, -1, 
     &   5,  1,  3, 
     &  -1,  4,  2, 
     &   7,  3, -1, 
     &   2,  6, -1, 
     &  -1,  5,  7, 
     &   4,  8,  6, 
     &  -1,  7, -1 /

      call r8mat_copy ( dim_num, node_num, node_xy_save, node_xy )

      call i4mat_copy ( 3, triangle_num, triangle_neighbor_save,
     &  triangle_neighbor )

      call i4mat_copy ( triangle_order, triangle_num, 
     &  triangle_node_save, triangle_node )

      return
      end
      subroutine triangulation_order6_example2_size ( node_num, 
     &  triangle_num, hole_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_EXAMPLE2_SIZE returns the size of an example.
c
c  Diagram:
c
c   21-22-23-24-25
c    |\  6 |\  8 |
c    | \   | \   |
c   16 17 18 19 20
c    |   \ |   \ |
c    | 5  \| 7  \|
c   11-12-13-14-15
c    |\  2 |\  4 |
c    | \   | \   |
c    6  7  8  9 10
c    | 1 \ | 3 \ |
c    |    \|    \|
c    1--2--3--4--5
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
c  Parameters
c
c    Output, integer NODE_NUM, the number of nodes.
c
c    Output, integer TRIANGLE_NUM, the number of triangles.
c
c    Output, integer HOLE_NUM, the number of holes.
c
      implicit none

      integer hole_num
      integer node_num
      integer triangle_num

      node_num = 25
      triangle_num = 8
      hole_num = 0

      return
      end
      subroutine triangulation_order6_neighbor ( triangle_num, 
     &  triangle_node, t1, s1, t2, s2 )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_NEIGHBOR determines a neighbor of a given triangle.
c
c  Discussion:
c
c    A set of nodes is given.  A triangulation of the nodes has been
c    defined and recorded in TRIANGLE.  The TRIANGLE data structure records
c    triangles as sets of six nodes, with the first three being the
c    vertices, in counterclockwise order.  The fourth node is the
c    midside node between nodes 1 and 2, and the other two are listed
c    logically.
c
c    The nodes of the triangle are listed in counterclockwise order.
c    This means that if two triangles share a side, then the nodes
c    defining that side occur in the order (N1,N2) for one triangle,
c    and (N2,N1) for the other.
c
c    The routine is given a triangle and a side, and asked to find
c    another triangle (if any) that shares that side.  The routine
c    simply searches the TRIANGLE structure for an occurrence of the
c    nodes in the opposite order.
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_ORDER(6,TRIANGLE_NUM), the nodes 
c    that define each triangle.
c
c    Input, integer T1, the index of the triangle.
c
c    Input, integer S1, the index of the triangle side.
c
c    Output, integer T2, the index of the triangle which is the
c    neighbor to T1 on side S1, or -1 if there is no such neighbor.
c
c    Output, integer S2, the index of the side of triangle T2 which
c    is shared with triangle T1, or -1 if there is no such neighbor.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer i4_wrap
      integer n1
      integer n2
      integer s
      integer s1
      integer s2
      integer ss
      integer t
      integer t1
      integer t2
      integer triangle_node(triangle_order,triangle_num)

      t2 = -1
      s2 = -1

      n1 = triangle_node(s1,t1)
      ss = s1 + 1
      ss = i4_wrap ( ss, 1, 3 )
      n2 = triangle_node(ss,t1)

      do t = 1, triangle_num
        do s = 1, 3
          if ( triangle_node(s,t) .eq. n1 ) then
            ss = s - 1
            ss = i4_wrap ( ss, 1, 3 )
            if ( triangle_node(ss,t) .eq. n2 ) then
              t2 = t
              s2 = ss
              return
            end if
          end if
        end do
      end do

      return
      end
      subroutine triangulation_order6_neighbor_triangles ( 
     &  triangle_num, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_NEIGHBOR_TRIANGLES determines triangle neighbors.
c
c  Discussion:
c
c    A triangulation of a set of nodes can be completely described by
c    the coordinates of the nodes, and the list of nodes that make up
c    each triangle.  However, in some cases, it is necessary to know
c    triangle adjacency information, that is, which triangle, if any,
c    is adjacent to a given triangle on a particular side.
c
c    This routine creates a data structure recording this information.
c
c    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
c    data items.
c
c    Note that COL is a work array allocated dynamically inside this
c    routine.  It is possible, for very large values of TRIANGLE_NUM,
c    that the necessary amount of memory will not be accessible, and the
c    routine will fail.  This is a limitation of the implementation of
c    dynamic arrays in FORTRAN.  One way to get around this would be
c    to require the user to declare COL in the calling routine
c    as an allocatable array, get the necessary memory explicitly with
c    an ALLOCATE statement, and then pass COL into this routine.
c
c    Of course, the point of dynamic arrays was to make it easy to
c    hide these sorts of temporary work arrays from the poor user!
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangles that are direct neighbors of a given triangle.  
c    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
c    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative
c    if there is no neighbor on that side.  In this case, that side of the 
c    triangle lies on the boundary of the triangulation.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer col(4,3*triangle_num)
      integer i
      integer icol
      integer j
      integer k
      integer side1
      integer side2
      integer triangle_neighbor(3,triangle_num)
      integer tri
      integer triangle_node(triangle_order,triangle_num)
      integer tri1
      integer tri2
c
c  Step 1.
c  From the list of nodes for triangle T, of the form: (I,J,K)
c  construct the three neighbor relations:
c
c    (I,J,1,T) or (J,I,1,T),
c    (J,K,2,T) or (K,J,2,T),
c    (K,I,3,T) or (I,K,3,T)
c
c  where we choose (I,J,1,T) if I .lt. J, or else (J,I,1,T)
c
      do tri = 1, triangle_num

        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)

        if ( i .lt. j ) then
          col(1,3*(tri-1)+1) = i
          col(2,3*(tri-1)+1) = j
          col(3,3*(tri-1)+1) = 1
          col(4,3*(tri-1)+1) = tri
        else
          col(1,3*(tri-1)+1) = j
          col(2,3*(tri-1)+1) = i
          col(3,3*(tri-1)+1) = 1
          col(4,3*(tri-1)+1) = tri
        end if

        if ( j .lt. k ) then
          col(1,3*(tri-1)+2) = j
          col(2,3*(tri-1)+2) = k
          col(3,3*(tri-1)+2) = 2
          col(4,3*(tri-1)+2) = tri
        else
          col(1,3*(tri-1)+2) = k
          col(2,3*(tri-1)+2) = j
          col(3,3*(tri-1)+2) = 2
          col(4,3*(tri-1)+2) = tri
        end if

        if ( k .lt. i ) then
          col(1,3*(tri-1)+3) = k
          col(2,3*(tri-1)+3) = i
          col(3,3*(tri-1)+3) = 3
          col(4,3*(tri-1)+3) = tri
        else
          col(1,3*(tri-1)+3) = i
          col(2,3*(tri-1)+3) = k
          col(3,3*(tri-1)+3) = 3
          col(4,3*(tri-1)+3) = tri
        end if

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c  We only intend to sort on columns 1 and 2; the routine we call here
c  sorts on columns 1 through 4 but that won't hurt us.
c
c  What we need is to find cases where two triangles share an edge.
c  Say they share an edge defined by the nodes I and J.  Then there are
c  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
c  we make sure that these two columns occur consecutively.  That will
c  make it easy to notice that the triangles are neighbors.
c
      call i4col_sort_a ( 4, 3*triangle_num, col )
c
c  Step 3. Neighboring triangles show up as consecutive columns with
c  identical first two entries.  Whenever you spot this happening,
c  make the appropriate entries in TRIANGLE_NEIGHBOR.
c
      do i = 1, 3
        do j = 1, triangle_num
          triangle_neighbor(i,j) = -1
        end do
      end do

      icol = 1

10    continue

        if ( 3 * triangle_num .le. icol ) then
          go to 20
        end if

        if ( col(1,icol) .ne. col(1,icol+1) .or. 
     &       col(2,icol) .ne. col(2,icol+1) ) then
          icol = icol + 1
          go to 10
        end if

        side1 = col(3,icol)
        tri1 = col(4,icol)
        side2 = col(3,icol+1)
        tri2 = col(4,icol+1)

        triangle_neighbor(side1,tri1) = tri2
        triangle_neighbor(side2,tri2) = tri1

        icol = icol + 2

      go to 10

20    continue

      return
      end
      subroutine triangulation_order6_plot ( file_name, node_num, 
     &  node_xy, triangle_num, triangle_node, node_show, 
     &  triangle_show )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
c
c  Discussion:
c
c    The triangulation is most usually a Delaunay triangulation,
c    but this is not necessary.
c
c    In a six node triangulation, it is assumed that nodes 1, 2, and 3
c    are the vertices of the triangles, and that nodes 4, 5, and 6
c    lie between 1 and 2, 2 and 3, and 3 and 1 respectively.
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
c  Parameters:
c
c    Input, character * ( * ) FILE_NAME, the name of the output file.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists, for 
c    each triangle, the indices of the nodes that form the vertices of the 
c    triangle.
c
c    Input, integer NODE_SHOW,
c    0, do not show nodes;
c    1, show nodes;
c    2, show nodes and label them.
c
c    Input, integer TRIANGLE_SHOW,
c    0, do not show triangles;
c    1, show triangles;
c    2, show triangles and label them.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      double precision ave_x
      double precision ave_y
      integer circle_size
      integer delta
      integer e
      character * ( * )  file_name
      integer file_unit
      integer i
      integer i4_wrap
      integer node
      integer node_show
      double precision node_xy(2,node_num)
      character * ( 40 ) string
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_show
      double precision x_max
      double precision x_min
      integer x_ps
      integer x_ps_max
      integer x_ps_max_clip
      integer x_ps_min
      integer x_ps_min_clip
      double precision x_scale
      double precision y_max
      double precision y_min
      integer y_ps
      integer y_ps_max
      integer y_ps_max_clip
      integer y_ps_min
      integer y_ps_min_clip
      double precision y_scale
c
c  We need to do some figuring here, so that we can determine
c  the range of the data, and hence the height and width
c  of the piece of paper.
c
      x_max = node_xy(1,1)
      x_min = node_xy(1,1)
      do node = 2, node_num
        x_max = max ( x_max, node_xy(1,node) )
        x_min = min ( x_min, node_xy(1,node) )
      end do
      x_scale = x_max - x_min

      x_max = x_max + 0.05D+00 * x_scale
      x_min = x_min - 0.05D+00 * x_scale
      x_scale = x_max - x_min

      y_max = node_xy(2,1)
      y_min = node_xy(2,1)
      do node = 2, node_num
        y_max = max ( y_max, node_xy(2,node) )
        y_min = min ( y_min, node_xy(2,node) )
      end do
      y_scale = y_max - y_min

      y_max = y_max + 0.05D+00 * y_scale
      y_min = y_min - 0.05D+00 * y_scale
      y_scale = y_max - y_min

      x_ps_max = 576
      x_ps_max_clip = 594
      x_ps_min = 36
      x_ps_min_clip = 18
      y_ps_max = 666
      y_ps_max_clip = 684
      y_ps_min = 126
      y_ps_min_clip = 108
      
      if ( x_scale .lt. y_scale ) then

        delta = nint ( dble ( x_ps_max - x_ps_min ) 
     &    * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

        x_ps_max = x_ps_max - delta
        x_ps_min = x_ps_min + delta

        x_ps_max_clip = x_ps_max_clip - delta
        x_ps_min_clip = x_ps_min_clip + delta

        x_scale = y_scale

      else if ( y_scale .lt. x_scale ) then

        delta = nint ( dble ( y_ps_max - y_ps_min ) 
     &    * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

        y_ps_max      = y_ps_max - delta
        y_ps_min      = y_ps_min + delta

        y_ps_max_clip = y_ps_max_clip - delta
        y_ps_min_clip = y_ps_min_clip + delta

        y_scale = x_scale

      end if
c
c  Open the file.
c
      call get_unit ( file_unit )

      open ( unit = file_unit, file = file_name, status = 'replace' )

      write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
      write ( file_unit, '(a)' ) 
     &  '%%Creator: triangulation_order6_plot.f'
      write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
      write ( file_unit, '(a)' ) '%%Pages: 1'
      write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ',
     &  x_ps_min, y_ps_min, x_ps_max, y_ps_max
      write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
      write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
      write ( file_unit, '(a)' ) '%%EndComments'
      write ( file_unit, '(a)' ) '%%BeginProlog'
      write ( file_unit, '(a)' ) '/inch {72 mul} def'
      write ( file_unit, '(a)' ) '%%EndProlog'
      write ( file_unit, '(a)' ) '%%Page: 1 1'
      write ( file_unit, '(a)' ) 'save'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Increase line width from default 0.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '2 setlinewidth'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Set the RGB line color to very light gray.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 
     &  '%  Draw a gray border around the page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_min, ' moveto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_max, y_ps_min, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_max, y_ps_max, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_max, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) 
     &  '  ', x_ps_min, y_ps_min, ' lineto'
      write ( file_unit, '(a)' ) 'stroke'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Set the font and its size.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '/Times-Roman findfont'
      write ( file_unit, '(a)' ) '0.50 inch scalefont'
      write ( file_unit, '(a)' ) 'setfont'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Print a title.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  210  702  moveto'
      write ( file_unit, '(a)' ) '%  (Triangulation)  show'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'newpath'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' moveto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_max_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_max_clip, ' lineto'
      write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', 
     &  x_ps_min_clip, y_ps_min_clip, ' lineto'
      write ( file_unit, '(a)' ) 'clip newpath'
c
c  Draw the nodes.
c
      if ( node_num .le. 200 ) then
        circle_size = 5
      else if ( node_num .le. 500 ) then
        circle_size = 4
      else if ( node_num .le. 1000 ) then
        circle_size = 3
      else if ( node_num .le. 5000 ) then
        circle_size = 2
      else
        circle_size = 1
      end if

      if ( 1 .le. node_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (         node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 
     &      'newpath ', x_ps, y_ps, 
     &      circle_size, '0 360 arc closepath fill'

        end do

      end if
c
c  Label the nodes.
c
      if ( 2 .le. node_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the nodes:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) 
     &    '%  Set the RGB color to darker blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) * dble ( x_ps_min )   
     &      + (       + node_xy(1,node) - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max                   - y_min ) )

          write ( string, '(i4)' ) node
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, 
     &      ' moveto (' // trim ( string ) // ') show'

        end do

      end if
c
c  Draw the triangles.
c
      if ( 1 .le. triangle_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw the triangles.'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          write ( file_unit, '(a)' ) 'newpath'

          node = triangle_node(6,triangle)

          x_ps = int ( 
     &      ( ( x_max - node_xy(1,node)         ) 
     &      * dble ( x_ps_min )   
     &      + (         node_xy(1,node) - x_min ) 
     &      * dble ( x_ps_max ) )
     &      / ( x_max                   - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - node_xy(2,node)         ) 
     &      * dble ( y_ps_min )   
     &      + (         node_xy(2,node) - y_min ) 
     &      * dble ( y_ps_max ) )
     &      / ( y_max                   - y_min ) )

          write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &      x_ps, y_ps, ' moveto'

          do i = 1, 3

            node = triangle_node(i,triangle)

            x_ps = int ( 
     &        ( ( x_max - node_xy(1,node)         ) 
     &        * dble ( x_ps_min )   
     &        + (         node_xy(1,node) - x_min ) 
     &        * dble ( x_ps_max ) )
     &        / ( x_max                   - x_min ) )

            y_ps = int ( 
     &        ( ( y_max - node_xy(2,node)         ) 
     &        * dble ( y_ps_min )   
     &        + (         node_xy(2,node) - y_min ) 
     &        * dble ( y_ps_max ) )
     &        / ( y_max                   - y_min ) )

            write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &        x_ps, y_ps, ' lineto'

            node = triangle_node(i+3,triangle)

            x_ps = int ( 
     &        ( ( x_max - node_xy(1,node)         ) 
     &        * dble ( x_ps_min )   
     &        + (         node_xy(1,node) - x_min ) 
     &        * dble ( x_ps_max ) )
     &        / ( x_max                   - x_min ) )

            y_ps = int ( 
     &        ( ( y_max - node_xy(2,node)         ) 
     &        * dble ( y_ps_min )   
     &        + (         node_xy(2,node) - y_min ) 
     &        * dble ( y_ps_max ) )
     &        / ( y_max                   - y_min ) )

            write ( file_unit, '(i3,2x,i3,2x,a)' ) 
     &        x_ps, y_ps, ' lineto'

          end do

          write ( file_unit, '(a)' ) 'stroke'

        end do

      end if
c
c  Label the triangles.
c
      if ( 2 .le. triangle_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the triangles:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.20 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, triangle_num

          ave_x = 0.0D+00
          ave_y = 0.0D+00

          do i = 1, 6

            node = triangle_node(i,triangle)

            ave_x = ave_x + node_xy(1,node)
            ave_y = ave_y + node_xy(2,node)

          end do

          ave_x = ave_x / 6.0D+00
          ave_y = ave_y / 6.0D+00

          x_ps = int ( 
     &      ( ( x_max - ave_x         ) * dble ( x_ps_min )   
     &      + (       + ave_x - x_min ) * dble ( x_ps_max ) ) 
     &      / ( x_max         - x_min ) )

          y_ps = int ( 
     &      ( ( y_max - ave_y         ) * dble ( y_ps_min )   
     &      + (         ave_y - y_min ) * dble ( y_ps_max ) ) 
     &      / ( y_max         - y_min ) )

          write ( string, '(i4)' ) triangle
          string = adjustl ( string )

          write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' 
     &      // trim ( string ) // ') show'

        end do

      end if

      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) 'restore  showpage'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%  End of page.'
      write ( file_unit, '(a)' ) '%'
      write ( file_unit, '(a)' ) '%%Trailer'
      write ( file_unit, '(a)' ) '%%EOF'
      close ( unit = file_unit )

      return
      end
      subroutine triangulation_order6_print ( node_num, triangle_num, 
     &  node_xy, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_PRINT prints information about a triangulation.
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
c    11 June 2009
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
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), the nodes 
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
      parameter ( triangle_order = 6 )

      integer boundary_num
      integer i
      integer i4_wrap
      integer j
      integer k
      integer n1
      integer n2
      integer n3
      double precision node_xy(dim_num,node_num)
      integer s
      logical skip
      integer t
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      integer vertex_list(triangle_order*triangle_num)
      integer vertex_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_ORDER6_PRINT'
      write ( *, '(a)' ) 
     &  '  Information defining an order 6 triangulation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of nodes is ', node_num

      call r8mat_transpose_print ( dim_num, node_num, node_xy, 
     &  '  Node coordinates' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of triangles is ', triangle_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Sets of six nodes are used as vertices of'
      write ( *, '(a)' ) 
     &  '  the triangles.  For each triangle, the vertices'
      write ( *, '(a)' ) '  are listed in counterclockwise order,'
      write ( *, '(a)' ) '  followed by the midside nodes.'

      call i4mat_transpose_print ( triangle_order, triangle_num, 
     &  triangle_node, '  Triangle nodes:' )

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
      write ( *, '(a)' ) '     #   Tri  Side    N1    N2    N3'
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
            n3 = triangle_node(s+3,t)
            n2 = triangle_node(i4_wrap(s+1,1,3),t)

            write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) 
     &        k, t, s, n1, n2, n3
          end if

        end do

        if ( skip ) then
          exit
        end if

      end do

      return
      end
      subroutine triangulation_order6_refine_compute ( node_num1, 
     &  triangle_num1, node_xy1, triangle_node1, node_num2, 
     &  triangle_num2, edge_data, node_xy2, triangle_node2 )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
c
c  Discussion:
c
c    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we 
c    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
c    quadratic subtriangles T1, T2, T3 and T4.
c
c    The task is more complicated by the fact that we are working with
c    a mesh of triangles, so that we want to create a node only once,
c    even though it may be shared by other triangles.  (In fact, only
c    the new nodes on the edges can be shared, and then only by at most
c    one other triangle.)
c
c            3
c           / \
c          36 35
c         / T3  \
c        6--56---5
c       / \ T4  / \
c      16 46  45  25
c     / T1  \ / T2  \
c    1--14---4--24---2
c
c    This routine is given sorted information defining the edges, and uses
c    it to build the new node and triangle arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NODE_NUM1, the number of nodes.
c
c    Input, integer TRIANGLE_NUM1, the number of triangles.
c
c    Input, double precision NODE_XY1(2,NODE_NUM1), the nodes.
c
c    Input, integer TRIANGLE_NODE1(6,TRIANGLE_NUM1), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise
c    order.
c
c    Input, integer NODE_NUM2, the number of nodes in the refined
c    mesh.
c
c    Input, integer TRIANGLE_NUM2, the number of triangles in 
c    the refined mesh.
c
c    Input, integer EDGE_DATA(5,3*TRIANGLE_NUM1), edge data.
c
c    Output, double precision NODE_XY2(2,NODE_NUM2), the refined nodes.
c
c    Output, integer TRIANGLE_NODE2(6,TRIANGLE_NUM2), the nodes 
c    that make up the triangles in the refined mesh.
c
      implicit none

      integer node_num1
      integer node_num2
      integer triangle_num1
      integer triangle_num2

      integer dim
      integer edge
      integer edge_data(5,3*triangle_num1)
      integer i
      integer l1
      integer l2
      integer l3
      integer n1
      integer n1_old
      integer n2
      integer n2_old
      integer node
      integer node_num
      double precision node_xy1(2,node_num1)
      double precision node_xy2(2,node_num2)
      integer t1
      integer t2
      integer t3
      integer t4
      integer triangle_node1(6,triangle_num1)
      integer triangle_node2(6,triangle_num2)
      integer triangle1
      integer triangle2
      integer v1
      integer v2
      integer v3
      integer v4
      integer v5
      integer v6
c
c  Step 1:
c  Copy old nodes.
c
      do node = 1, node_num1
        do dim = 1, 2
          node_xy2(dim,node) = node_xy1(dim,node)
        end do
      end do
c
c  Copy indices of existing nodes into new triangle array.
c
      do triangle2 = 1, triangle_num2
        do i = 1, 6
          triangle_node2(i,triangle2) = -1
        end do
      end do

      do triangle1 = 1, triangle_num1

        t1 = ( triangle1 - 1 ) * 4 + 1
        t2 = ( triangle1 - 1 ) * 4 + 2
        t3 = ( triangle1 - 1 ) * 4 + 3
        t4 = ( triangle1 - 1 ) * 4 + 4

        triangle_node2(1,t1) = triangle_node1(1,triangle1)
        triangle_node2(2,t1) = triangle_node1(4,triangle1)
        triangle_node2(3,t1) = triangle_node1(6,triangle1)

        triangle_node2(1,t2) = triangle_node1(4,triangle1)
        triangle_node2(2,t2) = triangle_node1(2,triangle1)
        triangle_node2(3,t2) = triangle_node1(5,triangle1)

        triangle_node2(1,t3) = triangle_node1(6,triangle1)
        triangle_node2(2,t3) = triangle_node1(5,triangle1)
        triangle_node2(3,t3) = triangle_node1(3,triangle1)

        triangle_node2(1,t4) = triangle_node1(5,triangle1)
        triangle_node2(2,t4) = triangle_node1(6,triangle1)
        triangle_node2(3,t4) = triangle_node1(4,triangle1)

      end do
c
c  Step 2.
c  Examine sorted edge information.  The first time an edge is encountered,
c  generate two new nodes, then assign them (usually) to the four subtriangles 
c  of the two triangles that share that edge.
c
      node = node_num1

      n1_old = -1
      n2_old = -1

      do edge = 1, 3 * triangle_num1

        n1 = edge_data(1,edge)
        n2 = edge_data(2,edge)

        l1 = edge_data(3,edge)
        l3 = edge_data(4,edge)

        if ( l1 .eq. 1 .and. l3 .eq. 2 ) then
          l2 = 4
        else if ( l1 .eq. 1 .and. l3 .eq. 3 ) then
          l2 = 6
        else if ( l1 .eq. 2 .and. l3 .eq. 3 ) then
          l2 = 5
        end if

        triangle1 = edge_data(5,edge)
c
c  If this is the first time we've encountered this edge,
c  create the new nodes.
c
        if ( n1 .ne. n1_old .or. n2 .ne. n2_old ) then

          n1_old = n1
          n2_old = n2

          v1 = triangle_node1(l1,triangle1)
          v2 = triangle_node1(l2,triangle1)
          v3 = triangle_node1(l3,triangle1)

          node = node + 1
          v4 = node
          do dim = 1, 2
            node_xy2(dim,node) = 
     &        0.5D+00 * ( node_xy1(dim,v1) + node_xy1(dim,v2) )
          end do

          node = node + 1
          v5 = node
          do dim = 1, 2
            node_xy2(dim,node) = 
     &        0.5D+00 * ( node_xy1(dim,v2) + node_xy1(dim,v3) )
          end do

        end if
     
        t1 = ( triangle1 - 1 ) * 4 + 1
        t2 = ( triangle1 - 1 ) * 4 + 2
        t3 = ( triangle1 - 1 ) * 4 + 3

        if ( l1 .eq. 1 .and. l3 .eq. 2 ) then

          if ( triangle_node1(1,triangle1) .eq. v1 ) then
            triangle_node2(4,t1) = v4
            triangle_node2(4,t2) = v5
          else
            triangle_node2(4,t1) = v5
            triangle_node2(4,t2) = v4
          end if

        else if ( l1 .eq. 1 .and. l3 .eq. 3 ) then

          if ( triangle_node1(1,triangle1) .eq. v1 ) then
            triangle_node2(6,t1) = v4
            triangle_node2(6,t3) = v5
          else
            triangle_node2(6,t1) = v5
            triangle_node2(6,t3) = v4
          end if

        else if ( l1 .eq. 2 .and. l3 .eq. 3 ) then

          if ( triangle_node1(2,triangle1) .eq. v1 ) then
            triangle_node2(5,t3) = v4
            triangle_node2(5,t2) = v5
          else
            triangle_node2(5,t3) = v5
            triangle_node2(5,t2) = v4
          end if

        end if

      end do
c
c  Step 3.
c  Each old triangle has a single central subtriangle, for which we now
c  need to generate three new "interior" nodes.
c
      do triangle1 = 1, triangle_num1

        v4 = triangle_node1(4,triangle1)
        v5 = triangle_node1(5,triangle1)
        v6 = triangle_node1(6,triangle1)

        t1 = ( triangle1 - 1 ) * 4 + 1
        t2 = ( triangle1 - 1 ) * 4 + 2
        t3 = ( triangle1 - 1 ) * 4 + 3
        t4 = ( triangle1 - 1 ) * 4 + 4

        node = node + 1
        do dim = 1, 2
          node_xy2(dim,node) = 
     &      0.5D+00 * ( node_xy1(dim,v5) + node_xy1(dim,v6) )
        end do
        triangle_node2(4,t4) = node
        triangle_node2(4,t3) = node

        node = node + 1
        do dim = 1, 2
          node_xy2(dim,node) = 
     &      0.5D+00 * ( node_xy1(dim,v6) + node_xy1(dim,v4) )
        end do
        triangle_node2(5,t4) = node
        triangle_node2(5,t1) = node

        node = node + 1
        do dim = 1, 2
          node_xy2(dim,node) = 
     &      0.5D+00 * ( node_xy1(dim,v4) + node_xy1(dim,v5) )
        end do
        triangle_node2(6,t4) = node
        triangle_node2(6,t2) = node

      end do

      return
      end
      subroutine triangulation_order6_refine_size ( node_num1, 
     &  triangle_num1, triangle_node1, node_num2, triangle_num2, 
     &  edge_data )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
c
c  Discussion:
c
c    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we 
c    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
c    quadratic subtriangles T1, T2, T3 and T4.
c
c    The task is more complicated by the fact that we are working with
c    a mesh of triangles, so that we want to create a node only once,
c    even though it may be shared by other triangles.  (In fact, only
c    the new nodes on the edges can be shared, and then only by at most
c    one other triangle.)
c
c            3
c           / \
c          36 35
c         / T3  \
c        6--56---5
c       / \ T4  / \
c      16 46  45  25
c     / T1  \ / T2  \
c    1--14---4--24---2
c
c    This routine determines the sizes of the resulting node and
c    triangles, and constructs an edge array that can be used to 
c    properly number the new nodes.
c
c    The primary work occurs in sorting a list related to the edges.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NODE_NUM1, the number of nodes.
c
c    Input, integer TRIANGLE_NUM1, the number of triangles.
c
c    Input, integer TRIANGLE_NODE1(6,TRIANGLE_NUM1), the nodes 
c    that make up the triangles.  These should be listed in counterclockwise
c    order.
c
c    Input, integer NODE_NUM2, the number of nodes in the refined
c    mesh.
c
c    Input, integer TRIANGLE_NUM2, the number of triangles in 
c    the refined mesh.
c
c    Output, integer EDGE_DATA(5,3*TRIANGLE_NUM1), edge data 
c    needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
c
      implicit none

      integer node_num1
      integer triangle_num1

      integer a
      integer b
      integer edge
      integer edge_data(5,3*triangle_num1)
      integer i
      integer j
      integer k
      integer n1
      integer n1_old
      integer n2
      integer n2_old
      integer node_num2
      integer triangle_num2
      integer triangle_node1(6,triangle_num1)
      integer triangle1
c
c  Step 1:
c  From the list of vertices for triangle T, of the form: (I,J,K),
c  construct the edge relations:
c
c    (I,J,1,2,T)
c    (I,K,1,3,T)
c    (J,K,2,3,T)
c
c  To make matching easier, we reorder each pair of nodes into
c  ascending order.
c
      do triangle1 = 1, triangle_num1

        i = triangle_node1(1,triangle1)
        j = triangle_node1(2,triangle1)
        k = triangle_node1(3,triangle1)

        call i4i4_sort_a ( i, j, a, b )
        edge_data(1,3*(triangle1-1)+1) = a
        edge_data(2,3*(triangle1-1)+1) = b
        edge_data(3,3*(triangle1-1)+1) = 1
        edge_data(4,3*(triangle1-1)+1) = 2
        edge_data(5,3*(triangle1-1)+1) = triangle1

        call i4i4_sort_a ( i, k, a, b )
        edge_data(1,3*(triangle1-1)+2) = a
        edge_data(2,3*(triangle1-1)+2) = b
        edge_data(3,3*(triangle1-1)+2) = 1
        edge_data(4,3*(triangle1-1)+2) = 3
        edge_data(5,3*(triangle1-1)+2) = triangle1

        call i4i4_sort_a ( j, k, a, b )
        edge_data(1,3*(triangle1-1)+3) = a
        edge_data(2,3*(triangle1-1)+3) = b
        edge_data(3,3*(triangle1-1)+3) = 2
        edge_data(4,3*(triangle1-1)+3) = 3
        edge_data(5,3*(triangle1-1)+3) = triangle1

      end do
c
c  Step 2: Perform an ascending dictionary sort on the relations.
c
      call i4col_sort_a ( 5, 3*triangle_num1, edge_data )
c
c  Step 3: Each shared edge will show up twice, consecutively,
c  in the EDGE_DATA array.  Each unique edge will generate
c  two new nodes, and each triangle will generate three new nodes.
c  
      node_num2 = node_num1

      n1_old = -1
      n2_old = -1

      do edge = 1, 3 * triangle_num1
        n1 = edge_data(1,edge)
        n2 = edge_data(2,edge)
        if ( n1 .ne. n1_old .or. n2 .ne. n2_old ) then
          node_num2 = node_num2 + 2
          n1_old = n1
          n2_old = n2
        end if
      end do

      node_num2 = node_num2 + 3 * triangle_num1

      triangle_num2 = 4 * triangle_num1

      return
      end
      subroutine triangulation_order6_to_order3 ( triangle_num1, 
     &  triangle_node1, triangle_num2, triangle_node2 )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_TO_ORDER3 linearizes a quadratic triangulation.
c
c  Discussion:
c
c    A quadratic triangulation is assumed to consist of 6-node triangles,
c    as in the following:
c
c    11-12-13-14-15
c     |\    |\    |
c     | \   | \   |
c     6  7  8  9 10
c     |   \ |   \ |
c     |    \|    \|
c     1--2--3--4--5
c
c   This routine rearranges information so as to define the 3-node
c   triangulation:
c
c    11-12-13-14-15
c     |\ |\ |\ |\ |
c     | \| \| \| \|
c     6--7--8--9-10
c     |\ |\ |\ |\ |
c     | \| \| \| \|
c     1--2--3--4--5
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
c  Parameters:
c
c    Input, integer TRIANGLE_NUM1, the number of triangles in 
c    the quadratic triangulation.
c
c    Input, integer TRIANGLE_NODE1(6,TRIANGLE_NUM1), the quadratic
c    triangulation.
c
c    Input, integer TRIANGLE_NUM2, the number of triangles in the 
c    linear triangulation.  TRIANGLE_NUM2 = 4 * TRIANGLE_NUM1.
c
c    Output, integer TRIANGLE_NODE2(3,TRIANGLE_NUM2), the linear
c    triangulation.
c
      implicit none

      integer triangle_num1
      integer triangle_num2
      integer triangle_order1
      parameter ( triangle_order1 = 6 )
      integer triangle_order2
      parameter ( triangle_order2 = 3 )

      integer n1
      integer n2
      integer n3
      integer n4
      integer n5
      integer n6
      integer tri1
      integer tri2
      integer triangle_node1(triangle_order1,triangle_num1)
      integer triangle_node2(triangle_order2,triangle_num2)

      tri2 = 0

      do tri1 = 1, triangle_num1

        n1 = triangle_node1(1,tri1)
        n2 = triangle_node1(2,tri1)
        n3 = triangle_node1(3,tri1)
        n4 = triangle_node1(4,tri1)
        n5 = triangle_node1(5,tri1)
        n6 = triangle_node1(6,tri1)

        tri2 = tri2 + 1
        triangle_node2(1,tri2) = n1
        triangle_node2(2,tri2) = n4
        triangle_node2(3,tri2) = n6
        tri2 = tri2 + 1
        triangle_node2(1,tri2) = n2
        triangle_node2(2,tri2) = n5
        triangle_node2(3,tri2) = n4
        tri2 = tri2 + 1
        triangle_node2(1,tri2) = n3
        triangle_node2(2,tri2) = n6
        triangle_node2(3,tri2) = n5
        tri2 = tri2 + 1
        triangle_node2(1,tri2) = n4
        triangle_node2(2,tri2) = n5
        triangle_node2(3,tri2) = n6

      end do

      return
      end
      subroutine triangulation_order6_vertex_count ( node_num, 
     &  triangle_num, triangle_node, vertex_num, midside_num )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_VERTEX_COUNT counts the vertex nodes in a triangulation.
c
c  Discussion:
c
c    In a triangulation of order 6, some nodes are midside nodes and some
c    nodes are vertex nodes.  
c
c    Especially when an order 6 triangulation is used to handle the
c    Navier Stokes equations, it is useful to know the number of
c    vertex and midside nodes.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  6   5  side 2
c       |    \
c    3  |     \
c       |      \
c       1---4---2
c
c         side 1
c
c    The local node numbering.  Local nodes 1, 2 and 3 are "vertex nodes",
c    while nodes 4, 5 and 6 are "midside nodes".
c
c
c   21-22-23-24-25
c    |\    |\    |
c    | \   | \   |
c   16 17 18 19 20
c    |   \ |   \ |
c    |    \|    \|
c   11-12-13-14-15
c    |\    |\    |
c    | \   | \   |
c    6  7  8  9 10
c    |   \ |   \ |
c    |    \|    \|
c    1--2--3--4--5
c
c    A sample grid, which contains 9 vertex nodes and 16 midside nodes.
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
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(6,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle.  The first three nodes are the vertices,
c    in counterclockwise order.  The fourth value is the midside
c    node between nodes 1 and 2; the fifth and sixth values are
c    the other midside nodes in the logical order.
c
c    Output, integer VERTEX_NUM, the number of nodes which are 
c    vertices.
c
c    Output, integer MIDSIDE_NUM, the number of nodes which are 
c    midsides.  This value is inferred from NODE_NUM - VERTEX_NUM, so the value
c    of NODE_NUM needs to be correct on inputc
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer midside_num
      integer node_num
      integer triangle
      integer triangle_node(triangle_order,triangle_num)
      integer vertex_num
      integer vertices(3*triangle_num)

      do triangle = 1, triangle_num
        vertices(triangle               ) = triangle_node(1,triangle)
        vertices(triangle+  triangle_num) = triangle_node(2,triangle)
        vertices(triangle+2*triangle_num) = triangle_node(3,triangle)
      end do

      call i4vec_sort_heap_a ( 3*triangle_num, vertices )

      call i4vec_sorted_unique ( 3*triangle_num, vertices, vertex_num )

      midside_num = node_num - vertex_num

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
      subroutine triangulation_search_naive ( node_num, node_xy, 
     &  element_order, element_num, element_node, p, triangle_index )

c*********************************************************************72
c
cc TRIANGULATION_SEARCH_NAIVE naively searches a triangulation.
c
c  Discussion:
c
c    The algorithm simply checks each triangle to see if point P is
c    contained in it.  Surprisingly, this is not the fastest way to
c    do the check, at least if the triangulation is Delaunay.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer ELEMENT_ORDER, the order of the triangles.
c
c    Input, integer ELEMENT_NUM, the number of triangles in
c    the triangulation.
c
c    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
c    the nodes that make up each triangle.
c
c    Input, double precision P(2), the coordinates of a point.
c
c    Output, integer TRIANGLE_INDEX, the index of the triangle
c    where the search ended, or -1 if no triangle was found containing
c    the point.
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
      double precision dxa
      double precision dxb
      double precision dxp
      double precision dya
      double precision dyb
      double precision dyp
      double precision gamma
      double precision node_xy(dim_num,node_num)
      double precision p(dim_num)
      integer triangle
      integer element_node(element_order,element_num)
      integer triangle_index

      triangle_index = -1

      do triangle = 1, element_num
c
c  Get the nodes of the triangle.
c
        a = element_node(1,triangle)
        b = element_node(2,triangle)
        c = element_node(3,triangle)
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
        beta  = ( dxa * dyp - dya * dxp ) / det
        gamma = 1.0D+00 - alpha - beta
c
c  If the barycentric coordinates are all positive, then the point
c  is inside the triangle and we're done.
c
        if ( 0.0D+00 .le. alpha .and. 
     &       0.0D+00 .le. beta  .and. 
     &       0.0D+00 .le. gamma ) then
          triangle_index = triangle
          return
        end if

      end do

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
      subroutine voronoi_polygon_area ( node, neighbor_num, 
     &  neighbor_index, node_num, node_xy, area )

c*********************************************************************72
c
cc VORONOI_POLYGON_AREA computes the area of a Voronoi polygon.
c
c  Discussion:
c
c    It is assumed that the Voronoi polygon is finite!  Every Voronoi
c    diagram includes some regions which are infinite, and for those,
c    this formula is not appropriate.
c
c    The routine is given the indices of the nodes that are neighbors of a
c    given "center" node.  A node is a neighbor of the center node if the
c    Voronoi polygons of the two nodes share an edge.  The triangles of the
c    Delaunay triangulation are formed from successive pairs of these neighbor
c    nodes along with the center node.
c
c    The assumption that the polygon is a Voronoi polygon is
c    used to determine the location of the boundaries of the polygon,
c    which are the perpendicular bisectors of the lines connecting
c    the center point to each of its neighbors.
c
c    The finiteness assumption is employed in part in the
c    assumption that the polygon is bounded by the finite
c    line segments from point 1 to 2, 2 to 3, ...,
c    M-1 to M, and M to 1, where M is the number of neighbors.
c
c    It is assumed that this subroutine is being called by a
c    process which has computed the Voronoi diagram of a large
c    set of nodes, so the arrays X and Y are dimensioned by
c    NODE_NUM, which may be much greater than the number of neighbor
c    nodes.
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
c  Reference:
c
c    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
c    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
c    Second Edition,
c    Wiley, 2000, page 485.
c
c  Parameters:
c
c    Input, integer NODE, the index of the node whose Voronoi
c    polygon is to be measured.  1 <= NODE <= NODE_NUM.
c
c    Input, integer NEIGHBOR_NUM, the number of neighbor nodes of
c    the given node.
c
c    Input, integer NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
c    of the neighbor nodes (used to access X and Y).  The neighbor
c    nodes should be listed in the (counter-clockwise) order in
c    which they occur as one circles the center node.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Output, double precision AREA, the area of the Voronoi polygon.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      integer node_num

      double precision a
      double precision area
      double precision b
      double precision c
      integer dim
      integer i
      integer j
      integer neighbor_index(neighbor_num)
      integer node
      double precision node_xy(dim_num,node_num)
      double precision pc(dim_num)
      double precision pi(dim_num)
      double precision pj(dim_num)
      double precision ui(dim_num)
      double precision uj(dim_num)

      area = 0.0D+00

      if ( node .lt. 1 .or. node_num .lt. node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VORONOI_POLYGON_AREA - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
        stop
      end if

      do dim = 1, dim_num
        pc(dim) = node_xy(dim,node)
      end do

      i = neighbor_num

      do dim = 1, dim_num
        pi(dim) = node_xy(dim,i)
      end do

      j = 1
      j = neighbor_index(j)

      do dim = 1, dim_num
        pj(dim) = node_xy(dim,j)
      end do

      a = ( pi(1)**2 + pi(2)**2 - pc(1)**2 - pc(2)**2 )
      b = ( pj(1)**2 + pj(2)**2 - pc(1)**2 - pc(2)**2 )
      c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) 
     &              - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
      uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
      uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

      do i = 1, neighbor_num

        do dim = 1, dim_num
          pi(dim) = pj(dim)
        end do

        do dim = 1, dim_num
          ui(dim) = uj(dim)
        end do

        j = i + 1
        if ( neighbor_num .lt. j ) then
          j = 1
        end if

        j = neighbor_index(j)

        do dim = 1, dim_num
          pj(dim) = node_xy(dim,j)
        end do

        a = ( pi(1)**2 + pi(2)**2 - pc(1)**2 - pc(2)**2 )
        b = ( pj(1)**2 + pj(2)**2 - pc(1)**2 - pc(2)**2 )
        c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) 
     &                - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
        uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
        uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

        area = area + uj(1) * ui(2) - ui(1) * uj(2)

      end do

      area = 0.5D+00 * area

      return
      end
      subroutine voronoi_polygon_centroid ( node, neighbor_num, 
     &  neighbor_index, node_num, node_xy, centroid )

c*********************************************************************72
c
cc VORONOI_POLYGON_CENTROID computes the centroid of a Voronoi polygon.
c
c  Discussion:
c
c    It is assumed that the Voronoi polygon is finite!  Every Voronoi
c    diagram includes some regions which are infinite, and for those,
c    this formula is not appropriate.
c
c    The routine is given the indices of the nodes that are neighbors of a
c    given "center" node.  A node is a neighbor of the center node if the
c    Voronoi polygons of the two nodes share an edge.  The triangles of the
c    Delaunay triangulation are formed from successive pairs of these neighbor
c    nodes along with the center node.
c
c    The assumption that the polygon is a Voronoi polygon is
c    used to determine the location of the boundaries of the polygon,
c    which are the perpendicular bisectors of the lines connecting
c    the center point to each of its neighbors.
c
c    The finiteness assumption is employed in part in the
c    assumption that the polygon is bounded by the finite
c    line segments from point 1 to 2, 2 to 3, ...,
c    M-1 to M, and M to 1, where M is the number of neighbors.
c
c    It is assumed that this subroutine is being called by a
c    process which has computed the Voronoi diagram of a large
c    set of nodes, so the arrays X and Y are dimensioned by
c    NODE_NUM, which may be much greater than the number of neighbor
c    nodes.
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
c  Reference:
c
c    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
c    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
c    Second Edition,
c    Wiley, 2000, page 490.
c
c  Parameters:
c
c    Input, integer NODE, the index of the node whose Voronoi
c    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
c
c    Input, integer NEIGHBOR_NUM, the number of neighbor nodes of
c    the given node.
c
c    Input, integer NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
c    of the neighbor nodes.  These indices are used to access the
c    X and Y arrays.  The neighbor nodes should be listed in the
c    (counter-clockwise) order in which they occur as one circles
c    the center node.
c
c    Input, integer NODE_NUM, the total number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Output, double precision CENTROID(2), the coordinates of the centroid
c    of the Voronoi polygon of node NODE.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      integer node_num

      double precision a
      double precision area
      double precision b
      double precision c
      double precision centroid(dim_num)
      integer dim
      integer i
      integer j
      integer neighbor_index(neighbor_num)
      integer node
      double precision node_xy(dim_num,node_num)
      double precision pc(dim_num)
      double precision pi(dim_num)
      double precision pj(dim_num)
      double precision ui(dim_num)
      double precision uj(dim_num)

      do dim = 1, dim_num
        centroid(dim) = 0.0D+00
      end do

      if ( node .lt. 1 .or. node_num .lt. node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VORONOI_POLYGON_CENTROID - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
        stop
      end if

      do dim = 1, dim_num
        pc(dim) = node_xy(dim,node)
      end do

      i = neighbor_num
      i = neighbor_index(i)

      do dim = 1, dim_num
        pi(dim) = node_xy(dim,i)
      end do

      j = 1
      j = neighbor_index(j)

      do dim = 1, dim_num
        pj(dim) = node_xy(dim,j)
      end do

      a = ( pi(1) * pi(1) + pi(2) * pi(2) 
     &    - pc(1) * pc(1) - pc(2) * pc(2) )

      b = ( pj(1) * pj(1) + pj(2) * pj(2) 
     &    - pc(1) * pc(1) - pc(2) * pc(2) )

      c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) 
     &              - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )
      uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
      uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

      do i = 1, neighbor_num

        do dim = 1, dim_num
          pi(dim) = pj(dim)
        end do

        do dim = 1, dim_num
          ui(dim) = uj(dim)
        end do

        j = i + 1
        if ( neighbor_num .lt. j ) then
          j = 1
        end if

        do dim = 1, dim_num
          pj(dim) = node_xy(dim,j)
        end do

        a = ( pi(1) * pi(1) + pi(2) * pi(2) 
     &      - pc(1) * pc(1) - pc(2) * pc(2) )

        b = ( pj(1) * pj(1) + pj(2) * pj(2) 
     &      - pc(1) * pc(1) - pc(2) * pc(2) )

        c = 2.0D+00 * ( ( pi(1) - pc(1) ) * ( pj(2) - pc(2) ) 
     &                - ( pj(1) - pc(1) ) * ( pi(2) - pc(2) ) )

        uj(1) = ( a * ( pj(2) - pc(2) ) - b * ( pi(2) - pc(2) )  ) / c
        uj(2) = ( a * ( pj(1) - pc(1) ) - b * ( pi(1) - pc(1) )  ) / c

        centroid(1) = centroid(1) + ( ui(2) - uj(2) ) 
     &    * ( ( uj(1) + ui(1) )**2 - uj(1) * ui(1) )

        centroid(2) = centroid(2) + ( ui(1) - uj(1) ) 
     &    * ( ( uj(2) + ui(2) )**2 - uj(2) * ui(2) )

      end do

      call voronoi_polygon_area ( node, neighbor_num, neighbor_index, 
     &  node_num, node_xy, area )

      do dim = 1, dim_num
        centroid(dim) = centroid(dim) / ( 6.0D+00 * area )
      end do

      return
      end
      subroutine voronoi_polygon_vertices ( node, neighbor_num, 
     &  neighbor_index, node_num, node_xy, v )

c*********************************************************************72
c
cc VORONOI_POLYGON_VERTICES computes the vertices of a Voronoi polygon.
c
c  Discussion:
c
c    This routine is only appropriate for Voronoi polygons that are finite.
c
c    The routine is given the indices of the nodes that are neighbors of a
c    given "center" node.  A node is a neighbor of the center node if the
c    Voronoi polygons of the two nodes share an edge.  The triangles of the
c    Delaunay triangulation are formed from successive pairs of these neighbor
c    nodes along with the center node.
c
c    Given only the neighbor node information, it is possible to determine
c    the location of the vertices of the polygonal Voronoi region by computing
c    the circumcenters of the Delaunay triangles.
c
c    It is assumed that this subroutine is being called by a process which has
c    computed the Voronoi diagram of a large set of nodes, so the arrays X and
c    Y are dimensioned by NODE_NUM, which may be much greater than the number
c    of neighbor nodes.
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
c  Reference:
c
c    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
c    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
c    Second Edition,
c    Wiley, 2000.
c
c  Parameters:
c
c    Input, integer NODE, the index of the node whose Voronoi
c    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
c
c    Input, integer NEIGHBOR_NUM, the number of neighbor nodes of
c    the given node.
c
c    Input, integer NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
c    of the neighbor nodes.  These indices are used to access the
c    X and Y arrays.  The neighbor nodes should be listed in the
c    (counter-clockwise) order in which they occur as one circles
c    the center node.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Output, double precision V(2,NEIGHBOR_NUM), the coordinates of
c    the vertices of the Voronoi polygon around node NODE.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer neighbor_num
      integer node_num

      integer dim
      integer i
      integer ip1
      integer neighbor_index(neighbor_num)
      integer node
      double precision node_xy(dim_num,node_num)
      double precision t(dim_num,3)
      double precision v(dim_num,neighbor_num)

      if ( node .lt. 1 .or. node_num .lt. node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VORONOI_POLYGON_VERTICES - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of input parameter NODE.'
        stop
      end if

      do dim = 1, dim_num
        t(dim,1) = node_xy(dim,node)
      end do

      ip1 = neighbor_index(1)

      do dim = 1, dim_num
        t(dim,3) = node_xy(dim,ip1)
      end do

      do i = 1, neighbor_num

        do dim = 1, dim_num
          t(dim,2) = t(dim,3)
        end do

        ip1 = i + 1
        if ( neighbor_num .lt. ip1 ) then
          ip1 = 1
        end if

        ip1 = neighbor_index(ip1)

        do dim = 1, dim_num
          t(dim,3) = node_xy(dim,ip1)
        end do

        call triangle_circumcenter_2d ( t, v(1,i) )

      end do

      return
      end
