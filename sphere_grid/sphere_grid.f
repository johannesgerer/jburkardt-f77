      function arc_cosine ( c )

c*********************************************************************72
c
cc ARC_COSINE computes the arc cosine function, with argument truncation.
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
c    02 December 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision C, the argument.
c
c    Output, double precision ARC_COSINE, an angle whose cosine is C.
c
      implicit none

      double precision arc_cosine
      double precision c
      double precision c2

      c2 = c
      c2 = max ( c2, -1.0D+00 )
      c2 = min ( c2, +1.0D+00 )

      arc_cosine = acos ( c2 )

      return
      end
      function arc_sine ( s )

c*********************************************************************72
c
cc ARC_SINE computes the arc sine function, with argument truncation.
c
c  Discussion:
c
c    If you call your system ASIN routine with an input argument that is
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
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision S, the argument.
c
c    Output, double precision ARC_SINE, an angle whose sine is S.
c
      implicit none

      double precision arc_sine
      double precision s
      double precision s2

      s2 = s
      s2 = max ( s2, -1.0D+00 )
      s2 = min ( s2, +1.0D+00 )

      arc_sine = asin ( s2 )

      return
      end
      function atan4 ( y, x )

c*********************************************************************72
c
cc ATAN4 computes the inverse tangent of the ratio Y / X.
c
c  Discussion:
c
c    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
c    the built in functions ATAN and ATAN2 already do.
c
c    However:
c
c    * ATAN4 always returns a positive angle, between 0 and 2 PI,
c      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
c      and [-PI,+PI] respectively;
c
c    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
c     function by contrast always returns an angle in the first or fourth
c     quadrants.
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
c    Input, double precision Y, X, two quantities which represent the 
c    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
c
c    Output, double precision ATAN4, an angle between 0 and 2 * PI, whose
c    tangent is (Y/X), and which lies in the appropriate quadrant so that 
c    the signs of its cosine and sine match those of X and Y.
c
      implicit none

      double precision abs_x
      double precision abs_y
      double precision atan4
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision theta_0
      double precision x
      double precision y
c
c  Special cases:
c
      if ( x .eq. 0.0D+00 ) then

        if ( 0.0D+00 .lt. y ) then
          theta = pi / 2.0D+00
        else if ( y .lt. 0.0D+00 ) then
          theta = 3.0D+00 * pi / 2.0D+00
        else if ( y .eq. 0.0D+00 ) then
          theta = 0.0D+00
        end if

      else if ( y .eq. 0.0D+00 ) then

        if ( 0.0D+00 .lt. x ) then
          theta = 0.0D+00
        else if ( x .lt. 0.0D+00 ) then
          theta = pi
        end if
c
c  We assume that ATAN2 is correct when both arguments are positive.
c
      else

        abs_y = abs ( y )
        abs_x = abs ( x )

        theta_0 = atan2 ( abs_y, abs_x )

        if ( 0.0D+00 .lt. x .and. 0.0D+00 .lt. y ) then
          theta = theta_0
        else if ( x .lt. 0.0D+00 .and. 0.0D+00 .lt. y ) then
          theta = pi - theta_0
        else if ( x .lt. 0.0D+00 .and. y .lt. 0.0D+00 ) then
          theta = pi + theta_0
        else if ( 0.0D+00 .lt. x .and. y .lt. 0.0D+00 ) then
          theta = 2.0D+00 * pi - theta_0
        end if

      end if

      atan4 = theta

      return
      end
      subroutine i4mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
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
c    Input, character * ( * ) TITLE, a title.
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
c    An I4MAT is an array of I4's.
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
c    Input, character * ( * ) TITLE, a title.
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )  title

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

          write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      subroutine icos_shape ( point_num, edge_num, face_num, 
     &  face_order_max, point_coord, edge_point, face_order, 
     &  face_point )

c*********************************************************************72
c
cc ICOS_SHAPE describes an icosahedron.
c
c  Discussion:
c
c    The input data required for this routine can be retrieved from ICOS_SIZE.
c
c    The vertices lie on the unit sphere.
c
c    The dual of an icosahedron is a dodecahedron.
c
c    The data has been rearranged from a previous assignment.  
c    The STRIPACK program refuses to triangulate data if the first
c    three nodes are "collinear" on the sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points (12).
c
c    Input, integer EDGE_NUM, the number of edges (30).
c
c    Input, integer FACE_NUM, the number of faces (20).
c
c    Input, integer FACE_ORDER_MAX, the maximum number of 
c    vertices per face (3).
c
c    Output, double precision POINT_COORD(3,POINT_NUM), the points.
c
c    Output, integer EDGE_POINT(2,EDGE_NUM), the points that 
c    make up each edge, listed in ascending order of their indexes.
c
c    Output, integer FACE_ORDER(FACE_NUM), the number of vertices
c    per face.
c
c    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
c    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
c    points are listed in the counter clockwise direction defined
c    by the outward normal at the face.  The nodes of each face are ordered 
c    so that the lowest index occurs first.  The faces are then sorted by
c    nodes.
c
      implicit none

      integer edge_num
      integer edge_order
      parameter ( edge_order = 2 )
      integer face_num
      integer face_order_max
      integer point_num

      double precision a
      double precision b
      integer edge_point(edge_order,edge_num)
      integer i
      integer face_order(face_num)
      integer face_point(face_order_max,face_num)
      double precision phi
      double precision point_coord(3,point_num)
      double precision z
c
c  Set the point coordinates.
c
      phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

      a = phi / sqrt ( 1.0D+00 + phi * phi )
      b = 1.0D+00 / sqrt ( 1.0D+00 + phi * phi )
      z = 0.0D+00
c
c  A*A + B*B + Z*Z = 1.
c
      point_coord(1,1) =   a
      point_coord(2,1) =   b
      point_coord(3,1) =   z
      point_coord(1,2) =   a
      point_coord(2,2) = - b
      point_coord(3,2) =   z
      point_coord(1,3) =   b
      point_coord(2,3) =   z
      point_coord(3,3) =   a
      point_coord(1,4) =   b
      point_coord(2,4) =   z
      point_coord(3,4) = - a
      point_coord(1,5) =   z
      point_coord(2,5) =   a
      point_coord(3,5) =   b
      point_coord(1,6) =   z
      point_coord(2,6) =   a
      point_coord(3,6) = - b
      point_coord(1,7) =   z
      point_coord(2,7) = - a
      point_coord(3,7) =   b
      point_coord(1,8) =   z
      point_coord(2,8) = - a
      point_coord(3,8) = - b
      point_coord(1,9) = - b
      point_coord(2,9) =   z
      point_coord(3,9) =   a
      point_coord(1,10) = - b
      point_coord(2,10) =   z
      point_coord(3,10) = - a
      point_coord(1,11) = - a
      point_coord(2,11) =   b
      point_coord(3,11) =   z
      point_coord(1,12) = - a
      point_coord(2,12) = - b
      point_coord(3,12) =   z
c
c  Set the edges.
c
      edge_point(1,1) = 1
      edge_point(2,1) = 2
      edge_point(1,2) = 1
      edge_point(2,2) = 3
      edge_point(1,3) = 1
      edge_point(2,3) = 4
      edge_point(1,4) = 1
      edge_point(2,4) = 5
      edge_point(1,5) = 1
      edge_point(2,5) = 6
      edge_point(1,6) = 2
      edge_point(2,6) = 3
      edge_point(1,7) = 2
      edge_point(2,7) = 4
      edge_point(1,8) = 2
      edge_point(2,8) = 7
      edge_point(1,9) = 2
      edge_point(2,9) = 8
      edge_point(1,10) = 3
      edge_point(2,10) = 5
      edge_point(1,11) = 3
      edge_point(2,11) = 7
      edge_point(1,12) = 3
      edge_point(2,12) = 9
      edge_point(1,13) = 4
      edge_point(2,13) = 6
      edge_point(1,14) = 4
      edge_point(2,14) = 8
      edge_point(1,15) = 4
      edge_point(2,15) = 10
      edge_point(1,16) = 5
      edge_point(2,16) = 6
      edge_point(1,17) = 5
      edge_point(2,17) = 9
      edge_point(1,18) = 5
      edge_point(2,18) = 11
      edge_point(1,19) = 6
      edge_point(2,19) = 10
      edge_point(1,20) = 6
      edge_point(2,20) = 11 
      edge_point(1,21) = 7
      edge_point(2,21) = 8
      edge_point(1,22) = 7
      edge_point(2,22) = 9
      edge_point(1,23) = 7
      edge_point(2,23) = 12
      edge_point(1,24) = 8
      edge_point(2,24) = 10
      edge_point(1,25) = 8
      edge_point(2,25) = 12
      edge_point(1,26) = 9
      edge_point(2,26) = 11
      edge_point(1,27) = 9
      edge_point(2,27) = 12
      edge_point(1,28) = 10
      edge_point(2,28) = 11
      edge_point(1,29) = 10
      edge_point(2,29) = 12
      edge_point(1,30) = 11
      edge_point(2,30) = 12
c
c  Set the face orders.
c
      do i = 1, face_num
        face_order(i) = 3
      end do
c
c  Set the faces.
c
      face_point(1,1) = 1
      face_point(2,1) = 2
      face_point(3,1) = 4
      face_point(1,2) = 1
      face_point(2,2) = 3
      face_point(3,2) = 2
      face_point(1,3) = 1
      face_point(2,3) = 4
      face_point(3,3) = 6
      face_point(1,4) = 1
      face_point(2,4) = 5
      face_point(3,4) = 3
      face_point(1,5) = 1
      face_point(2,5) = 6
      face_point(3,5) = 5
      face_point(1,6) = 2
      face_point(2,6) = 3
      face_point(3,6) = 7
      face_point(1,7) = 2
      face_point(2,7) = 7
      face_point(3,7) = 8
      face_point(1,8) = 2
      face_point(2,8) = 8
      face_point(3,8) = 4
      face_point(1,9) = 3
      face_point(2,9) = 5
      face_point(3,9) = 9
      face_point(1,10) = 3
      face_point(2,10) = 9
      face_point(3,10) = 7
      face_point(1,11) = 4
      face_point(2,11) = 8
      face_point(3,11) = 10
      face_point(1,12) = 4
      face_point(2,12) = 10
      face_point(3,12) = 6
      face_point(1,13) = 5
      face_point(2,13) = 6
      face_point(3,13) = 11
      face_point(1,14) = 5
      face_point(2,14) = 11
      face_point(3,14) = 9
      face_point(1,15) = 6
      face_point(2,15) = 10
      face_point(3,15) = 11
      face_point(1,16) = 7
      face_point(2,16) = 9
      face_point(3,16) = 12
      face_point(1,17) = 7
      face_point(2,17) = 12
      face_point(3,17) = 8
      face_point(1,18) = 8
      face_point(2,18) = 12
      face_point(3,18) = 10
      face_point(1,19) =  9
      face_point(2,19) = 11
      face_point(3,19) = 12
      face_point(1,20) = 10
      face_point(2,20) = 12
      face_point(3,20) = 11

      return
      end
      subroutine icos_size ( point_num, edge_num, face_num, 
     &  face_order_max )

c*********************************************************************72
c
cc ICOS_SIZE gives "sizes" for an icosahedron in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer POINT_NUM, the number of points.
c
c    Output, integer EDGE_NUM, the number of edges.
c
c    Output, integer FACE_NUM, the number of faces.
c
c    Output, integer FACE_ORDER_MAX, the maximum order of any face.
c
      implicit none

      integer edge_num
      integer face_num
      integer face_order_max
      integer point_num

      point_num = 12
      edge_num = 30
      face_num = 20
      face_order_max = 3

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
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
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
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
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
c    Input, character*(*) TITLE, a title.
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
c  Discussion:
c
c    An R8MAT is an array of R8's.
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
c    Input, character * ( * ) TITLE, a title.
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
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

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

          write ( *, '(2x,i8,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      function r8vec_diff_norm ( n, a, b )

c*********************************************************************72
c
cc R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), B(N), the vectors.
c
c    Output, double precision R8VEC_DIFF_NORM, the L2 norm of A - B.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision r8vec_diff_norm
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + ( a(i) + b(i) ) * ( a(i) - b(i) )
      end do
      value = sqrt ( value )

      r8vec_diff_norm = value

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
      function r8vec_norm ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM returns the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
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
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose L2 norm is desired.
c
c    Output, double precision R8VEC_NORM, the L2 norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r8vec_norm = value

      return
      end
      subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

c*********************************************************************72
c
cc R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The (nonzero) vector P defines a direction.
c
c    The vector A can be written as the sum
c
c      A = A_normal + A_parallel
c
c    where A_parallel is a linear multiple of P, and A_normal
c    is perpendicular to P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the vector to be polarized.
c
c    Input, double precision P(N), the polarizing direction.
c
c    Output, double precision A_NORMAL(N), A_PARALLEL(N), the normal
c    and parallel components of A.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_dot_p
      double precision a_normal(n)
      double precision a_parallel(n)
      integer i
      double precision p(n)
      double precision p_norm
      double precision r8vec_dot_product

      p_norm = 0.0D+00
      do i = 1, n
        p_norm = p_norm + p(i) * p(i)
      end do
      p_norm = sqrt ( p_norm )

      if ( p_norm .eq. 0.0D+00 ) then
        do i = 1, n
          a_normal(i) = a(i)
        end do
        do i = 1, n
          a_parallel(i) = 0.0D+00
        end do
        return
      end if

      a_dot_p = r8vec_dot_product ( n, a, p ) / p_norm

      do i = 1, n
        a_parallel(i) = a_dot_p * p(i) / p_norm
      end do

      do i = 1, n
        a_normal(i) = a(i) - a_parallel(i)
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
c    An R8VEC is a vector of R8's.
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
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine sphere_cubed_points ( n, ns, xyz )

c*********************************************************************72
c
cc SPHERE_CUBED_POINTS computes the points on a cubed sphere grid.
c
c  Discussion:
c
c    For a value of N = 3, for instance, each of the 6 cube faces will
c    be divided into 3 sections, so that a single cube face will have
c    (3+1)x(3+1) points:
c
c      X---X---X---X
c      | 1 | 4 | 7 |
c      X---X---X---X
c      | 2 | 5 | 8 |
c      X---X---X---X
c      | 3 | 6 | 9 |
c      X---X---X---X
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of sections into which each 
c    face of the cube is to be divided.
c
c    Input, integer NS, the number of points.
c
c    Output, double precision XYZ(3,NS), distinct points on the unit sphere
c    generated by a cubed sphere grid.
c
      implicit none

      integer ns

      integer n
      integer ns2
      double precision xyz(3,ns)

      ns2 = 0
c
c  Bottom full.
c
      call sphere_cubed_points_face ( n, 0, 0, 0, n, n, 0, ns2, xyz )
c
c  To avoid repetition, draw the middles as grids of n-2 x n-1 points.
c
      call sphere_cubed_points_face ( n, 0, 0, 1, 0,   n-1, n-1, 
     &  ns2, xyz )
      call sphere_cubed_points_face ( n, 0, n, 1, n-1, n,   n-1, 
     &  ns2, xyz )
      call sphere_cubed_points_face ( n, n, 1, 1, n,   n,   n-1, 
     &  ns2, xyz )
      call sphere_cubed_points_face ( n, 1, 0, 1, n,   0,   n-1, 
     &  ns2, xyz )
c
c  Top full.
c
      call sphere_cubed_points_face ( n, 0, 0, n, n, n, n, ns2, xyz )
      
      if ( ns2 .ne. ns ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPHERE_CUBED_POINTS - Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  Expected to generated NS = ', ns, ' points.'
        write ( *, '(a,i8,a)' ) '  Generated ', ns2, ' points.'
        stop
      end if

      return
      end
      subroutine sphere_cubed_points_face ( n, i1, j1, k1, i2, j2, 
     &  k2, ns, xyz )

c*********************************************************************72
c
cc SPHERE_CUBED_POINTS_FACE: points on one face of a cubed sphere grid.
c
c  Discussion:
c
c    This routine starts with NS = 0, and is called repeatedly to
c    add points for another face.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of sections into which each face 
c    of the cube is to be divided.
c
c    Input, integer I1, J1, K1, I2, J2, K2, the logical indices, 
c    between 0 and N, of two corners of the face grid.  It is guaranteed that 
c    I1 <= I2, J1 <= J2, and K1 <= K2.  
c
c    Input/output, integer NS, the number of points.
c
c    Input/output, real XYZ(3,NS), distinct points on the unit sphere
c    generated by a cubed sphere grid.
c
      implicit none

      integer i
      integer i1
      integer i2
      integer j
      integer j1
      integer j2
      integer k
      integer k1
      integer k2
      integer n
      integer ns
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision xyz(3,*)
      double precision xyzn
      double precision xc
      double precision yc
      double precision zc

      do i = i1, i2

        if ( i1 .lt. i2 ) then
          xc = tan ( dble ( 2 * i - n ) * 0.25D+00 * pi / dble ( n ) )
        else if ( i1 .eq. 0 ) then
          xc = -1.0D+00
        else if ( i1 .eq. n ) then
          xc = +1.0D+00
        else
          xc = 0.0D+00
        end if

        do j = j1, j2

          if ( j1 .lt. j2 ) then
            yc = tan ( dble ( 2 * j - n ) * 0.25D+00 * pi / dble ( n ) )
          else if ( j1 .eq. 0 ) then
            yc = -1.0D+00
          else if ( j1 .eq. n ) then
            yc = +1.0D+00
          else
            yc = 0.0D+00
          end if

          do k = k1, k2

            if ( k1 .lt. k2 ) then
              zc = tan ( dble ( 2 * k - n ) * 0.25D+00 * pi 
     &          / dble ( n ) )
            else if ( k1 .eq. 0 ) then
              zc = -1.0D+00
            else if ( k1 .eq. n ) then
              zc = +1.0D+00
            else
              zc = 0.0D+00
            end if

            xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )

            ns = ns + 1
            xyz(1,ns) = xc / xyzn
            xyz(2,ns) = yc / xyzn
            xyz(3,ns) = zc / xyzn

          end do
        end do
      end do

      return
      end
      subroutine sphere_cubed_points_size ( n, ns )

c*********************************************************************72
c
cc SPHERE_CUBED_POINTS_SIZE counts the points on a cubed sphere grid.
c
c  Discussion:
c
c    For a value of N = 3, for instance, each of the 6 cube faces will
c    be divided into 3 sections, so that a single cube face will have
c    (3+1)x(3+1) points:
c
c      X---X---X---X
c      | 1 | 4 | 7 |
c      X---X---X---X
c      | 2 | 5 | 8 |
c      X---X---X---X
c      | 3 | 6 | 9 |
c      X---X---X---X
c
c    The number of points is simply (N+1)^3 - (N-1)^3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of sections into which 
c    each face of the cube is to be divided.
c
c    Output, integer NS, the number of points.
c
      implicit none

      integer n
      integer ns

      ns = ( n + 1 ) ** 3 - ( n - 1 ) ** 3

      return
      end
      subroutine sphere_distance_xyz ( xyz1, xyz2, dist )

c*********************************************************************72
c
cc SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
c
c  Discussion:
c
c    XYZ coordinates are used.
c
c    We assume the points XYZ1 and XYZ2 lie on the same sphere.
c
c    This computation is a special form of the Vincenty formula.
c    It should be less sensitive to errors associated with very small 
c    or very large angular separations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    "Great-circle distance",
c    Wikipedia.
c
c  Parameters:
c
c    Input, double precision XYZ1(3), the coordinates of the first point.
c
c    Input, double precision XYZ2(3), the coordinates of the second point.
c
c    Output, double precision DIST, the great circle distance between
c    the points.
c
      implicit none

      double precision arc_sine
      double precision atan4
      double precision bot
      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r
      double precision r8vec_norm
      double precision top
      double precision xyz1(3)
      double precision xyz2(3)

      r = r8vec_norm ( 3, xyz1 )

      lat1 = arc_sine ( xyz1(3) )
      lon1 = atan4 ( xyz1(2), xyz1(1) )

      lat2 = arc_sine ( xyz2(3) )
      lon2 = atan4 ( xyz2(2), xyz2(1) )

      top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 
     &    + ( cos ( lat1 ) * sin ( lat2 ) 
     &    -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

      top = sqrt ( top )

      bot = sin ( lat1 ) * sin ( lat2 ) 
     &    + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

      dist = r * atan2 ( top, bot )

      return
      end
      subroutine sphere_grid_icos_size ( factor, node_num, edge_num, 
     &  triangle_num )

c*********************************************************************72
c
cc SPHERE_GRID_ICOS_SIZE sizes an icosahedral grid on a sphere.
c
c  Discussion:
c
c    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
c
c    With FACTOR = 2, each triangle of the icosahedron is subdivided into
c    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
c    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
c
c    With FACTOR = 3, each triangle of the icosahedron is subdivided into
c    3x3 subtriangles, resulting in 180 faces, 270 edges and 
c    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
c
c    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
c    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
c      12 
c    + 20 * 3          * (FACTOR-1) / 2 
c    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR, the subdivision factor, which must
c    be at least 1.
c
c    Output, integer NODE_NUM, the number of nodes.
c
c    Output, integer EDGE_NUM, the number of edges.
c
c    Output, integer TRIANGLE_NUM, the number of triangles.
c
      implicit none

      integer edge_num
      integer factor
      integer node_num
      integer triangle_num

      if ( factor .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPHERE_GRID_ICOS_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Input FACTOR < 1.'
        stop
      end if

      node_num = 12                                   
     &         + 10 * 3              * ( factor - 1 ) 
     &         + 10 * ( factor - 2 ) * ( factor - 1 )

      edge_num = 30 * factor * factor

      triangle_num = 20 * factor * factor

      return
      end
      subroutine sphere_grid_q4 ( lat_num, long_num, rectangle_node )

c*********************************************************************72
c
cc SPHERE_GRID_Q4: rectangular grid on a sphere.
c
c  Discussion:
c
c    The point numbering system is the same used in SPHERE_GRIDPOINTS,
c    and that routine may be used to compute the coordinates of the points.
c
c    A sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, the number of "rows" of rectangles to
c    be created.  LAT_NUM must be at least 2. 
c
c    Input, integer LONG_NUM, the number of "columns" of 
c    rectangles to be created.
c
c    Output, integer RECTANGLE_NODE(4,LAT_NUM*LONG_NUM), 
c    the indices of the nodes that make up each rectangle.
c
      implicit none

      integer lat_num
      integer long_num

      integer i
      integer j
      integer n
      integer n_max
      integer n_min
      integer ne
      integer nw
      integer s
      integer s_max
      integer s_min
      integer se
      integer sw
      integer rectangle_node(4,lat_num*long_num)
      integer rectangle_num

      rectangle_num = 0
c
c  The first row.
c
      n = 1

      sw = 2
      se = sw + 1

      s_min = 2
      s_max = long_num + 1

      do j = 1, long_num

        rectangle_num = rectangle_num + 1
        rectangle_node(1,rectangle_num) = sw
        rectangle_node(2,rectangle_num) = se
        rectangle_node(3,rectangle_num) = n
        rectangle_node(4,rectangle_num) = n

        sw = se

        if ( se .eq. s_max ) then
          se = s_min
        else
          se = se + 1
        end if

      end do
c
c  The intermediate rows.
c
      do i = 2, lat_num - 1

        n_max = s_max
        n_min = s_min

        s_max = s_max + long_num
        s_min = s_min + long_num

        nw = n_min
        ne = nw + 1
        sw = s_min
        se = sw + 1

        do j = 1, long_num

          rectangle_num = rectangle_num + 1
          rectangle_node(1,rectangle_num) = sw
          rectangle_node(2,rectangle_num) = se
          rectangle_node(3,rectangle_num) = ne
          rectangle_node(4,rectangle_num) = nw

          sw = se
          nw = ne

          if ( se .eq. s_max ) then
            se = s_min
          else
            se = se + 1
          end if

          if ( ne .eq. n_max ) then
            ne = n_min
          else
            ne = ne + 1
          end if

        end do

      end do
c
c  The last row.
c
      n_max = s_max
      n_min = s_min

      s = n_max + 1

      nw = n_min
      ne = nw + 1

      do j = 1, long_num

        rectangle_num = rectangle_num + 1
        rectangle_node(1,rectangle_num) = ne
        rectangle_node(2,rectangle_num) = nw
        rectangle_node(3,rectangle_num) = s
        rectangle_node(4,rectangle_num) = s

        nw = ne

        if ( ne .eq. n_max ) then
          ne = n_min
        else
          ne = ne + 1
        end if

      end do

      return
      end
      subroutine sphere_grid_t3 ( lat_num, long_num, triangle_node )

c*********************************************************************72
c
cc SPHERE_GRID_T3 produces a triangle grid on a sphere.
c
c  Discussion:
c
c    The point numbering system is the same used in SPHERE_GRIDPOINTS,
c    and that routine may be used to compute the coordinates of the points.
c
c    A sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude 
c    and longitude lines to draw.  The latitudes do not include the North 
c    and South poles, which will be included automatically, so LAT_NUM = 5, 
c    for instance, will result in points along 7 lines of latitude.
c
c    Output, integer TRIANGLE_NODE(3,2*(LAT_NUM+1)*LONG_NUM), the
c    triangle vertices.
c
      implicit none

      integer lat_num
      integer long_num

      integer i
      integer j
      integer n
      integer n_max
      integer n_min
      integer ne
      integer nw
      integer s
      integer s_max
      integer s_min
      integer se
      integer sw
      integer triangle_node(3,2*(lat_num+1)*long_num)
      integer triangle_num

      triangle_num = 0
c
c  The first row.
c
      n = 1

      sw = 2
      se = sw + 1

      s_min = 2
      s_max = long_num + 1

      do j = 0, long_num - 1

        triangle_num = triangle_num + 1
        triangle_node(1,triangle_num) = sw
        triangle_node(2,triangle_num) = se
        triangle_node(3,triangle_num) = n

        sw = se

        if ( se .eq. s_max ) then
          se = s_min
        else
          se = se + 1
        end if

      end do
c
c  The intermediate rows.
c
      do i = 1, lat_num

        n_max = s_max
        n_min = s_min

        s_max = s_max + long_num
        s_min = s_min + long_num

        nw = n_min
        ne = nw + 1
        sw = s_min
        se = sw + 1

        do j = 0, long_num - 1

          triangle_num = triangle_num + 1
          triangle_node(1,triangle_num) = sw
          triangle_node(2,triangle_num) = se
          triangle_node(3,triangle_num) = nw

          triangle_num = triangle_num + 1
          triangle_node(1,triangle_num) = ne
          triangle_node(2,triangle_num) = nw
          triangle_node(3,triangle_num) = se

          sw = se
          nw = ne

          if ( se .eq. s_max ) then
            se = s_min
          else
            se = se + 1
          end if

          if ( ne .eq. n_max ) then
            ne = n_min
          else
            ne = ne + 1
          end if

        end do

      end do
c
c  The last row.
c
      n_max = s_max
      n_min = s_min

      s = n_max + 1

      nw = n_min
      ne = nw + 1

      do j = 0, long_num - 1

        triangle_num = triangle_num + 1
        triangle_node(1,triangle_num) = ne
        triangle_node(2,triangle_num) = nw
        triangle_node(3,triangle_num) = s

        nw = ne

        if ( ne .eq. n_max ) then
          ne = n_min
        else
          ne = ne + 1
        end if

      end do

      return
      end
      subroutine sphere_gridpoints_icos1 ( factor, node_num, node_xyz )

c*********************************************************************72
c
cc SPHERE_GRIDPOINTS_ICOS1 returns icosahedral grid points on a sphere.
c
c  Discussion:
c
c    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
c
c    With FACTOR = 2, each triangle of the icosahedron is subdivided into
c    2x2 subtriangles, resulting in 80 faces and 
c    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
c
c    With FACTOR = 3, each triangle of the icosahedron is subdivided into
c    3x3 subtriangles, resulting in 180 faces and 
c    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
c
c    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
c    resulting in 20 * FACTOR * FACTOR faces and
c      12 
c    + 20 * 3          * (FACTOR-1) / 2 
c    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
c
c
c    This routine uses a simple, but only approximate, method of
c    carrying out the subdivision.  For each spherical triangle of the
c    face, we actually work in the planar triangle defined by those
c    three points.  All subdivisions are done on that planar triangle,
c    and the resulting points are then projected onto the sphere.
c    While these points are equally spaced on the planar triangle,
c    that is only approximately true on the sphere.
c
c    See SPHERE_GRID_POINTS_ICOS2 for a more accurate method of subdivision
c    that works on the spherical triangle itself.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR, the subdivision factor, which must
c    be at least 1.
c
c    Input, integer NODE_NUM, the number of nodes, as reported
c    by SPHERE_GRID_ICOS_SIZE.
c
c    Output, double precision NODE_XYZ(3,NODE_NUM), the node coordinates.
c
c  Local Parameters:
c
c    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
c    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
c    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
c    We need to refer to this data to generate the grid.
c
c    NODE counts the number of nodes we have generated so far.  At the
c    end of the routine, it should be equal to NODE_NUM.
c
      implicit none

      integer edge_num
      parameter ( edge_num = 30 )
      integer face_num
      parameter ( face_num = 20 )
      integer face_order_max
      parameter ( face_order_max = 3 )
      integer node_num
      integer point_num
      parameter ( point_num = 12 )

      integer a
      integer b
      integer c
      integer edge
      integer edge_point(2,edge_num)
      integer f
      integer f1
      integer f2
      integer face
      integer face_order(face_num)
      integer face_point(face_order_max,face_num)
      integer factor
      integer i
      integer j
      integer node
      double precision node_norm
      double precision node_xyz(3,node_num)
      double precision point_coord(3,point_num)
      double precision r8vec_norm

      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Generate the point coordinates.
c
c  A.  Points that are the icosahedral vertices.
c
      node = 0
      do j = 1, point_num
        do i = 1, 3
          node_xyz(i,j) = point_coord(i,j)
        end do
      end do
c
c  B. Points in the icosahedral edges, at 
c  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
c
      node = 12

      do edge = 1, edge_num

        a = edge_point(1,edge)
        b = edge_point(2,edge)

        do f = 1, factor - 1

          node = node + 1

          do i = 1, 3
            node_xyz(i,node) = 
     &        ( dble ( factor - f ) * point_coord(i,a)   
     &        + dble (          f ) * point_coord(i,b) ) 
     &        / dble ( factor     )
          end do

          node_norm = r8vec_norm ( 3, node_xyz(1,node) )

          do i = 1, 3
            node_xyz(i,node) = node_xyz(i,node) / node_norm
          end do

        end do
      end do
c
c  C.  Points in the icosahedral faces.
c
      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do f1 = 1, factor - 1
          do f2 = 1, factor - f1 - 1

            node = node + 1

            do i = 1, 3
              node_xyz(i,node) = 
     &          ( dble ( factor - f1 - f2 ) * point_coord(i,a)   
     &          + dble (          f1      ) * point_coord(i,b)   
     &          + dble (               f2 ) * point_coord(i,c) ) 
     &          / dble ( factor           )
            end do

            node_norm = r8vec_norm ( 3, node_xyz(1:3,node) )

            do i = 1, 3
              node_xyz(i,node) = node_xyz(i,node) / node_norm
            end do

          end do
        end do

      end do

      return
      end
      subroutine sphere_gridpoints_icos2 ( factor, node_num, node_xyz )

c*********************************************************************72
c
cc SPHERE_GRIDPOINTS_ICOS2 returns icosahedral grid points on a sphere.
c
c  Discussion:
c
c    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
c
c    With FACTOR = 2, each triangle of the icosahedron is subdivided into
c    2x2 subtriangles, resulting in 80 faces and 
c    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
c
c    With FACTOR = 3, each triangle of the icosahedron is subdivided into
c    3x3 subtriangles, resulting in 180 faces and 
c    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
c
c    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
c    resulting in 20 * FACTOR * FACTOR faces and
c      12 
c    + 20 * 3          * (FACTOR-1) / 2 
c    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR, the subdivision factor, which must
c    be at least 1.
c
c    Input, integer NODE_NUM, the number of nodes, as reported
c    by SPHERE_GRID_ICOS_SIZE.
c
c    Output, double precision NODE_XYZ(3,NODE_NUM), the node coordinates.
c
c  Local Parameters:
c
c    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
c    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
c    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
c    We need to refer to this data to generate the grid.
c
c    NODE counts the number of nodes we have generated so far.  At the
c    end of the routine, it should be equal to NODE_NUM.
c
      implicit none

      integer edge_num
      parameter ( edge_num = 30 )
      integer face_num
      parameter ( face_num = 20 )
      integer face_order_max
      parameter ( face_order_max = 3 )
      integer node_num
      integer point_num
      parameter ( point_num = 12 )

      integer a
      double precision angle
      double precision ab(3)
      double precision ac(3)
      double precision acn(3)
      double precision acn_norm
      double precision acp(3)
      integer b
      double precision bn(3)
      double precision bn_norm
      double precision bp(3)
      integer c
      double precision cn(3)
      double precision cn_norm
      double precision cp(3)
      integer edge
      integer edge_point(2,edge_num)
      integer f
      integer fa
      integer fbc
      integer face
      integer face_order(face_num)
      integer face_point(face_order_max,face_num)
      integer factor
      integer i
      integer j
      integer node
      double precision node_xyz(3,node_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision point_coord(3,point_num)
      double precision r8vec_norm
      double precision theta
      double precision theta_ab
      double precision theta_ac
      double precision theta_bc

      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Generate the point coordinates.
c
c  A.  Points that are the icosahedral vertices.
c
      node = 0
      do j = 1, point_num
        do i = 1, 3
          node_xyz(i,j) = point_coord(i,j)
        end do
      end do
c
c  B. Points in the icosahedral edges, at 
c  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
c
      node = 12

      do edge = 1, edge_num

        a = edge_point(1,edge)

        b = edge_point(2,edge)
c
c  Determine the "distance" = angle between points A and B.
c
        call sphere_distance_xyz ( point_coord(1,a), point_coord(1,b), 
     &    theta )
c
c  Polarize B into BP + BN and normalize BN.
c
        call r8vec_polarize ( 3, point_coord(1,b), point_coord(1,a), 
     &    bn, bp )

        bn_norm = r8vec_norm ( 3, bn )

        do i = 1, 3
          bn(i) = bn(i) / bn_norm
        end do
c
c  March from A to B, by taking equally spaced angles from 0 to THETA.
c  F = 0      => ANGLE = 0     => A
c  F = FACTOR => ANGLE = THETA => B
c
        do f = 1, factor - 1

          node = node + 1

          angle = ( dble ( f ) * theta ) / dble ( factor )

          do i = 1, 3
            node_xyz(i,node) = cos ( angle ) * point_coord(i,a) 
     &                       + sin ( angle ) * bn(i)
          end do

        end do
      end do
c
c  C.  Points in the icosahedral faces.
c
      do face = 1, face_num

        a = face_point(1,face)

        b = face_point(2,face)

        c = face_point(3,face)
c
c  Determine the "distance" = angle between points A and B, A and C.
c
        call sphere_distance_xyz ( point_coord(1,a), point_coord(1,b), 
     &    theta_ab )

        call sphere_distance_xyz ( point_coord(1,a), point_coord(1,c), 
     &    theta_ac )
c
c  Polarize B = BP + BN and normalize BN, C = CP + CN, and normalize CN.
c
        call r8vec_polarize ( 3, point_coord(1,b), point_coord(1,a), 
     &    bn, bp )
        bn_norm = r8vec_norm ( 3, bn )
        do i = 1, 3
          bn(i) = bn(i) / bn_norm
        end do

        call r8vec_polarize ( 3, point_coord(1,c), point_coord(1,a), 
     &    cn, cp )
        cn_norm = r8vec_norm ( 3, cn )
        do i = 1, 3
          cn(i) = cn(i) / cn_norm
        end do
c
c  March AB from A to B:
c    FA = 0      => ANGLE = 0        => AB = A
c    FA = FACTOR => ANGLE = THETA_AB => AB = B
c
c  March AC from A to C:
c    FA = 0      => ANGLE = 0        => AC = A
c    FA = FACTOR => ANGLE = THETA_AC => AC = C
c
        do fa = 2, factor - 1
c
c  Determine points AB and AC that use cos ( FA / FACTOR ) of A 
c  and cos ( ( FACTOR - FA ) / FACTOR ) of B or C.
c
          angle = ( dble ( fa ) * theta_ab ) / dble ( factor )
          do i = 1, 3
            ab(i) = cos ( angle ) * point_coord(i,a) 
     &            + sin ( angle ) * bn(i)
          end do

          angle = ( dble ( fa ) * theta_ac ) / dble ( factor )
          do i = 1, 3
            ac(i) = cos ( angle ) * point_coord(i,a) 
     &            + sin ( angle ) * cn(i)
          end do
c
c  Determine the "distance" = angle between points AB and AC.
c
          call sphere_distance_xyz ( ab, ac, theta_bc )
c
c  Polarize AC into ACP + ACN and normalize ACN.
c
          call r8vec_polarize ( 3, ac, ab, acn, acp )
          acn_norm = r8vec_norm ( 3, acn )
          do i = 1, 3
            acn(i) = acn(i) / acn_norm
          end do
c
c  The interval between AB and AC is broken into FA intervals.
c  Go from 1 to FA - 1.
c
          do fbc = 1, fa - 1

            node = node + 1

            angle = ( real ( fbc ) * theta_bc ) 
     &              / real ( fa )
         
            do i = 1, 3
              node_xyz(i,node) = cos ( angle ) * ab(i) 
     &                         + sin ( angle ) * acn(i)
            end do

          end do
        end do

      end do

      return
      end
      subroutine sphere_line_project ( r, pc, n, p, maxpnt2, n2, pp, 
     &  theta_min, theta_max )

c*********************************************************************72
c
cc SPHERE_LINE_PROJECT projects a line onto a sphere.
c
c  Discussion:
c
c    A sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
c
c    The line to be projected is specified as a sequence of points.
c    If two successive points subtend a small angle, then the second
c    point is essentially dropped.  If two successive points subtend
c    a large angle, then intermediate points are inserted, so that
c    the projected line stays closer to the sphere.
c
c    Note that if any P coincides with the center of the sphere, then
c    its projection is mathematically undefined.  PP will
c    be returned as the center.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.  If R is
c    zero, PP will be returned as the pc, and if R is
c    negative, points will end up diametrically opposite from where
c    you would expect them for a positive R.
c
c    Input, double precision PC(3), the center of the sphere.
c
c    Input, integer N, the number of points on the line that is
c    to be projected.
c
c    Input, double precision P(3,N), the coordinates of
c    the points on the line that is to be projected.
c
c    Input, integer MAXPNT2, the maximum number of points on the
c    projected line.  Even if the routine thinks that more points are needed,
c    no more than MAXPNT2 will be generated.
c
c    Output, integer N2, the number of points on the projected
c    line.  N2 can be zero, if the line has an angular projection of less
c    than THETA_MIN radians.
c
c    Output, double precision PP(3,N2), the coordinates
c    of the points representing the projected line.  These points lie on the
c    sphere.  Successive points are separated by at least THETA_MIN
c    radians, and by no more than THETA_MAX radians.
c
c    Input, double precision THETA_MIN, THETA_MAX, the minimum and maximum
c    angular projections allowed between successive projected points.
c    If two successive points on the original line have projections
c    separated by more than THETA_MAX radians, then intermediate points
c    will be inserted, in an attempt to keep the line closer to the
c    sphere.  If two successive points are separated by less than
c    THETA_MIN radians, then the second point is dropped, and the
c    line from the first point to the next point is considered.
c
      implicit none

      integer maxpnt2
      integer n

      double precision alpha
      double precision ang3d
      double precision arc_cosine
      double precision dot
      integer i
      integer j
      integer k
      integer nfill
      integer n2
      double precision p(3,n)
      double precision p1(3)
      double precision p2(3)
      double precision pc(3)
      double precision pd(3)
      double precision pp(3,maxpnt2)
      double precision r
      double precision r8vec_diff_norm
      double precision r8vec_norm
      double precision theta_max
      double precision theta_min
      double precision tnorm
c
c  Check the input.
c
      if ( r .eq. 0.0D+00 ) then
        n2 = 0
        return
      end if

      do i = 1, 3
        p1(i) = pc(i)
        p2(i) = pc(i)
      end do

      n2 = 0

      do i = 1, n

        if ( p(1,i) .eq. pc(1) .and.
     &       p(2,i) .eq. pc(2) .and.
     &       p(3,i) .eq. pc(3) ) then

        else

          do j = 1, 3
            p1(j) = p2(j)
          end do

          alpha = r8vec_diff_norm ( 3, p(1,i), pc )

          do j = 1, 3
            p2(j) = pc(j) + r * ( p(j,i) - pc(j) ) / alpha
          end do
c
c  If we have not gotten any points yet, take this point as our start.
c
          if ( n2 .eq. 0 ) then

            n2 = n2 + 1

            do j = 1, 3
              pp(j,n2) = p2(j)
            end do
c
c  Compute the angular projection of P1 to P2.
c
          else if ( 1 .le. n2 ) then

            dot = 0.0D+00
            do j = 1, 3
              dot = dot + ( p1(j) - pc(j) ) * ( p2(j) - pc(j) )
            end do

            ang3d = arc_cosine (  dot / ( r * r ) )
c
c  If the angle is at least THETA_MIN, (or it is the last point),
c  then we will draw a line segment.
c
            if ( theta_min .lt. abs ( ang3d ) .or. i .eq. n ) then
c
c  Now we check to see if the line segment is too long.
c
              if ( theta_max .lt. abs ( ang3d ) ) then

                nfill = int ( abs ( ang3d ) / theta_max )

                do j = 1, nfill - 1

                  do k = 1, 3
                   pd(k) = 
     &                 ( dble ( nfill - j ) * ( p1(k) - pc(k) ) 
     &                 + dble (         j ) * ( p2(k) - pc(k) ) )
                  end do

                  tnorm = r8vec_norm ( 3, pd )

                  if ( tnorm .ne. 0.0D+00 ) then
                    do k = 1, 3
                      pd(k) = pc(k) + r * pd(k) / tnorm
                    end do
                    n2 = n2 + 1
                    do k = 1, 3
                      pp(k,n2) = pd(k)
                    end do
                  end if

                end do

              end if
c
c  Now tack on the projection of point 2.
c
              n2 = n2 + 1

              do k = 1, 3
                pp(k,n2) = p2(k)
              end do

            end if

          end if

        end if

      end do

      return
      end
      subroutine sphere_ll_lines ( lat_num, long_num, line_num, line )

c*********************************************************************72
c
cc SPHERE_LL_LINES produces lines for a latitude/longitude grid.
c
c  Discussion:
c
c    The point numbering system is the same used in SPHERE_LL_POINTS,
c    and that routine may be used to compute the coordinates of the points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
c    longitude lines to draw.  The latitudes do not include the North and South
c    poles, which will be included automatically, so LAT_NUM = 5, for instance,
c    will result in points along 7 lines of latitude.
c
c    Input, integer LINE_NUM, the number of grid lines.
c
c    Output, integer LINE(2,LINE_NUM), contains pairs of point 
c    indices for line segments that make up the grid.
c
      implicit none

      integer line_num

      integer i
      integer j
      integer lat_num
      integer l
      integer line(2,line_num)
      integer long_num
      integer new
      integer newcol
      integer old

      l = 0
c
c  "Vertical" lines.
c
      do j = 0, long_num - 1

        old = 1
        new = j + 2

        l = l + 1
        line(1,l) = old
        line(2,l) = new

        do i = 1, lat_num - 1

          old = new
          new = old + long_num

          l = l + 1
          line(1,l) = old
          line(2,l) = new

        end do

        old = new

        l = l + 1
        line(1,l) = old
        line(2,l) = 1 + lat_num * long_num + 1

      end do
c
c  "Horizontal" lines.
c
      do i = 1, lat_num

        new = 1 + ( i - 1 ) * long_num + 1

        do j = 0, long_num - 2
          old = new
          new = old + 1
          l = l + 1
          line(1,l) = old
          line(2,l) = new
        end do

        old = new
        new = 1 + ( i - 1 ) * long_num + 1
        l = l + 1
        line(1,l) = old
        line(2,l) = new

      end do
c
c  "Diagonal" lines.
c
      do j = 0, long_num - 1

        old = 1
        new = j + 2
        newcol = j

        do i = 1, lat_num - 1

          old = new
          new = old + long_num + 1

          newcol = newcol + 1
          if ( long_num - 1 .lt. newcol ) then
            newcol = 0
            new = new - long_num
          end if

          l = l + 1
          line(1,l) = old
          line(2,l) = new

        end do

      end do

      return
      end
      subroutine sphere_ll_line_num ( lat_num, long_num, line_num )

c*********************************************************************72
c
cc SPHERE_LL_LINE_NUM counts lines for a latitude/longitude grid.
c
c  Discussion:
c
c    The number returned is the number of pairs of points to be connected.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
c    longitude lines to draw.  The latitudes do not include the North and South
c    poles, which will be included automatically, so LAT_NUM = 5, for instance,
c    will result in points along 7 lines of latitude.
c
c    Output, integer LINE_NUM, the number of grid lines.
c
      implicit none

      integer lat_num
      integer line_num
      integer long_num

      line_num = long_num * ( lat_num + 1 ) 
     &         + lat_num * long_num 
     &         + long_num * ( lat_num - 1 )

      return
      end
      subroutine sphere_ll_points ( r, pc, lat_num, long_num, 
     &  point_num, p )

c*********************************************************************72
c
cc SPHERE_LL_POINTS produces points on a latitude/longitude grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision PC(3), the center of the sphere.
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude 
c    and longitude lines to draw.  The latitudes do not include the North and 
c    South poles, which will be included automatically, so LAT_NUM = 5, for 
c    instance, will result in points along 7 lines of latitude.
c
c    Input, integer POINT_NUM, the number of points.
c
c    Output, double precision P(3,POINT_NUM), the grid points.
c
      implicit none

      integer lat_num
      integer long_num
      integer point_num

      integer lat
      integer long
      integer n
      double precision p(3,point_num)
      double precision pc(3)
      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      n = 0
c
c  The north pole.
c
      theta = 0.0D+00
      phi = 0.0D+00
      n = n + 1
      p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
      p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
      p(3,n) = pc(3) + r * cos ( phi )
c
c  Do each intermediate ring of latitude.
c
      do lat = 1, lat_num

        phi = dble ( lat         ) * pi 
     &      / dble ( lat_num + 1 )
c
c  Along that ring of latitude, compute points at various longitudes.
c
        do long = 0, long_num - 1

          theta = dble ( long     ) * 2.0D+00 * pi 
     &          / dble ( long_num )

          n = n + 1
          p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
          p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
          p(3,n) = pc(3) + r * cos ( phi )

        end do
      end do
c
c  The south pole.
c
      theta = 0.0D+00
      phi = pi
      n = n + 1
      p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
      p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
      p(3,n) = pc(3) + r * cos ( phi )

      return
      end
      subroutine sphere_ll_point_num ( lat_num, long_num, point_num )

c*********************************************************************72
c
cc SPHERE_LL_POINT_NUM counts points for a latitude/longitude grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude 
c    and longitude lines to draw.  The latitudes do not include the North and 
c    South poles, which will be included automatically, so LAT_NUM = 5, for 
c    instance, will result in points along 7 lines of latitude.
c
c    Output, integer ( kind = 4 ) POINT_NUM, the number of grid points.
c
      implicit none

      integer lat_num
      integer long_num
      integer point_num

      point_num = 2 + lat_num * long_num

      return
      end
      subroutine sphere_llq_lines ( lat_num, long_num, line_num, line )

c*********************************************************************72
c
cc SPHERE_LLQ_LINES: latitude/longitude quadrilateral grid lines.
c
c  Discussion:
c
c    The point numbering system is the same used in SPHERE_LL_POINTS,
c    and that routine may be used to compute the coordinates of the points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
c    longitude lines to draw.  The latitudes do not include the North and South
c    poles, which will be included automatically, so LAT_NUM = 5, for instance,
c    will result in points along 7 lines of latitude.
c
c    Input, integer LINE_NUM, the number of grid lines.
c
c    Output, integer LINE(2,LINE_NUM), contains pairs of point 
c    indices for line segments that make up the grid.
c
      implicit none

      integer line_num

      integer i
      integer j
      integer lat_num
      integer l
      integer line(2,line_num)
      integer long_num
      integer new
      integer newcol
      integer old

      l = 0
c
c  "Vertical" lines.
c
      do j = 0, long_num - 1

        old = 1
        new = j + 2

        l = l + 1
        line(1,l) = old
        line(2,l) = new

        do i = 1, lat_num - 1

          old = new
          new = old + long_num

          l = l + 1
          line(1,l) = old
          line(2,l) = new

        end do

        old = new

        l = l + 1
        line(1,l) = old
        line(2,l) = 1 + lat_num * long_num + 1

      end do
c
c  "Horizontal" lines.
c
      do i = 1, lat_num

        new = 1 + ( i - 1 ) * long_num + 1

        do j = 0, long_num - 2
          old = new
          new = old + 1
          l = l + 1
          line(1,l) = old
          line(2,l) = new
        end do

        old = new
        new = 1 + ( i - 1 ) * long_num + 1
        l = l + 1
        line(1,l) = old
        line(2,l) = new

      end do

      return
      end
      subroutine sphere_llq_line_num ( lat_num, long_num, line_num )

c*********************************************************************72
c
cc SPHERE_LLQ_LINE_NUM counts lines for a latitude/longitude quadrilateral grid.
c
c  Discussion:
c
c    The number returned is the number of pairs of points to be connected.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
c    longitude lines to draw.  The latitudes do not include the North and South
c    poles, which will be included automatically, so LAT_NUM = 5, for instance,
c    will result in points along 7 lines of latitude.
c
c    Output, integer LINE_NUM, the number of grid lines.
c
      implicit none

      integer lat_num
      integer line_num
      integer long_num

      line_num = long_num * ( lat_num + 1 ) 
     &         + lat_num * long_num

      return
      end
      subroutine sphere_spiralpoints ( r, pc, n, p )

c*********************************************************************72
c
cc SPHERE_SPIRALPOINTS: spiral points on a sphere.
c
c  Discussion:
c
c    The points should be arranged on the sphere in a pleasing design.
c
c    A sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Saff, Arno Kuijlaars,
c    Distributing Many Points on a Sphere,
c    The Mathematical Intelligencer,
c    Volume 19, Number 1, 1997, pages 5-11.
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision PC(3), the center of the sphere.
c
c    Input, integer N, the number of points to create.
c
c    Output, double precision P(3,N), the grid points.
c
      implicit none

      integer n

      double precision cosphi
      integer i
      double precision p(3,n)
      double precision pc(3)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision sinphi
      double precision theta

      do i = 1, n

        cosphi = ( dble ( n - i     ) * ( -1.0D+00 ) 
     &           + dble (     i - 1 ) * ( +1.0D+00 ) ) 
     &           / dble ( n     - 1 )

        sinphi = sqrt ( 1.0D+00 - cosphi * cosphi )

        if ( i == 1 .or. i == n ) then
          theta = 0.0D+00
        else
          theta = theta + 3.6D+00 / ( sinphi * sqrt ( dble ( n ) ) )
          theta = mod ( theta, 2.0D+00 * pi )
        end if

        p(1,i) = pc(1) + r * sinphi * cos ( theta )
        p(2,i) = pc(2) + r * sinphi * sin ( theta )
        p(3,i) = pc(3) + r * cosphi

      end do

      return
      end
      subroutine sphere_unit_sample ( n, seed, x )

c*********************************************************************72
c
cc SPHERE_UNIT_SAMPLE picks a random point on the unit sphere.
c
c  Discussion:
c
c    The unit sphere in 3D satisfies:
c
c      X * X + Y * Y + Z * Z = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of samples.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(3,N), the sample points.
c
      implicit none

      integer n

      double precision arc_cosine
      integer j
      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision theta
      double precision vdot
      double precision x(3,n)

      do j = 1, n
c
c  Pick a uniformly random VDOT, which must be between -1 and 1.
c  This represents the dot product of the random vector with the Z unit vector.
c
c  Note: this works because the surface area of the sphere between
c  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
c  a patch of area uniformly.
c
        vdot = r8_uniform_01 ( seed )
        vdot = 2.0D+00 * vdot - 1.0D+00

        phi = arc_cosine ( vdot )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the Z vector.
c
        theta = r8_uniform_01 ( seed )
        theta = 2.0D+00 * pi * theta

        x(1,j) = cos ( theta ) * sin ( phi )
        x(2,j) = sin ( theta ) * sin ( phi )
        x(3,j) = cos ( phi )

      end  do

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
