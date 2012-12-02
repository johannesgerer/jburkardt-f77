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
      subroutine icos_shape ( point_num, edge_num, face_num, 
     &  face_order_max, point_coord, edge_point, face_order, 
     &  face_point )

c*********************************************************************72
c
cc ICOS_SHAPE describes an icosahedron.
c
c  Discussion:
c
c    The input data required for this routine can be retrieved from
c    ICOS_SIZE.
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
c    23 September 2010
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
      integer dim_num
      parameter ( dim_num = 3 )
      integer point_num

      integer edge_point(edge_order,edge_num)
      integer edge_point_save(2,30)
      integer face_order(face_num)
      integer face_point(face_order_max,face_num)
      integer face_point_save(3,20)
      integer i
      integer j
      double precision point_coord(dim_num,point_num)
      double precision point_coord_save(3,12)

      save edge_point_save
      save face_point_save
      save point_coord_save

      data edge_point_save /
     &   1,  2, 
     &   1,  3, 
     &   1,  4,
     &   1,  5,
     &   1,  6, 
     &   2,  3, 
     &   2,  4, 
     &   2,  7, 
     &   2,  8, 
     &   3,  5, 
     &   3,  7, 
     &   3,  9, 
     &   4,  6, 
     &   4,  8, 
     &   4, 10, 
     &   5,  6, 
     &   5,  9, 
     &   5, 11, 
     &   6, 10, 
     &   6, 11, 
     &   7,  8, 
     &   7,  9, 
     &   7, 12, 
     &   8, 10, 
     &   8, 12, 
     &   9, 11, 
     &   9, 12, 
     &  10, 11, 
     &  10, 12, 
     &  11, 12 /

      data face_point_save /
     &   1,  2,  4, 
     &   1,  3,  2, 
     &   1,  4,  6, 
     &   1,  5,  3, 
     &   1,  6,  5, 
     &   2,  3,  7, 
     &   2,  7,  8, 
     &   2,  8,  4, 
     &   3,  5,  9, 
     &   3,  9,  7, 
     &   4,  8, 10, 
     &   4, 10,  6, 
     &   5,  6, 11, 
     &   5, 11,  9, 
     &   6, 10, 11, 
     &   7,  9, 12, 
     &   7, 12,  8, 
     &   8, 12, 10, 
     &   9, 11, 12, 
     &  10, 12, 11 /

      data point_coord_save /
     &    0.85065080835203999D+00,  0.52573111211913359D+00,  0.0D+00, 
     &    0.85065080835203999D+00, -0.52573111211913359D+00,  0.0D+00, 
     &    0.52573111211913359D+00,  0.0D+00,  0.85065080835203999D+00, 
     &    0.52573111211913359D+00,  0.0D+00, -0.85065080835203999D+00, 
     &    0.0D+00,  0.85065080835203999D+00,  0.52573111211913359D+00, 
     &    0.0D+00,  0.85065080835203999D+00, -0.52573111211913359D+00, 
     &    0.0D+00, -0.85065080835203999D+00,  0.52573111211913359D+00, 
     &    0.0D+00, -0.85065080835203999D+00, -0.52573111211913359D+00, 
     &   -0.52573111211913359D+00,  0.0D+00,  0.85065080835203999D+00, 
     &   -0.52573111211913359D+00,  0.0D+00, -0.85065080835203999D+00, 
     &   -0.85065080835203999D+00,  0.52573111211913359D+00,  0.0D+00, 
     &   -0.85065080835203999D+00, -0.52573111211913359D+00,  0.0D+00 /

      do j = 1, edge_num
        do i = 1, edge_order
          edge_point(i,j) = edge_point_save(i,j)
        end do
      end do

      do i = 1, face_num
        face_order(i) = 3
      end do

      do j = 1, face_num
        do i = 1, face_order_max
          face_point(i,j) = face_point_save(i,j)
        end do
      end do

      do j = 1, point_num
        do i = 1, dim_num
          point_coord(i,j) = point_coord_save(i,j)
        end do
      end do

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
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates Gamma(X) for a real argument.
c
c  Discussion:
c
c    This routine calculates the gamma function for a real argument X.
c    Computation is based on an algorithm outlined in reference 1.
c    The program uses rational functions that approximate the gamma
c    function to at least 20 significant decimal digits.  Coefficients
c    for the approximation over the interval (1,2) are unpublished.
c    Those for the approximation for 12 <= X are from reference 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    An Overview of Software Development for Special Functions,
c    in Numerical Analysis Dundee, 1975,
c    edited by GA Watson,
c    Lecture Notes in Mathematics 506,
c    Springer, 1976.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      integer i
      integer n
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
c
c  Mathematical constants
c
      data sqrtpi /0.9189385332046727417803297D+00/
      data pi /3.1415926535897932384626434D+00/
c
c  Machine dependent parameters
c
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (1,2).
c
      data p/
     & -1.71618513886549492533811d+00,
     &  2.47656508055759199108314d+01,
     & -3.79804256470945635097577d+02,
     &  6.29331155312818442661052d+02,
     &  8.66966202790413211295064d+02,
     & -3.14512729688483675254357d+04,
     & -3.61444134186911729807069d+04,
     &  6.64561438202405440627855d+04/

      data q/
     & -3.08402300119738975254353d+01,
     &  3.15350626979604161529144d+02,
     & -1.01515636749021914166146d+03,
     & -3.10777167157231109440444d+03,
     &  2.25381184209801510330112d+04,
     &  4.75584627752788110767815d+03,
     & -1.34659959864969306392456d+05,
     & -1.15132259675553483497211d+05/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     & -1.910444077728D-03,
     &  8.4171387781295D-04,
     & -5.952379913043012D-04,
     &  7.93650793500350248D-04,
     & -2.777777777777681622553D-03,
     &  8.333333333333333331554247D-02,
     &  5.7083835261D-03/

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
c
c  Argument is negative.
c
      if ( y .le. 0.0D+00 ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. 0.0D+00 ) then

          if ( y1 .ne. aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + 1.0D+00

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Argument is positive.
c
      if ( y .lt. eps ) then
c
c  Argument < EPS.
c
        if ( xminin .le. y ) then
          res = 1.0D+00 / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. 12.0D+00 ) then

        y1 = y
c
c  0.0 < argument < 1.0.
c
        if ( y .lt. 1.0D+00 ) then

          z = y
          y = y + 1.0D+00
c
c  1.0 < argument < 12.0.
c  Reduce argument if necessary.
c
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - 1.0D+00

        end if
c
c  Evaluate approximation for 1.0 < argument < 2.0.
c
        xnum = 0.0D+00
        xden = 1.0D+00
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + 1.0D+00
c
c  Adjust result for case  0.0 < argument < 1.0.
c
        if ( y1 .lt. y ) then

          res = res / y1
c
c  Adjust result for case 2.0 < argument < 12.0.
c
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + 1.0D+00
          end do

        end if

      else
c
c  Evaluate for 12.0 <= argument.
c
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - 0.5D+00 ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Final adjustments and return.
c
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. 1.0D+00 ) then
        res = fact / res
      end if

      r8_gamma = res

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
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
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
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
      subroutine sphere01_distance_xyz ( xyz1, xyz2, dist )

c*********************************************************************72
c
cc SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
c
c  Discussion:
c
c    XYZ coordinates are used.
c
c    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
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
c    23 September 2010
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
      double precision top
      double precision xyz1(3)
      double precision xyz2(3)

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

      dist = atan2 ( top, bot )

      return
      end
      subroutine sphere01_monomial_int ( e, integral )

c*********************************************************************72
c
cc SPHERE01_MONOMIAL_INT integrates a monomial on the unit sphere.
c
c  Discussion:
c
c    The integration region is 
c
c      X^2 + Y^2 + Z^2 = 1.
c
c    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Academic Press, 1984, page 263.
c
c  Parameters:
c
c    Input, integer E(3), the exponents of X, Y and Z in the 
c    monomial.  Each exponent must be nonnegative.
c
c    Output, double precision INTEGRAL, the integral.
c
      implicit none

      integer e(3)
      integer i
      double precision integral
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_gamma

      if ( e(1) .lt. 0 .or. e(2) .lt. 0 .or. e(3) .lt. 0 ) then
        integral = - huge ( integral )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INT - Fatal error!'
        write ( *, '(a)' ) '  All exponents must be nonnegative.'
        write ( *, '(a,i8)' ) '  E(1) = ', e(1)
        write ( *, '(a,i8)' ) '  E(2) = ', e(2)
        write ( *, '(a,i8)' ) '  E(3) = ', e(3)
        stop
      end if

      if ( e(1) .eq. 0 .and. e(2) .eq. 0 .and. e(3) .eq. 0 ) then

        integral = 2.0D+00 * sqrt ( pi**3 ) / r8_gamma ( 1.5D+00 )

      else if ( mod ( e(1), 2 ) .eq. 1 .or.
     &          mod ( e(2), 2 ) .eq. 1 .or.
     &          mod ( e(3), 2 ) .eq. 1 ) then

        integral = 0.0D+00

      else

        integral = 2.0D+00

        do i = 1, 3
          integral = integral * r8_gamma ( 0.5D+00 * dble ( e(i) + 1 ) )
        end do

        integral = integral 
     &    / r8_gamma ( 0.5D+00 * dble ( e(1) + e(2) + e(3) + 3 ) )

      end if

      return
      end
      subroutine sphere01_quad_icos1c ( factor, fun, node_num, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
c
c  Discussion:
c
c    This function estimates an integral over the surface of the unit sphere.
c
c    This function sets up an icosahedral grid, and subdivides each
c    edge of the icosahedron into FACTOR subedges.  These edges define a grid
c    within each triangular icosahedral face.  The centroids of these
c    triangles can be determined.  All of these calculations are done,
c    essentially, on the FLAT faces of the icosahedron.  Only then are
c    the triangle vertices and centroids projected to the sphere.  
c
c    The resulting grid of spherical triangles and projected centroids
c    is used to apply a centroid quadrature rule over the surface of
c    the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 September 2010
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
c    Input, external :: FUN, evaluates the integrand, of the form:
c      subroutine fun ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Output, integer NODE_NUM, the number of evaluation points.
c
c    Output, double precision RESULT, the estimated integral.
c
      implicit none

      integer node_num

      integer a
      double precision a_xyz(3)
      double precision a2_xyz(3)
      double precision area
      double precision area_total
      integer b
      double precision b_xyz(3)
      double precision b2_xyz(3)
      integer c
      double precision c_xyz(3)
      double precision c2_xyz(3)
      integer edge_num
      integer edge_point(2,30)
      integer f1
      integer f2
      integer f3
      integer face
      integer face_num
      integer face_order(20)
      integer face_point(3,20)
      integer face_order_max
      integer factor
      external fun
      integer i
      double precision node_xyz(3)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision point_coord(3,12)
      integer point_num
      double precision result
      double precision v
c
c  Size the icosahedron.
c
      call icos_size ( point_num, edge_num, face_num, face_order_max )
c
c  Set the icosahedron.
c
      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Initialize the integral data.
c
      result = 0.0D+00
      area_total = 0.0D+00
      node_num = 0
c
c  Pick a face of the icosahedron, and identify its vertices as A, B, C.
c
      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do i = 1, 3
          a_xyz(i) = point_coord(i,a)
          b_xyz(i) = point_coord(i,b)
          c_xyz(i) = point_coord(i,c)
        end do
c
c  Some subtriangles will have the same direction as the face.
c  Generate each in turn, by determining the barycentric coordinates
c  of the centroid (F1,F2,F3), from which we can also work out the barycentric
c  coordinates of the vertices of the subtriangle.
c
        do f3 = 1, 3 * factor - 2, 3
          do f2 = 1, 3 * factor - f3 - 1, 3

            f1 = 3 * factor - f3 - f2

            call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
     &        f2, f3, node_xyz )

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2, c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            call fun ( 1, node_xyz, v )    

            node_num = node_num + 1
            result = result + area * v
            area_total = area_total + area

          end do
        end do
c
c  The other subtriangles have the opposite direction from the face.
c  Generate each in turn, by determining the barycentric coordinates
c  of the centroid (F1,F2,F3), from which we can also work out the barycentric
c  coordinates of the vertices of the subtriangle.
c
        do f3 = 2, 3 * factor - 4, 3
          do f2 = 2, 3 * factor - f3 - 2, 3

            f1 = 3 * factor - f3 - f2

            call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
     &        f2, f3, node_xyz )

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            call fun ( 1, node_xyz, v )  

            node_num = node_num + 1  
            result = result + area * v
            area_total = area_total + area

          end do
        end do

      end do

      return
      end
      subroutine sphere01_quad_icos1m ( factor, fun, node_num, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
c
c  Discussion:
c
c    This function estimates an integral over the surface of the unit sphere.
c
c    This function sets up an icosahedral grid, and subdivides each
c    edge of the icosahedron into FACTOR subedges.  These edges define a grid
c    within each triangular icosahedral face.  The midsides of these
c    triangles can be determined.  All of these calculations are done,
c    essentially, on the FLAT faces of the icosahedron.  Only then are
c    the triangle vertices and midsides projected to the sphere.  
c
c    The resulting grid of spherical triangles and projected midsides
c    is used to apply a midside quadrature rule over the surface of
c    the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 September 2010
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
c    Input, external :: FUN, evaluates the integrand, of the form:
c      subroutine fun ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Output, integer NODE_NUM, the number of evaluation points.
c
c    Output, double precision RESULT, the estimated integral.
c
      implicit none

      integer node_num

      integer a
      double precision a_xyz(3)
      double precision a2_xyz(3)
      double precision a3_xyz(3)
      double precision area
      double precision area_total
      integer b
      double precision b_xyz(3)
      double precision b2_xyz(3)
      double precision b3_xyz(3)
      integer c
      double precision c_xyz(3)
      double precision c2_xyz(3)
      double precision c3_xyz(3)
      integer edge
      integer edge_num
      integer edge_point(2,30)
      integer f
      integer f1
      integer f2
      integer f3
      integer face
      integer face_num
      integer face_order(20)
      integer face_point(3,20)
      integer face_order_max
      integer factor
      external fun
      integer i
      integer j
      integer node
      double precision node_norm
      double precision node_xyz(3)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision point_coord(3,12)
      integer point_num
      double precision result
      double precision va
      double precision vb
      double precision vc
c
c  Size the icosahedron.
c
      call icos_size ( point_num, edge_num, face_num, face_order_max )
c
c  Set the icosahedron.
c
      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Initialize the integral data.
c
      result = 0.0D+00
      node_num = 0
      area_total = 0.0D+00
c
c  Consider each face.
c
      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do i = 1, 3
          a_xyz(i) = point_coord(i,a)
          b_xyz(i) = point_coord(i,b)
          c_xyz(i) = point_coord(i,c)
        end do
c
c  Deal with subtriangles that have same orientation as face.
c
        do f1 = 0, factor - 1
          do f2 = 0, factor - f1 - 1
            f3 = factor - f1 - f2

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1, 2 * f3 - 2, 
     &        a3_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 + 1, 2 * f3 - 1, 
     &        b3_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,     2 * f3 - 1, 
     &        c3_xyz )

            node_num = node_num + 3
            call fun ( 1, a3_xyz, va )
            call fun ( 1, b3_xyz, vb )   
            call fun ( 1, c3_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
c
c  Deal with subtriangles that have opposite orientation as face.
c
        do f3 = 0, factor - 2
          do f2 = 1, factor - f3 - 1
            f1 = factor - f2 - f3

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1, 2 * f3 + 2, 
     &        a3_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 - 1, 2 * f3 + 1, 
     &        b3_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,     2 * f3 + 1, 
     &        c3_xyz )

            node_num = node_num + 3
            call fun ( 1, a3_xyz, va )
            call fun ( 1, b3_xyz, vb )   
            call fun ( 1, c3_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
      end do

      return
      end
      subroutine sphere01_quad_icos1v ( factor, fun, node_num, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
c
c  Discussion:
c
c    This function estimates an integral over the surface of the unit sphere.
c
c    This function sets up an icosahedral grid, and subdivides each
c    edge of the icosahedron into FACTOR subedges.  These edges define a grid
c    within each triangular icosahedral face.  The vertices of these
c    triangles can be determined.  All of these calculations are done,
c    essentially, on the FLAT faces of the icosahedron.  Only then are
c    the triangle vertices projected to the sphere.  
c
c    The resulting grid of spherical triangles is used to apply a vertex
c    quadrature rule over the surface of the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
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
c    Input, external :: FUN, evaluates the integrand, of the form:
c      subroutine fun ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Output, integer NODE_NUM, the number of evaluation points.
c
c    Output, double precision RESULT, the estimated integral.
c
      implicit none

      integer node_num

      integer a
      double precision a_xyz(3)
      double precision a2_xyz(3)
      double precision area
      double precision area_total
      integer b
      double precision b_xyz(3)
      double precision b2_xyz(3)
      integer c
      double precision c_xyz(3)
      double precision c2_xyz(3)
      integer edge
      integer edge_num
      integer edge_point(2,30)
      integer f
      integer f1
      integer f2
      integer f3
      integer face
      integer face_num
      integer face_order(20)
      integer face_point(3,20)
      integer face_order_max
      integer factor
      external fun
      integer i
      integer j
      integer node
      double precision node_norm
      double precision node_xyz(3)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision point_coord(3,12)
      integer point_num
      double precision result
      double precision va
      double precision vb
      double precision vc
c
c  Size the icosahedron.
c
      call icos_size ( point_num, edge_num, face_num, face_order_max )
c
c  Set the icosahedron.
c
      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Initialize the integral data.
c
      result = 0.0D+00
      node_num = 0
      area_total = 0.0D+00
c
c  Consider each face.
c
      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do i = 1, 3
          a_xyz(i) = point_coord(i,a)
          b_xyz(i) = point_coord(i,b)
          c_xyz(i) = point_coord(i,c)
        end do
c
c  Deal with subtriangles that have same orientation as face.
c
        do f1 = 0, factor - 1
          do f2 = 0, factor - f1 - 1
            f3 = factor - f1 - f2

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            node_num = node_num + 3
            call fun ( 1, a2_xyz, va )
            call fun ( 1, b2_xyz, vb )   
            call fun ( 1, c2_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
c
c  Deal with subtriangles that have opposite orientation as face.
c
        do f3 = 0, factor - 2
          do f2 = 1, factor - f3 - 1
            f1 = factor - f2 - f3

            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
            call sphere01_triangle_project ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            node_num = node_num + 3
            call fun ( 1, a2_xyz, va )
            call fun ( 1, b2_xyz, vb )   
            call fun ( 1, c2_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
      end do

      return
      end
      subroutine sphere01_quad_icos2v ( factor, fun, node_num, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
c
c  Discussion:
c
c    This function estimates an integral over the surface of the unit sphere.
c
c    This function sets up an icosahedral grid, and subdivides each
c    edge of the icosahedron into FACTOR subedges.  These edges define a grid
c    within each triangular icosahedral face.  The vertices of these
c    triangles can be determined.  All of these calculations are done,
c    essentially, on the FLAT faces of the icosahedron.  Only then are
c    the triangle vertices projected to the sphere.  
c
c    The resulting grid of spherical triangles is used to apply a vertex
c    quadrature rule over the surface of the unit sphere.
c
c    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
c    sophisticated scheme to map points from the planar triangle to the surface
c    of the unit sphere.  Very little improvement to the estimated integral
c    was observed, so development of this scheme has been set aside for now.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 September 2010
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
c    Input, external :: FUN, evaluates the integrand, of the form:
c      subroutine fun ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Output, integer NODE_NUM, the number of evaluation points.
c
c    Output, double precision RESULT, the estimated integral.
c
      implicit none

      integer node_num

      integer a
      double precision a_xyz(3)
      double precision a2_xyz(3)
      double precision area
      double precision area_total
      integer b
      double precision b_xyz(3)
      double precision b2_xyz(3)
      integer c
      double precision c_xyz(3)
      double precision c2_xyz(3)
      integer edge
      integer edge_num
      integer edge_point(2,30)
      integer f
      integer f1
      integer f2
      integer f3
      integer face
      integer face_num
      integer face_order(20)
      integer face_point(3,20)
      integer face_order_max
      integer factor
      external fun
      integer i
      integer j
      integer node
      double precision node_norm
      double precision node_xyz(3)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision point_coord(3,12)
      integer point_num
      double precision result
      double precision va
      double precision vb
      double precision vc
c
c  Size the icosahedron.
c
      call icos_size ( point_num, edge_num, face_num, face_order_max )
c
c  Set the icosahedron.
c
      call icos_shape ( point_num, edge_num, face_num, face_order_max, 
     &  point_coord, edge_point, face_order, face_point )
c
c  Initialize the integral data.
c
      result = 0.0D+00
      node_num = 0
      area_total = 0.0D+00
c
c  Consider each face.
c
      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do i = 1, 3
          a_xyz(i) = point_coord(i,a)
          b_xyz(i) = point_coord(i,b)
          c_xyz(i) = point_coord(i,c)
        end do
c
c  Deal with subtriangles that have same orientation as face.
c
        do f1 = 0, factor - 1
          do f2 = 0, factor - f1 - 1
            f3 = factor - f1 - f2

            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            node_num = node_num + 3
            call fun ( 1, a2_xyz, va )
            call fun ( 1, b2_xyz, vb )   
            call fun ( 1, c2_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
c
c  Deal with subtriangles that have opposite orientation as face.
c
        do f3 = 0, factor - 2
          do f2 = 1, factor - f3 - 1
            f1 = factor - f2 - f3

            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
            call sphere01_triangle_project2 ( 
     &        a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

            call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, 
     &        c2_xyz, area )

            node_num = node_num + 3
            call fun ( 1, a2_xyz, va )
            call fun ( 1, b2_xyz, vb )   
            call fun ( 1, c2_xyz, vc )   
            result = result + area * ( va + vb + vc ) / 3.0D+00
            area_total = area_total + area

          end do
        end do
      end do

      return
      end
      subroutine sphere01_quad_llc ( f, h, n, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
c
c  Discussion:
c
c    The sphere is broken up into spherical triangles, whose sides
c    do not exceed the length H.  Then a centroid rule is used on
c    each spherical triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external :: F, evaluates the integrand, of the form:
c      subroutine f ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Input, double precision H, the maximum length of a side of the spherical
c    quadrilaterals.
c
c    Output, integer N, the number of points used.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision area
      external             f
      double precision h
      integer i
      integer j
      integer n
      double precision phi
      integer phi_num
      double precision phi1
      double precision phi2
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision result
      double precision sector_area
      double precision sphere_area
      double precision theta
      integer theta_num
      double precision theta1
      double precision theta2
      double precision v(1)
      double precision x(3)
      double precision x1(3)
      double precision x11(3)
      double precision x12(3)
      double precision x2(3)
      double precision x21(3)
      double precision x22(3)
c
c  Choose PHI and THETA counts that make short sides.
c
      phi_num = int ( pi / h )

      if ( h * dble ( phi_num ) .lt. pi ) then
        phi_num = phi_num + 1
      end if

      theta_num = int ( 2.0D+00 * pi / h )

      if ( h * dble ( theta_num ) .lt. pi ) then
        theta_num = theta_num + 1
      end if

      n = 0
      result = 0.0D+00
c
c  Only one THETA (and hence, only one PHI.)
c
      if ( theta_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi

        theta = 0.0D+00
        phi = pi / 2.0D+00
        call tp_to_xyz ( theta, phi, x )

        call f ( 1, x, v )
        n = n + 1
        result = sphere_area * v(1)
c
c  Several THETA's, one PHI.
c
      else if ( phi_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi
        sector_area = sphere_area / dble ( theta_num )

        result = 0.0D+00

        do j = 1, theta_num

          theta = dble ( j - 1 ) * 2 * pi / dble ( theta_num )
          phi = pi / 2.0D+00
          call tp_to_xyz ( theta, phi, x )
          call f ( 1, x, v )
          n = n + 1
          result = result + sector_area * v(1)

        end do
c
c  At least two PHI's.
c
      else

        result = 0.0D+00
c
c  Picture in top row, with V1 = north pole:
c
c        V1
c       /  \
c      /    \
c    V12----V22
c
        phi1 = 0.0D+00
        phi2 = pi / dble ( phi_num )

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x1 )
          call tp_to_xyz ( theta1, phi2, x12 )
          call tp_to_xyz ( theta2, phi2, x22 )

          call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )
          call sphere01_triangle_vertices_to_centroid ( x1, x12, x22, 
     &      x )
          call f ( 1, x, v )
          n = n + 1
          result = result + area * v(1)

        end do
c
c  Picture in all intermediate rows:
c
c    V11--V21
c     | \  |
c     |  \ |
c    V12--V22
c
        do i = 2, phi_num-1

          phi1 = dble ( i - 1 ) * pi / dble ( phi_num )
          phi2 = dble ( i ) * pi / dble ( phi_num )

          do j = 1, theta_num

            theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
            theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

            call tp_to_xyz ( theta1, phi1, x11 )
            call tp_to_xyz ( theta2, phi1, x21 )
            call tp_to_xyz ( theta1, phi2, x12 )
            call tp_to_xyz ( theta2, phi2, x22 )

            call sphere01_triangle_vertices_to_area ( x11, x12, x22, 
     &        area )
            call sphere01_triangle_vertices_to_centroid ( x11, x12, 
     &        x22, x )
            call f ( 1, x, v )
            n = n + 1
            result = result + area * v(1)

            call sphere01_triangle_vertices_to_area ( x22, x21, x11, 
     &        area )
            call sphere01_triangle_vertices_to_centroid ( x22, x21, 
     &        x11, x )
            call f ( 1, x, v )
            n = n + 1
            result = result + area * v(1)

          end do

        end do
c
c  Picture in last row, with V2 = south pole:
c
c    V11----V21
c      \    /
c       \  /
c        V2
c
        phi1 = dble ( phi_num - 1 ) * pi / dble ( phi_num )
        phi2 = pi

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x11 )
          call tp_to_xyz ( theta2, phi1, x21 )
          call tp_to_xyz ( theta2, phi2, x2 )

          call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )
          call sphere01_triangle_vertices_to_centroid ( x11, x2, x21, 
     &      x )
          call f ( 1, x, v )
          n = n + 1
          result = result + area * v(1)

        end do

      end if

      return
      end
      subroutine sphere01_quad_llm ( f, h, n, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
c
c  Discussion:
c
c    The sphere is broken up into spherical triangles, whose sides
c    do not exceed the length H.  Then the function is evaluated
c    at the midsides, and the average is multiplied by the area.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 June 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external :: F, evaluates the integrand, of the form:
c      subroutine f ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Input, double precision H, the maximum length of a side of the spherical
c    quadrilaterals.
c
c    Output, integer N, the number of points used.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision area
      external             f
      double precision h
      integer i
      integer j
      double precision m1(3)
      double precision m2(3)
      double precision m3(3)
      integer n
      double precision phi
      integer phi_num
      double precision phi1
      double precision phi2
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision result
      double precision sector_area
      double precision sphere_area
      double precision theta
      integer theta_num
      double precision theta1
      double precision theta2
      double precision v(1)
      double precision x(3)
      double precision x1(3)
      double precision x11(3)
      double precision x12(3)
      double precision x2(3)
      double precision x21(3)
      double precision x22(3)
c
c  Choose PHI and THETA counts that make short sides.
c
      phi_num = int ( pi / h )

      if ( h * dble ( phi_num ) .lt. pi ) then
        phi_num = phi_num + 1
      end if

      theta_num = int ( 2.0D+00 * pi / h )

      if ( h * dble ( theta_num ) .lt. pi ) then
        theta_num = theta_num + 1
      end if

      n = 0
      result = 0.0D+00
c
c  Only one THETA (and hence, only one PHI.)
c
      if ( theta_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi

        theta = 0.0D+00
        phi = pi / 2.0D+00
        call tp_to_xyz ( theta, phi, x )
        call f ( 1, x, v )
        n = n + 1
        result = sphere_area * v(1)
c
c  Several THETA's, one PHI.
c
      else if ( phi_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi
        sector_area = sphere_area / dble ( theta_num )

        result = 0.0D+00

        do j = 1, theta_num

          theta = dble ( ( j - 1 ) * 2 ) * pi / dble ( theta_num )
          phi = pi / 2.0D+00
          call tp_to_xyz ( theta, phi, x )
          call f ( 1, x, v )
          n = n + 1
          result = result + sector_area * v(1)

        end do
c
c  At least two PHI's.
c
      else

        result = 0.0D+00
c
c  Picture:
c
c        V1
c       /  \
c      /    \
c    V12----V22
c
        phi1 = 0.0D+00
        phi2 = pi / dble ( phi_num )

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x1 )
          call tp_to_xyz ( theta1, phi2, x12 )
          call tp_to_xyz ( theta2, phi2, x22 )

          call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

          call sphere01_triangle_vertices_to_midpoints ( x1, x12, x22, 
     &      m1, m2, m3 )

          call f ( 1, m1, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, m2, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, m3, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00

        end do
c
c  Picture:
c
c    V11--V21
c     | \  |
c     |  \ |
c    V12--V22
c
        do i = 2, phi_num-1

          phi1 = dble ( i - 1 ) * pi / dble ( phi_num )
          phi2 = dble ( i ) * pi / dble ( phi_num )

          do j = 1, theta_num

            theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
            theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

            call tp_to_xyz ( theta1, phi1, x11 )
            call tp_to_xyz ( theta2, phi1, x21 )
            call tp_to_xyz ( theta1, phi2, x12 )
            call tp_to_xyz ( theta2, phi2, x22 )

            call sphere01_triangle_vertices_to_area ( x11, x12, x22, 
     &        area )

            call sphere01_triangle_vertices_to_midpoints ( x11, x12, 
     &        x22, m1, m2, m3 )

            call f ( 1, m1, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, m2, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, m3, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00

            call sphere01_triangle_vertices_to_area ( x22, x21, x11, 
     &        area )

            call sphere01_triangle_vertices_to_midpoints ( x22, x21, 
     &        x11, m1, m2, m3 )

            call f ( 1, m1, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, m2, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, m3, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00

          end do

        end do
c
c  Picture:
c
c    V11----V21
c      \    /
c       \  /
c        V2
c
        phi1 = dble ( phi_num - 1 ) * pi / dble ( phi_num )
        phi2 = pi

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x11 )
          call tp_to_xyz ( theta2, phi1, x21 )
          call tp_to_xyz ( theta2, phi2, x2 )

          call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

          call sphere01_triangle_vertices_to_midpoints ( x11, x2, x21, 
     &      m1, m2, m3 )

          call f ( 1, m1, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, m2, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, m3, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00

        end do

      end if

      return
      end
      subroutine sphere01_quad_llv ( f, h, n, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
c
c  Discussion:
c
c    The sphere is broken up into spherical triangles, whose sides
c    do not exceed the length H.  Then the function is evaluated
c    at the vertices, and the average is multiplied by the area.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 June 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external :: F, evaluates the integrand, of the form:
c      subroutine f ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Input, double precision H, the maximum length of a side of the spherical
c    quadrilaterals.
c
c    Output, integer N, the number of points used.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision area
      external             f
      double precision h
      integer i
      integer j
      integer n
      double precision phi
      integer phi_num
      double precision phi1
      double precision phi2
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision result
      double precision sector_area
      double precision sphere_area
      double precision theta
      integer theta_num
      double precision theta1
      double precision theta2
      double precision v(1)
      double precision x(3)
      double precision x1(3)
      double precision x11(3)
      double precision x12(3)
      double precision x2(3)
      double precision x21(3)
      double precision x22(3)
c
c  Choose PHI and THETA counts that make short sides.
c
      phi_num = int ( pi / h )

      if ( h * dble ( phi_num ) .lt. pi ) then
        phi_num = phi_num + 1
      end if

      theta_num = int ( 2.0D+00 * pi / h )

      if ( h * dble ( theta_num ) .lt. pi ) then
        theta_num = theta_num + 1
      end if

      n = 0
      result = 0.0D+00
c
c  Only one THETA (and hence, only one PHI.)
c
      if ( theta_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi

        theta = 0.0D+00
        phi = pi / 2.0D+00
        call tp_to_xyz ( theta, phi, x )
        call f ( 1, x, v )
        result = sphere_area * v(1)
c
c  Several THETA's, one PHI.
c
      else if ( phi_num .eq. 1 ) then

        sphere_area = 4.0D+00 * pi
        sector_area = sphere_area / dble ( theta_num )

        result = 0.0D+00

        do j = 1, theta_num

          theta = dble ( ( j - 1 ) * 2 ) * pi / dble ( theta_num )
          phi = pi / 2.0D+00
          call tp_to_xyz ( theta, phi, x )
          call f ( 1, x, v )
          n = n + 1
          result = result + sector_area * v(1)

        end do
c
c  At least two PHI's.
c
      else

        result = 0.0D+00
c
c  Picture:
c
c        V1
c       /  \
c      /    \
c    V12----V22
c
        phi1 = 0.0D+00
        phi2 = pi / dble ( phi_num )

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x1 )
          call tp_to_xyz ( theta1, phi2, x12 )
          call tp_to_xyz ( theta2, phi2, x22 )

          call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

          call f ( 1, x1, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, x12, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, x22, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00

        end do
c
c  Picture:
c
c    V11--V21
c     | \  |
c     |  \ |
c    V12--V22
c
        do i = 2, phi_num-1

          phi1 = dble ( i - 1 ) * pi / dble ( phi_num )
          phi2 = dble ( i ) * pi / dble ( phi_num )

          do j = 1, theta_num

            theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
            theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

            call tp_to_xyz ( theta1, phi1, x11 )
            call tp_to_xyz ( theta2, phi1, x21 )
            call tp_to_xyz ( theta1, phi2, x12 )
            call tp_to_xyz ( theta2, phi2, x22 )

            call sphere01_triangle_vertices_to_area ( x11, x12, x22, 
     &        area )

            call f ( 1, x11, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, x12, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, x22, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00

            call sphere01_triangle_vertices_to_area ( x22, x21, x11, 
     &        area )

            call f ( 1, x22, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, x21, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00
            call f ( 1, x11, v )
            n = n + 1
            result = result + area * v(1) / 3.0D+00

          end do

        end do
c
c  Picture:
c
c    V11----V21
c      \    /
c       \  /
c        V2
c
        phi1 = dble ( phi_num - 1 ) * pi / dble ( phi_num )
        phi2 = pi

        do j = 1, theta_num

          theta1 = dble ( j - 1 ) * 2.0D+00 * pi / dble ( theta_num )
          theta2 = dble ( j ) * 2.0D+00 * pi / dble ( theta_num )

          call tp_to_xyz ( theta1, phi1, x11 )
          call tp_to_xyz ( theta2, phi1, x21 )
          call tp_to_xyz ( theta2, phi2, x2 )

          call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

          call f ( 1, x11, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, x2, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00
          call f ( 1, x21, v )
          n = n + 1
          result = result + area * v(1) / 3.0D+00

        end do

      end if

      return
      end
      subroutine sphere01_quad_mc ( f, h, seed, n, result )

c*********************************************************************72
c
cc SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
c
c  Discussion:
c
c    A number of points N are chosen at random on the sphere, with N
c    being determined so that, if the points were laid out on a regular
c    grid, the average spacing would be no more than H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external :: F, evaluates the integrand, of the form:
c      subroutine f ( n, x, v )
c      integer n
c      double precision v(n)
c      double precision x(3,n)
c
c    Input, double precision H, the maximum length of a side of the spherical
c    quadrilaterals.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Input, integer N, the number of points to use.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      external f
      double precision h
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8vec_sum
      double precision result
      integer seed
      double precision sphere_area
      double precision v(n)
      double precision x(3,n)
c
c  The sphere's area is 4 * PI.
c  Choose N so that we divide this area into N subareas of PI * H * H.
c
      sphere_area = 4.0D+00 * pi

      call sphere01_sample_3d ( n, seed, x )

      call f ( n, x, v )

      result = sphere_area * r8vec_sum ( n, v ) / dble ( n )

      return
      end
      subroutine sphere01_quad_mc_size ( h, n )

c*********************************************************************72
c
cc SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
c
c  Discussion:
c
c    A number of points N are chosen at random on the sphere, with N
c    being determined so that, if the points were laid out on a regular
c    grid, the average spacing would be no more than H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision H, the maximum length of a side of the spherical
c    quadrilaterals.
c
c    Output, integer N, the number of points to use.
c
      implicit none

      double precision h
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sphere_area
c
c  The sphere's area is 4 * PI.
c  Choose N so that we divide this area into N subareas of PI * H * H.
c
      sphere_area = 4.0D+00 * pi

      n = int ( sphere_area / h**2 )
      n = max ( n, 1 )

      return
      end
      subroutine sphere01_sample_3d ( n, seed, x )

c*********************************************************************72
c
cc SPHERE01_SAMPLE_3D picks random points on a sphere in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of samples.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision X(3,N), the sample points.
c
      implicit none

      integer n

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

        phi = acos ( vdot )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the Z vector.
c
        theta = r8_uniform_01 ( seed )
        theta = 2.0D+00 * pi * theta

        x(1,j) = cos ( theta ) * sin ( phi )
        x(2,j) = sin ( theta ) * sin ( phi )
        x(3,j) =                 cos ( phi )

      end do

      return
      end
      subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a spherical triangle on unit sphere.
c
c  Discussion:
c
c    A unit sphere centered at 0 in 3D satisfies the equation:
c
c      X*X + Y*Y + Z*Z = 1
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.
c
c    The area formula is known as Girard's formula.
c
c    The area of a spherical triangle on the unit sphere is:
c
c      AREA = ( A + B + C - PI )
c
c    where A, B and C are the (surface) angles of the triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the angles of the triangle.
c
c    Output, double precision AREA, the area of the spherical triangle.
c
      implicit none

      double precision a
      double precision area
      double precision b
      double precision c
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
c
c  Apply Girard's formula.
c
      area = a + b + c - pi

      return
      end
      subroutine sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
     &  f2, f3, node_xyz )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_PROJECT projects from a plane triangle to a spherical triangle.
c
c  Discussion:
c
c    We assume that points A, B and C lie on the unit sphere, and they
c    thus define a spherical triangle.
c
c    They also, of course, define a planar triangle.
c
c    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
c    planar triangle.
c
c    This function determines the coordinates of the point in the planar
c    triangle identified by the barycentric coordinates, and returns the
c    coordinates of the projection of that point onto the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
c    of the points A, B, and C.
c
c    Input, integer F1, F2, F3, the barycentric coordinates
c    of a point in the triangle ABC.  Normally, these coordinates would
c    be real numbers, and would sum to 1.  For convenience, we allow these
c    to be integers which must be divided by F1+F2+F3.
c
c    Output, double precision NODE_XYZ(3), the coordinates of the 
c    point on the unit sphere which is the projection of the point on the plane
c    whose barycentric coordinates with respect to A, B, and C is
c    (F1,F2,F3)/(F1+F2+F3).
c
      implicit none

      double precision a_xyz(3)
      double precision b_xyz(3)
      double precision c_xyz(3)
      integer f1
      integer f2
      integer f3
      integer i
      double precision node_norm
      double precision node_xyz(3)
      double precision r8vec_norm

      do i = 1, 3
        node_xyz(i) = 
     &    ( dble ( f1           ) * a_xyz(i)   
     &    + dble (      f2      ) * b_xyz(i)   
     &    + dble (           f3 ) * c_xyz(i) ) 
     &    / dble ( f1 + f2 + f3 )
      end do

      node_norm = r8vec_norm ( 3, node_xyz )

      do i = 1, 3
        node_xyz(i) = node_xyz(i) / node_norm
      end do

      return
      end
      subroutine sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, 
     &  f2, f3, node_xyz )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_PROJECT2 projects from a plane triangle to a spherical triangle.
c
c  Discussion:
c
c    We assume that points A, B and C lie on the unit sphere, and they
c    thus define a spherical triangle.
c
c    They also, of course, define a planar triangle.
c
c    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
c    planar triangle.
c
c    This function determines the coordinates of the point in the planar
c    triangle identified by the barycentric coordinates, and returns the
c    coordinates of the projection of that point onto the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
c    of the points A, B, and C.
c
c    Input, integer F1, F2, F3, the barycentric coordinates
c    of a point in the triangle ABC.  Normally, these coordinates would
c    be real numbers, and would sum to 1.  For convenience, we allow these
c    to be integers which must be divided by F1+F2+F3.
c
c    Output, double precision NODE_XYZ(3), the coordinates of the 
c    point on the unit sphere which is the projection of the point on the 
c    plane whose barycentric coordinates with respect to A, B, and C is
c    (F1,F2,F3)/(F1+F2+F3).
c
      implicit none

      double precision a_xyz(3)
      double precision ab(3)
      double precision ac(3)
      double precision acn(3)
      double precision acp(3)
      double precision angle
      double precision b_xyz(3)
      double precision bn(3)
      double precision bp(3)
      double precision c_xyz(3)
      double precision cn(3)
      double precision cp(3)
      integer f1
      integer f2
      integer f3
      integer i
      double precision node_xyz(3)
      double precision norm
      double precision r8vec_norm
      double precision theta_ab
      double precision theta_ac
      double precision theta_bc
c
c  This check avoids 0/0 calculations later.
c
      if ( f2 .eq. 0 .and. f3 .eq. 0 ) then
        do i = 1, 3
          node_xyz(i) = a_xyz(i)
        end do
        return
      else if ( f1 .eq. 0 .and. f3 .eq. 0 ) then
        do i = 1, 3
          node_xyz(i) = b_xyz(i)
        end do
        return
      else if ( f1 .eq. 0 .and. f2 .eq. 0 ) then
        do i = 1, 3
          node_xyz(i) = c_xyz(i)
        end do
        return
      end if
c
c  Determine the angular distances (A,B) and (A,C).
c
      call sphere01_distance_xyz ( a_xyz, b_xyz, theta_ab )

      call sphere01_distance_xyz ( a_xyz, c_xyz, theta_ac )
c
c  Polarize B = BP + BN
c  Normalize BN, 
c  Same for C.
c
      call r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp )
      norm = r8vec_norm ( 3, bn )
      do i = 1, 3
        bn(i) = bn(i) / norm
      end do

      call r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp )
      norm = r8vec_norm ( 3, cn )
      do i = 1, 3
        cn(i) = cn(i) / norm
      end do
c
c  Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
c  and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
c
      angle = ( dble ( f2 + f3 ) * theta_ab ) / dble ( f1 + f2 + f3 )
      do i = 1, 3
        ab(i) = cos ( angle ) * a_xyz(i) + sin ( angle ) * bn(i)
      end do

      angle = ( dble ( f2 + f3 ) * theta_ac ) / dble ( f1 + f2 + f3 )
      do i = 1, 3
        ac(i) = cos ( angle ) * a_xyz(i) + sin ( angle ) * cn(i)
      end do
c
c  Determine the angular distance between AB and AC.
c
      call sphere01_distance_xyz ( ab, ac, theta_bc )
c
c  Polarize AC = ACP + ACN, normalize ACN.
c
      call r8vec_polarize ( 3, ac, ab, acn, acp )
      norm = r8vec_norm ( 3, acn )
      do i = 1, 3
        acn(i) = acn(i) / norm
      end do
c
c  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
c
      angle = ( dble ( f3 ) * theta_bc ) / dble ( f2 + f3 )

      do i = 1, 3
        node_xyz(i) = cos ( angle ) * ab(i) + sin ( angle ) * acn(i)
      end do

      return
      end
      subroutine sphere01_triangle_sample ( n, v1, v2, v3, seed, x )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
c
c  Discussion:
c
c    The sphere has center 0 and radius 1.
c
c    A spherical triangle on the surface of the unit sphere contains those 
c    points with radius R = 1, bounded by the vertices V1, V2, V3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Arvo,
c    Stratified sampling of spherical triangles,
c    Computer Graphics Proceedings, Annual Conference Series, 
c    ACM SIGGRAPH '95, pages 437-438, 1995.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision V1(3), V2(3), V3(3), the XYZ coordinates of
c    the vertices of the spherical triangle.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision X(3,N), the XYZ coordinates of the 
c    sample points.
c
      implicit none

      integer n

      double precision a
      double precision alpha
      double precision area
      double precision area_hat
      double precision b
      double precision beta
      double precision c
      double precision gamma
      integer i
      integer j
      double precision norm
      double precision q
      double precision r8_uniform_01
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision s
      integer seed
      double precision t
      double precision u
      double precision v
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
      double precision v31(3)
      double precision v4(3)
      double precision v42(3)
      double precision w
      double precision x(3,n)
      double precision xsi1
      double precision xsi2
      double precision z

      call sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c )

      call sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, 
     &  gamma )

      call sphere01_triangle_angles_to_area ( alpha, beta, gamma, area )

      do j = 1, n
c
c  Select the new area.
c
        xsi1 = r8_uniform_01 ( seed )
        area_hat = xsi1 * area
c
c  Compute the sine and cosine of the angle phi.
c
        s = sin ( area_hat - alpha )
        t = cos ( area_hat - alpha )
c
c  Compute the pair that determines beta_hat.
c
        u = t - cos ( alpha )
        v = s + sin ( alpha ) * cos ( c )
c
c  Q is the cosine of the new edge length b_hat.
c
        q = ( ( v * t - u * s ) * cos ( alpha ) - v ) 
     &    / ( ( v * s + u * t ) * sin ( alpha ) )
c
c  Occasionally we get a Q value out of bounds.
c
        q = max ( q, - 1.0D+00 )
        q = min ( q, + 1.0D+00 )
c
c  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
c
        w = r8vec_dot_product ( 3, v3, v1 )
        do i = 1, 3
          v31(i) = v3(i) - w * v1(i)
        end do
        norm = r8vec_norm ( 3, v31 )
        do i = 1, 3
          v31(i) = v31(i) / norm
        end do
c
c  V4 is the third vertex of the subtriangle V1, V2, V4.
c
        do i = 1, 3
          v4(i) = q * v1(i) + sqrt ( 1.0D+00 - q * q ) * v31(i)
        end do
c
c  Select cos theta, which will sample along the edge from V2 to V4.
c
        xsi2 = r8_uniform_01 ( seed )
        z = 1.0D+00 - xsi2 
     &    * ( 1.0D+00 - r8vec_dot_product ( 3, v4, v2 ) )
c
c  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
c
        w = r8vec_dot_product ( 3, v4, v2 )
        do i = 1, 3
          v42(i) = v4(i) - w * v2(i)
        end do
        norm = r8vec_norm ( 3, v42 )
        do i = 1, 3
          v42(i) = v42(i) / norm
        end do
c
c  Construct the point.
c
        do i = 1, 3
          x(i,j) = z * v2(i) + sqrt ( 1.0D+00 - z * z ) * v42(i)
        end do

      end do

      return
      end
      subroutine sphere01_triangle_sides_to_angles ( as, bs, cs, 
     &  a, b, c )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_SIDES_TO_ANGLES; spherical triangle angles on the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision AS, BS, CS, the (geodesic) length of the 
c    sides of the triangle.
c
c    Output, double precision A, B, C, the spherical angles of the triangle.
c    Angle A is opposite the side of length AS, and so on.
c
      implicit none

      double precision a
      double precision as
      double precision asu
      double precision b
      double precision bs
      double precision bsu
      double precision c
      double precision cs
      double precision csu
      double precision ssu
      double precision tan_a2
      double precision tan_b2
      double precision tan_c2

      asu = as
      bsu = bs
      csu = cs
      ssu = ( asu + bsu + csu ) / 2.0D+00

      tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / 
     &                ( sin ( ssu ) * sin ( ssu - asu )     ) )

      a = 2.0D+00 * atan ( tan_a2 )

      tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / 
     &                ( sin ( ssu ) * sin ( ssu - bsu )     ) )

      b = 2.0D+00 * atan ( tan_b2 )

      tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / 
     &                ( sin ( ssu ) * sin ( ssu - csu )     ) )

      c = 2.0D+00 * atan ( tan_c2 )

      return
      end
      subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of a spherical triangle on the unit sphere.
c
c  Discussion:
c
c    A unit sphere centered at 0 in 3D satisfies the equation:
c
c      X * X + Y * Y + Z * Z = 1
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.
c
c    The area formula is known as Girard's formula.
c
c    The area of a spherical triangle on the unit sphere is:
c
c      AREA = ( A + B + C - PI ) * 1
c
c    where A, B and C are the (surface) angles of the triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Output, double precision AREA, the area of the spherical triangle.
c
      implicit none

      double precision a
      double precision area
      double precision as
      double precision b
      double precision bs
      double precision c
      double precision cs
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
c
c  Compute the lengths of the sides of the spherical triangle.
c
      call sphere01_triangle_vertices_to_sides ( v1, v2, v3, 
     &  as, bs, cs )
c
c  Get the spherical angles.
c
      call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )
c
c  Get the area.
c
      call sphere01_triangle_angles_to_area ( a, b, c, area )

      return
      end
      subroutine sphere01_triangle_vertices_to_centroid ( v1, v2, v3, 
     &  vs )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: spherical triangle centroid on unit sphere.
c
c  Discussion:
c
c    A unit sphere centered at 0 in 3D satisfies the equation:
c
c      X*X + Y*Y + Z*Z = 1
c
c    A spherical triangle is specified by three points on the sphere.
c
c    The (true) centroid of a spherical triangle is the point
c
c      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
c
c    Note that the true centroid does NOT, in general, lie on the sphere.  
c
c    The "flat" centroid VF is the centroid of the planar triangle defined by
c    the vertices of the spherical triangle.
c
c    The "spherical" centroid VS of a spherical triangle is computed by
c    the intersection of the geodesic bisectors of the triangle angles.
c    The spherical centroid lies on the sphere.
c
c    VF, VT and VS lie on a line through the center of the sphere.  We can
c    easily calculate VF by averaging the vertices, and from this determine
c    VS by normalizing.
c
c    Of course, we still will not have actually computed VT, which lies
c    somewhere between VF and VSc
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Output, double precision VS(3), the coordinates of the "spherical
c    centroid" of the spherical triangle.
c
      implicit none

      integer i
      double precision norm
      double precision r8vec_norm
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
      double precision vs(3)

      do i = 1, 3
        vs(i) = ( v1(i) + v2(i) + v3(i) ) / 3.0D+00
      end do

      norm = r8vec_norm ( 3, vs )

      do i = 1, 3
        vs(i) = vs(i) / norm
      end do

      return
      end
      subroutine sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, 
     &  m1, m2, m3 )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS gets the midsides of a spherical triangle.
c
c  Discussion:
c
c    The points are assumed to lie on the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 June 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Output, double precision M1(3), M2(3), M3(3), the coordinates of 
c    the midpoints of the sides of the spherical triangle.
c
      implicit none

      integer i
      double precision m1(3)
      double precision m2(3)
      double precision m3(3)
      double precision norm
      double precision r8vec_norm
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)

      do i = 1, 3
        m1(i) = ( v1(i) + v2(i) ) / 2.0D+00
      end do
      norm = r8vec_norm ( 3, m1 )
      do i = 1, 3
        m1(i) = m1(i) / norm
      end do

      do i = 1, 3
        m2(i) = ( v2(i) + v3(i) ) / 2.0D+00
      end do
      norm = r8vec_norm ( 3, m2 )
      do i = 1, 3
        m2(i) = m2(i) / norm
      end do

      do i = 1, 3
        m3(i) = ( v3(i) + v1(i) ) / 2.0D+00
      end do
      norm = r8vec_norm ( 3, m3 )
      do i = 1, 3
        m3(i) = m3(i) / norm
      end do

      return
      end
      subroutine sphere01_triangle_vertices_to_sides ( v1, v2, v3, 
     &  as, bs, cs )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_VERTICES_TO_SIDES: spherical triangle sides on unit sphere.
c
c  Discussion:
c
c    We can use the ACOS system call here, but the ARC_COSINE routine
c    will automatically take care of cases where the input argument is
c    (usually slightly) out of bounds.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the spherical
c    triangle.
c
c    Output, double precision AS, BS, CS, the (geodesic) length of the sides
c    of the triangle.
c
      implicit none

      double precision arc_cosine
      double precision as
      double precision bs
      double precision cs
      double precision r8vec_dot_product
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)

      as = arc_cosine ( r8vec_dot_product ( 3, v2, v3 ) )
      bs = arc_cosine ( r8vec_dot_product ( 3, v3, v1 ) )
      cs = arc_cosine ( r8vec_dot_product ( 3, v1, v2 ) )

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
      subroutine tp_to_xyz ( theta, phi, v )

c*********************************************************************72
c
cc TP_TO_XYZ converts unit spherical TP coordinates to XYZ coordinates.
c
c  Discussion:
c
c    The point is assume to lie on the unit sphere centered at the origin.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision THETA, PHI, the angular coordinates of a point
c    on the unit sphere.
c
c    Output, double precision V(3), the XYZ coordinates.
c
      implicit none

      double precision phi
      double precision theta
      double precision v(3)

      v(1) = cos ( theta ) * sin ( phi )
      v(2) = sin ( theta ) * sin ( phi )
      v(3) =                 cos ( phi )

      return
      end
