      subroutine piecewise_linear_product_integral ( a, b, f_num, f_x, 
     &  f_v, g_num, g_x, g_v, integral )

c*********************************************************************72
c
cc PIECEWISE_LINEAR_PRODUCT_INTEGRAL: piecewise linear product integral.
c
c  Discussion:
c
c    We are given two piecewise linear functions F(X) and G(X) and we wish
c    to compute the exact value of the integral
c
c      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
c
c    The functions F(X) and G(X) are defined as tables of coordinates X and
c    values V.  A piecewise linear function is evaluated at a point X by 
c    evaluating the interpolant to the data at the endpoints of the interval 
c    containing X.  
c
c    It must be the case that A <= B.
c
c    It must be the case that the node coordinates F_X(*) and G_X(*) are
c    given in ascending order.
c
c    It must be the case that:
c
c      F_X(1) <= A and B <= F_X(F_NUM)
c      G_X(1) <= A and B <= G_X(G_NUM)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer F_NUM, the number of nodes for F.
c
c    Input, double precision F_X(F_NUM), the node coordinates for F.
c
c    Input, double precision F_V(F_NUM), the nodal values for F.
c
c    Input, integer G_NUM, the number of nodes for G.
c
c    Input, double precision G_X(G_NUM), the node coordinates for G.
c
c    Input, double precision G_V(G_NUM), the nodal values for G.
c
c    Output, double precision INTEGRAL, the integral of F(X) * G(X)
c    from A to B.
c
      implicit none

      integer f_num
      integer g_num

      double precision a
      double precision b
      double precision bit
      integer f_left
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision f0
      double precision f1
      double precision fl
      double precision fr
      integer g_left
      double precision g_v(g_num)
      double precision g_x(g_num)
      double precision g0
      double precision g1
      double precision gl
      double precision gr
      double precision h0
      double precision h1
      double precision h2
      integer i
      double precision integral
      double precision r8_epsilon
      double precision xl
      double precision xr
      double precision xr_max

      integral = 0.0D+00

      if ( f_x(f_num) .le. a .or. g_x(g_num) .le. a ) then
        return
      end if

      if ( f_num .lt. 2 .or. g_num .lt. 2 ) then
        return
      end if

      xr = a

      f_left = 1
      call r8vec_bracket3 ( f_num, f_x, xr, f_left )
      fr = f_v(f_left) + ( xr - f_x(f_left) ) 
     &  * ( f_v(f_left+1) - f_v(f_left) ) 
     &  / ( f_x(f_left+1) - f_x(f_left) )

      g_left = 1
      call r8vec_bracket3 ( g_num, g_x, xr, g_left )
      gr = g_v(g_left) + ( xr - g_x(g_left) ) 
     &  * ( g_v(g_left+1) - g_v(g_left) ) 
     &  / ( g_x(g_left+1) - g_x(g_left) )

      xr_max = b
      xr_max = min ( xr_max, f_x(f_num) )
      xr_max = min ( xr_max, g_x(g_num) )

10    continue

      if ( xr .lt. xr_max ) then
c
c  Shift right values to left.
c
        xl = xr
        fl = fr
        gl = gr
c
c  Determine the new right values.
c  The hard part is figuring out how to advance XR some, but not too much.
c
        xr = xr_max

        do i = 1, 2
          if ( f_left + i .le. f_num ) then
            if ( xl .lt. f_x(f_left+i) .and. 
     &           f_x(f_left+i) .lt. xr ) then
              xr = f_x(f_left+i)
              go to 20
            end if
          end if
        end do

20      continue

        do i = 1, 2
          if ( g_left + i .le. g_num ) then
            if ( xl .lt. g_x(g_left+i) .and. 
     &           g_x(g_left+i) .lt. xr ) then
              xr = g_x(g_left+i)
              go to 30
            end if
          end if
        end do

30      continue

        call r8vec_bracket3 ( f_num, f_x, xr, f_left )
        fr = f_v(f_left) + ( xr - f_x(f_left) ) 
     &    * ( f_v(f_left+1) - f_v(f_left) ) 
     &    / ( f_x(f_left+1) - f_x(f_left) )

        call r8vec_bracket3 ( g_num, g_x, xr, g_left )
        gr = g_v(g_left) + ( xr - g_x(g_left) ) 
     &    * ( g_v(g_left+1) - g_v(g_left) ) 
     &    / ( g_x(g_left+1) - g_x(g_left) )
c
c  Form the linear polynomials for F(X) and G(X) over [XL,XR],
c  then the product H(X), integrate H(X) and add to the running total.
c
        if ( r8_epsilon ( ) .le. abs ( xr - xl ) ) then

          f1 = fl - fr
          f0 = fr * xl - fl * xr

          g1 = gl - gr
          g0 = gr * xl - gl * xr

          h2 = f1 * g1
          h1 = f1 * g0 + f0 * g1
          h0 = f0 * g0

          h2 = h2 / 3.0D+00
          h1 = h1 / 2.0D+00

          bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr 
     &        - ( ( h2 * xl + h1 ) * xl + h0 ) * xl

          integral = integral + bit / ( xr - xl ) / ( xr - xl )

        end if

        go to 10

      end if

      return
      end
      subroutine piecewise_linear_product_quad ( a, b, f_num, f_x, f_v,
     &  g_num, g_x, g_v, quad_num, quad )

c*********************************************************************72
c
c! PIECEWISE_LINEAR_PRODUCT_QUAD: estimate piecewise linear product integral.
c
c  Discussion:
c
c    We are given two piecewise linear functions F(X) and G(X) and we wish
c    to estimate the value of the integral
c
c      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
c
c    The functions F(X) and G(X) are defined as tables of coordinates X and
c    values V.  A piecewise linear function is evaluated at a point X by 
c    evaluating the interpolant to the data at the endpoints of the interval 
c    containing X.  
c
c    It must be the case that A <= B.
c
c    It must be the case that the node coordinates F_X(*) and G_X(*) are
c    given in ascending order.
c
c    It must be the case that:
c
c      F_X(1) <= A and B <= F_X(F_NUM)
c      G_X(1) <= A and B <= G_X(G_NUM)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer F_NUM, the number of nodes for F.
c
c    Input, double precision F_X(F_NUM), the node coordinates for F.
c
c    Input, double precision F_V(F_NUM), the nodal values for F.
c
c    Input, integer G_NUM, the number of nodes for G.
c
c    Input, double precision G_X(G_NUM), the node coordinates for G.
c
c    Input, double precision G_V(G_NUM), the nodal values for G.
c
c    Input, integer QUAD_NUM, the number of quadrature points.
c
c    Output, double precision QUAD, an estimate for the integral of F(X) * G(X)
c    from A to B.
c
      implicit none

      integer f_num
      integer g_num

      double precision a
      double precision a2
      double precision b
      double precision b2
      integer f_left
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision fq
      integer g_left
      double precision g_v(g_num)
      double precision g_x(g_num)
      double precision gq
      integer i
      double precision quad
      integer quad_num
      double precision xq

      quad = 0.0D+00

      f_left = 1
      g_left = 1

      a2 = a
      a2 = max ( a2, f_x(1) )
      a2 = max ( a2, g_x(1) )

      b2 = b
      b2 = min ( b2, f_x(f_num) )
      b2 = min ( b2, g_x(g_num) )

      do i = 1, quad_num

        xq =  ( dble (                2 * i - 1 ) * b2 
     &        + dble ( 2 * quad_num - 2 * i + 1 ) * a2 )  
     &        / dble ( 2 * quad_num             )

        call r8vec_bracket3 ( f_num, f_x, xq, f_left )

        fq = f_v(f_left) + ( xq - f_x(f_left) ) 
     &    * ( f_v(f_left+1) - f_v(f_left) ) 
     &    / ( f_x(f_left+1) - f_x(f_left) )

        call r8vec_bracket3 ( g_num, g_x, xq, g_left )

        gq = g_v(g_left) + ( xq - g_x(g_left) ) 
     &    * ( g_v(g_left+1) - g_v(g_left) ) 
     &    / ( g_x(g_left+1) - g_x(g_left) )

        quad = quad + fq * gq

      end do

      quad = quad * ( b - a ) / dble ( quad_num )

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
      subroutine r8vec_bracket3 ( n, t, tval, left )

c*********************************************************************72
c
cc R8VEC_BRACKET3 finds the interval containing or nearest a given value.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The routine always returns the index LEFT of the sorted array
c    T with the property that either
c    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
c    *  T .lt. T(LEFT) = T(1), or
c    *  T > T(LEFT+1) = T(N).
c
c    The routine is useful for interpolation problems, where
c    the abscissa must be located within an interval of data
c    abscissas for interpolation, or the "nearest" interval
c    to the (extreme) abscissa must be found so that extrapolation
c    can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of the input array.
c
c    Input, double precision T(N), an array that has been sorted
c    into ascending order.
c
c    Input, double precision TVAL, a value to be bracketed by entries of T.
c
c    Input/output, integer LEFT.
c    On input, if 1 .le. LEFT .le. N-1, LEFT is taken as a suggestion for the
c    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
c    is searched first, followed by the appropriate interval to the left
c    or right.  After that, a binary search is used.
c    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
c    is the closest to TVAL; it either contains TVAL, or else TVAL
c    lies outside the interval [ T(1), T(N) ].
c
      implicit none

      integer n

      integer high
      integer left
      integer low
      integer mid
      double precision t(n)
      double precision tval
c
c  Check the input data.
c
      if ( n .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
        write ( *, '(a)' ) '  N must be at least 2.'
        stop
      end if
c
c  If LEFT is not between 1 and N-1, set it to the middle value.
c
      if ( left .lt. 1 .or. n - 1 .lt. left ) then
        left = ( n + 1 ) / 2
      end if
c
c  CASE 1: TVAL .lt. T(LEFT):
c  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
c
      if ( tval .lt. t(left) ) then

        if ( left .eq. 1 ) then
          return
        else if ( left .eq. 2 ) then
          left = 1
          return
        else if ( t(left-1) .le. tval ) then
          left = left - 1
          return
        else if ( tval .le. t(2) ) then
          left = 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
c
        low = 2
        high = left - 2

10      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 10
c
c  CASE2: T(LEFT+1) .lt. TVAL:
c  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
c
      else if ( t(left+1) .lt. tval ) then

        if ( left .eq. n - 1 ) then
          return
        else if ( left .eq. n - 2 ) then
          left = left + 1
          return
        else if ( tval .le. t(left+2) ) then
          left = left + 1
          return
        else if ( t(n-1) .le. tval ) then
          left = n - 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
c
        low = left + 2
        high = n - 2

20      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 20
c
c  CASE3: T(LEFT) .le. TVAL .le. T(LEFT+1):
c  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
c
      else

      end if

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
