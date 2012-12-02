      function arc_sine ( s )

c*********************************************************************72
c
cc ARC_SINE computes the arc sine function, with argument truncation.
c
c  Discussion:
c
c    If you call your system ASIN routine with an input argument that is
c    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
c
c    In particular, you may get the value NaN returned.
c
c    This routine truncates arguments outside the range, avoiding the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
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
      subroutine ball_f1_nd ( func, n, center, r, result )

c*********************************************************************72
c
cc BALL_F1_ND approximates an integral inside a ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
c
c  Discussion:
c
c    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision CENTER(N), the center of the ball.
c
c    Input, double precision R, the radius of the ball.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision ball_volume_nd
      double precision center(n)
      double precision func
      external func
      integer i
      integer ihi
      integer itemp
      integer j
      integer k
      integer khi
      integer ktemp
      double precision quad
      double precision r
      double precision result
      double precision t
      double precision temp
      double precision u
      double precision u2
      double precision v
      double precision volume
      double precision w
      double precision x(n)
      double precision y

      if ( r .eq. 0.0D+00 ) then
        result = 0.0D+00
        return
      end if

      u2 = ( 1.0D+00 - 2.0D+00 * sqrt ( 1.0D+00 / dble ( n + 4 ) ) ) 
     &  / dble ( n + 2 )
      u = sqrt ( u2 )

      do i = 1,  n
        x(i) = center(i) - r * u
      end do

      w = 1.0D+00 / dble ( ( n + 1 ) * 2**n )

      quad = 0.0D+00
      ihi = 2**n

      do i = 1, ihi

        itemp = i - 1

        do j = 1, n

          u = ( center(j) - x(j) ) / r

          if ( mod ( itemp, 2 ) .eq. 1 ) then
            x(j) = center(j) - abs ( x(j) - center(j) )
          else
            x(j) = center(j) + abs ( x(j) - center(j) )
          end if

          itemp = itemp / 2

        end do

        quad = quad + w * func ( n, x )

      end do

      temp = sqrt ( dble ( n + 4 ) )

      t = sqrt ( 2.0D+00 * dble ( n + 1 ) / dble ( n + 2 ) ) 
     &  / ( dble ( n ) * temp )

      y = ( 1.0D+00 + 2.0D+00 / ( dble ( n ) * temp ) ) 
     &  / dble ( n + 2 )
      v = sqrt ( y - t )
      u = sqrt ( y + dble ( n - 1 ) * t )

      khi = 2**n

      do i = 1, n

        x(1:n) = center(1:n) - r * v

        x(i) = center(i) - r * u

        do k = 1, khi

          ktemp = k - 1

          do j = 1, n

            if ( mod ( ktemp, 2 ) .eq. 1 ) then
              x(j) = center(j) - abs ( x(j) - center(j) )
            else
              x(j) = center(j) + abs ( x(j) - center(j) )
            end if

            ktemp = ktemp / 2

          end do

          quad = quad + w * func ( n, x )

        end do

        x(i) = center(i) - r * v

      end do

      volume = ball_volume_nd ( n, r )
      result = quad * volume

      return
      end
      subroutine ball_f3_nd ( func, n, center, r, result )

c*********************************************************************72
c
cc BALL_F3_ND approximates an integral inside a ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
c
c  Discussion:
c
c    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision CENTER(N), the center of the ball.
c
c    Input, double precision R, the radius of the ball.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision ball_volume_nd
      double precision center(n)
      double precision func
      external func
      integer i
      integer j
      integer jtemp
      integer k
      double precision quad
      double precision r
      double precision result
      double precision ri
      double precision s
      double precision volume
      double precision weight
      double precision x(n)

      if ( r .eq. 0.0D+00 ) then
        result = 0.0D+00
        return
      end if

      quad = 0.0D+00
c
c  The first point is the center of the ball.
c
      do i = 1, n
        x(i) = center(i)
      end do

      weight = 4.0D+00 / dble ( ( n + 2 ) * ( n + 2 ) )
      quad = quad + weight * func ( n, x )

      s = 1.0D+00 / sqrt ( dble ( n + 4 ) )

      do i = 1, n

        ri = sqrt ( dble ( i + 2 ) / dble ( n + 4 ) )
c
c  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
c
        do j = 1, n

          if ( j .lt. i ) then
            x(j) = center(j)
          else if ( j .eq. i ) then
            x(j) = center(j) + r * ri
          else
            x(j) = center(j) + r * s
          end if

        end do

        weight = 2.0D+00**( i - n ) * dble ( n + 4 ) 
     &    / dble ( ( i + 1 ) * ( i + 2 ) * ( n + 2 ) )
c
c  Now go through all sign permutations of the basic point.
c
        do j = 1, 2**(n+1-i)

          jtemp = j - 1

          do k = i, n

            if ( mod ( jtemp, 2 ) .eq. 1 ) then
              x(k) = center(k) - abs ( x(k) - center(k) )
            else
              x(k) = center(k) + abs ( x(k) - center(k) )
            end if

            jtemp = jtemp / 2

          end do

          quad = quad + weight * func ( n, x )

        end do

      end do

      volume = ball_volume_nd ( n, r )
      result = quad * volume

      return
      end
      function ball_monomial_nd ( n, p, r )

c*********************************************************************72
c
cc BALL_MONOMIAL_ND integrates a monomial on a ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) <= R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gerald Folland,
c    How to Integrate a Polynomial Over a Sphere,
c    American Mathematical Monthly,
c    Volume 108, Number 5, May 2001, pages 446-448.
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, integer P(N), the exponents of X(1) through X(N) 
c    in the monomial.  The exponents P(N) must be nonnegative.
c
c    Input, double precision R, the radius of the ball.
c
c    Output, double precision BALL_MONOMIAL_ND, the integral of
c    X1^P(1) * X2^P(2) * ... * XN^P(N) over the ball.
c
      implicit none

      integer n

      double precision ball_monomial_nd
      integer p(n)
      double precision power
      double precision r
      double precision sphere_unit_monomial_nd

      power = dble ( sum ( p ) + n )

      ball_monomial_nd = sphere_unit_monomial_nd ( n, p ) 
     &  * r**power / power

      return
      end
      subroutine ball_unit_07_3d ( func, result )

c*********************************************************************72
c
cc BALL_UNIT_07_3D approximates an integral inside the unit ball in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z <= 1.
c
c  Discussion:
c
c    A 64 point 7-th degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 4 )

      double precision angle
      double precision ball_unit_volume_3d
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision result
      double precision volume
      double precision w
      double precision weight1(order)
      double precision weight2(order)
      double precision weight3(order)
      double precision x
      double precision xtab1(order)
      double precision xtab2(order)
      double precision xtab3(order)
      double precision y
      double precision z

      save xtab1
      save weight1
c
c  This is the 5 point Gauss-Legendre rule,
c  but with the midpoint deleted, and with different weights.
c
      data xtab1 /
     &  -0.906179845938663992797626878299D+00, 
     &  -0.538469310105683091036314420700D+00, 
     &   0.538469310105683091036314420700D+00, 
     &   0.906179845938663992797626878299D+00 /

      data weight1 /
     &  0.19455533421780251826D+00, 
     &  0.13877799911553081506D+00, 
     &  0.13877799911553081506D+00, 
     &  0.19455533421780251826D+00 /
c
c  Set XTAB2 and WEIGHT2.
c
      do j = 1, order
        angle = pi * dble ( 2 * j - 1 ) / dble ( 2 * order )
        xtab2(j) = cos ( angle )
      end do

      weight2(1:order) = 1.0D+00
c
c  Set XTAB3 and WEIGHT3 for the interval [-1,1].
c
      call legendre_set ( order, xtab3, weight3 )

      w = 3.0D+00 / 16.0D+00

      quad = 0.0D+00

      do i = 1, order
        do j = 1, order
          do k = 1, order

            x = xtab1(i) * sqrt ( 1.0D+00 - xtab2(j) * xtab2(j) ) 
     &                   * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
            y = xtab1(i) * xtab2(j) 
     &        * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
            z = xtab1(i) * xtab3(k)

            quad = quad + w * weight1(i) * weight2(j) * weight3(k) 
     &        * func ( x, y, z )

          end do
        end do
      end do

      volume = ball_unit_volume_3d ( )
      result = quad * volume

      return
      end
      subroutine ball_unit_14_3d ( func, result )

c*********************************************************************72
c
cc BALL_UNIT_14_3D approximates an integral inside the unit ball in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z <= 1.
c
c  Discussion:
c
c    A 288 point 14-th degree formula is used, Stroud number S3:14-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision ball_unit_volume_3d
      double precision func
      external func
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      double precision quad
      double precision r(4)
      double precision result
      double precision temp
      double precision volume
      double precision w1
      double precision w2
      double precision weight(4)
      double precision x
      double precision xtab(5)
      double precision y
      double precision ytab(5)
      double precision z
      double precision ztab(5)

      save r
      save weight
      save xtab
      save ytab
      save ztab

      data r /
     &  0.968160240D+00, 0.836031107D+00, 0.613371433D+00, 
     &  0.324253423D+00 /
      data weight /
     &  0.076181268D+00, 0.126263673D+00, 0.098048133D+00, 
     &  0.032840260D+00 /
      data xtab /
     &  -0.151108275D+00, 0.315838353D+00, 0.346307112D+00, 
     &  -0.101808787D+00, -0.409228403D+00 /
      data ytab /
     &  0.155240600D+00, 0.257049387D+00, 0.666277790D+00, 
     &  0.817386065D+00, 0.501547712D+00 /
      data ztab /
     &  0.976251323D+00, 0.913330032D+00, 0.660412970D+00, 
     &  0.567022920D+00, 0.762221757D+00 /

      quad = 0.0D+00

      do m = 1, 4

        w1 = 125.0D+00 * weight(m) / 3360.0D+00
        x = 0.525731112D+00 * r(m)
        y = 0.850650808D+00 * r(m)
        z = 0.0D+00

        do j = 1, 2
          x = -x
          do k = 1, 2
            y = -y
            do l = 1, 3
              call r8_swap3 ( x, y, z )
              quad = quad + w1 * func ( x, y, z )
            end do
          end do
        end do

        w2 = 143.0D+00 * weight(m) / 3360.0D+00

        do n = 1, 5

          x = xtab(n) * r(m)
          y = ytab(n) * r(m)
          z = ztab(n) * r(m)

          do i = 1, 3

            temp = x
            x = z
            z = -y
            y = -temp

            do j = 1, 3

              call r8_swap3 ( x, y, z )

              quad = quad + w2 * func ( x, y, z )

            end do

            y = -y
            z = -z
            quad = quad + w2 * func ( x, y, z )

          end do

        end do

      end do

      volume = ball_unit_volume_3d ( )
      result = quad * volume

      return
      end
      subroutine ball_unit_15_3d ( func, result )

c*********************************************************************72
c
cc BALL_UNIT_15_3D approximates an integral inside the unit ball in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z <= 1.
c
c  Discussion:
c
c    A 512 point 15-th degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order1
      parameter ( order1 = 4 )
      integer order2
      parameter ( order2 = 8 )

      double precision ball_unit_volume_3d
      double precision cj
      double precision ck
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision result
      double precision sj
      double precision sk
      double precision volume
      double precision w
      double precision weight1(order1)
      double precision weight2(order2)
      double precision x
      double precision xtab1(order1)
      double precision xtab2(order2)
      double precision y
      double precision z

      save weight1
      save xtab1

      data weight1 /
     &  0.0328402599D+00, 0.0980481327D+00, 0.1262636728D+00, 
     &  0.0761812678D+00 /
      data xtab1 /
     &  0.3242534234D+00, 0.6133714327D+00, 0.8360311073D+00, 
     &  0.9681602395D+00 /

      call legendre_set ( order2, xtab2, weight2 )

      w = 3.0D+00 / 32.0D+00

      quad = 0.0D+00

      do i = 1, order1

        do j = 1, order2

          sj = xtab2(j)
          cj = sqrt ( 1.0D+00 - sj * sj )

          do k = 1, 16
            sk = sin ( dble ( k ) * pi / 8.0D+00 )
            ck = cos ( dble ( k ) * pi / 8.0D+00 )
            x = xtab1(i) * cj * ck
            y = xtab1(i) * cj * sk
            z = xtab1(i) * sj
            quad = quad + w * weight1(i) * weight2(j) * func ( x, y, z )
          end do

        end do

      end do

      volume = ball_unit_volume_3d ( )
      result = quad * volume

      return
      end
      subroutine ball_unit_f1_nd ( func, n, result )

c*********************************************************************72
c
cc BALL_UNIT_F1_ND approximates an integral inside the unit ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) <= 1.
c
c  Discussion:
c
c    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision ball_unit_volume_nd
      double precision func
      external func
      integer i
      integer ihi
      integer itemp
      integer j
      integer k
      integer khi
      integer ktemp
      integer l
      double precision quad
      double precision result
      double precision t
      double precision temp
      double precision u
      double precision u2
      double precision v
      double precision volume
      double precision w
      double precision x(n)
      double precision y

      u2 = ( 1.0D+00 - 2.0D+00 * sqrt ( 1.0D+00 / dble ( n + 4 ) ) ) 
     &  / dble ( n + 2 )
      u = sqrt ( u2 )

      do i = 1, n
        x(i) = -u
      end do

      w = 1.0D+00 / dble ( ( n + 1 ) * 2**n )

      quad = 0.0D+00
      ihi = 2**n

      do i = 1, ihi

        itemp = i - 1

        do j = 1, n

          if ( mod ( itemp, 2 ) .eq. 1 ) then
            x(j) = -abs ( x(j) )
          else
            x(j) = abs ( x(j) )
          end if

          itemp = itemp / 2

        end do

        quad = quad + w * func ( n, x )

      end do

      temp = sqrt ( dble ( n + 4 ) )

      t = sqrt ( 2.0D+00 * dble ( n + 1 ) / dble ( n + 2 ) ) 
     &  / ( dble ( n ) * temp )

      y = ( 1.0D+00 + 2.0D+00 / ( dble ( n ) * temp ) ) 
     &  / dble ( n + 2 )
      v = sqrt ( y - t )
      u = sqrt ( y + dble ( n - 1 ) * t )

      khi = 2**n

      do i = 1, n

        do l = 1, n
          x(l) = -v
        end do

        x(i) = -u

        do k = 1, khi

          ktemp = k - 1

          do j = 1, n

            if ( mod ( ktemp, 2 ) .eq. 1 ) then
              x(j) = -abs ( x(j) )
            else
              x(j) = abs ( x(j) )
            end if

            ktemp = ktemp / 2

          end do

          quad = quad + w * func ( n, x )

        end do

        x(i) = -v

      end do

      volume = ball_unit_volume_nd ( n )
      result = quad * volume

      return
      end
      subroutine ball_unit_f3_nd ( func, n, result )

c*********************************************************************72
c
cc BALL_UNIT_F3_ND approximates an integral inside the unit ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) <= 1.
c
c  Discussion:
c
c    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision ball_unit_volume_nd
      double precision func
      external func
      integer i
      integer j
      integer jtemp
      integer k
      double precision quad
      double precision result
      double precision ri
      double precision s
      double precision volume
      double precision weight
      double precision x(n)

      quad = 0.0D+00
c
c  The first point is the center of the ball.
c
      do i = 1, n
        x(i) = 0.0D+00
      end do

      weight = 4.0D+00 / dble ( ( n + 2 ) * ( n + 2 ) )
      quad = quad + weight * func ( n, x )

      s = 1.0D+00 / sqrt ( dble ( n + 4 ) )

      do i = 1, n

        ri = sqrt ( dble ( i + 2 ) / dble ( n + 4 ) )
c
c  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
c
        do j = 1, n

          if ( j .lt. i ) then
            x(j) = 0.0D+00
          else if ( j .eq. i ) then
            x(j) = ri
          else
            x(j) = s
          end if

        end do

        weight = 2.0D+00**( i - n ) * dble ( n + 4 ) 
     &    / dble ( ( i + 1 ) * ( i + 2 ) * ( n + 2 ) )
c
c  Now go through all sign permutations of the basic point.
c
        do j = 1, 2**(n+1-i)

          jtemp = j - 1

          do k = i, n

            if ( mod ( jtemp, 2 ) .eq. 1 ) then
              x(k) = -abs ( x(k) )
            else
              x(k) = abs ( x(k) )
            end if

            jtemp = jtemp / 2

          end do

          quad = quad + weight * func ( n, x )

        end do

      end do

      volume = ball_unit_volume_nd ( n )
      result = quad * volume

      return
      end
      function ball_unit_volume_3d ( )

c*********************************************************************72
c
cc BALL_UNIT_VOLUME_3D computes the volume of the unit ball in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision BALL_UNIT_VOLUME_3D, the volume of the ball.
c
      implicit none

      double precision ball_unit_volume_3d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      ball_unit_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi

      return
      end
      function ball_unit_volume_nd ( n )

c*********************************************************************72
c
cc BALL_UNIT_VOLUME_ND computes the volume of the unit ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) <= 1.
c
c  Discussion:
c
c    N  Volume
c
c    2             PI
c    3  (4/3)    * PI
c    4  (1/2)    * PI^2
c    5  (8/15)   * PI^2
c    6  (1/6)    * PI^3
c    7  (16/105) * PI^3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision BALL_UNIT_VOLUME_ND, the volume of the ball.
c
      implicit none

      double precision ball_unit_volume_nd
      integer i
      integer m
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision volume

      if ( mod ( n, 2 ) .eq. 0 ) then
        m = n / 2
        volume = ( pi )**m
        do i = 1, m
          volume = volume / dble ( i )
        end do
      else
        m = ( n - 1 ) / 2
        volume = ( pi )**m * 2.0D+00**n
        do i = m+1, 2*m+1
          volume = volume / dble ( i )
        end do
      end if

      ball_unit_volume_nd = volume

      return
      end
      function ball_volume_3d ( r )

c*********************************************************************72
c
cc BALL_VOLUME_3D computes the volume of a ball in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z <= R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the ball.
c
c    Output, double precision BALL_VOLUME_3D, the volume of the ball.
c
      implicit none

      double precision ball_volume_3d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      ball_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi * r**3

      return
      end
      function ball_volume_nd ( n, r )

c*********************************************************************72
c
cc BALL_VOLUME_ND computes the volume of a ball in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) <= R * R
c
c  Discussion:
c
c    N  Volume
c
c    2             PI   * R^2
c    3  (4/3)    * PI   * R^3
c    4  (1/2)    * PI^2 * R^4
c    5  (8/15)   * PI^2 * R^5
c    6  (1/6)    * PI^3 * R^6
c    7  (16/105) * PI^3 * R^7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision R, the radius of the ball.
c
c    Output, double precision BALL_VOLUME_ND, the volume of the ball.
c
      implicit none

      double precision ball_unit_volume_nd
      double precision ball_volume_nd
      integer n
      double precision r

      ball_volume_nd = ball_unit_volume_nd ( n ) * r**n

      return
      end
      subroutine c1_geg_monomial_integral ( alpha, expon, value )

c*********************************************************************72
c
cc C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
c
c  Discussion:
c
c    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
c
c      w(alpha;x) = (1-x^2)^alpha
c
c    with -1.0 < alpha.
c
c    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the exponent of (1-X^2).
c    - 1.0 < ALPHA.
c
c    Input, integer EXPON, the exponent.
c    0 <= EXPON.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      double precision alpha
      double precision arg1
      double precision arg2
      double precision arg3
      double precision arg4
      double precision c
      integer expon
      double precision r8_gamma
      double precision value
      double precision value1

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C1_GEG_MONOMIAL_INTEGRAL - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      if ( mod ( expon, 2 ) .eq. 1 ) then
        value = 0.0D+00
        return
      end if

      c = dble ( expon )

      arg1 = - alpha
      arg2 =   1.0D+00 + c
      arg3 =   2.0D+00 + alpha + c
      arg4 = - 1.0D+00

      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

      value = 2.0D+00 * r8_gamma ( 1.0D+00 + c ) 
     &  * r8_gamma ( 1.0D+00 + alpha ) 
     &  * value1 / r8_gamma ( 2.0D+00 + alpha  + c )

      return
      end
      subroutine c1_jac_monomial_integral ( alpha, beta, expon, value )

c*********************************************************************72
c
cc C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
c
c  Discussion:
c
c    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the exponent of (1-X) in the weight factor.
c
c    Input, double precision BETA, the exponent of (1+X) in the weight factor.
c
c    Input, integer EXPON, the exponent.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      double precision alpha
      double precision arg1
      double precision arg2
      double precision arg3
      double precision arg4
      double precision beta
      double precision c
      integer expon
      double precision r8_gamma
      double precision s
      double precision value
      double precision value1
      double precision value2

      c = dble ( expon )

      if ( mod ( expon, 2 ) .eq. 0 ) then
        s = +1.0D+00
      else
        s = -1.0D+00
      end if

      arg1 = - alpha
      arg2 =   1.0D+00 + c
      arg3 =   2.0D+00 + beta + c
      arg4 = - 1.0D+00

      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

      arg1 = - beta
      arg2 =   1.0D+00 + c
      arg3 =   2.0D+00 + alpha + c
      arg4 = - 1.0D+00

      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

      value = r8_gamma ( 1.0D+00 + c ) * ( 
     &    s * r8_gamma ( 1.0D+00 + beta  ) * value1 
     &  / r8_gamma ( 2.0D+00 + beta  + c ) 
     &  +     r8_gamma ( 1.0D+00 + alpha ) * value2 
     &  / r8_gamma ( 2.0D+00 + alpha + c ) )

      return
      end
      subroutine c1_leg_monomial_integral ( expon, value )

c*********************************************************************72
c
cc C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
c
c  Discussion:
c
c    C1_LEG is the interval [-1,+1] with the Legendre weight function
c
c      w(x) = 1.
c
c    value = integral ( -1 <= x <= +1 ) x^expon dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer EXPON, the exponent.
c    0 <= EXPON.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer expon
      double precision value

      if ( mod ( expon, 2 ) .eq. 1 ) then
        value = 0.0D+00
        return
      end if

      if ( expon .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C1_LEG_MONOMIAL_INTEGRAL - Fatal error!'
        write ( *, '(a)' ) '  EXPON < 0.'
        stop
      end if

      value = 2.0D+00 / dble ( expon + 1 )

      return
      end
      subroutine circle_annulus ( func, center, radius1, radius2, nr, 
     &  result )

c*********************************************************************72
c
cc CIRCLE_ANNULUS approximates an integral in an annulus.
c
c  Integration region:
c
c    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Peirce,
c    Numerical Integration Over the Planar Annulus,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 5, Number 2, June 1957, pages 66-73.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of two
c    variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the circle.
c
c    Input, double precision RADIUS1, RADIUS2, the radii of the circles.
c
c    Input, integer NR, the order of the rule.  This quantity 
c    specifies the number of distinct radii to use.  The number of angles used 
c    will be 4*NR, for a total of 4*NR*NR points.
c
c    Output, double precision RESULT, the approximation to the integral.
c
      implicit none

      integer nr
      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision area
      double precision b
      double precision c
      double precision center(dim_num)
      double precision circle_annulus_area_2d
      double precision d
      double precision func
      external func
      integer i
      integer j
      integer nt
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision ra(nr)
      double precision radius1
      double precision radius2
      double precision result
      double precision rw(nr)
      double precision t
      double precision tw
      double precision x
      double precision y
c
c  Choose radial abscissas and weights.
c
      call legendre_set ( nr, ra, rw )
      a = -1.0D+00
      b = +1.0D+00
      c = radius1 * radius1
      d = radius2 * radius2
      call rule_adjust ( a, b, c, d, nr, ra, rw )
      do i = 1, nr
        ra(i) = sqrt ( ra(i) )
        rw(i) = rw(i) / ( radius2 - radius1 ) / ( radius2 + radius1 )
      end do
c
c  Set angular abscissas and weights.
c
      nt = 4 * nr

      tw = 1.0D+00 / dble ( nt )
c
c  Approximate the integral.
c
      quad = 0.0D+00
      do i = 1, nt
        t = 2.0D+00 * pi * dble ( i - 1 ) / dble ( nt )
        do j = 1, nr
          x = center(1) + ra(j) * cos ( t )
          y = center(2) + ra(j) * sin ( t )
          quad = quad + tw * rw(j) * func ( x, y )
        end do
      end do

      area = circle_annulus_area_2d ( radius1, radius2 )
      result = quad * area

      return
      end
      function circle_annulus_area_2d ( radius1, radius2 )

c*********************************************************************72
c
cc CIRCLE_ANNULUS_AREA_2D returns the area of a circular annulus in 2D.
c
c  Integration region:
c
c    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RADIUS1, RADIUS2, the radii of the circles.
c
c    Output, double precision CIRCLE_ANNULUS_AREA_2D, the area of the annulus.
c
      implicit none

      double precision circle_annulus_area_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision radius1
      double precision radius2

      circle_annulus_area_2d = pi * ( radius1 + radius2 ) 
     &  * ( radius2 - radius1 )

      return
      end
      subroutine circle_annulus_sector ( func, center, radius1, 
     &  radius2, theta1, theta2, nr, result )

c*********************************************************************72
c
cc CIRCLE_ANNULUS_SECTOR approximates an integral in a circular annulus sector.
c
c  Discussion:
c
c    A circular annulus sector comprises the area between two concentric
c    circles and two concentric rays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Peirce,
c    Numerical Integration Over the Planar Annulus,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 5, Number 2, June 1957, pages 66-73.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of two
c    variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the circle.
c
c    Input, double precision RADIUS1, RADIUS2, the radii of the circles.
c
c    Input, double precision THETA1, THETA2, the angles defining the sector.
c    The sector is measured from THETA1 to THETA2.
c
c    Input, integer NR, the order of the rule.  This quantity 
c    specifies the number of distinct radii to use.  The number of angles used 
c    will be 4*NR, for a total of 4*NR*NR points.
c
c    Output, double precision RESULT, the approximation to the integral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer nr

      double precision a
      double precision area
      double precision b
      double precision c
      double precision center(dim_num)
      double precision circle_annulus_sector_area_2d
      double precision d
      double precision func
      external func
      integer i
      integer j
      integer nt
      double precision quad
      double precision ra(nr)
      double precision radius1
      double precision radius2
      double precision result
      double precision rw(nr)
      double precision ta(4*nr)
      double precision theta1
      double precision theta2
      double precision tw(4*nr)
      double precision x
      double precision y
c
c  Set the radial abscissas and weights.
c
      call legendre_set ( nr, ra, rw )
      a = -1.0D+00
      b = +1.0D+00
      c = radius1 * radius1
      d = radius2 * radius2
      call rule_adjust ( a, b, c, d, nr, ra, rw )
      ra(1:nr) = sqrt ( ra(1:nr) )
      rw(1:nr) = rw(1:nr) / ( radius2 - radius1 ) 
     &  / ( radius2 + radius1 )
c
c  Pick angles evenly spaced between THETA1 and THETA2, but do not
c  include the endpoints, and use a half interval for the first and last.
c
      nt = 4 * nr

      call tvec_even_bracket3 ( nt, theta1, theta2, ta )
      tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )
c
c  Approximate the integral.
c
      quad = 0.0D+00
      do i = 1, nt
        do j = 1, nr
          x = center(1) + ra(j) * cos ( ta(i) )
          y = center(2) + ra(j) * sin ( ta(i) )
          quad = quad + tw(i) * rw(j) * func ( x, y )
        end do
      end do

      area = circle_annulus_sector_area_2d ( radius1, radius2, theta1, 
     &  theta2 )
      result = quad * area

      return
      end
      function circle_annulus_sector_area_2d ( radius1, radius2, 
     &  theta1, theta2 )

c*********************************************************************72
c
cc CIRCLE_ANNULUS_SECTOR_AREA_2D: area of a circular annulus sector in 2D.
c
c  Discussion:
c
c    A circular annulus sector comprises the area between two concentric
c    circles and two concentric rays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RADIUS1, RADIUS2, the radii of the circles.
c
c    Input, double precision THETA1, THETA2, the angles of the rays.
c    Ordinarily, (THETA2-THETA1) is between 0 and 2*PI.
c
c    Output, double precision CIRCLE_ANNULUS_SECTOR_AREA_2D, the area of the
c    circular annulus sector.
c
      implicit none

      double precision circle_annulus_sector_area_2d
      double precision radius1
      double precision radius2
      double precision theta1
      double precision theta2

      circle_annulus_sector_area_2d = 0.5D+00 * ( radius1 + radius2 ) 
     &  * ( radius2 - radius1 ) * ( theta2 - theta1 )

      return
      end
      function circle_area_2d ( r )

c*********************************************************************72
c
cc CIRCLE_AREA_2D returns the area of a circle in 2D.
c
c  Integration region:
c
c    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Output, double precision CIRCLE_AREA_2D, the area of the circle.
c
      implicit none

      double precision circle_area_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      circle_area_2d = pi * r * r

      return
      end
      subroutine circle_cap_area_2d ( r, h, area )

c*********************************************************************72
c
cc CIRCLE_CAP_AREA_2D computes the area of a circle cap in 2D.
c
c  Discussion:
c
c    Draw any radius R of the circle and denote as P the point where the
c    radius intersects the circle.  Now consider the point Q which lies
c    on the radius and which is H units from P.  The line which is
c    perpendicular to the radius R and passes through Q divides the
c    circle into two pieces.  The piece including the point P is the
c    circular cap of height (or thickness) H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision H, the "height" of the circle cap.  
c
c    Output, double precision AREA, the area of the circle cap.
c
      implicit none

      double precision arc_sine
      double precision area
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      if ( h .le. 0.0D+00 ) then

        area = 0.0D+00

      else if ( h .le. r ) then

        theta = 2.0D+00 
     &    * arc_sine ( sqrt ( h * ( 2.0D+00 * r - h ) ) / r )
        area = r * r * ( theta - sin ( theta ) ) / 2.0D+00

      else if ( h .le. 2.0D+00 * r ) then

        theta = 2.0D+00 
     &    * arc_sine ( sqrt ( h * ( 2.0D+00 * r - h ) ) / r )
        area = r * r * ( pi - ( theta - sin ( theta ) ) / 2.0D+00 )

      else if ( 2.0D+00 * r .le. h ) then

        area = pi * r * r

      end if

      return
      end
      subroutine circle_cum ( func, center, radius, order, result )

c*********************************************************************72
c
cc CIRCLE_CUM approximates an integral on the circumference of a circle in 2D.
c
c  Integration region:
c
c    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
c
c  Discussion:
c
c    An ORDER point, (ORDER-1)-th degree formula is used, 
c    Stroud number U2:M-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of two
c    variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the coordinates of the center of 
c    the circle.
c
c    Input, double precision RADIUS, the radius of the circle.
c
c    Input, integer ORDER, the number of points to use.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle
      double precision center(dim_num)
      double precision func
      external func
      integer i
      integer order
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision radius
      double precision result
      double precision volume
      double precision x
      double precision y

      quad = 0.0D+00

      do i = 1, order
        angle = dble ( 2 * i ) * pi / dble ( order )
        x = center(1) + radius * cos ( angle )
        y = center(2) + radius * sin ( angle )
        quad = quad + func ( x, y )
      end do

      quad = quad / dble ( order )

      volume = pi * radius * radius
      result = quad * volume

      return
      end
      subroutine circle_rt_set ( rule, nr, nt, nc, ra, rw, ta, tw, cw )

c*********************************************************************72
c
cc CIRCLE_RT_SET sets an R, THETA product quadrature rule in the unit circle.
c
c  Discussion:
c
c    For a given value of RULE, here are the number of points used at the
c    center (NC), the number of points along the radial direction (NR) and
c    the number of points along the theta direction (NT).  The total number
c    of points in the rule will be 
c
c      Total = NC + NR * NT.
c
c    The user, when choosing RULE, must allocate enough space in the arrays
c    RA, RW, TA and TW for the resulting values of NR and NT.
c
c    RULE  NC  NR  NT  Total
c    ----  --  --  --  -----
c       1   1   0   0      1
c       2   0   1   4      4
c       3   1   1   4      5
c       4   1   1   6      7
c       5   1   2   4      9
c       6   0   3   4     12
c       7   1   2  10     21
c       8   0   4  16     64
c       9   0   5  20    120
c
c    The integral of F(X,Y) over the unit circle is approximated by
c
c      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
c      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
c      = approximately
c        CW * F(0,0) 
c        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
c        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c
c    Input, integer NR, the number of R abscissas.
c
c    Input, integer NT, the number of Theta abscissas.
c
c    Input, integer NC, the number of center abscissas (0 or 1 ).
c
c    Output, double precision RA(NR), RW(NR), the R abscissas and weights.
c
c    Output, double precision TA(NT), TW(NT), the THETA abscissas and weights.
c
c    Output, double precision ZW, the weight to use for the center.
c
      implicit none

      integer nr
      integer nt

      double precision a
      double precision b
      double precision c
      double precision cw
      double precision d
      integer i
      integer nc
      double precision ra(nr)
      double precision rw(nr)
      integer rule
      double precision ta(nt)
      double precision tw(nt)
      double precision u
      double precision v
      double precision w

      if ( rule .eq. 1 ) then

        cw = 1.0D+00

      else if ( rule .eq. 2 ) then

        ra(1) = 0.5D+00
        rw(1) = 1.0D+00

        call tvec_even2 ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 0.0D+00

      else if ( rule .eq. 3 ) then

        ra(1) = 1.0D+00
        rw(1) = 1.0D+00

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 0.125D+00
        end do
        cw = 0.5D+00

      else if ( rule .eq. 4 ) then

        ra(1) = sqrt ( 2.0D+00 / 3.0D+00 )
        rw(1) = 1.0D+00

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 0.125D+00
        end do
        cw = 0.25D+00

      else if ( rule .eq. 5 ) then

        a = 1.0D+00
        b = sqrt ( 2.0D+00 ) / 2.0D+00
        u = 1.0D+00 / 6.0D+00
        v = 4.0D+00 / 6.0D+00

        ra(1) = a
        ra(2) = b
        rw(1) = u
        rw(2) = v

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 4.0D+00 / 24.0D+00

      else if ( rule .eq. 6 ) then

        a = sqrt ( 3.0D+00 ) / 2.0D+00
        b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &    / 52.0D+00 )
        c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &    / 52.0D+00 )

        u = 8.0D+00 / 27.0D+00
        v = ( 551.0D+00 + 41.0D+00 * sqrt ( 29.0D+00 ) ) / 1566.0D+00
        w = ( 551.0D+00 - 41.0D+00 * sqrt ( 29.0D+00 ) ) / 1566.0D+00

        ra(1) = a
        ra(2) = b
        ra(3) = c
        rw(1) = u
        rw(2) = v
        rw(3) = w

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 0.0D+00

      else if ( rule .eq. 7 ) then

        a = sqrt ( ( 6.0D+00 - sqrt ( 6.0D+00 ) ) / 10.0D+00 )
        b = sqrt ( ( 6.0D+00 + sqrt ( 6.0D+00 ) ) / 10.0D+00 )
        u = ( 16.0D+00 + sqrt ( 6.0D+00 ) ) / 36.0D+00
        v = ( 16.0D+00 - sqrt ( 6.0D+00 ) ) / 36.0D+00

        ra(1) = a
        ra(2) = b
        rw(1) = u
        rw(2) = v

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 1.0D+00 / 9.0D+00

      else if ( rule .eq. 8 ) then

        call legendre_set ( nr, ra, rw )
        a = -1.0D+00
        b = +1.0D+00
        c =  0.0D+00
        d = +1.0D+00
        call rule_adjust ( a, b, c, d, nr, ra, rw )
        do i = 1, nr
          ra(i) = sqrt ( ra(i) )
        end do

        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 0.0D+00

      else if ( rule .eq. 9 ) then

        call legendre_set ( nr, ra, rw )
        a = -1.0D+00
        b = +1.0D+00
        c =  0.0D+00
        d = +1.0D+00
        call rule_adjust ( a, b, c, d, nr, ra, rw )
        do i = 1, nr
          ra(i) = sqrt ( ra(i) )
        end do
        call tvec_even ( nt, ta )
        do i = 1, nt
          tw(i) = 1.0D+00 / dble ( nt )
        end do
        cw = 0.0D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CIRCLE_RT_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  There is no rule of index ', rule
        stop

      end if

      return
      end
      subroutine circle_rt_size ( rule, nr, nt, nc )

c*********************************************************************72
c
cc CIRCLE_RT_SIZE sizes an R, THETA product quadrature rule in the unit circle.
c
c  Discussion:
c
c    For a given value of RULE, here are the number of points used at the
c    center (NC), the number of points along the radial direction (NR) and
c    the number of points along the theta direction (NT).  The total number
c    of points in the rule will be 
c
c      Total = NC + NR * NT.
c
c    The user, when choosing RULE, must allocate enough space in the arrays
c    RA, RW, TA and TW for the resulting values of NR and NT.
c
c    RULE  NC  NR  NT  Total
c    ----  --  --  --  -----
c       1   1   0   0      1
c       2   0   1   4      4
c       3   1   1   4      5
c       4   1   1   6      7
c       5   1   2   4      9
c       6   0   3   4     12
c       7   1   2  10     21
c       8   0   4  16     64
c       9   0   5  20    120
c
c    The integral of F(X,Y) over the unit circle is approximated by
c
c      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
c      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
c      = approximately
c        ZW * F(0,0) 
c        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
c        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c
c    Output, integer NR, the number of R abscissas.
c    
c    Output, integer NT, the number of Theta abscissas.
c
c    Output, integer NC, the number of center abscissas (0 or 1).
c
      implicit none

      integer nc
      integer nr
      integer nt
      integer rule

      if ( rule .eq. 1 ) then

        nr = 0
        nt = 0
        nc = 1

      else if ( rule .eq. 2 ) then

        nr = 1
        nt = 4
        nc = 0

      else if ( rule .eq. 3 ) then

        nr = 1
        nt = 4
        nc = 1

      else if ( rule .eq. 4 ) then

        nr = 1
        nt = 6
        nc = 1

      else if ( rule .eq. 5 ) then

        nr = 2
        nt = 4
        nc = 1

      else if ( rule .eq. 6 ) then

        nr = 3
        nt = 4
        nc = 0

      else if ( rule .eq. 7 ) then

        nr = 2
        nt = 10
        nc = 1

      else if ( rule .eq. 8 ) then

        nr = 4
        nt = 16
        nc = 0

      else if ( rule .eq. 9 ) then

        nr = 5
        nt = 20
        nc = 0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CIRCLE_RT_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  There is no rule of index ', rule
        stop

      end if

      return
      end
      subroutine circle_rt_sum ( func, center, radius, nr, ra, rw, 
     &  nt, ta, tw, zw, result )

c*********************************************************************72
c
cc CIRCLE_RT_SUM applies an R, THETA product quadrature rule inside a circle.
c
c  Integration region:
c
c    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= RADIUS^2.
c
c  Discussion:
c
c    The product rule is assumed to be have the form:
c
c      Integral_Approx = ZW * F(CENTER(1),CENTER(2)) +
c        sum ( 1 <= IR <= NR ) Sum ( 1 <= IT <= NT )
c        RW(IR) * TW(IT) * F ( CENTER(1) + R(IR) * RADIUS * Cos ( TA(IT) ),
c                              CENTER(2) + R(IR) * RADIUS * Sin ( TA(IT) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of two variables which is to be integrated,
c    of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the circle.
c
c    Input, double precision RADIUS, the radius of the circle.
c
c    Input, integer NR, the number of R abscissas.
c
c    Input, double precision RA(NR), RW(NR), the R abscissas and weights.
c
c    Input, integer NT, the number of Theta abscissas.
c
c    Input, double precision TA(NT), TW(NT), the THETA abscissas and weights.
c
c    Input, double precision ZW, the weight to use for the center.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer nr
      integer nt

      double precision center(dim_num)
      double precision circle_area_2d
      double precision func
      external func
      integer ir
      integer it
      double precision quad
      double precision ra(nr)
      double precision radius
      double precision rct
      double precision result
      double precision rst
      double precision rw(nr)
      double precision ta(nt)
      double precision tw(nt)
      double precision volume
      double precision x
      double precision y
      double precision zw

      quad = 0.0D+00

      if ( zw .ne. 0.0D+00 ) then
        x = center(1)
        y = center(2)
        quad = quad + zw * func ( x, y )
      end if

      do it = 1, nt
        rct = radius * cos ( ta(it) )
        rst = radius * sin ( ta(it) )
        do ir = 1, nr
          x = center(1) + ra(ir) * rct
          y = center(2) + ra(ir) * rst
          quad = quad + tw(it) * rw(ir) * func ( x, y )
        end do
      end do

      volume = circle_area_2d ( radius )
      result = quad * volume

      return
      end
      subroutine circle_sector ( func, center, radius, theta1, theta2, 
     &  nr, result )

c*********************************************************************72
c
cc CIRCLE_SECTOR approximates an integral in a circular sector.
c
c  Discussion:
c
c    A sector is contained within a circular arc and the lines joining each
c    endpoint of the arc to the center of the circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of two
c    variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the circle.
c
c    Input, double precision RADIUS, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the sector.
c    The sector is measured from THETA1 to THETA2.
c
c    Input, integer NR, the number of radial values used in the 
c    approximation of the integral.  NR must be at least 1.  Higher values 
c    improve the accuracy of the integration, at the cost of more function 
c    evaluations.
c
c    Output, double precision RESULT, the approximation to the integral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer nr

      double precision a
      double precision area
      double precision b
      double precision c
      double precision center(dim_num)
      double precision circle_sector_area_2d
      double precision d
      double precision func
      external func
      integer i
      integer j
      integer nt
      double precision quad
      double precision ra(nr)
      double precision radius
      double precision result
      double precision rw(nr)
      double precision ta(4*nr)
      double precision theta1
      double precision theta2
      double precision tw(4*nr)
      double precision x
      double precision y
c
c  Set the radial abscissas and weights.
c
      call legendre_set ( nr, ra, rw )
      a = -1.0D+00
      b = +1.0D+00
      c =  0.0D+00
      d =  radius * radius
      call rule_adjust ( a, b, c, d, nr, ra, rw )
      do i = 1, nr
        ra(i) = sqrt ( ra(i) )
        rw(i) = rw(i) / radius / radius
      end do
c
c  Pick angles evenly spaced between THETA1 and THETA2, but do not
c  include the endpoints, and use a half interval for the first and last.
c
      nt = 4 * nr

      call tvec_even_bracket3 ( nt, theta1, theta2, ta )
      do i = 1, nt
        tw(i) = 1.0D+00 / dble ( nt )
      end do
c
c  Approximate the integral.
c
      quad = 0.0D+00

      do i = 1, nr
        do j = 1, nt
          x = center(1) + ra(i) * cos ( ta(j) )
          y = center(2) + ra(i) * sin ( ta(j) )
          quad = quad + rw(i) * tw(j) * func ( x, y )
        end do
      end do

      area = circle_sector_area_2d ( radius, theta1, theta2 )
      result = quad * area

      return
      end
      function circle_sector_area_2d ( r, theta1, theta2 )

c*********************************************************************72
c
cc CIRCLE_SECTOR_AREA_2D returns the area of a circular sector in 2D.
c
c  Discussion:
c
c    A sector is contained within a circular arc and the lines joining each
c    endpoint of the arc to the center of the circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles of the rays
c    that delimit the sector.
c
c    Output, double precision CIRCLE_SECTOR_AREA_2D, the area of the sector.
c
      implicit none

      double precision circle_sector_area_2d
      double precision r
      double precision theta1
      double precision theta2

      circle_sector_area_2d = 0.50D+00 * r * r * ( theta2 - theta1 )

      return
      end
      subroutine circle_sector_rule ( radius, theta1, theta2, nr, nt, 
     &  ra, rw, ta, tw )

c*********************************************************************72
c
cc CIRCLE_SECTOR_RULE approximates an integral in a circular sector.
c
c  Discussion:
c
c    A sector is contained within a circular arc and the lines joining each
c    endpoint of the arc to the center of the circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RADIUS, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the sector.
c    The sector is measured from THETA1 to THETA2.
c
c    Input, integer NR, the number of radial values.
c
c    Input, integer NT, the number of angular values.
c
c    Output, double precision RA(NR), RW(NR), the radial abscissas and weights.
c
c    Output, double precision TA(NT), TW(NT), the angular abscissas 
c    and weights.
c
      implicit none

      integer nr
      integer nt

      double precision a
      double precision area
      double precision b
      double precision c
      double precision d
      integer i
      integer j
      double precision ra(nr)
      double precision radius
      double precision rw(nr)
      double precision ta(nt)
      double precision theta1
      double precision theta2
      double precision tw(nt)
c
c  Set the radial abscissas and weights.
c
      call legendre_set ( nr, ra, rw )
      a = -1.0D+00
      b = +1.0D+00
      c =  0.0D+00
      d =  radius * radius
      call rule_adjust ( a, b, c, d, nr, ra, rw )
      do i = 1, nr
        ra(i) = sqrt ( ra(i) )
        rw(i) = rw(i) / radius / radius
      end do
c
c  Pick angles evenly spaced between THETA1 and THETA2, but do not
c  include the endpoints, and use a half interval for the first and last.
c
      call tvec_even_bracket3 ( nt, theta1, theta2, ta )
      do i = 1, nt
        tw(i) = 1.0D+00 / dble ( nt )
      end do

      return
      end
      function circle_triangle_area_2d ( r, theta1, theta2 )

c*********************************************************************72
c
cc CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
c
c  Discussion:
c
c    A circle triangle is formed by drawing a circular arc, and considering
c    the triangle formed by the endpoints of the arc plus the center of
c    the circle.
c
c    The normal situation is that 0 < ( THETA2 - THETA1 ) < PI.  Outside
c    this range, the triangle can actually have NEGATIVE area.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles of the rays that
c    delimit the arc.
c
c    Output, double precision CIRCLE_TRIANGLE_AREA_2D, the (signed) area
c    of the triangle.
c
      implicit none

      double precision circle_triangle_area_2d
      double precision r
      double precision theta1
      double precision theta2

      circle_triangle_area_2d = 
     &  0.5D+00 * r * r * sin ( theta2 - theta1 )

      return
      end
      subroutine circle_xy_set ( rule, order, xtab, ytab, weight )

c*********************************************************************72
c
cc CIRCLE_XY_SET sets an XY quadrature rule inside the unit circle in 2D.
c
c  Integration region:
c
c    X*X + Y*Y <= 1.0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Frank Lether,
c    A Generalized Product Rule for the Circle,
c    SIAM Journal on Numerical Analysis,
c    Volume 8, Number 2, June 1971, pages 249-253.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c      1, 1 point 1-st degree;
c      2, 4 point 3-rd degree, Stroud S2:3-1;
c      3, 4 point 3-rd degree, Lether #1;
c      4, 4 point 3-rd degree, Stroud S2:3-2;
c      5, 5 point 3-rd degree;
c      6, 7 point 5-th degree;
c      7, 9 point 5-th degree;
c      8, 9 point 5-th degree, Lether #2;
c      9, 12 point 7-th degree;
c     10, 16 point 7-th degree, Lether #3;
c     11, 21 point 9-th degree, Stroud S2:9-3;
c     12, 25 point 9-th degree, Lether #4 (after correcting error);
c     13, 64 point 15-th degree Gauss product rule.
c
c    Input, integer ORDER, the order of the desired rule.
c
c    Output, double precision XTAB(ORDER), YTAB(ORDER), the abscissas of 
c    the rule.
c
c    Output, double precision WEIGHT(ORDER), the ORDER weights of the rule.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision f
      double precision g
      double precision h
      integer i
      integer j
      integer k
      integer nr
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision ra(4)
      double precision rw(4)
      integer rule
      double precision s
      double precision u
      double precision v
      double precision w
      double precision w1
      double precision w2
      double precision w3
      double precision w4
      double precision w5
      double precision w6
      double precision w7
      double precision w8
      double precision w9
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)
      double precision z

      if ( rule .eq. 1 ) then

        xtab(1) = 0.0D+00
        ytab(1) = 0.0D+00
        weight(1) = 1.0D+00

      else if ( rule .eq. 2 ) then

        a = 0.5D+00
        b = 0.25D+00
        z = 0.0D+00

        xtab(1) =  a
        xtab(2) = -a
        xtab(3) =  z
        xtab(4) =  z

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  a
        ytab(4) = -a

        weight(1) = b
        weight(2) = b
        weight(3) = b
        weight(4) = b

      else if ( rule .eq. 3 ) then

        a = 0.5D+00
        b = 0.25D+00

        xtab(1) =  a
        xtab(2) = -a
        xtab(3) = -a
        xtab(4) =  a

        ytab(1) =  a
        ytab(2) =  a
        ytab(3) = -a
        ytab(4) = -a

        weight(1) = b
        weight(2) = b
        weight(3) = b
        weight(4) = b

      else if ( rule .eq. 4 ) then

        a = sqrt ( 2.0D+00 ) / 2.0D+00
        b = 0.25D+00

        xtab(1) =  a
        xtab(2) = -a
        xtab(3) = -a
        xtab(4) =  a

        ytab(1) =  a
        ytab(2) =  a
        ytab(3) = -a
        ytab(4) = -a

        weight(1) = b
        weight(2) = b
        weight(3) = b
        weight(4) = b

      else if ( rule .eq. 5 ) then

        a = 1.0D+00
        b = 0.5D+00
        c = 0.125D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  a
        xtab(3) =  z
        xtab(4) = -a
        xtab(5) =  z

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  a
        ytab(4) =  z
        ytab(5) = -a

        weight(1) = b
        weight(2) = c
        weight(3) = c
        weight(4) = c
        weight(5) = c

      else if ( rule .eq. 6 ) then

        a = sqrt ( 2.0D+00 / 3.0D+00 )
        b = sqrt ( 1.0D+00 / 6.0D+00 )
        c = sqrt ( 2.0D+00 ) / 2.0D+00
        d = 0.125D+00
        e = 0.25D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  a
        xtab(3) = -a
        xtab(4) =  b
        xtab(5) = -b
        xtab(6) =  b
        xtab(7) = -b

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  c
        ytab(5) =  c
        ytab(6) = -c
        ytab(7) = -c

        weight(1) = e
        weight(2) = d
        weight(3) = d
        weight(4) = d
        weight(5) = d
        weight(6) = d
        weight(7) = d

      else if ( rule .eq. 7 ) then

        a = 0.5D+00
        b = 1.0D+00
        c = 4.0D+00 / 24.0D+00
        d = 1.0D+00 / 24.0D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  b
        xtab(3) = -b
        xtab(4) =  z
        xtab(5) =  z
        xtab(6) =  a
        xtab(7) = -a
        xtab(8) = -a
        xtab(9) =  a

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  b
        ytab(5) = -b
        ytab(6) =  a
        ytab(7) =  a
        ytab(8) = -a
        ytab(9) = -a
 
        weight(1) = c
        weight(2) = d
        weight(3) = d
        weight(4) = d
        weight(5) = d
        weight(6) = c
        weight(7) = c
        weight(8) = c
        weight(9) = c

      else if ( rule .eq. 8 ) then

        a = sqrt ( 2.0D+00 ) / 2.0D+00
        b = sqrt ( 3.0D+00 / 5.0D+00 )
        c = sqrt ( 3.0D+00 / 10.0D+00 )

        w1 = 16.0D+00 / 72.0D+00
        w2 =  8.0D+00 / 72.0D+00
        w3 = 10.0D+00 / 72.0D+00
        w4 =  5.0D+00 / 72.0D+00

        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  a
        xtab(3) = -a
        xtab(4) =  z 
        xtab(5) =  z
        xtab(6) =  a
        xtab(7) =  a
        xtab(8) = -a
        xtab(9) = -a

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  b
        ytab(5) = -b
        ytab(6) =  c
        ytab(7) = -c
        ytab(8) =  c
        ytab(9) = -c

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w3
        weight(5) = w3
        weight(6) = w4
        weight(7) = w4
        weight(8) = w4
        weight(9) = w4

      else if ( rule .eq. 9 ) then

        a = sqrt ( 3.0D+00 ) / 2.0D+00
        b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &    / 104.0D+00 )
        c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &    / 104.0D+00 )
        u = 2.0D+00 / 27.0D+00
        v = ( 551.0D+00 + 41.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00
        w = ( 551.0D+00 - 41.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00
        z = 0.0D+00

        xtab(1) =  a
        xtab(2) = -a
        xtab(3) =  z
        xtab(4) =  z
        xtab(5) =  b
        xtab(6) = -b
        xtab(7) =  b
        xtab(8) = -b
        xtab(9) =  c
        xtab(10) = c
        xtab(11) = -c
        xtab(12) = -c
 
        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  a
        ytab(4) = -a
        ytab(5) =  b
        ytab(6) =  b
        ytab(7) = -b
        ytab(8) = -b
        ytab(9) =  c
        ytab(10) = -c
        ytab(11) =  c
        ytab(12) = -c

        weight(1) = u
        weight(2) = u
        weight(3) = u
        weight(4) = u
        weight(5) = v
        weight(6) = v
        weight(7) = v
        weight(8) = v
        weight(9) = w
        weight(10) = w
        weight(11) = w
        weight(12) = w

      else if ( rule .eq. 10 ) then

        a = sqrt ( ( 3.0D+00 - sqrt ( 5.0D+00 ) ) / 8.0D+00 )
        b = sqrt ( ( 15.0D+00 + 3.0D+00 * sqrt ( 5.0D+00 ) 
     &    - 2.0D+00 * sqrt ( 30.0D+00 ) 
     &    - 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
        c = sqrt ( ( 15.0D+00 + 3.0D+00 * sqrt ( 5.0D+00 ) 
     &    + 2.0D+00 * sqrt ( 30.0D+00 ) 
     &    + 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
        d = sqrt ( ( 3.0D+00 + sqrt ( 5.0D+00 ) ) / 8.0D+00 )
        e = sqrt ( ( 15.0D+00 - 3.0D+00 * sqrt ( 5.0D+00 ) 
     &    - 2.0D+00 * sqrt ( 30.0D+00 ) 
     &    + 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
        f = sqrt ( ( 15.0D+00 - 3.0D+00 * sqrt ( 5.0D+00 ) 
     &    + 2.0D+00 * sqrt ( 30.0D+00 ) 
     &    - 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
        w1 = ( 90.0D+00 + 5.0D+00 * sqrt ( 30.0D+00 ) 
     &     + 18.0D+00 * sqrt ( 5.0D+00 ) 
     &     + 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
        w2 = ( 90.0D+00 - 5.0D+00 * sqrt ( 30.0D+00 ) 
     &     + 18.0D+00 * sqrt ( 5.0D+00 ) 
     &     - 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
        w3 = ( 90.0D+00 + 5.0D+00 * sqrt ( 30.0D+00 ) 
     &     - 18.0D+00 * sqrt ( 5.0D+00 ) 
     &     - 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
        w4 = ( 90.0D+00 - 5.0D+00 * sqrt ( 30.0D+00 ) 
     &     - 18.0D+00 * sqrt ( 5.0D+00 ) 
     &     + 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00

        xtab(1) =  a
        xtab(2) =  a
        xtab(3) = -a
        xtab(4) = -a
        xtab(5) =  a
        xtab(6) =  a
        xtab(7) = -a
        xtab(8) = -a
        xtab(9) =  d
        xtab(10) =  d
        xtab(11) = -d
        xtab(12) = -d
        xtab(13) =  d
        xtab(14) =  d
        xtab(15) = -d
        xtab(16) = -d

        ytab(1) =  b
        ytab(2) = -b
        ytab(3) =  b
        ytab(4) = -b
        ytab(5) =  c
        ytab(6) = -c
        ytab(7) =  c
        ytab(8) = -c
        ytab(9) =  e
        ytab(10) = -e
        ytab(11) =  e
        ytab(12) = -e
        ytab(13) =  f
        ytab(14) = -f
        ytab(15) =  f
        ytab(16) = -f

        weight(1) = w1
        weight(2) = w1
        weight(3) = w1
        weight(4) = w1
        weight(5) = w2
        weight(6) = w2
        weight(7) = w2
        weight(8) = w2
        weight(9) = w3
        weight(10) = w3
        weight(11) = w3
        weight(12) = w3
        weight(13) = w4
        weight(14) = w4
        weight(15) = w4
        weight(16) = w4

      else if ( rule .eq. 11 ) then

        xtab(1) = 0.0D+00
        ytab(1) = 0.0D+00

        weight(1) = 1.0D+00 / 9.0D+00
        do i = 2, 11
          weight(i) = ( 16.0D+00 + sqrt ( 6.0D+00 ) ) / 360.0D+00
        end do
        do i = 12, 21
          weight(i) = ( 16.0D+00 - sqrt ( 6.0D+00 ) ) / 360.0D+00
        end do

        r = sqrt ( ( 6.0D+00 - sqrt ( 6.0D+00 ) ) / 10.0D+00 )

        do i = 1, 10
          a = 2.0D+00 * pi * real ( i, kind = 8 ) / 10.0D+00
          xtab(1+i) = r * cos ( a )
          ytab(1+i) = r * sin ( a )
        end do

        r = sqrt ( ( 6.0D+00 + sqrt ( 6.0D+00 ) ) / 10.0D+00 )

        do i = 1, 10
          a = 2.0D+00 * pi * real ( i, kind = 8 ) / 10.0D+00
          xtab(11+i) = r * cos ( a )
          ytab(11+i) = r * sin ( a )
        end do
c
c  There was apparently a misprint in the Lether paper.  The quantity
c  which here reads "322" was printed there as "332".
c
      else if ( rule .eq. 12 ) then

        a = 0.5D+00
        b = sqrt ( 3.0D+00 ) / 2.0D+00
        c = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 252.0D+00 )
        d = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 252.0D+00 )
        e = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 84.0D+00 )
        f = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 84.0D+00 )
        g = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 63.0D+00 )
        h = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) 
     &    / 63.0D+00 )

        w1 = 64.0D+00 / 675.0D+00
        w2 = 16.0D+00 / 225.0D+00
        w3 = 16.0D+00 / 675.0D+00
        w4 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 21600.0D+00
        w5 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 21600.0D+00
        w6 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 7200.0D+00
        w7 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 7200.0D+00
        w8 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 5400.0D+00
        w9 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 5400.0D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  a
        xtab(3) = -a
        xtab(4) =  b
        xtab(5) = -b
        xtab(6) =  b
        xtab(7) =  b
        xtab(8) = -b
        xtab(9) = -b
        xtab(10) =  b
        xtab(11) =  b
        xtab(12) = -b
        xtab(13) = -b
        xtab(14) =  a
        xtab(15) =  a
        xtab(16) = -a
        xtab(17) = -a
        xtab(18) =  a
        xtab(19) =  a
        xtab(20) = -a
        xtab(21) = -a
        xtab(22) =  z
        xtab(23) =  z
        xtab(24) =  z
        xtab(25) =  z

        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  z
        ytab(5) =  z
        ytab(6) =  c
        ytab(7) = -c
        ytab(8) =  c
        ytab(9) = -c
        ytab(10) =  d
        ytab(11) = -d
        ytab(12) =  d
        ytab(13) = -d
        ytab(14) =  e
        ytab(15) = -e
        ytab(16) =  e
        ytab(17) = -e
        ytab(18) =  f
        ytab(19) = -f
        ytab(20) =  f
        ytab(21) = -f
        ytab(22) =  g
        ytab(23) = -g
        ytab(24) =  h
        ytab(25) = -h

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w3
        weight(5) = w3
        weight(6) = w4
        weight(7) = w4
        weight(8) = w4
        weight(9) = w4
        weight(10) = w5
        weight(11) = w5
        weight(12) = w5
        weight(13) = w5
        weight(14) = w6
        weight(15) = w6
        weight(16) = w6
        weight(17) = w6
        weight(18) = w7
        weight(19) = w7
        weight(20) = w7
        weight(21) = w7
        weight(22) = w8
        weight(23) = w8
        weight(24) = w9
        weight(25) = w9

      else if ( rule .eq. 13 ) then

        nr = 4
        call legendre_set ( nr, ra, rw )
        a = -1.0D+00
        b = +1.0D+00
        c =  0.0D+00
        d = +1.0D+00
        call rule_adjust ( a, b, c, d, nr, ra, rw )
        do i = 1, nr
          ra(i) = sqrt ( ra(i) )
        end do
        i = 0

        do j = 1, 16

          c = cos ( pi * dble ( j ) / 8.0D+00 )
          s = sin ( pi * dble ( j ) / 8.0D+00 )

          do k = 1, 4

            i = i + 1
            xtab(i) = c * ra(k)
            ytab(i) = s * ra(k)
            weight(i) = rw(k) / 16.0D+00

          end do

        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CIRCLE_XY_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  There is no rule of index ', rule
        stop

      end if

      return
      end
      subroutine circle_xy_size ( rule, order )

c*********************************************************************72
c
cc CIRCLE_XY_SIZE sizes an XY quadrature rule inside the unit circle in 2D.
c
c  Integration region:
c
c    X*X + Y*Y <= 1.0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Frank Lether,
c    A Generalized Product Rule for the Circle,
c    SIAM Journal on Numerical Analysis,
c    Volume 8, Number 2, June 1971, pages 249-253.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c      1, 1 point 1-st degree;
c      2, 4 point 3-rd degree, Stroud S2:3-1;
c      3, 4 point 3-rd degree, Lether #1;
c      4, 4 point 3-rd degree, Stroud S2:3-2;
c      5, 5 point 3-rd degree;
c      6, 7 point 5-th degree;
c      7, 9 point 5-th degree;
c      8, 9 point 5-th degree, Lether #2;
c      9, 12 point 7-th degree;
c     10, 16 point 7-th degree, Lether #3;
c     11, 21 point 9-th degree, Stroud S2:9-3;
c     12, 25 point 9-th degree, Lether #4 (after correcting error);
c     13, 64 point 15-th degree Gauss product rule.
c
c    Output, integer ORDER, the order of the desired rule.
c
      implicit none

      integer order
      integer rule

      if ( rule == 1 ) then

        order = 1

      else if ( rule == 2 ) then

        order = 4

      else if ( rule == 3 ) then

        order = 4

      else if ( rule == 4 ) then

        order = 4

      else if ( rule == 5 ) then

        order = 5

      else if ( rule == 6 ) then

        order = 7

      else if ( rule == 7 ) then

        order = 9

      else if ( rule == 8 ) then

        order = 9

      else if ( rule == 9 ) then

        order = 12

      else if ( rule == 10 ) then

        order = 16

      else if ( rule == 11 ) then

        order = 21

      else if ( rule == 12 ) then

        order = 25

      else if ( rule == 13 ) then

        order = 64

      else

        order = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CIRCLE_XY_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  There is no rule of index ', rule
        stop

      end if

      return
      end
      subroutine circle_xy_sum ( func, center, r, order, xtab, ytab, 
     &  weight, result )

c*********************************************************************72
c
cc CIRCLE_XY_SUM applies an XY quadrature rule inside a circle in 2D.
c
c  Integration region:
c
c    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= R * R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of two variables which is to be integrated,
c    of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the coordinates of the center of 
c    the circle.
c
c    Input, double precision R, the radius of the circle.
c
c    Input, integer ORDER, the order of the rule.  The rule is
c    assumed to be defined on the unit circle.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the XY
c    coordinates of the abscissas of the quadrature rule for the unit circle.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer order

      double precision center(dim_num)
      double precision circle_area_2d
      double precision func
      external func
      integer i
      double precision quad
      double precision r
      double precision result
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y
      double precision ytab(order)

      quad = 0.0D+00
      do i = 1, order
        x = center(1) + r * xtab(i)
        y = center(2) + r * ytab(i)
        quad = quad + weight(i) * func ( x, y )
      end do

      volume = circle_area_2d ( r )
      result = quad * volume

      return
      end
      subroutine cn_geg_00_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc CN_GEG_00_1 implements the midpoint rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call c1_geg_monomial_integral ( alpha, expon, volume )
      volume = volume ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      w(k) = volume

      return
      end
      subroutine cn_geg_00_1_size ( n, alpha, o )

c*********************************************************************72
c
cc CN_GEG_00_1_SIZE sizes the midpoint rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_00_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine cn_geg_01_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision value1
      double precision value2
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha <= -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_01_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call c1_geg_monomial_integral ( alpha, expon, value1 )
      volume = value1 ** n

      expon = 1
      call c1_geg_monomial_integral ( alpha, expon, value2 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = value2 / value1
      end do
      w(k) = volume

      return
      end
      subroutine cn_geg_01_1_size ( n, alpha, o )

c*********************************************************************72
c
cc CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_01_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine cn_geg_02_xiu ( n, alpha, o, x, w )

c*********************************************************************72
c
cc CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision arg
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_02_XIU - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = 1.0D+00
      delta0 = 0.0D+00
      c1 = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      expon = 0
      call c1_geg_monomial_integral ( alpha, expon, volume_1d )
      volume = volume_1d ** n

      w(1:o) = volume / dble ( o )

      return
      end
      subroutine cn_geg_02_xiu_size ( n, alpha, o )

c*********************************************************************72
c
cc CN_GEG_02_XIU_SIZE sizes the Xiu rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_02_XIU_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = n + 1

      return
      end
      subroutine cn_geg_03_xiu ( n, alpha, o, x, w )

c*********************************************************************72
c
cc CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 January 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      double precision alpha
      double precision arg
      integer expon
      integer i
      integer j
      integer n
      integer o
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_03_XIU - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call c1_geg_monomial_integral ( alpha, expon, volume )
      volume = volume ** n

      do j = 1, o

        i = 0 
        do r = 1, n / 2
          arg = dble ( ( 2 * r - 1 ) * j ) * pi / dble ( n )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
     &      / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) 
     &      / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
        end do

        if ( i < n ) then
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) 
     &      / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
          if ( n .eq. 1 ) then
            x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
          end if
        end if

      end do

      w(1:o) = volume / dble ( o )

      return
      end
      subroutine cn_geg_03_xiu_size ( n, alpha, o )

c*********************************************************************72
c
cc CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the parameter.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_03_XIU_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 2 * n

      return
      end
      subroutine cn_geg_monomial_integral ( n, alpha, expon, value )

c*********************************************************************72
c
cc CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
c
c  Discussion:
c
c    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
c
c      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
c
c    with -1.0 < alpha.
c
c    value = integral ( CN ) 
c      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of (1-X).
c    -1.0 < ALPHA.
c
c    Input, integer EXPON(N), the exponents.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      double precision alpha
      integer expon(n)
      integer i
      double precision value
      double precision value2

      if ( alpha <= -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_GEG_MONOMIAL_INTEGRAL - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      value = 1.0D+00
      do i = 1, n
        call c1_geg_monomial_integral ( alpha, expon(i), value2 )
        value = value * value2
      end do

      return
      end
      subroutine cn_jac_00_1 ( n, alpha, beta, o, x, w )

c*********************************************************************72
c
cc CN_JAC_00_1 implements the midpoint rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha.
c
c    with -1 < alpha, -1 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision beta
      integer expon
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision volume
      double precision :: w(o)
      double precision :: x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      expon = 0
      call c1_jac_monomial_integral ( alpha, beta, expon, volume )
      volume = volume ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = volume

      return
      end
      subroutine cn_jac_00_1_size ( n, alpha, beta, o )

c*********************************************************************72
c
cc CN_JAC_00_1_SIZE sizes the midpoint rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha.
c
c    with -1 < alpha, -1 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      double precision beta
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine cn_jac_01_1 ( n, alpha, beta, o, x, w )

c*********************************************************************72
c
cc CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha. 
c
c    with -1 < alpha, -1 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision beta
      integer expon
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision value1
      double precision value2
      double precision volume
      double precision :: w(o)
      double precision :: x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      expon = 0
      call c1_jac_monomial_integral ( alpha, beta, expon, value1 )
      volume = value1 ** n

      expon = 1
      call c1_jac_monomial_integral ( alpha, beta, expon, value2 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = value2 / value1
      end do
      w(k) = volume

      return
      end
      subroutine cn_jac_01_1_size ( n, alpha, beta, o )

c*********************************************************************72
c
cc CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha. 
c
c    with -1 < alpha, -1 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      double precision beta
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine cn_jac_02_xiu ( n, alpha, beta, o, x, w )

c*********************************************************************72
c
cc CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha.
c
c    with -1.0 < alpha, -1.0 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision arg
      double precision beta
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = ( alpha + beta + 2.0D+00 ) / 2.0D+00
      delta0 = ( alpha - beta ) / 2.0D+00
      c1 = 2.0D+00 * ( alpha + 1.0D+00 ) * ( beta + 1.0D+00 ) 
     &  / ( alpha + beta + 3.0D+00 ) / ( alpha + beta + 2.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) 
     &      / gamma0
        end do
      end do

      expon = 0
      call c1_jac_monomial_integral ( alpha, beta, expon, volume_1d )
      volume = volume_1d ** n

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine cn_jac_02_xiu_size ( n, alpha, beta, o )

c*********************************************************************72
c
cc CN_JAC_02_XIU_SIZE sizes the Xiu rule for region CN_JAC.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
c
c      w(alpha,beta;x) = product ( 1 .le. i .le. n ) (1-x(i))^beta (1+x(i))^alpha.
c
c    with -1 < alpha, -1 < beta.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, BETA, the parameters.
c    -1.0 < ALPHA, -1.0 < BETA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      double precision beta
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
        write ( *, '(a)' ) '  ALPHA .le. -1.0'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
        write ( *, '(a)' ) '  BETA .le. -1.0'
        stop
      end if

      o = n + 1

      return
      end
      subroutine cn_jac_monomial_integral ( n, alpha, beta, expon, 
     &  value )

c*********************************************************************72
c
cc CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
c
c  Discussion:
c
c    value = integral ( CN ) 
c      product ( 1 .le. i .le. n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of (1-X) in the weight factor.
c
c    Input, double precision BETA, the exponent of (1+X) in the weight factor.
c
c    Input, integer EXPON(N), the exponents.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      double precision alpha
      double precision beta
      integer expon(n)
      integer i
      double precision value
      double precision value2

      value = 1.0D+00
      do i = 1, n
        call c1_jac_monomial_integral ( alpha, beta, expon(i), value2 )
        value = value * value2
      end do

      return
      end
      subroutine cn_leg_01_1 ( n, o, x, w )

c*********************************************************************72
c
cc CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      integer expon
      integer i
      integer j
      integer k
      double precision value1
      double precision volume
      double precision w(o)
      double precision x(n,o)

      expon = 0
      call c1_leg_monomial_integral ( expon, value1 )
      volume = value1 ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      w(k) = volume

      return
      end
      subroutine cn_leg_01_1_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 1

      return
      end
      subroutine cn_leg_02_xiu ( n, o, x, w )

c*********************************************************************72
c
cc CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision arg
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i < n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = 1.0D+00
      delta0 = 0.0D+00
      c1 = 1.0D+00 / 3.0D+00

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      expon = 0
      call c1_leg_monomial_integral ( expon, volume_1d )
      volume = volume_1d ** n

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine cn_leg_02_xiu_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_02_XIU_SIZE sizes the Xiu rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n + 1

      return
      end
      subroutine cn_leg_03_1 ( n, o, x, w )

c*********************************************************************72
c
cc CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c    The necessary treatment of the final coordinate of points when
c    N is odd seems to vary from what Stroud declares! 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit  none

      integer n
      integer o

      double precision arg
      integer expon
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision w(o)
      double precision x(n,o)

      expon = 0
      call c1_leg_monomial_integral ( expon, volume )
      volume = volume ** n

      do j = 1, o

        i = 0

        do r = 1, ( n / 2 )
          arg = dble ( ( 2 * r - 1 ) * j ) * pi / dble ( n )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )
        end do
c
c  The following code does not correspond to what Stroud declares.
c
        if ( i .lt. n ) then

          i = i + 1
          if ( n .eq. 1 ) then
            x(i,j) =                r8_mop ( j ) / sqrt ( 3.0D+00 )
          else
            x(i,j) = sqrt ( 2.0 ) * r8_mop ( j ) / sqrt ( 3.0D+00 )
          end if
        end if

      end do

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine cn_leg_03_1_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n

      return
      end
      subroutine cn_leg_03_xiu ( n, o, x, w )

c*********************************************************************72
c
cc CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision arg
      integer expon
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision w(o)
      double precision x(n,o)

      expon = 0
      call c1_leg_monomial_integral ( expon, volume )
      volume = volume ** n

      do j = 1, o

        i = 0 
        do r = 1, n / 2
          arg = dble ( ( 2 * r - 1 ) * j ) * pi / dble ( n )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )
        end do

        if ( i < n ) then
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) / sqrt ( 3.0D+00 )
          if ( n .eq. 1 ) then
            x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
          end if
        end if

      end do

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine cn_leg_03_xiu_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n

      return
      end
      subroutine cn_leg_05_1 ( n, option, o, x, w )

c*********************************************************************72
c
cc CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N^2 + N + 2.
c
c    The rule has precision P = 5.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    N must be 4, 5, or 6.
c
c    Input, integer OPTION, is only used in case N = 4 or 5.
c    In that case, OPTION should be 1 or 2 to select the
c    two available variants of the rule.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision arg
      double precision b
      double precision c
      double precision eta
      integer expon
      double precision gamma
      integer i
      integer i1
      integer i2
      integer k
      double precision lambda
      double precision mu
      integer option
      double precision volume
      double precision w(o)
      double precision x(n,o)
      double precision xsi

      if ( n .lt. 4 .or. 6 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
        write ( *, '(a)' ) '  The value of N must be 4, 5, or 6.'
        stop
      end if

      if ( n .eq. 4 .or. n .eq. 5 ) then
        if ( option .lt. 1 .or. 2 .lt. option ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
          write ( *, '(a)' ) '  When N = 4 or 5, OPTION must be 1 or 2.'
          stop
        end if
      end if

      expon = 0
      call c1_leg_monomial_integral ( expon, volume )
      volume = volume ** n

      if ( n .eq. 4 .and. option .eq. 1 ) then

        eta    =   0.778984505799815D+00
        lambda =   1.284565137874656D+00
        xsi =    - 0.713647298819253D+00
        mu =     - 0.715669761974162D+00
        gamma =    0.217089151000943D+00
        a =        0.206186096875899D-01 * volume
        b =        0.975705820221664D-02 * volume
        c =        0.733921929172573D-01 * volume

      else if ( n .eq. 4 .and. option .eq. 2 ) then

        eta    =   0.546190755827425D+00
        lambda =   0.745069130115661D+00
        xsi =    - 0.413927294508700D+00
        mu =     - 0.343989637454535D+00
        gamma =    1.134017894600344D+00
        a =        0.853094758323323D-01 * volume
        b =        0.862099000096395D-01 * volume
        c =        0.116418206881849D-01 * volume

      else if ( n .eq. 5 .and. option .eq. 1 ) then

        eta    =   0.522478547481276D+00
        lambda =   0.936135175985774D+00
        xsi =    - 0.246351362101519D+00
        mu =     - 0.496308106093758D+00
        gamma =    0.827180176822930D+00
        a =        0.631976901960153D-01 * volume
        b =        0.511464127430166D-01 * volume
        c =        0.181070246088902D-01 * volume

      else if ( n .eq. 5 .and. option .eq. 2 ) then

        eta    =   0.798317301388741D+00
        lambda =   0.637344273885728D+00
        xsi =    - 0.455245909918377D+00
        mu =     - 1.063446229997311D+00
        gamma =    0.354482076665770D+00
        a =        0.116952384292206D-01 * volume
        b =        0.701731258612708D-01 * volume
        c =        0.137439132264426D-01 * volume

      else if ( n .eq. 6 ) then

        eta    =   0.660225291773525D+00
        lambda =   1.064581294844754D+00
        xsi =      0.000000000000000D+00
        mu =     - 0.660225291773525D+00
        gamma =    0.660225291773525D+00
        a =        0.182742214532872D-01 * volume
        b =        0.346020761245675D-01 * volume
        c =        0.182742214532872D-01 * volume

      end if

      k = 0

      k = k + 1
      do i = 1, n
        x(i,k) = eta
      end do
      w(k) = a

      k = k + 1
      do i = 1, n
        x(i,k) = - eta
      end do
      w(k) = a

      do i1 = 1, n
        k = k + 1
        do i = 1, n
          x(i,k) = xsi
        end do
        x(i1,k) = lambda
        w(k) = b
      end do

      do i1 = 1, n
        k = k + 1
        do i = 1, n
          x(i,k) = - xsi
        end do
        x(i1,k) = - lambda
        w(k) = b
      end do

      do i1 = 1, n - 1
        do i2 = i1 + 1, n
          k = k + 1
          do i = 1, n
            x(i,k) = gamma
          end do
          x(i1,k) = mu
          x(i2,k) = mu
          w(k) = c
        end do
      end do

      do i1 = 1, n - 1
        do i2 = i1 + 1, n
          k = k + 1
          do i = 1, n
            x(i,k) = - gamma
          end do
          x(i1,k) = - mu
          x(i2,k) = - mu
          w(k) = c
        end do
      end do

      return
      end
      subroutine cn_leg_05_1_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N^2 + N + 2.
c
c    The rule has precision P = 5.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n * n + n + 2

      return o
      end
      subroutine cn_leg_05_2 ( n, o, x, w )

c*********************************************************************72
c
cc CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 N^2 + 1.
c
c    The rule has precision P = 5.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    N must be at least 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision b0
      double precision b1
      double precision b2
      integer expon
      integer i
      integer i1
      integer i2
      integer k
      double precision r
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CN_LEG_05_2 - Fatal error!'
        write ( *, '(a)' ) '  N must be at least 2.'
        stop
      end if

      expon = 0
      call c1_leg_monomial_integral ( expon, volume )
      volume = volume ** n

      b0 = dble ( 25 * n * n - 115 * n + 162 ) * volume / 162.0D+00
      b1 = dble ( 70 - 25 * n ) * volume / 162.0D+00
      b2 = 25.0D+00 * volume / 324.0D+00

      r = sqrt ( 3.0D+00 / 5.0D+00 )

      k = 0

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      w(k) = b0

      do i1 = 1, n

        k = k + 1
        do i = 1, n
          x(i,k) = 0.0D+00
        end do
        x(i1,k) = + r
        w(k) = b1

        k = k + 1
        do i = 1, n
          x(i,k) = 0.0D+00
        end do
        x(i1,k) = - r
        w(k) = b1

      end do

      do i1 = 1, n - 1
        do i2 = i1 + 1, n

          k = k + 1
          do i = 1, n
            x(i,k) = 0.0D+00
          end do
          x(i1,k) = + r
          x(i2,k) = + r
          w(k) = b2

          k = k + 1
          do i = 1, n
            x(i,k) = 0.0D+00
          end do
          x(i1,k) = + r
          x(i2,k) = - r
          w(k) = b2

          k = k + 1
          do i = 1, n
            x(i,k) = 0.0D+00
          end do
          x(i1,k) = - r
          x(i2,k) = + r
          w(k) = b2

          k = k + 1
          do i = 1, n
            x(i,k) = 0.0D+00
          end do
          x(i1,k) = - r
          x(i2,k) = - r
          w(k) = b2

        end do
      end do

      return
      end
      subroutine cn_leg_05_2_size ( n, o )

c*********************************************************************72
c
cc CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 N^2 + 1.
c
c    The rule has precision P = 5.
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n * n + 1

      return
      end
      subroutine cn_leg_monomial_integral ( n, expon, value )

c*********************************************************************72
c
cc CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
c
c  Discussion:
c
c    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
c
c      w(x) = 1.
c
c    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer EXPON(N), the exponents.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      integer expon(n)
      integer i
      double precision value
      double precision value2

      value = 1.0D+00
      do i = 1, n
        call c1_leg_monomial_integral ( expon(i), value2 )
        value = value * value2
      end do

      return
      end
      subroutine cone_unit_3d ( func, result )

c*********************************************************************72
c
cc CONE_UNIT_3D approximates an integral inside the unit cone in 3D.
c
c  Integration Region:
c
c      X*X + Y*Y <= 1 - Z  
c
c    and
c
c      0 <= Z <= 1.
c
c  Discussion:
c
c    An 48 point degree 7 formula, Stroud CN:S2:7-1, is used.
c
c    (There is a typographical error in the S2:7-1 formula for B3.)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision cone_volume_3d
      double precision func
      external func
      double precision h
      integer i
      double precision quad
      double precision r
      double precision result
      double precision u(4)
      double precision volume
      double precision w1(4)
      double precision w2(3)
      double precision x
      double precision y
      double precision z

      save u
      save w1

      data u /
     &  0.04850054945D+00, 0.2386007376D+00, 
     &  0.5170472951D+00,  0.7958514179D+00 /
      data w1 /
     &  0.1108884156D+00,  0.1434587878D+00, 
     &  0.06863388717D+00, 0.01035224075D+00 /

      a = sqrt ( 3.0D+00 ) / 2.0D+00
      b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &  / 104.0D+00 )
      c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) 
     &  / 104.0D+00 )
      w2(1:3) = 3.0D+00 * 2.0D+00 / 27.0D+00
      w2(2) = 3.0D+00 * ( 551.0D+00 + 4.0D+00 * sqrt ( 29.0D+00 ) ) 
     &  / 6264.0D+00
      w2(3) = 3.0D+00 * ( 551.0D+00 - 4.0D+00 * sqrt ( 29.0D+00 ) ) 
     &  / 6264.0D+00

      quad = 0.0D+00

      do i = 1, 4

        x = a * ( 1.0D+00 - u(i) )
        y = 0.0D+00
        z = u(i)
        quad = quad + w1(i) * w2(1) * func ( x, y, z )

        x = -a * ( 1.0D+00 - u(i) )
        y = 0.0D+00
        z = u(i)
        quad = quad + w1(i) * w2(1) * func ( x, y, z )

        x = 0.0D+00
        y = a * ( 1.0D+00 - u(i) )
        z = u(i)
        quad = quad + w1(i) * w2(1) * func ( x, y, z )

        x = 0.0D+00
        y = -a * ( 1.0D+00 - u(i) )
        z = u(i)
        quad = quad + w1(i) * w2(1) * func ( x, y, z )

      end do

      do i = 1, 4

        x =  b * ( 1.0D+00 - u(i) )
        y =  b * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(2) * func ( x, y, z )

        x = -b * ( 1.0D+00 - u(i) )
        y =  b * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(2) * func ( x, y, z )

        x = -b * ( 1.0D+00 - u(i) )
        y = -b * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(2) * func ( x, y, z )

        x =  b * ( 1.0D+00 - u(i) )
        y = -b * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(2) * func ( x, y, z )

        x =  c * ( 1.0D+00 - u(i) )
        y =  c * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(3) * func ( x, y, z )

        x = -c * ( 1.0D+00 - u(i) )
        y =  c * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(3) * func ( x, y, z )

        x = -c * ( 1.0D+00 - u(i) )
        y = -c * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(3) * func ( x, y, z )

        x =  c * ( 1.0D+00 - u(i) )
        y = -c * ( 1.0D+00 - u(i) )
        z =  u(i)
        quad = quad + w1(i) * w2(3) * func ( x, y, z )

      end do

      r = 1.0D+00
      h = 1.0D+00

      volume = cone_volume_3d ( r, h )
      result = quad * volume

      return
      end
      function cone_volume_3d ( r, h )

c*********************************************************************72
c
cc CONE_VOLUME_3D returns the volume of a cone in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the base of the cone.
c
c    Input, double precision H, the height of the cone.
c
c    Output, double precision CONE_VOLUME_3D, the volume of the cone.
c
      implicit none

      double precision cone_volume_3d
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      cone_volume_3d = ( pi / 3.0D+00 ) * h * r * r

      return
      end
      subroutine cube_shell_nd ( func, n, r1, r2, result )

c*********************************************************************72
c
cc CUBE_SHELL_ND approximates an integral inside a cubic shell in N dimensions.
c
c  Integration region:
c
c    R1 <= abs ( X(1:N) ) <= R2
c
c  Discussion:
c
c    An N*2^N point third degree formula is used, Stroud number CNSHELL:3-4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision R1, R2, the inner and outer radii of the cubical
c    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
c    2*R1.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision cube_shell_volume_nd
      logical done
      double precision func
      external func
      integer i
      integer j
      double precision quad
      double precision r1
      double precision r2
      double precision rmax
      double precision rmin
      double precision result
      double precision u
      double precision v
      double precision volume
      double precision x(n)

      if ( r1 .eq. r2 ) then
        result = 0.0D+00
        return
      end if

      rmax = max ( r1, r2 )
      rmin = min ( r1, r2 )    

      u = sqrt ( dble ( n ) * ( rmax**(n+2) - rmin**(n+2) ) 
     &  / ( dble ( n + 2 ) * ( rmax**n - rmin**n ) ) )
      v = u / sqrt ( 3.0D+00 )

      quad = 0.0D+00

      do i = 1, n

        do j = 1, n
          x(j) = v
        end do
        x(i) = u

        do

          quad = quad + func ( n, x )

          call r8vec_mirror_next ( n, x, done )

          if ( done ) then
            exit
          end if

        end do

      end do

      quad = quad / dble ( n * 2**n )

      volume = cube_shell_volume_nd ( n, r1, r2 )
      result = quad * volume
      
      return
      end
      function cube_shell_volume_nd ( n, r1, r2 )

c*********************************************************************72
c
cc CUBE_SHELL_VOLUME_ND computes the volume of a cubic shell in ND.
c
c  Integration region:
c
c    R1 <= abs ( X(1:N) ) <= R2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision R1, R2, the inner and outer radii of the cubic
c    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
c    2*R1.
c
c    Output, double precision CUBE_SHELL_VOLUME_ND, the volume of the cubic
c    shell.
c
      implicit none

      double precision cube_shell_volume_nd
      integer n
      double precision r1
      double precision r2

      cube_shell_volume_nd = ( r2**n - r1**n ) * 2**n

      return
      end
      subroutine cube_unit_3d ( func, result )

c*********************************************************************72
c
cc CUBE_UNIT_3D approximates an integral inside the unit cube in 3D.
c
c  Integration region:
c
c      -1 <= X <= 1,
c    and
c      -1 <= Y <= 1,
c    and
c      -1 <= Z <= 1.
c
c  Discussion:
c
c    An 8 point third degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision cube_unit_volume_nd
      double precision func
      external func
      double precision quad
      double precision result
      double precision s
      double precision volume
      double precision w
      double precision x
      double precision y
      double precision z

      s = 1.0D+00 / sqrt ( 3.0D+00 )
      w = 1.0D+00 / 8.0D+00

      x = s
      y = s
      z = s

      quad = w * ( 
     &    func (  x,  y,  z ) + func (  x,  y, -z ) 
     &  + func (  x, -y,  z ) + func (  x, -y, -z ) 
     &  + func ( -x,  y,  z ) + func ( -x,  y, -z ) 
     &  + func ( -x, -y,  z ) + func ( -x, -y, -z ) )

      volume = cube_unit_volume_nd ( 3 )
      result = quad * volume

      return
      end
      subroutine cube_unit_nd ( func, qa, qb, n, k )

c*********************************************************************72
c
cc CUBE_UNIT_ND approximates an integral inside the unit cube in ND.
c
c  Integration region:
c
c    -1 <= X(1:N) <= 1
c
c  Discussion:
c
c    A K^N point product formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Lyness, BJJ McHugh,
c    Integration Over Multidimensional Hypercubes, 
c    A Progressive Procedure,
c    The Computer Journal,
c    Volume 6, 1963, pages 264-270.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates the function, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Output, double precision QA(K), QB(K), two sets of estimates for
c    the integral.  The QB entries are obtained from the
c    QA entries by Richardson extrapolation, and QB(K) is
c    the best estimate for the integral.
c
c    Input, integer N, the dimension of the cube.
c
c    Input, integer K, the highest order of integration, and the 
c    order of Richardson extrapolation.  K can be no greater than 10.
c
      implicit none

      integer kmax
      parameter ( kmax = 10 )

      integer k
      integer n

      double precision func
      external func
      double precision g(kmax,kmax)
      integer i
      integer j
      double precision qa(k)
      double precision qb(k)

      do j = 1, kmax
        do i = 1, kmax
          g(i,j) = 0.0D+00
        end do
      end do

      g( 1, 1) =  1.0D+00
      g( 2, 1) = -0.3333333333333D+00
      g( 2, 2) =  0.1333333333333D+01
      g( 3, 1) =  0.4166666666667D-01
      g( 3, 2) = -0.1066666666667D+01
      g( 3, 3) =  0.2025000000000D+01
      g( 4, 1) = -0.2777777777778D-02
      g( 4, 2) =  0.3555555555556D+00
      g( 4, 3) = -0.2603571428571D+01
      g( 4, 4) =  0.3250793650794D+01
      g( 5, 1) =  0.1157407407407D-03
      g( 5, 2) = -0.6772486772487D-01
      g( 5, 3) =  0.1464508928571D+01
      g( 5, 4) = -0.5779188712522D+01
      g( 5, 5) =  0.5382288910935D+01
      g( 6, 1) = -0.3306878306878D-05
      g( 6, 2) =  0.8465608465608D-02
      g( 6, 3) = -0.4881696428571D+00
      g( 6, 4) =  0.4623350970018D+01
      g( 6, 5) = -0.1223247479758D+02
      g( 6, 6) =  0.9088831168831D+01
      g( 7, 1) =  0.6889329805996D-07
      g( 7, 2) = -0.7524985302763D-03
      g( 7, 3) =  0.1098381696429D+00
      g( 7, 4) = -0.2241624712736D+01
      g( 7, 5) =  0.1274216124748D+02
      g( 7, 6) = -0.2516907092907D+02
      g( 7, 7) =  0.1555944865432D+02
      g( 8, 1) = -0.1093544413650D-08
      g( 8, 2) =  0.5016656868509D-04
      g( 8, 3) = -0.1797351866883D-01
      g( 8, 4) =  0.7472082375786D+00
      g( 8, 5) = -0.8168052081717D+01
      g( 8, 6) =  0.3236023405166D+02
      g( 8, 7) = -0.5082753227079D+02
      g( 8, 8) =  0.2690606541646D+02
      g( 9, 1) =  0.1366930517063D-10
      g( 9, 2) = -0.2606055516108D-05
      g( 9, 3) =  0.2246689833604D-02
      g( 9, 4) = -0.1839281815578D+00
      g( 9, 5) =  0.3646451822195D+01
      g( 9, 6) = -0.2588818724133D+02
      g( 9, 7) =  0.7782965878964D+02
      g( 9, 8) = -0.1012934227443D+03
      g( 9, 9) =  0.4688718347156D+02
      g(10, 1) = -0.1380737896023D-12
      g(10, 2) =  0.1085856465045D-06
      g(10, 3) = -0.2222000934334D-03
      g(10, 4) =  0.3503393934435D-01
      g(10, 5) = -0.1215483940732D+01
      g(10, 6) =  0.1456210532325D+02
      g(10, 7) = -0.7477751530769D+02
      g(10, 8) =  0.1800771959898D+03
      g(10, 9) = -0.1998874663788D+03
      g(10,10) =  0.8220635246624D+02

      if ( kmax .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CUBE_UNIT_ND - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  K must be no greater than KMAX = ', kmax
        write ( *, '(a,i8)' ) '  but the input K is ', k
        stop
      end if

      do i = 1, k
        call qmdpt ( func, n, i, qa(i) )
      end do

      qb(1) = qa(1)

      do i = 2, k
        qb(i) = 0.0D+00
        do j = 1, i
          qb(i) = qb(i) + g(i,j) * qa(j)
        end do
      end do

      return
      end
      function cube_unit_volume_nd ( n )

c*********************************************************************72
c
cc CUBE_UNIT_VOLUME_ND returns the volume of the unit cube in ND.
c
c  Integration region:
c
c    -1 <= X(1:N) <= 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision CUBE_UNIT_VOLUME_ND, the volume of the unit
c    cube in ND.
c
      implicit none

      double precision cube_unit_volume_nd
      integer n

      cube_unit_volume_nd = 2.0D+00**n

      return
      end
      function ellipse_area_2d ( r1, r2 )

c*********************************************************************72
c
cc ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
c
c  Integration region:
c
c    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the major and minor semi-axes.
c
c    Output, double precision ELLIPSE_AREA_2D, the area of the ellipse.
c
      implicit none

      double precision ellipse_area_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2

      ellipse_area_2d = pi * r1 * r2

      return
      end
      function ellipse_circumference_2d ( r1, r2 )

c*********************************************************************72
c
cc ELLIPSE_CIRCUMFERENCE_2D returns the circumference of an ellipse in 2D.
c
c  Discussion:
c
c    There is no closed formula for the circumference of an ellipse.
c
c    Defining the eccentricity by
c
c      E = sqrt ( 1 - ( r2 / r1 )^2 )
c
c    where R1 and R2 are the major and minor axes, then
c
c      circumference
c        = 4 * R1 * E(K,2*PI)
c        = R1 * Integral ( 0 <= T <= 2*PI ) sqrt ( 1 - E * E * sin^2 ( T ) ) dT
c
c    This integral can be approximated by the Gauss-Kummer formula.
c
c  Integration region:
c
c    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Harris, Horst Stocker,
c    Handbook of Mathematics and Computational Science,
c    Springer, 1998,
c    ISBN: 0-387-94746-9,
c    LC: QA40.S76.
c
c  Parameters:
c
c    Input, double precision R1, R2, the major and minor semi-axes.
c
c    Output, double precision ELLIPSE_CIRCUMFERENCE_2D, the
c    circumference of the ellipse.
c
      implicit none

      double precision ellipse_circumference_2d
      double precision e
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_epsilon
      double precision term
      double precision test
      double precision value

      if ( r1 .eq. r2 ) then
        ellipse_circumference_2d = 2.0D+00 * pi * r1
        return
      end if

      test = r8_epsilon ( )
c
c  Compute the eccentricity of the ellipse.
c
      e = sqrt ( 1.0D+00 - ( min ( r1, r2 ) / max ( r1, r2 ) )**2 )

      value = 1.0D+00
      term = value
      i = 0

10    continue

        i = i + 1
        term = term * ( 2 * i - 3 ) * ( 2 * i - 1 ) * e * e 
     &    / dble ( 2 * 2 * i * i )

        if ( abs ( term ) .le. test * ( abs ( value ) + 1.0D+00 ) ) then
          go to 20
        end if

        value = value + term

      go to 10

20    continue

      ellipse_circumference_2d = 2.0D+00 * pi * max ( r1, r2 ) * value

      return
      end
      function ellipse_eccentricity_2d ( r1, r2 )

c*********************************************************************72
c
cc ELLIPSE_ECCENTRICITY_2D returns the eccentricity of an ellipse in 2D.
c
c  Integration region:
c
c    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the major and minor semi-axes.
c
c    Output, double precision ELLIPSE_ECCENTRICITY_2D, the eccentricity 
c    of the ellipse.
c
      implicit none

      double precision ellipse_eccentricity_2d
      double precision major
      double precision minor
      double precision r1
      double precision r2

      minor = min ( abs ( r1 ), abs ( r2 ) )
      major = max ( abs ( r1 ), abs ( r2 ) )

      if ( major .eq. 0.0D+00 ) then
        ellipse_eccentricity_2d = - huge ( r1 )
        return
      end if

      ellipse_eccentricity_2d = sqrt ( 1.0D+00 - ( minor / major )**2 )

      return
      end
      function ellipsoid_volume_3d ( r1, r2, r3 )

c*********************************************************************72
c
cc ELLIPSOID_VOLUME_3D returns the volume of an ellipsoid in 3d.
c
c  Discussion:
c
c    This is not a general ellipsoid, but one for which each of the 
c    axes lies along a coordinate axis.
c
c  Integration region:
c
c      ( ( X - CENTER(1) ) / R1 )^2 
c    + ( ( Y - CENTER(2) ) / R2 )^2
c    + ( ( Z - CENTER(3) ) / R3 )^2 <= 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, R3, the semi-axes of the ellipsoid.
c
c    Output, double precision ELLIPSOID_VOLUME_3D, the volume of the ellipsoid.
c
      implicit none

      double precision ellipsoid_volume_3d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r3

      ellipsoid_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi * r1 * r2 * r3

      return
      end
      subroutine en_r2_01_1 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_01_1 implements the Stroud rule 1.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi**n )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      w(k) = volume

      return
      end
      subroutine en_r2_01_1_size ( n, o )

c*********************************************************************72
c
cc EN_R2_01_1_SIZE sizes the Stroud rule 1.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 1

      return
      end
      subroutine en_r2_02_xiu ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_02_XIU implements the Xiu precision 2 rule for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision arg
      double precision c1
      double precision delta0
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = 2.0D+00
      delta0 = 0.0D+00
      c1 = 1.0D+00

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      volume_1d = sqrt ( pi )
      volume = volume_1d ** n

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine en_r2_02_xiu_size ( n, o )

c*********************************************************************72
c
cc EN_R2_02_XIU_SIZE sizes the Xiu rule for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n + 1

      return
      end
      subroutine en_r2_03_1 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_03_1 implements the Stroud rule 3.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2 * N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi**n )

      a = volume / dble ( o )
      r = sqrt ( dble ( n ) / 2.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - r
        w(k) = a
        k = k + 1
        x(i,k) = + r
        w(k) = a
      end do

      return
      end
      subroutine en_r2_03_1_size ( n, o )

c*********************************************************************72
c
cc EN_R2_03_1_SIZE sizes the Stroud rule 3.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2 * N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n

      return
      end
      subroutine en_r2_03_2 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_03_2 implements the Stroud rule 3.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi ** n )

      a = volume / dble ( o )
      r = sqrt ( 0.5D+00 );

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2^N points.
c
      k = k + 1
      do i = 1, n
        x(i,k) = - r
      end do
      w(k) = a
      more = .true.

10    continue

      if ( more ) then
        more = .false.
        do i = n, 1, -1
          if ( x(i,k) .lt. 0.0D+00 ) then
            k = k + 1
            do i2 = 1, i - 1
              x(i2,k) = x(i2,k-1)
            end do
            x(i,k)     = + r
            do i2 = i + 1, n
              x(i2,k) = - r
            end do
            w(k) = a
            more = .true.
            go to 20
          end if
        end do
20      continue
        go to 10
      end if

      return
      end
      subroutine en_r2_03_2_size ( n, o )

c*********************************************************************72
c
cc EN_R2_03_2_SIZE sizes the Stroud rule 3.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 ** n;

      return
      end
      subroutine en_r2_03_xiu ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_03_XIU implements the Xiu precision 3 rule for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      double precision arg
      integer i
      integer j
      integer n
      integer o
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi ** n )

      do j = 1, o

        i = 0 
        do r = 1, n / 2
          arg = dble ( ( 2 * r - 1 ) * j ) * pi / dble ( n )
          i = i + 1
          x(i,j) = cos ( arg )
          i = i + 1
          x(i,j) = sin ( arg )
        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j )
          if ( n .eq. 1 ) then
            x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
          end if
        end if

      end do

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine en_r2_03_xiu_size ( n, o )

c*********************************************************************72
c
cc EN_R2_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = 2 * N.
c
c    The rule has precision P = 3.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n

      return
      end
      subroutine en_r2_05_1 ( n, option, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_1 implements the Stroud rule 5.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = N^2 + N + 2.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
c    the OPTION variable to 1 or 2.
c
c    Versions of this rule are only available for N = 2 through 7.
c
c    There is a typographical error in the reference.
c    For the second version of the rule for N = 2, the line
c      gamma =    0.313300683022281D+00
c    should read
c      gamma =    0.312200683022281D+00
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    2 <= N <= 7.
c
c    Input, integer OPTION, selects option 1 or 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision b
      double precision c
      double precision eta
      double precision gamma
      integer i
      integer i2
      integer j
      integer k
      double precision lambda
      double precision mu
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision volume
      double precision w(o)
      double precision x(n,o)
      double precision xsi

      if ( n .lt. 2 .or. 7 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
        write ( *, '(a)' ) '  2 <= N <= 7 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      if ( option .eq. 2 ) then
        if ( n .ne. 3 .and. n .ne. 5 .and. n .ne. 6 ) then 
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
          write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
          stop
        end if
      end if

      volume = sqrt ( pi ** n )

      if ( n .eq. 2 ) then
        eta =      0.446103183094540D+00
        lambda =   0.136602540378444D+01
        xsi =    - 0.366025403784439D+00
        mu =       0.198167882945871D+01
        gamma =    0.000000000000000D+00
        a =        0.328774019778636D+00 * volume
        b =        0.833333333333333D-01 * volume
        c =        0.455931355469736D-02 * volume
      else if ( n .eq. 3 .and. option .eq. 1 ) then
        eta =      0.476731294622796D+00
        lambda =   0.935429018879534D+00
        xsi =    - 0.731237647787132D+00
        mu =       0.433155309477649D+00
        gamma =    0.266922328697744D+01
        a =        0.242000000000000D+00 * volume
        b =        0.810000000000000D-01 * volume
        c =        0.500000000000000D-02 * volume
c
c  The value of gamma that follows corrects an error in the reference.
c
      else if ( n .eq. 3 .and. option .eq. 2 ) then
        eta =      0.476731294622796D+00
        lambda =   0.128679320334269D+01
        xsi =    - 0.379873463323979D+00
        mu =     - 0.192386729447751D+01
        gamma =    0.312200683022281D+00
        a =        0.242000000000000D+00 * volume
        b =        0.810000000000000D-01 * volume
        c =        0.500000000000000D-02 * volume
      else if ( n .eq. 4 ) then
        eta =      0.523945658287507D+00
        lambda =   0.119433782552719D+01
        xsi =    - 0.398112608509063D+00
        mu =     - 0.318569372920112D+00
        gamma =    0.185675837424096D+01
        a =        0.155502116982037D+00 * volume
        b =        0.777510584910183D-01 * volume
        c =        0.558227484231506D-02 * volume
      else if ( n .eq. 5 .and. option .eq. 1 ) then
        eta =      0.214972564378798D+01
        lambda =   0.464252986016289D+01
        xsi =    - 0.623201054093728D+00
        mu =     - 0.447108700673434D+00
        gamma =    0.812171426076311D+00
        a =        0.487749259189752D-03 * volume
        b =        0.487749259189752D-03 * volume
        c =        0.497073504444862D-01 * volume
      else if ( n .eq. 5 .and. option .eq. 2 ) then
        eta =      0.615369528365158D+00
        lambda =   0.132894698387445D+01
        xsi =    - 0.178394363877324D+00
        mu =     - 0.745963266507289D+00
        gamma =    0.135503972310817D+01
        a =        0.726415024414905D-01 * volume
        b =        0.726415024414905D-01 * volume
        c =        0.641509853510569D-02 * volume
      else if ( n .eq. 6 .and. option .eq. 1 ) then
        eta =      0.100000000000000D+01
        lambda =   0.141421356237309D+01
        xsi =      0.000000000000000D+00
        mu =     - 0.100000000000000D+01
        gamma =    0.100000000000000D+01
        a =        0.781250000000000D-02 * volume
        b =        0.625000000000000D-01 * volume
        c =        0.781250000000000D-02 * volume
      else if ( n .eq. 6 .and. option .eq. 2 ) then
        eta =      0.100000000000000D+01
        lambda =   0.942809041582063D+00
        xsi =    - 0.471404520791032D+00
        mu =     - 0.166666666666667D+01
        gamma =    0.333333333333333D+00
        a =        0.781250000000000D-02 * volume
        b =        0.625000000000000D-01 * volume
        c =        0.781250000000000D-02 * volume
      else if ( n .eq. 7 ) then
        eta =      0.000000000000000D+00
        lambda =   0.959724318748357D+00
        xsi =    - 0.772326488820521D+00
        mu =     - 0.141214270131942D+01
        gamma =    0.319908106249452D+00
        a =        0.111111111111111D+00 * volume
        b =        0.138888888888889D-01 * volume
        c =        0.138888888888889D-01 * volume
      end if

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2 points.
c
      k = k + 1
      do i = 1, n
        x(i,k) = - eta
      end do
      w(k) = a
      k = k + 1
      do i = 1, n
        x(i,k) = + eta
      end do
      w(k) = a
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        do i2 = 1, n
          x(i2,k) = - xsi
        end do
        x(i,k) = - lambda
        w(k) = b
        k = k + 1
        do i2 = 1, n
          x(i2,k) = + xsi
        end do
        x(i,k) = + lambda
        w(k) = b
      end do
c
c  2 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          do i2 = 1, n
            x(i2,k) = - gamma
          end do
          x(i,k) = - mu
          x(j,k) = - mu
          w(k) = c
          k = k + 1
          do i2 = 1, n
            x(i2,k) = + gamma
          end do
          x(i,k) = + mu
          x(j,k) = + mu
          w(k) = c
        end do
      end do

      return
      end
      subroutine en_r2_05_1_size ( n, option, o )

c*********************************************************************72
c
cc EN_R2_05_1_SIZE sizes the Stroud rule 5.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = N^2 + N + 2.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
c    the OPTION variable to 1 or 2.
c
c    Versions of this rule are only available for N = 2 through 7.
c
c    There is a typographical error in the reference.
c    For the second version of the rule for N = 2, the line
c      gamma =    0.313300683022281D+00
c    should read
c      gamma =    0.312200683022281D+00
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    2 <= N <= 7.
c
c    Input, integer OPTION, selects option 1 or 2.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o
      integer option

      if ( n .lt. 2 .or. 7 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  2 <= N <= 7 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      if ( option .eq. 2 ) then
        if ( n .ne. 3 .and. n .ne. 5 .and. n .ne. 6 ) then 
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
          write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
          stop
        end if
      end if

      o = n * n + n + 2

      return
      end
      subroutine en_r2_05_2 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_2 implements the Stroud rule 5.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2 * N^2 + 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision b
      double precision c
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi ** n )

      a = 2.0D+00 * volume / dble ( n + 2 )
      b = dble ( 4 - n ) * volume / 2.0D+00 
     &  / dble ( ( n + 2 ) * ( n + 2 ) )
      c = volume / dble ( ( n + 2 ) * ( n + 2 ) )

      r = sqrt ( dble ( n + 2 ) / 2.0D+00 )
      s = sqrt ( dble ( n + 2 ) / 4.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = a
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - r
        w(k) = b
        k = k + 1
        x(i,k) = + r
        w(k) = b
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - s
          x(j,k) = - s
          w(k) = c
          k = k + 1
          x(i,k) = - s
          x(j,k) = + s
          w(k) = c
          k = k + 1
          x(i,k) = + s
          x(j,k) = - s
          w(k) = c
          k = k + 1
          x(i,k) = + s
          x(j,k) = + s
          w(k) = c
        end do
      end do

      return
      end
      subroutine en_r2_05_2_size ( n, o )

c*********************************************************************72
c
cc EN_R2_05_2_SIZE sizes the Stroud rule 5.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2 * N^2 + 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 * n * n + 1

      return
      end
      subroutine en_r2_05_3 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_3 implements the Stroud rule 5.3 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N + 2 * N.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 3 <= N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision b
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_3 - Fatal error!'
        write ( *, '(a)' ) '  3 <= N is required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      a = 4.0D+00 * volume / dble ( ( n + 2 ) * ( n + 2 ) )
      b = dble ( ( n - 2 ) * ( n - 2 ) ) * volume 
     &  / dble ( 2**n ) / dble ( ( n + 2 ) * ( n + 2 ) )
      r = sqrt ( dble ( n + 2 ) / 4.0D+00 )
      s = sqrt ( dble ( n + 2 ) / 2.0D+00 / dble ( n - 2 ) )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - r
        w(k) = a
        k = k + 1
        x(i,k) = + r
        w(k) = a
      end do
c
c  2^N points.
c
      k = k + 1
      do i = 1, n
        x(i,k) = - s
      end do
      w(k) = b
      more = .true.

10    continue

      if ( more ) then

        more = .false.
        do i = n, 1, -1
          if ( x(i,k) .lt. 0.0D+00 ) then
            k = k + 1
            do i2 = 1, i - 1
              x(i2,k) = x(i2,k-1)
            end do
            x(i,k)     = + s
            do i2 = i + 1, n
              x(i2,k) = - s
            end do
            w(k) = b
            more = .true.
            go to 20
          end if
        end do
20      continue
        go to 10
      end if

      return
      end
      subroutine en_r2_05_3_size ( n, o )

c*********************************************************************72
c
cc EN_R2_05_3_SIZE sizes the Stroud rule 5.3 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N + 2 * N.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 3 <= N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_3_SIZE - Fatal error!'
        write ( *, '(a)' ) '  3 <= N is required.'
        stop
      end if

      o = 2**n + 2 * n

      return
      end
      subroutine en_r2_05_4 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_4 implements the Stroud rule 5.4 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^(N+1) - 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision b
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi ** n )

      s = sqrt ( 0.5D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
c  but do the last point separately.
c
      do i = 1, n

        r = sqrt ( dble ( i + 2 ) / 2.0D+00 )
        b = 2.0D+00 ** ( i - n ) * volume / dble ( i + 1 ) 
     &    / dble ( i + 2 )

        k = k + 1
        x(i,k) = - r
        do i2 = i + 1, n
          x(i2,k) = - s
        end do
        w(k) = b
        more = .true.

10      continue

        if ( more ) then
          more = .false.
          do j = n, i, -1
            if ( x(j,k) .lt. 0.0D+00 ) then
              k = k + 1
              do i2 = 1, j - 1
                x(i2,k) = x(i2,k-1)
              end do
              x(j,k)     =   abs ( x(j,k) )
              do i2 = j + 1, n
                x(i2,k) = - abs ( x(i2,k) )
              end do
              w(k) = b
              more = .true.
              go to 20
            end if
          end do
20        continue
          go to 10
        end if

      end do
c
c  Last point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = 2.0D+00 * volume / dble ( n + 2 )

      return
      end
      subroutine en_r2_05_4_size ( n, o )

c*********************************************************************72
c
cc EN_R2_05_4_SIZE sizes the Stroud rule 5.4 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^(N+1) - 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 2 ** ( n + 1 ) - 1

      return
      end
      subroutine en_r2_05_5 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_5 implements the Stroud rule 5.5 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = N * 2^N + 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There is a second version of this rule however it results in
c    complex abscissas, and so it has been disabled.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision b
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision n_r8
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision volume
      double precision w(o)
      double precision x(n,o)

      volume = sqrt ( pi ** n )

      n_r8 = dble ( n )

      a = 2.0D+00 * volume / ( n_r8 + 2.0D+00 )
      b =           volume / ( n_r8 + 2.0D+00 ) / ( 2.0D+00 ** n )

      option = 1

      if ( option .eq. 1 ) then
        r = sqrt ( ( n_r8 + 2.0D+00 
     &    + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &    / 2.0D+00 / n_r8 )
        s = sqrt ( ( n_r8 + 2.0D+00 
     &    -                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &    / 2.0D+00 / n_r8 )
      else if ( option .eq. 2 ) then
        r = sqrt ( ( n_r8 + 2.0D+00 
     &    - ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &    / 2.0D+00 / n_r8 )
        s = sqrt ( ( n_r8 + 2.0D+00 
     &    +                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &    / 2.0D+00 / n_r8 )
      end if

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = a
c
c  N * 2^N points:
c  N choices for location of R, 2^N choices of sign pattern.
c
      do i = 1, n

        k = k + 1
        do i2 = 1, n
          x(i2,k) = - s
        end do
        x(i,k)   = - r
        w(k) = b

        more = .true.

10      continue

        if ( more ) then
          more = .false.
          do j = n, 1, -1
            if ( x(j,k) .lt. 0.0D+00 ) then
              k = k + 1
              do i2 = 1, j - 1
                x(i2,k) = x(i2,k-1)
              end do
              x(j,k)     =   abs ( x(j,k) )
              do i2 = j + 1, n
                x(i2,k) = - abs ( x(i2,k) )
              end do
              w(k) = b
              more = .true.
              go to 20
            end if
          end do
20        continue
          go to 10
        end if

      end do

      return
      end
      subroutine en_r2_05_5_size ( n, o )

c*********************************************************************72
c
cc EN_R2_05_5_SIZE sizes the Stroud rule 5.5 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = N * 2^N + 1.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There is a second version of this rule however it results in
c    complex abscissas, and so it has been disabled.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n * 2 ** n + 1

      return
      end
      subroutine en_r2_05_6 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_05_6 implements the Stroud rule 5.6 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( N + 1 ) * 2^N.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 5 <= N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    5 <= N.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision n_r8
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision t
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_6 - Fatal error!'
        write ( *, '(a)' ) '  5 <= N is required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      n_r8 = dble ( n )

      a = volume / ( 2.0D+00 ** n ) / ( n_r8 + 1.0D+00 )

      r = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) 
     &  + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) 
     &  / 2.0D+00 / n_r8 )
      s = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) 
     &  -                      sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) 
     &  / 2.0D+00 / n_r8 )
      t = sqrt ( ( 1.0D+00 + sqrt ( 2.0D+00 ) ) / 2.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  N * 2^N points.
c
      do i = 1, n

        k = k + 1
        do i2 = 1, n
          x(i2,k) = - s
        end do
        x(i,k)   = - r
        w(k) = a

        more = .true.

10      continue

        if ( more ) then
          more = .false.
          do j = n, 1, -1
            if ( x(j,k) .lt. 0.0D+00 ) then
              k = k + 1
              do i2 = 1, j - 1
                x(i2,k) = x(i2,k-1)
              end do
              x(j,k)     =   abs ( x(j,k) )
              do i2 = j + 1, n
                x(i2,k) = - abs ( x(i2,k) )
              end do
              w(k) = a
              more = .true.
              go to 20
            end if
          end do
20        continue
          go to 10
        end if

      end do
c
c  2^N points.
c
      k = k + 1
      do i = 1, n
        x(i,k) = - t
      end do
      w(k) = a
      more = .true.

30    continue

      if ( more ) then
        more = .false.
        do j = n, 1, -1
          if ( x(j,k) .lt. 0.0D+00 ) then
            k = k + 1
            do i2 = 1, j - 1
              x(i2,k) = x(i2,k-1)
            end do
            x(j,k)     =   abs ( x(j,k) )
            do i2 = j + 1, n
              x(i2,k) = - abs ( x(i2,k) )
            end do
            w(k) = a
            more = .true.
            go to 40
          end if
        end do
40      continue
        go to 30
      end if

      return
      end
      subroutine en_r2_05_6_size ( n, o )

c*********************************************************************72
c
cc EN_R2_05_6_SIZE sizes the Stroud rule 5.6 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( N + 1 ) * 2^N.
c
c    The rule has precision P = 5.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 5 <= N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    5 <= N.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      if ( n .lt. 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_05_6_SIZE - Fatal error!'
        write ( *, '(a)' ) '  5 <= N is required.'
        stop
      end if

      o = ( 2 ** n ) * ( n + 1 )

      return
      end
      subroutine en_r2_07_1 ( n, option, o, x, w )

c*********************************************************************72
c
cc EN_R2_07_1 implements the Stroud rule 7.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N + 2 * N^2 + 1.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of the rule, chosen by setting the
c    OPTION variable to 1 or 2.  
c
c    Option 1 is only valid for N = 3, 4, 6 or 7.
c    Option 2 is only valid for N = 3 or 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    N = 3, 4, 6 or 7.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision n_r8
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision s
      double precision t
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      if ( option .eq. 1 ) then
        if ( n .ne. 3 .and. n .ne. 4 .and. n .ne. 6 .and. n .ne. 7 ) 
     &    then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
          write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
          stop
        end if
      end if

      if ( option .eq. 2 ) then
        if ( n .ne. 3 .and. n .ne. 4 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
          write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
          stop
        end if
      end if

      volume = sqrt ( pi ** n )

      n_r8 = dble ( n )

      if ( option .eq. 1 ) then
        r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) - ( n_r8 - 2.0D+00 ) 
     &    * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 
     &    / ( 5.0D+00 - n_r8 ) )
        s = sqrt ( ( 3.0D+00 *             n_r8   -          2.0D+00   
     &    * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 
     &    / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
        t = sqrt ( ( 6.0D+00 + sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) 
     &    / 2.0D+00 )
      else if ( option .eq. 2 ) then
        r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) + ( n_r8 - 2.0D+00 ) 
     &    * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 
     &    / ( 5.0D+00 - n_r8 ) )
        s = sqrt ( ( 3.0D+00 *             n_r8   +          2.0D+00   
     &    * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 
     &    / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
        t = sqrt ( ( 6.0D+00 - sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) 
     &    / 2.0D+00 )
      end if

      b = ( 8.0D+00 - n_r8 ) * volume / 8.0D+00 / r ** 6
      c = volume / 2.0D+00 ** ( n + 3 ) / s ** 6
      d = volume / 16.0D+00 / t ** 6
      a = volume - 2.0D+00 * n_r8 * b - 2.0D+00 ** n * c - 2.0D+00 
     &  * n_r8 * ( n_r8 - 1.0D+00 ) * d

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = a
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - r
        w(k) = b
        k = k + 1
        x(i,k) = + r
        w(k) = b
      end do
c
c  2^N points.
c
      k = k + 1
      x(1:n,k) = - s
      w(k) = c
      more = .true.

10    continue

      if ( more ) then
        more = .false.
        do i = n, 1, -1
          if ( x(i,k) .lt. 0.0D+00 ) then
            k = k + 1
            do i2 = 1, i - 1
              x(i2,k) = x(i2,k-1)
            end do
            x(i,k)     =   abs ( x(i,k) )
            do i2 = i + 1, n
              x(i2,k) = - abs ( x(i2,k) )
            end do
            w(k) = c
            more = .true.
            go to 20
          end if
        end do
20      continue
        go to 10
      end if
c
c  2 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - t
          x(j,k) = - t
          w(k) = d
          k = k + 1
          x(i,k) = - t
          x(j,k) = + t
          w(k) = d
          k = k + 1
          x(i,k) = + t
          x(j,k) = - t
          w(k) = d
          k = k + 1
          x(i,k) = + t
          x(j,k) = + t
          w(k) = d
        end do
      end do

      return
      end
      subroutine en_r2_07_1_size ( n, option, o )

c*********************************************************************72
c
cc EN_R2_07_1_SIZE sizes the Stroud rule 7.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^N + 2 * N^2 + 1.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of the rule, chosen by setting the
c    OPTION variable to 1 or 2.  
c
c    Option 1 is only valid for N = 3, 4, 6 or 7.
c    Option 2 is only valid for N = 3 or 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    N = 3, 4, 6 or 7.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o
      integer option

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      if ( option .eq. 1 ) then
        if ( n .ne. 3 .and. n .ne. 4 .and. n .ne. 6 .and. n .ne. 7 ) 
     &    then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
          write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
          stop
        end if
      end if

      if ( option .eq. 2 ) then
        if ( n .ne. 3 .and. n .ne. 4 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
          write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
          stop
        end if
      end if

      o = 2 ** n + 2 * n ** 2 + 1

      return
      end
      subroutine en_r2_07_2 ( n, o, x, w )

c*********************************************************************72
c
cc EN_R2_07_2 implements the Stroud rule 7.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^(N+1) + 4 * N^2.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 3 <= N.
c
c    The reference has a typographical error in the description of this rule.
c    The formula:
c
c      (t,t,t,...,t)FS
c
c    should read
c
c      (t,t,0,...,0)FS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision a1
      double precision a2
      double precision b
      double precision c
      double precision d
      integer i
      integer i2
      integer j
      integer k
      logical more
      double precision n_r8
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision rho1
      double precision rho2
      double precision s
      double precision t
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_2 - Fatal error!'
        write ( *, '(a)' ) '  3 <= N is required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      n_r8 = dble ( n )

      rho1 = sqrt ( ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 
     &  * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 )
      rho2 = sqrt ( ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 
     &  * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 )
      a1 = ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &  / 2.0D+00 / ( n_r8 + 2.0D+00 )
      a2 = ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) 
     &  / 2.0D+00 / ( n_r8 + 2.0D+00 )

      r = 1.0D+00
      s = sqrt ( 1.0D+00 / n_r8 )
      t = sqrt ( 0.5D+00 )
      b = ( 8.0D+00 - n_r8 ) * volume / n_r8 / ( n_r8 + 2.0D+00 ) 
     &  / ( n_r8 + 4.0D+00 )
      c = n_r8 ** 3 * volume / 2.0D+00 ** n / n_r8 / ( n_r8 + 2.0D+00 ) 
     &  / ( n_r8 + 4.0D+00 )
      d = 4.0D+00 * volume / n_r8 / ( n_r8 + 2.0D+00 ) 
     &  / ( n_r8 + 4.0D+00 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  2 * 2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - rho1 * r
        w(k) = a1 * b
        k = k + 1
        x(i,k) = - rho2 * r
        w(k) = a2 * b
        k = k + 1
        x(i,k) = + rho1 * r
        w(k) = a1 * b
        k = k + 1
        x(i,k) = + rho2 * r
        w(k) = a2 * b
      end do
c
c  2 * 2^N points.
c
      k = k + 1
      x(1:n,k) = - rho1 * s
      w(k) = a1 * c
      k = k + 1
      x(1:n,k) = - rho2 * s
      w(k) = a2 * c
      more = .true.

10    continue

      if ( more ) then
        more = .false.
        do i = n, 1, -1
          if ( x(i,k) .lt. 0.0D+00 ) then
            k = k + 1
            do i2 = 1, i - 1
              x(i2,k) =     x(i2,k-2)
            end do
            x(i,k)     =   abs ( x(i,k) )
            do i2 = i + 1, n
              x(i2,k) = - abs ( x(i2,k) )
            end do
            w(k) = a1 * c
            k = k + 1
            do i2 = 1, i - 1
              x(i2,k) =     x(i2,k-2)
            end do
            x(i,k)     =   abs ( x(i,k) )
            do i2 = i + 1, n
              x(i2,k) = - abs ( x(i2,k) )
            end do
            w(k) = a2 * c
            more = .true.
            go to 20
          end if
        end do
20      continue
        go to 10
      end if
c
c  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - rho1 * t
          x(j,k) = - rho1 * t
          w(k) = a1 * d
          k = k + 1
          x(i,k) = - rho1 * t
          x(j,k) = + rho1 * t
          w(k) = a1 * d
          k = k + 1
          x(i,k) = + rho1 * t
          x(j,k) = - rho1 * t
          w(k) = a1 * d
          k = k + 1
          x(i,k) = + rho1 * t
          x(j,k) = + rho1 * t
          w(k) = a1 * d
          k = k + 1
          x(i,k) = - rho2 * t
          x(j,k) = - rho2 * t
          w(k) = a2 * d
          k = k + 1
          x(i,k) = - rho2 * t
          x(j,k) = + rho2 * t
          w(k) = a2 * d
          k = k + 1
          x(i,k) = + rho2 * t
          x(j,k) = - rho2 * t
          w(k) = a2 * d
          k = k + 1
          x(i,k) = + rho2 * t
          x(j,k) = + rho2 * t
          w(k) = a2 * d
        end do
      end do

      return
      end
      subroutine en_r2_07_2_size ( n, o )

c*********************************************************************72
c
cc EN_R2_07_2_SIZE sizes the Stroud rule 7.2 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = 2^(N+1) + 4 * N^2.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The rule requires 3 <= N.
c
c    The reference has a typographical error in the description of this rule.
c    The formula:
c
c      (t,t,t,...,t)FS
c
c    should read
c
c      (t,t,0,...,0)FS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_2_SIZE - Fatal error!'
        write ( *, '(a)' ) '  3 <= N is required.'
        stop
      end if

      o = 2 ** ( n + 1 ) + 4 * n * n

      return
      end
      subroutine en_r2_07_3 ( n, option, o, x, w )

c*********************************************************************72
c
cc EN_R2_07_3 implements the Stroud rule 7.3 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 6.
c
c     N    O
c    __  ___
c     3   45
c     4   97
c     5  181
c     6  305
c
c    The reference has a typographical error for N = 5, OPTION 1, B4:
c
c      -(1)0.736330882774831
c
c    should read
c
c      (-1)0.736330882774831
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 6.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision b0
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      integer i
      integer j
      integer k
      integer l
      logical more
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u
      double precision v
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 3 .or. 6 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 6 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      if ( n .eq. 3 .and. option .eq. 1 ) then
        u =    0.524647623275290D+00
        v =    0.165068012388578D+01
        b0 = - 0.166705761599566D+02
        b1 =   0.100296981655678D+02
        b2 =   0.161699246687754D+00
        b3 = - 0.604719151221535D+01
        b4 =   0.234381399489666D-01
        b5 =   0.417194501880647D+01
      else if ( n .eq. 3 .and. option .eq. 2 ) then
        u =    0.165068012388578D+01
        v =    0.524647623275290D+00
        b0 =   0.166705761599566D+02
        b1 =   0.178903161957074D+00
        b2 = - 0.665808190965810D+01
        b3 =   0.148361823143070D-01
        b4 =   0.229669852539758D+01
        b5 =   0.430097881732984D-02
      else if ( n .eq. 4 .and. option .eq. 1 ) then
        u  =   0.524647623275290D+00
        v  =   0.165068012388578D+01
        b0 = - 0.167539329651562D+03
        b1 =   0.687922329603575D+02
        b2 =   0.203518409659014D+00
        b3 = - 0.255075279116885D+02
        b4 =   0.415430214106084D-01
        b5 =   0.739458001434961D+01
      else if ( n .eq. 4 .and. option .eq. 2 ) then
        u =    0.165068012388578D+01
        v =    0.524647623275290D+00
        b0 =   0.688432856406677D+02
        b1 =   0.294997847268286D+00
        b2 = - 0.199427272118378D+02
        b3 =   0.110498755408511D-01
        b4 =   0.407079214570997D+01
        b5 =   0.762328646743931D-02
      else if ( n .eq. 5 .and. option .eq. 1 ) then
        u  =   0.524647623275290D+00
        v  =   0.165068012388578D+01
        b0 = - 0.826940846964452D+03
        b1 =   0.264779097660331D+03
        b2 =   0.213460812375320D+00
        b3 = - 0.714240197186780D+02
        b4 =   0.736330882774831D-01
        b5 =   0.131065518222629D+02
      else if ( n .eq. 5 .and. option .eq. 2 ) then
        u =    0.165068012388578D+01
        v =    0.524647623275290D+00
        b0 =   0.220502344940121D+03
        b1 =   0.537746975313769D+00
        b2 = - 0.497781460739792D+02
        b3 = - 0.743845245712926D-02
        b4 =   0.721529121489956D+01
        b5 =   0.135119234557687D-01
      else if ( n .eq. 6 .and. option .eq. 1 ) then
        u  =   0.524647623275290D+00
        v  =   0.165068012388578D+01
        b0 = - 0.309679578630802E+04
        b1 =   0.815423321880237D+03
        b2 =   0.117326937169073D+00
        b3 = - 0.173057295296448D+03
        b4 =   0.130511250871491D+00
        b5 =   0.232307582494626D+02
      else if ( n .eq. 6 .and. option .eq. 2 ) then
        u =    0.165068012388578D+01
        v =    0.524647623275290D+00
        b0 =   0.616293651884027D+03
        b1 =   0.107529736766179D+01
        b2 = - 0.113807008098269D+03
        b3 = - 0.610828352270520D-01
        b4 =   0.127887706992535D+02
        b5 =   0.239492607623178D-01
      end if

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = b0
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - u
        w(k) = b1
        k = k + 1
        x(i,k) = + u
        w(k) = b1
      end do
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - v
        w(k) = b2
        k = k + 1
        x(i,k) = + v
        w(k) = b2
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          w(k) = b3
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          w(k) = b3
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          w(k) = b3
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          w(k) = b3
        end do
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          w(k) = b4
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          w(k) = b4
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          w(k) = b4
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          w(k) = b4
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b5
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b5
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b5
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b5
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b5
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b5
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b5
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b5
          end do
        end do
      end do

      return
      end
      subroutine en_r2_07_3_size ( n, option, o )

c*********************************************************************72
c
cc EN_R2_07_3_SIZE sizes the Stroud rule 7.3 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
c
c    The rule has precision P = 7.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 6.
c
c     N    O
c    __  ___
c     3   45
c     4   97
c     5  181
c     6  305
c
c    The reference has a typographical error for N = 5, OPTION 1, B4:
c
c      -(1)0.736330882774831
c
c    should read
c
c      (-1)0.736330882774831
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 6.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o
      integer option

      if ( n .lt. 3 .or. 6 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 6 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      o = ( 4 * n ** 3 + 8 * n + 3 ) / 3

      return
      end
      subroutine en_r2_09_1 ( n, option, o, x, w )

c*********************************************************************72
c
cc EN_R2_09_1 implements the Stroud rule 9.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
c
c    The rule has precision P = 9.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the 
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 6.
c
c     N    O
c    __  ___
c     3   77
c     4  193
c     5  421
c     6  825
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 6.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision b0
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision b8
      integer i
      integer j
      integer k
      integer l
      integer m
      logical more
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u
      double precision v
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 3 .or. 6 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 6 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      if ( n .eq. 3 ) then
        u =    0.202018287045609D+01
        v =    0.958572464613819D+00
        b0 =   0.676448734429924D+00
        b1 =   0.511989106291551D-02
        b2 =   0.448595723493744D+00
        b3 =   0.235223454595606D-03
        b4 =   0.915390713080005D-01
        b5 =   0.139208199920793D-01
        b6 =   0.235223454595606D-03
        b7 =   0.915390713080008D-01
        b8 =   0.000000000000000D+00
      else if ( n .eq. 4 .and. option .eq. 1 ) then
        u =    0.202018287045609D+01
        v =    0.958572464613819D+00
        b0 = - 0.860452945007048D+00
        b1 = - 0.405511998533795D-01
        b2 =   0.107026475449715D+01
        b3 =   0.138974239307092D-03
        b4 = - 0.162248779448181D+00
        b5 =   0.246740110027234D-01
        b6 =   0.138974239307094D-03
        b7 =   0.162248779448181D+00
        b8 =   0.138974239307094D-03
      else if ( n .eq. 4 .and. option .eq. 2 ) then
        u =    0.958572464613819D+00
        v =    0.202018287045609D+01
        b0 =   0.265029088766810D-02
        b1 =   0.637601342635332D+00
        b2 = - 0.394394059389228D-01
        b3 =   0.540829264827264D-01
        b4 = - 0.416922717921281D-03
        b5 =   0.246740110027234D-01
        b6 =   0.540829264827270D-01
        b7 =   0.416922717921281D-03
        b8 =   0.540829264827269D-01
      else if ( n .eq. 5 .and. option .eq. 1 ) then
        u =    0.202018287045609D+01
        v =    0.958572464613819D+00
        b0 = - 0.827347006200826D+01
        b1 = - 0.160820174530905D+00
        b2 =   0.353499863758467D+01
        b3 =   0.738976276909564D-03
        b4 = - 0.862735421812943D+00
        b5 =   0.437335458190621D-01
        b6 = - 0.246325425636523D-03
        b7 =   0.287578473937648D+00
        b8 =   0.246325425636523D-03
      else if ( n .eq. 5 .and. option .eq. 2 ) then
        u =    0.958572464613819D+00
        v =    0.202018287045609D+01
        b0 = - 0.624416791055272D+00
        b1 =   0.467494915583104D+00
        b2 = - 0.152937760910536D+00
        b3 =   0.287578473937646D+00
        b4 = - 0.221692883072871D-02
        b5 =   0.437335458190621D-01
        b6 = - 0.958594913125490D-01
        b7 =   0.738976276909568D-03
        b8 =   0.958594913125492D-01
      else if ( n .eq. 6 .and. option .eq. 1 ) then
        u =    0.202018287045609D+01
        v =    0.958572464613819D+00
        b0 = - 0.361840434143098D+02
        b1 = - 0.447936529138517D+00
        b2 =   0.112077863004144D+02
        b3 =   0.392940404320855D-02
        b4 = - 0.254859786784158D+01
        b5 =   0.775156917007496D-01
        b6 = - 0.130980134773619D-02
        b7 =   0.509719573568315D+00
        b8 =   0.436600449245395D-03
      else if ( n .eq. 6 .and. option .eq. 2 ) then
        u =    0.958572464613819D+00
        v =    0.202018287045609D+01
        b0 =   0.448873836333650D+01
        b1 = - 0.238473566140736D+01
        b2 = - 0.413008493198885D+00
        b3 =   0.152915872070494D+01
        b4 = - 0.654900673868093D-02
        b5 =   0.775156917007496D-01
        b6 = - 0.509719573568314D+00
        b7 =   0.130980134773618D-02
        b8 =   0.169906524522772D+00
      end if

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = b0
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - u
        w(k) = b1
        k = k + 1
        x(i,k) = + u
        w(k) = b1
      end do
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - v
        w(k) = b2
        k = k + 1
        x(i,k) = + v
        w(k) = b2
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          w(k) = b3
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          w(k) = b3
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          w(k) = b3
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          w(k) = b3
        end do
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          w(k) = b4
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          w(k) = b4
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          w(k) = b4
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          w(k) = b4
        end do
      end do
c
c  4 * ( N * ( N - 1 ) ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - v
          w(k) = b5
          k = k + 1
          x(i,k) = - u
          x(j,k) = + v
          w(k) = b5
          k = k + 1
          x(i,k) = + u
          x(j,k) = - v
          w(k) = b5
          k = k + 1
          x(i,k) = + u
          x(j,k) = + v
          w(k) = b5
          k = k + 1
          x(i,k) = - v
          x(j,k) = - u
          w(k) = b5
          k = k + 1
          x(i,k) = - v
          x(j,k) = + u
          w(k) = b5
          k = k + 1
          x(i,k) = + v
          x(j,k) = - u
          w(k) = b5
          k = k + 1
          x(i,k) = + v
          x(j,k) = + u
          w(k) = b5
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b6
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b6
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b6
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b6
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b6
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b6
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b6
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b6
          end do
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - v
            x(j,k) = - v
            x(l,k) = - v
            w(k) = b7
            k = k + 1
            x(i,k) = - v
            x(j,k) = - v
            x(l,k) = + v
            w(k) = b7
            k = k + 1
            x(i,k) = - v
            x(j,k) = + v
            x(l,k) = - v
            w(k) = b7
            k = k + 1
            x(i,k) = - v
            x(j,k) = + v
            x(l,k) = + v
            w(k) = b7
            k = k + 1
            x(i,k) = + v
            x(j,k) = - v
            x(l,k) = - v
            w(k) = b7
            k = k + 1
            x(i,k) = + v
            x(j,k) = - v
            x(l,k) = + v
            w(k) = b7
            k = k + 1
            x(i,k) = + v
            x(j,k) = + v
            x(l,k) = - v
            w(k) = b7
            k = k + 1
            x(i,k) = + v
            x(j,k) = + v
            x(l,k) = + v
            w(k) = b7
          end do
        end do
      end do
c
c  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
c
      do i = 1, n - 3
        do j = i + 1, n - 2
          do l = j + 1, n - 1
            do m = l + 1, n
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b8
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b8
            end do
          end do
        end do
      end do

      return
      end
      subroutine en_r2_09_1_size ( n, option, o )

c*********************************************************************72
c
cc EN_R2_09_1_SIZE sizes the Stroud rule 9.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
c
c    The rule has precision P = 9.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the 
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 6.
c
c     N    O
c    __  ___
c     3   77
c     4  193
c     5  421
c     6  825
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 6.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o
      integer option

      if ( n .lt. 3 .or. 6 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 6 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      o = ( 2 * n ** 4 - 4 * n ** 3 + 22 * n ** 2 - 8 * n + 3 ) / 3

      return
      end
      subroutine en_r2_11_1 ( n, option, o, x, w )

c*********************************************************************72
c
cc EN_R2_11_1 implements the Stroud rule 11.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
c
c    The rule has precision P = 11.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 5.
c
c     N    O
c    __  ___
c     3  151
c     4  417
c     5  983
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 5.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision b0
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision b8
      double precision b9
      double precision b10
      double precision b11
      double precision b12
      double precision b13
      double precision b14
      double precision b15
      integer i
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer j
      integer k
      integer l
      integer m
      logical more
      integer option
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u
      double precision v
      double precision volume
      double precision w2
      double precision w(o)
      double precision x(n,o)

      if ( n .lt. 3 .or. 5 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 5 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      volume = sqrt ( pi ** n )

      if ( n .eq. 3 .and. option .eq. 1 ) then
        u =     0.235060497367449D+01
        v =     0.436077411927617D+00
        w2 =    0.133584907401370D+01
        b0 =  - 0.881591029957858D+01
        b1 =  - 0.751996143360650D-01
        b2 =    0.621743189471515D+01
        b3 =    0.241426451456494D+00
        b4 =  - 0.120709739276065D-02
        b5 =  - 0.427751221210138D+01
        b6 =    0.550169924840163D-01
        b7 =    0.237084999634707D-01
        b8 =  - 0.169791992887741D-02
        b9 =  - 0.252266276123350D-04
        b10 =   0.326777873717691D+01
        b11 =   0.968469949206802D-02
        b12 =   0.789754514877422D-03
        b13 =   0.000000000000000D+00
        b14 =   0.000000000000000D+00
        b15 =   0.000000000000000D+00
      else if ( n .eq. 3 .and. option .eq. 2 ) then
        u =     0.235060497367449D+01
        v =     0.133584907401370D+01
        w2 =    0.436077411927617D+00
        b0 =  - 0.141214037032900D+02
        b1 =  - 0.803730274707282D-01
        b2 =    0.235546545595906D+00
        b3 =    0.888123191556611D+01
        b4 =    0.142467131155533D-03
        b5 =    0.582993124006494D-01
        b6 =  - 0.561099173155661D+01
        b7 =  - 0.204028691521686D-02
        b8 =    0.252880089932256D-01
        b9 =  - 0.814378678627283D-04
        b10 =   0.804353953375146D-02
        b11 =   0.393451849690453D+01
        b12 =   0.171183493169724D-03
        b13 =   0.000000000000000D+00
        b14 =   0.000000000000000D+00
        b15 =   0.000000000000000D+00
      else if ( n .eq. 4 .and. option .eq. 1 ) then
        u =     0.235060497367449D+01
        v =     0.436077411927617D+00
        w2 =    0.133584907401370D+01
        b0 =    0.241502736147339D+03
        b1 =  - 0.196095938531478D+00
        b2 =  - 0.128675737999280D+03
        b3 =    0.307568784278696D+00
        b4 =  - 0.480908422319460D-02
        b5 =    0.698087019367085D+02
        b6 =    0.631837143743771D-01
        b7 =    0.392226151971179D-01
        b8 =  - 0.300948471646799D-02
        b9 =  - 0.650235306755170D-04
        b10 = - 0.386951974646715D+02
        b11 =   0.171656829095787D-01
        b12 =   0.139980343116450D-02
        b13 =   0.101552487093372D-04
        b14 =   0.222435922356439D+02
        b15 =   0.000000000000000D+00
      else if ( n .eq. 4 .and. option .eq. 2 ) then
        u =     0.235060497367449D+01
        v =     0.133584907401370D+01
        w2 =    0.436077411927617D+00
        b0 =  - 0.151944464736584D+03
        b1 =  - 0.223498438689039D+00
        b2 =    0.243574919068010D+00
        b3 =    0.634373877008693D+02
        b4 =  - 0.782065187814018D-04
        b5 =    0.911833754536616D-01
        b6 =  - 0.238927288245914D+02
        b7 =  - 0.422314408318853D-02
        b8 =    0.448218289217760D-01
        b9 =  - 0.138053374667391D-03
        b10 =   0.607473265800655D-02
        b11 =   0.697375246129742D+01
        b12 =   0.303414841680135D-03
        b13 = - 0.314574391771792D-05
        b14 =   0.409103498175100D-02
        b15 =   0.000000000000000D+00
      else if ( n .eq. 5 .and. option .eq. 1 ) then
        u =     0.235060497367449D+01
        v =     0.436077411927617D+00
        w2 =    0.133584907401370D+01
        b0 =    0.255885269311763E+04
        b1 =  - 0.439598677491526D+00
        b2 =  - 0.106541406144610E+04
        b3 =    0.453540909054264D+00
        b4 =  - 0.132100905623778D-01
        b5 =    0.418606568954203D+03
        b6 =    0.511394563043680D-01
        b7 =    0.645581013845604D-01
        b8 =  - 0.533417277494500D-02
        b9 =  - 0.137981626254496D-03
        b10 = - 0.147436933189884D+03
        b11 =   0.304253807765057D-01
        b12 =   0.248108698207828D-02
        b13 =   0.113652094546015D-04
        b14 =   0.394257407160391D+02
        b15 =   0.331725011358320D-05
      else if ( n .eq. 5 .and. option .eq. 2 ) then
        u =     0.235060497367449D+01
        v =     0.133584907401370D+01
        w2 =    0.436077411927617D+00
        b0 =  - 0.761305347548192D+03
        b1 =  - 0.536360805019297D+00
        b2 =    0.110669832078736D+00
        b3 =    0.246421088923968D+03
        b4 =  - 0.773649327968607D-03
        b5 =    0.169088641205970D+00
        b6 =  - 0.670700680243651D+02
        b7 =  - 0.856090560229205D-02
        b8 =    0.794446232770302D-01
        b9 =  - 0.220272863263544D-03
        b10 = - 0.373515812228225D-02
        b11 =   0.123606544052884D+02
        b12 =   0.537788804557843D-03
        b13 = - 0.122101861480881D-04
        b14 =   0.725117070759373D-02
        b15 =   0.331725011358320D-05
      end if

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
c     do i = 1, n
c       x(i,k) = 0.0D+00
c     end do
      w(k) = b0
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - u
        w(k) = b1
        k = k + 1
        x(i,k) = + u
        w(k) = b1
      end do
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - v
        w(k) = b2
        k = k + 1
        x(i,k) = + v
        w(k) = b2
      end do
c
c  2 * N points.
c
      do i = 1, n
        k = k + 1
        x(i,k) = - w2
        w(k) = b3
        k = k + 1
        x(i,k) = + w2
        w(k) = b3
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          w(k) = b4
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          w(k) = b4
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          w(k) = b4
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          w(k) = b4
        end do
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          w(k) = b5
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          w(k) = b5
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          w(k) = b5
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          w(k) = b5
        end do
      end do
c
c  4 * ( N * ( N - 1 ) / 2 ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - w2
          x(j,k) = - w2
          w(k) = b6
          k = k + 1
          x(i,k) = - w2
          x(j,k) = + w2
          w(k) = b6
          k = k + 1
          x(i,k) = + w2
          x(j,k) = - w2
          w(k) = b6
          k = k + 1
          x(i,k) = + w2
          x(j,k) = + w2
          w(k) = b6
        end do
      end do
c
c  4 * ( N * ( N - 1 ) ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - v
          w(k) = b7
          k = k + 1
          x(i,k) = - u
          x(j,k) = + v
          w(k) = b7
          k = k + 1
          x(i,k) = + u
          x(j,k) = - v
          w(k) = b7
          k = k + 1
          x(i,k) = + u
          x(j,k) = + v
          w(k) = b7
          k = k + 1
          x(i,k) = - v
          x(j,k) = - u
          w(k) = b7
          k = k + 1
          x(i,k) = - v
          x(j,k) = + u
          w(k) = b7
          k = k + 1
          x(i,k) = + v
          x(j,k) = - u
          w(k) = b7
          k = k + 1
          x(i,k) = + v
          x(j,k) = + u
          w(k) = b7
        end do
      end do
c
c  4 * ( N * ( N - 1 ) ) points.
c
      do i = 1, n - 1
        do j = i + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - w2
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + w2
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - w2
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + w2
          w(k) = b8
          k = k + 1
          x(i,k) = - w2
          x(j,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - w2
          x(j,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + w2
          x(j,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + w2
          x(j,k) = + u
          w(k) = b8
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b9
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b9
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b9
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b9
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b9
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b9
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b9
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b9
          end do
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - v
            x(j,k) = - v
            x(l,k) = - v
            w(k) = b10
            k = k + 1
            x(i,k) = - v
            x(j,k) = - v
            x(l,k) = + v
            w(k) = b10
            k = k + 1
            x(i,k) = - v
            x(j,k) = + v
            x(l,k) = - v
            w(k) = b10
            k = k + 1
            x(i,k) = - v
            x(j,k) = + v
            x(l,k) = + v
            w(k) = b10
            k = k + 1
            x(i,k) = + v
            x(j,k) = - v
            x(l,k) = - v
            w(k) = b10
            k = k + 1
            x(i,k) = + v
            x(j,k) = - v
            x(l,k) = + v
            w(k) = b10
            k = k + 1
            x(i,k) = + v
            x(j,k) = + v
            x(l,k) = - v
            w(k) = b10
            k = k + 1
            x(i,k) = + v
            x(j,k) = + v
            x(l,k) = + v
            w(k) = b10
          end do
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - w2
            x(j,k) = - w2
            x(l,k) = - w2
            w(k) = b11
            k = k + 1
            x(i,k) = - w2
            x(j,k) = - w2
            x(l,k) = + w2
            w(k) = b11
            k = k + 1
            x(i,k) = - w2
            x(j,k) = + w2
            x(l,k) = - w2
            w(k) = b11
            k = k + 1
            x(i,k) = - w2
            x(j,k) = + w2
            x(l,k) = + w2
            w(k) = b11
            k = k + 1
            x(i,k) = + w2
            x(j,k) = - w2
            x(l,k) = - w2
            w(k) = b11
            k = k + 1
            x(i,k) = + w2
            x(j,k) = - w2
            x(l,k) = + w2
            w(k) = b11
            k = k + 1
            x(i,k) = + w2
            x(j,k) = + w2
            x(l,k) = - w2
            w(k) = b11
            k = k + 1
            x(i,k) = + w2
            x(j,k) = + w2
            x(l,k) = + w2
            w(k) = b11
          end do
        end do
      end do
c
c  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
c
      do i = 1, n - 2
        do j = i + 1, n - 1
          do l = j + 1, n
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = - v
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = - u
            x(l,k) = + v
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = - v
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = + u
            x(l,k) = + v
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = - v
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = - u
            x(l,k) = + v
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = - v
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = + u
            x(l,k) = + v
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = - v
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = - v
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = + v
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = - u
            x(j,k) = + v
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = - v
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = - v
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = + v
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = + u
            x(j,k) = + v
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = - v
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = - v
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = - v
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = - v
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = + v
            x(j,k) = - u
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = + v
            x(j,k) = - u
            x(l,k) = + u
            w(k) = b12
            k = k + 1
            x(i,k) = + v
            x(j,k) = + u
            x(l,k) = - u
            w(k) = b12
            k = k + 1
            x(i,k) = + v
            x(j,k) = + u
            x(l,k) = + u
            w(k) = b12
          end do
        end do
      end do
c
c  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
c
      do i = 1, n - 3
        do j = i + 1, n - 2
          do l = j + 1, n - 1
            do m = l + 1, n
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = - u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = - u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = - u
              x(m,k) = + u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = - u
              w(k) = b13
              k = k + 1
              x(i,k) = + u
              x(j,k) = + u
              x(l,k) = + u
              x(m,k) = + u
              w(k) = b13
            end do
          end do
        end do
      end do
c
c  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
c
      do i = 1, n - 3
        do j = i + 1, n - 2
          do l = j + 1, n - 1
            do m = l + 1, n
              k = k + 1
              x(i,k) = - v
              x(j,k) = - v
              x(l,k) = - v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = - v
              x(l,k) = - v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = - v
              x(l,k) = + v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = - v
              x(l,k) = + v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = + v
              x(l,k) = - v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = + v
              x(l,k) = - v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = + v
              x(l,k) = + v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = - v
              x(j,k) = + v
              x(l,k) = + v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = - v
              x(l,k) = - v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = - v
              x(l,k) = - v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = - v
              x(l,k) = + v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = - v
              x(l,k) = + v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = + v
              x(l,k) = - v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = + v
              x(l,k) = - v
              x(m,k) = + v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = + v
              x(l,k) = + v
              x(m,k) = - v
              w(k) = b14
              k = k + 1
              x(i,k) = + v
              x(j,k) = + v
              x(l,k) = + v
              x(m,k) = + v
              w(k) = b14
            end do
          end do
        end do
      end do
c
c  All quintuples UUUUU with 32 sign combinations.
c
      do i1 = 1, n - 4
        do i2 = i1 + 1, n - 3
          do i3 = i2 + 1, n - 2
            do i4 = i3 + 1, n - 1
              do i5 = i4 + 1, n
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = - u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = - u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = - u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = - u
                x(i5,k) = + u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = - u
                w(k) = b15
                k = k + 1
                x(i1,k) = + u
                x(i2,k) = + u
                x(i3,k) = + u
                x(i4,k) = + u
                x(i5,k) = + u
                w(k) = b15
              end do
            end do
          end do
        end do
      end do

      return
      end
      subroutine en_r2_11_1_size ( n, option, o )

c*********************************************************************72
c
cc EN_R2_11_1_SIZE sizes the Stroud rule 11.1 for region EN_R2.
c
c  Discussion:
c
c    The rule has order 
c
c      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
c
c    The rule has precision P = 11.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    There are two versions of each rule, chosen by setting the
c    OPTION variable to 1 or 2.
c
c    The rule as tabulated by Stenger is available for N = 2 through 20.
c    This function accepts N = 3 through 5.
c
c     N    O
c    __  ___
c     3  151
c     4  417
c     5  983
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c    3 <= N <= 5.
c
c    Input, integer OPTION, chooses rule option 1 or 2.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o
      integer option

      if ( n .lt. 3 .or. 5 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  3 <= N <= 5 required.'
        stop
      end if

      if ( option .lt. 1 .or. 2 .lt. option ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
        write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
        stop
      end if

      o = ( 4 * n ** 5 - 20 * n ** 4 + 140 * n ** 3 - 130 * n ** 2 
     &  + 96 * n + 15 ) / 15


      return
      end
      subroutine en_r2_monomial_integral ( n, alpha, value )

c*********************************************************************72
c
cc EN_R2_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_R2.
c
c  Discussion:
c
c    ALPHA is the set of polynomial exponents.
c
c    EN_R2 is the entire N-dimensional space with weight function
c
c      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
c
c    The integral to be evaluated is
c
c      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
c        * x(n)^alpha(n) * w(x) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer ALPHA(N), the polynomial exponents.
c    0 <= ALPHA(*).
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      integer alpha(n)
      double precision arg
      integer i
      double precision r8_gamma
      double precision value

      do i = 1, n
        if ( alpha(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EN_R2_MONOMIAL_INTEGRAL - Fatal error!'
          write ( *, '(a)' ) '  Some ALPHA(I) .lt. 0.'
          stop
        end if
      end do

      do i = 1, n
        if ( mod ( alpha(i), 2 ) .eq. 1 ) then
          value = 0.0D+00
          return
        end if
      end do

      value = 1.0D+00
      do i = 1, n
        arg = ( dble ( alpha(i) + 1 ) ) / 2.0D+00
        value = value * r8_gamma ( arg )
      end do

      return
      end
      subroutine ep1_glg_monomial_integral ( expon, alpha, exact )

c*********************************************************************72
c
cc EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
c
c  Discussion:
c
c    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
c
c      w(alpha;x) = x^alpha exp ( - x )
c
c    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer EXPON, the exponent.
c    0 <= EXPON.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision arg
      double precision exact
      integer expon
      double precision r8_gamma

      arg = alpha + dble ( expon + 1 )

      exact = r8_gamma ( arg )

      return
      end
      subroutine ep1_lag_monomial_integral ( expon, value )

c*********************************************************************72
c
cc EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
c
c  Discussion:
c
c    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
c
c      w(x) = exp ( - x )
c
c    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer EXPON, the exponent.
c    0 <= EXPON.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer expon
      double precision r8_factorial
      double precision value

      value = r8_factorial ( expon )

      return
      end
      subroutine epn_glg_00_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, volume )
      volume = volume ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = 1.0D+00
      end do
      w(k) = volume

      return
      end
      subroutine epn_glg_00_1_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine epn_glg_01_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision value1
      double precision value2
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, value1 )
      volume = value1 ** n

      expon = 1
      call ep1_glg_monomial_integral ( expon, alpha, value2 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = value2 / value1
      end do
      w(k) = volume

      return
      end
      subroutine epn_glg_01_1_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine epn_glg_02_xiu ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision arg
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_02_XIU - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) * pi ) / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = - 1.0D+00
      delta0 = alpha + 1.0D+00
      c1 = - alpha - 1.0D+00

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
      volume = volume_1d ** n

      do j = 1, o
        w(j) = volume / dble ( o )
      end do

      return
      end
      subroutine epn_glg_02_xiu_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_02_XIU_SIZE sizes the Xiu rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = n + 1

      return
      end
      subroutine epn_glg_monomial_integral ( n, expon, alpha, value )

c*********************************************************************72
c
cc EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
c
c  Discussion:
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c    value = integral ( EPN ) 
c      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer EXPON(N), the exponents.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      double precision alpha
      integer expon(n)
      integer i
      double precision value
      double precision value2

      value = 1.0D+00
      do i = 1, n
        call ep1_glg_monomial_integral ( expon(i), alpha, value2 )
        value = value * value2
      end do

      return
      end
      subroutine ep1_lag_monomial_integral ( expon, value )

c*********************************************************************72
c
cc EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
c
c  Discussion:
c
c    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
c
c      w(x) = exp ( - x )
c
c    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer EXPON, the exponent.
c    0 <= EXPON.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer expon
      double precision r8_factorial
      double precision value

      value = r8_factorial ( expon )

      return
      end
      subroutine epn_glg_00_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, volume )
      volume = volume ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = 1.0D+00
      end do
      w(k) = volume

      return
      end
      subroutine epn_glg_00_1_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine epn_glg_01_1 ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      integer expon
      integer i
      integer j
      integer k
      double precision value1
      double precision value2
      double precision volume
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, value1 )
      volume = value1 ** n

      expon = 1
      call ep1_glg_monomial_integral ( expon, alpha, value2 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = value2 / value1
      end do
      w(k) = volume

      return
      end
      subroutine epn_glg_01_1_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = 1

      return
      end
      subroutine epn_glg_02_xiu ( n, alpha, o, x, w )

c*********************************************************************72
c
cc EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision alpha
      double precision arg
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_02_XIU - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = - 1.0D+00
      delta0 = alpha + 1.0D+00
      c1 = - alpha - 1.0D+00

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      expon = 0
      call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
      volume = volume_1d ** n

      do i = 1, o
        w(i) = volume / dble ( o )
      end do

      return
      end
      subroutine epn_glg_02_xiu_size ( n, alpha, o )

c*********************************************************************72
c
cc EPN_GLG_02_XIU_SIZE sizes the Xiu rule for region EPN_GLG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, integer O, the order.
c
      implicit none

      double precision alpha
      integer n
      integer o

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
        write ( *, '(a)' ) '  ALPHA <= -1.0'
        stop
      end if

      o = n + 1

      return
      end
      subroutine epn_glg_monomial_integral ( n, expon, alpha, value )

c*********************************************************************72
c
cc EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
c
c  Discussion:
c
c    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
c    Laguerre weight function:
c
c      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
c
c    value = integral ( EPN ) 
c      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer EXPON(N), the exponents.
c
c    Input, double precision ALPHA, the exponent of X in the weight function.
c    -1.0 < ALPHA.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      double precision alpha
      integer expon(n)
      integer i
      double precision value
      double precision value2

      value = 1.0D+00
      do i = 1, n
        call ep1_glg_monomial_integral ( expon(i), alpha, value2 )
        value = value * value2
      end do

      return
      end
      subroutine epn_lag_00_1 ( n, o, x, w )

c*********************************************************************72
c
cc EPN_LAG_00_1 implements the "midpoint rule" for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      integer expon
      integer i
      integer j
      integer k
      double precision volume
      double precision w(o)
      double precision x(n,o)

      expon = 0
      call ep1_lag_monomial_integral ( expon, volume )
      volume = volume ** n

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = 1.0D+00
      end do
      w(k) = volume

      return
      end
      subroutine epn_lag_00_1_size ( n, o )

c*********************************************************************72
c
cc EPN_LAG_00_1_SIZE sizes the midpoint rule for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 0.
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 1

      return
      end
      subroutine epn_lag_01_1 ( n, o, x, w )

c*********************************************************************72
c
cc EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      integer expon
      integer i
      integer j
      integer k
      double precision value1
      double precision value2
      double precision volume
      double precision w(o)
      double precision x(n,o)

      expon = 0
      call ep1_lag_monomial_integral ( expon, value1 )
      volume = value1 ** n

      expon = 1
      call ep1_lag_monomial_integral ( expon, value2 )

      do j = 1, o
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      k = 0
c
c  1 point.
c
      k = k + 1
      do i = 1, n
        x(i,k) = value2 / value1
      end do
      w(k) = volume

      return
      end
      subroutine epn_lag_01_1_size ( n, o )

c*********************************************************************72
c
cc EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order O = 1.
c
c    The rule has precision P = 1.
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = 1

      return
      end
      subroutine epn_lag_02_xiu ( n, o, x, w )

c*********************************************************************72
c
cc EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision arg
      double precision c1
      double precision coef
      double precision delta0
      integer expon
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      do j = 1, o

        i = 0 
        do r = 1, n / 2

          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )

          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do

      gamma0 = - 1.0D+00
      delta0 = 1.0D+00
      c1 = - 1.0D+00

      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do

      expon = 0
      call ep1_lag_monomial_integral ( expon, volume_1d )
      volume = volume_1d ** n

      do i = 1, o
        w(i) = volume / dble ( o )
      end do

      return
      end
      subroutine epn_lag_02_xiu_size ( n, o )

c*********************************************************************72
c
cc EPN_LAG_02_XIU_SIZE sizes the Xiu rule for region EPN_LAG.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n + 1

      return
      end
      subroutine epn_lag_monomial_integral ( n, expon, value )

c*********************************************************************72
c
cc EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
c
c  Discussion:
c
c    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
c    or Laguerre weight function:
c
c      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
c
c    value = integral ( EPN ) 
c      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer EXPON(N), the exponents.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      integer n

      integer expon(n)
      integer i
      double precision value
      double precision value2

      value = 1.0D+00
      do i = 1, n
        call ep1_lag_monomial_integral ( expon(i), value2 )
        value = value * value2
      end do

      return
      end
      subroutine gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )

c*********************************************************************72
c
cc GW_02_XIU implements the Golub-Welsch version of the Xiu rule.
c
c  Discussion:
c
c    The rule has order 
c
c      O = N + 1.
c
c    The rule has precision P = 2.
c
c    It is assumed that the integral is over an N-dimensional region,
c    and has the form
c
c      Integral f(x) w(x) dx
c
c    where w(x) is separable into identical and independent components:
c
c      w(x) = v(x1) * v(x2) * ... * v(xn)
c
c    Associated with the weight function v(x), we assume there is a
c    family of orthogonal polynomials satisfying a three-term recurrence
c    of the form:
c
c      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)
c
c    with P(0,x) = 1, and P(-1,x) = 0.
c
c    This routine can construct the desired quadrature rule by knowing
c    the values of C1, used in the definition of P2, the values
c    GAMMA0 = 1/A0 and DELTA0 = - B0/A0, for which it is the case that
c    P(1,X) = GAMMA0 * X + DELTA0, and the value of VOLUME_1D, that is,
c    the 1D integral of v(x) over the region.
c
c    Note the values for the following standard polynomial families:
c
c    Chebyshev Type 1
c      V(X) =      1 / sqrt ( 1 - X^2 )
c      Interval =  [-1,+1]
c      GAMMA0 =    1.0
c      DELTA0 =    0.0
c      C1 =        1/2
c      VOLUME_1D = PI
c
c    Chebyshev Type 2
c      V(X) =      sqrt ( 1 - X^2 )
c      Interval =  [-1,+1]
c      GAMMA0 =    2.0
c      DELTA0 =    0.0
c      C1 =        1/2
c      VOLUME_1D = PI / 2
c
c    Gegenbauer
c      V(X) =      ( 1 - X^2 )^A
c      Interval =  [-1,+1]
c      GAMMA0 =    2 * A + 1
c      DELTA0 =    0.0
c      C1 =        ( 2 * A + 1 ) / ( 2 A + 3 )
c      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
c
c    Gegenbauer* (Removes singularity at ALPHA = -0.5):
c      V(X) =      ( 1 - X^2 )^A
c      Interval =  [-1,+1]
c      GAMMA0 =    1
c      DELTA0 =    0.0
c      C1 =        1 / ( 2 A + 3 )
c      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
c
c    Generalized Hermite
c      V(X) = |x|^A exp ( - x^2 )
c      Interval = (-oo,+oo)
c      GAMMA0 =    2
c      DELTA0 =    0
c      C1 =        2+2A
c      VOLUME_1D = Gamma((A+1)/2)
c
c    Generalized Laguerre
c      V(X) =       x^A exp ( - x )
c      Interval =  [0,+oo)
c      GAMMA0 =    -1.0
c      DELTA0 =     A+1.0
c      C1 =        -A-1.0
c      VOLUME_1D =  Gamma(A+1)
c
c    Hermite (physicist)
c      V(X) =      exp ( - x^2 )
c      Interval =  (-oo,+oo)
c      GAMMA0 =    2.0
c      DELTA0 =    0.0
c      C1 =        1.0
c      VOLUME_1D = sqrt ( PI )
c
c    Hermite (probabilist)
c      V(X) =      exp ( - x^2 / 2 )
c      Interval =  (-oo,+oo)
c      GAMMA0 =    1.0
c      DELTA0 =    0.0
c      C1 =        1.0
c      VOLUME_1D = sqrt ( 2 PI )
c
c    Jacobi
c      V(X) =      (1-x)^A (1+x)^B
c      Interval =  [-1,+1]
c      GAMMA0 =    (A+B+2)/2  
c      DELTA0 =    (A-B)/2
c      C1 =        2(A+1)(B+1)/(A+B+3)/(A+B+2)
c      VOLUME_1D = 2^(A+B+1) * Gamma(A+1) * Gamma(B+1) / ( A+B+1) / Gamma(A+B+1)
c
c    Laguerre
c      V(X) =       exp ( - x )
c      Interval =  [0,+oo)
c      GAMMA0 =    -1.0
c      DELTA0 =     1.0
c      C1 =        -1.0
c      VOLUME_1D =  1.0
c
c    Legendre
c      V(X) =      1.0
c      Interval =  [-1,+1]
c      GAMMA0 =    1.0
c      DELTA0 =    0.0
c      C1 =        1/3
c      VOLUME_1D = 2.0
c                                  
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer O, the order.
c
c    Input, double precision GAMMA0, the ratio 1 / A0.
c
c    Input, double precision DELTA0, the ratio B0 / A0.
c
c    Input, double precision C1, the coefficient of P(0,X) in 
c    the definition of P(2,X).
c
c    Input, double precision VOLUME_1D, the 1D integral of V(X).
c
c    Output, double precision X(N,O), the abscissas.
c
c    Output, double precision W(O), the weights.
c
      implicit none

      integer n
      integer o

      double precision arg
      double precision c1
      double precision delta0
      double precision gamma0
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer r
      double precision r8_mop
      double precision volume_1d
      double precision w(o)
      double precision x(n,o)

      do j = 1, o

        i = 0
        do r = 1, ( n / 2 )
          arg = dble ( 2 * r * ( j - 1 ) ) * pi / dble ( n + 1 )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg )
          i = i + 1
          x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )
        end do

        if ( i .lt. n ) then
          i = i + 1
          x(i,j) = r8_mop ( j - 1 )
        end if

      end do
c
c  Adjust for the GW rule.
c
      do j = 1, o
        do i = 1, n
          x(i,j) = ( sqrt ( gamma0 * c1 ) * x(i,j) - delta0 ) / gamma0
        end do
      end do
c
c  The weights are equal.
c
      do j = 1, o
        w(j) = volume_1d ** n  / dble ( o )
      end do

      return
      end
      subroutine gw_02_xiu_size ( n, o )

c*********************************************************************72
c
cc GW_02_XIU_SIZE sizes the Golub Welsch version of the Xiu rule.
c
c  Discussion:
c
c    The rule has order O = N + 1;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical integration formulas of degree two,
c    Applied Numerical Mathematics,
c    Volume 58, 2008, pages 1515-1520.
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, integer O, the order.
c
      implicit none

      integer n
      integer o

      o = n + 1

      return
      end
      function hexagon_area_2d ( r )

c*********************************************************************72
c
cc HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
c
c  Discussion:
c
c    The formula for the area only requires the radius, and does
c    not depend on the location of the center, or the orientation
c    of the hexagon.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the hexagon.
c
c    Output, double precision HEXAGON_AREA_2D, the area of the hexagon.
c
      implicit none

      double precision hexagon_area_2d
      double precision hexagon_unit_area_2d
      double precision r

      hexagon_area_2d = r * r * hexagon_unit_area_2d ( )

      return
      end
      subroutine hexagon_sum ( func, center, r, order, xtab, ytab, 
     &  weight, result )

c*********************************************************************72
c
cc HEXAGON_SUM applies a quadrature rule inside a hexagon in 2D.
c
c  Discussion:
c
c    The input quadrature rule is assumed to be defined for a unit hexagon.
c
c    The input quadrature rule may be defined by calling HEXAGON_UNIT_SET.
c
c  Integration region:
c
c    The definition is given in terms of THETA, the angle in degrees of the
c    vector (X-CENTER(1),Y-CENTER(2)).  The following six conditions apply,
c    respectively, between the bracketing values of THETA of 0, 60, 120, 
c    180, 240, 300, and 360.
c
c      0 <= Y-CENTER(2) <= -SQRT(3) * (X-CENTER(1)) + R * SQRT(3)
c      0 <= Y-CENTER(2) <=                     R * SQRT(3)/2
c      0 <= Y-CENTER(2) <=  SQRT(3) * (X-CENTER(1)) + R * SQRT(3) 
c      -SQRT(3) * (X-CENTER(1)) - R * SQRT(3)	<= Y-CENTER(2) <= 0
c                        - R * SQRT(3)/2 <= Y-CENTER(2) <= 0
c       SQRT(3) * (X-CENTER(1)) - R * SQRT(3)   <= Y-CENTER(2) <= 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of two variables which is to be integrated,
c    of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the hexagon.
c
c    Input, double precision R, the radius of the hexagon.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer order

      double precision center(dim_num)
      double precision func
      external func
      double precision hexagon_area_2d
      integer i
      double precision quad
      double precision r
      double precision result
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y
      double precision ytab(order)

      quad = 0.0D+00
      do i = 1, order
        x = center(1) + r * xtab(i)
        y = center(2) + r * ytab(i)
        quad = quad + weight(i) * func ( x, y )
      end do

      volume = hexagon_area_2d ( r )
      result = quad * volume

      return
      end
      function hexagon_unit_area_2d ( )

c*********************************************************************72
c
cc HEXAGON_UNIT_AREA_2D returns the area of the unit regular hexagon in 2D.
c
c  Integration region:
c
c    The definition is given in terms of THETA, the angle in degrees of the
c    vector (X,Y).  The following six conditions apply, respectively,
c    between the bracketing values of THETA of 0, 60, 120, 180, 240,
c    300, and 360.
c
c                              0 <= Y <= -SQRT(3) * X + SQRT(3)
c                              0 <= Y <=                 SQRT(3)/2
c                              0 <= Y <=  SQRT(3) * X + SQRT(3)
c      - SQRT(3) * X - SQRT(3)   <= Y <= 0
c                    - SQRT(3)/2 <= Y <= 0
c        SQRT(3) * X - SQRT(3)   <= Y <= 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision HEXAGON_UNIT_AREA_2D, the area of the hexagon.
c
      implicit none

      double precision hexagon_unit_area_2d

      hexagon_unit_area_2d = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

      return
      end
      subroutine hexagon_unit_set ( rule, order, xtab, ytab, weight )

c*********************************************************************72
c
cc HEXAGON_UNIT_SET sets a quadrature rule inside the unit hexagon in 2D.
c
c  Integration region:
c
c    The definition is given in terms of THETA, the angle in degrees of the
c    vector (X,Y).  The following six conditions apply, respectively,
c    between the bracketing values of THETA of 0, 60, 120, 180, 240,
c    300, and 360.
c
c                              0 <= Y <= -SQRT(3) * X + SQRT(3)
c                              0 <= Y <=                 SQRT(3)/2
c                              0 <= Y <=  SQRT(3) * X + SQRT(3)
c       -SQRT(3) * X - SQRT(3)   <= Y <= 0
c                    - SQRT(3)/2 <= Y <= 0
c        SQRT(3) * X - SQRT(3)   <= Y <= 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 December 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c      1, 1 point,  degree 1;
c      2, 4 points, degree 3;
c      3, 7 points, degree 3;
c      4, 7 points, degree 5;
c
c    Input, integer ORDER, the order of the desired rule.
c
c    Output, double precision XTAB(*), YTAB(*), the abscissas of the rule.
c
c    Output, double precision WEIGHT(*), the ORDER weights of the rule.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      integer i
      integer rule
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)
      double precision z

      if ( rule .eq. 1 ) then

        xtab(1) = 0.0D+00
        ytab(1) = 0.0D+00
        weight(1) = 1.0D+00
c
c  Stroud rule H2:3-1.
c
      else if ( rule .eq. 2 ) then

        a = sqrt ( 5.0D+00 / 12.0D+00 )
        b = 1.0D+00 / 4.0D+00
        z = 0.0D+00

        xtab(1) =  a
        xtab(2) = -a
        xtab(3) =  z
        xtab(4) =  z
        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  a
        ytab(4) = -a
        weight(1) = b
        weight(2) = b
        weight(3) = b
        weight(4) = b
c
c  Stroud rule H2:3-2.
c
      else if ( rule .eq. 3 ) then

        a = sqrt ( 3.0D+00 ) / 2.0D+00
        b =  0.5D+00
        c =  1.0D+00
        d =  5.0D+00 / 72.0D+00
        e = 42.0D+00 / 72.0D+00
        z =  0.0D+00

        xtab(1) =  z
        xtab(2) =  c
        xtab(3) = -c
        xtab(4) =  b
        xtab(5) = -b
        xtab(6) =  b
        xtab(7) = -b
        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  a
        ytab(5) =  a
        ytab(6) = -a
        ytab(7) = -a
        weight(1) = e
        weight(2) = d
        weight(3) = d
        weight(4) = d
        weight(5) = d
        weight(6) = d
        weight(7) = d
c
c  Stroud rule H2:5-1.
c
      else if ( rule .eq. 4 ) then

        a = sqrt ( 14.0D+00 ) / 5.0D+00
        b = sqrt ( 14.0D+00 ) / 10.0D+00
        c = sqrt ( 42.0D+00 ) / 10.0D+00
        d = 125.0D+00 / 1008.0D+00
        e = 258.0D+00 / 1008.0D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  a
        xtab(3) = -a
        xtab(4) =  b
        xtab(5) = -b
        xtab(6) =  b
        xtab(7) = -b
        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  c
        ytab(5) =  c
        ytab(6) = -c
        ytab(7) = -c
        weight(1) = e
        weight(2) = d
        weight(3) = d
        weight(4) = d
        weight(5) = d
        weight(6) = d
        weight(7) = d

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HEXAGON_UNIT_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal input value of RULE = ', rule
        stop

      end if

      return
      end
      subroutine hexagon_unit_size ( rule, order )

c*********************************************************************72
c
cc HEXAGON_UNIT_SIZE sizes a quadrature rule inside the unit hexagon in 2D.
c
c  Integration region:
c
c    The definition is given in terms of THETA, the angle in degrees of the
c    vector (X,Y).  The following six conditions apply, respectively,
c    between the bracketing values of THETA of 0, 60, 120, 180, 240,
c    300, and 360.
c
c                              0 <= Y <= -SQRT(3) * X + SQRT(3)
c                              0 <= Y <=                 SQRT(3)/2
c                              0 <= Y <=  SQRT(3) * X + SQRT(3)
c       -SQRT(3) * X - SQRT(3)   <= Y <= 0
c                    - SQRT(3)/2 <= Y <= 0
c        SQRT(3) * X - SQRT(3)   <= Y <= 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule desired.
c      1, 1 point,  degree 1;
c      2, 4 points, degree 3;
c      3, 7 points, degree 3;
c      4, 7 points, degree 5;
c
c    Output, integer ORDER, the order of the desired rule.
c    If RULE is not legal, then ORDER is returned as -1.
c
      implicit none

      integer order
      integer rule

      if ( rule .eq. 1 ) then

        order = 1
c
c  Stroud rule H2:3-1.
c
      else if ( rule .eq. 2 ) then

        order = 4
c
c  Stroud rule H2:3-2.
c
      else if ( rule .eq. 3 ) then

        order = 7
c
c  Stroud rule H2:5-1.
c
      else if ( rule .eq. 4 ) then

        order = 7

      else

        order = -1

      end if

      return
      end
      function i4_factorial ( n )

c*********************************************************************72
c
cc I4_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c    0 <= N <= 13 is required.
c
c    Output, integer I4_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer i4_factorial
      integer n

      i4_factorial = 1

      if ( 13 .lt. n ) then
        i4_factorial = - 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
        write ( *, '(a)' )
     &  '  I4_FACTORIAL(N) cannot be computed as an integer'
        write ( *, '(a)' ) '  for 13 < N.'
        write ( *, '(a,i8)' ) '  Input value N = ', n
        stop
      end if

      do i = 1, n
        i4_factorial = i4_factorial * i
      end do

      return
      end
      function i4_factorial2 ( n )

c*********************************************************************72
c
cc I4_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial 
c    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
c
c    Output, integer I4_FACTORIAL2, the value of N!!.
c
      implicit none

      integer i4_factorial2
      integer n
      integer n_copy

      if ( n .lt. 1 ) then
        i4_factorial2 = 1
        return
      end if

      n_copy = n
      i4_factorial2 = 1

10    continue

      if ( 1 .lt. n_copy ) then
        i4_factorial2 = i4_factorial2 * n_copy
        n_copy = n_copy - 2
        go to 10
      end if

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    SUM function:
c
c      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
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
c    Output, integer I4VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_sum

      i4vec_sum = 0

      do i = 1, n
        i4vec_sum = i4vec_sum + a(i)
      end do

      return
      end
      subroutine ksub_next2 ( n, k, a, in, iout )

c*********************************************************************72
c
cc KSUB_NEXT2 generates the subsets of size K from a set of size N.
c
c  Discussion:
c
c    This routine uses the revolving door method.  It has no "memory".
c    It simply calculates the successor of the input set,
c    and will start from the beginning after the last set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
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
c    Input, integer N, the size of the set from which subsets are drawn.
c    N must be positive.
c
c    Input, integer K, the size of the desired subset.  K must be
c    between 0 and N.
c
c    Input/output, integer A(K).  On input, the user must
c    supply a subset of size K in A.  That is, A must
c    contain K unique numbers, in order, between 1 and N.  On
c    output, A(I) is the I-th element of the output subset.
c    The output array is also in sorted order.
c
c    Output, integer IN, the element of the output subset which
c    was not in the input set.  Each new subset differs from the
c    last one by adding one element and deleting another.
c
c    Output, integer IOUT, the element of the input subset which
c    is not in the output subset.
c
      implicit none

      integer k

      integer a(k)
      integer in
      integer iout
      integer j
      integer m
      integer n

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a)' ) '  but 0 .lt. N is required!'
        stop
      end if

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if

      j = 0

10    continue

        if ( 0 .lt. j .or. mod ( k, 2 ) .eq. 0 ) then

          j = j + 1

          if ( k .lt. j ) then
            a(k) = k
            in = k
            iout = n
            return
          end if

          if ( a(j) .ne. j ) then

            iout = a(j)
            in = iout - 1
            a(j) = in

            if ( j .ne. 1 ) then
              in = j - 1
              a(j-1) = in
            end if

            return

          end if

        end if

        j = j + 1
        m = n

        if ( j .lt. k ) then
          m = a(j+1) - 1
        end if

        if ( m .ne. a(j) ) then
          go to 20
        end if

      go to 10

20    continue

      in = a(j) + 1
      a(j) = in
      iout = in - 1

      if ( j .ne. 1 ) then
        a(j-1) = iout
        iout = j - 1
      end if

      return
      end
      subroutine legendre_set ( n, x, w )

c*********************************************************************72
c
cc LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
c
c  Discussion:
c
c    The integral:
c
c      integral ( -1 <= x <= 1 ) f(x) dx
c
c    The quadrature rule:
c
c      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
c
c    The quadrature rule is exact for polynomials through degree 2*N-1.
c
c    The abscissas are the zeroes of the Legendre polynomial P(N)(X).
c
c    Mathematica can compute the abscissas and weights of a Gauss-Legendre
c    rule of order N for the interval [A,B] with P digits of precision
c    by the commands:
c
c    Needs["NumericalDifferentialEquationAnalysis`"]
c    GaussianQuadratureWeights [ n, a, b, p ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Vladimir Krylov,
c    Approximate Calculation of Integrals,
c    Dover, 2006,
c    ISBN: 0486445798,
c    LC: QA311.K713.
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input, integer N, the order.
c    N must be between 1 and 33 or 63/64/65, 127/128/129,
c    255/256/257.
c
c    Output, double precision X(N), the abscissas.
c
c    Output, double precision W(N), the weights.
c
      implicit none

      integer n

      double precision w(n)
      double precision x(n)

      if ( n .eq. 1 ) then

        x(1) = 0.000000000000000000000000000000D+00

        w(1) = 2.000000000000000000000000000000D+00

      else if ( n .eq. 2 ) then

        x(1) = -0.577350269189625764509148780502D+00
        x(2) = 0.577350269189625764509148780502D+00

        w(1) = 1.000000000000000000000000000000D+00
        w(2) = 1.000000000000000000000000000000D+00

      else if ( n .eq. 3 ) then

        x(1) = -0.774596669241483377035853079956D+00
        x(2) = 0.000000000000000000000000000000D+00
        x(3) = 0.774596669241483377035853079956D+00

        w(1) = 0.555555555555555555555555555556D+00
        w(2) = 0.888888888888888888888888888889D+00
        w(3) = 0.555555555555555555555555555556D+00

      else if ( n .eq. 4 ) then

        x(1) = -0.861136311594052575223946488893D+00
        x(2) = -0.339981043584856264802665759103D+00
        x(3) = 0.339981043584856264802665759103D+00
        x(4) = 0.861136311594052575223946488893D+00

        w(1) = 0.347854845137453857373063949222D+00
        w(2) = 0.652145154862546142626936050778D+00
        w(3) = 0.652145154862546142626936050778D+00
        w(4) = 0.347854845137453857373063949222D+00

      else if ( n .eq. 5 ) then

        x(1) = -0.906179845938663992797626878299D+00
        x(2) = -0.538469310105683091036314420700D+00
        x(3) = 0.000000000000000000000000000000D+00
        x(4) = 0.538469310105683091036314420700D+00
        x(5) = 0.906179845938663992797626878299D+00

        w(1) = 0.236926885056189087514264040720D+00
        w(2) = 0.478628670499366468041291514836D+00
        w(3) = 0.568888888888888888888888888889D+00
        w(4) = 0.478628670499366468041291514836D+00
        w(5) = 0.236926885056189087514264040720D+00

      else if ( n .eq. 6 ) then

        x(1) = -0.932469514203152027812301554494D+00
        x(2) = -0.661209386466264513661399595020D+00
        x(3) = -0.238619186083196908630501721681D+00
        x(4) = 0.238619186083196908630501721681D+00
        x(5) = 0.661209386466264513661399595020D+00
        x(6) = 0.932469514203152027812301554494D+00

        w(1) = 0.171324492379170345040296142173D+00
        w(2) = 0.360761573048138607569833513838D+00
        w(3) = 0.467913934572691047389870343990D+00
        w(4) = 0.467913934572691047389870343990D+00
        w(5) = 0.360761573048138607569833513838D+00
        w(6) = 0.171324492379170345040296142173D+00

      else if ( n .eq. 7 ) then

        x(1) = -0.949107912342758524526189684048D+00
        x(2) = -0.741531185599394439863864773281D+00
        x(3) = -0.405845151377397166906606412077D+00
        x(4) = 0.000000000000000000000000000000D+00
        x(5) = 0.405845151377397166906606412077D+00
        x(6) = 0.741531185599394439863864773281D+00
        x(7) = 0.949107912342758524526189684048D+00

        w(1) = 0.129484966168869693270611432679D+00
        w(2) = 0.279705391489276667901467771424D+00
        w(3) = 0.381830050505118944950369775489D+00
        w(4) = 0.417959183673469387755102040816D+00
        w(5) = 0.381830050505118944950369775489D+00
        w(6) = 0.279705391489276667901467771424D+00
        w(7) = 0.129484966168869693270611432679D+00

      else if ( n .eq. 8 ) then

        x(1) = -0.960289856497536231683560868569D+00
        x(2) = -0.796666477413626739591553936476D+00
        x(3) = -0.525532409916328985817739049189D+00
        x(4) = -0.183434642495649804939476142360D+00
        x(5) = 0.183434642495649804939476142360D+00
        x(6) = 0.525532409916328985817739049189D+00
        x(7) = 0.796666477413626739591553936476D+00
        x(8) = 0.960289856497536231683560868569D+00

        w(1) = 0.101228536290376259152531354310D+00
        w(2) = 0.222381034453374470544355994426D+00
        w(3) = 0.313706645877887287337962201987D+00
        w(4) = 0.362683783378361982965150449277D+00
        w(5) = 0.362683783378361982965150449277D+00
        w(6) = 0.313706645877887287337962201987D+00
        w(7) = 0.222381034453374470544355994426D+00
        w(8) = 0.101228536290376259152531354310D+00

      else if ( n .eq. 9 ) then

        x(1) = -0.968160239507626089835576203D+00
        x(2) = -0.836031107326635794299429788D+00
        x(3) = -0.613371432700590397308702039D+00
        x(4) = -0.324253423403808929038538015D+00
        x(5) = 0.000000000000000000000000000D+00
        x(6) = 0.324253423403808929038538015D+00
        x(7) = 0.613371432700590397308702039D+00
        x(8) = 0.836031107326635794299429788D+00
        x(9) = 0.968160239507626089835576203D+00

        w(1) = 0.081274388361574411971892158111D+00
        w(2) = 0.18064816069485740405847203124D+00
        w(3) = 0.26061069640293546231874286942D+00
        w(4) = 0.31234707704000284006863040658D+00
        w(5) = 0.33023935500125976316452506929D+00
        w(6) = 0.31234707704000284006863040658D+00
        w(7) = 0.26061069640293546231874286942D+00
        w(8) = 0.18064816069485740405847203124D+00
        w(9) = 0.081274388361574411971892158111D+00

      else if ( n .eq. 10 ) then

        x(1) = -0.973906528517171720077964012D+00
        x(2) = -0.865063366688984510732096688D+00
        x(3) = -0.679409568299024406234327365D+00
        x(4) = -0.433395394129247190799265943D+00
        x(5) = -0.148874338981631210884826001D+00
        x(6) = 0.148874338981631210884826001D+00
        x(7) = 0.433395394129247190799265943D+00
        x(8) = 0.679409568299024406234327365D+00
        x(9) = 0.865063366688984510732096688D+00
        x(10) = 0.973906528517171720077964012D+00

        w(1) = 0.066671344308688137593568809893D+00
        w(2) = 0.14945134915058059314577633966D+00
        w(3) = 0.21908636251598204399553493423D+00
        w(4) = 0.26926671930999635509122692157D+00
        w(5) = 0.29552422471475287017389299465D+00
        w(6) = 0.29552422471475287017389299465D+00
        w(7) = 0.26926671930999635509122692157D+00
        w(8) = 0.21908636251598204399553493423D+00
        w(9) = 0.14945134915058059314577633966D+00
        w(10) = 0.066671344308688137593568809893D+00

      else if ( n .eq. 11 ) then

        x(1) = -0.978228658146056992803938001D+00
        x(2) = -0.887062599768095299075157769D+00
        x(3) = -0.730152005574049324093416252D+00
        x(4) = -0.519096129206811815925725669D+00
        x(5) = -0.269543155952344972331531985D+00
        x(6) = 0.000000000000000000000000000D+00
        x(7) = 0.269543155952344972331531985D+00
        x(8) = 0.519096129206811815925725669D+00
        x(9) = 0.730152005574049324093416252D+00
        x(10) = 0.887062599768095299075157769D+00
        x(11) = 0.978228658146056992803938001D+00

        w(1) = 0.055668567116173666482753720443D+00
        w(2) = 0.12558036946490462463469429922D+00
        w(3) = 0.18629021092773425142609764143D+00
        w(4) = 0.23319376459199047991852370484D+00
        w(5) = 0.26280454451024666218068886989D+00
        w(6) = 0.27292508677790063071448352834D+00
        w(7) = 0.26280454451024666218068886989D+00
        w(8) = 0.23319376459199047991852370484D+00
        w(9) = 0.18629021092773425142609764143D+00
        w(10) = 0.12558036946490462463469429922D+00
        w(11) = 0.055668567116173666482753720443D+00

      else if ( n .eq. 12 ) then

        x(1) = -0.981560634246719250690549090D+00
        x(2) = -0.904117256370474856678465866D+00
        x(3) = -0.769902674194304687036893833D+00
        x(4) = -0.587317954286617447296702419D+00
        x(5) = -0.367831498998180193752691537D+00
        x(6) = -0.125233408511468915472441369D+00
        x(7) = 0.125233408511468915472441369D+00
        x(8) = 0.367831498998180193752691537D+00
        x(9) = 0.587317954286617447296702419D+00
        x(10) = 0.769902674194304687036893833D+00
        x(11) = 0.904117256370474856678465866D+00
        x(12) = 0.981560634246719250690549090D+00

        w(1) = 0.047175336386511827194615961485D+00
        w(2) = 0.10693932599531843096025471819D+00
        w(3) = 0.16007832854334622633465252954D+00
        w(4) = 0.20316742672306592174906445581D+00
        w(5) = 0.23349253653835480876084989892D+00
        w(6) = 0.24914704581340278500056243604D+00
        w(7) = 0.24914704581340278500056243604D+00
        w(8) = 0.23349253653835480876084989892D+00
        w(9) = 0.20316742672306592174906445581D+00
        w(10) = 0.16007832854334622633465252954D+00
        w(11) = 0.10693932599531843096025471819D+00
        w(12) = 0.047175336386511827194615961485D+00

      else if ( n .eq. 13 ) then

        x(1) = -0.984183054718588149472829449D+00
        x(2) = -0.917598399222977965206547837D+00
        x(3) = -0.801578090733309912794206490D+00
        x(4) = -0.642349339440340220643984607D+00
        x(5) = -0.448492751036446852877912852D+00
        x(6) = -0.230458315955134794065528121D+00
        x(7) = 0.000000000000000000000000000D+00
        x(8) = 0.230458315955134794065528121D+00
        x(9) = 0.448492751036446852877912852D+00
        x(10) = 0.642349339440340220643984607D+00
        x(11) = 0.80157809073330991279420649D+00
        x(12) = 0.91759839922297796520654784D+00
        x(13) = 0.98418305471858814947282945D+00

        w(1) = 0.040484004765315879520021592201D+00
        w(2) = 0.092121499837728447914421775954D+00
        w(3) = 0.13887351021978723846360177687D+00
        w(4) = 0.17814598076194573828004669200D+00
        w(5) = 0.20781604753688850231252321931D+00
        w(6) = 0.22628318026289723841209018604D+00
        w(7) = 0.23255155323087391019458951527D+00
        w(8) = 0.22628318026289723841209018604D+00
        w(9) = 0.20781604753688850231252321931D+00
        w(10) = 0.17814598076194573828004669200D+00
        w(11) = 0.13887351021978723846360177687D+00
        w(12) = 0.092121499837728447914421775954D+00
        w(13) = 0.040484004765315879520021592201D+00

      else if ( n .eq. 14 ) then

        x(1) = -0.986283808696812338841597267D+00
        x(2) = -0.928434883663573517336391139D+00
        x(3) = -0.827201315069764993189794743D+00
        x(4) = -0.687292904811685470148019803D+00
        x(5) = -0.515248636358154091965290719D+00
        x(6) = -0.319112368927889760435671824D+00
        x(7) = -0.108054948707343662066244650D+00
        x(8) = 0.108054948707343662066244650D+00
        x(9) = 0.31911236892788976043567182D+00
        x(10) = 0.51524863635815409196529072D+00
        x(11) = 0.68729290481168547014801980D+00
        x(12) = 0.82720131506976499318979474D+00
        x(13) = 0.92843488366357351733639114D+00
        x(14) = 0.98628380869681233884159727D+00

        w(1) = 0.035119460331751863031832876138D+00
        w(2) = 0.08015808715976020980563327706D+00
        w(3) = 0.12151857068790318468941480907D+00
        w(4) = 0.15720316715819353456960193862D+00
        w(5) = 0.18553839747793781374171659013D+00
        w(6) = 0.20519846372129560396592406566D+00
        w(7) = 0.21526385346315779019587644332D+00
        w(8) = 0.21526385346315779019587644332D+00
        w(9) = 0.20519846372129560396592406566D+00
        w(10) = 0.18553839747793781374171659013D+00
        w(11) = 0.15720316715819353456960193862D+00
        w(12) = 0.12151857068790318468941480907D+00
        w(13) = 0.08015808715976020980563327706D+00
        w(14) = 0.035119460331751863031832876138D+00

      else if ( n .eq. 15 ) then

        x(1) = -0.987992518020485428489565719D+00
        x(2) = -0.937273392400705904307758948D+00
        x(3) = -0.848206583410427216200648321D+00
        x(4) = -0.724417731360170047416186055D+00
        x(5) = -0.570972172608538847537226737D+00
        x(6) = -0.394151347077563369897207371D+00
        x(7) = -0.201194093997434522300628303D+00
        x(8) = 0.00000000000000000000000000D+00
        x(9) = 0.20119409399743452230062830D+00
        x(10) = 0.39415134707756336989720737D+00
        x(11) = 0.57097217260853884753722674D+00
        x(12) = 0.72441773136017004741618605D+00
        x(13) = 0.84820658341042721620064832D+00
        x(14) = 0.93727339240070590430775895D+00
        x(15) = 0.98799251802048542848956572D+00

        w(1) = 0.030753241996117268354628393577D+00
        w(2) = 0.070366047488108124709267416451D+00
        w(3) = 0.107159220467171935011869546686D+00
        w(4) = 0.13957067792615431444780479451D+00
        w(5) = 0.16626920581699393355320086048D+00
        w(6) = 0.18616100001556221102680056187D+00
        w(7) = 0.19843148532711157645611832644D+00
        w(8) = 0.20257824192556127288062019997D+00
        w(9) = 0.19843148532711157645611832644D+00
        w(10) = 0.18616100001556221102680056187D+00
        w(11) = 0.16626920581699393355320086048D+00
        w(12) = 0.13957067792615431444780479451D+00
        w(13) = 0.107159220467171935011869546686D+00
        w(14) = 0.070366047488108124709267416451D+00
        w(15) = 0.030753241996117268354628393577D+00

      else if ( n .eq. 16 ) then

        x(1) = -0.989400934991649932596154173D+00
        x(2) = -0.944575023073232576077988416D+00
        x(3) = -0.865631202387831743880467898D+00
        x(4) = -0.755404408355003033895101195D+00
        x(5) = -0.617876244402643748446671764D+00
        x(6) = -0.458016777657227386342419443D+00
        x(7) = -0.281603550779258913230460501D+00
        x(8) = -0.09501250983763744018531934D+00
        x(9) = 0.09501250983763744018531934D+00
        x(10) = 0.28160355077925891323046050D+00
        x(11) = 0.45801677765722738634241944D+00
        x(12) = 0.61787624440264374844667176D+00
        x(13) = 0.75540440835500303389510119D+00
        x(14) = 0.86563120238783174388046790D+00
        x(15) = 0.94457502307323257607798842D+00
        x(16) = 0.98940093499164993259615417D+00

        w(1) = 0.027152459411754094851780572456D+00
        w(2) = 0.062253523938647892862843836994D+00
        w(3) = 0.09515851168249278480992510760D+00
        w(4) = 0.12462897125553387205247628219D+00
        w(5) = 0.14959598881657673208150173055D+00
        w(6) = 0.16915651939500253818931207903D+00
        w(7) = 0.18260341504492358886676366797D+00
        w(8) = 0.18945061045506849628539672321D+00
        w(9) = 0.18945061045506849628539672321D+00
        w(10) = 0.18260341504492358886676366797D+00
        w(11) = 0.16915651939500253818931207903D+00
        w(12) = 0.14959598881657673208150173055D+00
        w(13) = 0.12462897125553387205247628219D+00
        w(14) = 0.09515851168249278480992510760D+00
        w(15) = 0.062253523938647892862843836994D+00
        w(16) = 0.027152459411754094851780572456D+00

      else if ( n .eq. 17 ) then

        x(1) = -0.990575475314417335675434020D+00
        x(2) = -0.950675521768767761222716958D+00
        x(3) = -0.880239153726985902122955694D+00
        x(4) = -0.781514003896801406925230056D+00
        x(5) = -0.657671159216690765850302217D+00
        x(6) = -0.512690537086476967886246569D+00
        x(7) = -0.35123176345387631529718552D+00
        x(8) = -0.17848418149584785585067749D+00
        x(9) = 0.00000000000000000000000000D+00
        x(10) = 0.17848418149584785585067749D+00
        x(11) = 0.35123176345387631529718552D+00
        x(12) = 0.51269053708647696788624657D+00
        x(13) = 0.65767115921669076585030222D+00
        x(14) = 0.78151400389680140692523006D+00
        x(15) = 0.88023915372698590212295569D+00
        x(16) = 0.95067552176876776122271696D+00
        x(17) = 0.99057547531441733567543402D+00

        w(1) = 0.024148302868547931960110026288D+00
        w(2) = 0.055459529373987201129440165359D+00
        w(3) = 0.085036148317179180883535370191D+00
        w(4) = 0.111883847193403971094788385626D+00
        w(5) = 0.13513636846852547328631998170D+00
        w(6) = 0.15404576107681028808143159480D+00
        w(7) = 0.16800410215645004450997066379D+00
        w(8) = 0.17656270536699264632527099011D+00
        w(9) = 0.17944647035620652545826564426D+00
        w(10) = 0.17656270536699264632527099011D+00
        w(11) = 0.16800410215645004450997066379D+00
        w(12) = 0.15404576107681028808143159480D+00
        w(13) = 0.13513636846852547328631998170D+00
        w(14) = 0.111883847193403971094788385626D+00
        w(15) = 0.085036148317179180883535370191D+00
        w(16) = 0.055459529373987201129440165359D+00
        w(17) = 0.024148302868547931960110026288D+00

      else if ( n .eq. 18 ) then

        x(1) = -0.991565168420930946730016005D+00
        x(2) = -0.955823949571397755181195893D+00
        x(3) = -0.892602466497555739206060591D+00
        x(4) = -0.803704958972523115682417455D+00
        x(5) = -0.691687043060353207874891081D+00
        x(6) = -0.55977083107394753460787155D+00
        x(7) = -0.41175116146284264603593179D+00
        x(8) = -0.25188622569150550958897285D+00
        x(9) = -0.08477501304173530124226185D+00
        x(10) = 0.08477501304173530124226185D+00
        x(11) = 0.25188622569150550958897285D+00
        x(12) = 0.41175116146284264603593179D+00
        x(13) = 0.55977083107394753460787155D+00
        x(14) = 0.69168704306035320787489108D+00
        x(15) = 0.80370495897252311568241746D+00
        x(16) = 0.89260246649755573920606059D+00
        x(17) = 0.95582394957139775518119589D+00
        x(18) = 0.99156516842093094673001600D+00

        w(1) = 0.021616013526483310313342710266D+00
        w(2) = 0.049714548894969796453334946203D+00
        w(3) = 0.07642573025488905652912967762D+00
        w(4) = 0.10094204410628716556281398492D+00
        w(5) = 0.12255520671147846018451912680D+00
        w(6) = 0.14064291467065065120473130375D+00
        w(7) = 0.15468467512626524492541800384D+00
        w(8) = 0.16427648374583272298605377647D+00
        w(9) = 0.16914238296314359184065647013D+00
        w(10) = 0.16914238296314359184065647013D+00
        w(11) = 0.16427648374583272298605377647D+00
        w(12) = 0.15468467512626524492541800384D+00
        w(13) = 0.14064291467065065120473130375D+00
        w(14) = 0.12255520671147846018451912680D+00
        w(15) = 0.10094204410628716556281398492D+00
        w(16) = 0.07642573025488905652912967762D+00
        w(17) = 0.049714548894969796453334946203D+00
        w(18) = 0.021616013526483310313342710266D+00

      else if ( n .eq. 19 ) then

        x(1) = -0.992406843843584403189017670D+00
        x(2) = -0.960208152134830030852778841D+00
        x(3) = -0.903155903614817901642660929D+00
        x(4) = -0.822714656537142824978922487D+00
        x(5) = -0.72096617733522937861709586D+00
        x(6) = -0.60054530466168102346963816D+00
        x(7) = -0.46457074137596094571726715D+00
        x(8) = -0.31656409996362983199011733D+00
        x(9) = -0.16035864564022537586809612D+00
        x(10) = 0.00000000000000000000000000D+00
        x(11) = 0.16035864564022537586809612D+00
        x(12) = 0.31656409996362983199011733D+00
        x(13) = 0.46457074137596094571726715D+00
        x(14) = 0.60054530466168102346963816D+00
        x(15) = 0.72096617733522937861709586D+00
        x(16) = 0.82271465653714282497892249D+00
        x(17) = 0.90315590361481790164266093D+00
        x(18) = 0.96020815213483003085277884D+00
        x(19) = 0.99240684384358440318901767D+00

        w(1) = 0.019461788229726477036312041464D+00
        w(2) = 0.044814226765699600332838157402D+00
        w(3) = 0.069044542737641226580708258006D+00
        w(4) = 0.091490021622449999464462094124D+00
        w(5) = 0.111566645547333994716023901682D+00
        w(6) = 0.12875396253933622767551578486D+00
        w(7) = 0.14260670217360661177574610944D+00
        w(8) = 0.15276604206585966677885540090D+00
        w(9) = 0.15896884339395434764995643946D+00
        w(10) = 0.16105444984878369597916362532D+00
        w(11) = 0.15896884339395434764995643946D+00
        w(12) = 0.15276604206585966677885540090D+00
        w(13) = 0.14260670217360661177574610944D+00
        w(14) = 0.12875396253933622767551578486D+00
        w(15) = 0.111566645547333994716023901682D+00
        w(16) = 0.091490021622449999464462094124D+00
        w(17) = 0.069044542737641226580708258006D+00
        w(18) = 0.044814226765699600332838157402D+00
        w(19) = 0.019461788229726477036312041464D+00

      else if ( n .eq. 20 ) then

        x(1) = -0.993128599185094924786122388D+00
        x(2) = -0.963971927277913791267666131D+00
        x(3) = -0.912234428251325905867752441D+00
        x(4) = -0.83911697182221882339452906D+00
        x(5) = -0.74633190646015079261430507D+00
        x(6) = -0.63605368072651502545283670D+00
        x(7) = -0.51086700195082709800436405D+00
        x(8) = -0.37370608871541956067254818D+00
        x(9) = -0.22778585114164507808049620D+00
        x(10) = -0.07652652113349733375464041D+00
        x(11) = 0.07652652113349733375464041D+00
        x(12) = 0.22778585114164507808049620D+00
        x(13) = 0.37370608871541956067254818D+00
        x(14) = 0.51086700195082709800436405D+00
        x(15) = 0.63605368072651502545283670D+00
        x(16) = 0.74633190646015079261430507D+00
        x(17) = 0.83911697182221882339452906D+00
        x(18) = 0.91223442825132590586775244D+00
        x(19) = 0.96397192727791379126766613D+00
        x(20) = 0.99312859918509492478612239D+00

        w(1) = 0.017614007139152118311861962352D+00
        w(2) = 0.040601429800386941331039952275D+00
        w(3) = 0.062672048334109063569506535187D+00
        w(4) = 0.08327674157670474872475814322D+00
        w(5) = 0.10193011981724043503675013548D+00
        w(6) = 0.11819453196151841731237737771D+00
        w(7) = 0.13168863844917662689849449975D+00
        w(8) = 0.14209610931838205132929832507D+00
        w(9) = 0.14917298647260374678782873700D+00
        w(10) = 0.15275338713072585069808433195D+00
        w(11) = 0.15275338713072585069808433195D+00
        w(12) = 0.14917298647260374678782873700D+00
        w(13) = 0.14209610931838205132929832507D+00
        w(14) = 0.13168863844917662689849449975D+00
        w(15) = 0.11819453196151841731237737771D+00
        w(16) = 0.10193011981724043503675013548D+00
        w(17) = 0.08327674157670474872475814322D+00
        w(18) = 0.062672048334109063569506535187D+00
        w(19) = 0.040601429800386941331039952275D+00
        w(20) = 0.017614007139152118311861962352D+00

      else if ( n .eq. 21 ) then

        x( 1) =  -0.99375217062038950026024204D+00
        x( 2) =  -0.96722683856630629431662221D+00
        x( 3) =  -0.92009933415040082879018713D+00
        x( 4) =  -0.85336336458331728364725064D+00
        x( 5) =  -0.76843996347567790861587785D+00
        x( 6) =  -0.66713880419741231930596667D+00
        x( 7) =  -0.55161883588721980705901880D+00
        x( 8) =  -0.42434212020743878357366889D+00
        x( 9) =  -0.28802131680240109660079252D+00
        x(10) =  -0.14556185416089509093703098D+00
        x(11) =   0.00000000000000000000000000D+00
        x(12) =  +0.14556185416089509093703098D+00
        x(13) =  +0.28802131680240109660079252D+00
        x(14) =  +0.42434212020743878357366889D+00
        x(15) =  +0.55161883588721980705901880D+00
        x(16) =  +0.66713880419741231930596667D+00
        x(17) =  +0.76843996347567790861587785D+00
        x(18) =  +0.85336336458331728364725064D+00
        x(19) =  +0.92009933415040082879018713D+00
        x(20) =  +0.96722683856630629431662221D+00
        x(21) =  +0.99375217062038950026024204D+00

        w( 1) =   0.016017228257774333324224616858D+00
        w( 2) =   0.036953789770852493799950668299D+00
        w( 3) =   0.057134425426857208283635826472D+00
        w( 4) =   0.076100113628379302017051653300D+00
        w( 5) =   0.093444423456033861553289741114D+00
        w( 6) =   0.108797299167148377663474578070D+00
        w( 7) =   0.12183141605372853419536717713D+00
        w( 8) =   0.13226893863333746178105257450D+00
        w( 9) =   0.13988739479107315472213342387D+00
        w(10) =   0.14452440398997005906382716655D+00
        w(11) =   0.14608113364969042719198514768D+00
        w(12) =   0.14452440398997005906382716655D+00
        w(13) =   0.13988739479107315472213342387D+00
        w(14) =   0.13226893863333746178105257450D+00
        w(15) =   0.12183141605372853419536717713D+00
        w(16) =   0.108797299167148377663474578070D+00
        w(17) =   0.093444423456033861553289741114D+00
        w(18) =   0.076100113628379302017051653300D+00
        w(19) =   0.057134425426857208283635826472D+00
        w(20) =   0.036953789770852493799950668299D+00
        w(21) =   0.016017228257774333324224616858D+00

      else if ( n .eq. 22 ) then

        x(1) = -0.99429458548239929207303142D+00
        x(2) = -0.97006049783542872712395099D+00
        x(3) = -0.92695677218717400052069294D+00
        x(4) = -0.86581257772030013653642564D+00
        x(5) = -0.78781680597920816200427796D+00
        x(6) = -0.69448726318668278005068984D+00
        x(7) = -0.58764040350691159295887693D+00
        x(8) = -0.46935583798675702640633071D+00
        x(9) = -0.34193582089208422515814742D+00
        x(10) = -0.20786042668822128547884653D+00
        x(11) = -0.06973927331972222121384180D+00
        x(12) = 0.06973927331972222121384180D+00
        x(13) = 0.20786042668822128547884653D+00
        x(14) = 0.34193582089208422515814742D+00
        x(15) = 0.46935583798675702640633071D+00
        x(16) = 0.58764040350691159295887693D+00
        x(17) = 0.69448726318668278005068984D+00
        x(18) = 0.78781680597920816200427796D+00
        x(19) = 0.86581257772030013653642564D+00
        x(20) = 0.92695677218717400052069294D+00
        x(21) = 0.97006049783542872712395099D+00
        x(22) = 0.99429458548239929207303142D+00

        w(1) = 0.014627995298272200684991098047D+00
        w(2) = 0.033774901584814154793302246866D+00
        w(3) = 0.052293335152683285940312051273D+00
        w(4) = 0.06979646842452048809496141893D+00
        w(5) = 0.08594160621706772741444368137D+00
        w(6) = 0.10041414444288096493207883783D+00
        w(7) = 0.11293229608053921839340060742D+00
        w(8) = 0.12325237681051242428556098615D+00
        w(9) = 0.13117350478706237073296499253D+00
        w(10) = 0.13654149834601517135257383123D+00
        w(11) = 0.13925187285563199337541024834D+00
        w(12) = 0.13925187285563199337541024834D+00
        w(13) = 0.13654149834601517135257383123D+00
        w(14) = 0.13117350478706237073296499253D+00
        w(15) = 0.12325237681051242428556098615D+00
        w(16) = 0.11293229608053921839340060742D+00
        w(17) = 0.10041414444288096493207883783D+00
        w(18) = 0.08594160621706772741444368137D+00
        w(19) = 0.06979646842452048809496141893D+00
        w(20) = 0.052293335152683285940312051273D+00
        w(21) = 0.033774901584814154793302246866D+00
        w(22) = 0.014627995298272200684991098047D+00

      else if ( n .eq. 23 ) then

        x(1) = -0.99476933499755212352392572D+00
        x(2) = -0.97254247121811523195602408D+00
        x(3) = -0.93297108682601610234919699D+00
        x(4) = -0.87675235827044166737815689D+00
        x(5) = -0.80488840161883989215111841D+00
        x(6) = -0.71866136313195019446162448D+00
        x(7) = -0.61960987576364615638509731D+00
        x(8) = -0.50950147784600754968979305D+00
        x(9) = -0.39030103803029083142148887D+00
        x(10) = -0.26413568097034493053386954D+00
        x(11) = -0.13325682429846611093174268D+00
        x(12) = 0.00000000000000000000000000D+00
        x(13) = 0.13325682429846611093174268D+00
        x(14) = 0.26413568097034493053386954D+00
        x(15) = 0.39030103803029083142148887D+00
        x(16) = 0.50950147784600754968979305D+00
        x(17) = 0.61960987576364615638509731D+00
        x(18) = 0.71866136313195019446162448D+00
        x(19) = 0.80488840161883989215111841D+00
        x(20) = 0.87675235827044166737815689D+00
        x(21) = 0.93297108682601610234919699D+00
        x(22) = 0.97254247121811523195602408D+00
        x(23) = 0.99476933499755212352392572D+00

        w(1) = 0.013411859487141772081309493459D+00
        w(2) = 0.030988005856979444310694219642D+00
        w(3) = 0.048037671731084668571641071632D+00
        w(4) = 0.064232421408525852127169615159D+00
        w(5) = 0.079281411776718954922892524742D+00
        w(6) = 0.092915766060035147477018617370D+00
        w(7) = 0.104892091464541410074086185015D+00
        w(8) = 0.11499664022241136494164351293D+00
        w(9) = 0.12304908430672953046757840067D+00
        w(10) = 0.12890572218808214997859533940D+00
        w(11) = 0.13246203940469661737164246470D+00
        w(12) = 0.13365457218610617535145711055D+00
        w(13) = 0.13246203940469661737164246470D+00
        w(14) = 0.12890572218808214997859533940D+00
        w(15) = 0.12304908430672953046757840067D+00
        w(16) = 0.11499664022241136494164351293D+00
        w(17) = 0.104892091464541410074086185015D+00
        w(18) = 0.092915766060035147477018617370D+00
        w(19) = 0.079281411776718954922892524742D+00
        w(20) = 0.064232421408525852127169615159D+00
        w(21) = 0.048037671731084668571641071632D+00
        w(22) = 0.030988005856979444310694219642D+00
        w(23) = 0.013411859487141772081309493459D+00

      else if ( n .eq. 24 ) then

        x(1) = -0.99518721999702136017999741D+00
        x(2) = -0.97472855597130949819839199D+00
        x(3) = -0.93827455200273275852364900D+00
        x(4) = -0.88641552700440103421315434D+00
        x(5) = -0.82000198597390292195394987D+00
        x(6) = -0.74012419157855436424382810D+00
        x(7) = -0.64809365193697556925249579D+00
        x(8) = -0.54542147138883953565837562D+00
        x(9) = -0.43379350762604513848708423D+00
        x(10) = -0.31504267969616337438679329D+00
        x(11) = -0.19111886747361630915863982D+00
        x(12) = -0.06405689286260562608504308D+00
        x(13) = 0.06405689286260562608504308D+00
        x(14) = 0.19111886747361630915863982D+00
        x(15) = 0.31504267969616337438679329D+00
        x(16) = 0.43379350762604513848708423D+00
        x(17) = 0.54542147138883953565837562D+00
        x(18) = 0.64809365193697556925249579D+00
        x(19) = 0.74012419157855436424382810D+00
        x(20) = 0.82000198597390292195394987D+00
        x(21) = 0.88641552700440103421315434D+00
        x(22) = 0.93827455200273275852364900D+00
        x(23) = 0.97472855597130949819839199D+00
        x(24) = 0.99518721999702136017999741D+00

        w(1) = 0.012341229799987199546805667070D+00
        w(2) = 0.028531388628933663181307815952D+00
        w(3) = 0.044277438817419806168602748211D+00
        w(4) = 0.059298584915436780746367758500D+00
        w(5) = 0.07334648141108030573403361525D+00
        w(6) = 0.08619016153195327591718520298D+00
        w(7) = 0.09761865210411388826988066446D+00
        w(8) = 0.10744427011596563478257734245D+00
        w(9) = 0.11550566805372560135334448391D+00
        w(10) = 0.12167047292780339120446315348D+00
        w(11) = 0.12583745634682829612137538251D+00
        w(12) = 0.12793819534675215697405616522D+00
        w(13) = 0.12793819534675215697405616522D+00
        w(14) = 0.12583745634682829612137538251D+00
        w(15) = 0.12167047292780339120446315348D+00
        w(16) = 0.11550566805372560135334448391D+00
        w(17) = 0.10744427011596563478257734245D+00
        w(18) = 0.09761865210411388826988066446D+00
        w(19) = 0.08619016153195327591718520298D+00
        w(20) = 0.07334648141108030573403361525D+00
        w(21) = 0.059298584915436780746367758500D+00
        w(22) = 0.044277438817419806168602748211D+00
        w(23) = 0.028531388628933663181307815952D+00
        w(24) = 0.012341229799987199546805667070D+00

      else if ( n .eq. 25 ) then

        x(1) = -0.99555696979049809790878495D+00
        x(2) = -0.97666392145951751149831539D+00
        x(3) = -0.94297457122897433941401117D+00
        x(4) = -0.89499199787827536885104201D+00
        x(5) = -0.83344262876083400142102111D+00
        x(6) = -0.75925926303735763057728287D+00
        x(7) = -0.67356636847346836448512063D+00
        x(8) = -0.57766293024122296772368984D+00
        x(9) = -0.47300273144571496052218212D+00
        x(10) = -0.36117230580938783773582173D+00
        x(11) = -0.24386688372098843204519036D+00
        x(12) = -0.12286469261071039638735982D+00
        x(13) = 0.00000000000000000000000000D+00
        x(14) = 0.12286469261071039638735982D+00
        x(15) = 0.24386688372098843204519036D+00
        x(16) = 0.36117230580938783773582173D+00
        x(17) = 0.47300273144571496052218212D+00
        x(18) = 0.57766293024122296772368984D+00
        x(19) = 0.67356636847346836448512063D+00
        x(20) = 0.75925926303735763057728287D+00
        x(21) = 0.83344262876083400142102111D+00
        x(22) = 0.89499199787827536885104201D+00
        x(23) = 0.94297457122897433941401117D+00
        x(24) = 0.97666392145951751149831539D+00
        x(25) = 0.99555696979049809790878495D+00

        w(1) = 0.0113937985010262879479029641132D+00
        w(2) = 0.026354986615032137261901815295D+00
        w(3) = 0.040939156701306312655623487712D+00
        w(4) = 0.054904695975835191925936891541D+00
        w(5) = 0.068038333812356917207187185657D+00
        w(6) = 0.080140700335001018013234959669D+00
        w(7) = 0.091028261982963649811497220703D+00
        w(8) = 0.100535949067050644202206890393D+00
        w(9) = 0.108519624474263653116093957050D+00
        w(10) = 0.11485825914571164833932554587D+00
        w(11) = 0.11945576353578477222817812651D+00
        w(12) = 0.12224244299031004168895951895D+00
        w(13) = 0.12317605372671545120390287308D+00
        w(14) = 0.12224244299031004168895951895D+00
        w(15) = 0.11945576353578477222817812651D+00
        w(16) = 0.11485825914571164833932554587D+00
        w(17) = 0.108519624474263653116093957050D+00
        w(18) = 0.100535949067050644202206890393D+00
        w(19) = 0.091028261982963649811497220703D+00
        w(20) = 0.080140700335001018013234959669D+00
        w(21) = 0.068038333812356917207187185657D+00
        w(22) = 0.054904695975835191925936891541D+00
        w(23) = 0.040939156701306312655623487712D+00
        w(24) = 0.026354986615032137261901815295D+00
        w(25) = 0.0113937985010262879479029641132D+00

      else if ( n .eq. 26 ) then

        x(1) = -0.99588570114561692900321696D+00
        x(2) = -0.97838544595647099110058035D+00
        x(3) = -0.94715906666171425013591528D+00
        x(4) = -0.90263786198430707421766560D+00
        x(5) = -0.84544594278849801879750706D+00
        x(6) = -0.77638594882067885619296725D+00
        x(7) = -0.69642726041995726486381391D+00
        x(8) = -0.60669229301761806323197875D+00
        x(9) = -0.50844071482450571769570306D+00
        x(10) = -0.40305175512348630648107738D+00
        x(11) = -0.29200483948595689514283538D+00
        x(12) = -0.17685882035689018396905775D+00
        x(13) = -0.05923009342931320709371858D+00
        x(14) = 0.05923009342931320709371858D+00
        x(15) = 0.17685882035689018396905775D+00
        x(16) = 0.29200483948595689514283538D+00
        x(17) = 0.40305175512348630648107738D+00
        x(18) = 0.50844071482450571769570306D+00
        x(19) = 0.60669229301761806323197875D+00
        x(20) = 0.69642726041995726486381391D+00
        x(21) = 0.77638594882067885619296725D+00
        x(22) = 0.84544594278849801879750706D+00
        x(23) = 0.90263786198430707421766560D+00
        x(24) = 0.94715906666171425013591528D+00
        x(25) = 0.97838544595647099110058035D+00
        x(26) = 0.99588570114561692900321696D+00

        w(1) = 0.010551372617343007155651187685D+00
        w(2) = 0.024417851092631908789615827520D+00
        w(3) = 0.037962383294362763950303141249D+00
        w(4) = 0.050975825297147811998319900724D+00
        w(5) = 0.063274046329574835539453689907D+00
        w(6) = 0.07468414976565974588707579610D+00
        w(7) = 0.08504589431348523921044776508D+00
        w(8) = 0.09421380035591414846366488307D+00
        w(9) = 0.10205916109442542323841407025D+00
        w(10) = 0.10847184052857659065657942673D+00
        w(11) = 0.11336181654631966654944071844D+00
        w(12) = 0.11666044348529658204466250754D+00
        w(13) = 0.11832141527926227651637108570D+00
        w(14) = 0.11832141527926227651637108570D+00
        w(15) = 0.11666044348529658204466250754D+00
        w(16) = 0.11336181654631966654944071844D+00
        w(17) = 0.10847184052857659065657942673D+00
        w(18) = 0.10205916109442542323841407025D+00
        w(19) = 0.09421380035591414846366488307D+00
        w(20) = 0.08504589431348523921044776508D+00
        w(21) = 0.07468414976565974588707579610D+00
        w(22) = 0.063274046329574835539453689907D+00
        w(23) = 0.050975825297147811998319900724D+00
        w(24) = 0.037962383294362763950303141249D+00
        w(25) = 0.024417851092631908789615827520D+00
        w(26) = 0.010551372617343007155651187685D+00

      else if ( n .eq. 27 ) then

        x(1) = -0.99617926288898856693888721D+00
        x(2) = -0.97992347596150122285587336D+00
        x(3) = -0.95090055781470500685190803D+00
        x(4) = -0.90948232067749110430064502D+00
        x(5) = -0.85620790801829449030273722D+00
        x(6) = -0.79177163907050822714439734D+00
        x(7) = -0.71701347373942369929481621D+00
        x(8) = -0.63290797194649514092773464D+00
        x(9) = -0.54055156457945689490030094D+00
        x(10) = -0.44114825175002688058597416D+00
        x(11) = -0.33599390363850889973031903D+00
        x(12) = -0.22645936543953685885723911D+00
        x(13) = -0.11397258560952996693289498D+00
        x(14) = 0.00000000000000000000000000D+00
        x(15) = 0.11397258560952996693289498D+00
        x(16) = 0.22645936543953685885723911D+00
        x(17) = 0.33599390363850889973031903D+00
        x(18) = 0.44114825175002688058597416D+00
        x(19) = 0.54055156457945689490030094D+00
        x(20) = 0.63290797194649514092773464D+00
        x(21) = 0.71701347373942369929481621D+00
        x(22) = 0.79177163907050822714439734D+00
        x(23) = 0.85620790801829449030273722D+00
        x(24) = 0.90948232067749110430064502D+00
        x(25) = 0.95090055781470500685190803D+00
        x(26) = 0.97992347596150122285587336D+00
        x(27) = 0.99617926288898856693888721D+00

        w(1) = 0.0097989960512943602611500550912D+00
        w(2) = 0.022686231596180623196034206447D+00
        w(3) = 0.035297053757419711022578289305D+00
        w(4) = 0.047449412520615062704096710114D+00
        w(5) = 0.058983536859833599110300833720D+00
        w(6) = 0.069748823766245592984322888357D+00
        w(7) = 0.079604867773057771263074959010D+00
        w(8) = 0.088423158543756950194322802854D+00
        w(9) = 0.096088727370028507565652646558D+00
        w(10) = 0.102501637817745798671247711533D+00
        w(11) = 0.107578285788533187212162984427D+00
        w(12) = 0.111252488356845192672163096043D+00
        w(13) = 0.113476346108965148620369948092D+00
        w(14) = 0.11422086737895698904504573690D+00
        w(15) = 0.113476346108965148620369948092D+00
        w(16) = 0.111252488356845192672163096043D+00
        w(17) = 0.107578285788533187212162984427D+00
        w(18) = 0.102501637817745798671247711533D+00
        w(19) = 0.096088727370028507565652646558D+00
        w(20) = 0.088423158543756950194322802854D+00
        w(21) = 0.079604867773057771263074959010D+00
        w(22) = 0.069748823766245592984322888357D+00
        w(23) = 0.058983536859833599110300833720D+00
        w(24) = 0.047449412520615062704096710114D+00
        w(25) = 0.035297053757419711022578289305D+00
        w(26) = 0.022686231596180623196034206447D+00
        w(27) = 0.0097989960512943602611500550912D+00

      else if ( n .eq. 28 ) then

        x(1) = -0.99644249757395444995043639D+00
        x(2) = -0.98130316537087275369455995D+00
        x(3) = -0.95425928062893819725410184D+00
        x(4) = -0.91563302639213207386968942D+00
        x(5) = -0.86589252257439504894225457D+00
        x(6) = -0.80564137091717917144788596D+00
        x(7) = -0.73561087801363177202814451D+00
        x(8) = -0.65665109403886496121989818D+00
        x(9) = -0.56972047181140171930800328D+00
        x(10) = -0.47587422495511826103441185D+00
        x(11) = -0.37625151608907871022135721D+00
        x(12) = -0.27206162763517807767682636D+00
        x(13) = -0.16456928213338077128147178D+00
        x(14) = -0.05507928988403427042651653D+00
        x(15) = 0.05507928988403427042651653D+00
        x(16) = 0.16456928213338077128147178D+00
        x(17) = 0.27206162763517807767682636D+00
        x(18) = 0.37625151608907871022135721D+00
        x(19) = 0.47587422495511826103441185D+00
        x(20) = 0.56972047181140171930800328D+00
        x(21) = 0.65665109403886496121989818D+00
        x(22) = 0.73561087801363177202814451D+00
        x(23) = 0.80564137091717917144788596D+00
        x(24) = 0.86589252257439504894225457D+00
        x(25) = 0.91563302639213207386968942D+00
        x(26) = 0.95425928062893819725410184D+00
        x(27) = 0.98130316537087275369455995D+00
        x(28) = 0.99644249757395444995043639D+00

        w(1) = 0.009124282593094517738816153923D+00
        w(2) = 0.021132112592771259751500380993D+00
        w(3) = 0.032901427782304379977630819171D+00
        w(4) = 0.044272934759004227839587877653D+00
        w(5) = 0.055107345675716745431482918227D+00
        w(6) = 0.06527292396699959579339756678D+00
        w(7) = 0.07464621423456877902393188717D+00
        w(8) = 0.08311341722890121839039649824D+00
        w(9) = 0.09057174439303284094218603134D+00
        w(10) = 0.09693065799792991585048900610D+00
        w(11) = 0.10211296757806076981421663851D+00
        w(12) = 0.10605576592284641791041643700D+00
        w(13) = 0.10871119225829413525357151930D+00
        w(14) = 0.11004701301647519628237626560D+00
        w(15) = 0.11004701301647519628237626560D+00
        w(16) = 0.10871119225829413525357151930D+00
        w(17) = 0.10605576592284641791041643700D+00
        w(18) = 0.10211296757806076981421663851D+00
        w(19) = 0.09693065799792991585048900610D+00
        w(20) = 0.09057174439303284094218603134D+00
        w(21) = 0.08311341722890121839039649824D+00
        w(22) = 0.07464621423456877902393188717D+00
        w(23) = 0.06527292396699959579339756678D+00
        w(24) = 0.055107345675716745431482918227D+00
        w(25) = 0.044272934759004227839587877653D+00
        w(26) = 0.032901427782304379977630819171D+00
        w(27) = 0.021132112592771259751500380993D+00
        w(28) = 0.009124282593094517738816153923D+00

      else if ( n .eq. 29 ) then

        x(1) = -0.99667944226059658616319153D+00
        x(2) = -0.98254550526141317487092602D+00
        x(3) = -0.95728559577808772579820804D+00
        x(4) = -0.92118023295305878509375344D+00
        x(5) = -0.87463780492010279041779342D+00
        x(6) = -0.81818548761525244498957221D+00
        x(7) = -0.75246285173447713391261008D+00
        x(8) = -0.67821453760268651515618501D+00
        x(9) = -0.59628179713822782037958621D+00
        x(10) = -0.50759295512422764210262792D+00
        x(11) = -0.41315288817400866389070659D+00
        x(12) = -0.31403163786763993494819592D+00
        x(13) = -0.21135228616600107450637573D+00
        x(14) = -0.10627823013267923017098239D+00
        x(15) = 0.00000000000000000000000000D+00
        x(16) = 0.10627823013267923017098239D+00
        x(17) = 0.21135228616600107450637573D+00
        x(18) = 0.31403163786763993494819592D+00
        x(19) = 0.41315288817400866389070659D+00
        x(20) = 0.50759295512422764210262792D+00
        x(21) = 0.59628179713822782037958621D+00
        x(22) = 0.67821453760268651515618501D+00
        x(23) = 0.75246285173447713391261008D+00
        x(24) = 0.81818548761525244498957221D+00
        x(25) = 0.87463780492010279041779342D+00
        x(26) = 0.92118023295305878509375344D+00
        x(27) = 0.95728559577808772579820804D+00
        x(28) = 0.98254550526141317487092602D+00
        x(29) = 0.99667944226059658616319153D+00

        w(1) = 0.0085169038787464096542638133022D+00
        w(2) = 0.019732085056122705983859801640D+00
        w(3) = 0.030740492202093622644408525375D+00
        w(4) = 0.041402062518682836104830010114D+00
        w(5) = 0.051594826902497923912594381180D+00
        w(6) = 0.061203090657079138542109848024D+00
        w(7) = 0.070117933255051278569581486949D+00
        w(8) = 0.078238327135763783828144888660D+00
        w(9) = 0.085472257366172527545344849297D+00
        w(10) = 0.091737757139258763347966411077D+00
        w(11) = 0.096963834094408606301900074883D+00
        w(12) = 0.101091273759914966121820546907D+00
        w(13) = 0.104073310077729373913328471285D+00
        w(14) = 0.105876155097320941406591327852D+00
        w(15) = 0.10647938171831424424651112691D+00
        w(16) = 0.105876155097320941406591327852D+00
        w(17) = 0.104073310077729373913328471285D+00
        w(18) = 0.101091273759914966121820546907D+00
        w(19) = 0.096963834094408606301900074883D+00
        w(20) = 0.091737757139258763347966411077D+00
        w(21) = 0.085472257366172527545344849297D+00
        w(22) = 0.078238327135763783828144888660D+00
        w(23) = 0.070117933255051278569581486949D+00
        w(24) = 0.061203090657079138542109848024D+00
        w(25) = 0.051594826902497923912594381180D+00
        w(26) = 0.041402062518682836104830010114D+00
        w(27) = 0.030740492202093622644408525375D+00
        w(28) = 0.019732085056122705983859801640D+00
        w(29) = 0.0085169038787464096542638133022D+00

      else if ( n .eq. 30 ) then

        x(1) = -0.99689348407464954027163005D+00
        x(2) = -0.98366812327974720997003258D+00
        x(3) = -0.96002186496830751221687103D+00
        x(4) = -0.92620004742927432587932428D+00
        x(5) = -0.88256053579205268154311646D+00
        x(6) = -0.82956576238276839744289812D+00
        x(7) = -0.76777743210482619491797734D+00
        x(8) = -0.69785049479331579693229239D+00
        x(9) = -0.62052618298924286114047756D+00
        x(10) = -0.53662414814201989926416979D+00
        x(11) = -0.44703376953808917678060990D+00
        x(12) = -0.35270472553087811347103721D+00
        x(13) = -0.25463692616788984643980513D+00
        x(14) = -0.15386991360858354696379467D+00
        x(15) = -0.05147184255531769583302521D+00
        x(16) = 0.05147184255531769583302521D+00
        x(17) = 0.15386991360858354696379467D+00
        x(18) = 0.25463692616788984643980513D+00
        x(19) = 0.35270472553087811347103721D+00
        x(20) = 0.44703376953808917678060990D+00
        x(21) = 0.53662414814201989926416979D+00
        x(22) = 0.62052618298924286114047756D+00
        x(23) = 0.69785049479331579693229239D+00
        x(24) = 0.76777743210482619491797734D+00
        x(25) = 0.82956576238276839744289812D+00
        x(26) = 0.88256053579205268154311646D+00
        x(27) = 0.92620004742927432587932428D+00
        x(28) = 0.96002186496830751221687103D+00
        x(29) = 0.98366812327974720997003258D+00
        x(30) = 0.99689348407464954027163005D+00

        w(1) = 0.007968192496166605615465883475D+00
        w(2) = 0.018466468311090959142302131912D+00
        w(3) = 0.028784707883323369349719179611D+00
        w(4) = 0.038799192569627049596801936446D+00
        w(5) = 0.048402672830594052902938140423D+00
        w(6) = 0.057493156217619066481721689402D+00
        w(7) = 0.06597422988218049512812851512D+00
        w(8) = 0.07375597473770520626824385002D+00
        w(9) = 0.08075589522942021535469493846D+00
        w(10) = 0.08689978720108297980238753072D+00
        w(11) = 0.09212252223778612871763270709D+00
        w(12) = 0.09636873717464425963946862635D+00
        w(13) = 0.09959342058679526706278028210D+00
        w(14) = 0.10176238974840550459642895217D+00
        w(15) = 0.10285265289355884034128563671D+00
        w(16) = 0.10285265289355884034128563671D+00
        w(17) = 0.10176238974840550459642895217D+00
        w(18) = 0.09959342058679526706278028210D+00
        w(19) = 0.09636873717464425963946862635D+00
        w(20) = 0.09212252223778612871763270709D+00
        w(21) = 0.08689978720108297980238753072D+00
        w(22) = 0.08075589522942021535469493846D+00
        w(23) = 0.07375597473770520626824385002D+00
        w(24) = 0.06597422988218049512812851512D+00
        w(25) = 0.057493156217619066481721689402D+00
        w(26) = 0.048402672830594052902938140423D+00
        w(27) = 0.038799192569627049596801936446D+00
        w(28) = 0.028784707883323369349719179611D+00
        w(29) = 0.018466468311090959142302131912D+00
        w(30) = 0.007968192496166605615465883475D+00

      else if ( n .eq. 31 ) then

        x(1) = -0.99708748181947707405562655D+00
        x(2) = -0.98468590966515248400246517D+00
        x(3) = -0.96250392509294966178905240D+00
        x(4) = -0.93075699789664816495694576D+00
        x(5) = -0.88976002994827104337419201D+00
        x(6) = -0.83992032014626734008690454D+00
        x(7) = -0.78173314841662494040636002D+00
        x(8) = -0.71577678458685328390597087D+00
        x(9) = -0.64270672292426034618441820D+00
        x(10) = -0.56324916140714926272094492D+00
        x(11) = -0.47819378204490248044059404D+00
        x(12) = -0.38838590160823294306135146D+00
        x(13) = -0.29471806998170161661790390D+00
        x(14) = -0.19812119933557062877241300D+00
        x(15) = -0.09955531215234152032517479D+00
        x(16) = 0.00000000000000000000000000D+00
        x(17) = 0.09955531215234152032517479D+00
        x(18) = 0.19812119933557062877241300D+00
        x(19) = 0.29471806998170161661790390D+00
        x(20) = 0.38838590160823294306135146D+00
        x(21) = 0.47819378204490248044059404D+00
        x(22) = 0.56324916140714926272094492D+00
        x(23) = 0.64270672292426034618441820D+00
        x(24) = 0.71577678458685328390597087D+00
        x(25) = 0.78173314841662494040636002D+00
        x(26) = 0.83992032014626734008690454D+00
        x(27) = 0.88976002994827104337419201D+00
        x(28) = 0.93075699789664816495694576D+00
        x(29) = 0.96250392509294966178905240D+00
        x(30) = 0.98468590966515248400246517D+00
        x(31) = 0.99708748181947707405562655D+00

        w(1) = 0.0074708315792487758586968750322D+00
        w(2) = 0.017318620790310582463157996087D+00
        w(3) = 0.027009019184979421800608708092D+00
        w(4) = 0.036432273912385464024392010468D+00
        w(5) = 0.045493707527201102902315857895D+00
        w(6) = 0.054103082424916853711666259087D+00
        w(7) = 0.062174786561028426910343543687D+00
        w(8) = 0.069628583235410366167756126255D+00
        w(9) = 0.076390386598776616426357674901D+00
        w(10) = 0.082392991761589263903823367432D+00
        w(11) = 0.087576740608477876126198069695D+00
        w(12) = 0.091890113893641478215362871607D+00
        w(13) = 0.095290242912319512807204197488D+00
        w(14) = 0.097743335386328725093474010979D+00
        w(15) = 0.099225011226672307874875514429D+00
        w(16) = 0.09972054479342645142753383373D+00
        w(17) = 0.099225011226672307874875514429D+00
        w(18) = 0.097743335386328725093474010979D+00
        w(19) = 0.095290242912319512807204197488D+00
        w(20) = 0.091890113893641478215362871607D+00
        w(21) = 0.087576740608477876126198069695D+00
        w(22) = 0.082392991761589263903823367432D+00
        w(23) = 0.076390386598776616426357674901D+00
        w(24) = 0.069628583235410366167756126255D+00
        w(25) = 0.062174786561028426910343543687D+00
        w(26) = 0.054103082424916853711666259087D+00
        w(27) = 0.045493707527201102902315857895D+00
        w(28) = 0.036432273912385464024392010468D+00
        w(29) = 0.027009019184979421800608708092D+00
        w(30) = 0.017318620790310582463157996087D+00
        w(31) = 0.0074708315792487758586968750322D+00

      else if ( n .eq. 32 ) then

        x(1) = -0.99726386184948156354498113D+00
        x(2) = -0.98561151154526833540017504D+00
        x(3) = -0.96476225558750643077381193D+00
        x(4) = -0.93490607593773968917091913D+00
        x(5) = -0.89632115576605212396530724D+00
        x(6) = -0.84936761373256997013369300D+00
        x(7) = -0.79448379596794240696309730D+00
        x(8) = -0.73218211874028968038742667D+00
        x(9) = -0.66304426693021520097511517D+00
        x(10) = -0.58771575724076232904074548D+00
        x(11) = -0.50689990893222939002374747D+00
        x(12) = -0.42135127613063534536411944D+00
        x(13) = -0.33186860228212764977991681D+00
        x(14) = -0.23928736225213707454460321D+00
        x(15) = -0.14447196158279649348518637D+00
        x(16) = -0.04830766568773831623481257D+00
        x(17) = 0.04830766568773831623481257D+00
        x(18) = 0.14447196158279649348518637D+00
        x(19) = 0.23928736225213707454460321D+00
        x(20) = 0.33186860228212764977991681D+00
        x(21) = 0.42135127613063534536411944D+00
        x(22) = 0.50689990893222939002374747D+00
        x(23) = 0.58771575724076232904074548D+00
        x(24) = 0.66304426693021520097511517D+00
        x(25) = 0.73218211874028968038742667D+00
        x(26) = 0.79448379596794240696309730D+00
        x(27) = 0.84936761373256997013369300D+00
        x(28) = 0.89632115576605212396530724D+00
        x(29) = 0.93490607593773968917091913D+00
        x(30) = 0.96476225558750643077381193D+00
        x(31) = 0.98561151154526833540017504D+00
        x(32) = 0.99726386184948156354498113D+00

        w(1) = 0.007018610009470096600407063739D+00
        w(2) = 0.016274394730905670605170562206D+00
        w(3) = 0.025392065309262059455752589789D+00
        w(4) = 0.034273862913021433102687732252D+00
        w(5) = 0.042835898022226680656878646606D+00
        w(6) = 0.050998059262376176196163244690D+00
        w(7) = 0.058684093478535547145283637300D+00
        w(8) = 0.06582222277636184683765006371D+00
        w(9) = 0.07234579410884850622539935648D+00
        w(10) = 0.07819389578707030647174091883D+00
        w(11) = 0.08331192422694675522219907460D+00
        w(12) = 0.08765209300440381114277146275D+00
        w(13) = 0.09117387869576388471286857711D+00
        w(14) = 0.09384439908080456563918023767D+00
        w(15) = 0.09563872007927485941908200220D+00
        w(16) = 0.09654008851472780056676483006D+00
        w(17) = 0.09654008851472780056676483006D+00
        w(18) = 0.09563872007927485941908200220D+00
        w(19) = 0.09384439908080456563918023767D+00
        w(20) = 0.09117387869576388471286857711D+00
        w(21) = 0.08765209300440381114277146275D+00
        w(22) = 0.08331192422694675522219907460D+00
        w(23) = 0.07819389578707030647174091883D+00
        w(24) = 0.07234579410884850622539935648D+00
        w(25) = 0.06582222277636184683765006371D+00
        w(26) = 0.058684093478535547145283637300D+00
        w(27) = 0.050998059262376176196163244690D+00
        w(28) = 0.042835898022226680656878646606D+00
        w(29) = 0.034273862913021433102687732252D+00
        w(30) = 0.025392065309262059455752589789D+00
        w(31) = 0.016274394730905670605170562206D+00
        w(32) = 0.007018610009470096600407063739D+00

      else if ( n .eq. 33 ) then

        x(1) = -0.99742469424645521726616802D+00
        x(2) = -0.98645572623064248811037570D+00
        x(3) = -0.96682290968999276892837771D+00
        x(4) = -0.93869437261116835035583512D+00
        x(5) = -0.90231676774343358304053133D+00
        x(6) = -0.85800965267650406464306148D+00
        x(7) = -0.80616235627416658979620087D+00
        x(8) = -0.74723049644956215785905512D+00
        x(9) = -0.68173195996974278626821595D+00
        x(10) = -0.61024234583637902730728751D+00
        x(11) = -0.53338990478634764354889426D+00
        x(12) = -0.45185001727245069572599328D+00
        x(13) = -0.36633925774807334107022062D+00
        x(14) = -0.27760909715249702940324807D+00
        x(15) = -0.18643929882799157233579876D+00
        x(16) = -0.09363106585473338567074292D+00
        x(17) = 0.00000000000000000000000000D+00
        x(18) = 0.09363106585473338567074292D+00
        x(19) = 0.18643929882799157233579876D+00
        x(20) = 0.27760909715249702940324807D+00
        x(21) = 0.36633925774807334107022062D+00
        x(22) = 0.45185001727245069572599328D+00
        x(23) = 0.53338990478634764354889426D+00
        x(24) = 0.61024234583637902730728751D+00
        x(25) = 0.68173195996974278626821595D+00
        x(26) = 0.74723049644956215785905512D+00
        x(27) = 0.80616235627416658979620087D+00
        x(28) = 0.85800965267650406464306148D+00
        x(29) = 0.90231676774343358304053133D+00
        x(30) = 0.93869437261116835035583512D+00
        x(31) = 0.96682290968999276892837771D+00
        x(32) = 0.98645572623064248811037570D+00
        x(33) = 0.99742469424645521726616802D+00

        w(1) = 0.0066062278475873780586492352085D+00
        w(2) = 0.015321701512934676127945768534D+00
        w(3) = 0.023915548101749480350533257529D+00
        w(4) = 0.032300358632328953281561447250D+00
        w(5) = 0.040401541331669591563409790527D+00
        w(6) = 0.048147742818711695670146880138D+00
        w(7) = 0.055470846631663561284944495439D+00
        w(8) = 0.062306482530317480031627725771D+00
        w(9) = 0.068594572818656712805955073015D+00
        w(10) = 0.074279854843954149342472175919D+00
        w(11) = 0.079312364794886738363908384942D+00
        w(12) = 0.083647876067038707613928014518D+00
        w(13) = 0.087248287618844337607281670945D+00
        w(14) = 0.090081958660638577239743705500D+00
        w(15) = 0.092123986643316846213240977717D+00
        w(16) = 0.093356426065596116160999126274D+00
        w(17) = 0.09376844616020999656730454155D+00
        w(18) = 0.093356426065596116160999126274D+00
        w(19) = 0.092123986643316846213240977717D+00
        w(20) = 0.090081958660638577239743705500D+00
        w(21) = 0.087248287618844337607281670945D+00
        w(22) = 0.083647876067038707613928014518D+00
        w(23) = 0.079312364794886738363908384942D+00
        w(24) = 0.074279854843954149342472175919D+00
        w(25) = 0.068594572818656712805955073015D+00
        w(26) = 0.062306482530317480031627725771D+00
        w(27) = 0.055470846631663561284944495439D+00
        w(28) = 0.048147742818711695670146880138D+00
        w(29) = 0.040401541331669591563409790527D+00
        w(30) = 0.032300358632328953281561447250D+00
        w(31) = 0.023915548101749480350533257529D+00
        w(32) = 0.015321701512934676127945768534D+00
        w(33) = 0.0066062278475873780586492352085D+00

      else if ( n .eq. 63 ) then

        x(1) = -0.99928298402912378037893614D+00
        x(2) = -0.99622401277797010860219336D+00
        x(3) = -0.99072854689218946681089467D+00
        x(4) = -0.98280881059372723486251141D+00
        x(5) = -0.97248403469757002280196068D+00
        x(6) = -0.95977944975894192707035417D+00
        x(7) = -0.94472613404100980296637532D+00
        x(8) = -0.92736092062184320544703138D+00
        x(9) = -0.90772630277853155803695313D+00
        x(10) = -0.88587032850785342629029846D+00
        x(11) = -0.86184648236412371953961184D+00
        x(12) = -0.83571355431950284347180777D+00
        x(13) = -0.80753549577345676005146599D+00
        x(14) = -0.7773812629903723355633302D+00
        x(15) = -0.7453246483178474178293217D+00
        x(16) = -0.7114440995848458078514315D+00
        x(17) = -0.6758225281149860901311033D+00
        x(18) = -0.6385471058213653850003070D+00
        x(19) = -0.5997090518776252357390089D+00
        x(20) = -0.5594034094862850132676978D+00
        x(21) = -0.5177288132900332481244776D+00
        x(22) = -0.4747872479948043999222123D+00
        x(23) = -0.4306837987951116006620889D+00
        x(24) = -0.3855263942122478924776150D+00
        x(25) = -0.3394255419745844024688344D+00
        x(26) = -0.2924940585862514400361572D+00
        x(27) = -0.2448467932459533627484046D+00
        x(28) = -0.1966003467915066845576275D+00
        x(29) = -0.1478727863578719685698391D+00
        x(30) = -0.0987833564469452795297037D+00
        x(31) = -0.0494521871161596272342338D+00
        x(32) = 0.0000000000000000000000000D+00
        x(33) = 0.0494521871161596272342338D+00
        x(34) = 0.0987833564469452795297037D+00
        x(35) = 0.1478727863578719685698391D+00
        x(36) = 0.1966003467915066845576275D+00
        x(37) = 0.2448467932459533627484046D+00
        x(38) = 0.2924940585862514400361572D+00
        x(39) = 0.3394255419745844024688344D+00
        x(40) = 0.3855263942122478924776150D+00
        x(41) = 0.4306837987951116006620889D+00
        x(42) = 0.4747872479948043999222123D+00
        x(43) = 0.5177288132900332481244776D+00
        x(44) = 0.5594034094862850132676978D+00
        x(45) = 0.5997090518776252357390089D+00
        x(46) = 0.6385471058213653850003070D+00
        x(47) = 0.6758225281149860901311033D+00
        x(48) = 0.7114440995848458078514315D+00
        x(49) = 0.7453246483178474178293217D+00
        x(50) = 0.7773812629903723355633302D+00
        x(51) = 0.8075354957734567600514660D+00
        x(52) = 0.8357135543195028434718078D+00
        x(53) = 0.8618464823641237195396118D+00
        x(54) = 0.8858703285078534262902985D+00
        x(55) = 0.9077263027785315580369531D+00
        x(56) = 0.9273609206218432054470314D+00
        x(57) = 0.9447261340410098029663753D+00
        x(58) = 0.9597794497589419270703542D+00
        x(59) = 0.9724840346975700228019607D+00
        x(60) = 0.9828088105937272348625114D+00
        x(61) = 0.9907285468921894668108947D+00
        x(62) = 0.9962240127779701086021934D+00
        x(63) = 0.9992829840291237803789361D+00

        w(1) = 0.0018398745955770841170924455540D+00
        w(2) = 0.0042785083468637618660784110826D+00
        w(3) = 0.0067102917659601362519069307298D+00
        w(4) = 0.0091259686763266563540586454218D+00
        w(5) = 0.011519376076880041750750606149D+00
        w(6) = 0.013884612616115610824866086368D+00
        w(7) = 0.016215878410338338882283672975D+00
        w(8) = 0.018507464460161270409260545805D+00
        w(9) = 0.020753761258039090775341953421D+00
        w(10) = 0.022949271004889933148942319562D+00
        w(11) = 0.025088620553344986618630138068D+00
        w(12) = 0.027166574359097933225189839439D+00
        w(13) = 0.029178047208280526945551502154D+00
        w(14) = 0.031118116622219817508215988557D+00
        w(15) = 0.032982034883779341765683179672D+00
        w(16) = 0.034765240645355877697180504643D+00
        w(17) = 0.036463370085457289630452409788D+00
        w(18) = 0.038072267584349556763638324928D+00
        w(19) = 0.039587995891544093984807928149D+00
        w(20) = 0.041006845759666398635110037009D+00
        w(21) = 0.042325345020815822982505485403D+00
        w(22) = 0.043540267083027590798964315704D+00
        w(23) = 0.044648638825941395370332669517D+00
        w(24) = 0.045647747876292608685885992609D+00
        w(25) = 0.046535149245383696510395418747D+00
        w(26) = 0.047308671312268919080604988339D+00
        w(27) = 0.047966421137995131411052756195D+00
        w(28) = 0.048506789097883847864090099146D+00
        w(29) = 0.048928452820511989944709361549D+00
        w(30) = 0.049230380423747560785043116988D+00
        w(31) = 0.049411833039918178967039646117D+00
        w(32) = 0.04947236662393102088866936042D+00
        w(33) = 0.049411833039918178967039646117D+00
        w(34) = 0.049230380423747560785043116988D+00
        w(35) = 0.048928452820511989944709361549D+00
        w(36) = 0.048506789097883847864090099146D+00
        w(37) = 0.047966421137995131411052756195D+00
        w(38) = 0.047308671312268919080604988339D+00
        w(39) = 0.046535149245383696510395418747D+00
        w(40) = 0.045647747876292608685885992609D+00
        w(41) = 0.044648638825941395370332669517D+00
        w(42) = 0.043540267083027590798964315704D+00
        w(43) = 0.042325345020815822982505485403D+00
        w(44) = 0.041006845759666398635110037009D+00
        w(45) = 0.039587995891544093984807928149D+00
        w(46) = 0.038072267584349556763638324928D+00
        w(47) = 0.036463370085457289630452409788D+00
        w(48) = 0.034765240645355877697180504643D+00
        w(49) = 0.032982034883779341765683179672D+00
        w(50) = 0.031118116622219817508215988557D+00
        w(51) = 0.029178047208280526945551502154D+00
        w(52) = 0.027166574359097933225189839439D+00
        w(53) = 0.025088620553344986618630138068D+00
        w(54) = 0.022949271004889933148942319562D+00
        w(55) = 0.020753761258039090775341953421D+00
        w(56) = 0.018507464460161270409260545805D+00
        w(57) = 0.016215878410338338882283672975D+00
        w(58) = 0.013884612616115610824866086368D+00
        w(59) = 0.011519376076880041750750606149D+00
        w(60) = 0.0091259686763266563540586454218D+00
        w(61) = 0.0067102917659601362519069307298D+00
        w(62) = 0.0042785083468637618660784110826D+00
        w(63) = 0.0018398745955770841170924455540D+00

      else if ( n .eq. 64 ) then

        x(1) = -0.99930504173577213945690562D+00
        x(2) = -0.99634011677195527934692450D+00
        x(3) = -0.99101337147674432073938238D+00
        x(4) = -0.98333625388462595693129930D+00
        x(5) = -0.97332682778991096374185351D+00
        x(6) = -0.96100879965205371891861412D+00
        x(7) = -0.94641137485840281606248149D+00
        x(8) = -0.92956917213193957582149015D+00
        x(9) = -0.91052213707850280575638067D+00
        x(10) = -0.88931544599511410585340404D+00
        x(11) = -0.86599939815409281976078339D+00
        x(12) = -0.8406292962525803627516915D+00
        x(13) = -0.8132653151227975597419233D+00
        x(14) = -0.7839723589433414076102205D+00
        x(15) = -0.7528199072605318966118638D+00
        x(16) = -0.7198818501716108268489402D+00
        x(17) = -0.6852363130542332425635584D+00
        x(18) = -0.6489654712546573398577612D+00
        x(19) = -0.6111553551723932502488530D+00
        x(20) = -0.5718956462026340342838781D+00
        x(21) = -0.5312794640198945456580139D+00
        x(22) = -0.4894031457070529574785263D+00
        x(23) = -0.4463660172534640879849477D+00
        x(24) = -0.4022701579639916036957668D+00
        x(25) = -0.3572201583376681159504426D+00
        x(26) = -0.3113228719902109561575127D+00
        x(27) = -0.2646871622087674163739642D+00
        x(28) = -0.2174236437400070841496487D+00
        x(29) = -0.1696444204239928180373136D+00
        x(30) = -0.1214628192961205544703765D+00
        x(31) = -0.0729931217877990394495429D+00
        x(32) = -0.0243502926634244325089558D+00
        x(33) = 0.0243502926634244325089558D+00
        x(34) = 0.0729931217877990394495429D+00
        x(35) = 0.1214628192961205544703765D+00
        x(36) = 0.1696444204239928180373136D+00
        x(37) = 0.2174236437400070841496487D+00
        x(38) = 0.2646871622087674163739642D+00
        x(39) = 0.3113228719902109561575127D+00
        x(40) = 0.3572201583376681159504426D+00
        x(41) = 0.4022701579639916036957668D+00
        x(42) = 0.4463660172534640879849477D+00
        x(43) = 0.4894031457070529574785263D+00
        x(44) = 0.5312794640198945456580139D+00
        x(45) = 0.5718956462026340342838781D+00
        x(46) = 0.6111553551723932502488530D+00
        x(47) = 0.6489654712546573398577612D+00
        x(48) = 0.6852363130542332425635584D+00
        x(49) = 0.7198818501716108268489402D+00
        x(50) = 0.7528199072605318966118638D+00
        x(51) = 0.7839723589433414076102205D+00
        x(52) = 0.8132653151227975597419233D+00
        x(53) = 0.8406292962525803627516915D+00
        x(54) = 0.8659993981540928197607834D+00
        x(55) = 0.8893154459951141058534040D+00
        x(56) = 0.9105221370785028057563807D+00
        x(57) = 0.9295691721319395758214902D+00
        x(58) = 0.9464113748584028160624815D+00
        x(59) = 0.9610087996520537189186141D+00
        x(60) = 0.9733268277899109637418535D+00
        x(61) = 0.9833362538846259569312993D+00
        x(62) = 0.9910133714767443207393824D+00
        x(63) = 0.9963401167719552793469245D+00
        x(64) = 0.9993050417357721394569056D+00

        w(1) = 0.0017832807216964329472960791450D+00
        w(2) = 0.0041470332605624676352875357286D+00
        w(3) = 0.006504457968978362856117360400D+00
        w(4) = 0.008846759826363947723030914660D+00
        w(5) = 0.011168139460131128818590493019D+00
        w(6) = 0.013463047896718642598060766686D+00
        w(7) = 0.015726030476024719321965995298D+00
        w(8) = 0.017951715775697343085045302001D+00
        w(9) = 0.020134823153530209372340316729D+00
        w(10) = 0.022270173808383254159298330384D+00
        w(11) = 0.024352702568710873338177550409D+00
        w(12) = 0.026377469715054658671691792625D+00
        w(13) = 0.028339672614259483227511305200D+00
        w(14) = 0.030234657072402478867974059820D+00
        w(15) = 0.032057928354851553585467504348D+00
        w(16) = 0.033805161837141609391565482111D+00
        w(17) = 0.035472213256882383810693146715D+00
        w(18) = 0.037055128540240046040415101810D+00
        w(19) = 0.038550153178615629128962496947D+00
        w(20) = 0.039953741132720341386656926128D+00
        w(21) = 0.041262563242623528610156297474D+00
        w(22) = 0.042473515123653589007339767909D+00
        w(23) = 0.043583724529323453376827860974D+00
        w(24) = 0.044590558163756563060134710031D+00
        w(25) = 0.045491627927418144479770996971D+00
        w(26) = 0.046284796581314417295953249232D+00
        w(27) = 0.046968182816210017325326285755D+00
        w(28) = 0.047540165714830308662282206944D+00
        w(29) = 0.04799938859645830772812617987D+00
        w(30) = 0.04834476223480295716976952716D+00
        w(31) = 0.04857546744150342693479906678D+00
        w(32) = 0.04869095700913972038336539073D+00
        w(33) = 0.04869095700913972038336539073D+00
        w(34) = 0.04857546744150342693479906678D+00
        w(35) = 0.04834476223480295716976952716D+00
        w(36) = 0.04799938859645830772812617987D+00
        w(37) = 0.047540165714830308662282206944D+00
        w(38) = 0.046968182816210017325326285755D+00
        w(39) = 0.046284796581314417295953249232D+00
        w(40) = 0.045491627927418144479770996971D+00
        w(41) = 0.044590558163756563060134710031D+00
        w(42) = 0.043583724529323453376827860974D+00
        w(43) = 0.042473515123653589007339767909D+00
        w(44) = 0.041262563242623528610156297474D+00
        w(45) = 0.039953741132720341386656926128D+00
        w(46) = 0.038550153178615629128962496947D+00
        w(47) = 0.037055128540240046040415101810D+00
        w(48) = 0.035472213256882383810693146715D+00
        w(49) = 0.033805161837141609391565482111D+00
        w(50) = 0.032057928354851553585467504348D+00
        w(51) = 0.030234657072402478867974059820D+00
        w(52) = 0.028339672614259483227511305200D+00
        w(53) = 0.026377469715054658671691792625D+00
        w(54) = 0.024352702568710873338177550409D+00
        w(55) = 0.022270173808383254159298330384D+00
        w(56) = 0.020134823153530209372340316729D+00
        w(57) = 0.017951715775697343085045302001D+00
        w(58) = 0.015726030476024719321965995298D+00
        w(59) = 0.013463047896718642598060766686D+00
        w(60) = 0.011168139460131128818590493019D+00
        w(61) = 0.008846759826363947723030914660D+00
        w(62) = 0.006504457968978362856117360400D+00
        w(63) = 0.0041470332605624676352875357286D+00
        w(64) = 0.0017832807216964329472960791450D+00

      else if ( n .eq. 65 ) then

        x(1) = -0.99932609707541287726569361D+00
        x(2) = -0.99645094806184916305579494D+00
        x(3) = -0.99128527617680166872182118D+00
        x(4) = -0.98383981218703494137763778D+00
        x(5) = -0.97413153983355116907496789D+00
        x(6) = -0.96218275471805523771198375D+00
        x(7) = -0.94802092816840750637376974D+00
        x(8) = -0.93167862822874933796567699D+00
        x(9) = -0.91319344054284626173654692D+00
        x(10) = -0.89260788050473893142328554D+00
        x(11) = -0.8699692949264070361941320D+00
        x(12) = -0.8453297528999302839424500D+00
        x(13) = -0.8187459259226514534339191D+00
        x(14) = -0.7902789574921218430473804D+00
        x(15) = -0.7599943224419997868739828D+00
        x(16) = -0.7279616763294246790119737D+00
        x(17) = -0.6942546952139916335526225D+00
        x(18) = -0.6589509061936251330409408D+00
        x(19) = -0.6221315090854002415825996D+00
        x(20) = -0.5838811896604873133271545D+00
        x(21) = -0.5442879248622271385455725D+00
        x(22) = -0.5034427804550068823410431D+00
        x(23) = -0.4614397015691450576978341D+00
        x(24) = -0.4183752966234090092641990D+00
        x(25) = -0.3743486151220660120087939D+00
        x(26) = -0.3294609198374864076452867D+00
        x(27) = -0.2838154539022487306176554D+00
        x(28) = -0.2375172033464168065707124D+00
        x(29) = -0.1906726556261427697749124D+00
        x(30) = -0.1433895546989751711312496D+00
        x(31) = -0.0957766532091975056522186D+00
        x(32) = -0.0479434623531718575225298D+00
        x(33) = 0.0000000000000000000000000D+00
        x(34) = 0.0479434623531718575225298D+00
        x(35) = 0.0957766532091975056522186D+00
        x(36) = 0.1433895546989751711312496D+00
        x(37) = 0.1906726556261427697749124D+00
        x(38) = 0.2375172033464168065707124D+00
        x(39) = 0.2838154539022487306176554D+00
        x(40) = 0.3294609198374864076452867D+00
        x(41) = 0.3743486151220660120087939D+00
        x(42) = 0.4183752966234090092641990D+00
        x(43) = 0.4614397015691450576978341D+00
        x(44) = 0.5034427804550068823410431D+00
        x(45) = 0.5442879248622271385455725D+00
        x(46) = 0.5838811896604873133271545D+00
        x(47) = 0.6221315090854002415825996D+00
        x(48) = 0.6589509061936251330409408D+00
        x(49) = 0.6942546952139916335526225D+00
        x(50) = 0.7279616763294246790119737D+00
        x(51) = 0.7599943224419997868739828D+00
        x(52) = 0.7902789574921218430473804D+00
        x(53) = 0.8187459259226514534339191D+00
        x(54) = 0.8453297528999302839424500D+00
        x(55) = 0.8699692949264070361941320D+00
        x(56) = 0.8926078805047389314232855D+00
        x(57) = 0.9131934405428462617365469D+00
        x(58) = 0.9316786282287493379656770D+00
        x(59) = 0.9480209281684075063737697D+00
        x(60) = 0.9621827547180552377119837D+00
        x(61) = 0.9741315398335511690749679D+00
        x(62) = 0.9838398121870349413776378D+00
        x(63) = 0.9912852761768016687218212D+00
        x(64) = 0.9964509480618491630557949D+00
        x(65) = 0.9993260970754128772656936D+00

        w(1) = 0.0017292582513002508983395851463D+00
        w(2) = 0.0040215241720037363470786599528D+00
        w(3) = 0.0063079425789717545501888719039D+00
        w(4) = 0.0085801482668814598936358121592D+00
        w(5) = 0.0108326787895979686215140551272D+00
        w(6) = 0.013060311639994846336168342922D+00
        w(7) = 0.015257912146448310349265388145D+00
        w(8) = 0.017420421997670248495365759969D+00
        w(9) = 0.019542865836750062826837429313D+00
        w(10) = 0.021620361284934062841654274667D+00
        w(11) = 0.023648129691287236698780978994D+00
        w(12) = 0.025621506938037758214084978694D+00
        w(13) = 0.027535954088450343942499722327D+00
        w(14) = 0.029387067789310668062644859210D+00
        w(15) = 0.031170590380189142464431845777D+00
        w(16) = 0.032882419676368574984049638008D+00
        w(17) = 0.034518618398549058625221276859D+00
        w(18) = 0.036075423225565273932166270524D+00
        w(19) = 0.037549253448257709809772223198D+00
        w(20) = 0.038936719204051197616673806364D+00
        w(21) = 0.040234629273005533815446337743D+00
        w(22) = 0.041439998417240293022686299233D+00
        w(23) = 0.042550054246755802719217150803D+00
        w(24) = 0.043562243595800486532284821661D+00
        w(25) = 0.044474238395082974427323504000D+00
        w(26) = 0.045283941026300230657128240574D+00
        w(27) = 0.045989489146651696963893390818D+00
        w(28) = 0.046589259972233498302255136790D+00
        w(29) = 0.047081874010454522246006808290D+00
        w(30) = 0.047466198232885503152644458740D+00
        w(31) = 0.047741348681240621559038972227D+00
        w(32) = 0.047906692500495862031347289176D+00
        w(33) = 0.04796184939446661812070762137D+00
        w(34) = 0.047906692500495862031347289176D+00
        w(35) = 0.047741348681240621559038972227D+00
        w(36) = 0.047466198232885503152644458740D+00
        w(37) = 0.047081874010454522246006808290D+00
        w(38) = 0.046589259972233498302255136790D+00
        w(39) = 0.045989489146651696963893390818D+00
        w(40) = 0.045283941026300230657128240574D+00
        w(41) = 0.044474238395082974427323504000D+00
        w(42) = 0.043562243595800486532284821661D+00
        w(43) = 0.042550054246755802719217150803D+00
        w(44) = 0.041439998417240293022686299233D+00
        w(45) = 0.040234629273005533815446337743D+00
        w(46) = 0.038936719204051197616673806364D+00
        w(47) = 0.037549253448257709809772223198D+00
        w(48) = 0.036075423225565273932166270524D+00
        w(49) = 0.034518618398549058625221276859D+00
        w(50) = 0.032882419676368574984049638008D+00
        w(51) = 0.031170590380189142464431845777D+00
        w(52) = 0.029387067789310668062644859210D+00
        w(53) = 0.027535954088450343942499722327D+00
        w(54) = 0.025621506938037758214084978694D+00
        w(55) = 0.023648129691287236698780978994D+00
        w(56) = 0.021620361284934062841654274667D+00
        w(57) = 0.019542865836750062826837429313D+00
        w(58) = 0.017420421997670248495365759969D+00
        w(59) = 0.015257912146448310349265388145D+00
        w(60) = 0.013060311639994846336168342922D+00
        w(61) = 0.0108326787895979686215140551272D+00
        w(62) = 0.0085801482668814598936358121592D+00
        w(63) = 0.0063079425789717545501888719039D+00
        w(64) = 0.0040215241720037363470786599528D+00
        w(65) = 0.0017292582513002508983395851463D+00

      else if ( n .eq. 127 ) then

        x(1) = -0.9998221304153061462673512D+00
        x(2) = -0.9990629343553118951383159D+00
        x(3) = -0.9976975661898046210744170D+00
        x(4) = -0.9957265513520272266354334D+00
        x(5) = -0.9931510492545171473611308D+00
        x(6) = -0.9899726145914841576077867D+00
        x(7) = -0.9861931740169316667104383D+00
        x(8) = -0.9818150208038141100334631D+00
        x(9) = -0.9768408123430703268174439D+00
        x(10) = -0.9712735681615291922889469D+00
        x(11) = -0.9651166679452921210908251D+00
        x(12) = -0.9583738494252387711491029D+00
        x(13) = -0.9510492060778803105479076D+00
        x(14) = -0.9431471846248148273454496D+00
        x(15) = -0.9346725823247379685736349D+00
        x(16) = -0.9256305440562338491274647D+00
        x(17) = -0.9160265591914658093130886D+00
        x(18) = -0.9058664582618213828024613D+00
        x(19) = -0.8951564094170837089690438D+00
        x(20) = -0.8839029146800265699452579D+00
        x(21) = -0.8721128059985607114196375D+00
        x(22) = -0.8597932410977408098120313D+00
        x(23) = -0.8469516991340975984533393D+00
        x(24) = -0.8335959761548995143795572D+00
        x(25) = -0.8197341803650786741551191D+00
        x(26) = -0.8053747272046802146665608D+00
        x(27) = -0.7905263342398137999454500D+00
        x(28) = -0.7751980158702023824449628D+00
        x(29) = -0.7593990778565366715566637D+00
        x(30) = -0.7431391116709545129205669D+00
        x(31) = -0.7264279886740726855356929D+00
        x(32) = -0.7092758541221045609994446D+00
        x(33) = -0.6916931210077006701564414D+00
        x(34) = -0.6736904637382504853466825D+00
        x(35) = -0.6552788116554826302767651D+00
        x(36) = -0.6364693424002972413476082D+00
        x(37) = -0.6172734751268582838576392D+00
        x(38) = -0.5977028635700652293844120D+00
        x(39) = -0.5777693889706125800032517D+00
        x(40) = -0.5574851528619322329218619D+00
        x(41) = -0.5368624697233975674581664D+00
        x(42) = -0.5159138595042493572772773D+00
        x(43) = -0.4946520400227821173949402D+00
        x(44) = -0.4730899192454052416450999D+00
        x(45) = -0.4512405874502662273318986D+00
        x(46) = -0.4291173092801933762625441D+00
        x(47) = -0.4067335156897825634086729D+00
        x(48) = -0.3841027957915169357790778D+00
        x(49) = -0.3612388886058697060709248D+00
        x(50) = -0.3381556747203985013760003D+00
        x(51) = -0.3148671678628949814860148D+00
        x(52) = -0.2913875063937056207945188D+00
        x(53) = -0.2677309447223886208883435D+00
        x(54) = -0.2439118446539178579707132D+00
        x(55) = -0.2199446666696875424545234D+00
        x(56) = -0.1958439611486108515042816D+00
        x(57) = -0.1716243595336421650083449D+00
        x(58) = -0.1473005654490856693893293D+00
        x(59) = -0.1228873457740829717260337D+00
        x(60) = -0.0983995216776989707510918D+00
        x(61) = -0.0738519596210485452734404D+00
        x(62) = -0.0492595623319266303153793D+00
        x(63) = -0.0246372597574209446148971D+00
        x(64) = 0.0000000000000000000000000D+00
        x(65) = 0.0246372597574209446148971D+00
        x(66) = 0.0492595623319266303153793D+00
        x(67) = 0.0738519596210485452734404D+00
        x(68) = 0.0983995216776989707510918D+00
        x(69) = 0.1228873457740829717260337D+00
        x(70) = 0.1473005654490856693893293D+00
        x(71) = 0.1716243595336421650083449D+00
        x(72) = 0.1958439611486108515042816D+00
        x(73) = 0.2199446666696875424545234D+00
        x(74) = 0.2439118446539178579707132D+00
        x(75) = 0.2677309447223886208883435D+00
        x(76) = 0.2913875063937056207945188D+00
        x(77) = 0.3148671678628949814860148D+00
        x(78) = 0.3381556747203985013760003D+00
        x(79) = 0.3612388886058697060709248D+00
        x(80) = 0.3841027957915169357790778D+00
        x(81) = 0.4067335156897825634086729D+00
        x(82) = 0.4291173092801933762625441D+00
        x(83) = 0.4512405874502662273318986D+00
        x(84) = 0.4730899192454052416450999D+00
        x(85) = 0.4946520400227821173949402D+00
        x(86) = 0.5159138595042493572772773D+00
        x(87) = 0.5368624697233975674581664D+00
        x(88) = 0.5574851528619322329218619D+00
        x(89) = 0.5777693889706125800032517D+00
        x(90) = 0.5977028635700652293844120D+00
        x(91) = 0.6172734751268582838576392D+00
        x(92) = 0.6364693424002972413476082D+00
        x(93) = 0.6552788116554826302767651D+00
        x(94) = 0.6736904637382504853466825D+00
        x(95) = 0.6916931210077006701564414D+00
        x(96) = 0.7092758541221045609994446D+00
        x(97) = 0.7264279886740726855356929D+00
        x(98) = 0.7431391116709545129205669D+00
        x(99) = 0.7593990778565366715566637D+00
        x(100) = 0.7751980158702023824449628D+00
        x(101) = 0.7905263342398137999454500D+00
        x(102) = 0.8053747272046802146665608D+00
        x(103) = 0.8197341803650786741551191D+00
        x(104) = 0.8335959761548995143795572D+00
        x(105) = 0.8469516991340975984533393D+00
        x(106) = 0.8597932410977408098120313D+00
        x(107) = 0.8721128059985607114196375D+00
        x(108) = 0.8839029146800265699452579D+00
        x(109) = 0.8951564094170837089690438D+00
        x(110) = 0.9058664582618213828024613D+00
        x(111) = 0.9160265591914658093130886D+00
        x(112) = 0.9256305440562338491274647D+00
        x(113) = 0.9346725823247379685736349D+00
        x(114) = 0.9431471846248148273454496D+00
        x(115) = 0.9510492060778803105479076D+00
        x(116) = 0.9583738494252387711491029D+00
        x(117) = 0.965116667945292121090825D+00
        x(118) = 0.971273568161529192288947D+00
        x(119) = 0.976840812343070326817444D+00
        x(120) = 0.981815020803814110033463D+00
        x(121) = 0.986193174016931666710438D+00
        x(122) = 0.989972614591484157607787D+00
        x(123) = 0.993151049254517147361131D+00
        x(124) = 0.995726551352027226635433D+00
        x(125) = 0.997697566189804621074417D+00
        x(126) = 0.999062934355311895138316D+00
        x(127) = 0.999822130415306146267351D+00

        w(1) = 0.00045645726109586662791936519265D+00
        w(2) = 0.00106227668695384869596523598532D+00
        w(3) = 0.0016683488125171936761028862915D+00
        w(4) = 0.0022734860707492547802810840776D+00
        w(5) = 0.0028772587656289004082883197514D+00
        w(6) = 0.0034792893810051465908910894100D+00
        w(7) = 0.0040792095178254605327114733457D+00
        w(8) = 0.0046766539777779034772638165663D+00
        w(9) = 0.0052712596565634400891303815906D+00
        w(10) = 0.0058626653903523901033648343751D+00
        w(11) = 0.0064505120486899171845442463869D+00
        w(12) = 0.0070344427036681608755685893033D+00
        w(13) = 0.0076141028256526859356393930849D+00
        w(14) = 0.0081891404887415730817235884719D+00
        w(15) = 0.0087592065795403145773316804234D+00
        w(16) = 0.0093239550065309714787536985834D+00
        w(17) = 0.0098830429087554914716648010900D+00
        w(18) = 0.0104361308631410052256731719977D+00
        w(19) = 0.0109828830900689757887996573761D+00
        w(20) = 0.011522967656921087154811609735D+00
        w(21) = 0.012056056679400848183529562145D+00
        w(22) = 0.012581826520465013101514365424D+00
        w(23) = 0.013099957986718627426172681913D+00
        w(24) = 0.013610136522139249906034237534D+00
        w(25) = 0.014112052399003395774044161634D+00
        w(26) = 0.014605400905893418351737288079D+00
        w(27) = 0.015089882532666922992635733981D+00
        w(28) = 0.015565203152273955098532590263D+00
        w(29) = 0.016031074199309941802254151843D+00
        w(30) = 0.016487212845194879399346060358D+00
        w(31) = 0.016933342169871654545878815295D+00
        w(32) = 0.017369191329918731922164721250D+00
        w(33) = 0.017794495722974774231027912900D+00
        w(34) = 0.018208997148375106468721469154D+00
        w(35) = 0.018612443963902310429440419899D+00
        w(36) = 0.019004591238555646611148901045D+00
        w(37) = 0.019385200901246454628112623489D+00
        w(38) = 0.019754041885329183081815217323D+00
        w(39) = 0.020110890268880247225644623956D+00
        w(40) = 0.020455529410639508279497065713D+00
        w(41) = 0.020787750081531811812652137291D+00
        w(42) = 0.021107350591688713643523847922D+00
        w(43) = 0.021414136912893259295449693234D+00
        w(44) = 0.021707922796373466052301324695D+00
        w(45) = 0.021988529885872983756478409759D+00
        w(46) = 0.022255787825930280235631416460D+00
        w(47) = 0.022509534365300608085694429903D+00
        w(48) = 0.022749615455457959852242553241D+00
        w(49) = 0.022975885344117206754377437839D+00
        w(50) = 0.023188206663719640249922582982D+00
        w(51) = 0.023386450514828194170722043497D+00
        w(52) = 0.023570496544381716050033676844D+00
        w(53) = 0.023740233018760777777714726703D+00
        w(54) = 0.023895556891620665983864481754D+00
        w(55) = 0.024036373866450369675132086026D+00
        w(56) = 0.024162598453819584716522917711D+00
        w(57) = 0.024274154023278979833195063937D+00
        w(58) = 0.024370972849882214952813561907D+00
        w(59) = 0.024452996155301467956140198472D+00
        w(60) = 0.024520174143511508275183033290D+00
        w(61) = 0.024572466031020653286354137335D+00
        w(62) = 0.024609840071630254092545634003D+00
        w(63) = 0.024632273575707679066033370218D+00
        w(64) = 0.02463975292396109441957941748D+00
        w(65) = 0.024632273575707679066033370218D+00
        w(66) = 0.024609840071630254092545634003D+00
        w(67) = 0.024572466031020653286354137335D+00
        w(68) = 0.024520174143511508275183033290D+00
        w(69) = 0.024452996155301467956140198472D+00
        w(70) = 0.024370972849882214952813561907D+00
        w(71) = 0.024274154023278979833195063937D+00
        w(72) = 0.024162598453819584716522917711D+00
        w(73) = 0.024036373866450369675132086026D+00
        w(74) = 0.023895556891620665983864481754D+00
        w(75) = 0.023740233018760777777714726703D+00
        w(76) = 0.023570496544381716050033676844D+00
        w(77) = 0.023386450514828194170722043497D+00
        w(78) = 0.023188206663719640249922582982D+00
        w(79) = 0.022975885344117206754377437839D+00
        w(80) = 0.022749615455457959852242553241D+00
        w(81) = 0.022509534365300608085694429903D+00
        w(82) = 0.022255787825930280235631416460D+00
        w(83) = 0.021988529885872983756478409759D+00
        w(84) = 0.021707922796373466052301324695D+00
        w(85) = 0.021414136912893259295449693234D+00
        w(86) = 0.021107350591688713643523847922D+00
        w(87) = 0.020787750081531811812652137291D+00
        w(88) = 0.020455529410639508279497065713D+00
        w(89) = 0.020110890268880247225644623956D+00
        w(90) = 0.019754041885329183081815217323D+00
        w(91) = 0.019385200901246454628112623489D+00
        w(92) = 0.019004591238555646611148901045D+00
        w(93) = 0.018612443963902310429440419899D+00
        w(94) = 0.018208997148375106468721469154D+00
        w(95) = 0.017794495722974774231027912900D+00
        w(96) = 0.017369191329918731922164721250D+00
        w(97) = 0.016933342169871654545878815295D+00
        w(98) = 0.016487212845194879399346060358D+00
        w(99) = 0.016031074199309941802254151843D+00
        w(100) = 0.015565203152273955098532590263D+00
        w(101) = 0.015089882532666922992635733981D+00
        w(102) = 0.014605400905893418351737288079D+00
        w(103) = 0.014112052399003395774044161634D+00
        w(104) = 0.013610136522139249906034237534D+00
        w(105) = 0.013099957986718627426172681913D+00
        w(106) = 0.012581826520465013101514365424D+00
        w(107) = 0.012056056679400848183529562145D+00
        w(108) = 0.011522967656921087154811609735D+00
        w(109) = 0.0109828830900689757887996573761D+00
        w(110) = 0.0104361308631410052256731719977D+00
        w(111) = 0.0098830429087554914716648010900D+00
        w(112) = 0.0093239550065309714787536985834D+00
        w(113) = 0.0087592065795403145773316804234D+00
        w(114) = 0.0081891404887415730817235884719D+00
        w(115) = 0.0076141028256526859356393930849D+00
        w(116) = 0.0070344427036681608755685893033D+00
        w(117) = 0.0064505120486899171845442463869D+00
        w(118) = 0.0058626653903523901033648343751D+00
        w(119) = 0.0052712596565634400891303815906D+00
        w(120) = 0.0046766539777779034772638165663D+00
        w(121) = 0.0040792095178254605327114733457D+00
        w(122) = 0.0034792893810051465908910894100D+00
        w(123) = 0.0028772587656289004082883197514D+00
        w(124) = 0.0022734860707492547802810840776D+00
        w(125) = 0.0016683488125171936761028862915D+00
        w(126) = 0.00106227668695384869596523598532D+00
        w(127) = 0.00045645726109586662791936519265D+00

      else if ( n .eq. 128 ) then

        x(1) = -0.9998248879471319144736081D+00
        x(2) = -0.9990774599773758950119878D+00
        x(3) = -0.9977332486255140198821574D+00
        x(4) = -0.9957927585349811868641612D+00
        x(5) = -0.9932571129002129353034372D+00
        x(6) = -0.9901278184917343833379303D+00
        x(7) = -0.9864067427245862088712355D+00
        x(8) = -0.9820961084357185360247656D+00
        x(9) = -0.9771984914639073871653744D+00
        x(10) = -0.9717168187471365809043384D+00
        x(11) = -0.9656543664319652686458290D+00
        x(12) = -0.9590147578536999280989185D+00
        x(13) = -0.9518019613412643862177963D+00
        x(14) = -0.9440202878302201821211114D+00
        x(15) = -0.9356743882779163757831268D+00
        x(16) = -0.9267692508789478433346245D+00
        x(17) = -0.9173101980809605370364836D+00
        x(18) = -0.9073028834017568139214859D+00
        x(19) = -0.8967532880491581843864474D+00
        x(20) = -0.8856677173453972174082924D+00
        x(21) = -0.8740527969580317986954180D+00
        x(22) = -0.8619154689395484605906323D+00
        x(23) = -0.8492629875779689691636001D+00
        x(24) = -0.8361029150609068471168753D+00
        x(25) = -0.8224431169556438424645942D+00
        x(26) = -0.8082917575079136601196422D+00
        x(27) = -0.7936572947621932902433329D+00
        x(28) = -0.7785484755064119668504941D+00
        x(29) = -0.7629743300440947227797691D+00
        x(30) = -0.7469441667970619811698824D+00
        x(31) = -0.7304675667419088064717369D+00
        x(32) = -0.7135543776835874133438599D+00
        x(33) = -0.6962147083695143323850866D+00
        x(34) = -0.6784589224477192593677557D+00
        x(35) = -0.6602976322726460521059468D+00
        x(36) = -0.6417416925623075571535249D+00
        x(37) = -0.6228021939105849107615396D+00
        x(38) = -0.6034904561585486242035732D+00
        x(39) = -0.5838180216287630895500389D+00
        x(40) = -0.5637966482266180839144308D+00
        x(41) = -0.5434383024128103634441936D+00
        x(42) = -0.5227551520511754784539479D+00
        x(43) = -0.5017595591361444642896063D+00
        x(44) = -0.4804640724041720258582757D+00
        x(45) = -0.4588814198335521954490891D+00
        x(46) = -0.4370245010371041629370429D+00
        x(47) = -0.4149063795522750154922739D+00
        x(48) = -0.3925402750332674427356482D+00
        x(49) = -0.3699395553498590266165917D+00
        x(50) = -0.3471177285976355084261628D+00
        x(51) = -0.3240884350244133751832523D+00
        x(52) = -0.3008654388776772026671541D+00
        x(53) = -0.2774626201779044028062316D+00
        x(54) = -0.2538939664226943208556180D+00
        x(55) = -0.2301735642266599864109866D+00
        x(56) = -0.2063155909020792171540580D+00
        x(57) = -0.1823343059853371824103826D+00
        x(58) = -0.1582440427142249339974755D+00
        x(59) = -0.1340591994611877851175753D+00
        x(60) = -0.1097942311276437466729747D+00
        x(61) = -0.0854636405045154986364980D+00
        x(62) = -0.0610819696041395681037870D+00
        x(63) = -0.0366637909687334933302153D+00
        x(64) = -0.0122236989606157641980521D+00
        x(65) = 0.0122236989606157641980521D+00
        x(66) = 0.0366637909687334933302153D+00
        x(67) = 0.0610819696041395681037870D+00
        x(68) = 0.0854636405045154986364980D+00
        x(69) = 0.1097942311276437466729747D+00
        x(70) = 0.1340591994611877851175753D+00
        x(71) = 0.1582440427142249339974755D+00
        x(72) = 0.1823343059853371824103826D+00
        x(73) = 0.2063155909020792171540580D+00
        x(74) = 0.2301735642266599864109866D+00
        x(75) = 0.2538939664226943208556180D+00
        x(76) = 0.2774626201779044028062316D+00
        x(77) = 0.3008654388776772026671541D+00
        x(78) = 0.3240884350244133751832523D+00
        x(79) = 0.3471177285976355084261628D+00
        x(80) = 0.3699395553498590266165917D+00
        x(81) = 0.3925402750332674427356482D+00
        x(82) = 0.4149063795522750154922739D+00
        x(83) = 0.4370245010371041629370429D+00
        x(84) = 0.4588814198335521954490891D+00
        x(85) = 0.4804640724041720258582757D+00
        x(86) = 0.5017595591361444642896063D+00
        x(87) = 0.5227551520511754784539479D+00
        x(88) = 0.5434383024128103634441936D+00
        x(89) = 0.5637966482266180839144308D+00
        x(90) = 0.5838180216287630895500389D+00
        x(91) = 0.6034904561585486242035732D+00
        x(92) = 0.6228021939105849107615396D+00
        x(93) = 0.6417416925623075571535249D+00
        x(94) = 0.6602976322726460521059468D+00
        x(95) = 0.6784589224477192593677557D+00
        x(96) = 0.6962147083695143323850866D+00
        x(97) = 0.7135543776835874133438599D+00
        x(98) = 0.7304675667419088064717369D+00
        x(99) = 0.7469441667970619811698824D+00
        x(100) = 0.7629743300440947227797691D+00
        x(101) = 0.7785484755064119668504941D+00
        x(102) = 0.7936572947621932902433329D+00
        x(103) = 0.8082917575079136601196422D+00
        x(104) = 0.8224431169556438424645942D+00
        x(105) = 0.8361029150609068471168753D+00
        x(106) = 0.8492629875779689691636001D+00
        x(107) = 0.8619154689395484605906323D+00
        x(108) = 0.8740527969580317986954180D+00
        x(109) = 0.8856677173453972174082924D+00
        x(110) = 0.8967532880491581843864474D+00
        x(111) = 0.9073028834017568139214859D+00
        x(112) = 0.9173101980809605370364836D+00
        x(113) = 0.926769250878947843334625D+00
        x(114) = 0.935674388277916375783127D+00
        x(115) = 0.944020287830220182121111D+00
        x(116) = 0.951801961341264386217796D+00
        x(117) = 0.959014757853699928098919D+00
        x(118) = 0.965654366431965268645829D+00
        x(119) = 0.971716818747136580904338D+00
        x(120) = 0.977198491463907387165374D+00
        x(121) = 0.982096108435718536024766D+00
        x(122) = 0.986406742724586208871236D+00
        x(123) = 0.990127818491734383337930D+00
        x(124) = 0.993257112900212935303437D+00
        x(125) = 0.995792758534981186864161D+00
        x(126) = 0.997733248625514019882157D+00
        x(127) = 0.999077459977375895011988D+00
        x(128) = 0.999824887947131914473608D+00

        w(1) = 0.00044938096029209037639429223999D+00
        w(2) = 0.0010458126793403487793128516001D+00
        w(3) = 0.0016425030186690295387908755948D+00
        w(4) = 0.0022382884309626187436220542727D+00
        w(5) = 0.0028327514714579910952857346468D+00
        w(6) = 0.0034255260409102157743377846601D+00
        w(7) = 0.0040162549837386423131943434863D+00
        w(8) = 0.0046045842567029551182905419803D+00
        w(9) = 0.0051901618326763302050707671348D+00
        w(10) = 0.0057726375428656985893346176261D+00
        w(11) = 0.006351663161707188787214327826D+00
        w(12) = 0.006926892566898813563426670360D+00
        w(13) = 0.007497981925634728687671962688D+00
        w(14) = 0.008064589890486057972928598698D+00
        w(15) = 0.008626377798616749704978843782D+00
        w(16) = 0.009183009871660874334478743688D+00
        w(17) = 0.009734153415006805863548266094D+00
        w(18) = 0.010279479015832157133215340326D+00
        w(19) = 0.010818660739503076247659646277D+00
        w(20) = 0.011351376324080416693281668453D+00
        w(21) = 0.011877307372740279575891106926D+00
        w(22) = 0.012396139543950922968821728197D+00
        w(23) = 0.012907562739267347220442834004D+00
        w(24) = 0.013411271288616332314488951616D+00
        w(25) = 0.013906964132951985244288007396D+00
        w(26) = 0.014394345004166846176823892009D+00
        w(27) = 0.014873122602147314252385498520D+00
        w(28) = 0.015343010768865144085990853741D+00
        w(29) = 0.015803728659399346858965631687D+00
        w(30) = 0.016255000909785187051657456477D+00
        w(31) = 0.016696557801589204589091507954D+00
        w(32) = 0.017128135423111376830680987619D+00
        w(33) = 0.017549475827117704648706925634D+00
        w(34) = 0.017960327185008685940196927525D+00
        w(35) = 0.018360443937331343221289290991D+00
        w(36) = 0.018749586940544708650919548474D+00
        w(37) = 0.019127523609950945486518531668D+00
        w(38) = 0.019494028058706602823021918681D+00
        w(39) = 0.019848881232830862219944413265D+00
        w(40) = 0.020191871042130041180673158406D+00
        w(41) = 0.020522792486960069432284967788D+00
        w(42) = 0.020841447780751149113583948423D+00
        w(43) = 0.021147646468221348537019535180D+00
        w(44) = 0.021441205539208460137111853878D+00
        w(45) = 0.021721949538052075375260957768D+00
        w(46) = 0.021989710668460491434122106599D+00
        w(47) = 0.022244328893799765104629133607D+00
        w(48) = 0.022485652032744966871824603941D+00
        w(49) = 0.022713535850236461309712635923D+00
        w(50) = 0.022927844143686846920410987209D+00
        w(51) = 0.023128448824387027879297902403D+00
        w(52) = 0.023315229994062760122415671273D+00
        w(53) = 0.023488076016535913153025273282D+00
        w(54) = 0.023646883584447615143651392303D+00
        w(55) = 0.023791557781003400638780709885D+00
        w(56) = 0.023922012136703455672450408817D+00
        w(57) = 0.024038168681024052637587316820D+00
        w(58) = 0.024139957989019284997716653890D+00
        w(59) = 0.024227319222815248120093308442D+00
        w(60) = 0.024300200167971865323442606364D+00
        w(61) = 0.024358557264690625853268520246D+00
        w(62) = 0.024402355633849582093297989694D+00
        w(63) = 0.02443156909785004505484856143D+00
        w(64) = 0.02444618019626251821132585261D+00
        w(65) = 0.02444618019626251821132585261D+00
        w(66) = 0.02443156909785004505484856143D+00
        w(67) = 0.024402355633849582093297989694D+00
        w(68) = 0.024358557264690625853268520246D+00
        w(69) = 0.024300200167971865323442606364D+00
        w(70) = 0.024227319222815248120093308442D+00
        w(71) = 0.024139957989019284997716653890D+00
        w(72) = 0.024038168681024052637587316820D+00
        w(73) = 0.023922012136703455672450408817D+00
        w(74) = 0.023791557781003400638780709885D+00
        w(75) = 0.023646883584447615143651392303D+00
        w(76) = 0.023488076016535913153025273282D+00
        w(77) = 0.023315229994062760122415671273D+00
        w(78) = 0.023128448824387027879297902403D+00
        w(79) = 0.022927844143686846920410987209D+00
        w(80) = 0.022713535850236461309712635923D+00
        w(81) = 0.022485652032744966871824603941D+00
        w(82) = 0.022244328893799765104629133607D+00
        w(83) = 0.021989710668460491434122106599D+00
        w(84) = 0.021721949538052075375260957768D+00
        w(85) = 0.021441205539208460137111853878D+00
        w(86) = 0.021147646468221348537019535180D+00
        w(87) = 0.020841447780751149113583948423D+00
        w(88) = 0.020522792486960069432284967788D+00
        w(89) = 0.020191871042130041180673158406D+00
        w(90) = 0.019848881232830862219944413265D+00
        w(91) = 0.019494028058706602823021918681D+00
        w(92) = 0.019127523609950945486518531668D+00
        w(93) = 0.018749586940544708650919548474D+00
        w(94) = 0.018360443937331343221289290991D+00
        w(95) = 0.017960327185008685940196927525D+00
        w(96) = 0.017549475827117704648706925634D+00
        w(97) = 0.017128135423111376830680987619D+00
        w(98) = 0.016696557801589204589091507954D+00
        w(99) = 0.016255000909785187051657456477D+00
        w(100) = 0.015803728659399346858965631687D+00
        w(101) = 0.015343010768865144085990853741D+00
        w(102) = 0.014873122602147314252385498520D+00
        w(103) = 0.014394345004166846176823892009D+00
        w(104) = 0.013906964132951985244288007396D+00
        w(105) = 0.013411271288616332314488951616D+00
        w(106) = 0.012907562739267347220442834004D+00
        w(107) = 0.012396139543950922968821728197D+00
        w(108) = 0.011877307372740279575891106926D+00
        w(109) = 0.011351376324080416693281668453D+00
        w(110) = 0.010818660739503076247659646277D+00
        w(111) = 0.010279479015832157133215340326D+00
        w(112) = 0.009734153415006805863548266094D+00
        w(113) = 0.009183009871660874334478743688D+00
        w(114) = 0.008626377798616749704978843782D+00
        w(115) = 0.008064589890486057972928598698D+00
        w(116) = 0.007497981925634728687671962688D+00
        w(117) = 0.006926892566898813563426670360D+00
        w(118) = 0.006351663161707188787214327826D+00
        w(119) = 0.0057726375428656985893346176261D+00
        w(120) = 0.0051901618326763302050707671348D+00
        w(121) = 0.0046045842567029551182905419803D+00
        w(122) = 0.0040162549837386423131943434863D+00
        w(123) = 0.0034255260409102157743377846601D+00
        w(124) = 0.0028327514714579910952857346468D+00
        w(125) = 0.0022382884309626187436220542727D+00
        w(126) = 0.0016425030186690295387908755948D+00
        w(127) = 0.0010458126793403487793128516001D+00
        w(128) = 0.00044938096029209037639429223999D+00

      else if ( n .eq. 129 ) then

        x(1) = -0.9998275818477487191077441D+00
        x(2) = -0.9990916504696409986514389D+00
        x(3) = -0.9977681080525852721429460D+00
        x(4) = -0.9958574393142831982149111D+00
        x(5) = -0.9933607326210712814854011D+00
        x(6) = -0.9902794486488178389207689D+00
        x(7) = -0.9866153978313475022005761D+00
        x(8) = -0.9823707352517413115507418D+00
        x(9) = -0.9775479582993672474447814D+00
        x(10) = -0.9721499048427034297274163D+00
        x(11) = -0.9661797514202097197778763D+00
        x(12) = -0.9596410113101918904168119D+00
        x(13) = -0.9525375324342090471027732D+00
        x(14) = -0.9448734950776734726784764D+00
        x(15) = -0.9366534094216514605284616D+00
        x(16) = -0.9278821128840036204317296D+00
        x(17) = -0.9185647672698286252225115D+00
        x(18) = -0.9087068557320696331245539D+00
        x(19) = -0.8983141795436338850435985D+00
        x(20) = -0.8873928546826803665034968D+00
        x(21) = -0.8759493082329433892035217D+00
        x(22) = -0.8639902746011257878940216D+00
        x(23) = -0.8515227915535356930243826D+00
        x(24) = -0.8385541960742664442975407D+00
        x(25) = -0.8250921200473358809210133D+00
        x(26) = -0.8111444857653120742087717D+00
        x(27) = -0.7967195012670592680339606D+00
        x(28) = -0.7818256555073413245387500D+00
        x(29) = -0.7664717133611208816717785D+00
        x(30) = -0.7506667104654910227632368D+00
        x(31) = -0.7344199479022727047791516D+00
        x(32) = -0.7177409867244055767721220D+00
        x(33) = -0.7006396423293521790044710D+00
        x(34) = -0.6831259786828258512462248D+00
        x(35) = -0.6652103023962409818802202D+00
        x(36) = -0.6469031566613704719753373D+00
        x(37) = -0.6282153150457794374886895D+00
        x(38) = -0.6091577751526861909563306D+00
        x(39) = -0.5897417521489813916767844D+00
        x(40) = -0.5699786721652138894754096D+00
        x(41) = -0.5498801655714271702189358D+00
        x(42) = -0.5294580601328034000099406D+00
        x(43) = -0.5087243740491428186199463D+00
        x(44) = -0.4876913088822746111853066D+00
        x(45) = -0.4663712423755613514331869D+00
        x(46) = -0.4447767211697226217818454D+00
        x(47) = -0.4229204534192644388475065D+00
        x(48) = -0.4008153013138596117693121D+00
        x(49) = -0.3784742735090801012801265D+00
        x(50) = -0.3559105174709357969672656D+00
        x(51) = -0.3331373117387248575049982D+00
        x(52) = -0.3101680581107488341147318D+00
        x(53) = -0.2870162737574911929568755D+00
        x(54) = -0.2636955832669005409666949D+00
        x(55) = -0.2402197106264598167721148D+00
        x(56) = -0.2166024711467599103221439D+00
        x(57) = -0.1928577633313305998663880D+00
        x(58) = -0.1689995606975133227390302D+00
        x(59) = -0.1450419035531891084328306D+00
        x(60) = -0.1209988907342009817690539D+00
        x(61) = -0.0968846713073332753086909D+00
        x(62) = -0.0727134362437305599118207D+00
        x(63) = -0.0484994100676562986191764D+00
        x(64) = -0.0242568424855058415749954D+00
        x(65) = 0.0000000000000000000000000D+00
        x(66) = 0.0242568424855058415749954D+00
        x(67) = 0.0484994100676562986191764D+00
        x(68) = 0.0727134362437305599118207D+00
        x(69) = 0.0968846713073332753086909D+00
        x(70) = 0.1209988907342009817690539D+00
        x(71) = 0.1450419035531891084328306D+00
        x(72) = 0.1689995606975133227390302D+00
        x(73) = 0.1928577633313305998663880D+00
        x(74) = 0.2166024711467599103221439D+00
        x(75) = 0.2402197106264598167721148D+00
        x(76) = 0.2636955832669005409666949D+00
        x(77) = 0.2870162737574911929568755D+00
        x(78) = 0.3101680581107488341147318D+00
        x(79) = 0.3331373117387248575049982D+00
        x(80) = 0.3559105174709357969672656D+00
        x(81) = 0.3784742735090801012801265D+00
        x(82) = 0.4008153013138596117693121D+00
        x(83) = 0.4229204534192644388475065D+00
        x(84) = 0.4447767211697226217818454D+00
        x(85) = 0.4663712423755613514331869D+00
        x(86) = 0.4876913088822746111853066D+00
        x(87) = 0.5087243740491428186199463D+00
        x(88) = 0.5294580601328034000099406D+00
        x(89) = 0.5498801655714271702189358D+00
        x(90) = 0.5699786721652138894754096D+00
        x(91) = 0.5897417521489813916767844D+00
        x(92) = 0.6091577751526861909563306D+00
        x(93) = 0.6282153150457794374886895D+00
        x(94) = 0.6469031566613704719753373D+00
        x(95) = 0.6652103023962409818802202D+00
        x(96) = 0.6831259786828258512462248D+00
        x(97) = 0.7006396423293521790044710D+00
        x(98) = 0.7177409867244055767721220D+00
        x(99) = 0.7344199479022727047791516D+00
        x(100) = 0.7506667104654910227632368D+00
        x(101) = 0.7664717133611208816717785D+00
        x(102) = 0.7818256555073413245387500D+00
        x(103) = 0.7967195012670592680339606D+00
        x(104) = 0.8111444857653120742087717D+00
        x(105) = 0.8250921200473358809210133D+00
        x(106) = 0.8385541960742664442975407D+00
        x(107) = 0.8515227915535356930243826D+00
        x(108) = 0.8639902746011257878940216D+00
        x(109) = 0.875949308232943389203522D+00
        x(110) = 0.887392854682680366503497D+00
        x(111) = 0.898314179543633885043599D+00
        x(112) = 0.908706855732069633124554D+00
        x(113) = 0.918564767269828625222511D+00
        x(114) = 0.927882112884003620431730D+00
        x(115) = 0.936653409421651460528462D+00
        x(116) = 0.944873495077673472678476D+00
        x(117) = 0.952537532434209047102773D+00
        x(118) = 0.959641011310191890416812D+00
        x(119) = 0.966179751420209719777876D+00
        x(120) = 0.972149904842703429727416D+00
        x(121) = 0.977547958299367247444781D+00
        x(122) = 0.982370735251741311550742D+00
        x(123) = 0.986615397831347502200576D+00
        x(124) = 0.990279448648817838920769D+00
        x(125) = 0.993360732621071281485401D+00
        x(126) = 0.995857439314283198214911D+00
        x(127) = 0.997768108052585272142946D+00
        x(128) = 0.999091650469640998651439D+00
        x(129) = 0.999827581847748719107744D+00

        w(1) = 0.00044246794182939296923668005717D+00
        w(2) = 0.00102972844619622394463273519315D+00
        w(3) = 0.0016172530556785534682413679271D+00
        w(4) = 0.0022039015180966937075786419741D+00
        w(5) = 0.0027892681877797554940944677057D+00
        w(6) = 0.0033729979506246246117755709288D+00
        w(7) = 0.0039547444682113562172392974765D+00
        w(8) = 0.0045341644298525434513226874954D+00
        w(9) = 0.0051109164669246267289761565766D+00
        w(10) = 0.0056846609912469045788016012203D+00
        w(11) = 0.0062550602724461408889348709586D+00
        w(12) = 0.0068217785893519121070498527769D+00
        w(13) = 0.0073844824072454014447165055698D+00
        w(14) = 0.0079428405646668029041114107832D+00
        w(15) = 0.0084965244635723279730542832506D+00
        w(16) = 0.0090452082602137316404219313819D+00
        w(17) = 0.0095885690555104190787301294510D+00
        w(18) = 0.0101262870842733548093160774580D+00
        w(19) = 0.0106580459029055185304204093001D+00
        w(20) = 0.0111835325753305049735380697538D+00
        w(21) = 0.011702437856964778185746436834D+00
        w(22) = 0.012214456376582979416221105914D+00
        w(23) = 0.012719286815944623465099036330D+00
        w(24) = 0.013216632087061724231482387345D+00
        w(25) = 0.013706199506993971244060563234D+00
        w(26) = 0.014187700970062900419317230938D+00
        w(27) = 0.014660853117380060971041027493D+00
        w(28) = 0.015125377503587024690403432771D+00
        w(29) = 0.015581000760707523415881287558D+00
        w(30) = 0.016027454759014214436403950465D+00
        w(31) = 0.016464476764814667467169189640D+00
        w(32) = 0.016891809595063204177526208819D+00
        w(33) = 0.017309201768707240731293596444D+00
        w(34) = 0.017716407654678809269702031810D+00
        w(35) = 0.018113187616443980503999783812D+00
        w(36) = 0.018499308153024985727791918518D+00
        w(37) = 0.018874542036411948181617592169D+00
        w(38) = 0.019238668445283284085199492202D+00
        w(39) = 0.019591473094956024580283987216D+00
        w(40) = 0.019932748363489542089706675388D+00
        w(41) = 0.020262293413868438317104423081D+00
        w(42) = 0.020579914312192665948185517085D+00
        w(43) = 0.020885424141805311409990024684D+00
        w(44) = 0.021178643113290860912881038703D+00
        w(45) = 0.021459398670279205389981598196D+00
        w(46) = 0.021727525590993110687305178710D+00
        w(47) = 0.021982866085479386179554968899D+00
        w(48) = 0.022225269888466526554736910919D+00
        w(49) = 0.022454594347794176432066564511D+00
        w(50) = 0.022670704508362374313093970958D+00
        w(51) = 0.022873473191551169638592083492D+00
        w(52) = 0.023062781070063872924670495006D+00
        w(53) = 0.023238516738149892544490435771D+00
        w(54) = 0.023400576777165831146714346635D+00
        w(55) = 0.023548865816436258377269094263D+00
        w(56) = 0.023683296589378342897341543485D+00
        w(57) = 0.023803789984857314051325299744D+00
        w(58) = 0.023910275093742530302367230296D+00
        w(59) = 0.024002689250636756075547029720D+00
        w(60) = 0.024080978070754089272959634041D+00
        w(61) = 0.024145095481924836783843156014D+00
        w(62) = 0.024195003751708503129818111597D+00
        w(63) = 0.024230673509598936275508460625D+00
        w(64) = 0.024252083764308562906498864071D+00
        w(65) = 0.02425922191612154143202867472D+00
        w(66) = 0.024252083764308562906498864071D+00
        w(67) = 0.024230673509598936275508460625D+00
        w(68) = 0.024195003751708503129818111597D+00
        w(69) = 0.024145095481924836783843156014D+00
        w(70) = 0.024080978070754089272959634041D+00
        w(71) = 0.024002689250636756075547029720D+00
        w(72) = 0.023910275093742530302367230296D+00
        w(73) = 0.023803789984857314051325299744D+00
        w(74) = 0.023683296589378342897341543485D+00
        w(75) = 0.023548865816436258377269094263D+00
        w(76) = 0.023400576777165831146714346635D+00
        w(77) = 0.023238516738149892544490435771D+00
        w(78) = 0.023062781070063872924670495006D+00
        w(79) = 0.022873473191551169638592083492D+00
        w(80) = 0.022670704508362374313093970958D+00
        w(81) = 0.022454594347794176432066564511D+00
        w(82) = 0.022225269888466526554736910919D+00
        w(83) = 0.021982866085479386179554968899D+00
        w(84) = 0.021727525590993110687305178710D+00
        w(85) = 0.021459398670279205389981598196D+00
        w(86) = 0.021178643113290860912881038703D+00
        w(87) = 0.020885424141805311409990024684D+00
        w(88) = 0.020579914312192665948185517085D+00
        w(89) = 0.020262293413868438317104423081D+00
        w(90) = 0.019932748363489542089706675388D+00
        w(91) = 0.019591473094956024580283987216D+00
        w(92) = 0.019238668445283284085199492202D+00
        w(93) = 0.018874542036411948181617592169D+00
        w(94) = 0.018499308153024985727791918518D+00
        w(95) = 0.018113187616443980503999783812D+00
        w(96) = 0.017716407654678809269702031810D+00
        w(97) = 0.017309201768707240731293596444D+00
        w(98) = 0.016891809595063204177526208819D+00
        w(99) = 0.016464476764814667467169189640D+00
        w(100) = 0.016027454759014214436403950465D+00
        w(101) = 0.015581000760707523415881287558D+00
        w(102) = 0.015125377503587024690403432771D+00
        w(103) = 0.014660853117380060971041027493D+00
        w(104) = 0.014187700970062900419317230938D+00
        w(105) = 0.013706199506993971244060563234D+00
        w(106) = 0.013216632087061724231482387345D+00
        w(107) = 0.012719286815944623465099036330D+00
        w(108) = 0.012214456376582979416221105914D+00
        w(109) = 0.011702437856964778185746436834D+00
        w(110) = 0.0111835325753305049735380697538D+00
        w(111) = 0.0106580459029055185304204093001D+00
        w(112) = 0.0101262870842733548093160774580D+00
        w(113) = 0.0095885690555104190787301294510D+00
        w(114) = 0.0090452082602137316404219313819D+00
        w(115) = 0.0084965244635723279730542832506D+00
        w(116) = 0.0079428405646668029041114107832D+00
        w(117) = 0.0073844824072454014447165055698D+00
        w(118) = 0.0068217785893519121070498527769D+00
        w(119) = 0.0062550602724461408889348709586D+00
        w(120) = 0.0056846609912469045788016012203D+00
        w(121) = 0.0051109164669246267289761565766D+00
        w(122) = 0.0045341644298525434513226874954D+00
        w(123) = 0.0039547444682113562172392974765D+00
        w(124) = 0.0033729979506246246117755709288D+00
        w(125) = 0.0027892681877797554940944677057D+00
        w(126) = 0.0022039015180966937075786419741D+00
        w(127) = 0.0016172530556785534682413679271D+00
        w(128) = 0.00102972844619622394463273519315D+00
        w(129) = 0.00044246794182939296923668005717D+00

      else if ( n .eq. 255 ) then

        x(1) = -0.999955705317563751730191D+00
        x(2) = -0.999766621312000569367063D+00
        x(3) = -0.999426474680169959344386D+00
        x(4) = -0.998935241284654635142155D+00
        x(5) = -0.998292986136967889228248D+00
        x(6) = -0.997499804126615814044844D+00
        x(7) = -0.996555814435198617028738D+00
        x(8) = -0.995461159480026294089975D+00
        x(9) = -0.994216004616630164799381D+00
        x(10) = -0.992820538021989138984811D+00
        x(11) = -0.991274970630385567164523D+00
        x(12) = -0.989579536085920123498574D+00
        x(13) = -0.987734490699732356281248D+00
        x(14) = -0.985740113407419277752900D+00
        x(15) = -0.983596705724776358640192D+00
        x(16) = -0.981304591701017185126565D+00
        x(17) = -0.978864117869068155239121D+00
        x(18) = -0.976275653192735980815246D+00
        x(19) = -0.973539589010643617645393D+00
        x(20) = -0.970656338976880365477697D+00
        x(21) = -0.967626338998338798105523D+00
        x(22) = -0.964450047168726298761719D+00
        x(23) = -0.961127943699247839572910D+00
        x(24) = -0.957660530845962076295490D+00
        x(25) = -0.954048332833816317950921D+00
        x(26) = -0.950291895777368285733522D+00
        x(27) = -0.946391787598204251752103D+00
        x(28) = -0.942348597939064408301480D+00
        x(29) = -0.938162938074687317626793D+00
        x(30) = -0.933835440819386124349338D+00
        x(31) = -0.929366760431369935739045D+00
        x(32) = -0.924757572513824425220425D+00
        x(33) = -0.920008573912766315142721D+00
        x(34) = -0.915120482611686961035103D+00
        x(35) = -0.910094037623000801254172D+00
        x(36) = -0.904929998876314959753358D+00
        x(37) = -0.899629147103536800144342D+00
        x(38) = -0.894192283720836729335637D+00
        x(39) = -0.888620230707484040924981D+00
        x(40) = -0.882913830481574073645470D+00
        x(41) = -0.877073945772665439532627D+00
        x(42) = -0.871101459491346550796200D+00
        x(43) = -0.864997274595751144137121D+00
        x(44) = -0.858762313955042966785823D+00
        x(45) = -0.852397520209890250084237D+00
        x(46) = -0.845903855629951054143931D+00
        x(47) = -0.839282301968391021084600D+00
        x(48) = -0.832533860313455524647230D+00
        x(49) = -0.825659550937118650611534D+00
        x(50) = -0.818660413140831885432406D+00
        x(51) = -0.811537505098395829833580D+00
        x(52) = -0.804291903695978689734633D+00
        x(53) = -0.796924704369305728807154D+00
        x(54) = -0.789437020938044295117764D+00
        x(55) = -0.781829985437409458675147D+00
        x(56) = -0.774104747947015717207115D+00
        x(57) = -0.766262476417000644100858D+00
        x(58) = -0.758304356491446765092016D+00
        x(59) = -0.750231591329128358931528D+00
        x(60) = -0.742045401421610281838045D+00
        x(61) = -0.733747024408726316001889D+00
        x(62) = -0.725337714891464938687812D+00
        x(63) = -0.716818744242290800531501D+00
        x(64) = -0.708191400412930589382399D+00
        x(65) = -0.699456987739652339456557D+00
        x(66) = -0.690616826746067624571761D+00
        x(67) = -0.681672253943486448787259D+00
        x(68) = -0.672624621628855017806731D+00
        x(69) = -0.663475297680306939970658D+00
        x(70) = -0.654225665350358766508700D+00
        x(71) = -0.644877123056781136890077D+00
        x(72) = -0.635431084171177146547142D+00
        x(73) = -0.625888976805299900901619D+00
        x(74) = -0.616252243595141561442344D+00
        x(75) = -0.606522341482826526536576D+00
        x(76) = -0.596700741496341721653202D+00
        x(77) = -0.586788928527137300685706D+00
        x(78) = -0.576788401105631382036211D+00
        x(79) = -0.566700671174652760010815D+00
        x(80) = -0.556527263860855843833077D+00
        x(81) = -0.546269717244142383159817D+00
        x(82) = -0.535929582125124840335150D+00
        x(83) = -0.525508421790666565699453D+00
        x(84) = -0.515007811777534223035005D+00
        x(85) = -0.504429339634198197635551D+00
        x(86) = -0.493774604680816999489812D+00
        x(87) = -0.483045217767441948626854D+00
        x(88) = -0.472242801030478698742627D+00
        x(89) = -0.461368987647442418771401D+00
        x(90) = -0.450425421590043710043279D+00
        x(91) = -0.439413757375642589040685D+00
        x(92) = -0.428335659817108112494341D+00
        x(93) = -0.417192803771121462605751D+00
        x(94) = -0.405986873884960545511889D+00
        x(95) = -0.394719564341804385683361D+00
        x(96) = -0.383392578604595822734854D+00
        x(97) = -0.372007629158501235092510D+00
        x(98) = -0.360566437252006227074021D+00
        x(99) = -0.349070732636686422161576D+00
        x(100) = -0.337522253305692705554261D+00
        x(101) = -0.325922745230990453444769D+00
        x(102) = -0.314273962099392474845918D+00
        x(103) = -0.302577665047425574167140D+00
        x(104) = -0.290835622395070819082047D+00
        x(105) = -0.279049609378417768508970D+00
        x(106) = -0.267221407881273079721012D+00
        x(107) = -0.255352806165764071686080D+00
        x(108) = -0.243445598601977973686482D+00
        x(109) = -0.231501585396677734059116D+00
        x(110) = -0.219522572321135403508985D+00
        x(111) = -0.207510370438124240859625D+00
        x(112) = -0.195466795828110816293869D+00
        x(113) = -0.183393669314688508087976D+00
        x(114) = -0.171292816189293903533225D+00
        x(115) = -0.159166065935247723154292D+00
        x(116) = -0.147015251951161989456661D+00
        x(117) = -0.134842211273755257250625D+00
        x(118) = -0.122648784300117812092492D+00
        x(119) = -0.110436814509468826540991D+00
        x(120) = -0.098208148184447540736015D+00
        x(121) = -0.085964634131980604256000D+00
        x(122) = -0.073708123403767780288977D+00
        x(123) = -0.061440469016428270850728D+00
        x(124) = -0.049163525671349973093019D+00
        x(125) = -0.036879149474284021657652D+00
        x(126) = -0.024589197654727010541405D+00
        x(127) = -0.012295528285133320036860D+00
        x(128) = 0.000000000000000000000000D+00
        x(129) = 0.012295528285133320036860D+00
        x(130) = 0.024589197654727010541405D+00
        x(131) = 0.036879149474284021657652D+00
        x(132) = 0.049163525671349973093019D+00
        x(133) = 0.061440469016428270850728D+00
        x(134) = 0.073708123403767780288977D+00
        x(135) = 0.085964634131980604256000D+00
        x(136) = 0.098208148184447540736015D+00
        x(137) = 0.110436814509468826540991D+00
        x(138) = 0.122648784300117812092492D+00
        x(139) = 0.134842211273755257250625D+00
        x(140) = 0.147015251951161989456661D+00
        x(141) = 0.159166065935247723154292D+00
        x(142) = 0.171292816189293903533225D+00
        x(143) = 0.183393669314688508087976D+00
        x(144) = 0.195466795828110816293869D+00
        x(145) = 0.207510370438124240859625D+00
        x(146) = 0.219522572321135403508985D+00
        x(147) = 0.231501585396677734059116D+00
        x(148) = 0.243445598601977973686482D+00
        x(149) = 0.255352806165764071686080D+00
        x(150) = 0.267221407881273079721012D+00
        x(151) = 0.279049609378417768508970D+00
        x(152) = 0.290835622395070819082047D+00
        x(153) = 0.302577665047425574167140D+00
        x(154) = 0.314273962099392474845918D+00
        x(155) = 0.325922745230990453444769D+00
        x(156) = 0.337522253305692705554261D+00
        x(157) = 0.349070732636686422161576D+00
        x(158) = 0.360566437252006227074021D+00
        x(159) = 0.372007629158501235092510D+00
        x(160) = 0.383392578604595822734854D+00
        x(161) = 0.394719564341804385683361D+00
        x(162) = 0.405986873884960545511889D+00
        x(163) = 0.417192803771121462605751D+00
        x(164) = 0.428335659817108112494341D+00
        x(165) = 0.439413757375642589040685D+00
        x(166) = 0.450425421590043710043279D+00
        x(167) = 0.461368987647442418771401D+00
        x(168) = 0.472242801030478698742627D+00
        x(169) = 0.483045217767441948626854D+00
        x(170) = 0.493774604680816999489812D+00
        x(171) = 0.504429339634198197635551D+00
        x(172) = 0.515007811777534223035005D+00
        x(173) = 0.525508421790666565699453D+00
        x(174) = 0.535929582125124840335150D+00
        x(175) = 0.546269717244142383159817D+00
        x(176) = 0.556527263860855843833077D+00
        x(177) = 0.566700671174652760010815D+00
        x(178) = 0.576788401105631382036211D+00
        x(179) = 0.586788928527137300685706D+00
        x(180) = 0.596700741496341721653202D+00
        x(181) = 0.606522341482826526536576D+00
        x(182) = 0.616252243595141561442344D+00
        x(183) = 0.625888976805299900901619D+00
        x(184) = 0.635431084171177146547142D+00
        x(185) = 0.644877123056781136890077D+00
        x(186) = 0.654225665350358766508700D+00
        x(187) = 0.663475297680306939970658D+00
        x(188) = 0.672624621628855017806731D+00
        x(189) = 0.681672253943486448787259D+00
        x(190) = 0.690616826746067624571761D+00
        x(191) = 0.699456987739652339456557D+00
        x(192) = 0.708191400412930589382399D+00
        x(193) = 0.716818744242290800531501D+00
        x(194) = 0.725337714891464938687812D+00
        x(195) = 0.733747024408726316001889D+00
        x(196) = 0.742045401421610281838045D+00
        x(197) = 0.750231591329128358931528D+00
        x(198) = 0.758304356491446765092016D+00
        x(199) = 0.766262476417000644100858D+00
        x(200) = 0.774104747947015717207115D+00
        x(201) = 0.781829985437409458675147D+00
        x(202) = 0.789437020938044295117764D+00
        x(203) = 0.796924704369305728807154D+00
        x(204) = 0.804291903695978689734633D+00
        x(205) = 0.811537505098395829833580D+00
        x(206) = 0.818660413140831885432406D+00
        x(207) = 0.825659550937118650611534D+00
        x(208) = 0.832533860313455524647230D+00
        x(209) = 0.839282301968391021084600D+00
        x(210) = 0.845903855629951054143931D+00
        x(211) = 0.852397520209890250084237D+00
        x(212) = 0.858762313955042966785823D+00
        x(213) = 0.864997274595751144137121D+00
        x(214) = 0.871101459491346550796200D+00
        x(215) = 0.877073945772665439532627D+00
        x(216) = 0.882913830481574073645470D+00
        x(217) = 0.888620230707484040924981D+00
        x(218) = 0.894192283720836729335637D+00
        x(219) = 0.899629147103536800144342D+00
        x(220) = 0.904929998876314959753358D+00
        x(221) = 0.910094037623000801254172D+00
        x(222) = 0.915120482611686961035103D+00
        x(223) = 0.920008573912766315142721D+00
        x(224) = 0.924757572513824425220425D+00
        x(225) = 0.929366760431369935739045D+00
        x(226) = 0.933835440819386124349338D+00
        x(227) = 0.938162938074687317626793D+00
        x(228) = 0.942348597939064408301480D+00
        x(229) = 0.946391787598204251752103D+00
        x(230) = 0.950291895777368285733522D+00
        x(231) = 0.954048332833816317950921D+00
        x(232) = 0.957660530845962076295490D+00
        x(233) = 0.961127943699247839572910D+00
        x(234) = 0.964450047168726298761719D+00
        x(235) = 0.967626338998338798105523D+00
        x(236) = 0.970656338976880365477697D+00
        x(237) = 0.973539589010643617645393D+00
        x(238) = 0.976275653192735980815246D+00
        x(239) = 0.978864117869068155239121D+00
        x(240) = 0.981304591701017185126565D+00
        x(241) = 0.983596705724776358640192D+00
        x(242) = 0.985740113407419277752900D+00
        x(243) = 0.987734490699732356281248D+00
        x(244) = 0.989579536085920123498574D+00
        x(245) = 0.991274970630385567164523D+00
        x(246) = 0.992820538021989138984811D+00
        x(247) = 0.994216004616630164799381D+00
        x(248) = 0.995461159480026294089975D+00
        x(249) = 0.996555814435198617028738D+00
        x(250) = 0.997499804126615814044844D+00
        x(251) = 0.998292986136967889228248D+00
        x(252) = 0.998935241284654635142155D+00
        x(253) = 0.999426474680169959344386D+00
        x(254) = 0.999766621312000569367063D+00
        x(255) = 0.999955705317563751730191D+00

        w(1) = 0.00011367361999142272115645954414D+00
        w(2) = 0.00026459387119083065532790838855D+00
        w(3) = 0.00041569762526823913616284210066D+00
        w(4) = 0.00056675794564824918946626058353D+00
        w(5) = 0.00071773647800611087798371518325D+00
        w(6) = 0.00086860766611945667949717690640D+00
        w(7) = 0.00101934797642732530281229369360D+00
        w(8) = 0.0011699343729388079886897709773D+00
        w(9) = 0.0013203439900221692090523602144D+00
        w(10) = 0.0014705540427783843160097204304D+00
        w(11) = 0.0016205417990415653896921100325D+00
        w(12) = 0.0017702845706603213070421243905D+00
        w(13) = 0.0019197597117132050055085980675D+00
        w(14) = 0.0020689446195015801533643667413D+00
        w(15) = 0.0022178167367540171700373764020D+00
        w(16) = 0.0023663535543962867157201855305D+00
        w(17) = 0.0025145326145997073931298921370D+00
        w(18) = 0.0026623315139717112732749157331D+00
        w(19) = 0.0028097279068204407457332299361D+00
        w(20) = 0.0029566995084575002760043344138D+00
        w(21) = 0.0031032240985191112621977893133D+00
        w(22) = 0.0032492795242943133198690930777D+00
        w(23) = 0.0033948437040533928255056951665D+00
        w(24) = 0.0035398946303722552150296713510D+00
        w(25) = 0.0036844103734499176530742235517D+00
        w(26) = 0.0038283690844171626400743524999D+00
        w(27) = 0.0039717489986349171988699773906D+00
        w(28) = 0.0041145284389812475901826468094D+00
        w(29) = 0.0042566858191260658425395494472D+00
        w(30) = 0.0043981996467927779838546384780D+00
        w(31) = 0.0045390485270061921259394035112D+00
        w(32) = 0.0046792111653260640506279893190D+00
        w(33) = 0.0048186663710656988918572043815D+00
        w(34) = 0.0049573930604950563104281084148D+00
        w(35) = 0.0050953702600278273039420404117D+00
        w(36) = 0.0052325771093919661294970523234D+00
        w(37) = 0.0053689928647831724787741258653D+00
        w(38) = 0.0055045969020008281904902120813D+00
        w(39) = 0.0056393687195659001929970994675D+00
        w(40) = 0.0057732879418203275712033691864D+00
        w(41) = 0.0059063343220074160130475409466D+00
        w(42) = 0.0060384877453327676663371666884D+00
        w(43) = 0.0061697282320052788060812561217D+00
        w(44) = 0.0063000359402577418025981070425D+00
        w(45) = 0.0064293911693465917826140832500D+00
        w(46) = 0.0065577743625303421548456356354D+00
        w(47) = 0.0066851661100262568757892743568D+00
        w(48) = 0.0068115471519448109954345674817D+00
        w(49) = 0.0069368983812014946719507501243D+00
        w(50) = 0.0070612008464055194979848418291D+00
        w(51) = 0.0071844357547249896530757997058D+00
        w(52) = 0.0073065844747281040972736443146D+00
        w(53) = 0.0074276285391999597581348419714D+00
        w(54) = 0.0075475496479345294426435656724D+00
        w(55) = 0.0076663296705013920315933272426D+00
        w(56) = 0.0077839506489867963897419914623D+00
        w(57) = 0.0079003948007086443529587296692D+00
        w(58) = 0.0080156445209049821352946484008D+00
        w(59) = 0.0081296823853955935356080649925D+00
        w(60) = 0.0082424911532162924158504385939D+00
        w(61) = 0.0083540537692255160718568405530D+00
        w(62) = 0.0084643533666828253227353760036D+00
        w(63) = 0.0085733732697989214067758505840D+00
        w(64) = 0.0086810969962567940901133439612D+00
        w(65) = 0.0087875082597036197689825483144D+00
        w(66) = 0.0088925909722130327769834298578D+00
        w(67) = 0.0089963292467173975949700110383D+00
        w(68) = 0.0090987073994097142025303711406D+00
        w(69) = 0.0091997099521147934060534414075D+00
        w(70) = 0.0092993216346293436285393234867D+00
        w(71) = 0.0093975273870306153500305317074D+00
        w(72) = 0.0094943123619532541442165010292D+00
        w(73) = 0.0095896619268340180657610209655D+00
        w(74) = 0.0096835616661240200035669970076D+00
        w(75) = 0.0097759973834681605268499842249D+00
        w(76) = 0.0098669551038514217128483481814D+00
        w(77) = 0.0099564210757116974565448593910D+00
        w(78) = 0.0100443817730188408231888789497D+00
        w(79) = 0.0101308238973196141129538950955D+00
        w(80) = 0.0102157343797482324629939488415D+00
        w(81) = 0.0102991003830021970147153502911D+00
        w(82) = 0.0103809093032831189224876935085D+00
        w(83) = 0.0104611487722022407735015844669D+00
        w(84) = 0.0105398066586503673262517188088D+00
        w(85) = 0.0106168710706319228563864391054D+00
        w(86) = 0.0106923303570628578226139809571D+00
        w(87) = 0.0107661731095321330311788312990D+00
        w(88) = 0.0108383881640265149842990798832D+00
        w(89) = 0.0109089646026184216450603134401D+00
        w(90) = 0.0109778917551165634377595759712D+00
        w(91) = 0.0110451592006791299277436662993D+00
        w(92) = 0.0111107567693892782875426356195D+00
        w(93) = 0.0111746745437926853557086684962D+00
        w(94) = 0.0112369028603969308303734810332D+00
        w(95) = 0.0112974323111324849102690558722D+00
        w(96) = 0.0113562537447750795009464486204D+00
        w(97) = 0.011413358268329247942299599697D+00
        w(98) = 0.011468737248372824084374355981D+00
        w(99) = 0.011522382312362197440930930031D+00
        w(100) = 0.011574285349898127083439539046D+00
        w(101) = 0.011624438513951922901227922331D+00
        w(102) = 0.011672834222051808845465154244D+00
        w(103) = 0.011719465157429288794653489478D+00
        w(104) = 0.011764324270125341726399410909D+00
        w(105) = 0.011807404778056278953532930501D+00
        w(106) = 0.011848700168039102281222824051D+00
        w(107) = 0.011888204196776208064673282076D+00
        w(108) = 0.011925910891799288293359117699D+00
        w(109) = 0.011961814552372285996633285380D+00
        w(110) = 0.011995909750353268455989686823D+00
        w(111) = 0.012028191331015087920350431142D+00
        w(112) = 0.012058654413824705751531083631D+00
        w(113) = 0.012087294393181062176578184854D+00
        w(114) = 0.012114106939111380091025793650D+00
        w(115) = 0.012139087997925797641334635250D+00
        w(116) = 0.012162233792830230614908682534D+00
        w(117) = 0.012183540824497371981177306326D+00
        w(118) = 0.012203005871595742256331865516D+00
        w(119) = 0.012220625991276710706457005806D+00
        w(120) = 0.012236398519619413758040249691D+00
        w(121) = 0.012250321072033503350218104906D+00
        w(122) = 0.012262391543619664338660618398D+00
        w(123) = 0.012272608109487846445745237751D+00
        w(124) = 0.012280969225033162644659793962D+00
        w(125) = 0.012287473626169412265336919908D+00
        w(126) = 0.012292120329520193516690694701D+00
        w(127) = 0.012294908632567576531532225710D+00
        w(128) = 0.01229583811375831445681490730D+00
        w(129) = 0.012294908632567576531532225710D+00
        w(130) = 0.012292120329520193516690694701D+00
        w(131) = 0.012287473626169412265336919908D+00
        w(132) = 0.012280969225033162644659793962D+00
        w(133) = 0.012272608109487846445745237751D+00
        w(134) = 0.012262391543619664338660618398D+00
        w(135) = 0.012250321072033503350218104906D+00
        w(136) = 0.012236398519619413758040249691D+00
        w(137) = 0.012220625991276710706457005806D+00
        w(138) = 0.012203005871595742256331865516D+00
        w(139) = 0.012183540824497371981177306326D+00
        w(140) = 0.012162233792830230614908682534D+00
        w(141) = 0.012139087997925797641334635250D+00
        w(142) = 0.012114106939111380091025793650D+00
        w(143) = 0.012087294393181062176578184854D+00
        w(144) = 0.012058654413824705751531083631D+00
        w(145) = 0.012028191331015087920350431142D+00
        w(146) = 0.011995909750353268455989686823D+00
        w(147) = 0.011961814552372285996633285380D+00
        w(148) = 0.011925910891799288293359117699D+00
        w(149) = 0.011888204196776208064673282076D+00
        w(150) = 0.011848700168039102281222824051D+00
        w(151) = 0.011807404778056278953532930501D+00
        w(152) = 0.011764324270125341726399410909D+00
        w(153) = 0.011719465157429288794653489478D+00
        w(154) = 0.011672834222051808845465154244D+00
        w(155) = 0.011624438513951922901227922331D+00
        w(156) = 0.011574285349898127083439539046D+00
        w(157) = 0.011522382312362197440930930031D+00
        w(158) = 0.011468737248372824084374355981D+00
        w(159) = 0.011413358268329247942299599697D+00
        w(160) = 0.0113562537447750795009464486204D+00
        w(161) = 0.0112974323111324849102690558722D+00
        w(162) = 0.0112369028603969308303734810332D+00
        w(163) = 0.0111746745437926853557086684962D+00
        w(164) = 0.0111107567693892782875426356195D+00
        w(165) = 0.0110451592006791299277436662993D+00
        w(166) = 0.0109778917551165634377595759712D+00
        w(167) = 0.0109089646026184216450603134401D+00
        w(168) = 0.0108383881640265149842990798832D+00
        w(169) = 0.0107661731095321330311788312990D+00
        w(170) = 0.0106923303570628578226139809571D+00
        w(171) = 0.0106168710706319228563864391054D+00
        w(172) = 0.0105398066586503673262517188088D+00
        w(173) = 0.0104611487722022407735015844669D+00
        w(174) = 0.0103809093032831189224876935085D+00
        w(175) = 0.0102991003830021970147153502911D+00
        w(176) = 0.0102157343797482324629939488415D+00
        w(177) = 0.0101308238973196141129538950955D+00
        w(178) = 0.0100443817730188408231888789497D+00
        w(179) = 0.0099564210757116974565448593910D+00
        w(180) = 0.0098669551038514217128483481814D+00
        w(181) = 0.0097759973834681605268499842249D+00
        w(182) = 0.0096835616661240200035669970076D+00
        w(183) = 0.0095896619268340180657610209655D+00
        w(184) = 0.0094943123619532541442165010292D+00
        w(185) = 0.0093975273870306153500305317074D+00
        w(186) = 0.0092993216346293436285393234867D+00
        w(187) = 0.0091997099521147934060534414075D+00
        w(188) = 0.0090987073994097142025303711406D+00
        w(189) = 0.0089963292467173975949700110383D+00
        w(190) = 0.0088925909722130327769834298578D+00
        w(191) = 0.0087875082597036197689825483144D+00
        w(192) = 0.0086810969962567940901133439612D+00
        w(193) = 0.0085733732697989214067758505840D+00
        w(194) = 0.0084643533666828253227353760036D+00
        w(195) = 0.0083540537692255160718568405530D+00
        w(196) = 0.0082424911532162924158504385939D+00
        w(197) = 0.0081296823853955935356080649925D+00
        w(198) = 0.0080156445209049821352946484008D+00
        w(199) = 0.0079003948007086443529587296692D+00
        w(200) = 0.0077839506489867963897419914623D+00
        w(201) = 0.0076663296705013920315933272426D+00
        w(202) = 0.0075475496479345294426435656724D+00
        w(203) = 0.0074276285391999597581348419714D+00
        w(204) = 0.0073065844747281040972736443146D+00
        w(205) = 0.0071844357547249896530757997058D+00
        w(206) = 0.0070612008464055194979848418291D+00
        w(207) = 0.0069368983812014946719507501243D+00
        w(208) = 0.0068115471519448109954345674817D+00
        w(209) = 0.0066851661100262568757892743568D+00
        w(210) = 0.0065577743625303421548456356354D+00
        w(211) = 0.0064293911693465917826140832500D+00
        w(212) = 0.0063000359402577418025981070425D+00
        w(213) = 0.0061697282320052788060812561217D+00
        w(214) = 0.0060384877453327676663371666884D+00
        w(215) = 0.0059063343220074160130475409466D+00
        w(216) = 0.0057732879418203275712033691864D+00
        w(217) = 0.0056393687195659001929970994675D+00
        w(218) = 0.0055045969020008281904902120813D+00
        w(219) = 0.0053689928647831724787741258653D+00
        w(220) = 0.0052325771093919661294970523234D+00
        w(221) = 0.0050953702600278273039420404117D+00
        w(222) = 0.0049573930604950563104281084148D+00
        w(223) = 0.0048186663710656988918572043815D+00
        w(224) = 0.0046792111653260640506279893190D+00
        w(225) = 0.0045390485270061921259394035112D+00
        w(226) = 0.0043981996467927779838546384780D+00
        w(227) = 0.0042566858191260658425395494472D+00
        w(228) = 0.0041145284389812475901826468094D+00
        w(229) = 0.0039717489986349171988699773906D+00
        w(230) = 0.0038283690844171626400743524999D+00
        w(231) = 0.0036844103734499176530742235517D+00
        w(232) = 0.0035398946303722552150296713510D+00
        w(233) = 0.0033948437040533928255056951665D+00
        w(234) = 0.0032492795242943133198690930777D+00
        w(235) = 0.0031032240985191112621977893133D+00
        w(236) = 0.0029566995084575002760043344138D+00
        w(237) = 0.0028097279068204407457332299361D+00
        w(238) = 0.0026623315139717112732749157331D+00
        w(239) = 0.0025145326145997073931298921370D+00
        w(240) = 0.0023663535543962867157201855305D+00
        w(241) = 0.0022178167367540171700373764020D+00
        w(242) = 0.0020689446195015801533643667413D+00
        w(243) = 0.0019197597117132050055085980675D+00
        w(244) = 0.0017702845706603213070421243905D+00
        w(245) = 0.0016205417990415653896921100325D+00
        w(246) = 0.0014705540427783843160097204304D+00
        w(247) = 0.0013203439900221692090523602144D+00
        w(248) = 0.0011699343729388079886897709773D+00
        w(249) = 0.00101934797642732530281229369360D+00
        w(250) = 0.00086860766611945667949717690640D+00
        w(251) = 0.00071773647800611087798371518325D+00
        w(252) = 0.00056675794564824918946626058353D+00
        w(253) = 0.00041569762526823913616284210066D+00
        w(254) = 0.00026459387119083065532790838855D+00
        w(255) = 0.00011367361999142272115645954414D+00

      else if ( n .eq. 256 ) then

        x(1) = -0.999956050018992230734801D+00
        x(2) = -0.999768437409263186104879D+00
        x(3) = -0.999430937466261408240854D+00
        x(4) = -0.998943525843408856555026D+00
        x(5) = -0.998306266473006444055500D+00
        x(6) = -0.997519252756720827563409D+00
        x(7) = -0.996582602023381540430504D+00
        x(8) = -0.995496454481096356592647D+00
        x(9) = -0.994260972922409664962878D+00
        x(10) = -0.992876342608822117143534D+00
        x(11) = -0.991342771207583086922189D+00
        x(12) = -0.989660488745065218319244D+00
        x(13) = -0.987829747564860608916488D+00
        x(14) = -0.985850822286125956479245D+00
        x(15) = -0.983724009760315496166686D+00
        x(16) = -0.981449629025464405769303D+00
        x(17) = -0.979028021257622038824238D+00
        x(18) = -0.976459549719234155621011D+00
        x(19) = -0.973744599704370405266079D+00
        x(20) = -0.970883578480743029320923D+00
        x(21) = -0.967876915228489454909004D+00
        x(22) = -0.964725060975706430932612D+00
        x(23) = -0.961428488530732144006407D+00
        x(24) = -0.957987692411178129365790D+00
        x(25) = -0.954403188769716241764448D+00
        x(26) = -0.950675515316628276363852D+00
        x(27) = -0.946805231239127481372052D+00
        x(28) = -0.942792917117462443183076D+00
        x(29) = -0.938639174837814804981926D+00
        x(30) = -0.934344627502003094292477D+00
        x(31) = -0.929909919334005641180246D+00
        x(32) = -0.925335715583316202872730D+00
        x(33) = -0.920622702425146495505047D+00
        x(34) = -0.915771586857490384526670D+00
        x(35) = -0.910783096595065011890907D+00
        x(36) = -0.905657979960144647082682D+00
        x(37) = -0.900397005770303544771620D+00
        x(38) = -0.895000963223084577441223D+00
        x(39) = -0.889470661777610888828677D+00
        x(40) = -0.883806931033158284859826D+00
        x(41) = -0.878010620604706543986435D+00
        x(42) = -0.872082599995488289130046D+00
        x(43) = -0.866023758466554519297515D+00
        x(44) = -0.859835004903376350696173D+00
        x(45) = -0.853517267679502965073036D+00
        x(46) = -0.847071494517296207187072D+00
        x(47) = -0.840498652345762713895068D+00
        x(48) = -0.833799727155504894348444D+00
        x(49) = -0.826975723850812514289093D+00
        x(50) = -0.820027666098917067403478D+00
        x(51) = -0.812956596176431543136410D+00
        x(52) = -0.805763574812998623257389D+00
        x(53) = -0.798449681032170758782543D+00
        x(54) = -0.791016011989545994546707D+00
        x(55) = -0.783463682808183820750670D+00
        x(56) = -0.775793826411325739132053D+00
        x(57) = -0.768007593352445635975891D+00
        x(58) = -0.760106151642655454941907D+00
        x(59) = -0.752090686575492059587530D+00
        x(60) = -0.743962400549111568455683D+00
        x(61) = -0.735722512885917834620373D+00
        x(62) = -0.727372259649652126586894D+00
        x(63) = -0.718912893459971448372640D+00
        x(64) = -0.710345683304543313394566D+00
        x(65) = -0.701671914348685159406084D+00
        x(66) = -0.692892887742576960105342D+00
        x(67) = -0.684009920426075953124877D+00
        x(68) = -0.675024344931162763855919D+00
        x(69) = -0.665937509182048559906408D+00
        x(70) = -0.656750776292973221887500D+00
        x(71) = -0.647465524363724862617016D+00
        x(72) = -0.638083146272911368668689D+00
        x(73) = -0.628605049469014975432210D+00
        x(74) = -0.619032655759261219430968D+00
        x(75) = -0.609367401096333939522311D+00
        x(76) = -0.599610735362968321730388D+00
        x(77) = -0.589764122154454300785786D+00
        x(78) = -0.579829038559082944921832D+00
        x(79) = -0.569806974936568759057668D+00
        x(80) = -0.559699434694481145136907D+00
        x(81) = -0.549507934062718557042427D+00
        x(82) = -0.539234001866059181127936D+00
        x(83) = -0.528879179294822261951476D+00
        x(84) = -0.518445019673674476221662D+00
        x(85) = -0.507933088228616036231925D+00
        x(86) = -0.497344961852181477119512D+00
        x(87) = -0.486682228866890350103621D+00
        x(88) = -0.475946488786983306390738D+00
        x(89) = -0.465139352078479313645570D+00
        x(90) = -0.454262439917589998774455D+00
        x(91) = -0.443317383947527357216926D+00
        x(92) = -0.432305826033741309953441D+00
        x(93) = -0.421229418017623824976812D+00
        x(94) = -0.410089821468716550006434D+00
        x(95) = -0.398888707435459127713463D+00
        x(96) = -0.387627756194515583637985D+00
        x(97) = -0.376308656998716390283056D+00
        x(98) = -0.364933107823654018533465D+00
        x(99) = -0.353502815112969989537790D+00
        x(100) = -0.342019493522371636480730D+00
        x(101) = -0.330484865662416976229187D+00
        x(102) = -0.318900661840106275631683D+00
        x(103) = -0.307268619799319076258610D+00
        x(104) = -0.295590484460135614563787D+00
        x(105) = -0.283868007657081741799766D+00
        x(106) = -0.272102947876336609505245D+00
        x(107) = -0.260297069991942541978561D+00
        x(108) = -0.248452145001056666833243D+00
        x(109) = -0.236569949758284018477508D+00
        x(110) = -0.224652266709131967147878D+00
        x(111) = -0.212700883622625957937040D+00
        x(112) = -0.200717593323126670068001D+00
        x(113) = -0.188704193421388826461504D+00
        x(114) = -0.176662486044901997403722D+00
        x(115) = -0.164594277567553849829285D+00
        x(116) = -0.152501378338656395374607D+00
        x(117) = -0.140385602411375885913025D+00
        x(118) = -0.128248767270607094742050D+00
        x(119) = -0.116092693560332804940735D+00
        x(120) = -0.103919204810509403639197D+00
        x(121) = -0.091730127163519552031146D+00
        x(122) = -0.079527289100232965903227D+00
        x(123) = -0.067312521165716400242290D+00
        x(124) = -0.055087655694633984104561D+00
        x(125) = -0.042854526536379098381242D+00
        x(126) = -0.030614968779979029366279D+00
        x(127) = -0.018370818478813665117926D+00
        x(128) = -0.006123912375189529501170D+00
        x(129) = 0.006123912375189529501170D+00
        x(130) = 0.018370818478813665117926D+00
        x(131) = 0.030614968779979029366279D+00
        x(132) = 0.042854526536379098381242D+00
        x(133) = 0.055087655694633984104561D+00
        x(134) = 0.067312521165716400242290D+00
        x(135) = 0.079527289100232965903227D+00
        x(136) = 0.091730127163519552031146D+00
        x(137) = 0.103919204810509403639197D+00
        x(138) = 0.116092693560332804940735D+00
        x(139) = 0.128248767270607094742050D+00
        x(140) = 0.140385602411375885913025D+00
        x(141) = 0.152501378338656395374607D+00
        x(142) = 0.164594277567553849829285D+00
        x(143) = 0.176662486044901997403722D+00
        x(144) = 0.188704193421388826461504D+00
        x(145) = 0.200717593323126670068001D+00
        x(146) = 0.212700883622625957937040D+00
        x(147) = 0.224652266709131967147878D+00
        x(148) = 0.236569949758284018477508D+00
        x(149) = 0.248452145001056666833243D+00
        x(150) = 0.260297069991942541978561D+00
        x(151) = 0.272102947876336609505245D+00
        x(152) = 0.283868007657081741799766D+00
        x(153) = 0.295590484460135614563787D+00
        x(154) = 0.307268619799319076258610D+00
        x(155) = 0.318900661840106275631683D+00
        x(156) = 0.330484865662416976229187D+00
        x(157) = 0.342019493522371636480730D+00
        x(158) = 0.353502815112969989537790D+00
        x(159) = 0.364933107823654018533465D+00
        x(160) = 0.376308656998716390283056D+00
        x(161) = 0.387627756194515583637985D+00
        x(162) = 0.398888707435459127713463D+00
        x(163) = 0.410089821468716550006434D+00
        x(164) = 0.421229418017623824976812D+00
        x(165) = 0.432305826033741309953441D+00
        x(166) = 0.443317383947527357216926D+00
        x(167) = 0.454262439917589998774455D+00
        x(168) = 0.465139352078479313645570D+00
        x(169) = 0.475946488786983306390738D+00
        x(170) = 0.486682228866890350103621D+00
        x(171) = 0.497344961852181477119512D+00
        x(172) = 0.507933088228616036231925D+00
        x(173) = 0.518445019673674476221662D+00
        x(174) = 0.528879179294822261951476D+00
        x(175) = 0.539234001866059181127936D+00
        x(176) = 0.549507934062718557042427D+00
        x(177) = 0.559699434694481145136907D+00
        x(178) = 0.569806974936568759057668D+00
        x(179) = 0.579829038559082944921832D+00
        x(180) = 0.589764122154454300785786D+00
        x(181) = 0.599610735362968321730388D+00
        x(182) = 0.609367401096333939522311D+00
        x(183) = 0.619032655759261219430968D+00
        x(184) = 0.628605049469014975432210D+00
        x(185) = 0.638083146272911368668689D+00
        x(186) = 0.647465524363724862617016D+00
        x(187) = 0.656750776292973221887500D+00
        x(188) = 0.665937509182048559906408D+00
        x(189) = 0.675024344931162763855919D+00
        x(190) = 0.684009920426075953124877D+00
        x(191) = 0.692892887742576960105342D+00
        x(192) = 0.701671914348685159406084D+00
        x(193) = 0.710345683304543313394566D+00
        x(194) = 0.718912893459971448372640D+00
        x(195) = 0.727372259649652126586894D+00
        x(196) = 0.735722512885917834620373D+00
        x(197) = 0.743962400549111568455683D+00
        x(198) = 0.752090686575492059587530D+00
        x(199) = 0.760106151642655454941907D+00
        x(200) = 0.768007593352445635975891D+00
        x(201) = 0.775793826411325739132053D+00
        x(202) = 0.783463682808183820750670D+00
        x(203) = 0.791016011989545994546707D+00
        x(204) = 0.798449681032170758782543D+00
        x(205) = 0.805763574812998623257389D+00
        x(206) = 0.812956596176431543136410D+00
        x(207) = 0.820027666098917067403478D+00
        x(208) = 0.826975723850812514289093D+00
        x(209) = 0.833799727155504894348444D+00
        x(210) = 0.840498652345762713895068D+00
        x(211) = 0.847071494517296207187072D+00
        x(212) = 0.853517267679502965073036D+00
        x(213) = 0.859835004903376350696173D+00
        x(214) = 0.866023758466554519297515D+00
        x(215) = 0.872082599995488289130046D+00
        x(216) = 0.878010620604706543986435D+00
        x(217) = 0.883806931033158284859826D+00
        x(218) = 0.889470661777610888828677D+00
        x(219) = 0.895000963223084577441223D+00
        x(220) = 0.900397005770303544771620D+00
        x(221) = 0.905657979960144647082682D+00
        x(222) = 0.910783096595065011890907D+00
        x(223) = 0.915771586857490384526670D+00
        x(224) = 0.920622702425146495505047D+00
        x(225) = 0.925335715583316202872730D+00
        x(226) = 0.929909919334005641180246D+00
        x(227) = 0.934344627502003094292477D+00
        x(228) = 0.938639174837814804981926D+00
        x(229) = 0.942792917117462443183076D+00
        x(230) = 0.946805231239127481372052D+00
        x(231) = 0.950675515316628276363852D+00
        x(232) = 0.954403188769716241764448D+00
        x(233) = 0.957987692411178129365790D+00
        x(234) = 0.961428488530732144006407D+00
        x(235) = 0.964725060975706430932612D+00
        x(236) = 0.967876915228489454909004D+00
        x(237) = 0.970883578480743029320923D+00
        x(238) = 0.973744599704370405266079D+00
        x(239) = 0.976459549719234155621011D+00
        x(240) = 0.979028021257622038824238D+00
        x(241) = 0.981449629025464405769303D+00
        x(242) = 0.983724009760315496166686D+00
        x(243) = 0.985850822286125956479245D+00
        x(244) = 0.987829747564860608916488D+00
        x(245) = 0.989660488745065218319244D+00
        x(246) = 0.991342771207583086922189D+00
        x(247) = 0.992876342608822117143534D+00
        x(248) = 0.994260972922409664962878D+00
        x(249) = 0.995496454481096356592647D+00
        x(250) = 0.996582602023381540430504D+00
        x(251) = 0.997519252756720827563409D+00
        x(252) = 0.998306266473006444055500D+00
        x(253) = 0.998943525843408856555026D+00
        x(254) = 0.999430937466261408240854D+00
        x(255) = 0.999768437409263186104879D+00
        x(256) = 0.999956050018992230734801D+00

        w(1) = 0.00011278901782227217551253887725D+00
        w(2) = 0.00026253494429644590628745756250D+00
        w(3) = 0.00041246325442617632843218583774D+00
        w(4) = 0.00056234895403140980281523674759D+00
        w(5) = 0.0007121541634733206669089891511D+00
        w(6) = 0.0008618537014200890378140934163D+00
        w(7) = 0.0010114243932084404526058128414D+00
        w(8) = 0.0011608435575677247239705981135D+00
        w(9) = 0.0013100886819025044578316804271D+00
        w(10) = 0.0014591373333107332010883864996D+00
        w(11) = 0.0016079671307493272424499395690D+00
        w(12) = 0.0017565557363307299936069145295D+00
        w(13) = 0.0019048808534997184044191411746D+00
        w(14) = 0.0020529202279661431745487818492D+00
        w(15) = 0.0022006516498399104996848834189D+00
        w(16) = 0.0023480529563273120170064609087D+00
        w(17) = 0.0024951020347037068508395354372D+00
        w(18) = 0.0026417768254274905641208292516D+00
        w(19) = 0.0027880553253277068805747610763D+00
        w(20) = 0.0029339155908297166460123254142D+00
        w(21) = 0.0030793357411993375832053528316D+00
        w(22) = 0.0032242939617941981570107134269D+00
        w(23) = 0.0033687685073155510120191062489D+00
        w(24) = 0.0035127377050563073309710549844D+00
        w(25) = 0.0036561799581425021693892413052D+00
        w(26) = 0.0037990737487662579981170192082D+00
        w(27) = 0.0039413976414088336277290349840D+00
        w(28) = 0.0040831302860526684085997759212D+00
        w(29) = 0.0042242504213815362723565049060D+00
        w(30) = 0.0043647368779680566815684200621D+00
        w(31) = 0.0045045685814478970686417923159D+00
        w(32) = 0.0046437245556800603139790923525D+00
        w(33) = 0.0047821839258926913729317340448D+00
        w(34) = 0.0049199259218138656695587765655D+00
        w(35) = 0.0050569298807868423875578160762D+00
        w(36) = 0.0051931752508692809303287536296D+00
        w(37) = 0.0053286415939159303170811114788D+00
        w(38) = 0.0054633085886443102775705318566D+00
        w(39) = 0.0055971560336829100775514452572D+00
        w(40) = 0.005730163850601437177384417555D+00
        w(41) = 0.005862312086922653060661598801D+00
        w(42) = 0.005993580919115338221127696870D+00
        w(43) = 0.006123950655567932542389081187D+00
        w(44) = 0.006253401739542401272063645975D+00
        w(45) = 0.006381914752107880570375164275D+00
        w(46) = 0.006509470415053660267809899951D+00
        w(47) = 0.006636049593781065044590038355D+00
        w(48) = 0.006761633300173798780927861108D+00
        w(49) = 0.006886202695446320346713323775D+00
        w(50) = 0.007009739092969822621234436194D+00
        w(51) = 0.007132223961075390071672422986D+00
        w(52) = 0.007253638925833913783829137214D+00
        w(53) = 0.007373965773812346437572440695D+00
        w(54) = 0.007493186454805883358599761133D+00
        w(55) = 0.007611283084545659461618719618D+00
        w(56) = 0.007728237947381555631110194958D+00
        w(57) = 0.007844033498939711866810316151D+00
        w(58) = 0.007958652368754348353613161227D+00
        w(59) = 0.008072077362873499500946974804D+00
        w(60) = 0.008184291466438269935619761004D+00
        w(61) = 0.008295277846235225425171412553D+00
        w(62) = 0.008405019853221535756180301698D+00
        w(63) = 0.008513501025022490693838354790D+00
        w(64) = 0.008620705088401014305368838410D+00
        w(65) = 0.008726615961698807140336632217D+00
        w(66) = 0.008831217757248750025318272685D+00
        w(67) = 0.008934494783758207548408417085D+00
        w(68) = 0.009036431548662873680227775572D+00
        w(69) = 0.009137012760450806402000472219D+00
        w(70) = 0.009236223330956302687378716714D+00
        w(71) = 0.009334048377623269712466014486D+00
        w(72) = 0.009430473225737752747352764482D+00
        w(73) = 0.009525483410629284811829685754D+00
        w(74) = 0.009619064679840727857162164401D+00
        w(75) = 0.009711202995266279964249670496D+00
        w(76) = 0.009801884535257327825498800250D+00
        w(77) = 0.009891095696695828602630683809D+00
        w(78) = 0.009978823097034910124733949495D+00
        w(79) = 0.010065053576306383309460978930D+00
        w(80) = 0.010149774199094865654634066042D+00
        w(81) = 0.010232972256478219656954857160D+00
        w(82) = 0.010314635267934015068260713997D+00
        w(83) = 0.010394750983211728997101725205D+00
        w(84) = 0.010473307384170403003569566927D+00
        w(85) = 0.010550292686581481517533575536D+00
        w(86) = 0.010625695341896561133961681801D+00
        w(87) = 0.010699504038979785603048200583D+00
        w(88) = 0.010771707705804626636653631927D+00
        w(89) = 0.010842295511114795995293477058D+00
        w(90) = 0.010911256866049039700796847788D+00
        w(91) = 0.010978581425729570637988203448D+00
        w(92) = 0.011044259090813901263517571044D+00
        w(93) = 0.011108280009009843630460815451D+00
        w(94) = 0.011170634576553449462710881938D+00
        w(95) = 0.011231313439649668572656802083D+00
        w(96) = 0.011290307495875509508367594121D+00
        w(97) = 0.011347607895545491941625714297D+00
        w(98) = 0.011403206043039185964847059552D+00
        w(99) = 0.011457093598090639152334392298D+00
        w(100) = 0.011509262477039497958586392439D+00
        w(101) = 0.011559704854043635772668656950D+00
        w(102) = 0.011608413162253105722084706677D+00
        w(103) = 0.011655380094945242121298939730D+00
        w(104) = 0.011700598606620740288189823359D+00
        w(105) = 0.011744061914060550305376732759D+00
        w(106) = 0.011785763497343426181690117627D+00
        w(107) = 0.011825697100823977771160737958D+00
        w(108) = 0.011863856734071078731904572908D+00
        w(109) = 0.011900236672766489754287204237D+00
        w(110) = 0.011934831459563562255873201696D+00
        w(111) = 0.011967635904905893729007282670D+00
        w(112) = 0.011998645087805811934536710071D+00
        w(113) = 0.012027854356582571161267533498D+00
        w(114) = 0.012055259329560149814347085327D+00
        w(115) = 0.012080855895724544655975183976D+00
        w(116) = 0.012104640215340463097757829736D+00
        w(117) = 0.012126608720527321034718492205D+00
        w(118) = 0.012146758115794459815559837664D+00
        w(119) = 0.012165085378535502061307291839D+00
        w(120) = 0.012181587759481772174047585032D+00
        w(121) = 0.012196262783114713518180974196D+00
        w(122) = 0.012209108248037240407514094371D+00
        w(123) = 0.012220122227303969191708737227D+00
        w(124) = 0.012229303068710278904146266083D+00
        w(125) = 0.012236649395040158109242574767D+00
        w(126) = 0.012242160104272800769728083260D+00
        w(127) = 0.012245834369747920142463857550D+00
        w(128) = 0.01224767164028975590407032649D+00
        w(129) = 0.01224767164028975590407032649D+00
        w(130) = 0.012245834369747920142463857550D+00
        w(131) = 0.012242160104272800769728083260D+00
        w(132) = 0.012236649395040158109242574767D+00
        w(133) = 0.012229303068710278904146266083D+00
        w(134) = 0.012220122227303969191708737227D+00
        w(135) = 0.012209108248037240407514094371D+00
        w(136) = 0.012196262783114713518180974196D+00
        w(137) = 0.012181587759481772174047585032D+00
        w(138) = 0.012165085378535502061307291839D+00
        w(139) = 0.012146758115794459815559837664D+00
        w(140) = 0.012126608720527321034718492205D+00
        w(141) = 0.012104640215340463097757829736D+00
        w(142) = 0.012080855895724544655975183976D+00
        w(143) = 0.012055259329560149814347085327D+00
        w(144) = 0.012027854356582571161267533498D+00
        w(145) = 0.011998645087805811934536710071D+00
        w(146) = 0.011967635904905893729007282670D+00
        w(147) = 0.011934831459563562255873201696D+00
        w(148) = 0.011900236672766489754287204237D+00
        w(149) = 0.011863856734071078731904572908D+00
        w(150) = 0.011825697100823977771160737958D+00
        w(151) = 0.011785763497343426181690117627D+00
        w(152) = 0.011744061914060550305376732759D+00
        w(153) = 0.011700598606620740288189823359D+00
        w(154) = 0.011655380094945242121298939730D+00
        w(155) = 0.011608413162253105722084706677D+00
        w(156) = 0.011559704854043635772668656950D+00
        w(157) = 0.011509262477039497958586392439D+00
        w(158) = 0.011457093598090639152334392298D+00
        w(159) = 0.011403206043039185964847059552D+00
        w(160) = 0.011347607895545491941625714297D+00
        w(161) = 0.011290307495875509508367594121D+00
        w(162) = 0.011231313439649668572656802083D+00
        w(163) = 0.011170634576553449462710881938D+00
        w(164) = 0.011108280009009843630460815451D+00
        w(165) = 0.011044259090813901263517571044D+00
        w(166) = 0.010978581425729570637988203448D+00
        w(167) = 0.010911256866049039700796847788D+00
        w(168) = 0.010842295511114795995293477058D+00
        w(169) = 0.010771707705804626636653631927D+00
        w(170) = 0.010699504038979785603048200583D+00
        w(171) = 0.010625695341896561133961681801D+00
        w(172) = 0.010550292686581481517533575536D+00
        w(173) = 0.010473307384170403003569566927D+00
        w(174) = 0.010394750983211728997101725205D+00
        w(175) = 0.010314635267934015068260713997D+00
        w(176) = 0.010232972256478219656954857160D+00
        w(177) = 0.010149774199094865654634066042D+00
        w(178) = 0.010065053576306383309460978930D+00
        w(179) = 0.009978823097034910124733949495D+00
        w(180) = 0.009891095696695828602630683809D+00
        w(181) = 0.009801884535257327825498800250D+00
        w(182) = 0.009711202995266279964249670496D+00
        w(183) = 0.009619064679840727857162164401D+00
        w(184) = 0.009525483410629284811829685754D+00
        w(185) = 0.009430473225737752747352764482D+00
        w(186) = 0.009334048377623269712466014486D+00
        w(187) = 0.009236223330956302687378716714D+00
        w(188) = 0.009137012760450806402000472219D+00
        w(189) = 0.009036431548662873680227775572D+00
        w(190) = 0.008934494783758207548408417085D+00
        w(191) = 0.008831217757248750025318272685D+00
        w(192) = 0.008726615961698807140336632217D+00
        w(193) = 0.008620705088401014305368838410D+00
        w(194) = 0.008513501025022490693838354790D+00
        w(195) = 0.008405019853221535756180301698D+00
        w(196) = 0.008295277846235225425171412553D+00
        w(197) = 0.008184291466438269935619761004D+00
        w(198) = 0.008072077362873499500946974804D+00
        w(199) = 0.007958652368754348353613161227D+00
        w(200) = 0.007844033498939711866810316151D+00
        w(201) = 0.007728237947381555631110194958D+00
        w(202) = 0.007611283084545659461618719618D+00
        w(203) = 0.007493186454805883358599761133D+00
        w(204) = 0.007373965773812346437572440695D+00
        w(205) = 0.007253638925833913783829137214D+00
        w(206) = 0.007132223961075390071672422986D+00
        w(207) = 0.007009739092969822621234436194D+00
        w(208) = 0.006886202695446320346713323775D+00
        w(209) = 0.006761633300173798780927861108D+00
        w(210) = 0.006636049593781065044590038355D+00
        w(211) = 0.006509470415053660267809899951D+00
        w(212) = 0.006381914752107880570375164275D+00
        w(213) = 0.006253401739542401272063645975D+00
        w(214) = 0.006123950655567932542389081187D+00
        w(215) = 0.005993580919115338221127696870D+00
        w(216) = 0.005862312086922653060661598801D+00
        w(217) = 0.005730163850601437177384417555D+00
        w(218) = 0.0055971560336829100775514452572D+00
        w(219) = 0.0054633085886443102775705318566D+00
        w(220) = 0.0053286415939159303170811114788D+00
        w(221) = 0.0051931752508692809303287536296D+00
        w(222) = 0.0050569298807868423875578160762D+00
        w(223) = 0.0049199259218138656695587765655D+00
        w(224) = 0.0047821839258926913729317340448D+00
        w(225) = 0.0046437245556800603139790923525D+00
        w(226) = 0.0045045685814478970686417923159D+00
        w(227) = 0.0043647368779680566815684200621D+00
        w(228) = 0.0042242504213815362723565049060D+00
        w(229) = 0.0040831302860526684085997759212D+00
        w(230) = 0.0039413976414088336277290349840D+00
        w(231) = 0.0037990737487662579981170192082D+00
        w(232) = 0.0036561799581425021693892413052D+00
        w(233) = 0.0035127377050563073309710549844D+00
        w(234) = 0.0033687685073155510120191062489D+00
        w(235) = 0.0032242939617941981570107134269D+00
        w(236) = 0.0030793357411993375832053528316D+00
        w(237) = 0.0029339155908297166460123254142D+00
        w(238) = 0.0027880553253277068805747610763D+00
        w(239) = 0.0026417768254274905641208292516D+00
        w(240) = 0.0024951020347037068508395354372D+00
        w(241) = 0.0023480529563273120170064609087D+00
        w(242) = 0.0022006516498399104996848834189D+00
        w(243) = 0.0020529202279661431745487818492D+00
        w(244) = 0.0019048808534997184044191411746D+00
        w(245) = 0.0017565557363307299936069145295D+00
        w(246) = 0.0016079671307493272424499395690D+00
        w(247) = 0.0014591373333107332010883864996D+00
        w(248) = 0.0013100886819025044578316804271D+00
        w(249) = 0.0011608435575677247239705981135D+00
        w(250) = 0.0010114243932084404526058128414D+00
        w(251) = 0.0008618537014200890378140934163D+00
        w(252) = 0.0007121541634733206669089891511D+00
        w(253) = 0.00056234895403140980281523674759D+00
        w(254) = 0.00041246325442617632843218583774D+00
        w(255) = 0.00026253494429644590628745756250D+00
        w(256) = 0.00011278901782227217551253887725D+00

      else if ( n .eq. 257 ) then

        x(1) = -0.999956390712330402472857D+00
        x(2) = -0.999770232390338019056053D+00
        x(3) = -0.999435348366365078441838D+00
        x(4) = -0.998951714093223210129834D+00
        x(5) = -0.998319392445383847808766D+00
        x(6) = -0.997538475365520218731818D+00
        x(7) = -0.996609078365487004512326D+00
        x(8) = -0.995531339486830143483750D+00
        x(9) = -0.994305419008553630362377D+00
        x(10) = -0.992931499332908653172844D+00
        x(11) = -0.991409784923101705201254D+00
        x(12) = -0.989740502257507526030375D+00
        x(13) = -0.987923899788618253106809D+00
        x(14) = -0.985960247902290665366669D+00
        x(15) = -0.983849838875444644048531D+00
        x(16) = -0.981592986831381877693095D+00
        x(17) = -0.979190027692327124191591D+00
        x(18) = -0.976641319128992592610888D+00
        x(19) = -0.973947240507062326750976D+00
        x(20) = -0.971108192830542793021113D+00
        x(21) = -0.968124598681952354372943D+00
        x(22) = -0.964996902159337170373447D+00
        x(23) = -0.961725568810109767190665D+00
        x(24) = -0.958311085561711847074814D+00
        x(25) = -0.954753960649106318830855D+00
        x(26) = -0.951054723539105826691801D+00
        x(27) = -0.947213924851546682950881D+00
        x(28) = -0.943232136277318328151464D+00
        x(29) = -0.939109950493259404355123D+00
        x(30) = -0.934847981073932324370129D+00
        x(31) = -0.930446862400288909805510D+00
        x(32) = -0.925907249565240289235888D+00
        x(33) = -0.921229818276144817520964D+00
        x(34) = -0.916415264754228313295468D+00
        x(35) = -0.911464305630951423630955D+00
        x(36) = -0.906377677841339419411308D+00
        x(37) = -0.901156138514290206476301D+00
        x(38) = -0.895800464859876809085345D+00
        x(39) = -0.890311454053661045810287D+00
        x(40) = -0.884689923118035575018750D+00
        x(41) = -0.878936708800611938658765D+00
        x(42) = -0.873052667449672679799858D+00
        x(43) = -0.867038674886706051812473D+00
        x(44) = -0.860895626276042275514686D+00
        x(45) = -0.854624435991610735314055D+00
        x(46) = -0.848226037480837936478636D+00
        x(47) = -0.841701383125706473284556D+00
        x(48) = -0.835051444100995681967937D+00
        x(49) = -0.828277210229725073186687D+00
        x(50) = -0.821379689835822056081139D+00
        x(51) = -0.814359909594035880004229D+00
        x(52) = -0.807218914377120130552073D+00
        x(53) = -0.799957767100306523636066D+00
        x(54) = -0.792577548563093144962574D+00
        x(55) = -0.785079357288370682385816D+00
        x(56) = -0.777464309358910595129671D+00
        x(57) = -0.769733538251239556788216D+00
        x(58) = -0.761888194666924898264210D+00
        x(59) = -0.753929446361296162339238D+00
        x(60) = -0.745858477969628263337895D+00
        x(61) = -0.737676490830812123299244D+00
        x(62) = -0.729384702808539030149808D+00
        x(63) = -0.720984348110025333531072D+00
        x(64) = -0.712476677102304460118510D+00
        x(65) = -0.703862956126113592426171D+00
        x(66) = -0.695144467307402713168813D+00
        x(67) = -0.686322508366494071200553D+00
        x(68) = -0.677398392424920474813593D+00
        x(69) = -0.668373447809971163711735D+00
        x(70) = -0.659249017856974352220492D+00
        x(71) = -0.650026460709345873208532D+00
        x(72) = -0.640707149116433684724434D+00
        x(73) = -0.631292470229188329449219D+00
        x(74) = -0.621783825393689760680446D+00
        x(75) = -0.612182629942561267650033D+00
        x(76) = -0.602490312984301547488097D+00
        x(77) = -0.592708317190566281032495D+00
        x(78) = -0.582838098581430874902446D+00
        x(79) = -0.572881126308666332759406D+00
        x(80) = -0.562838882437060514424546D+00
        x(81) = -0.552712861723817332466074D+00
        x(82) = -0.542504571396066721967792D+00
        x(83) = -0.532215530926518500400434D+00
        x(84) = -0.521847271807293510797499D+00
        x(85) = -0.511401337321965712746629D+00
        x(86) = -0.500879282315849152005553D+00
        x(87) = -0.490282672964564000798817D+00
        x(88) = -0.479613086540916117008992D+00
        x(89) = -0.468872111180124821505728D+00
        x(90) = -0.458061345643433838720630D+00
        x(91) = -0.447182399080140586238810D+00
        x(92) = -0.436236890788079234603398D+00
        x(93) = -0.425226449972593188682213D+00
        x(94) = -0.414152715504032866791986D+00
        x(95) = -0.403017335673814873281489D+00
        x(96) = -0.391821967949078874408131D+00
        x(97) = -0.380568278725978696070941D+00
        x(98) = -0.369257943081644365255611D+00
        x(99) = -0.357892644524852014873858D+00
        x(100) = -0.346474074745438764010632D+00
        x(101) = -0.335003933362499872399782D+00
        x(102) = -0.323483927671405649204085D+00
        x(103) = -0.311915772389675771851948D+00
        x(104) = -0.300301189401748840754520D+00
        x(105) = -0.288641907502685160168097D+00
        x(106) = -0.276939662140840894253032D+00
        x(107) = -0.265196195159551900488370D+00
        x(108) = -0.253413254537865690008131D+00
        x(109) = -0.241592594130360106108882D+00
        x(110) = -0.229735973406087448117604D+00
        x(111) = -0.217845157186682897983880D+00
        x(112) = -0.205921915383676231351599D+00
        x(113) = -0.193968022735045913454182D+00
        x(114) = -0.181985258541054792946197D+00
        x(115) = -0.169975406399406713716337D+00
        x(116) = -0.157940253939763465806087D+00
        x(117) = -0.145881592557661591770148D+00
        x(118) = -0.133801217147868654144405D+00
        x(119) = -0.121700925837218653121859D+00
        x(120) = -0.109582519716966361063898D+00
        x(121) = -0.097447802574700412082119D+00
        x(122) = -0.085298580625855050603929D+00
        x(123) = -0.073136662244860502573600D+00
        x(124) = -0.060963857695971986730406D+00
        x(125) = -0.048781978863817431238958D+00
        x(126) = -0.036592838983704002816750D+00
        x(127) = -0.024398252371723591403953D+00
        x(128) = -0.012200034154697423345412D+00
        x(129) = 0.000000000000000000000000D+00
        x(130) = 0.012200034154697423345412D+00
        x(131) = 0.024398252371723591403953D+00
        x(132) = 0.036592838983704002816750D+00
        x(133) = 0.048781978863817431238958D+00
        x(134) = 0.060963857695971986730406D+00
        x(135) = 0.073136662244860502573600D+00
        x(136) = 0.085298580625855050603929D+00
        x(137) = 0.097447802574700412082119D+00
        x(138) = 0.109582519716966361063898D+00
        x(139) = 0.121700925837218653121859D+00
        x(140) = 0.133801217147868654144405D+00
        x(141) = 0.145881592557661591770148D+00
        x(142) = 0.157940253939763465806087D+00
        x(143) = 0.169975406399406713716337D+00
        x(144) = 0.181985258541054792946197D+00
        x(145) = 0.193968022735045913454182D+00
        x(146) = 0.205921915383676231351599D+00
        x(147) = 0.217845157186682897983880D+00
        x(148) = 0.229735973406087448117604D+00
        x(149) = 0.241592594130360106108882D+00
        x(150) = 0.253413254537865690008131D+00
        x(151) = 0.265196195159551900488370D+00
        x(152) = 0.276939662140840894253032D+00
        x(153) = 0.288641907502685160168097D+00
        x(154) = 0.300301189401748840754520D+00
        x(155) = 0.311915772389675771851948D+00
        x(156) = 0.323483927671405649204085D+00
        x(157) = 0.335003933362499872399782D+00
        x(158) = 0.346474074745438764010632D+00
        x(159) = 0.357892644524852014873858D+00
        x(160) = 0.369257943081644365255611D+00
        x(161) = 0.380568278725978696070941D+00
        x(162) = 0.391821967949078874408131D+00
        x(163) = 0.403017335673814873281489D+00
        x(164) = 0.414152715504032866791986D+00
        x(165) = 0.425226449972593188682213D+00
        x(166) = 0.436236890788079234603398D+00
        x(167) = 0.447182399080140586238810D+00
        x(168) = 0.458061345643433838720630D+00
        x(169) = 0.468872111180124821505728D+00
        x(170) = 0.479613086540916117008992D+00
        x(171) = 0.490282672964564000798817D+00
        x(172) = 0.500879282315849152005553D+00
        x(173) = 0.511401337321965712746629D+00
        x(174) = 0.521847271807293510797499D+00
        x(175) = 0.532215530926518500400434D+00
        x(176) = 0.542504571396066721967792D+00
        x(177) = 0.552712861723817332466074D+00
        x(178) = 0.562838882437060514424546D+00
        x(179) = 0.572881126308666332759406D+00
        x(180) = 0.582838098581430874902446D+00
        x(181) = 0.592708317190566281032495D+00
        x(182) = 0.602490312984301547488097D+00
        x(183) = 0.612182629942561267650033D+00
        x(184) = 0.621783825393689760680446D+00
        x(185) = 0.631292470229188329449219D+00
        x(186) = 0.640707149116433684724434D+00
        x(187) = 0.650026460709345873208532D+00
        x(188) = 0.659249017856974352220492D+00
        x(189) = 0.668373447809971163711735D+00
        x(190) = 0.677398392424920474813593D+00
        x(191) = 0.686322508366494071200553D+00
        x(192) = 0.695144467307402713168813D+00
        x(193) = 0.703862956126113592426171D+00
        x(194) = 0.712476677102304460118510D+00
        x(195) = 0.720984348110025333531072D+00
        x(196) = 0.729384702808539030149808D+00
        x(197) = 0.737676490830812123299244D+00
        x(198) = 0.745858477969628263337895D+00
        x(199) = 0.753929446361296162339238D+00
        x(200) = 0.761888194666924898264210D+00
        x(201) = 0.769733538251239556788216D+00
        x(202) = 0.777464309358910595129671D+00
        x(203) = 0.785079357288370682385816D+00
        x(204) = 0.792577548563093144962574D+00
        x(205) = 0.799957767100306523636066D+00
        x(206) = 0.807218914377120130552073D+00
        x(207) = 0.814359909594035880004229D+00
        x(208) = 0.821379689835822056081139D+00
        x(209) = 0.828277210229725073186687D+00
        x(210) = 0.835051444100995681967937D+00
        x(211) = 0.841701383125706473284556D+00
        x(212) = 0.848226037480837936478636D+00
        x(213) = 0.854624435991610735314055D+00
        x(214) = 0.860895626276042275514686D+00
        x(215) = 0.867038674886706051812473D+00
        x(216) = 0.873052667449672679799858D+00
        x(217) = 0.878936708800611938658765D+00
        x(218) = 0.884689923118035575018750D+00
        x(219) = 0.890311454053661045810287D+00
        x(220) = 0.895800464859876809085345D+00
        x(221) = 0.901156138514290206476301D+00
        x(222) = 0.906377677841339419411308D+00
        x(223) = 0.911464305630951423630955D+00
        x(224) = 0.916415264754228313295468D+00
        x(225) = 0.921229818276144817520964D+00
        x(226) = 0.925907249565240289235888D+00
        x(227) = 0.930446862400288909805510D+00
        x(228) = 0.934847981073932324370129D+00
        x(229) = 0.939109950493259404355123D+00
        x(230) = 0.943232136277318328151464D+00
        x(231) = 0.947213924851546682950881D+00
        x(232) = 0.951054723539105826691801D+00
        x(233) = 0.954753960649106318830855D+00
        x(234) = 0.958311085561711847074814D+00
        x(235) = 0.961725568810109767190665D+00
        x(236) = 0.964996902159337170373447D+00
        x(237) = 0.968124598681952354372943D+00
        x(238) = 0.971108192830542793021113D+00
        x(239) = 0.973947240507062326750976D+00
        x(240) = 0.976641319128992592610888D+00
        x(241) = 0.979190027692327124191591D+00
        x(242) = 0.981592986831381877693095D+00
        x(243) = 0.983849838875444644048531D+00
        x(244) = 0.985960247902290665366669D+00
        x(245) = 0.987923899788618253106809D+00
        x(246) = 0.989740502257507526030375D+00
        x(247) = 0.991409784923101705201254D+00
        x(248) = 0.992931499332908653172844D+00
        x(249) = 0.994305419008553630362377D+00
        x(250) = 0.995531339486830143483750D+00
        x(251) = 0.996609078365487004512326D+00
        x(252) = 0.997538475365520218731818D+00
        x(253) = 0.998319392445383847808766D+00
        x(254) = 0.998951714093223210129834D+00
        x(255) = 0.999435348366365078441838D+00
        x(256) = 0.999770232390338019056053D+00
        x(257) = 0.999956390712330402472857D+00

        w(1) = 0.00011191470145601756450862287886D+00
        w(2) = 0.00026049995580176964436806680831D+00
        w(3) = 0.00040926648283531339591138751432D+00
        w(4) = 0.00055799120546880640169677292533D+00
        w(5) = 0.00070663671051592291949335494247D+00
        w(6) = 0.00085517818446696565626595950963D+00
        w(7) = 0.00100359280467969441299468763292D+00
        w(8) = 0.0011518582377826677880963146741D+00
        w(9) = 0.0012999523174235227389668643832D+00
        w(10) = 0.0014478529559255120065233994722D+00
        w(11) = 0.0015955381166175133369701690235D+00
        w(12) = 0.0017429858051468299509941139300D+00
        w(13) = 0.0018901740676190104269878470891D+00
        w(14) = 0.0020370809914723626741694800322D+00
        w(15) = 0.0021836847075455253317921866057D+00
        w(16) = 0.0023299633927021828561308282641D+00
        w(17) = 0.0024758952727301488651840215879D+00
        w(18) = 0.0026214586253808109266552781372D+00
        w(19) = 0.0027666317834818283552560256501D+00
        w(20) = 0.0029113931380877846359302447381D+00
        w(21) = 0.0030557211416493711130936102459D+00
        w(22) = 0.0031995943111899437356540290142D+00
        w(23) = 0.0033429912314827618499065991316D+00
        w(24) = 0.0034858905582247143702551557840D+00
        w(25) = 0.0036282710212037760873102463983D+00
        w(26) = 0.0037701114274582873548537007645D+00
        w(27) = 0.0039113906644266662571543468015D+00
        w(28) = 0.0040520877030864825223229951262D+00
        w(29) = 0.0041921816010820254766367595011D+00
        w(30) = 0.0043316515058396297504806208252D+00
        w(31) = 0.0044704766576701092218388764046D+00
        w(32) = 0.0046086363928577081326523656522D+00
        w(33) = 0.0047461101467350184936945641585D+00
        w(34) = 0.0048828774567433411142588306018D+00
        w(35) = 0.0050189179654779878773297516544D+00
        w(36) = 0.0051542114237180378340642003713D+00
        w(37) = 0.0052887376934400710240953933529D+00
        w(38) = 0.0054224767508154127788846727083D+00
        w(39) = 0.0055554086891904284012033890901D+00
        w(40) = 0.0056875137220494140577838938236D+00
        w(41) = 0.0058187721859596348346566361185D+00
        w(42) = 0.0059491645434980654366600347567D+00
        w(43) = 0.0060786713861593931405204596709D+00
        w(44) = 0.0062072734372448464599330978665D+00
        w(45) = 0.0063349515547314166407936938524D+00
        w(46) = 0.0064616867341210426397202932350D+00
        w(47) = 0.0065874601112693336961737372300D+00
        w(48) = 0.0067122529651934070221351960200D+00
        w(49) = 0.0068360467208584215286561508406D+00
        w(50) = 0.0069588229519423919043121805236D+00
        w(51) = 0.0070805633835788707705149901066D+00
        w(52) = 0.0072012498950770900730828552207D+00
        w(53) = 0.0073208645226191563361371026044D+00
        w(54) = 0.0074393894619338979090297315972D+00
        w(55) = 0.0075568070709469658838993300454D+00
        w(56) = 0.0076730998724067939537782250476D+00
        w(57) = 0.0077882505564860261212726654404D+00
        w(58) = 0.0079022419833580248574070864277D+00
        w(59) = 0.0080150571857480760504667455353D+00
        w(60) = 0.0081266793714589108764118189068D+00
        w(61) = 0.0082370919258701685661946145361D+00
        w(62) = 0.0083462784144114279413811886655D+00
        w(63) = 0.0084542225850084395379670551258D+00
        w(64) = 0.0085609083705021941391459209280D+00
        w(65) = 0.0086663198910404675908861979240D+00
        w(66) = 0.0087704414564414858792445834744D+00
        w(67) = 0.0088732575685293586050755892934D+00
        w(68) = 0.0089747529234409331997949023068D+00
        w(69) = 0.0090749124139037264846862498962D+00
        w(70) = 0.0091737211314845944854270065178D+00
        w(71) = 0.0092711643688088057725325917169D+00
        w(72) = 0.0093672276217491880067391857021D+00
        w(73) = 0.0094618965915850218253881576301D+00
        w(74) = 0.0095551571871303607110514249099D+00
        w(75) = 0.0096469955268314600363329731559D+00
        w(76) = 0.0097373979408330030783691793250D+00
        w(77) = 0.0098263509730128164423854701706D+00
        w(78) = 0.0099138413829847720250916955489D+00
        w(79) = 0.0099998561480695773850435626986D+00
        w(80) = 0.0100843824652331611676814627839D+00
        w(81) = 0.0101674077529923650568895461852D+00
        w(82) = 0.0102489196532876585918958554047D+00
        w(83) = 0.0103289060333225980974485876288D+00
        w(84) = 0.0104073549873697559257355517893D+00
        w(85) = 0.0104842548385428511997370260353D+00
        w(86) = 0.0105595941405348182788823332058D+00
        w(87) = 0.0106333616793215542382761147904D+00
        w(88) = 0.0107055464748310917616231511294D+00
        w(89) = 0.0107761377825779489945556541150D+00
        w(90) = 0.0108451250952624130885928632830D+00
        w(91) = 0.0109124981443345193856719616965D+00
        w(92) = 0.0109782469015224934483083029166D+00
        w(93) = 0.0110423615803254284301924654946D+00
        w(94) = 0.0111048326374699756056269264803D+00
        w(95) = 0.0111656507743308312328559850485D+00
        w(96) = 0.0112248069383148083152535688671D+00
        w(97) = 0.0112822923242082872447042603128D+00
        w(98) = 0.0113380983754878447625379269120D+00
        w(99) = 0.011392216785593866154247619654D+00
        w(100) = 0.011444639499166951104119199270D+00
        w(101) = 0.011495358713246929174010288914D+00
        w(102) = 0.011544366878434306436012137033D+00
        w(103) = 0.011591656700013970380783131035D+00
        w(104) = 0.011637221139040985841125311445D+00
        w(105) = 0.011681053413388320313049670635D+00
        w(106) = 0.011723146998756342723302879656D+00
        w(107) = 0.011763495629643945382264331878D+00
        w(108) = 0.011802093300281144573421477037D+00
        w(109) = 0.011838934265523020964443424791D+00
        w(110) = 0.011874013041704866779344562066D+00
        w(111) = 0.011907324407458412445505183140D+00
        w(112) = 0.011938863404489011222535627643D+00
        w(113) = 0.011968625338313666131272065445D+00
        w(114) = 0.011996605778959789329711050159D+00
        w(115) = 0.012022800561624589927558893338D+00
        w(116) = 0.012047205787294992091420946532D+00
        w(117) = 0.012069817823327991167612855626D+00
        w(118) = 0.012090633303991361438266420912D+00
        w(119) = 0.012109649130964635027950450318D+00
        w(120) = 0.012126862473800277391553601370D+00
        w(121) = 0.012142270770344990738801546574D+00
        w(122) = 0.012155871727121082685623083829D+00
        w(123) = 0.012167663319667843366755737416D+00
        w(124) = 0.012177643792842880196606249581D+00
        w(125) = 0.012185811661083365425569178819D+00
        w(126) = 0.012192165708627157605870499188D+00
        w(127) = 0.012196704989693764053654538465D+00
        w(128) = 0.012199428828625117371582840212D+00
        w(129) = 0.01220033681998614507777289232D+00
        w(130) = 0.012199428828625117371582840212D+00
        w(131) = 0.012196704989693764053654538465D+00
        w(132) = 0.012192165708627157605870499188D+00
        w(133) = 0.012185811661083365425569178819D+00
        w(134) = 0.012177643792842880196606249581D+00
        w(135) = 0.012167663319667843366755737416D+00
        w(136) = 0.012155871727121082685623083829D+00
        w(137) = 0.012142270770344990738801546574D+00
        w(138) = 0.012126862473800277391553601370D+00
        w(139) = 0.012109649130964635027950450318D+00
        w(140) = 0.012090633303991361438266420912D+00
        w(141) = 0.012069817823327991167612855626D+00
        w(142) = 0.012047205787294992091420946532D+00
        w(143) = 0.012022800561624589927558893338D+00
        w(144) = 0.011996605778959789329711050159D+00
        w(145) = 0.011968625338313666131272065445D+00
        w(146) = 0.011938863404489011222535627643D+00
        w(147) = 0.011907324407458412445505183140D+00
        w(148) = 0.011874013041704866779344562066D+00
        w(149) = 0.011838934265523020964443424791D+00
        w(150) = 0.011802093300281144573421477037D+00
        w(151) = 0.011763495629643945382264331878D+00
        w(152) = 0.011723146998756342723302879656D+00
        w(153) = 0.011681053413388320313049670635D+00
        w(154) = 0.011637221139040985841125311445D+00
        w(155) = 0.011591656700013970380783131035D+00
        w(156) = 0.011544366878434306436012137033D+00
        w(157) = 0.011495358713246929174010288914D+00
        w(158) = 0.011444639499166951104119199270D+00
        w(159) = 0.011392216785593866154247619654D+00
        w(160) = 0.0113380983754878447625379269120D+00
        w(161) = 0.0112822923242082872447042603128D+00
        w(162) = 0.0112248069383148083152535688671D+00
        w(163) = 0.0111656507743308312328559850485D+00
        w(164) = 0.0111048326374699756056269264803D+00
        w(165) = 0.0110423615803254284301924654946D+00
        w(166) = 0.0109782469015224934483083029166D+00
        w(167) = 0.0109124981443345193856719616965D+00
        w(168) = 0.0108451250952624130885928632830D+00
        w(169) = 0.0107761377825779489945556541150D+00
        w(170) = 0.0107055464748310917616231511294D+00
        w(171) = 0.0106333616793215542382761147904D+00
        w(172) = 0.0105595941405348182788823332058D+00
        w(173) = 0.0104842548385428511997370260353D+00
        w(174) = 0.0104073549873697559257355517893D+00
        w(175) = 0.0103289060333225980974485876288D+00
        w(176) = 0.0102489196532876585918958554047D+00
        w(177) = 0.0101674077529923650568895461852D+00
        w(178) = 0.0100843824652331611676814627839D+00
        w(179) = 0.0099998561480695773850435626986D+00
        w(180) = 0.0099138413829847720250916955489D+00
        w(181) = 0.0098263509730128164423854701706D+00
        w(182) = 0.0097373979408330030783691793250D+00
        w(183) = 0.0096469955268314600363329731559D+00
        w(184) = 0.0095551571871303607110514249099D+00
        w(185) = 0.0094618965915850218253881576301D+00
        w(186) = 0.0093672276217491880067391857021D+00
        w(187) = 0.0092711643688088057725325917169D+00
        w(188) = 0.0091737211314845944854270065178D+00
        w(189) = 0.0090749124139037264846862498962D+00
        w(190) = 0.0089747529234409331997949023068D+00
        w(191) = 0.0088732575685293586050755892934D+00
        w(192) = 0.0087704414564414858792445834744D+00
        w(193) = 0.0086663198910404675908861979240D+00
        w(194) = 0.0085609083705021941391459209280D+00
        w(195) = 0.0084542225850084395379670551258D+00
        w(196) = 0.0083462784144114279413811886655D+00
        w(197) = 0.0082370919258701685661946145361D+00
        w(198) = 0.0081266793714589108764118189068D+00
        w(199) = 0.0080150571857480760504667455353D+00
        w(200) = 0.0079022419833580248574070864277D+00
        w(201) = 0.0077882505564860261212726654404D+00
        w(202) = 0.0076730998724067939537782250476D+00
        w(203) = 0.0075568070709469658838993300454D+00
        w(204) = 0.0074393894619338979090297315972D+00
        w(205) = 0.0073208645226191563361371026044D+00
        w(206) = 0.0072012498950770900730828552207D+00
        w(207) = 0.0070805633835788707705149901066D+00
        w(208) = 0.0069588229519423919043121805236D+00
        w(209) = 0.0068360467208584215286561508406D+00
        w(210) = 0.0067122529651934070221351960200D+00
        w(211) = 0.0065874601112693336961737372300D+00
        w(212) = 0.0064616867341210426397202932350D+00
        w(213) = 0.0063349515547314166407936938524D+00
        w(214) = 0.0062072734372448464599330978665D+00
        w(215) = 0.0060786713861593931405204596709D+00
        w(216) = 0.0059491645434980654366600347567D+00
        w(217) = 0.0058187721859596348346566361185D+00
        w(218) = 0.0056875137220494140577838938236D+00
        w(219) = 0.0055554086891904284012033890901D+00
        w(220) = 0.0054224767508154127788846727083D+00
        w(221) = 0.0052887376934400710240953933529D+00
        w(222) = 0.0051542114237180378340642003713D+00
        w(223) = 0.0050189179654779878773297516544D+00
        w(224) = 0.0048828774567433411142588306018D+00
        w(225) = 0.0047461101467350184936945641585D+00
        w(226) = 0.0046086363928577081326523656522D+00
        w(227) = 0.0044704766576701092218388764046D+00
        w(228) = 0.0043316515058396297504806208252D+00
        w(229) = 0.0041921816010820254766367595011D+00
        w(230) = 0.0040520877030864825223229951262D+00
        w(231) = 0.0039113906644266662571543468015D+00
        w(232) = 0.0037701114274582873548537007645D+00
        w(233) = 0.0036282710212037760873102463983D+00
        w(234) = 0.0034858905582247143702551557840D+00
        w(235) = 0.0033429912314827618499065991316D+00
        w(236) = 0.0031995943111899437356540290142D+00
        w(237) = 0.0030557211416493711130936102459D+00
        w(238) = 0.0029113931380877846359302447381D+00
        w(239) = 0.0027666317834818283552560256501D+00
        w(240) = 0.0026214586253808109266552781372D+00
        w(241) = 0.0024758952727301488651840215879D+00
        w(242) = 0.0023299633927021828561308282641D+00
        w(243) = 0.0021836847075455253317921866057D+00
        w(244) = 0.0020370809914723626741694800322D+00
        w(245) = 0.0018901740676190104269878470891D+00
        w(246) = 0.0017429858051468299509941139300D+00
        w(247) = 0.0015955381166175133369701690235D+00
        w(248) = 0.0014478529559255120065233994722D+00
        w(249) = 0.0012999523174235227389668643832D+00
        w(250) = 0.0011518582377826677880963146741D+00
        w(251) = 0.00100359280467969441299468763292D+00
        w(252) = 0.00085517818446696565626595950963D+00
        w(253) = 0.00070663671051592291949335494247D+00
        w(254) = 0.00055799120546880640169677292533D+00
        w(255) = 0.00040926648283531339591138751432D+00
        w(256) = 0.00026049995580176964436806680831D+00
        w(257) = 0.00011191470145601756450862287886D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of N = ', n
        write ( *, '(a)' ) '  Legal values are 1 through 33, ' //
     &    '63/64/65, 127/128/129 and 255/256/257.'
        stop

      end if

      return
      end
      subroutine legendre_set_x1 ( order, xtab, weight )

c*********************************************************************72
c
cc LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
c
c  Integration region:
c
c    [ -1, 1 ]
c
c  Weight function:
c
c    1 + X
c
c  Integral to approximate:
c
c    integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
c
c  Approximate integral:
c
c    sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4,G3S7.
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c    ORDER must be between 1 and 9.
c
c    Output, double precision XTAB(ORDER), the abscissas of the rule.
c
c    Output, double precision WEIGHT(ORDER), the weights of the rule.
c
      implicit none

      integer order

      double precision xtab(order)
      double precision weight(order)

      if ( order .eq. 1 ) then

        xtab(1) =  0.333333333333333333333333333333D+00

        weight(1) = 2.0D+00

      else if ( order .eq. 2 ) then

        xtab(1) = -0.289897948556635619639456814941D+00
        xtab(2) =  0.689897948556635619639456814941D+00

        weight(1) =  0.727834473024091322422523991699D+00
        weight(2) =  1.27216552697590867757747600830D+00

      else if ( order .eq. 3 ) then

        xtab(1) = -0.575318923521694112050483779752D+00
        xtab(2) =  0.181066271118530578270147495862D+00
        xtab(3) =  0.822824080974592105208907712461D+00

        weight(1) =  0.279307919605816490135525088716D+00
        weight(2) =  0.916964425438344986775682378225D+00
        weight(3) =  0.803727654955838523088792533058D+00

      else if ( order .eq. 4 ) then

        xtab(1) = -0.720480271312438895695825837750D+00
        xtab(2) = -0.167180864737833640113395337326D+00
        xtab(3) =  0.446313972723752344639908004629D+00
        xtab(4) =  0.885791607770964635613757614892D+00

        weight(1) =  0.124723883800032328695500588386D+00
        weight(2) =  0.519390190432929763305824811559D+00
        weight(3) =  0.813858272041085443165617903743D+00
        weight(4) =  0.542027653725952464833056696312D+00

      else if ( order .eq. 5 ) then

        xtab(1) = -0.802929828402347147753002204224D+00
        xtab(2) = -0.390928546707272189029229647442D+00
        xtab(3) =  0.124050379505227711989974959990D+00
        xtab(4) =  0.603973164252783654928415726409D+00
        xtab(5) =  0.920380285897062515318386619813D+00

        weight(1) =  0.0629916580867691047411692662740D+00
        weight(2) =  0.295635480290466681402532877367D+00
        weight(3) =  0.585547948338679234792151477424D+00
        weight(4) =  0.668698552377478261966702492391D+00
        weight(5) =  0.387126360906606717097443886545D+00

      else if ( order .eq. 6 ) then

        xtab(1) = -0.853891342639482229703747931639D+00
        xtab(2) = -0.538467724060109001833766720231D+00
        xtab(3) = -0.117343037543100264162786683611D+00
        xtab(4) =  0.326030619437691401805894055838D+00
        xtab(5) =  0.703842800663031416300046295008D+00
        xtab(6) =  0.941367145680430216055899446174D+00

        weight(1) =  0.0349532072544381270240692132496D+00
        weight(2) =  0.175820662202035902032706497222D+00
        weight(3) =  0.394644603562621056482338042193D+00
        weight(4) =  0.563170215152795712476307356284D+00
        weight(5) =  0.542169988926074467362761586552D+00
        weight(6) =  0.289241322902034734621817304499D+00

      else if ( order .eq. 7 ) then

        xtab(1) = -0.887474878926155707068695617935D+00
        xtab(2) = -0.639518616526215270024840114382D+00
        xtab(3) = -0.294750565773660725252184459658D+00
        xtab(4) =  0.0943072526611107660028971153047D+00
        xtab(5) =  0.468420354430821063046421216613D+00
        xtab(6) =  0.770641893678191536180719525865D+00
        xtab(7) =  0.955041227122575003782349000858D+00

        weight(1) =  0.0208574488112296163587654972151D+00
        weight(2) =  0.109633426887493901777324193433D+00
        weight(3) =  0.265538785861965879934591955055D+00
        weight(4) =  0.428500262783494679963649011999D+00
        weight(5) =  0.509563589198353307674937943100D+00
        weight(6) =  0.442037032763498409684482945478D+00
        weight(7) =  0.223869453693964204606248453720D+00

      else if ( order .eq. 8 ) then

        xtab(1) = -0.910732089420060298533757956283D+00
        xtab(2) = -0.711267485915708857029562959544D+00
        xtab(3) = -0.426350485711138962102627520502D+00
        xtab(4) = -0.0903733696068532980645444599064D+00
        xtab(5) =  0.256135670833455395138292079035D+00
        xtab(6) =  0.571383041208738483284917464837D+00
        xtab(7) =  0.817352784200412087992517083851D+00
        xtab(8) =  0.964440169705273096373589797925D+00

        weight(1) =  0.0131807657689951954189692640444D+00
        weight(2) =  0.0713716106239448335742111888042D+00
        weight(3) =  0.181757278018795592332221684383D+00
        weight(4) =  0.316798397969276640481632757440D+00
        weight(5) =  0.424189437743720042818124385645D+00
        weight(6) =  0.450023197883549464687088394417D+00
        weight(7) =  0.364476094545494505382889847132D+00
        weight(8) =  0.178203217446223725304862478136D+00

      else if ( order .eq. 9 ) then

        xtab(1) = -0.927484374233581078117671398464D+00
        xtab(2) = -0.763842042420002599615429776011D+00
        xtab(3) = -0.525646030370079229365386614293D+00
        xtab(4) = -0.236234469390588049278459503207D+00
        xtab(5) =  0.0760591978379781302337137826389D+00
        xtab(6) =  0.380664840144724365880759065541D+00
        xtab(7) =  0.647766687674009436273648507855D+00
        xtab(8) =  0.851225220581607910728163628088D+00
        xtab(9) =  0.971175180702246902734346518378D+00

        weight(1) =  0.00872338834309252349019620448007D+00
        weight(2) =  0.0482400171391415162069086091476D+00
        weight(3) =  0.127219285964216005046760427743D+00
        weight(4) =  0.233604781180660442262926091607D+00
        weight(5) =  0.337433287379681397577000079834D+00
        weight(6) =  0.401235236773473158616600898930D+00
        weight(7) =  0.394134968689382820640692081477D+00
        weight(8) =  0.304297020437232650320317215016D+00
        weight(9) =  0.145112014093119485838598391765D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_SET_X1 - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of ORDER = ', order
        stop

      end if

      return
      end
      subroutine legendre_set_x2 ( order, xtab, weight )

c*********************************************************************72
c
cc LEGENDRE_SET_X2 sets a Gauss-Legendre rule for ( 1 + X )^2 * F(X) on [-1,1].
c
c  Integration region:
c
c    [ -1, 1 ]
c
c  Weight function:
c
c    ( 1 + X )^2
c
c  Integral to approximate:
c
c    integral ( -1 <= X <= 1 ) ( 1 + X )^2 * F(X) dX
c
c  Approximate integral:
c
c    sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c    ORDER must be between 1 and 9.
c
c    Output, double precision XTAB(ORDER), the abscissas of the rule.
c
c    Output, double precision WEIGHT(ORDER), the weights of the rule.
c
      implicit none

      integer order

      double precision xtab(order)
      double precision weight(order)

      if ( order .eq. 1 ) then

        xtab(1) =  0.5D+00

        weight(1) =  2.66666666666666666666666666666D+00

      else if ( order .eq. 2 ) then

        xtab(1) = -0.0883036880224505775998524725910D+00
        xtab(2) =  0.754970354689117244266519139258D+00

        weight(1) =  0.806287056638603444666851075928D+00
        weight(2) =  1.86037961002806322199981559074D+00

      else if ( order .eq. 3 ) then

        xtab(1) = -0.410004419776996766244796955168D+00
        xtab(2) =  0.305992467923296230556472913192D+00
        xtab(3) =  0.854011951853700535688324041976D+00

        weight(1) =  0.239605624068645584091811926047D+00
        weight(2) =  1.16997015407892817602809616291D+00
        weight(3) =  1.25709088851909290654675857771D+00

      else if ( order .eq. 4 ) then

        xtab(1) = -0.591702835793545726606755921586D+00
        xtab(2) = -0.0340945902087350046811467387661D+00
        xtab(3) =  0.522798524896275389882037174551D+00
        xtab(4) =  0.902998901106005341405865485802D+00

        weight(1) =  0.0828179259993445222751812523731D+00
        weight(2) =  0.549071097383384602539010760334D+00
        weight(3) =  1.14767031839371367238662411421D+00
        weight(4) =  0.887107324890223869465850539752D+00

      else if ( order .eq. 5 ) then

        xtab(1) = -0.702108425894032836232448374820D+00
        xtab(2) = -0.268666945261773544694327777841D+00
        xtab(3) =  0.220227225868961343518209179230D+00
        xtab(4) =  0.653039358456608553790815164028D+00
        xtab(5) =  0.930842120163569816951085142737D+00

        weight(1) =  0.0329106016247920636689299329544D+00
        weight(2) =  0.256444805783695354037991444453D+00
        weight(3) =  0.713601289772720001490035944563D+00
        weight(4) =  1.00959169519929190423066348132D+00
        weight(5) =  0.654118274286167343239045863379D+00

      else if ( order .eq. 6 ) then

        xtab(1) = -0.773611232355123732602532012021D+00
        xtab(2) = -0.431362254623427837535325249187D+00
        xtab(3) = -0.0180728263295041680220798103354D+00
        xtab(4) =  0.395126163954217534500188844163D+00
        xtab(5) =  0.736872116684029732026178298518D+00
        xtab(6) =  0.948190889812665614490712786006D+00

        weight(1) =  0.0146486064549543818622276447204D+00
        weight(2) =  0.125762377479560410622810097040D+00
        weight(3) =  0.410316569036929681761034600615D+00
        weight(4) =  0.756617493988329628546336413760D+00
        weight(5) =  0.859011997894245060846045458784D+00
        weight(6) =  0.500309621812647503028212451747D+00

      else if ( order .eq. 7 ) then

        xtab(1) = -0.822366333126005527278634734418D+00
        xtab(2) = -0.547034493182875002223997992852D+00
        xtab(3) = -0.200043026557985860387937545780D+00
        xtab(4) =  0.171995710805880507163425502299D+00
        xtab(5) =  0.518891747903884926692601716998D+00
        xtab(6) =  0.793821941703901970495546427988D+00
        xtab(7) =  0.959734452453198985538996625765D+00

        weight(1) =  0.00714150426951365443207221475404D+00
        weight(2) =  0.0653034050584375560578544725498D+00
        weight(3) =  0.235377690316228918725962815880D+00
        weight(4) =  0.505171029671130381676271523850D+00
        weight(5) =  0.733870426238362032891332767175D+00
        weight(6) =  0.725590596901489156295739839779D+00
        weight(7) =  0.394212014211504966587433032679D+00

      else if ( order .eq. 8 ) then

        xtab(1) = -0.857017929919813794402037235698D+00
        xtab(2) = -0.631543407166567521509503573952D+00
        xtab(3) = -0.339104543648722903660229021109D+00
        xtab(4) = -0.0111941563689783438801237300122D+00
        xtab(5) =  0.316696017045595559454075475675D+00
        xtab(6) =  0.609049663022520165351466780939D+00
        xtab(7) =  0.834198765028697794599267293239D+00
        xtab(8) =  0.967804480896157932935972899807D+00

        weight(1) =  0.00374814227227757804631954025851D+00
        weight(2) =  0.0357961737041152639660521680263D+00
        weight(3) =  0.137974910241879862433949246199D+00
        weight(4) =  0.326515411108352185491692769217D+00
        weight(5) =  0.547577467373226177976217604887D+00
        weight(6) =  0.682278153375510121675529810121D+00
        weight(7) =  0.614544746137780998436053880546D+00
        weight(8) =  0.318231662453524478640851647411D+00

      else if ( order .eq. 9 ) then

        xtab(1) = -0.882491728426548422828684254270D+00
        xtab(2) = -0.694873684026474640346360850039D+00
        xtab(3) = -0.446537143480670863635920316400D+00
        xtab(4) = -0.159388112702326252531544826624D+00
        xtab(5) =  0.141092709224374414981503995427D+00
        xtab(6) =  0.428217823321559204544020866175D+00
        xtab(7) =  0.676480966471850715860378175342D+00
        xtab(8) =  0.863830940812464825046988286026D+00
        xtab(9) =  0.973668228805771018909618924364D+00

        weight(1) =  0.00209009877215570354392734918986D+00
        weight(2) =  0.0205951891648697848186537272448D+00
        weight(3) =  0.0832489326348178964194106978875D+00
        weight(4) =  0.210746247220398685903797568021D+00
        weight(5) =  0.388325022916052063676224499399D+00
        weight(6) =  0.554275165518437673725822282791D+00
        weight(7) =  0.621388553284444032628761363828D+00
        weight(8) =  0.523916296267173054255512857631D+00
        weight(9) =  0.262081160888317771694556320674D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_SET_X2 - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of ORDER = ', order
        stop

      end if

      return
      end
      function lens_half_2d ( func, center, r, theta1, theta2, order )

c*********************************************************************72
c
cc LENS_HALF_2D approximates an integral in a circular half lens in 2D.
c
c  Discussion:
c
c    A circular half lens is formed by drawing a circular arc,
c    and joining its endpoints.
c
c    This rule for a circular half lens simply views the region as 
c    a product region, with a coordinate "S" that varies along the
c    radial direction, and a coordinate "T" that varies in the perpendicular
c    direction, and whose extent varies as a function of S.
c
c    A Gauss-Legendre rule is used to construct a product rule that is
c    applied to the region.  The accuracy of the Gauss-Legendre rule,
c    which is valid for a rectangular product rule region, does not
c    apply straightforwardly to this region, since the limits in the
c    "T" coordinate are being handled implicitly.
c
c    This is simply an application of the QMULT_2D algorithm of Stroud.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of two
c    variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision CENTER(2), the center of the circle.
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles of the rays
c    that begin and end the arc.
c
c    Input, integer ORDER, the order of the Gauss-Legendre rule
c    to be used.  Legal values include 1 through 20, 32 or 64.
c
c    Output, double precision LENS_HALF_2D, the approximate value
c    of the integral of the function over the half lens.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer order

      double precision ax
      double precision ay
      double precision bx
      double precision by
      double precision center(dim_num)
      double precision cx
      double precision cy
      double precision dx
      double precision dy
      double precision func
      external func
      integer i
      integer j
      double precision lens_half_2d
      double precision quad
      double precision r
      double precision s_length
      double precision sx
      double precision sy
      double precision t_length
      double precision tdirx
      double precision tdiry
      double precision theta1
      double precision theta2
      double precision thi
      double precision tx
      double precision ty
      double precision w1
      double precision w2
      double precision weight(order)
      double precision xtab(order)
c
c  Retrieve the Legendre rule of the given order.
c
      call legendre_set ( order, xtab, weight )
c
c  Determine the points A (on the secant) and B (on the circumference)
c  that will form the "S" direction.
c
      ax = center(1) + r * 0.5D+00 * ( cos ( theta1 ) + cos ( theta2 ) )
      ay = center(2) + r * 0.5D+00 * ( sin ( theta1 ) + sin ( theta2 ) )

      bx = center(1) + r * cos ( 0.5D+00 * ( theta1 + theta2 ) )
      by = center(2) + r * sin ( 0.5D+00 * ( theta1 + theta2 ) )
c
c  Find the length of the line between A and B.
c
      s_length = sqrt ( ( ax - bx ) * ( ax - bx ) 
     &                + ( ay - by ) * ( ay - by ) )

      if ( s_length .eq. 0.0D+00 ) then
        lens_half_2d = 0.0D+00
        return
      end if
c
c  Determine the unit vector in the T direction.
c
      tdirx = ( ay - by ) / s_length
      tdiry = ( bx - ax ) / s_length

      quad = 0.0D+00

      do i = 1, order

        w1 = 0.5D+00 * s_length * weight(i)
c
c  Map the quadrature point to an S coordinate.
c
        sx = ( ( 1.0D+00 - xtab(i) ) * ax   
     &       + ( 1.0D+00 + xtab(i) ) * bx ) 
     &       /   2.0D+00
        sy = ( ( 1.0D+00 - xtab(i) ) * ay   
     &       + ( 1.0D+00 + xtab(i) ) * by ) 
     &       /   2.0D+00
c
c  Determine the length of the line in the T direction, from the
c  S axis to the circle circumference.
c
        thi = sqrt ( ( r - 0.25D+00 * ( 1.0D+00 - xtab(i) ) * s_length ) 0
     &                       *        ( 1.0D+00 - xtab(i) ) * s_length )
c 
c  Determine the maximum and minimum T coordinates by going
c  up and down in the T direction from the S axis.
c
        cx = sx + tdirx * thi
        cy = sy + tdiry * thi
        dx = sx - tdirx * thi
        dy = sy - tdiry * thi
c
c  Find the length of the T direction.
c
        t_length = sqrt ( ( cx - dx ) * ( cx - dx ) 
     &                  + ( cy - dy ) * ( cy - dy ) )

        do j = 1, order

          w2 = 0.5D+00 * t_length * weight(j)
c
c  Map the quadrature point to a T coordinate.
c
          tx = ( ( 1.0D+00 - xtab(j) ) * cx   
     &         + ( 1.0D+00 + xtab(j) ) * dx ) 
     &         /   2.0D+00
          ty = ( ( 1.0D+00 - xtab(j) ) * cy   
     &         + ( 1.0D+00 + xtab(j) ) * dy ) 
     &         /   2.0D+00

          quad = quad + w1 * w2 * func ( tx, ty )

        end do

      end do

      lens_half_2d = quad

      return
      end
      function lens_half_area_2d ( r, theta1, theta2 )

c*********************************************************************72
c
cc LENS_HALF_AREA_2D returns the area of a circular half lens in 2D.
c
c  Discussion:
c
c    A circular half lens is formed by drawing a circular arc, 
c    and joining its endpoints.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision THETA1, THETA2, the angles of the rays
c    that begin and end the arc.
c
c    Output, double precision LENS_HALF_AREA_2D, the area of the half lens.
c
      implicit none

      double precision circle_sector_area_2d
      double precision circle_triangle_area_2d
      double precision lens_half_area_2d
      double precision r
      double precision sector
      double precision theta1
      double precision theta2
      double precision triangle

      sector = circle_sector_area_2d ( r, theta1, theta2 )
      triangle = circle_triangle_area_2d ( r, theta1, theta2 )
      lens_half_area_2d = sector - triangle

      return
      end
      function lens_half_h_area_2d ( r, h )

c*********************************************************************72
c
cc LENS_HALF_H_AREA_2D returns the area of a circular half lens in 2D.
c
c  Discussion:
c
c    A circular half lens is formed by drawing a circular arc, and joining 
c    its endpoints.
c
c    This particular half lens is described by the "height" of the region.  
c    In other words, the half lens is the region that would be submerged 
c    if a circle of radius R were standing in water of depth H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision H, the height of the half lens region.
c
c    Output, double precision LENS_HALF_H_AREA_2D, the area of the half lens.
c
      implicit none

      double precision angle
      double precision area
      double precision h
      double precision half_width
      double precision lens_half_h_area_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision sector
      double precision triangle

      if ( h <= 0.0D+00 ) then

        area = 0.0D+00

      else if ( 2.0D+00 * r <= h ) then

        area = pi * r * r

      else

        half_width = sqrt ( h * ( 2.0D+00 * r - h ) )
        angle = 2.0D+00 * atan2 ( half_width, r - h )
        sector = r * r * angle / 2.0D+00
        triangle = ( r - h ) * half_width
        area = sector - triangle

      end if

      lens_half_h_area_2d = area

      return
      end
      function lens_half_w_area_2d ( r, w )

c*********************************************************************72
c
cc LENS_HALF_W_AREA_2D returns the area of a circular half lens in 2D.
c
c  Discussion:
c
c    A half lens is formed by drawing a circular arc, and joining its endpoints.
c    This half lens is described by the "width" of the region.  In other words,
c    it is the portion of the circle under water if the width
c    of the water surface is W.  There are two possible values for this
c    area, A and (PI*R*R-A).  The routine returns the smaller of the 
c    two values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision W, the width of the half lens region.
c
c    Output, double precision LENS_HALF_W_AREA_2D, the area of the half lens.
c
      implicit none

      double precision angle
      double precision area
      double precision h
      double precision half_width
      double precision lens_half_w_area_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision sector
      double precision triangle
      double precision w

      if ( w .le. 0.0D+00 ) then

        area = 0.0D+00

      else if ( 2.0D+00 * r .le. w ) then

        area = 0.5D+00 * pi * r * r

      else

        half_width = 0.5D+00 * w
        h = r - sqrt ( r * r - half_width * half_width )
        angle = 2.0D+00 * atan2 ( half_width, r - h )
        sector = r * r * angle / 2.0D+00
        triangle = ( r - h ) * half_width
        area = sector - triangle

      end if

      lens_half_w_area_2d = area

      return
      end
      subroutine monomial_value ( dim_num, point_num, x, expon, value )

c*********************************************************************72
c
cc MONOMIAL_VALUE evaluates a monomial.
c
c  Discussion:
c
c    This routine evaluates a monomial of the form
c
c      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
c
c    where the exponents are nonnegative integers.  Note that
c    if the combination 0^0 is encountered, it should be treated
c    as 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer POINT_NUM, the number of points at which the
c    monomial is to be evaluated.
c
c    Input, double precision X(DIM_NUM,POINT_NUM), the point coordinates.
c
c    Input, integer EXPON(DIM_NUM), the exponents.
c
c    Output, double precision VALUE(POINT_NUM), the value of the monomial.
c
      implicit none

      integer dim_num
      integer point_num

      integer dim
      integer expon(dim_num)
      integer j 
      double precision value(point_num)
      double precision x(dim_num,point_num)

      do j = 1, point_num
        value(j) = 1.0D+00
      end do

      do dim = 1, dim_num
        if ( 0 .ne. expon(dim) ) then
          do j = 1, point_num
            value(j) = value(j) * x(dim,j)**expon(dim)
          end do
        end if
      end do

      return
      end
      subroutine octahedron_unit_nd ( func, n, result )

c*********************************************************************72
c
cc OCTAHEDRON_UNIT_ND approximates integrals in the unit octahedron in ND.
c
c  Integration region:
c
c    sum ( abs ( X(1:N) ) ) <= 1.
c
c  Discussion:
c
c    A 2*N point 3rd degree formula is used, Stroud number GN:3-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which is to be integrated, of the form:
c
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the octahedron.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      integer j
      double precision octahedron_unit_volume_nd
      double precision quad
      double precision r
      double precision result
      double precision volume
      double precision w
      double precision x(n)

      w = 1.0D+00 / dble ( 2 * n )
      r = sqrt ( dble ( 2 * n ) / dble ( ( n + 1 ) * ( n + 2 ) ) )

      do i = 1, n
        x(i) = 0.0D+00
      end do

      quad = 0.0D+00
      do i = 1, n
        x(i) = r
        do j = 1, 2
          quad = quad + w * func ( n, x )
          x(i) = - x(i)
        end do
        x(i) = 0.0D+00
      end do

      volume = octahedron_unit_volume_nd ( n )
      result = quad * volume

      return
      end
      function octahedron_unit_volume_nd ( n )

!*****************************************************************************80
!
!! OCTAHEDRON_UNIT_VOLUME_ND returns the volume of the unit octahedron in ND.
!
!  Integration region:
!
!    sum ( abs ( X(1:N) ) ) <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, double precision OCTAHEDRON_UNIT_VOLUME_ND, the volume of
!    the unit octahedron.
!
      implicit none

      integer i
      integer n
      double precision octahedron_unit_volume_nd
      double precision volume

      volume = 1.0D+00
      do i = 1, n
        volume = volume * 2.0D+00 / dble ( i )
      end do

      octahedron_unit_volume_nd = volume

      return
      end
      function parallelipiped_volume_3d ( x, y, z )

c*********************************************************************72
c
cc PARALLELIPIPED_VOLUME_3D returns the volume of a parallelipiped in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(4), Y(4), Z(4), the coordinates of one corner
c    of the parallelipiped, and its 3 immediate neighbors.
c
c    Output, double precision PARALLELIPIPED_VOLUME_3D, the volume of
c    the parallelipiped.
c
      implicit none

      double precision parallelipiped_volume_3d
      double precision x(4)
      double precision y(4)
      double precision z(4)

      parallelipiped_volume_3d = abs ( 
     &  ( z(2) - z(1) ) * ( y(4) * x(3) - y(3) * x(4) ) + 
     &  ( z(3) - z(1) ) * ( x(4) * y(2) - x(2) * y(4) ) + 
     &  ( z(4) - z(1) ) * ( x(2) * y(3) - x(3) * y(2) ) + 
     &  ( z(3) - z(2) ) * ( y(4) * x(1) - y(1) * x(4) ) + 
     &  ( z(4) - z(2) ) * ( x(3) * y(1) - x(1) * y(3) ) + 
     &  ( z(4) - z(3) ) * ( x(1) * y(2) - x(2) * y(1) ) )

      return
      end
      function parallelipiped_volume_nd ( n, v )

c*********************************************************************72
c
cc PARALLELIPIPED_VOLUME_ND returns the volume of a parallelipiped in ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ( kind = 4 ) N, the dimension of the space.
c
c    Input, double precision V(N,N+1), the
c    N+1 columns of V contains the N coordinates of one of the
c    "corners" of the parallelipiped.
c
c    Output, double precision PARALLELIPIPED_VOLUME_ND, the volume of
c    the parallelipiped.
c
      implicit none

      integer n

      double precision det
      integer i
      integer info
      integer j
      double precision parallelipiped_volume_nd
      integer pivot(n)
      double precision v(n,n+1)
      double precision w(n,n)
c
c  Compute the volume of the N-dimensional parallelipiped.
c
      do j = 1, n
        do i = 1, n
          w(i,j) = v(i,j+1) - v(i,1)
        end do
      end do

      call r8ge_fa ( n, w, pivot, info )

      if ( info .ne. 0 ) then
        parallelipiped_volume_nd = 0.0D+00
        return
      end if

      call r8ge_det ( n, w, pivot, det )

      parallelipiped_volume_nd = abs ( det )

      return
      end
      subroutine polygon_1_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_1_2D integrates the function 1 over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
c      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c    The integral of 1 over a polygon is the area of the polygon.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_1_2D - Warningc'
        write ( *, '(a)' )
     &    '  The number of vertices must be at least 3.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result + 0.5D+00 * ( v(1,i) + v(1,im1) )
     &                            * ( v(2,i) - v(2,im1) )

      end do

      return
      end
      subroutine polygon_x_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_X_2D integrates the function X over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
c      ( X(I)*X(I) + X(I) * X(I-1) + X(I-1)*X(I-1) ) * ( Y(I) - Y(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_X_2D - Warning!'
        write ( *, '(a)' ) 
     &    '  The number of vertices must be at least 3.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result 
     &    + ( v(1,i)**2 + v(1,i) * v(1,im1) + v(1,im1)**2 )
     &    * ( v(2,i) - v(2,im1) )

      end do

      result = result / 6.0D+00

      return
      end
      subroutine polygon_xx_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
c      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
c      * ( Y(I) - Y(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in
c    counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_XX_2D - Warning!'
        write ( *, '(a)' ) 
     &    '  The number of vertices must be at least 3.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result + ( v(1,i)**3 + v(1,i)**2 * v(1,im1) 
     &    + v(1,i) * v(1,im1)**2 + v(1,im1)**3 ) * ( v(2,i) - v(2,im1) )

      end do

      result = result / 12.0D+00

      return
      end
      subroutine polygon_xy_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
c      ( Y(I)   * ( 3 * X(I)^2 + 2 * X(I) * X(I-1) +     X(I-1)^2 )
c      + Y(I-1) * (     X(I)^2 + 2 * X(I) * X(I-1) + 3 * X(I-1)^2 ) )
c      * ( Y(I) - Y(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in
c    counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_XY_2D - Warning!'
        write ( *, '(a)' ) 
     &    '  The number of vertices must be at least 3.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result + ( 
     &    v(2,i) * ( 3.0D+00 * v(1,i)**2 + 2.0D+00 * v(1,i) * v(1,im1) 
     &    + v(1,im1)**2 ) 
     &    + v(2,im1) * ( v(1,i)**2 + 2.0D+00 * v(1,i) * v(1,im1) 
     &    + 3.0D+00 * v(1,im1)**2 ) ) * ( v(2,i) - v(2,im1) )

      end do

      result = result / 24.0D+00

      return
      end
      subroutine polygon_y_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_Y_2D integrates the function Y over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
c      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in
c    counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_Y_2D - Warning!'
        write ( *, '(a)' ) 
     &    '  The number of vertices must be at least 3.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result 
     &    - ( v(2,i)**2 + v(2,i) * v(2,im1) + v(2,im1)**2 ) 
     &    * ( v(1,i) - v(1,im1) )

      end do

      result = result / 6.0D+00

      return
      end
      subroutine polygon_yy_2d ( n, v, result )

c*********************************************************************72
c
cc POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
c
c  Discussion:
c
c    The polygon is bounded by the points (X(1:N), Y(1:N)).
c
c    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
c      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
c      * ( X(I) - X(I-1) )
c
c    where X(0) and Y(0) should be replaced by X(N) and Y(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    SF Bockman,
c    Generalizing the Formula for Areas of Polygons to Moments,
c    American Mathematical Society Monthly,
c    1989, pages 131-132.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c    N should be at least 3 for a nonzero result.
c
c    Input, double precision V(2,N), the coordinates of the vertices
c    of the polygon.  These vertices should be given in
c    counter clockwise order.
c
c    Output, double precision RESULT, the value of the integral.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      integer im1
      double precision result
      double precision v(dim_num,n)

      result = 0.0D+00

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_YY_2D - Warning!'
        write ( *, '(a)' ) '  The number of polygonal vertices must be '
        write ( *, '(a,i8)' ) 
     &    '  at least 3, but the input polygon has N = ', n
        return
      end if

      do i = 1, n

        if ( i .eq. 1 ) then
          im1 = n
        else
          im1 = i - 1
        end if

        result = result - ( v(2,i)**3 + v(2,i)**2 * v(2,im1) 
     &    + v(2,i) * v(2,im1)**2 + v(2,im1)**3 ) * ( v(1,i) - v(1,im1) )

      end do

      result = result / 12.0D+00

      return
      end
      subroutine pyramid_unit_o01_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O01_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 1 point degree 1 formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision func
      external func
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w
      double precision x
      double precision y
      double precision z
c
c  Quadrature.
c
      quad = 0.0D+00

      x = 0.0D+00
      y = 0.0D+00
      z = 1.0D+00 / 4.0D+00
      w = 1.0D+00

      quad = quad + w * func ( x, y, z )
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o05_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O05_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 5 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 5 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w /
     &  0.21093750000000000000D+00, 
     &  0.21093750000000000000D+00, 
     &  0.21093750000000000000D+00, 
     &  0.21093750000000000000D+00, 
     &  0.15625000000000000000D+00 /
      data x /
     & -0.48686449556014765641D+00, 
     &  0.48686449556014765641D+00, 
     &  0.48686449556014765641D+00, 
     & -0.48686449556014765641D+00, 
     &  0.00000000000000000000D+00 /
      data y /
     & -0.48686449556014765641D+00, 
     & -0.48686449556014765641D+00, 
     &  0.48686449556014765641D+00, 
     &  0.48686449556014765641D+00, 
     &  0.00000000000000000000D+00 /
      data z /
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.70000000000000000000D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o06_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O06_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 6 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 6 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w /
     &  0.21000000000000000000D+00, 
     &  0.21000000000000000000D+00, 
     &  0.21000000000000000000D+00, 
     &  0.21000000000000000000D+00, 
     &  0.06000000000000000000D+00, 
     &  0.10000000000000000000D+00 /
      data x /
     & -0.48795003647426658968D+00, 
     &  0.48795003647426658968D+00, 
     &  0.48795003647426658968D+00, 
     & -0.48795003647426658968D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00 /
      data y /
     & -0.48795003647426658968D+00, 
     & -0.48795003647426658968D+00, 
     &  0.48795003647426658968D+00, 
     &  0.48795003647426658968D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00 /
      data z /
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.16666666666666666667D+00, 
     &  0.58333333333333333333D+00, 
     &  0.75000000000000000000D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o08_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O08_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    An 8 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 8 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w /
     &  0.075589411559869072938D+00, 
     &  0.075589411559869072938D+00, 
     &  0.075589411559869072938D+00, 
     &  0.075589411559869072938D+00, 
     &  0.17441058844013092706D+00, 
     &  0.17441058844013092706D+00, 
     &  0.17441058844013092706D+00, 
     &  0.17441058844013092706D+00 /
      data x /
     & -0.26318405556971359557D+00, 
     &  0.26318405556971359557D+00, 
     &  0.26318405556971359557D+00, 
     & -0.26318405556971359557D+00, 
     & -0.50661630334978742377D+00, 
     &  0.50661630334978742377D+00, 
     &  0.50661630334978742377D+00, 
     & -0.50661630334978742377D+00 /
      data y /
     & -0.26318405556971359557D+00, 
     & -0.26318405556971359557D+00, 
     &  0.26318405556971359557D+00, 
     &  0.26318405556971359557D+00, 
     & -0.50661630334978742377D+00, 
     & -0.50661630334978742377D+00, 
     &  0.50661630334978742377D+00, 
     &  0.50661630334978742377D+00 /
      data z /
     &  0.54415184401122528880D+00, 
     &  0.54415184401122528880D+00, 
     &  0.54415184401122528880D+00, 
     &  0.54415184401122528880D+00, 
     &  0.12251482265544137787D+00, 
     &  0.12251482265544137787D+00, 
     &  0.12251482265544137787D+00, 
     &  0.12251482265544137787D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o08b_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O08B_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    An 8 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 8 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w /
     &  0.16438287736328777572D+00, 
     &  0.16438287736328777572D+00, 
     &  0.16438287736328777572D+00, 
     &  0.16438287736328777572D+00, 
     &  0.085617122636712224276D+00, 
     &  0.085617122636712224276D+00, 
     &  0.085617122636712224276D+00, 
     &  0.085617122636712224276D+00 /
      data x /
     & -0.51197009372656270107D+00, 
     &  0.51197009372656270107D+00, 
     &  0.51197009372656270107D+00, 
     & -0.51197009372656270107D+00, 
     & -0.28415447557052037456D+00, 
     &  0.28415447557052037456D+00, 
     &  0.28415447557052037456D+00, 
     & -0.28415447557052037456D+00 /
      data y /
     & -0.51197009372656270107D+00, 
     & -0.51197009372656270107D+00, 
     &  0.51197009372656270107D+00, 
     &  0.51197009372656270107D+00, 
     & -0.28415447557052037456D+00, 
     & -0.28415447557052037456D+00, 
     &  0.28415447557052037456D+00, 
     &  0.28415447557052037456D+00 /
      data z /
     &  0.11024490204163285720D+00, 
     &  0.11024490204163285720D+00, 
     &  0.11024490204163285720D+00, 
     &  0.11024490204163285720D+00, 
     &  0.518326526529795714229D+00, 
     &  0.518326526529795714229D+00, 
     &  0.518326526529795714229D+00, 
     &  0.518326526529795714229D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o09_3d ( func, result )

c*********************************************************************72
c
cc PYRAMID_UNIT_O09_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 9 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 9 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w /
     &  0.13073389672275944791D+00, 
     &  0.13073389672275944791D+00, 
     &  0.13073389672275944791D+00, 
     &  0.13073389672275944791D+00, 
     &  0.10989110327724055209D+00, 
     &  0.10989110327724055209D+00, 
     &  0.10989110327724055209D+00, 
     &  0.10989110327724055209D+00, 
     &  0.03750000000000000000D+00 /
      data x /
     & -0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     & -0.52966422253852215131D+00, 
     & -0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     & -0.34819753825720418039D+00, 
     &  0.00000000000000000000D+00 /
      data y /
     & -0.52966422253852215131D+00, 
     & -0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     & -0.34819753825720418039D+00, 
     & -0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.00000000000000000000D+00 /
      data z /
     &  0.08176876558246862335D+00, 
     &  0.08176876558246862335D+00, 
     &  0.08176876558246862335D+00, 
     &  0.08176876558246862335D+00, 
     &  0.400374091560388519511D+00, 
     &  0.400374091560388519511D+00, 
     &  0.400374091560388519511D+00, 
     &  0.400374091560388519511D+00, 
     &  0.83333333333333333333D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o09_3d ( func, result )

c*****************************************************************************80
c
cc PYRAMID_UNIT_O09_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 9 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 9 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w / 
     & 0.13073389672275944791D+00, 
     & 0.13073389672275944791D+00, 
     & 0.13073389672275944791D+00, 
     & 0.13073389672275944791D+00, 
     & 0.10989110327724055209D+00, 
     & 0.10989110327724055209D+00, 
     & 0.10989110327724055209D+00, 
     & 0.10989110327724055209D+00, 
     & 0.03750000000000000000D+00 /
      data x / 
     & -0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     & -0.52966422253852215131D+00, 
     & -0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     & -0.34819753825720418039D+00, 
     &  0.00000000000000000000D+00 /
      data y / 
     & -0.52966422253852215131D+00, 
     & -0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     &  0.52966422253852215131D+00, 
     & -0.34819753825720418039D+00, 
     & -0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.34819753825720418039D+00, 
     &  0.00000000000000000000D+00 /
      data z / 
     & 0.08176876558246862335D+00, 
     & 0.08176876558246862335D+00, 
     & 0.08176876558246862335D+00, 
     & 0.08176876558246862335D+00, 
     & 0.400374091560388519511D+00, 
     & 0.400374091560388519511D+00, 
     & 0.400374091560388519511D+00, 
     & 0.400374091560388519511D+00, 
     & 0.83333333333333333333D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o13_3d ( func, result )

c*****************************************************************************80
c
cc PYRAMID_UNIT_O13_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 13 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 13 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w / 
     & 0.063061594202898550725D+00, 
     & 0.063061594202898550725D+00, 
     & 0.063061594202898550725D+00, 
     & 0.063061594202898550725D+00, 
     & 0.042101946815575556199D+00, 
     & 0.042101946815575556199D+00, 
     & 0.042101946815575556199D+00, 
     & 0.042101946815575556199D+00, 
     & 0.13172030707666776585D+00, 
     & 0.13172030707666776585D+00, 
     & 0.13172030707666776585D+00, 
     & 0.13172030707666776585D+00, 
     & 0.05246460761943250889D+00 /
      data x / 
     & -0.38510399211870384331D+00, 
     &  0.38510399211870384331D+00, 
     &  0.38510399211870384331D+00, 
     & -0.38510399211870384331D+00, 
     & -0.40345831960728204766D+00, 
     &  0.40345831960728204766D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     & -0.53157877436961973359D+00, 
     &  0.53157877436961973359D+00, 
     &  0.53157877436961973359D+00, 
     & -0.53157877436961973359D+00, 
     &  0.00000000000000000000D+00 /
      data y / 
     & -0.38510399211870384331D+00, 
     & -0.38510399211870384331D+00, 
     &  0.38510399211870384331D+00, 
     &  0.38510399211870384331D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     & -0.40345831960728204766D+00, 
     &  0.40345831960728204766D+00, 
     & -0.53157877436961973359D+00, 
     & -0.53157877436961973359D+00, 
     &  0.53157877436961973359D+00, 
     &  0.53157877436961973359D+00, 
     &  0.00000000000000000000D+00 /
      data z / 
     & 0.428571428571428571429D+00, 
     & 0.428571428571428571429D+00, 
     & 0.428571428571428571429D+00, 
     & 0.428571428571428571429D+00, 
     & 0.33928571428571428571D+00, 
     & 0.33928571428571428571D+00, 
     & 0.33928571428571428571D+00, 
     & 0.33928571428571428571D+00, 
     & 0.08496732026143790850D+00, 
     & 0.08496732026143790850D+00, 
     & 0.08496732026143790850D+00, 
     & 0.08496732026143790850D+00, 
     & 0.76219701803768503595D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o18_3d ( func, result )

c*****************************************************************************80
c
cc PYRAMID_UNIT_O18_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    An 18 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 18 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w / 
     & 0.023330065296255886709D+00, 
     & 0.037328104474009418735D+00, 
     & 0.023330065296255886709D+00, 
     & 0.037328104474009418735D+00, 
     & 0.059724967158415069975D+00, 
     & 0.037328104474009418735D+00, 
     & 0.023330065296255886709D+00, 
     & 0.037328104474009418735D+00, 
     & 0.023330065296255886709D+00, 
     & 0.05383042853090460712D+00, 
     & 0.08612868564944737139D+00, 
     & 0.05383042853090460712D+00, 
     & 0.08612868564944737139D+00, 
     & 0.13780589703911579422D+00, 
     & 0.08612868564944737139D+00, 
     & 0.05383042853090460712D+00, 
     & 0.08612868564944737139D+00, 
     & 0.05383042853090460712D+00 /
      data x / 
     & -0.35309846330877704481D+00, 
     &  0.00000000000000000000D+00, 
     &  0.35309846330877704481D+00, 
     & -0.35309846330877704481D+00, 
     &  0.00000000000000000000D+00, 
     &  0.35309846330877704481D+00, 
     & -0.35309846330877704481D+00, 
     &  0.00000000000000000000D+00, 
     &  0.35309846330877704481D+00, 
     & -0.67969709567986745790D+00, 
     &  0.00000000000000000000D+00, 
     &  0.67969709567986745790D+00, 
     & -0.67969709567986745790D+00, 
     &  0.00000000000000000000D+00, 
     &  0.67969709567986745790D+00, 
     & -0.67969709567986745790D+00, 
     &  0.00000000000000000000D+00, 
     &  0.67969709567986745790D+00 /
      data y / 
     & -0.35309846330877704481D+00, 
     & -0.35309846330877704481D+00, 
     & -0.35309846330877704481D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.35309846330877704481D+00, 
     &  0.35309846330877704481D+00, 
     &  0.35309846330877704481D+00, 
     & -0.67969709567986745790D+00, 
     & -0.67969709567986745790D+00, 
     & -0.67969709567986745790D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.67969709567986745790D+00, 
     &  0.67969709567986745790D+00, 
     &  0.67969709567986745790D+00 /
      data z / 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.544151844011225288800D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00, 
     & 0.12251482265544137787D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o27_3d ( func, result )

c*****************************************************************************80
c
cc PYRAMID_UNIT_O27_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    A 27 point formula is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carlos Felippa,
c    A compendium of FEM integration formulas for symbolic work,
c    Engineering Computation,
c    Volume 21, Number 8, 2004, pages 867-890.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 27 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w / 
     & 0.036374157653908938268D+00, 
     & 0.05819865224625430123D+00, 
     & 0.036374157653908938268D+00, 
     & 0.05819865224625430123D+00, 
     & 0.09311784359400688197D+00, 
     & 0.05819865224625430123D+00, 
     & 0.036374157653908938268D+00, 
     & 0.05819865224625430123D+00, 
     & 0.036374157653908938268D+00, 
     & 0.033853303069413431019D+00, 
     & 0.054165284911061489631D+00, 
     & 0.033853303069413431019D+00, 
     & 0.054165284911061489631D+00, 
     & 0.08666445585769838341D+00, 
     & 0.054165284911061489631D+00, 
     & 0.033853303069413431019D+00, 
     & 0.054165284911061489631D+00, 
     & 0.033853303069413431019D+00, 
     & 0.006933033103838124540D+00, 
     & 0.011092852966140999264D+00, 
     & 0.006933033103838124540D+00, 
     & 0.011092852966140999264D+00, 
     & 0.017748564745825598822D+00, 
     & 0.011092852966140999264D+00, 
     & 0.006933033103838124540D+00, 
     & 0.011092852966140999264D+00, 
     & 0.006933033103838124540D+00 /
      data x / 
     & -0.7180557413198889387D+00, 
     &  0.00000000000000000000D+00, 
     &  0.7180557413198889387D+00, 
     & -0.7180557413198889387D+00, 
     &  0.00000000000000000000D+00, 
     &  0.7180557413198889387D+00, 
     & -0.7180557413198889387D+00, 
     &  0.00000000000000000000D+00, 
     &  0.7180557413198889387D+00, 
     & -0.50580870785392503961D+00, 
     &  0.00000000000000000000D+00, 
     &  0.50580870785392503961D+00, 
     & -0.50580870785392503961D+00, 
     &  0.00000000000000000000D+00, 
     &  0.50580870785392503961D+00, 
     & -0.50580870785392503961D+00, 
     &  0.00000000000000000000D+00, 
     &  0.50580870785392503961D+00, 
     & -0.22850430565396735360D+00, 
     &  0.00000000000000000000D+00, 
     &  0.22850430565396735360D+00, 
     & -0.22850430565396735360D+00, 
     &  0.00000000000000000000D+00, 
     &  0.22850430565396735360D+00, 
     & -0.22850430565396735360D+00, 
     &  0.00000000000000000000D+00, 
     &  0.22850430565396735360D+00 /
      data y / 
     & -0.7180557413198889387D+00, 
     & -0.7180557413198889387D+00, 
     & -0.7180557413198889387D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.7180557413198889387D+00, 
     &  0.7180557413198889387D+00, 
     &  0.7180557413198889387D+00, 
     & -0.50580870785392503961D+00, 
     & -0.50580870785392503961D+00, 
     & -0.50580870785392503961D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.50580870785392503961D+00, 
     &  0.50580870785392503961D+00, 
     &  0.50580870785392503961D+00, 
     & -0.22850430565396735360D+00, 
     & -0.22850430565396735360D+00, 
     & -0.22850430565396735360D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.00000000000000000000D+00, 
     &  0.22850430565396735360D+00, 
     &  0.22850430565396735360D+00, 
     &  0.22850430565396735360D+00 /
      data z / 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00,  
     & 0.07299402407314973216D+00, 
     & 0.07299402407314973216D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.34700376603835188472D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00, 
     & 0.70500220988849838312D+00 /
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      subroutine pyramid_unit_o48_3d ( func, result )

c*****************************************************************************80
c
cc PYRAMID_UNIT_O48_3D approximates an integral inside the unit pyramid in 3D.
c
c  Discussion:
c
c    An 48 point degree 7 formula, Stroud CN:C2:7-1, is used.
c
c    The (X,Y,Z) integration region can be represented as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c    When Z is zero, the integration region is a square lying in the (X,Y) 
c    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
c    radius of the square diminishes, and when Z reaches 1, the square has 
c    contracted to the single point (0,0,1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 48 )

      double precision func
      external func
      integer i
      double precision pyramid_unit_volume_3d
      double precision quad
      double precision result
      double precision volume
      double precision w(order)
      double precision x(order)
      double precision y(order)
      double precision z(order)

      save w
      save x
      save y
      save z

      data w / 
     & 2.01241939442682455D-002, 
     & 2.01241939442682455D-002, 
     & 2.01241939442682455D-002, 
     & 2.01241939442682455D-002, 
     & 2.60351137043010779D-002, 
     & 2.60351137043010779D-002, 
     & 2.60351137043010779D-002, 
     & 2.60351137043010779D-002, 
     & 1.24557795239745531D-002, 
     & 1.24557795239745531D-002, 
     & 1.24557795239745531D-002, 
     & 1.24557795239745531D-002, 
     & 1.87873998794808156D-003, 
     & 1.87873998794808156D-003, 
     & 1.87873998794808156D-003, 
     & 1.87873998794808156D-003, 
     & 4.32957927807745280D-002, 
     & 4.32957927807745280D-002, 
     & 4.32957927807745280D-002, 
     & 4.32957927807745280D-002, 
     & 1.97463249834127288D-002, 
     & 1.97463249834127288D-002, 
     & 1.97463249834127288D-002, 
     & 1.97463249834127288D-002, 
     & 5.60127223523590526D-002, 
     & 5.60127223523590526D-002, 
     & 5.60127223523590526D-002, 
     & 5.60127223523590526D-002, 
     & 2.55462562927473852D-002, 
     & 2.55462562927473852D-002, 
     & 2.55462562927473852D-002, 
     & 2.55462562927473852D-002, 
     & 2.67977366291788643D-002, 
     & 2.67977366291788643D-002, 
     & 2.67977366291788643D-002, 
     & 2.67977366291788643D-002, 
     & 1.22218992265373354D-002, 
     & 1.22218992265373354D-002, 
     & 1.22218992265373354D-002, 
     & 1.22218992265373354D-002, 
     & 4.04197740453215038D-003, 
     & 4.04197740453215038D-003, 
     & 4.04197740453215038D-003, 
     & 4.04197740453215038D-003, 
     & 1.84346316995826843D-003, 
     & 1.84346316995826843D-003, 
     & 1.84346316995826843D-003, 
     & 1.84346316995826843D-003 /
      data x / 
     &  0.88091731624450909D+00,     
     & -0.88091731624450909D+00,     
     &  0.0000000000000000D+00,     
     &  0.0000000000000000D+00,     
     &  0.70491874112648223D+00,     
     & -0.70491874112648223D+00,     
     &  0.0000000000000000D+00,     
     &  0.0000000000000000D+00,     
     &  0.44712732143189760D+00,     
     & -0.44712732143189760D+00,     
     &  0.0000000000000000D+00,     
     &  0.0000000000000000D+00,     
     &  0.18900486065123448D+00,     
     & -0.18900486065123448D+00,     
     &  0.0000000000000000D+00,     
     &  0.0000000000000000D+00,     
     &  0.36209733410322176D+00,     
     & -0.36209733410322176D+00,     
     & -0.36209733410322176D+00,     
     &  0.36209733410322176D+00,     
     &  0.76688932060387538D+00,     
     & -0.76688932060387538D+00,     
     & -0.76688932060387538D+00,     
     &  0.76688932060387538D+00,     
     &  0.28975386476618070D+00,     
     & -0.28975386476618070D+00,     
     & -0.28975386476618070D+00,     
     &  0.28975386476618070D+00,     
     &  0.61367241226233160D+00,     
     & -0.61367241226233160D+00,     
     & -0.61367241226233160D+00,     
     &  0.61367241226233160D+00,     
     &  0.18378979287798017D+00,     
     & -0.18378979287798017D+00,     
     & -0.18378979287798017D+00,     
     &  0.18378979287798017D+00,     
     &  0.38925011625173161D+00,     
     & -0.38925011625173161D+00,     
     & -0.38925011625173161D+00,     
     &  0.38925011625173161D+00,     
     &  7.76896479525748113D-02, 
     & -7.76896479525748113D-02, 
     & -7.76896479525748113D-02, 
     &  7.76896479525748113D-02, 
     &  0.16453962988669860D+00,     
     & -0.16453962988669860D+00,     
     & -0.16453962988669860D+00,     
     &  0.16453962988669860D+00 /  
      data y / 
     & 0.0000000000000000D+00,      
     & 0.0000000000000000D+00,      
     & 0.88091731624450909D+00,      
     &-0.88091731624450909D+00,      
     & 0.0000000000000000D+00,      
     & 0.0000000000000000D+00,      
     & 0.70491874112648223D+00,      
     &-0.70491874112648223D+00,     
     & 0.0000000000000000D+00,      
     & 0.0000000000000000D+00,      
     & 0.44712732143189760D+00,      
     &-0.44712732143189760D+00,      
     & 0.0000000000000000D+00,      
     & 0.0000000000000000D+00,      
     & 0.18900486065123448D+00,      
     &-0.18900486065123448D+00,      
     & 0.36209733410322176D+00,      
     & 0.36209733410322176D+00,      
     &-0.36209733410322176D+00,      
     &-0.36209733410322176D+00,      
     & 0.76688932060387538D+00,      
     & 0.76688932060387538D+00,      
     &-0.76688932060387538D+00,      
     &-0.76688932060387538D+00,      
     & 0.28975386476618070D+00,      
     & 0.28975386476618070D+00,      
     &-0.28975386476618070D+00,      
     &-0.28975386476618070D+00,      
     & 0.61367241226233160D+00,      
     & 0.61367241226233160D+00,      
     &-0.61367241226233160D+00,      
     &-0.61367241226233160D+00,      
     & 0.18378979287798017D+00,      
     & 0.18378979287798017D+00,      
     &-0.18378979287798017D+00,      
     &-0.18378979287798017D+00,      
     & 0.38925011625173161D+00,      
     & 0.38925011625173161D+00,      
     &-0.38925011625173161D+00,      
     &-0.38925011625173161D+00,      
     & 7.76896479525748113D-02, 
     & 7.76896479525748113D-02, 
     &-7.76896479525748113D-02, 
     &-7.76896479525748113D-02, 
     & 0.16453962988669860D+00,      
     & 0.16453962988669860D+00,      
     &-0.16453962988669860D+00, 
     &-0.16453962988669860D+00 /
      data z / 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 4.85005494469969989D-02, 
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.23860073755186201D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.51704729510436798D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00,      
     & 0.79585141789677305D+00, 
     & 0.79585141789677305D+00 /     
c
c  Quadrature.
c
      quad = 0.0D+00
      do i = 1, order
        quad = quad + w(i) * func ( x(i), y(i), z(i) )
      end do
c
c  Volume.
c
      volume = pyramid_unit_volume_3d ( )
c
c  Result.
c
      result = quad * volume

      return
      end
      function pyramid_unit_monomial_3d ( alpha, beta, gamma )

c*****************************************************************************80
c
cc PYRAMID_UNIT_MONOMIAL_3D: monomial integral in a unit pyramid in 3D.
c
c  Discussion:
c
c    This routine returns the integral of X^ALPHA Y^BETA Z^GAMMA over
c    the unit pyramid.
c
c    The unit pyramid is defined as:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer ALPHA, BETA, GAMMA, the exponents of
c    X, Y and Z in the monomial.
c
c    Output, double precision PYRAMID_UNIT_MONOMIAL_3D, the volume of 
c    the pyramid.
c
      implicit none

      integer alpha
      integer beta
      integer gamma
      integer i
      integer i_hi
      double precision pyramid_unit_monomial_3d
      double precision r8_choose
      double precision r8_mop
      double precision value

      value = 0.0D+00

      if ( mod ( alpha, 2 ) .eq. 0 .and. mod ( beta, 2 ) .eq. 0 ) then    

        i_hi = 2 + alpha + beta

        do i = 0, i_hi
          value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) 
     &    / dble ( i + gamma + 1 )
        end do

        value = value 
     &        * 2.0D+00 / dble ( alpha + 1 ) 
     &        * 2.0D+00 / dble ( beta + 1 )

      end if

      pyramid_unit_monomial_3d = value

      return
      end
      function pyramid_unit_volume_3d ( )

c*****************************************************************************80
c
cc PYRAMID_UNIT_VOLUME_3D: volume of a unit pyramid with square base in 3D.
c
c  Integration region:
c
c    - ( 1 - Z ) <= X <= 1 - Z
c    - ( 1 - Z ) <= Y <= 1 - Z
c              0 <= Z <= 1.
c
c  Discussion:
c
c    The volume of this unit pyramid is 4/3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision PYRAMID_UNIT_VOLUME_3D, the volume of 
c    the pyramid.
c
      implicit none

      double precision pyramid_unit_volume_3d

      pyramid_unit_volume_3d = 4.0D+00 / 3.0D+00

      return
      end
      function pyramid_volume_3d ( r, h )

c*****************************************************************************80
c
cc PYRAMID_VOLUME_3D returns the volume of a pyramid with square base in 3D.
c
c  Integration region:
c
c    - ( H - Z ) * R <= X <= ( H - Z ) * R
c    - ( H - Z ) * R <= Y <= ( H - Z ) * R
c                  0 <= Z <= H.
c
c  Discussion:
c
c    A pyramid with square base can be regarded as the upper half of a
c    3D octahedron.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the "radius" of the pyramid, that is, half the
c    length of one of the sides of the square base.
c
c    Input, double precision H, the height of the pyramid.
c
c    Output, double precision PYRAMID_VOLUME_3D, the volume of the pyramid.
c
      implicit none

      double precision h
      double precision pyramid_volume_3d
      double precision r

      pyramid_volume_3d = ( 4.0D+00 / 3.0D+00 ) * h * r * r

      return
      end
      subroutine qmdpt ( func, n, nsub, result )

c*********************************************************************72
c
cc QMDPT carries out product midpoint quadrature for the unit cube in ND.
c
c  Integration region:
c
c    -1 <= X(1:N) <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates the function, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the cube.
c
c    Input, integer NSUB, the number of subdivisions 
c    (in each dimension).
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      integer ihi
      integer ix(n)
      integer j
      logical more
      integer nsub
      double precision quad
      double precision result
      double precision volume
      double precision w
      double precision x(n)

      w = 1.0D+00 / dble ( nsub ** n )
      quad = 0.0D+00

      more = .false.
      ihi = nsub ** n

      do i = 1, ihi

        call vec_lex_next ( n, nsub, ix, more )

        do j = 1, n
          x(j) = dble ( 2 * ix(j) + 1 - nsub ) / dble ( nsub )
        end do

        quad = quad + w * func ( n, x )

      end do

      volume = 2.0D+00 ** n
      result = quad * volume

      return
      end
      function qmult_1d ( func, a, b )

c*********************************************************************72
c
cc QMULT_1D approximates an integral over an interval in 1D.
c
c  Integration region:
c
c    A <= X <= B.
c
c  Discussion:
c
c    A 16 point 31-st degree Gauss-Legendre formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), of the form
c      function func ( x )
c      double precision func
c      double precision x
c
c    Input, double precision A, B, the lower and upper limits of integration.
c
c    Output, double precision QMULT_1D, the approximate integral of 
c    the function.
c
      implicit none

      integer order
      parameter ( order = 16 )

      double precision a
      double precision b
      double precision func
      external func
      integer i
      double precision quad
      double precision qmult_1d
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)

      call legendre_set ( order, xtab, weight )

      quad = 0.0D+00
      do i = 1, order
        x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( a + b )
        quad = quad + 0.5D+00 * weight(i) * func ( x )
      end do

      volume = b - a
      qmult_1d = quad * volume

      return
      end
      function qmult_2d ( func, a, b, fup, flo )

c*********************************************************************72
c
cc QMULT_2D approximates an integral with varying Y dimension in 2D.
c
c  Integration region:
c
c      A <= X <= B
c
c    and
c
c      FLO(X) <= Y <= FHI(X).
c
c  Discussion:
c
c    A 256 point product of two 16 point 31-st degree Gauss-Legendre
c    quadrature formulas is used.
c
c    This routine could easily be modified to use a different
c    order product rule by changing the value of ORDER.
c
c    Another easy change would allow the X and Y directions to
c    use quadrature rules of different orders.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y), of the form
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision A, B, the lower and upper limits of X integration.
c
c    Input, external FUP, FLO, the names of the user
c    supplied functions which evaluate the upper and lower
c    limits of the Y integration, of the form
c
c      function fup(x)
c      double precision fup
c      double precision x
c
c    and
c
c      function flo(x)
c      double precision flo
c      double precision x
c
c    Output, double precision QMULT_2D, the approximate integral of 
c    the function.
c
      implicit none

      integer order
      parameter ( order = 16 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision func
      external func
      double precision flo
      external flo
      double precision fup
      external fup
      integer i
      integer j
      double precision quad
      double precision qmult_2d
      double precision w1
      double precision w2
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y

      call legendre_set ( order, xtab, weight )

      quad = 0.0D+00

      do i = 1, order

        w1 = 0.5D+00 * ( b - a ) * weight(i)
        x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( b + a )
        c = flo ( x )
        d = fup ( x )

        do j = 1, order

          w2 = 0.5D+00 * ( d - c ) * weight(j)
          y = 0.5D+00 * ( d - c ) * xtab(j) + 0.5D+00 * ( d + c )
          quad = quad + w1 * w2 * func ( x, y )

        end do

      end do

      qmult_2d = quad

      return
      end
      function qmult_3d ( func, a, b, fup1, flo1, fup2, flo2 )

c*********************************************************************72
c
cc QMULT_3D approximates an integral with varying Y and Z dimension in 3D.
c
c  Integration region:
c
c      A         <= X <= B,
c    and
c      FLO(X)    <= Y <= FHI(X),
c    and
c      FLO2(X,Y) <= Z <= FHI2(X,Y).
c
c  Discussion:
c
c    A 4096 point product of three 16 point 31-st degree Gauss-Legendre
c    quadrature formulas is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision A, B, the lower and upper limits of X integration.
c
c    Input, external FUP1, FLO1, the names of the user
c    supplied functions which evaluate the upper and lower
c    limits of the Y integration, of the form
c
c      function fup1(x)
c      double precision fup1
c      double precision x
c
c    and
c
c      function flo1(x)
c      double precision flo1
c      double precision x
c
c    Input, external FUP2, FLO2, the names of the user
c    supplied functions which evaluate the upper and lower
c    limits of the Z integration, of the form
c
c      function fup2(x,y)
c      double precision fup2
c      double precision x
c      double precision y
c
c    and
c
c      function flo2(x,y)
c      double precision flo2
c      double precision x
c      double precision y
c
c    Output, double precision QMULT_3D, the approximate integral of 
c    the function.
c
      implicit none

      integer order
      parameter ( order = 16 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision f
      double precision func
      external func
      double precision flo1
      external flo1
      double precision flo2
      external flo2
      double precision fup1
      external fup1
      double precision fup2
      external fup2
      integer i
      integer j
      integer k
      double precision qmult_3d
      double precision quad
      double precision volume
      double precision w1
      double precision w2
      double precision w3
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y
      double precision z

      call legendre_set ( order, xtab, weight )

      quad = 0.0D+00

      do i = 1, order

        x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( b + a )
        w1 = 0.5D+00 * weight(i)
        c = flo1 ( x )
        d = fup1 ( x )

        do j = 1, order

          w2 = 0.5D+00 * ( d - c ) * weight(j)
          y = 0.5D+00 * ( d - c ) * xtab(j) + 0.5D+00 * ( d + c )
          e = flo2 ( x, y )
          f = fup2 ( x, y )

          do k = 1, order

            w3 = 0.5D+00 * ( f - e ) * weight(k)
            z = 0.5D+00 * ( f - e ) * xtab(k) + 0.5D+00 * ( f + e )
            quad = quad + w1 * w2 * w3 * func ( x, y, z )

          end do

        end do

      end do

      volume = b - a
      qmult_3d = quad * volume

      return
      end
      function r8_add ( x, y )

c*********************************************************************72
c
cc R8_ADD adds two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the numbers to be added.
c
c    Output, double precision R8_ADD, the sum of X and Y.
c
      implicit none

      double precision r8_add
      double precision x
      double precision y

      r8_add = x + y

      return
      end
      function r8_choose ( n, k )

c*********************************************************************72
c
cc R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R8 arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, double precision R8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      double precision r8_choose
      double precision value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0D+00

      else if ( mn .eq. 0 ) then

        value = 1.0D+00

      else

        mx = max ( k, n - k )
        value = dble ( mx + 1 )

        do i = 2, mn
          value = ( value * dble ( mx + i ) ) / dble ( i )
        end do

      end if

      r8_choose = value

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

      double precision one
      double precision r8_add
      double precision r8_epsilon
      double precision temp
      double precision test
      double precision value

      save value

      data value / 0.0D+00 /

      if ( value .ne. 0.0D+00 ) then
        r8_epsilon = value
        return
      end if

      one = dble ( 1 )

      value = one
      temp = value / 2.0D+00
      test = r8_add ( one, temp )

10    continue

      if ( one .lt. test ) then
        value = temp
        temp = value / 2.0D+00
        test = r8_add ( one, temp )
        go to 10
      end if

      r8_epsilon = value

      return
      end
      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

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
      function r8_gamma_log ( x )

c*********************************************************************72
c
cc R8_GAMMA_LOG evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
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
c    Output, double precision R8_GAMMA_LOG, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision eps
      double precision frtbig
      integer i
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision r8_gamma_log
      double precision res
      double precision sqrtpi
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
c
c  Mathematical constants
c
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( 0.0D+00 .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. 1.5D+00 ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = 0.0D+00
            xm1 = ( y - 0.5D+00 ) - 0.5D+00
          end if

          if ( y .le. 0.5D+00 .or. pnt68 .le. y ) then

            xden = 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - 0.5D+00 ) - 0.5D+00
            xden = 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do

            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. 4.0D+00 ) then

          xm2 = y - 2.0D+00
          xden = 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. 12.0D+00 ) then

          xm4 = y - 4.0D+00
          xden = - 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = 0.0D+00

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - 0.5D+00 * corr
          res = res + y * ( corr - 1.0D+00 )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      r8_gamma_log = res

      return
      end
      subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

c*********************************************************************72
c
cc R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
c
c  Discussion:
c
c    A minor bug was corrected.  The HW variable, used in several places as
c    the "old" value of a quantity being iteratively improved, was not
c    being initialized.  JVB, 11 February 2008.
c
c    The original version of this program allowed the input arguments to
c    be modified, although they were restored to their input values before exit.
c    This is unacceptable if the input arguments are allowed to be constants.
c    The code has been modified so that the input arguments are never modified.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
c    This FORTRAN77 version by John Burkardt.
c
c    The original FORTRAN77 version of this routine is copyrighted by
c    Shanjie Zhang and Jianming Jin.  However, they give permission to
c    incorporate this routine into a user program provided that the copyright
c    is acknowledged.
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45
c
c  Parameters:
c
c    Input, double precision A_INPUT, B_INPUT, C_INPUT, X_INPUT, 
c    the arguments of the function.  The user is allowed to pass these
c    values as constants or variables.
c    C_INPUT must not be equal to a nonpositive integer.
c    X_INPUT .lt. 1.
c
c    Output, double precision HF, the value of the function.
c
      implicit none

      double precision a
      double precision a_input
      double precision a0
      double precision aa
      double precision b
      double precision b_input
      double precision bb
      double precision c
      double precision c_input
      double precision c0
      double precision c1
      double precision el
      parameter ( el = 0.5772156649015329D+00 )
      double precision eps
      double precision f0
      double precision f1
      double precision g0
      double precision g1
      double precision g2
      double precision g3
      double precision ga
      double precision gabc
      double precision gam
      double precision gb
      double precision gbm
      double precision gc
      double precision gca
      double precision gcab
      double precision gcb
      double precision gm
      double precision hf
      double precision hw
      integer j
      integer k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      integer m
      integer nm
      double precision pa
      double precision pb
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r0
      double precision r1
      double precision r8_gamma
      double precision r8_psi
      double precision rm
      double precision rp
      double precision sm
      double precision sp
      double precision sp0
      double precision x
      double precision x_input
      double precision x1
c
c  Immediately copy the input argumentsc
c
      a = a_input
      b = b_input
      c = c_input
      x = x_input

      l0 = ( c .eq. aint ( c ) ) .and. ( c .lt. 0.0D+00 )
      l1 = ( 1.0D+00 - x .lt. 1.0D-15 ) .and. ( c - a - b .le. 0.0D+00 )
      l2 = ( a .eq. aint ( a ) ) .and. ( a .lt. 0.0D+00 )
      l3 = ( b .eq. aint ( b ) ) .and. ( b .lt. 0.0D+00 )
      l4 = ( c - a .eq. aint ( c - a ) ) .and. ( c - a .le. 0.0D+00 )
      l5 = ( c - b .eq. aint ( c - b ) ) .and. ( c - b .le. 0.0D+00 )

      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        return
      end if

      if ( 0.95D+00 .lt. x ) then
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if

      if ( x .eq. 0.0D+00 .or. a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then

        hf = 1.0D+00
        return

      else if ( 1.0D+00 - x .eq. eps .and. 0.0D+00 .lt. c - a - b ) then

        gc = r8_gamma ( c )
        gcab = r8_gamma ( c - a - b )
        gca = r8_gamma ( c - a )
        gcb = r8_gamma ( c - b )
        hf = gc * gcab / ( gca * gcb )
        return

      else if ( 1.0D+00 + x .le. eps .and. 
     &  abs ( c - a + b - 1.0D+00 ) .le. eps ) then

        g0 = sqrt ( pi ) * 2.0D+00**( - a )
        g1 = r8_gamma ( c )
        g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
        g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
        hf = g0 * g1 / ( g2 * g3 )
        return

      else if ( l2 .or. l3 ) then

        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if

        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if

        hf = 1.0D+00
        r = 1.0D+00

        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do

        return

      else if ( l4 .or. l5 ) then

        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if

        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if

        hf = 1.0D+00
        r  = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x )**( c - a - b ) * hf
        return

      end if

      aa = a
      bb = b
      x1 = x

      if ( x .lt. 0.0D+00 ) then
        x = x / ( x - 1.0D+00 )
        if ( a .lt. c .and. b .lt. a .and. 0.0D+00 .lt. b ) then
          a = bb
          b = aa
        end if
        b = c - b
      end if

      if ( 0.75D+00 .le. x ) then

        gm = 0.0D+00

        if ( abs ( c - a - b - aint ( c - a - b ) ) .lt. 1.0D-15 ) then

          m = int ( c - a - b )
          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gam = r8_gamma ( a + m )
          gbm = r8_gamma ( b + m )

          pa = r8_psi ( a )
          pb = r8_psi ( b )

          if ( m /= 0 ) then
            gm = 1.0D+00
          end if

          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do

          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do

          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00

          if ( 0 .le. m ) then

            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) 
     &          + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + ( 1.0D+00 - a ) 
     &            / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) 
     &            + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp + sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + m + k - 1.0D+00 ) 
     &          * ( b + m + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          else if ( m .lt. 0 ) then

            m = - m
            c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
            c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

            do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) 
     &          * ( b - m + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) 
     &          / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + 1.0D+00 / dble ( j + k )
              end do

              rp = pa + pb + 2.0D+00 * el + sp - sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          end if

        else

          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gca = r8_gamma ( c - a )
          gcb = r8_gamma ( c - b )
          gcab = r8_gamma ( c - a - b )
          gabc = r8_gamma ( a + b - c )
          c0 = gc * gcab / ( gca * gcb )
          c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
          hf = 0.0D+00
          hw = hf
          r0 = c0
          r1 = c1

          do k = 1, 250

            r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &        / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

            r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &        / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

            hf = hf + r0 + r1

            if ( abs ( hf - hw ) .lt. abs ( hf ) * eps ) then
              exit
            end if

            hw = hf

          end do

          hf = hf + c0 + c1

        end if

      else

        a0 = 1.0D+00

        if ( a .lt. c .and. c .lt. 2.0D+00 * a .and. 
     &       b .lt. c .and. c .lt. 2.0D+00 * b ) then

          a0 = ( 1.0D+00 - x )**( c - a - b )
          a = c - a
          b = c - b

        end if

        hf = 1.0D+00
        hw = hf
        r = 1.0D+00

        do k = 1, 250

          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x

          hf = hf + r

          if ( abs ( hf - hw ) .le. abs ( hf ) * eps ) then
            exit
          end if

          hw = hf

        end do

        hf = a0 * hf

      end if

      if ( x1 .lt. 0.0D+00 ) then
        x = x1
        c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
        hf = c0 * hf
      end if

      a = aa
      b = bb

      if ( 120 .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
        write ( *, '(a)' ) '  A large number of iterations were needed.'
        write ( *, '(a)' ) 
     &    '  The accuracy of the results should be checked.'
      end if

      return
      end
      function r8_mop ( i )

c*********************************************************************72
c
cc R8_MOP returns the I-th power of -1 as an R8.
c
c  Discussion:
c
c    An R8 is a double precision real value.
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
c    Input, integer I, the power of -1.
c
c    Output, double precision R8_MOP, the I-th power of -1.
c
      implicit none

      integer i
      double precision r8_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        r8_mop = + 1.0D+00
      else
        r8_mop = - 1.0D+00
      end if

      return
      end
      function r8_psi ( xx )

c*********************************************************************72
c
cc R8_PSI evaluates the function Psi(X).
c
c  Discussion:
c
c    This routine evaluates the logarithmic derivative of the
c    GAMMA function,
c
c      PSI(X) = d/dX (GAMMA(X)) / GAMMA(X) 
c             = d/dX LN ( GAMMA(X) )
c
c    for real X, where either
c
c      -XMAX1 < X < -XMIN  and X is not a negative integer), 
c
c    or
c
c      XMIN < X.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Anthony Strecok, Henry Thacher,
c    Chebyshev Approximations for the Psi Function,
c    Mathematics of Computation,
c    Volume 27, Number 121, January 1973, pages 123-127.
c
c  Parameters:
c
c    Input, double precision XX, the argument of the function.
c
c    Output, double precision R8_PSI, the value of the function.
c
      implicit none

      double precision aug
      double precision den
      integer i
      integer n
      integer nq
      double precision p1(9)
      double precision p2(7)
      double precision piov4
      double precision q1(8)
      double precision q2(6)
      double precision r8_psi
      double precision sgn
      double precision xlarge
      double precision upper
      double precision w
      double precision x
      double precision xinf
      double precision xmax1
      double precision xmin1
      double precision xsmall
      double precision x01
      double precision x01d
      double precision x02
      double precision xx
      double precision z
c
c  Mathematical constants.  PIOV4 = pi / 4
c
      data piov4 /7.8539816339744830962d-01/
c
c  Machine-dependent constants
c
      data xinf /1.70d+38/
      data xmin1 /5.89d-39/
      data xmax1 /3.60d+16/
      data xsmall /2.05d-09/
      data xlarge /2.04d+15/
c
c  Zero of psi(x)
c
      data x01 /187.0d0/
      data x01d /128.0d0/
      data x02 /6.9464496836234126266d-04/
c
c  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
c
      data p1/4.5104681245762934160d-03,5.4932855833000385356d+00,
     &        3.7646693175929276856d+02,7.9525490849151998065d+03,
     &        7.1451595818951933210d+04,3.0655976301987365674d+05,
     &        6.3606997788964458797d+05,5.8041312783537569993d+05,
     &        1.6585695029761022321d+05/
      data q1/9.6141654774222358525d+01,2.6287715790581193330d+03,
     &        2.9862497022250277920d+04,1.6206566091533671639d+05,
     &        4.3487880712768329037d+05,5.4256384537269993733d+05,
     &        2.4242185002017985252d+05,6.4155223783576225996d-08/
c
c  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
c  for 3.0 < x.
c
      data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01,
     &        -1.9784554148719218667d+01,-8.8100958828312219821d+00,
     &        -1.4479614616899842986d+00,-7.3689600332394549911d-02,
     &        -6.5135387732718171306d-21/
      data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02,
     &         2.4736979003315290057d+02, 1.0742543875702278326d+02,
     &         1.7463965060678569906d+01, 8.8427520398873480342d-01/

      x = xx
      w = abs ( x )
      aug = 0.0D+00
c
c  Check for valid arguments, then branch to appropriate algorithm.
c
      if ( - x .ge. xmax1 .or. w .lt. xmin1 ) then
        r8_psi = xinf
        if ( 0.0D+00 .lt. x ) then
          r8_psi = -xinf
        end if
        return
      end if

      if ( x .ge. 0.5D+00 ) then
        go to 200
c
c  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
c  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
c
      else if ( w .le. xsmall ) then
        aug = - 1.0D+00 / x
        go to 150
      end if
c
c  Argument reduction for cotangent.
c
  100 continue

      if ( x .lt. 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - aint ( w )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - dble ( nq ) * 0.25D+00 )
c
c  W is now related to the fractional part of 4.0 * X.
c  Adjust argument to correspond to values in the first
c  quadrant and determine the sign.
c
      n = nq / 2

      if ( n + n .ne. nq ) then
        w = 1.0D+00 - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) .ne. 0 ) then
        sgn = - sgn
      end if
c
c  Determine the final value for  -pi * cotan(pi*x).
c
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) .eq. 0 ) then
c
c  Check for singularity.
c
        if ( z .eq. 0.0D+00 ) then
          r8_psi = xinf
          if ( 0.0D+00 .lt. x ) then
            r8_psi = -xinf
          end if
          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else
        aug = sgn * ( 4.0D+00 * tan ( z ) )
      end if

  150 continue

      x = 1.0D+00 - x

  200 continue
c
c  0.5 <= X <= 3.0.
c
      if ( x .le. 3.0D+00 ) then

        den = x
        upper = p1(1) * x
        do i = 1, 7
          den = ( den + q1(i) ) * x
          upper = ( upper + p1(i+1) ) * x
        end do
        den = ( upper + p1(9) ) / ( den + q1(8) )
        x = ( x - x01 / x01d ) - x02
        r8_psi = den * x + aug
        return

      end if
c
c  3.0 < X.
c
      if ( x .lt. xlarge ) then
        w = 1.0D+00 / ( x * x )
        den = w
        upper = p2(1) * w
        do i = 1, 5
          den = ( den + q2(i) ) * w
          upper = ( upper + p2(i+1) ) * w
        end do
        aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
      end if

      r8_psi = aug + log ( x )

      return
      end
      subroutine r8_swap ( x, y )

c*********************************************************************72
c
cc R8_SWAP switches two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      double precision x
      double precision y
      double precision z

      z = x
      x = y
      y = z

      return
      end
      subroutine r8_swap3 ( x, y, z )

c*********************************************************************72
c
cc R8_SWAP3 swaps three R8's.
c
c  Example:
c
c    Input:
c
c      X = 1, Y = 2, Z = 3
c
c    Output:
c
c      X = 2, Y = 3, Z = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y, Z, three values to be swapped.
c
      implicit none

      double precision w
      double precision x
      double precision y
      double precision z

      w = x
      x = y
      y = z
      z = w

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
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
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

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8ge_det ( n, a_lu, pivot, det )

c*********************************************************************72
c
cc R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_LU(N,N), the LU factors from R8GE_FA or R8GE_TRF.
c
c    Input, integer PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision det
      integer i
      integer pivot(n)

      det = 1.0D+00

      do i = 1, n
        det = det * a_lu(i,i)
        if ( pivot(i) .ne. i ) then
          det = - det
        end if
      end do

      return
      end
      subroutine r8ge_fa ( n, a, pivot, info )

c*********************************************************************72
c
cc R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input/output, double precision A(N,N), the matrix to be factored.
c    On output, A contains an upper triangular matrix and the multipliers
c    which were used to obtain it.  The factorization can be written
c    A = L * U, where L is a product of permutation and unit lower
c    triangular matrices and U is upper triangular.
c
c    Output, integer PIVOT(N), a vector of pivot indices.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c 
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer info
      integer pivot(n)
      integer j
      integer k
      integer l
      double precision t

      info = 0

      do k = 1, n - 1
c
c  Find L, the index of the pivot row.
c
        l = k
        do i = k + 1, n
          if ( abs ( a(l,k) ) .lt. abs ( a(i,k) ) ) then
            l = i
          end if
        end do

        pivot(k) = l
c
c  If the pivot index is zero, the algorithm has failed.
c
        if ( a(l,k) .eq. 0.0D+00 ) then
          info = k
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Interchange rows L and K if necessary.
c
        if ( l .ne. k ) then
          t      = a(l,k)
          a(l,k) = a(k,k)
          a(k,k) = t
        end if
c
c  Normalize the values that lie below the pivot entry A(K,K).
c
        do i = k + 1, n
          a(i,k) = - a(i,k) / a(k,k)
        end do
c
c  Row elimination with column indexing.
c
        do j = k + 1, n

          if ( l .ne. k ) then
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end if

          do i = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do

        end do

      end do

      pivot(n) = n

      if ( a(n,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      end if

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
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
c    Input, character ( len = * ) TITLE, a title.
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
c    Input, character ( len = * ) TITLE, a title.
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
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

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

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
      subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

c*********************************************************************72
c
cc R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
c
c    Unless N = 1, X(1) = XLO and X(N) = XHI.
c
c    If N = 1, then X(1) = 0.5*(XLO+XHI).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision XLO, XHI, the low and high values.
c
c    Input, integer IVAL, the index of the desired point.
c    IVAL is normally between 1 and N, but may be any integer value.
c
c    Output, double precision XVAL, the IVAL-th of N evenly spaced values
c    between XLO and XHI.
c
      implicit none

      integer n

      integer ival
      double precision xhi
      double precision xlo
      double precision xval

      if ( n .eq. 1 ) then

        xval = 0.5D+00 * ( xlo + xhi )

      else

        xval = ( dble ( n - ival     ) * xlo
     &         + dble (     ival - 1 ) * xhi )
     &         / dble ( n        - 1 )

      end if

      return
      end
      subroutine r8vec_mirror_next ( n, a, done )

c*********************************************************************72
c
cc R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
c
c  Discussion:
c
c    In normal use, the user would set every element of A to be positive.
c    The routine will take the input value of A, and output a copy in
c    which the signs of one or more entries have been changed.  Repeatedly
c    calling the routine with the output from the previous call will generate
c    every distinct "variation" of A; that is, all possible sign variations.
c
c    When the output variable DONE is TRUE (or equal to 1), then the
c    output value of A_NEW is the last in the series.
c
c    Note that A may have some zero values.  The routine will essentially
c    ignore such entries; more exactly, it will not stupidly assume that -0
c    is a proper "variation" of 0c
c
c    Also, it is possible to call this routine with the signs of A set
c    in any way you like.  The routine will operate properly, but it
c    will nonethess terminate when it reaches the value of A in which
c    every nonzero entry has negative sign.
c
c
c    More efficient algorithms using the Gray code seem to require internal
c    memory in the routine, which is not one of MATLAB's strong points,
c    or the passing back and forth of a "memory array", or the use of
c    global variables, or unnatural demands on the user.  This form of
c    the routine is about as clean as I can make it.
c
c  Example:
c
c      Input         Output
c    ---------    --------------
c    A            A_NEW     DONE
c    ---------    --------  ----
c     1  2  3     -1  2  3  false
c    -1  2  3      1 -2  3  false
c     1 -2  3     -1 -2  3  false
c    -1 -2  3      1  2 -3  false
c     1  2 -3     -1  2 -3  false
c    -1  2 -3      1 -2 -3  false
c     1 -2 -3     -1 -2 -3  false
c    -1 -2 -3      1  2  3  true
c
c     1  0  3     -1  0  3  false
c    -1  0  3      1  0 -3  false
c     1  0 -3     -1  0 -3  false
c    -1  0 -3      1  0  3  true
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, double precision A(N), a vector of real numbers.  
c    On output, some signs have been changed.
c
c    Output, logical DONE, is TRUE if the input vector A was the last element
c    in the series (every entry was nonpositive); the output vector is reset 
c    so that all entries are nonnegative, but presumably the ride is overc
c
      implicit none

      integer n

      double precision a(n)
      logical done
      integer i
      integer positive
c
c  Seek the first strictly positive entry of A.
c
      positive = 0
      do i = 1, n
        if ( 0.0D+00 .lt. a(i) ) then
          positive = i
          go to 10
        end if
      end do

10    continue
c
c  If there is no strictly positive entry of A, there is no successor.
c
      if ( positive .eq. 0 ) then
        do i = 1, n
          a(i) = -a(i)
        end do
        done = .true.
        return
      end if
c
c  Otherwise, negate A up to the positive entry.
c
      do i = 1, positive
        a(i) = -a(i)
      end do
      done = .false.

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
      subroutine rectangle_3d ( func, a, b, result )

c*********************************************************************72
c
cc RECTANGLE_3D approximates an integral inside a rectangular block in 3D.
c
c  Integration region:
c
c      A(1) <= X <= B(1),
c    and
c      A(2) <= Y <= B(2),
c    and
c      A(3) <= Z <= B(3).
c
c  Discussion:
c
c    An 8 point third degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function which
c    evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision A(3), B(3), the lower and upper limits
c    for X, Y and Z.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision a(3)
      double precision b(3)
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision quad
      double precision result
      double precision sqr3
      double precision volume
      double precision w
      double precision x
      double precision y
      double precision z

      sqr3 = 1.0D+00 / sqrt ( 3.0D+00 )
      w = 1.0D+00 / 8.0D+00

      quad = 0.0D+00

      do i = 1, 2

        x = sqr3 * ( -1 )**i
        x = 0.5D+00 * ( ( 1.0D+00 - x ) * b(1) 
     &    + ( 1.0D+00 + x ) * a(1) )

        do j = 1, 2

          y = sqr3 * (  -1 )**j
          y = 0.5D+00 * ( ( 1.0D+00 - y ) * b(2) 
     &      + ( 1.0D+00 + y ) * a(2) )

          do k = 1, 2

            z = sqr3 * ( -1 )**k
            z = 0.5D+00 * ( ( 1.0D+00 - z ) * b(3) 
     &        + ( 1.0D+00 + z ) * a(3) )

            quad = quad + w * func ( x, y, z )

          end do

        end do

      end do

      volume = ( b(1) - a(1) ) * ( b(2) - a(2) ) * ( b(3) - a(3) )
      result = volume * quad

      return
      end
      subroutine rectangle_sub_2d ( func, xval, yval, nsub, order, xtab, 
     &  ytab, weight, result )

c*****************************************************************************80
c
cc RECTANGLE_SUB_2D carries out a composite quadrature over a rectangle in 2D.
c
c  Integration region:
c
c      XVAL(1) <= X <= XVAL(2),
c    and
c      YVAL(1) <= Y <= YVAL(2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the function to be
c    integrated.  The user must declare the name an EXTERNAL
c    parameter in the calling program, pass the name of the
c    function in FUNC, and write a function of the form
c      function func ( x, y )
c    which evaluates the function at the point (X,Y).
c
c    Input, double precision XVAL(2), the left and right X coordinates.
c
c    Input, double precision YVAL(2), the lower and upper Y coordinates.
c
c    Input, integer NSUB(2).
c    NSUB(1) is the number of subintervals to use in the X direction,
c    and NSUB(2) is the same thing for Y.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision a(2)
      double precision b(2)
      double precision func
      external func
      integer i
      integer j
      integer k
      integer nsub(2)
      double precision quad_sub
      double precision result
      double precision result_sub
      double precision volume
      double precision volume_sub
      double precision weight(order)
      double precision x
      double precision xhi
      double precision xlo
      double precision xtab(order)
      double precision xval(2)
      double precision y
      double precision yhi
      double precision ylo
      double precision ytab(order)
      double precision yval(2)

      a(1) = xval(1)
      a(2) = yval(1)
      b(1) = xval(2)
      b(2) = yval(2)

      do i = 1, 2
        if ( a(i) .eq. b(i) ) then
          result = 0.0D+00
          return
        end if
      end do

      do i = 1, 2
        if ( nsub(i) .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RECTANGLE_SUB_2D - Fatal error!'
          write ( *, '(a,i8)' ) 
     &      '  Nonpositive value of NSUB(I) = ', nsub(i)
          write ( *, '(a,i8)' ) '  for index I = ', i
          stop
        end if
      end do
c
c  Break up the X interval into NSUB(1) subintervals.
c
      volume = 0.0D+00
      result = 0.0D+00

      do i = 1, nsub(1)

        call r8vec_even_select ( nsub(1)+1, a(1), b(1), i, xlo )
        call r8vec_even_select ( nsub(1)+1, a(1), b(1), i+1, xhi )
c
c  Break up the Y interval into NSUB(2) subintervals.
c
        do j = 1, nsub(2)

          call r8vec_even_select ( nsub(2)+1, a(2), b(2), j,   ylo )
          call r8vec_even_select ( nsub(2)+1, a(2), b(2), j+1, yhi )

          quad_sub = 0.0D+00
          do k = 1, order

            x = xlo + 0.5D+00 * ( xtab(k) + 1.0D+00 ) * ( xhi - xlo )
            y = ylo + 0.5D+00 * ( ytab(k) + 1.0D+00 ) * ( yhi - ylo )

            quad_sub = quad_sub + weight(k) * func ( x, y ) / 4.0D+00

          end do

          volume_sub = ( xhi - xlo ) * ( yhi - ylo )
          result_sub = quad_sub * volume_sub

          volume = volume + volume_sub
          result = result + result_sub

        end do

      end do

      return
      end
      subroutine rule_adjust ( a, b, c, d, order, x, w )

c*********************************************************************72
c
cc RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
c
c  Discussion:
c
c    Most quadrature rules are defined on a special interval, like
c    [-1,1] or [0,1].  To integrate over an interval, the abscissas
c    and weights must be adjusted.  This can be done on the fly,
c    or by calling this routine.
c
c    If the weight function W(X) is not 1, then the W vector will
c    require further adjustment by the user.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the definition interval.
c
c    Input, double precision C, D, the endpoints of the integration interval.
c
c    Input, integer ORDER, the number of abscissas and weights.
c
c    Input/output, double precision X(ORDER), W(ORDER), the abscissas
c    and weights.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      double precision w(order)
      double precision x(order)

      do i = 1, order

        x(i) = ( ( b - x(i)     ) * c   
     &         + (     x(i) - a ) * d ) 
     &         / ( b              - a )

        w(i) = ( ( d - c ) / ( b - a ) ) * w(i)

      end do

      return
      end
      subroutine setsim ( n, v )

c*********************************************************************72
c
cc SETSIM defines a unit simplex.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      integer i
      integer j
      real ( kind = 8 ) v(n,n+1)

      do j = 1, n + 1
        do i = 1, n
          v(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        v(i,i+1) = 1.0D+00
      end do
 
      return
      end
      subroutine simplex_nd ( func, n, v, result )

c*********************************************************************72
c
cc SIMPLEX_ND approximates an integral inside a simplex in ND.
c
c  Discussion:
c
c    An N+1 point second degree formula is used.
c
c    The integration region is the simplex bounded by the origin and a 
c    convex combination of N points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input/output, double precision V(N,N+1).  On input, each of the
c    N+1 columns of V contains the N coordinates of one of the
c    "corners" of the simplex in entries 1 through N, with
c    the last column being left free.
c    On output, V has been overwritten in the process of
c    computing the volume of the simplex.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision c
      double precision func
      external func
      integer i
      integer j
      double precision quad
      double precision result
      double precision simplex_volume_nd
      double precision t
      double precision v(n,n+1)
      double precision volume
      double precision w
      double precision x(n)

      c = 1.0D+00 / sqrt ( dble ( n + 2 ) )
      w = 1.0D+00 / dble ( n + 1 )

      do j = 1, n
        t = 0.0D+00
        do i = 1, n + 1
          t = t + v(i,j)
        end do
        x(j) = w * ( 1.0D+00 - c ) * t
      end do

      quad = 0.0D+00

      do j = 1, n + 1

        do i = 1, n
          x(i) = x(i) + c * v(i,j)
        end do

        quad = quad + w * func ( n, x )

        do i = 1, n
          x(i) = x(i) - c * v(i,j)
        end do

      end do

      volume = simplex_volume_nd ( n, v )
      result = quad * volume

      return
      end
      subroutine simplex_unit_01_nd ( func, n, result )

c*********************************************************************72
c
cc SIMPLEX_UNIT_01_ND approximates an integral inside the unit simplex in ND.
c
c  Integration region:
c
c      0 <= X(1:N),
c    and
c      sum ( X(1:N) ) <= 1.
c
c  Discussion:
c
c    An 1 point formula of degree 1 is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.  
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision coef
      parameter ( coef = 1.0D+00 )
      double precision func
      external func
      double precision quad
      double precision result
      double precision simplex_unit_volume_nd
      double precision volume
      double precision x(n)

      quad = 0.0D+00

      x(1:n) = 1.0D+00 / dble ( n )
      quad = quad + coef * func ( n, x )

      volume = simplex_unit_volume_nd ( n )

      result = quad * volume

      return
      end
      subroutine simplex_unit_03_nd ( func, n, result )

c*********************************************************************72
c
cc SIMPLEX_UNIT_03_ND approximates an integral inside the unit simplex in ND.
c
c  Integration region:
c
c      0 <= X(1:N),
c    and
c      sum ( X(1:N) ) <= 1.
c
c  Discussion:
c
c    An N+2 point formula of degree 3 is used.  This is Stroud TN:3-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.  
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision coef
      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision simplex_unit_volume_nd
      double precision volume
      double precision x(n)

      quad = 0.0D+00

      do i = 1, n
        x(i) = 1.0D+00 / dble ( n + 1 )
      end do
      coef = -0.25D+00 * dble ( ( n + 1 ) * ( n + 1 ) ) 
     &     / dble ( n + 2 )
      quad = quad + coef * func ( n, x )

      a = 1.0D+00 / dble ( n + 3 )
      b = 3.0D+00 / dble ( n + 3 )

      do i = 1, n
        x(i) = a
      end do
      coef = 0.25D+00 * dble ( ( n + 3 ) * ( n + 3 ) ) 
     &  / dble ( ( n + 1 ) * ( n + 2 ) )
      quad = quad + coef * func ( n, x )

      do i = 1, n

        x(i) = b
        quad = quad + coef * func ( n, x )
        x(i) = a

      end do

      volume = simplex_unit_volume_nd ( n )

      result = quad * volume

      return
      end
      subroutine simplex_unit_05_nd ( func, n, result )

c*********************************************************************72
c
cc SIMPLEX_UNIT_05_ND approximates an integral inside the unit simplex in ND.
c
c  Integration region:
c
c      0 <= X(1:N),
c    and
c      sum ( X(1:N) ) <= 1.
c
c  Discussion:
c
c    An N^2 + 3 N + 3 point formula of degree 5 is used.  This is
c    Stroud formula TN:5-1.
c
c    (For N = 2, the number of points is actually only 7, and
c     for N = 3, the number of points is actually only 15.)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    A Fifth Degree Integration Formula for the N-Simplex,
c    SIAM Journal on Numerical Analysis,
c    Volume 6, Number 1, March 1969.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.  For this 
c    routine, it must be the case that 2 <= N <= 16.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision coef1(2:16)
      double precision coef21(2:16)
      double precision coef22(2:16)
      double precision coef31(2:16)
      double precision coef32(2:16)
      double precision func
      external func
      integer i
      integer j
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision s1
      double precision s2
      double precision simplex_unit_volume_nd
      double precision u1
      double precision u2
      double precision v1
      double precision v2
      double precision volume
      double precision x(n)

      save coef1
      save coef21
      save coef22
      save coef31
      save coef32

      data coef1 /
     &  0.225D+00, 
     &  0.118518518519D+00, 
     &  0.0631521898883D+00, 
     &  0.235714285714D+00, 
     &  0.791575476992D+00, 
     &  1.85798728021D+00, 
     &  3.53666958042D+00, 
     &  5.90844340844D+00, 
     &  9.03765432098D+00, 
     &  12.9758241758D+00, 
     &  17.7645108738D+00, 
     &  23.4375030259D+00, 
     &  30.0224941950D+00, 
     &  37.5423613501D+00, 
     &  46.0161454949D+00 /
      data coef21 /
     &  0.12593918054483D+00, 
     &  0.0719370837790D+00, 
     &  0.0470456145702D+00, 
     &  0.0333009774677D+00, 
     &  0.0248633014592D+00, 
     &  0.0192679696358D+00, 
     &  0.0153322153879D+00, 
     &  0.0124316229901D+00, 
     &  0.0102112988361D+00, 
     &  0.00845730697460D+00, 
     &  0.00703433430999D+00, 
     &  0.00585330520067D+00, 
     &  0.00485356735291D+00, 
     &  0.00399261092720D+00, 
     &  0.00323988713017D+00 /
      data coef22 /
     &  0.13239415278851D+00, 
     &  0.0690682072263D+00, 
     &  0.0371530185868D+00, 
     & -0.0719253160920D+00, 
     & -0.264323879461D+00, 
     & -0.537926779961D+00, 
     & -0.886895605701D+00, 
     & -1.30409181465D+00, 
     & -1.78227048964D+00, 
     & -2.31462336314D+00, 
     & -2.89499045158D+00, 
     & -3.51790849765D+00, 
     & -4.17858310668D+00, 
     & -4.87282884913D+00, 
     & -5.59699944261D+00 /
      data coef31 /
     &  0.0D+00, 
     &  0.0529100529100D+00, 
     &  0.0261368740713D+00, 
     &  0.0499020181331D+00, 
     &  0.0782233395867D+00, 
     &  0.109041040862D+00, 
     &  0.140874828568D+00, 
     &  0.172735353396D+00, 
     &  0.203992490408D+00, 
     &  0.234263814181D+00, 
     &  0.263332763315D+00, 
     &  0.291091849264D+00, 
     &  0.317504208212D+00, 
     &  0.342577872069D+00, 
     &  0.366348654344D+00 /
      data coef32 / 
     &  0.0D+00, 
     &  0.0D+00, 
     &  0.0254485903613D+00, 
     &  0.0165000982690D+00, 
     &  0.0115218303668D+00, 
     &  0.00850478779483D+00, 
     &  0.00655297510968D+00, 
     &  0.00522372456259D+00, 
     &  0.00428017828134D+00, 
     &  0.00358722367033D+00, 
     &  0.00306362964360D+00, 
     &  0.00265836687133D+00, 
     &  0.00233816221525D+00, 
     &  0.00208061510846D+00, 
     &  0.00187022027571D+00 /

      if ( n .lt. 2 .or. 16 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SIMPLEX_UNIT_05_ND - Fatal error!'
        write ( *, '(a)' ) '  Input spatial dimension N out of range.'
        write ( *, '(a,i8)' ) '  N = ', n
        result = 0.0D+00
        return
      end if

      quad = 0.0D+00
c
c  S1
c
      do i = 1, n
        x(i) = 1.0D+00 / dble ( n + 1 )
      end do
      quad = quad + coef1(n) * func ( n, x )
c
c  S21
c
      r1 = ( dble ( n + 4 ) - sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 8 * n + 1 )
      s1 = 1.0D+00 - dble ( n ) * r1

      do i = 1, n
        x(i) = r1
      end do

      do i = 1, n + 1

        quad = quad + coef21(n) * func ( n, x )

        if ( 1 .lt. i ) then
          x(i-1) = r1
        end if

        if ( i .lt. n + 1 ) then
          x(i) = s1
        end if

      end do
c
c  S22
c
      r2 = ( dble ( n + 4 ) + sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 8 * n + 1 )
      s2 = 1.0D+00 - dble ( n ) * r2

      do i = 1, n
        x(i) = r2
      end do

      do i = 1, n + 1

        quad = quad + coef22(n) * func ( n, x )

        if ( 1 .lt. i ) then
          x(i-1) = r2
        end if

        if ( i .lt. n + 1 ) then
          x(i) = s2
        end if

      end do
c
c  S31
c
      u1 = ( dble ( n + 7 ) + 2.0D+00 * sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 14 * n - 11 )
      v1 = ( dble ( 4 * n - 2 ) 
     &  - dble ( n - 1 ) * sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 14 * n - 11 )

      do i = 1, n

        do j = 1, n
          x(j) = u1
        end do
        x(i) = v1

        do j = i, n

          if ( i .lt. j - 1 ) then
            x(j-1) = u1
          end if

          x(j) = v1

          quad = quad + coef31(n) * func ( n, x )

        end do

      end do
c
c  S32
c
      u2 = ( dble ( n + 7 ) - 2.0D+00 * sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 14 * n - 11 )
      v2 = ( dble ( 4 * n - 2 ) 
     &  + dble ( n - 1 ) * sqrt ( 15.0D+00 ) ) 
     &  / dble ( n * n + 14 * n - 11 )

      do i = 1, n

        do j = 1, n
          x(j) = u2
        end do
        x(i) = v2

        do j = i, n

          if ( i .lt. j - 1 ) then
            x(j-1) = u2
          end if

          x(j) = v2

          quad = quad + coef32(n) * func ( n, x )

        end do

      end do

      volume = simplex_unit_volume_nd ( n )

      result = quad * volume

      return
      end
      subroutine simplex_unit_05_2_nd ( func, n, result )

c*********************************************************************72
c
cc SIMPLEX_UNIT_05_2_ND approximates an integral inside the unit simplex in ND.
c
c  Integration region:
c
c      0 <= X(1:N),
c    and
c      sum ( X(1:N) ) <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision coef
      double precision func
      external func
      integer i
      integer j
      double precision quad
      double precision result
      double precision simplex_unit_volume_nd
      double precision volume
      double precision x(n)

      quad = 0.0D+00
c
c  Group 1
c
      do i = 1, n
        x(i) = 1.0D+00 / dble ( n + 1 )
      end do
      coef = dble ( ( n + 1 )**4 ) 
     &  / dble ( 32 * ( n + 2 ) * ( n + 3 ) )
      quad = quad + coef * func ( n, x )
c
c  Group 2
c
      a = 1.0D+00 / dble ( n + 3 )
      b = 3.0D+00 / dble ( n + 3 )

      do i = 1, n
        x(i) = a
      end do
      coef = - dble ( ( n + 3 )**4 ) 
     &  / dble ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 4 ) )
      quad = quad + coef * func ( n, x )

      do i = 1, n

        x(i) = b
        quad = quad + coef * func ( n, x )
        x(i) = a

      end do
c
c  Group 3
c
      a = 1.0D+00 / dble ( n + 5 )
      b = 5.0D+00 / dble ( n + 5 )

      do i = 1, n
        x(i) = a
      end do
      coef = dble ( ( n + 5 )**4 ) 
     &  / dble ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) )
      quad = quad + coef * func ( n, x )

      do i = 1, n

        x(i) = b
        quad = quad + coef * func ( n, x )
        x(i) = a

      end do
c
c  Group 4
c
      a = 1.0D+00 / dble ( n + 5 )
      b = 3.0D+00 / dble ( n + 5 )

      coef = dble ( ( n + 5 )**4 ) 
     &  / dble ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) )

      do i = 1, n

        do j = 1, n
          x(j) = a
        end do
        x(i) = b
        quad = quad + coef * func ( n, x )

        do j = i + 1, n
          x(j) = b
          quad = quad + coef * func ( n, x )
          x(j) = a
        end do

      end do

      volume = simplex_unit_volume_nd ( n )

      result = quad * volume

      return
      end
      function simplex_unit_volume_nd ( n )

c*********************************************************************72
c
cc SIMPLEX_UNIT_VOLUME_ND returns the volume of the unit simplex in ND.
c
c  Integration region:
c
c      0 <= X(1:N),
c    and
c      sum ( X(1:N) ) <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision SIMPLEX_UNIT_VOLUME_ND, the volume of the
c    unit simplex.
c
      implicit none

      integer n
      integer i4_factorial
      double precision simplex_unit_volume_nd

      simplex_unit_volume_nd = 1.0D+00 / dble ( i4_factorial ( n ) )

      return
      end
      function simplex_volume_nd ( n, v )

c*********************************************************************72
c
cc SIMPLEX_VOLUME_ND returns the volume of a simplex in ND.
c
c  Integration region:
c
c    The simplex bounded by the origin and a convex combination of N points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision V(N,N+1), the coordinates of the
c    vertices.
c
c    Output, double precision SIMPLEX_VOLUME_ND, the volume of 
c    the unit simplex.
c
      implicit none

      integer n

      double precision det
      integer i
      integer info
      integer j
      integer pivot(n)
      double precision simplex_unit_volume_nd
      double precision simplex_volume_nd
      double precision v(n,n+1)
      double precision volume
      double precision w(n,n)

      do i = 1, n
        do j = 1, n
          w(i,j) = v(i,j+1) - v(i,1)
        end do
      end do

      call r8ge_fa ( n, w, pivot, info )

      call r8ge_det ( n, w, pivot, det )
c
c  Multiply by the volume of the unit simplex, which serves as a
c  conversion factor between a parallelipiped and the simplex.
c
      simplex_volume_nd = abs ( det ) * simplex_unit_volume_nd ( n )

      return
      end
      function sin_power_int ( a, b, n )

c*********************************************************************72
c
cc SIN_POWER_INT evaluates the sine power integral.
c
c  Discussion:
c
c    The function is defined by
c
c      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
c
c    The algorithm uses the following fact:
c
c      Integral sin^n ( t ) = (1/n) * (
c        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer N, the power of the sine function.
c
c    Output, double precision SIN_POWER_INT, the value of the integral.
c
      implicit none

      double precision a
      double precision b
      double precision ca
      double precision cb
      integer m
      integer mlo
      integer n
      double precision sa
      double precision sb
      double precision sin_power_int
      double precision value

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
        write ( *, '(a)' ) '  Power N < 0.'
        value = 0.0
        stop
      end if

      sa = sin ( a )
      sb = sin ( b )
      ca = cos ( a )
      cb = cos ( b )

      if ( mod ( n, 2 ) .eq. 0 ) then
        value = b - a
        mlo = 2
      else
        value = ca - cb
        mlo = 3
      end if

      do m = mlo, n, 2
        value = ( dble ( m - 1 ) * value 
     &            + sa**(m-1) * ca - sb**(m-1) * cb ) 
     &    / dble ( m )
      end do

      sin_power_int = value

      return
      end
      subroutine sphere_05_nd ( func, n, center, r, result )

c*********************************************************************72
c
cc SPHERE_05_ND approximates an integral on the surface of a sphere in ND.
c
c  Integration region:
c
c    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2
c
c  Discussion:
c
c    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision CENTER(N), the center of the sphere.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision center(n)
      double precision func
      external func
      integer i
      integer iadd
      integer ihi
      integer ix(n)
      logical more
      integer ncard
      double precision quad
      double precision r
      double precision result
      double precision sphere_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 1.0D+00
      x2 = 1.0D+00 / sqrt ( dble ( n ) )

      w1 = 1.0D+00 / dble ( n * ( n + 2 ) )
      w2 = dble ( n ) / dble ( ( n + 2 ) * 2**n )

      do i = 1, n
        x(i) = center(i)
      end do

      quad = 0.0D+00

      do i = 1, n
        x(i) = center(i) + r * x1
        quad = quad + w1 * func ( n, x )
        x(i) = center(i) - r * x1
        quad = quad + w1 * func ( n, x )
        x(i) = center(i)
      end do

      more = .false.
      ihi = 2**n

      do i = 1, n
        x(i) = center(i) - r * x2
      end do

      do i = 1, ihi

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = center(iadd) - ( x(iadd) - center(iadd) )
        end if

        quad = quad + w2 * func ( n, x )

      end do

      volume = sphere_area_nd ( n, r )
      result = quad * volume

      return
      end
      subroutine sphere_07_1_nd ( func, n, center, r, result )

c*********************************************************************72
c
cc SPHERE_07_1_ND approximates an integral on the surface of a sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N) - CENTER(1:N) )^2 = R * R.
c
c  Discussion:
c
c    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision CENTER(N), the center of the sphere.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision center(n)
      double precision func
      external func
      integer i
      integer iadd
      integer ix(n)
      integer j
      integer jhi
      logical more
      integer ncard
      double precision quad
      double precision r
      double precision result
      double precision sphere_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision w3
      double precision x(n)
      double precision x1
      double precision x2
      double precision x3

      do i = 1, n
        x(i) = center(i)
      end do

      w1 = dble ( 8 - n ) / dble ( n * ( n + 2 ) * ( n + 4 ) )
      w2 = dble ( n**3 ) 
     &  / dble ( 2**n * n * ( n + 2 ) * ( n + 4 ) )
      w3 = 4.0D+00 / dble ( n * ( n + 2 ) * ( n + 4 ) )

      x1 = 1.0D+00
      x2 = 1.0D+00 / sqrt ( dble ( n ) )
      x3 = 1.0D+00 / sqrt ( 2.0D+00 )

      quad = 0.0D+00
c
c  First term.
c
      do i = 1, n
        x(i) = center(i) + r * x1
        quad = quad + w1 * func ( n, x )
        x(i) = center(i) - r * x1
        quad = quad + w1 * func ( n, x )
        x(i) = center(i)
      end do
c
c  Second term.
c
      do i = 1, n
        x(i) = center(i) - r * x2
      end do

      more = .false.
      jhi = 2**n

      do j = 1, jhi

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = center(iadd) - ( x(iadd) - center(iadd) )
        end if

        quad = quad + w2 * func ( n, x )

      end do
c
c  Third term.
c
      do i = 1, n
        x(i) = center(i)
      end do

      do i = 1, n - 1
        do j = i + 1, n
          x(i) = center(i) + r * x3
          x(j) = center(j) + r * x3
          quad = quad + w3 * func ( n, x )
          x(i) = center(i) - r * x3
          x(j) = center(j) + r * x3
          quad = quad + w3 * func ( n, x )
          x(i) = center(i) + r * x3
          x(j) = center(j) - r * x3
          quad = quad + w3 * func ( n, x )
          x(i) = center(i) - r * x3
          x(j) = center(j) - r * x3
          quad = quad + w3 * func ( n, x )
          x(i) = center(i)
          x(j) = center(j)
        end do
      end do

      volume = sphere_area_nd ( n, r )
      result = quad * volume

      return
      end
      subroutine sphere_area_3d ( r, area )

c*********************************************************************72
c
cc SPHERE_AREA_3D computes the surface area of an implicit sphere in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision AREA, the area of the sphere.
c
      implicit none

      double precision area
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      area = 4.0D+00 * pi * r * r

      return
      end
      subroutine sphere_area_nd ( dim_num, r, area )

c*********************************************************************72
c
cc SPHERE_AREA_ND computes the surface area of an implicit sphere in ND.
c
c  Discussion:
c
c    An implicit sphere in ND satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c    DIM_NUM   Area
c
c    2      2       * PI   * R
c    3      4       * PI   * R^2
c    4      2       * PI^2 * R^3
c    5      (8/3)   * PI^2 * R^4
c    6                PI^3 * R^5
c    7      (16/15) * PI^3 * R^6
c
c    Sphere_Area ( DIM_NUM, R ) =
c      2 * PI^(DIM_NUM/2) * R^(DIM_NUM-1) / Gamma ( DIM_NUM / 2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision AREA, the area of the sphere.
c
      implicit none

      double precision area
      integer dim_num
      double precision r
      double precision sphere_unit_area_nd

      area = r**( dim_num -1  ) * sphere_unit_area_nd ( dim_num )

      return
      end
      subroutine sphere_cap_area_2d ( r, h, area )

c*********************************************************************72
c
cc SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
c
c  Discussion:
c
c    Draw any radius of the sphere and note the point P where the radius
c    intersects the sphere.  Consider the point on the radius line which is
c    H units from P.  Draw the circle that lies in the plane perpendicular to
c    the radius, and which intersects the sphere.  The circle divides the sphere
c    into two pieces, and the corresponding disk divides the solid sphere into
c    two pieces.  The spherical cap is the part of the solid sphere that
c    includes the point P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "height" of the spherical cap.
c    H must be between 0 and 2 * R.
c
c    Output, double precision AREA, the area of the spherical cap.
c
      implicit none

      double precision arc_sine
      double precision area
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      if ( h .le. 0.0D+00 ) then
        area = 0.0D+00
      else if ( 2.0D+00 * r .le. h ) then
        area = 2.0D+00 * pi * r
      else

        theta = 2.0D+00
     &    * arc_sine ( sqrt ( r * r - ( r - h )**2 ) / r )
        area = r * theta

        if ( r .le. h ) then
          area = 2.0D+00 * pi * r - area
        end if

      end if

      return
      end
      subroutine sphere_cap_area_3d ( r, h, area )

c*********************************************************************72
c
cc SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
c
c  Discussion:
c
c    Draw any radius of the sphere and note the point P where the radius
c    intersects the sphere.  Consider the point on the radius line which is
c    H units from P.  Draw the circle that lies in the plane perpendicular to
c    the radius, and which intersects the sphere.  The circle divides the sphere
c    into two pieces, and the corresponding disk divides the solid sphere into
c    two pieces.  The spherical cap is the part of the solid sphere that
c    includes the point P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "height" of the spherical cap.
c    H must be between 0 and 2 * R.
c
c    Output, double precision AREA, the area of the spherical cap.
c
      implicit none

      double precision area
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      if ( h .le. 0.0D+00 ) then
        area = 0.0D+00
      else if ( 2.0D+00 * r .le. h ) then
        area = 4.0D+00 * pi * r * r
      else
        area = 2.0D+00 * pi * r * h
      end if

      return
      end
      subroutine sphere_cap_area_nd ( dim_num, r, h, area )

c*********************************************************************72
c
cc SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
c
c  Discussion:
c
c    The spherical cap is a portion of the surface of the sphere:
c
c      sum ( X(1:N)^2 ) = R^2
c
c    which is no more than H units from the uppermost point on the sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Ericson, Victor Zinoviev,
c    Codes on Euclidean Spheres,
c    Elsevier, 2001, pages 439-441.
c    QA166.7 E75
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "thickness" of the spherical cap,
c    which is normally between 0 and 2 * R.
c
c    Output, double precision AREA, the area of the spherical cap.
c
      implicit none

      double precision arc_sine
      double precision area
      double precision area2
      double precision h
      double precision haver_sine
      integer i
      integer dim_num
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision sphere_k
      double precision theta
      double precision ti
      double precision tj
      double precision tk

      if ( h .le. 0.0D+00 ) then
        area = 0.0D+00
        return
      end if

      if ( 2.0D+00 * r .le. h ) then
        call sphere_area_nd ( dim_num, r, area )
        return
      end if
c
c  For cases where R < H < 2 * R, work with the complementary region.
c
      haver_sine = sqrt ( ( 2.0D+00 * r - h ) * h )

      theta = arc_sine ( haver_sine / r )

      if ( dim_num .lt. 1 ) then

        area = -1.0D+00
        return

      else if ( dim_num .eq. 1 ) then

        area = 0.0D+00

      else if ( dim_num .eq. 2 ) then

        area = 2.0D+00 * theta * r

      else

        ti = theta

        tj = ti
        ti = 1.0D+00 - cos ( theta )

        do i = 2, dim_num - 2
          tk = tj
          tj = ti
          ti = ( dble ( i - 1 ) * tk
     &      - cos ( theta ) * sin ( theta )**( i - 1 ) )
     &      / dble ( i )
        end do

        area = sphere_k ( dim_num-1 ) * ti * r**( dim_num - 1 )

      end if
c
c  Adjust for cases where R < H < 2R.
c
      if ( r .lt. h ) then
        call sphere_area_nd ( dim_num, r, area2 )
        area = area2 - area
      end if

      return
      end
      subroutine sphere_cap_volume_2d ( r, h, volume )

c*********************************************************************72
c
cc SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
c
c  Discussion:
c
c    Draw any radius R of the circle and denote as P the point where the
c    radius intersects the circle.  Now consider the point Q which lies
c    on the radius and which is H units from P.  The line which is
c    perpendicular to the radius R and passes through Q divides the
c    circle into two pieces.  The piece including the point P is the
c    spherical (circular) cap of height (or thickness) H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "height" of the spherical cap.  H must
c    be between 0 and 2 * R.
c
c    Output, double precision VOLUME, the volume (area) of the spherical cap.
c
      implicit none

      double precision arc_sine
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta
      double precision volume

      if ( h .le. 0.0D+00 ) then

        volume = 0.0D+00

      else if ( 2.0D+00 * r .le. h ) then

        volume = pi * r * r

      else

        theta = 2.0D+00 * arc_sine ( sqrt ( r * r - ( r - h )**2 ) / r )
        volume = r * r * ( theta - sin ( theta ) ) / 2.0D+00

        if ( r .lt. h ) then
          volume = pi * r * r - volume
        end if

      end if

      return
      end
      subroutine sphere_cap_volume_3d ( r, h, volume )

c*********************************************************************72
c
cc SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
c
c  Discussion:
c
c    Draw any radius of the sphere and note the point P where the radius
c    intersects the sphere.  Consider the point on the radius line which is
c    H units from P.  Draw the circle that lies in the plane perpendicular to
c    the radius, and which intersects the sphere.  The circle divides the sphere
c    into two pieces, and the corresponding disk divides the solid sphere into
c    two pieces.  The part of the solid sphere that includes the point P
c    is the spherical cap of height (or thickness) H.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "height" of the spherical cap.  H must
c    be between 0 and 2 * R.
c
c    Output, double precision VOLUME, the volume of the spherical cap.
c
      implicit none

      double precision h
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision r
      double precision volume

      if ( h .le. 0.0D+00 ) then
        volume = 0.0D+00
      else if ( 2.0D+00 * r .le. h ) then
        volume = ( 4.0D+00 / 3.0D+00 ) * pi * r * r * r
      else
        volume = ( 1.0D+00 / 3.0D+00 ) * pi * h * h
     &    * ( 3.0D+00 * r - h )
      end if

      return
      end
      subroutine sphere_cap_volume_nd ( dim_num, r, h, volume )

c*********************************************************************72
c
cc SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
c
c  Discussion:
c
c    The spherical cap is a portion of the surface and interior of the sphere:
c
c      sum ( X(1:N)^2 ) .le. R^2
c
c    which is no more than H units from some point P on the sphere.
c
c
c    The algorithm proceeds from the observation that the N-dimensional
c    sphere can be parameterized by a quantity RC that runs along the
c    radius from the center to the point P.  The value of RC at the
c    base of the spherical cap is (R-H) and at P it is R.  We intend to
c    use RC as our integration parameeter.
c
c    The volume of the spherical cap is then the integral, as RC goes
c    from (R-H) to R, of the N-1 dimensional volume of the sphere
c    of radius RS, where RC^2 + RS^2 = R^2.
c
c    The volume of the N-1 dimensional sphere of radius RS is simply
c    some constants times RS^(N-1).
c
c    After factoring out the constant terms, and writing RC = R * cos ( T ),
c    and RS = R * sin ( T ), and letting
c      T_MAX = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) ),
c    the "interesting part" of our integral becomes
c
c      constants * R^N * Integral ( T = 0 to T_MAX ) sin^N ( T ) dT
c
c    The integral of sin^N ( T ) dT can be handled by recursion.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H, the "thickness" of the spherical cap,
c    which is normally between 0 and 2 * R.
c
c    Output, double precision VOLUME, the volume of the spherical cap.
c
      implicit none

      double precision angle
      double precision arc_sine
      double precision factor1
      double precision factor2
      double precision h
      integer dim_num
      double precision r
      double precision sin_power_int
      double precision sphere_unit_volume_nd
      double precision volume
      double precision volume2

      if ( h .le. 0.0D+00 ) then
        volume = 0.0D+00
        return
      end if

      if ( 2.0D+00 * r .le. h ) then
        call sphere_volume_nd ( dim_num, r, volume )
        return
      end if

      if ( dim_num .lt. 1 ) then

        volume = - 1.0D+00

      else if ( dim_num .eq. 1 ) then

        volume = h

      else

        factor1 = sphere_unit_volume_nd ( dim_num - 1 )

        angle = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) )

        factor2 = sin_power_int ( 0.0D+00, angle, dim_num )

        volume = factor1 * factor2 * r**dim_num

        if ( r .lt. h ) then
          call sphere_volume_nd ( dim_num, r, volume2 )
          volume = volume2 - volume
        end if

      end if

      return
      end
      function sphere_k ( dim_num )

c*********************************************************************72
c
cc SPHERE_K computes a factor useful for spherical computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Ericson, Victor Zinoviev,
c    Codes on Euclidean Spheres,
c    Elsevier, 2001, pages 439-441.
c    QA166.7 E75
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Output, double precision SPHERE_K, the factor.
c
      implicit none

      integer i4_factorial2
      integer dim_num
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sphere_k

      if ( mod ( dim_num, 2 ) .eq. 0 ) then
        sphere_k = ( 2.0D+00 * pi )**( dim_num / 2 )
      else
        sphere_k = 2.0D+00 * ( 2.0D+00 * pi )**( ( dim_num - 1 ) / 2 )
      end if

      sphere_k = sphere_k / dble ( i4_factorial2 ( dim_num - 2 ) )

      return
      end
      subroutine sphere_monomial_int_nd ( n, r, e, integral )

c*********************************************************************72
c
cc SPHERE_MONOMIAL_INT_ND integrates a monomial on surface of a sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = R * R.
c
c  Discussion:
c
c    The sphere may have nonunit radius, but it must be centered at 0.
c
c    The monomial is F(X) = X(1)^E(1) * X(2)^E(2) * ... * X(N)^E(N).
c
c    This routine is useful for testing the accuracy of quadrature
c    rules on the sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
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
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, integer E(N), the exponents of X, Y and Z in 
c    the monomial.  Each exponent must be nonnegative.
c
c    Output, double precision INTEGRAL, the integral.
c
      implicit none

      integer n

      integer e(n)
      integer e_sum
      integer i
      integer i4vec_sum
      double precision integral
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_gamma

      do i = 1, n
        if ( e(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPHERE_MONOMIAL_INT_ND - Fatal error!'
          write ( *, '(a)' ) '  All exponents must be nonnegative.'
          stop
        end if
      end do

      e_sum = i4vec_sum ( n, e )

      if ( e_sum .eq. 0 ) then
        integral = 2.0D+00 * sqrt ( pi**n ) 
     &    / r8_gamma ( 0.5D+00 * dble ( n ) ) * r**2
        return
      end if

      do i = 1, n
        if ( mod ( e(i), 2 ) .eq. 1 ) then
          integral = 0.0D+00
          return
        end if
      end do

      integral = 2.0D+00

      do i = 1, n
        integral = integral * r8_gamma ( 0.5D+00 * dble ( e(i) + 1 ) )
      end do

      integral = integral / r8_gamma ( 0.5D+00 * dble ( e_sum + n ) )
     &  * r**( e_sum + 2 )

      return
      end
      subroutine sphere_shell_03_nd ( func, n, center, r1, r2, result )

c*********************************************************************72
c
cc SPHERE_SHELL_03_ND approximates an integral inside a spherical shell in ND.
c
c  Integration region:
c
c    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
c
c  Discussion:
c
c    An 2*N point 3-rd degree formula is used, Stroud number SN-Shell:3-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F at the N-vector X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision CENTER(N), the center of the spheres.
c
c    Input, double precision R1, R2, the inner and outer radiuses that
c    define the spherical shell.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision center(n)
      double precision func
      external func
      integer i
      double precision quad
      double precision r
      double precision r1
      double precision r2
      double precision result
      double precision rho
      double precision sphere_shell_volume_nd
      double precision volume
      double precision w
      double precision x(n)

      if ( r1 .eq. r2 ) then
        result = 0.0D+00
        return
      end if

      rho = r1 / r2

      r = dble ( n ) * ( 1.0D+00 - rho**(n+2) ) 
     &  / ( dble ( n + 2 ) * ( 1.0D+00 - rho**n ) )
      r = sqrt ( r )
      w = 1.0D+00 / dble ( 2 * n )

      do i = 1, n
        x(i) = center(i)
      end do

      quad = 0.0D+00
      do i = 1, n
        x(i) = center(i) + r * r2
        quad = quad + w * func ( n, x )
        x(i) = center(i) - r * r2
        quad = quad + w * func ( n, x )
        x(i) = center(i)
      end do

      volume = sphere_shell_volume_nd ( n, r1, r2 )
      result = quad * volume

      return
      end
      function sphere_shell_volume_nd ( n, r1, r2 )

c*********************************************************************72
c
cc SPHERE_SHELL_VOLUME_ND computes the volume of a spherical shell in ND.
c
c  Integration region:
c
c    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the space.
c
c    Input, double precision R1, R2, the radiuses of the inner and 
c    outer spheres.
c
c    Output, double precision SPHERE_SHELL_VOLUME_ND, the volume of the
c    spherical shell.
c
      implicit none

      double precision ball_volume_nd
      integer n
      double precision r1
      double precision r2
      double precision sphere_shell_volume_nd

      sphere_shell_volume_nd = ball_volume_nd ( n, r2 ) 
     &  - ball_volume_nd ( n, r1 )

      return
      end
      subroutine sphere_unit_03_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_03_ND approximates integral on surface of the unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1.
c
c  Discussion:
c
c    A 2*N point 3rd degree formula is used, Stroud number UN:3-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision sphere_unit_area_nd
      double precision volume
      double precision w
      double precision x(n)

      do i = 1, n
        x(i) = 0.0D+00
      end do

      w = 1.0D+00 / dble ( 2 * n )

      quad = 0.0D+00
      do i = 1, n
        x(i) = 1.0D+00
        quad = quad + w * func ( n, x )
        x(i) = -1.0D+00
        quad = quad + w * func ( n, x )
        x(i) = 0.0D+00
      end do

      volume = sphere_unit_area_nd ( n )
      result = quad * volume

      return
      end
      subroutine sphere_unit_04_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_04_ND approximates integral on surface of the unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1.
c
c  Discussion:
c
c    A 2*N*N point 5th degree formula is used, Stroud number UN:5-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      integer j
      double precision quad
      double precision result
      double precision s
      double precision sphere_unit_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision x(n)

      do i = 1, n
        x(i) = 0.0D+00
      end do

      w1 = dble ( 4 - n ) / dble ( 2 * n * ( n + 2 ) )

      quad = 0.0D+00

      do i = 1, n
        x(i) = 1.0D+00
        quad = quad + w1 * func ( n, x )
        x(i) = -1.0D+00
        quad = quad + w1 * func ( n, x )
        x(i) = 0.0D+00
      end do

      s = 1.0D+00 / sqrt ( 2.0D+00 )
      w2 = 1.0D+00 / dble ( n * ( n + 2 ) )

      do i = 1, n

        x(i) = s

        do j = i + 1, n
          x(j) = s
          quad = quad + w2 * func ( n, x )
          x(j) = -s
          quad = quad + w2 * func ( n, x )
          x(j) = 0.0D+00
        end do

        x(i) = -s

        do j = i + 1, n
          x(j) = s
          quad = quad + w2 * func ( n, x )
          x(j) = -s
          quad = quad + w2 * func ( n, x )
          x(j) = 0.0D+00
        end do

        x(i) = 0.0D+00

      end do

      volume = sphere_unit_area_nd ( n )
      result = quad * volume

      return
      end
      subroutine sphere_unit_05_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_05_ND approximates integral on surface of the unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1.
c
c  Discussion:
c
c    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      integer iadd
      integer ihi
      integer ix(n)
      logical more
      integer ncard
      double precision quad
      double precision result
      double precision sphere_unit_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 1.0D+00
      x2 = 1.0D+00 / sqrt ( dble ( n ) )

      w1 = 1.0D+00 / dble ( n * ( n + 2 ) )
      w2 = dble ( n ) / dble ( ( n + 2 ) * 2**n )

      do i = 1, n
        x(i) = 0.0D+00
      end do

      quad = 0.0D+00

      do i = 1, n
        x(i) = x1
        quad = quad + w1 * func ( n, x )
        x(i) = -x1
        quad = quad + w1 * func ( n, x )
        x(i) = 0.0D+00
      end do

      more = .false.
      ihi = 2**n

      do i = 1, n
        x(i) = -x2
      end do

      do i = 1, ihi

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = -x(iadd)
        end if

        quad = quad + w2 * func ( n, x )

      end do

      volume = sphere_unit_area_nd ( n )
      result = quad * volume

      return
      end
      subroutine sphere_unit_07_3d ( func, result )

c*********************************************************************72
c
cc SPHERE_UNIT_07_3D approximates integral on surface of the unit sphere in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z = 1.
c
c  Discussion:
c
c    A 32 point 7-th degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order1
      parameter ( order1 = 2 )
      integer order2
      parameter ( order2 = 4 )
      integer order3 
      parameter ( order3 = 4 )

      double precision angle
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision result
      double precision sphere_unit_area_3d
      double precision volume
      double precision weight1(order1)
      double precision weight2(order2)
      double precision weight3(order3)
      double precision x
      double precision xtab1(order1)
      double precision xtab2(order2)
      double precision xtab3(order3)
      double precision y
      double precision z
c
c  Set XTAB1 and WATE1.
c
      xtab1(1) = -1.0D+00
      xtab1(2) =  1.0D+00
      weight1(1) = 1.0D+00
      weight1(2) = 1.0D+00
c
c  Set XTAB2 and WATE2.
c
      do j = 1, order2
        angle = pi * dble ( 2 * j - 1 ) / dble ( 2 * order2 )
        xtab2(j) = cos ( angle )
      end do

      weight2(1:order2) = 1.0D+00 / dble ( 4 * order2 )
c
c  Set XTAB3 and WATE3.
c
      call legendre_set ( order3, xtab3, weight3 )

      quad = 0.0D+00
      do i = 1, order1
        do j = 1, order2
          do k = 1, order3

            x = xtab1(i) * sqrt ( 1.0D+00 - xtab2(j) * xtab2(j) ) 
     &                   * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
            y = xtab1(i) * xtab2(j) 
     &        * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
            z = xtab1(i) * xtab3(k)

            quad = quad + weight1(i) * weight2(j) * weight3(k) 
     &        * func ( x, y, z )

          end do
        end do
      end do

      volume = sphere_unit_area_3d ( )
      result = quad * volume

      return
      end
      subroutine sphere_unit_07_1_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_07_1_ND approximates integral on surface of unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1.
c
c  Discussion:
c
c    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer i
      integer iadd
      integer ix(n)
      integer j
      integer jhi
      logical more
      integer ncard
      double precision quad
      double precision result
      double precision sphere_unit_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision w3
      double precision x(n)
      double precision x1
      double precision x2
      double precision x3

      w1 = dble ( 8 - n ) / dble ( n * ( n + 2 ) * ( n + 4 ) )
      w2 = dble ( n**3 ) 
     &  / dble ( 2**n * n * ( n + 2 ) * ( n + 4 ) )
      w3 = 4.0D+00 / dble ( n * ( n + 2 ) * ( n + 4 ) )

      x1 = 1.0D+00
      x2 = 1.0D+00 / sqrt ( dble ( n ) )
      x3 = 1.0D+00 / sqrt ( 2.0D+00 )

      do i = 1, n
        x(i) = 0.0D+00
      end do

      quad = 0.0D+00
c
c  First term.
c
      do i = 1, n
        x(i) = x1
        quad = quad + w1 * func ( n, x )
        x(i) = -x1
        quad = quad + w1 * func ( n, x )
        x(i) = 0.0D+00
      end do
c
c  Second term.
c
      do i = 1, n
        x(i) = -x2
      end do

      more = .false.
      jhi = 2**n

      do j = 1, jhi

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = -x(iadd)
        end if

        quad = quad + w2 * func ( n, x )

      end do
c
c  Third term.
c
      do i = 1, n
        x(i) = 0.0D+00
      end do

      do i = 1, n - 1
        do j = i + 1, n
          x(i) = x3
          x(j) = x3
          quad = quad + w3 * func ( n, x )
          x(i) = -x3
          x(j) = x3
          quad = quad + w3 * func ( n, x )
          x(i) = x3
          x(j) = -x3
          quad = quad + w3 * func ( n, x )
          x(i) = -x3
          x(j) = -x3
          quad = quad + w3 * func ( n, x )
          x(i) = 0.0D+00
          x(j) = 0.0D+00
        end do
      end do

      volume = sphere_unit_area_nd ( n )
      result = quad * volume

      return
      end
      subroutine sphere_unit_07_2_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_07_2_ND approximates integral on surface of unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1.
c
c  Discussion:
c
c    A 2^N * ( N + 1 ) point 7th degree formula is used, Stroud number UN:7-2.
c
c    Some of the weights in this quadrature formula are negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X), at the N dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision func
      external func
      integer iadd
      integer i
      integer ix(n)
      integer j
      integer jhi
      logical more
      integer ncard
      double precision quad
      double precision result
      double precision sphere_unit_area_nd
      double precision volume
      double precision w1
      double precision w2
      double precision x(n)
      double precision x1
      double precision x2
      double precision x3

      do i = 1, n
        x(i) = 0.0D+00
      end do

      w1 = - dble ( n * n ) / dble ( 2**(n+3) * ( n + 2 ) )
      w2 = dble ( ( n + 4 ) * ( n + 4 ) ) 
     &  / dble ( 2**(n+3) * n * ( n + 2 ) )
      x1 = 1.0D+00 / sqrt ( dble ( n ) )
      x2 = sqrt ( 5.0D+00 / dble ( n + 4 ) )
      x3 = 1.0D+00 / sqrt ( dble ( n + 4 ) )

      quad = 0.0D+00

      do i = 1, n
        x(i) = - x1
      end do

      more = .false.
      jhi = 2**n

      do j = 1, jhi

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = - x(iadd)
        end if

        quad = quad + w1 * func ( n, x )

      end do

      do i = 1, n

        do j = 1, n
          x(j) = - x3
        end do

        x(i) = - x2
        more = .false.

        do j = 1, jhi

          call subset_gray_next ( n, ix, more, ncard, iadd )

          if ( iadd .ne. 0 ) then
            x(iadd) = - x(iadd)
          end if

          quad = quad + w2 * func ( n, x )

        end do

      end do

      volume = sphere_unit_area_nd ( n )
      result = quad * volume

      return
      end
      subroutine sphere_unit_11_3d ( func, result )

c*********************************************************************72
c
cc SPHERE_UNIT_11_3D approximates integral on surface of unit sphere in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z = 1.
c
c  Discussion:
c
c    A 50 point 11-th degree formula is used, Stroud number U3:11-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    AD McLaren,
c    Mathematics of Computation,
c    Volume 17, pages 361-383, 1963.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision func
      external func
      integer i
      integer j
      integer k
      integer l
      double precision quad
      double precision result
      double precision sphere_unit_area_3d
      double precision volume
      double precision w1
      double precision w2
      double precision w3
      double precision w4
      double precision x
      double precision y
      double precision z

      quad = 0.0D+00

      w1 = 9216.0D+00 / 725760.0D+00
      x = 1.0D+00
      y = 0.0D+00
      z = 0.0D+00
      do i = 1, 2
        x = -x
        do j = 1, 3
          call r8_swap3 ( x, y, z )
          quad = quad + w1 * func ( x, y, z )
        end do
      end do

      w2 = 16384.0D+00 / 725760.0D+00
      x = sqrt ( 0.5D+00 )
      y = sqrt ( 0.5D+00 )
      z = 0.0D+00
      do i = 1, 2
        x = -x
        do j = 1, 2
          y = -y
          do k = 1, 3
            call r8_swap3 ( x, y, z )
            quad = quad + w2 * func ( x, y, z )
          end do
        end do
      end do

      w3 = 15309.0D+00 / 725760.0D+00
      x = sqrt ( 1.0D+00 / 3.0D+00 )
      y = sqrt ( 1.0D+00 / 3.0D+00 )
      z = sqrt ( 1.0D+00 / 3.0D+00 )
      do i = 1, 2
        x = -x
        do j = 1, 2
          y = -y
          do k = 1, 2
            z = -z
            quad = quad + w3 * func ( x, y, z )
          end do
        end do
      end do

      w4 = 14641.0D+00 / 725760.0D+00
      x = sqrt ( 1.0D+00 / 11.0D+00 )
      y = sqrt ( 1.0D+00 / 11.0D+00 )
      z = 3.0D+00 * sqrt ( 1.0D+00 / 11.0D+00 )
      do i = 1, 2
        x = -x
        do j = 1, 2
          y = -y
          do k = 1, 2
            z = -z
            do l = 1, 3
              call r8_swap3 ( x, y, z )
              quad = quad + w4 * func ( x, y, z )
            end do
          end do
        end do
      end do

      volume = sphere_unit_area_3d ( )
      result = quad * volume

      return
      end
      subroutine sphere_unit_11_nd ( func, n, result )

c*********************************************************************72
c
cc SPHERE_UNIT_11_ND approximates integral on surface of unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) = 1
c
c  Discussion:
c
c    An 2^N * ( N^2 + N + 1 ) point formula of degree 5 is used.
c
c    (For N = 3, the number of points is actually only 56, and
c     for N = 4, the number of points is actually only 240.)
c
c    One element of COEF31 was changed from
c      0.0236339091329 to
c      0.0236639091329
c    by Stroud, when going from his paper to his later textbook.
c    This correction was pointed out by David Wright, 16 February 2010.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    A Fifth Degree Integration Formula for the N-Simplex,
c    SIAM Journal on Numerical Analysis,
c    Volume 6, Number 1, March 1969.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X) at the N-dimensional point
c    X, of the form
c      function func ( n, x )
c      integer n
c      double precision func
c      double precision x(n)
c
c    Input, integer N, the dimension of the space.  For this 
c    routine, it must be the case that 3 <= N <= 16.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer n

      double precision area
      double precision coef1(3:16)
      double precision coef21(3:16)
      double precision coef22(3:16)
      double precision coef31(3:16)
      double precision coef32(3:16)
      double precision func
      external func
      integer i
      integer iadd
      integer ix(n)
      integer j
      integer k
      logical more
      integer ncard
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision s1
      double precision s2
      double precision sphere_unit_area_nd
      double precision u1
      double precision u2
      double precision v1
      double precision v2
      double precision x(n)

      save coef1
      save coef21
      save coef22
      save coef31
      save coef32

      data coef1 /
     &  0.128571428571D+00, 
     &  0.0518518518518D+00, 
     &  0.0211979378646D+00, 
     &  0.281250000000D+00, 
     &  1.11934731935D+00, 
     &  2.82751322751D+00, 
     &  5.68266145619D+00, 
     &  9.93785824515D+00, 
     &  15.8196616478D+00, 
     &  23.5285714285D+00, 
     &  33.2409299392D+00, 
     &  45.1113811729D+00, 
     &  59.2754264177D+00, 
     &  75.8518518518D+00 /
      data coef21 /
     &  0.163795782462D+00, 
     &  0.0967270533860D+00, 
     &  0.0638253880175D+00, 
     &  0.0452340041459D+00, 
     &  0.0336329118818D+00, 
     &  0.0261275095270D+00, 
     &  0.0208331595340D+00, 
     &  0.0169937111647D+00, 
     &  0.0141147212492D+00, 
     &  0.0118949128383D+00, 
     &  0.0101424250926D+00, 
     &  0.00873046796644D+00, 
     &  0.00757257014768D+00, 
     &  0.00660813369775D+00 /
      data coef22 /
     &  0.126680408014D+00, 
     &  0.0514210947621D+00, 
     &  0.0213579471658D+00, 
     & -0.108726067638D+00, 
     & -0.371589499738D+00, 
     & -0.786048144448D+00, 
     & -1.36034060198D+00, 
     & -2.09547695631D+00, 
     & -2.98784764467D+00, 
     & -4.03107480702D+00, 
     & -5.21726499521D+00, 
     & -6.53783099707D+00, 
     & -7.98401677102D+00, 
     & -9.54722261180D+00 /
      data coef31 /
     &  0.0D+00, 
     &  0.0592592592592D+00, 
     &  0.0236639091329D+00, 
     &  0.0525940190875D+00, 
     &  0.0925052768546D+00, 
     &  0.141316953438D+00, 
     &  0.196818580052D+00, 
     &  0.257027634179D+00, 
     &  0.320299222258D+00, 
     &  0.385326226441D+00, 
     &  0.451098131789D+00, 
     &  0.516849445559D+00, 
     &  0.582010515746D+00, 
     &  0.646165210110D+00 /
      data coef32 /
     &  0.0D+00, 
     &  0.0D+00, 
     &  0.0316246294890D+00, 
     &  0.0207194729760D+00, 
     &  0.0144303800811D+00, 
     &  0.0105348984135D+00, 
     &  0.00798435122193D+00, 
     &  0.00623845929545D+00, 
     &  0.00499896882962D+00, 
     &  0.00409176297655D+00, 
     &  0.00341037426698D+00, 
     &  0.00288710646943D+00, 
     &  0.00247745182907D+00, 
     &  0.00215128820597D+00 /

      if ( n .lt. 3 .or. 16 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPHERE_UNIT_11_ND - Fatal error!'
        write ( *, '(a)' ) '  Input spatial dimension N out of range.'
        write ( *, '(a,i8)' ) '  N = ', n
        result = 0.0D+00
        return
      end if

      quad = 0.0D+00
c
c  S1
c
      do i = 1, n
        x(i) = 1.0D+00 / sqrt ( dble ( n ) )
      end do

      more = .false.

10    continue

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd .ne. 0 ) then
          x(iadd) = -x(iadd)
        end if

        quad = quad + coef1(n) * func ( n, x )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue
c
c  S21
c
      r1 = ( dble ( n + 6 ) - 4.0D+00 * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 12 * n - 12 )
      r1 = sqrt ( r1 )

      s1 = ( dble ( 7 * n - 6 ) 
     &  + dble ( 4 * ( n - 1 ) ) * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 12 * n - 12 )
      s1 = sqrt ( s1 )

      do i = 1, n

        do j = 1, n
          x(j) = r1
        end do
        x(i) = s1

        more = .false.

30      continue

          call subset_gray_next ( n, ix, more, ncard, iadd )

          if ( iadd .ne. 0 ) then
            x(iadd) = -x(iadd)
          end if

          quad = quad + coef21(n) * func ( n, x )

          if ( .not. more ) then
            go to 40
          end if

        go to 30

40      continue

      end do
c
c  S22
c
      r2 = ( dble ( n + 6 ) + 4.0D+00 * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 12 * n - 12 )
      r2 = sqrt ( r2 )

      s2 = ( dble ( 7 * n - 6 ) 
     &  - dble ( 4 * ( n - 1 ) ) * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 12 * n - 12 )
      s2 = sqrt ( s2 )

      do i = 1, n

        do j = 1, n
          x(j) = r2
        end do
        x(i) = s2

        more = .false.

50      continue

          call subset_gray_next ( n, ix, more, ncard, iadd )

          if ( iadd .ne. 0 ) then
            x(iadd) = -x(iadd)
          end if

          quad = quad + coef22(n) * func ( n, x )

          if ( .not. more ) then
            go to 60
          end if

        go to 50

60      continue

      end do
c
c  S31
c
      u1 = ( dble ( n + 12 ) + 8.0D+00 * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 24 * n - 48 )
      u1 = sqrt ( u1 )

      v1 = ( dble ( 7 * n - 12 ) 
     &  - dble ( 4 * n - 8 ) * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 24 * n - 48 )
      v1 = sqrt ( v1 )

      do i = 1, n

        do j = i + 1, n

          do k = 1, n
            x(k) = u1
          end do
          x(i) = v1
          x(j) = v1

          more = .false.

70        continue

            call subset_gray_next ( n, ix, more, ncard, iadd )

            if ( iadd .ne. 0 ) then
              x(iadd) = -x(iadd)
            end if

            quad = quad + coef31(n) * func ( n, x )

            if ( .not. more ) then
              go to 80
            end if

          go to 70

80        continue

        end do

      end do
c
c  S32
c
      u2 = ( dble ( n + 12 ) - 8.0D+00 * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 24 * n - 48 )
      u2 = sqrt ( u2 )

      v2 = ( dble ( 7 * n - 12 ) 
     &  + dble ( 4 * n - 8 ) * sqrt ( 3.0D+00 ) ) 
     &  / dble ( n * n + 24 * n - 48 )
      v2 = sqrt ( v2 )

      do i = 1, n

        do j = i + 1, n

          do k = 1, n
            x(k) = u2
          end do
          x(i) = v2
          x(j) = v2

          more = .false.

 90       continue

            call subset_gray_next ( n, ix, more, ncard, iadd )

            if ( iadd .ne. 0 ) then
              x(iadd) = -x(iadd)
            end if

            quad = quad + coef32(n) * func ( n, x )

            if ( .not. more ) then
              go to 100
            end if

          go to 90

100       continue

        end do

      end do

      area = sphere_unit_area_nd ( n )

      result = quad * area / 2.0D+00**n

      return
      end
      subroutine sphere_unit_14_3d ( func, result )

c*********************************************************************72
c
cc SPHERE_UNIT_14_3D approximates integral on surface of unit sphere in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z = 1.
c
c  Discussion:
c
c    A 72 point 14-th degree formula is used, Stroud number U3:14-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    AD McLaren,
c    Mathematics of Computation,
c    Volume 17, pages 361-383, 1963.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision func
      external func
      integer i
      integer j
      integer k
      double precision quad
      double precision result
      double precision sphere_unit_area_3d
      double precision temp
      double precision volume
      double precision w1
      double precision w2
      double precision x
      double precision xtab(5)
      double precision y
      double precision ytab(5)
      double precision z
      double precision ztab(5)

      save xtab
      save ytab
      save ztab

      data xtab /
     &  -0.151108275D+00, 0.315838353D+00, 0.346307112D+00, 
     &  -0.101808787D+00, -0.409228403D+00 /
      data ytab /
     &  0.155240600D+00, 0.257049387D+00, 0.666277790D+00,  
     &  0.817386065D+00, 0.501547712D+00 /
      data ztab /
     &  0.976251323D+00, 0.913330032D+00, 0.660412970D+00,  
     &  0.567022920D+00, 0.762221757D+00 /

      quad = 0.0D+00

      w1 = 125.0D+00 / 10080.0D+00
      x = 0.525731112D+00
      y = 0.850650808D+00
      z = 0.0D+00

      do i = 1, 2
        x = -x
        do j = 1, 2
          y = -y
          do k = 1, 3
            call r8_swap3 ( x, y, z )
            quad = quad + w1 * func ( x, y, z )
          end do
        end do
      end do

      w2 = 143.0D+00 / 10080.0D+00

      do i = 1, 5

        x = xtab(i)
        y = ytab(i)
        z = ztab(i)

        do j = 1, 3

          temp = x
          x = z
          z = -y
          y = -temp

          do k = 1, 3
            call r8_swap3 ( x, y, z )
            quad = quad + w2 * func ( x, y, z )
          end do

          y = -y
          z = -z
          quad = quad + w2 * func ( x, y, z )

        end do

      end do

      volume = sphere_unit_area_3d ( )
      result = quad * volume

      return
      end
      subroutine sphere_unit_15_3d ( func, result )

c*********************************************************************72
c
cc SPHERE_UNIT_15_3D approximates integral on surface of unit sphere in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z = 1.
c
c  Discussion:
c
c    A 128 point 15-th degree spherical product Gauss formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function which evaluates F(X,Y,Z), of the form
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 8 )

      double precision angle
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision result
      double precision sphere_unit_area_3d
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y
      double precision z

      call legendre_set ( order, xtab, weight )

      do i = 1, order
        weight(i) = weight(i) / 32.0D+00
      end do

      quad = 0.0D+00

      do j = 1, order

        do k = 1, 16

          angle = dble ( k ) * pi / 8.0D+00
          x = sqrt ( 1.0D+00 - xtab(j)**2 ) * cos ( angle )
          y = sqrt ( 1.0D+00 - xtab(j)**2 ) * sin ( angle )
          z = xtab(j)

          quad = quad + weight(j) * func ( x, y, z )

        end do
      end do

      volume = sphere_unit_area_3d ( )
      result = quad * volume

      return
      end
      function sphere_unit_area_3d ( )

c*********************************************************************72
c
cc SPHERE_UNIT_AREA_3D computes the surface area of the unit sphere in 3D.
c
c  Integration region:
c
c    X*X + Y*Y + Z*Z = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision SPHERE_UNIT_AREA_3D, the area of the sphere.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sphere_unit_area_3d

      sphere_unit_area_3d = 4.0D+00 * pi

      return
      end
      function sphere_unit_area_nd ( dim_num )

c*********************************************************************72
c
cc SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
c
c  Discussion:
c
c    The unit sphere in ND satisfies:
c
c      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
c
c    Results for the first few values of N are:
c
c    DIM_NUM   Area
c
c     2    2        * PI
c     3    4        * PI
c     4  ( 2 /   1) * PI^2
c     5  ( 8 /   3) * PI^2
c     6  ( 1 /   1) * PI^3
c     7  (16 /  15) * PI^3
c     8  ( 1 /   3) * PI^4
c     9  (32 / 105) * PI^4
c    10  ( 1 /  12) * PI^5
c
c    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
c
c    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI^(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Output, double precision SPHERE_UNIT_AREA_ND, the area of the sphere.
c
      implicit none

      double precision area
      integer dim_num
      integer i
      integer m
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sphere_unit_area_nd

      if ( mod ( dim_num, 2 ) .eq. 0 ) then
        m = dim_num / 2
        area = 2.0D+00 * ( pi )**m
        do i = 1, m-1
          area = area / dble ( i )
        end do
      else
        m = ( dim_num - 1 ) / 2
        area = ( pi )**m * 2.0D+00**dim_num
        do i = m+1, 2*m
          area = area / dble ( i )
        end do
      end if

      sphere_unit_area_nd = area

      return
      end
      subroutine sphere_unit_area_values ( n_data, n, area )

c*********************************************************************72
c
cc SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
c
c  Discussion:
c
c    The formula for the surface area of the unit sphere in N dimensions is:
c
c      Sphere_Unit_Area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
c
c    Some values of the function include:
c
c       N   Area
c
c       2    2        * PI
c       3  ( 4 /    ) * PI
c       4  ( 2 /   1) * PI^2
c       5  ( 8 /   3) * PI^2
c       6  ( 1 /   1) * PI^3
c       7  (16 /  15) * PI^3
c       8  ( 1 /   3) * PI^4
c       9  (32 / 105) * PI^4
c      10  ( 1 /  12) * PI^5
c
c    For the unit sphere, Area(N) = N * Volume(N)
c
c    In Mathematica, the function can be evaluated by:
c
c      2 * Pi^(n/2) / Gamma[n/2]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the spatial dimension.
c
c    Output, double precision AREA, the area of the unit sphere
c    in that dimension.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision area
      double precision area_vec(n_max)
      integer n_data
      integer n
      integer n_vec(n_max)

      save area_vec
      save n_vec

      data area_vec /
     &  0.2000000000000000D+01,
     &  0.6283185307179586D+01,
     &  0.1256637061435917D+02,
     &  0.1973920880217872D+02,
     &  0.2631894506957162D+02,
     &  0.3100627668029982D+02,
     &  0.3307336179231981D+02,
     &  0.3246969701133415D+02,
     &  0.2968658012464836D+02,
     &  0.2550164039877345D+02,
     &  0.2072514267328890D+02,
     &  0.1602315322625507D+02,
     &  0.1183817381218268D+02,
     &  0.8389703410491089D+01,
     &  0.5721649212349567D+01,
     &  0.3765290085742291D+01,
     &  0.2396678817591364D+01,
     &  0.1478625959000308D+01,
     &  0.8858104195716824D+00,
     &  0.5161378278002812D+00 /
      data n_vec /
     &   1,
     &   2,
     &   3,
     &   4,
     &   5,
     &   6,
     &   7,
     &   8,
     &   9,
     &  10,
     &  11,
     &  12,
     &  13,
     &  14,
     &  15,
     &  16,
     &  17,
     &  18,
     &  19,
     &  20 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        area = 0.0D+00
      else
        n = n_vec(n_data)
        area = area_vec(n_data)
      end if

      return
      end
      function sphere_unit_monomial_nd ( n, p )

c*********************************************************************72
c
cc SPHERE_UNIT_MONOMIAL_ND integrate monomial on surface of unit sphere in ND.
c
c  Integration region:
c
c    sum ( X(1:N)^2 ) == 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gerald Folland,
c    How to Integrate a Polynomial Over a Sphere,
c    American Mathematical Monthly,
c    Volume 108, May 2001, pages 446-448.
c
c  Parameters:
c
c    Input, integer  N, the dimension of the space.
c
c    Input, integer P(N), the exponents of X(1) through X(N) 
c    in the monomial.  The exponents P(N) must be nonnegative.
c
c    Output, double precision SPHERE_UNIT_MONOMIAL_ND, the integral of
c    X1^P(1) * X2^P(2) * ... * XN^P(N) over the unit sphere.
c
      implicit none

      integer n

      double precision arg1
      double precision arg2
      integer i
      integer p(n)
      double precision r8_gamma_log
      double precision sphere_unit_monomial_nd
      double precision temp

      sphere_unit_monomial_nd = 0.0D+00

      do i = 1, n
        if ( mod ( p(i), 2 ) .eq. 1 ) then
          return
        end if
      end do

      temp = 0.0D+00
      arg2 = 0.0D+00

      do i = 1, n
        arg1 = dble ( p(i) + 1 ) / 2.0D+00
        temp = temp + r8_gamma_log ( arg1 )
        arg2 = arg2 + arg1
      end do
      temp = temp - r8_gamma_log ( arg2 )
  
      sphere_unit_monomial_nd = 2.0D+00 * exp ( temp )

      return
      end
      function sphere_unit_volume_nd ( dim_num )

c*********************************************************************72
c
cc SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
c
c  Discussion:
c
c    The unit sphere in ND satisfies:
c
c      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
c
c    Results for the first few values of DIM_NUM are:
c
c     DIM_NUM  Volume
c
c     1    2
c     2    1        * PI
c     3  ( 4 /   3) * PI
c     4  ( 1 /   2) * PI^2
c     5  ( 8 /  15) * PI^2
c     6  ( 1 /   6) * PI^3
c     7  (16 / 105) * PI^3
c     8  ( 1 /  24) * PI^4
c     9  (32 / 945) * PI^4
c    10  ( 1 / 120) * PI^5
c
c    For the unit sphere, Volume(DIM_NUM) = 2 * PI * Volume(DIM_NUM-2)/ DIM_NUM
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Output, double precision SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
c
      implicit none

      integer dim_num
      integer i
      integer m
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sphere_unit_volume_nd
      double precision volume

      if ( mod ( dim_num, 2 ) .eq. 0 ) then
        m = dim_num / 2
        volume = ( pi )**m
        do i = 1, m
          volume = volume / dble ( i )
        end do
      else
        m = ( dim_num - 1 ) / 2
        volume = ( pi )**m * 2.0D+00**dim_num
        do i = m+1, 2*m+1
          volume = volume / dble ( i )
        end do
      end if

      sphere_unit_volume_nd = volume

      return
      end
      subroutine sphere_unit_volume_values ( n_data, n, volume )

c*********************************************************************72
c
cc SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
c
c  Discussion:
c
c    The formula for the volume of the unit sphere in N dimensions is
c
c      Volume(N) = 2 * PI^(N/2) / ( N * Gamma ( N / 2 ) )
c
c    This function satisfies the relationships:
c
c      Volume(N) = 2 * PI * Volume(N-2) / N
c      Volume(N) = Area(N) / N
c
c    Some values of the function include:
c
c       N  Volume
c
c       1    1
c       2    1        * PI
c       3  ( 4 /   3) * PI
c       4  ( 1 /   2) * PI^2
c       5  ( 8 /  15) * PI^2
c       6  ( 1 /   6) * PI^3
c       7  (16 / 105) * PI^3
c       8  ( 1 /  24) * PI^4
c       9  (32 / 945) * PI^4
c      10  ( 1 / 120) * PI^5
c
c    In Mathematica, the function can be evaluated by:
c
c      2 * Pi^(n/2) / ( n * Gamma[n/2] )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the spatial dimension.
c
c    Output, double precision VOLUME, the volume of the unit
c    sphere in that dimension.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer n_data
      integer n
      integer n_vec(n_max)
      double precision volume
      double precision volume_vec(n_max)

      save n_vec
      save volume_vec

      data n_vec /
     &   1,  2,
     &   3,  4,
     &   5,  6,
     &   7,  8,
     &   9, 10,
     &  11, 12,
     &  13, 14,
     &  15, 16,
     &  17, 18,
     &  19, 20 /
      data volume_vec /
     &  0.2000000000000000D+01,
     &  0.3141592653589793D+01,
     &  0.4188790204786391D+01,
     &  0.4934802200544679D+01,
     &  0.5263789013914325D+01,
     &  0.5167712780049970D+01,
     &  0.4724765970331401D+01,
     &  0.4058712126416768D+01,
     &  0.3298508902738707D+01,
     &  0.2550164039877345D+01,
     &  0.1884103879389900D+01,
     &  0.1335262768854589D+01,
     &  0.9106287547832831D+00,
     &  0.5992645293207921D+00,
     &  0.3814432808233045D+00,
     &  0.2353306303588932D+00,
     &  0.1409811069171390D+00,
     &  0.8214588661112823D-01,
     &  0.4662160103008855D-01,
     &  0.2580689139001406D-01 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        volume = 0.0D+00
      else
        n = n_vec(n_data)
        volume = volume_vec(n_data)
      end if

      return
      end
      subroutine sphere_volume_2d ( r, volume )

c*********************************************************************72
c
cc SPHERE_VOLUME_2D computes the volume of an implicit sphere in 2D.
c
c  Discussion:
c
c    An implicit sphere in 2D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision VOLUME, the volume of the sphere.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume

      volume = pi * r * r

      return
      end
      subroutine sphere_volume_3d ( r, volume )

c*********************************************************************72
c
cc SPHERE_VOLUME_3D computes the volume of an implicit sphere in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision VOLUME, the volume of the sphere.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume

      volume = ( 4.0D+00 / 3.0D+00 ) * pi * r * r * r

      return
      end
      subroutine sphere_volume_nd ( dim_num, r, volume )
c
c*********************************************************************72
c
cc SPHERE_VOLUME_ND computes the volume of an implicit sphere in ND.
c
c  Discussion:
c
c    An implicit sphere in ND satisfies the equation:
c
c      sum ( ( X(1:N) - PC(1:N) )^2 ) = R^2
c
c    where R is the radius and PC is the center.
c
c    Results for the first few values of N are:
c
c    DIM_NUM  Volume
c    -     -----------------------
c    2                PI   * R^2
c    3     (4/3)    * PI   * R^3
c    4     (1/2)    * PI^2 * R^4
c    5     (8/15)   * PI^2 * R^5
c    6     (1/6)    * PI^3 * R^6
c    7     (16/105) * PI^3 * R^7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision VOLUME, the volume of the sphere.
c
      implicit none

      integer dim_num
      double precision r
      double precision sphere_unit_volume_nd
      double precision volume

      volume = r**dim_num * sphere_unit_volume_nd ( dim_num )

      return
      end
      subroutine square_sum ( func, center, r, order, xtab, ytab, 
     &  weight, result )

c*********************************************************************72
c
cc SQUARE_SUM carries out a quadrature rule over a square.
c
c  Integration region:
c
c      abs ( X - CENTER(1) ) <= R 
c    and
c      abs ( Y - CENTER(2) ) <= R
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the function to be
c    integrated.  The user must declare the name an EXTERNAL
c    parameter in the calling program, pass the name of the
c    function in FUNC, and write a function of the form
c      function func(x,y)
c    which evaluates the function at the point (X,Y).
c
c    Input, double precision CENTER(2), the center of the square.
c
c    Input, double precision R, the radius of the square.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas of 
c    the rule.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer order

      double precision center(dim_num)
      double precision func
      external func
      integer i
      double precision quad
      double precision r
      double precision result
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision y
      double precision ytab(order)

      quad = 0.0D+00
      do i = 1, order
        x = center(1) + r * xtab(i)
        y = center(2) + r * ytab(i)
        quad = quad + 0.25D+00 * weight(i) * func ( x, y )
      end do

      volume = 4.0D+00 * r * r
      result = quad * volume

      return
      end
      subroutine square_unit_set ( rule, order, xtab, ytab, wtab )

c*********************************************************************72
c
cc SQUARE_UNIT_SET sets quadrature weights and abscissas in the unit square.
c
c  Discussion;
c
c    To get the value of ORDER associated with a given rule, 
c    call SQUARE_UNIT_SIZE first.
c
c  Integration region:
c
c      -1 <= X <= 1,
c    and
c      -1 <= Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gilbert Strang, George Fix,
c    An Analysis of the Finite Element Method,
c    Cambridge, 1973,
c    ISBN: 096140888X,
c    LC: TA335.S77.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule number.
c    1, order 1, degree 1 rule.
c    2, order 4, degree 3, rule.
c    3, order 9, degree 5 rule.
c    4, order 12 degree 7 rule, Stroud number C2:7-1.
c    5, order 13 degree 7 rule, Stroud number C2:7-3.
c    6, order 64 degree 15 product rule.
c
c    Input, integer ORDER, the order of the rule.
c
c    Output, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Output, double precision WTAB(ORDER), the weights.
c
      implicit none

      integer order
      integer order2
      parameter ( order2 = 8 )

      double precision a
      double precision c
      integer i
      integer j
      integer k
      double precision r
      integer rule
      double precision s
      double precision t
      double precision w1
      double precision w2
      double precision w3
      double precision wtab(order)
      double precision weight2(order2)
      double precision xtab(order)
      double precision xtab2(order2)
      double precision ytab(order)
      double precision z

      if ( rule .eq. 1 ) then

        wtab(1) = 4.0D+00

        xtab(1) = 0.0D+00
        ytab(1) = 0.0D+00

      else if ( rule .eq. 2 ) then

        a = 1.0D+00
        s = 1.0D+00 / sqrt ( 3.0D+00 )

        xtab(1) = -s
        xtab(2) = +s
        xtab(3) = -s
        xtab(4) = +s
        ytab(1) = -s
        ytab(2) = -s
        ytab(3) = +s
        ytab(4) = +s
        wtab(1) =  a
        wtab(2) =  a
        wtab(3) =  a
        wtab(4) =  a

      else if ( rule .eq. 3 ) then

        s = sqrt ( 0.6D+00 )
        z = 0.0D+00
        w1 = 64.0D+00 / 81.0D+00
        w2 = 25.0D+00 / 81.0D+00
        w3 = 40.0D+00 / 81.0D+00

        xtab(1) =  z
        xtab(2) = -s
        xtab(3) = +s
        xtab(4) = -s
        xtab(5) = +s
        xtab(6) =  z
        xtab(7) = -s
        xtab(8) = +s
        xtab(9) =  z 
        ytab(1) =  z
        ytab(2) = -s
        ytab(3) = -s
        ytab(4) = +s
        ytab(5) = +s
        ytab(6) = -s
        ytab(7) =  z
        ytab(8) =  z
        ytab(9) = +s
        wtab(1) = w1
        wtab(2) = w2
        wtab(3) = w2
        wtab(4) = w2
        wtab(5) = w2
        wtab(6) = w3
        wtab(7) = w3
        wtab(8) = w3
        wtab(9) = w3

      else if ( rule .eq. 4 ) then

        r = sqrt ( 6.0D+00 / 7.0D+00 )
        c = 3.0D+00 * sqrt ( 583.0D+00 )
        s = sqrt ( ( 114.0D+00 - c ) / 287.0D+00 )
        t = sqrt ( ( 114.0D+00 + c ) / 287.0D+00 )
        w1 = 4.0D+00 * 49.0D+00 / 810.0D+00
        w2 = 4.0D+00 * ( 178981.0D+00 + 923.0D+00 * c ) / 1888920.0D+00
        w3 = 4.0D+00 * ( 178981.0D+00 - 923.0D+00 * c ) / 1888920.0D+00
        z = 0.0D+00

        xtab(1) =  r
        xtab(2) =  z
        xtab(3) = -r
        xtab(4) =  z
        xtab(5) =  s
        xtab(6) = -s
        xtab(7) = -s
        xtab(8) =  s
        xtab(9) =  t
        xtab(10) = -t
        xtab(11) = -t
        xtab(12) =  t
        ytab(1) =  z
        ytab(2) =  r
        ytab(3) =  z
        ytab(4) = -r
        ytab(5) =  s
        ytab(6) =  s
        ytab(7) = -s
        ytab(8) = -s
        ytab(9) =  t
        ytab(10) =  t
        ytab(11) = -t
        ytab(12) = -t
        wtab(1) = w1
        wtab(2) = w1
        wtab(3) = w1
        wtab(4) = w1
        wtab(5) = w2
        wtab(6) = w2
        wtab(7) = w2
        wtab(8) = w2
        wtab(9) = w3
        wtab(10) = w3
        wtab(11) = w3
        wtab(12) = w3

      else if ( rule .eq. 5 ) then

        r = sqrt ( 12.0D+00 / 35.0D+00 )
        c = 3.0D+00 * sqrt ( 186.0D+00 )
        s = sqrt ( ( 93.0D+00 + c ) / 155.0D+00 )
        t = sqrt ( ( 93.0D+00 - c ) / 155.0D+00 )
        w1 =  8.0D+00 / 162.0D+00
        w2 = 98.0D+00 / 162.0D+00
        w3 = 31.0D+00 / 162.0D+00
        z = 0.0D+00

        xtab(1) =  z
        xtab(2) =  r
        xtab(3) = -r
        xtab(4) =  z
        xtab(5) =  z
        xtab(6) =  s
        xtab(7) =  s
        xtab(8) = -s
        xtab(9) = -s
        xtab(10) =  t
        xtab(11) =  t 
        xtab(12) = -t 
        xtab(13) = -t
        ytab(1) =  z
        ytab(2) =  z
        ytab(3) =  z
        ytab(4) =  r
        ytab(5) = -r
        ytab(6) =  t
        ytab(7) = -t
        ytab(8) =  t
        ytab(9) = -t
        ytab(10) =  s
        ytab(11) = -s
        ytab(12) =  s
        ytab(13) = -s
        wtab(1) = w1
        wtab(2) = w2
        wtab(3) = w2
        wtab(4) = w2
        wtab(5) = w2
        wtab(6) = w3
        wtab(7) = w3
        wtab(8) = w3
        wtab(9) = w3
        wtab(10) = w3
        wtab(11) = w3
        wtab(12) = w3
        wtab(13) = w3
 
      else if ( rule .eq. 6 ) then

        call legendre_set ( order2, xtab2, weight2 )

        k = 0

        do i = 1, order2

          do j = 1, order2

            k = k + 1
            xtab(k) = xtab2(i)
            ytab(k) = xtab2(j)
            wtab(k) = weight2(i) * weight2(j)

          end do

        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SQUARE_UNIT_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
        stop

      end if

      return
      end
      subroutine square_unit_size ( rule, order )

c*********************************************************************72
c
cc SQUARE_UNIT_SIZE sizes a quadrature rule in the unit square.
c
c  Integration region:
c
c      -1 <= X <= 1,
c    and
c      -1 <= Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gilbert Strang, George Fix,
c    An Analysis of the Finite Element Method,
c    Cambridge, 1973,
c    ISBN: 096140888X,
c    LC: TA335.S77.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the rule number.
c    1, a 1 point 1st degree rule.
c    2, a 4 point 3rd degree rule.
c    3, a 9 point 5th degree rule.
c    4, a 12 point 7-th degree rule, Stroud number C2:7-1.
c    5, a 13 point 7-th degree rule, Stroud number C2:7-3.
c    6, a 64 point 15-th degree product rule.
c
c    Output, integer ORDER, the order of the rule.
c
      implicit none

      integer order
      integer rule

      if ( rule .eq. 1 ) then

        order = 1

      else if ( rule .eq. 2 ) then

        order = 4

      else if ( rule .eq. 3 ) then

        order = 9

      else if ( rule .eq. 4 ) then

        order = 12

      else if ( rule .eq. 5 ) then

        order = 13

      else if ( rule .eq. 6 ) then

        order = 64

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SQUARE_UNIT_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
        stop

      end if

      return
      end
      subroutine square_unit_sum ( func, order, xtab, ytab, weight, 
     &  result )

c*********************************************************************72
c
cc SQUARE_UNIT_SUM carries out a quadrature rule over the unit square.
c
c  Integration region:
c
c      -1 <= X <= 1, 
c    and
c      -1 <= Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the function to be
c    integrated.  The user must declare the name an EXTERNAL
c    parameter in the calling program, pass the name of the
c    function in FUNC, and write a function of the form
c      function func ( x, y )
c    which evaluates the function at the point (X,Y).
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas of 
c    the rule.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision volume
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)

      quad = 0.0D+00
      do i = 1, order
        quad = quad + weight(i) * func ( xtab(i), ytab(i) ) / 4.0D+00
      end do

      volume = 1.0D+00
      result = quad * volume

      return
      end
      subroutine subset_gray_next ( n, a, more, ncard, iadd )

c*********************************************************************72
c
cc SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
c
c  Discussion:
c
c    It generates the subsets one at a time, by adding or subtracting
c    exactly one element on each step.
c
c    This uses a Gray code ordering of the subsets.
c
c    The user should set MORE = FALSE and the value of N before
c    the first call.  On return, the user may examine A which contains
c    the definition of the new subset, and must check MORE, because
c    as soon as it is FALSE on return, all the subsets have been
c    generated and the user probably should cease calling.
c
c    The first set returned is the empty set.
c
c  Example:
c
c    N = 4
c
c    0 0 0 0
c    1 0 0 0
c    1 1 0 0
c    0 1 0 0
c    0 1 1 0
c    1 1 1 0
c    1 0 1 0
c    0 0 1 0
c    0 0 1 1
c    1 0 1 1
c    1 1 1 1
c    0 1 1 1
c    0 1 0 1
c    1 1 0 1
c    1 0 0 1
c    0 0 0 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
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
c    Input, integer N, the order of the total set from which
c    subsets will be drawn.
c
c    Input/output, integer A(N).  On each return, the Gray code for the newly
c    generated subset.  A(I) = 0 if element I is in the subset, 1 otherwise.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but once
c    all the subsets have been generated, MORE will be
c    reset FALSE on return and you should stop calling the program.
c
c    Input/output, integer NCARD, the cardinality of the set returned,
c    which may be any value between 0 (the empty set) and N (the
c    whole set).
c
c    Output, integer IADD, the element which was added or removed to the
c    previous subset to generate the current one.  Exception:
c    the empty set is returned on the first call, and IADD is set to 0.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer iadd
      logical more
      integer ncard
c
c  The first set returned is the empty set.
c
      if ( .not. more ) then

        do i = 1, n
          a(i) = 0
        end do

        iadd = 0
        ncard = 0
        more = .true.

      else

        iadd = 1

        if ( mod ( ncard, 2 ) .ne. 0 ) then

10        continue

            iadd = iadd + 1
            if ( a(iadd-1) .ne. 0 ) then
              go to 20
            end if

          go to 10

20        continue

        end if

        a(iadd) = 1 - a(iadd)
        ncard = ncard + 2 * a(iadd) - 1
c
c  The last set returned is the singleton A(N).
c
        if ( ncard .eq. a(n) ) then
          more = .false.
        end if

      end if

      return
      end
      subroutine tetra_07 ( func, x, y, z, result )

c*********************************************************************72
c
cc TETRA_07 approximates an integral inside a tetrahedron in 3D.
c
c  Integration region:
c
c    Points inside a tetrahedron whose four corners are given.
c
c  Discussion:
c
c    A 64 point 7-th degree conical product Gauss formula is used,
c    Stroud number T3:7-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966, pages 42-43,
c    LC: QA299.4G3S7
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision X(4), Y(4), Z(4), the coordinates of 
c    the vertices.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 4 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision quad
      double precision result
      double precision t
      double precision tetra_volume
      double precision u
      double precision v
      double precision volume
      double precision w
      double precision weight1(order)
      double precision weight2(order)
      double precision weight3(order)
      double precision x(4)
      double precision xtab1(order)
      double precision xtab2(order)
      double precision xtab3(order)
      double precision xval
      double precision y(4)
      double precision yval
      double precision z(4)
      double precision zval

      save weight2
      save weight3
      save xtab2
      save xtab3

      data weight2 /
     &  0.1355069134D+00, 0.2034645680D+00, 0.1298475476D+00, 
     &  0.0311809709D+00 /
      data weight3 /
     &  0.1108884156D+00, 0.1434587898D+00, 0.0686338872D+00, 
     &  0.0103522407D+00 /
      data xtab2 /
     &  0.0571041961D+00, 0.2768430136D+00, 0.5835904324D+00, 
     &  0.8602401357D+00 /
      data xtab3 /
     &  0.0485005495D+00, 0.2386007376D+00, 0.5170472951D+00, 
     &  0.7958514179D+00 /
c
c  Get the Gauss-Legendre weights and abscissas for [-1,1].
c
      call legendre_set ( order, xtab1, weight1 )
c
c  Adjust the rule for the interval [0,1].
c
      a = -1.0D+00
      b = +1.0D+00

      c =  0.0D+00
      d =  1.0D+00

      call rule_adjust ( a, b, c, d, order, xtab1, weight1 )
c
c  Carry out the quadrature.
c
      quad = 0.0D+00

      do i = 1, order
        do j = 1, order
          do k = 1, order
c
c  Compute the barycentric coordinates of the point in the unit triangle.
c
            t = xtab3(k)
            u = xtab2(j) * ( 1.0D+00 - xtab3(k) )
            v = xtab1(i) * ( 1.0D+00 - xtab2(j) ) 
     &        * ( 1.0D+00 - xtab3(k) )
            w = 1.0D+00 - t - u - v
c
c  Compute the corresponding point in the triangle.
c
            xval = t * x(1) + u * x(2) + v * x(3) + w * x(4)
            yval = t * y(1) + u * y(2) + v * y(3) + w * y(4)
            zval = t * z(1) + u * z(2) + v * z(3) + w * z(4)

            quad = quad + 6.0D+00 * weight1(i) * weight2(j) 
     &        * weight3(k) * func ( xval, yval, zval )

          end do
        end do
      end do

      volume = tetra_volume ( x, y, z )
      result = quad * volume

      return
      end
      subroutine tetra_sum ( func, x, y, z, order, xtab, ytab, ztab, 
     &  weight, result )

c*********************************************************************72
c
cc TETRA_SUM carries out a quadrature rule in a tetrahedron in 3D.
c
c  Integration region:
c
c    A tetrahedron whose vertices are specified.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, name of the function, of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision X(4), Y(4), Z(4), the vertices.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER), the
c    abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision tetra_volume
      double precision volume
      double precision weight(order)
      double precision x(4)
      double precision xtab(order)
      double precision xval
      double precision y(4)
      double precision ytab(order)
      double precision yval
      double precision z(4)
      double precision ztab(order)
      double precision zval

      quad = 0.0D+00

      do i = 1, order

        xval =             xtab(i)                       * x(1) 
     &                             + ytab(i)             * x(2) 
     &                                       + ztab(i)   * x(3) 
     &       + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * x(4)

        yval =             xtab(i)                       * y(1) 
     &                             + ytab(i)             * y(2) 
     &                                       + ztab(i)   * y(3) 
     &       + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * y(4)

        zval =             xtab(i)                       * z(1) 
     &                             + ytab(i)             * z(2) 
     &                                       + ztab(i)   * z(3) 
     &       + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * z(4)

        quad = quad + weight(i) * func ( xval, yval, zval )

      end do

      volume = tetra_volume ( x, y, z )
      result = quad * volume

      return
      end
      subroutine tetra_tproduct ( func, order, x, y, z, result )

c*********************************************************************72
c
cc TETRA_TPRODUCT approximates an integral in a tetrahedron in 3D.
c
c  Discussion:
c
c    Integration is carried out over the points inside an arbitrary
c    tetrahedron whose four vertices are given.
c
c    An ORDER**3 point (2*ORDER-1)-th degree triangular product
c    Gauss-Legendre rule is used.
c
c    With ORDER = 8, this routine is equivalent to the routine TETR15
c    in the reference, page 367.
c
c    Thanks to Joerg Behrens, jbehren@gwdg.de, for numerous suggestions
c    and corrections.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, integer ORDER, the order of the basic quadrature rules.
c    ORDER should be between 1 and 9.
c
c    Input, double precision X(4), Y(4), Z(4), the vertices
c    of the tetrahedron.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      double precision func
      external func
      integer i
      integer j
      integer k
      double precision quad
      double precision result
      double precision tetra_volume
      double precision volume
      double precision weight0(order)
      double precision weight1(order)
      double precision weight2(order)
      double precision x(4)
      double precision xtab0(order)
      double precision xtab1(order)
      double precision xtab2(order)
      double precision xval
      double precision y(4)
      double precision yval
      double precision z(4)
      double precision zval

      if ( order < 1 .or. 9 < order ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TETRA_TPRODUCT - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The quadrature rule orders must be between 1 and 9.'
        write ( *, '(a,i8)' ) '  The input value was ORDER = ', order
        stop
      end if
c
c  Get the Gauss-Legendre ORDER point rules on [-1,1] for integrating
c    F(X),
c    X * F(X),
c    X * X * F(X).
c
      call legendre_set ( order, xtab0, weight0 )
      call legendre_set_x1 ( order, xtab1, weight1 )
      call legendre_set_x2 ( order, xtab2, weight2 )
c
c  Adjust the rules from [-1,1] to [0,1].
c
      a = -1.0D+00
      b = +1.0D+00
      c =  0.0D+00
      d =  1.0D+00

      call rule_adjust ( a, b, c, d, order, xtab0, weight0 )

      call rule_adjust ( a, b, c, d, order, xtab1, weight1 )

      call rule_adjust ( a, b, c, d, order, xtab2, weight2 )
c
c  For rules with a weight function that is not 1, the weight vectors
c  require further adjustment.
c
      do i = 1, order
        weight1(i) = weight1(i) / 2.0D+00
        weight2(i) = weight2(i) / 4.0D+00
      end do
c
c  Carry out the quadrature.
c
      quad = 0.0D+00

      do k = 1, order
        do j = 1, order
          do i = 1, order

            xval = x(1) + ( ( ( x(4) - x(3) )   * xtab0(i) 
     &                      + ( x(3) - x(2) ) ) * xtab1(j) 
     &                      + ( x(2) - x(1) ) ) * xtab2(k)

            yval = y(1) + ( ( ( y(4) - y(3) )   * xtab0(i) 
     &                      + ( y(3) - y(2) ) ) * xtab1(j) 
     &                      + ( y(2) - y(1) ) ) * xtab2(k)

            zval = z(1) + ( ( ( z(4) - z(3) )   * xtab0(i) 
     &                      + ( z(3) - z(2) ) ) * xtab1(j) 
     &                      + ( z(2) - z(1) ) ) * xtab2(k)

            quad = quad + 6.0D+00 * weight0(i) * weight1(j) 
     &        * weight2(k) * func ( xval, yval, zval )

          end do

        end do

      end do
c
c  Compute the volume of the tetrahedron.
c
      volume = tetra_volume ( x, y, z )
      result = quad * volume

      return
      end
      subroutine tetra_unit_set ( rule, order, xtab, ytab, ztab, 
     &  wtab )

c*********************************************************************72
c
cc TETRA_UNIT_SET sets quadrature weights and abscissas in the unit tetrahedron.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y,
c    and
c      0 <= Z, 
c    and
c      X + Y + Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hermann Engels,
c    Numerical Quadrature and Cubature,
c    Academic Press, 1980,
c    ISBN: 012238850X,
c    LC: QA299.3E5.
c
c    Patrick Keast,
c    Moderate Degree Tetrahedral Quadrature Formulas,
c    Computer Methods in Applied Mechanics and Engineering,
c    Volume 55, Number 3, May 1986, pages 339-348.
c
c    Olgierd Zienkiewicz,
c    The Finite Element Method,
c    Sixth Edition,
c    Butterworth-Heinemann, 2005,
c    ISBN: 0750663200,
c    LC: TA640.2.Z54
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
c     2, order 4, precision 1, Newton Cotes formula #1.
c     3, order 4, precision 2, Zienkiewicz #2.
c     4, order 10, precision 2, Newton Cotes formula #2
c     5, order 5, precision 3, Zienkiewicz #3.
c     6, order 8, precision 3, Newton Cotes formula #3.
c     7, order 35, precision 4, Newton Cotes formula #4.
c     8, order 11, precision 4, a Keast rule.
c
c    Input, integer ORDER, the order of the rule.
c
c    Output, double precision XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER),
c    the abscissas.
c
c    Output, double precision WTAB(ORDER), the weights.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision f
      double precision g
      double precision h
      integer rule
      double precision wtab(order)
      double precision xtab(order)
      double precision ytab(order)
      double precision z
      double precision ztab(order)
c
c  Newton Cotes #0.
c
      if ( rule .eq. 1 ) then

        xtab(1) = 1.0D+00 / 4.0D+00
        ytab(1) = 1.0D+00 / 4.0D+00
        ztab(1) = 1.0D+00 / 4.0D+00
        wtab(1) = 1.0D+00
c
c  Newton Cotes #1.
c
      else if ( rule .eq. 2 ) then

        a = 1.0D+00
        b = 1.0D+00 / 4.0D+00
        z = 0.0D+00

        xtab(1) = z
        xtab(2) = a
        xtab(3) = z
        xtab(4) = z
        ytab(1) = z
        ytab(2) = z
        ytab(3) = a
        ytab(4) = z
        ztab(1) = z
        ztab(2) = z
        ztab(3) = z
        ztab(4) = a
        wtab(1) = b
        wtab(2) = b
        wtab(3) = b
        wtab(4) = b
c
c  Zienkiewicz #2.
c
      else if ( rule .eq. 3 ) then

        a =  0.5854101966249685D+00
        b =  0.1381966011250105D+00
        c =  0.25D+00
 
        xtab(1) = a
        xtab(2) = b
        xtab(3) = b
        xtab(4) = b
        ytab(1) = b
        ytab(2) = a
        ytab(3) = b
        ytab(4) = b
        ztab(1) = b
        ztab(2) = b
        ztab(3) = a
        ztab(4) = b
        wtab(1) = c
        wtab(2) = c
        wtab(3) = c
        wtab(4) = c
c
c  Newton Cotes #2.
c
      else if ( rule .eq. 4 ) then

        a =  1.0D+00
        b =  0.5D+00
        c = -1.0D+00 / 20.0D+00
        d =  4.0D+00 / 20.0D+00
        z =  0.0D+00

        xtab(1) = z
        xtab(2) = a
        xtab(3) = z
        xtab(4) = z
        xtab(5) = b
        xtab(6) = z
        xtab(7) = z
        xtab(8) = b
        xtab(9) = b
        xtab(10) = z
        ytab(1) = z
        ytab(2) = z
        ytab(3) = a
        ytab(4) = z
        ytab(5) = z
        ytab(6) = b
        ytab(7) = z
        ytab(8) = b
        ytab(9) = z
        ytab(10) = b
        ztab(1) = z
        ztab(2) = z
        ztab(3) = z
        ztab(4) = a
        ztab(5) = z
        ztab(6) = z
        ztab(7) = b
        ztab(8) = z
        ztab(9) = b
        ztab(10) = b
        wtab(1) = c
        wtab(2) = c
        wtab(3) = c
        wtab(4) = c
        wtab(5) = d
        wtab(6) = d
        wtab(7) = d
        wtab(8) = d
        wtab(9) = d 
        wtab(10) = d
c
c  Zienkiewicz #3.
c
      else if ( rule .eq. 5 ) then

        a =  1.0D+00 / 6.0D+00
        b =  0.25D+00
        c =  0.5D+00
        d = -0.8D+00
        e =  0.45D+00

        xtab(1) = b
        xtab(2) = c
        xtab(3) = a
        xtab(4) = a
        xtab(5) = a
        ytab(1) = b
        ytab(2) = a
        ytab(3) = c
        ytab(4) = a
        ytab(5) = a
        ztab(1) = b
        ztab(2) = a
        ztab(3) = a
        ztab(4) = c
        ztab(5) = a
        wtab(1) = d
        wtab(2) = e
        wtab(3) = e
        wtab(4) = e
        wtab(5) = e
c
c  Newton Cotes #3.
c  (This is actually formally a 20 point rule, but with 12 zero coefficients!)
c
      else if ( rule .eq. 6 ) then

        a = 1.0D+00
        b = 1.0D+00 / 40.0D+00
        c = 1.0D+00 /  3.0D+00
        d = 9.0D+00 / 40.0D+00
        z = 0.0D+00

        xtab(1) = z
        xtab(2) = a
        xtab(3) = z
        xtab(4) = z
        xtab(5) = c
        xtab(6) = c
        xtab(7) = z
        xtab(8) = c
        ytab(1) = z
        ytab(2) = z
        ytab(3) = a
        ytab(4) = z
        ytab(5) = c
        ytab(6) = z
        ytab(7) = c
        ytab(8) = c
        ztab(1) = z
        ztab(2) = z
        ztab(3) = z
        ztab(4) = a
        ztab(5) = z
        ztab(6) = c
        ztab(7) = c
        ztab(8) = c
        wtab(1) = b
        wtab(2) = b
        wtab(3) = b
        wtab(4) = b
        wtab(5) = d
        wtab(6) = d
        wtab(7) = d
        wtab(8) = d
c
c  Newton Cotes #4.
c
      else if ( rule .eq. 7 ) then

        a =   0.25D+00
        b =   0.50D+00
        c =   0.75D+00
        d =   1.00D+00
        e =  -5.0D+00 / 420.0D+00
        f = -12.0D+00 / 420.0D+00
        g =  16.0D+00 / 420.0D+00
        h = 128.0D+00 / 420.0D+00
        z =   0.0D+00

        xtab(1) = z
        xtab(2) = d
        xtab(3) = z
        xtab(4) = z
        xtab(5) = a
        xtab(6) = z
        xtab(7) = z
        xtab(8) = c
        xtab(9) = c
        xtab(10) = c
        xtab(11) = z
        xtab(12) = a
        xtab(13) = z
        xtab(14) = z
        xtab(15) = a
        xtab(16) = z
        xtab(17) = b
        xtab(18) = z
        xtab(19) = z
        xtab(20) = b
        xtab(21) = b
        xtab(22) = z
        xtab(23) = a
        xtab(24) = b
        xtab(25) = a
        xtab(26) = a
        xtab(27) = b
        xtab(28) = z
        xtab(29) = b
        xtab(30) = z
        xtab(31) = a
        xtab(32) = a
        xtab(33) = z
        xtab(34) = a
        xtab(35) = a
        ytab(1) = z
        ytab(2) = z
        ytab(3) = d
        ytab(4) = z
        ytab(5) = z
        ytab(6) = a
        ytab(7) = z
        ytab(8) = z
        ytab(9) = a
        ytab(10) = z
        ytab(11) = c
        ytab(12) = c
        ytab(13) = c
        ytab(14) = z
        ytab(15) = z
        ytab(16) = a
        ytab(17) = z
        ytab(18) = b
        ytab(19) = z
        ytab(20) = b
        ytab(21) = z
        ytab(22) = b
        ytab(23) = a
        ytab(24) = a
        ytab(25) = b
        ytab(26) = z
        ytab(27) = z
        ytab(28) = a
        ytab(29) = a
        ytab(30) = b
        ytab(31) = b
        ytab(32) = z
        ytab(33) = a
        ytab(34) = a
        ytab(35) = a
        ztab(1) = z
        ztab(2) = z
        ztab(3) = z
        ztab(4) = d
        ztab(5) = z
        ztab(6) = z
        ztab(7) = a
        ztab(8) = z
        ztab(9) = z
        ztab(10) = a
        ztab(11) = z
        ztab(12) = z
        ztab(13) = a
        ztab(14) = c
        ztab(15) = c
        ztab(16) = c
        ztab(17) = z
        ztab(18) = z
        ztab(19) = b
        ztab(20) = z
        ztab(21) = b
        ztab(22) = b
        ztab(23) = z
        ztab(24) = z
        ztab(25) = z
        ztab(26) = a
        ztab(27) = a
        ztab(28) = a
        ztab(29) = a
        ztab(30) = a
        ztab(31) = a
        ztab(32) = b
        ztab(33) = b
        ztab(34) = b
        ztab(35) = a
        wtab(1) = e
        wtab(2) = e
        wtab(3) = e
        wtab(4) = e
        wtab(5) = g
        wtab(6) = g
        wtab(7) = g
        wtab(8) = g
        wtab(9) = g
        wtab(10) = g
        wtab(11) = g
        wtab(12) = g
        wtab(13) = g
        wtab(14) = g
        wtab(15) = g
        wtab(16) = g
        wtab(17) = f
        wtab(18) = f
        wtab(19) = f
        wtab(20) = f
        wtab(21) = f
        wtab(22) = f
        wtab(23) = g
        wtab(24) = g
        wtab(25) = g
        wtab(26) = g
        wtab(27) = g
        wtab(28) = g
        wtab(29) = g
        wtab(30) = g
        wtab(31) = g
        wtab(32) = g
        wtab(33) = g
        wtab(34) = g
        wtab(35) = h
c
c  Keast Rule of order 11
c
      else if ( rule .eq. 8 ) then

        a =  0.25D+00
        b =  11.0D+00 /    14.0D+00
        c =   1.0D+00 /    14.0D+00
        d =  0.25D+00 * ( 1.0D+00 + sqrt ( 5.0D+00 / 14.0D+00 ) )
        e =  0.25D+00 * ( 1.0D+00 - sqrt ( 5.0D+00 / 14.0D+00 ) )
        f = -74.0D+00 /  5625.0D+00
        g = 343.0D+00 / 45000.0D+00
        h =  56.0D+00 /  2250.0D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = c
        xtab(6) = d
        xtab(7) = d
        xtab(8) = d
        xtab(9) = e
        xtab(10) = e
        xtab(11) = e
        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = c
        ytab(6) = d
        ytab(7) = e
        ytab(8) = e
        ytab(9) = d
        ytab(10) = d
        ytab(11) = e
        ztab(1) = a
        ztab(2) = c
        ztab(3) = c
        ztab(4) = b
        ztab(5) = c
        ztab(6) = e
        ztab(7) = d
        ztab(8) = e
        ztab(9) = d
        ztab(10) = e
        ztab(11) = d
        wtab(1) = f
        wtab(2) = g
        wtab(3) = g
        wtab(4) = g
        wtab(5) = g
        wtab(6) = h
        wtab(7) = h
        wtab(8) = h
        wtab(9) = h
        wtab(10) = h
        wtab(11) = h

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TETRA_UNIT_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
        stop

      end if

      return
      end
      subroutine tetra_unit_size ( rule, order )

c*********************************************************************72
c
cc TETRA_UNIT_SIZE sizes quadrature rules in the unit tetrahedron.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y,
c    and
c      0 <= Z, 
c    and
c      X + Y + Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hermann Engels,
c    Numerical Quadrature and Cubature,
c    Academic Press, 1980,
c    ISBN: 012238850X,
c    LC: QA299.3E5.
c
c    Patrick Keast,
c    Moderate Degree Tetrahedral Quadrature Formulas,
c    Computer Methods in Applied Mechanics and Engineering,
c    Volume 55, Number 3, May 1986, pages 339-348.
c
c    Olgierd Zienkiewicz,
c    The Finite Element Method,
c    Sixth Edition,
c    Butterworth-Heinemann, 2005,
c    ISBN: 0750663200,
c    LC: TA640.2.Z54
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
c     2, order 4, precision 1, Newton Cotes formula #1.
c     3, order 4, precision 2, Zienkiewicz #2.
c     4, order 10, precision 2, Newton Cotes formula #2
c     5, order 5, precision 3, Zienkiewicz #3.
c     6, order 8, precision 3, Newton Cotes formula #3.
c     7, order 35, precision 4, Newton Cotes formula #4.
c     8, order 11, precision 4, a Keast rule.
c
c    Output, integer ORDER, the order of the rule.
c
      implicit none

      integer order
      integer rule
c
c  Newton Cotes #0.
c
      if ( rule .eq. 1 ) then

        order = 1
c
c  Newton Cotes #1.
c
      else if ( rule .eq. 2 ) then

        order = 4
c
c  Zienkiewicz #2.
c
      else if ( rule .eq. 3 ) then

        order = 4
c
c  Newton Cotes #2.
c
      else if ( rule .eq. 4 ) then

        order = 10
c
c  Zienkiewicz #3.
c
      else if ( rule .eq. 5 ) then

        order = 5
c
c  Newton Cotes #3.
c  (This is actually formally a 20 point rule, but with 12 zero coefficients!)
c
      else if ( rule .eq. 6 ) then

        order = 8
c
c  Newton Cotes #4.
c
      else if ( rule .eq. 7 ) then

        order = 35
c
c  Keast Rule of order 11
c
      else if ( rule .eq. 8 ) then

        order = 11

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TETRA_UNIT_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
        stop

      end if

      return
      end
      subroutine tetra_unit_sum ( func, order, xtab, ytab, ztab, 
     &  weight, result )

c*********************************************************************72
c
cc TETRA_UNIT_SUM carries out a quadrature rule in the unit tetrahedron in 3D.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y,
c    and
c      0 <= Z, 
c    and
c      X + Y + Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER), the
c    abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision tetra_unit_volume
      double precision volume
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)
      double precision ztab(order)

      quad = 0.0D+00

      do i = 1, order
        quad = quad + weight(i) * func ( xtab(i), ytab(i), ztab(i) )
      end do

      volume = tetra_unit_volume ( )
      result = quad * volume

      return
      end
      function tetra_unit_volume ( )

c*********************************************************************72
c
cc TETRA_UNIT_VOLUME returns the volume of the unit tetrahedron in 3D.
c
c  Discussion:
c
c    The integration region is:
c
c      0 <= X,
c      0 <= Y,
c      0 <= Z, 
c      X + Y + Z <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision TETRA_UNIT_VOLUME, the volume.
c
      implicit none

      double precision tetra_unit_volume

      tetra_unit_volume = 1.0D+00 / 6.0D+00

      return
      end
      function tetra_volume ( x, y, z )

c*********************************************************************72
c
cc TETRA_VOLUME computes the volume of a tetrahedron in 3D.
c
c  Integration region:
c
c    Points inside a tetrahedron whose four vertices are given.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(4), Y(4), Z(4), the vertices.
c
c    Output, double precision TETRA_VOLUME, the volume of the tetrahedron.
c
      implicit none

      double precision parallelipiped_volume_3d
      double precision tetra_unit_volume
      double precision tetra_volume
      double precision volume
      double precision x(4)
      double precision y(4)
      double precision z(4)

      volume = parallelipiped_volume_3d ( x, y, z )

      tetra_volume = volume * tetra_unit_volume ( )

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
      subroutine torus_1 ( func, r1, r2, n, result )

c*********************************************************************72
c
cc TORUS_1 approximates an integral on the surface of a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
c
c  Discussion:
c
c    An (N+1)*(N+2) point N-th degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Input, integer N, defines the degree of the formula
c    used to approximate the integral.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision angle
      double precision ct1
      double precision func
      external func
      integer i
      integer j
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision st1
      double precision torus_area_3d
      double precision u
      double precision volume
      double precision w
      double precision x
      double precision y
      double precision z

      w = 1.0D+00 / ( r1 * dble ( ( n + 1 ) * ( n + 2 ) ) )
      quad = 0.0D+00

      do i = 1, n + 1

        angle = 2.0D+00 * pi * dble ( i ) / dble ( n + 1 )
        ct1 = cos ( angle )
        st1 = sin ( angle )

        do j = 1, n + 2

          angle = 2.0D+00 * pi * dble ( j ) / dble ( n + 2 )
          u = r1 + r2 * cos ( angle )
          x = u * ct1
          y = u * st1
          z = r2 * sin ( angle )

          quad = quad + w * u * func ( x, y, z )

        end do

      end do

      volume = torus_area_3d ( r1, r2 )
      result = quad * volume

      return
      end
      subroutine torus_14s ( func, r1, r2, result )

c*********************************************************************72
c
cc TORUS_14S approximates an integral inside a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
c
c  Discussion:
c
c    A 960 point 14-th degree formula is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 4 )

      double precision angle
      double precision ct
      double precision cth
      double precision func
      external func
      integer i
      integer j
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r(order)
      double precision r1
      double precision r2
      double precision result
      double precision st
      double precision sth
      double precision torus_volume_3d
      double precision u
      double precision volume
      double precision weight(order)
      double precision x
      double precision y
      double precision z

      save r
      save weight

      data r /
     &  0.263499230D+00, 0.574464514D+00, 0.818529487D+00, 
     &  0.964659606D+00 /
      data weight /
     &  0.086963711D+00, 0.163036289D+00, 0.163036289D+00, 
     &  0.086963711D+00 /

      quad = 0.0D+00

      do n = 1, 15

        angle = 2.0D+00 * pi * dble ( n ) / 15.0D+00
        cth = cos ( angle )
        sth = sin ( angle )

        do i = 1, 16

          angle = 2.0D+00 * pi * dble ( i ) / 16.0D+00
          ct = cos ( angle )
          st = sin ( angle )

          do j = 1, order
            u = r1 + r(j) * ct * r2
            x = u * cth
            y = u * sth
            z = r(j) * st * r2
            quad = quad + u * weight(j) * func ( x, y, z ) 
     &        / ( 120.0D+00 * r1 )
          end do

        end do

      end do

      volume = torus_volume_3d ( r1, r2 )
      result = quad * volume

      return
      end
      subroutine torus_5s2 ( func, r1, r2, result )

c*********************************************************************72
c
cc TORUS_5S2 approximates an integral inside a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
c
c  Discussion:
c
c    A 24 point, 5-th degree formula is used, Stroud number TOR3-S2:5-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision angle
      double precision cs
      double precision func
      external func
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision sn
      double precision torus_volume_3d
      double precision u1
      double precision u2
      double precision u3
      double precision volume
      double precision w
      double precision x
      double precision y
      double precision z

      w = 1.0D+00 / 24.0D+00

      quad = 0.0D+00

      u1 = sqrt ( r1 * r1 + 0.5D+00 * r2 * r2 )
      u2 = sqrt ( r1 * r1 + sqrt ( 2.0D+00 ) * r1 * r2 + r2 * r2 )
      u3 = sqrt ( r1 * r1 - sqrt ( 2.0D+00 ) * r1 * r2 + r2 * r2 )

      do i = 1, 6

        angle = 2.0D+00 * pi * dble ( i ) / 6.0D+00
        cs = cos ( angle )
        sn = sin ( angle )

        x = u1 * cs
        y = u1 * sn
        z = r2 / sqrt ( 2.0D+00 )
        quad = quad + w * func ( x, y, z )

        x = u1 * cs
        y = u1 * sn
        z = -r2 / sqrt ( 2.0D+00 )
        quad = quad + w * func ( x, y, z )

        x = u2 * cs
        y = u2 * sn
        z = 0.0D+00
        quad = quad + w * func ( x, y, z )

        x = u3 * cs
        y = u3 * sn
        z = 0.0D+00
        quad = quad + w * func ( x, y, z )

      end do

      volume = torus_volume_3d ( r1, r2 )
      result = quad * volume

      return
      end
      subroutine torus_6s2 ( func, r1, r2, result )

c*********************************************************************72
c
cc TORUS_6S2 approximates an integral inside a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
c
c  Discussion:
c
c    An 84 point 6-th degree formula is used, Stroud number TOR3-S2:6-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 2 )

      double precision cth
      double precision func
      external func
      integer i
      integer j
      integer k
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision s(order)
      double precision sth
      double precision torus_volume_3d
      double precision u
      double precision v
      double precision volume
      double precision w
      double precision weight(order)
      double precision x
      double precision y
      double precision z

      save s
      save weight

      data s / 0.322914992D+00, 0.644171310D+00 /
      data weight / 0.387077796D+00, 0.165609800D+00 /

      w = 1.0D+00 / ( 7.0D+00 * r1 * pi )

      quad = 0.0D+00

      do n = 1, 7

        u = 0.5D+00 * sqrt ( 3.0D+00 ) * r2
        cth = cos ( 2.0D+00 * pi * dble ( n ) / 7.0D+00 )
        sth = sin ( 2.0D+00 * pi * dble ( n ) / 7.0D+00 )

        do i = 1, 2

          u = -u

          x = ( r1 + u ) * cth
          y = ( r1 + u ) * sth
          z = 0.0D+00
          quad = quad + 0.232710567D+00 * w * ( r1 + u ) 
     &      * func ( x, y, z )

          x = r1 * cth
          y = r1 * sth
          z = u
          quad = quad + 0.232710567D+00 * w * r1 * func ( x, y, z )

        end do

        do k = 1, order

          u = s(k) * r2
          v = u

          do i = 1, 2

            u = -u

            do j = 1, 2

              v = -v

              x = ( r1 + u ) * cth
              y = ( r1 + u ) * sth
              z = v
              quad = quad + weight(k) * w * ( r1 + u ) 
     &          * func ( x, y, z )

            end do
          end do
        end do
      end do

      volume = torus_volume_3d ( r1, r2 )
      result = quad * volume

      return
      end
      function torus_area_3d ( r1, r2 )

c*********************************************************************72
c
cc TORUS_AREA_3D returns the area of a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision TORUS_AREA_3D, the area of the torus.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision torus_area_3d

      torus_area_3d = 4.0D+00 * pi * pi * r1 * r2

      return
      end
      subroutine torus_square_14c ( func, r1, r2, result )

c*********************************************************************72
c
cc TORUS_SQUARE_14C approximates an integral in a "square" torus in 3D.
c
c  Discussion:
c
c    A 14-th degree 960 point formula is used.
c
c  Integration region:
c
c      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
c    and
c      -R2 <= Z <= R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated, of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, R2, the radii that define the torus.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order
      parameter ( order = 8 )

      double precision angle
      double precision cth
      double precision func
      external func
      integer i
      integer j
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision rtab(order)
      double precision sth
      double precision torus_square_volume_3d
      double precision u
      double precision volume
      double precision w
      double precision weight(order)
      double precision x
      double precision y
      double precision z

      call legendre_set ( order, rtab, weight )

      w = 1.0D+00 / ( 60.0D+00 * r1 )
      quad = 0.0D+00

      do n = 1, 15

        angle = 2.0D+00 * pi * dble ( n ) / 15.0D+00
        cth = cos ( angle )
        sth = sin ( angle )

        do i = 1, order

          u = r1 + rtab(i) * r2
          x = u * cth
          y = u * sth

          do j = 1, order
            z = rtab(j) * r2
            quad = quad + u * w * weight(i) * weight(j) 
     &        * func ( x, y, z )
          end do

        end do

      end do

      volume = torus_square_volume_3d ( r1, r2 )
      result = quad * volume

      return
      end
      subroutine torus_square_5c2 ( func, r1, r2, result )

c*********************************************************************72
c
cc TORUS_SQUARE_5C2 approximates an integral in a "square" torus in 3D.
c
c  Integration region:
c
c      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
c    and
c      -R2 <= Z <= R2.
c
c  Discussion:
c
c    A 24 point 5-th degree formula is used, Stroud number TOR3-C2:5-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of three variables which is to be integrated,
c    of the form:
c      function func ( x, y, z )
c      double precision func
c      double precision x
c      double precision y
c      double precision z
c
c    Input, double precision R1, the primary radius of the torus.
c
c    Input, double precision R2, one-half the length of a side of the
c    square cross-section.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      double precision b1 
      parameter ( b1 = 5.0D+00 / 108.0D+00 )
      double precision b2
      parameter ( b2 = 4.0D+00 / 108.0D+00 )
      double precision cs
      double precision func
      external func
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      double precision r1
      double precision r2
      double precision result
      double precision sn
      double precision torus_square_volume_3d
      double precision u1
      double precision u2
      double precision u3
      double precision v
      double precision volume
      double precision x
      double precision y
      double precision z

      quad = 0.0D+00

      u1 = sqrt ( r1 * r1 + r2 * r2 )

      v = r2 * sqrt ( 0.6D+00 )

      u2 = sqrt ( r1 * r1 - sqrt ( 3.0D+00 ) * r1 * r2 + r2 * r2 )

      u3 = sqrt ( r1 * r1 + sqrt ( 3.0D+00 ) * r1 * r2 + r2 * r2 )

      do i = 1, 6

        cs = cos ( dble ( i ) * pi / 3.0D+00 )
        sn = sin ( dble ( i ) * pi / 3.0D+00 )

        x = u1 * cs
        y = u1 * sn
        z = v
        quad = quad + b1 * func ( x, y, z )

        z = -v
        quad = quad + b1 * func ( x, y, z )

        x = u2 * cs
        y = u2 * sn
        z = 0.0D+00
        quad = quad + b2 * func ( x, y, z )

        x = u3 * cs
        y = u3 * sn
        z = 0.0D+00
        quad = quad + b2 * func ( x, y, z )

      end do

      volume = torus_square_volume_3d ( r1, r2 )
      result = quad * volume

      return
      end
      function torus_square_area_3d ( r1, r2 )

c*********************************************************************72
c
cc TORUS_SQUARE_AREA_3D returns the area of a square torus in 3D.
c
c  Integration region:
c
c      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
c    and
c      -R2 <= Z <= R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision TORUS_SQUARE_AREA_3D, the area of the torus.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision torus_square_area_3d

      torus_square_area_3d = 16.0D+00 * pi * r1 * r2

      return
      end
      function torus_square_volume_3d ( r1, r2 )

c*********************************************************************72
c
cc TORUS_SQUARE_VOLUME_3D returns the volume of a square torus in 3D.
c
c  Integration region:
c
c      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
c    and
c      -R2 <= Z <= R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision TORUS_SQUARE_VOLUME_3D, the volume of the torus.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision torus_square_volume_3d

      torus_square_volume_3d = 8.0D+00 * pi * r1 * r2 * r2

      return
      end
      function torus_volume_3d ( r1, r2 )

c*********************************************************************72
c
cc TORUS_VOLUME_3D returns the volume of a torus in 3D.
c
c  Integration region:
c
c    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the two radii that define the torus.
c
c    Output, double precision TORUS_VOLUME_3D, the volume of the torus.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision torus_volume_3d

      torus_volume_3d = 2.0D+00 * pi * pi * r1 * r2 * r2

      return
      end
      subroutine triangle_rule_adjust ( xval, yval, order, xtab, ytab, 
     &  weight, xtab2, ytab2, weight2 )

c*********************************************************************72
c
cc TRIANGLE_RULE_ADJUST adjusts a unit quadrature rule to an arbitrary triangle.
c
c  Integration region:
c
c      (X,Y) = ALPHA * (X1,Y1) + BETA * (X2,Y2) + ( 1 - ALPHA - BETA ) * (X3,Y3)
c    and
c      0 <= ALPHA <= 1 - BETA
c    and
c      0 <= BETA <= 1 - ALPHA
c
c  Discussion:
c
c    This routine accepts as input abscissas and weights appropriate for
c    quadrature in the unit triangle, and returns abscissas and weights
c    appropriate for quadrature in a given triangle.
c
c    Once this routine has been called, an integral over the given triangle
c    can be approximated as:
c
c      QUAD = sum ( 1 <= I <= ORDER ) WTAB2(I) * FUNC ( XTAB2(I), YTAB2(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision XVAL(3), YVAL(3), the coordinates of the nodes.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas for
c    the unit triangle.
c
c    Input, double precision WEIGHT(ORDER), the weights for the unit triangle.
c
c    Output, double precision XTAB2(ORDER), YTAB2(ORDER), the adjusted
c    abscissas.
c
c    Output, double precision WEIGHT2(ORDER), the adjusted weights.
c
      implicit none

      integer order

      integer i
      double precision triangle_volume
      double precision volume
      double precision weight(order)
      double precision weight2(order)
      double precision xtab(order)
      double precision xtab2(order)
      double precision xval(3)
      double precision ytab(order)
      double precision ytab2(order)
      double precision yval(3)

      volume = triangle_volume ( xval, yval )

      do i = 1, order

        xtab2(i) =             xtab(i)             * xval(1) 
     &           +                       ytab(i)   * xval(2) 
     &           + ( 1.0D+00 - xtab(i) - ytab(i) ) * xval(3)

        ytab2(i) =             xtab(i)             * yval(1) 
     &                                 + ytab(i)   * yval(2) 
     &           + ( 1.0D+00 - xtab(i) - ytab(i) ) * yval(3)

        weight2(i) = weight(i) * 2.0D+00 * volume

      end do

      return
      end
      subroutine triangle_sub ( func, xval, yval, nsub, order, xtab, 
     &  ytab, weight, result )

c*********************************************************************72
c
cc TRIANGLE_SUB carries out quadrature over subdivisions of a triangular region.
c
c  Integration region:
c
c      (X,Y) =       ALPHA          * ( XVAL(1), YVAL(1) )
c            +               BETA   * ( XVAL(2), YVAL(2) )
c            + ( 1 - ALPHA - BETA ) * ( XVAL(3), YVAL(3) )
c    and
c      0 <= ALPHA <= 1 - BETA
c    and
c      0 <= BETA <= 1 - ALPHA
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of
c    two variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision XVAL(3), YVAL(3), the coordinates of the 
c    triangle vertices.
c
c    Input, integer NSUB, the number of subdivisions of each side 
c    of the input triangle to be made.  NSUB = 1 means no subdivisions are made.
c    NSUB = 3 means that each side of the triangle is subdivided into
c    three portions, and that the original triangle is subdivided into
c    NSUB * NSUB triangles.  NSUB must be at least 1.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      integer j
      integer k
      integer nsub
      double precision quad
      double precision result
      double precision temp1
      double precision temp2
      double precision triangle_volume
      double precision volume
      double precision weight(order)
      double precision x
      double precision x1
      double precision x2
      double precision x3
      double precision xtab(order)
      double precision xval(3)
      double precision y
      double precision y1
      double precision y2
      double precision y3
      double precision ytab(order)
      double precision yval(3)
c
c  Initialize RESULT, the approximate integral.
c
      result = 0.0D+00
c
c  NSUB must be positive.
c
      if ( nsub .le. 0 ) then
        return
      end if
c
c  Initialize QUAD, the quadrature sum.
c
      quad = 0.0D+00
c
c  The sub-triangles can be grouped into NSUB strips.
c
      do i = 1, nsub

        temp1 = 0.0D+00
        temp2 = dble ( i ) / dble ( nsub )

        x2 = xval(2) + temp1 * ( xval(3) - xval(2) ) 
     &               + temp2 * ( xval(1) - xval(2) )

        y2 = yval(2) + temp1 * ( yval(3) - yval(2) ) 
     &               + temp2 * ( yval(1) - yval(2) )

        temp1 = 0.0D+00
        temp2 = dble ( i - 1 ) / dble ( nsub )

        x3 = xval(2) + temp1 * ( xval(3) - xval(2) ) 
     &               + temp2 * ( xval(1) - xval(2) )

        y3 = yval(2) + temp1 * ( yval(3) - yval(2) ) 
     &               + temp2 * ( yval(1) - yval(2) )
c
c  There are 2*I-1 triangles in strip number I.
c  The next triangle in the strip shares two nodes with the previous one.
c  Compute its corners, (X1,Y1), (X2,Y2), (X3,Y3).
c
        do j = 1, 2*i-1

          x1 = x2
          y1 = y2
          x2 = x3
          y2 = y3
          temp1 = dble ( ( ( j + 1 ) / 2 ) ) / dble ( nsub )
          temp2 = dble ( ( i - 1 - ( j / 2 ) ) ) / dble ( nsub )

          x3 = xval(2) + temp1 * ( xval(3) - xval(2) ) 
     &                 + temp2 * ( xval(1) - xval(2) )

          y3 = yval(2) + temp1 * ( yval(3) - yval(2) ) 
     &                 + temp2 * ( yval(1) - yval(2) )
c
c  Now integrate over the triangle, mapping the points ( XTAB(K), YTAB(K) )
c  into the triangle.
c
          do k = 1, order

            x = x2 + xtab(k) * ( x3 - x2 ) + ytab(k) * ( x1 - x2 )
            y = y2 + xtab(k) * ( y3 - y2 ) + ytab(k) * ( y1 - y2 )
            quad = quad + weight(k) * func ( x, y )

           end do

        end do

      end do

      volume = triangle_volume ( xval, yval ) / dble ( nsub * nsub )
      result = quad * volume

      return
      end
      subroutine triangle_sum ( func, xval, yval, order, xtab, ytab, 
     &  weight, result )

c*********************************************************************72
c
cc TRIANGLE_SUM carries out a unit quadrature rule in an arbitrary triangle.
c
c  Integration region:
c
c      (X,Y) =       ALPHA          * (X1,Y1) 
c            +               BETA   * (X2,Y2) 
c            + ( 1 - ALPHA - BETA ) * (X3,Y3)
c    and
c      0 <= ALPHA <= 1 - BETA
c    and
c      0 <= BETA <= 1 - ALPHA
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of
c    two variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, double precision XVAL(3), YVAL(3), the coordinates of the nodes.
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision triangle_volume
      double precision volume
      double precision weight(order)
      double precision x
      double precision xtab(order)
      double precision xval(3)
      double precision y
      double precision ytab(order)
      double precision yval(3)

      quad = 0.0D+00

      do i = 1, order

        x =             xtab(i)             * xval(1) 
     &    +                       ytab(i)   * xval(2) 
     &    + ( 1.0D+00 - xtab(i) - ytab(i) ) * xval(3)

        y =             xtab(i)             * yval(1) 
     &    +                       ytab(i)   * yval(2) 
     &    + ( 1.0D+00 - xtab(i) - ytab(i) ) * yval(3)

        quad = quad + weight(i) * func ( x, y )

      end do

      volume = triangle_volume ( xval, yval )
      result = quad * volume

      return
      end
      subroutine triangle_sum_adjusted ( func, order, xtab, ytab, 
     &  weight, result )

c*********************************************************************72
c
cc TRIANGLE_SUM_ADJUSTED carries out an adjusted quadrature rule in a triangle.
c
c  Integration region:
c
c      (X,Y) =       ALPHA          * (X1,Y1) 
c                          + BETA   * (X2,Y2) 
c            + ( 1 - ALPHA - BETA ) * (X3,Y3)
c    and
c      0 <= ALPHA <= 1 - BETA
c    and
c      0 <= BETA <= 1 - ALPHA
c
c  Discussion:
c
c    It is assumed that a quadrature rule approprate for the unit triangle
c    was generated, and then adjusted to a particular triangle by calling
c    TRIANGLE_RULE_ADJUST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied function of
c    two variables which is to be integrated, of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision result
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)

      result = 0.0D+00

      do i = 1, order
        result = result + weight(i) * func ( xtab(i), ytab(i) )
      end do

      return
      end
      subroutine triangle_unit_product_set ( rule, order, xtab, ytab, 
     &  weight )

c*********************************************************************72
c
cc TRIANGLE_UNIT_PRODUCT_SET sets a product rule on the unit triangle.
c
c  Discussion:
c
c    For a given order of accuracy, a product rule on a triangle usually
c    uses more points than necessary.  That is, there is usually a rule
c    of the same order that uses fewer points.
c
c    However, one advantage of product rules is that a rule of any
c    desired order can be generated automatically.
c   
c    The integration region is:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the order of the 1D rule.
c
c    Input, integer ORDER, the order of the rule.
c
c    Output, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Output, double precision WEIGHT(ORDER), the weights of the rule.
c
      implicit none

      integer order
      integer rule

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      integer j
      integer k
      integer order0
      integer order1
      double precision weight(order)
      double precision weight0(rule)
      double precision weight1(rule)
      double precision xtab(order)
      double precision xtab0(rule)
      double precision xtab1(rule)
      double precision ytab(order)

      a = -1.0D+00
      b = +1.0D+00
      c =  0.0D+00
      d = +1.0D+00

      order0 = rule
      call legendre_set ( order0, xtab0, weight0 )
      call rule_adjust ( a, b, c, d, order0, xtab0, weight0 )

      order1 = rule
      call legendre_set_x1 ( order1, xtab1, weight1 )
      call rule_adjust ( a, b, c, d, order1, xtab1, weight1 )

      k = 0
      do j = 1, order1
        do i = 1, order0
          k = k + 1
          xtab(k) = 1.0D+00 - xtab1(j)
          ytab(k) = xtab0(i) * xtab1(j)
          weight(k) = weight0(i) * weight1(j)
        end do
      end do

      return
      end
      subroutine triangle_unit_product_size ( rule, order )

c*********************************************************************72
c
cc TRIANGLE_UNIT_PRODUCT_SIZE sizes a product rule on the unit triangle.
c
c  Discussion:
c
c    The integration region is:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c  Parameters:
c
c    Input, integer RULE, the order of the 1D rule.
c
c    Output, integer ORDER, the order of the rule. 
c
      implicit none

      integer order
      integer rule

      order = rule * rule

      return
      end
      subroutine triangle_unit_set ( rule, order, xtab, ytab, weight )

c*********************************************************************72
c
cc TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
c
c  Discussion:
c
c    The user is responsible for determining the value of ORDER,
c    and appropriately dimensioning the arrays XTAB, YTAB and
c    WEIGHT so that they can accommodate the data.
c
c    The value of ORDER for each rule can be found by invoking
c    the function TRIANGLE_RULE_SIZE.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Graph:
c
c      ^
c    1 | *
c      | |\
c    Y | | \
c      | |  \
c    0 | *---*
c      +------->
c        0 X 1
c
c   The rules are accessed by an index number, RULE.  The indices,
c   and the descriptions of the corresponding rules, are:
c
c     1, ORDER =  1, precision 1, Zienkiewicz #1.
c     2, ORDER =  2, precision 1, (the "vertex rule").
c     3, ORDER =  3, precision 2, Strang and Fix formula #1.
c     4, ORDER =  3, precision 2, Strang and Fix formula #2,
c                                 Zienkiewicz #2.
c     5, ORDER =  4, precision 3, Strang and Fix formula #3,
c                                 Zienkiewicz #3.
c     6, ORDER =  6, precision 3, Strang and Fix formula #4.
c     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
c     8, ORDER =  6, precision 4, Strang and Fix formula #5.
c     9, ORDER =  7, precision 4, Strang and Fix formula #6.
c    10, ORDER =  7, precision 5, Strang and Fix formula #7,
c                                 Stroud formula T2:5-1, 
c                                 Zienkiewicz #4, 
c                                 Schwarz Table 2.2.
c    11, ORDER =  9, precision 6, Strang and Fix formula #8.
c    12, ORDER = 12, precision 6, Strang and Fix formula #9.
c    13, ORDER = 13, precision 7, Strang and Fix formula #10.
c        Note that there is a typographical error in Strang and Fix
c        which lists the value of the XSI(3) component of the
c        last generator point as 0.4869... when it should be 0.04869...
c    14, ORDER =  7, precision 3.
c    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
c    16, ORDER = 64, precision 15, triangular product Gauss rule.
c    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
c    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
c    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
c    20, ORDER = 37, precision 13, from ACM TOMS #706.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706,
c    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c    Dirk Laurie,
c    Algorithm 584,
c    CUBTRI, Automatic Cubature Over a Triangle,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 2, 1982, pages 210-218.
c
c    James Lyness, Dennis Jespersen,
c    Moderate Degree Symmetric Quadrature Rules for the Triangle,
c    Journal of the Institute of Mathematics and its Applications,
c    Volume 15, Number 1, February 1975, pages 19-32.
c
c    Hans Rudolf Schwarz,
c    Finite Element Methods,
c    Academic Press, 1988,
c    ISBN: 0126330107,
c    LC: TA347.F5.S3313.
c
c    Gilbert Strang, George Fix,
c    An Analysis of the Finite Element Method,
c    Cambridge, 1973,
c    ISBN: 096140888X,
c    LC: TA335.S77.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c    Olgierd Zienkiewicz,
c    The Finite Element Method,
c    Sixth Edition,
c    Butterworth-Heinemann, 2005,
c    ISBN: 0750663200,
c    LC: TA640.2.Z54
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c
c    Input, integer ORDER, the order of the rule.
c
c    Output, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Output, double precision WEIGHT(ORDER), the weights of the rule.
c
      implicit none

      integer order

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision f
      double precision g
      double precision h
      integer i
      integer j
      integer k
      integer order2
      double precision p
      double precision q
      double precision r
      integer rule
      double precision s
      double precision t
      double precision u
      double precision v
      double precision w
      double precision w1
      double precision w2
      double precision w3
      double precision w4
      double precision w5
      double precision w6
      double precision w7
      double precision w8
      double precision w9
      double precision weight(order)
      double precision weight1(8)
      double precision weight2(8)
      double precision wx
      double precision x
      double precision xtab(order)
      double precision xtab1(8)
      double precision xtab2(8)
      double precision y
      double precision ytab(order)
      double precision z
c
c  1 point, precision 1.
c
      if ( rule .eq. 1 ) then

        xtab(1)   = 0.33333333333333333333D+00

        ytab(1)   = 0.33333333333333333333D+00

        weight(1) = 1.00000000000000000000D+00
c
c  3 points, precision 1, the "vertex rule".
c
      else if ( rule .eq. 2 ) then

        xtab(1) =   1.00000000000000000000D+00
        xtab(2) =   0.00000000000000000000D+00
        xtab(3) =   0.00000000000000000000D+00

        ytab(1) =   0.00000000000000000000D+00
        ytab(2) =   1.00000000000000000000D+00
        ytab(3) =   0.00000000000000000000D+00

        weight(1) = 0.33333333333333333333D+00
        weight(2) = 0.33333333333333333333D+00
        weight(3) = 0.33333333333333333333D+00
c
c  3 points, precision 2, Strang and Fix formula #1.
c
      else if ( rule .eq. 3 ) then

        xtab(1)   = 0.66666666666666666667D+00
        xtab(2)   = 0.16666666666666666667D+00
        xtab(3)   = 0.16666666666666666667D+00

        ytab(1)   = 0.16666666666666666667D+00
        ytab(2)   = 0.66666666666666666667D+00
        ytab(3)   = 0.16666666666666666667D+00

        weight(1) = 0.33333333333333333333D+00
        weight(2) = 0.33333333333333333333D+00
        weight(3) = 0.33333333333333333333D+00
c
c  3 points, precision 2, Strang and Fix formula #2.
c
      else if ( rule .eq. 4 ) then

        xtab(1)   = 0.50000000000000000000D+00
        xtab(2)   = 0.50000000000000000000D+00
        xtab(3)   = 0.00000000000000000000D+00

        ytab(1)   = 0.00000000000000000000D+00
        ytab(2)   = 0.50000000000000000000D+00
        ytab(3)   = 0.50000000000000000000D+00

        weight(1) = 0.33333333333333333333D+00
        weight(2) = 0.33333333333333333333D+00
        weight(3) = 0.33333333333333333333D+00
c
c  4 points, precision 3, Strang and Fix formula #3.
c
      else if ( rule .eq. 5 ) then

        a =   6.0D+00
        b =  10.0D+00
        c =  18.0D+00
        d =  25.0D+00
        e = -27.0D+00
        f =  30.0D+00
        g =  48.0D+00

        xtab(1) = b / f 
        xtab(2) = c / f
        xtab(3) = a / f
        xtab(4) = a / f

        ytab(1) = b / f
        ytab(2) = a / f
        ytab(3) = c / f
        ytab(4) = a / f

        weight(1) = e / g
        weight(2) = d / g
        weight(3) = d / g
        weight(4) = d / g
c
c  6 points, precision 3, Strang and Fix formula #4.
c
      else if ( rule .eq. 6 ) then

        a = 0.659027622374092D+00
        b = 0.231933368553031D+00
        c = 0.109039009072877D+00

        xtab(1) = a
        xtab(2) = a
        xtab(3) = b
        xtab(4) = b
        xtab(5) = c
        xtab(6) = c

        ytab(1) = b
        ytab(2) = c
        ytab(3) = a
        ytab(4) = c
        ytab(5) = a
        ytab(6) = b

        weight(1) = 0.16666666666666666667D+00
        weight(2) = 0.16666666666666666667D+00
        weight(3) = 0.16666666666666666667D+00
        weight(4) = 0.16666666666666666667D+00
        weight(5) = 0.16666666666666666667D+00
        weight(6) = 0.16666666666666666667D+00
c
c  6 points, precision 3, Stroud T2:3-1.
c
      else if ( rule .eq. 7 ) then

        a = 0.0D+00
        b = 0.5D+00
        c = 2.0D+00 /  3.0D+00
        d = 1.0D+00 /  6.0D+00
        v = 1.0D+00 / 30.0D+00
        w = 3.0D+00 / 10.0D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = b
        xtab(4) = c
        xtab(5) = d 
        xtab(6) = d

        ytab(1) = b
        ytab(2) = a
        ytab(3) = b
        ytab(4) = d
        ytab(5) = c
        ytab(6) = d

        weight(1) = v
        weight(2) = v
        weight(3) = v
        weight(4) = w
        weight(5) = w 
        weight(6) = w
c
c  6 points, precision 4, Strang and Fix, formula #5.
c
      else if ( rule .eq. 8 ) then

        a = 0.816847572980459D+00
        b = 0.091576213509771D+00
        c = 0.108103018168070D+00
        d = 0.445948490915965D+00
        v = 0.109951743655322D+00
        w = 0.223381589678011D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = b
        xtab(4) = c
        xtab(5) = d
        xtab(6) = d

        ytab(1) = b
        ytab(2) = a
        ytab(3) = b
        ytab(4) = d
        ytab(5) = c
        ytab(6) = d

        weight(1) = v
        weight(2) = v
        weight(3) = v
        weight(4) = w
        weight(5) = w
        weight(6) = w
c
c  7 points, precision 4, Strang and Fix formula #6.
c
      else if ( rule .eq. 9 ) then

        a = 1.0D+00 / 3.0D+00
        c = 0.736712498968435D+00
        d = 0.237932366472434D+00
        e = 0.025355134551932D+00
        v = 0.375000000000000D+00
        w = 0.104166666666667D+00

        xtab(1) = a
        xtab(2) = c
        xtab(3) = c
        xtab(4) = d
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e

        ytab(1) = a
        ytab(2) = d
        ytab(3) = e
        ytab(4) = c
        ytab(5) = e
        ytab(6) = c
        ytab(7) = d

        weight(1) = v
        weight(2) = w
        weight(3) = w
        weight(4) = w
        weight(5) = w
        weight(6) = w
        weight(7) = w
c
c  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
c
      else if ( rule .eq. 10 ) then

        a = 1.0D+00 / 3.0D+00
        b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
        c = ( 6.0D+00 -           sqrt ( 15.0D+00 ) ) / 21.0D+00
        d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
        e = ( 6.0D+00 +           sqrt ( 15.0D+00 ) ) / 21.0D+00
        u = 0.225D+00
        v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
        w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e

        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = e
        ytab(6) = d
        ytab(7) = e

        weight(1) = u
        weight(2) = v
        weight(3) = v
        weight(4) = v
        weight(5) = w
        weight(6) = w
        weight(7) = w
c
c  9 points, precision 6, Strang and Fix formula #8.
c
      else if ( rule .eq. 11 ) then

        a = 0.124949503233232D+00
        b = 0.437525248383384D+00
        c = 0.797112651860071D+00
        d = 0.165409927389841D+00
        e = 0.037477420750088D+00

        u = 0.205950504760887D+00
        v = 0.063691414286223D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = b
        xtab(4) = c
        xtab(5) = c
        xtab(6) = d
        xtab(7) = d
        xtab(8) = e
        xtab(9) = e

        ytab(1) = b
        ytab(2) = a
        ytab(3) = b
        ytab(4) = d
        ytab(5) = e
        ytab(6) = c
        ytab(7) = e
        ytab(8) = c
        ytab(9) = d

        weight(1) = u
        weight(2) = u
        weight(3) = u
        weight(4) = v
        weight(5) = v
        weight(6) = v
        weight(7) = v
        weight(8) = v
        weight(9) = v
c
c  12 points, precision 6, Strang and Fix, formula #9.
c
      else if ( rule .eq. 12 ) then

        a = 0.873821971016996D+00
        b = 0.063089014491502D+00
        c = 0.501426509658179D+00
        d = 0.249286745170910D+00
        e = 0.636502499121399D+00
        f = 0.310352451033785D+00
        g = 0.053145049844816D+00

        u = 0.050844906370207D+00
        v = 0.116786275726379D+00
        w = 0.082851075618374D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = b
        xtab(4) = c
        xtab(5) = d
        xtab(6) = d
        xtab(7) = e
        xtab(8) = e
        xtab(9) = f
        xtab(10) = f
        xtab(11) = g
        xtab(12) = g

        ytab(1) = b
        ytab(2) = a
        ytab(3) = b
        ytab(4) = d
        ytab(5) = c
        ytab(6) = d
        ytab(7) = f
        ytab(8) = g
        ytab(9) = e
        ytab(10) = g
        ytab(11) = e
        ytab(12) = f

        weight(1) = u
        weight(2) = u
        weight(3) = u
        weight(4) = v
        weight(5) = v
        weight(6) = v
        weight(7) = w
        weight(8) = w
        weight(9) = w
        weight(10) = w
        weight(11) = w
        weight(12) = w
c
c  13 points, precision 7, Strang and Fix, formula #10.
c
c  Note that there is a typographical error in Strang and Fix
c  which lists the value of the XSI(3) component of the
c  last generator point as 0.4869... when it should be 0.04869...
c
      else if ( rule .eq. 13 ) then

        h = 1.0D+00 / 3.0D+00
        a = 0.479308067841923D+00
        b = 0.260345966079038D+00
        c = 0.869739794195568D+00
        d = 0.065130102902216D+00
        e = 0.638444188569809D+00
        f = 0.312865496004875D+00
        g = 0.048690315425316D+00

        w = -0.149570044467670D+00
        t =  0.175615257433204D+00
        u =  0.053347235608839D+00
        v =  0.077113760890257D+00

        xtab(1) = h
        xtab(2) = a
        xtab(3) = b
        xtab(4) = b
        xtab(5) = c
        xtab(6) = d
        xtab(7) = d
        xtab(8) = e
        xtab(9) = e
        xtab(10) = f
        xtab(11) = f
        xtab(12) = g
        xtab(13) = g

        ytab(1) = h
        ytab(2) = b
        ytab(3) = a
        ytab(4) = b
        ytab(5) = d
        ytab(6) = c
        ytab(7) = d
        ytab(8) = f
        ytab(9) = g
        ytab(10) = e
        ytab(11) = g
        ytab(12) = e
        ytab(13) = f

        weight(1) = w
        weight(2) = t
        weight(3) = t
        weight(4) = t
        weight(5) = u
        weight(6) = u
        weight(7) = u
        weight(8) = v
        weight(9) = v
        weight(10) = v
        weight(11) = v
        weight(12) = v
        weight(13) = v
c
c  7 points, precision 3.
c
      else if ( rule .eq. 14 ) then

        a = 1.0D+00 / 3.0D+00
        b = 1.0D+00
        c = 0.5D+00
        z = 0.0D+00

        u = 27.0D+00 / 60.0D+00
        v =  3.0D+00 / 60.0D+00
        w =  8.0D+00 / 60.0D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = z
        xtab(4) = z
        xtab(5) = z
        xtab(6) = c
        xtab(7) = c

        ytab(1) = a
        ytab(2) = z
        ytab(3) = b
        ytab(4) = z
        ytab(5) = c
        ytab(6) = z
        ytab(7) = c

        weight(1) = u
        weight(2) = v
        weight(3) = v
        weight(4) = v
        weight(5) = w
        weight(6) = w
        weight(7) = w
c
c  16 points, precision 5, Stroud T2:7-1.
c
      else if ( rule .eq. 15 ) then
c
c  Legendre rule of order 4.
c
        order2 = 4

        xtab(1) = -0.861136311594052575223946488893D+00
        xtab(2) = -0.339981043584856264802665759103D+00
        xtab(3) =  0.339981043584856264802665759103D+00
        xtab(4) =  0.861136311594052575223946488893D+00

        weight(1) = 0.347854845137453857373063949222D+00
        weight(2) = 0.652145154862546142626936050778D+00
        weight(3) = 0.652145154862546142626936050778D+00
        weight(4) = 0.347854845137453857373063949222D+00

        do i = 1, order2
          xtab1(i) = 0.5D+00 * ( xtab1(i) + 1.0D+00 )
        end do

        weight2(1) = 0.1355069134D+00
        weight2(2) = 0.2034645680D+00
        weight2(3) = 0.1298475476D+00
        weight2(4) = 0.0311809709D+00

        xtab2(1) = 0.0571041961D+00
        xtab2(2) = 0.2768430136D+00
        xtab2(3) = 0.5835904324D+00
        xtab2(4) = 0.8602401357D+00

        k = 0
        do i = 1, order2
          do j = 1, order2
            k = k + 1
            xtab(k) = xtab2(j)
            ytab(k) = xtab1(i) * ( 1.0D+00 - xtab2(j) )
            weight(k) = weight1(i) * weight2(j)
          end do
        end do
c
c  64 points, precision 15.
c
      else if ( rule .eq. 16 ) then
c
c  Legendre rule of order 8.
c
        order2 = 8

        xtab1(1) = -0.960289856497536231683560868569D+00
        xtab1(2) = -0.796666477413626739591553936476D+00
        xtab1(3) = -0.525532409916328985817739049189D+00
        xtab1(4) = -0.183434642495649804939476142360D+00
        xtab1(5) =  0.183434642495649804939476142360D+00
        xtab1(6) =  0.525532409916328985817739049189D+00
        xtab1(7) =  0.796666477413626739591553936476D+00
        xtab1(8) =  0.960289856497536231683560868569D+00

        weight1(1) = 0.101228536290376259152531354310D+00
        weight1(2) = 0.222381034453374470544355994426D+00
        weight1(3) = 0.313706645877887287337962201987D+00
        weight1(4) = 0.362683783378361982965150449277D+00
        weight1(5) = 0.362683783378361982965150449277D+00
        weight1(6) = 0.313706645877887287337962201987D+00
        weight1(7) = 0.222381034453374470544355994426D+00
        weight1(8) = 0.101228536290376259152531354310D+00

        weight2(1) = 0.00329519144D+00
        weight2(2) = 0.01784290266D+00
        weight2(3) = 0.04543931950D+00
        weight2(4) = 0.07919959949D+00
        weight2(5) = 0.10604735944D+00
        weight2(6) = 0.11250579947D+00
        weight2(7) = 0.09111902364D+00
        weight2(8) = 0.04455080436D+00

        xtab2(1) = 0.04463395529D+00
        xtab2(2) = 0.14436625704D+00
        xtab2(3) = 0.28682475714D+00
        xtab2(4) = 0.45481331520D+00
        xtab2(5) = 0.62806783542D+00
        xtab2(6) = 0.78569152060D+00
        xtab2(7) = 0.90867639210D+00
        xtab2(8) = 0.98222008485D+00

        k = 0
        do j = 1, order2
          do i = 1, order2
            k = k + 1
            xtab(k) = 1.0D+00 - xtab2(j)
            ytab(k) = 0.5D+00 * ( 1.0D+00 + xtab1(i) ) * xtab2(j)
            weight(k) = weight1(i) * weight2(j)
          end do
        end do
c
c  19 points, precision 8, from CUBTRI.
c
      else if ( rule .eq. 17 ) then

        a = 1.0D+00 / 3.0D+00
        b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
        c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
        d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
        e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
        f = ( 40.0D+00 - 10.0D+00 * sqrt ( 15.0D+00 ) 
     &    + 10.0D+00 * sqrt ( 7.0D+00 ) + 2.0D+00 * sqrt ( 105.0D+00 ) )
     &    / 90.0D+00
        g = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) 
     &    -  5.0D+00 * sqrt ( 7.0D+00 ) - sqrt ( 105.0D+00 ) ) 
     &    / 90.0D+00
        p = ( 40.0D+00 + 10.0D+00 * sqrt ( 15.0D+00 ) 
     &    + 10.0D+00 * sqrt ( 7.0D+00 ) - 2.0D+00 * sqrt ( 105.0D+00 ) )
     &   / 90.0D+00
        q = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) 
     &    -  5.0D+00 * sqrt ( 7.0D+00 ) + sqrt ( 105.0D+00 ) ) 
     &    / 90.0D+00
        r = ( 40.0D+00 + 10.0D+00 * sqrt ( 7.0D+00 ) ) / 90.0D+00
        s = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) 
     &    - 5.0D+00 * sqrt ( 7.0D+00 ) 
     &    - sqrt ( 105.0D+00 ) ) / 90.0D+00
        t = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) 
     &    - 5.0D+00 * sqrt ( 7.0D+00 ) 
     &    + sqrt ( 105.0D+00 ) ) / 90.0D+00

        w1 = ( 7137.0D+00 - 1800.0D+00 * sqrt ( 7.0D+00 ) ) 
     &    / 62720.0D+00
        w2 = -9301697.0D+00 / 4695040.0D+00 
     &    - 13517313.0D+00 * sqrt ( 15.0D+00 ) 
     &    / 23475200.0D+00 + 764885.0D+00 * sqrt ( 7.0D+00 ) 
     &    / 939008.0D+00 
     &    + 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
        w2 = w2 / 3.0D+00
        w3 = -9301697.0D+00 / 4695040.0D+00 + 13517313.0D+00 
     &    * sqrt ( 15.0D+00 ) 
     &    / 23475200.0D+00 
     &    + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 
     &    - 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
        w3 = w3 / 3.0D+00
        w4 = ( 102791225.0D+00 - 23876225.0D+00 * sqrt ( 15.0D+00 ) 
     &    - 34500875.0D+00 * sqrt ( 7.0D+00 ) 
     &    + 9914825.0D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
        w4 = w4 / 3.0D+00
        w5 = ( 102791225.0D+00 + 23876225.0D+00 * sqrt ( 15.0D+00 ) 
     &    - 34500875.0D+00 * sqrt ( 7.0D+00 ) 
     &    - 9914825D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
        w5 = w5 / 3.0D+00
        w6 = ( 11075.0D+00 - 3500.0D+00 * sqrt ( 7.0D+00 ) ) 
     &    / 8064.0D+00
        w6 = w6 / 6.0D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e
        xtab(8) = f
        xtab(9) = g
        xtab(10) = g
        xtab(11) = p
        xtab(12) = q
        xtab(13) = q
        xtab(14) = r
        xtab(15) = r
        xtab(16) = s
        xtab(17) = s
        xtab(18) = t
        xtab(19) = t

        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = e
        ytab(6) = d
        ytab(7) = e
        ytab(8) = g
        ytab(9) = f
        ytab(10) = g
        ytab(11) = q
        ytab(12) = p
        ytab(13) = q          
        ytab(14) = s
        ytab(15) = t
        ytab(16) = r
        ytab(17) = t
        ytab(18) = r
        ytab(19) = s

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w2
        weight(5) = w3
        weight(6) = w3
        weight(7) = w3
        weight(8) = w4
        weight(9) = w4
        weight(10) = w4
        weight(11) = w5
        weight(12) = w5
        weight(13) = w5
        weight(14) = w6
        weight(15) = w6
        weight(16) = w6
        weight(17) = w6
        weight(18) = w6
        weight(19) = w6
c
c  19 points, precision 9.
c  Lyness and Jesperson.
c
      else if ( rule .eq. 18 ) then

        a = 1.0D+00 / 3.0D+00
        b =  0.02063496160252593D+00
        c =  0.4896825191987370D+00
        d =  0.1258208170141290D+00
        e =  0.4370895914929355D+00
        f =  0.6235929287619356D+00
        g =  0.1882035356190322D+00
        r =  0.9105409732110941D+00
        s =  0.04472951339445297D+00
        t =  0.7411985987844980D+00
        u =  0.03683841205473626D+00
        v =  0.22196298916076574D+00

        w1 = 0.09713579628279610D+00
        w2 = 0.03133470022713983D+00
        w3 = 0.07782754100477543D+00
        w4 = 0.07964773892720910D+00
        w5 = 0.02557767565869810D+00
        w6 = 0.04328353937728940D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e
        xtab(8) = f
        xtab(9) = g
        xtab(10) = g
        xtab(11) = r
        xtab(12) = s
        xtab(13) = s
        xtab(14) = t
        xtab(15) = t
        xtab(16) = u
        xtab(17) = u
        xtab(18) = v
        xtab(19) = v

        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = e
        ytab(6) = d
        ytab(7) = e
        ytab(8) = g
        ytab(9) = f
        ytab(10) = g
        ytab(11) = s
        ytab(12) = r
        ytab(13) = s
        ytab(14) = u
        ytab(15) = v
        ytab(16) = t
        ytab(17) = v
        ytab(18) = t
        ytab(19) = u

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w2
        weight(5) = w3
        weight(6) = w3
        weight(7) = w3
        weight(8) = w4
        weight(9) = w4
        weight(10) = w4
        weight(11) = w5
        weight(12) = w5
        weight(13) = w5
        weight(14) = w6
        weight(15) = w6
        weight(16) = w6
        weight(17) = w6
        weight(18) = w6
        weight(19) = w6
c
c  28 points, precision 11.
c  Lyness and Jesperson.
c
      else if ( rule .eq. 19 ) then

        a = 1.0D+00 / 3.0D+00
        b = 0.9480217181434233D+00
        c = 0.02598914092828833D+00
        d = 0.8114249947041546D+00
        e = 0.09428750264792270D+00
        f = 0.01072644996557060D+00
        g = 0.4946367750172147D+00
        p = 0.5853132347709715D+00
        q = 0.2073433826145142D+00
        r = 0.1221843885990187D+00
        s = 0.4389078057004907D+00
        t = 0.6779376548825902D+00
        u = 0.04484167758913055D+00
        v = 0.27722066752827925D+00
        w = 0.8588702812826364D+00
        x = 0.0D+00
        y = 0.1411297187173636D+00

        w1 = 0.08797730116222190D+00
        w2 = 0.008744311553736190D+00
        w3 = 0.03808157199393533D+00
        w4 = 0.01885544805613125D+00
        w5 = 0.07215969754474100D+00
        w6 = 0.06932913870553720D+00
        w7 = 0.04105631542928860D+00
        w8 = 0.007362383783300573D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e
        xtab(8) = f
        xtab(9) = g
        xtab(10) = g
        xtab(11) = p
        xtab(12) = q
        xtab(13) = q
        xtab(14) = r
        xtab(15) = s
        xtab(16) = s
        xtab(17) = t
        xtab(18) = t
        xtab(19) = u
        xtab(20) = u
        xtab(21) = v
        xtab(22) = v
        xtab(23) = w
        xtab(24) = w
        xtab(25) = x
        xtab(26) = x
        xtab(27) = y
        xtab(28) = y

        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = e
        ytab(6) = d
        ytab(7) = e
        ytab(8) = g
        ytab(9) = f
        ytab(10) = g
        ytab(11) = q
        ytab(12) = p
        ytab(13) = q
        ytab(14) = s
        ytab(15) = r
        ytab(16) = s
        ytab(17) = u
        ytab(18) = v
        ytab(19) = t
        ytab(20) = v
        ytab(21) = t
        ytab(22) = u
        ytab(23) = x
        ytab(24) = y
        ytab(25) = w
        ytab(26) = y
        ytab(27) = w
        ytab(28) = x

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w2
        weight(5) = w3
        weight(6) = w3
        weight(7) = w3
        weight(8) = w4
        weight(9) = w4
        weight(10) = w4
        weight(11) = w5
        weight(12) = w5
        weight(13) = w5
        weight(14) = w6
        weight(15) = w6
        weight(16) = w6
        weight(17) = w7
        weight(18) = w7
        weight(19) = w7
        weight(20) = w7
        weight(21) = w7
        weight(22) = w7
        weight(23) = w8
        weight(24) = w8
        weight(25) = w8
        weight(26) = w8
        weight(27) = w8
        weight(28) = w8
c
c  37 points, precision 13.
c
      else if ( rule .eq. 20 ) then

        a = 1.0D+00 / 3.0D+00
        b = 0.950275662924105565450352089520D+00
        c = 0.024862168537947217274823955239D+00
        d = 0.171614914923835347556304795551D+00
        e = 0.414192542538082326221847602214D+00
        f = 0.539412243677190440263092985511D+00
        g = 0.230293878161404779868453507244D+00

        w1 = 0.051739766065744133555179145422D+00
        w2 = 0.008007799555564801597804123460D+00
        w3 = 0.046868898981821644823226732071D+00
        w4 = 0.046590940183976487960361770070D+00
        w5 = 0.031016943313796381407646220131D+00
        w6 = 0.010791612736631273623178240136D+00
        w7 = 0.032195534242431618819414482205D+00
        w8 = 0.015445834210701583817692900053D+00
        w9 = 0.017822989923178661888748319485D+00
        wx = 0.037038683681384627918546472190D+00

        xtab(1) = a
        xtab(2) = b
        xtab(3) = c
        xtab(4) = c
        xtab(5) = d
        xtab(6) = e
        xtab(7) = e
        xtab(8) = f
        xtab(9) = g
        xtab(10) = g

        ytab(1) = a
        ytab(2) = c
        ytab(3) = b
        ytab(4) = c
        ytab(5) = e
        ytab(6) = d
        ytab(7) = e
        ytab(8) = g
        ytab(9) = f
        ytab(10) = g

        weight(1) = w1
        weight(2) = w2
        weight(3) = w2
        weight(4) = w2
        weight(5) = w3
        weight(6) = w3
        weight(7) = w3
        weight(8) = w4
        weight(9) = w4
        weight(10) = w4
        weight(11) = w5
        weight(12) = w5
        weight(13) = w5
        weight(14) = w6
        weight(15) = w6
        weight(16) = w6
        weight(17) = w7
        weight(18) = w7
        weight(19) = w7
        weight(20) = w8
        weight(21) = w8
        weight(22) = w8
        weight(23) = w8
        weight(24) = w8
        weight(25) = w8
        weight(26) = w9
        weight(27) = w9
        weight(28) = w9
        weight(29) = w9
        weight(30) = w9
        weight(31) = w9
        weight(32) = wx
        weight(33) = wx
        weight(34) = wx
        weight(35) = wx
        weight(36) = wx
        weight(37) = wx

        a = 0.772160036676532561750285570113D+00
        b = 0.113919981661733719124857214943D+00

        xtab(11) = a
        ytab(11) = b

        xtab(12) = b
        ytab(12) = a

        xtab(13) = b
        ytab(13) = b

        a = 0.009085399949835353883572964740D+00
        b = 0.495457300025082323058213517632D+00

        xtab(14) = a
        ytab(14) = b

        xtab(15) = b
        ytab(15) = a

        xtab(16) = b
        ytab(16) = b

        a = 0.062277290305886993497083640527D+00
        b = 0.468861354847056503251458179727D+00

        xtab(17) = a
        ytab(17) = b

        xtab(18) = b
        ytab(18) = a

        xtab(19) = b
        ytab(19) = b

        a = 0.022076289653624405142446876931D+00
        b = 0.851306504174348550389457672223D+00
        c = 0.126617206172027096933163647918263D+00

        xtab(20) = a
        ytab(20) = b

        xtab(21) = a
        ytab(21) = c

        xtab(22) = b
        ytab(22) = a

        xtab(23) = b
        ytab(23) = c

        xtab(24) = c
        ytab(24) = a

        xtab(25) = c
        ytab(25) = b

        a = 0.018620522802520968955913511549D+00
        b = 0.689441970728591295496647976487D+00
        c = 0.291937506468887771754472382212953D+00

        xtab(26) = a
        ytab(26) = b

        xtab(27) = a
        ytab(27) = c

        xtab(28) = b
        ytab(28) = a

        xtab(29) = b
        ytab(29) = c

        xtab(30) = c
        ytab(30) = a

        xtab(31) = c
        ytab(31) = b

        a = 0.096506481292159228736516560903D+00
        b = 0.635867859433872768286976979827D+00
        c = 0.267625659273967961282458816185681D+00

        xtab(32) = a
        ytab(32) = b

        xtab(33) = a
        ytab(33) = c

        xtab(34) = b
        ytab(34) = a

        xtab(35) = b
        ytab(35) = c

        xtab(36) = c
        ytab(36) = a

        xtab(37) = c
        ytab(37) = b

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_UNIT_SET - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
        stop

      end if

      return
      end
      function triangle_unit_size ( rule )

c*********************************************************************72
c
cc TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706,
c    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c    DP Laurie,
c    Algorithm 584,
c    CUBTRI, Automatic Cubature Over a Triangle,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 2, 1982, pages 210-218.
c
c    James Lyness, Dennis Jespersen,
c    Moderate Degree Symmetric Quadrature Rules for the Triangle,
c    Journal of the Institute of Mathematics and its Applications,
c    Volume 15, Number 1, February 1975, pages 19-32.
c
c    Hans Rudolf Schwarz,
c    Methode der Finiten Elemente,
c    Teubner Studienbuecher, 1980,
c    ISBN: 3-519-02349-0.
c
c    Gilbert Strang, George Fix,
c    An Analysis of the Finite Element Method,
c    Prentice Hall, 1973, page 184,
c    ISBN: 096140888X,
c    LC: TA335.S77.
c
c    Arthur Stroud,
c    Approximate Calculation of Multiple Integrals,
c    Prentice Hall, 1971,
c    ISBN: 0130438936,
c    LC: QA311.S85.
c
c    Olgierd Zienkiewicz,
c    The Finite Element Method,
c    Sixth Edition,
c    Butterworth-Heinemann, 2005,
c    ISBN: 0750663200,
c    TA640.2.Z54
c
c  Parameters:
c
c    Input, integer RULE, the index of the rule.
c     1, ORDER =  1, precision 1, Zienkiewicz #1.
c     2, ORDER =  2, precision 1, (the "vertex rule").
c     3, ORDER =  3, precision 2, Strang and Fix formula #1.
c     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
c     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
c     6, ORDER =  6, precision 3, Strang and Fix formula #4.
c     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
c     8, ORDER =  6, precision 4, Strang and Fix formula #5.
c     9, ORDER =  7, precision 4, Strang and Fix formula #6.
c    10, ORDER =  7, precision 5, Strang and Fix formula #7,
c        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
c    11, ORDER =  9, precision 6, Strang and Fix formula #8.
c    12, ORDER = 12, precision 6, Strang and Fix formula #9.
c    13, ORDER = 13, precision 7, Strang and Fix formula #10.
c    14, ORDER =  7, precision ?.
c    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
c    16, ORDER = 64, precision 15, triangular product Gauss rule.
c    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
c    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
c    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
c    20, ORDER = 37, precision 13, from ACM TOMS #706.
c
c    Output, integer TRIANGLE_UNIT_SIZE, the order of the rule.
c
      implicit none

      integer rule
      integer triangle_unit_size

      if ( rule .eq. 1 ) then
        triangle_unit_size = 1
      else if ( rule .eq. 2 ) then
        triangle_unit_size = 3
      else if ( rule .eq. 3 ) then
        triangle_unit_size = 3
      else if ( rule .eq. 4 ) then
        triangle_unit_size = 3
      else if ( rule .eq. 5 ) then
        triangle_unit_size = 4
      else if ( rule .eq. 6 ) then
        triangle_unit_size = 6
      else if ( rule .eq. 7 ) then
        triangle_unit_size = 6
      else if ( rule .eq. 8 ) then
        triangle_unit_size = 6
      else if ( rule .eq. 9 ) then
        triangle_unit_size = 7
      else if ( rule .eq. 10 ) then
        triangle_unit_size = 7
      else if ( rule .eq. 11 ) then
        triangle_unit_size = 9
      else if ( rule .eq. 12 ) then
        triangle_unit_size = 12
      else if ( rule .eq. 13 ) then
        triangle_unit_size = 13
      else if ( rule .eq. 14 ) then
        triangle_unit_size = 7
      else if ( rule .eq. 15 ) then
        triangle_unit_size = 16
      else if ( rule .eq. 16 ) then
        triangle_unit_size = 64
      else if ( rule .eq. 17 ) then
        triangle_unit_size = 19
      else if ( rule .eq. 18 ) then
        triangle_unit_size = 19
      else if ( rule .eq. 19 ) then
        triangle_unit_size = 28
      else if ( rule .eq. 20 ) then
        triangle_unit_size = 37
      else
        triangle_unit_size = -1
      end if

      return
      end
      subroutine triangle_unit_sum ( func, order, xtab, ytab, weight, 
     &  result )

c*********************************************************************72
c
cc TRIANGLE_UNIT_SUM carries out a quadrature rule in the unit triangle.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external FUNC, the name of the user supplied
c    function of two variables which is to be integrated,
c    of the form:
c      function func ( x, y )
c      double precision func
c      double precision x
c      double precision y
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision XTAB(ORDER), YTAB(ORDER), the abscissas.
c
c    Input, double precision WEIGHT(ORDER), the weights of the rule.
c
c    Output, double precision RESULT, the approximate integral of the function.
c
      implicit none

      integer order

      double precision func
      external func
      integer i
      double precision quad
      double precision result
      double precision triangle_unit_volume
      double precision volume
      double precision weight(order)
      double precision xtab(order)
      double precision ytab(order)

      quad = 0.0D+00

      do i = 1, order
        quad = quad + weight(i) * func ( xtab(i), ytab(i) )
      end do

      volume = triangle_unit_volume ( )
      result = quad * volume

      return
      end
      function triangle_unit_volume ( )

c*********************************************************************72
c
cc TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Discussion:
c
c    The "volume" of a triangle is usually called its area.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision TRIANGLE_UNIT_VOLUME, the volume of the unit
c    triangle.
c
      implicit none

      double precision triangle_unit_volume

      triangle_unit_volume = 1.0D+00 / 2.0D+00

      return
      end
      function triangle_volume ( x, y )

c*********************************************************************72
c
cc TRIANGLE_VOLUME returns the "volume" of a triangle in 2D.
c
c  Integration region:
c
c      0 <= X,
c    and
c      0 <= Y, 
c    and
c      X + Y <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(3), Y(3), the vertices of the triangle.
c
c    Output, double precision TRIANGLE_VOLUME, the volume of the triangle.
c
      implicit none

      double precision triangle_volume
      double precision x(3)
      double precision y(3)

      triangle_volume = 0.5D+00 * abs ( 
     &  x(1) * ( y(2) - y(3) ) + 
     &  x(2) * ( y(3) - y(1) ) + 
     &  x(3) * ( y(1) - y(2) ) )

      return
      end
      subroutine tvec_even ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The computation realizes that 0 = 2 * PI.
c
c  Example:
c
c    NT = 4
c
c    T = ( 0, PI/2, PI, 3*PI/2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      do i = 1, nt
        t(i) = dble ( 2 * ( i - 1 ) ) * pi / dble ( nt )
      end do

      return
      end
      subroutine tvec_even2 ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The computation realizes that 0 = 2 * PI.  The values are equally
c    spaced in the circle, do not include 0, and are symmetric about 0.
c
c  Example:
c
c    NT = 4
c
c    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      do i = 1, nt
        t(i) = dble ( 2 * i - 1 ) * pi / dble ( nt )
      end do

      return
      end
      subroutine tvec_even3 ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN3 computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The angles begin with 0 and end with 2*PI.
c
c  Example:
c
c    NT = 4
c
c    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      if ( nt .eq. 1 ) then
        t(1) = pi
      else
        do i = 1, nt
          t(i) = dble ( 2 * ( i - 1 ) ) * pi / dble ( nt - 1 )
        end do
      end if

      return
      end
      subroutine tvec_even_bracket ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
c
c    The angles returned are the breakpoints of these subintervals,
c    including THETA1 and THETA2.
c
c  Example:
c
c    NT = 4
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 30, 50, 70, 90 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      if ( nt .eq. 1 ) then

        t(1) = ( theta1 + theta2 ) / 2.0D+00

      else

        do i = 1, nt
          t(i) = ( dble ( nt - i     ) * theta1   
     &           + dble (      i - 1 ) * theta2 ) 
     &           / dble ( nt     - 1 )
        end do

      end if

      return
      end
      subroutine tvec_even_bracket2 ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
c
c    The angles returned are the internal NT breakpoints of the subintervals.
c
c  Example:
c
c    NT = 5
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 40, 50, 60, 70, 80 )
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
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      do i = 1, nt
        t(i) = ( dble ( nt + 1 - i ) * theta1   
     &         + dble (          i ) * theta2 ) 
     &         / dble ( nt + 1     )
      end do

      return
      end
      subroutine tvec_even_bracket3 ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET3 computes evenly spaced angles between THETA1 and THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT subintervals.
c
c    The angles returned are the midpoints of each subinterval.
c
c  Example:
c
c    NT = 3
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 40, 60, 80 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      do i = 1, nt
        t(i) = ( dble ( 2 * nt - 2 * i + 1 ) * theta1   
     &         + dble (          2 * i - 1 ) * theta2 ) 
     &         / dble ( 2 * nt             )
      end do

      return
      end
      subroutine vec_lex_next ( dim_num, base, a, more )

c*********************************************************************72
c
cc VEC_LEX_NEXT generates vectors in lex order.
c
c  Discussion:
c
c    The vectors are produced in lexical order, starting with
c    (0,0,...,0),
c    (0,0,...,1),
c    ...
c    (BASE-1,BASE-1,...,BASE-1).
c
c  Example:
c
c    DIM_NUM = 2,
c    BASE = 3
c
c    0   0
c    0   1
c    0   2
c    1   0
c    1   1
c    1   2
c    2   0
c    2   1
c    2   2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the size of the vectors to be used.
c
c    Input, integer BASE, the base to be used.  BASE = 2 will
c    give vectors of 0's and 1's, for instance.
c
c    Input/output, integer A(DIM_NUM), the next vector.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  On return, MORE is TRUE if another vector has
c    been computed.  If MORE is returned FALSE, ignore the output
c    vector and stop calling the routine.
c
      implicit none

      integer dim_num

      integer a(dim_num)
      integer base
      integer i
      logical more

      if ( .not. more ) then

        do i = 1, dim_num
          a(i) = 0
        end do
        more = .true.

      else

        do i = dim_num, 1, -1

          a(i) = a(i) + 1

          if ( a(i) < base ) then
            return
          end if

          a(i) = 0

        end do

        more = .false.

      end if

      return
      end
