      subroutine angle_box_2d ( dist, p1, p2, p3, p4, p5 )

c*********************************************************************72
c
cc ANGLE_BOX_2D "boxes" an angle defined by three points in 2D.
c
c  Discussion:
c
c    The routine is given points P1, P2 and P3, determining the two lines:
c      P1 to P2
c    and
c      P2 to P3
c    and a nonnegative distance
c      DIST.
c
c    The routine returns a pair of "corner" points
c      P4 and P5
c    both of which are a distance DIST from both lines, and in fact,
c    both of which are a distance DIST from P2.
c
c                         /  P3
c                        /   /   /
c     - - - - - - - - -P4 - / -P6 - - -
c                      /   /   /
c    P1---------------/--P2-----------------
c                    /   /   /
c     - - - - - - -P7 - / -P5 - - - - -
c                  /   /   /
c
c    In the illustration, P1, P2 and P3 are the points defining the lines.
c
c    P4 and P5 represent the desired "corner points", which
c    are on the positive or negative sides of both lines.
c
c    P6 and P7 represent the undesired points, which
c    are on the positive side of one line and the negative of the other.
c
c    Special cases:
c
c    if P1 = P2, this is the same as extending the line from
c    P3 through P2 without a bend.
c
c    if P3 = P2, this is the same as extending the line from
c    P1 through P2 without a bend.
c
c    if P1 = P2 = P3 this is an error.
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
c    Input, double precision DIST, the nonnegative distance from P1
c    to the computed points P4 and P5.
c
c    Input, double precision P1(2), P2(2), P3(2).
c    P1 and P2 are distinct points that define a line.
c    P2 and P3 are distinct points that define a line.
c
c    Output, double precision P4(2), P5(2), points which lie DIST units from
c    the line between P1 and P2, and from the line between P2 and P3.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      integer i
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision p5(dim_num)
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision stheta
      double precision temp1
      double precision temp2
      double precision temp3
      double precision u(dim_num)
      double precision u1(dim_num)
      double precision u2(dim_num)
c
c  If DIST = 0, assume the user knows best.
c
      if ( dist .eq. 0.0D+00 ) then
        do i = 1, dim_num
          p4(i) = p2(i)
          p5(i) = p2(i)
        end do
        return
      end if
c
c  Fail if all three points are equal.
c
      if ( p1(1) .eq. p2(1) .and.
     &     p1(2) .eq. p2(2) .and.
     &     p2(1) .eq. p3(1) .and.
     &     p2(2) .eq. p3(2) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ANGLE_BOX_2D - Fatal error!'
        write ( *, '(a)' ) '  Input points P1 = P2 = P3.'
        write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1), p1(2)
        stop
      end if
c
c  If P1 = P2, extend the line through the doubled point.
c
      if ( p1(1) .eq. p2(1) .and.
     &     p1(2) .eq. p2(2) ) then
        u2(1) = p3(2) - p2(2)
        u2(2) = p2(1) - p3(1)
        temp1 = r8vec_norm ( dim_num, u2 )
        do i = 1, dim_num
          u2(i) = u2(i) / temp1
        end do
        do i = 1, dim_num
          p4(i) = p2(i) + dist * u2(i)
          p5(i) = p2(i) - dist * u2(i)
        end do
        return
      end if
c
c  If P2 = P3, extend the line through the doubled point.
c
      if ( p2(1) .eq. p3(1) .and.
     &     p2(2) .eq. p3(2) ) then
        u1(1) = p1(2) - p2(2)
        u1(2) = p2(1) - p1(1)
        temp1 = r8vec_norm ( dim_num, u1 )
        do i = 1, dim_num
          u1(i) = u1(i) / temp1
        end do
        do i = 1, dim_num
          p4(i) = p2(i) + dist * u1(i)
          p5(i) = p2(i) - dist * u1(i)
        end do
        return
      end if
c
c  Compute the unit normal vectors to each line.
c  We choose the sign so that the unit normal to line 1 has
c  a positive dot product with line 2.
c
      u1(1) = p1(2) - p2(2)
      u1(2) = p2(1) - p1(1)
      temp1 = r8vec_norm ( dim_num, u1 )
      do i = 1, dim_num
        u1(i) = u1(i) / temp1
      end do

      temp1 = 0.0D+00
      do i = 1, dim_num
        temp1 = temp1 + u1(i) * ( p3(i) - p2(i) )
      end do

      if ( temp1 .lt. 0.0D+00 ) then
        do i = 1, dim_num
          u1(i) = - u1(i)
        end do
      end if

      u2(1) = p3(2) - p2(2)
      u2(2) = p2(1) - p3(1)
      temp1 = r8vec_norm ( dim_num, u2 )
      do i = 1, dim_num
        u2(i) = u2(i) / temp1
      end do

      temp1 = 0.0D+00
      do i = 1, dim_num
        temp1 = temp1 + u2(i) * ( p1(i) - p2(i) )
      end do

      if ( temp1 .lt. 0.0D+00 ) then
        do i = 1, dim_num
          u2(i) = - u2(i)
        end do
      end if
c
c  Try to catch the case where we can not determine the
c  sign of U1, because both U1 and -U1 are perpendicular
c  to (P3-P2)...and similarly for U2 and (P1-P2).
c
      temp1 = 0.0D+00
      do i = 1, dim_num
        temp1 = temp1 + u1(i) * ( p3(i) - p2(i) )
      end do

      temp2 = 0.0D+00
      do i = 1, dim_num
        temp2 = temp2 + u2(i) * ( p1(i) - p2(i) )
      end do

      if ( temp1 .eq. 0.0D+00 .or. temp2 .eq. 0.0D+00 ) then

        if ( r8vec_dot_product ( dim_num, u1, u2 ) .lt. 0.0D+00 ) then
          do i = 1, dim_num
            u1(i) = - u1(i)
          end do
        end if

      end if
c
c  Try to catch a line turning back on itself, evidenced by
c    Cos(theta) = (P3-P2) dot (P2-P1) / ( norm(P3-P2) * norm(P2-P1) )
c  being -1, or very close to -1.
c
      temp1 = 0.0D+00
      temp2 = 0.0D+00
      temp3 = 0.0D+00
      do i = 1, dim_num
        temp1 = temp1 + ( p3(i) - p2(i) ) * ( p2(i) - p1(i) )
        temp2 = temp2 + ( p3(i) - p2(i) )**2
        temp3 = temp3 + ( p2(i) - p1(i) )**2
      end do
      temp2 = sqrt ( temp2 )
      temp3 = sqrt ( temp3 )

      temp1 = temp1 / temp2 / temp3

      if ( temp1 .lt. - 0.99D+00 ) then
        temp1 = 0.0D+00
        do i = 1, dim_num
          temp1 = temp1 + ( p2(i) - p1(i) )**2
        end do
        temp1 = sqrt ( temp1 )
        do i = 1, dim_num
          p4(i) = p2(i) + dist * ( p2(i) - p1(i) )
     &      / temp1 + dist * u1(i)
          p5(i) = p2(i) + dist * ( p2(i) - p1(i) )
     &      / temp1 - dist * u1(i)
         end do
        return
      end if
c
c  Compute the "average" unit normal vector.
c
c  The average of the unit normals could be zero, but only when
c  the second line has the same direction and opposite sense
c  of the first, and we have already checked for that case.
c
c  Well, check again.  This problem "bit" me in the case where
c  P1 = P2, which I now treat specially just to guarantee I
c  avoid this problemc
c
      if ( r8vec_dot_product ( dim_num, u1, u2 ) .lt. 0.0D+00 ) then
        do i = 1, dim_num
          u2(i) = - u2(i)
        end do
      end if

      do i = 1, dim_num
        u(i) = 0.5D+00 * ( u1(i) + u2(i) )
      end do
      temp1 = r8vec_norm ( dim_num, u )
      do i = 1, dim_num
        u(i) = u(i) / temp1
      end do
c
c  You must go DIST/STHETA units along this unit normal to
c  result in a distance DIST from line1 (and line2).
c
      stheta = r8vec_dot_product ( dim_num, u, u1 )

      do i = 1, dim_num
        p4(i) = p2(i) + dist * u(i) / stheta
        p5(i) = p2(i) - dist * u(i) / stheta
      end do

      return
      end
      subroutine angle_contains_point_2d ( p1, p2, p3, p, inside )

c*********************************************************************72
c
cc ANGLE_CONTAINS_POINT_2D determines if an angle contains a point, in 2D.
c
c  Discussion:
c
c    The angle is defined by the sequence of points P1, P2 and P3.
c
c    The point is "contained" by the angle if the ray P - P2
c    is between (in a counter clockwise sense) the rays P1 - P2
c    and P3 - P2.
c
c        P1
c        /
c       /   P
c      /  .
c     / .
c    P2--------->P3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), P3(2), the coordinates of
c    three points that define the angle.  The order of these points mattersc
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is inside the angle.
c
      implicit none

      double precision angle_rad_2d
      logical inside
      double precision p(2)
      double precision p1(2)
      double precision p2(2)
      double precision p3(2)

      if ( angle_rad_2d ( p1, p2, p ) .le.
     &     angle_rad_2d ( p1, p2, p3 ) ) then
        inside = .true.
      else
        inside = .false.
      end if

      return
      end
      function angle_deg_2d ( p1, p2, p3 )

c*********************************************************************72
c
cc ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
c
c  Discussion:
c
c    Except for the zero angle case, it should be true that
c
c      ANGLE_DEG_2D ( P1, P2, P3 ) + ANGLE_DEG_2D ( P3, P2, P1 ) = 360.0
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
c    25 June 2009
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
c    Output, double precision ANGLE_DEG_2D, the angle swept out by the
c    rays, measured in degrees.  0 .le. ANGLE_DEG_2D .lt. 360.  If either ray
c    has zero length, then ANGLE_DEG_2D is set to 0.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle_deg_2d
      double precision angle_rad_2d
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision radians_to_degrees
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)

      p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) )
     &     + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

      p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) )
     &     - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

      if ( p(1) .eq. 0.0D+00 .and. p(2) .eq. 0.0D+00 ) then
        angle_deg_2d = 0.0D+00
        return
      end if

      angle_rad_2d = atan2 ( p(2), p(1) )

      if ( angle_rad_2d .lt. 0.0D+00 ) then
        angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
      end if

      angle_deg_2d = radians_to_degrees ( angle_rad_2d )

      return
      end
      subroutine angle_half_2d ( p1, p2, p3, p4 )

c*********************************************************************72
c
cc ANGLE_HALF_2D finds half an angle in 2D.
c
c  Discussion:
c
c    The original angle is defined by the sequence of points P1, P2 and P3.
c
c    The point P4 is calculated so that:
c
c      (P1,P2,P4) = (P1,P2,P3) / 2
c
c        P1
c        /
c       /   P4
c      /  .
c     / .
c    P2--------->P3
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
c    Input, double precision P1(2), P2(2), P3(2), points defining the angle.
c
c    Input, double precision P4(2), a point defining the half angle.
c    The vector P4 - P2 will have unit norm.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      double precision norm
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)

      do i = 1, dim_num

        p4(i) = 0.5E+00 *
     &    ( ( p1(i) - p2(i) )
     &    / sqrt ( ( p1(1) - p2(1) )**2 + ( p1(2) - p2(2) )**2 )
     &    + ( p3(i) - p2(i) )
     &    / sqrt ( ( p3(1) - p2(1) )**2 + ( p3(2) - p2(2) )**2 ) )

      end do

      norm = sqrt ( p4(1)**2 + p4(2)**2 )

      do i = 1, dim_num
        p4(i) = p2(i) + p4(i) / norm
      end do

      return
      end
      function angle_rad_2d ( p1, p2, p3 )

c*********************************************************************72
c
cc ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
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
c    25 June 2009
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
c    in radians.  0 .le. ANGLE_RAD_2D .lt. 2 * PI.  If either ray has zero
c    length, then ANGLE_RAD_2D is set to 0.
c
      implicit none

      integer, parameter :: dim_num = 2

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

      if ( p(1) .eq. 0.0D+00 .and. p(2) .eq. 0.0D+00  ) then
        angle_rad_2d = 0.0D+00
        return
      end if

      angle_rad_2d = atan2 ( p(2), p(1) )

      if ( angle_rad_2d .lt. 0.0D+00 ) then
        angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
      end if

      return
      end
      function angle_rad_3d ( p1, p2, p3 )

c*********************************************************************72
c
cc ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
c
c  Discussion:
c
c    The routine always computes the SMALLER of the two angles between
c    two rays.  Thus, if the rays make an (exterior) angle of
c    1.5 pi radians, the (interior) angle of 0.5 pi radians will be reported.
c
c    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), points defining an angle.
c    The rays are P1 - P2 and P3 - P2.
c
c    Output, double precision ANGLE_RAD_3D, the angle between the two rays,
c    in radians.  This value will always be between 0 and PI.  If either ray has
c    zero length, then the angle is returned as zero.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision angle_rad_3d
      integer dim
      double precision dot
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision r8_acos
      double precision v1norm
      double precision v2norm

      v1norm = 0.0D+00
      do dim = 1, dim_num
        v1norm = v1norm + ( p1(dim) - p2(dim) )**2
      end do
      v1norm = sqrt ( v1norm )

      if ( v1norm .eq. 0.0D+00 ) then
        angle_rad_3d = 0.0D+00
        return
      end if

      v2norm = 0.0D+00
      do dim = 1, dim_num
        v2norm = v2norm + ( p3(dim) - p2(dim) )**2
      end do
      v2norm = sqrt ( v2norm )

      if ( v2norm .eq. 0.0D+00 ) then
        angle_rad_3d = 0.0D+00
        return
      end if

      dot = 0.0D+00
      do dim = 1, dim_num
        dot = dot + ( p1(dim) - p2(dim) ) * ( p3(dim) - p2(dim) )
      end do

      angle_rad_3d = r8_acos ( dot / ( v1norm * v2norm ) )

      return
      end
      function angle_rad_nd ( dim_num, v1, v2 )

c*********************************************************************72
c
cc ANGLE_RAD_ND returns the angle in radians between two rays in ND.
c
c  Discussion:
c
c    This routine always computes the SMALLER of the two angles between
c    two rays.  Thus, if the rays make an (exterior) angle of 1.5 PI,
c    then the (interior) angle of 0.5 PI is reported.
c
c    X dot Y = Norm(X) * Norm(Y) * Cos( Angle(X,Y) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, int DIM_NUM, the spatial dimension.
c
c    Input, double precision V1(DIM_NUM), V2(DIM_NUM), the two rays.
c
c    Output, double precision ANGLE_RAD_ND, the angle between the rays,
c    in radians.  This value will always be between 0 and PI.
c
      implicit none

      integer dim_num

      double precision angle_rad_nd
      double precision dot
      double precision r8_acos
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision v1(dim_num)
      double precision v1norm
      double precision v2(dim_num)
      double precision v2norm

      dot = r8vec_dot_product ( dim_num, v1, v2 )

      v1norm = r8vec_norm ( dim_num, v1 );

      if ( v1norm .eq. 0.0D+00 ) then
        angle_rad_nd = 0.0D+00
        return
      end if

      v2norm = r8vec_norm ( dim_num, v2 )

      if ( v2norm .eq. 0.0D+00 ) then
        angle_rad_nd = 0.0D+00
        return
      end if

      angle_rad_nd = r8_acos ( dot / ( v1norm * v2norm ) )

      return
      end
      subroutine angle_turn_2d ( p1, p2, p3, turn )

c*********************************************************************72
c
cc ANGLE_TURN_2D computes a turning angle in 2D.
c
c  Discussion:
c
c    This routine is most useful when considering the vertices of a
c    polygonal shape.  We wish to distinguish between angles that "turn
c    in" to the shape, (between 0 and 180 degrees) and angles that
c    "turn out" (between 180 and 360 degrees), as we traverse the boundary.
c
c    If we compute the interior angle and subtract 180 degrees, we get the
c    supplementary angle, which has the nice property that it is
c    negative for "in" angles and positive for "out" angles, and is zero if
c    the three points actually lie along a line.
c
c    Assuming P1, P2 and P3 define an angle, the TURN can be
c    defined to be either:
c
c    * the supplementary angle to the angle formed by P1=P2=P3, or
c
c    * the angle between the vector ( P3-P2) and the vector -(P1-P2),
c      where -(P1-P2) can be understood as the vector that continues
c      through P2 from the direction P1.
c
c    The turning will be zero if P1, P2 and P3 lie along a straight line.
c
c    It will be a positive angle if the turn from the previous direction
c    is counter clockwise, and negative if it is clockwise.
c
c    The turn is given in radians, and will lie between -PI and PI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), P3(2), the points that form
c    the angle.
c
c    Output, double precision TURN, the turn angle, between -PI and PI.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision p(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision r8_atan
      double precision turn

      p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) )
     &     + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

      p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) )
     &     - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

      if ( p(1) .eq. 0.0D+00 .and. p(2) .eq. 0.0D+00 ) then
        turn = 0.0D+00
      else
        turn = pi - r8_atan ( p(2), p(1) )
      end if

      return
      end
      subroutine annulus_area_2d ( r1, r2, area )

c*********************************************************************72
c
cc ANNULUS_AREA_2D computes the area of a circular annulus in 2D.
c
c  Discussion:
c
c    A circular annulus with center (XC,YC), inner radius R1 and
c    outer radius R2, is the set of points (X,Y) so that
c
c      R1**2 .le. (X-XC)**2 + (Y-YC)**2 .le. R2**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 August 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the inner and outer radii.
c
c    Output, double precision AREA, the area.
c
      implicit none

      double precision area
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2

      area = pi * ( r2 + r1 ) * ( r2 - r1 )

      return
      end
      subroutine annulus_sector_area_2d ( r1, r2, theta1, theta2, area )

c*********************************************************************72
c
cc ANNULUS_SECTOR_AREA_2D computes the area of an annular sector in 2D.
c
c  Discussion:
c
c    An annular sector with center PC, inner radius R1 and
c    outer radius R2, and angles THETA1, THETA2, is the set of points
c    P so that
c
c      R1**2 .le. (P(1)-PC(1))**2 + (P(2)-PC(2))**2 .le. R2**2
c
c    and
c
c      THETA1 .le. THETA ( P - PC ) .le. THETA2
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the inner and outer radii.
c
c    Input, double precision THETA1, THETA2, the angles.
c
c    Output, double precision AREA, the area.
c
      implicit none

      double precision area
      double precision r1
      double precision r2
      double precision theta1
      double precision theta2

      area = 0.5D+00 * ( theta2 - theta1 ) * ( r2 + r1 ) * ( r2 - r1 )

      return
      end
      subroutine annulus_sector_centroid_2d ( pc, r1, r2, theta1,
     &  theta2, centroid )

c*********************************************************************72
c
cc ANNULUS_SECTOR_CENTROID_2D computes the centroid of an annular sector in 2D.
c
c  Discussion:
c
c    An annular sector with center PC, inner radius R1 and
c    outer radius R2, and angles THETA1, THETA2, is the set of points
c    P so that
c
c      R1**2 <= (P(1)-PC(1))**2 + (P(2)-PC(2))**2 <= R2**2
c
c    and
c
c      THETA1 <= THETA ( P - PC ) <= THETA2
c
c    Thanks to Ed Segall for pointing out a mistake in the computation
c    of the angle THETA associated with the centroid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 December 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Harris, Horst Stocker,
c    Handbook of Mathematics and Computational Science,
c    Springer, 1998, QA40.S76
c
c  Parameters:
c
c    Input, double precision PC(2), the center.
c
c    Input, double precision R1, R2, the inner and outer radii.
c
c    Input, double precision THETA1, THETA2, the angles.
c
c    Output, double precision CENTROID(2), the centroid.
c
      implicit none

      double precision centroid(2)
      double precision pc(2)
      double precision r
      double precision r1
      double precision r2
      double precision theta
      double precision theta1
      double precision theta2

      theta = theta2 - theta1

      r = 4.0D+00 * sin ( theta / 2.0D+00 ) / ( 3.0D+00 * theta )
     &  * ( r1 * r1 + r1 * r2 + r2 * r2 ) / ( r1 + r2 )

      centroid(1) = pc(1) + r * cos ( theta1 + theta / 2.0D+00 )
      centroid(2) = pc(2) + r * sin ( theta1 + theta / 2.0D+00 )

      return
      end
      subroutine ball_unit_sample_2d ( seed, p )

c*********************************************************************72
c
cc BALL_UNIT_SAMPLE_2D picks a random point in the unit ball in 2D.
c
c  Discussion:
c
c    The unit ball is the set of points P such that
c
c      P(1) * P(1) + P(2) * P(2) <= 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precisionP(2), a random point in the unit ball.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision p(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      integer seed
      double precision theta
      double precision u(dim_num)

      call r8vec_uniform_01 ( dim_num, seed, u )

      r = sqrt ( u(1) )
      theta = 2.0D+00 * pi * u(2)

      p(1) = r * cos ( theta )
      p(2) = r * sin ( theta )

      return
      end
      subroutine ball_unit_sample_3d ( seed, p )

c*********************************************************************72
c
cc BALL_UNIT_SAMPLE_3D picks a random point in the unit ball in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision P(3), the sample point.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision p(dim_num)
      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_acos
      integer seed
      double precision theta
      double precision u(dim_num)
      double precision vdot

      call r8vec_uniform_01 ( dim_num, seed, u )
c
c  Pick a uniformly random VDOT, which must be between -1 and 1.
c  This represents the dot product of the random vector with the Z unit vector.
c
c  Note: this works because the surface area of the sphere between
c  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
c  a patch of area uniformly.
c
      vdot = 2.0D+00 * u(1) - 1.0D+00

      phi = r8_acos ( vdot )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the Z vector.
c
      theta = 2.0D+00 * pi * u(2)
c
c  Pick a random radius R.
c
      r = u(3)**( 1.0D+00 / 3.0D+00 )

      p(1) = r * cos ( theta ) * sin ( phi )
      p(2) = r * sin ( theta ) * sin ( phi )
      p(3) = r * cos ( phi )

      return
      end
      subroutine ball_unit_sample_nd ( dim_num, seed, p )

c*********************************************************************72
c
cc BALL_UNIT_SAMPLE_ND picks a random point in the unit ball in ND.
c
c  Discussion:
c
c    N-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
c
c    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
c    and has the form:
c
c     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
c     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
c
c    Finally, a scaling is applied to set the point at a distance R
c    from the origin, in a way that results in a uniform distribution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision P(N), the random point.
c
      implicit none

      integer dim_num

      double precision r8_uniform_01
      integer i
      double precision p(dim_num)
      double precision pi
      double precision r
      double precision random_cosine
      double precision random_sign
      double precision random_sine
      integer seed

      p(1) = 1.0D+00
      do i = 2, dim_num
        p(i) = 0.0D+00
      end do

      do i = 1, dim_num-1

        r = r8_uniform_01 ( seed )
        random_cosine = 2.0D+00 * r - 1.0D+00
        r = r8_uniform_01 ( seed )
        random_sign = dble ( 2 * int ( 2.0D+00 * r ) - 1 )
        r = r8_uniform_01 ( seed )
        random_sine = random_sign
     &    * sqrt ( 1.0D+00 - random_cosine * random_cosine )

        pi = p(i)
        p(i)   = random_cosine * pi
        p(i+1) = random_sine   * pi

      end do

      r = r8_uniform_01 ( seed )

      r = r**( 1.0D+00 / dble ( dim_num ) )

      do i = 1, dim_num
        p(i) = r * p(i)
      end do

      return
      end
      subroutine basis_map_3d ( u, v, a, ierror )

c*********************************************************************72
c
cc BASIS_MAP_3D computes the matrix which maps one basis to another in 3D.
c
c  Discussion:
c
c    As long as the column vectors U1, U2 and U3 are linearly independent,
c    a matrix A will be computed that maps U1 to V1, U2 to V2, and
c    U3 to V3, where V1, V2 and V3 are the columns of V.
c
c    Depending on the values of the vectors, A may represent a
c    rotation, reflection, dilation, projection, or a combination of these
c    basic linear transformations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision U(3,3), the columns of U are the three
c    "domain" or "preimage" vectors, which should be linearly independent.
c
c    Input, double precision V(3,3), the columns of V are the three
c    "range" or "image" vectors.
c
c    Output, double precision A(3,3), a matrix with the property that
c    A * U1 = V1, A * U2 = V2 and A * U3 = V3.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    nonzero, the matrix [ U1 | U2 | U3 ] is exactly singular.
c
      implicit none

      double precision a(3,3)
      double precision b(3,3)
      double precision c(3,3)
      double precision det
      integer ierror
      double precision u(3,3)
      double precision v(3,3)

      ierror = 0
c
c  Compute C = the inverse of [ U1 | U2 | U3 ].
c
      call r8mat_copy ( 3, 3, u, b )

      call r8mat_inverse_3d ( b, c, det )

      if ( det .eq. 0.0D+00 ) then
        ierror = 1
        return
      end if
c
c  A = [ V1 | V2 | V3 ] * inverse [ U1 | U2 | U3 ].
c
      call r8mat_mm ( 3, 3, 3, v, c, a )

      return
      end
      function box_01_contains_point_2d ( p )

c*********************************************************************72
c
cc BOX_01_CONTAINS_POINT_2D determines if a point is inside the unit box in 2D.
c
c  Discussion:
c
c    A unit box is assumed to be a rectangle with sides aligned on coordinate
c    axes.  It can be described as the set of points P satisfying:
c
c      0.0 <= P(1:DIM_NUM) <= 1.0
c
c      0.0 <= P(1:2) <= 1.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical BOX_01_CONTAINS_POINT_2D, is TRUE if the point is
c    inside the box.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical box_01_contains_point_2d
      integer dim
      double precision p(2)

      box_01_contains_point_2d = .false.

      do dim = 1, dim_num

        if ( p(dim) .lt. 0.0D+00 ) then
          return
        end if

        if ( 1.0D+00 .lt. p(dim) ) then
          return
        end if

      end do

      box_01_contains_point_2d = .true.

      return
      end
      function box_01_contains_point_nd ( dim_num, p )

c*********************************************************************72
c
cc BOX_01_CONTAINS_POINT_ND determines if a point is inside the unit box in ND.
c
c  Discussion:
c
c    A unit box is assumed to be a rectangle with sides aligned on coordinate
c    axes.  It can be described as the set of points P satisfying:
c
c      0.0 <= P(1:DIM_NUM) <= 1.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision P(DIM_NUM), the point to be checked.
c
c    Output, logical BOX_01_CONTAINS_POINT_ND, is TRUE if the point is
c    inside the box.
c
      implicit none

      integer dim_num

      logical box_01_contains_point_nd
      integer dim
      double precision p(dim_num)

      box_01_contains_point_nd = .false.

      do dim = 1, dim_num

        if ( p(dim) .lt. 0.0D+00 ) then
          return
        end if

        if ( 1.0D+00 .lt. p(dim) ) then
          return
        end if

      end do

      box_01_contains_point_nd = .true.

      return
      end
      function box_contains_point_2d ( p1, p2, p )

c*********************************************************************72
c
cc BOX_CONTAINS_POINT_2D determines if a point is inside a box in 2D.
c
c  Discussion:
c
c    A box in 2D is a rectangle with sides aligned on coordinate
c    axes.  It can be described by its low and high corners, P1 and P2
c    as the set of points P satisfying:
c
c      P1(1:2) <= P(1:2) <= P2(1:2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), the low and high
c    corners of the box.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical BOX_CONTAINS_POINT_2D, is TRUE if the point
c    is inside the box.
c
      implicit none

      logical box_contains_point_2d
      double precision p(2)
      double precision p1(2)
      double precision p2(2)

      if ( p(1)  .lt. p1(1) .or.
     &     p2(1) .lt. p(1)  .or.
     &     p(2)  .lt. p1(2) .or.
     &     p2(2) .lt. p(2) ) then
        box_contains_point_2d = .false.
      else
        box_contains_point_2d = .true.
      end if

      return
      end
      function box_contains_point_nd ( dim_num, p1, p2, p )

c*********************************************************************72
c
cc BOX_CONTAINS_POINT_ND determines if a point is inside a box in ND.
c
c  Discussion:
c
c    A box is a rectangle with sides aligned on coordinate
c    axes.  It can be described by its low and high corners, P1 and P2
c    as the set of points P satisfying:
c
c      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision P1(DIM_NUM), P2(DIM_NUM), the low and high
c    corners of the box.
c
c    Input, double precision P(DIM_NUM), the point to be checked.
c
c    Output, logical BOX_CONTAINS_POINT_ND, is TRUE if the point
c    is inside the box.
c
      implicit none

      integer dim_num

      logical box_contains_point_nd
      integer i
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)

      box_contains_point_nd = .false.

      do i = 1, dim_num
        if ( p(i) .lt. p1(i) .or. p2(i) .lt. p(i) ) then
          return
        end if
      end do

      box_contains_point_nd = .true.

      return
      end
      function box_contains_segment_nd ( dim_num, p1, p2, pa, pb  )

c*********************************************************************72
c
cc BOX_CONTAINS_SEGMENT_ND reports if a box contains a line segment in ND.
c
c  Discussion:
c
c    A box is assumed to be a rectangle with sides aligned on coordinate
c    axes.  It can be described by its low and high corners, P1 and P2
c    as the set of points P satisfying:
c
c      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
c
c    A line segment is the finite portion of a line that lies between
c    two points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision P1(DIM_NUM), P2(DIM_NUM), the low and high corners
c    of the box.
c
c    Input, double precision PA(DIM_NUM), PB(DIM_NUM), the endpoints of the
c    line segment.
c
c    Output, logical BOX_CONTAINS_SEGMENT_ND, is TRUE if the box contains
c    the line segment.
c
      implicit none

      integer dim_num

      logical box_contains_segment_nd
      logical box_contains_point_nd
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pa(dim_num)
      double precision pb(dim_num)

      box_contains_segment_nd = .false.

      if ( .not. box_contains_point_nd ( dim_num, p1, p2, pa ) ) then
        return
      end if

      if ( .not. box_contains_point_nd ( dim_num, p1, p2, pb ) ) then
        return
      end if

      box_contains_segment_nd = .true.

      return
      end
      subroutine box_ray_int_2d ( p1, p2, pa, pb, pint )

c*********************************************************************72
c
cc BOX_RAY_INT_2D: intersection ( box, ray ) in 2D.
c
c  Discussion:
c
c    A box in 2D is a rectangle with sides aligned on coordinate
c    axes.  It can be described by its low and high corners, P1 and P2
c    as the set of points P satisfying:
c
c      P1(1:2) <= P(1:2) <= P2(1:2).
c
c    The origin of the ray is assumed to be inside the box.  This
c    guarantees that the ray will intersect the box in exactly one point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), the low and high corners of the box.
c
c    Input, double precision PA(2), the origin of the ray, which should be
c    inside the box.
c
c    Input, double precision PB(2), a second point on the ray.
c
c    Output, double precision PINT(2), the point on the box intersected 
c    by the ray.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical inside
      integer ival
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pa(dim_num)
      double precision pb(dim_num)
      double precision pc(dim_num)
      double precision pd(dim_num)
      double precision pint(dim_num)
      integer side

      do side = 1, 4

        if ( side .eq. 1 ) then
          pd(1) = p1(1)
          pd(2) = p1(2)
          pc(1) = p2(1)
          pc(2) = p1(2)
        else if ( side .eq. 2 ) then
          pd(1) = p2(1)
          pd(2) = p1(2)
          pc(1) = p2(1)
          pc(2) = p2(2)
        else if ( side .eq. 3 ) then
          pd(1) = p2(1)
          pd(2) = p2(2)
          pc(1) = p1(1)
          pc(2) = p2(2)
        else if ( side .eq. 4 ) then
          pd(1) = p1(1)
          pd(2) = p2(2)
          pc(1) = p1(1)
          pc(2) = p1(2)
        end if

        call angle_contains_point_2d ( pc, pa, pd, pb, inside )

        if ( inside ) then
          go to 10
        end if

        if ( side .eq. 4 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BOX_RAY_INT_2D - Fatal error!'
          write ( *, '(a)' ) '  No intersection could be found.'
          stop
        end if

      end do

10    continue

      call lines_exp_int_2d ( pa, pb, pc, pd, ival, pint )

      return
      end
      subroutine box_segment_clip_2d ( p1, p2, pa, pb, ival )

c*********************************************************************72
c
cc BOX_SEGMENT_CLIP_2D uses a box to clip a line segment in 2D.
c
c  Discussion:
c
c    A box in 2D is a rectangle with sides aligned on coordinate
c    axes.  It can be described by its low and high corners, P1 and P2
c    as the set of points P satisfying:
c
c      P1(1:2) .le. P(1:2) .le. P2(1:2).
c
c    A line segment is the finite portion of a line that lies between
c    two points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), the low and high corners of the box.
c
c    Input/output, double precision PA(2), PB(2); on input, the endpoints 
c    of a line segment.  On output, the endpoints of the portion of the
c    line segment that lies inside the box.  However, if no part of the
c    initial line segment lies inside the box, the output value is the
c    same as the input value.
c
c    Output, integer IVAL:
c    -1, no part of the line segment is within the box.
c     0, no clipping was necessary.
c     1, PA was clipped.
c     2, PB was clipped.
c     3, PA and PB were clipped.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical clip_a
      logical clip_b
      integer ival
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pa(dim_num)
      double precision pb(dim_num)
      double precision q(dim_num)

      clip_a = .false.
      clip_b = .false.
c
c  Require that XMIN .le. X.
c
      if ( pa(1) .lt. p1(1) .and. pb(1) .lt. p1(1) ) then
        ival = -1
        return
      end if

      if ( pa(1) .lt. p1(1) .and. p1(1) .le. pb(1) ) then
        q(1) = p1(1)
        q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) 
     &    / ( pb(1) - pa(1) )
        pa(1) = q(1)
        pa(2) = q(2)
        clip_a = .true.
      else if ( p1(1) .le. pa(1) .and. pb(1) .lt. p1(1) ) then
        q(1) = p1(1)
        q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) 
     &    / ( pb(1) - pa(1) )
        pb(1) = q(1)
        pb(2) = q(2)
        clip_b = .true.
      end if
c
c  Require that X .le. XMAX.
c
      if ( p2(1) .lt. pa(1) .and. p2(1) .lt. pb(1) ) then
        ival = -1
        return
      end if

      if ( p2(1) .lt. pa(1) .and. pb(1) .le. p2(1) ) then
        q(1) = p2(1)
        q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) 
     &    / ( pb(1) - pa(1) )
        pa(1) = q(1)
        pa(2) = q(2)
        clip_a = .true.
      else if ( pa(1) .le. p2(1) .and. p2(1) .lt. pb(1) ) then
        q(1) = p2(1)
        q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) 
     &    / ( pb(1) - pa(1) )
        pb(1) = q(1)
        pb(2) = q(2)
        clip_b = .true.
      end if
c
c  Require that YMIN .le. Y.
c
      if ( pa(2) .lt. p1(2) .and. pb(2) .lt. p1(2) ) then
        ival = -1
        return
      end if

      if ( pa(2) .lt. p1(2) .and. p1(2) .le. pb(2) ) then
        q(2) = p1(2)
        q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) 
     &    / ( pb(2) - pa(2) )
        pa(1) = q(2)
        pa(2) = q(2)
        clip_a = .true.
      else if ( p1(2) .le. pa(2) .and. pb(2) .lt. p1(2) ) then
        q(2) = p1(2)
        q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) 
     &    / ( pb(2) - pa(2) )
        pb(1) = q(1)
        pb(2) = q(2)
        clip_b = .true.
      end if
c
c  Require that Y .le. YMAX.
c
      if ( p2(2) .lt. pa(2) .and. p2(2) .lt. pb(2) ) then
        ival = -1
        return
      end if

      if ( p2(2) .lt. pa(2) .and. pb(2) .le. p2(2) ) then
        q(2) = p2(2)
        q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) 
     &    / ( pb(2) - pa(2) )
        pa(1) = q(1)
        pa(2) = q(2)
        clip_a = .true.
      else if ( pa(2) .le. p2(2) .and. p2(2) .lt. pb(2) ) then
        q(2) = p2(2)
        q(1) = pa(1) + ( pb(1) - pa(1) ) * ( p2(2) - pa(2) ) 
     &    / ( pb(2) - pa(2) )
        pb(1) = q(1)
        pb(2) = q(2)
        clip_b = .true.
      end if

      ival = 0

      if ( clip_a ) then
        ival = ival + 1
      end if

      if ( clip_b ) then
        ival = ival + 2
      end if

      return
      end
      subroutine circle_arc_point_near_2d ( r, pc, theta1, theta2, 
     &  p, pn, dist )

c*********************************************************************72
c
cc CIRCLE_ARC_POINT_NEAR_2D : nearest point on a circular arc.
c
c  Discussion:
c
c    A circular arc is defined by the portion of a circle (R,C)
c    between two angles (THETA1,THETA2).
c
c    Thus, a point P on a circular arc satisfies
c
c      ( P(1) - PC(1) ) * ( P(1) - PC(1) ) 
c    + ( P(2) - PC(2) ) * ( P(2) - PC(2) ) = R * R
c
c    and
c
c      Theta1 <= Theta <= Theta2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision PN(2), a point on the circular arc which is
c    nearest to the point.
c
c    Output, double precision DIST, the distance to the nearest point.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision pn(dim_num)
      double precision r
      double precision r2
      double precision r8_atan
      double precision r8_modp
      double precision theta
      double precision theta1
      double precision theta2
c
c  Special case, the zero circle.
c
      if ( r .eq. 0.0D+00 ) then
        pn(1) = pc(1)
        pn(2) = pc(2)
        dist = sqrt ( ( p(1) - pn(1) )**2 + ( p(2) - pn(2) )**2 )
        return
      end if
c
c  Determine the angle made by the point.
c
      theta = r8_atan ( p(2) - pc(2), p(1) - pc(1) )
c
c  If the angle is between THETA1 and THETA2, then you can
c  simply project the point onto the arc.
c
      if ( r8_modp ( theta  - theta1,  2.0D+00 * pi ) .le. 
     &     r8_modp ( theta2 - theta1,  2.0D+00 * pi ) ) then

        r2 = sqrt ( ( p(1) - pc(1) )**2 + ( p(2) - pc(2) )**2 )

        pn(1) = pc(1) + ( p(1) - pc(1) ) * r / r2
        pn(2) = pc(2) + ( p(2) - pc(2) ) * r / r2
c
c  Otherwise, if the angle is less than the negative of the
c  average of THETA1 and THETA2, it's on the side of the arc
c  where the endpoint associated with THETA2 is closest.
c
      else if ( r8_modp ( theta - 0.5D+00 * ( theta1 + theta2 ), 
     &  2.0D+00 * pi ) .le. pi ) then

        pn(1) = pc(1) + r * cos ( theta2 )
        pn(2) = pc(2) + r * sin ( theta2 )
c
c  Otherwise, the endpoint associated with THETA1 is closest.
c
      else

        pn(1) = pc(1) + r * cos ( theta1 )
        pn(2) = pc(2) + r * sin ( theta1 )

      end if

      dist = sqrt ( ( p(1) - pn(1) )**2 + ( p(2) - pn(2) )**2 )

      return
      end
      subroutine circle_area_2d ( r, area )

c*********************************************************************72
c
cc CIRCLE_AREA_2D computes the area of a circle in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Output, double precision AREA, the area of the circle.
c
      implicit none

      double precision area
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      area = pi * r * r

      return
      end
      subroutine circle_dia2imp_2d ( p1, p2, r, pc )

c*********************************************************************72
c
cc CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
c
c  Discussion:
c
c    The diameter form of a circle is:
c
c      P1 and P2 are the endpoints of a diameter.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points that are the
c    endpoints of a diameter of the circle.
c
c    Output, double precision R, the radius of the circle.
c
c    Output, double precision PC(2), the center of the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      integer dim
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pc(dim_num)
      double precision r

      r = 0.0D+00
      do dim = 1, dim_num
        r = r + ( p2(dim) - p1(dim) )**2
      end do
      r = 0.5D+00 * sqrt ( r )

      do dim = 1, dim_num
        pc(dim) = 0.5D+00 * ( p1(dim) + p2(dim) )
      end do

      return
      end
      subroutine circle_exp_contains_point_2d ( p1, p2, p3, p, inside )

c*********************************************************************72
c
cc CIRCLE_EXP_CONTAINS_POINT_2D: explicit circle contains a point in 2D.
c
c  Discussion:
c
c    The explicit form of a circle in 2D is:
c
c      The circle passing through points P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), P3(2), three points on a circle.
c
c    Input, double precision P(2), the point to test.
c
c    Output, integer INSIDE, reports the result:
c   -1, the three points are distinct and noncolinear,
c    and P lies inside the circle.
c    0, the three points are distinct and noncolinear,
c    and P lies on the circle.
c    1, the three points are distinct and noncolinear,
c    and P lies outside the circle.
c    2, the three points are distinct and colinear,
c    and P lies on the line.
c    3, the three points are distinct and colinear,
c    and P does not lie on the line.
c    4, two points are distinct, and P lies on the line.
c    5, two points are distinct, and P does not lie on the line.
c    6, all three points are equal, and P is equal to them,
c    7, all three points are equal, and P is not equal to them.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a(4,4)
      double precision det
      double precision r8mat_det_4d
      logical r8vec_eq
      integer inside
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
c
c  P1 = P2?
c
      if ( r8vec_eq ( dim_num, p1, p2 ) ) then

        if ( r8vec_eq ( dim_num, p1, p3 ) ) then

          if ( r8vec_eq ( dim_num, p1, p ) ) then
            inside = 6
          else
            inside = 7
          end if

        else

          det = ( p1(1) - p3(1) ) * ( p(2)  - p3(2) ) 
     &        - ( p(1)  - p3(1) ) * ( p1(2) - p3(2) )

          if ( det .eq. 0.0D+00 ) then
            inside = 4
          else
            inside = 5
          end if
        end if

        return

      end if
c
c  P1 does not equal P2.  Does P1 = P3?
c
      if ( r8vec_eq ( dim_num, p1, p3 ) ) then

        det = ( p1(1) - p2(1) ) * ( p(2)  - p2(2) ) 
     &      - ( p(1)  - p2(1) ) * ( p1(2) - p2(2) )

        if ( det .eq. 0.0D+00 ) then
          inside = 4
        else
          inside = 5
        end if

        return

      end if
c
c  The points are distinct.  Are they colinear?
c
      det = ( p1(1) - p2(1) ) * ( p3(2) - p2(2) ) 
     &    - ( p3(1) - p2(1) ) * ( p1(2) - p2(2) )

      if ( det .eq. 0.0D+00 ) then

        det = ( p1(1) - p2(1) ) * ( p(2)  - p2(2) ) 
     &      - ( p(1)  - p2(1) ) * ( p1(2) - p2(2) )

        if ( det .eq. 0.0D+00 ) then
          inside = 2
        else
          inside = 3
        end if

        return

      end if
c
c  The points are distinct and non-colinear.
c
c  Compute the determinant
c
      a(1,1) = p1(1)
      a(1,2) = p1(2)
      a(1,3) = p1(1) * p1(1) + p1(2) * p1(2)
      a(1,4) = 1.0D+00

      a(2,1) = p2(1)
      a(2,2) = p2(2)
      a(2,3) = p2(1) * p2(1) + p2(2) * p2(2)
      a(2,4) = 1.0D+00

      a(3,1) = p3(1)
      a(3,2) = p3(2)
      a(3,3) = p3(1) * p3(1) + p3(2) * p3(2)
      a(3,4) = 1.0D+00

      a(4,1) = p(1)
      a(4,2) = p(2)
      a(4,3) = p(1) * p(1) + p(2) * p(2)
      a(4,4) = 1.0D+00

      det = r8mat_det_4d ( a )

      if ( det .lt. 0.0D+00 ) then
        inside = 1
      else if ( det .eq. 0.0D+00 ) then
        inside = 0
      else
        inside = -1
      end if

      return
      end
      subroutine circle_exp2imp_2d ( p1, p2, p3, r, pc )

c*********************************************************************72
c
cc CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
c
c  Discussion:
c
c    The explicit form of a circle in 2D is:
c
c      The circle passing through points P1, P2 and P3.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    Any three distinct points define a circle, as long as they don't lie 
c    on a straight line.  (If the points do lie on a straight line, we 
c    could stretch the definition of a circle to allow an infinite radius 
c    and a center at some infinite point.)
c
c    The diameter of the circle can be found by solving a 2 by 2 linear system.
c    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
c    and each forms a right triangle with the diameter.  Hence, the dot product
c    of P2 - P1 with the diameter is equal to the square of the length
c    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
c    diameter vector originating at P1.
c
c    If all three points are equal, return a circle of radius 0 and 
c    the obvious center.
c
c    If two points are equal, return a circle of radius half the distance
c    between the two distinct points, and center their average.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Joseph ORourke,
c    Computational Geometry,
c    Second Edition,
c    Cambridge, 1998,
c    ISBN: 0521649765,
c    LC: QA448.D38.
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), P3(2), three points on the circle.
c
c    Output, double precision R, the radius of the circle.  Normally, R will
c    be positive.  R will be (meaningfully) zero if all three points are 
c    equal.  If two points are equal, R is returned as the distance between
c    two nonequal points.  R is returned as -1 in the unlikely event that 
c    the points are numerically collinear; philosophically speaking, R 
c    should actually be "infinity" in this case.
c
c    Output, double precision PC(2), the center of the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision e
      double precision f
      double precision g
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pc(dim_num)
      double precision r
      logical r8vec_eq
      double precision r8vec_diff_norm
c
c  If all three points are equal, then the
c  circle of radius 0 and center P1 passes through the points.
c
      if ( r8vec_eq ( dim_num, p1, p2 ) .and. 
     &     r8vec_eq ( dim_num, p1, p3 ) ) then
        r = 0.0D+00
        pc(1) = p1(1)
        pc(2) = p1(2)
        return
      end if
c
c  If exactly two points are equal, then the circle is defined as
c  having the obvious radius and center.
c
           if ( r8vec_eq ( dim_num, p1, p2 ) ) then

        r = 0.5D+00 * r8vec_diff_norm ( dim_num, p1, p3 )
        pc(1) = 0.5D+00 * ( p1(1) + p3(1) )
        pc(2) = 0.5D+00 * ( p1(2) + p3(2) )
        return

      else if ( r8vec_eq ( dim_num, p1, p3 ) ) then

        r = 0.5D+00 * r8vec_diff_norm ( dim_num, p1, p2 )
        pc(1) = 0.5D+00 * ( p1(1) + p2(1) )
        pc(2) = 0.5D+00 * ( p1(2) + p2(2) )
        return

      else if ( r8vec_eq ( dim_num, p2, p3 ) ) then

        r = 0.5D+00 * r8vec_diff_norm ( dim_num, p1, p2 )
        pc(1) = 0.5D+00 * ( p1(1) + p2(1) )
        pc(2) = 0.5D+00 * ( p1(2) + p2(2) )
        return

      end if
c
c  We check for collinearity.  A more useful check would compare the
c  absolute value of G to a small quantity.
c
      e = ( p2(1) - p1(1) ) * ( p1(1) + p2(1) ) 
     &  + ( p2(2) - p1(2) ) * ( p1(2) + p2(2) )

      f = ( p3(1) - p1(1) ) * ( p1(1) + p3(1) ) 
     &  + ( p3(2) - p1(2) ) * ( p1(2) + p3(2) )

      g = ( p2(1) - p1(1) ) * ( p3(2) - p2(2) ) 
     &  - ( p2(2) - p1(2) ) * ( p3(1) - p2(1) )

      if ( g .eq. 0.0D+00 ) then
        pc(1) = 0.0D+00
        pc(2) = 0.0D+00
        r = -1.0D+00
        return
      end if
c
c  The center is halfway along the diameter vector from P1.
c
      pc(1) = 0.5D+00 * ( ( p3(2) - p1(2) ) * e 
     &                  - ( p2(2) - p1(2) ) * f ) / g
      pc(2) = 0.5D+00 * ( ( p2(1) - p1(1) ) * f 
     &                  - ( p3(1) - p1(1) ) * e ) / g
c
c  Knowing the center, the radius is now easy to compute.
c
      r = r8vec_diff_norm ( dim_num, p1, pc )

      return
      end
      subroutine circle_imp_contains_point_2d ( r, pc, p, inside )

c*********************************************************************72
c
cc CIRCLE_IMP_CONTAINS_POINT_2D: implicit circle contains a point in 2D?
c
c  Discussion:
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is inside or on the circle.
c
      implicit none

      integerdim_num
      parameter ( dim_num = 2 )

      logical inside
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision r

      if ( ( p(1) - pc(1) ) * ( p(1) - pc(1) )
     &   + ( p(2) - pc(2) ) * ( p(2) - pc(2) ) .le. r * r ) then
        inside = .true.
      else
        inside = .false.
      end if

      return
      end
      subroutine circle_imp_line_exp_dist_2d ( r, pc, p1, p2, dist )

c*********************************************************************72
c
cc CIRCLE_IMP_LINE_EXP_DIST_2D: distance ( implicit circle, explicit line ) in 2D.
c
c  Discussion:
c
c    The distance is zero if the line intersects the circle.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The distance between the circle and the line is zero if
c    and only if they intersect.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Output, double precision DIST, the distance of the line to the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pc(dim_num)
      double precision r

      call line_exp_point_dist_2d ( p1, p2, pc, dist )

      dist = dist - r

      if ( dist .lt. 0.0D+00 ) then
        dist = 0.0D+00
      end if

      return
      end
      subroutine circle_imp_line_par_int_2d ( r, pc, x0, y0, f, g,
     &  int_num, p )

c*********************************************************************72
c
cc CIRCLE_IMP_LINE_PAR_INT_2D: ( implicit circle, parametric line ) intersection in 2D.
c
c  Discussion:
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F**2 + G**2 = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision F, G, X0, Y0, the parametric parameters of
c    the line.
c
c    Output, integer INT_NUM, the number of intersecting points found.
c    INT_NUM will be 0, 1 or 2.
c
c    Output, double precision P(2,INT_NUM), the intersecting points.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision f
      double precision g
      integer int_num
      double precision p(dim_num,2)
      double precision pc(dim_num)
      double precision r
      double precision root
      double precision t
      double precision x0
      double precision y0

      root = r * r * ( f * f + g * g ) - ( f * ( pc(2) - y0 )
     &  - g * ( pc(1) - x0 ) )**2

      if ( root .lt. 0.0D+00 ) then

        int_num = 0

      else if ( root .eq. 0.0D+00 ) then

        int_num = 1

        t = ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) )
     &    / ( f * f + g * g )
        p(1,1) = x0 + f * t
        p(2,1) = y0 + g * t

      else if ( 0.0D+00 .lt. root ) then

        int_num = 2

        t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) )
     &    - sqrt ( root ) ) / ( f * f + g * g )

        p(1,1) = x0 + f * t
        p(2,1) = y0 + g * t

        t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) )
     &    + sqrt ( root ) ) / ( f * f + g * g )

        p(1,2) = x0 + f * t
        p(2,2) = y0 + g * t

      end if

      return
      end
      subroutine circle_imp_point_dist_2d ( r, pc, p, dist )

c*********************************************************************72
c
cc CIRCLE_IMP_POINT_DIST_2D: distance ( implicit circle, point ) in 2D.
c
c  Discussion:
c
c    The distance is zero if the point is on the circle.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision DIST, the distance of the point to the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision r
      double precision r2
      double precision r8vec_diff_norm

      r2 = r8vec_diff_norm ( 2, p, pc )
 
      dist = abs ( r2 - r )

      return
      end
      subroutine circle_imp_point_dist_signed_2d ( r, pc, p, dist )

c*********************************************************************72
c
cc CIRCLE_IMP_POINT_DIST_SIGNED_2D: signed distance ( imp circle, point ) in 2D.
c
c  Discussion:
c
c    The signed distance is zero if the point is on the circle.
c    The signed distance is negative if the point is inside the circle.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision DIST, the signed distance of the point
c    to the circle.  If the point is inside the circle, the signed distance
c    is negative.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision r
      double precision r2
      double precision r8vec_diff_norm

      r2 = r8vec_diff_norm ( dim_num, p, pc )

      dist = r2 - r

      return
      end
      subroutine circle_imp_point_near_2d ( r, pc, p, pn, dist )

c*********************************************************************72
c
cc CIRCLE_IMP_POINT_NEAR_2D: nearest ( implicit circle, point ) in 2D.
c
c  Discussion:
c
c    This routine finds the distance from a point to an implicitly
c    defined circle, and returns the point on the circle that is
c    nearest to the given point.
c
c    If the given point is the center of the circle, than any point
c    on the circle is "the" nearest.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision PN(2), the nearest point on the circle.
c
c    Output, double precision DIST, the distance of the point to the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision pn(dim_num)
      double precision r
      double precision r2
      double precision r8vec_diff_norm
      logical r8vec_eq

      if ( r8vec_eq ( dim_num, p, pc ) ) then
        dist = r
        pn(1) = pc(1) + r / sqrt ( dble ( dim_num ) )
        pn(2) = pc(2) + r / sqrt ( dble ( dim_num ) )
        return
      end if

      r2 = r8vec_diff_norm ( dim_num, p, pc )

      dist = abs (  r2 - r )

      pn(1) = pc(1) + r * ( p(1) - pc(1) ) / r2
      pn(2) = pc(2) + r * ( p(2) - pc(2) ) / r2

      return
      end
      subroutine circle_imp_points_2d ( r, pc, n, p )

c*********************************************************************72
c
cc CIRCLE_IMP_POINTS_2D returns points on an implicit circle in 2D.
c
c  Discussion:
c
c    The first point is always ( PC(1) + R, PC(2) ), and subsequent
c    points proceed counter clockwise around the circle.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, integer N, the number of points desired.  
c    N must be at least 1.
c
c    Output, double precision P(2,N), the coordinates of points 
c    on the circle.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer j
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      do j = 1, n
        theta = ( 2.0D+00 * pi * dble ( j - 1 ) ) / dble ( n )
        p(1,j) = pc(1) + r * cos ( theta )
        p(2,j) = pc(2) + r * sin ( theta )
      end do

      return
      end
      subroutine circle_imp_points_3d ( r, pc, nc, n, p )

c*********************************************************************72
c
cc CIRCLE_IMP_POINTS_3D returns points on an implicit circle in 3D.
c
c  Discussion:
c
c    Points P on an implicit circle in 3D satisfy the equations:
c
c      ( P(1) - PC(1) )**2 
c    + ( P(2) - PC(2) )**2 
c    + ( P(3) - PC(3) )**2 = R**2
c
c    and
c
c      ( P(1) - PC(1) ) * NC(1) 
c    + ( P(2) - PC(2) ) * NC(2) 
c    + ( P(3) - PC(3) ) * NC(3) = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(3), the center of the circle.
c
c    Input, double precision NC(3), a nonzero vector that is normal to
c    the plane of the circle.  It is customary, but not necessary,
c    that this vector have unit norm.
c
c    Input, integer N, the number of points desired.  
c    N must be at least 1.
c
c    Output, double precision P(3,N), the coordinates of points 
c    on the circle.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      integer j
      double precision n1(dim_num)
      double precision n2(dim_num)
      double precision nc(dim_num)
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta
c
c  Get two unit vectors N1 and N2 which are orthogonal to each other,
c  and to NC.
c
      call plane_normal_basis_3d ( pc, nc, n1, n2 )
c
c  Rotate R units away from PC in the plane of N1 and N2.
c
      do j = 1, n

        theta = ( 2.0D+00 * pi * dble ( j - 1 ) ) / dble ( n )

        do i = 1, dim_num
          p(i,j) = pc(i) 
     &      + r * ( cos ( theta ) * n1(i) 
     &            + sin ( theta ) * n2(i) )
        end do

      end do

      return
      end
      subroutine circle_imp_points_arc_2d ( r, pc, theta1, theta2, 
     &  n, p )

c*********************************************************************72
c
cc CIRCLE_IMP_POINTS_ARC_2D: N points on an arc of an implicit circle in 2D.
c
c  Discussion:
c
c    The first point is 
c      ( PC(1) + R * COS ( THETA1 ), PC(2) + R * SIN ( THETA1 ) );
c    The last point is
c      ( PC(1) + R * COS ( THETA2 ), PC(2) + R * SIN ( THETA2 ) );
c    and the intermediate points are evenly spaced in angle between these,
c    and in counter clockwise order.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angular coordinates of 
c    the first and last points to be drawn, in radians.
c
c    Input, integer N, the number of points desired.  
c    N must be at least 1.
c
c    Output, double precision P(2,N), the points on the circle.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_modp
      double precision theta
      double precision theta1
      double precision theta2
      double precision theta3
c
c  THETA3 is the smallest angle, no less than THETA1, which
c  coincides with THETA2.
c
      theta3 = theta1 + r8_modp ( theta2 - theta1, 2.0D+00 * pi )

      do i = 1, n

        if ( 1 .lt. n ) then
          theta = ( dble ( n - i     ) * theta1   
     &            + dble (     i - 1 ) * theta3 ) 
     &            / dble ( n     - 1 )
        else
          theta = 0.5D+00 * ( theta1 + theta3 )
        end if

        p(1,i) = pc(1) + r * cos ( theta )
        p(2,i) = pc(2) + r * sin ( theta )

      end do

      return
      end
      subroutine circle_imp_print_2d ( r, pc, title )

c*********************************************************************72
c
cc CIRCLE_IMP_PRINT_2D prints an implicit circle in 2D.
c
c  Discussion:
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision pc(dim_num)
      double precision r
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)'        ) ' '
      write ( *, '(a,g14.6)'  ) '  Radius = ', r
      write ( *, '(a,2g14.6)' ) '  Center = ', pc(1), pc(2)

      return
      end
      subroutine circle_imp_print_3d ( r, pc, nc, title )

c*********************************************************************72
c
cc CIRCLE_IMP_PRINT_3D prints an implicit circle in 3D.
c
c  Discussion:
c
c    Points P on an implicit circle in 3D satisfy the equations:
c
c      ( P(1) - PC(1) )**2 
c    + ( P(2) - PC(2) )**2 
c    + ( P(3) - PC(3) )**2 = R**2
c
c    and
c
c      ( P(1) - PC(1) ) * NC(1) 
c    + ( P(2) - PC(2) ) * NC(2) 
c    + ( P(3) - PC(3) ) * NC(3) = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(3), the center of the circle.
c
c    Input, double precision NC(3), the normal vector to the circle.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision nc(dim_num)
      double precision pc(dim_num)
      double precision r
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)'        ) ' '
      write ( *, '(a,g14.6)'  ) '  Radius = ', r
      write ( *, '(a,3g14.6)' ) '  Center = ', pc(1), pc(2), pc(3)
      write ( *, '(a,3g14.6)' ) '  Normal = ', nc(1), nc(2), nc(3)

      return
      end
      subroutine circle_imp2exp_2d ( r, pc, p1, p2, p3 )

c*********************************************************************72
c
cc CIRCLE_IMP2EXP_2D converts a circle from implicit to explicit form in 2D.
c
c  Discussion:
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    The explicit form of a circle in 2D is:
c
c      The circle passing through points P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Joseph ORourke,
c    Computational Geometry,
c    Second Edition,
c    Cambridge, 1998,
c    ISBN: 0521649765,
c    LC: QA448.D38.
c
c  Parameters:
c
c    Input, double precision R, PC(2), the radius and center of the circle.
c
c    Output, double precision P1(2), P2(2), P3(2), three points on the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pc(dim_num)
      double precision pi 
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      theta = 0.0D+00
      p1(1) = pc(1) + r * cos ( theta )
      p1(2) = pc(2) + r * sin ( theta )

      theta = 2.0D+00 * pi / 3.0D+00
      p2(1) = pc(1) + r * cos ( theta )
      p2(2) = pc(2) + r * sin ( theta )

      theta = 4.0D+00 * pi / 3.0D+00
      p3(1) = pc(1) + r * cos ( theta )
      p3(2) = pc(2) + r * sin ( theta )

      return
      end
      subroutine circle_llr2imp_2d ( p1, p2, q1, q2, r, pc )

c*********************************************************************72
c
cc CIRCLE_LLR2IMP_2D converts a circle from LLR to implicit form in 2D.
c
c  Discussion:
c
c    The LLR form of a circle in 2D is:
c
c      The circle of radius R tangent to the lines L1 and L2.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    Let S be the scaled distance of a point on L1 from P1 to P2, 
c    and let N1 be a unit normal vector to L1.  Then a point P that is
c    R units from L1 satisfies:
c
c      P = P1 + s * ( P2 - P1 ) + R * N1.
c
c    Let t be the scaled distance of a point on L2 from Q1 to Q2,
c    and let N2 be a unit normal vector to L2.  Then a point Q that is
c    R units from L2 satisfies:
c
c      Q = Q1 + t * ( Q2 - Q1 ) + R * N2.
c
c    For the center of the circle, then, we have P = Q, that is
c
c      ( P2 - P1 ) * s - ( Q2 - Q1 ) * t = - P1 + Q1 - R * N1 + R * N2 )
c
c    This is a linear system for ( s and t ) from which we can compute
c    the points of tangency, and the center.
c
c    Note that we have four choices for the circle based on the use
c    of plus or minus N1 and plus or minus N2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on line 1.
c
c    Input, double precision Q1(2), Q2(2), two points on line 2.
c
c    Input, double precision R, the radius of the circle.  
c
c    Output, double precision PC(2,4), the centers of the circles.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a(2,2)
      double precision b(2)
      double precision det
      integer i
      double precision n1(dim_num)
      double precision n2(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pc(dim_num,4)
      double precision q1(dim_num)
      double precision q2(dim_num)
      double precision r
      double precision x(dim_num)
c
c  Compute the normals N1 and N2.
c
      call line_exp_normal_2d ( p1, p2, n1 )

      call line_exp_normal_2d ( q1, q2, n2 )
c
c  Set the linear system.
c
      do i = 1, 2
        a(i,1) =   p2(i) - p1(i)
        a(i,2) = - q2(i) + q1(i)
      end do
c
c  Solve the 4 linear systems, using every combination of 
c  signs on the normal vectors.
c
      do i = 1, 2
        b(i) = - p1(i) + q1(i) + r * n1(i) + r * n2(i)
      end do

      call r8mat_solve_2d ( a, b, det, x )

      do i = 1, 2
        pc(i,1) = p1(i) + ( p2(i) - p1(i) ) * x(1) - r * n1(i) 
        b(i) = - p1(i) + q1(i) + r * n1(i) - r * n2(i)
      end do

      call r8mat_solve_2d ( a, b, det, x )

      do i = 1, 2
        pc(i,2) = p1(i) + ( p2(i) - p1(i) ) * x(1) - r * n1(i) 
        b(i) = - p1(i) + q1(i) - r * n1(i) + r * n2(i)
      end do

      call r8mat_solve_2d ( a, b, det, x )

      do i = 1, 2
        pc(i,3) = p1(i) + ( p2(i) - p1(i) ) * x(1) + r * n1(i) 
        b(i) = - p1(i) + q1(i) - r * n1(i) - r * n2(i)
      end do

      call r8mat_solve_2d ( a, b, det, x )

      do i = 1, 2
        pc(i,4) = p1(i) + ( p2(i) - p1(i) ) * x(1) + r * n1(i) 
      end do

      return
      end
      subroutine circle_lune_area_2d ( r, pc, theta1, theta2, area )

c*********************************************************************72
c
cc CIRCLE_LUNE_AREA_2D returns the area of a circular lune in 2D.
c
c  Discussion:
c
c    A lune is formed by drawing a circular arc, and joining its endpoints.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Output, double precision AREA, the area of the lune.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision area_sector
      double precision area_triangle
      double precision pc(dim_num)
      double precision r
      double precision theta1
      double precision theta2

      call circle_sector_area_2d ( r, pc, theta1, theta2, area_sector )

      call circle_triangle_area_2d ( r, pc, theta1, theta2, 
     &  area_triangle )

      area = area_sector - area_triangle

      return
      end
      subroutine circle_lune_centroid_2d ( r, pc, theta1, theta2, 
     &  centroid )

c*********************************************************************72
c
cc CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
c
c  Discussion:
c
c    A lune is formed by drawing a circular arc, and joining its endpoints.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Output, double precision CENTROID(2), the coordinates of the centroid
c    of the lune.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision centroid(dim_num)
      double precision d
      double precision pc(dim_num)
      double precision r
      double precision theta
      double precision theta1
      double precision theta2

      theta = theta2 - theta1

      if ( theta .eq. 0.0D+00 ) then
        d = r
      else
        d = 4.0D+00 * r * ( sin ( 0.5D+00 * theta ) )**3 / 
     &    ( 3.0D+00 * ( theta - sin ( theta ) ) )
      end if

      centroid(1) = pc(1) + d * cos ( theta )
      centroid(2) = pc(2) + d * sin ( theta )

      return
      end
      subroutine circle_pppr2imp_3d ( p1, p2, p3, r, pc, normal )

c*********************************************************************72
c
cc CIRCLE_PPPR2IMP_3D converts a circle from PPPR to implicit form in 3D.
c
c  Discussion:
c
c    The PPPR form of a circle in 3D is:
c
c      The circle of radius R passing through points P1 and P2,
c      and lying in the plane of P1, P2 and P3.
c
c    Points P on an implicit circle in 2D satisfy the equations:
c
c        ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 + ( P(3) - PC(3) )**2 = R**2
c      and 
c        ( P - PC ) dot NORMAL = 0.
c
c    There may be zero, one, or two circles that satisfy the
c    requirements of the PPPR form.
c
c    If there is no such circle, then PC(1:2,1) and PC(1:2,2)
c    are set to the midpoint of (P1,P2).
c
c    If there is one circle, PC(1:2,1) and PC(1:2,2) will be equal.
c
c    If there are two circles, then PC(1:2,1) is the first center,
c    and PC(1:2,2) is the second.
c
c    This calculation is equivalent to finding the intersections of
c    spheres of radius R at points P1 and P2, which lie in the plane
c    defined by P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), two points on the circle.
c
c    Input, double precision P3(3), a third point.
c
c    Input, double precision R, the radius of the circle.
c
c    Output, double precision PC(3,2), the centers of the two circles.
c
c    Output, double precision NORMAL(3), the normal to the circles.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dist
      double precision dot
      double precision h
      integer i
      integer j
      double precision length
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pc(dim_num,2)
      double precision r
      double precision r8vec_diff_norm
      double precision r8vec_norm
      double precision v(dim_num)
c
c  Compute the distance from P1 to P2.
c
      dist = r8vec_diff_norm ( dim_num, p2, p1 )
c
c  If R is smaller than DIST, we don't have a circle.
c
      if ( 2.0D+00 * r .lt. dist ) then
        do j = 1, 2
          do i = 1, dim_num
            pc(i,j) = 0.5D+00 * ( p1(i) + p2(i) )
          end do
        end do
        return
      end if
c
c  H is the distance from the midpoint of (P1,P2) to the center.
c
      h = sqrt ( ( r + 0.5D+00 * dist ) * ( r - 0.5D+00 * dist ) )
c
c  Define a unit direction V that is normal to P2-P1, and lying
c  in the plane (P1,P2,P3).
c
c  To do this, subtract from P3-P1 the component in the direction P2-P1.
c
      do i = 1, dim_num
        v(i) = p3(i) - p1(i)
      end do

      dot = 0.0D+00
      do i = 1, dim_num
        dot = dot + v(i) * ( p2(i) - p1(i) )
      end do
      dot = dot / dist

      do i = 1, dim_num
        v(i) = v(i) - dot * ( p2(i) - p1(i) ) / dist
      end do

      length = r8vec_norm ( dim_num, v )

      do i = 1, dim_num
        v(i) = v(i) / length
      end do
c
c  We can go with or against the given normal direction.
c
      do i = 1, dim_num
        pc(i,1) = 0.5D+00 * ( p2(i) + p1(i) ) + h * v(i)
      end do

      do i = 1, dim_num
        pc(i,2) = 0.5D+00 * ( p2(i) + p1(i) ) - h * v(i)
      end do

      call plane_exp_normal_3d ( p1, p2, p3, normal )

      return
      end
      subroutine circle_ppr2imp_2d ( p1, p2, r, pc )

c*********************************************************************72
c
cc CIRCLE_PPR2IMP_2D converts a circle from PPR to implicit form in 2D.
c
c  Discussion:
c
c    The PPR form of a circle in 2D is:
c
c      The circle of radius R passing through points P1 and P2.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    There may be zero, one, or two circles that satisfy the 
c    requirements of the PPR form.
c
c    If there is no such circle, then PC(1:2,1) and PC(1:2,2)
c    are set to the midpoint of (P1,P2).
c
c    If there is one circle, PC(1:2,1) and PC(1:2,2) will be equal.
c
c    If there are two circles, then PC(1:2,1) is the first center,
c    and PC(1:2,2) is the second.
c
c    This calculation is equivalent to finding the intersections of
c    circles of radius R at points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the circle.
c
c    Input, double precision R, the radius of the circle.  
c
c    Output, double precision PC(2,2), the centers of the two circles.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision h
      integer i
      integer j
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pc(dim_num,2)
      double precision r
      double precision r8vec_diff_norm
c
c  Compute the distance from P1 to P2.
c
      dist = r8vec_diff_norm ( dim_num, p2, p1 )
c
c  If R is smaller than DIST, we don't have a circle.
c
      if ( 2.0D+00 * r .lt. dist ) then
        do j = 1, 2
          do i = 1, dim_num
            pc(i,j) = 0.5D+00 * ( p1(i) + p2(i) )
          end do
        end do
        return
      end if
c
c  H is the distance from the midpoint of (P1,P2) to the center.
c
      h = sqrt ( ( r + 0.5D+00 * dist ) * ( r - 0.5D+00 * dist ) )
c
c  Determine the unit normal direction.
c
      normal(1) =   ( p2(2) - p1(2) ) / dist
      normal(2) = - ( p2(1) - p1(1) ) / dist
c
c  We can go with or against the given normal direction.
c
      do i = 1, dim_num
        pc(i,1) = 0.5D+00 * ( p2(i) + p1(i) ) + h * normal(i)
      end do

      do i = 1, dim_num
        pc(i,2) = 0.5D+00 * ( p2(i) + p1(i) ) - h * normal(i)
      end do

      return
      end
      subroutine circle_sector_area_2d ( r, pc, theta1, theta2, area )

c*********************************************************************72
c
cc CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
c
c  Discussion:
c
c    A circular sector is formed by a circular arc, and the two straight line
c    segments that join its ends to the center of the circle.
c
c    A circular sector is defined by the two conditions
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    and
c
c      Theta1 <= Theta <= Theta2
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
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the two angles defining the
c    sector, in radians.  Normally, THETA1 < THETA2.
c
c    Output, double precision AREA, the area of the circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision pc(dim_num)
      double precision r
      double precision theta1
      double precision theta2

      area = 0.5D+00 * r * r * ( theta2 - theta1 )

      return
      end
      subroutine circle_sector_centroid_2d ( r, pc, theta1, theta2,
     &  centroid )

c*********************************************************************72
c
cc CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
c
c  Discussion:
c
c    A circular sector is formed by a circular arc, and the two straight line
c    segments that join its ends to the center of the circle.
c
c    A circular sector is defined by
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    and
c
c      Theta1 <= Theta <= Theta2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Output, double precision CENTROID(2), the coordinates of the centroid
c    of the sector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision centroid(dim_num)
      double precision d
      double precision pc(dim_num)
      double precision r
      double precision theta
      double precision theta1
      double precision theta2

      theta = theta2 - theta1

      if ( theta .eq. 0.0D+00 ) then
        d = 2.0D+00 * r / 3.0D+00
      else
        d = 4.0D+00 * r * sin ( 0.5D+00 * theta ) / ( 3.0D+00 * theta )
      end if

      centroid(1) = pc(1) + d * cos ( theta )
      centroid(2) = pc(2) + d * sin ( theta )

      return
      end
      subroutine circle_sector_contains_point_2d ( r, pc, theta1, 
     &  theta2, p, inside )

c*********************************************************************72
c
cc CIRCLE_SECTOR_CONTAINS_POINT_2D : is a point inside a circular sector?
c
c  Discussion:
c
c    A circular sector is formed by a circular arc, and the two straight line 
c    segments that join its ends to the center of the circle.
c
c    A circular sector is defined by
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    and
c
c      Theta1 <= Theta <= Theta2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is inside or on the
c    circular sector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical inside
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_atan
      double precision r8_modp
      double precision theta
      double precision theta1
      double precision theta2

      inside = .false.
c
c  Is the point inside the (full) circle?
c
      if ( ( p(1) - pc(1) ) * ( p(1) - pc(1) ) 
     &   + ( p(2) - pc(2) ) * ( p(2) - pc(2) ) .le. r * r ) then
c
c  Is the point's angle within the arc's range?
c  Try to force the angles to lie between 0 and 2 * PI.
c
        theta = r8_atan ( p(2) - pc(2), p(1) - pc(1) )

        if ( r8_modp ( theta  - theta1,  2.0D+00 * pi ) .le. 
     &       r8_modp ( theta2 - theta1,  2.0D+00 * pi ) ) then

          inside = .true.

        end if

      end if

      return
      end
      subroutine circle_sector_print_2d ( r, pc, theta1, theta2 )

c*********************************************************************72
c
cc CIRCLE_SECTOR_PRINT_2D prints a circular sector in 2D.
c
c  Discussion:
c
c    A circular sector is formed by a circular arc, and the two straight line 
c    segments that join its ends to the center of the circle.
c
c    A circular sector is defined by
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c    and
c
c      Theta1 <= Theta <= Theta2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision pc(dim_num)
      double precision r
      double precision theta1
      double precision theta2

      write ( *, '(a)'        ) ' '
      write ( *, '(a)'        ) '  Circular sector definition:'
      write ( *, '(a)'        ) ' '
      write ( *, '(a,g14.6)'  ) '    Radius = ', r
      write ( *, '(a,2g14.6)' ) '    Center = ', pc(1:2)
      write ( *, '(a,2g14.6)' ) '    Theta  = ', theta1, theta2

      return
      end
      subroutine circle_triangle_area_2d ( r, pc, theta1, theta2, area )

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
c    Note that for angles greater than PI, the triangle will actually
c    have NEGATIVE area.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle.
c
c    Input, double precision PC(2), the center of the circle.
c
c    Input, double precision THETA1, THETA2, the angles defining the arc,
c    in radians.  Normally, THETA1 < THETA2.
c
c    Output, double precision AREA, the (signed) area of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision pc(dim_num)
      double precision r
      double precision theta1
      double precision theta2

      area = 0.5D+00 * r * r * sin ( theta2 - theta1 )

      return
      end
      subroutine circle_triple_angles_2d ( r1, r2, r3, angle1, angle2, 
     &  angle3 )

c*********************************************************************72
c
cc CIRCLE_TRIPLE_ANGLE_2D returns an angle formed by three circles in 2D.
c
c  Discussion:
c
c    A circle triple is a set of three tangent circles.  We assume
c    that no circle is contained in another.
c
c    We consider the triangle formed by joining the centers of the circles.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kenneth Stephenson,
c    Circle Packing, The Theory of Discrete Analytic Functions,
c    Cambridge, 2005.
c
c  Parameters:
c
c    Input, double precision R1, R2, R3, the radii of the circles.
c
c    Input, double precision ANGLE1, ANGLE2, ANGLE3, the angles
c    in the triangle.
c 
      implicit none

      double precision angle1
      double precision angle2
      double precision angle3
      double precision r1
      double precision r2
      double precision r3
      double precision r8_acos

      angle1 = r8_acos (
     &  ( r1 + r2 )**2 + ( r1 + r3 )**2 - ( r2 + r3 )**2 ) / 
     &  ( 2.0D+00 * ( r1 + r2 ) * ( r1 + r3 ) ) 

      angle2 = r8_acos ( 
     &  ( r2 + r3 )**2 + ( r2 + r1 )**2 - ( r3 + r1 )**2 ) / 
     &  ( 2.0D+00 * ( r2 + r3 ) * ( r2 + r1 ) ) 

      angle3 = r8_acos ( 
     &  ( r3 + r1 )**2 + ( r3 + r2 )**2 - ( r1 + r2 )**2 ) / 
     &  ( 2.0D+00 * ( r3 + r1 ) * ( r3 + r2 ) ) 

      return
      end
      subroutine circles_imp_int_2d ( r1, pc1, r2, pc2, int_num, p )

c*********************************************************************72
c
cc CIRCLES_IMP_INT_2D: finds the intersection of two implicit circles in 2D.
c
c  Discussion:
c
c    Two circles can intersect in 0, 1, 2 or infinitely many points.
c
c    The 0 and 2 intersection cases are numerically robust; the 1 and
c    infinite intersection cases are numerically fragile.  The routine
c    uses a tolerance to try to detect the 1 and infinite cases.
c
c    Points P on an implicit circle in 2D satisfy the equation:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 = R**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, the radius of the first circle.
c
c    Input, double precision PC1(2), the center of the first circle.
c
c    Input, double precision R2, the radius of the second circle.
c
c    Input, double precision PC2(2), the center of the second circle.
c
c    Output, integer INT_NUM, the number of intersecting points 
c    found.  INT_NUM will be 0, 1, 2 or 3.  3 indicates that there are an 
c    infinite number of intersection points.
c
c    Output, double precision P(2,2), if INT_NUM is 1 or 2,
c    the coordinates of the intersecting points.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision distsq
      integer i
      integer int_num
      integer j
      double precision p(dim_num,2)
      double precision pc1(dim_num)
      double precision pc2(dim_num)
      double precision r1
      double precision r2
      double precision r8_epsilon
      double precision root
      double precision sc1
      double precision sc2
      double precision t1
      double precision t2
      double precision tol

      tol = r8_epsilon ( )

      do j = 1, 2
        do i = 1, dim_num
          p(i,j) = 0.0D+00
        end do
      end do
c
c  Take care of the case in which the circles have the same center.
c
      t1 = ( abs ( pc1(1) - pc2(1) ) 
     &     + abs ( pc1(2) - pc2(2) ) ) / 2.0D+00

      t2 = ( abs ( pc1(1) ) + abs ( pc2(1) ) 
     &     + abs ( pc1(2) ) + abs ( pc2(2) ) + 1.0D+00 ) / 5.0D+00

      if ( t1 .le. tol * t2 ) then

        t1 = abs ( r1 - r2 )
        t2 = ( abs ( r1 ) + abs ( r2 ) + 1.0D+00 ) / 3.0D+00

        if ( t1 .le. tol * t2 ) then
          int_num = 3
        else
          int_num = 0
        end if

        return

      end if

      distsq = ( pc1(1) - pc2(1) )**2 + ( pc1(2) - pc2(2) )**2

      root = 2.0D+00 * ( r1**2 + r2**2 ) * distsq - distsq**2 
     &  - ( r1 - r2 )**2 * ( r1 + r2 )**2

      if ( root .lt. -tol ) then
        int_num = 0
        return
      end if

      sc1 = ( distsq - ( r2**2 - r1**2 ) ) / distsq

      if ( root .lt. tol ) then
        int_num = 1
        do i = 1, dim_num
          p(i,1) = pc1(i) + 0.5D+00 * sc1 * ( pc2(i) - pc1(i) )
        end do
        return
      end if

      sc2 = sqrt ( root ) / distsq

      int_num = 2

      p(1,1) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) 
     &                - 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
      p(2,1) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) 
     &                + 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

      p(1,2) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) 
     &                + 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
      p(2,2) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) 
     &                - 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

      return
      end
      subroutine combin2 ( n, k, icnk )

c*********************************************************************72
c
cc COMBIN2 computes the binomial coefficient C(N,K).
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in integer arithmetic.
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
c    01 April 2011
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
c    April, 1963.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer ICNK, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer icnk
      integer k
      integer mn
      integer mx
      integer n

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        icnk = 0

      else if ( mn .eq. 0 ) then

        icnk = 1

      else

        mx = max ( k, n - k )
        icnk = mx + 1

        do i = 2, mn
          icnk = ( icnk * ( mx + i ) ) / i
        end do

      end if

      return
      end
      subroutine cone_area_3d ( h, r, area )

c*********************************************************************72
c
cc CONE_AREA_3D computes the surface area of a right circular cone in 3D.
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
c    Input, double precision H, R, the height of the cone, and the radius
c    of the circle that forms the base of the cone.
c
c    Output, double precision AREA, the surface area of the cone.
c
      implicit none

      double precision area
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      area = pi * r * sqrt ( h * h + r * r )

      return
      end
      subroutine cone_centroid_3d ( r, pc, pt, centroid )

c*********************************************************************72
c
cc CONE_CENTROID_3D returns the centroid of a cone in 3D.
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
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision R, the radius of the circle at the base of
c    the cone.
c
c    Input, double precision PC(3), the center of the circle.
c
c    Input, double precision PT(3), the coordinates of the tip of the cone.
c
c    Output, double precision CENTROID(3), the coordinates of the centroid
c    of the cone.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision centroid(dim_num)
      integer dim
      double precision pc(dim_num)
      double precision pt(dim_num)
      double precision r

      do dim = 1, dim_num
        centroid(dim) = 0.75D+00 * pc(dim) + 0.25D+00 * pt(dim)
      end do

      return
      end
      subroutine cone_volume_3d ( h, r, volume )

c*********************************************************************72
c
cc CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
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
c    Input, double precision H, R, the height of the cone, and the radius
c    of the circle that forms the base of the cone.
c
c    Output, double precision VOLUME, the volume of the cone.
c
      implicit none

      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume

      volume = pi * r * r * h / 3.0D+00

      return
      end
      subroutine conv3d ( axis, theta, n, cor3, cor2 )

c*********************************************************************72
c
cc CONV3D converts 3D data to a 2D projection.
c
c  Discussion:
c
c    A "presentation angle" THETA is used to project the 3D point
c    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
c
c    If AXIS = 'X':
c
c      X2D = Y3D - sin ( THETA ) * X3D
c      Y2D = Z3D - sin ( THETA ) * X3D
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character AXIS, the coordinate axis to be projected.
c    AXIS should be 'X', 'Y', or 'Z'.
c
c    Input, double precision THETA, the presentation angle in degrees.
c
c    Input, integer N, the number of points.
c
c    Input, double precision COR3(3,N), the 3D points.
c
c    Output, double precision COR2(2,N), the 2D projections.
c
      implicit none

      integer n

      character axis
      double precision cor2(2,n)
      double precision cor3(3,n)
      double precision degrees_to_radians
      integer j
      double precision stheta
      double precision theta

      stheta = sin ( degrees_to_radians ( theta ) )

      if ( axis .eq. 'X' .or. axis .eq. 'x' ) then

        do j = 1, n
          cor2(1,j) = cor3(2,j) - stheta * cor3(1,j)
          cor2(2,j) = cor3(3,j) - stheta * cor3(1,j)
        end do

      else if ( axis .eq. 'Y' .or. axis .eq. 'y' ) then

        do j = 1, n
          cor2(1,j) = cor3(1,j) - stheta * cor3(2,j)
          cor2(2,j) = cor3(3,j) - stheta * cor3(2,j)
        end do

      else if ( axis .eq. 'Z' .or. axis .eq. 'z' ) then

        do j = 1, n
          cor2(1,j) = cor3(1,j) - stheta * cor3(3,j)
          cor2(2,j) = cor3(2,j) - stheta * cor3(3,j)
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CONV3D - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Illegal coordinate index = "' // axis // '".'
        stop

      end if

      return
      end
      function cos_deg ( angle_deg )

c*********************************************************************72
c
cc COS_DEG returns the cosine of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_DEG, the angle, in degrees.
c
c    Output, double precision COS_DEG, the cosine of the angle.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      double precision cos_deg
      double precision degrees_to_radians

      degrees_to_radians = 3.141592653589793D+00 / 180.0D+00

      angle_rad = degrees_to_radians * angle_deg

      cos_deg  = cos ( angle_rad )

      return
      end
      function cot_deg ( angle_deg )

c*********************************************************************72
c
cc COT_DEG returns the cotangent of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_DEG, the angle, in degrees.
c
c    Output, double precision COT_DEG, the cotangent of the angle.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      double precision cot_deg
      double precision degrees_to_radians 
      parameter ( degrees_to_radians = 
     &  3.141592653589793D+00 / 180.0D+00 )

      angle_rad = degrees_to_radians * angle_deg

      cot_deg  = cos ( angle_rad ) / sin ( angle_rad )

      return
      end
      function cot_rad ( angle_rad )

c*********************************************************************72
c
cc COT_RAD returns the cotangent of an angle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_RAD, the angle, in radians.
c
c    Output, double precision COT_RAD, the cotangent of the angle.
c
      implicit none

      double precision angle_rad
      double precision cot_rad

      cot_rad  = cos ( angle_rad ) / sin ( angle_rad )

      return
      end
      function csc_deg ( angle_deg )

c*********************************************************************72
c
cc CSC_DEG returns the cosecant of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_DEG, the angle, in degrees.
c
c    Output, double precision CSC_DEG, the cosecant of the angle.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      double precision csc_deg
      double precision degrees_to_radians 
      parameter ( degrees_to_radians = 
     &  3.141592653589793D+00 / 180.0D+00 )

      angle_rad = degrees_to_radians * angle_deg
      csc_deg  = 1.0D+00 / sin ( angle_rad )

      return
      end
      subroutine cube_shape_3d ( point_num, face_num, face_order_max, 
     &  point_coord, face_order, face_point )

c*********************************************************************72
c
cc CUBE_SHAPE_3D describes a cube in 3D.
c
c  Discussion:
c
c    The vertices lie on the unit sphere.
c
c    The dual of the cube is the octahedron.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer FACE_NUM, the number of faces.
c
c    Input, integer FACE_ORDER_MAX, the maximum number of vertices
c    in a face.
c
c    Output, double precision POINT_COORD(3,POINT_NUM),
c    the vertices.
c
c    Output, integer FACE_ORDER(FACE_NUM), the number of vertices
c    per face.
c
c    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
c    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  The
c    points are listed in the counter clockwise direction defined
c    by the outward normal at the face.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer face_num
      integer face_num_save
      parameter ( face_num_save = 6 )
      integer face_order_max
      integer face_order_max_save
      parameter ( face_order_max_save = 4 )
      integer point_num
      integer point_num_save
      parameter ( point_num_save = 8 )

      double precision a
      integer face_order(face_num)
      integer face_order_save(face_num_save)
      integer face_point(face_order_max,face_num)
      integer face_point_save(face_order_max_save,face_num_save)
      integer i
      integer j
      double precision point_coord(dim_num,point_num)
      double precision point_coord_save(dim_num,point_num_save)

      save face_order_save
      save face_point_save
      save point_coord_save

      data face_order_save / 4, 4, 4, 4, 4, 4 /

      data face_point_save /
     &   1, 4, 3, 2, 
     &   1, 2, 6, 5, 
     &   2, 3, 7, 6, 
     &   3, 4, 8, 7, 
     &   1, 5, 8, 4, 
     &   5, 6, 7, 8 /

      data point_coord_save /
     &   -0.577350269189626D+00, -0.577350269189626D+00,
     &   -0.577350269189626D+00, 
     &    0.577350269189626D+00, -0.577350269189626D+00,
     &   -0.577350269189626D+00, 
     &    0.577350269189626D+00,  0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &   -0.577350269189626D+00,  0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &   -0.577350269189626D+00, -0.577350269189626D+00,
     &    0.577350269189626D+00, 
     &    0.577350269189626D+00, -0.577350269189626D+00,   
     &    0.577350269189626D+00, 
     &    0.577350269189626D+00,  0.577350269189626D+00, 
     &    0.577350269189626D+00, 
     &   -0.577350269189626D+00,  0.577350269189626D+00, 
     &    0.577350269189626D+00 /

      do i = 1, face_num
        face_order(i) = face_order_save(i)
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
      subroutine cube_size_3d ( point_num, edge_num, face_num,
     &  face_order_max )

c*********************************************************************72
c
cc CUBE_SIZE_3D gives "sizes" for a cube in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 December 2008
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

      point_num = 8
      edge_num = 12
      face_num = 6
      face_order_max = 4

      return
      end
      subroutine cylinder_point_dist_3d ( p1, p2, r, p, distance )

c*********************************************************************72
c
cc CYLINDER_POINT_DIST_3D: distance from a cylinder to a point in 3D.
c
c  Discussion:
c
c    We are computing the distance to the SURFACE of the cylinder.
c
c    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
c    which is the line segment from point P1 to P2, and a radius R.  The points 
c    on the surface of the cylinder are:
c    * points at a distance R from the line through P1 and P2, and whose nearest
c      point on the line through P1 and P2 is strictly between P1 and P2, 
c    PLUS
c    * points at a distance less than or equal to R from the line through P1
c      and P2, whose nearest point on the line through P1 and P2 is either 
c      P1 or P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), the first and last points
c    on the axis line of the cylinder.
c
c    Input, double precision R, the radius of the cylinder.
c
c    Input, double precision P(3), the point.
c
c    Output, double precision DISTANCE, the distance from the point 
c    to the cylinder.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision axis(dim_num)
      double precision axis_length
      double precision distance
      integer i
      double precision off_axis_component
      double precision p(dim_num)
      double precision p_dot_axis
      double precision p_length
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r
      double precision r8_huge
      double precision r8vec_diff_norm
      double precision r8vec_norm

      do i = 1, dim_num
        axis(i) = p2(i) - p1(i)
      end do

      axis_length = r8vec_norm ( dim_num, axis )

      if ( axis_length .eq. 0.0D+00 ) then
        distance = - r8_huge ( )
        return
      end if

      do i = 1, dim_num
        axis(i) = axis(i) / axis_length
      end do

      p_dot_axis = 0.0D+00
      do i = 1, dim_num
        p_dot_axis = p_dot_axis + ( p(i) - p1(i) ) * axis(i)
      end do
c
c  Case 1: Below bottom cap.
c
      if ( p_dot_axis .le. 0.0D+00 ) then

        call disk_point_dist_3d ( p1, r, axis, p, distance )
c
c  Case 2: between cylinder planes.
c
      else if ( p_dot_axis .le. axis_length ) then

        p_length = r8vec_diff_norm ( dim_num, p, p1 )
        off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

        distance = abs ( off_axis_component - r )

        if ( off_axis_component .lt. r ) then
          distance = min ( distance, axis_length - p_dot_axis )
          distance = min ( distance, p_dot_axis )
        end if
c
c  Case 3: Above the top cap.
c  
      else if ( axis_length .lt. p_dot_axis ) then

        call disk_point_dist_3d ( p2, r, axis, p, distance )

      end if

      return
      end
      subroutine cylinder_point_inside_3d ( p1, p2, r, p, inside )

c*********************************************************************72
c
cc CYLINDER_POINT_INSIDE_3D determines if a cylinder contains a point in 3D.
c
c  Discussion:
c
c    The surface and interior of a (right) (finite) cylinder in 3D is defined 
c    by an axis, which is the line segment from point P1 to P2, and a 
c    radius R.  The points contained in the volume include:
c    * points at a distance less than or equal to R from the line through P1
c      and P2, whose nearest point on the line through P1 and P2 is, in fact,
c      P1, P2, or any point between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), the first and last points
c    on the axis line of the cylinder.
c
c    Input, double precision R, the radius of the cylinder.
c
c    Input, double precision P(3), the point.
c
c    Output, logical INSIDE, is TRUE if the point is inside the cylinder.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
     
      double precision axis(dim_num)
      double precision axis_length
      integer i
      logical inside
      double precision off_axis_component
      double precision p(dim_num)
      double precision p_dot_axis
      double precision p_length
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r
      double precision r8vec_diff_norm
      double precision r8vec_norm

      do i = 1, dim_num
        axis(i) = p2(i) - p1(i)
      end do

      axis_length = r8vec_norm ( dim_num, axis )

      if ( axis_length .eq. 0.0D+00 ) then
        inside = .false.
        return
      end if

      do i = 1, dim_num
        axis(i) = axis(i) / axis_length
      end do

      p_dot_axis = 0.0D+00
      do i = 1, dim_num
        p_dot_axis = p_dot_axis +  ( p(i) - p1(i) ) * axis(i)
      end do
c
c  If the point lies below or above the "caps" of the cylinder, we're done.
c
      if ( p_dot_axis .lt. 0.0D+00 .or. 
     &     axis_length .lt. p_dot_axis ) then

        inside = .false.
c
c  Otherwise, determine the distance from P to the axis.
c
      else

        p_length = r8vec_diff_norm ( dim_num, p, p1 )

        off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

        if ( off_axis_component .le. r ) then
          inside = .true.
        else
          inside = .false.
        end if

      end if

      return
      end
      subroutine cylinder_point_near_3d ( p1, p2, r, p, pn )

c*********************************************************************72
c
cc CYLINDER_POINT_NEAR_3D: nearest point on a cylinder to a point in 3D.
c
c  Discussion:
c
c    We are computing the nearest point on the SURFACE of the cylinder.
c
c    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
c    which is the line segment from point P1 to P2, and a radius R.  The points 
c    on the surface of the cylinder are:
c    * points at a distance R from the line through P1 and P2, and whose nearest
c      point on the line through P1 and P2 is strictly between P1 and P2, 
c    PLUS
c    * points at a distance less than or equal to R from the line through P1
c      and P2, whose nearest point on the line through P1 and P2 is either 
c      P1 or P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), the first and last points
c    on the axis line of the cylinder.
c
c    Input, double precision R, the radius of the cylinder.
c
c    Input, double precision P(3), the point.
c
c    Output, double precision PN(3), the nearest point on the cylinder.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision axial_component
      double precision axis(dim_num)
      double precision axis_length
      double precision distance
      integer i
      double precision off_axis(dim_num)
      double precision off_axis_component
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision r
      double precision r8vec_dot_product
      double precision r8vec_norm

      do i = 1, dim_num
        axis(i) = p2(i) - p1(i)
      end do

      axis_length = r8vec_norm ( dim_num, axis )

      do i = 1, dim_num
        axis(i) = axis(i) / axis_length
      end do

      axial_component = 0.0D+00
      do i = 1, dim_num
        axial_component = axial_component + ( p(i) - p1(i) ) * axis(i)
      end do

      do i = 1, dim_num
        off_axis(i) = p(i) - p1(i) - axial_component * axis(i)
      end do

      off_axis_component = r8vec_norm ( dim_num, off_axis )
c
c  Case 1: Below bottom cap.
c
      if ( axial_component .le. 0.0D+00 ) then

        if ( off_axis_component .le. r ) then
          do i = 1, dim_num
            pn(i) = p1(i) + off_axis(i)
          end do
        else
          do i = 1, dim_num
            pn(i) = p1(i) + ( r / off_axis_component ) * off_axis(i)
          end do
        end if
c
c  Case 2: between cylinder planes.
c
      else if ( axial_component .le. axis_length ) then

        if ( off_axis_component .eq. 0.0D+00 ) then

          call r8vec_any_normal ( dim_num, axis, off_axis )
          
          do i = 1, dim_num
            pn(i) = p(i) + r * off_axis(i)
          end do

        else

          distance = abs ( off_axis_component - r )

          do i = 1, dim_num
            pn(i) = p1(i) + axial_component * axis(i) 
     &        + ( r / off_axis_component ) * off_axis(i)
          end do

          if ( off_axis_component .lt. r ) then

            if ( axis_length - axial_component .lt. distance ) then
              distance = axis_length - axial_component
              do i = 1, dim_num
                pn(i) = p2(i) + off_axis(i)
              end do
            end if

            if ( axial_component .lt. distance ) then
              distance = axial_component
              do i = 1, dim_num
                pn(i) = p1(i) + off_axis(i)
              end do
            end if

          end if

        end if
c
c  Case 3: Above the top cap.
c  
      else if ( axis_length .lt. axial_component ) then

        if ( off_axis_component .le. r ) then
          do i = 1, dim_num
            pn(i) = p2(i) + off_axis(i)
          end do
        else
          do i = 1, dim_num
            pn(i) = p2(i) + ( r / off_axis_component ) * off_axis(i)
          end do
        end if

      end if

      return
      end
      subroutine cylinder_sample_3d ( p1, p2, r, n, seed, p )

c*********************************************************************72
c
cc CYLINDER_SAMPLE_3D samples a cylinder in 3D.
c
c  Discussion:
c
c    We are sampling the interior of a right finite cylinder in 3D.
c
c    The interior of a (right) (finite) cylinder in 3D is defined by an axis,
c    which is the line segment from point P1 to P2, and a radius R.  The points 
c    on or inside the cylinder are:
c    * points whose distance from the line through P1 and P2 is less than
c      or equal to R, and whose nearest point on the line through P1 and P2
c      lies (nonstrictly) between P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), the first and last points
c    on the axis line of the cylinder.
c
c    Input, double precision R, the radius of the cylinder.
c
c    Input, integer N, the number of sample points to compute.
c
c    Input/output, integer SEED, the random number seed.
c
c    Input, double precision P(3,N), the sample points.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer n

      double precision axis(dim_num)
      double precision axis_length
      integer i
      integer j
      double precision p(dim_num,n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r
      double precision r8vec_norm
      double precision radius(n)
      integer seed
      double precision theta(n)
      double precision v2(dim_num)
      double precision v3(dim_num)
      double precision z(n)
c
c  Compute the axis vector.
c
      do i = 1, dim_num
        axis(i) = p2(i) - p1(i)
      end do

      axis_length = r8vec_norm ( dim_num, axis )

      do i = 1, dim_num
        axis(i) = axis(i) / axis_length
      end do
c
c  Compute vectors V2 and V3 that form an orthogonal triple with AXIS.
c
      call plane_normal_basis_3d ( p1, axis, v2, v3 )
c
c  Assemble the randomized information.
c
      call r8vec_uniform_01 ( n, seed, radius )

      do i = 1, n
        radius(i) = r * sqrt ( radius(i) )
      end do

      call r8vec_uniform_01 ( n, seed, theta )

      do i = 1, n
        theta(i) = 2.0D+00 * pi * theta(i)
      end do

      call r8vec_uniform_01 ( n, seed, z )

      do i = 1, n
        z(i) = axis_length * z(i)
      end do

      do j = 1, n
        do i = 1, dim_num

          p(i,j) =                                     p1(i)   
     &              + z(j)                           * axis(i) 
     &              + radius(j) * cos ( theta(j) )   * v2(i)   
     &              + radius(j) * sin ( theta(j) )   * v3(i)
        end do
      end do


      return
      end
      subroutine cylinder_volume_3d ( p1, p2, r, volume )

c*********************************************************************72
c
cc CYLINDER_VOLUME_3D determines the volume of a cylinder in 3D.
c
c  Discussion:
c
c    The surface and interior of a (right) (finite) cylinder in 3D is defined by
c    an axis, which is the line segment from point P1 to P2, and a radius R.
c    The points contained in the volume include:
c    * points at a distance less than or equal to R from the line through P1
c      and P2, whose nearest point on the line through P1 and P2 is, in fact,
c      P1, P2, or any point between them.
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
c    Input, double precision P1(3), P2(3), the first and last points
c    on the axis line of the cylinder.
c
c    Input, double precision R, the radius of the cylinder.
c
c    Output, double precision VOLUME, the volume of the cylinder.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer dim
      double precision h
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume

      h = 0.0D+00
      do dim = 1, dim_num
        h = h + ( p2(dim) - p1(dim) )**2
      end do
      h = sqrt ( h )

      volume = pi * r * r * h

      return
      end
      function degrees_to_radians ( angle_deg )

c*********************************************************************72
c
cc DEGREES_TO_RADIANS converts an angle from degrees to radians.
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
c    Input, double precision ANGLE_DEG, an angle in degrees.
c
c    Output, double precision DEGREES_TO_RADIANS, the equivalent angle
c    in radians.
c
      implicit none

      double precision angle_deg
      double precision degrees_to_radians
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      degrees_to_radians = ( angle_deg / 180.0D+00 ) * pi

      return
      end
      subroutine direction_pert_3d ( sigma, vbase, seed, vran )

c*********************************************************************72
c
cc DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision SIGMA, determines the strength of the
c    perturbation.
c    SIGMA <= 0 results in a completely random direction.
c    1 <= SIGMA results in VBASE.
c    0 < SIGMA < 1 results in a perturbation from VBASE, which is
c    large when SIGMA is near 0, and small when SIGMA is near 1.
c
c    Input, double precision VBASE(3), the base direction vector, which
c    should have unit norm.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision VRAN(3), the perturbed vector, which will
c    have unit norm.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dphi
      integer i
      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision psi
      double precision r
      double precision r8_acos
      double precision r8_uniform_01
      integer seed
      double precision sigma
      double precision theta
      double precision v(dim_num)
      double precision vbase(dim_num)
      double precision vdot
      double precision vran(dim_num)
      double precision x
c
c  1 <= SIGMA, just use the base vector.
c
      if ( 1.0D+00 .le. sigma ) then

        do i = 1, dim_num
          vran(i) = vbase(i)
        end do

      else if ( sigma .le. 0.0D+00 ) then

        vdot = r8_uniform_01 ( seed )
        vdot = 2.0D+00 * vdot - 1.0D+00

        phi = r8_acos ( vdot )

        theta = r8_uniform_01 ( seed )
        theta = 2.0D+00 * pi * theta

        vran(1) = cos ( theta ) * sin ( phi )
        vran(2) = sin ( theta ) * sin ( phi )
        vran(3) = cos ( phi )

      else

        phi = r8_acos ( vbase(3) )
        theta = atan2 ( vbase(2), vbase(1) )
c
c  Pick VDOT, which must be between -1 and 1.  This represents
c  the dot product of the perturbed vector with the base vector.
c
c  R8_UNIFORM_01 returns a uniformly random value between 0 and 1.
c  The operations we perform on this quantity tend to bias it
c  out towards 1, as SIGMA grows from 0 to 1.
c
c  VDOT, in turn, is a value between -1 and 1, which, for large
c  SIGMA, we want biased towards 1.
c
        r = r8_uniform_01 ( seed )
        x = exp ( ( 1.0D+00 - sigma ) * log ( r ) )
        dphi = r8_acos ( 2.0D+00 * x - 1.0D+00 )
c
c  Now we know enough to write down a vector that is rotated DPHI
c  from the base vector.
c
        v(1) = cos ( theta ) * sin ( phi + dphi )
        v(2) = sin ( theta ) * sin ( phi + dphi )
        v(3) = cos ( phi + dphi )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the base vector.
c
        psi = r8_uniform_01 ( seed )
        psi = 2.0D+00 * pi * psi
c
c  Carry out the rotation.
c
        call rotation_axis_vector_3d ( vbase, psi, v, vran )

      end if

      return
      end
      subroutine direction_uniform_2d ( seed, vran )

c*********************************************************************72
c
cc DIRECTION_UNIFORM_2D picks a random direction vector in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
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
c    Output, double precision VRAN(2), the random direction vector, with
c    unit norm.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision r8_uniform_01
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seed
      double precision theta
      double precision vran(dim_num)

      theta = r8_uniform_01 ( seed )
      theta = 2.0D+00 * pi * theta

      vran(1) = cos ( theta )
      vran(2) = sin ( theta )

      return
      end
      subroutine direction_uniform_3d ( seed, vran )

c*********************************************************************72
c
cc DIRECTION_UNIFORM_3D picks a random direction vector in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
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
c    Output, double precision VRAN(3), the random direction vector,
c    with unit norm.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_acos
      double precision r8_uniform_01
      integer seed
      double precision theta
      double precision vdot
      double precision vran(dim_num)
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

      phi = r8_acos ( vdot )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the Z vector.
c
      theta = r8_uniform_01 ( seed )
      theta = 2.0D+00 * pi * theta

      vran(1) = cos ( theta ) * sin ( phi )
      vran(2) = sin ( theta ) * sin ( phi )
      vran(3) = cos ( phi )

      return
      end
      subroutine direction_uniform_nd ( dim_num, seed, w )

c*********************************************************************72
c
cc DIRECTION_UNIFORM_ND generates a random direction vector in ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision W(DIM_NUM), a random direction vector,
c    with unit norm.
c
      implicit none

      integer dim_num

      integer i
      double precision norm
      double precision r8vec_norm
      integer seed
      double precision w(dim_num)
c
c  Get N values from a standard normal distribution.
c
      call r8vec_normal_01 ( dim_num, seed, w )
c
c  Compute the length of the vector.
c
      norm = r8vec_norm ( dim_num, w )
c
c  Normalize the vector.
c
      do i = 1, dim_num
        w(i) = w(i) / norm
      end do

      return
      end
      subroutine disk_point_dist_3d ( pc, r, axis, p, dist )

c*********************************************************************72
c
cc DISK_POINT_DIST_3D determines the distance from a disk to a point in 3D.
c
c  Discussion:
c
c    A disk in 3D satisfies the equations:
c
c      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 + ( P(3) - PC(3) <= R**2
c
c    and
c
c      P(1) * AXIS(1) + P(2) * AXIS(2) + P(3) * AXIS(3) = 0
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
c    Input, double precision PC(3), the center of the disk.
c
c    Input, double precision R, the radius of the disk.
c
c    Input, double precision AXIS(3), the axis vector.
c
c    Input, double precision P(3), the point to be checked.
c
c    Output, double precision DIST, the distance of the point to the disk.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision axial_component
      double precision axis(dim_num)
      double precision axis_length
      double precision dist
      integer i
      double precision off_axis_component
      double precision off_axis(dim_num)
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision r
      double precision r8_huge
      double precision r8vec_diff_norm
      double precision r8vec_dot_product
      double precision r8vec_norm
c
c  Special case: the point is the center.
c
      do i = 1, dim_num
        if ( p(i) .eq. pc(i) ) then
          dist = 0.0D+00
          return
        end if
      end do

      axis_length = r8vec_norm ( dim_num, axis )

      if ( axis_length .eq. 0.0D+00 ) then
        dist = - r8_huge ( )
        return
      end if

      axial_component = 0.0D+00
      do i = 1, dim_num
        axial_component = axial_component + ( p(i) - pc(i) ) * axis(i)
      end do
      axial_component = axial_component / axis_length
c
c  Special case: the point satisfies the disk equation exactly.
c
      if ( r8vec_diff_norm ( dim_num, p, pc ) <= r .and.
     &      axial_component .eq. 0.0D+00 ) then
        dist = 0.0D+00
        return
      end if
c
c  Decompose P-PC into axis component and off-axis component.
c
      do i = 1, dim_num
        off_axis(i) = p(i) - pc(i) 
     &    - axial_component * axis(i) / axis_length
      end do

      off_axis_component = r8vec_norm ( dim_num, off_axis )
c
c  If the off-axis component has norm less than R, the nearest point is
c  the projection to the disk along the axial direction, and the distance 
c  is just the dot product of P-PC with unit AXIS.
c
      if ( off_axis_component .le. r ) then
        dist = abs ( axial_component )
        return
      end if
c
c  Otherwise, the nearest point is along the perimeter of the disk.
c
      dist = sqrt ( axial_component**2 + ( off_axis_component - r )**2 )

      return
      end
      subroutine dms_to_radians ( degrees, minutes, seconds, radians )

c*********************************************************************72
c
cc DMS_TO_RADIANS converts an angle from degrees/minutes/seconds to radians.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DEGREES, MINUTES, SECONDS, an angle in 
c    degrees, minutes, and seconds.
c
c    Output, double precision RADIANS, the equivalent angle in radians.
c
      implicit none

      double precision angle
      integer degrees
      integer minutes
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision radians
      integer seconds

      angle =   dble ( degrees ) 
     &      + ( dble ( minutes ) 
     &      + ( dble ( seconds ) / 60.0D+00 ) ) / 60.0D+00

      radians = ( angle / 180.0D+00 ) * pi

      return
      end
      subroutine dodec_shape_3d ( point_num, face_num, face_order_max, 
     &  point_coord, face_order, face_point )

c*********************************************************************72
c
cc DODEC_SHAPE_3D describes a dodecahedron in 3D.
c
c  Discussion:
c
c    The vertices lie on the unit sphere.
c
c    The dual of a dodecahedron is an icosahedron.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer FACE_NUM, the number of faces.
c
c    Input, integer FACE_ORDER_MAX, the maximum number of vertices
c    per face.
c
c    Output, double precision POINT_COORD(3,POINT_NUM), the vertices.
c
c    Output, integer FACE_ORDER[FACE_NUM], the number of vertices
c    per face.
c
c    Output, integer FACE_POINT(FACE_ORDER_MAX,POINT_NUM); 
c    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  
c    The points are listed in the counter clockwise direction defined
c    by the outward normal at the face.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer face_num
      integer face_num_save
      parameter ( face_num_save = 12 )
      integer face_order_max
      integer face_order_max_save
      parameter ( face_order_max_save = 5 )
      integer point_num
      integer point_num_save
      parameter ( point_num_save = 20 )

      integer i
      integer j
      integer face_order(face_num)
      integer face_order_save(face_num_save)
      integer face_point(face_order_max,face_num)
      integer face_point_save(face_order_max_save,face_num_save)
      double precision point_coord(dim_num,point_num)
      double precision point_coord_save(dim_num,point_num_save)

      save face_order_save
      save face_point_save
      save point_coord_save

      data face_order_save /
     &  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 /
      data face_point_save /
     &    2,  9,  1, 13, 14, 
     &    5, 10,  6, 16, 15, 
     &    3, 11,  4, 14, 13, 
     &    8, 12,  7, 15, 16, 
     &    3, 13,  1, 17, 18, 
     &    2, 14,  4, 20, 19, 
     &    5, 15,  7, 18, 17, 
     &    8, 16,  6, 19, 20, 
     &    5, 17,  1,  9, 10, 
     &    3, 18,  7, 12, 11, 
     &    2, 19,  6, 10,  9, 
     &    8, 20,  4, 11, 12 /
      data point_coord_save /
     &    0.577350269189626D+00,  0.577350269189626D+00, 
     &    0.577350269189626D+00, 
     &    0.577350269189626D+00,  0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &    0.577350269189626D+00, -0.577350269189626D+00, 
     &    0.577350269189626D+00, 
     &    0.577350269189626D+00, -0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &   -0.577350269189626D+00,  0.577350269189626D+00, 
     &    0.577350269189626D+00, 
     &   -0.577350269189626D+00,  0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &   -0.577350269189626D+00, -0.577350269189626D+00, 
     &    0.577350269189626D+00, 
     &   -0.577350269189626D+00, -0.577350269189626D+00, 
     &   -0.577350269189626D+00, 
     &    0.356822089773090D+00,  0.934172358962716D+00, 
     &    0.000000000000000D+00, 
     &   -0.356822089773090D+00,  0.934172358962716D+00, 
     &    0.000000000000000D+00, 
     &    0.356822089773090D+00, -0.934172358962716D+00, 
     &    0.000000000000000D+00, 
     &   -0.356822089773090D+00, -0.934172358962716D+00, 
     &    0.000000000000000D+00, 
     &    0.934172358962716D+00,  0.000000000000000D+00, 
     &    0.356822089773090D+00, 
     &    0.934172358962716D+00,  0.000000000000000D+00, 
     &   -0.356822089773090D+00, 
     &   -0.934172358962716D+00,  0.000000000000000D+00, 
     &    0.356822089773090D+00, 
     &   -0.934172358962716D+00,  0.000000000000000D+00, 
     &   -0.356822089773090D+00, 
     &    0.000000000000000D+00,  0.356822089773090D+00, 
     &    0.934172358962716D+00, 
     &    0.000000000000000D+00, -0.356822089773090D+00, 
     &    0.934172358962716D+00, 
     &    0.000000000000000D+00,  0.356822089773090D+00, 
     &   -0.934172358962716D+00, 
     &    0.000000000000000D+00, -0.356822089773090D+00,
     &   -0.934172358962716D+00 /

      do i = 1, face_num
        face_order(i) = face_order_save(i)
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
      subroutine dodec_size_3d ( point_num, edge_num, face_num,
     &  face_order_max )

c*********************************************************************72
c
cc DODEC_SIZE_3D gives "sizes" for a dodecahedron in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2007
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

      point_num = 20
      edge_num = 30
      face_num = 12
      face_order_max = 5

      return
      end
      subroutine dual_shape_3d ( point_num, face_num, face_order_max, 
     &  point_coord, face_order, face_point, point_num2, face_num2, 
     &  face_order_max2, point_coord2, face_order2, face_point2 )

c*********************************************************************72
c
cc DUAL_SHAPE_3D constructs the dual of a shape in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer FACE_NUM, the number of faces.
c
c    Input, integer FACE_ORDER_MAX, the maximum number of vertices
c    per face.
c
c    Input, double precision POINT_COORD(3,POINT_NUM), the points.
c
c    Input, integer FACE_ORDER(FACE_NUM), the number of vertices
c    per face.
c
c    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
c    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
c    points are listed in the counter clockwise direction defined
c    by the outward normal at the face.
c
c    Input, integer POINT_NUM2, the number of points in the dual.
c
c    Input, integer FACE_NUM2, the number of faces in the dual.
c
c    Input, integer FACE_ORDER_MAX2, the maximum number of 
c    vertices per face in the dual.
c
c    Output, double precision POINT_COORD2(3,POINT_NUM2), the point 
c    coordinates of the dual.
c
c    Output, integer FACE_ORDER2(FACE_NUM2), the number of 
c    vertices per face.
c
c    Output, integer FACE_POINT2(FACE_ORDER_MAX2,FACE_NUM2), 
c    the vertices of each face in the dual.
c
      implicit none

      integer face_num
      integer face_num2
      integer face_order_max
      integer face_order_max2
      integer dim_num
      parameter ( dim_num = 3 )
      integer point_num
      integer point_num2

      integer col
      integer face
      integer face_order(face_num)
      integer face_order2(face_num2)
      integer face_point(face_order_max,face_num)
      integer face_point2(face_order_max2,face_num2)
      integer i
      integer inext
      integer iprev
      integer istop
      integer j
      integer k
      double precision norm
      double precision p(dim_num)
      double precision point_coord(dim_num,point_num)
      double precision point_coord2(dim_num,point_num2)
      double precision r8vec_norm
      integer row
c
c  This computation should really compute the center of gravity
c  of the face, in the general case.
c
c  We'll also assume the vertices of the original and the dual
c  are to lie on the unit sphere, so we can normalize the
c  position vector of the vertex.
c
      do face = 1, face_num

        do i = 1, dim_num
          p(i) = 0.0D+00
        end do

        do j = 1, face_order(face)
          k = face_point(j,face)
          do i = 1, dim_num
            p(i) = p(i) + point_coord(i,k)
          end do
        end do

        norm = r8vec_norm ( dim_num, p )

        do i = 1, dim_num
          point_coord2(i,face) = p(i) / norm
        end do

      end do
c
c  Now build the face in the dual associated with each node FACE.
c
      do face = 1, face_num2
c
c  Initialize the order.
c
        face_order2(face) = 0
c
c  Find the first occurrence of FACE in an edge of polyhedron.
c
        call i4col_find_item ( face_order_max, face_num, face_point, 
     &    face, row, col )

        if ( row .le. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
          write ( *, '(a,i8)' ) 
     &      '  Could not find an edge using node ', face
          stop
        end if
c
c  Save the following node as ISTOP.
c  When we encounter ISTOP again, this will mark the end of our search.
c
        i = row + 1
        if ( face_order(col) .lt. i ) then
          i = 1
        end if

        istop = face_point(i,col)
c
c  Save the previous node as INEXT.
c
10      continue

          i = row - 1
          if ( i .lt. 1 ) then
            i = i + face_order(col)
          end if

          inext = face_point(i,col)

          face_order2(face) = face_order2(face) + 1

          face_point2(face_order2(face),face) = col
c
c  If INEXT =/= ISTOP, continue.
c
          if ( inext .eq. istop ) then
            go to 20
          end if
c
c  Set IPREV:= INEXT.
c
          iprev = inext
c
c  Search for the occurrence of the edge FACE-IPREV.
c
          call i4col_find_pair_wrap ( face_order_max, face_num, 
     &      face_point, face, iprev, row, col )

          if ( row .le. 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
            write ( *, '(a,i8)' ) '  No edge from node ', iprev
            write ( *, '(a,i8)' ) '  to node ', face
            stop
          end if

        go to 10

20      continue

      end do

      return
      end
      subroutine dual_size_3d ( point_num, edge_num, face_num,
     &  face_order_max, point_coord, face_order, face_point,
     &   point_num2, edge_num2, face_num2, face_order_max2 )

c*********************************************************************72
c
cc DUAL_SIZE_3D determines sizes for a dual of a shape in 3D.
c
c  Discussion:
c
c    We don't actually need FACE_POINT as input here.  But since the
c    three arrays occur together everywhere else, it seems unnecessarily
c    user-confusing to vary the usage herec
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
c    Input, integer POINT_NUM, the number of points.
c
c    Input, integer EDGE_NUM, the number of edges.
c
c    Input, integer FACE_NUM, the number of faces.
c
c    Input, integer  FACE_ORDER_MAX, the maximum number of vertices
c    per face.
c
c    Input, double precision POINT_COORD(3,POINT_NUM), the points.
c
c    Input, integer FACE_ORDER(FACE_NUM), the number of vertices
c    per face.
c
c    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM);
c    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
c    points are listed in the counter clockwise direction defined
c    by the outward normal at the face.
c
c    Output, integer POINT_NUM2, the number of points in the dual.
c
c    Output, integer EDGE_NUM2, the number of edges in the dual.
c
c    Output, integer FACE_NUM2, the number of faces in the dual.
c
c    Output, integer FACE_ORDER_MAX2, the maximum number of
c    vertices per face in the dual.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer face_num
      integer face_order_max
      integer point_num

      integer edge_num
      integer edge_num2
      integer face
      integer face_num2
      integer face_order(face_num)
      integer face_order2(point_num)
      integer face_order_max2
      integer face_point(face_order_max,face_num)
      integer face2
      integer i
      integer point_num2
      double precision point_coord(dim_num,point_num)
c
c  These values are easy to compute:
c
      point_num2 = face_num
      edge_num2 = edge_num
      face_num2 = point_num
c
c  To determine FACE_ORDER_MAX2 is not so easy.
c  You have to construct the FACE_ORDER array for the dual shape.
c  The order of a dual face is the number of edges that the vertex occurs in.
c  But then all we have to do is count how many times each item shows up
c  in the FACE_POINT array.
c
      do face2 = 1, face_num2
        face_order2(face2) = 0
      end do

      do face = 1, face_num
        do i = 1, face_order(face)
          face2 = face_point(i,face)
          face_order2(face2) = face_order2(face2) + 1
        end do
      end do

      face_order_max2 = 0
      do face2 = 1, face_num2
        face_order_max2 = max ( face_order_max2, face_order2(face2) )
      end do

      return
      end
      subroutine ellipse_area_2d ( r1, r2, area )

c*********************************************************************72
c
cc ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
c
c  Discussion:
c
c    An ellipse in standard position has a center at the origin, and
c    axes aligned with the coordinate axes.  Any point P on the ellipse
c    satisfies
c
c      (  P(1) / R1 )**2 + ( P(2) / R2 )**2 .eq. 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the "radius" of the ellipse in the major
c    and minor axis directions.  A circle has these values equal.
c
c    Output, double precision AREA, the area of the ellipse.
c
      implicit none

      double precision area
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2

      area = pi * r1 * r2

      return
      end
      subroutine ellipse_point_dist_2d ( r1, r2, p, dist )

c*********************************************************************72
c
cc ELLIPSE_POINT_DIST_2D finds the distance from a point to an ellipse in 2D.
c
c  Discussion:
c
c    An ellipse in standard position has a center at the origin, and
c    axes aligned with the coordinate axes.  Any point P on the ellipse
c    satisfies
c
c      (  P(1) / R1 )**2 + ( P(2) / R2 )**2 .eq. 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dianne O'Leary,
c    Elastoplastic Torsion: Twist and Stress,
c    Computing in Science and Engineering,
c    July/August 2004, pages 74-76.
c    September/October 2004, pages 63-65.
c
c  Parameters:
c
c    Input, double precision R1, R2, the ellipse parameters.  Normally,
c    these are both positive quantities.  Generally, they are also
c    distinct.
c
c    Input, double precision P(2), the point.
c
c    Output, double precision DIST, the distance to the ellipse.
c
      implicit none

      double precision dist
      double precision p(2)
      double precision pn(2)
      double precision r1
      double precision r2

      call ellipse_point_near_2d ( r1, r1, p, pn )

      dist = sqrt ( ( p(1) - pn(1) )**2 + ( p(2) - pn(2) )**2 )

      return
      end
      subroutine ellipse_point_near_2d ( r1, r2, p, pn )

c*********************************************************************72
c
cc ELLIPSE_POINT_NEAR_2D finds the nearest point on an ellipse in 2D.
c
c  Discussion:
c
c    An ellipse in standard position has a center at the origin, and
c    axes aligned with the coordinate axes.  Any point P on the ellipse
c    satisfies
c
c      (  P(1) / R1 )**2 + ( P(2) / R2 )**2 = 1
c
c    The nearest point PN on the ellipse has the property that the
c    line from PN to P is normal to the ellipse.  Points on the ellipse
c    can be parameterized by T, to have the form
c
c      ( R1 * cos ( T ), R2 * sin ( T ) ).
c
c    The tangent vector to the ellipse has the form
c
c      ( -R1 * sin ( T ), R2 * cos ( T ) ) 
c
c    At PN, the dot product of this vector with  ( P - PN ) must be
c    zero:
c
c      - R1 * sin ( T ) * ( X - R1 * cos ( T ) )
c      + R2 * cos ( T ) * ( Y - R2 * sin ( T ) ) = 0
c
c    This nonlinear equation for T can be solved by Newton's method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the ellipse parameters.  Normally,
c    these are both positive quantities.  Generally, they are also
c    distinct.
c
c    Input, double precision P(2), the point.
c
c    Output, double precision PN(2), the point on the ellipse which
c    is closest to P.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision ct
      double precision f
      double precision fp
      integer iteration
      integer iteration_max
      parameter ( iteration_max = 100 )
      double precision p(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision pn(dim_num)
      double precision r1
      double precision r2
      double precision r8_epsilon
      double precision st
      double precision t
      double precision x
      double precision y

      x = abs ( p(1) )
      y = abs ( p(2) )

      if ( y .eq. 0.0D+00 .and. r1 * r1 - r2 * r2 .le. r1 * x ) then

        t = 0.0D+00

      else if ( x .eq. 0.0D+00 .and. 
     &  r2 * r2 - r1 * r1 .le. r2 * y ) then

        t = pi / 2.0D+00

      else

        if ( y .eq. 0.0D+00 ) then
          y = sqrt ( r8_epsilon ( ) ) * abs ( r2 )
        end if

        if ( x .eq. 0.0D+00 ) then
          x = sqrt ( r8_epsilon ( ) ) * abs ( r1 )
        end if
c
c  Initial parameter T:
c
        t = atan2 ( y, x )

        iteration = 0

10      continue

          ct = cos ( t )
          st = sin ( t )

          f = ( x - abs ( r1 ) * ct ) * abs ( r1 ) * st 
     &      - ( y - abs ( r2 ) * st ) * abs ( r2 ) * ct

          if ( abs ( f ) .le. 100.0D+00 * epsilon ( f ) ) then
            go to 20
          end if

          if ( iteration_max .le. iteration ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ELLIPSE_POINT_NEAR_2D - Warning!'
            write ( *, '(a)' ) '  Reached iteration limit.'
            write ( *, '(a,f8.6)' ) '  T = ', t
            write ( *, '(a,g14.6)' ) '  F = ', f
            go to 20
          end if

          iteration = iteration + 1

          fp = r1 * r1 * st * st + r2 * r2 * ct * ct 
     &       + ( x - abs ( r1 ) * ct ) * abs ( r1 ) * ct 
     &       + ( y - abs ( r2 ) * st ) * abs ( r2 ) * st

          t = t - f / fp

        go to 10

20      continue

      end if
c
c  From the T value, we get the nearest point.
c
      pn(1) = abs ( r1 ) * cos ( t )
      pn(2) = abs ( r2 ) * sin ( t )
c
c  Take care of case where the point was in another quadrant.
c
      pn(1) = sign ( 1.0D+00, p(1) ) * pn(1)
      pn(2) = sign ( 1.0D+00, p(2) ) * pn(2)

      return
      end
      subroutine ellipse_points_2d ( pc, r1, r2, psi, n, p )

c*********************************************************************72
c
cc ELLIPSE_POINTS_2D returns N points on an tilted ellipse in 2D.
c
c  Discussion:
c
c    An ellipse in standard position has a center at the origin, and
c    axes aligned with the coordinate axes.  Any point P on the ellipse
c    satisfies
c
c      (  P(1) / R1 )**2 + ( P(2) / R2 )**2 .eq. 1
c
c    The points are "equally spaced" in the angular sense.  They are
c    not equally spaced along the perimeter of the ellipse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PC(2), the center of the ellipse.
c
c    Input, double precision R1, R2, the "radius" of the ellipse in the major
c    and minor axis directions.  A circle has these values equal.
c
c    Input, double precision PSI, the angle that the major axis of the ellipse
c    makes with the X axis.  A value of 0.0 means that the major and
c    minor axes of the ellipse will be the X and Y coordinate axes.
c
c    Input, integer N, the number of points desired.  N must 
c    be at least 1.
c
c    Output, double precision P(2,N), points on the ellipse.
c
      implicit none

      integer n

      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision psi
      double precision r1
      double precision r2
      double precision theta

      do i = 1, n

        theta = ( 2.0D+00 * pi * dble ( i - 1 ) ) / dble ( n )

        p(1,i) = pc(1) + r1 * cos ( psi ) * cos ( theta ) 
     &                 - r2 * sin ( psi ) * sin ( theta )

        p(2,i) = pc(2) + r1 * sin ( psi ) * cos ( theta ) 
     &                 + r2 * cos ( psi ) * sin ( theta )

      end do

      return
      end
      subroutine ellipse_points_arc_2d ( pc, r1, r2, psi, theta1, 
     &  theta2, n, p )

c*********************************************************************72
c
cc ELLIPSE_POINTS_ARC_2D returns N points on a tilted elliptical arc in 2D.
c
c  Discussion:
c
c    An ellipse in standard position has a center at the origin, and
c    axes aligned with the coordinate axes.  Any point P on the ellipse
c    satisfies
c
c      (  P(1) / R1 )**2 + ( P(2) / R2 )**2 == 1
c
c    The points are "equally spaced" in the angular sense.  They are
c    not equally spaced along the perimeter of the ellipse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PC(2), the coordinates of the center of
c    the ellipse.
c
c    Input, double precision R1, R2, the "radius" of the ellipse in the major
c    and minor axis directions.  A circle has these values equal.
c
c    Input, double precision PSI, the angle that the major axis of the ellipse
c    makes with the X axis.  A value of 0.0 means that the major and
c    minor axes of the ellipse will be the X and Y coordinate axes.
c
c    Input, double precision THETA1, THETA2, the angular coordinates of
c    the first and last points to be drawn, in radians.  This angle is measured
c    with respect to the (possibly tilted) major axis.
c
c    Input, integer N, the number of points desired.  N must 
c    be at least 1.
c
c    Output, double precision P(2,N), points on the ellipse.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      double precision r8_modp
      integer i
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision psi
      double precision r1
      double precision r2
      double precision theta
      double precision theta1
      double precision theta2
      double precision theta3
c
c  THETA3 is the smallest angle, no less than THETA1, which
c  coincides with THETA2.
c
      theta3 = theta1 + r8_modp ( theta2 - theta1, 2.0D+00 * pi )

      do i = 1, n

        if ( 1 .lt. n ) then
          theta = ( dble ( n - i     ) * theta1 
     &            + dble (     i - 1 ) * theta3 ) 
     &            / dble ( n     - 1 )
        else
          theta = 0.5D+00 * ( theta1 + theta3 )
        end if

        p(1,i) = pc(1) + r1 * cos ( psi ) * cos ( theta )
     &                 - r2 * sin ( psi ) * sin ( theta )

        p(2,i) = pc(2) + r1 * sin ( psi ) * cos ( theta )
     &                 + r2 * cos ( psi ) * sin ( theta )

      end do

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
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a valuer between 1 and 99, representing a
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

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine glob2loc_3d ( cospitch, cosroll, cosyaw, sinpitch, 
     &  sinroll, sinyaw, globas, glopts, locpts )

c*********************************************************************72
c
cc GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
c
c  Discussion:
c
c    A global coordinate system is given.
c
c    A local coordinate system has been translated to the point with
c    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
c    a roll.
c
c    A point has global coordinates GLOPTS, and it is desired to know
c    the point's local coordinates LOCPTS.
c
c    The transformation may be written as
c
c      LOC = M_ROLL * M_PITCH * M_YAW * ( GLOB - GLOBAS )
c
c    where
c
c               (       1            0            0      )
c    M_ROLL =   (       0        cos(Roll)    sin(Roll)  )
c               (       0      - sin(Roll)    cos(Roll)  )
c
c               (   cos(Pitch)       0      - sin(Pitch) )
c    M_PITCH =  (       0            1            0      )
c               (   sin(Pitch)       0        cos(Pitch) )
c
c               (   cos(Yaw)     sin(Yaw)         0      )
c    M_YAW    = ( - sin(Yaw)     cos(Yaw)         0      )
c               (       0            0            1      )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision COSPITCH, COSROLL, COSYAW, the cosines of
c    the pitch, roll and yaw angles.
c
c    Input, double precision SINPITCH, SINROLL, SINYAW, the sines of the pitch,
c    roll and yaw angles.
c
c    Input, double precision GLOBAS(3), the global base vector.
c
c    Input, double precision GLOPTS(3), the global coordinates
c    of the point whose coordinates are to be transformed.
c
c    Output, double precision LOCPTS(3), the local coordinates of the point
c    whose global coordinates were given in GLOPTS.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision cospitch
      double precision cosroll
      double precision cosyaw
      double precision globas(dim_num)
      double precision glopts(dim_num)
      double precision locpts(dim_num)
      double precision sinpitch
      double precision sinroll
      double precision sinyaw

      locpts(1) = ( cosyaw * cospitch ) * ( glopts(1) - globas(1) ) 
     &          + ( sinyaw * cospitch ) * ( glopts(2) - globas(2) ) 
     &          -   sinpitch * ( glopts(3) - globas(3) )

      locpts(2) = ( cosyaw * sinpitch * sinroll - sinyaw * cosroll ) 
     &  * ( glopts(1) - globas(1) ) 
     &  + ( sinyaw * sinpitch * sinroll + cosyaw * cosroll ) 
     &  * ( glopts(2) - globas(2) ) 
     &  +   cospitch * sinroll * ( glopts(3) - globas(3) )

      locpts(3) = ( cosyaw * sinpitch * cosroll + sinyaw * sinroll ) 
     &  * ( glopts(1) - globas(1) ) 
     &  + ( sinyaw * sinpitch * cosroll - cosyaw * sinroll  ) 
     &  * ( glopts(2) - globas(2) ) 
     &  + ( cospitch * cosroll ) * ( glopts(3) - globas(3) )

      return
      end
      function halfplane_contains_point_2d ( p1, p2, p )

c*********************************************************************72
c
cc HALFPLANE_CONTAINS_POINT_2D reports if a half-plane contains a point in 2d.
c
c  Discussion:
c
c    The halfplane is assumed to be all the points "to the left" of the
c    line that passes from P1 through P2.  Thus, one way to
c    understand where the point P is, is to compute the signed
c    area of the triangle ( P1, P2, P ).
c
c    If this area is
c      positive, the point is strictly inside the halfplane,
c      zero, the point is on the boundary of the halfplane,
c      negative, the point is strictly outside the halfplane.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two distinct points
c    on the line defining the half plane.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical HALFPLANE_CONTAINS_POINT_2D, is TRUE if the halfplane
c    contains the point.
c
      implicit none

      double precision area_signed
      logical halfplane_contains_point_2d
      double precision p(2)
      double precision p1(2)
      double precision p2(2)

      area_signed = 0.5D+00 *       
     &  ( p1(1) * ( p2(2) - p(2)  ) 
     &  + p2(1) * ( p(2)  - p1(2) ) 
     &  + p(1)  * ( p1(2) - p2(2) ) )

      halfplane_contains_point_2d = ( 0.0D+00 .le. area_signed )

      return
      end
      subroutine halfspace_imp_triangle_int_3d ( a, b, c, d, t, 
     &  int_num, pint )

c*********************************************************************72
c
cc HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( imp halfspace, triangle ).
c
c  Discussion:
c
c    The implicit form of a half-space in 3D may be described as the set
c    of points P on or "above" an implicit plane:
c
c      0 <= A * P(1) + B * P(2) + C * P(3) + D
c
c    The triangle is specified by listing its three vertices.
c
c    The intersection may be described by the number of vertices of the
c    triangle that are included in the halfspace, and by the location of
c    points between vertices that separate a side of the triangle into
c    an included part and an unincluded part.
c
c    0 vertices, 0 separators    (no intersection)
c    1 vertex, 0 separators      (point intersection)
c    2 vertices, 0 separators    (line intersection)
c    3 vertices, 0 separators    (triangle intersection)
c
c    1 vertex, 2 separators,     (intersection is a triangle)
c    2 vertices, 2 separators,   (intersection is a quadrilateral).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, D, the parameters that define the
c    implicit plane, which in turn define the implicit halfspace.
c
c    Input, double precision T(3,3), the vertices of the triangle.
c
c    Output, integer INT_NUM, the number of intersection points 
c    returned, which will always be between 0 and 4.
c
c    Output, double precision PINT(3,4), the coordinates of the INT_NUM
c    intersection points.  The points will lie in sequence on the triangle.
c    Some points will be vertices, and some may be separators.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision dist1
      double precision dist2
      double precision dist3
      integer int_num
      double precision pint(2,4)
      double precision t(2,3)
c
c  Compute the signed distances between the vertices and the plane.
c
      dist1 = a * t(1,1) + b * t(2,1) + c * t(3,1) + d
      dist2 = a * t(1,2) + b * t(2,2) + c * t(3,2) + d
      dist3 = a * t(1,3) + b * t(2,2) + c * t(3,3) + d
c
c  Now we can find the intersections.
c
      call halfspace_triangle_int_3d ( dist1, dist2, dist3, t, 
     &  int_num, pint )

      return
      end
      subroutine halfspace_normal_triangle_int_3d ( pp, normal, t, 
     &  int_num, pint )

c*********************************************************************72
c
cc HALFSPACE_NORMAL_TRIANGLE_INT_3D: intersection ( norm halfspace, triangle ).
c
c  Discussion:
c
c    The normal form of a halfspace in 3D may be described as the set
c    of points P on or "above" a plane described in normal form:
c
c      PP is a point on the plane,
c      NORMAL is the unit normal vector, pointing "out" of the
c      halfspace.
c
c    The triangle is specified by listing its three vertices.
c
c    The intersection may be described by the number of vertices of the
c    triangle that are included in the halfspace, and by the location of
c    points between vertices that separate a side of the triangle into
c    an included part and an unincluded part.
c
c    0 vertices, 0 separators    (no intersection)
c    1 vertex, 0 separators      (point intersection)
c    2 vertices, 0 separators    (line intersection)
c    3 vertices, 0 separators    (triangle intersection)
c
c    1 vertex, 2 separators,     (intersection is a triangle)
c    2 vertices, 2 separators,   (intersection is a quadrilateral).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the bounding plane
c    that defines the halfspace.
c
c    Input, double precision NORMAL(3), the components of the normal vector
c    to the bounding plane that defines the halfspace.  By convention, the
c    normal vector points "outwards" from the halfspace.
c
c    Input, double precision T(3,3), the vertices of the triangle.
c
c    Output, integer INT_NUM, the number of intersection points 
c    returned, which will always be between 0 and 4.
c
c    Output, double precision PINT(3,4), the coordinates of the INT_NUM
c    intersection points.  The points will lie in sequence on the triangle.
c    Some points will be vertices, and some may be separators.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision d
      double precision dist1
      double precision dist2
      double precision dist3
      double precision normal(dim_num)
      integer int_num
      double precision pp(dim_num)
      double precision pint(dim_num,4)
      double precision r8vec_dot_product
      double precision t(dim_num,3)
c
c  Compute the signed distances between the vertices and the plane.
c
      d = - r8vec_dot_product ( dim_num, normal, pp )
c
c  Compute the signed distances between the vertices and the plane.
c
      dist1 = d + r8vec_dot_product ( dim_num, normal, t(1,1) )
      dist2 = d + r8vec_dot_product ( dim_num, normal, t(1,2) )
      dist3 = d + r8vec_dot_product ( dim_num, normal, t(1,3) )
c
c  Now we can find the intersections.
c
      call halfspace_triangle_int_3d ( dist1, dist2, dist3, t, 
     &  int_num, pint )

      return
      end
      subroutine halfspace_triangle_int_3d ( dist1, dist2, dist3, t, 
     &  int_num, pint )

c*********************************************************************72
c
cc HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
c
c  Discussion:
c
c    The triangle is specified by listing its three vertices.
c
c    The halfspace is not described in the input data.  Rather, the
c    distances from the triangle vertices to the halfspace are given.
c
c    The intersection may be described by the number of vertices of the
c    triangle that are included in the halfspace, and by the location of
c    points between vertices that separate a side of the triangle into
c    an included part and an unincluded part.
c
c    0 vertices, 0 separators    (no intersection)
c    1 vertex, 0 separators      (point intersection)
c    2 vertices, 0 separators    (line intersection)
c    3 vertices, 0 separators    (triangle intersection)
c
c    1 vertex, 2 separators,     (intersection is a triangle)
c    2 vertices, 2 separators,   (intersection is a quadrilateral).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision DIST1, DIST2, DIST3, the distances from each of
c    the three vertices of the triangle to the halfspace.  The distance is
c    zero if a vertex lies within the halfspace, or on the plane that
c    defines the boundary of the halfspace.  Otherwise, it is the
c    distance from that vertex to the bounding plane.
c
c    Input, double precision T(3,3), the vertices of the triangle.
c
c    Output, integer INT_NUM, the number of intersection points
c    returned, which will always be between 0 and 4.
c
c    Output, double precision PINT(3,4), the coordinates of the INT_NUM
c    intersection points.  The points will lie in sequence on the triangle.
c    Some points will be vertices, and some may be separators.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist1
      double precision dist2
      double precision dist3
      integer i
      integer int_num
      double precision pint(dim_num,4)
      double precision t(dim_num,3)
c
c  Walk around the triangle, looking for vertices that are included,
c  and points of separation.
c
      int_num = 0

      if ( dist1 .le. 0.0D+00 ) then

        int_num = int_num + 1
        do i = 1, dim_num
          pint(i,int_num) = t(i,1)
        end do

      end if

      if ( dist1 * dist2 < 0.0D+00 ) then

        int_num = int_num + 1
        do i = 1, dim_num
          pint(i,int_num) = ( dist1 * t(i,2) - dist2 * t(i,1) ) 
     &    / ( dist1 - dist2 )
        end do

      end if

      if ( dist2 .le. 0.0D+00 ) then

        int_num = int_num + 1
        do i = 1, dim_num
          pint(i,int_num) = t(i,2)
        end do

      end if

      if ( dist2 * dist3 .lt. 0.0D+00 ) then

        int_num = int_num + 1

        do i = 1, dim_num
          pint(i,int_num) = ( dist2 * t(i,3) - dist3 * t(i,2) ) 
     &      / ( dist2 - dist3 )
        end do

      end if

      if ( dist3 .le. 0.0D+00 ) then

        int_num = int_num + 1
        do i = 1, dim_num
          pint(i,int_num) = t(i,3)
        end do

      end if

      if ( dist3 * dist1 .lt. 0.0D+00 ) then

        int_num = int_num + 1
        do i = 1, dim_num
          pint(i,int_num) = ( dist3 * t(i,1) - dist1 * t(i,3) ) 
     &      / ( dist3 - dist1 )
        end do

      end if

      return
      end
      function haversine ( a )

c*********************************************************************72
c
cc HAVERSINE computes the haversine of an angle.
c
c  Discussion:
c
c    haversine(A) = ( 1 - cos ( A ) ) / 2
c
c    The haversine is useful in spherical trigonometry.
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
c    Input, double precision A, the angle.
c
c    Output, double precision HAVERSINE, the haversine of the angle.
c
      implicit none

      double precision a
      double precision haversine

      haversine = ( 1.0D+00 - cos ( a ) ) / 2.0D+00

      return
      end
      subroutine helix_shape_3d ( a, n, r, theta1, theta2, p )

c*********************************************************************72
c
cc HELIX_SHAPE_3D computes points on a helix in 3D.
c
c  Discussion:
c
c    The user specifies the parameters A and R, the first and last
c    THETA values, and the number of equally spaced THETA values
c    at which point values are to be computed.
c
c    X = R * COS ( THETA )
c    Y = R * SIN ( THETA )
c    Z = A * THETA
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, the rate at which Z advances with THETA.
c
c    Input, integer N, the number of points to compute on 
c    the helix.
c
c    Input, double precision R, the radius of the helix.
c
c    Input, double precision THETA1, THETA2, the first and last THETA values at
c    which to compute points on the helix.  THETA is measured in
c    radians.
c
c    Output, double precision P(3,N), the coordinates of points on the helix.
c
      implicit none

      integer n

      double precision a
      integer i
      double precision p(3,n)
      double precision r
      double precision theta
      double precision theta1
      double precision theta2

      do i = 1, n

        if ( n .eq. 1 ) then
          theta = 0.5D+00 * ( theta1 + theta2 )
        else
          theta = ( dble ( n - i     ) * theta1 
     &            + dble (     i - 1 ) * theta2 ) 
     &            / dble ( n     - 1 )
        end if

        p(1,i) = r * cos ( theta )
        p(2,i) = r * sin ( theta )
        p(3,i) = a * theta

      end do

      return
      end
      function hexagon_area_2d ( r )

c*********************************************************************72
c
cc HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
c
c  Discussion:
c
c    The radius of a regular hexagon is the distance from the center
c    of the hexagon to any vertex.  This happens also to equal the
c    length of any side.
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
      function hexagon_unit_area_2d ( )

c*********************************************************************72
c
cc HEXAGON_UNIT_AREA_2D returns the area of a unit regular hexagon in 2D.
c
c  Discussion:
c
c    A "unit" regular hexagon has both a "radius" of 1 (distance
c    from the center to any vertex), and a side length of 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
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
      function i4_dedekind_factor ( p, q )

c*********************************************************************72
c
cc I4_DEDEKIND_FACTOR computes a function needed for a Dedekind sum. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hans Rademacher, Emil Grosswald,
c    Dedekind Sums,
c    Mathematics Association of America, 1972,
c    LC: QA241.R2.
c
c  Parameters:
c
c    Input, integer P, Q, two positive integers.
c
c    Input, double precision I4_DEDEKIND_FACTOR, the Dedekind factor of P / Q.
c
      implicit none

      double precision i4_dedekind_factor
      integer p
      integer q

      if ( mod ( p, q ) .eq. 0 ) then
        i4_dedekind_factor = 0.0D+00
      else
        i4_dedekind_factor = dble ( p ) / dble ( q ) 
     &    - dble ( ( p / q ) ) - 0.5D+00
      end if

      return
      end
      function i4_gcd ( i, j )

c*********************************************************************72
c
cc I4_GCD finds the greatest common divisor of I and J.
c
c  Discussion:
c
c    Only the absolute values of I and J are
c    considered, so that the result is always nonnegative.
c
c    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
c
c    If I and J have no common factor, I4_GCD is returned as 1.
c
c    Otherwise, using the Euclidean algorithm, I4_GCD is the
c    largest common factor of I and J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common divisor
c    is desired.
c
c    Output, integer I4_GCD, the greatest common divisor of I and J.
c
      implicit none

      integer i
      integer i4_gcd
      integer ip
      integer iq
      integer ir
      integer j

      i4_gcd = 1
c
c  Return immediately if either I or J is zero.
c
      if ( i .eq. 0 ) then
        i4_gcd = max ( 1, abs ( j ) )
        return
      else if ( j .eq. 0 ) then
        i4_gcd = max ( 1, abs ( i ) )
        return
      end if
c
c  Set IP to the larger of I and J, IQ to the smaller.
c  This way, we can alter IP and IQ as we go.
c
      ip = max ( abs ( i ), abs ( j ) )
      iq = min ( abs ( i ), abs ( j ) )
c
c  Carry out the Euclidean algorithm.
c
10    continue

        ir = mod ( ip, iq )

        if ( ir .eq. 0 ) then
          go to 20
        end if

        ip = iq
        iq = ir

      go to 10

20    continue

      i4_gcd = iq

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function i4_lcm ( i, j )

c*********************************************************************72
c
cc I4_LCM computes the least common multiple of two I4's.
c
c  Discussion:
c
c    The least common multiple may be defined as
c
c      LCM(I,J) = ABS( I * J ) / GCD(I,J)
c
c    where GCD(I,J) is the greatest common divisor of I and J.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the integers whose I4_LCM is desired.
c
c    Output, integer I4_LCM, the least common multiple of I and J.
c    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
c
      implicit none

      integer i
      integer i4_gcd
      integer j
      integer i4_lcm

      i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

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
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
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
c    12 November 2006
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
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

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
      subroutine i4col_find_item ( m, n, a, item, row, col )

c*********************************************************************72
c
cc I4COL_FIND_ITEM searches an I4COL for a given scalar value.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
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
c    Input, integer M, N, the number of rows and columns in
c    the table.
c
c    Input, integer A(M,N), an array of N columns of vectors
c    of length M.
c
c    Input, integer ITEM, the value to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM.  The search
c    is conducted by columns.  If the item is not found, then
c    ROW = COL = -1.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer item
      integer j
      integer row

      do j = 1, n
        do i = 1, m
          if ( a(i,j) .eq. item ) then
            row = i
            col = j
            return
          end if
        end do
      end do

      row = -1
      col = -1

      return
      end
      subroutine i4col_find_pair_wrap ( m, n, a, item1, item2, row,
     &  col )

c*********************************************************************72
c
cc I4COL_FIND_PAIR_WRAP searches an I4COL for a pair of items.
c
c  Discussion:
c
c    An I4COL is an M by N array of I4 values, regarded
c    as an array of N columns of length M.
c
c    The items (ITEM1, ITEM2) must occur consecutively.
c    However, wrapping is allowed, that is, if ITEM1 occurs
c    in the last row, and ITEM2 "follows" it in the first row
c    of the same column, a match is declared.
c
c    If the pair of items is not found, then ROW = COL = -1.
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
c    Input, integer M, N, the number of rows and columns
c    in the array.
c
c    Input, integer A(M,N), the array to search.
c
c    Input, integer ITEM1, ITEM2, the values to search for.
c
c    Output, integer ROW, COL, the row and column indices
c    of the first occurrence of the value ITEM1 followed immediately
c    by ITEM2.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer col
      integer i
      integer i2
      integer item1
      integer item2
      integer j
      integer row

      do j = 1, n
        do i = 1, m

          if ( a(i,j) .eq. item1 ) then

            i2 = i + 1

            if ( m .lt. i2 ) then
              i2 = 1
            end if

            if ( a(i2,j) .eq. item2 ) then
              row = i
              col = j
              return
            end if

          end if

        end do
      end do

      row = -1
      col = -1

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
c    Input, character*(*) TITLE, a title.
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
c    Input, character*(*) TITLE, a title.
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

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

          write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

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
      write ( *, '(a)' )  trim ( title )

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
      subroutine i4row_compare ( m, n, a, i, j, isgn )

c*********************************************************************72
c
cc I4ROW_COMPARE compares two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c  Example:
c
c    Input:
c
c    M = 3, N = 4, I = 2, J = 3
c
c    A = (
c      1  2  3  4
c      5  6  7  8
c      9 10 11 12 )
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
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of M rows of vectors
c    of length N.
c
c    Input, integer I, J, the rows to be compared.
c    I and J must be between 1 and M.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, row I .lt. row J,
c     0, row I = row J,
c    +1, row J .lt. row I.
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
c  Check that I and J are legal.
c
      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is less than 1.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      else if ( m .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is out of bounds.'
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      if ( j .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is less than 1.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      else if ( m .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is out of bounds.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      isgn = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. n ) then

        if ( a(i,k) .lt. a(j,k) ) then
          isgn = -1
          return
        else if ( a(j,k) .lt. a(i,k) ) then
          isgn = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4row_sort_a ( m, n, a )

c*********************************************************************72
c
cc I4ROW_SORT_A ascending sorts the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c    In lexicographic order, the statement "X .lt. Y", applied to two
c    vectors X and Y of length M, means that there is some index I, with
c    1 .le. I .le. M, with the property that
c
c      X(J) = Y(J) for J .lt. I,
c    and
c      X(I) .lt. Y(I).
c
c    In other words, X is less than Y if, at the first index where they
c    differ, the X value is less than the Y value.
c
c  Example:
c
c    Input:
c
c      M = 5, N = 3
c
c      A =
c        3  2  1
c        2  4  3
c        3  1  8
c        2  4  2
c        1  9  9
c
c    Output:
c
c      A =
c        1  9  9
c        2  4  2
c        2  4  3
c        3  1  8
c        3  2  1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of M rows of N-vectors.
c    On output, the rows of A have been sorted in ascending
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

      if ( m .le. 1 ) then
        return
      end if

      if ( n .le. 0 ) then
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

        call sort_heap_external ( m, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4row_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4row_compare ( m, n, a, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4row_sorted_unique_count ( m, n, a, unique_num )

c*********************************************************************72
c
cc I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of I4 values, regarded
c    as an array of M rows of length N.
c
c    The rows of the array may be ascending or descending sorted.
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
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), a sorted array, containing
c    M rows of data.
c
c    Output, integer UNIQUE_NUM, the number of unique rows.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      logical equal
      integer i1
      integer i2
      integer j
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if

      unique_num = 1
      i1 = 1

      do i2 = 2, m

        equal = .true.

        do j = 1, n
          if ( a(i1,j) .ne. a(i2,j) ) then
            equal = .false.
            go to 10
          end if
        end do

10      continue

        if ( .not. equal ) then
          unique_num = unique_num + 1
          i1 = i2
        end if

      end do

      return
      end
      subroutine i4row_swap ( m, n, a, i1, i2 )

c*********************************************************************72
c
cc I4ROW_SWAP swaps two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), an array of data.
c
c    Input, integer I1, I2, the two rows to swap.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i1
      integer i2
      integer row(n)
c
c  Check.
c
      if ( i1 .lt. 1 .or. m .lt. i1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I1 is out of range.'
        stop
      end if

      if ( i2 .lt. 1 .or. m .lt. i2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I2 is out of range.'
        stop
      end if

      if ( i1 .eq. i2 ) then
        return
      end if

      row(1:n)  = a(i1,1:n)
      a(i1,1:n) = a(i2,1:n)
      a(i2,1:n) = row(1:n)

      return
      end
      subroutine i4vec_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_D reorders an I4VEC into an descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
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
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
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
c    An I4VEC is a vector of I4's.
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
      function i4vec_lcm ( n, v )

c*********************************************************************72
c
cc I4VEC_LCM returns the least common multiple of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The value LCM returned has the property that it is the smallest integer
c    which is evenly divisible by every element of V.
c
c    The entries in V may be negative.
c
c    If any entry of V is 0, then LCM is 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of V.
c
c    Input, integer V(N), the vector.
c
c    Output, integer I4VEC_LCM, the least common multiple of V.
c
      implicit none

      integer n

      integer i
      integer i4_lcm
      integer i4vec_lcm
      integer lcm
      integer v(n)

      lcm = 1

      do i = 1, n

        if ( v(i) .eq. 0 ) then
          lcm = 0
          i4vec_lcm = lcm
          return
        end if

        lcm = i4_lcm ( lcm, v(i) )

      end do

      i4vec_lcm = lcm

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
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      function i4vec_product ( n, a )

c*********************************************************************72
c
cc I4VEC_PRODUCT returns the product of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c    In FORTRAN90, this facility is offered by the built in
c    PRODUCT function:
c
c      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
c
c    In MATLAB, this facility is offered by the built in
c    PROD function:
c
c      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2007
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
c    Output, integer I4VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_product

      i4vec_product = 1
      do i = 1, n
        i4vec_product = i4vec_product * a(i)
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
c    An I4VEC is a vector of I4's.
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
      subroutine i4vec_uniform ( n, a, b, seed, x )

c*********************************************************************72
c
cc I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The pseudorandom numbers should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer X(N), a vector of numbers between A and B.
c
      implicit none

      integer n

      integer a
      integer b
      integer i
      integer k
      real r
      integer seed
      integer value
      integer x(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
        r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
        value = nint ( r )

        value = max ( value, min ( a, b ) )
        value = min ( value, max ( a, b ) )

        x(i) = value

      end do

      return
      end
      subroutine i4vec_zero ( n, a )

c*********************************************************************72
c
cc I4VEC_ZERO sets the entries of an I4VEC to 0.
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
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Output, integer A(N), the vector, which has been set to zero.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = 0
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
      function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

c*********************************************************************72
c
cc LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
c
c  Discussion:
c
c    The explicit form of a line in ND is:
c
c      the line through the points P1 and P2.
c
c    An explicit line is degenerate if the two defining points are equal.
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
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision P1(DIM_NUM), P2(DIM_NUM), two points on the line.
c
c    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
c    is degenerate.
c
      implicit none

      integer dim_num

      integer i
      logical line_exp_is_degenerate_nd
      double precision p1(dim_num)
      double precision p2(dim_num)

      line_exp_is_degenerate_nd = .false.

      do i = 1, dim_num
        if ( p1(i) .ne. p2(i) ) then
          return
        end if
      end do

      line_exp_is_degenerate_nd = .true.

      return
      end
      subroutine line_exp_normal_2d ( p1, p2, normal )

c*********************************************************************72
c
cc LINE_EXP_NORMAL_2D computes a unit normal vector to a line in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The sign of the normal vector N is chosen so that the normal vector
c    points "to the left" of the direction of the line.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two distinct points on the line.
c
c    Output, double precision NORMAL(2), a unit normal vector to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical line_exp_is_degenerate_nd
      double precision norm
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r8vec_diff_norm

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        normal(1) = sqrt ( 2.0D+00 )
        normal(2) = sqrt ( 2.0D+00 )
        return
      end if

      norm = r8vec_diff_norm ( 2, p2, p1 )

      normal(1) = - ( p2(2) - p1(2) ) / norm
      normal(2) =   ( p2(1) - p1(1) ) / norm

      return
      end
      subroutine line_exp_perp_2d ( p1, p2, p3, p4, flag )

c*********************************************************************72
c
cc LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The input point P3 should NOT lie on the line (P1,P2).  If it
c    does, then the output value P4 will equal P3.
c
c    P1-----P4-----------P2
c            |
c            |
c           P3
c
c    P4 is also the nearest point on the line (P1,P2) to the point P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Input, double precision P3(2), a point (presumably not on the 
c    line (P1,P2)), through which the perpendicular must pass.
c
c    Output, double precision P4(2), a point on the line (P1,P2),
c    such that the line (P3,P4) is perpendicular to the line (P1,P2).
c
c    Output, logical FLAG, is TRUE if the value could not be computed.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision bot
      logical flag
      logical line_exp_is_degenerate_nd
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision r8_huge
      double precision r8vec_diff_norm_squared
      double precision t

      flag = .false.

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        flag = .true.
        p4(1) = r8_huge ( )
        p4(2) = r8_huge ( )
        return
      end if

      bot = r8vec_diff_norm_squared ( 2, p2, p1 )
c
c  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
c
c  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
c  of the projection of (P3-P1) onto (P2-P1).
c
      t = ( ( p1(1) - p3(1) ) * ( p1(1) - p2(1) )
     &    + ( p1(2) - p3(2) ) * ( p1(2) - p2(2) ) ) 
     &    / bot

      p4(1) = p1(1) + t * ( p2(1) - p1(1) )
      p4(2) = p1(2) + t * ( p2(2) - p1(2) )

      return
      end
      subroutine line_exp_point_dist_2d ( p1, p2, p, dist )

c*********************************************************************72
c
cc LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
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
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Input, double precision P(2), the point whose distance from the line is
c    to be measured.
c
c    Output, double precision DIST, the distance from the point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision bot
      double precision dist
      double precision dot
      logical line_exp_is_degenerate_nd
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision t

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

        pn(1) = p1(1)
        pn(2) = p1(2)
c
c  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
c
c  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
c  of the projection of (P-P1) onto (P2-P1).
c
      else

        dot = ( p(1) - p1(1) ) * ( p2(1) - p1(1) )
     &      + ( p(2) - p1(2) ) * ( p2(2) - p1(2) )

        bot = ( p2(1) - p1(1) )**2 + ( p2(2) - p1(2) )**2

        t = dot / bot

        pn(1) = p1(1) + t * ( p2(1) - p1(1) )
        pn(2) = p1(2) + t * ( p2(2) - p1(2) )

      end if

      dist = sqrt ( ( p(1) - pn(1) )**2 + ( p(2) - pn(2) )**2 )

      return
      end
      subroutine line_exp_point_dist_3d ( p1, p2, p, dist )

c*********************************************************************72
c
cc LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
c
c  Discussion:
c
c    The explicit form of a line in 3D is:
c
c      the line through the points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), two points on the line.
c
c    Input, double precision P(3), the point whose distance from the line is
c    to be measured.
c
c    Output, double precision DIST, the distance from the point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision bot
      double precision dist
      integer i
      logical line_exp_is_degenerate_nd
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision r8vec_diff_norm
      double precision r8vec_diff_norm_squared
      double precision t

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

        call r8vec_copy ( dim_num, p1, pn )
c
c  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
c
c  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
c  of the projection of (P-P1) onto (P2-P1).
c
      else

        bot = r8vec_diff_norm_squared ( dim_num, p2, p1 )

        t = 0.0D+00
        do i = 1, dim_num
          t = t + ( p(i) - p1(i) ) * ( p2(i) - p1(i) )
        end do
        t = t / bot

        do i = 1, dim_num
          pn(i) = p1(i) + t * ( p2(i) - p1(i) )
        end do

      end if
c
c  Now compute the distance between the projection point and P.
c
      dist = r8vec_diff_norm ( dim_num, p, pn )

      return
      end
      subroutine line_exp_point_dist_signed_2d ( p1, p2, p, 
     &  dist_signed )

c*********************************************************************72
c
cc LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( exp line, point ) in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The signed distance has two interesting properties:
c
c    *  The absolute value of the signed distance is the
c        usual (Euclidean) distance.
c
c    *  Points with signed distance 0 lie on the line,
c       points with a negative signed distance lie on one side
c         of the line,
c       points with a positive signed distance lie on the
c         other side of the line.
c
c    Assuming that C is nonnegative, then if a point is a positive
c    distance away from the line, it is on the same side of the
c    line as the point (0,0), and if it is a negative distance
c    from the line, it is on the opposite side from (0,0).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Input, double precision P(2), the point whose signed distance is desired.
c
c    Output, double precision DIST_SIGNED, the signed distance from the
c    point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision dist_signed
      logical line_exp_is_degenerate_nd
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r8vec_diff_norm
c
c  If the explicit line degenerates to a point, the computation is easy.
c
      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

        dist_signed = r8vec_diff_norm ( dim_num, p1, p )
c
c  Convert the explicit line to the implicit form A * P(1) + B * P(2) + C = 0.
c  This makes the computation of the signed distance to (X,Y) easy.
c
      else

        a = p2(2) - p1(2)
        b = p1(1) - p2(1)
        c = p2(1) * p1(2) - p1(1) * p2(2)

        dist_signed = ( a * p(1) + b * p(2) + c ) 
     &    / sqrt ( a * a + b * b )

      end if

      return
      end
      subroutine line_exp_point_near_2d ( p1, p2, p, pn, dist, t )

c*********************************************************************72
c
cc LINE_EXP_POINT_NEAR_2D: point on an explicit line nearest a point in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The nearest point PN = (XN,YN) has the form:
c
c      PN = (1-T) * P1 + T * P2.
c
c    If T is less than 0, PN is furthest from P2.
c    If T is between 0 and 1, PN is between P1 and P2.
c    If T is greater than 1, PN is furthest from P1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Input, double precision P(2), the point whose nearest neighbor on the
c    line is to be determined.
c
c    Output, double precision PN(2), the nearest point on the line to P.
c
c    Output, double precision DIST, the distance from the point to the line.
c
c    Output, double precision T, the relative position of the point
c    PN to the points P1 and P2.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision bot
      double precision dist
      integer i
      logical line_exp_is_degenerate_nd
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision r8vec_diff_dot_product
      double precision r8vec_diff_norm
      double precision r8vec_diff_norm_squared
      double precision t

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_POINT_NEAR_2D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if
c
c  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
c
c  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
c  of the projection of (P-P1) onto (P2-P1).
c
      bot = r8vec_diff_norm_squared ( dim_num, p2, p1 )

      t = r8vec_diff_dot_product ( p1, p, p1, p2 ) / bot

      do i = 1, dim_num
        pn(i) = p1(i) + t * ( p2(i) - p1(i) )
      end do

      dist = r8vec_diff_norm ( dim_num, pn, p )

      return
      end
      subroutine line_exp_point_near_3d ( p1, p2, p, pn, dist, t )

c*********************************************************************72
c
cc LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
c
c  Discussion:
c
c    The explicit form of a line in 3D is:
c
c      the line through the points P1 and P2.
c
c    The nearest point PN has the form:
c
c      PN = ( 1 - T ) * P1 + T * P2.
c
c    If T is less than 0, PN is furthest away from P2.
c    If T is between 0 and 1, PN is between P1 and P2.
c    If T is greater than 1, PN is furthest away from P1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), two points on the line.
c
c    Input, double precision P(3), the point whose nearest neighbor on
c    the line is to be determined.
c
c    Output, double precision PN(3), the point which is the nearest
c    point on the line to P.
c
c    Output, double precision DIST, the distance from the point to the 
c    nearest point on the line.
c
c    Output, double precision T, the relative position of the point
c    PN to P1 and P2.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision bot
      double precision dist
      integer i
      logical line_exp_is_degenerate_nd
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision r8vec_diff_dot_product
      double precision r8vec_diff_norm
      double precision r8vec_diff_norm_squared
      double precision t

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if
c
c  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
c
c  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
c  of the projection of (P-P1) onto (P2-P1).
c
      bot = r8vec_diff_norm_squared ( dim_num, p1, p2 )

      t = r8vec_diff_dot_product ( p, p1, p2, p1 ) / bot
c
c  Now compute the location of the projection point.
c
      do i = 1, dim_num
        pn(i) = p1(i) + t * ( p2(i) - p1(i) )
      end do
c
c  Now compute the distance between the projection point and P.
c
      dist = r8vec_diff_norm ( dim_num, pn, p )

      return
      end
      subroutine line_exp2imp_2d ( p1, p2, a, b, c )

c*********************************************************************72
c
cc LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The implicit form of a line in 2D is:
c
c      A * X + B * Y + C = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Output, double precision A, B, C, the implicit form of the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      logical line_exp_is_degenerate_nd
      double precision norm
      double precision p1(dim_num)
      double precision p2(dim_num)
c
c  Take care of degenerate cases.
c
      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
        write ( *, '(a)' ) '  The line is degenerate.'
      end if

      a = p2(2) - p1(2)
      b = p1(1) - p2(1)
      c = p2(1) * p1(2) - p1(1) * p2(2)

      norm = a * a + b * b + c * c

      if ( 0.0D+00 .lt. norm ) then
        a = a / norm
        b = b / norm
        c = c / norm
      end if

      if ( a .lt. 0.0D+00 ) then
        a = -a
        b = -b
        c = -c
      end if

      return
      end
      subroutine line_exp2par_2d ( p1, p2, f, g, x0, y0 )

c*********************************************************************72
c
cc LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F**2 + G**2 = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), two points on the line.
c
c    Output, double precision F, G, X0, Y0, the parametric parameters
c    of the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision f
      double precision g
      logical line_exp_is_degenerate_nd
      double precision norm
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision x0
      double precision y0

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_EXP2PAR_2D - Warning!'
        write ( *, '(a)' ) '  The line is degenerate.'
      end if

      x0 = p1(1)
      y0 = p1(2)

      f = p2(1) - p1(1)
      g = p2(2) - p1(2)

      norm = sqrt ( f * f + g * g )

      if ( norm .ne. 0.0D+00 ) then
        f = f / norm
        g = g / norm
      end if

      if ( f .lt. 0.0D+00 ) then
        f = -f
        g = -g
      end if

      return
      end
      subroutine line_exp2par_3d ( p1, p2, f, g, h, x0, y0, z0 )

c*********************************************************************72
c
cc LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
c
c  Discussion:
c
c    The explicit form of a line in 3D is:
c
c      the line through the points P1 and P2.
c
c    The parametric form of a line in 3D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c      Z = Z0 + H * T
c
c    We normalize by always choosing F*F + G*G + H*H = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), two points on the line.
c
c    Output, double precision F, G, H, X0, Y0, Z0, the parametric parameters
c    of the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision f
      double precision g
      double precision h
      logical line_exp_is_degenerate_nd
      double precision norm
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision x0
      double precision y0
      double precision z0

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_EXP2PAR_3D - Warning!'
        write ( *, '(a)' ) '  The line is degenerate.'
      end if

      x0 = p1(1)
      y0 = p1(2)
      z0 = p1(3)

      f = p2(1) - p1(1)
      g = p2(2) - p1(2)
      h = p2(3) - p1(3)

      norm = sqrt ( f * f + g * g + h * h )

      if ( norm .ne. 0.0D+00 ) then
        f = f / norm
        g = g / norm
        h = h / norm
      end if

      if ( f .lt. 0.0D+00 ) then
        f = -f
        g = -g
        h = -h
      end if

      return
      end
      function line_imp_is_degenerate_2d ( a, b, c )

c*********************************************************************72
c
cc LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
c
c  Discussion:
c
c    The implicit form of a line in 2D is:
c
c      A * X + B * Y + C = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit line parameters.
c
c    Output, logical LINE_IMP_IS_DEGENERATE_2D, is true if the
c    line is degenerate.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      logical line_imp_is_degenerate_2d

      line_imp_is_degenerate_2d = ( a * a + b * b .eq. 0.0D+00 )

      return
      end
      subroutine line_imp_point_dist_2d ( a, b, c, p, dist )

c*********************************************************************72
c
cc LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
c
c  Discussion:
c
c    The implicit form of a line in 2D is:
c
c      A * X + B * Y + C = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit line parameters.
c
c    Input, double precision P(2), the point whose distance from the line is
c    to be measured.
c
c    Output, double precision DIST, the distance from the point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision dist
      logical line_imp_is_degenerate_2d
      double precision p(dim_num)

      if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_2D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if

      dist = abs ( a * p(1) + b * p(2) + c ) / sqrt ( a * a + b * b )

      return
      end
      subroutine line_imp_point_dist_signed_2d ( a, b, c, p, 
     &  dist_signed )

c*********************************************************************72
c
cc LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( imp line, point ) in 2D.
c
c  Discussion:
c
c    The implicit form of a line in 2D is:
c
c      A * X + B * Y + C * Z + D = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit line parameters.
c
c    Input, double precision P(2), the coordinates of the point.
c
c    Output, double precision DIST_SIGNED, the signed distance from the
c    point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision dist_signed
      logical line_imp_is_degenerate_2d
      double precision p(dim_num)

      if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'LINE_IMP_POINT_DIST_SIGNED_2D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if

      dist_signed = - sign ( 1.0D+00, c ) 
     &  * ( a * p(1) + b * p(2) + c ) / 
     &  sqrt ( a * a + b * b )

      return
      end
      subroutine line_imp2exp_2d ( a, b, c, p1, p2 )

c*********************************************************************72
c
cc LINE_IMP2EXP_2D converts an implicit line to explicit form in 2D.
c
c  Discussion:
c
c    The implicit form of line in 2D is:
c
c      A * X + B * Y + C = 0
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit line parameters.
c
c    Output, double precision P1(2), P2(2), two points on the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      logical line_imp_is_degenerate_2d
      double precision normsq
      double precision p1(dim_num)
      double precision p2(dim_num)

      if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_IMP2EXP_2D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if

      normsq = a * a + b * b

      p1(1) = - a * c / normsq
      p1(2) = - b * c / normsq

      if ( abs ( b ) .lt. abs ( a ) ) then
        p2(1) = - ( a - b / a ) * c / normsq
        p2(2) = - ( b + 1.0D+00 ) * c / normsq
      else
        p2(1) = - ( a + 1.0D+00 ) * c / normsq
        p2(2) = - ( b - a / b ) * c / normsq
      end if

      return
      end
      subroutine line_imp2par_2d ( a, b, c, f, g, x0, y0 )

c*********************************************************************72
c
cc LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
c
c  Discussion:
c
c    The implicit form of line in 2D is:
c
c      A * X + B * Y + C = 0
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit line parameters.
c
c    Output, double precision F, G, X0, Y0, the parametric parameters of
c    the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision f
      logical line_imp_is_degenerate_2d
      double precision g
      double precision x0
      double precision y0

      if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_IMP2PAR_2D - Fatal error!'
        write ( *, '(a)' ) '  The line is degenerate.'
        stop
      end if

      x0 = - a * c / ( a * a + b * b )
      y0 = - b * c / ( a * a + b * b )

      f =   b / sqrt ( a * a + b * b )
      g = - a / sqrt ( a * a + b * b )

      if ( f .lt. 0.0D+00 ) then
        f = -f
        g = -g
      end if

      return
      end
      subroutine line_par_point_dist_2d ( f, g, x0, y0, p, dist )

c*********************************************************************72
c
cc LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
c
c  Discussion:
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision F, G, X0, Y0, the parametric line parameters.
c
c    Input, double precision P(2), the point whose distance from the line is
c    to be measured.
c
c    Output, double precision DIST, the distance from the point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dist
      double precision dx
      double precision dy
      double precision f
      double precision g
      double precision p(dim_num)
      double precision x0
      double precision y0

      dx =   g * g * ( p(1) - x0 ) - f * g * ( p(2) - y0 )
      dy = - f * g * ( p(1) - x0 ) + f * f * ( p(2) - y0 )

      dist = sqrt ( dx * dx + dy * dy ) / ( f * f + g * g )

      return
      end
      subroutine line_par_point_dist_3d ( f, g, h, x0, y0, z0, p, dist )

c*********************************************************************72
c
cc LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
c
c  Discussion:
c
c    The parametric form of a line in 3D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c      Z = Z0 + H * T
c
c    We normalize by always choosing F*F + G*G + H*H = 1, and F nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision F, G, H, X0, Y0, Z0, the parametric line
c    parameters.
c
c    Input, double precision P(3), the point whose distance from the line is
c    to be measured.
c
c    Output, double precision DIST, the distance from the point to the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dist
      double precision dx
      double precision dy
      double precision dz
      double precision f
      double precision g
      double precision h
      double precision p(dim_num)
      double precision x0
      double precision y0
      double precision z0

      dx =   g * ( f * ( p(2) - y0 ) - g * ( p(1) - x0 ) ) 
     &     + h * ( f * ( p(3) - z0 ) - h * ( p(1) - x0 ) )

      dy =   h * ( g * ( p(3) - z0 ) - h * ( p(2) - y0 ) )  
     &     - f * ( f * ( p(2) - y0 ) - g * ( p(1) - x0 ) )

      dz = - f * ( f * ( p(3) - z0 ) - h * ( p(1) - x0 ) ) 
     &     - g * ( g * ( p(3) - z0 ) - h * ( p(2) - y0 ) )

      dist = sqrt ( dx * dx + dy * dy + dz * dz ) 
     &  / ( f * f + g * g + h * h )

      return
      end
      subroutine line_par2exp_2d ( f, g, x0, y0, p1, p2 )

c*********************************************************************72
c
cc LINE_PAR2EXP_2D converts a parametric line to explicit form in 2D.
c
c  Discussion:
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision F, G, X0, Y0, the parametric line parameters.
c
c    Output, double precision P1(2), P2(2), two points on the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision f
      double precision g
      double precision x0
      double precision y0
      double precision p1(dim_num)
      double precision p2(dim_num)

      p1(1) = x0
      p1(2) = y0

      p2(1) = p1(1) + f
      p2(2) = p1(2) + g

      return
      end
      subroutine line_par2imp_2d ( f, g, x0, y0, a, b, c )

c*********************************************************************72
c
cc LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
c
c  Discussion:
c
c    The parametric form of a line in 2D is:
c
c      X = X0 + F * T
c      Y = Y0 + G * T
c
c    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
c
c    The implicit form of a line in 2D is:
c
c      A * X + B * Y + C = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision F, G, X0, Y0, the parametric line parameters.
c
c    Output, double precision A, B, C, the implicit line parameters.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision f
      double precision g
      double precision x0
      double precision y0

      a = -g
      b = f
      c = g * x0 - f * y0

      return
      end
      subroutine lines_exp_angle_3d ( p1, p2, q1, q2, angle )

c*********************************************************************72
c
cc LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
c
c  Discussion:
c
c    The explicit form of a line in 3D is:
c
c      the line through the points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), two points on the first line.
c
c    Input, double precision Q1(3), Q2(3), two points on the second line.
c
c    Output, double precision ANGLE, the angle in radians between the two
c    lines.  The angle is computed using the ACOS function, and so lies between
c    0 and PI.  But if one of the lines is degenerate, the angle is 
c    returned as -1.0.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision angle
      double precision ctheta
      logical line_exp_is_degenerate_nd
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pdotq
      double precision pnorm
      double precision q1(dim_num)
      double precision q2(dim_num)
      double precision qnorm
      double precision r8_acos
      double precision r8vec_diff_dot
      double precision r8vec_diff_norm

      if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
c       write ( *, '(a)' ) ' '
c       write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
c       write ( *, '(a)' ) '  The line (P1,P2) is degenerate!'
        angle = -1.0D+00
        return
      end if

      if ( line_exp_is_degenerate_nd ( dim_num, q1, q2 ) ) then
c       write ( *, '(a)' ) ' '
c       write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Warning!'
c       write ( *, '(a)' ) '  The line (Q1,Q2) is degenerate!'
        angle = -1.0D+00
        return
      end if

      pnorm = r8vec_diff_norm ( dim_num, p2, p1 )

      qnorm = r8vec_diff_norm ( dim_num, q2, q1 )

      pdotq = r8vec_diff_dot ( dim_num, p2, p1, q2, q1 )

      ctheta = pdotq / ( pnorm * qnorm )

      angle = r8_acos ( ctheta )

      return
      end
      function lines_exp_parallel_2d ( p1, p2, q1, q2 )

c*********************************************************************72
c
cc LINES_EXP_PARALLEL_2D determines if two lines are parallel in 2D.
c
c  Discussion:
c
c    The explicit form of a line in 2D is:
c
c      the line through the points P1 and P2.
c
c    The test is essentially a comparison of slopes, but should be
c    more accurate than an explicit slope comparison, and unfazed
c    by degenerate cases.
c
c    On the other hand, there is NO tolerance for error.  If the
c    slopes differ by a single digit in the last place, then the
c    lines are judged to be nonparallel.  A more robust test would
c    be to compute the angle between the lines, because then it makes
c    sense to say the lines are "almost" parallel: the angle is small.
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
c    Input, double precision P1(2), P2(2), two points on the first line.
c
c    Input, double precision Q1(2), Q2(2), two points on the second line.
c
c    Output, logical LINES_EXP_PARALLEL_2D is TRUE if the lines are parallel.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical lines_exp_parallel_2d
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision q1(dim_num)
      double precision q2(dim_num)

      lines_exp_parallel_2d = 
     &  ( p2(1) - p1(1) ) * ( q2(2) - q1(2) ) .eq. 
     &  ( q2(1) - q1(1) ) * ( p2(2) - p1(2) )

      return
      end
      function lines_exp_parallel_3d ( p1, p2, q1, q2 )

c*********************************************************************72
c
cc LINES_EXP_PARALLEL_3D determines if two lines are parallel in 3D.
c
c  Discussion:
c
c    The explicit form of a line in 3D is:
c
c      the line through the points P1 and P2.
c
c    The points P1, P2 define a direction (P2-P1).  Similarly, the
c    points (Q1,Q2) define a direction (Q2-Q1).  The quantity
c
c      (P2-P1) dot (Q2-Q1) = norm(P2-P1) * norm(Q2-Q1) * cos ( angle )
c
c    Therefore, the following value is between 0 and 1;
c
c      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) )
c
c    and the lines are parallel if
c
c      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) ) = 1
c
c    We can rephrase this as requiring:
c
c      ( (P2-P1)dot(Q2-Q1) )**2 = (P2-P1)dot(P2-P1) * (Q2-Q1)dot(Q2-Q1)
c
c    which avoids division and square roots.
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
c    Input, double precision P1(3), P2(3), two points on the first line.
c
c    Input, double precision Q1(3), Q2(3), two points on the second line.
c
c    Output, logical LINES_EXP_PARALLEL_3D is TRUE if the lines are parallel.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      logical lines_exp_parallel_3d
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pdotp
      double precision pdotq
      double precision q1(dim_num)
      double precision q2(dim_num)
      double precision qdotq
      double precision r8vec_diff_dot_product

      pdotq = r8vec_diff_dot_product ( dim_num, p2, p1, q2, q1 )
      pdotp = r8vec_diff_dot_product ( dim_num, p2, p1, p2, p1 )
      qdotq = r8vec_diff_dot_product ( dim_num, q2, q1, q2, q1 )

      lines_exp_parallel_3d = ( pdotq * pdotq .eq. pdotp * qdotq )

      return
      end
      subroutine octahedron_size_3d ( point_num, edge_num, face_num,
     &  face_order_max )

c*********************************************************************72
c
cc OCTAHEDRON_SIZE_3D returns size information for an octahedron in 3D.
c
c  Discussion:
c
c    This routine can be called before calling OCTAHEDRON_SHAPE_3D,
c    so that space can be allocated for the arrays.
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
c    Output, integer FACE_ORDER_MAX, the maximum number of
c    vertices per face.
c
      implicit none

      integer edge_num
      integer face_num
      integer face_order_max
      integer point_num

      point_num = 6
      edge_num = 12
      face_num = 8
      face_order_max = 3

      return
      end
      subroutine parallelogram_area_2d ( p, area )

c*********************************************************************72
c
cc PARALLELOGRAM_AREA_2D computes the area of a parallelogram in 2D.
c
c  Discussion:
c
c    A parallelogram is a polygon having four sides, with the property
c    that each pair of opposite sides is paralell.
c
c    Given the first three vertices of the parallelogram,
c    P1, P2, and P3, the fourth vertex must satisfy
c
c      P4 = P1 + ( P3 - P2 )
c
c    This routine uses the fact that the norm of the cross product
c    of two vectors is the area of the parallelogram they form:
c
c      Area = ( P3 - P2 ) x ( P1 - P2 ).
c
c        P4<-----P3
c        /       /
c       /       /
c      P1----->P2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P(2,4), the parallelogram vertices,
c    given in counterclockwise order.  The fourth vertex is ignored.
c
c    Output, double precision AREA, the (signed) area.
c
      implicit none

      double precision area
      double precision p(2,4)
c
c  Compute the cross product vector, which only has a single
c  nonzero component.
c
      area = ( p(1,2) - p(1,1) ) * ( p(2,3) - p(2,1) )
     &     - ( p(2,2) - p(2,1) ) * ( p(1,3) - p(1,1) )

      return
      end
      subroutine parallelogram_area_3d ( p, area )

c*********************************************************************72
c
cc PARALLELOGRAM_AREA_3D computes the area of a parallelogram in 3D.
c
c  Discussion:
c
c    A parallelogram is a polygon having four sides, with the property
c    that each pair of opposite sides is paralell.
c
c    A parallelogram in 3D must have the property that it is "really"
c    a 2D object, that is, that the four vertices that define it lie
c    in some plane.
c
c    Given the first three vertices of the parallelogram (in 2D or 3D),
c    P1, P2, and P3, the fourth vertex must satisfy
c
c      P4 = P1 + ( P3 - P2 )
c
c    This routine uses the fact that the norm of the cross product
c    of two vectors is the area of the parallelogram they form:
c
c      Area = ( P3 - P2 ) x ( P1 - P2 ).
c
c        P4<-----P3
c        /       /
c       /       /
c      P1----->P2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P(3,4), the parallelogram vertices,
c    given in counterclockwise order.  The fourth vertex is ignored.
c
c    Output, double precision AREA, the area
c
      implicit none

      double precision area
      double precision cross(3)
      double precision p(3,4)
c
c  Compute the cross product vector.
c
      cross(1) = ( p(2,2) - p(2,1) ) * ( p(3,3) - p(3,1) )
     &         - ( p(3,2) - p(3,1) ) * ( p(2,3) - p(2,1) )

      cross(2) = ( p(3,2) - p(3,1) ) * ( p(1,3) - p(1,1) )
     &         - ( p(1,2) - p(1,1) ) * ( p(3,3) - p(3,1) )

      cross(3) = ( p(1,2) - p(1,1) ) * ( p(2,3) - p(2,1) )
     &         - ( p(2,2) - p(2,1) ) * ( p(1,3) - p(1,1) )

      area = sqrt ( cross(1)**2 + cross(2)**2 + cross(3)**2 )

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
      call perm_check ( n, p, base, ierror )

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
      subroutine plane_exp_normal_3d ( p1, p2, p3, normal )

c*********************************************************************72
c
cc PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
c
c  Discussion:
c
c    The explicit form of a plane in 3D is
c
c      the plane through P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), three points on the plane.
c
c    Output, double precision NORMAL(3), the coordinates of the unit normal
c    vector to the plane containing the three points.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision normal(dim_num)
      double precision normal_norm
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
c
c  The cross product (P2-P1) x (P3-P1) is normal to (P2-P1) and (P3-P1).
c
      normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) )
     &          - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

      normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) )
     &          - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

      normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) )
     &          - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

      normal_norm = 0.0D+00
      do i = 1, 3
        normal_norm = normal_norm + normal(i)**2
      end do
      normal_norm = sqrt ( normal_norm )

      if ( normal_norm .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
        write ( *, '(a)' ) '  The plane is poorly defined.'
        stop
      end if

      do i = 1, 3
        normal(i) = normal(i) / normal_norm
      end do

      return
      end
      subroutine plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

c*********************************************************************72
c
cc PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
c
c  Discussion:
c
c    The explicit form of a plane in 3D is
c
c      the plane through P1, P2 and P3.
c
c    The implicit form of a plane in 3D is
c
c      A * X + B * Y + C * Z + D = 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), three points on the plane.
c
c    Output, double precision A, B, C, D, coefficients which describe
c    the plane.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)

      a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) )
     &  - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

      b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) )
     &  - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

      c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) )
     &  - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

      d = - p2(1) * a - p2(2) * b - p2(3) * c

      return
      end
      subroutine plane_exp2normal_3d ( p1, p2, p3, pp, normal )

c*********************************************************************72
c
cc PLANE_EXP2NORMAL_3D converts an explicit plane to normal form in 3D.
c
c  Discussion:
c
c    The explicit form of a plane in 3D is
c
c      the plane through P1, P2 and P3.
c
c    The normal form of a plane in 3D is
c
c      PP, a point on the plane, and
c      N, the unit normal to the plane.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), three points on the plane.
c
c    Output, double precision PP(3), a point on the plane.
c
c    Output, double precision NORMAL(3), a unit normal vector to the plane.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision norm
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pp(dim_num)

      do i = 1, 3
        pp(i) = p1(i)
      end do

      normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) )
     &          - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

      normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) )
     &          - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

      normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) )
     &          - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

      norm = 0.0D+00
      do i = 1, 3
        norm = norm + normal(i)**2
      end do
      norm = sqrt ( norm )

      if ( norm .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_EXP2NORMAL_3D - Fatal error!'
        write ( *, '(a)' ) '  The normal vector is null.'
        write ( *, '(a)' ) '  Two points coincide, or nearly so.'
        stop
      end if

      do i = 1, 3
        normal(i) = normal(i) / norm
      end do

      return
      end
      function plane_imp_is_degenerate_3d ( a, b, c )

c*********************************************************************72
c
cc PLANE_IMP_IS_DEGENERATE_3D is TRUE if an implicit plane is degenerate.
c
c  Discussion:
c
c    The implicit form of a plane in 3D is:
c
c      A * X + B * Y + C * Z + D = 0
c
c    The implicit plane is degenerate if A = B = C = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the implicit plane parameters.
c
c    Output, logical PLANE_IMP_IS_DEGENERATE_3D, is TRUE if the plane 
c    is degenerate.
c
      implicit none

      double precision a
      double precision b
      double precision c
      logical plane_imp_is_degenerate_3d

      if ( a .eq. 0.0D+00 .and. 
     &     b .eq. 0.0D+00 .and. 
     &     c .eq. 0.0D+00 ) then
        plane_imp_is_degenerate_3d = .true.
      else
        plane_imp_is_degenerate_3d = .false.
      end if

      return
      end
      subroutine plane_imp2exp_3d ( a, b, c, d, p1, p2, p3 )

c*********************************************************************72
c
cc PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
c
c  Discussion:
c
c    The implicit form of a plane in 3D is
c
c      A * X + B * Y + C * Z + D = 0.
c
c    The explicit form of a plane in 3D is
c
c      the plane through P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, D, the implicit plane parameters.
c
c    Output, double precision P1(3), P2(3), P3(3), three points on the plane.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pp(dim_num)

      call plane_imp2normal_3d ( a, b, c, d, pp, normal )

      call plane_normal2exp_3d ( pp, normal, p1, p2, p3 )

      return
      end
      subroutine plane_imp2normal_3d ( a, b, c, d, pp, normal )

c***********************************************************************72
c
cc PLANE_IMP2NORMAL_3D converts an implicit plane to normal form in 3D.
c
c  Discussion:
c
c    The implicit form of a plane in 3D is
c
c      A * X + B * Y + C * Z + D = 0.
c
c    The normal form of a plane in 3D is
c
c      PP, a point on the plane, and
c      N, the unit normal to the plane.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, D, the implicit plane parameters.
c
c    Output, double precision PP(3), a point on the plane.
c
c    Output, double precision NORMAL(3), the unit normal vector to the plane.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision norm
      double precision normal(dim_num)
      double precision pp(dim_num)

      norm = sqrt ( a * a + b * b + c * c )

      if ( norm .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP2NORMAL_3D - Fatal error!'
        write ( *, '(a)' ) '  The plane (A,B,C) has zero norm.'
        stop
      end if

      normal(1) = a / norm
      normal(2) = b / norm
      normal(3) = c / norm

      if ( a .ne. 0.0D+00 ) then
        pp(1) = - d / a
        pp(2) = 0.0D+00
        pp(3) = 0.0D+00
      else if ( b .ne. 0.0D+00 ) then
        pp(1) = 0.0D+00
        pp(2) = - d / b
        pp(3) = 0.0D+00
      else if ( c .ne. 0.0D+00 ) then
        pp(1) = 0.0D+00
        pp(2) = 0.0D+00
        pp(3) = - d / c
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP2NORMAL_3D - Fatal error!'
        write ( *, '(a)' ) '  The (A,B,C) vector is null.'
        stop
      end if

      return
      end
      subroutine plane_normal_basis_3d ( pp, normal, pq, pr )

c*********************************************************************72
c
cc PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
c
c  Discussion:
c
c    The normal form of a plane in 3D is:
c
c      PP is a point on the plane,
c      N is a normal vector to the plane.
c
c    The two vectors to be computed, PQ and PR, can be regarded as
c    the basis of a Cartesian coordinate system for points in the plane.
c    Any point in the plane can be described in terms of the "origin"
c    point PP plus a weighted sum of the two vectors PQ and PR:
c
c      P = PP + a * PQ + b * PR.
c
c    The vectors PQ and PR have unit length, and are perpendicular to N
c    and to each other.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.  (Actually,
c    we never need to know these values to do the calculationc)
c
c    Input, double precision NORMAL(3), a normal vector N to the plane.  The
c    vector must not have zero length, but it is not necessary for N
c    to have unit length.
c
c    Output, double precision PQ(3), a vector of unit length,
c    perpendicular to the vector N and the vector PR.
c
c    Output, double precision PR(3), a vector of unit length,
c    perpendicular to the vector N and the vector PQ.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision r8vec_norm
      double precision normal(dim_num)
      double precision normal_norm
      double precision pp(dim_num)
      double precision pq(dim_num)
      double precision pr(dim_num)
      double precision pr_norm
c
c  Compute the length of NORMAL.
c
      normal_norm = r8vec_norm ( dim_num, normal )

      if ( normal_norm .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_NORMAL_BASIS_3D - Fatal error!'
        write ( *, '(a)' ) '  The normal vector is 0.'
        stop
      end if
c
c  Find a vector PQ that is normal to NORMAL and has unit length.
c
      call r8vec_any_normal ( dim_num, normal, pq )
c
c  Now just take the cross product NORMAL x PQ to get the PR vector.
c
      call r8vec_cross_product_3d ( normal, pq, pr )

      pr_norm = r8vec_norm ( dim_num, pr )

      do i = 1, dim_num
        pr(i) = pr(i) / pr_norm
      end do

      return
      end
      subroutine plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr,
     &  xyz )

c*********************************************************************72
c
cc PLANE_NORMAL_QR_TO_XYZ: QR_TO_XYZ coordinates for a normal form plane.
c
c  Discussion:
c
c    The normal form of a plane in 3D is:
c
c      PP is a point on the plane,
c      NORMAL is a normal vector to the plane.
c
c    Two vectors PQ and PR can be computed with the properties that
c    * NORMAL, PQ and PR are pairwise orthogonal;
c    * PQ and PR have unit length;
c    * every point P in the plane has a "QR" representation
c      as P = PP + q * PQ + r * PR.
c
c    This function is given the QR coordinates of a set of points on the
c    plane, and returns the XYZ coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.
c
c    Input, double precision NORMAL(3), a normal vector N to the plane.  The
c    vector must not have zero length, but it is not necessary for N
c    to have unit length.
c
c    Input, double precision PQ(3), a vector of unit length,
c    perpendicular to the vector N and the vector PR.
c
c    Input, double precision PR(3), a vector of unit length,
c    perpendicular to the vector N and the vector PQ.
c
c    Input, integer N, the number of points on the plane.
c
c    Input, double precision QR(2,N), the QR coordinates of the points.
c
c    Output, double precision XYZ(3,N), the XYZ coordinates of the points.
c
      implicit none

      integer n

      integer i
      integer j
      double precision normal(3)
      double precision pp(3)
      double precision pq(3)
      double precision pr(3)
      double precision qr(2,n)
      double precision xyz(3,n)

      do j = 1, n
        do i = 1, 3
          xyz(i,j) = pp(i) + pq(i) * qr(1,j) + pr(i) * qr(2,j)
        end do
      end do

      return
      end
      subroutine plane_normal_tetrahedron_intersect ( pp, normal, t,
     &  int_num, pint )

c*********************************************************************72
c
cc PLANE_NORMAL_TETRAHEDRON_INTERSECT intersects a plane and a tetrahedron.
c
c  Discussion:
c
c    The intersection of a plane and a tetrahedron is one of:
c    0) empty
c    1) a single point
c    2) a single line segment
c    3) a triangle
c    4) a quadrilateral.
c
c    In each case, the region of intersection can be described by the
c    corresponding number of points.  In particular, cases 2, 3 and 4
c    are described by the vertices that bound the line segment, triangle,
c    or quadrilateral.
c
c    The normal form of a plane is:
c
c      PP is a point on the plane,
c      N is a normal vector to the plane.
c
c    The form of a tetrahedron is
c
c      T(1:3,1:4) contains the coordinates of the vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.
c
c    Input, double precision NORMAL(3), a normal vector to the plane.
c
c    Input, double precision T(3,4), the tetrahedron vertices.
c
c    Output, integer INT_NUM, the number of intersection
c    points returned.  This will be 0, 1, 2, 3 or 4.
c
c    Output, double precision PINT(3,4), the coordinates of the
c    intersection points.
c
      implicit none

      double precision area1
      double precision area2
      double precision d(4)
      double precision dn
      double precision dpp
      integer i
      integer int_num
      integer j
      integer j1
      integer j2
      double precision normal(3)
      double precision pint(3,4)
      double precision pp(3)
      logical r8_sign_opposite_strict
      double precision r8vec_dot_product
      logical r8vec_negative_strict
      logical r8vec_positive_strict
      double precision t(3,4)
      double precision temp

      int_num = 0
      do j = 1, 4
        do i = 1, 3
          pint(i,j) = 0.0D+00
        end do
      end do
c
c  DN is the length of the normal vector.
c
      dn = sqrt ( r8vec_dot_product ( 3, normal, normal ) )
c
c  DPP is the distance between the origin and the projection of the
c  point PP onto the normal vector.
c
      dpp = dn - r8vec_dot_product ( 3, normal, pp ) / dn
c
c  D(I) is positive, zero, or negative if vertex I is above,
c  on, or below the plane.
c
      do j = 1, 4
        d(j) = dn - dpp
        do i = 1, 3
          d(j) = d(j) - normal(i) * t(i,j)
        end do
      end do
c
c  If all D are positive or negative, no intersection.
c
      if ( r8vec_negative_strict ( 4, d ) .or.
     &     r8vec_positive_strict ( 4, d ) ) then
        int_num = 0
        return
      end if
c
c  Points with zero distance are automatically added to the list.
c
c  For each point with nonzero distance, seek another point
c  with opposite sign and higher index, and compute the intersection
c  of the line between those points and the plane.
c
      do j1 = 1, 4

        if ( d(j1) .eq. 0.0D+00 ) then
          int_num = int_num + 1
          do i = 1, 3
            pint(i,int_num) = t(i,j1)
          end do
        else
          do j2 = j1 + 1, 4
            if ( r8_sign_opposite_strict ( d(j1), d(j2) ) ) then
              int_num = int_num + 1
              do i = 1, 3
                pint(i,int_num) = ( d(j1)         * t(i,j2)
     &                                    - d(j2) * t(i,j1) )
     &                          / ( d(j1) - d(j2) )
              end do
            end if
          end do
        end if
      end do
c
c  If four points were found, try to order them properly.
c
      if ( int_num .eq. 4 ) then
        call quad_area_3d ( pint, area1 )
        do i = 1, 3
          temp      = pint(i,3)
          pint(i,3) = pint(i,4)
          pint(i,4) = temp
        end do
        call quad_area_3d ( pint, area2 )
        do i = 1, 3
          temp      = pint(i,3)
          pint(i,3) = pint(i,4)
          pint(i,4) = temp
        end do
      end if

      return
      end
      subroutine plane_normal_uniform_3d ( seed, pp, normal )

c*********************************************************************72
c
cc PLANE_NORMAL_UNIFORM_3D generates a random normal plane in 3D.
c
c  Discussion:
c
c    The normal form of a plane is:
c
c      PP is a point on the plane,
c      N is a normal vector to the plane.
c
c    The point PP will be chosen at random inside the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer  SEED, a seed for the random
c    number generator.
c
c    Output, double precision PP(3), a point on the plane.
c
c    Output, double precision NORMAL(3), the unit normal vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision norm
      double precision normal(dim_num)
      double precision pp(dim_num)
      integer seed
c
c  Pick PP as a random point inside the unit sphere in ND.
c
      call ball_unit_sample_3d ( seed, pp )
c
c  Get values from a standard normal distribution.
c
      call r8vec_normal_01 ( dim_num, seed, normal )
c
c  Compute the length of the vector.
c
      norm = 0.0D+00
      do i = 1, dim_num
        norm = norm + normal(i)**2
      end do
      norm = sqrt ( norm )
c
c  Normalize the vector.
c
      do i = 1, dim_num
        normal(i) = normal(i) / norm
      end do

      return
      end
      subroutine plane_normal_uniform_nd ( dim_num, seed, pp, normal )

c*********************************************************************72
c
cc PLANE_NORMAL_UNIFORM_ND generates a random normal plane in ND.
c
c  Discussion:
c
c    The normal form of a plane is:
c
c      PP is a point on the plane,
c      N is a normal vector to the plane.
c
c    The point PP will be chosen at random inside the unit sphere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input/output, integer  SEED, a seed for the random
c    number generator.
c
c    Output, double precision PP(DIM_NUM), a point on the plane.
c
c    Output, double precision NORMAL(DIM_NUM), the unit normal vector.
c
      implicit none

      integer dim_num

      integer i
      double precision norm
      double precision normal(dim_num)
      double precision pp(dim_num)
      integer seed
c
c  Pick PP as a random point inside the unit sphere in ND.
c
      call ball_unit_sample_nd ( dim_num, seed, pp )
c
c  Get values from a standard normal distribution.
c
      call r8vec_normal_01 ( dim_num, seed, normal )
c
c  Compute the length of the vector.
c
      norm = 0.0D+00
      do i = 1, dim_num
        norm = norm + normal(i)**2
      end do
      norm = sqrt ( norm )
c
c  Normalize the vector.
c
      do i = 1, dim_num
        normal(i) = normal(i) / norm
      end do

      return
      end
      subroutine plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz,
     &  qr )

c*********************************************************************72
c
cc PLANE_NORMAL_XYZ_TO_QR: XYZ to QR coordinates for a normal form plane.
c
c  Discussion:
c
c    The normal form of a plane in 3D is:
c
c      PP is a point on the plane,
c      NORMAL is a normal vector to the plane.
c
c    Two vectors PQ and PR can be computed with the properties that
c    * NORMAL, PQ and PR are pairwise orthogonal;
c    * PQ and PR have unit length;
c    * every point P in the plane has a "QR" representation
c      as P = PP + q * PQ + r * PR.
c
c    This function is given the XYZ coordinates of a set of points on the
c    plane, and returns the QR coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.
c
c    Input, double precision NORMAL(3), a normal vector N to the plane.  The
c    vector must not have zero length, but it is not necessary for N
c    to have unit length.
c
c    Input, double precision PQ(3), a vector of unit length,
c    perpendicular to the vector N and the vector PR.
c
c    Input, double precision PR(3), a vector of unit length,
c    perpendicular to the vector N and the vector PQ.
c
c    Input, integer N, the number of points on the plane.
c
c    Input, double precision XYZ(3,N), the XYZ coordinates of the points.
c
c    Output, double precision QR(2,N), the QR coordinates of the points.
c
      implicit none

      integer n

      integer j
      double precision normal(3)
      double precision pp(3)
      double precision pq(3)
      double precision pr(3)
      double precision qr(2,n)
      double precision xyz(3,n)

      do j = 1, n
        qr(1,j) = pq(1) * ( xyz(1,j) - pp(1) )
     &          + pq(2) * ( xyz(2,j) - pp(2) )
     &          + pq(3) * ( xyz(3,j) - pp(3) )
        qr(2,j) = pr(1) * ( xyz(1,j) - pp(1) )
     &          + pr(2) * ( xyz(2,j) - pp(2) )
     &          + pr(3) * ( xyz(3,j) - pp(3) )
      end do

      return
      end
      subroutine plane_normal2exp_3d ( pp, normal, p1, p2, p3 )

c*********************************************************************72
c
cc PLANE_NORMAL2EXP_3D converts a normal plane to explicit form in 3D.
c
c  Discussion:
c
c    The normal form of a plane in 3D is
c
c      PP, a point on the plane, and
c      N, the unit normal to the plane.
c
c    The explicit form of a plane in 3D is
c
c      the plane through P1, P2 and P3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.
c
c    Input, double precision NORMAL(3), a normal vector N to the plane.  The
c    vector must not have zero length, but it is not necessary for N
c    to have unit length.
c
c    Output, double precision P1(3), P2(3), P3(3), three points on the plane.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision normal(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pp(dim_num)
      double precision pq(dim_num)
      double precision pr(dim_num)

      call plane_normal_basis_3d ( pp, normal, pq, pr )

      do i = 1, dim_num
        p1(i) = pp(i)
        p2(i) = pp(i) + pq(i)
        p3(i) = pp(i) + pr(i)
      end do

      return
      end
      subroutine plane_normal2imp_3d ( pp, normal, a, b, c, d )

c*********************************************************************72
c
cc PLANE_NORMAL2IMP_3D converts a normal form plane to implicit form in 3D.
c
c  Discussion:
c
c    The normal form of a plane in 3D is
c
c      PP, a point on the plane, and
c      N, the unit normal to the plane.
c
c    The implicit form of a plane in 3D is
c
c      A * X + B * Y + C * Z + D = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PP(3), a point on the plane.
c
c    Input, double precision NORMAL(3), the unit normal vector to the plane.
c
c    Output, double precision A, B, C, D, the implicit plane parameters.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a
      double precision b
      double precision c
      double precision d
      double precision normal(dim_num)
      double precision pp(dim_num)

      a = normal(1)
      b = normal(2)
      c = normal(3)
      d = - a * pp(1) - b * pp(2) - c * pp(3)

      return
      end
      subroutine points_dist_nd ( dim_num, p1, p2, dist )

c*********************************************************************72
c
cc POINTS_DIST_ND finds the distance between two points in ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision P1(DIM_NUM), P2(DIM_NUM), the coordinates 
c    of two points.
c
c    Output, double precision DIST, the distance between the points.
c
      implicit none

      integer dim_num

      double precision dist
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r8vec_diff_norm

      dist = r8vec_diff_norm ( dim_num, p1, p2 )

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
c  Remember the starting point, so we know when to stopc
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
      subroutine points_point_near_naive_nd ( dim_num, set_num, pset,
     &  p, i_min, dist_min )

c*********************************************************************72
c
cc POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
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
c    03 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer SET_NUM, the number of points in the set.
c
c    Input, double precision PSET(DIM_NUM,SET_NUM), the points in the set.
c
c    Input, double precision P(DIM_NUM), the point whose nearest neighbor
c    is sought.
c
c    Output, integer I_MIN, the index of the nearest point in
c    PSET to P.
c
c    Output, double precision DIST_MIN, the distance between P(*)
c    and PSET(*,I_MIN).
c
      implicit none

      integer dim_num
      integer set_num

      double precision d
      double precision dist_min
      integer i
      integer i_min
      integer j
      double precision p(dim_num)
      double precision pset(dim_num,set_num)
      double precision r8_huge

      dist_min = r8_huge ( )
      i_min = -1

      do i = 1, set_num

        d = 0.0D+00
        do j = 1, dim_num
          d = d + ( p(j) - pset(j,i) )**2
        end do

        if ( d .lt. dist_min ) then
          dist_min = d
          i_min = i
        end if

      end do

      dist_min = sqrt ( dist_min )

      return
      end
      subroutine polar_to_xy ( r, t, xy )

c*********************************************************************72
c
cc POLAR_TO_XY converts polar coordinates to XY coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, T, the radius and angle (in radians).
c
c    Output, double precision XY(2), the Cartesian coordinates.
c
      implicit none

      double precision r
      double precision t
      double precision xy(2)

      xy(1) = r * cos ( t )
      xy(2) = r * sin ( t )

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
      subroutine polygon_angles_2d ( n, v, angle )

c*********************************************************************72
c
cc POLYGON_ANGLES_2D computes the interior angles of a polygon in 2D.
c
c  Discussion:
c
c    The vertices should be listed in counter clockwise order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision V(2,N), the vertices.
c
c    Output, double precision ANGLE(N), the angles of the polygon,
c    in radians.
c
      implicit none

      integer n
 
      double precision angle(n)
      double precision angle_rad_2d
      integer i
      integer i4_wrap
      integer im1
      integer ip1
      double precision v(2,n)

      if ( n .le. 2 ) then
        do i = 1, n
          angle(i) = 0.0D+00
        end do
        return
      end if
 
      do i = 1, n

        im1 = i4_wrap ( i - 1, 1, n )
        ip1 = i4_wrap ( i + 1, 1, n )

        angle(i) = angle_rad_2d ( v(1,im1), v(1,i), v(1,ip1) )

      end do

      return
      end
      subroutine polygon_area_2d ( n, v, area )

c*********************************************************************72
c
cc POLYGON_AREA_2D computes the area of a polygon in 2D.
c
c  Discussion:
c
c    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
c    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
c
c    If the vertices are given in counter clockwise order, the area
c    will be positive.  If the vertices are given in clockwise order,
c    the area will be negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision V(2,N), the vertices.
c
c    Output, double precision AREA, the absolute area of the polygon.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      integer i
      integer i4_wrap
      integer im1
      integer ip1
      double precision v(dim_num,n)

      area = 0.0D+00

      do i = 1, n

        im1 = i4_wrap ( i-1, 1, n )
        ip1 = i4_wrap ( i+1, 1, n )

        area = area + v(1,i) * ( v(2,ip1) - v(2,im1) )

      end do

      area = 0.5D+00 * area

      return
      end
      subroutine polygon_area_2d_2 ( n, v, area )

c*********************************************************************72
c
cc POLYGON_AREA_2D_2 computes the area of a polygon in 2D.
c
c  Discussion:
c
c    The area is the sum of the areas of the triangles formed by
c    node N with consecutive pairs of nodes.
c
c    If the vertices are given in counter clockwise order, the area
c    will be positive.  If the vertices are given in clockwise order,
c    the area will be negative.
c
c    Thanks to Martin Pineault for noticing that an earlier version
c    of this routine would not correctly compute the area of a nonconvex
c    polygon.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision V(2,N), the vertices.
c
c    Output, double precision AREA, the absolute area of the polygon.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision area_triangle
      integer j
      double precision t(dim_num,3)
      double precision v(dim_num,n)

      area = 0.0D+00

      do j = 1, n - 2

        t(1,1) = v(1,j)
        t(2,1) = v(2,j)
        t(1,2) = v(1,j+1)
        t(2,2) = v(2,j+1)
        t(1,3) = v(1,n)
        t(2,3) = v(2,n)

        call triangle_area_2d ( t, area_triangle )

        area = area + area_triangle

      end do

      return
      end
      subroutine polygon_area_3d ( n, v, area, normal )

c*********************************************************************72
c
cc POLYGON_AREA_3D computes the area of a polygon in 3D.
c
c  Discussion:
c
c    The computation is not valid unless the vertices of the polygon
c    lie in a plane, so that the polygon that is defined is "flat".
c
c    The polygon does not have to be "regular", that is, neither its
c    sides nor its angles need to be equal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Allen Van Gelder,
c    Efficient Computation of Polygon Area and Polyhedron Volume,
c    Graphics Gems V, 
c    edited by Alan Paeth,
c    AP Professional, 1995, T385.G6975.
c
c  Parameters:
c
c    Input, integer N, the number of vertices.
c
c    Input, double precision V(3,N), the coordinates of the vertices.
c    The vertices should be listed in neighboring order.
c
c    Output, double precision AREA, the area of the polygon.
c
c    Output, double precision NORMAL(3), the unit normal vector to the polygon.
c
      implicit none

      integer n

      integer dim_num
      parameter ( dim_num = 3 )

      double precision area
      double precision cross(dim_num)
      integer i
      integer j
      integer jp1
      double precision normal(dim_num)
      double precision r8vec_norm
      double precision v(dim_num,n)

      do i = 1, dim_num
        normal(i) = 0.0D+00
      end do

      do j = 1, n

        if ( j .lt. n ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if
c
c  Compute the cross product vector.
c
        cross(1) = v(2,j) * v(3,jp1) - v(3,j) * v(2,jp1)
        cross(2) = v(3,j) * v(1,jp1) - v(1,j) * v(3,jp1)
        cross(3) = v(1,j) * v(2,jp1) - v(2,j) * v(1,jp1)

        do i = 1, dim_num
          normal(i) = normal(i) + cross(i)
        end do

      end do

      area = r8vec_norm ( dim_num, normal )

      if ( area .ne. 0.0D+00 ) then
        do i = 1, dim_num
          normal(i) = normal(i) / area
        end do
      else
        do i = 1, dim_num
          normal(i) = 1.0D+00 / sqrt ( dble ( dim_num ) )
        end do
      end if

      area = 0.5D+00 * area

      return
      end
      subroutine polygon_area_3d_2 ( n, v, area )

c*********************************************************************72
c
cc POLYGON_AREA_3D_2 computes the area of a polygon in 3D.
c
c  Discussion:
c
c    The computation is not valid unless the vertices of the polygon
c    lie in a plane, so that the polygon that is defined is "flat".
c
c    The polygon does not have to be "regular", that is, neither its
c    sides nor its angles need to be equal.
c
c    The area is computed as the sum of the areas of the triangles 
c    formed by the last node with consecutive pairs of nodes (1,2),
c    (2,3), ..., and (N-2,N-1).
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
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision V(3,N), the coordinates of the vertices.
c
c    Output, double precision AREA, the area of the polygon.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 3 )

      double precision area
      double precision area_vector(dim_num)
      double precision area_vector_triangle(dim_num)
      integer i
      integer j
      double precision r8vec_norm
      double precision t(dim_num,3)
      double precision v(dim_num,n)

      do i = 1, dim_num
        area_vector(i) = 0.0D+00
      end do

      do j = 1, n - 2

        do i = 1, dim_num
          t(i,1) = v(i,j)
          t(i,2) = v(i,j+1)
          t(i,3) = v(i,n)
        end do

        call triangle_area_vector_3d ( t, area_vector_triangle )

        do i = 1, dim_num
          area_vector(i) = area_vector(i) + area_vector_triangle(i)
        end do

      end do

      area = 0.5D+00 * r8vec_norm ( dim_num, area_vector )

      return
      end
      subroutine polygon_inrad_data_2d ( n, radin, area, radout, side )

c*********************************************************************72
c
cc POLYGON_INRAD_DATA_2D determines polygonal data from its inner radius in 2D.
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
c  Parameters:
c
c    Input, integer N, the number of sides of the polygon.
c    N must be at least 3.
c
c    Input, double precision RADIN, the inner radius of the polygon, that is,
c    the radius of the largest circle that can be inscribed within
c    the polygon.
c
c    Output, double precision AREA, the area of the regular polygon.
c
c    Output, double precision RADOUT, the outer radius of the polygon, that is,
c    the radius of the smallest circle that can be described about
c    the polygon.
c
c    Output, double precision SIDE, the length of one side of the polygon.
c
      implicit none

      double precision angle
      double precision area
      integer n
      double precision pi 
      parameter ( pi = 3.141592653589793D+00 )
      double precision radin
      double precision radout
      double precision side

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_INRAD_DATA_2D - Fatal error!'
        write ( *, '(a)' ) '  Input value of N must be at least 3'
        write ( *, '(a,i8)' ) '  but your input value was N = ', n
        stop
      end if

      angle = pi / dble ( n )
      area = dble ( n ) * radin * radin * tan ( angle )
      side = 2.0D+00 * radin * tan ( angle )
      radout = 0.5D+00 * side / sin ( angle )

      return
      end
      subroutine polygon_lattice_area_2d ( i, b, area )

c*********************************************************************72
c
cc POLYGON_LATTICE_AREA_2D computes the area of a lattice polygon in 2D.
c
c  Discussion:
c
c    We define a lattice to be the 2D plane, in which the points
c    whose (X,Y) coordinates are both integers are given a special
c    status as "lattice points".
c
c    A lattice polygon is a polygon whose vertices are lattice points.
c
c    The area of a lattice polygon can be computed by Pick's Theorem:
c
c      Area = I + B / 2 - 1
c
c    where
c
c      I = the number of lattice points contained strictly inside the polygon;
c
c      B = the number of lattice points that lie exactly on the boundary.
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
c    Branko Gruenbaum, Geoffrey Shephard,
c    Pick's Theorem,
c    The American Mathematical Monthly,
c    Volume 100, Number 2, February 1993, pages 150-161.
c
c  Parameters:
c
c    Input, integer I, the number of interior lattice points.
c
c    Input, integer B, the number of boundary lattice points.
c
c    Output, double precision AREA, the area of the lattice polygon.
c
      implicit none

      double precision area
      integer b
      integer i

      area = dble ( i ) + dble ( b ) / 2.0D+00 - 1.0D+00

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
      subroutine polyhedron_area_3d ( coord, order_max, face_num, node, 
     &  node_num, order, area )

c*********************************************************************72
c
cc POLYHEDRON_AREA_3D computes the surface area of a polyhedron in 3D.
c
c  Discussion:
c
c    The computation is not valid unless the faces of the polyhedron
c    are planar polygons.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Allen Van Gelder,
c    Efficient Computation of Polygon Area and Polyhedron Volume,
c    in Graphics Gems V, 
c    edited by Alan Paeth,
c    AP Professional, 1995, T385.G6975
c
c  Parameters:
c
c    Input, double precision COORD(3,NODE_NUM), the coordinates of the
c    vertices.  The vertices may be listed in any order.
c
c    Input, integer ORDER_MAX, the maximum number of vertices 
c    that make up a face of the polyhedron.
c
c    Input, integer FACE_NUM, the number of faces of the 
c    polyhedron.
c
c    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined 
c    by the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
c    are listed in neighboring order.
c
c    Input, integer NODE_NUM, the number of points stored in COORD.
c
c    Input, integer ORDER(FACE_NUM), the number of vertices 
c    making up each face.
c
c    Output, double precision AREA, the total surface area of the polyhedron.
c
      implicit none

      integer face_num
      integer order_max
      integer dim_num
      parameter ( dim_num = 3 )
      integer node_num

      double precision ainc
      double precision area
      double precision coord(dim_num,node_num)
      integer face
      integer i
      integer j
      integer k1
      integer k2
      integer node(face_num,order_max)
      integer order(face_num)
      double precision r8vec_norm
      double precision v(dim_num)

      area = 0.0D+00
c
c  For each face
c
      do face = 1, face_num

        do i = 1, dim_num
          v(i) = 0.0D+00
        end do
c
c  For each triangle in the face, compute the normal vector.
c
        do j = 1, order(face)

          k1 = node(face,j)

          if ( j .lt. order(face) ) then
            k2 = node(face,j+1)
          else
            k2 = node(face,1)
          end if
c
c  Compute the cross product.
c
          v(1) = v(1) + coord(2,k1) * coord(3,k2) 
     &                - coord(3,k1) * coord(2,k2)
          v(2) = v(2) + coord(3,k1) * coord(1,k2) 
     &                - coord(1,k1) * coord(3,k2)
          v(3) = v(3) + coord(1,k1) * coord(2,k2) 
     &                - coord(2,k1) * coord(1,k2)

        end do
c
c  Add the magnitude of the normal vector to the sum.
c
        ainc = r8vec_norm ( dim_num, v )
        area = area + ainc

      end do

      area = 0.5D+00 * area

      return
      end
      subroutine polyhedron_volume_3d ( coord, order_max, face_num,
     &  node, node_num, order, volume )

c*********************************************************************72
c
cc POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision COORD(3,NODE_NUM), the coordinates of
c    the vertices.  The vertices may be listed in any order.
c
c    Input, integer ORDER_MAX, the maximum number of vertices
c    that make up a face of the polyhedron.
c
c    Input, integer FACE_NUM, the number of faces of the
c    polyhedron.
c
c    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined by
c    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
c    are listed in neighboring order.
c
c    Input, integer NODE_NUM, the number of points stored in COORD.
c
c    Input, integer ORDER(FACE_NUM), the number of vertices making
c    up each face.
c
c    Output, double precision VOLUME, the volume of the polyhedron.
c
      implicit none

      integer face_num
      integer order_max
      integer dim_num
      parameter ( dim_num = 3 )
      integer node_num

      double precision coord(dim_num,node_num)
      integer face
      integer n1
      integer n2
      integer n3
      integer node(face_num,order_max)
      integer order(face_num)
      integer v
      double precision volume

      volume = 0.0D+00
c
c  Triangulate each face.
c
      do face = 1, face_num

        n3 = node(face,order(face))

        do v = 1, order(face) - 2

          n1 = node(face,v)
          n2 = node(face,v+1)

          volume = volume
     &      + coord(1,n1)
     &      * ( coord(2,n2) * coord(3,n3) - coord(2,n3) * coord(3,n2) )
     &      + coord(1,n2)
     &      * ( coord(2,n3) * coord(3,n1) - coord(2,n1) * coord(3,n3) )
     &      + coord(1,n3)
     &      * ( coord(2,n1) * coord(3,n2) - coord(2,n2) * coord(3,n1) )

        end do

      end do

      volume = volume / 6.0D+00

      return
      end
      subroutine pyramid_volume_3d ( h, s, volume )

c*********************************************************************72
c
cc PYRAMID_VOLUME_3D computes the volume of a pyramid with square base in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision H, S, the height of the pyramid, and the
c    length of one side of the square base.
c
c    Output, double precision VOLUME, the volume of the pyramid.
c
      implicit none

      double precision h
      double precision s
      double precision volume

      volume = s * s * h / 3.0D+00

      return
      end
      subroutine quad_area_2d ( q, area )

c*********************************************************************72
c
cc QUAD_AREA_2D computes the area of a quadrilateral in 2D.
c
c  Discussion:
c
c    This algorithm should be able to handle nonconvex quadrilaterals.
c
c    The vertices of the quadrilateral should be listed in counter clockwise
c    order, so that the area is positive.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(2,4), the vertices of
c    the quadrilateral.  The corners should be specified in clockwise
c    or counter clockwise order.
c
c    Output, double precision AREA, the area of the quadrilateral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision area_triangle
      double precision q(dim_num,4)
      double precision t(dim_num,3)

      area = 0.0D+00

      t(1,1) = q(1,1)
      t(2,1) = q(2,1)

      t(1,2) = q(1,2)
      t(2,2) = q(2,2)

      t(1,3) = q(1,3)
      t(2,3) = q(2,3)

      call triangle_area_2d ( t, area_triangle )

      area = area + area_triangle

      t(1,1) = q(1,3)
      t(2,1) = q(2,3)

      t(1,2) = q(1,4)
      t(2,2) = q(2,4)

      t(1,3) = q(1,1)
      t(2,3) = q(2,1)

      call triangle_area_2d ( t, area_triangle )

      area = area + area_triangle

      return
      end
      subroutine quad_area2_2d ( q, area )

c*********************************************************************72
c
cc QUAD_AREA2_2D computes the area of a quadrilateral in 2D.
c
c  Discussion:
c
c    A quadrilateral is a polygon defined by 4 vertices.
c
c    This algorithm computes the area of the related
c    Varignon parallelogram first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(2,4), the vertices, specified in
c    counter clockwise order.
c
c    Output, double precision AREA, the area of the quadrilateral.
c
      implicit none

      double precision area
      integer i
      integer j
      double precision p(2,4)
      double precision q(2,4)
c
c  Define a parallelogram by averaging consecutive vertices.
c
      do j = 1, 3
        do i = 1, 2
          p(i,j) = ( q(i,j) + q(i,j+1) ) / 2.0D+00
        end do
      end do

      do i = 1, 2
        p(i,4) = ( q(i,4) + q(i,1) ) / 2.0D+00
      end do
c
c  Compute the area.
c
      call parallelogram_area_2d ( p, area )
c
c  The quadrilateral's area is twice that of the parallelogram.
c
      area = 2.0D+00 * area

      return
      end
      subroutine quad_area_3d ( q, area )

c*********************************************************************72
c
cc QUAD_AREA_3D computes the area of a quadrilateral in 3D.
c
c  Discussion:
c
c    A quadrilateral is a polygon defined by 4 vertices.
c
c    It is assumed that the four vertices of the quadrilateral
c    are coplanar.
c
c    This algorithm computes the area of the related
c    Varignon parallelogram first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(3,4), the vertices, specified in
c    counter clockwise order.
c
c    Output, double precision AREA, the area of the quadrilateral.
c
      implicit none

      double precision area
      integer i
      integer j
      double precision p(3,4)
      double precision q(3,4)
c
c  Define a parallelogram by averaging consecutive vertices.
c
      do j = 1, 3
        do i = 1, 3
          p(i,j) = ( q(i,j) + q(i,j+1) ) / 2.0D+00
        end do
      end do

      do i = 1, 3
        p(i,4) = ( q(i,4) + q(i,1) ) / 2.0D+00
      end do
c
c  Compute the area.
c
      call parallelogram_area_3d ( p, area )
c
c  The quadrilateral's area is twice that of the parallelogram.
c
      area = 2.0D+00 * area

      return
      end
      subroutine quad_contains_point_2d ( q, p, inside )

c*********************************************************************72
c
cc QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
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
c    Input, double precision Q(2,4), the vertices of the quadrilateral.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is in the quadrilateral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle_1
      double precision angle_2
      double precision angle_rad_2d
      logical inside
      double precision p(dim_num)
      double precision q(dim_num,4)
c
c  This will only handle convex quadrilaterals.
c
      inside = .false.

      angle_1 = angle_rad_2d ( q(1:2,1), q(1:2,2), q(1:2,3) )
      angle_2 = angle_rad_2d ( q(1:2,1), q(1:2,2), p(1:2) )

      if ( angle_1 .lt. angle_2 ) then
        return
      end if

      angle_1 = angle_rad_2d ( q(1:2,2), q(1:2,3), q(1:2,4) )
      angle_2 = angle_rad_2d ( q(1:2,2), q(1:2,3), p(1:2) )

      if ( angle_1 .lt. angle_2 ) then
        return
      end if

      angle_1 = angle_rad_2d ( q(1:2,3), q(1:2,4), q(1:2,1) )
      angle_2 = angle_rad_2d ( q(1:2,3), q(1:2,4), p(1:2) )

      if ( angle_1 .lt. angle_2 ) then
        return
      end if

      angle_1 = angle_rad_2d ( q(1:2,4), q(1:2,1), q(1:2,2) )
      angle_2 = angle_rad_2d ( q(1:2,4), q(1:2,1), p(1:2) )

      if ( angle_1 .lt. angle_2 ) then
        return
      end if

      inside = .true.

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
      subroutine quad_point_dist_2d ( q, p, dist )

c*********************************************************************72
c
cc QUAD_POINT_DIST_2D: distance ( quadrilateral, point ) in 2D.
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
c    Input, double precision Q(2,4), the quadrilateral vertices.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision DIST, the distance from the point to the
c    quadrilateral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer side_num
      parameter ( side_num = 4 )

      double precision dist
      double precision dist2
      integer i4_wrap
      integer j
      integer jp1
      double precision p(dim_num)
      double precision q(dim_num,side_num)
      double precision r8_huge
c
c  Find the distance to each of the line segments.
c
      dist = r8_huge ( )

      do j = 1, side_num

        jp1 = i4_wrap ( j+1, 1, side_num )

        call segment_point_dist_2d ( q(1,j), q(1,jp1), p, dist2 )

        if ( dist2 .lt. dist ) then
          dist = dist2
        end if

      end do

      return
      end
      subroutine quad_point_near_2d ( q, p, pn, dist )

c*********************************************************************72
c
cc QUAD_POINT_NEAR_2D computes the nearest point on a quadrilateral in 2D.
c
c  Discussion:
c
c    A quadrilateral is a polygon defined by 4 vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(2,4), the quadrilateral vertices.
c
c    Input, double precision P(2), the point whose nearest quadrilateral point
c    is to be determined.
c
c    Output, double precision PN(2), the nearest point to P.
c
c    Output, double precision DIST, the distance from the point to the
c    quadrilateral.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer side_num
      parameter ( side_num = 4 )

      double precision dist
      double precision dist2
      integer i
      integer i4_wrap
      integer j
      integer jp1
      double precision p(dim_num)
      double precision pn(dim_num)
      double precision pn2(dim_num)
      double precision q(dim_num,side_num)
      double precision r8_huge
      double precision tval
c
c  Find the distance to each of the line segments that make up the edges
c  of the quadrilateral.
c
      dist = r8_huge ( )

      do i = 1, dim_num
        pn(i) = 0.0D+00
      end do

      do j = 1, side_num

        jp1 = i4_wrap ( j+1, 1, side_num )

        call segment_point_near_2d ( q(1,j), q(1,jp1), p, 
     &    pn2, dist2, tval )

        if ( dist2 .lt. dist ) then
          dist = dist2
          do i = 1, dim_num
            pn(i) = pn2(i)
          end do
        end if

      end do

      return
      end
      subroutine quat_conj ( q, q2 )

c*********************************************************************72
c
cc QUAT_CONJ conjugates a quaternion.
c
c  Discussion:
c
c    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
c    may be written as
c
c      Q = A + Bi + Cj + Dk.
c
c    The conjugate of Q is
c
c      conj ( Q ) = A - Bi - Cj - Dk.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(4), the quaternion to be conjugated.
c
c    Output, double precision Q2(4), the conjugated quaternion.
c
      implicit none

      double precision q(4)
      double precision q2(4)

      q2(1) =   q(1)
      q2(2) = - q(2)
      q2(3) = - q(3)
      q2(4) = - q(4)

      return
      end
      subroutine quat_inv ( q, q2 )

c*********************************************************************72
c
cc QUAT_INV inverts a quaternion.
c
c  Discussion:
c
c    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
c    may be written as
c
c      Q = A + Bi + Cj + Dk.
c
c    The inverse of Q is
c
c      inverse ( Q ) = conjugate ( Q ) / ( norm ( Q ) )**2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(4), the quaternion to be inverted.
c
c    Output, double precision Q2(4), the inverse of the input quaternion.
c
      implicit none

      double precision bot
      integer i
      double precision q(4)
      double precision q2(4)

      bot = 0.0D+00
      do i = 1, 4
        bot = bot + q(i)**2
      end do

      q2(1) = q(1) / bot
      do i = 2, 4
        q2(1:4) = - q(1:4) / bot
      end do

      return
      end
      subroutine quat_mul ( q1, q2, q3 )

c*********************************************************************72
c
cc QUAT_MUL multiplies two quaternions.
c
c  Discussion:
c
c    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
c    may be written as
c
c      Q = A + Bi + Cj + Dk.
c
c    To multiply two quaternions, use the relationships:
c
c      i * j = -j * i = k
c      j * k = -k * j = i
c      k * i = -i * k = j
c      i * i =  j * j = k * k = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q1(4), Q2(4), the quaternions to be multiplied.
c
c    Output, double precision Q3(4), the product of the two quaternions.
c
      implicit none

      double precision q1(4)
      double precision q2(4)
      double precision q3(4)

      q3(1) = q1(1) * q2(1) - q1(2) * q2(2)
     &      - q1(3) * q2(3) - q1(4) * q2(4)
      q3(2) = q1(1) * q2(2) + q1(2) * q2(1)
     &      + q1(3) * q2(4) - q1(4) * q2(3)
      q3(3) = q1(1) * q2(3) - q1(2) * q2(4)
     &      + q1(3) * q2(1) + q1(4) * q2(2)
      q3(4) = q1(1) * q2(4) + q1(2) * q2(3)
     &      - q1(3) * q2(2) + q1(4) * q2(1)

      return
      end
      function quat_norm ( q )

c*********************************************************************72
c
cc QUAT_NORM computes the norm of a quaternion.
c
c  Discussion:
c
c    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
c    may be written as
c
c      Q = A + Bi + Cj + Dk.
c
c    The norm of Q is
c
c      norm ( Q ) = sqrt ( A * A + B * B + C * C + D * D ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(4), the quaternion.
c
c    Output, double precision QUAT_NORM, the norm of the quaternion.
c
      implicit none

      integer i
      double precision q(4)
      double precision quat_norm

      quat_norm = 0.0D+00
      do i = 1, 4
        quat_norm = quat_norm + q(i)**2
      end do
      quat_norm = sqrt ( quat_norm )

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r4_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
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
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r4_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
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
c  it generally cannot be represented exactly as a 32 bit real numberc
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

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
      function r8_asin ( s )

c*********************************************************************72
c
cc R8_ASIN computes the arc sine function, with argument truncation.
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
c    Output, double precision R8_ASIN, an angle whose sine is S.
c
      implicit none

      double precision r8_asin
      double precision s
      double precision s2

      s2 = s
      s2 = max ( s2, -1.0D+00 )
      s2 = min ( s2, +1.0D+00 )

      r8_asin = asin ( s2 )

      return
      end
      function r8_atan ( y, x )

c*********************************************************************72
c
cc R8_ATAN computes the inverse tangent of the ratio Y / X.
c
c  Discussion:
c
c    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
c    the built in functions ATAN and ATAN2 already do.
c
c    However:
c
c    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
c      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
c      and [-PI,+PI] respectively;
c
c    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
c     function by contrast always returns an angle in the first or fourth
c     quadrants.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2008
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
c    Output, double precision R8_ATAN, an angle between 0 and 2 * PI, whose
c    tangent is (Y/X), and which lies in the appropriate quadrant so that
c    the signs of its cosine and sine match those of X and Y.
c
      implicit none

      double precision abs_x
      double precision abs_y
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_atan
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

        abs_y = dabs ( y )
        abs_x = dabs ( x )

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

      r8_atan = theta

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
      function r8_is_int ( r )

c*********************************************************************72
c
cc R8_IS_INT determines if an R8 represents an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the number to be checked.
c
c    Output, logical R8_IS_INT, is TRUE if R is an integer value.
c
      implicit none

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      double precision r
      logical r8_is_int

      if ( dble ( i4_huge ) .lt. r ) then
        r8_is_int = .false.
      else if ( r .lt. - dble ( i4_huge ) ) then
        r8_is_int = .false.
      else if ( r .eq. dble ( int ( r ) ) ) then
        r8_is_int = .true.
      else
        r8_is_int = .false.
      end if

      return
      end
      function r8_modp ( x, y )

c*********************************************************************72
c
cc R8_MODP returns the nonnegative remainder of R8 division.
c
c  Formula:
c
c    If
c      REM = R8_MODP ( X, Y )
c      RMULT = ( X - REM ) / Y
c    then
c      X = Y * RMULT + REM
c    where REM is always nonnegative.
c
c  Discussion:
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360.0) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
c
c  Example:
c
c        X         Y     MOD R8_MODP  R8_MODP Factorization
c
c      107        50       7       7    107 =  2 *  50 + 7
c      107       -50       7       7    107 = -2 * -50 + 7
c     -107        50      -7      43   -107 = -3 *  50 + 43
c     -107       -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number to be divided.
c
c    Input, double precision Y, the number that divides X.
c
c    Output, double precision R8_MODP, the nonnegative remainder
c    when X is divided by Y.
c
      implicit none

      double precision r8_modp
      double precision x
      double precision y

      if ( y .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MODP - Fatal error!'
        write ( *, '(a,g14.6)' )
     &    '  R8_MODP ( X, Y ) called with Y = ', y
        stop
      end if

      r8_modp = mod ( x, y )

      if ( r8_modp .lt. 0.0D+00 ) then
        r8_modp = r8_modp + abs ( y )
      end if

      return
      end
      function r8_normal_01 ( seed )

c*********************************************************************72
c
cc R8_NORMAL_01 returns a unit pseudonormal R8.
c
c  Discussion:
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, the code can use the second
c    value that it calculated.
c
c    However, if the user has changed the SEED value between calls,
c    the routine automatically resets itself and discards the saved data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_uniform_01
      integer seed
      integer seed1
      integer seed2
      integer seed3
      integer used
      double precision v1
      double precision v2

      save seed1
      save seed2
      save seed3
      save used
      save v2

      data seed2 / 0 /
      data used / 0 /
      data v2 / 0.0D+00 /
c
c  If USED is odd, but the input SEED does not match
c  the output SEED on the previous call, then the user has changed
c  the seed.  Wipe out internal memory.
c
      if ( mod ( used, 2 ) .eq. 1 ) then

        if ( seed .ne. seed2 ) then
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        end if

      end if
c
c  If USED is even, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        seed1 = seed

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed

        r2 = r8_uniform_01 ( seed )

        seed3 = seed

        v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

        r8_normal_01 = v1
        seed = seed2
c
c  If USED is odd (and the input SEED matched the output value from
c  the previous call), return the second normal and its corresponding seed.
c
      else

        r8_normal_01 = v2
        seed = seed3

      end if

      used = used + 1

      return
      end
      function r8_pi ( )

c*********************************************************************72
c
cc R8_PI returns the value of pi as an R8.
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
c    Output, double precision R8_PI, the value of pi.
c
      implicit none

      double precision r8_pi

      r8_pi = 3.141592653589793D+00

      return
      end
      function r8_sign_opposite_strict ( r1, r2 )

c*********************************************************************72
c
cc R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.
c
c  Discussion:
c
c    This test could be coded numerically as
c
c      if ( r1 * r2 < 0.0 ) then ...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the values to check.
c
c    Output, logical R8_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
c    or ( R2 < 0 and 0 < R1 ).
c
      implicit none

      double precision r1
      double precision r2
      logical r8_sign_opposite_strict

      r8_sign_opposite_strict = 
     &  ( r1 .lt. 0.0D+00 .and. 0.0D+00 .lt. r2 ) .or.
     &  ( r2 .lt. 0.0D+00 .and. 0.0D+00 .lt. r1 )

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
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real numberc
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
c  it generally cannot be represented exactly as a 32 bit real numberc
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
      integer dim_num
      parameter ( dim_num = 2 )

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

      call perm_check ( n, p, base, ierror )

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
      subroutine r8ge_sl ( n, a_lu, pivot, b, job )

c*********************************************************************72
c
cc R8GE_SL solves a system factored by R8GE_FA.
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
c    R8GE_SL is a simplified version of the LINPACK routine SGESL.
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
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_LU(N,N), the LU factors from R8GE_FA.
c
c    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
c
c    Input/output, double precision B(N).
c    On input, the right hand side vector.
c    On output, the solution vector.
c
c    Input, integer JOB, specifies the operation.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision b(n)
      integer i
      integer job
      integer k
      integer l
      integer pivot(n)
      double precision t
c
c  Solve A * x = b.
c
      if ( job .eq. 0 ) then
c
c  Solve PL * Y = B.
c
        do k = 1, n - 1

          l = pivot(k)

          if ( l .ne. k ) then
            t    = b(l)
            b(l) = b(k)
            b(k) = t
          end if

          do i = k + 1, n
            b(i) = b(i) + a_lu(i,k) * b(k)
          end do

        end do
c
c  Solve U * X = Y.
c
        do k = n, 1, -1
          b(k) = b(k) / a_lu(k,k)
          do i = 1, k - 1
            b(i) = b(i) - a_lu(i,k) * b(k)
          end do
        end do
c
c  Solve A' * X = B.
c
      else
c
c  Solve U' * Y = B.
c
        do k = 1, n
          do i = 1, k - 1
            b(k) = b(k) - a_lu(i,k) * b(i)
          end do
          b(k) = b(k) / a_lu(k,k)
        end do
c
c  Solve ( PL )' * X = Y.
c
        do k = n - 1, 1, -1

          do i = k + 1, n
            b(k) = b(k) + a_lu(i,k) * b(i)
          end do

          l = pivot(k)

          if ( l .ne. k ) then
            t    = b(l)
            b(l) = b(k)
            b(k) = t
          end if

        end do

      end if

      return
      end
      subroutine r8mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc R8MAT_COPY copies an R8MAT.
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
      function r8mat_det_2d ( a )

c*********************************************************************72
c
cc R8MAT_DET_2D computes the determinant of a 2 by 2 matrix.
c
c  Discussion:
c
c    The determinant of a 2 by 2 matrix is
c
c      a11 * a22 - a12 * a21.
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
c    Input, double precision A(2,2), the matrix whose determinant is desired.
c
c    Output, double precision RMAT_DET_2D, the determinant of the matrix.
c
      implicit none

      double precision a(2,2)
      double precision r8mat_det_2d

      r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

      return
      end
      function r8mat_det_3d ( a )

c*********************************************************************72
c
cc R8MAT_DET_3D computes the determinant of a 3 by 3 matrix.
c
c  Discussion:
c
c    The determinant of a 3 by 3 matrix is
c
c        a11 * a22 * a33 - a11 * a23 * a32
c      + a12 * a23 * a31 - a12 * a21 * a33
c      + a13 * a21 * a32 - a13 * a22 * a31
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
c    Input, double precision A(3,3), the matrix whose determinant is desired.
c
c    Output, double precision RMAT_DET_3D, the determinant of the matrix.
c
      implicit none

      double precision a(3,3)
      double precision r8mat_det_3d

      r8mat_det_3d =
     &       a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) )
     &     + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) )
     &     + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

      return
      end
      function r8mat_det_4d ( a )

c*********************************************************************72
c
cc R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
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
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(4,4), the matrix whose determinant is desired.
c
c    Output, double precision R8MAT_DET_4D, the determinant of the matrix.
c
      implicit none

      double precision a(4,4)
      double precision r8mat_det_4d

      r8mat_det_4d =
     &       a(1,1) * (
     &           a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) )
     &         - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) )
     &         + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) )
     &     - a(1,2) * (
     &           a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) )
     &         - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) )
     &         + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) )
     &     + a(1,3) * (
     &           a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) )
     &         - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) )
     &         + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )
     &     - a(1,4) * (
     &           a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) )
     &         - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) )
     &         + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

      return
      end
      function r8mat_det_5d ( a )

c*********************************************************************72
c
cc R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
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
c    02 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(5,5), the matrix whose determinant is desired.
c
c    Output, double precision R8MAT_DET_5D, the determinant of the matrix.
c
      implicit none

      double precision a(5,5)
      double precision b(4,4)
      integer i
      integer inc
      integer j
      integer k
      double precision r8mat_det_4d
      double precision r8mat_det_5d
c
c  Expand the determinant into the sum of the determinants of the
c  five 4 by 4 matrices created by dropping row 1, and column k.
c
      r8mat_det_5d = 0.0D+00

      do k = 1, 5

        do i = 1, 4
          do j = 1, 4

            if ( j .lt. k ) then
              inc = 0
            else
              inc = 1
            end if

            b(i,j) = a(i+1,j+inc)

          end do
        end do

        r8mat_det_5d = r8mat_det_5d + (-1)**( k + 1 )
     &    * a(1,k) * r8mat_det_4d ( b )

      end do

      return
      end
      subroutine r8mat_inverse_2d ( a, b, det )

c*********************************************************************72
c
cc R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    If the determinant is zero, then A is singular, and does not have an
c    inverse.  In that case, B is simply set to zero, and a
c    message is printed.
c
c    If the determinant is nonzero, then its value is roughly an estimate
c    of how nonsingular the matrix A is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(2,2), the matrix to be inverted.
c
c    Output, double precision B(2,2), the inverse of the matrix A.
c
c    Output, double precision DET, the determinant of the matrix A.
c
      implicit none

      double precision a(2,2)
      double precision b(2,2)
      double precision det
      double precision r8mat_det_2d
c
c  Compute the determinant of A.
c
      det = r8mat_det_2d ( a )

      if ( det .eq. 0.0D+00 ) then

        b(1,1) = 0.0D+00
        b(1,2) = 0.0D+00
        b(2,1) = 0.0D+00
        b(2,2) = 0.0D+00

      else

        b(1,1) =  a(2,2) / det
        b(1,2) = -a(1,2) / det
        b(2,1) = -a(2,1) / det
        b(2,2) =  a(1,1) / det

      end if

      return
      end
      subroutine r8mat_inverse_3d ( a, b, det )

c*********************************************************************72
c
cc R8MAT_INVERSE_3D inverts a 3 by 3 R8MAT using Cramer's rule.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    If the determinant is zero, then A is singular, and does not have an
c    inverse.  In that case, B is simply set to zero, and a
c    message is printed.
c
c    If the determinant is nonzero, then its value is roughly an estimate
c    of how nonsingular the matrix A is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(3,3), the matrix to be inverted.
c
c    Output, double precision B(3,3), the inverse of the matrix A.
c
c    Output, double precision DET, the determinant of the matrix A.
c
      implicit none

      double precision a(3,3)
      double precision b(3,3)
      double precision det
      double precision r8mat_det_3d
c
c  Compute the determinant of A.
c
      det = r8mat_det_3d ( a )
c
c  If the determinant is zero, bail out.
c
      if ( det .eq. 0.0D+00 ) then
        call r8mat_zero ( 3, 3, b )
        return
      end if
c
c  Compute the entries of the inverse matrix using an explicit
c  formula.
c
      b(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
      b(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
      b(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

      b(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
      b(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
      b(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

      b(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
      b(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
      b(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

      return
      end
      subroutine r8mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MM multiplies two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n2)
      double precision b(n2,n3)
      double precision c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0D+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      return
      end
      subroutine r8mat_mv ( m, n, a, x, y )

c*********************************************************************72
c
cc R8MAT_MV multiplies a matrix times a vector.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c    In FORTRAN90, this operation can be more efficiently carried
c    out by the command
c
c      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of the matrix.
c
c    Input, double precision A(M,N), the M by N matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision Y(M), the product A*X.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x(n)
      double precision y(m)

      do i = 1, m
        y(i) = 0.0D+00
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
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
      subroutine r8mat_solve ( n, rhs_num, a, info )

c*********************************************************************72
c
cc R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer RHS_NUM, the number of right hand sides.
c    RHS_NUM must be at least 0.
c
c    Input/output, double precision A(N,N+rhs_num), contains in rows and
c    columns 1 to N the coefficient matrix, and in columns N+1 through
c    N+rhs_num, the right hand sides.  On output, the coefficient matrix
c    area has been destroyed, while the right hand sides have
c    been overwritten with the corresponding solutions.
c
c    Output, integer INFO, singularity flag.
c    0, the matrix was not singular, the solutions were computed;
c    J, factorization failed on step J, and the solutions could not
c    be computed.
c
      implicit none

      integer n
      integer rhs_num

      double precision a(n,n+rhs_num)
      double precision apivot
      double precision factor
      integer i
      integer info
      integer ipivot
      integer j
      integer k

      info = 0

      do j = 1, n
c
c  Choose a pivot row.
c
        ipivot = j
        apivot = a(j,j)

        do i = j+1, n
          if ( abs ( apivot ) .lt. abs ( a(i,j) ) ) then
            apivot = a(i,j)
            ipivot = i
          end if
        end do

        if ( apivot .eq. 0.0D+00 ) then
          info = j
          return
        end if
c
c  Interchange.
c
        do i = 1, n + rhs_num
          call r8_swap ( a(ipivot,i), a(j,i) )
        end do
c
c  A(J,J) becomes 1.
c
        a(j,j) = 1.0D+00
        do k = j + 1, n + rhs_num
          a(j,k) = a(j,k) / apivot
        end do
c
c  A(I,J) becomes 0.
c
        do i = 1, n

          if ( i .ne. j ) then

            factor = a(i,j)
            a(i,j) = 0.0D+00
            do k = j + 1, n + rhs_num
              a(i,k) = a(i,k) - factor * a(j,k)
            end do

          end if

        end do

      end do

      return
      end
      subroutine r8mat_solve_2d ( a, b, det, x )

c*********************************************************************72
c
cc R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
c
c  Discussion:
c
c    If DET is zero, then A is singular, and does not have an
c    inverse.  In that case, X is simply set to zero, and a
c    message is printed.
c
c    If DET is nonzero, then its value is roughly an estimate
c    of how nonsingular the matrix A is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(2,2), the matrix.
c
c    Input, double precision B(2), the right hand side.
c
c    Output, double precision DET, the determinant of the matrix A.
c
c    Output, double precision X(2), the solution of the system,
c    if DET is nonzero.
c
      implicit none

      double precision a(2,2)
      double precision b(2)
      double precision det
      double precision x(2)
c
c  Compute the determinant.
c
      det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
c
c  If the determinant is zero, bail out.
c
      if ( det .eq. 0.0D+00 ) then
        x(1) = 0.0D+00
        x(2) = 0.0D+00
        return
      end if
c
c  Compute the solution.
c
      x(1) = (  a(2,2) * b(1) - a(1,2) * b(2) ) / det
      x(2) = ( -a(2,1) * b(1) + a(1,1) * b(2) ) / det

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
      integer s_len_trim
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
      subroutine r8mat_uniform ( m, n, a, b, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM returns a scaled pseudorandom R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 February 2005
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      double precision a
      double precision b
      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

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
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8vec_angle_3d ( u, v, angle )

c*********************************************************************72
c
cc R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision U(3), V(3), the vectors.
c
c    Output, double precision ANGLE, the angle between the two vectors.
c
      implicit none

      double precision angle
      double precision angle_cos
      double precision r8_acos
      double precision r8vec_dot
      double precision u(3)
      double precision u_norm
      double precision uv_dot
      double precision v(3)
      double precision v_norm

      uv_dot = r8vec_dot ( 3, u, v )

      u_norm = sqrt ( r8vec_dot ( 3, u, u ) )

      v_norm = sqrt ( r8vec_dot ( 3, v, v ) )

      angle_cos = uv_dot / u_norm / v_norm

      angle = r8_acos ( angle_cos )

      return
      end
      subroutine r8vec_any_normal ( dim_num, v1, v2 )

c*********************************************************************72
c
cC R8VEC_ANY_NORMAL returns some normal vector to V1.
c
c  Discussion:
c
c    If DIM_NUM < 2, then no normal vector can be returned.
c
c    If V1 is the zero vector, then any unit vector will do.
c
c    No doubt, there are better, more robust algorithms.  But I will take
c    just about ANY reasonable unit vector that is normal to V1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision V1(DIM_NUM), the vector.
c
c    Output, double precision V2(DIM_NUM), a vector that is
c    normal to V2, and has unit Euclidean length.
c
      implicit none

      integer dim_num

      integer i
      integer j
      integer k
      double precision r8vec_norm
      double precision v1(dim_num)
      double precision v2(dim_num)
      double precision vj
      double precision vk

      if ( dim_num .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
        write ( *, '(a)' ) '  Called with DIM_NUM .lt. 2.'
        stop
      end if

      if ( r8vec_norm ( dim_num, v1 ) .eq. 0.0D+00 ) then
        v2(1) = 1.0D+00
        do i = 2, dim_num
          v2(i) = 0.0D+00
        end do
        return
      end if
c
c  Seek the largest entry in V1, VJ = V1(J), and the
c  second largest, VK = V1(K).
c
c  Since V1 does not have zero norm, we are guaranteed that
c  VJ, at least, is not zero.
c
      j = -1
      vj = 0.0D+00

      k = -1
      vk = 0.0D+00

      do i = 1, dim_num

        if ( abs ( vk ) .lt. abs ( v1(i) ) .or. k .lt. 1 ) then

          if ( abs ( vj ) .lt. abs ( v1(i) ) .or. j .lt. 1 ) then
            k = j
            vk = vj
            j = i
            vj = v1(i)
          else
            k = i
            vk = v1(i)
          end if

        end if

      end do
c
c  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
c  will just about do the trick.
c
      do i = 1, dim_num
        v2(i) = 0.0D+00
      end do

      v2(j) = -vk / sqrt ( vk * vk + vj * vj )
      v2(k) =  vj / sqrt ( vk * vk + vj * vj )

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
      function r8vec_cross_product_2d ( v1, v2 )

c*********************************************************************72
c
cc R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of vectors in 2D.
c
c  Discussion:
c
c    Strictly speaking, the vectors V1 and V2 should be considered
c    to lie in a 3D space, both having Z coordinate zero.  The cross
c    product value V3 then represents the standard cross product vector
c    (0,0,V3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(2), V2(2), the vectors.
c
c    Output, double precision R8VEC_CROSS_PRODUCT_2D, the cross product.
c
      implicit none

      double precision r8vec_cross_product_2d
      double precision v1(2)
      double precision v2(2)

      r8vec_cross_product_2d = v1(1) * v2(2) - v1(2) * v2(1)

      return
      end
      function r8vec_cross_product_affine_2d ( v0, v1, v2 )

c*********************************************************************72
c
cc R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
c
c  Discussion:
c
c    Strictly speaking, the vectors V1 and V2 should be considered
c    to lie in a 3D space, both having Z coordinate zero.  The cross
c    product value V3 then represents the standard cross product vector
c    (0,0,V3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V0(2), the base vector.
c
c    Input, double precision V1(2), V2(2), the vectors.
c
c    Output, double precision R8VEC_CROSS_PRODUCT_AFFINE_2D,
c    the cross product (V1-V0) x (V2-V0).
c
      implicit none

      double precision r8vec_cross_product_affine_2d
      double precision v0(2)
      double precision v1(2)
      double precision v2(2)

      r8vec_cross_product_affine_2d =
     &    ( v1(1) - v0(1) ) * ( v2(2) - v0(2) )
     &  - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

      return
      end
      subroutine r8vec_cross_product_3d ( v1, v2, v3 )

c*********************************************************************72
c
cc R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The cross product in 3D can be regarded as the determinant of the
c    symbolic matrix:
c
c          |  i  j  k |
c      det | x1 y1 z1 |
c          | x2 y2 z2 |
c
c      = ( y1 * z2 - z1 * y2 ) * i
c      + ( z1 * x2 - x1 * z2 ) * j
c      + ( x1 * y2 - y1 * x2 ) * k
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
c  Parameters:
c
c    Input, double precision V1(3), V2(3), the two vectors.
c
c    Output, double precision V3(3), the cross product vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision v1(dim_num)
      double precision v2(dim_num)
      double precision v3(dim_num)

      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

      return
      end
      subroutine r8vec_cross_product_affine_3d ( v0, v1, v2, v3 )

c*********************************************************************72
c
cc R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
c
c  Discussion:
c
c    The cross product in 3D can be regarded as the determinant of the
c    symbolic matrix:
c
c          |  i  j  k |
c      det | x1 y1 z1 |
c          | x2 y2 z2 |
c
c      = ( y1 * z2 - z1 * y2 ) * i
c      + ( z1 * x2 - x1 * z2 ) * j
c      + ( x1 * y2 - y1 * x2 ) * k
c
c    Here, we use V0 as the base of an affine system so we compute
c    the cross product of (V1-V0) and (V2-V0).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V0(3), the base vector.
c
c    Input, double precision V1(3), V2(3), the two vectors.
c
c    Output, double precision V3(3), the cross product vector
c    ( V1-V0) x (V2-V0).
c
      implicit none

      double precision v0(3)
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)

      v3(1) = ( v1(2) - v0(2) ) * ( v2(3) - v0(3) )
     &      - ( v2(2) - v0(2) ) * ( v1(3) - v0(3) )

      v3(2) = ( v1(3) - v0(3) ) * ( v2(1) - v0(1) )
     &      - ( v2(3) - v0(3) ) * ( v1(1) - v0(1) )

      v3(3) = ( v1(1) - v0(1) ) * ( v2(2) - v0(2) )
     &      - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

      return
      end
      function r8vec_diff_dot_product ( n, u1, v1, u2, v2 )

c*********************************************************************72
c
cc R8VEC_DIFF_DOT_PRODUCT: dot product of a pair of R8VEC differences.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
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
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision U1(N), V1(N), defines the vector U1-V1.
c
c    Input, double precision U2(N), V2(N), defines the vector U2-V2.
c
c    Output, double precision R8VEC_DIFF_DOT_PRODUCT, the dot product 
c    of (U1-V1)*(U2-V2).
c
      implicit none

      integer n

      integer i
      double precision r8vec_diff_dot_product
      double precision u1(n)
      double precision u2(n)
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + ( u1(i) - v1(i) ) * ( u2(i) - v2(i) )
      end do

      r8vec_diff_dot_product = value

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
      function r8vec_diff_norm_squared ( n, a, b )

c*********************************************************************72
c
cc R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A(I) - B(I) )**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2011
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
c    Output, double precision R8VEC_DIFF_NORM_SQUARED, the square
c    of the L2 norm of A - B.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision r8vec_diff_norm_squared
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + ( a(i) + b(i) ) * ( a(i) - b(i) )
      end do

      r8vec_diff_norm_squared = value

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
      function r8vec_dot_product_affine ( n, v0, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product V1-V0 * V2-V0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision V0(N), the base vector.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT_AFFINE, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product_affine
      double precision v0(n)
      double precision v1(n)
      double precision v2(n)

      r8vec_dot_product_affine = 0.0D+00

      do i = 1, n
        r8vec_dot_product_affine = r8vec_dot_product_affine
     &    + ( v1(i) - v0(i) ) * ( v2(i) - v0(i) )
      end do

      return
      end
      function r8vec_eq ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_EQ is true if every pair of entries in two vectors is equal.
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
cc R8VEC_GT: ( A1 > A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 > A2  <=>                                A1(1) > A2(1) or
c                 ( A1(1)     .eq. A2(1)     and A1(2) > A2(2) ) or
c                 ...
c                 ( A1(1:N-1) .eq. A2(1:N-1) and A1(N) > A2(N)
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
cc R8VEC_LT: ( A1 < A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 < A2  <=>                                A1(1) < A2(1) or
c                 ( A1(1)     .eq. A2(1)     and A1(2) < A2(2) ) or
c                 ...
c                 ( A1(1:N-1) .eq. A2(1:N-1) and A1(N) < A2(N)
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
      function r8vec_negative_strict ( n, a )

c*********************************************************************72
c
cc R8VEC_NEGATIVE_STRICT: every element of an R8VEC is strictly negative.
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
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N).
c
c    Output, logical R8VEC_NEGATIVE_STRICT, is TRUE every entry of the
c    vector is strictly negative.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      logical r8vec_negative_strict

      r8vec_negative_strict = .true.

      do i = 1, n
        if ( 0.0D+00 .le. a(i) ) then
          r8vec_negative_strict = .false.
          return
        end if
      end do

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
c    Output, double precision R8VEC_NORM, the L2 norm.
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
      function r8vec_norm_affine ( n, v0, v1 )

c*********************************************************************72
c
cc R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The affine vector L2 norm is defined as:
c
c      R8VEC_NORM_AFFINE(V0,V1)
c        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the vectors.
c
c    Input, double precision V0(N), the base vector.
c
c    Input, double precision V1(N), the vector whose affine norm is desired.
c
c    Output, double precision R8VEC_NORM_AFFINE, the affine L2.
c
      implicit none

      integer n

      integer i
      double precision r8vec_norm_affine
      double precision v0(n)
      double precision v1(n)

      r8vec_norm_affine = 0.0D+00
      do i = 1, n
        r8vec_norm_affine = r8vec_norm_affine
     &    + ( v0(i) - v1(i) )**2
      end do
      r8vec_norm_affine = sqrt ( r8vec_norm_affine )

      return
      end
      function r8vec_norm_squared ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM_SQUARED returns the square of the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    R8VEC_NORM_SQUARED = sum ( 1 <= I <= N ) A(I)**2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector.
c
c    Output, double precision R8VEC_NORM_SQUARED, the square of the L2 norm.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm_squared
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do

      r8vec_norm_squared = value

      return
      end
      subroutine r8vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
c
c    The Box-Muller method is used, which is efficient, but
c    generates an even number of values each time.  On any call
c    to this routine, an even number of new values are generated.
c    Depending on the situation, one value may be left over.
c    In that case, it is saved for the next call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values desired.  If N is negative,
c    then the code will flush its internal memory; in particular,
c    if there is a saved value to be used on the next call, it is
c    instead discarded.  This is useful if the user has reset the
c    random number seed, for instance.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision X(N), a sample of the standard normal PDF.
c
c  Local parameters:
c
c    Local, integer MADE, records the number of values that have
c    been computed.  On input with negative N, this value overwrites
c    the return value of N, so the user can get an accounting of
c    how much work has been done.
c
c    Local, integer SAVED, is 0 or 1 depending on whether there is a
c    single saved value left over from the previous call.
c
c    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
c    X that we need to compute.  This starts off as 1:N, but is adjusted
c    if we have a saved value that can be immediately stored in X(1),
c    and so on.
c
c    Local, double precision Y, the value saved from the previous call, if
c    SAVED is 1.
c
      implicit none

      integer n

      integer i
      integer m
      integer made
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(2)
      double precision r8_uniform_01
      integer saved
      integer seed
      double precision x(n)
      integer x_hi_index
      integer x_lo_index
      double precision y

      save made
      save saved
      save y

      data made / 0 /
      data saved / 0 /
      data y / 0.0D+00 /
c
c  I'd like to allow the user to reset the internal data.
c  But this won't work properly if we have a saved value Y.
c  I'm making a crock option that allows the user to signal
c  explicitly that any internal memory should be flushed,
c  by passing in a negative value for N.
c
      if ( n .lt. 0 ) then
        n = made
        made = 0
        saved = 0
        y = 0.0D+00
        return
      else if ( n .eq. 0 ) then
        return
      end if
c
c  Record the range of X we need to fill in.
c
      x_lo_index = 1
      x_hi_index = n
c
c  Use up the old value, if we have it.
c
      if ( saved .eq. 1 ) then
        x(1) = y
        saved = 0
        x_lo_index = 2
      end if
c
c  Maybe we don't need any more values.
c
      if ( x_hi_index - x_lo_index + 1 .eq. 0 ) then
c
c  If we need just one new value, do that here to avoid null arrays.
c
      else if ( x_hi_index - x_lo_index + 1 .eq. 1 ) then

        r(1) = r8_uniform_01 ( seed )

        if ( r(1) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        r(2) = r8_uniform_01 ( seed )

        x(x_hi_index) =
     &           sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * cos ( 2.0D+00 * pi * r(2) )
        y =      sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + 2
c
c  If we require an even number of values, that's easy.
c
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) .eq. 0 ) then

        do i = x_lo_index, x_hi_index, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        made = made + x_hi_index - x_lo_index + 1
c
c  If we require an odd number of values, we generate an even number,
c  and handle the last pair specially, storing one in X(N), and
c  saving the other for later.
c
      else

        do i = x_lo_index, x_hi_index - 1, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        call r8vec_uniform_01 ( 2, seed, r )

        x(n) = sqrt ( -2.0D+00 * log ( r(1) ) )
     &    * cos ( 2.0D+00 * pi * r(1) )

        y = sqrt ( -2.0D+00 * log ( r(2) ) )
     &    * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + x_hi_index - x_lo_index + 2

      end if

      return
      end
      function r8vec_normsq ( n, v )

c*********************************************************************72
c
cc R8VEC_NORMSQ returns the square of the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The square of the vector L2 norm is defined as:
c
c      R8VEC_NORMSQ = sum ( 1 <= I <= N ) V(I)^2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the vector dimension.
c
c    Input, double precision V(N), the vector.
c
c    Output, double precision R8VEC_NORMSQ, the squared L2 norm.
c
      implicit none

      integer n

      double precision r8vec_normsq
      double precision v(n)

      r8vec_normsq = sum ( v(1:n)**2 )

      return
      end
      function r8vec_normsq_affine ( n, v0, v1 )

c*********************************************************************72
c
cc R8VEC_NORMSQ_AFFINE returns the affine squared norm of an R8VEC.
c
c  Discussion:
c
c   An R8VEC is a vector of R8's.
c
c    The affine squared vector L2 norm is defined as:
c
c      R8VEC_NORMSQ_AFFINE(V0,V1)
c        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the vector dimension.
c
c    Input, double precision V0(N), the base vector.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_NORMSQ_AFFINE, the squared affine L2 norm.
c
      implicit none

      integer n

      integer i
      double precision r8vec_normsq_affine
      double precision v0(n)
      double precision v1(n)

      r8vec_normsq_affine = 0.0D+00
      do i = 1, n
        r8vec_normsq_affine = r8vec_normsq_affine + ( v0(i) - v1(i) )**2
      end do

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
      function r8vec_positive_strict ( n, a )

c*********************************************************************72
c
cc R8VEC_POSITIVE_STRICT: every element of an R8VEC is strictly positive.
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
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N).
c
c    Output, logical R8VEC_POSITIVE_STRICT, is TRUE every entry of the
c    vector is strictly positive.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      logical r8vec_positive_strict

      r8vec_positive_strict = .true.

      do i = 1, n
        if ( a(i) .le. 0.0D+00 ) then
          r8vec_positive_strict = .false.
          return
        end if
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
      function r8vec_scalar_triple_product ( v1, v2, v3 )

c*********************************************************************72
c
cc R8VEC_SCALAR_TRIPLE_PRODUCT finds the scalar triple product in 3D.
c
c  Discussion:
c
c    [A,B,C] = A dot ( B cross C )
c            = B dot ( C cross A )
c            = C dot ( A cross B )
c
c    The volume of a parallelepiped, whose sides are given by
c    vectors A, B, and C, is abs ( A dot ( B cross C ) ).
c
c    Three vectors are coplanar if and only if their scalar triple
c    product vanishes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    "Scalar Triple Product",
c    CRC Concise Encyclopedia of Mathematics,
c    CRC, 1999
c
c  Parameters:
c
c    Input, double precision V1(3), V2(3), V3(3), the vectors.
c
c    Output, double precision R8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
c    triple product.
c
      implicit none

      double precision r8vec_dot_product
      double precision r8vec_scalar_triple_product
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
      double precision v4(3)

      call r8vec_cross_product_3d ( v2, v3, v4 )

      r8vec_scalar_triple_product = r8vec_dot_product ( 3, v1, v4 )

      return
      end
      subroutine r8vec_swap ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_SWAP swaps the entries of two R8VEC's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
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
c    17 July 2006
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
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_uniform_ab ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
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
        write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
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
      subroutine r8vec_uniform_unit ( m, seed, w )

c*********************************************************************72
c
cc R8VEC_UNIFORM_UNIT generates a uniformly random unit vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision W(M), a random direction vector,
c    with unit norm.
c
      implicit none

      integer m

      integer i
      double precision norm
      double precision r8vec_norm_l2
      integer seed
      double precision w(m)
c
c  Get N values from a standard normal distribution.
c
      call r8vec_normal_01 ( m, seed, w )
c
c  Compute the length of the vector.
c
      norm = r8vec_norm_l2 ( m, w )
c
c  Normalize the vector.
c
      do i = 1, m
        w(i) = w(i) / norm
      end do

      return
      end
      subroutine radec_distance_3d ( ra1, dec1, ra2, dec2, theta )

c*********************************************************************72
c
cc RADEC_DISTANCE_3D - angular distance, astronomical units, sphere in 3D.
c
c  Discussion:
c
c    Right ascension is measured in hours, between 0 and 24, and
c    essentially measures longitude.
c
c    Declination measures the angle from the equator towards the north pole,
c    and ranges from -90 (South Pole) to 90 (North Pole).
c
c    On the unit sphere, the angular separation between two points is
c    equal to their geodesic or great circle distance.  On any other
c    sphere, multiply the angular separation by the radius of the
c    sphere to get the geodesic or great circle distance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RA1, DEC1, RA2, DEC2, the right ascension and
c    declination of the two points.
c
c    Output, double precision THETA, the angular separation between the points,
c    in radians.
c
      implicit none

      double precision cos_theta
      double precision dec1
      double precision dec2
      double precision degrees_to_radians
      double precision norm_v1
      double precision norm_v2
      double precision phi1
      double precision phi2
      double precision r8_acos
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision ra1
      double precision ra2
      double precision theta
      double precision theta1
      double precision theta2
      double precision v1(3)
      double precision v2(3)

      theta1 = degrees_to_radians ( 15.0D+00 * ra1 )
      phi1 = degrees_to_radians ( dec1 )

      v1(1) = cos ( theta1 ) * cos ( phi1 )
      v1(2) = sin ( theta1 ) * cos ( phi1 )
      v1(3) = sin ( phi1 )

      norm_v1 = r8vec_norm ( 3, v1 )

      theta2 = degrees_to_radians ( 15.0D+00 * ra2 )
      phi2 = degrees_to_radians ( dec2 )

      v2(1) = cos ( theta2 ) * cos ( phi2 )
      v2(2) = sin ( theta2 ) * cos ( phi2 )
      v2(3) = sin ( phi2 )

      norm_v2 = r8vec_norm ( 3, v2 )

      cos_theta = r8vec_dot_product ( 3, v1, v2 )
     &  / ( norm_v1 * norm_v2 )

      theta = r8_acos ( cos_theta )

      return
      end
      subroutine radec_to_xyz ( ra, dec, p )

c*********************************************************************72
c
cc RADEC_TO_XYZ converts right ascension/declination to (X,Y,Z) coordinates.
c
c  Discussion:
c
c    Right ascension is measured in hours, between 0 and 24, and
c    essentially measures longitude.
c
c    Declination measures the angle from the equator towards the north pole,
c    and ranges from -90 (South Pole) to 90 (North Pole).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RA, DEC, the right ascension and declination
c    of a point.
c
c    Output, double precision P(3), the corresponding coordinates of
c    a point with radius 1.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dec
      double precision degrees_to_radians
      double precision p(dim_num)
      double precision phi
      double precision ra
      double precision theta

      theta = degrees_to_radians ( 15.0D+00 * ra )
      phi = degrees_to_radians ( dec )

      p(1) = cos ( theta ) * cos ( phi )
      p(2) = sin ( theta ) * cos ( phi )
      p(3) = sin ( phi )

      return
      end
      function radians_to_degrees ( angle_rad )

c*********************************************************************72
c
cc RADIANS_TO_DEGREES converts an angle from radians to degrees.
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
c    Input, double precision ANGLE_RAD, an angle in radians.
c
c    Output, double precision RADIANS_TO_DEGREES, the equivalent angle
c    in degrees.
c
      implicit none

      double precision angle_rad
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision radians_to_degrees

      radians_to_degrees = ( angle_rad / pi ) * 180.0D+00

      return
      end
      subroutine radians_to_dms ( angle_rad, degrees, minutes, seconds )

c*********************************************************************72
c
cc RADIANS_TO_DMS converts an angle from radians to degrees/minutes/seconds.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_RAD, the angle in radians.
c
c    Output, integer DEGREES, MINUTES, SECONDS, the equivalent
c    angle in degrees, minutes, and seconds.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      integer degrees
      integer minutes
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seconds

      angle_deg = 180.0D+00 * abs ( angle_rad ) / pi

      degrees = int ( angle_deg )
      angle_deg = ( angle_deg - dble ( degrees ) ) * 60.0D+00
      minutes = int ( angle_deg )
      angle_deg = ( angle_deg - dble ( minutes ) ) * 60.0D+00
      seconds = nint ( angle_deg )

      if ( angle_rad .lt. 0.0D+00 ) then
        degrees = - degrees
        minutes = - minutes
        seconds = - seconds
      end if

      return
      end
      subroutine rotation_axis_vector_3d ( axis, angle, v, w )

c*********************************************************************72
c
cc ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
c
c  Discussion:
c
c    Thanks to Cody Farnell for correcting some mistakes in an earlier
c    version of this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision AXIS(3), the axis vector for the rotation.
c
c    Input, double precision ANGLE, the angle, in radians, of the rotation.
c
c    Input, double precision V(3), the vector to be rotated.
c
c    Output, double precision W(3), the rotated vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision angle
      double precision axis(dim_num)
      double precision axis_norm
      double precision dot
      integer i
      double precision norm
      double precision normal(dim_num)
      double precision normal_component
      double precision normal2(dim_num)
      double precision parallel(dim_num)
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision rot(dim_num)
      double precision u(dim_num)
      double precision v(dim_num)
      double precision w(dim_num)
c
c  Compute the length of the rotation axis.
c
      do i = 1, dim_num
        u(i) = axis(i)
      end do

      axis_norm = r8vec_norm ( dim_num, u )

      if ( axis_norm .eq. 0.0D+00 ) then
        w(1:dim_num) = 0.0D+00
        return
      end if

      u(1:dim_num) = u(1:dim_num) / axis_norm
c
c  Compute the dot product of the vector and the unit rotation axis.
c
      dot = r8vec_dot_product ( dim_num, u, v )
c
c  Compute the parallel component of the vector.
c
      do i = 1, dim_num
        parallel(i) = dot * u(i)
      end do
c
c  Compute the normal component of the vector.
c
      do i = 1, dim_num
        normal(i) = v(i) - parallel(i)
      end do

      normal_component = r8vec_norm ( dim_num, normal )

      if ( normal_component .eq. 0.0D+00 ) then
        do i = 1, dim_num
          w(i) = parallel(i)
        end do
        return
      end if

      do i = 1, dim_num
        normal(i) = normal(i) / normal_component
      end do
c
c  Compute a second vector, lying in the plane, perpendicular
c  to V, and forming a right-handed system, as the cross product
c  of the first two vectors.
c
      normal2(1) = u(2) * normal(3) - u(3) * normal(2)
      normal2(2) = u(3) * normal(1) - u(1) * normal(3)
      normal2(3) = u(1) * normal(2) - u(2) * normal(1)

      norm = r8vec_norm ( dim_num, normal2 )

      do i = 1, dim_num
        normal2(i) = normal2(i) / norm
      end do
c
c  Rotate the normal component by the angle.
c
      do i = 1, dim_num
        rot(i) = normal_component * (
     &      cos ( angle ) * normal(i)
     &    + sin ( angle ) * normal2(i) )
      end do
c
c  The rotated vector is the parallel component plus the rotated component.
c
      do i = 1, dim_num
        w(i) = parallel(i) + rot(i)
      end do

      return
      end
      subroutine rotation_axis2mat_3d ( axis, angle, a )

c*********************************************************************72
c
cc ROTATION_AXIS2MAT_3D converts a rotation from axis to matrix format in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision AXIS(3), the axis vector which remains 
c    unchanged by the rotation.
c
c    Input, double precision ANGLE, the angular measurement of the
c    rotation about the axis, in radians.
c
c    Output, double precision A(3,3), the rotation matrix.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a(dim_num,dim_num)
      double precision angle
      double precision axis(dim_num)
      double precision axis_norm
      double precision ca
      integer i
      integer j
      double precision r8vec_norm
      double precision sa
      double precision v1
      double precision v2
      double precision v3

      v1 = axis(1)
      v2 = axis(2)
      v3 = axis(3)

      axis_norm = r8vec_norm ( dim_num, axis )

      if ( axis_norm .eq. 0.0D+00 ) then
        do j = 1, dim_num
          do i = 1, dim_num
            a(i,j) = 0.0D+00
          end do
        end do
        return
      end if

      v1 = v1 / axis_norm
      v2 = v2 / axis_norm
      v3 = v3 / axis_norm

      ca = cos ( angle )
      sa = sin ( angle )

      a(1,1) =                    v1 * v1 + ca * ( 1.0D+00 - v1 * v1 )
      a(1,2) = ( 1.0D+00 - ca ) * v1 * v2 - sa * v3
      a(1,3) = ( 1.0D+00 - ca ) * v1 * v3 + sa * v2

      a(2,1) = ( 1.0D+00 - ca ) * v2 * v1 + sa * v3
      a(2,2) =                    v2 * v2 + ca * ( 1.0D+00 - v2 * v2 )
      a(2,3) = ( 1.0D+00 - ca ) * v2 * v3 - sa * v1

      a(3,1) = ( 1.0D+00 - ca ) * v3 * v1 - sa * v2
      a(3,2) = ( 1.0D+00 - ca ) * v3 * v2 + sa * v1
      a(3,3) =                    v3 * v3 + ca * ( 1.0D+00 - v3 * v3 )

      return
      end
      subroutine rotation_axis2quat_3d ( axis, angle, q )

c*********************************************************************72
c
cc ROTATION_AXIS2QUAT_3D converts rotation from axis to quaternion form in 3D.
c
c  Discussion:
c
c    A rotation quaternion Q has the form:
c
c      Q = A + Bi + Cj + Dk
c
c    where A, B, C and D are real numbers, and i, j, and k are to be regarded
c    as symbolic constant basis vectors, similar to the role of the "i"
c    in the representation of imaginary numbers.
c
c    A is the cosine of half of the angle of rotation.  (B,C,D) is a
c    unit vector pointing in the direction of the axis of rotation.
c    Rotation multiplication and inversion can be carried out using
c    this format and the usual rules for quaternion multiplication
c    and inversion.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision AXIS(3), the axis vector which remains 
c    unchanged by the rotation.
c
c    Input, double precision ANGLE, the angular measurement of the 
c    rotation about the axis, in radians.
c
c    Output, double precision Q(4), the quaternion representing the rotation.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision axis(dim_num)
      double precision axis_norm
      double precision angle
      integer i
      double precision q(4)
      double precision r8vec_norm

      axis_norm = r8vec_norm ( dim_num, axis )

      if ( axis_norm .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ROTATION_AXIS2QUAT_3D - Fatal error!'
        write ( *, '(a)' ) '  The axis vector is null.'
        q(1:4) = 0.0D+00
        stop
      end if

      q(1)   = cos ( 0.5D+00 * angle )
      do i = 1, 3
        q(i+1) = sin ( 0.5D+00 * angle ) * axis(i) / axis_norm
      end do

      return
      end
      subroutine rotation_mat_vector_3d ( a, v, w )

c*********************************************************************72
c
cc ROTATION_MAT_VECTOR_3D applies a marix rotation to a vector in 3d.
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
c  Parameters:
c
c    Input, double precision A(3,3), the matrix defining the rotation.
c
c    Input, double precision V(3), the vector to be rotated.
c
c    Output, double precision W(3), the rotated vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a(dim_num,dim_num)
      double precision v(dim_num)
      double precision w(dim_num)

      call r8mat_mv ( dim_num, dim_num, a, v, w )

      return
      end
      subroutine rtp_to_xyz ( r, theta, phi, xyz )

c*********************************************************************72
c
cc RTP_TO_XYZ converts (R,Theta,Phi) to (X,Y,Z) coordinates.
c
c  Discussion:
c
c    R measures the distance of the point to the origin.
c
c    Theta measures the "longitude" of the point, between 0 and 2 PI.
c
c    PHI measures the angle from the "north pole", between 0 and PI.
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
c    Input, double precision R, THETA, PHI, the radius, longitude, and
c    declination of a point.
c
c    Output, double precision XYZ(3), the corresponding Cartesian coordinates.
c
      implicit none

      double precision phi
      double precision r
      double precision theta
      double precision xyz(3)

      xyz(1) = r * cos ( theta ) * sin ( phi )
      xyz(2) = r * sin ( theta ) * sin ( phi )
      xyz(3) = r *                 cos ( phi )

      return
      end
      function sec_deg ( angle_deg )

c*********************************************************************72
c
cc SEC_DEG returns the secant of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_DEG, the angle, in degrees.
c
c    Output, double precision SEC_DEG, the secant of the angle.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      double precision degrees_to_radians 
      parameter ( degrees_to_radians 
     &  = 3.141592653589793D+00 / 180.0D+00 )
      double precision sec_deg

      angle_rad = degrees_to_radians * angle_deg
      sec_deg = 1.0D+00 / cos ( angle_rad )

      return
      end
      subroutine segment_contains_point_1d ( p1, p2, p, t )

c*********************************************************************72
c
cc SEGMENT_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
c
c  Discussion:
c
c    A line segment is the finite portion of a line that lies between
c    two points P1 and P2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1, P2, two points defining a line segment.
c    The line segment has T = 0 at P1, and T = 1 at P2.
c
c    Input, double precision P, a point to be tested.
c
c    Output, double precision T, the coordinate of P3 in units of (P2-P1).
c    The point P3 is contained in the line segment if 0 <= T <= 1.
c
      implicit none

      double precision p
      double precision p1
      double precision p2
      double precision r8_huge
      double precision t
      double precision unit

      unit = p2 - p1

      if ( unit .eq. 0.0D+00 ) then

        if ( p .eq. p1 ) then
          t = 0.5D+00
        else if ( p .lt. p1 ) then
          t = - r8_huge ( t )
        else if ( p1 .lt. p ) then
          t = r8_huge ( t )
        end if

      else

        t = ( p - p1 ) / unit

      end if

      return
      end
      subroutine segment_contains_point_2d ( p1, p2, p, u )

c*********************************************************************72
c
cc SEGMENT_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
c
c  Discussion:
c
c    A line segment is the finite portion of a line that lies between
c    two points P1 and P2.
c
c    In exact arithmetic, point P is on the line segment between
c    P1 and P2 if and only if 0 <= U <= 1 and V = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), the endpoints of the line segment.
c
c    Input, double precision P(2), a point to be tested.
c
c    Output, double precision U(2), the components of P, with the first
c    component measured along the axis with origin at P1 and unit at P2, 
c    and second component the magnitude of the off-axis portion of the
c    vector P-P1, measured in units of (P2-P1).
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision normsq
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision r8_huge
      double precision u(dim_num)

      normsq = ( p2(1) - p1(1) )**2 + ( p2(2) - p1(2) )**2

      if ( normsq .eq. 0.0D+00 ) then

        if ( p(1) .eq. p1(1) .and. p(2) .eq. p(2) ) then
          u(1) = 0.5D+00
          u(2) = 0.0D+00
        else
          u(1) = 0.5D+00
          u(2) = r8_huge ( )
        end if

      else

        u(1) = ( ( p(1) - p1(1) ) * ( p2(1) - p1(1) ) 
     &         + ( p(2) - p1(2) ) * ( p2(2) - p1(2) ) ) / normsq

        u(2) = sqrt ( ( ( u(1) - 1.0D+00 ) * p1(1) 
     &                  - u(1) * p2(1) + p(1) )**2 
     &              + ( ( u(1) - 1.0D+00 ) * p1(2) 
     &                  - u(1) * p2(2) + p(2) )**2 ) 
     &              / sqrt ( normsq )
 
      end if

      return
      end
      subroutine segment_point_coords_2d ( p1, p2, p, s, t )

c*********************************************************************72
c
cc SEGMENT_POINT_COORDS_2D: coordinates of a point on a line segment in 2D.
c
c  Discussion:
c
c    A line segment is the finite portion of a line that lies between
c    two points P1 and P2.
c
c    By the coordinates of a point P with respect to a line segment [P1,P2]
c    we mean numbers S and T such that S gives us the distance from the
c    point P to the nearest point PN on the line (not the line segment!), 
c    and T gives us the position of PN relative to P1 and P2.
c
c    If S is zero, then P lies on the line.
c
c    If 0 <= T <= 1, then PN lies on the line segment.
c
c    If both conditions hold, then P lies on the line segment.
c
c    If E is the length of the line segment, then the distance of the 
c    point to the line segment is:
c
c      sqrt ( S^2 +  T^2    * E^2 )     if      T <= 0;
c             S                         if 0 <= T <= 1
c      sqrt ( S^2 + (T-1)^2 * E^2 )     if 1 <= T
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), the endpoints of the line segment.
c
c    Input, double precision P(2), the point to be considered.
c
c    Output, double precision S, the distance of P to the nearest point PN
c    on the line through P1 and P2.  (S will always be nonnegative.)
c
c    Output, double precision T, the relative position of the point PN
c    to the points P1 and P2.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision bot
      integer i
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pn(dim_num)
      double precision s
      double precision t
c
c  If the line segment is actually a point, then the answer is easy.
c
      if ( p1(1) .eq. p2(1) .and. p1(2) .eq. p2(2) ) then

        t = 0.0D+00

      else

        bot = ( p2(1) - p1(1) )**2 + ( p2(2) - p1(2) )**2

        t = ( ( p(1) - p1(1) ) * ( p2(1) - p1(1) ) 
     &      + ( p(2) - p1(2) ) * ( p2(2) - p1(2) ) ) / bot

      end if

      do i = 1, dim_num
        pn(i) = p1(i) + t * ( p2(i) - p1(i) )
      end do

      s = sqrt ( ( p(1) - pn(1) )**2 + ( p(2) - pn(2) )**2 )

      return
      end
      subroutine segments_curvature_2d ( p1, p2, p3, curvature )

c*********************************************************************72
c
cc SEGMENTS_CURVATURE_2D computes the curvature of two line segments in 2D.
c
c  Discussion:
c
c    A line segment is the finite portion of a line that lies between
c    two points P1 and P2.
c
c    We assume that the segments are [P1,P2] and [P2,P3].
c
c    We compute the circle that passes through P1, P2 and P3.
c
c    The inverse of the radius of this circle is the local "curvature"
c    associated with the three points.
c
c    If curvature is 0, the two line segments have the same slope,
c    and the three points are collinear.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), P2(2), P3(2), the points.
c
c    Output, double precision CURVATURE, the local curvature.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision curvature
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision pc(dim_num)
      double precision r

      call circle_exp2imp_2d ( p1, p2, p3, r, pc )

      if ( 0.0D+00 .lt. r ) then
        curvature = 1.0D+00 / r
      else
        curvature = 0.0D+00
      end if

      return
      end
      subroutine simplex_lattice_layer_point_next ( n, c, v, more )

c*********************************************************************72
c
cc SIMPLEX_LATTICE_LAYER_POINT_NEXT: next simplex lattice layer point.
c
c  Discussion:
c
c    The simplex lattice layer L is bounded by the lines
c
c      0 <= X(1:N),
c      L - 1 < sum X(1:N) / C(1:N)  <= L.
c
c    In particular, layer L = 0 always contains just the origin.
c
c    This function returns, one at a time, the points that lie within
c    a given simplex lattice layer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer C(N+1), coefficients defining the
c    lattice layer in entries 1 to N, and the laver index in C(N+1).
c    The coefficients should be positive, and C(N+1) must be nonnegative.
c
c    Input/output, integer V(N).  On first call for a given layer,
c    the input value of V is not important.  On a repeated call for the same
c    layer, the input value of V should be the output value from the previous
c    call.  On output, V contains the next lattice layer point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given layer.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if the returned value V is a new point.
c    If the output value is FALSE, then no more points were found,
c    and V was reset to 0, and the lattice layer has been exhausted.
c
      implicit none

      integer n

      integer c(n+1)
      integer c1n
      integer i
      integer i4vec_lcm
      integer j
      integer lhs
      logical more
      integer rhs1
      integer rhs2
      integer v(n)
c
c  Treat layer C(N+1) = 0 specially.
c
      if ( c(n+1) .eq. 0 ) then
        if ( .not. more ) then
          do i = 1, n
            v(i) = 0
          end do
          more = .true.
        else
          more = .false.
        end if
        return
      end if
c
c  Compute the first point.
c
      if ( .not. more ) then

        v(1) = ( c(n+1) - 1 ) * c(1) + 1
        do i = 2, n
          v(i) = 0
        end do
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )

        rhs1 = c1n * ( c(n+1) - 1 )
        rhs2 = c1n *   c(n+1)
c
c  Try to increment component I.
c
        do i = 1, n

          v(i) = v(i) + 1

          do j = 1, i - 1
            v(j) = 0
          end do

          if ( 1 .lt. i ) then
            v(1) = rhs1
            do j = 2, n
              v(1) = v(1) - ( c1n / c(j) ) * v(j)
            end do
            v(1) = ( c(1) * v(1) ) / c1n
            v(1) = max ( v(1), 0 )
          end if

          lhs = 0
          do j = 1, n
            lhs = lhs + ( c1n / c(j) ) * v(j)
          end do

c         write (  *, * ) 'LHS = ', lhs, 'RHS1 = ', rhs1

          if ( lhs .le. rhs1 ) then
            v(1) = v(1) + 1
            lhs = lhs + c1n / c(1)
          end if

          if ( lhs .le. rhs2 ) then
            return
          end if

        end do

        do j = 1, n
          v(j) = 0
        end do
        more = .false.

      end if

      return
      end
      subroutine simplex_lattice_point_next ( n, c, v, more )

c*********************************************************************72
c
cc SIMPLEX_LATTICE_POINT_NEXT returns the next simplex lattice point.
c
c  Discussion:
c
c    The lattice simplex is defined by the vertices:
c
c      (0,0,...,0), (C(N+1)/C(1),0,...,0), (0,C(N+1)/C(2),...,0) ...
c      (0,0,...C(N+1)/C(N))
c
c    The lattice simplex is bounded by the lines
c
c      0 <= V(1:N),
c      V(1) / C(1) + V(2) / C(2) + ... + V(N) / C(N) <= C(N+1)
c
c    Lattice points are listed one at a time, starting at the origin,
c    with V(1) increasing first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, integer C(N+1), coefficients defining the
c    lattice simplex.  These should be positive.
c
c    Input/output, integer V(N).  On first call, the input
c    value is not important.  On a repeated call, the input value should
c    be the output value from the previous call.  On output, V contains
c    the next lattice point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given simplex.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if not only is the returned value V a lattice point,
c    but the routine can be called again for another lattice point.
c    If the output value is FALSE, then no more lattice points were found,
c    and V was reset to 0, and the routine should not be called further
c    for this simplex.
c
      implicit none

      integer n

      integer c(n+1)
      integer c1n
      integer i
      integer i4vec_lcm
      integer j
      integer lhs
      logical more
      integer rhs
      integer term
      integer v(n)

      if ( .not. more ) then

        call i4vec_zero ( n, v )
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )
        rhs = c1n * c(n+1)

        lhs = 0
        do i = 1, n
          term = 1
          do j = 1, n
            if ( i .eq. j ) then
              term = term * v(j)
            else
              term = term * c(j)
            end if
          end do
          lhs = lhs + term
        end do

        do i = 1, n
          if ( lhs + c1n / c(i) .le. rhs ) then
            v(i) = v(i) + 1
            more = .true.
            return
          end if
          lhs = lhs - c1n * v(i) / c(i)
          v(i) = 0
        end do

        more = .false.

      end if

      return
      end
      subroutine simplex_unit_lattice_point_num_nd ( d, s, n )

c*********************************************************************72
c
cc SIMPLEX_UNIT_LATTICE_POINT_NUM_ND: count lattice points.
c
c  Discussion:
c
c    The simplex is assumed to be the unit D-dimensional simplex:
c
c    ( (0,0,...,0), (1,0,...,0), (0,1,...,0), ... (0,,0,...,1) )
c
c    or a copy of this simplex scaled by an integer S:
c
c    ( (0,0,...,0), (S,0,...,0), (0,S,...,0), ... (0,,0,...,S) )
c
c    The routine returns the number of integer lattice points that appear
c    inside the simplex or on its boundary.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Matthias Beck, Sinai Robins,
c    Computing the Continuous Discretely,
c    Springer, 2006,
c    ISBN13: 978-0387291390,
c    LC: QA640.7.B43.
c
c  Parameters:
c
c    Input, integer D, the spatial dimension.
c
c    Input, integer S, the scale factor.
c
c    Output, integer N, the number of lattice points.
c
      implicit  none

      integer d
      integer i
      integer n
      integer s

      n = 1
      do i = 1, d
        n = ( n * ( s + i ) ) / i
      end do

      return
      end
      subroutine simplex_unit_volume_nd ( dim_num, volume )

c*********************************************************************72
c
cc SIMPLEX_UNIT_VOLUME_ND computes the volume of the unit simplex in ND.
c
c  Discussion:
c
c    The formula is simple: volume = 1/N!.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Output, double precision VOLUME, the volume.
c
      implicit none

      integer i
      integer dim_num
      double precision volume

      volume = 1.0D+00
      do i = 1, dim_num
        volume = volume / dble ( i )
      end do

      return
      end
      subroutine simplex_volume_nd ( dim_num, a, volume )

c*********************************************************************72
c
cc SIMPLEX_VOLUME_ND computes the volume of a simplex in ND.
c
c  Discussion:
c
c    The formula is: 
c
c      volume = 1/N! * det ( A )
c
c    where A is the N by N matrix obtained by subtracting one
c    vector from all the others.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Input, double precision A(DIM_NUM,DIM_NUM+1), the vertices.
c
c    Output, double precision VOLUME, the volume of the simplex.
c
      implicit none

      integer dim_num

      double precision a(dim_num,dim_num+1)
      double precision b(dim_num,dim_num)
      double precision det
      integer i
      integer info
      integer j
      integer pivot(dim_num)
      double precision volume

      do j = 1, dim_num
        do i = 1, dim_num
          b(i,j) = a(i,j)
        end do
      end do

      do j = 1, dim_num
        do i = 1, dim_num
          b(i,j) = b(i,j) - a(i,dim_num+1)
        end do
      end do

      call r8ge_fa ( dim_num, b, pivot, info )

      if ( info .ne. 0 ) then

        volume = -1.0D+00

      else

        call r8ge_det ( dim_num, b, pivot, det )

        volume = abs ( det )
        do i = 1, dim_num
          volume = volume / dble ( i )
        end do

      end if

      return
      end
      function sin_deg ( angle )

c*********************************************************************72
c
cc SIN_DEG returns the sine of an angle given in degrees.
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
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision SIN_DEG, the sine of the angle.
c
      implicit none

      double precision angle
      double precision angle_rad
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sin_deg

      angle_rad = pi * angle / 180.0D+00
      sin_deg  = sin ( angle_rad )

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
c    30 March 2011
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
        value = 0.0D+00
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
     &            + sa**( m - 1 ) * ca - sb**( m - 1 ) * cb ) 
     &    / dble ( m )
      end do

      sin_power_int = value

      return
      end
      subroutine soccer_size_3d ( point_num, edge_num, face_num,
     &  face_order_max )

c*********************************************************************72
c
cc SOCCER_SIZE_3D gives "sizes" for a truncated icosahedron in 3D.
c
c  Discussion:
c
c    The shape is a truncated icosahedron, which is the design used
c    on a soccer ball.  There are 12 pentagons and 20 hexagons.
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
c  Reference:
c
c    http://mathworld.wolfram.com/TruncatedIcosahedron.html
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

      point_num = 60
      edge_num = 90
      face_num = 32
      face_order_max = 6

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

      double precision area
      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_asin
      double precision theta

      if ( h .le. 0.0D+00 ) then
        area = 0.0D+00
      else if ( 2.0D+00 * r .le. h ) then
        area = 2.0D+00 * pi * r
      else

        theta = 2.0D+00
     &    * r8_asin ( sqrt ( r * r - ( r - h )**2 ) / r )
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

      double precision area
      double precision area2
      double precision h
      double precision haver_sine
      integer i
      integer dim_num
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_asin
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
        call sphere_imp_area_nd ( dim_num, r, area )
        return
      end if
c
c  For cases where R < H < 2 * R, work with the complementary region.
c
      haver_sine = sqrt ( ( 2.0D+00 * r - h ) * h )

      theta = r8_asin ( haver_sine / r )

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
        call sphere_imp_area_nd ( dim_num, r, area2 )
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

      double precision h
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r8_asin
      double precision theta
      double precision volume

      if ( h .le. 0.0D+00 ) then

        volume = 0.0D+00

      else if ( 2.0D+00 * r .le. h ) then

        volume = pi * r * r

      else

        theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h )**2 ) / r )
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
      double precision factor1
      double precision factor2
      double precision h
      integer dim_num
      double precision r
      double precision r8_asin
      double precision sin_power_int
      double precision sphere_unit_volume_nd
      double precision volume
      double precision volume2

      if ( h .le. 0.0D+00 ) then
        volume = 0.0D+00
        return
      end if

      if ( 2.0D+00 * r .le. h ) then
        call sphere_imp_volume_nd ( dim_num, r, volume )
        return
      end if

      if ( dim_num .lt. 1 ) then

        volume = - 1.0D+00

      else if ( dim_num .eq. 1 ) then

        volume = h

      else

        factor1 = sphere_unit_volume_nd ( dim_num - 1 )

        angle = r8_asin ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) )

        factor2 = sin_power_int ( 0.0D+00, angle, dim_num )

        volume = factor1 * factor2 * r**dim_num

        if ( r .lt. h ) then
          call sphere_imp_volume_nd ( dim_num, r, volume2 )
          volume = volume2 - volume
        end if

      end if

      return
      end
      subroutine sphere_dia2imp_3d ( p1, p2, r, pc )

c*********************************************************************72
c
cc SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) )^2 = R^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), are two points which form a
c    diameter of the sphere.
c
c    Output, double precision R, the computed radius of the sphere.
c
c    Output, double precision PC(3), the computed center of the sphere.
c
      implicit none

      integer i
      double precision p1(3)
      double precision p2(3)
      double precision pc(3)
      double precision r
      double precision r8vec_norm_affine

      r = 0.5D+00 * r8vec_norm_affine ( 3, p1, p2 )

      do i = 1, 3
        pc(i) = 0.5D+00 * ( p1(i) + p2(i) )
      end do

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

      double precision bot
      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r
      double precision r8_asin
      double precision r8_atan
      double precision r8vec_norm
      double precision top
      double precision xyz1(3)
      double precision xyz2(3)

      r = r8vec_norm ( 3, xyz1 )

      lat1 = r8_asin ( xyz1(3) )
      lon1 = r8_atan ( xyz1(2), xyz1(1) )

      lat2 = r8_asin ( xyz2(3) )
      lon2 = r8_atan ( xyz2(2), xyz2(1) )

      top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2
     &    + ( cos ( lat1 ) * sin ( lat2 )
     &    -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

      top = sqrt ( top )

      bot = sin ( lat1 ) * sin ( lat2 )
     &    + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

      dist = r * atan2 ( top, bot )

      return
      end
      subroutine sphere_distance1 ( lat1, lon1, lat2, lon2, r, dist )

c*********************************************************************72
c
cc SPHERE_DISTANCE1 computes great circle distances on a sphere.
c
c  Discussion:
c
c    This computation is based on the law of cosines for spheres.
c    This formula can suffer from rounding errors when the angular
c    distances are small.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
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
c    Input, double precision LAT1, LON1, the latitude and longitude of
c    the first point.
c
c    Input, double precision LAT2, LON2, the latitude and longitude of
c    the second point.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision DIST, the great circle distance between
c    the points, measured in the same units as R.
c
      implicit none

      double precision c
      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r

      c = cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )
     &  + sin ( lat1 ) * sin ( lat2 )

      dist = r * acos ( c )

      return
      end
      subroutine sphere_distance2 ( lat1, lon1, lat2, lon2, r, dist )

c*********************************************************************72
c
cc SPHERE_DISTANCE2 computes great circle distances on a sphere.
c
c  Discussion:
c
c    This computation is written in terms of haversines, and can be more
c    accurate when measuring small angular distances.  It can be somewhat
c    inaccurate when the two points are antipodal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2009
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
c    Input, double precision LAT1, LON1, the latitude and longitude of
c    the first point.
c
c    Input, double precision LAT2, LON2, the latitude and longitude of
c    the second point.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision DIST, the great circle distance between
c    the points, measured in the same units as R.
c
      implicit none

      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r
      double precision s

      s = ( sin ( ( lat1 - lat2 ) / 2.0D+00 ) )**2
     &  + cos ( lat1 ) * cos ( lat2 )
     &  * ( sin ( ( lon1 - lon2 ) / 2.0D+00 ) )**2

      s = sqrt ( s )

      dist = 2.0D+00 * r * asin ( s )

      return
      end
      subroutine sphere_distance3 ( lat1, lon1, lat2, lon2, r, dist )

c*********************************************************************72
c
cc SPHERE_DISTANCE3 computes great circle distances on a sphere.
c
c  Discussion:
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
c    22 February 2009
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
c    Input, double precision LAT1, LON1, the latitude and longitude of
c    the first point.
c
c    Input, double precision LAT2, LON2, the latitude and longitude of
c    the second point.
c
c    Input, double precision R, the radius of the sphere.
c
c    Output, double precision DIST, the great circle distance between
c    the points, measured in the same units as R.
c
      implicit none

      double precision bot
      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r
      double precision top

      top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2
     &    + ( cos ( lat1 ) * sin ( lat2 )
     &    -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

      top = sqrt ( top )

      bot = sin ( lat1 ) * sin ( lat2 )
     &    + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

      dist = r * atan2 ( top, bot )

      return
      end
      subroutine sphere_exp_contains_point_3d ( p1, p2, p3, p4, p, 
     &  inside )

c*********************************************************************72
c
cc SPHERE_EXP_CONTAINS_POINT_3D: does an explicit sphere contain a point in 3D.
c
c  Discussion:
c
c    An explicit sphere in 3D is determined by four points,
c    which should be distinct, and not coplanar.
c
c    The computation checks the determinant of the 5 by 5 matrix:
c
c      x1  y1  z1  x1^2+y1^2+z1^2  1
c      x2  y2  z2  x2^2+y2^2+z2^2  1
c      x3  y3  z3  x3^2+y3^2+z3^2  1
c      x4  y4  z4  x4^2+y4^2+z4^2  1
c      x   y   z   x^2 +y^2 +z^2   1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), P4(3),
c    four distinct noncoplanar points on the sphere.
c
c    Input, double precision P(3), the coordinates of a point, whose
c    position relative to the sphere is desired.
c
c    Output, logical INSIDE, is TRUE if the point is in the sphere.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a(5,5)
      double precision det
      logical inside
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision r8mat_det_5d
      double precision r8vec_norm_squared
c
c  Compute the determinant.
c
      a(1,1) = p1(1)
      a(1,2) = p1(2)
      a(1,3) = p1(3)
      a(1,4) = r8vec_norm_squared ( dim_num, p1 )
      a(1,5) = 1.0D+00

      a(2,1) = p2(1)
      a(2,2) = p2(2)
      a(2,3) = p2(3)
      a(2,4) = r8vec_norm_squared ( dim_num, p2 )
      a(2,5) = 1.0D+00

      a(3,1) = p3(1)
      a(3,2) = p3(2)
      a(3,3) = p3(3)
      a(3,4) = r8vec_norm_squared ( dim_num, p3 )
      a(3,5) = 1.0D+00

      a(4,1) = p4(1)
      a(4,2) = p4(2)
      a(4,3) = p4(3)
      a(4,4) = r8vec_norm_squared ( dim_num, p4 )
      a(4,5) = 1.0D+00

      a(5,1) = p(1)
      a(5,2) = p(2)
      a(5,3) = p(3)
      a(5,4) = r8vec_norm_squared ( dim_num, p )
      a(5,5) = 1.0D+00

      det = r8mat_det_5d ( a )

      if ( det .lt. 0.0D+00 ) then
        inside = .false.
      else if ( 0.0D+00 .le. det ) then
        inside = .true.
      end if

      return
      end
      subroutine sphere_exp_point_near_3d ( p1, p2, p3, p4, p, pn )

c*********************************************************************72
c
cc SPHERE_EXP_POINT_NEAR_3D: nearest point on explicit sphere to a point in 3D.
c
c  Discussion:
c
c    An explicit sphere in 3D is determined by four points,
c    which should be distinct, and not coplanar.
c
c    If the center of the sphere is PC, and the point is P, then
c    the desired point lies at a positive distance R along the vector 
c    P-PC unless P = PC in which case any point on the sphere is "nearest".
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
c    Input, double precision P1(3), P2(3), P3(3), P4(3),
c    four distinct noncoplanar points on the sphere.
c
c    Input, double precision P(3), a point whose nearest point on the 
c    sphere is desired.
c
c    Output, double precision PN(3), the nearest point on the sphere.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision norm
      double precision p(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision pc(dim_num)
      double precision pn(dim_num)
      double precision r
      double precision r8vec_diff_norm
c
c  Find the center.
c
      call sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )
c
c  If P = PC, bail out now.
c
      norm = r8vec_diff_norm ( dim_num, p, pc )

      if ( norm .eq. 0.0D+00 ) then
        do i = 1, dim_num
          pn(i) = pc(i) + r / sqrt ( dble ( dim_num ) )
        end do
        return
      end if
c
c  Compute the nearest point.
c
      do i = 1, dim_num
        pn(i) = pc(i) + r * ( p(i) - pc(i) ) / norm
      end do

      return
      end
      subroutine sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )

c*********************************************************************72
c
cc SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
c
c  Discussion:
c
c    An explicit sphere in 3D is determined by four points,
c    which should be distinct, and not coplanar.
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
c    30 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision P1(3), P2(3), P3(3), P4(3),
c    four distinct noncoplanar points on the sphere.
c
c    Output, double precision R, PC(3), the radius and the center
c    of the sphere.  If the linear system is
c    singular, then R = -1, PC(1:3) = 0.
c
      implicit none
  
      integer dim_num
      parameter ( dim_num = 3 )

      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision pc(dim_num)
      double precision r
      double precision tetra(dim_num,4)

      tetra(1,1) = p1(1)
      tetra(2,1) = p1(2)
      tetra(3,1) = p1(3)
      tetra(1,2) = p2(1)
      tetra(2,2) = p2(2)
      tetra(3,2) = p2(3)
      tetra(1,3) = p3(1)
      tetra(2,3) = p3(2)
      tetra(3,3) = p3(3)
      tetra(1,4) = p4(1)
      tetra(2,4) = p4(2)
      tetra(3,4) = p4(3)

      call tetrahedron_circumsphere_3d ( tetra, r, pc )

      return
      end
      subroutine sphere_exp2imp_nd ( n, p, r, pc )

c*********************************************************************72
c
cc SPHERE_EXP2IMP_ND finds an N-dimensional sphere through N+1 points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision P(N,N+1), the points.
c
c    Output, double precision R, PC(N), the radius and center of the
c    sphere.
c
      implicit none

      integer n

      double precision a(n,n+1)
      integer i
      integer info
      integer j
      double precision pc(n)
      double precision r
      double precision p(n,n+1)
c
c  Set up the linear system.
c
      do j = 1, n
        do i = 1, n
          a(i,j) = p(j,i+1)
        end do
      end do

      do j = 1, n
        do i = 1, n
          a(i,j) = a(i,j) - p(j,1)
        end do
      end do

      do i = 1, n
        r = 0.0D+00
        do j = 1, n
          r = r + a(i,j)**2
        end do
        a(i,n+1) = r
      end do
c
c  Solve the linear system.
c
      call r8mat_solve ( n, 1, a, info )
c
c  If the system was singular, return a consolation prize.
c
      if ( info .ne. 0 ) then
        r = -1.0D+00
        do i = 1, n
          pc(i) = 0.0D+00
        end do
        return
      end if
c
c  Compute the radius and center.
c
      r = 0.0D+00
      do i = 1, n
        r = r + a(i,n+1)**2
      end do
      r = 0.5D+00 * sqrt ( r )
      do i = 1, n
        pc(i) = p(i,1) + 0.5D+00 * a(i,n+1)
      end do

      return
      end
      subroutine sphere_imp_area_3d ( r, area )

c*********************************************************************72
c
cc SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
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
      subroutine sphere_imp_area_nd ( dim_num, r, area )

c*********************************************************************72
c
cc SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
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
      subroutine sphere_imp_contains_point_3d ( r, pc, p, inside )

c*********************************************************************72
c
cc SPHERE_IMP_CONTAINS_POINT_3D: point in implicit sphere in 3D?
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
c    Input, double precision PC(3), the center of the sphere.
c
c    Input, double precision P(3), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if the point is inside the sphere.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      logical inside
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision r
      double precision r8vec_diff_norm

      if ( r8vec_diff_norm ( dim_num, p, pc ) .le. r ) then
        inside = .true.
      else
        inside = .false.
      end if

      return
      end
      subroutine sphere_imp_line_project_3d ( r, pc, n, p, maxpnt2, n2, 
     &  pp, theta_min, theta_max )

c*********************************************************************72
c
cc SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
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
c    01 April 2011
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
      integer dim_num
      parameter ( dim_num = 3 )
      integer n

      double precision alpha
      double precision ang3d
      double precision dot
      integer i
      integer i2
      integer j
      integer nfill
      integer n2
      double precision p(dim_num,n)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pc(dim_num)
      double precision pd(dim_num)
      double precision pp(dim_num,maxpnt2)
      double precision r
      double precision r8_acos
      double precision r8vec_diff_dot_product
      double precision r8vec_diff_norm
      logical r8vec_eq
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

      call r8vec_copy ( dim_num, pc, p1 )
      call r8vec_copy ( dim_num, pc, p2 )

      n2 = 0

      do i = 1, n

        if ( r8vec_eq ( dim_num, p(1,i), pc ) ) then

        else

          call r8vec_copy ( dim_num, p2, p1 )

          alpha = r8vec_diff_norm ( dim_num, p, pc )

          do i2 = 1, dim_num
            p2(i2) = pc(i2) + r * ( p(i2,i) - pc(i2) ) / alpha
          end do
c
c  If we haven't gotten any points yet, take this point as our start.
c
          if ( n2 .eq. 0 ) then

            n2 = n2 + 1
            do i2 = 1, dim_num
              pp(i2,n2) = p2(i2)
            end do
c
c  Compute the angular projection of P1 to P2.
c
          else if ( 1 .le. n2 ) then

            dot = r8vec_diff_dot_product ( dim_num, p1, pc, p2, pc )

            ang3d = r8_acos (  dot / ( r * r ) )
c
c  If the angle is at least THETA_MIN, (or it's the last point),
c  then we will draw a line segment.
c
            if ( theta_min .lt. abs ( ang3d ) .or. i .eq. n ) then
c
c  Now we check to see if the line segment is too long.
c
              if ( theta_max .lt. abs ( ang3d ) ) then

                nfill = int ( abs ( ang3d ) / theta_max )

                do j = 1, nfill - 1

                  do i2 = 1, dim_num
                    pd(i2) = 
     &                ( dble ( nfill - j ) * ( p1(i2) - pc(i2) ) 
     &                + dble (         j ) * ( p2(i2) - pc(i2) ) )
                  end do

                  tnorm = r8vec_norm ( dim_num, pd )

                  if ( tnorm .ne. 0.0D+00 ) then
                    do i2 = 1, dim_num
                      pd(i2) = pc(i2) + r * pd(i2) / tnorm
                    end do
                    n2 = n2 + 1
                    do i2 = 1, dim_num
                      pp(i2,n2) = pd(i2)
                    end do
                  end if

                end do

              end if
c
c  Now tack on the projection of point 2.
c
              n2 = n2 + 1
              call r8vec_copy ( dim_num, p2, pp )

            end if

          end if

        end if

      end do

      return
      end
      subroutine sphere_imp_local2xyz_3d ( r, pc, theta, phi, p )

c*********************************************************************72
c
cc SPHERE_IMP_LOCAL2XYZ_3D: local to XYZ coordinates on implicit sphere in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c    The "local" spherical coordinates of a point are two angles, THETA and PHI.
c    PHI measures the angle that the vector from the origin to the point
c    makes with the positive Z axis.  THETA measures the angle that the
c    projection of the vector onto the XY plane makes with the positive X axis.
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
c    Input, double precision PC(3), the center of the sphere.
c
c    Input, double precision THETA, PHI, the local (THETA,PHI) spherical
c    coordinates of a point on the sphere.  THETA and PHI are angles,
c    measured in radians.  Usually, 0 <= THETA < 2 * PI, and 0 <= PHI <= PI.
c
c    Output, double precision P(3), the XYZ coordinates of the point.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision p(dim_num)
      double precision pc(dim_num)
      double precision phi
      double precision r
      double precision theta

      p(1) = pc(1) + r * sin ( phi ) * cos ( theta )
      p(2) = pc(2) + r * sin ( phi ) * sin ( theta )
      p(3) = pc(3) + r * cos ( phi )

      return
      end
      subroutine sphere_imp_point_near_3d ( r, pc, p, pn )

c*********************************************************************72
c
cc SPHERE_IMP_POINT_NEAR_3D: nearest point on implicit sphere to a point in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c    If the center of the sphere is PC, and the point is P, then
c    the desired point lies at a positive distance R along the vector 
c    P-PC unless P = PC, in which case any point on the sphere is "nearest".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
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
c    Input, double precision P(3), a point whose
c    nearest point on the sphere is desired.
c
c    Output, double precision PN(3), the nearest point on the sphere.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision norm
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision pn(dim_num)
      double precision r
      double precision r8vec_diff_norm
c
c  If P = PC, bail out now.
c
      norm = r8vec_diff_norm ( dim_num, p, pc )

      if ( norm .eq. 0.0D+00 ) then
        do i = 1, dim_num
          pn(i) = pc(i) + r / sqrt ( dble ( dim_num ) )
        end do
        return
      end if
c
c  Compute the nearest point.
c
      do i = 1, dim_num
        pn(i) = pc(i) + r * ( p(i) - pc(i) ) / norm
      end do

      return
      end
      subroutine sphere_imp_point_project_3d ( r, pc, p, pp )

c*********************************************************************72
c
cc SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere in 3D.
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
c    01 April 2011
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
c    Input, double precision P(3), a point.
c
c    Output, double precision PP(3), the projected point.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision norm
      double precision p(dim_num)
      double precision pc(dim_num)
      double precision pp(dim_num)
      double precision r
      double precision r8vec_diff_norm
      logical r8vec_eq

      if ( r .eq. 0.0D+00 ) then

        call r8vec_copy ( dim_num, pc, pp )

      else if ( r8vec_eq ( dim_num, p, pc ) ) then

        do i = 1, dim_num
          pp(i) = pc(i) + r / sqrt ( dble ( dim_num ) )
        end do

      else

        norm = r8vec_diff_norm ( dim_num, p, pc )

        do i = 1, dim_num
          pp(i) = pc(i) + r * ( p(i) - pc(i) ) / norm
        end do

      end if

      return
      end
      subroutine sphere_imp_volume_3d ( r, volume )

c*********************************************************************72
c
cc SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
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
      subroutine sphere_imp_volume_nd ( dim_num, r, volume )
c
c*********************************************************************72
c
cc SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
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
      subroutine sphere_imp_zone_area_3d ( r, h1, h2, area  )

c*********************************************************************72
c
cc SPHERE_IMP_ZONE_AREA_3D computes the surface area of a spherical zone in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c    Draw any radius of the sphere and note the point P where the radius
c    intersects the sphere.  Now choose two points on the radius line, a
c    distance H1 and H2 from the point P.  Consider all the points on or within
c    the sphere whose projection onto the radius lies between these two points.
c    These points constitute the spherical zone, which can also be considered
c    the difference of two spherical caps.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H1, H2, the distances that define the 
c    thickness of the zone.  H1 and H2 must be between 0 and 2 * R.
c
c    Output, double precision AREA, the area of the spherical zone.
c
      implicit none

      double precision area
      double precision h
      double precision h1
      double precision h2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r

      h = abs ( h1 - h2 )

      if ( h .le. 0.0D+00 ) then
        area = 0.0D+00
      else if ( 2.0D+00 * r .le. h ) then
        area = 4.0D+00 * pi * r * r
      else
        area = 2.0D+00 * pi * r * h
      end if

      return
      end
      subroutine sphere_imp_zone_volume_3d ( r, h1, h2, volume )

c*********************************************************************72
c
cc SPHERE_IMP_ZONE_VOLUME_3D computes the volume of a spherical zone in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R^2
c
c    Draw any radius of the sphere and note the point P where the radius
c    intersects the sphere.  Now choose two points on the radius line, a
c    distance H1 and H2 from the point P.  Consider all the points on or within
c    the sphere whose projection onto the radius lies between these two points.
c    These points constitute the spherical zone, which can also be considered
c    the difference of two spherical caps.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision H1, H2, the distances that define the 
c    thickness of the zone.  H1 and H2 must be between 0 and 2 * R.
c
c    Output, double precision VOLUME, the volume of the spherical zone
c
      implicit none

      double precision h1
      double precision h11
      double precision h2
      double precision h22
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision volume

      h11 = min ( h1, h2 )
      h11 = max ( h11, 0.0D+00 )

      if ( 2.0D+00 * r .le. h11 ) then
        volume = 0.0D+00
        return
      end if

      h22 = max ( h1, h2 )
      h22 = min ( h22, 2.0D+00 * r )

      if ( h22 .le. 0.0D+00 ) then
        volume = 0.0D+00
        return
      end if

      volume = ( 1.0D+00 / 3.0D+00 ) * pi * ( 
     &    h22 * h22 * ( 3.0D+00 * r - h22 ) 
     &  - h11 * h11 * ( 3.0D+00 * r - h11 ) )

      return
      end
      subroutine sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 )

c*********************************************************************72
c
cc SPHERE_IMP2EXP_3D converts a sphere from implicit to explicit form in 3D.
c
c  Discussion:
c
c    An implicit sphere in 3D satisfies the equation:
c
c      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
c
c    An explicit sphere in 3D is determined by four points,
c    which should be distinct, and not coplanar.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision R, PC(3), the radius and center of the sphere.
c
c    Output, double precision P1(3), P2(3), P3(3), P4(3),
c    four distinct noncoplanar points on the sphere.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision pc(dim_num)
      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta

      theta = 0.0D+00
      phi = 0.0D+00

      p1(1) = pc(1) + r * cos ( theta ) * sin ( phi )
      p1(2) = pc(2) + r * sin ( theta ) * sin ( phi )
      p1(3) = pc(3) + r                 * cos ( phi )

      theta = 0.0D+00
      phi = 2.0D+00 * pi / 3.0D+00

      p2(1) = pc(1) + r * cos ( theta ) * sin ( phi )
      p2(2) = pc(2) + r * sin ( theta ) * sin ( phi )
      p2(3) = pc(3) + r                 * cos ( phi )

      theta = 2.0D+00 * pi / 3.0D+00
      phi = 2.0D+00 * pi / 3.0D+00

      p3(1) = pc(1) + r * cos ( theta ) * sin ( phi )
      p3(2) = pc(2) + r * sin ( theta ) * sin ( phi )
      p3(3) = pc(3) + r                 * cos ( phi )

      theta = 4.0D+00 * pi / 3.0D+00
      phi = 2.0D+00 * pi / 3.0D+00

      p4(1) = pc(1) + r * cos ( theta ) * sin ( phi )
      p4(2) = pc(2) + r * sin ( theta ) * sin ( phi )
      p4(3) = pc(3) + r                 * cos ( phi )

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
      subroutine sphere_triangle_angles_to_area ( r, a, b, c, area )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
c
c  Discussion:
c
c    A sphere centered at 0 in 3D satisfies the equation:
c
c      X*X + Y*Y + Z*Z = R*R
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.
c
c    The area formula is known as Girard's formula.
c
c    The area of a spherical triangle is:
c
c      AREA = ( A + B + C - PI ) * R*R
c
c    where A, B and C are the (surface) angles of the triangle.
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
c    Input, double precision R, the radius of the sphere.
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
      double precision r
c
c  Apply Girard's formula.
c
      area = r * r * ( a + b + c - pi )

      return
      end
      subroutine sphere_triangle_contains_point ( v1, v2, v3, p,
     &  inside )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_CONTAINS_POINT determines if a spherical triangle contains a point.
c
c  Discussion:
c
c    A sphere centered at 0 in 3D satisfies the equation:
c
c      X*X + Y*Y + Z*Z = R*R
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.  The inside of the triangle is defined by the fact
c    that the three points are listed in counterclockwise order.
c    Here "counterclockwise" is with reference to an observer standing
c    outside the sphere.
c
c    If P is a point on the sphere, we say that the spherical triangle
c    "contains" P if P is in the interior of the spherical triangle.
c    We do not actually require that P be a point on the sphere.  Instead,
c    we consider the ray defined from the origin through P, which intersects
c    the sphere.  It is essentially this point of intersection we are
c    considering.
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
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Input, double precision P(3), a point on the sphere, or the point on
c    the sphere determined by the ray from the origin through P.  P must
c    not be zero.
c
c    Output, double precision INSIDE, is positive if the spherical triangle
c    contains P, zero if P is exactly on the boundary of the triangle, and
c    negative if P is outside the triangle.
c
      implicit none

      double precision inside
      integer i
      double precision p(3)
      double precision p_direction(3)
      double precision p_norm
      double precision normal_direction(3)
      double precision normal_norm
      double precision r8_huge
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
c
c  Determine the normal vector to (V1,V2,V3), which is (V2-V1)x(V3-V13)..
c
      normal_direction(1) = ( v2(2) - v1(2) ) * ( v3(3) - v1(3) )
     &                    - ( v2(3) - v1(3) ) * ( v3(2) - v1(2) )

      normal_direction(2) = ( v2(3) - v1(3) ) * ( v3(1) - v1(1) )
     &                    - ( v2(1) - v1(1) ) * ( v3(3) - v1(3) )

      normal_direction(3) = ( v2(1) - v1(1) ) * ( v3(2) - v1(2) )
     &                    - ( v2(2) - v1(2) ) * ( v3(1) - v1(1) )

      normal_norm = r8vec_norm ( 3, normal_direction )

      if ( normal_norm .eq. 0.0D+00 ) then
        inside = - r8_huge ( )
        return
      end if

      do i = 1, 3
        normal_direction(i) = normal_direction(i) / normal_norm
      end do
c
c  Determine the length of P.
c
      p_norm = r8vec_norm ( 3, p )

      if ( p_norm .eq. 0.0D+00 ) then
        inside = r8_huge ( )
        return
      end if

      do i = 1, 3
        p_direction(i) = p_direction(i) / p_norm
      end do
c
c  INSIDE is the dot product of the normal vector to (V1,V2,V3)
c  against the unit direction vector defined by P.
c
      inside = r8vec_dot_product ( 3, normal_direction, p_direction )

      return
      end
      subroutine sphere_triangle_sides_to_angles ( r, as, bs, cs,
     &  a, b, c )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles in 3D.
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
c    Input, double precision R, the radius of the sphere.
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
      double precision r
      double precision ssu
      double precision tan_a2
      double precision tan_b2
      double precision tan_c2

      asu = as / r
      bsu = bs / r
      csu = cs / r
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
      subroutine sphere_triangle_vertices_to_angles ( r, v1, v2, v3,
     &  a, b, c )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_VERTICES_TO_ANGLES: spherical triangle angles from vertices.
c
c  Discussion:
c
c    A sphere centered at 0 in 3D satisfies the equation:
c
c      X * X + Y * Y + Z * Z = R * R
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.
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
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Output, double precision A, B, C, the angles of the spherical triangle.
c
      implicit none

      double precision a
      double precision as
      double precision b
      double precision bs
      double precision c
      double precision cs
      double precision r
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
c
c  Compute the lengths of the sides of the spherical triangle.
c
      call sphere_triangle_vertices_to_sides ( r, v1, v2, v3,
     &  as, bs, cs )
c
c  Get the spherical angles.
c
      call sphere_triangle_sides_to_angles ( r, as, bs, cs,
     &  a, b, c )

      return
      end
      subroutine sphere_triangle_vertices_to_area ( r, v1, v2, v3,
     &  area )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle in 3D.
c
c  Discussion:
c
c    A sphere centered at 0 in 3D satisfies the equation:
c
c      X * X + Y * Y + Z * Z = R * R
c
c    A spherical triangle is specified by three points on the surface
c    of the sphere.
c
c    The area formula is known as Girard's formula.
c
c    The area of a spherical triangle is:
c
c      AREA = ( A + B + C - PI ) * R*R
c
c    where A, B and C are the (surface) angles of the triangle.
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
c    Input, double precision R, the radius of the sphere.
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
      double precision r
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
c
c  Compute the lengths of the sides of the spherical triangle.
c
      call sphere_triangle_vertices_to_sides ( r, v1, v2, v3,
     &  as, bs, cs )
c
c  Get the spherical angles.
c
      call sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c )
c
c  Get the area.
c
      call sphere_triangle_angles_to_area ( r, a, b, c, area )

      return
      end
      subroutine sphere_triangle_vertices_to_centroid ( r, v1, v2, v3,
     &  vs )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_VERTICES_TO_CENTROID gets a spherical triangle centroid in 3D.
c
c  Discussion:
c
c    A sphere centered at 0 in 3D satisfies the equation:
c
c      X*X + Y*Y + Z*Z = R*R
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
c    24 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
c
c    Output, double precision VS(3), the coordinates of the "spherical
c    centroid" of the spherical triangle.
c
      implicit none

      integer i
      double precision norm
      double precision r
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
        vs(i) = r * vs(i) / norm
      end do

      return
      end
      subroutine sphere_triangle_vertices_to_orientation ( a, b, c, o )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_VERTICES_TO_ORIENTATION seeks the orientation of a spherical triangle.
c
c  Discussion:
c
c    Three points on a sphere actually compute two triangles; typically
c    we are interested in the smaller of the two.
c
c    As long as our triangle is "small", we can define an orientation
c    by comparing the direction of the centroid against the normal
c    vector (C-B) x (A-B).  If the dot product of these vectors
c    is positive, we say the triangle has positive orientation.
c
c    By using information from the triangle orientation, we can correctly
c    determine the area of a Voronoi polygon by summing up the pieces
c    of Delaunay triangles, even in the case when the Voronoi vertex
c    lies outside the Delaunay triangle.  In that case, the areas of
c    some of the Delaunay triangle pieces must be formally negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A(3), B(3), C(3), three points on a sphere.
c
c    Output, integer O, is +1 if the spherical triangle is
c    judged to have positive orientation, and -1 otherwise.
c
      implicit none

      double precision a(3)
      double precision b(3)
      double precision c(3)
      double precision cd(3)
      double precision cp(3)
      integer i
      integer o
      double precision r8vec_dot_product
      double precision v1(3)
      double precision v2(3)
c
c  Centroid.
c
      do i = 1, 3
        cd(i) = ( a(i) + b(i) + c(i) ) / 3.0D+00
      end do
c
c  Cross product ( C - B ) x ( A - B );
c
      do i = 1, 3
        v1(i) = c(i) - b(i)
        v2(i) = a(i) - b(i)
      end do

      cp(1) = v1(2) * v2(3) - v1(3) * v2(2)
      cp(2) = v1(3) * v2(1) - v1(1) * v2(3)
      cp(3) = v1(1) * v2(2) - v1(2) * v2(1)
c
c  Compare the directions.
c
      if ( r8vec_dot_product ( 3, cp, cd ) .lt. 0.0D+00 ) then
        o = - 1
      else
        o = + 1
      end if

      return
      end
      subroutine sphere_triangle_vertices_to_sides ( r, v1, v2, v3,
     &  as, bs, cs )

c*********************************************************************72
c
cc SPHERE_TRIANGLE_VERTICES_TO_SIDES computes spherical triangle sides.
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
c    24 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the radius of the sphere.
c
c    Input, double precision V1(3), V2(3), V3(3), the vertices of the spherical
c    triangle.
c
c    Output, double precision AS, BS, CS, the (geodesic) length of the sides
c    of the triangle.
c
      implicit none

      double precision as
      double precision bs
      double precision cs
      double precision r
      double precision r8_acos
      double precision r8vec_dot_product
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)

      as = r * r8_acos ( r8vec_dot_product ( 3, v2, v3 ) / r**2 )
      bs = r * r8_acos ( r8vec_dot_product ( 3, v3, v1 ) / r**2 )
      cs = r * r8_acos ( r8vec_dot_product ( 3, v1, v2 ) / r**2 )

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
c      Sphere_Unit_Area ( N ) = 2 * PI**(N/2) / Gamma ( N / 2 )
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
      subroutine sphere_unit_sample_2d ( seed, x )

c*********************************************************************72
c
cc SPHERE_UNIT_SAMPLE_2D picks a random point on the unit sphere (circle) in 2D.
c
c  Discussion:
c
c    The unit sphere in 2D satisfies:
c
c      X * X + Y * Y = 1
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
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(2), a random point on the unit circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision r8_uniform_01
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seed
      double precision u
      double precision x(dim_num)

      u = r8_uniform_01 ( seed )

      x(1) = cos ( 2.0D+00 * pi * u )
      x(2) = sin ( 2.0D+00 * pi * u )

      return
      end
      subroutine sphere_unit_sample_3d ( seed, x )

c*********************************************************************72
c
cc SPHERE_UNIT_SAMPLE_3D picks a random point on the unit sphere in 3D.
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
c    29 May 2011
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
c    Output, double precision X(3), the sample point.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision phi
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_acos
      double precision r8_uniform_01
      integer seed
      double precision theta
      double precision vdot
      double precision x(dim_num)
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
 
      phi = r8_acos ( vdot )
c
c  Pick a uniformly random rotation between 0 and 2 Pi around the
c  axis of the Z vector.
c
      theta = r8_uniform_01 ( seed )
      theta = 2.0D+00 * pi * theta

      x(1) = cos ( theta ) * sin ( phi )
      x(2) = sin ( theta ) * sin ( phi )
      x(3) = cos ( phi )

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

      double precision bot
      double precision dist
      double precision lat1
      double precision lat2
      double precision lon1
      double precision lon2
      double precision r8_asin
      double precision r8_atan
      double precision top
      double precision xyz1(3)
      double precision xyz2(3)

      lat1 = r8_asin ( xyz1(3) )
      lon1 = r8_atan ( xyz1(2), xyz1(1) )

      lat2 = r8_asin ( xyz2(3) )
      lon2 = r8_atan ( xyz2(2), xyz2(1) )

      top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2
     &    + ( cos ( lat1 ) * sin ( lat2 )
     &    -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

      top = sqrt ( top )

      bot = sin ( lat1 ) * sin ( lat2 )
     &    + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

      dist = atan2 ( top, bot )

      return
      end
      function sphere01_polygon_area ( n, lat, lon )

c*********************************************************************72
c
cc SPHERE01_POLYGON_AREA returns the area of a spherical polygon.
c
c  Discussion:
c
c    On a unit sphere, the area of a spherical polygon with N sides
c    is equal to the spherical excess:
c
c      E = sum ( interior angles ) - ( N - 2 ) * pi.
c
c    On a sphere with radius R, the area is the spherical excess multiplied
c    by R * R.
c
c    The code was revised in accordance with suggestions in Carvalho and
c    Cavalcanti.
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
c    Original C version by Robert Miller.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
c    Point in Polyhedron Testing Using Spherical Polygons,
c    in Graphics Gems V,
c    edited by Alan Paeth,
c    Academic Press, 1995,
c    ISBN: 0125434553,
c    LC: T385.G6975.
c
c    Robert Miller,
c    Computing the Area of a Spherical Polygon,
c    Graphics Gems, Volume IV, pages 132-138,
c    Edited by Paul Heckbert,
c    Academic Press, 1994, T385.G6974.
c
c    Eric Weisstein,
c    "Spherical Polygon",
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 1999.
c
c  Parameters:
c
c    Input, integer N, the number of vertices.
c
c    Input, double precision LAT[N], LON[N], the latitudes and longitudes
c    of the vertices of the spherical polygon.
c
c    Output, double precision SPHERE01_POLYGON_AREA, the area of the
c    spherical polygon, measured in spherical radians.
c
      implicit none

      integer n

      double precision a
      double precision area
      double precision b
      double precision beta1
      double precision beta2
      double precision c
      double precision cos_b1
      double precision cos_b2
      double precision excess
      double precision hav_a
      double precision haversine
      integer j
      integer k
      double precision lam
      double precision lam1
      double precision lam2
      double precision lat(n)
      double precision lon(n)
      double precision pi_half
      parameter ( pi_half = 1.5707963267948966192313D+00 )
      double precision s
      double precision sphere01_polygon_area
      double precision t

      area = 0.0D+00

      do j = 1, n + 1

        if ( j .eq. 1 ) then
          lam1 = lon(j)
          beta1 = lat(j)
          lam2 = lon(j+1)
          beta2 = lat(j+1)
          cos_b1 = cos ( beta1 )
          cos_b2 = cos ( beta2 )
        else
          k = mod ( j + 1, n + 1 )
          lam1 = lam2
          beta1 = beta2
          lam2 = lon(k)
          beta2 = lat(k)
          cos_b1 = cos_b2
          cos_b2 = cos ( beta2 )
        end if

        if ( lam1 .ne. lam2 ) then

          hav_a = haversine ( beta2 - beta1 )
     &      + cos_b1 * cos_b2 * haversine ( lam2 - lam1 )
          a = 2.0D+00 * asin ( sqrt ( hav_a ) )

          b = pi_half - beta2
          c = pi_half - beta1
          s = 0.5D+00 * ( a + b + c )
c
c  Given the three sides of a spherical triangle, we can use a formula
c  to find the spherical excess.
c
          t = tan ( s / 2.0D+00 ) * tan ( ( s - a ) / 2.0D+00 )
     &      * tan ( ( s - b ) / 2.0D+00 ) * tan ( ( s - c ) / 2.0D+00 )

          excess = abs ( 4.0D+00 * atan ( sqrt ( abs ( t ) ) ) )

          if ( lam1 .lt. lam2 ) then
            lam = lam2 - lam1
          else
            lam = lam2 - lam1 + 4.0D+00 * pi_half
          end if

          if ( 2.0D+00 * pi_half .lt. lam ) then
            excess = -excess
          end if

          area = area + excess

        end if

      end do

      sphere01_polygon_area = abs ( area )

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
      subroutine sphere01_triangle_vertices_to_angles ( v1, v2, v3,
     &  a, b, c )

c*********************************************************************72
c
cc SPHERE01_TRIANGLE_VERTICES_TO_ANGLES: spherical triangle angles from vertices.
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
c    Output, double precision A, B, C, the angles of the spherical triangle.
c
      implicit none

      double precision a
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

      double precision as
      double precision bs
      double precision cs
      double precision r8_acos
      double precision r8vec_dot_product
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)

      as = r8_acos ( r8vec_dot_product ( 3, v2, v3 ) )
      bs = r8_acos ( r8vec_dot_product ( 3, v3, v1 ) )
      cs = r8_acos ( r8vec_dot_product ( 3, v1, v2 ) )

      return
      end
      subroutine super_ellipse_points_2d ( pc, r1, r2, expo, psi, n, p )

c*********************************************************************72
c
cc SUPER_ELLIPSE_POINTS_2D returns N points on a tilted superellipse in 2D.
c
c  Discussion:
c
c    The points are "equally spaced" in the angular sense.  They are
c    not equally spaced along the perimeter.
c
c    The parametric formula of the (untilted) superellipse is:
c
c      X = R1 * cos**EXPO ( THETA )
c      Y = R2 * sin**EXPO ( THETA )
c
c    An implicit form of the (untilted) superellipse is:
c
c      (X/R1)**(2/EXPO) + (Y/R2)**(2/EXPO) = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Martin Gardner,
c    The Mathematical Carnival,
c    Knopf, 1975, pages 240-254.
c
c  Parameters:
c
c    Input, double precision PC(2), the center of the superellipse.
c
c    Input, double precision R1, R2, the "radius" of the superellipse
c    in the major and minor axis directions.  A circle has these values equal.
c
c    Input, double precision EXPO, the exponent of the superellipse.
c    0 = a rectangle;
c    between 0 and 1, a "rounded" rectangle;
c    1.0 = an ellipse;
c    2.0 = a diamond;
c    > 2.0 a pinched shape.
c
c    Input, double precision PSI, the angle that the major axis of the
c    superellipse makes with the X axis.  A value of 0.0 means that the
c    major and minor axes of the superellipse will be the X and Y
c    coordinate axes.
c
c    Input, integer N, the number of points desired.  N must
c    be at least 1.
c
c    Output, double precision P(2,N), the coordinates of points
c    on the superellipse.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 2 )

      double precision act
      double precision ast
      integer i
      double precision expo
      double precision p(dim_num,n)
      double precision pc(dim_num)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision psi
      double precision r1
      double precision r2
      double precision sct
      double precision sst
      double precision theta

      do i = 1, n

        theta = ( 2.0D+00 * pi * dble ( i - 1 ) ) / dble ( n )

        act = abs ( cos ( theta ) )
        sct = sign ( 1.0D+00, cos ( theta ) )
        ast = abs ( sin ( theta ) )
        sst = sign ( 1.0D+00, sin ( theta ) )

        p(1,i) = pc(1) + r1 * cos ( psi ) * sct * ( act )**expo
     &                 - r2 * sin ( psi ) * sst * ( ast )**expo

        p(2,i) = pc(2) + r1 * sin ( psi ) * sct * ( act )**expo
     &                 + r2 * cos ( psi ) * sst * ( ast )**expo

      end do

      return
      end
      function tan_deg ( angle_deg )

c*********************************************************************72
c
cc TAN_DEG returns the tangent of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE_DEG, the angle, in degrees.
c
c    Output, double precision TAN_DEG, the tangent of the angle.
c
      implicit none

      double precision angle_deg
      double precision angle_rad
      double precision degrees_to_radians
      parameter ( degrees_to_radians =
     &  3.141592653589793D+00 / 180.0D+00 )
      double precision tan_deg

      angle_rad = degrees_to_radians * angle_deg
      tan_deg  = sin ( angle_rad ) / cos ( angle_rad )

      return
      end
      subroutine tetrahedron_barycentric_3d ( tetra, p, c )

c*********************************************************************72
c
cc TETRAHEDRON_BARYCENTRIC_3D: barycentric coordinates of a point in 3D.
c
c  Discussion:
c
c    The barycentric coordinates of a point P with respect to
c    a tetrahedron are a set of four values C(1:4), each associated
c    with a vertex of the tetrahedron.  The values must sum to 1.
c    If all the values are between 0 and 1, the point is contained
c    within the tetrahedron.
c
c    The barycentric coordinate of point P related to vertex A can be
c    interpreted as the ratio of the volume of the tetrahedron with 
c    vertex A replaced by vertex P to the volume of the original 
c    tetrahedron.
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
c    Input, double precision TETRA(3,4) the tetrahedron vertices.
c
c    Input, double precision P(3), the point to be checked.
c
c    Output, double precision C(4), the barycentric coordinates of P with
c    respect to the tetrahedron.
c
      implicit none

      integer dim_num 
      parameter ( dim_num = 3 )
      integer rhs_num
      parameter ( rhs_num = 1 )

      double precision a(dim_num,dim_num+rhs_num)
      double precision c(dim_num+1)
      integer i
      integer info
      integer j
      double precision p(dim_num)
      double precision tetra(dim_num,4)
c
c  Set up the linear system
c
c    ( X2-X1  X3-X1  X4-X1 ) C2    X - X1
c    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C3  = Y - Y1
c    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C4    Z - Z1
c
c  which is satisfied by the barycentric coordinates of P.
c
      do i = 1, dim_num
        do j = 1, 3
          a(i,j) = tetra(i,j+1) - tetra(i,1)
        end do
        a(i,4) = p(i) - tetra(i,1)
      end do
c
c  Solve the linear system.
c
      call r8mat_solve ( dim_num, rhs_num, a, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TETRAHEDRON_BARYCENTRIC_3D - Fatal error!'
        write ( *, '(a)' ) '  The linear system is singular.'
        write ( *, '(a)' ) 
     &    '  The input data does not form a proper tetrahedron.'
        stop
      end if

      do i = 1, dim_num
        c(i+1) = a(i,4)
      end do

      c(1) = 1.0D+00 - c(2) - c(3) - c(4)

      return
      end
      subroutine tetrahedron_centroid_3d ( tetra, centroid )

c*********************************************************************72
c
cc TETRAHEDRON_CENTROID_3D computes the centroid of a tetrahedron in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4) the tetrahedron vertices.
c
c    Output, double precision CENTROID(3), the coordinates of the centroid.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision centroid(dim_num)
      integer i
      integer j
      double precision tetra(dim_num,4)

      do i = 1, dim_num
        centroid(i) = 0.0D+00
        do j = 1, 4
          centroid(i) = centroid(i) + tetra(i,j)
        end do
        centroid(i) = centroid(i) / 4.0D+00
      end do

      return
      end
      subroutine tetrahedron_circumsphere_3d ( tetra, r, pc )

c*********************************************************************72
c
cc TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
c
c  Discussion:
c
c    The circumsphere, or circumscribed sphere, of a tetrahedron is the 
c    sphere that passes through the four vertices.  The circumsphere is
c    not necessarily the smallest sphere that contains the tetrahedron.
c
c    Surprisingly, the diameter of the sphere can be found by solving
c    a 3 by 3 linear system.  This is because the vectors P2 - P1,
c    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
c    right triangle with the diameter through P1.  Hence, the dot product of
c    P2 - P1 with that diameter is equal to the square of the length
c    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
c    the diameter vector originating at P1, and hence the radius and
c    center.
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
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision TETRA(3,4) the tetrahedron vertices.
c
c    Output, double precision R, PC(3), the center of the
c    circumscribed sphere, and its radius.  If the linear system is
c    singular, then R = -1, PC(1:3) = 0.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer rhs_num
      parameter ( rhs_num = 1 )

      double precision a(dim_num,dim_num+rhs_num)
      integer i
      integer info
      integer j
      double precision pc(dim_num)
      double precision r
      double precision tetra(dim_num,4)
c
c  Set up the linear system.
c
      do i = 1, dim_num
        do j = 1, dim_num
          a(i,j) = tetra(j,i+1) - tetra(j,1)
        end do
      end do

      do i = 1, 3
        a(i,4) = 0.0D+00
        do j = 1, 3
          a(i,4) = a(i,4) + a(i,j)**2
        end do
      end do
c
c  Solve the linear system.
c
      call r8mat_solve ( dim_num, rhs_num, a, info )
c
c  If the system was singular, return a consolation prize.
c
      if ( info .ne. 0 ) then
        r = -1.0D+00
        do i = 1, dim_num
          pc(i) = 0.0D+00
        end do
        return
      end if
c
c  Compute the radius and center.
c
      r = 0.0D+00
      do i = 1, dim_num
        r = r + a(i,4)**2
      end do
      r = 0.5D+00 * sqrt ( r )

      do i = 1, dim_num
        pc(i) = tetra(i,1) + 0.5D+00 * a(i,4)
      end do

      return
      end
      subroutine tetrahedron_contains_point_3d ( tetra, p, inside )

c*********************************************************************72
c
cc TETRAHEDRON_CONTAINS_POINT_3D finds if a point is inside a tetrahedron in 3D.
c
c  Discussion:
c
c    We compute the barycentric coordinates C(1:4) of the point, with respect
c    to the tetrahedron.  The point is inside the tetrahedron if and only
c    if each coordinate is nonnegative, and their sum is no greater than 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4), the tetrahedron vertices.
c
c    Input, double precision P(3), the point to be checked.
c
c    Output, logical INSIDE, is TRUE if P is inside the tetrahedron.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision c(dim_num+1)
      integer i
      logical inside
      double precision p(dim_num)
      double precision tetra(dim_num,4)

      call tetrahedron_barycentric_3d ( tetra, p, c )
c
c  If the point is in the tetrahedron, its barycentric coordinates
c  must be nonnegative.
c
      inside = .true.
      do i = 1, dim_num + 1
        if ( c(i) .lt. 0.0D+00 ) then
          inside = .false.
          return
        end if
      end do

      return
      end
      subroutine tetrahedron_dihedral_angles_3d ( tetra, angle )

c*********************************************************************72
c
cc TETRAHEDRON_DIHEDRAL_ANGLES_3D computes dihedral angles of a tetrahedron.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4), the vertices of the tetrahedron,
c    which can be labeled as A, B, C and D.
c
c    Output, double precision ANGLE(6), the dihedral angles along the
c    axes AB, AC, AD, BC, BD and CD, respectively.
c
      implicit none

      double precision ab(3)
      double precision abc_normal(3)
      double precision abd_normal(3)
      double precision ac(3)
      double precision acd_normal(3)
      double precision ad(3)
      double precision angle(6)
      double precision bc(3)
      double precision bcd_normal(3)
      double precision bd(3)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tetra(3,4)

      do i = 1, 3
        ab(i) = tetra(i,2) - tetra(i,1)
        ac(i) = tetra(i,3) - tetra(i,1)
        ad(i) = tetra(i,4) - tetra(i,1)
        bc(i) = tetra(i,3) - tetra(i,2)
        bd(i) = tetra(i,4) - tetra(i,2)
      end do

      call r8vec_cross_product_3d ( ac, ab, abc_normal )
      call r8vec_cross_product_3d ( ab, ad, abd_normal )
      call r8vec_cross_product_3d ( ad, ac, acd_normal )
      call r8vec_cross_product_3d ( bc, bd, bcd_normal )

      call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
      call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
      call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
      call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
      call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
      call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

      do i = 1, 6
        angle(i) = pi - angle(i)
      end do

      return
      end
      subroutine tetrahedron_edge_length_3d ( tetra, edge_length )

c*********************************************************************72
c
cc TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
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
c    Input, double precision TETRA(3,4), the tetrahedron vertices.
c
c    Output, double precision EDGE_LENGTH(6), the length of the edges.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision r8vec_norm
      double precision edge_length(6)
      integer i
      integer j1
      integer j2
      integer k
      double precision tetra(dim_num,4)

      k = 0
      do j1 = 1, 3
        do j2 = j1+1, 4
          k = k + 1
          edge_length(k) = r8vec_norm ( dim_num, 
     &      tetra(1,j2) - tetra(1,j1) )
        end do
      end do

      return
      end
      subroutine tetrahedron_face_angles_3d ( tetra, angles )

c*********************************************************************72
c
cc TETRAHEDRON_FACE_ANGLES_3D returns the 12 face angles of a tetrahedron 3D.
c
c  Discussion:
c
c    The tetrahedron has 4 triangular faces.  This routine computes the
c    3 planar angles associated with each face.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4) the tetrahedron vertices.
c
c    Output, double precision ANGLES(3,4), the face angles.
c
      implicit none

      double precision angles(3,4)
      integer i
      double precision tri(3,3)
      double precision tetra(3,4)
c
c  Face 123
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,2)
        tri(i,3) = tetra(i,3)
      end do
      call triangle_angles_3d ( tri, angles(1,1) )
c
c  Face 124
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,2)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_angles_3d ( tri, angles(1,2) )
c
c  Face 134
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,3)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_angles_3d ( tri, angles(1,3) )
c
c  Face 234
c
      do i = 1, 3
        tri(i,1) = tetra(i,2)
        tri(i,2) = tetra(i,3)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_angles_3d ( tri, angles(1,4) )

      return
      end
      subroutine tetrahedron_face_areas_3d ( tetra, areas )

c*********************************************************************72
c
cc TETRAHEDRON_FACE_AREAS_3D returns the 4 face areas of a tetrahedron 3D.
c
c  Discussion:
c
c    The tetrahedron has 4 triangular faces.  This routine computes the
c    area associated with each face.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4) the tetrahedron vertices.
c
c    Output, double precision AREAS(4), the face areas.
c
      implicit none

      double precision areas(4)
      integer i
      double precision tri(3,3)
      double precision tetra(3,4)
c
c  Face 123
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,2)
        tri(i,3) = tetra(i,3)
      end do
      call triangle_area_3d ( tri, areas(1) )
c
c  Face 124
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,2)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_area_3d ( tri, areas(2) )
c
c  Face 134
c
      do i = 1, 3
        tri(i,1) = tetra(i,1)
        tri(i,2) = tetra(i,3)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_area_3d ( tri, areas(3) )
c
c  Face 234
c
      do i = 1, 3
        tri(i,1) = tetra(i,2)
        tri(i,2) = tetra(i,3)
        tri(i,3) = tetra(i,4)
      end do
      call triangle_area_3d ( tri, areas(4) )

      return
      end
      subroutine tetrahedron_lattice_layer_point_next ( c, v, more )

c*********************************************************************72
c
cc TETRAHEDRON_LATTICE_LAYER_POINT_NEXT: next tetrahedron lattice layer point.
c
c  Discussion:
c
c    The tetrahedron lattice layer L is bounded by the lines
c
c      0 .le. X,
c      0 .le. Y,
c      0 .le. Z,
c      L - 1 < X / C(1) + Y / C(2) + Z/C(3) .le. L.
c
c    In particular, layer L = 0 always contains the single point (0,0).
c
c    This function returns, one at a time, the points that lie within
c    a given tetrahedron lattice layer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer C(4), coefficients defining the
c    lattice layer in entries 1 to 3, and the laver index in C(4).
c    The coefficients should be positive, and C(4) must be nonnegative.
c
c    Input/output, integer V(3).  On first call for a given layer,
c    the input value of V is not important.  On a repeated call for the same
c    layer, the input value of V should be the output value from the previous
c    call.  On output, V contains the next lattice layer point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given layer.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if the returned value V is a new point.
c    If the output value is FALSE, then no more points were found,
c    and V was reset to 0, and the lattice layer has been exhausted.
c
      implicit none

      integer c(4)
      integer c1n
      integer i4vec_lcm
      integer lhs
      logical more
      integer n
      parameter ( n = 3 )
      integer rhs1
      integer rhs2
      integer v(3)
c
c  Treat layer C(4) = 0 specially.
c
      if ( c(4) .eq. 0 ) then
        if ( .not. more ) then
          v(1) = 0
          v(2) = 0
          v(3) = 0
          more = .true.
        else
          more = .false.
        end if
        return
      end if
c
c  Compute the first point.
c
      if ( .not. more ) then

        v(1) = ( c(4) - 1 ) * c(1) + 1
        v(2) = 0
        v(3) = 0
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )

        rhs1 = c1n * ( c(4) - 1 )
        rhs2 = c1n *   c(4)
c
c  Can we simply increase X?
c
        v(1) = v(1) + 1

        lhs = ( c1n / c(1) ) * v(1)
     &      + ( c1n / c(2) ) * v(2)
     &      + ( c1n / c(3) ) * v(3)

        if ( lhs .le. rhs2 ) then
c
c  No.  Increase Y, and set X so we just exceed RHS1...if possible.
c
        else

          v(2) = v(2) + 1

          v(1) = ( c(1) * ( rhs1 - ( c1n / c(2) ) * v(2)
     &                           - ( c1n / c(3) ) * v(3) ) ) / c1n
          v(1) = max ( v(1), 0 )

          lhs = ( c1n / c(1) ) * v(1)
     &        + ( c1n / c(2) ) * v(2)
     &        + ( c1n / c(3) ) * v(3)

          if ( lhs .le. rhs1 ) then
            v(1) = v(1) + 1
            lhs = lhs + c1n / c(1)
          end if
c
c  We have increased Y by 1.  Have we stayed below the upper bound?
c
          if ( lhs .le. rhs2 ) then

          else
c
c  No.  Increase Z, and set X so we just exceed RHS1...if possible.
c
            v(3) = v(3) + 1
            v(2) = 0
            v(1) = ( c(1) * ( rhs1 - ( c1n / c(2) ) * v(2)
     &                             - ( c1n / c(3) ) * v(3) ) ) / c1n
            v(1) = max ( v(1), 0 )

            lhs = ( c1n / c(1) ) * v(1)
     &          + ( c1n / c(2) ) * v(2)
     &          + ( c1n / c(3) ) * v(3)

            if ( lhs .le. rhs1 ) then
              v(1) = v(1) + 1
              lhs = lhs + c1n / c(1)
            end if

            if ( lhs .le. rhs2 ) then

            else
              more = .false.
              v(1) = 0
              v(2) = 0
              v(3) = 0
            end if

          end if
        end if
      end if

      return
      end
      subroutine tetrahedron_lattice_point_next ( c, v, more )

c*********************************************************************72
c
cc TETRAHEDRON_LATTICE_POINT_NEXT returns the next tetrahedron lattice point.
c
c  Discussion:
c
c    The lattice tetrahedron is defined by the vertices:
c
c      (0,0,0), (C(4)/C(1),0,0), (0,C(4)/C(2),0) and (0,0,C(4)/C(3))
c
c    The lattice tetrahedron is bounded by the lines
c
c      0 <= X,
c      0 <= Y
c      0 <= Z,
c      X / C(1) + Y / C(2) + Z / C(3) <= C(4)
c
c    Lattice points are listed one at a time, starting at the origin,
c    with X increasing first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer C(4), coefficients defining the
c    lattice tetrahedron.  These should be positive.
c
c    Input/output, integer V(3).  On first call, the input
c    value is not important.  On a repeated call, the input value should
c    be the output value from the previous call.  On output, V contains
c    the next lattice point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given tetrahedron.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if not only is the returned value V a lattice point,
c    but the routine can be called again for another lattice point.
c    If the output value is FALSE, then no more lattice points were found,
c    and V was reset to 0, and the routine should not be called further
c    for this tetrahedron.
c
      implicit none

      integer c(4)
      integer c1n
      integer i4vec_lcm
      integer lhs
      logical more
      integer n
      parameter ( n = 3 )
      integer rhs
      integer v(3)

      if ( .not. more ) then

        v(1) = 0
        v(2) = 0
        v(3) = 0
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )
        rhs = c1n * c(4)

        lhs =        c(2) * c(3) * v(1)
     &      + c(1) *        c(3) * v(2)
     &      + c(1) * c(2)        * v(3)

        if ( lhs +  c1n / c(1) .le. rhs ) then

          v(1) = v(1) + 1

        else

          lhs = lhs - c1n * v(1) / c(1)
          v(1) = 0

          if ( lhs + c1n / c(2) .le. rhs ) then

            v(2) = v(2) + 1

          else

            lhs = lhs - c1n * v(2) / c(2)
            v(2) = 0

            if ( lhs + c1n / c(3) .le. rhs ) then

              v(3) = v(3) + 1

            else

              v(3) = 0
              more = .false.

            end if

          end if

        end if

      end if

      return
      end
      subroutine tetrahedron_quality1_3d ( tetra, quality )

c*********************************************************************72
c
cc TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
c
c  Discussion:
c
c    The quality of a tetrahedron is 3 times the ratio of the radius of 
c    the inscribed sphere divided by that of the circumscribed sphere.  
c
c    An equilateral tetrahredron achieves the maximum possible quality of 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 August 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4), the tetrahedron vertices.
c
c    Output, double precision QUALITY, the quality of the tetrahedron.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision pc(dim_num)
      double precision quality
      double precision r_in
      double precision r_out
      double precision tetra(dim_num,4)

      call tetrahedron_circumsphere_3d ( tetra, r_out, pc )

      call tetrahedron_insphere_3d ( tetra, r_in, pc )

      quality = 3.0D+00 * r_in / r_out

      return
      end
      subroutine tetrahedron_rhombic_size_3d ( point_num, edge_num,
     &  face_num, face_order_max )

c*********************************************************************72
c
cc TETRAHEDRON_RHOMBIC_SIZE_3D gives "sizes" for a rhombic tetrahedron in 3D.
c
c  Discussion:
c
c    Call this routine first, in order to learn the required dimensions
c    of arrays to be set up by TETRAHEDRON_RHOMBIC_SHAPE_3D.
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
c    Output, integer POINT_NUM, the number of vertices.
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

      point_num = 10
      edge_num = 6
      face_num = 4
      face_order_max = 6

      return
      end
      subroutine tetrahedron_size_3d ( point_num, edge_num, face_num,
     &  face_order_max )

c*********************************************************************72
c
cc TETRAHEDRON_SIZE_3D gives "sizes" for a tetrahedron in 3D.
c
c  Discussion:
c
c    Call this routine first, in order to learn the required dimensions
c    of arrays to be set up by TETRAHEDRON_SHAPE_3D.
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
c    Output, integer POINT_NUM, the number of vertices.
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

      point_num = 4
      edge_num = 6
      face_num = 4
      face_order_max = 3

      return
      end
      subroutine tetrahedron_solid_angles_3d ( tetra, angle )

c*********************************************************************72
c
cc TETRAHEDRON_SOLID_ANGLES_3D computes solid angles of a tetrahedron.
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
c    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
c
c    Output, double precision ANGLE(4), the solid angles.
c
      implicit none

      double precision angle(4)
      double precision dihedral_angles(6)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tetra(3,4)

      call tetrahedron_dihedral_angles_3d ( tetra, dihedral_angles )

      angle(1) = dihedral_angles(1) + dihedral_angles(2)
     &         + dihedral_angles(3) - pi
      angle(2) = dihedral_angles(1) + dihedral_angles(4)
     &         + dihedral_angles(5) - pi
      angle(3) = dihedral_angles(2) + dihedral_angles(4)
     &         + dihedral_angles(6) - pi
      angle(4) = dihedral_angles(3) + dihedral_angles(5)
     &         + dihedral_angles(6) - pi

      return
      end
      subroutine tetrahedron_unit_lattice_point_num_3d ( s, n )

c*********************************************************************72
c
cc TETRAHEDRON_UNIT_LATTICE_POINT_NUM_3D: count lattice points.
c
c  Discussion:
c
c    The tetrahedron is assumed to be the unit tetrahedron:
c
c    ( (0,0,0), (1,0,0), (0,1,0), (0,0,1) )
c
c    or a copy of this tetrahedron scaled by an integer S:
c
c    ( (0,0,0), (S,0,0), (0,S,0), (0,0,S) ).
c
c    The routine returns the number of integer lattice points that appear
c    inside the tetrahedron, or on its faces, edges or vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Matthias Beck, Sinai Robins,
c    Computing the Continuous Discretely,
c    Springer, 2006,
c    ISBN13: 978-0387291390,
c    LC: QA640.7.B43.
c
c  Parameters:
c
c    Input, integer S, the scale factor.
c
c    Output, integer N, the number of lattice points.
c
      implicit  none

      integer n
      integer s

      n = ( ( s + 3 ) * ( s + 2 ) * ( s + 1 ) ) / 6

      return
      end
      subroutine tetrahedron_volume_3d ( tetra, volume )

c*********************************************************************72
c
cc TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
c
c    Output, double precision VOLUME, the volume of the tetrahedron.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a(4,4)
      integer i
      integer j
      double precision r8mat_det_4d
      double precision tetra(dim_num,4)
      double precision volume

      do j = 1, 4
        do i = 1, dim_num
          a(i,j) = tetra(i,j)
        end do
        a(4,j) = 1.0D+00
      end do

      volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

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
      subroutine tmat_init ( a )

c*********************************************************************72
c
cc TMAT_INIT initializes the geometric transformation matrix.
c
c  Discussion:
c
c    The geometric transformation matrix can be thought of as a 4 by 4
c    matrix "A" having components:
c
c      r11 r12 r13 t1
c      r21 r22 r23 t2
c      r31 r32 r33 t3
c        0   0   0  1
c
c    This matrix encodes the rotations, scalings and translations that
c    are applied to graphical objects.
c
c    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as
c    PH = (x,y,z,1).  Then to apply the transformations encoded in A to
c    the point P, we simply compute A * PH.
c
c    Individual transformations, such as a scaling, can be represented
c    by simple versions of the transformation matrix.  If the matrix
c    A represents the current set of transformations, and we wish to
c    apply a new transformation B, then the original points are
c    transformed twice:  B * ( A * PH ).  The new transformation B can
c    be combined with the original one A, to give a single matrix C that
c    encodes both transformations: C = B * A.
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
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the geometric transformation matrix.
c
      implicit none

      double precision a(4,4)
      integer i
      integer j

      do i = 1, 4
        do j = 1, 4
          if ( i .eq. j ) then
            a(i,j) = 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine tmat_mxm ( a, b, c )

c*********************************************************************72
c
cc TMAT_MXM multiplies two geometric transformation matrices.
c
c  Discussion:
c
c    The product is accumulated in a temporary array, and then assigned
c    to the result.  Therefore, it is legal for any two, or all three,
c    of the arguments to share memory.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the first geometric transformation matrix.
c
c    Input, double precision B(4,4), the second geometric transformation
c    matrix.
c
c    Output, double precision C(4,4), the product A * B.
c
      implicit none

      double precision a(4,4)
      double precision b(4,4)
      double precision c(4,4)

      call r8mat_mm ( 4, 4, 4, a, b, c )

      return
      end
      subroutine tmat_mxp ( a, x, y )

c*********************************************************************72
c
cc TMAT_MXP multiplies a geometric transformation matrix times a point.
c
c  Discussion:
c
c    The matrix will normally have the form
c
c      xx xy xz tx
c      yx yy yz ty
c      zx zy zz tz
c       0  0  0  1
c
c    where the 3x3 initial block controls rotations and scalings,
c    and the values [ tx, ty, tz ] implement a translation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the geometric transformation matrix.
c
c    Input, double precision X(3), the point to be multiplied.  The fourth
c    component of X is implicitly assigned the value of 1.
c
c    Output, double precision Y(3), the result of A*X.  
c
      implicit none

      double precision a(4,4)
      integer i
      integer j
      double precision x(3)
      double precision y(3)

      do i = 1, 3
        y(i) = a(i,4)
        do j = 1, 3
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine tmat_mxp2 ( a, n, x, y )

c*********************************************************************72
c
cc TMAT_MXP2 multiplies a geometric transformation matrix times N points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the geometric transformation matrix.
c
c    Input, integer N, the number of points to be multiplied.
c
c    Input, double precision X(3,N), the points to be multiplied.
c
c    Output, double precision Y(3,N), the transformed points.  
c
      implicit none

      integer n

      double precision a(4,4)
      integer i
      integer j
      integer k
      double precision x(3,n)
      double precision y(3,n)

      do k = 1, n

        do i = 1, 3
          y(i,k) = a(i,4)
          do j = 1, 3
            y(i,k) = y(i,k) + a(i,j) * x(j,k)
          end do
        end do

      end do

      return
      end
      subroutine tmat_mxv ( a, x, y )

c*********************************************************************72
c
cc TMAT_MXV multiplies a geometric transformation matrix times a vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the geometric transformation matrix.
c
c    Input, double precision X(3), the vector to be multiplied.  The fourth
c    component of X is implicitly assigned the value of 1.
c
c    Output, double precision Y(3), the result of A*X.  The product is
c    accumulated in a temporary vector, and then assigned to the result. 
c    Therefore, it is legal for X and Y to share memory.
c
      implicit none

      double precision a(4,4)
      integer i
      integer j
      double precision x(3)
      double precision y(3)

      do i = 1, 3
        y(i) = a(i,4)
        do j = 1, 3
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine tmat_rot_axis ( a, angle, axis, b )

c*********************************************************************72
c
cc TMAT_ROT_AXIS: coordinate axis rotation to geometric transformation matrix.
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
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the current geometric transformation
c    matrix.
c
c    Input, double precision ANGLE, the angle, in degrees, of the rotation.
c
c    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
c    axis about which the rotation occurs.
c
c    Output, double precision B(4,4), the modified geometric 
c    transformation matrix.
c    A and B may share the same memory.
c
      implicit none

      double precision a(4,4)
      double precision angle
      double precision angle_rad
      character axis
      double precision b(4,4)
      double precision c(4,4)
      double precision degrees_to_radians

      angle_rad = degrees_to_radians ( angle )

      call tmat_init ( c )

      if ( axis .eq. 'X' .or. axis .eq. 'x' ) then
        c(2,2) =   cos ( angle_rad )
        c(2,3) = - sin ( angle_rad )
        c(3,2) =   sin ( angle_rad )
        c(3,3) =   cos ( angle_rad )
      else if ( axis .eq. 'Y' .or. axis .eq. 'y' ) then
        c(1,1) =   cos ( angle_rad )
        c(1,3) =   sin ( angle_rad )
        c(3,1) = - sin ( angle_rad )
        c(3,3) =   cos ( angle_rad )
      else if ( axis .eq. 'Z' .or. axis .eq. 'z' ) then
        c(1,1) =   cos ( angle_rad )
        c(1,2) = - sin ( angle_rad )
        c(2,1) =   sin ( angle_rad )
        c(2,2) =   cos ( angle_rad )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TMAT_ROT_AXIS - Fatal error!'
        write ( *, '(a)' ) '  Illegal rotation axis: ' // axis
        write ( *, '(a)' ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
        stop
      end if

      call r8mat_mm ( 4, 4, 4, c, a, b )

      return
      end
      subroutine tmat_rot_vector ( a, angle, axis, b )

c*********************************************************************72
c
cc TMAT_ROT_VECTOR: arbitrary axis rotation to geometric transformation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the current geometric transformation
c    matrix.
c
c    Input, double precision ANGLE, the angle, in degrees, of the rotation.
c
c    Input, double precision AXIS(3), the axis vector about which 
c    rotation occurs.  AXIS may not be the zero vector.
c
c    Output, double precision B(4,4), the modified geometric 
c    transformation matrix.
c    A and B may share the same memory.
c
      implicit none

      double precision a(4,4)
      double precision angle
      double precision angle_rad
      double precision axis(3)
      double precision b(4,4)
      double precision c(4,4)
      double precision ca
      double precision degrees_to_radians
      double precision norm
      double precision sa
      double precision v1
      double precision v2
      double precision v3

      v1 = axis(1)
      v2 = axis(2)
      v3 = axis(3)

      norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

      if ( norm .eq. 0.0D+00 ) then
        return
      end if

      v1 = v1 / norm
      v2 = v2 / norm
      v3 = v3 / norm

      angle_rad = degrees_to_radians ( angle )
      ca = cos ( angle_rad )
      sa = sin ( angle_rad )

      call tmat_init ( c )

      c(1,1) =                    v1 * v1 + ca * ( 1.0D+00 - v1 * v1 )
      c(1,2) = ( 1.0D+00 - ca ) * v1 * v2 - sa * v3
      c(1,3) = ( 1.0D+00 - ca ) * v1 * v3 + sa * v2

      c(2,1) = ( 1.0D+00 - ca ) * v2 * v1 + sa * v3
      c(2,2) =                    v2 * v2 + ca * ( 1.0D+00 - v2 * v2 )
      c(2,3) = ( 1.0D+00 - ca ) * v2 * v3 - sa * v1

      c(3,1) = ( 1.0D+00 - ca ) * v3 * v1 - sa * v2
      c(3,2) = ( 1.0D+00 - ca ) * v3 * v2 + sa * v1
      c(3,3) =                    v3 * v3 + ca * ( 1.0D+00 - v3 * v3 )

      call r8mat_mm ( 4, 4, 4, c, a, b )

      return
      end
      subroutine tmat_scale ( a, s, b )

c*********************************************************************72
c
cc TMAT_SCALE applies a scaling to the geometric transformation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the current geometric transformation
c    matrix.
c
c    Input, double precision S(3), the scalings to be applied to the 
c    X, Y and Z coordinates.
c
c    Output, double precision B(4,4), the modified geometric transformation
c    matrix.  A and B may share the same memory.
c
      implicit none

      double precision a(4,4)
      double precision b(4,4)
      double precision c(4,4)
      double precision s(3)

      call tmat_init ( c )

      c(1,1) = s(1)
      c(2,2) = s(2)
      c(3,3) = s(3)

      call r8mat_mm ( 4, 4, 4, c, a, b )

      return
      end


      subroutine tmat_shear ( a, axis, s, b )

c*********************************************************************72
c
cc TMAT_SHEAR applies a shear to the geometric transformation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the current geometric transformation
c    matrix.
c
c    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
c    specifying the shear equation:
c
c      XY:  x' = x + s * y;
c      XZ:  x' = x + s * z;
c      YX:  y' = y + s * x;
c      YZ:  y' = y + s * z;
c      ZX:  z' = z + s * x;
c      ZY:  z' = z + s * y.
c
c    Input, double precision S, the shear coefficient.
c
c    Output, double precision B(4,4), the modified geometric transformation
c    matrix.  A and B may share the same memory.
c
      implicit none

      double precision a(4,4)
      character * ( 2 ) axis
      double precision b(4,4)
      double precision c(4,4)
      double precision s

      call tmat_init ( c )

      if ( axis .eq. 'XY' .or. axis .eq. 'xy' ) then
        c(1,2) = s
      else if ( axis .eq. 'XZ' .or. axis .eq. 'xz' ) then
        c(1,3) = s
      else if ( axis .eq. 'YX' .or. axis .eq. 'yx' ) then
        c(2,1) = s
      else if ( axis .eq. 'YZ' .or. axis .eq. 'yz' ) then
        c(2,3) = s
      else if ( axis .eq. 'ZX' .or. axis .eq. 'zx' ) then
        c(3,1) = s
      else if ( axis .eq. 'ZY' .or. axis .eq. 'zy' ) then
        c(3,2) = s
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TMAT_SHEAR - Fatal error!'
        write ( *, '(a)' ) '  Illegal shear axis: "' // axis // '".'
        write ( *, '(a)' )
     &    '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
        stop
      end if

      call r8mat_mm ( 4, 4, 4, c, a, b )

      return
      end
      subroutine tmat_trans ( a, t, b )

c*********************************************************************72
c
cc TMAT_TRANS applies a translation to the geometric transformation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries van Dam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1990.
c
c  Parameters:
c
c    Input, double precision A(4,4), the current geometric transformation
c    matrix.
c
c    Input, double precision T(3), the translation.  This may be thought
c    of as the point that the origin moves to under the translation.
c
c    Output, double precision B(4,4), the modified transformation matrix.
c    A and B may share the same memory.
c
      implicit none

      double precision a(4,4)
      double precision b(4,4)
      integer i
      integer j
      double precision t(3)

      do j = 1, 4
        do i = 1, 4
          b(i,j) = a(i,j)
        end do
      end do

      do i = 1, 3
        b(i,4) = b(i,4) + t(i)
      end do

      return
      end
      function torus_area_3d ( r1, r2 )

c*********************************************************************72
c
cc TORUS_AREA_3D returns the area of a torus in 3D.
c
c  Discussion:
c
c    A torus with radii R1 and R2 is the set of points P satisfying:
c
c    ( sqrt ( P(1)^2 + P(2)^2 ) - R1 )^2 + P(3)^2 <= R2^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
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
      subroutine torus_volume_3d ( r1, r2, volume )

c*********************************************************************72
c
cc TORUS_VOLUME_3D computes the volume of a torus in 3D.
c
c  Discussion:
c
c    A torus with radii R1 and R2 is the set of points P satisfying:
c
c    ( sqrt ( P(1)^2 + P(2)^2 ) - R1 )^2 + P(3)^2 <= R2^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R1, R2, the "inner" and "outer" radii of the
c    torus.
c
c    Output, double precision VOLUME, the volume of the torus.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision volume

      volume = 2.0D+00 * pi * pi * r1 * r2 * r2

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
      subroutine triangle_angles_2d ( t, angle )

c*********************************************************************72
c
cc TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
c
c  Discussion:
c
c    The law of cosines is used:
c
c      C^2 = A^2 + B^2 - 2 * A * B * COS ( GAMMA )
c
c    where GAMMA is the angle opposite side C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 May 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision ANGLE(3), the angles opposite
c    sides P1-P2, P2-P3 and P3-P1, in radians.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision angle(3)
      double precision b
      double precision c
      integer dim
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_acos
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
      a = 0.0D+00
      b = 0.0D+00
      c = 0.0D+00
      do dim = 1, dim_num
        a = a + ( t(dim,1) - t(dim,2) )**2
        b = b + ( t(dim,2) - t(dim,3) )**2
        c = c + ( t(dim,3) - t(dim,1) )**2
      end do
      a = sqrt ( a )
      b = sqrt ( b )
      c = sqrt ( c )
c
c  Take care of ridiculous special cases.
c
      if ( a .eq. 0.0D+00 .and.
     &     b .eq. 0.0D+00 .and.
     &     c .eq. 0.0D+00 ) then
        do i = 1, 3
          angle(i) = 2.0D+00 * pi / 3.0D+00
        end do
        return
      end if

      if ( c .eq. 0.0D+00 .or. a .eq. 0.0D+00 ) then
        angle(1) = pi
      else
        angle(1) = r8_acos ( ( c * c + a * a - b * b )
     &    / ( 2.0D+00 * c * a ) )
      end if

      if ( a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then
        angle(2) = pi
      else
        angle(2) = r8_acos ( ( a * a + b * b - c * c )
     &    / ( 2.0D+00 * a * b ) )
      end if

      if ( b .eq. 0.0D+00 .or. c .eq. 0.0D+00 ) then
        angle(3) = pi
      else
        angle(3) = r8_acos ( ( b * b + c * c - a * a )
     &    / ( 2.0D+00 * b * c ) )
      end if

      return
      end
      subroutine triangle_angles_3d ( t, angle )

c*********************************************************************72
c
cc TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
c
c  Discussion:
c
c    The law of cosines is used:
c
c      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
c
c    where GAMMA is the angle opposite side C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision ANGLE(3), the angles opposite
c    sides P1-P2, P2-P3 and P3-P1, in radians.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision a
      double precision angle(3)
      double precision b
      double precision c
      integer dim
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_acos
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
c
c  Compute the length of each side.
c
      a = 0.0D+00
      b = 0.0D+00
      c = 0.0D+00
      do dim = 1, dim_num
        a = a + ( t(dim,1) - t(dim,2) )**2
        b = b + ( t(dim,2) - t(dim,3) )**2
        c = c + ( t(dim,3) - t(dim,1) )**2
      end do
      a = sqrt ( a )
      b = sqrt ( b )
      c = sqrt ( c )
c
c  Take care of a ridiculous special case.
c
      if ( a .eq. 0.0D+00 .and.
     &     b .eq. 0.0D+00 .and.
     &     c .eq. 0.0D+00 ) then
        do i = 1, 3
          angle(i) = 2.0D+00 * pi / 3.0D+00
        end do
        return
      end if

      if ( c .eq. 0.0D+00 .or. a .eq. 0.0D+00 ) then
        angle(1) = pi
      else
        angle(1) = r8_acos ( ( c * c + a * a - b * b )
     &    / ( 2.0D+00 * c * a ) )
      end if

      if ( a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then
        angle(2) = pi
      else
        angle(2) = r8_acos ( ( a * a + b * b - c * c )
     &    / ( 2.0D+00 * a * b ) )
      end if

      if ( b .eq. 0.0D+00 .or. c .eq. 0.0D+00 ) then
        angle(3) = pi
      else
        angle(3) = r8_acos ( ( b * b + c * c - a * a )
     &    / ( 2.0D+00 * b * c ) )
      end if

      return
      end
      subroutine triangle_area_2d ( t, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_2D computes the area of a triangle in 2D.
c
c  Discussion:
c
c    If the triangle's vertices are given in counter clockwise order,
c    the area will be positive.  If the triangle's vertices are given
c    in clockwise order, the area will be negativec
c
c    An earlier version of this routine always returned the absolute
c    value of the computed area.  I am convinced now that that is
c    a less useful resultc  For instance, by returning the signed
c    area of a triangle, it is possible to easily compute the area
c    of a nonconvex polygon as the sum of the (possibly negative)
c    areas of triangles formed by node 1 and successive pairs of vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision AREA, the area of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision area
      double precision t(dim_num,3)

      area = 0.5D+00 * (
     &    t(1,1) * ( t(2,2) - t(2,3) )
     &  + t(1,2) * ( t(2,3) - t(2,1) )
     &  + t(1,3) * ( t(2,1) - t(2,2) ) )

      return
      end
      subroutine triangle_area_3d ( t, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_3D computes the area of a triangle in 3D.
c
c  Discussion:
c
c    This routine uses the fact that the norm of the cross product
c    of two vectors is the area of the parallelogram they form.
c
c    Therefore, the area of the triangle is half of that value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision AREA, the area of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision area
      double precision cross(dim_num)
      double precision t(dim_num,3)
c
c  Compute the cross product vector.
c
      cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) )
     &         - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

      cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) )
     &         - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

      cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) )
     &         - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      area = 0.5D+00 * sqrt ( cross(1)**2 + cross(2)**2 + cross(3)**2 )

      return
      end
      subroutine triangle_area_3d_2 ( t, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_3D_2 computes the area of a triangle in 3D.
c
c  Discussion:
c
c    This routine computes the area "the hard way".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision AREA, the area of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision alpha
      double precision area
      double precision base
      double precision dot
      double precision height
      double precision t(dim_num,3)
c
c  Find the projection of (P3-P1) onto (P2-P1).
c
      dot = ( t(1,2) - t(1,1) ) * ( t(1,3) - t(1,1) )
     &    + ( t(2,2) - t(2,1) ) * ( t(2,3) - t(2,1) )
     &    + ( t(3,2) - t(3,1) ) * ( t(3,3) - t(3,1) )
c
c  Find the length of (P2-P1).
c
      base = sqrt ( ( t(1,2) - t(1,1) )**2
     &            + ( t(2,2) - t(2,1) )**2
     &            + ( t(3,2) - t(3,1) )**2 )
c
c  The height of the triangle is the length of (P3-P1) after its
c  projection onto (P2-P1) has been subtracted.
c
      if ( base .eq. 0.0D+00 ) then

        height = 0.0D+00

      else

        alpha = dot / ( base * base )

        height = sqrt (
     &      ( t(1,1) + alpha * ( t(1,2) - t(1,1) ) - t(1,3) )**2
     &    + ( t(2,1) + alpha * ( t(2,2) - t(2,1) ) - t(2,3) )**2
     &    + ( t(3,1) + alpha * ( t(3,2) - t(3,1) ) - t(3,3) )**2 )

      end if

      area = 0.5D+00 * base * height

      return
      end
      subroutine triangle_area_3d_3 ( t, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_3D_3 computes the area of a triangle in 3D.
c
c  Discussion:
c
c    This routine uses Heron's formula
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision AREA, the area of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision area
      integer i
      integer j
      integer jp1
      double precision s(3)
      double precision t(dim_num,3)

      do j = 1, 3
        jp1 = mod ( j, 3 ) + 1
        s(j) = 0.0D+00
        do i = 1, dim_num
          s(j) = s(j) + ( t(i,j) - t(i,jp1) )**2
        end do
        s(j) = sqrt ( s(j) )
      end do

      area = (   s(1) + s(2) + s(3) )
     &     * ( - s(1) + s(2) + s(3) )
     &     * (   s(1) - s(2) + s(3) )
     &     * (   s(1) + s(2) - s(3) )

      if ( area .lt. 0.0D+00 ) then
        area = -1.0D+00
        return
      end if

      area = 0.25D+00 * sqrt ( area )

      return
      end
      subroutine triangle_area_heron ( s, area )

c*********************************************************************72
c
cc TRIANGLE_AREA_HERON computes the area of a triangle using Heron's formula.
c
c  Discussion:
c
c    The formula is valid for any spatial dimension, depending only
c    on the lengths of the sides, and not the coordinates of the vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision S(3), the lengths of the three sides.
c
c    Output, double precision AREA, the area of the triangle, or -1.0 if the
c    sides cannot constitute a triangle.
c
      implicit none

      double precision area
      double precision s(3)

      area = (   s(1) + s(2) + s(3) )
     &     * ( - s(1) + s(2) + s(3) )
     &     * (   s(1) - s(2) + s(3) )
     &     * (   s(1) + s(2) - s(3) )

      if ( area .lt. 0.0D+00 ) then
        area = -1.0D+00
        return
      end if

      area = 0.25D+00 * sqrt ( area )

      return
      end
      subroutine triangle_area_vector_3d ( t, area_vector )

c*********************************************************************72
c
cc TRIANGLE_AREA_VECTOR_3D computes the area vector of a triangle in 3D.
c
c  Discussion:
c
c    The "area vector" of a triangle is simply a cross product of,
c    for instance, the vectors (V2-V1) and (V3-V1), where V1, V2
c    and V3 are the vertices of the triangle.
c
c    The norm of the cross product vector of two vectors is the area
c    of the parallelogram they form.
c
c    Therefore, the area of the triangle is half of the norm of the
c    area vector:
c
c      area = 0.5 * sqrt ( sum ( area_vector(1:3)^2 ) )
c
c    The reason for looking at the area vector rather than the area
c    is that this makes it possible to compute the area of a flat
c    polygon in 3D by summing the areas of the triangles that form
c    a decomposition of the polygon, while allowing for both positive
c    and negative areas.  (Sum the vectors, THEN take the norm and
c    multiply by 1/2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision AREA_VECTOR(3), the area vector of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision area_vector(dim_num)
      double precision t(dim_num,3)

      area_vector(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) )
     &               - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

      area_vector(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) )
     &               - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

      area_vector(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) )
     &               - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      return
      end
      subroutine triangle_barycentric_2d ( t, p, xsi )

c*********************************************************************72
c
cc TRIANGLE_BARYCENTRIC_2D finds the barycentric coordinates of a point in 2D.
c
c  Discussion:
c
c    The barycentric coordinate of point P related to vertex A can be
c    interpreted as the ratio of the area of the triangle with
c    vertex A replaced by vertex P to the area of the original
c    triangle.
c
c    This routine assumes that the triangle vertices are given in
c    counter clockwise order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c    The vertices should be given in counter clockwise order.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision XSI(3), the barycentric coordinates of P
c    with respect to the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer rhs_num
      parameter ( rhs_num = 1 )

      double precision a(dim_num,dim_num+rhs_num)
      integer info
      double precision p(dim_num)
      double precision t(dim_num,3)
      double precision xsi(dim_num+1)
c
c  Set up the linear system
c
c    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
c    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
c
c  which is satisfied by the barycentric coordinates of P.
c
      a(1,1) = t(1,2) - t(1,1)
      a(1,2) = t(1,3) - t(1,1)
      a(1,3) = p(1)   - t(1,1)

      a(2,1) = t(2,2) - t(2,1)
      a(2,2) = t(2,3) - t(2,1)
      a(2,3) = p(2)   - t(2,1)
c
c  Solve the linear system.
c
      call r8mat_solve ( dim_num, rhs_num, a, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_BARYCENTRIC_2D - Fatal error!'
        write ( *, '(a)' ) '  The linear system is singular.'
        write ( *, '(a)' )
     &    '  The input data does not form a proper triangle.'
        stop
      end if

      xsi(1) = a(1,3)
      xsi(2) = a(2,3)
      xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

      return
      end
      subroutine triangle_centroid_2d ( t, centroid )

c*********************************************************************72
c
cc TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
c
c  Discussion:
c
c    The centroid of a triangle can also be considered the
c    center of gravity, or center of mass, assuming that the triangle
c    is made of a thin uniform sheet of massy material.
c
c    The centroid of a triangle is the intersection of the medians.
c
c    A median of a triangle is a line connecting a vertex to the
c    midpoint of the opposite side.
c
c    In barycentric coordinates, in which the vertices of the triangle
c    have the coordinates (1,0,0), (0,1,0) and (0,0,1), the centroid
c    has coordinates (1/3,1/3,1/3).
c
c    In geometry, the centroid of a triangle is often symbolized by "G".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision CENTROID(2), the coordinates of the centroid.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision centroid(dim_num)
      integer i
      integer j
      double precision t(dim_num,3)

      do i = 1, dim_num
        centroid(i) = 0.0D+00
        do j = 1, 3
          centroid(i) = centroid(i) + t(i,j)
        end do
        centroid(i) = centroid(i) / 3.0D+00
      end do

      return
      end
      subroutine triangle_centroid_3d ( t, centroid )

c*********************************************************************72
c
cc TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
c
c  Discussion:
c
c    The centroid of a triangle can also be considered the
c    center of gravity or center of mass, assuming that the triangle
c    is made of a thin uniform sheet of massy material.
c
c    The centroid of a triangle is the intersection of the medians.
c    A median of a triangle is a line connecting any vertex to the
c    midpoint of the opposite side.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Output, double precision CENTROID(3), the coordinates of the centroid.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision centroid(dim_num)
      integer i
      integer j
      double precision t(dim_num,3)

      do i = 1, dim_num
        centroid(i) = 0.0D+00
        do j = 1, 3
          centroid(i) = centroid(i) + t(i,j)
        end do
        centroid(i) = centroid(i) / 3.0D+00
      end do

      return
      end
      subroutine triangle_circumcenter_2d ( t, pc )

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
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision PC(2), the circumcenter of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision det
      double precision f(2)
      double precision pc(dim_num)
      double precision t(dim_num,3)
      double precision top(dim_num)

      f(1) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
      f(2) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

      top(1) =    ( t(2,3) - t(2,1) ) * f(1)
     &          - ( t(2,2) - t(2,1) ) * f(2)
      top(2) =  - ( t(1,3) - t(1,1) ) * f(1)
     &          + ( t(1,2) - t(1,1) ) * f(2)

      det  =    ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )
     &        - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      pc(1) = t(1,1) + 0.5D+00 * top(1) / det
      pc(2) = t(2,1) + 0.5D+00 * top(2) / det

      return
      end
      subroutine triangle_circumcenter_2d_2 ( t, pc )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMCENTER_2D_2 computes the circumcenter of a triangle in 2D.
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
c    Surprisingly, the diameter of the circle can be found by solving
c    a 2 by 2 linear system.  If we label the vertices of the triangle
c    P1, P2 and P3, then the vectors P2 - P1 and P3 - P1 are secants of
c    the circle, and each forms a right triangle with the diameter
c    vector through P1.
c
c    Hence, the dot product of P2 - P1 with the diameter vector is equal
c    to the square of the length of P2 - P1, and similarly for P3 - P1.
c    This determines the diameter vector originating at P1.
c
c    In geometry, the circumcenter of a triangle is often symbolized by "O".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision PC(2), the circumcenter of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer rhs_num
      parameter ( rhs_num = 1 )

      double precision a(dim_num,dim_num+rhs_num)
      integer i
      integer info
      double precision pc(dim_num)
      double precision t(dim_num,3)
c
c  Set up the linear system.
c
      a(1,1) = t(1,2) - t(1,1)
      a(1,2) = t(2,2) - t(2,1)
      a(1,3) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2

      a(2,1) = t(1,3) - t(1,1)
      a(2,2) = t(2,3) - t(2,1)
      a(2,3) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
c
c  Solve the linear system.
c
      call r8mat_solve ( dim_num, rhs_num, a, info )
c
c  Compute the center
c
      if ( info .ne. 0 ) then
        do i = 1, dim_num
          pc(i) = 0.0D+00
        end do
      else
        do i = 1, dim_num
          pc(i) = t(i,1) + 0.5D+00 * a(i,dim_num+1)
        end do
      end if

      return
      end
      subroutine triangle_circumcenter ( n, t, p )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMCENTER computes the circumcenter of a triangle in ND.
c
c  Discussion:
c
c    Three ND points A, B and C lie on a circle.
c
c    The circumcenter P has the formula
c
c      P = ( Area ( PBC ) * A + Area ( APC) * B + Area ( ABP ) * C )
c        / ( Area ( PBC )     + Area ( APC )    + Area ( ABP ) )
c
c    The details of the formula rely on information supplied
c    by Oscar Lanzi III.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision T(N,3), the triangle vertices.
c
c    Output, double precision P(N), the circumcenter of the triangle.
c
      implicit none

      integer n

      double precision a
      double precision abp
      double precision apc
      double precision b
      double precision c
      integer i
      double precision p(n)
      double precision pbc
      double precision r8vec_normsq_affine
      double precision t(n,3)

      a = r8vec_normsq_affine ( n, t(1,2), t(1,3) )
      b = r8vec_normsq_affine ( n, t(1,3), t(1,1) )
      c = r8vec_normsq_affine ( n, t(1,1), t(1,2) )

      pbc = a * ( - a + b + c )
      apc = b * (   a - b + c )
      abp = c * (   a + b - c )

      do i = 1, n
        p(i) = ( pbc * t(i,1) + apc * t(i,2) + abp * t(i,3) )
     &       / ( pbc          + apc          + abp )
      end do

      return
      end
      subroutine triangle_circumcircle_2d ( t, r, pc )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
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
c    07 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision R, PC(2), the circumradius and circumcenter
c    of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision bot
      double precision c
      double precision det
      double precision f(2)
      double precision pc(dim_num)
      double precision r
      double precision top(dim_num)
      double precision t(dim_num,3)
c
c  Circumradius.
c
      a = sqrt ( ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2 )
      b = sqrt ( ( t(1,3) - t(1,2) )**2 + ( t(2,3) - t(2,2) )**2 )
      c = sqrt ( ( t(1,1) - t(1,3) )**2 + ( t(2,1) - t(2,3) )**2 )

      bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c )
     &  * (   a + b - c )

      if ( bot .le. 0.0D+00 ) then
        r = -1.0D+00
        pc(1) = 0.0D+00
        pc(2) = 0.0D+00
        return
      end if

      r = a * b * c / sqrt ( bot )
c
c  Circumcenter.
c
      f(1) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
      f(2) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

      top(1) =    ( t(2,3) - t(2,1) ) * f(1)
     &          - ( t(2,2) - t(2,1) ) * f(2)
      top(2) =  - ( t(1,3) - t(1,1) ) * f(1)
     &          + ( t(1,2) - t(1,1) ) * f(2)

      det  =    ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )
     &        - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      pc(1) = t(1,1) + 0.5D+00 * top(1) / det
      pc(2) = t(2,1) + 0.5D+00 * top(2) / det


      return
      end
      subroutine triangle_circumcircle_2d_2 ( t, r, pc )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcircle of a triangle in 2D.
c
c  Discussion:
c
c    The circumscribed circle of a triangle is the circle that passes through
c    the three vertices of the triangle.  The circumscribed circle contains
c    the triangle, but it is not necessarily the smallest triangle to do so.
c
c    Surprisingly, the diameter of the circle can be found by solving
c    a 2 by 2 linear system.  This is because the vectors P2 - P1
c    and P3 - P1 are secants of the circle, and each forms a right
c    triangle with the diameter.  Hence, the dot product of
c    P2 - P1 with the diameter is equal to the square of the length
c    of P2 - P1, and similarly for P3 - P1.  This determines the
c    diameter vector originating at P1.
c
c    If all angles of the triangle are no greater than 90 degrees, then
c    the center of the circumscribed circle will lie inside the triangle.
c    Otherwise, the center will lie outside the triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision R, PC(2), the circumradius and circumcenter.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer rhs_num
      parameter ( rhs_num = 1 )

      double precision a(dim_num,dim_num+rhs_num)
      integer info
      double precision pc(dim_num)
      double precision r
      double precision t(dim_num,3)
c
c  Set up the linear system.
c
      a(1,1) = t(1,2) - t(1,1)
      a(1,2) = t(2,2) - t(2,1)
      a(1,3) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2

      a(2,1) = t(1,3) - t(1,1)
      a(2,2) = t(2,3) - t(2,1)
      a(2,3) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
c
c  Solve the linear system.
c
      call r8mat_solve ( dim_num, rhs_num, a, info )

      if ( info .ne. 0 ) then
        r = -1.0D+00
        pc(1:dim_num) = 0.0D+00
      end if

      r = 0.5D+00 * sqrt ( a(1,dim_num+1) * a(1,dim_num+1)
     &                   + a(2,dim_num+1) * a(2,dim_num+1) )
      pc(1) = t(1,1) + 0.5D+00 * a(1,3)
      pc(2) = t(2,1) + 0.5D+00 * a(2,3)

      return
      end
      subroutine triangle_circumradius_2d ( t, r )

c*********************************************************************72
c
cc TRIANGLE_CIRCUMRADIUS_2D computes the circumradius of a triangle in 2D.
c
c  Discussion:
c
c    The circumscribed circle of a triangle is the circle that passes through
c    the three vertices of the triangle.  The circumscribed circle contains
c    the triangle, but it is not necessarily the smallest triangle to do so.
c
c    The circumradius of a triangle is the radius of the circumscribed
c    circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision R, the circumradius of the circumscribed circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision bot
      double precision c
      double precision r
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
      a = sqrt ( ( t(1,1) - t(1,2) )**2 + ( t(2,1) - t(2,2) )**2 )
      b = sqrt ( ( t(2,2) - t(2,3) )**2 + ( t(2,2) - t(2,3) )**2 )
      c = sqrt ( ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2 )

      bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c )
     &  * (   a + b - c )

      if ( bot .le. 0.0D+00 ) then
        r = -1.0D+00
        return
      end if

      r = a * b * c / sqrt ( bot )

      return
      end
      subroutine triangle_diameter_2d ( t, diameter )

c*********************************************************************72
c
cc TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
c
c  Discussion:
c
c    The diameter of a triangle is the diameter of the smallest circle
c    that can be drawn around the triangle.  At least two of the vertices
c    of the triangle will intersect the circle, but not necessarily
c    all three!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision DIAMETER, the diameter of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision asq
      double precision b
      double precision bsq
      double precision c
      double precision csq
      double precision diameter
      double precision r8vec_diff_norm_squared
      double precision t(dim_num,3)
c
c  Compute the squared length of each side.
c
      asq = r8vec_diff_norm_squared ( dim_num, t(1,1), t(1,2) )
      bsq = r8vec_diff_norm_squared ( dim_num, t(1,2), t(1,3) )
      csq = r8vec_diff_norm_squared ( dim_num, t(1,3), t(1,1) )
c
c  Take care of a zero side.
c
      if ( asq .eq. 0.0D+00 ) then
        diameter = sqrt ( bsq )
        return
      else if ( bsq .eq. 0.0D+00 ) then
        diameter = sqrt ( csq )
        return
      else if ( csq .eq. 0.0D+00 ) then
        diameter = sqrt ( asq )
        return
      end if
c
c  Make ASQ the largest.
c
      if ( asq .lt. bsq ) then
        call r8_swap ( asq, bsq )
      end if

      if ( asq .lt. csq ) then
        call r8_swap ( asq, csq )
      end if
c
c  If ASQ is very large...
c
      if ( bsq + csq .lt. asq ) then

        diameter = sqrt ( asq )

      else

        a = sqrt ( asq )
        b = sqrt ( bsq )
        c = sqrt ( csq )

        diameter = 2.0D+00 * a * b * c / sqrt ( ( a + b + c ) 
     &    * ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) )

      end if

      return
      end
      subroutine triangle_edge_length_2d ( t, edge_length )

c*********************************************************************72
c
cc TRIANGLE_EDGE_LENGTH_2D returns edge lengths of a triangle in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision EDGE_LENGTH(3), the length of the edges.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision edge_length(3)
      integer i4_wrap
      integer j1
      integer j2
      double precision r8vec_diff_norm
      double precision t(dim_num,3)

      do j1 = 1, 3
        j2 = i4_wrap ( j1 + 1, 1, 3 )
        edge_length(j1) = r8vec_diff_norm ( dim_num, t(1,j2), t(1,j1) )
      end do

      return
      end
      subroutine triangle_gridpoints_2d ( t, sub_num, grid_max, 
     &  grid_num, g )

c*********************************************************************72
c
cc TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
c
c  Discussion:
c
c    The gridpoints are computed by repeated halving of the triangle.
c    The 0-th set of grid points is the vertices themselves.
c    The first set of grid points is the midpoints of the sides.
c    These points can be used to draw 4 triangles that make up the original
c    triangle.  The second set of grid points is the side midpoints and centers
c    of these four triangles.
c
c     SUB_NUM                     GRID_NUM
c    -----                        -----
c        0      1                  =  1  (centroid)
c        1      1 + 2              =  3  (vertices)
c        2      1 + 2 + 3          =  6
c        3      1 + 2 + 3 + 4      = 10
c        4      1 + 2 + 3 + 4 + 5  = 15
c
c    GRID_NUM is the sum of the integers from 1 to SUB_NUM+1 or
c
c      GRID_NUM = (SUB_NUM+1) * (SUB_NUM+2) / 2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, integer SUB_NUM, the number of subdivisions.
c
c    Input, integer GRID_MAX, the maximum number of grid points.
c
c    Output, integer GRID_NUM, the number of grid points returned.
c
c    Output, double precision G(2,GRID_MAX), the grid points.
c
      implicit none

      integer grid_max
      integerdim_num
      parameter ( dim_num = 2 )

      double precision g(dim_num,grid_max)
      integer i
      integer j
      integer grid_num
      integer sub_num
      double precision t(dim_num,3)

      grid_num = 0
c
c  Special case, SUB_NUM = 0.
c
      if ( sub_num .eq. 0 ) then
        if ( 1 .le. grid_max ) then
          grid_num = 1
          g(1,1) = ( t(1,1) + t(1,2) + t(1,3) ) / 3.0D+00
          g(2,1) = ( t(2,1) + t(2,2) + t(2,3) ) / 3.0D+00
        end if
        return
      end if

      do i = 0, sub_num

        do j = 0, sub_num - i

          if ( grid_num .lt. grid_max ) then

            grid_num = grid_num + 1

            g(1,grid_num) = ( dble (           i     ) * t(1,1) 
     &                      + dble (               j ) * t(1,2) 
     &                      + dble ( sub_num - i - j ) * t(1,3) ) 
     &                      / dble ( sub_num         )

            g(2,grid_num) = ( dble (           i     ) * t(2,1) 
     &                      + dble (               j ) * t(2,2) 
     &                      + dble ( sub_num - i - j ) * t(2,3) ) 
     &                      / dble ( sub_num )
          end if

        end do
      end do

      return
      end
      subroutine triangle_incenter_2d ( t, pc )

c*********************************************************************72
c
cc TRIANGLE_INCENTER_2D computes the incenter of a triangle in 2D.
c
c  Discussion:
c
c    The incenter of a triangle is the center of the inscribed circle.
c
c    The inscribed circle of a triangle is the largest circle that can
c    be drawn inside the triangle.
c
c    The inscribed circle is tangent to all three sides of the triangle.
c
c    The angle bisectors of the triangle intersect at the center of the
c    inscribed circle.
c
c    In geometry, the incenter is often represented by "I".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision PC(2), the incenter.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      integer i
      double precision pc(dim_num)
      double precision perimeter
      double precision r8vec_diff_norm
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
      a = r8vec_diff_norm ( t(1,1), t(1,2) )
      b = r8vec_diff_norm ( t(1,2), t(1,3) )
      c = r8vec_diff_norm ( t(1,3), t(1,1) )

      perimeter = a + b + c

      if ( perimeter .eq. 0.0D+00 ) then

        pc(1) = t(1,1)
        pc(2) = t(2,1)

      else

        do i = 1, dim_num
          pc(i) = (  
     &        b * t(i,1) 
     &      + c * t(i,2) 
     &      + a * t(i,3) ) / perimeter
        end do

      end if

      return
      end
      subroutine triangle_incircle_2d ( t, r, pc )

c*********************************************************************72
c
cc TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
c
c  Discussion:
c
c    The inscribed circle of a triangle is the largest circle that can
c    be drawn inside the triangle.  It is tangent to all three sides,
c    and the lines from its center to the vertices bisect the angles
c    made by each vertex.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision R, PC(2), the radius and center of the
c    inscribed circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      integer i
      double precision pc(dim_num)
      double precision perimeter
      double precision r
      double precision r8vec_diff_norm
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
      a = r8vec_diff_norm ( t(1,1), t(1,2) )
      b = r8vec_diff_norm ( t(1,2), t(1,3) )
      c = r8vec_diff_norm ( t(1,3), t(1,1) )

      perimeter = a + b + c

      if ( perimeter .eq. 0.0D+00 ) then
        pc(1) = t(1,1)
        pc(2) = t(2,1)
        r = 0.0D+00
        return
      end if

      do i = 1, dim_num
        pc(i) = (  
     &      b * t(i,1) 
     &    + c * t(i,2) 
     &    + a * t(i,3) ) / perimeter
      end do

      r = 0.5D+00 * sqrt ( 
     &    ( - a + b + c )  
     &  * ( + a - b + c ) 
     &  * ( + a + b - c ) / perimeter )

      return
      end
      subroutine triangle_inradius_2d ( t, r )

c*********************************************************************72
c
cc TRIANGLE_INRADIUS_2D: radius of the inscribed circle of a triangle in 2D.
c
c  Discussion:
c
c    The inscribed circle of a triangle is the largest circle that can
c    be drawn inside the triangle.  It is tangent to all three sides,
c    and the lines from its center to the vertices bisect the angles
c    made by each vertex.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision R, the radius of the inscribed circle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision perimeter
      double precision r
      double precision r8vec_diff_norm
      double precision t(dim_num,3)
!
!  Compute the length of each side.
!
      a = r8vec_diff_norm ( t(1,1), t(1,2) )
      b = r8vec_diff_norm ( t(1,2), t(1,3) )
      c = r8vec_diff_norm ( t(1,3), t(1,1) )

      perimeter = a + b + c

      if ( perimeter .eq. 0.0D+00 ) then
        r = 0.0D+00
        return
      end if

      r = 0.5D+00 * sqrt ( 
     &    ( - a + b + c )  
     &  * ( + a - b + c ) 
     &  * ( + a + b - c ) / perimeter )

      return
      end
      function triangle_is_degenerate_nd ( dim_num, t )

c*********************************************************************72
c
cc TRIANGLE_IS_DEGENERATE_ND finds if a triangle is degenerate in ND.
c
c  Discussion:
c
c    A triangle in ND is described by the coordinates of its 3 vertices.
c
c    A triangle in ND is degenerate if any two vertices are equal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision T(DIM_NUM,3), the triangle vertices.
c
c    Output, logical TRIANGLE_IS_DEGENERATE_ND, is TRUE if the
c    triangle is degenerate.
c
      implicit none

      integer dim_num

      logical r8vec_eq
      double precision t(dim_num,3)
      logical triangle_is_degenerate_nd

      triangle_is_degenerate_nd = (
     &  r8vec_eq ( dim_num, t(1,1), t(1,2) ) .or.
     &  r8vec_eq ( dim_num, t(1,2), t(1,3) ) .or.
     &  r8vec_eq ( dim_num, t(1,3), t(1,1) ) )

      return
      end
      subroutine triangle_lattice_layer_point_next ( c, v, more )

c*********************************************************************72
c
cc TRIANGLE_LATTICE_LAYER_POINT_NEXT: next triangle lattice layer point.
c
c  Discussion:
c
c    The triangle lattice layer L is bounded by the lines
c
c      0 .le. X,
c      0 .le. Y,
c      L - 1 < X / C(1) + Y / C(2) .le. L.
c
c    In particular, layer L = 0 always contains the single point (0,0).
c
c    This function returns, one at a time, the points that lie within
c    a given triangle lattice layer.
c
c    Thus, if we set C(1) = 2, C(2) = 3, then we get the following layers:
c
c    L = 0: (0,0)
c    L = 1: (1,0), (2,0), (0,1), (1,1), (0,2), (0,3)
c    L = 2: (3,0), (4,0), (2,1), (3,1), (1,2), (2,2), (1,3), (2,3),
c           (0,4), (1,4), (0,5), (0,6).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer C(3), coefficients defining the
c    lattice layer.  Entry C(3) contains the layer index.
c    C(1) and C(2) should be positive, and C(3) must be nonnegative.
c
c    Input/output, integer V(2).  On first call for a given layer,
c    the input value of V is not important.  On a repeated call for the same
c    layer, the input value of V should be the output value from the previous
c    call.  On output, V contains the next lattice layer point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given layer.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if the returned value V is a new point.
c    If the output value is FALSE, then no more points were found,
c    and V was reset to 0, and the lattice layer has been exhausted.
c
      implicit none

      integer c(3)
      integer c1n
      integer i4vec_lcm
      logical more
      integer n
      parameter ( n = 2 )
      integer rhs1
      integer rhs2
      integer v(2)
c
c  Treat layer C(3) = 0 specially.
c
      if ( c(3) .eq. 0 ) then
        if ( .not. more ) then
          v(1) = 0
          v(2) = 0
          more = .true.
        else
          more = .false.
        end if
        return
      end if
c
c  Compute first point.
c
      if ( .not. more ) then

        v(1) = ( c(3) - 1 ) * c(1) + 1
        v(2) = 0
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )

        rhs1 = c1n * ( c(3) - 1 )
        rhs2 = c1n *   c(3)

        if ( c(2) * ( v(1) + 1 ) + c(1) * v(2) .le. rhs2 ) then
          v(1) = v(1) + 1
        else
          v(1) = ( rhs1 - c(1) * ( v(2) + 1 ) ) / c(2)
          v(1) = max ( v(1), 0 )
          v(2) = v(2) + 1
          if ( c(2) * v(1) + c(1) * v(2) .le. rhs1 ) then
            v(1) = v(1) + 1
          end if
          if ( c(2) * v(1) + c(1) * v(2) .le. rhs2 ) then

          else
            v(1) = 0
            v(2) = 0
            more = .false.
          end if
        end if
      end if

      return
      end
      subroutine triangle_lattice_point_next ( c, v, more )

c*********************************************************************72
c
cc TRIANGLE_LATTICE_POINT_NEXT returns the next triangle lattice point.
c
c  Discussion:
c
c    The lattice triangle is defined by the vertices:
c
c      (0,0), (C(3)/C(1), 0) and (0,C(3)/C(2))
c
c    The lattice triangle is bounded by the lines
c
c      0 <= X,
c      0 <= Y
c      X / C(1) + Y / C(2) <= C(3)
c
c    Lattice points are listed one at a time, starting at the origin,
c    with X increasing first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer C(3), coefficients defining the
c    lattice triangle.  These should be positive.
c
c    Input/output, integer V(2).  On first call, the input
c    value is not important.  On a repeated call, the input value should
c    be the output value from the previous call.  On output, V contains
c    the next lattice point.
c
c    Input/output, logical MORE.  On input, set MORE to FALSE to indicate
c    that this is the first call for a given triangle.  Thereafter, the input
c    value should be the output value from the previous call.  On output,
c    MORE is TRUE if the returned value V is a new lattice point.
c    If the output value is FALSE, then no more lattice points were found,
c    and V was reset to 0, and the routine should not be called further
c    for this triangle.
c
      implicit none

      integer c(3)
      integer c1n
      integer i4vec_lcm
      logical more
      integer n
      parameter ( n = 2 )
      integer rhs
      integer v(2)

      if ( .not. more ) then

        v(1) = 0
        v(2) = 0
        more = .true.

      else

        c1n = i4vec_lcm ( n, c )
        rhs = c1n * c(3)

        if ( c(2) * ( v(1) + 1 ) + c(1) * v(2) .le. rhs ) then
          v(1) = v(1) + 1
        else
          v(1) = 0
          if ( c(2) * v(1) + c(1) * ( v(2) + 1 ) .le. rhs ) then
            v(2) = v(2) + 1
          else
            v(2) = 0
            more = .false.
          end if
        end if

      end if

      return
      end
      subroutine triangle_line_imp_int_2d ( t, a, b, c, int_num, pint )

c*********************************************************************72
c
cc TRIANGLE_LINE_IMP_INT_2D: implicit line intersects a triangle in 2D.
c
c  Discussion:
c
c    An implicit line is the set of points ( X, Y ) satisfying
c
c      A * X + B * Y + C = 0
c
c    where at least one of A and B is not zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c   07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, double precision A, B, C, determine the equation of the line:
c    A*X + B*Y + C = 0.
c
c    Output, integer INT_NUM, the number of points of intersection
c    of the line with the triangle.  INT_NUM may be 0, 1, 2 or 3.
c
c    Output, double precision PINT(2,3), contains the intersection points.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision a1
      double precision b
      double precision b1
      double precision c
      double precision c1
      integer i
      integer i4_wrap
      integer int_num
      integer ival
      integer j
      double precision p(dim_num)
      double precision pint(dim_num,3)
      double precision t(dim_num,3)
      double precision test1
      double precision test2

      int_num = 0

      do i = 1, 3

        j = i4_wrap ( i+1, 1, 3 )
c
c  Get the implicit form of the line through vertices I and I+1.
c
        call line_exp2imp_2d ( t(1:2,i), t(1:2,j), a1, b1, c1 )
c
c  Seek an intersection with the original line.
c
        call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, p )
c
c  If there is an intersection, determine if it lies between the two vertices.
c
        if ( ival .eq. 1 ) then

          test1 = ( p(1)   - t(1,i) ) * ( t(1,j) - t(1,i) ) 
     &          + ( p(2)   - t(2,i) ) * ( t(2,j) - t(2,i) )

          test2 = ( t(1,j) - t(1,i) ) * ( t(1,j) - t(1,i) )
     &          + ( t(2,j) - t(2,i) ) * ( t(2,j) - t(2,i) )

          if ( 0 .le. test1 .and. test1 .le. test2 ) then
            int_num = int_num + 1
            pint(1,int_num) = p(1)
            pint(2,int_num) = p(2)
          end if

        end if

      end do

      return
      end
      function triangle_orientation_2d ( t )

c*********************************************************************72
c
cc TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
c
c  Discussion:
c
c    Three distinct non-colinear points in the plane define a circle.
c    If the points are visited in the order P1, P2, and then
c    P3, this motion defines a clockwise or counter clockwise
c    rotation along the circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, integer TRIANGLE_ORIENTATION_2D, reports if the 
c    three points lie clockwise on the circle that passes through them.  
c    The possible return values are:
c    0, the points are distinct, noncolinear, and lie counter clockwise
c    on their circle.
c    1, the points are distinct, noncolinear, and lie clockwise
c    on their circle.
c    2, the points are distinct and colinear.
c    3, at least two of the points are identical.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision det
      integer triangle_orientation_2d
      double precision t(dim_num,3)

      if ( ( t(1,1) .eq. t(1,2) .and. t(2,1) .eq. t(2,2) ) .or. 
     &     ( t(1,2) .eq. t(1,3) .and. t(2,2) .eq. t(2,3) ) .or. 
     &     ( t(1,3) .eq. t(1,1) .and. t(2,3) .eq. t(2,1) ) ) then
        triangle_orientation_2d = 3
        return
      end if

      det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) 
     &    - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

      if ( det .eq. 0.0D+00 ) then
        triangle_orientation_2d = 2
      else if ( det .lt. 0.0D+00 ) then
        triangle_orientation_2d = 1
      else if ( 0.0D+00 .lt. det ) then
        triangle_orientation_2d = 0
      end if

      return
      end
      subroutine triangle_orthocenter_2d ( t, pc )

c*********************************************************************72
c
cc TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
c
c  Discussion:
c
c    The orthocenter is defined as the intersection of the three altitudes
c    of a triangle.
c
c    An altitude of a triangle is the line through a vertex of the triangle
c    and perpendicular to the opposite side.
c
c    In geometry, the orthocenter of a triangle is often symbolized by "H".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision PC(2), the orthocenter of the triangle.
c
c    Output, logical FLAG, is TRUE if the value could not be computed.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      logical flag
      integer ival
      double precision p23(dim_num)
      double precision p31(dim_num)
      double precision pc(dim_num)
      double precision r8_huge
      double precision t(dim_num,3)
c
c  Determine a point P23 common to the line (P2,P3) and
c  its perpendicular through P1.
c
      call line_exp_perp_2d ( t(1,2), t(1,3), t(1,1), p23, flag )

      if ( flag ) then
        pc(1) = r8_huge ( )
        pc(2) = r8_huge ( )
        return
      end if
c
c  Determine a point P31 common to the line (P3,P1) and
c  its perpendicular through P2.
c
      call line_exp_perp_2d ( t(1,3), t(1,1), t(1,2), p31, flag )

      if ( flag ) then
        pc(1) = r8_huge ( )
        pc(2) = r8_huge ( )
        return
      end if
c
c  Determine PC, the intersection of the lines (P1,P23) and (P2,P31).
c
      call lines_exp_int_2d ( t(1,1), p23, t(1,2), p31, ival, pc )

      if ( ival .ne. 1 ) then
        pc(1) = r8_huge ( )
        pc(2) = r8_huge ( )
        flag = .true.
        return
      end if

      return
      end
      subroutine triangle_point_dist_2d ( t, p, dist )

c*********************************************************************72
c
cc TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, double precision P(2), the point to be checked.
c
c    Output, double precision DIST, the distance from the point to the
c    triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer side_num
      parameter ( side_num = 3 )

      double precision dist
      double precision dist2
      integer i4_wrap
      integer j
      integer jp1
      double precision p(dim_num)
      double precision r8_huge
      double precision t(dim_num,side_num)
c
c  Find the distance to each of the line segments.
c
      dist = r8_huge ( )

      do j = 1, side_num

        jp1 = i4_wrap ( j+1, 1, side_num )

        call segment_point_dist_2d ( t(1,j), t(1,jp1), p, dist2 )

        if ( dist2 .lt. dist ) then
          dist = dist2
        end if

      end do

      return
      end
      subroutine triangle_point_dist_3d ( t, p, dist )

c*********************************************************************72
c
cc TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(3,3), the triangle vertices.
c
c    Input, double precision P(3), the point which is to be checked.
c
c    Output, double precision DIST, the distance from the point to the
c    triangle.  DIST is zero if the point lies exactly on the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dist
      double precision dist2
      double precision p(dim_num)
      double precision t(dim_num,3)
c
c  Compute the distances from the point to each of the sides.
c
      call segment_point_dist_3d ( t(1,1), t(1,2), p, dist2 )

      dist = dist2

      call segment_point_dist_3d ( t(1,2), t(1,3), p, dist2 )

      dist = min ( dist, dist2 )

      call segment_point_dist_3d ( t(1,3), t(1,1), p, dist2 )

      dist = min ( dist, dist2 )

      return
      end
      subroutine triangle_point_dist_signed_2d ( t, p, dist_signed )

c*********************************************************************72
c
cc TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
c
c  Discussion:
c
c    If the signed distance is:
c    0, the point is on the boundary of the triangle;
c    negative, the point is in the triangle;
c    positive, the point is outside the triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c    These should be given in counter clockwise order.
c
c    Input, double precision P(2), the point which is to be checked.
c
c    Output, double precision DIST_SIGNED, the signed distance from the
c    point to the triangle.  
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision dis12
      double precision dis23
      double precision dis31
      double precision dist_signed
      double precision p(dim_num)
      double precision t(dim_num,3)
c
c  Compute the signed line distances to the point.
c
      call line_exp_point_dist_signed_2d ( t(1,1), t(1,2), p, dis12 )

      call line_exp_point_dist_signed_2d ( t(1,2), t(1,3), p, dis23 )

      call line_exp_point_dist_signed_2d ( t(1,3), t(1,1), p, dis31 )
c
c  If the point is inside the triangle, all the line distances are negative.
c  The largest (negative) line distance has the smallest magnitude,
c  and is the signed triangle distance.
c
      if ( dis12 .le. 0.0D+00 .and. 
     &     dis23 .le. 0.0D+00 .and. 
     &     dis31 .le. 0.0D+00 ) then
        dist_signed = max ( dis12, dis23, dis31 )
c
c  If the point is outside the triangle, then we have to compute
c  the (positive) line segment distances and take the minimum.
c
      else

        call segment_point_dist_2d ( t(1,1), t(1,2), p, dis12 )
        call segment_point_dist_2d ( t(1,2), t(1,3), p, dis23 )
        call segment_point_dist_2d ( t(1,3), t(1,1), p, dis31 )

        dist_signed = min ( dis12, dis23, dis31 )

      end if

      return
      end
      subroutine triangle_point_near_2d ( t, p, pn, dist )

c*********************************************************************72
c
cc TRIANGLE_POINT_NEAR_2D computes the nearest point on a triangle in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, double precision P(2), the point whose nearest triangle point
c    is to be determined.
c
c    Output, double precision PN(2), the nearest point to P.
c
c    Output, double precision DIST, the distance from the point to the
c    triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer side_num
      parameter ( side_num = 3 )

      double precision dist
      double precision dist2
      integer i4_wrap
      integer j
      integer jp1
      double precision p(dim_num)
      double precision pn(dim_num)
      double precision pn2(dim_num)
      double precision r8_huge
      double precision t(dim_num,side_num)
      double precision tval
c
c  Find the distance to each of the line segments that make up the edges
c  of the triangle.
c
      dist = r8_huge ( )
      pn(1) = 0.0D+00
      pn(2) = 0.0D+00

      do j = 1, side_num

        jp1 = i4_wrap ( j+1, 1, side_num )

        call segment_point_near_2d ( t(1,j), t(1,jp1), p, 
     &    pn2, dist2, tval )

        if ( dist2 .lt. dist ) then
          dist = dist2
          pn(1) = pn2(1)
          pn(2) = pn2(2)
        end if

      end do

      return
      end
      subroutine triangle_quality_2d ( t, quality )

c*********************************************************************72
c
cc TRIANGLE_QUALITY_2D: "quality" of a triangle in 2D.
c
c  Discussion:
c
c    The quality of a triangle is 2.0 times the ratio of the radius of 
c    the inscribed circle divided by that of the circumscribed circle.  
c    An equilateral triangle achieves the maximum possible quality of 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Adrian Bowyer, John Woodwark,
c    A Programmer's Geometry,
c    Butterworths, 1983.
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Output, double precision QUALITY, the quality of the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision a
      double precision b
      double precision c
      double precision quality
      double precision r8vec_diff_norm
      double precision t(dim_num,3)
c
c  Compute the length of each side.
c
      a = r8vec_diff_norm ( dim_num, t(1,1), t(1,2) )
      b = r8vec_diff_norm ( dim_num, t(1,2), t(1,3) )
      c = r8vec_diff_norm ( dim_num, t(1,3), t(1,1) )

      if ( a * b * c .eq. 0.0D+00 ) then
        quality = 0.0D+00
      else
        quality = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) 
     &    / ( a * b * c )
      end if

      return
      end
      subroutine triangle_right_lattice_point_num_2d ( a, b, n )

c*********************************************************************72
c
cc TRIANGLE_RIGHT_LATTICE_POINT_NUM_2D: count lattice points.
c
c  Discussion:
c
c    The triangle is assumed to be a right triangle which, without loss
c    of generality, has the coordinates:
c
c    ( (0,0), (a,0), (0,b) )
c
c    The routine returns the number of integer lattice points that appear
c    inside the triangle or on its edges or vertices.
c
c    The formula for this function occurred to me (JVB) after some thought,
c    on 06 July 2009.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, B, define the vertices.
c
c    Output, integer N, the number of lattice points.
c
      implicit  none

      integer a
      integer b
      integer i4_gcd
      integer n

      n = ( ( a + 1 ) * ( b + 1 ) + i4_gcd ( a, b ) + 1 ) / 2

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
c    28 March 2011
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
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision P(2,N), random points in the triangle.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer n

      double precision alpha(n)
      integer dim
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
      do j = 1, n
        alpha(j) = sqrt ( alpha(j) )
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
      subroutine triangle_unit_lattice_point_num_2d ( s, n )

c*********************************************************************72
c
cc TRIANGLE_UNIT_LATTICE_POINT_NUM_2D: count lattice points.
c
c  Discussion:
c
c    The triangle is assumed to be the unit triangle:
c
c    ( (0,0), (1,0), (0,1) )
c
c    or a copy of this triangle scaled by an integer S:
c
c    ( (0,0), (S,0), (0,S) ).
c
c    The routine returns the number of integer lattice points that appear
c    inside the triangle or on its edges or vertices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Matthias Beck, Sinai Robins,
c    Computing the Continuous Discretely,
c    Springer, 2006,
c    ISBN13: 978-0387291390,
c    LC: QA640.7.B43.
c
c  Parameters:
c
c    Input, integer S, the scale factor.
c
c    Output, integer N, the number of lattice points.
c
      implicit  none

      integer n
      integer s

      n = ( ( s + 2 ) * ( s + 1 ) ) / 2

      return
      end
      subroutine triangle_xsi_to_xy_2d ( t, xsi, p )

c*********************************************************************72
c
cc TRIANGLE_XSI_TO_XY_2D converts from barycentric to XY coordinates in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, double precision XSI(3), the barycentric coordinates of a point.
c    XSI(1) + XSI(2) + XSI(3) should equal 1, but this is not checked.
c
c    Output, double precision P(2), the XY coordinates of the point.
c
      implicit none

      integer i
      integer j
      double precision p(2)
      double precision t(2,3)
      double precision xsi(3)

      do i = 1, 2
        p(i) = 0.0D+00
        do j = 1, 3
          p(i) = p(i) + t(i,j) * xsi(j)
        end do
      end do

      return
      end
      subroutine triangle_xy_to_xsi_2d ( t, p, xsi )

c*********************************************************************72
c
cc TRIANGLE_XY_TO_XSI_2D converts from XY to barycentric in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the triangle vertices.
c
c    Input, double precision P(2), the XY coordinates of a point.
c
c    Output, double precision XSI(3), the barycentric coordinates of the point.
c    XSI1 + XSI2 + XSI3 should equal 1.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision det
      double precision p(dim_num)
      double precision t(dim_num,3)
      double precision xsi(3)

      det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) 
     &    - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

      xsi(1) = (   ( t(2,2) - t(2,3) ) * ( p(1) - t(1,3) ) 
     &           - ( t(1,2) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

      xsi(2) = ( - ( t(2,1) - t(2,3) ) * ( p(1) - t(1,3) ) 
     &           + ( t(1,1) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

      xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

      return
      end
      subroutine truncated_octahedron_size_3d ( point_num, edge_num,
     &  face_num, face_order_max )

c*********************************************************************72
c
cc TRUNCATED_OCTAHEDRON_SIZE_3D gives "sizes" for a truncated octahedron in 3D.
c
c  Discussion:
c
c    The truncated octahedron is "space-filling".
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

      point_num = 24
      edge_num = 36
      face_num = 14
      face_order_max = 6

      return
      end
      subroutine vector_directions_nd ( dim_num, v, angle )

c*********************************************************************72
c
cc VECTOR_DIRECTIONS_ND returns the direction angles of a vector in ND.
c
c  Discussion:
c
c    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
c    The I-th direction angle is the angle between V and E(I), which is
c    the angle whose cosine is equal to the direction cosine:
c
c      Direction_Cosine(I) = V dot E(I) / |V|.
c
c    If V is the null or zero vector, then the direction cosines and
c    direction angles are undefined, and this routine simply returns
c    zeroes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision V(DIM_NUM), the vector.
c
c    Output, double precision ANGLE(DIM_NUM), the direction angles, in radians,
c    that the vector V makes with the coordinate axes.
c
      implicit none

      integer dim_num

      double precision angle(dim_num)
      integer i
      double precision r8vec_norm
      double precision v(dim_num)
      double precision vnorm
c
c  Get the norm of the vector.
c
      vnorm = r8vec_norm ( dim_num, v )

      if ( vnorm .eq. 0.0D+00 ) then
        do i = 1, dim_num
          angle(i) = 0.0D+00
        end do
        return
      end if

      do i = 1, dim_num
        angle(i) = acos ( v(i) / vnorm )
      end do

      return
      end
      subroutine vector_rotate_2d ( v, angle, w )

c*********************************************************************72
c
cc VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
c
c  Discussion:
c
c    To see why this formula is so, consider that the original point
c    has the form ( R cos Theta, R sin Theta ), and the rotated point
c    has the form ( R cos ( Theta + Angle ), R sin ( Theta + Angle ) ).
c    Now use the addition formulas for cosine and sine to relate
c    the new point to the old one:
c
c      ( W1 ) = ( cos Angle  - sin Angle ) * ( V1 )
c      ( W2 )   ( sin Angle    cos Angle )   ( V2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V(2), the components of the vector to be
c    rotated.
c
c    Input, double precision ANGLE, the angle, in radians, of the rotation
c    to be carried out.  A positive angle rotates the vector in the
c    counter clockwise direction.
c
c    Output, double precision W(2), the rotated vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle
      double precision v(dim_num)
      double precision w(dim_num)

      w(1) = cos ( angle ) * v(1) - sin ( angle ) * v(2)
      w(2) = sin ( angle ) * v(1) + cos ( angle ) * v(2)

      return
      end
      subroutine vector_rotate_3d ( v1, axis, angle, v2 )

c*********************************************************************72
c
cc VECTOR_ROTATE_3D rotates a vector around an axis vector in 3D.
c
c  Discussion:
c
c    Thanks to Cody Farnell for correcting some errors in a previous
c    version of this routinec
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision V1(3), the vector to be rotated.
c
c    Input, double precision AXIS(3), the vector about which the
c    rotation is to be carried out.
c
c    Input, double precision ANGLE, the angle, in radians, of the rotation
c    to be carried out.
c
c    Output, double precision V2(3), the rotated vector.
c
      implicit none

      double precision angle
      double precision axis(3)
      double precision dot
      integer i
      double precision norm
      double precision norm_vn
      double precision normal2(3)
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision v1(3)
      double precision v2(3)
      double precision vn(3)
      double precision vp(3)
      double precision vr(3)
c
c  Compute the length of the rotation axis.
c
      norm = r8vec_norm ( 3, axis )

      if ( norm .eq. 0.0D+00 ) then
        do i = 1, 3
          v2(i) = v1(i)
        end do
        return
      end if
c
c  Compute the dot product of the vector and the (unit) rotation axis.
c
      dot = r8vec_dot_product ( 3, v1, axis ) / norm
c
c  Compute the parallel component of the vector.
c
      do i = 1, 3
        vp(i) = dot * axis(i) / norm
      end do
c
c  Compute the normal component of the vector.
c
      do i = 1, 3
        vn(i) = v1(i) - vp(i)
      end do

      norm_vn = r8vec_norm ( 3, vn )

      if ( norm .eq. 0.0D+00 ) then
        do i = 1, 3
          v2(i) = vp(i)
        end do
        return
      end if

      do i = 1, 3
        vn(i) = vn(i) / norm_vn
      end do
c
c  Compute a second vector, lying in the plane, perpendicular
c  to V1 and VN, and forming a right-handed system.
c
      normal2(1) = axis(2) * vn(3) - axis(3) * vn(2)
      normal2(2) = axis(3) * vn(1) - axis(1) * vn(3)
      normal2(3) = axis(1) * vn(2) - axis(2) * vn(1)

      norm = r8vec_norm ( 3, normal2 )

      do i = 1, 3
        normal2(i) = normal2(i) / norm
      end do
c
c  Rotate the normal component by the angle.
c
      do i = 1, 3
        vr(i) = norm_vn * ( cos ( angle ) * vn(i) 
     &                    + sin ( angle ) * normal2(i) )
      end do
c
c  The rotated vector is the parallel component plus the rotated component.
c
      do i = 1, 3
        v2(i) = vp(i) + vr(i)
      end do

      return
      end
      subroutine vector_rotate_base_2d ( p1, pb, angle, p2 )

c*********************************************************************72
c
cc VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
c
c  Discussion:
c
c    The original vector is assumed to be ( X1-XB, Y1-YB ), and the
c    rotated vector is ( X2-XB, Y2-YB ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision P1(2), the endpoint of the original vector.
c
c    Input, double precision PB(2), the location of the base point.
c
c    Input, double precision ANGLE, the angle, in radians, of the rotation
c    to be carried out.  A positive angle rotates the vector in the
c    counter clockwise direction.
c
c    Output, double precision P2(2), the endpoint of the rotated vector.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision pb(dim_num)

      p2(1) = pb(1) + cos ( angle ) * ( p1(1) - pb(1) ) 
     &              - sin ( angle ) * ( p1(2) - pb(2) )

      p2(2) = pb(2) + sin ( angle ) * ( p1(1) - pb(1) ) 
     &              + cos ( angle ) * ( p1(2) - pb(2) )

      return
      end
      subroutine vector_separation_nd ( dim_num, v1, v2, theta )

c*********************************************************************72
c
cc VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
c
c  Discussion:
c
c    Any two vectors lie in a plane, and are separated by a plane angle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision V1(DIM_NUM), V2(DIM_NUM), the two vectors.
c
c    Output, double precision THETA, the angle between the two vectors.
c
      implicit none

      integer dim_num

      double precision cos_theta
      double precision r8_acos
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision theta
      double precision v1(dim_num)
      double precision v1_norm
      double precision v2(dim_num)
      double precision v2_norm

      v1_norm = r8vec_norm ( dim_num, v1 )

      v2_norm = r8vec_norm ( dim_num, v2 )

      cos_theta = r8vec_dot_product ( dim_num, v1, v2 ) 
     &  / v1_norm / v2_norm

      theta = r8_acos ( cos_theta )

      return
      end
      subroutine vector_unit_nd ( dim_num, v )

c*********************************************************************72
c
cc VECTOR_UNIT_ND normalizes a vector in ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input/output, double precision V(DIM_NUM), the vector to be normalized.
c    On output, V should have unit Euclidean norm.  However, if the input vector
c    has zero Euclidean norm, it is not altered.
c
      implicit none

      integer dim_num

      integer i
      double precision norm
      double precision r8vec_norm
      double precision v(dim_num)

      norm = r8vec_norm ( dim_num, v )

      if ( norm .ne. 0.0D+00 ) then
        do i = 1, dim_num
          v(i) = v(i) / norm
        end do
      end if

      return
      end
      function voxels_dist_l1_nd ( dim_num, v1, v2 )

c*********************************************************************72
c
cc VOXELS_DIST_L1_ND computes the L1 distance between voxels in ND.
c
c  Discussion:
c
c    A voxel is generally a point in 3D space with integer coordinates.
c    There's no reason to stick with 3D, so this routine will handle
c    any dimension.
c
c    We can imagine that, in traveling from V1 to V2, we are allowed to 
c    increment or decrement just one coordinate at a time.  The minimum number 
c    of such changes required is the L1 distance. 
c
c    More formally,
c
c      DIST_L1 ( V1, V2 ) = sum ( 1 <= I <= N ) | V1(I) - V2(I) |
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer V1(DIM_NUM), the voxel that begins the line.
c
c    Input, integer V2(DIM_NUM), the voxel that ends the line.
c
c    Output, integer VOXELS_DIST_L1_ND, the L1 distance 
c    between the voxels.
c
      implicit none

      integer dim_num

      integer i
      integer v1(dim_num)
      integer v2(dim_num)
      integer value
      integer voxels_dist_l1_nd

      value = 0
      do i = 1, dim_num
        value = value + abs ( v1(i) - v2(i) )
      end do

      voxels_dist_l1_nd = value

      return
      end
      subroutine voxels_line_3d ( v1, v2, n, v )

c*********************************************************************72
c
cc VOXELS_LINE_3D computes voxels along a line in 3D.
c
c  Discussion:
c
c    The line itself is defined by two voxels.  The line will begin
c    at the first voxel, and move towards the second.  If the value of
c    N is equal to the L1 distance between the two voxels, then the
c    line will "almost" reach the second voxel.  Depending on the
c    direction, 1, 2 or 3 more steps may be needed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    Daniel Cohen
c    FORTRAN90 version by John Burkardt
c
c  Reference:
c
c    Daniel Cohen,
c    Voxel Traversal along a 3D Line,
c    in Graphics Gems IV,
c    edited by Paul Heckbert,
c    AP Professional, 1994,
c    LC: T385.G6974.
c
c  Parameters:
c
c    Input, integer V1(3), the voxel that begins the line.
c
c    Input, integer V2(3), the voxel that ends the line.
c
c    Input, integer N, the number of voxels to compute.
c
c    Output, integer V(3,N), a sequence of voxels, whose
c    first value is V1 and which proceeds towards V2.
c
      implicit none

      integer n
      integer dim_num
      parameter ( dim_num = 3 )

      integer a(dim_num)
      integer dim
      integer exy
      integer exz
      integer ezy
      integer i
      integer s(dim_num)
      integer v(dim_num,n)
      integer v1(dim_num)
      integer v2(dim_num)

      if ( n .le. 0 ) then
        return
      end if
c
c  Determine the number of voxels on the line.
c
      do dim = 1, dim_num
        s(dim) = sign ( 1, v2(dim) - v1(dim) )
      end do

      do dim = 1, dim_num
        a(dim) = abs ( v2(dim) - v1(dim) )
      end do

      exy = a(2) - a(1)
      exz = a(3) - a(1)
      ezy = a(2) - a(3)
c
c  We start at the starting point.
c
      do dim = 1, dim_num
        v(dim,1) = v1(dim)
      end do

      do i = 2, n

        do dim = 1, dim_num
          v(dim,i) = v(dim,i-1)
        end do

        if ( exy .lt. 0 ) then

          if ( exz .lt. 0 ) then
            v(1,i) = v(1,i) + s(1)
            exy = exy + 2 * a(2)
            exz = exz + 2 * a(3)
          else
            v(3,i) = v(3,i) + s(3)
            exz = exz - 2 * a(1)
            ezy = ezy + 2 * a(2)
          end if

        else if ( ezy .lt. 0 ) then

          v(3,i) = v(3,i) + s(3)
          exz = exz - 2 * a(1)
          ezy = ezy + 2 * a(2)

        else

          v(2,i) = v(2,i) + s(2)
          exy = exy - 2 * a(1)
          ezy = ezy - 2 * a(3)

        end if

      end do

      return
      end
      subroutine voxels_region_3d ( list_max, nx, ny, nz, ishow, 
     &  list_num, list, region_num )

c*********************************************************************72
c
cc VOXELS_REGION_3D arranges contiguous voxels into regions in 3D.
c
c  Discussion:
c
c    On input, the ISHOW array contains zero and nonzero values.  The nonzero
c    values are taken to be active voxels.  On output, the zero voxels remain
c    zero, and all the active voxels have been assigned a value which now
c    indicates membership in a region, or group of contiguous voxels.
c
c    On output, the array LIST contains information about the regions.
c    The last used element of LIST is LIST_NUM.
c
c    The number of elements in region REGION_NUM is NELEM = LIST(LIST_NUM).  
c    The (I,J,K) indices of the last element in this region are in
c    LIST(LIST_NUM-3) through LIST(LIST_NUM-1), and the first element is
c    listed in LIST(LIST_NUM-3*NELEM), LIST(LIST_NUM-3*NELEM+1),
c    LIST(LIST_NUM-3*NELEM+2).
c
c    The number of elements in REGION_NUM-1 is listed in
c    LIST(LIST_NUM-3*NELEM-1), 
c    and the (I,J,K) indices of the these elements are listed there.
c
c    Thanks to Emre Evren for pointing out a hard-to-spot error involving
c    a DO loop that mistakenly read "DO 1 = 1, N".
c
c  Picture:
c
c    Input:
c
c      0  2  0  0 17  0  3
c      0  0  3  0  1  0  4
c      1  0  4  8  8  0  7
c      3  0  6 45  0  0  0
c      3 17  0  5  9  2  5
c
c    Output:
c
c      0  1  0  0  2  0  3
c      0  0  2  0  2  0  3
c      4  0  2  2  2  0  3
c      4  0  2  2  0  0  0
c      4  4  0  2  2  2  2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LIST_MAX, the maximum length of the array 
c    used to list the elements of the regions.
c
c    Input, integer NX, NY, NZ, the number of voxels in the X, Y 
c    and Z directions.
c
c    Input/output, integer ISHOW(NX,NY,NZ).  On input, the only 
c    significance to the entries is whether they are zero or nonzero.  On 
c    output, the nonzero entries have now been revalued so that contiguous 
c    entries have the same value, indicating a grouping into a region.
c
c    Output, integer LIST_NUM, the number of entries of LIST that 
c    were used.  However, if LIST_MAX < LIST_NUM, then there was not enough 
c    space in LIST to store the data properly, and LIST should not be used,
c    although the data in ISHOW should be correct.
c
c    Output, integer LIST(LIST_MAX), contains, in stack form, a 
c    list of the indices of the elements in each region.
c
c    Output, integer REGION_NUM, the number of regions discovered.
c
      implicit none

      integer maxstack
      parameter ( maxstack = 100 )

      integer list_max
      integer nx
      integer ny
      integer nz

      integer i
      integer i2
      integer ibase
      integer ihi
      integer ilo
      integer ishow(nx,ny,nz)
      integer j
      integer j2
      integer jbase
      integer jhi
      integer jlo
      integer k
      integer k2
      integer kbase
      integer khi
      integer klo
      integer list(list_max)
      integer list_num
      integer nabes
      integer ncan
      integer nelements
      integer nstack
      integer region_num
      integer stack(maxstack)
c
c  Reset all nonzero entries of ISHOW to -1.
c
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            if ( ishow(i,j,k) .ne. 0 ) then
              ishow(i,j,k) = -1
            end if

          end do
        end do
      end do
c
c  Start the number of items in the region list at 0.
c
      list_num = 0
c
c  Start the number of regions at 0.
c
      region_num = 0
c
c  The stack begins empty.
c
      nstack = 0
c
c  Search for an unused "ON" voxel from which we can "grow" a new region.
c
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
c
c  We found a voxel that is "ON", and does not belong to any region.
c
            if ( ishow(i,j,k) .eq. -1 ) then
c
c  Increase the number of regions.
c
              region_num = region_num + 1
c
c  Add this voxel to the region.
c
              ishow(i,j,k) = region_num
c
c  Add this voxel to the stack.
c
              if ( maxstack .lt. nstack + 4 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'VOXELS_REGION - Fatal error!'
                write ( *, '(a)' ) '  The internal stack overflowed.'
                write ( *, '(a)' ) '  The algorithm has failed.'
                stop
              end if

              stack(nstack+1) = i
              stack(nstack+2) = j
              stack(nstack+3) = k
              stack(nstack+4) = 1

              nstack = nstack + 4
c
c  Add this voxel to the description of the region.
c
              nelements = 1

              if ( list_num + 3 .le. list_max ) then
                list(list_num+1) = i
                list(list_num+2) = j
                list(list_num+3) = k
              end if

              list_num = list_num + 3

10            continue
c
c  Find all neighbors of BASE that are "ON" but unused.
c  Mark them as belonging to this region, and stack their indices.
c
                ibase = stack(nstack-3)
                jbase = stack(nstack-2)
                kbase = stack(nstack-1)

                ilo = max ( ibase-1, 1 )
                ihi = min ( ibase+1, nx )
                jlo = max ( jbase-1, 1 )
                jhi = min ( jbase+1, ny )
                klo = max ( kbase-1, 1 )
                khi = min ( kbase+1, nz )

                nabes = 0

                do i2 = ilo, ihi
                  do j2 = jlo, jhi
                    do k2 = klo, khi
c
c  We found a neighbor to our current search point, which is "ON" and unused.
c
                      if ( ishow(i2,j2,k2) .eq. -1 ) then
c
c  Increase the number of neighbors.
c
                        nabes = nabes + 1
c
c  Mark the neighbor as belonging to the region.
c
                        ishow(i2,j2,k2) = region_num
c
c  Add the neighbor to the stack.
c
                        if ( maxstack .lt. nstack + 3 ) then
                          write ( *, '(a)' ) ' '
                          write ( *, '(a)' ) 
     &                      'VOXELS_REGION - Fatal error!'
                          write ( *, '(a)' ) 
     &                      '  The internal stack overflowed.'
                          write ( *, '(a)' ) 
     &                      '  The algorithm has failed.'
                          stop
                        end if

                        stack(nstack+1) = i2
                        stack(nstack+2) = j2
                        stack(nstack+3) = k2

                        nstack = nstack + 3
c
c  Add the neighbor to the description of the region.
c
                        nelements = nelements + 1

                        if ( list_num + 3 .le. list_max ) then
                          list(list_num+1) = i2
                          list(list_num+2) = j2
                          list(list_num+3) = k2
                        end if

                        list_num = list_num + 3

                      end if

                    end do
                  end do
                end do
c
c  If any new neighbors were found, take the last one as the basis
c  for a deeper search.
c
                if ( 0 .lt. nabes ) then

                  if ( maxstack .lt. nstack + 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'VOXELS_REGION - Fatal error!'
                    write ( *, '(a)' ) 
     &                '  The internal stack overflowed.'
                    write ( *, '(a)' ) '  The algorithm has failed.'
                    stop
                  end if

                  stack(nstack+1) = nabes
                  nstack = nstack + 1
                  go to 10

                end if
c
c  If the current search point had no new neighbors, drop it from the stack.
c
                ncan = stack(nstack) - 1
                nstack = nstack - 3
                stack(nstack) = ncan
c
c  If there are still any unused candidates at this level, take the
c  last one as the basis for a deeper search.
c
                if ( 0 .lt. stack(nstack) ) then
                  go to 10
                end if
c
c  If there are no more unused candidates at this level, then we need
c  to back up a level in the stack.  If there are any candidates at
c  that earlier level, then we can still do more searching.
c
                nstack = nstack - 1

                if ( nstack .le. 0 ) then
                  go to 20
                end if

              go to 10

20            continue
c
c  If we have exhausted the stack, we have completed this region.
c  Tag the number of elements to the end of the region description list.
c
              list_num = list_num + 1
              if ( list_num .le. list_max ) then
                list(list_num) = nelements
              end if

            end if

          end do
        end do
      end do
c
c  Print some warnings.
c
      if ( list_max .lt. list_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VOXELS_REGION - Warning!'
        write ( *, '(a)' ) 
     &    '  LIST_MAX was too small to list the regions.'
        write ( *, '(a)' ) '  Do not try to use the LIST array!'
        write ( *, '(a)' ) '  The ISHOW data is OK, however.'
      end if

      return
      end
      subroutine voxels_step_3d ( v1, v2, inc, jnc, knc, v3 )

c*********************************************************************72
c
cc VOXELS_STEP_3D computes voxels along a line from a given point in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer V1(3), the coordinates of the base voxel from
c    which the line begins.
c
c    Input, integer V2(3), the coordinates of the current voxel
c    on the line.  For the first call, these might be equal to V1.
c
c    Input, integer INC, JNC, KNC, the increments to the voxels.
c    These values define the direction along which the line proceeds.
c    However, the voxels on the line will typically be incremented
c    by a fractional value of the vector (INC,JNC,KNC), and the
c    result is essentially rounded.
c
c    Output, integer V3(3), the coordinates of the next voxel along
c    the line.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision alpha
      double precision alphai
      double precision alphaj
      double precision alphak
      integer i
      integer inc
      integer jnc
      integer knc
      double precision r8_huge
      integer v1(dim_num)
      integer v2(dim_num)
      integer v3(dim_num)

      do i = 1, dim_num
        v3(i) = v2(i)
      end do
c
c  Assuming for the moment that (I,J,K) can take on real values,
c  points on the line have the form:
c
c    I = V1(2) + alpha * inc
c    J = V1(2) + alpha * jnc
c    K = V1(3) + alpha * knc
c
      if ( inc == 0 .and. jnc == 0 .and. knc == 0 ) then
        return
      end if

      alpha = 0.0D+00
c
c  Compute the smallest ALPHA that will change one of V2(1:3) by +-0.5.
c
      if ( 0 .lt. inc ) then
        alphai = ( dble ( v2(1) - v1(1) ) + 0.5D+00 ) / dble ( inc )
      else if ( inc .lt. 0 ) then
        alphai = ( dble ( v2(1) - v1(1) ) - 0.5D+00 ) / dble ( inc )
      else
        alphai = huge ( alphai )
      end if

      if ( 0 .lt. jnc ) then
        alphaj = ( dble ( v2(2) - v1(2) ) + 0.5D+00 ) / dble ( jnc )
      else if ( jnc .lt. 0 ) then
        alphaj = ( dble ( v2(2) - v1(2) ) - 0.5D+00 ) / dble ( jnc )
      else
        alphaj = huge ( alphaj )
      end if

      if ( 0 .lt. knc ) then
        alphak = ( dble ( v2(3) - v1(3) ) + 0.5D+00 ) / dble ( knc )
      else if ( knc .lt. 0 ) then
        alphak = ( dble ( v2(3) - v1(3) ) - 0.5D+00 ) / dble ( knc )
      else
        alphaj = r8_huge ( )
      end if
c
c  The ALPHA of smallest positive magnitude represents the closest next voxel.
c
      alpha = r8_huge ( )

      if ( 0.0D+00 .lt. alphai ) then
        alpha = min ( alpha, alphai )
      end if

      if ( 0.0D+00 .lt. alphaj ) then
        alpha = min ( alpha, alphaj )
      end if

      if ( 0.0D+00 .lt. alphak ) then
        alpha = min ( alpha, alphak )
      end if
c
c  Move to the new voxel.  Whichever index just made the half
c  step must be forced to take a whole step.
c
      if ( alpha .eq. alphai ) then
        v3(1) = v2(1) + sign ( 1, inc )
        v3(2) = v1(2) + nint ( alpha * dble ( jnc ) )
        v3(3) = v1(3) + nint ( alpha * dble ( knc ) )
      else if ( alpha .eq. alphaj ) then
        v3(1) = v1(1) + nint ( alpha * dble ( inc ) )
        v3(2) = v2(2) + sign ( 1, jnc )
        v3(3) = v1(3) + nint ( alpha * dble ( knc ) )
      else if ( alpha .eq. alphak ) then
        v3(1) = v1(1) + nint ( alpha * dble ( inc ) )
        v3(2) = v1(2) + nint ( alpha * dble ( jnc ) )
        v3(3) = v2(3) + sign ( 1, knc )
      end if

      return
      end
      subroutine xy_to_polar ( xy, r, t )

c*********************************************************************72
c
cc XY_TO_POLAR converts XY coordinates to polar coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision XY(2), the Cartesian coordinates.
c
c    Output, double precision R, T, the radius and angle (in radians).
c
      implicit none

      double precision r
      double precision r8_atan
      double precision t
      double precision xy(2)
      double precision y

      r = sqrt ( xy(1) * xy(1) + xy(2) * xy(2) )

      if ( r .eq. 0.0D+00 ) then
        t = 0.0D+00
      else
        t = r8_atan ( xy(2), xy(1) )
      end if

      return
      end
      subroutine xyz_to_radec ( p, ra, dec )

c*********************************************************************72
c
cc XYZ_TO_RADEC converts (X,Y,Z) to right ascension/declination coordinates.
c
c  Discussion:
c
c    Given an XYZ point, compute its distance R from the origin, and
c    regard it as lying on a sphere of radius R, whose axis is the Z
c    axis.
c
c    The right ascension of the point is the "longitude", measured in hours,
c    between 0 and 24, with the X axis having right ascension 0, and the
c    Y axis having right ascension 6.
c
c    Declination measures the angle from the equator towards the north pole,
c    and ranges from -90 (South Pole) to 90 (North Pole).
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
c    Input, double precision P(3), the coordinates of a point in 3D.
c
c    Output, double precision RA, DEC, the corresponding right ascension
c    and declination.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      double precision dec
      double precision p(dim_num)
      double precision p_norm
      double precision phi
      double precision r8_asin
      double precision r8_atan
      double precision r8vec_norm
      double precision ra
      double precision radians_to_degrees
      double precision theta

      p_norm = r8vec_norm ( 3, p )

      if ( p_norm .eq. 0.0D+00 ) then
        dec = 0.0D+00
        ra = 0.0D+00
        return
      end if

      phi = r8_asin ( p(3) / p_norm )

      if ( cos ( phi ) .eq. 0.0D+00 ) then
        theta = 0.0D+00
      else
        theta = r8_atan ( p(2), p(1) )
      end if

      dec = radians_to_degrees ( phi )
      ra = radians_to_degrees ( theta ) / 15.0D+00

      return
      end
      subroutine xyz_to_rtp ( xyz, r, theta, phi )

c*********************************************************************72
c
cc XYZ_TO_RTP converts (X,Y,Z) to (R,Theta,Phi) coordinates.
c
c  Discussion:
c
c    Given an XYZ point, compute its distance R from the origin, and
c    regard it as lying on a sphere of radius R, whose axis is the Z
c    axis.
c
c    Theta measures the "longitude" of the point, between 0 and 2 PI.
c
c    PHI measures the angle from the "north pole", between 0 and PI.
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
c    Input, double precision XYZ(3), the coordinates of a point in 3D.
c
c    Output, double precision R, THETA, PHI, the radius, longitude and
c    declination of the point.
c
      implicit none

      double precision r
      double precision r8_acos
      double precision r8_atan
      double precision phi
      double precision theta
      double precision xyz(3)

      r = sqrt ( xyz(1) * xyz(1) + xyz(2) * xyz(2) + xyz(3) * xyz(3) )

      if ( r .eq. 0.0D+00 ) then
        theta = 0.0D+00
        phi = 0.0D+00
        return
      end if

      phi = r8_acos ( xyz(3) / r )

      theta = r8_atan ( xyz(2), xyz(1) )

      return
      end
      subroutine xyz_to_tp ( xyz, theta, phi )

c*********************************************************************72
c
cc XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates.
c
c  Discussion:
c
c    Given an XYZ point, regard it as lying on a sphere of radius R,
c    centered at the origin, whose axis is the Z axis.
c
c    We assume that the actual value of R is of no interest, and do
c    not report it.  This is especially appropriate if the point is
c    expected to lie on the unit sphere, for instance.
c
c    THETA measures the "longitude" of the point, between 0 and 2 PI.
c
c    PHI measures the angle from the "north pole", between 0 and PI.
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
c    Input, double precision XYZ(3), the coordinates of a point in 3D.
c
c    Output, double precision THETA, PHI, the longitude and declination
c    of the point.
c
      implicit none

      double precision r
      double precision r8_acos
      double precision r8_atan
      double precision phi
      double precision theta
      double precision xyz(3)

      r = sqrt ( xyz(1) * xyz(1) + xyz(2) * xyz(2) + xyz(3) * xyz(3) )

      if ( r .eq. 0.0D+00 ) then
        theta = 0.0D+00
        phi = 0.0D+00
        return
      end if

      phi = r8_acos ( xyz(3) / r )

      theta = r8_atan ( xyz(2), xyz(1) )

      return
      end
