      program main

c*********************************************************************72
c
cc MAIN is the main program for GEOMETRY_PRB.
c
c  Discussion:
c
c    GEOMETRY_PRB calls sample problems for the GEOMETRY library.
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
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GEOMETRY_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the GEOMETRY library.'

      call test0005 ( )
      call test002 ( )
      call test0477 ( )
      call test0478 ( )
      call test0616 ( )
      call test0617 ( )
      call test0685 ( )
      call test171 ( )
      call test1712 ( )
      call test1788 ( )
      call test1789 ( )
      call test1835 ( )
      call test1836 ( )
      call test203224 ( )
      call test203225 ( )
      call test2101 ( )
      call test2104 ( )
      call test2105 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GEOMETRY_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test0005 ( )

c*********************************************************************72
c
cc TEST0005 tests ANGLE_BOX_2D.
c
c  Discussion:
c
c    Test 1:
c
c      y = 0
c      y = 2x-6
c
c    Test 2:
c
c      y = 0
c      y = 2x-6
c
c    Test 3:
c
c      By setting P1 = P2, we are asking to extend the line
c      y = 2x-6
c      from P3 to P2 through to the other side.
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
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer test_num
      parameter ( test_num = 3 )

      double precision dist
      double precision dist_test(test_num)
      double precision p1(dim_num)
      double precision p1_test(dim_num,test_num)
      double precision p2(dim_num)
      double precision p2_test(dim_num,test_num)
      double precision p3(dim_num)
      double precision p3_test(dim_num,test_num)
      double precision p4(dim_num)
      double precision p5(dim_num)
      integer test

      save dist_test
      save p1_test
      save p2_test
      save p3_test

      data dist_test /
     &  1.0D+00, 1.0D+00, 1.0D+00 /
      data p1_test /
     &  0.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00, 
     &  3.0D+00, 0.0D+00 /
      data p2_test /
     &  3.0D+00, 0.0D+00, 
     &  3.0D+00, 0.0D+00, 
     &  3.0D+00, 0.0D+00 /
      data p3_test /
     &  4.0D+00,  2.0D+00, 
     &  2.0D+00, -2.0D+00, 
     &  2.0D+00, -2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0005'
      write ( *, '(a)' ) '  ANGLE_BOX_2D'
      write ( *, '(a)' ) '  Compute points P4 and P5, normal to '
      write ( *, '(a)' ) '  line through P1 and P2, and'
      write ( *, '(a)' ) '  line through P2 and P3, '
      write ( *, '(a)' ) '  and DIST units from P2.'

      do test = 1, test_num

        p1(1:dim_num) = p1_test(1:dim_num,test)
        p2(1:dim_num) = p2_test(1:dim_num,test)
        p3(1:dim_num) = p3_test(1:dim_num,test)
        dist = dist_test(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  DIST = ', dist
        write ( *, '(a,2g14.6)' ) '  P1:', p1(1:dim_num)
        write ( *, '(a,2g14.6)' ) '  P2:', p2(1:dim_num)
        write ( *, '(a,2g14.6)' ) '  P3:', p3(1:dim_num)
 
        call angle_box_2d ( dist, p1, p2, p3, p4, p5 )
 
        write ( *, '(a,2g14.6)' ) '  P4:', p4(1:dim_num)
        write ( *, '(a,2g14.6)' ) '  P5:', p5(1:dim_num)

      end do

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests ANGLE_DEG_2D and ANGLE_RAD_ND.
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
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )

      double precision angle_deg_2d
      integer angle_num
      parameter ( angle_num = 12 )
      double precision angle_rad_nd
      double precision degrees_to_radians
      integer i
      double precision radians_to_degrees
      double precision temp1
      double precision temp2
      double precision temp3
      double precision thetad
      double precision thetar
      double precision v1(dim_num)
      double precision v2(dim_num)
      double precision v3(dim_num)

      save v1
      save v3

      data v1 / 1.0D+00, 0.0D+00 /
      data v3 / 0.0D+00, 0.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' ) '  ANGLE_DEG_2D computes an angle;'
      write ( *, '(a)' ) '  ANGLE_RAD_ND computes an angle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X  Y  Theta  ATAN2(y, x), ' // 
     &  'ANGLE_RAD_ND, ANGLE_DEG_2D'
      write ( *, '(a)' ) ' '

      do i = 0, angle_num
 
        thetad = dble ( i ) * 360.0D+00 / dble ( angle_num )
        thetar = degrees_to_radians ( thetad )

        v2(1) = cos ( thetar )
        v2(2) = sin ( thetar )
    
        temp1 = radians_to_degrees ( atan2 ( v2(2), v2(1) ) )

        temp2 = angle_rad_nd ( dim_num, v1, v2 )

        temp3 = angle_deg_2d ( v1, v3, v2 )

        write ( *, '(2x,7f10.3)') v2(1:2), thetad, temp1, temp2, temp3
 
      end do
 
      return
      end
      subroutine test0477 ( )

c*********************************************************************72
c
cc TEST0477 tests PARALLELOGRAM_AREA_2D.
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
      implicit none

      double precision area
      double precision p(2,4)

      save p

      data p /
     &  2.0D+00, 7.0D+00, 
     &  5.0D+00, 7.0D+00, 
     &  6.0D+00, 9.0D+00, 
     &  3.0D+00, 9.0D+00  /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0477'
      write ( *, '(a)' ) '  PARALLELOGRAM_AREA_2D finds the area of a'
      write ( *, '(a)' ) '  parallelogram in 2D.'

      call r8mat_transpose_print ( 2, 4, p, '  Vertices:' )

      call parallelogram_area_2d ( p, area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  AREA = ', area

      return
      end
      subroutine test0478 ( )

c*********************************************************************72
c
cc TEST0478 tests PARALLELOGRAM_AREA_3D.
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
      implicit none

      double precision area
      double precision p(3,4)

      save p

      data p /
     & 1.0D+00,       2.0D+00,       3.0D+00, 
     & 2.4142137D+00, 3.4142137D+00, 3.0D+00, 
     & 1.7071068D+00, 2.7071068D+00, 4.0D+00, 
     & 0.2928931D+00, 0.2928931D+00, 4.0D+00  /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0478'
      write ( *, '(a)' ) '  PARALLELOGRAM_AREA_3D finds the area of a'
      write ( *, '(a)' ) '  parallelogram in 3D.'

      call r8mat_transpose_print ( 3, 4, p, '  Vertices:' )

      call parallelogram_area_3d ( p, area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  AREA = ', area

      return
      end
      subroutine test0616 ( )

c*********************************************************************72
c
cc TEST0616 tests PLANE_NORMAL_QR_TO_XYZ and PLANE_NORMAL_XYZ_TO_QR.
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
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 5 )

      double precision dif
      integer i
      integer j
      double precision normal(m)
      double precision pp(m)
      double precision pq(m)
      double precision pr(m)
      double precision qr1(m-1,n)
      double precision qr2(m-1,n)
      integer seed
      double precision t
      double precision xyz(m,n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0616'
      write ( *, '(a)' ) 
     &  '  For a normal plane, with point PP and NORMAL vector,'
      write ( *, '(a)' ) '  and in-plane basis vectors PQ and PR,'
      write ( *, '(a)' ) 
     &  '  PLANE_NORMAL_QR_TO_XYZ converts QR to XYZ coordinates;'
      write ( *, '(a)' ) 
     &  '  PLANE_NORMAL_XYZ_TO_QR converts XYZ to QR coordinates.'
c
c  Choose PP and NORMAL at random.
c
      call r8vec_uniform_01 ( m, seed, pp )

      call r8vec_uniform_01 ( m, seed, normal )
c
c  Compute in-plane basis vectors PQ and PR.
c
      call plane_normal_basis_3d ( pp, normal, pq, pr )
c
c  Choose random Q, R coordinates.
c
      call r8mat_uniform_01 ( m - 1, n, seed, qr1 )
c
c  Convert to XYZ.
c
      call plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr1, xyz )
c
c  Convert XYZ to QR.
c
      call plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz, qr2 )

      dif = 0.0D+00
      do j = 1, n
        t = 0.0
        do i = 1, m - 1
          t = t + ( qr1(i,j) - qr2(i,j) )**2
        end do
        t = sqrt ( t )
        dif = max ( dif, t )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Maximum difference was ', dif

      return
      end
      subroutine test0617 ( )

c*********************************************************************72
c
cc TEST0617 tests PLANE_NORMAL_TETRAHEDRON_INTERSECT.
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
      implicit none

      integer i
      integer j
      integer k
      integer l
      integer int_num
      double precision normal(3)
      double precision pint(3,4)
      double precision pp(3)
      double precision t(3,4)

      data t /
     &  0.0D+00, 0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 0.0D+00, 
     &  0.0D+00, 1.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0617'
      write ( *, '(a)' ) '  PLANE_NORMAL_TETRAHEDRON_INTERSECT'
      write ( *, '(a)' ) 
     &  '  determines the intersection of a plane and tetrahedron.'

      do k = 1, 2

        if ( k .eq. 1 ) then
          normal(1) = 0.0D+00
          normal(2) = 0.0D+00
          normal(3) = 1.0D+00
        else
          normal(1) = 1.0D+00 / sqrt ( 2.0D+00 )
          normal(2) = 1.0D+00 / sqrt ( 2.0D+00 )
          normal(3) = 0.0D+00
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Plane normal vector number ', k
        write ( *, '(a)' ) ' '
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    ( normal(i), i = 1, 3 )

        do l = 0, 6

          do i = 1, 3
            pp(i) = normal(i) * real ( l, kind = 8 ) / 5.0D+00
          end do

          call plane_normal_tetrahedron_intersect ( pp, normal, t, 
     &      int_num, pint )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Point on plane:'
          write ( *, '(a)' ) ' '
          write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      ( pp(i), i = 1, 3 )
          write ( *, '(a)' ) ' '
          write ( *, '(a,i4)' ) 
     &      '  Number of intersection points = ', int_num
          write ( *, '(a)' ) ' '
          do j = 1, int_num
            write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &        j, ( pint(i,j), i = 1, 3 )
          end do

        end do

      end do

      return
      end
      subroutine test0685 ( )

c*********************************************************************72
c
cc TEST0685 tests POLAR_TO_XY and XY_TO_POLAR.
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
      implicit none

      double precision b
      double precision c
      double precision r
      double precision r8_uniform
      integer seed
      double precision t
      integer test
      integer test_num
      parameter ( test_num = 10 )
      double precision xy1(2)
      double precision xy2(2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0685'
      write ( *, '(a)' ) '  POLAR_TO_XY converts (R,Theta) to (X,Y);'
      write ( *, '(a)' ) '  XY_TO_POLAR converts (X,Y) to (R,Theta).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '         X           Y     ',
     &  '===>  R           T   =>      X           Y'
      write ( *, '(a)' ) ' '

      b = -1.0D+00
      c = +1.0D+00
      seed = 123456789

      do test = 1, test_num

        xy1(1) = r8_uniform ( b, c, seed )
        xy1(2) = r8_uniform ( b, c, seed )

        call xy_to_polar ( xy1, r, t )
        call polar_to_xy ( r, t, xy2 )

        write ( *, '(2x,6f12.5)' ) 
     &  xy1(1), xy1(2), r, t, xy2(1), xy2(2)

      end do

      return
      end
      subroutine test171

c*********************************************************************72
c
cc TEST171 tests QUAD_AREA_2D, QUAD_AREA2_2D;
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
      implicit none

      double precision area
      double precision q(2,4)

      save q

      data q /
     &  0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 
     &  1.0D+00, 1.0D+00, 
     &  0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST171'
      write ( *, '(a)' ) '  For a quadrilateral in 2D:'
      write ( *, '(a)' ) '  QUAD_AREA_2D finds the area;'
      write ( *, '(a)' ) '  QUAD_AREA2_2D finds the area;'

      call r8mat_transpose_print ( 2, 4, q, '  The vertices:' )

      call quad_area_2d ( q, area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  QUAD_AREA_2D area is  ', area
  
      call quad_area2_2d ( q, area )

      write ( *, '(a,g14.6)' ) '  QUAD_AREA2_2D area is ', area
 
      return
      end
      subroutine test1712

c*********************************************************************72
c
cc TEST1712 tests QUAD_AREA_3D;
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
      implicit none

      double precision area
      double precision area1
      double precision area2
      integer i
      integer j
      double precision q(3,4)

      double precision t(3,3)

      save q

      data q /
     &  2.0D+00, 2.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00, 0.0D+00, 
     &  1.0D+00, 1.0D+00, 1.0D+00, 
     &  3.0D+00, 3.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1712'
      write ( *, '(a)' ) '  For a quadrilateral in 3D:'
      write ( *, '(a)' ) '  QUAD_AREA_3D finds the area.'

      call r8mat_transpose_print ( 3, 4, q, '  The vertices:' )

      call quad_area_3d ( q, area )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  QUAD_AREA_3D area is     ', area

      do j = 1, 3
        do i = 1, 3
          t(i,j) = q(i,j)
        end do
      end do

      call triangle_area_3d ( t, area1 )

      do j = 1, 2
        do i = 1, 3
          t(i,j) = q(i,j+2)
        end do
      end do

      do i = 1, 3
        t(i,3) = q(i,1)
      end do

      call triangle_area_3d ( t, area2 )
      write ( *, '(a,g14.6)' ) 
     &  '  Sum of TRIANGLE_AREA_3D: ', area1 + area2

      return
      end
      subroutine test1788 ( )

c*********************************************************************72
c
cc TEST17888 tests SIMPLEX_LATTICE_LAYER_POINT_NEXT.
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
      implicit none

      integer n_max
      parameter ( n_max = 4 )

      integer test_num
      parameter ( test_num = 4 )

      integer c(n_max+1)
      integer i
      integer layer
      logical more
      integer n
      integer n_test(test_num)
      integer test
      integer v(n_max)

      save n_test

      data n_test / 1, 2, 3, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1788'
      write ( *, '(a)' ) 
     &  '  SIMPLEX_LATTICE_LAYER_POINT_NEXT returns the next'
      write ( *, '(a)' ) 
     &  '  point in an N-dimensional simplex lattice layer defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '    C(N+1) - 1 <= X(1)/C(1) + X(2)/C(2) ',
     &  '+ ... + X(N)/C(N) <= C(N+1).'

      do test = 1, test_num

        n = n_test(test)

        do i = 1, n
          c(i) = i + 1
        end do
        do i = 1, n
          v(i) = 0
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  N = ', n
        write ( *, '(a,5(2x,i4))') 
     &    '  C =       ', ( c(i), i = 1, n )
        write ( *, '(a)' ) '  '

        do layer = 0, 2

          write ( *, '(a)' ) ' '
          write ( *, '(a,i4)' ) '  Layer ', layer
          write ( *, '(a)' ) ' '

          c(n+1) = layer
          more = .false.
          i = 0

10        continue

            call simplex_lattice_layer_point_next ( n, c, v, more )
            if ( .not. more ) then
              write ( *, '(a)' ) '  No more.'
              go to 20
            end if
            i = i + 1
            write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
          go to 10

20        continue

        end do

      end do

      return
      end
      subroutine test1789 ( )

c*********************************************************************72
c
cc TEST1789 tests SIMPLEX_LATTICE_POINT_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer n_max
      parameter ( n_max = 4 )

      integer c(n_max + 1)
      integer i
      integer j
      logical more
      integer n
      integer n_test(test_num)
      integer test
      integer v(n_max)

      save n_test
      data n_test / 1, 2, 3, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1789'
      write ( *, '(a)' ) 
     &  '  SIMPLEX_LATTICE_POINT_NEXT returns the next lattice'
      write ( *, '(a)' ) 
     &  '  point in an N-dimensional simplex defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    0 <= X(1)/C(1) + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).'

      do test = 1, test_num

        n = n_test(test)

        do i = 1, n + 1
          c(i) = n + 2 - i
        end do
        do i = 1, n
          v(i) = 0
        end do
        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  N = ', n
        write ( *, '(a,5(2x,i4))') 
     &    '  C =       ', ( c(i), i = 1, n + 1 )
        write ( *, '(a)' ) ' '

        i = 0

10      continue

          call simplex_lattice_point_next ( n, c, v, more )

          if ( .not. more ) then
            write ( *, '(a)' ) '  No more.'
            go to 20
          end if
          i = i + 1
          write ( *, '(2x,i4,6x,5(2x,i4))' ) i, ( v(j), j = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test1835 ( )

c*********************************************************************72
c
cc TEST1835 tests SPHERE_EXP2IMP_3D and SPHERE_IMP2EXP_3D.
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
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )

      integer i
      double precision pc(dim_num)
      double precision p1(dim_num)
      double precision p2(dim_num)
      double precision p3(dim_num)
      double precision p4(dim_num)
      double precision r

      save pc
      save p1
      save p2
      save p3
      save p4
      save r

      data pc / 1.0D+00, 2.0D+00, 3.0D+00 /
      data p1 / 4.0D+00, 2.0D+00, 3.0D+00 /
      data p2 / 1.0D+00, 5.0D+00, 3.0D+00 /
      data p3 / 1.0D+00, 2.0D+00, 6.0D+00 /
      data p4 / -2.0D+00, 2.0D+00, 3.0D+00 /
      data r / 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1835'
      write ( *, '(a)' ) 
     &  '  SPHERE_EXP2IMP_3D: explicit sphere => implicit form;'
      write ( *, '(a)' ) 
     &  '  SPHERE_IMP2EXP_3D: implicit sphere => explicit form.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial form of explicit sphere:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,3g14.6)' ) ( p1(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p2(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p3(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p4(i), i = 1, dim_num )

      call sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed form of implicit sphere:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Imputed radius = ', r

      call r8vec_print ( dim_num, pc, '  Imputed center' )

      call sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed form of explicit sphere:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,3g14.6)' ) ( p1(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p2(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p3(i), i = 1, dim_num )
      write ( *, '(2x,3g14.6)' ) ( p4(i), i = 1, dim_num )

      return
      end
      subroutine test1836 ( )

c*********************************************************************72
c
cc TEST1836 tests SPHERE_EXP2IMP_ND.
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
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      double precision p(n,n+1)
      double precision pc(n)
      double precision pc_true(n)
      double precision r
      double precision r_true

      save p
      save pc_true
      save r

      data pc_true / 1.0D+00, 2.0D+00, 3.0D+00 /
      data p / 
     &  4.0D+00, 2.0D+00, 3.0D+00, 
     &  1.0D+00, 5.0D+00, 3.0D+00, 
     &  1.0D+00, 2.0D+00, 6.0D+00, 
     & -2.0D+00, 2.0D+00, 3.0D+00 /
      data r_true / 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1836'
      write ( *, '(a)' ) 
     &  '  SPHERE_EXP2IMP_ND: explicit sphere => implicit form;'

      call r8mat_transpose_print ( n, n + 1, p, 
     &  '  Initial form of explicit sphere:' )

      call sphere_exp2imp_nd ( n, p, r, pc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computed form of implicit sphere:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Imputed radius = ', r
      write ( *, '(a,g14.6)' ) '  True radius =    ', r_true

      call r8vec_print ( n, pc, '  Imputed center' )

      call r8vec_print ( n, pc_true, '  True center' )

      return
      end
      subroutine test203224 ( )

c*********************************************************************72
c
cc TEST203224 tests TETRAHEDRON_LATTICE_LAYER_POINT_NEXT.
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
      implicit none

      integer n
      parameter ( n = 3 )

      integer c(n+1)
      integer i
      integer j
      integer layer
      logical more
      integer v(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST203224'
      write ( *, '(a)' ) 
     &  '  TETRAHEDRON_LATTICE_LAYER_POINT_NEXT returns the next'
      write ( *, '(a)' ) 
     &  '  point in a tetrahedron lattice layer defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    C(4) - 1 < X(1)/C(1) + X(2)/C(2) +X(3)/C(3) .le. C(4).'

      c(1) = 2
      c(2) = 3
      c(3) = 4
      v(1) = 0
      v(2) = 0
      v(3) = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  N = ', n
      write ( *, '(a,4(2x,i4))' ) '  C =       ', ( c(i), i = 1, n )

      do layer = 0, 2

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Layer ', layer
        write ( *, '(a)' ) ' '

        c(4) = layer
        more = .false.
        i = 0

        do
          call tetrahedron_lattice_layer_point_next ( c, v, more )
          if ( .not. more ) then
            write ( *, '(a)' ) '  No more.'
            exit
          end if
          i = i + 1
          write ( *, '(2x,i4,6x,10(2x,i4))' ) i, ( v(j), j = 1, n )

        end do

      end do

      return
      end
      subroutine test203225 ( )

c*********************************************************************72
c
cc TEST203225 tests TETRAHEDRON_LATTICE_POINT_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer c(n+1)
      integer i
      integer j
      logical more
      integer v(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST203225'
      write ( *, '(a)' ) 
     &  '  TETRAHEDRON_LATTICE_POINT_NEXT returns the next lattice'
      write ( *, '(a)' ) '  point in a tetrahedron defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    0 <= X(1)/C(1) + X(2)/C(2) + X(3)/C(3) <= C(4).'

      do i = 1, n + 1
        c(i) = n + 2 - i
      end do
      do i = 1, n
        v(i) = 0
      end do
      more = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  N = ', n
      write ( *, '(a,2x,i4,2x,i4,2x,i4,2x,i4)'  )
     &  '  C =       ', ( c(i), i = 1, n + 1 )
      write ( *, '(a)' ) ' '

      i = 0

10    continue

        call tetrahedron_lattice_point_next ( c, v, more )
        if ( .not. more ) then
          write ( *, '(a)' ) '  No more.'
          go to 20
        end if
        i = i + 1
        write ( *, '(2x,i4,6x,10(2x,i4))' ) i, ( v(j), j = 1, n )
      go to 10

20    continue

      return
      end
      subroutine test2101 ( )

c*********************************************************************72
c
cc TEST2101 tests TRIANGLE_CIRCUMCENTER_2D and others.
c
c  Discussion:
c
c    The functions tested include
c    * TRIANGLE_CIRCUMCENTER_2D;
c    * TRIANGLE_CIRCUMCENTER_2D_2;
c    * TRIANGLE_CIRCUMCENTER.
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
      implicit none

      integer m
      parameter ( m = 2 )
      integer test_num
      parameter ( test_num = 4 )

      integer i
      integer j
      double precision pc(m)
      double precision t(m,3)
      double precision t_test(m,3,test_num)
      integer test

      save t_test

      data t_test /
     &       10.0D+00,  5.0D+00, 
     &       11.0D+00,  5.0D+00, 
     &       10.0D+00,  6.0D+00, 
     &       10.0D+00,  5.0D+00, 
     &       11.0D+00,  5.0D+00, 
     &       10.5D+00,  5.86602539D+00, 
     &       10.0D+00,  5.0D+00, 
     &       11.0D+00,  5.0D+00, 
     &       10.5D+00, 15.0D+00, 
     &       10.0D+00,  5.0D+00, 
     &       11.0D+00,  5.0D+00, 
     &      20.0D+00,   7.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST2101'
      write ( *, '(a)' ) 
     &  '  For a triangle in 2D, the circumenter can be computed by:'
      write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER_2D;'
      write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER_2D_2;'
      write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER (any dimension);'

      do test = 1, test_num

        do j = 1, 3
          do i = 1, m
            t(i,j) = t_test(i,j,test)
          end do
        end do

        call r8mat_transpose_print ( m, 3, t, '  Triangle vertices:' )

        call triangle_circumcenter_2d ( t, pc )

        call r8vec_print ( m, pc, 
     &    '  Circumcenter by TRIANGLE_CIRCUMCENTER_2D:' )

        call triangle_circumcenter_2d_2 ( t, pc )

        call r8vec_print ( m, pc, 
     &    '  Circumcenter by TRIANGLE_CIRCUMCENTER_2D_2:' )

        call triangle_circumcenter ( m, t, pc )

        call r8vec_print ( m, pc, 
     &    '  Circumcenter by TRIANGLE_CIRCUMCENTER:' )

      end do

      return
      end
      subroutine test2104 ( )

c*********************************************************************72
c
cc TEST2104 tests TRIANGLE_LATTICE_LAYER_POINT_NEXT.
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
      implicit none

      integer n 
      parameter ( n = 2 )

      integer c(n+1)
      integer i
      integer j
      integer layer
      logical more
      integer v(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST2104'
      write ( *, '(a)' ) 
     &  '  TRIANGLE_LATTICE_LAYER_POINT_NEXT returns the next'
      write ( *, '(a)' ) 
     &  '  point in a triangle lattice layer defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    C(3) - 1 < X(1)/C(1) + X(2)/C(2) .le. C(3).'

      c(1) = 2
      c(2) = 3
      v(1) = 0
      v(2) = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  N = ', n
      write ( *, '(a,4(2x,i4))' ) '  C =       ', ( c(i), i = 1, n )

      do layer = 0, 4

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Layer ', layer
        write ( *, '(a)' ) ' '

        c(3) = layer
        more = .false.
        i = 0

        do
          call triangle_lattice_layer_point_next ( c, v, more )
          if ( .not. more ) then
            write ( *, '(a)' ) '  No more.'
            exit
          end if
          i = i + 1
          write ( *, '(2x,i4,6x,10(2x,i4))' ) i, ( v(j), j = 1, n )

        end do

      end do

      return
      end
      subroutine test2105 ( )

c*********************************************************************72
c
cc TEST2105 tests TRIANGLE_LATTICE_POINT_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      integer c(n+1)
      integer i
      integer j
      logical more
      integer v(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST2105'
      write ( *, '(a)' ) 
     &  '  TRIANGLE_LATTICE_POINT_NEXT returns the next lattice'
      write ( *, '(a)' ) '  point in a triangle defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    0 <= X(1)/C(1) + X(2)/C(2) <= C(3).'

      do i = 1, n + 1
        c(i) = n + 2 - i
      end do
      do i = 1, n
        v(i) = 0
      end do
      more = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  N = ', n
      write ( *, '(a,2x,i4,2x,i4,2x,i4)' ) 
     &  '  C =       ', ( c(i), i = 1, n + 1 )
      write ( *, '(a)' ) ' '

      i = 0

10    continue

        call triangle_lattice_point_next ( c, v, more )

        if ( .not. more ) then
          write ( *, '(a)' ) '  No more.'
          go to 20
        end if

        i = i + 1
        write ( *, '(2x,i4,6x,10(2x,i4))' ) i, ( v(j), j = 1, n )

      go to 10

20    continue

      return
      end
