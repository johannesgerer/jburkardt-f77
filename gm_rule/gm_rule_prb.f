      program main

c*********************************************************************72
c
cc MAIN is the main program for GM_RULE_PRB.
c
c  Discussion:
c
c    GM_RULE_PRB calls a set of problems for GM_RULE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GM_RULE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the GM_RULE library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GM_RULE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer vertex_num
      parameter ( vertex_num = dim_num + 1 )
      integer point_num
      parameter ( point_num = 10 )

      integer i
      integer j
      double precision phy(dim_num,point_num)
      double precision phy_unit(dim_num,dim_num+1)
      double precision ref(dim_num,point_num)
      integer seed
      double precision t(dim_num,vertex_num)
      double precision t_unit(dim_num,vertex_num)

      save t
      save t_unit

      data t /
     &  1.0D+00, 1.0D+00, 
     &  3.0D+00, 1.0D+00, 
     &  2.0D+00, 5.0D+00 /
      data t_unit /
     &  0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 
     &  0.0D+00, 1.0D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  SIMPLEX_UNIT_TO_GENERAL maps points in the unit'
      write ( *, '(a)' ) '  simplex to a general simplex.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Here we consider a simplex in 2D, a triangle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vertices of the general triangle are:'
      write ( *, '(a)' ) ' '
      do j = 1, vertex_num
        write ( *, '(2x,f8.4,2x,f8.4)' ) ( t(i,j), i = 1, dim_num )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   (  XSI     ETA )   ( X       Y  )'
      write ( *, '(a)' ) ' '

      call simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit, 
     &  phy_unit )

      do j = 1, dim_num + 1

        write ( *, '(2x,2f8.4,2x,2f8.4)' ) 
     &    ( t_unit(i,j), phy_unit(i,j), i = 1, dim_num )

      end do

      call simplex_unit_sample ( dim_num, point_num, seed, ref )

      call simplex_unit_to_general ( dim_num, point_num, t, ref, phy )

      do j = 1, point_num

        write ( *, '(2x,2f8.4,2x,2f8.4)' ) 
     &    ( ref(i,j), phy(i,j), i = 1, dim_num )

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests SIMPLEX_UNIT_TO_GENERAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer vertex_num
      parameter ( vertex_num = dim_num + 1 )
      integer point_num
      parameter ( point_num = 10 )

      integer i
      integer j
      double precision phy(dim_num,point_num)
      double precision phy_unit(dim_num,dim_num+1)
      double precision ref(dim_num,point_num)
      integer seed
      double precision t(dim_num,vertex_num)
      double precision t_unit(dim_num,vertex_num)

      save t
      save t_unit

      data t /
     &  1.0D+00, 1.0D+00, 1.0D+00, 
     &  3.0D+00, 1.0D+00, 1.0D+00, 
     &  1.0D+00, 4.0D+00, 1.0D+00, 
     &  1.0D+00, 1.0D+00, 5.0D+00 /
      data t_unit /
     &  0.0D+00, 0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 0.0D+00, 
     &  0.0D+00, 1.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00, 1.0D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  SIMPLEX_UNIT_TO_GENERAL'
      write ( *, '(a)' ) 
     &  '    maps points in the unit simplex to a general simplex.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Here we consider a simplex in 3D, a tetrahedron.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vertices are:'
      write ( *, '(a)' ) ' '
      do j = 1, vertex_num
        write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4)' ) 
     &    ( t(i,j), i = 1, dim_num )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   (  XSI     ETA     MU )    ( X       Y       Z )'
      write ( *, '(a)' ) ' '

      call simplex_unit_to_general ( dim_num, dim_num + 1, t, t_unit, 
     &  phy_unit )

      do j = 1, dim_num + 1

        write ( *, '(2x,3f8.4,2x,3f8.4)' ) 
     &    ( t_unit(i,j), phy_unit(i,j), i = 1, dim_num )

      end do

      call simplex_unit_sample ( dim_num, point_num, seed, ref )

      call simplex_unit_to_general ( dim_num, point_num, t, ref, phy )

      do j = 1, point_num

        write ( *, '(2x,3f8.4,2x,3f8.4)' ) 
     &    ( ref(i,j), phy(i,j), i = 1, dim_num )

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests GM_RULE_SIZE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer dim_num
      integer dim_num_test(test_num)
      integer degree
      integer point_num
      integer rule
      integer test

      save dim_num_test

      data dim_num_test / 2, 3, 5, 10 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  GM_RULE_SIZE returns POINT_NUM, the number of points'
      write ( *, '(a)' ) 
     &  '  associated with a Grundmann-Moeller quadrature rule'
      write ( *, '(a)' ) '  for the unit simplex of dimension DIM_NUM'
      write ( *, '(a)' ) '  with rule index RULE'
      write ( *, '(a)' ) '  and degree of exactness DEGREE = 2*RULE+1.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   DIM_NUM      RULE    DEGREE POINT_NUM'

      do test = 1, test_num

        dim_num = dim_num_test(test)

        write ( *, '(a)' ) ' '

        do rule = 0, 5

          call gm_rule_size ( rule, dim_num, point_num )
          degree = 2 * rule + 1

          write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) 
     &      dim_num, rule, degree, point_num

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests GM_RULE_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer point_num_max
      parameter ( point_num_max = 15 )

      integer i
      integer point
      integer point_num
      integer rule
      double precision w(point_num_max)
      double precision x(dim_num,point_num_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  GM_RULE_SET determines the weights and abscissas'
      write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
      write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
      write ( *, '(a)' ) '  using a rule of in index RULE,'
      write ( *, '(a)' ) 
     &  '  which will have degree of exactness 2*RULE+1.'

      rule = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num
      write ( *, '(a,i8)' ) '  RULE = ', rule
      write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

      call gm_rule_size ( rule, dim_num, point_num )

      call gm_rule_set ( rule, dim_num, point_num, w, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     POINT        W             X             Y             Z'
      write ( *, '(a)' ) ' '

      do point = 1, point_num
        write ( *, '(2x,i8,2x,f12.6,2x,f12.6,2x,f12.6,2x,f12.6)' ) 
     &    point, w(point), ( x(i,point), i = 1, dim_num )
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests GM_RULE_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )
      integer dim_max
      parameter ( dim_max = 10 )
      integer point_max
      parameter ( point_max = 4368 )

      integer dim_num
      integer dim_num_test(test_num)
      integer i
      integer point_num
      double precision r8vec_sum
      integer rule
      integer test
      double precision w(point_max)
      double precision w_sum
      double precision x(dim_max,point_max)

      save dim_num_test

      data dim_num_test / 2, 3, 5, 10 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  GM_RULE_SET determines the weights and abscissas'
      write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
      write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
      write ( *, '(a)' ) '  using a rule of in index RULE,'
      write ( *, '(a)' ) 
     &  '  which will have degree of exactness 2*RULE+1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  In this test, we compute various rules, and simply'
      write ( *, '(a)' ) 
     &  '  report the number of points, and the sum of weights.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   DIM_NUM      RULE    POINT_NUM  WEIGHT SUM'

      do test = 1, test_num

        dim_num = dim_num_test(test)

        write ( *, '(a)' ) ' '

        do rule = 0, 5

          call gm_rule_size ( rule, dim_num, point_num )

          call gm_rule_set ( rule, dim_num, point_num, w, x )

          w_sum = r8vec_sum ( point_num, w )

          write ( *, '(2x,i8,2x,i8,2x,i8,2x,g24.16)' ) 
     &      dim_num, rule, point_num, w_sum

        end do

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests GM_RULE_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer point_max
      parameter ( point_max = 15 )

      integer degree
      integer i
      integer point
      integer point_num
      integer rule
      double precision w(point_max)
      character * ( 12 ) w_file
      integer w_unit
      double precision x(dim_num,point_max)
      character * ( 12 ) x_file
      integer x_unit

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) 
     &  '  GM_RULE_SET determines the weights and abscissas'
      write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
      write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
      write ( *, '(a)' ) '  using a rule of in index RULE,'
      write ( *, '(a)' ) 
     &  '  which will have degree of exactness 2*RULE+1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this test, we write a rule to a file.'

      rule = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num
      write ( *, '(a,i8)' ) '  RULE = ', rule
      write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

      call gm_rule_size ( rule, dim_num, point_num )

      call gm_rule_set ( rule, dim_num, point_num, w, x )

      call get_unit ( w_unit )

      write ( w_file, '(a2,i1,a,i1,a7)' ) 
     &  'gm', rule, '_', dim_num, 'd_w.txt'

      open ( unit = w_unit, file = w_file, status = 'replace' )

      do point = 1, point_num
        write ( w_unit, '(f20.16)' ) w(point)
      end do

      close ( unit = w_unit )

      call get_unit ( x_unit )

      write ( x_file, '(a2,i1,a,i1,a7)' ) 
     &  'gm', rule, '_', dim_num, 'd_x.txt'

      open ( unit = x_unit, file = x_file, status = 'replace' )

      do point = 1, point_num
        write ( x_unit, '(3f20.16)' ) ( x(i,point), i = 1, dim_num )
      end do

      close ( unit = x_unit )

      write ( *, '(a,i2,a)' ) '  Wrote rule ', rule, ' to "' 
     &  // trim ( w_file ) // '" and "' // trim ( x_file ) // '".'

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests GM_RULE_SET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer point_max
      parameter ( point_max = 84 )

      integer degree
      integer degree_max
      parameter ( degree_max = 4 )
      integer expon(dim_num)
      integer h
      integer i
      double precision mono(point_max)
      logical more
      integer point
      integer point_num
      double precision quad
      double precision quad_error
      integer rule
      integer rule_max
      parameter ( rule_max = 3 )
      integer t
      double precision w(point_max)
      double precision x(dim_num,point_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  GM_RULE_SET determines the weights and abscissas'
      write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
      write ( *, '(a)' ) '  the DIM_NUM dimensional simplex,'
      write ( *, '(a)' ) '  using a rule of in index RULE,'
      write ( *, '(a)' ) 
     &  '  which will have degree of exactness 2*RULE+1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  In this test, look at all the monomials up to'
      write ( *, '(a)' ) 
     &  '  some maximum degree, choose a few low order rules'
      write ( *, '(a)' ) 
     &  '  and determine the quadrature error for each.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Here we use DIM_NUM = ', dim_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Rule     Order     Quad_Error'
      write ( *, '(a)' ) ' '

      do degree = 0, degree_max

        more = .false.

10      continue

          call comp_next ( degree, dim_num, expon, more, h, t )

          write ( *, '(a)' ) ' '
          write ( *, '(a,i1,a,i1,a,i1,a,i1,a,i1)' ) 
     &      '  F(X) = X1^', expon(1), ' * X2^', expon(2), 
     &      ' * X3^', expon(3), ' * X4^', expon(4), ' * X5^', expon(5)

          write ( *, '(a)' ) ' '

          do rule = 0, rule_max

            call gm_rule_size ( rule, dim_num, point_num )

            call gm_rule_set ( rule, dim_num, point_num, w, x )

            call simplex_unit_monomial_quadrature ( dim_num, expon, 
     &        point_num, x, w, quad_error )

            write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) 
     &        rule, point_num, quad_error

          end do

          if ( .not. more ) then
            go to 20
          end if

        go to 10

20      continue

      end do

      return
      end
