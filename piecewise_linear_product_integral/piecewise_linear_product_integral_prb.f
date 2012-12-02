      program main

c*********************************************************************72
c
cc MAIN is the main program for PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB.
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
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' )
     &  '  Test the PIECEWISE_LINEAR_PRODUCT_INTEGRAL library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
c
c  Discussion:
c
c    For the first test, we use the same single "piece" for both F and G.
c    Hence, we are actually integrating X^2.
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
      implicit none

      integer f_num
      integer g_num

      parameter ( f_num = 2 )
      parameter ( g_num = 2 )

      double precision a
      double precision b
      double precision exact
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision g_v(g_num)
      double precision g_x(g_num)
      integer i
      double precision integral

      save f_v
      save f_x
      save g_v
      save g_x

      data f_v / 0.0D+00, 5.0D+00 /
      data f_x / 0.0D+00, 5.0D+00 /
      data g_v / 0.0D+00, 5.0D+00 /
      data g_x / 0.0D+00, 5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a,a)' )
     &  '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL',
     &  ' on a very simple problem.'
      write ( *, '(a)' )
     &  '  F and G are both defined over a single common'
      write ( *, '(a)' ) '  interval, so that F(X) = G(X) = X.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '           A           B      Integral        Exact'
      write ( *, '(a)' ) ' '

      a = 1.0D+00
      do i = 1, 5
        b = dble ( i )
        call piecewise_linear_product_integral ( a, b, f_num, f_x,
     &    f_v, g_num, g_x, g_v, integral )
        exact = ( b * b * b - a * a * a ) / 3.0D+00
        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' )
     &    a, b, integral, exact
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
c
c  Discussion:
c
c    For this test, we use multiple "pieces" for both F and G,
c    but we define the values so that we are still actually integrating X^2.
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
      implicit none

      integer f_num
      integer g_num

      parameter ( f_num = 3 )
      parameter ( g_num = 4 )

      double precision a
      double precision b
      double precision exact
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision g_v(g_num)
      double precision g_x(g_num)
      integer i
      double precision integral

      save f_v
      save f_x
      save g_v
      save g_x

      data f_v / 0.0D+00, 2.0D+00, 5.0D+00 /
      data f_x / 0.0D+00, 2.0D+00, 5.0D+00 /
      data g_v / 0.0D+00, 1.5D+00, 3.0D+00, 5.0D+00 /
      data g_x / 0.0D+00, 1.5D+00, 3.0D+00, 5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.'
      write ( *, '(a)' )
     &  '  F and G are both defined over separate, multiple'
      write ( *, '(a)' )
     &  '  intervals, but still true that F(X) = G(X) = X.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '           A           B      Integral        Exact'
      write ( *, '(a)' ) ' '

      a = 1.0D+00
      do i = 1, 5
        b = dble ( i )
        call piecewise_linear_product_integral ( a, b, f_num, f_x,
     &    f_v, g_num, g_x, g_v, integral )
        exact = ( b * b * b - a * a * a ) / 3.0D+00
        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' )
     &    a, b, integral, exact
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
c
c  Discussion:
c
c    For this test, F(X) and G(X) are piecewise linear interpolants to
c    SIN(X) and 2 * COS(X), so we know the exact value of the integral
c    of the product of the original functions, but this is only an estimate
c    of the exact value of the integral of the product of the interpolants.
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
      implicit none

      integer f_num
      integer g_num

      parameter ( f_num = 11 )
      parameter ( g_num = 31 )

      double precision a
      double precision b
      double precision exact
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision g_v(g_num)
      double precision g_x(g_num)
      integer i
      double precision integral
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision quad
      integer quad_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.'
      write ( *, '(a)' ) '  F and G are defined over separate, multiple'
      write ( *, '(a)' ) '  intervals.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X) interpolates SIN(X),'
      write ( *, '(a)' ) '  G(X) interpolates 2*COS(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We compare:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  INTEGRAL, our value for the integral,'
      write ( *, '(a)' )
     &  '  QUAD, a quadrature estimate for the integral, and'
      write ( *, '(a)' )
     &  '  CLOSE, the value of the integral of 2*COS(X)*SIN(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '           A           B      Integral',
     &  '        Quad            Close'
      write ( *, '(a)' ) ' '

      do i = 1, f_num
        f_x(i) = ( dble ( f_num - i     ) * 0.0D+00
     &           + dble (         i - 1 ) * pi )
     &           / dble ( f_num     - 1 )
        f_v(i) = sin ( f_x(i) )
      end do

      do i = 1, g_num
        g_x(i) = ( dble ( g_num - i     ) * 0.0D+00
     &           + dble (         i - 1 ) * pi )
     &           / dble ( g_num     - 1 )
        g_v(i) = 2.0D+00 * cos ( g_x(i) )
      end do

      a = 0.0D+00
      do i = 0, 6
        b = dble ( i ) * pi / 6.0D+00
        call piecewise_linear_product_integral ( a, b, f_num, f_x,
     &    f_v, g_num, g_x, g_v, integral )
        exact = - ( cos ( 2.0D+00 * b ) - cos ( 2.0D+00 * a ) )
     &    / 2.0D+00
        quad_num = 2000
        call piecewise_linear_product_quad ( a, b, f_num, f_x, f_v,
     &    g_num, g_x, g_v, quad_num, quad )
        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &    a, b, integral, quad, exact
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
c
c  Discussion:
c
c    For this test, we compute the integrals of a hat function with itself,
c    and a hat function with its neighbor.
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
      implicit none

      integer f_num
      integer g_num

      parameter ( f_num = 3 )
      parameter ( g_num = 3 )

      double precision a
      double precision b
      double precision exact
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision g_v(g_num)
      double precision g_x(g_num)
      integer i
      double precision integral

      save f_v
      save f_x
      save g_v
      save g_x

      data f_v / 0.0D+00, 1.0D+00, 0.0D+00 /
      data f_x / 0.0D+00, 1.0D+00, 2.0D+00 /
      data g_v / 1.0D+00, 0.0D+00, 0.0D+00 /
      data g_x / 0.0D+00, 1.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL.'
      write ( *, '(a)' ) '  The nodes are at 0, 1, and 2.'
      write ( *, '(a)' ) '  F(X) = ( 0, 1, 0 ).'
      write ( *, '(a)' ) '  G(X) = ( 1, 0, 0 ).'
      write ( *, '(a)' ) ' '

      a = 0.0D+00
      b = 2.0D+00

      call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v,
     &  f_num, f_x, f_v, integral )

      write ( *, '(a,g14.6)' ) '  Integral F(X) * F(X) dx = ', integral

      call piecewise_linear_product_integral ( a, b, f_num, f_x, f_v,
     &  g_num, g_x, g_v, integral )

      write ( *, '(a,g14.6)' ) '  Integral F(X) * G(X) dx = ', integral

      call piecewise_linear_product_integral ( a, b, g_num, g_x, g_v,
     &  g_num, g_x, g_v, integral )

      write ( *, '(a,g14.6)' ) '  Integral G(X) * G(X) dx = ', integral

      return
      end
