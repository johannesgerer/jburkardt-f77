      program main

c*********************************************************************72
c
cc MAIN is the main program for HERMITE_CUBIC_PRB.
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
      implicit none

      call timestamp ( );
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the HERMITE_CUBIC library.'

      call hermite_cubic_test01 ( )
      call hermite_cubic_test02 ( )
      call hermite_cubic_test03 ( )
      call hermite_cubic_test04 ( )
      call hermite_cubic_test05 ( )
      call hermite_cubic_test06 ( )
      call hermite_cubic_test07 ( )
      call hermite_cubic_test08 ( )
      call hermite_cubic_test09 ( )

      call hermite_cubic_test10 ( )
      call hermite_cubic_test11 ( )
      call hermite_cubic_test12 ( )
      call hermite_cubic_test13 ( )
      call hermite_cubic_test14 ( )
      call hermite_cubic_test15 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( );

      stop
      end
      subroutine hermite_cubic_test01 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST01 tests HERMITE_CUBIC_VALUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1 )

      double precision d(n)
      double precision d1
      double precision d2
      double precision f(n)
      double precision f1
      double precision f2
      integer i
      integer j
      double precision s(n)
      double precision t(n)
      double precision x(n)
      integer x_interval
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST01:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.'
      write ( *, '(a)' ) '  Try out four sets of data:'
      write ( *, '(a)' ) 
     &  '  (F1,D1,F2,D2) = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)'
      write ( *, '(a)' ) '  on [0,1] and [1.0,-2.0] (interval reversed)'

      do x_interval = 1, 2

        if ( x_interval .eq. 1 ) then
          x1 = 0.0D+00
          x2 = 1.0D+00
        else
          x1 = 1.0D+00
          x2 = -2.0D+00
        end if

        do i = 1, 4

          f1 = 0.0D+00
          d1 = 0.0D+00
          f2 = 0.0D+00
          d2 = 0.0D+00

          if ( i .eq. 1 ) then
            f1 = 1.0D+00
          else if ( i .eq. 2 ) then
            d1 = 1.0D+00
          else if ( i .eq. 3 ) then
            f2 = 1.0D+00
          else if ( i .eq. 4 ) then
            d2 = 1.0D+00
          end if

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '    J      X           F           D'
          write ( *, '(a)' ) ' '

          do j = - 3, 12

            x = ( dble ( 10 - j ) * x1   
     &          + dble (      j ) * x2 ) 
     &          / dble ( 10     )

            call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, 1, x, 
     &        f, d, s, t )

            if ( j .eq. 0 ) then
              write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &          '*Data', x1, f1, d1
            end if
            write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &        j, x(1), f(1), d(1)
            if ( j .eq. 10 ) then
              write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &          '*Data', x2, f2, d2
            end if

          end do

        end do

      end do

      return
      end
      subroutine hermite_cubic_test02 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST02 tests HERMITE_CUBIC_VALUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1 )

      double precision d(n)
      double precision dc
      double precision d1
      double precision d2
      double precision f(n)
      double precision fc
      double precision f1
      double precision f2
      integer j
      double precision s(n)
      double precision s1
      double precision s2
      double precision sc
      double precision t(n)
      double precision t1
      double precision t2
      double precision tc
      double precision x(n)
      integer x_interval
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST02:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.'
      write ( *, '(a)' ) '  Try out data from a cubic function:'
      write ( *, '(a)' ) '  on [0,10] and [-1.0,1.0] and [0.5,0.75]'

      do x_interval = 1, 3

        if ( x_interval .eq. 1 ) then
          x1 = 0.0D+00
          x2 = 10.0D+00
        else if ( x_interval .eq. 2 ) then
          x1 = -1.0D+00
          x2 = +1.0D+00
        else if ( x_interval .eq. 3 ) then
          x1 = 0.5D+00
          x2 = 0.75D+00
        end if

        call cubic_value ( x1, f1, d1, s1, t1 )
        call cubic_value ( x2, f2, d2, s2, t2 )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '    J      X           F           D           S           T'

        do j = - 3, 12

          x(1) = ( ( 10 - j ) * x1   
     &           +        j   * x2 ) 
     &             / 10.0D+00

          call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, 
     &      d, s, t )

          call cubic_value ( x, fc, dc, sc, tc )

          write ( *, '(a)' ) ' '
          if ( j .eq. 0 ) then
            write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &        '*Data', x1, f1, d1
          end if
          write ( *, 
     &      '(a,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      'Exact', x(1), fc, dc, sc, tc
          write ( *, 
     &      '(2x,i3,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      j, x(1), f(1), d(1), s(1), t(1)
          if ( j .eq. 10 ) then
            write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &        '*Data', x2, f2, d2
          end if

        end do

      end do

      return
      end
      subroutine hermite_cubic_test03 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST03 tests HERMITE_CUBIC_INTEGRATE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision f1
      double precision f2
      integer j
      double precision q_computed
      double precision q_exact
      double precision s1
      double precision s2
      double precision t1
      double precision t2
      integer x_interval
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST03:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic'
      write ( *, '(a)' ) '  polynomial from A to B.'

      do x_interval = 1, 3

        if ( x_interval .eq. 1 ) then
          x1 = 0.0D+00
          x2 = 10.0D+00
        else if ( x_interval .eq. 2 ) then
          x1 = -1.0D+00
          x2 = +1.0D+00
        else if ( x_interval .eq. 3 ) then
          x1 = 0.5D+00
          x2 = 0.75D+00
        end if

        call cubic_value ( x1, f1, d1, s1, t1 )
        call cubic_value ( x2, f2, d2, s2, t2 )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '                                 Exact       Computed'
        write ( *, '(a)' ) 
     &    '    J      A           B         Integral    Integral'
        write ( *, '(a)' ) ' '

        a = x1 - 1.0D+00

        do j = - 3, 12

          b = ( dble ( 10 - j ) * x1   
     &        + dble (      j ) * x2 ) 
     &        / dble ( 10     )

          call cubic_integrate ( a, b, q_exact )

          call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, 
     &      q_computed )

          write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) 
     &      j, a, b, q_exact, q_computed

        end do

      end do

      return
      end
      subroutine hermite_cubic_test04 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST04 tests HERMITE_CUBIC_SPLINE_VALUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 51 )
      integer nn
      parameter ( nn = 11 )

      double precision d(n)
      double precision dn(nn)
      double precision f(n)
      double precision fn(nn)
      integer i
      double precision s(n)
      double precision t(n)
      double precision u
      double precision v
      double precision x1
      double precision x2
      double precision x(n)
      double precision xn(nn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST04:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.'

      x1 = 0.0D+00
      x2 = 10.0D+00

      call r8vec_even ( nn, x1, x2, xn )

      fn(1:nn) = sin ( xn(1:nn) )
      dn(1:nn) = cos ( xn(1:nn) )

      call r8vec_even ( n, x1, x2, x )

      call hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, f, d, 
     &  s, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I      X       F computed     F exact      Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        u = sin ( x(i) )
        v = abs ( f(i) - u )
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) 
     &    i, x(i), f(i), u, v
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I      X       D computed     D exact      Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        u = cos ( x(i) )
        v = abs ( d(i) - u )
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) 
     &    i, x(i), d(i), u, v
      end do

      return
      end
      subroutine hermite_cubic_test05 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST05 tests HERMITE_CUBIC_TO_POWER_CUBIC
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1 )

      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision d(n)
      double precision d1
      double precision d1r
      double precision d2
      double precision d2r
      double precision f(n)
      double precision f1
      double precision f1r
      double precision f2
      double precision f2r
      double precision fp
      integer j
      double precision s(n)
      double precision s1
      double precision s2
      double precision t(n)
      double precision t1
      double precision t2
      double precision x(n)
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST05:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_TO_POWER_CUBIC converts the Hermite data'
      write ( *, '(a)' ) 
     &  '  to the coefficients of the power form of the polynomial'
      write ( *, '(a)' ) 
     &  '  POWER_CUBIC_TO_HERMITE_CUBIC converts the power form'
      write ( *, '(a)' ) '  to Hermite form'

      x1 = -1.0D+00
      x2 = +1.0D+00

      call cubic_value ( x1, f1, d1, s1, t1 )
      call cubic_value ( x2, f2, d2, s2, t2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Hermite data:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &  '  X1, F1, D1:', x1, f1, d1
      write ( *, '(a,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &  '  X2, F2, D2:', x2, f2, d2

      call hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, d2, c0, 
     &  c1, c2, c3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Power form:'
      write ( *, '(a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)' ) 
     &  '    p(x) = ', c0, ' + ', c1, ' * x + ', c2, ' * x^2 + ',
     &  c3, ' * x^3'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X       F (Hermite)  F (power)'
      write ( *, '(a)' ) ' '

      do j = - 3, 12

        x(1) = ( dble ( 10 - j ) * x1   
     &         + dble (      j ) * x2 ) 
     &         / dble ( 10     )

        call hermite_cubic_value ( x1, f1, d1, x2, f2, d2, 1, x, f, 
     &    d, s, t )

        fp = c0 + x(1) * ( c1 + x(1) * ( c2 + x(1) * c3 ) )

        write ( *, '(2x,f10.4,2x,f10.4,2x,f10.4)' ) x, f, fp

      end do

      call power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, f1r, 
     &  d1r, f2r, d2r )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Use POWER_CUBIC_TO_HERMITE_CUBIC to recover the'
      write ( *, '(a)' ) '  original Hermite data:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         Original   Recovered'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  F1:  ', f1, f1r
      write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  D1:  ', d1, d1r
      write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  F2:  ', f2, f2r
      write ( *, '(a,2x,f10.4,2x,f10.4)' ) '  D2:  ', d2, d2r

      return
      end
      subroutine hermite_cubic_test06 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST06 tests HERMITE_CUBIC_INTEGRATE using vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision b_hi
      double precision b_lo
      double precision d1
      double precision d2
      double precision f1
      double precision f2
      integer i
      double precision q_computed
      double precision q_exact
      double precision s1
      double precision s2
      double precision t1
      double precision t2
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST06:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic'
      write ( *, '(a)' ) '  polynomial from A to B.'
      write ( *, '(a)' ) '  Use A, B vectors for the calculation.'

      x1 = 0.0D+00
      x2 = 10.0D+00

      call cubic_value ( x1, f1, d1, s1, t1 )
      call cubic_value ( x2, f2, d2, s2, t2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                                 Exact       Computed'
      write ( *, '(a)' ) 
     &  '    J      A           B         Integral    Integral'
      write ( *, '(a)' ) ' '

      do i = -3, 12

        a = x1 - 1.0D+00
        b = ( dble ( 10 - i ) * x1 
     &      + dble (      i ) * x2 ) 
     &      / dble ( 10     )

        call cubic_integrate ( a, b, q_exact )

        call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b, 
     &    q_computed )

        write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) 
     &    i, a, b, q_exact, q_computed

      end do

      return
      end
      subroutine hermite_cubic_test07 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST07 tests HERMITE_CUBIC_INTEGRAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision d1
      double precision d2
      double precision f1
      double precision f2
      double precision q_computed
      double precision q_exact
      double precision s1
      double precision s2
      double precision t1
      double precision t2
      integer x_interval
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST07:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_INTEGRAL integrates a Hermite cubic'
      write ( *, '(a)' ) 
     &  '  polynomial over the definition interval [X1,X2].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                            Exact       Computed'
      write ( *, '(a)' ) 
     &  '     X1          X2         Integral    Integral'
      write ( *, '(a)' ) ' '

      do x_interval = 1, 3

        if ( x_interval .eq. 1 ) then
          x1 = 0.0D+00
          x2 = 10.0D+00
        else if ( x_interval .eq. 2 ) then
          x1 = -1.0D+00
          x2 = +1.0D+00
        else if ( x_interval .eq. 3 ) then
          x1 = 0.5D+00
          x2 = 0.75D+00
        end if

        call cubic_value ( x1, f1, d1, s1, t1 )
        call cubic_value ( x2, f2, d2, s2, t2 )

        call cubic_integrate ( x1, x2, q_exact )

        call hermite_cubic_integral ( x1, f1, d1, x2, f2, d2, 
     &    q_computed )

        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) 
     &    x1, x2, q_exact, q_computed

      end do

      return
      end
      subroutine hermite_cubic_test08 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST08 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nn
      parameter ( nn = 11 )

      double precision a
      double precision b
      double precision dn(nn)
      double precision fn(nn)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision q_computed
      double precision q_exact
      integer test
      double precision xn(nn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST08:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite'
      write ( *, '(a)' ) 
     &  '  cubic spline over the definition interval [X1,XNN].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                            Exact       Computed'
      write ( *, '(a)' ) 
     &  '     X1          XNN        Integral    Integral'
      write ( *, '(a)' ) ' '

      do test = 1, 3

        if ( test .eq. 1 ) then

          a = 0.0D+00
          b = 1.0D+00

          call r8vec_even ( nn, a, b, xn )

          do i = 1, nn
            fn(i) = xn(i) * ( 4.0D+00 * xn(i) - 1.0D+00 ) 
     &        * ( xn(i) - 1.0D+00 )
            dn(i) = 1.0D+00 + xn(i) * ( - 10.0D+00 + xn(i) * 12.0D+00 )
          end do

          q_exact = 
     &      ( xn(nn) * xn(nn) * ( 0.5D+00 + xn(nn) 
     &      * ( - ( 5.0D+00 / 3.0D+00 ) + xn(nn) ) ) ) 
     &      - ( xn(1)  * xn(1)  * ( 0.5D+00 + xn(1)  
     &      * ( - ( 5.0D+00 / 3.0D+00 ) + xn(1)  ) ) )
c
c  Use variable spacing.
c
        else if ( test .eq. 2 ) then

          a = 0.0D+00
          b = 1.0D+00

          call r8vec_even ( nn, a, b, xn )

          do i = 1, nn
            xn(i) = sqrt ( xn(i) )
            fn(i) = xn(i) * ( 4.0D+00 * xn(i) - 1.0D+00 ) 
     &        * ( xn(i) - 1.0D+00 )
            dn(i) = 1.0D+00 + xn(i) * ( - 10.0D+00 + xn(i) * 12.0D+00 )
          end do

          q_exact = 
     &      ( xn(nn) * xn(nn) * ( 0.5D+00 + xn(nn)
     &      * ( - ( 5.0D+00 / 3.0D+00 ) + xn(nn) ) ) ) 
     &      - ( xn(1)  * xn(1)  * ( 0.5D+00 + xn(1)  
     &      * ( - ( 5.0D+00 / 3.0D+00 ) + xn(1)  ) ) )
c
c  Try a non-cubic.
c
        else if ( test .eq. 3 ) then

          a = 0.0D+00
          b = pi

          call r8vec_even ( nn, a, b, xn )

          do i = 1, nn
            fn(i) = sin ( xn(i) )
            dn(i) = cos ( xn(i) )
          end do

          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )

        end if

        call hermite_cubic_spline_integral ( nn, xn, fn, dn, 
     &    q_computed )

        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) 
     &    xn(1), xn(nn), q_exact, q_computed

      end do

      return
      end
      subroutine hermite_cubic_test09 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST09 tests HERMITE_CUBIC_SPLINE_INTEGRATE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nn
      parameter ( nn = 11 )
      integer n
      parameter ( n = 25 )

      double precision a(n)
      double precision b(n)
      double precision dn(n)
      double precision fn(n)
      integer i
      double precision q(n)
      double precision q_exact
      double precision sn(n)
      double precision tn(n)
      double precision x1
      double precision x2
      double precision xn(nn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST09:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite'
      write ( *, '(a)' ) '  cubic spline from A to B.'
c
c  Define the cubic spline.
c
      x1 = 0.0D+00
      x2 = 10.0D+00

      call r8vec_even ( nn, x1, x2, xn )

      do i = 1, nn
        call cubic_value ( xn(i), fn(i), dn(i), sn(i), tn(i) )
      end do

      a(1:n) = 2.5D+00
      call r8vec_even ( n, x1 - 1.0D+00, x2 + 1.0D+00, b )

      call hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, b, q )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                                 Exact       Computed'
      write ( *, '(a)' ) 
     &  '    I      A           B         Integral    Integral'
      write ( *, '(a)' ) ' '

      do i = 1, n

        call cubic_integrate ( a(i), b(i), q_exact )

        write ( *, '(2x,i3,2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6)' ) 
     &    i, a(i), b(i), q_exact, q(i)

      end do

      return
      end
      subroutine hermite_cubic_test10 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST10 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nn
      parameter ( nn = 11 )

      character * ( 100 ) comment
      double precision dn(nn)
      double precision fn(nn)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision q_computed
      double precision q_exact
      double precision r8_uniform_01
      integer seed
      integer test
      double precision xn(nn)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST10:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite'
      write ( *, '(a)' ) 
     &  '  cubic spline over the definition interval [X1,XNN].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  If the subintervals are equally spaced, the derivative'
      write ( *, '(a)' ) 
     &  '  information has no effect on the result, except for'
      write ( *, '(a)' ) 
     &  '  the first and last values, DN(1) and DN(NN).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                            Exact       Computed'
      write ( *, '(a)' ) 
     &  '     X1          XNN        Integral    Integral  Comment'
      write ( *, '(a)' ) ' '

      do test = 1, 5
c
c  Equal spacing.
c
        if ( test .eq. 1 ) then

          call r8vec_even ( nn, 0.0D+00, pi, xn )
          do i = 1, nn
            fn(i) = sin ( xn(i) )
            dn(i) = cos ( xn(i) )
          end do
          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
          comment = 'Equal spacing, correct DN'
c
c  Equal spacing, reset DN(2:NN-1) to random numbers.
c
        else if ( test .eq. 2 ) then

          call r8vec_even ( nn, 0.0D+00, pi, xn )
          fn(1:nn) = sin ( xn(1:nn) )
          dn(1:nn) = cos ( xn(1:nn) )
          do i = 2, nn - 1
            dn(i) = 1000.0D+00 * r8_uniform_01 ( seed )
          end do
          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
          comment = 'Equal spacing, DN(2:N-1) random'
c
c  Equal spacing, now reset all of DN to random numbers.
c
        else if ( test .eq. 3 ) then

          call r8vec_even ( nn, 0.0D+00, pi, xn )
          fn(1:nn) = sin ( xn(1:nn) )
          do i = 1, nn
            dn(i) = 1000.0D+00 * r8_uniform_01 ( seed )
          end do
          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
          comment = 'Equal spacing, DN(1:N) random'
c
c  Variable spacing, correct data.
c
        else if ( test .eq. 4 ) then

          call r8vec_even ( nn, 0.0D+00, pi**2, xn )
          do i = 1, nn
            xn(i) = sqrt ( xn(i) )
            fn(i) = sin ( xn(i) )
            dn(i) = cos ( xn(i) )
          end do
          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
          comment = 'Variable spacing, correct DN'
c
c  Variable spacing, change one entry in DN.
c
        else if ( test .eq. 5 ) then

          call r8vec_even ( nn, 0.0D+00, pi**2, xn )
          xn(1:nn) = sqrt ( xn(1:nn) )
          fn(1:nn) = sin ( xn(1:nn) )
          dn(1:nn) = cos ( xn(1:nn) )
          dn( ( nn + 1 ) / 2 ) = 1000.0D+00 * r8_uniform_01 ( seed )
          q_exact = - cos ( xn(nn) ) + cos ( xn(1) )
          comment = 'Variable spacing, a single internal DN randomized.'

        end if

        call hermite_cubic_spline_integral ( nn, xn, fn, dn, 
     &    q_computed )

        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6,2x,a)' ) 
     &    xn(1), xn(nn), q_exact, q_computed, trim ( comment )

      end do

      return
      end
      subroutine hermite_cubic_test11 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST11 tests HERMITE_CUBIC_LAGRANGE_VALUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision d(4,n)
      double precision f(4,n)
      integer j
      double precision s(4,n)
      double precision t(4,n)
      double precision x(n)
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST11:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_LAGRANGE_VALUE evaluates the four'
      write ( *, '(a)' ) 
     &  '  Lagrange basis functions associated with F1, D1,'
      write ( *, '(a)' ) '  F2 and D2 such that'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  P(X) = F1 * LF1(X) + D1 * LD1(X)'
      write ( *, '(a)' ) '       + F2 * LF2(X) + D2 * LD2(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The first, second and third derivatives of these four'
      write ( *, '(a)' ) '  Lagrange basis functions are also computed.'

      x1 = 1.0D+00
      x2 = 2.0D+00
      call r8vec_even ( n, 0.0D+00, 2.5D+00, x );

      call hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, s, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Lagrange basis functions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I        X           LF1         LD1         ' //
     &  'LF2         LD2'
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, 
     &    '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    j, x(j), f(1:4,j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The derivative of the Lagrange basis functions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I        X           LF1         LD1         ' //
     &  'LF2         LD2'
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, 
     &    '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    j, x(j), d(1:4,j)
      end do

      return
      end
      subroutine hermite_cubic_test12 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST12 tests HERMITE_CUBIC_LAGRANGE_INTEGRAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision q(4)
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST12:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_LAGRANGE_INTEGRAL returns the integrals'
      write ( *, '(a)' ) 
     &  '  of the four Lagrange basis functions associated'
      write ( *, '(a)' ) '  with F1, D1, F2 and D2 such that'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  P(X) = F1 * LF1(X) + D1 * LD1(X)'
      write ( *, '(a)' ) '       + F2 * LF2(X) + D2 * LD2(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Lagrange basis function integrals:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '        X1          X2          LF1         LD1' // 
     &  '         LF2         LD2'
      write ( *, '(a)' ) ' '

      x2 = 1.0D+00

      do i = -6, 2
        x1 = dble ( i )
        call hermite_cubic_lagrange_integral ( x1, x2, q )
        write ( *, '(6(2x,f10.4))' ) x1, x2, q(1:4)
      end do

      return
      end
      subroutine hermite_cubic_test13 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST13 tests HERMITE_CUBIC_LAGRANGE_INTEGRATE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision f1
      double precision f2
      integer j
      double precision p(4)
      double precision q(4)
      double precision x1
      double precision x2

      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  'HERMITE_CUBIC_TEST13:'
      write ( *, '(a)' )  
     &  '  HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates a Hermite cubic'
      write ( *, '(a)' )  '  Lagrange polynomial from A to B.'
      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  '  Compute each result TWICE:'
      write ( *, '(a)' )  
     &  '  First row computed using HERMITE_CUBIC_INTEGRATE.'
      write ( *, '(a)' )  
     &  '  Second row computed using HERMITE_CUBIC_LAGRANGE_INTEGRATE.'

      x1 = 0.0D+00
      x2 = 10.0D+00

      write ( *, '(a)' )  ' '
      write ( *, '(a)' )  
     &  '        A           B           LF1         LD1' // 
     &  '         LF2         LD2'
      write ( *, '(a)' )  ' '

      a = x1 - 1.0D+00

      do j = -3, 12

        b = ( dble ( 10 - j ) * x1   
     &      + dble (      j ) * x2 ) 
     &      / dble ( 10     )

        f1 = 1.0D+00
        d1 = 0.0D+00
        f2 = 0.0D+00
        d2 = 0.0D+00
        call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, 
     &    b, p(1) )

        f1 = 0.0D+00
        d1 = 1.0D+00
        f2 = 0.0D+00
        d2 = 0.0D+00
        call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, 
     &    b, p(2) )

        f1 = 0.0D+00
        d1 = 0.0D+00
        f2 = 1.0D+00
        d2 = 0.0D+00
        call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, 
     &    b, p(3) )

        f1 = 0.0D+00
        d1 = 0.0D+00
        f2 = 0.0D+00
        d2 = 1.0D+00
        call hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, 
     &    b, p(4) )

        call hermite_cubic_lagrange_integrate ( x1, x2, a, b, q )

        write ( *, '(6(2x,f10.4))' ) a, b, p(1:4)
        write ( *, '(24x,4(2x,f10.4))' )       q(1:4)

      end do

      return
      end
      subroutine hermite_cubic_test14 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST14 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision dn(n)
      double precision fn(n)
      integer i
      integer j
      integer k
      double precision q
      double precision r(n)
      integer seed
      double precision w(2,n)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST14:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule'
      write ( *, '(a)' ) '  for Hermite cubic splines.'

      do k = 1, 2

        write ( *, '(a)' ) ' '
        if ( k .eq. 1 ) then
          write ( *, '(a)' ) '  Case 1: Random spacing'

          seed = 123456789

          call r8vec_uniform_01 ( n, seed, r )

          x(1) = r(1)
          do i = 2, n
            x(i) = x(i-1) + r(i)
          end do
        else if ( k .eq. 2 ) then
          write ( *, '(a)' ) '  Case 2: Uniform spacing'
          write ( *, '(a)' ) '  F(2:N-1) have equal weight.'
          write ( *, '(a)' ) '  D(2:N-1) have zero weight.'
          do i = 1, n
            x(i) = dble ( 10 + i - 1 ) / 20.0D+00
          end do
        end if

        call hermite_cubic_spline_quad_rule ( n, x, w )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '   I   J        X         W                Q'
        write ( *, '(a)' ) ' '

        do i = 1, 2

          do j = 1, n

            fn(1:n) = 0.0D+00
            dn(1:n) = 0.0D+00

            if ( i .eq. 1 ) then
              fn(j) = 1.0D+00
            else
              dn(j) = 1.0D+00
            end if

            call hermite_cubic_spline_integral ( n, x, fn, dn, q )

            write ( *, '(2x,i2,2x,i2,2x,f10.4, 2x, g14.6,2x,g14.6)' ) 
     &        i, j, x(j), w(i,j), q

          end do

        end do

      end do

      return
      end
      subroutine hermite_cubic_test15 ( )

c*********************************************************************72
c
cc HERMITE_CUBIC_TEST15 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
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
      implicit none

      integer n
      parameter ( n = 11 )

      double precision dn(n)
      double precision fn(n)
      integer i
      integer j
      double precision q
      double precision q_exact
      double precision r(n)
      double precision s
      integer seed
      double precision t
      double precision w(2,n)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_CUBIC_TEST15:'
      write ( *, '(a)' ) 
     &  '  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule'
      write ( *, '(a)' ) '  for Hermite cubic splines.'

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, r )

      x(1) = r(1)
      do j = 2, n
        x(j) = x(j-1) + r(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Random spacing'
      write ( *, '(a,i8)' ) '  Number of points N = ', n
      write ( *, '(a,g14.6,a,g14.6,a)' ) 
     &  '  Interval = [', x(1), ',', x(n), ']'

      call hermite_cubic_spline_quad_rule ( n, x, w )

      do j = 1, n
        call cubic_value ( x(j), fn(j), dn(j), s, t )
      end do

      q = 0.0D+00
      do j = 1, n
        q = q + w(1,j) * fn(j) + w(2,j) * dn(j)
      end do

      call cubic_integrate ( x(1), x(n), q_exact )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Q         = ', q
      write ( *, '(a,g14.6)' ) '  Q (exact) = ', q_exact

      return
      end
      function cubic_antiderivative ( x )

c*********************************************************************72
c
cc CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 January 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision CUBIC_ANTIDERIVATIVE, the value.
c
      implicit none

      double precision cubic_antiderivative
      double precision x

      cubic_antiderivative = x * x * ( 5.0D+00 
     &  + x * ( - 7.0D+00 / 3.0D+00 + x * 1.0D+00 / 4.0D+00 ) )

      return
      end
      subroutine cubic_integrate ( a, b, q )

c*********************************************************************72
c
cc CUBIC_INTEGRATE integrates the cubic from A to B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the integration interval.
c
c    Output, double precision Q, the integral from A to B.
c
      implicit none

      double precision a
      double precision b
      double precision cubic_antiderivative
      double precision q

      q = cubic_antiderivative ( b ) - cubic_antiderivative ( a )

      return
      end
      subroutine cubic_value ( x, f, d, s, t )

c*********************************************************************72
c
cc CUBIC_VALUE evaluates a cubic function.
c
c  Discussion:
c
c    f(x) =   x^3 -  7 x^2 + 10 x
c    d(x) = 3 x^2 - 14 x   + 10
c    s(x) = 6 x   - 14
c    t(x) = 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision F, D, S, T, the value and first three
c    derivatives of the cubic function.
c
      implicit none

      double precision d
      double precision f
      double precision s
      double precision t
      double precision x

      f = 0.0D+00 + x * ( 10.0D+00 + x * (  - 7.0D+00 + x * 1.0D+00 ) )
      d =                 10.0D+00 + x * ( - 14.0D+00 + x * 3.0D+00 )
      s =                                  - 14.0D+00 + x * 6.0D+00
      t =                                                   6.0D+00

      return
      end
