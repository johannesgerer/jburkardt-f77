      program main

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_PRB tests the routines in FEM1D_BVP_LINEAR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FEM1D_BVP_LINEAR library.'

      call fem1d_bvp_linear_test01 ( )
      call fem1d_bvp_linear_test02 ( )
      call fem1d_bvp_linear_test03 ( )
      call fem1d_bvp_linear_test04 ( )
      call fem1d_bvp_linear_test05 ( )
      call fem1d_bvp_linear_test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine fem1d_bvp_linear_test01 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST01 carries out test case #1.
c
c  Discussion:
c
c    Use A1, C1, F1, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a1
      external a1
      double precision c1
      external c1
      double precision exact1
      external exact1
      double precision exact_ux1
      external exact_ux1
      double precision f1
      external f1
      integer i
      double precision l2_error
      double precision seminorm_error
      double precision u(n)
      double precision uexact
      double precision x(n)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST01'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
c
c  Geometry definitions.
c
      x_first = 0.0D+00
      x_last = 1.0D+00
      call r8vec_even ( n, x_first, x_last, x )

      call fem1d_bvp_linear ( n, a1, c1, f1, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      call compute_l2_error ( n, x, u, exact1, l2_error )
      call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
      write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

      return
      end
      subroutine fem1d_bvp_linear_test02 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST02 carries out test case #2.
c
c  Discussion:
c
c    Use A1, C2, F2, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a1
      external a1
      double precision c2
      external c2
      double precision exact1
      external exact1
      double precision exact_ux1
      external exact_ux1
      double precision f2
      external f2
      integer i
      double precision l2_error
      double precision seminorm_error
      double precision u(n)
      double precision uexact
      double precision x(n)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST02'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  C2(X)  = 2.0'
      write ( *, '(a)' ) '  F2(X)  = X * ( 5 - X ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
c
c  Geometry definitions.
c
      x_first = 0.0D+00
      x_last = 1.0D+00
      call r8vec_even ( n, x_first, x_last, x )

      call fem1d_bvp_linear ( n, a1, c2, f2, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      call compute_l2_error ( n, x, u, exact1, l2_error )
      call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
      write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

      return
      end
      subroutine fem1d_bvp_linear_test03 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST03 carries out test case #3.
c
c  Discussion:
c
c    Use A1, C3, F3, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a1
      external a1
      double precision c3
      external c3
      double precision exact1
      external exact1
      double precision exact_ux1
      external exact_ux1
      double precision f3
      external f3
      integer i
      double precision l2_error
      double precision seminorm_error
      double precision u(n)
      double precision uexact
      double precision x(n)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST03'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  C3(X)  = 2.0 * X'
      write ( *, '(a)' ) 
     &  '  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
c
c  Geometry definitions.
c
      x_first = 0.0D+00
      x_last = 1.0D+00
      call r8vec_even ( n, x_first, x_last, x )

      call fem1d_bvp_linear ( n, a1, c3, f3, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      call compute_l2_error ( n, x, u, exact1, l2_error )
      call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
      write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

      return
      end
      subroutine fem1d_bvp_linear_test04 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST04 carries out test case #4.
c
c  Discussion:
c
c    Use A2, C1, F4, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a2
      external a2
      double precision c1
      external c1
      double precision exact1
      external exact1
      double precision exact_ux1
      external exact_ux1
      double precision f4
      external f4
      integer i
      double precision l2_error
      double precision seminorm_error
      double precision u(n)
      double precision uexact
      double precision x(n)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST04'
      write ( *, '(a)' ) '  A2(X)  = 1.0 + X * X'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' ) 
     &  '  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
c
c  Geometry definitions.
c
      x_first = 0.0D+00
      x_last = 1.0D+00
      call r8vec_even ( n, x_first, x_last, x )

      call fem1d_bvp_linear ( n, a2, c1, f4, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      call compute_l2_error ( n, x, u, exact1, l2_error )
      call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
      write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

      return
      end
      subroutine fem1d_bvp_linear_test05 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST05 carries out test case #5.
c
c  Discussion:
c
c    Use A3, C1, F5, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a3
      external a3
      double precision c1
      external c1
      double precision exact1
      external exact1
      double precision exact_ux1
      external exact_ux1
      double precision f5
      external f5
      integer i
      double precision l2_error
      double precision seminorm_error
      double precision u(n)
      double precision uexact
      double precision x(n)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST05'
      write ( *, '(a)' ) '  A3(X)  = 1.0 + X * X for X .le. 1/3'
      write ( *, '(a)' ) '         = 7/9 + X     for      1/3 < X'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' ) 
     &  '  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
      write ( *, '(a)' ) '                       for X .le. 1/3'
      write ( *, '(a)' ) 
     &  '         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )'
      write ( *, '(a)' ) '                       for      1/3 .le. X'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
c
c  Geometry definitions.
c
      x_first = 0.0D+00
      x_last = 1.0D+00
      call r8vec_even ( n, x_first, x_last, x )

      call fem1d_bvp_linear ( n, a3, c1, f5, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      call compute_l2_error ( n, x, u, exact1, l2_error )
      call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
      write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

      return
      end
      subroutine fem1d_bvp_linear_test06 ( )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR_TEST06 does an error analysis.
c
c  Discussion:
c
c    Use A1, C1, F6, EXACT4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a1
      external a1
      double precision c1
      external c1
      double precision exact4
      external exact4
      double precision exact_ux4
      external exact_ux4
      double precision f6
      external f6
      integer i
      double precision l2_error
      integer n
      double precision seminorm_error
      double precision u(161)
      double precision uexact
      double precision x(161)
      double precision x_first
      double precision x_last

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST06'
      write ( *, '(a)' ) '  A1(X)  = 1.0 '
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' ) '  F6(X)  = pi*pi*sin(pi*X)'
      write ( *, '(a)' ) '  U4(X)  = sin(pi*x)'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Compute L2 norm and seminorm of error for various N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N        L2 error      Seminorm error'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '

      n = 11

      do i = 0, 4
!
!  Geometry definitions.
!
        x_first = 0.0D+00
        x_last = 1.0D+00

        call r8vec_even ( n, x_first, x_last, x )

        call fem1d_bvp_linear ( n, a1, c1, f6, x, u )

        call compute_l2_error ( n, x, u, exact4, l2_error )
        call compute_seminorm_error ( n, x, u, exact_ux4, 
     &    seminorm_error )

        write ( *, '(2x,i4,2x,f14.6,2x,f14.6)' ) 
     &    n, l2_error, seminorm_error

        n = 2 * ( n - 1 ) + 1

      end do

      return
      end
      function a1 ( x )

c*********************************************************************72
c
cc A1 evaluates A function #1.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A1, the value of A(X).
c
      implicit none

      double precision a1
      double precision x

      a1 = 1.0D+00

      return
      end
      function a2 ( x )

c*********************************************************************72
c
cc A2 evaluates A function #2.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A2, the value of A(X).
c
      implicit none

      double precision a2
      double precision x

      a2 = 1.0D+00 + x * x

      return
      end
      function a3 ( x )

c*********************************************************************72
c
cc A3 evaluates A function #3.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A3, the value of A(X).
c
      implicit none

      double precision a3
      double precision x

      if ( x .le. 1.0D+00 / 3.0D+00 ) then
        a3 = 1.0D+00 + x * x
      else
        a3 = x + 7.0D+00 / 9.0D+00
      end if

      return
      end
      function c1 ( x )

c*********************************************************************72
c
cc C1 evaluates C function #1.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision C1, the value of C(X).
c
      implicit none

      double precision c1
      double precision x

      c1 = 0.0D+00

      return
      end
      function c2 ( x )

c*********************************************************************72
c
cc C2 evaluates C function #2.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision C2, the value of C(X).
c
      implicit none

      double precision c2
      double precision x

      c2 = 2.0D+00

      return
      end
      function c3 ( x )

c*********************************************************************72
c
cc C3 evaluates C function #3.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision C3, the value of C(X).
c
      implicit none

      double precision c3
      double precision x

      c3 = 2.0D+00 * x

      return
      end
      function exact1 ( x )

c*********************************************************************72
c
cc EXACT1 evaluates exact solution #1.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT1, the value of U(X).
c
      implicit none

      double precision exact1
      double precision x

      exact1 = x * ( 1.0D+00 - x ) * exp ( x )

      return
      end
      function exact_ux1 ( x )

c*********************************************************************72
c
cc EXACT_UX1 evaluates the derivative of exact solution #1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT_UX1, the value of dUdX(X).
c
      implicit none

      double precision exact_ux1
      double precision x

      exact_ux1 = ( 1.0D+00 - x - x * x ) * exp ( x )

      return
      end
      function exact2 ( x )

c*********************************************************************72
c
cc EXACT2 returns exact solution #2.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT2, the value of U(X).
c
      implicit none

      double precision exact2
      double precision x

      if ( x .le. 2.0D+00 / 3.0D+00 ) then
        exact2 =  x * ( 1.0D+00 - x ) * exp ( x )
      else
        exact2 = x * ( 1.0D+00 - x )  * exp ( 2.0D+00 / 3.0D+00 )
      end if

      return
      end
      function exact3 ( x )

c*********************************************************************72
c
cc EXACT3 returns exact solution #3.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT3, the value of U(X).
c
      implicit none

      double precision exact3
      double precision x

      if ( x .le. 2.0D+00 / 3.0D+00 ) then
        exact3 = x * ( 1.0D+00 - x ) * exp ( x )
      else
        exact3 = x * ( 1.0D+00 - x )
      end if

      return
      end
      function exact4 ( x )

c*********************************************************************72
c
cc EXACT4 evaluates exact solution #4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT4, the value of U(X).
c
      implicit none

      double precision exact4
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      exact4 = sin ( pi * x )

      return
      end
      function exact_ux4 ( x )

c*********************************************************************72
c
cc EXACT_UX4 evaluates the derivative of exact solution #4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT_UX4, the value of dUdX(X).
c
      implicit none

      double precision exact_ux4
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      exact_ux4 = pi * cos ( pi * x )

      return
      end
      function f1 ( x )

c*********************************************************************72
c
cc F1 evaluates right hand side function #1.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F1, the value of F(X).
c
      implicit none

      double precision f1
      double precision x

      f1 = x * ( x + 3.0D+00 ) * exp ( x )

      return
      end
      function f2 ( x )

c*********************************************************************72
c
cc F2 evaluates right hand side function #2.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F2, the value of F(X).
c
      implicit none

      double precision f2
      double precision x

      f2 = x * ( 5.0D+00 - x ) * exp ( x )

      return
      end
      function f3 ( x )

c*********************************************************************72
c
cc F3 evaluates right hand side function #3.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F3, the value of F(X).
c
      implicit none

      double precision f3
      double precision x

      f3 = - x * ( 2.0D+00 * x * x - 3.0D+00 * x - 3.0D+00 ) * exp ( x )

      return
      end
      function f4 ( x )

c*********************************************************************72
c
cc F4 evaluates right hand side function #4.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F4, the value of F(X).
c
      implicit none

      double precision f4
      double precision x

      f4 = ( x + 3.0D+00 * x * x 
     &  + 5.0D+00 * x * x * x + x * x * x * x ) * exp ( x )

      return
      end
      function f5 ( x )

c*********************************************************************72
c
cc F5 evaluates right hand side function #5.
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
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F5, the value of F(X).
c
      implicit none

      double precision f5
      double precision x

      if ( x .le. 1.0D+00 / 3.0D+00 ) then
        f5 = ( x + 3.0D+00 * x * x + 5.0D+00 * x**3 + x**4 ) * exp ( x )
      else
        f5 = ( - 1.0D+00 + ( 10.0D+00 / 3.0D+00 ) * x 
     &    + ( 43.0D+00 / 9.0D+00 ) * x * x + x * x * x ) * exp ( x )
      end if

      return
      end
      function f6 ( x )

c*********************************************************************72
c
cc F6 evaluates right hand side function #6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision F6, the value of F(X).
c
      implicit none

      double precision f6
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      f6 = pi * pi * sin ( pi * x )

      return
      end
