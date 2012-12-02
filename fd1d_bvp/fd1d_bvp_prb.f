      program main

c*********************************************************************72
c
cc MAIN is the main program for FD1D_BVP_PRB.
c
c  Discussion:
c
c    FD1D_BVP_PRB tests the routines in FD1D_BVP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FD1D_BVP library.'

      call fd1d_bvp_test01 ( )
      call fd1d_bvp_test02 ( )
      call fd1d_bvp_test03 ( )
      call fd1d_bvp_test04 ( )
      call fd1d_bvp_test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine fd1d_bvp_test01 ( )

c*********************************************************************72
c
cc FD1D_BVP_TEST01 carries out test case #1.
c
c  Discussion:
c
c    Use A1, C1, F1, EXACT1.
c
c    Repeat for a nonuniform mesh.
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

      integer n
      parameter ( n = 21 )

      double precision a1
      external a1
      double precision a1prime
      external a1prime
      double precision c1
      external c1
      double precision f1
      external f1
      character * ( 80 ) filename
      integer i
      double precision u(n)
      double precision u2(2,n)
      double precision uexact(n)
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 0.0D+00
      x2 = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST01'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  A1''(X) = 0.0'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  X1 = ', x1
      write ( *, '(a,g14.6)' ) '  X2 = ', x2

      call r8vec_even ( n, x1, x2, x )

      call fd1d_bvp ( n, a1, a1prime, c1, f1, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat for a nonuniform mesh.'

      call r8vec_even ( n, x1, x2, x )

      do i = 1, n
        x(i) = sqrt ( x(i) )
      end do

      call fd1d_bvp ( n, a1, a1prime, c1, f1, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do
c
c  Write the data to files.
c
      filename = 'fd1d_bvp_test01_nodes.txt'
      call r8mat_write ( filename, 1, n, x );

      do i = 1, n
        u2(1,i) = u(i)
        u2(2,i) = uexact(i)
      end do

      filename = 'fd1d_bvp_test01_values.txt'
      call r8mat_write ( filename, 2, n, u2 )

      return
      end
      subroutine fd1d_bvp_test02 ( )

c*********************************************************************72
c
cc FD1D_BVP_TEST02 carries out test case #2.
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
c    15 February 2011
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
      double precision a1prime
      external a1prime
      double precision c2
      external c2
      double precision f2
      external f2
      character * ( 80 ) filename
      integer i
      double precision u(n)
      double precision u2(2,n)
      double precision uexact(n)
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 0.0D+00
      x2 = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST02'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  A1''(X) = 0.0'
      write ( *, '(a)' ) '  C2(X)  = 2.0'
      write ( *, '(a)' ) '  F2(X)  = X * ( 5 - X ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  X1 = ', x1
      write ( *, '(a,g14.6)' ) '  X2 = ', x2

      call r8vec_even ( n, x1, x2, x )

      call fd1d_bvp ( n, a1, a1prime, c2, f2, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do
c
c  Write the data to files.
c
      filename = 'fd1d_bvp_test02_nodes.txt'
      call r8mat_write ( filename, 1, n, x );

      do i = 1, n
        u2(1,i) = u(i)
        u2(2,i) = uexact(i)
      end do

      filename = 'fd1d_bvp_test02_values.txt'
      call r8mat_write ( filename, 2, n, u2 )

      return
      end
      subroutine fd1d_bvp_test03 ( )

c*********************************************************************72
c
cc FD1D_BVP_TEST03 carries out test case #3.
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
c    15 February 2011
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
      double precision a1prime
      external a1prime
      double precision c3
      external c3
      double precision f3
      external f3
      character * ( 80 ) filename
      integer i
      double precision u(n)
      double precision u2(2,n)
      double precision uexact(n)
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 0.0D+00
      x2 = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST03'
      write ( *, '(a)' ) '  A1(X)  = 1.0'
      write ( *, '(a)' ) '  A1''(X) = 0.0'
      write ( *, '(a)' ) '  C3(X)  = 2.0 * X'
      write ( *, '(a)' )
     &  '  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  X1 = ', x1
      write ( *, '(a,g14.6)' ) '  X2 = ', x2

      call r8vec_even ( n, x1, x2, x )

      call fd1d_bvp ( n, a1, a1prime, c3, f3, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do
c
c  Write the data to files.
c
      filename = 'fd1d_bvp_test03_nodes.txt'
      call r8mat_write ( filename, 1, n, x );

      do i = 1, n
        u2(1,i) = u(i)
        u2(2,i) = uexact(i)
      end do

      filename = 'fd1d_bvp_test03_values.txt'
      call r8mat_write ( filename, 2, n, u2 )

      return
      end
      subroutine fd1d_bvp_test04 ( )

c*********************************************************************72
c
cc FD1D_BVP_TEST04 carries out test case #4.
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
c    15 February 2011
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
      double precision a2prime
      external a2prime
      double precision c1
      external c1
      double precision f4
      external f4
      character * ( 80 ) filename
      integer i
      double precision u(n)
      double precision u2(2,n)
      double precision uexact(n)
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 0.0D+00
      x2 = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST04'
      write ( *, '(a)' ) '  A2(X)  = 1.0 + X * X'
      write ( *, '(a)' ) '  A2''(X) = 2.0 * X'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' )
     &  '  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  X1 = ', x1
      write ( *, '(a,g14.6)' ) '  X2 = ', x2

      call r8vec_even ( n, x1, x2, x )

      call fd1d_bvp ( n, a2, a2prime, c1, f4, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do
c
c  Write the data to files.
c
      filename = 'fd1d_bvp_test04_nodes.txt'
      call r8mat_write ( filename, 1, n, x );

      do i = 1, n
        u2(1,i) = u(i)
        u2(2,i) = uexact(i)
      end do

      filename = 'fd1d_bvp_test04_values.txt'
      call r8mat_write ( filename, 2, n, u2 )

      return
      end
      subroutine fd1d_bvp_test05 ( )

c*********************************************************************72
c
cc FD1D_BVP_TEST05 carries out test case #5.
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
c    15 February 2011
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
      double precision a3prime
      external a3prime
      double precision c1
      external c1
      double precision f5
      external f5
      character * ( 80 ) filename
      integer i
      double precision u(n)
      double precision u2(2,n)
      double precision uexact(n)
      double precision x(n)
      double precision x1
      double precision x2

      x1 = 0.0D+00
      x2 = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BVP_TEST05'
      write ( *, '(a)' ) '  A3(X)  = 1.0 + X * X for X <= 1/3'
      write ( *, '(a)' ) '         = 7/9 + X     for      1/3 < X'
      write ( *, '(a)' ) '  A3''(X) = 2.0 * X     for X <= 1/3'
      write ( *, '(a)' ) '           1           for      1/3 < X'
      write ( *, '(a)' ) '  C1(X)  = 0.0'
      write ( *, '(a)' )
     &  '  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
      write ( *, '(a)' ) '                       for X <= 1/3'
      write ( *, '(a)' )
     &  '         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) * exp ( X )'
      write ( *, '(a)' ) '                       for      1/3 <= X'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  X1 = ', x1
      write ( *, '(a,g14.6)' ) '  X2 = ', x2

      call r8vec_even ( n, x1, x2, x )

      call fd1d_bvp ( n, a3, a3prime, c1, f5, x, u )

      call exact1 ( n, x, uexact )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' )
     &    i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
      end do
c
c  Write the data to files.
c
      filename = 'fd1d_bvp_test05_nodes.txt'
      call r8mat_write ( filename, 1, n, x );

      do i = 1, n
        u2(1,i) = u(i)
        u2(2,i) = uexact(i)
      end do

      filename = 'fd1d_bvp_test05_values.txt'
      call r8mat_write ( filename, 2, n, u2 )

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
c    30 May 2009
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
      function a1prime ( x )

c*********************************************************************72
c
cc A1PRIME evaluates A' function #1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A1PRIME, the value of A'(X).
c
      implicit none

      double precision a1prime
      double precision x

      a1prime = 0.0D+00

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
c    30 May 2009
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
      function a2prime ( x )

c*********************************************************************72
c
cc A2PRIME evaluates A' function #2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A2PRIME, the value of A'(X).
c
      implicit none

      double precision a2prime
      double precision x

      a2prime = 2.0D+00 * x

      return
      end
      function a3 ( x )

c*****************************************************************************80
c
cc A3 evaluates A function #3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
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

      if ( x <= 1.0D+00 / 3.0D+00 ) then
        a3 = 1.0D+00 + x * x
      else
        a3 = x + 7.0D+00 / 9.0D+00
      end if

      return
      end
      function a3prime ( x )

c*********************************************************************72
c
cc A3PRIME evaluates A' function #3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision A3PRIME, the value of A'(X).
c
      implicit none

      double precision a3prime
      double precision x

      if ( x <= 1.0D+00 / 3.0D+00 ) then
        a3prime = 2.0D+00 * x
      else
        a3prime = 1.0D+00
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
c    30 May 2009
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
c    30 May 2009
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
c    30 May 2009
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
c    30 May 2009
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
c    30 May 2009
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
c    30 May 2009
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
c    30 May 2009
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

      f4 = ( x + 3.0D+00 * x**2 + 5.0D+00 * x**3 + x**4 ) * exp ( x )

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
c    30 May 2009
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

      if ( x <= 1.0D+00 / 3.0D+00 ) then
        f5 = ( x + 3.0D+00 * x**2 + 5.0D+00 * x**3 + x**4 ) * exp ( x )
      else
        f5 = ( - 1.0D+00 + ( 10.0D+00 / 3.0D+00 ) * x
     &    + ( 43.0D+00 / 9.0D+00 ) * x**2 + x**3 ) * exp ( x )
      end if

      return
      end
      subroutine exact1 ( n, x, uexact )

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
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision UEXACT(N), the values of U(X(1:N)).
c
      implicit none

      integer n

      integer i
      double precision uexact(n)
      double precision x(n)

      do i = 1, n
        uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
      end do

      return
      end
      subroutine exact2 ( n, x, uexact )

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
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision UEXACT(N), the values of U(X(1:N)).
c
      implicit none

      integer n

      integer i
      double precision uexact(n)
      double precision x(n)

      do i = 1, n
        if ( x(i) <= 2.0D+00 / 3.0D+00 ) then
          uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
        else
          uexact(i) = x(i) * ( 1.0D+00 - x(i) )
     &      * exp ( 2.0D+00 / 3.0D+00 )
        end if
      end do

      return
      end
      subroutine exact3 ( n, x, uexact )

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
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision UEXACT(N), the values of U(X(1:N)).
c
      implicit none

      integer n

      integer i
      double precision uexact(n)
      double precision x(n)

      do i = 1, n
        if ( x(i) <= 2.0D+00 / 3.0D+00 ) then
          uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
        else
          uexact(i) = x(i) * ( 1.0D+00 - x(i) )
        end if
      end do

      return
      end
