      program main

c*********************************************************************72
c
cc FEM1D_HEAT_STEADY_PRB tests the routines in FEM1D_HEAT_STEADY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FEM1D_HEAT_STEADY library.'

      call fem1d_heat_steady_test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine fem1d_heat_steady_test01 ( )

c*********************************************************************72
c
cc FEM1D_HEAT_STEADY_TEST01 carries out test case #1.
c
c  Discussion:
c
c    Use K1, F1, EXACT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision a
      double precision b
      double precision exact1
      double precision, external :: f1
      integer i
      double precision, external :: k1
      double precision u(n)
      double precision ua
      double precision ub
      double precision uexact
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_HEAT_STEADY_TEST01'
      write ( *, '(a)' ) '  K1(X)  = 1.0'
      write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
      write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
c
c  Geometry definitions.
c
      a = 0.0D+00
      b = 1.0D+00
      ua = 0.0D+00
      ub = 0.0D+00
      call r8vec_even ( n, a, b, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,f10.4)' ) '  Left endpoint A = ', a
      write ( *, '(a,f10.4)' ) '  Right endpoint B = ', b
      write ( *, '(a,f10.4)' ) '  Prescribed U(A) = ', ua
      write ( *, '(a,f10.4)' ) '  Prescribed U(B) = ', ub

      call fem1d_heat_steady ( n, a, b, ua, ub, k1, f1, x, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I    X         U         Uexact    Error'
      write ( *, '(a)' ) ' '

      do i = 1, n
        uexact = exact1 ( x(i) )
        write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), u(i), uexact, abs ( u(i) - uexact )
      end do

      return
      end
      function k1 ( x )

c*********************************************************************72
c
cc K1 evaluates K function #1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision K1, the value of K(X).
c
      implicit none

      double precision k1
      double precision x

      k1 = 1.0D+00

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
c    08 April 2011
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
c    08 April 2011
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

