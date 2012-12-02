      program main

c*********************************************************************72
c
cc MAIN is the main program for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
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
      double precision f1
      external f1
      double precision k1
      external k1
      double precision u(n)
      character * ( 80 ) u_file
      double precision ua
      double precision ub
      double precision x(n)
      character * ( 80 ) x_file

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROBLEM1:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'

      a = 0.0D+00
      b = 1.0D+00

      ua = 0.0D+00
      ub = 1.0D+00

      call fd1d_heat_steady ( n, a, b, ua, ub, k1, f1, x, u )

      x_file = 'problem1_nodes.txt'
      call r8mat_write ( x_file, 1, n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  X data written to "' // trim ( x_file ) // '".'

      u_file = 'problem1_values.txt'
      call r8mat_write ( u_file, 1, n, u )

      write ( *, '(a)' )
     &  '  U data written to "' // trim ( u_file ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROBLEM1:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function k1 ( x )

c*********************************************************************72
c
cc K1 evaluates the heat transfer coefficient K(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the position.
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
cc F1 evaluates the right hand side of the steady state heat equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the position.
c
c    Output, double precision F1, the value of F(X).
c
      implicit none

      double precision f1
      double precision x

      f1 = 0.0D+00

      return
      end


