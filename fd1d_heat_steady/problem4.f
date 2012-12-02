      program main

c*********************************************************************72
c
cc MAIN is the main program for problem 4.
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
      parameter ( n = 21 )

      double precision a
      double precision b
      double precision f4
      external f4
      double precision k4
      external k4
      double precision u(n)
      character * ( 80 ) u_file
      double precision ua
      double precision ub
      double precision x(n)
      character * ( 80 ) x_file

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROBLEM4:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'
      write ( *, '(a)' ) '  A heat source and a heat sink.'

      a = 0.0D+00
      b = 1.0D+00

      ua = 0.0D+00
      ub = 0.0D+00

      call fd1d_heat_steady ( n, a, b, ua, ub, k4, f4, x, u )

      x_file = 'problem4_nodes.txt'
      call r8mat_write ( x_file, 1, n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  X data written to "' // trim ( x_file ) // '".'

      u_file = 'problem4_values.txt'
      call r8mat_write ( u_file, 1, n, u )

      write ( *, '(a)' )
     &  '  U data written to "' // trim ( u_file ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROBLEM4:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function k4 ( x )

c*********************************************************************72
c
cc K4 evaluates the heat transfer coefficient K(X).
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
c    Output, double precision K4, the value of K(X).
c
      implicit none

      double precision k4
      double precision x

      k4 = 1.0D+00

      return
      end
      function f4 ( x )

c*********************************************************************72
c
cc F4 evaluates the right hand side of the steady state heat equation.
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
c    Output, double precision F4, the value of F(X).
c
      implicit none

      double precision f4
      double precision x

      if ( x < 0.15D+00 ) then
        f4 = 0.0D+00
      else if ( x < 0.35D+00 ) then
        f4 = 1.0D+00
      else if ( x < 0.75D+00 ) then
        f4 = 0.0D+00
      else if ( x < 0.85D+00 ) then
        f4 = - 2.0D+00
      else
        f4 = 0.0D+00
      end if

      return
      end


