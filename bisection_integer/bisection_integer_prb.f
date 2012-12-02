      program main

c*********************************************************************72
c
cc MAIN is the main program for BISECTION_INTEGER_PRB.
c
c  Discussion:
c
c    BISECTION_INTEGER_PRB tests the BISECTION_INTEGER library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BISECTION_INTEGER_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BISECTION_INTEGER library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BISECTION_INTEGER_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BISECTION_INTEGER;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer, external :: f01
      integer fc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  BISECTION_INTEGER attempts to locate an integer root C'
      write ( *, '(a)' ) '  of an equation F(C) = 0.'
      write ( *, '(a)' ) 
     &  '  The user supplies a change of sign interval [A,B].'
      write ( *, '(a)' ) 
     &  '  The function considered here has two real roots'
      write ( *, '(a)' ) 
     &  '  as well as an integer root, so the algorithm can'
      write ( *, '(a)' ) 
     &  '  fail depending on how the change of sign interval is chosen.'

      a = 4
      b = 100

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The initial change of sign interval is:'
      write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
      write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )

      call bisection_integer ( f01, a, b, c, fc )

      if ( fc .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  An exact root was found at C = ', c
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  An exact root was NOT found.'
        write ( *, '(a,i8)' ) '  The change of sign interval is now:'
        write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
        write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )
      end if

      a = -10
      b = 15

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The initial change of sign interval is:'
      write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
      write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )

      call bisection_integer ( f01, a, b, c, fc )

      if ( fc .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  An exact root was found at C = ', c
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  An exact root was NOT found.'
        write ( *, '(a,i8)' ) '  The change of sign interval is now:'
        write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', f01 ( a )
        write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', f01 ( b )
      end if


      return
      end
      function f01 ( n )

c*********************************************************************72
c
cc F01 is a test function.
c
c  Discussion:
c
c    The polynomial has roots 1/2, 7/2, and 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument.
c
c    Output, integer F01, the function value.
c
      implicit none

      integer f01
      integer n

      f01 = ( 2 * n - 7 ) * ( 2 * n - 1 ) * ( n - 10 )

      return
      end
