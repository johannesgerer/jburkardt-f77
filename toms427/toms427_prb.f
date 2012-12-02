      program main

c*********************************************************************72
c
cc TOMS427_PRB tests FRCOS.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS427_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS427 library.'
      write ( *, '(a)' ) ' '

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS427_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 handles F(X) = 1/(X*X+1).
c
      implicit none

      real et
      real exact
      external f1
      real frcos
      real g1
      real hl
      integer i
      real t
      real value
      real w

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  f(x) = 1/(X*X+1),'
      write ( *, '(a)' ) '  Look at effect of upper limit T.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       W                T             VALUE           EXACT'
      write ( *, '(a)' ) ' '

      t = 5.0E+00
      w = 1.0E+00
      et = 0.0001E+00
      hl = 1.0E+00

      do i = 1, 3

        value = frcos ( f1, w, t, et, hl )

        exact = g1 ( w, t )

        write ( *, '(4(2x,g14.6))' ) w, t, value, exact

        t = t * 4.0

      end do

      return
      end
      function f1 ( x )

c*********************************************************************72
c
cc F1 evaluates F(X) = 1/(X*X+1).
c
      implicit none

      real f1
      real x

      f1 = 1.0E+00 / ( x * x + 1.0E+00 )

      return
      end
      function g1 ( w, t )

c*********************************************************************72
c
cc G1 evaluates G(W,T) = Integral ( 0 <= X <= T ) F1(X) * cos ( W * X ) dX.
c
      implicit none

      real g1
      real pi
      real t
      real w

      pi = 3.141592653589793

      g1 = pi * exp ( -w ) / 2.0E+00

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 handles F(X) = EXP(-X).
c
      implicit none

      real et
      real exact
      external f2
      real frcos
      real g2
      real hl
      integer i
      real t
      real value
      real w

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  f(x) = exp(-x),'
      write ( *, '(a)' ) '  Look at effect of upper limit T.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       W                T             VALUE           EXACT'
      write ( *, '(a)' ) ' '

      t = 5.0E+00
      w = 1.0E+00
      et = 0.0001E+00
      hl = 1.0E+00

      do i = 1, 3

        value = frcos ( f2, w, t, et, hl )

        exact = g2 ( w, t )

        write ( *, '(4(2x,g14.6))' ) w, t, value, exact

        t = t * 4.0

      end do

      return
      end
      function f2 ( x )

c*********************************************************************72
c
cc F2 evaluates F(X) = 1/(X*X+1).
c
      implicit none

      real f2
      real x

      f2 = exp ( -x )

      return
      end
      function g2 ( w, t )

c*********************************************************************72
c
cc G2 evaluates G(W,T) = Integral ( 0 <= X <= T ) F2(X) * cos ( W * X ) dX.
c
      implicit none

      real g2
      real pi
      real t
      real w

      pi = 3.141592653589793E+00

      g2 = 1.0E+00 / ( w * w + 1.0E+00 )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end



