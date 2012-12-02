      program main

c*********************************************************************72
c
cc TOMS353_PRB tests FSPL2
c
c  Modified:
c
c    17 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS353_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS353 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS353_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests FSER1 with integrands of the form F(X)*COS(M*PI*X).
c
      implicit none

      real a
      real b
      real c
      real eps
      real exact
      external f1
      external f2
      external f3
      integer lc
      integer ls
      integer m
      integer max
      real pi
      real s

      pi = 3.141592653589793E+00
      a = 0.0E+00
      b = 1.0E+00

      eps = 0.00001E+00
      max = 8

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Use FSER1 to estimate the integrals of'
      write ( *, '(a)' ) '  the form F(X) * COS ( M * PI * X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Use integrand factors:'
      write ( *, '(a)' ) '  F(X) = 1, X, X*X.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '               M  Approximate         Exact'
      write ( *, '(a)' ) ' '

      m = 1

      lc = 1
      ls = 0
      call fser1 ( f1, eps, max, m, c, s, lc, ls )

      exact = ( sin ( m * pi * b ) - sin ( m * pi * a ) ) / ( m * pi )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    1:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f2, eps, max, m, c, s, lc, ls )

      exact = ( ( cos ( m * pi * b )
     &   + m * pi * b * sin ( m * pi * b ) )
     &   - ( cos ( m * pi * a ) + m * pi * a * sin ( m * pi * a ) ) )
     &   / ( m * pi )**2

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    X:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f3, eps, max, m, c, s, lc, ls )

      exact = ( ( 2.0E+00 * m * pi * b * cos ( m * pi * b )
     &   + ( ( m * pi )**2 * b**2 - 2.0E+00 ) * sin ( m * pi * b ) )
     &     - ( 2.0E+00 * m * pi * a * cos ( m * pi * a )
     &   + ( ( m * pi )**2 * a**2 - 2.0E+00 ) * sin ( m * pi * a ) ) )
     &   / ( m * pi )**3

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  X*X:  ', m, c, exact

      write ( *, '(a)' ) ' '

      m = 2

      lc = 1
      ls = 0
      call fser1 ( f1, eps, max, m, c, s, lc, ls )

      exact = ( sin ( m * pi * b ) - sin ( m * pi * a ) ) / ( m * pi )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    1:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f2, eps, max, m, c, s, lc, ls )

      exact = ( ( cos ( m * pi * b )
     &   + m * pi * b * sin ( m * pi * b ) )
     &   - ( cos ( m * pi * a ) + m * pi * a * sin ( m * pi * a ) ) )
     &   / ( m * pi )**2

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    X:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f3, eps, max, m, c, s, lc, ls )

      exact = ( ( 2.0E+00 * m * pi * b * cos ( m * pi * b )
     &   + ( ( m * pi )**2 * b**2 - 2.0E+00 ) * sin ( m * pi * b ) )
     &     - ( 2.0E+00 * m * pi * a * cos ( m * pi * a )
     &   + ( ( m * pi )**2 * a**2 - 2.0E+00 ) * sin ( m * pi * a ) ) )
     &   / ( m * pi )**3

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  X*X:  ', m, c, exact

      write ( *, '(a)' ) ' '

      m = 3

      lc = 1
      ls = 0
      call fser1 ( f1, eps, max, m, c, s, lc, ls )

      exact = ( sin ( m * pi * b ) - sin ( m * pi * a ) ) / ( m * pi )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    1:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f2, eps, max, m, c, s, lc, ls )

      exact = ( ( cos ( m * pi * b )
     &   + m * pi * b * sin ( m * pi * b ) )
     &   - ( cos ( m * pi * a ) + m * pi * a * sin ( m * pi * a ) ) )
     &   / ( m * pi )**2

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '    X:  ', m, c, exact

      lc = 1
      ls = 0
      call fser1 ( f3, eps, max, m, c, s, lc, ls )

      exact = ( ( 2.0E+00 * m * pi * b * cos ( m * pi * b )
     &   + ( ( m * pi )**2 * b**2 - 2.0E+00 ) * sin ( m * pi * b ) )
     &     - ( 2.0E+00 * m * pi * a * cos ( m * pi * a )
     &   + ( ( m * pi )**2 * a**2 - 2.0E+00 ) * sin ( m * pi * a ) ) )
     &   / ( m * pi )**3

      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  X*X:  ', m, c, exact

      return
      end
      function f1 ( x )

c*********************************************************************72
c
cc F1 evaluates the integrand factor F(X) = 1.
c
      implicit none

      real f1
      real x

      f1 = 1.0E+00

      return
      end
      function f2 ( x )

c*********************************************************************72
c
cc F2 evaluates the integrand factor F(X) = X.
c
      implicit none

      real f2
      real x

      f2 = x

      return
      end
      function f3 ( x )

c*********************************************************************72
c
cc F3 evaluates the integrand factor F(X) = X*X.
c
      implicit none

      real f3
      real x

      f3 = x*x

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
