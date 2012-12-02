      program main

c*********************************************************************72
c
cc CHEBYSHEV_TEST tests CHEBYSHEV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHEBYSHEV_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the CHEBYSHEV library.'

      call chebyshev_test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHEBYSHEV_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine chebyshev_test01 ( )

c*********************************************************************72
c
cc CHEBYSHEV_TEST01 tests CHEBYSHEV_COEFFICIENTS and CHEBYSHEV_INTERPOLANT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision a
      double precision b
      double precision c(n_max)
      double precision f1
      external f1
      double precision f2
      external f2
      double precision f3
      external f3
      double precision fc(n_max)
      integer i
      integer m
      integer n
      double precision x(n_max)
      
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHEBYSHEV_TEST01'
      write ( *, '(a)' ) 
     &  '  CHEBYSHEV_COEFFICIENTS computes the coefficients of the'
      write ( *, '(a)' ) '  Chebyshev interpolant.'
      write ( *, '(a)' ) 
     &  '  CHEBYSHEV_INTERPOLANT evaluates the interpolant.'

      n = 5
      a = -1.0D+00
      b = +1.0D+00

      call chebyshev_coefficients ( a, b, n, f1, c )

      call chebyshev_zeros ( n, x )
      x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

      m = n
      call chebyshev_interpolant ( a, b, n, c, m, x, fc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X) is a trig function:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X           C(I)        F(X)       C(F)(X)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) 
     &    x(i), c(i), f1( x(i) ), fc(i)
      end do
c
c  Try a variant interval.
c
      n = 5
      a = 0.0D+00
      b = +3.0D+00

      call chebyshev_coefficients ( a, b, n, f1, c )

      call chebyshev_zeros ( n, x )
      x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

      m = n
      call chebyshev_interpolant ( a, b, n, c, m, x, fc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Consider the same F(X), but now over [0,3]:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X           C(I)        F(X)       C(F)(X)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) 
     &    x(i), c(i), f1( x(i) ), fc(i)
      end do
c
c  Try a higher order.
c
      n = 10
      a = -1.0D+00
      b = +1.0D+00

      call chebyshev_coefficients ( a, b, n, f1, c )

      call chebyshev_zeros ( n, x )
      x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

      m = n
      call chebyshev_interpolant ( a, b, n, c, m, x, fc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Consider the same F(X), but now with higher order:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X           C(I)        F(X)       C(F)(X)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) 
     &    x(i), c(i), f1( x(i) ), fc(i)
      end do
c
c  Try a polynomial.
c
      n = 10
      a = -1.0D+00
      b = +1.0D+00

      call chebyshev_coefficients ( a, b, n, f3, c )

      call chebyshev_zeros ( n, x )
      x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

      m = n
      call chebyshev_interpolant ( a, b, n, c, m, x, fc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X) is a degree 4 polynomial:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X           C(I)        F(X)       C(F)(X)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) 
     &    x(i), c(i), f3( x(i) ), fc(i)
      end do
c
c  Try a function with decaying behavior.
c
      n = 10
      a = -1.0D+00
      b = +1.0D+00

      call chebyshev_coefficients ( a, b, n, f2, c )

      call chebyshev_zeros ( n, x )
      x(1:n) = 0.5D+00 * ( a + b ) + x(1:n) * 0.5D+00 * ( b - a )

      m = n
      call chebyshev_interpolant ( a, b, n, c, m, x, fc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The polynomial approximation to F(X) decays:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X           C(I)        F(X)       C(F)(X)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6)' ) 
     &    x(i), c(i), f2( x(i) ), fc(i)
      end do

      return
      end
      function f1 ( x )

c*********************************************************************72
c
cc F1 evaluates a function that can be used for Chebyshev interpolation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, a point where the function is to be evaluated.
c
c    Output, double precision F1, the function value.
c
      implicit none

      double precision f1
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      f1 = cos ( 2.0D+00 * pi * x ) * sin ( 3.0D+00 * pi * x )

      return
      end
      function f2 ( x )

c*********************************************************************72
c
cc F2 evaluates a function that can be used for Chebyshev interpolation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input,  double precision X, a point where the function is to be evaluated.
c
c    Output, double precision F2, the function value.
c
      implicit none

      double precision f2
      double precision x

      f2 = exp ( x )

      return
      end
      function f3 ( x )

c*********************************************************************72
c
cc F3 evaluates a function that can be used for Chebyshev interpolation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, a point where the function is to be evaluated.
c
c    Output, double precision F3, the function values.
c
      implicit none

      double precision f3
      double precision x

      f3 = ( x - 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 ) 
     &  * ( x + 2.0D+00 )

      return
      end
