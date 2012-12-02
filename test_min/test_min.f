      subroutine p00_f ( problem, x, f )

c*********************************************************************72
c
cc P00_F evaluates the function for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      integer problem
      double precision x

      if ( problem .eq. 1 ) then
        call p01_f ( x, f )
      else if ( problem .eq. 2 ) then
        call p02_f ( x, f )
      else if ( problem .eq. 3 ) then
        call p03_f ( x, f )
      else if ( problem .eq. 4 ) then
        call p04_f ( x, f )
      else if ( problem .eq. 5 ) then
        call p05_f ( x, f )
      else if ( problem .eq. 6 ) then
        call p06_f ( x, f )
      else if ( problem .eq. 7 ) then
        call p07_f ( x, f )
      else if ( problem .eq. 8 ) then
        call p08_f ( x, f )
      else if ( problem .eq. 9 ) then
        call p09_f ( x, f )
      else if ( problem .eq. 10 ) then
        call p10_f ( x, f )
      else if ( problem .eq. 11 ) then
        call p11_f ( x, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_F - Fatal error!'
        write ( *, '(a,i12)' ) 
     &    '  Illegal problem number PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_f1 ( problem, x, f1 )

c*********************************************************************72
c
cc P00_F1 evaluates the first derivative for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      integer problem
      double precision x

      if ( problem .eq. 1 ) then
        call p01_f1 ( x, f1 )
      else if ( problem .eq. 2 ) then
        call p02_f1 ( x, f1 )
      else if ( problem .eq. 3 ) then
        call p03_f1 ( x, f1 )
      else if ( problem .eq. 4 ) then
        call p04_f1 ( x, f1 )
      else if ( problem .eq. 5 ) then
        call p05_f1 ( x, f1 )
      else if ( problem .eq. 6 ) then
        call p06_f1 ( x, f1 )
      else if ( problem .eq. 7 ) then
        call p07_f1 ( x, f1 )
      else if ( problem .eq. 8 ) then
        call p08_f1 ( x, f1 )
      else if ( problem .eq. 9 ) then
        call p09_f1 ( x, f1 )
      else if ( problem .eq. 10 ) then
        call p10_f1 ( x, f1 )
      else if ( problem .eq. 11 ) then
        call p11_f1 ( x, f1 )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_F1 - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_f1_dif ( problem, x, f1_dif )

c*********************************************************************72
c
cc P00_F1_DIF approximates the first derivative via finite differences.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, double precision X, the point where the gradient is to 
c    be approximated.
c
c    Output, double precision F1_DIF, the approximated gradient vector.
c
      implicit none

      double precision dx
      double precision eps
      double precision f1_dif
      double precision fminus
      double precision fplus
      integer problem
      double precision r8_epsilon
      double precision x
      double precision xi

      eps = ( r8_epsilon ( ) )**0.33D+00

      if ( 0.0D+00 .le. x ) then
        dx = eps * ( x + 1.0D+00 )
      else
        dx = eps * ( x - 1.0D+00 )
      end if

      xi = x
      x = xi + dx
      call p00_f ( problem, x, fplus )

      x = xi - dx
      call p00_f ( problem, x, fminus )

      f1_dif = ( fplus - fminus ) / ( 2.0D+00 * dx )

      x = xi

      return
      end
      subroutine p00_f2 ( problem, x, f2 )

c*********************************************************************72
c
cc P00_F2 evaluates the second derivative for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      integer problem
      double precision x

      if ( problem .eq. 1 ) then
        call p01_f2 ( x, f2 )
      else if ( problem .eq. 2 ) then
        call p02_f2 ( x, f2 )
      else if ( problem .eq. 3 ) then
        call p03_f2 ( x, f2 )
      else if ( problem .eq. 4 ) then
        call p04_f2 ( x, f2 )
      else if ( problem .eq. 5 ) then
        call p05_f2 ( x, f2 )
      else if ( problem .eq. 6 ) then
        call p06_f2 ( x, f2 )
      else if ( problem .eq. 7 ) then
        call p07_f2 ( x, f2 )
      else if ( problem .eq. 8 ) then
        call p08_f2 ( x, f2 )
      else if ( problem .eq. 9 ) then
        call p09_f2 ( x, f2 )
      else if ( problem .eq. 10 ) then
        call p10_f2 ( x, f2 )
      else if ( problem .eq. 11 ) then
        call p11_f2 ( x, f2 )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_F2 - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_f2_dif ( problem, x, f2_dif )

c*********************************************************************72
c
cc P00_F2_DIF approximates the second derivative via finite differences.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F2_DIF, the approximate second derivative.
c
      implicit none

      double precision eps
      double precision f00
      double precision f2_dif
      double precision fmm
      double precision fpp
      integer problem
      double precision r8_epsilon
      double precision s
      double precision x
      double precision xi
c
c  Choose the stepsize.
c
      eps = ( r8_epsilon ( ) )**0.33D+00

      s = eps * ( abs ( x ) + 1.0D+00 )

      xi = x

      call p00_f ( problem, x, f00 )

      x = xi + s
      call p00_f ( problem, x, fpp )

      x = xi - s
      call p00_f ( problem, x, fmm )

      f2_dif = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s / s

      x = xi

      return
      end
      function p00_fmin ( a, b, problem, tol )

c*********************************************************************72
c
cc P00_FMIN seeks a minimizer of a scalar function of a scalar variable.
c
c  Discussion:
c
c    FMIN seeks an approximation to the point where F attains a minimum on
c    the interval (A,B).
c
c    The method used is a combination of golden section search and
c    successive parabolic interpolation.  Convergence is never much
c    slower than that for a Fibonacci search.  If F has a continuous
c    second derivative which is positive at the minimum (which is not
c    at A or B), then convergence is superlinear, and usually of the
c    order of about 1.324....
c
c    The function F is never evaluated at two points closer together
c    than EPS * ABS ( FMIN ) + (TOL/3), where EPS is approximately the
c    square root of the relative machine precision.  If F is a unimodal
c    function and the computed values of F are always unimodal when
c    separated by at least EPS * ABS ( XSTAR ) + (TOL/3), then FMIN
c    approximates the abcissa of the global minimum of F on the
c    interval [A, B] with an error less than 3 * EPS * ABS ( FMIN ) + TOL.
c    If F is not unimodal, then FMIN may approximate a local, but
c    perhaps non-global, minimum to the same accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization without Derivatives,
c    Prentice Hall, 1973.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1988.
c
c  Parameters
c
c    Input/output, double precision A, B.  On input, the left and right
c    endpoints of the initial interval.  On output, the lower and upper 
c    bounds for the minimizer.
c
c    Input, integer PROBLEM, the index of a problem.
c
c    Input, double precision TOL, the desired length of the interval of
c    uncertainty of the final result.  TOL must not be negative.
c
c    Output, double precision P00_FMIN, the abcissa approximating the 
c    minimizer of f.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision eps
      double precision fu
      double precision fv
      double precision fw
      double precision fx
      double precision midpoint
      integer problem
      double precision p
      double precision p00_fmin
      double precision q
      double precision r
      double precision r8_epsilon
      double precision tol
      double precision tol1
      double precision tol2
      double precision u
      double precision v
      double precision w
      double precision x

      c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )
c
c  C is the squared inverse of the golden ratio.
c
c  EPS is the square root of the relative machine precision.
c
      eps = sqrt ( r8_epsilon ( ) )
c
c  Initialization.
c
      v = a + c * ( b - a )
      w = v
      x = v
      e = 0.0D+00
      call p00_f ( problem, x, fx )
      fv = fx
      fw = fx
c
c  The main loop starts here.
c
10    continue

        midpoint = 0.5D+00 * ( a + b )
        tol1 = eps * abs ( x ) + tol / 3.0D+00
        tol2 = 2.0D+00 * tol1
c
c  Check the stopping criterion.
c
        if ( abs ( x - midpoint ) .le. 
     &    ( tol2 - 0.5D+00 * ( b - a ) ) ) then
          go to 20
        end if
c
c  Is golden-section necessary?
c
        if ( abs ( e ) .le. tol1 ) then
          if ( midpoint .le. x ) then
            e = a - x
          else
            e = b - x
          end if

          d = c * e
c
c  Consider fitting a parabola.
c
        else

          r = ( x - w ) * ( fx - fv )
          q = ( x - v ) * ( fx - fw )
          p = ( x - v ) * q - ( x - w ) * r
          q = 2.0D+00 * ( q - r )
          if ( 0.0D+00 .lt. q ) then
            p = -p
          end if
          q = abs ( q )
          r = e
          e = d
c
c  Choose a golden-section step if the parabola is not advised.
c
          if ( 
     &      ( abs ( 0.5D+00 * q * r ) .le. abs ( p ) ) .or. 
     &      ( p .le. q * ( a - x ) ) .or. 
     &      ( q * ( b - x ) .le. p ) ) then

            if ( midpoint .le. x ) then
              e = a - x
            else
              e = b - x
            end if

            d = c * e
c
c  Choose a parabolic interpolation step.
c
          else

            d = p / q
            u = x + d

            if ( ( u - a ) .lt. tol2 ) then
              d = sign ( tol1, midpoint - x )
            end if

            if ( ( b - u ) .lt. tol2 ) then
              d = sign ( tol1, midpoint - x )
            end if

         end if

       end if
c
c  F must not be evaluated too close to X.
c
        if ( tol1 .le. abs ( d ) ) then
          u = x + d
        end if

        if ( abs ( d ) .lt. tol1 ) then
          u = x + sign ( tol1, d )
        end if

        call p00_f ( problem, u, fu )
c
c  Update the data.
c
        if ( fu .le. fx ) then

          if ( x .le. u ) then
            a = x
          else
            b = x
          end if

          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
          go to 10

        end if

        if ( u .lt. x ) then
          a = u
        else
          b = u
        end if

        if ( fu .le. fw .or. w .eq. x ) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if ( fu .le. fv .or. v .eq. x .or. v .eq. w ) then
          v = u
          fv = fu
        end if

      go to 10

20    continue

      p00_fmin = x

      return
      end
      subroutine p00_interval ( problem, a, b )

c*********************************************************************72
c
cc P00_INTERVAL returns a bracketing interval for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Output, double precision A, B, two points, between which a local
c    minimizer should be sought.
c
      implicit none

      integer problem
      double precision a
      double precision b

      if ( problem .eq. 1 ) then
        call p01_interval ( a, b )
      else if ( problem .eq. 2 ) then
        call p02_interval ( a, b )
      else if ( problem .eq. 3 ) then
        call p03_interval ( a, b )
      else if ( problem .eq. 4 ) then
        call p04_interval ( a, b )
      else if ( problem .eq. 5 ) then
        call p05_interval ( a, b )
      else if ( problem .eq. 6 ) then
        call p06_interval ( a, b )
      else if ( problem .eq. 7 ) then
        call p07_interval ( a, b )
      else if ( problem .eq. 8 ) then
        call p08_interval ( a, b )
      else if ( problem .eq. 9 ) then
        call p09_interval ( a, b )
      else if ( problem .eq. 10 ) then
        call p10_interval ( a, b )
      else if ( problem .eq. 11 ) then
        call p11_interval ( a, b )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_INTERVAL - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_problem_num ( problem_num )

c*********************************************************************72
c
cc P00_PROBLEM_NUM returns the number of problems available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c   Output, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      integer problem_num

      problem_num = 11

      return
      end
      subroutine p00_sol ( problem, know, x )

c*********************************************************************72
c
cc P00_SOL returns the solution for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      integer problem
      double precision x

      if ( problem .eq. 1 ) then
        call p01_sol ( know, x )
      else if ( problem .eq. 2 ) then
        call p02_sol ( know, x )
      else if ( problem .eq. 3 ) then
        call p03_sol ( know, x )
      else if ( problem .eq. 4 ) then
        call p04_sol ( know, x )
      else if ( problem .eq. 5 ) then
        call p05_sol ( know, x )
      else if ( problem .eq. 6 ) then
        call p06_sol ( know, x )
      else if ( problem .eq. 7 ) then
        call p07_sol ( know, x )
      else if ( problem .eq. 8 ) then
        call p08_sol ( know, x )
      else if ( problem .eq. 9 ) then
        call p09_sol ( know, x )
      else if ( problem .eq. 10 ) then
        call p10_sol ( know, x )
      else if ( problem .eq. 11 ) then
        call p11_sol ( know, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_SOL - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_start ( problem, x )

c*********************************************************************72
c
cc P00_START returns a starting point for optimization for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      integer problem
      double precision x

      if ( problem .eq. 1 ) then
        call p01_start ( x )
      else if ( problem .eq. 2 ) then
        call p02_start ( x )
      else if ( problem .eq. 3 ) then
        call p03_start ( x )
      else if ( problem .eq. 4 ) then
        call p04_start ( x )
      else if ( problem .eq. 5 ) then
        call p05_start ( x )
      else if ( problem .eq. 6 ) then
        call p06_start ( x )
      else if ( problem .eq. 7 ) then
        call p07_start ( x )
      else if ( problem .eq. 8 ) then
        call p08_start ( x )
      else if ( problem .eq. 9 ) then
        call p09_start ( x )
      else if ( problem .eq. 10 ) then
        call p10_start ( x )
      else if ( problem .eq. 11 ) then
        call p11_start ( x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_START - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p00_title ( problem, title )

c*********************************************************************72
c
cc P00_TITLE returns a title for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      integer problem
      character * ( * ) title

      if ( problem .eq. 1 ) then
        call p01_title ( title )
      else if ( problem .eq. 2 ) then
        call p02_title ( title )
      else if ( problem .eq. 3 ) then
        call p03_title ( title )
      else if ( problem .eq. 4 ) then
        call p04_title ( title )
      else if ( problem .eq. 5 ) then
        call p05_title ( title )
      else if ( problem .eq. 6 ) then
        call p06_title ( title )
      else if ( problem .eq. 7 ) then
        call p07_title ( title )
      else if ( problem .eq. 8 ) then
        call p08_title ( title )
      else if ( problem .eq. 9 ) then
        call p09_title ( title )
      else if ( problem .eq. 10 ) then
        call p10_title ( title )
      else if ( problem .eq. 11 ) then
        call p11_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value of PROBLEM = ', problem
        stop
      end if

      return
      end
      subroutine p01_f ( x, f )

c*********************************************************************72
c
cc P01_F evaluates the objective function for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = ( x - 2.0D+00 ) * ( x - 2.0D+00 ) + 1.0D+00

      return
      end
      subroutine p01_f1 ( x, f1 )

c*********************************************************************72
c
cc P01_F1 evaluates the first derivative for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = 2.0D+00 * ( x - 2.0D+00 )

      return
      end
      subroutine p01_f2 ( x, f2 )

c*********************************************************************72
c
cc P01_F2 evaluates the second derivative for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 2.0D+00

      return
      end
      subroutine p01_interval ( a, b )

c*********************************************************************72
c
cc P01_INTERVAL returns a starting interval for optimization for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 3.141592653589793D+00

      return
      end
      subroutine p01_sol ( know, x )

c*********************************************************************72
c
cc P01_SOL returns the solution for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 2.0D+00

      return
      end
      subroutine p01_start ( x )

c*********************************************************************72
c
cc P01_START returns a starting point for optimization for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 3.141592653589793D+00

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns a title for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'Simple quadratic, (x-2)^2+1.'

      return
      end
      subroutine p02_f ( x, f )

c*********************************************************************72
c
cc P02_F evaluates the objective function for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = x * x + exp ( - x )

      return
      end
      subroutine p02_f1 ( x, f1 )

c*********************************************************************72
c
cc P02_F1 evaluates the first derivative for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = 2.0D+00 * x - exp ( -x )

      return
      end
      subroutine p02_f2 ( x, f2 )

c*********************************************************************72
c
cc P02_F2 evaluates the second derivative for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 2.0D+00 + exp ( -x )

      return
      end
      subroutine p02_interval ( a, b )

c*********************************************************************72
c
cc P02_INTERVAL returns a starting interval for optimization for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  0.0D+00
      b =  1.0D+00

      return
      end
      subroutine p02_sol ( know, x )

c*********************************************************************72
c
cc P02_SOL returns the solution for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 0.351734D+00

      return
      end
      subroutine p02_start ( x )

c*********************************************************************72
c
cc P02_START returns a starting point for optimization for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 0.8D+00

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns a title for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'Quadratic plus exponential, x^2 + e^(-x).'

      return
      end
      subroutine p03_f ( x, f )

c*********************************************************************72
c
cc P03_F evaluates the objective function for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = ( ( x * x + 2.0D+00 ) * x + 1.0D+00 ) * x + 3.0D+00

      return
      end
      subroutine p03_f1 ( x, f1 )

c*********************************************************************72
c
cc P03_F1 evaluates the first derivative for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = ( 4.0D+00 * x * x + 4.0D+00 ) * x + 1.0D+00

      return
      end
      subroutine p03_f2 ( x, f2 )

c*********************************************************************72
c
cc P03_F2 evaluates the second derivative for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 12.0D+00 * x * x + 4.0D+00

      return
      end
      subroutine p03_interval ( a, b )

c*********************************************************************72
c
cc P03_INTERVAL returns a starting interval for optimization for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  -2.0D+00
      b =  +2.0D+00

      return
      end
      subroutine p03_sol ( know, x )

c*********************************************************************72
c
cc P03_SOL returns the solution for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = -0.236733D+00

      return
      end
      subroutine p03_start ( x )

c*********************************************************************72
c
cc P03_START returns a starting point for optimization for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 1.5D+00

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns a title for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'Quartic, x^4 + 2x^2 + x + 3.'

      return
      end
      subroutine p04_f ( x, f )

c*********************************************************************72
c
cc P04_F evaluates the objective function for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = exp ( x ) + 0.01D+00 / x

      return
      end
      subroutine p04_f1 ( x, f1 )

c*********************************************************************72
c
cc P04_F1 evaluates the first derivative for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = exp ( x ) - 0.01D+00 / x / x

      return
      end
      subroutine p04_f2 ( x, f2 )

c*********************************************************************72
c
cc P04_F2 evaluates the second derivative for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = exp ( x ) + 0.02D+00 / x / x / x

      return
      end
      subroutine p04_interval ( a, b )

c*********************************************************************72
c
cc P04_INTERVAL returns a starting interval for optimization for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  0.0001D+00
      b =  1.0D+00

      return
      end
      subroutine p04_sol ( know, x )

c*********************************************************************72
c
cc P04_SOL returns the solution for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 0.0953446D+00

      return
      end
      subroutine p04_start ( x )

c*********************************************************************72
c
cc P04_START returns a starting point for optimization for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 0.95D+00

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns a title for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'Steep valley, e^x + 1/(100x).'

      return
      end
      subroutine p05_f ( x, f )

c*********************************************************************72
c
cc P05_F evaluates the objective function for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = exp ( x ) - 2.0D+00 * x + 0.01D+00 / x - 0.000001D+00 / x / x

      return
      end
      subroutine p05_f1 ( x, f1 )

c*********************************************************************72
c
cc P05_F1 evaluates the first derivative for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = exp ( x ) - 2.0D+00 - 0.01D+00 / x / x 
     &  + 0.000002D+00 / x / x / x

      return
      end
      subroutine p05_f2 ( x, f2 )

c*********************************************************************72
c
cc P05_F2 evaluates the second derivative for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = exp ( x ) + 0.02D+00 / x / x / x 
     &  - 0.000006D+00 / x / x / x / x

      return
      end
      subroutine p05_interval ( a, b )

c*********************************************************************72
c
cc P05_INTERVAL returns a starting interval for optimization for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  0.0002D+00
      b =  2.0D+00

      return
      end
      subroutine p05_sol ( know, x )

c*********************************************************************72
c
cc P05_SOL returns the solution for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 0.703206D+00

      return
      end
      subroutine p05_start ( x )

c*********************************************************************72
c
cc P05_START returns a starting point for optimization for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 1.5D+00

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns a title for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'Steep valley, e^x - 2x + 1/(100x) - 1/(1000000x^2).'

      return
      end
      subroutine p06_f ( x, f )

c*********************************************************************72
c
cc P06_F evaluates the objective function for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Prentice Hall 1973,
c    Reprinted Dover, 2002
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = 2.0D+00 - x

      return
      end
      subroutine p06_f1 ( x, f1 )

c*********************************************************************72
c
cc P06_F1 evaluates the first derivative for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = -1.0D+00

      return
      end
      subroutine p06_f2 ( x, f2 )

c*********************************************************************72
c
cc P06_F2 evaluates the second derivative for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    LE Scales,
c    Introduction to Non-Linear Optimization,
c    Springer, 1985.
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 0.0D+00

      return
      end
      subroutine p06_interval ( a, b )

c*********************************************************************72
c
cc P06_INTERVAL returns a starting interval for optimization for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  7.0D+00
      b =  9.0D+00

      return
      end
      subroutine p06_sol ( know, x )

c*********************************************************************72
c
cc P06_SOL returns the solution for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 9.0D+00

      return
      end
      subroutine p06_start ( x )

c*********************************************************************72
c
cc P06_START returns a starting point for optimization for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 7.2D+00

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns a title for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'line, 2 - x.'

      return
      end
      subroutine p07_f ( x, f )

c*********************************************************************72
c
cc P07_F evaluates the objective function for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Prentice Hall 1973,
c    Reprinted Dover, 2002
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = ( x + sin ( x ) ) * exp ( - x * x )

      return
      end
      subroutine p07_f1 ( x, f1 )

c*********************************************************************72
c
cc P07_F1 evaluates the first derivative for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = ( 1.0D+00 - 2.0D+00 * x * x + cos ( x ) 
     &       - 2.0D+00 * x * sin ( x ) ) * exp ( - x * x )

      return
      end
      subroutine p07_f2 ( x, f2 )

c*********************************************************************72
c
cc P07_F2 evaluates the second derivative for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = ( - 4.0D+00 - 2.0D+00 * x - 4.0D+00 * x * x * x 
     &  - 3.0D+00 * sin ( x ) - 4.0D+00 * x * cos ( x ) 
     &  + 4.0D+00 * x * x * sin ( x ) ) * exp ( - x * x )

      return
      end
      subroutine p07_interval ( a, b )

c*********************************************************************72
c
cc P07_INTERVAL returns a starting interval for optimization for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  -10.0D+00
      b =  +10.0D+00

      return
      end
      subroutine p07_sol ( know, x )

c*********************************************************************72
c
cc P07_SOL returns the solution for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = -0.6795786599525D+00

      return
      end
      subroutine p07_start ( x )

c*********************************************************************72
c
cc P07_START returns a starting point for optimization for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = -5.0D+00

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns a title for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'The dying snake, ( x + sin(x) ) * e^(-x^2).'

      return
      end
      subroutine p08_f ( x, f )

c*********************************************************************72
c
cc P08_F evaluates the objective function for problem 8.
c
c  Discussion:
c
c    This function looks positive, but has a pole at x = pi,
c    near which f -> negative infinity, and has two zeroes nearby.  
c    None of this will show up computationally.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Arnold Krommer, Christoph Ueberhuber,
c    Numerical Integration on Advanced Systems,
c    Springer, 1994, pages 185-186.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      if ( x .eq. pi ) then
        f = - 10000.0D+00
      else
        f = 3.0D+00 * x * x + 1.0D+00 
     &    + ( log ( ( x - pi ) * ( x - pi ) ) ) / pi**4
      end if

      return
      end
      subroutine p08_f1 ( x, f1 )

c*********************************************************************72
c
cc P08_F1 evaluates the first derivative for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      if ( x .eq. pi ) then
        f1 = 0.0D+00
      else
        f1 = 6.0D+00 * x + ( 2.0D+00 / ( x - pi ) ) / pi**4
      end if

      return
      end
      subroutine p08_f2 ( x, f2 )

c*********************************************************************72
c
cc P08_F2 evaluates the second derivative for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      if ( x .eq. pi ) then
        f2 = 1.0D+00
      else
        f2 = 6.0D+00 + ( - 2.0D+00 / ( x - pi ) / ( x - pi ) ) / pi**4
      end if

      return
      end
      subroutine p08_interval ( a, b )

c*********************************************************************72
c
cc P08_INTERVAL returns a starting interval for optimization for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  2.0D+00
      b =  4.0D+00

      return
      end
      subroutine p08_sol ( know, x )

c*********************************************************************72
c
cc P08_SOL returns the solution for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      know = 1

      x = pi

      return
      end
      subroutine p08_start ( x )

c*********************************************************************72
c
cc P08_START returns a starting point for optimization for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 3.1D+00

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns a title for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'The "Thin Pole", x^2+1+log((pi-x)^2)/pi^4'

      return
      end
      subroutine p09_f ( x, f )

c*********************************************************************72
c
cc P09_F evaluates the objective function for problem 9.
c
c  Discussion:
c
c    This function is oscillatory, with many local minima.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = x * x - 10.0D+00 * sin ( x * x - 3.0D+00 * x + 2.0D+00 )

      return
      end
      subroutine p09_f1 ( x, f1 )

c*********************************************************************72
c
cc P09_F1 evaluates the first derivative for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = 2.0D+00 * x 
     &  - 10.0D+00 * cos ( x * x - 3.0D+00 * x + 2.0D+00 ) 
     &  * ( 2.0D+00 * x - 3.0D+00 )

      return
      end
      subroutine p09_f2 ( x, f2 )

c*********************************************************************72
c
cc P09_F2 evaluates the second derivative for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 2.0D+00  
     &  + 10.0D+00 * sin ( x * x - 3.0D+00 * x + 2.0D+00 ) 
     &  * ( 2.0D+00 * x - 3.0D+00 ) * ( 2.0D+00 * x - 3.0D+00 ) 
     &  - 20.0D+00 * cos ( x * x - 3.0D+00 * x + 2.0D+00 )

      return
      end
      subroutine p09_interval ( a, b )

c*********************************************************************72
c
cc P09_INTERVAL returns a starting interval for optimization for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  -5.0D+00
      b =  +5.0D+00

      return
      end
      subroutine p09_sol ( know, x )

c*********************************************************************72
c
cc P09_SOL returns the solution for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 0.146621498932095D+00

      return
      end
      subroutine p09_start ( x )

c*********************************************************************72
c
cc P09_START returns a starting point for optimization for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = -2.0D+00

      return
      end
      subroutine p09_title ( title )

c*********************************************************************72
c
cc P09_TITLE returns a title for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'The oscillatory parabola'

      return
      end
      subroutine p10_f ( x, f )

c*********************************************************************72
c
cc P10_F evaluates the objective function for problem 10.
c
c  Discussion:
c
c    This function is oscillatory.
c
c    The function has a local minimum at 1.7922 whose function value is
c    very close to the minimum value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Isabel Beichl, Dianne O'Leary, Francis Sullivan,
c    Monte Carlo Minimization and Counting: One, Two, Too Many,
c    Computing in Science and Engineering,
c    Volume 9, Number 1, January/February 2007.
c
c    Dianne O'Leary,
c    Scientific Computing with Case Studies,
c    SIAM, 2008,
c    ISBN13: 978-0-898716-66-5,
c    LC: QA401.O44.
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f =           cos (           x ) 
     &  + 5.0D+00 * cos ( 1.6D+00 * x ) 
     &  - 2.0D+00 * cos ( 2.0D+00 * x ) 
     &  + 5.0D+00 * cos ( 4.5D+00 * x ) 
     &  + 7.0D+00 * cos ( 9.0D+00 * x )

      return
      end
      subroutine p10_f1 ( x, f1 )

c*********************************************************************72
c
cc P10_F1 evaluates the first derivative for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      f1 = -                     sin (           x )   
     &     - 5.0D+00 * 1.6D+00 * sin ( 1.6D+00 * x ) 
     &     + 2.0D+00 * 2.0D+00 * sin ( 2.0D+00 * x ) 
     &     - 5.0D+00 * 4.5D+00 * sin ( 4.5D+00 * x ) 
     &     - 7.0D+00 * 9.0D+00 * sin ( 9.0D+00 * x )

      return
      end
      subroutine p10_f2 ( x, f2 )

c*********************************************************************72
c
cc P10_F2 evaluates the second derivative for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = -                               cos (           x ) 
     &     - 5.0D+00 * 1.6D+00 * 1.6D+00 * cos ( 1.6D+00 * x ) 
     &     + 2.0D+00 * 2.0D+00 * 2.0D+00 * cos ( 2.0D+00 * x ) 
     &     - 5.0D+00 * 4.5D+00 * 4.5D+00 * cos ( 4.5D+00 * x ) 
     &     - 7.0D+00 * 9.0D+00 * 9.0D+00 * cos ( 9.0D+00 * x )

      return
      end
      subroutine p10_interval ( a, b )

c*********************************************************************72
c
cc P10_INTERVAL returns a starting interval for optimization for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  0.0D+00
      b =  7.0D+00

      return
      end
      subroutine p10_sol ( know, x )

c*********************************************************************72
c
cc P10_SOL returns the solution for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 5.975691087433868D+00

      return
      end
      subroutine p10_start ( x )

c*********************************************************************72
c
cc P10_START returns a starting point for optimization for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 0.5D+00

      return
      end
      subroutine p10_title ( title )

c*********************************************************************72
c
cc P10_TITLE returns a title for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = 'The cosine combo'

      return
      end
      subroutine p11_f ( x, f )

c*********************************************************************72
c
cc P11_F evaluates the objective function for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the objective function.
c
c    Output, double precision F, the value of the objective function.
c
      implicit none

      double precision f
      double precision x

      f = 1.0D+00 + abs ( 3.0D+00 * x - 1.0D+00 )

      return
      end
      subroutine p11_f1 ( x, f1 )

c*********************************************************************72
c
cc P11_F1 evaluates the first derivative for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value of the variable.
c
c    Output, double precision F1, the first derivative of the 
c    objective function.
c
      implicit none

      double precision f1
      double precision x

      if ( 3.0D+00 * x - 1.0D+00 < 0.0D+00 ) then
        f1 = - 3.0D+00
      else
        f1 = + 3.0D+00
      end if

      return
      end
      subroutine p11_f2 ( x, f2 )

c*********************************************************************72
c
cc P11_F2 evaluates the second derivative for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the values of the variables.
c
c    Output, double precision F2, the second derivative.
c
      implicit none

      double precision f2
      double precision x

      f2 = 0.0D+00

      return
      end
      subroutine p11_interval ( a, b )

c*********************************************************************72
c
cc P11_INTERVAL returns a starting interval for optimization for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, two points defining an interval in which
c    the local minimizer should be sought.
c
      implicit none

      double precision a
      double precision b

      a =  0.0D+00
      b =  1.0D+00

      return
      end
      subroutine p11_sol ( know, x )

c*********************************************************************72
c
cc P11_SOL returns the solution for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer KNOW.
c    If KNOW is 0, then the solution is not known.
c    If KNOW is positive, then the solution is known, and is returned in X.
c
c    Output, double precision X, the solution, if known.
c
      implicit none

      integer know
      double precision x

      know = 1

      x = 1.0D+00 / 3.0D+00

      return
      end
      subroutine p11_start ( x )

c*********************************************************************72
c
cc P11_START returns a starting point for optimization for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision X, a starting point for the optimization.
c
      implicit none

      double precision x

      x = 0.75D+00

      return
      end
      subroutine p11_title ( title )

c*********************************************************************72
c
cc P11_TITLE returns a title for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, a title for the problem.
c
      implicit none

      character * ( * ) title

      title = '1 + |3x-1|'

      return
      end
      function r8_add ( x, y )

c*********************************************************************72
c
cc R8_ADD adds two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the numbers to be added.
c
c    Output, double precision R8_ADD, the sum of X and Y.
c
      implicit none

      double precision r8_add
      double precision x
      double precision y

      r8_add = x + y

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision one
      double precision r8_add
      double precision r8_epsilon
      double precision temp
      double precision test
      double precision value

      save value

      data value / 0.0D+00 /

      if ( value .ne. 0.0D+00 ) then
        r8_epsilon = value
        return
      end if

      one = dble ( 1 )

      value = one
      temp = value / 2.0D+00
      test = r8_add ( one, temp )

10    continue

      if ( one .lt. test ) then
        value = temp
        temp = value / 2.0D+00
        test = r8_add ( one, temp )
        go to 10
      end if

      r8_epsilon = value

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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
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

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
