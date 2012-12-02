      subroutine crf ( zs, hs, hm, dm, func, ds, ze, he, de, n )

c*********************************************************************72
c
cc CRF seeks a root Z of a complex equation F(Z) = 0.
c
c  Discussion:
c
c    This routine determines a root of a transcendental complex
c    equation F(Z) = 0 by stepwise iteration, using the 
c    downhill method.
c
c    The function F(Z) must be analytic in the region where
c    the root is being sought.
c
c    The iteration continues until the stepsize has fallen below
c    the user's prescribed minimum value HM, or until the deviation
c    has fallen below the user prescribed tolerance DM (the preferred
c    condition).
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    Henning Bach
c
c  Reference:
c
c    Henning Bach,
c    On the Downhill Method,
c    Communications of the ACM,
c    Volume 121, Number 12, page 675-678, December 1969.
c
c    Henning Bach,
c    Algorithm 365: 
c    Complex Root Finding,
c    Communications of the ACM,
c    Volume 121, Number 12, page 686-687, December 1969.
c
c  Parameters:
c
c    Input, complex ZS, an initial estimate for the root.
c
c    Input, real HS, the length of the first step to try.
c
c    Input, real HM, the minimum stepsize.  The iteration will
c    terminate if the stepsize falls below this value.
c
c    Input, real DM, a tolerance for the deviation or error.
c    The iteration will terminate if the norm of the error
c    falls below this value.
c
c    Input, external function FUNC, the name of a FORTRAN function
c    which evaluates the function F, of the form:
c      function func ( z )
c      complex func
c      complex z
c
c    Output, real DS, the norm of the initial deviation, that
c    is, CABS ( FUNC ( ZS ) ).
c
c    Output, complex ZE, the estimate for the root.
c
c    Output, real HE, the step size used in the final iteration.
c
c    Output, real DE, the norm of the final deviation, that 
c    is, CABS ( FUNC ( ZE ) ).
c
c    Output, integer N, the number of iterations taken.
c
      implicit none

      complex a
      complex cw
      real de
      real dm
      real ds
      complex func
      real h
      real he
      real hm
      real hs
      integer i
      integer k
      integer n
      integer nr
      complex se
      complex u(7)
      complex v
      real w(3)
      real w0
      complex z(3)
      complex z0
      complex zd
      complex ze
      complex zs
      complex zz

      u(1) = (  1.0e+00,       0.0e+00 )
      u(2) = (  0.8660254e+00, 0.5000000e+00 )
      u(3) = (  0.0000000e+00, 1.0000000e+00 )
      u(4) = (  0.9659258e+00, 0.2588190e+00 )
      u(5) = (  0.7071068e+00, 0.7071068e+00 )
      u(6) = (  0.2588190e+00, 0.9659258e+00 )
      u(7) = ( -0.2588190e+00, 0.9659258e+00 )
      h = hs
      z0 = zs
      n = 0
c
c  Calculation of DS.
c
      cw = func ( z0 )
      w0 = abs ( real ( cw ) ) + abs ( aimag ( cw ) )
      ds = w0
      if ( w0 - dm ) 18, 18, 1

1     continue

      k = 1
      i = 0

2     continue

      v = ( -1.0e+00, 0.0e+00 )
c
c  Equilateral triangle walk pattern.
c
3     continue

      a = ( -0.5e+00, 0.866e+00 )
c
c  Calculation of deviations W in the new test points.
c
4     continue

      z(1) = z0 + h * v * a
      cw = func ( z(1) )
      w(1) = abs ( real ( cw ) ) + abs ( aimag ( cw ) )
      z(2) = z0 + h * v
      cw = func ( z(2) )
      w(2) = abs ( real ( cw ) ) + abs ( aimag ( cw ) )
      z(3) = z0 + h * conjg ( a ) * v
      cw = func ( z(3) )
      w(3) = abs ( real ( cw ) ) + abs ( aimag ( cw ) )
      n = n + 1
c
c  Determination of W(NR), the smallest of W(I).
c
      if ( w(1) - w(3) ) 5, 5, 6

5     continue

      if ( w(1) - w(2) ) 7, 8, 8

6     continue

      if ( w(2) - w(3) ) 8, 8, 9

7     continue

      nr = 1
      go to 10

8     continue

      nr = 2
      go to 10

9     continue

      nr = 3

10    continue

      if ( w0 - w(nr) ) 11, 12, 12

11    continue

      go to ( 13, 14, 15 ), k

12    continue

      k = 1
      i = 0
c
c  Forward directed walk pattern.
c
      a = ( 0.707e+00, 0.707e+00 )
      v = ( z(nr) - z0 ) / h
      w0 = w(nr)
      z0 = z(nr)
      if ( w0 - dm ) 18, 18, 4

13    continue

      k = 2
c
c  Reduction of step length.
c
      if ( h .lt. hm ) go to 18
      h = h * 0.25e+00
      go to 3

14    continue

      k = 3
c
c  Restoration of step length.
c
      h = h * 4.0e+00
      go to 2

15    continue

      i = i + 1
c
c  Rotation of walk pattern.
c
      if ( i - 7 ) 16, 16, 17

16    continue

      v = u(i)
      go to 3
c
c  Reduction of step length.
c
17    continue

      if ( h .lt. hm ) go to 18
      h = h * 0.25e+00
      i = 0
      go to 2
c
c  Terminate the iteration.
c
18    continue

      ze = z0
      he = h
      de = w0

      return
      end
