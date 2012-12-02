      subroutine bisection ( fatol, step_max, prob, xatol, xa, xb, 
     &  fxa, fxb )

c*********************************************************************72
c
cc BISECTION carries out the bisection method to seek a root of F(X) = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for 
c    the function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, an absolute error tolerance for the root.
c
c    Input/output, double precision XA, XB, two points at which the 
c    function differs in sign.  On output, these values have been adjusted 
c    to a smaller interval.
c
c    Input/output, double precision FXA, FXB, the value of the function 
c    at XA and XB.
c 
      implicit none

      double precision fatol
      double precision fxa
      double precision fxb
      double precision fxc
      integer prob
      integer step_max
      integer step_num
      double precision t
      double precision xa
      double precision xatol
      double precision xb
      double precision xc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BISECTION'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Step      XA            XB             F(XA)         F(XB)'
      write ( *, '(a)' ) ' '
c
c  Make A the root with negative F, B the root with positive F.
c
      if ( 0.0D+00 .lt. fxa ) then
        t = xa
        xa = xb
        xb = t
        t = fxa
        fxa = fxb
        fxb = t
      end if

      step_num = 0
c
c  Loop
c
10    continue

        write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) 
     &    step_num, xa, xb, fxa, fxb

        step_num = step_num + 1

        if ( step_max .lt. step_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '  Maximum number of steps taken without convergence.'
          go to 20
        end if

        if ( dabs ( xa - xb ) .lt. xatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Interval small enough for convergence.'
          go to 20
        end if

        if ( dabs ( fxa ) .le. fatol .or. dabs ( fxb ) .le. fatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function small enough for convergence.'
          go to 20
        end if
c
c  Compute the next iterate.
c
        xc = 0.5D+00 * ( xa + xb )
        call p00_fx ( prob, xc, fxc )
c
c  Replace one of the old points.
c
        if ( fxc .lt. 0.0D+00 ) then
          xa = xc
          fxa = fxc
        else
          xb = xc
          fxb = fxc
        end if

      go to 10

20    continue

      return
      end
      subroutine brent ( fatol, step_max, prob, xatol, xrtol, xa, 
     &  xb, fxa, fxb )

c*********************************************************************72
c
cc BRENT implements the Brent bisection-based zero finder.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2001
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
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for the
c    function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, XRTOL, absolute and relative error
c    tolerances for the root.
c
c    Input/output, double precision XA, XB, two points at which the 
c    function differs in sign.  On output, these values have been adjusted
c    to a smaller interval.
c
c    Input/output, double precision FXA, FXB, the value of the function 
c    at XA and XB.
c 
      implicit none

      double precision d
      double precision e
      double precision fatol
      double precision fxa
      double precision fxb
      double precision fxc
      integer step_max
      integer prob
      integer step_num
      double precision p
      double precision q
      double precision r
      double precision s
      double precision xa
      double precision xb
      double precision xc
      double precision xm
      double precision xatol
      double precision xrtol
      double precision xtol
c
c  Initialization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BRENT'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Step      XA            XB             F(XA)         F(XB)'
      write ( *, '(a)' ) ' '

      step_num = 0

      call p00_fx ( prob, xa, fxa )
      call p00_fx ( prob, xb, fxb )
c
c  Check that f(ax) and f(bx) have different signs.
c
      if ( ( fxa .lt. 0.0D+00 .and. fxb .lt. 0.0D+00 ) .or. 
     &     ( 0.0D+00 .lt. fxa .and. 0.0D+00 .lt. fxb ) ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(XA) and F(XB) have same sign.'
        return

      end if

      xc = xa
      fxc = fxa
      d = xb - xa
      e = d

10    continue

        write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) 
     &    step_num, xb, xc, fxb, fxc

        step_num = step_num + 1
     
        if ( step_max .lt. step_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Maximum number of steps taken.'
          go to 20
        end if

        if ( dabs ( fxc ) .lt. dabs ( fxb ) ) then
          xa = xb
          xb = xc
          xc = xa
          fxa = fxb
          fxb = fxc
          fxc = fxa
        end if

        xtol = 2.0D+00 * xrtol * dabs ( xb ) + 0.5D+00 * xatol
c
c  XM is the halfwidth of the current change-of-sign interval.
c
        xm = 0.5D+00 * ( xc - xb )

        if ( dabs ( xm ) .le. xtol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Interval small enough for convergence.'
          go to 20
        end if

        if ( dabs ( fxb ) .le. fatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function small enough for convergence.'
          go to 20
        end if
c
c  See if a bisection is forced.
c
        if ( dabs ( e ) .lt. xtol .or. 
     &       dabs ( fxa ) .le. dabs ( fxb ) ) then

          d = xm
          e = d
       
        else

          s = fxb / fxa
c 
c  Linear interpolation.
c
          if ( xa .eq. xc ) then

            p = 2.0D+00 * xm * s
            q = 1.0D+00 - s
c
c  Inverse quadratic interpolation.
c
          else

            q = fxa / fxc
            r = fxb / fxc
            p = s * ( 2.0D+00 * xm * q * ( q - r ) 
     &        - ( xb - xa ) * ( r - 1.0D+00 ) )
            q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

          end if

          if ( 0.0D+00 .lt. p ) then
            q = - q
          else
            p = - p
          end if

          s = e
          e = d

          if ( 
     &      3.0D+00 * xm * q - dabs ( xtol * q ) .le. 2.0D+00 * p .or. 
     &         dabs ( 0.5D+00 * s * q ) .le. p ) then
            d = xm
            e = d
          else
            d = p / q
          end if

        end if
c
c  Save in XA, FXA the previous values of XB, FXB.
c
        xa = xb
        fxa = fxb
c
c  Compute the new value of XB, and evaluate the function there.
c
        if ( xtol .lt. dabs ( d ) ) then
          xb = xb + d
        else if ( 0.0D+00 .lt. xm ) then
          xb = xb + xtol
        else if ( xm .le. 0.0D+00 ) then
          xb = xb - xtol
        end if

        call p00_fx ( prob, xb, fxb )
c
c  If the new FXB has the same sign as FXC, then replace XC by XA.
c
        if ( ( 0.0D+00 .lt. fxb .and. 0.0D+00 .lt. fxc ) .or. 
     &       ( fxb .lt. 0.0D+00 .and. fxc .lt. 0.0D+00 ) ) then
          xc = xa
          fxc = fxa
          d = xb - xa
          e = d
        end if

      go to 10

20    continue

      return
      end
      subroutine muller ( fatol, step_max, prob, xatol, xrtol, xa, 
     &  xb, xc, fxa, fxb, fxc )

c*********************************************************************72
c
cc MULLER carries out Muller's method for seeking a real root of a nonlinear function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for the
c    function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, XRTOL, absolute and relative error
c    tolerances  for the root.
c
c    Input/output, double precision XA, XB, XC, three points.
c
c    Input/output, double precision FXA, FXB, FXC, the value of the
c    function at XA, XB, and XC.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision fatol
      integer i
      integer step_max
      integer prob
      double precision xa
      double precision xatol
      double precision xb
      double precision xc
      double precision xd
      double precision xrtol
      double precision fxa
      double precision fxb
      double precision fxc
      double precision t
      double precision z1
      double precision z2
c
c  Initialization.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULLER'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Step      XA           XB           XC'
      write ( *, '(a)' ) '          F(XA)        F(XB)        F(XC)'
      write ( *, '(a)' ) ' '

      i = 0

      write ( *, '(2x,i4,3g14.6)' ) i,  xa,  xb,  xc
      write ( *, '(2x,4x,3g14.6)' )    fxa, fxb, fxc

      do i = 1, step_max
c
c  Determine the coefficients 
c    A, B, C 
c  of the polynomial 
c    Y(X) = A * (X-X2)**2 + B * (X-X2) + C
c  which goes through the data:
c    (X1,Y1), (X2,Y2), (X3,Y3).
c
        a = ( ( fxa - fxc ) * ( xb - xc ) 
     &      - ( fxb - fxc ) * ( xa - xc ) ) / 
     &        ( ( xa - xc ) * ( xb - xc ) * ( xa - xb ) )

        b = ( ( fxb - fxc ) * ( xa - xc )**2 
     &      - ( fxa - fxc ) * ( xb - xc )**2 ) / 
     &      ( ( xa - xc ) * ( xb - xc ) * ( xa - xb ) )

        c = fxc
c
c  Get the real roots of the polynomial, 
c  unless A = 0, in which case the algorithm is breaking down.
c
        if ( a .ne. 0.0D+00 ) then

          call r8poly2_rroot ( a, b, c, z1, z2 )

        else if ( b .ne. 0.0D+00 ) then

          z2 = - c / b

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Polynomial fitting has failed.'
          write ( *, '(a)' ) '  Muller''s algorithm breaks down.'
          return

        end if

        xd = xc + z2
c
c  Set XA, YA, based on which of XA and XB is closer to XD.
c
        if ( dabs ( xd - xb ) .lt. dabs ( xd - xa ) ) then
          t = xa
          xa = xb
          xb = t
          t = fxa
          fxa = fxb
          fxb = t
        end if
c
c  Set XB, YB, based on which of XB and XC is closer to XD.
c
        if ( dabs ( xd - xc ) .lt. dabs ( xd - xb ) ) then
          t = xb
          xb = xc
          xc = t
          t = fxb
          fxb = fxc
          fxc = t
        end if
c
c  Set XC, YC.
c
        xc = xd
        call p00_fx ( prob, xc, fxc )

        write ( *, '(2x,i4,3g14.6)' ) i,  xa,  xb,  xc
        write ( *, '(2x,4x,3g14.6)' )    fxa, fxb, fxc
c
c  Estimate the relative significance of the most recent correction.
c
        if ( dabs ( z2 ) .le. xrtol * dabs ( xc ) + xatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Stepsize small enough for convergence.'
          return
        end if

        if ( dabs ( fxc ) .lt. fatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function small enough for convergence.'
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Took maximum number of steps without convergence.'

      return
      end
      subroutine newton ( fatol, step_max, prob, xatol, xmin, xmax, 
     &  xa, fxa )

c*********************************************************************72
c
cc NEWTON carries out Newton's method to seek a root of F(X) = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for the
c    function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, an absolute error tolerance for the root.
c
c    Input, double precision RANGE(2), the interval in which the root should
c    be sought.
c
c    Input/output, double precision XA.  On input, the starting point for
c    the iteration.  On output, the current approximation to the root.
c
c    Input/output, double precision FXA, the function value at XA.
c 
      implicit none

      double precision fatol
      double precision fp
      double precision fxa
      integer step_max
      integer prob
      integer step_num
      double precision step
      double precision xa
      double precision xatol
      double precision xmax
      double precision xmin

      step = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWTON'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Step       X             F(X)      FP(X)'
      write ( *, '(a)' ) ' '

      step_num = 0
      call p00_fx1 ( prob, xa, fp )

      write ( *, '(2x,i4,2x,g16.8,2g14.6)' ) step_num, xa, fxa, fp

      do step_num = 1, step_max

        if ( xa .lt. xmin .or. xmax .lt. xa ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,g14.6)' ) '  The iterate X = ', xa
          write ( *, '(a)' ) '  has left the region [XMIN,XMAX].'
          return
        end if

        if ( dabs ( fxa ) .le. fatol )then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '  The function norm is small enough for convergence.'
          return
        end if

        if ( 1 .lt. step_num .and. dabs ( step ) .le. xatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '  The stepsize is small enough for convergence.'
          return
        end if

        if ( fp .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  F''(X)=0, the algorithm fails.'
          return
        end if

        step = fxa / fp

        xa = xa - step

        call p00_fx ( prob, xa, fxa )
        call p00_fx1 ( prob, xa, fp )

        write ( *, '(2x,i4,2x,g16.8,2g14.6)' ) step_num, xa, fxa, fp

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Took maximum number of steps without convergence.'

      return
      end
      subroutine p00_fx ( prob, x, fx )

c*********************************************************************72
c
cc P00_FX evaluates a function specified by problem number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      integer prob
      double precision x

      if ( prob .eq. 1 ) then
        call p01_fx ( x, fx )
      else if ( prob .eq. 2 ) then
        call p02_fx ( x, fx )
      else if ( prob .eq. 3 ) then
        call p03_fx ( x, fx )
      else if ( prob .eq. 4 ) then
        call p04_fx ( x, fx )
      else if ( prob .eq. 5 ) then
        call p05_fx ( x, fx )
      else if ( prob .eq. 6 ) then
        call p06_fx ( x, fx )
      else if ( prob .eq. 7 ) then
        call p07_fx ( x, fx )
      else if ( prob .eq. 8 ) then
        call p08_fx ( x, fx )
      else if ( prob .eq. 9 ) then
        call p09_fx ( x, fx )
      else if ( prob .eq. 10 ) then
        call p10_fx ( x, fx )
      else if ( prob .eq. 11 ) then
        call p11_fx ( x, fx )
      else if ( prob .eq. 12 ) then
        call p12_fx ( x, fx )
      else if ( prob .eq. 13 ) then
        call p13_fx ( x, fx )
      else if ( prob .eq. 14 ) then
        call p14_fx ( x, fx )
      else if ( prob .eq. 15 ) then
        call p15_fx ( x, fx )
      else if ( prob .eq. 16 ) then
        call p16_fx ( x, fx )
      else if ( prob .eq. 17 ) then
        call p17_fx ( x, fx )
      else if ( prob .eq. 18 ) then
        call p18_fx ( x, fx )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FX - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p00_fx1 ( prob, x, fx1 )

c*********************************************************************72
c
cc P00_FX1 evaluates the first derivative of a function specified by problem number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      integer prob
      double precision x

      if ( prob .eq. 1 ) then
        call p01_fx1 ( x, fx1 )
      else if ( prob .eq. 2 ) then
        call p02_fx1 ( x, fx1 )
      else if ( prob .eq. 3 ) then
        call p03_fx1 ( x, fx1 )
      else if ( prob .eq. 4 ) then
        call p04_fx1 ( x, fx1 )
      else if ( prob .eq. 5 ) then
        call p05_fx1 ( x, fx1 )
      else if ( prob .eq. 6 ) then
        call p06_fx1 ( x, fx1 )
      else if ( prob .eq. 7 ) then
        call p07_fx1 ( x, fx1 )
      else if ( prob .eq. 8 ) then
        call p08_fx1 ( x, fx1 )
      else if ( prob .eq. 9 ) then
        call p09_fx1 ( x, fx1 )
      else if ( prob .eq. 10 ) then
        call p10_fx1 ( x, fx1 )
      else if ( prob .eq. 11 ) then
        call p11_fx1 ( x, fx1 )
      else if ( prob .eq. 12 ) then
        call p12_fx1 ( x, fx1 )
      else if ( prob .eq. 13 ) then
        call p13_fx1 ( x, fx1 )
      else if ( prob .eq. 14 ) then
        call p14_fx1 ( x, fx1 )
      else if ( prob .eq. 15 ) then
        call p15_fx1 ( x, fx1 )
      else if ( prob .eq. 16 ) then
        call p16_fx1 ( x, fx1 )
      else if ( prob .eq. 17 ) then
        call p17_fx1 ( x, fx1 )
      else if ( prob .eq. 18 ) then
        call p18_fx1 ( x, fx1 )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FX1 - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p00_fx2 ( prob, x, fx2 )

c*********************************************************************72
c
cc P00_FX2 evaluates the second derivative of a function specified by problem number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      integer prob
      double precision x

      if ( prob .eq. 1 ) then
        call p01_fx2 ( x, fx2 )
      else if ( prob .eq. 2 ) then
        call p02_fx2 ( x, fx2 )
      else if ( prob .eq. 3 ) then
        call p03_fx2 ( x, fx2 )
      else if ( prob .eq. 4 ) then
        call p04_fx2 ( x, fx2 )
      else if ( prob .eq. 5 ) then
        call p05_fx2 ( x, fx2 )
      else if ( prob .eq. 6 ) then
        call p06_fx2 ( x, fx2 )
      else if ( prob .eq. 7 ) then
        call p07_fx2 ( x, fx2 )
      else if ( prob .eq. 8 ) then
        call p08_fx2 ( x, fx2 )
      else if ( prob .eq. 9 ) then
        call p09_fx2 ( x, fx2 )
      else if ( prob .eq. 10 ) then
        call p10_fx2 ( x, fx2 )
      else if ( prob .eq. 11 ) then
        call p11_fx2 ( x, fx2 )
      else if ( prob .eq. 12 ) then
        call p12_fx2 ( x, fx2 )
      else if ( prob .eq. 13 ) then
        call p13_fx2 ( x, fx2 )
      else if ( prob .eq. 14 ) then
        call p14_fx2 ( x, fx2 )
      else if ( prob .eq. 15 ) then
        call p15_fx2 ( x, fx2 )
      else if ( prob .eq. 16 ) then
        call p16_fx2 ( x, fx2 )
      else if ( prob .eq. 17 ) then
        call p17_fx2 ( x, fx2 )
      else if ( prob .eq. 18 ) then
        call p18_fx2 ( x, fx2 )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FX2 - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of problems available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROB_NUM, the number of problems available.
c
      implicit none

      integer prob_num

      prob_num = 18

      return
      end
      subroutine p00_range ( prob, range )

c*********************************************************************72
c
cc P00_RANGE returns an interval bounding the root for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      integer prob
      double precision range(2)

      if ( prob .eq. 1 ) then
        call p01_range ( range )
      else if ( prob .eq. 2 ) then
        call p02_range ( range )
      else if ( prob .eq. 3 ) then
        call p03_range ( range )
      else if ( prob .eq. 4 ) then
        call p04_range ( range )
      else if ( prob .eq. 5 ) then
        call p05_range ( range )
      else if ( prob .eq. 6 ) then
        call p06_range ( range )
      else if ( prob .eq. 7 ) then
        call p07_range ( range )
      else if ( prob .eq. 8 ) then
        call p08_range ( range )
      else if ( prob .eq. 9 ) then
        call p09_range ( range )
      else if ( prob .eq. 10 ) then
        call p10_range ( range )
      else if ( prob .eq. 11 ) then
        call p11_range ( range )
      else if ( prob .eq. 12 ) then
        call p12_range ( range )
      else if ( prob .eq. 13 ) then
        call p13_range ( range )
      else if ( prob .eq. 14 ) then
        call p14_range ( range )
      else if ( prob .eq. 15 ) then
        call p15_range ( range )
      else if ( prob .eq. 16 ) then
        call p16_range ( range )
      else if ( prob .eq. 17 ) then
        call p17_range ( range )
      else if ( prob .eq. 18 ) then
        call p18_range ( range )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_RANGE - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p00_root ( prob, i, x )

c*********************************************************************72
c
cc P00_ROOT returns a known root for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Input, integer I, the index of the root to return.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      integer prob
      double precision x

      if ( prob .eq. 1 ) then
        call p01_root ( i, x )
      else if ( prob .eq. 2 ) then
        call p02_root ( i, x )
      else if ( prob .eq. 3 ) then
        call p03_root ( i, x )
      else if ( prob .eq. 4 ) then
        call p04_root ( i, x )
      else if ( prob .eq. 5 ) then
        call p05_root ( i, x )
      else if ( prob .eq. 6 ) then
        call p06_root ( i, x )
      else if ( prob .eq. 7 ) then
        call p07_root ( i, x )
      else if ( prob .eq. 8 ) then
        call p08_root ( i, x )
      else if ( prob .eq. 9 ) then
        call p09_root ( i, x )
      else if ( prob .eq. 10 ) then
        call p10_root ( i, x )
      else if ( prob .eq. 11 ) then
        call p11_root ( i, x )
      else if ( prob .eq. 12 ) then
        call p12_root ( i, x )
      else if ( prob .eq. 13 ) then
        call p13_root ( i, x )
      else if ( prob .eq. 14 ) then
        call p14_root ( i, x )
      else if ( prob .eq. 15 ) then
        call p15_root ( i, x )
      else if ( prob .eq. 16 ) then
        call p16_root ( i, x )
      else if ( prob .eq. 17 ) then
        call p17_root ( i, x )
      else if ( prob .eq. 18 ) then
        call p18_root ( i, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p00_root_num ( prob, root_num )

c*********************************************************************72
c
cc P00_ROOT_NUM returns the number of known roots for a problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Output, integer ROOT_NUM, the number of known roots.
c    This value may be zero.
c
      implicit none

      integer prob
      integer root_num

      if ( prob == 1 ) then
        call p01_root_num ( root_num )
      else if ( prob == 2 ) then
        call p02_root_num ( root_num )
      else if ( prob == 3 ) then
        call p03_root_num ( root_num )
      else if ( prob == 4 ) then
        call p04_root_num ( root_num )
      else if ( prob == 5 ) then
        call p05_root_num ( root_num )
      else if ( prob == 6 ) then
        call p06_root_num ( root_num )
      else if ( prob == 7 ) then
        call p07_root_num ( root_num )
      else if ( prob == 8 ) then
        call p08_root_num ( root_num )
      else if ( prob == 9 ) then
        call p09_root_num ( root_num )
      else if ( prob == 10 ) then
        call p10_root_num ( root_num )
      else if ( prob == 11 ) then
        call p11_root_num ( root_num )
      else if ( prob == 12 ) then
        call p12_root_num ( root_num )
      else if ( prob == 13 ) then
        call p13_root_num ( root_num )
      else if ( prob == 14 ) then
        call p14_root_num ( root_num )
      else if ( prob == 15 ) then
        call p15_root_num ( root_num )
      else if ( prob == 16 ) then
        call p16_root_num ( root_num )
      else if ( prob == 17 ) then
        call p17_root_num ( root_num )
      else if ( prob == 18 ) then
        call p18_root_num ( root_num )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_ROOT_NUM - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_start ( prob, i, x )

c*********************************************************************72
c
cc P00_START returns a starting point for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      integer prob
      double precision x

      if ( prob .eq. 1 ) then
        call p01_start ( i, x )
      else if ( prob .eq. 2 ) then
        call p02_start ( i, x )
      else if ( prob .eq. 3 ) then
        call p03_start ( i, x )
      else if ( prob .eq. 4 ) then
        call p04_start ( i, x )
      else if ( prob .eq. 5 ) then
        call p05_start ( i, x )
      else if ( prob .eq. 6 ) then
        call p06_start ( i, x )
      else if ( prob .eq. 7 ) then
        call p07_start ( i, x )
      else if ( prob .eq. 8 ) then
        call p08_start ( i, x )
      else if ( prob .eq. 9 ) then
        call p09_start ( i, x )
      else if ( prob .eq. 10 ) then
        call p10_start ( i, x )
      else if ( prob .eq. 11 ) then
        call p11_start ( i, x )
      else if ( prob .eq. 12 ) then
        call p12_start ( i, x )
      else if ( prob .eq. 13 ) then
        call p13_start ( i, x )
      else if ( prob .eq. 14 ) then
        call p14_start ( i, x )
      else if ( prob .eq. 15 ) then
        call p15_start ( i, x )
      else if ( prob .eq. 16 ) then
        call p16_start ( i, x )
      else if ( prob .eq. 17 ) then
        call p17_start ( i, x )
      else if ( prob .eq. 18 ) then
        call p18_start ( i, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_start_num ( prob, start_num )

c*********************************************************************72
c
cc P00_START_NUM returns the number of starting points for a problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer prob
      integer start_num

      if ( prob == 1 ) then
        call p01_start_num ( start_num )
      else if ( prob == 2 ) then
        call p02_start_num ( start_num )
      else if ( prob == 3 ) then
        call p03_start_num ( start_num )
      else if ( prob == 4 ) then
        call p04_start_num ( start_num )
      else if ( prob == 5 ) then
        call p05_start_num ( start_num )
      else if ( prob == 6 ) then
        call p06_start_num ( start_num )
      else if ( prob == 7 ) then
        call p07_start_num ( start_num )
      else if ( prob == 8 ) then
        call p08_start_num ( start_num )
      else if ( prob == 9 ) then
        call p09_start_num ( start_num )
      else if ( prob == 10 ) then
        call p10_start_num ( start_num )
      else if ( prob == 11 ) then
        call p11_start_num ( start_num )
      else if ( prob == 12 ) then
        call p12_start_num ( start_num )
      else if ( prob == 13 ) then
        call p13_start_num ( start_num )
      else if ( prob == 14 ) then
        call p14_start_num ( start_num )
      else if ( prob == 15 ) then
        call p15_start_num ( start_num )
      else if ( prob == 16 ) then
        call p16_start_num ( start_num )
      else if ( prob == 17 ) then
        call p17_start_num ( start_num )
      else if ( prob == 18 ) then
        call p18_start_num ( start_num )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_START_NUM - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_title ( prob, title )

c*********************************************************************72
c
cc P00_TITLE returns the title for a given problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the problem.
c
c    Output, character * ( * ) TITLE, the title of the given problem.
c
      implicit none

      integer prob
      character * ( * ) title

      if ( prob .eq. 1 ) then
        call p01_title ( title )
      else if ( prob .eq. 2 ) then
        call p02_title ( title )
      else if ( prob .eq. 3 ) then
        call p03_title ( title )
      else if ( prob .eq. 4 ) then
        call p04_title ( title )
      else if ( prob .eq. 5 ) then
        call p05_title ( title )
      else if ( prob .eq. 6 ) then
        call p06_title ( title )
      else if ( prob .eq. 7 ) then
        call p07_title ( title )
      else if ( prob .eq. 8 ) then
        call p08_title ( title )
      else if ( prob .eq. 9 ) then
        call p09_title ( title )
      else if ( prob .eq. 10 ) then
        call p10_title ( title )
      else if ( prob .eq. 11 ) then
        call p11_title ( title )
      else if ( prob .eq. 12 ) then
        call p12_title ( title )
      else if ( prob .eq. 13 ) then
        call p13_title ( title )
      else if ( prob .eq. 14 ) then
        call p14_title ( title )
      else if ( prob .eq. 15 ) then
        call p15_title ( title )
      else if ( prob .eq. 16 ) then
        call p16_title ( title )
      else if ( prob .eq. 17 ) then
        call p17_title ( title )
      else if ( prob .eq. 18 ) then
        call p18_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if
     
      return
      end
      subroutine p01_fx ( x, fx )

c*********************************************************************72
c
cc P01_FX evaluates sin ( x ) - x / 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = dsin ( x ) - 0.5D+00 * x

      return
      end
      subroutine p01_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P01_FX1 evaluates the derivative of the function for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = dcos ( x ) - 0.5D+00

      return
      end
      subroutine p01_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P01_FX2 evaluates the second derivative of the function for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = - dsin ( x )

      return
      end
      subroutine p01_range ( range )

c*********************************************************************72
c
cc P01_RANGE returns an interval bounding the root for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -1000.0D+00
      range(2) =  1000.0D+00

      return
      end
      subroutine p01_root ( i, x )

c*********************************************************************72
c
cc P01_ROOT returns a root for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = - 1.895494267033981D+00
      else if ( i .eq. 2 ) then
        x = 0.0D+00
      else if ( i .eq. 3 ) then
        x = 1.895494267033981D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P01_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p01_root_num ( root_num )

c*********************************************************************72
c
cc P01_ROOT_NUM returns the number of known roots for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 3

      return
      end
      subroutine p01_start ( i, x )

c*********************************************************************72
c
cc P01_START returns a starting point for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision x

      if ( i == 1 ) then
        x = 0.5D+00 * pi
      else if ( i == 2 ) then
        x = pi
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P01_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p01_start_num ( start_num )

c*********************************************************************72
c
cc P01_START_NUM returns the number of starting points for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the title of problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = SIN(X) - 0.5 * X'

      return
      end
      subroutine p02_fx ( x, fx )

c*********************************************************************72
c
cc P02_FX evaluates 2 * x - exp ( - x ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = 2.0D+00 * x - dexp ( - x )

      return
      end
      subroutine p02_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P02_FX1 evaluates the derivative of the function for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = 2.0D+00 + dexp ( - x )

      return
      end
      subroutine p02_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P02_FX2 evaluates the second derivative of the function for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = - dexp ( - x )

      return
      end
      subroutine p02_range ( range )

c*********************************************************************72
c
cc P02_RANGE returns an interval bounding the root for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -10.0D+00
      range(2) = 100.0D+00

      return
      end
      subroutine p02_root ( i, x )

c*********************************************************************72
c
cc P02_ROOT returns a root for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.35173371124919584D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P02_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p02_root_num ( root_num )

c*********************************************************************72
c
cc P02_ROOT_NUM returns the number of known roots for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p02_start ( i, x )

c*********************************************************************72
c
cc P02_START returns a starting point for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.0D+00
      else if ( i == 2 ) then
        x = 1.0D+00
      else if ( i == 3 ) then
        x = -5.0D+00
      else if ( i == 4 ) then
        x = 10.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P02_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p02_start_num ( start_num )

c*********************************************************************72
c
cc P02_START_NUM returns the number of starting points for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 4

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title of problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = 2 * X - EXP ( - X )'

      return
      end
      subroutine p03_fx ( x, fx )

c*********************************************************************72
c
cc P03_FX evaluates x * exp ( - x ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = x * dexp ( - x )

      return
      end
      subroutine p03_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P03_FX1 evaluates the derivative of the function for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = dexp ( - x ) * ( 1.0D+00 - x )

      return
      end
      subroutine p03_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P03_FX2 evaluates the second derivative of the function for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = dexp ( - x ) * ( x - 2.0D+00 )

      return
      end
      subroutine p03_range ( range )

c*********************************************************************72
c
cc P03_RANGE returns an interval bounding the root for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -10.0D+00
      range(2) = 100.0D+00

      return
      end
      subroutine p03_root ( i, x )

c*********************************************************************72
c
cc P03_ROOT returns a root for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P03_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p03_root_num ( root_num )

c*********************************************************************72
c
cc P03_ROOT_NUM returns the number of known roots for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p03_start ( i, x )

c*********************************************************************72
c
cc P03_START returns a starting point for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = -1.0D+00
      else if ( i == 2 ) then
        x =  0.5D+00
      else if ( i == 3 ) then
        x =  2.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P03_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p03_start_num ( start_num )

c*********************************************************************72
c
cc P03_START_NUM returns the number of starting points for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title of problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = X * EXP ( - X )'

      return
      end
      subroutine p04_fx ( x, fx )

c*********************************************************************72
c
cc P04_FX evaluates exp ( x ) - 1 / ( 10 * x )^2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = dexp ( x ) - 1.0D+00 / ( 100.0D+00 * x * x )

      return
      end
      subroutine p04_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P04_FX1 evaluates the derivative of the function for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = dexp ( x ) + 2.0D+00 / ( 100.0D+00 * x**3 )

      return
      end
      subroutine p04_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P04_FX2 evaluates the second derivative of the function for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = dexp ( x ) - 6.0D+00 / ( 100.0D+00 * x**4 )

      return
      end
      subroutine p04_range ( range )

c*********************************************************************72
c
cc P04_RANGE returns an interval bounding the root for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) =  0.00001D+00
      range(2) = 20.0D+00

      return
      end
      subroutine p04_root ( i, x )

c*********************************************************************72
c
cc P04_ROOT returns a root for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.09534461720025875D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P04_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p04_root_num ( root_num )

c*********************************************************************72
c
cc P04_ROOT_NUM returns the number of known roots for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p04_start ( i, x )

c*********************************************************************72
c
cc P04_START returns a starting point for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.03D+00
      else if ( i == 2 ) then
        x = 1.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P04_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p04_start_num ( start_num )

c*********************************************************************72
c
cc P04_START_NUM returns the number of starting points for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title of problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = EXP ( X ) - 1 / ( 10 * X )^2'

      return
      end
      subroutine p05_fx ( x, fx )

c*********************************************************************72
c
cc P05_FX evaluates ( x + 3 ) * ( x - 1 )^2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = ( x + 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 )

      return
      end
      subroutine p05_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P05_FX1 evaluates the derivative of the function for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = ( 3.0D+00 * x + 5.0D+00 ) * ( x - 1.0D+00 )

      return
      end
      subroutine p05_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P05_FX2 evaluates the second derivative of the function for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 6.0D+00 * x + 2.0D+00

      return
      end
      subroutine p05_range ( range )

c*********************************************************************72
c
cc P05_RANGE returns an interval bounding the root for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -1000.0D+00
      range(2) =  1000.0D+00

      return
      end
      subroutine p05_root ( i, x )

c*********************************************************************72
c
cc P05_ROOT returns a root for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = - 3.0D+00
      else if ( i .eq. 2 ) then
        x = 1.0D+00
      else if ( i .eq. 3 ) then
        x = 1.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P05_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p05_root_num ( root_num )

c*********************************************************************72
c
cc P05_ROOT_NUM returns the number of known roots for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 3

      return
      end
      subroutine p05_start ( i, x )

c*********************************************************************72
c
cc P05_START returns a starting point for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x =  2.0D+00
      else if ( i == 2 ) then
        x = - 5.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P05_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p05_start_num ( start_num )

c*********************************************************************72
c
cc P05_START_NUM returns the number of starting points for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title of problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = ( X + 3 ) * ( X - 1 )^2'

      return
      end
      subroutine p06_fx ( x, fx )

c*********************************************************************72
c
cc P06_FX evaluates exp ( x ) - 2 - 1 / ( 10 * x )^2 + 2 / ( 100 * x )^3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = dexp ( x ) - 2.0D+00 - 1.0D+00 / ( 10.0D+00 * x )**2 
     &  + 2.0D+00 / ( 100.0D+00 * x )**3

      return
      end
      subroutine p06_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P06_FX1 evaluates the derivative of the function for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = dexp ( x ) + 2.0D+00 / ( 100.0D+00 * x**3 ) 
     &  - 6.0D+00 / ( 1000000.0D+00 * x**4 )

      return
      end
      subroutine p06_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P06_FX2 evaluates the second derivative of the function for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = dexp ( x ) - 6.0D+00 / ( 100.0D+00 * x**4 ) 
     &  + 24.0D+00 / ( 1000000.0D+00 * x**5 )

      return
      end
      subroutine p06_range ( range )

c*********************************************************************72
c
cc P06_RANGE returns an interval bounding the root for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) =  0.00001D+00
      range(2) = 20.0D+00

      return
      end
      subroutine p06_root ( i, x )

c*********************************************************************72
c
cc P06_ROOT returns a root for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.7032048403631358D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P06_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p06_root_num ( root_num )

c*********************************************************************72
c
cc P06_ROOT_NUM returns the number of known roots for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p06_start ( i, x )

c*********************************************************************72
c
cc P06_START returns a starting point for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.0002D+00
      else if ( i == 2 ) then
        x = 2.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P06_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p06_start_num ( start_num )

c*********************************************************************72
c
cc P06_START_NUM returns the number of starting points for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title of problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 
     &  'F(X) = EXP(X) - 2 - 1 / ( 10 * X )^2 - 2 / ( 100 * X )^3'

      return
      end
      subroutine p07_fx ( x, fx )

c*********************************************************************72
c
cc P07_FX evaluates x^3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = x**3

      return
      end
      subroutine p07_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P07_FX1 evaluates the derivative of the function for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = 3.0D+00 * x * x

      return
      end
      subroutine p07_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P07_FX2 evaluates the second derivative of the function for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 6.0D+00 * x

      return
      end
      subroutine p07_range ( range )

c*********************************************************************72
c
cc P07_RANGE returns an interval bounding the root for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -1000.0D+00
      range(2) =  1000.0D+00

      return
      end
      subroutine p07_root ( i, x )

c*********************************************************************72
c
cc P07_ROOT returns a root for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P07_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p07_root_num ( root_num )

c*********************************************************************72
c
cc P07_ROOT_NUM returns the number of known roots for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p07_start ( i, x )

c*********************************************************************72
c
cc P07_START returns a starting point for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 1.0D+00
      else if ( i == 2 ) then
        x = -1000.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P07_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p07_start_num ( start_num )

c*********************************************************************72
c
cc P07_START_NUM returns the number of starting points for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns the title of problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = X**3, only linear Newton convergence.'

      return
      end
      subroutine p08_fx ( x, fx )

c*********************************************************************72
c
cc P08_FX evaluates cos ( x ) - x.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = dcos ( x ) - x

      return
      end
      subroutine p08_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P08_FX1 evaluates the derivative of the function for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = - dsin ( x ) - 1.0D+00

      return
      end
      subroutine p08_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P08_FX2 evaluates the second derivative of the function for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = - dcos ( x )

      return
      end
      subroutine p08_range ( range )

c*********************************************************************72
c
cc P08_RANGE returns an interval bounding the root for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 10.0D+00
      range(2) =   10.0D+00

      return
      end
      subroutine p08_root ( i, x )

c*********************************************************************72
c
cc P08_ROOT returns a root for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.7390851332151607D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P08_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p08_root_num ( root_num )

c*********************************************************************72
c
cc P08_ROOT_NUM returns the number of known roots for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p08_start ( i, x )

c*********************************************************************72
c
cc P08_START returns a starting point for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 1.0D+00
      else if ( i == 2 ) then
        x = 0.5D+00
      else if ( i == 3 ) then
        x = - 1.6D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P08_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p08_start_num ( start_num )

c*********************************************************************72
c
cc P08_START_NUM returns the number of starting points for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns the title of problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'F(X) = COS(X) - X'

      return
      end
      subroutine p09_fx ( x, fx )

c*********************************************************************72
c
cc P09_FX evaluates the Newton Baffler.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x
      double precision x2

      x2 = x - 6.25D+00

      if ( x2 .lt. - 0.25D+00 ) then
        fx = 0.75D+00 * x2 - 0.3125D+00
      else if ( x2 .lt. 0.25D+00 ) then
        fx = 2.0D+00 * x2
      else 
        fx = 0.75D+00 * x2 + 0.3125D+00
      end if

      return
      end
      subroutine p09_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P09_FX1 evaluates the derivative of the function for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x
      double precision x2

      x2 = x - 6.25D+00

      if ( x2 .lt. - 0.25D+00 ) then
        fx1 = 0.75D+00
      else if ( x2 .lt. 0.25D+00 ) then
        fx1 = 2.0D+00
      else 
        fx1 = 0.75D+00
      end if

      return
      end
      subroutine p09_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P09_FX2 evaluates the second derivative of the function for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 0.0D+00

      return
      end
      subroutine p09_range ( range )

c*********************************************************************72
c
cc P09_RANGE returns an interval bounding the root for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = -4.0D+00
      range(2) = 16.0D+00

      return
      end
      subroutine p09_root ( i, x )

c*********************************************************************72
c
cc P09_ROOT returns a root for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 6.25D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p09_root_num ( root_num )

c*********************************************************************72
c
cc P09_ROOT_NUM returns the number of known roots for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p09_start ( i, x )

c*********************************************************************72
c
cc P09_START returns a starting point for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 6.25D+00 + 5.0D+00
      else if ( i == 2 ) then
        x = 6.25D+00 - 1.0D+00
      else if ( i == 3 ) then
        x = 6.25D+00 + 0.1D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p09_start_num ( start_num )

c*********************************************************************72
c
cc P09_START_NUM returns the number of starting points for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p09_title ( title )

c*********************************************************************72
c
cc P09_TITLE returns the title of problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'The Newton Baffler'

      return
      end
      subroutine p10_fx ( x, fx )

c*********************************************************************72
c
cc P10_FX evaluates the Repeller.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = 20.0D+00 * x / ( 100.0D+00 * x * x + 1.0D+00 )

      return
      end
      subroutine p10_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P10_FX1 evaluates the derivative of the function for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = ( 1.0D+00 - 10.0D+00 * x ) * ( 1.0D+00 + 10.0D+00 * x ) 
     &  / ( 100.0D+00 * x * x + 1.0D+00 )**2

      return
      end
      subroutine p10_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P10_FX2 evaluates the second derivative of the function for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 =  - 200.0D+00 * x * ( 3.0D+00 - 100.0D+00 * x * x ) 
     &  / ( 100.0D+00 * x * x + 1.0D+00 )**3

      return
      end
      subroutine p10_range ( range )

c*********************************************************************72
c
cc P10_RANGE returns an interval bounding the root for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 10.0D+00
      range(2) = + 10.0D+00

      return
      end
      subroutine p10_root ( i, x )

c*********************************************************************72
c
cc P10_ROOT returns a root for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P10_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p10_root_num ( root_num )

c*********************************************************************72
c
cc P10_ROOT_NUM returns the number of known roots for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p10_start ( i, x )

c*********************************************************************72
c
cc P10_START returns a starting point for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 1.0D+00
      else if ( i == 2 ) then
        x = - 0.14D+00
      else if ( i == 3 ) then
        x = 0.041D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P10_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p10_start_num ( start_num )

c*********************************************************************72
c
cc P10_START_NUM returns the number of starting points for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p10_title ( title )

c*********************************************************************72
c
cc P10_TITLE returns the title of problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'The Repeller'

      return
      end
      subroutine p11_fx ( x, fx )

c*********************************************************************72
c
cc P11_FX evaluates the Pinhead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision epsilon
      double precision fx
      double precision x

      epsilon = 0.00001D+00

      if ( epsilon .eq. 0.0D+00 ) then

        fx = ( 16.0D+00 - x**4 ) / ( 16.0D+00 * x**4 )

      else
      
        fx = ( 16.0D+00 - x**4 ) / ( 16.0D+00 * x**4 + epsilon )

      end if

      return
      end
      subroutine p11_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P11_FX1 evaluates the derivative of the function for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision epsilon
      double precision fx1
      double precision x

      epsilon = 0.00001D+00

      if ( epsilon .eq. 0.0D+00 ) then

        fx1 = - 4.0D+00 / x**5

      else

        fx1 = - 4.0D+00 * x**3 * ( epsilon + 256.0D+00 ) 
     &    / ( 16.0D+00 * x**4 + epsilon )**2

      end if

      return
      end
      subroutine p11_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P11_FX2 evaluates the second derivative of the function for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision epsilon
      double precision fx2
      double precision x

      epsilon = 0.00001D+00

      if ( epsilon .eq. 0.0D+00 ) then

        fx2 =  20.0D+00 / x**6

      else

        fx2 = - 4.0D+00 * ( epsilon + 256.0D+00 ) 
     &    * ( 3.0D+00 * epsilon - 80.0D+00 * x**4 ) * x * x 
     &    / ( 16.0D+00 * x**4 + epsilon )**3

      end if

      return
      end
      subroutine p11_range ( range )

c*********************************************************************72
c
cc P11_RANGE returns an interval bounding the root for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) =    0.0D+00
      range(2) = + 10.0D+00

      return
      end
      subroutine p11_root ( i, x )

c*********************************************************************72
c
cc P11_ROOT returns a root for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = - 2.0D+00
      else if ( i .eq. 2 ) then
        x = 2.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P11_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p11_root_num ( root_num )

c*********************************************************************72
c
cc P11_ROOT_NUM returns the number of known roots for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 2

      return
      end
      subroutine p11_start ( i, x )

c*********************************************************************72
c
cc P11_START returns a starting point for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.25D+00
      else if ( i == 2 ) then
        x = 5.0D+00
      else if ( i == 3 ) then
        x = 1.1D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P11_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p11_start_num ( start_num )

c*********************************************************************72
c
cc P11_START_NUM returns the number of starting points for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p11_title ( title )

c*********************************************************************72
c
cc P11_TITLE returns the title of problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'The Pinhead'

      return
      end
      subroutine p12_fx ( x, fx )

c*********************************************************************72
c
cc P12_FX evaluates Flat Stanley.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision factor
      double precision fx
      double precision s
      double precision x
      double precision y

      factor = 1000.0D+00

      if ( x .eq. 1.0D+00 ) then

        fx = 0.0D+00

      else

        y = x - 1.0D+00
        s = dsign ( 1.0D+00, y )
        
        fx = s * dexp ( dlog ( factor ) 
     &    + dlog ( dabs ( y ) ) - 1.0D+00 / y**2 )

      end if

      return
      end
      subroutine p12_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P12_FX1 evaluates the derivative of the function for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision factor
      double precision fx1
      double precision x
      double precision y

      factor = 1000.0D+00

      if ( x .eq. 1.0D+00 ) then
        fx1 = 0.0D+00
      else
        y = x - 1.0D+00
        fx1 = factor * dexp ( - 1.0D+00 / y**2 ) 
     &    * ( y**2 + 2.0D+00 ) / y**2
      end if

      return
      end
      subroutine p12_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P12_FX2 evaluates the second derivative of the function for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision factor
      double precision fx2
      double precision x
      double precision y

      factor = 1000.0D+00

      if ( x .eq. 1.0D+00 ) then
        fx2 = 0.0D+00
      else
        y = x - 1.0D+00
        fx2 = - 2.0D+00 * factor * dexp ( - 1.0D+00 / y**2 ) 
     &    * ( y**2 - 2.0D+00 ) / y**5 
      end if

      return
      end
      subroutine p12_range ( range )

c*********************************************************************72
c
cc P12_RANGE returns an interval bounding the root for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 4.0D+00
      range(2) =   4.0D+00

      return
      end
      subroutine p12_root ( i, x )

c*********************************************************************72
c
cc P12_ROOT returns a root for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 1.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P12_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p12_root_num ( root_num )

c*********************************************************************72
c
cc P12_ROOT_NUM returns the number of known roots for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p12_start ( i, x )

c*********************************************************************72
c
cc P12_START returns a starting point for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 2.0D+00
      else if ( i == 2 ) then
        x = 0.50D+00
      else if ( i == 3 ) then
        x = 4.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P12_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p12_start_num ( start_num )

c*********************************************************************72
c
cc P12_START_NUM returns the number of starting points for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p12_title ( title )

c*********************************************************************72
c
cc P12_TITLE returns the title of problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Flat Stanley (ALL derivatives are zero at the root.)'

      return
      end
      subroutine p13_fx ( x, fx )

c*********************************************************************72
c
cc P13_FX evaluates Lazy Boy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision slope
      double precision x

      slope = 0.00000000001D+00

      fx = slope * ( x - 100.0D+00 )

      return
      end
      subroutine p13_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P13_FX1 evaluates the derivative of the function for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision slope
      double precision x

      slope = 0.00000000001D+00
      fx1 = slope

      return
      end
      subroutine p13_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P13_FX2 evaluates the second derivative of the function for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 0.0D+00

      return
      end
      subroutine p13_range ( range )

c*********************************************************************72
c
cc P13_RANGE returns an interval bounding the root for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 10000000000000.0D+00
      range(2) =   10000000000000.0D+00

      return
      end
      subroutine p13_root ( i, x )

c*********************************************************************72
c
cc P13_ROOT returns a root for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 100.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p13_root_num ( root_num )

c*********************************************************************72
c
cc P13_ROOT_NUM returns the number of known roots for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p13_start ( i, x )

c*********************************************************************72
c
cc P13_START returns a starting point for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 100000000.0D+00
      else if ( i == 2 ) then
        x = 100000013.0D+00
      else if ( i == 3 ) then
        x = - 100000000000.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p13_start_num ( start_num )

c*********************************************************************72
c
cc P13_START_NUM returns the number of starting points for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p13_title ( title )

c*********************************************************************72
c
cc P13_TITLE returns the title of problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Lazy Boy (Linear function, almost flat.)'

      return
      end
      subroutine p14_fx ( x, fx )

c*********************************************************************72
c
cc P14_FX evaluates the Camel.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx =   1.0D+00 / ( ( x - 0.3D+00 )**2 + 0.01D+00 ) 
     &     + 1.0D+00 / ( ( x - 0.9D+00 )**2 + 0.04D+00 ) 
     &     + 2.0D+00 * x - 5.2D+00

      return
      end
      subroutine p14_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P14_FX1 evaluates the derivative of the function for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = - 2.0D+00 * ( x - 0.3D+00 ) 
     &  / ( ( x - 0.3D+00 )**2 + 0.01D+00 )**2 
     &      - 2.0D+00 * ( x - 0.9D+00 ) 
     &  / ( ( x - 0.9D+00 )**2 + 0.04D+00 )**2 
     &      + 2.0D+00

      return
      end
      subroutine p14_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P14_FX2 evaluates the second derivative of the function for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 0.0D+00

      return
      end
      subroutine p14_range ( range )

c*********************************************************************72
c
cc P14_RANGE returns an interval bounding the root for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 10.0D+00
      range(2) =   10.0D+00

      return
      end
      subroutine p14_root ( i, x )

c*********************************************************************72
c
cc P14_ROOT returns a root for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = - 0.1534804948126991D+00
      else if ( i .eq. 2 ) then
        x = 1.8190323925159182D+00
      else if ( i .eq. 3 ) then
        x = 2.1274329318603367D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p14_root_num ( root_num )

c*********************************************************************72
c
cc P14_ROOT_NUM returns the number of known roots for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 3

      return
      end
      subroutine p14_start ( i, x )

c*********************************************************************72
c
cc P14_START returns a starting point for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 3.0D+00
      else if ( i == 2 ) then
        x = - 0.5D+00
      else if ( i == 3 ) then
        x = 0.0D+00
      else if ( i == 4 ) then
        x = 2.12742D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p14_start_num ( start_num )

c*********************************************************************72
c
cc P14_START_NUM returns the number of starting points for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 4

      return
      end
      subroutine p14_title ( title )

c*********************************************************************72
c
cc P14_TITLE returns the title of problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'The Camel (double hump and some shallow roots.)'

      return
      end
      subroutine p15_fx ( x, fx )

c*********************************************************************72
c
cc P15_FX evaluates a pathological function for Newton's method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    George Donovan, Arnold Miller, Timothy Moreland,
c    Pathological Functions for Newton's Method,
c    American Mathematical Monthly, January 1993, pages 53-58.
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision r8_cube_root
      double precision x

      fx = r8_cube_root ( x ) * dexp ( - x * x )

      return
      end
      subroutine p15_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P15_FX1 evaluates the derivative of the function for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision r8_cube_root
      double precision x

      fx1 = ( 1.0D+00 - 6.0D+00 * x * x ) * r8_cube_root ( x ) 
     &  * dexp ( - x * x ) / ( 3.0D+00 * x )

      return
      end
      subroutine p15_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P15_FX2 evaluates the second derivative of the function for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision r8_cube_root
      double precision x

      fx2 = ( - 2.0D+00 - 30.0D+00 * x * x + 36.0D+00 * x**4 ) 
     &  * r8_cube_root ( x ) * dexp ( - x * x ) / ( 9.0D+00 * x * x )

      return
      end
      subroutine p15_range ( range )

c*********************************************************************72
c
cc P15_RANGE returns an interval bounding the root for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = - 10.0D+00
      range(2) =   10.0D+00

      return
      end
      subroutine p15_root ( i, x )

c*********************************************************************72
c
cc P15_ROOT returns a root for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 0.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P15_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p15_root_num ( root_num )

c*********************************************************************72
c
cc P15_ROOT_NUM returns the number of known roots for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p15_start ( i, x )

c*********************************************************************72
c
cc P15_START returns a starting point for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.01D+00
      else if ( i == 2 ) then
        x = - 0.25D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P15_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p15_start_num ( start_num )

c*********************************************************************72
c
cc P15_START_NUM returns the number of starting points for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p15_title ( title )

c*********************************************************************72
c
cc P15_TITLE returns the title of problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Donovan/Miller/Moreland Pathological Function'

      return
      end
      subroutine p16_fx ( x, fx )

c*********************************************************************72
c
cc P16_FX evaluates Kepler's Equation.
c
c  Discussion:
c
c    This is Kepler's equation.  The equation has the form:
c
c      X = M + E * sin ( X )
c
c    X represents the eccentric anomaly of a planet, the angle between the
c    perihelion (the point on the orbit nearest to the sun) 
c    through the sun to the center of the ellipse, and the
c    line from the center of the ellipse to the planet.
c
c    There are two parameters:
c
c    E is the eccentricity of the orbit, which should be between 0 and 1.0;
c
c    M is the angle from the perihelion made by a fictitious planet traveling
c    on a circular orbit centered at the sun, and traveling at a constant
c    angular velocity equal to the average angular velocity of the true planet.
c    M is usually between 0 and 180 degrees, but can have any value.
c
c    For convenience, X and M are measured in degrees.
c
c    Sample results:
c
c    E        M      X
c    -----  ---  ----------
c    0.100    5    5.554589
c    0.200    5    6.246908
c    0.300    5    7.134960
c    0.400    5    8.313903
c    0.500    5    9.950063
c    0.600    5   12.356653
c    0.700    5   16.167990
c    0.800    5   22.656579
c    0.900    5   33.344447
c    0.990    5   45.361023
c    0.990    1   24.725822
c    0.990   33   89.722155
c    0.750   70  110.302
c    0.990    2   32.361007
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Peter Colwell,
c    Solving Kepler's Equation Over Three Centuries,
c    Willmann-Bell, 1993
c
c    Jean Meeus,
c    Astronomical Algorithms,
c    Willman-Bell, Inc, 1991,
c    QB51.3.E43M42
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision e
      double precision fx
      double precision m
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision x

      e = 0.8D+00
      m = 5.0D+00

      fx = ( pi * ( x - m ) / 180.0D+00 ) 
     &  - e * dsin ( pi * x / 180.0D+00 )

      return
      end
      subroutine p16_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P16_FX1 evaluates the derivative of the function for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision e
      double precision fx1
      double precision m
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision x

      e = 0.8D+00
      m = 5.0D+00

      fx1 = ( pi / 180.0D+00 ) 
     &  - e * pi * dcos ( pi * x / 180.0D+00  ) / 180.0D+00

      return
      end
      subroutine p16_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P16_FX2 evaluates the second derivative of the function for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision e
      double precision fx2
      double precision m
      double precision, parameter :: pi = 3.141592653589793D+00
      double precision x

      e = 0.8D+00
      m = 5.0D+00

      fx2 = e * pi * pi * dsin ( pi * x / 180.0D+00  ) / 180.0D+00**2

      return
      end
      subroutine p16_range ( range )

c*********************************************************************72
c
cc P16_RANGE returns an interval bounding the root for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision e
      double precision m
      double precision range(2)

      e = 0.8D+00
      m = 5.0D+00

      range(1) = m - 180.0D+00
      range(2) = m + 180.0D+00

      return
      end
      subroutine p16_root ( i, x )

c*********************************************************************72
c
cc P16_ROOT returns a root for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_ROOT - Fatal error!'
      write ( *, '(a,i4)' ) '  Illegal root index = ', i

      stop
      end
      subroutine p16_root_num ( root_num )

c*********************************************************************72
c
cc P16_ROOT_NUM returns the number of known roots for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 0

      return
      end
      subroutine p16_start ( i, x )

c*********************************************************************72
c
cc P16_START returns a starting point for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision e
      double precision m
      double precision x

      e = 0.8D+00
      m = 5.0D+00

      if ( i == 1 ) then
        x = 0.0D+00
      else if ( i == 2 ) then
        x = m
      else if ( i == 3 ) then
        x = m + 180.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P16_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p16_start_num ( start_num )

c*********************************************************************72
c
cc P16_START_NUM returns the number of starting points for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 3

      return
      end
      subroutine p16_title ( title )

c*********************************************************************72
c
cc P16_TITLE returns the title of problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Kepler''s Eccentric Anomaly Equation, in degrees'

      return
      end
      subroutine p17_fx ( x, fx )

c*********************************************************************72
c
cc P17_FX evaluates the function for problem 17.
c
c  Discussion:
c
c    This simple example is of historical interest, since it was used
c    by Wallis to illustrate the use of Newton's method, and has been
c    a common example ever since.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = x**3 - 2.0D+00 * x - 5.0D+00

      return
      end
      subroutine p17_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P17_FX1 evaluates the derivative of the function for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = 3.0D+00 * x * x - 2.0D+00

      return
      end
      subroutine p17_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P17_FX2 evaluates the second derivative of the function for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 6.0D+00 * x

      return
      end
      subroutine p17_range ( range )

c*********************************************************************72
c
cc P17_RANGE returns an interval bounding the root for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = 2.0D+00
      range(2) = 3.0D+00

      return
      end
      subroutine p17_root ( i, x )

c*********************************************************************72
c
cc P17_ROOT returns a root for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 2.0945514815423265D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P17_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p17_root_num ( root_num )

c*********************************************************************72
c
cc P17_ROOT_NUM returns the number of known roots for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p17_start ( i, x )

c*********************************************************************72
c
cc P17_START returns a starting point for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 2.0D+00
      else if ( i == 2 ) then
        x = 3.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P17_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p17_start_num ( start_num )

c*********************************************************************72
c
cc P17_START_NUM returns the number of starting points for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p17_title ( title )

c*********************************************************************72
c
cc P17_TITLE returns the title of problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'The Wallis example, x^3-2x-5=0'

      return
      end
      subroutine p18_fx ( x, fx )

c*********************************************************************72
c
cc P18_FX evaluates the function for problem 18.
c
c  Discussion:
c
c    F(X) = 10^14 * (x-1)^7, but is written in term by term form.
c
c    This polynomial becomes difficult to evaluate accurately when 
c    written this way.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point at which F is to be evaluated.
c
c    Output, double precision FX, the value of the function at X.
c
      implicit none

      double precision fx
      double precision x

      fx = 10.0D+00**14 * ( 
     &                 x**7 
     &    -  7.0D+00 * x**6 
     &    + 21.0D+00 * x**5 
     &    - 35.0D+00 * x**4 
     &    + 35.0D+00 * x**3 
     &    - 21.0D+00 * x**2 
     &    +  7.0D+00 * x    
     &    -  1.0D+00 )

      return
      end
      subroutine p18_fx1 ( x, fx1 )

c*********************************************************************72
c
cc P18_FX1 evaluates the derivative of the function for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX1, the first derivative of the function at X.
c
      implicit none

      double precision fx1
      double precision x

      fx1 = 10.0D+00**14 * ( 
     &        7.0D+00 * x**6 
     &    -  42.0D+00 * x**5 
     &    + 105.0D+00 * x**4 
     &    - 140.0D+00 * x**3 
     &    + 105.0D+00 * x**2 
     &    -  42.0D+00 * x    
     &    +   7.0D+00 )

      return
      end
      subroutine p18_fx2 ( x, fx2 )

c*********************************************************************72
c
cc P18_FX2 evaluates the second derivative of the function for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the abscissa.
c
c    Output, double precision FX2, the second derivative of the function at X.
c
      implicit none

      double precision fx2
      double precision x

      fx2 = 10.0D+00**14 * ( 
     &       42.0D+00 * x**5 
     &    - 210.0D+00 * x**4 
     &    + 420.0D+00 * x**3 
     &    - 420.0D+00 * x**2 
     &    + 210.0D+00 * x    
     &    -  42.0D+00 )

      return
      end
      subroutine p18_range ( range )

c*********************************************************************72
c
cc P18_RANGE returns an interval bounding the root for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision RANGE(2), the minimum and maximum values of 
c    an interval containing the root.
c
      implicit none

      double precision range(2)

      range(1) = 0.988D+00
      range(2) = 1.012D+00

      return
      end
      subroutine p18_root ( i, x )

c*********************************************************************72
c
cc P18_ROOT returns a root for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the requested root.
c
c    Output, double precision X, the I-th root.
c
      implicit none

      integer i
      double precision x

      if ( i .eq. 1 ) then
        x = 1.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P18_ROOT - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal root index = ', i
        stop
      end if

      return
      end
      subroutine p18_root_num ( root_num )

c*********************************************************************72
c
cc P18_ROOT_NUM returns the number of known roots for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ROOT_NUM, the number of known roots.
c
      implicit none

      integer root_num

      root_num = 1

      return
      end
      subroutine p18_start ( i, x )

c*********************************************************************72
c
cc P18_START returns a starting point for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the starting point.
c
c    Output, double precision X, the starting point.
c
      implicit none

      integer i
      double precision x

      if ( i == 1 ) then
        x = 0.990D+00
      else if ( i == 2 ) then
        x = 1.013D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P18_START - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal start index = ', i
        stop
      end if

      return
      end
      subroutine p18_start_num ( start_num )

c*********************************************************************72
c
cc P18_START_NUM returns the number of starting points for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer START_NUM, the number of starting points.
c
      implicit none

      integer start_num

      start_num = 2

      return
      end
      subroutine p18_title ( title )

c*********************************************************************72
c
cc P18_TITLE returns the title of problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = '10^14 * (x-1)^7, written term by term.'

      return
      end
      function r8_cube_root ( x )

c*********************************************************************72
c
cc R8_CUBE_ROOT returns the cube root of an R8.
c
c  Discussion:
c
c    This routine is designed to avoid the possible problems that can occur
c    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose cube root is desired.
c
c    Output, double precision R8_CUBE_ROOT, the cube root of X.
c
      implicit none

      double precision r8_cube_root
      double precision value
      double precision x

      if ( 0.0D+00 .lt. x ) then
        value = x**(1.0D+00/3.0D+00)
      else if ( x .eq. 0.0D+00 ) then
        value = 0.0D+00
      else
        value = - ( dabs ( x ) )**(1.0D+00/3.0D+00)
      end if

      r8_cube_root = value

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_sign ( x )

c*********************************************************************72
c
cc R8_SIGN returns the sign of an R8.
c
c  Discussion:
c
c    value = -1 if X < 0;
c    value =  0 if X => 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose sign is desired.
c
c    Output, double precision R8_SIGN, the sign of X:
c
      implicit none

      double precision r8_sign
      double precision x

      if ( x .lt. 0.0D+00 ) then
        r8_sign = -1.0D+00
      else
        r8_sign = +1.0D+00
      end if

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8poly2_rroot ( a, b, c, r1, r2 )

c*********************************************************************72
c
cc R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
c
c  Example:
c
c    A    B    C       roots              R1   R2
c   --   --   --     ------------------   --   --
c    1   -4    3     1          3          1    3
c    1    0    4         2*i      - 2*i    0    0
c    2   -6    5     3 +   i    3 -   i    3    3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, C, the coefficients of the quadratic
c    polynomial A * X * X + B * X + C = 0 whose roots are desired.  
c    A must not be zero.
c
c    Output, double precision R1, R2, the real parts of the two roots 
c    of the polynomial.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision disc
      double precision q
      double precision r1
      double precision r2

      if ( a .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8POLY2_RROOT - Fatal error!'
        write ( *, '(a)' ) '  The coefficient A is zero.'
        stop
      end if

      disc = b * b - 4.0D+00 * a * c
      disc = dmax1 ( disc, 0.0D+00 )
         
      q = ( b + dsign ( 1.0D+00, b ) * dsqrt ( disc ) )
      r1 = - 0.5D+00 * q / a
      r2 = - 2.0D+00 * c / q

      return
      end
      subroutine regula_falsi ( fatol, step_max, prob, xatol, xa, 
     &  xb, fxa, fxb )

c*********************************************************************72
c
cc REGULA_FALSI carries out the Regula Falsi method to seek a root of F(X) = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for the
c    function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, an absolute error tolerance for the root.
c
c    Input/output, double precision XA, XB, two points at which the 
c    function differs in sign.  On output, these values have been adjusted
c    to a smaller interval.
c
c    Input/output, double precision FXA, FXB, the value of the function 
c    at XA and XB.
c 
      implicit none

      double precision fatol
      double precision fxa
      double precision fxb
      double precision fxc
      integer step_max
      integer prob
      integer step_num
      double precision t
      double precision xa
      double precision xb
      double precision xatol
      double precision xc
c
c  The method requires a change-of-sign interval.
c
      if ( dsign ( 1.0D+00, fxa ) .eq. dsign ( 1.0D+00, fxb ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'REGULA_FALSI - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Function does not change sign at endpoints.'
        stop
      end if
c
c  Make A the root with negative F, B the root with positive F.
c
      if ( 0.0D+00 .lt. fxa ) then
        t = xa
        xa = xb
        xb = t
        t = fxa
        fxa = fxb
        fxb = t
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGULA FALSI'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Step      XA            XB             F(XA)         F(XB)'
      write ( *, '(a)' ) ' '

      step_num = 0
      write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) 
     &  step_num, xa, xb, fxa, fxb

      do step_num = 1, step_max

        if ( dabs ( xa - xb ) .lt. xatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Interval small enough for convergence.'
          return
        end if

        if ( dabs ( fxa ) .le. fatol .or. dabs ( fxb ) .le. fatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function small enough for convergence.'
          return
        end if

        xc = ( fxa * xb - fxb * xa ) / ( fxa - fxb )
        call p00_fx ( prob, xc, fxc )

        if ( fxc .lt. 0.0D+00 ) then
          xa = xc
          fxa = fxc
        else
          xb = xc
          fxb = fxc
        end if

        write ( *, '(2x,i4,2x,2g16.8,2g14.6)' ) 
     &    step_num, xa, xb, fxa, fxb

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Took maximum number of steps without convergence.'

      return
      end
      subroutine secant ( fatol, step_max, prob, xatol, xmin, xmax, 
     &  xa, xb, fxa, fxb )

c*********************************************************************72
c
cc SECANT carries out the secant method to seek a root of F(X) = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FATOL, an absolute error tolerance for the
c    function value of the root.  If an approximate root X satisfies
c      ABS ( F ( X ) ) .le. FATOL, then X will be accepted as the
c    root and the iteration will be terminated.
c
c    Input, integer STEP_MAX, the maximum number of steps allowed
c    for an iteration.
c
c    Input, integer PROB, the index of the function whose root is
c    to be sought.
c
c    Input, double precision XATOL, an absolute error tolerance for the root.
c
c    Input, double precision XMAX, XMIN, the interval in which the root should
c    be sought.
c
c    Input/output, double precision XA, XB, two points at which the 
c    function differs in sign.  On output, these values have been adjusted
c    to a smaller interval.
c
c    Input/output, double precision FXA, FXB, the value of the function 
c    at XA and XB.
c 
      implicit none

      double precision fatol
      double precision fxa
      double precision fxb
      double precision fxc
      integer step_max
      integer prob
      integer step_num
      double precision xa
      double precision xatol
      double precision xb
      double precision xc
      double precision xmax
      double precision xmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SECANT'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Step       X             F(X)'
      write ( *, '(a)' ) ' '

      step_num = -1
      write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xa, fxa

      if ( dabs ( fxa ) .le. fatol ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Function small enough for convergence.'
        return
      end if

      step_num = 0
      write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xb, fxb

      do step_num = 1, step_max

        if ( dabs ( fxb ) .le. fatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function small enough for convergence.'
          return
        end if

        if ( dabs ( xa - xb ) .lt. xatol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Interval small enough for convergence.'
          return
        end if

        if ( xb .lt. xmin .or. xmax .lt. xb ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '  Iterate has left the region [XMIN,XMAX].'
          return
        end if

        if ( fxa .eq. fxb ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  F(A) = F(B), algorithm fails.'
          return
        end if

        xc = ( fxa * xb - fxb * xa ) / ( fxa - fxb )

        call p00_fx ( prob, xc, fxc )

        xa = xb
        fxa = fxb
        xb = xc
        fxb = fxc
        write ( *, '(2x,i4,2x,g16.8,g14.6)' ) step_num, xb, fxb

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Took maximum number of steps.'

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
