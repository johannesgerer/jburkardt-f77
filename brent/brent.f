      function glomin ( a, b, c, m, machep, e, t, f, x )

c*********************************************************************72
c
cc GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
c
c  Discussion:
c
c    This function assumes that F(X) is twice continuously differentiable
c    over [A,B] and that F''(X) <= M for all X in [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2008
c
c  Author:
c
c    Richard Brent
c    Modifications by John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the interval.
c    It must be the case that A < B.
c
c    Input, double precision C, an initial guess for the global
c    minimizer.  If no good guess is known, C = A or B is acceptable.
c
c    Input, double precision M, the bound on the second derivative.
c
c    Input, double precision MACHEP, an estimate for the relative machine
c    precision.
c
c    Input, double precision E, a positive tolerance, a bound for the
c    absolute error in the evaluation of F(X) for any X in [A,B].
c
c    Input, double precision T, a positive error tolerance.
c
c    Input, external double precision F, the name of a user-supplied
c    function, of the form "FUNCTION F ( X )", which evaluates the
c    function whose global minimum is being sought.
c
c    Output, double precision X, the estimated value of the abscissa
c    for which F attains its global minimum value in [A,B].
c
c    Output, double precision GLOMIN, the value F(X).
c
      implicit none

      double precision a
      double precision a0
      double precision a2
      double precision a3
      double precision b
      double precision c
      double precision d0
      double precision d1
      double precision d2
      double precision e
      double precision f
      double precision glomin
      double precision h
      integer k
      double precision m
      double precision m2
      double precision machep
      double precision p
      double precision q
      double precision qs
      double precision r
      double precision s
      double precision sc
      double precision t
      double precision x
      double precision y
      double precision y0
      double precision y1
      double precision y2
      double precision y3
      double precision yb
      double precision z0
      double precision z1
      double precision z2

      a0 = b
      x = a0
      a2 = a
      y0 = f ( b )
      yb = y0
      y2 = f ( a )
      y = y2

      if ( y0 < y ) then
        y = y0
      else
        x = a
      end if

      if ( m .le. 0.0D+00 .or. a .ge. b ) go to 140

      m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

      if ( sc .le. a .or. sc .ge. b ) then
        sc = 0.5D+00 * ( a + b )
      else
        sc = c
      end if

      y1 = f ( sc )
      k = 3
      d0 = a2 - sc
      h = 9.0D+00 / 11.0D+00

      if ( y1 .ge. y ) go to 30
      x = sc
      y = y1

30    continue

      d1 = a2 - a0
      d2 = sc - a0
      z2 = b - a2
      z0 = y2 - y1
      z1 = y2 - y0
      r = d1 * d1 * z0 - d0 * d0 * z1
      p = r
      qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
      q = qs
      if ( k .gt. 1000000 .and. y .lt. y2 ) go to 50

40    continue

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) .ge.
     &  z2 * m2 * r * ( z2 * q - r ) ) go to 50
      a3 = a2 + r / q
      y3 = f ( a3 )

      if ( y3 .lt. y ) then
        x = a3
        y = y3
      end if
c
c  Assume that 1611 * K does not overflow.
c
50    continue

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * dble ( k )
      if ( r .lt. z2 ) go to 40
      r = m2 * d0 * d1 * d2
      s = sqrt ( ( ( y2 - y ) + t ) / m2 )
      h = 0.5D+00 * ( 1.0D+00 + h )
      p = h * ( p + 2.0D+00 * r * s )
      q = q + 0.5D+00 * qs
      r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )
      if ( r .ge. s .and. d0 .ge. 0.0D+00 ) go to 60
      r = a2 + s
      go to 70

60    continue

      r = a2 + r

70    continue

      if ( p * q .le. 0.0D+00 ) go to 80
      a3 = a2 + p / q
      go to 90

80    continue

      a3 = r

90    continue

      if ( a3 .lt. r ) a3 = r
      if ( a3 .lt. b ) go to 100
      a3 = b
      y3 = yb
      go to 110

100   continue

      y3 = f ( a3 )

110   continue

      if ( y3 .lt. y ) then
        x = a3
        y = y3
      end if

120   continue

      d0 = a3 - a2
      if ( a3 .le. r ) go to 130
      p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )
      if ( abs ( p ) .ge. ( 1.0D+00 + 9.0D+00 * machep ) * d0 ) 
     &  go to 130
      if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) .le.
     &  ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) go to 130
      a3 = 0.5D+00 * ( a2 + a3 )
      h = 0.9D+00 * h
      go to 90

130   continue

      if ( a3 .ge. b ) go to 140
      a0 = sc
      sc = a2
      a2 = a3
      y0 = y1
      y1 = y2
      y2 = y3
      go to 30

140   continue

      glomin = y

      return
      end
      function local_min ( a, b, eps, t, f, x )

c*********************************************************************72
c
cc LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
c
c  Discussion:
c
c    The method used is a combination of golden section search and
c    successive parabolic interpolation.  Convergence is never much slower
c    than that for a Fibonacci search.  If F has a continuous second
c    derivative which is positive at the minimum (which is not at A or
c    B), then convergence is superlinear, and usually of the order of
c    about 1.324....
c
c    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
c    F is never evaluated at two points closer than TOL.  
c
c    If F is a unimodal function and the computed values of F are always
c    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
c    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
c    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  
c
c    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
c    perhaps non-global, minimum to the same accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2008
c
c  Author:
c
c    Richard Brent
c    Modifications by John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the interval.
c
c    Input, double precision EPS, a positive relative error tolerance.
c    EPS should be no smaller than twice the relative machine precision,
c    and preferably not much less than the square root of the relative
c    machine precision.
c
c    Input, double precision T, a positive absolute error tolerance.
c
c    Input, external double precision F, the name of a user-supplied
c    function, of the form "FUNCTION F ( X )", which evaluates the
c    function whose local minimum is being sought.
c
c    Output, double precision X, the estimated value of an abscissa
c    for which F attains a local minimum value in [A,B].
c
c    Output, double precision LOCAL_MIN, the value F(X).
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision eps
      double precision f
      double precision fu
      double precision fv
      double precision fw
      double precision fx
      double precision local_min
      double precision m
      double precision p
      double precision q
      double precision r
      double precision sa
      double precision sb
      double precision t
      double precision t2
      double precision tol
      double precision u
      double precision v
      double precision w
      double precision x
c
c  C is the square of the inverse of the golden ratio.
c
      c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

      sa = a
      sb = b
      x = sa + c * ( b - a )
      w = x
      v = w
      e = 0.0D+00
      fx = f ( x )
      fw = fx
      fv = fw

10    continue

      m = 0.5D+00 * ( sa + sb ) 
      tol = eps * abs ( x ) + t
      t2 = 2.0D+00 * tol
c
c  Check the stopping criterion.
c
      if ( abs ( x - m ) .le.  t2 - 0.5D+00 * ( sb - sa ) ) go to 190
      r = 0.0D+00
      q = r
      p = q
      if ( abs ( e ) .le. tol ) go to 40
c
c  Fit a parabola to the points.
c
      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )
      if ( q .le. 0.0D+00 ) go to 20
      p = - p
      go to 30

20    continue

      q = - q

30    continue

      r = e
      e = d

40    continue

      if ( abs ( p ) .ge. abs ( 0.5D+00 * q * r ) ) go to 60
      if ( p .le. q * ( sa - x ) .or. p .ge. q * ( sb - x ) ) go to 60
c
c  Take the parabolic interpolation step.
c
      d = p / q
      u = x + d
c
c  F must not be evaluated too close to A or B.
c
      if ( ( u - sa ) .ge. t2 .and. ( sb - u ) .ge. t2 ) go to 90
      if ( x .ge. m ) go to 50
      d = tol
      go to 90

50    continue

      d = - tol
      go to 90
c
c  A golden-section step.
c
60    continue

      if ( x .ge. m ) go to 70
      e = sb - x
      go to 80

70    continue

      e = a - x

80    continue

      d = c * e
c
c  F must not be evaluated too close to X.
c
90    continue

      if ( abs ( d ) .lt. tol ) go to 100
      u = x + d
      go to 120

100   continue

      if ( d .le. 0.0D+00 ) go to 110
      u = x + tol
      go to 120

110   continue

      u = x - tol

120   continue

      fu = f ( u )
c
c  Update A, B, V, W, and X.
c
      if ( fu .gt. fx ) go to 150
      if ( u .ge. x ) go to 130
      sb = x
      go to 140

130   continue

      sa = x

140   continue

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 10

150   continue

      if ( u .ge. x ) go to 160
      sa = u
      go to 170

160   continue

      sb = u

170   continue

      if ( fu .gt. fw .and. w .ne. x ) go to 180
      v = w
      fv = fw
      w = u
      fw = fu
      go to 10

180   continue

      if ( fu .gt. fv .and. v .ne. x .and. v .ne. w ) go to 20
      v = u
      fv = fu
      go to 10

190   continue

      local_min = fx

      return
      end
      subroutine local_min_rc ( a, b, arg, status, value )

c*********************************************************************72
c
cc LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
c
c  Discussion:
c
c    This routine seeks an approximation to the point where a function
c    F attains a minimum on the interval (A,B).
c
c    The method used is a combination of golden section search and
c    successive parabolic interpolation.  Convergence is never much
c    slower than that for a Fibonacci search.  If F has a continuous
c    second derivative which is positive at the minimum (which is not
c    at A or B), then convergence is superlinear, and usually of the
c    order of about 1.324...
c
c    The routine is a revised version of the Brent local minimization 
c    algorithm, using reverse communication.
c
c    It is worth stating explicitly that this routine will NOT be
c    able to detect a minimizer that occurs at either initial endpoint
c    A or B.  If this is a concern to the user, then the user must
c    either ensure that the initial interval is larger, or to check
c    the function value at the returned minimizer against the values
c    at either endpoint.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters
c
c    Input/output, double precision A, B.  On input, the left and right
c    endpoints of the initial interval.  On output, the lower and upper
c    bounds for an interval containing the minimizer.  It is required
c    that A < B.
c
c    Output, double precision ARG, the currently considered point.  The user
c    does not need to initialize this value.  On return with STATUS positive,
c    the user is requested to evaluate the function at ARG, and return
c    the value in VALUE.  On return with STATUS zero, ARG is the routine's
c    estimate for the function minimizer.
c
c    Input/output, integer STATUS, used to communicate between the user
c    and the routine.  The user only sets STATUS to zero on the first call,
c    to indicate that this is a startup call.  The routine returns STATUS
c    positive to request that the function be evaluated at ARG, or returns
c    STATUS as 0, to indicate that the iteration is complete and that
c    ARG is the estimated minimizer.
c
c    Input, double precision VALUE, the function value at ARG, as requested
c    by the routine on the previous call.
c
c  Local parameters:
c
c    C is the squared inverse of the golden ratio.
c
c    EPS is the square root of the relative machine precision.
c
      implicit none

      double precision a
      double precision arg
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
      double precision p
      double precision q
      double precision r
      double precision r8_epsilon
      integer status
      double precision tol
      double precision tol1
      double precision tol2
      double precision u
      double precision v
      double precision value
      double precision w
      double precision x

      save c
      save d
      save e
      save eps
      save fu
      save fv
      save fw
      save fx
      save midpoint
      save p
      save q
      save r
      save tol
      save tol1
      save tol2
      save u
      save v
      save w
      save x
c
c  STATUS (INPUT) = 0, startup.
c
      if ( status .eq. 0 ) then

        if ( b .le. a ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
          write ( *, '(a)' ) '  A < B is required, but'
          write ( *, '(a,g14.6)' ) '  A = ', a
          write ( *, '(a,g14.6)' ) '  B = ', b
          status = -1
          stop
        end if

        c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

        eps = sqrt ( r8_epsilon ( ) )
        tol = r8_epsilon ( )

        v = a + c * ( b - a )
        w = v
        x = v
        e = 0.0D+00

        status = 1
        arg = x

        return
c
c  STATUS (INPUT) = 1, return with initial function value of FX.
c
      else if ( status .eq. 1 ) then

        fx = value
        fv = fx
        fw = fx
c
c  STATUS (INPUT) = 2 or more, update the data.
c
      else if ( 2 .le. status ) then

        fu = value

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

        else

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

        end if

      end if
c
c  Take the next step.
c
      midpoint = 0.5D+00 * ( a + b )
      tol1 = eps * abs ( x ) + tol / 3.0D+00
      tol2 = 2.0D+00 * tol1
c
c  If the stopping criterion is satisfied, we can exit.
c
      if ( abs ( x - midpoint ) .le. 
     &   ( tol2 - 0.5D+00 * ( b - a ) ) ) then
        status = 0
        return
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
        if ( 0.0D+00 .le. q ) then
          p = - p
        end if
        q = abs ( q )
        r = e
        e = d
c
c  Choose a golden-section step if the parabola is not advised.
c
        if ( 
     &    ( abs ( 0.5D+00 * q * r ) .le. abs ( p ) ) .or. 
     &    ( p .le. q * ( a - x ) ) .or. 
     &    ( q * ( b - x ) .le. p ) ) then

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
c
c  Request value of F(U).
c
      arg = u
      status = status + 1

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

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

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
      function zero ( a, b, machep, t, f )

c*********************************************************************72
c
cc ZERO seeks the root of a function F(X) in an interval [A,B].
c
c  Discussion:
c
c    The interval [A,B] must be a change of sign interval for F.
c    That is, F(A) and F(B) must be of opposite signs.  Then
c    assuming that F is continuous implies the existence of at least
c    one value C between A and B for which F(C) = 0.
c
c    The location of the zero is determined to within an accuracy
c    of 6 * MACHEPS * abs ( C ) + 2 * T.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2008
c
c  Author:
c
c    Richard Brent
c    Modifications by John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the change of sign interval.
c
c    Input, double precision MACHEP, an estimate for the relative machine
c    precision.
c
c    Input, double precision T, a positive error tolerance.
c
c    Input, external double precision F, the name of a user-supplied
c    function, of the form "FUNCTION F ( X )", which evaluates the
c    function whose zero is being sought.
c
c    Output, double precision ZERO, the estimated value of a zero of
c    the function F.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision f
      double precision fa
      double precision fb
      double precision fc
      double precision m
      double precision machep
      double precision p
      double precision q
      double precision r
      double precision s
      double precision sa
      double precision sb
      double precision t
      double precision tol
      double precision zero
c
c  Make local copies of A and B.
c
      sa = a
      sb = b
      fa = f ( sa )
      fb = f ( sb )

10    continue

      c = sa
      fc = fa
      e = sb - sa
      d = e

20    continue

      if ( abs ( fc ) .lt. abs ( fb ) ) then
        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa
      end if

30    continue

      tol = 2.0D+00 * machep * abs ( sb ) + t
      m = 0.5D+00 * ( c - sb )
      if ( abs ( m ) .le. tol .or. fb .eq. 0.0D+00 ) go to 140
      if ( abs ( e ) .ge. tol .and. abs ( fa ) .gt. abs ( fb ) ) 
     &  go to 40

      e = m
      d = e
      go to 100

40    continue

      s = fb / fa
      if ( sa .ne. c ) go to 50

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s
      go to 60

50    continue

      q = fa / fc
      r = fb / fc
      p = s * 
     &  ( 2.0D+00 * m * a * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

60    continue

      if ( p .le. 0.0D+00 ) go to 70

      q = - q
      go to 80

70    continue

      p = - p

80    continue

      s = e
      e = d
      if ( 2.0D+00 * p .ge. 3.0D+00 * m * q - abs ( tol * q ) .or.
     &  p .ge. abs ( 0.5D+00 * s * q ) ) go to 90

      d = p / q
      go to 100

90    continue

      e = m
      d = e

100   continue

      sa = sb
      fa = fb
      if ( abs ( d ) .le. tol ) go to 110
      sb = sb + d
      go to 130

110   continue

      if ( m .le. 0.0D+00 ) go to 120
      sb = sb + tol
      go to 130

120   continue

      sb = sb - tol

130   continue

      fb = f ( sb )
      if ( fb .gt. 0.0D+00 .and. fc .gt. 0.0D+00 ) go to 10
      if ( fb .le. 0.0D+00 .and. fc .le. 0.0D+00 ) go to 10
      go to 20

140   continue

      zero = sb

      return
      end
      subroutine zero_rc ( a, b, t, arg, status, value )

c*********************************************************************72
c
cc ZERO_RC seeks the root of a function F(X) using reverse communication.
c
c  Discussion:
c
c    The interval [A,B] must be a change of sign interval for F.
c    That is, F(A) and F(B) must be of opposite signs.  Then
c    assuming that F is continuous implies the existence of at least
c    one value C between A and B for which F(C) = 0.
c
c    The location of the zero is determined to within an accuracy
c    of 6 * MACHEPS * abs ( C ) + 2 * T.
c
c    The routine is a revised version of the Brent zero finder 
c    algorithm, using reverse communication.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the change of sign interval.
c
c    Input, double precision T, a positive error tolerance.
c
c    Output, double precision ARG, the currently considered point.  The user
c    does not need to initialize this value.  On return with STATUS positive,
c    the user is requested to evaluate the function at ARG, and return
c    the value in VALUE.  On return with STATUS zero, ARG is the routine's
c    estimate for the function's zero.
c
c    Input/output, integer STATUS, used to communicate between 
c    the user and the routine.  The user only sets STATUS to zero on the first 
c    call, to indicate that this is a startup call.  The routine returns STATUS
c    positive to request that the function be evaluated at ARG, or returns
c    STATUS as 0, to indicate that the iteration is complete and that
c    ARG is the estimated zero
c
c    Input, double precision VALUE, the function value at ARG, as requested
c    by the routine on the previous call.
c
      implicit none

      double precision a
      double precision arg
      double precision b
      double precision c
      save c
      double precision d
      save d
      double precision e
      save e
      double precision fa
      save fa
      double precision fb
      save fb
      double precision fc
      save fc
      double precision m
      double precision machep
      save machep
      double precision p
      double precision q
      double precision r
      double precision r8_epsilon
      double precision s
      double precision sa
      save sa
      double precision sb
      save sb
      integer status
      double precision t
      double precision tol
      double precision value
c
c  Input STATUS = 0.
c  Initialize, request F(A).
c
      if ( status .eq. 0 ) then

        machep = r8_epsilon ( a )

        sa = a
        sb = b
        e = sb - sa
        d = e

        status = 1
        arg = a
        return
c
c  Input STATUS = 1.
c  Receive F(A), request F(B).
c
      else if ( status .eq. 1 ) then

        fa = value

        status = 2
        arg = sb
        return
c
c  Input STATUS = 2
c  Receive F(B).
c
      else if ( status .eq. 2 ) then

        fb = value

        if ( 0.0D+00 .lt. fa * fb ) then
          status = -1
          return
        end if

        c = sa
        fc = fa

      else

        fb = value

        if ( ( 0.0D+00 .lt. fb .and. 0.0D+00 .lt. fc ) .or. 
     &       ( fb .le. 0.0D+00 .and. fc .le. 0.0D+00 ) ) then
          c = sa
          fc = fa
          e = sb - sa
          d = e
        end if

      end if
c
c  Compute the next point at which a function value is requested.
c
      if ( abs ( fc ) .lt. abs ( fb ) ) then

        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa

      end if

      tol = 2.0D+00 * machep * abs ( sb ) + t
      m = 0.5D+00 * ( c - sb )

      if ( abs ( m ) .le. tol .or. fb .eq. 0.0D+00 ) then
        status = 0
        arg = sb
        return
      end if

      if ( abs ( e ) .lt. tol .or. abs ( fa ) .le. abs ( fb ) ) then

        e = m
        d = e

      else

        s = fb / fa

        if ( sa .eq. c ) then

          p = 2.0D+00 * m * s
          q = 1.0D+00 - s

        else

          q = fa / fc
          r = fb / fc
          p = s * ( 2.0D+00 * m * a * ( q - r ) 
     &      - ( sb - sa ) * ( r - 1.0D+00 ) )
          q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

        end if

        if ( 0.0D+00 .lt. p ) then
          q = - q
        else
          p = - p
        end if

        s = e
        e = d

        if ( 2.0D+00 * p .lt. 3.0D+00 * m * q - abs ( tol * q ) .and. 
     &    p .lt. abs ( 0.5D+00 * s * q ) ) then
          d = p / q
        else
          e = m
          d = e
        end if

      end if

      sa = sb
      fa = fb

      if ( tol .lt. abs ( d ) ) then
        sb = sb + d
      else if ( 0.0D+00 .lt. m ) then
        sb = sb + tol
      else
        sb = sb - tol
      end if

      arg = sb
      status = status + 1

      return
      end
