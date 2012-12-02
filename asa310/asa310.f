      function alngam ( xvalue, ifault )

c*********************************************************************72
c
cc ALNGAM computes the logarithm of the gamma function.
c
c  Modified:
c
c    30 March 1999
c
c  Author:
c
c    Allan Macleod
c    Modifications by John Burkardt
c
c  Reference:
c
c    Allan Macleod,
c    Algorithm AS 245,
c    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
c    Applied Statistics,
c    Volume 38, Number 2, 1989, pages 397-402.
c
c  Parameters:
c
c    Input, double precision XVALUE, the argument of the Gamma function.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    1, XVALUE is less than or equal to 0.
c    2, XVALUE is too big.
c
c    Output, double precision ALNGAM, the logarithm of the gamma function of X.
c
      implicit none

      double precision alngam
      double precision alr2pi
      parameter ( alr2pi = 0.918938533204673D+00 )
      integer ifault
      double precision r1(9)
      double precision r2(9)
      double precision r3(9)
      double precision r4(5)
      double precision x
      double precision x1
      double precision x2
      double precision xlge
      parameter ( xlge = 5.10D+06 )
      double precision xlgst
      parameter ( xlgst = 1.0D+30 )
      double precision xvalue
      double precision y

      data r1 /
     &  -2.66685511495D+00,
     &  -24.4387534237D+00,
     &  -21.9698958928D+00,
     &   11.1667541262D+00,
     &   3.13060547623D+00,
     &   0.607771387771D+00,
     &   11.9400905721D+00,
     &   31.4690115749D+00,
     &   15.2346874070D+00 /

      data r2 /
     &  -78.3359299449D+00,
     &  -142.046296688D+00,
     &   137.519416416D+00,
     &   78.6994924154D+00,
     &   4.16438922228D+00,
     &   47.0668766060D+00,
     &   313.399215894D+00,
     &   263.505074721D+00,
     &   43.3400022514D+00 /

      data r3 /
     &  -2.12159572323D+05,
     &   2.30661510616D+05,
     &   2.74647644705D+04,
     &  -4.02621119975D+04,
     &  -2.29660729780D+03,
     &  -1.16328495004D+05,
     &  -1.46025937511D+05,
     &  -2.42357409629D+04,
     &  -5.70691009324D+02 /

      data r4 / 
     &   0.279195317918525D+00, 
     &   0.4917317610505968D+00,
     &   0.0692910599291889D+00, 
     &   3.350343815022304D+00,
     &   6.012459259764103D+00 /

      x = xvalue
      alngam = 0.0D+00
c
c  Check the input.
c
      if ( xlgst .le. x ) then
        ifault = 2
        return
      end if

      if ( x .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0
c
c  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
c
      if ( x .lt. 1.5D+00 ) then

        if ( x .lt. 0.5D+00 ) then

          alngam = - dlog ( x )
          y = x + 1.0D+00
c
c  Test whether X < machine epsilon.
c
          if ( y .eq. 1.0D+00 ) then
            return
          end if

        else

          alngam = 0.0D+00
          y = x
          x = ( x - 0.5D+00 ) - 0.5D+00

        end if

        alngam = alngam + x * ((((
     &      r1(5)   * y 
     &    + r1(4) ) * y 
     &    + r1(3) ) * y
     &    + r1(2) ) * y 
     &    + r1(1) ) / ((((
     &                y 
     &    + r1(9) ) * y 
     &    + r1(8) ) * y
     &    + r1(7) ) * y 
     &    + r1(6) )

        return

      end if
c
c  Calculation for 1.5 <= X < 4.0.
c
      if ( x .lt. 4.0D+00 ) then

        y = ( x - 1.0D+00 ) - 1.0D+00

        alngam = y * ((((
     &      r2(5)   * x 
     &    + r2(4) ) * x 
     &    + r2(3) ) * x 
     &    + r2(2) ) * x
     &    + r2(1) ) / ((((
     &                x 
     &    + r2(9) ) * x 
     &    + r2(8) ) * x 
     &    + r2(7) ) * x
     &    + r2(6) )
c
c  Calculation for 4.0 <= X < 12.0.
c
      else if ( x .lt. 12.0D+00 ) then

        alngam = ((((
     &      r3(5)   * x 
     &    + r3(4) ) * x 
     &    + r3(3) ) * x 
     &    + r3(2) ) * x 
     &    + r3(1) ) / (((( 
     &                x 
     &    + r3(9) ) * x 
     &    + r3(8) ) * x 
     &    + r3(7) ) * x 
     &    + r3(6) )
c
c  Calculation for X >= 12.0.
c
      else

        y = dlog ( x )
        alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

        if ( x .le. xlge ) then

          x1 = 1.0D+00 / x
          x2 = x1 * x1

          alngam = alngam + x1 * ( ( 
     &           r4(3)   * 
     &      x2 + r4(2) ) * 
     &      x2 + r4(1) ) / ( ( 
     &      x2 + r4(5) ) * 
     &      x2 + r4(4) )

        end if

      end if

      return
      end
      function alnorm ( x, upper )

c*********************************************************************72
c
cc ALNORM computes the cumulative density of the standard normal distribution.
c
c  Modified:
c
c    28 March 1999
c
c  Author:
c
c    David Hill
c    Modifications by John Burkardt
c
c  Reference:
c
c    David Hill,
c    Algorithm AS 66:
c    The Normal Integral,
c    Applied Statistics,
c    Volume 22, Number 3, 1973, pages 424-427.
c
c  Parameters:
c
c    Input, double precision X, is one endpoint of the semi-infinite interval
c    over which the integration takes place.
c
c    Input, logical UPPER, determines whether the upper or lower
c    interval is to be integrated:
c    .TRUE.  => integrate from X to + Infinity;
c    .FALSE. => integrate from - Infinity to X.
c
c    Output, double precision ALNORM, the integral of the standard normal
c    distribution over the desired interval.
c
      implicit none

      double precision a1
      parameter ( a1 = 5.75885480458D+00 )
      double precision a2
      parameter ( a2 = 2.62433121679D+00 )
      double precision a3
      parameter ( a3 = 5.92885724438D+00 )
      double precision alnorm
      double precision b1
      parameter ( b1 = -29.8213557807D+00 )
      double precision b2
      parameter ( b2 = 48.6959930692D+00 )
      double precision c1
      parameter ( c1 = -0.000000038052D+00 )
      double precision c2
      parameter ( c2 = 0.000398064794D+00 )
      double precision c3
      parameter ( c3 = -0.151679116635D+00 )
      double precision c4
      parameter ( c4 = 4.8385912808D+00 )
      double precision c5
      parameter ( c5 = 0.742380924027D+00 )
      double precision c6
      parameter ( c6 = 3.99019417011D+00 )
      double precision con
      parameter ( con = 1.28D+00 )
      double precision d1
      parameter ( d1 = 1.00000615302D+00 )
      double precision d2
      parameter ( d2 = 1.98615381364D+00 )
      double precision d3
      parameter ( d3 = 5.29330324926D+00 )
      double precision d4
      parameter ( d4 = -15.1508972451D+00 )
      double precision d5
      parameter ( d5 = 30.789933034D+00 )
      double precision ltone
      parameter ( ltone = 7.0D+00 )
      double precision p
      parameter ( p = 0.398942280444D+00 )
      double precision q
      parameter ( q = 0.39990348504D+00 )
      double precision r
      parameter ( r = 0.398942280385D+00 )
      logical up
      logical upper
      double precision utzero
      parameter ( utzero = 18.66D+00 )
      double precision x
      double precision y
      double precision z

      up = upper
      z = x

      if ( z .lt. 0.0D+00 ) then
        up = .not. up
        z = - z
      end if

      if ( z .gt. ltone .and. 
     &  ( ( .not. up ) .or. utzero .lt. z ) ) then

        if ( up ) then
          alnorm = 0.0D+00
        else
          alnorm = 1.0D+00
        end if

        return

      end if

      y = 0.5D+00 * z * z

      if ( z .le. con ) then

        alnorm = 0.5D+00 - z * ( p - q * y
     &    / ( y + a1 + b1 
     &    / ( y + a2 + b2 
     &    / ( y + a3 ))))

      else

        alnorm = r * dexp ( - y )
     &    / ( z + c1 + d1
     &    / ( z + c2 + d2
     &    / ( z + c3 + d3
     &    / ( z + c4 + d4
     &    / ( z + c5 + d5
     &    / ( z + c6 ))))))

      end if

      if ( .not. up ) then
        alnorm = 1.0D+00 - alnorm
      end if

      return
      end
      subroutine beta_noncentral_cdf_values ( n_data, a, b, lambda, 
     &  x, fx )

c*********************************************************************72
c
cc BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
c
c  Discussion:
c
c    The values presented here are taken from the reference, where they
c    were given to a limited number of decimal places.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R Chattamvelli, R Shanmugam,
c    Algorithm AS 310:
c    Computing the Non-central Beta Distribution Function,
c    Applied Statistics,
c    Volume 46, Number 1, 1997, pages 146-156.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, B, the shape parameters.
c
c    Output, double precision LAMBDA, the noncentrality parameter.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      double precision lambda
      double precision lambda_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save lambda_vec
      save x_vec

      data a_vec /
     &   5.0D+00, 
     &   5.0D+00, 
     &   5.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  20.0D+00,
     &  10.0D+00,
     &  10.0D+00,
     &  15.0D+00,
     &  20.0D+00,
     &  20.0D+00,
     &  20.0D+00,
     &  30.0D+00,
     &  30.0D+00,
     &  10.0D+00,
     &  10.0D+00,
     &  10.0D+00,
     &  15.0D+00,
     &  10.0D+00,
     &  12.0D+00,
     &  30.0D+00,
     &  35.0D+00 /
      data b_vec /
     &   5.0D+00, 
     &   5.0D+00, 
     &   5.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  20.0D+00,
     &  20.0D+00,
     &  10.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  30.0D+00,
     &  50.0D+00,
     &  20.0D+00,
     &  40.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  30.0D+00,
     &  20.0D+00,
     &   5.0D+00,
     &  17.0D+00,
     &  30.0D+00,
     &  30.0D+00 /
      data fx_vec /
     &  0.4563021D+00,
     &  0.1041337D+00,
     &  0.6022353D+00,
     &  0.9187770D+00,
     &  0.6008106D+00,
     &  0.0902850D+00,
     &  0.9998655D+00,
     &  0.9925997D+00,
     &  0.9641112D+00,
     &  0.9376626573D+00,
     &  0.7306817858D+00,
     &  0.1604256918D+00,
     &  0.1867485313D+00,
     &  0.6559386874D+00,
     &  0.9796881486D+00,
     &  0.1162386423D+00,
     &  0.9930430054D+00,
     &  0.0506899273D+00,
     &  0.1030959706D+00,
     &  0.9978417832D+00,
     &  0.2555552369D+00,
     &  0.0668307064D+00,
     &  0.0113601067D+00,
     &  0.7813366615D+00,
     &  0.8867126477D+00 /
      data lambda_vec /
     &   54.0D+00, 
     &  140.0D+00,
     &  170.0D+00,
     &   54.0D+00,
     &  140.0D+00,
     &  250.0D+00,
     &   54.0D+00,
     &  140.0D+00,
     &  250.0D+00,
     &  150.0D+00,
     &  120.0D+00,
     &   80.0D+00,
     &  110.0D+00,
     &   65.0D+00,
     &  130.0D+00,
     &   80.0D+00,
     &  130.0D+00,
     &   20.0D+00,
     &   54.0D+00,
     &   80.0D+00,
     &  120.0D+00,
     &   55.0D+00,
     &   64.0D+00,
     &  140.0D+00,
     &   20.0D+00 /
      data x_vec /
     &  0.8640D+00,
     &  0.9000D+00,
     &  0.9560D+00,
     &  0.8686D+00,
     &  0.9000D+00,
     &  0.9000D+00,
     &  0.8787D+00,
     &  0.9000D+00,
     &  0.9220D+00,
     &  0.868D+00,
     &  0.900D+00,
     &  0.880D+00,
     &  0.850D+00,
     &  0.660D+00,
     &  0.720D+00,
     &  0.720D+00,
     &  0.800D+00,
     &  0.644D+00,
     &  0.700D+00,
     &  0.780D+00,
     &  0.760D+00,
     &  0.795D+00,
     &  0.560D+00,
     &  0.800D+00,
     &  0.670D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        lambda = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        lambda = lambda_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      function betain ( x, p, q, beta, ifault )

c*********************************************************************72
c
cc BETAIN computes the incomplete Beta function ratio.
c
c  Modified:
c
c    06 January 2008
c
c  Author:
c
c    KL Majumder, GP Bhattacharjee
c    Modifications by John Burkardt
c
c  Reference:
c
c    KL Majumder, GP Bhattacharjee,
c    Algorithm AS 63:
c    The incomplete Beta Integral,
c    Applied Statistics,
c    Volume 22, Number 3, 1973, pages 409-411.
c
c  Parameters:
c
c    Input, double precision X, the argument, between 0 and 1.
c
c    Input, double precision P, Q, the parameters, which
c    must be positive.
c
c    Input, double precision BETA, the logarithm of the complete
c    beta function.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    nonzero, an error occurred.
c
c    Output, double precision BETAIN, the value of the incomplete
c    Beta function ratio.
c
      implicit none

      double precision acu
      parameter ( acu = 0.1D-14 )
      double precision ai
      double precision beta
      double precision betain
      double precision cx
      integer ifault
      logical indx
      integer ns
      double precision p
      double precision pp
      double precision psq
      double precision q
      double precision qq
      double precision rx
      double precision temp
      double precision term
      double precision x
      double precision xx

      betain = x
      ifault = 0
c
c  Check the input arguments.
c
      if ( p .le. 0.0D+00 .or. q .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( x .lt. 0.0D+00 .or. 1.0D+00 .lt. x ) then
        ifault = 2
        return
      end if
c
c  Special cases.
c
      if ( x .eq. 0.0D+00 .or. x .eq. 1.0D+00 ) then
        return
      end if
c
c  Change tail if necessary and determine S.
c
      psq = p + q
      cx = 1.0D+00 - x

      if ( p .lt. psq * x ) then
        xx = cx
        cx = x
        pp = q
        qq = p
        indx = .true.
      else
        xx = x
        pp = p
        qq = q
        indx = .false.
      end if

      term = 1.0D+00
      ai = 1.0D+00
      betain = 1.0D+00
      ns = int ( qq + cx * psq )
c
c  Use Soper's reduction formula.
c
      rx = xx / cx
      temp = qq - ai
      if ( ns .eq. 0 ) then
        rx = xx
      end if

10    continue

      term = term * temp * rx / ( pp + ai )
      betain = betain + term
      temp = dabs ( term )

      if ( temp .le. acu .and. temp .le. acu * betain ) then

        betain = betain * dexp ( pp * dlog ( xx ) 
     &  + ( qq - 1.0D+00 ) * dlog ( cx ) - beta ) / pp

        if ( indx ) then
          betain = 1.0D+00 - betain
        end if

        return

      end if

      ai = ai + 1.0D+00
      ns = ns - 1

      if ( ns .ge. 0 ) then
        temp = qq - ai
        if ( ns .eq. 0 ) then
          rx = xx
        end if
      else
        temp = psq
        psq = psq + 1.0D+00
      end if

      go to 10

      return
      end
      function betanc ( x, a, b, lambda, ifault )

c*********************************************************************72
c
cc BETANC computes the tail of the noncentral Beta distribution.
c
c  Discussion:
c
c    This routine returns the cumulative probability of X for the non-central 
c    Beta distribution with parameters A, B and non-centrality LAMBDA.
c
c    Note that if LAMBDA = 0, the standard Beta distribution is defined.
c
c  Modified:
c
c    10 January 2008
c
c  Author:
c
c    Russell Lenth
c    Modifications by John Burkardt
c
c  Reference:
c
c    Russell Lenth,
c    Algorithm AS 226:
c    Computing Noncentral Beta Probabilities,
c    Applied Statistics,
c    Volume 36, Number 2, 1987, pages 241-244.
c
c  Parameters:
c
c    Input, double precision X, the value defining the cumulative
c    probability lower tail.  Normally, 0 <= X <= 1, but any value
c    is allowed.
c
c    Input, double precision A, B, the parameters of the distribution.
c    0 < A, 0 < B.
c
c    Input, double precision LAMBDA, the noncentrality parameter
c    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
c    accurate results for values of LAMBDA up to about 100.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    nonzero, an error occurred.
c
c    Output, double precision BETANC, the cumulative probability
c    of X.
c
      implicit none

      double precision a
      double precision a0
      double precision alnorm
      double precision alngam
      double precision ax
      double precision b
      double precision beta
      double precision betain
      double precision betanc
      double precision c
      double precision errbd
      double precision errmax
      parameter ( errmax = 1.0D-07 )
      double precision gx
      integer ifault
      integer itrmax
      parameter ( itrmax = 150 )
      double precision lambda
      double precision q
      double precision sumq
      double precision temp
      double precision ualpha
      parameter ( ualpha = 5.0D+00 )
      double precision x
      double precision x0
      double precision xj

      ifault = 0

      if ( lambda .lt. 0.0D+00 .or. 
     &     a .le. 0.0D+00 .or. 
     &     b .le. 0.0D+00 ) then
        ifault = 2
        betanc = -1.0D+00
        return
      end if

      if ( x .le. 0.0D+00 ) then
        betanc = 0.0D+00
        return
      end if

      if ( 1.0D+00 .le. x ) then
        betanc = 1.0D+00
        return
      end if

      c = 0.5D+00 * lambda
c
c  Initialize the series.
c
      beta = alngam ( a, ifault ) 
     &     + alngam ( b, ifault ) 
     &     - alngam ( a + b, ifault )

      temp = betain ( x, a, b, beta, ifault )

      gx = dexp ( a * dlog ( x ) + b * dlog ( 1.0D+00 - x ) 
     &  - beta - dlog ( a ) )

      q = dexp ( - c )

      xj = 0.0D+00
      ax = q * temp
      sumq = 1.0D+00 - q
      betanc = ax
c
c  Recur over subsequent terms until convergence is achieved.
c
   10 continue

      xj = xj + 1.0D+00
      temp = temp - gx
      gx = x * ( a + b + xj - 1.0D+00 ) * gx / ( a + xj )
      q = q * c / xj
      sumq = sumq - q
      ax = temp * q
      betanc = betanc + ax
c
c  Check for convergence and act accordingly.
c
      errbd = dabs ( ( temp - gx ) * sumq )

      if ( int ( xj ) .lt. itrmax .and. errbd .gt. errmax ) then
        go to 10
      end if

      if ( errbd .gt. errmax ) then
        ifault = 1
      end if

      return
      end
      function gammad ( x, p, ifault )

c*********************************************************************72
c
cc GAMMAD computes the Incomplete Gamma Integral
c
c  Auxiliary functions:
c
c    ALOGAM = logarithm of the gamma function, 
c    ALNORM = algorithm AS66
c
c  Modified:
c
c    31 March 1999
c
c  Author:
c
c    B Shea
c    Modifications by John Burkardt
c
c  Reference:
c
c    B Shea,
c    Algorithm AS 239:
c    Chi-squared and Incomplete Gamma Integral,
c    Applied Statistics,
c    Volume 37, Number 3, 1988, pages 466-473.
c
c  Parameters:
c
c    Input, double precision X, P, the parameters of the incomplete gamma ratio.
c    0 <= X, and 0 < P.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, X < 0 or P <= 0.
c
c    Output, double precision GAMMAD, the value of the incomplete Gamma integral.
c
      implicit none

      double precision a
      double precision alnorm
      double precision alngam
      double precision an
      double precision arg
      double precision b
      double precision c
      double precision elimit
      parameter ( elimit = - 88.0D+00 )
      double precision gammad
      integer ifault
      double precision oflo
      parameter ( oflo = 1.0D+37 )
      double precision p
      double precision plimit
      parameter ( plimit = 1000.0D+00 )
      double precision pn1
      double precision pn2
      double precision pn3
      double precision pn4
      double precision pn5
      double precision pn6
      double precision rn
      double precision tol
      parameter ( tol = 1.0D-14 )
      logical upper
      double precision x
      double precision xbig
      parameter ( xbig = 1.0D+08 )

      gammad = 0.0D+00
c
c  Check the input.
c
      if ( x .lt. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( p .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0

      if ( x .eq. 0.0D+00 ) then
        gammad = 0.0D+00
        return
      end if
c
c  If P is large, use a normal approximation.
c
      if ( plimit .lt. p ) then

        pn1 = 3.0D+00 * dsqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) 
     &  + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

        upper = .false.
        gammad = alnorm ( pn1, upper )
        return

      end if
c
c  If X is large set GAMMAD = 1.
c
      if ( xbig .lt. x ) then
        gammad = 1.0D+00
        return
      end if
c
c  Use Pearson's series expansion.
c  (Note that P is not large enough to force overflow in ALOGAM).
c  No need to test IFAULT on exit since P > 0.
c
      if ( x .le. 1.0D+00 .or. x .lt. p ) then

        arg = p * dlog ( x ) - x - alngam ( p + 1.0D+00, ifault )
        c = 1.0D+00
        gammad = 1.0D+00
        a = p

   40   continue

        a = a + 1.0D+00
        c = c * x / a
        gammad = gammad + c

        if ( c .gt. tol ) then
          go to 40
        end if

        arg = arg + dlog ( gammad )

        if ( elimit .le. arg ) then
          gammad = dexp ( arg )
        else
          gammad = 0.0D+00
        end if
c
c  Use a continued fraction expansion.
c
      else 

        arg = p * dlog ( x ) - x - alngam ( p, ifault )
        a = 1.0D+00 - p
        b = a + x + 1.0D+00
        c = 0.0D+00
        pn1 = 1.0D+00
        pn2 = x
        pn3 = x + 1.0D+00
        pn4 = x * b
        gammad = pn3 / pn4

   60   continue

        a = a + 1.0D+00
        b = b + 2.0D+00
        c = c + 1.0D+00
        an = a * c
        pn5 = b * pn3 - an * pn1
        pn6 = b * pn4 - an * pn2

        if ( pn6 .ne. 0.0D+00 ) then

          rn = pn5 / pn6

          if ( dabs ( gammad - rn ) .le. dmin1 ( tol, tol * rn ) ) then
            go to 80
          end if

          gammad = rn

        end if

        pn1 = pn3
        pn2 = pn4
        pn3 = pn5
        pn4 = pn6
c
c  Re-scale terms in continued fraction if terms are large.
c
        if ( dabs ( pn5 ) .ge. oflo ) then
          pn1 = pn1 / oflo
          pn2 = pn2 / oflo
          pn3 = pn3 / oflo
          pn4 = pn4 / oflo
        end if

        go to 60

   80   continue

        arg = arg + dlog ( gammad )

        if ( elimit .le. arg ) then
          gammad = 1.0D+00 - dexp ( arg )
        else
          gammad = 1.0D+00
        end if

      end if

      return
      end
      function ncbeta ( a, b, lambda, x, errmax, ifault )

c*********************************************************************72
c
cc NCBETA computes the noncentral Beta CDF.
c
c  Discussion:
c
c    Three corrections needed to be made to the text of this routine.
c    They are noted in the comments below.
c
c    Two of these corrections were errors in transcription made when
c    producing the online copy distributed by APSTAT.  
c
c    One error, an error of omission, occurred in the original printed 
c    copy of the routine, and was carried over into the online copy.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    R Chattamvelli, R Shanmugam
c    Modifications by John Burkardt
c
c  Reference:
c
c    R Chattamvelli, R Shanmugam,
c    Algorithm AS 310:
c    Computing the Non-central Beta Distribution Function,
c    Applied Statistics,
c    Volume 46, Number 1, 1997, pages 146-156.
c
c  Parameters:
c
c    Input, double precision A, B, the shape parameters.
c    0 <= A, 0 <= B.
c
c    Input, double precision LAMBDA, the noncentrality parameter.
c    0 <= LAMBDA.
c
c    Input, double precision X, the value at which the CDF is desired.
c
c    Input, double precision ERRMAX, the precision tolerance.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    1, X is 0 or 1.
c    2, X < 0 or 1 < X.
c    3, A, B or LAMBDA is less than 0.
c
c    Output, double precision NCBETA, the cumulative probability of X.
c
      implicit none

      double precision a
      double precision alngam
      external alngam
      double precision b
      double precision beta
      double precision betain
      external betain
      double precision betanc
      external betanc
      double precision c
      double precision ebd
      double precision errbd
      double precision errmax
      double precision ftemp
      double precision fx
      double precision gammad
      external gammad
      double precision gx
      integer i
      integer ifault
      integer iter1
      integer iter2
      integer iterhi
      integer iterlo
      integer j
      double precision lambda
      integer m
      double precision ncbeta
      double precision mr
      double precision psum
      double precision q
      double precision r
      double precision s
      double precision s0
      double precision s1
      double precision sum
      double precision t
      double precision t0
      double precision t1
      double precision temp
      double precision x
      integer xj

      ifault = 0
      ncbeta = x
c
c  Check parameters.
c
      if ( lambda .le. 0.0D+00 ) then
        ifault = 3
        return
      end if

      if ( a .le. 0.0D+00 ) then
        ifault = 3
        return
      end if

      if ( b .le. 0.0D+00 ) then
        ifault = 3
        return
      end if

      if ( x .le. 0.0D+00 ) then
        ncbeta = 0.0D+00
        return
      end if

      if ( 1.0D+00 .le. x ) then
        ncbeta = 1.0D+00
        return
      end if

      c = 0.5D+00 * lambda
      xj = 0.0D+00
c
c  AS 226 as it stands is sufficient in this situation.
c
      if ( lambda .lt. 54.0D+00 ) then

        ncbeta = betanc ( x, a, b, lambda, ifault )
        return

      else

        m = int ( c + 0.5D+00 )
        mr = dble ( m )
        iterlo = m - int ( 5.0D+00 * dsqrt ( mr ) )
        iterhi = m + int ( 5.0D+00 * dsqrt ( mr ) )
        t = - c + mr * dlog ( c ) - alngam ( mr + 1.0D+00, ifault )
        q = dexp ( t )
        r = q
        psum = q

        beta = alngam ( a + mr, ifault ) 
     &       + alngam ( b, ifault ) 
     &       - alngam ( a + mr + b, ifault )

        s1 = ( a + mr ) * dlog ( x ) 
     &     + b * dlog ( 1.0D+00 - x ) - dlog ( a + mr ) - beta
        gx = dexp ( s1 )
        fx = gx
        temp = betain ( x, a + mr, b, beta, ifault )
        ftemp = temp
        xj = xj + 1.0D+00
c
c  The online copy of AS 310 has "SUM = Q - TEMP" which is incorrect.
c
        sum = q * temp
        iter1 = m
c
c  The first set of iterations starts from M and goes downwards
c
   20   continue

        if ( iter1 .lt. iterlo ) then
          go to 30
        end if

        if ( q .lt. errmax ) then
          go to 30
        end if
c
c  The online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
c
        q = q * iter1 / c
        xj = xj + 1.0D+00
        gx = ( a + iter1 ) / ( x * ( a + b + iter1 - 1.0D+00 ) ) * gx
        iter1 = iter1 - 1
        temp = temp + gx
        psum = psum + q
        sum = sum + q * temp
        go to 20

   30   continue

        t0 = alngam ( a + b, ifault ) 
     &     - alngam ( a + 1.0D+00, ifault ) 
     &     - alngam ( b, ifault )

        s0 = a * dlog ( x ) + b * dlog ( 1.0D+00 - x )
c
c  Both the online copy of AS 310 and the text printed in the reference
c  did not initialize the variable S to zero, which is incorrect.
c  JVB, 12 January 2008.
c
        s = 0.0
        do i = 1, iter1
          j = i - 1
          s = s + dexp ( t0 + s0 + j * dlog ( x ) )
          t1 = dlog ( a + b + j ) - dlog ( a + 1.0D+00 + j ) + t0
          t0 = t1
        end do
c
c  Compute the first part of error bound.
c
        errbd = ( 1.0D+00 - gammad ( c, dble ( iter1 ), ifault ) ) 
     &  * ( temp + s )

        q = r
        temp = ftemp
        gx = fx
        iter2 = m

   50   continue

        ebd = errbd + ( 1.0D+00 - psum ) * temp

        if ( ebd .lt. errmax .or. iterhi .le. iter2 ) then
          go to 60
        end if

        iter2 = iter2 + 1
        xj = xj + 1.0D+00
        q = q * c / iter2
        psum = psum + q
        temp = temp - gx
        gx = x * ( a + b + iter2 - 1.0D+00 ) / ( a + iter2 ) * gx
        sum = sum + q * temp
        go to 50

   60   continue

      end if

      ncbeta = sum

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
