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
c    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
c    Algorithm AS 245,
c    Applied Statistics,
c    Volume 38, Number 2, pages 397-402, 1989.
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
      subroutine chi_square_cdf_values ( n_data, a, x, fx )

c*********************************************************************72
c
cc CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = ChiSquareDistribution [ df ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer A, the parameter of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      integer a
      integer a_vec(n_max)
      double precision fx
      double precision fx_vec(n_max) 
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save fx_vec
      save x_vec

      data a_vec /
     &   1,  2,  1,  2, 
     &   1,  2,  3,  4, 
     &   1,  2,  3,  4, 
     &   5,  3,  3,  3, 
     &   3,  3, 10, 10, 
     &  10 /
      data fx_vec /
     &  0.7965567455405796D-01, 
     &  0.4987520807317687D-02,  
     &  0.1124629160182849D+00, 
     &  0.9950166250831946D-02, 
     &  0.4729107431344619D+00,  
     &  0.1812692469220181D+00,  
     &  0.5975750516063926D-01,  
     &  0.1752309630642177D-01,  
     &  0.6826894921370859D+00,  
     &  0.3934693402873666D+00,  
     &  0.1987480430987992D+00,  
     &  0.9020401043104986D-01,  
     &  0.3743422675270363D-01,  
     &  0.4275932955291202D+00,  
     &  0.6083748237289110D+00,  
     &  0.7385358700508894D+00,  
     &  0.8282028557032669D+00,  
     &  0.8883897749052874D+00,  
     &  0.1721156299558408D-03,  
     &  0.3659846827343712D-02,  
     &  0.1857593622214067D-01 /
      data x_vec /
     &  0.01D+00,  
     &  0.01D+00,   
     &  0.02D+00,  
     &  0.02D+00,  
     &  0.40D+00,  
     &  0.40D+00,  
     &  0.40D+00,  
     &  0.40D+00,  
     &  1.00D+00,  
     &  1.00D+00,  
     &  1.00D+00,  
     &  1.00D+00,  
     &  1.00D+00,  
     &  2.00D+00,  
     &  3.00D+00,  
     &  4.00D+00,  
     &  5.00D+00,  
     &  6.00D+00,  
     &  1.00D+00,  
     &  2.00D+00,  
     &  3.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
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
c    Chi-squared and Incomplete Gamma Integral,
c    Algorithm AS 239,
c    Applied Statistics,
c    Volume 37, Number 3, pages 466-473, 1988.
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

      double precision elimit
      double precision oflo
      double precision plimit
      double precision tol
      double precision xbig

      parameter ( elimit = - 88.0D+00 )
      parameter ( oflo = 1.0D+37 )
      parameter ( plimit = 1000.0D+00 )
      parameter ( tol = 1.0D-14 )
      parameter ( xbig = 1.0D+08 )

      double precision a
      double precision alnorm
      double precision alngam
      double precision an
      double precision arg
      double precision b
      double precision c
      double precision gammad
      integer ifault
      double precision p
      double precision pn1
      double precision pn2
      double precision pn3
      double precision pn4
      double precision pn5
      double precision pn6
      double precision rn
      logical upper
      double precision x

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
      if ( p .gt. plimit ) then

        pn1 = 3.0D+00 * dsqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) 
     &  + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

        upper = .false.
        gammad = alnorm ( pn1, upper )
        return

      end if
c
c  If X is large set GAMMAD = 1.
c
      if ( x .gt. xbig ) then
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

        if ( arg .ge. elimit ) then
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

          if ( abs ( gammad - rn ) .le. min ( tol, tol * rn ) ) then
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
        if ( abs ( pn5 ) .ge. oflo ) then
          pn1 = pn1 / oflo
          pn2 = pn2 / oflo
          pn3 = pn3 / oflo
          pn4 = pn4 / oflo
        end if

        go to 60

   80   continue

        arg = arg + dlog ( gammad )

        if ( arg .ge. elimit ) then
          gammad = 1.0D+00 - dexp ( arg )
        else
          gammad = 1.0D+00
        end if

      end if

      return
      end
      function ppchi2 ( p, v, ifault )

c*********************************************************************72
c
cc PPCHI2 evaluates the percentage points of the chi-squared PDF.
c
c  Auxiliary routines:
c
c    PPND = AS 111 or AS 241;
c    GAMMAD = AS 239.
c
c  Modified:
c
c    30 March 1999
c
c  Author:
c
c    Best, Roberts
c    Modifications by John Burkardt
c
c  Reference:
c
c    Best, Roberts,
c    The Percentage Points of the Chi-Squared Distribution,
c    Algorithm AS 91,  
c    Applied Statistics,
c    Volume 24, Number ?, pages 385-390, 1975.
c
c  Parameters:
c
c    Input, double precision P, a value of the chi-squared cumulative probability
c    density function.
c    0.000002 <= P <= 0.999998.
c
c    Input, double precision V, the parameter of the chi-squared probability density
c    function.  V > 0.
c
c    Output, integer IFAULT, error flag.
c    0, no error detected.
c    1, P < PMIN or P > PMAX.
c    2, V <= 0.0.
c    3, an error occurred in the incomplete Gamma function routine.
c    4, the maximum number of iterations were taken without convergence.
c    5, an error occurred in the log Gamma routine.
c
c    Output, double precision PPCHI2, the value of the chi-squared random deviate
c    with the property that the probability that a chi-squared random
c    deviate with parameter V is less than or equal to PPCHI2 is P.
c
      implicit none

      double precision aa
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision c5
      double precision c6
      double precision c7
      double precision c8
      double precision c9
      double precision c10
      double precision c11
      double precision c12
      double precision c13
      double precision c14
      double precision c15
      double precision c16
      double precision c17
      double precision c18
      double precision c19
      double precision c20
      double precision c21
      double precision c22
      double precision c23
      double precision c24
      double precision c25
      double precision c26
      double precision c27
      double precision c28
      double precision c29
      double precision c30
      double precision c31
      double precision c32
      double precision c33
      double precision c34
      double precision c35
      double precision c36
      double precision c37
      double precision c38
      double precision e
      integer maxit
      double precision pmax
      double precision pmin

      parameter ( aa = 0.6931471806D+00 )
      parameter ( c1 = 0.01D+00 )
      parameter ( c2 = 0.222222D+00 )
      parameter ( c3 = 0.32D+00 )
      parameter ( c4 = 0.4D+00 )
      parameter ( c5 = 1.24D+00 )
      parameter ( c6 = 2.2D+00 )
      parameter ( c7 = 4.67D+00 )
      parameter ( c8 = 6.66D+00 )
      parameter ( c9 = 6.73D+00 )
      parameter ( c10 = 13.32D+00 )
      parameter ( c11 = 60.0D+00 )
      parameter ( c12 = 70.0D+00 )
      parameter ( c13 = 84.0D+00 )
      parameter ( c14 = 105.0D+00 )
      parameter ( c15 = 120.0D+00 )
      parameter ( c16 = 127.0D+00 )
      parameter ( c17 = 140.0D+00 )
      parameter ( c18 = 175.0D+00 )
      parameter ( c19 = 210.0D+00 )
      parameter ( c20 = 252.0D+00 )
      parameter ( c21 = 264.0D+00 )
      parameter ( c22 = 294.0D+00 )
      parameter ( c23 = 346.0D+00 )
      parameter ( c24 = 420.0D+00 )
      parameter ( c25 = 462.0D+00 )
      parameter ( c26 = 606.0D+00 )
      parameter ( c27 = 672.0D+00 )
      parameter ( c28 = 707.0D+00 )
      parameter ( c29 = 735.0D+00 )
      parameter ( c30 = 889.0D+00 )
      parameter ( c31 = 932.0D+00 )
      parameter ( c32 = 966.0D+00 )
      parameter ( c33 = 1141.0D+00 )
      parameter ( c34 = 1182.0D+00 )
      parameter ( c35 = 1278.0D+00 )
      parameter ( c36 = 1740.0D+00 )
      parameter ( c37 = 2520.0D+00 )
      parameter ( c38 = 5040.0D+00 )
      parameter ( e = 0.0000005D+00 )
      parameter ( maxit = 20 )
      parameter ( pmax = 0.999998D+00 )
      parameter ( pmin = 0.000002D+00 )

      double precision a
      double precision alngam
      double precision b
      double precision c
      double precision ch
      double precision g
      double precision gammad
      integer i
      integer ifault
      integer ifault2
      double precision p
      double precision p1
      double precision p2
      double precision ppchi2
      double precision ppnd
      double precision q
      double precision s1
      double precision s2
      double precision s3
      double precision s4
      double precision s5
      double precision s6
      double precision t
      double precision v
      double precision x
      double precision xx

      ifault2 = 0
c
c  Check the input.
c
      if ( p .lt. pmin .or. p .gt. pmax ) then
        ifault = 1
        ppchi2 = - 1.0D+00
        return
      end if

      if ( v .le. 0.0D+00 ) then
        ifault = 2
        ppchi2 = - 1.0D+00
        return
      end if

      ifault = 0
      xx = 0.5D+00 * v
      c = xx - 1.0D+00
c
c  Compute Log ( Gamma ( V/2 ) ).
c
      g = alngam ( v / 2.0D+00, ifault )

      if ( ifault .ne. 0 ) then
        ifault = 5
        return
      end if
c
c  Starting approximation for small chi-squared.
c
      if ( v .lt. - c5 * dlog ( p ) ) then

        ch = ( p * xx * dexp ( g + xx * aa ) )**( 1.0D+00 / xx )

        if ( ch .lt. e ) then
          ifault = 0
          ppchi2 = ch
          return
        end if
c
c  Starting approximation for V less than or equal to 0.32.
c
      else if ( v .le. c3 ) then

        ch = c4
        a = log ( 1.0D+00 - p )

10      continue

        q = ch
        p1 = 1.0D+00 + ch * ( c7 + ch )
        p2 = ch * ( c9 + ch * ( c8 + ch ) )

        t = - 0.5D+00 + ( c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 +
     &    3.0D+00 * ch ) ) / p2

        ch = ch - ( 1.0D+00 - exp ( a + g + 0.5D+00 * ch + c * aa ) *
     &    p2 / p1 ) / t

        if ( dabs ( q / ch - 1.0D+00 ) .gt. c1 ) then
          go to 10
        end if
c
c  Call to algorithm AS 111.
c  Note that P has been tested above.
c  AS 241 could be used as an alternative.
c
      else

        x = ppnd ( p, ifault2 )
c
c  Starting approximation using Wilson and Hilferty estimate.
c
        p1 = c2 / v
        ch = v * ( x * dsqrt ( p1 ) + 1.0D+00 - p1 )**3
c
c  Starting approximation for P tending to 1.
c
        if ( ch .gt. c6 * v + 6.0 ) then
          ch = - 2.0D+00 * ( log ( 1.0 - p ) 
     &    - c * dlog ( 0.5D+00 * ch ) + g )
        end if

      end if
c
c  Call to algorithm AS 239 and calculation of seven term Taylor series.
c
      do i = 1, maxit

        q = ch
        p1 = 0.5D+00 * ch
        p2 = p - gammad ( p1, xx, ifault2 )

        if ( ifault2 .ne. 0 ) then
          ppchi2 = - 1.0D+00
          ifault = 3
          return
        end if

        t = p2 * dexp ( xx * aa + g + p1 - c * dlog ( ch ) )
        b = t / ch
        a = 0.5D+00 * t - b * c

        s1 = 
     &      ( c19 + a 
     &    * ( c17 + a 
     &    * ( c14 + a 
     &    * ( c13 + a 
     &    * ( c12 + a 
     &    *   c11 ) ) ) ) ) / c24

        s2 = 
     &      ( c24 + a 
     &    * ( c29 + a 
     &    * ( c32 + a 
     &    * ( c33 + a
     &    *   c35 ) ) ) ) / c37

        s3 = 
     &      ( c19 + a 
     &    * ( c25 + a 
     &    * ( c28 + a * 
     &        c31 ) ) ) / c37

        s4 = 
     &      ( c20 + a 
     &    * ( c27 + a 
     &    *   c34 ) + c 
     &    * ( c22 + a 
     &    * ( c30 + a 
     &    *   c36 ) ) ) / c38

        s5 = ( c13 + c21 * a + c * ( c18 + c26 * a ) ) / c37

        s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38

        ch = ch + t * ( 1.0D+00 + 0.5D+00 * t * s1 - b * c 
     &    * ( s1 - b 
     &    * ( s2 - b 
     &    * ( s3 - b 
     &    * ( s4 - b 
     &    * ( s5 - b 
     &    *   s6 ) ) ) ) ) )

        if ( dabs ( q / ch - 1.0D+00 ) .gt. e ) then
          ifault = 0
          ppchi2 = ch
          return
        end if

      end do

      ifault = 4
      ppchi2 = ch

      return
      end
      function ppnd ( p, ifault )

c*********************************************************************72
c
cc PPND produces the normal deviate value corresponding to lower tail area = P.
c
c  Modified:
c
c    21 January 2008
c
c  Author:
c
c    J Beasley, S Springer
c    Modifications by John Burkardt
c
c  Reference:
c
c    J Beasley, S Springer,
c    Algorithm AS 111:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 26, Number 1, 1977, pages 118-121.
c
c  Parameters:
c
c    Input, double precision P, the value of the cumulative probability 
c    densitity function.  0 < P < 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, P <= 0 or P >= 1.  PPND is returned as 0.
c
c    Output, double precision PPND, the normal deviate value with the property that
c    the probability of a standard normal deviate being less than or
c    equal to PPND is P.
c
      implicit none

      double precision a0
      parameter ( a0 = 2.50662823884D+00 )
      double precision a1
      parameter ( a1 = -18.61500062529D+00 )
      double precision a2
      parameter ( a2 = 41.39119773534D+00 )
      double precision a3
      parameter ( a3 = -25.44106049637D+00 )
      double precision b1
      parameter ( b1 = -8.47351093090D+00 )
      double precision b2
      parameter ( b2 = 23.08336743743D+00 )
      double precision b3
      parameter ( b3 = -21.06224101826D+00 )
      double precision b4
      parameter ( b4 = 3.13082909833D+00 )
      double precision c0
      parameter ( c0 = -2.78718931138D+00 )
      double precision c1
      parameter ( c1 = -2.29796479134D+00 )
      double precision c2
      parameter ( c2 = 4.85014127135D+00 )
      double precision c3
      parameter ( c3 = 2.32121276858D+00 )
      double precision d1
      parameter ( d1 = 3.54388924762D+00 )
      double precision d2
      parameter ( d2 = 1.63706781897D+00 )
      integer ifault
      double precision p
      double precision ppnd
      double precision r
      double precision split
      parameter ( split = 0.42D+00 )
      double precision value

      ifault = 0
c
c  0.08 < P < 0.92
c
      if ( dabs ( p - 0.5D+00 ) .le. split ) then

        r = ( p - 0.5D+00 ) * ( p - 0.5D+00 )

        value = ( p - 0.5D+00 ) * ( ( ( 
     &      a3   * r 
     &    + a2 ) * r 
     &    + a1 ) * r 
     &    + a0 ) / ( ( ( (
     &      b4   * r 
     &    + b3 ) * r 
     &    + b2 ) * r
     &    + b1 ) * r 
     &    + 1.0D+00 )
c
c  P < 0.08 or P > 0.92, 
c  R = min ( P, 1-P )
c
      else if ( 0.0D+00 .lt. p .and. p .lt. 1.0D+00 ) then

        if ( p .gt. 0.5D+00 ) then
          r = dsqrt ( - dlog ( 1.0D+00 - p ) )
        else
          r = dsqrt ( - dlog ( p ) )
        end if

        value = ( ( (
     &      c3   * r 
     &    + c2 ) * r 
     &    + c1 ) * r 
     &    + c0 ) / ( ( 
     &      d2   * r 
     &    + d1 ) * r 
     &    + 1.0D+00 )

        if ( p .lt. 0.5D+00 ) then
          value = - value
        end if
c
c  P <= 0.0 or 1.0 <= P
c
      else

        ifault = 1
        value = 0.0D+00

      end if

      ppnd = value

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
