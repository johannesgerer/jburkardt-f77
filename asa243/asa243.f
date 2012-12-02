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
c    KL Majumder, GP Bhattacharjee,
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
c  Use the Soper reduction formula.
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
      subroutine student_noncentral_cdf_values ( n_data, df, lambda, 
     &  x, fx )

c*********************************************************************72
c
cc STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = NoncentralStudentTDistribution [ df, lambda ]
c      CDF [ dist, x ]
c
c    Mathematica seems to have some difficulty computing this function
c    to the desired number of digits.
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer DF, double precision LAMBDA, the parameters of the
c    function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 30 )

      integer df
      integer df_vec(n_max) 
      double precision fx
      double precision fx_vec(n_max) 
      double precision lambda
      double precision lambda_vec(n_max) 
      integer n_data
      double precision x
      double precision x_vec(n_max) 

      save df_vec
      save fx_vec
      save lambda_vec
      save x_vec

      data df_vec /
     &   1,  2,  3, 
     &   1,  2,  3, 
     &   1,  2,  3, 
     &   1,  2,  3, 
     &   1,  2,  3, 
     &  15, 20, 25, 
     &   1,  2,  3, 
     &  10, 10, 10, 
     &  10, 10, 10, 
     &  10, 10, 10 /
      data fx_vec /
     &  0.8975836176504333D+00, 
     &  0.9522670169D+00, 
     &  0.9711655571887813D+00, 
     &  0.8231218864D+00, 
     &  0.9049021510D+00, 
     &  0.9363471834D+00, 
     &  0.7301025986D+00, 
     &  0.8335594263D+00, 
     &  0.8774010255D+00, 
     &  0.5248571617D+00, 
     &  0.6293856597D+00, 
     &  0.6800271741D+00, 
     &  0.20590131975D+00, 
     &  0.2112148916D+00, 
     &  0.2074730718D+00, 
     &  0.9981130072D+00, 
     &  0.9994873850D+00, 
     &  0.9998391562D+00, 
     &  0.168610566972D+00, 
     &  0.16967950985D+00, 
     &  0.1701041003D+00, 
     &  0.9247683363D+00, 
     &  0.7483139269D+00, 
     &  0.4659802096D+00, 
     &  0.9761872541D+00, 
     &  0.8979689357D+00, 
     &  0.7181904627D+00, 
     &  0.9923658945D+00, 
     &  0.9610341649D+00, 
     &  0.8688007350D+00 /
      data lambda_vec /
     &  0.0D+00, 
     &  0.0D+00, 
     &  0.0D+00, 
     &  0.5D+00, 
     &  0.5D+00, 
     &  0.5D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  4.0D+00, 
     &  4.0D+00, 
     &  4.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  2.0D+00, 
     &  3.0D+00, 
     &  4.0D+00, 
     &  2.0D+00, 
     &  3.0D+00, 
     &  4.0D+00, 
     &  2.0D+00, 
     &  3.0D+00, 
     &  4.0D+00 /
      data x_vec /
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &   3.00D+00, 
     &  15.00D+00, 
     &  15.00D+00, 
     &  15.00D+00, 
     &   0.05D+00, 
     &   0.05D+00, 
     &   0.05D+00, 
     &   4.00D+00, 
     &   4.00D+00, 
     &   4.00D+00, 
     &   5.00D+00, 
     &   5.00D+00, 
     &   5.00D+00, 
     &   6.00D+00, 
     &   6.00D+00, 
     &   6.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        df = 0
        lambda = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        df = df_vec(n_data)
        lambda = lambda_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

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
      function tnc ( t, df, delta, ifault )

c*********************************************************************72
c
cc TNC computes the tail of the noncentral T distribution.
c
c  Discussion:
c
c    This routine computes the cumulative probability at T of the 
c    non-central T-distribution with DF degrees of freedom (which may 
c    be fractional) and non-centrality parameter DELTA.
c
c  Modified:
c
c    09 January 2008
c
c  Author:
c
c     Russell Lenth
c     Modifications by John Burkardt
c
c  Reference:
c
c    Russell Lenth,
c    Algorithm AS 243:
c    Cumulative Distribution Function of the Non-Central T Distribution,
c    Applied Statistics,
c    Volume 38, Number 1, 1989, pages 185-189.
c
c    William Guenther,
c    Evaluation of probabilities for the noncentral distributions and 
c    difference of two T-variables with a desk calculator,
c    Journal of Statistical Computation and Simulation, 
c    Volume 6, Number 3-4, 1978, pages 199-206.
c
c  Parameters:
c
c    Input, double precision T, the value whose cumulative density
c    is desired.
c
c    Input, double precision DF, the number of degrees of freedom.
c
c    Input, double precision DELTA, the noncentrality parameter.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    nonzero, an error occcurred.
c
c    Output, double precision TNC, the tail of the noncentral
c    T distribution.
c
      implicit none

      double precision a
      double precision albeta
      double precision alngam
      double precision alnorm
      double precision alnrpi
      parameter ( alnrpi = 0.57236494292470008707D+00 )
      double precision b
      double precision betain
      double precision del
      double precision delta
      double precision df
      double precision en
      double precision errbd
      double precision errmax
      parameter ( errmax = 1.0D-10 )
      double precision geven
      double precision godd
      double precision half
      integer ifault
      integer itrmax
      parameter ( itrmax = 100 )
      double precision lambda
      logical negdel
      double precision one
      double precision p
      double precision q
      double precision r2pi
      parameter ( r2pi = 0.79788456080286535588D+00 )
      double precision rxb
      double precision s
      double precision t
      double precision tnc
      double precision tt
      double precision two
      double precision x
      double precision xeven
      double precision xodd
      double precision zero

      tnc = 0.0D+00

      if ( df .le. 0.0D+00 ) then
        ifault = 2
        return
      end if

      ifault = 0

      tt = t
      del = delta
      negdel = .false.

      if ( t .lt. 0.0D+00 ) then
        negdel = .true.
        tt = - tt
        del = - del
      end if
c
c  Initialize twin series.
c
      en = 1.0D+00
      x = t * t / ( t * t + df )

      if ( x .le. 0.0D+00 ) then
        go to 20
      end if

      lambda = del * del
      p = 0.5D+00 * dexp ( - 0.5D+00 * lambda )
      q = r2pi * p * del
      s = 0.5D+00 - p
      a = 0.5D+00
      b = 0.5D+00 * df
      rxb = ( 1.0D+00 - x )**b
      albeta = alnrpi + alngam ( b, ifault ) - alngam ( a + b, ifault )
      xodd = betain ( x, a, b, albeta, ifault )
      godd = 2.0D+00 * rxb * dexp ( a * dlog ( x ) - albeta )
      xeven = 1.0D+00 - rxb
      geven = b * x * rxb
      tnc = p * xodd + q * xeven
c
c  Repeat until convergence.
c
   10 continue

      a = a + 1.0D+00
      xodd = xodd - godd
      xeven = xeven - geven
      godd = godd * x * ( a + b - 1.0D+00 ) / a
      geven = geven * x * ( a + b - 0.5D+00 ) / ( a + 0.5D+00 )
      p = p * lambda / ( 2.0D+00 * en )
      q = q * lambda / ( 2.0D+00 * en + 1.0D+00 )
      s = s - p
      en = en + 1.0D+00
      tnc = tnc + p * xodd + q * xeven
      errbd = 2.0D+00 * s * ( xodd - godd )

      if ( errmax .lt. errbd .and. en .le. itrmax ) then
        go to 10
      end if

   20 continue

      if ( itrmax .lt. en ) then
        ifault = 1
        return
      end if

      ifault = 0
      tnc = tnc + alnorm ( del, .true. )

      if ( negdel ) then
        tnc = 1.0D+00 - tnc
      end if

      return
      end
