      function alogam ( x, ifault )

c*********************************************************************72
c
cc ALOGAM computes the logarithm of the Gamma function.
c
c  Modified:
c
c    28 March 1999
c
c  Author:
c
c    Malcolm Pike,
c    David Hill.
c
c  Reference:
c
c    Malcolm Pike, David Hill,
c    Algorithm 291: 
c    Logarithm of Gamma Function,
c    Communications of the ACM,
c    Volume 9, Number 9, September 1966, page 684.
c
c  Parameters:
c
c    Input, double precision X, the argument of the Gamma function.
c    X should be greater than 0.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, X <= 0.
c
c    Output, double precision ALOGAM, the logarithm of the Gamma function of X.
c
      implicit none

      double precision alogam
      double precision f
      integer ifault
      double precision x
      double precision y
      double precision z

      if ( x .le. 0.0D+00 ) then
        ifault = 1
        alogam = 0.0D+00
        return
      end if

      ifault = 0
      y = x

      if ( x .lt. 7.0D+00 ) then

        f = 1.0D+00
        z = y

10      continue

        if ( z .lt. 7.0D+00 ) then
          f = f * z
          z = z + 1.0D+00
          go to 10
        end if

        y = z
        f = - dlog ( f )

      else

        f = 0.0D+00

      end if

      z = 1.0D+00 / y / y
        
      alogam = f + ( y - 0.5D+00 ) * dlog ( y ) - y
     &  + 0.918938533204673D+00 +
     &  ((( 
     &  - 0.000595238095238D+00   * z 
     &  + 0.000793650793651D+00 ) * z
     &  - 0.002777777777778D+00 ) * z 
     &  + 0.083333333333333D+00 ) / y

      return
      end
      subroutine beta_cdf_values ( n_data, a, b, x, fx )

c*********************************************************************72
c
cc BETA_CDF_VALUES returns some values of the Beta CDF.
c
c  Discussion:
c
c    The incomplete Beta function may be written
c
c      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
c                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
c
c    Thus,
c
c      BETA_INC(A,B,0.0) = 0.0
c      BETA_INC(A,B,1.0) = 1.0
c
c    The incomplete Beta function is also sometimes called the
c    "modified" Beta function, or the "normalized" Beta function
c    or the Beta CDF (cumulative density function.
c
c    In Mathematica, the function can be evaluated by:
c
c      BETA[X,A,B] / BETA[A,B]
c
c    The function can also be evaluated by using the Statistics package:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = BetaDistribution [ a, b ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Karl Pearson,
c    Tables of the Incomplete Beta Function,
c    Cambridge University Press, 1968.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, B, the parameters of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 42 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save x_vec

      data a_vec /
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   1.0D+00, 
     &   1.0D+00, 
     &   1.0D+00, 
     &   1.0D+00,
     &   1.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   5.5D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  30.0D+00, 
     &  30.0D+00, 
     &  40.0D+00, 
     &  0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.1D+01, 
     &   0.2D+01, 
     &   0.3D+01, 
     &   0.4D+01, 
     &   0.5D+01 /
      data b_vec /
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   1.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   5.0D+00, 
     &   0.5D+00, 
     &   5.0D+00, 
     &   5.0D+00, 
     &  10.0D+00, 
     &   5.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  20.0D+00, 
     &  20.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &  20.0D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.2D+01, 
     &   0.3D+01, 
     &   0.4D+01, 
     &   0.5D+01, 
     &   0.2D+01, 
     &   0.2D+01, 
     &   0.2D+01, 
     &   0.2D+01 /
      data fx_vec /
     &  0.6376856085851985D-01, 
     &  0.2048327646991335D+00, 
     &  0.1000000000000000D+01, 
     &  0.0000000000000000D+00, 
     &  0.5012562893380045D-02, 
     &  0.5131670194948620D-01, 
     &  0.2928932188134525D+00, 
     &  0.5000000000000000D+00, 
     &  0.2800000000000000D-01, 
     &  0.1040000000000000D+00, 
     &  0.2160000000000000D+00, 
     &  0.3520000000000000D+00, 
     &  0.5000000000000000D+00, 
     &  0.6480000000000000D+00, 
     &  0.7840000000000000D+00, 
     &  0.8960000000000000D+00, 
     &  0.9720000000000000D+00, 
     &  0.4361908850559777D+00, 
     &  0.1516409096347099D+00, 
     &  0.8978271484375000D-01, 
     &  0.1000000000000000D+01, 
     &  0.5000000000000000D+00, 
     &  0.4598773297575791D+00, 
     &  0.2146816102371739D+00, 
     &  0.9507364826957875D+00, 
     &  0.5000000000000000D+00, 
     &  0.8979413687105918D+00, 
     &  0.2241297491808366D+00, 
     &  0.7586405487192086D+00, 
     &  0.7001783247477069D+00, 
     &  0.5131670194948620D-01, 
     &  0.1055728090000841D+00, 
     &  0.1633399734659245D+00, 
     &  0.2254033307585166D+00, 
     &  0.3600000000000000D+00, 
     &  0.4880000000000000D+00, 
     &  0.5904000000000000D+00, 
     &  0.6723200000000000D+00, 
     &  0.2160000000000000D+00, 
     &  0.8370000000000000D-01, 
     &  0.3078000000000000D-01, 
     &  0.1093500000000000D-01 /
      data x_vec /
     &  0.01D+00, 
     &  0.10D+00, 
     &  1.00D+00, 
     &  0.00D+00, 
     &  0.01D+00, 
     &  0.10D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.10D+00, 
     &  0.20D+00, 
     &  0.30D+00, 
     &  0.40D+00, 
     &  0.50D+00, 
     &  0.60D+00, 
     &  0.70D+00, 
     &  0.80D+00, 
     &  0.90D+00, 
     &  0.50D+00, 
     &  0.90D+00, 
     &  0.50D+00, 
     &  1.00D+00, 
     &  0.50D+00, 
     &  0.80D+00, 
     &  0.60D+00, 
     &  0.80D+00, 
     &  0.50D+00, 
     &  0.60D+00, 
     &  0.70D+00, 
     &  0.80D+00, 
     &  0.70D+00, 
     &  0.10D+00, 
     &  0.20D+00, 
     &  0.30D+00, 
     &  0.40D+00, 
     &  0.20D+00, 
     &  0.20D+00, 
     &  0.20D+00, 
     &  0.20D+00, 
     &  0.30D+00, 
     &  0.30D+00, 
     &  0.30D+00, 
     &  0.30D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine gamma_log_values ( n_data, x, fx )

c*********************************************************************72
c
cc GAMMA_LOG_VALUES returns some values of the Log Gamma function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Log[Gamma[x]]
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.1524063822430784D+01, 
     &  0.7966778177017837D+00, 
     &  0.3982338580692348D+00, 
     &  0.1520596783998375D+00, 
     &  0.0000000000000000D+00, 
     & -0.4987244125983972D-01, 
     & -0.8537409000331584D-01, 
     & -0.1081748095078604D+00, 
     & -0.1196129141723712D+00, 
     & -0.1207822376352452D+00, 
     & -0.1125917656967557D+00, 
     & -0.9580769740706586D-01, 
     & -0.7108387291437216D-01, 
     & -0.3898427592308333D-01, 
     &  0.00000000000000000D+00, 
     &  0.69314718055994530D+00, 
     &  0.17917594692280550D+01, 
     &  0.12801827480081469D+02, 
     &  0.39339884187199494D+02, 
     &  0.71257038967168009D+02 /
      data x_vec /
     &  0.20D+00, 
     &  0.40D+00, 
     &  0.60D+00, 
     &  0.80D+00, 
     &  1.00D+00, 
     &  1.10D+00, 
     &  1.20D+00, 
     &  1.30D+00, 
     &  1.40D+00, 
     &  1.50D+00, 
     &  1.60D+00, 
     &  1.70D+00, 
     &  1.80D+00, 
     &  1.90D+00, 
     &  2.00D+00, 
     &  3.00D+00, 
     &  4.00D+00, 
     & 10.00D+00, 
     & 20.00D+00, 
     & 30.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine mdbeta ( x, p, q, prob, ier )

c*********************************************************************72
c
cc MDBETA evaluates the incomplete beta function.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    Oliver Ludwig
c    Modifications by John Burkardt
c
c  Reference:
c
c    Oliver Ludwig,
c    Algorithm 179:
c    Incomplete Beta Ratio,
c    Communications of the ACM,
c    Volume 6, Number 6, June 1963, page 314.
c
c  Parameters:
c
c    Input, double precision X, the value to which function is to be integrated.  
c    X must be in the range [0,1] inclusive.
c
c    Input, double precision P, the first parameter.  P must be greater than 0.0.
c
c    Input, double precision Q, the second parameter.  Q must be greater than 0.0.
c
c    Output, double precision PROB.  The probability that a random variable from a
c    Beta distribution having parameters P and Q will be less than or equal to X.
c
c    Output, integer IER, error parameter.
c    0, normal exit.
c    1, X is not in the range [0,1] inclusive.
c    2, P or Q is less than or equal to 0.
c
c  Local parameters:
c
c    Local, double precision ALEPS, the logarithm of EPS1.
c
c    Local, double precision EPS, the machine precision.
c
c    Local, double precision EPS1, the smallest representable number.
c
      implicit none

      double precision aleps
      parameter ( aleps = - 179.6016D+00 )
      double precision alogam
      double precision c
      double precision cnt
      double precision d4
      double precision dp
      double precision dq
      double precision eps
      parameter ( eps = 2.2D-16 )
      double precision eps1
      parameter ( eps1 = 1.0D-78 )
      double precision finsum
      integer ib
      integer ier
      integer ifault
      double precision infsum
      integer interval
      double precision p
      double precision p1
      double precision pq
      double precision prob
      double precision ps
      double precision px
      double precision q
      double precision temp
      double precision wh
      double precision x
      double precision xb
      double precision y
c
c  Check ranges of the arguments.
c
      y = x

      if ( x .lt. 0.0D+00 .or. 1.0D+00 .lt. x ) then
        ier = 1
        return
      end if

      if ( p .le. 0.0D+00 .or. q .le. 0.0D+00 ) then
        ier = 2
        return
      end if

      ier = 0

      if ( x .le. 0.5D+00 ) then
        interval = 0
      else
        interval = 1
        temp = p
        p = q
        q = temp
        y = 1.0D+00 - y
      end if

      if ( x .eq. 0.0D+00 .or. x .eq. 1.0D+00 ) then

        prob = 0.0D+00

        if ( interval .ne. 0 ) then
          prob = 1.0D+00 - prob
          temp = p
          p = q
          q = temp
        end if

        return
      end if

      ib = q
      temp = ib
      ps = q - dble ( ib )

      if ( q .eq. temp ) then
        ps = 1.0D+00
      end if

      dp = p
      dq = q
      px = dp * dlog ( y )
      pq = alogam ( dp + dq, ifault )
      p1 = alogam ( dp, ifault )
      c = alogam ( dq, ifault )
      d4 = dlog ( dp )
      xb = px + alogam ( ps + dp, ifault ) - alogam ( ps, ifault ) - d4 - p1
c
c  Scaling
c
      ib = int ( xb / aleps )
      infsum = 0.0D+00
c
c  First term of a decreasing series will underflow.
c
      if ( ib .eq. 0 ) then

        infsum = dexp ( xb )
        cnt = infsum * dp
c
c  CNT will equal dexp ( temp ) * ( 1.d0 - ps ) * i * p * y**i / factorial ( i ).
c
        wh = 0.0D+00

80      continue
  
        wh = wh + 1.0D+00
        cnt = cnt * ( wh - ps ) * y / wh
        xb = cnt / ( dp + wh )
        infsum = infsum + xb

        if ( infsum .lt. xb / eps ) then
          go to 80
        end if

      end if

      finsum = 0.0D+00

      if ( dq .le. 1.0D+00 ) then

        prob = finsum + infsum

        if ( interval .ne. 0 ) then
          prob = 1.0D+00 - prob
          temp = p
          p = q
          q = temp
        end if

        return
      end if

      xb = px + dq * dlog ( 1.0D+00 - y ) + pq - p1 - dlog ( dq ) - c
c
c  Scaling.
c
      ib = int ( xb / aleps )

      if ( ib .lt. 0 ) then
        ib = 0
      end if

      c = 1.0D+00 / ( 1.0D+00 - y )
      cnt = dexp ( xb - dble ( ib ) * aleps )
      ps = dq
      wh = dq

100   continue

      wh = wh - 1.0D+00

      if ( wh .le. 0.0D+00 ) then

        prob = finsum + infsum

        if ( interval .ne. 0 ) then
          prob = 1.0D+00 - prob
          temp = p
          p = q
          q = temp
        end if

        return
      end if

      px = ( ps * c ) / ( dp + wh )

      if ( px .le. 1.0D+00 ) then

        if ( cnt / eps .le. finsum .or. cnt .le. eps1 / px ) then

          prob = finsum + infsum

          if ( interval .ne. 0 ) then
            prob = 1.0D+00 - prob
            temp = p
            p = q
            q = temp
          end if

          return

        end if

      end if

      cnt = cnt * px
c
c  Rescale.
c
      if ( 1.0D+00 .lt. cnt ) then
        ib = ib - 1
        cnt = cnt * eps1
      end if

      ps = wh

      if ( ib .eq. 0 ) then
        finsum = finsum + cnt
      end if

      go to 100

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
