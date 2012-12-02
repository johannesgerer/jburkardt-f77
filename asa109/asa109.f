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
      subroutine beta_inc_values ( n_data, a, b, x, fx )

c*********************************************************************72
c
cc BETA_INC_VALUES returns some values of the incomplete Beta function.
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
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Karl Pearson,
c    Tables of the Incomplete Beta Function,
c    Cambridge University Press, 1968.
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
      function xinbta ( p, q, beta, alpha, ifault )

c*********************************************************************72
c
cc XINBTA computes inverse of the incomplete Beta function.
c
c  Modified:
c
c    10 January 2008
c
c  Author:
c
c    GW Cran, KJ Martin, GE Thomas
c    Modifications by John Burkardt
c
c  Reference:
c
c    GW Cran, KJ Martin, GE Thomas,
c    Remark AS R19 and Algorithm AS 109:
c    A Remark on Algorithms AS 63: The Incomplete Beta Integral
c    and AS 64: Inverse of the Incomplete Beta Integeral,
c    Applied Statistics,
c    Volume 26, Number 1, 1977, pages 111-114.
c
c  Parameters:
c
c    Input, double precision P, Q, the parameters of the incomplete
c    Beta function.
c
c    Input, double precision BETA, the logarithm of the value of
c    the complete Beta function.
c
c    Input, double precision ALPHA, the value of the incomplete Beta
c    function.  0 <= ALPHA <= 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    nonzero, an error occurred.
c
c    Output, double precision XINBTA, the argument of the incomplete
c    Beta function which produces the value ALPHA.
c
c  Local Parameters:
c
c    Local, double precision SAE, the most negative decimal exponent
c    which does not cause an underflow.
c
      implicit none

      double precision a
      double precision acu
      double precision adj
      double precision alpha
      double precision beta
      double precision betain
      double precision fpu
      double precision g
      double precision h
      integer iex
      integer ifault
      logical indx
      double precision p
      double precision pp
      double precision prev
      double precision q
      double precision qq
      double precision r
      double precision s
      double precision sae
      parameter ( sae = -37.0D+00 )
      double precision sq
      double precision t
      double precision tx
      double precision w
      double precision xin
      double precision xinbta
      double precision y
      double precision yprev

      fpu = 10.D+00**sae

      ifault = 0
      xinbta = alpha
c
c  Test for admissibility of parameters.
c
      if ( p .le. 0.0D+00 .or. q .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( alpha .lt. 0.0D+00 .or. 1.0D+00 .lt. alpha ) then
        ifault = 2
        return
      end if

      if ( alpha .eq. 0.0D+00 .or. alpha .eq. 1.0D+00 ) then
        return
      end if
c
c  Change tail if necessary.
c
      if ( 0.5D+00 .lt. alpha ) then
        a = 1.0D+00 - alpha
        pp = q
        qq = p
        indx = .true.
      else
        a = alpha
        pp = p
        qq = q
        indx = .false.
      end if
c
c  Calculate the initial approximation.
c
      r = dsqrt ( - dlog ( a * a ) )

      y = r - ( 2.30753D+00 + 0.27061D+00 * r ) 
     &  / ( 1.0D+00 + ( 0.99229D+00 + 0.04481D+00 * r ) * r )

      if ( 1.0D+00 .lt. p .and. 1.0D+00 .lt. q ) then

        r = ( y * y - 3.0D+00 ) / 6.0D+00
        s = 1.0D+00 / ( pp + pp - 1.0D+00 )
        t = 1.0D+00 / ( qq + qq - 1.0D+00 )
        h = 2.0D+00 / ( s + t )
        w = y * dsqrt ( h + r ) / h - ( t - s ) 
     &  * ( r + 5.0D+00 / 6.0D+00 - 2.0D+00 / ( 3.0D+00 * h ) )
        xinbta = pp / ( pp + qq * dexp ( w + w ) )

      else

        r = qq + qq
        t = 1.0D+00 / ( 9.0D+00 * qq )
        t = r * ( 1.0D+00 - t + y * dsqrt ( t ) )**3

        if ( t .le. 0.0D+00 ) then
          xinbta = 1.0D+00 
     &    - dexp ( ( dlog ( ( 1.0D+00 - a ) * qq ) + beta ) / qq )
        else

          t = ( 4.0D+00 * pp + r - 2.0D+00 ) / t
  
          if ( t .le. 1.0D+00 ) then
            xinbta = dexp ( ( dlog ( a * pp ) + beta ) / pp )
          else
            xinbta = 1.0D+00 - 2.0D+00 / ( t + 1.0D+00 )
          end if

        end if

      end if
c
c  Solve for X by a modified Newton-Raphson method,
c  using the function BETAIN.
c
      r = 1.0D+00 - pp
      t = 1.0D+00 - qq
      yprev = 0.0D+00
      sq = 1.0D+00
      prev = 1.0D+00

      if ( xinbta .lt. 0.0001D+00 ) then
        xinbta = 0.0001D+00
      end if

      if ( 0.9999D+00 .lt. xinbta ) then
        xinbta = 0.9999D+00
      end if

      iex = max ( - 5.0D+00 / pp**2 - 1.0D+00 / a**0.2D+00 - 13.0D+00, 
     &  sae )

      acu = 10.0D+00**iex

    7 continue

      y = betain ( xinbta, pp, qq, beta, ifault )

      if ( ifault .ne. 0 ) then
        ifault = 3
        return
      end if

      xin = xinbta
      y = ( y - a ) * dexp ( beta + r * dlog ( xin ) 
     &  + t * dlog ( 1.0D+00 - xin ) )

      if ( y * yprev .le. 0.0D+00 ) then
        prev = max ( sq, fpu )
      end if

      g = 1.0D+00

    9 continue

      adj = g * y
      sq = adj * adj

      if ( prev .le. sq ) then
        go to 10
      end if

      tx = xinbta - adj

      if ( 0.0D+00 .le. tx .and. tx .le. 1.0D+00 ) then
        go to 11
      end if

10    continue

      g = g / 3.0D+00
      go to 9

   11 continue

      if ( prev .le. acu ) then
        go to 12
      end if

      if ( y * y .le. acu ) then
        go to 12
      end if

      if ( tx .eq. 0.0D+00 .or. tx .eq. 1.0D+00 ) then
        go to 10
      end if

      if ( tx .eq. xinbta ) then
        go to 12
      end if

      xinbta = tx
      yprev = y
      go to 7

   12 continue

      if ( indx ) then 
        xinbta = 1.0D+00 - xinbta
      end if

      return
      end
