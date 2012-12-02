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
