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
c    Malcolm Pike, David Hill.
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
      subroutine beta_noncentral_cdf ( a, b, lambda, x, error_max,
     &  value )

c*********************************************************************72
c
cc BETA_NONCENTRAL_CDF evaluates the noncentral Beta CDF.
c
c  Discussion:
c
c    The reference mistakenly phrases the opposite of the correct
c    stopping criterion, that is, it says:
c
c      "stop when PSUM < 1 - ERROR"
c
c    but this must be:
c
c      "stop when 1 - ERROR < PSUM."
c
c  Modified:
c
c    29 January 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Harry Posten,
c    An Effective Algorithm for the Noncentral Beta Distribution Function,
c    The American Statistician,
c    Volume 47, Number 2, May 1993, pages 129-131.
c
c  Parameters:
c
c    Input, double precision A, B, the shape parameters.
c
c    Input, double precision LAMBDA, the noncentrality parameter.
c
c    Input, double precision X, the argument of the function.
c
c    Input, double precision ERROR_MAX, the error control.
c
c    Output, double precision VALUE, the value of the noncentral Beta CDF.
c
      implicit none

      double precision a
      double precision alogam
      double precision b
      double precision beta_log
      double precision betain
      double precision bi
      double precision bj
      double precision error_max
      integer i
      integer ifault
      double precision lambda
      double precision p_sum
      double precision pb_sum
      double precision pi
      double precision pj
      double precision si
      double precision sj
      double precision value
      double precision x

      i = 0
      pi = dexp ( - lambda / 2.0D+00 )

      beta_log = alogam ( a, ifault )
     &         + alogam ( b, ifault )
     &         - alogam ( a + b, ifault )

      bi = betain ( x, a, b, beta_log, ifault )

      si = dexp (
     &    a * dlog ( x )
     &  + b * dlog ( 1.0D+00 - x )
     &  - beta_log
     &  - dlog ( a ) )

      p_sum = pi
      pb_sum = pi * bi

10    continue

        if ( 1.0D+00 - error_max .le. p_sum ) then
          go to 20
        end if

        pj = pi
        bj = bi
        sj = si

        i = i + 1
        pi = 0.5D+00 * lambda * pj / dble ( i )
        bi = bj - sj
        si = x * ( a + b + i - 1 ) * sj / ( a + i )

        p_sum = p_sum + pi
        pb_sum = pb_sum + pi * bi

      go to 10

20    continue

      value = pb_sum

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
c    29 January 2008
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
c    Input/output, integer N_DATA.  The user sets N_DATA to 0
c    before the first call.  On each call, the routine increments N_DATA by 1,
c    and returns the corresponding data; when there is no more data, the
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
     &      5.0D+00,
     &      5.0D+00,
     &      5.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     15.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     30.0D+00,
     &     30.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     15.0D+00,
     &     10.0D+00,
     &     12.0D+00,
     &     30.0D+00,
     &     35.0D+00 /
      data b_vec /
     &      5.0D+00,
     &      5.0D+00,
     &      5.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     10.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     20.0D+00,
     &     10.0D+00,
     &      5.0D+00,
     &     10.0D+00,
     &     30.0D+00,
     &     50.0D+00,
     &     20.0D+00,
     &     40.0D+00,
     &      5.0D+00,
     &     10.0D+00,
     &     30.0D+00,
     &     20.0D+00,
     &      5.0D+00,
     &     17.0D+00,
     &     30.0D+00,
     &     30.0D+00 /
      data fx_vec /
     &     0.4563021D+00,
     &     0.1041337D+00,
     &     0.6022353D+00,
     &     0.9187770D+00,
     &     0.6008106D+00,
     &     0.0902850D+00,
     &     0.9998655D+00,
     &     0.9925997D+00,
     &     0.9641112D+00,
     &     0.9376626573D+00,
     &     0.7306817858D+00,
     &     0.1604256918D+00,
     &     0.1867485313D+00,
     &     0.6559386874D+00,
     &     0.9796881486D+00,
     &     0.1162386423D+00,
     &     0.9930430054D+00,
     &     0.0506899273D+00,
     &     0.1030959706D+00,
     &     0.9978417832D+00,
     &     0.2555552369D+00,
     &     0.0668307064D+00,
     &     0.0113601067D+00,
     &     0.7813366615D+00,
     &     0.8867126477D+00 /
      data lambda_vec /
     &      54.0D+00,
     &     140.0D+00,
     &     170.0D+00,
     &      54.0D+00,
     &     140.0D+00,
     &     250.0D+00,
     &      54.0D+00,
     &     140.0D+00,
     &     250.0D+00,
     &     150.0D+00,
     &     120.0D+00,
     &      80.0D+00,
     &     110.0D+00,
     &      65.0D+00,
     &     130.0D+00,
     &      80.0D+00,
     &     130.0D+00,
     &      20.0D+00,
     &      54.0D+00,
     &      80.0D+00,
     &     120.0D+00,
     &      55.0D+00,
     &      64.0D+00,
     &     140.0D+00,
     &      20.0D+00 /
      data x_vec /
     &     0.8640D+00,
     &     0.9000D+00,
     &     0.9560D+00,
     &     0.8686D+00,
     &     0.9000D+00,
     &     0.9000D+00,
     &     0.8787D+00,
     &     0.9000D+00,
     &     0.9220D+00,
     &     0.868D+00,
     &     0.900D+00,
     &     0.880D+00,
     &     0.850D+00,
     &     0.660D+00,
     &     0.720D+00,
     &     0.720D+00,
     &     0.800D+00,
     &     0.644D+00,
     &     0.700D+00,
     &     0.780D+00,
     &     0.760D+00,
     &     0.795D+00,
     &     0.560D+00,
     &     0.800D+00,
     &     0.670D+00 /

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
