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
      subroutine normal_01_cdf_values ( n_data, x, fx )

c*********************************************************************72
c
cc NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = NormalDistribution [ 0, 1 ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    24 March 2007
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
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 17 )

      double precision fx
      double precision fx_vec(n_max) 
      integer n_data
      double precision x
      double precision x_vec(n_max) 

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.5000000000000000D+00, 
     &  0.5398278372770290D+00, 
     &  0.5792597094391030D+00, 
     &  0.6179114221889526D+00, 
     &  0.6554217416103242D+00, 
     &  0.6914624612740131D+00, 
     &  0.7257468822499270D+00, 
     &  0.7580363477769270D+00, 
     &  0.7881446014166033D+00, 
     &  0.8159398746532405D+00, 
     &  0.8413447460685429D+00, 
     &  0.9331927987311419D+00, 
     &  0.9772498680518208D+00, 
     &  0.9937903346742239D+00, 
     &  0.9986501019683699D+00, 
     &  0.9997673709209645D+00, 
     &  0.9999683287581669D+00 /
      data x_vec /
     &  0.0000000000000000D+00,   
     &  0.1000000000000000D+00, 
     &  0.2000000000000000D+00, 
     &  0.3000000000000000D+00, 
     &  0.4000000000000000D+00, 
     &  0.5000000000000000D+00, 
     &  0.6000000000000000D+00, 
     &  0.7000000000000000D+00, 
     &  0.8000000000000000D+00, 
     &  0.9000000000000000D+00, 
     &  0.1000000000000000D+01, 
     &  0.1500000000000000D+01, 
     &  0.2000000000000000D+01, 
     &  0.2500000000000000D+01, 
     &  0.3000000000000000D+01, 
     &  0.3500000000000000D+01, 
     &  0.4000000000000000D+01 /

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
      subroutine normp ( z, p, q, pdf )

c*********************************************************************72
c
cc NORMP computes the cumulative density of the standard normal distribution.
c
c  Discussion:
c
c    This is algorithm 5666 from Hart, et al.
c
c  Modified:
c
c    27 March 1999
c
c  Author:
c
c    Alan Miller
c    Modifications by John Burkardt
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thacher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision Z, divides the double precision line into two semi-infinite
c    intervals, over each of which the standard normal distribution
c    is to be integrated.
c
c    Output, double precision P, Q, the integrals of the standard normal
c    distribution over the intervals ( - Infinity, Z] and 
c    [Z, + Infinity ), respectively.
c
c    Output, double precision PDF, the value of the standard normal distribution
c    at Z.
c
      implicit none

      double precision cutoff
      double precision p0
      double precision p1
      double precision p2
      double precision p3
      double precision p4
      double precision p5
      double precision p6
      double precision q0
      double precision q1
      double precision q2
      double precision q3
      double precision q4
      double precision q5
      double precision q6
      double precision q7
      double precision root2pi

      parameter ( cutoff = 7.071D+00 )
      parameter ( p0 = 220.2068679123761D+00 )
      parameter ( p1 = 221.2135961699311D+00 )
      parameter ( p2 = 112.0792914978709D+00 )
      parameter ( p3 = 33.91286607838300D+00 )
      parameter ( p4 = 6.373962203531650D+00 )
      parameter ( p5 = 0.7003830644436881D+00 )
      parameter ( p6 = 0.03526249659989109D+00 )
      parameter ( q0 = 440.4137358247522D+00 )
      parameter ( q1 = 793.8265125199484D+00 )
      parameter ( q2 = 637.3336333788311D+00 )
      parameter ( q3 = 296.5642487796737D+00 )
      parameter ( q4 = 86.78073220294608D+00 )
      parameter ( q5 = 16.06417757920695D+00 )
      parameter ( q6 = 1.755667163182642D+00 )
      parameter ( q7 = 0.08838834764831844D+00 )
      parameter ( root2pi = 2.506628274631001D+00 )

      double precision expntl
      double precision p
      double precision pdf
      double precision q
      double precision z
      double precision zabs

      zabs = dabs ( z )
c
c  37 < |Z|.
c
      if ( 37.0D+00 .lt. zabs ) then

        pdf = 0.0D+00
        p = 0.0D+00
c
c  |Z| <= 37.
c
      else

        expntl = dexp ( - 0.5D+00 * zabs * zabs )
        pdf = expntl / root2pi
c
c  |Z| < CUTOFF = 10 / sqrt(2).
c
        if ( zabs .lt. cutoff ) then

          p = expntl * ((((((
     &        p6   * zabs 
     &      + p5 ) * zabs 
     &      + p4 ) * zabs 
     &      + p3 ) * zabs
     &      + p2 ) * zabs 
     &      + p1 ) * zabs 
     &      + p0 ) / (((((((
     &        q7   * zabs 
     &      + q6 ) * zabs 
     &      + q5 ) * zabs 
     &      + q4 ) * zabs 
     &      + q3 ) * zabs 
     &      + q2 ) * zabs 
     &      + q1 ) * zabs 
     &      + q0 )
c
c  CUTOF <= |Z|.
c
        else

          p = pdf / ( 
     &      zabs + 1.0D+00 / ( 
     &      zabs + 2.0D+00 / ( 
     &      zabs + 3.0D+00 / ( 
     &      zabs + 4.0D+00 / (
     &      zabs + 0.65D+00 )))))

        end if

      end if

      if ( z .lt. 0.0D+00 ) then
        q = 1.0D+00 - p
      else
        q = p
        p = 1.0D+00 - q
      end if

      return
      end
      subroutine nprob ( z, p, q, pdf )

c*********************************************************************72
c
cc NPROB computes the cumulative density of the standard normal distribution.
c
c  Modified:
c
c    31 March 1999
c
c  Author:
c
c    AG Adams
c    Modifications by John Burkardt
c
c  Reference:
c
c    AG Adams,
c    Algorithm 39:
c    Areas Under the Normal Curve,
c    Computer Journal,
c    Volume 12, Number 2, May 1969, pages 197-198.
c
c  Parameters:
c
c    Input, double precision Z, divides the double precision line into two semi-infinite
c    intervals, over each of which the standard normal distribution
c    is to be integrated.
c
c    Output, double precision P, Q, the integrals of the standard normal
c    distribution over the intervals ( - Infinity, Z] and 
c    [Z, + Infinity ), respectively.
c
c    Output, double precision PDF, the value of the standard normal distribution at Z.
c
      implicit none

      double precision a0
      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision a6
      double precision a7
      double precision b0
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision b8
      double precision b9
      double precision b10
      double precision b11

      parameter ( a0 = 0.5D+00 )
      parameter ( a1 = 0.398942280444D+00 )
      parameter ( a2 = 0.399903438504D+00 )
      parameter ( a3 = 5.75885480458D+00 )
      parameter ( a4 = 29.8213557808D+00 )
      parameter ( a5 = 2.62433121679D+00 )
      parameter ( a6 = 48.6959930692D+00 )
      parameter ( a7 = 5.92885724438D+00 )
      parameter ( b0 = 0.398942280385D+00 )
      parameter ( b1 = 0.000000038052D+00 )
      parameter ( b2 = 1.00000615302D+00 )
      parameter ( b3 = 0.000398064794D+00 )
      parameter ( b4 = 1.98615381364D+00 )
      parameter ( b5 = 0.151679116635D+00 ) 
      parameter ( b6 = 5.29330324926D+00 )
      parameter ( b7 = 4.8385912808D+00 )
      parameter ( b8 = 15.1508972451D+00 )
      parameter ( b9 = 0.742380924027D+00 )
      parameter ( b10 = 30.789933034D+00 )
      parameter ( b11 = 3.99019417011D+00 )

      double precision p
      double precision pdf
      double precision q
      double precision y
      double precision z
      double precision zabs

      zabs = dabs ( z )
c
c  |Z| between 0 and 1.28
c
      if ( abs ( z ) .le. 1.28D+00 ) then

        y = a0 * z * z
        pdf = dexp ( - y ) * b0

        q = a0 - zabs * ( a1 - a2 * y 
     &    / ( y + a3 - a4 
     &    / ( y + a5 + a6 
     &    / ( y + a7 ))))
c
c  |Z| between 1.28 and 12.7
c
      else if ( abs ( z ) .le. 12.7D+00 ) then

        y = a0 * z * z
        pdf = dexp ( - y ) * b0

        q = pdf 
     &    / ( zabs - b1 + b2
     &    / ( zabs + b3 + b4
     &    / ( zabs - b5 + b6
     &    / ( zabs + b7 - b8
     &    / ( zabs + b9 + b10
     &    / ( zabs + b11 ))))))
c
c  Z far out in tail.
c
      else

        q = 0.0D+00
        pdf = 0.0D+00

      end if

      if ( z .lt. 0.0D+00 ) then
        p = q
        q = 1.0D+00 - p
      else
        p = 1.0D+00 - q
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
