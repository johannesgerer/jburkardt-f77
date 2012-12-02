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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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
      function ppnd7 ( p, ifault )

c*********************************************************************72
c
cc PPND7 produces the normal deviate value corresponding to lower tail area = P.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**7.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    Michael Wichura
c
c  Reference:
c
c    Michael Wichura,
c    Algorithm AS 241:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 37, Number 3, 1988, pages 477-484.
c
c  Parameters:
c
c    Input, real P, the value of the cumulative probability densitity function.
c    0 < P < 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, P <= 0 or P >= 1.
c
c    Output, real PPND7, the normal deviate value with the property that
c    the probability of a standard normal deviate being less than or
c    equal to PPND7 is P.
c
      implicit none

      real a0
      real a1
      real a2
      real a3
      real b1
      real b2
      real b3
      real c0
      real c1
      real c2
      real c3
      real const1
      real const2
      real d1
      real d2
      real e0
      real e1
      real e2
      real e3
      real f1
      real f2
      integer ifault
      real p
      real ppnd7
      real q
      real r
      real split1
      real split2

      parameter ( a0 = 3.3871327179E+00 )
      parameter ( a1 = 50.434271938E+00 )
      parameter ( a2 = 159.29113202E+00 )
      parameter ( a3 = 59.109374720E+00 )
      parameter ( b1 = 17.895169469E+00 )
      parameter ( b2 = 78.757757664E+00 )
      parameter ( b3 = 67.187563600E+00 )
      parameter ( c0 = 1.4234372777E+00 )
      parameter ( c1 = 2.7568153900E+00 )
      parameter ( c2 = 1.3067284816E+00 )
      parameter ( c3 = 0.17023821103E+00 )
      parameter ( const1 = 0.180625E+00 )
      parameter ( const2 = 1.6E+00 )
      parameter ( d1 = 0.73700164250E+00 )
      parameter ( d2 = 0.12021132975E+00 )
      parameter ( e0 = 6.6579051150E+00 )
      parameter ( e1 = 3.0812263860E+00 )
      parameter ( e2 = 0.42868294337E+00 )
      parameter ( e3 = 0.017337203997E+00 )
      parameter ( f1 = 0.24197894225E+00 )
      parameter ( f2 = 0.012258202635E+00 )
      parameter ( split1 = 0.425E+00 )
      parameter ( split2 = 5.0E+00 )

      ifault = 0
      q = p - 0.5E+00

      if ( abs ( q ) .le. split1 ) then

        r = const1 - q * q

        ppnd7 = q * (((
     &      a3   * r 
     &    + a2 ) * r 
     &    + a1 ) * r 
     &    + a0 ) / (((
     &      b3   * r 
     &    + b2 ) * r 
     &    + b1 ) * r 
     &    + 1.0E+00 )

      else

        if ( q .lt. 0.0E+00 ) then
          r = p
        else
          r = 1.0E+00 - p
        end if

        if ( r .le. 0.0E+00 ) then
          ifault = 1
          ppnd7 = 0.0E+00
          return
        end if

        r = sqrt ( - alog ( r ) )

        if ( r .le. split2 ) then

          r = r - const2

          ppnd7 = (((
     &      c3   * r 
     &    + c2 ) * r 
     &    + c1 ) * r
     &    + c0 ) / ((
     &      d2   * r 
     &    + d1 ) * r 
     &    + 1.0E+00 )

        else

          r = r - split2

          ppnd7 = (((
     &      e3   * r 
     &    + e2 ) * r 
     &    + e1 ) * r 
     &    + e0 ) / ((
     &      f2   * r 
     &    + f1 ) * r 
     &    + 1.0E+00 )

        end if

        if ( q .lt. 0.0E+00 ) then
          ppnd7 = - ppnd7
        end if

      end if

      return
      end
      function ppnd16 ( p, ifault )

c*********************************************************************72
c
cc PPND16 produces the normal deviate value corresponding to lower tail area = P.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**16.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    Michael Wichura
c
c  Reference:
c
c    Michael Wichura,
c    Algorithm AS 241:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 37, Number 3, 1988, pages 477-484.
c
c  Parameters:
c
c    Input, double precision P, the value of the cumulative probability 
c    densitity function.  0 < P < 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, P <= 0 or P >= 1.
c
c    Output, double precision PPND16, the normal deviate value with the 
c    property that the probability of a standard normal deviate being 
c    less than or equal to PPND16 is P.
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
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision c5
      double precision c6
      double precision c7
      double precision const1
      double precision const2
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      double precision d5
      double precision d6
      double precision d7
      double precision e0
      double precision e1
      double precision e2
      double precision e3
      double precision e4
      double precision e5
      double precision e6
      double precision e7
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision f5
      double precision f6
      double precision f7
      integer ifault
      double precision p
      double precision ppnd16
      double precision q
      double precision r
      double precision split1
      double precision split2

      parameter ( a0 = 3.3871328727963666080D+00 )
      parameter ( a1 = 1.3314166789178437745D+02 )
      parameter ( a2 = 1.9715909503065514427D+03 )
      parameter ( a3 = 1.3731693765509461125D+04 )
      parameter ( a4 = 4.5921953931549871457D+04 )
      parameter ( a5 = 6.7265770927008700853D+04 )
      parameter ( a6 = 3.3430575583588128105D+04 )
      parameter ( a7 = 2.5090809287301226727D+03 )
      parameter ( b1 = 4.2313330701600911252D+01 )
      parameter ( b2 = 6.8718700749205790830D+02 )
      parameter ( b3 = 5.3941960214247511077D+03 )
      parameter ( b4 = 2.1213794301586595867D+04 )
      parameter ( b5 = 3.9307895800092710610D+04 )
      parameter ( b6 = 2.8729085735721942674D+04 )
      parameter ( b7 = 5.2264952788528545610D+03 )
      parameter ( c0 = 1.42343711074968357734D+00 )
      parameter ( c1 = 4.63033784615654529590D+00 )
      parameter ( c2 = 5.76949722146069140550D+00 )
      parameter ( c3 = 3.64784832476320460504D+00 )
      parameter ( c4 = 1.27045825245236838258D+00 )
      parameter ( c5 = 2.41780725177450611770D-01 )
      parameter ( c6 = 2.27238449892691845833D-02 )
      parameter ( c7 = 7.74545014278341407640D-04 )
      parameter ( const1 = 0.180625D+00 )
      parameter ( const2 = 1.6D+00 )
      parameter ( d1 = 2.05319162663775882187D+00 )
      parameter ( d2 = 1.67638483018380384940D+00 )
      parameter ( d3 = 6.89767334985100004550D-01 )
      parameter ( d4 = 1.48103976427480074590D-01 )
      parameter ( d5 = 1.51986665636164571966D-02 )
      parameter ( d6 = 5.47593808499534494600D-04 )
      parameter ( d7 = 1.05075007164441684324D-09 )
      parameter ( e0 = 6.65790464350110377720D+00 )
      parameter ( e1 = 5.46378491116411436990D+00 )
      parameter ( e2 = 1.78482653991729133580D+00 )
      parameter ( e3 = 2.96560571828504891230D-01 )
      parameter ( e4 = 2.65321895265761230930D-02 )
      parameter ( e5 = 1.24266094738807843860D-03 )
      parameter ( e6 = 2.71155556874348757815D-05 )
      parameter ( e7 = 2.01033439929228813265D-07 )
      parameter ( f1 = 5.99832206555887937690D-01 )
      parameter ( f2 = 1.36929880922735805310D-01 )
      parameter ( f3 = 1.48753612908506148525D-02 )
      parameter ( f4 = 7.86869131145613259100D-04 )
      parameter ( f5 = 1.84631831751005468180D-05 )
      parameter ( f6 = 1.42151175831644588870D-07 )
      parameter ( f7 = 2.04426310338993978564D-15 )
      parameter ( split1 = 0.425D+00 )
      parameter ( split2 = 5.D+00 )

      ifault = 0
      q = p - 0.5D+00

      if ( dabs ( q ) .le. split1 ) then

        r = const1 - q * q

        ppnd16 = q * (((((((
     &      a7   * r 
     &    + a6 ) * r 
     &    + a5 ) * r 
     &    + a4 ) * r 
     &    + a3 ) * r 
     &    + a2 ) * r 
     &    + a1 ) * r 
     &    + a0 ) / (((((((
     &      b7   * r 
     &    + b6 ) * r 
     &    + b5 ) * r 
     &    + b4 ) * r 
     &    + b3 ) * r 
     &    + b2 ) * r 
     &    + b1 ) * r 
     &    + 1.0D+00 )

      else

        if ( q .lt. 0.0D+00 ) then
          r = p
        else
          r = 1.0D+00 - p
        end if

        if ( r .le. 0.0D+00 ) then
          ifault = 1
          ppnd16 = 0.0D+00
          return
        end if

        r = dsqrt ( - dlog ( r ) )

        if ( r .le. split2 ) then

          r = r - const2

          ppnd16 = (((((((
     &        c7   * r 
     &      + c6 ) * r 
     &      + c5 ) * r 
     &      + c4 ) * r 
     &      + c3 ) * r 
     &      + c2 ) * r 
     &      + c1 ) * r 
     &      + c0 ) / (((((((
     &        d7   * r 
     &      + d6 ) * r 
     &      + d5 ) * r 
     &      + d4 ) * r 
     &      + d3 ) * r 
     &      + d2 ) * r 
     &      + d1 ) * r 
     &      + 1.0D+00 )

        else

          r = r - split2

          ppnd16 = (((((((
     &        e7   * r 
     &      + e6 ) * r 
     &      + e5 ) * r 
     &      + e4 ) * r 
     &      + e3 ) * r 
     &      + e2 ) * r 
     &      + e1 ) * r 
     &      + e0 ) / (((((((
     &        f7   * r 
     &      + f6 ) * r 
     &      + f5 ) * r 
     &      + f4 ) * r 
     &      + f3 ) * r 
     &      + f2 ) * r 
     &      + f1 ) * r 
     &      + 1.0D+00 )

        end if

        if ( q .lt. 0.0D+00 ) then
          ppnd16 = - ppnd16
        end if

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
