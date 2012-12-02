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
