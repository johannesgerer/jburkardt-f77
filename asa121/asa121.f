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
      function trigam ( x, ifault )

c*********************************************************************72
c
cc TRIGAM calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
c
c  Modified:
c
c    28 March 1999
c
c  Author:
c
c    BE Schneider
c    Modifications by John Burkardt
c
c  Reference:
c
c    BE Schneider,
c    Algorithm AS 121:
c    Trigamma Function,
c    Applied Statistics, 
c    Volume 27, Number 1, pages 97-99, 1978.
c
c  Parameters:
c
c    Input, double precision X, the argument of the trigamma function.
c    0 < X.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, X <= 0.
c
c    Output, double precision TRIGAM, the value of the trigamma function at X.
c
      implicit none

      double precision a
      parameter ( a = 0.0001D+00 )
      double precision b
      parameter ( b = 5.0D+00 )
      double precision b2
      parameter ( b2 =  0.1666666667D+00 )
      double precision b4
      parameter ( b4 = -0.03333333333D+00 )
      double precision b6
      parameter ( b6 =  0.02380952381D+00 )
      double precision b8
      parameter ( b8 = -0.03333333333D+00 )
      integer ifault
      double precision trigam
      double precision x
      double precision y
      double precision z
c
c  Check the input.
c
      if ( x .le. 0.0D+00 ) then
        ifault = 1
        trigam = 0.0D+00
        return
      end if

      ifault = 0
      z = x
c
c  Use small value approximation if X <= A.
c
      if ( x .le. a ) then
        trigam = 1.0D+00 / x / x
        return
      end if
c
c  Increase argument to ( X + I ) >= B.
c
      trigam = 0.0D+00

10    continue

      if ( z .lt. b ) then
        trigam = trigam + 1.0D+00 / z / z
        z = z + 1.0D+00
        go to 10
      end if
c
c  Apply asymptotic formula if argument is B or greater.
c
      y = 1.0D+00 / z / z

      trigam = trigam + 0.5D+00 * 
     &    y + ( 1.0D+00
     &  + y * ( b2 
     &  + y * ( b4 
     &  + y * ( b6 
     &  + y *   b8 )))) / z

      return
      end
      subroutine trigamma_values ( n_data, x, fx )

c*********************************************************************72
c
cc TRIGAMMA_VALUES returns some values of the TriGamma function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      PolyGamma[1,x]
c
c    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
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
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision fx
      double precision fx_vec(n_max) 
      integer n_data
      double precision x
      double precision x_vec(n_max) 

      save fx_vec
      save x_vec

      data fx_vec /
     &   0.1644934066848226D+01, 
     &   0.1433299150792759D+01, 
     &   0.1267377205423779D+01, 
     &   0.1134253434996619D+01, 
     &   0.1025356590529597D+01, 
     &   0.9348022005446793D+00, 
     &   0.8584318931245799D+00, 
     &   0.7932328301639984D+00, 
     &   0.7369741375017002D+00, 
     &   0.6879720582426356D+00, 
     &   0.6449340668482264D+00 /
      data x_vec /
     &  1.0D+00, 
     &  1.1D+00, 
     &  1.2D+00, 
     &  1.3D+00, 
     &  1.4D+00, 
     &  1.5D+00, 
     &  1.6D+00, 
     &  1.7D+00, 
     &  1.8D+00, 
     &  1.9D+00, 
     &  2.0D+00 /

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
