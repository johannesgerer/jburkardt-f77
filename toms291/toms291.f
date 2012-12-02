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
