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
      function lngamma ( z, ier )

c*********************************************************************72
c
cc LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
c
c  Discussion:
c
c    This algorithm is not part of the Applied Statistics algorithms.   
c    It is slower but gives 14 or more significant decimal digits 
c    accuracy, except around X = 1 and X = 2.   The Lanczos series from 
c    which this algorithm is derived is interesting in that it is a 
c    convergent series approximation for the gamma function, whereas 
c    the familiar series due to De Moivre (and usually wrongly called 
c    Stirling's approximation) is only an asymptotic approximation, as
c    is the true and preferable approximation due to Stirling.
c
c  Modified:
c
c    30 March 1999
c
c  Author:
c
c    Alan Miller
c    CSIRO Division of Mathematics & Statistics
c
c  Reference:
c
c    Cornelius Lanczos,
c    A precision approximation of the gamma function, 
c    SIAM Journal on Numerical Analysis, B, 
c    Volume 1, 1964, pages 86-96.
c
c  Parameters:
c
c    Input, double precision Z, the argument of the Gamma function.
c
c    Output, integer IER, error flag.
c    0, no error occurred.
c    1, Z is less than or equal to 0.
c
c    Output, double precision LNGAMMA, the logarithm of the gamma function of Z.
c
      implicit none

      double precision a(9)
      integer ier
      integer j
      double precision lngamma
      double precision lnsqrt2pi
      parameter ( lnsqrt2pi = 0.9189385332046727D+00 )
      double precision tmp
      double precision z

      data a / 
     &       0.9999999999995183D+00, 
     &     676.5203681218835D+00,
     &  - 1259.139216722289D+00, 
     &     771.3234287757674D+00,
     &   - 176.6150291498386D+00, 
     &      12.50734324009056D+00,
     &     - 0.1385710331296526D+00, 
     &       0.9934937113930748D-05,
     &       0.1659470187408462D-06 /

      if ( z .le. 0.0D+00 ) then
        ier = 1
        lngamma = 0.0D+00
        return
      end if

      ier = 0

      lngamma = 0.0D+00
      tmp = z + 7.0D+00
      do j = 9, 2, -1
        lngamma = lngamma + a(j) / tmp
        tmp = tmp - 1.0D+00
      end do

      lngamma = lngamma + a(1)
      lngamma = dlog ( lngamma ) + lnsqrt2pi - ( z + 6.5D+00 ) +
     &  ( z - 0.5D+00 ) * dlog ( z + 6.5D+00 )

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
