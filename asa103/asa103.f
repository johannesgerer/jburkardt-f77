      function digama ( x, ifault )

c*********************************************************************72
c
cc DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    Jose Bernardo
c    Modifications by John Burkardt
c
c  Reference:
c
c    Jose Bernardo,
c    Algorithm AS 103:
c    Psi ( Digamma ) Function,
c    Applied Statistics,
c    Volume 25, Number 3, 1976, pages 315-317.
c
c  Parameters:
c
c    Input, double precision X, the argument of the digamma function.
c    0 < X.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, X <= 0.
c
c    Output, double precision DIGAMA, the value of the digamma function at X.
c
      implicit none


      double precision c
      parameter ( c = 8.5D+00 )
      double precision d1
      parameter ( d1 = -0.5772156649D+00 )
      double precision digama
      integer ifault
      double precision r
      double precision s
      parameter ( s = 0.00001D+00 )
      double precision s3
      parameter ( s3 = 0.08333333333D+00 )
      double precision s4
      parameter ( s4 = 0.0083333333333D+00 )
      double precision s5
      parameter ( s5 = 0.003968253968D+00 )
      double precision x
      double precision y
c
c  Check the input.
c
      if ( x .le. 0.0D+00 ) then
        digama = 0.0D+00
        ifault = 1
        return
      end if
c
c  Initialize.
c
      ifault = 0
      y = x
      digama = 0.0D+00
c
c  Use approximation if argument <= S.
c
      if ( y .le. s ) then
        digama = d1 - 1.0D+00 / y
        return
      end if
c
c  Reduce to DIGAMA(X + N) where (X + N) >= C.
c
10    continue

      if ( y .lt. c ) then
        digama = digama - 1.0D+00 / y
        y = y + 1.0D+00
        go to 10
      end if
c
c  Use Stirling's (actually de Moivre's) expansion if argument > C.
c
      r = 1.0D+00 / y
      digama = digama + log ( y ) - 0.5D+00 * r
      r = r * r
      digama = digama - r * ( s3 - r * ( s4 - r * s5 ) )

      return
      end
      subroutine psi_values ( n_data, x, fx )

c*********************************************************************72
c
cc PSI_VALUES returns some values of the Psi or Digamma function for testing.
c
c  Discussion:
c
c    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
c
c    PSI(1) = - Euler's constant.
c
c    PSI(X+1) = PSI(X) + 1 / X.
c
c  Modified:
c
c    31 March 2007
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
      double precision fxvec ( n_max )
      integer n_data
      double precision x
      double precision xvec ( n_max )

      data fxvec /
     &  -0.5772156649015329D+00, 
     &  -0.4237549404110768D+00, 
     &  -0.2890398965921883D+00, 
     &  -0.1691908888667997D+00, 
     &  -0.6138454458511615D-01, 
     &   0.3648997397857652D-01, 
     &   0.1260474527734763D+00, 
     &   0.2085478748734940D+00, 
     &   0.2849914332938615D+00, 
     &   0.3561841611640597D+00, 
     &   0.4227843350984671D+00 /

      data xvec /
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
        x = xvec(n_data)
        fx = fxvec(n_data)
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

