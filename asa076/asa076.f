      function alnorm ( x, upper )

c*********************************************************************72
c
cc ALNORM computes the cumulative density of the standard normal distribution.
c
c  Modified:
c
c    28 March 1999
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
      double precision a2
      double precision a3
      double precision b1
      double precision b2
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision c5
      double precision c6
      double precision con
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      double precision d5
      double precision ltone
      double precision p
      double precision q
      double precision r
      double precision utzero

      parameter ( a1 = 5.75885480458D+00 )
      parameter ( a2 = 2.62433121679D+00 )
      parameter ( a3 = 5.92885724438D+00 )
      parameter ( b1 = -29.8213557807D+00 )
      parameter ( b2 = 48.6959930692D+00 )
      parameter ( c1 = -0.000000038052D+00 )
      parameter ( c2 = 0.000398064794D+00 )
      parameter ( c3 = -0.151679116635D+00 )
      parameter ( c4 = 4.8385912808D+00 )
      parameter ( c5 = 0.742380924027D+00 )
      parameter ( c6 = 3.99019417011D+00 )
      parameter ( con = 1.28D+00 )
      parameter ( d1 = 1.00000615302D+00 )
      parameter ( d2 = 1.98615381364D+00 )
      parameter ( d3 = 5.29330324926D+00 )
      parameter ( d4 = -15.1508972451D+00 )
      parameter ( d5 = 30.789933034D+00 )
      parameter ( ltone = 7.0D+00 )
      parameter ( p = 0.398942280444D+00 )
      parameter ( q = 0.39990348504D+00 )
      parameter ( r = 0.398942280385D+00 )
      parameter ( utzero = 18.66D+00 )

      double precision alnorm
      logical up
      logical upper
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
     &  ( ( .not. up ) .or. z .gt. utzero ) ) then

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
      subroutine owen_values ( n_data, h, a, t )

c*********************************************************************72
c
cc OWEN_VALUES returns some values of Owen's T function.
c
c  Discussion:
c
c    Owen's T function is useful for computation of the bivariate normal
c    distribution and the distribution of a skewed normal distribution.
c
c    Although it was originally formulated in terms of the bivariate
c    normal function, the function can be defined more directly as
c
c      T(H,A) = 1 / ( 2 * pi ) *
c        Integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
c
c    In Mathematica, the function can be evaluated by:
c
c      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Mike Patefield, David Tandy,
c    Fast and Accurate Calculation of Owen's T Function,
c    Journal of Statistical Software,
c    Volume 5, Number 5, 2000, pages 1-25.
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
c    Output, double precision H, a parameter.
c
c    Output, double precision A, the upper limit of the integral.
c
c    Output, double precision T, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 28 )

      double precision a
      double precision a_vec(n_max)
      double precision h
      double precision h_vec(n_max)
      integer n_data
      double precision t
      double precision t_vec(n_max)

      save a_vec
      save h_vec
      save t_vec

      data a_vec /
     &  0.2500000000000000D+00,
     &  0.4375000000000000D+00,
     &  0.9687500000000000D+00,
     &  0.0625000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.9999975000000000D+00,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.3000000000000000D+01,
     &  0.1000000000000000D+02,
     &  0.1000000000000000D+03 /
      data h_vec /
     &  0.0625000000000000D+00,
     &  6.5000000000000000D+00,
     &  7.0000000000000000D+00,
     &  4.7812500000000000D+00,
     &  2.0000000000000000D+00,
     &  1.0000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.1250000000000000D+00,
     &  0.1250000000000000D+00,
     &  0.1250000000000000D+00,
     &  0.1250000000000000D+00,
     &  0.7812500000000000D-02,
     &  0.7812500000000000D-02,
     &  0.7812500000000000D-02,
     &  0.7812500000000000D-02,
     &  0.7812500000000000D-02,
     &  0.7812500000000000D-02 /
      data t_vec /
     &  3.8911930234701366D-02,
     &  2.0005773048508315D-11,
     &  6.3990627193898685D-13,
     &  1.0632974804687463D-07,
     &  8.6250779855215071D-03,
     &  6.6741808978228592D-02,
     &  0.4306469112078537D-01,
     &  0.6674188216570097D-01,
     &  0.7846818699308410D-01,
     &  0.7929950474887259D-01,
     &  0.6448860284750376D-01,
     &  0.1066710629614485D+00,
     &  0.1415806036539784D+00,
     &  0.1510840430760184D+00,
     &  0.7134663382271778D-01,
     &  0.1201285306350883D+00,
     &  0.1666128410939293D+00,
     &  0.1847501847929859D+00,
     &  0.7317273327500385D-01,
     &  0.1237630544953746D+00,
     &  0.1737438887583106D+00,
     &  0.1951190307092811D+00,
     &  0.7378938035365546D-01,
     &  0.1249951430754052D+00,
     &  0.1761984774738108D+00,
     &  0.1987772386442824D+00,
     &  0.2340886964802671D+00,
     &  0.2479460829231492D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        h = 0.0D+00
        a = 0.0D+00
        t = 0.0D+00
      else
        h = h_vec(n_data)
        a = a_vec(n_data)
        t = t_vec(n_data)
      end if

      return
      end
      function tfn ( x, fx )

c*********************************************************************72
c
cc TFN calculates the T-function of Owen.
c
c  Modified:
c
c    06 January 2008
c
c  Author:
c
c    JC Young, Christoph Minder
c    Modifications by John Burkardt
c
c  Reference:
c
c    MA Porter, DJ Winstanley,
c    Remark AS R30:
c    A Remark on Algorithm AS76:
c    An Integral Useful in Calculating Noncentral T and Bivariate
c    Normal Probabilities,
c    Applied Statistics,
c    Volume 28, Number 1, 1979, page 113.
c
c    JC Young, Christoph Minder,
c    Algorithm AS 76: 
c    An Algorithm Useful in Calculating Non-Central T and 
c    Bivariate Normal Distributions,
c    Applied Statistics,
c    Volume 23, Number 3, 1974, pages 455-457.
c
c  Parameters:
c
c    Input, double precision X, FX, the parameters of the function.
c
c    Output, double precision TFN, the value of the T-function.
c
      implicit none

      integer ng
      parameter ( ng = 5 )

      double precision fx
      double precision fxs
      integer i
      double precision r(ng)
      double precision r1
      double precision r2
      double precision rt
      double precision tfn
      double precision tp
      parameter ( tp = 0.159155D+00 )
      double precision tv1
      parameter ( tv1 = 1.0D-35 )
      double precision tv2
      parameter ( tv2 = 15.0D+00 )
      double precision tv3
      parameter ( tv3 = 15.0D+00 )
      double precision tv4
      parameter ( tv4 = 1.0D-05)
      double precision u(ng)
      double precision x
      double precision x1
      double precision x2
      double precision xs

      data u / 
     &  0.0744372D+00, 
     &  0.2166977D+00, 
     &  0.3397048D+00, 
     &  0.4325317D+00, 
     &  0.4869533D+00 /

      data r /
     &  0.1477621D+00, 
     &  0.1346334D+00, 
     &  0.1095432D+00, 
     &  0.0747257D+00, 
     &  0.0333357D+00 /
c
c  Test for X near zero.
c
      if ( dabs ( x ) .lt. tv1 ) then
        tfn = tp * atan ( fx )
        return
      end if
c
c  Test for large values of abs(X).
c
      if ( tv2 .lt. dabs ( x ) ) then
        tfn = 0.0D+00
        return
      end if
c
c  Test for FX near zero.
c
      if ( dabs ( fx ) .lt. tv1 ) then
        tfn = 0.0D+00
        return
      end if
c
c  Test whether abs ( FX ) is so large that it must be truncated.
c
      xs = - 0.5D+00 * x * x
      x2 = fx
      fxs = fx * fx

      if ( dlog ( 1.0D+00 + fxs ) - xs * fxs .lt. tv3 ) then
        go to 20 
      end if
c
c  Computation of truncation point by Newton iteration.
c
      x1 = 0.5D+00 * fx
      fxs = 0.25D+00 * fxs

10    continue

        rt = fxs + 1.0D+00

        x2 = x1 + ( xs * fxs + tv3 - dlog ( rt ) ) 
     &  / ( 2.0D+00 * x1 * ( 1.0D+00 / rt - xs ) )

        fxs = x2 * x2
        if ( dabs ( x2 - x1 ) .lt. tv4 ) then
          go to 20
        end if
        x1 = x2
      go to 10
c
c  Gaussian quadrature.
c
20    continue

      rt = 0.0D+00

      do i = 1, ng

        r1 = 1.0D+00 + fxs * ( 0.5D+00 + u(i) )**2
        r2 = 1.0D+00 + fxs * ( 0.5D+00 - u(i) )**2

        rt = rt + r(i) * ( dexp ( xs * r1 ) / r1 
     &  + dexp ( xs * r2 ) / r2 )

      end do

      tfn = rt * x2 * tp

      return
      end
      function tha ( h1, h2, a1, a2 )

c*********************************************************************72
c
cc THA computes Owen's T function.
c
c  Discussion:
c
c    This function computes T(H1/H2, A1/A2) for any real numbers H1, H2, 
c    A1 and A2.
c
c  Modified:
c
c    16 January 2008
c
c  Author:
c
c    JC Young, Christoph Minder
c    Modifications by John Burkardt
c
c  Reference:
c
c    Richard Boys,
c    Remark AS R80:
c    A Remark on Algorithm AS76:
c    An Integral Useful in Calculating Noncentral T and Bivariate
c    Normal Probabilities,
c    Applied Statistics,
c    Volume 38, Number 3, 1989, pages 580-582.
c
c    Youn-Min Chou,
c    Remark AS R55:
c    A Remark on Algorithm AS76:
c    An Integral Useful in Calculating Noncentral T and Bivariate
c    Normal Probabilities,
c    Applied Statistics,
c    Volume 34, Number 1, 1985, pages 100-101.
c
c    PW Goedhart, MJW Jansen,
c    Remark AS R89:
c    A Remark on Algorithm AS76:
c    An Integral Useful in Calculating Noncentral T and Bivariate
c    Normal Probabilities,
c    Applied Statistics,
c    Volume 41, Number 2, 1992, pages 496-497.
c
c    JC Young, Christoph Minder,
c    Algorithm AS 76: 
c    An Algorithm Useful in Calculating Noncentral T and 
c    Bivariate Normal Distributions,
c    Applied Statistics,
c    Volume 23, Number 3, 1974, pages 455-457.
c
c  Parameters:
c
c    Input, double precision H1, H2, A1, A2, define the arguments
c    of the T function.
c
c    Output, double precision THA, the value of Owen's T function.
c
      implicit none

      double precision a
      double precision a1
      double precision a2
      double precision absa
      double precision ah
      double precision alnorm
      double precision c1
      double precision c2
      double precision ex
      double precision g
      double precision gah
      double precision gh
      double precision h
      double precision h1
      double precision h2
      double precision lam
      double precision tfn
      double precision tha
      double precision twopi
      parameter ( twopi = 6.2831853071795864769D+00 )

      if ( h2 .eq. 0.0D+00 ) then
        tha = 0.0D+00
        return
      end if

      h = h1 / h2

      if ( a2 .eq. 0.0D+00 ) then

        g = alnorm ( h, .false. )

        if ( h .lt. 0.0D+00 ) then
          tha = g / 2.0D+00
        else
          tha = ( 1.0D+00 - g ) / 2.0D+00
        end if

        if ( a1 .lt. 0.0D+00 ) then
          tha = - tha
        end if

        return
      end if

      a = a1 / a2

      if ( dabs ( h ) .lt. 0.3D+00 .and. dabs ( a ) .gt. 7.0D+00 ) then

        lam = dabs ( a * h )
        ex = dexp ( - lam * lam / 2.0D+00 )
        g = alnorm ( lam, .false. )
        c1 = ( ex / lam + dsqrt ( twopi ) * ( g - 0.5D+00 ) ) / twopi
        c2 = ( ( lam * lam + 2.0D+00 ) * ex / lam**3 
     &  + dsqrt ( twopi ) * ( g - 0.5D+00 ) ) / ( 6.0D+00 * twopi )
        ah = dabs ( h )
        tha = 0.25D+00 - c1 * ah + c2 * ah**3
        tha = dsign ( tha, a )

      else

        absa = dabs ( a )

        if ( absa .le. 1.0D+00 ) then
          tha = tfn ( h, a )
          return
        end if

        ah = absa * h
        gh = alnorm ( h, .false. )
        gah = alnorm ( ah, .false. )
        tha = 0.5D+00 * ( gh + gah ) - gh * gah 
     &  - tfn ( ah, 1.0D+00 / absa )

        if ( a .lt. 0.0D+00 ) then
          tha = - tha
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
