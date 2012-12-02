      function alnfac ( n )

c*********************************************************************72
c
cc ALNFAC computes the logarithm of the factorial of N.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial.
c
c    Output, double precision ALNFAC, the logarithm of the factorial of N.
c
      implicit none

      double precision alnfac
      double precision alngam
      integer ier
      integer n
      
      alnfac = alngam ( dble ( n + 1 ), ier )

      return
      end
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
c    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
c    Algorithm AS 245,
c    Applied Statistics,
c    Volume 38, Number 2, pages 397-402, 1989.
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
      function chyper ( point, kk, ll, mm, nn, ifault )

c*********************************************************************72
c
cc CHYPER computes point or cumulative hypergeometric probabilities.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    Richard Lund
c    Modifications by John Burkardt
c
c  Reference:
c
c    PR Freeman,
c    Algorithm AS 59:
c    Hypergeometric Probabilities,
c    Applied Statistics, 
c    Volume 22, Number 1, 1973, pages 130-133.
c
c    Richard Lund,
c    Algorithm AS 152:
c    Cumulative hypergeometric probabilities,
c    Applied Statistics, 
c    Volume 29, Number 2, 1980, pages 221-223.
c
c    BL Shea,
c    Remark AS R77:
c    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
c    Applied Statistics,
c    Volume 38, Number 1, 1989, pages 199-204.
c
c  Parameters:
c
c    Input, logical POINT, is TRUE if the point probability is desired,
c    and FALSE if the cumulative probability is desired.
c
c    Input, integer KK, the sample size.
c    0 <= KK <= MM.
c
c    Input, integer LL, the number of successes in the sample.
c    0 <= LL <= KK.
c
c    Input, integer MM, the population size that was sampled.
c    0 <= MM.
c
c    Input, integer NN, the number of "successes" in the population.
c    0 <= NN <= MM.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    nonzero, an error occurred.
c
c    Output, double precision CHYPER, the PDF (point probability) of
c    exactly LL successes out of KK samples, or the CDF (cumulative 
c    probability) of up to LL successes out of KK samples.
c
      implicit none

      double precision alnfac
      double precision alnorm
      double precision arg
      double precision chyper
      logical dir
      double precision elimit
      parameter ( elimit = - 88.0D+00 )
      integer i
      integer ifault
      integer j
      integer k
      integer kk
      integer kl
      integer l
      integer ll
      integer m
      integer mbig
      parameter ( mbig = 600 )
      double precision mean
      integer mm
      integer mnkl
      integer mvbig
      parameter ( mvbig = 1000 )
      integer n
      integer nl
      integer nn
      double precision p
      logical point
      double precision pt
      double precision rootpi
      parameter ( rootpi = 2.506628274631001D+00 )
      double precision scale
      parameter ( scale = 1.0D+35 )
      double precision sig

      ifault = 0

      k = kk + 1
      l = ll + 1
      m = mm + 1
      n = nn + 1

      dir = .true.
c
c  Check arguments are within permitted limits.
c
      chyper = 0.0D+00

      if ( n .lt. 1 .or. m .lt. n .or. k .lt. 1 .or. m .lt. k ) then
        ifault = 1
        return
      end if

      if ( l .lt. 1 .or. m - n .lt. k - l ) then
        ifault = 2
        return
      end if

      if ( .not. point ) then
        chyper = 1.0D+00
      end if

      if ( n .lt. l .or. k .lt. l ) then
        ifault = 2
        return
      end if

      ifault = 0
      chyper = 1.0D+00

      if ( k .eq. 1 .or. k .eq. m .or. n .eq. 1 .or. n .eq. m ) then
        return
      end if

      if ( .not. point .and. ll .eq. min ( kk, nn ) ) then
        return
      end if

      p = dble ( nn ) / dble ( mm - nn )

      if ( 16.0D+00 * max ( p, 1.0D+00 / p ) .lt.
     &  dble ( min ( kk, mm - kk ) ) .and.
     &  mvbig .lt. mm .and. 
     &  - 100.0D+00 .lt. elimit ) then
c
c  Use a normal approximation.
c
        mean = dble ( kk * nn ) / dble ( mm )

        sig = dsqrt ( mean * ( dble ( mm - nn ) / dble ( mm ) ) 
     &  * ( dble ( mm - kk ) / ( dble ( mm - 1 ) ) ) )

        if ( point ) then

          arg = - 0.5D+00 * ((( dble ( ll ) - mean ) / sig )**2 )
          if ( elimit .le. arg ) then
            chyper = dexp ( arg ) / ( sig * rootpi )
          else
            chyper = 0.0D+00
          end if

        else

          chyper = alnorm ( ( dble ( ll ) + 0.5D+00 - mean ) / sig, 
     &    .false. )

        end if

      else
c
c  Calculate exact hypergeometric probabilities.
c  Interchange K and N if this saves calculations.
c
        if ( min ( n - 1, m - n ) .lt. min ( k - 1, m - k ) ) then
          i = k
          k = n
          n = i
        end if

        if ( m - k .lt. k - 1 ) then
          dir = .not. dir
          l = n - l + 1
          k = m - k + 1
        end if

        if ( mbig .lt. mm ) then
c
c  Take logarithms of factorials.
c
          p = alnfac ( nn ) 
     &      - alnfac ( mm ) 
     &      + alnfac ( mm - kk ) 
     &      + alnfac ( kk ) 
     &      + alnfac ( mm - nn ) 
     &      - alnfac ( ll ) 
     &      - alnfac ( nn - ll ) 
     &      - alnfac ( kk - ll )
     *      - alnfac ( mm - nn - kk + ll )

          if ( elimit .le. p ) then
            chyper = dexp ( p )
          else
            chyper = 0.0D+00
          end if

        else
c
c  Use Freeman/Lund algorithm.
c
          do i = 1, l - 1
            chyper = chyper * dble ( ( k - i ) * ( n - i ) )
     &      / dble ( ( l - i ) * ( m - i ) )
          end do

          if ( l .ne. k ) then
            j = m - n + l
            do i = l, k - 1
              chyper = chyper * dble ( j - i ) / dble ( m - i )
            end do

          end if

        end if

        if ( point ) then
          return
        end if

        if ( chyper .eq. 0.0D+00 ) then
c
c  We must recompute the point probability since it has underflowed.
c
          if ( mm .le. mbig ) then
            p = alnfac ( nn ) 
     &        - alnfac ( mm ) 
     &        + alnfac ( kk ) 
     &        + alnfac ( mm - nn ) 
     &        - alnfac ( ll ) 
     &        - alnfac ( nn - ll ) 
     &        - alnfac ( kk - ll ) 
     &        - alnfac ( mm - nn - kk + ll ) 
     &        + alnfac ( mm - kk )
          end if

          p = p + dlog ( scale )

          if ( p .lt. elimit ) then
            ifault = 3
            if ( 
     &      dble ( nn * kk + nn + kk + 1 ) / dble ( mm + 2 ) 
     &      .lt. dble ( ll ) ) then
              chyper = 1.0D+00
            end if
            return
          else
            p = dexp ( p )
          end if

        else
c
c  Scale up at this point.
c
          p = chyper * scale

        end if

        pt = 0.0D+00
        nl = n - l
        kl = k - l
        mnkl = m - n - kl + 1

        if ( l .le. kl ) then

          do i = 1, l - 1
            p = p * dble ( ( l - i ) * ( mnkl - i ) ) /
     &      dble ( ( nl + i ) * ( kl + i ) )
            pt = pt + p
          end do

        else

          dir = .not. dir
          do j = 0, kl - 1
            p = p * dble ( ( nl - j ) * ( kl - j ) )
     &      / dble ( ( l + j ) * ( mnkl + j ) )
            pt = pt + p
          end do

        end if

        if ( p .eq. 0.0D+00 ) then
          ifault = 3
        end if

        if ( dir ) then
          chyper = chyper + ( pt / scale )
        else
          chyper = 1.0D+00 - ( pt / scale )
        end if

      end if

      return
      end
      subroutine hypergeometric_cdf_values ( n_data, sam, suc, pop,
     &  n, fx )

c*********************************************************************72
c
cc HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
c
c  Discussion:
c
c    CDF(X)(A,B) is the probability of at most X successes in A trials,
c    given that the probability of success on a single trial is B.
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`DiscreteDistributions`]
c      dist = HypergeometricDistribution [ sam, suc, pop ]
c      CDF [ dist, n ]
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
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer SAM, integer SUC, integer POP, the sample size, 
c    success size, and population parameters of the function.
c
c    Output, integer N, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision fx
      double precision fx_vec(n_max) 
      integer n
      integer n_data
      integer n_vec(n_max) 
      integer pop
      integer pop_vec(n_max) 
      integer sam
      integer sam_vec(n_max) 
      integer suc
      integer suc_vec(n_max) 

      save fx_vec
      save n_vec
      save pop_vec
      save sam_vec
      save suc_vec

      data fx_vec /
     &  0.6001858177500578D-01, 
     &  0.2615284665839845D+00, 
     &  0.6695237889132748D+00, 
     &  0.1000000000000000D+01, 
     &  0.1000000000000000D+01, 
     &  0.5332595856827856D+00, 
     &  0.1819495964117640D+00, 
     &  0.4448047017527730D-01, 
     &  0.9999991751316731D+00, 
     &  0.9926860896560750D+00, 
     &  0.8410799901444538D+00, 
     &  0.3459800113391901D+00, 
     &  0.0000000000000000D+00, 
     &  0.2088888139634505D-02, 
     &  0.3876752992448843D+00, 
     &  0.9135215248834896D+00 /
      data n_vec /
     &   7,  8,  9, 10, 
     &   6,  6,  6,  6, 
     &   6,  6,  6,  6, 
     &   0,  0,  0,  0 /
      data pop_vec /
     &  100, 100, 100, 100, 
     &  100, 100, 100, 100, 
     &  100, 100, 100, 100, 
     &  90,  200, 1000, 10000 /
      data sam_vec /
     &  10, 10, 10, 10, 
     &   6,  7,  8,  9, 
     &  10, 10, 10, 10, 
     &  10, 10, 10, 10 /
      data suc_vec /
     &  90, 90, 90, 90, 
     &  90, 90, 90, 90, 
     &  10, 30, 50, 70, 
     &  90, 90, 90, 90 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if
     
      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        sam = 0
        suc = 0
        pop = 0
        n = 0
        fx = 0.0D+00
      else
        sam = sam_vec(n_data)
        suc = suc_vec(n_data)
        pop = pop_vec(n_data)
        n = n_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine hypergeometric_pdf_values ( n_data, sam, suc, pop,
     &  n, fx )

c*********************************************************************72
c
cc HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
c
c  Discussion:
c
c    PDF(X)(A,B) is the probability of X successes in A trials,
c    given that the probability of success on a single trial is B.
c
c    In Mathematica, the function can be evaluated by:
c
c      dist = HypergeometricDistribution [ sam, suc, pop ]
c      PDF [ dist, n ]
c
c  Modified:
c
c    08 January 2008
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
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer SAM, integer SUC, integer POP, the sample size, 
c    success size, and population parameters of the function.
c
c    Output, integer N, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision fx
      double precision fx_vec(n_max) 
      integer n
      integer n_data
      integer n_vec(n_max) 
      integer pop
      integer pop_vec(n_max) 
      integer sam
      integer sam_vec(n_max) 
      integer suc
      integer suc_vec(n_max) 

      save fx_vec
      save n_vec
      save pop_vec
      save sam_vec
      save suc_vec

      data fx_vec /
     &  0.05179370533242827D+00,
     &  0.2015098848089788D+00,
     &  0.4079953223292903D+00,
     &  0.3304762110867252D+00,
     &  0.5223047493549780D+00,
     &  0.3889503452643453D+00,
     &  0.1505614239732950D+00,
     &  0.03927689321042477D+00,
     &  0.00003099828465518108D+00,
     &  0.03145116093938197D+00,
     &  0.2114132170316862D+00,
     &  0.2075776621999210D+00,
     &  0.0000000000000000D+00,
     &  0.002088888139634505D+00,
     &  0.3876752992448843D+00,
     &  0.9135215248834896D+00 /
      data n_vec /
     &   7,  8,  9, 10, 
     &   6,  6,  6,  6, 
     &   6,  6,  6,  6, 
     &   0,  0,  0,  0 /
      data pop_vec /
     &  100, 100, 100, 100, 
     &  100, 100, 100, 100, 
     &  100, 100, 100, 100, 
     &  90,  200, 1000, 10000 /
      data sam_vec /
     &  10, 10, 10, 10, 
     &   6,  7,  8,  9, 
     &  10, 10, 10, 10, 
     &  10, 10, 10, 10 /
      data suc_vec /
     &  90, 90, 90, 90, 
     &  90, 90, 90, 90, 
     &  10, 30, 50, 70, 
     &  90, 90, 90, 90 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if
     
      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        sam = 0
        suc = 0
        pop = 0
        n = 0
        fx = 0.0D+00
      else
        sam = sam_vec(n_data)
        suc = suc_vec(n_data)
        pop = pop_vec(n_data)
        n = n_vec(n_data)
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
