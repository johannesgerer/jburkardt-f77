      program main

c*********************************************************************72
c
cc TOMS451_PRB tests CHISQD.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_num
      integer p_num

      parameter ( n_num = 10 )
      parameter ( p_num = 5 )

      real chisqd
      integer i
      integer j
      integer n
      integer n_vec(n_num)
      real p
      real p_vec(p_num)
      real row(p_num)
      real table(n_num,p_num)

      save n_vec
      save p_vec
      save table

      data n_vec / 1, 2, 3, 4, 5, 10, 15, 20, 50, 100 /
      data p_vec /
     &  0.9995E+00, 0.9950E+00, 0.500E+00, 0.0010E+00, 0.0001E+00 /
c
c  The data for TABLE is listed one column at a time.  That is,
c  a consecutive pair of rows of data represent a COLUMN of the array.
c
      data table /
     &   0.000000,  0.001000,  0.015312,  0.063955,   0.158168,
     &   1.264941,  3.107881,  5.398208, 23.460876,  59.895508,
     &   0.000039,  0.010025,  0.071641,  0.206904,   0.411690,
     &   2.155869,  4.601008,  7.433892, 27.990784,  67.327621,
     &   0.454933,  1.386293,  2.365390,  3.356400,   4.351295,
     &   9.341794, 14.338853, 19.337418, 49.334930,  99.334122,
     &  10.827576, 13.815512, 16.268982, 18.467987,  20.515503,
     &  39.589081, 37.697662, 45.314896, 86.660767, 149.449051,
     &  15.135827, 18.420670, 21.106873, 23.510040,  25.744583,
     &  35.565170, 44.267853, 52.387360, 95.969482, 161.319733 /

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS451_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS451 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Compute values of the Chi-square quantile,'
      write ( *, '(a)' ) '    CHISQD ( P, N )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We print out a table using both CHISQD'
      write ( *, '(a)' ) '  and the values printed in the writeup.'
      write ( *, '(a)' ) ' '
      write ( *, '(5x,5(2x,f12.6))' ) ( p_vec(j), j = 1, p_num )
      write ( *, '(a)' ) ' '

      do i = 1, n_num

        n = n_vec(i)

        do j = 1, p_num

          p = p_vec(j)
          row(j) = chisqd ( p, n )

        end do

        write ( *, '(2x,i3,5(2x,f12.6))' )
     &    n_vec(i), ( row(j), j = 1, p_num )
        write ( *, '(2x,i3,5(2x,f12.6))' )
     &    n_vec(i), ( table(i,j), j = 1, p_num )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS451_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function gaussd ( p )

c*********************************************************************72
c
cc GAUSSD inverts the standard normal CDF.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**16.
c
c  Modified:
c
c    27 December 2004
c
c  Author:
c
c    Michael Wichura
c
c  Reference:
c
c    Michael Wichura,
c    The Percentage Points of the Normal Distribution,
c    Algorithm AS 241,
c    Applied Statistics,
c    Volume 37, Number 3, pages 477-484, 1988.
c
c  Parameters:
c
c    Input, real P, the value of the cumulative probability
c    densitity function.  0 < P < 1.  If P is outside this range, an
c    "infinite" value will be returned.
c
c    Output, real GAUSSD, the normal deviate value
c    with the property that the probability of a standard normal deviate being
c    less than or equal to the value is P.
c
      implicit none

      real a(8)
      real b(8)
      real c(8)
      real const1
      real const2
      real d(8)
      real e(8)
      real f(8)
      real gaussd
      real p
      real q
      real r
      real rpoly_value
      real split1
      real split2

      save a
      save b
      save c
      save const1
      save const2
      save d
      save e
      save f
      save split1
      save split2

      data a /
     &  3.3871328727963666080E+00,
     &  1.3314166789178437745E+02,
     &  1.9715909503065514427E+03,
     &  1.3731693765509461125E+04,
     &  4.5921953931549871457E+04,
     &  6.7265770927008700853E+04,
     &  3.3430575583588128105E+04,
     &  2.5090809287301226727E+03 /
      data b /
     &  1.0E+00,
     &  4.2313330701600911252E+01,
     &  6.8718700749205790830E+02,
     &  5.3941960214247511077E+03,
     &  2.1213794301586595867E+04,
     &  3.9307895800092710610E+04,
     &  2.8729085735721942674E+04,
     &  5.2264952788528545610E+03 /
      data c /
     &  1.42343711074968357734E+00,
     &  4.63033784615654529590E+00,
     &  5.76949722146069140550E+00,
     &  3.64784832476320460504E+00,
     &  1.27045825245236838258E+00,
     &  2.41780725177450611770E-01,
     &  2.27238449892691845833E-02,
     &  7.74545014278341407640E-04 /
      data const1 / 0.180625E+00 /
      data const2 / 1.6E+00 /
      data d /
     &  1.0E+00,
     &  2.05319162663775882187E+00,
     &  1.67638483018380384940E+00,
     &  6.89767334985100004550E-01,
     &  1.48103976427480074590E-01,
     &  1.51986665636164571966E-02,
     &  5.47593808499534494600E-04,
     &  1.05075007164441684324E-09 /
      data e /
     &  6.65790464350110377720E+00,
     &  5.46378491116411436990E+00,
     &  1.78482653991729133580E+00,
     &  2.96560571828504891230E-01,
     &  2.65321895265761230930E-02,
     &  1.24266094738807843860E-03,
     &  2.71155556874348757815E-05,
     &  2.01033439929228813265E-07 /
      data f /
     &  1.0E+00,
     &  5.99832206555887937690E-01,
     &  1.36929880922735805310E-01,
     &  1.48753612908506148525E-02,
     &  7.86869131145613259100E-04,
     &  1.84631831751005468180E-05,
     &  1.42151175831644588870E-07,
     &  2.04426310338993978564E-15 /
      data split1 / 0.425E+00 /
      data split2 / 5.0E+00 /

      if ( p <= 0.0E+00 ) then
        gaussd = -1.0E+30
        return
      end if

      if ( 1.0D+00 <= p ) then
        gaussd = 1.0E+30
        return
      end if

      q = p - 0.5E+00

      if ( abs ( q ) <= split1 ) then

        r = const1 - q * q
        gaussd = q * rpoly_value ( 8, a, r ) / rpoly_value ( 8, b, r )

      else

        if ( q < 0.0E+00 ) then
          r = p
        else
          r = 1.0E+00 - p
        end if

        if ( r <= 0.0E+00 ) then
          gaussd = -1.0E+00
          stop
        end if

        r = sqrt ( -log ( r ) )

        if ( r <= split2 ) then

          r = r - const2
          gaussd = rpoly_value ( 8, c, r ) / rpoly_value ( 8, d, r )

        else

          r = r - split2
          gaussd = rpoly_value ( 8, e, r ) / rpoly_value ( 8, f, r )

        end if

        if ( q < 0.0E+00 ) then
          gaussd = -gaussd
        end if

      end if

      return
      end
      function rpoly_value ( n, a, x )

c*********************************************************************72
c
cc RPOLY_VALUE evaluates a double precision polynomial.
c
c  Discussion:
c
c    For sanity's sake, the value of N indicates the NUMBER of
c    coefficients, or more precisely, the ORDER of the polynomial,
c    rather than the DEGREE of the polynomial.  The two quantities
c    differ by 1, but cause a great deal of confusion.
c
c    Given N and A, the form of the polynomial is:
c
c      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
c
c  Modified:
c
c    13 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, real A(N), the coefficients of the polynomial.
c    A(1) is the constant term.
c
c    Input, real X, the point at which the polynomial is
c    to be evaluated.
c
c    Output, real RPOLY_VALUE, the value of the polynomial at X.
c
      implicit none

      integer n

      real a(n)
      integer i
      real rpoly_value
      real x

      rpoly_value = 0.0E+00
      do i = n, 1, -1
        rpoly_value = rpoly_value * x + a(i)
      end do

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
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
