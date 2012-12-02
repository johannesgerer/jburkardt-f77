      function besei0 ( x )

c*********************************************************************72
c
cc BESEI0 evaluates the exponentially scaled Bessel I0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the modified Bessel 
c    function of the first kind of order zero multiplied by EXP(-ABS(X)). 
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESEI0, the value of the function.
c
      implicit none

      double precision besei0
      integer jint
      double precision result
      double precision x

      jint = 2
      call calci0 ( x, result, jint )
      besei0 = result

      return
      end
      function besei1 ( x )

c*********************************************************************72
c
cc BESEI1 evaluates the exponentially scaled Bessel I1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the
c    modified Bessel function of the first kind of order one
c    multiplied by EXP(-ABS(X)).
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESEI1, the value of the function.
c
      implicit none

      double precision besei1
      integer jint
      double precision result
      double precision x

      jint = 2
      call calci1 ( x, result, jint )
      besei1 = result

      return
      end
      function besek0 ( x )

c*********************************************************************72
c
cc BESEK0 evaluates the exponentially scaled Bessel K0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the
c    modified Bessel function of the second kind of order zero
c    multiplied by the exponential function.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c    0 < X.
c
c    Output, double precision BESK0, the value of the function.
c
      implicit none

      double precision besek0
      integer jint
      double precision result
      double precision x

      jint = 2
      call calck0 ( x, result, jint )
      besek0 = result

      return
      end
      function besek1 ( x )

c*********************************************************************72
c
cc BESEK1 evaluates the exponentially scaled Bessel K1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the
c    modified Bessel function of the second kind of order one
c    multiplied by the exponential function, for arguments
c    XLEAST <= ARG <= XMAX.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESEK1, the value of the function.
c
      implicit none

      double precision besek1
      integer jint
      double precision result
      double precision x

      jint = 2
      call calck1 ( x, result, jint )
      besek1 = result

      return
      end
      function besi0 ( x )

c*********************************************************************72
c
cc BESI0 evaluates the Bessel I0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for
c    modified Bessel functions of the first kind of order zero for
c    arguments ABS(ARG) <= XMAX.
c
c    See comments heading CALCI0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESI0, the value of the function.
c
      implicit none

      double precision besi0
      integer jint
      double precision result
      double precision x

      jint = 1
      call calci0 ( x, result, jint )
      besi0 = result

      return
      end
      function besi1 ( x )

c*********************************************************************72
c
cc BESI1 evaluates the Bessel I1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for
c    modified Bessel functions of the first kind of order one for
c    arguments ABS(ARG) <= XMAX.
c
c    See comments heading CALCI1.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESI1, the value of the function.
c
      implicit none

      double precision besi1
      integer jint
      double precision result
      double precision x

      jint = 1
      call calci1 ( x, result, jint )
      besi1 = result

      return
      end
      function besj0 ( x )

c*********************************************************************72
c
cc BESJ0 evaluates the Bessel J0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the first kind of order zero for arguments  |X| <= XMAX
c
c    See comments heading CALJY0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESJ0, the value of the function.
c
      implicit none

      double precision besj0
      integer jint
      double precision result
      double precision x

      jint = 0
      call caljy0 ( x, result, jint )
      besj0 = result

      return
      end
      function besj1 ( x )

c*********************************************************************72
c
cc BESJ1 evaluates the Bessel J1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the first kind of order zero for arguments  |X| <= XMAX
c
c    See comments heading CALJY1.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESJ1, the value of the function.
c
      implicit none

      double precision besj1
      integer jint
      double precision result
      double precision x

      jint = 0
      call caljy1 ( x, result, jint )
      besj1 = result

      return
      end
      function besk0 ( x )

c*********************************************************************72
c
cc BESK0 evaluates the Bessel K0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the
c    modified Bessel function of the second kind of order zero
c    for arguments 0.0 < ARG <= XMAX.
c
c    See comments heading CALCK0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESK0, the value of the function.
c
      implicit none

      double precision besk0
      integer jint
      double precision result
      double precision x

      jint = 1
      call calck0 ( x, result, jint )
      besk0 = result

      return
      end
      function besk1 ( x )

c*********************************************************************72
c
cc BESK1 evaluates the Bessel K1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for the
c    modified Bessel function of the second kind of order one
c    for arguments XLEAST <= ARG <= XMAX.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESK1, the value of the function.
c
      implicit none

      double precision besk1
      integer jint
      double precision result
      double precision x

      jint = 1
      call calck1 ( x, result, jint )
      besk1 = result

      return
      end
      function besy0 ( x )

c*********************************************************************72
c
cc BESY0 evaluates the Bessel Y0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the second kind of order zero for arguments 0 < X <= XMAX.
c
c    See comments heading CALJY0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESY0, the value of the function.
c
      implicit none

      double precision besy0
      integer jint
      double precision result
      double precision x

      jint = 1
      call caljy0 ( x, result, jint )
      besy0 = result

      return
      end
      function besy1 ( x )

c*********************************************************************72
c
cc BESY1 evaluates the Bessel Y1(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the second kind of order zero for arguments 0 < X <= XMAX.
c
c    See comments heading CALJY1.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESY1, the value of the function.
c
      implicit none

      double precision besy1
      integer jint
      double precision result
      double precision x

      jint = 1
      call caljy1 ( x, result, jint )
      besy1 = result

      return
      end
      subroutine calcei ( arg, result, jint )

c*********************************************************************72
c
cc CALCEI computes various exponential integrals.
c
c  Discussion:
c
c    This routine computes the exponential integrals Ei(x),
c    E1(x), and  exp(-x)*Ei(x) for real arguments x where
c
c           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
c    Ei(x) =
c          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
c
c    and where the first integral is a principal value integral.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Henry Thacher,
c    Rational Chebyshev Approximations for the Exponential 
c    Integral E1(x),
c    Mathematics of Computation,
c    Volume 22, Number 103, July 1968, pages 641-649.
c
c    William Cody, Henry Thacher,
c    Chebyshev Approximations for the Exponential 
c    Integral Ei(x),
c    Mathematics of Computation,
c    Volume 23, Number 106, April 1969, pages 289-303.
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  The argument must not
c    be zero.  If JINT = 2, then the argument must be strictly positive.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    1, RESULT = EI ( ARG );
c    2, RESULT = EONE ( ARG );
c    3, RESULT = EXPEI ( ARG ).
c
c    Input, integer JINT, chooses the function to be computed.
c    1, Ei(x);
c    2, -Ei(-x);
c    3, exp(-x)*Ei(x).
c
      implicit none

      double precision a(7)
      double precision arg
      double precision b(6)
      double precision c(9)
      double precision d(9)
      double precision e(10)
      double precision ei
      double precision exp40
      double precision f(10)
      double precision four
      double precision fourty
      double precision frac
      double precision half
      integer i
      integer jint
      double precision one
      double precision p(10)
      double precision plg(4)
      double precision px(10)
      double precision p037
      double precision p1(10)
      double precision p2(10)
      double precision q(10)
      double precision qlg(4)
      double precision qx(10)
      double precision q1(9)
      double precision q2(9)
      double precision r(10)
      double precision result
      double precision s(9)
      double precision six
      double precision sump
      double precision sumq
      double precision t
      double precision three
      double precision twelve
      double precision two
      double precision two4
      double precision w
      double precision x
      double precision x0
      double precision x01
      double precision x02
      double precision x11
      double precision xbig
      double precision xinf
      double precision xmax
      double precision xmx0
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c  EXP40 = exp(40)
c  X0 = zero of Ei
c  X01/X11 + X02 = zero of Ei to extra precision
c
      data zero /0.0d0 /
      data p037 /0.037d0/
      data half / 0.5d0 /
      data one / 1.0d0 /
      data two / 2.0d0/
      data three /3.0d0 /
      data four / 4.0d0 /
      data six / 6.0d0 /
      data twelve / 12.0d0 /
      data two4 / 24.0d0/
      data fourty / 40.0d0 /
      data exp40 / 2.3538526683701998541d17/
      data x01 /381.5d0 /
      data x11 / 1024.0d0 /
      data x02 /-5.1182968633365538008d-5/
      data x0 / 3.7250741078136663466d-1 /
c
c  Machine-dependent constants
c
      data xinf /1.79d+308/
      data xmax /716.351d0/
      data xbig /701.84d0/
c
c  Coefficients  for -1.0 <= X < 0.0
c
      data a /1.1669552669734461083368d2, 2.1500672908092918123209d3,
     &        1.5924175980637303639884d4, 8.9904972007457256553251d4,
     &        1.5026059476436982420737d5,-1.4815102102575750838086d5,
     &        5.0196785185439843791020d0/

      data b /4.0205465640027706061433d1, 7.5043163907103936624165d2,
     &        8.1258035174768735759855d3, 5.2440529172056355429883d4,
     &        1.8434070063353677359298d5, 2.5666493484897117319268d5/
c
c  Coefficients for -4.0 <= X < -1.0
c
      data c /3.828573121022477169108d-1, 1.107326627786831743809d+1,
     &        7.246689782858597021199d+1, 1.700632978311516129328d+2,
     &        1.698106763764238382705d+2, 7.633628843705946890896d+1,
     &        1.487967702840464066613d+1, 9.999989642347613068437d-1,
     &        1.737331760720576030932d-8/
      data d /8.258160008564488034698d-2, 4.344836335509282083360d+0,
     &        4.662179610356861756812d+1, 1.775728186717289799677d+2,
     &        2.953136335677908517423d+2, 2.342573504717625153053d+2,
     &        9.021658450529372642314d+1, 1.587964570758947927903d+1,
     &        1.000000000000000000000d+0/
c
c  Coefficients for X < -4.0
c
      data e /1.3276881505637444622987d+2,3.5846198743996904308695d+4,
     &        1.7283375773777593926828d+5,2.6181454937205639647381d+5,
     &        1.7503273087497081314708d+5,5.9346841538837119172356d+4,
     &       1.0816852399095915622498d+4,1.0611777263550331766871d03,
     &        5.2199632588522572481039d+1,9.9999999999999999087819d-1/
      data f /3.9147856245556345627078d+4,2.5989762083608489777411d+5,
     &        5.5903756210022864003380d+5,5.4616842050691155735758d+5,
     &        2.7858134710520842139357d+5,7.9231787945279043698718d+4,
     &        1.2842808586627297365998d+4,1.1635769915320848035459d+3,
     &        5.4199632588522559414924d+1,1.0d0/
c
c  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
c
      data plg /-2.4562334077563243311d+01,2.3642701335621505212d+02,
     &          -5.4989956895857911039d+02,3.5687548468071500413d+02/
      data qlg /-3.5553900764052419184d+01,1.9400230218539473193d+02,
     &          -3.3442903192607538956d+02,1.7843774234035750207d+02/
c
c  Coefficients for  0.0 < X < 6.0,
c  ratio of Chebyshev polynomials
c
      data p /-1.2963702602474830028590d01,-1.2831220659262000678155d03,
     &        -1.4287072500197005777376d04,-1.4299841572091610380064d06,
     &        -3.1398660864247265862050d05,-3.5377809694431133484800d08,
     &         3.1984354235237738511048d08,-2.5301823984599019348858d10,
     &         1.2177698136199594677580d10,-2.0829040666802497120940d11/
      data q / 7.6886718750000000000000d01,-5.5648470543369082846819d03,
     &         1.9418469440759880361415d05,-4.2648434812177161405483d06,
     &         6.4698830956576428587653d07,-7.0108568774215954065376d08,
     &         5.4229617984472955011862d09,-2.8986272696554495342658d10,
     &         9.8900934262481749439886d10,-8.9673749185755048616855d10/
c
c  J-fraction coefficients for 6.0 <= X < 12.0
c
      data r/-2.645677793077147237806d00,-2.378372882815725244124d00,
     &       -2.421106956980653511550d01, 1.052976392459015155422d01,
     &        1.945603779539281810439d01,-3.015761863840593359165d01,
     &        1.120011024227297451523d01,-3.988850730390541057912d00,
     &        9.565134591978630774217d00, 9.981193787537396413219d-1/
      data s/ 1.598517957704779356479d-4, 4.644185932583286942650d00,
     &        3.697412299772985940785d02,-8.791401054875438925029d00,
     &        7.608194509086645763123d02, 2.852397548119248700147d01,
     &        4.731097187816050252967d02,-2.369210235636181001661d02,
     &        1.249884822712447891440d00/
c
c  J-fraction coefficients for 12.0 <= X < 24.0
c
      data p1/-1.647721172463463140042d00,-1.860092121726437582253d01,
     &        -1.000641913989284829961d01,-2.105740799548040450394d01,
     &        -9.134835699998742552432d-1,-3.323612579343962284333d01,
     &         2.495487730402059440626d01, 2.652575818452799819855d01,
     &        -1.845086232391278674524d00, 9.999933106160568739091d-1/
      data q1/ 9.792403599217290296840d01, 6.403800405352415551324d01,
     &         5.994932325667407355255d01, 2.538819315630708031713d02,
     &         4.429413178337928401161d01, 1.192832423968601006985d03,
     &         1.991004470817742470726d02,-1.093556195391091143924d01,
     &         1.001533852045342697818d00/
c
c  J-fraction coefficients for  24 <= X.
c
      data p2/ 1.75338801265465972390d02,-2.23127670777632409550d02,
     &        -1.81949664929868906455d01,-2.79798528624305389340d01,
     &        -7.63147701620253630855d00,-1.52856623636929636839d01,
     &        -7.06810977895029358836d00,-5.00006640413131002475d00,
     &        -3.00000000320981265753d00, 1.00000000000000485503d00/
      data q2/ 3.97845977167414720840d04, 3.97277109100414518365d00,
     &         1.37790390235747998793d02, 1.17179220502086455287d02,
     &         7.04831847180424675988d01,-1.20187763547154743238d01,
     &        -7.99243595776339741065d00,-2.99999894040324959612d00,
     &         1.99999999999048104167d00/

      x = arg

      if ( x .eq. zero ) then

        ei = -xinf

        if ( jint .eq. 2 ) then
          ei = -ei
        end if
c
c  Calculate EI for negative argument or for E1.
c
      else if ( x .lt. zero .or. jint .eq. 2 ) then

        y = abs ( x )

        if ( y .le. one ) then

          sump = a(7) * y + a(1)
          sumq = y + b(1)
          do i = 2, 6
            sump = sump * y + a(i)
            sumq = sumq * y + b(i)
          end do
          ei = log ( y ) - sump / sumq

          if ( jint .eq. 3 ) then
            ei = ei * exp ( y )
          end if

        else if ( y .le. four ) then

          w = one / y
          sump = c(1)
          sumq = d(1)
          do i = 2, 9
            sump = sump * w + c(i)
            sumq = sumq * w + d(i)
          end do

          ei = - sump / sumq

          if ( jint .ne. 3 ) then
            ei = ei * exp ( -y )
          end if

        else

          if ( xbig .lt. y .and. jint .lt. 3 ) then

            ei = zero

          else

            w = one / y
            sump = e(1)
            sumq = f(1)
            do i = 2, 10
              sump = sump * w + e(i)
              sumq = sumq * w + f(i)
            end do

            ei = -w * ( one - w * sump / sumq )

            if ( jint .ne. 3 ) then
              ei = ei * exp ( - y )
            end if

          end if

        end if

        if ( jint .eq. 2 ) then
          ei = -ei
        end if
c
c  To improve conditioning, rational approximations are expressed
c  in terms of Chebyshev polynomials for 0 <= X < 6, and in
c  continued fraction form for larger X.
c
      else if ( x .lt. six ) then

        t = x + x
        t = t / three - two
        px(1) = zero
        qx(1) = zero
        px(2) = p(1)
        qx(2) = q(1)
        do i = 2, 9
          px(i+1) = t * px(i) - px(i-1) + p(i)
          qx(i+1) = t * qx(i) - qx(i-1) + q(i)
        end do
        sump = half * t * px(10) - px(9) + p(10)
        sumq = half * t * qx(10) - qx(9) + q(10)
        frac = sump / sumq
        xmx0 = ( x - x01 / x11 ) - x02

        if ( abs ( xmx0 ) .ge. p037 ) then

          ei = log ( x / x0 ) + xmx0 * frac

          if ( jint .eq. 3 ) then
            ei = exp ( - x ) * ei
          end if
c
c  Special approximation to ln(X/X0) for X close to X0.
c
        else

          y = xmx0 / ( x + x0 )
          ysq = y * y
          sump = plg(1)
          sumq = ysq + qlg(1)
          do i = 2, 4
            sump = sump * ysq + plg(i)
            sumq = sumq * ysq + qlg(i)
          end do
          ei = ( sump / ( sumq * ( x + x0 ) ) + frac ) * xmx0

          if ( jint .eq. 3 ) then
            ei = exp ( - x ) * ei
          end if

        end if

      else if ( x .lt. twelve ) then

        frac = zero
        do i = 1, 9
          frac = s(i) / ( r(i) + x + frac )
        end do

        ei = ( r(10) + frac ) / x
        if ( jint .ne. 3 ) then
          ei = ei * exp ( x )
        end if

      else if ( x .le. two4 ) then

        frac = zero
        do i = 1, 9
          frac = q1(i) / ( p1(i) + x + frac )
        end do

        ei = ( p1(10) + frac ) / x

        if ( jint .ne. 3 ) then
          ei = ei * exp ( x )
        end if

      else

        if ( x .ge. xmax .and. jint .lt. 3 ) then

          ei = xinf

        else

          y = one / x
          frac = zero
          do i = 1, 9
            frac = q2(i) / ( p2(i) + x + frac )
          end do
          frac = p2(10) + frac
          ei = y + y * y * frac

          if ( jint .ne. 3 ) then

            if ( x .le. xmax - two4 ) then
              ei = ei * exp ( x )
c
c  Calculation reformulated to avoid premature overflow.
c
            else
              ei = ( ei * exp ( x - fourty ) ) * exp40
            end if

          end if
        end if
      end if

      result = ei

      return
      end
      subroutine calci0 ( arg, result, jint )

c*********************************************************************72
c
cc CALCI0 computes various I0 Bessel functions.
c
c  Discussion:
c
c    This routine computes modified Bessel functions of the first kind
c    and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
c    arguments X.  
c
c    The main computation evaluates slightly modified forms of
c    minimax approximations generated by Blair and Edwards, Chalk
c    River (Atomic Energy of Canada Limited) Report AECL-4928,
c    October, 1974. 
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT = 1, then
c    the argument must be less than XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    1, RESULT = I0(x);
c    2, RESULT = exp(-x) * I0(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    1, I0(x);
c    2, exp(-x) * I0(x);
c
      implicit none

      double precision a
      double precision arg
      double precision b
      double precision exp40
      double precision forty
      integer i
      integer jint
      double precision one
      double precision one5
      double precision p(15)
      double precision pp(8)
      double precision q(5)
      double precision qq(7)
      double precision result
      double precision rec15
      double precision sump
      double precision sumq
      double precision two25
      double precision x
      double precision xinf
      double precision xmax
      double precision xsmall
      double precision xx
c
c  Mathematical constants
c
      data one /1.0d0/
      data one5 /15.0d0/
      data exp40 /2.353852668370199854d17/
      data forty /40.0d0/
      data rec15 /6.6666666666666666666d-2/
      data two25 /225.0d0/
c
c  Machine-dependent constants
c
      data xsmall /5.55d-17/
      data xinf /1.79d308/
      data xmax /713.986d0/
c
c  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
c
      data  p/-5.2487866627945699800d-18,-1.5982226675653184646d-14,
     &        -2.6843448573468483278d-11,-3.0517226450451067446d-08,
     &        -2.5172644670688975051d-05,-1.5453977791786851041d-02,
     &        -7.0935347449210549190d+00,-2.4125195876041896775d+03,
     &        -5.9545626019847898221d+05,-1.0313066708737980747d+08,
     &        -1.1912746104985237192d+10,-8.4925101247114157499d+11,
     &        -3.2940087627407749166d+13,-5.5050369673018427753d+14,
     &        -2.2335582639474375249d+15/
      data  q/-3.7277560179962773046d+03, 6.5158506418655165707d+06,
     &        -6.5626560740833869295d+09, 3.7604188704092954661d+12,
     &        -9.7087946179594019126d+14/
c
c  Coefficients for 15.0 <= ABS(ARG)
c
      data pp/-3.9843750000000000000d-01, 2.9205384596336793945d+00,
     &        -2.4708469169133954315d+00, 4.7914889422856814203d-01,
     &        -3.7384991926068969150d-03,-2.6801520353328635310d-03,
     &         9.9168777670983678974d-05,-2.1877128189032726730d-06/
      data qq/-3.1446690275135491500d+01, 8.5539563258012929600d+01,
     &        -6.0228002066743340583d+01, 1.3982595353892851542d+01,
     &        -1.1151759188741312645d+00, 3.2547697594819615062d-02,
     &        -5.5194330231005480228d-04/

      x = abs ( arg )

      if ( x .lt. xsmall ) then

        result = one
c
c  XSMALL <= ABS(ARG) < 15.0.
c
      else if ( x .lt. one5 ) then

        xx = x * x
        sump = p(1)
        do i = 2, 15
          sump = sump * xx + p(i)
        end do
        xx = xx - two25

        sumq = ((((
     &      xx + q(1) ) 
     &    * xx + q(2) )
     &    * xx + q(3) )
     &    * xx + q(4) )
     &    * xx + q(5)

        result = sump / sumq

        if ( jint .eq. 2 ) then
          result = result * exp ( - x )
        end if

      else if ( x .ge. one5 ) then

        if ( jint .eq. 1 .and. xmax .lt. x ) then
          result = xinf
        else
c
c  15.0 <= ABS(ARG).
c
          xx = one / x - rec15

          sump = ((((((
     &             pp(1)
     &      * xx + pp(2) )
     &      * xx + pp(3) )
     &      * xx + pp(4) )
     &      * xx + pp(5) )
     &      * xx + pp(6) )
     &      * xx + pp(7) )
     &      * xx + pp(8)

          sumq = ((((((
     &        xx + qq(1) )
     &      * xx + qq(2) )
     &      * xx + qq(3) )
     &      * xx + qq(4) )
     &      * xx + qq(5) )
     &      * xx + qq(6) )
     &      * xx + qq(7)

          result = sump / sumq

          if ( jint .eq. 2 ) then
            result = ( result - pp(1) ) / sqrt ( x )
          else
c
c  Calculation reformulated to avoid premature overflow.
c
            if ( x .le.( xmax - one5 ) ) then
              a = exp ( x )
              b = one
            else
              a = exp ( x - forty )
              b = exp40
            end if

            result = ( ( result * a - pp(1) * a ) / sqrt ( x ) ) * b

          end if

        end if

      end if

      return
      end
      subroutine calci1 ( arg, result, jint )

c*********************************************************************72
c
cc CALCI1 computes various I1 Bessel functions.
c
c  Discussion:
c
c    This routine computes modified Bessel functioons of the first kind
c    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
c    arguments X.
c
c    The main computation evaluates slightly modified forms of
c    minimax approximations generated by Blair and Edwards, Chalk
c    River (Atomic Energy of Canada Limited) Report AECL-4928,
c    October, 1974.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT = 1, then
c    the argument must be less than XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    1, RESULT = I1(x);
c    2, RESULT = exp(-x) * I1(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    1, I1(x);
c    2, exp(-x) * I1(x);
c
      implicit none

      double precision a
      double precision arg
      double precision b
      double precision exp40
      double precision forty
      double precision half
      integer j
      integer jint
      double precision one
      double precision one5
      double precision p(15)
      double precision pbar
      double precision pp(8)
      double precision q(5)
      double precision qq(6)
      double precision rec15
      double precision result
      double precision sump
      double precision sumq
      double precision two25
      double precision x
      double precision xinf
      double precision xmax
      double precision xsmall
      double precision xx
      double precision zero
c
c  Mathematical constants
c
      data one /1.0d0/
      data one5 /15.0d0/
      data exp40 /2.353852668370199854d17/
      data forty /40.0d0/
      data rec15 /6.6666666666666666666d-2/
      data two25 /225.0d0/
      data half /0.5d0/
      data zero /0.0d0/
c
c  Machine-dependent constants
c
      data xsmall /5.55d-17/
      data xinf /1.79d308/
      data xmax /713.987d0/
c
c  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
c
      data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16,
     &       -1.1928788903603238754d-12,-1.4831904935994647675d-09,
     &       -1.3466829827635152875d-06,-9.1746443287817501309d-04,
     &       -4.7207090827310162436d-01,-1.8225946631657315931d+02,
     &       -5.1894091982308017540d+04,-1.0588550724769347106d+07,
     &       -1.4828267606612366099d+09,-1.3357437682275493024d+11,
     &       -6.9876779648010090070d+12,-1.7732037840791591320d+14,
     &       -1.4577180278143463643d+15/
      data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06,
     &       -8.0059518998619764991d+09, 4.8544714258273622913d+12,
     &       -1.3218168307321442305d+15/
c
c  Coefficients for 15.0 .LE. ABS(ARG)
c
      data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01,
     &        -4.2843766903304806403d-01, 9.7356000150886612134d-02,
     &        -3.2457723974465568321d-03,-3.6395264712121795296d-04,
     &         1.6258661867440836395d-05,-3.6347578404608223492d-07/
      data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00,
     &        -8.5017476463217924408d-01, 7.4212010813186530069d-02,
     &        -2.2835624489492512649d-03, 3.7510433111922824643d-05/
      data pbar/3.98437500d-01/

      x = abs ( arg )
c
c  Return for ABS(ARG) < XSMALL.
c
      if ( x .lt. xsmall ) then

        result = half * x
c
c  XSMALL <= ABS(ARG) < 15.0.
c
      else if ( x .lt. one5 ) then

        xx = x * x
        sump = p(1)
        do j = 2, 15
          sump = sump * xx + p(j)
        end do
        xx = xx - two25

        sumq = ((((
     &      xx + q(1) )
     &    * xx + q(2) )
     &    * xx + q(3) )
     &    * xx + q(4) )
     &    * xx + q(5)

        result = ( sump / sumq ) * x

        if ( jint .eq. 2 ) then
          result = result * exp ( -x )
        end if

      else if ( jint .eq. 1 .and. xmax .lt. x ) then

        result = xinf

      else
c
c  15.0 <= ABS(ARG).
c
        xx = one / x - rec15

        sump = ((((((
     &           pp(1)
     &    * xx + pp(2) ) 
     &    * xx + pp(3) )
     &    * xx + pp(4) )
     &    * xx + pp(5) )
     &    * xx + pp(6) )
     &    * xx + pp(7) )
     &    * xx + pp(8)

        sumq = (((((
     &      xx + qq(1) )
     &    * xx + qq(2) )
     &    * xx + qq(3) )
     &    * xx + qq(4) )
     &    * xx + qq(5) )
     &    * xx + qq(6)

        result = sump / sumq

        if ( jint .ne. 1 ) then
          result = ( result + pbar ) / sqrt ( x )
        else
c
c  Calculation reformulated to avoid premature overflow.
c
          if ( xmax - one5 .lt. x ) then
            a = exp ( x - forty )
            b = exp40
          else
            a = exp ( x )
            b = one
          end if

          result = ( ( result * a + pbar * a ) / sqrt ( x ) ) * b

        end if
      end if

      if ( arg .lt. zero ) then
        result = -result
      end if

      return
      end
      subroutine calck0 ( arg, result, jint )

c*********************************************************************72
c
cc CALCK0 computes various K0 Bessel functions.
c
c  Discussion:
c
c    This routine computes modified Bessel functions of the second kind
c    and order zero, K0(X) and EXP(X)*K0(X), for real
c    arguments X.  
c
c    The main computation evaluates slightly modified forms of near
c    minimax rational approximations generated by Russon and Blair,
c    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
c    1969.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  0 < ARG is
c    always required.  If JINT = 1, then the argument must also be 
c    less than XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    1, RESULT = K0(x);
c    2, RESULT = exp(x) * K0(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    1, K0(x);
c    2, exp(x) * K0(x);
c
      implicit none

      integer i
      integer jint
      double precision arg
      double precision f(4)
      double precision g(3)
      double precision one
      double precision p(6)
      double precision pp(10)
      double precision q(2)
      double precision qq(10)
      double precision result
      double precision sumf
      double precision sumg
      double precision sump
      double precision sumq
      double precision temp
      double precision x
      double precision xinf
      double precision xmax
      double precision xsmall
      double precision xx
      double precision zero
c
c  Mathematical constants
c
      data one /1.0d0/
      data zero /0.0d0/
c
c  Machine-dependent constants
c
      data xsmall /1.11d-16/
      data xinf /1.79d+308/
      data xmax /705.342d0/
c
c  Coefficients for XSMALL <= ARG <= 1.0
c
      data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01,
     &          1.1999463724910714109d+01, 4.6850901201934832188d+02,
     &          5.9169059852270512312d+03, 2.4708152720399552679d+03/
      data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
      data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02,
     &         -1.7733784684952985886d+04,-4.0320340761145482298d+05/
      data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04,
     &         -1.6128136304458193998d+06/
c
c  Coefficients for  1.0 < ARG
c
      data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03,
     &          3.1075408980684392399d+04, 1.0577068948034021957d+05,
     &          1.7398867902565686251d+05, 1.5097646353289914539d+05,
     &          7.1557062783764037541d+04, 1.8321525870183537725d+04,
     &          2.3444738764199315021d+03, 1.1600249425076035558d+02/
      data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03,
     &          3.1474655750295278825d+04, 9.7418829762268075784d+04,
     &          1.5144644673520157801d+05, 1.2689839587977598727d+05,
     &          5.8824616785857027752d+04, 1.4847228371802360957d+04,
     &          1.8821890840982713696d+03, 9.2556599177304839811d+01/

      x = arg
c
c  0.0 < ARG <= 1.0.
c
      if ( zero .lt. x ) then

        if ( x .le. one ) then

          temp = log ( x )

          if ( x .lt. xsmall ) then
c
c  Return for small ARG.
c
            result = p(6) / q(2) - temp

          else

            xx = x * x

            sump = ((((
     &               p(1)
     &        * xx + p(2) )
     &        * xx + p(3) )
     &        * xx + p(4) )
     &        * xx + p(5) )
     &        * xx + p(6)

            sumq = ( xx + q(1) ) * xx + q(2)
            sumf = ( ( 
     &               f(1) 
     &        * xx + f(2) )
     &        * xx + f(3) )
     &        * xx + f(4)

            sumg = ( ( xx + g(1) ) * xx + g(2) ) * xx + g(3)

            result = sump / sumq - xx * sumf * temp / sumg - temp

            if ( jint .eq. 2 ) then
              result = result * exp ( x )
            end if

          end if

        else if ( jint .eq. 1 .and. xmax .lt. x ) then
c
c  Error return for XMAX < ARG.
c
          result = zero

        else
c
c  1.0 < ARG.
c
          xx = one / x
          sump = pp(1)
          do i = 2, 10
            sump = sump * xx + pp(i)
          end do

          sumq = xx
          do i = 1, 9
            sumq = ( sumq + qq(i) ) * xx
          end do
          sumq = sumq + qq(10)
          result = sump / sumq / sqrt ( x )

          if ( jint .eq. 1 ) then
            result = result * exp ( -x )
          end if

        end if

      else
c
c  Error return for ARG <= 0.0.
c
        result = xinf

      end if

      return
      end
      subroutine calck1 ( arg, result, jint )

c*********************************************************************72
c
cc CALCK1 computes various K1 Bessel functions.
c
c  Discussion:
c
c    This routine computes modified Bessel functions of the second kind
c    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
c
c    The main computation evaluates slightly modified forms of near
c    minimax rational approximations generated by Russon and Blair,
c    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
c    1969.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  XLEAST < ARG is
c    always required.  If JINT = 1, then the argument must also be 
c    less than XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    1, RESULT = K1(x);
c    2, RESULT = exp(x) * K1(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    1, K1(x);
c    2, exp(x) * K1(x);
c
      implicit none

      double precision arg
      double precision f(5)
      double precision g(3)
      integer i
      integer jint
      double precision one
      double precision p(5)
      double precision pp(11)
      double precision q(3)
      double precision qq(9)
      double precision result
      double precision sumf
      double precision sumg
      double precision sump
      double precision sumq
      double precision x
      double precision xinf
      double precision xmax
      double precision xleast
      double precision xsmall
      double precision xx
      double precision zero
c
c  Mathematical constants
c
      data one /1.0d0/
      data zero /0.0d0/
c
c  Machine-dependent constants
c
      data xleast /2.23d-308/
      data xsmall /1.11d-16/
      data xinf /1.79d+308/
      data xmax /705.343d+0/
c
c  Coefficients for  XLEAST .LE.  ARG  .LE. 1.0
c
      data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1,
     &          7.1885382604084798576d+3, 1.7733324035147015630d+5,
     &          7.1938920065420586101d+5/
      data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4,
     &         -2.2149374878243304548d+6/
      data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1,
     &         -4.5051623763436087023d+3,-1.4758069205414222471d+5,
     &         -1.3531161492785421328d+6/
      data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4,
     &         -2.7062322985570842656d+6/
c
c  Coefficients for  1.0 < ARG
c
      data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0,
     &          1.3182609918569941308d+2, 8.1094256146537402173d+2,
     &          2.3123742209168871550d+3, 3.4540675585544584407d+3,
     &          2.8590657697910288226d+3, 1.3319486433183221990d+3,
     &          3.4122953486801312910d+2, 4.4137176114230414036d+1,
     &          2.2196792496874548962d+0/
      data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2,
     &          1.2082692316002348638d+3, 2.1181000487171943810d+3,
     &          1.9448440788918006154d+3, 9.6929165726802648634d+2,
     &          2.5951223655579051357d+2, 3.4552228452758912848d+1,
     &          1.7710478032601086579d+0/

      x = arg
c
c  Error return for ARG < XLEAST.
c
      if ( x .lt. xleast ) then

        result = xinf
c
c  XLEAST <= ARG <= 1.0.
c
      else if ( x .le. one ) then

        if ( x .lt. xsmall ) then
c
c  Return for small ARG.
c
          result = one / x

        else

          xx = x * x

          sump = (((( 
     &             p(1) 
     &      * xx + p(2) ) 
     &      * xx + p(3) ) 
     &      * xx + p(4) )
     &      * xx + p(5) ) 
     &      * xx + q(3)

          sumq = ((
     &        xx + q(1) )
     &      * xx + q(2) )
     &      * xx + q(3)

          sumf = (((
     &             f(1)
     &      * xx + f(2) )
     &      * xx + f(3) )
     &      * xx + f(4) )
     &      * xx + f(5)

          sumg = ((
     &        xx + g(1) )
     &      * xx + g(2) )
     &      * xx + g(3)

          result = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x

          if ( jint .eq. 2 ) then
            result = result * exp ( x )
          end if

        end if

      else if ( jint .eq. 1 .and. xmax .lt. x ) then
c
c  Error return for XMAX < ARG.
c
        result = zero

      else
c
c  1.0 < ARG.
c
        xx = one / x

        sump = pp(1)
        do i = 2, 11
          sump = sump * xx + pp(i)
        end do

        sumq = xx
        do i = 1, 8
          sumq = ( sumq + qq(i) ) * xx
        end do
        sumq = sumq + qq(9)

        result = sump / sumq / sqrt ( x )

        if ( jint .eq. 1 ) then
          result = result * exp ( -x )
        end if

      end if

      return
      end
      subroutine calerf ( arg, result, jint )

c*********************************************************************72
c
cc CALERF computes various forms of the error function.
c
c  Discussion:
c
c    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
c    for a real argument x.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody,
c    Rational Chebyshev Approximations for the Error Function,
c    Mathematics of Computation,
c    Volume 23, Number 107, July 1969, pages 631-638.
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT is 1, the
c    argument must be less than XBIG.  If JINT is 2, the argument
c    must lie between XNEG and XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    0, RESULT = erf(x);
c    1, RESULT = erfc(x) = 1 - erf(x);
c    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
c
c    Input, integer JINT, chooses the function to be computed.
c    0, erf(x);
c    1, erfc(x);
c    2, exp(x*x)*erfc(x).
c
      implicit none

      double precision a(5)
      double precision arg
      double precision b(4)
      double precision c(9)
      double precision d(8)
      double precision del
      double precision four
      double precision half
      integer i
      integer jint
      double precision one
      double precision p(6)
      double precision q(5)
      double precision result
      double precision sixten
      double precision sqrpi
      double precision two
      double precision thresh
      double precision x
      double precision xbig
      double precision xden
      double precision xhuge
      double precision xinf
      double precision xmax
      double precision xneg
      double precision xnum
      double precision xsmall
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c
      data four /4.0d0/
      data one /1.0d0/
      data half /0.5d0/
      data two /2.0d0/
      data zero /0.0d0/
      data sqrpi /5.6418958354775628695d-1/
      data thresh /0.46875d0/
      data sixten /16.0d0/
c
c  Machine-dependent constants
c
      data xinf /1.79d308 /
      data xneg / -26.628d0 /
      data xsmall /1.11d-16/
      data xbig /26.543d0 /
      data xhuge /6.71d7/
      data xmax /2.53d307/
c
c  Coefficients for approximation to  erf  in first interval
c
      data a/3.16112374387056560d00,1.13864154151050156d02,
     &       3.77485237685302021d02,3.20937758913846947d03,
     &       1.85777706184603153d-1/
      data b/2.36012909523441209d01,2.44024637934444173d02,
     &       1.28261652607737228d03,2.84423683343917062d03/
c
c  Coefficients for approximation to  erfc  in second interval
c
      data c/5.64188496988670089d-1,8.88314979438837594d0,
     &       6.61191906371416295d01,2.98635138197400131d02,
     &       8.81952221241769090d02,1.71204761263407058d03,
     &       2.05107837782607147d03,1.23033935479799725d03,
     &       2.15311535474403846d-8/
      data d/1.57449261107098347d01,1.17693950891312499d02,
     &       5.37181101862009858d02,1.62138957456669019d03,
     &       3.29079923573345963d03,4.36261909014324716d03,
     &       3.43936767414372164d03,1.23033935480374942d03/
c
c  Coefficients for approximation to  erfc  in third interval
c
      data p/3.05326634961232344d-1,3.60344899949804439d-1,
     &       1.25781726111229246d-1,1.60837851487422766d-2,
     &       6.58749161529837803d-4,1.63153871373020978d-2/
      data q/2.56852019228982242d00,1.87295284992346047d00,
     &       5.27905102951428412d-1,6.05183413124413191d-2,
     &       2.33520497626869185d-3/

      x = arg
      y = abs ( x )
c
c  Evaluate erf for |X| <= 0.46875.
c
      if ( y .le. thresh ) then

        ysq = zero
        if ( xsmall .lt. y ) then
          ysq = y * y
        end if

        xnum = a(5) * ysq
        xden = ysq

        do i = 1, 3
          xnum = ( xnum + a(i) ) * ysq
          xden = ( xden + b(i) ) * ysq
        end do

        result = x * ( xnum + a(4) ) / ( xden + b(4) )

        if ( jint .ne. 0 ) then
          result = one - result
        end if

        if ( jint .eq. 2 ) then
          result = exp ( ysq ) * result
        end if

        return
c
c  Evaluate erfc for 0.46875 <= |X| <= 4.0.
c
       else if ( y .le. four ) then

         xnum = c(9) * y
         xden = y

         do i = 1, 7
           xnum = ( xnum + c(i) ) * y
           xden = ( xden + d(i) ) * y
         end do

         result = ( xnum + c(8) ) / ( xden + d(8) )

         if ( jint .ne. 2 ) then
           ysq = aint ( y * sixten ) / sixten
           del = ( y - ysq ) * ( y + ysq )
           result = exp ( -ysq * ysq ) * exp ( -del ) * result
         end if
c
c  Evaluate erfc for 4.0 < |X|.
c
       else

         result = zero

         if ( y .ge. xbig ) then

           if ( jint .ne. 2 .or. y .ge. xmax ) then
             go to 300
           end if

           if ( y .ge. xhuge ) then
             result = sqrpi / y
             go to 300
           end if

         end if

         ysq = one / ( y * y )
         xnum = p(6) * ysq
         xden = ysq
         do i = 1, 4
           xnum = ( xnum + p(i) ) * ysq
           xden = ( xden + q(i) ) * ysq
          end do

          result = ysq * ( xnum + p(5) ) / ( xden + q(5) )
          result = ( sqrpi -  result ) / y

          if ( jint .ne. 2 ) then
            ysq = aint ( y * sixten ) / sixten
            del = ( y - ysq ) * ( y + ysq )
            result = exp ( -ysq * ysq ) * exp ( -del ) * result
          end if

      end if
c
c  Fix up for negative argument, erf, etc.
c
  300 continue

      if ( jint .eq. 0 ) then

        result = ( half - result ) + half
        if ( x .lt. zero ) then
          result = -result
        end if

      else if ( jint .eq. 1 ) then

        if ( x .lt. zero ) then
          result = two - result
        end if

      else

        if ( x .lt. zero ) then

          if ( x .lt. xneg ) then
            result = xinf
          else
            ysq = aint ( x * sixten ) / sixten
            del = ( x - ysq ) * ( x + ysq )
            y = exp ( ysq * ysq ) * exp ( del )
            result = ( y + y ) - result
          end if

        end if

      end if

      return
      end
      subroutine caljy0 ( arg, result, jint )

c*********************************************************************72
c
cc CALJY0 computes various J0 and Y0 Bessel functions.
c
c  Discussion:
c
c    This routine computes zero-order Bessel functions of the first and
c    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
c    for Y0, and |X| <= XMAX for J0.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT = 0, ARG
c    must satisfy 
c     -XMAX < ARG < XMAX;
c    If JINT = 1, then ARG must satisfy
c      0 < ARG < XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    0, RESULT = J0(x);
c    1, RESULT = Y0(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    0, J0(x);
c    1, Y0(x);
c
      implicit none

      integer i
      integer jint
      double precision arg
      double precision ax
      double precision cons
      double precision down
      double precision eight
      double precision five5
      double precision four
      double precision one
      double precision oneov8
      double precision pi2
      double precision pj0(7)
      double precision pj1(8)
      double precision plg(4)
      double precision prod
      double precision py0(6)
      double precision py1(7)
      double precision py2(8)
      double precision p0(6)
      double precision p1(6)
      double precision p17
      double precision qj0(5)
      double precision qj1(7)
      double precision qlg(4)
      double precision qy0(5)
      double precision qy1(6)
      double precision qy2(7)
      double precision q0(5)
      double precision q1(5)
      double precision resj
      double precision result
      double precision r0
      double precision r1
      double precision sixty4
      double precision three
      double precision twopi
      double precision twopi1
      double precision twopi2
      double precision two56
      double precision up
      double precision w
      double precision wsq
      double precision xden
      double precision xinf
      double precision xmax
      double precision xnum
      double precision xsmall
      double precision xj0
      double precision xj1
      double precision xj01
      double precision xj02
      double precision xj11
      double precision xj12
      double precision xy
      double precision xy0
      double precision xy01
      double precision xy02
      double precision xy1
      double precision xy11
      double precision xy12
      double precision xy2
      double precision xy21
      double precision xy22
      double precision z
      double precision zero
      double precision zsq
c
c  Mathematical constants
c  CONS = ln(.5) + Euler's gamma
c
      data zero / 0.0d0 /
      data one /1.0d0 /
      data three /3.0d0 /
      data four /4.0d0 /
      data eight /8.0d0/
      data five5 / 5.5d0 /
      data sixty4 /64.0d0 /
      data oneov8 /0.125d0 /
      data p17 /1.716d-1/
      data two56 /256.0d0/
      data cons / -1.1593151565841244881d-1/
      data pi2 /6.3661977236758134308d-1/
      data twopi /6.2831853071795864769d0/
      data twopi1 /6.28125d0 /
      data twopi2 / 1.9353071795864769253d-3/
c
c  Machine-dependent constants
c
      data xmax /1.07d+09/
      data xsmall /9.31d-10/
      data xinf /1.7d+38/
c
c  Zeroes of Bessel functions
c
      data xj0 /2.4048255576957727686d+0/
      data xj1 /5.5200781102863106496d+0/
      data xy0 /8.9357696627916752158d-1/
      data xy1 /3.9576784193148578684d+0/
      data xy2 /7.0860510603017726976d+0/
      data xj01 / 616.0d+0/
      data xj02 /-1.4244423042272313784d-03/
      data xj11 /1413.0d+0/
      data xj12 / 5.4686028631064959660d-04/
      data xy01 / 228.0d+0/
      data xy02 / 2.9519662791675215849d-03/
      data xy11 /1013.0d+0/
      data xy12 / 6.4716931485786837568d-04/
      data xy21 /1814.0d+0/
      data xy22 / 1.1356030177269762362d-04/
c
c  Coefficients for rational approximation to ln(x/a)
c
      data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02,
     &         -5.4989956895857911039d+02,3.5687548468071500413d+02/
      data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02,
     &         -3.3442903192607538956d+02,1.7843774234035750207d+02/
c
c  Coefficients for rational approximation of
c  J0(X) / (X**2 - XJ0**2),  XSMALL < |X| <= 4.0
c
      data pj0/6.6302997904833794242d+06,-6.2140700423540120665d+08,
     &         2.7282507878605942706d+10,-4.1298668500990866786d+11,
     &        -1.2117036164593528341d-01, 1.0344222815443188943d+02,
     &        -3.6629814655107086448d+04/
      data qj0/4.5612696224219938200d+05, 1.3985097372263433271d+08,
     &         2.6328198300859648632d+10, 2.3883787996332290397d+12,
     &         9.3614022392337710626d+02/
c
c  Coefficients for rational approximation of
c  J0(X) / (X**2 - XJ1**2), 4.0 < |X| <= 8.0
c
      data pj1/4.4176707025325087628d+03, 1.1725046279757103576d+04,
     &         1.0341910641583726701d+04,-7.2879702464464618998d+03,
     &        -1.2254078161378989535d+04,-1.8319397969392084011d+03,
     &         4.8591703355916499363d+01, 7.4321196680624245801d+02/
      data qj1/3.3307310774649071172d+02,-2.9458766545509337327d+03,
     &         1.8680990008359188352d+04,-8.4055062591169562211d+04,
     &         2.4599102262586308984d+05,-3.5783478026152301072d+05,
     &        -2.5258076240801555057d+01/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
c  XSMALL < |X| <= 3.0
c
      data py0/1.0102532948020907590d+04,-2.1287548474401797963d+06,
     &         2.0422274357376619816d+08,-8.3716255451260504098d+09,
     &         1.0723538782003176831d+11,-1.8402381979244993524d+01/
      data qy0/6.6475986689240190091d+02, 2.3889393209447253406d+05,
     &         5.5662956624278251596d+07, 8.1617187777290363573d+09,
     &         5.8873865738997033405d+11/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
c  3.0 < |X| <= 5.5
c
      data py1/-1.4566865832663635920d+04, 4.6905288611678631510d+06,
     &         -6.9590439394619619534d+08, 4.3600098638603061642d+10,
     &         -5.5107435206722644429d+11,-2.2213976967566192242d+13,
     &          1.7427031242901594547d+01/
      data qy1/ 8.3030857612070288823d+02, 4.0669982352539552018d+05,
     &          1.3960202770986831075d+08, 3.4015103849971240096d+10,
     &          5.4266824419412347550d+12, 4.3386146580707264428d+14/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
c  5.5 < |X| <= 8.0
c
      data py2/ 2.1363534169313901632d+04,-1.0085539923498211426d+07,
     &          2.1958827170518100757d+09,-1.9363051266772083678d+11,
     &         -1.2829912364088687306d+11, 6.7016641869173237784d+14,
     &         -8.0728726905150210443d+15,-1.7439661319197499338d+01/
      data qy2/ 8.7903362168128450017d+02, 5.3924739209768057030d+05,
     &          2.4727219475672302327d+08, 8.6926121104209825246d+10,
     &          2.2598377924042897629d+13, 3.9272425569640309819d+15,
     &          3.4563724628846457519d+17/
c
c  Coefficients for Hart,s approximation, 8.0 < |X|.
c
      data p0/3.4806486443249270347d+03, 2.1170523380864944322d+04,
     &        4.1345386639580765797d+04, 2.2779090197304684302d+04,
     &        8.8961548424210455236d-01, 1.5376201909008354296d+02/
      data q0/3.5028735138235608207d+03, 2.1215350561880115730d+04,
     &        4.1370412495510416640d+04, 2.2779090197304684318d+04,
     &        1.5711159858080893649d+02/
      data p1/-2.2300261666214198472d+01,-1.1183429920482737611d+02,
     &        -1.8591953644342993800d+02,-8.9226600200800094098d+01,
     &        -8.8033303048680751817d-03,-1.2441026745835638459d+00/
      data q1/1.4887231232283756582d+03, 7.2642780169211018836d+03,
     &        1.1951131543434613647d+04, 5.7105024128512061905d+03,
     &        9.0593769594993125859d+01/
c
c  Check for error conditions.
c
      ax = abs ( arg )

      if ( jint .eq. 1 .and. arg .le. zero ) then
        result = -xinf
        return
      else if ( xmax .lt. ax ) then
        result = zero
        return
      end if

      if ( eight .lt. ax ) then
        go to 800
      end if

      if ( ax .le. xsmall ) then
        if ( jint .eq. 0 ) then
          result = one
        else
          result = pi2 * ( log ( ax ) + cons )
        end if
        return
      end if
c
c  Calculate J0 for appropriate interval, preserving
c  accuracy near the zero of J0.
c
      zsq = ax * ax

      if ( ax .le. four ) then
        xnum = ( pj0(5) * zsq + pj0(6) ) * zsq + pj0(7)
        xden = zsq + qj0(5)
        do i = 1, 4
          xnum = xnum * zsq + pj0(i)
          xden = xden * zsq + qj0(i)
        end do
        prod = ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )
      else
        wsq = one - zsq / sixty4
        xnum = pj1(7) * wsq + pj1(8)
        xden = wsq + qj1(7)
        do i = 1, 6
          xnum = xnum * wsq + pj1(i)
          xden = xden * wsq + qj1(i)
        end do
        prod = ( ax + xj1 ) * ( ( ax - xj11 / two56 ) - xj12 )
      end if

      result = prod * xnum / xden

      if ( jint .eq. 0 ) then
        return
      end if
c
c  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
c  where xn is a zero of Y0.
c
      if ( ax .le. three ) then
        up = ( ax - xy01 / two56 ) - xy02
        xy = xy0
      else if ( ax .le. five5 ) then
        up = ( ax - xy11 / two56 ) - xy12
        xy = xy1
      else
        up = ( ax - xy21 / two56 ) - xy22
        xy = xy2
      end if

      down = ax + xy

      if ( abs ( up ) .lt. p17 * down ) then
        w = up / down
        wsq = w * w
        xnum = plg(1)
        xden = wsq + qlg(1)
        do i = 2, 4
          xnum = xnum * wsq + plg(i)
          xden = xden * wsq + qlg(i)
        end do
        resj = pi2 * result * w * xnum / xden
      else
        resj = pi2 * result * log ( ax / xy )
      end if
c
c  Now calculate Y0 for appropriate interval, preserving
c  accuracy near the zero of Y0.
c
      if ( ax .le. three ) then
        xnum = py0(6) * zsq + py0(1)
        xden = zsq + qy0(1)
        do i = 2, 5
          xnum = xnum * zsq + py0(i)
          xden = xden * zsq + qy0(i)
        end do
      else if ( ax .le. five5 ) then
        xnum = py1(7) * zsq + py1(1)
        xden = zsq + qy1(1)
        do i = 2, 6
          xnum = xnum * zsq + py1(i)
          xden = xden * zsq + qy1(i)
        end do
      else
        xnum = py2(8) * zsq + py2(1)
        xden = zsq + qy2(1)
        do i = 2, 7
          xnum = xnum * zsq + py2(i)
          xden = xden * zsq + qy2(i)
        end do
      end if

      result = resj + up * down * xnum / xden

      return
c
c  Calculate J0 or Y0 for 8.0 < |ARG|.
c
  800 continue

      z = eight / ax
      w = ax / twopi
      w = aint ( w ) + oneov8
      w = ( ax - w * twopi1 ) - w * twopi2
      zsq = z * z
      xnum = p0(5) * zsq + p0(6)
      xden = zsq + q0(5)
      up = p1(5) * zsq + p1(6)
      down = zsq + q1(5)

      do i = 1, 4
        xnum = xnum * zsq + p0(i)
        xden = xden * zsq + q0(i)
        up = up * zsq + p1(i)
        down = down * zsq + q1(i)
      end do

      r0 = xnum / xden
      r1 = up / down

      if ( jint .eq. 0 ) then
        result = sqrt ( pi2 / ax ) 
     &    * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
      else
        result = sqrt ( pi2 / ax ) 
     &    * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
      end if

      return
      end
      subroutine caljy1 ( arg, result, jint )

c*********************************************************************72
c
cc CALJY1 computes various J1 and Y1 Bessel functions.
c
c  Discussion:
c
c    This routine computes first-order Bessel functions of the first and
c    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
c    for Y1, and |X| <= XMAX for J1.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT = 0, ARG
c    must satisfy 
c     -XMAX < ARG < XMAX;
c    If JINT = 1, then ARG must satisfy
c      0 < ARG < XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    0, RESULT = J1(x);
c    1, RESULT = Y1(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    0, J1(x);
c    1, Y1(x);
c
      implicit none

      double precision arg
      double precision ax
      double precision down
      double precision eight
      double precision four
      double precision half
      integer i
      integer jint
      double precision pi2
      double precision pj0(7)
      double precision pj1(8)
      double precision plg(4)
      double precision prod
      double precision py0(7)
      double precision py1(9)
      double precision p0(6)
      double precision p1(6)
      double precision p17
      double precision qj0(5)
      double precision qj1(7)
      double precision qlg(4)
      double precision qy0(6)
      double precision qy1(8)
      double precision q0(6)
      double precision q1(6)
      double precision resj
      double precision result
      double precision rtpi2
      double precision r0
      double precision r1
      double precision throv8
      double precision twopi
      double precision twopi1
      double precision twopi2
      double precision two56
      double precision up
      double precision w
      double precision wsq
      double precision xden
      double precision xinf
      double precision xmax
      double precision xnum
      double precision xsmall
      double precision xj0
      double precision xj1
      double precision xj01
      double precision xj02
      double precision xj11
      double precision xj12
      double precision xy
      double precision xy0
      double precision xy01
      double precision xy02
      double precision xy1
      double precision xy11
      double precision xy12
      double precision z
      double precision zero
      double precision zsq
c
c  Mathematical constants
c
      data eight /8.0d0/
      data four /4.0d0/
      data half /0.5d0/
      data throv8 /0.375d0/
      data pi2 /6.3661977236758134308d-1/
      data p17 /1.716d-1/
      data twopi /6.2831853071795864769d+0/
      data zero /0.0d0/
      data twopi1 /6.28125d0/
      data twopi2 /1.9353071795864769253d-03/
      data two56 /256.0d+0/
      data rtpi2 /7.9788456080286535588d-1/
c
c  Machine-dependent constants
c
      data xmax /1.07d+09/
      data xsmall /9.31d-10/
      data xinf /1.7d+38/
c
c  Zeroes of Bessel functions
c
      data xj0 /3.8317059702075123156d+0/
      data xj1 /7.0155866698156187535d+0/
      data xy0 /2.1971413260310170351d+0/
      data xy1 /5.4296810407941351328d+0/
      data xj01 / 981.0d+0/
      data xj02 /-3.2527979248768438556d-04/
      data xj11 /1796.0d+0/
      data xj12 /-3.8330184381246462950d-05/
      data xy01 / 562.0d+0/
      data xy02 / 1.8288260310170351490d-03/
      data xy11 /1390.0d+0/
      data xy12 /-6.4592058648672279948d-06/
c
c  Coefficients for rational approximation to ln(x/a)
c
      data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02,
     &         -5.4989956895857911039d+02,3.5687548468071500413d+02/
      data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02,
     &         -3.3442903192607538956d+02,1.7843774234035750207d+02/
c
c  Coefficients for rational approximation of
c  J1(X) / (X * (X**2 - XJ0**2)), XSMALL < |X| <=  4.0
c
      data pj0/9.8062904098958257677d+05,-1.1548696764841276794d+08,
     &       6.6781041261492395835d+09,-1.4258509801366645672d+11,
     &      -4.4615792982775076130d+03, 1.0650724020080236441d+01,
     &      -1.0767857011487300348d-02/
      data qj0/5.9117614494174794095d+05, 2.0228375140097033958d+08,
     &       4.2091902282580133541d+10, 4.1868604460820175290d+12,
     &       1.0742272239517380498d+03/
c
c  Coefficients for rational approximation of
c  J1(X) / (X * (X**2 - XJ1**2)), 4.0 < |X| <= 8.0
c
      data pj1/4.6179191852758252280d+00,-7.1329006872560947377d+03,
     &       4.5039658105749078904d+06,-1.4437717718363239107d+09,
     &       2.3569285397217157313d+11,-1.6324168293282543629d+13,
     &       1.1357022719979468624d+14, 1.0051899717115285432d+15/
      data qj1/1.1267125065029138050d+06, 6.4872502899596389593d+08,
     &       2.7622777286244082666d+11, 8.4899346165481429307d+13,
     &       1.7128800897135812012d+16, 1.7253905888447681194d+18,
     &       1.3886978985861357615d+03/
c
c  Coefficients for rational approximation of
c  (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
c  XSMALL < |X| <=  4.0
c
      data py0/2.2157953222280260820d+05,-5.9157479997408395984d+07,
     &         7.2144548214502560419d+09,-3.7595974497819597599d+11,
     &         5.4708611716525426053d+12, 4.0535726612579544093d+13,
     &        -3.1714424660046133456d+02/
      data qy0/8.2079908168393867438d+02, 3.8136470753052572164d+05,
     &         1.2250435122182963220d+08, 2.7800352738690585613d+10,
     &         4.1272286200406461981d+12, 3.0737873921079286084d+14/
c
c  Coefficients for rational approximation of
c  (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
c  4.0 < |X| <= 8.0
c
      data py1/ 1.9153806858264202986d+06,-1.1957961912070617006d+09,
     &          3.7453673962438488783d+11,-5.9530713129741981618d+13,
     &          4.0686275289804744814d+15,-2.3638408497043134724d+16,
     &         -5.6808094574724204577d+18, 1.1514276357909013326d+19,
     &         -1.2337180442012953128d+03/
      data qy1/ 1.2855164849321609336d+03, 1.0453748201934079734d+06,
     &          6.3550318087088919566d+08, 3.0221766852960403645d+11,
     &          1.1187010065856971027d+14, 3.0837179548112881950d+16,
     &          5.6968198822857178911d+18, 5.3321844313316185697d+20/
c
c  Coefficients for Hart's approximation, 8.0 < |X|.
c
      data p0/-1.0982405543459346727d+05,-1.5235293511811373833d+06,
     &         -6.6033732483649391093d+06,-9.9422465050776411957d+06,
     &         -4.4357578167941278571d+06,-1.6116166443246101165d+03/
      data q0/-1.0726385991103820119d+05,-1.5118095066341608816d+06,
     &         -6.5853394797230870728d+06,-9.9341243899345856590d+06,
     &         -4.4357578167941278568d+06,-1.4550094401904961825d+03/
      data p1/ 1.7063754290207680021d+03, 1.8494262873223866797d+04,
     &          6.6178836581270835179d+04, 8.5145160675335701966d+04,
     &          3.3220913409857223519d+04, 3.5265133846636032186d+01/
      data q1/ 3.7890229745772202641d+04, 4.0029443582266975117d+05,
     &          1.4194606696037208929d+06, 1.8194580422439972989d+06,
     &          7.0871281941028743574d+05, 8.6383677696049909675d+02/
c
c  Check for error conditions.
c
      ax = abs ( arg )

      if ( jint .eq. 1 .and. ( arg .le. zero .or.
     &  ( arg .lt. half .and. ax * xinf .lt. pi2 ) ) ) then
        result = -xinf
        return
      else if ( xmax .lt. ax ) then
        result = zero
        return
      end if

      if ( eight .lt. ax ) then
        go to 800
      else if ( ax .le. xsmall ) then
        if ( jint .eq. 0 ) then
          result = arg * half
        else
          result = -pi2 / ax
        end if
        return
      end if
c
c  Calculate J1 for appropriate interval, preserving
c  accuracy near the zero of J1.
c
      zsq = ax * ax

      if ( ax .le. four ) then
        xnum = ( pj0(7) * zsq + pj0(6) ) * zsq + pj0(5)
        xden = zsq + qj0(5)
        do i = 1, 4
          xnum = xnum * zsq + pj0(i)
          xden = xden * zsq + qj0(i)
        end do
        prod = arg * ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )
      else
        xnum = pj1(1)
        xden = ( zsq + qj1(7) ) * zsq + qj1(1)
        do i = 2, 6
          xnum = xnum * zsq + pj1(i)
          xden = xden * zsq + qj1(i)
        end do
        xnum = xnum * ( ax - eight ) * ( ax + eight ) + pj1(7)
        xnum = xnum * ( ax - four ) * ( ax + four ) + pj1(8)
        prod = arg * ( ( ax - xj11 / two56 ) - xj12 ) * ( ax + xj1 )
      end if

      result = prod * ( xnum / xden )

      if ( jint .eq. 0 ) then
        return
      end if
c
c  Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
c  where xn is a zero of Y1.
c
      if ( ax .le. four ) then
        up = ( ax - xy01 / two56 ) - xy02
        xy = xy0
      else
        up = ( ax - xy11 / two56 ) - xy12
        xy = xy1
      end if

      down = ax + xy

      if ( abs ( up ) .lt. p17 * down ) then
        w = up / down
        wsq = w * w
        xnum = plg(1)
        xden = wsq + qlg(1)
        do i = 2, 4
          xnum = xnum * wsq + plg(i)
          xden = xden * wsq + qlg(i)
        end do
        resj = pi2 * result * w * xnum / xden
      else
        resj = pi2 * result * log ( ax / xy )
      end if
c
c  Now calculate Y1 for appropriate interval, preserving
c  accuracy near the zero of Y1.
c
      if ( ax .le. four ) then
        xnum = py0(7) * zsq + py0(1)
        xden = zsq + qy0(1)
        do i = 2, 6
          xnum = xnum * zsq + py0(i)
          xden = xden * zsq + qy0(i)
        end do
      else
        xnum = py1(9) * zsq + py1(1)
        xden = zsq + qy1(1)
        do i = 2, 8
          xnum = xnum * zsq + py1(i)
          xden = xden * zsq + qy1(i)
        end do
      end if

      result = resj + ( up * down / ax ) * xnum / xden
      return
c
c  Calculate J1 or Y1 for 8.0 < |ARG|.
c
  800 continue

      z = eight / ax
      w = aint ( ax / twopi ) + throv8
      w = ( ax - w * twopi1 ) - w * twopi2
      zsq = z * z
      xnum = p0(6)
      xden = zsq + q0(6)
      up = p1(6)
      down = zsq + q1(6)

      do i = 1, 5
        xnum = xnum * zsq + p0(i)
        xden = xden * zsq + q0(i)
        up = up * zsq + p1(i)
        down = down * zsq + q1(i)
      end do

      r0 = xnum / xden
      r1 = up / down

      if ( jint .eq. 0 ) then
        result = ( rtpi2 / sqrt ( ax ) ) 
     &    * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
      else
        result = ( rtpi2 / sqrt ( ax ) ) 
     &    * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
      end if

      if ( jint .eq. 0 .and. arg .lt. zero ) then
        result = -result
      end if

      return
      end
      function daw ( xx )

c*********************************************************************72
c
cc DAW evaluates Dawson's integral function.
c
c  Discussion:
c
c    This routine evaluates Dawson's integral,
c
c      F(x) = exp ( - x * x ) * Integral ( 0 <= t <= x ) exp ( t * t ) dt
c
c    for a real argument x.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Kathleen Paciorek, Henry Thacher,
c    Chebyshev Approximations for Dawson's Integral,
c    Mathematics of Computation,
c    Volume 24, Number 109, January 1970, pages 171-178.
c
c  Parameters:
c
c    Input, double precision XX, the argument of the function.
c
c    Output, double precision DAW, the value of the function.
c
      implicit none

      double precision daw
      double precision frac
      double precision half
      integer i
      double precision one
      double precision one225
      double precision p1(10)
      double precision p2(10)
      double precision p3(10)
      double precision p4(10)
      double precision q1(10)
      double precision q2(9)
      double precision q3(9)
      double precision q4(9)
      double precision six25
      double precision sump
      double precision sumq
      double precision two5
      double precision w2
      double precision x
      double precision xlarge
      double precision xmax
      double precision xsmall
      double precision xx
      double precision y
      double precision zero
c
c  Mathematical constants.
c
      data zero / 0.0D+00 /
      data half / 0.5D+00 /
      data one / 1.0D+00 /
      data six25 /6.25D+00 /
      data one225 /12.25d0 /
      data two5 /25.0d0/
c
c  Machine-dependent constants
c
      data xsmall /1.05d-08/
      data xlarge /9.49d+07/
      data xmax /2.24d+307/
c
c  Coefficients for R(9,9) approximation for  |x| < 2.5
c
      data p1/-2.69020398788704782410d-12, 4.18572065374337710778d-10,
     &        -1.34848304455939419963d-08, 9.28264872583444852976d-07,
     &        -1.23877783329049120592d-05, 4.07205792429155826266d-04,
     &        -2.84388121441008500446d-03, 4.70139022887204722217d-02,
     &        -1.38868086253931995101d-01, 1.00000000000000000004d+00/
      data q1/ 1.71257170854690554214d-10, 1.19266846372297253797d-08,
     &         4.32287827678631772231d-07, 1.03867633767414421898d-05,
     &         1.78910965284246249340d-04, 2.26061077235076703171d-03,
     &         2.07422774641447644725d-02, 1.32212955897210128811d-01,
     &         5.27798580412734677256d-01, 1.00000000000000000000d+00/
c
c  Coefficients for R(9,9) approximation in J-fraction form
c  for  x in [2.5, 3.5)
c
      data p2/-1.70953804700855494930d+00,-3.79258977271042880786d+01,
     &         2.61935631268825992835d+01, 1.25808703738951251885d+01,
     &        -2.27571829525075891337d+01, 4.56604250725163310122d+00,
     &        -7.33080089896402870750d+00, 4.65842087940015295573d+01,
     &        -1.73717177843672791149d+01, 5.00260183622027967838d-01/
      data q2/ 1.82180093313514478378d+00, 1.10067081034515532891d+03,
     &        -7.08465686676573000364d+00, 4.53642111102577727153d+02,
     &         4.06209742218935689922d+01, 3.02890110610122663923d+02,
     &         1.70641269745236227356d+02, 9.51190923960381458747d+02,
     &         2.06522691539642105009d-01/
c
c  Coefficients for R(9,9) approximation in J-fraction form
c  for  x in [3.5, 5.0]
c
      data p3/-4.55169503255094815112d+00,-1.86647123338493852582d+01,
     &        -7.36315669126830526754d+00,-6.68407240337696756838d+01,
     &         4.84507265081491452130d+01, 2.69790586735467649969d+01,
     &        -3.35044149820592449072d+01, 7.50964459838919612289d+00,
     &        -1.48432341823343965307d+00, 4.99999810924858824981d-01/
      data q3/ 4.47820908025971749852d+01, 9.98607198039452081913d+01,
     &         1.40238373126149385228d+01, 3.48817758822286353588d+03,
     &        -9.18871385293215873406d+00, 1.24018500009917163023d+03,
     &        -6.88024952504512254535d+01,-2.31251575385145143070d+00,
     &         2.50041492369922381761d-01/
c
c  Coefficients for R(9,9) approximation in J-fraction form
c  for 5.0 < |x|.
c
      data p4/-8.11753647558432685797d+00,-3.84043882477454453430d+01,
     &        -2.23787669028751886675d+01,-2.88301992467056105854d+01,
     &        -5.99085540418222002197d+00,-1.13867365736066102577d+01,
     &        -6.52828727526980741590d+00,-4.50002293000355585708d+00,
     &        -2.50000000088955834952d+00, 5.00000000000000488400d-01/
      data q4/ 2.69382300417238816428d+02, 5.04198958742465752861d+01,
     &         6.11539671480115846173d+01, 2.08210246935564547889d+02,
     &         1.97325365692316183531d+01,-1.22097010558934838708d+01,
     &        -6.99732735041547247161d+00,-2.49999970104184464568d+00,
     &         7.49999999999027092188d-01/

      x = xx

      if ( xlarge .lt. abs ( x ) ) then

        if ( abs ( x ) .le. xmax ) then
          daw = half / x
        else
          daw = zero
        end if

      else if ( abs ( x ) .lt. xsmall ) then

        daw = x

      else

        y = x * x
c
c  ABS(X) < 2.5.
c
        if ( y .lt. six25 ) then

          sump = p1(1)
          sumq = q1(1)
          do i = 2, 10
            sump = sump * y + p1(i)
            sumq = sumq * y + q1(i)
          end do

          daw = x * sump / sumq
c
c  2.5 <= ABS(X) < 3.5.
c
        else if ( y .lt. one225 ) then

          frac = zero
          do i = 1, 9
            frac = q2(i) / ( p2(i) + y + frac )
          end do

          daw = ( p2(10) + frac ) / x
c
c  3.5 <= ABS(X) < 5.0.
c
        else if ( y .lt. two5 ) then

          frac = zero
          do i = 1, 9
            frac = q3(i) / ( p3(i) + y + frac )
          end do

          daw = ( p3(10) + frac ) / x

        else
c
c  5.0 <= ABS(X) <= XLARGE.
c
          w2 = one / x / x

          frac = zero
          do i = 1, 9
            frac = q4(i) / ( p4(i) + y + frac )
          end do
          frac = p4(10) + frac

          daw = ( half + half * w2 * frac ) / x

        end if

      end if

      return
      end
      function dlgama ( x )

c*********************************************************************72
c
cc DLGAMA evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision DLGAMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision dlgama
      double precision eps
      double precision frtbig
      double precision four
      double precision half
      integer i
      double precision one
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision res
      double precision sqrtpi
      double precision thrhal
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c
      data one /1.0D+00/
      data half /0.5D+00/
      data twelve /12.0D+00/
      data zero /0.0D+00/
      data four /4.0D+00/
      data thrhal /1.5D+00/
      data two /2.0D+00/
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( zero .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. thrhal ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = zero
            xm1 = ( y - half ) - half
          end if

          if ( y .le. half .or. pnt68 .le. y ) then

            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - half ) - half
            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do

            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. four ) then

          xm2 = y - two
          xden = one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. twelve ) then

          xm4 = y - four
          xden = -one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = zero

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - half * corr
          res = res + y * ( corr - one )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      dlgama = res

      return
      end
      subroutine dsubn ( x, nmax, xmax, d )

c*********************************************************************72
c
cc DSUBN evaluates derivatives of Ei(X).
c
c  Discussion:
c
c    Translation of Gautschi CACM Algorithm 282 for derivatives of Ei(x).
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 282:
c    Derivatives of EXP(X)/X, COS(X)/X, and SIN(X)/X,
c    Communications of the ACM,
c    Volume 9, April 1966, page 272.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Input, integer NMAX, the number of derivatives to compute.
c
c    Input, double precision XMAX, the largest finite floating-point number,
c    as computed, for instance, by MACHAR.
c
c    Output, double precision D(0:NMAX), the value and first NMAX derivatives
c    of Ei(X).
c
      implicit none

      integer nmax

      double precision b0
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      logical bool1
      logical bool2
      double precision c0
      double precision c1
      double precision d(0:nmax)
      double precision e
      double precision en
      integer j
      integer lim
      integer mini
      integer n
      integer n0
      integer n1
      double precision one
      double precision p
      double precision q
      double precision t
      double precision ten
      double precision two
      double precision x
      double precision xmax
      double precision x1
      double precision z
      double precision zero

      data zero /0.0d0/
      data one /1.0d0/
      data two /2.0d0/
      data ten /10.0d0/
      data c0 /2.7183d0/
      data c1 /4.67452d1/
      data b0 /5.7941d-5/
      data b1 /-1.76148d-3/
      data b2 /2.08645d-2/
      data b3 /-1.29013d-1/
      data b4 /8.5777d-1/
      data b5 /1.0125d0/
      data b6 /7.75d-1/

      x1 = abs ( x )
      n0 = int ( x1 )
      e = exp ( x )
      d(0) = e / x
      bool1 = ( x .lt. zero ) .or. ( x1 .le. two )
      bool2 = n0 .lt. nmax
      mini = min ( n0, nmax )

      if ( bool1 ) then
        lim = nmax
      else
        lim = mini
      end if

      n = 1
      en = one

   50 continue

      d(n) = ( e - en * d(n-1) ) / x
      n = n + 1
      en = en + one

      if ( x1 .lt. one ) then

        if ( abs ( d(n-1) ) .lt. abs ( xmax * x / en ) .and. 
     &    n .le. lim ) then
          go to 50
        end if

      else

        if ( abs ( d(n-1) / x ) .lt. xmax / en .and. n .le. lim ) then
          go to 50
        end if

      end if

      do j = n, lim
        d(n) = zero
      end do

      if ( .not. bool1 .and. bool2 ) then

        t = ( x1 + c1 ) / ( c0 * x1 )

        if ( t .lt. ten ) then
          t = (((( 
     &            b0 
     &      * t + b1 )
     &      * t + b2 )
     &      * t + b3 )
     &      * t + b4 )
     &      * t + b5
        else
          z = log ( t ) - b6
          p = ( b6 - log ( z ) ) / ( one + z )
          p = one / ( one + p )
          t = t * p / z
        end if

        n1 = c0 * x1 * t - one
        if ( n1 .lt. nmax ) then
          n1 = nmax
        end if

        q = one / x
        en = one
        do n = 1, n1+1
          q = -en * q / x
          en = en + one
        end do

        do n = n1, n0+1, -1

          en = en - one
          q = ( e - x * q ) / en

          if ( n .le. nmax ) then
            d(n) = q
          end if

        end do

      end if

      return
      end
      function ei ( x )

c*********************************************************************72
c
cc EI evaluates the exponential integral Ei(X).
c
c  Discussion:
c
c    This routine computes approximate values for the
c    exponential integral Ei(x), where x is real.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision EI, the value of the function.
c
      implicit none

      double precision ei
      integer jint
      double precision result
      double precision x

      jint = 1
      call calcei ( x, result, jint )
      ei = result

      return
      end
      function eone ( x )

c*********************************************************************72
c
cc EONE evaluates the exponential integral E1(X).
c
c  Discussion:
c
c    This routine computes approximate values for the
c    exponential integral E1(x), where x is real.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision EONE, the value of the function.
c
      implicit none

      double precision eone
      integer jint
      double precision result
      double precision x

      jint = 2
      call calcei ( x, result, jint )
      eone = result

      return
      end
      function expei ( x )

c*********************************************************************72
c
cc EXPEI evaluates the scaled exponential integral exp(-X) * Ei(X).
c
c  Discussion:
c
c    This routine computes approximate values for the
c    function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
c    integral, and  x  is real.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision EXPEI, the value of the function.
c
      implicit none

      double precision expei
      integer jint
      double precision result
      double precision x

      jint = 3
      call calcei ( x, result, jint )
      expei = result

      return
      end
      subroutine machar ( ibeta, it, irnd, ngrd, machep, negep, iexp,
     &  minexp, maxexp, eps, epsneg, xmin, xmax )

c*********************************************************************72
c
cc MACHAR determines various machine arithmetic parameters.
c
c  Discussion:
c
c    This routine is intended to determine the parameters
c    of the floating-point arithmetic system specified below.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, William Waite,
c    Software Manual for the Elementary Functions,
c    Prentice Hall, 1980,
c    ISBN: 0138220646,
c    LC: QA331.C635.
c
c    Michael Malcolm,
c    Algorithms to Reveal Properties of Floating Point Arithmetic,
c    Communications of the ACM,
c    Volume 15, Number 11, November 1972, pages 949-951.
c
c  Parameters:
c
c    Output, integer IBETA, the radix for the floating-point representation
c
c    Output, integer IT, the number of base IBETA digits in the floating-point
c    significand.
c
c    Output, integer IRND.
c    * 0, if floating-point addition chops;
c    * 1, if floating-point addition rounds, but not in the IEEE style;
c    * 2, if floating-point addition rounds in the IEEE style;
c    * 3, if floating-point addition chops, and there is partial underflow;
c    * 4, if floating-point addition rounds, but not in the IEEE style, 
c      and there is partial underflow
c    * 5, if floating-point addition rounds in the IEEE style, and there 
c      is partial underflow
c
c    Output, integer NGRD, the number of guard digits for multiplication with
c    truncating arithmetic.  It is
c    * 0, if floating-point arithmetic rounds, or if it truncates and only
c      IT base IBETA digits participate in the post-normalization shift of the
c      floating-point significand in multiplication;
c    * 1, if floating-point arithmetic truncates and more than IT base IBETA
c      digits participate in the post-normalization shift of the floating-point
c      significand in multiplication.
c
c    Output, integer MACHEP. the largest negative integer such that
c    1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that MACHEP is bounded below 
c    by  -(IT+3).
c
c    Output, integer NEGEPS. the largest negative integer such that
c    1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that NEGEPS is bounded below 
c    by  -(IT+3).
c
c    Output, integer IEXP, the number of bits (decimal places if IBETA = 10)
c    reserved for the representation of the exponent, including the bias or 
c    sign, of a floating-point number.
c
c    Output, integer MINEXP, largest magnitude negative integer such that
c    FLOAT(IBETA)**MINEXP is positive and normalized.
c
c    Output, integer MAXEXP, the smallest positive power of BETA that 
c    overflows.
c
c    Output, double precision EPS, the value FLOAT(IBETA)**MACHEP.
c
c    Output, double precision EPSNEG, the value FLOAT(IBETA)**NEGEPS.
c
c    Output, double precision XMIN, the smallest non-vanishing normalized 
c    floating-point power of the radix, that is, XMIN = FLOAT(IBETA)**MINEXP.
c
c    Output, double precision XMAX, the largest finite floating-point number.
c    In particular XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP.
c    On some machines  XMAX  will be only the second, or perhaps third, 
c    largest number, being too small by 1 or 2 units in the last digit of
c    the significand.
c
      implicit none

      double precision a
      double precision b
      double precision beta
      double precision betain
      double precision betah
      double precision eps
      double precision epsneg
      integer i
      integer ibeta
      integer iexp
      integer irnd
      integer it
      integer itemp
      integer iz
      integer j
      integer k
      integer machep
      integer maxexp
      integer minexp
      integer mx
      integer negep
      integer ngrd
      integer nxres
      double precision one
      double precision t
      double precision temp
      double precision tempa
      double precision temp1
      double precision two
      double precision xmax
      double precision xmin
      double precision y
      double precision z
      double precision zero

      one = dble ( 1 )
      two = one + one
      zero = one - one
c
c  Determine IBETA, BETA ala Malcolm.
c
      a = one
   10 continue
      a = a + a
      temp = a + one
      temp1 = temp - a
      if ( temp1 - one .eq. zero ) then
        go to 10
      end if

      b = one
   20 continue
      b = b + b
      temp = a + b
      itemp = int ( temp - a )
      if ( itemp .eq. 0 ) then
        go to 20
      end if

      ibeta = itemp
      beta = dble ( ibeta )
c
c  Determine IT, IRND.
c
      it = 0
      b = one
  100 continue
      it = it + 1
      b = b * beta
      temp = b + one
      temp1 = temp - b
      if ( temp1 - one .eq. zero ) then
        go to 100
      end if

      irnd = 0
      betah = beta / two
      temp = a + betah

      if ( temp - a .ne. zero ) then
        irnd = 1
      end if

      tempa = a + beta
      temp = tempa + betah

      if ( irnd .eq. 0 .and. temp - tempa .ne. zero ) then
        irnd = 2
      end if
c
c  Determine NEGEP, EPSNEG.
c
      negep = it + 3
      betain = one / beta
      a = one
      do i = 1, negep
         a = a * betain
      end do

      b = a
  210 continue
      temp = one - a
      if ( temp - one .eq. zero ) then
        a = a * beta
        negep = negep - 1
        go to 210
      end if

      negep = -negep
      epsneg = a
c
c  Determine MACHEP, EPS.
c
      machep = -it - 3
      a = b
  300 continue
      temp = one + a
      if ( temp - one .eq. zero ) then
        a = a * beta
        machep = machep + 1
        go to 300
      end if

      eps = a
c
c  Determine NGRD.
c
      ngrd = 0
      temp = one + eps

      if ( irnd .eq. 0 .and. temp * one - one .ne. zero ) then
        ngrd = 1
      end if
c
c  Determine IEXP, MINEXP, XMIN.
c
c  Loop to determine largest I and K = 2**I such that
c    (1/BETA) ** (2**(I))
c  does not underflow.
c  Exit from loop is signaled by an underflow.
c
      i = 0
      k = 1
      z = betain
      t = one + eps
      nxres = 0
  400 continue
      y = z
      z = y * y
c
c  Check for underflow here.
c
      a = z * one
      temp = z * t

      if ( a + a .eq. zero .or. abs ( z ) .ge. y ) then
        go to 410
      end if

      temp1 = temp * betain
      if ( temp1 * beta .eq. z ) then
        go to 410
      end if

      i = i + 1
      k = k + k
      go to 400

  410 continue
      if ( ibeta .eq. 10 ) then
        go to 420
      end if

      iexp = i + 1
      mx = k + k
      go to 450
c
c  This segment is for decimal machines only.
c
  420 continue

      iexp = 2

      iz = ibeta

  430 continue

      if ( iz .le. k ) then
        iz = iz * ibeta
        iexp = iexp + 1
        go to 430
      end if

      mx = iz + iz - 1
c
c  Loop to determine MINEXP, XMIN.
c  Exit from loop is signaled by an underflow.
c
  450 continue

      xmin = y
      y = y * betain
c
c  Check for underflow here.
c
      a = y * one
      temp = y * t
      if ( a + a .eq. zero .or. abs ( y ) .ge. xmin ) then
        go to 460
      end if

      k = k + 1
      temp1 = temp * betain
      if ( temp1 * beta .ne. y .or. temp .eq. y ) then
        go to 450
      else
        nxres = 3
        xmin = y
      end if

  460 continue

      minexp = -k
c
c  Determine MAXEXP, XMAX.
c
      if ( k + k - 3 .lt. mx .or. ibeta .eq. 10 ) then
        go to 500
      end if

      mx = mx + mx
      iexp = iexp + 1

  500 continue

      maxexp = mx + minexp
c
c  Adjust IRND to reflect partial underflow.
c
      irnd = irnd + nxres
c
c  Adjust for IEEE-style machines.
c
      if ( irnd .ge. 2 ) then
        maxexp = maxexp - 2
      end if
c
c  Adjust for machines with implicit leading bit in binary
c  significand, and machines with radix point at extreme
c  right of significand.
c
      i = maxexp + minexp

      if ( ibeta .eq. 2 .and. i .eq. 0 ) then
        maxexp = maxexp - 1
      end if

      if ( 20 .lt. i ) then
        maxexp = maxexp - 1
      end if

      if ( a .ne. y ) then
        maxexp = maxexp - 2
      end if

      xmax = one - epsneg

      if ( xmax * one .ne. xmax ) then
        xmax = one - beta * epsneg
      end if

      xmax = xmax / ( beta * beta * beta * xmin )
      i = maxexp + minexp + 3

      do j = 1, i
        if ( ibeta .eq. 2 ) then
          xmax = xmax + xmax
        else
          xmax = xmax * beta
        end if
      end do

      return
      end
      function r8_erf ( x )

c*********************************************************************72
c
cc R8_ERF evaluates the error function.
c
c  Discussion:
c
c    This routine computes approximate values for erf(x).
c
c    This routine was renamed from "DERF" to "R8_ERF" to avoid being
c    overshadowed by routines of the same name supplied by some compilers.
c
c    See comments heading CALERF.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_ERF, the value of the function.
c
      implicit none

      integer jint
      double precision r8_erf
      double precision result
      double precision x

      jint = 0
      call calerf ( x, result, jint )
      r8_erf = result

      return
      end
      function r8_erfc ( x )

c*********************************************************************72
c
cc R8_ERFC evaluates the complementary error function.
c
c  Discussion:
c
c    This routine computes approximate values for erfc(x).
c
c    This routine was renamed from "DERFC" to "R8_ERFC" to avoid being
c    overshadowed by routines of the same name supplied by some compilers.
c
c    See comments heading CALERF.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_ERFC, the value of the function.
c
      implicit none

      integer jint
      double precision r8_erfc
      double precision result
      double precision x

      jint = 1
      call calerf ( x, result, jint )
      r8_erfc = result

      return
      end
      function r8_erfcx ( x )

c*********************************************************************72
c
cc R8_ERFCX evaluates the exponentially scaled complementary error function.
c
c  Discussion:
c
c    This routine computes approximate values for exp(x*x) * erfc(x).
c
c    This routine was renamed from "DERFCX" to "R8_ERFCX" to correspond to
c    the renamings of "DERF" and "DERFC".
c
c    See comments heading CALERF.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_ERFCX, the value of the function.
c
      implicit none

      integer jint
      double precision r8_erfcx
      double precision result
      double precision x

      jint = 2
      call calerf ( x, result, jint )
      r8_erfcx = result

      return
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates Gamma(X) for a real argument.
c
c  Discussion:
c
c    This function was originally named DGAMMA.
c
c    However, a number of Fortran compilers now include a library 
c    function of this name.  To avoid conflicts, this function was
c    renamed R8_GAMMA.
c
c    This routine calculates the GAMMA function for a real argument X.
c    Computation is based on an algorithm outlined in reference 1.
c    The program uses rational functions that approximate the GAMMA
c    function to at least 20 significant decimal digits.  Coefficients
c    for the approximation over the interval (1,2) are unpublished.
c    Those for the approximation for 12 <= X are from reference 2.
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Reference:
c
c    William Cody,
c    An Overview of Software Development for Special Functions,
c    in Numerical Analysis Dundee, 1975,
c    edited by GA Watson,
c    Lecture Notes in Mathematics 506,
c    Springer, 1976.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      double precision half
      integer i
      integer n
      double precision one
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
      double precision zero
c
c  Mathematical constants
c
      data one /1.0D+00 /
      data half /0.5D+00/
      data twelve /12.0D+00/
      data two /2.0D+00 /
      data zero /0.0D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
      data pi /3.1415926535897932384626434D+00/
c
c  Machine dependent parameters
c
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (1,2).
c
      data p/
     & -1.71618513886549492533811d+00,
     &  2.47656508055759199108314d+01,
     & -3.79804256470945635097577d+02,
     &  6.29331155312818442661052d+02,
     &  8.66966202790413211295064d+02,
     & -3.14512729688483675254357d+04,
     & -3.61444134186911729807069d+04,
     &  6.64561438202405440627855d+04/

      data q/
     & -3.08402300119738975254353d+01,
     &  3.15350626979604161529144d+02,
     & -1.01515636749021914166146d+03,
     & -3.10777167157231109440444d+03,
     &  2.25381184209801510330112d+04,
     &  4.75584627752788110767815d+03,
     & -1.34659959864969306392456d+05,
     & -1.15132259675553483497211d+05/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     & -1.910444077728D-03,
     &  8.4171387781295D-04,
     & -5.952379913043012D-04,
     &  7.93650793500350248D-04,
     & -2.777777777777681622553D-03,
     &  8.333333333333333331554247D-02,
     &  5.7083835261D-03/

      parity = .false.
      fact = one
      n = 0
      y = x
c
c  Argument is negative.
c
      if ( y .le. zero ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. zero ) then

          if ( y1 .ne. aint ( y1 * half ) * two ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + one

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Argument is positive.
c
      if ( y .lt. eps ) then
c
c  Argument < EPS.
c
        if ( xminin .le. y ) then
          res = one / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. twelve ) then

        y1 = y
c
c  0.0 < argument < 1.0.
c
        if ( y .lt. one ) then

          z = y
          y = y + one
c
c  1.0 < argument < 12.0.
c  Reduce argument if necessary.
c
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - one

        end if
c
c  Evaluate approximation for 1.0 < argument < 2.0.
c
        xnum = zero
        xden = one
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + one
c
c  Adjust result for case  0.0 < argument < 1.0.
c
        if ( y1 .lt. y ) then

          res = res / y1
c
c  Adjust result for case 2.0 < argument < 12.0.
c
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + one
          end do

        end if

      else
c
c  Evaluate for 12.0 <= argument.
c
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Final adjustments and return.
c
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. one ) then
        res = fact / res
      end if

      r8_gamma = res

      return
      end
      function r8_psi ( xx )

c*********************************************************************72
c
cc R8_PSI evaluates the function Psi(X).
c
c  Discussion:
c
c    This routine evaluates the logarithmic derivative of the
c    GAMMA function,
c
c      PSI(X) = d/dX (GAMMA(X)) / GAMMA(X) 
c             = d/dX LN ( GAMMA(X) )
c
c    for real X, where either
c
c      -XMAX1 < X < -XMIN  and X is not a negative integer), 
c
c    or
c
c      XMIN < X.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Anthony Strecok, Henry Thacher,
c    Chebyshev Approximations for the Psi Function,
c    Mathematics of Computation,
c    Volume 27, Number 121, January 1973, pages 123-127.
c
c  Parameters:
c
c    Input, double precision XX, the argument of the function.
c
c    Output, double precision R8_PSI, the value of the function.
c
      implicit none

      double precision aug
      double precision den
      double precision four
      double precision fourth
      double precision half
      integer i
      integer n
      integer nq
      double precision one
      double precision p1(9)
      double precision p2(7)
      double precision piov4
      double precision q1(8)
      double precision q2(6)
      double precision r8_psi
      double precision sgn
      double precision three
      double precision xlarge
      double precision upper
      double precision w
      double precision x
      double precision xinf
      double precision xmax1
      double precision xmin1
      double precision xsmall
      double precision x01
      double precision x01d
      double precision x02
      double precision xx
      double precision z
      double precision zero
c
c  Mathematical constants.  PIOV4 = pi / 4
c
      data zero /0.0d0/
      data fourth / 0.25d0/ 
      data half / 0.5d0 /
      data one / 1.0d0 /
      data three /3.0d0/
      data four /4.0d0/
      data piov4 /7.8539816339744830962d-01/
c
c  Machine-dependent constants
c
      data xinf /1.70d+38/
      data xmin1 /5.89d-39/
      data xmax1 /3.60d+16/
      data xsmall /2.05d-09/
      data xlarge /2.04d+15/
c
c  Zero of psi(x)
c
      data x01 /187.0d0/
      data x01d /128.0d0/
      data x02 /6.9464496836234126266d-04/
c
c  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
c
      data p1/4.5104681245762934160d-03,5.4932855833000385356d+00,
     &        3.7646693175929276856d+02,7.9525490849151998065d+03,
     &        7.1451595818951933210d+04,3.0655976301987365674d+05,
     &        6.3606997788964458797d+05,5.8041312783537569993d+05,
     &        1.6585695029761022321d+05/
      data q1/9.6141654774222358525d+01,2.6287715790581193330d+03,
     &        2.9862497022250277920d+04,1.6206566091533671639d+05,
     &        4.3487880712768329037d+05,5.4256384537269993733d+05,
     &        2.4242185002017985252d+05,6.4155223783576225996d-08/
c
c  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
c  for 3.0 < x.
c
      data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01,
     &        -1.9784554148719218667d+01,-8.8100958828312219821d+00,
     &        -1.4479614616899842986d+00,-7.3689600332394549911d-02,
     &        -6.5135387732718171306d-21/
      data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02,
     &         2.4736979003315290057d+02, 1.0742543875702278326d+02,
     &         1.7463965060678569906d+01, 8.8427520398873480342d-01/

      x = xx
      w = abs ( x )
      aug = zero
c
c  Check for valid arguments, then branch to appropriate algorithm.
c
      if ( -x .ge. xmax1 .or. w .lt. xmin1 ) then
        r8_psi = xinf
        if ( zero .lt. x ) then
          r8_psi = -xinf
        end if
        return
      end if

      if ( x .ge. half ) then
        go to 200
c
c  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
c  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
c
      else if ( w .le. xsmall ) then
        aug = -one / x
        go to 150
      end if
c
c  Argument reduction for cotangent.
c
  100 continue

      if ( x .lt. zero ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - aint ( w )
      nq = int ( w * four )
      w = four * ( w - dble ( nq ) * fourth )
c
c  W is now related to the fractional part of 4.0 * X.
c  Adjust argument to correspond to values in the first
c  quadrant and determine the sign.
c
      n = nq / 2

      if ( n + n .ne. nq ) then
        w = one - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) .ne. 0 ) then
        sgn = - sgn
      end if
c
c  Determine the final value for  -pi * cotan(pi*x).
c
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) .eq. 0 ) then
c
c  Check for singularity.
c
        if ( z .eq. zero ) then
          r8_psi = xinf
          if ( zero .lt. x ) then
            r8_psi = -xinf
          end if
          return
        end if

        aug = sgn * ( four / tan ( z ) )

      else
        aug = sgn * ( four * tan ( z ) )
      end if

  150 continue

      x = one - x

  200 continue
c
c  0.5 <= X <= 3.0.
c
      if ( x .le. three ) then

        den = x
        upper = p1(1) * x
        do i = 1, 7
          den = ( den + q1(i) ) * x
          upper = ( upper + p1(i+1) ) * x
        end do
        den = ( upper + p1(9) ) / ( den + q1(8) )
        x = ( x - x01 / x01d ) - x02
        r8_psi = den * x + aug
        return

      end if
c
c  3.0 < X.
c
      if ( x .lt. xlarge ) then
        w = one / ( x * x )
        den = w
        upper = p2(1) * w
        do i = 1, 5
          den = ( den + q2(i) ) * w
          upper = ( upper + p2(i+1) ) * w
        end do
        aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug
      end if

      r8_psi = aug + log ( x )

      return
      end
      function ren ( k )

c*********************************************************************72
c
cc REN is a random number generator.
c
c  Discussion:
c
c    This routine is intended for use on computers with
c    fixed point wordlength of at least 29 bits.  It is
c    best if the floating-point significand has at most
c    29 bits.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    Malcolm Pike, David Hill,
c    Algorithm 266:
c    Pseudo-Random Numbers,
c    Communications of the ACM,
c    Volume 8, Number 10, October 1965, page 605.
c
c  Parameters:
c
c    Input, integer K, a seed for the random number generator. 
c    (Not actually used.)
c
c    Output, double precision REN, a pseudorandom value.
c
      implicit none

      double precision c1
      parameter ( c1 = 2796203.0D+00 )
      double precision c2
      parameter ( c2 = 1.0D-06 )
      double precision c3
      parameter ( c3 = 1.0D-12 )
      integer iy
      integer k
      double precision ren

      save iy

      data iy / 100001 /

      iy = iy * 125
      iy = iy - ( iy / 2796203 ) * 2796203
      ren = dble ( iy ) / c1 * ( 1.0D+00 + c2 + c3 )

      return
      end
      subroutine ribesl ( x, alpha, nb, ize, b, ncalc )

c*********************************************************************72
c
cc RIBESL calculates I Bessel function with non-integer orders.
c
c  Discussion:
c
c    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
c    for non-negative argument X, and non-negative order N+ALPHA,
c    with or without exponential scaling.
c
c    This program is based on a program written by David 
c    Sookne that computes values of the Bessel functions J or
c    I of real argument and integer order.  Modifications include
c    the restriction of the computation to the I Bessel function
c    of non-negative real argument, the extension of the computation
c    to arbitrary positive order, the inclusion of optional
c    exponential scaling, and the elimination of most underflow.
c
c    In case of an error,  NCALC .NE. NB, and not all I's are
c    calculated to the desired accuracy.
c
c    NCALC .LT. 0:  An argument is out of range. For example,
c    NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG.
c    In this case, the B-vector is not calculated, and NCALC is
c    set to MIN0(NB,0)-1 so that NCALC .NE. NB.
c
c    NB .GT. NCALC .GT. 0: Not all requested function values could
c    be calculated accurately.  This usually occurs because NB is
c    much larger than ABS(X).  In this case, B(N) is calculated
c    to the desired accuracy for N .LE. NCALC, but precision
c    is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
c    for N .GT. NCALC (because it is too small to be represented),
c    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
c    significant figures of B(N) can be trusted.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    Frank Olver, David Sookne,
c    A Note on Backward Recurrence Algorithms,
c    Mathematics of Computation,
c    Volume 26, 1972, pages 941-947.
c
c    David Sookne,
c    Bessel Functions of Real Argument and Integer Order,
c    NBS Journal of Research B,
c    Volume 77B, 1973, pages 125-132.
c
c    William Cody,
c    Algorithm 597:
c    Sequence of Modified Bessel Functions of the First Kind,
c    ACM Transactions of Mathematical Software,
c    Volume 9, Number 2, June 1983, pages 242-245.
c
c  Parameters:
c
c    Input, double precision X, the argument for which the functions
c    are to be calculated.
c
c    Input, double precision ALPHA,the fractional part of the order 
c    for which the functions are to be calculated.
c    0 <= ALPHA < 1.0.
c
c    Input, integer NB, the number of functions to be calculated.
c    The first function calculated is of order ALPHA, and the
c    last is of order (NB - 1 + ALPHA).  1 <= NB.
c
c    Input, integer IZE, scaling option.
c    1, unscaled I's are to calculated,
c    2, exponentially scaled I's are to be calculated.
c
c    Output, double precision B(NB), the values of the functions 
c    I(ALPHA,X) through I(NB-1+ALPHA,X), with scaling if requested.
c
c    Output, integer NCALC, error indicator.  If NCALC = NB, then all the 
c    requested values were calculated to the desired accuracy.  
c
c  Local Parameeters:
c
c    BETA, radix for the floating-point system
c
c    MINEXP, the smallest representable power of beta
c
c    MAXEXP, the smallest power of BETA that overflows
c
c    IT, the number of bits in the mantissa of a working precision
c    variable.
c
c    NSIG, the decimal significance desired.  Should be set to
c    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
c    in decreased accuracy while setting NSIG higher will
c    increase CPU time without increasing accuracy.  The
c    truncation error is limited to a relative error of
c    T=.5*10**(-NSIG).
c
c    ENTEN = 10.0**K, where K is the largest integer such that
c    ENTEN is machine-representable in working precision
c
c    ENSIG = 10.0**NSIG
c
c    RTNSIG = 10.0**(-K) for the smallest integer K such that
c    NSIG/4 <= K.
c
c    ENMTEN, the smallest ABS(X) such that X/4 does not underflow
c
c    XLARGE, the upper limit on the magnitude of X when IZE=2.  Bear
c    in mind that if ABS(X)=N, then at least N iterations
c    of the backward recursion will be executed.  The value
c    of 10000.0 is used on every machine.
c
c    EXPARG, the largest working precision argument that the library
c    EXP routine can handle and upper limit on the
c    magnitude of X when IZE=1; approximately LOG(BETA**MAXEXP).
c
      implicit none

      integer nb

      double precision alpha
      double precision b(nb)
      double precision const
      double precision em
      double precision empal
      double precision emp2al
      double precision en
      double precision enmten
      double precision ensig
      double precision enten
      double precision exparg
      double precision half
      double precision halfx
      integer ize
      integer k
      integer l
      integer magx
      integer n
      integer nbmx
      integer ncalc
      integer nend
      integer nsig
      integer nstart
      double precision one
      double precision p
      double precision plast
      double precision pold
      double precision psave
      double precision psavel
      double precision r8_gamma
      double precision rtnsig
      double precision sum
      double precision tempa
      double precision tempb
      double precision tempc
      double precision test
      double precision tover
      double precision two
      double precision x
      double precision xlarge
      double precision zero

      data const / 1.585D+00 /
      data enmten / 8.9D-308 /
      data ensig / 1.0D+16 /
      data enten / 1.0D+308 /
      data exparg / 709.0D+00 /
      data half / 0.5D+00 /
      data nsig / 16 /
      data one / 1.0D+00 /
      data rtnsig / 1.0D-04 /
      data two / 2.0D+00 /
      data xlarge / 1.0D+04 /
      data zero / 0.0D+00 /
c
c  Check the input.
c
      if ( nb .le. 0 ) then
        ncalc = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RIBESL - Fatal error!'
        write ( *, '(a)' ) '  NB <= 0.'
        return
      end if

      if ( x .lt. 0.0D+00 ) then
        ncalc = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RIBESL - Fatal error!'
        write ( *, '(a)' ) '  X < 0.0.'
        return
      end if

      if ( alpha .lt. 0.0D+00 ) then
        ncalc = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RIBESL - Fatal error!'
        write ( *, '(a)' ) '  ALPHA < 0.'
        return
      end if

      if ( 1.0D+00 < alpha ) then
        ncalc = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RIBESL - Fatal error!'
        write ( *, '(a)' ) '  1 < ALPHA.'
        return
      end if
c
c  Check for X, NB, OR IZE out of range.
c
      if (
     &    ( ize .eq. 1 .and. x .le. exparg ) .or.
     &    ( ize .eq. 2 .and. x .le. xlarge ) ) then
c
c  Use 2-term ascending series for small X.
c
        ncalc = nb
        magx = int ( x )
c
c  Initialize the forward sweep, the P-sequence of Olver.
c
        if ( rtnsig .le. x ) then

          nbmx = nb - magx
          n = magx + 1
          en = dble ( n + n ) + ( alpha + alpha )
          plast = one
          p = en / x
c
c  Calculate general significance test.
c
          test = ensig + ensig

          if ( 5 * nsig .lt. 2 * magx ) then
            test = sqrt ( test * p )
          else
            test = test / const**magx
          end if

          if ( 3 .le. nbmx ) then
c
c  Calculate P-sequence until N = NB-1.  Check for possible overflow.
c
            tover = enten / ensig
            nstart = magx + 2
            nend = nb - 1

            do k = nstart, nend

              n = k
              en = en + two
              pold = plast
              plast = p
              p = en * plast / x + pold
c
c  To avoid overflow, divide P-sequence by TOVER.  Calculate
c  P-sequence until 1 < ABS(P).
c
              if ( tover .lt. p ) then

                tover = enten
                p = p / tover
                plast = plast / tover
                psave = p
                psavel = plast
                nstart = n + 1

   60           continue

                n = n + 1
                en = en + two
                pold = plast
                plast = p
                p = en * plast / x + pold

                if ( p .le. one ) then
                  go to 60
                end if

                tempb = en / x
c
c  Calculate backward test, and find NCALC, the highest N
c  such that the test is passed.
c
                test = pold * plast / ensig
                test = test * ( half - half / ( tempb * tempb ) )
                p = plast * tover
                n = n - 1
                en = en - two
                nend = min0 ( nb, n )

                do l = nstart, nend

                  ncalc = l
                  pold = psavel
                  psavel = psave
                  psave = en * psavel / x + pold

                  if ( test .lt. psave * psavel ) then
                    go to 90
                  end if

                end do

                ncalc = nend + 1

   90           continue

                ncalc = ncalc - 1
                go to 120

              end if

            end do

            n = nend
            en = dble ( n + n ) + ( alpha + alpha )
c
c  Calculate special significance test for 2 < NBMX.
c
            test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

          end if
c
c  Calculate P-sequence until significance test passed.
c
  110     continue

          n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x + pold
          if ( p .lt. test ) then
            go to 110
          end if
c
c  Initialize the backward recursion and the normalization sum.
c
  120     continue

          n = n + 1
          en = en + two
          tempb = zero
          tempa = one / p
          em = dble ( n ) - one
          empal = em + alpha
          emp2al = ( em - one ) + ( alpha + alpha )
          sum = tempa * empal * emp2al / em
          nend = n - nb
c
c  N < NB, so store B(N) and set higher orders to zero.
c
          if ( nend .lt. 0 ) then

            b(n) = tempa
            nend = -nend
            do l = 1, nend
              b(n+l) = zero
            end do

          else
c
c  Recur backward via difference equation, calculating (but
c  not storing) B(N), until N = NB.
c
            if ( 0 .lt. nend ) then

              do l = 1, nend

                n = n - 1
                en = en - two
                tempc = tempb
                tempb = tempa
                tempa = ( en * tempb ) / x + tempc
                em = em - one
                emp2al = emp2al - one

                if ( n .eq. 1 ) then
                  go to 150
                end if

                if ( n .eq. 2 ) then
                  emp2al = one
                end if

                empal = empal - one
                sum = ( sum + tempa * empal ) * emp2al / em

              end do

            end if
c
c  Store B(NB).
c
150         continue

            b(n) = tempa

            if ( nb .le. 1 ) then
              sum = ( sum + sum ) + tempa
              go to 230
            end if
c
c  Calculate and Store B(NB-1).
c
            n = n - 1
            en = en - two
            b(n)  = ( en * tempa ) / x + tempb

            if ( n .eq. 1 ) then
              go to 220
            end if

            em = em - one
            emp2al = emp2al - one

            if ( n .eq. 2 ) then
              emp2al = one
            end if

            empal = empal - one
            sum = ( sum + b(n) * empal ) * emp2al / em

          end if

          nend = n - 2
c
c  Calculate via difference equation and store B(N), until N = 2.
c
          if ( 0 .lt. nend ) then

            do l = 1, nend
              n = n - 1
              en = en - two
              b(n) = ( en * b(n+1) ) / x + b(n+2)
              em = em - one
              emp2al = emp2al - one
              if ( n .eq. 2 ) then
                emp2al = one
              end if
              empal = empal - one
              sum = ( sum + b(n) * empal ) * emp2al / em
            end do

          end if
c
c  Calculate B(1).
c
          b(1) = two * empal * b(2) / x + b(3)

  220     continue

          sum = ( sum + sum ) + b(1)
c
c  Normalize.  Divide all B(N) by sum.
c
  230     continue

          if ( alpha .ne. zero ) then
            sum = sum * r8_gamma ( one + alpha ) 
     &      * ( x * half )**(-alpha)
          end if

          if ( ize .eq. 1 ) then
            sum = sum * exp ( -x )
          end if

          tempa = enmten

          if ( one .lt. sum ) then
            tempa = tempa * sum
          end if

          do n = 1, nb
            if ( b(n) .lt. tempa ) then
              b(n) = zero
            end if
            b(n) = b(n) / sum
          end do

          return
c
c  Two-term ascending series for small X.
c
        else

          tempa = one
          empal = one + alpha
          halfx = zero

          if ( enmten .lt. x ) then
            halfx = half * x
          end if

          if ( alpha .ne. zero ) then
            tempa = halfx**alpha / r8_gamma ( empal )
          end if

          if ( ize .eq. 2 ) then
            tempa = tempa * exp ( -x )
          end if

          tempb = zero

          if ( one .lt. x + one ) then
            tempb = halfx * halfx
          end if

          b(1) = tempa + tempa * tempb / empal

          if ( x .ne. zero .and. b(1) .eq. zero ) then
            ncalc = 0
          end if

          if ( 1 .lt. nb ) then

            if ( x .eq. zero ) then

              do n = 2, nb
                b(n) = zero
              end do

            else
c
c  Calculate higher-order functions.
c
              tempc = halfx
              tover = ( enmten + enmten ) / x

              if ( tempb .ne. zero ) then
                tover = enmten / tempb
              end if

              do n = 2, nb

                tempa = tempa / empal
                empal = empal + one
                tempa = tempa * tempc

                if ( tempa .le. tover * empal ) then
                  tempa = zero
                end if

                b(n) = tempa + tempa * tempb / empal

                if ( b(n) .eq. zero .and. n .lt. ncalc ) then
                  ncalc = n - 1
                end if

              end do
            end if
          end if
        end if

      else

        ncalc = min ( nb, 0 ) - 1

      end if

      return
      end
      subroutine rjbesl ( x, alpha, nb, b, ncalc )

c*********************************************************************72
c
cc RJBESL calculates J Bessel function with non-integer orders.
c
c  Discussion:
c
c    This routine calculates Bessel functions J sub(N+ALPHA) (X)
c    for non-negative argument X, and non-negative order N+ALPHA.
c
c    This program is based on a program written by David Sookne
c    that computes values of the Bessel functions J or I of real
c    argument and integer order.  Modifications include the restriction
c    of the computation to the J Bessel function of non-negative real
c    argument, the extension of the computation to arbitrary positive
c    order, and the elimination of most underflow.
c
c    In case of an error,  NCALC .NE. NB, and not all J's are
c    calculated to the desired accuracy.
c
c    NCALC .LT. 0:  An argument is out of range. For example,
c    NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
c    In this case, B(1) is set to zero, the remainder of the
c    B-vector is not calculated, and NCALC is set to
c    MIN(NB,0)-1 so that NCALC .NE. NB.
c
c    NB .GT. NCALC .GT. 0: Not all requested function values could
c    be calculated accurately.  This usually occurs because NB is
c    much larger than ABS(X).  In this case, B(N) is calculated
c    to the desired accuracy for N .LE. NCALC, but precision
c    is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
c    for N .GT. NCALC (because it is too small to be represented),
c    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
c    significant figures of B(N) can be trusted.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    Frank Olver, David Sookne,
c    A Note on Backward Recurrence Algorithms,
c    Mathematics of Computation,
c    Volume 26, 1972, pages 941-947.
c
c    David Sookne,
c    Bessel Functions of Real Argument and Integer Order,
c    NBS Journal of Res. B,
c    Volume 77B, 1973, pages 125-132.
c
c  Parameters:
c
c    Input, double precision X, the argument for which the
c    J's are to be calculated.
c
c    Input, double precision ALPHA, the fractional part of order for which
c    the J's or exponentially scaled J's (J*exp(X)) are to be calculated.  
c    0 <= ALPHA < 1.0.
c
c    Input, integer NB, the number of functions to be calculated.
c    0 < NB.  The first function calculated is of order ALPHA, and the
c    last is of order (NB - 1 + ALPHA).
c
c    Output, double precision B(NB).  If RJBESL terminates normally, with
c    NCALC = NB, then B contains the functions J/ALPHA/(X) through 
c    J/NB-1+ALPHA/(X), or the corresponding exponentially scaled functions.
c
c    Output, integer NCALC, error indicator.  If NCALC = NB, then all the 
c    requested values were calculated to the desired accuracy.  
c
c  Local Parameters:
c
c    IT, the number of bits in the mantissa of a working precision
c    variable.
c
c    NSIG, the decimal significance desired.  Should be set to
c    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
c    in decreased accuracy while setting NSIG higher will
c    increase CPU time without increasing accuracy.  The
c    truncation error is limited to a relative error of
c    T=.5*10**(-NSIG).
c
c    ENTEN = 10.0**K, where K is the largest integer such that
c    ENTEN is machine-representable in working precision
c
c    ENSIG = 10.0**NSIG
c
c    RTNSIG = 10.0**(-K) for the smallest integer K such that
c    K .GE. NSIG/4
c
c    ENMTEN, the smallest ABS(X) such that X/4 does not underflow
c
c    XLARGE, the upper limit on the magnitude of X.  If ABS(X)=N,
c    then at least N iterations of the backward recursion
c    will be executed.  The value of 10000.0 is used on
c    every machine.
c
      implicit none

      integer nb

      double precision alpha
      double precision alpem
      double precision alp2em
      double precision b(nb)
      double precision capp
      double precision capq
      double precision eighth
      double precision em
      double precision en
      double precision enmten
      double precision ensig
      double precision enten
      double precision fact(25)
      double precision four
      double precision gnu
      double precision half
      double precision halfx
      integer i
      integer j
      integer k
      integer l
      integer m
      integer magx
      integer n
      integer nbmx
      integer ncalc
      integer nend
      integer nstart
      double precision one
      double precision one30
      double precision p
      double precision pi2
      double precision plast
      double precision pold
      double precision psave
      double precision psavel
      double precision r8_gamma
      double precision rtnsig
      double precision s
      double precision sum
      double precision t
      double precision t1
      double precision tempa
      double precision tempb
      double precision tempc
      double precision test
      double precision three
      double precision three5
      double precision tover
      double precision two
      double precision twofiv
      double precision twopi1
      double precision twopi2
      double precision x
      double precision xc
      double precision xin
      double precision xk
      double precision xlarge
      double precision xm
      double precision vcos
      double precision vsin
      double precision z
      double precision zero
c
c  Mathematical constants
c
c   PI2    - 2 / PI
c   TWOPI1 - first few significant digits of 2 * PI
c   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
c            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
c
      data pi2 / 0.636619772367581343075535d0 /
      data twopi1 / 6.28125d0 /
      data twopi2 / 1.935307179586476925286767d-3 /
      data zero /0.0d0 /
      data eighth / 0.125d0 /
      data half / 0.5d0 /
      data one / 1.0d0/
      data two /2.0d0 / 
      data three / 3.0d0 /
      data four / 4.0d0 /
      data twofiv /25.0d0/
      data one30 /130.0d0 /
      data three5 / 35.0d0/
c
c  Machine-dependent parameters
c
      data enten /1.0d38 /
      data ensig / 1.0d17 /
      data rtnsig / 1.0d-4/
      data enmten /1.2d-37 /
      data xlarge / 1.0d4/
c
c  Factorial(N)
c
      data fact / 
     & 1.0d0,
     & 1.0d0,
     & 2.0d0,
     & 6.0d0,
     & 24.0d0,
     & 1.2d2,
     & 7.2d2,
     & 5.04d3,
     & 4.032d4,
     & 3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9,
     & 8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14,
     & 6.402373705728d15,1.21645100408832d17,2.43290200817664d18,
     & 5.109094217170944d19,1.12400072777760768d21,
     & 2.585201673888497664d22,
     & 6.2044840173323943936d23/
c
c  Check for out of range arguments.
c
      magx = int ( x )

      if ( 
     &  0 .lt. nb .and. 
     &  zero .le. x .and. 
     &  x .le. xlarge .and. 
     &  zero .le. alpha .and. 
     &  alpha .lt. one ) then
c
c  Initialize result array to zero.
c
        ncalc = nb
        do i = 1, nb
          b(i) = zero
        end do
c
c  Branch to use 2-term ascending series for small X and asymptotic
c  form for large X when NB is not too large.
c
        if ( x .lt. rtnsig ) then
c
c  Two-term ascending series for small X.
c
          tempa = one
          alpem = one + alpha
          halfx = zero

          if ( enmten .lt. x ) then
            halfx = half * x
          end if

          if ( alpha .ne. zero ) then
            tempa = halfx**alpha / ( alpha * r8_gamma ( alpha ) )
          end if

          tempb = zero

          if ( one .lt. x + one ) then
            tempb = -halfx * halfx
          end if

          b(1) = tempa + tempa * tempb / alpem

          if ( x .ne. zero .and. b(1) .eq. zero ) then
            ncalc = 0
          end if

          if ( nb .ne. 1 ) then

            if ( x .le. zero ) then

              do n = 2, nb
                b(n) = zero
              end do
c
c  Calculate higher order functions.
c
            else

              tempc = halfx
              tover = ( enmten + enmten ) / x

              if ( tempb .ne. zero ) then
                tover = enmten / tempb
              end if

              do n = 2, nb

                tempa = tempa / alpem
                alpem = alpem + one
                tempa = tempa * tempc

                if ( tempa .le. tover * alpem ) then
                  tempa = zero
                end if

                b(n) = tempa + tempa * tempb / alpem

                if ( b(n) .eq. zero .and. n .lt. ncalc ) then
                  ncalc = n - 1
                end if

              end do

            end if
          end if
c
c  Asymptotic series for 21 < X.
c
        else if ( twofiv .lt. x .and. nb .le. magx + 1 ) then

          xc = sqrt ( pi2 / x )
          xin = ( eighth / x )**2
          m = 11

          if ( x .ge. three5 ) then
            m = 8
          end if

          if ( x .ge. one30 ) then
            m = 4
          end if

          xm = four * dble ( m )
c
c  Argument reduction for SIN and COS routines.
c
          t = aint ( x / ( twopi1 + twopi2 ) + half )
          z = ( ( x - t * twopi1 ) - t * twopi2 ) 
     &      - ( alpha + half ) / pi2
          vsin = sin ( z )
          vcos = cos ( z )
          gnu = alpha + alpha

          do i = 1, 2

            s = ( ( xm - one ) - gnu ) * ( ( xm - one ) + gnu ) 
     &        * xin * half
            t = ( gnu - ( xm - three ) ) * ( gnu + ( xm - three ) )
            capp = s * t / fact(2*m+1)
            t1 = ( gnu - ( xm + one ) ) * ( gnu + ( xm + one ) )
            capq = s * t1 / fact(2*m+2)
            xk = xm
            k = m + m
            t1 = t

            do j = 2, m
              xk = xk - four
              s = ( ( xk - one ) - gnu ) * ( ( xk - one ) + gnu )
              t = ( gnu - ( xk - three ) ) * ( gnu + ( xk - three ) )
              capp = ( capp + one / fact(k-1) ) * s * t * xin
              capq = ( capq + one / fact(k) ) * s * t1 * xin
              k = k - 2
              t1 = t
            end do

            capp = capp + one
            capq = ( capq + one ) * ( gnu * gnu - one ) * ( eighth / x )
            b(i) = xc * ( capp * vcos - capq * vsin )

            if ( nb .eq. 1 ) then
              return
            end if

            t = vsin
            vsin = -vcos
            vcos = t
            gnu = gnu + two

          end do
c
c  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1.
c
          if ( 2 .lt. nb ) then
            gnu = alpha + alpha + two
            do j = 3, nb
              b(j) = gnu * b(j-1) / x - b(j-2)
              gnu = gnu + two
            end do
          end if
c
c  Use recurrence to generate results.  First initialize the
c  calculation of P's.
c
        else

          nbmx = nb - magx
          n = magx + 1
          en = dble ( n + n ) + ( alpha + alpha )
          plast = one
          p = en / x
c
c  Calculate general significance test.
c
          test = ensig + ensig
c
c  Calculate P's until N = NB-1.  Check for possible overflow.
c
          if ( 3 .le. nbmx ) then

            tover = enten / ensig
            nstart = magx + 2
            nend = nb - 1
            en = dble ( nstart + nstart ) - two + ( alpha + alpha )

            do k = nstart, nend

              n = k
              en = en + two
              pold = plast
              plast = p
              p = en * plast / x - pold
c
c  To avoid overflow, divide P's by TOVER.  Calculate P's until
c  1 < ABS(P).
c
              if ( tover .lt. p ) then

                tover = enten
                p = p / tover
                plast = plast / tover
                psave = p
                psavel = plast
                nstart = n + 1

  100           continue

                n = n + 1
                en = en + two
                pold = plast
                plast = p
                p = en * plast / x - pold

                if ( p .le. one ) then
                  go to 100
                end if

                tempb = en / x
c
c  Calculate backward test and find NCALC, the highest N such that
c  the test is passed.
c
                test = pold * plast 
     &            * ( half - half / ( tempb * tempb ) )
                test = test / ensig
                p = plast * tover
                n = n - 1
                en = en - two
                nend = min ( nb, n )

                do l = nstart, nend
                  pold = psavel
                  psavel = psave
                  psave = en * psavel / x - pold
                  if ( test .lt. psave * psavel ) then
                    ncalc = l - 1
                    go to 190
                  end if
                end do

                ncalc = nend
                go to 190

              end if

            end do

            n = nend
            en = dble ( n + n ) + ( alpha + alpha )
c
c  Calculate special significance test for 2 < NBMX.
c
            test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

          end if
c
c  Calculate P's until significance test passes.
c
  140     continue

          n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x - pold
          if ( p .lt. test ) then
            go to 140
          end if
c
c  Initialize the backward recursion and the normalization sum.
c
  190     continue

          n = n + 1
          en = en + two
          tempb = zero
          tempa = one / p
          m = 2 * n - 4 * ( n / 2 )
          sum = zero
          em = dble ( n / 2 )
          alpem = ( em - one ) + alpha
          alp2em = ( em + em ) + alpha

          if ( m .ne. 0 ) then
            sum = tempa * alpem * alp2em / em
          end if

          nend = n - nb
c
c  Recur backward via difference equation, calculating (but not
c  storing) B(N), until N = NB.
c
          if ( 0 .lt. nend ) then

            do l = 1, nend

              n = n - 1
              en = en - two
              tempc = tempb
              tempb = tempa
              tempa = ( en * tempb ) / x - tempc
              m = 2 - m

              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                if ( n .eq. 1 ) then
                  go to 210
                end if
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + tempa * alp2em ) * alpem / em
              end if

            end do

          end if
c
c  Store B(NB).
c
  210     continue

          b(n) = tempa

          if ( nend .ge. 0 ) then

            if ( nb .le. 1 ) then

              alp2em = alpha
              if ( alpha + one .eq. one ) then
                alp2em = one
              end if
              sum = sum + b(1) * alp2em
              go to 250

            else
c
c  Calculate and store B(NB-1).
c
              n = n - 1
              en = en - two
              b(n) = ( en * tempa ) / x - tempb

              if ( n .eq. 1 ) then
                go to 240
              end if

              m = 2 - m

              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + b(n) * alp2em ) * alpem / em
              end if

            end if

          end if

          nend = n - 2
c
c  Calculate via difference equation and store B(N), until N = 2.
c
          if ( nend .ne. 0 ) then

            do l = 1, nend
              n = n - 1
              en = en - two
              b(n) = ( en * b(n+1) ) / x - b(n+2)
              m = 2 - m
              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + b(n) * alp2em ) * alpem / em
              end if
            end do

          end if
c
c  Calculate B(1).
c
          b(1) = two * ( alpha + one ) * b(2) / x - b(3)

  240     continue

          em = em - one
          alp2em = ( em + em ) + alpha

          if ( alp2em .eq. zero ) then
            alp2em = one
          end if

          sum = sum + b(1) * alp2em
c
c  Normalize.  Divide all B(N) by sum.
c
  250     continue

          if ( alpha + one .ne. one ) then
            sum = sum * r8_gamma ( alpha ) * ( x * half )**( -alpha )
          end if

          tempa = enmten
          if ( one .lt. sum ) then
            tempa = tempa * sum
          end if

          do n = 1, nb
            if ( abs ( b(n) ) .lt. tempa ) then
              b(n) = zero
            end if
            b(n) = b(n) / sum
          end do

        end if
c
c  Error return: X, NB, or ALPHA is out of range.
c
      else
        b(1) = zero
        ncalc = min ( nb, 0 ) - 1
      end if

      return
      end
      subroutine rkbesl ( x, alpha, nb, ize, bk, ncalc )

c*********************************************************************72
c
cc RKBESL calculates K Bessel function with non-integer orders.
c
c  Discussion:
c
c    This routine calculates modified Bessel functions of the second 
c    kind, K SUB(N+ALPHA) (X), for non-negative argument X, and 
c    non-negative order N+ALPHA, with or without exponential scaling.
c
c    This program is based on a program written by J. B. Campbell
c    that computes values of the Bessel functions K of real
c    argument and real order.  Modifications include the addition
c    of non-scaled functions, parameterization of machine
c    dependencies, and the use of more accurate approximations
c    for SINH and SIN.
c
c    In case of an error, NCALC .NE. NB, and not all K's are
c    calculated to the desired accuracy.
c
c    NCALC .LT. -1:  An argument is out of range. For example,
c    NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
c    XMAX.  In this case, the B-vector is not calculated,
c    and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
c
c    NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
c    K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
c    the B-vector is not calculated.  Note that again
c    NCALC .NE. NB.
c
c    0 .LT. NCALC .LT. NB: Not all requested function values could
c    be calculated accurately.  BK(I) contains correct function
c    values for I .LE. NCALC, and contains the ratios
c    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Reference:
c
c    JB Campbell,
c    On Temme's Algorithm for the Modified Bessel Functions of the 
c    Third Kind,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 581-586.
c
c    JB Campbell,
c    A FORTRAN IV Subroutine for the Modified Bessel Functions of
c    the Third Kind of Real Order and Real Argument,
c    Report NRC/ERB-925,
c    National Research Council, Canada.
c
c  Parameters:
c
c    Input, double precision X, the non-negative argument for which
c    K's or exponentially scaled K's (K*EXP(X))
c    are to be calculated.  If K's are to be calculated,
c    X must not be greater than XMAX.
c
c    Input, double precision ALPHA, the fractional part of order for which
c    K's or exponentially scaled K's (K*EXP(X)) are to be calculated.  
c    0 <= ALPHA < 1.0.
c
c    Input, integer NB, the number of functions to be calculated, NB .GT. 0.
c    The first function calculated is of order ALPHA, and the
c    last is of order (NB - 1 + ALPHA).
c
c    Input, integer IZE, scaling option.
c    1, unscaled functions are to calculated,
c    2, exponentially scaled functions are to be calculated.
c
c    Output, double precision BK(NB), the results.  If the routine
c    terminates normally, with NCALC = NB, the vector BK contains the 
c    functions K(ALPHA,X), ... , K(NB-1+ALPHA,X), or the corresponding 
c    exponentially scaled functions.
c    If (0 .LT. NCALC .LT. NB), BK(I) contains correct function
c    values for I .LE. NCALC, and contains the ratios
c    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
c
c    Output, integer NCALC, error indicator.  If NCALC = NB, then all the 
c    requested values were calculated to the desired accuracy.  
c
      implicit none

      double precision a
      double precision alpha
      double precision blpha
      double precision bk(1)
      double precision bk1
      double precision bk2
      double precision c
      double precision d
      double precision dm
      double precision d1
      double precision d2
      double precision d3
      double precision enu
      double precision eps
      double precision estf(7)
      double precision estm(6)
      double precision ex
      double precision four
      double precision f0
      double precision f1
      double precision f2
      double precision half
      integer i
      integer iend
      integer itemp
      integer ize
      integer j
      integer k
      integer m
      integer mplus1
      integer nb
      integer ncalc
      double precision one
      double precision p(8)
      double precision p0
      double precision q(7)
      double precision q0
      double precision r(5)
      double precision ratio
      double precision s(4)
      double precision sqxmin
      double precision t(6)
      double precision tinyx
      double precision two
      double precision twonu
      double precision twox
      double precision t1
      double precision t2
      double precision wminf
      double precision x
      double precision xinf
      double precision xmax
      double precision xmin
      double precision x2by4
      double precision zero
c
c  Mathematical constants
c    A = LOG(2.D0) - Euler's constant
c    D = SQRT(2.D0/PI)
c
      data half / 0.5d0 /
      data one / 1.0d0 /
      data two / 2.0d0 /
      data zero /0.0d0/
      data four /4.0d0 /
      data tinyx / 1.0d-10/
      data a / 0.11593151565841244881d0/
      data d /0.797884560802865364d0/
c
c  Machine dependent parameters
c
      data eps /2.22d-16/
      data sqxmin /1.49d-154/
      data xinf /1.79d+308/
      data xmin /2.23d-308/
      data xmax /705.342d0/
c
c  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
c                                         + Euler's constant
c  Coefficients converted from hex to decimal and modified
c  by W. J. Cody, 2/26/82
c
c  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
c  T    - Approximation for SINH(Y)/Y
c
      data p/ 0.805629875690432845d00,    0.204045500205365151d02,
     &        0.157705605106676174d03,    0.536671116469207504d03,
     &        0.900382759291288778d03,    0.730923886650660393d03,
     &        0.229299301509425145d03,    0.822467033424113231d00/
      data q/ 0.294601986247850434d02,    0.277577868510221208d03,
     &        0.120670325591027438d04,    0.276291444159791519d04,
     &        0.344374050506564618d04,    0.221063190113378647d04,
     &        0.572267338359892221d03/
      data r/-0.48672575865218401848d+0,  0.13079485869097804016d+2,
     &       -0.10196490580880537526d+3,  0.34765409106507813131d+3,
     &        0.34958981245219347820d-3/
      data s/-0.25579105509976461286d+2,  0.21257260432226544008d+3,
     &       -0.61069018684944109624d+3,  0.42269668805777760407d+3/
      data t/ 0.16125990452916363814d-9, 0.25051878502858255354d-7,
     &        0.27557319615147964774d-5, 0.19841269840928373686d-3,
     &        0.83333333333334751799d-2, 0.16666666666666666446d+0/
      data estm/5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2,
     &          9.3715d0/
      data estf/4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0,
     &          8.45096d1, 2.0d1/

      ex = x
      enu = alpha
      ncalc = min ( nb, 0 ) - 2

      if ( 0 .lt. nb .and. 
     &  ( enu .ge. zero .and. enu .lt. one ) .and. 
     &  ( ize .ge. 1 .and. ize .le. 2 ) .and.
     &  ( ize .ne. 1 .or. ex .le. xmax ) .and.
     &  zero .lt. ex )  then

        k = 0
        if ( enu .lt. sqxmin ) then
          enu = zero
        end if

        if ( half .lt. enu ) then
          k = 1
          enu = enu - one
        end if

        twonu = enu + enu
        iend = nb + k - 1
        c = enu * enu
        d3 = -c

        if ( ex .le. one ) then
c
c  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA,
c                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA.
c
          d1 = zero
          d2 = p(1)
          t1 = one
          t2 = q(1)

          do i = 2, 7, 2
            d1 = c * d1 + p(i)
            d2 = c * d2 + p(i+1)
            t1 = c * t1 + q(i)
            t2 = c * t2 + q(i+1)
          end do

          d1 = enu * d1
          t1 = enu * t1
          f1 = log ( ex )
          f0 = a + enu * ( p(8) 
     &       - enu * ( d1 + d2 ) / ( t1 + t2 ) ) - f1
          q0 = exp ( -enu * ( a - enu * 
     &       ( p(8) + enu * ( d1 - d2 ) / ( t1 - t2 ) ) - f1 ) )
          f1 = enu * f0
          p0 = exp ( f1 )
c
c  Calculation of F0.
c
          d1 = r(5)
          t1 = one
          do i = 1, 4
            d1 = c * d1 + r(i)
            t1 = c * t1 + s(i)
          end do

          if ( abs ( f1 ) .le. half ) then
            f1 = f1 * f1
            d2 = zero
            do i = 1, 6
              d2 = f1 * d2 + t(i)
            end do
            d2 = f0 + f0 * f1 * d2
          else
            d2 = sinh ( f1 ) / enu
          end if

          f0 = d2 - enu * d1 / ( t1 * p0 )
c
c  X <= 1.0E-10.
c
c  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X).
c
          if ( ex .le. tinyx ) then

            bk(1) = f0 + ex * f0

            if ( ize .eq. 1 ) then
              bk(1) = bk(1) - ex * bk(1)
            end if

            ratio = p0 / f0
            c = ex * xinf
c
c  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
c  1/2 <= ALPHA.
c
            if ( k .ne. 0 ) then

              ncalc = -1

              if ( bk(1) .ge. c / ratio ) then
                return
              end if

              bk(1) = ratio * bk(1) / ex
              twonu = twonu + two
              ratio = twonu

            end if

            ncalc = 1

            if ( nb .eq. 1 ) then
              return
            end if
c
c  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1.
c
            ncalc = -1
            do i = 2, nb
              if ( ratio .ge. c ) then
                return
              end if
              bk(i) = ratio / ex
              twonu = twonu + two
              ratio = twonu
            end do

            ncalc = 1
            go to 420
c
c  1.0E-10 < X <= 1.0.
c
          else

            c = one
            x2by4 = ex * ex / four
            p0 = half * p0
            q0 = half * q0
            d1 = -one
            d2 = zero
            bk1 = zero
            bk2 = zero
            f1 = f0
            f2 = p0

  100       continue

            d1 = d1 + two
            d2 = d2 + one
            d3 = d1 + d3
            c = x2by4 * c / d2
            f0 = ( d2 * f0 + p0 + q0 ) / d3
            p0 = p0 / ( d2 - enu )
            q0 = q0 / ( d2 + enu )
            t1 = c * f0
            t2 = c * ( p0 - d2 * f0 )
            bk1 = bk1 + t1
            bk2 = bk2 + t2

            if ( eps .lt. abs ( t1 / ( f1 + bk1 ) ) .or.
     &           eps .lt. abs ( t2 / ( f2 + bk2 ) ) )  then
              go to 100
            end if

            bk1 = f1 + bk1
            bk2 = two * ( f2 + bk2 ) / ex

            if ( ize .eq. 2 ) then
              d1 = exp ( ex )
              bk1 = bk1 * d1
              bk2 = bk2 * d1
            end if

            wminf = estf(1) * ex + estf(2)

          end if
c
c  1/EPS < X.
c
        else if ( one .lt. eps * ex ) then

          ncalc = nb
          bk1 = one / ( d * sqrt ( ex ) )
          do i = 1, nb
            bk(i) = bk1
          end do

          return

        else
c
c  1 < X.
c
          twox = ex + ex
          blpha = zero
          ratio = zero

          if ( ex .le. four ) then
c
c  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0.
c
            d2 = aint ( estm(1) / ex + estm(2) )
            m = int ( d2 )
            d1 = d2 + d2
            d2 = d2 - half
            d2 = d2 * d2
            do i = 2, m
              d1 = d1 - two
              d2 = d2 - d1
              ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
            end do
c
c  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
c  recurrence and K(ALPHA,X) from the Wronskian.
c
            d2 = aint ( estm(3) * ex + estm(4) )
            m = int ( d2 )
            c = abs ( enu )
            d3 = c + c
            d1 = d3 - one
            f1 = xmin
            f0 = ( two * ( c + d2 ) / ex 
     &         + half * ex / ( c + d2 + one ) ) * xmin

            do i = 3, m
              d2 = d2 - one
              f2 = ( d3 + d2 + d2 ) * f0
              blpha = ( one + d1 / d2 ) * ( f2 + blpha )
              f2 = f2 / ex + f1
              f1 = f0
              f0 = f2
            end do

            f1 = ( d3 + two ) * f0 / ex + f1
            d1 = zero
            t1 = one
            do i = 1, 7
              d1 = c * d1 + p(i)
              t1 = c * t1 + q(i)
            end do

            p0 = exp ( c * ( a + c * ( p(8) 
     &         - c * d1 / t1 ) - log ( ex ) ) ) / ex
            f2 = ( c + half - ratio ) * f1 / ex
            bk1 = p0 + ( d3 * f0 - f2 + f0 + blpha ) 
     &        / ( f2 + f1 + f0 ) * p0

            if ( ize .eq. 1 ) then
              bk1 = bk1 * exp ( - ex )
            end if

            wminf = estf(3) * ex + estf(4)

          else
c
c  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
c  recurrence, for 4 < X.
c
            dm = aint ( estm(5) / ex + estm(6) )
            m = int ( dm )
            d2 = dm - half
            d2 = d2 * d2
            d1 = dm + dm

            do i = 2, m
              dm = dm - one
              d1 = d1 - two
              d2 = d2 - d1
              ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
              blpha = ( ratio + ratio * blpha ) / dm
            end do

            bk1 = one / ( ( d + d * blpha ) * sqrt ( ex ) )

            if ( ize .eq. 1 ) then
              bk1 = bk1 * exp ( - ex )
            end if

            wminf = estf(5) * ( ex - abs ( ex - estf(7) ) ) + estf(6)

          end if
c
c  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
c  K(ALPHA+1,X)/K(ALPHA,X).
c
          bk2 = bk1 + bk1 * ( enu + half - ratio ) / ex

        end if
c
c  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
c  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1.
c
        ncalc = nb
        bk(1) = bk1

        if ( iend .eq. 0 ) then
          return
        end if

        j = 2 - k 

        if ( 0 .lt. j ) then
          bk(j) = bk2
        end if

        if ( iend .eq. 1 ) then
          return
        end if

        m = min ( int ( wminf - enu ), iend )

        do i = 2, m

          t1 = bk1
          bk1 = bk2
          twonu = twonu + two

          if ( ex .lt. one ) then

            if ( ( xinf / twonu ) * ex .le. bk1 ) then
              go to 195
            end if
            go to 187

          else

            if ( bk1 / ex .ge. xinf / twonu ) then
              go to 195
            end if

          end if

  187     continue

          bk2 = twonu / ex * bk1 + t1
          itemp = i
          j = j + 1

          if ( 0 .lt. j ) then
            bk(j) = bk2
          end if

        end do

  195   continue

        m = itemp

        if ( m .eq. iend ) then
          return
        end if

        ratio = bk2 / bk1
        mplus1 = m + 1
        ncalc = -1

        do i = mplus1, iend

          twonu = twonu + two
          ratio = twonu / ex + one / ratio
          j = j + 1

          if ( 1 .lt. j ) then
            bk(j) = ratio
          else
            if ( bk2 .ge. xinf / ratio ) then
                       return
            end if
            bk2 = ratio * bk2
          end if

        end do

        ncalc = max ( mplus1 - k, 1 )

        if ( ncalc .eq. 1 ) then
          bk(1) = bk2
        end if

        if ( nb .eq. 1 ) then
          return
        end if

  420   continue

        j = ncalc + 1

        do i = j, nb
          if ( bk(ncalc) .ge. xinf / bk(i) ) then
            return
          end if
          bk(i) = bk(ncalc) * bk(i)
          ncalc = i
        end do

      end if

      return
      end
      subroutine rybesl ( x, alpha, nb, by, ncalc )

c*********************************************************************72
c
cc RYBESL calculates Y Bessel function with non-integer orders.
c
c  Discussion:
c
c    This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
c    for non-negative argument X, and non-negative order N+ALPHA.
c
c    This program draws heavily on Temme's Algol program for Y(a,x)
c    and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
c    scheme is used for x < THRESH, and Campbell's scheme is used
c    in the asymptotic region.  Segments of code from both sources
c    have been translated into Fortran77, merged, and heavily modified.
c    Modifications include parameterization of machine dependencies,
c    use of a new approximation for ln(gamma(x)), and built-in
c    protection against over/underflow.
c
c    In case of an error, NCALC .NE. NB, and not all Y's are
c    calculated to the desired accuracy.
c
c    NCALC .LT. -1:  An argument is out of range. For example,
c    NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
c    XMAX.  In this case, BY(1) = 0.0, the remainder of the
c    BY-vector is not calculated, and NCALC is set to
c    MIN0(NB,0)-2  so that NCALC .NE. NB.
c
c    NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
c    values are set to 0.0.
c
c    1 .LT. NCALC .LT. NB: Not all requested function values could
c    be calculated accurately.  BY(I) contains correct function
c    values for I .LE. NCALC, and and the remaining NB-NCALC
c    array elements contain 0.0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    JB Campbell,
c    Bessel functions J_nu(x) and Y_nu(x) of real order and real argument,
c    Computational Physics Communications,
c    Volume 18, 1979, pages 133-142.
c
c    NM Temme,
c    On the numerical evaluation of the ordinary Bessel function 
c    of the second kind,
c    Journal of Computational Physics,
c    Volume 21, 1976, pages 343-350.
c
c  Parameters:
c
c    Input, double precision X, the argument.  0 <= X.
c
c    Input, double precision ALPHA, the fractional part of the order 
c    for which the Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
c
c    Input, integer NB, the number of functions to be calculated, NB .GT. 0.
c    The first function calculated is of order ALPHA, and the
c    last is of order (NB - 1 + ALPHA).
c
c    Output, double precision BY(NB).  If the routine terminates normally
c    (NCALC=NB), the vector BY contains the functions Y(ALPHA,X) through 
c    Y(NB-1+ALPHA,X),  If (0 .LT. NCALC .LT. NB), BY(I) contains correct 
c    function values for I .LE. NCALC, and contains the ratios
c    Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
c
c    Output, integer NCALC, error flag.  Before using the vector BY, the 
c    user should check that NCALC=NB, i.e., all orders have been calculated
c    to the desired accuracy.
c
      implicit none

      integer nb

      double precision alfa
      double precision alpha
      double precision aye
      double precision b
      double precision by(nb)
      double precision c
      double precision ch(21)
      double precision cosmu
      double precision d
      double precision del
      double precision den
      double precision ddiv
      double precision div
      double precision dmu
      double precision d1
      double precision d2
      double precision e
      double precision eight
      double precision en
      double precision enu
      double precision en1
      double precision eps
      double precision even
      double precision ex
      double precision f
      double precision fivpi
      double precision g
      double precision gamma
      double precision h
      double precision half
      integer i
      integer k
      integer na
      integer ncalc
      double precision odd
      double precision onbpi
      double precision one
      double precision one5
      double precision p
      double precision pa
      double precision pa1
      double precision pi
      double precision piby2
      double precision pim5
      double precision q
      double precision qa
      double precision qa1
      double precision q0
      double precision r
      double precision s
      double precision sinmu
      double precision sq2bpi
      double precision ten9
      double precision term
      double precision three
      double precision thresh
      double precision two
      double precision twobyx
      double precision x
      double precision xinf
      double precision xlarge
      double precision xmin
      double precision xna
      double precision x2
      double precision ya
      double precision ya1
      double precision zero
c
c  Mathematical constants
c    FIVPI = 5*PI
c    PIM5 = 5*PI - 15
c    ONBPI = 1/PI
c    PIBY2 = PI/2
c    SQ2BPI = SQUARE ROOT OF 2/PI
c
      data zero / 0.0d0 /
      data half / 0.5d0 /
      data one / 1.0d0 /
      data two / 2.0d0 /
      data three / 3.0d0 /
      data eight /8.0d0 /
      data one5 / 15.0d0 /
      data ten9 / 1.9d1/
      data fivpi /1.5707963267948966192d1 /
      data piby2 / 1.5707963267948966192d0/
      data pi /3.1415926535897932385d0 /
      data sq2bpi / 7.9788456080286535588d-1/
      data pim5 /7.0796326794896619231d-1/
      data onbpi / 3.1830988618379067154d-1/
c
c  Machine-dependent constants
c
      data del / 1.0d-8 /
      data xmin / 4.46d-308 /
      data xinf / 1.79d308 /
      data eps / 1.11d-16 /
      data thresh / 16.0d0 /
      data xlarge / 1.0d8 /
c
c  Coefficients for Chebyshev polynomial expansion of
c  1/gamma(1-x), abs(x) <= .5
c
      data ch/-0.67735241822398840964d-23,-0.61455180116049879894d-22,
     &         0.29017595056104745456d-20, 0.13639417919073099464d-18,
     &         0.23826220476859635824d-17,-0.90642907957550702534d-17,
     &        -0.14943667065169001769d-14,-0.33919078305362211264d-13,
     &        -0.17023776642512729175d-12, 0.91609750938768647911d-11,
     &         0.24230957900482704055d-09, 0.17451364971382984243d-08,
     &        -0.33126119768180852711d-07,-0.86592079961391259661d-06,
     &        -0.49717367041957398581d-05, 0.76309597585908126618d-04,
     &         0.12719271366545622927d-02, 0.17063050710955562222d-02,
     &        -0.76852840844786673690d-01,-0.28387654227602353814d+00,
     &         0.92187029365045265648d+00/

      ex = x
      enu = alpha

      if ( 0 .lt. nb .and. 
     &  x .ge. xmin .and. 
     &  ex .lt. xlarge .and.
     &  enu .ge. zero .and. 
     &  enu .lt. one )  then

        xna = aint ( enu + half )
        na = int ( xna )

        if ( na .eq. 1 ) then
          enu = enu - xna
        end if

        if ( enu .eq. -half ) then

          p = sq2bpi / sqrt ( ex )
          ya = p * sin ( ex )
          ya1 = -p * cos ( ex )
c
c  Use Temme's scheme for small X.
c
        else if ( ex .lt. three ) then

          b = ex * half
          d = - log ( b )
          f = enu * d
          e = b**( -enu )

          if ( abs ( enu ) .lt. del ) then
            c = onbpi
          else
            c = enu / sin ( enu * pi )
          end if
c
c  Computation of sinh(f)/f.
c
          if ( abs ( f ) .lt. one ) then
            x2 = f * f
            en = ten9
            s = one
            do i = 1, 9
              s = s * x2 / en / ( en - one ) + one
              en = en - two
            end do
          else
            s = ( e - one / e ) * half / f
          end if
c
c  Computation of 1/gamma(1-a) using Chebyshev polynomials.
c
          x2 = enu * enu * eight
          aye = ch(1)
          even = zero
          alfa = ch(2)
          odd = zero

          do i = 3, 19, 2
            even = - ( aye + aye + even )
            aye = - even * x2 - aye + ch(i)
            odd = - ( alfa + alfa + odd )
            alfa = - odd * x2 - alfa + ch(i+1)
          end do

          even = ( even * half + aye ) * x2 - aye + ch(21)
          odd = ( odd + alfa ) * two
          gamma = odd * enu + even
c
c  End of computation of 1/gamma(1-a).
c
          g = e * gamma
          e = ( e + one / e ) * half
          f = two * c * ( odd * e + even * s * d )
          e = enu * enu
          p = g * c
          q = onbpi / g
          c = enu * piby2

          if ( abs ( c ) .lt. del ) then
            r = one
          else
            r = sin ( c ) / c
          end if

          r = pi * c * r * r
          c = one
          d = - b * b
          h = zero
          ya = f + r * q
          ya1 = p
          en = zero

  100     continue

          en = en + one

          if ( eps .lt. abs ( g / ( one + abs ( ya ) ) )
     &      + abs ( h / ( one + abs ( ya1 ) ) ) ) then
            f = ( f * en + p + q ) / ( en * en - e )
            c = c * d / en
            p = p / ( en - enu )
            q = q / ( en + enu )
            g = c * ( f + r * q )
            h = c * p - en * g
            ya = ya + g
            ya1 = ya1 + h
            go to 100
          end if

          ya = -ya
          ya1 = -ya1 / b

        else if ( ex .lt. thresh ) then
c
c  Use Temme's scheme for moderate X.
c
          c = ( half - enu ) * ( half + enu )
          b = ex + ex
          e = ( ex * onbpi * cos ( enu * pi ) / eps )
          e = e * e
          p = one
          q = -ex
          r = one + ex * ex
          s = r
          en = two

  200     continue

          if ( r * en * en .lt. e ) then
            en1 = en + one
            d = ( en - one + c / en ) / s
            p = ( en + en - p * d ) / en1
            q = ( -b + q * d ) / en1
            s = p * p + q * q
            r = r * s
            en = en1
            go to 200
          end if

          f = p / s
          p = f
          g = -q / s
          q = g

  220     continue

          en = en - one

          if ( zero .lt. en ) then
            r = en1 * ( two - p ) - two
            s = b + en1 * q
            d = ( en - one + c / en ) / ( r * r + s * s )
            p = d * r
            q = d * s
            e = f + one
            f = p * e - g * q
            g = q * e + p * g
            en1 = en
            go to 220
          end if

          f = one + f
          d = f * f + g * g
          pa = f / d
          qa = -g / d
          d = enu + half -p
          q = q + ex
          pa1 = ( pa * q - qa * d ) / ex
          qa1 = ( qa * q + pa * d ) / ex
          b = ex - piby2 * ( enu + half )
          c = cos ( b )
          s = sin ( b )
          d = sq2bpi / sqrt ( ex )
          ya = d * ( pa * s + qa * c )
          ya1 = d * ( qa1 * s - pa1 * c )
        else
c
c  Use Campbell's asymptotic scheme.
c
          na = 0
          d1 = aint ( ex / fivpi )
          i = int ( d1 )
          dmu = (( ex - one5 * d1 ) - d1 * pim5 ) 
     &      - ( alpha + half ) * piby2

          if ( i - 2 * ( i / 2 ) .eq. 0 ) then
            cosmu = cos ( dmu )
            sinmu = sin ( dmu )
          else
            cosmu = -cos ( dmu )
            sinmu = -sin ( dmu )
          end if

          ddiv = eight * ex
          dmu = alpha
          den = sqrt ( ex )

          do k = 1, 2

            p = cosmu
            cosmu = sinmu
            sinmu = -p
            d1 = ( two * dmu - one ) * ( two * dmu + one )
            d2 = zero
            div = ddiv
            p = zero
            q = zero
            q0 = d1 / div
            term = q0

            do i = 2, 20

              d2 = d2 + eight
              d1 = d1 - d2
              div = div + ddiv
              term = -term * d1 / div
              p = p + term
              d2 = d2 + eight
              d1 = d1 - d2
              div = div + ddiv
              term = term * d1 / div
              q = q + term

              if ( abs ( term ) .le. eps ) then
                go to 320
              end if

            end do

  320       continue

            p = p + one
            q = q + q0

            if ( k .eq. 1 ) then
              ya = sq2bpi * ( p * cosmu - q * sinmu ) / den
            else
              ya1 = sq2bpi * ( p * cosmu - q * sinmu ) / den
            end if

            dmu = dmu + one

          end do

        end if

        if ( na .eq. 1 ) then
          h = two * ( enu + one ) / ex
          if ( one .lt. h ) then
            if ( xinf / h .lt. abs ( ya1 ) ) then
              h = zero
              ya = zero
            end if
          end if
          h = h * ya1 - ya
          ya = ya1
          ya1 = h
        end if
c
c  Now have first one or two Y's.
c
        by(1) = ya
        by(2) = ya1

        if ( ya1 .eq. zero ) then

          ncalc = 1

        else

          aye = one + alpha
          twobyx = two / ex
          ncalc = 2

          do i = 3, nb

            if ( twobyx .lt. one ) then
              if ( abs ( by(i-1) ) * twobyx .ge. xinf / aye ) then
                go to 450
              end if
            else
              if ( abs ( by(i-1) ) .ge. xinf / aye / twobyx ) then
                go to 450
              end if
            end if

            by(i) = twobyx * aye * by(i-1) - by(i-2)
            aye = aye + one
            ncalc = ncalc + 1

          end do

        end if

  450   continue

        do i = ncalc + 1, nb
          by(i) = zero
        end do

      else

        by(1) = zero
        ncalc = min ( nb, 0 ) - 1

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
