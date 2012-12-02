      subroutine bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

c*********************************************************************72
c
cc BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
c
c  Discussion:
c
c    FXY is the probability that two variables A and B, which are
c    related by a bivariate normal distribution with correlation R,
c    respectively satisfy A .lt.= X and B .lt.= Y.
c
c    Mathematica can evaluate the bivariate normal CDF via the commands:
c
c      <<MultivariateStatisticsc
c      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    National Bureau of Standards,
c    Tables of the Bivariate Normal Distribution and Related Functions,
c    NBS, Applied Mathematics Series, Number 50, 1959.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, Y, the parameters of the function.
c
c    Output, double precision R, the correlation value.
c
c    Output, double precision FXY, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 41 )

      double precision fxy
      double precision fxy_vec(n_max)
      integer n_data
      double precision r
      double precision r_vec(n_max)
      double precision x
      double precision x_vec(n_max)
      double precision y
      double precision y_vec(n_max)

      save fxy_vec
      save r_vec
      save x_vec
      save y_vec

      data fxy_vec /
     &  0.02260327218569867D+00,
     &  0.1548729518584100D+00,
     &  0.4687428083352184D+00,
     &  0.7452035868929476D+00,
     &  0.8318608306874188D+00,
     &  0.8410314261134202D+00,
     &  0.1377019384919464D+00,
     &  0.1621749501739030D+00,
     &  0.1827411243233119D+00,
     &  0.2010067421506235D+00,
     &  0.2177751155265290D+00,
     &  0.2335088436446962D+00,
     &  0.2485057781834286D+00,
     &  0.2629747825154868D+00,
     &  0.2770729823404738D+00,
     &  0.2909261168683812D+00,
     &  0.3046406378726738D+00,
     &  0.3183113449213638D+00,
     &  0.3320262544108028D+00,
     &  0.3458686754647614D+00,
     &  0.3599150462310668D+00,
     &  0.3742210899871168D+00,
     &  0.3887706405282320D+00,
     &  0.4032765198361344D+00,
     &  0.4162100291953678D+00,
     &  0.6508271498838664D+00,
     &  0.8318608306874188D+00,
     &  0.0000000000000000D+00,
     &  0.1666666666539970D+00,
     &  0.2500000000000000D+00,
     &  0.3333333333328906D+00,
     &  0.5000000000000000D+00,
     &  0.7452035868929476D+00,
     &  0.1548729518584100D+00,
     &  0.1548729518584100D+00,
     &  0.06251409470431653D+00,
     &  0.7452035868929476D+00,
     &  0.1548729518584100D+00,
     &  0.1548729518584100D+00,
     &  0.06251409470431653D+00,
     &  0.6337020457912916D+00 /
      data r_vec /
     &   0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,
     &   0.500D+00, -0.900D+00, -0.800D+00, -0.700D+00, -0.600D+00,
     &  -0.500D+00, -0.400D+00, -0.300D+00, -0.200D+00, -0.100D+00,
     &   0.000D+00,  0.100D+00,  0.200D+00,  0.300D+00,  0.400D+00,
     &   0.500D+00,  0.600D+00,  0.700D+00,  0.800D+00,  0.900D+00,
     &   0.673D+00,  0.500D+00, -1.000D+00, -0.500D+00,  0.000D+00,
     &   0.500D+00,  1.000D+00,  0.500D+00,  0.500D+00,  0.500D+00,
     &   0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,
     &   0.500D+00 /
      data x_vec /
     &  -2.0D+00, -1.0D+00,  0.0D+00,  1.0D+00,  2.0D+00,
     &   3.0D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00,
     &  -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00,
     &  -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00,
     &  -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00,
     &   1.0D+00,  2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00,
     &   0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, -1.0D+00,
     &  -1.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, -1.0D+00,
     &   0.7071067811865475D+00 /
      data y_vec /
     &   1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,
     &   1.0D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,
     &   0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,
     &   0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,
     &   0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,
     &   0.5D+00,  1.0D+00,  0.0D+00,  0.0D+00,  0.0D+00,
     &   0.0D+00,  0.0D+00,  1.0D+00, -1.0D+00,  1.0D+00,
     &  -1.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, -1.0D+00,
     &   0.7071067811865475D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        r = 0.0D+00
        x = 0.0D+00
        y = 0.0D+00
        fxy = 0.0D+00
      else
        r = r_vec(n_data)
        x = x_vec(n_data)
        y = y_vec(n_data)
        fxy = fxy_vec(n_data)
      end if

      return
      end
      DOUBLE PRECISION FUNCTION BIVNOR ( AH, AK, R )

c*********************************************************************72
c
cc BIVNOR computes the bivariatenormal probability.
c
c  Discussion:
c
c    BIVNOR computes the probability for two normal variates X and Y
c    whose correlation is R, that X .gt. AH and Y .gt. AK.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Donnelly.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Thomas Donnelly,
c    Algorithm 462: Bivariate Normal Distribution,
c    Communications of the ACM,
c    October 1973, Volume 16, Number 10, page 638.
c
c  Parameters:
c
c    Input, double precision AH, AK, the limits of integration.
c
c    Input, double precision R, the correlation.
c
c    Output, double precision BIVNOR, the bivariate normal CDF.
c
      implicit none

      integer i
      integer idig
      integer is
      DOUBLE PRECISION TWOPI, B, AH, AK, R, GH, GK, RR, GAUSS,
     &  DERF, H2, A2, H4, DEXP, EX, W2, AP, S2, SP, S1, SN, SQR,
     &  DSQRT, CON, DATAN, WH, WK, GW, SGN, T, DABS, G2, CONEX,
     &  CN

      GAUSS ( T ) = ( 1.0D0 + DERF ( T / DSQRT ( 2.0D0 ) ) ) / 2.0D0
C
C  GAUSS IS A UNIVARIATE LOWER NORMAL
C  TAIL AREA CALCULATED HERE FROM THE
C  CENTRAL ERROR FUNCTION DERF.
C  IT MAY BE REPLACED BY THE ALGORITHM IN
C  HILL,I.D. AND JOYCE,S.A. ALGORITHM 304,
C  NORMAL CURVE INTEGRAL(515), COMM.A.C.M.(10)
C  (JUNE,1967),P.374.
C  SOURCE: OWERN, D.B. ANN.MATH.STAT.
C  VOL. 27(1956), P.1075.
C  TWOPI = 2. * PI
      TWOPI = 6.283185307179587D0
      B = 0.0D0
      IDIG = 15
C  THE PARAMETER 'IDIG' GIVES THE
C  NUMBER OF SIGNIFICANT DIGITS
C  TO THE RIGHT OF THE DECIMAL POINT
C  DESIRED IN THE ANSWER, IF
C  IT IS WITHIN THE COMPUTER'S
C  CAPACITY OF COURSE.
      GH = GAUSS ( -AH ) / 2.0D0
      GK = GAUSS ( -AK ) / 2.0D0
      IF ( R ) 10, 30, 10
10    RR = 1.0D0 - R * R
      IF ( RR ) 20, 40, 100
20    WRITE ( 6, 99999 ) R
C  ERROR EXIT FOR ABS(R) .GT. 1.0D=
99999 FORMAT ( 12H BIVNOR R IS, D26.16 )
      STOP
30    B = 4.0D0 * GH * GK
      GO TO 350
40    IF ( R ) 50, 70, 70
50    IF ( AH + AK ) 60, 350, 350
60    B = 2.0D0 * ( GH + GK ) - 1.0D0
      GO TO 350
70    IF ( AH - AK ) 80, 90, 90
80    B = 2.0D0 * GK
      GO TO 350
90    B = 2.0D0 * GH
      GO TO 350
100   SQR = DSQRT ( RR )
      IF ( IDIG - 15 ) 120, 110, 120
110   CON = TWOPI * 1.D-15 / 2.0D0
      GO TO 140
120   CON = TWOPI / 2.0D0
      DO 130 I = 1, IDIG
        CON = CON / 10.0D0
130   CONTINUE
140   IF ( AH ) 170, 150, 170
150   IF ( AK ) 190, 160, 190
c
c  AH = AK = 0.0
c  The original code had a mistake here, JVB, 23 November 2010.
c
160   continue
      b = 0.25D+00 + dasin ( r ) / twopi
      GO TO 350
170   B = GH
      IF ( AH * AK ) 180, 200, 190
180   B = B - 0.5D0
190   B = B + GK
      IF ( AH ) 200, 340, 200
200   WH = -AH
      WK = ( AK / AH - R ) / SQR
      GW = 2.0D0 * GH
      IS = -1
210   SGN = -1.0D0
      T = 0.0D0
      IF ( WK ) 220, 320, 220
220   IF ( DABS ( WK ) - 1.0D0 ) 270, 230, 240
230   T = WK * GW * ( 1.0D0 - GW ) / 2.0D0
      GO TO 310
240   SGN = -SGN
      WH = WH * WK
      G2 = GAUSS ( WH )
      WK = 1.0D0 / WK
      IF ( WK ) 250, 260, 260
250   B = B + 0.5D0
260   B = B - ( GW + G2 ) / 2.0D0 + GW * G2
270   H2 = WH * WH
      A2 = WK * WK
      H4 = H2 / 2.0D0
      EX = DEXP ( - H4 )
      W2 = H4 * EX
      AP = 1.0D0
      S2 = AP - EX
      SP = AP
      S1 = 0.0D0
      SN = S1
      CONEX = DABS ( CON / WK )
      GO TO 290
280   SN = SP
      SP = SP + 1.0D0
      S2 = S2 - W2
      W2 = W2 * H4 / SP
      AP = - AP * A2
290   CN = AP * S2 / ( SN + SP )
      S1 = S1 + CN
      IF ( DABS ( CN ) - CONEX ) 300, 300, 280
300   T = ( DATAN ( WK ) - WK * S1 ) / TWOPI
310   B = B + SGN * T
320   IF ( IS ) 330, 350, 350
330   IF ( AK ) 340, 350, 340
340   WH = -AK
      WK = ( AH / AK - R ) / SQR
      GW = 2.0D0 * GK
      IS = 1
      GO TO 210
350   IF ( B ) 360, 370, 370
360   B = 0.0D0
370   IF ( B - 1.0D0 ) 390, 390, 380
380   B = 1.0D0
390   BIVNOR = B
      RETURN
      END
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
      function q ( h, ah ) 

c*********************************************************************72
c 
cc Q computes (1/2) * p(H<Z) - T(H,AH).
c
c  Discussion:
c
c    The routine computes Q = (1/2) * P( H < Z ) - T ( H, AH ).
c
c    The result for Q is non-negative. 
c
c    Warning : Q is computed as the difference between two terms; 
c    When the two terms are of similar value this may produce 
c    error in Q. 
c
c  Modified:
c
c    16 January 2008
c
c  Author:
c
c    Original version by Mike Patefield, David Tandy.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Mike Patefield, David Tandy,
c    Fast and Accurate Calculation of Owen's T Function,
c    Journal of Statistical Software,
c    Volume 5, Number 5, 2000, pages 1-25.
c
c  Parameters:
c
c    Input, double precision H, the lower limit for Z.
c    0 < H.
c
c    Input, double precision AH, one of the arguments for the T function.
c
c    Output, double precision Q, the desired quantity.
c
      implicit none

      double precision ah
      double precision ahh
      double precision h
      integer ifail
      double precision q
      double precision rroot2
      parameter ( rroot2 = 0.70710678118654752440D+00 )
      double precision t
      double precision tfun
      double precision x
      double precision znorm1
      double precision znorm2

      if ( 1.0D+00 .lt. ah ) then 
        ahh = ah * h 
        q = tfun ( ahh, 1.0D+00 / ah, h) - znorm2 ( ahh ) * znorm1 ( h ) 
      else 
        q = 0.5D+00 * znorm2 ( h ) - t ( h, ah ) 
      end if 

      return 
      end
      function t ( h, a ) 

c*********************************************************************72
c
cc T computes Owen's T function for arbitrary H and A.
c
c  Modified:
c
c    13 April 2012
c
c  Author:
c
c    Original version by Mike Patefield, David Tandy.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Mike Patefield, David Tandy,
c    Fast and Accurate Calculation of Owen's T Function,
c    Journal of Statistical Software,
c    Volume 5, Number 5, 2000, pages 1-25.
c
c  Parameters:
c
c    Input, double precision H, A, the arguments.
c
c    Output, double precision T, the value of Owen's T function.
c
      implicit none

      double precision a
      double precision absa
      double precision absh
      double precision ah
      double precision cut
      parameter ( cut = 0.67D+00 )
      double precision h
      integer ifail
      double precision normah
      double precision normh
      double precision rroot2
      parameter ( rroot2 = 0.70710678118654752440D+00 )
      double precision t
      double precision tfun
      double precision x
      double precision znorm1
      double precision znorm2

      absh = dabs ( h )
      absa = dabs ( a ) 
      ah = absa * absh

      if ( absa .le. 1.0D+00 ) then

        t = tfun ( absh, absa, ah )
c
c  In the printed paper, the formula for T that follows was incorrect.
c
      else if ( absh .le. cut ) then

        t = 0.25D+00 - znorm1 ( absh ) * znorm1 ( ah )
     &    - tfun ( ah, 1.0D+00 / absa, absh ) 

      else

        normh = znorm2 ( absh ) 
        normah = znorm2 ( ah )
        t = 0.5D+00 * ( normh + normah ) - normh * normah 
     &  - tfun ( ah, 1.0D+00 / absa, absh ) 

      end if
 
      if ( a .lt. 0.0D+00 ) then
        t = - t 
      end if

      return 
      end 
      function tfun ( h, a, ah ) 

c*********************************************************************72
c
cc TFUN computes Owen's T function for a restricted range of parameters.
c
c  Discussion:
c
c    This routine computes Owen's T-function of H and A.
c
c    Originally named "TF", this routine was renamed "TFUN" to avoid a
c    conflict with a MATLAB built in function.
c
c    Thanks to Benjamin Sobotta for pointing out a missing factor of
c    0.5 that occurred where ZNORM1 was used, 15 December 2011.
c
c  Modified:
c
c    15 December 2011
c
c  Author:
c
c    Original version by Mike Patefield, David Tandy.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Mike Patefield, David Tandy,
c    Fast and Accurate Calculation of Owen's T Function,
c    Journal of Statistical Software,
c    Volume 5, Number 5, 2000, pages 1-25.
c
c  Parameters:
c
c    Input, double precision H, the H argument of the function.
c    0 <= H.
c
c    Input, double precision A, the A argument of the function.
c    0 <= A <= 1.
c
c    Input, double precision AH, the value of A*H.
c
c    Output, double precision TF, the value of Owen's T function.
c
      implicit none

      double precision a
      double precision ah
      double precision ai
      double precision aj
      double precision arange(7)
      double precision as
      double precision c2(21)
      double precision dhs
      double precision dj
      double precision gj
      double precision h
      double precision hrange(14)
      double precision hs
      integer i
      integer iaint
      integer icode
      integer ifail
      integer ihint
      integer ii
      integer j
      integer jj
      integer m
      integer maxii
      integer meth(18)
      double precision normh
      integer ord(18)
      double precision pts(13)
      double precision r
      double precision rroot2
      parameter ( rroot2 = 0.70710678118654752440D+00 )
      double precision rrtpi
      parameter ( rrtpi = 0.39894228040143267794D+00 )
      double precision rtwopi
      parameter ( rtwopi = 0.15915494309189533577D+00 )
      integer select(15,8)
      double precision tfun
      double precision vi
      double precision wts(13)
      double precision x
      double precision y
      double precision yi
      double precision z
      double precision zi
      double precision znorm1
      double precision znorm2

      data arange /
     &  0.025D+00, 0.09D+00, 0.15D+00, 0.36D+00, 0.5D+00, 
     &  0.9D+00, 0.99999D+00 / 

      data c2 /                         0.99999999999999987510D+00, 
     &    -0.99999999999988796462D+00,  0.99999999998290743652D+00, 
     &    -0.99999999896282500134D+00,  0.99999996660459362918D+00, 
     &    -0.99999933986272476760D+00,  0.99999125611136965852D+00, 
     &    -0.99991777624463387686D+00,  0.99942835555870132569D+00, 
     &    -0.99697311720723000295D+00,  0.98751448037275303682D+00, 
     &    -0.95915857980572882813D+00,  0.89246305511006708555D+00,
     &    -0.76893425990463999675D+00,  0.58893528468484693250D+00, 
     &    -0.38380345160440256652D+00,  0.20317601701045299653D+00, 
     &    -0.82813631607004984866D-01,  0.24167984735759576523D-01, 
     &    -0.44676566663971825242D-02,  0.39141169402373836468D-03 / 

      data hrange /
     &  0.02D+00, 0.06D+00, 0.09D+00, 0.125D+00, 0.26D+00, 
     &  0.4D+00,  0.6D+00,  1.6D+00,  1.7D+00,   2.33D+00,  
     &  2.4D+00,  3.36D+00, 3.4D+00,  4.8D+00 / 

      data meth / 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6 / 

      data ord /  2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0 / 

      data pts /                        0.35082039676451715489D-02, 
     &     0.31279042338030753740D-01,  0.85266826283219451090D-01, 
     &     0.16245071730812277011D+00,  0.25851196049125434828D+00, 
     &     0.36807553840697533536D+00,  0.48501092905604697475D+00, 
     &     0.60277514152618576821D+00,  0.71477884217753226516D+00, 
     &     0.81475510988760098605D+00,  0.89711029755948965867D+00, 
     &     0.95723808085944261843D+00,  0.99178832974629703586D+00 / 

      data select / 1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9, 
     &              1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9, 
     &              2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10, 
     &              2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10, 
     &              2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11, 
     &              2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12, 
     &              2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12, 
     &              2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12 / 

      data wts /                        0.18831438115323502887D-01, 
     &     0.18567086243977649478D-01,  0.18042093461223385584D-01, 
     &     0.17263829606398753364D-01,  0.16243219975989856730D-01, 
     &     0.14994592034116704829D-01,  0.13535474469662088392D-01, 
     &     0.11886351605820165233D-01,  0.10070377242777431897D-01, 
     &     0.81130545742299586629D-02,  0.60419009528470238773D-02, 
     &     0.38862217010742057883D-02,  0.16793031084546090448D-02 / 
c 
c  Determine appropriate method from t1...t6 
c
      ihint = 15

      do i = 1, 14 
        if ( h .le. hrange(i) ) then
          ihint = i
          go to 20 
        end if
      end do

   20 continue

      iaint = 8

      do i = 1, 7 
        if ( a .le. arange(i) ) then
          iaint = i
          go to 40 
        end if
      end do 

   40 continue

      icode = select(ihint,iaint) 
      m = ord(icode) 
c 
c  t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18 
c  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j! 
c  aj = a**(2j-1) / (2*pi) 
c 
      if ( meth(icode) .eq. 1 ) then

        hs = - 0.5D+00 * h * h 
        dhs = exp ( hs ) 
        as = a * a 
        j = 1 
        jj = 1 
        aj = rtwopi * a
        tfun = rtwopi * atan ( a ) 
        dj = dhs - 1.0D+00 
        gj = hs * dhs
 
  110   continue

        tfun = tfun + dj * aj / dble ( jj ) 

        if ( m .le. j ) then
          return
        end if 

        j = j + 1 
        jj = jj + 2 
        aj = aj * as 
        dj = gj - dj 
        gj = gj * hs / dble ( j ) 
        go to 110 
c 
c  t2(h, a, m) ; m = 10, 20 or 30 
c  z = (-1)**(i-1) * zi 
c  ii = 2i - 1 
c  vi = (-1)**(i-1) * a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi) 
c 
      else if ( meth(icode) .eq. 2 ) then

        maxii = m + m + 1 
        ii = 1 
        tfun = 0.0D+00 
        hs = h * h 
        as = - a * a 
        vi = rrtpi * a * exp ( - 0.5D+00 * ah * ah ) 
        z = 0.5D+00 * znorm1 ( ah ) / h 
        y = 1.0D+00 / hs 

  210   continue

        tfun = tfun + z 

        if ( maxii .le. ii ) then
          tfun = tfun * rrtpi * exp ( - 0.5D+00 * hs ) 
          return 
        end if

        z = y * ( vi - dble ( ii ) * z ) 
        vi = as * vi 
        ii = ii + 2 
        go to 210 
c 
c  t3(h, a, m) ; m = 20 
c  ii = 2i - 1 
c  vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi) 
c 
      else if ( meth(icode) .eq. 3 ) then

        i = 1 
        ii = 1 
        tfun = 0.0D+00 
        hs = h * h 
        as = a * a 
        vi = rrtpi * a * exp ( - 0.5D+00 * ah * ah ) 
        zi = 0.5D+00 * znorm1 ( ah ) / h 
        y = 1.0D+00 / hs 

  310   continue

        tfun = tfun + zi * c2(i) 

        if ( m .lt. i ) then
          tfun = tfun * rrtpi * exp ( - 0.5D+00 * hs ) 
          return 
        end if

        zi = y  * ( dble ( ii ) * zi - vi ) 
        vi = as * vi 
        i = i + 1 
        ii = ii + 2 
        go to 310 
c 
c  t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1 
c  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi) 
c 
      else if ( meth(icode) .eq. 4 ) then

        maxii = m + m + 1
        ii = 1 
        hs = h * h
        as = - a * a 
        tfun = 0.0D+00 
        ai = rtwopi * a * exp ( - 0.5D+00 * hs * ( 1.0D+00 - as ) )  
        yi = 1.0D+00

  410   continue

        tfun = tfun + ai * yi 

        if ( maxii .le. ii ) then
          return
        end if
 
        ii = ii + 2 
        yi = ( 1.0D+00 - hs * yi ) / dble ( ii ) 
        ai = ai * as 
        go to 410 
c 
c  t5(h, a, m) ; m = 13 
c  2m - point gaussian quadrature 
c 
      else if ( meth(icode) .eq. 5 ) then

        tfun = 0.0D+00 
        as = a * a 
        hs = - 0.5D+00 * h * h 
        do i = 1, m 
          r = 1.0D+00 + as * pts(i) 
          tfun = tfun + wts(i) * exp ( hs * r ) / r 
        end do
        tfun = a * tfun 
c 
c  t6(h, a);  approximation for a near 1, (a.le.1) 
c 
      else if ( meth(icode) .eq. 6 ) then

        normh = znorm2 ( h ) 
        tfun = 0.5D+00 * normh * ( 1.0D+00 - normh ) 
        y = 1.0D+00 - a 
        r = atan ( y / ( 1.0D+00 + a ) ) 

        if ( r .ne. 0.0D+00 ) then
          tfun = tfun - rtwopi * r * exp ( - 0.5D+00 * y * h * h / r ) 
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
      function znorm1 ( z )

c*********************************************************************72
c
cc ZNORM1 evaluates the normal CDF from 0 to Z.
c
c  Modified:
c
c    13 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Z, the upper limit.
c
c    Output, double precision ZNORM1, the probability that a standard
c    normal variable will lie between 0 and Z.
c
      implicit none

      double precision z
      double precision znorm1

      znorm1 = 0.5D+00 * erf ( z / sqrt ( 2.0D+00 ) )

      return
      end
      function znorm2 ( z )

c*********************************************************************72
c
cc ZNORM2 evaluates the normal CDF from Z to +oo.
c
c  Modified:
c
c    13 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Z, the lower limit.
c
c    Output, double precision ZNORM2, the probability that a standard
c    normal variable will lie between Z and +oo.
c
      implicit none

      double precision z
      double precision znorm2

      znorm2 = 0.5D+00 * erfc ( z / sqrt ( 2.0D+00 ) )

      return
      end
