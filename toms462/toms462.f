      DOUBLE PRECISION FUNCTION BIVNOR ( AH, AK, R )

c*********************************************************************72
c
cc BIVNOR computes the bivariate normal CDF.
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
c    Input, double precision R, the correlation between X and Y.
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
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
