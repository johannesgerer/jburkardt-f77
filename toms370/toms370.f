      SUBROUTINE RANDG ( L, X, R, Y )

c*********************************************************************72
C
C  SUBROUTINE RANDG
C
C  PURPOSE:
C
C    COMPUTE RANDOM NUMBERS FROM ANY GENERAL DISTRIBUTION.
C
C  USAGE:
C
C    CALL RANDG ( L, X, R, Y )
C
C  DESCRIPTION OF PARAMETERS:
C
C    INPUT
C
C      L  A NONZERO ODD RANDOM INTEGER.
C
C      X  VECTOR OF LENGTH 257 CONTAINING ORDINATE POINTS
C         SEPARATED BY EQUAL PROBABILITY ON DESIRED DISTRIBUTION.
C         (CAN BE CALCULATED IN RANDGI.)
C
C      R  VECTOR OF LENGTH 256 CONTAINING RATIOS OF DERIVATIVE*DX
C         TO AREA/DX FOR EACH ORDINATE POINT IN X.
C         (CAN BE CALCULATED IN RANDGI.)
C
C    OUTPUT
C
c      L  has been updated so that it may be used as a seed for
c      a subsequent call.
c
C      Y  RANDOM NUMBER.
C
C  REMARKS:
C
C    QUADRATIC APPROXIMATION OF CDF (CUMULATIVE DENSITY FUNCTION)
C    WHICH IMPLIES LINEAR APPROXIMATION OF PDF (PROBABILITY DENSITY
C    FUNCTION).
C
C  SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED:
C
C    NONE DIRECTLY.  RANDGI MAY BE USED FOR INITIALIZATION.
C
C  METHOD:
C
C    TABLE LOOKUP PLUS UNIFORM AND TRIANGULAR DISTRIBUTION
C    VARIABLES ARE USED.
C
      IMPLICIT NONE

      REAL AK2
      INTEGER K1
      INTEGER L
      INTEGER L1
      INTEGER L2
      REAL R(256)
      REAL X(257)
      REAL Y
C
C  GENERATE TWO UNIFORM RANDOM NUMBERS ON INTERVAL ( 1, 2**31 ).
C  ANY GOOD GENERATOR MAY BE SUBSTITUTED.
C
      L1 = IABS ( 65539 * L )
      L = IABS ( 65539 * L1 )
      L2 = L
C
C  CALCULATE TWO UNIFORM RANDOM NUMBERS
C  K1 INTEGER ON INTERVAL ( 1, 256 ),
C  AK2 REAL ON INTERVAL ( 0.0, 1.0 ).
C
      K1 = L1 / 8388608 + 1
      AK2 = FLOAT ( MOD ( L1, 8388608 ) ) * 1.192093E-07
      IF ( AK2 - ABS ( R(K1) ) ) 8, 8, 30
8     IF ( R(K1) ) 20, 10, 10
C
C  CALCULATE TRIANGULAR RANDOM SKEWED LEFT.
C
10    Y = X(K1) + ( X(K1+1) - X(K1) ) * AMAX0 ( L1, L2 ) * 4.656613E-10
      RETURN
C
C  CALCULATE TRIANGULAR RANDOM SKEWED RIGHT.
C
20    Y = X(K1) + ( X(K1+1) - X(K1) ) * AMIN0 ( L1, L2 ) * 4.656613E-10
      RETURN
C
C  CALCULATE UNIFORM RANDOM.
C
30    Y = X(K1) + ( X(K1+1) - X(K1) ) * FLOAT ( L2 ) * 4.656613E-10

      RETURN
      END
      SUBROUTINE RANDGI ( N, X, Y, P, Q, R, K, IER )

c*********************************************************************72
C
C  SUBROUTINE RANDGI
C
C  PURPOSE:
C
C    COMPUTE INITIALIZING VECTORS FOR RANDG.
C
C  USAGE:
C
C    CALL RANDGI ( N, X, Y, P, Q, R, K, IER )
C
C  DESCRIPTION OF PARAMETERS:
C
C    INPUT
C
C      N   NUMBER OF (X,Y) POINTS OF APPROXIMATION TO PDF
C          (PROBABILITY DENSITY FUNCTION) OR CDF
C          (CUMULATIVE DENSITY FUNCTION), N.GE.3.
c
C      X   VECTOR OF LENGTH N CONTAINING ORDINATE OF PDF OR CDF.
c
C      Y   VECTOR OF LENGTH N CONTAINING ABSCISSA OF PDF OR CDF.
C
C    OUTPUT
C
C      P   WORK VECTOR OF LENGTH N.
c
C      Q   VECTOR OR LENGTH 257 CONTAINING ORDINATE POINTS
C          SEPARATED BY EQUAL PROBABILITY ON DESIRED DISTRIBUTION.
C      R   VECTOR OF LENGTH 256 CONTAINING RATIOS OF DERIVATIVES
C          TO AREA/DX FOR EACH ORDINATE POINT IN Q.
c
C      K   K SHOULD BE SET TO 1 IF Y REPRESENTS A CDF.  OTHERWISE,
C          Y WILL BE INTERPRETED AS A PDF.
c
C      IER ERROR INDICATOR
C          1, ERROR IN SCALING.  I.E., TOTAL AREA OF PDF NOT EQUAL TO
C             1.  ASSUMING ESTIMATION ERRORS, A FUDGE FACTOR IS USED
C             TO SCALE A RESULT.
C          2, DENSITY NOT POSITIVE.  I.E., SOME Y(I) .LT. 0.  ABORT.
C          3, NOT IN SORT.  I.E., SOME X(I) .LE. X(I-1).  ABORT.
C          4, SEARCH ERROR.  SHOULD NEVER OCCUR BECAUSE OF FUDGE
C             FACTOR USED.  THIS MEANS SOME P(I) NOT LARGE ENOUGH
C             FOR SEARCH OF PROPER Q.  INVESTIGATION IS NEEDED.
c          5, N < 3.  Abort.
C
C  REMARKS:
C
C    NONE
C
C  SUBROUTINES AND FUNCTION SUBPROGRAMS NEEDED:
C
C    NONE
C
C  METHOD:
C
C    LINEAR APPROXIMATION OF PDF TO FIND CDF (CUMULATIVE DENSITY
C    FUNCTION) AND ICDF (INVERSE CDF).  QUADRATIC INTERPOLATION ON
C    ICDF TO FIND Q AND R.
C
      IMPLICIT NONE

      REAL F
      INTEGER I
      INTEGER IER
      INTEGER J
      INTEGER J1
      INTEGER K
      INTEGER N
      REAL P(300)
      REAL Q(257)
      REAL R(256)
      REAL T1
      REAL T2
      REAL V
      REAL X(300)
      REAL XT1
      REAL XT2
      REAL XT3
      REAL XV1
      REAL XV2
      REAL XV3
      REAL Y(300)
C
C  CALCULATE CUMULATIVE PROBABILITIES.
C
      IER = 0

      if ( n < 3 ) then
        ier = 5
        return
      end if

      IF ( Y(1) ) 5, 10, 10
C
C  ERROR 2.
C
5     IER = 2
      RETURN
10    P(1) = 0.0E+00

      DO 15 I = 2, N

        IF ( Y(I) ) 5, 11, 11
11      IF ( X(I) - X(I-1) ) 6, 6, 12
C
C  ERROR 3.
C
6       IER = 3
        RETURN
12      IF ( K .EQ. 1 ) GO TO 13
C
C  Y IS A PDF.
C
        P(I) = ( Y(I) + Y(I-1) ) * ( X(I) - X(I-1) ) * 0.5E+00 + P(I-1)
        GO TO 15
C
C  Y IS A CDF.
C
13      P(I) = Y(I)

15    CONTINUE

      IF ( P(N) - 0.996094E+00 ) 7, 7, 16
16    IF ( P(N) - 1.003906E+00 ) 3, 7, 7
C
C  ERROR 1.
C
7     IER = 1
3     F = 1.0E+00 / P(N)
      DO 4 I = 2, N
4       P(I) = P(I) * F
C
C  CALCULATE X POINTS FOR EQUAL-DISTANT CUMULATIVE PROBABILITIES.
C
      V = 0.0E+00
      Q(1) = X(1)
      T1 = Y(1)
      J1 = 2
100   DO 150 I = 2, 257
        IF ( I - 257 ) 102, 103, 103
102     V = V + 3.90625E-03
C
C  LOCATE BEST POINT FOR INTERPOLATION.
C
        DO 101 J = J1, N
          IF ( P(J) - V ) 101, 104, 105
101     CONTINUE
C
C  ERROR 4.
C
        IER = 4
103     J = N
104     Q(I) = X(J)
        T2 = Y(J)
        GO TO 125

105     IF ( J .LE. 3 ) GO TO 113
        J1 = J - 2
        GO TO 120
113     J1 = 1
C
C  QUADRATIC INTERPOLATION OF P INVERSE FOR Q.
C
120     XT2 = P(J1+2) - P(J1)
        XT3 = P(J1+2) - P(J1+1)
        XT1 = P(J1+1) - P(J1)
        XV1 = V - P(J1)
        XV2 = V - P(J1+1)
        XV3 = V - P(J1+2)
        Q(I) = ( XV3 * XV2 * X(J1)   ) / ( XT1 * XT2 )
     &       - ( XV3 * XV1 * X(J1+1) ) / ( XT1 * XT3 )
     &       + ( XV2 * XV1 * X(J1+2) ) / ( XT2 * XT3 )
C
C  LINEAR INTERPOLATION OF Y FOR T2 AND R.
C
        T2 = ( Y(J) - Y(J-1) ) * ( Q(I) - X(J-1) )
     &    / ( X(J) - X(J-1) ) + Y(J-1)
125     R(I-1) = ( T2 - T1 ) / ( T2 + T1 )
        T1 = T2
        J1 = J
150   CONTINUE

      RETURN
      END
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
