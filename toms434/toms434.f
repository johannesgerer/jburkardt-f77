      SUBROUTINE CONP ( MATRIX, NR, NC, PT, PS, PC )

c*********************************************************************72
C
C  INPUT ARGUMENTS.
C
C    MATRIX = SPECIFICATION OF THE CONTINGENCY TABLE.
C      THIS MATRIX IS PARTIONED AS FOLLOWS:
C
C      X(11).....X(IC)  R(1)
C        .  .....  .      .
C        .  .....  .      .
C      X(R1).....X(RC)  R(R)
C       C(1)..... C(C)    N
C
C      WHERE X(IJ) ARE THE OBSERVED CELL FREQUENCIES,
C      R(I) ARE THE ROW TOTALS, C(J) ARE THE COLUMN
C      TOTALS, AND N IS THE TOTAL SAMPLE SIZE.
C      NOTE THAT THE ORIGINAL CELL FREQUENCIES ARE
C      DESTROYED BY THIS SUBROUTINE.
C
C    NR = THE NUMBER OF ROWS IN MATRIX (R=NR-1).
C
C    NC = THE NUMBER OF COLUMNS IN MATRIX (C=NC-1).
C
C  OUTPUT ARGUMENTS.
C
C    PT = THE PROBABILITY OF OBTAINING THE GIVEN TABLE.
C
C    PS = THE PROBABILITY OF OBTAINING A TABLE AS PROBABLE
C      AS, OR LESS PROBABLE THAN, THE GIVEN TABLE.
C
C    PC = THE PROBABILITY OF OBTAINING SOME OF THE
C      TABLES POSSIBLE WITHIN THE CONSTRAINTS OF THE
C      MARGINAL TOTALS.  (THIS SHOULD BE 1.0.  DEVIATIONS
C      FROM 1.0 REFLECT THE ACCURACY OF THE COMPUTATION.)
C
C  EXTERNALS
C
C    FACLOG(N) = FUNCTION TO RETURN THE FLOATING POINT
C      VALUE OF LOG BASE 10 OF N FACTORIAL.
C
      DIMENSION MATRIX(NR,NC)
      INTEGER R,C,TEMP
C
      R = NR - 1
      C = NC - 1
C
C  COMPUTE LOG OF CONSTANT NUMERATOR.
C
      QXLOG = -FACLOG ( MATRIX(NR,NC) )
      DO 10 I = 1, R
10      QXLOG = QXLOG + FACLOG ( MATRIX(I,NC) )
      DO 20 J = 1, C
20      QXLOG = QXLOG + FACLOG ( MATRIX(NR,J) )
C
C  COMPUTE PROBABILITY OF GIVEN TABLE.
C
      RXLOG = 0.0
      DO 50 I = 1, R
        DO 50 J = 1, C
50        RXLOG = RXLOG + FACLOG ( MATRIX(I,J) )
      PT = 10.0**( QXLOG - RXLOG )
C
      PS = 0.0
      PC = 0.0
C
C  FILL LOWER RIGHT (R-1) X (C-1) CELLS WITH
C  MINIMIMUM OF ROW AND COLUMN TOTALS.
C
      DO 100 I = 2, R
        DO 100 J = 2, C
100       MATRIX(I,J) = MIN0 ( MATRIX(I,NC), MATRIX(NR,J) )
      GO TO 300
C
C  OBTAIN A NEW SET OF FREQUENCIES IN
C  LOWER RIGHT (R-1) X (C-1) CELLS.
C
200   DO 220 I = 2, R
        DO 220 J = 2, C
          MATRIX(I,J) = MATRIX(I,J) - 1
          IF ( MATRIX(I,J) .GE. 0 ) GO TO 300
220       MATRIX(I,J) = MIN0 ( MATRIX(I,NC), MATRIX(NR,J) )
      RETURN
C
C  FILL REMAINDER OF OBSERVED CELLS
C  ...COMPLETE COLUMN 1.
C
300   DO 320 I = 2, R
        TEMP = MATRIX(I,NC)
        DO 310 J = 2, C
310       TEMP = TEMP - MATRIX(I,J)
        IF ( TEMP .LT. 0 ) GO TO 200
320     MATRIX(I,1) = TEMP
C
C  ...COMPLETE ROW 1.
C
      DO 340 J = 1, C
        TEMP = MATRIX(NR,J)
        DO 330 I = 2, R
330       TEMP = TEMP - MATRIX(I,J)
        IF ( TEMP .LT. 0 ) GO TO 200
340     MATRIX(1,J) = TEMP
C
C  COMPUTE LOG OF THE DENOMINATOR.
C
      RXLOG = 0.0
      DO 350 I = 1, R
        DO 350 J = 1, C
350       RXLOG = RXLOG + FACLOG ( MATRIX(I,J) )
C
C  COMPUTE PX.  ADD TO PS IF PX .LE. PT
C  (ALLOW FOR ROUND-OFF ERROR).
C
      PX = 10.0** ( QXLOG - RXLOG )
      PC = PC + PX
      IF ( ( PT / PX ) .GT. 0.99999 ) PS = PS + PX
      GO TO 200
      END
      FUNCTION FACLOG ( N )

c*********************************************************************72
C
C  INPUT ARGUMENT.
C
C    N = AN INTEGER GREATER THAN OR EQUAL TO ZERO.
C
C  FUNCTION RESULT.
C
C    FACLOG = THE LOG TO THE BASE 10 OF N FACTORIAL.
C
      DIMENSION TABLE(101)
      SAVE ELOG,IFLAG,TABLE,TPILOG
      DATA TPILOG / 0.3990899342 /
      DATA ELOG / 0.4342944819 /
      DATA IFLAG / 0 /
C
C  USE STIRLINGS APPROXIMATION IF N GT 100.
C
      IF ( N .GT. 100 ) GO TO 50
C
C  LOOKUP UP ANSWER IF TABLE WAS GENERATED.
C
      IF ( IFLAG .EQ. 0 ) GO TO 100
10    FACLOG = TABLE(N+1)
      RETURN
C
C  HERE FOR STIRLINGS APPROXIMATION.
C
50    X = FLOAT ( N )
      FACLOG = ( X + 0.5 ) * ALOG10 ( X ) - X * ELOG + TPILOG
     &  + ELOG / ( 12.0 * X ) - ELOG / ( 360.0 * X * X * X )
      RETURN
C
C  HERE TO GENERATE LOG FACTORIAL TABLE.
C
100   TABLE(1) = 0.0
      DO 120 I = 2, 101
        X = FLOAT ( I - 1 )
120     TABLE(I) = TABLE(I-1) + ALOG10 ( X )
      IFLAG = 1
      GO TO 10
      END
      function factorial ( n )

c*********************************************************************72
c
cc FACTORIAL evaluates the factorial function.
c
      real factorial

      factorial = 1.0E+00
      do i = 1, n
        factorial = factorial * i
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
