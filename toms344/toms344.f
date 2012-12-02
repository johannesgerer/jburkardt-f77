      SUBROUTINE TTEST ( T, DF, ANS, KERR )

c*********************************************************************72
c
cc TTEST
c
      IMPLICIT NONE

      REAL ANS
      REAL D1
      REAL D2
      INTEGER DF
      REAL F1
      REAL F2
      INTEGER I
      INTEGER KERR
      INTEGER N
      REAL T
      REAL T1
      REAL T2

      SAVE D1

      DATA D1 / 0.63661977236758134E+00 /
C
C  0.63661977236758134... = 2 / PI
C
      KERR = 0
C
      IF ( DF .GT. 0 ) GO TO 1
C
C  ERROR RETURN IF DF NOT POSITIVE.
C
      KERR = 1
      ANS = 0.0E+00
      RETURN
C
C  BEGIN COMPUTATION OF SERIES.
C
1     T2 = T * T / FLOAT ( DF )
      T1 = SQRT ( T2 )
      T2 = 1.0E+00 / ( 1.0E+00 + T2 )
C
      IF ( ( DF / 2 ) * 2 .EQ. DF ) GO TO 5
C
C  DF IS AN ODD INTEGER.
C
      ANS = 1.0E+00 - D1 * ATAN ( T1 )
C
      IF ( DF .EQ. 1 ) GO TO 4
C
      D2 = D1 * T1 * T2
      ANS = ANS - D2
C
      IF ( DF .EQ. 3 ) GO TO 4
C
      F1 = 0.0E+00
2     N = ( DF - 2 ) / 2
      DO 3 I = 1, N
        F2 = 2.0E+00 * FLOAT ( I ) - F1
        D2 = D2 * T2 * F2 / ( F2 + 1.0E+00 )
3       ANS = ANS - D2
C
C  COMMON RETURN AFTER COMPUTATION.
C
4      IF ( ANS .LT. 0.0E+00 ) ANS = 0.0E+00
       RETURN
C
C  DF IS AN EVEN INTEGER.
C
5      D2 = T1 * SQRT ( T2 )
       ANS = 1.0E+00 - D2
C
       IF ( DF .EQ. 2 ) GO TO 4
C
       F1 = 1.0E+00
       GO TO 2
       END
