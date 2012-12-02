      SUBROUTINE CHARAC ( DATA, M, IFAIL )

c*********************************************************************72
C
C  THIS ROUTINE ADVANCES ONE GRID STEP THE SOLUTION OF THE SYSTEM
C
C    A1 * U(X) + A2 * U(Y) + A3 * V(X) + A4 * V(Y) = H1
C    B1 * U(X) + B2 * U(Y) + B3 * V(X) + B4 * V(Y) = H2
C
C  WHERE U(X) MEANS PARTIAL DERIVATIVE OF U WITH RESPECT TO X, AND
C  SO ON, AND A1 = A1(X,Y,U,V), AND SO ON.
C
C  THE INITIAL DATA IS GIVEN IN THE MATRIX "DATA", EACH COLUMN OF
C  FOUR ELEMENTS CONTAINING A VALUE OF X, Y, U, V.
C
C  M IS THE NUMBER OF DATA POINTS ON THE INITIAL CURVE.
C
      IMPLICIT NONE

      INTEGER M

      REAL D0(4)
      REAL D1(4)
      REAL D2(4)
      REAL D3(4)
      REAL D4(4)
      REAL DATA(4,M)
      REAL DX
      REAL DX1
      REAL DX2
      REAL DY
      INTEGER FAIL
      INTEGER FLAG
      INTEGER IFAIL
      INTEGER I
      INTEGER J
      REAL S1(4)
      REAL S2(4)
      REAL S3(4)
      REAL T1(4)
      REAL T2(4)
      REAL T3(4)
      REAL T4(4)
      REAL T5(4)
      REAL T6(4)
      REAL T7(4)
      REAL V1(8)
      REAL V2(8)
      REAL V3(8)
      REAL V4(8)
      REAL V5(8)
      REAL V6(8)
      REAL V7(8)
      REAL V8(8)
      REAL V9(8)
      REAL V10(8)

      COMMON / CHFAIL / FAIL, FLAG

      FAIL = 0
      M = M - 4
      IF ( M .LE. 0 ) RETURN

      DO 145 J = 1, 4
        D1(J) = DATA(J,1)
        D2(J) = DATA(J,2)
        D3(J) = DATA(J,3)
145     D4(J) = DATA(J,4)

      CALL CHVAR ( D1, V2 )
      IF ( IFAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( D2, V3 )
      IF ( IFAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( D3, V4 )
      IF ( IFAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( D4, V5 )
      IF ( IFAIL .NE. 0 ) GO TO 250

      FLAG = 0
      DY = D2(2) - D1(2)
      DX = D2(1) - D1(1)
      IF ( DY ) 148, 149, 150
148   DY = -DY
      DX = -DX
      GO TO 150
149   DY = 1.0E+00
      DX = 1.0E+30
150   DX1 = DY / V2(1)
      DX2 = DY / V2(2)
      IF ( DX1 .LT. DX .AND. DX2 .GE. DX ) GO TO 170
      IF ( DX2 .LT. DX .AND. DX1 .GE. DX ) GO TO 175
      IF ( DX2 .GE. DX1 ) GO TO 175

170   FLAG = 1
175   CONTINUE

      CALL CHSTEP ( D1, V2, D3, V4, T1 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T1, V6 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( D2, V3, D4, V5, T6 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T6, V10 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( D3, V4, D4, V5, T2 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T2, V7 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( D2, V3, D3, V4, T4 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T4, V9 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( T4, V9, T2, V7, T3 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T3, V8 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( D1, V2, D2, V3, T5 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T5, V1 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( T5, V1, T4, V9, T4 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T4, V1 )
      IF ( FAIL .NE. 0 ) GO TO 250

      CALL CHSTEP ( T4, V1, T3, V8, T4 )
      IF ( FAIL .NE. 0 ) GO TO 250
      CALL CHVAR ( T4, V9 )
      IF ( FAIL .NE. 0 ) GO TO 250

      DO 100 I = 1, M

        DO 200 J = 1, 8
          V1(J) = V2(J)
          V2(J) = V3(J)
          V3(J) = V4(J)
          V4(J) = V5(J)
200     CONTINUE

        DO 201 J = 1, 4
          D0(J) = D1(J)
          D1(J) = D2(J)
          D2(J) = D3(J)
          D3(J) = D4(J)
201       D4(J) = DATA(J,I+4)

        CALL CHVAR ( D4, V5 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHSTEP ( D0, V1, D4, V5, S1 )
        IF ( FAIL .NE. 0 ) GO TO 250
        CALL CHSTEP ( D2, V3, D4, V5, T5 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHVAR ( T5, V1 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHSTEP ( T1, V6, T5, V1, S2 )
        IF ( FAIL .NE. 0 ) GO TO 250
        CALL CHSTEP ( D3, V4, D4, V5, T1 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHVAR ( T1, V6 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHSTEP ( T2, V7, T1, V6, T7 )
        IF ( FAIL .NE. 0 ) GO TO 250

        DO 210 J = 1, 4
          T2(J)    = T1(J)
          V7(J)    = V6(J)
          V7(J+4)  = V6(J+4)
          T1(J)    = T6(J)
          V6(J)    = V10(J)
          V6(J+4)  = V10(J+4)
          T6(J)    = T5(J)
          V10(J)   = V1(J)
          V10(J+4) = V1(J+4)
210     CONTINUE

        CALL CHVAR ( T7, V1 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHSTEP ( T3, V8, T7, V1, T5 )
        IF ( FAIL .NE. 0 ) GO TO 250

        DO 220 J = 1, 4
          T3(J) = T7(J)
          V8(J) = V1(J)
          V8(J+4) = V1(J+4)
220     CONTINUE

        CALL CHVAR ( T5, V1 )
        IF ( FAIL .NE. 0 ) GO TO 250

        CALL CHSTEP ( T4, V9, T5, V1, S3 )
        IF ( FAIL .NE. 0 ) GO TO 250

        DO 230 J = 1, 4
          T4(J) = T5(J)
          V9(J) = V1(J)
          V9(J+4) = V1(J+4)
230     CONTINUE
C
C  EXTRAPOLATE THE THREE SUCCESSIVE APPROXIMATIONS.
C
        DO 300 J = 1, 4
300       DATA(J,I) = 0.3333333333E+00 
     &      * ( 8.0E+00 * S3(J) - 6.0E+00 * S2(J) + S1(J) )

100   CONTINUE

      IFAIL = 0
      RETURN
C
C  ERROR EXIT.
C
250   IFAIL = FAIL
      RETURN
      END
      SUBROUTINE CHVAR ( XYUV, VAR )

c*********************************************************************72
C
C  COMPUTES THE VALUES SIGMA1, SIGMA2, R1, R2, S1, S2, T1, T2,
C  STORED IN THE LIST VAR, FROM THE COEFFICIENT FUNCTIONS AND
C  THE VALUES X, Y, U, V IN THE LIST XYUV.
C
      IMPLICIT NONE

      REAL A
      REAL B
      REAL C
      REAL D
      INTEGER I
      INTEGER FAIL
      INTEGER FLAG
      REAL T(10)
      REAL VAR(8)
      REAL XYUV(4)

      COMMON / CHFAIL / FAIL, FLAG

      CALL CHCOEF ( T, XYUV )
C
C  COEFFICIENTS OF SYSTEM ARE STORED IN THE LIST T.
C
      A = T(1) * T(8) - T(3) * T(6)
      B = 0.5E+00 
     &  * ( T(1) * T(9) - T(4) * T(6) - T(3) * T(7) + T(2) * T(8) )
      C = T(2) * T(9) - T(4) * T(7)

      IF ( A .NE. 0.0E+00 ) GO TO 150
      IF ( B .EQ. 0.0E+00 ) GO TO 500
      VAR(1) = 1.0E+15
      VAR(2) = 0.5E+00 * C / B
      IF ( B .GT. 0.0E+00 ) GO TO 400
      VAR(1) = VAR(2)
      VAR(2) = 1.0E+15
      GO TO 400

150   D = B * B - A * C
      IF ( D .LE. 0.0E+00 ) GO TO 500
      D = SQRT ( D )
      IF ( B .LT. 0.0E+00 ) GO TO 300
      VAR(1) = ( B + D ) / A
      VAR(2) = C / ( A * VAR(1) )
      GO TO 400

300   VAR(2) = ( B - D ) / A
      VAR(1) = C / ( A * VAR(2) )

400   DO 100 I = 1, 2
        T(4) = T(1) * VAR(I) - T(2)
        T(9) = T(6) * VAR(I) - T(7)
        VAR(I+2) = T(1) * T(9) - T(6) * T(4)
        VAR(I+4) = T(3) * T(9) - T(8) * T(4)
        VAR(I+6) = T(5) * T(9) - T(10) * T(4)
100   CONTINUE

      RETURN
C
C  ERROR EXIT.
C
500   FAIL = 1
      RETURN
      END
      SUBROUTINE CHSTEP ( D1, V1, D2, V2, D3 )

c*********************************************************************72
C
C  THIS ROUTINE COMPUTES THE VALUES OF X3, Y3, U3, V3 AT THE POINT
C  DETERMINED BY THE INTERSECTION OF THE CHARACTERISTICS THROUGH (X1,Y1)
C  AND (X2,Y2).
C
      IMPLICIT NONE

      REAL AA
      REAL BB
      REAL D1(4)
      REAL D2(4)
      REAL D3(4)
      REAL D31
      REAL D32
      REAL DEM1
      REAL DEM2
      REAL EPS
      INTEGER FAIL
      INTEGER FLAG
      REAL R1
      REAL R2
      REAL RHO
      REAL S1
      REAL S2
      REAL SIG
      REAL T1
      REAL T2
      REAL TA
      REAL TB
      REAL TC
      REAL TD
      REAL V1(8)
      REAL V2(8)

      COMMON / CHFAIL / FAIL, FLAG

      EPS = 1.0E-05
      IF ( FLAG .NE. 0 ) GO TO 150
100   SIG = V1(2)
      R1 = V1(4)
      S1 = V1(6)
      T1 = V1(8)
      GO TO 180

150   SIG = V1(1)
      R1 = V1(3)
      S1 = V1(5)
      T1 = V1(7)

180   CONTINUE

      IF ( FLAG .NE. 0 ) GO TO 250

200   RHO = V2(1)
      R2 = V2(3)
      S2 = V2(5)
      T2 = V2(7)
      GO TO 280

250   RHO = V2(2)
      R2 = V2(4)
      S2 = V2(6)
      T2 = V2(8)

280   CONTINUE
C
C  COMPUTE X3, Y3, U3, V3.
C
      DEM1 = SIG - RHO
      IF ( ABS ( DEM1 ) .LT. EPS * ABS ( SIG ) ) GO TO 900
      AA = D1(2) - SIG * D1(1)
      BB = D2(2) - RHO * D2(1)
      D32 = ( SIG * BB - RHO * AA ) / DEM1
      D31 = ( BB - AA ) / DEM1
      TA = T1 * ( D31 - D1(1) ) + R1 * D1(3) + S1 * D1(4)
      TB = T2 * ( D31 - D2(1) ) + R2 * D2(3) + S2 * D2(4)
      TC = R1 * S2
      TD = R2 * S1
      DEM2 = TC - TD
      TC = AMAX1 ( ABS ( TC ), ABS ( TD ) )
      IF ( ABS ( DEM2 ) .LE. EPS * TC ) GO TO 900
      D3(3) = ( TA * S2 - TB * S1 ) / DEM2
      D3(4) = ( TB * R1 - TA * R2 ) / DEM2
      D3(2) = D32
      D3(1) = D31
      RETURN

900   FAIL = 2
      RETURN
      END
