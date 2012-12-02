      SUBROUTINE MFCVAL ( N, R, QQ, CV, J )

c*********************************************************************72
C
CC MFCVAL computes characteristic values of Mathieu's differential equation.
C
C  Modified:
C
C    06 October 2006
C
C  Author:
C
C    Donald Clemm
C
C  Reference:
C
C    Donald Clemm,
C    Algorithm 352: Characteristic Values and Associated 
C    Solutions of Matthieu's Differential Equation,
C    Communications of the ACM,
C    Volume 12, Number 7, pages 399-407, June 1969.
C
C  Parameters:
C
C    Input, integer N, the number of characteristic values desired.
C
C    Input, integer R, given as N-1 or N according as the characteristic
C    values are to be associated with the even or odd solutions, respectively.
C
C    Input, double precision QQ, the nonnegative parameter.
C
C    Output, double precision CV(6,N), the computed array of characteristic
C    values and bounds.
C
C    Output, integer J, the number of characteristic values successfully
C    computed.  If J is not equal to N then an error occurred while trying 
C    to compute value J+1.
C
      IMPLICIT NONE

      INTEGER N

      DOUBLE PRECISION A
      DOUBLE PRECISION A0
      DOUBLE PRECISION A1
      DOUBLE PRECISION A2
      DOUBLE PRECISION CV(6,N)
      DOUBLE PRECISION DL
      DOUBLE PRECISION DR
      DOUBLE PRECISION DTM
      INTEGER J
      INTEGER K
      INTEGER KK
      INTEGER L
      INTEGER M
      INTEGER M0
      INTEGER M1
      INTEGER M2S
      INTEGER MF
      DOUBLE PRECISION Q
      DOUBLE PRECISION QQ
      INTEGER R
      DOUBLE PRECISION T
      DOUBLE PRECISION TM
      DOUBLE PRECISION TOL
      DOUBLE PRECISION TOLA
      INTEGER TYPE

      EQUIVALENCE ( DL, DR, T )

      COMMON / MF1 / Q, TOL, TYPE, M1, M0, M2S, MF

      COMMON / MF2 / A0, A2, A1

      TOL = 1.0D-13

      IF ( N - R ) 10, 10, 20

10    L = 1
      GO TO 30

20    L = 2

30    Q = QQ

      DO 500 K = 1, N

        J = K

        IF ( Q .EQ. 0.0D+00 ) THEN
          GO TO 490
        END IF

        IF ( Q .LT. 0.0D+00 ) THEN
960       WRITE ( 6, 961 )
          Q = -Q
        END IF

40      KK = MIN0 ( K, 4 )
        TYPE = 2 * MOD ( L, 2 ) + MOD ( K - L + 1, 2 )
C
C  FIRST APPROXIMATION.
C
        GO TO ( 100, 200, 300, 400 ), KK

100     IF ( Q - 1.0D+00 ) 110, 140, 140

110     GO TO ( 120, 130 ), L

120     A = 1.0D+00 - Q - 0.125D+00 * Q * Q
        GO TO 420

130     A = Q * Q
        A = A * ( -0.5D+00 + 0.0546875D+00 * A )
        GO TO 420

140     IF ( Q - 2.0D+00 ) 150, 180, 180

150     GO TO ( 160, 170 ), L

160     A = 1.033D+00 - 1.0746D+00 * Q - 0.0688D+00 * Q * Q
        GO TO 420

170     A = 0.23D+00 - 0.495D+00 * Q - 0.191D+00 * Q * Q
        GO TO 420

180     A = -0.25D+00 - 2.0D+00 * Q + 2.0D+00 * DSQRT ( Q )
        GO TO 420

200     DL = L
        IF ( Q * DL - 6.0D+00 ) 210, 350, 350

210     GO TO ( 220, 230 ), L

220     A = 4.01521D+00 - Q * ( 0.046D+00 + 0.0667857D+00 * Q )
        GO TO 420

230     A = 1.0D+00 + 1.05007D+00 * Q - 0.180143D+00 * Q * Q
        GO TO 420

300     IF ( Q - 8.0D+00 ) 310, 350, 350

310     GO TO ( 320, 330 ), L

320     A = 8.93867D+00 + 0.178156D+00 * Q - 0.0252132D+00 * Q * Q
        GO TO 420

330     A = 3.70017D+00 + 0.953485D+00 * Q - 0.0475065D+00 * Q * Q
        GO TO 420

350     DR = K - 1
        A = CV(1,K-1) - DR + 4.0D+00 * DSQRT ( Q )
        GO TO 420

400     A = CV(1,K-1) - CV(1,K-2)
        A = 3.0D+00 * A + CV(1,K-3)

420     IF ( Q .GE. 1.0D+00 ) GO TO 440
        IF ( K .NE. 1 ) GO TO 430

425     TOLA = DMAX1 ( DMIN1 ( TOL, DABS ( A ) ), 1.0D-14 )
        GO TO 450

430     TOLA = TOL * DABS ( A )
        GO TO 450

440     TOLA = TOL * DMAX1 ( Q, DABS ( A ) )

445     TOLA = DMAX1 ( 
     &    DMIN1 ( TOLA, DABS ( A ), 0.4D+00 * DSQRT ( Q ) ), 
     &    1.0D-14 )
C
C  CRUDE UPPER AND LOWER BOUNDS.
C
450     CONTINUE

        CALL BOUNDS ( K, A, TOLA, CV, N, M )
        IF ( M .NE. 0 ) IF ( M - 1 ) 470, 910, 900
C
C  ITERATE.
C
        CALL MFITR8 ( TOLA, CV(1,K), CV(2,K), M )
        IF ( M .GT. 0 ) GO TO 920
C
C  FINAL BOUNDS AND FUNCTIONS.
C
470     T = CV(1,K) - TOLA
        CALL TMOFA ( T, TM, DTM, M )

        IF ( M .GT. 0 ) THEN
          WRITE ( 6, 941 ) K
          CV(3,K) = 0.0D+00
          CV(4,K) = 0.0D+00
        ELSE
          CV(3,K) = T
          CV(4,K) = -TM / DTM
        END IF

480     T = CV(1,K) + TOLA
        CALL TMOFA ( T, TM, DTM, M )
        IF ( M .GT. 0 ) GO TO 950
        CV(5,K) = T
        CV(6,K) = -TM / DTM
        GO TO 500
C
C  Q EQUALS ZERO.
C
490     CV(1,K) = ( K - L + 1 )**2
        CV(2,K) = 0.0D+00
        CV(3,K) = CV(1,K)
        CV(4,K) = 0.0D+00
        CV(5,K) = CV(1,K)
        CV(6,K) = 0.0D+00

500   CONTINUE

550   RETURN
C
C  PRINT ERROR MESSAGES.
C
900   WRITE ( 6, 901 ) K
901   FORMAT ( 20H0CRUDE BOUNDS CANNOT,
     &  22H BE LOCATED, NO OUTPUT,
     &  7H FOR K=, I2 )
      GO TO 930

910   WRITE ( 6, 911 ) K
911   FORMAT ( ' ERROR IN SUBPROGRAM TMOFA, VIA SUBPROGRAM',
     &  ' BOUNDS, NO OUTPUT, FOR K = ', I2 )
      GO TO 930

920   WRITE ( 6, 921 ) K
921   FORMAT ( ' ERROR IN SUBPROGRAM TMOFA, VIA SUBPROGRAM',
     &  ' MFITR8, NO OUTPUT FOR K = ', I2 )

930   J = J - 1
      GO TO 550
941   FORMAT ( ' ERROR IN SUBPROGRAM TMOFA, NO LOWER BOUND',
     &  ' FOR K = ', I2 )

950   WRITE ( 6, 951 ) K
951   FORMAT ( ' ERROR IN SUBPROGRAM TMOFA, NO UPPER BOUND',
     &  ' FOR K = ', I2 )
      CV(5,K) = 0.0D+00
      CV(6,K) = 0.0D+00
      GO TO 500
961   FORMAT ( 20H0Q GIVEN NEGATIVELY,,
     &  20H USED ABSOLUTE VALUE )
      END
      SUBROUTINE MATH ( XX, QQ, R, CV, SOL, FNC, NORM, F, K, M )
 
c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A
      DOUBLE PRECISION AB(200)
      DOUBLE PRECISION CV
      EXTERNAL DC
      DOUBLE PRECISION DC
      EXTERNAL DDC
      DOUBLE PRECISION DDC
      EXTERNAL DDS
      DOUBLE PRECISION DDS
      DOUBLE PRECISION DLAST
      DOUBLE PRECISION DMAX
      EXTERNAL DPC
      DOUBLE PRECISION DPC
      EXTERNAL DPS
      DOUBLE PRECISION DPS
      EXTERNAL DS
      DOUBLE PRECISION DS
      REAL DUM1(578)
      REAL DUM2(6)
      DOUBLE PRECISION F(3)
      INTEGER FNC
      DOUBLE PRECISION G
      INTEGER I
      DOUBLE PRECISION J(250)
      INTEGER K(2)
      INTEGER KLAST
      INTEGER KMAX
      INTEGER L
      INTEGER LL
      INTEGER M
      INTEGER MF
      INTEGER ML
      INTEGER MM
      INTEGER M0
      INTEGER M1
      INTEGER M2S
      INTEGER N
      INTEGER NORM
      INTEGER P
      EXTERNAL PC
      EXTERNAL PS
      DOUBLE PRECISION Q
      DOUBLE PRECISION QQ
      INTEGER R
      INTEGER S
      INTEGER SOL
      DOUBLE PRECISION T
      DOUBLE PRECISION TOL
      INTEGER TYPE
      DOUBLE PRECISION U1
      DOUBLE PRECISION U2
      DOUBLE PRECISION X
      DOUBLE PRECISION XX
      DOUBLE PRECISION Y(250)

      COMMON J, Y, U1, U2, N, P, S, L, X, T, I, LL, G, DMAX, DLAST,
     &  KMAX, KLAST, DUM1, A, DUM2, MM, ML, AB

      COMMON / MF1 / Q, TOL, TYPE, M1, M0, M2S, MF

      M = 0
      IF ( SOL .LT. 1 .OR. 
     &  SOL .GT. 3 .OR. 
     &  FNC .LT. 1 .OR. 
     &  FNC .GT. 4 ) GO TO 400

      A = CV
      Q = QQ
      TOL = 1.0D-13
      TYPE = 2 * MOD ( FNC, 2 ) + MOD ( R, 2 )
      CALL COEF ( M )

      IF ( M ) 410, 10, 420

10    N = R / 2
      P = MOD ( R, 2 )
      S = MM / 2
      L = ML / 2
      X = XX
      T = 1.0D+00

      IF ( SOL .EQ. 3 ) GO TO ( 150, 160, 170, 180 ), FNC

      U1 = DSQRT ( Q ) * DEXP ( -X )
      U2 = Q / U1
      LL = L + S + P
C
C  COMPUTE BESSEL FUNCTIONS.
C
      CALL BESSEL ( 1, U1, J, LL )
      CALL BESSEL ( SOL, U2, Y, LL )
C
C  EVALUATE SELECTED FUNCTION.
C
      GO TO ( 50, 60, 70, 80 ), FNC

50    CALL SUM_SERIES ( DS )
      GO TO 300
60    CALL SUM_SERIES ( DC )
      GO TO 300
70    CALL SUM_SERIES ( DDS )
      GO TO 300
80    CALL SUM_SERIES ( DDC )
      GO TO 300
150   CALL SUM_SERIES ( PS )
      GO TO 200
160   CALL SUM_SERIES ( PC )
      GO TO 200
170   CALL SUM_SERIES ( DPS )
      GO TO 200
180   CALL SUM_SERIES ( DPC )

200   IF ( NORM - 2 ) 300, 210, 250
C
C  INCE NORMALIZATION.
C
210   T = AB(1)**2
      IF ( TYPE .EQ. 0 ) T = T + T
      DO 220 I = 1, L
        T = T + AB(I+1)**2
220   CONTINUE
      T = DSQRT ( T )
      I = M0 / 2
      IF ( AB(I) .LT. 0.0D+00 ) T = -T
      GO TO 300
C
C  STRATTON NORMALIZATION.
C
250   IF ( TYPE .GT. 1 ) GO TO 270
      T = AB(1)
      DO 260 I = 1, L
        T = T + AB(I+1)
260   CONTINUE
      GO TO 300
270   T = DBLE ( FLOAT ( P ) ) * AB(1)
      DO 280 I = 1, L
        T = T + AB(I+1) * DBLE ( FLOAT ( 2 * I + P ) )
280   CONTINUE

300   F(1) = G / T
      F(2) = DMAX / T
      F(3) = DLAST / T
      K(1) = KMAX
      K(2) = KLAST
350   RETURN
C
C  PRINT ERROR MESSAGES.
C
400   WRITE ( 6, 401 )
401   FORMAT ( 18H0SOL OR FNC OUT OF,
     &  17H RANGE, NO OUTPUT )
      GO TO 450
410   WRITE ( 6, 411 )
411   FORMAT ( 15H0MORE THAN 200 ,
     &  22HCOEFFICIENTS REQUIRED,,
     &  20H QQ AND R TOO LARGE,,
     &  10H NO OUTPUT )
      GO TO 450
420   WRITE ( 6, 421 )
421   FORMAT ( 20H0ERROR IN SUBPROGRAM,
     &  22H TMOFA, VIA SUBPROGRAM,
     &  13H COEF, VERIFY,
     &  21H ARGUMENTS, NO OUTPUT )

450   M = 1
      F(1) = 0.0D+00
      F(2) = 0.0D+00
      F(3) = 0.0D+00
      K(1) = 0
      K(2) = 0
      GO TO 350
      END
      SUBROUTINE BOUNDS ( K, APPROX, TOLA, CV, N, MM )

c*********************************************************************72

      IMPLICIT NONE

      INTEGER N

      DOUBLE PRECISION A
      DOUBLE PRECISION A0
      DOUBLE PRECISION A1
      DOUBLE PRECISION APPROX
      DOUBLE PRECISION CV(6,N)
      DOUBLE PRECISION D0
      DOUBLE PRECISION D1
      DOUBLE PRECISION DTM
      INTEGER K
      INTEGER KA
      INTEGER M
      INTEGER M0
      INTEGER M1
      INTEGER M2S
      INTEGER MF
      INTEGER MM
      DOUBLE PRECISION Q
      DOUBLE PRECISION TM
      DOUBLE PRECISION TOL
      DOUBLE PRECISION TOLA
      INTEGER TYPE

      COMMON / MF1 / Q, TOL, TYPE, M1, M0, M2S, MF

      COMMON / MF2 / A0, A, A1

      KA = 0
      IF ( K .EQ. 1 ) GO TO 20

      IF ( APPROX - CV(1,K-1) ) 10, 10, 20

10    A0 = CV(1,K-1) + 1.0D+00
      GO TO 30

20    A0 = APPROX

30    CALL TMOFA ( A0, TM, DTM, M )
      IF ( M .GT. 0 ) GO TO 250
      D0 = -TM / DTM
      IF ( D0 ) 100, 300, 50
C
C  A0 IS LOWER BOUND, SEARCH FOR UPPER BOUND.
C
50    A1 = A0 + D0 + 0.1D+00
      CALL TMOFA ( A1, TM, DTM, M )

      IF ( M .GT. 0 ) GO TO 250
      D1 = -TM / DTM

      IF ( D1 ) 200, 350, 60

60    A0 = A1
      D0 = D1
      KA = KA + 1

      IF ( KA - 4 ) 50, 400, 400
C
C  A1 IS UPPER BOUND, SEARCH FOR LOWER BOUND.
C
100   A1 = A0
      D1 = D0
      A0 = DMAX1 ( A1 + D1 - 0.1D+00, -2.0D+00 * Q )
      IF ( K .EQ. 1 ) GO TO 110
      IF ( A0 - CV(1,K-1) ) 150, 150, 110

110   CALL TMOFA ( A0, TM, DTM, M )
      IF ( M .GT. 0 ) GO TO 250
      D0 = -TM / DTM
      IF ( D0 ) 120, 300, 200

120   KA = KA + 1
      IF ( KA - 4 ) 100, 400, 400

150   KA = KA + 1
      IF ( KA - 4 ) 160, 400, 400

160   A0 = A1 + DMAX1 ( TOLA, DABS ( D1 ) )
      GO TO 30

200   A = 0.5D+00 * ( A0 + D0 + A1 + D1 )
      IF ( A .LE. A0 .OR. A .GE. A1 ) A = 0.5D+00 * ( A0 + A1 )

250   MM = M
      RETURN

300   CV(1,K) = A0

310   CV(2,K) = 0.0D+00
      M = -1
      GO TO 250

350   CV(1,K) = A1
      GO TO 310

400   M = 2
      GO TO 250
      END
      SUBROUTINE MFITR8 ( TOLA, CV, DCV, MM )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A
      DOUBLE PRECISION A0
      DOUBLE PRECISION A1
      DOUBLE PRECISION A2
      DOUBLE PRECISION CV
      DOUBLE PRECISION D
      DOUBLE PRECISION DCV
      DOUBLE PRECISION DTM
      LOGICAL LAST
      INTEGER M
      INTEGER MM
      INTEGER N
      DOUBLE PRECISION TM
      DOUBLE PRECISION TOLA

      COMMON / MF2 / A0, A, A1

      N = 0
      LAST = .FALSE.

50    N = N + 1

      CALL TMOFA ( A, TM, DTM, M )
      IF ( M .GT. 0 ) GO TO 400
      D = -TM / DTM
C
C  IS TOLERANCE MET?
C
      IF ( N .EQ. 40 .OR. 
     &  A - A0 .LE. TOLA .OR.
     &  A1 - A .LE. TOLA .OR.
     &  DABS ( D ) .LT. TOLA ) LAST = .TRUE.
      IF ( D ) 110, 100, 120

100   CV = A
      DCV = 0.0D+00
      GO TO 320
C
C  REPLACE UPPER BOUND BY A.
C
110   A1 = A
      GO TO 200
C
C  REPLACE LOWER BOUND BY A.
C
120   A0 = A

200   A2 = A + D

      IF ( LAST ) GO TO 300
      IF ( A2 .GT. A0 .AND. A2 .LT. A1 ) GO TO 250

      A = 0.5D+00 * ( A0 + A1 )
      GO TO 50

250   A = A2
      GO TO 50

300   IF ( A2 .LE. A0 .OR. A2 .GE. A1 ) GO TO 350

      CALL TMOFA ( A2, TM, DTM, M )

      IF ( M .GT. 0 ) GO TO 400

      D = - TM / DTM
      CV = A2

310   DCV = D

320   MM = M
      RETURN

350   CV = A
      GO TO 310

400   CV = 0.0D+00
      DCV = 0.0D+00
      GO TO 320
      END
      SUBROUTINE TMOFA ( ALFA, TM, DTM, ND )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A(3)
      DOUBLE PRECISION AA
      DOUBLE PRECISION ALFA
      DOUBLE PRECISION B(3)
      DOUBLE PRECISION DG(200,2)
      DOUBLE PRECISION DTM
      DOUBLE PRECISION DTYPE
      DOUBLE PRECISION F
      DOUBLE PRECISION FL
      DOUBLE PRECISION G(200,2)
      DOUBLE PRECISION H(200)
      DOUBLE PRECISION HP
      INTEGER K
      INTEGER KK
      INTEGER KT
      INTEGER L
      INTEGER MF
      INTEGER M0
      INTEGER M1
      INTEGER M2S
      INTEGER ND
      INTEGER TYPE
      DOUBLE PRECISION Q
      DOUBLE PRECISION Q1
      DOUBLE PRECISION Q2
      DOUBLE PRECISION QINV
      DOUBLE PRECISION T
      DOUBLE PRECISION TM
      DOUBLE PRECISION TOL
      DOUBLE PRECISION TT
      DOUBLE PRECISION V

      COMMON G, DG, AA, A, B, DTYPE, QINV, Q1, Q2, T, TT, K, L, KK, KT

      COMMON / MF1 / Q, TOL, TYPE, M1, M0, M2S, MF

      EQUIVALENCE ( H(1), G(1,1) )
      EQUIVALENCE ( Q1, HP )
      EQUIVALENCE ( Q2, F )

      DATA FL / 1.0D+30 /
C
C  STATEMENT FUNCTION.
C
      V ( K ) = ( AA - DBLE ( FLOAT ( K ) )**2 ) / Q

      ND = 0
      KT = 0
      AA = ALFA
      DTYPE = TYPE
      QINV = 1.0D+00 / Q
      DO 10 L = 1, 2
        DO 5 K = 1, 200
          G(K,L) = 0.0D+00
          DG(K,L) = 0.0D+00
5       CONTINUE
10    CONTINUE

      IF ( MOD ( TYPE, 2 ) ) 20, 30, 20

20    M0 = 3
      GO TO 40

30    M0 = TYPE + 2
40    K = 0.5D+00 * DSQRT ( DMAX1 ( 3.0D+00 * Q + AA, 0.0D+00 ) )
      M2S = MIN0 ( 2 * K + M0 + 4, 398 + MOD ( M0, 2 ) )
C
C  EVALUATION OF THE TAIL OF A CONTINUED FRACTION.
C
      A(1) = 1.0D+00
      A(2) = V(M2S+2)
      B(1) = V(M2S)
      B(2) = A(2) * B(1) - 1.0D+00
      Q1 = A(2) / B(2)

      DO 50 K = 1, 200
        MF = M2S + 2 + 2 * K
        T = V(MF)
        A(3) = T * A(2) - A(1)
        B(3) = T * B(2) - B(1)
        Q2 = A(3) / B(3)

        IF ( DABS ( Q1 - Q2 ) .LT. TOL ) GO TO 70

        Q1 = Q2
        A(1) = A(2)
        A(2) = A(3)
        B(1) = B(2)
        B(2) = B(3)
50    CONTINUE

      KT = 1

70    T = 1.0D+00 / T
      TT = - T * T * QINV
      L = MF - M2S

      DO 80 K = 2, L, 2
        T = 1.0D+00 / ( V(MF-K) - T )
        TT = T * T * ( TT - QINV )
80    CONTINUE

      KK = M2S / 2 + 1

      IF ( KT .EQ. 1 ) Q2 = T

      G(KK,2) = 0.5D+00 * ( Q2 + T )
      DG(KK,2) = TT
C
C  STAGE 1.
C
      G(2,1) = 1.0D+00

      DO 140 K = M0, M2S, 2

        KK = K / 2 + 1
        IF ( K .LT. 5 ) IF ( K - 3 ) 100, 110, 120
        G(KK,1) = V(K-2) - 1.0D+00 / G(KK-1,1)
        DG(KK,1) = QINV + DG(KK-1,1) / G(KK-1,1)**2
        GO TO 130

100     G(2,1) = V ( 0 )
        DG(2,1) = QINV
        GO TO 130

110     G(2,1) = V ( 1 ) + DTYPE - 2.0D+00
        DG(2,1) = QINV
        GO TO 130


120     G(3,1) = V ( 2 ) + ( DTYPE - 2.0D+00 ) / G(2,1)
        DG(3,1) = QINV + ( 2.0D+00 - DTYPE ) * DG(2,1) / G(2,1)**2

        IF ( TYPE .EQ. 2 ) G(2,1) = 0.0D+00

130     IF ( DABS ( G(KK,1) ) .LT. 1.0D+00 ) GO TO 200

140   CONTINUE
C
C  BACKTRACK.
C
      TM = G(KK,2) - G(KK,1)
      DTM = DG(KK,2) - DG(KK,1)
      M1 = M2S
      KT = M2S - M0

      DO 180 L = 2, KT, 2

        K = M2S - L
        KK = K / 2 + 1
        G(KK,2) = 1.0D+00 / ( V(K) - G(KK+1,2) )
        DG(KK,2) = -G(KK,2)**2 * ( QINV - DG(KK+1,2) )
        IF ( K - 2 ) 150, 150, 160

150     G(2,2) = 2.0D+00 * G(2,2)
        DG(2,2) = 2.0D+00 * DG(2,2)

160     T = G(KK,2) - G(KK,1)
        IF ( DABS ( T ) - DABS ( TM ) ) 170, 180, 180

170     TM = T
        DTM = DG(KK,2) - DG(KK,1)
        M1 = K

180   CONTINUE

      GO TO 320
C
C  STAGE 2.
C
200   CONTINUE

      M1 = K
      K = M2S
      KK = K / 2 + 1

210   CONTINUE

      IF ( K .EQ. M1 ) IF ( K - 2 ) 300, 300, 310
      K = K - 2
      KK = KK - 1
      T = V(K) - G(KK+1,2)
      IF ( DABS ( T ) - 1.0D+00 ) 250, 220, 220

220   CONTINUE

      G(KK,2) = 1.0D+00 / T
      DG(KK,2) = ( DG(KK+1,2) - QINV ) / T**2
      GO TO 210
C
C  STAGE 3.
C
250   CONTINUE

      IF ( K .EQ. M1 ) IF ( T ) 220, 290, 220
      HP = DG(KK+1,2) - QINV

260   CONTINUE

      G(KK,2) = FL
      H(KK) = T
      K = K - 2
      KK = KK - 1
      F = V(K) * T - 1.0D+00
      IF ( K .EQ. M1 ) IF ( F ) 280, 290, 280
      IF ( DABS ( F ) - DABS ( T ) ) 270, 280, 280

270   CONTINUE

      HP = HP / T**2 - QINV
      T = F / T
      GO TO 270

280   G(KK,2) = T / F
      DG(KK,2) = ( HP - QINV * T * T ) / F**2
      GO TO 210

290   ND = 1
      GO TO 320
C
C  CHAINING M EQUALS 3.
C
300   G(2,2) = 2.0D+00 * G(2,2)
      DG(2,2) = 2.0D+00 * DG(2,2)

310   TM = G(KK,2) - G(KK,1)
      DTM = DG(KK,2) - DG(KK,1)

320   RETURN
      END
      SUBROUTINE COEF ( M )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A
      DOUBLE PRECISION AB(200)
      REAL DUM1(800)
      DOUBLE PRECISION FL
      DOUBLE PRECISION G(200,2)
      DOUBLE PRECISION H(200)
      INTEGER K
      INTEGER KA
      INTEGER KB
      INTEGER KK
      INTEGER M
      INTEGER MF
      INTEGER ML
      INTEGER MM
      INTEGER M0
      INTEGER M1
      INTEGER M2S
      DOUBLE PRECISION Q
      DOUBLE PRECISION T
      DOUBLE PRECISION TOL
      INTEGER TYPE
      DOUBLE PRECISION V
      DOUBLE PRECISION V2

      SAVE FL
      SAVE V2

      COMMON G, DUM1, A, T, K, KA, KB, KK, MM, ML, AB

      COMMON / MF1 / Q, TOL, TYPE, M1, M0, M2S, MF

      EQUIVALENCE ( H(1), G(1,1) )

      DATA FL / 1.0D+30 /
      DATA V2 / 1.0D-15 /
C
C  STATEMENT FUNCTION:
C
      V ( K ) = ( A - DBLE ( FLOAT ( K ) )**2 ) / Q

      CALL TMOFA ( A, T, T, M )

      IF ( M .NE. 0 ) GO TO 300

      DO 60 K = 1, 200
        AB(K) = 0.0D+00
60    CONTINUE

      KA = M1 - M0 + 2

      DO 90 K = 2, KA, 2
        KK = ( M1 - K ) / 2 + 1
        IF ( K - 2 ) 70, 70, 80
70      AB(KK) = 1.0D+00
        GO TO 90
80      AB(KK) = AB(KK+1) / G(KK+1,1)
90    CONTINUE

      KA = 0

      DO 130 K = M1, M2S, 2

        KK = K / 2 + 1
        ML = K
        IF ( G(KK,2) .EQ. FL ) GO TO 100
        AB(KK) = AB(KK-1) * G(KK,2)
        GO TO 110

100     T = AB(KK-2)
        IF ( K .EQ. 4 .AND. M1 .EQ. 2 ) T = T + T
        AB(KK) = T / ( V(K-2) * H(KK) - 1.0D+00 )

110     IF ( DABS ( AB(KK) ) .GE. 1.0D-17 ) KA = 0
        IF ( KA .EQ. 5 ) GO TO 260
        KA = KA + 1

130   CONTINUE

      T = DLOG ( DABS ( AB(KK) ) / V2 ) 
     &  / DLOG ( 1.0D+00 / DABS ( G(KK,2) ) )
      KA = 2 * IDINT ( T )
      ML = KA + 2 + M2S
      IF ( ML .GT. 399 ) GO TO 400
      KB = KA + 2 + MF
      T = 1.0D+00 / V(KB)
      KK = MF - M2S

      DO 150 K = 2, KK, 2
        T = 1.0D+00 / ( V(KB-K) - T )
150   CONTINUE

      KK = ML / 2 + 1
      G(KK,2) = T

      DO 200 K = 2, KA, 2
        KK = ( ML - K ) / 2 + 1
        G(KK,2) = 1.0D+00 / ( V(ML-K) - G(KK+1,2) )
200   CONTINUE

      KA = M2S + 2
      DO 250 K = KA, ML, 2
        KK = K / 2 + 1
        AB(KK) = AB(K-1) * G(KK,2)
250   CONTINUE
C
C  NEUTRAL NORMALIZATION.
C
260   T = AB(1)
      MM = MOD ( TYPE, 2 )
      KA = MM + 2
      DO 280 K = KA, ML, 2
        KK = K / 2 + 1
        IF ( DABS ( T ) - DABS ( AB(KK) ) ) 270, 280, 280
270     T = AB(KK)
        MM = K
280   CONTINUE

      DO 290 K = 1, KK
        AB(K) = AB(K) / T
290   CONTINUE

300   RETURN
400   M = -1
      GO TO 300
      END
      SUBROUTINE SUM_SERIES ( DUM )

c*********************************************************************72
c
c  To avoid conflict with the FORTRAN90 intrinsic SUM function,
c  this routine, originally named "SUM", has been renamed "SUM_SERIES".
c
      IMPLICIT NONE

      DOUBLE PRECISION DLAST
      DOUBLE PRECISION DMAX
      DOUBLE PRECISION DUM
      REAL DUM1(1006)
      REAL DUM2(6)
      DOUBLE PRECISION F
      INTEGER K
      INTEGER KLAST
      INTEGER KMAX
      INTEGER L
      INTEGER S
      DOUBLE PRECISION T

      COMMON DUM1, S, L, DUM2, F, DMAX, DLAST, KMAX, KLAST, T

      K = 0
      F = DUM ( 0 )
      DMAX = F
      T = DABS ( F )
      KMAX = 0
      DO 30 KLAST = 1, L
        DLAST = DUM ( KLAST )
        F = F + DLAST
        IF ( T - DABS ( DLAST ) ) 10, 10, 20
10      DMAX = DLAST
        T = DABS ( DMAX )
        KMAX = KLAST
20      IF ( KLAST .LE. S ) GO TO 30
        IF ( DABS ( DLAST ) / T .GT. 1.0D-13 ) K = 0
        K = K + 1
        IF ( K .EQ. 3 ) GO TO 40
30    CONTINUE

      KLAST = L

40    RETURN
      END
      SUBROUTINE BESSEL ( SOL, U, JY, N )

c*********************************************************************72
C
      IMPLICIT NONE

      INTEGER N

      DOUBLE PRECISION JY(N)
      INTEGER K
      INTEGER NN
      INTEGER SOL
      DOUBLE PRECISION U

      NN = MIN0 ( N, 249 )

      IF ( U .EQ. 0.0D+00 .AND. SOL .EQ. 2 ) GO TO 80
      IF ( U .GE. 8.0D+00 ) GO TO 30
      GO TO ( 10, 20 ), SOL

10    CONTINUE

      CALL J0J1 ( U, JY )
      GO TO 40

20    CONTINUE

      CALL Y0Y1 ( U, JY )
      GO TO 40

30    CONTINUE

      CALL LUKE ( U, SOL, JY )

40    CONTINUE

      IF ( N .LT. 2 ) GO TO 100
      GO TO ( 50, 60 ), SOL

50    CONTINUE

      CALL JNS ( JY, U, NN )
      GO TO 100
C
C  RECURRENCE FORMULA.
C
60    CONTINUE

      DO 70 K = 2, NN
        JY(K+1) = 2.0D+00 
     &    * DBLE ( FLOAT ( K - 1 ) ) * JY(K) / U - JY(K-1)
70    CONTINUE

      GO TO 100

80    CONTINUE

      NN = NN + 1

      DO 90 K = 1, NN
        JY(K) = -1.0D+37
90    CONTINUE

100   CONTINUE

      RETURN
      END
      SUBROUTINE J0J1 ( X, J )

c*********************************************************************72

      IMPLICIT NONE

      REAL DUM(1014)
      DOUBLE PRECISION J(2)
      DOUBLE PRECISION T(5)
      DOUBLE PRECISION X

      COMMON DUM, T

      T(1) = X / 2.0D+00
      J(1) = 1.0D+00
      J(2) = T(1)
      T(2) = -T(1)**2
      T(3) = 1.0D+00
      T(4) = 1.0D+00

10    T(4) = T(4) * T(2) / T(3)**2
      J(1) = J(1) + T(4)
      T(5) = T(4) * T(1) / ( T(3) + 1.0D+00 )
      J(2) = J(2) + T(5)

      IF ( DMAX1 ( DABS ( T(4) ), DABS ( T(5) ) ) .LT. 1.0D-15 ) RETURN

      T(3) = T(3) + 1.0D+00
      GO TO 10
      END
      SUBROUTINE Y0Y1 ( X, Y )

      IMPLICIT NONE

      REAL DUM(1014)
      DOUBLE PRECISION T(10)
      DOUBLE PRECISION X
      DOUBLE PRECISION Y(2)

      COMMON DUM, T

      T(1) = X / 2.0D+00
      T(2) = -T(1)**2
      Y(1) = 1.0D+00
      Y(2) = T(1)
      T(7) = 0.0D+00
      T(10) = -T(1)
      T(3) = 0.0D+00
      T(4) = 0.0D+00
      T(5) = 1.0D+00

10    CONTINUE

      T(3) = T(3) + 1.0D+00
      T(4) = T(4) + 1.0D+00 / T(3)
      T(5) = T(5) * T(2) / T(3)**2
      Y(1) = Y(1) + T(5)
      T(6) = -T(5) * T(4)
      T(7) = T(7) + T(6)
      T(8) = T(5) * T(1) / ( T(3) + 1.0D+00 )
      Y(2) = Y(2) + T(8)
      T(9) = -T(8) * ( 2.0D+00 * T(4) + 1.0D+00 / ( T(3) + 1.0D+00 ) )
      T(10) = T(10) + T(9)

      IF ( DMAX1 ( DABS ( T(6) ), DABS ( T(9) ) ) .GE. 1.0D-15 ) 
     &  GO TO 10

      T(2) = 0.57721566490153286D+00 + DLOG ( T(1) )
      Y(1) = 0.63661977236758134D+00 * ( Y(1) * T(2) + T(7) )
      Y(2) = 0.63661977236758134D+00 * ( Y(2) * T(2) - 1.0D+00 / X ) 
     &  + T(10) / 3.1415926535897932D+00

      RETURN
      END
      SUBROUTINE LUKE ( U, KIND, JY )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A(19)
      DOUBLE PRECISION B(19)
      DOUBLE PRECISION C(19)
      DOUBLE PRECISION CS
      DOUBLE PRECISION D(19)
      REAL DUM(1014)
      DOUBLE PRECISION G(3)
      DOUBLE PRECISION JY(2)
      INTEGER K
      INTEGER KIND
      DOUBLE PRECISION R(2)
      DOUBLE PRECISION S(2)
      DOUBLE PRECISION SN
      DOUBLE PRECISION T
      DOUBLE PRECISION U
      DOUBLE PRECISION X

      SAVE A
      SAVE B
      SAVE C
      SAVE D

      COMMON DUM, R, S, G, X, T, SN, CS

      DATA A /
     &   0.99959506476867287416D+00,
     &  -0.53807956139606913D-03,
     &  -0.13179677123361570D-03,
     &   0.151422497048644D-05,
     &   0.15846861792063D-06,
     &  -0.856069553946D-08,
     &  -0.29572343355D-09,
     &   0.6573556254D-10,
     &  -0.223749703D-11,
     &  -0.44821140D-12,
     &   0.6954827D-13,
     &  -0.151340D-14,
     &  -0.92422D-15,
     &   0.15558D-15,
     &  -0.476D-17,
     &  -0.274D-17,
     &   0.61D-18,
     &  -0.4D-19,
     &  -0.1D-19 /
      DATA B /
     &  -0.776935569420532136D-02,
     &  -0.774803230965447670D-02,
     &   0.2536541165430796D-04,
     &   0.394273598399711D-05,
     &  -0.10723498299129D-06,
     &  -0.721389799328D-08,
     &   0.73764602893D-09,
     &   0.150687811D-11,
     &  -0.574589537D-11,
     &   0.45996574D-12,
     &   0.2270323D-13,
     &  -0.887890D-14,
     &   0.74497D-15,
     &   0.5847D-16,
     &  -0.2410D-16,
     &   0.265D-17,
     &   0.13D-18,
     &  -0.10D-18,
     &   0.2D-19 /

      DATA C /
     &   1.00067753586591346234D+00,
     &   0.90100725195908183D-03,
     &   0.22172434918599454D-03,
     &  -0.196575946319104D-05,
     &  -0.20889531143270D-06,
     &   0.1028144350894D-07,
     &   0.37597054789D-10,
     &  -0.7638891358D-10,
     &   0.238734670D-11,
     &   0.51825489D-12,
     &  -0.7693969D-13,
     &   0.144008D-14,
     &   0.103294D-14,
     &  -0.16821D-15,
     &   0.459D-17,
     &   0.302D-17,
     &  -0.65D-18,
     &   0.4D-19,
     &   0.1D-19 /
      DATA D /
     &   0.2337682998628580328D-01,
     &   0.2334680122354557533D-01,
     &  -0.3576010590901382D-04,
     &  -0.560863149492627D-05,
     &   0.13273894084340D-06,
     &   0.916975845066D-08,
     &  -0.86838880371D-09,
     &  -0.378073005D-11,
     &   0.663145586D-11,
     &  -0.50584390D-12,
     &  -0.2720782D-13,
     &   0.985381D-14,
     &  -0.79398D-15,
     &  -0.6757D-16,
     &   0.2625D-16,
     &  -0.280D-17,
     &  -0.15D-18,
     &   0.10D-18,
     &  -0.2D-19 /

      X = 8.0D+00 / U
      G(1) = 1.0D+00
      G(2) = 2.0D+00 * X - 1.0D+00
      R(1) = A(1) + A(2) * G(2)
      S(1) = B(1) + B(2) * G(2)
      R(2) = C(1) + C(2) * G(2)
      S(2) = D(1) + D(2) * G(2)
      DO 10 K = 3, 19
        G(3) = ( 4.0D+00 * X - 2.0D+00 ) * G(2) - G(1)
        R(1) = R(1) + A(K) * G(3)
        S(1) = S(1) + B(K) * G(3)
        R(2) = R(2) + C(K) * G(3)
        S(2) = S(2) + D(K) * G(3)
        G(1) = G(2)
        G(2) = G(3)
10    CONTINUE

      T = 0.7978845608028654D+00 / DSQRT ( U )
      SN = DSIN ( U - 0.7853981633974483D+00 )
      CS = DCOS ( U - 0.7853981633974483D+00 )

      GO TO ( 20, 30 ), KIND

20    JY(1) = T * ( R(1) * CS - S(1) * SN )
      JY(2) = T * ( R(2) * SN + S(2) * CS )
      GO TO 40

30    JY(1) = T * ( S(1) * CS + R(1) * SN )
      JY(2) = T * ( S(2) * SN - R(2) * CS )

40    RETURN
      END
      SUBROUTINE JNS ( JJ, U, M )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION D(2)
      DOUBLE PRECISION DM
      REAL DUM(1014)
      DOUBLE PRECISION G(249)
      DOUBLE PRECISION JJ(250)
      INTEGER K
      INTEGER KA
      INTEGER KK
      INTEGER M
      INTEGER M2
      DOUBLE PRECISION P(3)
      DOUBLE PRECISION Q(3)
      DOUBLE PRECISION U

      EQUIVALENCE ( A, G(1) )
      EQUIVALENCE ( D(1), G(2) )
      EQUIVALENCE ( P(1), G(4) )
      EQUIVALENCE ( Q(1), G(7) )
      EQUIVALENCE ( DM, G(10) )
      EQUIVALENCE ( B, G(11) )

      COMMON DUM, G, M2, K, KK, KA

      M2 = M

      DM = 2.0D+00 * M
      P(1) = 0.0D+00
      Q(1) = 1.0D+00
      P(2) = 1.0D+00
      Q(2) = DM / U
      D(1) = 1.0D+00 / Q(2)
      A = 2.0D+00

10    B = ( DM + A ) / U
      P(3) = B * P(2) - P(1)
      Q(3) = B * Q(2) - Q(1)
      D(2) = P(3) / Q(3)

      IF ( DABS ( D(1) - D(2) ) .LT. 1.0D-15 ) GO TO 20

      P(1) = P(2)
      P(2) = P(3)
      Q(1) = Q(2)
      Q(2) = Q(3)
      D(1) = D(2)
      A = A + 2.0D+00
      GO TO 10

20    G(M) = D(2)
      KA = M - 2

      DO 30 K = 1, KA
        KK = M - K
        A = 2 * KK
        G(KK) = U / ( A - U * G(KK+1) )
        IF ( G(KK) .EQ. 0.0D+00 ) G(KK) = 1.0D-35
30    CONTINUE

      DO 40 K = 2, M
        JJ(K+1) = G(K) * JJ(K)
40    CONTINUE

      RETURN
      END
      DOUBLE PRECISION FUNCTION DS ( KK )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE RADIAL SOLUTION ASSOCIATED WITH B(Q).
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1004)
      REAL DUM2(17)
      REAL DUM3(583)
      DOUBLE PRECISION FJ
      DOUBLE PRECISION FY
      INTEGER K
      INTEGER KK
      INTEGER N
      INTEGER N1
      INTEGER N2
      INTEGER P
      INTEGER S

      COMMON DUM1, N, P, S, DUM2,
     &  K, N1, N2, DUM3, AB

      K = KK
      N1 = K - S
      N2 = K + S + P
      DS = AB(K+1) * ( FJ ( N1 ) * FY ( N2 ) - FJ ( N2 ) * FY ( N1 ) )
      IF ( MOD ( K + N, 2 ) .NE. 0 ) DS = -DS

      RETURN
      END
      DOUBLE PRECISION FUNCTION DC ( KK )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE RADIAL SOLUTION ASSOCIATED WITH A(Q).
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1004)
      REAL DUM2(17)
      REAL DUM3(583)
      DOUBLE PRECISION FJ
      DOUBLE PRECISION FY
      INTEGER K
      INTEGER KK
      INTEGER N
      INTEGER N1
      INTEGER N2
      INTEGER P
      INTEGER S     

      COMMON DUM1, N, P, S, DUM2,
     &  K, N1, N2, DUM3, AB

      K = KK
      N1 = K - S
      N2 = K + S + P     
      DC = AB(K+1) * ( FJ ( N1 ) * FY ( N2 ) - FJ ( N2 ) * FY ( N1 ) )
      IF ( MOD ( K + N, 2 ) .NE. 0 ) DC = -DC
      IF ( S + P .EQ. 0 ) DC = 0.5D+00 * DC

      RETURN
      END
      DOUBLE PRECISION FUNCTION DDS ( KK )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE DERIVATIVE OF THE RADIAL SOLUTION
C  ASSOCIATED WITH B(Q).
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1000)
      REAL DUM2(17)
      REAL DUM3(583)
      DOUBLE PRECISION DJ
      DOUBLE PRECISION DY
      DOUBLE PRECISION FJ
      DOUBLE PRECISION FY
      INTEGER K
      INTEGER KK
      INTEGER N
      INTEGER N1
      INTEGER N2
      INTEGER P
      INTEGER S
      DOUBLE PRECISION U1
      DOUBLE PRECISION U2

      COMMON DUM1, U1, U2, N, P, S, DUM2,
     &  K, N1, N2, DUM3, AB

      K = KK
      N1 = K - S
      N2 = K + S + P
      DDS = AB(K+1) * ( 
     &    U2 * ( FJ ( N1 ) * DY ( N2 ) 
     &         - FJ ( N2 ) * DY ( N1 ) )
     &  - U1 * ( FY ( N2 ) * DJ ( N1 ) 
     &         - FY ( N1 ) * DJ ( N2 ) ) )
 
      IF ( MOD ( K + N, 2 ) .NE. 0 ) DDS = -DDS

      RETURN
      END
      DOUBLE PRECISION FUNCTION DDC ( KK )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE DERIVATIVE OF THE RADIAL SOLUTION 
C  ASSOCIATED WITH A(Q).
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1000)
      REAL DUM2(17)
      REAL DUM3(583)
      DOUBLE PRECISION DJ
      DOUBLE PRECISION DY
      DOUBLE PRECISION FJ
      DOUBLE PRECISION FY
      INTEGER K
      INTEGER KK
      INTEGER N
      INTEGER N1
      INTEGER N2
      INTEGER P
      INTEGER S
      DOUBLE PRECISION U1
      DOUBLE PRECISION U2

      COMMON DUM1, U1, U2, N, P, S, DUM2,
     &  K, N1, N2, DUM3, AB

      K = KK
      N1 = K - S
      N2 = K + S + P     
      DDC = AB(K+1) * ( 
     &    U2 * ( FJ ( N1 ) * DY ( N2 ) 
     &         + FJ ( N2 ) * DY ( N1 ) )
     &  - U1 * ( FY ( N2 ) * DJ ( N1 ) 
     &         + FY ( N1 ) * DJ ( N2 ) ) )

      IF ( MOD ( K + N, 2 ) .NE. 0 ) DDC = -DDC
      IF ( S + P .EQ. 0 ) DDC = 0.5D+00 * DDC

      RETURN
      END
      DOUBLE PRECISION FUNCTION PS ( K )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE ODD PERIODIC SOLUTION.
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1005)
      REAL DUM2(2)
      REAL DUM3(600)
      INTEGER K
      INTEGER P
      DOUBLE PRECISION X

      COMMON DUM1, P, DUM2, X, DUM3, AB

      PS = AB(K+1) * DSIN ( DBLE ( FLOAT ( 2 * K + P ) ) * X )

      RETURN
      END
      DOUBLE PRECISION FUNCTION PC ( K )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE EVEN PERIODIC SOLUTION.
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1005)
      REAL DUM2(2)
      REAL DUM3(600)
      INTEGER K
      INTEGER P
      DOUBLE PRECISION X

      COMMON DUM1, P, DUM2, X, DUM3, AB

      PC = AB(K+1) * DCOS ( DBLE ( FLOAT ( 2 * K + P ) ) * X )

      RETURN
      END
      DOUBLE PRECISION FUNCTION DPS ( K )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE DERIVATIVE OF THE ODD PERIODIC SOLUTION.
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1005)
      REAL DUM2(2)
      REAL DUM3(14)
      REAL DUM4(584)
      INTEGER K
      INTEGER P
      DOUBLE PRECISION T
      DOUBLE PRECISION X

      COMMON DUM1, P, DUM2, X, DUM3, T, DUM4, AB

      T = 2 * K + P
      DPS = AB(K+1) * T * DCOS ( T * X )

      RETURN
      END
      DOUBLE PRECISION FUNCTION DPC ( K )

c*********************************************************************72
C
C  EVALUATES ONE TERM OF THE DERIVATIVE OF THE EVEN PERIODIC SOLUTION.
C
      IMPLICIT NONE

      DOUBLE PRECISION AB(200)
      REAL DUM1(1005)
      REAL DUM2(2)
      REAL DUM3(14)
      REAL DUM4(584)
      INTEGER K
      INTEGER P
      DOUBLE PRECISION T
      DOUBLE PRECISION X

      COMMON DUM1, P, DUM2, X, DUM3, T, DUM4, AB

      T = 2 * K + P
      DPC = -AB(K+1) * T * DSIN ( T * X )

      RETURN
      END
      DOUBLE PRECISION FUNCTION FJ ( N )

c*********************************************************************72
C
CC FJ produces Bessel functions of the first kind.
C
      IMPLICIT NONE

      REAL DUM(527)
      DOUBLE PRECISION J(250)
      INTEGER K
      INTEGER N

      COMMON J, DUM, K

      K = IABS ( N )

      IF ( K .GE. 250 ) GO TO 20
      FJ = J(K+1)
      IF ( MOD ( N, 2 ) .LT. 0 ) FJ = -FJ
10    RETURN

20    FJ = 0.0D+00
      WRITE ( 6, 99 ) N
99    FORMAT ( 2H0J, I3, 7H NEEDED )
      GO TO 10
      END
      DOUBLE PRECISION FUNCTION FY ( N )

c*********************************************************************72
C
CC FY produces Bessel functions of the second kind.
C
      IMPLICIT NONE

      REAL DUM1(500)
      REAL DUM2(27)
      INTEGER K
      INTEGER N
      DOUBLE PRECISION Y(250)

      COMMON DUM1, Y, DUM2, K

      K = IABS ( N )

      IF ( K .GE. 250 ) GO TO 20
      FY = Y(K+1)
      IF ( MOD ( N, 2 ) .LT. 0 ) FY = -FY
10    RETURN

20    FY = 0.0D+00
      WRITE ( 6, 99 ) N
99    FORMAT ( 2H0Y, I3, 7H NEEDED )
      GO TO 10
      END
      DOUBLE PRECISION FUNCTION DJ ( N )

c*********************************************************************72
C
CC DJ computes derivatives of Bessel functions of the first kind.
C
      IMPLICIT NONE

      REAL DUM1(1000)
      REAL DUM2(26)
      DOUBLE PRECISION FJ
      DOUBLE PRECISION FN
      INTEGER N
      DOUBLE PRECISION U1

      COMMON DUM1, U1, DUM2, FN

      FN = N

      IF ( N - 249 ) 10, 20, 40

10    DJ = FN * FJ ( N ) / U1 - FJ ( N+1 )
      GO TO 30

20    DJ = FJ ( N-1 ) - FN * FJ ( N ) / U1

30    RETURN

40    DJ = 0.0D+00
      WRITE ( 6, 99 ) N
99    FORMAT ( 3H0J@, I3, 7H NEEDED )
      GO TO 30
      END
      DOUBLE PRECISION FUNCTION DY ( N )

c*********************************************************************72
C
C  DERIVATIVES OF BESSEL FUNCTIONS OF THE SECOND KIND.
C
      IMPLICIT NONE

      REAL DUM1(1002)
      REAL DUM2(24)
      DOUBLE PRECISION FN
      DOUBLE PRECISION FY
      INTEGER N
      DOUBLE PRECISION U2

      COMMON DUM1, U2, DUM2, FN

      IF ( N .GE. 250 ) GO TO 20

      FN = N
      DY = FY ( N-1 ) * FN * FY ( N ) / U2
10    RETURN

20    DY = 0.0D+00
      WRITE ( 6, 99 ) N
99    FORMAT ( 3H0Y@, I3, 7H NEEDED )
      GO TO 10
      END
