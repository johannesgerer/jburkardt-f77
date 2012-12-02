      SUBROUTINE SYMQR ( A, D, E, K0, N, NA, EPS, ABSCNV, VEC, TRD, 
     &  FAIL )

c*********************************************************************72
C
C  EXPLANATION OF THE PARAMETERS IN THE CALLING SEQUENCE:
C
C    A      A DOUBLE DIMENSIONED ARRAY.  IF THE MATRIX IS NOT
C           INITIALLY TRIDIAGONAL, IT IS CONTAINED IN THE LOWER
C           TRIANGLE OF A.  IF EIGENVECTORS ARE NOT REQUESTED,
C           THE LOWER TRIANGLE OF A IS DESTROYED WHILE THE
C           ELEMENTS ABOVE THE DIAGONAL ARE LEFT UNDISTURBED.
C           IF EIGENVECTORS ARE REQUESTED, THEY ARE RETURNED IN THE
C           COLUMNS OF A.
C
C    D      A SINGLY SUBSCRIPTED ARRAY.  IF THE MATRIX IS
C           INITIALLY TRIDIAGONAL, D CONTAINS ITS DIAGONAL 
C           ELEMENTS.  ON RETURN, D CONTAINS THE EIGENVALUES OF
C           THE MATRIX.
C
C    E      A SINGLY SUBSCRIPTED ARRAY.  IF THE MATRIX IS
C           INITIALLY TRIDIAGONAL, E CONTAINS ITS OFF-DIAGONAL
C           ELEMENTS.  UPON RETURN, E(I) CONTAINS THE NUMBER OF
C           ITERATIONS REQUIRED TO COMPUTE THE APPROXIMATE
C           EIGENVALUE D(I).
C
C    K0     A REAL VARIABLE CONTAINING AN INITIAL ORIGIN SHIFT TO
C           BE USED UNTIL THE COMPUTED SHIFTS SETTLE DOWN.
C
C    N      AN INTEGER VARIABLE CONTAINING THE FIRST DIMENSION
C           OF THE ARRAY A.
C
C    EPS    A REAL VARIABLE CONTAINING A CONVERGENCE TOLERANCE.
C
C    ABSCNV A LOGICAL VARIABLE CONTAINING THE VALUE .TRUE. IF
C           THE ABSOLUTE CONVERGENCE CRITERION IS TO BE USED
C           OR THE VALUE .FALSE. IF THE RELATIVE CRITERION
C           IS TO BE USED.
C
C    VEC    A LOGICAL VARIABLE CONTAINING THE VALUE .TRUE. IF
C           EIGENVECTORS ARE TO BE COMPUTED AND RETURNED IN
C           THE ARRAY A AND OTHERWISE CONTAINING THE VALUE
C           .FALSE..
C
C    TRD    A LOGICAL VARIABLE CONTAINING THE VALUE .TRUE.
C           IF THE MATRIX IS TRIDIAGONAL AND LOCATED IN THE ARRAYS
C           D AND E AND OTHERWISE CONTAINING THE VALUE .FALSE..
C
C    FAIL   AN INTEGER VARIABLE CONTAINING AN ERROR SIGNAL.
C           ON RETURN, THE EIGENVALUES IN D(FAIL+1), ..., D(N)
C           AND THEIR CORRESPONDING EIGENVECTORS MAY BE PRESUMED
C           ACCURATE.
C
      IMPLICIT NONE

      INTEGER N
      INTEGER NA

      REAL A(NA,N)
      LOGICAL ABSCNV
      REAL C
      REAL CB
      REAL CC
      REAL CD
      REAL CON
      REAL D(N)
      REAL E(N)
      REAL EPS
      INTEGER FAIL
      INTEGER I
      INTEGER I1
      INTEGER J
      REAL K
      REAL K0
      REAL K1
      REAL K2
      INTEGER L
      INTEGER L1
      INTEGER LL
      INTEGER LL1
      REAL MAX
      REAL NINF
      INTEGER NL
      INTEGER NM1
      INTEGER NM2
      INTEGER NNL
      REAL NORM
      INTEGER NU
      INTEGER NUM1
      REAL P
      REAL PP
      REAL Q
      REAL QQ
      REAL R
      INTEGER RETURN
      REAL S
      REAL S2
      LOGICAL SHFT
      REAL SUM
      REAL SUM1
      REAL TEMP
      REAL TEST
      REAL TITTER
      LOGICAL TRD
      LOGICAL VEC

      TITTER = 50.0E+00
      NM1 = N - 1
      NM2 = N - 2
      NINF = 0.0E+00
      FAIL = 0
C
C  SIGNAL ERROR IF N IS NOT POSITIVE.
C
      IF ( N .GT. 0 ) GO TO 1
      FAIL = -1
      RETURN
C
C  SPECIAL TREATMENT FOR A MATRIX OF ORDER ONE.
C
1     IF ( N .GT. 1 ) GO TO 5
      IF ( .NOT. TRD ) D(1) = A(1,1)
      IF ( VEC ) A(1,1) = 1.0E+00
      FAIL = 0
      RETURN
C
C  IF THE MATRIX IS TRIDIAGONAL, SKIP THE REDUCTION.C
C
5     IF ( TRD ) GO TO 100
      IF ( N .EQ. 2 ) GO TO 80
C
C  REDUCE THE MATRIX TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD.
C
      DO 70 L = 1, NM2

        L1 = L + 1
        D(L) = A(L,L)
        MAX = 0.0E+00

        DO 10 I = L1, N
10        MAX = AMAX1 ( MAX, ABS ( A(I,L) ) )

        IF ( MAX .NE. 0.0 ) GO TO 13
        E(L) = 0.0E+00
        A(L,L) = 1.0E+00
        GO TO 70
13      SUM = 0.0E+00

        DO 17 I = L1, N
          A(I,L) = A(I,L) / MAX
17        SUM = SUM + A(I,L)**2

        S2 = SUM
        S2 = SQRT ( S2 )
        IF ( A(L1,L) .LT. 0.0 ) S2 = -S2
        E(L) = -S2 * MAX
        A(L1,L) = A(L1,L) + S2
        A(L,L) = S2 * A(L1,L)

        SUM1 = 0.0E+00
        DO 50 I = L1, N
          SUM = 0.0E+00
          DO 20 J = L1, I
20          SUM = SUM + A(I,J) * A(J,L)
          IF ( I .EQ. N ) GO TO 40
          I1 = I + 1
          DO 30 J = I1, N
30          SUM = SUM + A(J,L) * A(J,I)
40        E(I) = SUM / A(L,L)
50        SUM1 = SUM1 + A(I,L) * E(I)

        CON = 0.5E+00 * SUM1 / A(L,L)
        DO 60 I = L1, N
          E(I) = E(I) - CON * A(I,L)
          DO 60 J = L1, I
60          A(I,J) = A(I,J) - A(I,L) * E(J) - A(J,L) * E(I)

70    CONTINUE

80    D(NM1) = A(NM1,NM1)
      D(N) = A(N,N)
      E(NM1) = A(N,NM1)
C
C  IF EIGENVECTORS ARE REQUIRED, INITIALIZE A.
C
100   IF ( .NOT. VEC ) GO TO 180
C
C  IF THE MATRIX WAS TRIDIAGONAL, SET A EQUAL TO THE IDENTITY MATRIX.
C
      IF ( .NOT. TRD .AND. N .NE. 2 ) GO TO 130

      DO 120 I = 1, N
        DO 110 J = 1, N
110       A(I,J) = 0.0E+00
120     A(I,I) = 1.0E+00

      GO TO 180
C
C  IF THE MATRIX WAS NOT TRIDIAGONAL, MULTIPLY OUT THE
C  TRANSFORMATIONS OBTAINED IN THE HOUSEHOLDER REDUCTION.C
C
130   A(N,N) = 1.0E+00
      A(NM1,NM1) = 1.0E+00
      A(NM1,N) = 0.0E+00

      DO 170 L = 1, NM2
        LL = NM2 - L + 1
        LL1 = LL + 1

        DO 140 I = LL1, N
          SUM = 0.0E+00
          DO 135 J = LL1, N
135         SUM = SUM + A(J,LL) * A(J,I)
140       A(LL,I) = SUM / A(LL,LL)

          DO 150 I = LL1, N
            DO 150 J = LL1, N
150           A(I,J) = A(I,J) - A(I,LL) * A(LL,J)

          DO 160 I = LL1, N
            A(I,LL) = 0.0E+00
160         A(LL,I) = 0.0E+00

170       A(LL,LL) = 1.0E+00
C
C  IF AN ABSOLUTE CONVERGENCE CRITERION IS REQUESTED,
C  ( ABSCNV = .TRUE. ), COMPUTE THE INFINITY NORM OF THE MATRIX.
C
180   IF ( .NOT. ABSCNV ) GO TO 200
      NINF = AMAX1 ( 
     &  NINF, ABS ( D(1) ) + ABS ( E(1) ) + ABS ( E(I-1) ) )
      IF ( N .EQ. 2 ) GO TO 200

      DO 190 I = 2, NM1
190     NINF = AMAX1 ( NINF, 
     &    ABS ( D(I) ) + ABS ( E(I) ) + ABS ( E(I-1) ) )
C
C  START THE QR ITERATION.
C
200   NU = N
      NUM1 = N - 1
      SHFT = .FALSE.
      K1 = K0
      TEST = NINF * EPS
      E(N) = 0.0E+00
C
C  CHECK FOR CONVERGENCE AND LOCATE THE SUBMATRIX IN WHICH THE
C  QR STEP IS TO BE PERFORMED.
C
210   DO 220 NNL = 1, NUM1
        NL = NUM1 - NNL + 1
        IF ( .NOT. ABSCNV ) 
     &    TEST = EPS * AMIN1 ( ABS ( D(NL) ), ABS ( D(NL+1) ) )
        IF ( ABS ( E(NL) ) .LE. TEST ) GO TO 230
220   CONTINUE

      GO TO 240
230   E(NL) = 0.0E+00
      NL = NL + 1
      IF ( NL .NE. NU ) GO TO 240
      IF ( NUM1 .EQ. 1 ) RETURN
      NU = NUM1
      NUM1 = NU - 1
      GO TO 210

240   E(NU) = E(NU) + 1.0E+00
      IF ( E(NU) .LE. TITTER ) GO TO 250
      FAIL = NU
      RETURN
C
C  CALCULATE THE SHIFT.
C
250   CB = ( D(NUM1) - D(NU) ) / 2.0E+00
      MAX = AMAX1 ( ABS ( CB ), ABS ( E(NUM1) ) )
      CB = CB / MAX
      CC = ( E(NUM1) / MAX )**2
      CD = SQRT ( CB**2 + CC )
      IF ( CB .NE. 0.0E+00 ) CD = SIGN ( CD, CB )
      K2 = D(NU) - MAX * CC / ( CB + CD )
      IF ( SHFT ) GO TO 270
      IF ( ABS ( K2 - K1 ) .LT. 0.5E+00 * ABS ( K2 ) ) GO TO 260
      K1 = K2
      K = K0
      GO TO 300
260   SHFT = .TRUE.
270   K = K2
C
C  PERFORM ONE QR STEP WITH SHIFT K ON ROWS AND COLUMNS
C  NL THROUGH NU.
C
300   P = D(NL) - K
      Q = E(NL)
      CALL SINCOS ( P, Q, C, S, NORM )

310   DO 380 I = NL, NUM1
C
C  IF REQUIRED, ROTATE THE EIGENVECTORS.
C
        IF ( .NOT. VEC ) GO TO 330

        DO 320 J = 1, N
          TEMP = C * A(J,I) + S * A(J,I+1)
          A(J,I+1) = -S * A(J,I) + C*A(J,I+1)
320       A(J,I) = TEMP
C
C  PERFORM THE SIMILARITY TRANSFORMATION AND CALCULATE THE NEXT
C  ROTATION.
C
330     D(I) = C * D(I) + S * E(I)
        TEMP = C * E(I) + S * D(I+1)
        D(I+1) = -S * E(I) + C * D(I+1)
        E(I) = -S * K
        D(I) = C * D(I) + S * TEMP
        IF ( I .EQ. NUM1 ) GO TO 380
        IF ( ABS ( S ) .GT. ABS ( C ) ) GO TO 350
        R = S / C
        D(I+1) = -S * E(I) + C * D(I+1)
        P = D(I+1) - K
        Q = C * E(I+1)
        CALL SINCOS ( P, Q, C, S, NORM )
340     E(I) = R * NORM
        E(I+1) = Q
        GO TO 380
350     P = C * E(I) + S * D(I+1)
        Q = S * E(I+1)
        D(I+1) = C * P / S + K
        E(I+1) = C * E(I+1)
        CALL SINCOS ( P, Q, C, S, NORM )
360     E(I) = NORM

380   CONTINUE

      TEMP = C * E(NUM1) + S * D(NU)
      D(NU) = -S * E(NUM1) + C * D(NU)
      E(NUM1) = TEMP
      GO TO 210

      END
      SUBROUTINE SINCOS ( P, Q, C, S, NORM )

c*********************************************************************72
C
C  SINCOS CALCULATES THE ROTATION CORRESPONDING TO THE VECTOR (P,Q).
C
      IMPLICIT NONE

      REAL C
      REAL NORM
      REAL P
      REAL PP
      REAL Q
      REAL QQ
      REAL S

500   PP = ABS ( P )
      QQ = ABS ( Q )
      IF ( QQ .GT. PP ) GO TO 510
      NORM = PP * SQRT ( 1.0E+00 + ( QQ / PP )**2 )
      GO TO 520
510   IF ( QQ .EQ. 0.0E+00 ) GO TO 530
      NORM = QQ * SQRT ( 1.0E+00 + ( PP / QQ )**2 )
520   C = P / NORM
      S = Q / NORM
      RETURN
530   C = 1.0E+00
      S = 0.0E+00
      NORM = 0.0E+00
      RETURN
      END
