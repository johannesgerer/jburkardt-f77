      SUBROUTINE ROMINT ( VAL, ERR, EPS, A, B, N, MAXE, F )

c*********************************************************************72
c
cc ROMINT
c
      IMPLICIT NONE

      REAL A
      REAL B
      REAL BB
      REAL EPS
      REAL ERR
      EXTERNAL F
      REAL F
      REAL H
      INTEGER J
      INTEGER K
      INTEGER K0
      INTEGER K1
      INTEGER K2
      INTEGER KK
      INTEGER KKK
      INTEGER L
      INTEGER MAXE
      INTEGER N
      INTEGER N0
      INTEGER N1
      REAL R
      REAL RM(16)
      REAL S
      REAL S0
      REAL S1
      REAL T
      REAL VAL
C
C  INITIAL TRAPEZOID RULE.
C
      T = ( B - A ) * ( F(A) + F(B) ) * 0.5E+00
C
C  INITIAL RECTANGLE VALUE.
C
      RM(1) = ( B - A ) * F ( ( A + B ) * 0.5E+00 )

      N = 2
      R = 4.0E+00

      DO 11 K = 1, MAXE

        BB = ( R * 0.5E+00 - 1.0E+00 ) / ( R - 1.0E+00 )
C
C  IMPROVED TRAPEZOID VALUE.
C
        T = RM(1) + BB * ( T - RM(1) )
C
C  DOUBLE NUMBER OF SUBDIVISIONS OF (A,B).
C
        N = 2 * N
        S = 0.0E+00
        H = ( B - A ) / FLOAT ( N )
C
C  CALCULATE RECTANGLE VALUE.
C
        IF ( N - 32 ) 1, 1, 2
1       N0 = N
        GO TO 4
2       N0 = 32
3       IF ( N - 512 ) 4, 4, 5
4       N1 = N
        GO TO 6
5       N1 = 512
6       DO 9 K2 = 1, N, 512
          S1 = 0.0E+00
          KK = K2 + N1 - 1
          DO 8 K1 = K2, KK, 32
            S0 = 0.0E+00
            KKK = K1 + N0 - 1
            DO 7 K0 = K1, KKK, 2
              S0 = S0 + F ( A + FLOAT ( K0 ) * H )
7           CONTINUE
            S1 = S1 + S0
8         CONTINUE
          S = S + S1
9       CONTINUE
        RM(K+1) = 2.0E+00 * H * S
C
C  END CALCULATION OF RECTANGLE VALUE.
C
        R = 4.0E+00
C
C  FORM ROMBERG TABLE FROM RECTANGLE VALUES.
C
        DO 10 J = 1, K
          L = K + 1 - J
          RM(L) = RM(L+1) + ( RM(L+1) - RM(L) ) / ( R - 1.0E+00 )
          R = 4.0E+00 * R
10      CONTINUE

        ERR = ABS ( T - RM(1) ) * 0.5E+00
C
C  CONVERGENCE TEST.
C
        IF ( ERR - EPS ) 12, 12, 11

11    CONTINUE

      K = 0

12    VAL = ( T + RM(1) ) * 0.5E+00
      N = N + 1
      MAXE = K

      RETURN
      END
