      SUBROUTINE JACOBI ( DEGREE, ALFA, BETA, X, F, FD, E, ED,
     &  FLAGF, FLAGFD )

c*********************************************************************72

      IMPLICIT NONE

      DOUBLE PRECISION A
      DOUBLE PRECISION ALF
      DOUBLE PRECISION ALFA
      DOUBLE PRECISION B
      DOUBLE PRECISION BET
      DOUBLE PRECISION BETA
      DOUBLE PRECISION C
      DOUBLE PRECISION D
      INTEGER DEGREE
      REAL E
      REAL ED
      REAL EG
      REAL E1
      REAL E2
      DOUBLE PRECISION F
      DOUBLE PRECISION FD
      INTEGER FLAGF
      INTEGER FLAGFD
      DOUBLE PRECISION G
      DOUBLE PRECISION H
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER M
      INTEGER N
      DOUBLE PRECISION P
      DOUBLE PRECISION PD
      DOUBLE PRECISION Q
      DOUBLE PRECISION QD
      REAL S
      DOUBLE PRECISION T1
      DOUBLE PRECISION T2
      DOUBLE PRECISION U
      DOUBLE PRECISION V
      DOUBLE PRECISION W
      DOUBLE PRECISION X
      REAL Y

      DIMENSION P(25)
      DIMENSION PD(25)
      DIMENSION Q(25)
      DIMENSION QD(25)
      DIMENSION U(25)
      DIMENSION V(25)
      DIMENSION W(25)

      SAVE A
      SAVE ALF
      SAVE B
      SAVE BET
      SAVE M
      SAVE U
      SAVE V
      SAVE W
      SAVE Y

      DATA ALF / -2.0D+00 /
      DATA BET / -2.0D+00 /
      DATA M / -2 /
      DATA Y / 2.2D-16 /

      IF ( DEGREE .EQ. 0 ) GO TO 10

      IF ( ( ALFA .NE. ALF ) 
     & .OR. ( BETA .NE. BET ) ) GO TO 1

      IF ( DEGREE .LE. M ) GO TO 5

      I = M
      K = DEGREE - 1
      M = DEGREE

      IF ( I - 2 ) 2, 3, 3
C
C  CALCULATE THE U(J), V(J), W(J) IN 
C  THE RECURRENCE RELATION
C  P(J) = P(J-1) * ( U(J) + V(J) * X ) - P(J-2) * W(J)
C
1     M = DEGREE
      ALF = ALFA
      BET = BETA
      A = ALF + BET
      B = ALF - BET
      U(1) = B / 2.0D+00
      V(1) = 1.0D+00 + A / 2.0D+00
      W(1) = 0.0D+00

      IF ( DEGREE .EQ. 1 ) GO TO 5

2     U(2) = A * B * ( A + 3.0D+00 ) 
     &  / ( 4.0D+00 * ( A + 2.0D+00 )**2 )
      V(2) = ( A + 3.0D+00 ) * ( A + 4.0D+00 ) 
     &  / ( 4.0D+00 * ( A + 2.0D+00 ) )
      W(2) = ( 1.0D+00 + ALF ) * ( 1.0D+00 + BET ) * ( A + 4.0D+00 )
      W(2) = W(2) / ( 2.0D+00 * ( A + 2.0D+00 )**2 )
      I = 2
      K = DEGREE - 1

3     IF ( ( DEGREE .EQ. 2 )
     &  .OR. ( I .GT. K ) ) GO TO 5

      DO 4 J = I, K
        A = 2 * J + 2
        D = ALF + BET
        A = A + D
        B = D * ( A - 1.0D+00 ) * ( ALF - BET )
        C = J + 1
        C = 2.0D+00 * C * ( A - 2.0D+00 ) * ( C + D )
        U(J+1) = B / C
        D = A * ( A - 1.0D+00 ) * ( A - 2.0D+00 )
        V(J+1) = D / C
        D = J
        A = 2.0D+00 * ( D + ALF ) * ( D + BET ) * A
        W(J+1) = A / C
4     CONTINUE
C
C  FIND THE STARTING VALES FOR J = 1 
C  AND J = 2 FOR USE IN THE RECURSION.
C
5     T1 = V(1) * X
      P(1) = U(1) + T1
      S = Y * DMAX1 ( DABS ( U(1) ), DABS ( T1 ) )
      Q(1) = P(1) + S
      PD(1) = V(1)
      QD(1) = V(1)

      IF ( DEGREE .EQ. 1 ) GO TO 7
      T1 = V(2) * X
      G = U(2) + T1
      EG = Y * DMAX1 ( DABS ( U(2) ), DABS ( T1 ) )
      H = G + EG
      T1 = G * P(1)
      E1 = DABS ( EG * P(1) )
      P(2) = T1 - W(2)
      S = Y * DABS ( W(2) )
      S = AMAX1 ( E1, S )
      Q(2) = H * Q(1) - W(2) + S
      PD(2) = G * PD(1) + V(2) * P(1)
      QD(2) = H * QD(1) + V(2) * Q(1)

      IF ( DEGREE .EQ. 2 ) GO TO 7
C
C  USE THE RECURSION.
C
      DO 6 J = 3, DEGREE
        T2 = V(J) * X
        G = U(J) + T2
        EG = Y * DMAX1 ( DABS ( U(J) ), DABS ( T2 ) )
        H = G + EG
        T1 = G * P(J-1)
        T2 = W(J) * P(J-2)
        E1 = DABS ( EG * P(J-1) )
        E2 = DABS ( T2 ) * Y
        P(J) = T1 - T2
        S = AMAX1 ( E1, E2 )
        Q(J)  = H * Q(J-1)  - W(J) * Q(J-2) + S
        PD(J) = G * PD(J-1) - W(J) * PD(J-2)
        QD(J) = H * QD(J-1) - W(J) * QD(J-2)
        PD(J) = PD(J) + V(J) * P(J-1)
        QD(J) = QD(J) + V(J) * Q(J-1)
6     CONTINUE
C
C  PREPARE THE OUTPUT.
C
7     N = DEGREE
      F = P(N)
      IF ( DABS ( F ) .LT. Y ) GO TO 8
      FLAGF = 0
      E = DABS ( 1.0D+00 - Q(N) / F )
      GO TO 9
8     E = DABS ( F - Q(N) )
      FLAGF = 1
9     FD = PD(N)
      IF ( DABS ( FD ) .LT. Y ) GO TO 11
      FLAGFD = 0
      ED = DABS ( 1.0D+00 - QD(N) / FD )
      GO TO 12
10    F = 1.0D+00
      E = 0.0D+00
      FD = 0.0D+00
      ED = 0.0D+00
      FLAGF = 2
      FLAGFD = 2
      GO TO 12
11    ED = DABS ( FD - QD(N) )
      FLAGFD = 1
12    RETURN
      END
