      SUBROUTINE DECOMP ( N, NDIM, A, IP )

c*********************************************************************72
c
c  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
c
c  INPUT..
c
c    N = ORDER OF MATRIX.
c    NDIM = DECLARED DIMENSION OF ARRAY A.
c    A = MATRIX TO BE TRIANGULARIZED.
c
c  OUTPUT..
c
c    A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR U.
c    A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR
c                     FACTOR, I-L.
c    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
c    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
c
c  USE 'SOLVE' TO OBTAIN SOLUTION OF LINEAR SYSTEM.
c  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
c  IF IP(N) = 0, A IS SINGULAR, SOLVE WILL DIVIDE BY ZERO.
c  INTERCHANGES FINISHED IN U, ONLY PARTLY IN L.
c
      REAL A(NDIM,NDIM), T
      INTEGER IP(NDIM)

      IP(N) = 1
      DO 6 K = 1, N
        IF ( K .EQ. N ) GO TO 5
        KP1 = K + 1
        M = K
        DO 1 I = KP1, N
          IF ( ABS ( A(I,K)) .GT. ABS ( A(M,K) ) ) M = I
1       CONTINUE
        IP(K) = M
        IF ( M .NE. K ) IP(N) = -IP(N)
        T = A(M,K)
        A(M,K) = A(K,K)
        A(K,K) = T
        IF ( T .EQ. 0. ) GO TO 5
        DO 2 I = KP1, N
2         A(I,K) = -A(I,K) / T
        DO 4 J = KP1, N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF ( T .EQ. 0. ) GO TO 4
          DO 3 I = KP1, N
3           A(I,J) = A(I,J) + A(I,K) * T
4       CONTINUE
5       IF ( A(K,K) .EQ. 0. ) IP(N) = 0
6     CONTINUE

      RETURN
      END
      SUBROUTINE SOLVE ( N, NDIM, A, B, IP )

c*********************************************************************72
c
c  SOLUTION OF LINEAR SYSTEM, A*X = B.
c
c  INPUT..
c
c    N = ORDER OF MATRIX.
c    NDIM = DECLARED DIMENSION OF ARRAY A.
c    A = TRIANGULARIZED MATRIX OBTAINED FROM 'DECOMP'.
c    B = RIGHT HAND SIDE VECTOR.
c    IP = PIVOT VECTOR OBTAINED FROM 'DECOMP'.
c  DO NOT USE IF DECOMP HAS SET IP(N) = 0.
c
c  OUTPUT..
c
c    B = SOLUTION VECTOR, X.
c
      REAL A(NDIM,NDIM), B(NDIM), T
      INTEGER IP(NDIM)

      IF ( N .EQ. 1 ) GO TO 9
      NM1 = N - 1
      DO 7 K = 1, NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 7 I = KP1, N
7         B(I) = B(I) + A(I,K) * T
      DO 8 KB = 1, NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K) / A(K,K)
        T = -B(K)
        DO 8 I = 1,KM1
8         B(I) = B(I) + A(I,K) * T
9       B(1) = B(1) / A(1,1)

      RETURN
      END
