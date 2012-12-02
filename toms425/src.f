C Incorporates the changes suggested in the remark by Page
C CACM 17(6) June 1974 p325
      SUBROUTINE RNVR ( X, A, NV, NI, IENT, IARG )
C  GENERATES A RANDOM NORMAL VECTOR (M,S)
C  A      INPUT COVARIANCE MATRIX, CONDITIONAL MOMENTS RETURN
C  X,B,C  WORK ARRAYS, RETURN VECTOR OF RANDOM NORMAL VARIATES IN X.
C  NV,NI  ORDER OF COVARIANCE MATRIX, ORDER OF ARRAY.
C  IENT   -1, INITIAL ENTRY.
C          0, RETURN IF NOT POSITIVE DEFINITE.
C          1, RETURN IF POSITIVE DEFINITE.
C  IARG    ARGUMENT FOR RANDOM NUMBER GENERATOR.
C     .. Scalar Arguments ..
      INTEGER IARG,IENT,NI,NV
C     ..
C     .. Array Arguments ..
      REAL A(NI,NI),X(NI)
C     ..
C     .. Local Scalars ..
      REAL T
      INTEGER I,J,K,NA,NB
C     ..
C     .. External Functions ..
      REAL RNOR
      EXTERNAL RNOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT

      IF ( IENT ) 1, 9, 6
C  ***  COMPUTE CONDITIONAL MOMENTS
1     NA = NV - 1
      DO 4 K = 1, NA
        T = A(K,K)
        IF ( T ) 10, 10, 2
2       NB = K + 1
        A(K,K) = SQRT ( T )
        DO 3 I = NB, NV
          A(I,K) = A(K,I) / T
3       CONTINUE
        DO 42 I = NB, NB
          DO 41 J = I, NV
            A(I,J) = A(I,J) - A(I,K) * A(K,J)
41        CONTINUE
42      CONTINUE
4     CONTINUE
      IF ( A(NV,NV) ) 10, 10, 5
5     IENT = 1
      A(NV,NV) = SQRT ( A(NV,NV) )
C  ***  COMPUTE A RANDOM VECTOR
6     DO 7 I = 1, NV
        X(I) = RNOR ( IARG ) * A(I,I)
7     CONTINUE
      DO 8 I = 2, NV
        NB = NV - I + 1
        DO 81 J = 1, NB
          X(NB+1) = X(NB+1) + A(NB+1,J) * X(J)
81      CONTINUE
8     CONTINUE
9     RETURN
10    IENT = 0
      RETURN
      END
      REAL FUNCTION RNOR ( IR )
C  GENERATES A RANDOM NORMAL NUMBER (0,1)
C  IARG IS A LARGE ODD INTEGER FOR A BEGINNING ARGUMENT.
C  REQUIRES FUNCTION RN WHICH GENERATES A UNIFORM RANDOM NUMBER 0-1.
C     .. Scalar Arguments ..
      INTEGER IR
C     ..
C     .. Local Scalars ..
      REAL GO2,S,X,Y
      INTEGER I
C     ..
C     .. External Functions ..
      REAL RN
      EXTERNAL RN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,SQRT
      DATA I /0/,GO2/0.0/
      IF ( I .GT. 0 ) GO TO 30
10    X = 2.0 * RN ( IR ) - 1.0
      Y = 2.0 * RN ( IR ) - 1.0
      S = X * X + Y * Y
      IF ( S .GE. ( 1.0 ) ) GO TO 10
      S = SQRT ( -2.0 * LOG ( S ) / S )
      RNOR = X * S
      GO2 = Y * S
      I = 1
      GO TO 40
30    RNOR = GO2
      I = 0
40    RETURN
      END

