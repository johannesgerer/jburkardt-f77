      SUBROUTINE POLYAN ( C, CM, N )

c*********************************************************************72
c
C  POLYAN OBTAINS INFORMATION ABOUT THE LOCATION
C  OF THE ROOTS OF A POLYNOMIAL BY USING
C  BOUND, RADIUS AND HRWTZR.
C  C IS AN N ELEMENT ARRAY CONTAINING THE COEFFICIENTS
C  NORMALIZED SO THAT THE LEADING COEFFICIENT (WHICH
C  IS NOT INCLUDED IN C) IS +1.0.
C  CM IS A WORKING ARRAY THE SAME SIZE AS C.
C  N = DEGREE OF POLYNOMIAL.
C
      DIMENSION C(N), CM(N)
      LOGICAL HRWTZR
C  TEST FOR ZERO ROOT
      IF ( C(N) .EQ. 0.0 ) GO TO 50
C  COEFFICIENTS FOR RECIPROCAL POLYNOMIAL ARE PUT IN CM.
      CM(N) = 1. / C(N)
      NM1 = N - 1
      DO 5 I = 1, NM1
        NI = N - I
5       CM(I) = CM(N) * C(NI)
      ROUT = RADIUS ( C, N )
      RIN = 1. / RADIUS ( CM, N )
      WRITE ( 6, 201 ) RIN, ROUT
201   FORMAT ( 40H ROOTS ARE IN AN ANNULUS OF INNER RADIUS,
     & E10.3, 17H AND OUTER RADIUS, E10.3 )
      RPU = BOUND ( C, N )
      IF ( RPU .NE. 0.0 ) GO TO 10
      WRITE ( 6, 202 )
202   FORMAT ( 33H THERE ARE NO REAL POSITIVE ROOTS )
      GO TO 20
10    RPL = 1. / BOUND ( CM, N )
      WRITE ( 6, 203 ) RPL, RPU
203   FORMAT 
     &  ( 40H THE POSITIVE ROOTS(IF ANY) ARE BETWEEN,
     &  E10.3, 4H AND, E10.3 )
C  COEFFICIENTS FOR NEGATIVE RECIPROCAL ARE PUT IN CM.
20    DO 25 I = 1, N, 2
25      CM(I) = - CM(I)
      RNU = BOUND ( CM, N )
      IF ( RNU .NE. 0.0 ) GO TO 30
      WRITE ( 6, 204 )
204   FORMAT
     &  ( 33H THERE ARE NO NEGATIVE REAL ROOTS )
      GO TO 40
C  COEFFICIENTS FOR NEGATIVE ROOTS ARE PUT IN CM.
30    X = -1.0
      DO 35 I = 1, N
        CM(I) = X * C(I)
35      X = -X
      RNU = -1. / RNU
      RNL = - BOUND ( CM, N )
      WRITE ( 6, 205 ) RNU, RNL
205   FORMAT 
     &  ( 44H THE REAL NEGATIVE ROOTS(IF ANY) ARE BETWEEN,
     &  E10.3, 4H AND, E10.3 )
40    IF ( HRWTZR ( C, N ) ) WRITE ( 6, 206 )
206   FORMAT 
     &  ( 44H THERE ARE NO ROOTS WITH POSITIVE REAL PARTS )
      IF ( HRWTZR ( CM, N ) ) WRITE ( 6, 207 )
207   FORMAT 
     &  ( 44H THERE ARE NO ROOTS WITH NEGATIVE REAL PARTS )
      RETURN
50    WRITE ( 6, 208 )
208   FORMAT
     &  ( 41H POLYNOMIAL HAS A ZERO ROOT-REDUCE DEGREE )
      RETURN
      END
      FUNCTION RADIUS ( C, N )

c*********************************************************************72
c
C  RADIUS RETURNS AN UPPER LIMIT FOR THE MODULUS
C  OF THE ROOTS OF AN N DEGREE POLYNOMIAL.
C
      DIMENSION C(N)
      RADIUS = ABS ( C(1) )
      DO 10 I = 2, N
10      IF ( ABS ( C(I) ) .GT. RADIUS ) RADIUS = ABS ( C(I) )
      RADIUS = 1. + RADIUS

      RETURN
      END
      FUNCTION BOUND ( C, N )

c*********************************************************************72
c
C  BOUND RETURNS AN UPPER LIMIT FOR THE
C  POSITIVE REAL ROOTS OF AN N DEGREE POLYNOMIAL.
C
      DIMENSION C(N)
      M = 0
      BOUND = 0.0
      DO 10 I = 1, N
        IF ( M .GT. 0 ) GO TO 10
        IF ( C(I) .LT. 0.0 ) M = I
10      IF ( C(I) .LT. BOUND ) BOUND = C(I)
      IF ( M .EQ. 0 ) RETURN
      BOUND = 1. + ( -BOUND )** ( 1. / FLOAT ( M ) )

      RETURN
      END
      LOGICAL FUNCTION HRWTZR ( C, N )

c*********************************************************************72
c
C  HRWTZR RETURNS .TRUE. IF ALL THE ROOTS HAVE
C  NEGATIVE REAL PARTS, OTHERWISE .FALSE. IS RETURNED.
C  IF A REAL PART IS ZERO, THEN .FALSE. IS RETURNED.
C
      DIMENSION C(N)
      HRWTZR = .FALSE.
C%%
C%% Modifications as suggested by Driessen and Hunt,
C%% CACM Volume 16, Number 9, page 579, September 1973.
C%%
      IF ( C(1) .LE. 0.0 .OR. C(N) .LE. 0. ) RETURN
      C1 = C(1)
      M = N - 1
      DO 30 I = 2, M
        DO 20 K = I, M, 2
20        C(K) = C(K) - C(K+1) / C1
        C1 = C(I) / C1
        IF ( C1 .LE. 0.0 ) RETURN
30    CONTINUE
      HRWTZR = .TRUE.

      RETURN
      END
