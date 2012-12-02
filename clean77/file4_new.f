      SUBROUTINE DOBAD (X,Y)
      COMMON /XYZ/ A(10)
C
      DO 10 I = 1, 10
         IF (A(I).GT.0.0) RETURN
         A(I) = I
         Y = A(I)
   10 CONTINUE
   10 X = Y
C
      DO 20 I = 1, 10
         IF (A(I).GT.10.0) GO TO 20
         A(I) = 2*I
   20    Y = A(I)
   20 X = Y
C
      RETURN
C
      END
