      COMMON/XYZ/A(30,50),B(30,50)
C
      DO 1 I=1,30
      DO 1 J=1,50
      DO 2 K=1,3
2     A(I,J)=A(I,J)*FUNCT(K)
1     CONTINUE
C
      DO 3 I=1,30
      DO 3 J=1,50
      IF (B(I,J).LT.0.0)GOTO3
      B(I,J)=A(I,J)
3     B(I,J)=B(I,J)F*UNCT2(83.7)
C
      RETURN
C
      END
