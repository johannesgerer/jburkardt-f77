      SUBROUTINE PSIMP ( A, B, N, FN, GN, VINT )

c*********************************************************************72
C
C  THIS SUBROUTINE USES THE PRODUCT TYPE SIMPSON RULE
C  COMPOUNDED N TIMES TO APPROXIMATE THE INTEGRAL FROM A TO B
C  OF THE FUNCTION FN(X) * GN(X).  FN AND GN ARE FUNCTION
C  SUBPROGRAMS WHICH MUST BE SUPPLIED BY THE USER.  THE
C  RESULT IS STORED IN VINT.
C
      DOUBLE PRECISION A, AG, AM(3,3), B, F(3), FN, G(3),
     &  GN, H, VINT, X(2), DBLE

      DATA AM(1,1), AM(3,3) /2 * 4.D0 /, AM(1,2), AM(2,1),
     &  AM(2,3), AM(3,2) / 4 * 2.D0/, AM(1,3), AM(3,1)
     &  /2 * -1.D0/, AM(2,2) / 16.0D0 /

      H = ( B - A ) / DBLE ( FLOAT ( N ) )
      X(1) = A + H / 2.D0
      X(2) = A + H
      VINT = 0.D0
      F(3) = FN ( A )
      G(3) = GN ( A )
      DO 3 I = 1, N
        F(1) = F(3)
        G(1) = G(3)
        DO 1 J = 1, 2
          F(J+1) = FN ( X(J) )
          G(J+1) = GN ( X(J) )
1         X(J) = X(J) + H
        DO 3 J = 1, 3
          AG = 0.D0
          DO 2 K = 1, 3
2           AG = AG + AM(J,K) * G(K)
3       VINT = VINT + F(J) * AG
      VINT = H * VINT / 30.D0

      RETURN
      END
