      SUBROUTINE PTRAP ( A, B, N, FN, GN, VINT )

c*********************************************************************72
c
cc PTRAP carries out product-type trapezoidal quadrature.
C
C  THIS SUBROUTINE USES THE PRODUCT TYPE TRAPEZOIDAL RULE
C  COMPOUNDED N TIMES TO APPROXIMATE THE INTEGRAL FROM A TO B
C  OF THE FUNCTION FN(X) * GN(X).  FN AND GN ARE FUNCTION
C  SUBPROGRAMS WHICH MUST BE SUPPLIED BY THE USER.  THE
C  RESULT IS STORED IN VINT.
C
      DOUBLE PRECISION A, AG, AM(2,2), B, F(2), FN, G(2),
     &  GN, H, VINT, X, DBLE
      DATA AM(1,1), AM(2,2) /2 * 2.D0 /, AM(1,2), AM(2,1)
     &  / 2 * 1.D0/
      H = ( B - A ) / DBLE ( FLOAT ( N ) )
      VINT = 0.D0
      X = A
      F(2) = FN ( A )
      G(2) = GN ( A )
      DO 2 I = 1, N
        F(1) = F(2)
        G(1) = G(2)
        X = X + H
        F(2) = FN ( X )
        G(2) = GN ( X )
        DO 2 J = 1, 2
          AG = 0.D0
          DO 1 K = 1, 2
1           AG = AG + AM(J,K) * G(K)
2         VINT = VINT + F(J) * AG
      VINT = H * VINT / 6.D0
      RETURN
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
