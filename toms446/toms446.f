      SUBROUTINE CHEBY ( NF, NPL, NPLMAX, N2, FUNCTN, X, FXJ, GC )

c*********************************************************************72
c
C  SIMULTANEOUS CHEBYSHEV ANALYSIS OF NF FUNCTIONS
C  COMPUTES A MATRIX, X, CONTAINING ONE CHEBYSHEV SERIES PER
C  COLUMN FOR A GIVEN NUMBER OF FUNCTIONS, NF.  INPUT NFL,
C  THE NUMBER OF TERMS IN ALL SERIES, NPLMAX, THE ROW
C  DIMENSION OF X IN THE CALLING PROGRAM (MUST BE.GE.NPL),
C  N2, DIMENSION OF GC (MUST BE.GE.2*(NPL-1)), AND FUNCTN,
C  THE NAME OF USER SUBROUTINE WHICH DEFINES THE NF
C  FUNCTIONS.  FXJ AND GC ARE WORK SPACE.
C  AN EXAMPLE OF SUCH A SUBROUTINE IS AS FOLLOWS
C  SUBROUTINE FUNCTN ( A, VAL )
C  DOUBLE PRECISION A, VAL(2)
C  VAL(1) = DSIN(A)
C  VAL(2) = DCOS(A)
C  RETURN
C  END
c
      DOUBLE PRECISION X(NPLMAX,NF), FXJ(NF), GC(N2), ENW, XJ,
     &  FK, PEN, FAC

      DO 20 K = 1, NPL
        DO 10 J = 1, NF
          X(K,J) = 0.D0
10      CONTINUE
20    CONTINUE
      N = NPL - 1
      ENN = N
      PEN = 3.1415926535897932D0 / ENN
      DO 30 K = 1, N2
        FK = K - 1
        GC(K) = DCOS ( FK * PEN )
30    CONTINUE
      DO 80 J = 1, NPL
        XJ = GC(J)
        CALL FUNCTN ( XJ, FXJ )
        IF ( J .NE. 1 .AND. J .NE. NPL ) GO TO 50
        DO 40 K = 1, NF
          FXJ(K) = .5D0 * FXJ(K)
40      CONTINUE
50      DO 70 L = 1, NPL
          LM = MOD ( ( L - 1 ) * ( J - 1 ), N2 ) + 1
          DO 60 K = 1, NF
            X(L,K) = X(L,K) + FXJ(K) * GC(LM)
60        CONTINUE
70      CONTINUE
80    CONTINUE
      FAC = 2.0D0 / ENN
      DO 100 K = 1, NPL
        DO 90 J = 1, NF
          X(K,J) = FAC * X(K,J)
90      CONTINUE
100   CONTINUE
      RETURN
      END
      SUBROUTINE MLTPLY ( XX, X2, NPL, X3 )

c*********************************************************************72
c
C  MULTIPLIES TWO GIVEN CHEBYSHEV SERIES, XX AND X2, WITH
C  NPL TERMS TO PRODUCE AN OUTPUT CHEBYSHEV SERIES, X3.
c
      DOUBLE PRECISION XX(NPL), X2(NPL), X3(NPL), EX

      DO 10 K = 1, NPL
        X3(K) = 0.0D0
10    CONTINUE
      DO 30 K = 1, NPL
        EX = 0.0D0
        MM = NPL - K + 1
        DO 20 M = 1, MM
          L = M + K - 1
          EX = EX + XX(M) * X2(L) + XX(L) * X2(M)
20      CONTINUE
        X3(K) = 0.5D0 * EX
30    CONTINUE
      X3(1) = X3(1) - 0.5D0 * XX(1) * X2(1)
      DO 50 K = 3, NPL
        EX = 0.0D0
        MM = K - 1
        DO 40 M = 2, MM
          L = K - M + 1
          EX = EX + XX(M) * X2(L)
40      CONTINUE
        X3(K) = 0.5D0 * EX + X3(K)
50    CONTINUE
      RETURN
      END
      SUBROUTINE ECHEB ( X, COEF, NPL, FX )

c*********************************************************************72
c
C  EVALUATES THE VALUE FX(X) OF A GIVEN CHEBYSHEV SERIES,
C  COEF, WITH NPL TERMS AT A GIVEN VALUE OF X BETWEEN
C  -1. AND 1.
c
      DOUBLE PRECISION COEF(NPL), X, FX, BR, BRPP, BRP2

      BR = 0.0D0
      BRPP = 0.0D0
      DO 10 K = 1, NPL
        J = NPL - K + 1
        BRP2 = BRPP
        BRPP = BR
        BR = 2.0D00 * X * BRPP - BRP2 + COEF(J)
10    CONTINUE
      FX = 0.5D0 * ( BR - BRP2 )
      RETURN
      END
      SUBROUTINE EDCHEB ( X, COEF, NPL, FX )

c*********************************************************************72
c
C  EVALUATES THE VALUE FX(X) OF THE DERIVATIVE OF A
C  CHEBYSHEV SERIES, COEF, WITH NPL TERMS AT A GIVEN
C  VALUE OF X BETWEEN -1. AND 1.
c
      DOUBLE PRECISION COEF(NPL), X, FX, XJP2, XJPL, XJ, BJP2,
     &  BJPL, BJ, BF, DJ

      XJP2 = 0.0D0
      XJPL = 0.0D0
      BJP2 = 0.0D0
      BJPL = 0.0D0
      N = NPL - 1
      DO 10 K = 1, N
        J = NPL - K
        DJ = J
        XJ = 2.D0 * COEF(J+1) * DJ + XJP2
        BJ = 2.D0 * X * BJPL - BJP2 + XJ
        BF = BJP2
        BJP2 = BJPL
        BJPL = BJ
        XJP2 = XJPL
        XJPL = XJ
10    CONTINUE
      FX = .5D0 * ( BJ - BF )
      RETURN
      END
      SUBROUTINE DFRNT ( XX, NPL, X2 )

c*********************************************************************72
c
C  COMPUTES THE DERIVATIVE CHEBYSHEV SERIES, X2, OF A GIVEN
C  CHEBYSHEV SERIES, XX, WITH NPL TERMS.
C  TO REPLACE A SERIES X BY ITS DERIVATIVE, USE
C  CALL DFRNT ( X, NPL, X )
c
      DOUBLE PRECISION XX(NPL), XXN, XXL, XN, DL, X2(NPL)

      DN = NPL - 1
      XXN = XX(NPL-1)
      X2(NPL-1) = 2.D0 * XX(NPL) * DN
      X2(NPL) = 0.D0
      DO 10 K = 3, NPL
        L = NPL - K + 1
        DL = L
        XXL = XX(L)
        X2(L) = X2(L+2) + 2.D0 * XXN * DL
        XXN = XXL
10    CONTINUE
      RETURN
      END
      SUBROUTINE NTGRT ( XX, NPL, X2 )

c*********************************************************************72
c
C  COMPUTES THE INTEGRAL CHEBYSHEV SERIES, X2, OF A GIVEN
C  CHEBYSHEV SERIES, XX, WITH NPL TERMS.
C  TO REPLACE A SERIES X BY ITS INTEGRAL, USE
C  CALL NTGRT ( X, NPL, X )
c
      DOUBLE PRECISION XX(NPL), XPR, TERM, DK, X2(NPL)

      XPR = XX(1)
      X2(1) = 0.0D0
      N = NPL - 1
      DO 10 K = 2, N
        DK = K - 1
        TERM = ( XPR - XX(K+1) ) / ( 2.D0 * DK )
        XPR = XX(K)
        X2(K) = TERM
10    CONTINUE
      DK = N
      X2(NPL) = XPR / ( 2.D0 * DK )
      RETURN
      END
      SUBROUTINE INVERT ( X, XX, NPL, NET, XNVSE, WW, W2 )

c*********************************************************************72
c
C  COMPUTES THE INVERSE CHEBYSHEV SERIES, XNVSE, GIVEN A
C  CHEBYSHEV SERIES, X, A FIRST APPROXIMATION CHEBYSHEV
C  SERIES, XX, WITH NPL TERMS, AND THE NUMBER OF
C  ITERATIONS, NET.  THE SUBROUTINE USES THE EULER METHOD
C  AND COMPUTES ALL POWERS EPS**K UP TO K=2**(NET+1),
C  WHERE EPS = 1 - X * ( XX INVERSE ).  WW AND W2 ARE WORK SPACE.
C  SUBROUTINES USED - MLTPLY
c
      DOUBLE PRECISION X(NPL), XX(NPL), XNVSE(NPL), WW(NPL),
     &  W2(NPL)

      CALL MLTPLY ( X, XX, NPL, WW )
      WW(1) = 2.D0 - WW(1)
      DO 10 K = 2, NPL
        WW(K) = -WW(K)
10    CONTINUE
      CALL MLTPLY ( WW, WW, NPL, W2 )
      WW(1) = 2.D0 * WW(1)
      DO 40 K = 1, NET
        CALL MLTPLY ( WW, W2, NPL, XNVSE )
        DO 20 J = 1, NPL
          WW(J) = WW(J) + XNVSE(J)
20      CONTINUE
        CALL MLTPLY ( W2, W2, NPL, XNVSE )
        DO 30 J = 1, NPL
          W2(J) = XNVSE(J)
30      CONTINUE
40    CONTINUE
      CALL MLTPLY ( WW, XX, NPL, XNVSE )
      RETURN
      END
      SUBROUTINE BINOM ( X, XX, NPL, M, NT, XA, WW, W2, W3 )

c*********************************************************************72
c
C  COMPUTES THE BINOMIAL EXPANSION SERIES, XA, FOR (-1/M)
C  POWER OF A GIVEN CHEBYSHEV SERIES, X, WITH NPL TERMS,
C  WHERE M IS A POSITIVE INTEGER.  XX IS A GIVEN INITIAL
C  APPROXIMATION TO X**(-1/M).  NT IS A GIVEN NUMBER OF
C  TERMS IN BINOMIAL SERIES.  WW, W2 AND W3 ARE WORK SPACE.
C  SUBROUTINES USED - MLTPLY.
c
      DOUBLE PRECISION X(NPL), XX(NPL), XA(NPL), WW(NPL),
     &  W2(NPL), W3(NPL), ALFA, COEF, DM, DKMM, DKM2

      DM = M
      ALFA = -1.D0 / DM
      DO 10 J = 1, NPL
        WW(J) = X(J)
10    CONTINUE
      DO 30 K = 1, M
        CALL MLTPLY ( WW, XX, NPL, W2 )
        DO 20 J = 1, NPL
          WW(J) = W2(J)
20      CONTINUE
30    CONTINUE
      WW(1) = WW(1) - 2.D0
      XA(1) = 2.D0
      DO 40 J = 2, NPL
        XA(J) = 0.0D0
        W3(J) = 0.D0
40    CONTINUE
      W3(1) = 2.D0
      DO 60 K = 2, NT
        DKMM = K - 1
        DKM2 = K - 2
        COEF = ( ALFA - DKM2 ) / DKMM
        CALL MLTPLY ( W3, WW, NPL, W2 )
        DO 50 J = 1, NPL
          W3(J) = W2(J) * COEF
          XA(J) = XA(J) + W3(J)
50      CONTINUE
60    CONTINUE
      CALL MLTPLY ( XA, XX, NPL, W2 )
      DO 70 J = 1, NPL
        XA(J) = W2(J)
70    CONTINUE
      RETURN
      END
      SUBROUTINE XALFA2 ( X, XX, NPL, M, MAXET, EPSLN, NET, WW,
     &  W2 )

c*********************************************************************72
c
C  REPLACES A GIVEN INITIAL APPROXIMATION CHEBYSHEV SERIES,
C  XX, BY A GIVEN CHEBYSHEV SERIES, X, WITH NPL TERMS,
C  RAISED TO THE (-1/M) POWER, WHERE M IS AN INTEGER.
C  INPUT MAXET, MAXIMUM ALLOWED NUMBER OF ITERATIONS, AND
C  EPSLN, REQUIRED PRECISION EPSILON.  OUTPUT ARGUMENT,
C  NET, IS NUMBER OF ITERATIONS PERFORMED.  IF MAXET=NET,
C  MAY BE DIVERGENCE.  WW AND W2 ARE WORK SPACE.
C  CONVERGENCE IS QUADRATIC.
C  SUBROUTINES USED - MULTIPLY.
c
      DOUBLE PRECISION X(NPL), XX(NPL), WW(NPL), W2(NPL),
     & EPSLN, DALFA, DM, S, TDMM

      DM = M
      DALFA = 1.D0 / DM
      TDMM = 2.D0 * ( DM + 1.D0 )
      DO 60 JX = 1, MAXET
        DO 10 L = 1, NPL
          WW(L) = X(L)
10      CONTINUE
        DO 30 K = 1, M
          CALL MLTPLY ( WW, XX, NPL, W2 )
          DO 20 L = 1, NPL
            WW(L) = W2(L)
20        CONTINUE
30      CONTINUE
        S = -2.D0
        DO 40 L = 1, NPL
          S = S + DABS ( WW(L) )
          WW(L) = -WW(L)
40      CONTINUE
        WW(1) = WW(1) + TDMM
        CALL MLTPLY ( WW, XX, NPL, W2 )
        DO 50 L = 1, NPL
          XX(L) = W2(L) * DALFA
50      CONTINUE
        NET = JX
        IF ( DABS ( S ) .LT. EPSLN ) RETURN
60    CONTINUE
      RETURN
      END
      SUBROUTINE XALFA3 ( X, XX, NPL, M, MAXET, EPSLN, NET, WW,
     &  W2 )

c*********************************************************************72
c
C  REPLACES A GIVEN INITIAL APPROXIMATION CHEBYSHEV SERIES,
C  XX, BY A GIVEN CHEBYSHEV SERIES, X, WITH NPL TERMS,
C  RAISED TO THE (-1/M) POWER, WHERE M IS AN INTEGER.
C  INPUT MAXET, MAXIMUM ALLOWED NUMBER OF ITERATIONS, AND
C  EPSLN, REQUIRED PRECISION EPSILON.  OUTPUT ARGUMENT
C  NET, IS NUMBER OF ITERATIONS PERFORMED.  IF MAXET = NET,
C  MAY BE DIVERGENCE.  WW AND W2 ARE WORK SPACE.
C  CONVERGENCE IS OF ORDER THREE.
C  SUBROUTINES USED - MLTPLY
c
      DOUBLE PRECISION X(NPL), XX(NPL), WW(NPL), W2(NPL),
     &  EPSLN, DALFA, DM, S, TDMM, P5DML

      DM = M
      DALFA = 1.D0 / DM
      TDMM = 2.D0 * ( DM + 1.D0 )
      P5DML = .5D0 * ( DM + 1.D0 )
      DO 90 JX = 1, MAXET
        DO 10 L = 1, NPL
          WW(L) = X(L)
10      CONTINUE
        DO 30 K = 1, M
          CALL MLTPLY ( WW, XX, NPL, W2 )
          DO 20 L = 1, NPL
            WW(L) = W2(L)
20        CONTINUE
30      CONTINUE
        S = -2.D0
        DO 40 L = 1, NPL
          S = S + DABS ( WW(L) )
40      CONTINUE
        WW(1) = WW(1) - 2.D0
        DO 50 L = 1, NPL
          WW(L) = WW(L) * DALFA
50      CONTINUE
        CALL MLTPLY ( WW, WW, NPL, W2 )
        DO 60 L = 1, NPL
          WW(L) = -WW(L)
          W2(L) = W2(L) * P5DML
60      CONTINUE
        WW(1) = WW(1) + 2.D0
        DO 70 L = 1, NPL
          W2(L) = W2(L) + WW(L)
70      CONTINUE
        CALL MLTPLY ( W2, XX, NPL, WW )
        DO 80 L = 1, NPL
          XX(L) = WW(L)
80      CONTINUE
        NET = JX
        IF ( DABS ( S ) .LT. EPSLN ) RETURN
90    CONTINUE

      RETURN
      END
