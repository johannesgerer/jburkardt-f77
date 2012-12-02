      SUBROUTINE SPNBVP ( A, B, NP, NK, X, XDP, EP, GT, KG, VP,
     &  VQ, VR, VG, MAT, VM )

c*********************************************************************72
C
C  THIS ALGORITHM COMPUTES ITERATIVELY A CUBIC SPLINE
C  APPROXIMATION TO THE SOLUTION OF THE DIFFERENTIAL EQUATION
C  X''(T)=P(T)X(T)+Q(T)X(G(T))+R(T) ON THE INTERVAL (A,B)
C  WITH BOUNDARY CONDITIONS GIVEN BY U(T) IF T.LE.A AND
C  V(T) IF T.GE.B.
C
C  A AND B ARE TWO REAL VARIABLES DEFINED AS ABOVE.
C
C  NP   AN INTEGER VARIABLE SPECIFYING THE NUMBER OF KNOTS
C       ON THE INTERVAL (A,B).
C
C  NK   AN INTEGER VARIABLE SPECIFYING THE NUMBER OF INTERIOR
C       KNOTS.  THUS NK=NP-2.  IT IS USED TO ESTABLISH THE
C       DIMENSION OF CERTAIN ARRAYS MENTIONED BELOW.
C
C  X    ON RETURN TO THE CALLING PROGRAM X WILL CONTAIN THE
C       VALUES OF THE APPROXIMATION TO THE SOLUTION AT THE
C       NK INTERIOR KNOTS.  THIS IS AN ARRAY OF DIMENSION
C       NK AND TYPE REAL.
C
C  XDP  ON RETURN, XDP CONTAINS THE QUANTITIES H*H/6.0
C       MULTIPLIED BY THE SECOND DERIVATIVE VALUES AT ALL THE
C       KNOTS.  XDP IS A REAL ARRAY OF DIMENSION NP.
C
C  EP   THIS REAL VARIABLE IS SET TO THE VALUE 1.0E-M IF WE
C       REQUIRE M IDENTICAL FIGURES IN SUCCESSIVE ITERATES.
C
C  GT   AN INTEGER ARRAY OF LENGTH NP WHICH ASSIGNS TO EACH
C       KNOT T SUB J AN INTEGER VALUE BETWEEN 1 AND 6.  THIS
C       VALUE DESIGNATES RESPECTIVELY THE CASES WHEN
C       G(T SUB J) IS 1) .LE. A, 2) .GE. B, 3) WITHIN EP OF
C       SOME KNOT VALUE, 4) IN THE FIRST SUBINTERVAL,
C       5) IN THE LAST SUBINTERVAL, AND 6) IN ANY OTHER
C       SUBINTERVAL.  GT(I+1) CORRESPONDS TO KNOT T SUB I.
C
C  KG   AN INTEGER ARRAY OF LENGTH NP WHICH ASSIGNS TO EACH
C       KNOT AN INTEGER BETWEEN 2 AND NP-1.  IF GT(I+1) = 3,
C       THEN KG(I+1) CONTAINS THE SUBSCRIPT OF THE KNOT
C       AT THE POINT G(T SUB I).  IF GT(I+1) = 6, THEN KG(I+1)
C       CONTAINS THE SUBSCRIPT OF THE KNOW AT THE RIGHT HAND
C       ENDPOINT OF THE SUBINTERVAL IN WHICH G(T SUB I) LIES.
C
C  VP, VQ, VR, AND VG ARE ALL REAL ARRAYS OF DIMENSION NP AND
C       CONTAINS THE VALUES OF THE FUNCTIONS P, Q, R AND G
C       RESPECTIVELY, EACH EVALUATED AT THE NP KNOTS.
C
C  MAT  IS A REAL NK BY NK ARRAY USED IN THE MATRIX EQUATION
C       (MAT)(X)=(VM) SET UP TO SOLVE FOR THE X SUB J VALUES
C       STORED IN ARRAY X.
C
C  VM   AN ARRAY OF LENGTH NK AND TYPE REAL USED AS
C       DESCRIBED ABOVE.
C
C  THE USER MUST SUPPLY REAL FUNCTION SUBPROGRAMS TO COMPUTE
C  THE FUNCTIONS U(T),V(T),P(T),Q(T),R(T) AND G(T) DEFINED AS
C  ABOVE.  HE MUST ALSO SUPPLY SUBPROGRAMS WHICH SOLVE THE
C  SYSTEM (MAT)*(X)=(VM).  THE ROUTINE LUDCMP(MAT,NK) IS TO
C  REAPLCE MAT BY ITS LU DECOMPOSITION.  THE ROUTINE
C  LUSUB(VM,MAT,X,NK) IS TO COMPUTE X WHEN VM AND THE LU
C  FORM OF MAT IS GIVEN.
C
      IMPLICIT NONE

      INTEGER NK
      INTEGER NP

      REAL A
      REAL APLSE
      REAL B
      REAL BMINE
      REAL C
      REAL CCC
      REAL COEF
      REAL EP
      REAL G
      INTEGER GT1
      INTEGER GTE
      INTEGER GT(NP)
      INTEGER GTNP
      INTEGER GTYP
      REAL H
      REAL HR
      REAL HS
      INTEGER J
      INTEGER JF
      INTEGER JJ
      INTEGER JJZ
      INTEGER JZ
      INTEGER K
      INTEGER KG(NP)
      INTEGER KNOT
      INTEGER KPR
      REAL MAT(NK,NK)
      INTEGER N
      INTEGER NK1
      INTEGER NNN
      REAL P
      REAL Q
      REAL R
      REAL RK
      REAL RKNOT
      REAL RN
      REAL T
      REAL TJ
      REAL TM
      REAL TSTVL1
      REAL TSTVL2
      REAL U
      REAL V
      REAL VDH
      REAL VG(NP)
      REAL VGE
      REAL VM(NK)
      REAL VP(NP)
      REAL VPA
      REAL VPB
      REAL VQ(NP)
      REAL VR(NP)
      REAL X(NK)
      REAL XA
      REAL XB
      REAL XDP(NP)
      REAL XDPA
      REAL XDPB
      REAL XGAA
      REAL XGAB
C
C  KPR IS PRINTER DEVICE.
C
      DATA KPR / 6 /

      C ( T ) = T * ( T * T - 1.0 )
C
C  INITIALIZATION.
C
      N = NP - 1
      RN = N
C
C  NK is presumably set by the user!
C
C     NK = N - 1

      DO 20 K = 1, NK
        DO 10 J = 1, NK
          MAT(K,J) = 0.0
10      CONTINUE
20    CONTINUE

      XA = U ( A )
      XB = V ( B )
C
C  INITIALIZE XDP TO ZERO (INITIAL SPLINE).
C
      DO 30 K = 1, NP
        XDP(K) = 0.0
30    CONTINUE
C
C  SET UP P, Q, R, G VECTORS.
C
      H = ( B - A ) / RN
      HS = H * H / 6.0
      HR = 1.0 / H

      DO 40 K = 1, NP
        RK = K - 1
        TM = A + RK * H
        VP(K) = P ( TM )
        VQ(K) = Q ( TM )
        VR(K) = R ( TM )
        VG(K) = G ( TM )
40    CONTINUE
C
C  SET UP *TYPE OF G VALUE* ARRAY AND KG ARRAY.
C
      APLSE = A + EP * ABS ( A )
      BMINE = B - EP * ABS ( B )
      DO 70 K = 1, NP
        GTE = 6
        VGE = VG(K)
        IF ( VGE .LT. A + H ) GTE = 4
        IF ( VGE .GT. B - H ) GTE = 5
        IF ( VGE .LE. APLSE ) GTE = 1
        IF ( VGE .GE. BMINE ) GTE = 2
        VDH = ( VGE - A ) / H
        KNOT = VDH + EP
        RKNOT = KNOT
        IF ( ( KNOT .LT. 1 ) .OR. ( KNOT .GT. NK ) ) GO TO 50
        IF ( ABS ( VDH - RKNOT ) .GT. EP ) GO TO 50
        GTE = 3
        KG(K) = KNOT
        GO TO 60
50      KG(K) = KNOT + 1
60      GT(K) = GTE
70    CONTINUE
C
C  PUT XSUBJ COEFFICIENTS INTO (MAT) AND INITIALIZE X TO ZERO.
C
      DO 90 J = 1, NK
        X(J) = 0.0
        IF ( J .EQ. 1 ) GO TO 80
        MAT(J,J-1) = 1.0 - HS * VP(J)
80      MAT(J,J) = -2.0 * ( 1.0 + 2.0 * HS * VP(J+1) )
        IF ( J .EQ. NK ) GO TO 90
        MAT(J,J+1) = 1.0 - HS * VP(J+2)
90    CONTINUE
C
C  ADD INTO (MAT) X SUB G SUB T COEFFICIENTS.
C
      DO 150 J = 1, NK
        DO 140 JJ = 1, 3
          JZ = JJ - 1
          JJZ = J + JZ
          COEF = HS * VQ(JJZ)
          IF ( JZ .EQ. 1 ) COEF = COEF * 4.0
          GTYP = GT(JJZ)
          GO TO ( 140, 140, 100, 110, 120, 130 ), GTYP
100       KNOT = KG(JJZ)
          MAT(J,KNOT) = MAT(J,KNOT) - COEF
          GO TO 140
110       MAT(J,1) = MAT(J,1) - COEF * ( VG(JJZ) - A ) * HR
          GO TO 140
120       MAT(J,NK) = MAT(J,NK) - COEF * ( B - VG(JJZ) ) * HR
          GO TO 140
130       KNOT = KG(JJZ)
          RKNOT = KNOT
          CCC = RKNOT + ( A - VG(JJZ) ) * HR
          MAT(J,KNOT-1) = MAT(J,KNOT-1) - COEF * CCC
          CCC = ( VG(JJZ) - A ) * HR - RKNOT + 1.0
          MAT(J,KNOT) = MAT(J,KNOT) - COEF * CCC
140     CONTINUE
150   CONTINUE
C
C  REPLACE (MAT) BY ITS LU DECOMPOSITION.
C
      CALL LUDCMP ( MAT, NK )
C
C  A SEQUENCE OF SPLINES (UP TO 20) IS NOT GENERATED.
C
      VPA = VP(1)
      VPB = VP(NP)

      DO 250 NNN = 1, 20
C
C  VECTOR VM I SNOW SET UP.
C
        DO 200 J = 1, NK
          VM(J) = ( VR(J) + 4.0 * VR(J+1) + VR(J+2) ) * HS
          DO 190 JJ = 1, 3
            JZ = JJ - 1
            JJZ = J + JZ
            GTYP = GT(JJZ)
            COEF = HS * VQ(JJZ)
            IF ( JZ .EQ. 1 ) COEF = COEF * 4.
            IF ( GTYP .EQ. 1 ) VM(J) = VM(J) + COEF * U(VG(JJZ))
            IF ( GTYP .EQ. 2 ) VM(J) = VM(J) + COEF * V(VG(JJZ))
            GO TO ( 190, 190, 190, 160, 170, 180 ), GTYP
160         TM = ( VG(JJZ) - A ) * HR
            TJ = 1.0 - TM
            CCC = TJ * XA + C ( TM ) * XDP(2) + C ( TJ ) * XDP(1)
            VM(J) = VM(J) + COEF * CCC
            GO TO 190
170         TJ = ( B - VG(JJZ) ) * HR
            TM = 1. - TJ
            CCC = TM * XB + C ( TM ) * XDP(NK+1) + C ( TJ ) * XDP(NK+1)
            VM(J) = VM(J) + COEF * CCC
            GO TO 190
180         KNOT = KG(JJZ)
            RKNOT = KNOT
            TJ = ( A - VG(JJZ) ) * HR + RKNOT
            TM = 1.0 - TJ
            CCC = C ( TM ) * XDP(KNOT+1) - C ( TJ ) * XDP(KNOT)
            VM(J) = VM(J) + COEF * CCC
190       CONTINUE
200     CONTINUE
        VM(1) = VM(1) - ( 1.0 - HS * VPA ) * U(A)
        VM(NK) = VM(NK) - ( 1.0 - HS * VPB ) * V(B)
C
C  THE NEW X ARRAY IS NOW COMPUTED.
C  THE ARRAY VP SERVES AS A WORK AREA.
C
        DO 210 JF = 1, NK
          VP(JF) = X(JF)
210     CONTINUE

        CALL LUSUB ( VM, MAT, X, NK )

        TSTVL1 = 0.0
        TSTVL2 = 0.0
        DO 220 JF = 1, NK
          TSTVL1 = TSTVL1 + ABS ( VP(JF) - X(JF) )
          TSTVL2 = TSTVL2 + ABS ( X(JF) )
220     CONTINUE
C
C  CALCULATION OF XDP AT A AND B.
C
        GT1 = GT(1)
        GTNP = GT(NP)
        IF ( GT1 .EQ. 1 ) XGAA = U ( VG(1) )
        IF ( GT1 .EQ. 2 ) XGAA = V ( VG(1) )
        IF ( GTNP .EQ. 1 ) XGAB = U ( VG(NP) )
        IF ( GTNP .EQ. 2 ) XGAB = V ( VG(NP) )
        CALL GAGB ( GT1, XGAA, KG(1), VG(1), X, XDP, A, B, NP, NK )
        CALL GAGB ( GTNP, XGAB, KG(NP), VG(NP), X, XDP, A, B, NP,
     &    NK )
        XDPA = ( VPA * XA + VG(1) * XGAA + VR(1) ) * HS
        XDPB = ( VPB * XB + VG(NP) * XGAB + VR(NP) ) * HS
C
C  SOLVE FOR XDP VALUES OF CURRENT SPLINE USING CONTINUITY
C  EQUATIONS.  VM AND VP ARE USED AS WORKING AREAS.
C
        VM(1) = XA + X(2) - 2.0 * X(1) - XDPA
        NK1 = NK - 1
        VM(NK) = XB + X(NK1) - 2.0 * X(NK) - XDPB
        DO 230 J = 2, NK1
          VM(J) = X(J-1) + X(J+1) - 2.0 * X(J)
230     CONTINUE

        CALL SOLVE ( VM, NK, VP, NP )

        XDP(1) = XDPA
        XDP(NP) = XDPB
        DO 240 J = 1, NK
          XDP(J+1) = VM(J)
240     CONTINUE

        IF ( TSTVL1 .LE. TSTVL2 * EP ) RETURN
        IF ( NNN .EQ. 20 ) WRITE ( KPR, 1000 )
1000    FORMAT ( 32H NO CONVERGENCE IN 20 ITERATIONS )

250   CONTINUE

      RETURN
      END
      SUBROUTINE GAGB ( GTYP, ANS, K, GV, X, XDP, A, B, NP, NK )

c*********************************************************************72

      IMPLICIT NONE

      INTEGER NK
      INTEGER NP

      REAL A
      REAL ANS
      REAL B
      REAL C
      INTEGER GTYP
      REAL GV
      REAL H
      INTEGER K
      REAL RK
      REAL RNKD
      REAL T
      REAL TJ
      REAL TM
      REAL U
      REAL V
      REAL X(NK)
      REAL XA
      REAL XB
      REAL XDP(NP)

      C ( T ) = T * ( T * T - 1. )
      RNKD = NK + 1
      RK = K
      XA = U ( A )
      XB = V ( B )
      H = ( B - A ) / RNKD
      GO TO ( 10, 20, 30, 40, 50, 60 ), GTYP
10    RETURN
20    RETURN
30    ANS = X(K)
      RETURN

40    TM = ( GV - A ) / H
      TJ = 1. - TM
      ANS = TM * X(1) + TJ * XA + C ( TM ) * XDP(2) + C ( TJ ) * XDP(1)
      RETURN

50    TJ = ( B - GV ) / H
      TM = 1. - TJ
      ANS = TM * XB + TJ * X(NK) + C(TM) * XDP(NK+2) + C(TJ) * XDP(NK+1)
      RETURN

60    TJ = ( A - GV ) / H + RK
      TM = 1. - TJ
      ANS = TM * X(K) + TJ * X(K-1) + C(TM) * XDP(K+1) + C(TJ) * XDP(K)

      RETURN
      END
      SUBROUTINE SOLVE ( D, NK, M, NP )

c*********************************************************************72

      IMPLICIT NONE

      INTEGER NK
      INTEGER NP

      REAL D(NK)
      REAL M(NP)
      INTEGER I
      INTEGER J
      INTEGER NK1

      NK1 = NK - 1

      M(NK) = 0.25
      DO 10 I = 1, NK1
        J = NK - I
        M(J) = 1.0 / ( 4.0 - M(J+1) )
        D(J) = D(J) - D(J+1) * M(J+1)
10    CONTINUE

      D(1) = D(1) * M(1)
      DO 20 I = 2, NK
        D(I) = ( D(I) - D(I-1) ) * M(I)
20    CONTINUE

      RETURN
      END
