      double precision FUNCTION PRAXIS(T0,MACHEP,H0,N,PRIN,X,F,FMIN)

c*********************************************************************72
C                             LAST MODIFIED 3/1/73
      double precision T0,MACHEP,H0,X(N),F,FMIN
      EXTERNAL F
      INTEGER PRIN
C
C     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
C     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
C     NOT REQUIRED.
C
C     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
C     "ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
C     CALCULATING DERIVATIVES" BY RICHARD P BRENT.
C
C     THE PARAMETERS ARE:
C     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
C              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
C              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
C     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
C              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
C              2.22D-16) FOR double precision ARITHMETIC ON THE IBM 360.
C     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
C              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
C              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
C              CONVERGENCE MAY BE SLOW.)
C     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
C              THE FUNCTION DEPENDS.
C     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
C              IF PRIN=0, NOTHING IS PRINTED.
C              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
C              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
C              PRINTED ONLY IF N IS AT MOST 4.
C              IF PRIN=2, THE SCALE FACTORS AND THE PRINCIPAL VALUES OF
C              THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
C              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
C              MINIMIZATIONS.
C              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
C              QUADRATIC FORM ARE ALSO PRINTED.
C     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
C              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
C     F(X,N)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE A double precision
C              FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
C     FMIN     IS AN ESTIMATE OF THE MINIMUM, USED ONLY IN PRINTING
C              INTERMEDIATE RESULTS.
C     THE APPROXIMATING QUADRATIC FORM IS
C              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
C     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
C              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
C     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
C     OF SECOND DIFFERENCES).  IF F HAS CONTINUOUS SECOND DERIVATIVES
C     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.
C
C     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
C     TO ZERO.
C     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
C     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.
C
      LOGICAL ILLC
      INTEGER NL,NF,KL,KT,KTM
      double precision S,SL,DN,DMIN,FX,F1,LDS,LDT
      double precision T,H,SF,DF,QF1,QD0,QD1,QA,QB,QC
      double precision M2,M4,SMALL,VSMALL,LARGE
      double precision VLARGE,SCBD,LDFAC,T2,DNI,VALUE
      double precision RANDOM,DSQRT,DABS
C
C.....IF N>20 OR IF N<20 AND YOU NEED MORE SPACE, CHANGE '20' TO THE
C     LARGEST VALUE OF N IN THE NEXT CARD, IN THE CARD 'IDIM=20', AND
C     IN THE DIMENSION STATEMENTS IN SUBROUTINES MINFIT,MIN,FLIN,QUAD.
C
      double precision D(20),Y(20),Z(20),Q0(20),Q1(20),V(20,20)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
C
C.....INITIALIZATION.....
C     MACHINE DEPENDENT NUMBERS:
C
      SMALL=MACHEP*MACHEP
      VSMALL=SMALL*SMALL
      LARGE=1.D0/SMALL
      VLARGE=1.D0/VSMALL
      M2=DSQRT(MACHEP)
      M4=DSQRT(M2)
C
C     HEURISTIC NUMBERS:
C     IF THE AXES MAY BE BADLY SCALED (WHICH IS TO BE AVOIDED IF
C     POSSIBLE), THEN SET SCBD=10.  OTHERWISE SET SCBD=1.
C     IF THE PROBLEM IS KNOWN TO BE ILL-CONDITIONED, SET ILLC=TRUE.
C     OTHERWISE SET ILLC=FALSE.
C     KTM IS THE NUMBER OF ITERATIONS WITHOUT IMPROVEMENT BEFORE THE
C     ALGORITHM TERMINATES.  KTM=4 IS VERY CAUTIOUS; USUALLY KTM=1
C     IS SATISFACTORY.
C
      SCBD=1.D0
      ILLC=.FALSE.
      KTM=1
C
      LDFAC=0.01D0
      IF (ILLC) LDFAC=0.1D0
      KT=0
      NL=0
      NF=1
      FX=F(X,N)
      QF1=FX
      T=SMALL+DABS(T0)
      T2=T
      DMIN=SMALL
      H=H0
      IF (H.LT.100*T) H=100*T
      LDT=H
C.....THE FIRST SET OF SEARCH DIRECTIONS V IS THE IDENTITY MATRIX.....
      DO 20 I=1,N
           DO 10 J=1,N
10              V(I,J)=0.0D0
20         V(I,I)=1.D0
      D(1)=0.D0
      QD0=0.D0
      DO 30 I=1,N
           Q0(I)=X(I)
30         Q1(I)=X(I)
      IF (PRIN.GT.0) CALL PRINT(N,X,PRIN,FMIN)
C
C.....THE MAIN LOOP STARTS HERE.....
40    SF=D(1)
      D(1)=0.D0
      S=0.D0
C
C.....MINIMIZE ALONG THE FIRST DIRECTION V(*,1).
C     FX MUST BE PASSED TO MIN BY VALUE.
      VALUE=FX
      CALL MIN(N,1,2,D(1),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
      IF (S.GT.0.D0) GO TO 50
           DO 45 I=1,N
45              V(I,1)=-V(I,1)
50    IF (SF.GT.0.9D0*D(1).AND.0.9D0*SF.LT.D(1)) GO TO 70
           DO 60 I=2,N
60              D(I)=0.D0
C
C.....THE INNER LOOP STARTS HERE.....
70    DO 170 K=2,N
           DO 75 I=1,N
75              Y(I)=X(I)
           SF=FX
           IF (KT.GT.0) ILLC=.TRUE.
80         KL=K
           DF=0.D0
C
C.....A RANDOM STEP FOLLOWS (TO AVOID RESOLUTION VALLEYS).
C     PRAXIS ASSUMES THAT RANDOM RETURNS A RANDOM NUMBER UNIFORMLY
C     DISTRIBUTED IN (0,1).
C
           IF(.NOT.ILLC) GO TO 95
                DO 90 I=1,N
                     S=(0.1D0*LDT+T2*(10**KT))*(RANDOM(N)-0.5D0)
                     Z(I)=S
                     DO 85 J=1,N
85                        X(J)=X(J)+S*V(J,I)
90              CONTINUE
                FX=F(X,N)
                NF=NF+1
C
C.....MINIMIZE ALONG THE "NON-CONJUGATE" DIRECTIONS V(*,K),...,V(*,N)
C
95         DO 105 K2=K,N
                SL=FX
                S=0.D0
                VALUE=FX
                CALL MIN(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
                IF (ILLC) GO TO 97
                     S=SL-FX
                     GO TO 99
97              S=D(K2)*((S+Z(K2))**2)
99              IF (DF.GT.S) GO TO 105
                     DF=S
                     KL=K2
105        CONTINUE
           IF (ILLC.OR.(DF.GE.DABS((100*MACHEP)*FX))) GO TO 110
C
C.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
C     ILLC=TRUE AND START THE INNER LOOP AGAIN.....
C
           ILLC=.TRUE.
           GO TO 80
110        IF (K.EQ.2.AND.PRIN.GT.1) CALL VCPRNT(1,D,N)
C
C.....MINIMIZE ALONG THE "CONJUGATE" DIRECTIONS V(*,1),...,V(*,K-1)
C
           KM1=K-1
           DO 120 K2=1,KM1
           S=0
           VALUE=FX
           CALL MIN(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
120        CONTINUE
           F1=FX
           FX=SF
           LDS=0
           DO 130 I=1,N
                SL=X(I)
                X(I)=Y(I)
                SL=SL-Y(I)
                Y(I)=SL
130             LDS=LDS+SL*SL
           LDS=DSQRT(LDS)
           IF (LDS.LE.SMALL) GO TO 160
C
C.....DISCARD DIRECTION V(*,KL).
C     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE "NON-CONJUGATE"
C     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
C
           KLMK=KL-K
           IF (KLMK.LT.1) GO TO 141
           DO 140 II=1,KLMK
                I=KL-II
                DO 135 J=1,N
135                  V(J,I+1)=V(J,I)
140             D(I+1)=D(I)
141        D(K)=0
           DO 145 I=1,N
145             V(I,K)=Y(I)/LDS
C
C.....MINIMIZE ALONG THE NEW "CONJUGATE" DIRECTION V(*,K), WHICH IS
C     THE NORMALIZED VECTOR:  (NEW X) - (0LD X).....
C
           VALUE=F1
           CALL MIN(N,K,4,D(K),LDS,VALUE,.TRUE.,F,X,T,MACHEP,H)
           IF (LDS.GT.0.D0) GO TO 160
                LDS=-LDS
                DO 150 I=1,N
150                  V(I,K)=-V(I,K)
160        LDT=LDFAC*LDT
           IF (LDT.LT.LDS) LDT=LDS
           IF (PRIN.GT.0) CALL PRINT(N,X,PRIN,FMIN)
           T2=0.D0
           DO 165 I=1,N
165             T2=T2+X(I)**2
           T2=M2*DSQRT(T2)+T
C
C.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
C     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
C
           IF (LDT.GT.(0.5*T2)) KT=-1
           KT=KT+1
           IF (KT.GT.KTM) GO TO 400
170   CONTINUE
C.....THE INNER LOOP ENDS HERE.
C
C     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.
C
171   CALL QUAD(N,F,X,T,MACHEP,H)
      DN=0.D0
      DO 175 I=1,N
           D(I)=1.D0/DSQRT(D(I))
           IF (DN.LT.D(I)) DN=D(I)
175   CONTINUE
      IF (PRIN.GT.3) CALL MAPRNT(1,V,IDIM,N)
      DO 180 J=1,N
           S=D(J)/DN
           DO 180 I=1,N
180             V(I,J)=S*V(I,J)
C
C.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....
C
      IF (SCBD.LE.1.D0) GO TO 200
           S=VLARGE
           DO 185 I=1,N
                SL=0.D0
                DO 182 J=1,N
182                  SL=SL+V(I,J)*V(I,J)
                Z(I)=DSQRT(SL)
                IF (Z(I).LT.M4) Z(I)=M4
                IF (S.GT.Z(I)) S=Z(I)
185        CONTINUE
           DO 195 I=1,N
                SL=S/Z(I)
                Z(I)=1.D0/SL
                IF (Z(I).LE.SCBD) GO TO 189
                     SL=1.D0/SCBD
                     Z(I)=SCBD
189             DO 190 J=1,N
190                  V(I,J)=SL*V(I,J)
195        CONTINUE
C
C.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
C     THE MAIN LOOP.
C     FIRST TRANSPOSE V FOR MINFIT:
C
200   DO 220 I=2,N
           IM1=I-1
           DO 210 J=1,IM1
                S=V(I,J)
                V(I,J)=V(J,I)
210             V(J,I)=S
220   CONTINUE
C
C.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
C     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
C     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
C     NUMBER.....
C
      CALL MINFIT(IDIM,N,MACHEP,VSMALL,V,D)
C
C.....UNSCALE THE AXES.....
C
      IF (SCBD.LE.1.D0) GO TO 250
           DO 230 I=1,N
                S=Z(I)
                DO 225 J=1,N
225                  V(I,J)=S*V(I,J)
230        CONTINUE
           DO 245 I=1,N
                S=0.D0
                DO 235 J=1,N
235                  S=S+V(J,I)**2
                S=DSQRT(S)
                D(I)=S*D(I)
                S=1/S
                DO 240 J=1,N
240                  V(J,I)=S*V(J,I)
245        CONTINUE
C
250   DO 270 I=1,N
           DNI=DN*D(I)
           IF (DNI.GT.LARGE) GO TO 265
                IF (DNI.LT.SMALL) GO TO 260
                     D(I)=1/(DNI*DNI)
                     GO TO 270
260             D(I)=VLARGE
                GO TO 270
265        D(I)=VSMALL
270   CONTINUE
C
C.....SORT THE EIGENVALUES AND EIGENVECTORS.....
C
      CALL SORT(IDIM,N,D,V)
      DMIN=D(N)
      IF (DMIN.LT.SMALL) DMIN=SMALL
      ILLC=.FALSE.
      IF (M2*D(1).GT.DMIN) ILLC=.TRUE.
      IF (PRIN.GT.1.AND.SCBD.GT.1.D0) CALL VCPRNT(2,Z,N)
      IF (PRIN.GT.1) CALL VCPRNT(3,D,N)
      IF (PRIN.GT.3) CALL MAPRNT(2,V,IDIM,N)
C.....THE MAIN LOOP ENDS HERE.....
C
      GO TO 40
C
C.....RETURN.....
C
400   IF (PRIN.GT.0) CALL VCPRNT(4,X,N)
      PRAXIS=FX
      RETURN
      END
      double precision FUNCTION FLIN (N,J,L,F,X,NF)

c*********************************************************************72
      IMPLICIT double precision (A-H,O-Z)
      double precision L,X(N)
      DIMENSION V(20,20),Q0(20),Q1(20)
C...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
C   BY THE SUBROUTINE MIN...
      COMMON /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
      DIMENSION T(20)
      IF (J .EQ. 0) GO TO 2
C...THE SEARCH IS LINEAR...
      DO 1 I=1,N
1        T(I) = X(I) + L*V(I,J)
      GO TO 4
C...THE SEARCH IS ALONG A PARABOLIC SPACE CURVE...
2     QA = (L*(L - QD1))/(QD0*(QD0 + QD1))
      QB = ((L + QD0)*(QD1 - L))/(QD0*QD1)
      QC = (L*(L + QD0))/(QD1*(QD0 + QD1))
      DO 3 I=1,N
3        T(I) = (QA*Q0(I) + QB*X(I)) + QC*Q1(I)
C...THE FUNCTION EVALUATION COUNTER NF IS INCREMENTED...
4     NF = NF + 1
      FLIN = F(T,N)
      RETURN
      END
      SUBROUTINE MAPRNT(OPTION,V,M,N)

c*********************************************************************72
      double precision V(M,N)
      INTEGER OPTION,LOW,UPP
C...THE SUBROUTINE MAPRNT PRINTS THE COLUMNS OF THE NXN MATRIX V
C   WITH A HEADING AS SPECIFIED BY OPTION.
C   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM...
      LOW = 1
      UPP = 5
      GO TO (1,2),OPTION
1     WRITE (6,101)
101   FORMAT (/' THE NEW DIRECTIONS ARE:')
      GO TO 3
2     WRITE (6,102)
102   FORMAT (' AND THE PRINCIPAL AXES:')
3     IF (N.LT.UPP) UPP = N
      DO 4 I=1,N
4        WRITE (6,104) (V(I,J),J=LOW,UPP)
      LOW = LOW + 5
      IF (N.LT.LOW) RETURN
      UPP = UPP + 5
      WRITE (6,103)
      GO TO 3
103   FORMAT (' ')
104   FORMAT (E32.14,4E25.14)
      END
      SUBROUTINE MIN(N,J,NITS,D2,X1,F1,FK,F,X,T,MACHEP,H)

c*********************************************************************72
      IMPLICIT double precision (A-H,O-Z)
      EXTERNAL F
      LOGICAL FK
      double precision MACHEP,X(N),LDT
      DIMENSION V(20,20),Q0(20),Q1(20)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
C...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
C   J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
C   DEFINED BY Q0,Q1,X.
C   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
C   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
C   ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
C   FOUND.
C   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
C   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
C   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
C   THE INTERVAL.
      LOGICAL DZ
      double precision M2,M4
      SMALL = MACHEP**2
      M2 = DSQRT(MACHEP)
      M4 = DSQRT(M2)
      SF1 = F1
      SX1 = X1
      K = 0
      XM = 0.D0
      FM = FX
      F0 = FX
      DZ = D2.LT.MACHEP
C...FIND THE STEP SIZE...
      S = 0.D0
      DO 1 I=1,N
1        S = S + X(I)**2
      S = DSQRT(S)
      TEMP = D2
      IF (DZ) TEMP = DMIN
      T2 = M4*DSQRT(DABS(FX)/TEMP + S*LDT) + M2*LDT
      S = M4*S + T
      IF (DZ.AND.T2.GT.S) T2 = S
      T2 = DMAX1(T2,SMALL)
      T2 = DMIN1(T2,.01D0*H)
      IF (.NOT.FK.OR.F1.GT.FM) GO TO 2
      XM = X1
      FM = F1
2     IF (FK.AND.DABS(X1).GE.T2) GO TO 3
      TEMP=1.D0
      IF (X1.LT.0.D0) TEMP=-1.D0
      X1=TEMP*T2
      F1 = FLIN(N,J,X1,F,X,NF)
3     IF (F1.GT.FM) GO TO 4
      XM = X1
      FM = F1
4     IF (.NOT.DZ) GO TO 6
C...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
      X2 = -X1
      IF (F0.GE.F1) X2 = 2.D0*X1
      F2 = FLIN(N,J,X2,F,X,NF)
      IF (F2.GT.FM) GO TO 5
         XM = X2
         FM = F2
5     D2 = (X2*(F1 - F0)-X1*(F2 - F0))/((X1*X2)*(X1 - X2))
C...ESTIMATE THE FIRST DERIVATIVE AT 0...
6     D1 = (F1 - F0)/X1 - X1*D2
      DZ = .TRUE.
C...PREDICT THE MINIMUM...
      IF (D2.GT.SMALL) GO TO 7
         X2 = H
         IF (D1.GE.0.D0) X2 = -X2
         GO TO 8
7        X2 = (-.5D0*D1)/D2
8     IF (DABS(X2).LE.H) GO TO 11
         IF (X2) 9,9,10
9        X2 = -H
         GO TO 11
10       X2 = H
C...EVALUATE F AT THE PREDICTED MINIMUM...
11    F2 = FLIN(N,J,X2,F,X,NF)
      IF (K.GE.NITS.OR.F2.LE.F0) GO TO 12
C...NO SUCCESS, SO TRY AGAIN...
         K = K + 1
         IF (F0.LT.F1.AND.(X1*X2).GT.0.D0) GO TO 4
         X2 = 0.5D0*X2
         GO TO 11
C...INCREMENT THE ONE-DIMENSIONAL SEARCH COUNTER...
12    NL = NL + 1
      IF (F2.LE.FM) GO TO 13
      X2 = XM
      GO TO 14
13    FM = F2
C...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
14    IF (DABS(X2*(X2 - X1)).LE.SMALL) GO TO 15
         D2 = (X2*(F1-F0) - X1*(FM-F0))/((X1*X2)*(X1 - X2))
         GO TO 16
15       IF (K.GT.0) D2 = 0.D0
16    IF (D2.LE.SMALL) D2 = SMALL
      X1 = X2
      FX = FM
      IF (SF1.GE.FX) GO TO 17
         FX = SF1
         X1 = SX1
C...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
17    IF (J.EQ.0) RETURN
      DO 18 I=1,N
18       X(I) = X(I) + X1*V(I,J)
      RETURN
      END
      SUBROUTINE MINFIT(M,N,MACHEP,TOL,AB,Q)

c*********************************************************************72
      IMPLICIT double precision (A-H,O-Z)
      double precision MACHEP
      DIMENSION AB(M,N),Q(N)
C...AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
C   RESTRICTED TO M=N,P=0.
C   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
C   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
C   WHERE U IS ANOTHER ORTHOGONAL MATRIX.
      DIMENSION E(20)
C...HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM...
      IF (N.EQ.1) GO TO 200
      EPS = MACHEP
      G = 0.D0
      X = 0.D0
      DO 11 I=1,N
         E(I) = G
         S = 0.D0
         L = I + 1
         DO 1 J=I,N
1           S = S + AB(J,I)**2
         G = 0.D0
         IF (S.LT.TOL) GO TO 4
            F = AB(I,I)
            G = DSQRT(S)
            IF (F.GE.0.D0) G = -G
            H = F*G - S
            AB(I,I)=F-G
            IF (L.GT.N) GO TO 4
            DO 3 J=L,N
               F = 0.D0
               DO 2 K=I,N
2                 F = F + AB(K,I)*AB(K,J)
               F = F/H
               DO 3 K=I,N
3                 AB(K,J) = AB(K,J) + F*AB(K,I)
4        Q(I) = G
         S = 0.D0
         IF (I.EQ.N) GO TO 6
         DO 5 J=L,N
5           S = S + AB(I,J)*AB(I,J)
6        G = 0.D0
         IF (S.LT.TOL) GO TO 10
            IF (I.EQ.N) GO TO 16
            F = AB(I,I+1)
16          G = DSQRT(S)
            IF (F.GE.0.D0) G = -G
            H = F*G - S
            IF (I.EQ.N) GO TO 10
            AB(I,I+1) = F - G
            DO 7 J=L,N
7              E(J) = AB(I,J)/H
            DO 9 J=L,N
               S = 0.D0
               DO 8 K=L,N
8                 S = S + AB(J,K)*AB(I,K)
               DO 9 K=L,N
9                 AB(J,K) = AB(J,K) + S*E(K)
10       Y = DABS(Q(I)) + DABS(E(I))
11       IF (Y.GT.X) X = Y
C...ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS...
      AB(N,N) = 1.D0
      G = E(N)
      L = N
      DO 25 II=2,N
         I = N - II + 1
         IF (G.EQ.0.D0) GO TO 23
         H = AB(I,I+1)*G
         DO 20 J=L,N
20          AB(J,I) = AB(I,J)/H
         DO 22 J=L,N
            S = 0.D0
            DO 21 K=L,N
21             S = S + AB(I,K)*AB(K,J)
            DO 22 K=L,N
22             AB(K,J) = AB(K,J) + S*AB(K,I)
23       DO 24 J=L,N
            AB(I,J) = 0.D0
24          AB(J,I) = 0.D0
         AB(I,I) = 1.D0
         G = E(I)
25       L = I
C...DIAGONALIZATION OF THE BIDIAGONAL FORM...
100   EPS = EPS*X
      DO 150 KK=1,N
         K = N - KK + 1
         KT = 0
101      KT = KT + 1
         IF (KT.LE.30) GO TO 102
            E(K) = 0.D0
            WRITE (6,1000)
1000        FORMAT (' QR FAILED')
102      DO 103 LL2=1,K
            L2 = K - LL2 + 1
            L = L2
            IF (DABS(E(L)).LE.EPS) GO TO 120
            IF (L.EQ.1) GO TO 103
            IF (DABS(Q(L-1)).LE.EPS) GO TO 110
103         CONTINUE
C...CANCELLATION OF E(L) IF L>1...
110      C = 0.D0
         S = 1.D0
         DO 116 I=L,K
            F = S*E(I)
            E(I) = C*E(I)
            IF (DABS(F).LE.EPS) GO TO 120
            G = Q(I)
C...Q(I) = H = DSQRT(G*G + F*F)...
            IF (DABS(F).LT.DABS(G)) GO TO 113
            IF (F) 112,111,112
111         H = 0.D0
            GO TO 114
112         H = DABS(F)*DSQRT(1 + (G/F)**2)
            GO TO 114
113         H = DABS(G)*DSQRT(1 + (F/G)**2)
114         Q(I) = H
            IF (H.NE.0.D0) GO TO 115
               G = 1.D0
               H = 1.D0
115         C = G/H
116         S = -F/H
C...TEST FOR CONVERGENCE...
120      Z = Q(K)
         IF (L.EQ.K) GO TO 140
C...SHIFT FROM BOTTOM 2*2 MINOR...
         X = Q(L)
         Y = Q(K-1)
         G = E(K-1)
         H = E(K)
         F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(2*H*Y)
         G = DSQRT(F*F + 1.0D0)
         TEMP = F - G
         IF (F.GE.0.D0) TEMP = F + G
         F = ((X - Z)*(X + Z) + H*(Y/TEMP - H))/X
C...NEXT QR TRANSFORMATION...
         C = 1.D0
         S = 1.D0
         LP1 = L + 1
         IF (LP1.GT.K) GO TO 133
         DO 132 I=LP1,K
            G = E(I)
            Y = Q(I)
            H = S*G
            G = G*C
            IF (DABS(F).LT.DABS(H)) GO TO 123
            IF (F) 122,121,122
121         Z = 0.D0
            GO TO 124
122         Z = DABS(F)*DSQRT(1 + (H/F)**2)
            GO TO 124
123         Z = DABS(H)*DSQRT(1 + (F/H)**2)
124         E(I-1) = Z
            IF (Z.NE.0.D0) GO TO 125
               F = 1.D0
               Z = 1.D0
125         C = F/Z
            S = H/Z
            F = X*C + G*S
            G = -X*S + G*C
            H = Y*S
            Y = Y*C
            DO 126 J=1,N
               X = AB(J,I-1)
               Z = AB(J,I)
               AB(J,I-1) = X*C + Z*S
126            AB(J,I) = -X*S + Z*C
            IF (DABS(F).LT.DABS(H)) GO TO 129
            IF (F) 128,127,128
127         Z = 0.D0
            GO TO 130
128         Z = DABS(F)*DSQRT(1 + (H/F)**2)
            GO TO 130
129         Z = DABS(H)*DSQRT(1 + (F/H)**2)
130         Q(I-1) = Z
            IF (Z.NE.0.D0) GO TO 131
               F = 1.D0
               Z = 1.D0
131         C = F/Z
            S = H/Z
            F = C*G + S*Y
132         X = -S*G + C*Y
133      E(L) = 0.D0
         E(K) = F
         Q(K) = X
         GO TO 101
C...CONVERGENCE:  Q(K) IS MADE NON-NEGATIVE...
140      IF (Z.GE.0.D0) GO TO 150
         Q(K) = -Z
         DO 141 J=1,N
141         AB(J,K) = -AB(J,K)
150      CONTINUE
      RETURN
200   Q(1) = AB(1,1)
      AB(1,1) = 1.D0
      RETURN
      END
      SUBROUTINE PRINT(N,X,PRIN,FMIN)

c*********************************************************************72
      IMPLICIT double precision (A-H,O-Z)
      INTEGER PRIN
      double precision X(N),LN,LDT
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
      WRITE (6,101) NL,NF,FX
      IF (FX.LE.FMIN) GO TO 1
      LN = DLOG10(FX-FMIN)
      WRITE (6,102) FMIN,LN
      GO TO 2
1     WRITE (6,103) FMIN
2     IF (N.GT.4.AND.PRIN.LE.2) RETURN
      WRITE (6,104) (X(I),I=1,N)
      RETURN
101   FORMAT (/' AFTER',I6,
     . ' LINEAR SEARCHES, THE FUNCTION HAS BEEN EVALUATED',I6,
     . ' TIMES.  THE SMALLEST VALUE FOUND IS F(X) = ',E21.14)
102   FORMAT (' LOG (F(X)-',E21.14,') = ',E21.14)
103   FORMAT (' LOG (F(X)-',E21.14,') IS UNDEFINED.')
104   FORMAT (' X IS:',E26.14/(E32.14))
      END
      SUBROUTINE QUAD(N,F,X,T,MACHEP,H)

c*********************************************************************72
      IMPLICIT double precision (A-H,O-Z)
      EXTERNAL F
C...QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X...
      double precision X(N),MACHEP,LDT,L
      DIMENSION V(20,20),Q0(20),Q1(20)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
      S = FX
      FX = QF1
      QF1 = S
      QD1 = 0.D0
      DO 1 I=1,N
         S = X(I)
         L = Q1(I)
         X(I) = L
         Q1(I) = S
1        QD1 = QD1 + (S-L)**2
      QD1 = DSQRT(QD1)
      L = QD1
      S = 0.D0
      IF (QD0. LE. 0.D0 .OR. QD1 .LE. 0.D0 .OR. NL .LT. 3*N*N) GO TO 2
      VALUE=QF1
      CALL MIN(N,0,2,S,L,VALUE,.TRUE.,F,X,T,MACHEP,H)
      QA = (L*(L-QD1))/(QD0*(QD0+QD1))
      QB = ((L+QD0)*(QD1-L))/(QD0*QD1)
      QC = (L*(L+QD0))/(QD1*(QD0+QD1))
      GO TO 3
2     FX = QF1
      QA = 0.D0
      QB = QA
      QC = 1.D0
3     QD0 = QD1
      DO 4 I=1,N
         S = Q0(I)
         Q0(I) = X(I)
4        X(I) = (QA*S + QB*X(I)) + QC*Q1(I)
      RETURN
      END
      function random ( naught )

c*********************************************************************72
c
cc RANDOM is a random number generator.
c
      implicit none

      real*8 half
      integer i
      logical init
      integer j
      integer naught
      integer q
      integer r
      double precision ran1
      integer ran2
      double precision ran3(127)
      double precision random

      save init
      save ran2
      save ran3

      data init /.false./

      if ( .not. init) then
        r = mod(naught,8190) + 1
        ran2 = 128
        do i=1,127
          ran2 = ran2 - 1
          ran1 = -2.d0**55
          do j=1,7
            r = mod(1756*r,8191)
            q = r/32
            ran1 = (ran1 + q)*(1.0d0/256)
          end do
          ran3(ran2) = ran1
        end do
        init = .true.
      end if

      if (ran2.eq.1) then
        ran2 = 128
      end if

      ran2 = ran2 - 1
      ran1 = ran1 + ran3(ran2)
      half = .5d0

      if (ran1.ge.0.d0) then
        half = -half
      end if

      ran1 = ran1 + half
      ran3(ran2) = ran1
      random = ran1 + .5d0

      return
      end
      SUBROUTINE SORT(M,N,D,V)

c*********************************************************************72
      double precision D(N),V(M,N)
C...SORTS THE ELEMENTS OF D(N) INTO DESCENDING ORDER AND MOVES THE
C   CORRESPONDING COLUMNS OF V(N,N).
C   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM.
      double precision S
      IF (N.EQ.1) RETURN
      NM1 = N - 1
      DO 3 I = 1,NM1
         K=I
         S = D(I)
         IP1 = I + 1
         DO 1 J = IP1,N
            IF (D(J) .LE. S) GO TO 1
            K = J
            S = D(J)
1           CONTINUE
         IF (K .LE. I) GO TO 3
         D(K) = D(I)
         D(I) = S
         DO 2 J = 1,N
            S = V(J,I)
            V(J,I) = V(J,K)
2           V(J,K) = S
3        CONTINUE
      RETURN
      END
      SUBROUTINE VCPRNT(OPTION,V,N)

c*********************************************************************72
      double precision V(N)
      INTEGER OPTION
      GO TO (1,2,3,4),OPTION
1     WRITE (6,101) (V(I),I=1,N)
      RETURN
2     WRITE (6,102) (V(I),I=1,N)
      RETURN
3     WRITE (6,103) (V(I),I=1,N)
      RETURN
4     WRITE (6,104) (V(I),I=1,N)
      RETURN
101   FORMAT (/' THE SECOND DIFFERENCE ARRAY D(*) IS:'/
     .        (E32.14,4E25.14))
102   FORMAT (/' THE SCALE FACTORS ARE:'/(E32.14,4E25.14))
103   FORMAT (/' THE APPROXIMATING QUADRATIC FORM HAS THE PRINCIPAL VALU
     .ES:'/(E32.14,4E25.14))
104   FORMAT (/' X IS:',E26.14/(E32.14))
      END
